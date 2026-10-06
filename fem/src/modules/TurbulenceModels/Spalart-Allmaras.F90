!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 16 Nov 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for the Spalart-Allmaras turbulence model.
!> This is the vectorized/threaded implementation: Basis/dBasisdx are computed
!> for all Gauss points of an element at once via ElementInfoVec, the
!> nonlinear closure coefficients become per-Gauss-point arrays via plain
!> array syntax, and each bilinear/linear form maps onto one LinearForms_*
!> call (same approach as HeatSolve.F90 and IncompressibleNS.F90). Bulk
!> assembly is threaded over elements, each thread resolving its own material
!> ValueHandle_t slot (mirroring IncompressibleNS's NSHandles_t).
!>
!> Axisymmetric/cylindrical coordinates go through LocalMatrixScalar instead
!> of LocalMatrixVec -- a scalar, per-Gauss-point fallback carrying the same
!> metric-tensor math as Spalart-AllmarasLegacy.F90's own LocalMatrix (no
!> LinearForms equivalent for that), called serially -- same role as
!> HeatSolve.F90's own AxiSymmetric branch between its LocalMatrixVec and
!> LocalMatrix. LocalMatrixScalar also carries the legacy "Bubbles = True"
!> per-node scheme (doubling nd to 2*n), for the rare sif that sets a plain
!> nodal "Element" explicitly instead of letting SpalartAllmaras_Init0 default
!> it to a p-element bubble.
!>
!> The original scalar-element solver lives on in Spalart-AllmarasLegacy.F90
!> (subroutine SpalartAllmarasLegacy), reachable either directly by that name
!> or via "Legacy Assembly = Logical True" here (see SpalartAllmarasFront
!> below). With axisymmetric and the per-node scheme both covered above, that
!> path is now only needed by a sif that must reproduce the legacy solver's
!> exact historical numbers.
!> \ingroup Solvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Whether this solver should run the original scalar-element assembly
!> (Spalart-AllmarasLegacy.F90) instead of this file's own implementation.
!------------------------------------------------------------------------------
MODULE SpalartAllmarasFront
  USE DefUtils
  USE LoadMod, ONLY: ExecSolver
  IMPLICIT NONE IMPLICIT_EXTERNAL

CONTAINS

  FUNCTION LegacyAssembly( Solver ) RESULT( Legacy )
    TYPE(Solver_t) :: Solver
    LOGICAL :: Legacy, Found

    Legacy = ListGetLogical( Solver % Values, 'Legacy Assembly', Found )
  END FUNCTION LegacyAssembly

!------------------------------------------------------------------------------
!> Call one of SpalartAllmarasLegacy's entry points with this solver. The
!> name is resolved at run time, as the core resolves any solver, so this
!> file and Spalart-AllmarasLegacy.so stay independent of one another.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateToSpalartAllmarasLegacy( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'Spalart-AllmarasLegacy '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'SpalartAllmaras', &
        '"Legacy Assembly" was requested but "'//TRIM(Entry)//'" could not be found. '// &
        'Is Spalart-AllmarasLegacy.so installed beside this solver?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateToSpalartAllmarasLegacy

END MODULE SpalartAllmarasFront


!------------------------------------------------------------------------------
MODULE SpalartAllmarasLocalForms

  USE DefUtils
  USE LinearForms

  IMPLICIT NONE IMPLICIT_EXTERNAL

  ! Per-element bubble history, needed by Default1stOrderTime's Nb path
  ! (DefUtils.F90) to form a consistent BDF(1) time derivative for a
  ! condensed p-bubble (or, here, LocalMatrixScalar's legacy per-node
  ! "Bubbles" scheme, Nb=n) -- lives on Solver % Variable's own BubbleValues/
  ! BubblePrevValues (Types.F90), not a separate type; see the matching
  ! bx/bxprev comment in Spalart-AllmarasLegacy.F90, which this mirrors.

  ! Per-thread ValueHandle_t storage for LocalMatrixVec's material lookups.
  ! NOT THREADPRIVATE -- see the matching comment on IncompressibleNS.F90's
  ! NSHandles_t for the Windows/GCC emutls hazard that rules that out. Instead
  ! a plain module-level array indexed by omp_get_thread_num()+1; each thread
  ! only ever touches its own slot.
  TYPE :: SAHandles_t
    TYPE(ValueHandle_t) :: Visc_h, Dens_h
  END TYPE SAHandles_t
  TYPE(SAHandles_t), ALLOCATABLE, SAVE :: SAHandles(:)

CONTAINS

!------------------------------------------------------------------------------
!> Assemble and glue local matrix/RHS for one bulk element. Vectorized over
!> Gauss points, safe to call concurrently from multiple threads (each thread
!> passes its own InitHandles and only touches SAHandles(tid)).
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, dt, Transient, GlobalBubbles, Stabilize, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE IMPLICIT_EXTERNAL
    TYPE(Element_t), POINTER :: Element
    INTEGER, INTENT(IN) :: n, nd, nb
    REAL(KIND=dp), INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: Transient, GlobalBubbles, Stabilize
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    REAL(KIND=dp), ALLOCATABLE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:)
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), TimeForce(:)

    REAL(KIND=dp), POINTER :: RhoVec(:), MuVec(:)

    REAL(KIND=dp), ALLOCATABLE :: VeloNodal(:,:), TmuNodal(:), DistNodal(:)
    REAL(KIND=dp), ALLOCATABLE :: VeloVec(:,:), dVelodxVec(:,:,:), TmuVec(:), &
        GradTmuVec(:,:), DistVec(:), StrainVec(:,:,:), VorticityVec(:,:,:), &
        StrainMeasureVec(:), VorticityMeasureVec(:), XiVec(:), fw1Vec(:), &
        fw2Vec(:), StVec(:), rVec(:), gVec(:), fwVec(:), EffmuVec(:), &
        ReactCoeffVec(:), LoadVec(:), EffVeloVec(:,:), StreamVec(:,:), &
        TauVec(:), TmpVec(:), RadiusVec(:)

    REAL(KIND=dp) :: Cb1,Cb2,Cv1,Sigma,Cw1,Cw2,Cw3,hK,mK,VNorm
    INTEGER :: i,j,k,p,ngp,dim,allocstat,tid,ntot
    LOGICAL :: Stat, Found, IsAxiSymmetric
!------------------------------------------------------------------------------
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Visc_h => SAHandles(tid) % Visc_h, Dens_h => SAHandles(tid) % Dens_h )

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity' )
      CALL ListInitElementKeyword( Dens_h,'Material','Density' )
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ! nd is the RETAINED (nodal/edge/face) dof count, exactly
    ! GetElementNOFDOFs()'s own meaning -- nb (p-bubble) dofs are separate.
    ! The test/trial basis actually used in the local matrix spans both, so
    ! everything below sized off the bubble-augmented total uses ntot, while
    ! CondensateP/CondensatePTransient (which need to know where the retained
    ! block ends) still get nd and nb apart. Mirrors the legacy driver's own
    ! "nd+nb" passed into LocalMatrix's "nd" parameter
    ! (Spalart-AllmarasLegacy.F90) and IncompressibleNS's ntot=nd+nb passed
    ! alongside its own plain nd.
    ntot = nd + nb

    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ALLOCATE( BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        MASS(ntot,ntot), STIFF(ntot,ntot), FORCE(ntot), TimeForce(ntot), &
        VeloNodal(3,n), TmuNodal(n), DistNodal(n), &
        VeloVec(ngp,3), dVelodxVec(ngp,3,3), TmuVec(ngp), GradTmuVec(ngp,3), &
        DistVec(ngp), StrainVec(ngp,3,3), VorticityVec(ngp,3,3), &
        StrainMeasureVec(ngp), VorticityMeasureVec(ngp), XiVec(ngp), &
        fw1Vec(ngp), fw2Vec(ngp), StVec(ngp), rVec(ngp), gVec(ngp), &
        fwVec(ngp), EffmuVec(ngp), ReactCoeffVec(ngp), LoadVec(ngp), &
        EffVeloVec(ngp,3), StreamVec(ngp,ntot), TauVec(ngp), TmpVec(ngp), &
        RadiusVec(ngp), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('SpalartAllmaras','Local storage allocation failed')

    CALL GetElementNodesVec( Nodes, UElement=Element )

    MASS = 0._dp; STIFF = 0._dp; FORCE = 0._dp

    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, DetJVec, &
        SIZE(BasisVec,2), BasisVec, dBasisdxVec )
    DetJVec(1:ngp) = DetJVec(1:ngp) * IP % s(1:ngp)

    ! Axisymmetric (no swirl): the integration measure is r dr dz instead of
    ! dr dz -- the diffusion/reaction operator itself needs no further metric
    ! correction, since CoordinateSystemInfo's Cylindrical Metric stays
    ! diagonal-identity in the r,z block (Metric(3,3)=1/r^2 only matters for a
    ! genuine theta/swirl derivative, which a 2D mesh's dBasisdx never has).
    ! This much mirrors HeatSolve.F90's own "Weight = Weight * r" exactly, and
    ! genuine swirl ("Cylindric Symmetric") still goes through
    ! LocalMatrixScalar.
    !
    ! The strain-rate measure below is a SEPARATE story: unlike an ordinary
    ! partial derivative, the covariant strain e_theta_theta = u_r/r is
    ! nonzero even with no swirl and no theta-dependence -- see
    ! SecondInvariant's own dedicated AxisSymmetric branch
    ! (MaterialModels.F90), whose last term "(2*Velo(1)*symb(1,3,3))**2" is
    ! exactly this, since symb(1,3,3) = 1/r. RadiusVec is kept around to add
    ! that same term to StrainMeasureVec further down.
    IsAxiSymmetric = ( CurrentCoordinateSystem() == AxisSymmetric )
    IF( IsAxiSymmetric ) THEN
      RadiusVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), Nodes % x(1:n) )
      DetJVec(1:ngp) = DetJVec(1:ngp) * RadiusVec(1:ngp)
    END IF

    ! Nodal input fields, at the n element corner nodes -- exactly the arrays
    ! the legacy LocalMatrix samples (UX/UY/UZ/Tviscosity/Distance are all
    ! declared (n), not (nd), there too): a p-bubble mode never participates
    ! in the closure-law evaluation, only in the test/trial basis.
    VeloNodal = 0._dp
    CALL GetScalarLocalSolution( VeloNodal(1,1:n), 'Velocity 1', UElement=Element )
    CALL GetScalarLocalSolution( VeloNodal(2,1:n), 'Velocity 2', UElement=Element )
    IF( dim == 3 ) CALL GetScalarLocalSolution( VeloNodal(3,1:n), 'Velocity 3', UElement=Element )

    CALL GetScalarLocalSolution( TmuNodal, UElement=Element )
    CALL GetScalarLocalSolution( DistNodal, 'Wall Distance', UElement=Element )

    RhoVec => ListGetElementRealVec( Dens_h, ngp, BasisVec, Element, Found )
    MuVec  => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element, Found )

    VeloVec = 0._dp
    DO i=1,dim
      VeloVec(1:ngp,i) = MATMUL( BasisVec(1:ngp,1:n), VeloNodal(i,1:n) )
    END DO

    dVelodxVec = 0._dp
    DO i=1,dim
      DO k=1,dim
        dVelodxVec(1:ngp,i,k) = MATMUL( dBasisdxVec(1:ngp,1:n,k), VeloNodal(i,1:n) )
      END DO
    END DO

    TmuVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), TmuNodal(1:n) )

    GradTmuVec = 0._dp
    DO k=1,dim
      GradTmuVec(1:ngp,k) = MATMUL( dBasisdxVec(1:ngp,1:n,k), TmuNodal(1:n) )
    END DO

    DistVec(1:ngp) = MAX( MATMUL( BasisVec(1:ngp,1:n), DistNodal(1:n) ), 1.0d-10 )

    ! Strain-rate / vorticity tensors and their (Frobenius) measures, Cartesian
    ! only -- mirrors the legacy CurrentCoordinateSystem()==Cartesian branch.
    StrainVec = 0._dp; VorticityVec = 0._dp
    DO i=1,dim
      DO k=1,dim
        StrainVec(1:ngp,i,k)    = 0.5_dp*( dVelodxVec(1:ngp,i,k) + dVelodxVec(1:ngp,k,i) )
        VorticityVec(1:ngp,i,k) = 0.5_dp*( dVelodxVec(1:ngp,i,k) - dVelodxVec(1:ngp,k,i) )
      END DO
    END DO

    StrainMeasureVec = 0._dp; VorticityMeasureVec = 0._dp
    DO i=1,dim
      DO k=1,dim
        StrainMeasureVec(1:ngp)    = StrainMeasureVec(1:ngp)    + StrainVec(1:ngp,i,k)**2
        VorticityMeasureVec(1:ngp) = VorticityMeasureVec(1:ngp) + VorticityVec(1:ngp,i,k)**2
      END DO
    END DO

    ! Axisymmetric (no swirl) hoop strain e_theta_theta = u_r/r: a genuine
    ! extra diagonal strain component (see the comment on RadiusVec above),
    ! entering StrainMeasureVec's sum of squares exactly like any other
    ! diagonal Strain(i,i)**2 term above. Vorticity needs no such addition --
    ! its diagonal is zero by antisymmetry, and there is no swirl velocity to
    ! give it an off-diagonal r-theta/z-theta component either.
    IF( IsAxiSymmetric ) THEN
      StrainMeasureVec(1:ngp) = StrainMeasureVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF

    StrainMeasureVec(1:ngp)    = SQRT( 2._dp*StrainMeasureVec(1:ngp) )
    VorticityMeasureVec(1:ngp) = SQRT( 2._dp*VorticityMeasureVec(1:ngp) )

    ! Spalart-Allmaras closure coefficients -- same constants and formulas as
    ! Spalart-AllmarasLegacy.F90's LocalMatrix, just evaluated for all ngp
    ! points at once via array syntax instead of one point at a time. The
    ! rotation and curvature correction is left out here exactly as it is
    ! there (r=1, "NOT IN USE").
    Cb1 = 0.1355_dp; Cb2 = 0.6220_dp; Cv1 = 7.1_dp; Sigma = 2._dp/3._dp
    Cw1 = Cb1/0.41_dp**2 + (1._dp+Cb2)/Sigma
    Cw2 = 0.3_dp; Cw3 = 2.0_dp

    XiVec(1:ngp)  = TmuVec(1:ngp) / ( MuVec(1:ngp)/RhoVec(1:ngp) )
    fw1Vec(1:ngp) = XiVec(1:ngp)**3 / ( XiVec(1:ngp)**3 + Cv1**3 )
    fw2Vec(1:ngp) = 1._dp - XiVec(1:ngp) / ( 1._dp + XiVec(1:ngp)*fw1Vec(1:ngp) )

    StVec(1:ngp) = VorticityMeasureVec(1:ngp) + &
        2._dp*MIN( 0._dp, StrainMeasureVec(1:ngp) - VorticityMeasureVec(1:ngp) )
    StVec(1:ngp) = StVec(1:ngp) + TmuVec(1:ngp)/DistVec(1:ngp)**2/0.41_dp**2*fw2Vec(1:ngp)

    rVec(1:ngp) = TmuVec(1:ngp) / MAX(StVec(1:ngp),1.0d-10) / 0.41_dp**2 / DistVec(1:ngp)**2
    gVec(1:ngp) = rVec(1:ngp) + Cw2*( rVec(1:ngp)**6 - rVec(1:ngp) )
    fwVec(1:ngp) = gVec(1:ngp) * ( (1._dp+Cw3**6) / (gVec(1:ngp)**6+Cw3**6) )**(1._dp/6._dp)

    EffmuVec(1:ngp) = ( MuVec(1:ngp) + RhoVec(1:ngp)*TmuVec(1:ngp) ) / Sigma

    ! Reaction term (production - destruction, linearized about the previous
    ! iterate) and the RHS load, exactly as the legacy A/LoadAtIp terms.
    ReactCoeffVec(1:ngp) = RhoVec(1:ngp) * ( -Cb1*StVec(1:ngp)/4._dp + &
        Cw1*fwVec(1:ngp)*TmuVec(1:ngp)/DistVec(1:ngp)**2 )

    LoadVec(1:ngp) = 3._dp*RhoVec(1:ngp)*Cb1*StVec(1:ngp)*TmuVec(1:ngp)/4._dp

    ! Effective convection velocity: Velo - (Cb2/Sigma)*grad(Tmu).
    EffVeloVec = 0._dp
    DO i=1,dim
      EffVeloVec(1:ngp,i) = VeloVec(1:ngp,i) - (Cb2/Sigma)*GradTmuVec(1:ngp,i)
    END DO

    ! MASS(p,q) = (rho*Basis(q), Basis(p))
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, MASS, RhoVec )

    ! STIFF(p,q) += (reaction coefficient * Basis(q), Basis(p))
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, STIFF, ReactCoeffVec )

    ! STIFF(p,q) += (Effmu * grad Basis(q), grad Basis(p))  -- diffusion
    CALL LinearForms_GradUdotGradU( ngp, ntot, dim, dBasisdxVec, DetJVec, STIFF, EffmuVec )

    ! STIFF(p,q) += (rho*EffVelo . grad Basis(q), Basis(p))  -- convection
    CALL LinearForms_GradUdotU( ngp, ntot, dim, dBasisdxVec, BasisVec, DetJVec, STIFF, &
        RhoVec, EffVeloVec )

    ! FORCE(p) += (LoadAtIp, Basis(p))
    CALL LinearForms_UdotF( ngp, ntot, BasisVec, DetJVec, LoadVec, FORCE )

    !------------------------------------------------------------------------
    ! SUPG (equal-order) stabilization, opt-in via "Stabilize"/"Stabilization
    ! Method" (see GetStabilizeFlag). Same Franca et al. tau construction and
    ! streamline-weighted test function as HeatSolve.F90's own Vec SUPG
    ! block, adapted to this equation's coefficients: rho multiplies both the
    ! time derivative and the convection term here (rho*cp does in Heat), and
    ! Effmu is the diffusion coefficient (Heat Conductivity there). Like
    ! HeatSolve's version, the reaction term (ReactCoeffVec) is left out of
    ! the streamline residual -- it is not needed for the pure convection-
    ! diffusion stability mechanism SUPG targets, and keeping it out matches
    ! the precedent (HeatSolve drops its own C0/reaction piece for the same
    ! reason).
    !------------------------------------------------------------------------
    IF( Stabilize ) THEN
      hK = Element % hK
      mK = Element % StabilizationMK

      ! Streamline-weighted trial/test "basis": rho*(EffVelo . grad Basis(p))
      StreamVec(1:ngp,1:ntot) = 0._dp
      DO i=1,dim
        DO p=1,ntot
          StreamVec(1:ngp,p) = StreamVec(1:ngp,p) + &
              RhoVec(1:ngp) * EffVeloVec(1:ngp,i) * dBasisdxVec(1:ngp,p,i)
        END DO
      END DO

      DO j=1,ngp
        VNorm = SQRT( SUM( EffVeloVec(j,1:dim)**2 ) )
        IF( VNorm > 0._dp .AND. EffmuVec(j) /= 0._dp ) THEN
          TauVec(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(EffmuVec(j))) )
          TauVec(j) = hK * TauVec(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauVec(j) = 0._dp
        END IF
      END DO

      ! STIFF(p,q) += tau*(StreamVec(q), StreamVec(p))
      TmpVec(1:ngp) = TauVec(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, STIFF )

      ! FORCE(p) += tau*(LoadAtIp, StreamVec(p))
      CALL LinearForms_UdotF( ngp, ntot, StreamVec, TmpVec, LoadVec, FORCE )

      IF( Transient ) THEN
        ! MASS(p,q) += tau*rho*(Basis(q), StreamVec(p))
        TmpVec(1:ngp) = TauVec(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MASS )
      END IF
    END IF

    !------------------------------------------------------------------------
    ! Time discretization and p-bubble condensation -- mirrors the nb>0
    ! branch of Spalart-AllmarasLegacy.F90's driver exactly (no legacy
    ! "Bubbles" per-node scheme here, see the file header).
    !------------------------------------------------------------------------
    TimeForce = 0._dp
    IF( nb > 0 ) THEN
      IF( Transient .AND. .NOT. GlobalBubbles ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=nb )
      ELSE
        IF( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
        CALL CondensateP( nd, nb, STIFF, FORCE, TimeForce )
      END IF
    ELSE
      IF( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element, VecAssembly=.TRUE. )

    END ASSOCIATE
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Scalar, coordinate-system-aware fallback for what LocalMatrixVec cannot
!> do: axisymmetric/cylindrical coordinates (the ElementKernel below carries
!> the same metric-tensor branch as Spalart-AllmarasLegacy.F90's own
!> LocalMatrix -- verbatim, no new math), and the legacy "Bubbles = True"
!> per-node scheme for a plain nodal "Element" set explicitly. Always called
!> serially (see the AxiSymmetric branch in SpalartAllmaras below), so it
!> uses the classic GetMaterial()/GetReal()/argument-less GetElementNOF*()
!> accessors exactly as the legacy driver does, relying on the same-thread
!> "current element" state its own preceding GetActiveElement() call set --
!> not safe to call from inside an OMP parallel region, unlike LocalMatrixVec.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixScalar( Element, dt, Transient, GlobalBubbles, BubblesDefault )
!------------------------------------------------------------------------------
    IMPLICIT NONE IMPLICIT_EXTERNAL
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: Transient, GlobalBubbles, BubblesDefault
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:), &
        TimeForce(:)
    LOGICAL :: Bubbles
    INTEGER :: n, nd, nb, allocstat
!------------------------------------------------------------------------------
    Bubbles = BubblesDefault .AND. .NOT. ASSOCIATED( Element % PDefs )
    Material => GetMaterial()

    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    IF ( Bubbles ) nd = 2*n
    nb = GetElementNOFBDOFs()
    CALL GetElementNodes( ElementNodes )

    ALLOCATE( MASS(nd+nb,nd+nb), STIFF(nd+nb,nd+nb), FORCE(nd+nb), LOAD(1,n), &
        TimeForce(nd+nb), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('SpalartAllmaras','Local storage allocation failed')

    CALL ElementKernel( MASS, STIFF, FORCE, LOAD, Element, n, nd+nb, ElementNodes )

    TimeForce = 0.0_dp
    IF ( Bubbles ) THEN
      IF ( Transient ) THEN
        ! Same convention as the nb>0 branch below, just with "as many
        ! bubbles as nodes" (Nb=n). DOFs=1 here, so no interleaving is needed.
        ! Never gated on Solver % GlobalBubbles -- legacy bubbles are always
        ! locally condensed, same reasoning as Spalart-AllmarasLegacy.F90's
        ! driver.
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=n )
      ELSE
        CALL Condensate( n, STIFF, FORCE, TimeForce )
      END IF
    ELSE IF ( nb > 0 ) THEN
      IF ( Transient .AND. .NOT. GlobalBubbles ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=nb )
      ELSE
        IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
        CALL CondensateP( nd, nb, STIFF, FORCE, TimeForce )
      END IF
    ELSE
      IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )

  CONTAINS

!------------------------------------------------------------------------------
!> Verbatim port of Spalart-AllmarasLegacy.F90's nested LocalMatrix: same
!> per-Gauss-point scalar math, same metric-tensor (axisymmetric/cylindrical)
!> branch, same closure-law formulas. "Bubbles" and "Material" come from the
!> host (LocalMatrixScalar above), exactly as they did from the legacy
!> driver's own host scope.
!------------------------------------------------------------------------------
    SUBROUTINE ElementKernel( MASS,STIFF,FORCE, LOAD, Element,n,nd,Nodes )
!------------------------------------------------------------------------------
      USE MaterialModels

      IMPLICIT NONE IMPLICIT_EXTERNAL

      REAL(KIND=dp), DIMENSION(:)   :: FORCE
      REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

      INTEGER :: n, nd

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: ddBasisddx(nd,3,3)
      REAL(KIND=dp) :: Basis(nd)
      REAL(KIND=dp) :: dBasisdx(nd,3),detJ

      REAL(KIND=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3),Tviscosity(n), &
                       Distance(n), Density(n), Viscosity(n)

      REAL(KIND=dp) :: A,M,Prod,div
      INTEGER :: i,j,c,p,q,t,dim,NBasis
      REAL(KIND=dp) :: LoadatIp,Cmu,Rho,mu,Tmu,Effmu

      REAL(KIND=dp) :: s,u,v,w,Strain(3,3), Vorticity(3,3), dist

      REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,Sigma, &
         GradTmu(3), Cw1,Cw2,Cw3,fw,fw1,fw2,Cb1,Cb2,Cb3,St,Xi,Cv1,g,r

      REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

      LOGICAL :: stat
      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      Viscosity(1:n) = GetReal( Material, 'Viscosity' )
      Density(1:n) = GetReal( Material, 'Density' )

      CALL GetScalarLocalSolution( UX, 'Velocity 1' )
      CALL GetScalarLocalSolution( UY, 'Velocity 2' )
      CALL GetScalarLocalSolution( UZ, 'Velocity 3' )

      CALL GetScalarLocalSolution( TViscosity )
      CALL GetScalarLocalSolution( Distance, 'Wall Distance' )

      FORCE = 0.0_dp
      STIFF = 0.0_dp
      MASS  = 0.0_dp

      NBasis = nd

      IF ( Bubbles ) THEN
         IntegStuff = GaussPoints( element, element % Type % GaussPoints2 )
      ELSE
         IntegStuff = GaussPoints( element )
      END IF

      DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)
        stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
              Basis,dBasisdx,Bubbles=Bubbles )

        s = detJ * IntegStuff % s(t)
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( Nodes % x(1:n)*Basis(1:n) )
          y = SUM( Nodes % y(1:n)*Basis(1:n) )
          z = SUM( nodes % z(1:n)*Basis(1:n) )
          CALL CoordinateSystemInfo(Metric,SqrtMetric,Symb,dSymb,X,Y,Z)
          s = s * SqrtMetric
        END IF

        Velo = 0.0_dp
        Velo(1) = SUM( UX(1:n)*Basis(1:n) )
        Velo(2) = SUM( UY(1:n)*Basis(1:n) )
        Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

        dVelodx = 0.0_dp
        DO i=1,dim
          dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
          dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
          dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
        END DO

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
          Strain  = 0.5_dp * (dVelodx + TRANSPOSE(dVelodx))
          StrainMeasure = SQRT(2 * SUM(Strain * Strain))

          Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
          VorticityMeasure = SQRT(2 * SUM(Vorticity * Vorticity))
        ELSE
          StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
        END IF

        Tmu = SUM( Tviscosity(1:n) * Basis(1:n) )
        DO i=1,dim
          GradTmu(i) = SUM( Tviscosity(1:n) * dBasisdx(1:n,i) )
        END DO

        mu   = SUM( Viscosity(1:n) * Basis(1:n) )
        rho  = SUM( Density(1:n) * Basis(1:n) )
        dist = MAX( SUM( Distance(1:n) * Basis(1:n) ), 1.0d-10 )

        Cb1 = 0.1355_dp
        Cb2 = 0.6220_dp
        Cv1 = 7.1_dp
        Sigma = 2._dp/3._dp

        ! Rotation and curvature correction of Schweighofer & Helsten; NOT IN USE
        r = 1._dp
        Cw1 = r*(Cb1/0.41_dp**2 + (1+Cb2)/Sigma)

        Cw2 = 0.3_dp
        Cw3 = 2.0_dp

        Xi  = Tmu/(mu/rho)
        fw1 = Xi**3 / (Xi**3 + Cv1**3)
        fw2 = 1 - Xi / ( 1+Xi*fw1 )

        St = VorticityMeasure + 2 * MIN(0.0_dp, StrainMeasure-VorticityMeasure)
        St = St + Tmu / dist**2 / 0.41_dp**2 * fw2

        r  = Tmu / MAX( St, 1.0d-10 ) / 0.41_dp**2 / dist**2
        g  = r + Cw2 * (r**6-r)
        fw = g*((1+Cw3**6)/(g**6+Cw3**6))**(1._dp/6._dp)

        Effmu = (mu + rho*Tmu)/Sigma

        DO p=1,NBasis
        DO q=1,NBasis
           M = 0.0d0
           A = 0.0d0

           M = rho * Basis(q) * Basis(p)
           A = A - rho * Cb1 * St * Basis(q) * Basis(p)/4
           A = A + rho * Cw1 * fw * Tmu / dist**2 * Basis(q) * Basis(p)

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO i=1,dim
                A = A + Effmu * dBasisdx(q,i) * dBasisdx(p,i)
              END DO
           ELSE
              DO i=1,dim
                DO j=1,dim
                   A = A + Metric(i,j) * Effmu * dBasisdx(q,i) * dBasisdx(p,j)
                END DO
              END DO
           END IF

           DO i=1,dim
             A = A + rho * (Velo(i)-Cb2*GradTmu(i)/Sigma) * dBasisdx(q,i) * Basis(p)
           END DO

           MASS(p,q)  = MASS(p,q)  + s*M
           STIFF(p,q) = STIFF(p,q) + s*A
         END DO
         END DO

         LoadAtIp = 3._dp*rho * Cb1 * St * Tmu/4

         DO p=1,NBasis
           FORCE(p) = FORCE(p)+s*LoadAtIp*Basis(p)
         END DO
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE ElementKernel
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixScalar
!------------------------------------------------------------------------------

END MODULE SpalartAllmarasLocalForms


!------------------------------------------------------------------------------
!> Vectorized/threaded Spalart-Allmaras driver. See the file header above for
!> what this does and does not support, and "Legacy Assembly" for the fallback.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SpalartAllmaras( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SpalartAllmarasLocalForms
  USE SpalartAllmarasFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: KE
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC
  INTEGER :: i,k,n,nb,nd,Active,iter,NonlinearIter,nthr
  LOGICAL :: GotIt, InitHandles, GlobalBubbles, Stabilize, AxiSymmetric, BubblesDefault
  REAL(KIND=dp) :: Norm, KVal
  CHARACTER(*), PARAMETER :: Caller = 'SpalartAllmaras'
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToSpalartAllmarasLegacy( 'SpalartAllmarasLegacy', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  KE => Solver % Variable
  IF ( .NOT. ASSOCIATED(KE) ) RETURN
  IF ( COUNT( KE % Perm > 0 ) <= 0 ) RETURN

  ! LocalMatrixVec now carries the plain "Axi Symmetric" (no swirl) case
  ! itself (r-weighted measure plus the hoop strain term -- see its own
  ! comments); genuine swirl ("Cylindric Symmetric") and general "Cylindric"
  ! still need the full metric/Christoffel treatment only LocalMatrixScalar
  ! has, so those still go through it, serially -- same branch HeatSolve.F90
  ! makes between its own LocalMatrixVec/LocalMatrix, just with one more
  ! coordinate system (AxisSymmetric) now on the fast side of it.
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian .AND. &
                   CurrentCoordinateSystem() /= AxisSymmetric )

  IF (.NOT. ALLOCATED(SAHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(SAHandles(nthr))
  END IF

  GlobalBubbles = Solver % GlobalBubbles
  Stabilize = GetStabilizeFlag( Solver % Values, GotIt )

  ! Only LocalMatrixScalar's legacy per-node branch uses this -- see the
  ! matching BubblesDefault resolution in Spalart-AllmarasLegacy.F90.
  BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  IF ( .NOT.GotIt ) BubblesDefault = .TRUE.

  ! A single scalar field, Dofs=1 -- exactly as Spalart-AllmarasLegacy.F90
  ! sizes its own bx/bxprev.
  IF ( TransientSimulation ) CALL DefaultBubbleHistoryUpdate( Dofs=1 )

  NonlinearIter = ListGetInteger( Solver % Values, 'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT.GotIt ) NonlinearIter = 1

  DO i=1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    IF ( ListGetLogical( BC, 'Noslip wall BC', gotit ) ) THEN
      CALL ListAddConstReal( BC, 'Turbulent Viscosity', 0.0_dp )
    END IF
  END DO

  DO iter=1,NonlinearIter
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)
    CALL Info(Caller,'Spalart-Allmaras iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)

    CALL DefaultInitialize()

    Active = GetNOFActive()

    IF ( AxiSymmetric ) THEN
      ! Serial only: LocalMatrixScalar relies on the classic argument-less
      ! Get*() "current element" accessors, exactly as
      ! Spalart-AllmarasLegacy.F90's own driver does.
      DO i=1,Active
        Element => GetActiveElement(i)
        CALL LocalMatrixScalar( Element, dt, TransientSimulation, GlobalBubbles, BubblesDefault )
      END DO
    ELSE
      ! Element 1 serially, so a Fatal/Warn triggered while resolving handles
      ! for the first time surfaces cleanly before any thread starts -- same
      ! reasoning as IncompressibleNS's own element-1 pass.
      InitHandles = .TRUE.
      Element => GetActiveElement(1)
      n  = GetElementNOFNodes(Element)
      nb = GetElementNOFBDOFs(Element, Update=.TRUE.)
      nd = GetElementNOFDOFs(Element)
      CALL LocalMatrixVec( Element, n, nd, nb, dt, TransientSimulation, GlobalBubbles, Stabilize, InitHandles )

      ! Each thread resolves its own handles rather than inheriting slot 1's --
      ! see the matching comment in IncompressibleNS.F90 on why a shallow copy
      ! of a ValueHandle_t (whose buffers are POINTERs) across threads is unsafe.
      InitHandles = .TRUE.
      !$OMP PARALLEL SHARED(Active, dt, TransientSimulation, GlobalBubbles, Stabilize) &
      !$OMP          PRIVATE(Element, i, n, nd, nb) FIRSTPRIVATE(InitHandles) DEFAULT(NONE)
      !$OMP DO
      DO i=2,Active
        Element => GetActiveElement(i)
        n  = GetElementNOFNodes(Element)
        nb = GetElementNOFBDOFs(Element, Update=.TRUE.)
        nd = GetElementNOFDOFs(Element)
        CALL LocalMatrixVec( Element, n, nd, nb, dt, TransientSimulation, GlobalBubbles, Stabilize, InitHandles )
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    END IF

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    ! Turbulent viscosity should stay positive.
    DO i=1,SIZE(Solver % Variable % Perm)
      k = Solver % Variable % Perm(i)
      IF ( k <= 0 ) CYCLE
      KVal = Solver % Variable % Values(k)
      Solver % Variable % Values(k) = MAX( KVal, 1.0d-12 )
    END DO

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE SpalartAllmaras
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Element-type resolution, called before the mesh/basis functions are
!> finalized -- same role and timing as HeatSolver_Init0. Legacy's default
!> ("Bubbles = True" whenever "Bubbles"/"Element" is left unset) stabilizes
!> with a per-node bubble that needs no p-element setup at all; this file's
!> only stabilization mechanisms (the p-element bubble, or SUPG) both DO need
!> the mesh built with the right "Element" string, so that default has to be
!> picked here, this early, rather than in SpalartAllmaras_Init below. Skipped
!> under "Legacy Assembly", since that mode never uses either.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SpalartAllmaras_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SpalartAllmarasFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, Serendipity, Stabilize
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) RETURN

  Params => GetSolverParams()

  ! Resolve whether the transport equation will be stabilized by SUPG
  ! (equal-order, no bubble) or by a residual-free bubble (the default) --
  ! same choice and same keyword as HeatSolver_Init0's.
  Stabilize = GetStabilizeFlag( Params )

  IF( .NOT. ListCheckPresent( Params,'Element' ) ) THEN
    IF( Stabilize ) THEN
      ! SUPG is the equal-order alternative to the bubble: a plain linear
      ! nodal element on every family, no bubble to condense at all.
      CALL ListAddNewString( Params,'Element','n:1' )
    ELSE
      Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found )
      IF(.NOT.Found) Serendipity = .TRUE.
      IF( Serendipity ) THEN
        CALL ListAddString( Params,'Element', &
            'p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4' )
      ELSE
        CALL ListAddString( Params,'Element', &
            'p:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:4 -pyramid b:4' )
      END IF
      CALL ListAddNewLogical( Params,'Bubbles in Global System',.FALSE. )
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE SpalartAllmaras_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: SpalartAllmaras. Under "Legacy
!> Assembly", delegates to SpalartAllmarasLegacy_Init and returns -- the
!> p-bubble setup below is specific to this file's own implementation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SpalartAllmaras_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SpalartAllmarasFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, PBubble
  CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToSpalartAllmarasLegacy( 'SpalartAllmarasLegacy_Init', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  SolverParams => GetSolverParams()

  IF ( ListGetLogical( SolverParams, 'Bubbles', Found ) ) THEN
    CALL Warn('SpalartAllmaras_Init', &
        '"Bubbles = True" (the legacy per-node scheme) is not used here -- '// &
        'SpalartAllmaras_Init0 already defaulted "Element" to an equivalent '// &
        'p-element bubble unless SUPG ("Stabilize"/"Stabilization Method") was '// &
        'requested instead. Use "Legacy Assembly = Logical True" for the original '// &
        'per-node scheme.')
  END IF

  ! Same reasoning as Spalart-AllmarasLegacy.F90's own Init: condense a
  ! p-bubble out locally by default (ListAddNew, so an explicit sif setting
  ! still wins), and a transient condensed bubble needs at least two solves
  ! per timestep to recover its history -- see the bubble history comment above.
  str = ListGetString( SolverParams,'Element', Found )
  PBubble = .FALSE.
  IF ( Found ) PBubble = INDEX( str, 'b:' ) > 0

  IF ( PBubble ) THEN
    CALL ListAddNewLogical(SolverParams, 'Bubbles in Global System', .FALSE.)

    IF ( TransientSimulation .AND. &
         .NOT. ListGetLogical(SolverParams,'Bubbles in Global System',Found) ) THEN
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Min Iterations', 2)
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Max Iterations', 2)
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE SpalartAllmaras_Init
!------------------------------------------------------------------------------
