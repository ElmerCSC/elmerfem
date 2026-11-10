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
! *  Original Date: 10 Nov 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for the K-omega turbulence model. Vectorized/threaded implementation
!> following the same pattern as Spalart-Allmaras.F90: Basis/dBasisdx computed
!> for all Gauss points at once via ElementInfoVec, closure coefficients as
!> per-Gauss-point arrays, each bilinear/linear form mapped onto one
!> LinearForms_* call. K and omega have no matrix cross-coupling in the legacy
!> LocalMatrix (only their nonlinear coefficients and the production term
!> couple the two) -- see A(1,1)/A(2,2) and M(1,1)/M(2,2) there, A(1,2)/A(2,1)
!> are never written -- so this assembles K and omega as two independent
!> scalar problems (LocalMatrixVec's "K"/"O" suffixed blocks) and interleaves
!> them into the final 2*ntot-sized MASS/STIFF/FORCE right before
!> condensation/glue, rather than needing IncompressibleNS's full block-
!> coupling machinery.
!>
!> Axisymmetric/cylindrical coordinates go through LocalMatrixScalar instead
!> of LocalMatrixVec -- a scalar, per-Gauss-point fallback carrying the same
!> metric-tensor math as KomegaLegacy.F90's own LocalMatrix (already
!> interleaved directly via STIFF(2*(p-1)+i,2*(q-1)+j), so ported verbatim),
!> called serially -- same role as HeatSolve.F90's own AxiSymmetric branch.
!> LocalMatrixScalar also carries the legacy "Bubbles = True" per-node scheme,
!> for a sif that sets a plain nodal "Element" explicitly, and OmegaWall, the
!> near-wall Dirichlet omega treatment (ported unchanged bar taking Solver as
!> an explicit argument, since it is no longer nested inside the driver).
!>
!> The original scalar-element solver lives on in KomegaLegacy.F90 (subroutine
!> KOmegaLegacy), reachable either directly by that name or via "Legacy
!> Assembly = Logical True" here (see KOmegaFront below). With axisymmetric
!> and the per-node scheme both covered above, that path is now only needed
!> by a sif that must reproduce the legacy solver's exact historical numbers.
!> \ingroup Solvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Whether this solver should run the original scalar-element assembly
!> (KomegaLegacy.F90) instead of this file's own implementation.
!------------------------------------------------------------------------------
MODULE KOmegaFront
  USE DefUtils
  USE LoadMod, ONLY: ExecSolver
  IMPLICIT NONE

CONTAINS

  FUNCTION LegacyAssembly( Solver ) RESULT( Legacy )
    TYPE(Solver_t) :: Solver
    LOGICAL :: Legacy, Found

    Legacy = ListGetLogical( Solver % Values, 'Legacy Assembly', Found )
  END FUNCTION LegacyAssembly

!------------------------------------------------------------------------------
!> Call one of KOmegaLegacy's entry points with this solver. The name is
!> resolved at run time, as the core resolves any solver, so this file and
!> KomegaLegacy.so stay independent of one another.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateToKOmegaLegacy( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'KomegaLegacy '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'KOmega', &
        '"Legacy Assembly" was requested but "'//TRIM(Entry)//'" could not be found. '// &
        'Is KomegaLegacy.so installed beside this solver?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateToKOmegaLegacy

END MODULE KOmegaFront


!------------------------------------------------------------------------------
MODULE KOmegaLocalForms

  USE DefUtils
  USE LinearForms

  IMPLICIT NONE

  ! Per-element bubble history, used by Default1stOrderTime's Nb path
  ! (DefUtils.F90) -- lives on Solver % Variable's own BubbleValues/
  ! BubblePrevValues (Types.F90), not a separate type; see the matching
  ! bx/bxprev comment in KomegaLegacy.F90, which this mirrors.

  ! Per-thread ValueHandle_t storage for LocalMatrixVec's material lookups.
  ! NOT THREADPRIVATE -- see the matching comment on IncompressibleNS.F90's
  ! NSHandles_t for the Windows/GCC emutls hazard that rules that out.
  TYPE :: KOHandles_t
    TYPE(ValueHandle_t) :: Visc_h, Dens_h
  END TYPE KOHandles_t
  TYPE(KOHandles_t), ALLOCATABLE, SAVE :: KOHandles(:)

CONTAINS

!------------------------------------------------------------------------------
!> Assemble and glue local matrix/RHS for one bulk element. Vectorized over
!> Gauss points, safe to call concurrently from multiple threads (each thread
!> passes its own InitHandles and only touches KOHandles(tid)).
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, dt, Transient, GlobalBubbles, Stabilize, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
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
    REAL(KIND=dp), ALLOCATABLE :: MassK(:,:), StiffK(:,:), ForceK(:), &
        MassO(:,:), StiffO(:,:), ForceO(:)

    REAL(KIND=dp), POINTER :: RhoVec(:), MuVec(:)

    REAL(KIND=dp), ALLOCATABLE :: VeloNodal(:,:), KNodal(:), ONodal(:)
    REAL(KIND=dp), ALLOCATABLE :: VeloVec(:,:), dVelodxVec(:,:,:), KVec(:), OVec(:), &
        StrainVec(:,:,:), VorticityVec(:,:,:), StrainMeasureVec(:), &
        VorticityMeasureVec(:), TmuVec(:), Effmu1Vec(:), Effmu2Vec(:), &
        ProdVec(:), ReactK(:), ReactO(:), LoadK(:), LoadO(:), &
        StreamVec(:,:), TauK(:), TauO(:), TmpVec(:), RadiusVec(:)

    REAL(KIND=dp) :: Beta,SigmaK,SigmaO,rGamma,hK,mK,VNorm
    INTEGER :: i,j,k,p,ngp,dim,allocstat,tid,ntot
    LOGICAL :: Stat, Found, IsAxiSymmetric
!------------------------------------------------------------------------------
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Visc_h => KOHandles(tid) % Visc_h, Dens_h => KOHandles(tid) % Dens_h )

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity' )
      CALL ListInitElementKeyword( Dens_h,'Material','Density' )
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ! nd is the RETAINED dof count (GetElementNOFDOFs()'s own meaning); ntot
    ! is the bubble-augmented total the test/trial basis actually spans --
    ! see the matching comment in Spalart-Allmaras.F90's own LocalMatrixVec.
    ntot = nd + nb

    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ALLOCATE( BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        MassK(ntot,ntot), StiffK(ntot,ntot), ForceK(ntot), &
        MassO(ntot,ntot), StiffO(ntot,ntot), ForceO(ntot), &
        MASS(2*ntot,2*ntot), STIFF(2*ntot,2*ntot), FORCE(2*ntot), TimeForce(2*ntot), &
        VeloNodal(3,n), KNodal(n), ONodal(n), &
        VeloVec(ngp,3), dVelodxVec(ngp,3,3), KVec(ngp), OVec(ngp), &
        StrainVec(ngp,3,3), VorticityVec(ngp,3,3), StrainMeasureVec(ngp), &
        VorticityMeasureVec(ngp), TmuVec(ngp), Effmu1Vec(ngp), Effmu2Vec(ngp), &
        ProdVec(ngp), ReactK(ngp), ReactO(ngp), LoadK(ngp), LoadO(ngp), &
        StreamVec(ngp,ntot), TauK(ngp), TauO(ngp), TmpVec(ngp), &
        RadiusVec(ngp), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('KOmega','Local storage allocation failed')

    CALL GetElementNodesVec( Nodes, UElement=Element )

    MassK = 0._dp; StiffK = 0._dp; ForceK = 0._dp
    MassO = 0._dp; StiffO = 0._dp; ForceO = 0._dp

    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, DetJVec, &
        SIZE(BasisVec,2), BasisVec, dBasisdxVec )
    DetJVec(1:ngp) = DetJVec(1:ngp) * IP % s(1:ngp)

    ! Axisymmetric (no swirl): r-weighted measure, plus a hoop strain
    ! correction further down -- see the matching (more detailed) comments in
    ! Spalart-Allmaras.F90's own LocalMatrixVec. Genuine swirl ("Cylindric
    ! Symmetric") still goes through LocalMatrixScalar.
    IsAxiSymmetric = ( CurrentCoordinateSystem() == AxisSymmetric )
    IF( IsAxiSymmetric ) THEN
      RadiusVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), Nodes % x(1:n) )
      DetJVec(1:ngp) = DetJVec(1:ngp) * RadiusVec(1:ngp)
    END IF

    ! Nodal input fields, at the n element corner nodes -- as in
    ! Spalart-Allmaras.F90's own LocalMatrixVec, a p-bubble mode never
    ! participates in the closure-law evaluation, only in the test/trial basis.
    VeloNodal = 0._dp
    CALL GetScalarLocalSolution( VeloNodal(1,1:n), 'Velocity 1', UElement=Element )
    CALL GetScalarLocalSolution( VeloNodal(2,1:n), 'Velocity 2', UElement=Element )
    IF( dim == 3 ) CALL GetScalarLocalSolution( VeloNodal(3,1:n), 'Velocity 3', UElement=Element )

    CALL GetScalarLocalSolution( KNodal, 'Kinetic energy', UElement=Element )
    CALL GetScalarLocalSolution( ONodal, 'Kinetic Dissipation', UElement=Element )

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

    KVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), KNodal(1:n) )
    OVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), ONodal(1:n) )

    ! Strain-rate / vorticity tensors and their (Frobenius) measures, Cartesian
    ! only -- mirrors the legacy CurrentCoordinateSystem()==Cartesian branch.
    ! StrainMeasure is floored at 1e-10 here, exactly as the legacy LocalMatrix
    ! does (unlike Spalart-Allmaras, which does not floor it).
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
    ! extra diagonal strain component (covariant, not an ordinary partial
    ! derivative -- see SecondInvariant's own dedicated AxisSymmetric branch
    ! in MaterialModels.F90), entering both the strain measure and the
    ! production term below exactly like any other diagonal Strain(i,i)**2
    ! term. Vorticity needs no such addition (zero diagonal by antisymmetry,
    ! no swirl velocity to give it an off-diagonal r/z-theta component).
    IF( IsAxiSymmetric ) THEN
      StrainMeasureVec(1:ngp) = StrainMeasureVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF

    StrainMeasureVec(1:ngp)    = MAX( SQRT( 2._dp*StrainMeasureVec(1:ngp) ), 1.0d-10 )
    VorticityMeasureVec(1:ngp) = SQRT( 2._dp*VorticityMeasureVec(1:ngp) )

    ! Production, exactly as legacy's Prod = 2*Tmu*SUM(Strain*dVelodx) (the
    ! full, unsymmetrized velocity gradient dotted with the strain tensor --
    ! kept literal rather than simplified to Tmu*StrainMeasure**2, which it
    ! is mathematically equal to, since Strain:Vorticity vanishes).
    Beta = 0.075_dp; SigmaK = 2.000_dp; SigmaO = 2.000_dp; rGamma = 5._dp/9._dp

    TmuVec(1:ngp) = RhoVec(1:ngp)*KVec(1:ngp) / MAX(OVec(1:ngp),1.0d-10)
    Effmu1Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaK
    Effmu2Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaO

    ProdVec = 0._dp
    DO i=1,dim
      DO j=1,dim
        ProdVec(1:ngp) = ProdVec(1:ngp) + StrainVec(1:ngp,i,j)*dVelodxVec(1:ngp,i,j)
      END DO
    END DO
    IF( IsAxiSymmetric ) THEN
      ProdVec(1:ngp) = ProdVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF
    ProdVec(1:ngp) = 2._dp*TmuVec(1:ngp)*ProdVec(1:ngp)

    ! K-equation: destruction rho*0.09*Omega, diffusion Effmu1, production Prod.
    ReactK(1:ngp) = RhoVec(1:ngp) * 0.09_dp * OVec(1:ngp)
    LoadK(1:ngp)  = ProdVec(1:ngp)

    ! Omega-equation: destruction rho*Beta*Omega, diffusion Effmu2,
    ! production rGamma*Prod*Omega/K (no floor on K here, exactly as legacy).
    ReactO(1:ngp) = RhoVec(1:ngp) * Beta * OVec(1:ngp)
    LoadO(1:ngp)  = rGamma * ProdVec(1:ngp) * OVec(1:ngp) / KVec(1:ngp)

    ! Both equations share the same convection velocity (no extra grad term,
    ! unlike Spalart-Allmaras's Cb2 correction).
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, MassK, RhoVec )
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffK, ReactK )
    CALL LinearForms_GradUdotGradU( ngp, ntot, dim, dBasisdxVec, DetJVec, StiffK, Effmu1Vec )
    CALL LinearForms_GradUdotU( ngp, ntot, dim, dBasisdxVec, BasisVec, DetJVec, StiffK, &
        RhoVec, VeloVec )
    CALL LinearForms_UdotF( ngp, ntot, BasisVec, DetJVec, LoadK, ForceK )

    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, MassO, RhoVec )
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffO, ReactO )
    CALL LinearForms_GradUdotGradU( ngp, ntot, dim, dBasisdxVec, DetJVec, StiffO, Effmu2Vec )
    CALL LinearForms_GradUdotU( ngp, ntot, dim, dBasisdxVec, BasisVec, DetJVec, StiffO, &
        RhoVec, VeloVec )
    CALL LinearForms_UdotF( ngp, ntot, BasisVec, DetJVec, LoadO, ForceO )

    !------------------------------------------------------------------------
    ! SUPG (equal-order) stabilization, opt-in via "Stabilize"/"Stabilization
    ! Method" -- same Franca et al. tau/streamline construction as
    ! Spalart-Allmaras.F90's own SUPG block. The streamline-weighted test
    ! function StreamVec is shared (both equations convect with the same
    ! Velo), only Tau differs (it depends on each equation's own Effmu).
    !------------------------------------------------------------------------
    IF( Stabilize ) THEN
      hK = Element % hK
      mK = Element % StabilizationMK

      StreamVec(1:ngp,1:ntot) = 0._dp
      DO i=1,dim
        DO p=1,ntot
          StreamVec(1:ngp,p) = StreamVec(1:ngp,p) + &
              RhoVec(1:ngp) * VeloVec(1:ngp,i) * dBasisdxVec(1:ngp,p,i)
        END DO
      END DO

      DO j=1,ngp
        VNorm = SQRT( SUM( VeloVec(j,1:dim)**2 ) )
        IF( VNorm > 0._dp .AND. Effmu1Vec(j) /= 0._dp ) THEN
          TauK(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(Effmu1Vec(j))) )
          TauK(j) = hK * TauK(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauK(j) = 0._dp
        END IF
        IF( VNorm > 0._dp .AND. Effmu2Vec(j) /= 0._dp ) THEN
          TauO(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(Effmu2Vec(j))) )
          TauO(j) = hK * TauO(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauO(j) = 0._dp
        END IF
      END DO

      TmpVec(1:ngp) = TauK(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, StiffK )
      CALL LinearForms_UdotF( ngp, ntot, StreamVec, TmpVec, LoadK, ForceK )

      TmpVec(1:ngp) = TauO(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, StiffO )
      CALL LinearForms_UdotF( ngp, ntot, StreamVec, TmpVec, LoadO, ForceO )

      IF( Transient ) THEN
        TmpVec(1:ngp) = TauK(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MassK )
        TmpVec(1:ngp) = TauO(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MassO )
      END IF
    END IF

    ! Interleave the two independent scalar blocks into the final K/omega
    ! system, odd rows/columns for K, even for omega -- the off-diagonal
    ! (K,omega) blocks stay zero, matching legacy's own A(1,2)/A(2,1) (never
    ! written).
    MASS = 0._dp; STIFF = 0._dp; FORCE = 0._dp
    MASS(1:2*ntot-1:2,1:2*ntot-1:2) = MassK
    MASS(2:2*ntot:2,  2:2*ntot:2)   = MassO
    STIFF(1:2*ntot-1:2,1:2*ntot-1:2) = StiffK
    STIFF(2:2*ntot:2,  2:2*ntot:2)   = StiffO
    FORCE(1:2*ntot-1:2) = ForceK
    FORCE(2:2*ntot:2)   = ForceO

    !------------------------------------------------------------------------
    ! Time discretization and p-bubble condensation -- mirrors the nb>0
    ! branch of KomegaLegacy.F90's driver exactly (DOFs=2, K/omega
    ! interleaved; no legacy "Bubbles" per-node scheme here, see the file
    ! header).
    !------------------------------------------------------------------------
    TimeForce = 0._dp
    IF( nb > 0 ) THEN
      IF( Transient .AND. .NOT. GlobalBubbles ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=nb )
      ELSE
        IF( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
        CALL CondensateP( 2*nd, 2*nb, STIFF, FORCE, TimeForce )
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
!> the same metric-tensor branch as KomegaLegacy.F90's own LocalMatrix --
!> verbatim, no new math, and already interleaves K/omega directly via
!> STIFF(2*(p-1)+i,2*(q-1)+j)), and the legacy "Bubbles = True" per-node
!> scheme for a plain nodal "Element" set explicitly. Always called serially
!> (see the AxiSymmetric branch in KOmega below), so it uses the classic
!> GetMaterial()/GetReal()/argument-less GetElementNOF*() accessors exactly
!> as the legacy driver does -- not safe to call from inside an OMP parallel
!> region, unlike LocalMatrixVec.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixScalar( Element, dt, Transient, GlobalBubbles, BubblesDefault )
!------------------------------------------------------------------------------
    IMPLICIT NONE
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

    ALLOCATE( MASS(2*(nd+nb),2*(nd+nb)), STIFF(2*(nd+nb),2*(nd+nb)), &
        FORCE(2*(nd+nb)), LOAD(2,n), TimeForce(2*(nd+nb)), &
        STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('KOmega','Local storage allocation failed')

    CALL ElementKernel( MASS, STIFF, FORCE, LOAD, Element, n, nd+nb, ElementNodes )

    TimeForce = 0.0_dp
    IF ( Bubbles ) THEN
      IF ( Transient ) THEN
        ! Nb=n: "as many bubbles as nodes", same convention as the nb>0
        ! branch below. Never gated on Solver % GlobalBubbles -- legacy
        ! bubbles are always locally condensed.
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=n )
      ELSE
        CALL Condensate( 2*n, STIFF, FORCE, TimeForce )
      END IF
    ELSE IF ( nb > 0 ) THEN
      IF ( Transient .AND. .NOT. GlobalBubbles ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element, &
            Nb=nb )
      ELSE
        IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
        CALL CondensateP( 2*nd, 2*nb, STIFF, FORCE, TimeForce )
      END IF
    ELSE
      IF ( Transient ) CALL Default1stOrderTime( MASS, STIFF, FORCE, UElement=Element )
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )

  CONTAINS

!------------------------------------------------------------------------------
!> Verbatim port of KomegaLegacy.F90's nested LocalMatrix: same per-Gauss-
!> point scalar math, same metric-tensor (axisymmetric/cylindrical) branch,
!> same closure-law formulas, same direct K/omega interleaving. "Bubbles" and
!> "Material" come from the host (LocalMatrixScalar above), exactly as they
!> did from the legacy driver's own host scope.
!------------------------------------------------------------------------------
    SUBROUTINE ElementKernel( MASS,STIFF,FORCE, LOAD, Element,n,nd,Nodes )
!------------------------------------------------------------------------------
      USE MaterialModels

      IMPLICIT NONE

      REAL(KIND=dp), DIMENSION(:)   :: FORCE
      REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

      INTEGER :: n, nd

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: ddBasisddx(nd,3,3)
      REAL(KIND=dp) :: Basis(nd)
      REAL(KIND=dp) :: dBasisdx(nd,3),detJ

      REAL(KIND=dp) :: UX(n), UY(n), UZ(n), Velo(3), dVelodx(3,3), Energy(n), &
                       Dissipation(n), Distance(n), Density(n), Viscosity(n)

      REAL(KIND=dp) :: A(2,2),M(2,2),Prod,div,ProdTensor(3,3)
      INTEGER :: i,j,c,p,q,t,dim,NBasis
      REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2)

      REAL(KIND=dp) :: s,u,v,w, K,Omega,Strain(3,3), Vorticity(3,3), dist

      REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,SigmaK, &
              SigmaO,Beta,CD,F1,F2,F3,F4,rGamma, GradK(3), GradO(3)

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

      CALL GetScalarLocalSolution( Energy, 'Kinetic energy' )
      CALL GetScalarLocalSolution( Dissipation, 'Kinetic Dissipation' )

      FORCE = 0.0D0
      STIFF = 0.0D0
      MASS  = 0.0D0

      NBasis = nd

      IF ( Bubbles ) THEN
         IntegStuff = GaussPoints( element, element % TYPE % GaussPoints2 )
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
          X = SUM( Nodes % x(1:n)*Basis(1:n) )
          Y = SUM( Nodes % y(1:n)*Basis(1:n) )
          Z = SUM( nodes % z(1:n)*Basis(1:n) )
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
          StrainMeasure = MAX( SQRT(2 * SUM(Strain * Strain)), 1.0d-10 )

          Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
          VorticityMeasure = SQRT(2 * SUM(Vorticity * Vorticity))
        ELSE
          StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
        END IF

        K = SUM( Energy(1:n) * Basis(1:n) )
        Omega = SUM( Dissipation(1:n) * Basis(1:n) )

        mu   = SUM( Viscosity(1:n) * Basis(1:n) )
        rho  = SUM( Density(1:n) * Basis(1:n) )

        Beta   = 0.075_dp
        SigmaK = 2.000_dp
        SigmaO = 2.000_dp
        rGamma = 5._dp/9._dp

        Tmu = rho*K/MAX(Omega,1.0d-10)
        Effmu(1) = mu + Tmu / SigmaK
        Effmu(2) = mu + Tmu / SigmaO

        Prod = 2*Tmu*SUM(Strain*dVelodx)

        DO p=1,NBasis
        DO q=1,NBasis
           M = 0.0d0
           A = 0.0d0

           M(1,1) = rho * Basis(q) * Basis(p)
           M(2,2) = rho * Basis(q) * Basis(p)

           A(1,1) = A(1,1) + rho * 0.09_dp * Omega * Basis(q) * Basis(p)
           A(2,2) = A(2,2) + rho * Beta * Omega * Basis(q) * Basis(p)

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO i=1,dim
                A(1,1) = A(1,1) + Effmu(1) * dBasisdx(q,i) * dBasisdx(p,i)
                A(2,2) = A(2,2) + Effmu(2) * dBasisdx(q,i) * dBasisdx(p,i)
              END DO
           ELSE
              DO i=1,dim
                DO j=1,dim
                   A(1,1) = A(1,1) + Metric(i,j) * Effmu(1) * &
                        dBasisdx(q,i) * dBasisdx(p,j)

                   A(2,2) = A(2,2) + Metric(i,j) * Effmu(2) * &
                        dBasisdx(q,i) * dBasisdx(p,j)
                END DO
              END DO
           END IF

           DO i=1,dim
             A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
             A(2,2) = A(2,2) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO

           DO i=1,2
              DO j=1,2
                STIFF(2*(p-1)+i,2*(q-1)+j) = STIFF(2*(p-1)+i,2*(q-1)+j)+s*A(i,j)
                MASS(2*(p-1)+i,2*(q-1)+j)  = MASS(2*(p-1)+i,2*(q-1)+j) +s*M(i,j)
              END DO
           END DO
         END DO
         END DO

         LoadAtIP(1) = Prod
         LoadAtIP(2) = rGamma * Prod * Omega / K

         DO p=1,NBasis
           FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1)+s*LoadAtIp(1)*Basis(p)
           FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2)+s*LoadAtIp(2)*Basis(p)
         END DO
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE ElementKernel
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixScalar
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Wall law for the k-omega turbulence model, ported unchanged from
!> KomegaLegacy.F90's own nested OmegaWall bar taking Solver explicitly (it is
!> a module procedure here, not nested inside the driver, so it no longer
!> gets Solver/Model by host association) and CurrentModel % Nodes in place
!> of Model % Nodes.
!------------------------------------------------------------------------------
  SUBROUTINE OmegaWall( Element, n, Solver )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: omega_wall,dist,mu(32),rho(32),x0(n),y0(n),z0(n),x,y,z
    INTEGER :: i,j,np
    TYPE(Element_t), POINTER :: Parent
!------------------------------------------------------------------------------
    Parent => Element % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) &
      Parent => Element % BoundaryInfo % Right
    IF ( .NOT. ASSOCIATED(Parent) ) RETURN

    np = GetElementNOFNodes(Parent)

    rho(1:np)= GetReal( GetMaterial(Parent), 'Density', UElement=Parent )
    mu(1:np) = GetReal( GetMaterial(Parent), 'Viscosity', UElement=Parent )

    x0(1:n) = CurrentModel % Nodes % x(Element % NodeIndexes)
    y0(1:n) = CurrentModel % Nodes % y(Element % NodeIndexes)
    z0(1:n) = CurrentModel % Nodes % z(Element % NodeIndexes)

    omega_wall = 1.d10
    DO i=1,np
      j = Parent % NodeIndexes(i)
      IF ( ANY( j==Element % NodeIndexes(1:n) ) ) CYCLE

      x = CurrentModel % Nodes % x(j)
      y = CurrentModel % Nodes % y(j)
      z = CurrentModel % Nodes % z(j)

      dist = MINVAL( (x-x0(1:n))**2 + (y-y0(1:n))**2 + (z-z0(1:n))**2 )
      IF ( dist < AEPS ) CYCLE

      omega_wall = 6*mu(i)/rho(i)/0.075_dp/dist

      j = 2*Solver % Variable % Perm(j)

      CALL UpdateDirichletDof( Solver % Matrix, j, omega_wall )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE OmegaWall
!------------------------------------------------------------------------------

END MODULE KOmegaLocalForms


!------------------------------------------------------------------------------
!> Vectorized/threaded K-omega driver. See the file header above for what
!> this does and does not support, and "Legacy Assembly" for the fallback.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KOmega( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KOmegaLocalForms
  USE KOmegaFront
  IMPLICIT NONE
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
  REAL(KIND=dp) :: Norm, KVal, EVal, KMax, EMax
  CHARACTER(*), PARAMETER :: Caller = 'KOmega'
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToKOmegaLegacy( 'KOmegaLegacy', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  KE => Solver % Variable
  IF ( .NOT. ASSOCIATED(KE) ) RETURN
  IF ( COUNT( KE % Perm > 0 ) <= 0 ) RETURN

  ! LocalMatrixVec now carries the plain "Axi Symmetric" (no swirl) case
  ! itself; genuine swirl ("Cylindric Symmetric") and general "Cylindric"
  ! still need LocalMatrixScalar's full metric/Christoffel treatment, so
  ! those still go through it, serially -- same branch HeatSolve.F90 makes,
  ! and the same treatment Spalart-Allmaras.F90 now has.
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian .AND. &
                   CurrentCoordinateSystem() /= AxisSymmetric )

  IF (.NOT. ALLOCATED(KOHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(KOHandles(nthr))
  END IF

  GlobalBubbles = Solver % GlobalBubbles
  Stabilize = GetStabilizeFlag( Solver % Values, GotIt )

  ! Only LocalMatrixScalar's legacy per-node branch uses this -- see the
  ! matching BubblesDefault resolution in KomegaLegacy.F90.
  BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  IF ( .NOT.GotIt ) BubblesDefault = .TRUE.

  ! K and Omega are interleaved, Dofs=2 -- exactly as KomegaLegacy.F90 sizes
  ! its own bx/bxprev.
  IF ( TransientSimulation ) CALL DefaultBubbleHistoryUpdate( Dofs=2 )

  NonlinearIter = ListGetInteger( Solver % Values, 'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT.GotIt ) NonlinearIter = 1

  DO i=1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    IF ( ListGetLogical( BC, 'Noslip wall BC', gotit ) ) THEN
      CALL ListAddConstReal( BC, 'Kinetic Energy', 0.0_dp )
    END IF
  END DO

  DO iter=1,NonlinearIter
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)
    CALL Info(Caller,'Komega iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)

    CALL DefaultInitialize()

    Active = GetNOFActive()

    IF ( AxiSymmetric ) THEN
      ! Serial only: LocalMatrixScalar relies on the classic argument-less
      ! Get*() "current element" accessors, exactly as KomegaLegacy.F90's
      ! own driver does.
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

    CALL DefaultFinishBulkAssembly()

    DO i=1,Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(i)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE
      n = GetElementNOFNodes()
      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE
      IF (ListGetLogical(BC, 'Omega Wall BC', gotIt ) .OR. &
          ListGetLogical(BC, 'Noslip Wall BC',  gotIt)) CALL OmegaWall(Element,n,Solver)
    END DO

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    ! K and omega should stay positive.
    DO i=1,SIZE(Solver % Variable % Perm)
      k = Solver % Variable % Perm(i)
      IF ( k <= 0 ) CYCLE
      KVal = Solver % Variable % Values(2*k-1)
      EVal = Solver % Variable % Values(2*k-0)
      Solver % Variable % Values(2*k-1) = MAX( KVal, 1.0d-12 )
      Solver % Variable % Values(2*k-0) = MAX( EVal, 1.0d-12 )
    END DO

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE KOmega
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Element-type resolution, called before the mesh/basis functions are
!> finalized -- same role and timing as HeatSolver_Init0 and
!> Spalart-Allmaras.F90's own SpalartAllmaras_Init0.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KOmega_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KOmegaFront
  IMPLICIT NONE
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

  Stabilize = GetStabilizeFlag( Params )

  IF( .NOT. ListCheckPresent( Params,'Element' ) ) THEN
    IF( Stabilize ) THEN
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
END SUBROUTINE KOmega_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: KOmega. Under "Legacy Assembly",
!> delegates to KOmegaLegacy_Init and returns -- the p-bubble setup below is
!> specific to this file's own implementation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KOmega_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KOmegaFront
  IMPLICIT NONE
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
    CALL DelegateToKOmegaLegacy( 'KOmegaLegacy_Init', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  SolverParams => GetSolverParams()

  IF ( ListGetLogical( SolverParams, 'Bubbles', Found ) ) THEN
    CALL Warn('KOmega_Init', &
        '"Bubbles = True" (the legacy per-node scheme) is not used here -- '// &
        'KOmega_Init0 already defaulted "Element" to an equivalent p-element '// &
        'bubble unless SUPG ("Stabilize"/"Stabilization Method") was requested '// &
        'instead. Use "Legacy Assembly = Logical True" for the original per-node '// &
        'scheme.')
  END IF

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
END SUBROUTINE KOmega_Init
!------------------------------------------------------------------------------
