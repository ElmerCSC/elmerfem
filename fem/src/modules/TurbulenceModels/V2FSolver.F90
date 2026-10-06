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
! *  Module containing a solver for the V2-F-turbulence model.
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 08 Jun 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solve the V2-F (LDM) turbulence model equations. Vectorized/threaded
!> implementation following the same pattern as KESolver.F90/Komega.F90/
!> SSTKomega.F90/Spalart-Allmaras.F90.
!>
!> Unlike any of those, V2 and F are genuinely coupled BOTH ways: the legacy
!> LocalMatrix writes A(1,2) = -Rho*K*Basis(q)*Basis(p) (V2-row/F-column) and
!> A(2,1) = (C1-6)/K/TimeScale*Basis(q)*Basis(p) (F-row/V2-column). This
!> still needs no block-coupling machinery beyond LinearForms: each
!> off-diagonal block is just one more LinearForms_UdotU call (the same
!> "alpha=coefficient" weighted mass-type form the diagonal reaction terms
!> already use), interleaved into STIFF(1:2*ntot-1:2, 2:2*ntot:2) and
!> STIFF(2:2*ntot:2, 1:2*ntot-1:2) instead of a diagonal slot. F itself has
!> no mass/time-derivative term at all (legacy's own MASS(2,2) is never
!> written) and no convection (F is a purely elliptic relaxation equation),
!> so its own diagonal block skips LinearForms_UdotU/GradUdotU for those.
!>
!> Otherwise simpler than KESolver/SSTKomega: a single closure ("V2-F Model"
!> is read but never branched on in the legacy LocalMatrix either -- there
!> is only the LDM formulation), no buoyancy/Mach-number terms, and Density
!> is a plain material property (never ElementDensity), so there is no
!> "Compressibility Model" case to route elsewhere.
!>
!> Axisymmetric/cylindrical coordinates go through LocalMatrixScalar instead
!> of LocalMatrixVec -- a scalar, per-Gauss-point fallback carrying the same
!> metric-tensor math as V2FSolverLegacy.F90's own LocalMatrix (already
!> interleaved directly via STIFF(2*(p-1)+i,2*(q-1)+j)), called serially --
!> same role as HeatSolve.F90's own AxiSymmetric branch. LocalMatrixScalar
!> also carries the legacy "Bubbles = True" per-node scheme for a plain
!> nodal "Element" set explicitly. There is no wall-function/weak-form
!> boundary treatment here (unlike KESolver/SSTKomega): V2/F just get a
!> plain zero Dirichlet value on "Noslip wall BC", exactly as legacy.
!>
!> The original scalar-element solver lives on in V2FSolverLegacy.F90
!> (subroutine V2F_LDM_Legacy), reachable either directly by that name or
!> via "Legacy Assembly = Logical True" here (see V2FSolverFront below).
!> With axisymmetric and the per-node scheme both covered above, that path
!> is now only needed by a sif that must reproduce the legacy solver's exact
!> historical numbers.
!> \ingroup Solvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Whether this solver should run the original scalar-element assembly
!> (V2FSolverLegacy.F90) instead of this file's own implementation.
!------------------------------------------------------------------------------
MODULE V2FSolverFront
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
!> Call one of V2F_LDM_Legacy's entry points with this solver. The name is
!> resolved at run time, as the core resolves any solver, so this file and
!> V2FSolverLegacy.so stay independent of one another.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateToV2FSolverLegacy( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'V2FSolverLegacy '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'V2F_LDM', &
        '"Legacy Assembly" was requested but "'//TRIM(Entry)//'" could not be found. '// &
        'Is V2FSolverLegacy.so installed beside this solver?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateToV2FSolverLegacy

END MODULE V2FSolverFront


!------------------------------------------------------------------------------
MODULE V2FSolverLocalForms

  USE DefUtils
  USE LinearForms

  IMPLICIT NONE IMPLICIT_EXTERNAL

  ! Per-element bubble history, used by Default1stOrderTime's Nb path
  ! (DefUtils.F90) -- lives on Solver % Variable's own BubbleValues/
  ! BubblePrevValues (Types.F90), not a separate type; see the matching
  ! bx/bxprev comment in V2FSolverLegacy.F90, which this mirrors.

  ! Per-thread ValueHandle_t storage for LocalMatrixVec's material lookups.
  ! NOT THREADPRIVATE -- see the matching comment on IncompressibleNS.F90's
  ! NSHandles_t for the Windows/GCC emutls hazard that rules that out.
  TYPE :: V2FHandles_t
    TYPE(ValueHandle_t) :: Visc_h, Dens_h, Sigma_h, Cmu_h, C1_h, C2_h, CT_h, CL_h, Cnu_h
  END TYPE V2FHandles_t
  TYPE(V2FHandles_t), ALLOCATABLE, SAVE :: V2FHandles(:)

CONTAINS

!------------------------------------------------------------------------------
!> Assemble and glue local matrix/RHS for one bulk element. Vectorized over
!> Gauss points, safe to call concurrently from multiple threads (each thread
!> passes its own InitHandles and only touches V2FHandles(tid)).
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
    REAL(KIND=dp), ALLOCATABLE :: MassV2(:,:), StiffV2(:,:), StiffV2F(:,:), &
        StiffFV2(:,:), StiffFF(:,:), ForceF(:)

    REAL(KIND=dp), POINTER :: RhoVec(:), MuVec(:), SigmaVec(:), CmuVec(:), &
        C1Vec(:), C2Vec(:), CTVec(:), CLVec(:), CnuVec(:)

    REAL(KIND=dp), ALLOCATABLE :: VeloNodal(:,:), KNodal(:), ONodal(:), V2Nodal(:), FNodal(:)
    REAL(KIND=dp), ALLOCATABLE :: VeloVec(:,:), dVelodxVec(:,:,:), KVec(:), OVec(:), &
        V2Vec(:), FVec(:), StrainVec(:,:,:), SecInvVec(:), TimeScaleVec(:), &
        LengthScale2Vec(:), TmuVec(:), EffViscVec(:), ProdVec(:), &
        ReactV2(:), CrossV2F(:), CrossFV2(:), LoadF(:), &
        StreamVec(:,:), TauVec(:), TmpVec(:), RadiusVec(:)

    REAL(KIND=dp) :: hK, mK, VNorm
    INTEGER :: i,j,k,p,ngp,dim,allocstat,tid,ntot
    LOGICAL :: Stat, Found, IsAxiSymmetric
!------------------------------------------------------------------------------
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Visc_h => V2FHandles(tid) % Visc_h, Dens_h => V2FHandles(tid) % Dens_h, &
               Sigma_h => V2FHandles(tid) % Sigma_h, Cmu_h => V2FHandles(tid) % Cmu_h, &
               C1_h => V2FHandles(tid) % C1_h, C2_h => V2FHandles(tid) % C2_h, &
               CT_h => V2FHandles(tid) % CT_h, CL_h => V2FHandles(tid) % CL_h, &
               Cnu_h => V2FHandles(tid) % Cnu_h )

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity' )
      CALL ListInitElementKeyword( Dens_h,'Material','Density' )
      CALL ListInitElementKeyword( Sigma_h,'Material','V2-F Sigma' )
      CALL ListInitElementKeyword( Cmu_h,'Material','KE Cmu' )
      CALL ListInitElementKeyword( C1_h,'Material','V2-F C1' )
      CALL ListInitElementKeyword( C2_h,'Material','V2-F C2' )
      CALL ListInitElementKeyword( CT_h,'Material','V2-F CT' )
      CALL ListInitElementKeyword( CL_h,'Material','V2-F CL' )
      CALL ListInitElementKeyword( Cnu_h,'Material','V2-F Cnu' )
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ntot = nd + nb

    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ALLOCATE( BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        MassV2(ntot,ntot), StiffV2(ntot,ntot), StiffV2F(ntot,ntot), &
        StiffFV2(ntot,ntot), StiffFF(ntot,ntot), ForceF(ntot), &
        MASS(2*ntot,2*ntot), STIFF(2*ntot,2*ntot), FORCE(2*ntot), TimeForce(2*ntot), &
        VeloNodal(3,n), KNodal(n), ONodal(n), V2Nodal(n), FNodal(n), &
        VeloVec(ngp,3), dVelodxVec(ngp,3,3), KVec(ngp), OVec(ngp), &
        V2Vec(ngp), FVec(ngp), StrainVec(ngp,3,3), SecInvVec(ngp), TimeScaleVec(ngp), &
        LengthScale2Vec(ngp), TmuVec(ngp), EffViscVec(ngp), ProdVec(ngp), &
        ReactV2(ngp), CrossV2F(ngp), CrossFV2(ngp), LoadF(ngp), &
        StreamVec(ngp,ntot), TauVec(ngp), TmpVec(ngp), &
        RadiusVec(ngp), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('V2FSolver','Local storage allocation failed')

    CALL GetElementNodesVec( Nodes, UElement=Element )

    MassV2 = 0._dp; StiffV2 = 0._dp
    StiffV2F = 0._dp; StiffFV2 = 0._dp; StiffFF = 0._dp; ForceF = 0._dp

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

    VeloNodal = 0._dp
    CALL GetScalarLocalSolution( VeloNodal(1,1:n), 'Velocity 1', UElement=Element )
    CALL GetScalarLocalSolution( VeloNodal(2,1:n), 'Velocity 2', UElement=Element )
    IF( dim == 3 ) CALL GetScalarLocalSolution( VeloNodal(3,1:n), 'Velocity 3', UElement=Element )

    CALL GetScalarLocalSolution( KNodal, 'Kinetic Energy', UElement=Element )
    CALL GetScalarLocalSolution( ONodal, 'Kinetic Dissipation', UElement=Element )
    CALL GetScalarLocalSolution( V2Nodal, 'V2', UElement=Element )
    CALL GetScalarLocalSolution( FNodal, 'F', UElement=Element )

    RhoVec => ListGetElementRealVec( Dens_h, ngp, BasisVec, Element, Found )
    MuVec  => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element, Found )

    SigmaVec => ListGetElementRealVec( Sigma_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) SigmaVec = 1.0_dp

    CmuVec => ListGetElementRealVec( Cmu_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) CmuVec = 0.22_dp

    C1Vec => ListGetElementRealVec( C1_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) C1Vec = 1.4_dp

    C2Vec => ListGetElementRealVec( C2_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) C2Vec = 0.3_dp

    CTVec => ListGetElementRealVec( CT_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) CTVec = 6.0_dp

    CLVec => ListGetElementRealVec( CL_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) CLVec = 0.23_dp

    CnuVec => ListGetElementRealVec( Cnu_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) CnuVec = 70.0_dp

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

    KVec(1:ngp) = MAX( MATMUL( BasisVec(1:ngp,1:n), KNodal(1:n) ), 1.0d-10 )
    OVec(1:ngp) = MAX( MATMUL( BasisVec(1:ngp,1:n), ONodal(1:n) ), 1.0d-10 )
    V2Vec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), V2Nodal(1:n) )
    FVec(1:ngp)  = MATMUL( BasisVec(1:ngp,1:n), FNodal(1:n) )

    ! Strain tensor and SecInv = 2*Strain:Strain, Cartesian only. Not floored
    ! -- legacy's own LocalMatrix does not floor SecInv either.
    StrainVec = 0._dp
    DO i=1,dim
      DO k=1,dim
        StrainVec(1:ngp,i,k) = 0.5_dp*( dVelodxVec(1:ngp,i,k) + dVelodxVec(1:ngp,k,i) )
      END DO
    END DO
    SecInvVec = 0._dp
    DO i=1,dim
      DO j=1,dim
        SecInvVec(1:ngp) = SecInvVec(1:ngp) + StrainVec(1:ngp,i,j)**2
      END DO
    END DO

    ! Axisymmetric (no swirl) hoop strain e_theta_theta = u_r/r: a genuine
    ! extra diagonal strain component (covariant, not an ordinary partial
    ! derivative -- see SecondInvariant's own dedicated AxisSymmetric branch
    ! in MaterialModels.F90, and the matching comment in Spalart-Allmaras.F90).
    IF( IsAxiSymmetric ) THEN
      SecInvVec(1:ngp) = SecInvVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF

    SecInvVec(1:ngp) = 2._dp*SecInvVec(1:ngp)

    TimeScaleVec(1:ngp) = MAX( KVec(1:ngp)/OVec(1:ngp), &
        CTVec(1:ngp)*SQRT((MuVec(1:ngp)/RhoVec(1:ngp))/OVec(1:ngp)) )
    LengthScale2Vec(1:ngp) = CLVec(1:ngp)**2 * MAX( KVec(1:ngp)**3/OVec(1:ngp)**2, &
        CnuVec(1:ngp)**2*SQRT((MuVec(1:ngp)/RhoVec(1:ngp))**3/OVec(1:ngp)) )

    TmuVec(1:ngp) = RhoVec(1:ngp)*CmuVec(1:ngp)*V2Vec(1:ngp)*TimeScaleVec(1:ngp)
    EffViscVec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaVec(1:ngp)

    ProdVec(1:ngp) = TmuVec(1:ngp) * SecInvVec(1:ngp) / RhoVec(1:ngp)

    ! V2-equation reaction (diagonal) and the V2/F cross coefficients --
    ! exactly legacy's A(1,1) reaction term, A(1,2) and A(2,1).
    ReactV2(1:ngp)  = 6._dp * RhoVec(1:ngp) / TimeScaleVec(1:ngp)
    CrossV2F(1:ngp) = -RhoVec(1:ngp) * KVec(1:ngp)
    CrossFV2(1:ngp) = ( C1Vec(1:ngp) - 6._dp ) / KVec(1:ngp) / TimeScaleVec(1:ngp)

    ! F-equation load -- legacy's LoadAtIP(2); LoadAtIP(1) (V2's own force)
    ! is always zero, so ForceV2 is never built.
    LoadF(1:ngp) = ( (C1Vec(1:ngp)-1._dp)*(2._dp/3._dp) + C2Vec(1:ngp)*ProdVec(1:ngp)/OVec(1:ngp) ) &
        / TimeScaleVec(1:ngp)

    ! V2's own diagonal block: mass, reaction, diffusion, convection.
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, MassV2, RhoVec )
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffV2, ReactV2 )
    CALL LinearForms_GradUdotGradU( ngp, ntot, dim, dBasisdxVec, DetJVec, StiffV2, EffViscVec )
    CALL LinearForms_GradUdotU( ngp, ntot, dim, dBasisdxVec, BasisVec, DetJVec, StiffV2, &
        RhoVec, VeloVec )

    ! V2/F and F/V2 cross blocks (one LinearForms_UdotU call each -- same
    ! weighted-mass form the diagonal reaction terms use, just placed off
    ! the interleaved diagonal below).
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffV2F, CrossV2F )
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffFV2, CrossFV2 )

    ! F's own diagonal block: reaction coefficient 1 (no alpha needed) plus
    ! diffusion -- no mass, no convection (F has neither in the legacy
    ! LocalMatrix: MASS(2,2) is never written, and A(2,2)'s convection loop
    ! is absent).
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffFF )
    CALL LinearForms_GradUdotGradU( ngp, ntot, dim, dBasisdxVec, DetJVec, StiffFF, LengthScale2Vec )
    CALL LinearForms_UdotF( ngp, ntot, BasisVec, DetJVec, LoadF, ForceF )

    !------------------------------------------------------------------------
    ! SUPG (equal-order) stabilization for V2's own convection-diffusion,
    ! opt-in via "Stabilize"/"Stabilization Method" -- same Franca et al.
    ! tau/streamline construction as KESolver.F90's own SUPG block. F has no
    ! convection, so it is never stabilized. Legacy hard-disables any
    ! stabilization for this solver (Stabilize is set to .FALSE.
    ! unconditionally, the ListGetLogical read commented out) -- this is a
    ! genuinely new, opt-in capability, off by default like everywhere else.
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
        IF( VNorm > 0._dp .AND. EffViscVec(j) /= 0._dp ) THEN
          TauVec(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(EffViscVec(j))) )
          TauVec(j) = hK * TauVec(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauVec(j) = 0._dp
        END IF
      END DO

      TmpVec(1:ngp) = TauVec(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, StiffV2 )

      IF( Transient ) THEN
        TmpVec(1:ngp) = TauVec(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MassV2 )
      END IF
    END IF

    ! Interleave: odd rows/columns for V2, even for F. Both cross blocks are
    ! real here (unlike Komega/SSTKomega/KESolver), matching legacy's own
    ! A(1,2) and A(2,1).
    MASS = 0._dp; STIFF = 0._dp; FORCE = 0._dp
    MASS(1:2*ntot-1:2,1:2*ntot-1:2)  = MassV2
    STIFF(1:2*ntot-1:2,1:2*ntot-1:2) = StiffV2
    STIFF(2:2*ntot:2,  2:2*ntot:2)   = StiffFF
    STIFF(1:2*ntot-1:2,2:2*ntot:2)   = StiffV2F
    STIFF(2:2*ntot:2,  1:2*ntot-1:2) = StiffFV2
    FORCE(2:2*ntot:2) = ForceF

    !------------------------------------------------------------------------
    ! Time discretization and p-bubble condensation -- mirrors the nb>0
    ! branch of V2FSolverLegacy.F90's driver exactly (DOFs=2, V2/F
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
!> the same metric-tensor math as V2FSolverLegacy.F90's own LocalMatrix --
!> verbatim, no new math, and already interleaves V2/F directly via
!> STIFF(2*(p-1)+i,2*(q-1)+j)), and the legacy "Bubbles = True" per-node
!> scheme for a plain nodal "Element" set explicitly. Always called serially
!> (see the AxiSymmetric branch in V2F_LDM below), so it uses the classic
!> GetMaterial()/GetReal()/argument-less GetElementNOF*() accessors exactly
!> as the legacy driver does -- not safe to call from inside an OMP parallel
!> region, unlike LocalMatrixVec.
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
    CHARACTER(LEN=MAX_NAME_LEN) :: V2FModel
    REAL(KIND=dp) :: Clip
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:), &
        TimeForce(:), LocalV2(:), LocalF(:), &
        U(:), V(:), W(:), Density(:), Viscosity(:), &
        KECmu(:), V2FC1(:), V2FC2(:), V2FCL(:), V2FCT(:), V2FCnu(:), V2FSigma(:)
    LOGICAL :: Bubbles, GotIt
    INTEGER :: n, nd, nb, allocstat
!------------------------------------------------------------------------------
    Bubbles = BubblesDefault .AND. .NOT. ASSOCIATED( Element % PDefs )
    Material => GetMaterial()

    Clip = GetConstReal( Material, 'KE Clip', GotIt )
    IF ( .NOT.GotIt ) Clip = 1.0d-6

    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    IF ( Bubbles ) nd = 2*n
    nb = GetElementNOFBDOFs()

    CALL GetElementNodes( ElementNodes )

    ALLOCATE( MASS(2*(nd+nb),2*(nd+nb)), STIFF(2*(nd+nb),2*(nd+nb)), &
        FORCE(2*(nd+nb)), LOAD(2,n), TimeForce(2*(nd+nb)), &
        LocalV2(nd+nb), LocalF(nd+nb), &
        U(n), V(n), W(n), Density(n), Viscosity(n), &
        KECmu(n), V2FC1(n), V2FC2(n), V2FCL(n), V2FCT(n), V2FCnu(n), V2FSigma(n), &
        STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('V2FSolver','Local storage allocation failed')

    V2FModel = GetString( Material, 'V2-F Model', GotIt )
    IF ( .NOT. GotIt ) V2FModel = 'ldm'

    V2FSigma(1:n) = GetReal( Material, 'V2-F Sigma', GotIt )
    IF ( .NOT. GotIt ) V2FSigma(1:n) = 1.0d0

    KECmu(1:n) = ListGetConstReal( Material, 'KE Cmu', GotIt )
    IF ( .NOT. GotIt ) KECmu(1:n) = 0.22_dp

    V2FC1(1:n) = GetReal( Material, 'V2-F C1', GotIt )
    IF ( .NOT. GotIt ) V2FC1(1:n) = 1.4_dp

    V2FC2(1:n) = GetReal( Material, 'V2-F C2', GotIt )
    IF ( .NOT. GotIt ) V2FC2(1:n) = 0.3_dp

    V2FCT(1:n) = GetReal( Material, 'V2-F CT', GotIt )
    IF ( .NOT. GotIt ) V2FCT(1:n) = 6.0_dp

    V2FCL(1:n) = GetReal( Material, 'V2-F CL', GotIt )
    IF ( .NOT. GotIt ) V2FCL(1:n) = 0.23_dp

    V2FCnu(1:n) = GetReal( Material, 'V2-F Cnu', GotIt )
    IF ( .NOT. GotIt ) V2FCnu(1:n) = 70.0d0

    Density(1:n)   = GetReal( Material,'Density' )
    Viscosity(1:n) = GetReal( Material,'Viscosity' )

    CALL GetScalarLocalSolution( LocalV2, 'V2' )
    CALL GetScalarLocalSolution( LocalF, 'F' )

    CALL GetScalarLocalSolution( U, 'Velocity 1' )
    CALL GetScalarLocalSolution( V, 'Velocity 2' )
    CALL GetScalarLocalSolution( W, 'Velocity 3' )

    CALL ElementKernel( MASS, STIFF, FORCE, LOAD, U, V, W, Element, n, nd+nb, ElementNodes )

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
!> Verbatim port of V2FSolverLegacy.F90's nested LocalMatrix. "Bubbles" and
!> "Material" come from the host (LocalMatrixScalar above), exactly as they
!> did from the legacy driver's own host scope.
!------------------------------------------------------------------------------
    SUBROUTINE ElementKernel( MASS,STIFF,FORCE, &
             LOAD,UX,UY,UZ,Element,n,nd,Nodes )
!------------------------------------------------------------------------------
      USE MaterialModels
      IMPLICIT NONE IMPLICIT_EXTERNAL

      REAL(KIND=dp), DIMENSION(:)   :: FORCE,UX,UY,UZ
      REAL(KIND=dp), DIMENSION(:,:) :: MASS,STIFF,LOAD

      INTEGER :: n,nd

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: ddBasisddx(nd,3,3)
      REAL(KIND=dp) :: Basis(nd)
      REAL(KIND=dp) :: dBasisdx(nd,3),detJ

      REAL(KIND=dp) :: Velo(3),dVelodx(3,3)

      REAL(KIND=dp) :: A(2,2),M(2,2)
      REAL(KIND=dp) :: LoadatIp(2),Rho,mu,Tmu,EffVisc
      INTEGER :: i,j,c,p,q,t,dim,N_Integ,NBasis

      REAL(KIND=dp) :: s,u,v,w, K,E,Eta,LV2,Lf,Strain(3,3),ProdTensor(3,3), &
                      Vorticity(3,3), nu, Cv, Prod,div

      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

      REAL(KIND=dp) :: C1,C2,CT,CL,Cnu,Cmu,TimeScale,LengthScale2,aparm(nd)
      REAL(KIND=dp) :: SecInv,X,Y,Z,Re_T
      REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

      REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

      LOGICAL :: stat,Convection
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      FORCE = 0.0_dp
      STIFF = 0.0_dp
      MASS  = 0.0_dp

      IF ( Bubbles ) THEN
         IntegStuff = GaussPoints( element, element % TYPE % GaussPoints2 )
      ELSE
         IntegStuff = GaussPoints( element )
      END IF

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n

      DO t=1,N_Integ
        u = U_Integ(t)
        v = V_Integ(t)
        w = W_Integ(t)
        stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                Basis,dBasisdx,Bubbles=Bubbles )

        s = detJ * S_Integ(t)
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          X = SUM( Nodes % x(1:n)*Basis(1:n) )
          Y = SUM( Nodes % y(1:n)*Basis(1:n) )
          Z = SUM( nodes % z(1:n)*Basis(1:n) )
          CALL CoordinateSystemInfo(Metric,SqrtMetric,Symb,dSymb,X,Y,Z)

          s = s * SqrtMetric
        END IF

        Velo = 0.0d0
        Velo(1) = SUM( UX(1:n)*Basis(1:n) )
        Velo(2) = SUM( UY(1:n)*Basis(1:n) )
        Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

        dVelodx = 0.0d0
        DO i=1,dim
          dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
          dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
          dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
        END DO

        Strain    = 0.5_dp * ( dVelodx + TRANSPOSE(dVelodx) )
        Vorticity = 0.5_dp * ( dVelodx - TRANSPOSE(dVelodx) )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           Secinv = 2*SUM(Strain * Strain)
        ELSE
           SecInv = SecondInvariant(Velo,dVelodx,Metric,Symb)/2
        END IF

        K   = MAX( SUM( LocalV2(1:n) * Basis(1:n) ), 1.0d-10 )
        E   = MAX( SUM( LocalF(1:n)  * Basis(1:n) ), 1.0d-10 )
        LF  = SUM( LocalF(1:n)  * Basis(1:n) )
        LV2 = SUM( LocalV2(1:n) * Basis(1:n) )

        Cnu = SUM( V2FCnu(1:n) * Basis(1:n) )
        Cmu = SUM( KECMu(1:n)  * Basis(1:n) )
        CT  = SUM( V2FCT(1:n)  * Basis(1:n) )
        CL  = SUM( V2FCL(1:n)  * Basis(1:n) )
        C1  = SUM( V2FC1(1:n)  * Basis(1:n) )
        C2  = SUM( V2FC2(1:n)  * Basis(1:n) )

        rho = SUM( Basis(1:n) * Density(1:n) )
        mu  = SUM( Basis(1:n) * Viscosity(1:n) )

        Re_T = K**2 / ((mu/Rho)*E)
        TimeScale = MAX( K/E, CT*SQRT((mu/rho)/E))
        Lengthscale2 = CL**2 * MAX(K**3 /E**2, Cnu**2*SQRT((mu/rho)**3/E))

        Tmu = rho*Cmu*LV2*TimeScale
        EffVisc = mu + Tmu / SUM(V2FSigma(1:n)*Basis(1:n))

        Prod = Tmu * SecInv / Rho

        DO p=1,nd
        DO q=1,nd
           M = 0.0d0
           A = 0.0d0

           M(1,1) = rho * Basis(q) * Basis(p)

           A(1,1) = A(1,1) + 6 * Rho / TimeScale  * Basis(q) * Basis(p)
           A(1,2) = A(1,2) - Rho * K * Basis(q) * Basis(p)

           A(2,1) = A(2,1) + (C1-6) / K / TimeScale * Basis(q) * Basis(p)
           A(2,2) = A(2,2) + Basis(q) * Basis(p)

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO i=1,dim
                A(1,1) = A(1,1) + EffVisc * dBasisdx(q,i) * dBasisdx(p,i)
                A(2,2) = A(2,2) + LengthScale2*dBasisdx(q,i) * dBasisdx(p,i)
              END DO
           ELSE
              DO i=1,dim
                DO j=1,dim
                   A(1,1) = A(1,1) + Metric(i,j) * EffVisc * &
                        dBasisdx(q,i) * dBasisdx(p,j)

                   A(2,2) = A(2,2) + Metric(i,j) * LengthScale2 * &
                        dBasisdx(q,i) * dBasisdx(p,j)
                END DO
              END DO
           END IF

           DO i=1,dim
             A(1,1) = A(1,1) + rho * Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO

           DO i=1,2
              DO j=1,2
                STIFF(2*(p-1)+i,2*(q-1)+j) = STIFF(2*(p-1)+i,2*(q-1)+j)+s*A(i,j)
                MASS(2*(p-1)+i,2*(q-1)+j)  = MASS(2*(p-1)+i,2*(q-1)+j) +s*M(i,j)
              END DO
           END DO
        END DO
        END DO

        LoadAtIP = 0.0_dp
        LoadAtIP(2) = ((C1-1)*2/3.0_dp+C2*Prod/E)/TimeScale

        DO p=1,nd
           FORCE(2*(p-1)+1) = FORCE(2*(p-1)+1) + s*LoadAtIp(1)*Basis(p)
           FORCE(2*(p-1)+2) = FORCE(2*(p-1)+2) + s*LoadAtIp(2)*Basis(p)
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE ElementKernel
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixScalar
!------------------------------------------------------------------------------

END MODULE V2FSolverLocalForms


!------------------------------------------------------------------------------
!> Vectorized/threaded V2-F driver. See the file header above for what this
!> does and does not support, and "Legacy Assembly" for the fallback.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE V2F_LDM( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE V2FSolverLocalForms
  USE V2FSolverFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC
  INTEGER :: i,n,nb,nd,Active,iter,NonlinearIter,nthr
  LOGICAL :: GotIt, InitHandles, GlobalBubbles, Stabilize, AxiSymmetric, BubblesDefault
  REAL(KIND=dp) :: Norm
  CHARACTER(*), PARAMETER :: Caller = 'V2F_LDM'
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToV2FSolverLegacy( 'V2F_LDM_Legacy', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN

  ! LocalMatrixVec now carries the plain "Axi Symmetric" (no swirl) case
  ! itself; genuine swirl ("Cylindric Symmetric") and general "Cylindric"
  ! still need LocalMatrixScalar's full metric/Christoffel treatment, so
  ! those still go through it, serially -- same branch HeatSolve.F90 makes,
  ! and the same treatment the other turbulence solvers in this directory
  ! now have. There is no "Compressibility Model" or alternate "V2-F Model"
  ! to route here -- Density is a plain material property and there is only
  ! the LDM formulation, in both the legacy LocalMatrix and here.
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian .AND. &
                   CurrentCoordinateSystem() /= AxisSymmetric )

  IF (.NOT. ALLOCATED(V2FHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(V2FHandles(nthr))
  END IF

  GlobalBubbles = Solver % GlobalBubbles
  Stabilize = GetStabilizeFlag( Solver % Values, GotIt )

  ! Only LocalMatrixScalar's legacy per-node branch uses this -- see the
  ! matching BubblesDefault resolution in V2FSolverLegacy.F90.
  BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  IF ( .NOT.GotIt ) BubblesDefault = .TRUE.

  ! K and Epsilon are interleaved, Dofs=2 -- exactly as V2FSolverLegacy.F90
  ! sizes its own bx/bxprev.
  IF ( TransientSimulation ) CALL DefaultBubbleHistoryUpdate( Dofs=2 )

  NonlinearIter = ListGetInteger( Solver % Values, 'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT.GotIt ) NonlinearIter = 1

  DO i=1,Model % NumberOfBCs
    BC => Model % BCs(i) % Values
    IF ( ListGetLogical( BC, 'Noslip wall BC', gotit ) ) THEN
      CALL ListAddConstReal( BC, 'F',  0.0_dp )
      CALL ListAddConstReal( BC, 'V2', 0.0_dp )
    END IF
  END DO

  DO iter=1,NonlinearIter
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,' ', Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)
    CALL Info(Caller,'V2-F iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)

    CALL DefaultInitialize()

    Active = GetNOFActive()

    IF ( AxiSymmetric ) THEN
      ! Serial only: LocalMatrixScalar relies on the classic argument-less
      ! Get*() "current element" accessors, exactly as
      ! V2FSolverLegacy.F90's own driver does.
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
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    ! V2 positive; F is left alone (legacy only clips V2, Values(1::2)).
    Solver % Variable % Values(1::2) = &
      MAX( Solver % Variable % Values(1::2), 1.0d-9 )

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE V2F_LDM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Element-type resolution, called before the mesh/basis functions are
!> finalized -- same role and timing as HeatSolver_Init0 and the other
!> turbulence solvers' own _Init0 in this directory.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE V2F_LDM_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE V2FSolverFront
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
END SUBROUTINE V2F_LDM_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: V2F_LDM. Under "Legacy Assembly",
!> delegates to V2F_LDM_Legacy_Init and returns -- the p-bubble setup below
!> is specific to this file's own implementation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE V2F_LDM_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE V2FSolverFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found, PBubble, LegacyBubbles
  CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToV2FSolverLegacy( 'V2F_LDM_Legacy_Init', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  SolverParams => GetSolverParams()

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
  ELSE IF ( TransientSimulation ) THEN
    LegacyBubbles = ListGetLogical( SolverParams, 'Bubbles', Found )
    IF ( .NOT. Found ) LegacyBubbles = .TRUE.

    IF ( LegacyBubbles ) THEN
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Min Iterations', 2)
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Max Iterations', 2)
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE V2F_LDM_Init
!------------------------------------------------------------------------------
