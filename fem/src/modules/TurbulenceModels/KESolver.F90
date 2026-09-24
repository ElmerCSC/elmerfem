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
! *  Original Date: 08 Jun 1997
! *
! ****************************************************************************/

!------------------------------------------------------------------------------
!> Solver for the k-epsilon turbulence model ("standard" and "v2-f"; "rng"
!> falls back to LocalMatrixScalar, see below). Vectorized/threaded
!> implementation following the same pattern as Komega.F90/SSTKomega.F90 and
!> Spalart-Allmaras.F90.
!>
!> Unlike Komega/SSTKomega, K and epsilon here are NOT matrix-diagonal only:
!> the legacy LocalMatrix linearizes the K equation's destruction term
!> implicitly in epsilon's own trial function -- A(1,2) = C0(1)*Basis(q)*
!> Basis(p), which becomes a genuine K-row/epsilon-column STIFF block
!> (A(2,1) is never written, so the coupling is one-directional). This still
!> needs no block-coupling machinery beyond what LinearForms already gives:
!> the cross block is built with one more LinearForms_UdotU call (weight =
!> rho, the same "C0(1) = Rho" legacy uses) and interleaved into
!> STIFF(1:2*ntot-1:2, 2:2*ntot:2) instead of a diagonal slot.
!>
!> Genuine swirl ("Cylindric Symmetric"), general "Cylindric", any non-default
!> "Compressibility Model", and "KE Model = RNG" (whose alpha root-solve isn't
!> worth vectorizing for a variant none of this codebase's tests exercise) go
!> through LocalMatrixScalar instead -- a scalar, per-Gauss-point fallback
!> carrying the same math as KESolverLegacy.F90's own LocalMatrix (already
!> interleaved directly via STIFF(2*(p-1)+i,2*(q-1)+j)), called serially --
!> same role as HeatSolve.F90's own AxiSymmetric branch. LocalMatrixScalar
!> also carries the legacy "Bubbles = True" per-node scheme. EpsilonWall
!> (the weak-form near-wall epsilon BC) and the "Wall Law" wall-function
!> treatment (via the core KEWall routine in Walls.F90) are boundary
!> treatments shared by both bulk paths, called once from the driver's own
!> boundary loop -- ported unchanged bar taking Solver/Model explicitly,
!> since they are module procedures here, not nested inside the driver.
!>
!> The original scalar-element solver lives on in KESolverLegacy.F90
!> (subroutine KESolverLegacy), reachable either directly by that name or
!> via "Legacy Assembly = Logical True" here (see KESolverFront below). With
!> axisymmetric, non-default compressibility, RNG and the per-node scheme
!> all covered above, that path is now only needed by a sif that must
!> reproduce the legacy solver's exact historical numbers.
!> \ingroup Solvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Whether this solver should run the original scalar-element assembly
!> (KESolverLegacy.F90) instead of this file's own implementation.
!------------------------------------------------------------------------------
MODULE KESolverFront
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
!> Call one of KESolverLegacy's entry points with this solver. The name is
!> resolved at run time, as the core resolves any solver, so this file and
!> KESolverLegacy.so stay independent of one another.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateToKESolverLegacy( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'KESolverLegacy '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'KESolver', &
        '"Legacy Assembly" was requested but "'//TRIM(Entry)//'" could not be found. '// &
        'Is KESolverLegacy.so installed beside this solver?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateToKESolverLegacy

END MODULE KESolverFront


!------------------------------------------------------------------------------
MODULE KESolverLocalForms

  USE DefUtils
  USE LinearForms

  IMPLICIT NONE IMPLICIT_EXTERNAL

  ! Per-element bubble history, used by Default1stOrderTime's Nb path
  ! (DefUtils.F90) -- lives on Solver % Variable's own BubbleValues/
  ! BubblePrevValues (Types.F90), not a separate type; see the matching
  ! bx/bxprev comment in KESolverLegacy.F90, which this mirrors.

  ! Per-thread ValueHandle_t storage for LocalMatrixVec's material lookups.
  ! NOT THREADPRIVATE -- see the matching comment on IncompressibleNS.F90's
  ! NSHandles_t for the Windows/GCC emutls hazard that rules that out.
  TYPE :: KEHandles_t
    TYPE(ValueHandle_t) :: Visc_h, Dens_h, SigmaK_h, SigmaE_h, Cmu_h, C1_h, C2_h, V2FCT_h
  END TYPE KEHandles_t
  TYPE(KEHandles_t), ALLOCATABLE, SAVE :: KEHandles(:)

CONTAINS

!------------------------------------------------------------------------------
!> Assemble and glue local matrix/RHS for one bulk element. Vectorized over
!> Gauss points, safe to call concurrently from multiple threads (each thread
!> passes its own InitHandles and only touches KEHandles(tid)).
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
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp), ALLOCATABLE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:)
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), TimeForce(:)
    REAL(KIND=dp), ALLOCATABLE :: MassK(:,:), StiffK(:,:), ForceK(:), &
        MassO(:,:), StiffO(:,:), ForceO(:), StiffKO(:,:)

    REAL(KIND=dp), POINTER :: RhoVec(:), MuVec(:), SigmaKVec(:), SigmaOVec(:), &
        CmuVec(:), C1Vec(:), C2Vec(:), V2FCTVec(:)

    REAL(KIND=dp), ALLOCATABLE :: VeloNodal(:,:), KNodal(:), ONodal(:), &
        V2Nodal(:), PressureNodal(:)
    REAL(KIND=dp), ALLOCATABLE :: VeloVec(:,:), dVelodxVec(:,:,:), KVec(:), OVec(:), &
        V2Vec(:), PressureVec(:), StrainVec(:,:,:), SecInvVec(:), &
        TmuVec(:), TimeScaleVec(:), Effmu1Vec(:), Effmu2Vec(:), SoundSpeedSqVec(:), &
        MachSqVec(:), ProdKVec(:), ProdEVec(:), ReactO(:), &
        LoadK(:), LoadO(:), StreamVec(:,:), TauVec(:), TmpVec(:), RadiusVec(:)

    CHARACTER(LEN=MAX_NAME_LEN) :: KEModelStr
    LOGICAL :: UseV2F
    REAL(KIND=dp) :: hK, mK, VNorm, SpecificHeatRatio, ReferencePressure, V2FCp
    INTEGER :: i,j,k,p,ngp,dim,allocstat,tid,ntot
    LOGICAL :: Stat, Found, IsAxiSymmetric
!------------------------------------------------------------------------------
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Visc_h => KEHandles(tid) % Visc_h, Dens_h => KEHandles(tid) % Dens_h, &
               SigmaK_h => KEHandles(tid) % SigmaK_h, SigmaE_h => KEHandles(tid) % SigmaE_h, &
               Cmu_h => KEHandles(tid) % Cmu_h, C1_h => KEHandles(tid) % C1_h, &
               C2_h => KEHandles(tid) % C2_h, V2FCT_h => KEHandles(tid) % V2FCT_h )

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity' )
      CALL ListInitElementKeyword( Dens_h,'Material','Density' )
      CALL ListInitElementKeyword( SigmaK_h,'Material','KE SigmaK' )
      CALL ListInitElementKeyword( SigmaE_h,'Material','KE SigmaE' )
      CALL ListInitElementKeyword( Cmu_h,'Material','KE Cmu' )
      CALL ListInitElementKeyword( C1_h,'Material','KE C1' )
      CALL ListInitElementKeyword( C2_h,'Material','KE C2' )
      CALL ListInitElementKeyword( V2FCT_h,'Material','V2-F CT' )
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ntot = nd + nb

    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ALLOCATE( BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        MassK(ntot,ntot), StiffK(ntot,ntot), ForceK(ntot), &
        MassO(ntot,ntot), StiffO(ntot,ntot), ForceO(ntot), StiffKO(ntot,ntot), &
        MASS(2*ntot,2*ntot), STIFF(2*ntot,2*ntot), FORCE(2*ntot), TimeForce(2*ntot), &
        VeloNodal(3,n), KNodal(n), ONodal(n), V2Nodal(n), PressureNodal(n), &
        VeloVec(ngp,3), dVelodxVec(ngp,3,3), KVec(ngp), OVec(ngp), V2Vec(ngp), &
        PressureVec(ngp), StrainVec(ngp,3,3), SecInvVec(ngp), &
        TmuVec(ngp), TimeScaleVec(ngp), Effmu1Vec(ngp), Effmu2Vec(ngp), SoundSpeedSqVec(ngp), &
        MachSqVec(ngp), ProdKVec(ngp), ProdEVec(ngp), ReactO(ngp), &
        LoadK(ngp), LoadO(ngp), StreamVec(ngp,ntot), TauVec(ngp), TmpVec(ngp), &
        RadiusVec(ngp), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('KESolver','Local storage allocation failed')

    CALL GetElementNodesVec( Nodes, UElement=Element )

    MassK = 0._dp; StiffK = 0._dp; ForceK = 0._dp
    MassO = 0._dp; StiffO = 0._dp; ForceO = 0._dp; StiffKO = 0._dp

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
    CALL GetScalarLocalSolution( PressureNodal, 'Pressure', UElement=Element )

    RhoVec => ListGetElementRealVec( Dens_h, ngp, BasisVec, Element, Found )
    MuVec  => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element, Found )

    SigmaKVec => ListGetElementRealVec( SigmaK_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) SigmaKVec = 1.0_dp

    ! KE Model resolves the defaults for SigmaE/Cmu/C1/C2 below and picks the
    ! K/epsilon closure (standard vs. v2-f). "rng" (and anything unrecognized)
    ! is refused by the driver before this is ever called -- see
    ! UseScalarFallback in KESolver.
    Material => GetMaterial( Element )
    CALL GetStringThreadSafe( Material, 'KE Model', KEModelStr, Found )
    IF( .NOT. Found ) KEModelStr = 'standard'
    UseV2F = ( KEModelStr == 'v2-f' )

    SigmaOVec => ListGetElementRealVec( SigmaE_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) THEN
      IF( UseV2F ) THEN
        SigmaOVec = 1.3_dp
      ELSE
        SigmaOVec = 1.3_dp
      END IF
    END IF

    CmuVec => ListGetElementRealVec( Cmu_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) THEN
      IF( UseV2F ) THEN
        CmuVec = 0.22_dp
      ELSE
        CmuVec = 0.09_dp
      END IF
    END IF

    C1Vec => ListGetElementRealVec( C1_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) THEN
      IF( UseV2F ) THEN
        C1Vec = 1.4_dp
      ELSE
        C1Vec = 1.44_dp
      END IF
    END IF

    C2Vec => ListGetElementRealVec( C2_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) THEN
      IF( UseV2F ) THEN
        C2Vec = 1.9_dp
      ELSE
        C2Vec = 1.92_dp
      END IF
    END IF

    ! Per-material constants (not spatially varying) -- same GetCReal calls
    ! and same (Elmer default) 0.0 fallback as the legacy LocalMatrix, so
    ! Mach_number_sq below comes out zero exactly as it does there whenever
    ! "Specific Heat Ratio" is unset. Buoyancy (rho_g) is exactly zero here
    ! -- the vectorized path only supports a spatially uniform Density, see
    ! the file header -- so "Turbulent Prandtl Number"/"Dissipation buoyancy
    ! coefficient" (which only ever multiply rho_g) are not needed here.
    SpecificHeatRatio = GetCReal( Material, 'Specific Heat Ratio', Found )
    ReferencePressure = GetCReal( Material, 'Reference Pressure', Found )

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
    PressureVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), PressureNodal(1:n) )

    ! Strain tensor and SecInv = 2*Strain:Strain, Cartesian only. Not floored
    ! -- legacy's own LocalMatrix does not floor SecInv either (unlike
    ! Komega/SSTKomega's StrainMeasure).
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

    IF( UseV2F ) THEN
      CALL GetScalarLocalSolution( V2Nodal, 'V2', UElement=Element )
      V2Vec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), V2Nodal(1:n) )

      V2FCTVec => ListGetElementRealVec( V2FCT_h, ngp, BasisVec, Element, Found )
      IF( .NOT. Found ) V2FCTVec = 6.0_dp

      V2FCp = GetCReal( Material, 'V2-F Cp', Found )
      IF( .NOT. Found ) V2FCp = 0.05_dp

      ! C1 is overridden for v2-f (legacy: LC1 = 1.4*(1+Cp*SQRT(K/V2))),
      ! regardless of whatever "KE C1" resolved to above.
      C1Vec(1:ngp) = 1.4_dp * ( 1._dp + V2FCp*SQRT(KVec(1:ngp)/V2Vec(1:ngp)) )

      TimeScaleVec(1:ngp) = MAX( KVec(1:ngp)/OVec(1:ngp), &
          V2FCTVec(1:ngp)*SQRT(MuVec(1:ngp)/RhoVec(1:ngp)/OVec(1:ngp)) )
      TmuVec(1:ngp) = RhoVec(1:ngp) * CmuVec(1:ngp) * V2Vec(1:ngp) * TimeScaleVec(1:ngp)
    ELSE
      TmuVec(1:ngp) = RhoVec(1:ngp) * CmuVec(1:ngp) * KVec(1:ngp)**2 / OVec(1:ngp)
    END IF

    ! rho_g omitted (zero, see above): ProdK/ProdE reduce to Tmu*SecInv/rho.
    ProdKVec(1:ngp) = TmuVec(1:ngp) * SecInvVec(1:ngp) / RhoVec(1:ngp)
    ProdEVec(1:ngp) = ProdKVec(1:ngp)

    Effmu1Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaKVec(1:ngp)
    Effmu2Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaOVec(1:ngp)

    SoundSpeedSqVec(1:ngp) = (PressureVec(1:ngp)+ReferencePressure) * SpecificHeatRatio / RhoVec(1:ngp)
    MachSqVec = 0._dp
    WHERE( SoundSpeedSqVec(1:ngp) > 0._dp ) MachSqVec(1:ngp) = KVec(1:ngp) / SoundSpeedSqVec(1:ngp)

    ! K-equation destruction (rho*epsilon) is linearized IMPLICITLY in
    ! epsilon's own trial function -- StiffKO below, not a reaction
    ! coefficient here -- so ReactK does not exist; only epsilon's own
    ! destruction (C0(2)) is a plain diagonal reaction.
    IF( UseV2F ) THEN
      ReactO(1:ngp) = RhoVec(1:ngp) * C2Vec(1:ngp) / TimeScaleVec(1:ngp)
      LoadO(1:ngp)  = RhoVec(1:ngp) * C1Vec(1:ngp) * ProdEVec(1:ngp) / TimeScaleVec(1:ngp)
    ELSE
      ReactO(1:ngp) = RhoVec(1:ngp) * C2Vec(1:ngp) * OVec(1:ngp) / KVec(1:ngp)
      LoadO(1:ngp)  = RhoVec(1:ngp) * C1Vec(1:ngp) * ProdEVec(1:ngp) * OVec(1:ngp) / KVec(1:ngp)
    END IF

    LoadK(1:ngp) = RhoVec(1:ngp) * ( ProdKVec(1:ngp) - 2._dp*OVec(1:ngp)*MachSqVec(1:ngp) )

    ! MASS/STIFF diagonal blocks: mass, diffusion, convection -- same plain
    ! Velo for both equations (no extra cross-diffusion correction here,
    ! unlike SSTKomega's omega equation).
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, MassK, RhoVec )
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

    ! K-row/epsilon-column cross block: C0(1) = Rho, exactly legacy's A(1,2).
    CALL LinearForms_UdotU( ngp, ntot, dim, BasisVec, DetJVec, StiffKO, RhoVec )

    !------------------------------------------------------------------------
    ! SUPG (equal-order) stabilization, opt-in via "Stabilize"/"Stabilization
    ! Method" -- same Franca et al. tau/streamline construction as
    ! Komega.F90's own SUPG block (K and epsilon share the same plain Velo,
    ! so one StreamVec/Tau pair per equation's own Effmu, as there). Only
    ! the diagonal blocks are stabilized -- SUPG targets the convection-
    ! diffusion operator, not the implicit reaction cross term above.
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
          TauVec(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(Effmu1Vec(j))) )
          TauVec(j) = hK * TauVec(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauVec(j) = 0._dp
        END IF
      END DO
      TmpVec(1:ngp) = TauVec(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, StiffK )
      CALL LinearForms_UdotF( ngp, ntot, StreamVec, TmpVec, LoadK, ForceK )
      IF( Transient ) THEN
        TmpVec(1:ngp) = TauVec(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MassK )
      END IF

      DO j=1,ngp
        VNorm = SQRT( SUM( VeloVec(j,1:dim)**2 ) )
        IF( VNorm > 0._dp .AND. Effmu2Vec(j) /= 0._dp ) THEN
          TauVec(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(Effmu2Vec(j))) )
          TauVec(j) = hK * TauVec(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauVec(j) = 0._dp
        END IF
      END DO
      TmpVec(1:ngp) = TauVec(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, StreamVec, TmpVec, StiffO )
      CALL LinearForms_UdotF( ngp, ntot, StreamVec, TmpVec, LoadO, ForceO )
      IF( Transient ) THEN
        TmpVec(1:ngp) = TauVec(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVec, BasisVec, TmpVec, MassO )
      END IF
    END IF

    ! Interleave: odd rows/columns for K, even for epsilon. The (K,epsilon)
    ! cross block goes into STIFF(1:2*ntot-1:2, 2:2*ntot:2); (epsilon,K)
    ! stays zero, matching legacy's own A(2,1) (never written).
    MASS = 0._dp; STIFF = 0._dp; FORCE = 0._dp
    MASS(1:2*ntot-1:2,1:2*ntot-1:2) = MassK
    MASS(2:2*ntot:2,  2:2*ntot:2)   = MassO
    STIFF(1:2*ntot-1:2,1:2*ntot-1:2) = StiffK
    STIFF(2:2*ntot:2,  2:2*ntot:2)   = StiffO
    STIFF(1:2*ntot-1:2,2:2*ntot:2)   = StiffKO
    FORCE(1:2*ntot-1:2) = ForceK
    FORCE(2:2*ntot:2)   = ForceO

    !------------------------------------------------------------------------
    ! Time discretization and p-bubble condensation -- mirrors the nb>0
    ! branch of KESolverLegacy.F90's driver exactly (DOFs=2, K/epsilon
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
!> do: axisymmetric/cylindrical coordinates, any non-default "Compressibility
!> Model", and "KE Model = RNG" (the ElementKernel below carries the same
!> metric-tensor, ElementDensity and RNG alpha-root-solve as
!> KESolverLegacy.F90's own LocalMatrix -- verbatim, no new math, and already
!> interleaves K/epsilon directly via STIFF(2*(p-1)+i,2*(q-1)+j)), and the
!> legacy "Bubbles = True" per-node scheme for a plain nodal "Element" set
!> explicitly. Always called serially (see the UseScalarFallback branch in
!> KESolver below), so it uses the classic GetMaterial()/GetReal()/
!> argument-less GetElementNOF*() accessors exactly as the legacy driver
!> does -- not safe to call from inside an OMP parallel region, unlike
!> LocalMatrixVec.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixScalar( Element, dt, Transient, GlobalBubbles, BubblesDefault )
!------------------------------------------------------------------------------
    IMPLICIT NONE IMPLICIT_EXTERNAL
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: Transient, GlobalBubbles, BubblesDefault
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    TYPE(ValueList_t), POINTER :: Material, Equation
    CHARACTER(LEN=MAX_NAME_LEN) :: KEModel, V2FModel
    REAL(KIND=dp) :: Clip
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:), &
        TimeForce(:), LocalKinEnergy(:), LocalDissipation(:), &
        U(:), V(:), W(:), Density(:), Viscosity(:), &
        KESigmaK(:), KESigmaE(:), KECmu(:), KEC1(:), KEC2(:), LocalV2(:), V2FCT(:)
    REAL(KIND=dp) :: V2FCp
    LOGICAL :: Bubbles, GotIt
    INTEGER :: n, nd, nb, allocstat
!------------------------------------------------------------------------------
    Bubbles = BubblesDefault .AND. .NOT. ASSOCIATED( Element % PDefs )
    Material => GetMaterial()
    Equation => GetEquation()

    Clip = GetConstReal( Material, 'KE Clip', GotIt )
    IF ( .NOT.GotIt ) Clip = 1.0d-6

    KEModel = GetString( Material, 'KE Model', GotIt )
    IF ( .NOT. GotIt ) KEModel = 'standard'

    V2FModel = GetString( Material, 'V2-F Model', GotIt )

    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
    IF ( Bubbles ) nd = 2*n
    nb = GetElementNOFBDOFs()

    CALL GetElementNodes( ElementNodes )

    ALLOCATE( MASS(2*(nd+nb),2*(nd+nb)), STIFF(2*(nd+nb),2*(nd+nb)), &
        FORCE(2*(nd+nb)), LOAD(2,n), TimeForce(2*(nd+nb)), &
        LocalKinEnergy(nd+nb), LocalDissipation(nd+nb), &
        U(n), V(n), W(n), Density(n), Viscosity(n), &
        KESigmaK(n), KESigmaE(n), KECmu(n), KEC1(n), KEC2(n), LocalV2(n), V2FCT(n), &
        STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('KESolver','Local storage allocation failed')

    CALL GetScalarLocalSolution( LocalV2, 'V2' )
    CALL GetScalarLocalSolution( LocalKinEnergy, 'Kinetic Energy' )
    CALL GetScalarLocalSolution( LocalDissipation, 'Kinetic Dissipation' )

    CALL GetScalarLocalSolution( U, 'Velocity 1' )
    CALL GetScalarLocalSolution( V, 'Velocity 2' )
    CALL GetScalarLocalSolution( W, 'Velocity 3' )

    KESigmaK(1:n) = GetReal( Material, 'KE SigmaK', GotIt )
    IF ( .NOT. GotIt ) KESigmaK = 1.0d0

    KESigmaE(1:n) = GetReal( Material, 'KE SigmaE', GotIt )
    IF ( .NOT. GotIt ) THEN
      SELECT CASE( KEModel )
      CASE( 'standard','v2-f' )
        KESigmaE = 1.3_dp
      CASE( 'rng' )
        KESigmaE = 1.0_dp
      CASE DEFAULT
        CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
      END SELECT
    END IF

    KECmu(1:n) = ListGetConstReal( Material, 'KE Cmu', GotIt )
    IF ( .NOT. GotIt ) THEN
      SELECT CASE( KEModel )
      CASE( 'standard' )
        KECmu = 0.09_dp
      CASE( 'v2-f')
        KECmu = 0.22_dp
      CASE( 'rng' )
        KECmu = 0.0845_dp
      CASE DEFAULT
        CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
      END SELECT
    END IF

    KEC1(1:n) = GetReal( Material, 'KE C1', GotIt )
    IF ( .NOT. GotIt ) THEN
      SELECT CASE( KEModel )
      CASE( 'standard' )
        KEC1 = 1.44_dp
      CASE( 'v2-f' )
        KEC1(1:n) = 1.4_dp
      CASE( 'rng' )
        KEC1 = 1.42_dp
      CASE DEFAULT
        CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
      END SELECT
    END IF

    KEC2(1:n) = GetReal( Material, 'KE C2', GotIt )
    IF ( .NOT. GotIt ) THEN
      SELECT CASE( KEModel )
      CASE( 'standard' )
        KEC2 = 1.92_dp
      CASE( 'v2-f' )
        KEC2 = 1.9_dp
      CASE( 'rng' )
        KEC2 = 1.68_dp
      CASE DEFAULT
        CALL Fatal( 'KESolver', 'Unknown K-Epsilon model' )
      END SELECT
    END IF

    IF ( KEModel == 'v2-f' ) THEN
      V2FCT(1:n) = GetReal( Material, 'V2-F CT', GotIt )
      IF ( .NOT. GotIt ) V2FCT(1:n) = 6.0_dp

      V2FCp = GetCReal( Material, 'V2-F Cp', GotIt )
      IF ( .NOT. GotIt ) V2FCp = 0.05_dp
    END IF

    Density(1:n)   = GetReal( Material, 'Density' )
    Viscosity(1:n) = GetReal( Material, 'Viscosity' )

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
!> Verbatim port of KESolverLegacy.F90's nested LocalMatrix. "Bubbles",
!> "Material", "KEModel", the "KE*"/V2F nodal arrays etc. come from the host
!> (LocalMatrixScalar above), exactly as they did from the legacy driver's
!> own host scope.
!------------------------------------------------------------------------------
    SUBROUTINE ElementKernel( MASS,STIFF,FORCE, &
          LOAD, UX,UY,UZ,Element,n,nd,Nodes )
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

      REAL(KIND=dp) :: A(2,2),M(2,2),ProdK,ProdE,div,ProdTensor(3,3)
      INTEGER :: i,j,c,p,q,t,dim,NBasis
      REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2),TimeScale,Re_T,LC1,LC2,LC3

      REAL(KIND=dp) :: s,u,v,w, K,E,Eta,Strain(3,3), GradAsymm(3,3), nu, Cv,&
              alpha,oldalpha,dalpha,err,ww,olderr,derr

      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

      REAL(KIND=dp) :: SpecificHeatRatio, Pressure(n), ReferencePressure, &
           Sound_speed_sq, Mach_number_sq
      REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
               SqrtMetric,Gravity(3),Pr_rho(n),Pr,rho_g, c3(n)
      REAL(KIND=dp) :: C0(2),C1,CT,C2(2),dC2dx(3),SecInv,X,Y,Z,LV2,LCT,SigmaK,SigmaE

      REAL(KIND=dp), POINTER :: gWork(:,:)

      LOGICAL :: stat,Convection, UseRNGModel, GotIt
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      FORCE = 0.0D0
      STIFF = 0.0D0
      MASS  = 0.0D0

      NBasis = nd
      IF ( Bubbles ) NBasis = 2*n

      UseRNGModel = KEModel == 'rng'

      ! CurrentModel in place of Model: ElementKernel is nested inside
      ! LocalMatrixScalar, a module procedure, not inside the driver the way
      ! the legacy nested LocalMatrix was.
      gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',GotIt)
      IF ( GotIt ) THEN
        Gravity = gWork(1:3,1)*gWork(4,1)
      ELSE
        Gravity    =  0.00_dp
        Gravity(2) = -9.81_dp
      END IF

      Pr_rho(1:n) = GetReal( Material, 'Turbulent Prandtl Number', stat )
      IF ( .NOT. stat ) Pr_rho(1:n) = 0.85_dp

      c3(1:n) = GetReal( Material, 'Dissipation buoyancy coefficient', stat )
      IF ( .NOT. stat ) c3(1:n) = 0.0_dp

      CALL ElementDensity( Density, n )
      SpecificheatRatio = GetCReal( Material, 'Specific Heat Ratio', stat )
      CALL getScalarLocalSolution( Pressure, 'Pressure' )
      ReferencePressure = GetCReal( Material, 'Reference Pressure', stat )

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

        Velo = 0.0D0
        Velo(1) = SUM( UX(1:n)*Basis(1:n) )
        Velo(2) = SUM( UY(1:n)*Basis(1:n) )
        Velo(3) = SUM( UZ(1:n)*Basis(1:n) )

        dVelodx = 0.0d0
        DO i=1,dim
          dVelodx(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
          dVelodx(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
          dVelodx(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
        END DO

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
          Strain = 0.5d0 * ( dVelodx + TRANSPOSE(dVelodx) )
          Secinv = 2 * SUM( Strain * Strain )
        ELSE
          SecInv = SecondInvariant( Velo,dVelodx,Metric,Symb ) / 2
        END IF

        K = SUM( LocalKinEnergy(1:n) * Basis(1:n) )
        E = SUM( LocalDissipation(1:n) * Basis(1:n) )
        Eta =  SQRT(SecInv) * K / E

        Pr = SUM(Pr_rho(1:n) * Basis(1:n) )
        rho_g = 0._dp
        DO i=1,dim
          rho_g = rho_g + SUM(Density(1:n) * dBasisdx(1:n,i)) * Gravity(i)
        END DO

        mu  = SUM( Viscosity(1:n) * Basis(1:n) )
        rho = SUM( Density(1:n) * Basis(1:n) )

        SigmaK = SUM( KESigmaK(1:n) * Basis(1:n) )
        SigmaE = SUM( KESigmaE(1:n) * Basis(1:n) )

        Cmu = SUM( KECMu(1:n) * Basis(1:n) )

        LC1 = SUM( KEC1(1:n) * Basis(1:n) )
        LC2 = SUM( KEC2(1:n) * Basis(1:n) )
        LC3 = SUM( c3(1:n) * Basis(1:n) )

        Sound_speed_sq = (SUM(Basis(1:n)*Pressure(1:n))+ReferencePressure) * &
                      SpecificHeatRatio / rho
        Mach_number_sq = 0._dp
        IF ( Sound_speed_sq > 0._dp ) Mach_number_sq = K/Sound_speed_sq

        IF ( KEModel=='v2-f' ) THEN
          LV2 = SUM( LocalV2(1:n) * Basis(1:n) )
          LCT = SUM( V2FCT(1:n) * Basis(1:n) )
          LC1 = 1.4_dp * (1+V2FCp*SQRT(K/LV2))

          Timescale = MAX(K/E,LCT*SQRT(mu/rho/E))
          Tmu = Rho * Cmu * LV2 * TimeScale
        ELSE
          Tmu  = Rho * Cmu*K**2  / E
        END IF
        ProdK = Tmu * (SecInv-rho_g/(rho*Pr)) / rho
        ProdE = Tmu * (SecInv-LC3*rho_g/(rho*Pr)) / rho

        Effmu(1) = mu + Tmu / SigmaK
        Effmu(2) = mu + Tmu / SigmaE

        C0(1) = Rho
        IF ( KEModel == 'v2-f' ) THEN
          C0(2) = Rho * LC2 / TimeScale
        ELSE
          C0(2) = Rho * LC2 * E / K
        END IF

        C1 = Rho
        CT = Rho

        Alpha = 1.0d0

        IF ( UseRNGModel ) THEN
           ww = mu / Effmu(1)
           alpha = 1.3929d0
           oldalpha = 1

           olderr = ABS((oldalpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
           olderr = olderr * ABS((oldalpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
           olderr = olderr - ww

           DO i=1,100
              err = ABS((alpha-1.3929d0)/(1.0d0-1.3929d0))**0.6321d0
              err = err * ABS((alpha+2.3929d0)/(1.0d0+2.3929d0))**0.3679d0
              err = err - ww
              derr = olderr - err
              olderr = err
              dalpha = oldalpha - alpha
              oldalpha = alpha
              alpha = alpha - 0.5 * err * dalpha / derr
              IF ( ABS(err) < 1.0d-8 ) EXIT
           END DO

           IF ( ABS(err) > 1.0d-8 ) THEN
              alpha = 1.3929d0
           END IF
        END IF

        C2(1) = Alpha * Effmu(1)
        C2(2) = Alpha * Effmu(2)

        DO p=1,NBasis
        DO q=1,NBasis
           M = 0.0d0
           A = 0.0d0

           M(1,1) = CT * Basis(q) * Basis(p)
           M(2,2) = CT * Basis(q) * Basis(p)

           A(1,2) = C0(1) * Basis(q) * Basis(p)
           A(2,2) = C0(2) * Basis(q) * Basis(p)

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO i=1,dim
                A(1,1) = A(1,1) + C2(1) * dBasisdx(q,i) * dBasisdx(p,i)
                A(2,2) = A(2,2) + C2(2) * dBasisdx(q,i) * dBasisdx(p,i)
              END DO
           ELSE
              DO i=1,dim
                DO j=1,dim
                   A(1,1) = A(1,1) + Metric(i,j) * C2(1) * &
                        dBasisdx(q,i) * dBasisdx(p,j)

                   A(2,2) = A(2,2) + Metric(i,j) * C2(2) * &
                        dBasisdx(q,i) * dBasisdx(p,j)
                END DO
              END DO
           END IF

           DO i=1,dim
             A(1,1) = A(1,1) + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
             A(2,2) = A(2,2) + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
           END DO

           DO i=1,2
              DO j=1,2
                STIFF(2*(p-1)+i,2*(q-1)+j) = STIFF(2*(p-1)+i,2*(q-1)+j)+s*A(i,j)
                MASS(2*(p-1)+i,2*(q-1)+j)  = MASS(2*(p-1)+i,2*(q-1)+j) +s*M(i,j)
              END DO
           END DO
         END DO
         END DO

         LoadAtIP(1) = rho*(ProdK-2*E*Mach_number_sq)
         IF ( KEModel=='v2-f' ) THEN
           LoadAtIP(2) = Rho*LC1*ProdE/TimeScale
         ELSE
           LoadAtIP(2) = Rho*LC1*ProdE*E/K
         END IF

         IF ( UseRNGModel ) &
             LoadatIP(2) = LoadatIP(2) - Cmu*Rho*Eta**3*(1-Eta/4.38d0) / &
                        (1.0d0 + 0.012d0*Eta**3) * E**2 / K

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
!> Weak-form near-wall epsilon BC, ported unchanged from KESolverLegacy.F90's
!> own nested EpsilonWall bar taking Solver/Transient explicitly (module
!> procedure here, not nested in the driver) and self-contained local
!> Density/Viscosity (the legacy version shared the driver's own N-sized
!> SAVE arrays; np is bounded by the same 32-node assumption BasisB etc.
!> already make).
!------------------------------------------------------------------------------
  SUBROUTINE EpsilonWall( Element, n, Solver, Transient )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
    TYPE(Solver_t) :: Solver
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: BasisB(32), Basis(32), BasisK(32), dBasisdx(32,3), &
                EVals(32), KVals(32), Density(32), Viscosity(32), detJ
    REAL(KIND=dp) :: MASS(2*32,2*32), STIFF(2*32,2*32), FORCE(2*32)

    INTEGER :: i,j,c,p,q,t,np, dim,N_Integ,NBasis

    REAL(KIND=dp) :: s,u,v,w,E,K,Kder,X,Y,Z,Normal(3),mu,rho,Relax

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    LOGICAL :: stat
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes, ParentNodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    FORCE = 0.0_dp
    STIFF = 0.0_dp

    BC => GetBC( Element )

    IntegStuff = GaussPoints( Element )

    Relax = GetCReal( BC, 'Epsilon Relax', stat )
    IF (.NOT. stat ) Relax = 1

    Parent => Element % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) &
      Parent => Element % BoundaryInfo % Right
    IF(.NOT.ASSOCIATED(Parent))RETURN

    np = GetElementNOFDOFs(Parent)
    CALL GetElementNodes( Nodes, Element )
    CALL GetElementNodes( ParentNodes, Parent )

    Density(1:np)   = GetReal( GetMaterial(Parent), 'Density',UElement=Parent )
    Viscosity(1:np) = GetReal( GetMaterial(Parent), 'Viscosity',UElement=Parent )

    CALL GetScalarLocalSolution( KVals, 'Kinetic Energy', Parent )
    CALL GetScalarLocalSolution( EVals, 'Kinetic dissipation', Parent )

    DO t=1,IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)

      stat = ElementInfo( Element,Nodes,u,v,w,detJ, BasisB )

      s = detJ * IntegStuff % s(t)
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        x = SUM( Nodes % x(1:n)*BasisB(1:n) )
        y = SUM( Nodes % y(1:n)*BasisB(1:n) )
        z = SUM( Nodes % z(1:n)*BasisB(1:n) )
        s = s *  CoordinateSqrtMetric(x,y,z)
      END IF

      Normal = NormalVector( Element, Nodes, U, V, .TRUE. )

      CALL GetParentUVW( Element,n,Parent,np,U,V,W,BasisB )
      stat = ElementInfo( Parent,ParentNodes,U,V,W,detJ, &
            Basis,dBasisdx )

      IF ( ABS(u)>0.999_dp ) THEN
        u = -u
      ELSE IF ( ABS(v)>0.999_dp ) THEN
        v = -v
      ELSE IF ( ABS(w)>0.999_dp ) THEN
        w = -w
      END IF
      stat = ElementInfo( Parent,ParentNodes,U,V,W,detJ,BasisK )

      rho = SUM( Basis(1:np) * Density(1:np) )
      mu  = SUM( Basis(1:np) * Viscosity(1:np) )

      E = SUM( Basis(1:np) *Evals(1:np) )
      K = SUM( BasisK(1:np)*Kvals(1:np) )
      KVals(1:np) = SQRT(KVals(1:np))
      Kder = 0.0_dp
      DO i=1,3
        Kder = Kder+SUM(dBasisdx(1:np,i)*Kvals(1:np))*Normal(i)
      END DO

      DO p=1,np
        DO q=1,np
          STIFF(2*p,2*q)   = STIFF(2*p,2*q) + s*Basis(q)*Basis(p)
          STIFF(2*p,2*q-1) = STIFF(2*p,2*q-1) - s*Relax*2*mu/rho*Kder**2*BasisK(q)/MAX(K,AEPS)*Basis(p)
        END DO
        FORCE(2*p) = FORCE(2*p) + s*(1-Relax)*E*Basis(p)
      END DO
    END DO

    IF ( Transient ) THEN
      MASS = 0.0_dp
      CALL Default1stOrderTime(MASS, STIFF, FORCE, UElement=Parent)
    END IF
    CALL DefaultUpdateEquations(STIFF, FORCE, UElement=Parent)
!------------------------------------------------------------------------------
  END SUBROUTINE EpsilonWall
!------------------------------------------------------------------------------

END MODULE KESolverLocalForms


!------------------------------------------------------------------------------
!> Vectorized/threaded k-epsilon driver. See the file header above for what
!> this does and does not support, and "Legacy Assembly" for the fallback.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KESolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KESolverLocalForms
  USE KESolverFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
  EXTERNAL :: KEWALL
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: KE, FlowSol
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: BC, Material
  INTEGER, POINTER :: FlowPerm(:), KinPerm(:)
  REAL(KIND=dp), POINTER :: FlowSolution(:)
  INTEGER :: i,j,k,n,nb,nd,Active,iter,NonlinearIter,nthr,NSDOFs
  LOGICAL :: GotIt, InitHandles, GlobalBubbles, Stabilize, UseScalarFallback, BubblesDefault
  REAL(KIND=dp) :: Norm, KVal, EVal, KMax, EMax, Clip
  REAL(KIND=dp) :: U(32), V(32), W(32), Density(32), Viscosity(32), &
      SurfaceRoughness(32), LayerThickness(32), Work(3)
  REAL(KIND=dp), ALLOCATABLE, SAVE :: NodalDensity(:), NodalViscosity(:), NodalCmu(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: KEModelStr
  CHARACTER(*), PARAMETER :: Caller = 'KESolver'
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToKESolverLegacy( 'KESolverLegacy', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  KE => Solver % Variable
  IF ( .NOT. ASSOCIATED(KE) ) RETURN
  IF ( COUNT( KE % Perm > 0 ) <= 0 ) RETURN
  KinPerm => KE % Perm

  FlowSol => VariableGet( Model % Variables, 'Flow Solution' )
  IF ( ASSOCIATED( FlowSol ) ) THEN
    FlowPerm     => FlowSol % Perm
    NSDOFs       =  FlowSol % DOFs
    FlowSolution => FlowSol % Values
  END IF

  ! LocalMatrixVec has no metric tensor, only supports a spatially uniform
  ! Density (its buoyancy term is exactly zero, see the file header), and
  ! its closure coefficients only cover "standard" and "v2-f" -- "rng" needs
  ! LocalMatrixScalar's own root-solve for alpha. Resolve "KE Model" the
  ! same way LocalMatrixScalar's Init-time check would (ListCheckPresentAny
  ! Material can't distinguish which model string a material has, only
  ! whether the keyword is present at all, so check every material's own
  ! value here).
  UseScalarFallback = ( CurrentCoordinateSystem() /= Cartesian .AND. &
      CurrentCoordinateSystem() /= AxisSymmetric ) .OR. &
      ListCheckPresentAnyMaterial( Model, 'Compressibility Model' )
  IF( .NOT. UseScalarFallback ) THEN
    DO i=1,Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      KEModelStr = ListGetString( Material, 'KE Model', GotIt )
      IF( GotIt .AND. KEModelStr /= 'standard' .AND. KEModelStr /= 'v2-f' ) THEN
        UseScalarFallback = .TRUE.
        EXIT
      END IF
    END DO
  END IF

  IF (.NOT. ALLOCATED(KEHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(KEHandles(nthr))
  END IF

  GlobalBubbles = Solver % GlobalBubbles
  Stabilize = GetStabilizeFlag( Solver % Values, GotIt )

  ! Only LocalMatrixScalar's legacy per-node branch uses this -- see the
  ! matching BubblesDefault resolution in KESolverLegacy.F90 (which also
  ! consults "Stabilization method" as a synonym for the literal "Bubbles"
  ! keyword).
  BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  IF ( .NOT. GotIt ) BubblesDefault = GetString(GetSolverParams(Solver), &
            'Stabilization method', GotIt ) == 'bubbles'
  IF ( .NOT. GotIt ) BubblesDefault = .TRUE.

  IF ( .NOT. ALLOCATED(NodalDensity) ) THEN
    ALLOCATE( NodalDensity(Solver % Mesh % NumberOfNodes), &
        NodalViscosity(Solver % Mesh % NumberOfNodes), &
        NodalCmu(Solver % Mesh % NumberOfNodes) )
  END IF

  ! K and Epsilon are interleaved, Dofs=2 -- exactly as KESolverLegacy.F90
  ! sizes its own bx/bxprev.
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
    CALL Info(Caller,'KEpsilon iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)

    CALL DefaultInitialize()

    Active = GetNOFActive()

    IF ( UseScalarFallback ) THEN
      ! Serial only: LocalMatrixScalar relies on the classic argument-less
      ! Get*() "current element" accessors, exactly as KESolverLegacy.F90's
      ! own driver does. It also needs Clip and the nodal Density/Viscosity/
      ! Cmu arrays (for the positivity clip below) recorded per element --
      ! mirrored here from its own Material lookups.
      Clip = 1.0d-6
      DO i=1,Active
        Element => GetActiveElement(i)
        CALL LocalMatrixScalar( Element, dt, TransientSimulation, GlobalBubbles, BubblesDefault )
        Material => GetMaterial( Element )
        n = GetElementNOFNodes( Element )
        Clip = GetConstReal( Material, 'KE Clip', GotIt )
        IF ( .NOT. GotIt ) Clip = 1.0d-6
        Density(1:n)   = GetReal( Material, 'Density', UElement=Element )
        Viscosity(1:n) = GetReal( Material, 'Viscosity', UElement=Element )
        NodalDensity( Element % NodeIndexes(1:n) )   = Density(1:n)
        NodalViscosity( Element % NodeIndexes(1:n) ) = Viscosity(1:n)
        NodalCmu( Element % NodeIndexes(1:n) )       = GetCmu( Material )
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
      Clip = RecordNodalMatProps( Element, n )

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

      ! Nodal Density/Viscosity/Cmu for the positivity clip below, and Clip
      ! itself: a second, serial, read-only sweep (cheap: one GetReal/element)
      ! -- see RecordNodalMatProps.
      DO i=2,Active
        Element => GetActiveElement(i)
        n = GetElementNOFNodes( Element )
        Clip = RecordNodalMatProps( Element, n )
      END DO
    END IF

    CALL DefaultFinishBulkAssembly()

    IF( ListGetLogicalAnyBC(Model,'Epsilon Wall BC') .OR. &
        ListGetLogicalAnyBC(Model,'Noslip Wall BC') ) THEN
      DO i=1,Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(i)
        n = GetElementNOFNodes( Element )

        BC => GetBC( Element )
        IF ( ASSOCIATED( BC ) ) THEN
          IF ( ListGetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
               ListGetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
            DO j=1,n
              k = KinPerm(Element % NodeIndexes(j))
              CALL ZeroRow( Solver % Matrix, 2*k )
              Solver % Matrix % RHS(2*k) = 0.0_dp
            END DO
          END IF
        END IF
      END DO
    END IF

    DO i=1,Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(i)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE

      n = GetElementNOFNodes( Element )
      BC => GetBC( Element )

      IF ( ASSOCIATED( BC ) ) THEN
        IF ( ListGetLogical( BC, 'Wall Law',gotIt ) ) THEN

          Density(1:n)   = GetParentMatProp( 'Density', Element )
          Viscosity(1:n) = GetParentMatProp( 'Viscosity', Element )

          SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness', gotIt )
          LayerThickness(1:n)   = GetReal( BC, 'Boundary Layer Thickness' )

          DO j=1,n
            k = FlowPerm(Element % NodeIndexes(j))
            IF ( k > 0 ) THEN
              SELECT CASE( NSDOFs )
                CASE(3)
                  U(j) = FlowSolution( NSDOFs*k-2 )
                  V(j) = FlowSolution( NSDOFs*k-1 )
                  W(j) = 0.0D0

                CASE(4)
                  U(j) = FlowSolution( NSDOFs*k-3 )
                  V(j) = FlowSolution( NSDOFs*k-2 )
                  W(j) = FlowSolution( NSDOFs*k-1 )
              END SELECT
            ELSE
              U(j) = 0.0d0
              V(j) = 0.0d0
              W(j) = 0.0d0
            END IF
          END DO

          DO j=1,n
            CALL KEWall( Work(1), Work(2), Work(3), SQRT(U(j)**2+V(j)**2+W(j)**2), &
             LayerThickness(j), SurfaceRoughness(j), Viscosity(j), &
               Density(j) )

            k = 2*(KinPerm(Element % NodeIndexes(j))-1)

            CALL UpdateDirichletDof( Solver % Matrix, k+1, Work(1) )
            CALL UpdateDirichletDof( Solver % Matrix, k+2, Work(2) )
          END DO
        END IF

        IF ( ListGetLogical( BC, 'Epsilon Wall BC', gotIt ) .OR. &
             ListGetLogical( BC, 'Noslip Wall BC',  gotIt ) ) THEN
          CALL EpsilonWall( Element, n, Solver, TransientSimulation )
        END IF
      END IF
    END DO

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    ! K positive; epsilon floored, using the viscous-sublayer bound
    ! eps_min = rho*Cmu*k^2/mu with THIS node's own material properties --
    ! exactly as KESolverLegacy.F90's own positivity clip.
    n = Solver % Mesh % NumberOfNodes
    Kmax = MAXVAL( Solver % Variable % Values(1::2) )
    Emax = MAXVAL( Solver % Variable % Values(2::2) )
    DO i=1,n
      k = Solver % Variable % Perm(i)
      IF ( k <= 0 ) CYCLE

      KVal = Solver % Variable % Values(2*k-1)
      EVal = Solver % Variable % Values(2*k-0)

      IF ( KVal < Clip*Kmax ) KVal = Clip*Kmax

      IF ( EVal < Clip*Emax ) THEN
        EVal = MAX( NodalDensity(i)*NodalCmu(i)*KVal**2/NodalViscosity(i), Clip*Emax )
      END IF

      Solver % Variable % Values(2*k-1) = MAX( KVal, 1.0d-10 )
      Solver % Variable % Values(2*k-0) = MAX( EVal, 1.0d-10 )
    END DO

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  n = SIZE( Solver % Variable % Values )
  KE => VariableGet( Solver % Mesh % Variables, 'Kinetic Energy' )
  IF (ASSOCIATED(KE)) KE % Values = Solver % Variable % Values(1:n:2)

  KE => VariableGet( Solver % Mesh % Variables, 'Kinetic Dissipation' )
  IF (ASSOCIATED(KE)) KE % Values = Solver % Variable % Values(2:n:2)

CONTAINS

!------------------------------------------------------------------------------
!> Records this element's own Density/Viscosity/Cmu into the nodal arrays the
!> positivity clip above reads, and returns "KE Clip" -- same GetReal/
!> GetConstReal calls LocalMatrixVec's own handles resolve, just read again
!> here (cheap, once per element, serially) since the positivity clip needs
!> them per NODE, not per Gauss point.
!------------------------------------------------------------------------------
  FUNCTION RecordNodalMatProps( Element, n ) RESULT( Clip )
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    REAL(KIND=dp) :: Clip
    TYPE(ValueList_t), POINTER :: Mat
    LOGICAL :: Found
    REAL(KIND=dp) :: DensityN(n), ViscosityN(n)

    Mat => GetMaterial( Element )
    Clip = GetConstReal( Mat, 'KE Clip', Found )
    IF ( .NOT. Found ) Clip = 1.0d-6

    DensityN(1:n)   = GetReal( Mat, 'Density', UElement=Element )
    ViscosityN(1:n) = GetReal( Mat, 'Viscosity', UElement=Element )
    NodalDensity( Element % NodeIndexes(1:n) )   = DensityN(1:n)
    NodalViscosity( Element % NodeIndexes(1:n) ) = ViscosityN(1:n)
    NodalCmu( Element % NodeIndexes(1:n) )       = GetCmu( Mat )
  END FUNCTION RecordNodalMatProps

!------------------------------------------------------------------------------
!> "KE Cmu", with the same KE-Model-dependent default LocalMatrixVec's own
!> Cmu_h handle falls back to.
!------------------------------------------------------------------------------
  FUNCTION GetCmu( Mat ) RESULT( Cmu )
    TYPE(ValueList_t), POINTER :: Mat
    REAL(KIND=dp) :: Cmu
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: Model

    Cmu = ListGetConstReal( Mat, 'KE Cmu', Found )
    IF ( Found ) RETURN

    CALL GetStringThreadSafe( Mat, 'KE Model', Model, Found )
    IF ( .NOT. Found ) Model = 'standard'
    SELECT CASE( Model )
    CASE( 'v2-f' )
      Cmu = 0.22_dp
    CASE( 'rng' )
      Cmu = 0.0845_dp
    CASE DEFAULT
      Cmu = 0.09_dp
    END SELECT
  END FUNCTION GetCmu
!------------------------------------------------------------------------------
END SUBROUTINE KESolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Element-type resolution, called before the mesh/basis functions are
!> finalized -- same role and timing as HeatSolver_Init0 and
!> Spalart-Allmaras.F90's/Komega.F90's/SSTKomega.F90's own _Init0.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KESolver_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KESolverFront
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
END SUBROUTINE KESolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: KESolver. Under "Legacy Assembly",
!> delegates to KESolverLegacy_Init and returns -- the "Variable"/p-bubble
!> setup below is specific to this file's own implementation (identical to
!> KESolverLegacy.F90's own, since that part was never scalar-specific).
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE KESolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE KESolverFront
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
    CALL DelegateToKESolverLegacy( 'KESolverLegacy_Init', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  SolverParams => GetSolverParams()

  str = GetString( SolverParams,'Variable', Found )
  IF ( .NOT. Found ) str = 'K-eps'
  IF ( INDEX( str, '[' ) <= 0 ) THEN
    CALL ListAddString( SolverParams, 'Variable', &
          TRIM(str) // '[Kinetic Energy:1 Kinetic Dissipation:1]' )
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
  ELSE IF ( TransientSimulation ) THEN
    LegacyBubbles = ListGetLogical( SolverParams, 'Bubbles', Found )
    IF ( .NOT. Found ) LegacyBubbles = ListGetString( SolverParams, &
        'Stabilization method', Found ) == 'bubbles'
    IF ( .NOT. Found ) LegacyBubbles = .TRUE.

    IF ( LegacyBubbles ) THEN
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Min Iterations', 2)
      CALL ListAddNewInteger(SolverParams, 'Nonlinear System Max Iterations', 2)
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE KESolver_Init
!------------------------------------------------------------------------------
