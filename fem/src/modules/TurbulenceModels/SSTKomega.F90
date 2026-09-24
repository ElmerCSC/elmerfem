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
!> Solver for the (RC)SST K-omega turbulence model. Vectorized/threaded
!> implementation following the same pattern as Komega.F90 and
!> Spalart-Allmaras.F90. Like plain K-omega, the SST blending only adds
!> per-Gauss-point nonlinear coefficients (F1,F2,CD,Tmu,SigmaK,SigmaO,Beta,
!> rGamma) and, for the omega equation only, an extra "effective velocity"
!> cross-diffusion correction (-2*rho*(1-F1)/1.168/Omega*GradK, folded into
!> the same rho*EffVelo.grad(Basis) form the plain convection term already
!> has) -- see A(1,1)/A(2,2) in the legacy LocalMatrix: there is still no
!> genuine A(1,2)/A(2,1) matrix coupling. So this assembles K and omega as
!> two independent scalar problems, each with its own EffVelo, and
!> interleaves them right before condensation/glue, exactly as Komega.F90
!> does.
!>
!> The buoyancy term rho_g = grad(rho).Gravity is exactly zero here: the
!> vectorized path only supports a spatially uniform Density (see
!> UseScalarFallback below), whose gradient is zero regardless of Gravity.
!> A case that needs a real (non-uniform) density -- any "Compressibility
!> Model" other than the default -- routes to LocalMatrixScalar instead,
!> which computes it properly via ElementDensity, exactly as legacy does.
!>
!> Axisymmetric/cylindrical coordinates go through LocalMatrixScalar instead
!> of LocalMatrixVec -- a scalar, per-Gauss-point fallback carrying the same
!> metric-tensor math as SSTKomegaLegacy.F90's own LocalMatrix (already
!> interleaved directly via STIFF(2*(p-1)+i,2*(q-1)+j), so ported verbatim),
!> called serially -- same role as HeatSolve.F90's own AxiSymmetric branch.
!> LocalMatrixScalar also carries the legacy "Bubbles = True" per-node scheme
!> and any non-default Compressibility Model. OmegaWall (near-wall Dirichlet
!> omega) and KomegaWallLaw (the "Wall Law" wall-function BC, via the core
!> KEWall routine in Walls.F90) are boundary treatments shared by both bulk
!> paths, called once from the driver's own boundary loop -- ported unchanged
!> bar taking Solver/Model explicitly, since they are module procedures here,
!> not nested inside the driver.
!>
!> The original scalar-element solver lives on in SSTKomegaLegacy.F90
!> (subroutine SSTKOmegaLegacy), reachable either directly by that name or
!> via "Legacy Assembly = Logical True" here (see SSTKOmegaFront below). With
!> axisymmetric, non-default compressibility and the per-node scheme all
!> covered above, that path is now only needed by a sif that must reproduce
!> the legacy solver's exact historical numbers.
!> \ingroup Solvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Whether this solver should run the original scalar-element assembly
!> (SSTKomegaLegacy.F90) instead of this file's own implementation.
!------------------------------------------------------------------------------
MODULE SSTKOmegaFront
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
!> Call one of SSTKOmegaLegacy's entry points with this solver. The name is
!> resolved at run time, as the core resolves any solver, so this file and
!> SSTKomegaLegacy.so stay independent of one another.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateToSSTKOmegaLegacy( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'SSTKomegaLegacy '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'SSTKOmega', &
        '"Legacy Assembly" was requested but "'//TRIM(Entry)//'" could not be found. '// &
        'Is SSTKomegaLegacy.so installed beside this solver?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateToSSTKOmegaLegacy

END MODULE SSTKOmegaFront


!------------------------------------------------------------------------------
MODULE SSTKOmegaLocalForms

  USE DefUtils
  USE LinearForms

  IMPLICIT NONE

  ! Per-element bubble history, used by Default1stOrderTime's Nb path
  ! (DefUtils.F90) -- lives on Solver % Variable's own BubbleValues/
  ! BubblePrevValues (Types.F90), not a separate type; see the matching
  ! bx/bxprev comment in SSTKomegaLegacy.F90, which this mirrors.

  ! Per-thread ValueHandle_t storage for LocalMatrixVec's material lookups.
  ! NOT THREADPRIVATE -- see the matching comment on IncompressibleNS.F90's
  ! NSHandles_t for the Windows/GCC emutls hazard that rules that out.
  TYPE :: SSTHandles_t
    TYPE(ValueHandle_t) :: Visc_h, Dens_h, PrRho_h, C3Omega_h
  END TYPE SSTHandles_t
  TYPE(SSTHandles_t), ALLOCATABLE, SAVE :: SSTHandles(:)

CONTAINS

!------------------------------------------------------------------------------
!> Assemble and glue local matrix/RHS for one bulk element. Vectorized over
!> Gauss points, safe to call concurrently from multiple threads (each thread
!> passes its own InitHandles and only touches SSTHandles(tid)).
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
    TYPE(ValueList_t), POINTER :: Material

    REAL(KIND=dp), ALLOCATABLE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:)
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), FORCE(:), TimeForce(:)
    REAL(KIND=dp), ALLOCATABLE :: MassK(:,:), StiffK(:,:), ForceK(:), &
        MassO(:,:), StiffO(:,:), ForceO(:)

    REAL(KIND=dp), POINTER :: RhoVec(:), MuVec(:), PrVec(:), C3Vec(:)

    REAL(KIND=dp), ALLOCATABLE :: VeloNodal(:,:), KNodal(:), ONodal(:), &
        DistNodal(:), PressureNodal(:)
    REAL(KIND=dp), ALLOCATABLE :: VeloVec(:,:), dVelodxVec(:,:,:), KVec(:), OVec(:), &
        GradKVec(:,:), GradOVec(:,:), DistVec(:), PressureVec(:), &
        StrainVec(:,:,:), VorticityVec(:,:,:), StrainMeasureVec(:), &
        VorticityMeasureVec(:), StrainDotGradUVec(:), DivVelVec(:), &
        CDVec(:), F1Vec(:), F2Vec(:), BetaVec(:), SigmaKVec(:), SigmaOVec(:), &
        rGammaVec(:), TmuVec(:), Effmu1Vec(:), Effmu2Vec(:), SoundSpeedSqVec(:), &
        MachSqVec(:), ProdKVec(:), ProdOVec(:), ReactK(:), ReactO(:), &
        LoadK(:), LoadO(:), EffVeloO(:,:), &
        StreamVecK(:,:), StreamVecO(:,:), TauK(:), TauO(:), TmpVec(:), RadiusVec(:)

    REAL(KIND=dp) :: hK, mK, VNorm, SpecificHeatRatio, ReferencePressure
    INTEGER :: i,j,k,p,ngp,dim,allocstat,tid,ntot
    LOGICAL :: Stat, Found, IsAxiSymmetric
!------------------------------------------------------------------------------
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Visc_h => SSTHandles(tid) % Visc_h, Dens_h => SSTHandles(tid) % Dens_h, &
               PrRho_h => SSTHandles(tid) % PrRho_h, C3Omega_h => SSTHandles(tid) % C3Omega_h )

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Visc_h,'Material','Viscosity' )
      CALL ListInitElementKeyword( Dens_h,'Material','Density' )
      CALL ListInitElementKeyword( PrRho_h,'Material','Turbulent Prandtl Number' )
      CALL ListInitElementKeyword( C3Omega_h,'Material','Dissipation buoyancy coefficient' )
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ntot = nd + nb

    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ALLOCATE( BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        MassK(ntot,ntot), StiffK(ntot,ntot), ForceK(ntot), &
        MassO(ntot,ntot), StiffO(ntot,ntot), ForceO(ntot), &
        MASS(2*ntot,2*ntot), STIFF(2*ntot,2*ntot), FORCE(2*ntot), TimeForce(2*ntot), &
        VeloNodal(3,n), KNodal(n), ONodal(n), DistNodal(n), PressureNodal(n), &
        VeloVec(ngp,3), dVelodxVec(ngp,3,3), KVec(ngp), OVec(ngp), &
        GradKVec(ngp,3), GradOVec(ngp,3), DistVec(ngp), PressureVec(ngp), &
        StrainVec(ngp,3,3), VorticityVec(ngp,3,3), StrainMeasureVec(ngp), &
        VorticityMeasureVec(ngp), StrainDotGradUVec(ngp), DivVelVec(ngp), &
        CDVec(ngp), F1Vec(ngp), F2Vec(ngp), BetaVec(ngp), SigmaKVec(ngp), SigmaOVec(ngp), &
        rGammaVec(ngp), TmuVec(ngp), Effmu1Vec(ngp), Effmu2Vec(ngp), SoundSpeedSqVec(ngp), &
        MachSqVec(ngp), ProdKVec(ngp), ProdOVec(ngp), ReactK(ngp), ReactO(ngp), &
        LoadK(ngp), LoadO(ngp), EffVeloO(ngp,3), &
        StreamVecK(ngp,ntot), StreamVecO(ngp,ntot), TauK(ngp), TauO(ngp), TmpVec(ngp), &
        RadiusVec(ngp), STAT=allocstat )
    IF( allocstat /= 0 ) CALL Fatal('SSTKomega','Local storage allocation failed')

    CALL GetElementNodesVec( Nodes, UElement=Element )

    MassK = 0._dp; StiffK = 0._dp; ForceK = 0._dp
    MassO = 0._dp; StiffO = 0._dp; ForceO = 0._dp

    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, DetJVec, &
        SIZE(BasisVec,2), BasisVec, dBasisdxVec )
    DetJVec(1:ngp) = DetJVec(1:ngp) * IP % s(1:ngp)

    ! Axisymmetric (no swirl): r-weighted measure, plus the hoop strain and
    ! hoop divergence corrections further down -- see the matching (more
    ! detailed) comments in Spalart-Allmaras.F90's own LocalMatrixVec. Genuine
    ! swirl ("Cylindric Symmetric") still goes through LocalMatrixScalar.
    IsAxiSymmetric = ( CurrentCoordinateSystem() == AxisSymmetric )
    IF( IsAxiSymmetric ) THEN
      RadiusVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), Nodes % x(1:n) )
      DetJVec(1:ngp) = DetJVec(1:ngp) * RadiusVec(1:ngp)
    END IF

    VeloNodal = 0._dp
    CALL GetScalarLocalSolution( VeloNodal(1,1:n), 'Velocity 1', UElement=Element )
    CALL GetScalarLocalSolution( VeloNodal(2,1:n), 'Velocity 2', UElement=Element )
    IF( dim == 3 ) CALL GetScalarLocalSolution( VeloNodal(3,1:n), 'Velocity 3', UElement=Element )

    CALL GetScalarLocalSolution( KNodal, 'Kinetic energy', UElement=Element )
    CALL GetScalarLocalSolution( ONodal, 'Kinetic Dissipation', UElement=Element )
    CALL GetScalarLocalSolution( DistNodal, 'Wall Distance', UElement=Element )
    CALL GetScalarLocalSolution( PressureNodal, 'Pressure', UElement=Element )

    RhoVec => ListGetElementRealVec( Dens_h, ngp, BasisVec, Element, Found )
    MuVec  => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element, Found )

    PrVec => ListGetElementRealVec( PrRho_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) PrVec = 0.85_dp

    C3Vec => ListGetElementRealVec( C3Omega_h, ngp, BasisVec, Element, Found )
    IF( .NOT. Found ) C3Vec = 0.0_dp

    ! Per-material constants (not spatially varying, unlike the handles
    ! above) -- same GetCReal calls and same (Elmer default) 0.0 fallback as
    ! the legacy LocalMatrix, so Mach_number_sq below comes out zero exactly
    ! as it does there whenever "Specific Heat Ratio" is left unset.
    Material => GetMaterial( Element )
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

    GradKVec = 0._dp; GradOVec = 0._dp
    DO i=1,dim
      GradKVec(1:ngp,i) = MATMUL( dBasisdxVec(1:ngp,1:n,i), KNodal(1:n) )
      GradOVec(1:ngp,i) = MATMUL( dBasisdxVec(1:ngp,1:n,i), ONodal(1:n) )
    END DO

    DistVec(1:ngp) = MAX( MATMUL( BasisVec(1:ngp,1:n), DistNodal(1:n) ), 1.0d-10 )
    PressureVec(1:ngp) = MATMUL( BasisVec(1:ngp,1:n), PressureNodal(1:n) )

    ! Strain-rate / vorticity tensors and their (Frobenius) measures, Cartesian
    ! only. StrainMeasure floored at 1e-10, exactly as legacy.
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
    ! in MaterialModels.F90, and the matching comment in Spalart-Allmaras.F90).
    ! Vorticity needs no such addition (zero diagonal by antisymmetry).
    IF( IsAxiSymmetric ) THEN
      StrainMeasureVec(1:ngp) = StrainMeasureVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF

    StrainMeasureVec(1:ngp)    = MAX( SQRT( 2._dp*StrainMeasureVec(1:ngp) ), 1.0d-10 )
    VorticityMeasureVec(1:ngp) = SQRT( 2._dp*VorticityMeasureVec(1:ngp) )

    ! SUM(Strain*dVelodx) and trace(dVelodx) = div(Velo), shared by ProdK/ProdO.
    StrainDotGradUVec = 0._dp
    DO i=1,dim
      DO j=1,dim
        StrainDotGradUVec(1:ngp) = StrainDotGradUVec(1:ngp) + StrainVec(1:ngp,i,j)*dVelodxVec(1:ngp,i,j)
      END DO
    END DO
    IF( IsAxiSymmetric ) THEN
      StrainDotGradUVec(1:ngp) = StrainDotGradUVec(1:ngp) + ( VeloVec(1:ngp,1) / RadiusVec(1:ngp) )**2
    END IF

    DivVelVec = 0._dp
    DO i=1,dim
      DivVelVec(1:ngp) = DivVelVec(1:ngp) + dVelodxVec(1:ngp,i,i)
    END DO
    ! The true divergence in axisymmetric (no swirl) coordinates is
    ! du_r/dr + u_r/r + du_z/dz -- trace(dVelodx) above is only the first and
    ! third terms.
    IF( IsAxiSymmetric ) THEN
      DivVelVec(1:ngp) = DivVelVec(1:ngp) + VeloVec(1:ngp,1) / RadiusVec(1:ngp)
    END IF

    ! Cross-diffusion CD and the F1/F2 blending functions.
    CDVec = 0._dp
    DO i=1,dim
      CDVec(1:ngp) = CDVec(1:ngp) + GradKVec(1:ngp,i)*GradOVec(1:ngp,i)
    END DO
    CDVec(1:ngp) = MAX( 2._dp*RhoVec(1:ngp)/1.168_dp/OVec(1:ngp)*CDVec(1:ngp), 1.0d-10 )

    F1Vec(1:ngp) = SQRT(KVec(1:ngp)) / 0.09_dp / OVec(1:ngp) / DistVec(1:ngp)
    F1Vec(1:ngp) = MAX( F1Vec(1:ngp), 500._dp*MuVec(1:ngp)/RhoVec(1:ngp)/OVec(1:ngp)/DistVec(1:ngp)**2 )
    F1Vec(1:ngp) = MIN( F1Vec(1:ngp), 4._dp*RhoVec(1:ngp)/1.168_dp*KVec(1:ngp)/CDVec(1:ngp)/DistVec(1:ngp)**2 )
    F1Vec(1:ngp) = TANH( F1Vec(1:ngp)**4 )

    F2Vec(1:ngp) = MAX( 2._dp*SQRT(KVec(1:ngp))/0.09_dp/OVec(1:ngp)/DistVec(1:ngp), &
        500._dp*MuVec(1:ngp)/RhoVec(1:ngp)/OVec(1:ngp)/DistVec(1:ngp)**2 )
    F2Vec(1:ngp) = TANH( F2Vec(1:ngp)**2 )
    ! F3 (rough-wall correction) and F4 (rotation/curvature correction) are
    ! "NOT IN USE" in the legacy LocalMatrix too -- both fixed at 1.

    BetaVec(1:ngp)   = 0.075_dp*F1Vec(1:ngp) + 0.0828_dp*(1._dp-F1Vec(1:ngp))
    SigmaKVec(1:ngp) = 1.176_dp*F1Vec(1:ngp) + 1.0000_dp*(1._dp-F1Vec(1:ngp))
    SigmaOVec(1:ngp) = 2.000_dp*F1Vec(1:ngp) + 1.1680_dp*(1._dp-F1Vec(1:ngp))
    rGammaVec(1:ngp) = (5._dp/9._dp)*F1Vec(1:ngp) + 0.44_dp*(1._dp-F1Vec(1:ngp))

    TmuVec(1:ngp) = 0.31_dp*RhoVec(1:ngp)*KVec(1:ngp) / &
        MAX( 0.31_dp*OVec(1:ngp), VorticityMeasureVec(1:ngp)*F2Vec(1:ngp) )
    Effmu1Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaKVec(1:ngp)
    Effmu2Vec(1:ngp) = MuVec(1:ngp) + TmuVec(1:ngp)/SigmaOVec(1:ngp)

    ! Buoyancy (rho_g) is exactly zero here -- see the file header. Mach
    ! number correction kept for fidelity, evaluates to zero whenever
    ! "Specific Heat Ratio" is unset, exactly as legacy.
    SoundSpeedSqVec(1:ngp) = (PressureVec(1:ngp)+ReferencePressure) * SpecificHeatRatio / RhoVec(1:ngp)
    MachSqVec = 0._dp
    WHERE( SoundSpeedSqVec(1:ngp) > 0._dp ) MachSqVec(1:ngp) = KVec(1:ngp) / SoundSpeedSqVec(1:ngp)

    ProdKVec(1:ngp) = 2._dp*TmuVec(1:ngp)*StrainDotGradUVec(1:ngp) &
        - (2._dp/3._dp)*RhoVec(1:ngp)*KVec(1:ngp)*DivVelVec(1:ngp) &
        - 2._dp*RhoVec(1:ngp)*0.09_dp*KVec(1:ngp)*OVec(1:ngp)*MachSqVec(1:ngp)

    ProdOVec(1:ngp) = 2._dp*RhoVec(1:ngp)*StrainDotGradUVec(1:ngp) &
        - (2._dp/3._dp)*RhoVec(1:ngp)*OVec(1:ngp)*DivVelVec(1:ngp)

    ReactK(1:ngp) = RhoVec(1:ngp) * 0.09_dp * OVec(1:ngp)
    ReactO(1:ngp) = RhoVec(1:ngp) * BetaVec(1:ngp) * OVec(1:ngp)

    LoadK(1:ngp) = ProdKVec(1:ngp)
    LoadO(1:ngp) = rGammaVec(1:ngp) * ProdOVec(1:ngp)

    ! K convects with the plain velocity; omega additionally carries the SST
    ! cross-diffusion correction, folded into an effective velocity exactly
    ! as Spalart-Allmaras folds its own Cb2 correction.
    EffVeloO = 0._dp
    DO i=1,dim
      EffVeloO(1:ngp,i) = VeloVec(1:ngp,i) - &
          2._dp*(1._dp-F1Vec(1:ngp))/1.168_dp*GradKVec(1:ngp,i)/OVec(1:ngp)
    END DO

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
        RhoVec, EffVeloO )
    CALL LinearForms_UdotF( ngp, ntot, BasisVec, DetJVec, LoadO, ForceO )

    !------------------------------------------------------------------------
    ! SUPG (equal-order) stabilization, opt-in via "Stabilize"/"Stabilization
    ! Method" -- same Franca et al. tau/streamline construction as
    ! Spalart-Allmaras.F90's and Komega.F90's own SUPG blocks. Unlike plain
    ! Komega, K and omega convect with different effective velocities here,
    ! so each equation gets its own streamline-weighted test function.
    !------------------------------------------------------------------------
    IF( Stabilize ) THEN
      hK = Element % hK
      mK = Element % StabilizationMK

      StreamVecK(1:ngp,1:ntot) = 0._dp
      StreamVecO(1:ngp,1:ntot) = 0._dp
      DO i=1,dim
        DO p=1,ntot
          StreamVecK(1:ngp,p) = StreamVecK(1:ngp,p) + &
              RhoVec(1:ngp) * VeloVec(1:ngp,i) * dBasisdxVec(1:ngp,p,i)
          StreamVecO(1:ngp,p) = StreamVecO(1:ngp,p) + &
              RhoVec(1:ngp) * EffVeloO(1:ngp,i) * dBasisdxVec(1:ngp,p,i)
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

        VNorm = SQRT( SUM( EffVeloO(j,1:dim)**2 ) )
        IF( VNorm > 0._dp .AND. Effmu2Vec(j) /= 0._dp ) THEN
          TauO(j) = MIN( 1._dp, mK*hK*RhoVec(j)*VNorm / (2._dp*ABS(Effmu2Vec(j))) )
          TauO(j) = hK * TauO(j) / ( 2._dp * RhoVec(j) * VNorm )
        ELSE
          TauO(j) = 0._dp
        END IF
      END DO

      TmpVec(1:ngp) = TauK(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVecK, StreamVecK, TmpVec, StiffK )
      CALL LinearForms_UdotF( ngp, ntot, StreamVecK, TmpVec, LoadK, ForceK )

      TmpVec(1:ngp) = TauO(1:ngp) * DetJVec(1:ngp)
      CALL LinearForms_UdotV( ngp, ntot, dim, StreamVecO, StreamVecO, TmpVec, StiffO )
      CALL LinearForms_UdotF( ngp, ntot, StreamVecO, TmpVec, LoadO, ForceO )

      IF( Transient ) THEN
        TmpVec(1:ngp) = TauK(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVecK, BasisVec, TmpVec, MassK )
        TmpVec(1:ngp) = TauO(1:ngp) * RhoVec(1:ngp) * DetJVec(1:ngp)
        CALL LinearForms_UdotV( ngp, ntot, dim, StreamVecO, BasisVec, TmpVec, MassO )
      END IF
    END IF

    ! Interleave the two independent scalar blocks -- odd rows/columns for K,
    ! even for omega, off-diagonal (K,omega) blocks stay zero, matching
    ! legacy's own A(1,2)/A(2,1) (never written).
    MASS = 0._dp; STIFF = 0._dp; FORCE = 0._dp
    MASS(1:2*ntot-1:2,1:2*ntot-1:2) = MassK
    MASS(2:2*ntot:2,  2:2*ntot:2)   = MassO
    STIFF(1:2*ntot-1:2,1:2*ntot-1:2) = StiffK
    STIFF(2:2*ntot:2,  2:2*ntot:2)   = StiffO
    FORCE(1:2*ntot-1:2) = ForceK
    FORCE(2:2*ntot:2)   = ForceO

    !------------------------------------------------------------------------
    ! Time discretization and p-bubble condensation -- mirrors the nb>0
    ! branch of SSTKomegaLegacy.F90's driver exactly (DOFs=2, K/omega
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
!> do: axisymmetric/cylindrical coordinates and any non-default
!> "Compressibility Model" (the ElementKernel below carries the same
!> metric-tensor and ElementDensity calls as SSTKomegaLegacy.F90's own
!> LocalMatrix -- verbatim, no new math, and already interleaves K/omega
!> directly via STIFF(2*(p-1)+i,2*(q-1)+j)), and the legacy "Bubbles = True"
!> per-node scheme for a plain nodal "Element" set explicitly. Always called
!> serially (see the UseScalarFallback branch in SSTKOmega below), so it uses
!> the classic GetMaterial()/GetReal()/argument-less GetElementNOF*()
!> accessors exactly as the legacy driver does -- not safe to call from
!> inside an OMP parallel region, unlike LocalMatrixVec.
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
    IF( allocstat /= 0 ) CALL Fatal('SSTKomega','Local storage allocation failed')

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
!> Verbatim port of SSTKomegaLegacy.F90's nested LocalMatrix. "Bubbles" and
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

      REAL(KIND=dp) :: A(2,2),M(2,2),Prod,div,ProdK,ProdO,ProdTensor(3,3),Ident(3,3)
      INTEGER :: i,j,c,p,q,t,dim,NBasis
      REAL(KIND=dp) :: LoadatIp(2),Cmu,Rho,mu,Tmu,Effmu(2)

      REAL(KIND=dp) :: s,u,v,w, K,Omega,Strain(3,3), Vorticity(3,3), dist, &
             Mach_number_sq, Sound_speed_sq

      REAL(KIND=dp) :: StrainMeasure,VorticityMeasure,X,Y,Z,SigmaK,  &
              SigmaO,Beta,CD,F1,F2,F3,F4,rGamma, GradK(3), GradO(3), &
              Gravity(3), rho_g, Pr_rho(n), Pr, c3_omega(n), c3, Pressure(n), &
              ReferencePressure, SpecificHeatRatio

      REAL(KIND=dp), POINTER :: gWork(:,:)

      REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),SqrtMetric

      LOGICAL :: stat, GotIt
      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()

      ! Model isn't in scope here by host association (ElementKernel is
      ! nested inside LocalMatrixScalar, a module procedure, not inside the
      ! driver the way the legacy nested LocalMatrix was) -- CurrentModel is
      ! the same model regardless.
      gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',GotIt)
      IF ( GotIt ) THEN
        Gravity = gWork(1:3,1)*gWork(4,1)
      ELSE
        Gravity    =  0.00_dp
        Gravity(2) = -9.81_dp
      END IF

      Viscosity(1:n) = GetReal( Material, 'Viscosity' )
      CALL ElementDensity( Density, n )

      Pr_rho(1:n) = GetReal( Material, 'Turbulent Prandtl Number', stat )
      IF ( .NOT. stat ) Pr_rho(1:n) = 0.85_dp

      c3_omega(1:n) = GetReal( Material, 'Dissipation buoyancy coefficient', stat )
      IF ( .NOT. stat ) c3_omega(1:n) = 0.0_dp

      SpecificheatRatio = GetCReal( Material, 'Specific Heat Ratio', stat )
      CALL getScalarLocalSolution( Pressure, 'Pressure' )
      ReferencePressure = GetCReal( Material, 'Reference Pressure', stat )

      CALL GetScalarLocalSolution( UX, 'Velocity 1' )
      CALL GetScalarLocalSolution( UY, 'Velocity 2' )
      CALL GetScalarLocalSolution( UZ, 'Velocity 3' )

      CALL GetScalarLocalSolution( Energy, 'Kinetic energy' )
      CALL GetScalarLocalSolution( Distance, 'Wall Distance' )
      CALL GetScalarLocalSolution( Dissipation, 'Kinetic Dissipation' )

      FORCE = 0.0D0
      STIFF = 0.0D0
      MASS  = 0.0D0

      NBasis = nd

      Ident = 0._dp
      DO i=1,3
        Ident(i,i) = 1._dp
      END DO

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
          StrainMeasure = MAX( SQRT(2*SUM(Strain*Strain)), 1.0d-10 )

          Vorticity = 0.5_dp * (dVelodx - TRANSPOSE(dVelodx))
          VorticityMeasure = SQRT(2*SUM(Vorticity*Vorticity))
        ELSE
          StrainMeasure = SQRT(SecondInvariant( Velo,dVelodx,Metric,Symb )/2)
        END IF

        K = SUM( Energy(1:n) * Basis(1:n) )
        Omega = SUM( Dissipation(1:n) * Basis(1:n) )

        DO i=1,dim
          GradK(i) = SUM( dBasisdx(1:n,i) * Energy(1:n) )
          GradO(i) = SUM( dBasisdx(1:n,i) * Dissipation(1:n) )
        END DO

        mu   = SUM( Viscosity(1:n) * Basis(1:n) )
        rho  = SUM( Density(1:n) * Basis(1:n) )
        rho_g = 0._dp
        DO i=1,dim
          rho_g = rho_g + SUM(Density(1:n) * dBasisdx(1:n,i)) * Gravity(i)
        END DO
        dist = MAX( SUM( Distance(1:n) * Basis(1:n) ), 1.0d-10 )

        Sound_speed_sq = (SUM(Basis(1:n)*Pressure(1:n))+ReferencePressure) * &
                      SpecificHeatRatio / rho

        Mach_number_sq = 0._dp
        IF ( Sound_speed_sq > 0._dp ) Mach_number_sq = K / Sound_speed_sq

        Pr = SUM( Basis(1:n) * Pr_rho(1:n) )
        c3 = SUM( Basis(1:n) * c3_omega(1:n) )

        CD = MAX(2*rho/1.168_dp/Omega*SUM(GradK(1:dim)*GradO(1:dim)),1.d-10)

        F1 = SQRT(K) / 0.09_dp / Omega / dist
        F1 = MAX( F1, 500 * mu / rho / Omega / dist**2 )
        F1 = MIN( F1, 4 * rho / 1.168_dp* K / CD / dist**2 )
        F1 = TANH( F1**4 )

        F2 = MAX( 2*SQRT(K)/0.09_dp/Omega/dist,500*mu/rho/omega/dist**2 )
        F2 = TANH(F2**2)

        F3 = 1
        F4 = 1

        Beta   = 0.075_dp*F1 + 0.0828_dp*(1-F1)
        SigmaK = 1.176_dp*F1 + 1.0000_dp*(1-F1)
        SigmaO = 2.000_dp*F1 + 1.1680_dp*(1-F1)

        rGamma = 5._dp/9._dp * F1 + 0.44_dp * (1-F1)

        Tmu = 0.31_dp*rho*K/MAX(0.31_dp*Omega,VorticityMeasure*F2*F3)
        Effmu(1) = mu + Tmu / SigmaK
        Effmu(2) = mu + Tmu / SigmaO

        ProdK = SUM((2*Tmu*Strain-2/3._dp*rho*Ident*K)*dVelodx) - &
              Tmu * rho_g / (rho * Pr) - 2*rho*0.09_dp*K*Omega*Mach_number_sq

        ProdO = SUM((2*rho*Strain-2/3._dp*rho*Ident*Omega)*dVelodx) - &
              c3 * rho_g / Pr

        DO p=1,NBasis
        DO q=1,NBasis
           M = 0.0d0
           A = 0.0d0

           M(1,1) = rho * Basis(q) * Basis(p)
           M(2,2) = rho * Basis(q) * Basis(p)

           A(1,1) = A(1,1) + rho * 0.09_dp * Omega * Basis(q) * Basis(p)
           A(2,2) = A(2,2) + rho * F4 * Beta * Omega * Basis(q) * Basis(p)

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO i=1,dim
                A(1,1) = A(1,1) + Effmu(1) * dBasisdx(q,i) * dBasisdx(p,i)
                A(2,2) = A(2,2) + Effmu(2) * dBasisdx(q,i) * dBasisdx(p,i)
                A(2,2) = A(2,2) - 2*rho*(1-F1)/1.168_dp/Omega*GradK(i)*dBasisdx(q,i)*Basis(p)
              END DO
           ELSE
              DO i=1,dim
                DO j=1,dim
                   A(1,1) = A(1,1) + Metric(i,j) * Effmu(1) * &
                        dBasisdx(q,i) * dBasisdx(p,j)

                   A(2,2) = A(2,2) + Metric(i,j) * Effmu(2) * &
                        dBasisdx(q,i) * dBasisdx(p,j)

                   A(2,2) = A(2,2) - 2*rho*(1-F1)/1.168_dp/Omega*Metric(i,j)* &
                        GradK(i)*dBasisdx(q,i)*Basis(p)
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

         LoadAtIP(1) = ProdK
         LoadAtIP(2) = rGamma*ProdO

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
!> Near-wall Dirichlet omega from the local Wall Distance field, ported
!> unchanged from SSTKomegaLegacy.F90's own nested OmegaWall bar taking
!> Solver explicitly (module procedure here, not nested in the driver).
!------------------------------------------------------------------------------
  SUBROUTINE OmegaWall( Element, n, Solver )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Distance(32), omega_wall, dist, mu(32), rho(32)
    INTEGER :: i,j,np
    TYPE(Element_t), POINTER :: Parent
!------------------------------------------------------------------------------
    Parent => Element % BoundaryInfo % Left
    IF ( .NOT. ASSOCIATED(Parent) ) THEN
      Parent => Element % BoundaryInfo % Right
    ELSE
      IF ( .NOT. ALL(Solver % Variable % Perm(Parent % NodeIndexes)>0) ) &
        Parent => Element % BoundaryInfo % Right
    END IF
    IF(.NOT.ASSOCIATED(Parent))RETURN

    np = GetElementNOFNodes(Parent)

    rho(1:np)= GetReal( GetMaterial(Parent), 'Density', UElement=Parent )
    mu(1:np) = GetReal( GetMaterial(Parent), 'Viscosity', UElement=Parent )

    CALL GetScalarLocalSolution( Distance, 'Wall distance', UElement=Parent )

    omega_wall = 1.d10
    DO i=1,np
      j = Parent % NodeIndexes(i)
      IF ( Distance(i) < AEPS ) CYCLE
      IF ( ANY( j==Element % NodeIndexes(1:n) ) ) CYCLE

      omega_wall = 6*mu(i)/rho(i)/0.075_dp/Distance(i)**2

      j = 2*Solver % Variable % Perm(j)

      CALL UpdateDirichletDof( Solver % Matrix, j, omega_wall )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE OmegaWall
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> "Wall Law" wall-function BC: sets K/omega Dirichlet values from the
!> near-wall analytical profile (KEWall, the same core Walls.F90 routine
!> KESolver.F90 uses). Ported unchanged from SSTKomegaLegacy.F90's own nested
!> KomegaWallLaw bar taking Element/BC/Solver explicitly.
!------------------------------------------------------------------------------
  SUBROUTINE KomegaWallLaw( Element, n, BC, Solver )

    EXTERNAL :: KEWALL
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    LOGICAL :: GotIt
    INTEGER :: i,j,k,DOFs
    REAL(KIND=dp) :: Density(n),Viscosity(n),SurfaceRoughness(n),LayerThickness(n), &
            U(n), V(n), W(n), Kin, Eps, Omega
!------------------------------------------------------------------------------
    Density(1:n)   = GetParentMatProp( 'Density', Element )
    Viscosity(1:n) = GetParentMatProp( 'Viscosity', Element )

    SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness', gotIt )
    LayerThickness(1:n)   = GetReal( BC, 'Boundary Layer Thickness' )

    CALL GetScalarLocalSolution(U, 'Velocity 1', UElement=Element)
    CALL GetScalarLocalSolution(V, 'Velocity 2', UElement=Element)
    CALL GetScalarLocalSolution(W, 'Velocity 3', UElement=Element)

    DOFs = Solver % Variable % DOFs
    DO j=1,n
     CALL KEWall( Kin, Eps, Omega, SQRT(U(j)**2+V(j)**2+W(j)**2), &
      LayerThickness(j), SurfaceRoughness(j), Viscosity(j), &
        Density(j) )

      k = DOFs*(Solver % Variable % Perm(Element % NodeIndexes(j))-1)

      CALL UpdateDirichletDof( Solver % Matrix, k+1, Kin )
      CALL UpdateDirichletDof( Solver % Matrix, k+2, Omega )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE KomegaWallLaw
!------------------------------------------------------------------------------

END MODULE SSTKOmegaLocalForms


!------------------------------------------------------------------------------
!> Vectorized/threaded SST K-omega driver. See the file header above for what
!> this does and does not support, and "Legacy Assembly" for the fallback.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SSTKOmega( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SSTKOmegaLocalForms
  USE SSTKOmegaFront
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
  LOGICAL :: GotIt, InitHandles, GlobalBubbles, Stabilize, UseScalarFallback, BubblesDefault
  REAL(KIND=dp) :: Norm, KVal, EVal
  CHARACTER(*), PARAMETER :: Caller = 'SSTKomega'
!------------------------------------------------------------------------------
  IF ( LegacyAssembly( Solver ) ) THEN
    CALL DelegateToSSTKOmegaLegacy( 'SSTKomegaLegacy', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  KE => Solver % Variable
  IF ( .NOT. ASSOCIATED(KE) ) RETURN
  IF ( COUNT( KE % Perm > 0 ) <= 0 ) RETURN

  ! LocalMatrixVec now carries the plain "Axi Symmetric" (no swirl) case
  ! itself; it still only supports a spatially uniform Density (its buoyancy
  ! term is exactly zero, see the file header), so genuine swirl ("Cylindric
  ! Symmetric"), general "Cylindric", and any non-default "Compressibility
  ! Model" still go through the scalar LocalMatrixScalar fallback, serially
  ! -- same branch HeatSolve.F90 makes, and the same treatment
  ! Spalart-Allmaras.F90/Komega.F90 now have.
  UseScalarFallback = ( CurrentCoordinateSystem() /= Cartesian .AND. &
      CurrentCoordinateSystem() /= AxisSymmetric ) .OR. &
      ListCheckPresentAnyMaterial( Model, 'Compressibility Model' )

  IF (.NOT. ALLOCATED(SSTHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(SSTHandles(nthr))
  END IF

  GlobalBubbles = Solver % GlobalBubbles
  Stabilize = GetStabilizeFlag( Solver % Values, GotIt )

  ! Only LocalMatrixScalar's legacy per-node branch uses this -- see the
  ! matching BubblesDefault resolution in SSTKomegaLegacy.F90.
  BubblesDefault = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  IF ( .NOT.GotIt ) BubblesDefault = .TRUE.

  ! K and Omega are interleaved, Dofs=2 -- exactly as SSTKomegaLegacy.F90
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
    CALL Info(Caller,'SSTKomega iteration: '//I2S(iter), Level=4)
    CALL Info(Caller,'-------------------------------------', Level=4)

    CALL DefaultInitialize()

    Active = GetNOFActive()

    IF ( UseScalarFallback ) THEN
      ! Serial only: LocalMatrixScalar relies on the classic argument-less
      ! Get*() "current element" accessors, exactly as
      ! SSTKomegaLegacy.F90's own driver does.
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

      IF ( ListGetLogical( BC, 'Wall Law',gotIt ) ) THEN
        CALL KomegaWallLaw(Element,n,BC,Solver)
      END IF
    END DO

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
END SUBROUTINE SSTKOmega
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Element-type resolution, called before the mesh/basis functions are
!> finalized -- same role and timing as HeatSolver_Init0 and
!> Spalart-Allmaras.F90's/Komega.F90's own _Init0.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SSTKOmega_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SSTKOmegaFront
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
END SUBROUTINE SSTKOmega_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: SSTKOmega. Under "Legacy Assembly",
!> delegates to SSTKOmegaLegacy_Init and returns -- the p-bubble setup below
!> is specific to this file's own implementation.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SSTKOmega_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE SSTKOmegaFront
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
    CALL DelegateToSSTKOmegaLegacy( 'SSTKomegaLegacy_Init', Model, Solver, dt, TransientSimulation )
    RETURN
  END IF

  SolverParams => GetSolverParams()

  IF ( ListGetLogical( SolverParams, 'Bubbles', Found ) ) THEN
    CALL Warn('SSTKOmega_Init', &
        '"Bubbles = True" (the legacy per-node scheme) is not used here -- '// &
        'SSTKOmega_Init0 already defaulted "Element" to an equivalent p-element '// &
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
END SUBROUTINE SSTKOmega_Init
!------------------------------------------------------------------------------
