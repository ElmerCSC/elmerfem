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
! *  Module for the solution of incompressible Navier-Stokes equation.
! *  Utilizes multithreading and vectorization features initially introduced by Mikko Byckling.
! *  Replaces partly the legacy solver FlowSolve which is not optimized.
! *
! *  Authors: Mika Malinen, Juhani Kataja, Juha Ruokolainen, Peter Råback, Thomas Zwinger
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 28.01.2019
! *
!/*****************************************************************************/

MODULE IncompressibleLocalForms

  USE DefUtils

  REAL(KIND=dp), ALLOCATABLE, SAVE :: bx(:), bxprev(:)

  ! Per-thread handle/cache storage for LocalBulkMatrix and EffectiveViscosityVec.
  ! NOT THREADPRIVATE (Windows/GCC emutls bug inherits master's ALLOCATABLE/POINTER
  ! THREADPRIVATE data into workers instead of giving independent copies). Instead,
  ! this is a plain module-level array indexed by omp_get_thread_num()+1; each thread
  ! only ever touches its own slot in the parallel bulk-assembly loop. Slot 1 is filled
  ! by the serial InitHandles=.TRUE. pass (element 1); IncompressibleNSSolver then
  ! broadcasts slot 1 to all other slots before the parallel region starts, so every
  ! thread enters with identical resolved handle configuration.
  TYPE :: NSHandles_t
    TYPE(ValueHandle_t) :: Dens_h, Load_h(3)
    TYPE(ValueHandle_t) :: Visc_h, ViscModel_h, ViscExp_h, ViscCritical_h, &
        ViscNominal_h, ViscDiff_h, ViscTrans_h, ViscYasuda_h, ViscGlenExp_h, ViscGlenFactor_h, &
        ViscArrSet_h, ViscArr_h, ViscTLimit_h, ViscRate1_h, ViscRate2_h, ViscEne1_h, ViscEne2_h, &
        ViscTemp_h
    REAL(KIND=dp) :: R = 8.314_dp, NewtonRelax = 0.0_dp
    LOGICAL :: ConstantVisc = .FALSE., Visited = .FALSE., GotRelax = .FALSE.
    LOGICAL :: SaveShear = .FALSE., SaveVisc = .FALSE., SaveWeight = .FALSE.
    TYPE(Variable_t), POINTER :: ShearVar => NULL(), ViscVar => NULL(), WeightVar => NULL()
  END TYPE NSHandles_t
  TYPE(NSHandles_t), ALLOCATABLE :: NSHandles(:)

CONTAINS

!------------------------------------------------------------------------------
! Assemble local finite element matrix for a single bulk element and glue
! it to the global matrix. Routine is vectorized and multithreaded.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(Element, n, nd, ntot, dim, &
       DivCurlForm, GradPVersion, SpecificLoad, StokesFlow, &
       dt, LinearAssembly, nb, Newton, Transient, InitHandles, SchurSolver, &
       PStab, PStabCoeff, PStabHmode )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, ntot, dim, nb
    LOGICAL, INTENT(IN) :: DivCurlForm, GradPVersion, SpecificLoad, StokesFlow
    REAL(KIND=dp), INTENT(IN) :: dt   
    LOGICAL, INTENT(IN) :: LinearAssembly, Newton, Transient, InitHandles
    TYPE(Solver_t), POINTER :: SchurSolver
    LOGICAL, INTENT(IN) :: PStab
    REAL(KIND=dp), INTENT(IN) :: PStabCoeff
    INTEGER, INTENT(IN) :: PStabHmode
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    LOGICAL :: Stat, Found

    REAL(KIND=dp), TARGET :: MASS(ntot*(dim+1),ntot*(dim+1)), &
        STIFF(ntot*(dim+1),ntot*(dim+1)), FORCE(ntot*(dim+1))
    REAL(KIND=dp) :: NodalSol(dim+1,ntot)
    REAL(KIND=dp) :: PrevNodalSol(dim+1,ntot)
    REAL(KIND=dp) :: s, rho

    REAL(KIND=dp), ALLOCATABLE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:), &
        rhoVec(:), VeloPresVec(:,:), loadAtIpVec(:,:), VelocityMass(:,:), &
        PressureMass(:,:), ForcePart(:), &
        weight_a(:), weight_b(:), weight_c(:), tauVec(:), PrevTempVec(:), PrevPressureVec(:), &
        VeloVec(:,:), PresVec(:), GradVec(:,:,:), ConvVec(:,:)
    REAL(KIND=dp), POINTER :: muVec(:), LoadVec(:)
    REAL(KIND=dp), ALLOCATABLE :: muDerVec0(:),g(:,:,:),StrainRateVec(:,:,:)
    ! Caller-owned scratch for the viscosity vector; see EffectiveViscosityVec.
    REAL(KIND=dp), ALLOCATABLE, TARGET :: ViscWork(:)
    REAL(kind=dp) :: stifford(ntot,ntot,dim+1,dim+1), jacord(ntot,ntot,dim+1,dim+1), &
        JAC(ntot*(dim+1),ntot*(dim+1) )

    INTEGER :: t, i, j, k, p, q, ngp, allocstat, elemdim
    INTEGER :: DOFs, tid
    REAL(KIND=dp) :: Kfam, hStab, VNorm, PeStab

!DIR$ ATTRIBUTES ALIGN:64 :: BasisVec, dBasisdxVec, DetJVec, rhoVec, VeloPresVec, loadAtIpVec
!DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE, weight_a, weight_b, weight_c
!------------------------------------------------------------------------------

    ! Dens_h/Load_h and (inside EffectiveViscosityVec) the viscosity handles live in
    ! NSHandles(tid) — a per-thread slot, not SAVE/THREADPRIVATE. See NSHandles_t above.
    tid = 1
    !$ tid = OMP_GET_THREAD_NUM() + 1

    ASSOCIATE( Dens_h => NSHandles(tid) % Dens_h, Load_h => NSHandles(tid) % Load_h )

    CALL GetElementNodesVec( Nodes )
    STIFF = 0._dp
    MASS  = 0._dp
    FORCE = 0._dp
    JAC   = 0._dp
    JacOrd = 0._dp
    stifford = 0._dp

    DOFs = dim + 1

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt(Element)
    ngp = IP % n

    ElemDim = Element % Type % Dimension
    
    ! Storage size depending ngp
    !-------------------------------------------------------------------------------

    ALLOCATE(BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
        rhoVec(ngp), VeloVec(ngp, dim), PresVec(ngp), velopresvec(ngp,dofs), LoadAtIpVec(ngp,dim+1), &
        weight_a(ngp), weight_b(ngp), weight_c(ngp), tauVec(ngp), PrevTempVec(ngp), &
        PrevPressureVec(ngp), GradVec(ngp,dim,dim), ConvVec(ngp,ntot), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('IncompressibleNSSolver::LocalBulkMatrix','Local storage allocation failed')

    ALLOCATE(VelocityMass(ntot,ntot), PressureMass(ntot, ntot), ForcePart(ntot))
           
    IF (Newton) THEN
      ALLOCATE(muDerVec0(ngp), g(ngp,ntot,dim), StrainRateVec(ngp,dim,dim))
      muDerVec0 = 0._dp
    END IF

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Dens_h,'Material','Density')      
      CALL ListInitElementKeyword( Load_h(1),'Body Force','Flow Bodyforce 1')      
      CALL ListInitElementKeyword( Load_h(2),'Body Force','Flow Bodyforce 2')      
      CALL ListInitElementKeyword( Load_h(3),'Body Force','Flow Bodyforce 3')      
    END IF

    ! We assume constant density so far:
    !-----------------------------------
    rho = ListGetElementReal( Dens_h, Element = Element ) 

    ! Get the previous elementwise velocity-pressure iterate:
    !--------------------------------------------------------
    IF ( LinearAssembly ) THEN
      VeloPresVec = 0._dp
    ELSE
      CALL GetLocalSolution( NodalSol )
      IF (nb > 0 .AND. Transient .AND. .NOT. StokesFlow) & 
         CALL GetLocalSolution(PrevNodalSol, tStep=-1)
    END IF

    VelocityMass = 0.0d0

    ! Vectorized basis functions
    stat = ElementInfoVec(Element, nodes, ngp, IP % U, IP % V, &
        IP % W, detJvec, SIZE(basisVec, 2), BasisVec, dBasisdxVec)

    ! Weights at integration points
    DO t = 1, ngp
      DetJVec(t) = DetJVec(t) * IP % s(t)
    END DO

    ! Velocity and pressure from previous iteration at integration points
    IF(.NOT. LinearAssembly) THEN
      CALL DGEMM('N', 'T', ngp, DOFs, n, 1._dp, BasisVec, ngp, NodalSol, &
          dofs, 0._dp, VeloPresVec, ngp)
    END IF

    ! Return the effective viscosity. Currently only non-newtonian models supported.
    muvec => EffectiveViscosityVec( ngp, BasisVec, dBasisdxVec, Element, NodalSol, &
              muDerVec0, Newton,  InitHandles, DetJVec, ViscWork )        

    ! Rho 
    rhovec(1:ngp) = rho

    ! Flow bodyforce if present
    LoadAtIpVec = 0._dp
    DO i=1,dim
      LoadVec => ListGetElementRealVec( Load_h(i), ngp, BasisVec, Element, Found ) 
      IF( Found ) THEN
        IF (SpecificLoad) THEN
          LoadAtIpVec(1:ngp,i) = LoadVec(1:ngp)
        ELSE
          LoadAtIpVec(1:ngp,i) = rho * LoadVec(1:ngp)
        END IF
      END IF
    END DO

    IF ( Newton ) THEN

      DO i = 1,dim
        DO j = 1,dim
          GradVec(1:ngp, i, j) = MATMUL(dBasisdxVec(1:ngp,1:ntot,j),nodalsol(i,1:ntot))
        END DO
      END DO

      IF( .NOT. StokesFlow ) THEN
        DO i = 1,dim
          LoadAtIpVec(1:ngp, i) = LoadAtIpVec(1:ngp, i) + rhovec(1:ngp) * &
               SUM(gradvec(1:ngp,i,1:dim)*velopresvec(1:ngp,1:dim),2)
        END DO
      END IF

      IF (ANY(muDerVec0(1:ngp)/=0)) THEN
        DO i = 1,dim
          DO j = 1,dim
            StrainRateVec(1:ngp,i,j) = ( GradVec(1:ngp,i,j) + GradVec(1:ngp,j,i) ) / 2
          END DO
        END DO

        muDerVec0(1:ngp) = muderVec0(1:ngp)*detJVec(1:ngp)*8
        DO i=1,dim
          DO q = 1,ntot
            g(1:ngp,q,i) = SUM(StrainRateVec(1:ngp,i,:)*dBasisdxvec(1:ngp,q,1:dim),2)
          END DO
        END DO

        DO i=1,dim
          DO j=1,dim
            CALL LinearForms_udotv(ngp,ntot,dim,g(:,:,j),g(:,:,i),mudervec0,jacord(:,:,j,i))
          END DO
        END DO
      END IF
    END IF

    weight_a(1:ngp) = muVec(1:ngp) * detJVec(1:ngp)
    IF (DivCurlForm) THEN
      weight_b(1:ngp) = -weight_a(1:ngp)
      ! The following assumes that the bulk viscosity of the fluid vanishes:
      weight_c(1:ngp) =  4.0_dp / 3.0_dp * weight_a(1:ngp)

      DO i = 1, dim
        DO j = 1, dim
          ! Div-Div term
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dbasisdxvec(:,:,i), dbasisdxvec(:,:,j), weight_c, stifford(:,:,i,j))

          IF (i == j) THEN
            ! Diagonal terms for curl-curl 
            DO k = 1, dim
              IF (k /= i) THEN
                CALL LinearForms_UdotV(ngp, ntot, elemdim, &
                    dbasisdxvec(:,:,k), dbasisdxvec(:,:,k), weight_a, stifford(:,:,i,j))
              END IF
            END DO
          ELSE
            ! Off-diagonal terms: Cross derivatives for curl-curl 
            CALL LinearForms_UdotV(ngp, ntot, elemdim, &
                dbasisdxvec(:,:,j), dbasisdxvec(:,:,i), weight_b, stifford(:,:,i,j))
          END IF
        END DO
      END DO

    ELSE ! Standard Grad-U Form (Non-DivCurlForm)
      DO i=1,dim
        DO j=1,dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(:,:,j), dBasisdxVec(:,:,j), weight_a, stifford(:,:,i,i))

          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(:,:,j), dBasisdxVec(:,:,i), weight_a, stifford(:,:,i,j))
        END DO
      END DO
    END IF

    IF (GradPVersion) THEN
       ! b(u,q) = (u, grad q) part
      DO i = 1, dim
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            BasisVec, dbasisdxvec(:,:,i), detJVec, stifford(:,:,i,dofs))
        StiffOrd(:,:,dofs,i) = transpose(stifford(:,:,i,dofs))
      END DO
    ELSE
       DO i = 1, dim
         CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dBasisdxVec(:, :, i), BasisVec, -detJVec, StiffOrd(:,:,i,dofs))
        StiffOrd(:,:,dofs,i) = transpose(stifford(:,:,i,dofs))
      END DO
    END IF

    !---------------------------------------------------------------------------
    ! PRESSURE STABILISATION, the equal-order alternative to the bubble.
    !
    ! The pressure/pressure block is empty in every inf-sup stable pair above --
    ! that is what a saddle point matrix looks like -- and here it is filled with
    !
    !   -tau * (grad q, grad p)
    !
    ! the Brezzi-Pitkaranta term, plus the load part of the momentum residual on
    ! the right. The sign is negative because the constraint row assembled just
    ! above is -(div u, q), so the system stays [A B^T; B -C] with C positive
    ! semi-definite rather than turning the saddle point into something indefinite
    ! in the wrong direction.
    !
    ! tau is built AT THE INTEGRATION POINTS, not once per element, because muVec
    ! is a shear-rate dependent viscosity here: a non-Newtonian element can span a
    ! decade in mu, and a single element value would over-stabilise the stiff end
    ! of it and under-stabilise the soft end.
    !
    ! The load term is kept, which is where this departs from the same scheme in
    ! ElasticSolve. There the load is a surface traction and dropping the volume
    ! residual is the classical Brezzi-Pitkaranta simplification, consistent to
    ! O(h^2). Here the load is typically gravity, and the whole pressure field is
    ! the response to it -- grad p balances rho*g to leading order -- so dropping
    ! it would leave the stabilisation fighting the hydrostatic gradient itself.
    ! With it, the term vanishes identically on the exact hydrostatic solution.
    !
    ! What is deliberately NOT in the residual is -mu*div(grad u), which for the
    ! equal-order linear basis this path forces is zero on simplices anyway, and
    ! the convective and transient terms -- see the Stokes-only guard in the
    ! solver. Adding those is what turns this into a full PSPG/SUPG scheme.
    !---------------------------------------------------------------------------
    IF ( PStab ) THEN
      ! Per element family. The starting point is the same scheme in ElasticSolve,
      ! whose incompressible limit IS Stokes with the shear modulus in place of
      ! the viscosity, so the constants transfer in principle. Simplices need far
      ! less than tensor elements of the same size, consistent with hK being a
      ! diameter and so overstating a simplex.
      !
      ! MEASURED HERE, on the backward facing step of Step_stokes_vec, as the
      ! value at which the stabilised pair reproduces the default MINI answer.
      ! The optimum is a real minimum -- the error changes sign through it -- and
      ! it is mesh independent, which is what says the h^2/mu scaling is right:
      !
      !   quad      500 / 2000 / 8000 elements, optimum at 0.63x the elastic
      !             constant on all three; at 1x the error is -0.434 / -0.128 /
      !             -0.033 %, so it also SHRINKS under refinement rather than
      !             converging to something merely close
      !   triangle  1062 / 4012 elements, optimum at 0.92x, i.e. the elastic
      !             constant already
      !
      ! So the quadrilateral is restated at 0.0417*0.63 and the triangle at
      ! 0.0104*0.92. That the two families moved by different factors is the
      ! point: the discrepancy is NOT a global difference between the elastic and
      ! flow references that could be absorbed into the coefficient. It is the
      ! problem sensitivity ElasticSolve's own note flags as untested, now
      ! measured once -- tens of percent on a quadrilateral, nothing on a
      ! triangle.
      !
      ! The hexahedron is measured too, on FlowResNoslip -- flow around a hole,
      ! a native 3D brick mesh -- where the optimum is 1.1x, so the inherited
      ! value stands as it is.
      !
      ! A glacier slab, Friction_WeertmanNewton3D, appears to want 4.6x instead --
      ! but that is a RESOLUTION artefact, not a property of the problem, and the
      ! way it was first measured is the trap. Refining that case by hand gave
      ! 4.8, 4.2 and 4.9 at 10x10x3, 20x20x3 and 40x40x5, which reads as mesh
      ! independent. It is not: those meshes refine the HORIZONTAL divisions while
      ! the vertical stays at three levels and then five. The slab is thin and its
      ! velocity gradient is vertical, so the vertical spacing is what limits it,
      ! and that was very nearly held fixed across all three.
      !
      ! "Mesh Levels", which splits the same mesh uniformly, refines the direction
      ! that matters and the crossing walks straight down:
      !
      !   Mesh Levels    1      2      3
      !   optimum       4.2    2.0    1.6
      !
      ! heading for the same order as every other family. The error at coefficient
      ! one falls with it, 0.32 -> 0.12 -> 0.086 %, so a refined glacier case needs
      ! no special coefficient at all.
      !
      ! The lesson for anyone measuring one of these constants again: refine with
      ! "Mesh Levels" rather than by regenerating meshes. Regenerating changes
      ! element shape and resolution together, and a constant fitted that way can
      ! look settled across three meshes while tracking whichever dimension was
      ! not being refined.
      !
      ! Much of that factor is the LENGTH SCALE rather than the physics, and it is
      ! worth knowing that hK is not one measure. ElementDiameter returns the
      ! SHORTEST EDGE for every family in its default branch -- hexahedron,
      ! tetrahedron, prism, pyramid -- an area based quantity for the triangle,
      ! and a harmonic mean of two edges for the quadrilateral. So part of the
      ! spread between the per-family constants above is the measure changing
      ! under them, not the elements differing.
      !
      ! "Pressure Stabilization Element Size" exists because of what happens when
      ! the measure is varied, at coefficient one:
      !
      !                    long/short   shortest    longest     volume
      !   glacier slab        2.0        +0.323 %   +0.015 %   +0.147 %
      !   flow past hole      2.95       +0.038 %   -0.828 %   -0.326 %
      !
      ! Each problem is best on a DIFFERENT measure, and the more anisotropic of
      ! the two is the one happy with the shortest edge -- so this is not a
      ! geometric correction waiting to be derived. Taking the longest edge on
      ! the slab reproduces MINI outright, which is where its 4.6x went: long
      ! over short is 2 there and tau goes as h^2. The same substitution on the
      ! hole would predict 8.7x and the truth is 1.1x.
      !
      ! Hence a keyword and not a cleverer default. The default stays hK, which
      ! is what every constant above was measured against.
      !
      ! THE TETRAHEDRON AND PYRAMID REMAIN UNMEASURED for this solver -- there is
      ! no incompressible case of either to take a measurement from, exactly as
      ! in ElasticSolve. Note also that extruding a 2D case is not a shortcut to
      ! measuring anything: see the hK guard in the solver.
      SELECT CASE( Element % TYPE % ElementCode / 100 )
      CASE( 3 )
        Kfam = 0.0096_dp
      CASE( 4 )
        Kfam = 0.0263_dp
      CASE( 5 )
        Kfam = 0.0029_dp
      CASE( 6 )
        Kfam = 0.005_dp
      CASE( 7 )
        Kfam = 0.0090_dp
      CASE DEFAULT
        Kfam = 0.0232_dp
      END SELECT

      ! The length scale. hK is the default and is what the constants above were
      ! measured with, but it is not one measure -- see the note above -- so the
      ! keyword lets a case pick a consistent one instead. The volume form is the
      ! only genuinely family independent choice of the three.
      SELECT CASE( PStabHmode )
      CASE( 1 )
        hStab = ElementDiameter( Element, Nodes, UseLongEdge = .TRUE. )
      CASE( 2 )
        hStab = SUM( DetJVec(1:ngp) ) ** ( 1.0_dp / elemdim )
      CASE DEFAULT
        hStab = Element % hK
      END SELECT

      ! An element size of zero makes tau vanish, and that does not fail -- it
      ! quietly returns the unstabilised equal-order system, whose pressure has a
      ! null space. The answer is then wherever round-off points, and it does not
      ! move when the coefficient is changed, which is the signature to know it
      ! by. So it is checked here, against the value actually used, rather than
      ! once at start-up against element one.
      !
      ! That distinction is the whole point. MeshStabParams fills hK when a mesh
      ! is loaded, but meshes get BUILT later too, and the low level builders do
      ! not call it: CreateExtrudedMesh does not, nor do RemeshMMG3D,
      ! SequentialRemeshParMMG or ParallelRemesh. Their callers all do today --
      ! the calving and adaptive paths call MeshStabParams on the way out -- so
      ! nothing is broken, but a mesh replaced MID RUN by adaptation would sail
      ! past a start-up test, and a new caller that forgot would too. It does not
      ! follow that every extruded or remeshed case is affected either: anything
      ! that later moves the mesh repairs hK on the way past, DisplaceMesh
      ! recomputing it per element, which is why FixTangentVelo never sees this.
      IF( hStab <= 0.0_dp ) CALL Fatal('IncompressibleNSSolver::LocalBulkMatrix', &
          'Stabilization element size is zero, so the stabilisation would '// &
          'silently vanish. This mesh has not been through MeshStabParams -- a '// &
          'freshly extruded or remeshed one is how that happens')

      ! tau. The diffusive expression Kfam*h^2/mu is what the constants were
      ! measured against and is exact for Stokes. With advection it becomes the
      ! interpolating form NavierStokes uses, scaled so that the weak advection
      ! limit returns that same expression rather than merely resembling it --
      ! the identical construction as in ElasticSolve, so a case moving between
      ! the two solvers sees one scheme.
      !
      ! Per integration point, since both |u| and mu vary within an element.
      IF( StokesFlow ) THEN
        tauVec(1:ngp) = PStabCoeff * Kfam * hStab**2 / muVec(1:ngp)
      ELSE
        DO t = 1, ngp
          VNorm = SQRT( SUM( VeloPresVec(t,1:dim)**2 ) )
          IF( rhoVec(t) * VNorm * hStab > 0.0_dp ) THEN
            PeStab = MIN( 1.0_dp, 2.0_dp * Kfam * rhoVec(t) * VNorm * hStab / muVec(t) )
            tauVec(t) = hStab * PeStab / ( 2.0_dp * rhoVec(t) * VNorm )
          ELSE
            tauVec(t) = Kfam * hStab**2 / muVec(t)
          END IF
        END DO
        tauVec(1:ngp) = PStabCoeff * tauVec(1:ngp)
      END IF

      ! The convective operator applied to each basis function, rho*(u.grad)phi_q.
      ! It is both the convective part of the residual SU and, tested against
      ! itself, the SUPG weight SW -- so it is built once and used four ways.
      ConvVec = 0.0_dp
      IF( .NOT. StokesFlow ) THEN
        DO k = 1, dim
          DO q = 1, ntot
            ConvVec(1:ngp,q) = ConvVec(1:ngp,q) + &
                rhoVec(1:ngp) * VeloPresVec(1:ngp,k) * dBasisdxVec(1:ngp,q,k)
          END DO
        END DO
      END IF

      ! In LinearForms_UdotV the FIRST basis argument indexes the row (the test
      ! function) and the second the column, which matters for every term below
      ! but not for the symmetric one this started as.
      !
      ! -tau * (grad q, grad p) into the pressure/pressure block
      weight_a(1:ngp) = -tauVec(1:ngp) * detJVec(1:ngp)
      DO i = 1, dim
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dBasisdxVec(:,:,i), dBasisdxVec(:,:,i), weight_a, StiffOrd(:,:,dofs,dofs))
      END DO

      IF( .NOT. StokesFlow ) THEN
        ! The rest of the residual, and the SUPG weight. Writing R for
        ! rho*(u.grad)u + grad p - f, the scheme is
        !
        !   momentum_i  +=  tau * ( rho*(u.grad)v_i , R_i )      SUPG
        !   constraint  += -tau * ( grad q          , R   )      PSPG
        !
        ! and the four matrix blocks below are those two weights against the two
        ! operators in R. The load half of each goes in with the other forces.
        !
        ! The viscous part of R, -div(mu grad u), is NOT included. On the linear
        ! equal-order basis this path forces it vanishes on simplices, and the
        ! vectorized ElementInfoVec returns no second derivatives to build it from
        ! on the others -- which is the classical simplification, consistent to
        ! the order of the discretisation, not an oversight.

        ! PSPG against convection: -tau * (grad q, rho*(u.grad)u)
        DO j = 1, dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(:,:,j), ConvVec, weight_a, StiffOrd(:,:,dofs,j))
        END DO

        weight_b(1:ngp) = tauVec(1:ngp) * detJVec(1:ngp)

        ! SUPG against convection, diagonal in the velocity components
        DO i = 1, dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              ConvVec, ConvVec, weight_b, StiffOrd(:,:,i,i))
        END DO

        ! SUPG against the pressure gradient
        DO i = 1, dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              ConvVec, dBasisdxVec(:,:,i), weight_b, StiffOrd(:,:,i,dofs))
        END DO

        ! The Newton term of the residual, rho*(v.grad)u linearised, so that the
        ! Jacobian matches the operator the stabilisation actually applies.
        IF( Newton ) THEN
          DO i = 1, dim
            DO j = 1, dim
              CALL LinearForms_UdotV(ngp, ntot, elemdim, &
                  ConvVec, BasisVec, weight_b, StiffOrd(:,:,i,j), &
                  rhoVec*GradVec(:,i,j))
              CALL LinearForms_UdotV(ngp, ntot, elemdim, &
                  dBasisdxVec(:,:,i), BasisVec, weight_a, StiffOrd(:,:,dofs,j), &
                  rhoVec*GradVec(:,i,j))
            END DO
          END DO
        END IF
      END IF
    END IF

    ! Masses (use symmetry)
    ! Compute bilinear form G=G+(alpha u, u) = u .dot. (grad u) 
    IF ( .NOT. StokesFlow ) THEN
      CALL LinearForms_UdotU(ngp, ntot, elemdim, BasisVec, DetJVec, VelocityMass, rhovec)

      ! Scatter to the usual local mass matrix
      DO i = 1, dim
        mass(i::dofs, i::dofs) = mass(i::dofs, i::dofs) + VelocityMass(1:ntot, 1:ntot)
      END DO

      ! Accumulate convection block once, then scatter to diagonal velocity blocks
      VelocityMass = 0._dp
      DO k = 1, dim
        weight_a(1:ngp) = rhovec(1:ngp) * veloPresVec(1:ngp,k)
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            basisvec, dbasisdxvec(:,:,k), detJvec, VelocityMass, weight_a)
      END DO
      DO i = 1, dim
        stifford(1:ntot,1:ntot,i,i) = stifford(1:ntot,1:ntot,i,i) + VelocityMass(1:ntot,1:ntot)
      END DO

      IF ( Newton ) THEN
        DO i = 1, dim
          DO j = 1, dim
            CALL LinearForms_UdotV(ngp, ntot, elemdim, &
                basisvec, basisvec, detJvec, stifford(:,:,i,j), rhovec*gradvec(:,i,j))
          END DO
        END DO
      END IF
    END IF

    ! add loads
    DO i = 1,dim+1
      ForcePart = 0._dp
      CALL LinearForms_UdotF(ngp, ntot, basisVec, detJVec, LoadAtIpVec(:,i), ForcePart)
      FORCE(i::dofs) = ForcePart(1:ntot)
    END DO

    ! The load half of the stabilised constraint equation, -tau*(grad q, f). The
    ! matrix half sits in the pressure/pressure block far above; this is the same
    ! term with grad p replaced by the load it balances, so the two cancel on a
    ! hydrostatic field and the stabilisation costs nothing there.
    !
    ! weight_a is rebuilt rather than carried down from that block: the convection
    ! assembly between here and there overwrites it. tauVec is untouched.
    IF ( PStab ) THEN
      weight_a(1:ngp) = -tauVec(1:ngp) * detJVec(1:ngp)
      ForcePart = 0._dp
      DO i = 1,dim
        CALL LinearForms_UdotF(ngp, ntot, dBasisdxVec(:,:,i), weight_a, &
            LoadAtIpVec(:,i), ForcePart)
      END DO
      FORCE(dofs::dofs) = FORCE(dofs::dofs) + ForcePart(1:ntot)

      ! The SUPG half of the same load, this one per velocity component and with
      ! the opposite sign, the momentum rows carrying the load on the right.
      IF( .NOT. StokesFlow ) THEN
        weight_b(1:ngp) = tauVec(1:ngp) * detJVec(1:ngp)
        DO i = 1,dim
          ForcePart = 0._dp
          CALL LinearForms_UdotF(ngp, ntot, ConvVec, weight_b, &
              LoadAtIpVec(:,i), ForcePart)
          FORCE(i::dofs) = FORCE(i::dofs) + ForcePart(1:ntot)
        END DO
      END IF
    END IF

    DO i = 1, DOFS
      DO j = 1, DOFS
        Stiff(i::DOFS, j::DOFS) = StiffOrd(1:ntot, 1:ntot, i,j)
      END DO
    END DO

    IF ( Newton ) THEN
      BLOCK
        REAL(KIND=dp) :: SOL(ntot*(dim+1))

        SOL=0._dp
        DO i = 1, DOFS
          DO j = 1, DOFS
            JAC(i::DOFS, j::DOFS) = JacOrd(1:ntot, 1:ntot, i,j)
          END DO
          SOL(i::DOFs) = NodalSol(i,1:ntot)
        END DO

        STIFF = STIFF + JAC
        FORCE = FORCE + MATMUL(JAC,SOL)
      END BLOCK
    END IF

    IF(StokesFlow) THEN
      IF ( nb>0 ) THEN
        CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE)
      ELSE
        DO p = n+1,ntot
          i = DOFs * p
          FORCE(i)   = 0._dp
          STIFF(i,:) = 0._dp
          STIFF(:,i) = 0._dp
          STIFF(i,i) = 1._dp
        END DO
      END IF

    ELSE IF (nb > 0 .AND. nd==n .AND. Transient) THEN
      !-------------------------------------------------------------------------
      ! This branch is primarily intended to handle the (enhanced) MINI element 
      ! approximation together with the static condensation for the velocity
      ! bubbles. The subroutine LCondensate constructs the time derivative of 
      ! the bubble-augmented part and performs the static condensation.
      !-------------------------------------------------------------------------
      CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE, PrevNodalSol, &
                   NodalSol, Element % ElementIndex)
    ELSE
      !-------------------------------------------------------------------------
      ! The cases handled here include the MINI element approximation with the 
      ! velocity bubbles left in the global system and P2/Q2-P1/Q1 approximation.
      ! First, enforce P1/Q1 pressure approximation by setting Dirichlet 
      ! constraints for unused dofs: 
      !-------------------------------------------------------------------------
      DO p = n+1,ntot
        i = DOFs * p
        FORCE(i)   = 0._dp
        MASS(:,i)  = 0._dp
        MASS(i,:)  = 0._dp
        STIFF(i,:) = 0._dp
        STIFF(:,i) = 0._dp
        STIFF(i,i) = 1._dp
      END DO

      IF ( Transient ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF
      IF (nb > 0) THEN
        IF (Transient) THEN
          CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE, &
              PrevNodalSol, NodalSol, Element % ElementIndex)
        ELSE
          CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE)
        END IF
      END IF
    END IF

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element, VecAssembly=.TRUE.)

    IF( ASSOCIATED( SchurSolver ) ) THEN
      ! Preconditioner for pressure block when using block preconditioning               
      weight_a(1:ngp) = -1.0_dp / muvec(1:ngp) * detJVec(1:ngp)
      PressureMass = 0.0_dp
      FORCE = 0.0_dp
      CALL LinearForms_UdotU(ngp, nd, elemdim, BasisVec, weight_a, PressureMass)
      CALL DefaultUpdateEquations( PressureMass, FORCE, UElement=Element, &
                 Usolver = SchurSolver, VecAssembly = .TRUE.)
    END IF

    END ASSOCIATE

!------------------------------------------------------------------------------

  CONTAINS


    FUNCTION EffectiveViscosityVec( ngp, BasisVec, dBasisdxVec, Element, NodalSol, &
        ViscDerVec, ViscNewton, InitHandles, DetJVec, ViscWork ) RESULT ( EffViscVec ) 

      INTEGER :: ngp
      REAL(KIND=dp) :: BasisVec(:,:), dBasisdxVec(:,:,:)
      TYPE(Element_t), POINTER :: Element
      REAL(KIND=dp) :: NodalSol(:,:)
      REAL(KIND=dp), ALLOCATABLE :: ViscDerVec(:)
      LOGICAL :: InitHandles , ViscNewton
      REAL(KIND=dp), POINTER  :: EffViscVec(:)
      REAL(KIND=dp), ALLOCATABLE :: DetJVec(:)
      ! ViscWork is supplied by the caller so that the vector this function returns a
      ! pointer to is owned by the caller's frame: a local of a routine called inside
      ! the parallel region, hence per-thread by construction and released on return.
      ! It was POINTER+SAVE+THREADPRIVATE before; as a plain local pointer it leaked
      ! ngp reals per call, and routing it through a module-level per-thread struct
      ! ended up shared between threads.
      REAL(KIND=dp), ALLOCATABLE, TARGET :: ViscWork(:)

      LOGICAL :: Found
      CHARACTER(LEN=MAX_NAME_LEN) :: ViscModel
      REAL(KIND=dp) :: c1, c2, c3, c4, Ehf, Tlimit, ArrheniusFactor, A1, A2, Q1, Q2, ViscCond
      REAL(KIND=dp), ALLOCATABLE :: ss(:), s(:), ArrheniusFactorVec(:)
      REAL(KIND=dp), POINTER :: ViscVec0(:), ViscVec(:), TempVec(:), EhfVec(:)
      CHARACTER(*), PARAMETER :: Caller = 'EffectiveViscosityVec'

      ! All handle/cache state lives in NSHandles(tid) (see NSHandles_t, module scope)
      ! — not SAVE/THREADPRIVATE. "tid" is host-associated from LocalBulkMatrix.
      ASSOCIATE( &
          Visc_h => NSHandles(tid) % Visc_h, ViscModel_h => NSHandles(tid) % ViscModel_h, &
          ViscExp_h => NSHandles(tid) % ViscExp_h, ViscCritical_h => NSHandles(tid) % ViscCritical_h, &
          ViscNominal_h => NSHandles(tid) % ViscNominal_h, ViscDiff_h => NSHandles(tid) % ViscDiff_h, &
          ViscTrans_h => NSHandles(tid) % ViscTrans_h, ViscYasuda_h => NSHandles(tid) % ViscYasuda_h, &
          ViscGlenExp_h => NSHandles(tid) % ViscGlenExp_h, ViscGlenFactor_h => NSHandles(tid) % ViscGlenFactor_h, &
          ViscArrSet_h => NSHandles(tid) % ViscArrSet_h, ViscArr_h => NSHandles(tid) % ViscArr_h, &
          ViscTLimit_h => NSHandles(tid) % ViscTLimit_h, ViscRate1_h => NSHandles(tid) % ViscRate1_h, &
          ViscRate2_h => NSHandles(tid) % ViscRate2_h, ViscEne1_h => NSHandles(tid) % ViscEne1_h, &
          ViscEne2_h => NSHandles(tid) % ViscEne2_h, ViscTemp_h => NSHandles(tid) % ViscTemp_h, &
          R => NSHandles(tid) % R, NewtonRelax => NSHandles(tid) % NewtonRelax, &
          ConstantVisc => NSHandles(tid) % ConstantVisc, Visited => NSHandles(tid) % Visited, &
          GotRelax => NSHandles(tid) % GotRelax, &
          SaveShear => NSHandles(tid) % SaveShear, SaveVisc => NSHandles(tid) % SaveVisc, &
          SaveWeight => NSHandles(tid) % SaveWeight )
      ! ShearVar/ViscVar/WeightVar are reassigned (=>) below, so they cannot be ASSOCIATE
      ! names (Fortran does not let an associate-name for a POINTER component be
      ! re-pointed) — referenced via NSHandles(tid)%... directly instead.

      IF(InitHandles ) THEN
        CALL Info(Caller,'Initializing handles for viscosity models',Level=8)

        CALL ListInitElementKeyword( Visc_h,'Material','Viscosity')      
        CALL ListInitElementKeyword( ViscModel_h,'Material','Viscosity Model')      

        IF( ListGetElementSomewhere( ViscModel_h) ) THEN
          ViscCond = ListGetCReal( CurrentModel % Solver % Values,&
              'Newtonian Viscosity Condition',Found )      
          ConstantVisc = ( Found .AND. ViscCond > 0.0_dp ) 
          
          IF( ListGetLogical( CurrentModel % Solver % Values,&
              'Constant-Viscosity Start', Found) ) ConstantVisc = (.NOT. Visited ) 
          
          CALL ListInitElementKeyword( ViscExp_h,'Material','Viscosity Exponent')      
          CALL ListInitElementKeyword( ViscCritical_h,'Material','Critical Shear Rate')      
          CALL ListInitElementKeyword( ViscNominal_h,'Material','Nominal Shear Rate')      
          CALL ListInitElementKeyword( ViscDiff_h,'Material','Viscosity Difference')      
          CALL ListInitElementKeyword( ViscTrans_h,'Material','Viscosity Transition')      
          CALL ListInitElementKeyword( ViscYasuda_h,'Material','Yasuda Exponent')      

          ! Do these initializations for glen's model only
          IF ( ListCompareElementAnyString( ViscModel_h,'glen') ) THEN
            CALL ListInitElementKeyword( ViscGlenExp_h,'Material','Glen Exponent',DefRValue=3.0_dp)
            CALL ListInitElementKeyword( ViscGlenFactor_h,'Material','Glen Enhancement Factor',DefRValue=1.0_dp)           
            CALL ListInitElementKeyword( ViscArrSet_h,'Material','Set Arrhenius Factor',DefLValue=.FALSE.)
            CALL ListInitElementKeyword( ViscArr_h,'Material','Arrhenius Factor')            
            CALL ListInitElementKeyword( ViscTLimit_h,'Material','Limit Temperature',DefRValue=-10.0_dp)
            CALL ListInitElementKeyword( ViscRate1_h,'Material','Rate Factor 1',DefRValue=3.985d-13)
            CALL ListInitElementKeyword( ViscRate2_h,'Material','Rate Factor 2',DefRValue=1.916d3)
            CALL ListInitElementKeyword( ViscEne1_h,'Material','Activation Energy 1',DefRValue=60.0d03)
            CALL ListInitElementKeyword( ViscEne2_h,'Material','Activation Energy 2',DefRValue=139.0d03)       
            CALL ListInitElementKeyword( ViscTemp_h,'Material','Relative Temperature')            

            IF (.NOT.ListCheckPresentAnyMaterial( CurrentModel,'Glen Allow Old Keywords')) THEN
              IF( ListCheckPresentAnyMaterial( CurrentModel,'Constant Temperature') ) THEN                
                CALL Warn(Caller,'Replace >Constant Temperature< with >Relative Temperature<')
              END IF
              IF( ListCheckPresentAnyMaterial( CurrentModel,'Temperature Field Variable') ) THEN             
                CALL Warn(Caller,'Replace >Temperature Field Variable< with >Relative Temperature = Equals ...<')
              END IF
            END IF
            IF (ViscArrSet_h % NotPresentAnywhere .AND. ViscTemp_h % NotPresentAnywhere ) THEN
              CALL Fatal(Caller,'>Relative Temperature< not given for viscosity model "glen"')
            END IF            
            
            IF( ListCheckPresentAnyMaterial( CurrentModel,'Glen Enhancement Factor Function')  ) THEN
              CALL Fatal(Caller,'No Glen function API yet!')
            END IF
            R = GetConstReal( CurrentModel % Constants,'Gas Constant',Found)
            IF (.NOT.Found) R = 8.314_dp

            NewtonRelax = ListGetCReal( CurrentModel % Solver % Values,&
                'Viscosity Newton Relaxation Factor',GotRelax )
          END IF
        END IF

        NSHandles(tid) % ShearVar => VariableGet( CurrentModel % Mesh % Variables,'Shearrate',ThisOnly=.TRUE.)
        SaveShear = ASSOCIATED(NSHandles(tid) % ShearVar)
        Found = .FALSE.
        IF(SaveShear) THEN
          IF(NSHandles(tid) % ShearVar % TYPE == Variable_on_gauss_points ) THEN
            CALL Info(Caller,'Saving "Shearrate" on ip points!',Level=10)
          ELSE IF( NSHandles(tid) % ShearVar % TYPE == Variable_on_elements ) THEN
            CALL Info(Caller,'Saving "Shearrate" on elements!',Level=10)
          ELSE IF( NSHandles(tid) % ShearVar % TYPE == Variable_on_nodes ) THEN
            CALL Info(Caller,'Saving "Shearrate" on nodes!',Level=10)
            Found = .TRUE.
          ELSE
            CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
          END IF
          NSHandles(tid) % ShearVar % Values = 0.0_dp
        END IF

        NSHandles(tid) % ViscVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
        SaveVisc = ASSOCIATED(NSHandles(tid) % ViscVar)        
        IF(SaveVisc) THEN
          IF(NSHandles(tid) % ViscVar % TYPE == Variable_on_gauss_points ) THEN
            CALL Info(Caller,'Saving "Viscosity" on ip points!',Level=10)
          ELSE IF( NSHandles(tid) % ViscVar % TYPE == Variable_on_elements ) THEN
            CALL Info(Caller,'Saving "Viscosity" on elements!',Level=10)
          ELSE IF( NSHandles(tid) % ViscVar % TYPE == Variable_on_nodes ) THEN
            CALL Info(Caller,'Saving "Viscosity" on nodes!',Level=10)
            Found = .TRUE.
          ELSE
            CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
          END IF
          NSHandles(tid) % ViscVar % Values = 0.0_dp
        END IF

        SaveWeight = .FALSE.
        IF( Found ) THEN
          NSHandles(tid) % WeightVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity Weight',ThisOnly=.TRUE.)
          IF( ASSOCIATED( NSHandles(tid) % WeightVar ) ) THEN
            IF( NSHandles(tid) % WeightVar % TYPE /= Variable_on_nodes ) THEN
              CALL Fatal(Caller,'Invalid field type for "Viscosity Weight"!')
            END IF
            SaveWeight = .TRUE.
            NSHandles(tid) % WeightVar % Values = 0.0_dp
          END IF
        END IF        
          
        Visited = .TRUE.
      END IF

      ViscVec0 => ListGetElementRealVec( Visc_h, ngp, BasisVec, Element )

      ViscModel = ListGetElementString( ViscModel_h, Element, Found ) 
      IF( .NOT. Found ) THEN
        ! Return the plain viscosity
        EffViscVec => ViscVec0
        RETURN
      END IF

      ! Initialize derivative of viscosity for when newtonian linearization is used
      IF( ViscNewton ) THEN
        ViscDerVec(1:ngp) = 0.0_dp
      END IF

      ! This reverts the viscosity model to linear 
      IF( ConstantVisc ) THEN
        EffViscVec => ViscVec0        
        RETURN      
      END IF
        
      ALLOCATE(ss(ngp), s(ngp), ArrheniusFactorVec(ngp))
      IF( .NOT. ALLOCATED( ViscWork ) ) THEN
        ALLOCATE( ViscWork(ngp) )
      ELSE IF( SIZE( ViscWork ) < ngp ) THEN
        DEALLOCATE( ViscWork )
        ALLOCATE( ViscWork(ngp) )
      END IF
      ViscVec => ViscWork(1:ngp)

      ! For non-newtonian models compute the viscosity here
      EffViscVec => ViscVec

      ! Calculate the strain rate velocity at all integration points
      ss(1:ngp) = 0.0_dp
      DO i=1,dim
        DO j=1,dim
          s(1:ngp) = MATMUL( dBasisdxVec(1:ngp,1:ntot,i), nodalsol(j,1:ntot) ) + &
              MATMUL( dBasisdxVec(1:ngp,1:ntot,j), nodalsol(i,1:ntot) )
          ss(1:ngp) = ss(1:ngp) + s(1:ngp)**2
        END DO
      END DO
      ss(1:ngp) = 0.5_dp * ss(1:ngp)

      IF(SaveShear) THEN
        IF( NSHandles(tid) % ShearVar % TYPE == Variable_on_nodes ) THEN
          DO i=1,n
            NSHandles(tid) % ShearVar % Values(NSHandles(tid) % ShearVar % Perm(Element % NodeIndexes(i))) = &
                NSHandles(tid) % ShearVar % Values(NSHandles(tid) % ShearVar % Perm(Element % NodeIndexes(i))) + &
                SUM(BasisVec(1:ngp,i)*ss(1:ngp)*DetJVec(1:ngp))
          END DO
        ELSE
          i = Element % ElementIndex
          IF( NSHandles(tid) % ShearVar % TYPE == Variable_on_gauss_points ) THEN
            j = NSHandles(tid) % ShearVar % Perm(i+1) - NSHandles(tid) % ShearVar % Perm(i)
            IF(j /= ngp) THEN
              CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Shearrate" got '//I2S(ngp))
            END IF
            NSHandles(tid) % ShearVar % Values( &
                NSHandles(tid) % ShearVar % Perm(i)+1 : NSHandles(tid) % ShearVar % Perm(i+1) ) = ss(1:ngp)
          ELSE IF( NSHandles(tid) % ShearVar % TYPE == Variable_on_elements ) THEN
            NSHandles(tid) % ShearVar % Values(NSHandles(tid) % ShearVar % Perm(i)) = SUM(ss(1:ngp)) / ngp
          END IF
        END IF
      END IF
            
      
      SELECT CASE( ViscModel )       

      CASE('glen')
        c2 = ListGetElementReal( ViscGlenExp_h,Element=Element,Found=Found)

        ! the second invariant is not taken from the strain rate tensor,
        ! but rather 2*strain rate tensor (that's why we divide by 4 = 2**2)        
        s(1:ngp) = ss(1:ngp)/4.0_dp

        c3 = ListGetElementReal( ViscCritical_h,Element=Element,Found=Found)
        IF( Found ) THEN
          c3 = c3**2
          WHERE( s(1:ngp) < c3 ) s(1:ngp) = c3
        END IF

        IF( ListGetElementLogical( ViscArrSet_h,Element,Found=Found) ) THEN
          ArrheniusFactor = ListGetElementReal( ViscArr_h,Element=Element)
          ViscVec(1:ngp) = 0.5_dp * (ArrheniusFactor)**(-1.0_dp/c2) * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp);                    
          
          IF( ViscNewton ) THEN
            WHERE( s(1:ngp) > c3 ) ViscDerVec(1:ngp) = 0.5_dp * ArrheniusFactor**(-1.0_dp/c2) &
                * ((1.0_dp/c2)-1.0_dp)/2.0_dp * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp - 1.0_dp)/4.0_dp
          END IF
        ELSE         
          ! lets for the time being have this hardcoded
          Tlimit = ListGetElementReal( ViscTlimit_h,Element=Element)
          A1 = ListGetElementReal( ViscRate1_h,Element=Element)
          A2 = ListGetElementReal( ViscRate2_h,Element=Element)
          Q1 = ListGetElementReal( ViscEne1_h,Element=Element)
          Q2 = ListGetElementReal( ViscEne2_h,Element=Element)

          ! WHERE is faster than DO + IF
          TempVec => ListGetElementRealVec( ViscTemp_h, ngp, BasisVec, Element )
          
          WHERE( TempVec(1:ngp ) < Tlimit )
            ArrheniusFactorVec(1:ngp) = A1 * EXP( -Q1/(R * (273.15_dp + TempVec(1:ngp))))
          ELSE WHERE( TempVec(1:ngp) > 0.0_dp ) 
            ArrheniusFactorVec(1:ngp) = A2 * EXP( -Q2/(R * (273.15_dp)))
          ELSE WHERE
            ArrheniusFactorVec(1:ngp) = A2 * EXP( -Q2/(R * (273.15_dp + TempVec(1:ngp))))
          END WHERE
          
          EhfVec => ListGetElementRealVec( ViscGlenFactor_h, ngp, BasisVec,Element=Element )
          ViscVec(1:ngp) = 0.5_dp * (EhFVec(1:ngp) * ArrheniusFactorVec(1:ngp))**(-1.0_dp/c2) * &
              s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp);
          
          IF( ViscNewton ) THEN
            WHERE( s(1:ngp) > c3 ) 
              ViscDerVec(1:ngp) = 0.5_dp * (  EhFVec(1:ngp) * ArrheniusFactorVec(1:ngp))**(-1.0_dp/c2) &
                    * ((1.0_dp/c2)-1.0_dp)/2.0_dp * s(1:ngp)**(((1.0_dp/c2)-1.0_dp)/2.0_dp - 1.0_dp)/4.0_dp
            END WHERE
          END IF                      
        END IF
        

      CASE('power law')
        c2 = ListGetElementReal( ViscExp_h,Element=Element)

        c3 = ListGetElementReal( ViscCritical_h,Element=Element,Found=Found)       
        IF( Found ) THEN
          c3 = c3**2
          WHERE( ss(1:ngp) < c3 ) ss(1:ngp) = c3
        END IF
        
        ViscVec(1:ngp) = ViscVec0(1:ngp) * ss(1:ngp)**((c2-1)/2)
       
        IF (ViscNewton ) THEN
          WHERE(ss(1:ngp) /= 0) ViscDerVec(1:ngp) = &
              ViscVec0(1:ngp) * (c2-1)/2 * ss(1:ngp)**((c2-1)/2-1)
        END IF
        
        c4 = ListGetElementReal( ViscNominal_h,Element=Element,Found=Found)
        IF( Found ) THEN
          ViscVec(1:ngp) = ViscVec(1:ngp) / c4**(c2-1)
          IF (ViscNewton ) THEN
            ViscDerVec(1:ngp) = ViscDerVec(1:ngp) / c4**(c2-1)
          END IF
        END IF

      CASE('power law too')
        c2 = ListGetElementReal( ViscExp_h,Element=Element)           
        ViscVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)* ss(1:ngp)**(-(c2-1)/(2*c2)) / 2

        IF (ViscNewton ) THEN
          ViscDerVec(1:ngp) = ViscVec0(1:ngp)**(-1/c2)*(-(c2-1)/(2*c2))*ss(1:ngp)*(-(c2-1)/(2*c2)-1) / 2
        END IF
                
      CASE ('carreau')      
        c1 = ListGetElementReal( ViscDiff_h,Element=Element)
        c2 = ListGetElementReal( ViscExp_h,Element=Element)
        c3 = ListGetElementReal( ViscTrans_h,Element=Element)
        c4 = ListGetElementReal( ViscYasuda_h,Element=Element,Found=Found)
        IF( Found ) THEN
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4) 
          
          IF( ViscNewton ) THEN
            ViscDerVec(1:ngp) = c1*(1+c3**c4*ss(1:ngp)**(c4/2))**((c2-1)/c4-1)*(c2-1)/2*c3**c4*&
                ss(1:ngp)**(c4/2-1)
          END IF
        ELSE
          ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * (1 + c3*c3*ss(1:ngp))**((c2-1)/2) 

          IF( ViscNewton ) THEN
            ViscDerVec(1:ngp) = c1*(c2-1)/2*c3**2*(1+c3**2*ss(1:ngp))**((c2-1)/2-1)
          END IF
        END IF

      CASE ('cross')
        c1 = ListGetElementReal( ViscDiff_h,Element=Element)
        c2 = ListGetElementReal( ViscExp_h,Element=Element)
        c3 = ListGetElementReal( ViscTrans_h,Element=Element)

        ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 / (1 + c3*ss(1:ngp)**(c2/2))

        IF( ViscNewton ) THEN
          ViscDerVec(1:ngp) = -c1*c3*ss(1:ngp)**(c2/2)*c2 / (2*(1+c3*ss(1:ngp)**(c2/2))**2*ss(1:ngp))
        END IF
          
      CASE ('powell eyring')
        c1 = ListGetElementReal( ViscDiff_h,Element=Element)
        c2 = ListGetElementReal( ViscTrans_h,Element=Element)

        s(1:ngp) = SQRT(ss(1:ngp))

        IF( ViscNewton ) THEN          
          WHERE( c2*s(1:ngp) < 1.0d-5 )
            ViscVec(1:ngp) = ViscVec0(1:ngp) + c1
            ViscDerVec(1:ngp) = 0.0_dp
          ELSE WHERE
            ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * LOG(c2*s(1:ngp)+SQRT(c2*c2*ss(1:ngp)+1))/(c2*ss(1:ngp))            
            ViscDerVec(1:ngp) = c1*(c2/(2*s(1:ngp))+c2**2/(2*SQRT(c2**2*ss(1:ngp)+1)))/ &
                ((c2*s(1:ngp)+SQRT(c2*ss(1:ngp)+1))*c2*s(1:ngp)) - &
                c1*LOG(c2*s(1:ngp)+SQRT(c2**2*ss(1:ngp)+1))/(c2*s(1:ngp)**3)/2
          END WHERE
        ELSE
          WHERE( c2*s(1:ngp) < 1.0d-5 )
            ViscVec(1:ngp) = ViscVec0(1:ngp) + c1
          ELSE WHERE
            ViscVec(1:ngp) = ViscVec0(1:ngp) + c1 * LOG(c2*s(1:ngp)+SQRT(c2*c2*ss(1:ngp)+1))/(c2*ss(1:ngp))            
          END WHERE
        END IF

      CASE DEFAULT 
        CALL Fatal(Caller,'Unknown material model')

      END SELECT

      IF( ViscNewton ) THEN
        IF(GotRelax) ViscDerVec(1:ngp) = NewtonRelax * ViscDerVec(1:ngp)
      END IF
      
      ! If requested, save viscosity field (on nodes, ip points or elements). 
      IF(SaveVisc) THEN
        IF( NSHandles(tid) % ViscVar % TYPE == Variable_on_nodes ) THEN
          DO i=1,n
            NSHandles(tid) % ViscVar % Values(NSHandles(tid) % ViscVar % Perm(Element % NodeIndexes(i))) = &
                NSHandles(tid) % ViscVar % Values(NSHandles(tid) % ViscVar % Perm(Element % NodeIndexes(i))) + &
                SUM(BasisVec(1:ngp,i)*ViscVec(1:ngp)*detJVec(1:ngp))
          END DO
        ELSE
          i = Element % ElementIndex
          IF( NSHandles(tid) % ViscVar % TYPE == Variable_on_gauss_points ) THEN
            j = NSHandles(tid) % ViscVar % Perm(i+1) - NSHandles(tid) % ViscVar % Perm(i) 
            IF(j /= ngp) THEN
              CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Viscosity" got '//I2S(ngp))
            END IF
            NSHandles(tid) % ViscVar % Values( &
                NSHandles(tid) % ViscVar % Perm(i)+1 : NSHandles(tid) % ViscVar % Perm(i+1) ) = ViscVec(1:ngp)
          ELSE
            NSHandles(tid) % ViscVar % Values(NSHandles(tid) % ViscVar % Perm(i)) = SUM(ViscVec(1:ngp)) / ngp
          END IF
        END IF
      END IF

      ! If requested, save normalization weight associated to viscosity (and shearrate).
      IF(SaveWeight) THEN
        DO i=1,n
          NSHandles(tid) % WeightVar % Values(NSHandles(tid) % WeightVar % Perm(Element % NodeIndexes(i))) = &
              NSHandles(tid) % WeightVar % Values(NSHandles(tid) % WeightVar % Perm(Element % NodeIndexes(i))) + &
              SUM(BasisVec(1:ngp,i)*detJVec(1:ngp))
        END DO
      END IF

      END ASSOCIATE

    END FUNCTION EffectiveViscosityVec
      

    
    !------------------------------------------------------------------------------
    ! A special subroutine for performing static condensation when the velocity
    ! (the first dim components of the solution) is augmented by a bubble part.
    ! The speciality is that the bubble part at the previous time level
    ! is retrieved to approximate the first time derivative in a consistent manner.
    ! This version works only with the BDF(1) method, as higher-order versions have
    ! not yet been implemented.    
    !------------------------------------------------------------------------------
    SUBROUTINE LCondensate( N, nb, dim, M, K, F, xprev, x, Element_id )
      !------------------------------------------------------------------------------
      USE LinearAlgebra
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N   ! The number of retained DOFs per scalar field
      INTEGER, INTENT(IN) :: nb  ! The number of eliminated DOFs per scalar field
      INTEGER, INTENT(IN) :: dim ! The number of bubble-augmented fields 
      REAL(KIND=dp), INTENT(IN) :: M(:,:)
      REAL(KIND=dp), INTENT(INOUT) :: K(:,:), F(:)
      REAL(KIND=dp), INTENT(IN), OPTIONAL :: x(:,:)     ! The solution without 
      ! bubbles
      REAL(KIND=dp), INTENT(IN), OPTIONAL :: xprev(:,:) ! The previous solution 
      ! without bubbles
      INTEGER, INTENT(IN), OPTIONAL :: Element_id       ! The element identifier

      ! This subroutine accesses also the real-valued arrays bx(:) and bxprev(:)
      ! that contain coefficients of the bubble basis functions

  !------------------------------------------------------------------------------
      LOGICAL :: ComputeBubblePart
      REAL(KIND=dp) :: Kbb(nb*dim,nb*dim), Fb(nb*dim)
      REAL(KIND=dp) :: Kbl(nb*dim,n*(dim+1)), Klb(n*(dim+1), nb*dim)

      REAL(KIND=dp) :: xl(n*(dim+1)), xlprev((n+nb)*(dim+1))

      INTEGER :: DOFs
      INTEGER :: p, q, Cdofs((dim+1)*n), Bdofs(dim*nb)
  !------------------------------------------------------------------------------
      ComputeBubblePart = PRESENT(x) .AND. PRESENT(xprev) .AND. PRESENT(Element_id)

      DOFs = dim + 1 

      ! Vectorize the input array x and
      ! create xlprev that contains the full previous solution including 
      ! the bubble DOFs. First insert the DOFs that are retained:
      q = 0
      DO p = 1,n
        DO i = 1,DOFs
          q = q + 1
          Cdofs(q) = DOFs*(p-1) + i
        END DO
      END DO

      ! Then the DOFs of the bubble part:
      q = 0
      DO p = 1,nb
        DO i = 1,dim
          q = q + 1
          Bdofs(q) = DOFs*(p-1) + i + n*DOFs
        END DO
      END DO

      ! The following only works for the BDF(1) method: 
      IF (ComputeBubblePart) THEN
        xlprev = 0
        q = 0
        DO p = 1,n
          DO i = 1,DOFs
            q = q + 1
            xl(q) = x(i,p)
            xlprev(cdofs(q)) = xprev(i,p) ! cdofs identity mapping?
          END DO
        END DO

        q = 0
        DO p = 1,nb
          DO i = 1,dim
            q = q + 1
            xlprev(bdofs(q)) = bxprev((Element_id-1)*dim*nb+q)
          END DO
        END DO

        K = K + M / dt
        F = F + MATMUL(M,xlprev) / dt
      END IF

      Kbb = K(Bdofs,Bdofs)
      Kbl = K(Bdofs,Cdofs)
      Klb = K(Cdofs,Bdofs)
      Fb  = F(Bdofs)

      CALL InvertMatrix( Kbb,Nb*dim )

      F(cdofs) = F(cdofs) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
      K(cdofs,cdofs) = &
          K(cdofs,cdofs) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )

      ! The bubble part evaluated for the current solution candidate: 
      IF (ComputeBubblePart) bx((Element_id-1)*dim*nb+1:Element_id*dim*nb) = &
          MATMUL(Kbb,Fb-MATMUL(Kbl,xl))
      !------------------------------------------------------------------------------
    END SUBROUTINE LCondensate
    !------------------------------------------------------------------------------

  END SUBROUTINE LocalBulkMatrix


!------------------------------------------------------------------------------
! Assemble local finite element matrix for a single boundary element and glue
! it to the global matrix.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix( Element, n, nd, dim, dt, SpecificLoad, InitHandles, &
      FrictionNewton)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, dim
    REAL(KIND=dp), INTENT(IN) :: dt
    LOGICAL, INTENT(IN) :: SpecificLoad, FrictionNewton
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------    
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), TARGET :: STIFF(nd*(dim+1),nd*(dim+1)), FORCE(nd*(dim+1))
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str, FSSAFlag
    INTEGER :: c,i,j,k,l,p,q,t,ngp,norm_comp,no_slip_comp
    LOGICAL :: NormalTangential, HaveSlip, HaveForce, HavePres, HaveFrictionW, HaveFrictionU, &
        HaveFriction, HaveNormal, FrictionNormal, Found, Stat, HaveFSSA, &
        FoundLoad, GotRelax, LocalNewton, HaveNormalSlip
    REAL(KIND=dp) :: ExtPressure, s, detJ, FSSAtheta, wut0, wexp, wcoeff, un, ut, rho
    REAL(KIND=dp) :: SlipCoeff(3), NormalSlipCoeff, SurfaceTraction(3), Normal(3), Tangent(3), Tangent2(3), &
        Vect(3), Velo(3), TanFrictionCoeff, DummyVals(1), LoadVec(dim), FSSAaccum, &
        NewtonRelax 
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: ExtPressure_h, SurfaceTraction_h, SlipCoeff_h, NormalSlipCoeff_h, &
        NormalTangential_h, NormalTangentialVelo_h, WeertmanCoeff_h, WeertmanExp_h, &
        FrictionUt0_h, FrictionNormal_h, FrictionCoeff_h, &
        FSSAtheta_h, Dens_h, Load_h(3), FSSAaccum_h
    TYPE(VariableHandle_t), SAVE :: Velo_v
    TYPE(Variable_t), POINTER, SAVE :: NrmSol, VeloSol
    TYPE(ValueList_t), POINTER :: BC    
    REAL(KIND=dp) :: TanFder,JAC(nd*(dim+1),nd*(dim+1)),SOL(nd*(dim+1)),NodalSol(dim+1,nd)
    TYPE(Variable_t), POINTER, SAVE :: SlipCoeffVar, SlipSpeedVar, SlipWeightVar
    LOGICAL, SAVE :: SaveSlipSpeed, SaveSlipCoeff, SaveSlipWeight
    
    SAVE Basis, HaveNormal, GotRelax, NewtonRelax
    
!------------------------------------------------------------------------------
    
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','Normal Surface Traction')
      IF( .NOT. ListGetElementSomewhere( ExtPressure_h) ) THEN
        CALL ListInitElementKeyword( ExtPressure_h,'Boundary Condition','External Pressure')      
      END IF
      CALL ListInitElementKeyword( SurfaceTraction_h,'Boundary Condition','Surface Traction',InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( SlipCoeff_h,'Boundary Condition','Slip Coefficient',InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( NormalSlipCoeff_h,'Boundary Condition','Normal Slip Coefficient')
      
      CALL ListInitElementKeyword( FrictionCoeff_h,'Boundary Condition','Friction Coefficient',&
          EvaluateAtIp=.TRUE., DummyCount=1)     
      CALL ListInitElementKeyword( FrictionNormal_h,'Boundary Condition','Friction Normal Velocity Zero')     
      CALL ListInitElementKeyword( FrictionUt0_h,'Boundary Condition','Friction Linear Velocity')
      IF(FrictionUt0_h % NotPresentAnywhere) THEN
        CALL ListInitElementKeyword( FrictionUt0_h,'Boundary Condition','Weertman Linear Velocity')
      END IF
        
      CALL ListInitElementKeyword( WeertmanCoeff_h,'Boundary Condition','Weertman Friction Coefficient')
      CALL ListInitElementKeyword( WeertmanExp_h,'Boundary Condition','Weertman Exponent')
      
      CALL ListInitElementKeyword( NormalTangentialVelo_h,'Boundary Condition',&
          'Normal-Tangential Velocity' )
      CALL ListInitElementKeyword( NormalTangential_h,'Boundary Condition',&
          'Normal-Tangential '//GetVarName(CurrentModel % Solver % Variable) )
      CALL ListInitElementKeyword( FSSAtheta_h, 'Boundary Condition',&
           'FSSA Theta')
      str = ListGetString( CurrentModel % Solver % Values,'Normal Vector Name',Found )
      IF(.NOT. Found) str = 'Normal Vector'
      NrmSol => VariableGet( CurrentModel % Solver % Mesh % Variables, str, ThisOnly = .TRUE.) 
      
      VeloSol => CurrentModel % Solver % Variable
      
      !CALL ListInitElementVariable( Normal_v, str, Found=HaveNormal)

      CALL ListInitElementVariable( Velo_v )
      CALL ListInitElementKeyword( Dens_h,'Material','Density')
      DO i=1,dim 
        CALL ListInitElementKeyword( Load_h(i),'Body Force','Flow Bodyforce '//I2S(i))
      END DO
      CALL ListInitElementKeyword( FSSAaccum_h,'Boundary Condition','FSSA Accumulation')

      NewtonRelax = ListGetCReal( CurrentModel % Solver % Values,&
          'Friction Newton Relaxation Factor',GotRelax )
      IF(.NOT. GotRelax) NewtonRelax = 1.0_dp
      
      SlipCoeffVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Coefficient',ThisOnly=.TRUE.)
      SaveSlipCoeff = ASSOCIATED(SlipCoeffVar)
      IF(SaveSlipCoeff) SlipCoeffVar % Values = 0.0_dp

      SaveSlipSpeed = .FALSE.
      IF( SaveSlipCoeff ) THEN
        SlipSpeedVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Speed',ThisOnly=.TRUE.)
        SaveSlipSpeed = ASSOCIATED(SlipSpeedVar)
        IF(SaveSlipSpeed) SlipSpeedVar % Values = 0.0_dp
        
        SlipWeightVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Weight',ThisOnly=.TRUE.)
        SaveSlipWeight = ASSOCIATED(SlipWeightVar)
        IF(SaveSlipWeight) SlipWeightVar % Values = 0.0_dp        
      END IF
      
      InitHandles = .FALSE.
    END IF
    
    BC => GetBC()
    IF( ALLOCATED( Basis ) ) THEN
      IF( SIZE( Basis ) < nd ) THEN
        DEALLOCATE( Basis ) 
      END IF
    END IF

    IF( .NOT. ALLOCATED( Basis ) ) THEN
      ALLOCATE( Basis(nd) )
    END IF
          
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    JAC = 0.0d0
    FORCE = 0.0d0
    c = dim + 1

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    ngp = IP % n
    
    LocalNewton = .FALSE.
    NormalTangential = ListGetElementLogical( NormalTangentialVelo_h, Element, Found )
    IF (.NOT.Found) THEN
      NormalTangential = ListGetElementLogical( NormalTangential_h, Element, Found )
    END IF

    FrictionNormal = .FALSE.
    TanFder=0._dp
    no_slip_comp = 0
    
    ! There is no elemental routine for this.
    ! So whereas this breaks the beauty it does not cost too much.
    FSSAFlag = GetString(BC, 'FSSA Flag', Found)
    IF (.NOT.Found) FSSAFlag = 'none'
    
    HaveFrictionW = ListCheckPresent( BC,'Weertman Friction Coefficient') 
    HaveFrictionU = ListCheckPresent( BC,'Friction Coefficient')
    HaveFriction = HaveFrictionU .OR. HaveFrictionW
    
    IF( HaveFriction ) THEN
      wut0 = ListGetElementReal( FrictionUt0_h, Element = Element )
      FrictionNormal = ListGetElementLogical( FrictionNormal_h, Element ) 
    END IF
    
    DO t=1,ngp
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t), detJ, Basis )

      s = detJ * IP % s(t)

      ! Given force on a boundary componentwise
      !----------------------------------------
      SurfaceTraction = ListGetElementReal3D( SurfaceTraction_h, Basis, Element, HaveForce, GaussPoint = t )

      ! Given force to the normal direction
      !------------------------------------
      ExtPressure = ListGetElementReal( ExtPressure_h, Basis, Element, HavePres, GaussPoint = t )

      ! Slip coefficient
      !----------------------------------
      SlipCoeff = ListGetElementReal3D( SlipCoeff_h, Basis, Element, HaveSlip, GaussPoint = t )
      NormalSlipCoeff = ListGetElementReal( NormalSlipCoeff_h, Basis, Element, HaveNormalSlip, GaussPoint = t )

      FSSAtheta = ListGetElementReal( FSSAtheta_h, Basis, Element, HaveFSSA, GaussPoint = t )
      IF (HaveFSSA) THEN
        IF (FSSAtheta == 0.0) HaveFSSA = .FALSE.
        rho = ListGetElementRealParent( Dens_h, Basis, Element, Found )
        IF (.NOT.Found .AND. t==1) THEN
          CALL WARN('IncompressibleNSSolver (FSSA)','"Density" in Parent element not found!')
          HaveFSSA = .FALSE.
        END IF
      END IF

      IF (HaveFSSA) THEN
        ! Flow bodyforce if present
        LoadVec = 0._dp
        FoundLoad = .FALSE.
        DO i=1,dim
          LoadVec(i) = ListGetElementRealParent( Load_h(i), Basis, Element, Found )
          FoundLoad = FoundLoad .OR. Found
          IF( Found .AND. .NOT.SpecificLoad) THEN
            LoadVec(i) = rho * LoadVec(i)
          END IF
        END DO
        IF (.NOT.FoundLoad) THEN
          CALL WARN('IncompressibleNSSolver (FSSA)','No component of "Flow Body Force" in Parent element found!')
          HaveFSSA = .FALSE.
        END IF
      END IF
      
      ! Nothing to do, exit the routine
      !---------------------------------
      IF(.NOT. (HaveForce .OR. HavePres .OR. HaveSlip .OR. HaveNormalSlip .OR. HaveFriction .OR. HaveFSSA)) RETURN

      ! Calculate normal vector only if needed, which is almost always...
      IF( HavePres .OR. NormalTangential .OR. HaveFriction  .OR. HaveFSSA .OR. HaveNormalSlip ) THEN
        Normal = ConsistentNormalVector( CurrentModel % Solver, NrmSol, Element, Found, Basis = Basis )
        IF(.NOT. Found) Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t),.TRUE. )

        ! Define normal component so we can differentiate the treatment of normal and tangent directions
        ! if needed. If n-t coordinates are not used we take the dominating direction assuming cartesian system.
        norm_comp = 1
        IF( .NOT. NormalTangential) THEN
          DO i=2,dim
            IF(ABS(Normal(i)) > ABS(Normal(norm_comp))) norm_comp = i
          END DO
        END IF

        ! If we have normal slip given we now now to which component it affects.
        ! This overwrites potential "slip coefficient" if given for this component!
        IF(HaveNormalSlip) THEN
          SlipCoeff(norm_comp) = NormalSlipCoeff
        END IF
      END IF

      ! Create friction model if "slip coefficient" is not given.
      ! Because of backward compatibility we cannot allow for both slip and
      ! friction at the same time.
      !---------------------------------------------------------------------
      IF( HaveFriction  .AND. .NOT. HaveSlip ) THEN
        ! Velocity at integration point for nonlinear friction laws
        Velo = ListGetElementVectorSolution( Velo_v, Basis, Element, dofs = dim )

        ! It seems futile to take normal component away as it is usually zero by construction!
        ! We include here the option not to remove it. 
        IF(.NOT. FrictionNormal ) THEN
          un = SUM( Normal(1:dim) * Velo(1:dim) )
          velo(1:dim) = velo(1:dim)-un*normal(1:dim)
        END IF
        ut = MAX(wut0, SQRT(SUM(Velo(1:dim)**2)))
                     
        IF( HaveFrictionW ) THEN
          ! Weertman friction law computed internally
          wcoeff = ListGetElementReal( WeertmanCoeff_h, Basis, Element, GaussPoint = t )
          wexp = ListGetElementReal( WeertmanExp_h, Basis, Element, GaussPoint = t )
          TanFrictionCoeff = MIN(wcoeff * ut**(wexp-1.0_dp),1.0e20)
          ! dTanFrictionCoeff/dut for Newton
          IF(FrictionNewton ) THEN
            TanFder=0._dp
            IF ((ut > wut0).AND.(TanFrictionCoeff < 1.0e20)) &
                TanFder = (wexp-1.0_dp) * wcoeff * ut**(wexp-2.0_dp)
            LocalNewton = .TRUE.
          END IF
        ELSE
          ! Else, user defined friction law
          DummyVals(1) = ut          
          TanFrictionCoeff = ListGetElementReal( FrictionCoeff_h, Basis, Element, &
              GaussPoint = t, DummyVals = DummyVals )             
        END IF

        ! We Do not set the slip in normal direction if given by friction model only. 
        DO i=1,dim
          IF(i /= norm_comp) SlipCoeff(i) = TanFrictionCoeff
        END DO
         
        IF(.NOT. HaveNormalSlip) THEN
          no_slip_comp = norm_comp
        END IF

        ! We initially have friction given and translate it to terms of slip coefficient.
        ! The code is then continued as if we would have slip coefficients. 
        HaveSlip = .TRUE.
        
        IF(SaveSlipSpeed) THEN
          DO i=1,n
            j = SlipSpeedVar % Perm(Element % NodeIndexes(i))
            IF(j>0) THEN
              SlipSpeedVar % Values(j) = SlipSpeedVar % Values(j) + &
                  Basis(i) * IP % s(t) * ut
            END IF
          END DO
        END IF
      END IF

      IF(SaveSlipCoeff) THEN
        DO i=1,n
          j = SlipCoeffVar % Perm(Element % NodeIndexes(i))
          IF(j>0) THEN
            SlipCoeffVar % Values(j) = SlipCoeffVar % Values(j) + &
                Basis(i) * IP % s(t) * MAXVAL(SlipCoeff(1:dim))
          END IF
        END DO
        IF(SaveSlipWeight) THEN
          DO i=1,n
            j = SlipWeightVar % Perm(Element % NodeIndexes(i))
            IF(j>0) SlipWeightVar % Values(j) = SlipWeightVar % Values(j) + Basis(i) * IP % s(t)
          END DO
        END IF
      END IF

      ! Project external pressure to the normal direction
      IF( HavePres ) THEN
        IF( NormalTangential ) THEN
          SurfaceTraction(1) = SurfaceTraction(1) + ExtPressure
        ELSE
          SurfaceTraction = SurfaceTraction + ExtPressure * Normal
        END IF
        HaveForce = .TRUE. 
      END IF
      
      ! Calculate directions for N-T system
      IF( NormalTangential ) THEN       
        SELECT CASE( dim ) 
        CASE(2)
          Tangent(1) =  Normal(2)
          Tangent(2) = -Normal(1)
          Tangent(3) =  0.0_dp
          Tangent2   =  0.0_dp
        CASE(3)
          CALL TangentDirections( Normal, Tangent, Tangent2 ) 
        END SELECT
      END IF

      ! Assemble the slip coefficients to the stiffness matrix
      IF( HaveSlip ) THEN               
        IF ( NormalTangential ) THEN
          DO i=1,dim
            ! We know by construction that if only friction models are given there is no normal slip. 
            IF(i==no_slip_comp) CYCLE
            SELECT CASE(i)
            CASE(1)
              Vect = Normal
            CASE(2)
              Vect = Tangent
            CASE(3)
              Vect = Tangent2
            END SELECT            
            
            DO p=1,nd
              DO q=1,nd               
                DO j=1,dim
                  DO k=1,dim
                    STIFF( (p-1)*c+j,(q-1)*c+k ) = &
                        STIFF( (p-1)*c+j,(q-1)*c+k ) + &
                        s * SlipCoeff(i) * Basis(q) * Basis(p) * Vect(j) * Vect(k)

                    IF(LocalNewton) THEN
                      ! Only tangential directions have Newton linearization
                      IF(i==norm_comp) CYCLE
                      JAC((p-1)*c+j,(q-1)*c+k ) = JAC((p-1)*c+j,(q-1)*c+k ) + &
                          s * TanFder * Basis(q) * Basis(p) * Vect(j) * velo(k) * SUM(velo(1:dim)*Vect(1:dim))/ut
                    END IF

                  END DO
                END DO
              END DO
            END DO
          END DO
        ELSE       
          DO p=1,nd
            DO q=1,nd
              DO i=1,dim
                IF(i == no_slip_comp) CYCLE
                STIFF( (p-1)*c+i,(q-1)*c+i ) = &
                    STIFF( (p-1)*c+i,(q-1)*c+i ) + &
                    s * SlipCoeff(i) * Basis(q) * Basis(p)

               IF(LocalNewton) THEN
                 DO j=1,dim
                  IF(i==norm_comp) CYCLE
                  JAC((p-1)*c+i,(q-1)*c+j ) = JAC((p-1)*c+i,(q-1)*c+j ) + &
                        s * TanFder * Basis(q) * Basis(p) * velo(i) * velo(j) / ut
                 END DO
               END IF

              END DO
            END DO
          END DO
        END IF
      END IF

           
      ! Assemble given forces to r.h.s.
      IF( HaveForce .OR. HavePres .OR. HaveFSSA ) THEN
        FSSAaccum = ListGetElementReal( FSSAaccum_h, Basis, Element, GaussPoint = t )        
        IF ( NormalTangential ) THEN
          DO i=1,dim
            SELECT CASE(i)
            CASE(1)
              Vect = Normal
            CASE(2)
              Vect = Tangent
            CASE(3)
              Vect = Tangent2
            END SELECT

            DO q=1,nd
              DO j=1,dim
                l = (q-1)*c + j
                FORCE(l) = FORCE(l) + s * Basis(q) * SurfaceTraction(i) * Vect(j)
                IF(HaveFSSA) THEN
                  IF (i==1) FORCE(l) = FORCE(l) + s * FSSAtheta * dt * FSSAaccum * Basis(q) * LoadVec(j) 
                END IF
              END DO
            END DO
          END DO
        ELSE       
          DO i=1,dim
            DO q=1,nd
              k = (q-1)*c + i
              FORCE(k) = FORCE(k) + s * Basis(q) * SurfaceTraction(i)
              IF(HaveFSSA) THEN
                FORCE(k) = FORCE(k) + s * FSSAtheta * dt *  FSSAaccum * Basis(q) * LoadVec(i)
              END IF
            END DO
          END DO
        END IF
      END IF
      
      !FSSA stabilization: sum_i <u_i*n_i*v_z> * StabCoeff, i=1,..,dim
      !---------------------------------------------------------------
     
      IF ( HaveFSSA ) THEN
        SELECT CASE(FSSAFlag)
          ! version 1,  approximation with normal pointing into z-direction
        CASE ('normal')
          DO p=1,nd
            DO q=1,nd
              STIFF( (p-1)*c+dim,(q-1)*c+dim ) = &
                  STIFF( (p-1)*c+dim,(q-1)*c+dim )  &
                  - s * FSSAtheta * dt * LoadVec(dim) * Basis(q) * Basis(p) * Normal(dim)
            END DO
          END DO
          ! version 2, transposed
        CASE ('transposed')
          DO p=1,nd
            DO q=1,nd
              DO i=1,dim
                STIFF( (p-1)*c+i,(q-1)*c+dim ) = & 
                     STIFF( (p-1)*c+i,(q-1)*c+dim )  &
                     - s * FSSAtheta * dt * LoadVec(dim) * Basis(q) * Basis(p) * Normal(i)
              END DO
            END DO
          END DO        
        CASE ('full') ! full entry matrix FSSA
          DO p=1,nd
            DO q=1,nd
              DO i=1,dim
                DO j=1,dim
                  STIFF( (p-1)*c+j,(q-1)*c+i ) = & 
                     STIFF( (p-1)*c+j,(q-1)*c+i )  &
                     - s * FSSAtheta * dt * LoadVec(j) * Basis(q) * Basis(p) * Normal(i)
                !PRINT *, "K(",p,q,i,")=", s * FSSAtheta * Basis(q) * Basis(p) * Normal(i), STIFF( (p-1)*c+dim,(q-1)*c+i )
                END DO
              END DO
            END DO
          END DO          
        CASE DEFAULT

        END SELECT

      END IF
    END DO


    IF(LocalNewton) THEN
      CALL GetLocalSolution( NodalSol )
      SOL=0._dp
      DO i = 1, c
        SOL(i::c) = NodalSol(i,1:nd)
      END DO

      IF(GotRelax) JAC = NewtonRelax * JAC
      
      STIFF=STIFF+JAC
      FORCE=FORCE + MATMUL(JAC,SOL)
    END IF
      
    CALL DefaultUpdateEquations( STIFF, FORCE )
        
  END SUBROUTINE LocalBoundaryMatrix

  
END MODULE IncompressibleLocalForms


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  LOGICAL :: Found, Serendipity

  ! Equal-order stabilisation asks for the opposite element to everything below:
  ! a plain "p:1" with no bubble at all, the velocity and the pressure carried on
  ! the same nodal basis and the pressure mode held by the Brezzi-Pitkaranta term
  ! in LocalBulkMatrix rather than by an inf-sup stable pair.
  !
  ! It has to be settled HERE and not in _Init. Def_Dofs is filled from the
  ! "Element" string while the sif is read, before the mesh is loaded, and _Init0
  ! is the only entry point that runs earlier (ModelDescription: the _Init0 sweep
  ! precedes PrepareMesh). Adding the bubble string and removing it later would be
  ! too late -- the bubble dofs are already counted into the mesh by then, which
  ! is visible as EnlargeCoordinates growing the node array by nb per element.
  !
  ! ListAddNewString throughout, so an explicit sif "Element" still wins over
  ! either branch. What that can produce -- stabilisation on top of a bubble -- is
  ! refused in the solver, where the element is known.
  !
  ! "n:1" and not "p:1". Both give one dof per node and no bubble, but "p:1"
  ! makes it a p-element, and then every basis evaluation goes through the
  ! hierarchic path and GaussPointsAdapt through the p-element rule -- paid on
  ! every element, every iteration, for a degree-one basis that the nodal path
  ! already describes exactly. Reaching "n:1" needs the PDefs guard on the
  ! GaussPointsAdapt call in the solver; without it this segfaults.
  IF( GetLogical( GetSolverParams(), 'Pressure Stabilization', Found ) ) THEN
    CALL ListAddNewString(GetSolverParams(),'Element','n:1')
    RETURN
  END IF

  Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found)
  IF(.NOT.Found) Serendipity = .TRUE.
  
  IF(Serendipity) THEN
    CALL ListAddNewString(GetSolverParams(),'Element', &
      'p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4')
  ELSE
    CALL ListAddNewString(GetSolverParams(),'Element', &
      'p:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:4 -pyramid b:4')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE IncompressibleNSSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver_init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params 
  LOGICAL :: Found
  INTEGER :: dim
  CHARACTER(:), ALLOCATABLE :: str
  CHARACTER(*), PARAMETER :: Caller = 'IncompressibleNSSolver_init'
!------------------------------------------------------------------------------ 
  Params => GetSolverParams() 

  IF( ListCheckPresentAnyBC( Model, 'Pressure 1' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 1< instead of >Pressure 1<')
  END IF
  IF( ListCheckPresentAnyBC( Model, 'Pressure 2' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 3< instead of >Pressure 2<')
  END IF
  IF( ListCheckPresentAnyBC( Model, 'Pressure 3' ) ) THEN
    CALL Fatal( Caller,'Use >Surface Traction 3< instead of >Pressure 3<')
  END IF
  
  dim = CoordinateSystemDimension()
  
  IF ( dim == 2 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow Solution[Velocity:2 Pressure:1]')
  ELSE
    CALL ListAddNewString(Params, 'Variable', &
        'Flow Solution[Velocity:3 Pressure:1]')
  END IF

  ! Study only velocity components in linear system
  CALL ListAddNewInteger(Params, 'Nonlinear System Norm DOFs', dim )

  ! This should be true to incompressible flows where pressure level is not uniquely determined
  CALL ListAddNewLogical(Params, 'Relative Pressure Relaxation', .TRUE. )

  ! Automate the choice for the variational formulation:
  CALL ListAddNewLogical(Params, 'GradP Discretization', .FALSE.)
  CALL ListAddNewLogical(Params, 'Div-Curl Discretization', .FALSE.)

  
  ! It makes sense to eliminate the bubbles to save memory and time
  CALL ListAddNewLogical(Params, 'Bubbles in Global System', .FALSE.)

  ! The recovery of transient bubble DOFs is done such that at least two 
  ! iterations within the same time step are needed. The multiple solutions are
  ! ensured by making at least two nonlinear iterations. However, if "steady 
  ! state" iterations are performed, computing just one nonlinear iterate 
  ! is sufficient.
  IF ( .NOT. ListGetLogical(Params, 'Bubbles In Global System', Found) ) THEN
    CALL ListAddNewInteger(Params, 'Nonlinear System Min Iterations', 2)
  END IF

  ! Backward compatibility with old FlowSolver
  str = GetString( Params, 'Flow Model', Found )
  IF( Found ) THEN
    SELECT CASE(str)
    CASE('no convection')
      CALL Warn(Caller,'Option "Flow Model = no convection" not used in this Solver!')
    CASE('stokes')
      CALL ListAddNewLogical( Params,'Stokes Flow',.TRUE.)
    CASE DEFAULT
    END SELECT
  END IF

  IF( GetLogical( Params,'Save Viscosity', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Viscosity')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Shearrate')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-nooutput Viscosity Weight')
  END IF

  IF( GetLogical( Params,'Save Slip', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Slip Coefficient')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'Slip Speed')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-nooutput Slip Weight')
  END IF

  CALL ListAddNewLogical(Params,'schur: Variable Output',.FALSE.)

  
!------------------------------------------------------------------------------ 
END SUBROUTINE IncompressibleNSSolver_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE IncompressibleLocalForms
  USE MainUtils

  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(GaussIntegrationPoints_t) :: IP

  INTEGER :: Element_id
  INTEGER :: i, n, nb, nd, nbdofs, dim, Active, maxiter, iter, nthr
  INTEGER :: stimestep = -1
 
  REAL(KIND=dp) :: Norm

  LOGICAL :: AllocationsDone = .FALSE., Found, StokesFlow, BlockPrec, Converged
  LOGICAL :: GradPVersion, DivCurlForm, SpecificLoad, InitBCHandles
  LOGICAL :: PStab
  REAL(KIND=dp) :: PStabCoeff
  INTEGER :: PStabHmode
  CHARACTER(LEN=MAX_NAME_LEN) :: PStabStr

  TYPE(Solver_t), POINTER, SAVE :: SchurSolver => Null()
  
  CHARACTER(*), PARAMETER :: Caller = 'IncompressibleNSSolver'

  SAVE AllocationsDone, stimestep

!------------------------------------------------------------------------------
! Local variables to be accessed by the contained subroutines:
!------------------------------------------------------------------------------
  LOGICAL :: LinearAssembly, Newton
!------------------------------------------------------------------------------ 

  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()

  ! Per-thread ValueHandle_t/cache storage for LocalBulkMatrix and
  ! EffectiveViscosityVec (see NSHandles_t, IncompressibleLocalForms module scope).
  ! Not called concurrently from more than one thread, so no locking needed here.
  IF (.NOT. ALLOCATED(NSHandles)) THEN
    nthr = 1
    !$ nthr = OMP_GET_MAX_THREADS()
    ALLOCATE(NSHandles(nthr))
  END IF

  !-----------------------------------------------------------------------------
  ! Allocate some permanent storage, this is done first time only:
  !-----------------------------------------------------------------------------
  IF (.NOT. AllocationsDone .AND. Transient .AND. &
      Mesh % MaxBDOFs > 0) THEN
    !
    ! Allocate arrays having a sufficient size for listing all bubble entries of
    ! the velocity solution (current and previous). These are needed in order to
    ! evaluate the time derivative of the bubble part.
    !
    nbdofs = Mesh % MaxBDOFs*dim*GetNOFActive()
    ALLOCATE(bx(nbdofs), bxprev(nbdofs)); 
    bx=0.0_dp; bxprev=0.0_dp

    AllocationsDone = .TRUE.
  END IF

  ! Check if the previous bubble part (bxprev) needs to be updated:
  IF (Transient .AND. GetTimestep() /= stimestep .AND. &
      Mesh % MaxBDOFs > 0) THEN
    bxprev = bx
    stimestep = GetTimestep()
  END IF

  Params => GetSolverParams() 

  !-----------------------------------------------------------------------------
  ! Output the number of integration points as information.
  ! This in not fully informative if several element types are present.
  !-----------------------------------------------------------------------------
  Element => Mesh % Elements( Solver % ActiveElements(1) ) 
  ! PReferenceElement only where there is a p-element to take the reference from.
  ! Asking for it unconditionally reads Element % PDefs, which a plain nodal
  ! element does not have, and segfaults in GaussPoints. It never showed because
  ! _Init0 always handed this solver a "p:1 ..." definition; an explicit
  ! "Element = n:1" in a sif is enough to reach it.
  IP = GaussPointsAdapt( Element, PReferenceElement = isActivePElement(Element) )
  CALL Info('IncompressibleNSSolver', &
      'Number of 1st integration points: '//I2S(IP % n), Level=5)
  
  !-----------------------------------------------------------------------------
  ! Set the flags/parameters which define how the system is assembled: 
  !-----------------------------------------------------------------------------
  LinearAssembly = GetLogical(Params, 'Linear Equation', Found )
  StokesFlow = GetLogical(Params, 'Stokes Flow', Found )
  GradPVersion = GetLogical(Params, 'GradP Discretization', Found)
  DivCurlForm = GetLogical(Params, 'Div-Curl Discretization', Found)
  SpecificLoad = GetLogical(Params,'Specific Load',Found)

  ! This solver assembles the Cartesian equations and only those. There is no
  ! metric anywhere in LocalBulkMatrix, no r weighting of the integration and no
  ! hoop term, so an axisymmetric or cylindrical case would not be approximated
  ! here, it would be answered with the equations of a different problem. FlowSolve
  ! handles those by dispatching to NavierStokesCylindricalCompose instead, a
  ! separate assembly with the metric in it, and this solver has no counterpart.
  !
  ! Nothing about the sif makes that visible -- the run simply produces a number --
  ! so it is refused here rather than discovered later.
  IF( CurrentCoordinateSystem() /= Cartesian ) CALL Fatal(Caller, &
      'This solver assembles the Cartesian equations only. Axisymmetric and '// &
      'cylindrical coordinates need the metric terms, which are not here: use '// &
      'FlowSolve, which dispatches to NavierStokesCylindrical for them')

  !-----------------------------------------------------------------------------
  ! Equal-order velocity/pressure held by a stabilised pressure block instead of
  ! an inf-sup stable pair. Opt-in per solver and never selected by default: it
  ! is a different discretisation, not an optimisation, and the element it needs
  ! was already chosen back in _Init0 on the strength of this same keyword.
  !-----------------------------------------------------------------------------
  PStab = GetLogical(Params, 'Pressure Stabilization', Found)
  PStabCoeff = GetConstReal(Params, 'Pressure Stabilization Coefficient', Found)
  IF (.NOT. Found) PStabCoeff = 1.0_dp

  ! Which element size enters tau ~ h^2/mu. The default reproduces the constants
  ! as they were measured; the other two exist because hK is not one quantity
  ! across families, so a case spanning several of them -- or one whose elements
  ! are far from cubic -- may want a measure it can reason about.
  PStabHmode = 0
  PStabStr = ListGetString(Params, 'Pressure Stabilization Element Size', Found)
  IF (Found) THEN
    SELECT CASE( PStabStr )
    CASE( 'shortest edge' )
      PStabHmode = 0
    CASE( 'longest edge' )
      PStabHmode = 1
    CASE( 'volume' )
      PStabHmode = 2
    CASE DEFAULT
      CALL Fatal(Caller, 'Unknown "Pressure Stabilization Element Size": '// &
          TRIM(PStabStr)//'. Use "shortest edge", "longest edge" or "volume"')
    END SELECT
    IF (PStab) CALL Info(Caller, 'Stabilization element size: '//TRIM(PStabStr), Level=5)
  END IF

  IF (PStab) THEN
    CALL Info(Caller, 'Using pressure stabilisation instead of bubbles', Level=5)
    ! Both together would stabilise an already stable pair, which is not a
    ! sharper answer but a softer one -- the bubble supplies the inf-sup
    ! condition and the tau then only adds a consistent but unnecessary
    ! perturbation. It can only arise from a sif overriding the _Init0 element.
    IF ( GetElementNOFBDOFs( Mesh % Elements( Solver % ActiveElements(1) ) ) > 0 ) &
        CALL Fatal(Caller, '"Pressure Stabilization" is the equal-order alternative '// &
        'to the bubble: remove the "b:" from "Element", or the keyword')
    ! Stage-limited on purpose. The residual assembled in LocalBulkMatrix carries
    ! the pressure gradient and the load, not the convective or transient terms,
    ! so it is a consistent scheme for Stokes and an inconsistent one for
    ! Navier-Stokes -- and it stabilises the pressure only, never the advection.
    ! Navier-Stokes is carried by the SUPG/PSPG terms in LocalBulkMatrix. What is
    ! still missing from the residual is rho*du/dt, so a run that actually carries
    ! that term would be stabilising the wrong operator -- refused rather than
    ! assembled quietly.
    !
    ! The test is NOT on Transient alone. Under "Stokes Flow" the assembly takes
    ! the branch that never reaches Default1stOrderTime and never fills MASS, so
    ! the momentum equation has no time derivative even inside a transient
    ! simulation: each step poses a quasi-static Stokes problem and the residual
    ! here is complete. That is the normal shape of a prognostic glaciological
    ! run -- ice creeps, the free surface evolves -- and FSSA_perlin2d is one, so
    ! refusing on Transient would have locked out the case this is most wanted
    ! for.
    IF (Transient .AND. .NOT. StokesFlow) CALL Fatal(Caller, &
        '"Pressure Stabilization" needs "Stokes Flow" in a transient run: '// &
        'rho*du/dt is not part of the stabilised residual')
  END IF
  
  Maxiter = GetInteger(Params, 'Nonlinear system max iterations', Found)
  IF (.NOT.Found) Maxiter = 1
  !-----------------------------------------------------------------------------

  IF (DivCurlForm) CALL Info(Caller, 'The div-curl form is used for the viscous terms')
  IF (GradPVersion) CALL Info(Caller, 'The pressure gradient is not integrated by parts')

  BlockPrec = GetLogical(Params,'Block Preconditioner',Found )
  IF(BlockPrec) THEN
    ! Do not create the Schur complement approximation if some other means are used. 
    IF( ListGetLogical( Params,'Create Schur Matrix Approximation',Found ) ) THEN
      CALL Info(Caller,'Schur complement not created, matrix approximation used instead!')
      BlockPrec = .FALSE.
    END IF
    IF( ListCheckPrefix( Params,'Schur Operator' ) ) THEN
      CALL Info(Caller,'Schur complement not created, matrix operations used instead!')
      BlockPrec = .FALSE.
    END IF
  END IF
  IF(BlockPrec ) THEN
    CALL Info(Caller,'Creating pressure block for block preconditioner')

    ! Create solver related to variable "schur" when using block preconditioning
    ! These keywords ensure that the matrix is truly used in the library version of the
    ! block solver.
    CALL ListAddNewString( Params,'Block Matrix Schur Variable','schur')

    ! Create solver that only acts as a container for the schur complement
    ! matrix used in the block preconditioning solver of the library.
    IF( .NOT. ASSOCIATED( SchurSolver ) ) THEN
      SchurSolver => CreateChildSolver( Solver,'schur', 1,'schur:') 
    END IF
  END IF
  
  
  
  DO iter=1,maxiter

    CALL Info(Caller,'--------------------------------------------------------', Level=4)
    WRITE( Message,'(A,I4)') 'Nonlinear iteration:', Iter
    CALL Info(Caller, Message, Level=4)
    CALL Info(Caller,'--------------------------------------------------------', Level=4)

100 CONTINUE
    
    Active = GetNOFActive()
    CALL DefaultInitialize()
    IF (ASSOCIATED(SchurSolver)) THEN
      CALL DefaultInitialize(USolver=SchurSolver)
    END IF
    
    Newton = GetNewtonActive( Solver )

    DO Element_id=1,1
      Element => GetActiveElement(Element_id)
      n  = GetElementNOFNodes(Element)
      !
      ! When the number of bubbles is obtained with the Update=.TRUE. flag,
      ! we need to call GetElementNOFBDOFs before calling GetElementNOFDOFs.
      !
      nb = GetElementNOFBDOFs(Element)
      nd = GetElementNOFDOFs(Element)
      
      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(Element, n, nd, nd+nb, dim,  DivCurlForm, GradPVersion, &
          SpecificLoad, StokesFlow, dt, LinearAssembly, nb, Newton, Transient, .TRUE., &
          SchurSolver, PStab, PStabCoeff, PStabHmode )
    END DO

    ! The serial call above (element 1) resolved NSHandles(1) — the material/BF
    ! keyword lookups and the ShearVar/ViscVar/WeightVar zeroing. Broadcast that
    ! resolved state to every other thread's slot before the parallel bulk loop,
    ! since each thread there calls LocalBulkMatrix with InitHandles=.FALSE. and
    ! only ever touches its own NSHandles(tid).
    IF (nthr > 1) NSHandles(2:nthr) = NSHandles(1)

    !$OMP PARALLEL SHARED(Active, dim, SpecificLoad, StokesFlow, &
    !$OMP                 DivCurlForm, GradPVersion, &
    !$OMP                 dt, LinearAssembly, Newton, Transient, SchurSolver, &
    !$OMP                 PStab, PStabCoeff, PStabHmode ) &
    !$OMP PRIVATE(Element, Element_id, n, nd, nb)  DEFAULT(None)
    !$OMP DO    
    DO Element_id=2,Active
      Element => GetActiveElement(Element_id)
      n  = GetElementNOFNodes(Element)
      nb = GetElementNOFBDOFs(Element, Update=.TRUE.)
      nd = GetElementNOFDOFs(Element)
      
      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(Element, n, nd, nd+nb, dim,  DivCurlForm, GradPVersion, &
          SpecificLoad, StokesFlow, dt, LinearAssembly, nb, Newton, Transient, .FALSE.,&
          SchurSolver, PStab, PStabCoeff, PStabHmode )
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    
    CALL DefaultFinishBulkAssembly()
    
    Active = GetNOFBoundaryElements()
    InitBCHandles = .TRUE.  
    
    DO Element_id=1,Active
      Element => GetBoundaryElement(Element_id)
      IF (ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()

        ! Skip 101 elements in 2D, and additionally 202's in 3D.
        IF ( GetElementFamily() < dim ) CYCLE        
        
        ! Get element local matrix and rhs vector:
        !-----------------------------------------
        CALL LocalBoundaryMatrix(Element, n, nd, dim, dt, SpecificLoad, InitBCHandles, Newton)
        InitBCHandles = .FALSE.
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    
    ! This is a matrix level routine for setting friction such that tangential
    ! traction is the normal traction multiplied by a coefficient.
    CALL SetImplicitFriction(Model, Solver,'Implicit Friction Coefficient')
    
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    IF(ASSOCIATED(SchurSolver)) CALL DefaultDirichletBCs(USolver=SchurSolver)

    ! Check stepsize for nonlinear iteration
    !------------------------------------------------------------------------------
    IF( DefaultLinesearch( Converged ) ) GOTO 100
    IF( Converged ) EXIT
    
    Norm = DefaultSolve()

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  BLOCK
    TYPE(Variable_t), POINTER, SAVE :: pVar, wVar
    REAL(KIND=dp) :: minw
    minw = 1.0e-20
    DO i=1,4
      SELECT CASE(i)
      CASE( 1 ) 
        wVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Weight',ThisOnly=.TRUE.)
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Coefficient',ThisOnly=.TRUE.)
      CASE( 2 ) 
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Slip Speed',ThisOnly=.TRUE.)
      CASE( 3 )         
        wVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity Weight',ThisOnly=.TRUE.)
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
      CASE( 4 ) 
        pVar => VariableGet( CurrentModel % Mesh % Variables,'Strainrate',ThisOnly=.TRUE.)
      END SELECT
      IF(ASSOCIATED(wVar) .AND. ASSOCIATED(pVar) ) THEN     
        CALL Info('IncompressibleNSSolver','Normalizing field number: '//I2S(i),Level=15)
        WHERE(wVar % Values > minw )
          pVar % Values = pVar % Values / wVar % Values
        END WHERE
      END IF
    END DO
  END BLOCK

  
  CALL DefaultFinish()
  
  CALL Info( Caller,'All done',Level=10)
!------------------------------------------------------------------------------
END SUBROUTINE IncompressibleNSSolver
!------------------------------------------------------------------------------
