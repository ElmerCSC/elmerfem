!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
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

CONTAINS

!------------------------------------------------------------------------------
! Assemble local finite element matrix for a single bulk element and glue
! it to the global matrix. Routine is vectorized and multithreaded.
!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(Element, n, nd, ntot, dim, &
       DivCurlForm, GradPVersion, SpecificLoad, StokesFlow, &
       dt, LinearAssembly, nb, Newton, Transient, InitHandles, SchurSolver )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, ntot, dim, nb
    LOGICAL, INTENT(IN) :: DivCurlForm, GradPVersion, SpecificLoad, StokesFlow
    REAL(KIND=dp), INTENT(IN) :: dt   
    LOGICAL, INTENT(IN) :: LinearAssembly, Newton, Transient, InitHandles
    TYPE(Solver_t), POINTER :: SchurSolver
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    LOGICAL :: Stat, Found

    REAL(KIND=dp), TARGET :: MASS(ntot*(dim+1),ntot*(dim+1)), &
        STIFF(ntot*(dim+1),ntot*(dim+1)), FORCE(ntot*(dim+1))
    REAL(KIND=dp) :: NodalSol(dim+1,ntot)
    REAL(KIND=dp) :: PrevNodalSol(dim+1,ntot)
    REAL(KIND=dp) :: s, rho

    REAL(KIND=dp), ALLOCATABLE, SAVE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:), &
        rhoVec(:), VeloPresVec(:,:), loadAtIpVec(:,:), VelocityMass(:,:), &
        PressureMass(:,:), ForcePart(:), &
        weight_a(:), weight_b(:), weight_c(:), tauVec(:), PrevTempVec(:), PrevPressureVec(:), &
        VeloVec(:,:), PresVec(:), GradVec(:,:,:)
    REAL(KIND=dp), POINTER :: muVec(:), LoadVec(:)
    REAL(KIND=dp), ALLOCATABLE :: muDerVec0(:),g(:,:,:),StrainRateVec(:,:,:)
    REAL(kind=dp) :: stifford(ntot,ntot,dim+1,dim+1), jacord(ntot,ntot,dim+1,dim+1), &
        JAC(ntot*(dim+1),ntot*(dim+1) )

    INTEGER :: t, i, j, k, p, q, ngp, allocstat
    INTEGER, SAVE :: elemdim
    INTEGER :: DOFs

    TYPE(ValueHandle_t), SAVE :: Dens_h, Load_h(3)
    
!DIR$ ATTRIBUTES ALIGN:64 :: BasisVec, dBasisdxVec, DetJVec, rhoVec, VeloPresVec, loadAtIpVec
!DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE, weight_a, weight_b, weight_c
!$OMP THREADPRIVATE(BasisVec, dBasisdxVec, DetJVec, rhoVec, VeloPresVec, loadAtIpVec, ElemDim )
!$OMP THREADPRIVATE(VelocityMass, PressureMass, ForcePart, Weight_a, weight_b, weight_c)
!$OMP THREADPRIVATE(tauVec, PrevTempVec, PrevPressureVec, VeloVec, PresVec, GradVec, Nodes)

    SAVE Nodes
!------------------------------------------------------------------------------

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

    ! Deallocate storage if needed 
    IF (ALLOCATED(BasisVec)) THEN
      IF (SIZE(BasisVec,1) < ngp .OR. SIZE(BasisVec,2) < ntot) &
          DEALLOCATE(BasisVec,dBasisdxVec, DetJVec, rhoVec, VeloVec, PresVec, &
          LoadAtIpVec, weight_a, weight_b, weight_c, tauVec, PrevTempVec, &
          PrevPressureVec, VeloPresVec, GradVec)
    END IF
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(BasisVec)) THEN
      ALLOCATE(BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
          rhoVec(ngp), VeloVec(ngp, dim), PresVec(ngp), velopresvec(ngp,dofs), LoadAtIpVec(ngp,dim+1), &
          weight_a(ngp), weight_b(ngp), weight_c(ngp), tauVec(ngp), PrevTempVec(ngp), &
          PrevPressureVec(ngp), GradVec(ngp,dim,dim), &
          STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('IncompressibleNSSolver::LocalBulkMatrix','Local storage allocation failed')
      END IF
    END IF


    ! Deallocate storage (ntot) if needed
    IF (ALLOCATED(VelocityMass)) THEN
      IF(SIZE(VelocityMass,1) < ntot ) THEN
        DEALLOCATE(VelocityMass, PressureMass, ForcePart)
      END IF
    END IF

    ! Allocate storage (ntot) if needed
    IF(.NOT. ALLOCATED(VelocityMass)) THEN
      ALLOCATE(VelocityMass(ntot,ntot), PressureMass(ntot, ntot), &
          ForcePart(ntot))
    END IF
           
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
              muDerVec0, Newton,  InitHandles, DetJVec )        

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

    IF (DivCurlForm) THEN
      weight_a(1:ngp) = muVec(1:ngp) * detJVec(1:ngp)
      weight_b(1:ngp) = -muVec(1:ngp) * detJVec(1:ngp) 

      ! The following assumes that the bulk viscosity of the fluid vanishes:
      weight_c(1:ngp) =  4.0_dp / 3.0_dp * muVec(1:ngp) * detJVec(1:ngp)

      ! curl-curl part and div-div parts
      SELECT CASE(dim)
      CASE(3) ! {{{
        i = 1; j = 1 
        ! curl-curl
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,2), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:,3), weight_a, stifford(:,:,i,j))
        ! div-div
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:, 1), weight_c, stifford(:,:,i,j))

        i=1;j=2
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,1), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:, 2), weight_c, stifford(:,:,i,j))

        i=1;j=3
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:,1), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:, 3), weight_c, stifford(:,:,i,j))

        i=2;j=1
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,2), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:, 1), weight_c, stifford(:,:,i,j))

        i=2;j=2
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,1), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:,3), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:, 2), weight_c, stifford(:,:,i,j))

        i=2;j=3
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:,2), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:, 3), weight_c, stifford(:,:,i,j))

        i=3;j=1
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,3), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:, 1), weight_c, stifford(:,:,i,j))

        i=3;j=2
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,3), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:, 2), weight_c, stifford(:,:,i,j))

        i=3;j=3
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,1), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,2), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 3), dbasisdxvec(:,:, 3), weight_c, stifford(:,:,i,j))
        ! }}}
      CASE(2)  ! {{{
        i = 1; j = 1
        ! curl-curl
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,2), weight_a, stifford(:,:,i,j))
        ! div-div
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:, 1), weight_c, stifford(:,:,i,j))

        i = 1; j = 2
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:,1), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:, 2), weight_c, stifford(:,:,i,j))

        i = 2; j = 1
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,2), weight_b, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:, 1), weight_c, stifford(:,:,i,j))

        i = 2; j = 2
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 1), dbasisdxvec(:,:,1), weight_a, stifford(:,:,i,j))
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            dbasisdxvec(:, :, 2), dbasisdxvec(:,:, 2), weight_c, stifford(:,:,i,j))
        ! }}}
      END SELECT

    ELSE !DivCurlForm

      weight_a(1:ngp) = muVec(1:ngp) * detJvec(1:ngp)      

      ! The following assumes that the bulk viscosity of the fluid vanishes:
!     weight_c(1:ngp) = -2.0_dp / 3.0_dp * muVec(1:ngp) * detJVec(1:ngp)

      DO i=1,dim
        DO j=1,dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(1:ngp,1:ntot,j), dBasisdxVec(1:ngp,1:ntot,j), weight_a, stifford(1:ntot,1:ntot,i,i))

          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(1:ngp,1:ntot,j), dBasisdxVec(1:ngp,1:ntot,i), weight_a, stifford(1:ntot,1:ntot,i,j))

!         CALL LinearForms_UdotV(ngp, ntot, elemdim, &
!             dBasisdxVec(1:ngp,1:ntot,i), dBasisdxVec(1:ngp,1:ntot,j), weight_c, stifford(1:ntot,1:ntot,i,j))
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

    ! Masses (use symmetry)
    ! Compute bilinear form G=G+(alpha u, u) = u .dot. (grad u) 
    IF ( .NOT. StokesFlow ) THEN
      CALL LinearForms_UdotU(ngp, ntot, elemdim, BasisVec, DetJVec, VelocityMass, rhovec)

      ! Scatter to the usual local mass matrix
      DO i = 1, dim
        mass(i::dofs, i::dofs) = mass(i::dofs, i::dofs) + VelocityMass(1:ntot, 1:ntot)
      END DO

      !mass(dofs::dofs, dofs::dofs) = mass(dofs::dofs, dofs::dofs) + PressureMass(1:ntot,1:ntot)

      ! These loop unrolls look bad, maybe do nicer weight precomputation?
      weight_a(1:ngp) = rhovec(1:ngp) * veloPresVec(1:ngp,1)
      weight_b(1:ngp) = rhovec(1:ngp) * veloPresVec(1:ngp,2)
      weight_c(1:ngp) = rhovec(1:ngp) * veloPresVec(1:ngp,3)
      DO i = 1, dim
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            basisvec, dbasisdxvec(:,:,1), detJvec, stifford(:,:,i,i), weight_a)
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            basisvec, dbasisdxvec(:,:,2), detJvec, stifford(:,:,i,i), weight_b)
        CALL LinearForms_UdotV(ngp, ntot, elemdim, &
            basisvec, dbasisdxvec(:,:,3), detJvec, stifford(:,:,i,i), weight_c)
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



!------------------------------------------------------------------------------

  CONTAINS


    FUNCTION EffectiveViscosityVec( ngp, BasisVec, dBasisdxVec, Element, NodalSol, &
        ViscDerVec, ViscNewton, InitHandles, DetJVec ) RESULT ( EffViscVec ) 

      INTEGER :: ngp
      REAL(KIND=dp) :: BasisVec(:,:), dBasisdxVec(:,:,:)
      TYPE(Element_t), POINTER :: Element
      REAL(KIND=dp) :: NodalSol(:,:)
      REAL(KIND=dp), ALLOCATABLE :: ViscDerVec(:)
      LOGICAL :: InitHandles , ViscNewton
      REAL(KIND=dp), POINTER  :: EffViscVec(:)
      REAL(KIND=dp), ALLOCATABLE :: DetJVec(:)
      
      LOGICAL :: Found     
      CHARACTER(LEN=MAX_NAME_LEN) :: ViscModel
      TYPE(ValueHandle_t), SAVE :: Visc_h, ViscModel_h, ViscExp_h, ViscCritical_h, &
          ViscNominal_h, ViscDiff_h, ViscTrans_h, ViscYasuda_h, ViscGlenExp_h, ViscGlenFactor_h, &
          ViscArrSet_h, ViscArr_h, ViscTLimit_h, ViscRate1_h, ViscRate2_h, ViscEne1_h, ViscEne2_h, &
          ViscTemp_h
      REAL(KIND=dp), SAVE :: R, NewtonRelax
      REAL(KIND=dp) :: c1, c2, c3, c4, Ehf, Tlimit, ArrheniusFactor, A1, A2, Q1, Q2, ViscCond
      LOGICAL, SAVE :: ConstantVisc = .FALSE., Visited = .FALSE., GotRelax = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, SAVE :: ss(:), s(:), ArrheniusFactorVec(:)
      REAL(KIND=dp), POINTER, SAVE :: ViscVec0(:), ViscVec(:), TempVec(:), EhfVec(:) 
      TYPE(Variable_t), POINTER, SAVE :: ShearVar, ViscVar, WeightVar
      LOGICAL, SAVE :: SaveShear, SaveVisc, SaveWeight
      CHARACTER(*), PARAMETER :: Caller = 'EffectiveViscosityVec'
     
!$OMP THREADPRIVATE(ss,s,ViscVec0,ViscVec,ArrheniusFactorVec)
     
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

        ShearVar => VariableGet( CurrentModel % Mesh % Variables,'Shearrate',ThisOnly=.TRUE.)
        SaveShear = ASSOCIATED(ShearVar)
        Found = .FALSE.
        IF(SaveShear) THEN
          IF(ShearVar % TYPE == Variable_on_gauss_points ) THEN
            CALL Info(Caller,'Saving "Shearrate" on ip points!',Level=10)
          ELSE IF( ShearVar % TYPE == Variable_on_elements ) THEN
            CALL Info(Caller,'Saving "Shearrate" on elements!',Level=10)
          ELSE IF( ShearVar % TYPE == Variable_on_nodes ) THEN
            CALL Info(Caller,'Saving "Shearrate" on nodes!',Level=10)
            Found = .TRUE.
          ELSE
            CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
          END IF
          ShearVar % Values = 0.0_dp
        END IF

        ViscVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity',ThisOnly=.TRUE.)
        SaveVisc = ASSOCIATED(ViscVar)        
        IF(SaveVisc) THEN
          IF(ViscVar % TYPE == Variable_on_gauss_points ) THEN
            CALL Info(Caller,'Saving "Viscosity" on ip points!',Level=10)
          ELSE IF( ViscVar % TYPE == Variable_on_elements ) THEN
            CALL Info(Caller,'Saving "Viscosity" on elements!',Level=10)
          ELSE IF( ViscVar % TYPE == Variable_on_nodes ) THEN
            CALL Info(Caller,'Saving "Viscosity" on nodes!',Level=10)
            Found = .TRUE.
          ELSE
            CALL Fatal(Caller,'Invalid field type for "Shearrate"!')
          END IF
          ViscVar % Values = 0.0_dp
        END IF

        SaveWeight = .FALSE.
        IF( Found ) THEN
          WeightVar => VariableGet( CurrentModel % Mesh % Variables,'Viscosity Weight',ThisOnly=.TRUE.)
          IF( ASSOCIATED( WeightVar ) ) THEN
            IF( WeightVar % TYPE /= Variable_on_nodes ) THEN
              CALL Fatal(Caller,'Invalid field type for "Viscosity Weight"!')
            END IF
            SaveWeight = .TRUE.
            WeightVar % Values = 0.0_dp
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
        
      ! Deallocate too small storage if needed 
      IF (ALLOCATED(ss)) THEN
        IF (SIZE(ss) < ngp ) DEALLOCATE(ss, s, ViscVec, ArrheniusFactorVec )
      END IF

      ! Allocate storage if needed
      IF (.NOT. ALLOCATED(ss)) THEN
        ALLOCATE(ss(ngp),s(ngp),ViscVec(ngp),ArrheniusFactorVec(ngp),STAT=allocstat)
        IF (allocstat /= 0) THEN
          CALL Fatal(Caller,'Local storage allocation failed')
        END IF
      END IF

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
        IF( ShearVar % TYPE == Variable_on_nodes ) THEN
          DO i=1,n
            ShearVar % Values(ShearVar % Perm(Element % NodeIndexes(i))) = &
                ShearVar % Values(ShearVar % Perm(Element % NodeIndexes(i))) + &
                SUM(BasisVec(1:ngp,i)*ss(1:ngp)*DetJVec(1:ngp))
          END DO
        ELSE
          i = Element % ElementIndex
          IF( ShearVar % TYPE == Variable_on_gauss_points ) THEN
            j = ShearVar % Perm(i+1) - ShearVar % Perm(i)
            IF(j /= ngp) THEN
              CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Shearrate" got '//I2S(ngp))
            END IF
            ShearVar % Values(ShearVar % Perm(i)+1:ShearVar % Perm(i+1)) = ss(1:ngp)
          ELSE IF( ShearVar % TYPE == Variable_on_elements ) THEN
            ShearVar % Values(ShearVar % Perm(i)) = SUM(ss(1:ngp)) / ngp
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
        IF( ViscVar % TYPE == Variable_on_nodes ) THEN
          DO i=1,n
            ViscVar % Values(ViscVar % Perm(Element % NodeIndexes(i))) = &
                ViscVar % Values(ViscVar % Perm(Element % NodeIndexes(i))) + &
                SUM(BasisVec(1:ngp,i)*ViscVec(1:ngp)*detJVec(1:ngp))
          END DO
        ELSE
          i = Element % ElementIndex
          IF( ViscVar % TYPE == Variable_on_gauss_points ) THEN
            j = ViscVar % Perm(i+1) - ViscVar % Perm(i) 
            IF(j /= ngp) THEN
              CALL Fatal(Caller,'Expected '//I2S(j)//' gauss point for "Viscosity" got '//I2S(ngp))
            END IF
            ViscVar % Values(ViscVar % Perm(i)+1:ViscVar % Perm(i+1)) = ViscVec(1:ngp)
          ELSE
            ViscVar % Values(ViscVar % Perm(i)) = SUM(ViscVec(1:ngp)) / ngp
          END IF
        END IF
      END IF

      ! If requested, save normalization weight associated to viscosity (and shearrate).
      IF(SaveWeight) THEN
        DO i=1,n
          WeightVar % Values(WeightVar % Perm(Element % NodeIndexes(i))) = &
              WeightVar % Values(WeightVar % Perm(Element % NodeIndexes(i))) + &
              SUM(BasisVec(1:ngp,i)*detJVec(1:ngp))
        END DO
      END IF
      
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
    LOGICAL, INTENT(INOUT) :: SpecificLoad, InitHandles, FrictionNewton 
!------------------------------------------------------------------------------    
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), TARGET :: STIFF(nd*(dim+1),nd*(dim+1)), FORCE(nd*(dim+1))
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str, FSSAFlag
    INTEGER :: c,i,j,k,l,p,q,t,ngp,norm_comp
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
    norm_comp = 0

    ! There is no elemental routine for this.
    ! So whereas this breaks the beuty it does not cost too much.
    FSSAFlag = GetString(BC, 'FSSA Flag', Found)
    IF (.NOT.Found) THEN
      WRITE (FSSAFlag,*) "none"
    END IF
    
    HaveFrictionW = ListCheckPresent( BC,'Weertman Friction Coefficient') 
    HaveFrictionU = ListCheckPresent( BC,'Friction Coefficient')
    HaveFriction = HaveFrictionU .OR. HaveFrictionW
    
    IF( HaveFriction ) THEN
      wut0 = ListGetElementReal( FrictionUt0_h, Element = Element )
      FrictionNormal = ListGetElementLogical( FrictionNormal_h, Element ) 
    END IF
    
    FSSAtheta = ListGetElementReal( FSSAtheta_h,  Basis, Element, HaveFSSA, GaussPoint = t )
    IF (HaveFSSA) THEN
      IF (FSSAtheta == 0.0) HaveFSSA=.FALSE.
      rho = ListGetElementRealParent( Dens_h, Basis, Element, Found )
      IF (.NOT.Found) THEN
        CALL WARN('IncompressibleNSSolver (FSSA)','"Density" in Parent element not found!')          
        HaveFSSA = .FALSE.
      END IF
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
          CALL WARN('IncompressibleNSSolver (FSSA)','No component of "Flow Body Force" in Parent element not found!')
          HaveFSSA = .FALSE.
        END IF
      END IF
      
      ! Nothing to do, exit the routine
      !---------------------------------
      IF(.NOT. (HaveForce .OR. HavePres .OR. HaveSlip .OR. HaveNormalSlip .OR. HaveFriction .OR. HaveFSSA)) RETURN

      ! Calculate normal vector only if needed
      IF( HavePres .OR. NormalTangential .OR. HaveFriction  .OR. HaveFSSA .OR. HaveNormalSlip ) THEN
        Normal = ConsistentNormalVector( CurrentModel % Solver, NrmSol, Element, Found, Basis = Basis )
        IF(.NOT. Found) Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t),.TRUE. )
      END IF
            
      !-----------------------------------------------------------------
      IF( HaveFriction  .AND. .NOT. HaveSlip ) THEN
        ! Velocity at integration point for nonlinear friction laws
        Velo = ListGetElementVectorSolution( Velo_v, Basis, Element, dofs = dim )

        ! It seems futile to take normal component away as it is usually zero by construction!
        ! We include here the option not to remove it. 
        IF(.NOT. FrictionNormal ) THEN
          un = SUM( Normal(1:dim) * Velo(1:dim) )
          velo(1:dim) = velo(1:dim)-un*normal(1:dim)

          ! Define normal component to skip for Weertman friction law!
          ! This is just the dominating direction that assumes rectangular directions.
          norm_comp = 1
          IF( .NOT. NormalTangential) THEN
            DO i=2,dim
              IF(ABS(Normal(i)) > ABS(Normal(norm_comp))) norm_comp = i
            END DO
          END IF
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

        DO i=1,dim
          IF(i==norm_comp) THEN
            SlipCoeff(i) = NormalSlipCoeff
          ELSE
            SlipCoeff(i) = TanFrictionCoeff
          END IF
        END DO
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

      IF(t==0 ) THEN
        PRINT *,'Normal:',Element % ElementIndex, Normal,un,ut,norm_comp
        PRINT *,'wexp:',wexp,wcoeff,wut0,HaveFrictionW,HaveFrictionU,HaveSlip,HaveNormalSlip
        PRINT *,'Velo:',dim,Velo,NormalTangential, LocalNewton
        PRINT *,'Slip:',SlipCoeff,NormalSlipCoeff,TanFrictionCoeff
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
            IF(i==norm_comp .AND. .NOT. HaveNormalSlip ) CYCLE
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
                IF(i == norm_comp .AND. .NOT. HaveNormalSlip ) CYCLE
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
              DO i=dim,dim
                STIFF( (p-1)*c+dim,(q-1)*c+i ) = & 
                    STIFF( (p-1)*c+dim,(q-1)*c+i )  &
                    - s * FSSAtheta * dt * LoadVec(dim) * Basis(q) * Basis(p) * Normal(i)
              END DO
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

  ! Create solver related to variable "schur" when using block preconditioning
  ! These keywords ensure that the matrix is truly used in the library version of the
  ! block solver.
  IF( ListGetLogical( Params,'Block Preconditioner',Found ) ) THEN
    CALL ListAddNewString( Params,'Block Matrix Schur Variable','schur')
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
  INTEGER :: i, n, nb, nd, nbdofs, dim, Active, maxiter, iter
  INTEGER :: stimestep = -1
 
  REAL(KIND=dp) :: Norm

  LOGICAL :: AllocationsDone = .FALSE., Found, StokesFlow, BlockPrec, Converged
  LOGICAL :: GradPVersion, DivCurlForm, SpecificLoad, InitBCHandles

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
  IP = GaussPointsAdapt( Element, PReferenceElement = .TRUE. )
  CALL Info('IncompressibleNSSolver', &
      'Number of 1st integration points: '//I2S(IP % n), Level=5)
  
  !-----------------------------------------------------------------------------
  ! Set the flags/parameters which define how the system is assembled: 
  !-----------------------------------------------------------------------------
  LinearAssembly = GetLogical(Params, 'Linear Equation', Found )
  StokesFlow = GetLogical(Params, 'Stokes Flow', Found )
  GradPVersion = GetLogical(Params, 'GradP Discretization', Found)
  DivCurlForm = GetLogical(Params, 'Div-Curl Discretization', Found)
  BlockPrec = GetLogical(Params,'Block Preconditioner',Found )
  SpecificLoad = GetLogical(Params,'Specific Load',Found)
  
  Maxiter = GetInteger(Params, 'Nonlinear system max iterations', Found)
  IF (.NOT.Found) Maxiter = 1
  !-----------------------------------------------------------------------------

  IF (DivCurlForm) CALL Info(Caller, 'The div-curl form is used for the viscous terms')
  IF (GradPVersion) CALL Info(Caller, 'The pressure gradient is not integrated by parts')
  IF (BlockPrec) CALL Info(Caller,'Creating pressure block for block preconditioner')

  IF(BlockPrec ) THEN
    ! Create solver that only acts as a container for the shcur complement
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
          SchurSolver )
    END DO
    
    !$OMP PARALLEL SHARED(Active, dim, SpecificLoad, StokesFlow, &
    !$OMP                 DivCurlForm, GradPVersion, &
    !$OMP                 dt, LinearAssembly, Newton, Transient, SchurSolver ) &
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
          SchurSolver )
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
