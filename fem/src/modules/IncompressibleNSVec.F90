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
! *  Authors: Mika Malinen, Juhani Kataja, Juha Ruokolainen, Peter Råback
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
  !$OMP THREADPRIVATE (bx, bxprev)

CONTAINS

  SUBROUTINE LocalMatrix(Element, n, nd, ntot, dim, ReadBodyForce, &
      ThermalCorrection, DivCurlForm, GradPVersion, &
      dt, RelOrder, LinearAssembly, nb, Newton, TransientSimulation)
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE

    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, ntot, dim, nb
    LOGICAL, INTENT(IN) :: ReadBodyForce
    LOGICAL, INTENT(IN) :: ThermalCorrection
    LOGICAL, INTENT(IN) :: DivCurlForm, GradPVersion
    REAL(KIND=dp), INTENT(IN) :: dt   
    INTEGER, INTENT(IN) :: RelOrder
    LOGICAL, INTENT(IN) :: LinearAssembly, Newton, TransientSimulation 
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material, BodyForce
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes

    LOGICAL :: Stat, GotLoad, Found

    REAL(KIND=dp), TARGET :: MASS(ntot*(dim+1),ntot*(dim+1)), &
        STIFF(ntot*(dim+1),ntot*(dim+1)), FORCE(ntot*(dim+1))
    REAL(KIND=dp), TARGET :: K(ntot*(dim+1),ntot*(dim+1))
    REAL(KIND=dp) :: s, LOAD(dim+1,n), LoadAtIP(dim+1), Velo(dim)    
    REAL(KIND=dp) :: NodalSol(dim+1,ntot), NodalTemp(nd)
    REAL(KIND=dp) :: PrevNodalSol(dim+1,ntot), PrevNodalTemp(nd)
    REAL(KIND=dp) :: mu, kappa, rho
    REAL(KIND=dp) :: Pressure, PrevPressure
    REAL(KIND=dp) :: Temp
    
    REAL(KIND=dp), ALLOCATABLE, SAVE :: BasisVec(:,:), dBasisdxVec(:,:,:), DetJVec(:), &
        rhoVec(:), kappaVec(:), VeloPresVec(:,:), loadAtIpVec(:,:), VelocityMass(:,:), &
        PressureMass(:,:), ForcePart(:), &
        weight_a(:), weight_b(:), weight_c(:), tauVec(:), PrevTempVec(:), PrevPressureVec(:), &
        VeloVec(:,:), PresVec(:), GradVec(:,:,:)

    REAL(kind=dp) :: stifford(ntot,ntot,dim+1,dim+1)

    INTEGER :: t, i, j, ii, jj, p, q, ngp, allocstat
    INTEGER, SAVE :: elemdim
    INTEGER :: DOFs

!DIR$ ATTRIBUTES ALIGN:64 :: BasisVec, dBasisdxVec, DetJVec, rhoVec, kappaVec, VeloPresVec, loadAtIpVec
!DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE, weight_a, weight_b, weight_c
!$OMP THREADPRIVATE(BasisVec, dBasisdxVec, DetJVec, rhoVec, kappaVec, VeloPresVec, loadAtIpVec, ElemDim )
!$OMP THREADPRIVATE(VelocityMass, PressureMass, ForcePart, Weight_a, weight_b, weight_c)
!$OMP THREADPRIVATE(tauVec, PrevTempVec, PrevPressureVec, VeloVec, PresVec, GradVec, Nodes)

    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodesVec( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0
    stifford = 0._dp

    DOFs = dim + 1

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    ngp = IP % n

    ElemDim = Element % Type % Dimension
    
    ! Storage size depending ngp
    !-------------------------------------------------------------------------------

    ! Deallocate storage if needed 
    IF (ALLOCATED(BasisVec)) THEN
       IF (SIZE(BasisVec,1) < ngp .OR. SIZE(BasisVec,2) < ntot) &
            DEALLOCATE(BasisVec,dBasisdxVec, DetJVec, rhoVec, VeloVec, PresVec, kappaVec, &
            LoadAtIpVec, weight_a, weight_b, weight_c, tauVec, PrevTempVec, &
           PrevPressureVec, VeloPresVec, GradVec)
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(BasisVec)) THEN
       ALLOCATE(BasisVec(ngp,ntot), dBasisdxVec(ngp,ntot,3), DetJVec(ngp), &
           rhoVec(ngp), VeloVec(ngp, dim), PresVec(ngp), velopresvec(ngp,dofs), kappaVec(ngp), LoadAtIpVec(ngp,dim+1), &
           weight_a(ngp), weight_b(ngp), weight_c(ngp), tauVec(ngp), PrevTempVec(ngp), &
           PrevPressureVec(ngp), GradVec(ngp,dim,dim), &
           STAT=allocstat)
       IF (allocstat /= 0) THEN
          CALL Fatal('IncompressibleNSSolver::LocalMatrix','Local storage allocation failed')
       END IF
    END IF

    ! Storage size depending ntot
    !-------------------------------------------------------------------------------

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

    ! Volume forces:
    !---------------
    GotLoad = .FALSE.

    IF ( ReadBodyForce ) THEN
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) THEN
        LOAD = 0._dp

        Load(1, 1:n) = GetReal( BodyForce, 'Flow Source 1', Found )
        GotLoad = GotLoad .OR. Found
        
        Load(2, 1:n) = GetReal( BodyForce, 'Flow Source 2', Found )
        GotLoad = GotLoad .OR. Found
        
        Load(3, 1:n) = GetReal( BodyForce, 'Flow Source 3', Found )
        GotLoad = GotLoad .OR. Found
      END IF
    END IF
      
    ! Material parameters:
    !---------------------
    Material => GetMaterial()
    mu = ListGetCReal(Material, 'Viscosity')
    rho = ListGetCReal(Material, 'Density')
    
    ! Get the previous elementwise velocity-pressure iterate:
    !--------------------------------------------------------
    IF ( LinearAssembly ) THEN
      VeloPresVec = 0._dp
    ELSE
      CALL GetLocalSolution( NodalSol )
      IF (nb > 0 .AND. TransientSimulation .OR. ThermalCorrection) &
          CALL GetLocalSolution(PrevNodalSol, tStep=-1)
    END IF

    Temp = 0._dp 

    VelocityMass = 0.0d0
    !PressureMass = 0.0d0

    stat = ElementInfoVec(Element, nodes, ngp, IP % U, IP % V, &
        IP % W, detJvec, SIZE(basisVec, 2), BasisVec, dBasisdxVec)
    
    DO t = 1, ngp
      DetJVec(t) = DetJVec(t) * IP % s(t)
    END DO
    
    IF(.NOT. LinearAssembly) THEN
      CALL DGEMM('N', 'T', ngp, DOFs, n, 1._dp, &
          BasisVec, ngp, NodalSol, dofs, 0._dp, &
          VeloPresVec, ngp)
    END IF
    
    ! Rho and Kappa
    rhovec(1:ngp) = rho
    !kappavec(1:ngp) = 0.0d0

    ! Load
    LoadAtIpVec = 0._dp
    IF (GotLoad) THEN
      CALL DGEMM('N', 'T', ngp, DOFs, n, 1._dp, &
          BasisVec, ngp, load, dofs, 0._dp, &
          LoadAtIpVec, ngp)
      ! The preceding DGEMM call is equivalent to
      ! LoadAtIpVec(1:ngp, 1:dofs) = MATMUL(basisvec(1:ngp, 1:n), transpose(load(1:dofs, 1:n)))
    END IF

    IF ( Newton ) THEN
      do i = 1,dim
        do j = 1,dim
          GradVec(1:ngp, i, j) = MATMUL(dBasisdxVec(1:ngp,1:ntot,j),nodalsol(i,1:ntot))
        end do
        LoadAtIpVec(1:ngp, i) = LoadAtIpVec(1:ngp, i) + rhovec(1:ngp)*sum(gradvec(1:ngp,i,1:dim)*velopresvec(1:ngp,1:dim),2)
      end do
    END IF

    
    IF (DivCurlForm) THEN
      weight_a = mu * detJVec
      weight_b = -mu * detJVec

      ! The following assumes that the bulk viscosity of the fluid vanishes:
      weight_c =  4.0_dp / 3.0_dp * mu * detJVec

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
        i=2;j=1
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
      
      weight_a = mu*detJvec      

      ! The following assumes that the bulk viscosity of the fluid vanishes:
      weight_c = -2.0_dp / 3.0_dp * mu * detJVec

      DO i=1,dim
        DO j=1,dim
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(1:ngp,1:ntot,j), dBasisdxVec(1:ngp,1:ntot,j), weight_a, stifford(1:ntot,:ntot,i,i))
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(1:ngp,1:ntot,j), dBasisdxVec(1:ngp,1:ntot,i), weight_a, stifford(1:ntot,1:ntot,i,j))
          CALL LinearForms_UdotV(ngp, ntot, elemdim, &
              dBasisdxVec(1:ngp,1:ntot,i), dBasisdxVec(1:ngp,1:ntot,j), weight_c, stifford(1:ntot,:ntot,i,j))
        END DO
      END DO

  END IF

  ! Masses (use symmetry)
  ! Compute bilinear form G=G+(alpha u, u) = u .dot. (grad u) 
  CALL LinearForms_UdotU(ngp, ntot, elemdim, BasisVec, DetJVec, VelocityMass, rhovec)

  ! Scatter to the usual local mass matrix
  DO i = 1, dim
    mass(i::dofs, i::dofs) = mass(i::dofs, i::dofs) + VelocityMass(1:ntot, 1:ntot)
  END DO
  !CALL LinearForms_UdotU(ngp, ntot, elemdim, BasisVec, DetJVec, PressureMass, -kappavec)

  !mass(dofs::dofs, dofs::dofs) = mass(dofs::dofs, dofs::dofs) + PressureMass(1:ntot,1:ntot)


  IF (GradPVersion) THEN
    ! b(u,q) = (u, grad q) part
    DO i = 1, dim
      CALL LinearForms_UdotV(ngp, ntot, elemdim, &
          BasisVec, dbasisdxvec(:,:,i), detJVec, stifford(:,:,i,dofs))
      stifford(:,:,dofs,i) = transpose(stifford(:,:,i,dofs))
    END DO
  ELSE
    DO i = 1, dim
      CALL LinearForms_UdotV(ngp, ntot, elemdim, &
          dBasisdxVec(:, :, i), BasisVec, -detJVec, StiffOrd(:,:,i,dofs))
      StiffOrd(:,:,dofs,i) = transpose(stifford(:,:,i,dofs))
    END DO
  END IF

  ! These loop unrolls look bad, maybe do nicer weight precomputation?
  weight_a = rhovec*veloPresVec(1:ngp,1)
  weight_b = rhovec*veloPresVec(1:ngp,2)
  weight_c = rhovec*veloPresVec(1:ngp,3)
  do i = 1, dim
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        basisvec, dbasisdxvec(:,:,1), detJvec, stifford(:,:,i,i), weight_a)
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        basisvec, dbasisdxvec(:,:,2), detJvec, stifford(:,:,i,i), weight_b)
    CALL LinearForms_UdotV(ngp, ntot, elemdim, &
        basisvec, dbasisdxvec(:,:,3), detJvec, stifford(:,:,i,i), weight_c)
  end do

  do i = 1, dim
    do j = 1, dim
      IF ( Newton ) THEN
        CALL LinearForms_UdotV(ngp, ntot, Element%TYPE%DIMENSION, &
            basisvec, basisvec, detJvec, stifford(:,:,i,j), rhovec*gradvec(:,i,j))
      end if
    end do
  end do

  !weight_a = -kappavec*veloPresVec(1:ngp,1)
  !weight_b = -kappavec*veloPresVec(1:ngp,2)
  !weight_c = -kappavec*veloPresVec(1:ngp,3)
  !CALL LinearForms_UdotV(ngp, ntot, elemdim, &
  !    basisvec, dbasisdxvec(:,:,1), detJvec, stifford(1:ntot,1:ntot,dofs,dofs), weight_a)
  !CALL LinearForms_UdotV(ngp, ntot, elemdim, &
  !    basisvec, dbasisdxvec(:,:,2), detJvec, stifford(1:ntot,1:ntot,dofs,dofs), weight_b)
  !CALL LinearForms_UdotV(ngp, ntot, elemdim, &
  !    basisvec, dbasisdxvec(:,:,3), detJvec, stifford(1:ntot,1:ntot,dofs,dofs), weight_c)

  ! convection and loads
  DO ii = 1,dim+1
    ForcePart = 0._dp
    CALL LinearForms_UdotF(ngp, ntot, basisVec, detJVec, LoadAtIpVec(:,ii), ForcePart)
    FORCE(ii::dofs) = ForcePart(1:ntot)
  END DO

  DO i = 1, DOFS
    DO j = 1, DOFS
      Stiff(i::DOFS, j::DOFS) = StiffOrd(1:ntot, 1:ntot, i,j)
    END DO
  END DO


  IF (nb > 0 .AND. nd==n .AND. TransientSimulation) THEN
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
      s = STIFF(i,i)
      STIFF(i,:) = 0._dp
      STIFF(:,i) = 0._dp
      STIFF(i,i) = ABS(s)
    END DO

    IF ( TransientSimulation ) THEN
      CALL Default1stOrderTime( MASS, STIFF, FORCE )
    END IF
    IF (nb > 0) THEN
      IF (TransientSimulation) THEN
        CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE, &
            PrevNodalSol, NodalSol, Element % ElementIndex)
      ELSE
        CALL LCondensate(nd, nb, dim, MASS, STIFF, FORCE)
      END IF
    END IF
  END IF

  CALL DefaultUpdateEquationsR( STIFF, FORCE, UElement=Element, VecAssembly = .TRUE.)
  !------------------------------------------------------------------------------

  CONTAINS

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
      xlprev = 0
      q = 0
      DO p = 1,n
        DO i = 1,DOFs
          q = q + 1
          Cdofs(q) = DOFs*(p-1) + i
          IF (ComputeBubblePart) THEN
            xl(q) = x(i,p)
            xlprev(cdofs(q)) = xprev(i,p) ! cdofs identity mapping?
          END IF
        END DO
      END DO
        
      ! Then the DOFs of the bubble part:
      q = 0
      DO p = 1,nb
        DO i = 1,dim
          q = q + 1
          Bdofs(q) = DOFs*(p-1) + i + n*DOFs
          IF (ComputeBubblePart) xlprev(bdofs(q)) = &
              bxprev((Element_id-1)*dim*nb+q)
        END DO
      END DO

      ! The following only works for the BDF(1) method: 
      IF (ComputeBubblePart) THEN
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

END SUBROUTINE LocalMatrix

END MODULE IncompressibleLocalForms


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver_init(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params 
  LOGICAL :: Found
  INTEGER :: dim
!------------------------------------------------------------------------------ 
  Params => GetSolverParams() 

  dim = CoordinateSystemDimension()

  IF ( dim == 2 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow Solution[Velocity:2 Pressure:1]')
  ELSE
    CALL ListAddNewString(Params, 'Variable', &
        'Flow Solution[Velocity:3 Pressure:1]')
  END IF

  ! Study only velocity components in linear system
  CALL ListAddNewInteger( Solver % Values, 'Nonlinear System Norm DOFs', dim )
  
  ! Automate the choice for the variational formulation:
  CALL ListAddNewLogical(Params, 'GradP Discretization', .FALSE.)
  CALL ListAddNewLogical(Params, 'Div-Curl Form', .TRUE.)
  
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
!------------------------------------------------------------------------------ 
END SUBROUTINE IncompressibleNSSolver_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  USE IncompressibleLocalForms
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: TempVar => NULL()
  TYPE(ValueList_t), POINTER :: Material
  TYPE(GaussIntegrationPoints_t) :: IP

  INTEGER :: Element_id
  INTEGER :: i, n, nb, nd, nbdofs, dim, Active, maxiter, iter, RelOrder
  INTEGER :: stimestep = -1

  REAL(KIND=dp) :: Norm, pres

  LOGICAL :: AllocationsDone = .FALSE., Found
  LOGICAL :: ReadBodyForce
  LOGICAL :: ThermalCorrection, GradPVersion, DivCurlForm

  CHARACTER(*), PARAMETER :: Caller = 'IncompressibleNSSolver'

  
  SAVE AllocationsDone, stimestep
  !$OMP THREADPRIVATE(AllocationsDone, stimestep)

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
  !$OMP PARALLEL DEFAULT(shared) PRIVATE(nbdofs) 
  IF (.NOT. AllocationsDone .AND. TransientSimulation .AND. &
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
  IF (TransientSimulation .AND. GetTimestep() /= stimestep .AND. &
      Mesh % MaxBDOFs > 0) THEN
    bxprev = bx
    stimestep = GetTimestep()
  END IF
  !$OMP END PARALLEL

  Params => GetSolverParams() 


  !-----------------------------------------------------------------------------
  ! Output the number of integration points as information.
  ! This in not fully informative if several element types are present.
  !-----------------------------------------------------------------------------
  Element => Mesh % Elements( Solver % ActiveElements(1) ) 
  IP = GaussPointsAdapt( Element )
  CALL Info('IncompressibleNSSolver', &
      'Number of 1st integration points: '//TRIM(I2S(IP % n)), Level=5)
  !-----------------------------------------------------------------------------

  Material => GetMaterial( Element )

  !TempVar => VariableGet(Mesh % Variables, 'Temperature')
  !ThermalCorrection = ASSOCIATED( TempVar )
  !
  ! Disable the temperature coupling:
  ThermalCorrection = .FALSE.

  !-----------------------------------------------------------------------------
  ! Set the flags/parameters which define how the system is assembled: 
  !-----------------------------------------------------------------------------
  LinearAssembly = GetLogical(Params, 'Linear Equation', Found )
  GradPVersion = GetLogical(Params, 'GradP Discretization', Found)
  DivCurlForm = GetLogical(Params, 'Div-Curl Form', Found)
  ReadBodyForce = ListCheckPrefixAnyBodyForce(Model, 'Body Force')

  Maxiter = GetInteger(Params, 'Nonlinear system max iterations', Found)
  IF (.NOT.Found) Maxiter = 1
  !-----------------------------------------------------------------------------

  IF (DivCurlForm) CALL Info(Caller, 'The div-curl form is used for the viscous terms')
  IF (GradPVersion) CALL Info(Caller, 'The pressure gradient is not integrated by parts')


  DO iter=1,maxiter

    CALL Info(Caller,'--------------------------------------------------------', Level=4)
    WRITE( Message,'(A,I4)') 'Nonlinear iteration:', Iter
    CALL Info(Caller, Message, Level=4)
    CALL Info(Caller,'--------------------------------------------------------', Level=4)

    Active = GetNOFActive()
    CALL DefaultInitialize()
    call ResetTimer('IncompressibleNSBulkAssembly')

    Newton = GetNewtonActive( Solver )
    
    !$OMP PARALLEL SHARED(ReadBodyForce, ThermalCorrection, Active, dim, &
    !$OMP                 DivCurlForm, GradPVersion, &
    !$OMP                 dt, RelOrder, LinearAssembly, Newton, TransientSimulation) &
    !$OMP          PRIVATE(Element, Element_id, n, nd, nb) &
    !$OMP          DEFAULT(None)
    !$OMP DO 
    DO Element_id=1,Active
      Element => GetActiveElement(Element_id)
      n  = GetElementNOFNodes(Element)
      !
      ! When the number of bubbles is obtained with the Update=.TRUE. flag,
      ! we need to call GetElementNOFBDOFs before calling GetElementNOFDOFs.
      !
      nb = GetElementNOFBDOFs(Element, Update=.TRUE.)
      nd = GetElementNOFDOFs(Element)

      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalMatrix(Element, n, nd, nd+nb, dim, ReadBodyForce, &
          ThermalCorrection, DivCurlForm, GradPVersion, &
          dt, RelOrder, LinearAssembly, nb, Newton, TransientSimulation)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL DefaultFinishBulkAssembly()
    call CheckTimer('IncompressibleNSBulkAssembly', level=5)


    Active = GetNOFBoundaryElements()
    DO Element_id=1,Active
      Element => GetBoundaryElement(Element_id)
      IF (ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        !
        ! Elements defining the boundary can never have DOFs that could be
        ! eliminated via the static consensation, so here the function call
        ! nb = GetElementNOFBDOFs() would be worthless
        !

        ! A ROUTINE FOR THE LOCAL MATRIX ASSEMBLY SHOULD BE CALLED HERE 

      END IF
    END DO

    !CALL DefaultFinishBoundaryAssembly()

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    Norm = DefaultSolve()

    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO
 
  CALL DefaultFinish()

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE IncompressibleNSSolver
!------------------------------------------------------------------------------
