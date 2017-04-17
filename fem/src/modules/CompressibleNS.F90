!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
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
! *  A monolithic compressible Navier-Stokes solver 
! *
! ******************************************************************************
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           P.O. Box 405
! *           FI-02101 Espoo, Finland 
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> A monolithic compressible Navier-Stokes solver. The equation includes
!> velocity components, pressure and temperature. Ideal gas law is assumed. 
!------------------------------------------------------------------------------
SUBROUTINE CompressibleNS( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation    !< Steady state or transient simulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: AMatrix, SMatrix, A1Matrix, A2Matrix, &
       A3Matrix, SMatrix2, MMatrix, HMatrix

  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found, Convect, &
       OptimizeBW, GotIt, ExplicitStabilization, &
       ComponentwisePreconditioning =.TRUE., BlockPreconditioning = .FALSE., &
       NormalTractionBoundary, SlipBoundary

  TYPE(Element_t),POINTER :: Element, Parent
  INTEGER, POINTER :: NodeIndexes(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: OuterIterationMethod, NonlinearIterationMethod
  INTEGER :: i,j,k,n, nb, nd, t, istat, dim, m, p, q, &
       MaxIterations, NumberOfNormalTractionNodes, BDOFs=0, NormalDir
  REAL(KIND=dp) :: Norm = 0, PrevNorm, Tolerance, ToleranceRatio, &
       atime, stime, at0, SlipCoefficient, &
       gamma, cv, kcoeff, rho0, T0, lambda, bulkvisc
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime, RealTime
#endif

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), Mass(:,:), &
       FORCE(:), rho(:), mu(:), Velocity(:,:), Temperature(:), &
       ALocal(:,:), SLocal(:,:), AnLocal(:,:), BoundaryVelocities(:), &
       MLocal(:,:), HLocal(:,:), Pressure(:)

  REAL(KIND=dp), ALLOCATABLE :: ForceVector(:)
  INTEGER, ALLOCATABLE :: TractionBCIndeces(:)


  LOGICAL :: NewTimeStep, OutflowBC, PicardIteration, FirstSolve
  INTEGER :: CurrentDoneTime = 0, NonlinearIter, iter, CoordSys
  REAL(KIND=dp) :: NonlinearTol

  INTEGER, ALLOCATABLE :: Indexes(:)

  SAVE STIFF, LOAD, FORCE, rho, Temperature, mu, Velocity, AllocationsDone, &
       AMatrix, SMatrix, SMatrix2, ALocal, SLocal, &
       A1Matrix, A2Matrix, A3Matrix, AnLocal, Mass, TractionBCIndeces, &
       BoundaryVelocities, CurrentDoneTime, Indexes, &
       MMatrix, HMatrix, MLocal, HLocal, Pressure, Norm
  !------------------------------------------------------------------------------

  IF ( CurrentDoneTime == 0 ) THEN
     FirstSolve = .TRUE.
  ELSE 
     FirstSolve = .FALSE.
  END IF

  NewTimeStep = .FALSE.
  IF (CurrentDoneTime /= Solver % DoneTime) THEN
     NewTimeStep = .TRUE.
     ! IF (NewTimeStep) PRINT *, 'New time step begins...'
     CurrentDoneTime = CurrentDoneTime + 1
  END IF

  dim = CoordinateSystemDimension()
  CoordSys = CurrentCoordinateSystem()  
  ! Check whether convection stabilization is used...
  !ExplicitStabilization = ListGetLogical( Solver % Values, 'Stabilize', GotIt )
  !IF (.NOT. GotIt) ExplicitStabilization = .false.
  ExplicitStabilization = .FALSE.
  NonlinearIterationMethod = ListGetString(Solver % Values, 'Nonlinear Iteration Method', GotIt)
  IF ( .NOT. GotIt )  NonlinearIterationMethod = 'picard'
  IF ( NonlinearIterationMethod /= 'picard') THEN
     PicardIteration = .FALSE.
  ELSE
     PicardIteration = .TRUE.     
  END IF

  !--------------------------------------------------------------------------
  ! Check whether the block preconditioning is to be used.
  !---------------------------------------------------------------------------
  BlockPreconditioning = ListGetLogical( Solver % Values, 'Block Preconditioning', GotIt )
  BlockPreconditioning = .FALSE.

  !----------------------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !----------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     p = Solver % Mesh % MaxElementDOFs       ! The size of the one-component stiffness matrix block
     n = (dim+2)*p                            ! The size of the complete stiffness matrix
     m = dim * Solver % Mesh % MaxElementDOFs ! The size of the (1,1)-block (the velocity part)

     ALLOCATE( FORCE(n), LOAD(p,5), STIFF(n,n), Mass(n,n),&
          rho(p), Temperature(p), mu(p), Velocity(dim+1,p), Pressure(p), &
          ALocal(m+1,m+1), SLocal(2*p,2*p), Indexes(p), &
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'NavierStokesSolver', 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
  END IF

  Velocity = 0.0d0
  !-------------------------------------------------------------------------------
  ! Initialization and assembly for the primary system and preconditioning systems
  !-------------------------------------------------------------------------------
  !CALL DefaultInitialize()
  !IF (TransientSimulation) CALL InitializeTimestep(Solver)

  atime = CPUTime()
  at0 = RealTime()

  NonlinearIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Max Iterations', minv=0 )
  NonlinearTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance',minv=0.0d0 )

  
  CALL DefaultStart()

  DO iter=1, NonlinearIter
     ! Initialize the system matrices and vectrors...
     CALL DefaultInitialize()

     DO t=1,Solver % NumberOfActiveElements
        IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                (Solver % NumberOfActiveElements - t) / &
                (Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'NavierStokesSolver', Message, Level=5 )
           at0 = RealTime()
        END IF
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()

        nd = GetElementDOFs( Indexes )
        nb = GetElementNOFBDOFs()

        !-----------------------------------------------
        ! Volume forces:
        !-----------------------------------------------
        BodyForce => GetBodyForce()
        LOAD = 0.0d0
        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1:n,1) = GetReal( BodyForce, 'Body Force 1', Found )
           Load(1:n,2) = GetReal( BodyForce, 'Body Force 2', Found )
           Load(1:n,3) = GetReal( BodyForce, 'Body Force 3', Found )
        END IF
        !-----------------------------------------------
        ! Material parameters:
        !-----------------------------------------------
        Material => GetMaterial()
        !rho(1:n) = GetReal( Material, 'Density' )
        mu(1:n)  = GetReal( Material, 'Viscosity' )
        cv = ListGetConstReal( Material, 'Specific Heat')
        kcoeff = ListGetConstReal( Material, 'Heat Conductivity')
        gamma = ListGetConstReal( Material, 'Specific Heat Ratio')
        bulkvisc = ListGetConstReal( Material, 'Bulk Viscosity', Found)
        IF (.NOT. Found) bulkvisc = 0.0d0
        lambda = bulkvisc - 2.0d0/3.0d0 * mu(1)
        T0 = ListGetConstReal( Material, 'Equilibrium Temperature')
        rho0 = ListGetConstReal( Material, 'Equilibrium Density')          

        !---------------------------------------------------------------------
        ! Get previous elementwise velocity, temperature and density iterates:
        !---------------------------------------------------------------------
        IF (FirstSolve .AND. iter == 1) THEN
           Temperature(1:nd) = 0.0d0
           rho(1:nd) = 0.0d0
           DO i=1,dim
              Velocity(i,1:nd) = 0.0d0
           END DO
        ELSE
           IF ( NewTimeStep .AND. iter==1 ) THEN
              !print *, 'Using previous velo'
              DO i=1,dim
                 Velocity(i,1:nd) = Solver % Variable % PrevValues( &
                      Solver % Variable % DOFs*(Solver % Variable % &
                      Perm(Indexes(1:nd))-1)+i,1) 
              END DO
              Temperature(1:nd) = Solver % Variable % PrevValues( &
                   Solver % Variable % DOFs*(Solver % Variable % &
                   Perm(Indexes(1:nd))-1)+dim+1,1)
              rho(1:nd) = Solver % Variable % PrevValues( &
                   Solver % Variable % DOFs*(Solver % Variable % &
                   Perm(Indexes(1:nd))-1)+dim+2,1)
           ELSE
              !print *, 'Using current velo'
              DO i=1,dim
                 Velocity(i,1:nd) = Solver % Variable % Values( &
                      Solver % Variable % DOFs*(Solver % Variable % &
                      Perm(Indexes(1:nd))-1)+i)
              END DO
              Temperature(1:nd) = Solver % Variable % Values( &
                   Solver % Variable % DOFs*(Solver % Variable % &
                   Perm(Indexes(1:nd))-1)+dim+1)
              rho(1:nd) = Solver % Variable % Values( &
                   Solver % Variable % DOFs*(Solver % Variable % &
                   Perm(Indexes(1:nd))-1)+dim+2)
           END IF
        END IF

        !------------------------------------------------------
        ! Get element local matrix and rhs vector:
        !-------------------------------------------------------
        CALL LocalMatrix(  STIFF, Mass, FORCE, LOAD, Temperature, rho, mu, &
             Velocity, Element, n, nd, dim, ExplicitStabilization, &
             PicardIteration, gamma, cv, kcoeff, lambda, T0, rho0)

        IF (TransientSimulation) THEN
           CALL Default1stOrderTime( Mass, STIFF, FORCE )
        END IF
        !-----------------------------------------------------------------
        ! Update global matrix and rhs vector from local matrix & vector:
        !----------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )

     END DO

     CALL DefaultFinishBulkAssembly()

     ! No flux BCs
     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     atime = CPUTime() - atime

     !---------------------------------------------------------------------------
     ! The solution of the linear system...
     !---------------------------------------------------------------------------
     stime = CPUTime()
     PrevNorm = Norm

     Norm = DefaultSolve()

     stime = CPUTime() - stime
     WRITE(Message, '(a,F8.2)') ' Assembly:  (s)', atime
     CALL Info( 'NavierStokesSolver', Message, Level=4)
     WRITE(Message, '(a,F8.2)') ' Solution:  (s)', stime    
     CALL Info( 'NavierStokesSolver', Message, Level=4)

     IF ( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()
  
  
CONTAINS



  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, Mass, FORCE, LOAD, NodalT, Nodalrho, &
       Nodalmu, NodalVelo, Element, n, nd, dim, Stabilization, &
       PicardIteration, gamma, cv, k, la, T0, rho0)
    !------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), Mass(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), NodalT(:), Nodalrho(:), NodalVelo(:,:), &
         gamma, cv, k, la, T0, rho0
    INTEGER :: dim, n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Stabilization, PicardIteration
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), &
         DetJ,LoadAtIP(dim+2),Velo(dim), Grad(dim,dim), AK, w3
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c, c2, ch, rotterm, PrevT, R, &
         GradV1(2), GradV2(2), rpos, TestTerm

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    PicardIteration = .TRUE.
    Stabilization = .FALSE.
 
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    Mass = 0.0d0

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )  

    !AK = 0.0d0
    !ch = (sqrt( (Nodes % x(2) - Nodes % x(1))**2 + (Nodes % y(2) - Nodes % y(1))**2) + & 
    !     sqrt( (Nodes % x(2) - Nodes % x(3))**2 + (Nodes % y(2) - Nodes % y(3))**2) + &
    !     sqrt( (Nodes % x(1) - Nodes % x(3))**2 + (Nodes % y(1) - Nodes % y(3))**2) )/ 3.0d0 
    !rotterm = 0.0d0

    R = (gamma - 1.0d0) * cv

    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ
       IF (CoordSys == AxisSymmetric) THEN
          rpos = SUM( Basis(1:n) * Nodes % x(1:n) )
          s = s * rpos
       END IF

       !AK = AK + s

       Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
       GradV1(1) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,1) )
       GradV1(2) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )
       GradV2(1) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) )
       GradV2(2) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,2) )
        
       !w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
       !     SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )     

       !----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:nd) * Nodalrho(1:nd) ) + rho0
       PrevT = SUM( Basis(1:nd) * NodalT(1:nd) ) + T0
       !--------------------------------------------
       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:dim+2) = MATMUL( Basis(1:n), LOAD(1:n,1:dim+2) )

       !rotterm = rotterm + s * ( w3*w3 * (Velo(1) * Velo(1) + Velo(2) * Velo(2)) )

       !----------------------------------------------------------------------------------
       ! The system matrix with only the velocity space augmented by bubbles  
       !---------------------------------------------------------------------------------
       DO p=1,nd
          DO q=1,nd
             i = (dim+2) * (p-1) + 1
             j = (dim+2) * (q-1) + 1
             A => STIFF(i:i+dim+1,j:j+dim+1)
             M => Mass(i:i+dim+1,j:j+dim+1)
             DO i=1,dim

                IF ( (CoordSys == AxisSymmetric) .AND. i==1 ) THEN
                   TestTerm = dBasisdx(p,i) + Basis(p)/rpos
                ELSE
                   TestTerm = dBasisdx(p,i)
                END IF

                DO j = 1,dim
                   A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                   A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
                   ! A(i,j) = A(i,j) + s * la * dBasisdx(q,j) * dBasisdx(p,i)
                   A(i,j) = A(i,j) + s * la * dBasisdx(q,j) * TestTerm

                   IF ( (CoordSys == AxisSymmetric) .AND. j==1 ) THEN
                      A(i,j) = A(i,j) + s * la * 1/rpos * Basis(q) * TestTerm
                   END IF

                END DO
                
                IF ( (CoordSys == AxisSymmetric) .AND. i==1 ) THEN
                   A(i,i) = A(i,i) + s * 2.0d0 * mu * 1/rpos**2 * Basis(q) * Basis(p)
                END IF
                   

                M(i,i) = M(i,i) + s * rho * Basis(p) * Basis(q)
                
                ! Temperature and density bubbles are handled elsewhere... 
                IF (q <= n) THEN
                   ! Testing w.r.t. all velocity test functions and omitting temperature/density bubbles...
                   A(i,dim+1) = A(i,dim+1) - s * R * rho * Basis(q) * dBasisdx(p,i)
                   A(i,dim+2) = A(i,dim+2) - s * R * PrevT * Basis(q) * dBasisdx(p,i)

                   IF ( (CoordSys == AxisSymmetric) .AND. i==1 ) THEN
                      A(i,dim+1) = A(i,dim+1) - s * R * rho * Basis(q) * Basis(p)/rpos
                      A(i,dim+2) = A(i,dim+2) - s * R * PrevT * Basis(q) * Basis(p)/rpos
                   END IF

                END IF
                IF (p <= n) THEN
                   ! Testing w.r.t  standard temperature/density test functions...
                   A(dim+1,i) = A(dim+1,i) + s * R * rho * dBasisdx(q,i) * Basis(p)
                   A(dim+2,i) = A(dim+2,i) + s * rho * dBasisdx(q,i) * Basis(p)

                   IF ( (CoordSys == AxisSymmetric) .AND. i==1 ) THEN
                      A(dim+1,i) = A(dim+1,i) + s * R * rho * Basis(q)/rpos * Basis(p)
                      A(dim+2,i) = A(dim+2,i) + s * rho * Basis(q)/rpos * Basis(p)
                   END IF

                END IF
             END DO

             M(dim+1,dim+1) = M(dim+1,dim+1) + s * rho * cv / PrevT * Basis(p) * Basis(q)  
             M(dim+2,dim+2) = M(dim+2,dim+2) + s * Basis(p) * Basis(q) 

             A(dim+1,dim+1) = A(dim+1,dim+1) + s * k / PrevT * SUM( dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )
             A(dim+1,dim+1) = A(dim+1,dim+1) + s * rho * cv / PrevT * &
                  ( Velo(1) * dBasisdx(q,1) + Velo(2) * dBasisdx(q,2) ) * Basis(p)
             A(dim+2,dim+2) = A(dim+2,dim+2) + s * ( Velo(1) * dBasisdx(q,1) + Velo(2) * dBasisdx(q,2) ) * &
                  Basis(p)

             IF (dim==3) THEN
                A(dim+1,dim+1) = A(dim+1,dim+1) + s * rho * cv / PrevT * &
                      Velo(3) * dBasisdx(q,3) * Basis(p)
                A(dim+2,dim+2) = A(dim+2,dim+2) + s * Velo(3) * dBasisdx(q,3) * &
                  Basis(p)                
             END IF

 
             IF (.FALSE.) THEN
!             IF ( PicardIteration ) THEN 
                Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) = Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) + &
                     rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+2 ) = Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+2 ) - &
                     rho * s * Velo(2) * dBasisdx(q,1) * Basis(p) 

                Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+1 ) = Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+1 ) - &
                     rho * s * Velo(1) * dBasisdx(q,2) * Basis(p) 
                Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) = Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) + &
                     rho * s * Velo(1) * dBasisdx(q,1) * Basis(p)
             ELSE
                Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) = Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) + &
                     rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) = Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) + &
                     rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 

                Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) = Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) + &
                     rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) = Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) + &
                     rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)

                IF (dim > 2) THEN

                   Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) = Stiff( (dim+2)*(p-1)+1, (dim+2)*(q-1)+1 ) + &
                        rho * s * Velo(3) * dBasisdx(q,3) * Basis(p) 
                   Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) = Stiff( (dim+2)*(p-1)+2, (dim+2)*(q-1)+2 ) + &
                        rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                   Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) = Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) + &
                        rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                   Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) = Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) + &
                        rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)                  
                   Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) = Stiff( (dim+2)*(p-1)+3, (dim+2)*(q-1)+3 ) + &
                        rho * s * Velo(3) * dBasisdx(q,3) * Basis(p)                       

                END IF

             END IF


             

          END DO

          i = (dim+2) * (p-1) + 1
          F => FORCE(i:i+dim+1)
          F = F + s * LoadAtIP * Basis(p)

          ! The explicit treatment of the convection term grad(v*v) 
          !F(1) = F(1) - s * rho * ( Velo(1) * GradV1(1) + Velo(2) * GradV2(1) )
          !F(2) = F(2) - s * rho * ( Velo(1) * GradV1(2) + Velo(2) * GradV2(2) )


       END DO
    END DO

    !print *, 'Average rot = ', sqrt(rotterm)/sqrt(AK)   
    !rotterm = sqrt(rotterm)/sqrt(AK)  
    !ch = ch * sqrt(rotterm)
    
    !-------------------------------------------------------------------------------------------
    ! The system matrix may have been allocated for the case where all approximation spaces 
    ! are augmented by bubbles. This nullifies the effect of the temperature/density bubbles.
    !-------------------------------------------------------------------------------------------
    DO p = n+1,nd
      i = (dim+2) * p
      FORCE(i)   = 0.0d0
      MASS(:,i)  = 0.0d0
      MASS(i,:)  = 0.0d0
      STIFF(i,:) = 0.0d0
      STIFF(:,i) = 0.0d0
      FORCE(i-1)   = 0.0d0
      MASS(:,i-1)  = 0.0d0
      MASS(i-1,:)  = 0.0d0
      STIFF(i-1,:) = 0.0d0
      STIFF(:,i-1) = 0.0d0

      IF ( .NOT. Stabilization) THEN
         STIFF(i,i) = 1.0d0
         STIFF(i-1,i-1) = 1.0d0
      END IF
    END DO



  !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
  !------------------------------------------------------------------------------










!------------------------------------------------------------------------------
END SUBROUTINE CompressibleNS
!------------------------------------------------------------------------------
