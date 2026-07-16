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
! *  Authors: Juha Ruokolainen, Mika  Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 12 Dec 2003
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> A solver for the incompressible Navier-Stokes equations in rotation form.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StokesSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: AMatrix, MMatrix, PMatrix

  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Hybrid = .FALSE., Found, Convect, &
       GotIt, ConvectionStabilization, BlockPreconditioning = .FALSE., &
       NormalTractionBoundary, GradDivStabilization

  TYPE(Element_t),POINTER :: Element, Parent
  TYPE(Solver_t),POINTER :: PressureSolver, ProjectionSolver, VelocitySolver
  INTEGER, POINTER :: NodeIndexes(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: OuterIterationMethod, NonlinearIterationMethod, eq
  INTEGER :: i, j, k, n, nb, nd, t, istat, dim, m, p, q, &
       MaxIterations, NumberOfNormalTractionNodes
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, Tolerance, ToleranceRatio, &
       atime, stime, at0, NonLinError

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), Mass(:,:), &
       FORCE(:), rho(:), mu(:), Velocity(:,:), &
       ALocal(:,:), MLocal(:,:), PLocal(:,:), PrevSol(:), NonLinRes(:)

  INTEGER, ALLOCATABLE :: TractionBCIndeces(:)

  LOGICAL :: NewTimeStep, OutflowBC, PicardIteration
  INTEGER :: CurrentDoneTime = 0, NonlinearIter, iter, nlen, &
       MaxPicardIterations
  REAL(KIND=dp) :: NonlinearTol, RelaxationFactor, GradDivParam, PresErr, EK, StNumber, &
       NewtonThreshold

  INTEGER, ALLOCATABLE :: Indexes(:)

  SAVE STIFF, LOAD, FORCE, rho, mu, Velocity, AllocationsDone, &
       ALocal, Mass, TractionBCIndeces, NumberOfNormalTractionNodes, &
       CurrentDoneTime, Indexes, MLocal, Norm, &
       PLocal, PrevSol, NonLinRes
  !------------------------------------------------------------------------------

  
  IF (CurrentDoneTime /= Solver % DoneTime) THEN
     CurrentDoneTime = CurrentDoneTime + 1
  END IF

  dim = CoordinateSystemDimension()
  ! Check whether convection stabilization is used...
  ConvectionStabilization = ListGetLogical( Solver % Values, 'Stabilize', GotIt )
  IF (.NOT. GotIt) ConvectionStabilization = .FALSE.
  ConvectionStabilization = .FALSE.

  ! The option for using grad-div stabilization is not available currently... 
  GradDivStabilization = ListGetLogical( Solver % Values, 'Grad-Div Stabilization', GotIt )
  IF ( GotIt .AND. GradDivStabilization) THEN
     GradDivParam =  ListGetConstReal( Solver % Values, &
          'Grad-Div Stabilization Parameter', GotIt )
     IF ( .NOT. GotIt ) GradDivParam = 0.0d0
  ELSE
     GradDivStabilization = .FALSE.
     GradDivParam = 0.0d0
  END IF

  NonlinearIterationMethod = ListGetString(Solver % Values, 'Nonlinear Iteration Method', GotIt)
  IF ( .NOT. GotIt )  NonlinearIterationMethod = 'picard'

  IF ( NonlinearIterationMethod == 'hybrid') THEN
     Newton = .TRUE.
     Hybrid = .TRUE.
     PicardIteration = .FALSE.
  ELSE
     IF ( NonlinearIterationMethod /= 'picard') THEN
        PicardIteration = .FALSE.
        Newton = .TRUE.
        Hybrid = .FALSE.        
     ELSE
        PicardIteration = .TRUE.
        Newton = .FALSE.
        Hybrid = .FALSE.
     END IF
  END IF

  !print *, NonlinearIterationMethod
  !print *, 'Newton = ', Newton
  !print *, 'Hybrid = ', Hybrid
  !print *, 'Picard = ', PicardIteration   


  !--------------------------------------------------------------------------
  ! Check whether the block preconditioning is to be used.
  !---------------------------------------------------------------------------
  BlockPreconditioning = ListGetLogical( Solver % Values, 'Block Preconditioning', GotIt )
  IF (BlockPreconditioning) PRINT *, 'OuterIteration residuals for time level', CurrentDoneTime

  IF (BlockPreconditioning) THEN

     !OuterIterationMethod = ListGetString(Solver % Values, 'Outer Iteration Method', GotIt)
     !IF ( .NOT. GotIt ) OuterIterationMethod = 'gcr'   

     !ToleranceRatio = ListGetConstReal( Solver % Values, &
     !     'Ratio of Convergence Tolerances', GotIt )
     !IF (GotIt) THEN
     !   Tolerance = ToleranceRatio * ListGetConstReal( Solver % Values, &
     !        'Linear System Convergence Tolerance' )
     !ELSE
     !   Tolerance = ListGetConstReal( Solver % Values, &
     !        'Linear System Convergence Tolerance' )
     !END IF

     Tolerance = ListGetConstReal( Solver % Values, &
          'Linear System Convergence Tolerance' )     
     !MaxIterations = ListGetInteger( Solver % Values, &
     !     'Max Outer Iterations', GotIt )
     MaxIterations = ListGetInteger( Solver % Values, &
          'Linear System Max Iterations')    
  END IF


  !----------------------------------------------------------------------------
  !Allocate some permanent storage, this is done first time only:
  !----------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN

     p = Solver % Mesh % MaxElementDOFs       ! The size of the one-component stiffness matrix
     n = (dim+1)*p                            ! The size of the complete stiffness matrix
     m = dim * Solver % Mesh % MaxElementDOFs ! The size of the (1,1)-block (the velocity part)
     !print *, Solver % Matrix % NumberOfRows

     ALLOCATE( &
          FORCE(n), &
          STIFF(n,n), &
          Mass(n,n),&
          LOAD(p,4), &
          rho(p), &
          mu(p), &
          Velocity(dim+1,p), &
          ALocal(m,m), & 
          Indexes(p), &
          MLocal(p,p), &
          PLocal(p,p), &
          PrevSol(Solver % Matrix % NumberOfRows), &
          NonLinRes(Solver % Matrix % NumberOfRows), &          
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'NavierStokesSolver', 'Memory allocation error.' )
     END IF

     IF (BlockPreconditioning) THEN
        !---------------------------------------------------------------------------------
        ! Allocate additional arrays for handling data on the outflow boundary
        !---------------------------------------------------------------------------------
        NumberOfNormalTractionNodes = 0
        i = 0
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           n  = GetElementNOFNodes()
           nd = GetElementDOFs( Indexes )

           IF ( GetElementFamily() == 1 ) CYCLE

           BC => GetBC()
           IF ( ASSOCIATED( BC ) ) THEN
              OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found )

              IF ( .NOT. OutflowBC ) CYCLE
             
              IF ( ANY( Solver % Variable % Perm( Indexes(1:nd) ) == 0 ) ) CYCLE
              i = i + 1
              NumberOfNormalTractionNodes = NumberOfNormalTractionNodes + nd
           END IF
        END DO

        IF (NumberOfNormalTractionNodes > 0) THEN
           ALLOCATE( TractionBCIndeces(NumberOfNormalTractionNodes), &
                STAT=istat )
           IF ( istat /= 0 ) &
                CALL Fatal( 'NavierStokesSolver', 'Memory allocation error.' )
        END IF
     END IF

     AllocationsDone = .TRUE.
  END IF


  !----------------------------------------------------------------------------------
  ! Find the preconditioning matrices and solvers generated before simulation
  !----------------------------------------------------------------------------------
  IF (BlockPreconditioning) THEN
     DO i=1,Model % NumberOfSolvers
        eq = ListGetString( Model % Solvers(i) % Values,'Equation', Found )
        IF (Found) THEN
           nlen = LEN_TRIM(eq)
           IF ( eq(1:nlen) == 'pressure preconditioning' ) THEN
              PressureSolver => Model % Solvers(i)
              PMatrix => PressureSolver % Matrix
           ELSE
              IF ( eq(1:nlen) == 'divergence projection' ) THEN
                 ProjectionSolver => Model % Solvers(i)
                 MMatrix => ProjectionSolver % Matrix
              ELSE
                 IF ( eq(1:nlen) == 'velocity preconditioning' ) THEN
                    VelocitySolver => Model % Solvers(i)
                    AMatrix => VelocitySolver % Matrix                 
                 END IF
              END IF
           END IF
        END IF
     END DO
 
     IF ( .NOT. ASSOCIATED(PMatrix) ) CALL Fatal( 'NavierStokesSolver', &
          'Pressure Preconditioning matrix (S) was not found' )
     IF ( .NOT. ASSOCIATED(MMatrix) ) CALL Fatal( 'NavierStokesSolver', &
          'Pressure Mass matrix (M) was not found' )     
     IF ( .NOT. ASSOCIATED(AMatrix) ) CALL Fatal( 'NavierStokesSolver', &
          'Velocity Preconditioning matrix (A) was not found' )

     n = SIZE(Solver % Variable % Perm)
     Found = ALL(PressureSolver % Variable % Perm(1:n) == Solver % Variable % Perm(1:n))
     IF (.NOT. Found) CALL Fatal( 'NavierStokesSolver', &
          'Incompatible permutations for pressure preconditioning and primary unknown')
     Found = ALL(ProjectionSolver % Variable % Perm(1:n) == Solver % Variable % Perm(1:n))
     IF (.NOT. Found) CALL Fatal( 'NavierStokesSolver', &
          'Incompatible permutations for divergence projection and primary unknown')
     Found = ALL(VelocitySolver % Variable % Perm(1:n) ==  Solver % Variable % Perm(1:n))
     IF (.NOT. Found) CALL Fatal( 'NavierStokesSolver', &
        'Incompatible permutations for divergence projection and primary unknown')

  END IF


  Convect = GetLogical( GetSolverParams(), 'Convective', Found )
  IF ( .NOT. Found ) Convect = .TRUE.
  IF ( .NOT. Convect ) ConvectionStabilization = .FALSE.
  
  atime = CPUTime()
  at0 = RealTime()

  NonlinearIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Max Iterations', minv=0 )
  NonlinearTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance',minv=0.0d0 )
  IF ( .NOT. Convect ) NonlinearIter = 0

  DO iter=1, NonlinearIter+1
     !if (.true.) PRINT *, 'NonlinearIteration residuals for level', iter
     
     NewtonThreshold = ListGetConstReal( Solver % Values, &
          'Nonlinear System Newton After Tolerance', Found )

     IF ( Found .AND. (iter > 1) ) THEN
       IF (NonLinError < NewtonThreshold ) THEN
         Newton = .TRUE.
         PicardIteration = .FALSE.
         Hybrid = .FALSE.
       END IF
     END IF

     MaxPicardIterations = ListGetInteger( Solver % Values, &
          'Nonlinear System Newton After Iterations', Found )    
     IF ( Found ) THEN
        IF ( iter > MaxPicardIterations ) THEN
           Newton = .TRUE.
           PicardIteration = .FALSE.
           Hybrid = .FALSE.  
        END IF
     END IF

     ! Initialize the system matrices and vectors...
     CALL DefaultInitialize()

     IF (BlockPreconditioning) THEN
        CALL CRS_ZeroMatrix( AMatrix )
        CALL CRS_ZeroMatrix( MMatrix )
        CALL CRS_ZeroMatrix( PMatrix )      
     END IF


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
        ! Body forces:
        !-----------------------------------------------
        BodyForce => GetBodyForce()
        LOAD = 0.0d0
        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1:n,1) = GetReal( BodyForce, 'Force 1', Found )
           Load(1:n,2) = GetReal( BodyForce, 'Force 2', Found )
           Load(1:n,3) = GetReal( BodyForce, 'Force 3', Found )
        END IF

        !-----------------------------------------------
        ! Material parameters:
        !-----------------------------------------------
        Material => GetMaterial()
        !rho(1:n) = GetReal( Material, 'Strouhal' )
        !StNumber = rho(1)
        rho(1:n) = GetReal( Material, 'Density' )
        mu(1:n)  = GetReal( Material, 'Viscosity' )
        !-------------------------------------------------------
        ! Get previous elementwise velocity iterate:
        !-------------------------------------------------------
        Velocity = 0.0d0
        IF ( .NOT. Convect ) THEN
           velocity = 0.0d0
        ELSE
           IF (TransientSimulation) THEN
              IF ( iter==1 ) THEN
                 DO i=1,dim
                    Velocity(i,1:nd) = Solver % Variable % PrevValues( &
                         Solver % Variable % DOFs*(Solver % Variable % &
                         Perm(Indexes(1:nd))-1)+i,1)
                 END DO
              ELSE
                 DO i=1,dim
                    Velocity(i,1:nd) = Solver % Variable % Values( &
                         Solver % Variable % DOFs*(Solver % Variable % &
                         Perm(Indexes(1:nd))-1)+i)
                 END DO
              END IF
           ELSE
              DO i=1,dim
                 Velocity(i,1:nd) = Solver % Variable % Values( &
                      Solver % Variable % DOFs*(Solver % Variable % &
                      Perm(Indexes(1:nd))-1)+i)
              END DO
           END IF
        END IF

        !------------------------------------------------------
        ! Get element local matrix and rhs vector:
        !-------------------------------------------------------
        CALL LocalMatrix(  STIFF, Mass, FORCE, LOAD, rho, mu, &
             Velocity, Element, n, nd, dim, ConvectionStabilization, &
             Convect, GradDivParam)

        IF (TransientSimulation) THEN
           CALL Default1stOrderTime( Mass, STIFF, FORCE )
        END IF
        !-----------------------------------------------------------------
        ! Update global matrix and rhs vector from local matrix & vector:
        !----------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )

        !-----------------------------------------------------------------------------
        ! Do the assembly for the preconditioners...
        !----------------------------------------------------------------------------
        IF (BlockPreconditioning) THEN
           Alocal = 0.0d0
           DO p=1,nd
              DO i=1,dim
                 DO q=1,nd
                    DO j=1,dim
                       Alocal( (p-1)*dim+i, (q-1)*dim+j ) = &
                            stiff( (p-1)*(dim+1)+i, (q-1)*(dim+1)+j)
                    END DO
                 END DO
              END DO
           END DO
           CALL UpdateGlobalPreconditioner( AMatrix, ALocal, nd, dim, &
                Solver % Variable % Perm( Indexes(1:nd) ) )
 
           CALL SchurComplementMatrix( rho, mu, Element, n, nd, dim, Solver % dt, &
                MLocal, PLocal)
           CALL UpdateGlobalPreconditioner( MMatrix, MLocal, nd, 1, &
                Solver % Variable % Perm( Indexes(1:nd) ) )
           CALL UpdateGlobalPreconditioner( PMatrix, PLocal, nd, 1, &
                Solver % Variable % Perm( Indexes(1:nd) ) )

        END IF
     END DO
     CALL DefaultFinishBulkAssembly()

     ! Terms on outflow boundary
     !--------------------------------------------------------------
     IF (Convect) THEN
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           n  = GetElementNOFNodes()
           nd = GetElementDOFs( Indexes )

           IF ( GetElementFamily() == 1 ) CYCLE

           BC => GetBC()
           IF ( ASSOCIATED( BC ) ) THEN
              OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found )

              IF (OutflowBC) THEN

                 IF (TransientSimulation) THEN
                    IF ( iter ==1 ) THEN
                       DO i=1,dim
                          Velocity(i,1:nd) = Solver % Variable % PrevValues( &
                               Solver % Variable % DOFs*(Solver % Variable % &
                               Perm(Indexes(1:nd))-1)+i,1)
                       END DO
                    ELSE
                       DO i=1,dim
                          Velocity(i,1:nd) = Solver % Variable % Values( &
                               Solver % Variable % DOFs*(Solver % Variable % &
                               Perm(Indexes(1:nd))-1)+i)
                       END DO
                    END IF
                 ELSE
                    DO i=1,dim
                       Velocity(i,1:nd) = Solver % Variable % Values( &
                            Solver % Variable % DOFs*(Solver % Variable % &
                            Perm(Indexes(1:nd))-1)+i)
                    END DO
                 END IF


                 CALL LocalMatrixBoundary( STIFF, FORCE, &
                      Velocity, rho(1), Element, nd, dim, BlockPreconditioning, ALocal )
                 CALL DefaultUpdateEquations( STIFF, FORCE )

                 IF (BlockPreconditioning) THEN
                    CALL UpdateGlobalPreconditioner( AMatrix, ALocal, nd, dim, &
                         Solver % Variable % Perm( Indexes(1:nd) ) )
                 END IF

              END IF

           END IF
        END DO
     END IF


     CALL DefaultFinishAssembly()
     !CALL DefaultDirichletBCs()


     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 1), &
          1, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)
     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 2), &
          2, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)     
     CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 3), &
          3, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)          
     IF ( dim > 2 ) CALL SetBoundaryConditions(Model, Solver % Matrix, ComponentName(Solver % Variable % name, 4), &
          4, dim+1,  Solver % Variable % Perm, Solver % Matrix % RHS)          
     


     IF (BlockPreconditioning) THEN
        !-------------------------------------------------------------------------
        ! Set boundary conditions for the preconditioners...
        !-------------------------------------------------------------------------
        CALL SetBoundaryConditions(Model, AMatrix, ComponentName(Solver % Variable % name, 1), &
             1, dim,  Solver % Variable % Perm)
        CALL SetBoundaryConditions(Model, AMatrix, ComponentName(Solver % Variable % name, 2), &
             2, dim,  Solver % Variable % Perm)
        IF( dim > 2 ) &
             CALL SetBoundaryConditions(Model, AMatrix, ComponentName(Solver % Variable % name, 3), &
             3, dim,  Solver % Variable % Perm)


        !-----------------------------------------------------------------------------
        ! Boundary conditions for the pressure preconditioner on the part of the boundary where 
        ! the normal component of the surface force vector is given  
        !-----------------------------------------------------------------------------
        IF (NumberOfNormalTractionNodes > 0) THEN
           m = 0

           DO t=1, Solver % Mesh % NumberOfBoundaryElements
              Element => GetBoundaryElement(t)
              IF ( .NOT. ActiveBoundaryElement() ) CYCLE

              n  = GetElementNOFNodes()
              nd = GetElementDOFs( Indexes )
              
              IF ( GetElementFamily() == 1 ) CYCLE

              BC => GetBC()
              IF ( ASSOCIATED( BC ) ) THEN
                 OutflowBC= ListGetLogical( BC, 'Outflow Boundary', Found )

                 IF ( .NOT. OutflowBC ) CYCLE
                 IF ( ANY( Solver % Variable % Perm( Indexes(1:nd) ) == 0 ) ) CYCLE

                 DO i=1,nd
                    k = Solver % Variable % Perm(Indexes(i))
                    TractionBCIndeces(m+i) = k
                    CALL ZeroRow( PMatrix, k)
                    CALL SetMatrixElement( PMatrix, k, k, 1.0d0 )
                    !print *, 'MANIPULATING THE PRESSURE LAPLACIAN...'
                 END DO
                 m = m + nd
              END IF
           END DO
        END IF
     END IF


     atime = CPUTime() - atime

     !---------------------------------------------------------------------------
     ! The solution of the linear system...
     !---------------------------------------------------------------------------
     stime = CPUTime()
     n = Solver % Matrix % NumberOfRows

     IF ( iter == 1 .AND. TransientSimulation ) THEN
        PrevSol(1:n) = Solver % Variable % PrevValues(1:n,1)
     ELSE
        PrevSol(1:n) = Solver % Variable % Values(1:n)
     END IF

     ! Test whether the nonlinear residual is small enough to terminate the nonlinear iteration
     CALL CRS_MatrixVectorMultiply( Solver % Matrix, PrevSol, NonLinRes )
     NonLinRes(1:n) = Solver % Matrix % RHS(1:n) - NonLinRes(1:n)
     NonLinError = PNorm(n,NonLinRes)/PNorm(n,Solver % Matrix % RHS)




     ! IF ( NonLinError < 1.0d-2 ) Newton = .TRUE.
     WRITE(Message, '(a,E12.4)') 'Nonlinear iteration residual: ', NonLinError
     CALL Info( 'NavierStokesSolver', Message, Level=4)

     IF ( (NonLinError < NonLinearTol .OR. iter == (NonlinearIter+1)) .AND. Convect ) THEN
       EXIT
     ELSE
       IF ( BlockPreconditioning ) THEN
         
         PRINT *, 'OuterIteration residuals for nonlinear iteration: ', iter

         CALL GCROuterIteration( Solver % Matrix % NumberOfRows, Solver % Matrix, &
           AMatrix % NumberOfRows, AMatrix, &
           2*(MMatrix % NumberOfRows), MMatrix, PMatrix, &
           Solver % Variable % Values, Solver % Matrix % RHS, &
           MaxIterations, Tolerance, dim, TractionBCIndeces, mu(1), dt ) !       Tolerance
         
       ELSE
         Norm = DefaultSolve()
       END IF

       Norm = PNorm(n, Solver % Variable % Values)
       WRITE(Message, '(a,E12.4)') 'Solution norm: ', Norm
       CALL Info( 'NavierStokesSolver', Message, Level=4)

       stime = CPUTime() - stime
       WRITE(Message, '(a,F8.2)') ' Assembly:  (s)', atime
       CALL Info( 'NavierStokesSolver', Message, Level=4)
       WRITE(Message, '(a,F8.2)') ' Solution:  (s)', stime    
       CALL Info( 'NavierStokesSolver', Message, Level=4)
     END IF

     RelaxationFactor = ListGetCReal( Solver % Values, &
          'Nonlinear System Relaxation Factor', Found )
     IF ( Found ) THEN
        Solver % Variable % Values(1:n) = (1.0d0 - RelaxationFactor) * PrevSol(1:n) + &
             RelaxationFactor * Solver % Variable % Values(1:n)
     END IF

   END DO


   ! Computation of pressure error in the case of test case...
   IF (.FALSE.) THEN
      stime = GetTime()
      PresErr = 0.0d0
      DO t=1,Solver % NumberOfActiveElements

         Element => GetActiveElement(t)
         n  = GetElementNOFNodes()

         nd = GetElementDOFs( Indexes )
         nb = GetElementNOFBDOFs()

         Load(1:nd,4) = Solver % Variable % Values( &
              Solver % Variable % DOFs*(Solver % Variable % &
              Perm(Indexes(1:nd))-1)+dim+1)

         EK = ComputePressureError(Load, Element, n, nd, stime, mu(1))

         PresErr = PresErr + EK

      END DO

      PresErr = SQRT(PresErr)

      PRINT * , 'L2 error of pressure = ', PresErr

      ! Computation of velocity error in the case of test case...
      PresErr = 0.0d0
      DO t=1,Solver % NumberOfActiveElements

         Element => GetActiveElement(t)
         n  = GetElementNOFNodes()

         nd = GetElementDOFs( Indexes )
         nb = GetElementNOFBDOFs()

         DO i=1,dim
            Load(1:nd,i) = Solver % Variable % Values( &
                 Solver % Variable % DOFs*(Solver % Variable % &
                 Perm(Indexes(1:nd))-1)+i)
         END DO

         EK = ComputeVelocityError(Load, Element, n, nd, stime, mu(1))

         PresErr = PresErr + EK

      END DO

      PresErr = SQRT(PresErr)
      PRINT * , 'L2 error of velocity = ', PresErr

   END IF




  
CONTAINS

  !----------------------------------------------------------------------------
  FUNCTION ComputePressureError(ElementSol, Element, n, nd, tp, mu) RESULT(Err)
  !------------------------------------------------------------------------------
    REAL(KIND=dp) :: ElementSol(:,:), Err, tp, mu
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
         DetJ, LoadAtIP(dim+1), Velo(dim), Grad(dim,dim), AK, w1, w2, w3
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: rho = 1.0d0, s, c, c2, ch, rotterm, EK, &
         xp, yp, pref

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    Err = 0.0d0
    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       xp = SUM( Nodes % x(1:n) * Basis(1:n) )
       yp = SUM( Nodes % y(1:n) * Basis(1:n) )
 
       pref = 2 * COS(xp) * SIN( yp + tp )

       !if (.false.) then
       IF (Convect) THEN
          pref = 1.0d-0 * pref + 0.5d0 * ( ( SIN(xp) * SIN(yp+tp) )**2 + ( COS(xp) * COS(yp+tp))**2 )  
       END IF
       
       Err = Err + s * ( pref - SUM( Basis(1:nd) * ElementSol(1:nd,4) ) )**2

    END DO

  !-------------------------------
  END FUNCTION ComputePressureError
  !--------------------------------


  !----------------------------------------------------------------------------
  FUNCTION ComputeVelocityError(ElementSol, Element, n, nd, tp, mu) RESULT(Err)
  !------------------------------------------------------------------------------
    REAL(KIND=dp) :: ElementSol(:,:), Err, tp, mu
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
         DetJ, LoadAtIP(dim+1), Velo(dim), Grad(dim,dim), AK, w1, w2, w3
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: rho = 1.0d0, s, c, c2, ch, rotterm, EK, &
         xp, yp, uref, vref

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    Err = 0.0d0
    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       xp = SUM( Nodes % x(1:n) * Basis(1:n) )
       yp = SUM( Nodes % y(1:n) * Basis(1:n) )
 
       uref = SIN(xp) * SIN( yp + tp )
       vref = COS(xp) * COS( yp + tp )
       
       Err = Err + s * ( uref - SUM( Basis(1:nd) * ElementSol(1:nd,1) ) )**2 + &
            s * ( vref - SUM( Basis(1:nd) * ElementSol(1:nd,2) ) )**2
 
    END DO

  !-------------------------------
  END FUNCTION ComputeVelocityError
  !--------------------------------






  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, Mass, FORCE, LOAD, Nodalrho, &
       Nodalmu, NodalVelo, Element, n, nd, dim, Stabilization, &
       Convect, GradDivParam)
  !------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), Mass(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalVelo(:,:), GradDivParam
    INTEGER :: dim, n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Stabilization, Convect
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
         DetJ, LoadAtIP(dim+1), Velo(dim), Grad(dim,dim), AK, w1, w2, w3
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c, c2, ch, rotterm, mK, hK, Re, VNorm, &
         a1, a2

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    Mass = 0.0d0

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    AK = 0.0d0
    ch = ElementDiameter(Element, Nodes)
    rotterm = 0.0d0

    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       AK = AK + s

       Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
       w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
            SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )

       IF ( (dim > 2) .AND. ( Newton .OR. Stabilization ) ) THEN
          w1 = SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,2) ) - &
               SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,3) )             
          w2 = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,3) ) - &
               SUM( NodalVelo(3,1:nd) * dBasisdx(1:nd,1) )          
       END IF
       
       IF (Stabilization ) THEN
          IF ( dim > 2 ) THEN
             rotterm = rotterm + s * ( ( w2*Velo(3) - w3*Velo(2) )**2 + &
                  ( w3*Velo(1) - w1*Velo(3) )**2 + &
                  ( w1*Velo(2) - w2*Velo(1) )**2 )
          ELSE
             rotterm = rotterm + s * ( w3*w3 * (Velo(1) * Velo(1) + Velo(2) * Velo(2)) )
          END IF
       END IF

       Grad(1,1) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,1) )
       Grad(1,2) = SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )
       Grad(2,1) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) )
       Grad(2,2) = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,2) )      
     
       !----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       !--------------------------------------------
       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:dim+1) = MATMUL( Basis(1:n), LOAD(1:n,1:dim+1) )

       !----------------------------------------------------------------------------------
       ! The system matrix with only the velocity space augmented by bubbles  
       !---------------------------------------------------------------------------------
       DO p=1,nd
          DO q=1,nd
             i = (dim+1) * (p-1) + 1
             j = (dim+1) * (q-1) + 1
             A => STIFF(i:i+dim,j:j+dim)
             M => Mass(i:i+dim,j:j+dim)
             DO i=1,dim
                DO j = 1,dim
                   A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                   A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
                   
                   IF (.FALSE.) THEN
                      A(i,j) = A(i,j) + s * 1.0d-0 * dBasisdx(q,j) * dBasisdx(p,i)
                   END IF

                END DO
                M(i,i) = M(i,i) + s * rho * Basis(p) * Basis(q)            
                !M(i,i) = M(i,i) + s * StNumber * Basis(p) * Basis(q)    
                IF ( Stabilization ) THEN
                   A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)
                ELSE
                   ! Testing w.r.t. all velocity test functions and omitting pressure bubbles...
                   IF (q <= n) &
                        A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)
                END IF

                IF (p <= n) &   
                     ! Testing w.r.t  standard pressure test functions...                
                     ! Pressure bubbles are constructed elsewhere... 
                     A(dim+1,i) = A(dim+1,i) - s * dBasisdx(q,i) * Basis(p)
             END DO


             IF ( Convect ) THEN


                IF ( .FALSE. ) THEN
                   ! The standard convection form in 2d for testing purposes 

                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                        rho * s * Velo(2) * dBasisdx(q,2) * Basis(p) 
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                        rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                   
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                        rho * s * Velo(1) * dBasisdx(q,1) * Basis(p) 
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                        rho * s * Velo(2) * dBasisdx(q,2) * Basis(p)

                ELSE
                   ! The convection term in the rotation form

                   IF ( dim > 2 ) THEN

                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                           s * rho * Velo(3) * dBasisdx(q,3) * Basis(p)
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,1) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) - &
                           s * rho * Velo(3) * dBasisdx(q,1) * Basis(p)                   

                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,2) * Basis(p)
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p) + &
                           s * rho * Velo(3) * dBasisdx(q,3) * Basis(p)
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) - &
                           s * rho * Velo(3) * dBasisdx(q,2) * Basis(p) 


                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,3) * Basis(p)

                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,3) * Basis(p) 

                      Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+3 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p)                  


                      IF ( Newton ) THEN

                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              s * rho * w3 * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+3 ) + &
                              s * rho * w2 * Basis(q) * Basis(p) 

                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) + &
                              s * rho * w3 * Basis(q) * Basis(p)
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+3 ) - &
                              s * rho * w1 * Basis(q) * Basis(p) 

                         Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+1 ) - &
                              s * rho * w2 * Basis(q) * Basis(p)
                         Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+3, (dim+1)*(q-1)+2 ) + &
                              s * rho * w1 * Basis(q) * Basis(p) 


                      END IF


                   ELSE

                      
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) + &
                           s * rho * Velo(2) * dBasisdx(q,2) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                           s * rho * Velo(2) * dBasisdx(q,1) * Basis(p) 

                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                           s * rho * Velo(1) * dBasisdx(q,2) * Basis(p) 
                      Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) + &
                           s * rho * Velo(1) * dBasisdx(q,1) * Basis(p) 

                      
                      IF ( Newton ) THEN

                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              s * rho * w3 * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) + &
                              s * rho * w3 * Basis(q) * Basis(p)                      
                      END IF


                      ! Testing a hybrid linearization strategy...

                      IF ( Hybrid ) THEN
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1 ) - &
                              rho * s * Grad(1,1) * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2 ) - &
                              rho * s * Grad(1,2) * Basis(q) * Basis(p) 
                   
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1 ) - &
                              rho * s * Grad(2,1) * Basis(q) * Basis(p) 
                         Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2 ) - &
                              rho * s * Grad(2,2) * Basis(q) * Basis(p) 
                      END IF



                   END IF

                END IF

             END IF

          END DO

          i = (dim+1) * (p-1) + 1
          F => FORCE(i:i+dim)
          F = F + s * LoadAtIP * Basis(p)

          IF ( Newton ) THEN

             F(1) = F(1) - s * rho * Basis(p) * Velo(2) * w3
             F(2) = F(2) + s * rho * Basis(p) * Velo(1) * w3                    

             ! Testing a hybrid linearization strategy...
             IF ( Hybrid .AND. (dim==2) ) THEN
                F(1) = F(1) - s * rho * Basis(p) * ( Grad(1,1) * Velo(1) + Grad(1,2) * Velo(2) )
                F(2) = F(2) - s * rho * Basis(p) * ( Grad(2,1) * Velo(1) + Grad(2,2) * Velo(2) )
             END IF
             
             IF ( dim > 2 ) THEN
                
                F(1) = F(1) + s * rho * Basis(p) * Velo(3) * w2
                F(2) = F(2) - s * rho * Basis(p) * Velo(3) * w1                                    
                F(3) = F(3) - s * rho * Basis(p) * Velo(1) * w2
                F(3) = F(3) + s * rho * Basis(p) * Velo(2) * w1        

             END IF
          END IF

       END DO
    END DO


    
    !-------------------------------------------------------------------------------------------
    ! The system matrix may have been allocated for the case where both velocities and pressure 
    ! are augmented by bubbles. This nullifies the effect of the pressure bubbles.
    !-------------------------------------------------------------------------------------------
    DO p = n+1,nd
       i = (dim+1) * p
       FORCE(i)   = 0.0d0
       !MASS(:,i)  = 0.0d0
       !MASS(i,:)  = 0.0d0
       STIFF(i,:) = 0.0d0
       !STIFF(:,i) = 0.0d0
       IF ( .NOT. Stabilization ) STIFF(i,i) = 1.0d0
    END DO


    IF ( Stabilization ) THEN
       ! This is the convection stabilization part...
       rotterm = SQRT(rotterm)/SQRT(AK)  
       ch = ch/3.0d0 * SQRT(rotterm)


       hK = element % hK
       mK = element % StabilizationMK
       

       DO t=1,IP % n
          !--------------------------------------------------------------
          ! Basis function values & derivatives at the integration point:
          !--------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis, dBasisdx )

          s = IP % s(t) * detJ

          Velo = MATMUL( NodalVelo(1:dim,1:nd), Basis(1:nd) )
          w3 = SUM( NodalVelo(2,1:nd) * dBasisdx(1:nd,1) ) - &
               SUM( NodalVelo(1,1:nd) * dBasisdx(1:nd,2) )     
          !c = -Velo(2)*dBasisdx(nd,1)+Velo(1)*dBasisdx(nd,2)

           
          !VNorm = MAX( abs(Velo(1)), abs(Velo(2)), 1.0d-12 )
          a1 = MAXVAL( ABS( NodalVelo(1,1:nd) ) )
          a2 = MAXVAL( ABS( NodalVelo(2,1:nd) ) )
          VNorm = MAX( a1, a2 )
          Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
          Vnorm = 1.0d0   ! The velocity norm over the entire domain 
          ch = SQRT(3.0d0/100.0d0) * SQRT( hK * Re * rho / VNorm )


          !a1 = maxval( abs( NodalVelo(1,1:nd) ) )
          !a2 = maxval( abs( NodalVelo(2,1:nd) ) )
          !VNorm = MAX( a1, a2 )
          !Re = hK * VNorm / mu
          !ch = sqrt( 2.0d-2 * hK * 2.0d0 * Re / (1 + Re) )
          

          !----------------------------------------------
          ! The construction of the pressure bubbles 
          !------------------------------------------------
          DO p=n+1,nd
             DO q=n+1,nd
                Stiff( (dim+1)*p, (dim+1)*q ) = Stiff( (dim+1)*p, (dim+1)*q ) + 1.0d0 * &                 
                     s * SUM( dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )                
             END DO

             IF ( .TRUE. ) THEN
                ! The robust version of Picard linearization                
                c = -Velo(2)*dBasisdx(p,1)+Velo(1)*dBasisdx(p,2)
                DO q=1,nd

                   Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) + rho * &          
                        s * dBasisdx(q,1) * c

                   Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) - rho * &          
                        s * dBasisdx(q,2) * c
                
                END DO
             ELSE
                ! The alternative Picard
                DO q=1,nd
                   
                   Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+2 ) - w3 * &          
                        s * Basis(q) * dBasisdx(p,1)
                   
                   Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*p, (dim+1)*(q-1)+1 ) + w3 * &          
                        s * Basis(q) * dBasisdx(p,2)
                   
                END DO
             END IF


          END DO

          !---------------------------------------------------------
          ! The projection part of the equation of motion...
          !--------------------------------------------------------

          DO p=1,nd

             IF ( .TRUE. ) THEN
                ! The robust version of Picard linearization   
                DO q=n+1,nd
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) + &   
                        ch**2 * s * dBasisdx(q,1) * dBasisdx(p,2) * Velo(2)
             
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) - &   
                        ch**2 * s * dBasisdx(q,2) * dBasisdx(p,2) * Velo(1)
             
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) - &   
                        ch**2 * s * dBasisdx(q,1) * dBasisdx(p,1) * Velo(2)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) + &   
                        ch**2 * s * dBasisdx(q,2) * dBasisdx(p,1) * Velo(1)

                END DO

                c2 = rho * ch**2 * dBasisdx(p,2) * ( Velo(1)*Velo(1) + Velo(2)*Velo(2) )
                DO q = 1,nd

                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) + &
                        s * c2 * dBasisdx(q,2)
                
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+2) - &
                        s * c2 * dBasisdx(q,1)   

                END DO

                c2 = rho * ch**2 * dBasisdx(p,1) * ( Velo(1)*Velo(1) + Velo(2)*Velo(2) )
                DO q = 1,nd
                
                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+1) - &
                        s * c2 * dBasisdx(q,2)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) + &
                        s * c2 * dBasisdx(q,1)   
                
                END DO
                
             ELSE
                ! The alternative Picard

                DO q=n+1,nd       
        
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*q) = Stiff( (dim+1)*(p-1)+1, (dim+1)*q) + &   
                        ch**2 * w3 * s * dBasisdx(q,2) * Basis(p)

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*q) = Stiff( (dim+1)*(p-1)+2, (dim+1)*q) - &   
                        ch**2 * w3 * s * dBasisdx(q,1) * Basis(p)

                END DO

                DO q = 1,nd
                   
                   Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) = Stiff( (dim+1)*(p-1)+1, (dim+1)*(q-1)+1) + &
                        w3 * w3 * ch**2 * s * Basis(p) * Basis(q) 

                   Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) = Stiff( (dim+1)*(p-1)+2, (dim+1)*(q-1)+2) + &
                        w3 * w3 * ch**2 * s * Basis(p) * Basis(q) 
                        
                END DO

             END IF
                
          END DO


       END DO


       IF (.FALSE.) THEN
          
          ! Condensate the pressure stabilization bubble dof
          CALL LCondensateStabilizationBubble( n, nd, dim, Stiff)          

          ! Finally nullify the effect of pressure bubbles
          
          DO q=n+1,nd
             i = (dim+1) * q
             STIFF(i,:) = 0.0d0
             STIFF(:,i) = 0.0d0
             STIFF(i,i) = 1.0d0
          END DO

       END IF
       
    END IF
  !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LCondensateStabilizationBubble( n, nd, dim, K )
!------------------------------------------------------------------------------
    USE LinearAlgebra
!------------------------------------------------------------------------------
    INTEGER :: n, nd, dim
    REAL(KIND=dp) :: K(:,:), Kbb(nd-n,nd-n), Kbl(nd-n,(dim+1)*nd-nd+n), &
         Klb((dim+1)*nd-nd+n,nd-n)

    INTEGER :: i, Ldofs((dim+1)*nd-nd+n), Bdofs

    Ldofs = (/ (i, i=1,(dim+1)*n-1) /)
    Bdofs = n*(dim+1)

    Kbb(1,1) = K(Bdofs,Bdofs)
    !Fb(1) = F(Bdofs)
    Kbl(1,1:dim*n-1) = K(Bdofs,Ldofs)
    Klb(1:dim*n-1,1) = K(Ldofs,Bdofs)

    CALL InvertMatrix( Kbb,1 )

    !F(1:dim*n) = F(1:dim*n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:dim*n-1,1:dim*n-1) = &
         K(1:dim*n-1,1:dim*n-1) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensateStabilizationBubble
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, &
      Velocity, rho, Element, nd, dim, BlockPreconditioning, ABlock )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:), ABlock(:,:)
    REAL(KIND=dp) :: Velocity(:,:), rho
    INTEGER :: dim, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: BlockPreconditioning
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), &
        DetJ, LoadAtIP(dim), Velo(dim), Normal(3), SquaredVelo
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    STIFF = 0.0d0
    FORCE = 0.0d0
    ABlock = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    IP = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx )
      s = IP % s(t) * detJ
      
      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      Velo = MATMUL( Velocity(1:dim,1:nd), Basis(1:nd) )
      SquaredVelo = rho * SUM( Velo(1:dim) * Velo(1:dim) ) 

      IF ( Newton ) THEN
         c = 1.0d0
      ELSE
         c = 0.5d0
      END IF

      IF ( .TRUE. ) THEN
         
         DO p=1,nd
            DO i=1,dim
               
               ! The force for the Newton iteration
               IF (Newton) &
                    FORCE( (dim+1)*(p-1)+i ) = FORCE( (dim+1)*(p-1)+i) + s * Normal(i) * Basis(p) * &
                    0.5d0 * SquaredVelo

               DO q=1,nd

                  Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+1 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+1 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(1) * Basis(q) * rho

                  Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+2 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+2 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(2) * Basis(q) * rho             

                  IF ( dim > 2 ) &
                       Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+3 ) = Stiff( (dim+1)*(p-1)+i, (dim+1)*(q-1)+3 )  + &
                       c * s * Basis(p) * Normal(i) * Velo(3) * Basis(q) * rho                 
                  
               END DO

            END DO
         END DO

         ! Change this so that entries are copied only!  
         IF ( BlockPreconditioning) THEN
            DO p=1,nd
               DO q=1,nd
                  DO i=1,dim            
                     
                     ABlock( dim*(p-1)+i, dim*(q-1)+1 ) = ABlock( dim*(p-1)+i, dim*(q-1)+1 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(1) * Basis(q) * rho
                     
                     ABlock( dim*(p-1)+i, dim*(q-1)+2 ) = ABlock( dim*(p-1)+i, dim*(q-1)+2 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(2) * Basis(q) * rho            

                     IF ( dim > 2 ) &
                          ABlock( dim*(p-1)+i, dim*(q-1)+3 ) = ABlock( dim*(p-1)+i, dim*(q-1)+3 )  + &
                          c * s * Basis(p) * Normal(i) * Velo(3) * Basis(q) * rho               

                  END DO
               END DO
            END DO
         END IF



      ELSE
         ! This is the explicit approximation (not in use)
         DO p=1,nd
            DO i=1,dim
               FORCE( (dim+1)*(p-1)+i ) = FORCE( (dim+1)*(p-1)+i) - s * Normal(i) * Basis(p) * &
                    0.5d0 * SquaredVelo
            END DO
         END DO
      END IF


   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE GCROuterIteration( n, A, l, PA, q, PM, PP, x, b, &
      Rounds, TOL, dim, TractionBCIndeces, mu, dt )
!------------------------------------------------------------------------------
!    This is the preconditioned GCR(m) iteration for the linear system Ax=b.
!    The preconditioning strategy is based on the idea of block factorization.
!    The preconditioners are constructed by using inner iterations for systems
!    PA*y=z  and PS*y=z. 
!------------------------------------------------------------------------------  
    TYPE(Matrix_t), POINTER :: A, PA, PM, PP
    INTEGER :: n, l, q, Rounds, dim, TractionBCIndeces(:)
    REAL(KIND=dp) :: x(n), b(n), TOL, mu, dt
    
!-------------------------------------------------------------------------------
    INTEGER :: i, j, k, m, IluOrder, InnerRounds, MaxRestarts
    REAL(KIND=dp) :: r(n), P(n), V(n,Rounds), T(n), T1(n), T2(n), S(n,Rounds), &
        alpha, beta, omega, rho, oldrho, res, norm, stime, tottime, &
        y(l), f(l), InnerTol, u(q/2), g(q/2), h(q), w(q)  
    LOGICAL :: Condition, TrivialPressure, GotIt, ConvergedSol = .FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------

    InnerTol = ListGetConstReal( Solver % Values, &
         'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
         'Linear System Max Iterations') 
    MaxRestarts = ListGetInteger( Solver % Values, &
         'Max GCR Restarts', GotIt) 
    IF ( .NOT. GotIt ) MaxRestarts = 1

    
    !--------------------------------------------------------------------------------
    !    The start of the GCR iteration... 
    !--------------------------------------------------------------------------------     
    tottime = CPUTime()
    CALL Pmv( A, x, r )
    r(1:n) = b(1:n) - r(1:n)

    res = Pnorm(n,r)/Pnorm(n,b)   
    WRITE(Message,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
            0, res  
    CALL Info('NavierStokesSolver',Message,Level=4)
 
    ConvergedSol = ( res < TOL)


    V(1:n,1:Rounds) = 0.0d0
    S(1:n,1:Rounds) = 0.0d0

    IF (.NOT. ConvergedSol) THEN 

       DO k=1,Rounds
          !----------------------------------------------------------
          ! Perform the preconditioning...
          !---------------------------------------------------------------
          T1(1:n) = r(1:n)
          CALL PreconditioningIteration( n, x, A, l, PA, q, PM, PP, T1, dim, InnerTol, &
               InnerRounds, TractionBCIndeces,  mu, dt )
          CALL Pmv( A, T1, T2 )

          !--------------------------------------------------------------
          ! Perform the orthogonalisation of the search directions....
          !--------------------------------------------------------------
          DO i=1,k-1
             beta = PDot( n, V(1:n,i), T2(1:n) )
             T1(1:n) = T1(1:n) - beta * S(1:n,i)
             T2(1:n) = T2(1:n) - beta * V(1:n,i)        
          END DO

          alpha = PNorm(n,T2)
          T1(1:n) = 1.0d0/alpha * T1(1:n)
          T2(1:n) = 1.0d0/alpha * T2(1:n)

          !-------------------------------------------------------------
          ! The update of the solution and save the search data...
          !------------------------------------------------------------- 
          beta = PDot(n, T2, r)
          x(1:n) = x(1:n) + beta * T1(1:n)      
          r(1:n) = r(1:n) - beta * T2(1:n)
          S(1:n,k) = T1(1:n)
          V(1:n,k) = T2(1:n) 

          !--------------------------------------------------------------
          ! Check whether the convergence criterion is met 
          !--------------------------------------------------------------
          res = Pnorm(n,r)/Pnorm(n,b)
          WRITE(Message,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
               k, res     !, CPUTime() - tottime 
          CALL Info('NavierStokesSolver',Message,Level=4)
          ConvergedSol = ( res < TOL)
          IF (ConvergedSol) EXIT

          !print * , 'Velo recompute', ListGetLogical(Solver % Values,  'No Precondition Recompute', stat 
          ! Some performance optimization
          IF (k==1) THEN
             str = ListGetString( VelocitySolver % Values, &
                  'Linear System Preconditioning', GotIt )          
             IF ( GotIt .AND. SEQL(str, 'ilu') ) &
                  CALL ListAddLogical(VelocitySolver % Values, 'No Precondition Recompute', .TRUE.)

             str = ListGetString( ProjectionSolver % Values, &
                  'Linear System Preconditioning', GotIt )          
             IF ( GotIt .AND. SEQL(str, 'ilu') ) &
                  CALL ListAddLogical(ProjectionSolver % Values, 'No Precondition Recompute', .TRUE.)

             str = ListGetString( PressureSolver % Values, &
                  'Linear System Preconditioning', GotIt )          
             IF ( GotIt .AND. SEQL(str, 'ilu') ) &
                  CALL ListAddLogical(PressureSolver % Values, 'No Precondition Recompute', .TRUE.)
          END IF


       END DO

       ! These are related to the performance optimization
       str = ListGetString( VelocitySolver % Values, &
            'Linear System Preconditioning', GotIt )          
       IF ( GotIt .AND. SEQL(str, 'ilu') ) &
            CALL ListAddLogical(VelocitySolver % Values, 'No Precondition Recompute', .FALSE.)     
       str = ListGetString( ProjectionSolver % Values, &
            'Linear System Preconditioning', GotIt )          
       IF ( GotIt .AND. SEQL(str, 'ilu') ) &
            CALL ListAddLogical(ProjectionSolver % Values, 'No Precondition Recompute', .FALSE.)     
       str = ListGetString( PressureSolver % Values, &
            'Linear System Preconditioning', GotIt )          
       IF ( GotIt .AND. SEQL(str, 'ilu') ) &
            CALL ListAddLogical(PressureSolver % Values, 'No Precondition Recompute', .FALSE.)     
 


    END IF


!----------------------------------------------------------------------------------------
  END SUBROUTINE GCROuterIteration
!-----------------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE PreconditioningIteration( n, x, A, l, PA, q, PM, PP, r, dim, &
      TOL, Rounds, TractionBCIndeces, mu, dt )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A, PA, PM, PP
    INTEGER :: n, l, q, Rounds, dim, TractionBCIndeces(:)
    REAL(KIND=dp) :: x(n), r(n), TOL, mu, dt
!-------------------------------------------------------------------------------
    INTEGER :: i, j, m
    REAL(KIND=dp) :: f(l), y(l), z(n), u(q/2), g(q/2), Norm, &
         dvec(q/2), fvec(q/2), tmpvec(q/2), prevpres(q/2)
!----------------------------------------------------------------------------------       
    z(1:n) = x(1:n)

    DO j=1,n/(dim+1)
       DO m=1,dim
          f( dim*(j-1)+m ) = r( (dim+1)*(j-1)+m )
       END DO
    END DO

    y = 0.0d0
    PRINT *, 'Velocity solve...'
    !CALL BiCGStab( l, PA, y, f, Rounds, Tol )    
    CurrentModel % Solver => VelocitySolver
    CALL SolveLinearSystem( PA, f, y, Norm, dim, VelocitySolver )
    CurrentModel % Solver => Solver

    DO j=1,n/(dim+1)
       DO m=1,dim
          z( (dim+1)*(j-1)+m ) = z( (dim+1)*(j-1)+m ) + y( dim*(j-1)+m )
       END DO
    END DO

    CALL CRS_MatrixDivMultiply( A,z,g,dim )

    DO j=1,n/(dim+1)
       fvec(j) = 2.0d0 * mu * g(j)
       prevpres(j) = x((dim+1)*j)
    END DO


    ! Solve the continuous approximation
    dvec = 0.0d0
    PRINT *, 'Pressure mass matrix solve...'
    !CALL BiCGStab( q/2, PM, dvec, fvec, Rounds, Tol )
    CurrentModel % Solver => ProjectionSolver
    CALL SolveLinearSystem( PM, fvec, dvec, Norm, 1, ProjectionSolver )
    CurrentModel % Solver => Solver    

    ! Create rhs for the pressure equation
    tmpvec = 0.0d0
    CALL CRS_MatrixVectorMultiply( PP, dvec, tmpvec )
 
    IF (TransientSimulation) THEN
       DO j=1,q/2
          fvec(j) = -tmpvec(j) - 1.0d0/dt * g(j)
       END DO
    ELSE
       DO j=1,q/2
          fvec(j) = -tmpvec(j)
       END DO
    END IF

    ! Pressure BC update...
    j = SIZE(TractionBCIndeces)
    DO m=1,j
       i = TractionBCIndeces(m)
       fvec(i) = -dvec(i)
    END DO


    dvec = 0.0d0
    PRINT *, 'Pressure Laplace matrix solve...'       
    !CALL BiCGStab( q/2, PP, dvec, fvec, Rounds, Tol )

    CurrentModel % Solver => PressureSolver
    CALL SolveLinearSystem( PP, fvec, dvec, Norm, 1, PressureSolver )
    CurrentModel % Solver => Solver

    DO j=1,n/(dim+1)
      z( (dim+1)*j ) = prevpres(j) + dvec(j)
    END DO
    
    r(1:n) = z(1:n) - x(1:n)
!---------------------------------------------------------------------------------------
  END SUBROUTINE PreconditioningIteration
!---------------------------------------------------------------------------------------







!------------------------------------------------------------------------------
  SUBROUTINE SchurComplementMatrix(  Nodalrho, Nodalmu, Element, &
       n, nd, dim, dt, MLocal, Plocal)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MLocal(:,:), Plocal(:,:) 
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:)
    INTEGER :: dim, n, nd
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), OPTIONAL :: dt    
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), ddBasisddx(nd,3,3), DetJ
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    MLocal = 0.0d0
    PLocal = 0.0d0
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = IP % s(t) * detJ
       !-----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       !---------------------------------------------
       ! the stiffness matrix...
       !----------------------------------------
       DO p=1,n
          DO q=1,n
             MLocal(p,q) = MLocal(p,q) + s * Basis(q) * Basis(p)
             PLocal(p,q) = PLocal(p,q) + s * SUM( dBasisdx(q,1:dim) * dBasisdx(p,1:dim) )
          END DO
       END DO

    END DO

    DO p=n+1,nd
      MLocal(p,p) = 1.0d0
      PLocal(p,p) = 1.0d0
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SchurComplementMatrix
!------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
  SUBROUTINE CRS_MatrixDivMultiply( A,u,v,dim )
!------------------------------------------------------------------------------
!
!  DESCRIPTION:
!    Specific matrix vector product for preconditioning purposes...
!
!  ARGUMENTS:
!
!    TYPE(Matrix_t), POINTER :: A
!    REAL(KIND=dp) :: u(*),v(*)
!
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(*) :: u,v
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)
    
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------

     n = A % NumberOfRows
     Rows   => A % Rows
     Cols   => A % Cols
     Diag   => A % Diag
     Values => A % Values
     
     DO i=1,n
        k = Diag(i)
        IF ( MOD(i,dim+1)==0 ) THEN        
           s = 0.0d0
           DO j=Rows(i),Rows(i+1)-1,dim+1           
              s = s - u(Cols(j)) * Values(j)
           END DO
           DO j=Rows(i)+1,Rows(i+1)-1,dim+1           
              s = s - u(Cols(j)) * Values(j)
           END DO
           IF ( dim > 2 ) THEN
              DO j=Rows(i)+2,Rows(i+1)-1,dim+1           
                 s = s - u(Cols(j)) * Values(j)
              END DO
           END IF

           v( i/(dim+1) ) = s
        END IF
     END DO

!-----------------------------------------------------------------------------
   END SUBROUTINE CRS_MatrixDivMultiply
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     This is an internal ILU preconditioned BiCGStab method.
!------------------------------------------------------------------------------
    SUBROUTINE BiCGStab( n, A, x, b, Rounds, TOL )
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A,M
      INTEGER :: n,Rounds
      REAL(KIND=dp) :: TOL,x(n),b(n),r(n),Ri(n),P(n),V(n),T(n),T1(n),T2(n),S(n)
!------------------------------------------------------------------------------
      INTEGER :: i
      LOGICAL :: Condition
      REAL(KIND=dp) :: alpha,beta,omega,rho,oldrho, res, tottime
!------------------------------------------------------------------------------
      IF ( ALL(x == 0.0d0 ) ) x = b      

      CALL Pmv( A, x, r )
      r(1:n) = b(1:n) - r(1:n)

      Ri(1:n) = r(1:n)
      P(1:n) = 0
      V(1:n) = 0
      omega  = 1
      alpha  = 0
      oldrho = 1
      tottime = CPUTime()

      DO i=1,Rounds
        rho = Pdot( n, r, Ri )
        beta = alpha * rho / ( oldrho * omega )
        P(1:n) = r(1:n) + beta * (P(1:n) - omega*V(1:n))
        V(1:n) = P(1:n)
        CALL CRS_LUSolve( n, A, V )
        T1(1:n) = V(1:n)
        CALL Pmv( A, T1, V )
        alpha = rho / Pdot( n, Ri, V )
        S(1:n) = r(1:n) - alpha * V(1:n)
        T(1:n) = S(1:n)
        CALL CRS_LUSolve( n, A, T )
        T2(1:n) = T(1:n)
        CALL Pmv( A, T2, T )
        omega = Pdot( n,T,S ) / Pdot( n,T,T )
        oldrho = rho
        r(1:n) = S(1:n) - omega*T(1:n)
        x(1:n) = x(1:n) + alpha*T1(1:n) + omega*T2(1:n)
        res = Pnorm(n,r)/Pnorm(n,b)
        WRITE(*,'(I4,ES12.3)') i, res
        IF ( res < TOL ) THEN
          EXIT
        END IF
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE BiCGStab
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE Pmv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       CALL CRS_MatrixVectorMultiply( A, x, b )
!------------------------------------------------------------------------------
    END SUBROUTINE Pmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION Pnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s,x(:)
!------------------------------------------------------------------------------
       s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
!------------------------------------------------------------------------------
    END FUNCTION Pnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION Pdot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s,x(:),y(:)
!------------------------------------------------------------------------------
       s = DOT_PRODUCT( x(1:n), y(1:n) )
!------------------------------------------------------------------------------
    END FUNCTION Pdot
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE UpdateGlobalPreconditioner( StiffMatrix, LocalStiffMatrix, &
      n, NDOFs, NodeIndexes )
!------------------------------------------------------------------------------
! 
! Add element matrices to global matrices
!
! TYPE(Matrix_t), POINTER :: StiffMatrix
!   INOUT: The global matrix
!
! REAL(KIND=dp) :: LocalStiffMatrix(:,:)
!   INPUT: Local matrix to be added to the global matrix
!
! INTEGER :: n, NDOFs
!   INPUT :: number of nodes / element and number of DOFs / node
!
! INTEGER :: NodeIndexes(:)
!   INPUT: Element node to global node numbering mapping
! 
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix

     REAL(KIND=dp) :: LocalStiffMatrix(:,:)

     INTEGER :: n, NDOFs, NodeIndexes(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
!------------------------------------------------------------------------------
!    Update global matrix .
!------------------------------------------------------------------------------
     SELECT CASE( StiffMatrix % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_GlueLocalMatrix( StiffMatrix,n,NDOFs, NodeIndexes, &
                        LocalStiffMatrix )

     CASE( MATRIX_BAND,MATRIX_SBAND )
       CALL Band_GlueLocalMatrix( StiffMatrix,n,NDOFs, NodeIndexes, &
                        LocalStiffMatrix )

     END SELECT

!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalPreconditioner
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sets internally Dirichlet boundary condition for given dof.
!------------------------------------------------------------------------------
   SUBROUTINE SetBoundaryConditions( Model, StiffMatrix, &
                   Name, DOF, NDOFs, Perm, rhs )
!------------------------------------------------------------------------------
!
! TYPE(Model_t) :: Model
!   INPUT: the current model structure
!
! TYPE(Matrix_t), POINTER :: StiffMatrix
!   INOUT: The global matrix
!
! CHARACTER(LEN=*) :: Name
!   INPUT: name of the dof to be set
!
! INTEGER :: DOF, NDOFs
!   INPUT: The order number of the dof and the total number of DOFs for
!          this equation
!
! INTEGER :: Perm(:)
!   INPUT: The node reordering info, this has been generated at the
!          beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Matrix_t), POINTER :: StiffMatrix

    CHARACTER(LEN=*) :: Name 
    INTEGER :: DOF, NDOFs, Perm(:)
    REAL(KIND=dp), OPTIONAL :: rhs(:)    
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,t,k1,k2,nd
    LOGICAL :: GotIt, periodic
    REAL(KIND=dp) :: Work(Model % MaxElementNodes),s

!------------------------------------------------------------------------------

    DO t = Model % NumberOfBulkElements + 1, &
        Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

      CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
!      Set the current element pointer in the model structure to
!      reflect the element being processed
!------------------------------------------------------------------------------
      Model % CurrentElement => Model % Elements(t)
!------------------------------------------------------------------------------
      !n = CurrentElement % TYPE % NumberOfNodes

      n  = GetElementNOFNodes()
      nd = GetElementDOFs( Indexes )
      NodeIndexes => CurrentElement % NodeIndexes(1:n)

      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            Model % BCs(i) % Tag ) THEN

          Work(1:n) = ListGetReal( Model % BCs(i) % Values, &
              Name,n,NodeIndexes, gotIt )
          IF ( gotIt ) THEN
            DO j=1,n
              k = Perm(NodeIndexes(j))
              IF ( k > 0 ) THEN
                k = NDOFs * (k-1) + DOF
                CALL ZeroRow( StiffMatrix,k )
                CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 )
                IF ( PRESENT(rhs) ) rhs(k) = work(j) 
              END IF
            END DO

            DO j=n+1,nd
               k = Perm(Indexes(j))
               k = NDOFs * (k-1) + DOF              
               CALL ZeroRow( StiffMatrix,k )
               CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 )
               IF ( PRESENT(rhs) ) rhs(k) = 0.0d0  
            END DO

          END IF
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetBoundaryConditions
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE StokesSolver
!------------------------------------------------------------------------------
