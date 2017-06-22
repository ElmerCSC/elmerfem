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
! *  Authors: Juha Ruokolainen, Mika Malinen
! *  Email:   mika.malinen@csc.fi & Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2012-01-30
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE StokesSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  IF ( .NOT. ListCheckPresent(SolverParams, 'Bubbles in Global System') ) &
      CALL ListAddLogical(SolverParams, 'Bubbles in Global System', .FALSE.)  

  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Solver') ) &
      CALL ListAddString(SolverParams, 'Linear System Solver', 'Iterative')
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Iterative Method') ) &
      CALL ListAddString(SolverParams, 'Linear System Iterative Method', 'GCR')
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System GCR Restart') ) &
      CALL ListAddInteger(SolverParams, 'Linear System GCR Restart', 50)
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Max Iterations') ) &
      CALL ListAddInteger(SolverParams, 'Linear System Max Iterations', 200)
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Row Equilibration') ) &
      CALL ListAddLogical(SolverParams, 'Linear System Row Equilibration', .TRUE.)
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Convergence Tolerance') ) &
      CALL ListAddConstReal(SolverParams, 'Linear System Convergence Tolerance', 1.0d-6)
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Base Tolerance') ) &
      CALL ListAddConstReal(SolverParams, 'Linear System Base Tolerance', 1.0d-3)
  IF ( .NOT. ListCheckPresent(SolverParams, 'Linear System Relative Tolerance') ) &
      CALL ListAddConstReal(SolverParams, 'Linear System Relative Tolerance', 1.0d-2)

!------------------------------------------------------------------------------
END SUBROUTINE StokesSolver_Init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE StokesSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  A parallel solver that uses two-level iterations to solve the discrete
!  Stokes model. Inner iterations can be associated with preconditioning and 
!  they provide search directions for the outer iterative method (GCR) applied 
!  to the primitive Stokes problem. 
!
!  A key design choice here has been that the inner iterations are performed
!  via calling DefaultSolve routine, so that the full range of standard parallel
!  methods can be applied for this computation and a clean interface is obtained
!  with this piece of code. The outer iterative method is however implemented
!  locally into this solver.
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils
  USE SolverUtils
  USE ElementUtils
  USE MaterialModels

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found, Convect, GotForceBC, NormalTangential, &
       SkipPowerLaw=.TRUE., Newton = .FALSE.
  TYPE(Element_t), POINTER :: Element

  REAL(KIND=dp) :: Norm, NonlinearTol, NewtonThreshold, NonLinError, res
  INTEGER, POINTER :: EdgeMap(:,:)
  INTEGER :: BrickFaceMap(6,4)
  INTEGER :: n, m, p, q, nb, nd, t, istat, active, dim, &
       iter, NonlinearIter, MaxPicardIterations, MidEdgeNodes(12)
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), FORCE(:), rho(:), mu(:), &
       Velocity(:,:), LoadVector(:,:), SlipCoeff(:,:), ExtPressure(:)

  TYPE(Nodes_t) :: ElementNodes
  CHARACTER(LEN=MAX_NAME_LEN) :: NormalTangentialName

  SAVE STIFF, LOAD, FORCE, rho, mu, Velocity, AllocationsDone, ElementNodes, &
      LoadVector, SlipCoeff, ExtPressure


!-----------------------------------------------------------------------------
! Variables for two-level iterations... 
!----------------------------------------------------------------------------- 
  LOGICAL :: BlockPreconditioning, Parallel, UpdateMatrix, UseTrueResidual, Timing, &
       OptimizeBW, GlobalBubbles, P2P1
  LOGICAL :: DoScaling, DoEquilibration, BlockDiagonalA, UseVeloLaplacian, AdaptiveTols
  CHARACTER(LEN=MAX_NAME_LEN) :: Eq
  INTEGER :: i, j, k, nlen, Round, MaxIterations, RestartM
  INTEGER, ALLOCATABLE :: Indexes(:)
  INTEGER, POINTER :: NodeIndexes(:), Permutation(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: MinTolerance, t0, rt0, st, rst, ct, maxU=0.0d0, maxV=0.0d0, &
       BaseTolerance, TargetTol, RelTolerance, PrecondTol
#else
  REAL(KIND=dp) :: MinTolerance, t0, rt0, st, rst, ct, CPUTime, RealTime, maxU=0.0d0, maxV=0.0d0, &
       BaseTolerance, TargetTol, RelTolerance, PrecondTol
#endif
  TYPE(Solver_t), POINTER :: PrecondSolver, PressureSolver, VelocitySolver, SubSolver, &
       GradSolver1, GradSolver2, GradSolver3
  TYPE(Matrix_t), POINTER :: PrecondMatrix, MMatrix, PMatrix, AMatrix, B1Matrix, B2Matrix, B3Matrix
  REAL(KIND=dp), ALLOCATABLE :: Snew(:), R(:), S(:,:), V(:,:), &
       ALocal(:,:), PLocal(:,:), B1Local(:,:), B2Local(:,:), B3Local(:,:), Vx(:), Vy(:), Vz(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: Residual(:), RotatedVarValues(:), TempRes(:), TempS(:), &
      TempMx(:), TempMb(:), TempMr(:), TempRHS(:)
  REAL(KIND=dp), POINTER :: Mx(:), Mb(:), Mr(:)
  SAVE Snew, R, S, V, Indexes, ALocal, PLocal, B1Local, B2Local, B3Local, Vx, Vy, Vz
  REAL(KIND=dp), ALLOCATABLE :: RowScalingA(:), RowScalingP(:)
  SAVE RowScalingA, RowScalingP
  SAVE MaxIterations,RestartM
!------------------------------------------------------------------------------
  NonLinError=2.0_dp
  Solver % Variable % NonLinChange = NonLinError

  Parallel = ParEnv % Pes > 1

  dim = CoordinateSystemDimension()

  P2P1 = ListGetLogical( Solver % Values, 'P2-P1 Approximation', Found)
  BlockPreconditioning = ListGetLogical( Solver % Values, 'Block Preconditioning', Found)
  IF ( .NOT. Found) BlockPreconditioning = .FALSE.

  Convect = GetLogical( GetSolverParams(), 'Convective', Found )
  IF ( .NOT. Found ) Convect = .FALSE. 

  SkipPowerLaw = ListGetLogical( Solver % Values, 'Constant-Viscosity Start', Found)
  IF ( .NOT. Found) SkipPowerLaw = .TRUE.

  DoScaling=ListGetLogical(Solver % Values, 'Linear System Scaling', Found) 
  DoEquilibration=ListGetLogical(Solver % Values, 'Linear System Row Equilibration', Found)
  DoScaling = DoScaling .OR. DoEquilibration

  BlockDiagonalA = ListGetLogical( Solver % Values, 'Block Diagonal A', Found)
  UseVeloLaplacian = ListGetLogical( Solver % Values, 'Use Velocity Laplacian', Found)

  AdaptiveTols = ListGetLogical(Solver % Values, 'Linear System Adaptive Tolerance', Found)
  IF (AdaptiveTols) THEN
     TargetTol = ListGetConstReal(Solver % Values, 'Linear System Convergence Tolerance')
     IF (.NOT. ListCheckPresent(Solver % Values, 'Linear System Relative Tolerance')) THEN
        RelTolerance = 1.0d-2
     ELSE
        RelTolerance = ListGetConstReal(Solver % Values, 'Linear System Relative Tolerance')
     END IF
     
     IF (.NOT. ListCheckPresent(Solver % Values, 'Linear System Base Tolerance')) THEN
        BaseTolerance = 1.0d-4
     ELSE
        BaseTolerance = ListGetConstReal(Solver % Values, 'Linear System Base Tolerance')
     END IF
  END IF

  !----------------------------------------------------------------------------------
  ! Find the preconditioning matrices and solvers generated before simulation
  !----------------------------------------------------------------------------------
  IF (BlockPreconditioning) THEN
     DO i=1,Model % NumberOfSolvers
        Eq = ListGetString( Model % Solvers(i) % Values, 'Equation', Found )
        IF (Found) THEN
           nlen = LEN_TRIM(eq)
           SELECT CASE( Eq(1:nlen) )
           CASE ('pressure preconditioning')
              PressureSolver => Model % Solvers(i)
              PMatrix => PressureSolver % Matrix
           CASE ('velocity preconditioning')
              VelocitySolver => Model % Solvers(i)
              AMatrix => VelocitySolver % Matrix
           END SELECT
        END IF
     END DO

     IF ( .NOT. ASSOCIATED(PressureSolver) ) &
          CALL Fatal( 'ParStokes', 'Pressure preconditioning operation missing...' )     

     IF ( .NOT. ASSOCIATED(VelocitySolver) ) &
          CALL Fatal( 'ParStokes', 'Velocity preconditioning operation missing...' ) 

     !---------------------------------------------------------------------------
     ! Perform permutation check...
     !---------------------------------------------------------------------------
     n = SIZE(Solver % Variable % Perm)
     IF ( ANY(Solver % Variable % Perm(1:n) /= PressureSolver % Variable % Perm(1:n)) ) &
          CALL Fatal( 'ParStokes', 'Nonmatching variable permutations, &
          use Optimize Bandwidth in pressure preconditioning' )
     IF ( ANY(Solver % Variable % Perm(1:n) /= VelocitySolver % Variable % Perm(1:n)) ) &
          CALL Fatal( 'ParStokes', 'Nonmatching variable permutations, &
          use Optimize Bandwidth in velocity preconditioning' )
  END IF


  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     p = Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     n = p * (dim+1)
     m = p * dim
     ALLOCATE( &
          FORCE(n), &
          LOAD(p,dim+1),  &
          STIFF(n,n), &
          Indexes(p), &
          rho(p), &
          mu(p), &
          Velocity(dim,p), &
          Vx(p), &
          Vy(p), &
          Vz(p), &
          LoadVector(4,p), &
          SlipCoeff(3,p), &
          ExtPressure(p), &
          ALocal(m,m), & 
          PLocal(p,p), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'ParStokes', 'Memory allocation error.' )
     END IF

     IF (BlockPreconditioning) THEN
        MaxIterations = ListGetInteger(Solver % Values, 'Linear System Max Iterations')
        ReStartM = ListGetInteger(Solver % Values, 'Linear System GCR Restart', Found)
        IF ( .NOT. Found) THEN
           ReStartM = MaxIterations
        END IF
        ALLOCATE( &
             Snew( Solver % Matrix % NumberOfRows ), &    ! Search direction variable
             R( Solver % Matrix % NumberOfRows ),   &     ! For saving the residual generated by the iterative method
             S( Solver % Matrix % NumberOfRows, RestartM), &      ! For saving search directions
             V( Solver % Matrix % NumberOfRows, RestartM), &      ! For saving the range of coefficient matrix
             STAT=istat )
        IF ( istat /= 0 ) THEN
           CALL Fatal( 'ParStokes', 'Memory allocation error.' )
        END IF
     END IF

     AllocationsDone = .TRUE.
  END IF

  !-------------------------------------------------------------------
  ! System assembly & solve for each system arising from linearization
  !--------------------------------------------------------------------
  NonlinearIter = ListGetInteger( Solver % Values, &
       'Nonlinear System Max Iterations', minv=0 )
  NonlinearTol = ListGetConstReal( Solver % Values, &
       'Nonlinear System Convergence Tolerance', minv=0.0d0 )

  IF (BlockPreconditioning) THEN
     !Solver % Variable % Values = 0.0d0 !Fab: No need to reset to 0
     Solver % Matrix % RHS = 0.0d0
     NonlinearIter = NonlinearIter + 1   ! To enable the computation of true nonlinear residual
  END IF

  Active = GetNOFActive()

  DO iter=1, NonlinearIter
     !-----------------------------------------------------------------------
     ! Check whether it is the time to switch to Newton linearization...
     !-----------------------------------------------------------------------
     Newton = .FALSE.
     IF (iter > 2) THEN
       NewtonThreshold = ListGetConstReal( Solver % Values, &
           'Nonlinear System Newton After Tolerance', Found )
       IF ( Found .AND. (NonLinError < NewtonThreshold) ) THEN
         Newton = .TRUE.
       ELSE
         MaxPicardIterations = ListGetInteger( Solver % Values, &
             'Nonlinear System Newton After Iterations', Found )    
         IF ( Found ) THEN
           IF ( iter > MaxPicardIterations ) THEN
             Newton = .TRUE.
           END IF
         END IF
       END IF
     END IF

     IF (Newton) THEN
        WRITE(Message,'(a)') 'The Newton linearization is used...'
        CALL Info('ParSolver', Message, Level=4)
     END IF
 
     !------------------------------------------------------------------
     ! Initialize matrix structures...
     !------------------------------------------------------------------
     CALL DefaultInitialize()

     IF (BlockPreconditioning) THEN
       CALL InitializeToZero( AMatrix, AMatrix % RHS ) 
       CALL InitializeToZero( PMatrix, PMatrix % RHS )        
     END IF

     !------------------------------------------------------------
     ! We need to make a copy of the rotated solution vector to 
     ! handle normal-tangential bcs.
     !------------------------------------------------------------
     IF ( BlockPreconditioning .AND. Iter==1 ) THEN     
       n =  Solver % Matrix % NumberOfRows
       ALLOCATE( RotatedVarValues(n) )
       RotatedVarValues = 0.0d0          
     END IF

     !---------------------------------------------------------------
     ! We always generate the initial guess by adjusting material law
     ! exponent to obtain linear behaviour
     !---------------------------------------------------------------
     IF ((Iter==1) .AND. (SkipPowerlaw)) THEN !SkipPowerlaw only 1st time
        SkipPowerlaw = .TRUE.
     ELSE
        SkipPowerlaw = .FALSE.
     END IF

     CALL StartAdvanceOutput( 'StokesSolver', 'Assembly:' )
     DO t=1,Active
        CALL AdvanceOutput(t, Active)        

        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementDOFs( Indexes )
        nb = GetElementNOFBDOFs()

        !-----------------------------------------------
        ! Body forces:
        !-----------------------------------------------
        LOAD = 0.0d0
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1:n,1) = GetReal( BodyForce, 'Flow BodyForce 1', Found )
           Load(1:n,2) = GetReal( BodyForce, 'Flow BodyForce 2', Found )
           IF (dim > 2) &
                Load(1:n,3) = GetReal( BodyForce, 'Flow BodyForce 3', Found )
        END IF

        !-----------------------------------------------
        ! Material parameters:
        !-----------------------------------------------
        Material => GetMaterial()
        rho(1:n) = GetReal( Material, 'Density' )
        mu(1:n)  = GetReal( Material, 'Viscosity' )
        !-------------------------------------------------------    

        !--------------------------------------------------------------------------
        ! Get previous elementwise velocity iterate for approximating nonlinearity
        !--------------------------------------------------------------------------
        Vx = 0.0d0
        Vy = 0.0d0
        Vz = 0.0d0
        CALL GetScalarLocalSolution( Vx, ComponentName(Solver % Variable % Name,1) )
        CALL GetScalarLocalSolution( Vy, ComponentName(Solver % Variable % Name,2) )      
        IF( dim > 2 ) THEN
          CALL GetScalarLocalSolution( Vz, ComponentName(Solver % Variable % Name,3) )
        END IF

        !---------------------------------------------------------------------
        ! Get element local matrix and rhs vector:
        !---------------------------------------------------------------------
        CALL LocalMatrix( STIFF, ALocal, PLocal, FORCE, LOAD, rho, mu, &
             Vx, Vy, Vz, Element, n, nd+nb, dim, Convect, P2P1, SkipPowerLaw, Newton, &
             BlockDiagonalA)

        IF ( nb > 0 ) THEN
           CALL StaticCondensation( nd, nb, dim, STIFF, FORCE )
        END IF

        IF (BlockPreconditioning) THEN
          IF (.NOT. BlockDiagonalA) THEN
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
          ELSE
            ! The block diagonal approximation of the velocity preconditioning matrix: 
	    IF (.NOT. UseVeloLaplacian) THEN
              Alocal = 0.0d0
              DO p=1,nd
                DO i=1,dim
                  DO q=1,nd
                    Alocal( (p-1)*dim+i, (q-1)*dim+i ) = &
                        stiff( (p-1)*(dim+1)+i, (q-1)*(dim+1)+i)
                  END DO
                END DO
              END DO
	    END IF
          END IF
        END IF

        CALL DefaultUpdateEquations( STIFF, FORCE )

        IF (BlockPreconditioning) THEN
          Force = 0.0d0
          CALL DefaultUpdateEquations( ALocal, FORCE, USolver=VelocitySolver)
          CALL DefaultUpdateEquations( PLocal, FORCE, USolver=PressureSolver)
        END IF
 
     END DO
     
     CALL DefaultFinishBulkAssembly()
     
     !------------------------------------------------------------------------------
     ! Surface force and slip boundary conditions
     !------------------------------------------------------------------------------
     DO t = 1, GetNOFBoundaryElements()

        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE
        
        n = GetElementNOFNodes()
        nd = GetElementDOFs( Indexes )
      
        IF ( GetElementFamily() == 1 ) CYCLE

        CALL GetElementNodes( ElementNodes )
        NodeIndexes => Element % NodeIndexes

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

        GotForceBC = GetLogical( BC, 'Flow Force BC', Found )
        IF ( .NOT. Found ) GotForceBC = .TRUE.

        IF ( GotForceBC ) THEN
           STIFF = 0.0d0
           FORCE = 0.0d0

           ExtPressure = 0.0d0
           ExtPressure(1:n) = GetReal( BC, 'Normal Surface Traction', GotForceBC )
           IF (.NOT. GotForceBC) ExtPressure(1:n) = GetReal( BC, 'External Pressure', GotForceBC )

           !------------------------------------------------------------------------------
           ! Surface force in given direction BC: \sigma n = F
           !------------------------------------------------------------------------------
           LoadVector = 0.0d0
           LoadVector(1,1:n) = GetReal( BC, 'Surface Traction 1', Found )
           IF (.NOT. Found) LoadVector(1,1:n) = GetReal( BC, 'Pressure 1', Found )
           LoadVector(2,1:n) = GetReal( BC, 'Surface Traction 2', Found )
           IF (.NOT. Found) LoadVector(2,1:n) = GetReal( BC, 'Pressure 2', Found )
           LoadVector(3,1:n) = GetReal( BC, 'Surface Traction 3', Found )
           IF (.NOT. Found) LoadVector(3,1:n) = GetReal( BC, 'Pressure 3', Found )
           LoadVector(4,1:n) =  0.0d0

           !------------------------------------------------------------------------------
           ! Slip boundary condition BC: \sigma n = R_k u_k
           !------------------------------------------------------------------------------
           SlipCoeff = 0.0d0
           SlipCoeff(1,1:n) =  GetReal( BC, 'Slip Coefficient 1', Found )
           SlipCoeff(2,1:n) =  GetReal( BC, 'Slip Coefficient 2', Found )
           SlipCoeff(3,1:n) =  GetReal( BC, 'Slip Coefficient 3', Found )


           NormalTangentialName = 'Normal-Tangential'
           NormalTangentialName = TRIM(NormalTangentialName) // ' ' // &
                GetVarName(Solver % Variable)

           NormalTangential = GetLogical( BC, &
                NormalTangentialName, Found )
         
           CALL NavierStokesBoundary( STIFF, FORCE, &
                LoadVector, ExtPressure, SlipCoeff, NormalTangential,   &
                Element, n, nd, ElementNodes )

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
           END IF

           CALL DefaultUpdateEquations( STIFF, FORCE )

           IF (BlockPreconditioning) THEN
             Force = 0.0d0
             CALL DefaultUpdateEquations( ALocal, FORCE, USolver=VelocitySolver)
           END IF

        END IF
     END DO

     CALL DefaultFinishAssembly()

     CALL DefaultDirichletBCs()

     IF (BlockPreconditioning) THEN
        ! CALL DefaultFinishAssembly(VelocitySolver)
        ! CALL DefaultFinishAssembly(PressureSolver)
        CALL DefaultDirichletBCs(USolver=PressureSolver)
        CALL DefaultDirichletBCs(USolver=VelocitySolver)
     END IF


     !---------------------------
     ! Perform linear solve...
     !----------------------------
     IF (.NOT. BlockPreconditioning) THEN
        Norm = DefaultSolve()
        IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
     ELSE
        !-------------------------------------------------------------------------------
        ! This branch performs a locally implemented block preconditioned GCR iteration.
        ! First, perform linear system scaling by row equilibration and make some
        ! initializations...
        !-------------------------------------------------------------------------------
        IF ( DoScaling ) THEN
           IF ( DoEquilibration ) THEN
              CALL RowEquilibration(Solver % Matrix,Solver % Matrix % RHS,Parallel)
           ELSE
              CALL ScaleLinearSystem(Solver,Solver % Matrix,Solver % Matrix % RHS)
           END IF
           ! store the scaling for the blocks A and Q
           j = PressureSolver % Matrix % NumberOfRows
           IF (.NOT. ALLOCATED(RowScalingP)) THEN
              ALLOCATE(RowScalingP(j))
           END IF
           RowScalingP=1.0_dp
           IF (.NOT. ALLOCATED(RowScalingA)) THEN
              ALLOCATE(RowScalingA(AMatrix%NumberOfRows))
           END IF
           RowScalingA=1.0_dp
           ! for row equilibration we do not scale the A and Q matrices
           ! so that their symmetry is not destroyed, instead we scale 
           ! the RHS back to the symmetric scaling
           IF (DoEquilibration) THEN
              DO i=1,j
                 RowScalingP(i) = 1.0_dp/Solver%Matrix%DiagScaling((dim+1)*i)
                 DO p=1,dim
                    RowScalingA( dim*(i-1)+p ) = &
                         1.0_dp/Solver%Matrix%DiagScaling( (dim+1)*(i-1)+p )
                 END DO
              END DO
           ELSE
              IF (.NOT. ASSOCIATED(PMatrix%DiagScaling)) THEN
                 ALLOCATE(PMatrix%DiagScaling(j))
              END IF
              PMatrix%DiagScaling(:)=1.0_dp
              IF (.NOT. ASSOCIATED(AMatrix % DiagScaling)) THEN
                 ALLOCATE(AMatrix%DiagScaling(AMatrix%NumberOfRows))
              END IF
              AMatrix%DiagScaling(:)=1.0_dp
              DO i=1,j
                 PMatrix%DiagScaling(i) = Solver%Matrix%DiagScaling((dim+1)*i)
                 DO p=1,dim
                    AMatrix%DiagScaling( dim*(i-1)+p ) = &
                         Solver%Matrix%DiagScaling( (dim+1)*(i-1)+p )
                 END DO
              END DO
              CALL ScaleLinearSystem( VelocitySolver, AMatrix, DiagScaling=AMatrix%DiagScaling)
              CALL ScaleLinearSystem( PressureSolver, PMatrix, DiagScaling=PMatrix%DiagScaling)
           END IF

        ELSE
           j = PressureSolver % Matrix % NumberOfRows
           IF (.NOT. ALLOCATED(RowScalingP)) THEN
             ALLOCATE(RowScalingP(j))
           END IF
           RowScalingP=1.0_dp
           IF (.NOT. ALLOCATED(RowScalingA)) THEN
             ALLOCATE(RowScalingA(AMatrix%NumberOfRows))
           END IF
           RowScalingA=1.0_dp

           IF (.NOT. ASSOCIATED(PMatrix%DiagScaling)) THEN
             ALLOCATE(PMatrix%DiagScaling(j))
           END IF
           PMatrix%DiagScaling(:)=1.0_dp
           IF (.NOT. ASSOCIATED(AMatrix % DiagScaling)) THEN
             ALLOCATE(AMatrix%DiagScaling(AMatrix%NumberOfRows))
           END IF
           AMatrix%DiagScaling(:)=1.0_dp

           IF (.NOT. ASSOCIATED(Solver % Matrix % DiagScaling)) THEN
             ALLOCATE(Solver % Matrix % DiagScaling(Solver%Matrix%NumberOfRows))
           END IF
           Solver % Matrix % DiagScaling(:)=1.0_dp

        END IF

        !-------------------------------------------------------------------------------
        ! Switch to using the rotated variables (again)
        !-------------------------------------------------------------------------------
        IF (Iter > 1) THEN
          j = Solver % Matrix % NumberOfRows
          Solver % Variable % Values(1:j) = RotatedVarValues(1:j)
        END IF

        IF (Parallel) THEN
           IF ( .NOT. ASSOCIATED(Solver % Matrix % ParMatrix) ) &
                CALL ParallelInitMatrix( Solver, Solver % Matrix )

           n =  Solver % Matrix % NumberOfRows
           IF ( .NOT. ALLOCATED(Residual) ) THEN
              ALLOCATE( Residual(n) )
              Residual = 0.0d0
           END IF
           IF ( .NOT. ALLOCATED(TempRes) ) THEN           
              ALLOCATE( TempRes(n) )
              TempRes = 0.0d0
           END IF
           IF ( .NOT. ALLOCATED(TempS) ) THEN           
              ALLOCATE( TempS(n) )
              TempS = 0.0d0
           END IF
           IF ( .NOT. ALLOCATED(TempRHS) ) THEN           
              ALLOCATE( TempRHS(n) )
              TempRHS = 0.0d0
           END IF


           UpdateMatrix = .TRUE.
           CALL ParallelInitSolve( Solver % Matrix, Solver % Variable % Values, &
                Solver % Matrix % RHS, Residual, UpdateMatrix )

           MMatrix => ParallelMatrix( Solver % Matrix, Mx, Mb, Mr )
           n = MMatrix % NumberOfRows

           IF ( .NOT. ALLOCATED(TempMx) ) THEN           
             ALLOCATE( TempMx(n) )
           END IF
           IF ( .NOT. ALLOCATED(TempMb) ) THEN           
             ALLOCATE( TempMb(n) )
           END IF
           IF ( .NOT. ALLOCATED(TempMr) ) THEN           
             ALLOCATE( TempMr(n) )
           END IF

        ELSE
           n =  Solver % Matrix % NumberOfRows
           IF ( .NOT. ALLOCATED(Residual) ) THEN
              ALLOCATE( Residual(n) )
              Residual = 0.0d0
           END IF
           IF ( .NOT. ALLOCATED(TempRes) ) THEN           
              ALLOCATE( TempRes(n) )
              TempRes = 0.0d0
           END IF
           IF ( .NOT. ALLOCATED(TempS) ) THEN           
              ALLOCATE( TempS(n) )
              TempS = 0.0d0
           END IF

           MMatrix => Solver % Matrix
           Mx => Solver % Variable % Values
           Mb => Solver % Matrix % RHS
           Mr => Residual

        END IF

        !------------------------------------------------------------------------------------------
        ! Test whether the nonlinear residual is small enough to terminate the nonlinear iteration
        !------------------------------------------------------------------------------------------
        IF (Iter >= 1) THEN
           !--------------------------------------------------------------------------------------
           ! The rotated variables are used to compute the current residual.
           !--------------------------------------------------------------------------------------
           j = Solver % Matrix % NumberOfRows
           CALL Mymv( Solver % Matrix, Solver % Variable % Values, Residual, .TRUE.)
           Residual(1:j) = Solver % Matrix % RHS(1:j) - Residual(1:j)

           IF ( Parallel ) THEN
             DO j=1,SIZE(Mr)
               Mr(j) = Mb(j) - Mr(j)
             END DO
           END IF
           NonLinError = Mynorm( n, Mr )/Mynorm( n, Mb )
           Solver % Variable % NonLinChange = NonLinError

           WRITE(Message,'(a,I4,ES12.3)') 'Residual for nonlinear iterate', &
              Iter-1, NonLinError
           CALL Info('StokesSolver', Message, Level=3)            

           IF ( NonLinError < NonlinearTol .OR. Iter==NonlinearIter ) THEN
             DEALLOCATE( Residual )
             DEALLOCATE( TempRes )
             DEALLOCATE( TempS )
             IF (Parallel) THEN
               DEALLOCATE( TempMx )
               DEALLOCATE( TempMr )
               DEALLOCATE( TempMb )             
               DEALLOCATE( TempRHS )
             END IF
             CALL BackRotateNTSystem( Solver % Variable % Values, Solver % Variable % Perm, &
                  Solver % Variable % DOFs )
             EXIT
           END IF
        !---------------------------------------------------------------------
        END IF ! ...Nonlinear convergence check
        !--------------------------------------------------------------------

        !------------------------------------------------------------------
        ! Initialize the variable containing the new search direction...
        !------------------------------------------------------------------
        Snew = 0.0d0 

        IF (AdaptiveTols) THEN
           !---------------------------------------------------------------------------------------
           ! Adapt the linear system convergence tolerance in terms of the current nonlinear error
           !---------------------------------------------------------------------------------------
           IF (Iter >= 1 ) THEN
              MinTolerance = RelTolerance * NonLinError
              IF (MinTolerance > BaseTolerance) MinTolerance = BaseTolerance
              IF (MinTolerance < TargetTol) MinTolerance = TargetTol
           ELSE
              MinTolerance = BaseTolerance
           END IF
        ELSE
           MinTolerance = ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance')
        END IF

        UseTrueResidual = .TRUE.
        !----------------------------------------------------------------------
        ! Perform solution timing for the actual GCR loop if desired...
        !--------------------------------------------------------------------- 
        Timing = ListGetLogical(Solver % Values,'Linear System Timing', Found)

        IF( Timing ) THEN
          t0 = CPUTime(); rt0 = RealTime()
        END IF

        !---------------------------------------------------------------------------
        ! The actual GCR iteration loop...
        !---------------------------------------------------------------------------
        DO Round=1, MaxIterations
           !-----------------------------------------------------------------------
           ! Generate a new search direction by applying the preconditioner. That is,
           ! find an approximation to the error by solving Pe=r, with P the
           ! preconditioner and r=b-Ax the current residual of the primary system.
           ! In principle, we can employ the true residual obtained by performing
           ! the matrix-vector multiplication or utilize the residual vector
           ! generated by the iterative algorithm. The iterative residual cannot
           ! be employed currently, since the residual entries corresponding to nodes
           ! which are not owned by the partition are not received in parallel 
           ! computations...
           !-----------------------------------------------------------------------
           IF ( UseTrueResidual ) THEN
              !-------------------------------------------------------------------------
              ! Compute the current true residual...
              !-------------------------------------------------------------------------
              CALL Mymv( Solver % Matrix, Solver % Variable % Values, Residual, .TRUE. )
              j = Solver % Matrix % NumberOfRows
              Residual(1:j) = Solver % Matrix % RHS(1:j) - Residual(1:j)
           END IF

           !----------------------------------------------------------------------------
           ! Update the preconditioner system rhs by assuming that the preconditioning
           ! system and the primary system have the same permutation...
           !--------------------------------------------------------------------------
           ! we multiply by the inverse of the row scaling of the global matrix, for  
           ! symmetric scaling we have set the arrays to 1 so nothing happens here.   
           ! For row equilibration this makes sure that the correct system is solved  
           ! without having to scale the A and Q matrices.
           j = PressureSolver % Matrix % NumberOfRows
           DO i=1,j
              PressureSolver % Matrix % RHS(i) = Residual((dim+1)*i) * RowScalingP(i)
              DO p=1,dim
                 VelocitySolver % Matrix % RHS( dim*(i-1)+p ) = Residual( (dim+1)*(i-1)+p )  &
                        * RowScalingA( dim*(i-1)+p )
              END DO
           END DO

           !--------------------------------------------------------------------------
           ! Solve the preconditioning systems sequentially...
           !--------------------------------------------------------------------------
           IF (Round > 1) THEN 
             CALL ListAddLogical(PressureSolver % Values, &
               'No Precondition Recompute', .TRUE.)
             CALL ListAddLogical(PressureSolver % Values, &
               'Linear System Refactorize', .FALSE.)
           END IF
           CurrentModel % Solver => PressureSolver
           CALL Info( 'Pressure Preconditioning', ' ', Level=4 )
           CALL Info( 'Pressure Preconditioning', '-------------------------------------', Level=4 )
           CALL Info( 'Pressure Preconditioning', 'Performing linear solve', Level=4 )
           CALL Info( 'Pressure Preconditioning', '-------------------------------------', Level=4 )
           CALL Info( 'Pressure Preconditioning', ' ', Level=4 )
           Norm = DefaultSolve(PressureSolver)
           CurrentModel % Solver => Solver
           CALL ListAddLogical(PressureSolver % Values, &
               'No Precondition Recompute', .FALSE.)
           CALL ListAddLogical(PressureSolver % Values, &
               'Linear System Refactorize', .TRUE.)               
           
           !-----------------------------------------------------------------------------
           ! Perform matrix-vector product to produce upper triangular preconditioner...
           !-----------------------------------------------------------------------------
           IF (.NOT. Parallel) THEN
             TempS = 0.0d0
             j = Solver % Matrix % NumberOfRows/(dim+1)
             DO i=1,j
               TempS( i*(dim+1) ) = PressureSolver % Variable % Values(i)
             END DO
             CALL Mymv( Solver % Matrix, TempS, TempRes )
             DO i=1,j
               DO p=1,dim
                 VelocitySolver % Matrix % RHS( dim*(i-1)+p ) = &
                     VelocitySolver % Matrix % RHS( dim*(i-1)+p ) - &
                     TempRes( (dim+1)*(i-1)+p ) * RowScalingA( dim*(i-1)+p )
               END DO
             END DO
           ELSE
             !--------------------------------------------------------
             ! Save the current status of parallel matrix structure
             !-------------------------------------------------------
             j = MMatrix % NumberOfRows
             TempMx(1:j) = Mx(1:j)
             TempMb(1:j) = Mb(1:j)
             TempMr(1:j) = Mr(1:j)
             j=Solver % Matrix % NumberOfRows
             TempS(1:j) = Solver % Variable % Values(1:j)
             TempRHS(1:j) = Solver % Matrix % RHS(1:j)

             TempRes = 0.0d0
             Solver % Matrix % RHS = 0.0d0
             Solver % Variable % Values = 0.0d0
             j = Solver % Matrix % NumberOfRows/(dim+1)
             DO i=1,j
               Solver % Variable % Values( i*(dim+1) ) = PressureSolver % Variable % Values(i)
             END DO

             !-------------------------------------------------------------------------------
             ! Update the vectors of the parallel matrix structure   
             !-------------------------------------------------------------------------------
             CALL ParallelUpdateSolve( Solver % Matrix, Solver % Variable % Values, TempRes )
             !-----------------------------------------------------------------------
             ! Perform mv-product and insert the result into the argument vectors...
             !-----------------------------------------------------------------------
             CALL Mymv( Solver % Matrix, Solver % Variable % Values, TempRes, .TRUE.)
           
             !-------------------------------------------------------------------------
             ! Update the velocity solver rhs by adding the result of the gradient mv...
             !-------------------------------------------------------------------------
             DO i=1,j
               DO p=1,dim
                 VelocitySolver % Matrix % RHS( dim*(i-1)+p ) = &
                     VelocitySolver % Matrix % RHS( dim*(i-1)+p ) - &
                     TempRes( (dim+1)*(i-1)+p ) * RowScalingA(dim*(i-1)+p )
                       
               END DO
             END DO
             !-----------------------------------------------------------------------------
             ! Finally retrieve the original status of the parallel matrix structure....
             !----------------------------------------------------------------------------
             j=Solver % Matrix % NumberOfRows
             Solver % Variable % Values(1:j) = TempS(1:j)
             Solver % Matrix % RHS(1:j) = TempRHS(1:j)
 
             j = MMatrix % NumberOfRows
             Mx(1:j) = TempMx(1:j)
             Mb(1:j) = TempMb(1:j)
             Mr(1:j) = TempMr(1:j)
           END IF

           IF (Round > 1) THEN
             CALL ListAddLogical(VelocitySolver % Values, &
               'No Precondition Recompute', .TRUE.)
             CALL ListAddLogical(VelocitySolver % Values, &
               'Linear System Refactorize', .FALSE.)
           END IF
           CurrentModel % Solver => VelocitySolver
           !----------------------------------------------------------------------------
           ! Adapt the preconditioning system convergence tolerance in terms of the 
           ! current convergence tolerance for the primary system
           !----------------------------------------------------------------------------
           !IF (AdaptiveTols) THEN
           IF (.FALSE.) THEN ! Disabled: Using a fixed tolerance may actually be a better strategy
              PrecondTol = 1/RelTolerance * MinTolerance
              IF (PrecondTol > BaseTolerance) PrecondTol = BaseTolerance
              CALL ListAddConstReal(VelocitySolver % Values, &
                   'Linear System Convergence Tolerance', BaseTolerance)
           END IF
           CALL Info( 'Velocity Preconditioning', ' ', Level=4 )
           CALL Info( 'Velocity Preconditioning', '-------------------------------------', Level=4 )
           CALL Info( 'Velocity Preconditioning', 'Performing linear solve', Level=4 )
           CALL Info( 'Velocity Preconditioning', '-------------------------------------', Level=4 )
           CALL Info( 'Velocity Preconditioning', ' ', Level=4 )
           Norm = DefaultSolve(VelocitySolver)      
           CurrentModel % Solver => Solver
           CALL ListAddLogical(VelocitySolver % Values, &
               'No Precondition Recompute', .FALSE.)  
           CALL ListAddLogical(VelocitySolver % Values, &
               'Linear System Refactorize', .TRUE.)  

           !-------------------------------------------------------------------------
           ! Insert the computed approximation of the error into the search direction 
           ! variable which will be given to the GCR routine...
           !--------------------------------------------------------------------------
           DO t=1, Active
              Element => GetActiveElement(t)
              nd = GetElementDOFs( Indexes )

              DO i=1,nd
                 j = Solver % Variable % Perm( Indexes(i) )
                 k = PressureSolver % Variable % Perm( Indexes(i) )
                 Snew( j*(dim+1) ) = PressureSolver % Variable % Values(k)

                 k = VelocitySolver % Variable % Perm( Indexes(i) )   
                 DO p=1,dim          
                    Snew( (j-1)*(dim+1)+p ) = VelocitySolver % Variable % Values( (j-1)*dim+p )
                 END DO

              END DO
           END DO

           !-----------------------------------------------------------------------------
           ! Perform GCR update for solving the primitive linear system... 
           !-----------------------------------------------------------------------------
           CALL GCRUpdate(n, Solver % Matrix, MMatrix, Mx, Mb, Mr, Snew, S, V, R, Round, &
                Norm, RestartM)

           IF (Parallel) &
                CALL ParallelUpdateResult( Solver % Matrix, Solver % Variable % Values, Residual )

           IF (Norm < MinTolerance) EXIT

        END DO
        
        ! unscale solution if column scaling was used
        IF (DoScaling .AND. (.NOT. DoEquilibration)) THEN
          CALL BackScaleLinearSystem(Solver,Solver%Matrix,Solver%Matrix%RHS,Solver%Variable%Values)
        END IF

        IF( Timing ) THEN
          st  = CPUTime() - t0;
          rst = RealTime() - rt0

          CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys cpu time '&
              //GetVarName(Solver % Variable),st)
          CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys real time '&
              //GetVarName(Solver % Variable),rst)
          WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time (CPU,REAL) for '&
              //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
          CALL Info('SolveSystem',Message)    

          IF( ListGetLogical(Solver % Values,'Linear System Timing Cumulative',Found)) THEN
            ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
                //GetVarName(Solver % Variable),Found)
            st = st + ct
            ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
                //GetVarName(Solver % Variable),Found)
            rst = rst + ct
            CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
                //GetVarName(Solver % Variable),st)
            CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
                //GetVarName(Solver % Variable),rst)
          END IF

        END IF

        DEALLOCATE( Residual )
        DEALLOCATE( TempRes )        
        DEALLOCATE( TempS )          

        ! Compute Var Loads if needed
        IF ( ListGetLogical( Solver % Values,'Calculate Loads', Found ) ) &
                         CALL ComputeVarLoads(Solver)


        n =  Solver % Matrix % NumberOfRows
        RotatedVarValues(1:n) = Solver % Variable % Values(1:n)
        CALL BackRotateNTSystem( Solver % Variable % Values, Solver % Variable % Perm, &
            Solver % Variable % DOFs )


      END IF

   END DO

   IF (BlockPreconditioning) &
       DEALLOCATE(RotatedVarValues)

   IF (P2P1) THEN
     !----------------------------------------------------------------------------------------
     ! Replace the zero pressure solution at the nodes which are not needed in the linear
     ! pressure approximation by the interpolated values for right visualization:
     !----------------------------------------------------------------------------------------
     DO t=1,Active
       ! First the midedge nodes:
       Element => GetActiveElement(t)
       nd = GetElementDOFs( Indexes )
       k = GetElementFamily()
       EdgeMap => GetEdgeMap(k)
       SELECT CASE(k)
       CASE (3)
         MidEdgeNodes(1:3) = (/ 4, 5, 6 /)
       CASE (4)
         MidEdgeNodes(1:4) = (/ 5, 6, 7, 8 /)
       CASE (5)
         MidEdgeNodes(1:6) = (/ 5, 6, 7, 8, 9, 10 /) 
       CASE (6)
         MidEdgeNodes(1:8) = (/ 6, 7, 8, 9, 10, 11, 12, 13 /)
       CASE (7)
         MidEdgeNodes(1:9) = (/ 7, 8, 9, 10, 11, 12, 13, 14, 15 /)
       CASE (8)
         MidEdgeNodes(1:12) = (/ 9, 10, 11, 12, 17, 18, 19, 20, 13, 14, 15, 16 /)
       END SELECT

       DO q=1,SIZE(EdgeMap,1)
         m = (dim+1) * Solver % Variable % Perm(Indexes(MidEdgeNodes(q)))
         i = (dim+1) * Solver % Variable % Perm(Indexes(EdgeMap(q,1)))
         j = (dim+1) * Solver % Variable % Perm(Indexes(EdgeMap(q,2)))
         Solver % Variable % Values(m) = 0.5d0 * ( Solver % Variable % Values(i) + &
             Solver % Variable % Values(j) )   
       END DO

       ! The pressure at the midface nodes for 409 elements: 
       IF (k==4 .AND. nd==9) THEN
         res = 0.0d0
         DO q=1,4
           i = (dim+1) * Solver % Variable % Perm(Indexes(q))
           res = res + Solver % Variable % Values(i)
         END DO
         m = (dim+1) * Solver % Variable % Perm(Indexes(9))
         Solver % Variable % Values(m) = 0.25d0 * res
       END IF

       ! The pressure at the midpoint and at the midface nodes for 827 elements:
       IF (k==8 .AND. nd==27) THEN
         BrickFaceMap(1,:) = (/ 1,2,6,5 /)         
         BrickFaceMap(2,:) = (/ 2,3,7,6 /)
         BrickFaceMap(3,:) = (/ 4,3,7,8 /)
         BrickFaceMap(4,:) = (/ 1,4,8,5 /)
         BrickFaceMap(5,:) = (/ 1,2,3,4 /)          
         BrickFaceMap(6,:) = (/ 5,6,7,8 /)
         DO j=1,6
           res = 0.0d0
           DO q=1,4
             i = (dim+1) * Solver % Variable % Perm(Indexes(BrickFaceMap(j,q)))
             res = res + Solver % Variable % Values(i)
           END DO
           m = (dim+1) * Solver % Variable % Perm(Indexes(20+j))
           Solver % Variable % Values(m) = 0.25d0 * res
         END DO           

         res = 0.0d0
         DO q=1,8
           i = (dim+1) * Solver % Variable % Perm(Indexes(q))
           res = res + Solver % Variable % Values(i)
         END DO
         m = (dim+1) * Solver % Variable % Perm(Indexes(27))      
         Solver % Variable % Values(m) = 0.125d0 * res
       END IF

     END DO
   END IF

CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE GCRUpdate( n, A, M, x, b, r, Snew, S, V, RR, Round, res, MParam)
!------------------------------------------------------------------------------
    INTEGER :: n, Round, MParam
    TYPE(Matrix_t), POINTER :: A, M
    REAL(KIND=dp) :: x(:), b(:), r(:), Snew(:)  
    REAL(KIND=dp) :: S(:,:), V(:,:), RR(:), res
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: T1(n), T2(n), beta, alpha
    INTEGER :: i,j,k
!--------------------------------------------------------------------------------

    IF ( Parallel ) CALL ParallelVector(A,Snew)   

    !----------------------------------------------
    ! Check for restarting
    !--------------------------------------------- 
    IF ( MOD(Round,MParam)==0 ) THEN
       j = MParam
    ELSE
       j = MOD(Round,MParam)
    END IF

    IF ( j == 1) THEN
      CALL Mymv( A, x, r ) 
      r(1:n) = b(1:n)-r(1:n)
      res = MyNorm(n,r)/Mynorm(n,b)
      IF (Round==1) THEN
         WRITE(Message,'(a,I4,ES12.3)') 'GCR residual for iterate', &
              0, res
         CALL Info('Outer Iteration',Message,Level=4)
      END IF
     ELSE
      r(1:n) = RR(1:n)
    END IF

    T1(1:n) = Snew(1:n)
    CALL Mymv( A, T1, T2 )
    !--------------------------------------------------------------
    ! Perform the orthogonalisation of the search directions...
    !--------------------------------------------------------------
    DO i=1,j-1
       beta = Mydot( n, V(1:n,i), T2(1:n) )
       T1(1:n) = T1(1:n) - beta * S(1:n,i)
       T2(1:n) = T2(1:n) - beta * V(1:n,i)    
    END DO

    alpha = Mynorm(n,T2)
    T1(1:n) = 1.0d0/alpha * T1(1:n)
    T2(1:n) = 1.0d0/alpha * T2(1:n)

    !-------------------------------------------------------------
    ! The update of the solution and save the search data...
    !------------------------------------------------------------- 
    beta = Mydot(n, T2, r)

    x(1:n) = x(1:n) + beta * T1(1:n)
    r(1:n) = r(1:n) - beta * T2(1:n)

    IF ( j /= MParam ) THEN
       S(1:n,j) = T1(1:n)
       V(1:n,j) = T2(1:n)
    END IF

    RR(1:n) = r(1:n)
    res = MyNorm(n,r)/Mynorm(n,b)

    WRITE(Message,'(a,I4,ES12.3)') 'GCR residual for iterate', &
        Round, res
    CALL Info('Outer Iteration',Message,Level=4)

!-------------------------------------------
  END SUBROUTINE GCRUpdate
!---------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, VeloBlock, Mass, FORCE, LOAD, Nodalrho, &
       Nodalmu, Vx, Vy, Vz, Element, n, nd, dim, Convect, P2P1, SkipPowerLaw, &
       Newton, DiagonalA)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), VeloBlock(:,:), Mass(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), Vx(:), Vy(:), Vz(:), GradDivParam
    INTEGER :: dim, n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Stabilization, Convect, P2P1, SkipPowerLaw, Newton, DiagonalA
    !------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: PressureElement => NULL()
    REAL(KIND=dp) :: Basis(nd), dBasisdx(nd,3), &
         DetJ, LoadAtIP(dim+1), Velo(dim), Grad(dim,dim), AK, w1, w2, w3, ViscAtIp, RhoAtIP
    REAL(KIND=dp) :: LinBasis(nd)
    REAL(KIND=dp), POINTER :: A(:,:), F(:), M(:,:), Jac(:,:)
    REAL(KIND=dp), TARGET :: JacM(nd*(dim+1),nd*(dim+1)), Sol(nd*(dim+1)) 
    LOGICAL :: Stat, ViscNewtonLin
    INTEGER :: t, i, j, k, l, p, q
    INTEGER :: Code, nlin
    INTEGER :: LinearCode(3:8) = (/ 303,404,504,605,706,808 /)
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c, c2, ch, rotterm, mK, hK, Re, VNorm, &
         a1, a2, muder, muder0, Strain(dim,dim)

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes, PressureElement
    !------------------------------------------------------------------------------
    
    IF (P2P1) THEN
      IF ( .NOT. ASSOCIATED(PressureElement) ) PressureElement => AllocateElement()
      k = GetElementFamily()
      PressureElement % Type => GetElementType( LinearCode(k), .FALSE. )
      nlin = PressureElement % Type % NumberOfNodes
    END IF

    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0
    Mass = 0.0d0
    VeloBlock = 0.0d0
    IF (Newton) JacM = 0.0d0
    ViscNewtonLin = .FALSE.

    !----------------------
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

       !--------------------------------------------------------------
       ! Linear pressure basis functions for P2-P1 approximation:
       !--------------------------------------------------------------
       IF (P2P1) CALL NodalBasisFunctions( nlin, LinBasis, PressureElement, &
           IP % U(t), IP % V(t), IP % W(t) )

       s = IP % s(t) * detJ

       Grad = 0.0d0
       DO i=1,dim
          Grad(1,i) = SUM( Vx(1:n) * dBasisdx(1:n,i) )
          Grad(2,i) = SUM( Vy(1:n) * dBasisdx(1:n,i) )
          IF ( DIM > 2 ) Grad(3,i) = SUM( Vz(1:n) * dBasisdx(1:n,i) )
       END DO

       !----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       ViscAtIP  = SUM( Basis(1:n) * Nodalmu(1:n) )
       RhoAtIP = SUM( Basis(1:n) * Nodalrho(1:n) )
       
       IF ( SkipPowerLaw ) THEN
          mu = ViscAtIP
       ELSE
          IF( ListCheckPresent( Material, 'Viscosity Model' ) ) THEN
             mu = EffectiveViscosity( ViscAtIP, RhoAtIp, Vx, Vy, Vz, &
                  Element, Nodes, n, n, IP % U(t), IP % V(t), &
                  IP % W(t), muder0 )
             ViscNewtonLin = Newton .AND. muder0/= 0.0d0
             IF ( ViscNewtonLin )  Strain = (Grad+TRANSPOSE(Grad))/2
          ELSE
             mu = ViscAtIP
          END IF
       END IF

       IF (Convect) THEN
          w1 = SUM( Vx(1:n) * Basis(1:n) ) 
          w2 = SUM( Vy(1:n) * Basis(1:n) )           
          IF (dim > 2) w3 = SUM( Vz(1:n) * Basis(1:n) )  
       END IF

       !--------------------------------------------
       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:dim+1) = MATMUL( Basis(1:n), LOAD(1:n,1:dim+1) )

       !----------------------------------------------------------------------------------
       ! The system matrix with linear pressure approximation
       !---------------------------------------------------------------------------------
       DO p=1,nd
          DO q=1,nd
             i = (dim+1) * (p-1) + 1
             j = (dim+1) * (q-1) + 1
             A => STIFF(i:i+dim,j:j+dim)
             IF ( ViscNewtonLin ) Jac => JacM( i:i+dim,j:j+dim )

             IF ( ViscNewtonLin ) THEN
               DO i=1,dim
                 muder = 4.0d0*muder0*SUM(Strain(i,1:dim)*dBasisdx(q,1:dim))
                 DO j=1,dim
                   a1 = 0.0d0
                   a2 = 0.0d0
                   DO k=1,dim
                     a1 = a1 + Grad(j,k)*dBasisdx(p,k)
                     a2 = a2 + Grad(k,j)*dBasisdx(p,k)
                   END DO
                   Jac(j,i) = Jac(j,i) + s * muder * (a1+a2)
                   !  Jac(i,i)=Jac(i,i) + s*muder*Grad(i,j)*dBasisdx(p,j)
                   !   Jac(j,i)=Jac(j,i) + s*muder*Grad(i,j)*dBasisdx(p,i)
                 END DO
               END DO
             END IF

             DO i=1,dim
               DO j = 1,dim
                 A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                 A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
               END DO

               ! Only linear pressure approximation is used...
               IF (P2P1) THEN
                 IF (q <= nlin) &
                     A(i,dim+1) = A(i,dim+1) - s * LinBasis(q) * dBasisdx(p,i)
                 IF (p <= nlin) &
                     A(dim+1,i) = A(dim+1,i) - s * dBasisdx(q,i) * LinBasis(p)
               ELSE
                 IF (q <= n) &
                     A(i,dim+1) = A(i,dim+1) - s * Basis(q) * dBasisdx(p,i)

                 IF (p <= n) &   
                     A(dim+1,i) = A(dim+1,i) - s * dBasisdx(q,i) * Basis(p)
               END IF
             END DO

             IF (DiagonalA) THEN
               DO i=1,dim
                 DO j = 1,dim
                   
                   VeloBlock( dim*(p-1)+i, dim*(q-1)+i ) = VeloBlock( dim*(p-1)+i, dim*(q-1)+i ) + &
                       s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                   
                 END DO
               END DO
             END IF

             IF (Convect) THEN
                A(1,1) = A(1,1) + s * RhoAtIP * w1 * dBasisdx(q,1) * Basis(p)
                A(1,1) = A(1,1) + s * RhoAtIP * w2 * dBasisdx(q,2) * Basis(p)
                A(2,2) = A(2,2) + s * RhoAtIP * w1 * dBasisdx(q,1) * Basis(p)
                A(2,2) = A(2,2) + s * RhoAtIP * w2 * dBasisdx(q,2) * Basis(p)

                IF (dim > 2) THEN
                   A(1,1) = A(1,1) + s * RhoAtIP * w3 * dBasisdx(q,3) * Basis(p)
                   A(2,2) = A(2,2) + s * RhoAtIP * w3 * dBasisdx(q,3) * Basis(p)
                   A(3,3) = A(3,3) + s * RhoAtIP * w1 * dBasisdx(q,1) * Basis(p)
                   A(3,3) = A(3,3) + s * RhoAtIP * w2 * dBasisdx(q,2) * Basis(p)
                   A(3,3) = A(3,3) + s * RhoAtIP * w3 * dBasisdx(q,3) * Basis(p)
                END IF

             END IF
             

          END DO

          i = (dim+1) * (p-1) + 1
          F => FORCE(i:i+dim)
          F = F + s * RhoAtIP * LoadAtIP * Basis(p)

       END DO

       IF (P2P1) THEN
         DO p=1,nlin
           DO q=1,nlin       
             Mass(p,q) = Mass(p,q) - s * 1.0d0/mu * Basis(p) * Basis(q)            
           END DO
         END DO
       ELSE
         DO p=1,n
           DO q=1,n       
             Mass(p,q) = Mass(p,q) - s * 1.0d0/mu * Basis(p) * Basis(q)            
           END DO
         END DO
       END IF
    END DO

    IF ( ViscNewtonLin ) THEN
       SOL=0.0d0
       SOL(1:(dim+1)*n:dim+1) = Vx(1:n)
       SOL(2:(dim+1)*n:dim+1) = Vy(1:n)
       IF ( dim>2 ) SOL(3:(dim+1)*n:dim+1) = Vz(1:n)
       p = (dim+1)*n
       Stiff(1:p,1:p) = Stiff(1:p,1:p)+JacM(1:p,1:p)
       FORCE(1:p)=FORCE(1:p)+MATMUL(JacM(1:p,1:p),SOL(1:p))
    END IF

    ! Omit dofs that do not correspond to linear pressure approximation...
    IF (P2P1) THEN
      DO p = nlin+1,nd
        i = (dim+1) * p
        FORCE(i)   = 0.0d0
        STIFF(i,:) = 0.0d0
        STIFF(:,i) = 0.0d0       
        STIFF(i,i) = 1.0d0
        Mass(p,p) = 1.0d0
      END DO
    ELSE
      DO p = n+1,nd
        i = (dim+1) * p
        FORCE(i)   = 0.0d0
        STIFF(i,:) = 0.0d0
        STIFF(:,i) = 0.0d0       
        STIFF(i,i) = 1.0d0
        Mass(p,p) = 1.0d0
      END DO
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE StaticCondensation( N, Nb, dim, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb, dim
    REAL(KIND=dp) :: K(:,:), F(:), Kbb(nb*dim,nb*dim), &
         Kbl(nb*dim,n*(dim+1)),Klb(n*(dim+1),nb*dim),Fb(nb*dim)

    INTEGER :: m, i, j, l, p, Cdofs((dim+1)*n), Bdofs(dim*nb)

    m = 0
    DO p = 1,n
      DO i = 1,dim+1
        m = m + 1
        Cdofs(m) = (dim+1)*(p-1) + i
      END DO
    END DO

    m = 0
    DO p = 1,nb
      DO i = 1,dim
        m = m + 1
        Bdofs(m) = (dim+1)*(p-1) + i + n*(dim+1)
      END DO
    END DO

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Cdofs)
    Klb = K(Cdofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb*dim )

    F(1:(dim+1)*n) = F(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    K(1:(dim+1)*n,1:(dim+1)*n) = &
    K(1:(dim+1)*n,1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )

!------------------------------------------------------------------------------
  END SUBROUTINE StaticCondensation
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
      FUNCTION Mydot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp) :: s,x(:),y(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = ParallelDot( n, x, y )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION Mydot
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
    SUBROUTINE Mymv( A, x, b, Update )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       IF ( .NOT. Parallel ) THEN
         CALL CRS_MatrixVectorMultiply( A, x, b )
       ELSE
         IF ( PRESENT( Update ) ) THEN
           CALL ParallelMatrixVector( A,x,b,Update )
         ELSE
           CALL ParallelMatrixVector( A,x,b )
         END IF
       END IF
!------------------------------------------------------------------------------
     END SUBROUTINE Mymv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      FUNCTION Mynorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp) :: s,x(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
        ELSE
          s = ParallelNorm( n, x )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION Mynorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE SetBoundaryConditions( Model, StiffMatrix, &
                   Name, DOF, NDOFs, Perm, rhs )
!------------------------------------------------------------------------------
!
! Set dirichlet boundary condition for given dof
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
!******************************************************************************
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

      !------------------------------------------------------------------------------
      ! If standard element types are used, there is no need for these modifications
      !------------------------------------------------------------------------------
      IF (n==nd) CYCLE

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
 SUBROUTINE NavierStokesBoundary( BoundaryMatrix, BoundaryVector, LoadVector,   &
     NodalExtPressure, NodalSlipCoeff, NormalTangential, Element, n, nd, Nodes )
             
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Return element local matrices and RSH vector for Navier-Stokes-equations
!  boundary conditions.
!
!  ARGUMENTS:
!
!  REAL(KIND=dp) :: BoundaryMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: BoundaryVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:,:)
!     INPUT: Nodal values force in coordinate directions
!
!  REAL(KIND=dp) :: NodalAlpha(:,:)
!     INPUT: Nodal values of force in normal direction
!
!  REAL(KIND=dp) :: NodalBeta(:,:)
!     INPUT: Nodal values of something which will be taken derivative in
!            tangential direction and added to force...
!
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of boundary element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!******************************************************************************
!------------------------------------------------------------------------------
   USE ElementUtils

   IMPLICIT NONE

   REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), LoadVector(:,:), &
       NodalSlipCoeff(:,:), NodalExtPressure(:)

   INTEGER :: n, nd, pn

   TYPE(Element_t),POINTER  :: Element, Parent
   TYPE(Nodes_t) :: Nodes, ParentNodes

   LOGICAL :: NormalTangential

!------------------------------------------------------------------------------
!  Local variables
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
   REAL(KIND=dp) :: detJ,FlowStress(3,3),SlipCoeff

   REAL(KIND=dp) :: u,v,w,ParentU,ParentV,ParentW,s,x(n),y(n),z(n)
   REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)
   REAL(KIND=dp) :: TangentForce(3),Force(3),Normal(3),Tangent(3),Tangent2(3), &
               Vect(3), Alpha, mu,Grad(3,3),Velo(3)

   REAL(KIND=dp) :: xx, yy, ydot, ydotdot, MassFlux

   INTEGER :: i,j,k,l,k1,k2,t,q,p,c,dim,N_Integ

   LOGICAL :: stat, Found

   TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------
!   NormalTangential = .TRUE.
   dim = CoordinateSystemDimension()
   c = dim + 1
!
!------------------------------------------------------------------------------
!  Integration stuff
!------------------------------------------------------------------------------
   IntegStuff = GaussPoints( element )
   U_Integ => IntegStuff % u
   V_Integ => IntegStuff % v
   W_Integ => IntegStuff % w
   S_Integ => IntegStuff % s
   N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!  Now we start integrating
!------------------------------------------------------------------------------
   DO t=1,N_Integ

     u = U_Integ(t)
     v = V_Integ(t)
     w = W_Integ(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                 Basis, dBasisdx )

     s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!    Add to load: tangetial derivative of something
!------------------------------------------------------------------------------
     DO i=1,dim
       TangentForce(i) = 0.0d0
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in coordinate directions
!------------------------------------------------------------------------------
     Force = 0.0d0
     DO i=1,dim
       Force(i) = Force(i) + SUM( LoadVector(i,1:n)*Basis )
     END DO

!------------------------------------------------------------------------------
!    Add to load: given force in normal direction
!------------------------------------------------------------------------------
     Normal = NormalVector( Element, Nodes, u,v,.TRUE. )

     Alpha = SUM( NodalExtPressure(1:n) * Basis )
     IF ( NormalTangential ) THEN
       Force(1) = Force(1) + Alpha
     ELSE
        DO i=1,dim
           Force(i) = Force(i) + Alpha * Normal(i)
        END DO
     END IF

!------------------------------------------------------------------------------

     Alpha = 0.0d0
     MassFlux = 0.0d0

!------------------------------------------------------------------------------

     SELECT CASE( Element % TYPE % DIMENSION )
     CASE(1)
        Tangent(1) =  Normal(2)
        Tangent(2) = -Normal(1)
        Tangent(3) =  0.0_dp
        Tangent2   =  0.0_dp
     CASE(2)
        CALL TangentDirections( Normal, Tangent, Tangent2 ) 
     END SELECT

     IF ( ANY( NodalSlipCoeff(:,:) /= 0.0d0 ) ) THEN
       DO p=1,nd
         DO q=1,nd
           DO i=1,dim
             SlipCoeff = SUM( NodalSlipCoeff(i,1:n) * Basis(1:n) )

             IF ( NormalTangential ) THEN
                SELECT CASE(i)
                   CASE(1)
                     Vect = Normal
                   CASE(2)
                     Vect = Tangent
                   CASE(3)
                     Vect = Tangent2
                END SELECT

                DO j=1,dim
                   DO k=1,dim
                      BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) = &
                         BoundaryMatrix( (p-1)*c+j,(q-1)*c+k ) + &
                          s * SlipCoeff * Basis(q) * Basis(p) * Vect(j) * Vect(k)
                   END DO
                END DO
             ELSE
                 BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) = &
                     BoundaryMatrix( (p-1)*c+i,(q-1)*c+i ) + &
                          s * SlipCoeff * Basis(q) * Basis(p)
             END IF
           END DO
         END DO
       END DO
     END IF

     DO q=1,nd
       DO i=1,dim
         k = (q-1)*c + i
         IF ( NormalTangential ) THEN
            SELECT CASE(i)
               CASE(1)
                 Vect = Normal
               CASE(2)
                 Vect = Tangent
               CASE(3)
                 Vect = Tangent2
            END SELECT

            DO j=1,dim
               l = (q-1)*c + j
               BoundaryVector(l) = BoundaryVector(l) + &
                 s * Basis(q) * Force(i) * Vect(j)
            END DO
         ELSE
            BoundaryVector(k) = BoundaryVector(k) + s*Basis(q)*Force(i)
         END IF
         BoundaryVector(k) = BoundaryVector(k) - s * Alpha * dBasisdx(q,i)
         BoundaryVector(k) = BoundaryVector(k) + s * TangentForce(i)*Basis(q)
       END DO
       BoundaryVector(q*c) = BoundaryVector(q*c) + s * MassFlux * Basis(q)
     END DO

   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE NavierStokesBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE ComputeVarLoads(Solver)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: Aaid, Projector
    TYPE(Variable_t), POINTER ::  NodalLoads

    REAL(KIND=dp), POINTER :: SaveValues(:)
    REAL(KIND=dp), ALLOCATABLE :: x(:),TempVector(:), TempRHS(:)

    INTEGER :: DOFs
    INTEGER :: i,ii,j,l,DOF,This

!------------------------------------------------------------------------------


    DOFs=Solver % Variable % DOFs

    NodalLoads => VariableGet( Solver % Mesh % Variables, &
       GetVarName(Solver % Variable) // ' Loads' )

    Aaid => Solver % Matrix

    IF ( ASSOCIATED(NodalLoads) .AND. ASSOCIATED(Aaid % BulkValues) ) THEN
      ALLOCATE(x(SIZE(Solver % Variable % Values)))
      x(:)=Solver % Variable % Values(:)
      ALLOCATE( TempVector(Aaid % NumberOfRows) )

      SaveValues => Aaid % Values
      Aaid % Values => Aaid % BulkValues

      IF ( ParEnv % PEs > 1 ) THEN
        ALLOCATE(TempRHS(SIZE(AAid % BulkRHS)))
        TempRHS = Aaid % BulkRHS
        CALL ParallelInitSolve( Aaid, x, TempRHS, Tempvector )
        CALL ParallelMatrixVector( Aaid, x, TempVector, .TRUE. )
      ELSE
        CALL MatrixVectorMultiply( Aaid, x, TempVector )
      END IF


      Aaid % Values => SaveValues
      IF ( ParEnv % PEs>1 ) THEN
        DO i=1,Aaid % NumberOfRows
          IF ( AAid % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % Mype ) THEN
            TempVector(i) = TempVector(i)-TempRHS(i)
          ELSE
            TempVector(i) = 0
          END IF
        END DO
        DEALLOCATE(TempRHS)
        CALL ParallelSumVector( AAid, Tempvector )
      ELSE
        TempVector = TempVector - Aaid % BulkRHS
      END IF

      DO This=1,CurrentModel % NumberOfBCs
        Projector=>CurrentModel  % BCs(This) % PMatrix
        IF (ASSOCIATED(Projector))THEN
          DO DOF=1,DOFs
            DO i=1,Projector % NumberOfRows
              ii = Projector % InvPerm(i)
              k = Solver % Variable % Perm(ii)
              IF(k<=0) CYCLE
              k = DOFs * (k-1) + DOF
              TempVector(k)=0

              DO l = Projector % Rows(i), Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = Solver % Variable % Perm( Projector % Cols(l) )
                IF ( m > 0 ) THEN
                  m = DOFs * (m-1) + DOF
                  TempVector(k) = TempVector(k) + Projector % Values(l)*TempVector(m)
                 END IF
              END DO
            END DO
          END DO
        END IF
      END DO

      DO i=1,SIZE( NodalLoads % Perm )
        IF ( NodalLoads % Perm(i)>0 .AND. Solver % Variable % Perm(i)>0 ) THEN
           DO j=1,DOFs
             NodalLoads % Values(DOFs*(NodalLoads % Perm(i)-1)+j) =  &
                TempVector(DOFs*(Solver % Variable % Perm(i)-1)+j)
           END DO
         END IF
      END DO
      DEALLOCATE( x,TempVector )

      CALL BackRotateNTSystem(NodalLoads % Values,NodalLoads % Perm,DOFs)

    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE ComputeVarLoads

!------------------------------------------------------------------------------

END SUBROUTINE StokesSolver
!------------------------------------------------------------------------------
