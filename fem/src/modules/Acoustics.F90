!/*****************************************************************************
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
! ****************************************************************************/
!
!/*****************************************************************************
! *
! *****************************************************************************
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           P.O. Box 405
! *           FI-02101 Espoo, Finland 
! *
! ****************************************************************************/



SUBROUTINE AcousticsSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  INTEGER :: dim

  CALL Info('AcousticsSolver','Initialization the solver')
  
  Params => GetSolverParams()
  CALL ListAddNewLogical( Params,'Linear System Complex',.TRUE.)

  dim = CoordinateSystemDimension() 

  IF( ListCheckPresent( Params,'Variable' ) ) THEN
    CALL Warn('AcousticsSolver','Redefining variable name from the given one!')
  END IF

  ! Leave for now since the strings are too short
  IF(.FALSE.) THEN
    IF( dim == 2 ) THEN
      CALL ListAddString( Params,'Variable',&
          'Flow[Re Velocity 1:1 Im Velocity 1:1 Re Velocity 2:1 Im Velocity 2:1 '&
          //' Re Temperature:1 Im Temperature:1 Re Pressure:1 Im Pressure]')
    ELSE
      CALL ListAddString( Params,'Variable',&
          'Flow[Re Velocity 1:1 Im Velocity 1:1 Re Velocity 2:1 Im Velocity 2:1 Re Velocity 3:1 Im Velocity 3:1 '&
          //' Re Temperature:1 Im Temperature:1 Re Pressure:1 Im Pressure:1]')
    END IF
  END IF
    
  
END SUBROUTINE AcousticsSolver_init
  
  

!-----------------------------------------------------------------------------
!>  Solve the time-harmonic, generalized NS-equations assuming ideal gas law.
!------------------------------------------------------------------------------
SUBROUTINE AcousticsSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model            !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt               !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation    !< Steady state or transient simulation
  !------------------------------------------------------------------------------
  ! Local variables related to the block preconditioning...
  !------------------------------------------------------------------------------
  TYPE(Matrix_t),POINTER  :: SMatrix, AMatrix, MMatrix, LMatrix
  LOGICAL :: BlockPreconditioning = .FALSE., OptimizeBW, VelocityAssembly = .TRUE.
  CHARACTER(LEN=MAX_NAME_LEN) :: OuterIterationMethod, IterationMethod
  INTEGER :: MaxIterations, dim
  REAL(KIND=dp) :: Tolerance, ToleranceRatio
  REAL(KIND=dp), ALLOCATABLE :: SLocal(:,:), ALocal(:,:), AnLocal(:,:)
  SAVE SMatrix, AMatrix, MMatrix, LMatrix, SLocal, ALocal, AnLocal
  !-----------------------------------------------------------------------------
  ! Local variables for performing analyses with acoustic interfaces
  !-----------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: HSol, tmp
  LOGICAL :: AcousticInterface, ImpedanceBCIteration
  REAL(KIND=dp), POINTER :: H(:)
  INTEGER, POINTER ::  HPerm(:)
  CHARACTER(3) :: str
  CHARACTER(30) :: str2
  CHARACTER(LEN=MAX_NAME_LEN) :: ElmerStatusFile, CompaStatusFile, CompaMeshPrefix, DataFile
  INTEGER :: MaxCoupledIterations, SolverCalls = 0
  REAL(KIND=dp) :: CoupledTolerance
  SAVE SolverCalls 
  !------------------------------------------------------------------------------
  ! Other local variables
  !------------------------------------------------------------------------------
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(Nodes_t) :: ElementNodes, ParentNodes
  TYPE(Element_t),POINTER :: CurrentElement, Parent
  TYPE(ValueList_t), POINTER :: Material

  INTEGER, POINTER :: NodeIndexes(:), FlowPerm(:)
  INTEGER :: i, j, k, m, n, t, istat, LocalNodes, Dofs, VelocityComponents, &
      VelocityDofs, CoordSys, pn, AcousticI, MaxNodesOnBoundary, np, nlen, nb, RelOrder
  INTEGER, ALLOCATABLE :: Bndries(:), BemElementIndeces(:), AcousticInterfaceNodes(:), &
      NodesOnBoundary(:)

  LOGICAL :: AllocationsDone = .FALSE., Bubbles, MiniBubbles, GotIt, GotIt2, stat, &
      VisitedNodes(Model % NumberOfNodes), SlipBoundary, UtilizePreviousSolution, &
      FirstVisit = .TRUE., BEMCoupling, BEMBoundary, BEMNodesCreated = .FALSE., &
      NodeOnBoundary(10), CalculateMoment, TimeEstimated

  LOGICAL, ALLOCATABLE :: MomentOutput(:)
  SAVE MomentOutput

  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName, VariableName, BoundaryName
!  CHARACTER(LEN=MAX_NAME_LEN) :: VersionID = "$Id: Acoustics.f90,v 1.4 2006/12/12 10:32:03 mmalinen Exp $"

  COMPLEX(KIND=dp) :: A, AverVel, ZAppr, ZAppr2, C1, C2, C3, C4

  REAL(KIND=dp), POINTER :: Flow(:), ForceVector(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: Norm, PrevNorm = 0.0d0, RelativeChange, AngularFrequency, &
      SlipCoefficient1, SlipCoefficient2, SlipCoefficient3, &
      Rea, Ima, Reb, Imb, Rec, Imc, ReP, ImP, &
      NonlinearTol, at, at0, totat, st, totst, &
      MomentAbout(3), Moment(6), Traction(6), Area
#else
  REAL(KIND=dp) :: Norm, PrevNorm = 0.0d0, RelativeChange, AngularFrequency, &
      SlipCoefficient1, SlipCoefficient2, SlipCoefficient3, &
      Rea, Ima, Reb, Imb, Rec, Imc, ReP, ImP, &
      NonlinearTol, at, at0, totat, st, totst, CPUTime, RealTime, &
      MomentAbout(3), Moment(6), Traction(6), Area
#endif
  REAL(KIND=dp), ALLOCATABLE :: LocalStiffMatrix(:,:), Load(:,:), &
      HeatSource(:,:), temp(:), &
      LocalForce(:), SpecificHeat(:), HeatRatio(:), &
      Density(:), Pressure(:), Temperature(:), Conductivity(:), &
      Viscosity(:), Lambda(:), BulkViscosity(:), Impedance(:,:), &
      WallTemperature(:), WallVelocity(:,:), AcImpedances(:,:), &
      AcousticInterfaceResults(:,:), TotalForce(:,:), TotalMoment(:,:), TotalArea(:)
  INTEGER, ALLOCATABLE :: GapIndexes(:)

  SAVE BemElementIndeces, AcousticInterfaceNodes, AcousticInterfaceResults, NodesOnBoundary
  SAVE LocalStiffMatrix, temp, Load, HeatSource,  LocalForce, ElementNodes, ParentNodes, &
       SpecificHeat, HeatRatio, Density, Pressure, &
       Temperature, Conductivity, Viscosity, Lambda, BulkViscosity, &
       Impedance, WallTemperature, WallVelocity, AllocationsDone, FirstVisit, &
       BEMNodesCreated, PrevNorm, TotalForce, TotalMoment, TotalArea, &
       GapIndexes

  TYPE(Variable_t), POINTER :: TimeVar
  LOGICAL :: NewTimeStep, ScanningUsed, FirstTimeStep
  INTEGER :: CurrentDoneTime = 0
  CHARACTER(LEN=MAX_NAME_LEN) :: SimulationType
  REAL(KIND=dp) :: ReIterationCoeff, ImIterationCoeff, RePotentialCoeff, ImPotentialCoeff 
  
  SAVE CurrentDoneTime, ReIterationCoeff, ImIterationCoeff

  LOGICAL ::  PotentialFlowBC, DDPreconditioning, Found
  TYPE(GaussIntegrationPoints_t) :: IP


  !------------------------------------------------------------------------------
  !    Check if version number output is requested
  !------------------------------------------------------------------------------
  !IF ( .NOT. AllocationsDone ) THEN
  !  IF ( ListGetLogical( GetSimulation(), 'Output Version Numbers', GotIt ) ) THEN
  !    CALL Info( 'AcousticsSolver', 'Acoustics version:', Level = 0 ) 
  !    CALL Info( 'AcousticsSolver', VersionID, Level = 0 ) 
  !    CALL Info( 'AcousticsSolver', ' ', Level = 0 ) 
  !  END IF
  !END IF

  !------------------------------------------------------------------------------
  ! Get variables needed for solution
  !------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  !Solver % Matrix % COMPLEX = .TRUE.

  ! Nullify pointer
  HSol => Null()
  Flow     => Solver % Variable % Values
  FlowPerm => Solver % Variable % Perm
  LocalNodes = COUNT( FlowPerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN
  
  ! PRINT *, 'The systems size is ', LocalNodes

  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS
  Norm = Solver % Variable % Norm
  Dofs = Solver % Variable % DOFs
  VelocityComponents = CoordinateSystemDimension()
  VelocityDofs = VelocityComponents*2
  IF (Dofs /= VelocityDofs + 4) THEN
    CALL Warn('AcousticsSolver', 'Inconsistent number of Variable Dofs')
  END IF


  !------------------------------------------------------------------
  ! Check whether scanning option for coupled simulations is used...
  !------------------------------------------------------------------
  ScanningUsed = .FALSE.     
  SimulationType = ListGetString(Model % Simulation, 'Simulation Type', GotIt)
  IF (SimulationType == 'scanning') THEN
     ScanningUsed = .TRUE.
     NewTimeStep = .FALSE.
     FirstTimeStep = .FALSE.
     TimeVar => VariableGet( Model % Mesh  % Variables, 'Time' ) 
     IF (.NOT. ASSOCIATED(TimeVar) ) & 
          CALL Fatal( 'AcousticsSolver',  'Error in reading time variable.' )

     IF ( NINT(TimeVar % Values(1)) == 1) THEN
        FirstTimeStep = .TRUE.          
     ELSE
        IF ( (NINT(TimeVar % Values(1)) - CurrentDoneTime) == 2) THEN 
           NewTimeStep = .TRUE.
           CurrentDoneTime = CurrentDoneTime + 1
        END IF
     END IF
  END IF

  !-----------------------------------------------------------------------------
  ! If an external BEM solver is used, examine the contents of the COMPA_STATUS
  ! file and modify the ELMER_STATUS file
  !-----------------------------------------------------------------------------
  BEMCoupling = ListGetLogical( Model % Simulation, 'BEM Coupling', GotIt )
  IF (.NOT. GotIt) BEMCoupling = .FALSE.
  SolverCalls = SolverCalls + 1
  IF ( BEMCoupling ) THEN

    ElmerStatusFile = ListGetString(Model % Simulation, 'Elmer Status File', GotIt)
    IF ( .NOT. GotIt ) ElmerStatusFile = 'ELMER_STATUS'
    CompaStatusFile = ListGetString(Model % Simulation, 'Compa Status File', GotIt)
    IF ( .NOT. GotIt ) CompaStatusFile = 'COMPA_STATUS'

    CompaMeshPrefix = ListGetString(Model % Simulation, 'Compa Mesh File Prefix', GotIt)
    IF ( .NOT. GotIt ) CompaMeshPrefix = 'mesh'

    IF (.NOT. FirstVisit) THEN
      OPEN( 10, FILE = CompaStatusFile, status='OLD', IOSTAT = istat )
      IF ( istat /= 0) &
          CALL Fatal( 'AcousticsSolver', 'Cannot open Compa Status File' )           
      READ(10,'(A)') str2
      CLOSE(10)
      i = INDEX(str2,'$') + 1    
      IF ( str2(i:i+3) == 'STOP' .OR. str2(i:i+4) == 'PANIC' ) &
          CALL Fatal( 'AcousticsSolver', 'CompaSolver is not willing to continue' )           
    END IF

    OPEN( 10, FILE = ElmerStatusFile, status='REPLACE', IOSTAT = istat )
    IF ( istat /= 0) THEN
      CALL Fatal( 'AcousticsSolver', 'Cannot open Elmer Status File' )    
    ELSE
      WRITE( 10, '(A)', ADVANCE='NO') '$WAIT'
    END IF
    CLOSE(10)

  ELSE

     DDPreconditioning = ListGetLogical( Model % Simulation, 'DD Preconditioning', GotIt )
     IF (DDPreconditioning) THEN
        HSol => VariableGet( Model % Variables, 'Pres') 
        IF( ASSOCIATED(HSol) ) THEN 
           H => HSol % Values
           HPerm => HSol % Perm
        ELSE
           CALL Fatal('AcousticsSolver', 'Helmholtz solution Pres is not available')
        END IF

     ELSE
        !------------------------------------------------------------------------------
        ! Get (Elmer) Helmholtz solution for acoustic interface boundary conditions...
        !------------------------------------------------------------------------------
        Tmp => Solver % Mesh % Variables
        DO WHILE( ASSOCIATED(tmp) )
           IF ( tmp % NameLen == 4 ) THEN
              IF ( tmp % Name == 'pres' ) THEN
                 IF ( Tmp % Valid ) THEN
                    HSol => Tmp
                 END IF
                 EXIT
              END IF
           END IF
           tmp => tmp % Next
        END DO

        IF (ASSOCIATED(HSol)) THEN
           H => HSol % Values
           HPerm => HSol % Perm
        ELSE
           CALL Info('AcousticsSolver', 'Helmholtz solution is not available')
           CALL Info('AcousticsSolver', 'Acoustic interface boundary conditions cannot be applied')
        END IF
     END IF
  END IF

  RelOrder = ListGetInteger( Solver % Values,'Relative Integration Order',Found )
  IF( Found ) THEN
    IF( ParEnv % MyPe == 0 ) THEN
      CurrentElement => Solver % Mesh % Elements( Solver % ActiveElements(1) )
      IP = GaussPoints( CurrentElement )
      i = IP % n 
      IP = GaussPoints( CurrentElement, RelOrder = RelOrder )
      j = IP % n 
      CALL Info('AcousticsSolver','Number of Gauss Points: '&
          //TRIM(I2S(j))//' (vs. '//TRIM(I2S(i))//')',Level=6)
    END IF
  END IF

  
  !-------------------------------------------------------------------------------
  ! If the Helmholtz solution is done with an external BEM solver, create an array
  ! for the node indices on the acoustic interface.
  !-------------------------------------------------------------------------------
  IF ( BEMCoupling .AND.  (.NOT. FirstVisit) ) THEN

    WRITE(DataFile,'(A,A)') TRIM(CompaMeshPrefix), '.P'
    OPEN( 10, FILE = DataFile, status='OLD', IOSTAT = istat )
    IF ( istat /= 0) THEN
      WRITE( Message, * ) 'Cannot open file ', DataFile     
      CALL Fatal( 'AcousticsSolver', Message )
    END IF

    j = 0
    DO
      READ( 10, *, IOSTAT = istat) i
      IF (istat == 0) THEN
        j = j + 1
      ELSE
        EXIT
      END IF
    END DO
    !------------------------------------------------------------------------------
    IF (.NOT. BEMNodesCreated) THEN
      ALLOCATE( AcousticInterfaceNodes(j), AcousticInterfaceResults(j,2), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
      AcousticInterfaceNodes(1:j) = 0
    END IF
    !------------------------------------------------------------------------------
    !PRINT *, j,' nodes were found' 
    REWIND 10
    
    OPEN( 20, FILE = 'mesh.nodes', status='OLD', IOSTAT = istat )
    IF ( istat /= 0) &
        CALL Fatal( 'AcousticsSolver', 'Cannot open file mesh.nodes' )    
    DO i = 1, j
      READ( 10, *) t, ReP, Imp
      !-----------------------------------------
      ! Find the row where this node is defined
      !-----------------------------------------
      k = 0
      DO
        READ( 20, *, IOSTAT = istat) m 
        IF (istat == 0) THEN
          k = k + 1
          IF ( m==t ) THEN
            AcousticInterfaceNodes(i) = k
            AcousticInterfaceResults(i,1) = ReP
            AcousticInterfaceResults(i,2) = -1.0d0*ImP            
            EXIT
          END IF
        ELSE
          WRITE( Message, * ) 'Inconsistent node numbering in the file ', DataFile
          CALL Fatal( 'AcousticsSolver', Message )          
        END IF
      END DO
      REWIND 20
    END DO
    !----------------------------------
    CLOSE(20)
    CLOSE(10)     
    BEMNodesCreated = .TRUE.
  END IF



  !------------------------------------------------------------------------------
  ! Find out whether the block preconditioning is used...
  !------------------------------------------------------------------------------
  dim = CoordinateSystemDimension() 
  BlockPreconditioning = ListGetLogical( Solver % Values, 'Block Preconditioning', GotIt )
  
  IF (BlockPreconditioning) THEN
    CALL Info('AcousticsSolver', 'Block preconditioning will be used.')

    OuterIterationMethod = ListGetString(Solver % Values, 'Outer Iteration Method', GotIt)
    IF ( .NOT. GotIt ) OuterIterationMethod = 'nested gcr'     
    ToleranceRatio = ListGetConstReal( Solver % Values, &
        'Ratio of Convergence Tolerances', GotIt )
    IF (GotIt) THEN
      Tolerance = ToleranceRatio * ListGetConstReal( Solver % Values, &
          'Linear System Convergence Tolerance' )
    ELSE
      Tolerance = ListGetConstReal( Solver % Values, &
          'Linear System Convergence Tolerance' )
    END IF
    MaxIterations = ListGetInteger( Solver % Values, &
        'Max Outer Iterations', GotIt )
    IF ( .NOT. GotIt ) MaxIterations = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations')

    VelocityAssembly = ListGetLogical(Solver % Values, 'Velocity Assembly', GotIt)
    IF (.NOT. GotIt) VelocityAssembly = .FALSE.
    
  END IF

  !------------------------------------------------------------------------------
  ! Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
    !N = Solver % Mesh % MaxElementNodes
    N = Solver % Mesh % MaxElementDOFs
    
    IF ( AllocationsDone ) THEN
      DEALLOCATE(                 &
          ElementNodes % x,      &
          ElementNodes % y,      &
          ElementNodes % z,      &
          ParentNodes % x,      &
          ParentNodes % y,      &
          ParentNodes % z,      &         
          LocalForce,            &
          temp,                  &
          LocalStiffMatrix,      &
          Load,                  &
          HeatSource,            &
          SpecificHeat,          &
          HeatRatio,             &
          Density,               &
          Pressure,              &
          Temperature,           &
          Conductivity,          &
          Viscosity,             &
          Lambda,                &
          BulkViscosity,         &
          Impedance,             &
          WallTemperature,       &
          WallVelocity,          &
          TotalForce,            &
          TotalMoment,           &
          TotalArea,             &
          MomentOutput,          &
          GapIndexes,            &
          NodesOnBoundary)           
      
    END IF

    !-----------------------------------------------------------------
    ! Create new matrices for the block preconditioning...
    !------------------------------------------------------------------
    IF (BlockPreconditioning) THEN
      !---------------------------------------------------------------------------
      ! Make sure that the new matrices are optimized as the primary one.
      !---------------------------------------------------------------------------
      OptimizeBW = ListGetLogical(Solver % Values,'Optimize Bandwidth',gotIt)
      IF(.NOT. GotIt) OptimizeBW = .TRUE.

      IF (VelocityAssembly) &
          AMatrix => CreateMatrix( Model, Solver, Solver % Mesh, FlowPerm, 2*dim, &
          MATRIX_CRS, OptimizeBW, ListGetString( Solver % Values, 'Equation' ) )

      SMatrix => CreateMatrix( Model, Solver, Solver % Mesh, Solver % Variable % Perm, &
          4, MATRIX_CRS, OptimizeBW, ListGetString( Solver % Values, 'Equation' ) )
      
      !------------------------------------------------------
      ! Matrices to construct consistent boundary conditions
      !------------------------------------------------------
      MMatrix => CreateMatrix( Model, Solver, Solver % Mesh, FlowPerm, &
          2, MATRIX_CRS, OptimizeBW, ListGetString( Solver % Values, 'Equation' ) )

      LMatrix => CreateMatrix( Model, Solver, Solver % Mesh, FlowPerm, &
          2, MATRIX_CRS, OptimizeBW, ListGetString( Solver % Values, 'Equation' ) )

      ALLOCATE( SLocal(4*n,4*n), ALocal(2*dim*n,2*dim*n), AnLocal(2*n,2*n), &
          STAT=istat )
      IF ( istat /= 0 ) THEN
        CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
      END IF
    END IF
    
    ALLOCATE( ElementNodes % x( N ), &
        ElementNodes % y( N ),       &
        ElementNodes % z( N ),       &
        ParentNodes % x( N ),        &
        ParentNodes % y( N ),       &
        ParentNodes % z( N ),       &
        LocalForce( Dofs*N ),        &
        temp( N ),                   &
        LocalStiffMatrix( Dofs*N,Dofs*N ), &
        Load( 6,N ),   &
        HeatSource(2,N),               &
        SpecificHeat(N),             &
        HeatRatio(N),                &
        Density(N),                  &
        Pressure(N),        &
        Temperature(N),              &
        Conductivity(N),             &
        Viscosity(N),                &
        Lambda(N),                   &
        BulkViscosity(N),            &
        Impedance( 4,N ),            &
        WallTemperature(N),          &
        WallVelocity(6,N),           &
        TotalForce(Model % NumberOfBCs,6), &
        TotalMoment(Model % NumberOfBCs,6), &
        TotalArea(Model % NumberOfBCs), &
        MomentOutput(Model % NumberOfBCs), &
        ! NodesOnBoundary(N*Solver % Mesh % NumberOfBoundaryElements), &
        GapIndexes(N), &
        STAT=istat )

    ElementNodes % x = 0._dp
    ElementNodes % y = 0._dp
    ElementNodes % z = 0._dp

    IF ( istat /= 0 ) THEN
      CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
    END IF

    AllocationsDone = .TRUE.
  END IF

  !---------------------------------------------------------------------------------
  ! Create the array which contains the indeces of the nodes located on the boundary. 
  ! This is for testing an explicit stabilisation technique 
  !---------------------------------------------------------------------------------
  IF (.FALSE.) THEN
    NodesOnBoundary = 0
    i = 1
    DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
        Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
      
      CurrentElement => Solver % Mesh % Elements(t)
      Model % CurrentElement => CurrentElement
      !------------------------------------------------------------------------------
      ! Check that the dimension of element is suitable for BCs
      !------------------------------------------------------------------------------
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      IF (ANY(FlowPerm(NodeIndexes(1:n)) == 0)) CYCLE
      NodesOnBoundary(i:i+n-1) = NodeIndexes(1:n)
      i = i + n
    END DO
    MaxNodesOnBoundary = i - 1
  END IF

  !---------------------------------------------------------------------------
  ! Initialization related to the block preconditioning
  !---------------------------------------------------------------------------
  IF (BlockPreconditioning) THEN
    IF (VelocityAssembly) CALL CRS_ZeroMatrix( AMatrix ) 
    CALL CRS_ZeroMatrix( SMatrix )
    CALL CRS_ZeroMatrix( MMatrix )
    CALL CRS_ZeroMatrix( LMatrix )
  END IF

  !------------------------------------------------------------------------------
  ! Do some additional initialization and start assembly
  !------------------------------------------------------------------------------
  EquationName = ListGetString( Solver % Values, 'Equation' )
  Bubbles = ListGetLogical( Solver % Values, 'Bubbles', GotIt )
  MiniBubbles = ListGetLogical( Solver % Values, 'Mini Bubbles', GotIt )
  CoordSys = CurrentCoordinateSystem()
  IF (.NOT. (CoordSys==Cartesian .OR. CoordSys==AxisSymmetric) ) THEN
    CALL Fatal('AcousticsSolver',&
        'Currently only Cartesian or cylindrical coordinates are allowed')
  END IF

  !------------------------------------------------------------------------------
  ! Figure out angular frequency 
  !------------------------------------------------------------------------------
  AngularFrequency = GetCReal( Model % Simulation, 'Angular Frequency', GotIt )
  IF(.NOT. GotIt ) THEN
    AngularFrequency = 2*PI*GetCReal( Model % Simulation, 'Frequency',GotIT )
  END IF

  !------------------------------------------------------------------------------
  ! The nodal values of the approximations of 
  ! div(v)-type term and the scaled temperature may be overwritten in such a way 
  ! that the previous solution is used as an initial guess...
  !----------------------------------------------------------------------------- 
  UtilizePreviousSolution = ListGetLogical( Solver % Values, &
      'Utilize Previous Solution', GotIt )

  IF (UtilizePreviousSolution) THEN
    VisitedNodes = .FALSE.
    DO t=1,Solver % Mesh % NumberOfBulkElements
      CurrentElement => Solver % Mesh % Elements(t)
      IF ( .NOT. CheckElementEquation( Model, &
          CurrentElement, EquationName ) ) CYCLE
      Model % CurrentElement => CurrentElement
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      !------------------------------------------------------------------------------
      !   Get equation & material parameters
      !------------------------------------------------------------------------------
      k = ListGetInteger( Model % Bodies( CurrentElement % &
          Bodyid ) % Values, 'Material' )
      Material => Model % Materials(k) % Values

      SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
          n, NodeIndexes )
      HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
          n, NodeIndexes )
      Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
          n, NodeIndexes )    
      Temperature(1:n) = ListGetReal( Material, 'Equilibrium Temperature', &
          n, NodeIndexes )        
      Viscosity(1:n) = ListGetReal( Material, 'Viscosity', &
          n, NodeIndexes )   
      Lambda(1:n) = -2.0d0/3.0d0 * Viscosity(1:n)
      BulkViscosity(1:n) = ListGetReal( Material, 'Bulk Viscosity', &
          n, NodeIndexes, GotIt )
      IF (GotIt) Lambda(1:n) = BulkViscosity(1:n) - 2.0d0/3.0d0 * Viscosity(1:n)

      DO i=1,n
        j = FlowPerm(NodeIndexes(i))
        IF ( .NOT. VisitedNodes(j) ) THEN
          VisitedNodes(j) = .TRUE.
          !-------------------------------------
          ! Rescaling of the temperature...
          !--------------------------------------
          A = (HeatRatio(i)-1.0d0) * SpecificHeat(i) / AngularFrequency * &
              CMPLX( Flow( (j-1)*Dofs+VelocityDofs+1 ), &
              Flow( (j-1)*Dofs+VelocityDofs+2 ), kind=dp )
          Flow( (j-1)*Dofs+VelocityDofs+1 ) = REAL(A)
          Flow( (j-1)*Dofs+VelocityDofs+2 ) = AIMAG(A)
          !-------------------------------------
          ! Rescaling of the pressure...
          !--------------------------------------
          A = CMPLX( 1.0d0, Lambda(i)*AngularFrequency/  &
              ( (HeatRatio(i)-1.0d0) * SpecificHeat(i) * Density(i) * Temperature(i) ), kind=dp ) * &
              (1.0d0/(AngularFrequency * Density(i)) * &
              CMPLX( Flow( (j-1)*Dofs+VelocityDofs+3 ), Flow( (j-1)*Dofs+VelocityDofs+4 ), kind=dp ) - &
              CMPLX( Flow( (j-1)*Dofs+VelocityDofs+1 ), Flow( (j-1)*Dofs+VelocityDofs+2 ), kind=dp ) )
          Flow( (j-1)*Dofs+VelocityDofs+3 ) = REAL(A)
          Flow( (j-1)*Dofs+VelocityDofs+4 ) = AIMAG(A)
        END IF
      END DO
    END DO
  END IF

  totat = 0.0d0
  totst = 0.0d0
  at  = CPUTime()
  at0 = RealTime()
  TimeEstimated = .FALSE.

  CALL Info( 'AcousticsSolver', ' ', Level=4 )
  CALL Info( 'AcousticsSolver', '-------------------------------------', &
      Level=4 )
  WRITE( Message, * ) 'Frequency (Hz): ', AngularFrequency/(2*PI)
  CALL Info( 'AcousticsSolver', Message, Level=4 )
  CALL Info( 'AcousticsSolver', '-------------------------------------', &
      Level=4 )
  CALL Info( 'AcousticsSolver', ' ', Level=4 )
  CALL Info( 'AcousticsSolver', 'Starting Assembly', Level=4 )

  CALL DefaultStart()
  CALL DefaultInitialize()
  !InitializeToZero( StiffMatrix, ForceVector )

  !------------------------------------------------------------------------------
  DO t=1,Solver % NumberOfActiveElements
  !------------------------------------------------------------------------------

    IF ( RealTime() - at0 > 1.0 ) THEN
      IF( .NOT. TimeEstimated ) THEN
        IF( CPUTime() - at > 1.0_dp .AND. t < Solver % NumberOfActiveElements / 2) THEN
          WRITE(Message,'(a,F8.2)' ) 'Estimated assembly time (s): ', &
              (CPUTime()-at)*(Solver % NumberOfActiveElements) / t
          CALL Info( 'AcousticsSolver', Message, Level=5 )        
          TimeEstimated = .TRUE.
        END IF
      END IF

      IF( MOD( t, (Solver % NumberOfActiveElements / 12)) == 0) THEN
        WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
            (Solver % NumberOfActiveElements)), ' % done'
        CALL Info( 'AcousticsSolver', Message, Level=5 )
        at0 = RealTime()
      END IF
    END IF


    CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))

    Model % CurrentElement => CurrentElement
    n = CurrentElement % TYPE % NumberOfNodes
    nb = GetElementNOFBDOFs()
    IF ( (n /= GetElementNOFDOFs()) .AND. (nb >0) ) THEN
       WRITE( Message, * ) 'The keyword Bubbles in Global System should have the value False'
       CALL Fatal( 'AcousticsSolver', Message)
    END  IF

    NodeIndexes => CurrentElement % NodeIndexes
 
    ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

    !------------------------------------------------------------------------------
    ! Get equation & material parameters
    !------------------------------------------------------------------------------
    k = ListGetInteger( Model % Bodies( CurrentElement % &
        Bodyid ) % Values, 'Material' )
    Material => Model % Materials(k) % Values

    SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
        n, NodeIndexes )
    HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
        n, NodeIndexes )
    Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
        n, NodeIndexes )
    Temperature(1:n) = ListGetReal( Material, 'Equilibrium Temperature', &
        n, NodeIndexes )        
    Conductivity(1:n) = ListGetReal( Material, 'Heat Conductivity', &
        n, NodeIndexes )   
    Viscosity(1:n) = ListGetReal( Material, 'Viscosity', &
        n, NodeIndexes )   
    Lambda(1:n) = -2.0d0/3.0d0 * Viscosity(1:n)
    BulkViscosity(1:n) = ListGetReal( Material, ' Bulk Viscosity', &
        n, NodeIndexes, GotIt )
    IF (GotIt) Lambda(1:n) = BulkViscosity(1:n) - 2.0d0/3.0d0 * Viscosity(1:n)
    Pressure(1:n) = (HeatRatio(1:n)-1.0d0)* SpecificHeat(1:n) * Density(1:n) * Temperature(1:n)

    !------------------------------------------------------------------------------
    !   The heat source and body force at nodes
    !------------------------------------------------------------------------------
    Load = 0.0d0
    HeatSource = 0.0d0
    k = ListGetInteger( Model % Bodies( CurrentElement % BodyId ) % &
                       Values, 'Body Force', GotIt )
    IF ( k > 0 ) THEN
      HeatSource(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Re Heat Source', n, NodeIndexes, GotIt )
      HeatSource(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Im Heat Source', n, NodeIndexes, GotIt )
      Load(1,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Re Body Force 1', n, NodeIndexes, GotIt )
      Load(2,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Im Body Force 1', n, NodeIndexes, GotIt )
      Load(3,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Re Body Force 2', n, NodeIndexes, GotIt )
      Load(4,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Im Body Force 2', n, NodeIndexes, GotIt )
      Load(5,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Re Body Force 3', n, NodeIndexes, GotIt )
      Load(6,1:n) = ListGetReal( Model % BodyForces(k) % Values, &
          'Im Body Force 3', n, NodeIndexes, GotIt )
    END IF

    !------------------------------------------------------------------------------
    !   Get element local matrix and rhs vector
    !------------------------------------------------------------------------------
    CALL LocalMatrix(  LocalStiffMatrix, LocalForce, AngularFrequency, &
        SpecificHeat, HeatRatio, Density,                     &
        Temperature, Conductivity, Viscosity, Lambda,                   &
        HeatSource, Load, Bubbles, MiniBubbles, CurrentElement, n, ElementNodes,     &
        Dofs, nb)  

    !------------------------------------------------------------------------------
    !   Update global matrix and rhs vector from local matrix & vector
    !------------------------------------------------------------------------------
    CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
        ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )


    !-----------------------------------------------------------------------------
    ! This is for testing an explicit stabilisation technique
    !-----------------------------------------------------------------
    IF (.FALSE.) THEN
      !------------------------------------------------------------------
      ! Check whether the element has nodes located on boundary...
      !------------------------------------------------------------------
      NodeOnBoundary(1:10) = .FALSE.
      DO i = 1,n
        DO j = 1, MaxNodesOnBoundary
          IF ( NodeIndexes(i) == NodesOnBoundary(j) ) THEN
            NodeOnBoundary(i) = .TRUE.
            EXIT
          END IF
        END DO
      END DO

      !----------------------------------
      ! Compute the stabilisation matrix 
      !----------------------------------
      CALL ExplicitStabilisationMatrix(  LocalStiffMatrix, LocalForce, AngularFrequency, &
          SpecificHeat, HeatRatio, Density,                    &
          Temperature, Conductivity, Viscosity, Lambda,                  &
          HeatSource, Load, Bubbles, CurrentElement, n, ElementNodes,    &
          Dofs, NodeOnBoundary)

      !-------------------
      ! Do the assembly 
      !-------------------
      CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
          ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )
    END IF
    !----------------- the end of testing ---------------------------------------



    !-----------------------------------------------------------------------------
    ! Do the assembly for the preconditioners...
    !----------------------------------------------------------------------------
    IF (BlockPreconditioning) THEN
      IF (VelocityAssembly) THEN
        CALL CoupledVelocityMatrix( ALocal, Viscosity, AngularFrequency, Density, &
            CurrentElement, n, dim)
        CALL UpdateGlobalPreconditioner( AMatrix, ALocal, n, 2*dim, &
            FlowPerm( CurrentElement % NodeIndexes ) )
      END IF

      CALL PressureLaplaceMatrix( AnLocal, Viscosity, AngularFrequency, Density, &
          CurrentElement, n, dim)      
      CALL UpdateGlobalPreconditioner( LMatrix, AnLocal, n, 2, &
          FlowPerm( CurrentElement % NodeIndexes ) )

      CALL PressureMassMatrix( AnLocal, CurrentElement, n, dim)      
      CALL UpdateGlobalPreconditioner( MMatrix, AnLocal, n, 2, &
          FlowPerm( CurrentElement % NodeIndexes ) )

      CALL SchurComplementMatrix( SLocal, AngularFrequency, SpecificHeat, &
          HeatRatio, Density, Pressure, Temperature, Conductivity, Viscosity, Lambda, &
          CurrentElement, n, dim)
      CALL UpdateGlobalPreconditioner( SMatrix, SLocal, n, 4, &
          FlowPerm( CurrentElement % NodeIndexes ) )

    END IF
  !------------------------------------------------------------------------------
  END DO
  !------------------------------------------------------------------------------

  CALL DefaultFinishBulkAssembly()


  !------------------------------------------------------------------------------
  ! Compute the matrix corresponding to the integral over boundaries:
  !------------------------------------------------------------------------------
  DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
      Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements

    CurrentElement => Solver % Mesh % Elements(t)
    Model % CurrentElement => CurrentElement

    IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel % Solver) ) CYCLE    

    !------------------------------------------------------------------------------
    ! Extract the parent element to find its material parameters... 
    !----------------------------------------------------------------------------- 
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    IF (ANY(FlowPerm(NodeIndexes(1:n)) == 0)) CYCLE
    
    DO i=1,Model % NumberOfBCs
      IF ( CurrentElement % BoundaryInfo % Constraint == &
          Model % BCs(i) % Tag ) THEN

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
        
        Parent => CurrentELement % BoundaryInfo % Left
        stat = ASSOCIATED( Parent )
        IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

        IF ( .NOT. stat) THEN
          Parent => CurrentELement % BoundaryInfo % Right
          stat = ASSOCIATED( Parent )
          IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

          IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
              'No parent element can be found for given boundary element' )
        END IF
        ! IF ( .NOT. CheckElementEquation( Model, Parent, EquationName ) ) CYCLE
    
        k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
            'Material' )
        Material => Model % Materials(k) % Values
        
        SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
            n, NodeIndexes )
        HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
            n, NodeIndexes )
        Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
            n, NodeIndexes )
        Temperature(1:n) = ListGetReal( Material, &
            'Equilibrium Temperature', n, NodeIndexes )        
        Conductivity(1:n) = ListGetReal( Material, 'Heat Conductivity', &
            n, NodeIndexes )   
        Pressure(1:n) = (HeatRatio(1:n)-1.0d0)* SpecificHeat(1:n) * Density(1:n) * Temperature(1:n)


        AcousticInterface = ListGetLogical( Model % BCs(i) % Values, &
            'Acoustic Interface', GotIt )
        ImpedanceBCIteration = ListGetLogical( Model % BCs(i) % Values, &
            'Impedance Update', GotIt )
        IF (.NOT. GotIt) ImpedanceBCIteration = .TRUE.

        PotentialFlowBC = ListGetLogical( Model % BCs(i) % Values, &
            'Potential Flow BC', GotIt )


        IF ( (AcousticInterface .AND. (.NOT. FirstVisit)) ) THEN
          IF (ImpedanceBCIteration) THEN
            !----------------------------------------------------------------
            ! Compute new approximations to the impedances. First
            ! compute the average of the normal velocity over the element.
            !----------------------------------------------------------------
            DO j=1,n
              k = FlowPerm( NodeIndexes(j) ) 
              WallVelocity(1,j) =  Flow( (k-1)*(dim*2+4)+1 )
              WallVelocity(2,j) =  Flow( (k-1)*(dim*2+4)+2 )
              WallVelocity(3,j) =  Flow( (k-1)*(dim*2+4)+3 )
              WallVelocity(4,j) =  Flow( (k-1)*(dim*2+4)+4 )
              IF (dim > 2) THEN
                WallVelocity(5,j) =  Flow( (k-1)*(dim*2+4)+5 )
                WallVelocity(6,j) =  Flow( (k-1)*(dim*2+4)+6 )
              ELSE
                WallVelocity(5,j) =  0.0d0
                WallVelocity(6,j) =  0.0d0
              END IF
            END DO

            CALL ComputeAverageVelocity( CurrentElement, n, WallVelocity, AverVel)

            IF (BEMCoupling) THEN
              DO j=1,n
                istat = 1
                DO k = 1, SIZE(AcousticInterfaceNodes)
                  IF ( AcousticInterfaceNodes(k) == NodeIndexes(j) )THEN
                    Load(1,j) = AcousticInterfaceResults(k,1)
                    Load(2,j) = AcousticInterfaceResults(k,2)
                    !ZAppr = -CMPLX( AcousticInterfaceResults(k,1), & 
                    !    AcousticInterfaceResults(k,2), kind=dp) / AverVel
                    !Impedance(1,j) = REAL(ZAppr)
                    !Impedance(2,j) = AIMAG(ZAppr)
                    !ZAppr = CMPLX( 0.0d0, AngularFrequency*Density(j), kind=dp) / ZAppr
                    !Impedance(3,j) = REAL(ZAppr)
                    !Impedance(4,j) = AIMAG(ZAppr)                      
                    istat = 0
                    EXIT
                  END IF
                END DO
                IF (istat /= 0) & 
                    CALL Fatal('AcousticsSolver', 'Helmholtz solution is not available on the node')             
              END DO
              ! Impedance calculation by pressure averaging
              ZAppr = -CMPLX( SUM( Load(1,1:n) )/n, & 
                  SUM( Load(2,1:n) )/n, kind=dp) / AverVel
              DO j=1,n
                Impedance(1,j) = REAL(ZAppr)
                Impedance(2,j) = AIMAG(ZAppr)
                ZAppr2 = CMPLX( 0.0d0, AngularFrequency*Density(j), kind=dp) / ZAppr
                Impedance(3,j) = REAL(ZAppr2)
                Impedance(4,j) = AIMAG(ZAppr2)
              END DO

            ELSE
              IF ( ALL( HPerm( NodeIndexes(1:n) ) > 0 ) ) THEN
                DO j=1,n
                  k = HPerm( NodeIndexes(j) ) 
                  ZAppr = -CMPLX( H(2*k-1), H(2*k), kind=dp) / AverVel
                  Impedance(1,j) = REAL(ZAppr)
                  Impedance(2,j) = AIMAG(ZAppr)
                  ZAppr = CMPLX( 0.0d0, AngularFrequency*Density(j), kind=dp) / ZAppr
                  Impedance(3,j) = REAL(ZAppr)
                  Impedance(4,j) = AIMAG(ZAppr)                
                END DO
              END IF
            END IF

          ELSE
            Impedance(1,1:n) = 0.0d0
            Impedance(2,1:n) = 0.0d0
            Impedance(3,1:n) = 0.0d0
            Impedance(4,1:n) = 0.0d0
          END IF

        ELSE
          Impedance(1,1:n) = ListGetReal( Model % BCs(i) % Values, &
              'Re Specific Acoustic Impedance', n, NodeIndexes, GotIt )
          Impedance(2,1:n) = ListGetReal( Model % BCs(i) % Values, &
              'Im Specific Acoustic Impedance', n, NodeIndexes, GotIt )
          Impedance(3,1:n) = ListGetReal( Model % BCs(i) % Values, &
              'Re Specific Thermal Impedance', n, NodeIndexes, GotIt )
          Impedance(4,1:n) = ListGetReal( Model % BCs(i) % Values, &
              'Im Specific Thermal Impedance', n, NodeIndexes, GotIt )
        END IF
          
        Load(1,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Re Surface Traction 1', n, NodeIndexes, GotIt )
        Load(2,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Im Surface Traction 1', n, NodeIndexes, GotIt )
        Load(3,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Re Surface Traction 2', n, NodeIndexes, GotIt )
        Load(4,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Im Surface Traction 2', n, NodeIndexes, GotIt )
        Load(5,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Re Surface Traction 3', n, NodeIndexes, GotIt )
        Load(6,1:n) = ListGetReal( Model % BCs(i) % Values, &
            'Im Surface Traction 3', n, NodeIndexes, GotIt )

        !------------------------------------------------------------------------------
        ! Get element local matrix and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrixBoundary(  LocalStiffMatrix, LocalForce, &
            AngularFrequency , SpecificHeat, HeatRatio, Density, &
            Pressure, Temperature, Conductivity,     &
            Impedance, Load, CurrentElement, n, ElementNodes, Dofs)
        
        !------------------------------------------------------------------------------
        ! Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
            ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )
        !------------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        ! Boundary conditions on acoustic interface.
        ! The iteration based on updating the impedance seems to perform better,
        ! so the default is that this branch is never executed. However,
        ! in the case of domain decomposition preconditioning this is performed.
        !---------------------------------------------------------------------------
        IF ( (AcousticInterface .AND. (.NOT. FirstVisit) .AND. (.NOT. ImpedanceBCIteration) ) .OR. &
             PotentialFlowBC ) THEN
           IF (BEMCoupling) THEN
              DO j=1,n
                 istat = 1
                 DO k = 1, SIZE(AcousticInterfaceNodes)
                    IF ( AcousticInterfaceNodes(k) == NodeIndexes(j) )THEN
                       Load(1,j) = AcousticInterfaceResults(k,1)
                       Load(2,j) = AcousticInterfaceResults(k,2)
                       istat = 0
                       EXIT
                    END IF
                 END DO
                 IF (istat /= 0) & 
                      CALL Fatal('AcousticsSolver', 'Helmholtz solution is not available on the node')             
              END DO

              CALL LocalInterfaceMatrix( LocalStiffMatrix, LocalForce, &
                   AngularFrequency , SpecificHeat, HeatRatio, Density, &
                   Pressure, Temperature, Conductivity,     &
                   Impedance, Load, CurrentElement, n, ElementNodes, Dofs)
            
              CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
                   ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )
            
              !--------------------------------------------------------------------
              ! Set Dirichlet BCs for temperature...
              !--------------------------------------------------------------------
              DO j=1,n
                 m = Solver % Variable % Perm( NodeIndexes(j) )
                 CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+1 )
                 CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+1, &
                      (m-1)*(dim*2+4)+dim*2+1, 1.0d0 )
                 CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+2 )
                 CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+2, &
                      (m-1)*(dim*2+4)+dim*2+2, 1.0d0 )  
                 Solver % Matrix % RHS( (m-1)*(dim*2+4)+dim*2+1 ) = &
                      (HeatRatio(j)-1.0d0)/(Density(j)*HeatRatio(j)*AngularFrequency) * Load(1,j)
                 Solver % Matrix % RHS( (m-1)*(dim*2+4)+dim*2+2 ) = &
                      (HeatRatio(j)-1.0d0)/(Density(j)*HeatRatio(j)*AngularFrequency) * Load(2,j)
              END DO

           ELSE
              IF ( ALL( HPerm( NodeIndexes(1:n) ) > 0 ) ) THEN
                 DO j=1,n
                    k = HPerm( NodeIndexes(j) ) 
                    Load(1,j) = H(2*k-1)
                    Load(2,j) = H(2*k)
                 END DO
                 
                 CALL LocalInterfaceMatrix( LocalStiffMatrix, LocalForce, &
                      AngularFrequency , SpecificHeat, HeatRatio, Density, &
                      Pressure, Temperature, Conductivity,     &
                      Impedance, Load, CurrentElement, n, ElementNodes, Dofs)
            
                 CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
                      ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )
            
                 !--------------------------------------------------------------------
                 ! Set Dirichlet BCs for temperature...
                 !--------------------------------------------------------------------
                 DO j=1,n
                    k = HPerm( NodeIndexes(j) )
                    m = Solver % Variable % Perm( NodeIndexes(j) )
                    CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+1 )
                    CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+1, &
                         (m-1)*(dim*2+4)+dim*2+1, 1.0d0 )
                    CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+2 )
                    CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+2, &
                         (m-1)*(dim*2+4)+dim*2+2, 1.0d0 )  
                    Solver % Matrix % RHS( (m-1)*(dim*2+4)+dim*2+1 ) = &
                         (HeatRatio(j)-1.0d0)/(Density(j)*HeatRatio(j)*AngularFrequency) * H(2*k-1)
                    Solver % Matrix % RHS( (m-1)*(dim*2+4)+dim*2+2 ) = &
                         (HeatRatio(j)-1.0d0)/(Density(j)*HeatRatio(j)*AngularFrequency) * H(2*k)
                 END DO

              ELSE
                 CALL Fatal('AcousticsSolver', 'Helmholtz solution is not available on boundary')
              END IF
           END IF
        END IF

        !-----------------------------------------------------------------------------
        ! Do additional assembly for the preconditioners in the case of impedance BC
        !----------------------------------------------------------------------------
        IF (BlockPreconditioning .AND. (ANY(Impedance(1:4,1:n) /= 0.0d0)) ) THEN
          IF (VelocityAssembly) THEN
            CALL VelocityImpedanceMatrix(  ALocal, AngularFrequency, Density, &
                Impedance, CurrentElement, n, ElementNodes, dim)
            CALL UpdateGlobalPreconditioner( AMatrix, ALocal, n, 2*dim, &
                Solver % Variable % Perm( CurrentElement % NodeIndexes ) )
          END IF
          !------------------------------------------------------------------------
          ! Impedance BC's for the Schur complement
          !------------------------------------------------------------------------
          CALL SchurComplementImpedanceMatrix( SLocal, Impedance, AngularFrequency, &
              SpecificHeat, HeatRatio, Density, Conductivity, &
              CurrentElement, n, ElementNodes, dim)
          CALL UpdateGlobalPreconditioner( SMatrix, SLocal, n, 4, &
              Solver % Variable % Perm( CurrentElement % NodeIndexes ) )
        END IF

      END IF
    END DO
  !------------------------------------------------------------------------------
  END DO
  !------------------------------------------------------------------------------

 

  !------------------------------------------------------------------------------
  !    Slip boundary conditions
  !------------------------------------------------------------------------------
  DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
      Solver % Mesh % NumberOfBulkElements +  &
      Solver % Mesh % NumberOfBoundaryElements
    CurrentElement => Solver % Mesh % Elements(t)
    Model % CurrentElement => CurrentElement
    IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel% Solver) ) CYCLE    

    !------------------------------------------------------------------------------
    ! Extract the parent element to find its material parameters... 
    !----------------------------------------------------------------------------- 
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    IF (ANY(FlowPerm(NodeIndexes(1:n)) == 0)) CYCLE

    !------------------------------------------------------------------------------
    DO i=1,Model % NumberOfBCs
      IF ( CurrentElement % BoundaryInfo % Constraint == &
          Model % BCs(i) % Tag ) THEN
        !------------------------------------------------------------------------------
        SlipBoundary = ListGetLogical( Model % BCs(i) % Values, &
            'Slip Boundary', GotIt )
        IF ( .NOT. SlipBoundary) CYCLE

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

        Parent => CurrentELement % BoundaryInfo % Left
        stat = ASSOCIATED( Parent )
        IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

        IF ( .NOT. stat) THEN
          Parent => CurrentELement % BoundaryInfo % Right
          stat = ASSOCIATED( Parent )
          IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

          IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
              'No parent element can be found for given boundary element' )
        END IF

        k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
            'Material' )
        Material => Model % Materials(k) % Values
        SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
            n, NodeIndexes )
        HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
            n, NodeIndexes )
        Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
            n, NodeIndexes )
        Conductivity(1:n) = ListGetReal( Material, 'Heat Conductivity', &
            n, NodeIndexes )   
        Temperature(1:n) = ListGetReal( Material, &
            'Equilibrium Temperature', n, NodeIndexes )
        Pressure(1:n) = (HeatRatio(1:n)-1.0d0)* SpecificHeat(1:n) * Density(1:n) * Temperature(1:n)
    
        WallTemperature(1:n) = ListGetReal(Model % BCs(i) % Values, &
            'Reference Wall Temperature', n, NodeIndexes )
        WallVelocity(1,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Re Reference Wall Velocity 1', n, NodeIndexes, GotIt)
        WallVelocity(2,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Im Reference Wall Velocity 1', n, NodeIndexes, GotIt )
        WallVelocity(3,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Re Reference Wall Velocity 2', n, NodeIndexes, GotIt )
        WallVelocity(4,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Im Reference Wall Velocity 2', n, NodeIndexes, GotIt )
        WallVelocity(5,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Re Reference Wall Velocity 3', n, NodeIndexes, GotIt )
        WallVelocity(6,1:n) = ListGetReal(Model % BCs(i) % Values, & 
            'Im Reference Wall Velocity 3', n, NodeIndexes, GotIt ) 

        
        SlipCoefficient1 = ListGetConstReal( Model % BCs(i) % Values, &
            'Tangential Momentum Accommodation Coefficient')
        SlipCoefficient2 = ListGetConstReal( Model % BCs(i) % Values, &
            'Energy Accommodation Coefficient')
        SlipCoefficient3 = ListGetConstReal( Model % BCs(i) % Values, &
            'Normal Momentum Accommodation Coefficient')

        !------------------------------------------------------------------------------
        !  Get element local matrix and rhs vector
        !------------------------------------------------------------------------------
        IF (CoordSys==Cartesian .OR. CoordSys==AxisSymmetric) THEN
          CALL SlipMatrix(  LocalStiffMatrix, LocalForce, SpecificHeat, &
              HeatRatio, Density, Conductivity, Pressure, Temperature, &
              AngularFrequency, WallTemperature, WallVelocity, SlipCoefficient1, &
              SlipCoefficient2, SlipCoefficient3, CurrentElement, n, ElementNodes, Dofs)
        END IF
        !------------------------------------------------------------------------------
        !          Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
            ForceVector, LocalForce, n, Dofs, FlowPerm(NodeIndexes) )
        !------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------
        ! Do an additional assembly for the preconditioners
        !----------------------------------------------------------------------------
        IF (BlockPreconditioning) THEN
          IF (VelocityAssembly) THEN
            CALL VelocitySlipMatrix(  ALocal, SpecificHeat, &
                HeatRatio, Density, Temperature, &
                AngularFrequency, WallTemperature, SlipCoefficient1, &
                CurrentElement, n, ElementNodes, dim )
            CALL UpdateGlobalPreconditioner( AMatrix, ALocal, n, 2*dim, &
                Solver % Variable % Perm( CurrentElement % NodeIndexes ) )
          END IF
          !------------------------------------------------------------------------
          ! Slip BC for the Schur complement
          !------------------------------------------------------------------------
          CALL SchurComplementSlipMatrix( SLocal, SpecificHeat, &
              HeatRatio, Density, Temperature, AngularFrequency, Conductivity, & 
              WallTemperature, SlipCoefficient2, &
              CurrentElement, n, ElementNodes, dim)
          CALL UpdateGlobalPreconditioner( SMatrix, SLocal, n, 4, &
              Solver % Variable % Perm( CurrentElement % NodeIndexes ) )
        END IF
      END IF
    END DO
  !------------------------------------------------------------------------------
  END DO
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !     CALL FinishAssembly( Solver, ForceVector )
  !------------------------------------------------------------------------------
  !    Dirichlet BCs:
  !------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  ! This is for testing the effect of initial BC's in the case of Helmholtz coupling
  !----------------------------------------------------------------------------------
  !IF ( FirstVisit ) THEN
  !  CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
  !      'Re Initialvelo 1', 1, Dofs, FlowPerm )
  !  CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
  !      'Im Initialvelo 1', 2, Dofs, FlowPerm )   
  !END IF
  !---------- The end of testing ----------------------------------------------------

  DO i = 1, VelocityComponents
    WRITE(VariableName,'(A,A,I1)') 'Re Velocity',' ', i
    CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
        VariableName, (i-1)*2+1, Dofs, FlowPerm )
    WRITE(VariableName,'(A,A,I1)') 'Im Velocity',' ', i
    CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
        VariableName, (i-1)*2+2, Dofs, FlowPerm )
  END DO

  CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, & 
      'Re Temperature', Dofs-3, Dofs, FlowPerm )
  CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, & 
      'Im Temperature', Dofs-2, Dofs, FlowPerm )
  
  IF(.TRUE.) CALL AcousticShellInterface()

  CALL DefaultDirichletBCs()
  
  CALL Info( 'AcousticsSolver', 'Assembly done', Level=4 )

  !------------------------------------------------------------------------------
  ! Set Dirichlet BCs for the normal velocity on the slip boundary:
  ! Nothing is done in the current implementation as it is assumed that 
  ! the user specifies explicitly the boundary condition for the normal 
  ! velocity. 
  !----------------------------------------------------------------------------- 

  !-------------------------------------------------------------------------
  ! Set boundary conditions for the preconditioners...
  !-------------------------------------------------------------------------
  IF (BlockPreconditioning) THEN
    IF (VelocityAssembly) THEN
      CALL SetBoundaryConditions(Model, AMatrix, 'Re Velocity 1', 1, dim*2, &
          FlowPerm)        
      CALL SetBoundaryConditions(Model, AMatrix, 'Im Velocity 1', 2, dim*2,  &
          Solver % Variable % Perm)        
      CALL SetBoundaryConditions(Model, AMatrix, 'Re Velocity 2', 3, dim*2,  &
          FlowPerm)        
      CALL SetBoundaryConditions(Model, AMatrix, 'Im Velocity 2', 4, dim*2,  &
          FlowPerm)      
      IF (dim > 2) THEN
        CALL SetBoundaryConditions(Model, AMatrix, 'Re Velocity 3', 5, dim*2,  &
            FlowPerm)        
        CALL SetBoundaryConditions(Model, AMatrix, 'Im Velocity 3', 6, dim*2,  &
            FlowPerm)  
      END IF
    END IF

    CALL SetBoundaryConditions(Model, SMatrix, 'Re Temperature', 1, 4,  &
        FlowPerm)
    CALL SetBoundaryConditions(Model, SMatrix, 'Im Temperature', 2, 4,  &
        FlowPerm)


    !------------------------------------------------------
    ! ? Testing boundary conditions on traction boundary 
    !------------------------------------------------------
    DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
         Solver % Mesh % NumberOfBulkElements +  &
         Solver % Mesh % NumberOfBoundaryElements
       CurrentElement => Solver % Mesh % Elements(t)
       Model % CurrentElement => CurrentElement
       IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel% Solver) ) CYCLE    

       !------------------------------------------------------------------------------
       ! Extract the parent element to find its material parameters... 
       !----------------------------------------------------------------------------- 
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes
       IF (ANY(FlowPerm(NodeIndexes(1:n)) == 0)) CYCLE

       !------------------------------------------------------------------------------
       DO i=1,Model % NumberOfBCs
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN
             !------------------------------------------------------------------------------
             SlipBoundary = ListGetLogical( Model % BCs(i) % Values, &
                  'Traction Boundary', GotIt )
             IF ( .NOT. SlipBoundary) CYCLE

             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

             Parent => CurrentELement % BoundaryInfo % Left
             stat = ASSOCIATED( Parent )
             IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

             IF ( .NOT. stat) THEN
                Parent => CurrentELement % BoundaryInfo % Right
                stat = ASSOCIATED( Parent )
                IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )
                
                IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
                     'No parent element can be found for given boundary element' )
             END IF

             k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
                  'Material' )
             Material => Model % Materials(k) % Values
             
             
             SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
                  n, NodeIndexes )
             HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
                  n, NodeIndexes )
             Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
                  n, NodeIndexes )
             Temperature(1:n) = ListGetReal( Material, 'Equilibrium Temperature', &
                  n, NodeIndexes )        
             Conductivity(1:n) = ListGetReal( Material, 'Heat Conductivity', &
                  n, NodeIndexes )   
             Viscosity(1:n) = ListGetReal( Material, 'Viscosity', &
                  n, NodeIndexes )   
             Lambda(1:n) = -2.0d0/3.0d0 * Viscosity(1:n)
             BulkViscosity(1:n) = ListGetReal( Material, ' Bulk Viscosity', &
                  n, NodeIndexes, GotIt )
             IF (GotIt) Lambda(1:n) = BulkViscosity(1:n) - 2.0d0/3.0d0 * Viscosity(1:n)
             Pressure(1:n) = (HeatRatio(1:n)-1.0d0)* SpecificHeat(1:n) * Density(1:n) * Temperature(1:n)


             C1 = CMPLX(0.0d0, -0.5d0*Density(1)*AngularFrequency/Viscosity(1), kind=dp)

             C2 = C1 * CMPLX( 1.0d0,AngularFrequency/Pressure(1)*(2.0d0*Viscosity(1)+Lambda(1)), kind=dp ) / &
                  CMPLX( 1.0d0,AngularFrequency/Pressure(1)*Lambda(1), kind=dp )  

             
             reb = REAL(C1)
             imb = AIMAG(C1)
             rec = REAL(C2)
             imc = AIMAG(C2)


             DO j=1,n
                m = Solver % Variable % Perm( NodeIndexes(j) )
                CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+3 )
                CALL ZeroRow( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+4 )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+3, &
                  (m-1)*(dim*2+4)+dim*2+1, reb )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+3, &
                  (m-1)*(dim*2+4)+dim*2+2, -imb )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+3, &
                  (m-1)*(dim*2+4)+dim*2+3, rec )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+3, &
                  (m-1)*(dim*2+4)+dim*2+4, -imc )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+4, &
                  (m-1)*(dim*2+4)+dim*2+1, imb )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+4, &
                  (m-1)*(dim*2+4)+dim*2+2, reb )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+4, &
                  (m-1)*(dim*2+4)+dim*2+3, imc )
                CALL SetMatrixElement( Solver % Matrix, (m-1)*(dim*2+4)+dim*2+4, &
                  (m-1)*(dim*2+4)+dim*2+4, rec )
             END DO

          END IF
       END DO
    !------------------------------------------------------------------------------
    END DO
    !------------------------------------------------------------------------------

  END IF
 

  !------------------------------------------------------------------------------
  !    Solve the linear system...
  !------------------------------------------------------------------------------
  !PrevNorm = Norm
  at = CPUTime() - at
  st = CPUTime()

  IF ( BlockPreconditioning ) THEN
    m = MMatrix % NumberOfRows

    SELECT CASE(OuterIterationMethod)
    CASE ('gcr')
      CALL GCROuterIteration( Solver % Matrix % NumberOfRows, Solver % Matrix, &
          m, AMatrix, SMatrix, MMatrix, LMatrix, Solver % Variable % Values, &
          Solver % Matrix % RHS, MaxIterations, Tolerance, dim, Norm )

    CASE('bicgstab')
      CALL BiCGStabOuterIteration( Solver % Matrix % NumberOfRows, Solver % Matrix, &
          m, AMatrix, SMatrix, MMatrix, LMatrix, Solver % Variable % Values, &
          Solver % Matrix % RHS, MaxIterations, Tolerance, dim, Norm )

    CASE('bicgstabl')
      j = ListGetInteger(Solver % Values, 'BiCGStab Polynomial Degree', GotIt)
      IF ( .NOT. GotIt) j = 4
      CALL BiCGStablOuterIteration( j, Solver % Matrix % NumberOfRows, Solver % Matrix, &
          m, AMatrix, SMatrix, MMatrix, LMatrix, Solver % Variable % Values, &
          Solver % Matrix % RHS, MaxIterations, Tolerance, dim, Norm )

    CASE DEFAULT
      CALL InnerOuterIteration( Solver % Matrix % NumberOfRows, Solver % Matrix, &
          m, AMatrix, SMatrix, MMatrix, LMatrix,&
          Solver % Variable % Values, Solver % Matrix % RHS, MaxIterations, &
          Tolerance, dim, Norm )

    END SELECT
    Solver % Variable % Norm = Norm

  ELSE
    Norm = DefaultSolve()      
  END IF
  st = CPUTime() - st


  IF ( PrevNorm + Norm /= 0.0d0 ) THEN
    RelativeChange = 2*ABS(PrevNorm - Norm) / (PrevNorm + Norm)
  ELSE
    RelativeChange = 0.0d0
  END IF

  CALL Info( 'AcousticsSolver', ' ', Level=4 )
  WRITE( Message, * ) 'Result Norm    : ', Norm
  CALL Info( 'AcousticsSolver', Message, Level=4 )
  WRITE( Message, * ) 'Relative Change: ', RelativeChange
  CALL Info( 'AcousticsSolver', Message, Level=4 )

  WRITE(Message,'(a,F8.2)') ' Assembly: (s)', at  
  CALL Info( 'AcousticsSolver', Message, Level=4 )
  WRITE(Message,'(a,F8.2)') ' Solve:    (s)', st
  CALL Info( 'AcousticsSolver', Message, Level=4 )

  PrevNorm = Norm  

  !------------------------------------------------------------------------------
  ! Overwrite the nodal values of the approximations of the div(v)-type term 
  ! and the scaled temperature in such a way that the resulting nodal values 
  ! are approximations to the true pressure and the temperature.
  !------------------------------------------------------------------------------
  VisitedNodes = .FALSE.
  DO t=1,Solver % Mesh % NumberOfBulkElements
    CurrentElement => Solver % Mesh % Elements(t)
    IF ( .NOT. CheckElementEquation( Model, &
        CurrentElement, EquationName ) ) CYCLE
    Model % CurrentElement => CurrentElement
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    !------------------------------------------------------------------------------
    ! Get equation & material parameters
    !------------------------------------------------------------------------------
    k = ListGetInteger( Model % Bodies( CurrentElement % &
        Bodyid ) % Values, 'Material' )
    Material => Model % Materials(k) % Values

    SpecificHeat(1:n) = ListGetReal( Material, 'Specific Heat', &
        n, NodeIndexes )
    HeatRatio(1:n) = ListGetReal( Material, 'Specific Heat Ratio', &
        n, NodeIndexes )
    Temperature(1:n) = ListGetReal( Material, 'Equilibrium Temperature', &
        n, NodeIndexes )        
    Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
        n, NodeIndexes )
    Viscosity(1:n) = ListGetReal( Material, 'Viscosity', &
        n, NodeIndexes )   
    Lambda(1:n) = -2.0d0/3.0d0 * Viscosity(1:n)
    BulkViscosity(1:n) = ListGetReal( Material, ' Bulk Viscosity', &
        n, NodeIndexes, GotIt )
    IF (GotIt) Lambda(1:n) = BulkViscosity(1:n) - 2.0d0/3.0d0 * Viscosity(1:n)
    Pressure(1:n) = (HeatRatio(1:n)-1.0d0)* SpecificHeat(1:n) * Density(1:n) * Temperature(1:n)
    
    DO i=1,n
      j = FlowPerm(NodeIndexes(i))
      IF ( .NOT. VisitedNodes(j) ) THEN
        VisitedNodes(j) = .TRUE.
        !-------------------------------------
        ! Rescaling of the pressure...
        !--------------------------------------
        A = AngularFrequency * Density(i) * ( CMPLX( Flow( (j-1)*Dofs+VelocityDofs+1 ), &
            Flow( (j-1)*Dofs+VelocityDofs+2 ), kind=dp ) + &
            1.0d0/CMPLX( 1.0d0, Lambda(i)*AngularFrequency/Pressure(i), kind=dp) * &
            CMPLX( Flow( (j-1)*Dofs+VelocityDofs+3 ), &
            Flow( (j-1)*Dofs+VelocityDofs+4 ), kind=dp ) )
        Flow( (j-1)*Dofs+VelocityDofs+3 ) = REAL(A)
        Flow( (j-1)*Dofs+VelocityDofs+4 ) = AIMAG(A)
        !-------------------------------------
        ! Rescaling of the temperature...
        !--------------------------------------
        A = AngularFrequency * Density(i) * Temperature(i)/Pressure(i) * &
            CMPLX( Flow( (j-1)*Dofs+VelocityDofs+1 ), &
            Flow( (j-1)*Dofs+VelocityDofs+2 ), kind=dp )
        Flow( (j-1)*Dofs+VelocityDofs+1 ) = REAL(A)
        Flow( (j-1)*Dofs+VelocityDofs+2 ) = AIMAG(A)
      END IF
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Write an input file for an external BEM solver if it is required 
  !-----------------------------------------------------------------------------
  IF ( BEMCoupling ) THEN
    CALL Info('AcousticsSolver', 'Writing input file for BEM solver')
    
    IF (FirstVisit) THEN
      !----------------------------------------------
      ! Create an array for boundary element indeces 
      !----------------------------------------------
      ALLOCATE( BemElementIndeces(Solver % Mesh % NumberOfBoundaryElements), &
          STAT=istat )
      IF ( istat /= 0 ) CALL Fatal( 'AcousticsSolver', 'Memory allocation error.' )
 
      OPEN( 10, FILE = 'mesh.boundary', status='OLD')
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements    
        READ( 10, *) j
        BemElementIndeces(t) = j
      END DO
      CLOSE(10)
    END IF

    WRITE(DataFile,'(A,A,A)') TRIM(CompaMeshPrefix), '.', 'param'
    OPEN( 10, FILE = DataFile, status='REPLACE')
  
    DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
        Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements

      CurrentElement => Solver % Mesh % Elements(t)
      Model % CurrentElement => CurrentElement

      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      
      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
            Model % BCs(i) % Tag ) THEN
          
          BEMBoundary = ListGetLogical( Model % BCs(i) % Values, &
              'BEM Boundary', GotIt )
          AcousticInterface = ListGetLogical( Model % BCs(i) % Values, &
              'Acoustic Interface', GotIt )

          IF ( BEMBoundary .OR. AcousticInterface) THEN
            IF (AcousticInterface) THEN

              ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
              ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
              ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
        
              Parent => CurrentELement % BoundaryInfo % Left
              stat = ASSOCIATED( Parent )
              IF (stat) stat = ALL(FlowPerm(Parent % NodeIndexes(1:n)) > 0)
              
              IF ( .NOT. stat) THEN
                Parent => CurrentELement % BoundaryInfo % Right
                stat = ASSOCIATED( Parent )
                IF (stat) stat = ALL(FlowPerm(Parent % NodeIndexes(1:n)) > 0)
                IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
                    'No parent element can be found for given boundary element' )
              END IF
              IF ( .NOT. CheckElementEquation( Model, Parent, EquationName ) ) CYCLE
              
              k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
                  'Material' )
              Material => Model % Materials(k) % Values

              Density(1:n) = ListGetReal( Material, 'Equilibrium Density', &
                  n, NodeIndexes )
        
              !------------------------------------------------------------
              ! Compute the average of the normal velocity over the element
              !------------------------------------------------------------
              DO j=1,n
                k = FlowPerm( NodeIndexes(j) ) 
                WallVelocity(1,j) =  Flow( (k-1)*(dim*2+4)+1 )
                WallVelocity(2,j) =  Flow( (k-1)*(dim*2+4)+2 )
                WallVelocity(3,j) =  Flow( (k-1)*(dim*2+4)+3 )
                WallVelocity(4,j) =  Flow( (k-1)*(dim*2+4)+4 )
                IF (dim > 2) THEN
                  WallVelocity(5,j) =  Flow( (k-1)*(dim*2+4)+5 )
                  WallVelocity(6,j) =  Flow( (k-1)*(dim*2+4)+6 )
                ELSE
                  WallVelocity(5,j) =  0.0d0
                  WallVelocity(6,j) =  0.0d0
                END IF
              END DO

              CALL ComputeAverageVelocity( CurrentElement, n, WallVelocity, AverVel)
              AverVel = AverVel * CMPLX(0.0d0, -1.0d0 * AngularFrequency * Density(1), kind=dp )
              AverVel = CONJG(AverVel)

              WRITE( 10, '(I9,6e23.15)',ADVANCE='NO') BemElementIndeces( CurrentElement % ElementIndex - &
                  Solver % Mesh % NumberOfBulkElements ), 0.0d+0, 0.0d+0, 1.0d+0, 0.0d+0, &
                  REAL(AverVel), AIMAG(AverVel)
              WRITE( 10,* ) ''             

            ELSE
              !---------------------------------------------------------  
              ! Handling bem boundaries not shared with the fem domain
              !---------------------------------------------------------
              Load(1,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Re a', n, NodeIndexes, GotIt )
              Load(2,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Im a', n, NodeIndexes, GotIt )
              Load(3,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Re b', n, NodeIndexes, GotIt )
              Load(4,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Im b', n, NodeIndexes, GotIt )
              Load(5,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Re c', n, NodeIndexes, GotIt )
              Load(6,1:n) = ListGetReal( Model % BCs(i) % Values, &
                  'Im c', n, NodeIndexes, GotIt )             
              
              WRITE( 10, '(I9,6e23.15)',ADVANCE='NO') BemElementIndeces( CurrentElement % ElementIndex - &
                  Solver % Mesh % NumberOfBulkElements ), Load(1,1), Load(2,1), Load(3,1), Load(4,1), & 
                  Load(5,1), Load(6,1)
              WRITE( 10,* ) ''   
              
            END IF
          END IF
        END IF
      END DO
    END DO

    CLOSE(10)

  END IF


  !---------------------------------------------------------------------------------
  ! Compute and save surface force on boundary
  !---------------------------------------------------------------------------------
  TotalForce = 0.0d0
  TotalMoment = 0.0d0
  TotalArea = 0.0d0
  MomentOutput = .FALSE.

  DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
      Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements

    CurrentElement => Solver % Mesh % Elements(t)
    Model % CurrentElement => CurrentElement

    IF ( .NOT. ActiveBoundaryElement(CurrentElement, CurrentModel % Solver) ) CYCLE    

    !------------------------------------------------------------------------------
    ! Extract the parent element to find its material parameters... 
    !----------------------------------------------------------------------------- 
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    IF (ANY(FlowPerm(NodeIndexes(1:n)) == 0)) CYCLE
    
    DO i=1,Model % NumberOfBCs
      IF ( CurrentElement % BoundaryInfo % Constraint == &
          Model % BCs(i) % Tag ) THEN

        IF ( .NOT. ListGetLogical( Model % BCs(i) % Values, &
            'Calculate Fluidic Force', GotIt ) ) CYCLE

        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
        
        Parent => CurrentELement % BoundaryInfo % Left
        stat = ASSOCIATED( Parent )
        IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

        IF ( .NOT. stat) THEN
          Parent => CurrentELement % BoundaryInfo % Right
          stat = ASSOCIATED( Parent )
          IF (stat) stat = ALL( FlowPerm(Parent % NodeIndexes(1: Parent % TYPE % NumberOfNodes)) > 0 )

          IF ( .NOT. stat )  CALL Fatal( 'AcousticsSolver', &
              'No parent element can be found for given boundary element' )
        END IF
    
        np = Parent % TYPE % NumberOfNodes
        ParentNodes % x(1:np) = Solver % Mesh % Nodes % x(Parent % NodeIndexes)
        ParentNodes % y(1:np) = Solver % Mesh % Nodes % y(Parent % NodeIndexes)
        ParentNodes % z(1:np) = Solver % Mesh % Nodes % z(Parent % NodeIndexes)

        k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
            'Material' )
        Material => Model % Materials(k) % Values
        
        Viscosity(1:np) = ListGetReal( Material, 'Viscosity', &
            np, Parent % NodeIndexes )   
        Lambda(1:np) = -2.0d0/3.0d0 * Viscosity(1:np)
        BulkViscosity(1:np) = ListGetReal( Material, 'Bulk Viscosity', &
            np, Parent % NodeIndexes, GotIt )
        IF (GotIt) Lambda(1:np) = BulkViscosity(1:np) - 2.0d0/3.0d0 * Viscosity(1:np)
        
        MomentAbout(1) = ListGetConstReal( Model % BCs(i) % Values, &
            'Moment About 1', GotIt )

        MomentAbout(2) = ListGetConstReal( Model % BCs(i) % Values, &
            'Moment About 2', CalculateMoment )
        CalculateMoment = GotIt .OR. CalculateMoment
        
        MomentAbout(3) = ListGetConstReal( Model % BCs(i) % Values, &
            'Moment About 3', GotIt)
        CalculateMoment = GotIt .OR. CalculateMoment
        MomentOutput(i) = CalculateMoment

        DO j=1,np
          k = FlowPerm( Parent % NodeIndexes(j) ) 
          WallVelocity(1,j) =  Flow( (k-1)*(dim*2+4)+1 )
          WallVelocity(2,j) =  Flow( (k-1)*(dim*2+4)+2 )
          WallVelocity(3,j) =  Flow( (k-1)*(dim*2+4)+3 )
          WallVelocity(4,j) =  Flow( (k-1)*(dim*2+4)+4 )
          ! Load will contain the pressure solution...
          IF (dim > 2) THEN
            WallVelocity(5,j) =  Flow( (k-1)*(dim*2+4)+5 )
            WallVelocity(6,j) =  Flow( (k-1)*(dim*2+4)+6 )
            Load(1,j) =  Flow( (k-1)*(dim*2+4)+9 )
            Load(2,j) =  Flow( (k-1)*(dim*2+4)+10 )           
          ELSE
            WallVelocity(5,j) =  0.0d0
            WallVelocity(6,j) =  0.0d0
            Load(1,j) =  Flow( (k-1)*(dim*2+4)+7 )
            Load(2,j) =  Flow( (k-1)*(dim*2+4)+8 )                       
          END IF
        END DO
        
        Traction = 0.0d0
        Moment = 0.0d0
  
        CALL SurfaceForceIntegration(CurrentElement, Parent, Traction, Moment, &
            MomentAbout, Area, CalculateMoment, WallVelocity, Load, Viscosity, &
            Lambda, ElementNodes, ParentNodes, np)

        TotalForce(i,1:6) = TotalForce(i,1:6) + Traction(1:6)
        TotalArea(i) = TotalArea(i) + Area
        IF (CalculateMoment) THEN
          TotalMoment(i,1:6) = TotalMoment(i,1:6) + Moment(1:6)          
        END IF

      END IF
    END DO
  END DO



  CALL ListAddConstReal( Model % Simulation, 'res: frequency', &
      AngularFrequency / (2*PI) )

  DO k=1, Model % NumberOfBCs
    IF ( .NOT. ListGetLogical(Model % BCs(k) % Values,'Calculate Fluidic Force',stat ) ) CYCLE
    IF(Model % NumberOfBCs < 10) THEN
      WRITE( BoundaryName, '("bc ",I1)') k
    ELSE
      WRITE( BoundaryName, '("bc ",I2)') k
    END IF

    nlen = LEN_TRIM(BoundaryName)

    CALL Info('ForceCompute','Forces on Boundary '//BoundaryName(1:nlen),Level=4 )
    WRITE( Message, '("Real Fluidic Force (X,Y,Z):", 3ES17.6E2)') TotalForce(k,1:3)
    CALL Info( 'ForceCompute', Message, Level=4 )    
    WRITE( Message, '("Imaginary Fluidic Force (X,Y,Z):", 3ES17.6E2)') TotalForce(k,4:6)
    CALL Info( 'ForceCompute', Message, Level=4 ) 
    WRITE( Message, '("Contact Area:   ", ES17.6E2)') TotalArea(k)
    CALL Info( 'ForceCompute', Message, Level=4 )


    CALL ListAddConstReal( Model % Simulation, &
        'res: re contact force 1 '//BoundaryName(1:nlen), TotalForce(k,1) )
    CALL ListAddConstReal( Model % Simulation, &
        'res: re contact force 2 '//BoundaryName(1:nlen), TotalForce(k,2) )
    IF ( DIM > 2 )  CALL ListAddConstReal( Model % Simulation, &
        'res: re contact force 3 '//BoundaryName(1:nlen), TotalForce(k,3) )
    CALL ListAddConstReal( Model % Simulation, &
        'res: im contact force 1 '//BoundaryName(1:nlen), TotalForce(k,4) )
    CALL ListAddConstReal( Model % Simulation, &
        'res: im contact force 2 '//BoundaryName(1:nlen), TotalForce(k,5) )
    IF ( DIM > 2 )  CALL ListAddConstReal( Model % Simulation, &
        'res: im contact force 3 '//BoundaryName(1:nlen), TotalForce(k,6) )

    IF ( MomentOutput(k) ) THEN
      WRITE( Message, &
          '("Real Moment about (",ES9.3E1,",",ES9.3E1,",",ES9.3E1,"):",3Es14.6E2)') &
          MomentAbout(1:3), TotalMoment(k,1:3)
      CALL Info( 'ForceCompute', Message, Level=4 )
      WRITE( Message, &
          '("Imaginary Moment about (",ES9.3E1,",",ES9.3E1,",",ES9.3E1,"):",3Es14.6E2)') &
          MomentAbout(1:3), TotalMoment(k,4:6)
      CALL Info( 'ForceCompute', Message, Level=4 )
      CALL ListAddConstReal( Model % Simulation, &
          'res: re contact moment 1 '//BoundaryName(1:nlen), TotalMoment(k,1) )
      CALL ListAddConstReal( Model % Simulation, &
          'res: re contact moment 2 '//BoundaryName(1:nlen), TotalMoment(k,2) )
      CALL ListAddConstReal( Model % Simulation, &
          'res: re contact moment 3 '//BoundaryName(1:nlen), TotalMoment(k,3) )
      CALL ListAddConstReal( Model % Simulation, &
          'res: im contact moment 1 '//BoundaryName(1:nlen), TotalMoment(k,4) )
      CALL ListAddConstReal( Model % Simulation, &
          'res: im contact moment 2 '//BoundaryName(1:nlen), TotalMoment(k,5) )
      CALL ListAddConstReal( Model % Simulation, &
          'res: im contact moment 3 '//BoundaryName(1:nlen), TotalMoment(k,6) )
    END IF

    CALL ListAddConstReal( Model % Simulation, & 
        'res: contact force area '//BoundaryName(1:nlen), TotalArea(k) )

  END DO


  IF (.FALSE.) THEN
     ! print *, 'Contact area is = ', Area
     ZAppr = CMPLX(0.0d0, 1.0d0, kind=dp) * AngularFrequency * C1/C2
     ZAppr2 = CMPLX(0.0d0, 1.0d0, kind=dp) * AngularFrequency * C3/C4
  END IF
  
  IF (.FALSE.) THEN
  !if (ScanningUsed) then
     !if (FirstTimeStep) then
     !   CALL ListAddConstReal( Model % Simulation, & 
     !        'res: Re Iteration Coefficient', 0.0d0)
     !   CALL ListAddConstReal( Model % Simulation, & 
     !        'res: Im Iteration Coefficient', 0.0d0)        
     !else
     !   if (NewTimeStep) then
     !ReIterationCoeff = 0.0d0
     !ImIterationCoeff = 0.0d0
     ReIterationCoeff = REAL(ZAppr2)
     ImIterationCoeff = AIMAG(ZAppr2)
     !   end if
     CALL ListAddConstReal( Model % Simulation, & 
          'res: Re Iteration Coefficient', ReIterationCoeff)
     CALL ListAddConstReal( Model % Simulation, & 
          'res: Im Iteration Coefficient', ImIterationCoeff)
     !end if
  END IF
        
  
  !WRITE( Message, '("Iteration Coefficient:   ", 2ES17.6E2)') real(ZAppr), aimag(ZAppr)
  !CALL Info( 'IterationCoefficientCompute', Message, Level=4 )

  !WRITE( Message, '("Iteration Coefficient 2:   ", 2ES17.6E2)') real(ZAppr2), aimag(ZAppr2)
  !CALL Info( 'IterationCoefficientCompute', Message, Level=4 )



  !------------------------------------------------------------------------------
  ! Check if computing acoustic impedance required
  !------------------------------------------------------------------------------
  ALLOCATE( Bndries( Model % NumberOfBCs ) )
  Bndries = 0
  
  AcousticI = 0
  DO i = 1, Model % NumberOfBCs
    IF ( ListGetLogical( Model % BCs(i) % Values, &
        'Calculate Acoustic Impedance', GotIt ) ) THEN
      AcousticI = 1
      Bndries(1) = i
      EXIT
    END IF
  END DO

  IF ( AcousticI > 0 ) THEN
    j = 1
    DO i = 1, Model % NumberOfBCs
      IF ( ListGetLogical( Model % BCs(i) % Values, &
          'Impedance Target Boundary', GotIt ) ) THEN
        AcousticI = AcousticI + 1
        j = j + 1
        Bndries(j) = i
      END IF
    END DO

    ALLOCATE( AcImpedances( AcousticI, 2 ) )
    CALL ComputeAcousticImpedance( AcousticI, Bndries, AcImpedances )

    WRITE( Message, * ) 'Self specific acoustic impedance on bc ', &
        Bndries(1)
    CALL INFO( 'Acoustics', Message, LEVEL=5 )
    WRITE( Message, * ) '  In-phase with velocity:     ', &
        AcImpedances(1,1)
    CALL INFO( 'Acoustics', Message, LEVEL=5 )
    WRITE( Message, * ) '  Out-of-phase with velocity: ', &
        AcImpedances(1,2)
    CALL INFO( 'Acoustics', Message, LEVEL=5 )

    DO i = 2, AcousticI
      WRITE( Message, * ) 'Cross specific acoustic impedance, bcs ', &
          Bndries(i), ' , ', Bndries(1)
      CALL INFO( 'Acoustics', Message, LEVEL=5 )
      WRITE( Message, * ) '  In-phase with velocity:     ', &
          AcImpedances(i,1)
      CALL INFO( 'Acoustics', Message, LEVEL=5 )
      WRITE( Message, * ) '  Out-of-phase with velocity: ', &
          AcImpedances(i,2)
      CALL INFO( 'Acoustics', Message, LEVEL=5 )
    END DO

    DO i = 2, AcousticI
      WRITE( Message, '(A,I1,A,I1)' ) &
          'res: out-of-phase cross acoustic impedance ', &
          Bndries(AcousticI - i + 2), ' , ', Bndries(1)
      CALL ListAddConstReal( Model % Simulation, Message, &
          AcImpedances(AcousticI - i + 2,2) )
      WRITE( Message, '(A,I1,A,I1)' ) &
          'res: in-phase cross acoustic impedance ', &
          Bndries(AcousticI - i + 2), ' , ', Bndries(1)
      CALL ListAddConstReal( Model % Simulation, Message, &
          AcImpedances(AcousticI - i + 2,1) )
    END DO

    WRITE( Message, '(A,I1)' ) &
        'res: out-of-phase self acoustic impedance ', &
        Bndries(1)
    CALL ListAddConstReal( Model % Simulation, Message, &
        AcImpedances(1,2) )
    WRITE( Message, '(A,I1)' ) &
        'res: in-phase self acoustic impedance ', &
        Bndries(1)
    CALL ListAddConstReal( Model % Simulation, Message, &
        AcImpedances(1,1) )
    
    DEALLOCATE( AcImpedances )
  END IF

  DEALLOCATE( Bndries )

!------------------------------------------------------------------------------
!  CALL ListAddConstReal( Model % Simulation, 'res: frequency', &
!      AngularFrequency / (2*PI) )
!------------------------------------------------------------------------------


  IF ( BEMCoupling ) THEN
    OPEN( 10, FILE = ElmerStatusFile, status='REPLACE', IOSTAT = istat )
    IF ( istat /= 0) THEN
      CALL Fatal( 'AcousticsSolver', 'Cannot open Elmer Status File' )    
    ELSE
      CoupledTolerance = ListGetConstReal( Solver % Values, &
          'Steady State Convergence Tolerance' )
      MaxCoupledIterations =  ListGetInteger( Model % Simulation, &
          'Steady State Max Iterations' )
      IF ( (RelativeChange < CoupledTolerance) .OR. (SolverCalls >= MaxCoupledIterations) ) THEN
        WRITE( 10, '(A)', ADVANCE='NO') '$STOP'  
      ELSE  
        WRITE(str(1:3),'(I3)') SolverCalls
        WRITE( 10, '(A)', ADVANCE='NO') '$CONTINUE #' // ADJUSTL(str(1:3))
      END IF
    END IF
    CLOSE(10)
  END IF

  FirstVisit = .FALSE.
  CALL Info( 'AcousticsSolver', 'Exiting the solver...', Level=4 )  


CONTAINS



!------------------------------------------------------------------------------
  SUBROUTINE OptimalScaling( n, A, s )
!------------------------------------------------------------------------------
!   This subroutine equilibrates the rows of the complex-coefficient matrix A 
!   to minimize the condition number. Only the coefficients necessary for obtaining 
!   Re(A*z), with z a complex vector, are modified. The vector s will contain 
!   the complex scaling factors.
!-----------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(kind=dp) :: s(n)
!-----------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    !-------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp     
      s(i) = CMPLX( 1.0d0,0.0d0, kind=dp) / tmp
    END DO

!------------------------------------------------------------------------------
    END SUBROUTINE OptimalScaling
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE OptimalMatrixScaling( n, A, s )
!------------------------------------------------------------------------------
!   This subroutine equilibrates the rows of the complex-coefficient matrix A 
!   to minimize the condition number. Only the coefficients necessary for obtaining 
!   Re(A*z), with z a complex vector, are modified. The vector s will contain 
!   the complex scaling factors.
!-----------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(kind=dp) :: s(n)
!-----------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    !-------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp     
      s(i) = CMPLX( 1.0d0,0.0d0, kind=dp) / tmp
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        Values(j) = Values(j) * s(i)
        Values(j+1) = Values(j+1) * s(i)
      END DO
    END DO
    !--------------------------------------------------------------
    WRITE( Message, * ) 'Unscaled matrix norm: ', norm    
    CALL Info( 'AcousticsSolver', Message, Level=5 )
    !--------------------------------------------------------------
    norm = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
    END DO
    !---------------------------------------------------------------
    WRITE( Message, * ) 'Scaled matrix norm: ', norm
    CALL Info( 'AcousticsSolver', Message, Level=5 )
    !---------------------------------------------------------------

    ! recovering original system
    !DO i=1,n
    !  b(i) = b(i)/s(i)
    !  DO j=Rows(2*i-1),Rows(2*i)-1,2
    !    Values(j) = Values(j)/s(i) 
    !    Values(j+1) = Values(j+1)/s(i) 
    !  END DO
    !END DO


!------------------------------------------------------------------------------
    END SUBROUTINE OptimalMatrixScaling
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
  SUBROUTINE GCROuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!    This is the preconditioned GCR iteration for the complex linear system Ax=b.
!    The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------  
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp) :: x(n), b(n), TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, k, m, InnerRounds, IluOrder, VelocityPrecond = 1, &
        MaxRestarts        
    LOGICAL :: Condition, GotIt, SystemScaling, ConvergedSol = .FALSE.
    REAL(KIND=dp) :: res, tottime, res0, const, stime, InnerTol, alpha
    COMPLEX(KIND=dp) :: r(n/2),T1(n/2),T2(n/2), &
        S(n/2,Rounds), V(n/2,Rounds), y(n/2),f(n/2), e(n/2), &
        Vel(dim*q/2), VelRhs(dim*q/2), da(n/2), dps(q), dpa(dim*q/2)
    COMPLEX(KIND=dp) :: beta 
!------------------------------------------------------------------------------

    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number 
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)     
    IF (.NOT. GotIt) SystemScaling = .FALSE.
    IF (SystemScaling) THEN
      CALL OptimalMatrixScaling( n/2, A, da ) 
      CALL OptimalMatrixScaling( q, PS, dps )
      CALL OptimalMatrixScaling( dim*q/2, PA, dpa )   
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices. 
    ! This needs to be done only once. 
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement') 

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt) 

    IF (GotIt) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !-------------------------------------------------------------------------------- 
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations') 
    MaxRestarts = ListGetInteger( Solver % Values, &
        'Max GCR Restarts', GotIt) 
    IF (.NOT. GotIt) MaxRestarts = 1
 
    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution cannot be used as an initial guess
    !--------------------------------------------------------------------------------------------
    x(1:n) = 0.0d0 

    Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
    DO i=1,dim
      DO j=1,q/2
        VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
            b(2*(dim+2)*(j-1)+2*i), kind=dp )
      END DO
    END DO

    IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
      !-----------------------------------
      ! Scale the right-hand side vector 
      !-----------------------------------
      IF (SystemScaling) THEN 
        DO i=1,dim*q/2
          VelRhs(i) = dpa(i) * VelRhs(i)       
        END DO
      END IF
      !-----------------------------------------------------------------
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', &
          'Solving initial guess for velocities', Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      !------------------------------------------------------------------
      CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

      DO j=1,q/2
        DO i=1,dim
          x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
          x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
        END DO
      END DO
    END IF
 
    !--------------------------------------------------------------------------------
    ! The start of the GCR iteration... 
    !--------------------------------------------------------------------------------     
    m = n/2
    !----------------------------------------------------------------------------
    ! Transform the solution vector x and the right-hand side vector b to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------
    DO i=1,m
      y(i) = CMPLX( x(2*i-1), x(2*i), kind=dp )
      f(i) = CMPLX( b(2*i-1), b(2*i), kind=dp )
    END DO

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)       
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    ! PRINT *,'OuterIteration ',0, ComplexNorm(m,r)/res0, StoppingCriterion(m,A,y,f,r), &
    !     CPUTime() - tottime        
    WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, CPUTime() - tottime 

    DO j=1,MaxRestarts
      IF (ConvergedSol) EXIT
      V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      DO k=1,Rounds
        !----------------------------------------------------------
        ! Perform the preconditioning...
        !---------------------------------------------------------------
        T1(1:m) = r(1:m)
        CALL PreconditioningIteration(n, A, q, PA, PS, T1, dim, &
            da, dpa, dps, SystemScaling, MMatrix, LMatrix)
        CALL ComplexMatrixVectorProduct( A, T1, T2 )  
      
        !--------------------------------------------------------------
        ! Perform the orthogonalisation of the search directions...
        !--------------------------------------------------------------
        DO i=1,k-1
          beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
          T1(1:m) = T1(1:m) - beta * S(1:m,i)
          T2(1:m) = T2(1:m) - beta * V(1:m,i)        
        END DO
        alpha = ComplexNorm(m,T2)
        T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
        T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

        !-------------------------------------------------------------
        ! The update of the solution and save the search data...
        !------------------------------------------------------------- 
        beta = ComplexDotProduct(m, T2, r)
        y(1:m) = y(1:m) + beta * T1(1:m)      
        r(1:m) = r(1:m) - beta * T2(1:m)
        S(1:m,k) = T1(1:m)
        V(1:m,k) = T2(1:m) 

        !----------------------------------------------------------------
        ! Check the accuracy of the residual, if desired...
        !----------------------------------------------------------------
        IF (.FALSE.) THEN
          e(1:m) = CMPLX( 0.0d0,0.0d0, kind=dp )
          CALL ComplexMatrixVectorProduct( A, y, e )
          e(1:m) = f(1:m) - e(1:m)
          norm = ComplexNorm(m,e(1:m)-r(1:m))/ComplexNorm(m,e(1:m))
          PRINT *, 'Relative error of the residual: ', norm
        END IF

        !---------------------------------------------------- 
        ! Check whether the convergence criterion is met 
        !----------------------------------------------------
        res = ComplexNorm(m,r)/res0

        !res = StoppingCriterion( m, A, y, f, r ) 
        !PRINT *,'OuterIteration ',i,res, StoppingCriterion(m,A,y,f,r),CPUTime() - tottime

        WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
            k + (j-1) * Rounds, res, CPUTime() - tottime 
        ConvergedSol = ( res < TOL)
        IF (ConvergedSol) EXIT      
      END DO
    END DO

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
        ConditionEstimate( m, A, y, f)

    !----------------------------------------------
    ! Return the solution as a real vector...
    !----------------------------------------------
    DO i=1,m
      x( 2*i-1 ) = REAL( y(i) )
      x( 2*i ) = AIMAG( y(i) )
    END DO
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))  

!------------------------------------------------------------------------------
  END SUBROUTINE GCROuterIteration
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE BiCGStabOuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!   This is the preconditioned BiCGStab iteration for the complex linear 
!   system Ax=b. The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------  
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp) :: x(n), b(n), TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, m, InnerRounds, IluOrder, VelocityPrecond = 1
    LOGICAL :: Condition, GotIt, SystemScaling
    REAL(KIND=dp) :: res, tottime, res0, InnerTol
    COMPLEX(KIND=dp) :: r(n/2),Ri(n/2),P(n/2),V(n/2),T(n/2),T1(n/2),T2(n/2),S(n/2), &
        y(n/2), f(n/2), Vel(dim*q/2), VelRhs(dim*q/2), da(n/2), dps(q), dpa(dim*q/2)
    COMPLEX(KIND=dp) :: alpha, beta, omega, rho, oldrho
!------------------------------------------------------------------------------

    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number 
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)     
    IF (.NOT. GotIt) SystemScaling = .FALSE.
    IF (SystemScaling) THEN
      CALL OptimalMatrixScaling( n/2, A, da ) 
      CALL OptimalMatrixScaling( q, PS, dps )
      CALL OptimalMatrixScaling( dim*q/2, PA, dpa )   
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices. 
    ! This needs to be done only once. 
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement') 

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt) 

    IF (GotIt) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !-------------------------------------------------------------------------------- 
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations') 

    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution cannot be used as an initial guess
    !--------------------------------------------------------------------------------------------
    x(1:n) = 0.0d0 

    Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
    DO i=1,dim
      DO j=1,q/2
        VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
            b(2*(dim+2)*(j-1)+2*i), kind=dp )
      END DO
    END DO

    IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
      !-----------------------------------
      ! Scale the right-hand side vector 
      !-----------------------------------
      IF (SystemScaling) THEN 
        DO i=1,dim*q/2
          VelRhs(i) = dpa(i) * VelRhs(i)       
        END DO
      END IF
      !-----------------------------------------------------------------
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', &
          'Solving initial guess for velocities', Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      !------------------------------------------------------------------
      CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

      DO j=1,q/2
        DO i=1,dim
          x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
          x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
        END DO
      END DO
    END IF
 
    !--------------------------------------------------------------------------------
    ! The start of the BiCGStab iteration... 
    !--------------------------------------------------------------------------------     
    m = n/2
    !----------------------------------------------------------------------------
    ! Transform the solution vector x and the right-hand side vector b to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------
    DO i=1,m
      y(i) = CMPLX( x(2*i-1), x(2*i), kind=dp )
      f(i) = CMPLX( b(2*i-1), b(2*i), kind=dp )
    END DO

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)       
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    ! PRINT *,'OuterIteration ',0, ComplexNorm(m,r)/res0, StoppingCriterion(m,A,y,f,r), &
    !     CPUTime() - tottime
        
    WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, CPUTime() - tottime 

    Ri(1:m) = r(1:m)
    P(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    V(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    omega  = CMPLX( 1.0d0, 0.0d0, kind=dp)
    alpha  = CMPLX( 0.0d0, 0.0d0, kind=dp)
    oldrho = CMPLX( 1.0d0, 0.0d0, kind=dp)

    DO i=1,Rounds
      rho = ComplexDotProduct( m, Ri, r )
      beta = alpha * rho / ( oldrho * omega )
      P(1:m) = r(1:m) + beta * (P(1:m) - omega*V(1:m))
      V(1:m) = P(1:m)
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !---------------------------------------------------------------
      CALL PreconditioningIteration(n, A, q, PA, PS, V, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      !--------------------------------------------------------------
      T1(1:m) = V(1:m)
      CALL ComplexMatrixVectorProduct( A, T1, V )
      alpha = rho / ComplexDotProduct( m, Ri, V )
      S(1:m) = r(1:m) - alpha * V(1:m)
      
      !---------------------------------------------------------------------------------
      ! The update of the solution and the computation of the residual-based error indicator  
      !---------------------------------------------------------------------------------
      y(1:m) = y(1:m) + alpha*T1(1:m)

      res = ComplexNorm(m,S)/res0
      !res = StoppingCriterion( m, A, y, f, S )
      IF ( res < TOL ) THEN
        WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', i, res, &
            StoppingCriterion( m, A, y, f, S )
        WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
            ConditionEstimate( m, A, y, f)
        EXIT
      END IF

      T(1:m) = S(1:m)
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !-----------------------------------------------------------------         
      CALL PreconditioningIteration(n, A, q, PA, PS, T, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      !-----------------------------------------------------------------
      T2(1:m) = T(1:m)
      CALL ComplexMatrixVectorProduct( A, T2, T )
      omega = ComplexDotProduct( m,T,S ) / ComplexDotProduct( m,T,T )
      oldrho = rho
      r(1:m) = S(1:m) - omega*T(1:m)
      y(1:m) = y(1:m) + omega*T2(1:m)

      res = ComplexNorm(m,r)/res0
      ! res = StoppingCriterion( m, A, y, f, r ) 

      WRITE(*,'(a,I4,ES12.3,ES12.3)') 'OuterIteration residual for iterate', i, res, &
          StoppingCriterion( m, A, y, f, r )

      IF ( res < TOL .OR. i==Rounds) THEN
        WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
            ConditionEstimate( m, A, y, f)
        EXIT
      END IF
    END DO

    !----------------------------------------------
    ! Return the solution as a real vector...
    !----------------------------------------------
    DO i=1,m
      x( 2*i-1 ) = REAL( y(i) )
      x( 2*i ) = AIMAG( y(i) )
    END DO
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))  


!------------------------------------------------------------------------------
  END SUBROUTINE BiCGStabOuterIteration
!------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------
  SUBROUTINE BiCGStablOuterIteration( l, m, A, q, PA, PS, MMatrix, LMatrix, &
      v, f, MaxRounds, Tol, dim, Norm)
!-----------------------------------------------------------------------------------
!  This is the preconditioned BiCGStab(l) iteration for the complex linear 
!  system Ax=b. The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------  
   INTEGER :: l, m, q, MaxRounds, dim   
   TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix  
   REAL(KIND=dp) :: v(m), f(m), Tol, Norm
!----------------------------------------------------------------------------------- 
   LOGICAL :: SystemScaling, GotIt, Condition
   COMPLEX(KIND=dp) :: da(m/2), dps(q), dpa(dim*q/2), Vel(dim*q/2), VelRhs(dim*q/2), &
       x(m/2), b(m/2)
   INTEGER :: IluOrder, InnerRounds, VelocityPrecond = 1
   REAL(KIND=dp) :: InnerTol, tottime

   COMPLEX(KIND=dp) :: zzero, zone, t(m/2), kappa0, kappal
   REAL(KIND=dp) :: dznrm2, rnrm0, bnrm, rnrm, mxnrmx, mxnrmr, errorind, &
       delta = 1.0d-2
   INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round
   COMPLEX(KIND=dp) :: work(m/2,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
       alpha, beta, omega, rho0, rho1, sigma, zdotc, varrho, hatgamma
   LOGICAL rcmp, xpdt 
!------------------------------------------------------------------------------

   tottime = CPUTime()
   !-----------------------------------------------------------------------------
   ! Scaling to minimize the condition number 
   !-----------------------------------------------------------------------------
   SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)     
   IF (.NOT. GotIt) SystemScaling = .FALSE.
   IF (SystemScaling) THEN
     CALL OptimalMatrixScaling( m/2, A, da ) 
     CALL OptimalMatrixScaling( q, PS, dps )
     CALL OptimalMatrixScaling( dim*q/2, PA, dpa )   
   END IF
 
   !-------------------------------------------------------------------------------
   ! Compute ILU factorizations for the preconditioner matrices. 
   ! This needs to be done only once. 
   !-------------------------------------------------------------------------------
   IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement') 

   CALL Info( 'AcousticsSolver', ' ', Level=4)
   CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
       Level=4)
   CALL Info( 'AcousticsSolver', ' ', Level=4)
   Condition = CRS_ComplexIncompleteLU( PS, IluOrder )

   IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt) 

   IF (GotIt) THEN
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
         Level=4)
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
   ELSE
     VelocityPrecond = 0
   END IF

   !--------------------------------------------------------------------------------
   !   Some initializations
   !-------------------------------------------------------------------------------- 
   InnerTol = ListGetConstReal( Solver % Values, &
       'Linear System Convergence Tolerance' )
   InnerRounds = ListGetInteger( Solver % Values, &
       'Linear System Max Iterations') 

   !--------------------------------------------------------------------------------------------
   ! The solution of an initial guess, the previous solution cannot be used as an initial guess
   !--------------------------------------------------------------------------------------------
   v(1:m) = 0.0d0

   Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
   DO i=1,dim
     DO j=1,q/2
       VelRhs((j-1)*dim+i) = CMPLX( f(2*(dim+2)*(j-1)+2*i-1), &
           f(2*(dim+2)*(j-1)+2*i), kind=dp )
     END DO
   END DO

   IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN   
     !-----------------------------------
     ! Scale the right-hand side vector 
     !-----------------------------------
     IF (SystemScaling) THEN
       DO i=1,dim*q/2
         VelRhs(i) = dpa(i) * VelRhs(i)
       END DO
     END IF
     !-----------------------------------------------------------------
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     CALL Info( 'AcousticsSolver', &
         'Solving initial guess for velocities', Level=4)
     CALL Info( 'AcousticsSolver', ' ', Level=4)
     !------------------------------------------------------------------
     CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )

     DO j=1,q/2
       DO i=1,dim
         v(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
         v(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
       END DO
     END DO
   END IF

   !--------------------------------------------------------------------------------
   ! The start of the BiCGStabl iteration... 
   !--------------------------------------------------------------------------------     
   n = m/2
   !----------------------------------------------------------------------------
   ! Transform the solution vector v and the right-hand side vector f to 
   ! complex-valued vectors x and b
   !---------------------------------------------------------------------------
   DO i=1,n
     x(i) = CMPLX( v(2*i-1), v(2*i), kind=dp )
     b(i) = CMPLX( f(2*i-1), f(2*i), kind=dp )
   END DO

   zzero = CMPLX( 0.0d0,0.0d0, kind=dp)
   zone =  CMPLX( 1.0d0,0.0d0, kind=dp)
   work = CMPLX( 0.0d0, 0.0d0, kind=dp )
   rwork = CMPLX( 0.0d0, 0.0d0, kind=dp )
    
   rr = 1
   r = rr+1
   u = r+(l+1)
   xp = u+(l+1)
   bp = xp+1

   z = 1
   zz = z+(l+1)
   y0 = zz+(l+1)
   yl = y0+1
   y = yl+1

   CALL ComplexMatrixVectorProduct( A, x, work(1:n,r) )
   work(1:n,r) = b(1:n) - work(1:n,r)
   bnrm = dznrm2(n, b(1:n), 1)
    
   work(1:n,rr) = work(1:n,r) 
   work(1:n,bp) = work(1:n,r)
   work(1:n,xp) = x(1:n)
   x(1:n) = zzero
   rnrm0 = dznrm2(n, work(1:n,r), 1)
   rnrm = rnrm0
   mxnrmx = rnrm0
   mxnrmr = rnrm0
    
   alpha = zzero
   omega = zone
   sigma = zone
   rho0 = zone

   Round = 0
   errorind = rnrm/bnrm

   WRITE(*,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', &
        0, errorind 

   DO WHILE ( errorind > Tol .AND. Round < MaxRounds) 
     Round = Round + 1 
     !-------------------------
     ! --- The BiCG part ---
     !-------------------------
     rho0 = -omega*rho0

     DO k=1,l
       rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
       IF (rho0 == zzero) THEN
         CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
       ENDIF
       beta = alpha*(rho1/rho0)
       rho0 = rho1
       DO j=0,k-1
         work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
       ENDDO
       t(1:n) = work(1:n,u+k-1)
       !----------------------------------------------------------
       ! Perform the preconditioning...
       !---------------------------------------------------------------
       CALL PreconditioningIteration(m, A, q, PA, PS, t, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
       !--------------------------------------------------------------
       
       CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) )      
       sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
       IF (sigma == zzero) THEN
         CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
       ENDIF
       alpha = rho1/sigma
       x(1:n) = x(1:n) + alpha * work(1:n,u)
       DO j=0,k-1
         work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
       ENDDO
       t(1:n) = work(1:n,r+k-1)
       !-----------------------------------------------------------------
       ! Preconditioning...
       !-----------------------------------------------------------------
       CALL PreconditioningIteration(m, A, q, PA, PS, t, dim, &
           da, dpa, dps, SystemScaling, MMatrix, LMatrix)
       !------------------------------------------------------------------
       CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )  
       rnrm = dznrm2(n, work(1:n,r), 1)
       mxnrmx = MAX (mxnrmx, rnrm)
       mxnrmr = MAX (mxnrmr, rnrm)
     ENDDO

     !--------------------------------------
     ! --- The convex polynomial part ---
     !--------------------------------------

     DO i=1,l+1
       DO j=1,i
         rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 ) 
       END DO
     END DO
     DO j=2,l+1
       rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
     END DO

     rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
     CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
         iwork, stat)

     ! --- tilde r0 and tilde rl (small vectors)
     
     rwork(1,y0) = -zone
     rwork(2:l,y0) = rwork(2:l,z) 
     CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
         rwork(2:l,y0), l-1, stat)
     rwork(l+1,y0) = zzero

     rwork(1,yl) = zzero
     rwork(2:l,yl) = rwork(2:l,z+l) 
     CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
         rwork(2:l,yl), l-1, stat)
     rwork(l+1,yl) = -zone

     ! --- Convex combination

     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
     kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
         (kappa0*kappal)
     hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
         kappa0/kappal
     rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

     !  --- Update
     
     omega = rwork(l+1,y0)
     DO j=1,l
       work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
       x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
       work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
     ENDDO


     CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
         rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
     rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------
     IF (.FALSE.) THEN
       IF (.TRUE.) THEN
         mxnrmx = MAX (mxnrmx, rnrm)
         mxnrmr = MAX (mxnrmr, rnrm)
         xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
         rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
         IF (rcmp) THEN
           PRINT *, 'Performing residual update...'
           t(1:n) = x(1:n)
           CALL CRS_ComplexLUSolve2( n, A, t )         
           CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
           work(1:n,r) = work(1:n,bp) - work(1:n,r)
           mxnrmr = rnrm
           IF (xpdt) THEN
             PRINT *, 'Performing solution update...'
             !CALL CRS_ComplexLUSolve2( n, A, x )
             work(1:n,xp) = work(1:n,xp) + t(1:n)
             x(1:n) = zzero
             work(1:n,bp) = work(1:n,r)
             mxnrmx = rnrm
           ENDIF
         ENDIF
       ENDIF

       IF (rcmp) THEN
         IF (xpdt) THEN       
           t(1:n) = work(1:n,xp)
         ELSE
           t(1:n) = t(1:n) + work(1:n,xp)  
         END IF
       ELSE
         t(1:n) = x(1:n)
         CALL CRS_ComplexLUSolve2( n, A, t ) 
         t(1:n) =  t(1:n) + work(1:n,xp)
       END IF
     END IF

     !errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )

     errorind = rnrm/bnrm
     WRITE(*,'(a,I4,ES12.3)') 'OuterIteration residual for iterate', Round, errorind

   END DO

   !------------------------------------------------------------
   ! We have solved z = P*x, so finally solve the true unknown x
   !------------------------------------------------------------
   !
   ! This is a nontrivial computation ...
   !
   !CALL CRS_ComplexLUSolve2( n, A, x )
   !x(1:n) = x(1:n) + work(1:n,xp)
   STOP

   WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
       ConditionEstimate( m, A, x, b)
   !----------------------------------------------
   ! Return the solution as a real vector...
   !----------------------------------------------
   DO i=1,n
     v( 2*i-1 ) = REAL( x(i) )
     v( 2*i ) = AIMAG( x(i) )
   END DO

   Norm = SQRT(DOT_PRODUCT( v(1:2*n), v(1:2*n) )/(2*n))

!------------------------------------------------------------------------------
   END SUBROUTINE BiCGStablOuterIteration
!------------------------------------------------------------------------------









!------------------------------------------------------------------------------
  SUBROUTINE InnerOuterIteration( n, A, q, PA, PS, MMatrix, LMatrix, &
      x, b, Rounds, TOL, dim, Norm )
!------------------------------------------------------------------------------
!   This is the nested GCR iteration for the complex linear system Ax=b.
!   The new search direction is solved from the residual equation As = r. 
!   The residual equations is solved using the block-preconditioned GCR(m) method
!   ideally with a small m.
!------------------------------------------------------------------------------  
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp) :: x(n), b(n), TOL, Norm
!-------------------------------------------------------------------------------
    INTEGER :: i, j, k, m, InnerRounds, IluOrder, VelocityPrecond = 1, InnerRestart, &
        MaxRestarts
    LOGICAL :: Condition, GotIt, SystemScaling, Truncation, &
        ConvergedSol, FirstVisit = .TRUE. 
    REAL(KIND=dp) :: res, tottime, res0, InnerTol, alpha, &
        ResidualReductionRatio, bw_error
    COMPLEX(KIND=dp) :: r(n/2),T1(n/2),T2(n/2), &
        S(n/2,Rounds), V(n/2,Rounds), y(n/2), f(n/2), &
        Vel(dim*q/2), VelRhs(dim*q/2), Sol(n/2), da(n/2), dps(q), dpa(dim*q/2)
    COMPLEX(KIND=dp) :: beta 
    SAVE FirstVisit
!------------------------------------------------------------------------------
    ConvergedSol = .FALSE.
    tottime = CPUTime()
    !-----------------------------------------------------------------------------
    ! Scaling to minimize the condition number 
    !-----------------------------------------------------------------------------
    SystemScaling = ListGetLogical( Solver % Values, 'Linear System Scaling', GotIt)     

    IF (.NOT. GotIt) SystemScaling = .TRUE.
    IF (SystemScaling) THEN
      WRITE( Message, * ) 'Scaling the system matrix...'
      CALL Info( 'AcousticsSolver', Message, Level=5 )
      CALL OptimalMatrixScaling( n/2, A, da ) 
      WRITE( Message, * ) 'Scaling the preconditioning matrix for Schur complement...'
      CALL Info( 'AcousticsSolver', Message, Level=5 )
      CALL OptimalMatrixScaling( q, PS, dps )
      IF ( ASSOCIATED(PA) ) THEN
        WRITE( Message, * ) 'Scaling the preconditioning matrix for velocities...'    
        CALL Info( 'AcousticsSolver', Message, Level=5 )
        CALL OptimalMatrixScaling( dim*q/2, PA, dpa )
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! Compute ILU factorizations for the preconditioner matrices. 
    ! This needs to be done only once. 
    !-------------------------------------------------------------------------------
    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Schur Complement') 

    CALL Info( 'AcousticsSolver', ' ', Level=4)
    CALL Info( 'AcousticsSolver', 'ILU factorization for the Schur complement preconditioner', &
        Level=4)
    CALL Info( 'AcousticsSolver', ' ', Level=4)
    Condition = CRS_ComplexIncompleteLU( PS, IluOrder )
    !Condition = CRS_ComplexILUT(PS,1.0d-2)

    IluOrder = ListGetInteger( Solver % Values, 'ILU Order for Velocities', GotIt) 

    IF (GotIt .AND. ASSOCIATED(PA) ) THEN
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      CALL Info( 'AcousticsSolver', 'ILU factorization for the Velocity preconditioner', &
          Level=4)
      CALL Info( 'AcousticsSolver', ' ', Level=4)
      Condition = CRS_ComplexIncompleteLU( PA, IluOrder )
    ELSE
      VelocityPrecond = 0
    END IF

    !--------------------------------------------------------------------------------
    !   Some initializations
    !-------------------------------------------------------------------------------- 
    InnerTol = ListGetConstReal( Solver % Values, &
        'Linear System Convergence Tolerance' )
    InnerRounds = ListGetInteger( Solver % Values, &
        'Linear System Max Iterations') 
    MaxRestarts = ListGetInteger( Solver % Values, &
        'Max Outer GCR Cycles', GotIt)
    IF (.NOT. GotIt) MaxRestarts = 1   
    InnerRestart = ListGetInteger( Solver % Values, 'Max Inner GCR Iterations', GotIt)
    IF (.NOT. GotIt) InnerRestart = 5
    ResidualReductionRatio = ListGetConstReal( Solver % Values, 'Residual Reduction Ratio', GotIt)
    IF (.NOT. GotIt) ResidualReductionRatio = 1.0d-1
    Truncation = ListGetLogical( Solver % Values, &
        'Use Truncation', GotIt)
    IF (.NOT. GotIt) Truncation = .FALSE.

    

    !--------------------------------------------------------------------------------------------
    ! The solution of an initial guess, the previous solution can be used as an initial guess
    !--------------------------------------------------------------------------------------------
    IF ( FirstVisit .OR. (.NOT. UtilizePreviousSolution) ) THEN
      x(1:n) = 0.0d0 

      Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      DO i=1,dim
        DO j=1,q/2
          VelRhs((j-1)*dim+i) = CMPLX( b(2*(dim+2)*(j-1)+2*i-1), &
              b(2*(dim+2)*(j-1)+2*i), kind=dp )
        END DO
      END DO

      IF ( ANY( VelRhs /= CMPLX(0.0d0, 0.0d0, kind=dp) ) ) THEN
        IF ( ASSOCIATED(PA) ) THEN 
          !-----------------------------------
          ! Scale the right-hand side vector 
          !-----------------------------------
          IF (SystemScaling) THEN 
            DO i=1,dim*q/2
              VelRhs(i) = dpa(i) * VelRhs(i)       
            END DO
          END IF
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Solving initial guess for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          !CALL ComplexBiCGStab( q*dim, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0 )      
          CALL ComplexBiCGStabl( 4, q*dim/2, PA, Vel, VelRhs, InnerRounds, 1.0d-7, 0)
        ELSE
          !-----------------------------------
          ! Scale the right-hand side vector 
          !-----------------------------------
          IF (SystemScaling) THEN
            DO i=1,dim
              DO j=1,q/2
                VelRhs((j-1)*dim+i) = da((j-1)*(dim+2)+i) * VelRhs((j-1)*dim+i)   
              END DO
            END DO
          END IF
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Solving initial guess for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          CALL VelocitySolve( 4, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, 1.0d-7, 0)         
        END IF

        DO j=1,q/2
          DO i=1,dim
            x(2*(dim+2)*(j-1)+2*i-1) = REAL(Vel((j-1)*dim+i))
            x(2*(dim+2)*(j-1)+2*i) = AIMAG(Vel((j-1)*dim+i))
          END DO
        END DO
      END IF
      FirstVisit = .FALSE.
    END IF
    
    !--------------------------------------------------------------------------------
    ! The start of the GCR iteration... 
    !--------------------------------------------------------------------------------     
    m = n/2
    !----------------------------------------------------------------------------
    ! Transform the solution vector x and the right-hand side vector b to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------
    DO i=1,m
      y(i) = CMPLX( x(2*i-1), x(2*i), kind=dp )
      f(i) = CMPLX( b(2*i-1), b(2*i), kind=dp )
    END DO

    CALL ComplexMatrixVectorProduct( A, y, r )
    !--------------------------------------------------------
    ! Scale the right-hand side vector and form the residual
    !--------------------------------------------------------
    IF (SystemScaling) THEN
      DO i=1,m
        f(i) = da(i) * f(i)       
      END DO
    END IF
    r(1:m) = f(1:m) - r(1:m)
    res0 = ComplexNorm(m,f)
    res = ComplexNorm(m,r)/res0

    WRITE(*,'(a,I4,ES12.3,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
        0, res, StoppingCriterion(m,A,y,f,r), CPUTime() - tottime 

    DO j=1,MaxRestarts
      IF (ConvergedSol) EXIT
      V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      DO k=1,Rounds
        !----------------------------------------------------------
        ! Perform the preconditioning...
        !---------------------------------------------------------------
        T1(1:m) = r(1:m)
        Sol(1:m) = CMPLX(0.0d0, 0.0d0, kind=dp)
        CALL InnerIteration( n, A, q, PA, PS, Sol, T1, &
            InnerRestart, ResidualReductionRatio, dim, MMatrix, LMatrix, &
            SystemScaling, da, dpa, dps)
        T1(1:m) = Sol(1:m)
        CALL ComplexMatrixVectorProduct( A, T1, T2 )  
      
        !--------------------------------------------------------------
        ! Perform the orthogonalisation of the search directions....
        !--------------------------------------------------------------
        DO i=1,k-1
          beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
          T1(1:m) = T1(1:m) - beta * S(1:m,i)
          T2(1:m) = T2(1:m) - beta * V(1:m,i)        
        END DO
        alpha = ComplexNorm(m,T2)
        T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
        T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

        !-------------------------------------------------------------
        ! The update of the solution and save the search data...
        !------------------------------------------------------------- 
        beta = ComplexDotProduct(m, T2, r)
        y(1:m) = y(1:m) + beta * T1(1:m)      
        r(1:m) = r(1:m) - beta * T2(1:m)
        S(1:m,k) = T1(1:m)
        V(1:m,k) = T2(1:m) 

        !--------------------------------------------------------------
        ! Check whether the convergence criterion is met 
        !--------------------------------------------------------------
        res = ComplexNorm(m,r)/res0

        bw_error = StoppingCriterion( m, A, y, f, r ) 
        !PRINT *,'OuterIteration residual',i,res, StoppingCriterion(m,A,y,f,r),CPUTime() - tottime

        WRITE(*,'(a,I4,ES12.3,ES12.3,ES12.3)') 'OuterIteration residual for iterate', &
            k + (j-1) * Rounds, res, bw_error, CPUTime() - tottime 
        ConvergedSol = ( bw_error < TOL)
        IF (ConvergedSol) EXIT              
      END DO
    END DO

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number: ', &
        ConditionEstimate( m, A, y, f )

    !----------------------------------------------
    ! Return the solution as a real vector...
    !----------------------------------------------
    DO i=1,m
      x( 2*i-1 ) = REAL( y(i) )
      x( 2*i ) = AIMAG( y(i) )
    END DO
    Norm = SQRT(DOT_PRODUCT( x(1:2*m), x(1:2*m) )/(2*m))          

!------------------------------------------------------------------------------
   END SUBROUTINE InnerOuterIteration
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE InnerIteration( n, A, q, PA, PS, x, b, Rounds, TOL, dim, &
      MMatrix, LMatrix, SystemScaling, da, dpa, dps)
!------------------------------------------------------------------------------
!   This is the preconditioned GCR(Rounds) iteration for the complex linear system Ax=b.
!   The Schur complement preconditioning strategy is used.
!------------------------------------------------------------------------------  
    TYPE(Matrix_t), POINTER :: A, PA, PS, MMatrix, LMatrix
    INTEGER :: n, q, Rounds, dim
    REAL(KIND=dp) :: TOL
    COMPLEX(KIND=dp) :: x(n/2), b(n/2)
    LOGICAL, OPTIONAL :: SystemScaling 
    COMPLEX(KIND=dp), OPTIONAL :: da(n/2), dpa(dim*q/2), dps(q)
!-------------------------------------------------------------------------------
    INTEGER :: k, m
    REAL(KIND=dp) :: res, tottime, res0, alpha
    COMPLEX(KIND=dp) :: r(n/2), T1(n/2), T2(n/2), S(n/2,Rounds), V(n/2,Rounds) 
    COMPLEX(KIND=dp) :: beta 
!------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!   The start of the GCR iteration... 
!--------------------------------------------------------------------------------     
    tottime = CPUTime()
    m = n/2

    CALL ComplexMatrixVectorProduct( A, x, r )
    r(1:m) = b(1:m) - r(1:m)
    res0 = ComplexNorm(m,r)

    V(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)
    S(1:m,1:Rounds) = CMPLX( 0.0d0, 0.0d0, kind=dp)

    DO k=1,Rounds
      !----------------------------------------------------------
      ! Perform the preconditioning...
      !---------------------------------------------------------------
      T1(1:m) = r(1:m)
      CALL PreconditioningIteration( n, A, q, PA, PS, T1, dim, &
          da, dpa, dps, SystemScaling, MMatrix, LMatrix)
      CALL ComplexMatrixVectorProduct( A, T1, T2 )  
      
      !--------------------------------------------------------------
      ! Perform the orthogonalisation of the search directions....
      !--------------------------------------------------------------
      DO i=1,k-1
        beta = ComplexDotProduct( m, V(1:m,i), T2(1:m) )
        T1(1:m) = T1(1:m) - beta * S(1:m,i)
        T2(1:m) = T2(1:m) - beta * V(1:m,i)        
      END DO
      alpha = ComplexNorm(m,T2)
      T1(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T1(1:m)
      T2(1:m) = CMPLX( 1.0d0, 0.0d0, kind=dp)/CMPLX( alpha, 0.0d0, kind=dp) * T2(1:m)

      !-------------------------------------------------------------
      ! The update of the solution and save the search data...
      !------------------------------------------------------------- 
      beta = ComplexDotProduct(m, T2, r)
      x(1:m) = x(1:m) + beta * T1(1:m)      
      r(1:m) = r(1:m) - beta * T2(1:m)
      S(1:m,k) = T1(1:m)
      V(1:m,k) = T2(1:m) 

      !--------------------------------------------------------------
      ! Check whether the convergence criterion is met 
      !--------------------------------------------------------------
      res = ComplexNorm(m,r)/res0
      WRITE(*,'(a,I4,ES12.3,ES12.3)') 'InnerIteration residual for iterate', &
          k, res, CPUTime() - tottime 
      IF ( res < TOL) EXIT
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE InnerIteration
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
     SUBROUTINE PreconditioningIteration( n, A, q, PA, PS, V, dim, &
         da, dpa, dps, Scaling, MMatrix, LMatrix)
!------------------------------------------------------------------------------
!     This subroutine solves iteratively the upper triangular preconditioning 
!     system P*z = V. The vector V is overwritten by the solution z.
!------------------------------------------------------------------------------  
      TYPE(Matrix_t), POINTER :: A, PA, PS
      INTEGER :: n, q, dim
      COMPLEX(KIND=dp) :: V(n/2)
      COMPLEX(KIND=dp), OPTIONAL :: da(n/2), dpa(dim*q/2), dps(q)
      LOGICAL, OPTIONAL :: Scaling
      TYPE(Matrix_t), POINTER, OPTIONAL :: MMatrix, LMatrix    
!--------------------------------------------------------------------------------
      REAL(kind=dp) :: InnerTol, VelocityTol, SchurTol
      INTEGER :: i, j, InnerRounds
      LOGICAL :: SystemScaling, ConsistentSplitting, VelocityCriterion, &
          SchurCriterion, GotIt
      COMPLEX(kind=dp) :: z(n/2), y(q/2), f(q/2), w(q), g(q), & 
          Vel(dim*q/2), VelRhs(dim*q/2)
!------------------------------------------------------------------------------
      InnerTol = ListGetConstReal( Solver % Values, &
          'Linear System Convergence Tolerance', GotIt )
      InnerRounds = ListGetInteger( Solver % Values, &
          'Linear System Max Iterations', GotIt)
      VelocityTol =  ListGetConstReal( Solver % Values, &
          'Velocity Convergence Tolerance', VelocityCriterion)
      SchurTol = ListGetConstReal( Solver % Values, &
          'Schur Complement Convergence Tolerance', SchurCriterion)

      ConsistentSplitting = PRESENT(MMatrix) .AND. PRESENT(LMatrix)

      IF ( PRESENT(Scaling) ) THEN
        SystemScaling = Scaling
      ELSE
        SystemScaling = .FALSE.
      END IF

      z(1:n/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      !-------------------------------------------------------------------------
      ! If the consistent splitting approach is used, 
      ! compute the continuous approximation of the continuity equation residual
      !-------------------------------------------------------------------------
      IF ( ConsistentSplitting ) THEN  
        y(1:q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
        IF (SystemScaling) THEN
          !-------------------------------
          ! Recover the unscaled residual
          !--------------------------------
          DO j=1,q/2 
            f(j) = 1.0d0/da( (dim+2)*j ) * V( (dim+2)*j )
          END DO
        ELSE
          DO j=1,q/2
            f(j) = V( (dim+2)*j )
          END DO
        END IF
        !-----------------------------------------------------------------
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        CALL Info( 'AcousticsSolver', &
            'Iteration for the continuity equation residual', Level=4)
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        !------------------------------------------------------------------
        CALL ComplexBiCGStab( q, MMatrix, y, f, InnerRounds, 1.0d-7, 0 )
        CALL ComplexMatrixVectorProduct( LMatrix, y, f )
      END IF

      !-----------------------------------------------------------------
      ! The solution of the Schur complement equation...
      !------------------------------------------------------------------ 
      w(1:q) = CMPLX( 0.0d0,0.0d0, kind=dp )
      !------------------------------------------------------------------
      ! Construct the right-hand side vector g for the Schur complement
      ! system. The vector g consists of the unscaled residuals and the
      ! modification term which arises from the consistent splitting.
      !------------------------------------------------------------------ 
      IF ( ConsistentSplitting ) THEN 

        IF ( SystemScaling ) THEN
          DO j=1,q/2
            g(2*j-1) = 1.0d0 / da( (dim+2)*j-1 ) * V((dim+2)*j-1)
            g(2*j) = 1.0d0 / da( (dim+2)*j ) * V((dim+2)*j) + f(j)
          END DO
        ELSE
          DO j=1,q/2
            g(2*j-1) = V((dim+2)*j-1)
            g(2*j) = V((dim+2)*j) + f(j)
          END DO
        END IF

      ELSE

        IF ( SystemScaling ) THEN
          g(2*j-1) = 1.0d0 / da( (dim+2)*j-1 ) * V((dim+2)*j-1)
          g(2*j) = 1.0d0 / da( (dim+2)*j ) * V((dim+2)*j)         
        ELSE
          DO j=1,q/2
            g(2*j-1) = V((dim+2)*j-1)
            g(2*j) = V((dim+2)*j)
          END DO
        END IF

      END IF

      !--------------------------------------------------------------------
      ! If scaling is used, we must scale g according to the scaling of 
      ! the Schur complement matrix 
      !--------------------------------------------------------------------
      IF ( SystemScaling ) THEN
        DO j=1,q/2
          g(2*j-1) = dps( 2*j-1 ) * g(2*j-1)
          g(2*j) = dps( 2*j ) * g(2*j)         
        END DO
      END IF
      !-------------------------------------------------------------------
      ! Solve the Schur complement system...
      !-------------------------------------------------------------------
      IF ( ANY( g /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
        !-----------------------------------------------------------------
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        CALL Info( 'AcousticsSolver', &
            'Preconditioning iteration for the Schur complement system', Level=4)
        CALL Info( 'AcousticsSolver', ' ', Level=4)
        !-------------------------------------------------------------------
        IF (SchurCriterion) THEN 
          CALL ComplexBiCGStabl( 4, q, PS, w, g, InnerRounds, SchurTol, 0 )
        ELSE
          CALL ComplexBiCGStabl( 4, q, PS, w, g, InnerRounds, InnerTol, 0 )
          !CALL ComplexBiCGStabl( 2, q, PS, w, g, 5, 1.0d-6, 0 )
        END IF
      END IF

      DO j=1,q/2
        V((dim+2)*j-1) = w(2*j-1) 
        V((dim+2)*j) = w(2*j)
      END DO

      !---------------------------------------------------------------------------------
      ! The computation of the specific matrix vector product related to preconditining
      !---------------------------------------------------------------------------------
      CALL ComplexMatrixVectorProduct2( A, V, z, dim )

      !--------------------------------------------------------
      ! The solution of the velocity preconditioning equation
      !--------------------------------------------------------
      Vel(1:dim*q/2) = CMPLX( 0.0d0,0.0d0, kind=dp )
      DO i=1,dim
        DO j=1,q/2
          VelRhs((j-1)*dim+i) = z((dim+2)*(j-1)+i)
        END DO
      END DO


      IF ( ASSOCIATED(PA) ) THEN
        IF ( SystemScaling ) THEN
          DO i=1,dim
            DO j=1,q/2
              VelRhs((j-1)*dim+i) = dpa( (j-1)*dim+i ) / da( (dim+2)*(j-1)+i ) * &
                  VelRhs((j-1)*dim+i)
            END DO
          END DO
        END IF

        IF ( ANY( VelRhs /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Preconditioning iteration for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          IF (VelocityCriterion) THEN
            CALL ComplexBiCGStabl( 2, q*dim/2, PA, Vel, VelRhs, InnerRounds, VelocityTol, 0)            
          ELSE
            CALL ComplexBiCGStabl( 2, q*dim/2, PA, Vel, VelRhs, InnerRounds, InnerTol, 0)
          END IF
        END IF

      ELSE

        IF ( ANY( VelRhs /= CMPLX( 0.0d0,0.0d0, KIND=dp ) ) ) THEN
          !-----------------------------------------------------------------
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          CALL Info( 'AcousticsSolver', &
              'Preconditioning iteration for velocities', Level=4)
          CALL Info( 'AcousticsSolver', ' ', Level=4)
          !------------------------------------------------------------------
          IF (VelocityCriterion) THEN
            CALL VelocitySolve( 2, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, VelocityTol, 0)            
          ELSE
            CALL VelocitySolve( 2, q*dim/2, A, Vel, VelRhs, dim, InnerRounds, InnerTol, 0)
          END IF
        END IF
      END IF

      DO j=1,q/2
        DO i=1,dim
          z((dim+2)*(j-1)+i) = Vel((j-1)*dim+i)
        END DO
      END DO

      V(1:n/2)=z(1:n/2)
!------------------------------------------------------------------------------
    END SUBROUTINE PreconditioningIteration
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  FUNCTION StoppingCriterion( n, A, x, b, res, CriterionType ) RESULT(err)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n), res(n)
    INTEGER, OPTIONAL :: CriterionType
    REAL(kind=dp) :: err
!-----------------------------------------------------------------------------
    REAL(kind=dp) :: norm, tmp, normb, normres, normx
    INTEGER :: i, j
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
!--------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    normb = 0.0d0
    normres = 0.0d0
    normx = 0.0d0
    
    IF (PRESENT(CriterionType)) THEN
      DO i=1,n
        normb = MAX( normb, CDABS(b(i)) )
        normres = MAX( normres, CDABS(res(i)) )
      END DO
      err = normres/normb
    ELSE
      DO i=1,n
        tmp = 0.0d0
        DO j=Rows(2*i-1),Rows(2*i)-1,2
          tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
        END DO
        IF (tmp > norm) norm = tmp
        normb = MAX( normb, CDABS(b(i)) )
        normres = MAX( normres, CDABS(res(i)) )
        normx = MAX( normx, CDABS(x(i)) )
      END DO
      err = normres / (norm * normx + normb)
    END IF

!------------------------------------------------------------------------------
  END FUNCTION StoppingCriterion
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ConditionEstimate( n, A, x, b ) RESULT(err)
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(kind=dp) :: err
!------------------------------------------------------------------------------
    INTEGER :: i, j
    REAL(kind=dp) :: norm, tmp, normb, normx
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
!-------------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    norm = 0.0d0
    normb = 0.0d0
    normx = 0.0d0
    DO i=1,n
      tmp = 0.0d0
      DO j=Rows(2*i-1),Rows(2*i)-1,2
        tmp = tmp + CDABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
      END DO
      IF (tmp > norm) norm = tmp
      normb = MAX( normb, CDABS(b(i)) )
      normx = MAX( normx, CDABS(x(i)) )
    END DO

    err = norm*normx/normb 
!------------------------------------------------------------------------------
  END FUNCTION ConditionEstimate
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVectorProduct( A,u,v )
!------------------------------------------------------------------------------
!   Matrix vector product (v = Au) for a matrix given in CRS format.
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: u,v
    TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)

    INTEGER :: i,j,n
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = A % NumberOfRows / 2
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    v(1:n) = CMPLX( 0.0d0, 0.0d0, kind=dp )
    DO i=1,n
       DO j=Rows(2*i-1),Rows(2*i)-1,2
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) + s * u((Cols(j)+1)/2)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ComplexMatrixVectorProduct
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ComplexInfinityNorm( m, x ) RESULT(s)
!------------------------------------------------------------------------------
    INTEGER :: m
    COMPLEX(KIND=dp) :: x(m)
    REAL(KIND=dp) :: s  
!--------------------------------------------------------------------
    s = 0.0d0
    DO i=1,m
      s = MAX( s, ABS(x(i)) )
    END DO
!--------------------------------------
 END FUNCTION ComplexInfinityNorm
!--------------------------------------




!------------------------------------------------------------------------------
    FUNCTION ComplexNorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n, i
       REAL(KIND=dp) :: s
       COMPLEX(KIND=dp) :: r, x(:)
!------------------------------------------------------------------------------
       s =  SQRT( REAL( DOT_PRODUCT( x(1:n), x(1:n) ) ) )
!------------------------------------------------------------------------------
    END FUNCTION ComplexNorm
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
    FUNCTION ComplexDotProduct( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       COMPLEX(KIND=dp) :: s, x(:), y(:)
!------------------------------------------------------------------------------
       s = DOT_PRODUCT( x(1:n), y(1:n) )
!------------------------------------------------------------------------------
    END FUNCTION ComplexDotProduct
!------------------------------------------------------------------------------





!-----------------------------------------------------------------------------------
  SUBROUTINE ComplexBiCGStabl( l, n, A, x, b, MaxRounds, Tol, StoppingCriterionType )
!-----------------------------------------------------------------------------------
!   This subroutine solves complex linear systems by using the BiCGStab(l) algorithm 
!   with l >= 2. It has been developed by using as a starting point the work of D.R. Fokkema 
!   (subroutine zbistbl v1.1 1998). Dr. Fokkema has given the right to distribute
!   the derived work under GPL and hence the original more conservative 
!   copyright notice of the subroutine has been removed accordingly. 
!
!   This version uses right-oriented ILU(n) preconditioning.
!----------------------------------------------------------------------------------- 
    INTEGER :: l   ! polynomial degree
    INTEGER :: n, MaxRounds   
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(KIND=dp) :: Tol
    INTEGER, OPTIONAL :: StoppingCriterionType 
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: zzero, zone, t(n), kappa0, kappal 
    REAL(KIND=dp) :: dznrm2, rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
        delta = 1.0d-2, bnrm, bw_errorind, tottime
    INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
        IluOrder
    COMPLEX(KIND=dp) :: work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
        alpha, beta, omega, rho0, rho1, sigma, zdotc, varrho, hatgamma
    LOGICAL rcmp, xpdt, GotIt, BackwardError  
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    tottime = CPUTime()
    BackwardError = .TRUE.
    IF ( PRESENT(StoppingCriterionType) ) THEN
      IF (StoppingCriterionType == 0) BackwardError = .FALSE.
    END IF
    IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b

    zzero = dcmplx( 0.0d0,0.0d0)
    zone =  dcmplx( 1.0d0,0.0d0)
    work = dcmplx( 0.0d0, 0.0d0 )
    rwork = dcmplx( 0.0d0, 0.0d0 )
    
    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1

    z = 1
    zz = z+(l+1)
    y0 = zz+(l+1)
    yl = y0+1
    y = yl+1

    CALL ComplexMatrixVectorProduct( A, x, work(1:n,r) )
    work(1:n,r) = b(1:n) - work(1:n,r)
    bnrm = dznrm2(n, b(1:n), 1)

    work(1:n,rr) = work(1:n,r) 
    work(1:n,bp) = work(1:n,r)
    work(1:n,xp) = x(1:n)
    x(1:n) = zzero
    rnrm0 = dznrm2(n, work(1:n,r), 1)
    rnrm = rnrm0
    mxnrmx = rnrm0
    mxnrmr = rnrm0
    
    alpha = zzero
    omega = zone
    sigma = zone
    rho0 = zone

    Round = 0
    errorind = 1.0d0
    DO WHILE ( errorind > Tol .AND. Round < MaxRounds) 
      Round = Round + 1 
      !-------------------------
      ! --- The BiCG part ---
      !-------------------------
      rho0 = -omega*rho0

      DO k=1,l
        rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
        IF (rho0 == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        beta = alpha*(rho1/rho0)
        rho0 = rho1
        DO j=0,k-1
          work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
        ENDDO
        t(1:n) = work(1:n,u+k-1)
        CALL CRS_ComplexLUSolve2( n, A, t )
        CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) )      
        sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
        IF (sigma == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        alpha = rho1/sigma
        x(1:n) = x(1:n) + alpha * work(1:n,u)
        DO j=0,k-1
          work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
        ENDDO
        t(1:n) = work(1:n,r+k-1)
        CALL CRS_ComplexLUSolve2( n, A, t ) 
        CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )  
        rnrm = dznrm2(n, work(1:n,r), 1)
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
      ENDDO

      !--------------------------------------
      ! --- The convex polynomial part ---
      !--------------------------------------

      DO i=1,l+1
        DO j=1,i
          rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 ) 
        END DO
      END DO
      DO j=2,l+1
        rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
      END DO

      rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
      CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
          iwork, stat)

      ! --- tilde r0 and tilde rl (small vectors)

      rwork(1,y0) = -zone
      rwork(2:l,y0) = rwork(2:l,z) 
      CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,y0), l-1, stat)
      rwork(l+1,y0) = zzero

      rwork(1,yl) = zzero
      rwork(2:l,yl) = rwork(2:l,z+l) 
      CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,yl), l-1, stat)
      rwork(l+1,yl) = -zone

      ! --- Convex combination

      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
      kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
          (kappa0*kappal)
      hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
          kappa0/kappal
      rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

      !  --- Update

      omega = rwork(l+1,y0)
      DO j=1,l
        work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
        x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
        work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
      ENDDO


      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------

      IF (.TRUE.) THEN
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          t(1:n) = x(1:n)
          CALL CRS_ComplexLUSolve2( n, A, t )         
          CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
          work(1:n,r) = work(1:n,bp) - work(1:n,r)
          mxnrmr = rnrm
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            work(1:n,xp) = work(1:n,xp) + t(1:n)
            x(1:n) = zzero
            work(1:n,bp) = work(1:n,r)
            mxnrmx = rnrm
          ENDIF
        ENDIF
      ENDIF

      IF (rcmp) THEN
        IF (xpdt) THEN       
          t(1:n) = work(1:n,xp)
        ELSE
          t(1:n) = t(1:n) + work(1:n,xp)  
        END IF
      ELSE
        t(1:n) = x(1:n)
        CALL CRS_ComplexLUSolve2( n, A, t ) 
        t(1:n) =  t(1:n) + work(1:n,xp)
      END IF

      bw_errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )
      errorind = rnrm/bnrm
      WRITE(*,'(I4,ES12.3,ES12.3)') Round, errorind, bw_errorind
      
      IF (BackwardError) errorind = bw_errorind

    END DO

    !------------------------------------------------------------
    ! We have solved z = P*x, so finally solve the true unknown x
    !------------------------------------------------------------
    CALL CRS_ComplexLUSolve2( n, A, x )
    x(1:n) = x(1:n) + work(1:n,xp)

    WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
        ConditionEstimate( n, A, x, b )
    WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
        ComplexNorm( n, x )  

!------------------------------------------------------------------------------
  END SUBROUTINE ComplexBiCGStabl
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
    SUBROUTINE ComplexBiCGStab( n, A, x, b, Rounds, TOL, StoppingCriterionType)
!------------------------------------------------------------------------------
!   This is the ILU or Jacobi preconditioned BiCGStab method for the complex 
!   linear system Ax = b. The optional argument StoppingCriterionType can be 
!   used to define the type of the stopping criterion.
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: n, Rounds
      REAL(KIND=dp) :: TOL
      COMPLEX(kind=dp) :: x(n/2), b(n/2)
      INTEGER, OPTIONAL :: StoppingCriterionType 
!------------------------------------------------------------------------------
      INTEGER :: i, m, k
      LOGICAL :: BackwardError 
      REAL(KIND=dp) :: res, tottime, res0, const
      COMPLEX(KIND=dp) :: r(n/2),Ri(n/2),P(n/2),V(n/2),T(n/2),T1(n/2),T2(n/2),&
          S(n/2)
      COMPLEX(KIND=dp) :: alpha,beta,omega,rho,oldrho
!------------------------------------------------------------------------------

      BackwardError = .TRUE.
      IF ( PRESENT(StoppingCriterionType) ) THEN
        IF (StoppingCriterionType == 0) BackwardError = .FALSE.
      END IF

      IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b
      m = n/2

      CALL ComplexMatrixVectorProduct( A, x, r )
      r(1:m) = b(1:m) - r(1:m)
      IF ( .NOT. BackwardError) res0 = ComplexNorm(m,b)

      Ri(1:m) = r(1:m)
      P(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      V(1:m) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      omega  = CMPLX( 1.0d0, 0.0d0, kind=dp)
      alpha  = CMPLX( 0.0d0, 0.0d0, kind=dp)
      oldrho = CMPLX( 1.0d0, 0.0d0, kind=dp)
      tottime = CPUTime()

      DO i=1,Rounds
        rho = ComplexDotProduct( m, Ri, r )
        beta = alpha * rho / ( oldrho * omega )
        P(1:m) = r(1:m) + beta * (P(1:m) - omega*V(1:m))
        V(1:m) = P(1:m)

        CALL CRS_ComplexLUSolve2( m, A, V )

        T1(1:m) = V(1:m)
        CALL ComplexMatrixVectorProduct( A, T1, V )
        alpha = rho / ComplexDotProduct( m, Ri, V )
        S(1:m) = r(1:m) - alpha * V(1:m)
        x(1:m) = x(1:m) + alpha*T1(1:m)

        IF (BackwardError) THEN
          res = StoppingCriterion( m, A, x, b, S)
        ELSE
          res = ComplexNorm(m,S)/res0
        END IF

        IF ( res < TOL ) THEN
          WRITE(*,'(I4,ES12.3)') i, res
          WRITE(*,'(a,F8.2)') 'Solution time (s):    ', CPUTime() - tottime
          WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
              ConditionEstimate( m, A, x, b )
          WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
              ComplexNorm( m, x )          
          EXIT
        END IF

        T(1:m) = S(1:m)

        CALL CRS_ComplexLUSolve2( m, A, T )
           
        T2(1:m) = T(1:m)
        CALL ComplexMatrixVectorProduct( A, T2, T )
        omega = ComplexDotProduct( m,T,S ) / ComplexDotProduct( m,T,T )
        oldrho = rho
        r(1:m) = S(1:m) - omega*T(1:m)
        x(1:m) = x(1:m) + omega*T2(1:m)

        IF (BackwardError) THEN
          res = StoppingCriterion( m, A, x, b, r )
        ELSE
          res = ComplexNorm(m,r)/res0
        END IF

        WRITE(*,'(I4,ES12.3)') i, res
        IF ( res < TOL ) THEN
          WRITE(*,'(a,F8.2)') 'Solution time (s):    ', CPUTime() - tottime
          WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
              ConditionEstimate( m, A, x, b )
          WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
              ComplexNorm( m, x )              
          EXIT
        END IF
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE ComplexBiCGStab
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexLUSolve2( N,A,b )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  DESCRIPTION:
!    Solve a (complex9 system (Ax=b) after factorization A=LUD has been
!    done. This routine is meant as a part of  a preconditioner for an
!    iterative solver.
!
!  ARGUMENTS:
!
!    INTEGER :: N
!      INPUT: Size of the system
!
!    TYPE(Matrix_t) :: A
!      INPUT: Structure holding input matrix
!
!    DOUBLE PRECISION :: b
!      INOUT: on entry the RHS vector, on exit the solution vector.
!
!******************************************************************************
!------------------------------------------------------------------------------
 
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: N
    COMPLEX(KIND=dp) :: b(N)

!------------------------------------------------------------------------------

    COMPLEX(KIND=dp), POINTER :: Values(:)
    INTEGER :: i,j
    COMPLEX(KIND=dp) :: x, s
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    
!------------------------------------------------------------------------------

    Diag => A % ILUDiag
    Rows => A % ILURows
    Cols => A % ILUCols
    Values => A % CILUValues

!
!   if no ilu provided do diagonal solve:
!   -------------------------------------
    IF ( .NOT. ASSOCIATED( Values ) ) THEN
       Diag => A % Diag

       !DO i=1,n/2    ! THIS IS PROBABLY ERRATIC IN CVS-VERSION
       DO i=1,n
          x = CMPLX( A % Values(Diag(2*i-1)), -A % Values(Diag(2*i-1)+1), kind=dp )
          b(i) = b(i) / x
       END DO
       RETURN
    END IF

    IF (.FALSE.) THEN
      CALL ComplexLUSolve( n,SIZE(Cols),Rows,Cols,Diag,Values,b )
    ELSE
      ! Forward substitute
      DO i=1,n
        s = b(i)
        DO j=Rows(i),Diag(i)-1
          s = s - Values(j) * b(Cols(j))
        END DO
        b(i) = s
      END DO
      !
      ! Backward substitute
      DO i=n,1,-1
        s = b(i)
        DO j=Diag(i)+1,Rows(i+1)-1
          s = s - Values(j) * b(Cols(j))
        END DO
        b(i) = Values(Diag(i)) * s
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_ComplexLUSolve2
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE ComplexLUSolve( n,m,Rows,Cols,Diag,Values,b )
!-----------------------------------------------------------------------------
    INTEGER :: n,m,Rows(n+1),Cols(m),Diag(n)
    COMPLEX(KIND=dp) :: Values(m),b(n)
    INTEGER :: i,j


    ! Forward substitute
    DO i=1,n
      DO j=Rows(i),Diag(i)-1
        b(i) = b(i) - Values(j) * b(Cols(j))
      END DO
    END DO

    ! Backward substitute
    DO i=n,1,-1
      DO j=Diag(i)+1,Rows(i+1)-1
        b(i) = b(i) - Values(j) * b(Cols(j))
      END DO
      b(i) = Values(Diag(i)) * b(i)
    END DO
!-----------------------------------------------------------------------------
  END SUBROUTINE ComplexLUSolve
!------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVectorProduct2( A,u,v,dim )
!------------------------------------------------------------------------------
!
!   The computation of a specific matrix-vector product for preconditioning  
!
!------------------------------------------------------------------------------

    COMPLEX(KIND=dp), DIMENSION(*) :: u,v
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)

    INTEGER :: i,j,k,n,p,q
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = A % NumberOfRows / 2
    q = n/(dim+2)

    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag 
    Values => A % Values

    v(1:n) = u(1:n) 

    DO k=1,q
      DO p=1,dim
        i = (k-1)*(dim+2)+p
        DO j = Rows(2*i-1)+2*dim, Rows(2*i)-1, 2*(dim+2)       
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) - s * u((Cols(j)+1)/2)
        END DO
        DO j = Rows(2*i-1)+2*(dim+1), Rows(2*i)-1, 2*(dim+2)       
          s = CMPLX( Values(j), -Values(j+1), kind=dp )
          v(i) = v(i) - s * u((Cols(j)+1)/2)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------------
   END SUBROUTINE ComplexMatrixVectorProduct2
!------------------------------------------------------------------------------






!-----------------------------------------------------------------------------------
  SUBROUTINE VelocitySolve( l, n, A, x, b, dim, MaxRounds, Tol, StoppingCriterionType )
!-----------------------------------------------------------------------------------
!   This subroutine solves the velocity preconditioning system without requiring 
!   that the assemby has been made for the preconditioning system. The required
!   matrix-vector products are performed by extracting the required entries from
!   the primary coefficient matrix A. This version uses diagonal preconditioning.
!----------------------------------------------------------------------------------- 
    INTEGER :: l   ! polynomial degree
    INTEGER :: n, MaxRounds, dim  
    TYPE(Matrix_t), POINTER :: A
    COMPLEX(KIND=dp) :: x(n), b(n)
    REAL(KIND=dp) :: Tol
    INTEGER, OPTIONAL :: StoppingCriterionType 
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: zzero, zone, t(n), kappa0, kappal 
    REAL(KIND=dp) :: dznrm2, rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
        delta = 1.0d-2, bnrm, bw_errorind, tottime
    INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
        IluOrder
    COMPLEX(KIND=dp) :: work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)), &
        alpha, beta, omega, rho0, rho1, sigma, zdotc, varrho, hatgamma
    LOGICAL rcmp, xpdt, GotIt, BackwardError  
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    tottime = CPUTime()
    BackwardError = .FALSE.
    IF ( PRESENT(StoppingCriterionType) ) THEN
      IF (StoppingCriterionType == 0) BackwardError = .FALSE.
    END IF
    IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b

    zzero = dcmplx( 0.0d0,0.0d0)
    zone =  dcmplx( 1.0d0,0.0d0)
    work = dcmplx( 0.0d0, 0.0d0 )
    rwork = dcmplx( 0.0d0, 0.0d0 )
    
    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1

    z = 1
    zz = z+(l+1)
    y0 = zz+(l+1)
    yl = y0+1
    y = yl+1

    CALL ComplexMatrixVelocityVectorProduct( A, x, work(1:n,r), dim )
    work(1:n,r) = b(1:n) - work(1:n,r)
    bnrm = dznrm2(n, b(1:n), 1)

    work(1:n,rr) = work(1:n,r) 
    work(1:n,bp) = work(1:n,r)
    work(1:n,xp) = x(1:n)
    x(1:n) = zzero
    rnrm0 = dznrm2(n, work(1:n,r), 1)
    rnrm = rnrm0
    mxnrmx = rnrm0
    mxnrmr = rnrm0
    
    alpha = zzero
    omega = zone
    sigma = zone
    rho0 = zone

    Round = 0
    errorind = 1.0d0
    DO WHILE ( errorind > Tol .AND. Round < MaxRounds) 
      Round = Round + 1 
      !-------------------------
      ! --- The BiCG part ---
      !-------------------------
      rho0 = -omega*rho0

      DO k=1,l
        rho1 = zdotc(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
        IF (rho0 == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        beta = alpha*(rho1/rho0)
        rho0 = rho1
        DO j=0,k-1
          work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
        ENDDO
        t(1:n) = work(1:n,u+k-1)
        !CALL CRS_ComplexLUSolve2( n, A, t )
        CALL DiagonalVelocityPreconditioning( A, t, dim)
        !CALL ComplexMatrixVectorProduct( A, t, work(1:n,u+k) ) 
        CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,u+k), dim ) 
        sigma = zdotc(n, work(1:n,rr), 1, work(1:n,u+k), 1)
        IF (sigma == zzero) THEN
          CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
        ENDIF
        alpha = rho1/sigma
        x(1:n) = x(1:n) + alpha * work(1:n,u)
        DO j=0,k-1
          work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
        ENDDO
        t(1:n) = work(1:n,r+k-1)
        !CALL CRS_ComplexLUSolve2( n, A, t )
        CALL DiagonalVelocityPreconditioning( A, t, dim) 
        !CALL ComplexMatrixVectorProduct( A, t, work(1:n,r+k) )
        CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,r+k), dim )  
        rnrm = dznrm2(n, work(1:n,r), 1)
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
      ENDDO

      !--------------------------------------
      ! --- The convex polynomial part ---
      !--------------------------------------

      DO i=1,l+1
        DO j=1,i
          rwork(i,j) = zdotc(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 ) 
        END DO
      END DO
      DO j=2,l+1
        rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
      END DO

      rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
      CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
          iwork, stat)

      ! --- tilde r0 and tilde rl (small vectors)

      rwork(1,y0) = -zone
      rwork(2:l,y0) = rwork(2:l,z) 
      CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,y0), l-1, stat)
      rwork(l+1,y0) = zzero

      rwork(1,yl) = zzero
      rwork(2:l,yl) = rwork(2:l,z+l) 
      CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
          rwork(2:l,yl), l-1, stat)
      rwork(l+1,yl) = -zone

      ! --- Convex combination

      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
      kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )
      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
          (kappa0*kappal)
      hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
          kappa0/kappal
      rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

      !  --- Update

      omega = rwork(l+1,y0)
      DO j=1,l
        work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
        x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
        work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
      ENDDO


      CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
      rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

      !---------------------------------------
      !  --- The reliable update part ---
      !---------------------------------------

      IF (.TRUE.) THEN
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          t(1:n) = x(1:n)
          !CALL CRS_ComplexLUSolve2( n, A, t )
          CALL DiagonalVelocityPreconditioning( A, t, dim)         
          !CALL ComplexMatrixVectorProduct( A, t, work(1:n,r) )
          CALL ComplexMatrixVelocityVectorProduct( A, t, work(1:n,r), dim )
          work(1:n,r) = work(1:n,bp) - work(1:n,r)
          mxnrmr = rnrm
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            work(1:n,xp) = work(1:n,xp) + t(1:n)
            x(1:n) = zzero
            work(1:n,bp) = work(1:n,r)
            mxnrmx = rnrm
          ENDIF
        ENDIF
      ENDIF

      IF (rcmp) THEN
        IF (xpdt) THEN       
          t(1:n) = work(1:n,xp)
        ELSE
          t(1:n) = t(1:n) + work(1:n,xp)  
        END IF
      ELSE
        t(1:n) = x(1:n)
        !CALL CRS_ComplexLUSolve2( n, A, t ) 
        CALL DiagonalVelocityPreconditioning( A, t, dim)
        t(1:n) =  t(1:n) + work(1:n,xp)
      END IF

      ! bw_errorind = StoppingCriterion( n, A, t, b, work(1:n,r) )
      errorind = rnrm/bnrm
      WRITE(*,'(I4,ES12.3,ES12.3)') Round, errorind
      
      !IF (BackwardError) errorind = bw_errorind

    END DO

    !------------------------------------------------------------
    ! We have solved z = P*x, so finally solve the true unknown x
    !------------------------------------------------------------
    !CALL CRS_ComplexLUSolve2( n, A, x )
    CALL DiagonalVelocityPreconditioning( A, x, dim)
    x(1:n) = x(1:n) + work(1:n,xp)

    !WRITE(*,'(a,ES12.3)') 'An approximate lower bound for the condition number ', &
    !    ConditionEstimate( n, A, x, b )
    WRITE(*,'(a,ES12.3)') 'The 2-norm of the solution: ', &
        ComplexNorm( n, x )  

!------------------------------------------------------------------------------
  END SUBROUTINE VelocitySolve
!------------------------------------------------------------------------------
















!-------------------------------------------------------------------------------
  SUBROUTINE ComplexMatrixVelocityVectorProduct( K,u,v,dim )
!------------------------------------------------------------------------------
!
!   The computation of a specific matrix-vector product needed in preconditioning.
!   The subroutine computes the matrix-vector product v = Au where A is
!   the coefficient matrix for unknown velocities occupying
!   the (1,1) block of the coeffiecient matrix K. 
!
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: u, v
    TYPE(Matrix_t), POINTER :: K
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: i, j, l, m, n, p, q, t
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = K % NumberOfRows / 2
    q = n/(dim+2)

    Rows   => K % Rows
    Cols   => K % Cols
    Values => K % Values

    v(1:dim*q) = CMPLX( 0.0d0, 0.0d0, kind=dp) 

    DO m=1,q
      DO p=1,dim
        i = (m-1)*(dim+2)+p
        DO l = 1,dim
          DO j = Rows(2*i-1)+2*(l-1), Rows(2*i)-1, 2*(dim+2)       
            s = CMPLX( Values(j), -Values(j+1), kind=dp )
            t = (Cols(j)+1)/2
            v((m-1)*dim+p) = v((m-1)*dim+p) + s * u( (t-l)/(dim+2)*dim+p )  
          END DO
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------------
  END SUBROUTINE ComplexMatrixVelocityVectorProduct
!------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
  SUBROUTINE DiagonalVelocityPreconditioning( K, v, dim )
!------------------------------------------------------------------------------
!
!   This subroutine performs a specific diagonal preconditioning needed in
!   solving the velocities. The subroutine computes the matrix-vector product 
!   v = inv(P)* v where P is the diagonal matrix consisting of the diagonal 
!   entries of the coefficient matrix for unknown velocities occupying
!   the (1,1) block of the coefficient matrix K. 
!
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*) :: v
    TYPE(Matrix_t), POINTER :: K
    INTEGER :: dim
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: i, j, n, p, q
    COMPLEX(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = K % NumberOfRows / 2
    q = n/(dim+2)

    Diag   => K % Diag 
    Values => K % Values

    DO j=1,q
      DO p=1,dim
        i = (j-1)*(dim+2)+p
        s = CMPLX( Values(Diag(2*i-1)), -Values(Diag(2*i-1)+1), kind=dp )
        v((j-1)*dim+p) = v((j-1)*dim+p) / s  
      END DO
    END DO

!-----------------------------------------------------------------------------
  END SUBROUTINE DiagonalVelocityPreconditioning
!------------------------------------------------------------------------------







!------------------------------------------------------------------------------
   SUBROUTINE CoupledVelocityMatrix(  StiffMatrix, Viscosity, AngularFrequency, &
       Density, Element, n, dim)
!------------------------------------------------------------------------------
!    This subroutine computes the element matrix for the preconditioning
!    equation of velocities. The coordinate system can be either orthogonal
!    Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
     REAL(KIND=dp) :: Viscosity(:), AngularFrequency, Density(:)
     INTEGER :: dim, n
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
     COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
         mu, rho0, s
     LOGICAL :: Stat
     INTEGER :: t, i, j, p, q
     TYPE(GaussIntegrationPoints_t) :: IP
     TYPE(Nodes_t) :: Nodes
     SAVE Nodes
!------------------------------------------------------------------------------
     CALL GetElementNodes( Nodes )
     StiffMatrix = 0.0d0
     CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
     !-------------------------
     ! Numerical integration:
     !-------------------------
     IP = GaussPoints( Element,n)
     DO t=1,IP % n
       !--------------------------------------------------------------
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       s = IP % s(t) * detJ
       IF (CoordSys == AxisSymmetric) THEN
         r = SUM( Basis * Nodes % x(1:n) )
         s = r * s
       END IF

       !-----------------------------------------------
       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Viscosity(1:n) )
       rho0 = SUM( Density(1:n) * Basis(1:n) )
       !---------------------------------------------
       ! the stiffness matrix...
       !---------------------------------------------
       DO i=1,dim
         DO p=1,n
           DO q=1,n
             DO j = 1,dim
               CStiff((p-1)*dim+i,(q-1)*dim+i) = CStiff((p-1)*dim+i,(q-1)*dim+i) + &
                   s * CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), KIND=dp ) * &
                   dBasisdx(q,j) * dBasisdx(p,j)
               CStiff((p-1)*dim+i, (q-1)*dim+j) = &
                   CStiff((p-1)*dim+i, (q-1)*dim+j) + &
                   CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), KIND=dp ) * &
                   dBasisdx(q,i) * dBasisdx(p,j) * s
             END DO

             IF ( (i==1) .AND. (CoordSys == AxisSymmetric) ) THEN
               CStiff((p-1)*dim+i, (q-1)*dim+i) = &
                   CStiff((p-1)*dim+i, (q-1)*dim+i) + &
                   CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), KIND=dp ) * 1/r**2 * &
                   Basis(q) * Basis(p) * s  
             END IF

             CStiff((p-1)*dim+i,(q-1)*dim+i) = CStiff((p-1)*dim+i,(q-1)*dim+i) + &
                 s * CMPLX(1.0d0, 0.0d0, KIND=dp ) * Basis(q) * Basis(p)             
           END DO
         END DO
       END DO
     END DO


     DO p=1,n
       DO i=1,DIM
         DO q=1,n
           DO j=1,DIM
             StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                 REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                 -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                 AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
             StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                 REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
           END DO
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CoupledVelocityMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE PressureLaplaceMatrix(  StiffMatrix, Viscosity, AngularFrequency, &
       Density, Element, n, dim)
!------------------------------------------------------------------------------
!   This subroutine computes the element matrix for the Laplacian type 
!   term arising in the consistent splitting approach. The coordinate system can 
!   be either orthogonal Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: Viscosity(:), AngularFrequency, Density(:)
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(n,n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        mu, rho0, s
    LOGICAL :: Stat
    INTEGER :: t, i, j, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
      
      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      mu  = SUM( Basis(1:n) * Viscosity(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO q=1,n
          DO j = 1,dim
            CStiff(p,q) = CStiff(p,q) - &
                s * CMPLX( 0.0d0, 2.0d0*mu/(AngularFrequency*rho0), KIND=dp ) * &
                dBasisdx(q,j) * dBasisdx(p,j)
          END DO
        END DO
      END DO
    END DO

    DO i=1,n
      DO j=1,n
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( CStiff(i,j) )
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE PressureLaplaceMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE PressureMassMatrix(  StiffMatrix, Element, n, dim)
!------------------------------------------------------------------------------
!   This subroutine computes the element matrix for the mass matrix which is 
!   needed in the consistent splitting approach to compute the continuous 
!   approximation of residual. The coordinate system can be either orthogonal 
!   Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(n,n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, s
    LOGICAL :: Stat
    INTEGER :: t, i, j, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, KIND=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO q=1,n
          CStiff(p,q) = CStiff(p,q) + s * CMPLX(1.0d0, 0.0d0, KIND=dp) * &
              Basis(q) * Basis(p)             
        END DO
      END DO
    END DO

    DO i=1,n
      DO j=1,n
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+1 ) =  REAL( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+1, 2*(j-1)+2 ) = -AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+1 ) =  AIMAG( CStiff(i,j) )
        StiffMatrix( 2*(i-1)+2, 2*(j-1)+2 ) =  REAL( CStiff(i,j) )
      END DO
    END DO
!------------------------------------------------------------------------------
    END SUBROUTINE PressureMassMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE SchurComplementMatrix(  StiffMatrix, AngularFrequency , &
        SpecificHeat, HeatRatio, Density, Pressure,               &
        Temperature, Conductivity, Viscosity, Lambda,             &
        Element, n, dim)
!------------------------------------------------------------------------------
!     This subroutine computes the element matrix for the preconditioning
!     equation of temperature and pressure. The coordinate system can be either 
!     orthogonal Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
      REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
      REAL(KIND=dp) :: AngularFrequency, SpecificHeat(:), HeatRatio(:), Density(:), &
          Pressure(:),  Temperature(:), Conductivity(:), Viscosity(:), Lambda(:) 
      TYPE(Element_t), POINTER :: Element
      INTEGER :: n, dim
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
          CV, kappa, mu, rho0, P0, gamma, la, K1, K2, s
      COMPLEX(kind=dp) :: CStiff(2*n,2*n), C1, C2, C3, C4
      LOGICAL :: Stat
      INTEGER :: t, i, j, p, q
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t) :: Nodes
      SAVE Nodes
!------------------------------------------------------------------------------
      CALL GetElementNodes( Nodes )
      StiffMatrix = 0.0d0
      CStiff = CMPLX( 0.0d0,0.0d0, kind=dp )
      !-------------------------
      ! Numerical integration:
      !-------------------------
      IP = GaussPoints( Element,n)
      DO t=1,IP % n
        !--------------------------------------------------------------
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

        s = IP % s(t) * detJ
        IF (CoordSys == AxisSymmetric) THEN
          r = SUM( Basis * Nodes % x(1:n) )
          s = r * s
        END IF

        !-----------------------------------------------
        ! Material parameters at the integration point:
        !----------------------------------------------
        CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
        kappa = SUM( Conductivity(1:n) * Basis(1:n) )
        mu  = SUM( Basis(1:n) * Viscosity(1:n) )
        rho0  = SUM( Basis(1:n) * Density(1:n) )
        P0  = SUM( Basis(1:n) * Pressure(1:n) )
        gamma  = SUM( Basis(1:n) * HeatRatio(1:n) )
        la = SUM( Basis(1:n) * Lambda(1:n) )

        C1 = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
            CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  

        C2 = CMPLX( rho0*AngularFrequency,0.0d0, kind=dp) / &
            CMPLX( P0/AngularFrequency, la, kind=dp )
        C3 = CMPLX( 1.0d0, 0.0d0, kind=dp)
        C4 = CMPLX( 1.0d0,0.0d0, kind=dp) / &
            CMPLX( 1.0d0, AngularFrequency/P0*(la), kind=dp )

        K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)
        K2 = 1.0d0/(gamma-1.0d0)

        !---------------------------------------------
        ! the stiffness matrix...
        !----------------------------------------
        DO p=1,n
          DO q=1,n
            DO i=1,dim

              CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
                  s * CMPLX( 0.0d0, K1, kind=dp) * dBasisdx(q,i) * dBasisdx(p,i)

              CStiff((p-1)*2+2,(q-1)*2+1) = CStiff((p-1)*2+2,(q-1)*2+1) + &
                  s * C3 * dBasisdx(q,i) * dBasisdx(p,i)

              CStiff((p-1)*2+2,(q-1)*2+2) = CStiff((p-1)*2+2,(q-1)*2+2) + &
                  s * C1 * dBasisdx(q,i) * dBasisdx(p,i)
            END DO
            CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
                s * CMPLX(-K2, 0.0d0, kind=dp) * Basis(q) * Basis(p)

            CStiff((p-1)*2+1,(q-1)*2+2) = CStiff((p-1)*2+1,(q-1)*2+2) + &
                s * C4 * Basis(q) * Basis(p)

            CStiff((p-1)*2+2,(q-1)*2+2) = CStiff((p-1)*2+2,(q-1)*2+2) - &
                s * C2 * Basis(q) * Basis(p)
          END DO
        END DO
      END DO

      DO p=1,n
        DO i=1,2
          DO q=1,n
            DO j=1,2
              StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                  REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                  -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                  AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
              StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                  REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            END DO
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE SchurComplementMatrix
!------------------------------------------------------------------------------






!------------------------------------------------------------------------------
  SUBROUTINE VelocityImpedanceMatrix(  StiffMatrix, AngularFrequency, &
      Density, Impedance, Element, n, Nodes, dim)
!------------------------------------------------------------------------------
!   This subroutine computes the contribution of the impedance boundary
!   condition to the preconditioning equation of velocities. The coordinate 
!   system can be either orthogonal Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: AngularFrequency, Density(:), Impedance(:,:)
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: rho0, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
      !-----------------------------------------------
      ! Material parameters etc. at the integration point:
      !----------------------------------------------
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      Impedance1 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(1,1:n) * Basis(1:n) )
      Impedance2 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(2,1:n) * Basis(1:n) ) 
      !---------------------------------------------
      ! the stiffness matrix...
      !---------------------------------------------
      DO p=1,n
        DO i=1,dim
          DO q=1,n
            DO j=1,dim
              CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                  CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &  
                  CMPLX(-Impedance2, Impedance1, kind=dp) * &
                  Basis(q) * Normal(j) * Basis(p) * Normal(i) * s
            END DO
          END DO
        END DO
      END DO
    END DO   ! Loop over integration points

    DO p=1,n
      DO i=1,DIM
        DO q=1,n
          DO j=1,DIM
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE VelocityImpedanceMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SchurComplementImpedanceMatrix( StiffMatrix, Impedance, &
      AngularFrequency, SpecificHeat, HeatRatio, Density, Conductivity, &
      Element, n, Nodes, dim)
!------------------------------------------------------------------------------
!   This subroutine computes the contribution of the impedance boundary
!   conditions to the preconditioning equation of temperature and pressure. 
!   The coordinate system can be either orthogonal Cartesian or axially symmetric.       
!-------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: Impedance(:,:), AngularFrequency, SpecificHeat(:), &
        HeatRatio(:), Density(:), Conductivity(:)
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(2*n,2*n), K1, ZT
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: CV, kappa, rho0, gamma, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      rho0  = SUM( Basis(1:n) * Density(1:n) )
      gamma  = SUM( Basis(1:n) * HeatRatio(1:n) )

      K1 = CMPLX( 0.0d0, kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV), kind=dp )
      ZT = CMPLX( SUM( Impedance(3,1:n) * Basis(1:n) ), SUM( Impedance(4,1:n) * Basis(1:n) ), kind=dp )
      
      DO p=1,n
        DO q=1,n
          CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) - &
              s * K1 * ZT * Basis(q) * Basis(p)
        END DO
      END DO
    END DO
   
    DO p=1,n
      DO i=1,2
        DO q=1,n
          DO j=1,2
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SchurComplementImpedanceMatrix
!------------------------------------------------------------------------------






!------------------------------------------------------------------------------
  SUBROUTINE VelocitySlipMatrix(  StiffMatrix, SpecificHeat, &
      HeatRatio, Density, Temperature, &
      AngularFrequency, WallTemperature, SlipCoefficient1, &
      Element, n, Nodes, dim)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: SpecificHeat(:), HeatRatio(:), Density(:), &
        Temperature(:), AngularFrequency, WallTemperature(:), SlipCoefficient1
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(dim*n,dim*n)    
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        CV, gamma, rho0, T0, WallT0, C1, s, &
        Normal(3), Tangent1(3), Tangent2(3)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0,kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )
      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
      !-----------------------------------------------
      ! Material parameters etc. at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )

      C1 = -SlipCoefficient1/(2.0d0-SlipCoefficient1) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      Normal = Normalvector(Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      CALL TangentDirections(Normal, Tangent1, Tangent2)

      DO p=1,n
        DO i=1,dim
          DO j=1,dim
            DO q=1,n
              CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                  CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  Basis(q) * Tangent1(j) * Basis(p) * Tangent1(i) * s
            END DO
          END DO
        END DO
      END DO

      IF (dim > 2) THEN
        DO p=1,n
          DO i=1,dim
            DO j=1,dim
              DO q=1,n
                CStiff( (p-1)*DIM+i, (q-1)*DIM+j) = &
                    CStiff( (p-1)*DIM+i, (q-1)*DIM+j) + &
                    CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    Basis(q) * Tangent2(j) * Basis(p) * Tangent2(i) * s
              END DO
            END DO
          END DO
        END DO
      END IF
    END DO

    DO p=1,n
      DO i=1,DIM
        DO q=1,n
          DO j=1,DIM
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j-1 ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i-1, 2*DIM*(q-1)+2*j ) =  &
                -AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j-1 ) =  &
                AIMAG( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
            StiffMatrix( 2*DIM*(p-1)+2*i, 2*DIM*(q-1)+2*j ) =  &
                REAL( CSTIFF(DIM*(p-1)+i,DIM*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE VelocitySlipMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SchurComplementSlipMatrix( StiffMatrix, SpecificHeat, &
      HeatRatio, Density, Temperature, AngularFrequency, Conductivity, & 
      WallTemperature, SlipCoefficient2, &
      Element, n, Nodes, dim)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: StiffMatrix(:,:)
    REAL(KIND=dp) :: SpecificHeat(:), HeatRatio(:), Density(:), Temperature(:), &
        AngularFrequency, Conductivity(:), WallTemperature(:), SlipCoefficient2
    INTEGER :: dim, n
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    COMPLEX(kind=dp) :: CStiff(2*n,2*n)
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), DetJ, r, &
        Normal(3), Impedance1, Impedance2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: CV, kappa, rho0, gamma, T0, WallT0, K1, C2, s
!------------------------------------------------------------------------------
    StiffMatrix = 0.0d0
    CStiff = CMPLX(0.0d0,0.0d0, kind=dp)
    !-------------------------
    ! Numerical integration:
    !-------------------------
    IP = GaussPoints( Element,n)
    DO t=1,IP % n
      !--------------------------------------------------------------
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = IP % s(t) * detJ
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF

      !-----------------------------------------------
      ! Material parameters at the integration point:
      !----------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )

      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)
      C2 = 1/kappa*SlipCoefficient2/(2.0d0-SlipCoefficient2) * &
          (gamma+1.0d0)/2.0d0 * rho0 * CV * &
          SQRT(2.0d0*(gamma-1.0d0)*CV*(T0+WallT0)/PI)
      
      DO p=1,n
        DO q=1,n
          CStiff((p-1)*2+1,(q-1)*2+1) = CStiff((p-1)*2+1,(q-1)*2+1) + &
              s * CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * Basis(q) * Basis(p)
        END DO
      END DO
    END DO
   
    DO p=1,n
      DO i=1,2
        DO q=1,n
          DO j=1,2
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j-1 ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i-1, 4*(q-1)+2*j ) =  &
                -AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j-1 ) =  &
                AIMAG( CStiff(2*(p-1)+i,2*(q-1)+j) )
            StiffMatrix( 4*(p-1)+2*i, 4*(q-1)+2*j ) =  &
                REAL( CStiff(2*(p-1)+i,2*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SchurComplementSlipMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  StiffMatrix, Force, AngularFrequency , &
      SpecificHeat, HeatRatio, Density,               &
      Temperature, Conductivity, Viscosity, Lambda,             &
      HeatSource, Load, Bubbles, Mini_Bubbles, Element, n, Nodes, Dofs, nb )
!------------------------------------------------------------------------------
!    This subroutine computes the element stiffness matrix. The coordinate 
!    system can be either orthogonal Cartesian or axially symmetric.       
!    The bubble basis can be specified by the user.
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, &
        SpecificHeat(:), HeatRatio(:), Density(:),    &       
        Temperature(:), Conductivity(:), Viscosity(:), Lambda(:),  &
        HeatSource(:,:), Load(:,:)
    LOGICAL :: Bubbles, Mini_Bubbles
    INTEGER :: n, Dofs, nb
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*n), dBasisdx(2*n,3), ddBasisddx(n,3,3)
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, L(6), &
        CV, gamma, rho0, P0, T0, kappa, mu, la, f1, f2, K1, K2, K3, r  
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2)), A, &
        SchurConst, C1

    INTEGER :: i, j, p, q, t, DIM, NBasis, VelocityDofs, & 
        VelocityComponents
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp) :: X, Y, Z, Metric(3,3), SqrtMetric, Symb(3,3,3), &
        dSymb(3,3,3,3), StabTerms(n), AK, k, BubbleLaplacian
    LOGICAL :: AxialSymmetry, MyBubbles = .FALSE.
!------------------------------------------------------------------------------
    AxialSymmetry = (CoordSys == AxisSymmetric)

    StabTerms = 0.0d0
    AK = ElementArea(Solver % Mesh, Element, n)

    DIM = CoordinateSystemDimension()

    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
    !------------------------------------------------------------------------------
    !   Numerical integration
    !------------------------------------------------------------------------------
    IF ( Bubbles .OR. Mini_Bubbles) THEN
      IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
      IF (Mini_Bubbles) THEN
        NBasis = n+1 
      ELSE
        NBasis = 2*n
      ENDIF
    ELSE       
      NBasis = n + nb
      IntegStuff = GaussPoints( Element, RelOrder = RelOrder )
    END IF

    !------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
       !------------------------------------------------------------------------------
       !  Basis function values & derivatives at the integration point
       !------------------------------------------------------------------------------
       IF (Mini_Bubbles .OR. MyBubbles) THEN
         stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
             Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
        
         dBasisdx(n+1:,:) = 0._dp
         SELECT CASE( Element % TYPE % ElementCode ) 
           
         CASE(504)
           Basis(n+1) = 1.0d0
           DO i=1,n
             Basis(n+1) = Basis(n+1)*Basis(i)
           END DO
       
           DO j=1,dim
             dBasisdx(n+1,j) = dBasisdx(1,j)*Basis(2)*Basis(3)*Basis(4) + &
                 dBasisdx(2,j)*Basis(1)*Basis(3)*Basis(4) + &
                 dBasisdx(3,j)*Basis(1)*Basis(2)*Basis(4) + &
                 dBasisdx(4,j)*Basis(1)*Basis(2)*Basis(3)  
           END DO
          
         CASE(404)
           Basis(n+1) = Basis(1)*Basis(3)           
           DO j=1,dim
             dBasisdx(n+1,j) = dBasisdx(1,j)*Basis(3) + &
                 dBasisdx(3,j)*Basis(1)
           END DO
          
         CASE(303)
           Basis(n+1) = Basis(1)*Basis(2)*Basis(3)           
           DO j=1,dim
             dBasisdx(n+1,j) = dBasisdx(1,j)*Basis(2)*Basis(3) + &
                 dBasisdx(2,j)*Basis(1)* Basis(3) + &
                 dBasisdx(3,j)*Basis(1)*Basis(2)
           END DO

         CASE(706)
           Basis(n+1) = ( Basis(1)+Basis(2)+Basis(3) ) * ( Basis(4)+Basis(5)+Basis(6) ) * &
               ( Basis(1)+Basis(4) ) * ( Basis(2)+Basis(5) ) * ( Basis(3)+Basis(6) )
           DO j=1,dim
             dBasisdx(n+1,j) = ( dBasisdx(1,j)+dBasisdx(2,j)+dBasisdx(3,j) ) * &
                 ( Basis(4)+Basis(5)+Basis(6) ) * &
                 ( Basis(1)+Basis(4) ) * ( Basis(2)+Basis(5) ) * ( Basis(3)+Basis(6) ) + &
                 ( dBasisdx(4,j)+dBasisdx(5,j)+dBasisdx(6,j) ) * &
                 ( Basis(1)+Basis(2)+Basis(3) ) * &
                 ( Basis(1)+Basis(4) ) * ( Basis(2)+Basis(5) ) * ( Basis(3)+Basis(6) ) + &  
                 ( Basis(1)+Basis(2)+Basis(3) ) * ( Basis(4)+Basis(5)+Basis(6) ) * &
                 ( dBasisdx(1,j)+dBasisdx(4,j) ) * ( Basis(2)+Basis(5) ) * ( Basis(3)+Basis(6) ) + &
                 ( Basis(1)+Basis(2)+Basis(3) ) * ( Basis(4)+Basis(5)+Basis(6) ) * &
                 ( dBasisdx(2,j)+dBasisdx(5,j) ) * ( Basis(1)+Basis(4) ) * ( Basis(3)+Basis(6) ) + &
                 ( Basis(1)+Basis(2)+Basis(3) ) * ( Basis(4)+Basis(5)+Basis(6) ) * &
                 ( dBasisdx(3,j)+dBasisdx(6,j) ) * ( Basis(1)+Basis(4) ) * ( Basis(2)+Basis(5) ) 
           END DO

         CASE(808)
           Basis(n+1) = Basis(1)*Basis(7)
           DO j=1,dim
             dBasisdx(n+1,j) = dBasisdx(1,j)*Basis(7) + &
                 dBasisdx(7,j)*Basis(1)
           END DO
          
         CASE DEFAULT
           WRITE( Message, '(a,i4,a)' ) 'Mini-bubbles for element ', &
               Element % TYPE % ElementCode, ' are not implemented.'
           CALL Error( 'Acoustics', Message )
           
         END SELECT
        
         IF (MyBubbles) THEN
            ! Second derivatives of bubble function... 
            SELECT CASE( Element % TYPE % ElementCode )
            
            CASE(303)
               BubbleLaplacian = 0.0d0
               DO j = 1,dim
                  BubbleLaplacian = BubbleLaplacian + &
                       2.0d0 * dBasisdx(1,j) * dBasisdx(2,j) * Basis(3) + &
                       2.0d0 * dBasisdx(1,j) * dBasisdx(3,j) * Basis(2) + & 
                       2.0d0 * dBasisdx(2,j) * dBasisdx(3,j) * Basis(1)                       
               END DO
                
            CASE DEFAULT
               WRITE( Message, '(a,i4,a)' ) 'Second derivatives for element ', &
                    Element % TYPE % ElementCode, ' are not implemented.'
               CALL Error( 'Acoustics', Message ) 
               
            END SELECT
         END IF

       ELSE

          IF ( nb > 0) THEN
             stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
                  Basis, dBasisdx )
          ELSE
             stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
                  Basis, dBasisdx, ddBasisddx, .FALSE., Bubbles )
          END IF

       ENDIF

       s = s * SqrtElementMetric 
       IF (AxialSymmetry) THEN
         r = SUM( Basis(1:n) * Nodes % x(1:n) )
         s = r * s
       END IF

       !------------------------------------------------------------------------------
       !  Problem parameters and the real and imaginary part of the 
       !  load at the integration point
       !------------------------------------------------------------------------------
       CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
       gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
       rho0 = SUM( Density(1:n) * Basis(1:n) )
       !       P0 = SUM( Pressure(1:n) * Basis(1:n) )
       T0 = SUM( Temperature(1:n) * Basis(1:n) )
       kappa = SUM( Conductivity(1:n) * Basis(1:n) )
       mu = SUM( Viscosity(1:n) * Basis(1:n) )
       la = SUM( Lambda(1:n) * Basis(1:n) )

       P0 = (gamma-1.0d0)*CV*rho0*T0
       !       K1 = kappa*AngularFrequency/(CV*(gamma-1.0d0)*P0)
       K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)
       !       K2 = AngularFrequency**2*rho0/(P0*(gamma-1.0d0))
       K2 = 1.0d0/(gamma-1.0d0)
       K3 = 1.0d0/AngularFrequency**2

       SchurConst = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
          CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  
       k = 1.0d0/SQRT( gamma*P0/(rho0*AngularFrequency**2) )   

       C1 = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
            CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  

       f1 = K3*SUM( HeatSource(1,1:n) * Basis(1:n) )
       f2 = K3*SUM( HeatSource(2,1:n) * Basis(1:n) )
       DO i=1,2*DIM
         L(i) = 1.0d0/AngularFrequency*SUM( Load(i,1:n) * Basis(1:n) )
       END DO

       !------------------------------------------------------------------------------
       !  The stiffness matrix and load vector...
       !  The following is the contribution from the heat equation and pressure
       !  equation, i.e. the part arising from the loop over the test functions 
       !  for temperature and pressure 
       !------------------------------------------------------------------------------
       DO p=1,N
         StabTerms(p) = StabTerms(p) + s * Basis(p)
         DO q=1,N
           A = CMPLX( -K2, 0.0d0, kind=dp ) * Basis(q) * Basis(p)
           DO i=1,DIM
             !-----------------------------------------------
             !  Coefficients for the nodal velocities...
             !-----------------------------------------------
             LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+i) =  &
                 LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+i) + &
                 CMPLX(0.0d0, 1.0d0, kind=dp) * dBasisdx(q,i) * Basis(p) * s 

             IF ((i==1) .AND. AxialSymmetry) THEN
               LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+i) =  &
                   LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+i) + &
                   CMPLX(0.0d0, 1.0d0, kind=dp) * 1/r * Basis(q) * Basis(p) * s 
             END IF

             A = A + CMPLX( 0.0d0, K1, kind=dp ) * &
                   dBasisdx(q,i) * dBasisdx(p,i)
           END DO
           !------------------------------------------------------------------------------
           !  Coefficients for the nodal temperatures
           !------------------------------------------------------------------------------
           LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1 ) = &
               LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1 ) + s*A
           !------------------------------------------------------------------------------
           !  Coefficients for the nodal pressures
           !------------------------------------------------------------------------------
           LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+2) =  &
               LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+2) + &
                s * 1.0d0/CMPLX( 1.0d0, AngularFrequency/P0*la, kind=dp ) * &
                Basis(q) * Basis(p)

           LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
                LSTIFF( p*(DIM+2), q*(DIM+2) ) - &
                s * rho0*AngularFrequency/CMPLX( P0/AngularFrequency, la, kind=dp ) * &
                Basis(q) * Basis(p)

           !--------------------------------------------------------------------------------
           !  Try stabilisation without bubbles...
           !-------------------------------------------------------------------------------
           IF (.FALSE.) THEN  !(.NOT. Bubbles) THEN
             LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
                 LSTIFF( p*(DIM+2), q*(DIM+2) ) - CMPLX( 1.0d0/k**2, 0.0d0, kind=dp) * & 
                 CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
                 CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  * &
                 s * Basis(q) * Basis(p) 
            
             LSTIFF( p*(DIM+2), q*(DIM+2)-1 ) = &
                 LSTIFF( p*(DIM+2), q*(DIM+2)-1 ) - CMPLX( 1.0d0/k**2, 0.0d0, kind=dp) * & 
                 s * Basis(q) * Basis(p)
           END IF

         END DO
         LFORCE( (p-1)*(DIM+2)+DIM+1 ) = &
             LFORCE( (p-1)*(DIM+2)+DIM+1 ) + s * Basis(p) * CMPLX( -f2,f1, kind=dp )
       END DO

       IF ( Bubbles .OR. Mini_Bubbles .OR. (nb > 0) ) THEN
         DO p=1,n
           DO q=n+1,NBasis
             DO i=1,DIM
               !------------------------------------------------------------------------------
               !  Coefficients for the nodal velocities
               !------------------------------------------------------------------------------
               LSTIFF( p*(DIM+2), (q-1)*DIM+2*n+i) = &
                   LSTIFF( p*(DIM+2), (q-1)*DIM+2*n+i) + &
                   CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(q,i) * Basis(p) * s

               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF( p*(DIM+2), (q-1)*DIM+2*n+i) = &
                     LSTIFF( p*(DIM+2), (q-1)*DIM+2*n+i) + &
                     CMPLX( 0.0d0, 1.0d0, kind=dp ) * 1/r * Basis(q) * Basis(p) * s  
               END IF
             END DO
           END DO
         END DO
       END IF

       !------------------------------------------------------------------------------
       !  The following is the contribution from the NS-equation, i.e. the 
       !  part arising from the loop over the test functions for velocity 
       !------------------------------------------------------------------------------
       DO i=1,DIM
         DO p=1,n
           DO q=1,n
             !------------------------------------------------------------------------------
             !  Coefficients for the nodal temperatures...               
             !------------------------------------------------------------------------------
             LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+DIM+1 ) = &
                 LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+DIM+1 ) + &
                 CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(p,i) * Basis(q) * s
             IF ((i==1) .AND. AxialSymmetry) THEN
               LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+DIM+1 ) = &
                   LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+DIM+1 ) + &
                   CMPLX( 0.0d0, 1.0d0, kind=dp ) * 1/r * Basis(p) * Basis(q) * s
             END IF
             !------------------------------------------------------------------------------
             !  Coefficients for the nodal pressures...               
             !------------------------------------------------------------------------------
             LSTIFF( (p-1)*(DIM+2)+i, q*(DIM+2) ) = &
                 LSTIFF( (p-1)*(DIM+2)+i, q*(DIM+2) ) + &
                 CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(p,i) * Basis(q) * s
             IF ((i==1) .AND. AxialSymmetry) THEN
               LSTIFF( (p-1)*(DIM+2)+i, q*(DIM+2) ) = &
                   LSTIFF( (p-1)*(DIM+2)+i, q*(DIM+2) ) + &
                   CMPLX( 0.0d0, 1.0d0, kind=dp ) * 1/r * Basis(p) * Basis(q) * s  
             END IF
             !------------------------------------------------------------------------------
             !  Coefficients for the nodal velocities...             
             !------------------------------------------------------------------------------
             LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+i) = &
                 LSTIFF( (p-1)*(DIM+2)+i,(q-1)*(DIM+2)+i) + &
                 CMPLX( 1.0d0, 0.0d0, kind=dp ) * &
                 Basis(q) * Basis(p) * s

             !------------------------------------------------------------------------------
             !  grad(v)grav(w)-type terms
             !------------------------------------------------------------------------------                
             DO j=1,DIM
               LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+i) = &
                   LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+i) + &
                   CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * &
                   dBasisdx(q,j) * dBasisdx(p,j) * s
               LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                   LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                   CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * &
                   dBasisdx(q,i) * dBasisdx(p,j) * s
             END DO

             IF ((i==1) .AND. AxialSymmetry) THEN
               LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+i) = &
                   LSTIFF((p-1)*(DIM+2)+i, (q-1)*(DIM+2)+i) + &
                   CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), kind=dp ) * 1/r**2 * Basis(q) * Basis(p) * s  
             END IF
           END DO
           !------------------------------------------------------------------------------
           !  The components of force vector...
           !------------------------------------------------------------------------------
           LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i ) - &
               s * Basis(p) * CMPLX( -L((i-1)*2+2), L((i-1)*2+1), kind=dp  )
         END DO
       END DO

       IF ( Bubbles .OR. Mini_Bubbles .OR. ( nb>0 ) ) THEN
         DO i=1,DIM
           DO p=1,n
             DO q=n+1,NBasis
               !------------------------------------------------------------------------------
               !  coefficients for the nodal velocities             
               !------------------------------------------------------------------------------
               LSTIFF( (p-1)*(DIM+2)+i, (q-1)*DIM+2*n+i) = &
                   LSTIFF( (p-1)*(DIM+2)+i,(q-1)*DIM+2*n+i) + &
                   CMPLX( 1.0d0, 0.0d0, kind=dp ) * &
                   Basis(q) * Basis(p) * s

               !------------------------------------------------------------------------------
               !  grad(v)grav(w)-type terms
               !------------------------------------------------------------------------------                
               DO j=1,DIM
                 LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+i) = &
                     LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+i) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * &
                     dBasisdx(q,j) * dBasisdx(p,j) * s
                 LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+j) = &
                     LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+j) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * &
                     dBasisdx(q,i) * dBasisdx(p,j) * s
               END DO

               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+i) = &
                     LSTIFF((p-1)*(DIM+2)+i, (q-1)*DIM+2*n+i) + &
                     CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), kind=dp ) * 1/r**2 * Basis(q) * Basis(p) * s 
               END IF
             END DO
           END DO
         END DO

         DO i=1,DIM
           DO p=n+1,NBasis
             DO q=1,n
               !------------------------------------------------------------------------------
               !  Coefficients for the nodal temperatures               
               !------------------------------------------------------------------------------
               LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+DIM+1 ) = &
                   LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+DIM+1 ) + &
                   CMPLX( 0.0d0, 1.0d0, kind=dp) * dBasisdx(p,i) * Basis(q) * s
               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+DIM+1 ) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+DIM+1 ) + &
                     CMPLX( 0.0d0, 1.0d0, kind=dp) * 1/r * Basis(p) * Basis(q) * s
               END IF
               !------------------------------------------------------------------------------
               !  Coefficients for the nodal pressures...               
               !------------------------------------------------------------------------------
               LSTIFF( (p-1)*DIM+2*n+i, q*(DIM+2) ) = &
                   LSTIFF( (p-1)*DIM+2*n+i, q*(DIM+2) ) + &
                   CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(p,i) * Basis(q) * s
               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF( (p-1)*DIM+2*n+i, q*(DIM+2) ) = &
                     LSTIFF( (p-1)*DIM+2*n+i, q*(DIM+2) ) + &
                     CMPLX( 0.0d0, 1.0d0, kind=dp ) * 1/r * Basis(p) * Basis(q) * s  
               END IF
               !------------------------------------------------------------------------------
               !  coefficients for the nodal velocities 
               !------------------------------------------------------------------------------
               LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+i) = &
                   LSTIFF( (p-1)*DIM+2*n+i,(q-1)*(DIM+2)+i) + &
                   CMPLX( 1.0d0, 0.0d0, kind=dp ) * &
                   Basis(q) * Basis(p) * s

               !------------------------------------------------------------------------------
               !  grad(v)grav(w)-type terms
               !------------------------------------------------------------------------------                
               DO j=1,DIM
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+i) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+i) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp) * dBasisdx(q,j) * dBasisdx(p,j) * s
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+j) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+j) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp) * dBasisdx(q,i) * dBasisdx(p,j) * s
               END DO

               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+i) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*(DIM+2)+i) + &
                     CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), kind=dp) * 1/r**2 * Basis(q) * Basis(p) * s
               END IF
             END DO
             !------------------------------------------------------------------------------
             !  The components of force vector...
             !------------------------------------------------------------------------------
             LFORCE( (p-1)*DIM+2*n+i ) = LFORCE( (p-1)*DIM+2*n+i ) - &
                 s * Basis(p) * CMPLX( -L((i-1)*2+2), L((i-1)*2+1), kind=dp )
           END DO
         END DO

         DO i=1,DIM
           DO p=n+1,NBasis
             DO q=n+1,NBasis
               !------------------------------------------------------------------------------
               !  coefficients for the nodal velocities 
               !------------------------------------------------------------------------------
               LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+i) = &
                   LSTIFF( (p-1)*DIM+2*n+i,(q-1)*DIM+2*n+i) + &
                   CMPLX( 1.0d0, 0.0d0, kind=dp ) * &
                   Basis(q) * Basis(p) * s

               !------------------------------------------------------------------------------
               !  grad(v)grav(w)-type terms
               !------------------------------------------------------------------------------                
               DO j=1,DIM
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+i) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+i) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * dBasisdx(q,j) * dBasisdx(p,j) * s
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+j) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+j) + &
                     CMPLX( 0.0d0, -mu/(AngularFrequency*rho0), kind=dp ) * dBasisdx(q,i) * dBasisdx(p,j) * s
               END DO
               IF ((i==1) .AND. AxialSymmetry) THEN
                 LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+i) = &
                     LSTIFF( (p-1)*DIM+2*n+i, (q-1)*DIM+2*n+i) + &
                     CMPLX( 0.0d0, -2*mu/(AngularFrequency*rho0), kind=dp) * 1/r**2 * Basis(q) * Basis(p) * s
               END IF
            END DO
           END DO
         END DO         
       END IF   ! IF (Bubbles)...


       IF (MyBubbles) THEN

          DO j=1,dim
             LSTIFF( (dim+2)*n+1, (dim+2)*n+1 ) = LSTIFF( (dim+2)*n+1, (dim+2)*n+1 ) + &
                  dBasisdx(n+1,j) * dBasisdx(n+1,j) * s                
          END DO

          DO q=1,n
             DO j=1,dim
                LSTIFF( (dim+2)*n+1, dim*n+1 ) = LSTIFF( (dim+2)*n+1, dim*n+1) + &
                     CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(n+1,j) * dBasisdx(q,j) * s

                LSTIFF( (dim+2)*n+1, dim*n+2 ) = LSTIFF( (dim+2)*n+1, dim*n+2) + &
                     C1 * CMPLX( 0.0d0, 1.0d0, kind=dp ) * dBasisdx(n+1,j) * dBasisdx(q,j) * s
                
                LSTIFF( (dim+2)*n+1, (q-1)*(dim+2)+j ) = LSTIFF( (dim+2)*n+1, (q-1)*(dim+2)+j ) - &
                  Basis(q) * dBasisdx(n+1,j)

             END DO
          END DO

       END IF



    !------------------------------------------------------------------------------
    END DO      ! Loop over integration points
    !------------------------------------------------------------------------------
    IF (.FALSE.) THEN    !(.NOT. Bubbles) THEN
      DO p=1,n
        DO q=1,n

          LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
              LSTIFF( p*(DIM+2), q*(DIM+2) ) + CMPLX( 1.0d0/k**2, 0.0d0, kind=dp ) * &
              CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
              CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  * &
              1.0d0/AK * StabTerms(p) * StabTerms(q)
          
          LSTIFF( p*(DIM+2), q*(DIM+2)-1 ) = &
              LSTIFF( p*(DIM+2), q*(DIM+2)-1 ) + CMPLX( 1.0d0/k**2, 0.0d0, kind=dp ) * &
              1.0d0/AK * StabTerms(p) * StabTerms(q)
          
        END DO
      END DO
    END IF

    IF ( nb > 0 ) THEN
       CALL LCondensateBubbles( n, dim, LSTIFF, LFORCE, nb )
    ELSE
       IF ( Bubbles .OR. Mini_Bubbles) THEN
          IF (Mini_Bubbles) THEN
             CALL LCondensateBubbles( n, dim, LSTIFF, LFORCE, 1 )  
          ELSE
             CALL LCondensate( n, dim, LSTIFF, LFORCE ) 
          END IF
       END IF
    END IF

!    if (MyBubbles) then
!        CALL LCondensateMyBubble( n, dim, LSTIFF, LFORCE )  
!    end if



    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LCondensateBubbles( n, dim, K, F, nb )
!------------------------------------------------------------------------------
    USE LinearAlgebra
!------------------------------------------------------------------------------
    INTEGER :: n, dim, nb
    COMPLEX(KIND=dp) :: K(:,:), F(:), Kbb(dim*nb,dim*nb), &
         Kbl(dim*nb,(dim+2)*n), Klb((dim+2)*n,dim*nb), Fb(nb*dim)

    INTEGER :: i, Ldofs((dim+2)*n), Bdofs(dim*nb)

    Ldofs = (/ (i, i=1,(dim+2)*n) /)
    Bdofs = (/ ((dim+2)*n+i, i=1, dim*nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL ComplexInvertMatrix( Kbb,nb*dim )

    F(1:(dim+2)*n) = F(1:(dim+2)*n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:(dim+2)*n,1:(dim+2)*n) = &
         K(1:(dim+2)*n,1:(dim+2)*n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensateBubbles
!------------------------------------------------------------------------------









!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( n, dim, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
!------------------------------------------------------------------------------
    INTEGER :: n, dim
    COMPLEX(KIND=dp) :: K(:,:), F(:), Kbb(dim*n,dim*n), &
         Kbl(dim*n,(dim+2)*n), Klb((dim+2)*n,dim*n), Fb(n*dim)

    INTEGER :: i, Ldofs((dim+2)*n), Bdofs(dim*n)

    Ldofs = (/ (i, i=1,(dim+2)*n) /)
    Bdofs = (/ ((dim+2)*n+i, i=1, dim*n) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL ComplexInvertMatrix( Kbb,n*dim )

    F(1:(dim+2)*n) = F(1:(dim+2)*n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:(dim+2)*n,1:(dim+2)*n) = &
         K(1:(dim+2)*n,1:(dim+2)*n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary(  StiffMatrix, Force, AngularFrequency, &
      SpecificHeat, HeatRatio, Density, Pressure,               &
      Temperature, Conductivity, Impedance, Load, &
      Element, n, Nodes, Dofs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, & 
        SpecificHeat(:), HeatRatio(:), Density(:), Pressure(:),    &
        Temperature(:), Conductivity(:), Impedance(:,:), &
        Load(:,:)
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, Impedance1, Impedance2, &
        Impedance3, Impedance4, CV, gamma, rho0, P0, T0, kappa, K1, L(6),   & 
        Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), X, Y, Z, Normal(3), r
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2))
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, DIM
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
      s = s * SqrtElementMetric
      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
!------------------------------------------------------------------------------
!     Problem parameters at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      P0 = SUM( Pressure(1:n) * Basis(1:n) )
      T0 = SUM( Temperature(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
!      K1 = kappa*AngularFrequency/(CV*(gamma-1.0d0)*P0)
      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

      Impedance1 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(1,1:n) * Basis(1:n) )
      Impedance2 = 1.0d0/(AngularFrequency*rho0) * SUM( Impedance(2,1:n) * Basis(1:n) ) 
      Impedance3 = SUM( Impedance(3,1:n) * Basis(1:n) ) 
      Impedance4 = SUM( Impedance(4,1:n) * Basis(1:n) ) 

      DO i=1,2*DIM
        L(i) = 1.0d0/(AngularFrequency*rho0) * SUM( Load(i,1:n) * Basis(1:n) )
      END DO

      DO p=1,n
        DO i=1,dim
          DO q=1,n
            DO j=1,dim
              LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                  LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                  CMPLX(-Impedance2, Impedance1, kind=dp) * &
                  Basis(q) * Normal(j) * Basis(p) * Normal(i) * s
            END DO
          END DO
          LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i ) - &
              s * Basis(p) * CMPLX( -L((i-1)*2+2), L((i-1)*2+1), kind=dp  )
        END DO
      END DO

      DO p=1,n
        DO q=1,n
          LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) = &
              LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) - &
              CMPLX(0.0d0, K1, kind=dp) * CMPLX(Impedance3, Impedance4, kind=dp) * &
              Basis(q) * Basis(p) * s
        END DO
      END DO

    END DO   ! Loop over integration points

    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------









!------------------------------------------------------------------------------
  SUBROUTINE LocalInterfaceMatrix(  StiffMatrix, Force, AngularFrequency, &
      SpecificHeat, HeatRatio, Density, Pressure,               &
      Temperature, Conductivity, Impedance, Load, &
      Element, n, Nodes, Dofs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, & 
        SpecificHeat(:), HeatRatio(:), Density(:), Pressure(:),    &
        Temperature(:), Conductivity(:), Impedance(:,:), &
        Load(:,:)
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, Impedance1, Impedance2, &
        Impedance3, Impedance4, CV, gamma, rho0, P0, T0, kappa, K1, L(6),   & 
        Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), X, Y, Z, Normal(3), r
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2))
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, DIM
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()
    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
      s = s * SqrtElementMetric

      IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
      END IF
!------------------------------------------------------------------------------
!     Problem parameters at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      P0 = SUM( Pressure(1:n) * Basis(1:n) )
      T0 = SUM( Temperature(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
!      K1 = kappa*AngularFrequency/(CV*(gamma-1.0d0)*P0)
      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

      DO j=1,dim
        L(2*j-1) = -1.0d0/(AngularFrequency*rho0) * SUM( Load(1,1:n) * Basis(1:n) ) * Normal(j)
        L(2*j) = -1.0d0/(AngularFrequency*rho0) * SUM( Load(2,1:n) * Basis(1:n) ) * Normal(j)
      END DO

      DO p=1,n
        DO i=1,dim
          LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i ) - &
              s * Basis(p) * CMPLX( -L((i-1)*2+2), L((i-1)*2+1), kind=dp  )
        END DO
      END DO

    END DO   ! Loop over integration points

    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalInterfaceMatrix
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE SlipMatrix(  StiffMatrix, Force, SpecificHeat, HeatRatio, &
      Density, Conductivity, Pressure, Temperature, &
      AngularFrequency, WallTemperature, &
      WallVelocity, SlipCoefficient1, SlipCoefficient2, SlipCoefficient3, &
      Element, n, Nodes, Dofs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), & 
        SpecificHeat(:), HeatRatio(:), Density(:), Conductivity(:), &
        Pressure(:), Temperature(:), AngularFrequency, WallTemperature(:), &
        WallVelocity(:,:), SlipCoefficient1, SlipCoefficient2, SlipCoefficient3
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, Impedance1, Impedance2, &
        Impedance3, Impedance4, CV, gamma, rho0, P0, T0, ReV0(3), ImV0(3), &
        WallT0, kappa, K1, L(6),C1, C2, C3, & 
        Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), X, Y, Z, Normal(3), r,  &
        Tangent1(3), Tangent2(3)
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2))
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, DIM
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()
    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1, IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
      IF (CoordSys == Cartesian) THEN
        s = s * SqrtElementMetric
      ELSE
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s * SqrtElementMetric
      END IF
!------------------------------------------------------------------------------
!     Problem parameters at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      P0 =  SUM( Pressure(1:n) * Basis(1:n) )
      T0 =  SUM( Temperature(1:n) * Basis(1:n) )
      WallT0 = SUM( WallTemperature(1:n) * Basis(1:n) )
      
      DO i=1,3
        ReV0(i) = SUM( WallVelocity((i-1)*2+1,1:n) * Basis(1:n) )
        ImV0(i) = SUM( WallVelocity((i-1)*2+2,1:n) * Basis(1:n) )
      END DO

      K1 = kappa/(rho0*AngularFrequency*(gamma-1.0d0)*CV)

      C1 = -SlipCoefficient1/(2.0d0-SlipCoefficient1) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      C2 = 1/kappa*SlipCoefficient2/(2.0d0-SlipCoefficient2) * &
          (gamma+1.0d0)/2.0d0 * rho0 * CV * &
          SQRT(2.0d0*(gamma-1.0d0)*CV*(T0+WallT0)/PI)

      C3 = -SlipCoefficient3/(2.0d0-SlipCoefficient3) * &
          rho0*SQRT(2.0d0*(gamma-1.0d0)*CV* &
          (T0+WallT0)/PI)

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

      CALL TangentDirections(Normal, Tangent1, Tangent2)

      DO p=1,n
        DO i=1,dim
          DO j=1,dim
            DO q=1,n
              LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                  LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  Basis(q) * Tangent1(j) * Basis(p) * Tangent1(i) * s
            END DO

            LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                Tangent1(j) * Basis(p) * Tangent1(i) * s 
          END DO
        END DO
      END DO

      IF (dim > 2) THEN
        DO p=1,n
          DO i=1,dim
            DO j=1,dim
              DO q=1,n
                LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                    LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                    CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    Basis(q) * Tangent2(j) * Basis(p) * Tangent2(i) * s
              END DO

              LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                  CMPLX(0.0d0, C1, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                  CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                  Tangent2(j) * Basis(p) * Tangent2(i) * s 
            END DO
          END DO
        END DO
      END IF

      !Normal direction

      DO p=1,n
         DO i=1,dim
            DO j=1,dim
               DO q=1,n
                  LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) = &
                       LSTIFF( (p-1)*(DIM+2)+i, (q-1)*(DIM+2)+j) + &
                       CMPLX(0.0d0, C3, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                       Basis(q) * Normal(j) * Basis(p) * Normal(i) * s
               END DO

               LFORCE( (p-1)*(DIM+2)+i) = LFORCE( (p-1)*(DIM+2)+i) + &
                    CMPLX(0.0d0, C3, kind=dp) * 1.0d0/(AngularFrequency*rho0) * &
                    CMPLX(ReV0(j),ImV0(j),kind=dp) * &
                    Normal(j) * Basis(p) * Normal(i) * s 
            END DO
         END DO
      END DO




      DO p=1,n
        DO q=1,n
          LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) = &
              LSTIFF( (p-1)*(DIM+2)+DIM+1, (q-1)*(DIM+2)+DIM+1) + &
              CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * &
              Basis(q) * Basis(p) * s
        END DO
        LFORCE( (p-1)*(DIM+2)+DIM+1) = LFORCE( (p-1)*(DIM+2)+DIM+1) + &
            CMPLX(0.0d0, K1, kind=dp) * CMPLX(C2, 0.0d0, kind=dp) * &
            WallT0 * Basis(p) * s
      END DO

    END DO  ! Loop over integration points

    DO p=1,n
      DO i=1,DIM+2
        Force( 2*(DIM+2)*(p-1)+2*i-1 ) = REAL(LFORCE( (p-1)*(DIM+2)+i ))
        Force( 2*(DIM+2)*(p-1)+2*i ) = AIMAG(LFORCE( (p-1)*(DIM+2)+i ))
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SlipMatrix
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
     END SELECT

!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalPreconditioner
!------------------------------------------------------------------------------






!------------------------------------------------------------------------------
   SUBROUTINE SetBoundaryConditions( Model, StiffMatrix, &
                   Name, DOF, NDOFs, Perm )
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
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,t,k1,k2
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
      n = CurrentElement % TYPE % NumberOfNodes
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
                s = 1.0d0 
                CALL ZeroRow( StiffMatrix,k )
                CALL SetMatrixElement( StiffMatrix,k,k, 1.0d0 * s )
              END IF
            END DO
          END IF
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetBoundaryConditions
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
  SUBROUTINE ComputeAcousticImpedance( AcousticI, Bndries, AcImpedances )
!------------------------------------------------------------------------------

    INTEGER :: AcousticI, Bndries(:)
    REAL(KIND=dp) :: AcImpedances(:,:)
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    REAL(KIND=dp), ALLOCATABLE :: PressureData(:,:), VelocityData(:)
    REAL(KIND=dp) :: ElementPresRe, ElementPresIm
    REAL(KIND=dp) :: ElementVeloRe, ElementVeloIm
    REAL(KIND=dp) :: u, v, w, s, xpos, ypos, zpos
    REAL(KIND=dp) :: SqrtMetric, SqrtElementMetric
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes), Metric(3,3), &
         dBasisdx(Model % MaxElementNodes,3), Symb(3,3,3), &
         ddBasisddx(Model % MaxElementNodes,3,3), dSymb(3,3,3,3)
    REAL(KIND=dp) :: Normal(3), VReal(3), VIm(3), Preal, PIm
    REAL(KIND=dp) :: AbsNormalVelo, InPhasePres, OutPhasePres
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i, j, k, n, t
    LOGICAL :: stat
!------------------------------------------------------------------------------

    ALLOCATE( PressureData( AcousticI, 2 ) )
    PressureData = 0.0d0
    ALLOCATE( VelocityData( 2 ) )
    VelocityData = 0.0d0

    VReal = 0.0d0
    VIm = 0.0d0

    DO t = Model % NumberOfBulkElements +1, &
         Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

       CurrentElement => Solver % Mesh % Elements(t)
       n = CurrentElement % TYPE % NumberOfNodes
       NodeIndexes => CurrentElement % NodeIndexes

       DO i = 1, AcousticI
          IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(Bndries(i)) % Tag ) THEN

             ElementPresRe = 0.0d0
             ElementPresIm = 0.0d0
             ElementVeloRe = 0.0d0
             ElementVeloIm = 0.0d0

             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------

             IntegStuff = GaussPoints( CurrentElement )
             DO j = 1, IntegStuff % n
                u = IntegStuff % u(j)
                v = IntegStuff % v(j)
                w = IntegStuff % w(j)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
               
                stat = ElementInfo( CurrentElement, ElementNodes, u, v, w, &
                     SqrtElementMetric, Basis, dBasisdx, ddBasisddx, &
                     .FALSE., .FALSE. )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
                s = 1.0d0
                IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                   xpos = SUM( ElementNodes % x(1:n) * Basis(1:n) )
                   ypos = SUM( ElementNodes % y(1:n) * Basis(1:n) )
                   zpos = SUM( ElementNodes % z(1:n) * Basis(1:n) )
                   s = 2*PI
                END IF
         
                CALL CoordinateSystemInfo( Metric, SqrtMetric, Symb, dSymb, &
                     xpos, ypos, zpos)
 
                s = s * SqrtMetric * SqrtElementMetric * IntegStuff % s(j)
        
!------------------------------------------------------------------------------

                PReal = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes) ) &
                     - 1 ) * Basis(1:n) )
                PIm = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes) ) ) &
                     * Basis(1:n) )


                ElementPresRe = ElementPresRe + s * PReal
                ElementPresIm = ElementPresIm + s * PIm

!------------------------------------------------------------------------------
!   The normal velocity is needed only for the transmitting boundary i = 1
!------------------------------------------------------------------------------

                IF ( i == 1 ) THEN
                   Normal = Normalvector( CurrentElement, ElementNodes, &
                        u, v, .TRUE. )
                   Normal = -Normal   ! use inward normal vector

                   DO k = 1, VelocityComponents
                      VReal(k) = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes)-1 ) &
                           + 2 * k - 1 ) * Basis(1:n) )
                      VIm(k) = SUM( Flow( DOFs * ( FlowPerm(NodeIndexes)-1 ) &
                           + 2 * k ) * Basis(1:n) )
                   END DO

                   ElementVeloRe = ElementVeloRe + s * &
                        SUM( VReal(1:3) * Normal(1:3) )
                   ElementVeloIm = ElementVeloIm + s * &
                        SUM( VIm(1:3) * Normal(1:3) )
                END IF

!------------------------------------------------------------------------------

             END DO

          PressureData(i,1) = PressureData(i,1) + ElementPresRe
          PressureData(i,2) = PressureData(i,2) + ElementPresIm
          VelocityData(1) = VelocityData(1) + ElementVeloRe
          VelocityData(2) = VelocityData(2) + ElementVeloIm

          END IF

       END DO
    END DO

!------------------------------------------------------------------------------

    DO i = 1, AcousticI
       AcImpedances(i,1) = ( PressureData(i,1) * VelocityData(1) + &
            PressureData(i,2) * VelocityData(2) ) / &
            SUM( VelocityData**2 )
       AcImpedances(i,2) = ( PressureData(i,2) * VelocityData(1) - &
            PressureData(i,1) * VelocityData(2) ) / &
            SUM( VelocityData**2 )
    END DO

    DEALLOCATE( PressureData )
    DEALLOCATE( VelocityData )

!------------------------------------------------------------------------------
  END SUBROUTINE ComputeAcousticImpedance
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ComputeAverageVelocity( Element, n, WallVelo, res)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    REAL(KIND=dp) :: WallVelo(:,:)
    COMPLEX(KIND=dp) :: res   
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    INTEGER :: t
    REAL(KIND=dp) :: AK, SqrtElementMetric,U,V,W,S 
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: Normal(3)
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    COMPLEX(KIND=dp) :: Velo(3), NormVelo
!------------------------------------------------------------------------------
    res = CMPLX(0.0d0,0.0d0, kind=dp)
    AK = ElementArea(Solver % Mesh, Element, n)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )

      s = s * SqrtElementMetric

      Normal = Normalvector(Element, Nodes, U, V, .TRUE.)
      !PRINT *, Normal

      Velo(1) = CMPLX( SUM( WallVelo(1,1:n) * Basis(1:n) ), SUM( WallVelo(2,1:n) * Basis(1:n) ), kind=dp )
      Velo(2) = CMPLX( SUM( WallVelo(3,1:n) * Basis(1:n) ), SUM( WallVelo(4,1:n) * Basis(1:n) ), kind=dp )
      Velo(3) = CMPLX( SUM( WallVelo(5,1:n) * Basis(1:n) ), SUM( WallVelo(6,1:n) * Basis(1:n) ), kind=dp )
      NormVelo = Velo(1) * CMPLX(Normal(1), 0.0d0, kind=dp) + Velo(2) * CMPLX(Normal(2), 0.0d0, kind=dp) + &
          Velo(3) * CMPLX(Normal(3), 0.0d0, kind=dp)
      res = res + s * NormVelo
!------------------------------------------------------------------------------
    END DO
    res = res / CMPLX(AK, 0.0d0, kind=dp)
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeAverageVelocity
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
  SUBROUTINE ExplicitStabilisationMatrix(  StiffMatrix, Force, AngularFrequency , &
      SpecificHeat, HeatRatio, Density,               &
      Temperature, Conductivity, Viscosity, Lambda,             &
      HeatSource, Load, Bubbles, Element, n, Nodes, Dofs, NodeOnBoundary )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), AngularFrequency, &
        SpecificHeat(:), HeatRatio(:), Density(:),    &       
        Temperature(:), Conductivity(:), Viscosity(:), Lambda(:),  &
        HeatSource(:,:), Load(:,:)
    LOGICAL :: Bubbles, NodeOnBoundary(:)
    INTEGER :: n, Dofs
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*n), dBasisdx(2*n,3), ddBasisddx(n,3,3)
    REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, L(6), &
        CV, gamma, rho0, P0, T0, kappa, mu, la, f1, f2, K1, K2, K3, k  
    COMPLEX(KIND=dp) :: LSTIFF(n*(Dofs-2),n*(Dofs-2)), LFORCE(n*(Dofs-2)), A, &
        SchurConst, C1, C2, C3

    INTEGER :: i, j, p, q, t, DIM, NBasis, CoordSys, VelocityDofs, & 
        VelocityComponents
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    REAL(KIND=dp) :: X, Y, Z, Metric(3,3), SqrtMetric, Symb(3,3,3), &
        dSymb(3,3,3,3), StabTerms(n), AK
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
   
    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

    LSTIFF = CMPLX(0.0d0,0.0d0, kind=dp)
    LFORCE = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    NBasis = n
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., Bubbles )
      s = s * SqrtElementMetric
!------------------------------------------------------------------------------
!     Problem parameters and the real and imaginary part of the 
!     load at the integration point
!------------------------------------------------------------------------------
      CV = SUM( SpecificHeat(1:n) * Basis(1:n) )       
      gamma = SUM( HeatRatio(1:n) * Basis(1:n) )
      rho0 = SUM( Density(1:n) * Basis(1:n) )
      T0 = SUM( Temperature(1:n) * Basis(1:n) )
      kappa = SUM( Conductivity(1:n) * Basis(1:n) )
      mu = SUM( Viscosity(1:n) * Basis(1:n) )
      la = SUM( Lambda(1:n) * Basis(1:n) )

      P0 = (gamma-1.0d0)*CV*rho0*T0
      C1 = CMPLX( 1.0d0,AngularFrequency/P0*(2.0d0*mu+la), kind=dp ) / &
          CMPLX( 1.0d0,AngularFrequency/P0*(la), kind=dp )  
      C2 = CMPLX( rho0*AngularFrequency,0.0d0, kind=dp) / &
          CMPLX( P0/AngularFrequency, la, kind=dp )
      C3 = CMPLX( 1.0d0, 0.0d0, kind=dp)

      DO p=1,N
        IF ( NodeOnBoundary(p) ) CYCLE
        DO q=1,N
          DO i=1,dim
            
            LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+DIM+1 ) = &
                LSTIFF( p*(DIM+2), (q-1)*(DIM+2)+DIM+1 ) - &
                s * C3 * dBasisdx(q,i) * dBasisdx(p,i)

            LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
                LSTIFF( p*(DIM+2), q*(DIM+2) ) - &
                s * C1 * dBasisdx(q,i) * dBasisdx(p,i)

          END DO
          
          LSTIFF( p*(DIM+2), q*(DIM+2) ) = &
              LSTIFF( p*(DIM+2), q*(DIM+2) ) + &
              s * C2 * Basis(q) * Basis(p)
        END DO
      END DO

    END DO


    DO p=1,n
      DO i=1,DIM+2
        DO q=1,n
          DO j=1,DIM+2
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i-1, 2*(DIM+2)*(q-1)+2*j ) =  &
                -AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j-1 ) =  &
                AIMAG( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
            StiffMatrix( 2*(DIM+2)*(p-1)+2*i, 2*(DIM+2)*(q-1)+2*j ) =  &
                REAL( LSTIFF((DIM+2)*(p-1)+i,(DIM+2)*(q-1)+j) )
          END DO
        END DO
      END DO
    END DO
    Force = 0.0d0
!--------------------------------------------------------------------------------   
  END SUBROUTINE ExplicitStabilisationMatrix
!--------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
  SUBROUTINE SurfaceForceIntegration(Element, Parent, Traction, Moment, MomentAbout, Area, &
     CalculateMoment, Velo, Pressure, Viscosity, BulkViscosity, Nodes, ParentNodes, np)
!----------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element, Parent  
  REAL(kind=dp) :: Traction(6), Moment(6), MomentAbout(3), Area, Velo(:,:), &
      Pressure(:,:), Viscosity(:), BulkViscosity(:)
  LOGICAL :: CalculateMoment
  TYPE(Nodes_t) :: Nodes, ParentNodes
  INTEGER :: np
!----------------------------------------------------------------------------------
  LOGICAL :: stat
  INTEGER :: N_Integ, t, i, dim
  REAL(kind=dp) :: u, v, w, detJ, s, Basis(np), dBasisdx(np,3), Normal(3), &
      Visc, Lambda, ReGrad(3,3), ImGrad(3,3), ReDiv, ImDiv, ReD(3,3), ImD(3,3), &
      ReP, ImP, tmpmat(3,np), tmp = 0.0d0, r1, r2, r3, r
  REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  SAVE tmp
!----------------------------------------------------------------------------------
  dBasisdx = 0.0d0
  ReDiv = 0.0d0
  ImDiv = 0.0d0
  Traction = 0.0d0
  Area = 0.0d0
  dim = CoordinateSystemDimension()
!----------------------------------------------------------------------------------
! The normal vector to the boundary element
!----------------------------------------------------------------------------------
  SELECT CASE( Element % TYPE % NumberOfNodes )
  CASE(2)
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.)    
  CASE( 3 )
    Normal = Normalvector(Element, Nodes, 0.3d0, 0.3d0, .TRUE.)          
  CASE( 4 )
    Normal = Normalvector(Element, Nodes, 0.0d0, 0.0d0, .TRUE.) 
  END SELECT


!----------------------------------------------------------------------------------
! Compute field derivatives at the integration points of the parent element
!----------------------------------------------------------------------------------
  IntegStuff = Gausspoints(Parent)
  U_Integ => IntegStuff % u
  V_Integ => IntegStuff % v
  W_Integ => IntegStuff % w
  S_Integ => IntegStuff % s
  N_Integ =  IntegStuff % n 

!------------------------------------------------------------------------------
  DO t=1,N_Integ
!------------------------------------------------------------------------------
    u = U_Integ(t)
    v = V_Integ(t)
    w = W_Integ(t)
!------------------------------------------------------------------------------
!   Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
    stat = ElementInfo( Parent, ParentNodes, u, v, w, &
        detJ, Basis, dBasisdx )

    s = detJ * S_Integ(t)
 
    IF (CoordSys == AxisSymmetric) THEN
       r = SUM( Basis(1:np) * ParentNodes % x(1:np) )
       s = r * s
    END IF

    Visc = SUM( Viscosity(1:np) * Basis(1:np) ) 
!    Visc = 0.0d0
    Lambda = SUM( BulkViscosity(1:np) * Basis(1:np) )
!    Lambda = 0.0d0
    ReP = SUM(Pressure(1,1:np) * Basis(1:np))
    ImP = SUM(Pressure(2,1:np) * Basis(1:np))
    tmpmat(1,1:np) = Velo(1,1:np)
    tmpmat(2,1:np) = Velo(3,1:np)
    tmpmat(3,1:np) = Velo(5,1:np)    
    ReGrad = MATMUL( tmpmat, dBasisdx )
    tmpmat(1,1:np) = Velo(2,1:np)
    tmpmat(2,1:np) = Velo(4,1:np)
    tmpmat(3,1:np) = Velo(6,1:np)    
    ImGrad = MATMUL( tmpmat, dBasisdx )

    ReD = Visc * ( ReGrad + TRANSPOSE(ReGrad) )
    ImD = Visc * ( ImGrad + TRANSPOSE(ImGrad) )    
    
    ReDiv = 0.0d0
    ImDiv = 0.0d0
    DO i = 1,dim
      ReDiv = ReDiv + ReGrad(i,i)
      ImDiv = ImDiv + ImGrad(i,i)
    END DO

    IF (CoordSys == AxisSymmetric) THEN
       ReDiv = ReDiv + 1/r * SUM( Velo(1,1:np) * Basis(1:np) )
       ImDiv = ImDiv + 1/r * SUM( Velo(2,1:np) * Basis(1:np) )
    END IF

    Traction(1:3) = Traction(1:3) - ReP * Normal(1:3) + Lambda * &
        ReDiv * Normal(1:3) + MATMUL( ReD, Normal)

    Traction(4:6) = Traction(4:6) - ImP * Normal(1:3) + Lambda * &
        ImDiv * Normal(1:3) + MATMUL( ImD, Normal)

  END DO
  
  Area = ElementArea(Solver % Mesh, Element, Element % TYPE % NumberOfNodes)

  IF (CoordSys == AxisSymmetric) Area = 2*pi*Area

  Traction = -Area/N_Integ * Traction


  IF (CalculateMoment) THEN
    r1 = SUM( Nodes % x(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(1)
    r2 = SUM( Nodes % y(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(2)
    r3 = SUM( Nodes % z(1:Element % TYPE % NumberOfNodes) )/Element % TYPE % NumberOfNodes - MomentAbout(3)
    Moment(1) = r2 * Traction(3) - r3 * Traction(2)
    Moment(2) = r3 * Traction(1) - r1 * Traction(3)
    Moment(3) = r1 * Traction(2) - r2 * Traction(1)
    Moment(4) = r2 * Traction(6) - r3 * Traction(5)
    Moment(5) = r3 * Traction(4) - r1 * Traction(6)
    Moment(6) = r1 * Traction(5) - r2 * Traction(4)
  END IF

!----------------------------------------------------------------------------------
END SUBROUTINE SurfaceForceIntegration
!-----------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
SUBROUTINE SurfaceImpedanceIntegration(Element, Velo, Pressure, Nodes, n, Impedance1, Impedance2, &
     C3, C4, Area)
!---------------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element  
  REAL(kind=dp) :: Velo(:,:), Pressure(:,:), Area
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n
  COMPLEX(kind=dp) :: Impedance1, Impedance2, C3, C4
!------------------------------------------------------------------------------
  REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, &
       Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), Normal(3), r, ReP, ImP, &
       u1, v1, u2, v2, u3, v3, un, vn

  LOGICAL :: Stat
  INTEGER :: i, j, p, q, t, DIM
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
!  Impedance1 = CMPLX(0.0d0,0.0d0, kind=dp)
!  Impedance2 = CMPLX(0.0d0,0.0d0, kind=dp)
!------------------------------------------------------------------------------
! Numerical integration
!------------------------------------------------------------------------------
  IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
  DO t=1, IntegStuff % n
     U = IntegStuff % u(t)
     V = IntegStuff % v(t)
     W = IntegStuff % w(t)
     S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
     s = s * SqrtElementMetric
     IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
     END IF

     Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

     ReP = SUM(Pressure(1,1:n) * Basis(1:n))
     ImP = SUM(Pressure(2,1:n) * Basis(1:n))
      
     u1 = SUM(Velo(1,1:n) * Basis(1:n))
     v1 = SUM(Velo(2,1:n) * Basis(1:n))
      
     u2 = SUM(Velo(3,1:n) * Basis(1:n))
     v2 = SUM(Velo(4,1:n) * Basis(1:n))
      
     u3 = SUM(Velo(5,1:n) * Basis(1:n))
     v3 = SUM(Velo(6,1:n) * Basis(1:n))

     un = u1 * Normal(1) + u2 * Normal(2) + u3 * Normal(3)
     vn = v1 * Normal(1) + v2 * Normal(2) + v3 * Normal(3)

     Impedance1 = Impedance1 + CMPLX(ReP, ImP, kind=dp) * CMPLX(ReP, ImP, kind=dp) * s
     Impedance2 = Impedance2 - CMPLX(ReP, ImP, kind=dp) * CMPLX(un, vn, kind=dp) * s

     C3 = C3 - CMPLX(ReP, ImP, kind=dp) * CMPLX(un, vn, kind=dp) * s
     C4 = C4 + CMPLX(un, vn, kind=dp) * CMPLX(un, vn, kind=dp) * s

    
     Area = Area + s
   END DO

!-----------------------------------------------------------------------------
END SUBROUTINE SurfaceImpedanceIntegration
!-----------------------------------------------------------------------------



!------------------------------------------------------------------------------
FUNCTION DiscontIndexes( Solver, Element, n ) RESULT(GapIndexes)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    INTEGER :: GapIndexes(n)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent, Left, Right
    INTEGER :: i,j,k,mini
    REAL(KIND=dp) :: x0,y0,z0,dx,dy,dz,dist,mindist
!------------------------------------------------------------------------------
    Left  => Element % BoundaryInfo % Left
    Right => Element % BoundaryInfo % Right

    GapIndexes(1:n) = 0
    IF ( .NOT. ASSOCIATED(Left) .OR. .NOT. ASSOCIATED(Right) ) RETURN
    
    DO i=1,n
      Parent => Left
      k = Element % NodeIndexes(i)
      IF ( ANY( Parent % NodeIndexes == k ) ) Parent => Right
      
      x0 = Solver % Mesh % Nodes % x(k) 
      y0 = Solver % Mesh % Nodes % y(k) 
      z0 = Solver % Mesh % Nodes % z(k) 
      
      mini = 0
      mindist = HUGE(mindist)
      DO j=1,Parent % TYPE % NumberOfNodes
        k = Parent % NodeIndexes(j)
        dx = Solver % Mesh % Nodes % x(k) - x0
        dy = Solver % Mesh % Nodes % y(k) - y0
        dz = Solver % Mesh % Nodes % z(k) - z0
        dist = dx**2 + dy**2 + dz**2 
        IF( dist < mindist) THEN
          mini = k
          mindist = dist
        END IF
      END DO

      IF(mini /= Element % NodeIndexes(i)) GapIndexes(i) = mini
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DiscontIndexes
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE SetDirichletPoints( StiffMatrix, ForceVector, DOF, NDOFs, &
      Perm, n, NodeIndexes, NodeValues) 
!------------------------------------------------------------------------------

    INTEGER :: n
    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:), NodeValues(:)
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndexes(:)
!------------------------------------------------------------------------------

    INTEGER :: i,PermIndex
    REAL(KIND=dp) :: s

!------------------------------------------------------------------------------

    DO i=1,n
      PermIndex = Perm(NodeIndexes(i))
      IF ( PermIndex == 0 ) CYCLE
      
      PermIndex = NDOFs * (PermIndex-1) + DOF
      
      IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN
        
        CALL SBand_SetDirichlet( StiffMatrix,ForceVector,PermIndex,NodeValues(i) )
        
      ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. &
          StiffMatrix % Symmetric ) THEN
        
        CALL CRS_SetSymmDirichlet(StiffMatrix,ForceVector,PermIndex,NodeValues(i) )
        
      ELSE                  
        
        s = StiffMatrix % Values(StiffMatrix % Diag(PermIndex))
        ForceVector(PermIndex) = NodeValues(i) * s
        CALL ZeroRow( StiffMatrix,PermIndex )
        CALL SetMatrixElement( StiffMatrix,PermIndex,PermIndex,1.0d0*s )
        
      END IF
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoints
!------------------------------------------------------------------------------


SUBROUTINE AcousticShellInterface()

  INTEGER :: BC,t,n
  TYPE(Variable_t), POINTER :: Var
  REAL (KIND=dp), POINTER :: Deflection(:)
  INTEGER, POINTER :: DeflectionPerm(:)
  INTEGER :: Indexes(MAX_ELEMENT_NODES)
  REAL (KIND=dp) :: NodeValues(MAX_ELEMENT_NODES)
  INTEGER :: sign
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: AcousticShell, AcousticNoslip, AcousticIsothermal, DoneDiscont

  Var => VariableGet( Model % Variables, 'U' )
  IF(.NOT. ASSOCIATED(Var)) THEN
    CALL Warn('Acoustics','Could not find variable U for FSI simulation')
    RETURN
  END IF
  Deflection => Var % Values
  DeflectionPerm => Var % Perm

  DO BC=1,Model % NumberOfBCs

    AcousticShell = GetLogical( Model % BCs(BC) % Values,'Solid Interaction',GotIt)
    AcousticNoslip = GetLogical( Model % BCs(BC) % Values,'Acoustic Wall Noslip',GotIt)
    AcousticIsothermal = GetLogical( Model % BCs(BC) % Values,'Acoustic Wall Isothermal',GotIt)

    IF(.NOT. (AcousticShell .OR. AcousticNoslip .OR. AcousticIsothermal)) CYCLE

    DO t = Model % NumberOfBulkElements + 1, &
        Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
      
      Element => Model % Elements(t)
      IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
      
      Model % CurrentElement => Element
      n = Element % TYPE % NumberOfNodes
      Indexes(1:n) = Element % NodeIndexes
      
      IF( AcousticShell) THEN
        GapIndexes(1:n) = Indexes(1:n)       
        !DoneDiscont = .FALSE.
        
        DO i=1,VelocityComponents
          ! Real part 
          NodeValues(1:n) = -AngularFrequency * Deflection( 6 * (DeflectionPerm(Indexes(1:n)) - 1) + 2*i)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i-1, 10, FlowPerm, n, GapIndexes, NodeValues) 
          ! Im part 
          NodeValues(1:n) = AngularFrequency * Deflection( 6 * (DeflectionPerm(Indexes(1:n)) - 1) + 2*i-1)
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i, 10, FlowPerm, n, GapIndexes, NodeValues) 
        END DO
        
        IF(AcousticIsothermal) THEN
          NodeValues(1:n) = 0.0
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-3, Dofs, FlowPerm, n, GapIndexes, NodeValues) 
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-2, Dofs, FlowPerm, n, GapIndexes, NodeValues) 
        END IF
        
        !IF( .NOT. DoneDiscont ) THEN
        !  DoneDiscont = .TRUE.
        !  GapIndexes(1:n) = DiscontIndexes(Solver, Element, n)
        !  IF( ALL(GapIndexes(1:n) /= 0)) THEN
        !    IF(.FALSE.) PRINT *,'found discontinuity',n,GapIndexes(1:n)
        !    GOTO 100
        !  END IF
        !END IF

      ELSE IF (AcousticNoslip) THEN
        NodeValues(1:n) = 0.0

        DO i=1, VelocityComponents
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i-1, Dofs, FlowPerm, n, Indexes, NodeValues) 
          CALL SetDirichletPoints( StiffMatrix, ForceVector, 2*i, Dofs, FlowPerm, n, Indexes, NodeValues) 
        END DO
        IF(AcousticIsothermal) THEN
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-3, Dofs, FlowPerm, n, Indexes, NodeValues) 
          CALL SetDirichletPoints( StiffMatrix, ForceVector, Dofs-2, Dofs, FlowPerm, n, Indexes, NodeValues) 
        END IF
      END IF
    END DO
  END DO

END SUBROUTINE AcousticShellInterface




!---------------------------------------------------------------------------------------
SUBROUTINE FSIIntegration(Element, Velo, Pressure, Nodes, n, uf, wf, Area)
!---------------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element  
  REAL(kind=dp) :: Velo(:,:), Pressure(:,:), Area
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n
  COMPLEX(kind=dp) :: uf, wf
!------------------------------------------------------------------------------
  REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S, &
       Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), Normal(3), r, ReP, ImP, &
       u1, v1, u2, v2, u3, v3, un, vn

  LOGICAL :: Stat
  INTEGER :: i, j, p, q, t, DIM
  TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Numerical integration
!------------------------------------------------------------------------------
  IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
  DO t=1, IntegStuff % n
     U = IntegStuff % u(t)
     V = IntegStuff % v(t)
     W = IntegStuff % w(t)
     S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!    Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
     stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, ddBasisddx, .FALSE., .FALSE. )
     s = s * SqrtElementMetric
     IF (CoordSys == AxisSymmetric) THEN
        r = SUM( Basis * Nodes % x(1:n) )
        s = r * s
     END IF

     Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

     ReP = SUM(Pressure(1,1:n) * Basis(1:n))
     ImP = SUM(Pressure(2,1:n) * Basis(1:n))
      
     u1 = SUM(Velo(1,1:n) * Basis(1:n))
     v1 = SUM(Velo(2,1:n) * Basis(1:n))
      
     u2 = SUM(Velo(3,1:n) * Basis(1:n))
     v2 = SUM(Velo(4,1:n) * Basis(1:n))
      
     u3 = SUM(Velo(5,1:n) * Basis(1:n))
     v3 = SUM(Velo(6,1:n) * Basis(1:n))

     un = u1 * Normal(1) + u2 * Normal(2) + u3 * Normal(3) 
     vn = v1 * Normal(1) + v2 * Normal(2) + v3 * Normal(3) 


     ! U = \int (u.u*)
     uf = uf + (s / AngularFrequency ** 2) * ( CMPLX(u1, v1, kind=dp) * CMPLX(u1, -v1, kind=dp) + &
          CMPLX(u2, v2, kind=dp) * CMPLX(u2, -v2, kind=dp) + &
          CMPLX(u3, v3, kind=dp) * CMPLX(u3, -v3, kind=dp) )         
     
     ! Kf = \int p.u* 
     wf = wf + CMPLX(0.0d0, 1.0d0, kind=dp) * (s / AngularFrequency) * CMPLX(ReP, ImP, kind=dp) * &
          CMPLX(un, -vn, kind=dp)
    
     Area = Area + s
   END DO

!-----------------------------------------------------------------------------
END SUBROUTINE FSIIntegration
!-----------------------------------------------------------------------------









!------------------------------------------------------------------------------
END SUBROUTINE AcousticsSolver
!------------------------------------------------------------------------------
