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
! *  ELMER/FEM Data Interpolation Main Program (File to file)
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Ville Savolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \ingroup Programs 
!> \{

!> \defgroup ResultToResult Program ResultToResult
!> \{

!------------------------------------------------------------------------------
!>  Stand-alone program for Elmer results file to Elmer results processing 
!> file conversion. The meshes may be different for the two result files.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
PROGRAM ResultToResult

!------------------------------------------------------------------------------
     USE Types
     USE Lists
     USE Integration
     USE Interpolation

     USE MainUtils

     USE CoordinateSystems
     USE ModelDescription
     USE ElementDescription

     USE SParIterSolve
     USE SParIterComm
     USE SParIterGlobals

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------


     INTEGER :: i,j,k,n,l,t,k1,k2,iter,Ndeg,Time,istat,nproc

     REAL(KIND=dp) :: s,dt

     REAL(KIND=dp), TARGET :: SimulationTime(1)

     REAL(KIND=dp), POINTER ::  TimeVariable(:)

     LOGICAL :: GotIt,TransientSimulation,LastSaved

     INTEGER :: TimeIntervals,interval,timestep,SavedSteps,TimeStepMaxIter

     REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
     INTEGER, POINTER :: Timesteps(:),OutputIntervals(:)

     LOGICAL :: KEModelSolved = .FALSE., SteadyStateReached = .FALSE.

     TYPE(ElementType_t),POINTER :: elmt

     TYPE(ParEnv_t), POINTER :: ParallelEnv

     CHARACTER(LEN=MAX_NAME_LEN) :: OldModelName
     CHARACTER(:), ALLOCATABLE :: NewModelName, eq, OutputName, OutputFile, &
            OldOutputFile, PostFile, RestartFile


     TYPE(Model_t), POINTER :: OldModel, NewModel
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: Solver

! Note CurrentModel is a global variable defined in Types.f90
!    TYPE(Model_t),  POINTER :: CurrentModel

     INTEGER :: dim

     REAL(KIND=dp), DIMENSION(6) :: BoundingBox
     TYPE(Quadrant_t), POINTER :: RootQuadrant
     LOGICAL :: QuadrantTreeExists=.FALSE.

!------------------------------------------------------------------------------
!    Read input file name (the old model .sif file)
!------------------------------------------------------------------------------

     OPEN( 1,file='ELMERSOLVER_STARTINFO')
       READ(1,'(a)') OldModelName
       READ(1,*) nproc
     CLOSE(1)

     NewModelName = 'Interpolation.sif'

!------------------------------------------------------------------------------
!    If parallel execution requested, initialize parallel environment
!------------------------------------------------------------------------------
     IF ( nproc > 1 ) THEN
       ParallelEnv => ParCommInit()
     ELSE
       ParEnv % PEs  = 1 
       ParEnv % MyPE = 0
     END IF

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     CALL InitializeElementDescriptions
!------------------------------------------------------------------------------
!    Read Model from Elmer mesh data base
!------------------------------------------------------------------------------
     IF ( ParEnv % MyPE == 0 ) THEN
       PRINT*, ' '
       PRINT*, ' '
       PRINT*, '-----------------------'
       PRINT*,'Reading Old Model ...       '
     END IF

     OldModel => LoadModel( OldModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE )
     CALL SetCurrentMesh( OldModel, OldModel % Meshes )

     IF ( ParEnv % MyPE == 0 ) THEN
       PRINT*,'... Done               '
       PRINT*, '-----------------------'
     END IF

     PRINT*, ' '
     PRINT*, ' '
     PRINT*, '-----------------------'
     PRINT*,'Reading New Model ...       '
     NewModel => LoadModel( NewModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE )
     PRINT*,'... Done               '
     PRINT*, '-----------------------'
     CALL SetCurrentMesh( NewModel, NewModel % Meshes )

!------------------------------------------------------------------------------
! Old Model
!------------------------------------------------------------------------------
     eq = ListGetString( OldModel % Simulation, 'Simulation Type' )
     TransientSimulation = .FALSE.

     IF ( eq == 'transient' ) THEN
       TransientSimulation= .TRUE.

       elmt => GetElementType( 303 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 3

       elmt => GetElementType( 504 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 4
     END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Figure out requested coordinate system
!------------------------------------------------------------------------------
     eq = ListGetString( OldModel % Simulation, 'Coordinate System', GotIt )

     Coordinates = Cartesian
     IF ( eq == 'cartesian 2d' ) THEN

       OldModel % Dimension = 2
       Coordinates = Cartesian

     ELSE IF ( eq == 'cartesian 3d' ) THEN

       OldModel % Dimension = 3
       Coordinates = Cartesian

     ELSE IF ( eq == 'axi symmetric' ) THEN

       OldModel % Dimension = 2
       Coordinates = AxisSymmetric

     ELSE IF( eq == 'cylindric symmetric' ) THEN

       OldModel % Dimension = 2
       Coordinates = CylindricSymmetric

     ELSE IF( eq == 'cylindrical' ) THEN

       OldModel % Dimension = 3
       Coordinates = Cylindric

     ELSE IF( eq == 'polar 2d' ) THEN

       OldModel % Dimension = 2
       Coordinates = Polar

     ELSE IF( eq == 'polar 3d' ) THEN

       OldModel % Dimension = 3
       Coordinates = Polar

     ELSE

       PRINT*,'Solver: ERROR: Unknown global coordinate system: ',TRIM(eq),' Aborting'
       STOP

     END IF

!------------------------------------------------------------------------------
!   Figure out what (flow,heat,stress,...) should be computed, and get
!   memory for the dofs
!------------------------------------------------------------------------------
     DO i=1,OldModel % NumberOfSolvers
       eq = ListGetString( OLdModel % Solvers(i) % Values,'Equation' )

       Solver => OldModel % Solvers(i)
       CALL AddEquationBasics( Solver, eq, TransientSimulation )
       CALL AddEquationSolution( Solver, TransientSimulation )
     END DO

!------------------------------------------------------------------------------
!    Add coordinates to list of variables so that coordinate dependent
!    parameter computing routines can ask for them...
!------------------------------------------------------------------------------
     Mesh => OldModel % Meshes 
     NULLIFY( Solver )
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 1',1,Mesh % Nodes % x )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 2',1,Mesh % Nodes % y )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 3',1,Mesh % Nodes % z )
       Mesh => Mesh % Next
     END DO

!------------------------------------------------------------------------------
!    Time integration and/or steady state steps
!------------------------------------------------------------------------------
     IF ( TransientSimulation ) THEN
       Timesteps => ListGetIntegerArray( OldModel % Simulation, &
                 'Timestep Intervals', GotIt )

       IF ( .NOT.GotIt ) THEN
         PRINT*,'Solver Input error: Time step intervals must be defined.'
         STOP
       END IF 

       TimestepSizes => ListGetConstRealArray( OldModel % Simulation, &
                             'Timestep Sizes', GotIt )

       IF ( .NOT.GotIt ) THEN
         PRINT*,'Solver Input error: Time step sizes must be defined.'
         STOP
       END IF 

       TimeIntervals = SIZE(Timesteps)

       TimestepMaxIter = ListGetInteger( OldModel % Simulation, &
                  'Steady State Max Iterations' )
!------------------------------------------------------------------------------
     ELSE
!------------------------------------------------------------------------------
!      Steady state
!------------------------------------------------------------------------------
       ALLOCATE( Timesteps(1) )

       Timesteps(1) = ListGetInteger( OldModel % Simulation, &
                  'Steady State Max Iterations' )

       ALLOCATE( TimestepSizes(1,1) )
       TimestepSizes(1,1) = 1.0D0

       TimeIntervals   = 1
       TimeStepMaxIter = 1
     END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Get Output File Options
!------------------------------------------------------------------------------
     OldOutputFile = ListGetString(OldModel % Simulation,'Output File',GotIt)
     IF ( .NOT.GotIt ) OldOutputFile = 'Result.dat'

     OutputIntervals => ListGetIntegerArray( OldModel % Simulation, &
                     'Output Intervals' )


!------------------------------------------------------------------------------
!    Add simulation time to list of system variables so that parameters
!    depending on time  may be set.
!------------------------------------------------------------------------------
     TimeVariable => SimulationTime
     Mesh => OldModel % Meshes
     TimeVariable => SimulationTime
     DO WHILE( ASSOCIATED(Mesh) )
       CALL VariableAdd( Mesh % Variables,Mesh,Solver,'Time',1,TimeVariable )
       Mesh => Mesh % Next
     END DO

!------------------------------------------------------------------------------
!    Check if we are restarting
!------------------------------------------------------------------------------

     RestartFile = ListGetString( OldModel % Simulation, &
                 'Restart File', GotIt )

     IF ( GotIt ) THEN
       k = ListGetInteger( OldModel % Simulation,'Restart Position',GotIt )

       Mesh => OldModel % Meshes
       DO WHILE( ASSOCIATED(Mesh) ) 
         IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
           OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
         ELSE
           OutputName = TRIM(RestartFile)
         END IF
         CALL SetCurrentMesh( OldModel, Mesh )
         CALL LoadRestartFile( OutputName,k,Mesh )
         Mesh => Mesh % Next
       END DO
     END IF

!------------------------------------------------------------------------------
! New Model
!------------------------------------------------------------------------------
     eq = ListGetString( NewModel % Simulation, 'Simulation Type' )
     TransientSimulation = .FALSE.

     IF ( eq == 'transient' ) THEN
       TransientSimulation= .TRUE.

       elmt => GetElementType( 303 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 3

       elmt => GetElementType( 504 )
       IF ( elmt % GaussPoints == 1 ) elmt % GaussPoints = 4
     END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Figure out requested coordinate system
!------------------------------------------------------------------------------
     eq = ListGetString( NewModel % Simulation, 'Coordinate System', GotIt )

     Coordinates = Cartesian
     IF ( eq == 'cartesian 2d' ) THEN

       NewModel % Dimension = 2
       Coordinates = Cartesian

     ELSE IF ( eq == 'cartesian 3d' ) THEN

       NewModel % Dimension = 3
       Coordinates = Cartesian

     ELSE IF ( eq == 'axi symmetric' ) THEN

       NewModel % Dimension = 2
       Coordinates = AxisSymmetric

     ELSE IF( eq == 'cylindric symmetric' ) THEN

       NewModel % Dimension = 2
       Coordinates = CylindricSymmetric

     ELSE IF( eq == 'cylindrical' ) THEN

       NewModel % Dimension = 3
       Coordinates = Cylindric

     ELSE IF( eq == 'polar 2d' ) THEN

       NewModel % Dimension = 2
       Coordinates = Polar

     ELSE IF( eq == 'polar 3d' ) THEN

       NewModel % Dimension = 3
       Coordinates = Polar

     ELSE

       PRINT*,'Solver: ERROR: Unknown global coordinate system: ',TRIM(eq),' Aborting'
       STOP

     END IF

!------------------------------------------------------------------------------
!   Figure out what (flow,heat,stress,...) should be computed, and get
!   memory for the dofs
!------------------------------------------------------------------------------
     DO i=1,NewModel % NumberOfSolvers
       eq = ListGetString( NewModel % Solvers(i) % Values,'Equation' )

       Solver => NewModel % Solvers(i)
       CALL AddEquationBasics( Solver, eq, TransientSimulation )
       CALL AddEquationSolution( Solver, TransientSimulation )
     END DO

!------------------------------------------------------------------------------
!    Add coordinates to list of variables so that coordinate dependent
!    parameter computing routines can ask for them...
!------------------------------------------------------------------------------
     Mesh => NewModel % Meshes 
     NULLIFY( Solver )
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 1',1,Mesh % Nodes % x )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 2',1,Mesh % Nodes % y )
       CALL VariableAdd(Mesh % Variables,Mesh,Solver,'Coordinate 3',1,Mesh % Nodes % z )
       Mesh => Mesh % Next
     END DO

!------------------------------------------------------------------------------
!    Time integration and/or steady state steps
!------------------------------------------------------------------------------
     IF ( TransientSimulation ) THEN
       Timesteps => ListGetIntegerArray( NewModel % Simulation, &
                 'Timestep Intervals', GotIt )

       IF ( .NOT.GotIt ) THEN
         PRINT*,'Solver Input error: Time step intervals must be defined.'
         STOP
       END IF 

       TimestepSizes => ListGetConstRealArray( NewModel % Simulation, &
                             'Timestep Sizes', GotIt )

       IF ( .NOT.GotIt ) THEN
         PRINT*,'Solver Input error: Time step sizes must be defined.'
         STOP
       END IF 

       TimeIntervals = SIZE(Timesteps)

       TimestepMaxIter = ListGetInteger( NewModel % Simulation, &
                  'Steady State Max Iterations' )
!------------------------------------------------------------------------------
     ELSE
!------------------------------------------------------------------------------
!      Steady state
!------------------------------------------------------------------------------
       ALLOCATE( Timesteps(1) )

       Timesteps(1) = ListGetInteger( NewModel % Simulation, &
                  'Steady State Max Iterations' )

       ALLOCATE( TimestepSizes(1,1) )
       TimestepSizes(1,1) = 1.0D0

       TimeIntervals   = 1
       TimeStepMaxIter = 1
     END IF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Get Output File Options
!------------------------------------------------------------------------------
     OutputFile = ListGetString(NewModel % Simulation,'Output File',GotIt)
     IF ( .NOT.GotIt ) OutputFile = 'Result.dat'

     OutputIntervals => ListGetIntegerArray( NewModel % Simulation, &
                     'Output Intervals' )


     CALL SaveToPost

!------------------------------------------------------------------------------
!    Add simulation time to list of system variables so that parameters
!    depending on time  may be set.
!------------------------------------------------------------------------------
     TimeVariable => SimulationTime
     Mesh => NewModel % Meshes
     TimeVariable => SimulationTime
     DO WHILE( ASSOCIATED(Mesh) )
       CALL VariableAdd( Mesh % Variables,Mesh,Solver,'Time',1,TimeVariable )
       Mesh => Mesh % Next
     END DO

! Build QuadrantTree for the old mesh
!    BoundingBox = (/ 0.0d0, 0.0d0, 0.0d0, 0.5d0, 2.11d0, 0.d0 /)
!    Mesh => OldModel % Meshes
!    CALL BuildQuadrantTree(Mesh, BoundingBox, RootQuadrant)
!    QuadrantTreeExists = .TRUE.

     dim = CoordinateSystemDimension()

! Interpolate all variables (that apply to the new mesh)
! from old mesh to new mesh
! If you want to interpolate only a subset of variables,
! pass it as the fourth argument instead of NewModel % Variables

     CALL InterpolateMeshToMesh( OldModel % Meshes, NewModel % Meshes, &
      OldModel % Meshes % Variables, NewModel % Meshes % Variables )

!------------------------------------------------------------------------------
     PRINT*,'*** Interpolation: ALL DONE ***'
!------------------------------------------------------------------------------
! Write the new result file (TimeStep=1, SimulationTime=1.d0)
     PRINT*,'Saving file:',OutputFile
!------------------------------------------------------------------------------
     PRINT*,'*** New result file written ***'
!------------------------------------------------------------------------------
     CALL SaveCurrent
!------------------------------------------------------------------------------
     PRINT*,'*** New post file written ***'
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE SaveCurrent
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Var
     INTEGER :: i
     LOGICAL :: BinaryOutput, SaveAll
     CHARACTER(:), ALLOCATABLE :: Simul, OutputName
 
     Simul = ListGetString( CurrentModel % Simulation, &
                     'Simulation Type' )

     BinaryOutput = ListGetLogical( CurrentModel % Simulation,'Binary Output',GotIt )
     IF ( .NOT.GotIt ) BinaryOutput = .FALSE.

     SaveAll = .NOT.ListGetLogical( CurrentModel % Simulation,&
                                    'Omit unchanged variables in output',GotIt )
     IF ( .NOT.GotIt ) SaveAll = .TRUE.

     Mesh => CurrentModel % Meshes
     DO WHILE( ASSOCIATED( Mesh ) ) 
       IF ( LEN_TRIM(Mesh % Name )>0 ) THEN
         OutputName = TRIM(Mesh % Name) // '/' // TRIM(OutputFile)
       ELSE
         OutputName = OutputFile
       END IF
       IF ( Simul == 'eigen analysis' ) THEN
         DO i=1,CurrentModel % Solvers(1) % NOFEigenValues
           Var => Mesh % Variables
           DO WHILE( ASSOCIATED( Var ) ) 
             IF ( Var % Name  /= 'time' .AND. &
                     Var % Name /= 'coordinate' ) THEN
               Var % Values = REAL(Var % EigenVectors(i,:))
               simulationtime(1) = real(var % eigenvalues(i))
             END IF
             Var => Var % Next
           END DO
           SavedSteps = SaveResult( OutputName,Mesh, &
                   i,SimulationTime(1),BinaryOutput,SaveAll )
         END DO
       ELSE
         SavedSteps = SaveResult( OutputName,Mesh,Time,SimulationTime(1),BinaryOutput,SaveAll )
       END IF
       Mesh => Mesh % Next
     END DO
     CALL SaveToPost
!------------------------------------------------------------------------------
   END SUBROUTINE SaveCurrent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE SaveToPost 
!------------------------------------------------------------------------------
!    Convert results file to post processing file, if requested
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: Var
     INTEGER :: i, TotalTimeSteps = 1
     CHARACTER(:), ALLOCATABLE :: Simul, PostFile, PostName, OutputName
 
     Simul = ListGetString( CurrentModel % Simulation, &
                     'Simulation Type' )

     PostFile = ListGetString( CurrentModel % Simulation,'Post File',GotIt )
     IF ( GotIt ) THEN
       IF ( ParEnv % PEs > 1 ) THEN
         PostFile = PostFile // '.' // i2s(ParEnv % MyPE)
       END IF
       Mesh => CurrentModel % Meshes
       DO WHILE( ASSOCIATED( Mesh ) )
         IF ( LEN_TRIM( Mesh % Name ) > 0 )  THEN
           OutputName = TRIM(Mesh % Name) // '/' // TRIM(OutputFile)
           PostName   = TRIM(Mesh % Name) // '/' // TRIM(PostFile)
         ELSE
           PostName   = PostFile
           OutputName = OutputFile
         END IF
         CALL SetCurrentMesh( CurrentModel, Mesh )
         IF ( Simul == 'eigen analysis' ) THEN
           DO i=1,CurrentModel % Solvers(1) % NOFEigenValues
             Var => Mesh % Variables
             DO WHILE( ASSOCIATED( Var ) ) 
               IF ( Var % Name  /= 'time' .AND. &
                       Var % Name /= 'coordinate' ) THEN
                 Var % Values = REAL(Var % EigenVectors(i,:))
                 simulationtime(1) = real(var % eigenvalues(i))
               END IF
               Var => Var % Next
             END DO
             mesh % savesdone = i
             CALL WritePostFile( PostName,OutputName,CurrentModel, &
              CurrentModel % Solvers(1) % NOFEigenValues, .TRUE. )
           END DO
         ELSE
           CALL WritePostFile( PostName,OutputName,CurrentModel, &
                     TotalTimesteps, .TRUE. )
         END IF
         Mesh => Mesh % Next
       END DO
     END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SaveToPost
!------------------------------------------------------------------------------


   END PROGRAM ResultToResult
!------------------------------------------------------------------------------

!> \}
!> \}
