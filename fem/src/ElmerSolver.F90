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
!
!/******************************************************************************
! *
! *  ELMER/FEM Solver main program
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! *****************************************************************************/

!> \defgroup Solvers Dynamically linked solvers

!> \defgroup UDF Dynamically linked functions

!> \defgroup Programs Utility programs

!> \degroup ElmerLib Elmer library routines

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> The main program for Elmer. Solves the equations as defined by the input files.
!------------------------------------------------------------------------------
   SUBROUTINE ElmerSolver(initialize)
!------------------------------------------------------------------------------

     USE MainUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------

     INTEGER :: Initialize

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     INTEGER :: i,j,k,n,l,t,k1,k2,iter,Ndeg,istat,nproc,tlen,nthreads
     CHARACTER(LEN=MAX_STRING_LEN) :: threads, CoordTransform

     REAL(KIND=dp) :: s,dt,dtfunc
     REAL(KIND=dP), POINTER :: WorkA(:,:,:) => NULL()
     REAL(KIND=dp), POINTER, SAVE :: sTime(:), sStep(:), sInterval(:), sSize(:), &
           steadyIt(:),nonlinIt(:),sPrevSizes(:,:),sPeriodic(:),sPar(:)

     TYPE(Element_t),POINTER :: CurrentElement

     LOGICAL :: GotIt,Transient,Scanning,LastSaved, MeshMode = .FALSE.

     INTEGER :: TimeIntervals,interval,timestep, &
       TotalTimesteps,SavedSteps,CoupledMaxIter,CoupledMinIter

     INTEGER, POINTER, SAVE :: Timesteps(:),OutputIntervals(:), ActiveSolvers(:)
     REAL(KIND=dp), POINTER, SAVE :: TimestepSizes(:,:)

     INTEGER(KIND=AddrInt) :: ControlProcedure

     LOGICAL :: InitDirichlet, ExecThis

     TYPE(ElementType_t),POINTER :: elmt

     TYPE(ParEnv_t), POINTER :: ParallelEnv

     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName, eq, ExecCommand, ExtrudedMeshName
     CHARACTER(LEN=MAX_STRING_LEN) :: OutputFile, PostFile, RestartFile, &
                OutputName=' ',PostName=' ', When, OptionString

     TYPE(Variable_t), POINTER :: Var
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: Solver

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: CT0,RT0,tt
#else
     REAL(KIND=dp) :: RealTime,CPUTime,CT0,RT0,tt
#endif

     LOGICAL :: FirstLoad = .TRUE., FirstTime=.TRUE., Found
     LOGICAL :: Silent, Version, GotModelName

     INTEGER :: iargc, NoArgs

     INTEGER :: ExtrudeLevels, MeshIndex
     TYPE(Mesh_t), POINTER :: ExtrudedMesh

     INTEGER :: omp_get_max_threads

#ifdef HAVE_TRILINOS
INTERFACE
      SUBROUTINE TrilinosCleanup() BIND(C,name='TrilinosCleanup')
      IMPLICIT NONE
      END SUBROUTINE TrilinosCleanup
END INTERFACE
#endif

     ! Start the watches, store later
     !--------------------------------
     RT0 = RealTime()
     CT0 = CPUTime()


     ! If parallel execution requested, initialize parallel environment:
     !------------------------------------------------------------------
     IF(FirstTime)  ParallelEnv => ParallelInit()

     OutputPE = -1
     IF( ParEnv % MyPe == 0 ) THEN
       OutputPE = 0 
     END IF
     
     IF ( FirstTime ) THEN

       !
       ! Print banner to output:
#include "../config.h"
       ! -----------------------
#ifdef USE_ISO_C_BINDINGS
       NoArgs = COMMAND_ARGUMENT_COUNT()
#else
       NoArgs = IARGC()
#endif
       ! Info Level is always true until the model has been read!
       ! This makes it possible to cast something 
       Silent = .FALSE.
       Version = .FALSE.
       IF( NoArgs > 0 ) THEN 
         DO i = 1, NoArgs 
#ifdef USE_ISO_C_BINDINGS
           CALL GET_COMMAND_ARGUMENT(i, OptionString)
#else
           CALL getarg( i,OptionString )
#endif
           Silent = Silent .OR. &
               ( OptionString=='-s' .OR. OptionString=='--silent' ) 
           Version = Version .OR. &
               ( OptionString=='-v' .OR. OptionString == '--version' )
         END DO
       END IF

       ! Set number of OpenMP threads
       nthreads = 1
       !$ nthreads = omp_get_max_threads()
       IF (nthreads > 1) THEN
         ! Check if OMP_NUM_THREADS environment variable is set
#if USE_ISO_C_BINDINGS
         CALL envir( 'OMP_NUM_THREADS', threads, tlen )
#else
         CALL envir( 'OMP_NUM_THREADS'//CHAR(0), threads, tlen )
#endif
         IF (tlen==0) THEN
           CALL Warn('MAIN','OMP_NUM_THREADS not set. Using only 1 thread.')
           nthreads = 1
           ! Set number of threads to 1
           !$ CALL omp_set_num_threads(nthreads)
#ifdef HAVE_MKL
           CALL mkl_set_num_threads(nthreads)
#endif
         END IF
       END IF


       IF( .NOT. Silent ) THEN
         CALL Info( 'MAIN', ' ')
         CALL Info( 'MAIN', '=============================================================')
         CALL Info( 'MAIN', 'ElmerSolver finite element software, Welcome!                ')
         CALL Info( 'MAIN', 'This program is free software licensed under (L)GPL          ')
         CALL Info( 'MAIN', 'Copyright 1st April 1995 - , CSC - IT Center for Science Ltd.')
         CALL Info( 'MAIN', 'Webpage http://www.csc.fi/elmer, Email elmeradm@csc.fi       ')
         CALL Info( 'MAIN', 'Version: ' // VERSION &
#ifdef REVISION
             // ' (Rev: ' // REVISION  &
#endif
#ifdef COMPILATIONDATE
             // ', Compiled: ' // COMPILATIONDATE // ')' &
#else
             // ')' &
#endif
         )

         IF ( ParEnv % PEs > 1 ) &
             CALL Info( 'MAIN', ' Running in parallel using ' // &
             TRIM(i2s(ParEnv % PEs)) // ' tasks.')

         ! Print out number of threads in use
         IF ( nthreads > 1 ) &
             CALL Info('MAIN', ' Running in parallel with ' // &
                       TRIM(i2s(nthreads)) // ' threads per task.')

#ifdef HAVE_HYPRE
         CALL Info( 'MAIN', ' HYPRE library linked in.')
#endif
#ifdef HAVE_TRILINOS
         CALL Info( 'MAIN', ' Trilinos library linked in.')
#endif
#ifdef HAVE_MUMPS
         CALL Info( 'MAIN', ' MUMPS library linked in.')
#endif
#ifdef HAVE_CHOLMOD
         CALL Info( 'MAIN', ' CHOLMOD library linked in.')
#endif
#ifdef HAVE_SUPERLU
         CALL Info( 'MAIN', ' SUPERLU library linked in.')
#endif
#ifdef HAVE_PARDISO
         CALL Info( 'MAIN', ' PARDISO library linked in.')
#endif
#ifdef HAVE_MKL
         CALL Info( 'MAIN', ' Intel MKL linked in.' )
#endif
         CALL Info( 'MAIN', '=============================================================')
       END IF

       IF( Version ) RETURN
       
       CALL InitializeElementDescriptions()
     END IF

     ! Read input file name either as an argument, or from the default file:
     !----------------------------------------------------------------------
     GotModelName = .FALSE.
     IF ( ParEnv % PEs <= 1 .AND. NoArgs > 0 ) THEN
#ifdef USE_ISO_C_BINDINGS
       CALL GET_COMMAND_ARGUMENT(1, ModelName)
#else
       CALL getarg( 1,ModelName )
#endif
       IF( ModelName(1:1) /= '-') THEN 
         GotModelName = .TRUE.
#ifdef USE_ISO_C_BINDINGS
         IF (NoArgs > 1) CALL GET_COMMAND_ARGUMENT(2, eq)
#else
         IF ( NoArgs > 1 ) CALL getarg( 2,eq )
#endif
       END IF
     END IF
         
     IF( .NOT. GotModelName ) THEN
       OPEN( 1, File='ELMERSOLVER_STARTINFO', STATUS='OLD', ERR=10 )
       READ(1,'(a)') ModelName
       CLOSE(1)
     END IF

!------------------------------------------------------------------------------
!    Read element definition file, and initialize element types
!------------------------------------------------------------------------------
     IF ( Initialize==1 ) THEN
       CALL FreeModel(CurrentModel)
       FirstLoad=.TRUE.
     END IF

!------------------------------------------------------------------------------
!    Read Model and mesh from Elmer mesh data base
!------------------------------------------------------------------------------
     MeshIndex = 0
     DO WHILE( .TRUE. )

       IF ( initialize==2 ) GOTO 1

       IF(MeshMode) THEN
         CALL FreeModel(CurrentModel)
         MeshIndex = MeshIndex + 1
         FirstLoad=.TRUE.
       END IF

       IF ( FirstLoad ) THEN
         IF( .NOT. Silent ) THEN
           CALL Info( 'MAIN', ' ')
           CALL Info( 'MAIN', ' ')
           CALL Info( 'MAIN', '-------------------------------------')
           CALL Info( 'MAIN', 'Reading Model: '//TRIM( ModelName) )
         END IF

         INQUIRE(Unit=InFileUnit, Opened=GotIt)
         IF ( gotIt ) CLOSE(inFileUnit)
         
         OPEN( Unit=InFileUnit, Action='Read',File=ModelName,Status='OLD',ERR=20 )
         CurrentModel => LoadModel(ModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE,MeshIndex)
         IF(.NOT.ASSOCIATED(CurrentModel)) EXIT

         MeshMode = ListGetLogical( CurrentModel % Simulation, 'Mesh Mode', Found)

         !------------------------------------------------------------------------------
         ! Some keywords automatically require other keywords to be set
         ! We could complain on the missing keywords later on, but sometimes 
         ! it may be just as simple to add them directly. 
         !------------------------------------------------------------------------------
         CALL CompleteModelKeywords( )


         ! Optionally perform simple extrusion to increase the dimension of the mesh
         !----------------------------------------------------------------------------------
         ExtrudeLevels=GetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found)
         IF (Found) THEN
            IF(ExtrudeLevels>1) THEN
               ExtrudedMeshName = GetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
               IF (Found) THEN
                  ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, ExtrudeLevels-2, ExtrudedMeshName)
               ELSE
                  ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, ExtrudeLevels-2)
               END IF
               DO i=1,CurrentModel % NumberOfSolvers
                  IF(ASSOCIATED(CurrentModel % Solvers(i) % Mesh,CurrentModel % Meshes)) &
                       CurrentModel % Solvers(i) % Mesh => ExtrudedMesh 
               END DO
               ExtrudedMesh % Next => CurrentModel % Meshes % Next
               CurrentModel % Meshes => ExtrudedMesh
               
               ! If periodic BC given, compute boundary mesh projector:
               ! ------------------------------------------------------
               DO i = 1,CurrentModel % NumberOfBCs
                  IF(ASSOCIATED(CurrentModel % Bcs(i) % PMatrix)) &
                       CALL FreeMatrix( CurrentModel % BCs(i) % PMatrix )
                  CurrentModel % BCs(i) % PMatrix => NULL()
                  k = ListGetInteger( CurrentModel % BCs(i) % Values, 'Periodic BC', GotIt )
                  IF( GotIt ) THEN
                     CurrentModel % BCs(i) % PMatrix =>  PeriodicProjector( CurrentModel, ExtrudedMesh, i, k )
                  END IF
               END DO
            END IF
         END IF

         ! If requested perform coordinate transformation directly after is has been obtained.
         ! Don't maintain the original mesh. 
         !----------------------------------------------------------------------------------
         CoordTransform = ListGetString(CurrentModel % Simulation,'Coordinate Transformation',GotIt)
         IF( GotIt ) THEN
           CALL CoordinateTransformation( CurrentModel % Meshes, CoordTransform, &
               CurrentModel % Simulation, .TRUE. )
         END IF

         IF(.NOT. Silent ) THEN
            CALL Info( 'MAIN', '-------------------------------------')
          END IF 
       ELSE
         IF ( Initialize==3 ) THEN
           INQUIRE(Unit=InFileUnit, Opened=GotIt)
           IF ( gotIt ) CLOSE(inFileUnit)
           OPEN( Unit=InFileUnit, Action='Read', & 
                    File=ModelName,Status='OLD',ERR=20 )
         END IF

         IF ( .NOT.ReloadInputFile(CurrentModel) ) EXIT

         Mesh => CurrentModel % Meshes
         DO WHILE( ASSOCIATED(Mesh) )
           Mesh % SavesDone = 0
           Mesh => Mesh % Next
         END DO
       END IF

1      CONTINUE

       CALL ListAddLogical( CurrentModel % Simulation, &
             'Initialization Phase', .TRUE. )

       ! Save the start time to the list so that it may be retrieved when necessary
       ! This could perhaps also be a global variable etc, but this will do for now.
       !-------------------------------------------------------------------------
       IF( ListGetLogical( CurrentModel % Simulation,'Simulation Timing',GotIt) ) THEN
         CALL ListAddConstReal( CurrentModel % Simulation,'cputime0',ct0 )
         CALL ListAddConstReal( CurrentModel % Simulation,'realtime0',rt0 )
       END IF

!------------------------------------------------------------------------------
!      Check for transient case
!------------------------------------------------------------------------------
       eq = ListGetString( CurrentModel % Simulation, 'Simulation Type', GotIt )
       Scanning  = eq == 'scanning'
       Transient = eq == 'transient'

!------------------------------------------------------------------------------
!      To more conveniently support the use of VTK based visualization there 
!      is a hack that recognizes VTU suffix and creates a instance of output solver.
!      Note that this is really quite a dirty hack, and is not a good example.
!-----------------------------------------------------------------------------
       CALL AddVtuOutputSolverHack()

!------------------------------------------------------------------------------
!      Figure out what (flow,heat,stress,...) should be computed, and get
!      memory for the dofs
!------------------------------------------------------------------------------
       CALL AddSolvers()

!------------------------------------------------------------------------------
!      Time integration and/or steady state steps
!------------------------------------------------------------------------------
       IF ( Transient .OR. Scanning ) THEN
         Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
                       'Timestep Intervals', GotIt )

         IF ( .NOT.GotIt ) THEN
           CALL Fatal( ' ', 'Keyword > Timestep Intervals < MUST be ' //  &
                   'defined for transient and scanning simulations' )
         END IF 
         TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation, &
                               'Timestep Sizes', GotIt )

         IF ( .NOT.GotIt ) THEN
           IF( Scanning .OR. ListCheckPresent( CurrentModel % Simulation,'Timestep Size') ) THEN
             ALLOCATE(TimestepSizes(SIZE(Timesteps),1))
             TimestepSizes = 1.0_dp
           ELSE
             CALL Fatal( ' ', 'Keyword [Timestep Sizes] MUST be ' //  &
                 'defined for time dependent simulations' )
             STOP
           END IF
         END IF 

         TimeIntervals = SIZE(Timesteps)

         CoupledMaxIter = ListGetInteger( CurrentModel % Simulation, &
               'Steady State Max Iterations', GotIt, minv=1 )
         IF ( .NOT. GotIt ) CoupledMaxIter = 1
!------------------------------------------------------------------------------
       ELSE
!------------------------------------------------------------------------------
!        Steady state
!------------------------------------------------------------------------------
         ALLOCATE(Timesteps(1))

         Timesteps(1) = ListGetInteger( CurrentModel % Simulation, &
               'Steady State Max Iterations', GotIt,minv=1 )
         IF ( .NOT. GotIt ) Timesteps(1)=1
  
         ALLOCATE(TimestepSizes(1,1))
         TimestepSizes(1,1) = 1.0D0

         CoupledMaxIter = 1
         TimeIntervals  = 1
       END IF

       IF ( FirstLoad ) &
         ALLOCATE( sTime(1), sStep(1), sInterval(1), sSize(1), &
             steadyIt(1), nonLinit(1), sPrevSizes(1,5), sPeriodic(1), sPar(1) )

       dt   = 0._dp

       sTime = 0._dp
       sStep = 0
       sPeriodic = 0._dp

       sSize = dt
       sPrevSizes = 0_dp

       sInterval = 0._dp

       steadyIt = 0
       nonlinIt = 0
       sPar = 0

       CoupledMinIter = ListGetInteger( CurrentModel % Simulation, &
                  'Steady State Min Iterations', GotIt )

!------------------------------------------------------------------------------
!      Add coordinates and simulation time to list of variables so that
!      coordinate dependent parameter computing routines can ask for
!      them...
!------------------------------------------------------------------------------
       IF ( FirstLoad ) CALL AddMeshCoordinatesAndTime()

!------------------------------------------------------------------------------
!      Get Output File Options
!------------------------------------------------------------------------------

       OutputIntervals => ListGetIntegerArray( CurrentModel % Simulation, &
                       'Output Intervals', GotIt )
       IF ( .NOT. GotIt ) THEN
         ALLOCATE( OutputIntervals(SIZE(TimeSteps)) )
         OutputIntervals = 1
       END IF


       ! Initial Conditions:
       ! -------------------
       IF ( FirstLoad ) CALL SetInitialConditions()


       ! Compute the total number of steps that will be saved to the files
       ! Particularly look if the last step will be saved, or if it has
       ! to be saved separately.
       !------------------------------------------------------------------
       TotalTimesteps = 0
       LastSaved = .TRUE.
       DO interval=1,TimeIntervals
         DO timestep = 1,Timesteps(interval)
           IF( OutputIntervals(Interval) == 0 ) CYCLE
           LastSaved = .FALSE.
           IF ( MOD(Timestep-1, OutputIntervals(Interval))==0 ) THEN
              LastSaved = .TRUE.
              TotalTimesteps = TotalTimesteps + 1
           END IF
         END DO
       END DO

       DO i=1,CurrentModel % NumberOfSolvers
          Solver => CurrentModel % Solvers(i)
          IF(.NOT.ASSOCIATED(Solver % Variable)) CYCLE
          IF(.NOT.ASSOCIATED(Solver % Variable % Values)) CYCLE
          When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
          IF ( GotIt ) THEN
             IF ( When == 'after simulation' .OR. When == 'after all' ) THEN
                LastSaved = .FALSE.
             END IF
          ELSE
           IF ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_ALL ) THEN
              LastSaved = .FALSE.
           END IF
          END IF
       END DO

       IF ( .NOT.LastSaved ) TotalTimesteps = TotalTimesteps + 1
       IF( TotalTimesteps == 0 ) TotalTimesteps = 1

       CALL ListAddLogical( CurrentModel % Simulation,  &
            'Initialization Phase', .FALSE. )

       FirstLoad = .FALSE.
       IF ( Initialize == 1 ) EXIT

!------------------------------------------------------------------------------
!      Here we actually start the simulation ....
!      First go trough timeintervals
!------------------------------------------------------------------------------
       ExecCommand = ListGetString( CurrentModel % Simulation, &
                 'Control Procedure', GotIt )
       IF ( GotIt ) THEN
          ControlProcedure = GetProcAddr( ExecCommand )
          CALL ExecSimulationProc( ControlProcedure, CurrentModel )
       ELSE
          CALL ExecSimulation( TimeIntervals, CoupledMinIter, &
              CoupledMaxIter, OutputIntervals, Transient, Scanning)
       END IF
!------------------------------------------------------------------------------
!    Always save the last step to output
!------------------------------------------------------------------------------
       IF ( .NOT.LastSaved ) THEN
         DO i=1,CurrentModel % NumberOfSolvers
            Solver => CurrentModel % Solvers(i)
            IF ( Solver % PROCEDURE == 0 ) CYCLE
            ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
            When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
            IF ( GotIt ) ExecThis = ( When == 'before saving') 
            IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
         END DO

         CALL SaveToPost(0)
         CALL SaveCurrent(Timestep)

         DO i=1,CurrentModel % NumberOfSolvers
            Solver => CurrentModel % Solvers(i)
            IF ( Solver % PROCEDURE == 0 ) CYCLE
            ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_SAVE)
            When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
            IF ( GotIt ) ExecThis = ( When == 'after saving') 
            IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
         END DO
       END IF

       CALL CompareToReferenceSolution( )

       IF ( Initialize >= 2 ) EXIT
     END DO


     CALL CompareToReferenceSolution( Finalize = .TRUE. )


!------------------------------------------------------------------------------
!    THIS IS THE END (...,at last, the end, my friend,...)
!------------------------------------------------------------------------------
     IF ( Initialize /= 1 ) CALL Info( 'ElmerSolver', '*** Elmer Solver: ALL DONE ***',Level=3 )

     IF ( Initialize <= 0 ) CALL FreeModel(CurrentModel)

#ifdef HAVE_TRILINOS
  CALL TrilinosCleanup()
#endif

     IF ( FirstTime ) CALL ParallelFinalize()
     FirstTime = .FALSE.
     CALL Info('ElmerSolver','The end',Level=3)

     RETURN

10   CONTINUE
     CALL Fatal( 'ElmerSolver', 'Unable to find ELMERSOLVER_STARTINFO, can not execute.' )

20   CONTINUE
     CALL Fatal( 'ElmerSolver', 'Unable to find input file [' // &
              TRIM(Modelname) // '], can not execute.' )

   CONTAINS 
     

     ! The user may request unit tests to be performed. 
     ! This will be done if any reference norm is given.
     ! The success will be written to file TEST.PASSED as 0/1. 
     !--------------------------------------------------------
     SUBROUTINE CompareToReferenceSolution( Finalize ) 
       LOGICAL, OPTIONAL :: Finalize 

       INTEGER :: i, j, k, n, solver_id, TestCount=0, PassCount=0, FailCount, Dofs
       REAL(KIND=dp) :: Norm, RefNorm, Tol, Err, val, refval
       TYPE(Solver_t), POINTER :: Solver
       TYPE(Variable_t), POINTER :: Var
       LOGICAL :: Found, Success = .TRUE., FinalizeOnly, CompareNorm, CompareSolution, AbsoluteErr
       CHARACTER(LEN=MAX_STRING_LEN) :: PassedMsg

       SAVE TestCount, PassCount 

       IF( ParEnv % MyPe > 0 ) RETURN

       ! Write the success to a file for further use e.g. by cmake
       !----------------------------------------------------------
       FinalizeOnly = .FALSE.
       IF( PRESENT( Finalize ) ) FinalizeOnly = .TRUE.

       IF( FinalizeOnly ) THEN

         ! Nothing tested
         IF( TestCount == 0 ) RETURN
         
         Success = ( PassCount == TestCount )  
         FailCount = TestCount - PassCount
         
         IF( Success ) THEN
           CALL Info('CompareToReferenceSolution',&
               'PASSED all '//TRIM(I2S(TestCount))//' tests!',Level=4)
         ELSE         
           CALL Warn('CompareToReferenceSolution','FAILED '//TRIM(I2S(FailCount))//&
               ' tests out of '//TRIM(I2S(TestCount))//'!')
         END IF
         
         IF( FinalizeOnly ) THEN
           IF( ParEnv % MyPe == 0 ) THEN
             IF( ParEnv % PEs > 1 ) THEN
               ! Parallel test, add the number of tasks as a suffix
               WRITE(PassedMsg, '("TEST.PASSED_",I0)') ParEnv % PEs
               OPEN( 10, FILE = PassedMsg )
             ELSE
               OPEN( 10, FILE = 'TEST.PASSED' )
             END IF
             IF( Success ) THEN
               WRITE( 10,'(I1)' ) 1
             ELSE
               WRITE( 10,'(I1)' ) 0
             END IF
             CLOSE( 10 )
           END IF
         END IF

         RETURN
       END IF


       DO solver_id=1,CurrentModel % NumberOfSolvers
         Solver => CurrentModel % Solvers(solver_id)

         RefNorm = ListGetConstReal( Solver % Values,'Reference Norm', CompareNorm )
         CompareSolution = ListCheckPrefix( Solver % Values,'Reference Solution')
         
         IF(.NOT. ( CompareNorm .OR. CompareSolution ) ) CYCLE
         
         Var => Solver % Variable
         IF( .NOT. ASSOCIATED( Var ) ) THEN
           CALL Warn('CompareToReferenceSolution','Variable in Solver '&
               //TRIM(I2S(i))//' not associated, cannot compare')
           CYCLE
         END IF

         TestCount = TestCount + 1
         Success = .TRUE.

         ! Compare either to existing norm (ensures consistancy) 
         ! or to existing solution (may also be used to directly verify)
         ! Usually only either of these is given but for the sake of completeness
         ! both may be used at the same time. 
         IF( CompareNorm ) THEN
           Tol = ListGetConstReal( Solver % Values,'Reference Norm Tolerance', Found )
           IF(.NOT. Found ) Tol = 1.0d-5
           Norm = Var % Norm 
           AbsoluteErr = ListGetLogical( Solver % Values,'Reference Norm Absolute', Found ) 
           Err = ABS( Norm - RefNorm ) 

           IF(.NOT. AbsoluteErr ) THEN
             IF( RefNorm < TINY( RefNorm ) ) THEN
               CALL Warn('CompareToReferenceSolution','Refenrece norm too small for relative error')
               AbsoluteErr = .TRUE.
             ELSE
               Err = Err / RefNorm 
             END IF
           END IF

           ! Compare to given reference norm
           IF( Err > Tol ) THEN
             ! Warn only in the main core
             IF( ParEnv % MyPe == 0 ) THEN
               WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                   'Solver ',solver_id,' FAILED:  Norm =',Norm,'  RefNorm =',RefNorm
               CALL Warn('CompareToReferenceSolution',Message)
               IF( AbsoluteErr ) THEN
                 WRITE( Message,'(A,ES13.6)') 'Absolute Error to reference norm:',Err
               ELSE
                 WRITE( Message,'(A,ES13.6)') 'Relative Error to reference norm:',Err
               END IF
               CALL Info('CompareToReferenceSolution',Message, Level = 4 )
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                 'Solver ',solver_id,' PASSED:  Norm =',Norm,'  RefNorm =',RefNorm
             CALL Info('CompareToReferenceSolution',Message,Level=4)
           END IF
         END IF

         IF( CompareSolution ) THEN
           Tol = ListGetConstReal( Solver % Values,'Reference Solution Tolerance', Found )
           IF(.NOT. Found ) Tol = 1.0d-5
           Dofs = Var % Dofs
           n = 0 
           RefNorm = 0.0_dp
           Norm = 0.0_dp
           Err = 0.0_dp
           DO i=1,Solver % Mesh % NumberOfNodes
             j = Var % Perm(i)
             IF( j == 0 ) CYCLE
             DO k=1,Dofs
               IF( Dofs == 1 ) THEN
                 refval = ListGetRealAtNode( Solver % Values,'Reference Solution',i,Found ) 
               ELSE
                 refval = ListGetRealAtNode( Solver % Values,'Reference Solution '//TRIM(I2S(k)),i,Found ) 
               END IF
               IF( Found ) THEN
                 val = Var % Values( Dofs*(j-1)+k)
                 RefNorm = RefNorm + refval**2
                 Norm = Norm + val**2
                 Err = Err + (refval-val)**2
                 n = n + 1
               END IF
             END DO
           END DO
           IF( ParEnv % PEs > 1 ) CALL Warn('CompareToReferefenSolution','Not implemented in parallel!')
           IF( n == 0 ) CALL Fatal('CompareToReferenceSolution','Could not find any reference solution')
           RefNorm = SQRT( RefNorm / n ) 
           Norm = SQRT( Norm / n )
           Err = SQRT( Err / n ) 

           IF( Err > Tol ) THEN
             ! Normally warning is done for every partition but this time it is the same for all
             IF( ParEnv % MyPe == 0 ) THEN
               WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                   'Solver ',solver_id,' FAILED:  Solution = ',Norm,'  RefSolution =',RefNorm
               CALL Warn('CompareToReferenceSolution',Message)
               WRITE( Message,'(A,ES13.6)') 'Relative Error to reference solution:',Err
               CALL Info('CompareToReferenceSolution',Message, Level = 4 )
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                 'Solver ',solver_id,' PASSED:  Solution =',Norm,'  RefSolution =',RefNorm
             CALL Info('CompareToReferenceSolution',Message,Level=4)
           END IF
         END IF

         IF( Success ) PassCount = PassCount + 1
       END DO
      

     END SUBROUTINE CompareToReferenceSolution

     


     ! This is a dirty hack that adds an instance of ResultOutputSolver to the list of Solvers.
     ! The idea is that it is much easier for the end user to take into use the vtu output this way.
     ! The solver itself has limited set of parameters needed and is therefore approapriate for this
     ! kind of hack. It can of course be also added as a regular solver also.
     !----------------------------------------------------------------------------------------------
     SUBROUTINE AddVtuOutputSolverHack()     
       TYPE(Solver_t), POINTER :: ABC(:), PSolver
       CHARACTER(LEN=MAX_NAME_LEN) :: str
       INTEGER :: i,j,j2,j3,k,n
       TYPE(ValueList_t), POINTER :: Params
       LOGICAL :: gotIt, VtuFormat

       str = ListGetString( CurrentModel % Simulation,'Post File',GotIt) 
       IF(.NOT. GotIt) RETURN
       k = INDEX( str,'.vtu' )
       VtuFormat = ( k /= 0 ) 

       IF(.NOT. VtuFormat ) RETURN
       
       CALL Info('AddVtuOutputSolverHack','Adding ResultOutputSolver to write VTU output in file: '&
           //TRIM(str(1:k-1)))
     
       CALL ListRemove( CurrentModel % Simulation,'Post File')
       n = CurrentModel % NumberOfSolvers+1
       ALLOCATE( ABC(n) )
       DO i=1,n-1
         ! Def_Dofs is the only allocatable structure within Solver_t:
         IF( ALLOCATED( CurrentModel % Solvers(i) % Def_Dofs ) ) THEN
           j = SIZE(CurrentModel % Solvers(i) % Def_Dofs,1)
           j2 = SIZE(CurrentModel % Solvers(i) % Def_Dofs,2)
           j3 = SIZE(CurrentModel % Solvers(i) % Def_Dofs,3)
           ALLOCATE( ABC(i) % Def_Dofs(j,j2,j3) )
         END IF

         ! Copy the content of the Solver structure
         ABC(i) = CurrentModel % Solvers(i)

         ! Nullify the old structure since otherwise bad things may happen at deallocation
         NULLIFY( CurrentModel % Solvers(i) % ActiveElements )
         NULLIFY( CurrentModel % Solvers(i) % Mesh )
         NULLIFY( CurrentModel % Solvers(i) % BlockMatrix )
         NULLIFY( CurrentModel % Solvers(i) % Matrix )
         NULLIFY( CurrentModel % Solvers(i) % Variable )
       END DO

       ! Deallocate the old structure and set the pointer to the new one
       DEALLOCATE( CurrentModel % Solvers )
       CurrentModel % Solvers => ABC
       CurrentModel % NumberOfSolvers = n

       ! Now create the ResultOutputSolver instance on-the-fly
       CurrentModel % Solvers(n) % PROCEDURE = 0
       NULLIFY( CurrentModel % Solvers(n) % Matrix )
       NULLIFY( CurrentModel % Solvers(n) % BlockMatrix )
       NULLIFY( CurrentModel % Solvers(n) % Variable )
       NULLIFY( CurrentModel % Solvers(n) % ActiveElements )
       CurrentModel % Solvers(n) % NumberOfActiveElements = 0
       j = CurrentModel % NumberOfBodies
       ALLOCATE( CurrentModel % Solvers(n) % Def_Dofs(10,j,6))
       CurrentModel % Solvers(n) % Def_Dofs(:,1:j,6) = -1
       
       ! Add some keywords to the list
       CurrentModel % Solvers(n) % Values => ListAllocate()
       CALL ListAddString(CurrentModel % Solvers(n) % Values,&
           'Procedure', 'ResultOutputSolve ResultOutputSolver',.FALSE.)
       CALL ListAddString(CurrentModel % Solvers(n) % Values,'Output Format','vtu')
       CALL ListAddString(CurrentModel % Solvers(n) % Values,'Output File Name',str(1:k-1))
       CALL ListAddString(CurrentModel % Solvers(n) % Values,'Exec Solver','after saving')
       CALL ListAddString(CurrentModel % Solvers(n) % Values,'Equation','InternalVtuOutputSolver')
       CALL ListAddLogical(CurrentModel % Solvers(n) % Values,'Save Geometry IDs',.TRUE.)

     END SUBROUTINE AddVtuOutputSolverHack


!------------------------------------------------------------------------------
!> Adds flags for active solvers. 
!------------------------------------------------------------------------------
    SUBROUTINE AddSolvers()
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,nlen
      LOGICAL :: InitSolver, Found
!------------------------------------------------------------------------------

      CALL Info('AddSolvers','Setting up '//TRIM(I2S(CurrentModel % NumberOfSolvers))//&
          ' solvers',Level=10)

      ! This is a hack that sets Equation flags True for the "Active Solvers".
      ! The Equation flag is the legacy way of setting a Solver active and is still
      ! used internally.
      !----------------------------------------------------------------------------
      DO i=1,CurrentModel % NumberOfSolvers

        eq = ListGetString( CurrentModel % Solvers(i) % Values,'Equation', Found )
     
        IF ( Found ) THEN
          nlen = LEN_TRIM(eq)
          DO j=1,CurrentModel % NumberOFEquations
             ActiveSolvers => ListGetIntegerArray( CurrentModel % Equations(j) % Values, &
                                'Active Solvers', Found )
             IF ( Found ) THEN
                DO k=1,SIZE(ActiveSolvers)
                   IF ( ActiveSolvers(k) == i ) THEN
                      CALL ListAddLogical( CurrentModel % Equations(j) % Values, eq(1:nlen), .TRUE. )
                      EXIT
                   END IF
                END DO
             END IF
          END DO
       END IF
     END DO

     ! Add the dynamically linked solver to be called later
     !---------------------------------------------------------------------
     DO i=1,CurrentModel % NumberOfSolvers
        eq = ListGetString( CurrentModel % Solvers(i) % Values,'Equation', Found )
        CALL Info('AddSolvers','Setting up solver '//TRIM(I2S(i))//': '//TRIM(eq),Level=10)

        Solver => CurrentModel % Solvers(i)
        InitSolver = ListGetLogical( Solver % Values, 'Initialize', Found )
        IF ( Found .AND. InitSolver ) THEN
          CALL FreeMatrix( Solver % Matrix )
          CALL ListAddLogical( Solver % Values, 'Initialize', .FALSE. )
        END IF

        IF ( Solver % PROCEDURE == 0 .OR. InitSolver ) THEN
          IF ( .NOT. ASSOCIATED( Solver % Mesh ) ) THEN
            Solver % Mesh => CurrentModel % Meshes
          END IF
          CurrentModel % Solver => Solver
          CALL AddEquationBasics( Solver, eq, Transient )
          CALL AddEquationSolution( Solver, Transient )
        END IF
     END DO

     CALL Info('AddSolvers','Setting up solvers done',Level=12)

!------------------------------------------------------------------------------
  END SUBROUTINE AddSolvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds coordinate and time variables to the current mesh structure. 
!------------------------------------------------------------------------------
  SUBROUTINE AddMeshCoordinatesAndTime()
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: DtVar

     CALL Info('AddMeshCoordinatesAndTime','Setting mesh coordinates and time',Level=10)

     NULLIFY( Solver )

     Mesh => CurrentModel % Meshes 
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd( Mesh % Variables, Mesh,Solver, &
             'Coordinate 1',1,Mesh % Nodes % x )

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, &
             'Coordinate 2',1,Mesh % Nodes % y )

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, &
             'Coordinate 3',1,Mesh % Nodes % z )

       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Time', 1, sTime )
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Periodic Time', 1, sPeriodic )
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Timestep', 1, sStep )
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Timestep size', 1, sSize )
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Timestep interval', 1, sInterval )

       ! Save some previous timesteps for variable timestep multistep methods
       DtVar => VariableGet( Mesh % Variables, 'Timestep size' )
       DtVar % PrevValues => sPrevSizes

       CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
               'nonlin iter', 1, nonlinIt )
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, &
               'coupled iter', 1, steadyIt )

       sPar(1) = 1.0_dp * ParEnv % MyPe 
       CALL VariableAdd( Mesh % Variables, Mesh, Solver, 'Partition', 1, sPar ) 

       Mesh => Mesh % Next
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddMeshCoordinatesAndTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sets initial conditions for the fields. 
!------------------------------------------------------------------------------
   SUBROUTINE SetInitialConditions()
!------------------------------------------------------------------------------
     USE DefUtils
     INTEGER :: DOFs
     CHARACTER(LEN=MAX_NAME_LEN) :: str
     LOGICAL :: Found, NamespaceFound
     TYPE(Solver_t), POINTER :: Solver
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp),ALLOCATABLE :: Work(:)

     INTEGER :: i,j,k,l,m,vect_dof,real_dof,dim

     REAL(KIND=dp) :: nrm(3),t1(3),t2(3),vec(3),tmp(3),udot
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Nodes_t), SAVE :: Nodes
     LOGICAL :: nt_boundary
     TYPE(Element_t), POINTER :: Element
     TYPE(Variable_t), POINTER :: var, vect_var

     CALL Info('SetInitialConditions','Setting up initial conditions (if any)',Level=10)


     dim = CoordinateSystemDimension()

     IF (GetLogical(GetSimulation(),'Restart Before Initial Conditions',Found)) THEN
       CALL Restart
       CALL InitCond
     ELSE
       CALL InitCond
       CALL Restart
     END IF

!------------------------------------------------------------------------------
!    Make sure that initial values at boundaries are set correctly.
!    NOTE: This overrides the initial condition setting for field variables!!!!
!-------------------------------------------------------------------------------
     InitDirichlet = ListGetLogical( CurrentModel % Simulation, &
            'Initialize Dirichlet Conditions', GotIt ) 
     IF ( .NOT. GotIt ) InitDirichlet = .TRUE.

     vect_var => NULL()
     IF ( InitDirichlet ) THEN
       Mesh => CurrentModel % Meshes
       DO WHILE( ASSOCIATED(Mesh) )
         ALLOCATE( Work(Mesh % MaxElementDOFs) )
         CALL SetCurrentMesh( CurrentModel, Mesh )

         DO t = Mesh % NumberOfBulkElements + 1, &
                 Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           Element => Mesh % Elements(t)

           ! Set also the current element pointer in the model structure to
           ! reflect the element being processed:
           ! ---------------------------------------------------------------
           CurrentModel % CurrentElement => Element
           n = Element % TYPE % NumberOfNodes

           BC => GetBC()
           IF(.NOT.ASSOCIATED(BC)) CYCLE

           Var => Mesh % Variables
           DO WHILE( ASSOCIATED(Var) )
             Solver => Var % Solver
             IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver

             str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
             IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))


             IF ( Var % DOFs <= 1 ) THEN
               Work(1:n) = GetReal( BC,Var % Name, gotIt )
               IF ( GotIt ) THEN

                 nt_boundary = .FALSE.
                 IF ( GetElementFamily() /= 1 ) THEN
                   k = LEN_TRIM(var % name)
                   vect_dof = ICHAR(Var % Name(k:k))-ICHAR('0');
                   IF ( vect_dof>=1 .AND. vect_dof<= 3 ) THEN
                     nt_boundary =  GetLogical( BC, &
                        'normal-tangential '//var % name(1:k-2), gotIt)

                     IF ( nt_boundary ) THEN
                       nt_boundary = .FALSE.
                       Vect_var => Mesh % Variables
                       DO WHILE( ASSOCIATED(Vect_var) )
                         IF ( Vect_var % Dofs>=dim ) THEN
                           DO real_dof=1,Vect_var % Dofs
                             nt_boundary = ASSOCIATED(Var % Values, &
                               Vect_var % Values(real_dof::Vect_var % Dofs))
                             IF ( nt_boundary ) EXIT
                           END DO
                           IF ( nt_boundary ) EXIT
                         END IF
                         Vect_var => Vect_var % Next
                       END DO
                     END IF

                     IF ( nt_boundary ) THEN
                       CALL GetElementNodes(Nodes)
                       nrm = NormalVector(Element,Nodes,0._dp,0._dp,.TRUE.)
                       SELECT CASE(Element % TYPE % DIMENSION)
                       CASE(1)
                         t1(1) =  nrm(2)
                         t1(2) = -nrm(1)
                       CASE(2)
                         CALL TangentDirections(nrm,t1,t2)
                       END SELECT

                       SELECT CASE(vect_dof)
                         CASE(1)
                           vec = nrm
                         CASE(2)
                           vec = t1
                         CASE(3)
                           vec = t2
                       END SELECT
                     END IF
                   END IF
                 END IF

                 DO j=1,n
                   k = Element % NodeIndexes(j)
                   IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                   IF ( k>0 ) THEN
                     IF ( nt_boundary ) THEN
                       DO l=1,dim
                         m = l+real_dof-vect_dof
                         tmp(l)=Vect_var % Values(Vect_var % Dofs*(k-1)+m)
                       END DO
                       udot = SUM(vec(1:dim)*tmp(1:dim))
                       tmp(1:dim)=tmp(1:dim)+(work(j)-udot)*vec(1:dim)
                       DO l=1,dim
                         m = l+real_dof-vect_dof
                         Vect_var % Values(Vect_var % Dofs*(k-1)+m)=tmp(l)
                       END DO
                     ELSE
                       Var % Values(k) = Work(j)
                     END IF
                   END IF
                 END DO
               END IF

               IF ( Transient .AND. Solver % TimeOrder==2 ) THEN
                  Work(1:n) = GetReal( BC, TRIM(Var % Name) // ' Velocity', GotIt )
                  IF ( GotIt ) THEN
                    DO j=1,n
                      k = Element % NodeIndexes(j)
                      IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                      IF ( k>0 ) Var % PrevValues(k,1) = Work(j)
                    END DO
                  END IF
                  Work(1:n) = GetReal( BC, TRIM(Var % Name) // ' Acceleration', GotIt )
                  IF ( GotIt ) THEN
                    DO j=1,n
                      k = Element % NodeIndexes(j)
                      IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                      IF ( k>0 ) Var % PrevValues(k,2) = Work(j)
                    END DO
                  END IF
               END IF
             ELSE
               CALL ListGetRealArray( BC, &
                 Var % Name, WorkA, n, Element % NodeIndexes, gotIt )
               IF ( GotIt ) THEN
                 DO j=1,n
                   k = Element % NodeIndexes(j)
                   IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                   IF(k>0) THEN
                     DO l=1,MIN(SIZE(WorkA,1),Var % DOFs)
                       Var % Values(Var % DOFs*(k-1)+l) = WorkA(l,1,j)
                     END DO
                   END IF
                 END DO
               ELSE
               END IF
             END IF

             IF(NamespaceFound) CALL ListPopNamespace()
             Var => Var % Next
           END DO
         END DO
         DEALLOCATE( Work )
         Mesh => Mesh % Next
       END DO
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE SetInitialConditions
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE InitCond()
!------------------------------------------------------------------------------
     USE DefUtils
     TYPE(Element_t), POINTER :: Edge
     INTEGER :: DOFs,i,j,k,l
     CHARACTER(LEN=MAX_NAME_LEN) :: str
     LOGICAL :: Found, ThingsToDO, NamespaceFound
     TYPE(Solver_t), POINTER :: Solver
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp) :: Val
     REAL(KIND=dp),ALLOCATABLE :: Work(:)
     TYPE(ValueList_t), POINTER :: IC
!------------------------------------------------------------------------------

     Mesh => CurrentModel % Meshes
     DO WHILE( ASSOCIATED( Mesh ) )
       ALLOCATE( Indexes(Mesh % MaxElementDOFs), Work(Mesh % MaxElementDOFs) )
       CALL SetCurrentMesh( CurrentModel, Mesh )

       ! First set the global variables and check whether there is anything left to do
       ThingsToDo = .FALSE.
       DO j=1,CurrentModel % NumberOfICs

         IC => CurrentModel % ICs(j) % Values
         
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 
           
           Solver => Var % Solver
           IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver
           
           str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
           IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
           
           ! global variable
           IF( SIZE( Var % Values ) == Var % DOFs ) THEN
             Val = ListGetCReal( IC, Var % Name, GotIt )
             IF( GotIt ) THEN
               WRITE( Message,'(A,ES12.3)') 'Initializing global variable: > '&
                   //TRIM(Var % Name)//' < to :',Val
               CALL Info('InitCond', Message,Level=8)
               Var % Values = Val
             END IF
           ELSE
             ThingsToDo = ThingsToDo .OR. &
                 ListCheckPresent( IC, TRIM(Var % Name) )
             ThingsToDo = ThingsToDo .OR. &
                 ListCheckPresent( IC, TRIM(Var % Name)//' Velocity' )
             ThingsToDo = ThingsToDo .OR. &
                 ListCheckPresent( IC, TRIM(Var % Name)//' Acceleration' )
             ThingsToDo = ThingsToDo .OR. &
                  ListCheckPresent( IC, TRIM(Var % Name)//' {e}' )
           END IF
           IF (NamespaceFound) CALL ListPopNamespace()
           Var => Var % Next
         END DO
       END DO

       ! And now do the ordinary fields
       !--------------------------------
       IF( ThingsToDo ) THEN
         DO t=1, Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
           
           CurrentElement =>  Mesh % Elements(t)
           
           i = CurrentElement % BodyId 
           IF( i == 0 ) CYCLE
           
           j = ListGetInteger(CurrentModel % Bodies(i) % Values, &
               'Initial Condition',GotIt, 1, CurrentModel % NumberOfICs )           
           IF ( .NOT. GotIt ) CYCLE
           
           IC => CurrentModel % ICs(j) % Values
           CurrentModel % CurrentElement => CurrentElement
           n = GetElementNOFNodes()
           
           Var => Mesh % Variables
           DO WHILE( ASSOCIATED(Var) ) 
             
             
             Solver => Var % Solver
             IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver

             str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
             IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
             
             ! global variables were already set
             IF( SIZE( Var % Values ) == Var % DOFs ) THEN
               CONTINUE
               
             ELSE IF ( Var % DOFs <= 1 ) THEN
               
               Work(1:n) = GetReal( IC, Var % Name, GotIt )
               IF ( GotIt ) THEN
                 DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                 DO k=1,n
                   k1 = Indexes(k)
                   IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                   IF ( k1>0 ) Var % Values(k1) = Work(k)
                 END DO
               END IF
               
               IF ( Transient .AND. Solver % TimeOrder==2 ) THEN
                 Work(1:n) = GetReal( IC, TRIM(Var % Name) // ' Velocity', GotIt )
                 IF ( GotIt ) THEN
                   DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   DO k=1,n
                     k1 = Indexes(k)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % PrevValues(k1,1) = Work(k)
                   END DO
                 END IF
                 Work(1:n) = GetReal( IC, TRIM(Var % Name) // ' Acceleration', GotIt )
                 IF ( GotIt ) THEN
                   DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   DO k=1,n
                     k1 = Indexes(k)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % PrevValues(k1,2) = Work(k)
                   END DO
                 END IF
               END IF
               
               IF(ASSOCIATED(Mesh % Edges)) THEN
                 IF ( i<=Mesh % NumberOfBulkElements) THEN
                   Gotit = ListCheckPresent( IC, TRIM(Var % Name)//' {e}' )
                   IF ( Gotit ) THEN
                     DO k=1,CurrentElement % TYPE % NumberOfedges
                       Edge => Mesh % Edges(CurrentElement % EdgeIndexes(k))
                       l = Var % Perm(CurrentElement % EdgeIndexes(k)+Mesh % NumberOfNodes)
                       IF ( l>0 ) THEN
                         CALL LocalBcIntegral( IC, &
                             Edge, Edge % TYPE % NumberOfNodes, CurrentElement, n, &
                             TRIM(Var % Name)//' {e}', Work(1) )
                         Var % Values(l) = Work(1)
                       END IF
                     END DO
                   END IF
                 END IF
               END IF
               
             ELSE
               CALL ListGetRealArray( IC, &
                   Var % Name, WorkA, n, CurrentElement % NodeIndexes, gotIt )
               
               IF ( GotIt ) THEN
                 DO k=1,n
                   k1 = Indexes(k)
                   IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                   IF(k1>0) THEN
                     DO l=1,MIN(SIZE(WorkA,1),Var % DOFs)
                       IF ( k1>0 ) Var % Values(Var % DOFs*(k1-1)+l) = WorkA(l,1,k)
                     END DO
                   END IF
                 END DO
               END IF
             END IF
             IF(NamespaceFound) CALL ListPopNamespace()
             Var => Var % Next
           END DO
         END DO
       END IF

       DEALLOCATE( Indexes, Work )
       Mesh => Mesh % Next
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE InitCond
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if we are restarting are if yes, read in field values.
!------------------------------------------------------------------------------
   SUBROUTINE Restart()
!------------------------------------------------------------------------------
     USE DefUtils
     LOGICAL :: Gotit
     INTEGER :: k
     REAL(KIND=dp) :: StartTime
!------------------------------------------------------------------------------


     RestartFile = ListGetString( CurrentModel % Simulation, &
         'Restart File', GotIt )

     IF ( GotIt ) THEN
       k = ListGetInteger( CurrentModel % Simulation,'Restart Position',GotIt, &
                  minv=0 )

       Mesh => CurrentModel % Meshes
       DO WHILE( ASSOCIATED(Mesh) ) 

         IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
           OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
         ELSE
           OutputName = TRIM(RestartFile)
         END IF
         IF ( ParEnv % PEs > 1 ) &
           OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))

         CALL SetCurrentMesh( CurrentModel, Mesh )
         CALL LoadRestartFile( OutputName,k,Mesh )

         StartTime = ListGetConstReal( CurrentModel % Simulation,'Restart Time',GotIt)
         IF( GotIt ) THEN
	   Var  => VariableGet( Mesh % Variables, 'Time' )
           IF ( ASSOCIATED( Var ) )  Var % Values(1)  = StartTime
         END IF

         Mesh => Mesh % Next

       END DO


     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE Restart
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Execute the individual solvers in defined sequence. 
!------------------------------------------------------------------------------
   SUBROUTINE ExecSimulation(TimeIntervals,  CoupledMinIter, &
              CoupledMaxIter, OutputIntervals, Transient, Scanning)
      INTEGER :: TimeIntervals,CoupledMinIter, CoupledMaxIter,OutputIntervals(:)
      LOGICAL :: Transient,Scanning
!------------------------------------------------------------------------------
     INTEGER :: interval, timestep, i, j, k, n
     REAL(KIND=dp) :: dt, ddt, dtfunc, timeleft
     INTEGER :: cum_timestep
     INTEGER, SAVE ::  stepcount=0, RealTimestep
     LOGICAL :: ExecThis,SteadyStateReached=.FALSE.

     REAL(KIND=dp) :: CumTime, MaxErr, AdaptiveLimit, &
           AdaptiveMinTimestep, AdaptiveMaxTimestep, timePeriod
     INTEGER :: SmallestCount, AdaptiveKeepSmallest, StepControl=-1
     LOGICAL :: AdaptiveTime = .TRUE., Found

     TYPE(Solver_t), POINTER :: Solver

#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: newtime, prevtime=0, maxtime, exitcond
#else
     REAL(KIND=dp) :: RealTime, newtime, prevtime=0, maxtime, exitcond
#endif
     REAL(KIND=dp), ALLOCATABLE :: xx(:,:), xxnrm(:), yynrm(:), PrevXX(:,:,:)

!$omp parallel
!$   IF(.NOT.GaussPointsInitialized()) CALL GaussPointsInit
!$omp end parallel

     DO i=1,CurrentModel % NumberOfSolvers
        Solver => CurrentModel % Solvers(i)
        IF ( Solver % PROCEDURE==0 ) CYCLE
        IF ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL ) THEN
          ! solver to be called prior to time looping can never be transient
           CALL SolverActivate( CurrentModel,Solver,dt,.FALSE. )
        END IF
     END DO

     IF( ListGetLogical( CurrentModel % Simulation,'Calculate Mesh Pieces',Found ) ) THEN
       CALL CalculateMeshPieces( CurrentModel % Mesh ) 
     END IF
     

     DO interval = 1, TimeIntervals
        stepcount = stepcount + Timesteps(interval)
     END DO 

     cum_Timestep = 0
     ddt = 0.0d0
     DO interval = 1,TimeIntervals

!------------------------------------------------------------------------------
       IF ( Transient .OR. Scanning ) THEN
         dt = TimestepSizes(interval,1)
       ELSE
         dt = 1
       END IF
!------------------------------------------------------------------------------
!      go trough number of timesteps within an interval
!------------------------------------------------------------------------------
       timePeriod = ListGetCReal(CurrentModel % Simulation, 'Time Period',gotIt)
       IF(.NOT.GotIt) timePeriod = HUGE(timePeriod)


       RealTimestep = 1
       DO timestep = 1,Timesteps(interval)

         cum_Timestep = cum_Timestep + 1
         sStep(1) = cum_Timestep

         dtfunc = ListGetConstReal( CurrentModel % Simulation, &
                  'Timestep Function', gotIt)
         IF(GotIt) THEN
	   CALL Warn('ExecSimulation','Obsolite keyword > Timestep Function < , use > Timestep Size < instead')
         ELSE	
           dtfunc = ListGetCReal( CurrentModel % Simulation, &
                  'Timestep Size', gotIt)
         END IF
         IF(GotIt) dt = dtfunc

!------------------------------------------------------------------------------
         sTime(1) = sTime(1) + dt
         sPeriodic(1) = sTime(1)
         DO WHILE(sPeriodic(1) > timePeriod)
           sPeriodic(1) = sPeriodic(1) - timePeriod 
         END DO

         ! Move the old timesteps one step down the ladder
         IF(timestep > 1 .OR. interval > 1) THEN
           DO i = SIZE(sPrevSizes,2),2,-1
             sPrevSizes(1,i) = sPrevSizes(1,i-1)
           END DO
           sPrevSizes(1,1) = sSize(1)
         END IF 
         sSize(1) = dt

         sInterval(1) = interval
         IF (.NOT. Transient ) steadyIt(1) = steadyIt(1) + 1
!------------------------------------------------------------------------------
         IF ( ParEnv % MyPE == 0 ) THEN
           CALL Info( 'MAIN', ' ', Level=3 )
           CALL Info( 'MAIN', '-------------------------------------', Level=3 )

           IF ( Transient .OR. Scanning ) THEN
             WRITE( Message, * ) 'Time: ',TRIM(i2s(cum_Timestep)),'/', &
                   TRIM(i2s(stepcount)), sTime(1)
             CALL Info( 'MAIN', Message, Level=3 )

             newtime= RealTime()

             IF( cum_Timestep > 1 ) THEN
               maxtime = ListGetConstReal( CurrentModel % Simulation,'Real Time Max',GotIt)
               IF( GotIt ) THEN
                 WRITE( Message,'(A,F8.3)') 'Fraction of real time left: ',&
                     1.0_dp-RealTime() / maxtime
               END IF

               ! Compute estimated time left in seconds
               timeleft = (stepcount-(cum_Timestep-1))*(newtime-prevtime)
               
               ! No sense to show too short estimated times
               IF( timeleft > 1 ) THEN
                 IF (timeleft >= 24 * 3600) THEN ! >24 hours
                   WRITE( Message,'(A)') 'Estimated time left: '//I2S(NINT(timeleft/3600))//' hours'
                 ELSE IF (timeleft >= 3600) THEN   ! 1 to 20 hours
                   WRITE( Message,'(A,F5.1,A)') 'Estimated time left:',timeleft/3600,' hours'
                 ELSE IF(timeleft >= 60) THEN ! 1 to 60 minutes
                   WRITE( Message,'(A,F5.1,A)') 'Estimated time left:',timeleft/60,' minutes'
                 ELSE                         ! 1 to 60 seconds
                   WRITE( Message,'(A,F5.1,A)') 'Estimated time left:',timeleft,' seconds'
                 END IF
                 CALL Info( 'MAIN', Message, Level=3 )
               END IF
               
             END IF
             prevtime = newtime
           ELSE
             WRITE( Message, * ) 'Steady state iteration: ',cum_Timestep
             CALL Info( 'MAIN', Message, Level=3 )
           END IF

           CALL Info( 'MAIN', '-------------------------------------', Level=3 )
           CALL Info( 'MAIN', ' ', Level=3 )
         END IF

!------------------------------------------------------------------------------
!        Solve any and all governing equations in the system
!------------------------------------------------------------------------------
         AdaptiveTime = ListGetLogical( CurrentModel % Simulation, &
                  'Adaptive Timestepping', GotIt )

         IF ( Transient .AND. AdaptiveTime ) THEN 
            AdaptiveLimit = ListGetConstReal( CurrentModel % Simulation, &
                        'Adaptive Time Error', GotIt )
 
            IF ( .NOT. GotIt ) THEN 
               WRITE( Message, * ) 'Adaptive Time Limit must be given for' // &
                        'adaptive stepping scheme.'
               CALL Fatal( 'ElmerSolver', Message )
            END IF

            AdaptiveMaxTimestep = ListGetConstReal( CurrentModel % Simulation, &
                     'Adaptive Max Timestep', GotIt )
            IF ( .NOT. GotIt ) AdaptiveMaxTimestep =  dt
            AdaptiveMaxTimestep =  MIN(AdaptiveMaxTimeStep, dt)

            AdaptiveMinTimestep = ListGetConstReal( CurrentModel % Simulation, &
                     'Adaptive Min Timestep', GotIt )

            AdaptiveKeepSmallest = ListGetInteger( CurrentModel % Simulation, &
                       'Adaptive Keep Smallest', GotIt, minv=0  )

            n = CurrentModel % NumberOfSolvers
            j = 0
            k = 0
            DO i=1,n
               Solver => CurrentModel % Solvers(i)
               IF ( ASSOCIATED( Solver % Variable  % Values ) ) THEN
                  IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                     j = MAX( j, SIZE( Solver % Variable % PrevValues,2 ) )
                  END IF
                  k = MAX( k, SIZE( Solver % Variable % Values ) )
               END IF
            END DO
            ALLOCATE( xx(n,k), yynrm(n), xxnrm(n), prevxx( n,k,j ) )

            CumTime = 0.0d0
            IF ( ddt == 0.0d0 .OR. ddt > AdaptiveMaxTimestep ) ddt = AdaptiveMaxTimestep

            s = sTime(1) - dt
            SmallestCount = 0
            DO WHILE( CumTime < dt-1.0d-12 )
               ddt = MIN( dt - CumTime, ddt )

               DO i=1,CurrentModel % NumberOFSolvers
                  Solver => CurrentModel % Solvers(i)
                  IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                     n = SIZE( Solver % Variable % Values )
                     xx(i,1:n) = Solver % Variable % Values
                     xxnrm(i) = Solver % Variable % Norm
                     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                        DO j=1,SIZE( Solver % Variable % PrevValues,2 )
                           prevxx(i,1:n,j) = Solver % Variable % PrevValues(:,j)
                        END DO
                     END IF
                  END IF
               END DO

               sTime(1) = s + CumTime + ddt
               sSize(1) = ddt
               CALL SolveEquations( CurrentModel, ddt, Transient, &
                 CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )


               MaxErr = ListGetConstReal( CurrentModel % Simulation, &
                          'Adaptive Error Measure', GotIt )

               DO i=1,CurrentModel % NumberOFSolvers
                  Solver => CurrentModel % Solvers(i)
                  IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                     n = SIZE(Solver % Variable % Values)
                     yynrm(i) = Solver % Variable % Norm
                     Solver % Variable % Values = xx(i,1:n)
                     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                        DO j=1,SIZE( Solver % Variable % PrevValues,2 )
                           Solver % Variable % PrevValues(:,j) = prevxx(i,1:n,j)
                        END DO
                     END IF
                  END IF
               END DO

               sStep(1) = ddt / 2
               sTime(1) = s + CumTime + ddt/2
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                  CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )
               sTime(1) = s + CumTime + ddt
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                  CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )

               MaxErr = ABS( MaxErr - ListGetConstReal( CurrentModel % Simulation, &
                           'Adaptive Error Measure', GotIt ) )

               IF ( .NOT. GotIt ) THEN
                  MaxErr = 0.0d0
                  DO i=1,CurrentModel % NumberOFSolvers
                     Solver => CurrentModel % Solvers(i)
                     IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                        IF ( yynrm(i) /= Solver % Variable % Norm ) THEN
                           Maxerr = MAX(Maxerr,ABS(yynrm(i)-Solver % Variable % Norm)/yynrm(i))
                        END IF
                     END IF
                  END DO
               END IF

               IF ( MaxErr < AdaptiveLimit .OR. ddt <= AdaptiveMinTimestep ) THEN
                 CumTime = CumTime + ddt
                 RealTimestep = RealTimestep+1
                 IF ( SmallestCount >= AdaptiveKeepSmallest .OR. StepControl > 0 ) THEN
                    ddt = MIN( 2*ddt, AdaptiveMaxTimeStep )
                    StepControl   = 1
                    SmallestCount = 0
                  ELSE
                    StepControl   = 0
                    SmallestCount = SmallestCount + 1
                  END IF
               ELSE
                  DO i=1,CurrentModel % NumberOFSolvers
                     Solver => CurrentModel % Solvers(i)
                     IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                        n = SIZE(Solver % Variable % Values)
                        Solver % Variable % Norm = xxnrm(i)
                        Solver % Variable % Values = xx(i,1:n)
                        IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                           DO j=1,SIZE( Solver % Variable % PrevValues,2 )
                              Solver % Variable % PrevValues(:,j) = prevxx(i,1:n,j)
                           END DO
                        END IF
                     END IF
                  END DO
                  ddt = ddt / 2
                  StepControl = -1
               END IF
               WRITE(*,'(a,3e20.12)') 'Adaptive(cum,ddt,err): ', cumtime, ddt, maxerr
            END DO
            sSize(1) = dt
            sTime(1) = s + dt
  
            DEALLOCATE( xx, xxnrm, yynrm, prevxx )
         ELSE ! Adaptive timestepping
            CALL SolveEquations( CurrentModel, dt, Transient, &
              CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )
            RealTimestep = RealTimestep+1
         END IF
!------------------------------------------------------------------------------
!        Save results to disk, if requested
!------------------------------------------------------------------------------

         LastSaved = .FALSE.
         IF( OutputIntervals(Interval) /= 0 ) THEN

           CALL SaveToPost(0)
           k = MOD( Timestep-1, OutputIntervals(Interval) )
           IF ( k == 0 .OR. SteadyStateReached ) THEN
            
             DO i=1,CurrentModel % NumberOfSolvers
               Solver => CurrentModel % Solvers(i)
               IF ( Solver % PROCEDURE == 0 ) CYCLE
               ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
               When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
               IF ( GotIt ) ExecThis = ( When == 'before saving') 
               IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
             END DO 

             CALL SaveCurrent(Timestep)
             LastSaved = .TRUE.

             DO i=1,CurrentModel % NumberOfSolvers
               Solver => CurrentModel % Solvers(i)
               IF ( Solver % PROCEDURE == 0 ) CYCLE
               ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_SAVE)
               When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
               IF ( GotIt ) ExecThis = ( When == 'after saving') 
               IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
             END DO 
           END IF
         END IF
!------------------------------------------------------------------------------

         maxtime = ListGetCReal( CurrentModel % Simulation,'Real Time Max',GotIt)
         IF( GotIt .AND. RealTime() - RT0 > maxtime ) THEN
            CALL Info('ElmerSolver','Reached allowed maximum real time, exiting...')
            GOTO 100
         END IF

	 exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',GotIt)
	 IF( GotIt .AND. exitcond > 0.0_dp ) THEN
            CALL Info('ElmerSolver','Found a positive exit condition, exiting...')
            GOTO 100
         END IF
	 
!------------------------------------------------------------------------------

         IF ( SteadyStateReached .AND. .NOT. (Transient .OR. Scanning) ) THEN
            IF ( Timestep >= CoupledMinIter ) EXIT
         END IF

!------------------------------------------------------------------------------
       END DO ! timestep within an iterval
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     END DO ! timestep intervals, i.e. the simulation
!------------------------------------------------------------------------------

100   DO i=1,CurrentModel % NumberOfSolvers
        Solver => CurrentModel % Solvers(i)
        IF ( Solver % PROCEDURE == 0 ) CYCLE
        When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
        IF ( GotIt ) THEN
           IF ( When == 'after simulation' .OR. When == 'after all' ) THEN
              CALL SolverActivate( CurrentModel,Solver,dt,Transient )
              IF (ASSOCIATED(Solver % Variable % Values) ) LastSaved = .FALSE.
           END IF
        ELSE
           IF ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_ALL ) THEN
              CALL SolverActivate( CurrentModel,Solver,dt,Transient )
              IF (ASSOCIATED(Solver % Variable % Values) ) LastSaved = .FALSE.
           END IF
        END IF
     END DO

     IF( .NOT. LastSaved ) THEN
        DO i=1,CurrentModel % NumberOfSolvers
           Solver => CurrentModel % Solvers(i)
           IF ( Solver % PROCEDURE == 0 ) CYCLE
           ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
           When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
           IF ( GotIt ) ExecThis = ( When == 'before saving') 
           IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
        END DO
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE ExecSimulation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Saves current timestep to external files.
!------------------------------------------------------------------------------
  SUBROUTINE SaveCurrent( CurrentStep )
!------------------------------------------------------------------------------
    INTEGER :: i, j,k,l,n,q,CurrentStep,nlen
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: EigAnal, GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: Simul
    LOGICAL :: BinaryOutput, SaveAll
    
    Simul = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
    
    OutputFile = ListGetString(CurrentModel % Simulation,'Output File',GotIt)
    IF ( GotIt ) THEN
      IF ( ParEnv % PEs > 1 ) THEN
        DO i=1,MAX_NAME_LEN
          IF ( OutputFile(i:i) == ' ' ) EXIT
        END DO
        OutputFile(i:i) = '.'
        WRITE( OutputFile(i+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
      END IF
      
      BinaryOutput = ListGetLogical( CurrentModel % Simulation,'Binary Output',GotIt )
      IF ( .NOT.GotIt ) BinaryOutput = .FALSE.
      
      SaveAll = .NOT.ListGetLogical( CurrentModel % Simulation,&
          'Omit unchanged variables in output',GotIt )
      IF ( .NOT.GotIt ) SaveAll = .TRUE.
      
      Mesh => CurrentModel % Meshes
      DO WHILE( ASSOCIATED( Mesh ) ) 
        IF ( Mesh % OutputActive ) THEN
          nlen = LEN_TRIM(Mesh % Name )
          IF ( nlen > 0 ) THEN
            OutputName = Mesh % Name(1:nlen) // '/' // TRIM(OutputFile)
          ELSE
            OutputName = OutputFile
          END IF
          
          EigAnal = .FALSE.
          DO i=1,CurrentModel % NumberOfSolvers
            IF ( ASSOCIATED( CurrentModel % Solvers(i) % Mesh, Mesh ) ) THEN
              EigAnal = ListGetLogical( CurrentModel % Solvers(i) % Values, &
                  'Eigen Analysis', GotIt )
              EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % Solvers(i) % Values, &
                  'Harmonic Analysis', GotIt )
              
              IF ( EigAnal ) THEN
                Var => CurrentModel % Solvers(i) % Variable
                IF ( ASSOCIATED(Var % EigenValues) ) THEN
                  IF ( TotalTimesteps == 1 ) THEN
                    DO j=1,CurrentModel % Solvers(i) % NOFEigenValues
                      IF ( CurrentModel % Solvers(i) % Matrix % COMPLEX ) THEN

                        n = SIZE(Var % Values)/Var % DOFs
                        DO k=1,n
                          DO l=1,Var % DOFs/2
                            q = Var % DOFs*(k-1)
                            Var % Values(q+l) = REAL(Var % EigenVectors(j,q/2+l))
                            Var % Values(q+l+Var % DOFs/2) = AIMAG(Var % EigenVectors(j,q/2+l))
                          END DO
                        END DO
                      ELSE
                        Var % Values = REAL( Var % EigenVectors(j,:) )
                      END IF
                      SavedSteps = SaveResult( OutputName, Mesh, &
                          j, sTime(1), BinaryOutput, SaveAll )
                    END DO
                  ELSE
                    j = MIN( CurrentStep, SIZE( Var % EigenVectors,1 ) )
                    IF ( CurrentModel % Solvers(i) % Matrix % COMPLEX ) THEN
                      n = SIZE(Var % Values)/Var % DOFs
                      DO k=1,n
                        DO l=1,Var % DOFs/2
                          q = Var % DOFs*(k-1)
                          Var % Values(q+l) = REAL(Var % EigenVectors(j,q/2+l))
                          Var % Values(q+l+Var % DOFs/2) = AIMAG(Var % EigenVectors(j,q/2+l))
                        END DO
                      END DO
                    ELSE
                      Var % Values = REAL(Var % EigenVectors(j,:))
                    END IF
                    SavedSteps = SaveResult( OutputName, Mesh, &
                        CurrentStep, sTime(1), BinaryOutput, SaveAll )
                  END IF
                  Var % Values = 0.0d0
                END IF
              END IF
            END IF
          END DO
          
          IF ( .NOT. EigAnal ) THEN
            SavedSteps = SaveResult( OutputName,Mesh, NINT(sStep(1)), &
                sTime(1), BinaryOutput, SaveAll )
          END IF
        END IF
        Mesh => Mesh % Next
      END DO
    ELSE
      Mesh => CurrentModel % Meshes
      DO WHILE( ASSOCIATED( Mesh ) ) 
        IF ( Mesh % OutputActive ) &
            Mesh % SavesDone=Mesh % SavesDone+1
        Mesh => Mesh % Next
      END DO
    END IF
    CALL SaveToPost(CurrentStep)
!------------------------------------------------------------------------------
  END SUBROUTINE SaveCurrent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Saves results file to post processing file of ElmerPost format, if requested.
!------------------------------------------------------------------------------
  SUBROUTINE SaveToPost(CurrentStep)
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: EigAnal = .FALSE., Found
    INTEGER :: i, j,k,l,n,q,CurrentStep,nlen,timesteps,SavedEigenValues
    CHARACTER(LEN=MAX_NAME_LEN) :: Simul, SaveWhich
    
    Simul = ListGetString( CurrentModel % Simulation,  'Simulation Type' )

    OutputFile = ListGetString( CurrentModel % Simulation,'Output File',GotIt )
    IF ( Gotit ) THEN
      IF ( ParEnv % PEs > 1 ) THEN
        DO i=1,MAX_NAME_LEN
          IF ( OutputFile(i:i) == ' ' ) EXIT
        END DO
        OutputFile(i:i) = '.'
        WRITE( OutputFile(i+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
      END IF
    END IF
    
    PostFile = ListGetString( CurrentModel % Simulation,'Post File',GotIt )
    IF( .NOT. GotIt ) RETURN

    IF ( ParEnv % PEs > 1 ) THEN
      DO i=1,MAX_NAME_LEN
        IF ( PostFile(i:i) == ' ' ) EXIT
      END DO
      PostFile(i:i) = '.'
      WRITE( PostFile(i+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
    END IF

    ! Loop over all meshes
    !---------------------------------------------------
    Mesh => CurrentModel % Meshes
    DO WHILE( ASSOCIATED( Mesh ) )
      IF ( CurrentStep == 0 .AND. Mesh % SavesDone > 0) THEN
        Mesh => Mesh % Next
        CYCLE
      END IF

      ! Check whether this mesh is active for saving
      !--------------------------------------------------
      IF ( Mesh % OutputActive ) THEN
        nlen = LEN_TRIM(Mesh % Name)
        IF ( nlen==0 .OR. FileNameQualified(OutputFile) ) THEN
          OutputName = OutputFile
        ELSE
          OutputName = Mesh % Name(1:nlen)//'/'//TRIM(OutputFile)
        END IF
        
        IF ( nlen==0 .OR. FileNameQualified(PostFile) ) THEN
          PostName = PostFile
        ELSE
          Postname = Mesh % Name(1:nlen)//'/'//TRIM(PostFile)
        END IF
        
        IF ( ListGetLogical( CurrentModel % Simulation,'Filename Numbering',GotIt) ) THEN
          IF( CurrentStep == 0 ) THEN
            PostName = NextFreeFilename(PostName)
          ELSE 
            PostName = NextFreeFilename(PostName,LastExisting = .TRUE. ) 
          END IF
        END IF

        ! Set the Mesh pointer in the CurrentModel 
        CALL SetCurrentMesh( CurrentModel, Mesh )

        ! Use number of timesteps or number of eigenmodes
        !------------------------------------------------
        EigAnal = .FALSE.
        timesteps = TotalTimeSteps
        DO i=1,CurrentModel % NumberOfSolvers
          IF (ASSOCIATED(CurrentModel % Solvers(i) % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( CurrentModel % &
                Solvers(i) % Values, 'Eigen Analysis', GotIt )
            
            EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % &
                Solvers(i) % Values, 'Harmonic Analysis', GotIt )
            
            IF ( EigAnal ) timesteps = MAX( timesteps, &
                CurrentModel % Solvers(i) % NOFEigenValues )
          END IF
        END DO

        DO i=1,CurrentModel % NumberOfSolvers
          IF (ASSOCIATED(CurrentModel % Solvers(i) % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( CurrentModel % &
                Solvers(i) % Values, 'Eigen Analysis', GotIt )
            
            EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % &
                Solvers(i) % Values, 'Harmonic Analysis', GotIt )
            
            IF ( EigAnal ) THEN
              SaveWhich = ListGetString( CurrentModel % Solvers(i) % Values, &
                  'Eigen and Harmonic Solution Output', Found )
              
              SavedEigenValues = CurrentModel % Solvers(i) % NOFEigenValues
              IF( TotalTimesteps > 1 ) THEN
!                SavedEiegnValues = MIN( CurrentStep, SIZE( Var % EigenVectors,1 ) )
              END IF

              DO j=1, SavedEigenValues
                Var => Mesh % Variables
                DO WHILE(ASSOCIATED(Var))
                  IF ( .NOT. ASSOCIATED(Var % EigenValues) ) THEN
                    Var => Var % Next
                    CYCLE
                  END IF
                  
                  IF ( CurrentModel % Solvers(i) % Matrix % COMPLEX ) THEN
                    IF(Var % DOFs==1) THEN
                      Var => Var % Next
                      CYCLE
                    END IF

                    n = SIZE(Var % Values)/Var % DOFs
                    DO k=1,n
                      DO l=1,Var % DOFs/2
                        q = Var % DOFs*(k-1)
                        Var % Values(q+l) = REAL(Var % EigenVectors(j,q/2+l))
                        Var % Values(q+l+Var % DOFs/2) = AIMAG(Var % EigenVectors(j,q/2+l))
                      END DO
                    END DO

                  ELSE
                    SELECT CASE( SaveWhich )
                    CASE('real part')
                      Var % Values = Var % EigenVectors(j,:)
                    CASE('imag part')
                      Var % Values = AIMAG(Var % EigenVectors(j,:))
                    CASE('abs value')
                      Var % Values = ABS(Var % EigenVectors(j,:))
                    CASE('phase angle')
                      Var % Values = ATAN2(AIMAG(Var % EigenVectors(j,:)), &
                          REAL(Var % EigenVectors(j,:)))
                    CASE DEFAULT
                      Var % CValues => Var % EigenVectors(j,:)
                    END SELECT
                  END IF
                  Var => Var % Next
                END DO
                
                IF ( CurrentStep > 0 ) THEN
                  IF ( Mesh % SavesDone /= 0 ) THEN
                    IF( TotalTimeSteps == 1 ) THEN
                      Mesh % SavesDone = j
                    ELSE
                      Mesh % SavesDone = CurrentStep
                    END IF
                  END IF
                  CALL WritePostFile( PostName,OutputName, CurrentModel, &
                      CurrentModel % Solvers(i) % NOFEigenValues, .TRUE. )
                END IF
              END DO
              EXIT
            END IF
          END IF
        END DO
        
        ! If this mesh has not been saved, then do so
        !----------------------------------------------------------------------------
        IF ( .NOT. EigAnal .OR. CurrentStep == 0 ) THEN
          CALL WritePostFile( PostName, OutputName, CurrentModel, timesteps, .TRUE. )
        END IF

        Var => Mesh % Variables
        DO WHILE(ASSOCIATED(Var))
          IF ( ASSOCIATED(Var % EigenValues) ) THEN
            Var % Values = 0._dp
            Var % CValues => NULL()
          END IF
          Var => Var % Next
        END DO
      END IF
      Mesh => Mesh % Next
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SaveToPost
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE ElmerSolver
!------------------------------------------------------------------------------

!> \} ElmerLib
