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

MODULE ElmerSolver_mod
  
  USE MainUtils
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: ElmerSolver, ElmerSolver_init, ElmerSolver_run, ElmerSolver_runAll, ElmerSolver_finalize



#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp)      :: CT0  
  REAL(KIND=dp),SAVE :: RT0
#else
  REAL(KIND=dp)      :: CT0  
  REAL(KIND=dp),SAVE :: RT0
  INTERFACE
     FUNCTION RealTime()
       USE MainUtils
       REAL(KIND=dp) :: RealTime
     END FUNCTION RealTime
  END INTERFACE
  INTERFACE
     FUNCTION CPUTime()
       USE MainUtils
       REAL(KIND=dp) :: CPUTime
     END FUNCTION CPUTime
  END INTERFACE
#endif

  INTEGER            :: NoArgs

  INTEGER,SAVE       :: iter,Ndeg,istat,nproc,tlen,nthreads

  ! note: this declaration list needs checking for which variables actually
  ! require the SAVE attribute (note also the SAVE attribute is imposed by 
  ! setting a default).  It aslo needs checking for variables that can be 
  ! removed to a subroutine.
  CHARACTER(LEN=MAX_STRING_LEN),SAVE :: threads, CoordTransform  
  REAL(KIND=dp),SAVE            :: ss,dt,dtfunc
  REAL(KIND=dP), POINTER, SAVE  :: WorkA(:,:,:) => NULL()
  REAL(KIND=dp), POINTER, SAVE  :: sTime(:), sStep(:), sInterval(:), sSize(:), &
       steadyIt(:),nonlinIt(:),sPrevSizes(:,:),sPeriodic(:)  
  TYPE(Element_t),POINTER, SAVE :: CurrentElement  
  LOGICAL, SAVE                 :: GotIt,Transient,Scanning,LastSaved  
  INTEGER, SAVE                 :: TimeIntervals,interval,timestep, &
       TotalTimesteps,SavedSteps,CoupledMaxIter,CoupledMinIter  
  INTEGER, POINTER, SAVE        :: Timesteps(:),OutputIntervals(:)
  INTEGER, POINTER, SAVE        :: ActiveSolvers(:)
  REAL(KIND=dp), POINTER, SAVE  :: TimestepSizes(:,:)  
  INTEGER(KIND=AddrInt), SAVE   :: ControlProcedure
  LOGICAL, SAVE                 :: InitDirichlet, ExecThis  
  TYPE(ElementType_t),POINTER,SAVE :: elmt  
  TYPE(ParEnv_t),POINTER,SAVE   :: ParallelEnv  
  CHARACTER(LEN=MAX_NAME_LEN),SAVE :: ModelName, eq, ExecCommand
  CHARACTER(LEN=MAX_STRING_LEN),SAVE :: OutputFile, PostFile, RestartFile, &
       OutputName=' ',PostName=' ', When, OptionString  
  TYPE(Variable_t),POINTER,SAVE :: Var
  TYPE(Mesh_t), POINTER, SAVE   :: Mesh
  TYPE(Solver_t), POINTER, SAVE :: Solver  
  LOGICAL, SAVE                 :: FirstLoad = .TRUE., FirstTime=.TRUE., Found
  LOGICAL, SAVE                 :: Silent, Version, GotModelName  
  INTEGER, SAVE                 :: ExtrudeLevels
  TYPE(Mesh_t), POINTER, SAVE   :: ExtrudedMesh
  INTEGER, SAVE                 :: omp_get_max_threads
  INTEGER, SAVE                 :: Initialize = 0
  
CONTAINS
  
  !------------------------------------------------------------------------------
  ! The top level subroutine for Elmer.  For default stand alone simulations it 
  ! is to be called from Solver.src with Initialize set to 0.  The Elmer API calls 
  ! this routine with different values of initialize.  The ESMF framework calls 
  ! the  ElmerSolver_int, _run and _finalize subroutines directly.
  !------------------------------------------------------------------------------
  SUBROUTINE ElmerSolver(init)
    
    INTEGER :: Init
    
    Initialize = init
    
    CALL ElmerSolver_init()
    IF ( Initialize /= 1 ) THEN       
       CALL ElmerSolver_runAll()
       CALL ElmerSolver_finalize()
    END IF
    
  END SUBROUTINE ElmerSolver
  
  !------------------------------------------------------------------------------
  SUBROUTINE ElmerSolver_init(meshFootprint,ParEnvInitialised)
  
#ifdef HAVE_TRILINOS
    INTERFACE
       SUBROUTINE TrilinosCleanup() BIND(C,name='TrilinosCleanup')
         IMPLICIT NONE
       END SUBROUTINE TrilinosCleanup
    END INTERFACE
#endif

    TYPE(mesh_t),INTENT(INOUT),OPTIONAL :: meshFootprint
    LOGICAL,INTENT(IN),OPTIONAL         :: ParEnvInitialised

    INTEGER :: ii, kk


    ! Start the watches, store later
    !--------------------------------
    RT0 = RealTime()
    CT0 = CPUTime()
    
    
    ! If parallel execution requested, initialize parallel environment, 
    ! unless the calling program already initialised ParEnv.
    !------------------------------------------------------------------
    IF (PRESENT(ParEnvInitialised)) THEN
       IF (ParEnvInitialised) THEN
          WRITE( Message, * ) 'ParEnv initialised by calling program'
          CALL Info( 'ParCommInit', Message, Level=4 )
          ParallelEnv => ParEnv
       ELSE
          ParallelEnv => ParallelInit()
       END IF
    ELSE
       ParallelEnv => ParallelInit()
    END IF
    OutputPE = ParEnv % MyPE
    
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
          DO ii = 1, NoArgs 
#ifdef USE_ISO_C_BINDINGS
             CALL GET_COMMAND_ARGUMENT(ii, OptionString)
#else
             CALL getarg( ii,OptionString )
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
          CALL Info( 'MAIN', 'Library version: ' // VERSION &
#ifdef REVISION
               // ' (Rev: ' // REVISION // ')' )
#else
          )
#endif
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
       
       CALL InitializeElementDescriptions
       FirstTime = .FALSE.
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
       
    IF ( initialize==2 ) GOTO 1
    
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
       CurrentModel => LoadModel( ModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE )
       
       ! Optionally perform simple extrusion to increase the dimension of the mesh
       !----------------------------------------------------------------------------------
       ExtrudeLevels=GetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found)
       IF(ExtrudeLevels>1) THEN
          IF (PRESENT(meshFootprint)) THEN
             meshFootprint = CurrentModel % Meshes
          END IF
          ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, ExtrudeLevels-2)
          DO ii=1,CurrentModel % NumberOfSolvers
             IF(ASSOCIATED(CurrentModel % Solvers(ii) % Mesh,CurrentModel % Meshes)) &
                  CurrentModel % Solvers(ii) % Mesh => ExtrudedMesh 
          END DO
          ExtrudedMesh % Next => CurrentModel % Meshes % Next
          CurrentModel % Meshes => ExtrudedMesh
          
          ! If periodic BC given, compute boundary mesh projector:
          ! ------------------------------------------------------
          DO ii = 1,CurrentModel % NumberOfBCs
             IF(ASSOCIATED(CurrentModel % Bcs(ii) % PMatrix)) &
                  CALL FreeMatrix( CurrentModel % BCs(ii) % PMatrix )
             CurrentModel % BCs(ii) % PMatrix => NULL()
             kk = ListGetInteger( CurrentModel % BCs(ii) % Values, 'Periodic BC', GotIt )
             IF( GotIt ) THEN
                CurrentModel % BCs(ii) % PMatrix =>  PeriodicProjector( CurrentModel, ExtrudedMesh, ii, kk )
             END IF
          END DO
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

       IF ( .NOT.ReloadInputFile(CurrentModel) ) RETURN

       Mesh => CurrentModel % Meshes
       DO WHILE( ASSOCIATED(Mesh) )
          Mesh % SavesDone = 0
          Mesh => Mesh % Next
       END DO
    END IF
       
1   CONTINUE
       
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
         steadyIt(1), nonLinit(1), sPrevSizes(1,5), sPeriodic(1) )
    
    dt   = 0._dp
    
    sTime = 0._dp
    sStep = 0
    sPeriodic = 0._dp
    
    sSize = dt
    sPrevSizes = 0_dp
    
    sInterval = 0._dp
    
    steadyIt = 0
    nonlinIt = 0
    
    CoupledMinIter = ListGetInteger( CurrentModel % Simulation, &
         'Steady State Min Iterations', GotIt )
    
    !------------------------------------------------------------------------------
    !      Add coordinates and simulation time to list of variables so that
    !      coordinate dependent parameter computing routines can ask for
    !      them...
    !------------------------------------------------------------------------------
    IF ( FirstLoad ) CALL AddMeshCoordinatesAndTime
    
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
    
    DO ii=1,CurrentModel % NumberOfSolvers
       Solver => CurrentModel % Solvers(ii)
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

    RETURN

10  CONTINUE
    CALL Fatal( 'ElmerSolver', 'Unable to find ELMERSOLVER_STARTINFO, can not execute.' )

20  CONTINUE
    CALL Fatal( 'ElmerSolver', 'Unable to find input file [' // &
         TRIM(Modelname) // '], can not execute.' )
  
  END SUBROUTINE ElmerSolver_init
    
  !------------------------------------------------------------------------------
  SUBROUTINE  ElmerSolver_runAll()

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

  END SUBROUTINE ElmerSolver_runAll
  
  !------------------------------------------------------------------------------
  SUBROUTINE  ElmerSolver_run()

    CALL ExecSimulation( 1, CoupledMinIter, &
            CoupledMaxIter, OutputIntervals, Transient, Scanning)
  
  END SUBROUTINE ElmerSolver_run
  
  !------------------------------------------------------------------------------
  SUBROUTINE ElmerSolver_finalize()
    
    INTEGER :: ii

    !------------------------------------------------------------------------------
    !    Always save the last step to output
    !------------------------------------------------------------------------------
    IF ( .NOT.LastSaved ) THEN
       DO ii=1,CurrentModel % NumberOfSolvers
          Solver => CurrentModel % Solvers(ii)
          IF ( Solver % PROCEDURE == 0 ) CYCLE
          ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
          When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
          IF ( GotIt ) ExecThis = ( When == 'before saving') 
          IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
       END DO
       
       CALL SaveToPost(0)
       CALL SaveCurrent(Timestep)
       
       DO ii=1,CurrentModel % NumberOfSolvers
          Solver => CurrentModel % Solvers(ii)
          IF ( Solver % PROCEDURE == 0 ) CYCLE
          ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_SAVE)
          When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
          IF ( GotIt ) ExecThis = ( When == 'after saving') 
          IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
       END DO
    END IF
    
    IF (Initialize==0) THEN 
       IF ( .NOT.ReloadInputFile(CurrentModel) ) CONTINUE
    END IF


     CALL CompareToReferenceSolution( Finalize = .TRUE. )

    !------------------------------------------------------------------------------
    !    THIS IS THE END (...,at last, the end, my friend,...)
    !------------------------------------------------------------------------------
    IF ( Initialize /= 1 ) CALL Info( 'ElmerSolver', '*** Elmer Solver: ALL DONE ***',Level=3 )
    
    IF ( Initialize <=0 ) CALL FreeModel(CurrentModel)
    
#ifdef HAVE_TRILINOS
    CALL TrilinosCleanup()
#endif
    
    CALL ParallelFinalize()
!    IF ( ParEnv % PEs>1 )  CALL ParallelFinalize()
    CALL Info('ElmerSolver','The end',Level=3)
            
  END SUBROUTINE ElmerSolver_finalize

  
  

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
       LOGICAL :: Found, Success = .TRUE., FinalizeOnly, CompareNorm, CompareSolution

       SAVE TestCount, PassCount 


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
             OPEN( 10, FILE = 'TEST.PASSED' )
             IF( Success ) THEN
               WRITE( 10,'(I1)' ) 1
             ELSE
               WRITE( 10,'(I1)' ) 0
             END IF
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
           Err = ABS( Norm - RefNorm ) / RefNorm 

           ! Compare to given reference norm
           IF( Err > Tol ) THEN
             ! Warn only in the main core
             IF( ParEnv % MyPe == 0 ) THEN
               WRITE( Message,'(A,I0,A,ES12.6,A,ES12.6)') &
                   'Solver ',solver_id,' FAILED:  Norm = ',Norm,'  RefNorm = ',RefNorm
               CALL Warn('CompareToReferenceSolution',Message)
               WRITE( Message,'(A,ES12.6)') 'Relative Error to reference norm: ',Err
               CALL Info('CompareToReferenceSolution',Message, Level = 4 )
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES12.6,A,ES12.6)') &
                 'Solver ',solver_id,' PASSED:  Norm = ',Norm,'  RefNorm = ',RefNorm
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
               WRITE( Message,'(A,I0,A,ES12.6,A,ES12.6)') &
                   'Solver ',solver_id,' FAILED:  Solution = ',Norm,'  RefSolution = ',RefNorm
               CALL Warn('CompareToReferenceSolution',Message)
               WRITE( Message,'(A,ES12.6)') 'Relative Error to reference solution: ',Err
               CALL Info('CompareToReferenceSolution',Message, Level = 4 )
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES12.6,A,ES12.6)') &
                 'Solver ',solver_id,' PASSED:  Solution = ',Norm,'  RefSolution = ',RefNorm
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
       INTEGER :: ii,jj,j2,j3,kk,nn
       TYPE(ValueList_t), POINTER :: Params
       LOGICAL :: gotIt, VtuFormat

       str = ListGetString( CurrentModel % Simulation,'Post File',GotIt) 
       IF(.NOT. GotIt) RETURN
       kk = INDEX( str,'.vtu' )
       VtuFormat = ( kk /= 0 ) 

       IF(.NOT. VtuFormat ) RETURN
       
       CALL Info('AddVtuOutputSolverHack','Adding ResultOutputSolver to write VTU output in file: '&
           //TRIM(str(1:kk-1)))
     
       CALL ListRemove( CurrentModel % Simulation,'Post File')
       nn = CurrentModel % NumberOfSolvers+1
       ALLOCATE( ABC(nn) )
       DO ii=1,nn-1
         ! Def_Dofs is the only allocatable structure within Solver_t:
         IF( ALLOCATED( CurrentModel % Solvers(ii) % Def_Dofs ) ) THEN
           jj = SIZE(CurrentModel % Solvers(ii) % Def_Dofs,1)
           j2 = SIZE(CurrentModel % Solvers(ii) % Def_Dofs,2)
           j3 = SIZE(CurrentModel % Solvers(ii) % Def_Dofs,3)
           ALLOCATE( ABC(ii) % Def_Dofs(jj,j2,j3) )
         END IF

         ! Copy the content of the Solver structure
         ABC(ii) = CurrentModel % Solvers(ii)

         ! Nullify the old structure since otherwise bad things may happen at deallocation
         NULLIFY( CurrentModel % Solvers(ii) % ActiveElements )
         NULLIFY( CurrentModel % Solvers(ii) % Values )
         NULLIFY( CurrentModel % Solvers(ii) % Mesh )
         NULLIFY( CurrentModel % Solvers(ii) % BlockMatrix )
         NULLIFY( CurrentModel % Solvers(ii) % Matrix )
         NULLIFY( CurrentModel % Solvers(ii) % Variable )
       END DO

       ! Deallocate the old structure and set the pointer to the new one
       DEALLOCATE( CurrentModel % Solvers )
       CurrentModel % Solvers => ABC
       CurrentModel % NumberOfSolvers = nn

       ! Now create the ResultOutputSolver instance on-the-fly
       NULLIFY( CurrentModel % Solvers(nn) % Values )
       CurrentModel % Solvers(nn) % PROCEDURE = 0
       NULLIFY( CurrentModel % Solvers(nn) % Matrix )
       NULLIFY( CurrentModel % Solvers(nn) % BlockMatrix )
       NULLIFY( CurrentModel % Solvers(nn) % Values )
       NULLIFY( CurrentModel % Solvers(nn) % Variable )
       NULLIFY( CurrentModel % Solvers(nn) % ActiveElements )
       CurrentModel % Solvers(nn) % NumberOfActiveElements = 0
       NULLIFY( CurrentModel % Solvers(nn) % Values )
       jj = CurrentModel % NumberOfBodies
       ALLOCATE( CurrentModel % Solvers(nn) % Def_Dofs(10,jj,6))
       CurrentModel % Solvers(nn) % Def_Dofs(:,1:jj,6) = -1
       
       ! Add some keywords to the list
       CALL ListAddString(CurrentModel % Solvers(nn) % Values,&
           'Procedure', 'ResultOutputSolve ResultOutputSolver',.FALSE.)
       CALL ListAddString(CurrentModel % Solvers(nn) % Values,'Output Format','vtu')
       CALL ListAddString(CurrentModel % Solvers(nn) % Values,'Output File Name',str(1:kk-1))
       CALL ListAddString(CurrentModel % Solvers(nn) % Values,'Exec Solver','after saving')
       CALL ListAddString(CurrentModel % Solvers(nn) % Values,'Equation','InternalVtuOutputSolver')
       CALL ListAddLogical(CurrentModel % Solvers(nn) % Values,'Save Geometry IDs',.TRUE.)

     END SUBROUTINE AddVtuOutputSolverHack


!------------------------------------------------------------------------------
!> Adds flags for active solvers. 
!------------------------------------------------------------------------------
    SUBROUTINE AddSolvers()
!------------------------------------------------------------------------------
      INTEGER :: ii,jj,kk,nlen
      LOGICAL :: InitSolver, Found
!------------------------------------------------------------------------------

      ! This is a hack that sets Equation flags True for the "Active Solvers".
      ! The Equation flag is the legacy way of setting a Solver active and is still
      ! used internally.
      !----------------------------------------------------------------------------
      DO ii=1,CurrentModel % NumberOfSolvers

        eq = ListGetString( CurrentModel % Solvers(ii) % Values,'Equation', Found )
     
        IF ( Found ) THEN
          nlen = LEN_TRIM(eq)
          DO jj=1,CurrentModel % NumberOFEquations
             ActiveSolvers => ListGetIntegerArray( CurrentModel % Equations(jj) % Values, &
                                'Active Solvers', Found )
             IF ( Found ) THEN
                DO kk=1,SIZE(ActiveSolvers)
                   IF ( ActiveSolvers(kk) == ii ) THEN
!------------------------------------------------------------------------------------------
                      CALL ListAddLogical( CurrentModel % Equations(jj) % Values, eq(1:nlen), .TRUE. )
                      EXIT
                   END IF
                END DO
             END IF
          END DO
       END IF
     END DO

     DO ii=1,CurrentModel % NumberOfSolvers
        eq = ListGetString( CurrentModel % Solvers(ii) % Values,'Equation', Found )

        Solver => CurrentModel % Solvers(ii)
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
!------------------------------------------------------------------------------
  END SUBROUTINE AddSolvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds coordinate and time variables to the current mesh structure. 
!------------------------------------------------------------------------------
  SUBROUTINE AddMeshCoordinatesAndTime()
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: DtVar

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
     LOGICAL :: Found
     TYPE(Solver_t), POINTER :: Solver
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp),ALLOCATABLE :: Work(:)

     INTEGER :: ii,jj,kk,ll,nn,mm,tt,vect_dof,real_dof,dim

     REAL(KIND=dp) :: nrm(3),t1(3),t2(3),vec(3),tmp(3),udot
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Nodes_t), SAVE :: Nodes
     LOGICAL :: nt_boundary
     TYPE(Element_t), POINTER :: Element
     TYPE(Variable_t), POINTER :: var, vect_var

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

         DO tt = Mesh % NumberOfBulkElements + 1, &
                 Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           Element => Mesh % Elements(tt)

           ! Set also the current element pointer in the model structure to
           ! reflect the element being processed:
           ! ---------------------------------------------------------------
           CurrentModel % CurrentElement => Element
           nn = Element % TYPE % NumberOfNodes

           BC => GetBC()

           Var => Mesh % Variables
           DO WHILE( ASSOCIATED(Var) )
             Solver => Var % Solver
             IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver

             str = ListGetString( Solver % Values, 'Namespace', Found )
             IF (Found) CALL ListSetNamespace(TRIM(str))


             IF ( Var % DOFs <= 1 ) THEN
               Work(1:nn) = GetReal( BC,Var % Name, gotIt )
               IF ( GotIt ) THEN

                 nt_boundary = .FALSE.
                 IF ( GetElementFamily() /= 1 ) THEN
                   kk = LEN_TRIM(var % name)
                   vect_dof = ICHAR(Var % Name(kk:kk))-ICHAR('0');
                   IF ( vect_dof>=1 .AND. vect_dof<= 3 ) THEN
                     nt_boundary =  GetLogical( BC, &
                        'normal-tangential '//var % name(1:kk-2), gotIt)

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

                 DO jj=1,nn
                   kk = Element % NodeIndexes(jj)
                   IF ( ASSOCIATED(Var % Perm) ) kk = Var % Perm(kk)
                   IF ( kk>0 ) THEN
                     IF ( nt_boundary ) THEN
                       DO ll=1,dim
                         mm = ll+real_dof-vect_dof
                         tmp(ll)=Vect_var % Values(Vect_var % Dofs*(kk-1)+mm)
                       END DO
                       udot = SUM(vec(1:dim)*tmp(1:dim))
                       tmp(1:dim)=tmp(1:dim)+(work(jj)-udot)*vec(1:dim)
                       DO ll=1,dim
                         mm = ll+real_dof-vect_dof
                         Vect_var % Values(Vect_var % Dofs*(kk-1)+mm)=tmp(ll)
                       END DO
                     ELSE
                       Var % Values(kk) = Work(jj)
                     END IF
                   END IF
                 END DO
               END IF

               IF ( Transient .AND. Solver % TimeOrder==2 ) THEN
                  Work(1:nn) = GetReal( BC, TRIM(Var % Name) // ' Velocity', GotIt )
                  IF ( GotIt ) THEN
                    DO jj=1,nn
                      kk = Element % NodeIndexes(jj)
                      IF ( ASSOCIATED(Var % Perm) ) kk = Var % Perm(kk)
                      IF ( kk>0 ) Var % PrevValues(kk,1) = Work(jj)
                    END DO
                  END IF
                  Work(1:nn) = GetReal( BC, TRIM(Var % Name) // ' Acceleration', GotIt )
                  IF ( GotIt ) THEN
                    DO jj=1,nn
                      kk = Element % NodeIndexes(jj)
                      IF ( ASSOCIATED(Var % Perm) ) kk = Var % Perm(kk)
                      IF ( kk>0 ) Var % PrevValues(kk,2) = Work(jj)
                    END DO
                  END IF
               END IF
             ELSE
               CALL ListGetRealArray( BC, &
                 Var % Name, WorkA, nn, Element % NodeIndexes, gotIt )
               IF ( GotIt ) THEN
                 DO jj=1,nn
                   kk = Element % NodeIndexes(jj)
                   DO ll=1,MIN(SIZE(WorkA,1),Var % DOFs)
                     IF ( ASSOCIATED(Var % Perm) ) kk = Var % Perm(kk)
                     IF ( kk>0 ) Var % Values(Var % DOFs*(kk-1)+ll) = WorkA(ll,1,jj)
                   END DO
                 END DO
               ELSE
               END IF
             END IF

             CALL ListSetNamespace('')
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
     INTEGER :: DOFs,bid,jj,kk,ll,k1,nn,tt
     CHARACTER(LEN=MAX_NAME_LEN) :: str
     LOGICAL :: Found, ThingsToDO
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
       DO jj=1,CurrentModel % NumberOfICs

         IC => CurrentModel % ICs(jj) % Values
         
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 
           
           Solver => Var % Solver
           IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver
           
           str = ListGetString( Solver % Values, 'Namespace', Found )
           IF (Found) CALL ListSetNamespace(TRIM(str))
           
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
           Var => Var % Next
         END DO
       END DO

       ! And now do the ordinary fields
       !--------------------------------
       IF( ThingsToDo ) THEN
         DO tt=1, Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
           
           CurrentElement =>  Mesh % Elements(tt)
           
           bid = CurrentElement % BodyId 
           IF( bid == 0 ) CYCLE
           
           jj = ListGetInteger(CurrentModel % Bodies(bid) % Values, &
               'Initial Condition',GotIt, 1, CurrentModel % NumberOfICs )           
           IF ( .NOT. GotIt ) CYCLE
           
           IC => CurrentModel % ICs(jj) % Values
           CurrentModel % CurrentElement => CurrentElement
           nn = GetElementNOFNodes()
           
           Var => Mesh % Variables
           DO WHILE( ASSOCIATED(Var) ) 
             
             
             Solver => Var % Solver
             IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver

             str = ListGetString( Solver % Values, 'Namespace', Found )
             IF (Found) CALL ListSetNamespace(TRIM(str))
             
             ! global variables were already set
             IF( SIZE( Var % Values ) == Var % DOFs ) THEN
               CONTINUE
               
             ELSE IF ( Var % DOFs <= 1 ) THEN
               
               Work(1:nn) = GetReal( IC, Var % Name, GotIt )
               IF ( GotIt ) THEN
                 DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                 DO kk=1,nn
                   k1 = Indexes(kk)
                   IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                   IF ( k1>0 ) Var % Values(k1) = Work(kk)
                 END DO
               END IF
               
               IF ( Transient .AND. Solver % TimeOrder==2 ) THEN
                 Work(1:nn) = GetReal( IC, TRIM(Var % Name) // ' Velocity', GotIt )
                 IF ( GotIt ) THEN
                   DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   DO kk=1,nn
                     k1 = Indexes(kk)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % PrevValues(k1,1) = Work(kk)
                   END DO
                 END IF
                 Work(1:nn) = GetReal( IC, TRIM(Var % Name) // ' Acceleration', GotIt )
                 IF ( GotIt ) THEN
                   DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   DO kk=1,nn
                     k1 = Indexes(kk)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % PrevValues(k1,2) = Work(kk)
                   END DO
                 END IF
               END IF
               
               IF(ASSOCIATED(Mesh % Edges)) THEN
                 IF ( bid<=Mesh % NumberOfBulkElements) THEN
                   Gotit = ListCheckPresent( IC, TRIM(Var % Name)//' {e}' )
                   IF ( Gotit ) THEN
                     DO kk=1,CurrentElement % TYPE % NumberOfedges
                       Edge => Mesh % Edges(CurrentElement % EdgeIndexes(kk))
                       ll = Var % Perm(CurrentElement % EdgeIndexes(kk)+Mesh % NumberOfNodes)
                       IF ( ll>0 ) THEN
                         CALL LocalBcIntegral( IC, &
                             Edge, Edge % TYPE % NumberOfNodes, CurrentElement, nn, &
                             TRIM(Var % Name)//' {e}', Work(1) )
                         Var % Values(ll) = Work(1)
                       END IF
                     END DO
                   END IF
                 END IF
               END IF
               
             ELSE
               CALL ListGetRealArray( IC, &
                   Var % Name, WorkA, nn, CurrentElement % NodeIndexes, gotIt )
               
               IF ( GotIt ) THEN
                 DO kk=1,nn
                   k1 = Indexes(kk)
                   DO ll=1,MIN(SIZE(WorkA,1),Var % DOFs)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % Values(Var % DOFs*(k1-1)+ll) = WorkA(ll,1,kk)
                   END DO
                 END DO
               END IF
             END IF
             CALL ListSetNamespace('')
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
     INTEGER :: kk,StartTime
!------------------------------------------------------------------------------


     RestartFile = ListGetString( CurrentModel % Simulation, &
         'Restart File', GotIt )

     IF ( GotIt ) THEN
       kk = ListGetInteger( CurrentModel % Simulation,'Restart Position',GotIt, &
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
         CALL LoadRestartFile( OutputName,kk,Mesh )

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
     INTEGER :: interval, timestep, ii, jj, kk, nn
     REAL(KIND=dp) :: dt, ddt, dtfunc
     INTEGER :: timeleft,cum_timestep
     INTEGER, SAVE ::  stepcount=0, RealTimestep
     LOGICAL :: ExecThis,SteadyStateReached=.FALSE.

     REAL(KIND=dp) :: CumTime, MaxErr, AdaptiveLimit, &
           AdaptiveMinTimestep, AdaptiveMaxTimestep, timePeriod
     INTEGER :: SmallestCount, AdaptiveKeepSmallest, StepControl=-1
     LOGICAL :: AdaptiveTime = .TRUE.

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

     DO ii=1,CurrentModel % NumberOfSolvers
        Solver => CurrentModel % Solvers(ii)
        IF ( Solver % PROCEDURE==0 ) CYCLE
        IF ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL ) THEN
           CALL SolverActivate( CurrentModel,Solver,dt,Transient )
        END IF
     END DO

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
           DO ii = SIZE(sPrevSizes,2),2,-1
             sPrevSizes(1,ii) = sPrevSizes(1,ii-1)
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
               ELSE             
                 timeleft = NINT((stepcount-(cum_Timestep-1))*(newtime-prevtime)/60._dp);
                 IF (timeleft > 120) THEN
                   WRITE( Message, *) 'Estimated time left: ', &
                     TRIM(i2s(timeleft/60)),' hours.'
                 ELSE IF(timeleft > 60) THEN
                   WRITE( Message, *) 'Estimated time left: 1 hour ', &
                     TRIM(i2s(MOD(timeleft,60))), ' minutes.'
                 ELSE IF(timeleft >= 1) THEN
                   WRITE( Message, *) 'Estimated time left: ', &
                     TRIM(i2s(timeleft)),' minutes.'
                 ELSE
                   WRITE( Message, *) 'Estimated time left: less than a minute.'
                 END IF
               END IF
               CALL Info( 'MAIN', Message, Level=3 )
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

            nn = CurrentModel % NumberOfSolvers
            jj = 0
            kk = 0
            DO ii=1,nn
               Solver => CurrentModel % Solvers(ii)
               IF ( ASSOCIATED( Solver % Variable  % Values ) ) THEN
                  IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                     jj = MAX( jj, SIZE( Solver % Variable % PrevValues,2 ) )
                  END IF
                  kk = MAX( kk, SIZE( Solver % Variable % Values ) )
               END IF
            END DO
            ALLOCATE( xx(nn,kk), yynrm(nn), xxnrm(nn), prevxx( nn,kk,jj ) )

            CumTime = 0.0d0
            IF ( ddt == 0.0d0 .OR. ddt > AdaptiveMaxTimestep ) ddt = AdaptiveMaxTimestep

            ss = sTime(1) - dt
            SmallestCount = 0
            DO WHILE( CumTime < dt-1.0d-12 )
               ddt = MIN( dt - CumTime, ddt )

               DO ii=1,CurrentModel % NumberOFSolvers
                  Solver => CurrentModel % Solvers(ii)
                  IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                     nn = SIZE( Solver % Variable % Values )
                     xx(ii,1:nn) = Solver % Variable % Values
                     xxnrm(ii) = Solver % Variable % Norm
                     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                        DO jj=1,SIZE( Solver % Variable % PrevValues,2 )
                           prevxx(ii,1:nn,jj) = Solver % Variable % PrevValues(:,jj)
                        END DO
                     END IF
                  END IF
               END DO

               sTime(1) = ss + CumTime + ddt
               sSize(1) = ddt
               CALL SolveEquations( CurrentModel, ddt, Transient, &
                 CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )


               MaxErr = ListGetConstReal( CurrentModel % Simulation, &
                          'Adaptive Error Measure', GotIt )

               DO ii=1,CurrentModel % NumberOFSolvers
                  Solver => CurrentModel % Solvers(ii)
                  IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                     nn = SIZE(Solver % Variable % Values)
                     yynrm(ii) = Solver % Variable % Norm
                     Solver % Variable % Values = xx(ii,1:nn)
                     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                        DO jj=1,SIZE( Solver % Variable % PrevValues,2 )
                           Solver % Variable % PrevValues(:,jj) = prevxx(ii,1:nn,jj)
                        END DO
                     END IF
                  END IF
               END DO

               sStep(1) = ddt / 2
               sTime(1) = ss + CumTime + ddt/2
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                  CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )
               sTime(1) = ss + CumTime + ddt
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                  CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep )

               MaxErr = ABS( MaxErr - ListGetConstReal( CurrentModel % Simulation, &
                           'Adaptive Error Measure', GotIt ) )

               IF ( .NOT. GotIt ) THEN
                  MaxErr = 0.0d0
                  DO ii=1,CurrentModel % NumberOFSolvers
                     Solver => CurrentModel % Solvers(ii)
                     IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                        IF ( yynrm(ii) /= Solver % Variable % Norm ) THEN
                           Maxerr = MAX(Maxerr,ABS(yynrm(ii)-Solver % Variable % Norm)/yynrm(ii))
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
                  DO ii=1,CurrentModel % NumberOFSolvers
                     Solver => CurrentModel % Solvers(ii)
                     IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                        nn = SIZE(Solver % Variable % Values)
                        Solver % Variable % Norm = xxnrm(ii)
                        Solver % Variable % Values = xx(ii,1:nn)
                        IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                           DO jj=1,SIZE( Solver % Variable % PrevValues,2 )
                              Solver % Variable % PrevValues(:,jj) = prevxx(ii,1:nn,jj)
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
            sTime(1) = ss + dt
  
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
           kk = MOD( Timestep-1, OutputIntervals(Interval) )
           IF ( kk == 0 .OR. SteadyStateReached ) THEN
            
             DO ii=1,CurrentModel % NumberOfSolvers
               Solver => CurrentModel % Solvers(ii)
               IF ( Solver % PROCEDURE == 0 ) CYCLE
               ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
               When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
               IF ( GotIt ) ExecThis = ( When == 'before saving') 
               IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
             END DO 

             CALL SaveCurrent(Timestep)
             LastSaved = .TRUE.

             DO ii=1,CurrentModel % NumberOfSolvers
               Solver => CurrentModel % Solvers(ii)
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

100   DO ii=1,CurrentModel % NumberOfSolvers
        Solver => CurrentModel % Solvers(ii)
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
        DO ii=1,CurrentModel % NumberOfSolvers
           Solver => CurrentModel % Solvers(ii)
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
    INTEGER :: ii, jj,kk,ll,nn,qq,CurrentStep,nlen
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: EigAnal, GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: Simul
    LOGICAL :: BinaryOutput, SaveAll
    
    Simul = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
    
    OutputFile = ListGetString(CurrentModel % Simulation,'Output File',GotIt)
    IF ( GotIt ) THEN
      IF ( ParEnv % PEs > 1 ) THEN
        DO ii=1,MAX_NAME_LEN
          IF ( OutputFile(ii:ii) == ' ' ) EXIT
        END DO
        OutputFile(ii:ii) = '.'
        WRITE( OutputFile(ii+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
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
          DO ii=1,CurrentModel % NumberOfSolvers
            IF ( ASSOCIATED( CurrentModel % Solvers(ii) % Mesh, Mesh ) ) THEN
              EigAnal = ListGetLogical( CurrentModel % Solvers(ii) % Values, &
                  'Eigen Analysis', GotIt )
              EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % Solvers(ii) % Values, &
                  'Harmonic Analysis', GotIt )
              
              IF ( EigAnal ) THEN
                Var => CurrentModel % Solvers(ii) % Variable
                IF ( ASSOCIATED(Var % EigenValues) ) THEN
                  IF ( TotalTimesteps == 1 ) THEN
                    DO jj=1,CurrentModel % Solvers(ii) % NOFEigenValues
                      IF ( CurrentModel % Solvers(ii) % Matrix % COMPLEX ) THEN

                        nn = SIZE(Var % Values)/Var % DOFs
                        DO kk=1,nn
                          DO ll=1,Var % DOFs/2
                            qq = Var % DOFs*(kk-1)
                            Var % Values(qq+ll) = REAL(Var % EigenVectors(jj,qq/2+ll))
                            Var % Values(qq+ll+Var % DOFs/2) = AIMAG(Var % EigenVectors(jj,qq/2+ll))
                          END DO
                        END DO
                      ELSE
                        Var % Values = REAL( Var % EigenVectors(jj,:) )
                      END IF
                      SavedSteps = SaveResult( OutputName, Mesh, &
                          jj, sTime(1), BinaryOutput, SaveAll )
                    END DO
                  ELSE
                    jj = MIN( CurrentStep, SIZE( Var % EigenVectors,1 ) )
                    IF ( CurrentModel % Solvers(ii) % Matrix % COMPLEX ) THEN
                      nn = SIZE(Var % Values)/Var % DOFs
                      DO kk=1,nn
                        DO ll=1,Var % DOFs/2
                          qq = Var % DOFs*(kk-1)
                          Var % Values(qq+ll) = REAL(Var % EigenVectors(jj,qq/2+ll))
                          Var % Values(qq+ll+Var % DOFs/2) = AIMAG(Var % EigenVectors(jj,qq/2+ll))
                        END DO
                      END DO
                    ELSE
                      Var % Values = REAL(Var % EigenVectors(jj,:))
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
    INTEGER :: ii, jj,kk,ll,nn,qq,CurrentStep,nlen,timesteps,SavedEigenValues
    CHARACTER(LEN=MAX_NAME_LEN) :: Simul, SaveWhich
    
    Simul = ListGetString( CurrentModel % Simulation,  'Simulation Type' )

    OutputFile = ListGetString( CurrentModel % Simulation,'Output File',GotIt )
    IF ( Gotit ) THEN
      IF ( ParEnv % PEs > 1 ) THEN
        DO ii=1,MAX_NAME_LEN
          IF ( OutputFile(ii:ii) == ' ' ) EXIT
        END DO
        OutputFile(ii:ii) = '.'
        WRITE( OutputFile(ii+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
      END IF
    END IF
    
    PostFile = ListGetString( CurrentModel % Simulation,'Post File',GotIt )
    IF( .NOT. GotIt ) RETURN

    IF ( ParEnv % PEs > 1 ) THEN
      DO ii=1,MAX_NAME_LEN
        IF ( PostFile(ii:ii) == ' ' ) EXIT
      END DO
      PostFile(ii:ii) = '.'
      WRITE( PostFile(ii+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
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
        DO ii=1,CurrentModel % NumberOfSolvers
          IF (ASSOCIATED(CurrentModel % Solvers(ii) % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( CurrentModel % &
                Solvers(ii) % Values, 'Eigen Analysis', GotIt )
            
            EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % &
                Solvers(ii) % Values, 'Harmonic Analysis', GotIt )
            
            IF ( EigAnal ) timesteps = MAX( timesteps, &
                CurrentModel % Solvers(ii) % NOFEigenValues )
          END IF
        END DO

        DO ii=1,CurrentModel % NumberOfSolvers
          IF (ASSOCIATED(CurrentModel % Solvers(ii) % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( CurrentModel % &
                Solvers(ii) % Values, 'Eigen Analysis', GotIt )
            
            EigAnal = EigAnal .OR. ListGetLogical( CurrentModel % &
                Solvers(ii) % Values, 'Harmonic Analysis', GotIt )
            
            IF ( EigAnal ) THEN
              SaveWhich = ListGetString( CurrentModel % Solvers(ii) % Values, &
                  'Eigen and Harmonic Solution Output', Found )
              
              SavedEigenValues = CurrentModel % Solvers(ii) % NOFEigenValues
              IF( TotalTimesteps > 1 ) THEN
!                SavedEiegnValues = MIN( CurrentStep, SIZE( Var % EigenVectors,1 ) )
              END IF

              DO jj=1, SavedEigenValues
                Var => Mesh % Variables
                DO WHILE(ASSOCIATED(Var))
                  IF ( .NOT. ASSOCIATED(Var % EigenValues) ) THEN
                    Var => Var % Next
                    CYCLE
                  END IF
                  
                  IF ( CurrentModel % Solvers(ii) % Matrix % COMPLEX ) THEN
                    IF(Var % DOFs==1) THEN
                      Var => Var % Next
                      CYCLE
                    END IF

                    nn = SIZE(Var % Values)/Var % DOFs
                    DO kk=1,nn
                      DO ll=1,Var % DOFs/2
                        qq = Var % DOFs*(kk-1)
                        Var % Values(qq+ll) = REAL(Var % EigenVectors(jj,qq/2+ll))
                        Var % Values(qq+ll+Var % DOFs/2) = AIMAG(Var % EigenVectors(jj,qq/2+ll))
                      END DO
                    END DO

                  ELSE
                    SELECT CASE( SaveWhich )
                    CASE('real part')
                      Var % Values = Var % EigenVectors(jj,:)
                    CASE('imag part')
                      Var % Values = AIMAG(Var % EigenVectors(jj,:))
                    CASE('abs value')
                      Var % Values = ABS(Var % EigenVectors(jj,:))
                    CASE('phase angle')
                      Var % Values = ATAN2(AIMAG(Var % EigenVectors(jj,:)), &
                          REAL(Var % EigenVectors(jj,:)))
                    CASE DEFAULT
                      Var % CValues => Var % EigenVectors(jj,:)
                    END SELECT
                  END IF
                  Var => Var % Next
                END DO
                
                IF ( CurrentStep > 0 ) THEN
                  IF ( Mesh % SavesDone /= 0 ) THEN
                    IF( TotalTimeSteps == 1 ) THEN
                      Mesh % SavesDone = jj
                    ELSE
                      Mesh % SavesDone = CurrentStep
                    END IF
                  END IF
                  CALL WritePostFile( PostName,OutputName, CurrentModel, &
                      CurrentModel % Solvers(ii) % NOFEigenValues, .TRUE. )
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
END MODULE ElmerSolver_mod
!------------------------------------------------------------------------------

!> \} ElmerLib
