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

#include "../config.h"

!------------------------------------------------------------------------------
!> The main program for Elmer. Solves the equations as defined by the input files.
!------------------------------------------------------------------------------
   SUBROUTINE ElmerSolver(initialize)
!------------------------------------------------------------------------------

     USE Lists
     USE MainUtils
     USE OptimizationUtils
     USE SolverUtils, ONLY: GetControlValue
    
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------

     INTEGER :: Initialize

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     INTEGER :: i,j,k,n,l,t,k1,k2,iter,Ndeg,istat,nproc,tlen,nthreads
     CHARACTER(LEN=MAX_STRING_LEN) :: threads
     CHARACTER(:), ALLOCATABLE :: CoordTransform

     REAL(KIND=dp) :: s,dt,dtfunc
     REAL(KIND=dP), POINTER :: WorkA(:,:,:) => NULL()
     REAL(KIND=dp), POINTER, SAVE :: sTime(:), sStep(:), sInterval(:), sSize(:), &
         steadyIt(:),nonlinIt(:),sPrevSizes(:,:),sPeriodicTime(:),sPeriodicCycle(:),&
         sScan(:),sSweep(:),sPar(:),sFinish(:),sProduce(:),sSlice(:),sSliceRatio(:),&
         sSliceWeight(:), sAngle(:), sAngleVelo(:),sSector(:)

     LOGICAL :: GotIt,Transient,Scanning, LastSaved, MeshMode = .FALSE.

     INTEGER :: TimeIntervals,interval,timestep, &
       TotalTimesteps,SavedSteps,CoupledMaxIter,CoupledMinIter

     INTEGER, POINTER, SAVE :: Timesteps(:),OutputIntervals(:) => NULL(), ActiveSolvers(:)
     REAL(KIND=dp), POINTER, SAVE :: TimestepSizes(:,:),TimestepRatios(:,:)

     INTEGER(KIND=AddrInt) :: ControlProcedure

     LOGICAL :: InitDirichlet, ExecThis, GotTimestepRatios = .FALSE.

     TYPE(ElementType_t),POINTER :: elmt

     TYPE(ParEnv_t), POINTER :: ParallelEnv

     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName, eq
     CHARACTER(LEN=MAX_STRING_LEN) :: OptionString

     CHARACTER(:), ALLOCATABLE :: str, PostFile, ExecCommand, OutputFile, RestartFile, &
          OutputName, PostName, When

     TYPE(Variable_t), POINTER :: Var
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: Solver, iSolver

     REAL(KIND=dp) :: CT0,RT0,tt

     LOGICAL :: Silent=.FALSE., Version=.FALSE., GotModelName, FinishEarly=.FALSE.
     LOGICAL :: FirstLoad = .TRUE., FirstTime=.TRUE., Found

     INTEGER :: iargc, NoArgs
     INTEGER :: iostat, iSweep = 1, OptimIters
     LOGICAL :: GotOptimIters
     INTEGER :: MeshIndex
     TYPE(Mesh_t), POINTER :: ExtrudedMesh

     TYPE(Model_t), POINTER, SAVE :: Control
     LOGICAL :: DoControl=.FALSE., ProcControl=.FALSE., GotParams=.FALSE., DoIt
     INTEGER :: nr,ni,ExtMethod
     INTEGER, ALLOCATABLE :: ipar(:)
     REAL(KIND=dp), ALLOCATABLE :: rpar(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: MeshDir, MeshName
     
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
       ! -----------------------
       NoArgs = COMMAND_ARGUMENT_COUNT()
       ! Info Level is always true until the model has been read!
       ! This makes it possible to cast something 
       Silent = .FALSE.
       Version = .FALSE.

       IF( NoArgs > 0 ) THEN 
         i = 0
         DO WHILE( i < NoArgs )
           i = i + 1 
           CALL GET_COMMAND_ARGUMENT(i, OptionString)
           IF( OptionString=='-rpar' ) THEN
             ! Followed by number of parameters + the parameter values
             i = i + 1
             CALL GET_COMMAND_ARGUMENT(i, OptionString)
             READ( OptionString,*) nr             
             ALLOCATE( rpar(nr) )
             DO j=1,nr
               i = i + 1
               CALL GET_COMMAND_ARGUMENT(i, OptionString)
               READ( OptionString,*) rpar(j)
             END DO
             CALL Info('MAIN','Read '//I2S(nr)//' real parameters from command line!')
             CALL SetRealParametersMATC(nr,rpar)
           END IF

           IF( OptionString=='-ipar' ) THEN
             ! Followed by number of parameters + the parameter values
             i = i + 1
             CALL GET_COMMAND_ARGUMENT(i, OptionString)
             READ( OptionString,*) ni             
             ALLOCATE( ipar(nr) )
             DO j=1,ni
               i = i + 1
               CALL GET_COMMAND_ARGUMENT(i, OptionString)
               READ( OptionString,*) ipar(j)
             END DO
             CALL Info('MAIN','Read '//I2S(ni)//' integer parameters from command line!')
             CALL SetIntegerParametersMATC(ni,ipar)
           END IF

           Silent = Silent .OR. &
               ( OptionString=='-s' .OR. OptionString=='--silent' ) 
           Version = Version .OR. &
               ( OptionString=='-v' .OR. OptionString == '--version' )
         END DO
       END IF

       ! Set number of OpenMP threads
       nthreads = 1
       !$ nthreads = omp_get_max_threads()
       IF (nthreads > 1 ) THEN
         ! Check if OMP_NUM_THREADS environment variable is set
         CALL envir( 'OMP_NUM_THREADS', threads, tlen )

         IF (tlen==0) THEN
           CALL Info('MAIN','OMP_NUM_THREADS not set. Using only 1 thread per task.',Level=6)
           nthreads = 1
           ! Set number of threads to 1
           !$ CALL omp_set_num_threads(nthreads)
#ifdef HAVE_MKL
           CALL mkl_set_num_threads(nthreads)
#endif
         END IF
       END IF
       ParEnv % NumberOfThreads = nthreads

       
       IF( .NOT. Silent ) THEN
         CALL Info( 'MAIN', ' ')
         CALL Info( 'MAIN', '=============================================================')
         CALL Info( 'MAIN', 'ElmerSolver finite element software, Welcome!                ')
         CALL Info( 'MAIN', 'This program is free software licensed under (L)GPL          ')
         CALL Info( 'MAIN', 'Copyright 1st April 1995 - , CSC - IT Center for Science Ltd.')
         CALL Info( 'MAIN', 'Webpage http://www.csc.fi/elmer, Email elmeradm@csc.fi       ')
         CALL Info( 'MAIN', 'Version: ' // GetVersion() // ' (Rev: ' // GetRevision() // &
                            ', Compiled: ' // GetCompilationDate() // ')' )

         IF ( ParEnv % PEs > 1 ) THEN
           CALL Info( 'MAIN', ' Running in parallel using ' // &
               i2s(ParEnv % PEs) // ' tasks.')
         ELSE
           CALL Info('MAIN', ' Running one task without MPI parallelization.',Level=10)
         END IF
         
         ! Print out number of threads in use
         IF ( nthreads > 1 ) THEN
           CALL Info('MAIN', ' Running in parallel with ' // &
               i2s(nthreads) // ' threads per task.')
         ELSE
           CALL Info('MAIN', ' Running with just one thread per task.',Level=10)
         END IF
         
#ifdef HAVE_FETI4I
         CALL Info( 'MAIN', ' FETI4I library linked in.')
#endif
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
#ifdef HAVE_MMG
         CALL Info( 'MAIN', ' MMG library linked in.')
#endif
#ifdef HAVE_PARMMG
         CALL Info( 'MAIN', ' ParMMG library linked in.')
#endif         
#ifdef HAVE_MKL
         CALL Info( 'MAIN', ' Intel MKL linked in.' )
#endif
#ifdef HAVE_LUA
         CALL Info( 'MAIN', ' Lua interpreter linked in.' )
#endif
#ifdef HAVE_EXTOPTIM
         CALL Info( 'MAIN', ' External optimization routines linked in.' )
#endif
#ifdef HAVE_ZOLTAN
         CALL Info( 'MAIN', ' Zoltan library linked in.' )
#endif
#ifdef HAVE_AMGX
         CALL Info( 'MAIN', ' AMGX library linked in.' )
#endif
#ifdef HAVE_ROCALUTION
         CALL Info( 'MAIN', ' ROCALUTION library linked in.' )
#endif
         CALL Info( 'MAIN', '=============================================================')
       END IF

       IF( Version ) RETURN
       
       CALL InitializeElementDescriptions()
     END IF

     ! Read input file name either as an argument, or from the default file:
     !----------------------------------------------------------------------
     GotModelName = .FALSE.
     IF ( NoArgs > 0 ) THEN
       CALL GET_COMMAND_ARGUMENT(1, ModelName)
       IF( ModelName(1:1) /= '-') THEN 
         GotModelName = .TRUE.
         IF (NoArgs > 1) CALL GET_COMMAND_ARGUMENT(2, eq)
       END IF
     END IF

     IF( .NOT. GotModelName ) THEN
       OPEN( 1, File='ELMERSOLVER_STARTINFO', STATUS='OLD', IOSTAT=iostat )       
       IF( iostat /= 0 ) THEN
         CALL Fatal( 'MAIN', 'Unable to find ELMERSOLVER_STARTINFO, can not execute.' )
       END IF
       READ(1,'(a)') ModelName
       CLOSE(1)
     END IF

     ! This sets optionally some internal parameters for doing scanning
     ! over a parameter space / optimization. 
     !-----------------------------------------------------------------
     IF( FirstTime ) THEN
       OPEN( Unit=InFileUnit, Action='Read',File=ModelName,Status='OLD',IOSTAT=iostat)         
       IF( iostat /= 0 ) THEN
         CALL Fatal( 'MAIN', 'Unable to find input file [' // &
             TRIM(Modelname) // '], can not execute.' )
       END IF
       ALLOCATE( Control )          
       ! Read only the "Run Control" section of the sif file.
       CALL LoadInputFile( Control,InFileUnit,ModelName,MeshDir,MeshName, &
           .FALSE., .TRUE., ControlOnly = .TRUE.)
       DoControl =  ASSOCIATED( Control % Control )

       IF( DoControl ) THEN
         CALL Info('MAIN','Run Control section active!')
         OptimIters = ListGetInteger( Control % Control,'Run Control Iterations', GotOptimIters )
         IF(.NOT. GotOptimIters) OptimIters = 1              
         
         ! If there are no parameters this does nothing
         CALL ControlParameters(Control % Control,1,GotParams,FinishEarly)
       ELSE
         OptimIters = 1 
       END IF
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
         FirstLoad = .TRUE.
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

         ! Here we read the whole model including command file and default mesh file
         !---------------------------------------------------------------------------------
         OPEN( Unit=InFileUnit, Action='Read',File=ModelName,Status='OLD',IOSTAT=iostat)
         IF( iostat /= 0 ) THEN
           CALL Fatal( 'MAIN', 'Unable to find input file [' // &
               TRIM(Modelname) // '], can not execute.' )
         END IF
         
         CurrentModel => LoadModel(ModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE,MeshIndex)
         IF(.NOT.ASSOCIATED(CurrentModel)) EXIT


         IF( nthreads > 1 ) THEN
           MaxOutputThread = ListGetInteger( CurrentModel % Simulation,'Max Output Thread',GotIt)
           IF(.NOT. GotIt) MaxOutputThread = 1
         END IF
         
         !----------------------------------------------------------------------------------
         ! Set namespace searching mode
         !----------------------------------------------------------------------------------
         CALL SetNamespaceCheck( ListGetLogical( CurrentModel % Simulation, &
                    'Additive namespaces', Found ) )

         !----------------------------------------------------------------------------------
         MeshMode = ListGetLogical( CurrentModel % Simulation, 'Mesh Mode', Found)

         !------------------------------------------------------------------------------
         ! Some keywords automatically require other keywords to be set
         ! We could complain on the missing keywords later on, but sometimes 
         ! it may be just as simple to add them directly. 
         !------------------------------------------------------------------------------
         CALL CompleteModelKeywords( )
         
         !----------------------------------------------------------------------------------
         ! Optionally perform simple extrusion to increase the dimension of the mesh
         !----------------------------------------------------------------------------------
         CALL CreateExtrudedMesh() 

         !----------------------------------------------------------------------------------
         ! If requested perform coordinate transformation directly after is has been obtained.
         ! Don't maintain the original mesh. 
         !----------------------------------------------------------------------------------
         CoordTransform=ListGetString(CurrentModel % Simulation,'Coordinate Transformation',GotIt)
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
               File=ModelName,Status='OLD',IOSTAT=iostat)
           IF( iostat /= 0 ) THEN
             CALL Fatal( 'MAIN', 'Unable to find input file [' // &
                 TRIM(Modelname) // '], can not execute.' )
           END IF                               
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
       Scanning  = ( eq == 'scanning' )
       Transient = ( eq == 'transient' )

!------------------------------------------------------------------------------
!      To more conveniently support the use of VTK based visualization there 
!      is a hack that recognizes VTU suffix and creates a instance of output solver.
!      Note that this is really quite a dirty hack, and is not a good example.
!-----------------------------------------------------------------------------
       CALL AddVtuOutputSolverHack()

! Support easily saving scalars to a file activated by "Scalars File" keyword.
!-----------------------------------------------------------------------------
       CALL AddSaveScalarsHack()


!------------------------------------------------------------------------------
!      Add coordinates such that if there is a solver that is run on creation
!      the coordinates are already usable then.
!------------------------------------------------------------------------------
       IF ( FirstLoad ) CALL AddMeshCoordinates()

!------------------------------------------------------------------------------
!      Figure out what (flow,heat,stress,...) should be computed, and get
!      memory for the dofs
!------------------------------------------------------------------------------
       CALL AddSolvers()

!------------------------------------------------------------------------------
!      Time integration and/or steady state steps
!------------------------------------------------------------------------------
       CALL InitializeIntervals()

!------------------------------------------------------------------------------
!      Add time and other global variables so that we can have dependence on these.
!------------------------------------------------------------------------------
       IF ( FirstLoad ) CALL AddTimeEtc()

!------------------------------------------------------------------------------
!      Initialize the random seeds so that all simulation depending on that
!      give consistent results.
!------------------------------------------------------------------------------      
       IF( FirstLoad ) CALL InitializeRandomSeed()
       
!------------------------------------------------------------------------------
!      Get Output File Options
!------------------------------------------------------------------------------
       
       ! Initial Conditions:
       ! -------------------
       IF ( FirstLoad ) THEN
         CALL SetInitialConditions()     

         DO i=1,CurrentModel % NumberOfSolvers
           Solver => CurrentModel % Solvers(i)
           IF( ListGetLogical( Solver % Values, 'Initialize Exported Variables', GotIt ) ) THEN
             CurrentModel % Solver => Solver
             CALL UpdateExportedVariables( Solver )	 
           END IF
         END DO
       END IF
         
       ! Compute the total number of steps that will be saved to the files
       ! Particularly look if the last step will be saved, or if it has
       ! to be saved separately.
       !------------------------------------------------------------------
       CALL CountSavedTimesteps() 
       
       CALL ListAddLogical( CurrentModel % Simulation,  &
            'Initialization Phase', .FALSE. )

       FirstLoad = .FALSE.
       IF ( Initialize == 1 ) EXIT

       ! Check whether we are using external optimization routine that
       ! needs to have basically Elmer given as a function that returns
       ! the cost function. Hence this is treated separately from the internal
       ! optimization methods.
       !---------------------------------------------------------------------
       ExtMethod = 0
       str = ListGetString( CurrentModel % Control,'Optimization Method',Found)       
       IF( Found ) THEN
         IF( SEQL(str, 'hybrd') ) THEN
           ExtMethod = 1
         ELSE IF( SEQL(str,'newuoa') ) THEN
           ExtMethod = 2
         ELSE IF( SEQL(str,'bobyqa') ) THEN 
           ExtMethod = 3
         END IF
       END IF

       ExecCommand = ListGetString( CurrentModel % Simulation, &
           'Control Procedure', ProcControl )       
       IF ( ProcControl ) THEN
         ControlProcedure = GetProcAddr( ExecCommand )
         CALL ExecSimulationProc( ControlProcedure, CurrentModel )
         
       ELSE IF( ExtMethod > 0 ) THEN
#ifdef HAVE_EXTOPTIM
         SELECT CASE( ExtMethod ) 
         CASE(1)
           CALL ExternalOptimization_minpack(ExecSimulationFunVec)           
         CASE(2)
           CALL ExternalOptimization_newuoa(ExecSimulationFunCost)                    
         CASE(3)
           CALL ExternalOptimization_bobyqa(ExecSimulationFunCost)                    
         END SELECT
#else
         CALL Fatal('MAIN','Compile WITH_EXTOPTIM to activate method: '//TRIM(str))
#endif
           
       ELSE IF( DoControl ) THEN

         ! This sets optionally some internal parameters for doing scanning
         ! over a parameter space / optimization. 
         !-----------------------------------------------------------------
         iSweep = 0
         DO WHILE (.TRUE.) 
           iSweep = iSweep + 1
           CALL Info('MAIN','========================================================',Level=5)
           CALL Info('MAIN','Control Loop '//I2S(iSweep))
           CALL Info('MAIN','========================================================',Level=5)
           
           sSweep = 1.0_dp * iSweep
           ! If there are no parameters this does nothing                  
           CALL ControlResetMesh(CurrentModel % Control, iSweep )            
           IF( iSweep > 1 ) THEN
             CALL ControlParameters(CurrentModel % Control,iSweep,&
                 GotParams,FinishEarly)           
             IF( FinishEarly ) EXIT
             Found = ReloadInputFile(CurrentModel,RewindFile=.TRUE.)
             CALL InitializeIntervals()
           END IF
           
           ! This is another calling slot as here we have formed the model structure and
           ! may toggle with the keyword coefficients. 
           CALL ControlParameters(CurrentModel % Control,iSweep,&
               GotParams,FinishEarly,SetCoeffs=.TRUE.)
           
           IF( iSweep > 1 ) THEN
             IF( ListGetLogical( CurrentModel % Control,'Reset Adaptive Mesh',Found ) ) THEN
               CALL ResetAdaptiveMesh()
             END IF
             IF( ListGetLogical( CurrentModel % Control,'Reset Initial Conditions',Found ) ) THEN
               CALL SetInitialConditions()
             END IF
           END IF

           CALL ExecSimulation( TimeIntervals, CoupledMinIter, &
               CoupledMaxIter, OutputIntervals, Transient, Scanning) 
           
           ! This evaluates the cost function and saves the results of control
           CALL ControlParameters(CurrentModel % Control, &
               iSweep,GotParams,FinishEarly,.TRUE.)

           IF( iSweep == 1 ) THEN
             DO i=1,CurrentModel % NumberOfSolvers 
               iSolver => CurrentModel % Solvers(i)
               j = iSolver % NumberOfConstraintModes
               IF( j <= 0 ) CYCLE
               IF( ListGetLogical( iSolver % Values,'Run Control Constraint Modes', Found ) .OR. &
                   ListGetLogical( CurrentModel % Control,'Constraint Modes Analysis',Found ) ) THEN
                 IF( GotOptimIters ) THEN
                   IF( OptimIters /= j ) THEN
                     CALL Warn('MAIN','Incompatible number of run control iterations and constraint modes!')
                   END IF
                 ELSE
                   CALL Info('MAIN','Setting run control iterations to constraint modes count: '//I2S(j))
                   OptimIters = j
                 END IF
                 EXIT
               END IF
             END DO
           END IF
           
           ! We use this type of condition so that OptimIters can be changed on-the-fly
           IF(iSweep == OptimIters) EXIT
         END DO

         DO i=1,CurrentModel % NumberOfSolvers 
           iSolver => CurrentModel % Solvers(i)
           IF( iSolver % NumberOfConstraintModes > 0 ) THEN
             IF( ListGetLogical( iSolver % Values,'Run Control Constraint Modes', Found ) .OR. &
                 ListGetLogical( CurrentModel % Control,'Constraint Modes Analysis',Found ) ) THEN
               CALL FinalizeLumpedMatrix( iSolver )            
             END IF
           END IF
         END DO

         DO i=1,CurrentModel % NumberOfSolvers 
           iSolver => CurrentModel % Solvers(i)
           IF ( iSolver % PROCEDURE == 0 ) CYCLE
           When = ListGetString( iSolver % Values, 'Exec Solver', Found )
           IF ( Found ) THEN
             DoIt = ( When == 'after control' ) 
           ELSE
             DoIt = ( iSolver % SolverExecWhen == SOLVER_EXEC_AFTER_CONTROL )
           END IF
           IF(DoIt) CALL SolverActivate( CurrentModel,iSolver,dt,Transient )
         END DO
       ELSE
         CALL ExecSimulation( TimeIntervals, CoupledMinIter, &
             CoupledMaxIter, OutputIntervals, Transient, Scanning) 
       END IF
       
       ! Comparison to reference is done to enable consistency test under ctest.
       !-------------------------------------------------------------------------
       CALL CompareToReferenceSolution( )

       IF ( Initialize >= 2 ) EXIT
     END DO

     IF( ListGetLogical( CurrentModel % Simulation,'Echo Keywords at End', GotIt ) ) THEN
       CALL ListEchoKeywords( CurrentModel )        
     END IF
      
     CALL CompareToReferenceSolution( Finalize = .TRUE. )
     
#ifdef DEVEL_LISTUSAGE
     IF(InfoActive(10) .AND. ParEnv % MyPe == 0 ) THEN
       CALL Info('MAIN','Reporting unused list entries for sif improvement!')
       CALL Info('MAIN','If you do not want these lines undefine > DEVEL_LISTUSAGE < !')       
       CALL ReportListCounters( CurrentModel, 2 )
     END IF
#endif
#ifdef DEVEL_LISTCOUNTER
     IF( ParEnv % MyPe == 0 ) THEN
       CALL Info('MAIN','Reporting list counters for code optimization purposes only!')
       CALL Info('MAIN','If you get these lines with production code undefine > DEVEL_LISTCOUNTER < !')
       CALL ReportListCounters( CurrentModel, 3 )
     END IF
#endif
          
!------------------------------------------------------------------------------
!    THIS IS THE END (...,at last, the end, my friend,...)
!------------------------------------------------------------------------------
     IF ( Initialize /= 1 ) CALL Info( 'MAIN', '*** Elmer Solver: ALL DONE ***',Level=3 )

     ! This may be used to study problems at the finish
     IF( ListGetLogical( CurrentModel % Simulation,'Dirty Finish', GotIt ) ) THEN
       CALL Info('MAIN','Skipping freeing of the Model structure',Level=4)
       RETURN
     END IF
     
     IF ( Initialize <= 0 ) CALL FreeModel(CurrentModel)

#ifdef HAVE_TRILINOS
  CALL TrilinosCleanup()
#endif

     IF ( FirstTime ) CALL ParallelFinalize()
     FirstTime = .FALSE.

     CALL Info('MAIN','The end',Level=3)

     RETURN

   CONTAINS 


     ! If we want to start a new adaptive simulation with the original mesh
     ! call this subroutine.
     !---------------------------------------------------------------------
     SUBROUTINE ResetAdaptiveMesh()

       TYPE(Mesh_t), POINTER :: pMesh, pMesh0
       TYPE(Solver_t), POINTER :: iSolver
       TYPE(Variable_t), POINTER :: pVar
       LOGICAL :: GB, BO
       CHARACTER(*), PARAMETER :: Caller = 'ResetAdaptiveMesh'
       
       ! Find the 1st mesh
       pMesh0 => CurrentModel % Mesh 
       DO WHILE( ASSOCIATED(pMesh0 % Parent) )
         pMesh0 => pMesh0 % Parent
       END DO
       !PRINT *,'First mesh:',pMesh0 % AdaptiveDepth, TRIM(pMesh0 % Name)

       ! Find the last mesh
       pMesh => CurrentModel % Mesh 
       DO WHILE( ASSOCIATED(pMesh % Child) )
         pMesh => pMesh % Child
       END DO
       !PRINT *,'Last mesh:',pMesh % AdaptiveDepth, TRIM(pMesh % Name)

       ! Move point to the 1st mesh and related fields
       CALL SetCurrentMesh( CurrentModel, pMesh0 )

       DO i=1,CurrentModel % NumberOfSolvers 
         iSolver => CurrentModel % Solvers(i)

         IF(.NOT. ASSOCIATED(iSolver % Variable)) THEN
           CALL Info(Caller,'No Variable in this mesh for solver index '//I2S(i),Level=10)
           CYCLE
         END IF

         ! Set Solver mesh 
         IF(ASSOCIATED(iSolver % Mesh)) iSolver % Mesh => pMesh0

         ! Set Solver variable point to the field in the original mesh
         IF(ASSOCIATED(iSolver % Variable)) THEN
           pVar => VariableGet(pMesh0 % Variables, &
               iSolver % Variable % Name, ThisOnly = .TRUE.)  
           IF(.NOT. ASSOCIATED(pVar)) THEN
             CALL Info(Caller,'No Variable in coarsest mesh for solver index '//I2S(i),Level=10)
             CYCLE
           END IF
           iSolver % Variable => pVar
         END IF

         ! Reset active element table
         IF( iSolver % NumberOfActiveElements == 0 ) THEN
           CALL Info(Caller,'No active elements for solver index '//I2S(i),Level=10)
           CYCLE
         END IF
           
         iSolver % NumberOfActiveElements = 0
         CALL SetActiveElementsTable( CurrentModel, iSolver )                   

         ! Create the matrix related to the original mesh 
         IF( .NOT. ASSOCIATED( iSolver % Matrix ) ) THEN
           CALL Info(Caller,'No matrix for solver index '//I2S(i),Level=10)
           CYCLE
         END IF
           
         CALL FreeMatrix( iSolver % Matrix)

         GB = ListGetLogical( iSolver % Values,'Bubbles in Global System', Found )
         IF ( .NOT. Found ) GB = .TRUE.

         BO = ListGetLogical( iSolver % Values,'Optimize Bandwidth', Found )
         IF ( .NOT. Found ) BO = .TRUE.

         iSolver % Matrix => CreateMatrix( CurrentModel, iSolver, iSolver % Mesh,  &
             iSolver % Variable % Perm, iSolver % Variable % DOFs, MATRIX_CRS, &
             BO, ListGetString( iSolver % Values, 'Equation' ), GlobalBubbles=GB )
         ALLOCATE( iSolver % Matrix % rhs(iSolver % Matrix % NumberOfRows ) )
         iSolver % Matrix % rhs = 0.0_dp
       END DO

       ! Release the old adaptive meshes
       DO WHILE( ASSOCIATED(pMesh % Parent))
         pMesh => pMesh % Parent             
         CALL ReleaseMesh( pMesh % Child ) 
       END DO
       pMesh % Child => NULL()
                 
     END SUBROUTINE ResetAdaptiveMesh


     
     SUBROUTINE InitializeRandomSeed()
       INTEGER :: i,n
       INTEGER, ALLOCATABLE :: seeds(:)
       
       CALL RANDOM_SEED() ! initialize with system generated seed

       i = ListGetInteger( CurrentModel % Simulation,'Random Number Seed',Found ) 
       IF( .NOT. Found ) i = 314159265

       CALL RANDOM_SEED(size=j) ! find out size of seed
       ALLOCATE(seeds(j))
       !CALL RANDOM_SEED(get=seeds) ! get system generated seed
       !WRITE(*,*) seeds            ! writes system generated seed
       seeds = i
       CALL RANDOM_SEED(put=seeds) ! set current seed
       !CALL RANDOM_SEED(get=seeds) ! get current seed
       !WRITE(*,*) seeds            ! writes 314159265
       DEALLOCATE(seeds)           
       
       CALL Info('MAIN','Random seed initialized to: '//I2S(i),Level=10)
     END SUBROUTINE InitializeRandomSeed

     
     ! Optionally create extruded mesh on-the-fly.
     !--------------------------------------------------------------------
     SUBROUTINE CreateExtrudedMesh()

       LOGICAL :: SliceVersion

       IF(.NOT. ListCheckPrefix(CurrentModel % Simulation,'Extruded Mesh') ) RETURN
       
       SliceVersion = GetLogical(CurrentModel % Simulation,'Extruded Mesh Slices',Found )              
       IF( SliceVersion ) THEN
         ExtrudedMesh => MeshExtrudeSlices(CurrentModel % Meshes, CurrentModel % Simulation )
       ELSE
         ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, CurrentModel % Simulation)
       END IF
         
       ! Make the solvers point to the extruded mesh, not the original mesh
       !-------------------------------------------------------------------
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

     END SUBROUTINE CreateExtrudedMesh
     

     ! Given geometric ratio of timesteps redistribute them so that the ratio
     ! is met as closely as possible, maintaining total time and sacrificing
     ! number of timesteps. 
     !----------------------------------------------------------------------
     SUBROUTINE GeometricTimesteps(m,n0,dt0,r)
       INTEGER :: m
       INTEGER :: n0(:)
       REAL(KIND=dp) :: dt0(:),r(:)

       INTEGER :: i,n
       REAL(KIND=dp) :: q
       LOGICAL :: Visited = .FALSE.

       ! Only do this once since it tampers stuff in lists.
       IF(Visited) RETURN
       Visited = .TRUE.

       CALL Info('MAIN','Creating geometric timestepping strategy',Level=6)
       
       DO i=1,m
         ! Some users may give zero ratio, assume that they mean one.
         IF(ABS(r(i)) < EPSILON(q) ) r(i) = 1.0_dp
         ! Ratio one means even distribution.
         IF(ABS(r(i)-1.0) < EPSILON(q) ) CYCLE

         q = 1 + (r(i)-1)/n0(i)
         n = NINT( LOG(1+(q-1)*n0(i)) / LOG(q) )
         dt = n0(i)*dt0(i)*(1-q)/(1-q**n)
         
         !PRINT *,'ratio:',i,n0(i),dt0(i),r(i),n,dt,q

         ! Replace the new distribution
         r(i) = q
         dt0(i) = dt
         n0(i) = n
       END DO
       
     END SUBROUTINE GeometricTimesteps

     
     ! Initialize intervals for steady state, transient and scanning types.
     !----------------------------------------------------------------------
     SUBROUTINE InitializeIntervals()

       IF ( Transient .OR. Scanning ) THEN
         Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
             'Timestep Intervals', GotIt )

         IF ( .NOT.GotIt ) THEN
           CALL Fatal('MAIN', 'Keyword > Timestep Intervals < MUST be ' //  &
               'defined for transient and scanning simulations' )
         END IF
         
         TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation, &
             'Timestep Sizes', GotIt )
         IF ( .NOT.GotIt ) THEN
           IF( Scanning .OR. ListCheckPresent( CurrentModel % Simulation,'Timestep Size') ) THEN
             ALLOCATE(TimestepSizes(SIZE(Timesteps),1))
             TimestepSizes = 1.0_dp
           ELSE
             CALL Fatal( 'MAIN', 'Keyword [Timestep Sizes] MUST be ' //  &
                 'defined for time dependent simulations' )
           END IF
         END IF
         
         CoupledMaxIter = ListGetInteger( CurrentModel % Simulation, &
             'Steady State Max Iterations', GotIt, minv=1 )
         IF ( .NOT. GotIt ) CoupledMaxIter = 1
         
         TimeIntervals = SIZE(Timesteps)

         TimestepRatios => ListGetConstRealArray( CurrentModel % Simulation, &
             'Timestep Ratios', GotTimestepRatios )

         
         IF ( GotTimestepRatios ) THEN           
           CALL GeometricTimesteps(TimeIntervals,Timesteps,TimestepSizes(:,1),TimestepRatios(:,1))
         END IF
         

       ELSE
         ! Steady state
         !------------------------------------------------------------------------------
         IF( .NOT. ASSOCIATED(Timesteps) ) THEN
           ALLOCATE(Timesteps(1))
         END IF
         Timesteps(1) = ListGetInteger( CurrentModel % Simulation, &
             'Steady State Max Iterations', GotIt,minv=1 )
         IF ( .NOT. GotIt ) Timesteps(1) = 1
         CoupledMaxIter = 1 
         
         ALLOCATE(TimestepSizes(1,1))
         TimestepSizes(1,1) = 1.0_dp
         TimeIntervals  = 1
       END IF

      
       CoupledMinIter = ListGetInteger( CurrentModel % Simulation, &
           'Steady State Min Iterations', GotIt )
       IF( .NOT. GotIt ) CoupledMinIter = 1 
       
       OutputIntervals => ListGetIntegerArray( CurrentModel % Simulation, &
           'Output Intervals', GotIt )
       IF( GotIt ) THEN
         IF( SIZE(OutputIntervals) /= SIZE(TimeSteps) ) THEN
           CALL Fatal('MAIN','> Output Intervals < should have the same size as > Timestep Intervals < !')
         END IF
       ELSE 
         IF( .NOT. ASSOCIATED( OutputIntervals ) ) THEN
           ALLOCATE( OutputIntervals(SIZE(TimeSteps)) )
           OutputIntervals = 1
         END IF
       END IF 


       IF ( FirstLoad ) &
           ALLOCATE( sTime(1), sStep(1), sInterval(1), sSize(1), &
           steadyIt(1), nonLinit(1), sPrevSizes(1,5), sPeriodicTime(1), &
           sPeriodicCycle(1), sPar(1), sScan(1), sSweep(1), sFinish(1), &
           sProduce(1),sSlice(1), sSector(1), sSliceRatio(1), sSliceWeight(1), sAngle(1), &
           sAngleVelo(1) )
       
       dt = 0._dp       
       sTime = 0._dp
       sStep = 0
       sPeriodicTime = 0._dp
       sPeriodicCycle = 0._dp
       sScan = 0._dp       
       sSize = dt
       sPrevSizes = 0_dp       
       sInterval = 0._dp       
       steadyIt = 0
       nonlinIt = 0
       sPar = 0
       sFinish = -1.0_dp
       sProduce = -1.0_dp
       sSlice = 0._dp
       sSector = 0._dp
       sSliceRatio = 0._dp
       sSliceWeight = 1.0_dp
       sAngle = 0.0_dp
       sAngleVelo = 0.0_dp
       
     END SUBROUTINE InitializeIntervals
       

     ! Calculate the number of timesteps for saving.
     ! This is needed for some legacy formats.
     !-----------------------------------------------------------------
     SUBROUTINE CountSavedTimesteps()

       INTEGER :: i
       
       TotalTimesteps = 0
       LastSaved = .FALSE.
       DO i=1,TimeIntervals
         DO timestep = 1,Timesteps(i)
           IF( OutputIntervals(i) == 0 ) CYCLE
           LastSaved = .FALSE.
           IF ( MOD(Timestep-1, OutputIntervals(i))==0 ) THEN
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
       
       CALL Info('MAIN','Number of timesteps to be saved: '//I2S(TotalTimesteps))
       
     END SUBROUTINE CountSavedTimesteps
     
     
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
       CHARACTER(:), ALLOCATABLE :: PassedMsg

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
               'PASSED all '//I2S(TestCount)//' tests!',Level=3)
         ELSE         
           CALL Warn('CompareToReferenceSolution','FAILED '//I2S(FailCount)//&
               ' tests out of '//I2S(TestCount)//'!')
         END IF
         
         IF( FinalizeOnly ) THEN
           IF( ParEnv % MyPe == 0 ) THEN
             IF( ParEnv % PEs > 1 ) THEN
               ! Parallel test, add the number of tasks as a suffix
               PassedMsg = "TEST.PASSED_"//I2S(ParEnv % PEs)
               OPEN( 10, FILE = PassedMsg )
             ELSE
               OPEN( 10, FILE = 'TEST.PASSED' )
             END IF
             IF( Success ) THEN
               WRITE( 10,'(I1)' ) 1
             ELSE
               WRITE( 10,'(I1)' ) 0
             END IF
             CALL FLUSH( 10 )
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
               //I2S(i)//' not associated, cannot compare')
           CYCLE
         END IF

         TestCount = TestCount + 1
         Success = .TRUE.

         ! Compare either to existing norm (ensures consistency)
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
           IF( Err > Tol .OR. Err /= Err ) THEN
             ! Warn only in the main core
             IF( ParEnv % MyPe == 0 ) THEN
               WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                   'Solver ',solver_id,' FAILED:  Norm =',Norm,'  RefNorm =',RefNorm
               CALL Warn('CompareToReferenceSolution',Message)
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                 'Solver ',solver_id,' PASSED:  Norm =',Norm,'  RefNorm =',RefNorm
             CALL Info('CompareToReferenceSolution',Message,Level=3)
           END IF
           IF( AbsoluteErr ) THEN
             WRITE( Message,'(A,ES13.6)') 'Absolute Error to reference norm:',Err
           ELSE
             WRITE( Message,'(A,ES13.6)') 'Relative Error to reference norm:',Err
           END IF
           CALL Info('CompareToReferenceSolution',Message, Level = 3 )
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
                 refval = ListGetRealAtNode( Solver % Values,'Reference Solution '//I2S(k),i,Found ) 
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
               CALL Info('CompareToReferenceSolution',Message, Level = 3 )
             END IF
             Success = .FALSE.
           ELSE         
             WRITE( Message,'(A,I0,A,ES15.8,A,ES15.8)') &
                 'Solver ',solver_id,' PASSED:  Solution =',Norm,'  RefSolution =',RefNorm
             CALL Info('CompareToReferenceSolution',Message,Level=3)
           END IF
         END IF

         IF( Success ) PassCount = PassCount + 1
       END DO
      
     END SUBROUTINE CompareToReferenceSolution


     SUBROUTINE AppendNewSolver(Model,pSolver)
       TYPE(Model_t) :: Model
       TYPE(Solver_t), POINTER :: pSolver
       
       TYPE(Solver_t), POINTER :: OldSolvers(:),NewSolvers(:)
       INTEGER :: i, j,j2,j3,n, AllocStat
       
       n = Model % NumberOfSolvers+1
       ALLOCATE( NewSolvers(n), STAT = AllocStat )
       IF( AllocStat /= 0 ) CALL Fatal('AppendNewSolver','Allocation error 1')

       OldSolvers => Model % Solvers
       
       CALL Info('AppendNewSolver','Increasing number of solvers to: '&
           //I2S(n),Level=8)
       DO i=1,n-1
         ! Def_Dofs is the only allocatable structure within Solver_t:
         IF( ALLOCATED( OldSolvers(i) % Def_Dofs ) ) THEN
           j = SIZE(OldSolvers(i) % Def_Dofs,1)
           j2 = SIZE(OldSolvers(i) % Def_Dofs,2)
           j3 = SIZE(OldSolvers(i) % Def_Dofs,3)
           ALLOCATE( NewSolvers(i) % Def_Dofs(j,j2,j3), STAT = AllocStat )
           IF( AllocStat /= 0 ) CALL Fatal('AppendNewSolver','Allocation error 2')           
         END IF

         ! Copy the content of the Solver structure
         NewSolvers(i) = OldSolvers(i)

         ! Nullify the old structure since otherwise bad things may happen at deallocation
         NULLIFY( OldSolvers(i) % ActiveElements )
         NULLIFY( OldSolvers(i) % Mesh )
         NULLIFY( OldSolvers(i) % BlockMatrix )
         NULLIFY( OldSolvers(i) % Matrix )
         NULLIFY( OldSolvers(i) % Variable )
       END DO

       ! Deallocate the old structure and set the pointer to the new one
       DEALLOCATE( Model % Solvers )
       Model % Solvers => NewSolvers
       Model % NumberOfSolvers = n

       pSolver => NewSolvers(n)

       NULLIFY( pSolver % Matrix )
       NULLIFY( pSolver % Mesh ) 
       NULLIFY( pSOlver % BlockMatrix )
       NULLIFY( pSolver % Variable )
       NULLIFY( pSolver % ActiveElements )
       
       pSolver % PROCEDURE = 0
       pSolver % NumberOfActiveElements = 0
       j = CurrentModel % NumberOfBodies
       ALLOCATE( pSolver % Def_Dofs(10,j,6),STAT=AllocStat)       
       IF( AllocStat /= 0 ) CALL Fatal('AppendNewSolver','Allocation error 3')
       pSolver % Def_Dofs = -1
       pSolver % Def_Dofs(:,:,1) =  1
       
       ! Create empty list to add some keywords to 
       pSolver % Values => ListAllocate()
       
     END SUBROUTINE AppendNewSolver
     

     !------------------------------------------------------------------------
     !> Given name of a solver module find it among the active solvers and
     !> upon success return its index. 
     !------------------------------------------------------------------------
     FUNCTION FindSolverByProcName(Model,ProcName) RESULT (solver_id)
       IMPLICIT NONE

       TYPE(Model_t), POINTER :: Model
       CHARACTER(*) :: ProcName
       INTEGER :: solver_id
       
       LOGICAL :: Found
       INTEGER :: i,j
       CHARACTER(:), ALLOCATABLE :: str
       TYPE(Solver_t), POINTER :: pSolver

       solver_id = 0       
       Found = .FALSE.

       DO i=1, Model % NumberOfSolvers
         pSolver => CurrentModel % Solvers(i)
         str = ListGetString(pSolver % Values,'Procedure',Found)
         IF(.NOT. Found) CYCLE

         j = INDEX(str,ProcName)         
         IF( j > 0 ) THEN
           solver_id = i
           EXIT
         END IF
       END DO
       
     END FUNCTION FindSolverByProcName
     !------------------------------------------------------------------------------


     
     ! This is a dirty hack that adds an instance of ResultOutputSolver to the list of Solvers.
     ! The idea is that it is much easier for the end user to take into use the vtu output this way.
     ! The solver itself has limited set of parameters needed and is therefore approapriate for this
     ! kind of hack. It can of course be also added as a regular solver also.
     !----------------------------------------------------------------------------------------------
     SUBROUTINE AddVtuOutputSolverHack()     
       TYPE(Solver_t), POINTER :: pSolver
       CHARACTER(:), ALLOCATABLE :: str
       INTEGER :: j,k
       TYPE(ValueList_t), POINTER :: Params, Simu
       LOGICAL :: Found, VtuFormat
       INTEGER :: AllocStat
       LOGICAL, SAVE :: Visited = .FALSE.
       
       Simu => CurrentModel % Simulation
       str = ListGetString( Simu,'Post File',Found) 
       IF(.NOT. Found) RETURN
       
       k = INDEX( str,'.vtu' )
       VtuFormat = ( k /= 0 ) 
       IF(.NOT. VtuFormat ) RETURN

       ! No use to create the same solver twice
       IF( Visited ) RETURN
       Visited = .TRUE.
       
       CALL Info('AddVtuOutputSolverHack','Adding ResultOutputSolver to write VTU output in file: '&
           //TRIM(str(1:k-1)))

       j = FindSolverByProcName(CurrentModel,'ResultOutputSolver')
       IF(j>0) THEN
         pSolver => CurrentModel % Solvers(j)
         IF( ListGetLogical( pSolver % Values,'Vtu Format',Found) .OR. &
             ListGetString(pSolver % Values,'Output Format',Found) == 'vtu' ) THEN
           CALL Warn('AddVtuOutputSolverHack','ResultOutputSolver instance with VTU format already exists, doing nothing!')
           RETURN
         END IF
       END IF
              
       ! Remove the post file from the simulation list as it will be dealt by the solver section
       CALL ListRemove( Simu,'Post File')

       ! Allocate one new solver to the end of list and get pointer to it
       CALL AppendNewSolver(CurrentModel,pSolver)
       
       ! Now create the ResultOutputSolver instance on-the-fly
       Params => pSolver % Values
       CALL ListAddString(Params,'Procedure', 'ResultOutputSolve ResultOutputSolver',.FALSE.)
       CALL ListAddString(Params,'Equation','InternalVtuOutputSolver')
       CALL ListAddString(Params,'Output Format','vtu')
       CALL ListAddString(Params,'Output File Name',str(1:k-1),.FALSE.)
       CALL ListAddLogical(Params,'No Matrix',.TRUE.)
       CALL ListAddNewString(Params,'Variable','-global vtu_internal_dummy')
       CALL ListAddString(Params,'Exec Solver','after saving')
       CALL ListAddLogical(Params,'Save Geometry IDs',.TRUE.)
       
       ! Add a few often needed keywords also if they are given in simulation section
       CALL ListCopyPrefixedKeywords( Simu, Params, 'vtu:' )

       ! It makes sense to inherit global ascii/binary flags if not given
       CALL ListCompareAndCopy(CurrentModel % Simulation, Params, 'ascii output', nooverwrite = .TRUE. )
       CALL ListCompareAndCopy(CurrentModel % Simulation, Params, 'binary output', nooverwrite = .TRUE. )
       
       CALL Info('AddVtuOutputSolverHack','Finished appending VTU output solver',Level=12)
       
     END SUBROUTINE AddVtuOutputSolverHack


     ! This is a dirty hack that adds an instance of SaveScalars to the list of Solvers.
     ! The idea is that it is much easier for the end user to add a basic instance.
     !----------------------------------------------------------------------------------------------
     SUBROUTINE AddSaveScalarsHack()     
       TYPE(Solver_t), POINTER :: ABC(:), PSolver
       CHARACTER(:), ALLOCATABLE :: str
       INTEGER :: k
       TYPE(ValueList_t), POINTER :: Params, Simu
       LOGICAL :: Found, VtuFormat
       INTEGER :: AllocStat
       LOGICAL, SAVE :: Visited = .FALSE.
       
       Simu => CurrentModel % Simulation
       str = ListGetString( Simu,'Scalars File',Found) 
       IF(.NOT. Found) RETURN
       
       ! No use to create the same solver twice
       IF( Visited ) RETURN
       Visited = .TRUE.
       
       CALL Info('AddSaveScalarsHack','Adding SaveScalars solver to write scalars into file: '&
           //TRIM(str))

       k = FindSolverByProcName(CurrentModel,'SaveScalars')
       IF(k>0) THEN
         CALL Warn('AddSaveScalarsHack','SaveScalars instance already exists, doing nothing!')
         RETURN
       END IF
       
       ! Allocate one new solver to the end of list and get pointer to it
       CALL AppendNewSolver(CurrentModel,pSolver)
       
       ! Add some keywords to the list
       Params => pSolver % Values
       CALL ListAddString(Params,'Procedure', 'SaveData SaveScalars',.FALSE.)
       CALL ListAddString(Params,'Equation','InternalSaveScalars')
       CALL ListAddString(Params,'Filename',TRIM(str),.FALSE.)
       CALL ListAddString(Params,'Exec Solver','after saving')

       ! Add a few often needed keywords also if they are given in simulation section
       CALL ListCopyPrefixedKeywords( Simu, Params, 'scalars:' )

       CALL Info('AddSaveScalarsHack','Finished appending SaveScalars solver',Level=12)
       
     END SUBROUTINE AddSaveScalarsHack


     
!------------------------------------------------------------------------------
!> Adds flags for active solvers. 
!------------------------------------------------------------------------------
    SUBROUTINE AddSolvers()
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,n,nlen
      LOGICAL :: InitSolver, Found, DoTiming
      TYPE(Mesh_t), POINTER :: pMesh
!------------------------------------------------------------------------------

      CALL Info('AddSolvers','Setting up '//I2S(CurrentModel % NumberOfSolvers)//&
          ' solvers',Level=10)

      ! This is a hack that sets Equation flags True for the "Active Solvers".
      ! The Equation flag is the legacy way of setting a Solver active and is still
      ! used internally. Also set WhenExec flag since we might want to call some
      ! solvers immediately.
      !----------------------------------------------------------------------------
      DO i=1,CurrentModel % NumberOfSolvers
        Solver => CurrentModel % Solvers(i)
        
        eq = ListGetString( Solver % Values,'Equation', Found )     
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

        CALL AddExecWhenFlag( Solver ) 
     END DO


     ! Add the dynamically linked solver to be called later
     ! First do the initialization for solvers that other solvers except
     ! before initialization having the "when created" slot.  
     !---------------------------------------------------------------------
     DO i=1,CurrentModel % NumberOfSolvers
       Solver => CurrentModel % Solvers(i)
       IF ( Solver % SolverExecWhen /= SOLVER_EXEC_WHENCREATED ) CYCLE

       eq = ListGetString( Solver % Values,'Equation', Found )
       CALL Info('AddSolvers','Setting up solver '//I2S(i)//': '//TRIM(eq),Level=10)

       InitSolver = ListGetLogical( Solver % Values, 'Initialize', Found )
       IF ( Found .AND. InitSolver ) THEN
         CALL FreeMatrix( Solver % Matrix )
         CALL ListAddLogical( Solver % Values, 'Initialize', .FALSE. )
       END IF

       IF ( Solver % PROCEDURE == 0 .OR. InitSolver ) THEN
         IF ( .NOT. ASSOCIATED( Solver % Mesh ) ) THEN
           Solver % Mesh => CurrentModel % Meshes
         END IF

         n = ListGetInteger( Solver % Values,'Relative Mesh Level',Found )
         IF( Found .AND. n /= 0) THEN
           IF( n > 0 ) CALL Fatal('AddSolvers','Relative Mesh Level should ne negative!')
           j = 0
           DO WHILE(j > n)             
             pMesh => pMesh % Parent
             IF(.NOT. ASSOCIATED(pMesh)) EXIT
             j = j-1
           END DO
           IF(ASSOCIATED(pMesh) ) THEN
             CALL Info('AddSolvers','Using relative mesh level '//I2S(n),Level=10)
           ELSE
             CALL Fatal('AddSolvers','Could not find relative mesh level: '//I2S(n))
           END IF
           Solver % Mesh => pMesh
         END IF
         
         
         DoTiming = ListGetLogical( Solver % Values,'Solver Timing', Found ) 
         IF( DoTiming ) CALL ResetTimer('SolverInitialization')
         
         CurrentModel % Solver => Solver
         CALL AddEquationBasics( Solver, eq, Transient )
         CALL AddEquationSolution( Solver, Transient )

         IF( DoTiming ) CALL CheckTimer('SolverInitialization',Level=7,Delete=.TRUE.)
         
         CALL Info('AddSolvers','Executing solver '//I2S(i)//' immediately when created!,Level=5')
         CALL SetCurrentMesh( CurrentModel, Solver % Mesh )
         CALL SingleSolver( CurrentModel, Solver, 0.0_dp, .FALSE. )
       END IF
     END DO

     
     ! And now initialize the other solvers
     !---------------------------------------------------------------------
     DO i=1,CurrentModel % NumberOfSolvers
       Solver => CurrentModel % Solvers(i)
       IF ( Solver % SolverExecWhen == SOLVER_EXEC_WHENCREATED ) CYCLE

       eq = ListGetString( Solver % Values,'Equation', Found )
       CALL Info('AddSolvers','Setting up solver '//I2S(i)//': '//TRIM(eq),Level=10)
       
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

         DoTiming = ListGetLogical( Solver % Values,'Solver Timing', Found ) 
         IF( DoTiming ) CALL ResetTimer('SolverInitialization')
         
         CALL AddEquationBasics( Solver, eq, Transient )
         CALL AddEquationSolution( Solver, Transient )

         IF( DoTiming ) CALL CheckTimer('SolverInitialization',Level=7,Delete=.TRUE.)
       END IF
     END DO

     CALL Info('AddSolvers','Setting up solvers done',Level=12)

!------------------------------------------------------------------------------
  END SUBROUTINE AddSolvers
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Adds coordinate as variables to the current mesh structure. 
!------------------------------------------------------------------------------
  SUBROUTINE AddMeshCoordinates()
!------------------------------------------------------------------------------
     CALL Info('AddMeshCoordinates','Setting mesh coordinates and time',Level=10)

     Mesh => CurrentModel % Meshes 
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd( Mesh % Variables, Mesh, &
           Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )
       
       CALL VariableAdd(Mesh % Variables,Mesh, &
           Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )
       
       CALL VariableAdd(Mesh % Variables,Mesh, &
           Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )
       Mesh => Mesh % Next        
     END DO
     
!------------------------------------------------------------------------------
   END SUBROUTINE AddMeshCoordinates
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Adds coordinate and time variables to the current mesh structure. 
!------------------------------------------------------------------------------
  SUBROUTINE AddTimeEtc()
!------------------------------------------------------------------------------
     TYPE(Variable_t), POINTER :: DtVar
     
     CALL Info('AddTimeEtc','Setting time and other global variables',Level=10)

     NULLIFY( Solver )

     Mesh => CurrentModel % Meshes 
     DO WHILE( ASSOCIATED( Mesh ) )
       CALL VariableAdd( Mesh % Variables, Mesh, Name='Time',DOFs=1, Values=sTime )
       CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep', DOFs=1, Values=sStep )
       CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep size', DOFs=1, Values=sSize )
       CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep interval', DOFs=1, Values=sInterval )

       ! Save some previous timesteps for variable timestep multistep methods
       DtVar => VariableGet( Mesh % Variables, 'Timestep size' )
       DtVar % PrevValues => sPrevSizes

       CALL VariableAdd( Mesh % Variables, Mesh, &
               Name='nonlin iter', DOFs=1, Values=nonlinIt )
       CALL VariableAdd( Mesh % Variables, Mesh, &
               Name='coupled iter', DOFs=1, Values=steadyIt )


       IF( ListCheckPrefix( CurrentModel % Simulation,'Periodic Time') .OR. &
           ListCheckPresent( CurrentModel % Simulation,'Time Period') ) THEN
         ! For periodic systems we may do several cycles.
         CALL VariableAdd( Mesh % Variables, Mesh, Name='Periodic Time',DOFs=1, Values=sPeriodicTime )
         CALL VariableAdd( Mesh % Variables, Mesh, Name='Periodic Cycle',DOFs=1, Values=sPeriodicCycle )      
         
         ! After convergence is reached we may start producing the results.         
         CALL VariableAdd( Mesh % Variables, Mesh, Name='Finish',DOFs=1, Values=sFinish )
         CALL VariableAdd( Mesh % Variables, Mesh, Name='Produce',DOFs=1, Values=sProduce )         
       END IF

       IF( ListCheckPresent( CurrentModel % Simulation,'Rotor Angle') ) THEN
         CALL VariableAdd( Mesh % Variables, Mesh, Name='rotor angle',DOFs=1, Values=sAngle )
         CALL VariableAdd( Mesh % Variables, Mesh, Name='rotor velo',DOFs=1, Values=sAngleVelo )
       END IF
              
       IF( ListCheckPresentAnySolver( CurrentModel,'Scanning Loops') ) THEN
         CALL VariableAdd( Mesh % Variables, Mesh, Name='scan', DOFs=1, Values=sScan )
       END IF

       IF( DoControl ) THEN
         sSweep = 1.0_dp * iSweep
         CALL VariableAdd( Mesh % Variables, Mesh, Name='run', DOFs=1, Values=sSweep )
       END IF
              
       sPar(1) = 1.0_dp * ParEnv % MyPe 
       CALL VariableAdd( Mesh % Variables, Mesh, Name='Partition', DOFs=1, Values=sPar ) 

       IF( ListCheckPresent( CurrentModel % Simulation,'Parallel Slices') ) THEN
         CALL VariableAdd( Mesh % Variables, Mesh, Name='slice', DOFs=1, Values=sSlice )
         CALL VariableAdd( Mesh % Variables, Mesh, Name='slice ratio', DOFs=1, Values=sSliceRatio )
         CALL VariableAdd( Mesh % Variables, Mesh, Name='slice weight', DOFs=1, Values=sSliceWeight )
       END IF
       
       IF( ListCheckPresent( CurrentModel % Simulation,'Parallel Timestepping') ) THEN
         CALL VariableAdd( Mesh % Variables, Mesh, Name='time sector', DOFs=1, Values=sSector )
       END IF
             
       ! Add partition as a elemental field in case we have just one partition
       ! and have asked still for partitioning into many.
       IF( ParEnv % PEs == 1 .AND. ASSOCIATED( Mesh % Repartition ) ) THEN
         BLOCK
           REAL(KIND=dp), POINTER :: PartField(:)
           INTEGER, POINTER :: PartPerm(:)
           INTEGER :: i,n
           
           CALL Info('AddTimeEtc','Adding partitioning also as a field')
           
           n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           NULLIFY( PartField, PartPerm )
           ALLOCATE( PartField(n), PartPerm(n) )
           DO i=1,n
             PartPerm(i) = i
             PartField(i) = 1.0_dp * Mesh % RePartition(i)
           END DO
           
           CALL VariableAdd( Mesh % Variables, Mesh, Name='PartField',DOFs=1, &
               Values=PartField, Perm=PartPerm, TYPE=Variable_on_elements)
         END BLOCK
       END IF

       Mesh => Mesh % Next
        
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE AddTimeEtc
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sets initial conditions for the fields. 
!------------------------------------------------------------------------------
   SUBROUTINE SetInitialConditions()
!------------------------------------------------------------------------------
     USE DefUtils
     INTEGER :: DOFs
     CHARACTER(:), ALLOCATABLE :: str
     LOGICAL :: Found, NamespaceFound
     TYPE(Solver_t), POINTER :: Solver
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp),ALLOCATABLE :: Work(:)

     INTEGER :: i,j,k,l,m,vect_dof,real_dof,dim

     REAL(KIND=dp) :: nrm(3),t1(3),t2(3),vec(3),tmp(3),udot
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Nodes_t), SAVE :: Nodes
     LOGICAL :: nt_boundary, DG
     TYPE(Element_t), POINTER :: Element
     TYPE(Variable_t), POINTER :: var, vect_var
     LOGICAL :: AnyNameSpace
     TYPE(Element_t), POINTER :: p
     
     CALL Info('SetInitialConditions','Setting up initial conditions (if any)',Level=10)


     dim = CoordinateSystemDimension()

     IF (GetLogical(GetSimulation(),'Restart Before Initial Conditions',Found)) THEN
       CALL Restart()
       CALL InitCond()
     ELSE
       CALL InitCond()
       CALL Restart()
     END IF

         
!------------------------------------------------------------------------------
!    Make sure that initial values at boundaries are set correctly.
!    NOTE: This overrides the initial condition setting for field variables!!!!
!-------------------------------------------------------------------------------
     InitDirichlet = ListGetLogical( CurrentModel % Simulation, &
            'Initialize Dirichlet Conditions', GotIt ) 
     IF ( .NOT. GotIt ) InitDirichlet = .TRUE.

     AnyNameSpace = ListCheckPresentAnySolver( CurrentModel,'Namespace')
     NamespaceFound = .FALSE.
     
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

             IF( AnyNameSpace ) THEN
               str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
               IF (NamespaceFound) CALL ListPushNamespace(str)
             END IF               

             ! This seems to be a more robust marker for DG type
             DG = ( Var % Type == Variable_on_nodes_on_elements ) 
             
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
                   IF ( DG ) THEN
                     k = 0
                     p => Element % BoundaryInfo % Left                       
                     IF( ASSOCIATED( p ) ) THEN
                       DO i=1,p % TYPE % NumberOfNodes
                         IF(p % NodeIndexes(i) == Element % NodeIndexes(j) ) THEN
                           k = p % DGIndexes(i); EXIT
                         END IF
                       END DO
                       IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                     END IF
                     ! The active BC could be on either side!
                     ! If this is an internal BC this may really be poorly defined.
                     IF( k == 0 ) THEN
                       p => Element % BoundaryInfo % Right                                                
                       IF( ASSOCIATED( p ) ) THEN
                         DO i=1,p % TYPE % NumberOfNodes
                           IF(p % NodeIndexes(i) == Element % NodeIndexes(j) ) THEN
                             k = p % DGIndexes(i); EXIT
                           END IF
                         END DO
                         IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                       END IF
                     END IF
                   ELSE
                     k = Element % NodeIndexes(j)
                     IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
                   END IF

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
                      IF ( k>0 ) THEN
                        Var % PrevValues(k,2) = Work(j)
                        Var % PrevValues(k,6) = Work(j)
                      END IF
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
     INTEGER :: DOFs,i,i2,j,k,k1,k2,l,n,n2,m,nsize
     CHARACTER(:), ALLOCATABLE :: str, VarName
     LOGICAL :: Found, ThingsToDO, NamespaceFound, AnyNameSpace
     TYPE(Solver_t), POINTER :: Solver, CSolver
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp) :: Val
     REAL(KIND=dp),ALLOCATABLE :: Work(:)
     TYPE(ValueList_t), POINTER :: IC
     TYPE(Nodes_t) :: Nodes
     TYPE(GaussIntegrationPoints_t) :: IP
     REAL(KIND=dp), ALLOCATABLE :: Basis(:)
     REAL(KIND=dp) :: DetJ
     TYPE(Element_t), POINTER :: Element, p
     TYPE(ValueHandle_t) :: LocalSol_h
     LOGICAL :: Stat, FoundIC, PrevFoundIC
     INTEGER :: VarOrder, PrevBodyId
     !------------------------------------------------------------------------------

     AnyNameSpace = ListCheckPresentAnySolver( CurrentModel,'namespace')
     NameSpaceFound = .FALSE.
     
     Mesh => CurrentModel % Meshes
     DO WHILE( ASSOCIATED( Mesh ) )
       
       CALL SetCurrentMesh( CurrentModel, Mesh )

       IF( InfoActive( 30 ) ) THEN
         CALL Info('InitCond','Initial conditions for '//I2S(Mesh % MeshDim)//'D mesh:'//TRIM(Mesh % Name))
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 
           IF( ListCheckPresentAnyIC( CurrentModel, Var % Name ) ) THEN
             CALL VectorValuesRange(Var % Values,SIZE(Var % Values),'PreInit: '//TRIM(Var % Name))       
           END IF
           Var => Var % Next
         END DO
       END IF
       
       m = Mesh % MaxElementDofs
       n = Mesh % MaxElementNodes      
       ALLOCATE( Indexes(m), Work(m) , Basis(m), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

       ! First set the global variables and check whether there is anything left to do
       ThingsToDo = .FALSE.
       DO j=1,CurrentModel % NumberOfICs

         IC => CurrentModel % ICs(j) % Values
         
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 

           IF( .NOT. ASSOCIATED( Var % Values ) ) THEN
             Var => Var % Next
             CYCLE
           END IF

           IF( SIZE( Var % Values ) == 0 ) THEN
             Var => Var % Next
             CYCLE
           END IF
          
           Solver => Var % Solver
           IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver

           IF( AnyNameSpace ) THEN
             str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
             IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
           END IF
             
           ! global variable
           IF( SIZE( Var % Values ) == Var % DOFs .OR. Var % Type == Variable_global ) THEN
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
           
           Element =>  Mesh % Elements(t)
           
           i = Element % BodyId 
           IF( i == 0 ) CYCLE
           
           j = ListGetInteger(CurrentModel % Bodies(i) % Values, &
               'Initial Condition',GotIt, 1, CurrentModel % NumberOfICs )           
           IF ( .NOT. GotIt ) CYCLE
           
           IC => CurrentModel % ICs(j) % Values
           CurrentModel % CurrentElement => Element
           n = GetElementNOFNodes()
           
           Var => Mesh % Variables
           DO WHILE( ASSOCIATED(Var) ) 
             
             IF( .NOT. ASSOCIATED( Var % Values ) ) THEN
               Var => Var % Next
               CYCLE
             END IF
             
             IF( SIZE( Var % Values ) == 0 ) THEN
               Var => Var % Next
               CYCLE
             END IF

             IF( t == 1 ) THEN
               CALL Info('InitCond','Trying to initialize variable: '//TRIM(Var % Name),Level=20)
             END IF
               
             Solver => Var % Solver
             IF ( .NOT. ASSOCIATED(Solver) ) THEN
              Solver => CurrentModel % Solver
            END IF
            CSolver => CurrentModel % Solver
            CurrentModel % Solver => Solver

             IF( AnyNameSpace ) THEN
               str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
               IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
             END IF
               
             ! global variables were already set
             IF( SIZE( Var % Values ) == Var % DOFs .OR. Var % Type == Variable_global ) THEN
               CONTINUE
               
             ELSE IF( Var % TYPE == Variable_on_elements ) THEN
               IF( Var % DOFs > 1 ) THEN
                 CALL Fatal('InitCond','Initialization only for scalar elements fields!')
               END IF
               
               Work(1:n) = GetReal( IC, Var % Name, GotIt )
               IF ( GotIt ) THEN
                 k1 = Element % ElementIndex 
                 IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                 IF ( k1>0 ) Var % Values(k1) = SUM( Work(1:n) ) / n
               END IF               
               
             ELSE IF( Var % TYPE == Variable_on_gauss_points ) THEN
               ! We do this elsewhere in a more efficient manner
               CONTINUE
               
             ELSE IF ( Var % DOFs == 1 ) THEN
                
               Work(1:n) = ListGetReal( IC, Var % Name, n, Element % NodeIndexes, GotIt )

               IF ( GotIt ) THEN
                 ! Sometimes you may have both DG and bubbles,
                 ! this way DG always has priority. 
                 IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN 
                   DO k=1,n
                     IF( ASSOCIATED( Element % DGIndexes) ) THEN
                       ! DG variable has always a permutation associated to it!
                       k1 = Var % Perm(Element % DgIndexes(k))
                     ELSE                                            
                       k1 = 0
                       p => Element % BoundaryInfo % Left                       
                       IF( ASSOCIATED( p ) ) THEN
                         DO i=1,p % TYPE % NumberOfNodes
                           IF(p % NodeIndexes(i) == Element % NodeIndexes(k) ) THEN
                             k1 = Var % Perm(p % DGIndexes(i)); EXIT
                           END IF
                         END DO
                       END IF
                       IF( k1 == 0 ) THEN
                         p => Element % BoundaryInfo % Right                       
                         IF( ASSOCIATED( p ) ) THEN
                           DO i=1,p % TYPE % NumberOfNodes
                             IF(p % NodeIndexes(i) == Element % NodeIndexes(k) ) THEN
                               k1 = Var % Perm(p % DGIndexes(i)); EXIT
                             END IF
                           END DO
                         END IF
                       END IF
                     END IF
                     
                     IF ( k1>0 ) Var % Values(k1) = Work(k)
                   END DO
 
                 ELSE
                   DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   DO k=1,n
                     k1 = Indexes(k)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % Values(k1) = Work(k)
                   END DO
                 END IF
                   
               END IF

               IF ( Transient .AND. Solver % TimeOrder==2 ) THEN
                 Work(1:n) = GetReal( IC, TRIM(Var % Name) // ' Velocity', GotIt )
                 IF ( GotIt ) THEN
                   IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN 
                     Indexes(1:n) = Element % DgIndexes(1:n)
                   ELSE
                     DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   END IF

                   DO k=1,n
                     k1 = Indexes(k)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) Var % PrevValues(k1,1) = Work(k)
                   END DO
                 END IF
                 Work(1:n) = GetReal( IC, TRIM(Var % Name) // ' Acceleration', GotIt )
                 IF ( GotIt ) THEN
                   IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN 
                     Indexes(1:n) = Element % DgIndexes(1:n)
                   ELSE
                     DOFs = GetElementDOFs( Indexes, USolver=Var % Solver )
                   END IF

                   DO k=1,n
                     k1 = Indexes(k)
                     IF ( ASSOCIATED(Var % Perm) ) k1 = Var % Perm(k1)
                     IF ( k1>0 ) THEN
                       Var % PrevValues(k1,2) = Work(k)
                       Var % PrevValues(k1,6) = Work(k)
                     END IF
                   END DO
                 END IF
               END IF
               
               IF(ASSOCIATED(Mesh % Edges)) THEN
                 IF ( i<=Mesh % NumberOfBulkElements) THEN
                   Gotit = ListCheckPresent( IC, TRIM(Var % Name)//' {e}' )
                   IF ( Gotit ) THEN
                     DO k=1,Element % TYPE % NumberOfedges
                       Edge => Mesh % Edges(Element % EdgeIndexes(k))
                       l = Var % Perm(Element % EdgeIndexes(k)+Mesh % NumberOfNodes)
                       IF ( l>0 ) THEN
                         CALL VectorElementEdgeDOFs( IC, &
                             Edge, Edge % TYPE % NumberOfNodes, Element, n, &
                             TRIM(Var % Name)//' {e}', Work )
                         Var % Values(l) = Work(1)
                       END IF
                     END DO
                   END IF
                 END IF
               END IF
               
             ELSE
               CALL ListGetRealArray( IC, &
                   Var % Name, WorkA, n, Element % NodeIndexes, gotIt )
               
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
           CSolver => CurrentModel % Solver
         END DO
         
         ! Here we do just the gauss point values for now.
         ! It would really make sense to do the ICs in this order since probably
         ! there are quite few IC variables to set but quite many elements.
         
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 

           VarName = TRIM(Var % Name)
           Solver => Var % Solver
           IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver
           
           IF( AnyNameSpace ) THEN
             str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
             IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
           END IF
           
           VarOrder = -1
           DO VarOrder = 0, 2
             IF( VarOrder == 0 ) THEN
               VarName = TRIM(Var % Name)
             ELSE IF( VarOrder == 1 ) THEN
               VarName = TRIM(Var % Name )//' Velocity'
             ELSE IF( VarOrder == 2 ) THEN
               VarName = TRIM(Var % Name )//' Acceleration'
             END IF

             !CALL ListInitElementKeyword( LocalSol_h,'Initial Condition',VarName, &
             !    FoundSomewhere = Found )
             Found = ListCheckPresentAnyIC( CurrentModel, VarName )

             IF( Found .AND. VarOrder > 0 ) THEN              
               IF ( .NOT. ( Transient .AND. Solver % TimeOrder==2 ) ) THEN
                 CALL Warn('InitCond','We can only set timederivative for transients')
                 Found = .FALSE.
               END IF
             END IF

             IF( Found ) THEN
               
               CALL ListInitElementKeyword( LocalSol_h,'Initial Condition',VarName )
               
               IF( Var % TYPE /= Variable_on_gauss_points ) CYCLE

               CALL Info('InitCond','Initializing gauss point field: '//TRIM(VarName),Level=7)               
               IF( Var % DOFs > 1 ) THEN
                 CALL Fatal('InitCond','Initialization only for scalar elemental fields!')
               END IF
               
100            PrevBodyId = -1 
               DO t=1, Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
                 
                 Element => Mesh % Elements(t)

                 i = Element % BodyId 
                 IF( i == 0 ) CYCLE         

                 IF( i == PrevBodyId ) THEN
                   FoundIC = PrevFoundIC
                 ELSE
                   j = ListGetInteger(CurrentModel % Bodies(i) % Values, &
                       'Initial Condition',FoundIC, 1, CurrentModel % NumberOfICs )           
                   IF ( FoundIC ) THEN                       
                     IC => CurrentModel % ICs(j) % Values
                     FoundIC = ListCheckPresent( IC, VarName )
                   END IF
                   PrevFoundIC = FoundIC 
                 END IF

                 IF( .NOT. FoundIC ) CYCLE

                 CurrentModel % CurrentElement => Element
                 n = GetElementNOFNodes()                 
                 
                 k1 = Var % Perm( Element % ElementIndex )
                 k2 = Var % Perm( Element % ElementIndex + 1 )

                 IF( k2- k1 > 0 ) THEN
                   
                   IP = GaussPointsAdapt( Element, Solver )
                   
                   IF( k2 - k1 /= Ip % n ) THEN
                     CALL Info('InitCond','Number of Gauss points has changed, redoing permutations!',Level=8)
                     CALL UpdateIpPerm( Solver, Var % Perm )
                     nsize = MAXVAL( Var % Perm )
                     
                     CALL Info('InitCond','Total number of new IP dofs: '//I2S(nsize))
                     
                     IF( SIZE( Var % Values ) /= Var % Dofs * nsize ) THEN
                       DEALLOCATE( Var % Values )
                       ALLOCATE( Var % Values( nsize * Var % Dofs ) )
                     END IF
                     Var % Values = 0.0_dp
                     GOTO 100 
                   END IF

                   Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
                   Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
                   Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)

                   DO k=1,IP % n
                     stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
                         IP % W(k), detJ, Basis )

                     val = ListGetElementReal( LocalSol_h,Basis,Element,Found,GaussPoint=k)

                     IF( VarOrder == 0 ) THEN
                       Var % Values(k1+k) = val
                     ELSE 
                       Var % PrevValues(k1+k,VarOrder) = val
                     END IF

                   END DO
                   
                 END IF 
               END DO
             END IF
           END DO
           IF(NamespaceFound) CALL ListPopNamespace()
           Var => Var % Next
         END DO
       END IF

       DEALLOCATE( Indexes, Work, Basis, Nodes % x, Nodes % y, Nodes % z)
       
       IF( InfoActive( 20 ) ) THEN
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 
           IF( ListCheckPresentAnyIC( CurrentModel, Var % Name ) ) THEN
             CALL VectorValuesRange(Var % Values,SIZE(Var % Values),'PostInit: '//TRIM(Var % Name))       
           END IF
           Var => Var % Next
         END DO
       END IF
     
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
     LOGICAL :: Gotit, DoIt
     INTEGER :: i, j, k, l
     REAL(KIND=dp) :: StartTime
     TYPE(Mesh_t), POINTER :: Mesh, pMesh
     TYPE(ValueList_t), POINTER :: RestartList
     LOGICAL, ALLOCATABLE :: MeshDone(:)
     INTEGER, POINTER :: MeshesToRestart(:)
     LOGICAL :: CheckMesh, DoMesh, isParallel
!------------------------------------------------------------------------------

     
     ! Count the number of meshes first so that we can identify them
     j = 0
     pMesh => CurrentModel % Meshes       
     DO WHILE( ASSOCIATED(pMesh) ) 
       j = j + 1
       pMesh => pMesh % Next
     END DO
     ALLOCATE( MeshDone( j ) )
     MeshDone = .FALSE.

     
     ! Do Solver-mesh specific restart only
     !-----------------------------------------------------------------
     IF ( ListCheckPresentAnySolver( CurrentModel,'Restart File') ) THEN
       DO i=1, CurrentModel % NumberOfSolvers
         RestartList => CurrentModel % Solvers(i) % Values 
         
         RestartFile = ListGetString( RestartList, 'Restart File', GotIt )
         IF ( GotIt ) THEN

           Mesh => CurrentModel % Solvers(i) % Mesh 
           IF( .NOT. ASSOCIATED(Mesh) ) THEN
             CALL Warn('Restart','Solver has no mesh associated!')
             CYCLE
           END IF
           
           DoIt = .TRUE.
           pMesh => CurrentModel % Meshes       
           j = 0
           DO WHILE( ASSOCIATED(pMesh) ) 
             j = j + 1 
             pMesh => pMesh % Next
             IF( ASSOCIATED( Mesh, pMesh ) ) THEN
               IF( MeshDone(j) ) THEN
                 DoIt = .FALSE.
               ELSE
                 MeshDone(j) = .TRUE.
               END IF
             END IF
           END DO

           ! Variables for this mesh already done!
           IF(.NOT. DoIt ) CYCLE
           
           CALL Info('Restart','Perfoming solver specific Restart for: '//TRIM(Mesh % Name),Level=6)
           IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
             OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
           ELSE
             OutputName = RestartFile
           END IF

           ! If we have single mesh we have the luxury of using either parallel or serial restart
           isParallel = ParEnv % PEs > 1
           IF(isParallel .AND. Mesh % SingleMesh ) THEN
             isParallel = ListGetLogical( RestartList,'Restart Parallel',Found )
           END IF                        
           IF(isParallel) OutputName = OutputName // '.' // i2s(ParEnv % MyPe)

           CALL SetCurrentMesh( CurrentModel, Mesh )
           
           k = ListGetInteger( RestartList,'Restart Position',GotIt, minv=0 )
           CALL LoadRestartFile( OutputName, k, Mesh, SolverId = i )
           
           StartTime = ListGetConstReal( RestartList ,'Restart Time',GotIt)
           IF( GotIt ) THEN
             Var  => VariableGet( Mesh % Variables, 'Time' )
             IF ( ASSOCIATED( Var ) )  Var % Values(1) = StartTime
           END IF
         END IF
         
       END DO
     END IF

     ! Do the standard global restart
     !-----------------------------------------------------------------
     RestartList => CurrentModel % Simulation

     ! We may suppress restart from certain meshes.
     ! This was initially only related to calving, but no need to limit to that.
     l = 0
     MeshesToRestart => ListGetIntegerArray(RestartList,&
         'Meshes To Restart', CheckMesh )

     RestartFile = ListGetString( RestartList, 'Restart File', GotIt )
     IF ( GotIt ) THEN
       k = ListGetInteger( RestartList,'Restart File Number',GotIt)
       IF( GotIt ) RestartFile = TRIM(RestartFile)//'_'//I2S(k)//'nc'
              
       k = ListGetInteger( RestartList,'Restart Position',GotIt, minv=0 )
       Mesh => CurrentModel % Meshes

       j = 0
       DO WHILE( ASSOCIATED(Mesh) ) 
         j = j + 1

         ! Make sure that if a mesh has already been restarted 
         ! it is not being done again. 
         IF( MeshDone(j) ) THEN
           CALL Info('Restart','Already done mesh: '//TRIM(Mesh % Name))
           Mesh => Mesh % Next
           CYCLE
         END IF
         MeshDone(j) = .TRUE.
         
         CALL Info('Restart','Perfoming global Restart for: '//TRIM(Mesh % Name),Level=6)
         
         IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
           OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
         ELSE
           OutputName = TRIM(RestartFile)
         END IF

         isParallel = ParEnv % PEs > 1
         IF(isParallel .AND. Mesh % SingleMesh ) THEN
           isParallel = ListGetLogical( RestartList,'Restart Parallel',Found )
         END IF                  
         IF(isParallel ) OutputName = TRIM(OutputName) // '.' // i2s(ParEnv % MyPe)
         
         l = l+1

         DoMesh = .TRUE.
         IF(CheckMesh) THEN
           DoMesh = (ANY(MeshesToRestart == l))
         END IF

         IF( DoMesh ) THEN
           CALL SetCurrentMesh( CurrentModel, Mesh )
           CALL LoadRestartFile( OutputName, k, Mesh )
           
           StartTime = ListGetConstReal( RestartList ,'Restart Time',GotIt)
           IF( GotIt ) THEN
             Var  => VariableGet( Mesh % Variables, 'Time' )
             IF ( ASSOCIATED( Var ) )  Var % Values(1) = StartTime
           END IF
         END IF           
         Mesh => Mesh % Next
       END DO
     ELSE
       StartTime = ListGetConstReal( CurrentModel % Simulation  ,'Start Time',GotIt)
       Mesh => CurrentModel % Meshes

       j = 0
       IF( GotIt ) THEN
         DO WHILE( ASSOCIATED(Mesh) ) 
           j = j + 1

           ! Make sure that if a mesh has already been restarted 
           ! it is not being done again. 
           IF( MeshDone(j) ) THEN
             CALL Info('Restart','Already done mesh: '//TRIM(Mesh % Name))
             Mesh => Mesh % Next
             CYCLE
           END IF
           MeshDone(j) = .TRUE.
           Var  => VariableGet( Mesh % Variables, 'Time' )
           IF ( ASSOCIATED( Var ) )  Var % Values(1) = StartTime
           WRITE(Message,*) 'Found "Start Time" (without restart) and setting time-variable to t=',StartTime
           CALL INFO("Restart",Message,Level=1)
         END DO
       END IF
     END IF
 
     
!------------------------------------------------------------------------------
   END SUBROUTINE Restart
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Execute the individual solvers in defined sequence. 
!------------------------------------------------------------------------------
   SUBROUTINE ExecSimulation(TimeIntervals,  CoupledMinIter, &
              CoupledMaxIter, OutputIntervals, Transient, Scanning)
     IMPLICIT NONE
      INTEGER :: TimeIntervals,CoupledMinIter, CoupledMaxIter,OutputIntervals(:)
      LOGICAL :: Transient,Scanning
!------------------------------------------------------------------------------
     INTEGER :: interval, timestep, i, j, k, n
     REAL(KIND=dp) :: dt, ddt, dtfunc, timeleft
     INTEGER :: cum_timestep
     INTEGER, SAVE ::  stepcount, RealTimestep
     LOGICAL :: ExecThis,SteadyStateReached=.FALSE.,PredCorrControl, &
         DivergenceControl, HaveDivergence
     REAL(KIND=dp) :: CumTime, MaxErr, AdaptiveLimit, &
         AdaptiveMinTimestep, AdaptiveMaxTimestep, timePeriod
     INTEGER :: SmallestCount, AdaptiveKeepSmallest, StepControl=-1, nSolvers
     LOGICAL :: AdaptiveTime = .TRUE., AdaptiveRough, AdaptiveSmart, Found, DoIt
     INTEGER :: AllocStat
     REAL(KIND=dp) :: AdaptiveIncrease, AdaptiveDecrease     
     TYPE(Solver_t), POINTER :: Solver    
     TYPE AdaptiveVariables_t 
       TYPE(Variable_t) :: Var
       REAL(KIND=dp) :: Norm
     END TYPE AdaptiveVariables_t
     TYPE(AdaptiveVariables_t), ALLOCATABLE, SAVE :: AdaptVars(:)     
     REAL(KIND=dp) :: newtime, prevtime=0, maxtime, exitcond
     INTEGER, SAVE :: PrevMeshI = 0
     INTEGER :: nPeriodic, nSlices, nTimes, iSlice, iTime
     LOGICAL :: ParallelTime, ParallelSlices, IsPeriodic
     CHARACTER(*), PARAMETER :: Caller = 'ExecSimulation'
     
     !$OMP PARALLEL
     IF(.NOT.GaussPointsInitialized()) CALL GaussPointsInit()
     !$OMP END PARALLEL

     nSolvers = CurrentModel % NumberOfSolvers
     DO i=1,nSolvers
        Solver => CurrentModel % Solvers(i)
        IF ( Solver % PROCEDURE==0 ) CYCLE
        DoIt = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL )
        IF(.NOT. DoIt) THEN
          DoIt = ListGetLogical( Solver % Values,'Before All',Found ) .OR. &
              ListGetLogical( Solver % Values,'Before Simulation',Found )
        END IF
        
        IF( DoIt ) THEN
          ! solver to be called prior to time looping can never be transient
          dt = 1.0_dp
          CALL SolverActivate( CurrentModel,Solver,dt,.FALSE. )
        END IF
     END DO

     IF( ListGetLogical( CurrentModel % Simulation,'Calculate Mesh Pieces',Found ) .OR. &
         ListCheckPresent( CurrentModel % Simulation,'Desired Mesh Pieces') ) THEN
       CALL CalculateMeshPieces( CurrentModel % Mesh ) 
     END IF

     ! Predictor-Corrector time stepping control 
     PredCorrControl = ListGetLogical( CurrentModel % Simulation, &
         'Predictor-Corrector Control', gotIt)

     ! Divergence control 
     DivergenceControl = ListGetLogical( CurrentModel % Simulation, &
         'Convergence Control', gotIt)     

     AdaptiveTime = ListGetLogical( CurrentModel % Simulation, &
         'Adaptive Timestepping', GotIt )

     stepcount = 0
     DO interval = 1, TimeIntervals
        stepcount = stepcount + Timesteps(interval)
     END DO 

     dt = 1.0_dp
     cum_Timestep = 0
     ddt = -1.0_dp

     ! For parallel timestepping we need to divide the periodic timesteps for each partition. 
     ParallelTime = ListGetLogical( CurrentModel % Simulation,'Parallel Timestepping', GotIt ) &
         .AND. ( ParEnv % PEs > 1 ) 

     ! For parallel slices we need to introduce the slices
     ! Let this be active even for serial case to have consistent setup
     ParallelSlices = ListGetLogical( CurrentModel % Simulation,'Parallel Slices',GotIt ) 
         !.AND. ( ParEnv % PEs > 1 )

     IF( ParallelTime .OR. ParallelSlices ) THEN
       IF( ParEnv % PEs > 1 ) THEN
         IF(.NOT. ListGetLogical( CurrentModel % Simulation,'Single Mesh',GotIt ) ) THEN
           CALL Fatal(Caller,'Parallel time and slices only available with "Single Mesh"')
         END IF
       END IF
     END IF

     
     nSlices = 1
     nTimes = 1
     iTime = 0
     iSlice = 0
     
     IF( ParallelTime .AND. ParallelSlices ) THEN
       nSlices = ListGetInteger( CurrentModel % Simulation,'Number Of Slices',GotIt)
       IF(GotIt) THEN
         IF( nSlices > ParEnv % PEs ) THEN
           CALL Fatal(Caller,'"Number Of Slices" cannot be be larger than #np')
         END IF
       ELSE
         IF( ParEnv % PEs == 1 ) THEN
           CALL ListAddInteger( CurrentModel % Simulation,'Number Of Slices',nSlices)
         ELSE
           CALL Fatal(Caller,'We need "Number Of Slices" with parallel timestepping')
         END IF
       END IF
       IF( MODULO( ParEnv % PEs, nSlices ) /= 0 ) THEN
         CALL Fatal(Caller,'For hybrid parallellism #np must be divisible with "Number of Slices"')
       END IF
       nTimes = ParEnv % PEs / nSlices 
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Times',nTimes )
       iSlice = MODULO( ParEnv % MyPe, nSlices ) 
       iTime = ParEnv % MyPe / nSlices
     ELSE IF( ParallelTime ) THEN
       nTimes = ParEnv % PEs
       iTime = ParEnv % MyPe
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Times',nTimes )
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Slices',nSlices )
       CALL Info(Caller,'Setting one time sector for each partition!')
     ELSE IF( ParallelSlices ) THEN
       nSlices = ParEnv % PEs
       iSlice = ParEnv % MyPe
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Times',nTimes )
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Slices',nSlices )
       CALL Info(Caller,'Setting one slice for each partition!')
     END IF

     IF( nTimes > 1 ) THEN
       DO i=1,SIZE(Timesteps,1)
         IF( MODULO( Timesteps(i), nTimes ) /= 0 ) THEN
           CALL Fatal(Caller,'"Timestep Intervals" should be divisible by nTimes: '//I2S(nTimes))
         END IF
         Timesteps(i) = Timesteps(i) / nTimes
       END DO
       CALL Info(Caller,'Divided timestep intervals equally for each partition!',Level=4)
     END IF

     IsPeriodic = ListCheckPrefix( CurrentModel % Simulation,'Periodic Time') .OR. &
         ListCheckPresent( CurrentModel % Simulation,'Time Period' )
     
     nPeriodic = ListGetInteger( CurrentModel % Simulation,'Periodic Timesteps',GotIt )
     IF( ParallelTime ) THEN
       IF( nPeriodic <= 0 ) THEN
         CALL Fatal(Caller,'Parallel timestepping requires "Periodic Timesteps"')
       END IF
       IF( MODULO( nPeriodic, nTimes ) /= 0 ) THEN
         CALL Fatal(Caller,'For parallel timestepping "Periodic Timesteps" must be divisible by #np')
       END IF
       nPeriodic = nPeriodic / nTimes
     END IF
     
     IF( ListGetLogical( CurrentModel % Simulation,'Parallel Slices',GotIt ) ) THEN
       IF( nSlices <= 1 ) THEN
         sSlice = 0.0_dp
         sSliceRatio = 0.0_dp
         sSliceWeight = 1.0_dp
       ELSE
         sSlice = 1.0_dp * iSlice
         sSliceRatio = ( iSlice + 0.5_dp ) / nSlices - 0.5_dp
         sSliceWeight = 1.0_dp / nSlices 
       END IF

       BLOCK
         REAL :: z_min, z_max, z_mid
         z_min = ListGetConstReal(CurrentModel % Simulation,'Extruded Min Coordinate',GotIt)
         z_max = ListGetConstReal(CurrentModel % Simulation,'Extruded Max Coordinate',GotIt)
         
         IF( GotIt .AND. nSlices > 1) THEN
           z_mid = 0.5_dp * ( z_min + z_max ) 
           CALL Info(Caller,'Moving parallel slices in z-direction!',Level=6)
           i = CurrentModel % Mesh % NumberOfNodes 
           CurrentModel % Mesh % Nodes % z(1:i) = z_mid + (z_max-z_min) * sSliceRatio(1)  
         END IF
       END BLOCK
     END IF

     IF( ListGetLogical( CurrentModel % Simulation,'Parallel Timestepping',GotIt ) ) THEN
       IF( nTimes <= 1 ) THEN
         sSector = 0.0_dp
       ELSE         
         sSector = 1.0_dp * iTime 
       END IF
     END IF       
     
     DO interval = 1,TimeIntervals
       
!------------------------------------------------------------------------------
!      go through number of timesteps within an interval
!------------------------------------------------------------------------------       
       timePeriod = ListGetCReal(CurrentModel % Simulation, 'Time Period',gotIt)       
       IF(.NOT.GotIt) timePeriod = HUGE(timePeriod)
         
       IF(GetNameSpaceCheck()) THEN
         IF(Scanning) THEN
           CALL ListPushNamespace('scan:')
         ELSE IF(Transient) THEN
           CALL ListPushNamespace('time:')
         ELSE
           CALL ListPushNamespace('steady:')
         END IF
       END IF

       RealTimestep = 1

       timestep = 1
       DO WHILE(timestep <= Timesteps(interval))

         cum_Timestep = cum_Timestep + 1
         sStep(1) = cum_Timestep

         IF ( GetNamespaceCheck() ) THEN
           IF( Scanning ) THEN
             CALL ListPushNamespace('scan '//i2s(cum_Timestep)//':')
           ELSE IF ( Transient ) THEN
             CALL ListPushNamespace('time '//i2s(cum_Timestep)//':')
           ELSE
             CALL ListPushNamespace('steady '//i2s(cum_Timestep)//':')
           END IF
         END IF

         ! Sometimes when timestep depends on time we need to have first timestep size
         ! given separately to avoid problems. 
         GotIt = .FALSE.
         IF( cum_Timestep == 1 ) THEN
           dtfunc = ListGetCReal( CurrentModel % Simulation,'First Timestep Size',GotIt )
           IF(GotIt) dt = dtfunc
         END IF
         
         IF ( ( Transient .OR. Scanning ) .AND. .NOT. GotIt ) THEN
           dtfunc = ListGetCReal( CurrentModel % Simulation,'Timestep Function',GotIt )
           IF(GotIt) THEN
             CALL Warn(Caller,'Obsolete keyword > Timestep Function < , use > Timestep Size < instead')
           ELSE           
             dtfunc = ListGetCReal( CurrentModel % Simulation,'Timestep Size', gotIt)
           END IF
           IF( GotIt ) THEN             
             BLOCK                 
               TYPE(Variable_t), POINTER :: tVar
               INTEGER :: dtiter 
               REAL(KIND=dp) :: t0
               dtiter = ListGetInteger( CurrentModel % Simulation,'Timestep Size Iterations',Found)

               ! If timestep size depends on time i.e. dt=dt(t) we may sometimes want to iterate
               ! so that the timestep is consistent with the new time t+dt i.e. dt=dt(t+dt).
               IF( dtiter > 0 ) THEN
                 tVar => VariableGet( CurrentModel % Mesh % Variables, 'time' )
                 IF(ASSOCIATED(tVar)) THEN
                   t0 = tVar % Values(1)
                   ! Iterate for consistent time + timestep size. 
                   DO i=1,dtiter  
                     tVar % Values(1) = t0 + dtfunc
                     dtfunc = ListGetCReal( CurrentModel % Simulation,'Timestep Size' )
                   END DO
                   ! Revert back to original time, this will be added later on again.
                   tVar % Values(1) = t0
                 END IF
               END IF
             END BLOCK
           END IF
           
           IF(GotIt) THEN
             dt = dtfunc
             IF(dt < EPSILON(dt) ) THEN
               WRITE(Message,'(A,ES12.3)') 'Timestep smaller than epsilon: ',dt
               CALL Fatal(Caller, Message)
             END IF             
           ELSE 
             dt = TimestepSizes(interval,1)
             IF( GotTimestepRatios ) THEN
               BLOCK 
                 REAL(KIND=dp) :: q
                 q = TimestepRatios(interval,1)
                 IF( ABS(1-q) > EPSILON(q) ) dt = dt * q**(timestep-1)
               END BLOCK
             END IF
             
           END IF
         END IF

!------------------------------------------------------------------------------
         ! Predictor-Corrector time stepping control 
         IF ( PredCorrControl ) THEN 
           CALL PredictorCorrectorControl( CurrentModel, dt, timestep )
         END IF

!------------------------------------------------------------------------------
         
         IF( AdaptiveTime ) THEN
           AdaptiveLimit = ListGetConstReal( CurrentModel % Simulation, &
               'Adaptive Time Error', GotIt )       
           IF ( .NOT. GotIt ) THEN 
             CALL Fatal('MAIN','Adaptive Time Error must be given for ' // &
                 'adaptive stepping scheme.')
           END IF
           AdaptiveKeepSmallest = ListGetInteger( CurrentModel % Simulation, &
               'Adaptive Keep Smallest', GotIt, minv=0  )         
        END IF
         
         IF( AdaptiveTime .OR. DivergenceControl ) THEN
           AdaptiveMaxTimestep = ListGetConstReal( CurrentModel % Simulation, &
               'Adaptive Max Timestep', GotIt )
           IF ( GotIt ) THEN
             AdaptiveMaxTimestep =  MIN(AdaptiveMaxTimeStep, dt)
           ELSE
             AdaptiveMaxTimestep =  dt
           END IF
                        
           AdaptiveMinTimestep = ListGetConstReal( CurrentModel % Simulation, &
               'Adaptive Min Timestep', GotIt )
           IF(.NOT. GotIt) AdaptiveMinTimestep = 1.0e-8 * AdaptiveMaxTimestep
           
           AdaptiveIncrease =  ListGetConstReal( CurrentModel % Simulation, &
               'Adaptive Increase Coefficient', GotIt )
           IF(.NOT. GotIt) AdaptiveIncrease = 2.0_dp
           
           AdaptiveDecrease =  ListGetConstReal( CurrentModel % Simulation, &
               'Adaptive Decrease Coefficient', GotIt )
           IF(.NOT. GotIt) AdaptiveDecrease = 0.5_dp
           
           AdaptiveRough = ListGetLogical( CurrentModel % Simulation, &
               'Adaptive Rough Timestep', GotIt )           

           AdaptiveSmart = ListGetLogical( CurrentModel % Simulation, &
               'Adaptive Smart Timestep', GotIt )           
         END IF

         
!------------------------------------------------------------------------------
         IF(cum_Timestep == 1 .AND. ListGetLogical( CurrentModel % Simulation,'Timestep Start Zero',GotIt) ) THEN
           CALL Info(Caller,'Not advancing the 1st timestep!')
         ELSE
           sTime(1) = sTime(1) + dt
         END IF
         
         IF( nPeriodic > 0 ) THEN
           IF( ParallelTime ) THEN
             timePeriod = nTimes * nPeriodic * dt                        
             IF( cum_Timestep == 1 ) THEN
               sTime(1) = sTime(1) + iTime * nPeriodic * dt
             ELSE IF( MODULO( cum_Timestep, nPeriodic ) == 1 ) THEN
               CALL Info(Caller,'Making jump in time-parallel scheme!')
               sTime(1) = sTime(1) + nPeriodic * (nTimes - 1) * dt
             END IF
           ELSE
             timePeriod = nPeriodic * dt           
           END IF
         END IF

         IF( isPeriodic ) THEN
           sPeriodicCycle(1) = sTime(1) / timePeriod 
           sPeriodicTime(1) = MODULO( sTime(1), timePeriod )
         END IF
           
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

!-----------------------------------------------------------------------------
         IF( ListCheckPresent( CurrentModel % Simulation,'Rotor Angle') ) THEN
           BLOCK
             REAL(KIND=dp) :: PrevAngle
             PrevAngle = sAngle(1)
             sAngle(1) = ListGetCReal( CurrentModel % Simulation,'Rotor Angle')
             sAngleVelo(1) = (sAngle(1)-PrevAngle)/dt
             WRITE(Message,'(A,ES12.3)') '"Rotor Angle" set to value: ',sAngle(1)
             CALL Info(Caller,Message,Level=6)
             WRITE(Message,'(A,ES12.3)') '"Rotor Velo" set to value: ',sAngleVelo(1)
             CALL Info(Caller,Message,Level=6)
           END BLOCK
         END IF       
         
!------------------------------------------------------------------------------

         BLOCK
           TYPE(Mesh_t), POINTER :: Mesh
           REAL(KIND=dp) :: MeshR
           CHARACTER(:), ALLOCATABLE :: MeshStr
           
           IF( ListCheckPresent( GetSimulation(), 'Mesh Name Index') ) THEN
             ! we cannot have mesh depend on "time" or "timestep" if they are not available as
             ! variables. 
             Mesh => CurrentModel % Meshes             
             CALL VariableAdd( Mesh % Variables, Mesh, Name='Time',DOFs=1, Values=sTime )
             CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep', DOFs=1, Values=sStep )

             MeshR = GetCREal( GetSimulation(), 'Mesh Name Index',gotIt )            
             i = NINT( MeshR )

             IF( i > 0 .AND. i /= PrevMeshI ) THEN                             
               MeshStr = ListGetString( GetSimulation(),'Mesh Name '//I2S(i),GotIt)
               IF( GotIt ) THEN
                 CALL Info(Caller,'Swapping mesh to: '//MeshStr,Level=5)
               ELSE
                 CALL Fatal(Caller,'Could not find >Mesh Name '//I2S(i)//'<')
               END IF
               CALL SwapMesh( CurrentModel, Mesh, MeshStr )
               PrevMeshI = i
             END IF
           END IF
         END BLOCK


!------------------------------------------------------------------------------
         IF ( ParEnv % MyPE == 0 ) THEN
           CALL Info( 'MAIN', ' ', Level=3 )
           CALL Info( 'MAIN', '-------------------------------------', Level=3 )

           IF ( Transient .OR. Scanning ) THEN
             WRITE( Message,'(A,ES12.3)') 'Time: '//i2s(cum_Timestep)//'/'// &
                   i2s(stepcount)//':', sTime(1)
             CALL Info( 'MAIN', Message, Level=3 )

             newtime= RealTime()

             IF( cum_Timestep > 1 ) THEN
               maxtime = ListGetConstReal( CurrentModel % Simulation,'Real Time Max',GotIt)
               IF( GotIt ) THEN
                 WRITE( Message,'(A,F8.3)') 'Fraction of real time left: ',&
                     1.0_dp-(RealTime()-RT0)/ maxtime
               END IF

               ! This gives elapsed time in same format as estimated time
               timeleft = (newtime-RT0)
               IF( timeleft > 1 ) THEN
                 IF (timeleft >= 10 * 3600) THEN 
                   WRITE( Message,'(A)') 'Elapsed time: '//I2S(NINT(timeleft/3600))//' hours'
                 ELSE IF (timeleft >= 3600) THEN   
                   WRITE( Message,'(A,F5.1,A)') 'Elapsed time:',timeleft/3600,' hours'
                 ELSE IF(timeleft >= 60) THEN 
                   WRITE( Message,'(A,F5.1,A)') 'Elapsed time:',timeleft/60,' minutes'
                 ELSE                         
                   WRITE( Message,'(A,F5.1,A)') 'Elapsed time:',timeleft,' seconds'
                 END IF
                 CALL Info( 'MAIN', Message, Level=6 )
               END IF
               
               ! Compute estimated time left in seconds
               timeleft = (stepcount-(cum_Timestep-1))*(newtime-prevtime)
               
               ! No sense to show too short estimated times
               IF( timeleft > 1 ) THEN
                 IF (timeleft >= 10 * 3600) THEN ! >10 hours
                   WRITE( Message,'(A)') 'Estimated time left: '//I2S(NINT(timeleft/3600))//' hours'
                 ELSE IF (timeleft >= 3600) THEN   ! >1 hours
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

         IF ( Transient .AND. AdaptiveTime ) THEN 

            IF(.NOT. ALLOCATED( AdaptVars ) ) THEN
              ALLOCATE( AdaptVars( nSolvers ), STAT = AllocStat )
              IF( AllocStat /= 0 ) CALL Fatal(Caller,'Allocation error for AdaptVars')
              
              DO i=1,nSolvers
                Solver => CurrentModel % Solvers(i)

                NULLIFY( AdaptVars(i) % Var % Values )
                NULLIFY( AdaptVars(i) % Var % PrevValues )

                IF( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE
                IF( .NOT. ASSOCIATED( Solver % Variable  % Values ) ) CYCLE
                CALL Info(Caller,'Allocating adaptive work space for: '//I2S(i),Level=12)
                j = SIZE( Solver % Variable % Values )
                ALLOCATE( AdaptVars(i) % Var % Values( j ), STAT=AllocStat )
                IF( AllocStat /= 0 ) CALL Fatal(Caller,'Allocation error AdaptVars Values')

                IF( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                  k = SIZE( Solver % Variable % PrevValues, 2 )
                  ALLOCATE( AdaptVars(i) % Var % PrevValues( j, k ), STAT=AllocStat)
                  IF( AllocStat /= 0 ) CALL Fatal(Caller,'Allocation error for AdaptVars PrevValues')
                END IF
              END DO
            END IF
            
            CumTime = 0.0d0
            IF ( ddt < 0.0_dp .OR. ddt > AdaptiveMaxTimestep ) ddt = AdaptiveMaxTimestep
            
            s = sTime(1) - dt
            SmallestCount = 0
            DO WHILE( CumTime < dt-1.0d-12 )
              IF( .NOT. AdaptiveRough ) THEN
                ddt = MIN( dt - CumTime, ddt )
                IF( AdaptiveSmart ) THEN
                  ! If the next timestep will not get us home but the next one would
                  ! then split the timestep equally into two parts.
                  IF( dt - CumTime - ddt > 1.0d-12 ) THEN
                    CALL Info(Caller,'Splitted timestep into two equal parts',Level=12)
                    ddt = MIN( ddt, ( dt - CumTime ) / 2.0_dp )
                  END IF
                END IF
                
              END IF
                
               ! Store the initial values before the start of the step
               DO i=1,nSolvers
                 Solver => CurrentModel % Solvers(i)
                 IF ( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE
                 IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) CYCLE
                 AdaptVars(i) % Var % Values = Solver % Variable % Values
                 AdaptVars(i) % Var % Norm = Solver % Variable % Norm
                 IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                   AdaptVars(i) % Var % PrevValues = Solver % Variable % PrevValues
                 END IF
               END DO
               
               sTime(1) = s + CumTime + ddt
               sSize(1) = ddt

               ! Solve with full timestep
               CALL SolveEquations( CurrentModel, ddt, Transient, &
                   CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                   BeforeTime = .TRUE., AtTime = .TRUE., AfterTime = .FALSE.)           

               ! External adaptive error given 
               MaxErr = ListGetConstReal( CurrentModel % Simulation, &
                   'Adaptive Error Measure', GotIt )

               DO i=1,nSolvers
                 Solver => CurrentModel % Solvers(i)
                 IF ( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE 
                 IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) CYCLE
                 Solver % Variable % Values = AdaptVars(i) % Var % Values 
                 AdaptVars(i) % Norm = Solver % Variable % Norm
                 IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                   Solver % Variable % PrevValues = AdaptVars(i) % Var % PrevValues 
                 END IF
               END DO

               ! Test the error for half the timestep
               sStep(1) = ddt / 2
               sTime(1) = s + CumTime + ddt/2
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                   CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                   BeforeTime = .TRUE., AtTime = .TRUE., AfterTime = .FALSE.)           

               sTime(1) = s + CumTime + ddt
               CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                   CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                   BeforeTime = .TRUE., AtTime = .TRUE., AfterTime = .FALSE.)           

               MaxErr = ABS( MaxErr - ListGetConstReal( CurrentModel % Simulation, &
                   'Adaptive Error Measure', GotIt ) )

               ! If not measure given then use the maximum change in norm as the measure
               IF ( .NOT. GotIt ) THEN
                 MaxErr = 0.0d0
                 DO i=1,nSolvers
                   Solver => CurrentModel % Solvers(i)
                   IF ( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE
                   IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) CYCLE
                   IF ( AdaptVars(i) % norm /= Solver % Variable % Norm ) THEN
                     Maxerr = MAX(Maxerr,ABS(AdaptVars(i) % norm - Solver % Variable % Norm)/&
                         AdaptVars(i) % norm )
                   END IF
                 END DO
               END IF
               
               ! We have a success no need to redo this step
               IF ( MaxErr < AdaptiveLimit .OR. ddt <= AdaptiveMinTimestep ) THEN
                 CumTime = CumTime + ddt
                 RealTimestep = RealTimestep+1
                 IF ( SmallestCount >= AdaptiveKeepSmallest .OR. StepControl > 0 ) THEN
                   ddt = MIN( AdaptiveIncrease * ddt, AdaptiveMaxTimeStep )
                   StepControl   = 1
                   SmallestCount = 0
                 ELSE
                   StepControl   = 0
                   SmallestCount = SmallestCount + 1
                 END IF

                 ! Finally solve only the postprocessing solver after the timestep has been accepted
                 CALL SolveEquations( CurrentModel, ddt/2, Transient, &
                     CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                     BeforeTime = .FALSE., AtTime = .FALSE., AfterTime = .TRUE.)           
               ELSE
                 DO i=1,nSolvers
                   Solver => CurrentModel % Solvers(i)
                   IF ( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE
                   IF ( .NOT. ASSOCIATED( Solver % Variable % Values ) ) CYCLE
                   Solver % Variable % Norm = AdaptVars(i) % Var % Norm 
                   Solver % Variable % Values = AdaptVars(i) % Var % Values 
                   IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                     Solver % Variable % PrevValues = AdaptVars(i) % Var % PrevValues
                   END IF
                 END DO
                 ddt = AdaptiveDecrease * ddt 
                 StepControl = -1
               END IF

               WRITE(*,'(a,3e20.12)') 'Adaptive(cum,ddt,err): ', cumtime, ddt, maxerr
             END DO
            sSize(1) = dt
            sTime(1) = s + dt

          ELSE IF( DivergenceControl ) THEN
            ! This is still tentative 
            CALL Info(Caller,'Solving equations with divergence control',Level=6)
            
            CumTime = 0.0d0
            ddt = AdaptiveIncrease * ddt
            IF ( ddt < 0.0_dp .OR. ddt > AdaptiveMaxTimestep ) ddt = AdaptiveMaxTimestep            
            s = sTime(1) - dt
            !StepControl = 0
            
            DO WHILE( CumTime < dt-1.0d-12 )              

              IF( .NOT. AdaptiveRough ) THEN
                ddt = MIN( dt - CumTime, ddt )
                IF( AdaptiveSmart ) THEN
                  IF( dt - CumTime - ddt > 1.0d-12 ) THEN
                    CALL Info(Caller,'Splitted timestep into two equal parts',Level=12)
                    ddt = MIN( ddt, ( dt - CumTime ) / 2.0_dp )
                  END IF
                END IF
              END IF

              sTime(1) = s + CumTime + ddt
              sSize(1) = ddt

              CALL SolveEquations( CurrentModel, ddt, Transient, &
                  CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                  BeforeTime = .TRUE., AtTime = .TRUE., AfterTime = .FALSE.)

              HaveDivergence = .FALSE.
              DO i=1,nSolvers
                Solver => CurrentModel % Solvers(i) 
                IF( ASSOCIATED( Solver % Variable ) ) THEN
                  IF( Solver % Variable % NonlinConverged > 1 ) THEN
                    CALL Info(Caller,'Solver '//I2S(i)//' has diverged',Level=8)
                    HaveDivergence = .TRUE.
                    EXIT
                  END IF
                END IF
              END DO
              
              IF( .NOT. HaveDivergence ) THEN
                CALL Info(Caller,'No solver has diverged',Level=8)

                ! Finally solve only the postprocessing solver after the timestep has been accepted
                CALL SolveEquations( CurrentModel, ddt, Transient, &
                    CoupledMinIter, CoupledMaxIter, SteadyStateReached, RealTimestep, &
                    BeforeTime = .FALSE., AtTime = .FALSE., AfterTime = .TRUE.)           
                
                ! If step control was active for this interval then we can
                ! again start to increase timestep, otherwise not

                CumTime = CumTime + ddt
                RealTimestep = RealTimestep+1
                
                !IF ( StepControl > 0 ) THEN
                  ddt = AdaptiveIncrease * ddt
                  IF( ddt > AdaptiveMaxTimestep ) THEN
                    ddt = AdaptiveMaxTimestep
                    StepControl = 0
                  END IF
                !END IF
             ELSE
               IF( ddt < AdaptiveMinTimestep * (1+1.0e-8) ) THEN
                 CALL Fatal(Caller,'Could not find stable timestep above given minimum')
               END IF

               CALL Info(Caller,'Reducing timestep due to divergence problems!',Level=6)

               ddt = MAX( AdaptiveDecrease * ddt, AdaptiveMinTimestep ) 
               StepControl = 1

               CALL Info(Caller,'Reverting to previous timestep as initial guess',Level=8)
               DO i=1,nSolvers
                 Solver => CurrentModel % Solvers(i)
                 IF ( ASSOCIATED( Solver % Variable % Values ) ) THEN
                   IF( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                     Solver % Variable % Values = Solver % Variable % PrevValues(:,1)
                     Solver % Variable % Norm = Solver % Variable % PrevNorm 
                   END IF
                 END IF
               END DO
             END IF
            END DO
          ELSE
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
             DO i=1,nSolvers
               Solver => CurrentModel % Solvers(i)
               IF ( Solver % PROCEDURE == 0 ) CYCLE
               ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_SAVE)
               When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
               IF ( GotIt ) ExecThis = ( When == 'before saving') 
               IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
             END DO 

             ! Output file is used to Save the results for restart.
             ! Optionally we may save just the final stage which saves disk space and time.
             IF( .NOT. ListGetLogical( CurrentModel % Simulation,'Output File Final Only',GotIt) ) THEN               
               CALL SaveCurrent(Timestep)
             END IF

             CALL SaveToPost(TimeStep)
             LastSaved = .TRUE.

             DO i=1,nSolvers
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
         CALL ListPopNameSpace()
!------------------------------------------------------------------------------

         maxtime = ListGetCReal( CurrentModel % Simulation,'Real Time Max',GotIt)
         IF( GotIt .AND. RealTime() - RT0 > maxtime ) THEN
            CALL Info('MAIN','Reached allowed maximum real time, exiting...',Level=3)
            GOTO 100
         END IF

	 exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',GotIt)
	 IF( GotIt .AND. exitcond > 0.0_dp ) THEN
            CALL Info('MAIN','Found a positive exit condition, exiting...',Level=3)
            GOTO 100
         END IF

         IF( sFinish(1) > 0.0_dp ) THEN
           CALL Info('MAIN','Finishing condition "finish" found to be positive, exiting...',Level=3)
           GOTO 100
         END IF
           
!------------------------------------------------------------------------------

         IF ( SteadyStateReached .AND. .NOT. (Transient .OR. Scanning) ) THEN
            IF ( Timestep >= CoupledMinIter ) EXIT
         END IF

         ! if extra steps have been added need to check if the loop needs extended
         ! (used for calving algorithm)
         IF(Transient) THEN
            Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
            'Timestep Intervals', GotIt )
            IF ( .NOT.GotIt ) THEN
                CALL Fatal('ElmerSolver', 'Keyword > Timestep Intervals < MUST be ' //  &
                    'defined for transient and scanning simulations' )
            END IF
         END IF
         Timestep = Timestep + 1

         stepcount = 0
         DO i = 1, TimeIntervals
            stepcount = stepcount + Timesteps(i)
         END DO

!------------------------------------------------------------------------------
       END DO ! timestep within an iterval
!------------------------------------------------------------------------------
       
!------------------------------------------------------------------------------
     END DO ! timestep intervals, i.e. the simulation
!------------------------------------------------------------------------------

     BLOCK
       TYPE(Solver_t), POINTER :: iSolver

       DO i=1,CurrentModel % NumberOfSolvers 
         iSolver => CurrentModel % Solvers(i)
         IF( iSolver % NumberOfConstraintModes > 0 ) THEN
           IF( ListGetLogical( iSolver % Values,'Steady State Constraint Modes', Found ) ) THEN
             CALL FinalizeLumpedMatrix( iSolver )            
           END IF
         END IF
       END DO
     END BLOCK
     
100  CONTINUE

     CALL ListPopNamespace()

     DO i=1,nSolvers
        Solver => CurrentModel % Solvers(i)
        IF ( Solver % PROCEDURE == 0 ) CYCLE
        When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
        IF ( GotIt ) THEN
           IF ( When == 'after simulation' .OR. When == 'after all' ) THEN
              CALL SolverActivate( CurrentModel,Solver,dt,Transient )
              !IF( ASSOCIATED(Solver % Variable) ) THEN
              ! This construct seems to be for cases when we solve something "after all"
              ! that affects results elsewhere. Hence we set "LastSaved" to false even
              ! if it would be true before.                 
              !IF (ASSOCIATED(Solver % Variable % Values) ) LastSaved = .FALSE.
              !END IF
           END IF
        ELSE
           IF ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_ALL ) THEN
              CALL SolverActivate( CurrentModel,Solver,dt,Transient )
              !IF( ASSOCIATED(Solver % Variable) ) THEN
              !  IF (ASSOCIATED(Solver % Variable % Values) ) LastSaved = .FALSE.
              !END IF
           END IF
        END IF
     END DO

!------------------------------------------------------------------------------
!    Always save the last step to output
!-----------------------------------------------------------------------------
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
       CALL SaveToPost(TimeStep)
       
       IF( .NOT. ListGetLogical( CurrentModel % Simulation,'Output File Final Only',GotIt) ) THEN               
         CALL SaveCurrent(Timestep)
       END IF

       DO i=1,CurrentModel % NumberOfSolvers
         Solver => CurrentModel % Solvers(i)
         IF ( Solver % PROCEDURE == 0 ) CYCLE
         ExecThis = ( Solver % SolverExecWhen == SOLVER_EXEC_AFTER_SAVE)
         When = ListGetString( Solver % Values, 'Exec Solver', GotIt )
         IF ( GotIt ) ExecThis = ( When == 'after saving') 
         IF( ExecThis ) CALL SolverActivate( CurrentModel,Solver,dt,Transient )
       END DO
     ELSE IF( ListGetLogical( CurrentModel % Simulation,'Output File Final Only',GotIt) ) THEN               
       CALL SaveCurrent(Timestep)
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE ExecSimulation
!------------------------------------------------------------------------------


#ifdef HAVE_EXTOPTIM
!------------------------------------------------------------------------------
!> This is a handle for doing parametric simulation using minpack optimization
!> library. The function returns a set of function values. 
!------------------------------------------------------------------------------
   SUBROUTINE ExecSimulationFunVec(NoParam,Param,Fvec,iflag ) 
     INTEGER, INTENT(in) :: NoParam
     REAL(KIND=dp), INTENT(in) :: Param(NoParam)
     REAL(KIND=dp), INTENT(out) :: Fvec(NoParam)
     INTEGER, INTENT(inout) :: iflag

     INTEGER :: i, cnt, iSweep = 0
     LOGICAL :: Found
     
     iSweep = iSweep + 1

     CALL Info('ExecSimulationFunVec','Calling Elmer as a cost function: '//I2S(iSweep))
     
     IF(iSweep==1) THEN
       CONTINUE
     ELSE
       CALL ControlResetMesh(Control % Control, iSweep )            
     END IF

     ! Optionally reset the mesh if it has been modified
     CALL ControlResetMesh(Control % Control, iSweep )            

     ! Set parameters to be accessible to the MATC preprocessor when reading sif file.
     CALL SetRealParametersMATC(NoParam,Param)
     ! Reread the sif file for MATC changes to take effect
     Found = ReloadInputFile(CurrentModel,RewindFile=.TRUE.)

     ! Update the parameters also as coefficient as we don't know which one we are using
     CALL SetRealParametersKeywordCoeff(NoParam,Param,cnt)
     CALL Info('ExecSimulationFunVec','Set '//I2S(cnt)//&
         ' coefficients with parameter tags!',Level=10)

     CALL InitializeIntervals()
     CALL SetInitialConditions()

     CALL ExecSimulation( TimeIntervals, CoupledMinIter, &
         CoupledMaxIter, OutputIntervals, Transient, Scanning)

     DO i=1,NoParam     
       Fvec(i) = GetControlValue(CurrentModel % Mesh,CurrentModel % Control,i)
     END DO

     PRINT *,'Fvec:',iSweep, Fvec

   END SUBROUTINE ExecSimulationFunVec


!------------------------------------------------------------------------------
!> This is a handle for doing parametric simulation with optimizers that
!> except one single value of cost function. 
!------------------------------------------------------------------------------

   SUBROUTINE ExecSimulationFunCost(NoParam,Param,Cost) 
     INTEGER, INTENT(in) :: NoParam
     REAL(KIND=dp), INTENT(in) :: Param(NoParam)
     REAL(KIND=dp), INTENT(out) :: Cost

     INTEGER :: i, cnt, iSweep = 0
     LOGICAL :: Found
     
     iSweep = iSweep + 1

     CALL Info('ExecSimulationFunCost','Calling Elmer as a cost function: '//I2S(iSweep))
     
     IF(iSweep==1) THEN
       CONTINUE
     ELSE
       CALL ControlResetMesh(Control % Control, iSweep )            
     END IF

     ! Optionally reset the mesh if it has been modified
     CALL ControlResetMesh(Control % Control, iSweep )            

     ! Set parameters to be accessible to the MATC preprocessor when reading sif file.
     CALL SetRealParametersMATC(NoParam,Param)
     ! Reread the sif file for MATC changes to take effect
     Found = ReloadInputFile(CurrentModel,RewindFile=.TRUE.)

     ! Update the parameters also as coefficient as we don't know which one we are using
     CALL SetRealParametersKeywordCoeff(NoParam,Param,cnt)
     CALL Info('ExecSimulationFunCost','Set '//I2S(cnt)//&
         ' coefficients with parameter tags!',Level=10)

     CALL InitializeIntervals()
     CALL SetInitialConditions()

     CALL ExecSimulation( TimeIntervals, CoupledMinIter, &
         CoupledMaxIter, OutputIntervals, Transient, Scanning)

     CALL GetCostFunction(CurrentModel % Control,Cost,Found)
     IF(.NOT. Found ) THEN
       CALL Fatal('ExecSimulationFunCost','Could not find cost function!')
     END IF
     
     PRINT *,'Cost:',iSweep, Cost

   END SUBROUTINE ExecSimulationFunCost
!------------------------------------------------------------------------------
#endif   

!------------------------------------------------------------------------------
!> Saves current timestep to external files.
!------------------------------------------------------------------------------
  SUBROUTINE SaveCurrent( CurrentStep )
!------------------------------------------------------------------------------
    INTEGER :: i, j,k,l,n,m,q,CurrentStep,nlen,Time
    TYPE(Variable_t), POINTER :: Var, TimeVar
    LOGICAL :: EigAnal, GotIt, BinaryOutput, SaveAll, OutputActive, EveryTime
    TYPE(ValueList_t), POINTER :: vList
    TYPE(Solver_t), POINTER :: pSolver
    
    CALL Info('SaveCurrent','Saving information on current step',Level=20)
    
    ! There are currently global definitions that apply also for solver specific meshes
    vList => CurrentModel % Simulation       
    BinaryOutput = ListGetLogical( vList,'Binary Output',GotIt )      
    SaveAll = .NOT. ListGetLogical( vList,'Omit unchanged variables in output',GotIt )
    IF ( .NOT.GotIt ) SaveAll = .TRUE.
    

    Mesh => CurrentModel % Meshes
    DO WHILE( ASSOCIATED( Mesh ) ) 
      
      ! We want to be able to give the "output file" a different name for solver-specific
      ! meshes. Otherwise we might write data on the same file. 
      vList => CurrentModel % Simulation
      OutputActive = Mesh % OutputActive 
      m = Mesh % SolverId
      IF( m > 0 ) THEN
        pSolver => CurrentModel % Solvers(m)
        IF( ListCheckPresent( pSolver % Values,'Output File') ) THEN
          vList => pSolver % Values
          OutputActive = ListGetLogical( vList,'Mesh Output',GotIt)
          IF(.NOT. GotIt) OutputActive = Mesh % OutputActive
        END IF
      END IF

      IF(.NOT. OutputActive ) THEN
        Mesh => Mesh % Next
        CYCLE
      END IF
      
      OutputFile = ListGetString( vList,'Output File',GotIt)                   
      IF(GotIt) THEN
        i = INDEX( OutputFile,'/')
        IF( i > 0 ) THEN
          CALL Warn('SaveCurrent','> Output File < for restart should not include directory: '&
              //TRIM(OutputFile))
        END IF

        ! This was added for needs of calving where remeshing is applied and can go wrong but
        ! we want to study the results. 
        EveryTime = ListGetLogical( vList,'Output File Each Timestep',GotIt)
        IF(EveryTime) THEN
          TimeVar => VariableGet( CurrentModel % Variables, 'Timestep' )
          Time = INT(TimeVar % Values(1))

          OutputFile = OutputFile // '.' // i2s(Time)

          ! set saves to zero. This will insure new save file even if remeshing fails
          Mesh % SavesDone = 0
        END IF

        !IF ( ParEnv % PEs > 1 ) THEN
        !  DO i=1,MAX_NAME_LEN
        !    IF ( OutputFile(i:i) == ' ' ) EXIT
        !  END DO
        !  OutputFile(i:i) = '.'
        !  WRITE( OutputFile(i+1:), '(a)' ) i2s(ParEnv % MyPE)
        !END IF

        ! Always write the output with respect to mesh file
        nlen = LEN_TRIM(Mesh % Name )
        IF ( nlen > 0 ) THEN
          OutputName = Mesh % Name(1:nlen) // '/' // TRIM(OutputFile)
        END IF
        
        EigAnal = .FALSE.
        DO i=1,CurrentModel % NumberOfSolvers
          pSolver => CurrentModel % Solvers(i)
          IF ( ASSOCIATED( pSolver % Mesh, Mesh ) ) THEN
            EigAnal = ListGetLogical( pSolver % Values,'Eigen Analysis', GotIt ) .OR. &
                ListGetLogical( pSolver % Values,'Harmonic Analysis', GotIt )

            IF ( EigAnal ) THEN
              Var => pSolver % Variable
              IF ( ASSOCIATED(Var % EigenValues) ) THEN
                IF ( TotalTimesteps == 1 ) THEN
                  DO j=1,pSolver % NOFEigenValues
                    IF ( pSolver % Matrix % COMPLEX ) THEN

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
                        j, sTime(1), BinaryOutput, SaveAll, vList = vList )
                  END DO
                ELSE
                  j = MIN( CurrentStep, SIZE( Var % EigenVectors,1 ) )
                  IF ( pSolver % Matrix % COMPLEX ) THEN
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
                      CurrentStep, sTime(1), BinaryOutput, SaveAll, vList = vList )
                END IF
                Var % Values = 0.0d0
              END IF
            END IF
          END IF
        END DO
          
        IF ( .NOT. EigAnal ) THEN
          SavedSteps = SaveResult( OutputName,Mesh, NINT(sStep(1)), &
              sTime(1), BinaryOutput, SaveAll, vList = vList )
        END IF
      ELSE
        Mesh % SavesDone = Mesh % SavesDone+1
      END IF
      Mesh => Mesh % Next
    END DO
    
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
    INTEGER :: i, j,k,l,n,q,CurrentStep,nlen,nlen2,timesteps,SavedEigenValues
    CHARACTER(LEN=MAX_NAME_LEN) :: Simul, SaveWhich
    CHARACTER(MAX_NAME_LEN) :: OutputDirectory
    TYPE(Solver_t), POINTER :: pSolver
    
    Simul = ListGetString( CurrentModel % Simulation,'Simulation Type' )

    OutputFile = ListGetString( CurrentModel % Simulation,'Output File',GotIt )

    IF ( Gotit ) THEN
      IF ( ParEnv % PEs > 1 ) THEN
        OutputFile = OutputFile // '.' // i2s(ParEnv % MyPE)
      END IF
    END IF
    
    PostFile = ListGetString( CurrentModel % Simulation,'Post File',GotIt )
    IF( .NOT. GotIt ) RETURN

    IF ( ParEnv % PEs > 1 ) THEN
      PostFile = PostFile // '.' // i2s(ParEnv % MyPE)
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
        
        nlen2 = LEN_TRIM(OutputPath)
        IF(nlen2 == 1) THEN
          IF(OutputPath(1:1) == '.') nlen2 = 0
        END IF
        
        ! If "Results Directory" is given (nlen2>0) we want to give that
        ! priority over mesh directory. 
        IF ( FileNameQualified(PostFile) .OR. nlen2 > 0 .OR. nlen==0 ) THEN
          PostName = PostFile
        ELSE
          PostName = Mesh % Name(1:nlen)//'/'//TRIM(PostFile)
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
          pSolver => CurrentModel % Solvers(i)
          IF (ASSOCIATED(pSolver % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( pSolver % Values, 'Eigen Analysis', GotIt ) .OR. &
                ListGetLogical( pSolver % Values, 'Harmonic Analysis', GotIt )            
            IF ( EigAnal ) timesteps = MAX( timesteps, pSolver % NOFEigenValues )
          END IF
        END DO

        DO i=1,CurrentModel % NumberOfSolvers
          pSolver => CurrentModel % Solvers(i)
          IF (ASSOCIATED(pSolver % Mesh, Mesh)) THEN
            EigAnal = ListGetLogical( pSolver % Values, 'Eigen Analysis', GotIt ) .OR. &
                ListGetLogical( pSolver % Values, 'Harmonic Analysis', GotIt )
            
            IF ( EigAnal ) THEN
              SaveWhich = ListGetString( pSolver % Values, &
                  'Eigen and Harmonic Solution Output', Found )              
              SavedEigenValues = pSolver % NOFEigenValues

              DO j=1, SavedEigenValues
                Var => Mesh % Variables
                DO WHILE(ASSOCIATED(Var))
                  IF ( .NOT. ASSOCIATED(Var % EigenValues) ) THEN
                    Var => Var % Next
                    CYCLE
                  END IF
                  
                  IF ( pSolver % Matrix % COMPLEX ) THEN
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
                      pSolver % NOFEigenValues, .TRUE. )
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
