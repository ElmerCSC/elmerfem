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

     CHARACTER(LEN=MAX_NAME_LEN) :: ModelName, eq, ExecCommand, ExtrudedMeshName
     CHARACTER(LEN=MAX_STRING_LEN) :: OutputFile, PostFile, RestartFile, &
                OutputName=' ',PostName=' ', When, OptionString

     TYPE(Variable_t), POINTER :: Var
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Solver_t), POINTER :: Solver

     REAL(KIND=dp) :: CT0,RT0,tt

     LOGICAL :: FirstLoad = .TRUE., FirstTime=.TRUE., Found
     LOGICAL :: Silent, Version, GotModelName, FinishEarly

     INTEGER :: iargc, NoArgs
     INTEGER :: iostat, iSweep = 1, OptimIters
     
     INTEGER :: MeshIndex
     TYPE(Mesh_t), POINTER :: ExtrudedMesh

     TYPE(Model_t), POINTER, SAVE :: Control
     CHARACTER(LEN=MAX_NAME_LEN) :: MeshDir, MeshName
     LOGICAL :: DoControl, GotParams
     INTEGER :: nr,ni
     REAL(KIND=dp), ALLOCATABLE :: rpar(:)
     INTEGER, ALLOCATABLE :: ipar(:)
     
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
             CALL Info('MAIN','Read '//TRIM(I2S(nr))//' real parameters from command line!')
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
             CALL Info('MAIN','Read '//TRIM(I2S(ni))//' integer parameters from command line!')
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
               TRIM(i2s(ParEnv % PEs)) // ' tasks.')
         ELSE
           CALL Info('MAIN', ' Running one task without MPI parallelization.',Level=10)
         END IF
         
         ! Print out number of threads in use
         IF ( nthreads > 1 ) THEN
           CALL Info('MAIN', ' Running in parallel with ' // &
               TRIM(i2s(nthreads)) // ' threads per task.')
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
#ifdef HAVE_MKL
         CALL Info( 'MAIN', ' Intel MKL linked in.' )
#endif
#ifdef HAVE_LUA
         CALL Info( 'MAIN', ' Lua interpreter linked in.' )
#endif
#ifdef HAVE_ZOLTAN
         CALL Info( 'MAIN', ' Zoltan library linked in.' )
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
         CALL Fatal( 'ElmerSolver', 'Unable to find ELMERSOLVER_STARTINFO, can not execute.' )
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
         CALL Fatal( 'ElmerSolver', 'Unable to find input file [' // &
             TRIM(Modelname) // '], can not execute.' )
       END IF
       ALLOCATE( Control )          
       ! Read only the "Run Control" section of the sif file.
       CALL LoadInputFile( Control,InFileUnit,ModelName,MeshDir,MeshName, &
           .FALSE., .TRUE., ControlOnly = .TRUE.)
       DoControl =  ASSOCIATED( Control % Control )
       IF( DoControl ) THEN
         CALL Info('ElmerSolver','Run Control section active!')
         OptimIters = ListGetInteger( Control % Control,'Run Control Iterations', Found )
         IF(.NOT. Found) OptimIters = 1              

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

         ! Here we read the whole model including command file and detault mesh file
         !---------------------------------------------------------------------------------
         OPEN( Unit=InFileUnit, Action='Read',File=ModelName,Status='OLD',IOSTAT=iostat)         
         IF( iostat /= 0 ) THEN
           CALL Fatal( 'ElmerSolver', 'Unable to find input file [' // &
               TRIM(Modelname) // '], can not execute.' )
         END IF
         

         CurrentModel => LoadModel(ModelName,.FALSE.,ParEnv % PEs,ParEnv % MyPE,MeshIndex)
         IF(.NOT.ASSOCIATED(CurrentModel)) EXIT

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
               File=ModelName,Status='OLD',IOSTAT=iostat)
           IF( iostat /= 0 ) THEN
             CALL Fatal( 'ElmerSolver', 'Unable to find input file [' // &
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
!      Figure out what (flow,heat,stress,...) should be computed, and get
!      memory for the dofs
!------------------------------------------------------------------------------
       CALL AddSolvers()

!------------------------------------------------------------------------------
!      Time integration and/or steady state steps
!------------------------------------------------------------------------------
       CALL InitializeIntervals()

!------------------------------------------------------------------------------
!      Add coordinates and simulation time to list of variables so that
!      coordinate dependent parameter computing routines can ask for
!      them...
!------------------------------------------------------------------------------
       IF ( FirstLoad ) CALL AddMeshCoordinatesAndTime()

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

       ! This sets optionally some internal parameters for doing scanning
       ! over a parameter space / optimization. 
       !-----------------------------------------------------------------
       DO iSweep = 1, OptimIters
         sSweep = 1.0_dp * iSweep
         ! If there are no parameters this does nothing                  
         IF( DoControl ) THEN
           CALL ControlResetMesh(Control % Control, iSweep )            
           IF( iSweep > 1 ) THEN
             CALL ControlParameters(Control % Control,iSweep,&
                 GotParams,FinishEarly)           
             IF( FinishEarly ) EXIT
             Found = ReloadInputFile(CurrentModel,RewindFile=.TRUE.)
             CALL InitializeIntervals()
           END IF

           ! This is another calling slot as here we have formed the model structure and
           ! may toggle with the keyword coefficients. 
           CALL ControlParameters(Control % Control,iSweep,&
               GotParams,FinishEarly,SetCoeffs=.TRUE.)

           IF( iSweep > 1 ) THEN
             IF( ListGetLogical( Control % Control,'Reset Initial Conditions',Found ) ) THEN
               CALL SetInitialConditions()
             END IF
           END IF
           
         END IF
           
         !------------------------------------------------------------------------------
         ! Here we actually perform the simulation: ExecSimulation does it all ....
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

         ! This evaluates the cost function and saves the results of control
         IF( DoControl ) THEN
           CALL ControlParameters(Control % Control, &
               iSweep,GotParams,FinishEarly,.TRUE.)
         END IF
       END DO
       
       ! Comparison to reference is done to enable consistency test underc ctest.
       !-------------------------------------------------------------------------
       CALL CompareToReferenceSolution( )

       IF ( Initialize >= 2 ) EXIT
     END DO

 
     CALL CompareToReferenceSolution( Finalize = .TRUE. )


#ifdef DEVEL_LISTCOUNTER
     CALL Info('ElmerSolver','Reporting list counters for code optimization purposes only!')
     CALL Info('ElmerSolver','If you get these lines with production code undefine > DEVEL_LISTCOUNTER < !')
     CALL ReportListCounters( CurrentModel )
#endif
     

     
!------------------------------------------------------------------------------
!    THIS IS THE END (...,at last, the end, my friend,...)
!------------------------------------------------------------------------------
     IF ( Initialize /= 1 ) CALL Info( 'ElmerSolver', '*** Elmer Solver: ALL DONE ***',Level=3 )

     ! This may be used to study problems at the finish
     IF( ListGetLogical( CurrentModel % Simulation,'Dirty Finish', GotIt ) ) THEN
       CALL Info('ElmerSolver','Skipping freeing of the Model structure',Level=4)
       RETURN
     END IF
     
     IF ( Initialize <= 0 ) CALL FreeModel(CurrentModel)

#ifdef HAVE_TRILINOS
  CALL TrilinosCleanup()
#endif

     IF ( FirstTime ) CALL ParallelFinalize()
     FirstTime = .FALSE.

     CALL Info('ElmerSolver','The end',Level=3)

     RETURN

   CONTAINS 


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
       
       CALL Info('ElmerSolver','Random seed initialized to: '//TRIM(I2S(i)),Level=10)
     END SUBROUTINE InitializeRandomSeed

     
     ! Optionally create extruded mesh on-the-fly.
     !--------------------------------------------------------------------
     SUBROUTINE CreateExtrudedMesh()

       INTEGER :: ExtrudeLayers
       
       ExtrudeLayers = GetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found) - 1 
       IF( .NOT. Found ) THEN
         ExtrudeLayers = GetInteger(CurrentModel % Simulation,'Extruded Mesh Layers',Found)
       END IF
       IF(.NOT. Found ) RETURN

       IF(ExtrudeLayers < 2) THEN
         CALL Fatal('ElmerSolver','There must be at least two "Extruded Mesh Layers"!')
       END IF

       ExtrudedMeshName = GetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
       IF (Found) THEN
         ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, ExtrudeLayers-1, ExtrudedMeshName)
       ELSE
         ExtrudedMesh => MeshExtrude(CurrentModel % Meshes, ExtrudeLayers-1)
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

       CALL Info('ElmerSolver','Creating geometric timestepping strategy',Level=6)
       
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
           CALL Fatal('ElmerSolver', 'Keyword > Timestep Intervals < MUST be ' //  &
               'defined for transient and scanning simulations' )
         END IF
         
         TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation, &
             'Timestep Sizes', GotIt )
         IF ( .NOT.GotIt ) THEN
           IF( Scanning .OR. ListCheckPresent( CurrentModel % Simulation,'Timestep Size') ) THEN
             ALLOCATE(TimestepSizes(SIZE(Timesteps),1))
             TimestepSizes = 1.0_dp
           ELSE
             CALL Fatal( 'ElmerSolver', 'Keyword [Timestep Sizes] MUST be ' //  &
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
           CALL Fatal('ElmerSolver','> Output Intervals < should have the same size as > Timestep Intervals < !')
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
       
       CALL Info('ElmerSolver','Number of timesteps to be saved: '//TRIM(I2S(TotalTimesteps)))
       
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
               'PASSED all '//TRIM(I2S(TestCount))//' tests!',Level=3)
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
               //TRIM(I2S(i))//' not associated, cannot compare')
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
           //TRIM(I2S(n)),Level=8)
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
     

     !------------------------------------------------------------------------------
     FUNCTION FindSolverByProcName(Model,ProcName) RESULT (solver_id)
       IMPLICIT NONE

       TYPE(Model_t), POINTER :: Model
       CHARACTER(*) :: ProcName
       INTEGER :: solver_id
       
       LOGICAL :: Found
       INTEGER :: i,j
       TYPE(Solver_t), POINTER :: pSolver
       CHARACTER(LEN=MAX_NAME_LEN) :: str

       solver_id = 0       
       Found = .FALSE.

       !PRINT *,'procname:',TRIM(ProcName)
       
       DO i=1, Model % NumberOfSolvers
         pSolver => CurrentModel % Solvers(i)
         str = ListGetString(pSolver % Values,'Procedure',Found)
         IF(.NOT. Found) CYCLE

         !PRINT *,'str:',i,TRIM(str)
         
         j = INDEX(str,ProcName)         
         IF( j > 0 ) THEN
           solver_id = i
           EXIT
         END IF
       END DO

       !PRINT *,'j:',j,solver_id
       
     END FUNCTION FindSolverByProcName
     !------------------------------------------------------------------------------


     
     ! This is a dirty hack that adds an instance of ResultOutputSolver to the list of Solvers.
     ! The idea is that it is much easier for the end user to take into use the vtu output this way.
     ! The solver itself has limited set of parameters needed and is therefore approapriate for this
     ! kind of hack. It can of course be also added as a regular solver also.
     !----------------------------------------------------------------------------------------------
     SUBROUTINE AddVtuOutputSolverHack()     
       TYPE(Solver_t), POINTER :: pSolver
       CHARACTER(LEN=MAX_NAME_LEN) :: str
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
         CALL Warn('AddVtuOutputSolverHack','ResultOutputSolver instance already exists, doing nothing!')
         RETURN
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
       CALL ListAddString(Params,'Exec Solver','after saving')
       CALL ListAddLogical(Params,'Save Geometry IDs',.TRUE.)
       CALL ListAddLogical(Params,'Check Simulation Keywords',.TRUE.)

       ! Add a few often needed keywords also if they are given in simulation section
       CALL ListCopyPrefixedKeywords( Simu, Params, 'vtu:' )

       CALL Info('AddVtuOutputSolverHack','Finished appending VTU output solver',Level=12)
       
     END SUBROUTINE AddVtuOutputSolverHack


     ! This is a dirty hack that adds an instance of SaveScalars to the list of Solvers.
     ! The idea is that it is much easier for the end user to add a basic instance.
     !----------------------------------------------------------------------------------------------
     SUBROUTINE AddSaveScalarsHack()     
       TYPE(Solver_t), POINTER :: ABC(:), PSolver
       CHARACTER(LEN=MAX_NAME_LEN) :: str
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
       CALL VariableAdd( Mesh % Variables, Mesh, &
             Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )

       CALL VariableAdd(Mesh % Variables,Mesh, &
             Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )

       CALL VariableAdd(Mesh % Variables,Mesh, &
             Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )

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
           
           CALL Info('AddMeshCoordinatesAndTime','Adding partitioning also as a field')
           
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
               IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
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
     INTEGER :: DOFs,i,i2,j,k,k1,k2,l,n,n2,m,nsize
     CHARACTER(LEN=MAX_NAME_LEN) :: str, VarName
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

       IF( InfoActive( 20 ) ) THEN
         PRINT *,'InitCond mesh:',TRIM(Mesh % Name), Mesh % MeshDim, Mesh % NumberOfNodes 
         Var => Mesh % Variables
         DO WHILE( ASSOCIATED(Var) ) 
           IF( ListCheckPresentAnyIC( CurrentModel, Var % Name ) ) THEN
             PRINT *,'InitCond pre range:',TRIM(Var % Name),MINVAL(Var % Values),MAXVAL( Var % Values)
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
                     IF ( k1>0 ) Var % PrevValues(k1,2) = Work(k)
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

           VarName = Var % Name 
           Solver => Var % Solver
           IF ( .NOT. ASSOCIATED(Solver) ) Solver => CurrentModel % Solver
           
           IF( AnyNameSpace ) THEN
             str = ListGetString( Solver % Values, 'Namespace', NamespaceFound )
             IF (NamespaceFound) CALL ListPushNamespace(TRIM(str))
           END IF
           
           VarOrder = -1
           DO VarOrder = 0, 2
             IF( VarOrder == 0 ) THEN
               VarName = Var % Name
             ELSE IF( VarOrder == 1 ) THEN
               VarName = TRIM( Var % Name )//' Velocity'
             ELSE IF( VarOrder == 2 ) THEN
               VarName = TRIM( Var % Name )//' Acceleration'
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
                     
                     CALL Info('InitCond','Total number of new IP dofs: '//TRIM(I2S(nsize)))
                     
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
             PRINT *,'InitCond post range:',TRIM(Var % Name),&
                 MINVAL(Var % Values),MAXVAL( Var % Values),SUM(Var % Values)/SIZE(Var % Values)
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
     LOGICAL :: CheckMesh, DoMesh
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
             OutputName = TRIM(RestartFile)
           END IF
                                 
           IF ( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh ) &
               OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))
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
       IF( GotIt ) RestartFile = TRIM(RestartFile)//'_'//TRIM(I2S(k))//'nc'
              
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
         IF ( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh ) &
           OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))
         
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
     LOGICAL :: AdaptiveTime = .TRUE., AdaptiveRough, AdaptiveSmart, Found
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
     
     !$OMP PARALLEL
     IF(.NOT.GaussPointsInitialized()) CALL GaussPointsInit()
     !$OMP END PARALLEL

     nSolvers = CurrentModel % NumberOfSolvers
     DO i=1,nSolvers
        Solver => CurrentModel % Solvers(i)
        IF ( Solver % PROCEDURE==0 ) CYCLE
        IF ( Solver % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL ) THEN
          ! solver to be called prior to time looping can never be transient
          dt = 1.0_dp
          CALL SolverActivate( CurrentModel,Solver,dt,.FALSE. )
        END IF
     END DO

     IF( ListGetLogical( CurrentModel % Simulation,'Calculate Mesh Pieces',Found ) ) THEN
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
     ParallelSlices = ListGetLogical( CurrentModel % Simulation,'Parallel Slices',GotIt ) &
         .AND. ( ParEnv % PEs > 1 )

     IF( ParallelTime .OR. ParallelSlices ) THEN
       IF(.NOT. ListGetLogical( CurrentModel % Simulation,'Single Mesh',GotIt ) ) THEN
         CALL Fatal('ExecSimulation','Parallel time and slices only available with "Single Mesh"')
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
           CALL Fatal('ExecSimulation','"Number Of Slices" cannot be be larger than #np')
         END IF
       ELSE
         IF( ParEnv % PEs == 1 ) THEN
           CALL ListAddInteger( CurrentModel % Simulation,'Number Of Slices',nSlices)
         ELSE
           CALL Fatal('ExecSimulation','We need "Number Of Slices" with parallel timestepping')
         END IF
       END IF
       IF( MODULO( ParEnv % PEs, nSlices ) /= 0 ) THEN
         CALL Fatal('ExecSimulation','For hybrid parallellism #np must be divisible with "Number of Slices"')
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
       CALL Info('ExecSimulation','Setting one time sector for each partition!')
     ELSE IF( ParallelSlices ) THEN
       nSlices = ParEnv % PEs
       iSlice = ParEnv % MyPe
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Times',nTimes )
       CALL ListAddInteger( CurrentModel % Simulation,'Number Of Slices',nSlices )
       CALL Info('ExecSimulation','Setting one slice for each partition!')
     END IF

     IF( nTimes > 1 ) THEN
       DO i=1,SIZE(Timesteps,1)
         IF( MODULO( Timesteps(i), nTimes ) /= 0 ) THEN
           CALL Fatal('ExecSimulation','"Timestep Intervals" should be divisible by nTimes: '//TRIM(I2S(nTimes)))
         END IF
         Timesteps(i) = Timesteps(i) / nTimes
       END DO
       CALL Info('ExecSimulation','Divided timestep intervals equally for each partition!',Level=4)
     END IF

     IsPeriodic = ListCheckPrefix( CurrentModel % Simulation,'Periodic Time') .OR. &
         ListCheckPresent( CurrentModel % Simulation,'Time Period' )
     
     nPeriodic = ListGetInteger( CurrentModel % Simulation,'Periodic Timesteps',GotIt )
     IF( ParallelTime ) THEN
       IF( nPeriodic <= 0 ) THEN
         CALL Fatal('ExecSimulation','Parallel timestepping requires "Periodic Timesteps"')
       END IF
       IF( MODULO( nPeriodic, nTimes ) /= 0 ) THEN
         CALL Fatal('ExecSimulation','For parallel timestepping "Periodic Timesteps" must be divisible by #np')
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
              
       DO timestep = 1,Timesteps(interval)
         
         cum_Timestep = cum_Timestep + 1
         sStep(1) = cum_Timestep

         IF ( GetNamespaceCheck() ) THEN
           IF( Scanning ) THEN
             CALL ListPushNamespace('scan '//TRIM(i2s(cum_Timestep))//':')
           ELSE IF ( Transient ) THEN
             CALL ListPushNamespace('time '//TRIM(i2s(cum_Timestep))//':')
           ELSE
             CALL ListPushNamespace('steady '//TRIM(i2s(cum_Timestep))//':')
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
             CALL Warn('ExecSimulation','Obsolete keyword > Timestep Function < , use > Timestep Size < instead')
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
               CALL Fatal('ExecSimulation', Message)
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
             CALL Fatal('ElmerSolver','Adaptive Time Error must be given for ' // &
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
         sTime(1) = sTime(1) + dt

         IF( nPeriodic > 0 ) THEN
           IF( ParallelTime ) THEN
             timePeriod = nTimes * nPeriodic * dt                        
             IF( cum_Timestep == 1 ) THEN
               sTime(1) = sTime(1) + iTime * nPeriodic * dt
             ELSE IF( MODULO( cum_Timestep, nPeriodic ) == 1 ) THEN
               CALL Info('ExecSimulation','Making jump in time-parallel scheme!')
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
             CALL Info('ExecSimulation',Message,Level=6)
             WRITE(Message,'(A,ES12.3)') '"Rotor Velo" set to value: ',sAngleVelo(1)
             CALL Info('ExecSimulation',Message,Level=6)
           END BLOCK
         END IF       
         
!------------------------------------------------------------------------------

         BLOCK
           TYPE(Mesh_t), POINTER :: Mesh
           REAL(KIND=dp) :: MeshR
           CHARACTER(LEN=MAX_NAME_LEN) :: MeshStr
           
           IF( ListCheckPresent( GetSimulation(), 'Mesh Name Index') ) THEN
             ! we cannot have mesh depend on "time" or "timestep" if they are not available as
             ! variables. 
             Mesh => CurrentModel % Meshes             
             CALL VariableAdd( Mesh % Variables, Mesh, Name='Time',DOFs=1, Values=sTime )
             CALL VariableAdd( Mesh % Variables, Mesh, Name='Timestep', DOFs=1, Values=sStep )

             MeshR = GetCREal( GetSimulation(), 'Mesh Name Index',gotIt )            
             i = NINT( MeshR )

             IF( i > 0 .AND. i /= PrevMeshI ) THEN                             
               MeshStr = ListGetString( GetSimulation(),'Mesh Name '//TRIM(I2S(i)),GotIt)
               IF( GotIt ) THEN
                 CALL Info('ExecSimulation','Swapping mesh to: '//TRIM(MeshStr),Level=5)
               ELSE
                 CALL Fatal('ExecSimulation','Could not find >Mesh Name '//TRIM(I2S(i))//'<')
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
             WRITE( Message,'(A,ES12.3)') 'Time: '//TRIM(i2s(cum_Timestep))//'/'// &
                   TRIM(i2s(stepcount))//':', sTime(1)
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

         IF ( Transient .AND. AdaptiveTime ) THEN 

            IF(.NOT. ALLOCATED( AdaptVars ) ) THEN
              ALLOCATE( AdaptVars( nSolvers ), STAT = AllocStat )
              IF( AllocStat /= 0 ) CALL Fatal('ExecSimulation','Allocation error for AdaptVars')
              
              DO i=1,nSolvers
                Solver => CurrentModel % Solvers(i)

                NULLIFY( AdaptVars(i) % Var % Values )
                NULLIFY( AdaptVars(i) % Var % PrevValues )

                IF( .NOT. ASSOCIATED( Solver % Variable ) ) CYCLE
                IF( .NOT. ASSOCIATED( Solver % Variable  % Values ) ) CYCLE
                CALL Info('ExecSimulation','Allocating adaptive work space for: '//TRIM(I2S(i)),Level=12)
                j = SIZE( Solver % Variable % Values )
                ALLOCATE( AdaptVars(i) % Var % Values( j ), STAT=AllocStat )
                IF( AllocStat /= 0 ) CALL Fatal('ExecSimulation','Allocation error AdaptVars Values')

                IF( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
                  k = SIZE( Solver % Variable % PrevValues, 2 )
                  ALLOCATE( AdaptVars(i) % Var % PrevValues( j, k ), STAT=AllocStat)
                  IF( AllocStat /= 0 ) CALL Fatal('ExecSimulation','Allocation error for AdaptVars PrevValues')
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
                    CALL Info('ExecSimulation','Splitted timestep into two equal parts',Level=12)
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
            CALL Info('ExecSimulation','Solving equations with divergence control',Level=6)
            
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
                    CALL Info('ExecSimulation','Splitted timestep into two equal parts',Level=12)
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
                    CALL Info('ExecSimulation','Solver '//TRIM(I2S(i))//' has diverged',Level=8)
                    HaveDivergence = .TRUE.
                    EXIT
                  END IF
                END IF
              END DO
              
              IF( .NOT. HaveDivergence ) THEN
                CALL Info('ExecSimulation','No solver has diverged',Level=8)

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
                 CALL Fatal('ExecSimulation','Could not find stable timestep above given minimum')
               END IF

               CALL Info('ExecSimulation','Reducing timestep due to divergence problems!',Level=6)

               ddt = MAX( AdaptiveDecrease * ddt, AdaptiveMinTimestep ) 
               StepControl = 1

               CALL Info('ExecSimulation','Reverting to previous timestep as initial guess',Level=8)
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
            CALL Info('ElmerSolver','Reached allowed maximum real time, exiting...',Level=3)
            GOTO 100
         END IF

	 exitcond = ListGetCReal( CurrentModel % Simulation,'Exit Condition',GotIt)
	 IF( GotIt .AND. exitcond > 0.0_dp ) THEN
            CALL Info('ElmerSolver','Found a positive exit condition, exiting...',Level=3)
            GOTO 100
         END IF

         IF( sFinish(1) > 0.0_dp ) THEN
           CALL Info('ElmerSolver','Finishing condition "finish" found to be positive, exiting...',Level=3)
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

100  CONTINUE

     CALL ListPopNamespace()

     DO i=1,nSolvers
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
      i = INDEX( OutputFile,'/')
      IF( i > 0 ) THEN
        CALL Warn('SaveCurrent','> Output File < for restart should not include directory: '&
            //TRIM(OutputFile))
      END IF
      
      !IF ( ParEnv % PEs > 1 ) THEN
      !  DO i=1,MAX_NAME_LEN
      !    IF ( OutputFile(i:i) == ' ' ) EXIT
      !  END DO
      !  OutputFile(i:i) = '.'
      !  WRITE( OutputFile(i+1:), '(a)' ) TRIM(i2s(ParEnv % MyPE))
      !END IF
      
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
! We want to separate saving of ElmerPost file and Result file.
!    CALL SaveToPost(CurrentStep)
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
    
    Simul = ListGetString( CurrentModel % Simulation,'Simulation Type' )

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
