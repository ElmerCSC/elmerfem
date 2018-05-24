!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Joe Todd
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *
! *
! *****************************************************************************

!Routines for dealing with the Remeshing of the 3D calving model

SUBROUTINE CheckFlowConvergence( Model, Solver, dt, Transient )

  USE CalvingGeometry

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !-------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t) :: RemeshSolver
  TYPE(Variable_t), POINTER :: FlowVar, TimeVar
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Parallel, Found, CheckFlowDiverge=.TRUE., CheckFlowMax, FirstTime=.TRUE.,&
       NSDiverge, NSFail, NSTooFast
  REAL(KIND=dp) :: SaveNewtonTol, MaxNSDiverge, MaxNSValue, FirstMaxNSValue, FlowMax,&
  SaveFlowMax, Mag, NSChange, SaveDt, SaveRelax,SaveMeshMinLC,SaveMeshRmLC,SaveMeshRmThresh
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
  INTEGER :: i,j,SaveNewtonIter,Num, ierr, FailCount
  CHARACTER(MAX_NAME_LEN) :: FlowVarName, SolverName, EqName, RemeshEqName

  SAVE ::SaveNewtonTol, SaveNewtonIter, SaveFlowMax, SaveDt, FirstTime, FailCount,&
       SaveRelax,SaveMeshMinLC,SaveMeshRmLC,SaveMeshRmThresh

  Mesh => Solver % Mesh
  SolverName = 'CheckFlowConvergence'
  Params => Solver % Values
  Parallel = (ParEnv % PEs > 1)

  FlowVarName = ListGetString(Params,'Flow Solver Name',Found)
  IF(.NOT. Found) FlowVarName = "Flow Solution"
  FlowVar => VariableGet(Mesh % Variables, FlowVarName, .TRUE., UnfoundFatal=.TRUE.)

  RemeshEqName = ListGetString(Params,'Remesh Equation Name',Found)
  IF(.NOT. Found) RemeshEqName = "remesh"

  !Get a handle to the remesh solver
  Found = .FALSE.
  DO j=1,Model % NumberOfSolvers
    IF(ListGetString(Model % Solvers(j) % Values, "Equation") == RemeshEqName) THEN
      RemeshSolver = Model % Solvers(j)
      Found = .TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT. Found) CALL Fatal(SolverName, "Failed to get handle to Remesh Solver.")

  IF(FirstTime) THEN

    FailCount = 0
    TimestepSizes => ListGetConstRealArray( Model % Simulation, &
         'Timestep Sizes', Found, UnfoundFatal=.TRUE.)
    IF(SIZE(TimestepSizes,1) > 1) CALL Fatal(SolverName,&
         "Calving solver requires a single constant 'Timestep Sizes'")
    SaveDt = TimestepSizes(1,1)

    SaveNewtonTol = ListGetConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", Found)
    IF(.NOT. Found) SaveNewtonTol = 1.0E-3
    SaveNewtonIter = ListGetInteger(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", Found)
    IF(.NOT. Found) SaveNewtonIter = 15

    SaveRelax = ListGetConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Relaxation Factor", Found)
    IF(.NOT. Found) SaveRelax = 1.0

    SaveMeshMinLC = ListGetConstReal(RemeshSolver % Values, &
         "Remesh Min Characteristic Length", Found, UnfoundFatal=.TRUE.)

    SaveMeshRmLC = ListGetConstReal(RemeshSolver % Values, &
         "Remesh Remove Nodes Closer Than", Found, UnfoundFatal=.TRUE.)

    SaveMeshRmThresh = ListGetConstReal(RemeshSolver % Values, &
         "Remesh Remove Nodes Deviation Threshold", Found, UnfoundFatal=.TRUE.)
  END IF

  !Get current simulation time
  TimeVar => VariableGet(Model % Variables, 'Time', .TRUE.)

  MaxNSDiverge = ListGetConstReal(Params, "Maximum Flow Solution Divergence", CheckFlowDiverge)
  MaxNSValue = ListGetConstReal(Params, "Maximum Velocity Magnitude", CheckFlowMax)
  FirstMaxNSValue = ListGetConstReal(Params, "First Time Max Expected Velocity", Found)
  IF(.NOT. Found .AND. CheckFlowDiverge) THEN
    CALL Info(SolverName, "'First Time Max Expected Velocity' not found, setting to 1.0E4")
    FirstMaxNSValue = 1.0E4
  END IF

  !====================================================!
  !---------------- DO THINGS! ------------------------!
  !====================================================!

  NSFail=.FALSE.
  NSDiverge=.FALSE.
  NSTooFast=.FALSE.

  !In addition to checking for absolute failure (% NonlinConverged = 0), we can check
  !for suspiciously large shift in the max variable value (this usually indicates a problem)
  !and we can also check for unphysically large velocity values
  IF(CheckFlowDiverge .OR. CheckFlowMax) THEN

    FlowMax = 0.0_dp
    DO i=1,Mesh % NumberOfNodes
      Mag = 0.0_dp

      DO j=1,FlowVar % DOFs-1
        Mag = Mag + (FlowVar % Values( (FlowVar % Perm(i)-1)*FlowVar % DOFs + j ) ** 2.0_dp)
      END DO
      Mag = Mag ** 0.5_dp
      FlowMax = MAX(FlowMax, Mag)
    END DO

    CALL MPI_AllReduce(MPI_IN_PLACE, FlowMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)
  END IF

  IF(CheckFlowDiverge) THEN
    !First time, there's no previous velocity against which to check divergence.
    !This is somewhat messy because of the separate 'Maximum Velocity Magnitude'
    IF(FirstTime) SaveFlowMax = MIN(FlowMax,FirstMaxNSValue)

    NSChange = FlowMax / SaveFlowMax
    PRINT *,'Debug, Flow Max (old/new): ',SaveFlowMax, FlowMax,' NSChange: ',NSChange

    IF(NSChange > MaxNSDiverge) THEN
      NSDiverge = .TRUE.
      CALL Info(SolverName,"Large change in maximum velocity suggests dodgy&
           &Navier-Stokes solution.")
    END IF
    IF(.NOT. NSDiverge) SaveFlowMax = FlowMax
  END IF

  IF(CheckFlowMax) THEN
    IF(FlowMax > MaxNSValue) THEN
      NSTooFast = .TRUE.
      CALL Info(SolverName,"Large maximum velocity suggests dodgy&
           &Navier-Stokes solution.")
    END IF
  END IF

  NSFail = FlowVar % NonlinConverged < 1 .OR. NSDiverge .OR. NSTooFast

  IF(NSFail) THEN
    CALL Info(SolverName, "Skipping solvers except Remesh because NS failed to converge.")

    FailCount = FailCount + 1
    IF(FailCount >= 4) THEN
       CALL Fatal(SolverName, "Don't seem to be able to recover from NS failure, giving up...")
    END IF

    !Set the clock back one second less than a timestep size.
    !This means next time we are effectively redoing the same timestep
    !but without any solvers which are dependent on (t > told) to reallocate
    TimeVar % Values(1) = TimeVar % Values(1) - SaveDt + (1.0/(365.25 * 24 * 60 * 60.0_dp))


    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", 0.0_dp)
    CALL ListAddInteger( FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", 10000)

    !If this is the second failure in a row, fiddle with the mesh
    IF(FailCount >= 2) THEN
      CALL Info(SolverName,"NS failed twice, fiddling with the mesh...")
      CALL ListAddConstReal(RemeshSolver % Values, &
           "Remesh Min Characteristic Length", 150.0_dp)
      CALL ListAddConstReal(RemeshSolver % Values, &
           "Remesh Remove Nodes Closer Than", 120.0_dp)
      CALL ListAddConstReal(RemeshSolver % Values, &
           "Remesh Remove Nodes Deviation Threshold", 50.0_dp)
    END IF

    IF( .NOT. (NSTooFast .OR. NSDiverge)) THEN
      !---Not quite converging---!

      CALL ListAddConstReal(FlowVar % Solver % Values, &
           "Nonlinear System Relaxation Factor", 0.9_dp)

    ELSE
      !---Solution blowing up----!

      !Set var values to zero so it doesn't mess up viscosity next time
      FlowVar % Values = 0.0_dp

      !TODO: What else? different linear method? more relaxation?
    END IF
  ELSE

    FailCount = 0
    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", SaveNewtonTol)
    CALL ListAddInteger( FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", SaveNewtonIter)
    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Relaxation Factor", SaveRelax)

    CALL ListAddConstReal(RemeshSolver % Values, &
         "Remesh Min Characteristic Length", SaveMeshMinLC)
    CALL ListAddConstReal(RemeshSolver % Values, &
         "Remesh Remove Nodes Closer Than", SaveMeshRmLC)
    CALL ListAddConstReal(RemeshSolver % Values, &
         "Remesh Remove Nodes Deviation Threshold", SaveMeshRmThresh)

  END IF

  !Set a simulation switch to be picked up by Remesher
  CALL ListAddLogical( Model % Simulation, 'Flow Solution Failed', NSFail )

  !Switch off solvers
  DO Num = 1,999
    WRITE(Message,'(A,I0)') 'Switch Off Equation ',Num
    EqName = ListGetString( Params, Message, Found)
    IF( .NOT. Found) EXIT

    Found = .FALSE.
    DO j=1,Model % NumberOfSolvers
      IF(ListGetString(Model % Solvers(j) % Values, "Equation") == EqName) THEN
        Found = .TRUE.
        !Turn off (or on) the solver
        !If NS failed to converge, (switch) off = .true.
        CALL SwitchSolverExec(Model % Solvers(j), NSFail)
        EXIT
      END IF
    END DO

    IF(.NOT. Found) THEN
      WRITE (Message,'(A,A,A)') "Failed to find Equation Name: ",EqName,&
           " to switch off after calving."
      CALL Fatal(SolverName,Message)
    END IF
  END DO

  FirstTime = .FALSE.

END SUBROUTINE CheckFlowConvergence

SUBROUTINE Remesher( Model, Solver, dt, Transient )

  USE DefUtils
  USE GeneralUtils
  USE ElementDescription
  USE MeshUtils  
  USE SParIterComm
  USE CalvingGeometry

  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Transient
  !-------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  TYPE(Variable_t), POINTER :: Var, RefVar, TimeVar, CalvingVar, TangledVar
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) ::FrontOrientation(3), RotationMatrix(3,3), UnRotationMatrix(3,3), NodeHolder(3)
  REAL(KIND=dp), POINTER :: PArray(:,:) => NULL(), TimestepSizes(:,:)
  REAL(KIND=dp) :: time, dt, PseudoSSdt, SaveDt, LastRemeshTime, TimeSinceRemesh, ForceRemeshTime,&
       ZThresh, global_eps, local_eps
  LOGICAL :: Debug, Parallel, CalvingOccurs, RemeshOccurs, PauseSolvers, Found, &
       RotFS, FirstTime = .TRUE.,CalvingLastTime, PauseAfterCalving, FrontalBecomeBasal, &
       TangleOccurs, CheckTangled, NSFail, CheckFlowConvergence
  LOGICAL, POINTER :: NewBasalNode(:)=>NULL()
  CHARACTER(MAX_NAME_LEN) :: SolverName, VarName, EqName, CalvingVarName,&
       FrontMaskName,InMaskName,TopMaskName,BotMaskName,LeftMaskName,RightMaskName, &
       TangledVarName
  INTEGER :: Num, i,j, SaveSSiter, PauseTimeCount=0, PauseTimeMax

  SAVE :: FirstTime, SaveDt, SaveSSiter, PseudoSSdt, LastRemeshTime, ForceRemeshTime, &
       CalvingLastTime,FrontMaskName,InMaskName,TopMaskName,BotMaskName,LeftMaskName,&
       RightMaskName, ZThresh, CalvingVarName, PauseAfterCalving, global_eps, local_eps,&
       PauseTimeCount

  Debug = .FALSE.

  Mesh => Solver % Mesh
  SolverName = 'Remesher'
  Params => Solver % Values
  Parallel = (ParEnv % PEs > 1)

  !-----------------------------------------------------------
  ! Some notes on linking this solver to Calving3D
  !
  ! Input: -CalvingOccurs logical
  !        -Pseudo mesh update for each 
  !         calving front point.
  !
  !-----------------------------------------------------------

  !Trying out new functionality:
  !  CALL Assert(1==0, SolverName, "1 does not equal zero!")
  !  CALL NumericalError( SolverName, "Couldn't do pretend convergence")

  IF(FirstTime) THEN

     TopMaskName = "Top Surface Mask"
     BotMaskName = "Bottom Surface Mask"
     LeftMaskName = "Left Sidewall Mask"
     RightMaskName = "Right Sidewall Mask"
     FrontMaskName = "Calving Front Mask"
     InMaskName = "Inflow Mask"

     LastRemeshTime = GetTime()

     !Need this for temporarily stopping simulation clock when calving occurs,
     ! to recheck for multiple calving events triggered in the same timestep
     TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation, &
          'Timestep Sizes', Found, UnfoundFatal=.TRUE.)
     IF(SIZE(TimestepSizes,1) > 1) CALL Fatal(SolverName,&
          "Calving solver requires a single constant 'Timestep Sizes'")

     SaveDt = TimestepSizes(1,1)

     SaveSSiter = ListGetInteger(Model % Simulation, "Steady State Max Iterations", Found)
     IF(.NOT. Found) SaveSSiter = 1

     PseudoSSdt = ListGetConstReal( Params, 'Pseudo SS dt', Found)
     IF(.NOT. Found) THEN
        CALL Warn(SolverName,"No value specified for 'Pseudo SS dt', taking 1.0e-10")
        PseudoSSdt = 1.0e-10
     END IF

     ForceRemeshTime = ListGetConstReal(Params, 'Force Remesh After Time', Found)
     IF(.NOT. Found) THEN
        CALL Warn(SolverName, 'No "Force Remesh After Time" found, defaulting to 1 month...')
        ForceRemeshTime = 1.0/12.0
     END IF

     !Get the calving variable
     CalvingVarName = ListGetString(Params,"Calving Variable Name", Found)
     IF(.NOT. Found) THEN
        CALL Info(SolverName, "Can't find Calving Variable Name in solver section, &
             & assuming 'Calving'")
        CalvingVarName = "Calving"
     END IF

     !The threshold for converting overhanging frontal elements to basal elements.
     ZThresh = ListGetConstReal(Params, "Front Normal Z Threshold", Found)
     IF(.NOT. Found) THEN
        CALL Warn(SolverName, "Couldn't find Front Normal Z Threshold, setting to -0.9.")
        ZThresh = -0.9
     ELSE
        IF(ZThresh > 0) CALL Fatal(SolverName, "Front Normal Z Threshold controls the &
             &conversion of overhanging frontal elements to basal elements. Sensible values &
             &lie in the range -0.5 to -0.95")
     END IF

     PauseAfterCalving = ListGetLogical(Params, "Pause After Calving Event", Found)
     IF(.NOT. Found) THEN
        CALL Info(SolverName, "Can't find 'Pause After Calving Event' logical in Solver section, &
             & assuming True")
        PauseAfterCalving = .TRUE.
     END IF

     global_eps = 1.0E-2_dp
     local_eps = 1.0E-2_dp
  END IF !FirstTime

  !Get the orientation of the front and compute rotation matrices
  PArray => ListGetConstRealArray( Model % Constants,'Front Orientation', &
       Found, UnfoundFatal=.TRUE.)
  DO i=1,3
     FrontOrientation(i) = PArray(i,1)
  END DO
  RotationMatrix = ComputeRotationMatrix(FrontOrientation)
  UnRotationMatrix = TRANSPOSE(RotationMatrix)

  CALL Info( SolverName, ' ---- Front Rotation Matrix ---- ')
  DO i=1,3
     WRITE(Message, '(f10.7,2x,f10.7,2x,f10.7)') &
          RotationMatrix(i,1),&
          RotationMatrix(i,2),&
          RotationMatrix(i,3)
     CALL Info(SolverName, Message)
  END DO

  !Get current simulation time
  TimeVar => VariableGet(Model % Variables, 'Time', .TRUE.)
  time = TimeVar % Values(1)

  !Is there a calving event?
  CalvingOccurs = ListGetLogical(Model % Simulation, 'CalvingOccurs', Found)
  IF(.NOT.Found) CALL Warn(SolverName, "Unable to find CalvingOccurs logical, &
       & assuming no calving event.")

  !Note - two switches: PauseAfterCalving is the rule in the sif (i.e. do or don't)
  !PauseSolvers is marked by Calving3D if calving event exceeding certain size occurred
  PauseSolvers = ListGetLogical(Model % Simulation, 'Calving Pause Solvers', Found)
  IF(.NOT.Found) THEN
    CALL Warn(SolverName, "Unable to find 'Calving Pause Solvers' logical, &
         & assuming true.")
    PauseSolvers = CalvingOccurs
  END IF

  PauseTimeMax = ListGetInteger(Params, "Calving Pause Max Steps", Found)
  IF(.NOT. Found) THEN
    PauseTimeMax = 15
  END IF

  IF(PauseSolvers) THEN
    PauseTimeCount = PauseTimeCount + 1
    IF(PauseTimeCount > PauseTimeMax) THEN
      PauseSolvers = .FALSE.
      PauseTimeCount = 0
      CALL Info(SolverName,"Calving paused steps exceeded given threshold, moving on...")
    END IF
  ELSE
    PauseTimeCount = 0
  END IF

  CalvingVar => VariableGet(Mesh % Variables, CalvingVarName, .TRUE.)
  IF(.NOT. ASSOCIATED(CalvingVar)) &
       CALL Fatal(SolverName, "Couldn't get Calving variable.")

  !-----------------------------------------
  ! Action!
  !-----------------------------------------

  !Initialize
  RemeshOccurs = .FALSE.

  !If FlowSolver failed to converge (usually a result of weird mesh), large unphysical
  !calving events can be predicted. So, turn off CalvingOccurs, and ensure a remesh
  !Also undo this iterations mesh update
  NSFail = ListGetLogical(Model % Simulation, "Flow Solution Failed",CheckFlowConvergence)
  IF(CheckFlowConvergence) THEN
    IF(NSFail) THEN
      CalvingOccurs = .FALSE.
      RemeshOccurs = .TRUE.
      CALL Info(SolverName, "Remeshing but not calving because NS failed to converge.")
    ELSE

    END IF
  END IF

  !The lagrangian front advance method can sometimes (though rarely) result
  !in columns getting tangled up. FrontAdvance fixes this tangling by restoring
  !the tangled region to a very thin pinnacle (or rift) and then marks them via TangledVar > 0.5
  !When Remeshing, these columns should be removed from the new mesh, and we also ensure that new
  !nodes in FootprintMesh are not in the region where they may be interpolated from the tangled columns
  TangledVarName = ListGetString(Params, "Tangled Variable Name", CheckTangled)
  IF(.NOT. CheckTangled) THEN
    CALL Info(SolverName, "No 'Tangled Variable Name' found, not checking for tangled nodes")
    TangleOccurs = .FALSE.
  ELSE
    CALL Info(SolverName, "Checking for tangled nodes")

    TangledVar => VariableGet(Mesh % Variables, TangledVarName, .TRUE., UnfoundFatal=.TRUE.)

    TangleOccurs = ANY(TangledVar % Values > 0.5)
    IF(Parallel) CALL SParIterAllReduceOR(TangleOccurs)
    IF(TangleOccurs) CALL Info(SolverName, "Some front columns are tangled, remeshing...")
  END IF

  !-------------------------------------------------------------------------

  IF(CalvingOccurs) THEN
    !--Move front nodes--
    CALL DisplaceCalvingFront(Mesh, CalvingVar, 1)
  ELSE
    !-- or ensure set to zero
    CalvingVar % Values = 0.0_dp
  END IF
  !----------------------------------------------------
  ! Check front elements normal vectors
  ! convert downward pointing into basal elements, and nodes
  !----------------------------------------------------
  CALL ConvertFrontalToBasal(Model, Mesh, FrontMaskName, BotMaskName, ZThresh, &
       NewBasalNode, FrontalBecomeBasal)

  IF(CalvingOccurs) LastRemeshTime = time
  TimeSinceRemesh = time - LastRemeshTime

  !Force remesh every so often to maintain mesh quality
  IF( (TimeSinceRemesh > ForceRemeshTime) .OR. FrontalBecomeBasal) THEN
     RemeshOccurs = .TRUE.
     LastRemeshTime = time
  END IF

  IF(CalvingOccurs .OR. RemeshOccurs .OR. TangleOccurs) THEN
    !TODO: there are evidently some large buffered sends in this subroutine
    !Find them and calculate this properly...
    CALL CheckBuffer(104857600)

     CALL Info( SolverName, ' ',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     IF(CalvingOccurs) THEN
        WRITE(Message,'(A,f8.4)') " Calving Event at time: ",time
     ELSE IF(TangleOccurs) THEN
        WRITE(Message,'(A,f8.4)') " Tangled columns, forcing glacier remesh at time: ",time
     ELSE
        WRITE(Message,'(A,f8.4)') " Forcing glacier remesh at time: ",time
     END IF
     CALL Info( SolverName, Message,Level=4 )
     CALL Info( SolverName, ' ',Level=4 )
     CALL Info( SolverName, ' Remeshing Glacier',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ',Level=4 )

     CALL CalvingRemesh(Model, Solver, Mesh, NewMesh, Parallel, Transient)

     NewMesh % Name = TRIM(Mesh % Name)
     NewMesh % OutputActive = .TRUE.
     NewMesh % Changed = .TRUE. 

     CALL SwitchMesh(Model, Solver, Mesh, NewMesh)
     CALL MeshStabParams( Model % Mesh )

     CALL Info( SolverName, ' ',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' Remeshing Complete',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ',Level=4 )

  ELSE

     CALL Info( SolverName, ' ',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' No calving or remesh, doing nothing...',Level=4 )
     CALL Info( SolverName, '-------------------------------------',Level=4 )
     CALL Info( SolverName, ' ',Level=4 )

     !Mesh % Changed forces solvers to reallocate their internal arrays
     !Needs to be .TRUE. for first timestep after a calving event, as 
     !Mesh Update and FreeSurface solvers haven't been called in the mean time
     IF(.NOT. CalvingLastTime) Mesh % Changed = .FALSE.
  END IF

  !Reset listed mesh update variable values to zero
  !Regardless of whether calving occurs.
  DO Num = 1,999
     WRITE(Message,'(A,I0)') 'Mesh Update Variable ',Num
     VarName = ListGetString( Params, Message, Found)
     IF( .NOT. Found) EXIT

     Var => VariableGet( Model % Mesh % Variables, VarName, .TRUE. )
     IF(.NOT. ASSOCIATED(Var)) THEN
        WRITE(Message,'(A,A)') "Listed mesh update variable but cant find: ",VarName
        CALL Fatal(SolverName, Message)
     END IF
     Var % Values = 0.0_dp

     !Turn off (or on) the solver
     !If CalvingOccurs, (switch) off = .true.
     IF(PauseAfterCalving) CALL SwitchSolverExec(Var % Solver, (CalvingOccurs .AND. PauseSolvers))
  END DO

  !Turn off free surface solvers for next timestep
  !And set values equal to z (or rotated) coordinate
  DO Num = 1,999
     WRITE(Message,'(A,I0)') 'FreeSurface Variable ',Num
     VarName = ListGetString( Params, Message, Found)
     IF( .NOT. Found) EXIT

     Var => VariableGet( Model % Mesh % Variables, VarName, .TRUE. )
     IF(.NOT. ASSOCIATED(Var)) THEN
        WRITE(Message,'(A,A)') "Listed FreeSurface variable but cant find: ",VarName
        CALL Fatal(SolverName, Message)
     END IF

     RefVar => VariableGet( Model % Mesh % Variables, "Reference "//TRIM(VarName), .TRUE. )
     IF(.NOT. ASSOCIATED(RefVar)) THEN
        WRITE(Message,'(A,A)') "Listed FreeSurface variable but cant find: ",&
             "Reference "//TRIM(VarName)
        CALL Fatal(SolverName, Message)
     END IF

     WRITE(Message, '(A,A)') TRIM(Message) // " Rotated"
     RotFS = ListGetLogical(Params, Message, Found)
     IF(.NOT. Found) RotFS = .FALSE.

     IF(RotFS) THEN !calving "zs" front
        DO i=1,Model % Mesh % NumberOfNodes
           IF(Var % Perm(i) <= 0) CYCLE
           NodeHolder(1) = Model % Mesh % Nodes % x(i)
           NodeHolder(2) = Model % Mesh % Nodes % y(i)
           NodeHolder(3) = Model % Mesh % Nodes % z(i)
           NodeHolder = MATMUL(RotationMatrix, NodeHolder)
           Var % Values(Var % Perm(i)) = NodeHolder(3)
           RefVar % Values(RefVar % Perm(i)) = NodeHolder(3)
        END DO
     ELSE
        DO i=1,Model % Mesh % NumberOfNodes
           IF(Var % Perm(i) <= 0) CYCLE
           Var % Values(Var % Perm(i)) = Model % Mesh % Nodes % z(i)
           RefVar % Values(RefVar % Perm(i)) = Model % Mesh % Nodes % z(i)
        END DO
     END IF
     
     !Turn off (or on) the solver
     !If CalvingOccurs, (switch) off = .true.
     IF(PauseAfterCalving) CALL SwitchSolverExec(Var % Solver, (CalvingOccurs .AND. PauseSolvers))
  END DO

  IF(PauseAfterCalving) THEN

     IF(CalvingOccurs .AND. PauseSolvers) THEN
        CALL ListAddConstReal( Model % Simulation, 'Timestep Size', PseudoSSdt)
        CALL ListAddInteger( Model % Simulation, 'Steady State Max Iterations', 1)
     ELSE
        CALL ListAddConstReal( Model % Simulation, 'Timestep Size', SaveDt)
        CALL ListAddInteger( Model % Simulation, 'Steady State Max Iterations', SaveSSiter)
     END IF

     DO Num = 1,999
        WRITE(Message,'(A,I0)') 'Switch Off Equation ',Num
        EqName = ListGetString( Params, Message, Found)
        IF( .NOT. Found) EXIT

        Found = .FALSE.
        DO j=1,Model % NumberOfSolvers
           IF(ListGetString(Model % Solvers(j) % Values, "Equation") == EqName) THEN
              Found = .TRUE.
              !Turn off (or on) the solver
              !If CalvingOccurs, (switch) off = .true.
              CALL SwitchSolverExec(Model % Solvers(j), (CalvingOccurs .AND. PauseSolvers))
              EXIT
           END IF
        END DO

        IF(.NOT. Found) THEN
           WRITE (Message,'(A,A,A)') "Failed to find Equation Name: ",EqName,&
                " to switch off after calving."
           CALL Fatal(SolverName,Message)
        END IF
     END DO
  END IF

  !Zero TangledVar on new mesh - this avoids the problem of interpolated values not being
  !updated next timestep, leading to erroneous tangling.
  IF(CheckTangled) THEN
    TangledVar => VariableGet(Model % Mesh % Variables, TangledVarName, .TRUE., UnfoundFatal=.TRUE.)
    TangledVar % Values = 0.0_dp
  END IF

  FirstTime = .FALSE.
  CalvingLastTime = CalvingOccurs

  IF(ASSOCIATED(NewBasalNode)) DEALLOCATE(NewBasalNode)

CONTAINS

  !--------------------------------------------------------------
  ! Takes an existing extruded mesh with a non-vertical face, 
  ! remeshes the footprint, extrudes the mesh, mesh-updates the
  ! calving front to match the old one, and returns a handle
  ! to the new mesh : NewMesh
  !--------------------------------------------------------------
  SUBROUTINE CalvingRemesh(Model, Solver, Mesh, NewMesh, Parallel, Transient)

    USE CalvingGeometry
    USE MainUtils
    USE InterpVarToVar

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh, FootprintMesh, ExtrudedMesh
    LOGICAL :: Parallel, Transient
    !-------------------------------------------
    TYPE(Mesh_t), POINTER :: OldMesh
    TYPE(Solver_t), POINTER :: MUSolver=>NULL()
    TYPE(Matrix_t), POINTER :: StiffMatrix
    TYPE(Element_t), POINTER :: Element, CurrentElement
    TYPE(ValueList_t), POINTER :: Params, Material
    TYPE(Variable_t), POINTER :: TopVar=>NULL(), BottomVar=>NULL(), OldGLVar, NewGLVar, &
         WorkVar, HeightVar, Var, TimestepVar
    TYPE(Nodes_t), TARGET :: FaceNodesT, LeftNodes, RightNodes, BackNodes, FrontNodes, WorkNodes, Nodes
    TYPE(Nodes_t), POINTER :: WriteNodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: Boss, Found, Debug, FirstTime = .TRUE., BadMesh, &
         TriedMetis(5), ThisBC, DoGL, AllFail, AllFailCatch, RemovedOne, InGroup, First, &
         FixDegenerate, MoveMesh=.FALSE.,CalvingColumn, AnyDegenerate,does_intersect
    LOGICAL, ALLOCATABLE :: RemoveNode(:), IsCalvingNode(:), Degenerate(:)
    LOGICAL, POINTER :: UnfoundNodes(:)=>NULL(), OldElemMask(:)
    REAL(KIND=dp) :: MeshEdgeLC, MeshMinDist, MeshMaxDist, MeshMinLC, MeshMaxLC, MeshRmLC,&
         MeshRmThresh, extrude_localeps, extrude_globaleps, Norm, MuStretchZ, NodeHolder(3), &
         ColMin(3), ColMax(3), p1(2), p2(2), q1(2), q2(2), intersect(2), &
         BotDisplacement, TopDisplacement, Displacement, prop, x,dx,maxdz,maxdzdx,DzDxThresh,&
         DzDxMaxDev,ThisDzDxMaxDev,dist,detJ,ShiftBuffer,&
#ifdef USE_ISO_C_BINDINGS
        rt0, rt
#else
        rt0, rt, RealTime
#endif

    REAL(KIND=dp), POINTER :: TopVarValues(:), BottomVarValues(:), ZeroOneHeight(:),&
         ActualHeight(:), WorkReal(:), WorkReal2(:), ForceVector(:),Basis(:)
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), BedHeight(:), RemovalDeviation(:),&
         RemovalDistance(:), ColumnRange(:),TangledZone(:,:),LocalCalvingVar(:),&
         FrontCalving1Values(:), MyDegenerateCoords(:,:), DegenerateCoords(:,:)
    CHARACTER(MAX_NAME_LEN) :: SolverName, Str, NonVertBCName, MeshDir, MeshName, filename, &
         filename_root, MaskName, MuVarName, TopVarName, BottomVarName,&
         NameSuffix, FrontLineMethod, GLVarName, VarName, MoveMeshDir,MoveMeshFullPath,&
         WorkName
    INTEGER :: i,j,k,n, counter, NoNodes, dummyint, FaceNodeCount, &
         ExtrudedLevels, ExtrudeLevels, NodesPerLevel, start, fin, stride, next, WriteNodeCount, &
         MeshBC,col, dim, MetisMethod, MetisMethodDefault, active, NextBasalPerm, &
         FrontBCtag, GroupCount, GroupEnd, GroupStart, TangledGroups

    INTEGER :: comm, ierr, Me, PEs, TotalNodes, DegenCount
    INTEGER, PARAMETER :: GeoUnit = 10
    INTEGER, ALLOCATABLE :: MyFaceNodeNums(:), PFaceNodeCount(:), FNColumns(:), disps(:), &
         WritePoints(:), LocalTangledNode(:), TangledNode(:), WorkInt(:), TangledColumn(:),&
         PartCountDegenerate(:)
    INTEGER, POINTER :: TopPerm(:)=>NULL(), BotPerm(:)=>NULL(), FrontPerm(:)=>NULL(), &
         ExtrudedFrontPerm(:)=>NULL(),  WorkPerm(:),OrderPerm(:),BList(:), &
         InterpDim(:)=>NULL(), ColumnPerm(:), TopVarPerm(:), FootprintFrontPerm(:)=>NULL(),&
         BottomVarPerm(:), TopVarOldPerm(:), BottomVarOldPerm(:)
    INTEGER, POINTER :: LeftLineNodeNums(:), RightLineNodeNums(:), &
         BackLineNodeNums(:), FrontLineNodeNums(:), NodeNums(:),FaceNodeNums(:)=>NULL()

    SAVE FirstTime, MoveMesh

    INTERFACE
       SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
         !------------------------------------------------------------------------------
         USE Lists
         USE SParIterComm
         USE Interpolation
         USE CoordinateSystems
         !-------------------------------------------------------------------------------
         TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
         TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
         LOGICAL, OPTIONAL :: UseQuadrantTree
         LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
         TYPE(Projector_t), POINTER, OPTIONAL :: Projector
         CHARACTER(LEN=*),OPTIONAL :: MaskName
       END SUBROUTINE InterpolateMeshToMesh
    END INTERFACE

    Debug = .FALSE.

    CALL Info( 'Remesher', '-----------------------------------------------',Level=4 )
    CALL Info( 'Remesher', ' Using calving glacier remeshing implementation',Level=4 )
    CALL Info( 'Remesher', '-----------------------------------------------',Level=4 )

    rt0 = RealTime()

    dim = CoordinateSystemDimension()
    SolverName = "Calving Remesh"
    Params => Solver % Values
    Boss = (ParEnv % MyPE == 0) .OR. (.NOT. Parallel)
    TriedMetis = .FALSE.

    extrude_globaleps = global_eps
    extrude_localeps = local_eps

    OldMesh => Mesh
    DegenCount = 0

    NonVertBCName = ListGetString(Params, "Non-Vertical Face Name", Found, UnfoundFatal=.TRUE.)
    NameSuffix = ListGetString(Params, "Remesh Append Name", Found, UnfoundFatal=.TRUE.)

    !Produce mesh and gmsh .geo filenames
    WRITE(filename_root,'(A,A)') "Remesh_temp",TRIM(NameSuffix)
    WRITE(filename,'(A,A)') TRIM(filename_root), ".geo"

    MoveMeshDir = ListGetString(Params, "Remesh Move Mesh Dir", Found)
    IF(Found) THEN
      MoveMesh = .TRUE.
      CALL Info(SolverName, "Moving temporary mesh files after done")
    END IF

    MetisMethodDefault = ListGetInteger(Params, "Metis Algorithm", Found)
    IF(.NOT. Found ) MetisMethodDefault = 4
    MetisMethod = MetisMethodDefault

    !NOTE: don't first time these, they can be altered if meshing fails
    !TODO: MeshEdgeLC isn't really necessary because we use the Distance field...
    MeshEdgeLC = ListGetConstReal(Params, "Remesh Default Characteristic Length", Found, UnfoundFatal=.TRUE.)
    MeshMinLC = ListGetConstReal(Params, "Remesh Min Characteristic Length", Found, UnfoundFatal=.TRUE.)
    MeshMaxLC = ListGetConstReal(Params, "Remesh Max Characteristic Length", Found, UnfoundFatal=.TRUE.)
    MeshMinDist = ListGetConstReal(Params, "Remesh Min Distance Threshold", Found, UnfoundFatal=.TRUE.)
    MeshMaxDist = ListGetConstReal(Params, "Remesh Max Distance Threshold", Found, UnfoundFatal=.TRUE.)
    MeshRmLC = ListGetConstReal(Params, "Remesh Remove Nodes Closer Than", Found, UnfoundFatal=.TRUE.)
    MeshRmThresh = ListGetConstReal(Params, "Remesh Remove Nodes Deviation Threshold", Found, UnfoundFatal=.TRUE.)

    DzDxThresh = ListGetConstReal(Params, "Remesh Max Displacement Gradient", FixDegenerate)
    IF(FixDegenerate) THEN
      CALL Info(SolverName, &
           "Attempting to prevent degeneracy by limiting node displacement gradient.")

      DzDxMaxDev = ListGetConstReal(Params, "Remesh Displacement Deviation Limit", Found)
      IF(.NOT. Found) THEN
        DzDxMaxDev = 50.0_dp
        CALL Info(SolverName, &
             "No deviation limit found for displacement gradient limitation, setting to 50m.")
      END IF
    END IF

    !How the footprint mesh's calving front is computed:
    FrontLineMethod = ListGetString(Params, "Vertical Front Computation", Found)
    IF(.NOT. Found) FrontLineMethod = "midrange"

    GLVarName = ListGetString(Params, "Grounding Line Variable Name", Found)
    IF(.NOT. Found) THEN
       CALL Info(SolverName, "No 'Grounding Line Variable Name' found, assuming GroundedMask")
       GLVarName = "GroundedMask"
    END IF

    !Return here if new mesh has degenerate elements
8989 CONTINUE

    OldGLVar => VariableGet(OldMesh % Variables, GLVarName, .TRUE.)
    IF(ASSOCIATED(OldGLVar)) THEN
       DoGL = .TRUE.
    ELSE
       DoGL = .FALSE.
       IF(Found) THEN
          CALL Fatal(SolverName, "Specified 'Grounding Line Variable Name' but variable not found!")
       ELSE
          CALL Info(SolverName, "Didn't find GroundedMask, not accounting for Grounding Line in remeshing.")
       END IF
    END IF

    !------------------------------------------
    ! Notes on new strategy:
    !
    ! - deformed by calving above to calc new footprint and 0-1
    !
    ! - when cycling front columns, determine and store 0-1 z height
    !
    ! - extrude footprint (already implemented)
    !
    ! - Change front to 0-1 domain, Rotate, and InterpVarToVar
    !
    ! - Mesh Update NewMesh based on interped values
    !
    ! - Undeform from 0-1 back to real mesh
    !
    ! - InterpVarToVarReduced on top and bottom (post deformed new mesh)
    !
    ! - Poisson eq / Mesh Update?
    !
    !------------------------------------------

    !-------------------------------------------
    ! Create nodal BC perms to make lookup simpler
    !-------------------------------------------
    NoNodes = OldMesh % NumberOfNodes
    ALLOCATE( TopPerm(NoNodes), BotPerm(NoNodes), FrontPerm(NoNodes))

    !Generate perms to quickly get nodes on each boundary
    CALL MakePermUsingMask( Model, Solver, OldMesh, TopMaskName, &
         .FALSE., TopPerm, dummyint)
    CALL MakePermUsingMask( Model, Solver, OldMesh, BotMaskName, &
         .FALSE., BotPerm, dummyint)
    CALL MakePermUsingMask( Model, Solver, OldMesh, FrontMaskName, &
         .FALSE., FrontPerm, FaceNodeCount) !<- Number of nodes in this part on calving face

    IF(FrontalBecomeBasal .AND. Debug) &
         PRINT *,ParEnv % MyPe,'Frontal become basal!'

    
    !---------------------------------------------------
    !            FOOTPRINT GENERATION
    !
    ! Get the global node
    ! numbers, connectivity etc of the footprint. Needs
    ! special care in parallel
    !
    !---------------------------------------------------

    !Get each of the 4 edges of the top surface
    !The result of these calls is only valid in Boss part
    CALL GetDomainEdge(Model, OldMesh, TopPerm, LeftMaskName, &
         LeftNodes, LeftLineNodeNums, Parallel, Simplify=.TRUE.)
    IF(Debug) CALL Info(SolverName, "Done left domain edge")
    CALL GetDomainEdge(Model, OldMesh, TopPerm, RightMaskName, &
         RightNodes, RightLineNodeNums, Parallel, Simplify=.TRUE.)
    IF(Debug) CALL Info(SolverName, "Done right domain edge")
    CALL GetDomainEdge(Model, OldMesh, TopPerm, InMaskName, &
         BackNodes, BackLineNodeNums, Parallel, Simplify=.TRUE.)
    IF(Debug) CALL Info(SolverName, "Done back domain edge")
    CALL GetDomainEdge(Model, OldMesh, TopPerm, FrontMaskName, &
         FrontNodes, FrontLineNodeNums, Parallel, Simplify=.FALSE.)
    IF(Debug) CALL Info(SolverName, "Done front domain edge")

    !Call this again later to remove too close nodes from result
    IF(Boss .AND. Debug) PRINT *, 'Debug CalvingRemesh, FrontLineNodeNums: ',FrontLineNodeNums
    IF(Boss .AND. Debug) PRINT *, 'Debug CalvingRemesh, BackLineNodeNums: ',BackLineNodeNums
    IF(Boss .AND. Debug) PRINT *, 'Debug CalvingRemesh, LeftLineNodeNums: ',LeftLineNodeNums

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken for variable loading, making perms, GetDomainEdge: ', rt
    rt0 = RealTime()

    !---------------------------------------------
    !
    ! Get average calving line
    ! 
    ! TODO: front to basal important here
    !---------------------------------------------

    !Cycle all boundary elements in current partition, get calving face nodes
    ALLOCATE(MyFaceNodeNums(FaceNodeCount))
    j = 0
    DO i=1, NoNodes
       IF(FrontPerm(i) <= 0) CYCLE
       j = j + 1
       MyFaceNodeNums(j) = i
    END DO
    !Now MyFaceNodeNums is a list of all nodes on the calving boundary.

    !Send calving nodes to boss partition
    !This sacrifices a little performance for potentially unneeded generality
    !Because the mesh is extruded, we could simply average each column of points
    !and then send result to boss part.  But let's not, in case we want to do
    !some other kind of smoothing in the future.
    IF(Parallel) THEN

       Me = ParEnv % MyPe
       PEs = ParEnv % PEs
       comm = ParEnv % ActiveComm

       !Send node COUNT
       IF(Boss) ALLOCATE(PFaceNodeCount(PEs))


       CALL MPI_GATHER(FaceNodeCount,1,MPI_INTEGER,PFaceNodeCount,&
            1,MPI_INTEGER, 0, comm, ierr)

       IF(Boss) THEN
          FaceNodesT % NumberOfNodes = SUM(PFaceNodeCount)
          n = FaceNodesT % NumberOfNodes
          ALLOCATE(FaceNodeNums(n),&
               FaceNodesT % x(n),&
               FaceNodesT % y(n),&
               FaceNodesT % z(n),&
               disps(PEs)) !variable to hold array offset from each proc
          !work out where in array to put data from each proc
          !how to deal with zero size?
          disps(1) = 0
          DO i=2,PEs
             disps(i) = disps(i-1) + PFaceNodeCount(i-1)
          END DO
       END IF

       !Global NodeNumbers
       CALL MPI_GATHERV(OldMesh % ParallelInfo % GlobalDOFs(MyFaceNodeNums),&
            FaceNodeCount,MPI_INTEGER,&
            FaceNodeNums,PFaceNodeCount,&
            disps,MPI_INTEGER,0,comm, ierr)
       !X coords
       CALL MPI_GATHERV(OldMesh % Nodes % x(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % x,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       !Y coords
       CALL MPI_GATHERV(OldMesh % Nodes % y(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % y,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       !Z coords
       CALL MPI_GATHERV(OldMesh % Nodes % z(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % z,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)

       !Gather tangled columns info to boss
       IF(TangleOccurs .AND. NSFail) &
            TangleOccurs=.FALSE.

       IF(TangleOccurs) THEN
         IF(Boss) ALLOCATE(TangledNode(n), FrontCalving1Values(n))
         ALLOCATE(LocalTangledNode(FaceNodeCount),&
              LocalCalvingVar(FaceNodeCount))

         LocalTangledNode = NINT(TangledVar % Values(TangledVar % Perm(MyFaceNodeNums)))
         LocalCalvingVar = &
              CalvingVar % Values((CalvingVar % Perm(MyFaceNodeNums)-1)*CalvingVar % DOFs + 1)

         CALL MPI_GATHERV(LocalTangledNode,&
              FaceNodeCount,MPI_INTEGER,&
              TangledNode,PFaceNodeCount,&
              disps,MPI_INTEGER,0,comm, ierr)

         CALL MPI_GATHERV(LocalCalvingVar,&
              FaceNodeCount,MPI_DOUBLE_PRECISION,&
              FrontCalving1Values,PFaceNodeCount,&
              disps,MPI_DOUBLE_PRECISION,0,comm, ierr)

         IF(Boss) THEN
           TangledGroups = MAXVAL(TangledNode)
         END IF

         CALL MPI_BCAST( TangledGroups, 1, MPI_INTEGER, 0, comm, ierr)
         ALLOCATE(TangledZone(TangledGroups,2))
       END IF

       IF(Boss) THEN
          IF(Debug) THEN
             PRINT *, 'Debug Remesh, pre removal nodes: '
             DO i=1,FaceNodesT % NumberOfNodes
                PRINT *, FaceNodesT % x(i),FaceNodesT % y(i),FaceNodesT % z(i)
             END DO
          END IF

          !Remove duplicates!!!
          ALLOCATE(RemoveNode(FaceNodesT % NumberOfNodes))
          RemoveNode = .FALSE.
          DO i=2,FaceNodesT % NumberOfNodes
             IF(ANY(FaceNodeNums(1:i-1) == FaceNodeNums(i))) THEN
                RemoveNode(i) = .TRUE.
             END IF
          END DO

          !Remove duplicate values from TangledNode
          IF(TangleOccurs) THEN
            ALLOCATE(WorkInt(COUNT(.NOT. RemoveNode)))
            WorkInt = PACK(TangledNode, .NOT. RemoveNode)
            DEALLOCATE(TangledNode)
            ALLOCATE(TangledNode(SIZE(WorkInt)))
            TangledNode = WorkInt
            DEALLOCATE(WorkInt)

            ALLOCATE(WorkReal(COUNT(.NOT. RemoveNode)))
            WorkReal = PACK(FrontCalving1Values, .NOT. RemoveNode)
            DEALLOCATE(FrontCalving1Values)
            ALLOCATE(FrontCalving1Values(SIZE(WorkReal)))
            FrontCalving1Values = WorkReal
            DEALLOCATE(WorkReal)
          END IF

          CALL RemoveNodes(FaceNodesT, RemoveNode, FaceNodeNums)
          DEALLOCATE(RemoveNode)

          IF(Debug) THEN
             PRINT *, 'Size of FaceNodeNums: ', SIZE(FaceNodeNums)
             PRINT *, 'Debug Remesh, post removal nodes: '
             DO i=1,FaceNodesT % NumberOfNodes
                PRINT *, FaceNodesT % x(i),FaceNodesT % y(i),FaceNodesT % z(i)
             END DO
          END IF
       END IF

    ELSE !Serial
       n = FaceNodeCount
       ALLOCATE(FaceNodeNums(n),&
            FaceNodesT % x(n),&
            FaceNodesT % y(n),&
            FaceNodesT % z(n))

       !This seems a little redundant but it saves
       !some lines of code later on...
       FaceNodeNums = MyFaceNodeNums
       FaceNodesT % x = OldMesh % Nodes % x(MyFaceNodeNums)
       FaceNodesT % y = OldMesh % Nodes % y(MyFaceNodeNums)
       FaceNodesT % z = OldMesh % Nodes % z(MyFaceNodeNums)
    END IF

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken for collecting front node, removing duplicates: ', rt
    rt0 = RealTime()

    !--------------------------------------------------------------------

    !Need global mesh structure info
    IF(Parallel) THEN
       !Rather than summing NoNodes from each part, we simply find
       !the maximum global node number
       CALL MPI_Reduce(MAXVAL(OldMesh % ParallelInfo % GlobalDOFs), TotalNodes, &
            1, MPI_INTEGER, MPI_MAX, 0,comm,ierr)
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    ELSE
       TotalNodes = NoNodes
    END IF

    ExtrudedLevels = GetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found)
    IF(.NOT. Found) ExtrudedLevels = &
         GetInteger(CurrentModel % Simulation,'Remesh Extruded Mesh Levels',Found)
    IF(.NOT. Found) CALL Fatal("Remesh",&
         "Unable to find 'Extruded Mesh Levels' or 'Remesh Extruded Mesh Levels'")

    IF(Boss) THEN !Master/Slave problem (or serial)

       IF(Debug) PRINT *, 'Debug remesh, Total nodes: ', TotalNodes

       IF(MOD(TotalNodes, ExtrudedLevels) /= 0) &
            CALL Fatal("Remesh","Total number of nodes isn't divisible by number&
            &of mesh levels. WHY?")

       NodesPerLevel = TotalNodes / ExtrudedLevels

       PRINT *, 'Debug CalvingRemesh, NodesPerLevel: ', NodesPerLevel

       !note, if ExtrudedLevels is messed with in remeshing, it'll need to be
       !pushed back to simulation

       !Now columns *should* have the same integer FNColumns
       !NOTE: Think about how messing with BCs following undercutting affects this
       ALLOCATE(FNColumns(FaceNodesT % NumberOfNodes), &
            ColumnRange(SIZE(FrontLineNodeNums)))

       FNColumns = MOD(FaceNodeNums, NodesPerLevel)

       !Check for tangled columns
       IF(TangleOccurs) THEN
         ALLOCATE(TangledColumn(SIZE(FrontLineNodeNums)))
         TangledColumn = 0
         IF(Debug) PRINT *,'tangled nodes: ',COUNT(TangledNode > 0)

         DO i=1,TangledGroups
           TangledZone(i,1) = HUGE(TangledZone(i,1))
           TangledZone(i,2) = -HUGE(TangledZone(i,2))
         END DO

         DO i=1,SIZE(FrontLineNodeNums) !cycle columns
           j = MOD(FrontLineNodeNums(i), NodesPerLevel)
           CalvingColumn = .FALSE.
           !If any node in column is tangled, mark column
           !   (in fact all nodes in column will be)
           DO k=1,SIZE(FNColumns)
             IF(FNColumns(k)==j) THEN
               IF(TangledNode(k) > 0) THEN
                 TangledColumn(i) = TangledNode(k)
                 IF(FrontCalving1Values(k) /= 0.0_dp) CalvingColumn = .TRUE.
               END IF
             END IF
           END DO

           !Mark the extent of tangling
           IF((TangledColumn(i) > 0) .AND. .NOT. CalvingColumn) THEN
             NodeHolder(1) = FrontNodes % x(i)
             NodeHolder(2) = FrontNodes % y(i)
             NodeHolder(3) = FrontNodes % z(i)
             NodeHolder = MATMUL(RotationMatrix, NodeHolder)
             TangledZone(TangledColumn(i),1) = MIN(TangledZone(TangledColumn(i),1), NodeHolder(2))
             TangledZone(TangledColumn(i),2) = MAX(TangledZone(TangledColumn(i),2), NodeHolder(2))
           END IF
         END DO
         PRINT *,'Remesh, identified ',COUNT(TangledColumn > 0), ' tangled columns.'
         IF(COUNT(TangledColumn > 0) == 0) CALL Fatal(SolverName, &
              "TangleOccurs is true, but found no tangled columns...")
       END IF

       !---------------------------------------------------------------------
       ! Cycle FrontLineNodeNums, finding column averages and updating locations
       !---------------------------------------------------------------------
       FrontNodes % x = 0
       FrontNodes % y = 0
       FrontNodes % z = 0

       DO i=1,SIZE(FrontLineNodeNums) !cycle columns
          j = MOD(FrontLineNodeNums(i), NodesPerLevel)
          WorkNodes % NumberOfNodes = COUNT(FNColumns == j)
          n = WorkNodes % NumberOfNodes

          IF(n < 2) CALL Fatal("CalvingRemesh",&
               "Found fewer than 2 nodes for a column of calving face nodes.")
          IF(Debug) THEN
             PRINT *, 'Debug CalvingRemesh, number of worknodes: ',n
          END IF

          ALLOCATE(WorkNodes % x(n),&
               WorkNodes % y(n),&
               WorkNodes % z(n))

          !----------------------------------------------
          ! Compute and store column horizontal displacement range
          ColMin(3) = HUGE(ColMin(3))
          ColMax(3) = -HUGE(ColMax(3))
          DO k=1,SIZE(FNColumns)
            IF(FNColumns(k)==j) THEN
              NodeHolder(1) = FaceNodesT % x(k)
              NodeHolder(2) = FaceNodesT % y(k)
              NodeHolder(3) = FaceNodesT % z(k)

              NodeHolder = MATMUL(RotationMatrix, NodeHolder)

              IF(NodeHolder(3) < ColMin(3)) ColMin = NodeHolder

              IF(NodeHolder(3) > ColMax(3)) ColMax = NodeHolder
            END IF
          END DO
          ColumnRange(i) = ColMax(3) - ColMin(3)

          !----------------------------------------------
          ! Compute front points for plane mesh

          SELECT CASE(FrontLineMethod)
          CASE("mean")

             counter=1
             DO k=1,SIZE(FNColumns)
                IF(FNColumns(k)==j) THEN
                   WorkNodes % x(counter) = FaceNodesT % x(k)
                   WorkNodes % y(counter) = FaceNodesT % y(k)
                   WorkNodes % z(counter) = FaceNodesT % z(k)
                   counter = counter + 1
                END IF
             END DO

             !Order by ascending WorkNodes % z
             ALLOCATE(OrderPerm(n))
             DO k=1,n; OrderPerm(k) = k;
             END DO

             CALL SortD( n, WorkNodes % z, OrderPerm )
             CALL MySortF( n, OrderPerm, WorkNodes % x )
             CALL MySortF( n, OrderPerm, WorkNodes % y )

             !Zs aren't necessarily equidistant, better: mean = integral/total_z
             !**Trapezoid rule over nodes in column**
             DO k=2,n
                FrontNodes % x(i) = FrontNodes % x(i) + &
                     ((WorkNodes % x(k) + WorkNodes % x(k-1))/2.0_dp) * &
                     (WorkNodes % z(k) - WorkNodes % z(k-1))

                FrontNodes % y(i) = FrontNodes % y(i) + &
                     ((WorkNodes % y(k) + WorkNodes % y(k-1))/2.0_dp) * &
                     (WorkNodes % z(k) - WorkNodes % z(k-1))
             END DO

             FrontNodes % x(i) = FrontNodes % x(i) / (WorkNodes % z(n) - WorkNodes % z(1))
             FrontNodes % y(i) = FrontNodes % y(i) / (WorkNodes % z(n) - WorkNodes % z(1))

             IF(Debug) PRINT *,'Debug, calving column ',i,' has points :',&
                  FrontNodes % x(i), FrontNodes % y(i)

             DEALLOCATE(OrderPerm)

          CASE("minimum")

             ColMin(3) = HUGE(ColMin(3))
             DO k=1,SIZE(FNColumns)
                IF(FNColumns(k)==j) THEN
                   NodeHolder(1) = FaceNodesT % x(k)
                   NodeHolder(2) = FaceNodesT % y(k)
                   NodeHolder(3) = FaceNodesT % z(k)

                   NodeHolder = MATMUL(RotationMatrix, NodeHolder)

                   IF(NodeHolder(3) < ColMin(3)) ColMin = NodeHolder
                END IF
             END DO

             NodeHolder = MATMUL(UnRotationMatrix, ColMin)

             FrontNodes % x(i) = NodeHolder(1)
             FrontNodes % y(i) = NodeHolder(2)

             IF(Debug) PRINT *,'Debug, calving column ',i,' has points :',&
                  FrontNodes % x(i), FrontNodes % y(i)

          CASE("midrange")

             ColMin(3) = HUGE(ColMin(3))
             ColMax(3) = -HUGE(ColMax(3))
             DO k=1,SIZE(FNColumns)
                IF(FNColumns(k)==j) THEN
                   NodeHolder(1) = FaceNodesT % x(k)
                   NodeHolder(2) = FaceNodesT % y(k)
                   NodeHolder(3) = FaceNodesT % z(k)

                   NodeHolder = MATMUL(RotationMatrix, NodeHolder)

                   IF(NodeHolder(3) < ColMin(3)) ColMin = NodeHolder

                   IF(NodeHolder(3) > ColMax(3)) ColMax = NodeHolder
                END IF
             END DO

             ! NodeHolder(1) = MinColX + 0.5*(MaxColX - MinColX)
             ! NodeHolder(2) = MinColY + 0.5*(MaxColY - MinColY)
             ! NodeHolder(3) = MinColZ + 0.5*(MaxColZ - MinColZ)
             NodeHolder = ColMin + 0.5*(ColMax - ColMin)
             NodeHolder = MATMUL(UnRotationMatrix, NodeHolder)

             FrontNodes % x(i) = NodeHolder(1)
             FrontNodes % y(i) = NodeHolder(2)

             IF(Debug) PRINT *,'Debug, calving column ',i,' has points :',&
                  FrontNodes % x(i), FrontNodes % y(i)
          CASE DEFAULT
             CALL Fatal(SolverName, "Invalid choice given for 'Vertical Front Computation'")
          END SELECT

          DEALLOCATE(WorkNodes % x, WorkNodes % y, WorkNodes % z)
       END DO

       DEALLOCATE(FNColumns)
       IF(Debug) THEN
          PRINT *, 'Debug CalvingRemesh, FrontNodes: '
          DO i=1,FrontNodes % NumberOfNodes
             PRINT *, 'node: ',i,' x: ',FrontNodes % x(i),' y: ', FrontNodes % y(i)
          END DO
       END IF

    END IF !Boss

    !------------------------------------------------------------
    ! Calculate ZeroOneHeight for all calving front nodes
    !------------------------------------------------------------

    !Send NodesPerLevel (to get column mod)
    CALL MPI_BCAST( NodesPerLevel, 1, MPI_INTEGER, 0, comm, ierr)

    ALLOCATE(FNColumns(OldMesh % NumberOfNodes),&
         ZeroOneHeight(FaceNodeCount))

    FNColumns = MOD(OldMesh % ParallelInfo % GlobalDOFs, NodesPerLevel)
    ZeroOneHeight = -1.0_dp

    DO WHILE(.TRUE.) !cycle columns

       !Find a new column
       col = -1
       DO j=1, OldMesh % NumberOfNodes
          IF(FrontPerm(j) <= 0) CYCLE
          IF(ZeroOneHeight(FrontPerm(j)) == -1.0_dp) THEN
             col = FNColumns(j)
             EXIT
          END IF
       END DO
       IF(col == -1) EXIT !All done

       !Gather front nodes in specified column
       WorkNodes % NumberOfNodes = COUNT((FrontPerm > 0) .AND. (FNColumns == col))
       n = WorkNodes % NumberOfNodes
       ALLOCATE(WorkNodes % z(n), ColumnPerm(n))

       counter = 1
       DO j=1, OldMesh % NumberOfNodes
          IF(FrontPerm(j) <= 0) CYCLE
          IF(FNColumns(j) == col) THEN
             WorkNodes % z(counter) = OldMesh % Nodes % z(j)
             ColumnPerm(counter) = j
             counter = counter + 1
          END IF
       END DO

       !Calc ZeroOneHeight
       !Take advantage of the fact that nodes are numbered up from the bottom level, so 
       !worknodes(1) is the lowest, and worknodes(n) is the top.
       DO k=1,n !Node elevation in 0-1
          ZeroOneHeight(FrontPerm(ColumnPerm(k))) = (WorkNodes % z(k) - WorkNodes % z(1)) / &
               (WorkNodes % z(n) - WorkNodes % z(1))
       END DO

       DEALLOCATE(WorkNodes % z, ColumnPerm)
    END DO

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken for calculating zero one height: ', rt
    rt0 = RealTime()

    !Remove very close nodes & those which will lead to mesh degeneracy
    !Also remove those which were marked by FrontAdvance3D.F90 as being tangled
    IF(Boss) THEN
       !------------------------------------------------------------
       ! Check no calving nodes too close (too far will be dealt with by gmsh)
       !
       ! Strategy:
       !   For each node, calculate:
       !          - Mean dist from 2 neighbours
       !          - Distance from point to line which would result from its removal
       !            Effectively, this is the resulting 'deviation' of the front
       !
       !   Each cycle, remove the best candidate:
       !          lowest distance & below deviation threshold
       !
       ! This isn't terribly efficient, everything is recomputed each time a
       ! node is removed, but it's not very computationally expensive...
       !------------------------------------------------------------
       ALLOCATE(RemovalDeviation(FrontNodes % NumberOfNodes), &
            RemovalDistance(FrontNodes % NumberOfNodes),&
            RemoveNode(FrontNodes % NumberOfNodes))

       RemoveNode = .FALSE.

       !Remove tangled columns
       IF(TangleOccurs) THEN
         DO i=2,FrontNodes % NumberOfNodes-1
           IF(TangledColumn(i) > 0) RemoveNode(i) = .TRUE.
         END DO
       END IF

       RemovedOne = .TRUE.
       DO WHILE(RemovedOne)
          RemovedOne = .FALSE.

          RemovalDeviation = HUGE(0.0_dp)
          RemovalDistance  = HUGE(0.0_dp)

          !find deviation of point in front direction if removed...
          !and average distance between each node and its neighbours
          !NOTE: code duplication here from above
          DO i=2,FrontNodes % NumberOfNodes-1

             IF(RemoveNode(i)) CYCLE !already got

             j = i - 1
             k = i + 1

             !If neighbours are marked for removal, look for next available
             !neighbour in each direction
             DO WHILE(RemoveNode(j))
                j = j-1
             END DO
             DO WHILE(RemoveNode(k))
                k = k+1
             END DO

             p1(1) = FrontNodes % x(j) !prev point
             p1(2) = FrontNodes % y(j)
             p2(1) = FrontNodes % x(k) !next point
             p2(2) = FrontNodes % y(k)
             q1(1) = FrontNodes % x(i) !this point
             q1(2) = FrontNodes % y(i)
             q2(1) = FrontNodes % x(i) + FrontOrientation(1)
             q2(2) = FrontNodes % y(i) + FrontOrientation(2)

             !Find where this node would theoretically sit between its
             !neighbours, were it removed. Then find the distance
             !p1,p2 are the surrounding points on front
             !q1 is this point, q2 is another point on the line from
             !this point, normal to the front direction
             CALL LinesIntersect ( p1, p2, q1, q2, intersect, does_intersect )

             IF(does_intersect) THEN
                RemovalDeviation(i) = &
                     ( ((q1(1) - intersect(1))**2)  +  ((q1(2) - intersect(2))**2) ) ** 0.5
             ELSE
                !Rare (impossible?) case where line between the two neighbours
                !is parallel to front orientation.
                RemovalDeviation(i) = 0.0_dp
             END IF

             RemovalDistance(i) = (NodeDist2D(FrontNodes,i,j) + NodeDist2D(FrontNodes,i,k)) / 2.0_dp

             !Scale RemovalDeviation by removal distance.
             !For a given front deviation, the larger the RemovalDistance, the larger
             !the volume artificially missing from the front.
             !So we scale dev = dev * (maxdist / dist)
             !Thus, note that given value for MeshRmThresh is for dist = MeshRmLC
             RemovalDistance(i) = MAX(RemovalDistance(i),0.0001_dp) !avoid /0
             RemovalDeviation(i) = RemovalDeviation(i) / (MeshRmLC / RemovalDistance(i))
          END DO

          IF(Debug) THEN
             PRINT *,'Debug minval, minloc, removaldistance', MINVAL(RemovalDistance), &
                  MINLOC(RemovalDistance,1)
             PRINT *,'Debug minval, minloc, removaldeviation', MINVAL(RemovalDeviation), &
                  MINLOC(RemovalDeviation,1)
          END IF

          !If the closest node on the front has a sufficiently low removal deviation, remove it
          !else, look for the second closest, etc.
          DO WHILE(.TRUE.)

             IF(MINVAL(RemovalDistance) > MeshRmLC) EXIT
             IF(MINVAL(RemovalDeviation) > MeshRmThresh) EXIT

             IF(RemovalDeviation(MINLOC(RemovalDistance,1)) < MeshRmThresh) THEN
                IF(Debug) THEN
                   PRINT *, 'Debug CalvingRemesh, MinDist, removing node ',MINLOC(RemovalDistance)&
                        ,' dist: ', MINVAL(RemovalDistance),' deviation: ',&
                        RemovalDeviation(MINLOC(RemovalDistance,1))
                END IF

                RemoveNode(MINLOC(RemovalDistance)) = .TRUE.
                RemovedOne = .TRUE.
                EXIT
             END IF
             RemovalDeviation(MINLOC(RemovalDistance)) = HUGE(0.0_dp)
             RemovalDistance(MINLOC(RemovalDistance)) = HUGE(0.0_dp)
          END DO
       END DO

       !Look for likely degenerates
       IF(FixDegenerate) THEN

         RemovedOne = .TRUE.

         DO WHILE(RemovedOne)
           RemovedOne = .FALSE.

           !find deviation of point in front direction if removed...
           DO i=2,FrontNodes % NumberOfNodes-1

             IF(RemoveNode(i)) CYCLE !already got

             j = i - 1
             k = i + 1

             !If neighbours are marked for removal, look for next available
             !neighbour in each direction
             DO WHILE(RemoveNode(j))
               j = j-1
             END DO

             DO WHILE(RemoveNode(k))
               k = k+1
             END DO

             p1(1) = FrontNodes % x(j) !prev point
             p1(2) = FrontNodes % y(j)
             p2(1) = FrontNodes % x(k) !next point
             p2(2) = FrontNodes % y(k)
             q1(1) = FrontNodes % x(i) !this point
             q1(2) = FrontNodes % y(i)
             q2(1) = FrontNodes % x(i) + FrontOrientation(1)
             q2(2) = FrontNodes % y(i) + FrontOrientation(2)

             !Find where this node would theoretically sit between its
             !neighbours, were it removed. Then find the distance
             !p1,p2 are the surrounding points on front
             !q1 is this point, q2 is another point on the line from
             !this point, normal to the front direction
             CALL LinesIntersect ( p1, p2, q1, q2, intersect, does_intersect )

             IF(does_intersect) THEN
               RemovalDeviation(i) = &
                    ( ((q1(1) - intersect(1))**2)  +  ((q1(2) - intersect(2))**2) ) ** 0.5
             ELSE
               !Rare (impossible?) case where line between the two neighbours
               !is parallel to front orientation.
               RemovalDeviation(i) = 0.0_dp
             END IF
           END DO

           !Look for columns where there is a large amount of horizontal
           !deviation (i.e. melt undercutting), and remove close nodes
           DO i=1,FrontNodes % NumberOfNodes-1
             IF(RemoveNode(i)) CYCLE

             j = i + 1
             DO WHILE(RemoveNode(j))
               j = j+1
             END DO

             !Find front-horizontal distance between nodes (dx)
             NodeHolder(1) = FrontNodes % x(i)
             NodeHolder(2) = FrontNodes % y(i)
             NodeHolder(3) = 0.0_dp
             NodeHolder = MATMUL(RotationMatrix, NodeHolder)

             x = NodeHolder(2)

             NodeHolder(1) = FrontNodes % x(j)
             NodeHolder(2) = FrontNodes % y(j)
             NodeHolder(3) = 0.0_dp
             NodeHolder = MATMUL(RotationMatrix, NodeHolder)

             dx = ABS(x - NodeHolder(2))

             !Compute maximum possible 'deviation' of plane nodes
             maxdz = MAX(ColumnRange(i), ColumnRange(j))
             maxdzdx = maxdz / dx

             dist = NodeDist2D(FrontNodes,i,j)

             ThisDzDxMaxDev = DzDxMaxDev / (DzDxThresh / maxdzdx)

             IF(Debug) PRINT *,'Debug, Remesh, i, dx, maxdz, maxdzdx: ', i, dx, maxdz, maxdzdx

             !Two conditions: 1) high gradient, 2) columns actually close in xy distance
             IF(maxdzdx > DzDxThresh .AND. dist < maxdz) THEN

               !Check removal deviation here
               !If lower columnrange node has low deviation, remove it
               !Else if higher columnrange node has low deviation, remove it
               !Else nothing
               IF(j == FrontNodes % NumberOfNodes) THEN
                 !Special case at end of line
                 IF(RemovalDeviation(i) < ThisDzDxMaxDev) THEN
                   RemoveNode(i) = .TRUE.
                   RemovedOne = .TRUE.
                   EXIT !go back and recalc
                 END IF
               ELSE IF(i == 1) THEN
                 !Special case at start of line
                 IF(RemovalDeviation(j) < ThisDzDxMaxDev) THEN
                   RemoveNode(j) = .TRUE.
                   RemovedOne = .TRUE.
                   EXIT !go back and recalc
                 END IF
               ELSE
                 !Remove whichever node has smallest displacement range
                 IF((ColumnRange(i) < ColumnRange(j)) ) THEN
                   IF(RemovalDeviation(i) < ThisDzDxMaxDev) THEN
                     RemoveNode(i) = .TRUE.
                     RemovedOne = .TRUE.
                     EXIT !go back and recalc
                   ELSE IF(RemovalDeviation(j) < ThisDzDxMaxDev) THEN
                     RemoveNode(j) = .TRUE.
                     RemovedOne = .TRUE.
                     EXIT !go back and recalc
                   END IF
                 ELSE
                   IF(RemovalDeviation(j) < ThisDzDxMaxDev) THEN
                     RemoveNode(j) = .TRUE.
                     RemovedOne = .TRUE.
                     EXIT !go back and recalc
                   ELSE IF(RemovalDeviation(i) < ThisDzDxMaxDev) THEN
                     RemoveNode(i) = .TRUE.
                     RemovedOne = .TRUE.
                     EXIT !go back and recalc
                   END IF
                 END IF
               END IF
             END IF
           END DO

           IF(RemovedOne) PRINT *,'Remesh, removing a column to prevent degeneracy'
         END DO
       END IF !FixDegenerate

       IF(COUNT(RemoveNode) > 0) CALL RemoveNodes(FrontNodes, RemoveNode, FrontLineNodeNums)
       DEALLOCATE(RemovalDeviation, RemovalDistance, RemoveNode)

    END IF

1280 CONTINUE
    !------------------------------------------------------------
    ! Write and generate footprint mesh
    !------------------------------------------------------------
    IF(Boss) THEN
       AllFail = .FALSE.

       OPEN( UNIT=GeoUnit, FILE=filename, STATUS='UNKNOWN') 

       WriteNodeCount = SIZE(FrontLineNodeNums) + &
            SIZE(LeftLineNodeNums) + &
            SIZE(BackLineNodeNums) + &
            SIZE(RightLineNodeNums) - 4

       !-------------Write Points---------------
       ALLOCATE(WritePoints(WriteNodeCount))
       counter = 1
       DO i=1,4 !cycle boundaries
          SELECT CASE(i)
          CASE(1) !left
             WriteNodes => LeftNodes
             NodeNums => LeftLineNodeNums
          CASE(2) !back
             WriteNodes => BackNodes
             NodeNums => BackLineNodeNums
          CASE(3) !right
             WriteNodes => RightNodes             
             NodeNums => RightLineNodeNums
          CASE(4) !front
             WriteNodes => FrontNodes
             NodeNums => FrontLineNodeNums
          END SELECT

          n = WriteNodes % NumberOfNodes
          IF(n /= SIZE(NodeNums)) CALL Fatal("CalvingRemesh","Size mismatch in perm size")

          !Determine order
          !TODO - MMigrate issue here with providing all but last node from each boundary
          ! however, this isn't the root cause...
          IF(i==1) THEN !left edge, find which end neighbours calving front
             IF(ANY(FrontLineNodeNums == LeftLineNodeNums(1))) THEN
                start=1;fin=n-1;stride=1
                next = NodeNums(n)
             ELSE IF(ANY(FrontLineNodeNums == LeftLineNodeNums(SIZE(LeftLineNodeNums)))) THEN
                start=n;fin=2;stride=-1
                next = NodeNums(1)
             ELSE
                CALL Fatal("CalvingRemesh","Problem joining up a closed loop for footprint mesh.")
             END IF
          ELSE
             IF(NodeNums(1)==next) THEN
                start=1;fin=n-1;stride=1
                next = NodeNums(n)
             ELSE IF(NodeNums(n)==next) THEN
                start=n;fin=2;stride=-1
                next = NodeNums(1)
             ELSE
                PRINT *, 'i, NodeNums(1), (n), next: ',i,NodeNums(1), NodeNums(n), next
                CALL Fatal("CalvingRemesh","Problem joining up a closed loop for footprint mesh.")
             END IF
          END IF

          DO j=start,fin,stride !cycle nodes except last, these are overlap
             WRITE( GeoUnit,'(A,I0,A,ES20.11,A,ES20.11,A,ES20.11,A,ES20.11,A)')&
                  'Point(',NodeNums(j),') = {',&
                  WriteNodes % x(j),',',&
                  WriteNodes % y(j),',',&
                  0.0,',',& !don't need z coord for footprint
                  MeshEdgeLC,'};'
             WritePoints(counter) = NodeNums(j)
             counter = counter + 1
          END DO
          WRITE(GeoUnit,'(A)') ''
       END DO

       !---------------Write lines------------------
       DO i=1,WriteNodeCount-1
          WRITE( GeoUnit,'(A,i0,A,i0,A,i0,A)') 'Line(',WritePoints(i),') = {'&
               ,WritePoints(i),',',WritePoints(i+1),'};'
       END DO
       WRITE( GeoUnit,'(A,i0,A,i0,A,i0,A)') 'Line(',WritePoints(WriteNodeCount),') = {',&
            WritePoints(WriteNodeCount),',',WritePoints(1),'};'

       !------------Write physical lines-------------
       counter = 1
       DO i=1,4 !cycle boundaries
          SELECT CASE(i)
          CASE(1) !left
             NodeNums => LeftLineNodeNums
             MaskName = LeftMaskName
          CASE(2) !back
             NodeNums => BackLineNodeNums
             MaskName = InMaskName
          CASE(3) !right
             NodeNums => RightLineNodeNums
             MaskName = RightMaskName
          CASE(4) !front
             NodeNums => FrontLineNodeNums
             MaskName = FrontMaskName
          END SELECT

          !Find BC number for physical line
          DO j=1,Model % NumberOfBCs
             Found = ListCheckPresent(Model % BCs(j) % Values,MaskName)
             IF(Found) THEN
                BList => ListGetIntegerArray( Model % BCs(j) % Values, &
                     'Target Boundaries', Found )
                IF(SIZE(BList)>1) CALL Fatal(SolverName,&
                     "Could not uniquely determine target BC")
                MeshBC = BList(1)
                EXIT
             END IF
          END DO
          IF(Debug) THEN
             PRINT *, 'Debug CalvingRemesh, BC number for ',MaskName,' is: ',MeshBC
          END IF

          WRITE(GeoUnit,'(A,i0,A)') 'Physical Line(',MeshBC,') = {'
          DO j=1,SIZE(NodeNums)-2
             WRITE(GeoUnit,'(i0,A)') WritePoints(counter),','
             counter = counter + 1
          END DO
          !Last line
          WRITE(GeoUnit,'(i0,A)') WritePoints(counter),'};'
          counter = counter + 1
       END DO

       !--------------Write Line Loop-----------------
       WRITE(GeoUnit, '(A)') 'Line Loop(1) = {'
       DO i=1,WriteNodeCount-1
          WRITE(GeoUnit,'(i0,A)') WritePoints(i),','
       END DO
       WRITE(GeoUnit,'(i0,A)') WritePoints(WriteNodeCount),'};'

       WRITE(GeoUnit,'(A)') 'Plane Surface(1)={1};'
       WRITE(GeoUnit,'(A)') 'Physical Surface(1)={1};'

       !-------------Write attractor etc--------------
       WRITE(GeoUnit,'(A)') 'Field[1] = Attractor;'
       WRITE(GeoUnit,'(A)') 'Field[1].NNodesByEdge = 100.0;'
       WRITE(GeoUnit,'(A)') 'Field[1].NodesList = {'
       DO i=1,SIZE(FrontLineNodeNums)-1
          WRITE(GeoUnit,'(I0,A)') FrontLineNodeNums(i),','
       END DO
       WRITE(GeoUnit,'(I0,A)') FrontLineNodeNums(SIZE(FrontLineNodeNums)),'};'

       WRITE(GeoUnit, '(A)') 'Field[2] = Threshold;'
       WRITE(GeoUnit, '(A)') 'Field[2].IField = 1;'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].LcMin = ',MeshMinLC,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].LcMax = ',MeshMaxLC,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].DistMin = ',MeshMinDist,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].DistMax = ',MeshMaxDist,';'

       WRITE(GeoUnit, '(A)') 'Background Field = 2;'
       WRITE(GeoUnit, '(A)') 'Mesh.CharacteristicLengthExtendFromBoundary = 0;'

       rt = RealTime() - rt0
       IF(ParEnv % MyPE == 0) &
            PRINT *, 'Remesh, Time taken to write footprint mesh: ', rt
       rt0 = RealTime()

       !-----------system call gmsh------------------
       !'env -i' obscures all environment variables, so gmsh doesn't see any
       !MPI stuff and break down.
       CALL EXECUTE_COMMAND_LINE( "env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH gmsh -2 "// filename, .TRUE., ierr )
       IF(ierr > 1) THEN
         IF(ierr == 127) THEN
           CALL Fatal(SolverName, "The 3D Calving implementation depends on GMSH, but this has not been found.")
         END IF
         WRITE(Message, '(A,i0)') "Error executing gmsh, error code: ",ierr
         CALL Fatal(SolverName,Message)
       END IF

       rt = RealTime() - rt0
       IF(ParEnv % MyPE == 0) &
            PRINT *, 'Remesh, Time taken to execute gmsh: ', rt
       rt0 = RealTime()

    END IF

999 CONTINUE

    IF(Boss) THEN
       !-----------system call ElmerGrid------------------
       IF(Parallel) THEN
          WRITE(Message,'(A,A,A,i0,A,i0)') "ElmerGrid 14 2 ",TRIM(filename_root),".msh -metis ",&
               ParEnv % PEs,' ',MetisMethod
       ELSE
          WRITE(Message, '(A,A)') "ElmerGrid 14 2 ",TRIM(filename_root)
       END IF

       CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )
       IF(ierr /= 0) THEN
          WRITE(Message, '(A,i0)') "Error executing ElmerGrid, error code: ",ierr
          CALL Fatal(SolverName,Message)
       END IF

       rt = RealTime() - rt0
       IF(ParEnv % MyPE == 0) &
            PRINT *, 'Remesh, Time taken to execute ElmerGrid: ', rt
       rt0 = RealTime()

    END IF !Boss only

    !Ensure all partitions wait until boss has remeshed
    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    !------------------------------------------------------------
    ! Read in our section of footprint mesh
    ! This is copied from above, and should be modified probably...
    !------------------------------------------------------------

    MeshName = TRIM(filename_root)  !will eventually be set internally, so not important for testing
    MeshDir = ""

    !Read our section of the mesh
    !Note: MeshDir isn't used by LoadMesh2, so keep meshes in WD
    FootprintMesh => LoadMesh2( Model, MeshDir, MeshName, .FALSE., &
         ParEnv % PEs, ParEnv % MyPE)
    FootprintMesh % Name = TRIM(OldMesh % Name)//'_footprint'
    FootprintMesh % OutputActive = .TRUE.
    FootprintMesh % Changed = .TRUE. 

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to read new mesh: ', rt
    rt0 = RealTime()

    !-------------------------------------------------
    ! Check for 101 elements and remesh if found
    !-------------------------------------------------
    BadMesh = .FALSE.
    DO i=FootprintMesh % NumberOfBulkElements+1, &
         FootprintMesh % NumberOfBulkElements + FootprintMesh % NumberOfBoundaryElements

       Element => FootprintMesh % Elements(i)
       IF(GetElementFamily(Element) == 1) THEN
          BadMesh = .TRUE.
          PRINT *, 'PE: ', ParEnv % MyPE,' BAD MESH!' 
          EXIT
       END IF
    END DO

    CALL SParIterAllReduceOR(BadMesh)


    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to check bad mesh: ', rt
    rt0 = RealTime()

    IF(BadMesh .AND. .FALSE.) THEN
       CALL ReleaseMesh(FootprintMesh)
       IF(Boss) THEN
          TriedMetis(MetisMethod+1) = .TRUE.
          IF(.NOT. ALL(TriedMetis)) THEN !try another metis algo
             DO i=5,1,-1
                IF(.NOT. TriedMetis(i)) THEN
                   MetisMethod = i-1
                   EXIT
                END IF
             END DO
             WRITE(Message, '(A,i0)' ) "Resulting mesh contained 101 elements, &
                  &trying again with metis method: ",MetisMethod
             CALL Info(SolverName, Message)

             WRITE(Message,'(A,A)') "rm -r ",TRIM(filename_root)

             CALL EXECUTE_COMMAND_LINE( Message, .FALSE., ierr )
          ELSE
             TriedMetis = .FALSE. 
             AllFail = .TRUE.
             DEALLOCATE(WritePoints)
             MetisMethod = MetisMethodDefault

             WRITE(Message, '(A,i0)' ) "All metis algorithms produced orphan nodes, &
                  &perturbing gmsh parameters."
             CALL Info(SolverName, Message)

             WRITE(Message,'(A,A,A,A,A,A)') "rm -r ",TRIM(filename)," ",&
                  TRIM(filename_root),".msh ",TRIM(filename_root)
             CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )

             MeshMinLC = MeshMinLC * 1.05
             MeshMaxLC = MeshMaxLC * 1.05
             MeshMinDist = MeshMinDist * 0.95
             MeshMaxDist = MeshMaxDist * 1.05
          END IF
       END IF

       CALL MPI_Scatter(AllFail, 1, MPI_LOGICAL, &
            AllFailCatch, 1, MPI_LOGICAL,0,ELMER_COMM_WORLD, ierr)

       IF(AllFailCatch) THEN
          GO TO 1280
       ELSE
          GO TO 999
       END IF
    END IF !bad mesh

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    !Check to ensure no nodes are within the region of tangled
    !columns - if they are, shift them!
    IF(TangleOccurs) THEN

      ALLOCATE(FootprintFrontPerm(FootprintMesh % NumberOfNodes))

      CALL MakePermUsingMask( Model, Solver, FootprintMesh, FrontMaskName, &
           .FALSE., FootprintFrontPerm, dummyint) !<- Number of nodes in this part on calving face

      !Broadcast tangled zones to all partitions
      DO i=1,TangledGroups
        CALL MPI_BCAST( TangledZone(i,:), 2, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        IF(Debug) PRINT *,ParEnv % MyPE, ' tangled zone: ',i,TangledZone(i,:)
      END DO

      DO i=1,FootprintMesh % NumberOfNodes
        IF(FootprintFrontPerm(i) <= 0) CYCLE

        NodeHolder(1) = FootprintMesh % Nodes % x(i)
        NodeHolder(2) = FootprintMesh % Nodes % y(i)
        NodeHolder(3) = FootprintMesh % Nodes % z(i)

        NodeHolder = MATMUL(RotationMatrix, NodeHolder)
        DO j=1,TangledGroups
          IF((NodeHolder(2) <= TangledZone(j,2)) .AND. (NodeHolder(2) >= TangledZone(j,1))) THEN
            PRINT *,ParEnv % MyPE, 'Shifted node ',i,&
                 ' in footprint mesh because its in the tangled zone: ',NodeHolder(2)

            NodeHolder(2) = TangledZone(j,1) - 0.1_dp !TODO: eps param
            NodeHolder = MATMUL(UnRotationMatrix, NodeHolder)
            FootprintMesh % Nodes % x(i) = NodeHolder(1)
            FootprintMesh % Nodes % y(i) = NodeHolder(2)
            FootprintMesh % Nodes % z(i) = NodeHolder(3)
          END IF
        END DO
      END DO

      DEALLOCATE(FootprintFrontPerm)
    END IF

    !----------------------------------------------
    !   Extrude new mesh
    !----------------------------------------------

    ExtrudeLevels = ListGetInteger(Model % Simulation, "Remesh Extruded Mesh Levels", Found, UnfoundFatal=.TRUE.)
    ExtrudedMesh => NULL()
    ExtrudedMesh => MeshExtrude(FootprintMesh, ExtrudeLevels-2)

    !----------------------------------------------------
    ! Interp front position from squished front nodes
    ! onto freshly extruded footprint mesh
    !----------------------------------------------------

    !Save Z and set to ZeroOne for all front nodes
    ALLOCATE(ActualHeight(FaceNodeCount))
    DO i=1, OldMesh % NumberOfNodes
       IF(FrontPerm(i) <= 0) CYCLE
       ActualHeight(FrontPerm(i)) = OldMesh % Nodes % z(i)
       OldMesh % Nodes % z(i) = ZeroOneHeight(FrontPerm(i))
    END DO
    ALLOCATE(WorkPerm(SIZE(FrontPerm)))
    WorkPerm = FrontPerm

    CALL VariableRemove(OldMesh % Variables, "ActualHeight")
    CALL VariableAdd(OldMesh % Variables, OldMesh, Solver, "ActualHeight", 1,&
         ActualHeight, WorkPerm)

    CALL RotateMesh(OldMesh, RotationMatrix)

    NULLIFY(WorkReal, WorkPerm)
    ALLOCATE(WorkReal(FaceNodeCount), WorkPerm(OldMesh % NumberOfNodes))
    WorkPerm = FrontPerm

    DO i=1, OldMesh % NumberOfNodes
       IF(WorkPerm(i) <= 0) CYCLE
       WorkReal(WorkPerm(i)) = OldMesh % Nodes % z(i)
    END DO

    CALL VariableRemove(OldMesh % Variables, "FrontExtent")
    CALL VariableAdd(OldMesh % Variables, OldMesh, Solver, "FrontExtent", 1,&
         WorkReal, WorkPerm)

    ALLOCATE(ExtrudedFrontPerm(ExtrudedMesh % NumberOfNodes))
    CALL MakePermUsingMask( Model, Solver, ExtrudedMesh, FrontMaskName, &
         .FALSE., ExtrudedFrontPerm, dummyint) !<- Number of nodes in this part on calving face

    CALL RotateMesh(ExtrudedMesh, RotationMatrix)

    n = ExtrudedMesh % NumberOfNodes
    NULLIFY(WorkReal, WorkPerm)
    ALLOCATE(WorkReal(n), WorkPerm(n))

    WorkPerm = ExtrudedFrontPerm
    WorkReal = 0.0_dp

    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, "ActualHeight", 1, &
         WorkReal, WorkPerm, .FALSE.)

    NULLIFY(WorkReal, WorkPerm)

    ! --- Add front extent (rotated height) var ---
    ALLOCATE(WorkReal(n), WorkPerm(n))

    WorkPerm = ExtrudedFrontPerm
    WorkReal = 0.0_dp
    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, "FrontExtent", 1, &
         WorkReal, WorkPerm, .FALSE.)
    NULLIFY(WorkReal, WorkPerm)


    ! --- Add Calving variable so we can determine which nodes on --- !
    ! ---       the new front can be shifted vertically           --- !
    ALLOCATE(WorkPerm(n))
    WorkPerm = ExtrudedFrontPerm

    IF(COUNT(WorkPerm > 0) > 0) THEN
      ALLOCATE(WorkReal(COUNT(WorkPerm>0)*3))
    ELSE
      ALLOCATE(WorkReal(SIZE(WorkPerm)*3))
    END IF
    WorkReal = 0.0_dp

    VarName = TRIM(CalvingVarName)

    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, VarName, &
         CalvingVar % DOFs, WorkReal, WorkPerm, .FALSE.)
    
    Var => VariableGet(ExtrudedMesh % Variables,VarName,.TRUE.)
    ALLOCATE(Var % PrevValues(SIZE(WorkReal),SIZE(CalvingVar % PrevValues, 2)))
    Var % PrevValues = 0.0_dp

    DO i=1,3 !add each calving DOF
      WorkReal2 => WorkReal( i::3 )
      WorkName = ComponentName(TRIM(Var % Name),i)

      NULLIFY(WorkPerm); ALLOCATE(WorkPerm(SIZE(ExtrudedFrontPerm)))
      WorkPerm = ExtrudedFrontPerm

      CALL VariableAdd( ExtrudedMesh % Variables, ExtrudedMesh, &
           Solver, WorkName, &
           1, WorkReal2, WorkPerm, .FALSE.)

      WorkVar => VariableGet( ExtrudedMesh % Variables, WorkName, .TRUE. )
      IF(.NOT. ASSOCIATED(WorkVar)) CALL Fatal(SolverName, &
           "Error allocating calving PrevValues.")

      NULLIFY(WorkVar % PrevValues)
      WorkVar % PrevValues => Var % PrevValues( i::3, : )
    END DO

    NULLIFY(WorkReal, WorkPerm)
    
    !Add rotated front height as var to both
    !InterpVarToVarReduced
    ALLOCATE(InterpDim(1)); InterpDim = (/3/);
    CALL ParallelActive(.TRUE.)

    !Create an element mask to mask out elements not on the front.
    !  Find the right BC
    DO i=1,Model % NumberOfBCs
       ThisBC = ListGetLogical(Model % BCs(i) % Values,FrontMaskName,Found)
       IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
       FrontBCtag =  Model % BCs(i) % Tag
       EXIT
    END DO

    !  Create the mask
    ALLOCATE(OldElemMask(OldMesh % NumberOfBulkElements+&
         OldMesh % NumberOfBoundaryElements))

    OldElemMask = .TRUE.
    DO i=OldMesh % NumberOfBulkElements+1,&
         OldMesh % NumberOfBulkElements+OldMesh % NumberOfBoundaryElements
       IF(OldMesh % Elements(i) % BoundaryInfo % Constraint /= FrontBCtag) CYCLE !not on front
       OldElemMask(i) = .FALSE.
    END DO

    CALL InterpolateVarToVarReduced(OldMesh, ExtrudedMesh, "ActualHeight", &
         InterpDim, UnfoundNodes, OldElemMask=OldElemMask, Variables=OldMesh % Variables, &
         GlobalEps=extrude_globaleps, LocalEps=extrude_localeps)

    IF(ANY(UnfoundNodes)) THEN
       DO i=1, SIZE(UnfoundNodes)
          IF(UnfoundNodes(i)) THEN
             PRINT *,ParEnv % MyPE,' Didnt find point: ', i, ' x:', ExtrudedMesh % Nodes % x(i),&
                  ' y:', ExtrudedMesh % Nodes % y(i),&
                  ' z:', ExtrudedMesh % Nodes % z(i)
             CALL InterpolateUnfoundPoint( i, ExtrudedMesh, "ActualHeight", InterpDim,&
                  Variables=ExtrudedMesh % Variables )
          END IF
       END DO
       CALL Warn(SolverName,"Failed to find all nodes in calving front interp")
    END IF

    DEALLOCATE(InterpDim, OldElemMask)

    Var => VariableGet(ExtrudedMesh % Variables, VarName, .TRUE.)
    IF(.NOT. ASSOCIATED(Var)) CALL Fatal(SolverName, &
         "Couldn't get calving var on new mesh to determining calving nodes")

    ALLOCATE(IsCalvingNode(ExtrudedMesh % NumberOfNodes))
    IsCalvingNode = .FALSE.
    DO i=1,ExtrudedMesh % NumberOfNodes
       IF(Var % Perm(i) <= 0) CYCLE
       IF(ANY(Var % Values((Var % Perm(i)*3)-2:(Var % Perm(i)*3)) /= 0.0_dp)) THEN
         IsCalvingNode(i) = .TRUE.
       END IF
    END DO

    ! ---------------------------------------------------
    !           Internal Mesh Update New Mesh
    !         (to ensure proper height interp)
    ! ---------------------------------------------------
    !
    !  Mesh Update 3 = dist between squished non-vertical 
    !                  front and vertical extruded front
    !  Mesh Update 1,2 = 0
    ! ---------------------------------------------------

    ! Find mesh update solver
    MuVarName = ListGetString(Params,"Mesh Update Helper Variable",Found, UnfoundFatal=.TRUE.)

    DO i=1,Model % NumberOfSolvers
       IF(.NOT. ASSOCIATED(Model % Solvers(i) % Variable)) CYCLE
       IF(Model % Solvers(i) % Variable % Name == MuVarName) THEN
          MUSolver => Model % Solvers(i)
          EXIT
       END IF
    END DO
    IF(.NOT. ASSOCIATED(MUSolver)) &
         CALL Fatal("Calving Remesh","Couldn't find Remesh Update variable")

    !Point mesh update to new mesh
    MUSolver % Mesh => ExtrudedMesh
    CALL CopyIntrinsicVars(OldMesh, ExtrudedMesh)

    !Add mesh velocity variable to extruded mesh (this shouldn't be necessary really...)
    n = ExtrudedMesh % NumberOfNodes
    NULLIFY(WorkPerm, WorkReal)
    ALLOCATE(WorkPerm(n), WorkReal(n*dim))
    WorkPerm = [(i,i=1,n)]
    WorkReal = 0.0_dp
    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, MUSolver, "Mesh Velocity",&
         dim, WorkReal, WorkPerm, .FALSE.)

    !--------------------------------------------
    ! Reallocate variable and reconstruct matrix 
    ! based on new mesh (ExtrudedMesh)
    !--------------------------------------------

    NULLIFY(WorkPerm, WorkReal)
    ALLOCATE(WorkPerm(n), WorkReal(n*dim))
    WorkPerm = [(i,i=1,n)]
    WorkReal = 0.0_dp
    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, MUSolver, MuVarName,&
         DIM, WorkReal, WorkPerm)
    NULLIFY(WorkReal, WorkPerm)
    MuSolver % Variable => VariableGet(ExtrudedMesh % Variables, MuVarName, .TRUE.)

    !Reconstruct matrix
    IF(ASSOCIATED(MUSolver % Matrix)) CALL FreeMatrix(MUSolver % Matrix)
    MUSolver % Matrix => CreateMatrix(Model, MUSolver, ExtrudedMesh, &
         MUSolver % Variable % Perm, dim, MATRIX_CRS, .TRUE., &
         ListGetString(MUSolver % Values, "Equation"))

    MUSolver % Matrix % Perm => MUSolver % Variable % Perm
    CALL AllocateVector( MUSolver % Matrix % RHS, MUSolver % Matrix % NumberOfRows )

    Model % Solver => MUSolver
    CALL SetCurrentMesh( Model, ExtrudedMesh)

    !------------------------------------------------
    !  Issue here: 'thin plate' style glacier meshes 
    !    aren't great for mesh update in their long axes.
    !
    !  We want to make it behave better (deformations 
    !  transmitted further upstream). Ideally, specify 
    !  anisotropy in mesh update, but easier to 
    !  temporarily stretch the mesh vertically to improve
    !  its aspect ratio...
    !------------------------------------------------

    !Stretch new mesh vertically to induce pseudo-anisotropy
    MuStretchZ = ListGetConstReal(Params, "Remesh Vertical Stretch", Found, UnfoundFatal=.TRUE.)

    NULLIFY(WorkReal)
    ALLOCATE(WorkReal(ExtrudedMesh % NumberOfNodes))

    CALL RotateMesh(ExtrudedMesh, UnRotationMatrix)
    WorkReal = ExtrudedMesh % Nodes % z !save to avoid FP error
    DO i=1, ExtrudedMesh % NumberOfNodes
       ExtrudedMesh % Nodes % z(i) = ExtrudedMesh % Nodes % z(i) * MuStretchZ
    END DO
    CALL RotateMesh(ExtrudedMesh, RotationMatrix)

    !---------- Call mesh update solver ------------
    CALL SingleSolver( Model, MUSolver, MUSolver % dt, Transient) 

    !Unstretch the mesh
    CALL RotateMesh(ExtrudedMesh, UnRotationMatrix)
    ExtrudedMesh % Nodes % Z = WorkReal
    CALL RotateMesh(ExtrudedMesh, RotationMatrix)
    DEALLOCATE(WorkReal)

    !Put things back
    CALL SetCurrentMesh( Model, OldMesh)
    Model % Solver => Solver
    MuSolver % Variable => VariableGet(OldMesh % Variables, MuVarName, .TRUE.)

    !Unrotate meshes
    CALL RotateMesh(ExtrudedMesh, UnRotationMatrix)
    CALL RotateMesh(OldMesh, UnRotationMatrix)

    !Put nodes back
    DO i=1, OldMesh % NumberOfNodes
       IF(FrontPerm(i) <= 0) CYCLE
       OldMesh % Nodes % z(i) = ActualHeight(FrontPerm(i))
    END DO

    !Now that we've set FrontExtent for Mesh Update dirichlet,
    !  put the nodes back to pre-calving geometry
    CALL DisplaceCalvingFront(OldMesh, CalvingVar, -1)


    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to Extrude Mesh and do front interp: ', rt
    rt0 = RealTime()

    !----------------------------------------------
    !   Interp Top and Bottom Height
    !----------------------------------------------

    !Get pointer to top and bottom vars in old mesh 
    TopVarName = "RemeshTopSurf"
    BottomVarName = "RemeshBottomSurf"

    TopVar => VariableGet(OldMesh % Variables, TopVarName, .TRUE.)
    IF(.NOT.ASSOCIATED(TopVar)) CALL Fatal(SolverName, "Couldn't get variable:&
         &RemeshTopSurf")
    BottomVar => VariableGet(OldMesh % Variables, BottomVarName, .TRUE.)
    IF(.NOT.ASSOCIATED(BottomVar)) CALL Fatal(SolverName, "Couldn't get variable:&
         &RemeshBottomSurf")

    !When the model initialises, these two exported variables
    !share a perm with other vars, so we don't want to deallocate
    !However, once variables are copied to a new mesh, they have 
    !their own perm, and so we deallocate it here to avoid memory leak
    IF(FirstTime) THEN
       NULLIFY(BottomVar % Perm, TopVar % Perm)
    ELSE
       DEALLOCATE(BottomVar % Perm, TopVar % Perm)
    END IF

    n = OldMesh % NumberOfNodes
    ALLOCATE(BottomVarOldPerm(n), TopVarOldPerm(n))

    BottomVar % Perm => BottomVarOldPerm
    TopVar % Perm => TopVarOldPerm

    !mess with botperm before this point
    !and bottomvar % values
    BottomVar % Perm = BotPerm
    TopVar % Perm = TopPerm

    IF(FrontalBecomeBasal) THEN
      NextBasalPerm = MAXVAL(BottomVar % Perm)
      DO i=1,OldMesh % NumberOfNodes
        IF(BottomVar % Perm(i) > 0) CYCLE
        IF(.NOT. NewBasalNode(i)) CYCLE

        NextBasalPerm = NextBasalPerm + 1
        BottomVar % Perm(i) = NextBasalPerm
        IF(Debug) PRINT *,ParEnv % MyPE,' debug, adding basalperm to node: ',i
      END DO

      DEALLOCATE(BottomVar % Values)
      ALLOCATE(BottomVar % Values(MAXVAL(BottomVar % Perm)))
    END IF

    DO i=1,OldMesh % NumberOfNodes
       IF(TopVar % Perm(i) > 0) THEN
          TopVar % Values(TopVar % Perm(i)) = OldMesh % Nodes % z(i)
       END IF
       IF(BottomVar % Perm(i) > 0) THEN
          BottomVar % Values(BottomVar % Perm(i)) = OldMesh % Nodes % z(i)
       END IF
    END DO

    !Add to ExtrudedMesh 
    n = ExtrudedMesh % NumberOfNodes
    ALLOCATE(TopVarValues(n),BottomVarValues(n),TopVarPerm(n),BottomVarPerm(n))
    TopVarPerm = 0; BottomVarPerm = 0;
    NodesPerLevel = n / ExtrudeLevels

    DO i=1,NodesPerLevel
       BottomVarPerm(i) = i
       TopVarPerm(n - NodesPerLevel + i) = i
    END DO

    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, TopVarName, 1, &
         TopVarValues, TopVarPerm, .TRUE.)

    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, BottomVarName, 1, &
         BottomVarValues, BottomVarPerm, .TRUE.)

    !If required, add Grounding Line variable to new mesh to be interpolated.
    !We do this here, instead of in SwitchMesh, because nodes which will be
    !grounded read their z coordinate from the specified bed, rather than the 
    !old mesh, to avoid problems with GL on the next timestep
    IF(DoGL) THEN
      ALLOCATE(WorkPerm(ExtrudedMesh % NumberOfNodes), &
            WorkReal(COUNT(BottomVarPerm > 0)))
       WorkPerm = BottomVarPerm
       WorkReal = 0.0_dp
       CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, GLVarName, 1, &
         WorkReal, WorkPerm, .TRUE.)
       NULLIFY(WorkReal, WorkPerm)
       NewGLVar => VariableGet(ExtrudedMesh % Variables, GLVarName, .TRUE.)
       IF(ASSOCIATED(OldGLVar % PrevValues)) THEN
          ALLOCATE(NewGLVar % PrevValues(SIZE(NewGLVar % Values), SIZE(OldGLVar % PrevValues,2)))
       END IF
    END IF

    IF(ASSOCIATED(InterpDim)) DEALLOCATE(InterpDim)
    ALLOCATE(InterpDim(1)); InterpDim(1) = 3

    CALL InterpolateVarToVarReduced(OldMesh, ExtrudedMesh, TopVarName, InterpDim, UnfoundNodes,&
         GlobalEps=extrude_globaleps, LocalEps=extrude_localeps) 

    IF(ANY(UnfoundNodes)) THEN
      DO i=1, SIZE(UnfoundNodes)
        IF(UnfoundNodes(i)) THEN
          PRINT *,ParEnv % MyPE,' Missing interped point: ', i, &
               ' x:', ExtrudedMesh % Nodes % x(i),&
               ' y:', ExtrudedMesh % Nodes % y(i),&
               ' z:', ExtrudedMesh % Nodes % z(i)
          CALL InterpolateUnfoundPoint( i, ExtrudedMesh, TopVarName, InterpDim )
        END IF
      END DO
      WRITE(Message,'(a,i0,a,i0,a)') "Failed to find ",COUNT(UnfoundNodes),' of ',&
           SIZE(UnfoundNodes),' nodes on top surface for mesh extrusion.'
      CALL Warn(SolverName, Message)
    END IF

    !Avoid accidentally interpolating any other variables (specifically, ActualHeight in
    !cases where a basal element has all 3 nodes on the front)
    !So, we temporarily chop OldGLVar out the variable linked list
    IF(DoGL) THEN
       WorkVar => OldMesh % Variables
       IF(ASSOCIATED(WorkVar, OldGLVar)) THEN
          OldMesh % Variables => OldMesh % Variables % Next
          First = .TRUE.
       ELSE

          DO WHILE(ASSOCIATED(WorkVar))
             IF(ASSOCIATED(WorkVar % Next, OldGLVar)) EXIT
             WorkVar => WorkVar % Next
          END DO

          WorkVar % Next => OldGLVar % Next
          First = .FALSE.
       END IF
       NULLIFY(OldGLVar % Next)
    END IF

    CALL InterpolateVarToVarReduced(OldMesh, ExtrudedMesh, BottomVarName, InterpDim, UnfoundNodes,&
         Variables=OldGLVar, GlobalEps=extrude_globaleps, LocalEps=extrude_localeps)

    !Put GL var back in the linked list
    IF(DoGL) THEN
       IF(First) THEN
          OldGLVar % Next => OldMesh % Variables
          OldMesh % Variables => OldGLVar
       ELSE
          OldGLVar % Next => WorkVar % Next
          WorkVar % Next => OldGLVar
       END IF
    END IF

    IF(ANY(UnfoundNodes)) THEN
       DO i=1, SIZE(UnfoundNodes)
          IF(UnfoundNodes(i)) THEN
             PRINT *,ParEnv % MyPE,'Didnt find point: ', i,&
                  ' frontperm: ',ExtrudedFrontPerm(i),&
                  ' x:', ExtrudedMesh % Nodes % x(i),&
                  ' y:', ExtrudedMesh % Nodes % y(i),&
                  ' z:', ExtrudedMesh % Nodes % z(i)
             CALL InterpolateUnfoundPoint( i, ExtrudedMesh, BottomVarName, &
                  InterpDim, Variables=NewGLVar )
          END IF
       END DO
       WRITE(Message,'(a,i0,a,i0,a)') "Failed to find ",COUNT(UnfoundNodes),' of ',&
            SIZE(UnfoundNodes),' nodes on bottom surface for mesh extrusion.'
       CALL Warn(SolverName, Message)
    END IF

    TopVar => NULL(); BottomVar => NULL()
    TopVar => VariableGet(ExtrudedMesh % Variables, TopVarName, .TRUE.)
    IF(.NOT. ASSOCIATED(TopVar)) CALL Fatal(SolverName, &
         "Couldn't find top surface variable on extruded mesh.")
    BottomVar => VariableGet(ExtrudedMesh % Variables,BottomVarName, .TRUE.)
    IF(.NOT. ASSOCIATED(BottomVar)) CALL Fatal(SolverName, &
         "Couldn't find bottom surface variable on extruded mesh.")

    !Grounded nodes get BottomVar (i.e. nodes % z) from the bed function
    IF(DoGL) THEN
       NewGLVar => VariableGet(ExtrudedMesh % Variables, GLVarName, .TRUE.)
       IF(.NOT. ASSOCIATED(NewGLVar)) CALL Fatal(SolverName,&
            "Trying to account for the grounding line, but can't find GL var on new mesh.")

       ALLOCATE(BedHeight(ExtrudedMesh % NumberOfNodes))
       BedHeight = 0.0_dp
       Material => GetMaterial(ExtrudedMesh % Elements(1)) !TODO, this is not generalised

       CALL SetCurrentMesh( Model, ExtrudedMesh)

       DO i=ExtrudedMesh % NumberOfBulkElements+1, &
            ExtrudedMesh % NumberOfBulkElements+ExtrudedMesh % NumberOfBoundaryElements

          Element => ExtrudedMesh % Elements(i)
          n = Element % TYPE % NumberOfNodes

          IF(ANY(BottomVarPerm(Element % NodeIndexes) <= 0)) CYCLE

          BedHeight(Element % Nodeindexes(1:N)) = &
              ListGetReal(Material,'Min Zs Bottom',n,Element % NodeIndexes, Found, UnfoundFatal=.TRUE.)

       END DO

       CALL SetCurrentMesh( Model, OldMesh)

       PRINT *, ParEnv % MyPE, ' Remesh, max/min bedheight: ', MAXVAL(BedHeight), MINVAL(BedHeight)

       DO i=1,ExtrudedMesh % NumberOfNodes
          IF(NewGLVar % Perm(i) <= 0) CYCLE
          IF(NewGLVar % Values(NewGLVar % Perm(i)) < -0.5) CYCLE !floating, so interped bed is fine
          BottomVar % Values(BottomVar % Perm(i)) = BedHeight(i)
       END DO
    END IF

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to interp top and bottom: ', rt
    rt0 = RealTime()

    !--------------------------------------------------------------------------
    ! Adjust 'ActualHeight' where the removal of nodes on the front has resulted
    !  in a discrepancy between 'ActualHeight' and 'Top (or Bottom) Var'
    ! Three possible cases:
    !   1) Full depth calving - vertical shift constrained by top and bottom node
    !   2) Partial calving reaching top or bottom - 
    !        vertical shift constrained by 1 of top and bottom, and 
    !        the node below/above the last calving node
    !   3) Partial calving reaches neither top nor bottom - ignore, doesn't matter
    !--------------------------------------------------------------------------

    HeightVar => VariableGet(ExtrudedMesh % Variables, "ActualHeight", .TRUE.)
    IF(.NOT.ASSOCIATED(HeightVar)) &
         CALL Fatal(SolverName, "Unable to get ActualHeight var from new mesh")

    DEALLOCATE(FNColumns)
    ALLOCATE(FNColumns(ExtrudedMesh % NumberOfNodes))
    FNColumns = [(i,i=1,ExtrudedMesh % NumberOfNodes)]
    !NodesPerLevel is currently for ExtrudedMesh, see above

    FNColumns = MOD(FNColumns, NodesPerLevel)

    DO WHILE(.TRUE.) !cycle columns

       !Find a new column
       col = -1
       DO j=1, ExtrudedMesh % NumberOfNodes
          IF(ExtrudedFrontPerm(j) <= 0) CYCLE
          IF(FNColumns(j) /= -1) THEN !not done
             col = FNColumns(j)
             EXIT
          END IF
       END DO
       IF(col == -1) EXIT !All done

       !Gather front nodes in specified column
       WorkNodes % NumberOfNodes = COUNT((ExtrudedFrontPerm > 0) .AND. (FNColumns == col))

       IF(WorkNodes % NumberOfNodes /= ExtrudedLevels ) THEN
          WRITE(Message,'(A,i0,A)') "Error in FNColumns, only found ",&
               WorkNodes % NumberOfNodes," nodes in column."
          CALL Fatal(SolverName,Message)
       END IF

       n = WorkNodes % NumberOfNodes
       ALLOCATE(WorkNodes % z(n), ColumnPerm(n))

       counter = 1
       DO j=1, ExtrudedMesh % NumberOfNodes
          IF(ExtrudedFrontPerm(j) <= 0) CYCLE
          IF(FNColumns(j) == col) THEN
             WorkNodes % z(counter) = ExtrudedMesh % Nodes % z(j)
             ColumnPerm(counter) = j
             counter = counter + 1
          END IF
       END DO

       !Order by ascending WorkNodes % z
       ALLOCATE(OrderPerm(n))
       OrderPerm = [(i,i=1,n)]

       CALL SortD( n, WorkNodes % z, OrderPerm )

       start = 1
       DO WHILE(.TRUE.)

          !The difference between front interpolated and top interpolated height
          ! at the top and bottom of the column - set these inside this loop because
          ! they'll be modified in some scenarios
          BotDisplacement = BottomVar % Values(BottomVar % Perm(ColumnPerm(OrderPerm(1)))) - &
               HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(1))))

          TopDisplacement = TopVar % Values(TopVar % Perm(ColumnPerm(OrderPerm(n)))) - &
               HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(n))))

          InGroup = .FALSE.
          GroupCount = 0
          GroupEnd = 0

          !Cycle nodes upwards
          DO k=start,n

             !Add node to group
             IF(IsCalvingNode(ColumnPerm(OrderPerm(k)))) THEN

                IF(.NOT. InGroup) GroupStart = k
                InGroup = .TRUE.

                GroupEnd = k
                GroupCount = GroupCount + 1

                !top node, forces empty DO loop next time, InGroup = .FALSE.
                IF(k == n) start = k + 1

             !Not in group
             ELSE
                IF(InGroup) THEN
                   start = k + 1
                   EXIT
                END IF
             END IF

          END DO

          IF(.NOT. InGroup) EXIT !didn't find any more calving nodes this round, done

          !Which case?:
          IF((GroupStart==1) .AND. (GroupEnd==n)) THEN
             !Top and BotDisplacement are already set above, good to go...
             CONTINUE
          ELSE IF(GroupStart==1) THEN
             TopDisplacement = 0.0_dp
          ELSE IF(GroupEnd==n) THEN
             BotDisplacement = 0.0_dp
          ELSE !do nothing
             CYCLE
          END IF

          DO k=GroupStart, GroupEnd
             IF(GroupCount > 1) THEN
                !How far through the column are we? Drowning in permutations...
                prop = (HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(k)))) - &
                     HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(GroupStart))))) / &
                     (HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(GroupEnd)))) - &
                     HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(GroupStart)))))
             ELSE
                !Special case - only top or bottom node moves
                IF(GroupStart == 1) THEN
                   prop = 0.0_dp
                ELSE
                   prop = 1.0_dp
                END IF
             END IF

             Displacement = (prop * TopDisplacement) + ((1 - prop) * BotDisplacement)

             !Displace this node's heightvar
             HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(k)))) = &
                  HeightVar % Values(HeightVar % Perm(ColumnPerm(OrderPerm(k)))) + &
                  Displacement
          END DO
       END DO

       FNColumns(ColumnPerm) = -1 !mark column nodes as done already

       DEALLOCATE(WorkNodes % z, ColumnPerm, OrderPerm)
    END DO
    DEALLOCATE(FNColumns)

    !------------------------------------------------
    !
    ! Height Solver
    !
    !------------------------------------------------

    n = ExtrudedMesh % NumberOfNodes
    NULLIFY(WorkPerm, WorkReal)
    ALLOCATE(WorkPerm(n), WorkReal(n))
    WorkPerm = [(i,i=1,n)]
    WorkReal = 0.0_dp
    CALL VariableAdd(ExtrudedMesh % Variables, ExtrudedMesh, Solver, &
         Solver % Variable % Name, 1, WorkReal, WorkPerm)

    Solver % Mesh => ExtrudedMesh
    Solver % Variable => VariableGet(ExtrudedMesh % Variables, Solver % Variable % Name, .TRUE.)

    !Reconstruct matrix
    IF(ASSOCIATED(Solver % Matrix)) CALL FreeMatrix(Solver % Matrix)
    Solver % Matrix => CreateMatrix(Model, Solver, ExtrudedMesh, &
         Solver % Variable % Perm, 1, MATRIX_CRS, .TRUE., &
         ListGetString(Solver % Values, "Equation"))

    Solver % Matrix % Perm => Solver % Variable % Perm
    CALL AllocateVector( Solver % Matrix % RHS, Solver % Matrix % NumberOfRows )

    Solver % Matrix % RHS = 0.0_dp
    Solver % Matrix % RHS_im => NULL()

    ALLOCATE(Solver % Matrix % Force(Solver % Matrix % NumberOfRows, Solver % TimeOrder+1))
    Solver % Matrix % Force = 0.0_dp

    ParEnv % ActiveComm = Solver % Matrix % Comm

    n = ExtrudedMesh % MaxElementNodes
    ALLOCATE( FORCE(n), STIFF(n,n))

    Solver % NumberOfActiveElements = 0
    DEALLOCATE(Solver % ActiveElements)
    ALLOCATE(Solver % ActiveElements(Solver % Mesh % NumberOfBulkElements + &
         Solver % Mesh % NumberOfBoundaryElements))

    DO i=1,Solver % Mesh % NumberOfBulkElements+Solver % Mesh % NumberOFBoundaryElements
       CurrentElement => Solver % Mesh % Elements(i)
       IF( CurrentElement % PartIndex /= ParEnv % myPE ) CYCLE

       IF ( CheckElementEquation( Model, CurrentElement, &
            ListGetString(Solver % Values, "Equation")) ) THEN
          Solver % NumberOfActiveElements = Solver % NumberOFActiveElements + 1
          Solver % ActiveElements( Solver % NumberOFActiveElements ) = i
       END IF

    END DO

    CALL DefaultInitialize()

    active = GetNOFActive()
    DO i=1,active
       Element => GetActiveElement(i)
       n = GetElementNOFNodes(Element)
       CALL LocalMatrix(  STIFF, FORCE, Element, n )
       CALL DefaultUpdateEquations( STIFF, FORCE, Element )
    END DO

    CALL DefaultFinishBulkAssembly()
    CALL DefaultFinishAssembly()

    StiffMatrix => Solver % Matrix
    ForceVector => StiffMatrix % RHS

    !-----------------------------------------
    ! Set dirichlet points for height solution
    !   1) Top and Bottom from respective vars
    !   2) Calving front from interpolated "ActualHeight"
    !-----------------------------------------
    !TODO: Possibility to remove other dirichlets from SIF and implement them like this:
    DO i=1, ExtrudedMesh % NumberOfNodes
       IF(HeightVar % Perm(i)>0) THEN
          CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
               Solver % Variable % Perm, i, HeightVar % Values(HeightVar % Perm(i)))
       ELSE IF(TopVar % Perm(i) > 0) THEN
          CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
               Solver % Variable % Perm, i, TopVar % Values(TopVar % Perm(i)))
       ELSE IF(BottomVar % Perm(i) > 0) THEN
          CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
               Solver % Variable % Perm, i, BottomVar % Values(BottomVar % Perm(i)))
       END IF
    END DO

    Norm = DefaultSolve()

    DO i=1,ExtrudedMesh % NumberOfNodes
       ExtrudedMesh % Nodes % z(i) = Solver % Variable % Values(Solver % Variable % Perm(i))
    END DO

    !Put var back to avoid screwing up SwitchMesh
    Solver % Variable => VariableGet(OldMesh % Variables, Solver % Variable % Name, .TRUE.)

    NewMesh => ExtrudedMesh

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to adjust front nodes and HeightSolver: ', rt
    rt0 = RealTime()

    IF(.FALSE.) THEN
      !------------------------------------------------------
      ! The idea here is to identify degenerate elements in the new mesh. If found,
      ! go back to the start of the subroutine, but using the NewMesh as the OldMesh,
      ! and generating a new one. This ISN'T properly implemented - I gave up on it for
      ! now because it's a very rare issue and I have stuff to be getting on with!
      !------------------------------------------------------
      ! If pursuing this idea, it is essential to consider:
      !    DisplaceMesh by calving
      !    GL var, etc for interp...
      !    Var removal/deallocation?
      !    Other deallocations
      !------------------------------------------------------

      !Counter to make sure it doesn't keep trying forever...
      DegenCount = DegenCount + 1
      IF(DegenCount <= 3) THEN

        n = NewMesh % NumberOfBulkElements + NewMesh % NumberOfBoundaryElements
        ALLOCATE(Basis(NewMesh % MaxElementNodes))
        ALLOCATE(Degenerate(n))
        Degenerate = .FALSE.
        AnyDegenerate = .FALSE.

        DO i=1,NewMesh % NumberOfBulkElements + NewMesh % NumberOfBoundaryElements

          Element => NewMesh % Elements(i)
          CALL GetElementNodes(Nodes, Element)

          !Cycle Gauss points, checking for degenerate element.
          !Not sure if/why this needs to be done at the integration
          !points...
          IP = GaussPoints( Element )
          DO j=1,IP % n
            IF(.NOT. ElementInfo( Element, Nodes, IP % U(j), IP % V(j), &
                 IP % W(j),  detJ, Basis )) THEN
              Degenerate(i) = .TRUE.
              EXIT
            END IF
          END DO


        END DO

        AnyDegenerate = COUNT(Degenerate) > 0
        CALL SParIterAllReduceOR(AnyDegenerate)

        IF(AnyDegenerate) THEN

          ALLOCATE(MyDegenerateCoords(COUNT(Degenerate),2))

          counter = 0
          DO i=1,NewMesh % NumberOfBulkElements + NewMesh % NumberOfBoundaryElements
            !Cycle nodes, getting range of effect
            IF(Degenerate(i)) THEN

              counter = counter + 1
              MyDegenerateCoords(counter,1) = HUGE(1.0_dp)
              MyDegenerateCoords(counter,2) = -HUGE(1.0_dp)

              n = Element % TYPE % NumberOfNodes
              DO j=1,n
                NodeHolder(1) = NewMesh % Nodes % x(Element % NodeIndexes(j))
                NodeHolder(2) = NewMesh % Nodes % y(Element % NodeIndexes(j))
                NodeHolder(3) = NewMesh % Nodes % z(Element % NodeIndexes(j))
                NodeHolder = MATMUL(RotationMatrix, NodeHolder)

                PRINT *,ParEnv % MyPE, 'Debug, found degenerate element: ',&
                     i,' with rotated y: ',NodeHolder(2)

                MyDegenerateCoords(counter,1) = MIN(MyDegenerateCoords(counter,1),NodeHolder(2))
                MyDegenerateCoords(counter,2) = MAX(MyDegenerateCoords(counter,2),NodeHolder(2))
              END DO
            END IF
          END DO

          IF(Parallel) THEN
            !MPI Gather count in each part
            IF(Boss) ALLOCATE(PartCountDegenerate(ParEnv % PEs))
            CALL MPI_GATHER(COUNT(Degenerate),1,MPI_INTEGER,PartCountDegenerate,&
                 1,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)


            !MPI GatherV degenerate coords
            IF(Boss) THEN
              ALLOCATE(DegenerateCoords(SUM(PartCountDegenerate),2))

              disps(1) = 0
              DO i=2,PEs
                disps(i) = disps(i-1) + PartCountDegenerate(i-1)
              END DO
            END IF

            CALL MPI_GATHERV(MyDegenerateCoords(:,1),&
                 counter,MPI_DOUBLE_PRECISION,&
                 DegenerateCoords(:,1),PartCountDegenerate,&
                 disps,MPI_DOUBLE_PRECISION,0,ELMER_COMM_WORLD, ierr)

            CALL MPI_GATHERV(MyDegenerateCoords(:,2),&
                 counter,MPI_DOUBLE_PRECISION,&
                 DegenerateCoords(:,2),PartCountDegenerate,&
                 disps,MPI_DOUBLE_PRECISION,0,ELMER_COMM_WORLD, ierr)
          ELSE
            ALLOCATE(DegenerateCoords(SIZE(MyDegenerateCoords,1),2))
            DegenerateCoords = MyDegenerateCoords
          ENDIF

          IF(Boss) THEN
            !Organise the gathered data, including a buffer around either size, probably
            ShiftBuffer = 50.0 !50m buffer either side, TODO parameterise

            !Remove columns from line etc
          END IF
        END IF

        DEALLOCATE(Basis)

        IF(AnyDegenerate) THEN
          !TODO: Take care of:
          !   DisplaceMesh by calving
          !   GL var, etc for interp...
          !   Var removal/deallocation?
          !   Other deallocations

          DEALLOCATE(MyDegenerateCoords)

          IF(Boss) THEN
            DEALLOCATE(PartCountDegenerate)
          END IF

          CALL Warn(SolverName, "Redoing mesh because of degenerate elements.")
          GO TO 8989
        END IF

      ELSE

        CALL Warn(SolverName, "Tried to fix mesh degeneracy 3 times and failed, continuing!")
      END IF !DegenCount <= 3

    END IF !.FALSE. <- not used

    !---------------------------------
    !
    ! Deallocations etc
    !
    !--------------------------------

    !Delete unnecessary meshes
    CALL ReleaseMesh(FootprintMesh)
    DEALLOCATE(FootprintMesh)

    !Remove variables from NewMesh
    CALL ReleaseVariableList(NewMesh % Variables)
    NULLIFY(NewMesh % Variables)

    !free the dummy matrix
    CALL FreeMatrix(Solver % Matrix)

    !---------------------------------------------------------------


    DEALLOCATE( TopPerm, BotPerm, FrontPerm, STIFF, FORCE, &
         ExtrudedFrontPerm, UnfoundNodes )

    DEALLOCATE(MyFaceNodeNums, BedHeight, IsCalvingNode)
    IF(TangleOccurs) DEALLOCATE(LocalTangledNode,LocalCalvingVar)

    IF(Boss) THEN
      ! clear up mesh files
      ! don't think this does anything because system calls above flush
      CLOSE(GeoUnit)
      IF(MoveMesh) THEN

	TimestepVar => VariableGet( OldMesh % Variables, "Timestep", .TRUE. )

	!Make the directory
	WRITE(MoveMeshFullPath,'(A,A,I4.4)') TRIM(MoveMeshDir), &
	     TRIM(filename_root),INT(TimestepVar % Values(1))

	WRITE(Message,'(A,A)') "mkdir -p ",TRIM(MoveMeshFullPath)
	CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )

	WRITE(Message,'(A,A,A,A)') "mv ",TRIM(filename_root),"* ",&
	     TRIM(MoveMeshFullPath)
	CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )

      ELSE
	WRITE(Message,'(A,A,A,A,A,A)') "rm -r ",TRIM(filename)," ",&
	     TRIM(filename_root),".msh ",TRIM(filename_root)
	CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )
      END IF

      !Deallocations
      DEALLOCATE(FaceNodeNums,&
	   BackNodes % x, BackNodes % y, BackNodes % z,&
	   LeftNodes % x, LeftNodes % y, LeftNodes % z,&
	   RightNodes % x, RightNodes % y, RightNodes % z,&
	   FaceNodesT % x, FaceNodesT % y, FaceNodesT % z,&
	   FrontNodes % x, FrontNodes % y, FrontNodes % z)
      IF(TangleOccurs) THEN
	DEALLOCATE(TangledNode,TangledColumn, TangledZone, FrontCalving1Values)
      END IF
    END IF

    !---------------------------------------------------------------

    FirstTime = .FALSE.
    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr) !Wait for the boss

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Remesh, Time taken to tidy up etc: ', rt
    rt0 = RealTime()

  END SUBROUTINE CalvingRemesh



  ! Takes two meshes which are assumed to represent the same domain
  ! and interpolates variables between them. Uses full dimension 
  ! interpolation (InterpolateMeshToMesh) for all nodes, then picks
  ! up missing boundary nodes using reduced dim 
  ! (InterpolateVarToVarReduced)
  SUBROUTINE SwitchMesh(Model, Solver, OldMesh, NewMesh)

    USE CalvingGeometry

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
    !-------------------------------------------------
    TYPE(Solver_t), POINTER :: WorkSolver
    TYPE(Variable_t), POINTER :: Var=>NULL(), NewVar=>NULL(), WorkVar=>NULL()
    TYPE(Valuelist_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: WorkMatrix=>NULL()
    LOGICAL :: Found, Global, GlobalBubbles, Debug, DoPrevValues, &
         NoMatrix, DoOptimizeBandwidth, PrimaryVar, HasValuesInPartition, &
         PrimarySolver
    LOGICAL, POINTER :: UnfoundNodes(:)=>NULL()
    INTEGER :: i,j,k,DOFs, nrows,n
    INTEGER, POINTER :: WorkPerm(:)=>NULL()
    REAL(KIND=dp), POINTER :: WorkReal(:)=>NULL(), WorkReal2(:)=>NULL(), PArray(:,:) => NULL()
    REAL(KIND=dp) :: FrontOrientation(3), RotationMatrix(3,3), UnRotationMatrix(3,3), &
         globaleps, localeps
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, WorkName

    INTERFACE
       SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
            NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
         !------------------------------------------------------------------------------
         USE Lists
         USE SParIterComm
         USE Interpolation
         USE CoordinateSystems
         !-------------------------------------------------------------------------------
         TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
         TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
         LOGICAL, OPTIONAL :: UseQuadrantTree
         LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
         TYPE(Projector_t), POINTER, OPTIONAL :: Projector
         CHARACTER(LEN=*),OPTIONAL :: MaskName
       END SUBROUTINE InterpolateMeshToMesh
    END INTERFACE

    SolverName = "SwitchMesh"
    Debug = .FALSE.
    Params => Solver % Values
    CALL Info( 'Remesher', ' ',Level=4 )
    CALL Info( 'Remesher', '-------------------------------------',Level=4 )
    CALL Info( 'Remesher', ' Switching from old to new mesh...',Level=4 )
    CALL Info( 'Remesher', '-------------------------------------',Level=4 )
    CALL Info( 'Remesher', ' ',Level=4 )

    IF(ASSOCIATED(NewMesh % Variables)) CALL Fatal(SolverName,&
         "New mesh already has variables associated!")

    !interpolation epsilons
    globaleps = global_eps
    localeps = local_eps

    !----------------------------------------------
    ! Get the orientation of the calving front
    ! & compute rotation matrix
    !----------------------------------------------
    PArray => ListGetConstRealArray( Model % Constants,'Front Orientation', &
         Found, UnfoundFatal=.TRUE.)
    DO i=1,3
       FrontOrientation(i) = PArray(i,1)
    END DO
    RotationMatrix = ComputeRotationMatrix(FrontOrientation)
    UnRotationMatrix = TRANSPOSE(RotationMatrix)

    !----------------------------------------------
    !               Action
    !----------------------------------------------

    CALL CopyIntrinsicVars(OldMesh, NewMesh)

    !----------------------------------------------
    ! Add Variables to NewMesh
    !----------------------------------------------

    Var => OldMesh % Variables
    DO WHILE( ASSOCIATED(Var) )

       DoPrevValues = ASSOCIATED(Var % PrevValues)
       WorkSolver => Var % Solver
       HasValuesInPartition = .TRUE.

       !Do nothing if it already exists
       NewVar => VariableGet( NewMesh % Variables, Var % Name, ThisOnly = .TRUE.)
       IF(ASSOCIATED(NewVar)) THEN
          NULLIFY(NewVar)
          Var => Var % Next
          CYCLE
       END IF

       DOFs = Var % DOFs
       Global = (SIZE(Var % Values) .EQ. DOFs)

       !Allocate storage for values and perm
       IF(Global) THEN 
          ALLOCATE(WorkReal(DOFs))
          WorkReal = Var % Values

          CALL VariableAdd( NewMesh % Variables, NewMesh, &
               Var % Solver, TRIM(Var % Name), &
               Var % DOFs, WorkReal)

       ELSE !Regular field variable
          ALLOCATE(WorkPerm(NewMesh % NumberOfNodes))

          IF(.NOT. ASSOCIATED(WorkSolver)) THEN
             WRITE(Message, '(a,a,a)') "Variable ",Var % Name," has no solver, unexpected."
             CALL Fatal(SolverName, Message)
          END IF

          PrimaryVar = ASSOCIATED(WorkSolver % Variable, Var)

          IF(PrimaryVar) THEN !Take care of the matrix
             NoMatrix = ListGetLogical( WorkSolver % Values, 'No matrix',Found)
             !Issue here, this will recreate matrix for every variable associated w/ solver.

             IF(.NOT. NoMatrix) THEN
                IF(ParEnv % MyPE == 0) PRINT *, 'Computing matrix for variable: ',TRIM(Var % Name)

                DoOptimizeBandwidth = ListGetLogical( WorkSolver % Values, &
                     'Optimize Bandwidth', Found )
                IF ( .NOT. Found ) DoOptimizeBandwidth = .TRUE.

                GlobalBubbles = ListGetLogical( WorkSolver % Values, &
                     'Bubbles in Global System', Found )
                IF ( .NOT. Found ) GlobalBubbles = .TRUE.

                WorkMatrix => CreateMatrix(Model, WorkSolver, &
                     NewMesh, WorkPerm, DOFs, MATRIX_CRS, DoOptimizeBandwidth, &
                     ListGetString( WorkSolver % Values, 'Equation' ), &
                     GlobalBubbles = GlobalBubbles )

                IF(ASSOCIATED(WorkMatrix)) THEN
                   WorkMatrix % Comm = ELMER_COMM_WORLD

                   WorkMatrix % Symmetric = ListGetLogical( WorkSolver % Values, &
                        'Linear System Symmetric', Found )

                   WorkMatrix % Lumped = ListGetLogical( WorkSolver % Values, &
                        'Lumped Mass Matrix', Found )

                   CALL AllocateVector( WorkMatrix % RHS, WorkMatrix % NumberOfRows )
                   WorkMatrix % RHS = 0.0_dp
                   WorkMatrix % RHS_im => NULL()

                   ALLOCATE(WorkMatrix % Force(WorkMatrix % NumberOfRows, WorkSolver % TimeOrder+1))
                   WorkMatrix % Force = 0.0_dp
                ELSE
                   !No nodes in this partition now
                   NoMatrix = .TRUE.
                END IF
             END IF

             IF ( ASSOCIATED(Var % EigenValues) ) THEN
                n = SIZE(Var % EigenValues)

                IF ( n > 0 ) THEN
                   WorkSolver % NOFEigenValues = n
                   CALL AllocateVector( NewVar % EigenValues,n )
                   CALL AllocateArray( NewVar % EigenVectors, n, &
                        SIZE(NewVar % Values) ) 

                   NewVar % EigenValues  = 0.0d0
                   NewVar % EigenVectors = 0.0d0
                   IF(.NOT.NoMatrix) THEN
                      CALL AllocateVector( WorkMatrix % MassValues, SIZE(WorkMatrix % Values) )
                      WorkMatrix % MassValues = 0.0d0
                   END IF
                END IF
             END IF

             IF(ASSOCIATED(WorkSolver % Matrix)) CALL FreeMatrix(WorkSolver % Matrix)
             WorkSolver % Matrix => WorkMatrix

             !Check for duplicate solvers with same var
             DO j=1,Model % NumberOfSolvers
                IF(ASSOCIATED(WorkSolver, Model % Solvers(j))) CYCLE
                IF(.NOT. ASSOCIATED(Model % Solvers(j) % Variable)) CYCLE
                IF( TRIM(Model % Solvers(j) % Variable % Name) /= TRIM(Var % Name)) CYCLE
                !Ideally, the solver's old matrix would be freed here, but apart from the 
                !first timestep, it'll be a duplicate
                IF(ASSOCIATED(Model % Solvers(j) % Matrix, WorkMatrix)) CYCLE
                CALL FreeMatrix(Model % Solvers(j) % Matrix)
                Model % Solvers(j) % Matrix => WorkMatrix
             END DO

             NULLIFY(WorkMatrix)

             !NOTE: We don't switch Solver % Variable here, because
             !Var % Solver % Var doesn't necessarily point to self
             !if solver has more than one variable. We do this below.
          ELSE
             k = InitialPermutation(WorkPerm, Model, WorkSolver, &
                  NewMesh, ListGetString(WorkSolver % Values,'Equation'))
          END IF !Primary var

          HasValuesInPartition = COUNT(WorkPerm>0) > 0
          IF(HasValuesInPartition) THEN
             ALLOCATE(WorkReal(COUNT(WorkPerm>0)*DOFs))
          ELSE
             !this is silly but it matches AddEquationBasics
             ALLOCATE(WorkReal(NewMesh % NumberOfNodes * DOFs))
          END IF

          WorkReal = 0.0_dp
          CALL VariableAdd( NewMesh % Variables, NewMesh, &
               Var % Solver, TRIM(Var % Name), &
               Var % DOFs, WorkReal, WorkPerm, &
               Var % Output, Var % Secondary, Var % TYPE )

       END IF !Not global

       NewVar => VariableGet( NewMesh % Variables, Var % Name, ThisOnly = .TRUE. )
       IF(.NOT.ASSOCIATED(NewVar)) CALL Fatal(SolverName,&
            "Problem creating variable on new mesh.")

       IF(DoPrevValues) THEN 
          ALLOCATE(NewVar % PrevValues( SIZE(NewVar % Values), SIZE(Var % PrevValues,2) ))
       END IF

       !Add the components of variables with more than one DOF
       !NOTE, this implementation assumes the vector variable
       !comes before the scalar components in the list.
       !e.g., we add Mesh Update and so here we add MU 1,2,3
       !SO: next time round, new variable (MU 1) already exists
       !and so it's CYCLE'd
       IF((DOFs > 1) .AND. (.NOT.Global)) THEN
          nrows = SIZE(WorkReal)
          DO i=1,DOFs

             WorkReal2 => WorkReal( i:nrows-DOFs+i:DOFs )
             WorkName = ComponentName(TRIM(Var % Name),i)
             CALL VariableAdd( NewMesh % Variables, NewMesh, &
                  Var % Solver, WorkName, &
                  1, WorkReal2, WorkPerm, &
                  Var % Output, Var % Secondary, Var % TYPE )

             IF(DoPrevValues) THEN
                WorkVar => VariableGet( NewMesh % Variables, WorkName, .TRUE. )
                IF(.NOT. ASSOCIATED(WorkVar)) CALL Fatal(SolverName, &
                     "Error allocating Remesh Update PrevValues.")

                NULLIFY(WorkVar % PrevValues)
                WorkVar % PrevValues => NewVar % PrevValues(i:nrows-DOFs+i:DOFs,:)
             END IF

             NULLIFY(WorkReal2)
          END DO
       END IF

       NULLIFY(WorkReal, WorkPerm)
       Var => Var % Next
    END DO

    !set partitions to active, so variable can be -global -nooutput
    CALL ParallelActive(.TRUE.) 
    !MPI_BSend buffer issue in this call to InterpolateMeshToMesh
    CALL InterpolateMeshToMesh( OldMesh, NewMesh, OldMesh % Variables, UnfoundNodes=UnfoundNodes)
    IF(ANY(UnfoundNodes)) THEN
       PRINT *, ParEnv % MyPE, ' missing ', COUNT(UnfoundNodes),' out of ',SIZE(UnfoundNodes),&
            ' nodes in SwitchMesh.'
    END IF

    !---------------------------------------------------------
    ! For top, bottom and calving front BC, do reduced dim 
    ! interpolation to avoid epsilon problems
    !---------------------------------------------------------

    CALL InterpMaskedBCReduced(Model, Solver, OldMesh, NewMesh, OldMesh % Variables, &
         "Top Surface Mask",globaleps=globaleps,localeps=localeps)
    CALL InterpMaskedBCReduced(Model, Solver, OldMesh, NewMesh, OldMesh % Variables, &
         "Bottom Surface Mask",globaleps=globaleps,localeps=localeps)

    CALL RotateMesh(OldMesh, RotationMatrix)
    CALL RotateMesh(NewMesh, RotationMatrix)

    CALL InterpMaskedBCReduced(Model, Solver, OldMesh, NewMesh, OldMesh % Variables, &
         "Calving Front Mask", UnfoundNodes,globaleps=globaleps,localeps=localeps)

    !NOTE: InterpMaskedBCReduced on the calving front will most likely fail to
    ! find a few points, due to vertical adjustment to account for GroundedSolver.
    ! Briefly, the 'DoGL' sections of CalvingRemesh adjust the Z coordinate of
    ! basal nodes which are grounded, to ensure they match the bed dataset.
    ! Thus, it's not impossible for points on the new mesh to sit slightly outside
    ! the old.
    ! However, these points should sit behind or on the old calving front, so
    ! InterpMaskedBC... on the bed should get them. Thus the only thing that may
    ! be missed would be variables defined solely on the front. Currently, none
    ! of these are important for the next timestep, so this should be fine.

    CALL RotateMesh(NewMesh, UnrotationMatrix)
    CALL RotateMesh(OldMesh, UnrotationMatrix)

    !-----------------------------------------------
    ! Point solvers at the correct mesh and variable
    !-----------------------------------------------

    DO i=1,Model % NumberOfSolvers
       WorkSolver => Model % Solvers(i)

       WorkSolver % Mesh => NewMesh !note, assumption here that there's only one active mesh

       !hack to get SingleSolver to recompute
       !should be taken care of by Mesh % Changed, but
       !this is reset by CoupledSolver for some reason
       WorkSolver % NumberOfActiveElements = -1 

       IF(.NOT. ASSOCIATED(WorkSolver % Variable)) CYCLE
       IF(WorkSolver % Variable % NameLen == 0) CYCLE !dummy  !invalid read

       !Check for multiple solvers with same var:
       !If one of the duplicate solvers is only executed before the simulation (or never),
       !then we don't point the variable at this solver. (e.g. initial groundedmask).
       !If both solvers are executed during each timestep, we have a problem.
       !If neither are, it doesn't matter, and so the the later occurring solver will have
       !the variable pointed at it (arbitrary).
       PrimarySolver = .TRUE.
       DO j=1,Model % NumberOfSolvers
          IF(j==i) CYCLE
          IF(.NOT. ASSOCIATED(Model % Solvers(j) % Variable)) CYCLE
          IF(TRIM(Model % Solvers(j) % Variable % Name) == WorkSolver % Variable % Name) THEN

             IF( (WorkSolver % SolverExecWhen == SOLVER_EXEC_NEVER) .OR. &
                  (WorkSolver % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL) ) THEN
                IF((Model % Solvers(j) % SolverExecWhen == SOLVER_EXEC_NEVER) .OR. &
                     (Model % Solvers(j) % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL) ) THEN
                   PrimarySolver = .TRUE.
                ELSE
                   PrimarySolver = .FALSE.
                   WorkSolver % Matrix => NULL()
                   EXIT
                END IF
             ELSE
                IF( (Model % Solvers(j) % SolverExecWhen == SOLVER_EXEC_NEVER) .OR. &
                     (Model % Solvers(j) % SolverExecWhen == SOLVER_EXEC_AHEAD_ALL) ) THEN
                   PrimarySolver = .TRUE.
                   EXIT
                ELSE
                   WRITE(Message, '(A,A)') "Unable to determine main solver for variable: ", &
                        TRIM(WorkSolver % Variable % Name)
                   CALL Fatal(SolverName, Message)
                END IF
             END IF

          END IF
       END DO

       WorkVar => VariableGet(NewMesh % Variables, &
            WorkSolver % Variable % Name, .TRUE.) !invalid read

       IF(ASSOCIATED(WorkVar)) THEN
          WorkSolver % Variable => WorkVar
          IF(PrimarySolver) WorkVar % Solver => WorkSolver
       ELSE
          WRITE(Message, '(a,a,a)') "Variable ",WorkSolver % Variable % Name," wasn't &
               &correctly switched to the new mesh." !invalid read
          PRINT *, i,' debug, solver equation: ', ListGetString(WorkSolver % Values, "Equation")
          CALL Fatal(SolverName, Message)
       END IF

    END DO


    NewMesh % Next => OldMesh % Next
    Model % Meshes => NewMesh
    Model % Mesh => NewMesh
    Model % Variables => NewMesh % Variables

    !Free old mesh and associated variables
    CALL ReleaseMesh(OldMesh)
    DEALLOCATE(OldMesh)
    DEALLOCATE(UnfoundNodes)

    OldMesh => Model % Meshes

  END SUBROUTINE SwitchMesh

  SUBROUTINE InterpMaskedBCReduced(Model, Solver, OldMesh, NewMesh, Variables, MaskName, &
       SeekNodes, globaleps, localeps)

    USE InterpVarToVar

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
    TYPE(Variable_t), POINTER :: Variables
    INTEGER, POINTER :: OldMaskPerm(:)=>NULL(), NewMaskPerm(:)=>NULL()
    INTEGER, POINTER  :: InterpDim(:)
    INTEGER :: dummyint
    REAL(KIND=dp), OPTIONAL :: globaleps,localeps
    REAL(KIND=dp) :: geps,leps
    LOGICAL, POINTER :: OldMaskLogical(:), NewMaskLogical(:), UnfoundNodes(:)=>NULL()
    LOGICAL, POINTER, OPTIONAL :: SeekNodes(:)
    CHARACTER(LEN=*) :: MaskName

    CALL MakePermUsingMask( Model, Solver, NewMesh, MaskName, &
         .FALSE., NewMaskPerm, dummyint)

    CALL MakePermUsingMask( Model, Solver, OldMesh, MaskName, &
         .FALSE., OldMaskPerm, dummyint)

    ALLOCATE(OldMaskLogical(SIZE(OldMaskPerm)),&
         NewMaskLogical(SIZE(NewMaskPerm)))

    OldMaskLogical = (OldMaskPerm <= 0)
    NewMaskLogical = (NewMaskPerm <= 0)
    IF(PRESENT(SeekNodes)) NewMaskLogical = &
         NewMaskLogical .OR. .NOT. SeekNodes

    IF(PRESENT(globaleps)) THEN
      geps = globaleps
    ELSE
      geps = 1.0E-4
    END IF

    IF(PRESENT(localeps)) THEN
      leps = localeps
    ELSE
      leps = 1.0E-4
    END IF

    IF(Debug) PRINT *, ParEnv % MyPE,'Debug, on boundary: ',TRIM(MaskName),' seeking ',&
         COUNT(.NOT. NewMaskLogical),' of ',SIZE(NewMaskLogical),' nodes.'

    ALLOCATE(InterpDim(1))
    InterpDim(1) = 3

    CALL ParallelActive(.TRUE.)
    CALL InterpolateVarToVarReduced(OldMesh, NewMesh, "remesh update 1", InterpDim, &
         UnfoundNodes, OldMaskLogical, NewMaskLogical, Variables=OldMesh % Variables, &
         GlobalEps=geps, LocalEps=leps)

    IF(ANY(UnfoundNodes)) THEN
      !NewMaskLogical changes purpose, now it masks supporting nodes
      NewMaskLogical = (NewMaskPerm <= 0)

      DO i=1, SIZE(UnfoundNodes)
          IF(UnfoundNodes(i)) THEN
             PRINT *,ParEnv % MyPE,'Didnt find point: ', i, &
                  ' x:', NewMesh % Nodes % x(i),&
                  ' y:', NewMesh % Nodes % y(i),&
                  ' z:', NewMesh % Nodes % z(i)

             CALL InterpolateUnfoundPoint( i, NewMesh, "remesh update 1", InterpDim, &
                  NodeMask=NewMaskLogical, Variables=NewMesh % Variables )
          END IF
       END DO

       WRITE(Message, '(i0,a,a,a,i0,a,i0,a)') ParEnv % MyPE,&
            ' Failed to find all points on face: ',MaskName, ', ',&
            COUNT(UnfoundNodes),' of ',COUNT(.NOT. NewMaskLogical),' missing points.'
       CALL Warn("InterpMaskedBCReduced", Message)
    END IF

    DEALLOCATE(OldMaskLogical, &
         NewMaskLogical, NewMaskPerm, &
         OldMaskPerm, UnfoundNodes)

  END SUBROUTINE InterpMaskedBCReduced


  ! Sets the value of coordinate variables from 
  ! a given mesh.
  SUBROUTINE SetCoordVar(Var, Mesh)
    IMPLICIT NONE

    TYPE(Variable_t) :: Var
    TYPE(Mesh_t) :: Mesh
    INTEGER :: int
    CHARACTER(MAX_NAME_LEN) :: str
    str = Var % Name(12:12)
    read( str, '(I1)') int

    SELECT CASE(int) !1,2,3
    CASE(1)
       Var % Values = Mesh % Nodes % x
    CASE(2)
       Var % Values = Mesh % Nodes % y
    CASE(3)
       Var % Values = Mesh % Nodes % z
    CASE DEFAULT
       CALL FATAL("Remesh","Problem setting coordinate variable.  This shoudn't have happened.")
    END SELECT

  END SUBROUTINE SetCoordVar

  SUBROUTINE DisplaceCalvingFront(Mesh, CalvingVar, Sign)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: CalvingVar
    INTEGER :: Sign, k, i

    DO i=1,Mesh % NumberOfNodes
       k = CalvingVar % Perm(i)
       IF(k <= 0) CYCLE
       k = k * CalvingVar % DOFs

       Mesh % Nodes % x(i) = Mesh % Nodes % x(i) + &
            Sign * CalvingVar % Values(k-2)

       Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + &
            Sign * CalvingVar % Values(k-1)

       Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + &
            Sign * CalvingVar % Values(k)
    END DO

  END SUBROUTINE DisplaceCalvingFront

  !Constructs the local matrix for the "d2U/dz2 = 0" Equation
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element)

    FORCE = 0.0_dp
    STIFF = 0.0_dp

    dim = CoordinateSystemDimension()

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx )

       DO p=1,n
          DO q=1,n
             STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
          END DO
       END DO
    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix

  !Taken from TwoMeshes
  !------------------------------------------------------------------------------
  SUBROUTINE SetDirichtletPoint( StiffMatrix, ForceVector,DOF, NDOFs, &
       Perm, NodeIndex, NodeValue) 
    !------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:), NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
    !------------------------------------------------------------------------------

    INTEGER :: PermIndex
    REAL(KIND=dp) :: s

    !------------------------------------------------------------------------------

    PermIndex = Perm(NodeIndex)

    IF ( PermIndex > 0 ) THEN
       PermIndex = NDOFs * (PermIndex-1) + DOF

       IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN        
          CALL SBand_SetDirichlet( StiffMatrix,ForceVector,PermIndex,NodeValue )        
       ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. &
            StiffMatrix % Symmetric ) THEN        
          CALL CRS_SetSymmDirichlet(StiffMatrix,ForceVector,PermIndex,NodeValue)        
       ELSE                          
          s = StiffMatrix % Values(StiffMatrix % Diag(PermIndex))
          ForceVector(PermIndex) = NodeValue * s
          CALL ZeroRow( StiffMatrix,PermIndex )
          CALL SetMatrixElement( StiffMatrix,PermIndex,PermIndex,1.0d0*s )        
       END IF
    END IF

    !------------------------------------------------------------------------------
  END SUBROUTINE SetDirichtletPoint
  !------------------------------------------------------------------------------

END SUBROUTINE Remesher
