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

!Subroutine for predicting calving in 3D

 SUBROUTINE Find_Calving3D_LSet ( Model, Solver, dt, TransientSimulation )

   USE CalvingGeometry
   USE MainUtils
   USE InterpVarToVar
   USE MeshUtils

   IMPLICIT NONE

!-----------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!-----------------------------------------------
   TYPE(ValueList_t), POINTER :: Params, MeshParams => NULL()
   TYPE(Variable_t), POINTER :: CalvingVar, DistVar, CrevVar, &
        SignDistVar, HitCountVar, IsolineIDVar, MeshCrossVar
   TYPE(Solver_t), POINTER :: PCSolver => NULL(), &
        VTUOutputSolver => NULL(), IsoSolver => NULL()
   TYPE(Mesh_t), POINTER :: Mesh, PlaneMesh, IsoMesh, WorkMesh, WorkMesh2
   TYPE(Element_t), POINTER :: Element, WorkElements(:),IceElement
   TYPE(Nodes_t), TARGET :: FaceNodesT
   INTEGER :: i,j,jmin,k,n,dim, dummyint, NoNodes, ierr, PEs, &
        FaceNodeCount, DOFs, PathCount, LeftConstraint, RightConstraint, &
        FrontConstraint, NoCrevNodes, NoPaths, IMBdryCount, &
        node, nodecounter, CurrentNodePosition, StartNode, &
        directions, Counter, ClosestCrev, NumEdgeNodes, UnFoundLoops, EdgeLength,&
        status(MPI_STATUS_SIZE), NUnfoundConstraint
   INTEGER, POINTER :: CalvingPerm(:), TopPerm(:)=>NULL(), BotPerm(:)=>NULL(), &
        LeftPerm(:)=>NULL(), RightPerm(:)=>NULL(), FrontPerm(:)=>NULL(), InflowPerm(:)=>NULL(),&
        FrontNodeNums(:), FaceNodeNums(:)=>NULL(), DistPerm(:), WorkPerm(:), &
        SignDistPerm(:), NodeIndexes(:),IceNodeIndexes(:),&
        EdgeMap(:,:)
   INTEGER, ALLOCATABLE :: CrevEnd(:),CrevStart(:),IMBdryConstraint(:),IMBdryENums(:),&
         PolyStart(:), PolyEnd(:), EdgeLine(:,:), EdgeCount(:), Nodes(:), StartNodes(:,:),&
         WorkInt(:), WorkInt2D(:,:), PartCount(:), ElemsToAdd(:), PartElemsToAdd(:), &
         EdgeLineNodes(:), NodePositions(:), FrontToLateralConstraint(:), UnfoundConstraints(:)
#ifdef ELMER_BROKEN_MPI_IN_PLACE
   INTEGER, ALLOCATABLE :: buffer2(:)
#endif
   REAL(KIND=dp) :: FrontOrientation(3), &
        RotationMatrix(3,3), UnRotationMatrix(3,3), NodeHolder(3), MaxMeshDist,&
        y_coord(2), TempDist,MinDist, xl,xr,yl, yr, xx,yy,&
        angle,angle0,a1(2),a2(2),b1(2),b2(2),a2a1(2),isect(2),front_extent(4), &
        buffer, gridmesh_dx, FrontLeft(2), FrontRight(2), ElemEdge(2,5), &
        CalvingLimit, CrevPenetration, PrevValue, PartMinDist, &
        Calv_delta_t, Calv_dmax, Calv_k, RandomNumber, x, maxprob, Mu, Prob, calv_f, &
        rt0, rt

   REAL(KIND=dp), POINTER :: DistValues(:), SignDistValues(:), WorkReal(:), &
        CalvingValues(:)
   REAL(KIND=dp), ALLOCATABLE :: CrevX(:),CrevY(:),IMBdryNodes(:,:), Polygon(:,:), PathPoly(:,:), &
        EdgeX(:), EdgeY(:), EdgePoly(:,:), CrevOrient(:,:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, DistVarname, &
        FrontMaskName,TopMaskName,BotMaskName,LeftMaskName,RightMaskName,InflowMaskName, &
        PC_EqName, Iso_EqName, VTUSolverName, EqName, FileName
   LOGICAL :: Found, Parallel, Boss, Debug, FirstTime = .TRUE., CalvingOccurs=.FALSE., &
        SaveParallelActive, LeftToRight, inside, Complete, SaveToFile, &
        LatCalvMargins, FullThickness, UnfoundConstraint, RandomCalving, FileCreated
   LOGICAL, ALLOCATABLE :: RemoveNode(:), IMNOnFront(:), IMOnMargin(:), &
        IMNOnLeft(:), IMNOnRight(:), IMElmONFront(:), IMElmOnLeft(:), IMElmOnRight(:), FoundNode(:), &
        IMElemOnMargin(:), DeleteMe(:), IsCalvingNode(:), PlaneEdgeElem(:), EdgeNode(:), UsedElem(:), &
        CrevLR(:), DeleteNode(:), DeleteElement(:)

   TYPE(CrevassePath_t), POINTER :: CrevassePaths, CurrentPath

   SAVE :: FirstTime, SolverName, Params, Parallel, Boss, dim, Debug, &
        DistVarName, PC_EqName, Iso_EqName, LeftConstraint, &
        RightConstraint, FrontConstraint,TopMaskName, BotMaskName, &
        LeftMaskName, RightMaskName, FrontMaskName, InflowMaskName, FileCreated

!---------------Get Variables and Parameters for Solution-----------

   rt0 = RealTime()

   IF(FirstTime) THEN
      SolverName = "Find_Calving3D"
      Params => Solver % Values
      Parallel = (ParEnv % PEs > 1)
      Boss = (ParEnv % MyPE == 0) .OR. (.NOT. Parallel)
      Debug = .TRUE.

      TopMaskName = "Top Surface Mask"
      BotMaskName = "Bottom Surface Mask"
      LeftMaskName = "Left Sidewall Mask"
      RightMaskName = "Right Sidewall Mask"
      FrontMaskName = "Calving Front Mask"
      InflowMaskName = "Inflow Mask"

      dim = CoordinateSystemDimension()
      IF(dim /= 3) CALL Fatal(SolverName, "Solver only works in 3D!")

      DistVarName = ListGetString(Params,"Distance Variable Name", Found)
      IF(.NOT. Found) DistVarName = "Distance"

      PC_EqName = ListGetString(Params,"Project Calving Equation Name",Found, UnfoundFatal=.TRUE.)
      Iso_EqName = ListGetString(Params,"Isosurface Equation Name",Found, UnfoundFatal=.TRUE.)

      !Get the boundary constraint of the left, right & front boundaries
      LeftConstraint = 0; RightConstraint = 0; FrontConstraint = 0
      DO i=1,Model % NumberOfBCs
        IF(ListCheckPresent(Model % BCs(i) % Values,FrontMaskName)) &
             FrontConstraint = Model % BCs(i) % Tag
        IF(ListCheckPresent(Model % BCs(i) % Values,LeftMaskName)) &
             LeftConstraint = Model % BCs(i) % Tag
        IF(ListCheckPresent(Model % BCs(i) % Values,RightMaskName)) &
             RightConstraint = Model % BCs(i) % Tag
      END DO
      IF(Debug) PRINT *,'Front, Left, Right constraints: ',FrontConstraint,LeftConstraint,RightConstraint
   END IF !FirstTime

   RandomCalving = ListGetLogical(Params,"Random Calving")
   IF(RandomCalving) THEN
    IF(Boss) THEN

      Calv_k = ListGetConstReal(Params, "Calving k",Found)
      IF(.NOT. Found) CALL FATAL(SolverName, "Random calving specified but no 'Calving k' value")
      Calv_delta_t = ListGetConstReal(Params, "Calving delta t",Found)
      IF(.NOT. Found) CALL FATAL(SolverName, "Random calving specified but no 'Calving delta t' value")
      Calv_dmax = ListGetConstReal(Params, "Calving dmax",Found)
      IF(.NOT. Found) CALL FATAL(SolverName, "Random calving specified but no 'Calving dmax' value")
      Calv_f = ListGetConstReal(Params, "Calving f",Found)
      IF(.NOT. Found) CALL FATAL(SolverName, "Random calving specified but no 'Calving f' value")

      SaveToFile = ListGetLogical(Params,"Save Probability to File")

      ! calculate dcrev from calving probability
      !x = (dcrev - dmin) / (1 - dmin)
      !x[dcrev < dmin] = 0
      !mu = (x**k) / ( (x**k) + (1-x)**k )
      !prob = 1 - np.exp(-mu * delta_t)

      maxprob = 1.0_dp - EXP(-1.0_dp * Calv_delta_t * Calv_f)

      ! one random number per timestep
      CALL Random_Seed()
      CALL Random_Number(RandomNumber)

      IF(RandomNumber >= maxprob + AEPS) THEN
        CALL INFO(SolverName, 'shortening random number')
        Prob = maxprob
      ELSE
        Prob = RandomNumber
      END IF

      Mu = -LOG(1.0_dp-Prob)/(Calv_delta_t*calv_f)

      CrevPenetration = Calv_dmax - (LOG(1/mu) / Calv_k)

      IF(SaveToFile) THEN
        Filename = ListGetString(Params,"Probability File Name", Found)
        IF(.NOT. Found) THEN
          CALL WARN(SolverName, 'Output file name not given so using Probability.txt')
          Filename = "Probability.txt"
        END IF

        ! write to file
        IF(FileCreated) THEN
          OPEN( 47, FILE=filename, STATUS='UNKNOWN', ACCESS='APPEND')
        ELSE
          OPEN( 47, FILE=filename, STATUS='UNKNOWN')
          WRITE(47, '(A)') "Calving Probability Output File"
          WRITE(47, '(A)') "TimeStep, Time, RandomNumber, Probability, Mu, DCrev"
        END IF

        WRITE(47, *) GetTimestep(), GetTime(), RandomNumber, Prob, Mu, CrevPenetration
        CLOSE(47)
        FileCreated = .TRUE.
      END IF
    ELSE
      CrevPenetration = 0.0_dp
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, CrevPenetration, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)

    PRINT*, 'CrevPenetration: ', CrevPenetration

    WRITE (Message,'(A,F1.2,A)') 'Calving occuring using crevasse penetration:', CrevPenetration, &
        'due to probability.'
    CALL INFO(SolverName, Message)
   ELSE
    CrevPenetration = ListGetConstReal(Params, "Crevasse Penetration",Found, DefValue = 1.0_dp)
    IF(.NOT. Found) CALL Info(SolverName, "No Crevasse Penetration specified so assuming full thickness")
    PRINT*, 'CrevPenetration: ', CrevPenetration
    IF(CrevPenetration > 1 .OR. CrevPenetration < 0) CALL FATAL(SolverName, "Crevasse Penetraion must be between 0-1")
   END IF

   CalvingOccurs = .FALSE.

   Mesh => Model % Mesh

   ! addition of the lateral boundaries when calculating constrictions for crevasses
   ! on the lateral boundaries
   LatCalvMargins = ListGetLogical(Params,"Lateral Calving Margins", DefValue=.TRUE.)

   MaxMeshDist = ListGetConstReal(Params, "Calving Search Distance",Found, UnfoundFatal=.TRUE.)

   CrevPenetration = ListGetConstReal(Params, "Crevasse Penetration",Found, DefValue = 1.0_dp)
   IF(.NOT. Found) CALL Info(SolverName, "No Crevasse Penetration specified so assuming full thickness")
   FullThickness = CrevPenetration == 1.0_dp
   PRINT*, 'CrevPenetration: ', CrevPenetration
   IF(CrevPenetration > 1 .OR. CrevPenetration < 0) CALL FATAL(SolverName, "Crevasse Penetraion must be between 0-1")

   DistVar => VariableGet(Model % Variables, DistVarName, .TRUE., UnfoundFatal=.TRUE.)
   DistValues => DistVar % Values
   DistPerm => DistVar % Perm

   !This solver's variable - holds the levelset value for
   ! CalvingRemeshMMG - negative where calving occurs
   CalvingVar => Solver % Variable
   IF(.NOT.ASSOCIATED(CalvingVar)) &
      CALL Fatal(SolverName, "Find_Calving3D_Lset has no variable!")
   IF(CalvingVar % DOFs /= 1) &
        CALL Fatal(SolverName,"Solver variable has more than one DOF")

   CalvingValues => CalvingVar % Values
   CalvingPerm => CalvingVar % Perm
   DOFs = CalvingVar % DOFs

   NoNodes = Mesh % NumberOfNodes
   ALLOCATE( TopPerm(NoNodes), BotPerm(NoNodes), LeftPerm(NoNodes),&
        RightPerm(NoNodes), FrontPerm(NoNodes), InflowPerm(NoNodes))

   !Generate perms to quickly get nodes on each boundary
   CALL MakePermUsingMask( Model, Solver, Mesh, TopMaskName, &
        .FALSE., TopPerm, dummyint)
   CALL MakePermUsingMask( Model, Solver, Mesh, BotMaskName, &
        .FALSE., BotPerm, dummyint)
   CALL MakePermUsingMask( Model, Solver, Mesh, LeftMaskName, &
        .FALSE., LeftPerm, dummyint)
   CALL MakePermUsingMask( Model, Solver, Mesh, RightMaskName, &
        .FALSE., RightPerm, dummyint)
   CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
        .FALSE., FrontPerm, FaceNodeCount)
   CALL MakePermUsingMask( Model, Solver, Mesh, InflowMaskName, &
        .FALSE., InflowPerm, dummyint)

   !Get the orientation of the calving front, compute rotation matrix
   FrontOrientation = GetFrontOrientation(Model)
   RotationMatrix = ComputeRotationMatrix(FrontOrientation)
   UnRotationMatrix = TRANSPOSE(RotationMatrix)


   DO i=1, Mesh % NumberOfNodes
      IF(InflowPerm(i) > 0) THEN
        IF(ABS(DistValues(DistPerm(i)) * FrontOrientation(2)) <= MaxMeshDist) &
          CALL FATAL(SolverName, "Reduce Calving Search Distance as parts of the inflow &
            boundary have a lower distance.")
      END IF
   END DO

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Time taken for variable loading, making perms etc: ', rt
   rt0 = RealTime()

   !-----------------------------------------------------------
   !                  Action!
   !-----------------------------------------------------------

   !-----------------------------------------------------------
   !
   ! The CIndex "calving criterion":
   ! At its core, the Nye criterion says, in 2D, crevasses exist
   ! if sigma_xx > 0.
   ! Water filled extension: sigma_xx + water_pressure > 0
   !
   ! 3D extension: If any ( planar principal stress + water_pressure) > 0, crevasse exists
   !-----------------------------------------------------------

   !TODO here - rotate the main mesh, or determine the rotated extent
   ! to pass for Grid Min X, etc
   ! Use MaxMeshDist...

   gridmesh_dx = ListGetConstReal(Params, "PlaneMesh Grid Size",Found, DefValue=20.0_dp)
   IF(.NOT. Found) CALL WARN(SolverName, "No PlaneMesh Grid Size set in sif so setting to 20")
   buffer = gridmesh_dx * 4

   front_extent = GetFrontExtent(Mesh, RotationMatrix, DistVar, MaxMeshDist, buffer)

   CALL ListAddConstReal(MeshParams,"Grid Mesh Min X",front_extent(1))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Max X",front_extent(2))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Min Y",front_extent(3))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Max Y",front_extent(4))
   CALL ListAddConstReal(MeshParams,"Grid Mesh dx",gridmesh_dx)


   PRINT *,ParEnv % MyPE,' front extent: ',front_extent

   PlaneMesh => CreateRectangularMesh(MeshParams)

   CALL SetMeshMaxDOFs(PlaneMesh)

   PlaneMesh % Nodes % z = PlaneMesh % Nodes % y
   PlaneMesh % Nodes % y = PlaneMesh % Nodes % x
   PlaneMesh % Nodes % x = 0.0
   CALL RotateMesh(PlaneMesh, UnrotationMatrix)

   PlaneMesh % Name = "calving_plane"
   PlaneMesh % OutputActive = .TRUE.
   PlaneMesh % Changed = .TRUE.

    ! CALL WriteMeshToDisk2(Model, PlaneMesh, ".")

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken to create rectangular mesh: ', rt
    rt0 = RealTime()

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    WorkMesh => Mesh % Next
    Mesh % Next => PlaneMesh

    CALL CopyIntrinsicVars(Mesh, PlaneMesh)

    !----------------------------------------------------
    ! Project Calving Solver to get % fractured
    !----------------------------------------------------

    ! Locate ProjectCalving Solver
    DO i=1,Model % NumberOfSolvers
       IF(GetString(Model % Solvers(i) % Values, 'Equation') == PC_EqName) THEN
          PCSolver => Model % Solvers(i)
          EXIT
       END IF
    END DO
    IF(.NOT. ASSOCIATED(PCSolver)) &
         CALL Fatal("Calving Remesh","Couldn't find Project Calving Solver")
    PCSolver % Mesh => PlaneMesh
    PCSolver % dt = dt
    PCSolver % NumberOfActiveElements = 0 !forces recomputation of matrix, is this necessary?

    n = PlaneMesh % NumberOfNodes
    ALLOCATE(WorkPerm(n), WorkReal(n))
    WorkReal = 0.0_dp
    WorkPerm = [(i,i=1,n)]

    IF(ASSOCIATED(PCSolver % Matrix)) CALL FreeMatrix(PCSolver % Matrix) !shouldn't be required...
    PCSolver % Matrix => CreateMatrix(Model, PCSolver, PlaneMesh, &
         WorkPerm, 1, MATRIX_CRS, .TRUE.)

    !NOTE: ave_cindex is a misnomer. It actually contains the proprtion (0-1) of intact ice...
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "ave_cindex", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal)

    !Surf cindex (surface crevasses reaching waterline)
    ALLOCATE(WorkReal(n))
    WorkReal = 0.0_dp
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "surf_cindex", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal)

    !Basal cindex (surface and basal crevasses meeting)
    ALLOCATE(WorkReal(n))
    WorkReal = 0.0_dp
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "basal_cindex", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal)

    !Helper var for determining edges in isoline
    ALLOCATE(WorkReal(n))
    WorkReal = 0.0_dp
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "isoline id", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal)

    !Variable to hold the number of hits in ProjectCalving
    ALLOCATE(WorkReal(n))
    WorkReal = 0.0_dp
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "hitcount", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal)

    !Variable for edge of Solver % Mesh for polygon creation for levelset sign
    ALLOCATE(WorkReal(n))
    WorkReal = 0.0_dp
    CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "meshcrossover", &
         1, WorkReal, WorkPerm)
    NULLIFY(WorkReal, WorkPerm)

    !Solver variable - crevasse penetration in range 0-1
    CrevVar => VariableGet(PlaneMesh % Variables, "ave_cindex", .TRUE.)
    PCSolver % Variable => CrevVar
    PCSolver % Matrix % Perm => CrevVar % Perm

    !----------------------------------------------------
    ! Run Project Calving solver
    !----------------------------------------------------
    SaveParallelActive = ParEnv % Active(ParEnv % MyPE+1)
    CALL ParallelActive(.TRUE.)

    Model % Solver => PCSolver
    CALL SingleSolver( Model, PCSolver, PCSolver % dt, .FALSE. )
    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken up to project cindex into plane: ', rt
    rt0 = RealTime()

    PlaneMesh % OutputActive = .TRUE.

    Model % Solver => Solver
    Mesh % Next => WorkMesh !Probably null...

    !--------------------------------------------------------------
    ! Isosurface solver to find cindex 0 contour
    !--------------------------------------------------------------

    IF(Boss) THEN

      HitCountVar => VariableGet(PlaneMesh % Variables, "hitcount", .TRUE.,UnfoundFatal=.TRUE.)
      IsolineIdVar => VariableGet(PlaneMesh % Variables, "isoline id", .TRUE.,UnfoundFatal=.TRUE.)
      MeshCrossVar => VariableGet(PlaneMesh % Variables, "meshcrossover", .TRUE.,UnfoundFatal=.TRUE.)

      ! meshcrossover - 2 = glacier hit, 1 = edge, first node outside glacier domain, 0=outside glacier domain
      MeshCrossVar % Values = 2.0
      !Set PlaneMesh exterior nodes to ave_cindex = 1 (no calving)
      !                             and IsolineID = 1 (exterior)
      DO i=1,PlaneMesh % NumberOfNodes
        IF(HitCountVar % Values(HitCountVar % Perm(i)) <= 0.0_dp) THEN
          CrevVar % Values(CrevVar % Perm(i)) = 1.0
          IsolineIDVar % Values(IsolineIDVar % Perm(i)) = 1.0
          MeshCrossVar % Values(MeshCrossVar % Perm(i)) = 0.0
        END IF
      END DO

      !PlaneMesh nodes in elements which cross the 3D domain boundary
      !get ave_cindex = 0 (to enable CrevassePaths/isolines to reach the edge)
      !Also mark elements which cross the edge of the 3D domain (TODO - not used!)
      ALLOCATE(PlaneEdgeElem(PlaneMesh % NumberOfBulkElements), EdgeNode(PlaneMesh % NumberOfNodes))
      PlaneEdgeElem = .FALSE.
      EdgeNode = .FALSE.

      DO i=1, PlaneMesh % NumberOfBulkElements
        NodeIndexes => PlaneMesh % Elements(i) % NodeIndexes
        n = PlaneMesh % Elements(i) % TYPE % NumberOfNodes

        !Some but not all element nodes have hitcount == 0
        IF((.NOT. ALL(HitCountVar % Values(HitCountVar % Perm(NodeIndexes(1:n))) <= 0.0_dp)) .AND. &
             ANY(HitCountVar % Values(HitCountVar % Perm(NodeIndexes(1:n))) <= 0.0_dp)) THEN

          PlaneEdgeElem(i) = .TRUE.

          !Mark nodes just outside the edge to CrevVar = 0.0
          DO j=1,n
            IF(HitCountVar % Values(HitCountVar % Perm(NodeIndexes(j))) <= 0.0_dp) THEN
                 CrevVar % Values(CrevVar % Perm(NodeIndexes(j))) = 0.0
                 MeshCrossVar % Values(MeshCrossVar % Perm(NodeIndexes(j))) = 1.0
                 EdgeNode(NodeIndexes(j)) = .TRUE.
            END IF
          END DO

        END IF
      END DO
      PRINT *,'Debug - PlaneEdgeElem: ',COUNT(PlaneEdgeElem), 'EdgeNode: ', COUNT(EdgeNode)

      ! get planmesh edgenodes in a line- nodes just outside the solver%mesh domain-
      ! for calving lset sign. These nodes used to make calving polygon
      !! Elem 404 Structure
      !  4____1
      !  |    |
      !  |____|
      !  3    2
      ! following code assumes this structure in order to find edge nodes in order
      ! create edgeline so we can stitch this with crevasses

      ALLOCATE(UsedElem(PlaneMesh % NumberOfBulkElements), StartNodes(2, 5), NodePositions(3))
      UsedElem=.FALSE.
      StartNodes=0
      DO startnode=1, PlaneMesh % NumberOfNodes
        IF(.NOT. EdgeNode(startnode)) CYCLE
        !StartNode = !random start point
        directions=0

        ! find how many starting directions
        DO i=1, PlaneMesh % NumberOfBulkElements
          IF(.NOT. PlaneEdgeElem(i)) CYCLE
          NodeIndexes => PlaneMesh % Elements(i) % NodeIndexes
          n = PlaneMesh % Elements(i) % TYPE % NumberOfNodes
          IF(.NOT. ANY(NodeIndexes(1:n)==startnode)) CYCLE
          NodePositions=0
          nodecounter=0
          DO j=1,n
            IF(.NOT. EdgeNode(NodeIndexes(j))) CYCLE ! not a edge node
            IF(NodeIndexes(j) == startnode) THEN
              CurrentNodePosition = j
              CYCLE ! this node
            END IF
            nodecounter=nodecounter+1
            NodePositions(nodecounter)=j
          END DO

          ! assign nodes found so we can create edgeline in two directions
          IF(nodecounter == 0) THEN
            CYCLE !only contains startnode
          ELSE IF(nodecounter == 1) THEN
            directions=directions+1
            StartNodes(directions, 1) = NodeIndexes(NodePositions(1))
          ELSE IF(nodecounter == 2) THEN
            IF(NodePositions(2)-NodePositions(1) == 2) THEN
              directions=directions+2
              StartNodes(1, 1)=NodeIndexes(NodePositions(1))
              StartNodes(2, 1)=NodeIndexes(NodePositions(2))
            ELSE
              directions=directions+1
              IF(ABS(NodePositions(1)-CurrentNodePosition) /= 2) THEN ! closest node
                StartNodes(directions, 1)=NodeIndexes(NodePositions(1))
                StartNodes(directions, 2)=NodeIndexes(NodePositions(2))
              ELSE
                StartNodes(directions, 1)=NodeIndexes(NodePositions(2))
                StartNodes(directions, 2)=NodeIndexes(NodePositions(1))
              END IF
            END IF
          END IF
          UsedElem(i)=.TRUE.
        END DO

        NumEdgeNodes = COUNT(EdgeNode)
        ! add startnodes to edgeline
        ALLOCATE(Edgeline(directions, NumEdgeNodes), nodes(directions), &
                EdgeCount(directions), FoundNode(directions))
        EdgeCount = 0
        EdgeLine = 0
        DO i=1, directions
          EdgeCount(i) = COUNT(StartNodes(i,:) /= 0)
          IF(i==1) THEN ! for initial startnode
            EdgeCount(i) = EdgeCount(i) + 1
            EdgeLine(i, 1) = startnode
            EdgeLine(i, 2:EdgeCount(i)) = PACK(StartNodes(i,:), StartNodes(i,:) /= 0)
          ELSE
            EdgeLine(i, 1:EdgeCount(i)) = PACK(StartNodes(i,:), StartNodes(i,:) /= 0)
          END IF

          nodes(i) = EdgeLine(i, EdgeCount(i))
        END DO

        ! loop through starting with startnodes to find remaining edgenodes
        ! in order
        UnFoundLoops = 0
        Complete = .FALSE.
        DO WHILE(.NOT. Complete)
          FoundNode=.FALSE.
          DO i=1, PlaneMesh % NumberOfBulkElements
            IF(.NOT. PlaneEdgeElem(i)) CYCLE
            IF(UsedElem(i)) CYCLE ! already used elem
            NodeIndexes => PlaneMesh % Elements(i) % NodeIndexes
            n = PlaneMesh % Elements(i) % TYPE % NumberOfNodes
            DO j=1, directions
              node = Nodes(j)
              IF(.NOT. ANY(NodeIndexes(1:n)==node)) CYCLE
              nodecounter=0
              NodePositions=0
              DO k=1,n
                IF(.NOT. EdgeNode(NodeIndexes(k))) CYCLE ! not a edge node
                IF(NodeIndexes(k) == node) THEN
                  CurrentNodePosition = k
                  CYCLE ! this node
                END IF
                nodecounter=nodecounter+1
                NodePositions(nodecounter) = k
              END DO
              IF(nodecounter == 1) THEN
                IF(ABS(CurrentNodePosition - NodePositions(1)) == 2) CYCLE ! nodes not connected
                IF(ANY(EdgeLine(j, :) == NodeIndexes(NodePositions(1)))) THEN
                  ! increase capactiy by one
                  ALLOCATE(WorkInt2D(directions, NumEdgeNodes))
                  WorkInt2D = EdgeLine
                  DEALLOCATE(EdgeLine)
                  NumEdgeNodes=NumEdgeNodes+1
                  ALLOCATE(EdgeLine(directions, NumEdgeNodes))
                  EdgeLine=0
                  EdgeLine(:,1:NumEdgeNodes-1) = WorkInt2D
                  DEALLOCATE(WorkInt2D)
                END IF
                EdgeCount(j) = EdgeCount(j) + 1
                EdgeLine(j, EdgeCount(j)) = NodeIndexes(NodePositions(1))
              ELSE IF(nodecounter == 2) THEN
                counter=0
                DO k=1, nodecounter
                  IF(ANY(EdgeLine(j, :) == NodeIndexes(NodePositions(k)))) THEN
                    ! increase capactiy by one
                    ALLOCATE(WorkInt2D(directions, NumEdgeNodes))
                    WorkInt2D = EdgeLine
                    DEALLOCATE(EdgeLine)
                    NumEdgeNodes=NumEdgeNodes+1
                    ALLOCATE(EdgeLine(directions, NumEdgeNodes))
                    EdgeLine=0
                    EdgeLine(:,1:NumEdgeNodes-1) = WorkInt2D
                    DEALLOCATE(WorkInt2D)
                  END IF
                  IF(ABS(CurrentNodePosition - NodePositions(k)) /= 2) THEN ! first node
                    EdgeLine(j, EdgeCount(j) + 1 + Counter) = NodeIndexes(NodePositions(k))
                  ELSE ! second node
                    EdgeLine(j, EdgeCount(j) + 1 + Counter) = NodeIndexes(NodePositions(k))
                  END IF
                  counter=counter+1
                END DO
                EdgeCount(j) = EdgeCount(j) + Counter
              ELSE IF(nodecounter == 3) THEN
                CALL FATAL('PlaneMesh', '4 edge nodes in one element!')
              END IF
              Nodes(j)=Edgeline(j, EdgeCount(j))
              IF(nodecounter>0) THEN
                UsedElem(i)=.TRUE.
                FoundNode(j)=.TRUE.
              END IF
            END DO
          END DO
          IF(.NOT. ANY(FoundNode)) THEN
            DO i=1, directions
              IF(FoundNode(i)) CYCLE
              Nodes(i)=Edgeline(i, EdgeCount(i)-1)
            END DO
            UnFoundLoops=UnFoundLoops+1
          ELSE
            UnFoundLoops=0
          END IF
          IF(UnFoundLoops == 2) THEN
            ! if can't find all edge nodes
            ! check if unfound node lies in poly of edgeline
            ALLOCATE(EdgePoly(2, SUM(EdgeCount)+1))
            nodecounter=0
            IF(directions > 1) THEN
              DO i=EdgeCount(2),1, -1 ! backwards iteration
                nodecounter=nodecounter+1
                EdgePoly(1,nodecounter) = PlaneMesh % Nodes % x(EdgeLine(2, i))
                EdgePoly(2,nodecounter) = PlaneMesh % Nodes % y(EdgeLine(2, i))
              END DO
            END IF
            DO i=1, EdgeCount(1)
              nodecounter=nodecounter+1
              EdgePoly(1,nodecounter) = PlaneMesh % Nodes % x(EdgeLine(1, i))
              EdgePoly(2,nodecounter) = PlaneMesh % Nodes % y(EdgeLine(1, i))
            END DO
            EdgePoly(1, SUM(EdgeCount)+1) = EdgePoly(1,1)
            EdgePoly(2, SUM(EdgeCount)+1) = EdgePoly(2,1)

            counter=0
            DO i=1, PlaneMesh % NumberOfNodes
              inside=.FALSE.
              IF(.NOT. EdgeNode(i)) CYCLE
              IF(ANY(EdgeLine == i)) CYCLE
              xx = PlaneMesh % Nodes % x(i)
              yy = PlaneMesh % Nodes % y(i)
              inside = PointInPolygon2D(EdgePoly, (/xx, yy/))
              IF(inside) counter=counter+1
            END DO
            IF(NumEdgeNodes == SUM(EdgeCount) + counter) Complete = .TRUE.
            DEALLOCATE(EdgePoly)
          END IF
          IF(NumEdgeNodes == SUM(EdgeCount)) Complete = .TRUE.
          IF(UnFoundLoops == 10) CALL FATAL('Calving_lset', 'Cannot find all edge nodes')
        END DO
        EXIT ! initial loop only to get random start point
      END DO
      DEALLOCATE(UsedElem, NodePositions)

      EdgeLength = SUM(EdgeCount)
      ! translate edge nodes into x and y and join both strings
      ALLOCATE(EdgeX(EdgeLength), EdgeY(EdgeLength), EdgeLineNodes(EdgeLength))
      nodecounter=0
      IF(directions > 1) THEN
        DO i=EdgeCount(2),1, -1 ! backwards iteration
          nodecounter=nodecounter+1
          EdgeX(nodecounter) = PlaneMesh % Nodes % x(EdgeLine(2, i))
          EdgeY(nodecounter) = PlaneMesh % Nodes % y(EdgeLine(2, i))
          EdgeLineNodes(nodecounter) = EdgeLine(2,i)
        END DO
      END IF
      DO i=1, EdgeCount(1)
        nodecounter=nodecounter+1
        EdgeX(nodecounter) = PlaneMesh % Nodes % x(EdgeLine(1, i))
        EdgeY(nodecounter) = PlaneMesh % Nodes % y(EdgeLine(1, i))
        EdgeLineNodes(nodecounter) = EdgeLine(1,i)
      END DO

      ! see if all mesh nodes lie with edgeline if not adjust
      ALLOCATE(EdgePoly(2, EdgeLength+1))
      EdgePoly(1,1:EdgeLength) = EdgeX
      EdgePoly(1,EdgeLength+1) = EdgeX(1)
      EdgePoly(2,1:EdgeLength) = EdgeY
      EdgePoly(2,EdgeLength+1) = EdgeY(1)
      ! end of finding planmesh edge
    END IF ! BOSS

    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
    CALL MPI_BCAST(EdgeLength, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
    IF(.NOT. Boss) ALLOCATE(EdgePoly(2,EdgeLength+1))
    CALL MPI_BCAST(EdgePoly, (EdgeLength+1)*2, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)

    ! check all nodes lie within edgeline
    ALLOCATE(ElemsToAdd(Mesh % NumberOfNodes))
    counter=0; ElemsToAdd = 0
    DO i=1, Mesh % NumberOfNodes
      IF(DistValues(DistPerm(i)) > MaxMeshDist/2) CYCLE
      inside = PointInPolygon2D(EdgePoly, (/Mesh % Nodes % x(i), Mesh % Nodes % y(i)/))
      IF(.NOT. inside) THEN
        DO j=1, PlaneMesh % NumberOfBulkElements
          NodeIndexes => PlaneMesh % Elements(j) % NodeIndexes
          DO k=1,4
            ElemEdge(1,k) = PlaneMesh % Nodes % x(NodeIndexes(k))
            ElemEdge(2,k) = PlaneMesh % Nodes % y(NodeIndexes(k))
          END DO
          ElemEdge(1,5) = PlaneMesh % Nodes % x(NodeIndexes(1))
          ElemEdge(2,5) = PlaneMesh % Nodes % y(NodeIndexes(1))
          inside = PointInPolygon2D(ElemEdge, (/Mesh % Nodes % x(i), Mesh % Nodes % y(i)/))
          IF(inside) THEN
            IF(ANY(ElemsToAdd == j)) CYCLE
            counter=counter+1
            ElemsToAdd(counter) = j
          END IF
        END DO
      END IF
    END DO

    ALLOCATE(PartCount(ParEnv % PEs))
    PartCount=0
    CALL MPI_GATHER(counter, 1, MPI_INTEGER, PartCount, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

    IF(.NOT. Boss .AND. Counter > 0) THEN
      CALL MPI_BSEND(ElemsToAdd, counter, MPI_INTEGER, 0, &
          8000, ELMER_COMM_WORLD, ierr)
    END IF
    IF(Boss.AND. SUM(PartCount) > 0) THEN
      ALLOCATE(PartElemsToAdd(SUM(PartCount)))
      IF(counter > 0) PartElemsToAdd(1:Counter) = ElemsToAdd(1:counter)
      DO i=2, ParEnv % PEs ! boss is zero
        IF(PartCount(i) == 0) CYCLE
        PRINT*, i-1, counter, PartCount(i)
        CALL MPI_RECV(PartElemsToAdd(counter+1:counter+PartCount(i)), PartCount(i), &
            MPI_INTEGER, i-1, 8000, ELMER_COMM_WORLD, status, ierr)
        counter=counter+PartCount(i)
      END DO

    END IF

    IF(Boss .AND. SUM(PartCount) > 0) THEN
      !remove duplicates
      ALLOCATE(WorkInt(SUM(PartCount)), RemoveNode(SUM(PartCount)))
      WorkInt = PartElemsToAdd
      RemoveNode = .FALSE.
      DO i=1,SUM(PartCount)
        IF(RemoveNode(i)) CYCLE
        DO j=1,SUM(PartCount)
          IF(i==j) CYCLE
          IF(WorkInt(i) == WorkInt(j)) RemoveNode(j) = .TRUE.
        END DO
      END DO
      DEALLOCATE(PartElemsToAdd)
      PartElemsToAdd = PACK(WorkInt, .NOT. RemoveNode)
      DEALLOCATE(WorkInt, RemoveNode)

      !if nodes lie outwith edgeline expand edgeline
      ALLOCATE(UsedElem(SIZE(PartElemsToAdd)), NodePositions(4))
      UsedElem=.FALSE.
      DO WHILE(.NOT. ALL(UsedElem))
        DO i=1, SIZE(PartElemsToAdd)
          IF(UsedElem(i)) CYCLE
          NodeIndexes => PlaneMesh % Elements(PartElemsToAdd(i)) % NodeIndexes
          counter=0
          DO j=1,4
            IF(ANY(EdgeLineNodes == NodeIndexes(j))) counter=counter+1
          END DO
          IF(counter == 1) THEN
            ! add three missing nodes in order
            DO j=1,EdgeLength
              IF(ANY(NodeIndexes == EdgeLineNodes(j))) THEN ! found first node
                DO k=1,4
                  IF(NodeIndexes(k) == EdgeLineNodes(j)) NodePositions(1) = k
                END DO
                ALLOCATE(WorkInt(EdgeLength+4))
                WorkInt(1:j) = EdgeLineNodes(1:j)
                DO k=1,4 ! 4 as need to add original node back in at end
                  n = NodePositions(1) + k
                  IF(n > 4) n = n - 4
                  WorkInt(j+k) = NodeIndexes(n)
                END DO
                WorkInt(j+5:EdgeLength+4) = EdgeLineNodes(j+1:EdgeLength)
                DEALLOCATE(EdgeLineNodes)
                EdgeLength = EdgeLength + 4
                ALLOCATE(EdgeLineNodes(EdgeLength))
                EdgeLineNodes = WorkInt
                DEALLOCATE(WorkInt)
                EXIT
              END IF
            END DO
          ELSE IF(Counter == 2) THEN
            ! add two missing nodes in order
            DO j=1,EdgeLength
              IF(ANY(NodeIndexes == EdgeLineNodes(j))) THEN ! found first node
                DO k=1,4
                  IF(NodeIndexes(k) == EdgeLineNodes(j)) NodePositions(1) = k
                  IF(NodeIndexes(k) == EdgeLineNodes(j+1)) NodePositions(2) = k
                END DO
                IF(NodePositions(1) == NodePositions(2)) CYCLE
                IF(ABS(NodePositions(1) - NodePositions(2)) == 2) &
                  CALL FATAL('Calving3D_lset', 'Error building edgeline')
                ! add in other nodes
                ALLOCATE(WorkInt(EdgeLength+2))
                WorkInt(1:j) = EdgeLineNodes(1:j)
                ! first node is opposite NodePostiions(1)
                DO k=1,4
                  IF(ANY(NodePositions(1:2) == k)) CYCLE
                  IF(ABS(k-NodePositions(1)) /= 2) WorkInt(j+1) = NodeIndexes(k)
                  IF(ABS(k-NodePositions(2)) /= 2) WorkInt(j+2) = NodeIndexes(k)
                END DO
                WorkInt(j+3:EdgeLength+2) = EdgeLineNodes(j+1:EdgeLength)
                DEALLOCATE(EdgeLineNodes)
                EdgeLength = EdgeLength + 2
                ALLOCATE(EdgeLineNodes(EdgeLength))
                EdgeLineNodes = WorkInt
                DEALLOCATE(WorkInt)
                EXIT
              END IF
            END DO
          ELSE IF(counter == 3) THEN !swap middle node for unsed node
            DO j=1,EdgeLength
              IF(ANY(NodeIndexes == EdgeLineNodes(j))) THEN ! found first node
                IF(ANY(EdgeLineNodes(j+1:j+2) == EdgeLineNodes(j))) CYCLE ! loops back on itself
                IF(.NOT. ANY(NodeIndexes == EdgeLineNodes(j+1))) CYCLE ! moves out of element
                IF(.NOT. ANY(NodeIndexes == EdgeLineNodes(j+2))) CYCLE ! moves out of element
                DO k=1,4
                  IF(NodeIndexes(k) == EdgeLineNodes(j)) NodePositions(1) = k
                  IF(NodeIndexes(k) == EdgeLineNodes(j+1)) NodePositions(2) = k
                  IF(ANY(EdgeLineNodes(j:j+2) == NodeIndexes(k))) CYCLE
                  NodePositions(4) = k
                END DO
                ! swap node 2 with 4
                EdgeLineNodes(j+1) = NodeIndexes(NodePositions(4))
                EXIT
              END IF
            END DO
          ELSE IF(counter == 4) THEN ! remove two middle nodes
            DO j=1,EdgeLength
              IF(ANY(NodeIndexes == EdgeLineNodes(j))) THEN ! found first node
                IF(ANY(EdgeLineNodes(j+1:j+3) == EdgeLineNodes(j))) CYCLE ! loops back on itself
                IF(ANY(EdgeLineNodes(j+2:j+3) == EdgeLineNodes(j+1))) CYCLE ! loops back on itself
                IF(EdgeLineNodes(j+2) == EdgeLineNodes(j+3)) CYCLE ! loops back on itself
                IF(.NOT. ANY(NodeIndexes == EdgeLineNodes(j+1))) CYCLE ! moves out of element
                IF(.NOT. ANY(NodeIndexes == EdgeLineNodes(j+2))) CYCLE ! moves out of element
                IF(.NOT. ANY(NodeIndexes == EdgeLineNodes(j+3))) CYCLE ! moves out of element
                DO k=1,4
                  IF(NodeIndexes(k) == EdgeLineNodes(j)) NodePositions(1) = k
                  IF(NodeIndexes(k) == EdgeLineNodes(j+3)) NodePositions(4) = k
                END DO
                IF(ABS(NodePositions(1) - NodePositions(4)) == 2) &
                  CALL FATAL('Calving3D_lset', 'Error building edgeine')
                ALLOCATE(WorkInt(EdgeLength-2))
                WorkInt(1:j) = EdgeLineNodes(1:j)
                WorkInt(j+1:) = EdgeLineNodes(j+3:)
                DEALLOCATE(EdgeLineNodes)
                EdgeLength = EdgeLength - 2
                ALLOCATE(EdgeLineNodes(EdgeLength))
                EdgeLineNodes = WorkInt
                DEALLOCATE(WorkInt)
                EXIT
              END IF
            END DO
          END IF
          UsedElem(i)=.TRUE.
          IF(counter == 0) CALL WARN(SolverName, 'Element to be added not connected')
        END DO
      END DO
      ! update edge line
      DEALLOCATE(EdgeX, EdgeY)
      ALLOCATE(EdgeX(EdgeLength), EdgeY(EdgeLength))
      DO i=1, EdgeLength
        EdgeX(i) = PlaneMesh % Nodes % x(EdgeLineNodes(i))
        EdgeY(i) = PlaneMesh % Nodes % y(EdgeLineNodes(i))
      END DO
    END IF

    IF(Debug .AND. Boss) THEN
      PRINT*, 'EdgeLine ----'
      DO i=1, EdgeLength
        PRINT*, EdgeX(i), ',', EdgeY(i), ','
      END DO
    END IF

    IF(Boss) THEN

        IF(.NOT. FullThickness) THEN
          ! save original ave_cindex
          n = PlaneMesh % NumberOfNodes
          ALLOCATE(WorkPerm(n), WorkReal(n))
          WorkReal = CrevVar % Values(CrevVar % Perm)
          WorkPerm = [(i,i=1,n)]
          CALL VariableAdd(PlaneMesh % Variables, PlaneMesh, PCSolver, "ave_cindex_fullthickness", &
          1, WorkReal, WorkPerm)
          NULLIFY(WorkReal, WorkPerm)

          CalvingLimit = 1.0_dp - CrevPenetration

          !alter ave_cindix so relflect %crevasses indicated in sif
          ! 0 full
          DO i=1, n
            IF(CrevVar % Values(CrevVar % Perm(i)) <= CalvingLimit) THEN
              CrevVar % Values(CrevVar % Perm(i)) = 0.0_dp
            ELSE
              PrevValue = CrevVar % Values(CrevVar % Perm(i))
              CrevVar % Values(CrevVar % Perm(i)) = PrevValue - CalvingLimit
            END IF
          END DO
        END IF



       ! Locate Isosurface Solver
       DO i=1,Model % NumberOfSolvers
          IF(GetString(Model % Solvers(i) % Values, 'Equation') == Iso_EqName) THEN
             IsoSolver => Model % Solvers(i)
             EXIT
          END IF
       END DO
       IF(.NOT. ASSOCIATED(IsoSolver)) &
            CALL Fatal("Calving Remesh","Couldn't find Isosurface Solver")

       IsoSolver % Mesh => PlaneMesh
       IsoSolver % dt = dt
       IsoSolver % NumberOfActiveElements = 0 !forces recomputation of matrix, is this necessary?

       PEs = ParEnv % PEs
       ParEnv % PEs = 1

       Model % Solver => IsoSolver
       CALL SingleSolver( Model, IsoSolver, IsoSolver % dt, .FALSE. )

       Model % Solver => Solver
       ParEnv % PEs = PEs

       !Immediately following call to Isosurface solver, the resulting
       !isoline mesh is the last in the list. We want to remove it from the
       !Model % Mesh linked list and attach it to PlaneMesh
       WorkMesh => Model % Mesh
       DO WHILE(ASSOCIATED(WorkMesh % Next))
          WorkMesh2 => WorkMesh
          WorkMesh => WorkMesh % Next
       END DO
       IsoMesh => WorkMesh
       WorkMesh2 % Next => NULL() !break the list
       PlaneMesh % Next => IsoMesh !add to planemesh

       ALLOCATE(DeleteNode(PlaneMesh % NumberOfNodes), DeleteElement(PlaneMesh % NumberOfBulkElements))
       DeleteNode = .FALSE.; DeleteElement = .FALSE.
       !remove unwanted elements
       DO i=1, PlaneMesh % NumberOfNodes
         IF(CrevVar % Values(CrevVar % Perm(i)) /= CrevPenetration) CYCLE !all nodes outside should equal 1
         inside = PointInPolygon2D(EdgePoly, (/PlaneMesh % Nodes % x(i), PlaneMesh % Nodes % y(i)/))
         IF(inside) CYCLE
         DeleteNode(i) = .TRUE.
       END DO

       DO i=1, PlaneMesh % NumberOfBulkElements
         IF(ANY(DeleteNode(PlaneMesh % Elements(i) % NodeIndexes))) &
           DeleteElement(i) = .TRUE.
       END DO

       CALL CutPlaneMesh(PlaneMesh, DeleteNode, DeleteElement)

       !-------------------------------------------------
       !
       ! Compare rotated isomesh to current calving front
       ! position to see if and where calving occurs
       !
       !-------------------------------------------------

       !-------------------------------------------------
       ! Need to map boundary info from the main
       ! 3D mesh to the IsoMesh. We ask Isosurface Solver
       ! to interpolate the 'Isoline ID' var, which is
       ! 1 where the PlaneMesh didn't hit the 3D mesh, and
       ! 0 elsewhere. This does not allow us to distinguish
       ! between frontal and lateral Isomesh elements/nodes
       ! but it identifies candidate boundary elements in
       ! Isomesh.
       !-------------------------------------------------

       ALLOCATE(IMNOnFront(IsoMesh % NumberOfNodes), &
            IMNOnLeft(IsoMesh % NumberOfNodes),&
            IMNOnRight(IsoMesh % NumberOfNodes),&
            IMOnMargin(IsoMesh % NumberOfNodes))
       IMNOnFront=.FALSE.; IMOnMargin=.FALSE.
       IMNOnLeft=.FALSE.; IMNOnRight=.FALSE.

       IsolineIdVar => VariableGet(IsoMesh % Variables, "isoline id", .TRUE.,UnfoundFatal=.TRUE.)
       DO i=1, IsoMesh % NumberOfNodes
         ImOnMargin(i) = IsolineIDVar % Values(IsolineIDVar % Perm(i)) > 0.0_dp
       END DO

       !Count the elements which cross through the ice domain boundary
       !and fill IMBdryNodes to send to other partitions
       ALLOCATE(IMElemOnMargin(IsoMesh % NumberOfBulkElements))
       IMElemOnMargin = .FALSE.

       DO k=1,2 !2 passes, count then fill
         IMBdryCount = 0
         DO i=1, IsoMesh % NumberOfBulkElements
           NodeIndexes => Isomesh % Elements(i) % NodeIndexes
           IF(IMOnMargin(NodeIndexes(1)) .EQV. IMOnMargin(NodeIndexes(2))) CYCLE
           IMElemOnMargin(i) = .TRUE.
           IMBdryCount = IMBdryCount + 1

           IF(k==2) THEN
             IMBdryNodes(IMBdryCount,1) = Isomesh % Nodes % x(NodeIndexes(1))
             IMBdryNodes(IMBdryCount,2) = Isomesh % Nodes % y(NodeIndexes(1))
             IMBdryNodes(IMBdryCount,3) = Isomesh % Nodes % x(NodeIndexes(2))
             IMBdryNodes(IMBdryCount,4) = Isomesh % Nodes % y(NodeIndexes(2))
             IMBdryENums(IMBdryCount) = i
           END IF
         END DO
         IF(k==1) ALLOCATE(IMBdryNodes(IMBdryCount,4),IMBdryENums(IMBdryCount))
       END DO

     END IF

     IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
     ! check front boundary is connected returns FrontToLateralConstraint which
     ! reassigns unconnected elems to nearest lateral margin in IMBdryConstraints
     CALL CheckFrontBoundary(Model, FrontConstraint, RightConstraint, LeftConstraint, FrontToLateralConstraint)

     !Pass isoline boundary elements to all partitions to get info
     !about which boundary they cross

     CALL MPI_BCAST(IMBdryCount,1,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     IF(.NOT. Boss) ALLOCATE(IMBdryNodes(IMBdryCount,4))

     CALL MPI_BCAST(IMBdryNodes,IMBdryCount*4, MPI_DOUBLE_PRECISION,&
          0, ELMER_COMM_WORLD, ierr)

     ALLOCATE(IMBdryConstraint(IMBdryCount))
     IMBdryConstraint = 0
#ifdef ELMER_BROKEN_MPI_IN_PLACE
     ALLOCATE(buffer2(IMBdryCount))
#endif

     !Now cycle elements: for those with a node either side
     !of domain boundary, cycle 3d mesh boundary elements
     !looking for a 2D intersection. This identifies the
     !physical boundary for the isomesh node.
     DO i=1, IMBdryCount

       a1(1) = IMBdryNodes(i,1)
       a1(2) = IMBdryNodes(i,2)
       a2(1) = IMBdryNodes(i,3)
       a2(2) = IMBdryNodes(i,4)

       !Slightly extend the length of the element to account for floating point
       !precision issues:
       a2a1 = (a2 - a1) * 0.1
       a1 = a1 - a2a1
       a2 = a2 + a2a1

       DO j=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
            Mesh % NumberOfBoundaryElements

         IceElement => Mesh % Elements(j)
         IF(IceElement % BoundaryInfo % Constraint /= FrontConstraint .AND. &
              IceElement % BoundaryInfo % Constraint /= LeftConstraint .AND. &
              IceElement % BoundaryInfo % Constraint /= RightConstraint) CYCLE

         IceNodeIndexes => IceElement  % NodeIndexes

         EdgeMap => GetEdgeMap(IceElement % TYPE % ElementCode / 100)

         !set K = element number of edges
         ! cycle K, setting edge to b1 b2, checking
         n = SIZE(EdgeMap,1)
         DO k=1,n
           b1(1) = Mesh % Nodes % x(IceNodeIndexes(EdgeMap(k,1)))
           b1(2) = Mesh % Nodes % y(IceNodeIndexes(EdgeMap(k,1)))
           b2(1) = Mesh % Nodes % x(IceNodeIndexes(EdgeMap(k,2)))
           b2(2) = Mesh % Nodes % y(IceNodeIndexes(EdgeMap(k,2)))

           CALL LineSegmentsIntersect ( a1, a2, b1, b2, isect, Found)

           IF(Found) THEN
             !set the element BC id here!
             CALL Assert(ASSOCIATED(IceElement % BoundaryInfo), SolverName,&
                  "Boundary element has no % boundaryinfo!")
             IMBdryConstraint(i) = IceElement % BoundaryInfo % Constraint
             IF(Debug) PRINT *,ParEnv % MyPE,' isomesh element ',i,' hit ',j,'constraint: ', &
                  IceElement % BoundaryInfo % Constraint,' xy: ',a1
             EXIT
           END IF
         END DO
         ! when lateral margin advances it doesn't follow the z axis so we want to determine
         ! the further point the lateral margins have advanced
         IF(Found .AND. FrontToLateralConstraint(j) /= 0) &
            IMBdryConstraint(i) = FrontToLateralConstraint(j)
         IF(Found .AND. IMBdryConstraint(i) /= FrontConstraint) EXIT
       END DO
     END DO

     !Send info back to boss
     IF(Boss) THEN
#ifdef ELMER_BROKEN_MPI_IN_PLACE
       buffer2 = IMBdryConstraint
       CALL MPI_Reduce(buffer2, &
#else
       CALL MPI_Reduce(MPI_IN_PLACE, &
#endif
            IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
            MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
     ELSE
       CALL MPI_Reduce(IMBdryConstraint, IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
            MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
     END IF

     CALL GetFrontCorners(Model, Solver, FrontLeft, FrontRight)

     IF(Boss) THEN

       UnfoundConstraint = .FALSE.
       IF(ANY(IMBdryConstraint == 0)) THEN
         PRINT *,'Debug constraints: ', IMBdryConstraint
         CALL WARN(SolverName,"Failed to identify boundary constraint for all isomesh edge elements")
         UnfoundConstraints = PACK((/ (i,i=1,IMBdryCount) /),IMBdryConstraint == 0)
         UnfoundConstraint = .TRUE.
       END IF
     END IF

     CALL MPI_BCAST(UnfoundConstraint, 1, MPI_LOGICAL, 0, ELMER_COMM_WORLD, ierr)

     IF(UnfoundConstraint) THEN
       CALL INFO(SolverName, "Unfound boundary constraints so searching for nearest boundary element")

       IF(Boss) NUnfoundConstraint = SIZE(UnfoundConstraints)
       CALL MPI_BCAST(NUnfoundConstraint, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
       IF(.NOT. Boss) ALLOCATE(UnfoundConstraints(NUnfoundConstraint))
       CALL MPI_BCAST(UnfoundConstraints, NUnfoundConstraint, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

       DO node=1, NUnfoundConstraint

         MinDist = HUGE(1.0_dp)
         i = UnfoundConstraints(node)

         a1(1) = IMBdryNodes(i,1)
         a1(2) = IMBdryNodes(i,2)
         a2(1) = IMBdryNodes(i,3)
         a2(2) = IMBdryNodes(i,4)

         !Slightly extend the length of the element to account for floating point
         !precision issues:
         a2a1 = (a2 - a1) * 0.1
         a1 = a1 - a2a1
         a2 = a2 + a2a1

         DO j=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
              Mesh % NumberOfBoundaryElements

           IceElement => Mesh % Elements(j)
           IF(IceElement % BoundaryInfo % Constraint /= FrontConstraint .AND. &
                IceElement % BoundaryInfo % Constraint /= LeftConstraint .AND. &
                IceElement % BoundaryInfo % Constraint /= RightConstraint) CYCLE

           IceNodeIndexes => IceElement  % NodeIndexes

           n = IceElement % TYPE % NumberOfNodes

           DO k=1,n
             b1(1) = Mesh % Nodes % x(IceNodeIndexes(k))
             b1(2) = Mesh % Nodes % y(IceNodeIndexes(k))

             tempdist = PointLineSegmDist2D ( a1, a2, b1)

             IF(TempDist < MinDist) THEN
               MinDist = TempDist
               jmin = j
             END IF
           END DO
         END DO

         CALL MPI_AllReduce(MinDist, PartMinDist, 1, MPI_DOUBLE_PRECISION, &
            MPI_MIN, ELMER_COMM_WORLD, ierr)

         IF(MinDist == PartMinDist) THEN
            IMBdryConstraint(i) = Mesh % Elements(jmin) % BoundaryInfo % Constraint
            PRINT*, ParEnv % MyPE, 'Assigning constraint', IMBdryConstraint(i), &
               'based off mindist:', PartMinDist
         END IF
       END DO

       IF(Boss) THEN
#ifdef ELMER_BROKEN_MPI_IN_PLACE
         buffer2 = IMBdryConstraint
         CALL MPI_Reduce(buffer2, &
#else
         CALL MPI_Reduce(MPI_IN_PLACE, &
#endif
              IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
              MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
       ELSE
         CALL MPI_Reduce(IMBdryConstraint, IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
              MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
       END IF
     END IF


     IF(Boss) THEN
       !-----------------------------------------------------------------
       ! Cycle elements, deleting any which lie wholly on a margin
       !-----------------------------------------------------------------

       ALLOCATE(DeleteMe( IsoMesh % NumberOfBulkElements ))
       DeleteMe = .FALSE.
       DO i=1, IsoMesh % NumberOfBulkElements
          Element => IsoMesh % Elements(i)
          N = Element % TYPE % NumberOfNodes
          IF(ALL(IMOnMargin(Element % NodeIndexes(1:N)))) THEN
             !delete the element, we don't need it
             DeleteMe(i) = .TRUE.
          END IF
       END DO

       PRINT *, 'debug, ', COUNT(DeleteMe), ' elements marked for deletion from IsoMesh.'

       ALLOCATE(WorkElements(COUNT(.NOT. DeleteMe)))
       WorkElements = PACK(IsoMesh % Elements, (.NOT. DeleteMe))

       IF(Debug) THEN
         DO i=1, SIZE(WorkElements)
           PRINT *, 'Workelements', i, ' nodeindexes: ', &
                WorkElements(i) % NodeIndexes
         END DO
       END IF

       DO i=1, IsoMesh % NumberOfBulkElements
          IF(DeleteMe(i)) CALL DeallocateElement(IsoMesh % Elements(i))
       END DO

       DEALLOCATE(DeleteMe)
       IF(ASSOCIATED(IsoMesh % Elements)) DEALLOCATE(IsoMesh % Elements)
       IsoMesh % Elements => WorkElements
       IsoMesh % NumberOfBulkElements = SIZE(WorkElements)
       NULLIFY(WorkElements)

       Counter=0
       DO i=1,IsoMesh % NumberOfBulkElements
          NodeIndexes => Isomesh % Elements(i) % NodeIndexes
          IF(IMOnMargin(NodeIndexes(1)) .EQV. IMOnMargin(NodeIndexes(2))) CYCLE
          counter = counter + 1
          IMBdryENums(Counter) = i
       END DO
       IF(counter /= IMBdryCount) THEN
          IMBdryCount = Counter
          ALLOCATE(WorkInt(IMBdryCount))
          WorkInt = IMBdryENums(1:IMBdryCount)
          DEALLOCATE(IMBdryENums)
          ALLOCATE(IMBdryENums(IMBdryCount))
          IMBdryENums = WorkInt
          DEALLOCATE(WorkInt)
       END IF

       ALLOCATE(IMElmOnFront(IsoMesh % NumberOfBulkElements), &
       IMElmOnLeft(IsoMesh % NumberOfBulkElements),&
       IMElmOnRight(IsoMesh % NumberOfBulkElements))
       IMElmOnFront=.FALSE.; IMElmOnLeft=.FALSE.; IMElmOnRight=.FALSE.

       ! remember IMOnMargin - nodes, IMBdryConstrain - Elems
       DO i=1,IMBdryCount
         k = IMBdryENums(i)
         PRINT*, i, k, IMBdryConstraint(i)
         IF(k==0) CALL FATAL(SolverName, 'IMBdryConstraint = zero')
         IF(IMBdryConstraint(i) == FrontConstraint) IMElmOnFront(k) = .TRUE.
         IF(IMBdryConstraint(i) == LeftConstraint) IMElmOnLeft(k) = .TRUE.
         IF(IMBdryConstraint(i) == RightConstraint) IMElmOnRight(k) = .TRUE.
         NodeIndexes => Isomesh % Elements(k) % NodeIndexes
         DO j=1,2
          IF(.NOT. IMOnMargin(NodeIndexes(j))) CYCLE
          IF(IMBdryConstraint(i) == FrontConstraint) IMNOnFront(NodeIndexes(j)) = .TRUE.
          IF(IMBdryConstraint(i) == LeftConstraint) IMNOnLeft(NodeIndexes(j)) = .TRUE.
          IF(IMBdryConstraint(i) == RightConstraint) IMNOnRight(NodeIndexes(j)) = .TRUE.
         END DO
       END DO

       IF(Debug) THEN
          PRINT *, 'debug, count IMNOnFront: ', COUNT(IMNOnFront), 'IMElmOnFront', COUNT(IMElmOnFront)
          PRINT *, 'debug, count IMNOnLeft: ', COUNT(IMNOnLeft), 'IMElmOnLeft', COUNT(IMElmOnLeft)
          PRINT *, 'debug, count IMNOnRight: ', COUNT(IMNOnRight), 'IMElmOnRight', COUNT(IMElmOnRight)
          PRINT *, 'debug, count IMOnMargin: ', COUNT(IMOnMargin)
          PRINT *, 'debug, isomesh bulkelements,', IsoMesh % NumberOfBulkElements
          PRINT *, 'debug, isomesh boundaryelements,', IsoMesh % NumberOfBoundaryElements
          PRINT *, 'debug, size isomesh elements: ', SIZE(IsoMesh % Elements)
       END IF

       !Find chains which make contact with front twice
       !===============================================
       ! NOTES:
       ! Because intersections lie exactly on nodes, IsoSurface
       ! creates several nodes per real node, and several 202 elems.
       ! This probably isn't a problem, but we'll see...
       ! Also, we have deleted unnecessary elements, but their nodes
       ! remain, but this also isn't a problem as InterpVarToVar
       ! cycles elements, not nodes.
       CALL FindCrevassePaths(IsoMesh, IMOnMargin, CrevassePaths, PathCount)
       CALL CheckCrevasseNodes(IsoMesh, CrevassePaths, IMNOnLeft, IMNOnRight)
       !CALL ValidateCrevassePaths(IsoMesh, CrevassePaths, FrontOrientation, PathCount,&
       !     IMOnLeft, IMOnRight, .FALSE.)
       IF(Debug) THEN
          PRINT*, 'Debug: Crevs PreChecks'
          counter=0
          CurrentPath => CrevassePaths
          DO WHILE(ASSOCIATED(CurrentPath))
              counter=counter+1
              PRINT*, counter, ',',counter+CurrentPath % NumberOfNodes-1, ','
              counter=counter+CurrentPath % NumberOfNodes-1
              CurrentPath => CurrentPath % Next
          END DO
          CurrentPath => CrevassePaths
          DO WHILE(ASSOCIATED(CurrentPath))
              DO i=1, CurrentPath % NumberOfNodes
                PRINT*, IsoMesh % Nodes % x(CurrentPath % NodeNumbers(i)),',',&
                  IsoMesh % Nodes % y(CurrentPath % NodeNumbers(i)), ','
              END DO
              CurrentPath => CurrentPath % Next
          END DO
       END IF

       ! do not remove crevs inside each other as the bigger crev may be severely constricted
       CALL RemoveInvalidCrevs(IsoMesh, CrevassePaths, EdgeX, EdgeY, .FALSE., .FALSE., &
                      IMNOnleft, IMNOnRight, IMNOnFront, gridmesh_dx)
       CALL ValidateNPCrevassePaths(IsoMesh, CrevassePaths, IMNOnLeft, IMNOnRight, &
                      FrontLeft, FrontRight, EdgeX, EdgeY, LatCalvMargins, gridmesh_dx)
       ! this call to remove crevs within other crevs
       CALL RemoveInvalidCrevs(IsoMesh, CrevassePaths, EdgeX, EdgeY, .TRUE., .TRUE., GridSize=gridmesh_dx)

       IF(Debug) THEN
          PRINT*, 'Debug: Crevs PostChecks'
          counter=0
          CurrentPath => CrevassePaths
          DO WHILE(ASSOCIATED(CurrentPath))
              counter=counter+1
              PRINT*, counter, ',',counter+CurrentPath % NumberOfNodes-1, ','
              counter=counter+CurrentPath % NumberOfNodes-1
              CurrentPath => CurrentPath % Next
          END DO
          CurrentPath => CrevassePaths
          DO WHILE(ASSOCIATED(CurrentPath))
              DO i=1, CurrentPath % NumberOfNodes
                PRINT*, IsoMesh % Nodes % x(CurrentPath % NodeNumbers(i)),',',&
                  IsoMesh % Nodes % y(CurrentPath % NodeNumbers(i)), ','
              END DO
              CurrentPath => CurrentPath % Next
          END DO
       END IF

       IF(Debug) THEN
          PRINT *,'Crevasse Path Count: ', PathCount
          CurrentPath => CrevassePaths
          DO WHILE(ASSOCIATED(CurrentPath))
             PRINT *, 'New Crevasse Path:'
             PRINT *, 'Current path elems:', CurrentPath % ElementNumbers
             PRINT *, 'Current path nodes:', CurrentPath % NodeNumbers
             DO i=1, CurrentPath % NumberOfNodes
                PRINT *,'path ID', CurrentPath % ID, 'path node: ',i,&
                     ' x: ',IsoMesh % Nodes % x(CurrentPath % NodeNumbers(i)),&
                     ' y: ',IsoMesh % Nodes % y(CurrentPath % NodeNumbers(i))
             END DO
             PRINT *, ''
             CurrentPath => CurrentPath % Next
          END DO
       END IF

       ALLOCATE(DeleteMe(IsoMesh % NumberOfBulkElements))
       DeleteMe = .FALSE.

       DO i=1, IsoMesh % NumberOfBulkElements
          IF(ElementPathID(CrevassePaths, i) == 0) DeleteMe(i) = .TRUE.
       END DO

       IF(Debug) THEN
          WRITE (Message,'(A,i0,A)') "Deleting ",COUNT(DeleteMe)," elements which &
               &don't lie on valid crevasse paths."
          CALL Info(SolverName, Message)
       END IF

       ALLOCATE(WorkElements(COUNT(.NOT. DeleteMe)))
       WorkElements = PACK(IsoMesh % Elements, (.NOT. DeleteMe))

       DO i=1, IsoMesh % NumberOfBulkElements
          IF(DeleteMe(i)) CALL DeallocateElement(IsoMesh % Elements(i))
       END DO

       DEALLOCATE(DeleteMe)
       DEALLOCATE(IsoMesh % Elements)
       IsoMesh % Elements => WorkElements
       IsoMesh % NumberOfBulkElements = SIZE(WorkElements)
       NULLIFY(WorkElements)

       !InterpVarToVar only looks at boundary elements. In reality, the
       !two are equivalent here because all are 202 elements.
       IsoMesh % NumberOfBoundaryElements = IsoMesh % NumberOfBulkElements
       IsoMesh % NumberOfBulkElements = 0
       IsoMesh % MaxElementNodes = 2
    ELSE
       !not boss, allocate dummy IsoMesh
       !all front partitions will ask root for interp
       IsoMesh => AllocateMesh()
       ALLOCATE(IsoMesh % Nodes % x(0),&
            IsoMesh % Nodes % y(0),&
            IsoMesh % Nodes % z(0))
    END IF !Boss


    ! calculate signed distance here
    ALLOCATE(WorkReal(NoNodes))
    WorkReal = 0.0_dp
    CALL VariableAdd(Mesh % Variables, Mesh, Solver, &
         "SignedDistance", 1, WorkReal, CalvingPerm) ! TO DO check this is correct perm
    NULLIFY(WorkReal)
    SignDistVar => VariableGet(Mesh % Variables, "SignedDistance", .TRUE.)
    SignDistValues => SignDistVar % Values
    SignDistPerm => SignDistVar % Perm
    CurrentPath => CrevassePaths

    IF(Debug) THEN
       PRINT *, ' IsoMesh % NumberOfNodes',  IsoMesh % NumberOfNodes
       PRINT *, 'IsoMesh % NumberOfBoundaryElements', IsoMesh % NumberOfBoundaryElements
       PRINT *, 'PathCount', PathCount
    END IF

    ! boss contains isomesh and crevpaths 
    IF (Boss) THEN
       IF(IsoMesh % NumberOfNodes > 0 ) THEN ! check if there is any valid paths
          CurrentPath => CrevassePaths
          NoCrevNodes=0 ! Number of Crevasse nodes
          NoPaths=0

          DO WHILE(ASSOCIATED(CurrentPath))
             NoPaths=NoPaths+1
             NoCrevNodes=NoCrevNodes+CurrentPath % NumberOfNodes
             CurrentPath => CurrentPath % Next
          END DO

          IF(NoPaths > 0) THEN

            ALLOCATE(CrevX(NoCrevNodes),CrevY(NoCrevNodes),&
                 CrevEnd(NoPaths),CrevStart(NoPaths),CrevOrient(NoPaths,2),&
                 CrevLR(NoPaths))

            j=0
            CurrentPath => CrevassePaths
            DO WHILE(ASSOCIATED(CurrentPath))
              j=j+1
              CrevOrient(j,:) = CurrentPath % Orientation
              CrevLR(j) = CurrentPath % LeftToRight
              CurrentPath => CurrentPath % Next
            END DO

            CurrentPath => CrevassePaths
            j=1;k=1
            n=CurrentPath % ID ! first ID may not be 1..
            CrevStart(1)=1
            DO WHILE(ASSOCIATED(CurrentPath))
              DO i=1,CurrentPath % NumberOfNodes
                CrevX(j)=IsoMesh % Nodes % x(CurrentPath % NodeNumbers(i))
                CrevY(j)=IsoMesh % Nodes % y(CurrentPath % NodeNumbers(i))
                IF(n < CurrentPath % ID .AND. CurrentPath % ID>0) THEN ! non valid paths set to 0?
                  CrevEnd(k)=j-1
                  IF(k < NoPaths) CrevStart(k+1)=j
                  n=CurrentPath % ID
                  k=k+1
                END IF
                j=j+1
              END DO
              CurrentPath => CurrentPath % Next
            END DO
            IF(j/=NoCrevNodes+1) PRINT *, 'programming error'
            CrevEnd(NoPaths)=NoCrevNodes
            IF (Debug) THEN
              PRINT *, 'number of crevasse nodes', NoCrevNodes
              PRINT *, 'crevasse start numbers', CrevStart
              PRINT *, 'crevasse end numbers',CrevEnd
            END IF
          END IF
       END IF

       NodeHolder(3)=0.0_dp
       DO i=1,NoPaths
          NodeHolder(1) = CrevX(CrevStart(i))
          NodeHolder(2) = CrevY(CrevStart(i))
          NodeHolder = MATMUL(RotationMatrix, NodeHolder)
          y_coord(1) = NodeHolder(2)
          NodeHolder(1) = CrevX(CrevEnd(i))
          NodeHolder(2) = CrevY(CrevEnd(i))
          ! NodeHolder(3) = FrontNodes % z(FrontLineCount)
          !TODO - Ask Eef about this ---^
          NodeHolder = MATMUL(RotationMatrix, NodeHolder)
          y_coord(2) = NodeHolder(2)
          LeftToRight = y_coord(2) > y_coord(1) ! TODO check if this doesn't break for special cases
          IF(LeftToRight) THEN
             CrevX(CrevStart(i):CrevEnd(i))=CrevX(CrevEnd(i):CrevStart(i):-1)
             CrevY(CrevStart(i):CrevEnd(i))=CrevY(CrevEnd(i):CrevStart(i):-1)
          END IF
       END DO
       IF (Debug) Print *, CrevX
     END IF


     CALL MPI_BCAST(NoCrevNodes,1,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(NoPaths,1,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     IF (.NOT. Boss) ALLOCATE(CrevX(NoCrevNodes),CrevY(NoCrevNodes),&
          CrevEnd(NoPaths),CrevStart(NoPaths), CrevOrient(NoPaths,2),&
          CrevLR(NoPaths))! (because already created on boss)
     CALL MPI_BCAST(CrevX,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevY,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevEnd,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevStart,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevOrient,NoPaths*2,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevLR,NoPaths,MPI_LOGICAL, 0, ELMER_COMM_WORLD, ierr)

     ! if there are no crevasses no need to calculate signed distance
     IF(NoCrevNodes == 0) THEN
        CalvingOccurs = .FALSE.

        CALL SParIterAllReduceOR(CalvingOccurs)
        CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )

        !EqName = ListGetString( Params, "Remesh Equation Name", Found, UnfoundFatal = .TRUE.)
        !DO j=1,Model % NumberOfSolvers
        !  IF(ListGetString(Model % Solvers(j) % Values, "Equation") == EqName) THEN
        !    Found = .TRUE.
            !Turn off (or on) the solver
            !If CalvingOccurs, (switch) off = .true.
        !    CALL SwitchSolverExec(Model % Solvers(j), .NOT. CalvingOccurs)
        !    IF(.NOT. CalvingOccurs) CALL ResetMeshUpdate(Model, Model % Solvers(j))

            ! if remeshing if skipped need to ensure calving solvers are not paused
        !    IF(.NOT. CalvingOccurs) CALL PauseCalvingSolvers(Model, Model % Solvers(j) % Values, .FALSE.)

        !    EXIT
        !  END IF
        !END DO
        CalvingValues(CalvingPerm) = DistValues(DistPerm) + 1
        CALL WARN(SolverName, 'No crevasses so not calculating signed distance')
        !RETURN
     ELSE
        !above IF(Parallel) but that's assumed implicitly
        !make sure that code works for empty isomesh as well!!

        ! get calving polygons
        IF(Boss) THEN
          CALL GetCalvingPolygons(IsoMesh, CrevassePaths, EdgeX, EdgeY, Polygon, PolyStart, PolyEnd, gridmesh_dx)
          IF(Debug) THEN
            PRINT*, 'Calving Polygons ----'
            DO i=1, SIZE(PolyStart)
              PRINT*, PolyStart(i), ',', PolyEnd(i), ','
            END DO
            DO i=1, SIZE(Polygon(1,:))
              PRINT*, Polygon(1,i), ',', Polygon(2,i), ','
            END DO
          END IF
          ! release crevasse paths
          IF(ASSOCIATED(CrevassePaths)) CALL ReleaseCrevassePaths(CrevassePaths)
        END IF

        ! send bdrynode info to all procs
        CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
        IF(Boss) counter=SIZE(Polygon(1,:))
        CALL MPI_BCAST(counter, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
        IF(.NOT. Boss) ALLOCATE(Polygon(2, counter), PolyStart(NoPaths), &
                              PolyEnd(NoPaths))
        CALL MPI_BCAST(Polygon, Counter*2, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
        CALL MPI_BCAST(PolyEnd, NoPaths, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
        CALL MPI_BCAST(PolyStart, NoPaths, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

        CALL CheckLateralCalving(Mesh, Params, FrontPerm, CrevX,CrevY,CrevStart,CrevEnd, CrevOrient,&
          CrevLR, Polygon, PolyStart, PolyEnd)
        NoPaths = SIZE(CrevStart)

        ALLOCATE(IsCalvingNode(Mesh % NumberOfNodes))
        IsCalvingNode=.FALSE.

        IF (NoPaths > 0 ) THEN
          DO i=1,Solver % Mesh % NumberOfNodes
            IF (Debug) PRINT *, ParEnv % MyPE, 'For node i', i,' out of ',Solver % Mesh % NumberOfNodes
            xx = Solver % Mesh % Nodes % x(i)
            yy = Solver % Mesh % Nodes % y(i)
            MinDist = HUGE(1.0_dp)

            inside=.FALSE.
            ClosestCrev=0
            DO j=1, NoPaths
              ALLOCATE(PathPoly(2, PolyEnd(j)-PolyStart(j)+1))
              PathPoly = Polygon(:, PolyStart(j):PolyEnd(j))
              inside = PointInPolygon2D(PathPoly, (/xx, yy/))
              DEALLOCATE(PathPoly)
              IF(inside) THEN
                ClosestCrev = j
                EXIT
              END IF
            END DO

            ! TO DO; brute force here, checking all crevasse segments, better to find closest crev first
            DO j=1, NoPaths
              IF(inside .AND. j/=ClosestCrev) CYCLE
              DO k=CrevStart(j), CrevEnd(j)-1
                xl=CrevX(k);yl=CrevY(k)
                xr=CrevX(k+1);yr=CrevY(k+1)
                TempDist=PointLineSegmDist2D( (/xl,yl/),(/xr,yr/), (/xx,yy/))
                IF(TempDist <= (ABS(MinDist)+AEPS) ) THEN ! as in ComputeDistanceWithDirection
                  ! updated so no longer based off an angle. This caused problems when the front was
                  ! no longer projectable. Instead the node is marked to calve if it is inside the calving polygon.
                  ! PointInPolygon2D based of the winding number algorithm
                  IF(j==ClosestCrev) THEN ! inside calved area
                    IsCalvingNode(i) = .TRUE.
                    MinDist = -TempDist
                  ELSE
                    MinDist = TempDist
                  END IF
                  jmin=k ! TO DO rm jmin? not necessary anymore
                END IF
              END DO
            END DO
            SignDistValues(SignDistPerm(i)) =  MinDist
            IF(Debug)     PRINT *, ParEnv % MyPE, i, 'Shortest distance to closest segment in crevassepath ',MinDist
            IF(Debug)     PRINT *, ParEnv % MyPE, i, 'jmin is ',jmin
          END DO
        ELSE
            SignDistValues = 0.0_dp
        END IF

        Debug = .FALSE.
        CalvingValues(CalvingPerm) = SignDistValues(SignDistPerm)
        IF(MINVAL(SignDistValues) < - AEPS) CalvingOccurs = .TRUE.


        ! TO DO check in each partition whether calving sizes are sufficient and whether calving occurs by
        !  IF( integrate negative part of signdistance < -MinCalvingSize) CalvingOccurs = .TRUE.
        ! IsCalvingNode = (SignDistValues < 0 ) if they have same perm, but not actually used?    
        ! then parallel reduce with MPI_LOR on 'CalvingOccurs'
        ! NOTE; what happens if calving event is divided over 2 partitions with insignificant in each
        ! but significant in total? should something be done to avoid that they're filtered out?
        !Pass CalvingOccurs to all processes, add to simulation
        CALL SParIterAllReduceOR(CalvingOccurs)
        CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )
        
        IF(CalvingOccurs) THEN
          CALL Info( SolverName, 'Calving Event',Level=1)
        ELSE
          !If only insignificant calving events occur, reset everything
          !CalvingValues = 0.0_dp
          IsCalvingNode = .FALSE.
        END IF

        !if calving doesn't occur then no need to run remeshing solver
        !EqName = ListGetString( Params, "Remesh Equation Name", Found, UnfoundFatal = .TRUE.)
        !DO j=1,Model % NumberOfSolvers
        !  IF(ListGetString(Model % Solvers(j) % Values, "Equation") == EqName) THEN
        !    Found = .TRUE.
            !Turn off (or on) the solver
            !If CalvingOccurs, (switch) off = .true.
        !    CALL SwitchSolverExec(Model % Solvers(j), .NOT. CalvingOccurs)
        !    IF(.NOT. CalvingOccurs) CALL ResetMeshUpdate(Model, Model % Solvers(j))

            ! if remeshing if skipped need to ensure calving solvers are not paused
        !    IF(.NOT. CalvingOccurs) CALL PauseCalvingSolvers(Model, Model % Solvers(j) % Values, .FALSE.)

        !    EXIT
        !  END IF
        !END DO

        !IF(.NOT. Found) THEN
        !  WRITE (Message,'(A,A,A)') "Failed to find Equation Name: ",EqName,&
        !      " to switch off after calving."
        !  CALL Fatal(SolverName,Message)
        !END IF
    END IF

    ! because isomesh has no bulk?
    IsoMesh % NumberOfBulkElements =   IsoMesh % NumberOfBoundaryElements
    IsoMesh % NumberOfBoundaryElements = 0


    !=========================================================
    !  DONE, just plane VTU and deallocations left
    !=========================================================

    !-----------------------------------------------------------
    ! Call ResultOutputSolver in faux-serial to output calving
    ! plane for debugging
    !-----------------------------------------------------------
    IF(Boss) THEN
       VTUSolverName = "calvingresultoutput"
       DO i=1,Model % NumberOfSolvers
          IF(GetString(Model % Solvers(i) % Values, 'Equation') == VTUSolverName) THEN
             VTUOutputSolver => Model % Solvers(i)
             EXIT
          END IF
       END DO
       IF(.NOT. ASSOCIATED(VTUOutputSolver)) &
            CALL Fatal("Calving Remesh","Couldn't find VTUSolver")

       PEs = ParEnv % PEs
       ParEnv % PEs = 1

       WorkMesh => Model % Mesh
       WorkMesh2 => Model % Meshes

       Model % Mesh => PlaneMesh
       Model % Meshes => PlaneMesh
       VTUOutputSolver % Mesh => PlaneMesh
       PlaneMesh % Next => IsoMesh

       VTUOutputSolver % dt = dt
       VTUOutputSolver % NumberOfActiveElements = 0
       Model % Solver => VTUOutputSolver

       CALL SingleSolver( Model, VTUOutputSolver, VTUOutputSolver % dt, TransientSimulation )

       Model % Solver => Solver
       Model % Mesh => WorkMesh
       Model % Meshes => WorkMesh2
       PlaneMesh % Next => NULL()
       VTUOutputSolver % Mesh => WorkMesh

       ParEnv % PEs = PEs

     END IF

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken up to finish up: ', rt
    rt0 = RealTime()

    !-----------------------------------------------------------
    ! Tidy up and deallocate
    !-----------------------------------------------------------

    FirstTime = .FALSE.

    PCSolver % Variable => NULL()
    PCSolver % Matrix % Perm => NULL()
    CALL FreeMatrix(PCSolver % Matrix)
    CALL ReleaseMesh(PlaneMesh)
    ! isomesh doesn't seem to released as part of planemesh
    CALL ReleaseMesh(IsoMesh)

    DEALLOCATE(TopPerm, BotPerm, LeftPerm, RightPerm, FrontPerm, InflowPerm)

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

CONTAINS

  !Returns extent: min_x, max_x, min_y, max_y of front in rotated coord system
  FUNCTION GetFrontExtent(Mesh, RotationMatrix, DistVar, SearchDist, Buffer) RESULT(Extent)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: DistVar
    REAL(KIND=dp) :: SearchDist, RotationMatrix(3,3), Extent(4), Buffer
#ifdef ELMER_BROKEN_MPI_IN_PLACE
    REAL(KIND=dp) :: buffer2
#endif
    !------------------------------
    REAL(KIND=dp) :: NodeHolder(3)
    INTEGER :: i,ierr

    !Get the min,max x and y in rotated coords:
    extent(1) = HUGE(extent(1))
    extent(2) = -HUGE(extent(1))
    extent(3) = HUGE(extent(1))
    extent(4) = -HUGE(extent(1))

    DO i=1,Mesh % NumberOfNodes

      IF(DistVar % Values(DistVar % Perm(i)) > SearchDist) CYCLE

      NodeHolder(1) = Mesh % Nodes % x(i)
      NodeHolder(2) = Mesh % Nodes % y(i)
      NodeHolder(3) = Mesh % Nodes % z(i)
      NodeHolder = MATMUL(RotationMatrix, NodeHolder)

      extent(1) = MIN(extent(1), NodeHolder(2))
      extent(2) = MAX(extent(2), NodeHolder(2))
      extent(3) = MIN(extent(3), NodeHolder(3))
      extent(4) = MAX(extent(4), NodeHolder(3))

    END DO

#ifdef ELMER_BROKEN_MPI_IN_PLACE
    buffer2 = extent(1)
    CALL MPI_AllReduce(buffer2, &
#else
    CALL MPI_AllReduce(MPI_IN_PLACE, &
#endif
         extent(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

#ifdef ELMER_BROKEN_MPI_IN_PLACE
    buffer2 = extent(2)
    CALL MPI_AllReduce(buffer2, &
#else
    CALL MPI_AllReduce(MPI_IN_PLACE, &
#endif
         extent(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)

#ifdef ELMER_BROKEN_MPI_IN_PLACE
    buffer2 = extent(3)
    CALL MPI_AllReduce(buffer2, &
#else
    CALL MPI_AllReduce(MPI_IN_PLACE, &
#endif
         extent(3), 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)

#ifdef ELMER_BROKEN_MPI_IN_PLACE
    buffer2 = extent(4)
    CALL MPI_AllReduce(buffer2, &
#else
    CALL MPI_AllReduce(MPI_IN_PLACE, &
#endif
         extent(4), 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)

    extent(1) = extent(1) - buffer
    extent(2) = extent(2) + buffer
    extent(3) = extent(3) - buffer
    extent(4) = extent(4) + buffer
  END FUNCTION GetFrontExtent

  !subroutine to check if lateral calving would be evacuated or would jam on fjord walls
  SUBROUTINE CheckLateralCalving(Mesh, SolverParams, FrontPerm, CrevX,CrevY,CrevStart,CrevEnd, CrevOrient,&
              CrevLR, Polygon, PolyStart, PolyEnd)

    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: SolverParams
    INTEGER, POINTER :: FrontPerm(:)
    REAL(KIND=dp),ALLOCATABLE :: CrevX(:),CrevY(:),CrevOrient(:,:),Polygon(:,:)
    INTEGER, ALLOCATABLE :: CrevStart(:),CrevEnd(:),PolyStart(:),PolyEnd(:)
    LOGICAL, ALLOCATABLE :: CrevLR(:)
    !--------------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: AdvSolver
    TYPE(Valuelist_t), POINTER :: AdvParams
    TYPE(Variable_t), POINTER :: SignDistVar
    INTEGER, POINTER :: SignDistPerm(:)
    REAL(KIND=dp), POINTER :: SignDistValues(:)
    INTEGER :: i,j,NoPaths,ClosestCrev,Naux,Nl,Nr,ok,RemovePoints
    REAL(KIND=dp) :: xx,yy,MinDist,a1(2),a2(2),b1(2),b2(2),Orientation(2),intersect(2), crevdist,&
      TempDist, SecDist
    INTEGER, ALLOCATABLE :: NodeClosestCrev(:), WorkInt(:)
    REAL(KIND=dp), ALLOCATABLE :: PathPoly(:,:),xL(:),yL(:),xR(:),yR(:),WorkReal(:),WorkReal2(:,:)
    LOGICAL :: inside,does_intersect,FoundIntersect
    LOGICAL, ALLOCATABLE :: RemoveCrev(:), WorkLogical(:)
#ifdef ELMER_BROKEN_MPI_IN_PLACE
    LOGICAL, ALLOCATABLE :: buffer(:)
#endif
    CHARACTER(MAX_NAME_LEN) :: FuncName="CheckLateralCalving", Adv_EqName, LeftRailFName, RightRailFName
    INTEGER, PARAMETER :: io=20

    SignDistVar => VariableGet(Mesh % Variables, "SignedDistance", .TRUE.)
    SignDistValues => SignDistVar % Values
    SignDistPerm => SignDistVar % Perm

    NoPaths = SIZE(CrevStart)

    ALLOCATE(NodeClosestCrev(Mesh % NumberOfNodes))

    DO i=1, Mesh % NumberOfNodes

      IF(FrontPerm(i) == 0) CYCLE

      xx = Solver % Mesh % Nodes % x(i)
      yy = Solver % Mesh % Nodes % y(i)
      MinDist = HUGE(1.0_dp)

      inside=.FALSE.
      ClosestCrev=0
      DO j=1, NoPaths
        ALLOCATE(PathPoly(2, PolyEnd(j)-PolyStart(j)+1))
        PathPoly = Polygon(:, PolyStart(j):PolyEnd(j))
        inside = PointInPolygon2D(PathPoly, (/xx, yy/))
        DEALLOCATE(PathPoly)
        IF(inside) THEN
          ClosestCrev = j
          NodeClosestCrev(i) = j
          EXIT
        END IF
      END DO

      ! TO DO; brute force here, checking all crevasse segments, better to find closest crev first
      DO j=1, NoPaths
        IF(inside .AND. j/=ClosestCrev) CYCLE
        DO k=CrevStart(j), CrevEnd(j)-1
          a1(1)=CrevX(k);a1(2)=CrevY(k)
          a2(1)=CrevX(k+1);a2(2)=CrevY(k+1)
          TempDist=PointLineSegmDist2D( a1,a2, (/xx,yy/))
          IF(TempDist <= (ABS(MinDist)+AEPS) ) THEN ! as in ComputeDistanceWithDirection
            ! updated so no longer based off an angle. This caused problems when the front was
            ! no longer projectable. Instead the node is marked to calve if it is inside the calving polygon.
            ! PointInPolygon2D based of the winding number algorithm
            IF(j==ClosestCrev) THEN ! inside calved area
              MinDist = -TempDist
            ELSE
              MinDist = TempDist
            END IF
            jmin=k ! TO DO rm jmin? not necessary anymore
          END IF
        END DO
      END DO

      !hold in dist values temporarily
      SignDistValues(SignDistPerm(i)) =  MinDist
    END DO

    Adv_EqName = ListGetString(SolverParams,"Front Advance Solver", UnfoundFatal=.TRUE.)
    ! Locate CalvingAdvance Solver
    DO i=1,Model % NumberOfSolvers
      IF(GetString(Model % Solvers(i) % Values, 'Equation') == Adv_EqName) THEN
         AdvSolver => Model % Solvers(i)
         EXIT
      END IF
    END DO
    AdvParams => AdvSolver % Values

    LeftRailFName = ListGetString(AdvParams, "Left Rail File Name", Found)
    IF(.NOT. Found) THEN
      CALL Info(FuncName, "Left Rail File Name not found, assuming './LeftRail.xy'")
      LeftRailFName = "LeftRail.xy"
    END IF
    Nl = ListGetInteger(AdvParams, "Left Rail Number Nodes", Found)
    IF(.NOT.Found) THEN
      WRITE(Message,'(A,A)') 'Left Rail Number Nodes not found'
      CALL FATAL(FuncName, Message)
    END IF
    !TO DO only do these things if firsttime=true?
    OPEN(unit = io, file = TRIM(LeftRailFName), status = 'old',iostat = ok)
    IF (ok /= 0) THEN
      WRITE(message,'(A,A)') 'Unable to open file ',TRIM(LeftRailFName)
      CALL FATAL(Trim(FuncName),Trim(message))
    END IF
    ALLOCATE(xL(Nl), yL(Nl))

    ! read data
    DO i = 1, Nl
      READ(io,*,iostat = ok, end=200) xL(i), yL(i)
    END DO
200   Naux = Nl - i
    IF (Naux > 0) THEN
      WRITE(Message,'(I0,A,I0,A,A)') Naux,' out of ',Nl,' datasets in file ', TRIM(LeftRailFName)
      CALL INFO(Trim(FuncName),Trim(message))
    END IF
    CLOSE(io)
    RightRailFName = ListGetString(AdvParams, "Right Rail File Name", Found)
    IF(.NOT. Found) THEN
      CALL Info(FuncName, "Right Rail File Name not found, assuming './RightRail.xy'")
      RightRailFName = "RightRail.xy"
    END IF

    Nr = ListGetInteger(AdvParams, "Right Rail Number Nodes", Found)
    IF(.NOT.Found) THEN
      WRITE(Message,'(A,A)') 'Right Rail Number Nodes not found'
      CALL FATAL(FuncName, Message)
    END IF
    !TO DO only do these things if firsttime=true?
    OPEN(unit = io, file = TRIM(RightRailFName), status = 'old',iostat = ok)

    IF (ok /= 0) THEN
      WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(RightRailFName)
      CALL FATAL(Trim(FuncName),Trim(message))
    END IF
    ALLOCATE(xR(Nr), yR(Nr))

    ! read data
    DO i = 1, Nr
       READ(io,*,iostat = ok, end=100) xR(i), yR(i)
    END DO
100   Naux = Nr - i
    IF (Naux > 0) THEN
      WRITE(message,'(I0,A,I0,A,A)') Naux,' out of ',Nl,' datasets in file ', TRIM(RightRailFName)
      CALL INFO(Trim(FuncName),Trim(message))
    END IF
    CLOSE(io)

    ALLOCATE(RemoveCrev(NoPaths))
#ifdef ELMER_BROKEN_MPI_IN_PLACE
    ALLOCATE(buffer(NoPaths))
#endif
    RemoveCrev = .FALSE.
    DO i=1, Mesh % NumberOfNodes

      IF(FrontPerm(i) == 0) CYCLE
      IF(SignDistValues(SignDistPerm(i)) >= 0) CYCLE
      IF(RemoveCrev(NodeClosestCrev(i))) CYCLE

      b1(1) = Solver % Mesh % Nodes % x(i)
      b1(2) = Solver % Mesh % Nodes % y(i)
      MinDist = HUGE(1.0_dp)

      Orientation = CrevOrient(NodeClosestCrev(i),:)

      IF(CrevLR(NodeClosestCrev(i))) THEN
        b2 = b1 - 10*Orientation
      ELSE
        b2 = b1 + 10*Orientation
      END IF

      FoundIntersect = .FALSE.
      DO j=1, Nl-1
        a1 = (/xL(j), yL(j)/)
        a2 = (/xL(j+1), yL(j+1)/)
        IF(PointInPolygon2D(Polygon(:, PolyStart(NodeClosestCrev(i)):PolyEnd(NodeClosestCrev(i))),a1)) THEN
          DO k=CrevStart(NodeClosestCrev(i)), CrevEnd(NodeClosestCrev(i))-1
            CALL LineSegmentsIntersect (a1, a2, (/CrevX(k), CrevY(k)/), (/CrevX(k+1), CrevY(k+1)/),&
                           intersect, does_intersect )
            IF(does_intersect) a1 = intersect
          END DO
        END IF
        IF(PointInPolygon2D(Polygon(:, PolyStart(NodeClosestCrev(i)):PolyEnd(NodeClosestCrev(i))),a2)) THEN
          DO k=CrevStart(NodeClosestCrev(i)), CrevEnd(NodeClosestCrev(i))-1
            CALL LineSegmentsIntersect (a1, a2, (/CrevX(k), CrevY(k)/), (/CrevX(k+1), CrevY(k+1)/),&
                           intersect, does_intersect )
            IF(does_intersect) a2 = intersect
          END DO
        END IF
        IF(a1(1)==a2(1) .AND. a1(2)==a2(2)) CYCLE
        CALL LineSegmLineIntersect (a1, a2, b1, b2, intersect, does_intersect )

        IF(.NOT. does_intersect) CYCLE
        tempdist = PointDist2D(b1,intersect)
        secdist = PointDist2D(b2, intersect)
        IF(secdist > tempdist) CYCLE
        IF(tempdist < mindist) THEN
          mindist = tempdist
          FoundIntersect = .TRUE.
        END IF
      END DO

      IF(.NOT. FoundIntersect) THEN
        DO j=1, Nr-1
          a1 = (/xR(j), yR(j)/)
          a2 = (/xR(j+1), yR(j+1)/)
          IF(PointInPolygon2D(Polygon(:, PolyStart(NodeClosestCrev(i)):PolyEnd(NodeClosestCrev(i))),a1)) THEN
            DO k=CrevStart(NodeClosestCrev(i)), CrevEnd(NodeClosestCrev(i))-1
              CALL LineSegmentsIntersect (a1, a2, (/CrevX(k), CrevY(k)/), (/CrevX(k+1), CrevY(k+1)/),&
                             intersect, does_intersect )
              IF(does_intersect) a1 = intersect
            END DO
          END IF
          IF(PointInPolygon2D(Polygon(:, PolyStart(NodeClosestCrev(i)):PolyEnd(NodeClosestCrev(i))),a2)) THEN
            DO k=CrevStart(NodeClosestCrev(i)), CrevEnd(NodeClosestCrev(i))-1
              CALL LineSegmentsIntersect (a1, a2, (/CrevX(k), CrevY(k)/), (/CrevX(k+1), CrevY(k+1)/),&
                             intersect, does_intersect )
              IF(does_intersect) a2 = intersect
            END DO
          END IF
          IF(a1(1)==a2(1) .AND. a1(2)==a2(2)) CYCLE
          CALL LineSegmLineIntersect (a1, a2, b1, b2, intersect, does_intersect )
          IF(.NOT. does_intersect) CYCLE
          tempdist = PointDist2D(b1,intersect)
          secdist = PointDist2D(b2, intersect)
          IF(secdist > tempdist) CYCLE
          IF(tempdist < mindist) THEN
            mindist = tempdist
            FoundIntersect = .TRUE.
          END IF
        END DO
      END IF

      crevdist = 0.0_dp
      IF(FoundIntersect) THEN
        FoundIntersect = .FALSE.
        crevdist = HUGE(1.0_dp)
        DO j=CrevStart(NodeClosestCrev(i)), CrevEnd(NodeClosestCrev(i))-1
          a1 = (/CrevX(j), CrevY(j)/)
          a2 = (/CrevX(j+1), CrevY(j+1)/)
          CALL LineSegmLineIntersect (a1, a2, b1, b2, intersect, does_intersect )
          IF(.NOT. does_intersect) CYCLE
          tempdist = PointDist2D(b1,intersect)
          IF(tempdist < crevdist) THEN
            crevdist = tempdist
          END IF
          FoundIntersect = .TRUE.
        END DO
      END IF

      IF(MinDist < crevdist .AND. FoundIntersect) THEN
        CALL WARN(FuncName, 'Removing lateral calving event as it would jam on lateral margins')
        RemoveCrev(NodeClosestCrev(i)) = .TRUE.
      END IF
    END DO

#ifdef ELMER_BROKEN_MPI_IN_PLACE
         buffer = RemoveCrev
         CALL MPI_ALLREDUCE(buffer, &
#else
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, &
#endif
        RemoveCrev, NoPaths, MPI_LOGICAL, MPI_LOR, ELMER_COMM_WORLD, ierr)

    DO WHILE(ANY(RemoveCrev))
      DO i=1, NoPaths
        IF(.NOT. RemoveCrev(i)) CYCLE
        !CrevX,CrevY,CrevStart,CrevEnd, CrevOrient,&
        !Polygon, PolyStart, PolyEnd)

        !change crevasse allocations
        RemovePoints = CrevEnd(i) - CrevStart(i) + 1
        ALLOCATE(WorkReal(CrevEnd(NoPaths) - RemovePoints))
        ! alter crevx
        IF(i > 1) THEN
          WorkReal(1:CrevEnd(i-1)) = CrevX(1:CrevEnd(i-1))
          IF(i < NoPaths) WorkReal(CrevEnd(i-1)+1:) = CrevX(CrevStart(i+1):)
        ELSE IF (i < NoPaths) THEN
          WorkReal(1:) = CrevX(CrevStart(i+1):)
        END IF
        DEALLOCATE(CrevX)
        ALLOCATE(CrevX((CrevEnd(NoPaths) - RemovePoints)))
        CrevX = WorkReal
        ! alter crevy
        IF(i > 1) THEN
          WorkReal(1:CrevEnd(i-1)) = CrevY(1:CrevEnd(i-1))
          IF(i < NoPaths) WorkReal(CrevEnd(i-1)+1:) = CrevY(CrevStart(i+1):)
        ELSE IF (i < NoPaths) THEN
          WorkReal(1:) = CrevY(CrevStart(i+1):)
        END IF
        DEALLOCATE(CrevY)
        ALLOCATE(CrevY((CrevEnd(NoPaths) - RemovePoints)))
        CrevY = WorkReal
        DEALLOCATE(WorkReal)

        ALLOCATE(WorkInt(NoPaths-1))
        !alter crevstart
        IF(i > 1) THEN
          WorkInt(1:i-1) = CrevStart(1:i-1)
          IF(i < NoPaths) WorkInt(i:) = CrevStart(i+1:) - RemovePoints
        ELSE IF (i < NoPaths) THEN
          WorkInt = CrevStart(i+1:) - RemovePoints
        END IF
        DEALLOCATE(CrevStart)
        ALLOCATE(CrevStart(NoPaths-1))
        CrevStart = WorkInt
        !alter crevend
        IF(i > 1) THEN
          WorkInt(1:i-1) = CrevEnd(1:i-1)
          IF(i < NoPaths) WorkInt(i:) = CrevEnd(i+1:) - RemovePoints
        ELSE IF (i < NoPaths) THEN
          WorkInt = CrevEnd(i+1:) - RemovePoints
        END IF
        DEALLOCATE(CrevEnd)
        ALLOCATE(CrevEnd(NoPaths-1))
        CrevEnd = WorkInt

        !now polygons
        RemovePoints = PolyEnd(i) - PolyStart(i) + 1
        ALLOCATE(WorkReal2(2, PolyEnd(NoPaths) - RemovePoints))
        IF(i > 1) THEN
          WorkReal2(:,1:PolyEnd(i-1)) = Polygon(:,1:PolyEnd(i-1))
          IF(i < NoPaths) WorkReal2(:,PolyEnd(i-1)+1:) = Polygon(:,PolyStart(i+1):)
        ELSE IF (i < NoPaths) THEN
          WorkReal2(:,1:) = Polygon(:,PolyStart(i+1):)
        END IF
        DEALLOCATE(Polygon)
        ALLOCATE(Polygon(2,PolyEnd(NoPaths) - RemovePoints))
        Polygon = WorkReal2
        DEALLOCATE(WorkReal2)

        !alter polystart
        IF(i > 1) THEN
          WorkInt(1:i-1) = PolyStart(1:i-1)
          IF(i < NoPaths) WorkInt(i:) = PolyStart(i+1:) - RemovePoints
        ELSE IF (i < NoPaths) THEN
          WorkInt = PolyStart(i+1:) - RemovePoints
        END IF
        DEALLOCATE(PolyStart)
        ALLOCATE(PolyStart(NoPaths-1))
        PolyStart = WorkInt
        !alter polyend
        IF(i > 1) THEN
          WorkInt(1:i-1) = PolyEnd(1:i-1)
          IF(i < NoPaths) WorkInt(i:) = PolyEnd(i+1:) - RemovePoints
        ELSE IF (i < NoPaths) THEN
          WorkInt = PolyEnd(i+1:) - RemovePoints
        END IF
        DEALLOCATE(PolyEnd)
        ALLOCATE(PolyEnd(NoPaths-1))
        PolyEnd = WorkInt
        DEALLOCATE(WorkInt)

        !finally remove crev marker
        ALLOCATE(WorkLogical(NoPaths-1))
        IF(i > 1) THEN
          WorkLogical(1:i-1) = RemoveCrev(1:i-1)
          IF(i < NoPaths) WorkLogical(i:) = RemoveCrev(i+1:)
        ELSE IF (i < NoPaths) THEN
          WorkLogical = RemoveCrev(i+1:)
        END IF
        DEALLOCATE(RemoveCrev)
        ALLOCATE(RemoveCrev(NoPaths-1))
        RemoveCrev = WorkLogical
        DEALLOCATE(WorkLogical)

        EXIT
      END DO

      NoPaths = SIZE(RemoveCrev)
    END DO

  END SUBROUTINE CheckLateralCalving

  !Cleanly removes elements & nodes from a mesh based on mask
  !Any element containing a removed node (RmNode) will be deleted
  !May optionally specify which elements to remove
  !If only RmElem is provided, no nodes are removed
  !(should this be changed? i.e. detect orphaned nodes?)
  SUBROUTINE CutPlaneMesh(Mesh, RmNode, RmElem)

    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: RmNode(:)
    LOGICAL, OPTIONAL, TARGET :: RmElem(:)
    !--------------------------------
    TYPE(Element_t), POINTER :: Element, Work_Elements(:)
    TYPE(Variable_t), POINTER :: Variables, Var
    TYPE(Nodes_t), POINTER :: Nodes
    TYPE(NeighbourList_t), POINTER :: work_neighlist(:)
    REAL(KIND=dp), ALLOCATABLE :: work_xyz(:,:)
    REAL(KIND=dp), POINTER CONTIG :: work_x(:),work_y(:), work_z(:)
    REAL(KIND=dp), POINTER :: WorkReal(:)=>NULL(), ptr(:)
    INTEGER :: i,j,counter,NNodes,NBulk, NBdry,NewNNodes, NewNElems, NewNBulk,&
         NewNbdry, ElNNodes
    INTEGER, ALLOCATABLE :: Nodeno_map(:),EIdx_map(:)
    INTEGER, POINTER :: NodeIndexes(:), work_pInt(:),WorkPerm(:)=>NULL()
    LOGICAL, POINTER :: RmElement(:),work_logical(:)
    CHARACTER(*), PARAMETER :: FuncName="CutMesh"
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    NNodes = Mesh % NumberOfNodes
    NBulk = Mesh % NumberOfBulkElements
    NBdry = Mesh % NumberOfBoundaryElements
    Variables => Mesh % Variables

    IF(.NOT. (PRESENT(RmElem) .OR. PRESENT(RmNode))) &
         CALL Fatal(FuncName,"Need to provide at least one of RmElem, RmNode!")

    IF(PRESENT(RmElem)) THEN
      RmElement => RmElem
    ELSE
      !Mark elements containing deleted nodes for deletion
      ALLOCATE(RmElement(NBulk + Nbdry))
      RmElement = .FALSE.
      DO i=1,NBulk + NBdry
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes
        ElNNodes = Element % TYPE % NumberOfNodes
        IF(ANY(RmNode(NodeIndexes(1:ElNNodes)))) RmElement(i) = .TRUE.
        !Get rid of orphans
        IF(i > NBulk) THEN
          IF(ASSOCIATED(Element % BoundaryInfo)) THEN
            !If the main body element is removed, remove this BC elem
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % ElementIndex
              IF(RmElement(j)) RmElement(i) = .TRUE.
            END IF
            !Point % Right => NULL() if % Right element is removed
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j = Element % BoundaryInfo % Right % ElementIndex
              IF(RmElement(j)) Element % BoundaryInfo % Right => NULL()
            END IF
          END IF
        END IF
      END DO
    END IF

    IF(PRESENT(RmNode)) THEN
      !Removing nodes implies shifting element nodeindexes
      !Map pre -> post deletion node nums
      ALLOCATE(Nodeno_map(NNodes))
      Nodeno_map = 0
      counter = 0
      DO i=1,NNodes
        IF(RmNode(i)) CYCLE
        counter = counter + 1
        Nodeno_map(i) = counter
      END DO

      !Update the element nodeindexes
      DO i=1,NBulk+NBdry
        Element => Mesh % Elements(i)
        IF(RmElement(i)) CYCLE
        DO j=1,Element % TYPE % NumberOfNodes
          Element % NodeIndexes(j) = Nodeno_map(Element % NodeIndexes(j))
          IF(Element % NodeIndexes(j) == 0) CALL Fatal(FuncName, &
              "Programming error: mapped nodeno = 0")
        END DO
      END DO

      !Clear out deleted nodes
      Nodes => Mesh % Nodes
      NewNNodes = COUNT(.NOT. RmNode)

      ALLOCATE(work_x(NewNNodes),&
          work_y(NewNNodes),&
          work_z(NewNNodes))

      counter = 0
      DO i=1,NNodes
        IF(RmNode(i)) CYCLE
        counter = counter + 1
        work_x(counter) = Nodes % x(i)
        work_y(counter) = Nodes % y(i)
        work_z(counter) = Nodes % z(i)
      END DO

      DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
      Nodes % x => work_x
      Nodes % y => work_y
      Nodes % z => work_z

      Nodes % NumberOfNodes = NewNNodes
      Mesh % NumberOfNodes = NewNNodes

      !update perms
      Var => Variables
      DO WHILE(ASSOCIATED(Var))
        IF(SIZE(Var % Values) == Var % DOFs) THEN
          Var => Var % Next
          CYCLE
        ELSE IF(LEN(Var % Name) > 10) THEN
          IF(Var % Name(1:10)=='coordinate') THEN
            Var => Var % Next
            CYCLE
          END IF
        END IF

        ALLOCATE(WorkPerm(NewNNodes))
        WorkPerm = [(i,i=1,NewNNodes)]
        IF(.NOT. ASSOCIATED(WorkReal)) ALLOCATE(WorkReal(NewNNodes))

        VarName = Var % Name

        counter = 0
        DO i=1,NNodes
          IF(RmNode(i)) CYCLE
          counter = counter + 1
          WorkReal(counter) = Var % Values(Var % Perm(i))
        END DO

        ptr => Var % Values(i::Var % DOFs)

        !DEALLOCATE(Var % Perm)
        IF ( ASSOCIATED( Var % Values ) ) &
           DEALLOCATE( Var % Values )

        Var % Perm => WorkPerm
        Var % Values => WorkReal
        NULLIFY(WorkPerm, WorkReal)
        Var => Var % Next
      END DO

      !Clear out ParallelInfo
      IF(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs)) THEN
        ALLOCATE(work_pInt(NewNNodes))
        counter = 0
        DO i=1,NNodes
          IF(RmNode(i)) CYCLE
          counter = counter + 1
          work_pInt(counter) = Mesh % ParallelInfo % GlobalDOFs(i)
        END DO
        DEALLOCATE(Mesh % ParallelInfo % GlobalDOFs)
        Mesh % ParallelInfo % GlobalDOFs => work_pInt
        work_pInt => NULL()
      END IF

      !Get rid of NeighbourList
      IF(ASSOCIATED(Mesh % ParallelInfo % NeighbourList)) THEN
        ALLOCATE(work_neighlist(NewNNodes))
        DO i=1,NNodes
          IF(.NOT. RmNode(i)) CYCLE
          IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
              DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
        END DO

        counter = 0
        DO i=1,NNodes
          IF(RmNode(i)) CYCLE
          counter = counter + 1
          work_neighlist(counter) % Neighbours => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
        END DO
        DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
        Mesh % ParallelInfo % NeighbourList => work_neighlist
        work_neighlist => NULL()
      END IF

      !Get rid of ParallelInfo % NodeInterface
      IF(ASSOCIATED(Mesh % ParallelInfo % GInterface)) THEN
        ALLOCATE(work_logical(NewNNodes))
        counter = 0
        DO i=1,NNodes
          IF(RmNode(i)) CYCLE
          counter = counter + 1
          work_logical(counter) = Mesh % ParallelInfo % GInterface(i)
        END DO
        DEALLOCATE(Mesh % ParallelInfo % GInterface)
        Mesh % ParallelInfo % GInterface => work_logical
        work_logical => NULL()
      END IF
    END IF

    !TODO - Mesh % Edges - see ReleaseMeshEdgeTables
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Edges not yet implemented.")
    END IF

    !TODO - Mesh % Faces - see ReleaseMeshFaceTables
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Faces not yet implemented.")
    END IF

    !TODO - Mesh % ViewFactors -  see ReleaseMeshFactorTables
    IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % ViewFactors not yet implemented.")
    END IF

    !TODO - Mesh % Projector - see FreeMatrix in ReleaseMesh
    IF (ASSOCIATED(Mesh % Projector) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % Projector not yet implemented.")
    END IF

    !TODO - Mesh % RootQuadrant - see FreeQuadrantTree
    IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
      CALL Fatal(FuncName,"Clearing out Mesh % RootQuadrant not yet implemented.")
    END IF

    NewNElems = COUNT(.NOT. RmElement)


    !Clear out deleted elements structures - taken from ReleaseMesh
    DO i=1,NBulk+NBdry
      IF(.NOT. RmElement(i)) CYCLE

      !          Boundaryinfo structure for boundary elements
      !          ---------------------------------------------
      IF ( Mesh % Elements(i) % Copy ) CYCLE

      IF ( i > NBulk ) THEN
        IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo ) ) THEN
          DEALLOCATE( Mesh % Elements(i) % BoundaryInfo )
        END IF
      END IF

      IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
      Mesh % Elements(i) % NodeIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
      Mesh % Elements(i) % EdgeIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
      Mesh % Elements(i) % FaceIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % DGIndexes )
      Mesh % Elements(i) % DGIndexes => NULL()

      IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
           DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
      Mesh % Elements(i) % BubbleIndexes => NULL()

      ! This creates problems later on!!!
      !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
      !   DEALLOCATE( Mesh % Elements(i) % PDefs )

      Mesh % Elements(i) % PDefs => NULL()
    END DO

    !Construct a map of old element indexes => new element indexes
    ALLOCATE(EIdx_map(NBdry+NBulk))
    EIdx_map = 0
    counter = 0
    DO i=1,NBulk+NBdry
      IF(RmElement(i)) CYCLE
      counter = counter + 1
      IF(Mesh % Elements(i) % ElementIndex /= i) CALL Warn(FuncName,&
           "Assumption Elements(i) % ElementIndex == i not valid! Expect memory corruption...")
      EIdx_map(Mesh % Elements(i) % ElementIndex) = counter
    END DO

    !Repoint elements
    ALLOCATE(work_elements(NewNElems))
    counter = 0
    DO i=1,Nbulk+NBdry
      IF(RmElement(i)) CYCLE
      counter = counter + 1

      Element => Mesh % Elements(i)
      work_elements(counter) = Element

      !Repoint BoundaryInfo % Left, % Right
      IF(i > NBdry) THEN
        IF(ASSOCIATED(Element % BoundaryInfo)) THEN
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            j = Element % BoundaryInfo % Left % ElementIndex
            work_elements(counter) % BoundaryInfo % Left => work_elements(EIdx_map(j))
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            j = Element % BoundaryInfo % Right % ElementIndex
            work_elements(counter) % BoundaryInfo % Right => work_elements(EIdx_map(j))
          END IF
        END IF
      END IF
    END DO

    !Update ElementIndexes
    DO i=1,NewNElems
      work_elements(i) % ElementIndex = i
    END DO

    DEALLOCATE(Mesh % Elements)
    Mesh % Elements => work_elements
    work_elements => NULL()

    Mesh % NumberOfBulkElements = COUNT(.NOT. RmElement(1:nbulk))
    Mesh % NumberOfBoundaryElements = COUNT(.NOT. RmElement(nbulk+1:nbulk+nbdry))

    IF(.NOT. PRESENT(RmElem)) DEALLOCATE(RmElement)

  END SUBROUTINE CutPlaneMesh

END SUBROUTINE Find_Calving3D_LSet
