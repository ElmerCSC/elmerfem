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
        node, nodecounter, CurrentNodePosition, StartNode, NodePositions(3), &
        directions, Counter, ClosestCrev, NumEdgeNodes, UnFoundLoops
   INTEGER, POINTER :: CalvingPerm(:), TopPerm(:)=>NULL(), BotPerm(:)=>NULL(), &
        LeftPerm(:)=>NULL(), RightPerm(:)=>NULL(), FrontPerm(:)=>NULL(), &
        FrontNodeNums(:), FaceNodeNums(:)=>NULL(), DistPerm(:), WorkPerm(:), &
        SignDistPerm(:), NodeIndexes(:),IceNodeIndexes(:),&
        EdgeMap(:,:)
   INTEGER, ALLOCATABLE :: CrevEnd(:),CrevStart(:),IMBdryConstraint(:),IMBdryENums(:),&
         PolyStart(:), PolyEnd(:), EdgeLine(:,:), EdgeCount(:), Nodes(:), StartNodes(:,:),&
         WorkInt(:), IMBdryNNums(:), WorkInt2D(:,:)
   REAL(KIND=dp) :: FrontOrientation(3), &
        RotationMatrix(3,3), UnRotationMatrix(3,3), NodeHolder(3), &
        MaxMeshDist, MeshEdgeMinLC, MeshEdgeMaxLC, MeshLCMinDist, MeshLCMaxDist,&
        CrevasseThreshold, MinCalvingSize, PauseVolumeThresh, &
        y_coord(2), TempDist,MinDist, xl,xr,yl, yr, xx,yy,&
        angle,angle0,a1(2),a2(2),b1(2),b2(2),a2a1(2),isect(2),front_extent(4), &
        buffer, gridmesh_dx, FrontLeft(2), FrontRight(2),&
#ifdef USE_ISO_C_BINDINGS
        rt0, rt
#else
        rt0, rt, RealTime
#endif

   REAL(KIND=dp), POINTER :: DistValues(:), SignDistValues(:), WorkReal(:), &
        CalvingValues(:)
   REAL(KIND=dp), ALLOCATABLE :: CrevX(:),CrevY(:),IMBdryNodes(:,:), Polygon(:,:), PathPoly(:,:), &
        EdgeX(:), EdgeY(:), EdgePoly(:,:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, DistVarname, &
        filename_root, &
        FrontMaskName,TopMaskName,BotMaskName,LeftMaskName,RightMaskName, &
        PC_EqName, Iso_EqName, VTUSolverName, NameSuffix,&
        MoveMeshDir
   LOGICAL :: Found, Parallel, Boss, Debug, FirstTime = .TRUE., CalvingOccurs=.FALSE., &
        SaveParallelActive, PauseSolvers, LeftToRight, MoveMesh=.FALSE., inside, Complete
   LOGICAL, ALLOCATABLE :: RemoveNode(:), IMOnFront(:), IMOnSide(:), IMOnMargin(:), &
        IMOnLeft(:), IMOnRight(:), FoundNode(:), &
        IMElemOnMargin(:), DeleteMe(:), IsCalvingNode(:), PlaneEdgeElem(:), EdgeNode(:), UsedElem(:)

   TYPE(CrevassePath_t), POINTER :: CrevassePaths, CurrentPath

   SAVE :: FirstTime, SolverName, Params, Parallel, Boss, dim, Debug, &
        DistVarName, PC_EqName, Iso_EqName, &
        MinCalvingSize, PauseVolumeThresh, MoveMesh,LeftConstraint, &
        RightConstraint, FrontConstraint,TopMaskName, BotMaskName, &
        LeftMaskName, RightMaskName, FrontMaskName

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

      dim = CoordinateSystemDimension()
      IF(dim /= 3) CALL Fatal(SolverName, "Solver only works in 3D!")

      DistVarName = ListGetString(Params,"Distance Variable Name", Found)
      IF(.NOT. Found) DistVarName = "Distance"

      PC_EqName = ListGetString(Params,"Project Calving Equation Name",Found, UnfoundFatal=.TRUE.)
      Iso_EqName = ListGetString(Params,"Isosurface Equation Name",Found, UnfoundFatal=.TRUE.)

      MinCalvingSize = ListGetConstReal(Params, "Minimum Calving Event Size", Found)
      IF(.NOT. Found) THEN
         MinCalvingSize = 0.0_dp
         CALL Warn(SolverName,"No 'Minimum Calving Event Size' given, simulation may run slowly&
              & due to remeshing after tiny events!")
      END IF

      PauseVolumeThresh = ListGetConstReal(Params, "Pause Solvers Minimum Iceberg Volume", Found)
      IF(.NOT. Found) PauseVolumeThresh = 0.0_dp

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

   Mesh => Model % Mesh

   !TODO - these 4 are defunct
   MeshEdgeMinLC = ListGetConstReal(Params, "Calving Mesh Min LC",Found, UnfoundFatal=.TRUE.)
   MeshEdgeMaxLC = ListGetConstReal(Params, "Calving Mesh Max LC",Found, UnfoundFatal=.TRUE.)
   MeshLCMinDist = ListGetConstReal(Params, "Calving Mesh LC Min Dist",Found, UnfoundFatal=.TRUE.)
   MeshLCMaxDist = ListGetConstReal(Params, "Calving Mesh LC Max Dist",Found, UnfoundFatal=.TRUE.)


   MaxMeshDist = ListGetConstReal(Params, "Calving Search Distance",Found, UnfoundFatal=.TRUE.)
   CrevasseThreshold = ListGetConstReal(Params, "Crevasse Penetration Threshold", Found, &
        UnfoundFatal=.TRUE.)

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
        RightPerm(NoNodes), FrontPerm(NoNodes))

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

   !Get the orientation of the calving front, compute rotation matrix
   FrontOrientation = GetFrontOrientation(Model)
   RotationMatrix = ComputeRotationMatrix(FrontOrientation)
   UnRotationMatrix = TRANSPOSE(RotationMatrix)

   NameSuffix = ListGetString(Params, "Calving Append Name", Found, UnfoundFatal=.TRUE.)
   MoveMeshDir = ListGetString(Params, "Calving Move Mesh Dir", Found)
   IF(Found) THEN
     CALL Info(SolverName, "Moving temporary mesh files after done")
     MoveMesh = .TRUE.
   END IF

   WRITE(filename_root,'(A,A)') "Calving_temp",TRIM(NameSuffix)

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

   gridmesh_dx = 20.0_dp
   buffer = gridmesh_dx * 4

   front_extent = GetFrontExtent(Mesh, RotationMatrix, DistVar, MaxMeshDist, buffer)

   CALL ListAddConstReal(MeshParams,"Grid Mesh Min X",front_extent(1))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Max X",front_extent(2))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Min Y",front_extent(3))
   CALL ListAddConstReal(MeshParams,"Grid Mesh Max Y",front_extent(4))
   CALL ListAddConstReal(MeshParams,"Grid Mesh dx",gridmesh_dx)


   PRINT *,ParEnv % MyPE,' front extent: ',front_extent

   PlaneMesh => CreateRectangularMesh(MeshParams)
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

      ALLOCATE(UsedElem(PlaneMesh % NumberOfBulkElements), StartNodes(2, 5))
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

      ! translate edge nodes into x and y and join both strings
      ALLOCATE(EdgeX(SUM(EdgeCount)), EdgeY(SUM(EdgeCount)))
      nodecounter=0
      IF(directions > 1) THEN
        DO i=EdgeCount(2),1, -1 ! backwards iteration
          nodecounter=nodecounter+1
          EdgeX(nodecounter) = PlaneMesh % Nodes % x(EdgeLine(2, i))
          EdgeY(nodecounter) = PlaneMesh % Nodes % y(EdgeLine(2, i))
        END DO
      END IF
      DO i=1, EdgeCount(1)
        nodecounter=nodecounter+1
        EdgeX(nodecounter) = PlaneMesh % Nodes % x(EdgeLine(1, i))
        EdgeY(nodecounter) = PlaneMesh % Nodes % y(EdgeLine(1, i))
      END DO

      ! end of finding planmesh edge

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

       ALLOCATE(IMOnFront(IsoMesh % NumberOfNodes), &
            IMOnLeft(IsoMesh % NumberOfNodes),&
            IMOnRight(IsoMesh % NumberOfNodes),&
            IMOnSide(IsoMesh % NumberOfNodes),&
            IMOnMargin(IsoMesh % NumberOfNodes))
       IMOnFront=.FALSE.; IMOnSide=.FALSE.; IMOnMargin=.FALSE.
       IMOnLeft=.FALSE.; IMOnRight=.FALSE.

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
             WorkInt(1+(IMBdryCount-1)*2: IMBdryCount*2) = NodeIndexes ! becomes IMBdryNNums
           END IF
         END DO
         IF(k==1) ALLOCATE(IMBdryNodes(IMBdryCount,4),IMBdryENums(IMBdryCount), WorkInt(IMBdryCount*2))
       END DO

     END IF

     IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

     !Pass isoline boundary elements to all partitions to get info
     !about which boundary they cross

     CALL MPI_BCAST(IMBdryCount,1,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     IF(.NOT. Boss) ALLOCATE(IMBdryNodes(IMBdryCount,4))

     CALL MPI_BCAST(IMBdryNodes,IMBdryCount*4, MPI_DOUBLE_PRECISION,&
          0, ELMER_COMM_WORLD, ierr)
     ALLOCATE(IMBdryConstraint(IMBdryCount))
     IMBdryConstraint = 0

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
         IF(Found) EXIT
       END DO
     END DO

     !Send info back to boss
     IF(Boss) THEN
       CALL MPI_Reduce(MPI_IN_PLACE, IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
            MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
     ELSE
       CALL MPI_Reduce(IMBdryConstraint, IMBdryConstraint, IMBdryCount, MPI_INTEGER, &
            MPI_MAX, 0, ELMER_COMM_WORLD, ierr)
     END IF

     CALL GetFrontCorners(Model, FrontLeft, FrontRight)

     IF(Boss) THEN

       IF(ANY(IMBdryConstraint == 0)) THEN
         PRINT *,'Debug constraints: ', IMBdryConstraint
         CALL Fatal(SolverName,"Failed to identify boundary constraint for all isomesh edge elements")
       END IF

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

       ALLOCATE(IMBdryNNums(IMBdryCount))
       nodecounter=0
       IMBdryNNums=0
       DO i=1, IMBdryCount*2
        counter=0
        DO j=1, IsoMesh % NumberOfBulkElements
          NodeIndexes => Isomesh % Elements(j) % NodeIndexes
          IF(.NOT. ANY(NodeIndexes == WorkInt(i))) CYCLE
          counter=counter+1
        END DO
        IF(counter==1) THEN
          Nodecounter=nodecounter+1
          IMBdryNNums(nodecounter) = WorkInt(i)
        END IF
       END DO
       DEALLOCATE(WorkInt)

       DO i=1,IMBdryCount
         k = IMBdryNNums(i)
         PRINT*, i, k, IMBdryConstraint(i)
         IF(k==0) CYCLE
         IF(IMBdryConstraint(i) == FrontConstraint) IMOnFront(k) = .TRUE.
         IF(IMBdryConstraint(i) == LeftConstraint) IMOnLeft(k) = .TRUE.
         IF(IMBdryConstraint(i) == RightConstraint) IMOnRight(k) = .TRUE.
       END DO

       IF(Debug) THEN
          PRINT *, 'debug, count IMOnFront: ', COUNT(IMOnFront)
          PRINT *, 'debug, count IMOnSide: ', COUNT(IMOnSide)
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
       CALL CheckCrevasseNodes(IsoMesh, CrevassePaths, IMOnLeft, IMOnRight)
       !CALL ValidateCrevassePaths(IsoMesh, CrevassePaths, FrontOrientation, PathCount,&
       !     IMOnLeft, IMOnRight, .FALSE.)
       CALL RemoveInvalidCrevs(IsoMesh, CrevassePaths, EdgeX, EdgeY)
       CALL ValidateNPCrevassePaths(IsoMesh, CrevassePaths, IMOnLeft, IMOnRight, FrontLeft, FrontRight)

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
                 CrevEnd(NoPaths),CrevStart(NoPaths))
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
          LeftToRight = y_coord(2) > y_coord(1) ! TO DO check if this doesnt break for special cases
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
          CrevEnd(NoPaths),CrevStart(NoPaths))! (because already created on boss)
     CALL MPI_BCAST(CrevX,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevY,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevEnd,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     CALL MPI_BCAST(CrevStart,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
     !above IF(Parallel) but that's assumed implicitly
     !make sure that code works for empty isomesh as well!!

     ! get calving polygons
     IF(Boss) THEN
      CALL GetCalvingPolygons(IsoMesh, CrevassePaths, EdgeX, EdgeY, Polygon, PolyStart, PolyEnd, gridmesh_dx)
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
           IF(inside) CYCLE ! already found
           ALLOCATE(PathPoly(2, PolyEnd(j)-PolyStart(j)+1))
           PathPoly = Polygon(:, PolyStart(j):PolyEnd(j))
           inside = PointInPolygon2D(PathPoly, (/xx, yy/))
           IF(inside) ClosestCrev = j
           DEALLOCATE(PathPoly)
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
     END IF

     Debug = .FALSE.
     CalvingValues = SignDistValues
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
       CalvingValues = 0.0_dp
       IsCalvingNode = .FALSE.
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

    !Output some information
    !CALL CalvingStats(MaxBergVolume)
    ! IF(MaxBergVolume > PauseVolumeThresh) THEN
    !   PauseSolvers = .TRUE.
    ! ELSE
    !   PauseSolvers = .FALSE.
    ! END IF

    PRINT *,'not calculating maxbergvolume now, depends on columns!'

    CALL SParIterAllReduceOR(PauseSolvers) !Really this should just be Boss sending to all
    CALL ListAddLogical( Model % Simulation, 'Calving Pause Solvers', PauseSolvers )
    IF(Debug) PRINT *,ParEnv % MyPE, ' Calving3D, pause solvers: ', PauseSolvers

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken up to output calving stats: ', rt
    rt0 = RealTime()
    !-----------------------------------------------------------
    ! Tidy up and deallocate
    !-----------------------------------------------------------

    FirstTime = .FALSE.

    PCSolver % Variable => NULL()
    PCSolver % Matrix % Perm => NULL()
    CALL FreeMatrix(PCSolver % Matrix)
    CALL ReleaseMesh(PlaneMesh)

    DEALLOCATE(TopPerm, BotPerm, LeftPerm, RightPerm, FrontPerm)

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

CONTAINS

  !Returns extent: min_x, max_x, min_y, max_y of front in rotated coord system
  FUNCTION GetFrontExtent(Mesh, RotationMatrix, DistVar, SearchDist, Buffer) RESULT(Extent)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: DistVar
    REAL(KIND=dp) :: SearchDist, RotationMatrix(3,3), Extent(4), Buffer
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

    CALL MPI_AllReduce(MPI_IN_PLACE, extent(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
    CALL MPI_AllReduce(MPI_IN_PLACE, extent(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)
    CALL MPI_AllReduce(MPI_IN_PLACE, extent(3), 1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD, ierr)
    CALL MPI_AllReduce(MPI_IN_PLACE, extent(4), 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)

    extent(1) = extent(1) - buffer
    extent(2) = extent(2) + buffer
    extent(3) = extent(3) - buffer
    extent(4) = extent(4) + buffer
  END FUNCTION GetFrontExtent

END SUBROUTINE Find_Calving3D_LSet
