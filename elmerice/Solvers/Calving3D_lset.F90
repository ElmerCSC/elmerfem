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
 
   IMPLICIT NONE

!-----------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!-----------------------------------------------
   TYPE(ValueList_t), POINTER :: Params
   TYPE(Variable_t), POINTER :: CalvingVar, &  ! TO DO rm many unused variables here
        DistVar, CIndexVar, HeightVar, CrevVar, TimestepVar,SignDistVar
   TYPE(Solver_t), POINTER :: PCSolver => NULL(), &
        VTUOutputSolver => NULL(), IsoSolver => NULL()
   TYPE(Matrix_t), POINTER :: StiffMatrix
   TYPE(Mesh_t), POINTER :: Mesh, PlaneMesh, IsoMesh, WorkMesh, WorkMesh2
   TYPE(Element_t), POINTER :: Element, WorkElements(:)
   TYPE(Nodes_t), TARGET :: WorkNodes, FaceNodesT, LeftNodes, RightNodes, FrontNodes
   TYPE(Nodes_t), POINTER :: WriteNodes
   INTEGER :: i,j,jmin,k,n,counter, dim, dummyint, TotalNodes, NoNodes, &
!        comm, ierr, Me, PEs, ExtrudedLevels, FaceNodeCount, start, fin, &
!        stride, MaxNN, Next, NodesPerLevel, LeftTgt, RightTgt, &
        comm, ierr, Me, PEs, FaceNodeCount, start, fin, &
        stride, MaxNN, Next, LeftTgt, RightTgt, &   
        county, PMeshBCNums(3), DOFs, PathCount, ValidPathCount, active,&
        WriteNodeCount, MeshBC, GroupCount, GroupStart, GroupEnd, col, &
        FrontLineCount, ShiftIdx,NoCrevNodes, NoPaths
   INTEGER, PARAMETER :: GeoUnit = 11
   INTEGER, POINTER :: CalvingPerm(:), TopPerm(:)=>NULL(), BotPerm(:)=>NULL(), &
        LeftPerm(:)=>NULL(), RightPerm(:)=>NULL(), FrontPerm(:)=>NULL(), &
        PlaneFrontPerm(:)=>NULL(), PlaneLeftPerm(:)=>NULL(), &
        PlaneRightPerm(:)=>NULL(), NodeNums(:), FrontNodeNums(:), RightNodeNums(:), &
        LeftNodeNums(:), FaceNodeNums(:)=>NULL(), DistPerm(:), ColumnPerm(:), &
        CIndexPerm(:), BList(:), WorkPerm(:), InterpDim(:),&
        OrderPerm(:), SignDistPerm(:)
   INTEGER, ALLOCATABLE :: MyFaceNodeNums(:), PFaceNodeCount(:),&
        !FNColumns(:),
        disps(:), WritePoints(:), CrevEnd(:),CrevStart(:),FrontColumnList(:) ! TO DO rm frontcolumnlist
   REAL(KIND=dp) :: FrontOrientation(3), MaxHolder, &
        RotationMatrix(3,3), UnRotationMatrix(3,3), NodeHolder(3), &
        MaxMeshDist, MeshEdgeMinLC, MeshEdgeMaxLC, MeshLCMinDist, MeshLCMaxDist,&
        Projection, CrevasseThreshold, search_eps, Norm, MinCalvingSize,&
        PauseVolumeThresh, BotZ, TopZ, prop, MaxBergVolume, dy, dz, dzdy, dxdy,&
        gradLimit, Displace, y_coord(2), ShiftTo, TempDist,MinDist, xl,xr,yl, yr, xx,yy,&
        angle,angle0,&
#ifdef USE_ISO_C_BINDINGS
        rt0, rt
#else
        rt0, rt, RealTime
#endif

   REAL(KIND=dp), POINTER :: PArray(:,:) => NULL(), &
        DistValues(:), SignDistValues(:), CIndexValues(:), WorkReal(:), &
        CalvingValues(:), ForceVector(:)
   REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), HeightDirich(:), &
        Rot_y_coords(:,:), Rot_z_coords(:,:), CrevX(:),CrevY(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, DistVarname, &
        CIndexVarName, filename_root, filename,MaskName,&
        FrontMaskName,TopMaskName,BotMaskName,LeftMaskName,RightMaskName, &
        MeshDir, PC_EqName, Iso_EqName, VTUSolverName, NameSuffix,&
        MoveMeshDir,MoveMeshFullPath
   LOGICAL :: Found, Parallel, Boss, Debug, FirstTime = .TRUE., CalvingOccurs=.FALSE., &
        SaveParallelActive, InGroup, PauseSolvers, LeftToRight, MovedOne, ShiftSecond, &
        MoveMesh=.FALSE.
   LOGICAL, POINTER :: UnfoundNodes(:)=>NULL(), PWorkLogical(:)
   LOGICAL, ALLOCATABLE :: RemoveNode(:), IMOnFront(:), IMOnSide(:), &
        DeleteMe(:), IsCalvingNode(:), WorkLogical(:)

   TYPE(CrevassePath_t), POINTER :: CrevassePaths, CurrentPath, ClosestPath

   SAVE :: FirstTime, SolverName, Params, Parallel, Boss, dim, Debug, &
        DistVarName, CIndexVarName, PC_EqName, Iso_EqName, &
        MinCalvingSize, PauseVolumeThresh, gradLimit, MoveMesh

!---------------Get Variables and Parameters for Solution-----------

   rt0 = RealTime()

   IF(FirstTime) THEN
      SolverName = "Find_Calving3D"
      Params => Solver % Values
      Parallel = (ParEnv % PEs > 1)
      Boss = (ParEnv % MyPE == 0) .OR. (.NOT. Parallel)
      Debug = .FALSE.

      dim = CoordinateSystemDimension()
      IF(dim /= 3) CALL Fatal(SolverName, "Solver only works in 3D!")

      DistVarName = ListGetString(Params,"Distance Variable Name", Found)
      IF(.NOT. Found) DistVarName = "Distance"

      CIndexVarName = ListGetString(Params,"CIndex Variable Name", Found)
      IF(.NOT. Found) CIndexVarName = "CIndex"

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

      gradLimit = ListGetConstReal(Params, "Front Gradient Limit", Found)
      IF(.NOT. Found) gradLimit = 100.0_dp

   END IF !FirstTime

   Mesh => Model % Mesh

   !TODO, could take default value
   MeshEdgeMinLC = ListGetConstReal(Params, "Calving Mesh Min LC",Found, UnfoundFatal=.TRUE.)
   MeshEdgeMaxLC = ListGetConstReal(Params, "Calving Mesh Max LC",Found, UnfoundFatal=.TRUE.)
   MeshLCMinDist = ListGetConstReal(Params, "Calving Mesh LC Min Dist",Found, UnfoundFatal=.TRUE.)
   MeshLCMaxDist = ListGetConstReal(Params, "Calving Mesh LC Max Dist",Found, UnfoundFatal=.TRUE.)
   MaxMeshDist = ListGetConstReal(Params, "Calving Search Distance",Found, UnfoundFatal=.TRUE.)
   CrevasseThreshold = ListGetConstReal(Params, "Crevasse Penetration Threshold", Found, UnfoundFatal=.TRUE.)

   DistVar => VariableGet(Model % Variables, DistVarName, .TRUE., UnfoundFatal=.TRUE.)
   DistValues => DistVar % Values
   DistPerm => DistVar % Perm

   CIndexVar => VariableGet(Model % Variables, CIndexVarName, .TRUE., UnfoundFatal=.TRUE.)
   CIndexValues => CIndexVar % Values
   CIndexPerm => CIndexVar % Perm

   !This solver's variable, two dofs, holds x and y pseudo mesh update
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

   TopMaskName = "Top Surface Mask"
   BotMaskName = "Bottom Surface Mask"
   LeftMaskName = "Left Sidewall Mask"
   RightMaskName = "Right Sidewall Mask"
   FrontMaskName = "Calving Front Mask"

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
   PArray => ListGetConstRealArray( Model % Constants,'Front Orientation', &
        Found, UnfoundFatal=.TRUE.)
   DO i=1,3
      FrontOrientation(i) = PArray(i,1)
   END DO
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

!--------------Generate 2D plane mesh-------------------

   !Use GetDomainEdge for Front, Left and Right
   !Returns full domain edge to Boss partition only
   CALL GetDomainEdge(Model, Mesh, TopPerm, LeftMaskName, &
        LeftNodes, LeftNodeNums, Parallel, Simplify=.FALSE.) ! TO DO add MinDist here?
   CALL GetDomainEdge(Model, Mesh, TopPerm, RightMaskName, &
        RightNodes, RightNodeNums, Parallel, Simplify=.FALSE.)  ! TO DO add MinDist here?
   CALL GetDomainEdge(Model, Mesh, TopPerm, FrontMaskName, &
        FrontNodes, FrontNodeNums, Parallel, Simplify=.FALSE.)  ! TO DO add MinDist here?

   !Determine whether front columns are arranged
   !left to right, and reorder if not. useful later...
   IF(Boss) THEN
     FrontLineCount = SIZE(FrontNodeNums)

     NodeHolder(1) = FrontNodes % x(1)
     NodeHolder(2) = FrontNodes % y(1)
     NodeHolder(3) = FrontNodes % z(1)
     NodeHolder = MATMUL(RotationMatrix, NodeHolder)
     y_coord(1) = NodeHolder(2)

     NodeHolder(1) = FrontNodes % x(FrontLineCount)
     NodeHolder(2) = FrontNodes % y(FrontLineCount)
     NodeHolder(3) = FrontNodes % z(FrontLineCount)
     NodeHolder = MATMUL(RotationMatrix, NodeHolder)
     y_coord(2) = NodeHolder(2)

     LeftToRight = y_coord(2) > y_coord(1)

     IF(.NOT. LeftToRight) THEN
       IF(Debug) PRINT *,'Debug, switching to LeftToRight'
       FrontNodeNums = FrontNodeNums(FrontLineCount:1:-1)
       FrontNodes % x = FrontNodes % x(FrontLineCount:1:-1)
       FrontNodes % y = FrontNodes % y(FrontLineCount:1:-1)
       FrontNodes % z = FrontNodes % z(FrontLineCount:1:-1)
     END IF
   END IF

   IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

   IF(Boss) THEN
      !Remove Left and Right nodes beyond the calving search distance
      !NOTE: technically this doesn't guarantee that resulting nodes
      !will be a contiguous subsection of the original, but it almost
      !certainly will be the case

      IF(ANY(FrontNodeNums == LeftNodeNums(1))) THEN
         LeftTgt = 1
      ELSE
         LeftTgt = LeftNodes % NumberOfNodes
      END IF

      IF(ANY(FrontNodeNums == RightNodeNums(1))) THEN
         RightTgt = 1
      ELSE
         RightTgt = RightNodes % NumberOfNodes
      END IF

      ALLOCATE(RemoveNode(LeftNodes % NumberOfNodes))
      RemoveNode = .FALSE.
      DO i=1,LeftNodes % NumberOfNodes
         IF(NodeDist2D(LeftNodes, LeftTgt, i) > MaxMeshDist) &
              RemoveNode(i) = .TRUE.
      END DO
      CALL RemoveNodes(LeftNodes, RemoveNode, LeftNodeNums)
      DEALLOCATE(RemoveNode)

      ALLOCATE(RemoveNode(RightNodes % NumberOfNodes))
      RemoveNode = .FALSE.
      DO i=1,RightNodes % NumberOfNodes
         IF(NodeDist2D(RightNodes, RightTgt, i) > MaxMeshDist) &
              RemoveNode(i) = .TRUE.
      END DO
      CALL RemoveNodes(RightNodes, RemoveNode, RightNodeNums)
      DEALLOCATE(RemoveNode)

   END IF

   !---------------------------------------------
   !
   ! Get maximum coverage of calving front (i.e. vertical shadow of front)
   !
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

    !Send calving front nodes to boss partition
    IF(Parallel) THEN

       Me = ParEnv % MyPe
       PEs = ParEnv % PEs
       comm = ELMER_COMM_WORLD

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
       CALL MPI_GATHERV(Mesh % ParallelInfo % GlobalDOFs(MyFaceNodeNums),&
            FaceNodeCount,MPI_INTEGER,&
            FaceNodeNums,PFaceNodeCount,&
            disps,MPI_INTEGER,0,comm, ierr)
       !X coords
       CALL MPI_GATHERV(Mesh % Nodes % x(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % x,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       !Y coords
       CALL MPI_GATHERV(Mesh % Nodes % y(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % y,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       !Z coords
       CALL MPI_GATHERV(Mesh % Nodes % z(MyFaceNodeNums),&
            FaceNodeCount,MPI_DOUBLE_PRECISION,&
            FaceNodesT % z,PFaceNodeCount,&
            disps,MPI_DOUBLE_PRECISION,0,comm, ierr)

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

          CALL RemoveNodes(FaceNodesT, RemoveNode, FaceNodeNums)

          IF(Debug) THEN
             PRINT *, 'Size of FaceNodeNums: ', SIZE(FaceNodeNums)
             PRINT *, 'Debug Calving3D, post removal nodes: '
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
       FaceNodesT % x = Mesh % Nodes % x(MyFaceNodeNums)
       FaceNodesT % y = Mesh % Nodes % y(MyFaceNodeNums)
       FaceNodesT % z = Mesh % Nodes % z(MyFaceNodeNums)
    END IF

    !--------------------------------------------------------------------

    !Need global mesh structure info
    IF(Parallel) THEN
       !Rather than summing NoNodes from each part, we simply find
       !the maximum global node number
       CALL MPI_AllReduce(MAXVAL(Mesh % ParallelInfo % GlobalDOFs), TotalNodes, &
            1, MPI_INTEGER, MPI_MAX, comm,ierr)
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    ELSE
       TotalNodes = NoNodes
    END IF
 
    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken up to mesh creation: ', rt
    rt0 = RealTime()

    IF(Boss) THEN

       !Produce gmsh .geo file
       WRITE(filename,'(A,A)') TRIM(filename_root), ".geo"

       OPEN( UNIT=GeoUnit, FILE=filename, STATUS='UNKNOWN')

       WriteNodeCount = SIZE(FrontNodeNums) + &
            SIZE(LeftNodeNums) + &
            SIZE(RightNodeNums) - 2

       !-------------Write Points---------------
       ALLOCATE(WritePoints(WriteNodeCount))
       counter = 1
       DO i=1,3 !cycle boundaries
          SELECT CASE(i)
          CASE(1) !left
             WriteNodes => LeftNodes
             NodeNums => LeftNodeNums
          CASE(2) !front
             WriteNodes => FrontNodes
             NodeNums => FrontNodeNums
          CASE(3) !right
             WriteNodes => RightNodes
             NodeNums => RightNodeNums
          END SELECT

          n = WriteNodes % NumberOfNodes
          IF(n /= SIZE(NodeNums)) CALL Fatal("Calving3D","Size mismatch in perm size")

          !Determine order
          IF(i==1) THEN !left edge, find which end neighbours calving front
             IF(ANY(FrontNodeNums == LeftNodeNums(1))) THEN
                start=n;fin=2;stride=-1
                next = NodeNums(1)
             ELSE IF(ANY(FrontNodeNums == LeftNodeNums(SIZE(LeftNodeNums)))) THEN
                start=1;fin=n-1;stride=1
                next = NodeNums(n)
             ELSE
                CALL Fatal(SolverName,"Problem joining up a closed loop for footprint mesh.")
             END IF
          ELSE
             IF(NodeNums(1)==next) THEN
                start=1;fin=n-1;stride=1
                IF(i==3) fin = n
                next = NodeNums(n)
             ELSE IF(NodeNums(n)==next) THEN
                start=n;fin=2;stride=-1
                IF(i==3) fin = 1
                next = NodeNums(1)
             ELSE
                PRINT *, 'i, NodeNums(1), (n), next: ',i,NodeNums(1), NodeNums(n), next
                CALL Fatal(SolverName,"Problem joining up a closed loop for footprint mesh.")
             END IF
          END IF

          DO j=start,fin,stride !cycle nodes except last, these are overlap
             WRITE( GeoUnit,'(A,I0,A,ES20.11,A,ES20.11,A,ES20.11,A,ES20.11,A)')&
                  'Point(',NodeNums(j),') = {',&
                  WriteNodes % x(j),',',&
                  WriteNodes % y(j),',',&
                  0.0,',',& !don't need z coord for footprint
                  MeshEdgeMinLC,'};'
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
       DO i=1,3 !cycle boundaries
          SELECT CASE(i)
          CASE(1) !left
             NodeNums => LeftNodeNums
             MaskName = LeftMaskName
          CASE(2) !front
             NodeNums => FrontNodeNums
             MaskName = FrontMaskName
          CASE(3) !right
             NodeNums => RightNodeNums
             MaskName = RightMaskName
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
                PMeshBCNums(i) = j !Use this later
                EXIT
             END IF
          END DO
          IF(Debug) THEN
             PRINT *, 'Debug Calving3D, BC number for ',TRIM(MaskName),' is: ',MeshBC
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
       WRITE(GeoUnit,'(A)') 'Physical Surface(3)={1};'
       !TODO, check existing number of bodies, write next, instead of '3'

       !-------------Write attractor etc--------------
       WRITE(GeoUnit,'(A)') 'Field[1] = Attractor;'
       WRITE(GeoUnit,'(A)') 'Field[1].NNodesByEdge = 100.0;'
       WRITE(GeoUnit,'(A)') 'Field[1].NodesList = {'
       DO i=1,SIZE(FrontNodeNums)-1
          WRITE(GeoUnit,'(I0,A)') FrontNodeNums(i),','
       END DO
       WRITE(GeoUnit,'(I0,A)') FrontNodeNums(SIZE(FrontNodeNums)),'};'

       WRITE(GeoUnit, '(A)') 'Field[2] = Threshold;'
       WRITE(GeoUnit, '(A)') 'Field[2].IField = 1;'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].LcMin = ',MeshEdgeMinLC,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].LcMax = ',MeshEdgeMaxLC,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].DistMin = ',MeshLCMinDist,';'
       WRITE(GeoUnit, '(A,F9.1,A)') 'Field[2].DistMax = ',MeshLCMaxDist,';'

       WRITE(GeoUnit, '(A)') 'Background Field = 2;'
       WRITE(GeoUnit, '(A)') 'Mesh.CharacteristicLengthExtendFromBoundary = 0;'

       rt = RealTime() - rt0
       IF(ParEnv % MyPE == 0) &
            PRINT *, 'Time taken to write mesh file: ', rt
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
            PRINT *, 'Time taken to execute gmsh: ', rt
       rt0 = RealTime()

       !-----------system call ElmerGrid------------------
       WRITE(Message, '(A,A,A)') "ElmerGrid 14 2 ",TRIM(filename_root),".msh"

       CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )
       IF(ierr /= 0) THEN
          WRITE(Message, '(A,i0)') "Error executing ElmerGrid, error code: ",ierr
          CALL Fatal(SolverName,Message)
       END IF


    END IF !Boss only
    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken to execute ElmerGrid: ', rt
    rt0 = RealTime()

    !Load the mesh
    MeshDir = ""

    CurrentModel % DIMENSION = 2
    PlaneMesh => LoadMesh2( Model, MeshDir, filename_root, .FALSE., 1, 0 )
    CurrentModel % DIMENSION = 3
    !NOTE: checked that planemesh exists on every PE, seems fine

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken to load mesh: ', rt
    rt0 = RealTime()

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    IF(Boss) THEN
      ! clear up mesh files
      CLOSE(GeoUnit)
      IF(MoveMesh) THEN

        !Make the directory
        TimestepVar => VariableGet( Mesh % Variables, "Timestep", .TRUE. )
        WRITE(MoveMeshFullPath,'(A,A,I4.4)') TRIM(MoveMeshDir), &
             TRIM(filename_root),INT(TimestepVar % Values(1))

        WRITE(Message,'(A,A)') "mkdir -p ",TRIM(MoveMeshFullPath)
        CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )

        !Move the files

        WRITE(Message,'(A,A,A,A)') "mv ",TRIM(filename_root),"* ",&
             TRIM(MoveMeshFullPath)
        CALL EXECUTE_COMMAND_LINE( Message, .TRUE., ierr )

      ELSE
        WRITE(Message,'(A,A,A,A,A,A)') "rm -r ",TRIM(filename)," ",&
             TRIM(filename_root),".msh ",TRIM(filename_root)
        CALL EXECUTE_COMMAND_LINE( Message, .FALSE., ierr )
      END IF
    END IF

    PlaneMesh % Name = "calving_plane"
    PlaneMesh % OutputActive = .TRUE.
    PlaneMesh % Changed = .TRUE.

    WorkMesh => Mesh % Next
    Mesh % Next => PlaneMesh

    CALL CopyIntrinsicVars(Mesh, PlaneMesh)

    rt = RealTime() - rt0
    IF(ParEnv % MyPE == 0) &
         PRINT *, 'Time taken to tidy mesh files etc: ', rt
    rt0 = RealTime()

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
    NULLIFY(WorkReal, WorkPerm)

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

       !Generate PlaneMesh perms to quickly get nodes on each boundary
       CALL MakePermUsingMask( Model, Solver, PlaneMesh, FrontMaskName, &
            .FALSE., PlaneFrontPerm, dummyint)
       CALL MakePermUsingMask( Model, Solver, PlaneMesh, LeftMaskName, &
            .FALSE., PlaneLeftPerm, dummyint)
       CALL MakePermUsingMask( Model, Solver, PlaneMesh, RightMaskName, &
            .FALSE., PlaneRightPerm, dummyint)

       ! Set ave_cindex values to 0.0 on front
       ! In fact, right at the very front, ave_cindex is undefined,
       ! but it ensures that isolines will contact the front
       DO i=1, PlaneMesh % NumberOfNodes
          IF(PlaneFrontPerm(i) > 0) CrevVar % Values(CrevVar % Perm(i)) = 0.0_dp
       END DO


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
       !isoline mesh is the last the list. We want to remove it from the
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

       ALLOCATE(IMOnFront(IsoMesh % NumberOfNodes), &
            IMOnSide(IsoMesh % NumberOfNodes))
       IMOnFront = .FALSE.; IMOnSide = .FALSE.
       search_eps = EPSILON(PlaneMesh % Nodes % x(1))

       DO i=1, IsoMesh % NumberOfNodes
          Found = .FALSE.
          DO j=1, PlaneMesh % NumberOfNodes
             IF(ABS(PlaneMesh % Nodes % x(j) - &
                  IsoMesh % Nodes % x(i)) < search_eps) THEN
                IF(ABS(PlaneMesh % Nodes % y(j) - &
                     IsoMesh % Nodes % y(i)) < search_eps) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END IF
          END DO
          IF(.NOT. Found) THEN
             WRITE(Message,'(A,i0,A)') "Unable to locate isomesh node ",i," in PlaneMesh."
             CALL Fatal(SolverName, Message)
          END IF

          IF(PlaneFrontPerm(j) > 0) IMOnFront(i) = .TRUE.
          IF((PlaneLeftPerm(j) > 0) .OR. &
               (PlaneRightPerm(j) > 0)) IMOnSide(i) = .TRUE.
       END DO

       IF(Debug) THEN
          PRINT *, 'debug, count IMOnFront: ', COUNT(IMOnFront)
          PRINT *, 'debug, count IMOnSide: ', COUNT(IMOnSide)
          PRINT *, 'debug, isomesh bulkelements,', IsoMesh % NumberOfBulkElements
          PRINT *, 'debug, isomesh boundaryelements,', IsoMesh % NumberOfBoundaryElements
          PRINT *, 'debug, size isomesh elements: ', SIZE(IsoMesh % Elements)
       END IF

       !-----------------------------------------------------------------
       ! Cycle elements, deleting any which lie wholly on the front (or side)
       !-----------------------------------------------------------------

       ALLOCATE(DeleteMe( IsoMesh % NumberOfBulkElements ))
       DeleteMe = .FALSE.
       DO i=1, IsoMesh % NumberOfBulkElements
          Element => IsoMesh % Elements(i)
          N = Element % TYPE % NumberOfNodes
          IF(ALL(IMOnFront(Element % NodeIndexes(1:N))) .OR. &
             ALL(IMOnSide(Element % NodeIndexes(1:N)))) THEN
             !delete the element, we don't need it
             DeleteMe(i) = .TRUE.
          END IF
       END DO

       PRINT *, 'debug, ', COUNT(DeleteMe), ' elements marked for deletion from IsoMesh.'

       ALLOCATE(WorkElements(COUNT(.NOT. DeleteMe)))
       WorkElements = PACK(IsoMesh % Elements, (.NOT. DeleteMe))

       PRINT *, 'debug, size of workelements: ', SIZE(WorkElements)

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

       !Find chains which make contact with front twice
       !===============================================
       ! NOTES:
       ! Because intersections lie exactly on nodes, IsoSurface
       ! creates several nodes per real node, and several 202 elems.
       ! This probably isn't a problem, but we'll see...
       ! Also, we have deleted unnecessary elements, but their nodes
       ! remain, but this also isn't a problem as InterpVarToVar
       ! cycles elements, not nodes.

       CALL FindCrevassePaths(IsoMesh, IMOnFront, CrevassePaths, PathCount)
       CALL CheckCrevasseNodes(IsoMesh, CrevassePaths)
       CALL ValidateCrevassePaths(IsoMesh, CrevassePaths, FrontOrientation, PathCount)
!Debug = .TRUE.    !Debugging statements
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
       NodeHolder(3)=0.0_dp
       DO i=1,NoPaths
          NodeHolder(1) = CrevX(CrevStart(i))
          NodeHolder(2) =CrevY(CrevStart(i))
          NodeHolder = MATMUL(RotationMatrix, NodeHolder)
          y_coord(1) = NodeHolder(2)
          NodeHolder(1) = CrevX(CrevEnd(i))
          NodeHolder(2) = CrevY(CrevEnd(i))
          NodeHolder(3) = FrontNodes % z(FrontLineCount)
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
         CrevEnd(PathCount),CrevStart(PathCount))! (because already created on boss)
CALL MPI_BCAST(CrevX,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
CALL MPI_BCAST(CrevY,NoCrevNodes,MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
CALL MPI_BCAST(CrevEnd,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr) 
CALL MPI_BCAST(CrevStart,NoPaths,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)  
!above IF(Parallel) but that's assumed implicitly 
    !make sure that code works for empty isomesh as well!!

ALLOCATE(IsCalvingNode(Mesh % NumberOfNodes))
IsCalvingNode=.FALSE.

IF (NoPaths > 0 ) THEN
   DO i=1,Solver % Mesh % NumberOfNodes !
      IF (Debug) PRINT *, 'For node i', i,' out of ',Solver % Mesh % NumberOfNodes 
      xx = Solver % Mesh % Nodes % x(i)
      yy = Solver % Mesh % Nodes % y(i)
      MinDist = MAXVAL(DistValues) ! NOTE dependency on Dist here
      ! TO DO; brute force here, checking all crevasse segments, better to find closest crev first
      DO j=1, NoPaths
         DO k=CrevStart(j), CrevEnd(j)-1
            xl=CrevX(k);yl=CrevY(k)
            xr=CrevX(k+1);yr=CrevY(k+1)
            TempDist=PointLineSegmDist2D( (/xl,yl/),(/xr,yr/), (/xx,yy/)) 
            IF(TempDist <= (ABS(MinDist)+AEPS) ) THEN ! as in ComputeDistanceWithDirection
               angle = (xr-xl)*(yy-yl)-(xx-xl)*(yr-yl) ! angle to get sign of distance function
               ! In case distance is equal, favor segment with clear angles
               ! angle0 doesn't need to be initialized as TempDist < (ABS(MinDist) - AEPS) holds
               ! at some point, angle only important if TempDist ~= MinDist
               IF( (TempDist < (ABS(MinDist) - AEPS)) .OR. (ABS(angle) > ABS(angle0)) ) THEN
                  IF( angle < 0.0) THEN ! TO DO possibly direction(i) * angle if not always same
                     MinDist = -TempDist
                     IsCalvingNode(i) = .TRUE.  
                  ELSE
                     MinDist = TempDist
                  END IF
                  angle0 = angle
                  jmin=k ! TO DO rm jmin? not necessary anymore
               END IF
            END IF
         END DO
      END DO
      SignDistValues(SignDistPerm(i)) =  MinDist
 IF(Debug)     PRINT *, 'Shortest distance to closest segment in crevassepath ',MinDist
 IF(Debug)     PRINT *, 'jmin is ',jmin
   END DO
END IF
Debug = .FALSE.
CalvingValues = SignDistValues
 
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
       !TODO, check we can just RETURN now...
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

    IF(Parallel) DEALLOCATE(MyFaceNodeNums)

    IF(Boss) THEN
       DEALLOCATE(FaceNodeNums, &
                  FaceNodesT % x, &
                  FaceNodesT % y, &
                  FaceNodesT % z, &
                  RemoveNode,&
                  PlaneFrontPerm,&
                  PlaneLeftPerm,&
                  PlaneRightPerm&
                  )

       IF(Parallel) DEALLOCATE(disps, PFaceNodeCount)

    END IF

    IF(Parallel) CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

CONTAINS

 !Subroutine to print iceberg information to a file, to be processed in python.
 !Also calculates the size of the largest iceberg and returns it
 SUBROUTINE CalvingStats(MaxBergVol)
   IMPLICIT NONE

   TYPE(Nodes_t) :: ElementNodes
   TYPE(Element_t), POINTER :: CalvingElements(:), Element
   TYPE(GaussIntegrationPoints_t) :: IntegStuff

   INTEGER :: i,j,k,county, sendcount, status(MPI_STATUS_SIZE), NoIcebergs,n,&
        col, row, elemcorners(4), countcalve, ElemBergID, start, fin,LeftIndex, RightIndex
   INTEGER, ALLOCATABLE :: disps(:), FNColumns(:), FNRows(:), FNColumnOrder(:),&
        WorkInt(:), IDVector(:)
   INTEGER, POINTER :: IcebergID(:), NodeIndexes(:)
   INTEGER, PARAMETER :: FileUnit = 57
   REAL(KIND=dp), ALLOCATABLE :: AllCalvingValues(:), MyOrderedCalvingValues(:), &
        WorkReal(:), CalvingMagnitude(:)
   REAL(KIND=dp) :: LeftMost, RightMost, s, U, V, W, Basis(Mesh % MaxElementNodes), &
        SqrtElementMetric, MaxBergVol, BergVolume, ElemVolume
   LOGICAL, POINTER :: NodesAreNeighbours(:,:), CalvingNeighbours(:,:), BergBoundaryNode(:,:)
   LOGICAL :: Visited = .FALSE., IcebergCondition(4), Debug, stat
   CHARACTER(MAX_NAME_LEN) :: FileName

   SAVE :: Visited

   Debug = .FALSE.
   MaxBergVol = 0.0_dp

   !Boss has:
   !   FaceNodesT   -  which is all the frontal nodes
   !   FaceNodeNums -  global front node numbers from all parts
   !   FNColumns    -  info about the structure of the mesh
   !   NodesPerLevel, ExtrudedLevels

   !BUT, doesn't have info about neighbours...
   !although this could be ascertained from FrontNodeNums, which are in order
! TO DO fix for unstructured
   IF(Boss) THEN

      FileName = TRIM(NameSuffix)//"_IcebergStats.txt"

      ALLOCATE(WorkReal(SUM(PFaceNodeCount)*DOFs),&
           AllCalvingValues(FaceNodesT % NumberOfNodes * DOFs),&
           FNRows(FaceNodesT % NumberOfNodes),&
           FNColumns(FaceNodesT % NumberOfNodes),&
          ! FNColumnOrder(NodesPerLevel),&
           NodesAreNeighbours(FaceNodesT % NumberOfNodes, FaceNodesT % NumberOfNodes),&
           CalvingNeighbours(FaceNodesT % NumberOfNodes, FaceNodesT % NumberOfNodes),&
           WorkInt(SIZE(FrontNodeNums)),&
           IDvector(FaceNodesT % NumberOfNodes),&
           disps(ParEnv % PEs), STAT=ierr)

      IDvector = [(i,i=1,FaceNodesT % NumberOfNodes)]

      disps(1) = 0
      DO i=2,ParEnv % PEs
         disps(i) = disps(i-1) + (PFaceNodeCount(i-1)*DOFs)
      END DO


      !FrontNodeNums (from GetDomainEdge) are ordered, and thus so are the columns
   !   WorkInt = MOD(FrontNodeNums, NodesPerLevel)
   !   DO i=1, SIZE(WorkInt)
       !  FNColumnOrder(WorkInt(i)) = i
    !  END DO

     FNRows = (FaceNodeNums - 1) 
   !   FNColumns = MOD(FaceNodeNums, NodesPerLevel)
      NodesAreNeighbours = .FALSE.

      DO i=1,FaceNodesT % NumberOfNodes
         DO j=1,FaceNodesT % NumberOfNodes
            IF(i==j) CYCLE
            IF( ABS(FNColumnOrder(FNColumns(j)) - FNColumnOrder(FNColumns(i))) > 1) CYCLE
            IF( ABS(FNRows(j) - FNRows(i)) > 1) CYCLE
            !Neighbour must be in either same column or row... i.e. no diag neighbours
            IF( (FNRows(j) /= FNRows(i)) .AND. &
                 (FNColumnOrder(FNColumns(j)) /= FNColumnOrder(FNColumns(i)))) CYCLE
            NodesAreNeighbours(i,j) = .TRUE.
         END DO
      END DO

      IF(Debug) PRINT *,'Debug CalvingStats, FaceNodes: ',FaceNodesT % NumberOfNodes,&
           ' neighbourships: ', COUNT(NodesAreNeighbours)
   END IF

   ALLOCATE(MyOrderedCalvingValues(COUNT(CalvingPerm>0)*DOFs), STAT=ierr)
   MyOrderedCalvingValues = 0.0_dp

   !Order the calving values to match front node numbers
   county = 0
   DO i=1,NoNodes
      IF(CalvingPerm(i) <= 0) CYCLE

      county = county + 1

      MyOrderedCalvingValues((county*DOFs)-2) = CalvingValues((CalvingPerm(i)*DOFs)-2)
      MyOrderedCalvingValues((county*DOFs)-1) = CalvingValues((CalvingPerm(i)*DOFs)-1)
      MyOrderedCalvingValues(county*DOFs) = CalvingValues(CalvingPerm(i)*DOFs)
   END DO

   !Gather calving var values
   sendcount = COUNT(CalvingPerm>0)*DOFs
   IF(Debug) PRINT *,ParEnv % MyPE,'send count: ',sendcount

   IF(sendcount > 0) THEN
      CALL MPI_BSEND(MyOrderedCalvingValues,sendcount, MPI_DOUBLE_PRECISION,0,&
           1000+ParEnv % MyPE, ELMER_COMM_WORLD, ierr)
   END IF

   IF(BOSS) THEN
      IF(Debug) PRINT *,'Debug, size workreal:',SIZE(WorkReal)
      DO i=1,ParEnv % PEs
         IF(PFaceNodeCount(i) <= 0) CYCLE

         start = 1+disps(i)
         fin = (disps(i)+PFaceNodeCount(i)*DOFs)
         IF(Debug) PRINT *,'Debug, ',i,' start, end',start, fin

         CALL MPI_RECV(WorkReal(start:fin), &
              PFaceNodeCount(i)*DOFs, MPI_DOUBLE_PRECISION, i-1, &
              1000+i-1, ELMER_COMM_WORLD, status, ierr)
      END DO
   END IF

   CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

   IF(Boss) THEN
      !Remove duplicates, using previously computed duplicate positions
      county = 0
      DO i=1,SIZE(RemoveNode)
         IF(RemoveNode(i)) CYCLE
         county = county+1
         AllCalvingValues((county*DOFs)-2) = WorkReal((i*DOFs)-2)
         AllCalvingValues((county*DOFs)-1) = WorkReal((i*DOFs)-1)
         AllCalvingValues((county*DOFs)) = WorkReal((i*DOFs))
      END DO

      IF(Debug) THEN
         PRINT *,'Debug, AllCalvingValues: '
         DO i=1,FaceNodesT % NumberOfNodes
            PRINT *, 'Node: ',i, 'x,y,z: ',AllCalvingValues((i*3)-2), &
                 AllCalvingValues((i*3)-1), AllCalvingValues(i*3)
         END DO
      END IF

      CalvingNeighbours = NodesAreNeighbours
      DO i=1,SIZE(CalvingNeighbours,1)
         k = i * DOFs
         IF(ALL(AllCalvingValues(k-2:k) == 0.0_dp)) THEN
            CalvingNeighbours(i,:) = .FALSE.
            CalvingNeighbours(:,i) = .FALSE.
         END IF
      END DO

      !Mark connected calving neighbours with a unique iceberg ID
      ALLOCATE(IcebergID(FaceNodesT % NumberOfNodes))
      IcebergID = 0

      NoIcebergs = 0
      DO i=1,FaceNodesT % NumberOfNodes
         k = i * DOFs
         !pretty sure no perm between CalvingValues and NodesAreNeighbours
         IF(ALL(AllCalvingValues(k-2:k) == 0.0_dp)) CYCLE
         IF(IcebergID(i) > 0) CYCLE !Got already

         !new group
         NoIcebergs = NoIcebergs + 1
         IcebergID(i) = NoIcebergs
         CALL MarkNeighbours(i, CalvingNeighbours, IcebergID, NoIcebergs)

      END DO

      !Now cycle icebergs
      DEALLOCATE(WorkInt)
      ALLOCATE(BergBoundaryNode(NoIcebergs, FaceNodesT % NumberOfNodes))
      BergBoundaryNode = .FALSE.

      !Cycle icebergs, cycle nodes in iceberg, mark all neighbours (edges of iceberg)
      DO i=1,NoIcebergs
         ALLOCATE(WorkInt(COUNT(IcebergID == i)))
         WorkInt = PACK(IDVector, (IcebergID == i))

         DO j=1,SIZE(WorkInt)
            DO k=1,SIZE(NodesAreNeighbours,1)
               IF(.NOT. NodesAreNeighbours(WorkInt(j),k)) CYCLE

               !Not a boundary node
               IF(IcebergID(k) /= 0) THEN
                  IF(IcebergID(k) /= i) CALL Fatal("CalvingStats",&
                       "This shouldn't happen - two adjacent nodes in different icebergs...")
                  CYCLE
               END IF

               BergBoundaryNode(i,k) = .TRUE.
            END DO
         END DO
         DEALLOCATE(WorkInt)
      END DO

      !------------------------------------------
      ! Scan along/down front, constructing calving elements
      !
      ! Strategy: cycle front nodes, looking for a calving element which
      !           has node i as a top left corner. NB, it needn't necessarily
      !           be a calving node itself, so long as three of the 4 are...
      !
      ! elem corners go (tl, tr, bl, br)
      !------------------------------------------
      ALLOCATE(CalvingElements(FaceNodesT % NumberOfNodes)) !<- excessive, but meh
      county = 0

      DO i=1,FaceNodesT % NumberOfNodes
         !NB use NodesAreNeighbours instead of CalvingNeighbours, check both
         elemcorners = 0
         elemcorners(1) = i
         col = FNColumnOrder(FNColumns(i)) !<-- pay attention
         row = FNRows(i)

         !gather 3 other nodes
         DO j=1,FaceNodesT % NumberOfNodes
            !note, this doesn't guarantee left, or right,
            !but as long as consistent, doesn't matter
!            IF(.NOT. NodesAreNeighbours(i,j)) CYCLE
            SELECT CASE(FNColumnOrder(FNColumns(j)) - col)
            CASE(1)
               !Next column
               SELECT CASE(row - FNRows(j))
               CASE(1)
                  !Next Row
                  elemcorners(4) = j
               CASE(0)
                  !Same Row
                  elemcorners(2) = j
               CASE DEFAULT
                  CYCLE
               END SELECT
            CASE(0)

               !Same column
               SELECT CASE(row - FNRows(j))
               CASE(1)
                  !Next row
                  elemcorners(3) = j
               CASE(0)
                  !Same row
                  !same node!!
                  CYCLE
               CASE DEFAULT
                  CYCLE
               END SELECT
            CASE DEFAULT
               CYCLE
            END SELECT
         END DO

         IF(ANY(ElemCorners == 0)) CYCLE !edge of domain

         ElemBergID = MAXVAL(IcebergID(elemcorners))

         IF(ElemBergID == 0) CYCLE

         IcebergCondition = (IcebergID(elemcorners) == ElemBergID) &
              .OR. BergBoundaryNode(ElemBergID,elemcorners)

         countcalve = COUNT(IcebergCondition)

         IF( countcalve < 3) CYCLE

         county = county + 1
         ALLOCATE(CalvingElements(county) % NodeIndexes(countcalve))
         CalvingElements(county) % TYPE => GetElementType(countcalve*100+countcalve, .FALSE.)
         CalvingElements(county) % BodyID = ElemBergID

         countcalve = 0
         DO j=1,4
            IF(IcebergCondition(j)) THEN
               countcalve = countcalve+1
               CalvingElements(county) % NodeIndexes(countcalve) = elemcorners(j)
            END IF
         END DO

      END DO

      !------------------------------------------
      ! Write info to file
      !------------------------------------------

      !Find left and rightmost nodes for info
      LeftMost = HUGE(0.0_dp)
      RightMost = -HUGE(0.0_dp) !??

      DO i=1,FaceNodesT % NumberOfNodes
         Nodeholder(1) = FaceNodesT % x(i)
         Nodeholder(2) = FaceNodesT % y(i)
         Nodeholder(3) = FaceNodesT % z(i)
         Nodeholder = MATMUL(RotationMatrix, NodeHolder)

         IF(Nodeholder(2) < LeftMost) THEN
            LeftIndex = i
            LeftMost = NodeHolder(2)
         END IF

         IF(Nodeholder(2) > RightMost) THEN
            RightIndex = i
            RightMost = NodeHolder(2)
         END IF
      END DO

      IF(Visited) THEN
         OPEN( UNIT=FileUnit, FILE=filename, STATUS='UNKNOWN', ACCESS='APPEND')
      ELSE
         OPEN( UNIT=FileUnit, FILE=filename, STATUS='UNKNOWN')
         WRITE(FileUnit, '(A,ES20.11,ES20.11,ES20.11)') "FrontOrientation: ",FrontOrientation
      END IF

      !Write out the left and rightmost points
      WRITE(FileUnit, '(A,i0,ES30.21)') 'Time: ',GetTimestep(),GetTime()
      WRITE(FileUnit, '(A,ES20.11,ES20.11)') 'Left (xy): ',&
           FaceNodesT % x(LeftIndex),&
           FaceNodesT % y(LeftIndex)

      WRITE(FileUnit, '(A,ES20.11,ES20.11)') 'Right (xy): ',&
           FaceNodesT % x(RightIndex),&
           FaceNodesT % y(RightIndex)

      !Write the iceberg count
      WRITE(FileUnit, '(A,i0)') 'Icebergs: ',NoIcebergs

      !TODO, write element count
      !  would need to modify Icebergs.py too
      DO i=1,NoIcebergs
         county = 0
         !count elements

         WRITE(FileUnit, '(A,i0)') 'Iceberg ',i
         WRITE(FileUnit, '(i0,A,i0)') COUNT(BergBoundaryNode(i,:))," ", COUNT(IceBergID == i)
         WRITE(FileUnit, '(A)') "Boundary nodes"
         DO j=1,FaceNodesT % NumberOfNodes
            IF(BergBoundaryNode(i,j)) THEN

               county = county + 1
               WRITE(FileUnit,'(i0,A,i0,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11)') &
                    county," ",j,&
                    FaceNodesT % x(j), FaceNodesT % y(j), FaceNodesT % z(j),&
                    0.0_dp, 0.0_dp, 0.0_dp

            END IF
         END DO

         WRITE(FileUnit, '(A)') "Calving nodes"
         DO j=1,FaceNodesT % NumberOfNodes
            IF(IcebergID(j) == i) THEN
               county = county + 1
               k = j*DOFs
               WRITE(FileUnit,'(i0,A,i0,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11)')&
                    county," ",j,&
               FaceNodesT % x(j),&
               FaceNodesT % y(j),&
               FaceNodesT % z(j),&
               AllCalvingValues(k-2), &
               AllCalvingValues(k-1), &
               AllCalvingValues(k)
            END IF
         END DO

         WRITE(FileUnit, '(A)') "Elements"
         DO j=1,SIZE(CalvingElements)
            IF(CalvingElements(j) % BodyID /= i) CYCLE
            DO k=1,CalvingElements(j) % TYPE % NumberOfNodes
               WRITE(FileUnit,'(i0,A)', ADVANCE="NO") CalvingElements(j) % NodeIndexes(k),'  '
            END DO
            WRITE(FileUnit,'(A)') ''
         END DO

      END DO

      CLOSE(FileUnit)

      !------------------------------------------
      ! Compute berg volumes. Largest berg size determines
      ! whether the pause the timestep.
      !------------------------------------------
      n = Mesh % MaxElementNodes
      ALLOCATE(ElementNodes % x(n),&
           ElementNodes % y(n),&
           ElementNodes % z(n),&
           CalvingMagnitude(FaceNodesT % NumberOfNodes))

      DO i=1,SIZE(CalvingMagnitude)
        k = i*DOFs
        CalvingMagnitude(i) = ((AllCalvingValues(k-2) ** 2) + &
             (AllCalvingValues(k-1) ** 2) + &
             (AllCalvingValues(k) ** 2)) ** 0.5
      END DO

      MaxBergVol = 0.0_dp
      DO i=1,NoIcebergs
        BergVolume = 0.0_dp

        DO j=1,SIZE(CalvingElements)
          Element => CalvingElements(j)

          IF(Element % BodyID /= i) CYCLE

          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes

          ElementNodes % x(1:n) = FaceNodesT % x(NodeIndexes(1:n))
          ElementNodes % y(1:n) = FaceNodesT % y(NodeIndexes(1:n))
          ElementNodes % z(1:n) = FaceNodesT % z(NodeIndexes(1:n))

          IntegStuff = GaussPoints( Element )

          ElemVolume = 0.0_dp
          DO k=1,IntegStuff % n

            U = IntegStuff % u(k)
            V = IntegStuff % v(k)
            W = IntegStuff % w(k)

            stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
                 Basis )

            !assume cartesian here
            s = SqrtElementMetric * IntegStuff % s(k)

            ElemVolume = ElemVolume + s * SUM(CalvingMagnitude(NodeIndexes(1:n)) * Basis(1:n))
          END DO

          BergVolume = BergVolume + ElemVolume
        END DO
        IF(Debug) PRINT *,'Berg ',i,' volume: ', BergVolume
        MaxBergVol = MAX(BergVolume, MaxBergVol)
      END DO
      IF(Debug) PRINT *,'Max berg volume: ',MaxBergVol
    END IF


    Visited = .TRUE.
    DEALLOCATE(MyOrderedCalvingValues)
    IF(Boss) THEN
      DEALLOCATE(WorkReal,&
           AllCalvingValues, &
           FNRows,&
           FNColumns,&
           FNColumnOrder,&
           NodesAreNeighbours,&
           CalvingNeighbours,&
           disps,&
           BergBoundaryNode,&
           IDVector,&
           IcebergID,&
           ElementNodes % x,&
           ElementNodes % y,&
           ElementNodes % z,&
           CalvingMagnitude&
           )

      !Cycle and deallocate element % Nodeindexes, and elements
      DO i=1,SIZE(CalvingElements)
        IF(ASSOCIATED(CalvingElements(i) % NodeIndexes)) &
             DEALLOCATE(CalvingElements(i) % NodeIndexes)
      END DO
      DEALLOCATE(CalvingElements)
    END IF

  END SUBROUTINE CalvingStats

  
END SUBROUTINE Find_Calving3D_LSet
