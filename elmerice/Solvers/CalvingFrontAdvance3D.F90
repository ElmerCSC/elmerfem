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

!A routine for getting frontal advance in 3D, given velocity and melt

 SUBROUTINE FrontAdvance3D ( Model, Solver, dt, TransientSimulation )

   USE CalvingGeometry
   USE DefUtils
   IMPLICIT NONE

!-----------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!-----------------------------------------------
   TYPE(Mesh_t), POINTER :: Mesh
   TYPE(Nodes_t), TARGET :: FrontNodes, ColumnNodes
   TYPE(ValueList_t), POINTER :: Params
   TYPE(Variable_t), POINTER :: Var, VeloVar, MeltVar, NormalVar, TangledVar
   INTEGER :: i, j, k, m, n, DOFs, TotalNodes, NodesPerLevel, ExtrudedLevels,&
        FaceNodeCount, DummyInt, col, hits, ierr, FrontLineCount, county,&
        ShiftIdx, ShiftToIdx, NoTangledGroups, PivotIdx, CornerIdx, &
        SecondIdx, FirstTangleIdx, LastTangleIdx
   INTEGER, POINTER :: Perm(:), FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        FrontNodeNums(:)=>NULL()
   INTEGER, ALLOCATABLE :: FNColumns(:), FrontColumnList(:), FrontLocalNodeNumbers(:), &
        NodeNumbers(:), TangledGroup(:), TangledPivotIdx(:), UpdatedDirection(:)
   REAL(KIND=dp) :: FrontOrientation(3),RotationMatrix(3,3),&
        NodeVelo(3), NodeMelt(3), NodeNormal(3), ColumnNormal(3), CornerDirection(2),&
        MeltRate, Displace(3), NodeHolder(3), DangerGrad, ShiftTo, direction, &
        ShiftDist, y_coord(2), epsShift, ShiftToY, LongRangeLimit, MaxDisplacement, LimitZ, &
        p1(2),p2(2),q1(2),q2(2),intersect(2), LeftY, RightY, EpsTangle,thisEps,Shift, thisY
   REAL(KIND=dp), POINTER :: PArray(:,:) => NULL(), Advance(:)
   REAL(KIND=dp), ALLOCATABLE :: Rot_y_coords(:,:), Rot_z_coords(:,:), ColumnNormals(:,:), &
        TangledShiftTo(:)
   LOGICAL :: Found, Debug, Parallel, Boss, ShiftLeft, LeftToRight, MovedOne, ShiftSecond, &
        Protrusion, SqueezeLeft, SqueezeRight, FirstTime=.TRUE., intersect_flag, FrontMelting
   LOGICAL, ALLOCATABLE :: DangerZone(:), WorkLogical(:), UpdatedColumn(:),&
        Tangled(:), DontMove(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VeloVarName, MeltVarName, &
        NormalVarName, FrontMaskName, TopMaskName, TangledVarName

   SAVE :: FirstTime

   !-----------------------------------------------
   ! Initialisation
   !-----------------------------------------------

   Debug = .FALSE.
   Parallel = (ParEnv % PEs > 1)
   Boss = ParEnv % MyPE == 0

   SolverName = "FrontAdvance3D"
   Params => Solver % Values
   Mesh => Solver % Mesh

   !The main solver var contains the magnitude of front advance
   Var => Solver % Variable
   Advance => Var % Values
   Perm => Var % Perm
   DOFs = Var % DOFs
   IF(Var % DOFs /= 3) CALL Fatal(SolverName, "Variable should have 3 DOFs...")

   !Get the flow solution
   VeloVarName = ListGetString(Params, "Flow Solution Variable Name", Found)
   IF(.NOT. Found) THEN
     CALL Info(SolverName, "Flow Solution Variable Name not found, assuming 'Flow Solution'")
     VeloVarName = "Flow Solution"
   END IF
   VeloVar => VariableGet(Mesh % Variables, VeloVarName, .TRUE., UnfoundFatal=.TRUE.)

   !Get melt rate
   MeltVarName = ListGetString(Params, "Melt Variable Name", Found)
   IF(.NOT. Found) THEN
     CALL Info(SolverName, "Melt Variable Name not found, assuming no frontal melting")
     FrontMelting = .FALSE.
   ELSE
     FrontMelting = .TRUE.
     MeltVar => VariableGet(Mesh % Variables, MeltVarName, .TRUE., UnfoundFatal=.TRUE.)
   END IF

   TangledVarName = "Tangled"
   TangledVar => VariableGet(Mesh % Variables, TangledVarName, .TRUE., UnfoundFatal=.TRUE.)
   TangledVar % Values = 0.0_dp

   !Get front normal vector
   NormalVarName = ListGetString(Params, "Normal Vector Variable Name", UnfoundFatal=.TRUE.)
   NormalVar => VariableGet(Mesh % Variables, NormalVarName, .TRUE., UnfoundFatal=.TRUE.)

   !Get the orientation of the calving front, compute rotation matrix
   !TODO: generalize and link
   PArray => ListGetConstRealArray( Model % Constants,'Front Orientation', Found, UnfoundFatal=.TRUE.)
   DO i=1,3
      FrontOrientation(i) = PArray(i,1)
   END DO
   RotationMatrix = ComputeRotationMatrix(FrontOrientation)

   DangerGrad = ListGetConstReal(Params, "Front Gradient Threshold", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Front Gradient Threshold' not found, setting to 5.0")
     DangerGrad = 5.0
   END IF

   !When correcting for unprojectability, to prevent interp problems,
   !ensure that it is SLIGHTLY beyond beyond parallel, by adding this
   !value (metres) to the displacement
   epsShift = ListGetConstReal(Params, "Front Projectability Epsilon", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Front Projectability Epsilon' not found, setting to 1.0")
     epsShift = 1.0
   END IF

   epsTangle = ListGetConstReal(Params, "Front Tangle Epsilon", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Front Tangle Epsilon' not found, setting to 1.0")
     epsTangle = 1.0
   END IF

   LongRangeLimit = ListGetConstReal(Params, "Column Max Longitudinal Range", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Column Max Longitudinal Range' not found, setting to 300.0")
     LongRangeLimit = 300.0_dp
   END IF

   MaxDisplacement = ListGetConstReal(Params, "Maximum Node Displacement", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Maximum Node Displacement' not found, setting to 1.0E4.")
     MaxDisplacement = 1.0E4_dp
   END IF

   !-------------------------------------------
   ! Find FNColumns
   CALL MPI_AllReduce(MAXVAL(Mesh % ParallelInfo % GlobalDOFs), TotalNodes, &
        1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD,ierr)

   ExtrudedLevels = ListGetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found)
   IF(.NOT. Found) ExtrudedLevels = &
        ListGetInteger(CurrentModel % Simulation,'Remesh Extruded Mesh Levels',Found)
   IF(.NOT. Found) CALL Fatal(SolverName,&
        "Unable to find 'Extruded Mesh Levels' or 'Remesh Extruded Mesh Levels'")
   NodesPerLevel = TotalNodes / ExtrudedLevels

   ALLOCATE(FNColumns(Mesh % NumberOfNodes))
   FNColumns = MOD(Mesh % ParallelInfo % GlobalDOFs, NodesPerLevel)

   !Get the front line
   FrontMaskName = "Calving Front Mask"
   TopMaskName = "Top Surface Mask"

   CALL MakePermUsingMask( Model, Solver, Mesh, TopMaskName, &
        .FALSE., TopPerm, dummyint)

   CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
        .FALSE., FrontPerm, FaceNodeCount)

   CALL GetDomainEdge(Model, Mesh, TopPerm, FrontMaskName, &
        FrontNodes, FrontNodeNums, Parallel, Simplify=.FALSE.)

   !Pass FrontNodeNums to all CPUs
   IF(Boss) FrontLineCount = SIZE(FrontNodeNums)
   CALL MPI_BCAST(FrontLineCount , 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

   !Determine whether front columns are arranged
   !left to right, and reorder if not
   IF(Boss) THEN

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

   IF(.NOT. Boss) ALLOCATE(FrontNodeNums(FrontLineCount))
   CALL MPI_BCAST(FrontNodeNums , FrontLineCount, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

   !FrontNodeNums is the GlobalDOFs of the nodes on the front
   !So this can be easily converted into a FNColumn like list
   ALLOCATE(FrontColumnList(FrontLineCount), &
        FrontLocalNodeNumbers(FrontLineCount), &
        DangerZone(FrontLineCount), &
        ColumnNormals(FrontLineCount,3))

   FrontColumnList = MOD(FrontNodeNums, NodesPerLevel)
   FrontLocalNodeNumbers = -1
   DangerZone = .FALSE.
   ColumnNormals = 0.0_dp

   DO i=1,FrontLineCount
     n = FrontNodeNums(i)
     DO j=1,Mesh % NumberOfNodes
       IF(Mesh % ParallelInfo % GlobalDOFs(j) == n) THEN
         FrontLocalNodeNumbers(i) = j
       END IF
     END DO
   END DO

   !--------------------------------------
   ! Action: Compute lagrangian displacement for all nodes
   !         This is the main function of the Solver.
   !         Everything following this DO loop is taking care of
   !         unprojectability/high gradient etc
   !--------------------------------------

   DO i=1,Mesh % NumberOfNodes
     IF(Perm(i) <= 0) CYCLE

     IF(FrontMelting) THEN
       IF(MeltVar % Perm(i) <= 0) &
            CALL Fatal(SolverName, "Permutation error on front node!")


       !Scalar melt value from Plume solver
       MeltRate = MeltVar % Values(MeltVar % Perm(i))

       NodeNormal(1) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 1)
       NodeNormal(2) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 2)
       NodeNormal(3) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 3)

       NodeMelt = NodeNormal * MeltRate
     ELSE
       NodeMelt = 0.0_dp
     END IF

     !Compute front normal component of velocity
     NodeVelo(1) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 1)
     NodeVelo(2) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 2)
     NodeVelo(3) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 3)


     Displace = 0.0

     Displace(1) =  NodeVelo(1) - NodeMelt(1)
     Displace(2) =  NodeVelo(2) - NodeMelt(2)
     Displace(3) =  NodeVelo(3) - NodeMelt(3)

     Displace = Displace * dt

     IF(MAXVAL(Displace) > MaxDisplacement) THEN
       WRITE(Message,'(A,i0,A)') "Maximum allowable front displacement exceeded for node ",i,". Limiting..."
       CALL Warn(SolverName, Message)
       Displace = Displace * (MaxDisplacement/MAXVAL(Displace))
     END IF

     Advance((Perm(i)-1)*DOFs + 1) = Displace(1)
     Advance((Perm(i)-1)*DOFs + 2) = Displace(2)
     Advance((Perm(i)-1)*DOFs + 3) = Displace(3)

     !First and last (i.e. lateral margin) columns don't move
     IF((FNColumns(i) == FrontColumnList(1)) .OR. &
          (FNColumns(i) == FrontColumnList(FrontLineCount))) THEN
       Advance((Perm(i)-1)*DOFs + 1) = 0.0_dp
       Advance((Perm(i)-1)*DOFs + 2) = 0.0_dp
       Advance((Perm(i)-1)*DOFs + 3) = 0.0_dp
     END IF
   END DO


   !----------------------------------------------------------
   ! Now we need to look for and limit  regions on the
   ! where high gradient + 'straightness' makes
   ! remeshing into vertical columns troublesome.
   !----------------------------------------------------------
   !
   ! Five things:
   !  - Limit horizontal range (otherwise, Remesh will smear a column of nodes
   !    all over the place)
   !  - Limit undercut range (longitudinal range) to prevent mesh degeneracy
   !  - Prevent retreat at node near lateral margins
   !  - Check for 'tangled' regions, where fixing unprojectability one way
   !    causes unprojectability the other way.
   !  - Prevent Unprojectability - requires knowledge of neighbouring columns
   !----------------------------------------------------------

   !-------------------------------------------
   ! Cycle columns, looking at average normal for the column
   !-------------------------------------------
   ! Normal vector is limited by projectability, so if average is high,
   ! all are high
   !
   ! Note that we are cycling the parallel gathered line, so each part
   ! will only have a few columns (if any)
   ! So, mark troublesome columns where present, then communicate
   !-------------------------------------------

   DO i=1,FrontLineCount
     col = FrontColumnList(i)

     hits = COUNT(FNColumns == col)
     IF(hits == 0) CYCLE

     ColumnNormal = 0.0_dp

     !Gather normal vector
     DO j=1,Mesh % NumberOfNodes
       IF(FNColumns(j) == col) THEN
         ColumnNormal(1) = ColumnNormal(1) + &
              NormalVar % Values(((NormalVar % Perm(j)-1)*NormalVar % DOFs) + 1)
         ColumnNormal(2) = ColumnNormal(2) + &
              NormalVar % Values(((NormalVar % Perm(j)-1)*NormalVar % DOFs) + 2)
         ColumnNormal(3) = ColumnNormal(3) + &
              NormalVar % Values(((NormalVar % Perm(j)-1)*NormalVar % DOFs) + 3)
       END IF
     END DO

     !unit vector
     ColumnNormal = ColumnNormal / hits
     !Convert to front coordinate system
     ColumnNormal = MATMUL( RotationMatrix, ColumnNormal )
     !Save for later
     ColumnNormals(i,1:3) = ColumnNormal(1:3)

     !We're concerned with the lateral (rather than vertical) component
     IF( ABS(ColumnNormal(2) / ColumnNormal(3)) > DangerGrad ) THEN
       DangerZone(i) = .TRUE.
     END IF
   END DO

   !Communicate between partitions
   !This is an 'OR' condition insofar as only one partition
   !will have actually computed the gradient
   IF(Parallel) THEN
     DO i=1,FrontLineCount
       CALL SParIterAllReduceOR(DangerZone(i))
     END DO
   END IF

   !Mark before and after dangerzone too
   ALLOCATE(WorkLogical(FrontLineCount))
   WorkLogical = .FALSE.

   DO i=1,FrontLineCount
     IF(DangerZone(i)) WorkLogical(i) = .TRUE.

     IF(i > 1) THEN
       IF(DangerZone(i-1)) WorkLogical(i) = .TRUE.
     END IF

     IF(i < FrontLineCount) THEN
       IF(DangerZone(i+1)) WorkLogical(i) = .TRUE.
     END IF
   END DO

   DangerZone = WorkLogical
   DEALLOCATE(WorkLogical)

   !Find and store leftmost and rightmost (rotated) coordinate of
   !each column for checking projectability later
   !Rot_y_coords(:,1) is leftmost, and (:,2) is right most y coord of
   !each column's nodes, and Rot_z_coords(:,:) is the min(1)/max(2) z
   ALLOCATE(Rot_y_coords(FrontLineCount,2),&
        Rot_z_coords(FrontLineCount,2))
   Rot_y_coords(:,1) = HUGE(0.0_dp)
   Rot_y_coords(:,2) = -HUGE(0.0_dp)
   Rot_z_coords(:,1) = HUGE(0.0_dp)
   Rot_z_coords(:,2) = -HUGE(0.0_dp)

   DO i=1,FrontLineCount
     col = FrontColumnList(i)
     IF(COUNT(FNColumns == col) == 0) CYCLE

     DO j=1,Mesh % NumberOfNodes
       IF(FNColumns(j) /= col) CYCLE

       NodeHolder(1) = Mesh % Nodes % x(j) + Advance((Perm(j)-1)*DOFs + 1)
       NodeHolder(2) = Mesh % Nodes % y(j) + Advance((Perm(j)-1)*DOFs + 2)
       NodeHolder(3) = Mesh % Nodes % z(j) + Advance((Perm(j)-1)*DOFs + 3)
       NodeHolder = MATMUL(RotationMatrix, NodeHolder)

       Rot_y_coords(i,1) = MIN(Rot_y_coords(i,1), NodeHolder(2))
       Rot_y_coords(i,2) = MAX(Rot_y_coords(i,2), NodeHolder(2))

       Rot_z_coords(i,1) = MIN(Rot_z_coords(i,1), NodeHolder(3))
       Rot_z_coords(i,2) = MAX(Rot_z_coords(i,2), NodeHolder(3))
     END DO
   END DO

   !-----------------------------------
   ! Limit lateral range (to 0) in DangerZone
   !-----------------------------------
   ! The physical interpretation of this is that, when melt undercutting occurs
   ! in a region of high gradient, the unmelted parts of that column also shift
   ! with the melt.
   ! This is a minor limitation, but better than Free Surface Equation!

   DO i=1,FrontLineCount
     IF(.NOT. DangerZone(i)) CYCLE

     col = FrontColumnList(i)
     hits = COUNT(FNColumns == col)

     IF(hits == 0) CYCLE

     ALLOCATE(NodeNumbers(hits), &
          ColumnNodes % x(hits),&
          ColumnNodes % y(hits),&
          ColumnNodes % z(hits))

     ColumnNodes % NumberOfNodes = hits

     !Gather nodenumbers in column
     county = 0
     DO j=1,Mesh % NumberOfNodes
       IF(FNColumns(j) /= col) CYCLE
       county = county + 1

       NodeNumbers(county) = j
       ColumnNodes % x(county) = Mesh % Nodes % x(j) + Advance((Perm(j)-1)*DOFs + 1)
       ColumnNodes % y(county) = Mesh % Nodes % y(j) + Advance((Perm(j)-1)*DOFs + 2)
       ColumnNodes % z(county) = Mesh % Nodes % z(j) + Advance((Perm(j)-1)*DOFs + 3)
     END DO

     !direction - which way is this part of the front pointing?
     ShiftLeft = ColumnNormals(i,2) > 0
     IF(ShiftLeft) THEN
       ShiftTo = HUGE(ShiftTo)
     ELSE
       ShiftTo = -HUGE(ShiftTo)
     END IF

     !Rotate points and find furthest left (or right)
     DO j=1,ColumnNodes % NumberOfNodes
       NodeHolder(1) = ColumnNodes % x(j)
       NodeHolder(2) = ColumnNodes % y(j)
       NodeHolder(3) = ColumnNodes % z(j)
       NodeHolder = MATMUL(RotationMatrix, NodeHolder)

       IF(ShiftLeft) THEN
         IF(NodeHolder(2) < ShiftTo) ShiftTo = NodeHolder(2)
       ELSE
         IF(NodeHolder(2) > ShiftTo) ShiftTo = NodeHolder(2)
       END IF
     END DO

     !Now, for each node in column, compute the displacement (in rotated y coordinate)
     DO j=1,ColumnNodes % NumberOfNodes
       NodeHolder(1) = ColumnNodes % x(j)
       NodeHolder(2) = ColumnNodes % y(j)
       NodeHolder(3) = ColumnNodes % z(j)
       NodeHolder = MATMUL(RotationMatrix, NodeHolder)

       ShiftDist = ShiftTo - NodeHolder(2)

       Displace = 0.0_dp
       Displace(2) = ShiftDist

       Displace = MATMUL(TRANSPOSE(RotationMatrix), Displace)

       !Adjust the variable values to shift nodes in line.
       Advance((Perm(NodeNumbers(j))-1)*DOFs + 1) = &
            Advance((Perm(NodeNumbers(j))-1)*DOFs + 1) + Displace(1)
       Advance((Perm(NodeNumbers(j))-1)*DOFs + 2) = &
            Advance((Perm(NodeNumbers(j))-1)*DOFs + 2) + Displace(2)

       IF(Debug) PRINT *,'node: ',NodeNumbers(j),' shifting: ',ShiftDist,' to ',ShiftTo,' point: ', &
            ColumnNodes % x(j),ColumnNodes % y(j),ColumnNodes % z(j)
     END DO

     DEALLOCATE(NodeNumbers, &
          ColumnNodes % x, &
          ColumnNodes % y, &
          ColumnNodes % z)
   END DO

   !-----------------------------------
   ! Limit longitudinal range everywhere
   !-----------------------------------
   ! Two options for how to do this:
   ! Drag nodes back with melt, or effectively
   ! limit melt (keeping nodes forward)
   ! Opt for the latter, or else we're
   ! forcing melting to have an effect
   !-----------------------------------
   DO i=1,FrontLineCount

     col = FrontColumnList(i)
     hits = COUNT(FNColumns == col)

     IF(hits == 0) CYCLE

     IF((Rot_z_coords(i,2) - Rot_z_coords(i,1)) > LongRangeLimit) THEN

       LimitZ = Rot_z_coords(i,2) - LongRangeLimit
       Rot_z_coords(i,1) = LimitZ

       DO j=1,Mesh % NumberOfNodes
         IF(FNColumns(j) /= col) CYCLE

         NodeHolder(1) = Mesh % Nodes % x(j) + Advance((Perm(j)-1)*DOFs + 1)
         NodeHolder(2) = Mesh % Nodes % y(j) + Advance((Perm(j)-1)*DOFs + 2)
         NodeHolder(3) = Mesh % Nodes % z(j) + Advance((Perm(j)-1)*DOFs + 3)
         NodeHolder = MATMUL(RotationMatrix, NodeHolder)

         IF(NodeHolder(3) >= LimitZ) CYCLE

         Displace = 0.0_dp
         Displace(3) = LimitZ - NodeHolder(3)

         Displace = MATMUL(TRANSPOSE(RotationMatrix), Displace)

         !Adjust the variable values to shift nodes in line.
         Advance((Perm(j)-1)*DOFs + 1) = &
              Advance((Perm(j)-1)*DOFs + 1) + Displace(1)
         Advance((Perm(j)-1)*DOFs + 2) = &
              Advance((Perm(j)-1)*DOFs + 2) + Displace(2)
         IF(Debug) PRINT *,ParEnv % MyPE,' Debug, limiting node range: col ',&
              i,' node ',j,' limit: ',LimitZ,' displace: ',Displace
       END DO

     END IF
   END DO


   !-----------------------------------
   ! Now check projectability
   !-----------------------------------
   ! If an unprojectable portion of the front is found,
   ! the inland nodes remain in place, while those further
   ! advanced shift sideways to maintain projectability
   !
   ! Typical case is melt plume at end of recently calved
   ! straight edge.
   !
   ! requires MPI comms

   !Gather rotated y coord of all columns
   DO i=1,FrontLineCount
     CALL MPI_AllReduce(MPI_IN_PLACE, Rot_y_coords(i,1), &
          1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD,ierr)
     CALL MPI_AllReduce(MPI_IN_PLACE, Rot_y_coords(i,2), &
          1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD,ierr)

     CALL MPI_AllReduce(MPI_IN_PLACE, Rot_z_coords(i,1), &
          1, MPI_DOUBLE_PRECISION, MPI_MIN, ELMER_COMM_WORLD,ierr)
     CALL MPI_AllReduce(MPI_IN_PLACE, Rot_z_coords(i,2), &
          1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD,ierr)

     IF(Boss .AND. Debug) PRINT *,'Debug, rot_y_coords: ',i,rot_y_coords(i,:)
     IF(Boss .AND. Debug) PRINT *,'Debug, rot_z_coords: ',i,rot_z_coords(i,:)
   END DO

   ALLOCATE(UpdatedColumn(FrontLineCount), &
        UpdatedDirection(FrontLineCount), &
        Tangled(FrontLineCount),&
        TangledGroup(FrontLineCount),&
        TangledPivotIdx(FrontLineCount),&
        TangledShiftTo(FrontLineCount),&
        DontMove(FrontLineCount))

   UpdatedColumn = .FALSE.
   UpdatedDirection = 0
   Tangled = .FALSE.
   TangledGroup = 0
   TangledPivotIdx = 0
   TangledShiftTo = 0.0_dp
   DontMove = .FALSE.

   DontMove(1) = .TRUE.
   DontMove(FrontLineCount) = .TRUE.

   !TODO: Deal with DontMove and tangling properly...
   ! Ensure that comms is not necessary for DontMove

   !Test for thin sliver at lateral margins - potential to get
   !stuck like this, so impose an angle here
   !We require that the 2nd node from the end not retreat beyond
   !a 45 degree angle inland. PYTHAGORAS
   !Shift the 2nd column in the lateral direction only,
   !to the point where it makes a 45 degree angle w/
   !the corner.
   !         *
   !       / |
   !     /   |
   !    * <- *

   DO i=1,2

     IF(i==1) THEN
       CornerIdx = 1
       SecondIdx = 2
       direction = 1.0_dp
     ELSE
       CornerIdx = FrontLineCount
       SecondIdx = FrontLineCount-1
       direction = -1.0_dp
     END IF

     !Conveniently, first time round we want the left most coord (Rot_y_coords(:,1))
     !whereas second time we want the right (Rot_y_cords(:,2))
     CornerDirection(1) = Rot_y_coords(SecondIdx,i) - Rot_y_coords(CornerIdx,1)
     CornerDirection(2) = Rot_z_coords(SecondIdx,i) - Rot_z_coords(CornerIdx,1)
     !Unit vector
     CornerDirection = CornerDirection / (SUM(CornerDirection ** 2.0)**0.5)

     IF(Debug) PRINT *,ParEnv % MyPE, 'CornerDirection ',i,': ',CornerDirection

     IF( CornerDirection(2) < (-1.0_dp / (2.0_dp ** 0.5)) ) THEN
       !Shift 2nd node and mark DontMove

       ShiftToY = Rot_y_coords(CornerIdx,1) + &
            direction * ((Rot_z_coords(CornerIdx,1) - Rot_z_coords(SecondIdx,1)))

       Rot_y_coords(SecondIdx,1) = ShiftToY
       Rot_y_coords(SecondIdx,2) = ShiftToY

       DontMove(SecondIdx) = .TRUE.
       IF(Debug) PRINT *,ParEnv % MyPE,' debug, side ',i,' of front is too inwards'

       DO k=1,Mesh % NumberOfNodes
         IF(FNColumns(k) /= FrontColumnList(SecondIdx)) CYCLE

         NodeHolder(1) = Mesh % Nodes % x(k) + Advance((Perm(k)-1)*DOFs + 1)
         NodeHolder(2) = Mesh % Nodes % y(k) + Advance((Perm(k)-1)*DOFs + 2)
         NodeHolder(3) = Mesh % Nodes % z(k) + Advance((Perm(k)-1)*DOFs + 3)

         NodeHolder = MATMUL(RotationMatrix, NodeHolder)

         Displace = 0.0_dp
         Displace(2) = ShiftToY - NodeHolder(2)

         Displace = MATMUL(TRANSPOSE(RotationMatrix), Displace)

         IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, shifting next to corner node ',k,&
              NodeHolder,' by ',Displace

         !Adjust the variable values to shift nodes in line.
         Advance((Perm(k)-1)*DOFs + 1) = Advance((Perm(k)-1)*DOFs + 1) + Displace(1)
         Advance((Perm(k)-1)*DOFs + 2) = Advance((Perm(k)-1)*DOFs + 2) + Displace(2)
       END DO

     END IF
   END DO

   !-----------------------------------------------
   ! Cycle the line, looking for unprojectable
   ! (but not tangled) regions.
   ! e.g. plume melting behind a vertical portion
   !-----------------------------------------------

   !In case shifting nodes affects another partition, we loop so long as
   !at least one partition has MovedOne
   MovedOne = .TRUE.
   county = 0
   DO WHILE(MovedOne)
     county = county + 1
     IF(county > 100) CALL Fatal(SolverName, "Infinite loop!")
     IF(Debug) PRINT *,ParEnv % MyPE, 'Debug, iterating projectability ', county

     MovedOne = .FALSE.
     UpdatedColumn = .FALSE.

     DO i=2,FrontLineCount

       !If already detangled, skip
       !IF(Tangled(i)) CYCLE

       !epsShift here
       IF((Rot_y_coords(i,1) - Rot_y_coords(i-1,2)) < 0.9*epsShift) THEN

         IF(Debug) PRINT *, 'Debug diff ',i,i-1, &
              (Rot_y_coords(i,1) - Rot_y_coords(i-1,2)), epsShift

         IF(DontMove(i) .AND. DontMove(i-1)) THEN
           CALL Warn(SolverName,&
                "Both nodes are marked Dont Move, stopping shifting.")
           EXIT
         END IF
         !Work out which to shift
         !If one is marked DontMove (corner node, or 2nd inland and weird shape)
         ! then move the other
         !Otherwise, shift whichever is further forward
         IF(DontMove(i)) THEN
           ShiftSecond = .FALSE.
           DontMove(i-1) = .TRUE.
           IF(Debug) PRINT *,ParEnv % MyPE,'Debug, dont move second: ',i
         ELSE IF(DontMove(i-1)) THEN
           ShiftSecond = .TRUE.
           DontMove(i) = .TRUE.
           IF(Debug) PRINT *,ParEnv % MyPE,'Debug, dont move first: ',i-1
         ELSE
           IF(SUM(Rot_z_coords(i,:)) > SUM(Rot_z_coords(i-1,:))) THEN
             ShiftSecond = .TRUE.
           ELSE
             ShiftSecond = .FALSE.
           END IF
         END IF

         IF(ShiftSecond) THEN
           ShiftIdx = i
           ShiftToIdx = i-1
         ELSE
           ShiftIdx = i-1
           ShiftToIdx = i
         END IF

         !Calculate ShiftTo, depends on ShiftSecond
         !Also update shifted Rot_y_coords
         IF(ShiftSecond) THEN
           !If this node has already been moved the other way, leave it alone
           IF((UpdatedDirection(ShiftIdx) < 0) .AND. .NOT. DontMove(ShiftToIdx)) CYCLE

           ShiftTo = Rot_y_coords(i-1,2) + epsShift
           IF(Rot_y_coords(i,2) < ShiftTo) UpdatedDirection(ShiftIdx) = 1

           Rot_y_coords(i,1) = MAX(Rot_y_coords(i,1), ShiftTo)
           Rot_y_coords(i,2) = MAX(Rot_y_coords(i,2), ShiftTo)
         ELSE
           !If this node has already been moved the other way, leave it alone
           IF((UpdatedDirection(ShiftIdx) > 0) .AND. .NOT. DontMove(ShiftToIdx)) CYCLE

           ShiftTo = Rot_y_coords(i,1) - epsShift
           IF(Rot_y_coords(i,1) > ShiftTo) UpdatedDirection(ShiftIdx) = -1

           Rot_y_coords(i-1,1) = MIN(Rot_y_coords(i-1,1), ShiftTo)
           Rot_y_coords(i-1,2) = MIN(Rot_y_coords(i-1,2), ShiftTo)
         END IF

         IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, updating col ',ShiftIdx,&
              ' rot_y_coords: ',Rot_y_coords(ShiftIdx,:)

         UpdatedColumn(ShiftIdx) = .TRUE.

         !Shift all the nodes in this column
         DO j=1,Mesh % NumberOfNodes
           IF(FNColumns(j) /= FrontColumnList(ShiftIdx)) CYCLE

           NodeHolder(1) = Mesh % Nodes % x(j) + Advance((Perm(j)-1)*DOFs + 1)
           NodeHolder(2) = Mesh % Nodes % y(j) + Advance((Perm(j)-1)*DOFs + 2)
           NodeHolder(3) = Mesh % Nodes % z(j) + Advance((Perm(j)-1)*DOFs + 3)
           NodeHolder = MATMUL(RotationMatrix, NodeHolder)

           Displace = 0.0_dp
           Displace(2) = ShiftTo - NodeHolder(2)

           !Check if node already projectable
           IF(ShiftSecond) THEN
             IF(Displace(2) < 0 ) CYCLE
           ELSE
             IF(Displace(2) > 0) CYCLE
           END IF

           Displace = MATMUL(TRANSPOSE(RotationMatrix), Displace)

           IF(Debug) PRINT *,ParEnv % MyPE, 'Debug, shifting node ',j,' col ',ShiftIdx,'xyz: ',&
                Mesh % Nodes % x(j)+Advance((Perm(j)-1)*DOFs + 1),&
                Mesh % Nodes % y(j)+Advance((Perm(j)-1)*DOFs + 2),&
                Mesh % Nodes % z(j)+Advance((Perm(j)-1)*DOFs + 3),&
                ' by ',Displace,' to ensure projectability. 1'


           !Adjust the variable values to shift nodes in line.
           Advance((Perm(j)-1)*DOFs + 1) = Advance((Perm(j)-1)*DOFs + 1) + Displace(1)
           Advance((Perm(j)-1)*DOFs + 2) = Advance((Perm(j)-1)*DOFs + 2) + Displace(2)

         END DO
       END IF
     END DO

     IF(ANY(UpdatedColumn)) MovedOne = .TRUE.

     IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, MovedOne: ',MovedOne
     IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, UpdatedColumn: ',UpdatedColumn

   END DO


   !-----------------------------------------------
   ! Check for pinnacles or rifts which can't be made
   ! projectable. Set the columns in a line and mark them
   ! for removal by Remesh.F90
   !-----------------------------------------------
   !TODO: At least 3 separate issues here:
   !
   !     *This one should be fixed now, at the end we check for groups disappearing*
   ! 1 : If there's a tangled group, then another, superset tangled
   !     group, this causes "programming error, tangled group has zero nodes"
   !     potential solution here is to just remove this error message, because
   !     its not really a problem *I DON'T THINK*. However, if so, need to
   !     check for it, and shift TangledGroups down one, so Remesh behaves
   !
   !     *This one should now be fixed, because we do all the unprojectable shifting first:*
   ! 2 : Here we iteratively: check tangled, then check unprojectable
   !     If correcting unprojectability causes further tangling, this
   !     isn't caught. Potential solutions:
   !            Improve parallel unprojectability checker (dirty fix, won't catch every poss case)
   !            Better: recheck previously tangled regions
   !
   ! 3 : Long straight sides which also happen to tangle with the end of their headlands
   !     are not well handled. All the nodes along these sides are shifted out the way
   !     in the new mesh


   NoTangledGroups = 0
   DO i=1,FrontLineCount-2
     IF(Tangled(i)) CYCLE
     j = i+2

     !If the column two columns away is less than 2*epsShift away, and there
     !is a change of direction, problems...
     IF((Rot_y_coords(j,1) - Rot_y_coords(i,2)) < 2*epsShift) THEN

       IF(((Rot_z_coords(i,1) < Rot_z_coords(i+1,1)) .NEQV. &
            (Rot_z_coords(i+1,1) < Rot_z_coords(j,1))) .OR. &
            ((Rot_z_coords(i,2) < Rot_z_coords(i+1,2)) .NEQV. &
            (Rot_z_coords(i+1,2) < Rot_z_coords(j,2)))) THEN

         !Either a pinnacle or a slit (i.e protrusion or rift)
         Protrusion = SUM(Rot_z_coords(i+1,:)) > SUM(Rot_z_coords(i,:))
         IF(Debug) PRINT *,'Debug, protrusion: ',Protrusion

         NoTangledGroups = NoTangledGroups + 1
         PivotIdx = i+1

         DO k=2,PivotIdx
           p1(1) = Rot_y_coords(k,2)
           p1(2) = Rot_z_coords(k,2)

           p2(1) = Rot_y_coords(k-1,2)
           p2(2) = Rot_z_coords(k-1,2)

           DO m=FrontLineCount-1,PivotIdx,-1
             IF(k==m) CYCLE !first two will always intersect by definition, not what we want
             q1(1) = Rot_y_coords(m,1)
             q1(2) = Rot_z_coords(m,1)

             q2(1) = Rot_y_coords(m+1,1)
             q2(2) = Rot_z_coords(m+1,1)

             CALL LineSegmentsIntersect ( p1, p2, q1, q2, intersect, intersect_flag )

             IF(intersect_flag) EXIT
           END DO
           IF(intersect_flag) EXIT
         END DO

         IF(intersect_flag) THEN

           !Found an intersection point (intersect)
           PRINT *,'Debug, found tangle intersection ',intersect,' leaving last tangled nodes: ',k,m

           Tangled(k:m) = .TRUE.
           TangledPivotIdx(k:m) = PivotIdx
           TangledShiftTo(k:m) = intersect(1)

           IF(Protrusion) THEN
             TangledGroup(k:m) = NoTangledGroups
           ELSE
             TangledGroup(k:m) = -NoTangledGroups
           END IF

         ELSE
           PRINT *,'Debug, found no intersection, so nodes arent QUITE tangled: ',i,j

           Tangled(i:j) = .TRUE.
           TangledPivotIdx(i:j) = PivotIdx
           TangledShiftTo(i:j) = (SUM(rot_y_coords(i,:)) + SUM(rot_y_coords(j,:))) / 4.0_dp

           IF(Protrusion) THEN
             TangledGroup(i:j) = NoTangledGroups
           ELSE
             TangledGroup(i:j) = -NoTangledGroups
           END IF

         END IF
       END IF
     END IF
   END DO

   !Check for a tangled group 'swallowing' another
   county = 0
   DO i=1,NoTangledGroups
     IF(COUNT(TangledGroup == i) > 0) CYCLE
     county = county + 1
     !Shift group numbers down one
     DO j=1,SIZE(TangledGroup)
       IF(TangledGroup(j) > i) TangledGroup(j) = TangledGroup(j) - 1
     END DO
   END DO
   NoTangledGroups = NoTangledGroups - county

   !Strategy:
   ! Shift all nodes to near (offset) the y-coordinate
   ! of the tangle pivot (i+1, above)
   DO i=1,NoTangledGroups
     n = COUNT(TangledGroup == i)
     IF(n==0) THEN
       n = COUNT(TangledGroup == -i)
       Protrusion = .FALSE.
     ELSE
       Protrusion = .TRUE.
     END IF
     IF(n==0) CALL Fatal(SolverName, "Programming error: tangled group has 0 nodes?")

     !Get the pivot node index, ShiftToY, and tangled node range
     FirstTangleIdx = 0
     LastTangleIdx = 0
     DO j=1,FrontLineCount
       IF(TangledGroup(j) == i .OR. TangledGroup(j) == -i) THEN
         IF(FirstTangleIdx == 0) FirstTangleIdx = j
         LastTangleIdx = j
         PivotIdx = TangledPivotIdx(j)
         ShiftToY = TangledShiftTo(j)
         !EXIT
       END IF
     END DO

     IF(LastTangleIdx - FirstTangleIdx /= n-1) CALL Fatal(SolverName, &
          "Programming error: wrong number of nodes in TangledGroup")
     IF(Debug) PRINT *,'Debug, tangled pivot index, ShiftToY: ', PivotIdx, ShiftToY

     LeftY =  ShiftToY - (epsTangle * ((n-1) / 2.0_dp))
     RightY = ShiftToY + (epsTangle * ((n-1) / 2.0_dp))

     !Potentially 'squeezed' on either side by neighbour columns
     SqueezeLeft = .FALSE.
     SqueezeRight = .FALSE.

     IF(LeftY < Rot_y_coords(FirstTangleIdx-1,2) + epsTangle) THEN
       SqueezeLeft = .TRUE.
       IF( ( (Rot_y_coords(LastTangleIdx+1,1) - epsTangle) - (epsTangle * (n-1)) ) < &
            (Rot_y_coords(FirstTangleIdx-1,2) + epsTangle)) THEN
         SqueezeRight = .TRUE.
       END IF
     END IF
     IF(RightY > (Rot_y_coords(LastTangleIdx+1,1) - epsTangle)) THEN
       SqueezeRight = .TRUE.
       IF( ( (Rot_y_coords(LastTangleIdx+1,1) - epsTangle) - (epsTangle * (n-1)) ) < &
            (Rot_y_coords(FirstTangleIdx-1,2) + epsTangle)) THEN
         SqueezeLeft = .TRUE.
       END IF
     END IF

     IF(Debug) PRINT *,'Debug, LeftY: ',LeftY,' RightY: ',RightY,' SqueezeL: ',&
          SqueezeLeft,' SqueezeR: ',SqueezeRight,' prev: ',Rot_y_coords(FirstTangleIdx-1,2),&
          ' next: ',Rot_y_coords(LastTangleIdx+1,1)

     !If squeezed, adjust the new y coord of the tangled nodes
     IF(SqueezeLeft .AND. SqueezeRight) THEN
       thisEps = (Rot_y_coords(LastTangleIdx+1,1) - Rot_y_coords(FirstTangleIdx-1,2)) / (n+1)
       LeftY = Rot_y_coords(FirstTangleIdx-1,2) + thisEps
       RightY = Rot_y_coords(LastTangleIdx+1,1) - thisEps
       ShiftToY = (RightY + LeftY) / 2.0_dp !midpoint
     ELSE IF(SqueezeLeft) THEN
       Shift = (Rot_y_coords(FirstTangleIdx-1,2) + epsTangle) - LeftY
       LeftY = LeftY + Shift
       RightY = RightY + Shift
       ShiftToY = ShiftToY + Shift
     ELSE IF(SqueezeRight) THEN
       Shift = (Rot_y_coords(LastTangleIdx+1,1) - epsTangle) - RightY
       LeftY = LeftY + Shift
       RightY = RightY + Shift
       ShiftToY = ShiftToY + Shift
     END IF

     !Where tangling occurs, we shift the tangled columns to be
     !1m apart in a series. Then Remesh gets rid of them
     DO j=FirstTangleIdx,LastTangleIdx
       IF(.NOT. (TangledGroup(j) == i .OR. TangledGroup(j) == -i)) &
            CALL Fatal(SolverName, "Programming error: node in specified idx range not tangled?")

       thisY = LeftY + ((REAL(j - FirstTangleIdx)/(n-1)) * (RightY - LeftY))
       IF(Debug) PRINT *,'Debug, thisY: ',thisY, j

       DO k=1,Mesh % NumberOfNodes
         IF(FNColumns(k) /= FrontColumnList(j)) CYCLE

         NodeHolder(1) = Mesh % Nodes % x(k) + Advance((Perm(k)-1)*DOFs + 1)
         NodeHolder(2) = Mesh % Nodes % y(k) + Advance((Perm(k)-1)*DOFs + 2)
         NodeHolder(3) = Mesh % Nodes % z(k) + Advance((Perm(k)-1)*DOFs + 3)
         IF(Debug) PRINT *, ParEnv % MyPE, ' Debug, pre shift tangled node ',k,': ',NodeHolder

         NodeHolder = MATMUL(RotationMatrix, NodeHolder)

         Displace = 0.0_dp
         Displace(2) = thisY - NodeHolder(2)

         Displace = MATMUL(TRANSPOSE(RotationMatrix), Displace)

         IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, shifting node ',k, ' col ',j,&
              ' by ',Displace,' to detangle.'

         !Adjust the variable values to shift nodes in line.
         Advance((Perm(k)-1)*DOFs + 1) = Advance((Perm(k)-1)*DOFs + 1) + Displace(1)
         Advance((Perm(k)-1)*DOFs + 2) = Advance((Perm(k)-1)*DOFs + 2) + Displace(2)
         !Set this variable so that Remesh knows
         TangledVar % Values(TangledVar % Perm(k)) = 1.0_dp * ABS(TangledGroup(j))
       END DO
     END DO
   END DO


   !---------------------------------------
   !Done, just deallocations

   FirstTime = .FALSE.

   DEALLOCATE(FrontPerm, &
        TopPerm, &
        FrontNodeNums, &
        Rot_y_coords, &
        UpdatedColumn, &
        UpdatedDirection, &
        Tangled, &
        TangledGroup, &
        TangledPivotIdx, &
        TangledShiftTo, &
        FrontColumnList, &
        FrontLocalNodeNumbers)

 END SUBROUTINE FrontAdvance3D
