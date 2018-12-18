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

!This moduled, loosely named 'CalvingGeometry' is for basically any
!reusable routines for the 3D calving model.

MODULE CalvingGeometry

  USE Types
  USE SParIterComm
  USE MainUtils

  IMPLICIT NONE

  INTERFACE DoubleIntVectorSize
     MODULE PROCEDURE DoubleIntVectorSizeP, DoubleIntVectorSizeA
  END INTERFACE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived type for 3D crevasse group info
  ! 
  ! Using the => Next, => Prev format like 
  ! variables_t, because there's no way of 
  ! knowing, a priori, how many we need.
  ! 
  ! Actually the only use of this is borrowed by BasalMelt3D.F90, so its misnamed...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE CrevasseGroup3D_t
     INTEGER :: NumberOfNodes = 0, ID = 0
     INTEGER, POINTER :: NodeNumbers(:) => NULL()
     INTEGER, POINTER :: BoundaryNodes(:) => NULL(), FrontNodes(:) => NULL() !allocatable too?
     REAL(KIND=dp) :: BoundingBox(4) !min_x, max_x, min_y, max_y
     
     LOGICAL :: FrontConnected !Does the group touch the terminus?
     TYPE(CrevasseGroup3D_t), POINTER :: Next => NULL(), Prev => NULL()
  END TYPE CrevasseGroup3D_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Derived type for a calving path defined by
  ! the IsoSurface/Line solver.
  ! (doubly linked list)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE CrevassePath_t
     INTEGER :: NumberOfNodes = 0, NumberOfElements = 0, ID = 0
     INTEGER, POINTER :: NodeNumbers(:) => NULL(), ElementNumbers(:)=>NULL()
!     INTEGER :: Ends(2)
     REAL(KIND=dp) :: Left, Right, Extent
     TYPE(CrevassePath_t), POINTER :: Next => NULL(), Prev => NULL()
     LOGICAL :: Valid = .TRUE.
  END TYPE  CrevassePath_t

CONTAINS
  

  !Returns the neighbours of a specified node using the matrix 
  !provided.
  !Note the current definition of neighbours:
  !Two nodes are neighbours if they are in the same bulk element
  !NOT ONLY if they are joined by a bar... 
  !User may provide an inverse perm (InvPerm_in), or else this will recomputed
  !each time (which would be pretty inefficient)
  FUNCTION FindNodeNeighbours(NodeNumber, Matrix, Perm, DOFs, InvPerm_in) RESULT (Neighbours)
    INTEGER :: NodeNumber, NoNeighbours, i, count, DOFs !<---!!!
    TYPE(Matrix_t), POINTER :: Matrix
    INTEGER, POINTER :: Perm(:), Neighbours(:), InvPerm(:)
    INTEGER, POINTER, OPTIONAL, INTENT(IN) :: InvPerm_in(:)
    LOGICAL :: Debug
    Debug = .FALSE.

    IF(PRESENT(InvPerm_in)) THEN
       InvPerm => InvPerm_in
    ELSE
       IF(Debug) PRINT *, 'Debug FindNodeNeighbours, creating InvPerm'
       InvPerm => CreateInvPerm(Perm)
    END IF

    NoNeighbours = Matrix % Rows((Perm(NodeNumber)*DOFs)+1) &
         - Matrix % Rows(Perm(NodeNumber)*DOFs)

    IF(MOD(NoNeighbours, DOFs).NE. 0) &
         CALL FATAL("Geometry","This shouldn't have happened...")

    !Each neighbour appears once per DOF, and there's also the current node thus: (x/DOFS) - 1...
    NoNeighbours = (NoNeighbours / DOFs) - 1

    ALLOCATE(Neighbours(NoNeighbours))
    Neighbours = 0

    count = 0

    DO i=Matrix % Rows(Perm(NodeNumber)*DOFs),&
         (Matrix % Rows((Perm(NodeNumber)*DOFs)+1)-1) !move along the row
       IF(MOD(i,DOFs) /= 0) CYCLE !Stored DOF1, DOF2, DOF3, only need every DOFth
       IF(MOD(Matrix % Cols(i), DOFs) /= 0) CALL Fatal("Geometry:FindNodeNeighbours", &
            "This is a programming error, Matrix structure is not what was expected.")

       IF(InvPerm(Matrix % Cols(i)/DOFs) == NodeNumber) CYCLE !Not our own neighbour
       count = count + 1
       Neighbours(count) = &
            InvPerm(Matrix % Cols(i)/DOFs)
    END DO

    IF(.NOT. PRESENT(InvPerm_in)) DEALLOCATE(InvPerm)

  END FUNCTION FindNodeNeighbours


  !-----------------------------------------------------------------------------
  !Returns the 2D (x,y) Cartesian distance between two nodes
  !NOTE: This isn't well programmed, should probably pass nodes...
  FUNCTION NodeDist2D(Nodes, NodeNum1, NodeNum2 ) RESULT (dist)
    TYPE(Nodes_t) :: Nodes
    INTEGER :: NodeNum1, NodeNum2
    REAL(KIND=dp) :: xdist,ydist,dist
    !Pythagoras in 2D
    xdist = Nodes % x(NodeNum1)&
         - Nodes % x(NodeNum2)
    ydist = Nodes % y(NodeNum1)&
         - Nodes % y(NodeNum2)
    !TODO: Can this be simplified?  See Interpolation.f90
    dist = ((xdist**2) + (ydist**2))**0.5
  END FUNCTION NodeDist2D

  !-----------------------------------------------------------------------------
  !Returns the 3D Cartesian distance between two nodes
  !NOTE: This isn't well programmed, should probably pass nodes...
  FUNCTION NodeDist3D( Nodes, Node1, Node2 ) RESULT (dist)
    TYPE(Nodes_t) :: Nodes
    INTEGER :: Node1, Node2
    REAL(KIND=dp) :: xdist,ydist,zdist,xydist,dist
    !Pythagoras in 3D
    xdist = Nodes % x(Node1)&
         - Nodes % x(Node2)
    ydist = Nodes % y(Node1)&
         - Nodes % y(Node2)
    zdist = Nodes % z(Node1)&
         - Nodes % z(Node2)
    !TODO: Can this be simplified?  See Interpolation.f90
    xydist = ((xdist**2) + (ydist**2))**0.5
    dist = ((xydist**2) + (zdist**2))**0.5
  END FUNCTION NodeDist3D

  !-----------------------------------------------------------------------------
  !Returns the inverse permutation table for a given perm and DOFs
  !NOTE, differs from the definition of InvPerm currently used in
  !Calving.F90
  FUNCTION CreateInvPerm(Perm) RESULT(InvPerm)
    INTEGER, POINTER :: Perm(:), InvPerm(:)
    INTEGER :: i, j

    ALLOCATE(InvPerm(MAXVAL(Perm)))

    j = 0
    DO i=1,SIZE(Perm)
       IF(Perm(i) == 0) CYCLE
       j = j + 1
       InvPerm( Perm(i) ) = j
    END DO

  END FUNCTION CreateInvPerm

  !-----------------------------------------------------------------------------
  !Returns dx/dy for two given nodes
  FUNCTION NodesGradXY(Nodes, Node1, Node2)RESULT(dxdy)
    INTEGER :: Node1, Node2
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: dx,dy,dxdy

    dx = Nodes % x(Node1) - Nodes % x(Node2)
    dy = Nodes % y(Node1) - Nodes % y(Node2)
    dxdy = dx/dy
  END FUNCTION NodesGradXY

  !-----------------------------------------------------------------------------
  !Returns the number of decimal places of a real number
  !which has been read from a text file (.e.g mesh.nodes)
  !this differs from intrinsic PRECISION() because these
  !numbers often have trailing 000s or 999s
  FUNCTION RealAeps(in)RESULT(myaeps)
    REAL(KIND=dp) :: in, toler, x, myaeps
    INTEGER :: sigs, mag, decs

    !Find how many decimal places
    mag = FLOOR(LOG10(ABS(in))) + 1 !Order of magnitude of number
    decs = PRECISION(in) - mag  !total digits - magnitude = decimal places

    toler = 10.0_dp**(-decs)
    sigs = 0
    x = in

    DO WHILE (.TRUE.)
       IF(ABS(x - NINT(x)) < toler) THEN !found the precision limit
          EXIT
       ELSE
          sigs = sigs + 1
          x = x * 10 !move the decimal point along
          x = x - FLOOR(x) !convert number to O(1) so FLOOR doesn't reach integer limit
          toler = toler * 10.0_dp !1 fewer remaining decimal places
       END IF
    END DO
    myaeps = 10.0**(-sigs)
  END FUNCTION RealAeps

  !-----------------------------------------------------------------------------
  ! Constructs paths of connected isoline (202) elements which intersect the 
  ! front. Each path will begin and end with a node where OnFront=.TRUE.
  !-----------------------------------------------------------------------------
  SUBROUTINE FindCrevassePaths(IsoMesh, OnFront, CrevassePaths, PathCount)
    IMPLICIT NONE
    TYPE(Mesh_t), POINTER :: IsoMesh
    LOGICAL, ALLOCATABLE :: OnFront(:)
    TYPE(CrevassePath_t), POINTER :: CrevassePaths
    INTEGER :: PathCount
    !----------------------------------------------
    TYPE(CrevassePath_t), POINTER :: CurrentPath
    LOGICAL :: Found, Debug
    INTEGER :: i,j,NodeCount,ElemCount, NextElem
    INTEGER, ALLOCATABLE :: WorkElems(:), WorkNodes(:)

    Debug = .FALSE.
    PathCount = 1

    !TODO assert all 202 elements

    ALLOCATE(CrevassePaths)
    CurrentPath => CrevassePaths

    ALLOCATE(WorkElems(100), WorkNodes(100))
    WorkElems = 0; WorkNodes = 0

    DO i=1, IsoMesh % NumberOfBulkElements

       IF(ANY(OnFront(Isomesh % Elements(i) % NodeIndexes))) THEN
          !Found an element with one node on calving front

          IF(ElementPathID(CrevassePaths, i) /= 0) CYCLE !already in a path

          !Starting a new group...
          CurrentPath % ID = PathCount
          IF(Debug) PRINT *, 'Potential calving isomesh element: ',i

          ElemCount = 1
          NextElem = i

          !Identify which of the two nodes are on the front...
          DO j=1,2
             IF(OnFront(IsoMesh % Elements(i) % NodeIndexes(j))) EXIT
          END DO
          IF(j==3) CALL Fatal("FindCrevassePaths", "Couldn't find node on boundary")

          !... and put it first in the list
          WorkNodes(1) = IsoMesh % Elements(i) % NodeIndexes(j)
          NodeCount = 2

          !Follow the chain
          DO WHILE(.TRUE.) 

             WorkElems(ElemCount) = NextElem
             ElemCount = ElemCount + 1
             !Put the other node into the list
             DO j=1,2
                IF(ANY(WorkNodes == IsoMesh % Elements(NextElem) % NodeIndexes(j))) CYCLE
                WorkNodes(NodeCount) = IsoMesh % Elements(NextElem) % NodeIndexes(j)
                NodeCount = NodeCount + 1
                EXIT
             END DO

             !Look for element which contains previous element's node
             Found = .FALSE.
             DO j=1,IsoMesh % NumberOfBulkElements
                IF(ANY(IsoMesh % Elements(j) % NodeIndexes == WorkNodes(NodeCount-1))) THEN

                   !already in another group (is this possible?)
                   IF(ElementPathID(CrevassePaths, j ) /= 0) CYCLE 
                   !Already in current group
                   IF(ANY(WorkElems == j)) CYCLE 
                   
                   NextElem = j
                   Found = .TRUE.
                   EXIT
                END IF
             END DO

             IF(.NOT. Found) EXIT

             IF(ElemCount > SIZE(WorkElems)) THEN
                IF(Debug) PRINT *, 'FindCrevassePaths, doubling size of element array.'
                CALL DoubleIntVectorSize(WorkElems)
             END IF
             IF(NodeCount > SIZE(WorkNodes)) THEN
                IF(Debug) PRINT *, 'FindCrevassePaths, doubling size of node array.'
                CALL DoubleIntVectorSize(WorkNodes)
             END IF
          END DO
          
          ElemCount = ElemCount - 1
          NodeCount = NodeCount - 1

          CurrentPath % NumberOfNodes = NodeCount
          CurrentPath % NumberOfElements = ElemCount

          ALLOCATE(CurrentPath % ElementNumbers(ElemCount), &
               CurrentPath % NodeNumbers(NodeCount))
          
          CurrentPath % NodeNumbers = WorkNodes(1:NodeCount)
          CurrentPath % ElementNumbers = WorkElems(1:ElemCount)

          WorkNodes = 0
          WorkElems = 0

          ALLOCATE(CurrentPath % Next)
          CurrentPath % Next % Prev => CurrentPath
          CurrentPath => CurrentPath % Next
          PathCount = PathCount + 1
       END IF
    END DO

    !We always overshoot by one
    PathCount = PathCount - 1

    IF(PathCount > 0) THEN
       PRINT *,'Number of crevasse paths: ', PathCount
       CurrentPath % Prev % Next => NULL()
       DEALLOCATE(CurrentPath)
    ELSE
       PRINT *,'No crevasse paths'
       DEALLOCATE(CrevassePaths)
    END IF

    DEALLOCATE(WorkNodes, WorkElems)

  END SUBROUTINE FindCrevassePaths

  !Removes a CrevassePath from a linked list of CrevassePaths
  SUBROUTINE RemoveCrevassePath(Path)
    IMPLICIT NONE
    TYPE(CrevassePath_t), POINTER :: Path
    !------------------------------------------------
    IF(ASSOCIATED(Path % Prev)) Path % Prev % Next => Path % Next
    IF(ASSOCIATED(Path % Next)) Path % Next % Prev => Path % Prev

    IF(ASSOCIATED(Path % NodeNumbers)) DEALLOCATE(Path % NodeNumbers)
    IF(ASSOCIATED(Path % ElementNumbers)) DEALLOCATE(Path % ElementNumbers)
    DEALLOCATE(Path)

  END SUBROUTINE RemoveCrevassePath
  
  !--------------------------------------------------------------------
  ! 'Tidies up' isomesh and the CrevassePaths found by FindCrevassePaths
  !--------------------------------------------------------------------
  ! This involves removing duplicate nodes, taking care to replace node
  ! indexes in affected elements. This then allows easy removal of 
  ! 202 elements with zero length.
  !
  ! Closed loops are removed from crevasse paths
  !--------------------------------------------------------------------

  SUBROUTINE CheckCrevasseNodes(Mesh, CrevassePaths)
    IMPLICIT NONE
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(CrevassePath_t), POINTER :: CrevassePaths
    !-------------------------------------------------
    TYPE(CrevassePath_t), POINTER :: CurrentPath,WorkPath
    INTEGER :: i,j,ElNo,counter, ElementNumbers(2)
    INTEGER, ALLOCATABLE :: ReplaceWithNode(:),WorkInt(:)
    LOGICAL :: Debug
    LOGICAL, ALLOCATABLE :: RemoveElement(:), RemoveNode(:), PathRemoveElement(:)

    Debug = .FALSE.

    ALLOCATE(RemoveNode(Mesh % NumberOfNodes),&
         ReplaceWithNode(Mesh % NumberOfNodes),&
         RemoveElement(Mesh % NumberOfBulkElements))
    RemoveNode = .FALSE.
    RemoveElement = .FALSE.
    ReplaceWithNode = 0

    !Cycle mesh NODES, looking for duplicates and marking them for deletion
    DO i=1,Mesh % NumberOfNodes
       IF(RemoveNode(i)) CYCLE
       DO j=1,Mesh % NumberOfNodes
          IF(i==j .OR. RemoveNode(j)) CYCLE
          IF(Mesh % Nodes % x(i) == Mesh % Nodes % x(j) .AND.&
               Mesh % Nodes % y(i) == Mesh % Nodes % y(j) .AND.&
               Mesh % Nodes % z(i) == Mesh % Nodes % z(j)) THEN
             RemoveNode(j) = .TRUE.
             ReplaceWithNode(j) = i
          END IF
       END DO
    END DO

    !Replace element nodeindexes where nodes are removed
    DO i=1,Mesh % NumberOfBulkElements
       DO j=1,SIZE(Mesh % Elements(i) % NodeIndexes)
          IF(RemoveNode(Mesh % Elements(i) % NodeIndexes(j))) &
               Mesh % Elements(i) % NodeIndexes(j) = &
               ReplaceWithNode(Mesh % Elements(i) % NodeIndexes(j))
       END DO
    END DO

    !Mark elements with zero length (duplicate node indexes) for removal
    DO i=1,Mesh % NumberOfBulkElements
       IF(Mesh % Elements(i) % NodeIndexes(1) == Mesh % Elements(i) % NodeIndexes(2)) THEN
          RemoveElement(i) = .TRUE.
          IF(Debug) PRINT *,'debug, removing element: ',i,' with identical nodes: ',&
               Mesh % Elements(i) % NodeIndexes(1)
       END IF
    END DO

    IF(Debug) PRINT *,'Debug, removing ',COUNT(RemoveElement),' of ',SIZE(RemoveElement),' elements'

    !Cycle paths, looking for nodes which are identical and removing them, joining up elements etc
    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))

       IF(Debug) PRINT *,'Debug, Path: ',CurrentPath % ID,'initial no elems: ',&
            CurrentPath % NumberOfElements,&
            ' no nodes: ', CurrentPath % NumberOfNodes

       ALLOCATE(WorkInt(CurrentPath % NumberOfElements))
       WorkInt = 0
       counter = 0

       !Mark pairs of duplicate elements in path for removal
       ALLOCATE(PathRemoveElement(CurrentPath % NumberOfElements))
       PathRemoveElement = .FALSE.

       IF(CurrentPath % NumberOfElements == 1) THEN
         !Only has one element, remove
         PathRemoveElement = .TRUE.
       ELSE
         DO i=1,CurrentPath % NumberOfElements-1

           IF(PathRemoveElement(i)) CYCLE
           ElementNumbers(1) = CurrentPath % ElementNumbers(i)
           IF(RemoveElement(ElementNumbers(1))) CYCLE

           j = i+1
           IF(PathRemoveElement(j)) CYCLE
           ElementNumbers(2) = CurrentPath % ElementNumbers(j)
           IF(RemoveElement(ElementNumbers(2))) CYCLE

           IF( ANY(Mesh % Elements(ElementNumbers(1)) % NodeIndexes == &
                Mesh % Elements(ElementNumbers(2)) % NodeIndexes(1)) .AND. &
                ANY(Mesh % Elements(ElementNumbers(1)) % NodeIndexes == &
                Mesh % Elements(ElementNumbers(2)) % NodeIndexes(2)) ) THEN
             PathRemoveElement(j) = .TRUE.
             PathRemoveElement(i) = .TRUE.
             IF(Debug) PRINT *,'Path: ',CurrentPath % ID,' removing identical elements: ',i,' ',j
           END IF

         END DO


         !Check if entire crevasse group is a closed loop
         ElementNumbers(1) = CurrentPath % ElementNumbers(1)
         ElementNumbers(2) = CurrentPath % ElementNumbers(CurrentPath % NumberOfElements)
         DO i=1,2
           IF(.NOT. ANY(Mesh % Elements(CurrentPath % ElementNumbers(2)) % NodeIndexes == &
                Mesh % Elements(ElementNumbers(1)) % NodeIndexes(i))) EXIT
         END DO
         IF(i==3) CALL Fatal("CheckCrevassePaths","Programming error: unable to determine first node")
         IF(ANY(Mesh % Elements(ElementNumbers(2)) % NodeIndexes == &
              Mesh % Elements(ElementNumbers(1)) % NodeIndexes(i))) THEN
           PathRemoveElement = .TRUE.
           IF(Debug) PRINT *,'Debug, removing path ',CurrentPath % ID,' because its entirely closed.'
         END IF

         !For each element 'i' in turn, cycle backwards through element list looking
         !for element(i)'s nodes. If found, this indicates a closed loop which should
         !be removed.
         DO i=1,CurrentPath % NumberOfElements
           IF(PathRemoveElement(i)) CYCLE
           IF(RemoveElement(CurrentPath % ElementNumbers(i))) CYCLE
           ElementNumbers(1) = CurrentPath % ElementNumbers(i)

           DO j=CurrentPath % NumberOfElements,i+1,-1 !cycle backwards from end to i+1
             IF(PathRemoveElement(j)) CYCLE
             IF(RemoveElement(CurrentPath % ElementNumbers(j))) CYCLE
             ElementNumbers(2) = CurrentPath % ElementNumbers(j)

             IF( ANY(Mesh % Elements(ElementNumbers(1)) % NodeIndexes == &
                  Mesh % Elements(ElementNumbers(2)) % NodeIndexes(1)) .OR. &
                  ANY(Mesh % Elements(ElementNumbers(1)) % NodeIndexes == &
                  Mesh % Elements(ElementNumbers(2)) % NodeIndexes(2)) ) THEN
               PathRemoveElement(i+1:j-1) = .TRUE.
               IF(Debug) PRINT *,'CheckCrevasseNodes, &
                    &Removing a closed loop from ',i+1,' to ',j-1
             END IF

           END DO

         END DO
       END IF

       !Replace CrevassePath % ElementNumbers based on previous removals
       DO i=1,CurrentPath % NumberOfElements
          IF(.NOT.RemoveElement(CurrentPath % ElementNumbers(i)) .AND. &
               .NOT.PathRemoveElement(i)) THEN
             counter = counter + 1
             WorkInt(counter) = CurrentPath % ElementNumbers(i)
             IF(Debug) THEN
                PRINT *,'Debug, keeping element: ',i,' from path: ',CurrentPath % ID
                PRINT *,'Debug, element global: ',CurrentPath % ElementNumbers(i),' and nodes :',&
                     Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes
             END IF
          ELSE
             IF(Debug) THEN
                PRINT *,'Debug, removing element: ',i,' from path: ',CurrentPath % ID
                PRINT *,'Debug, element global: ',CurrentPath % ElementNumbers(i),' and nodes :',&
                     Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes
             END IF
          END IF
       END DO
       IF(counter < CurrentPath % NumberOfElements) THEN
          IF(Debug) PRINT *,'debug, path loses ',CurrentPath % NumberOfElements - counter,&
               ' of ',CurrentPath % NumberOfElements,' elements.'

          CurrentPath % NumberOfElements = counter
          DEALLOCATE(CurrentPath % ElementNumbers)
          ALLOCATE(CurrentPath % ElementNumbers(counter))

          CurrentPath % ElementNumbers = WorkInt(1:counter)
       END IF
       DEALLOCATE(WorkInt,PathRemoveElement)

       IF (CurrentPath % NumberOfElements <= 0) THEN
          WorkPath => CurrentPath % Next

          IF(ASSOCIATED(CurrentPath,CrevassePaths)) CrevassePaths => WorkPath
          CALL RemoveCrevassePath(CurrentPath)
          IF(Debug) CALL Info("CheckCrevasseNodes",&
               "Removing a crevasse path with no elements")
          CurrentPath => WorkPath
          CYCLE
       END IF

       !Now reconstruct node list for path:
       DEALLOCATE(CurrentPath % NodeNumbers)
       CurrentPath % NumberOfNodes = CurrentPath % NumberOfElements + 1
       ALLOCATE(CurrentPath % NodeNumbers(CurrentPath % NumberOfNodes))
       CurrentPath % NodeNumbers = 0

       !First node
       IF(CurrentPath % NumberOfElements >= 2) THEN
          DO i=1,2
             IF( ANY(Mesh % Elements(CurrentPath % ElementNumbers(2)) % NodeIndexes == &
                  Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(i))) CYCLE
             CurrentPath % NodeNumbers(1) = &
                  Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(i)

             IF(i==2) THEN !Reorder so that nodeindexes(1) and (2) are in chain order
               Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(2) = &
                    Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(1)
               Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(1) = &
                    CurrentPath % NodeNumbers(1)
             END IF
             EXIT
          END DO
       ELSE !Rare, single element path, choice of first node is arbitrary...
             CurrentPath % NodeNumbers(1) = &
                  Mesh % Elements(CurrentPath % ElementNumbers(1)) % NodeIndexes(1)
       END IF

       IF(Debug) PRINT *,'Path ',CurrentPath % ID,' has first node: ',CurrentPath % NodeNumbers(1)

       !Follow the chain...
       DO i=1,CurrentPath % NumberOfElements
          ElNo = CurrentPath % ElementNumbers(i)
          DO j=1,2
             IF(ANY(CurrentPath % NodeNumbers == Mesh % Elements(ElNo) % NodeIndexes(j))) CYCLE
             CurrentPath % NodeNumbers(i+1) = Mesh % Elements(ElNo) % NodeIndexes(j)

             IF(j==1) THEN !Reorder so that nodeindexes(1) and (2) are in chain order
               Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes(1) = &
                    Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes(2)
               Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes(2) = &
                    CurrentPath % NodeNumbers(i+1)
             END IF

             EXIT
          END DO
       END DO

       IF(Debug) PRINT *,'Debug, path ',CurrentPath % ID,' has nodes: ',CurrentPath % NodeNumbers
       IF(ANY(CurrentPath % NodeNumbers == 0)) CALL Fatal("CheckCrevasseNodes","Failed to fill node indexes")
       CurrentPath => CurrentPath % Next
    END DO

  END SUBROUTINE CheckCrevasseNodes

  !----------------------------------------------------
  ! Checks paths for projectability and overlap
  ! In case of overlap, smaller enclosed path is deleted
  ! In case of unprojectability, nodes are moved laterally
  ! to restore projectability.
  !----------------------------------------------------
  ! NOTE: if this breaks, it could be due to two paths
  !       sharing a node. Thinking about it, I see no reason
  !       this should be an issue, but we'll see...
  SUBROUTINE ValidateCrevassePaths(Mesh, CrevassePaths, FrontOrientation, PathCount, ValidPathCount)
    IMPLICIT NONE
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(CrevassePath_t), POINTER :: CrevassePaths
    REAL(KIND=dp) :: FrontOrientation(3)
    INTEGER :: PathCount, ValidPathCount, First, Last, LeftIdx, RightIdx
    !---------------------------------------------------
    REAL(KIND=dp) :: RotationMatrix(3,3), UnRotationMatrix(3,3), FrontDist, MaxDist, &
         ShiftTo, Dir1(2), Dir2(2)
    REAL(KIND=dp), ALLOCATABLE :: ConstrictDirection(:,:)
    TYPE(CrevassePath_t), POINTER :: CurrentPath, OtherPath, WorkPath, LeftPath, RightPath
    INTEGER :: i,j,k,n,ElNo,ShiftToMe, NodeNums(2),A,B,FirstIndex, LastIndex,Start
    INTEGER, ALLOCATABLE :: WorkInt(:)
    LOGICAL :: Debug, Shifted, CCW, ToLeft, Snakey, OtherRight, ShiftRightPath
    LOGICAL, ALLOCATABLE :: PathMoveNode(:), DeleteElement(:), BreakElement(:), &
         FarNode(:), DeleteNode(:), Constriction(:)

    Debug = .FALSE.
    Snakey = .TRUE.

    RotationMatrix = ComputeRotationMatrix(FrontOrientation)
    UnRotationMatrix = TRANSPOSE(RotationMatrix)

    ! Temporarily rotate the mesh
    CALL RotateMesh(Mesh, RotationMatrix)

    ! Find path %left, %right, %extent (width)
    CALL ComputePathExtent(CrevassePaths, Mesh % Nodes, .TRUE.)

    IF(Snakey) THEN
      !-----------------------------------------------------
      ! Paths should not 'snake' inwards in a narrow slit...
      !-----------------------------------------------------

      !it's insufficent to require that no nodes be
      !further away than the two edge nodes.
      !Instead, must ensure that no nodes are further away than any
      !surrounding nodes.

      !First need to determine path orientation
      !with respect to front....

      CurrentPath => CrevassePaths
      DO WHILE(ASSOCIATED(CurrentPath))

        !First and last node on path
        First = CurrentPath % NodeNumbers(1)
        Last = CurrentPath % NodeNumbers(CurrentPath % NumberOfNodes)

        !if ToLeft, the crevasse path goes from right to left, from the
        !perspective of someone sitting in the fjord, looking at the front
        ToLeft = Mesh % Nodes % y(Last) > Mesh % Nodes % y(First)

        IF(Debug) THEN
          FrontDist = NodeDist3D(Mesh % Nodes,First, Last)
          PRINT *,'PATH: ', CurrentPath % ID, ' FrontDist: ',FrontDist
          PRINT *,'PATH: ', CurrentPath % ID, &
               ' nonodes: ',CurrentPath % NumberOfNodes,&
               ' noelems: ',CurrentPath % NumberOfElements
        END IF

        !Cycle path nodes, finding those which are too far away
        ALLOCATE(FarNode(CurrentPath % NumberOfNodes), &
             Constriction(CurrentPath % NumberOfNodes),&
             ConstrictDirection(CurrentPath % NumberOfNodes,2))
        FarNode = .FALSE.
        Constriction = .FALSE.
        ConstrictDirection = 0.0_dp

        !Determine which nodes have the potential to be constriction (based on angle)
        !and compute constriction direction (i.e. which way the 'pointy bit' points...')
        DO i=2,CurrentPath % NumberOfNodes-1
          First = CurrentPath % NodeNumbers(i-1)
          Last = CurrentPath % NodeNumbers(i+1)
          n = CurrentPath % NodeNumbers(i)

          CCW = ((Mesh % Nodes % y(n) - Mesh % Nodes % y(First)) * &
               (Mesh % Nodes % z(Last) - Mesh % Nodes % z(First))) > &
               ((Mesh % Nodes % z(n) - Mesh % Nodes % z(First)) * &
               (Mesh % Nodes % y(Last) - Mesh % Nodes % y(First)))

          IF(CCW .NEQV. ToLeft) THEN
            Constriction(i) = .TRUE.
            !Calculate constriction direction

            Dir1(1) = Mesh % Nodes % y(n) - Mesh % Nodes % y(First)
            Dir1(2) = Mesh % Nodes % z(n) - Mesh % Nodes % z(First)
            Dir1 = Dir1 / ((Dir1(1)**2.0 + Dir1(2)**2.0) ** 0.5)

            Dir2(1) = Mesh % Nodes % y(n) - Mesh % Nodes % y(Last)
            Dir2(2) = Mesh % Nodes % z(n) - Mesh % Nodes % z(Last)
            Dir2 = Dir2 / ((Dir2(1)**2.0 + Dir2(2)**2.0) ** 0.5)

            ConstrictDirection(i,1) = Dir1(1) + Dir2(1)
            ConstrictDirection(i,2) = Dir1(2) + Dir2(2)
            ConstrictDirection(i,:) = ConstrictDirection(i,:) / &
                 ((ConstrictDirection(i,1)**2.0 + ConstrictDirection(i,2)**2.0) ** 0.5)

            IF(Debug) PRINT *, 'Debug, node ',i,' dir1,2: ',Dir1, Dir2
            IF(Debug) PRINT *, 'Debug, node ',i,' constriction direction: ',ConstrictDirection(i,:)
          END IF
        END DO

        !First and last can always be constriction
        Constriction(1) = .TRUE.
        Constriction(SIZE(Constriction)) = .TRUE.

        !Compute constriction direction for first and last
        !We don't have info about the third node, so take orthogonal to 2
        Last = CurrentPath % NodeNumbers(2)
        n = CurrentPath % NodeNumbers(1)
        Dir1(1) = Mesh % Nodes % y(n) - Mesh % Nodes % y(Last)
        Dir1(2) = Mesh % Nodes % z(n) - Mesh % Nodes % z(Last)
        Dir1 = Dir1 / ((Dir1(1)**2.0 + Dir1(2)**2.0) ** 0.5)

        !Depending on which end of the chain we are,
        !we take either the right or left orthogonal vector
        IF(ToLeft) THEN
          ConstrictDirection(1,1) = Dir1(2)
          ConstrictDirection(1,2) = -1.0 * Dir1(1)
        ELSE
          ConstrictDirection(1,1) = -1.0 * Dir1(2)
          ConstrictDirection(1,2) = Dir1(1)
        END IF
        IF(Debug) PRINT *, 'Debug, node 1 constriction direction: ',ConstrictDirection(1,:)

        Last = CurrentPath % NodeNumbers(CurrentPath % NumberOfNodes - 1)
        n = CurrentPath % NodeNumbers(CurrentPath % NumberOfNodes)

        Dir1(1) = Mesh % Nodes % y(n) - Mesh % Nodes % y(Last)
        Dir1(2) = Mesh % Nodes % z(n) - Mesh % Nodes % z(Last)
        Dir1 = Dir1 / ((Dir1(1)**2.0 + Dir1(2)**2.0) ** 0.5)

        IF(.NOT. ToLeft) THEN
          ConstrictDirection(CurrentPath % NumberOfNodes,1) = Dir1(2)
          ConstrictDirection(CurrentPath % NumberOfNodes,2) = -1.0 * Dir1(1)
        ELSE
          ConstrictDirection(CurrentPath % NumberOfNodes,1) = -1.0 * Dir1(2)
          ConstrictDirection(CurrentPath % NumberOfNodes,2) = Dir1(1)
        END IF
        IF(Debug) PRINT *, 'Debug, node last constriction direction: ',&
             ConstrictDirection(CurrentPath % NumberOfNodes,:)

        !---------------------------------------
        ! Now that we have constrictions marked and directions computed, cycle nodes

        DO i=1,CurrentPath % NumberOfNodes
          IF(.NOT. Constriction(i)) CYCLE

          DO j=CurrentPath % NumberOfNodes,i+1,-1
            IF(.NOT. Constriction(j)) CYCLE


            First = CurrentPath % NodeNumbers(i)
            Last = CurrentPath % NodeNumbers(j)

            !Check that these constrictions 'face' each other via dot product
            Dir1(1) = Mesh % Nodes % y(Last) - Mesh % Nodes % y(First)
            Dir1(2) = Mesh % Nodes % z(Last) - Mesh % Nodes % z(First)
            Dir2(1) = -Dir1(1)
            Dir2(2) = -Dir1(2)

            !If the two constrictions aren't roughly facing each other:
            ! <  > rather than    > <
            ! then skip this combo
            IF(SUM(ConstrictDirection(i,:)*Dir1) < 0) THEN
              IF(Debug) PRINT *,'Constrictions ',i,j,' dont face each other 1: ',&
                   SUM(ConstrictDirection(i,:)*Dir1)
              CYCLE
            END IF

            IF(SUM(ConstrictDirection(j,:)*Dir2) < 0) THEN
              IF(Debug) PRINT *,'Constrictions ',j,i,' dont face each other 2: ',&
                   SUM(ConstrictDirection(j,:)*Dir2)
              CYCLE
            END IF

            MaxDist = NodeDist3D(Mesh % Nodes,First, Last)

            DO k=i+1,j-1
              IF(FarNode(k)) CYCLE

              n = CurrentPath % NodeNumbers(k)

              IF((NodeDist3D(Mesh % Nodes, First, n) <= MaxDist) .AND. &
                   (NodeDist3D(Mesh % Nodes, Last, n) <= MaxDist)) CYCLE !within range

              FarNode(k) = .TRUE.

            END DO
          END DO
        END DO

        !Cycle elements, marking those which need to be adjusted
        ALLOCATE(BreakElement(CurrentPath % NumberOfElements),&
             DeleteElement(CurrentPath % NumberOfElements))
        BreakElement = .FALSE.
        DeleteElement = .FALSE.

        DO i=1,CurrentPath % NumberOfElements
          IF(ANY(FarNode(i:i+1))) THEN
            IF(ALL(FarNode(i:i+1))) THEN
              DeleteElement(i) = .TRUE.
              IF(Debug) PRINT  *,'PATH: ', CurrentPath % ID, ' element: ',i,' is deleted.'
            ELSE
              BreakElement(i) = .TRUE.
              IF(Debug) PRINT  *,'PATH: ', CurrentPath % ID, ' element: ',i,' is broken.'
            END IF
          END IF
        END DO

        DO i=1,CurrentPath % NumberOfElements
          IF((.NOT. BreakElement(i)) .OR. DeleteElement(i)) CYCLE
          !This element needs to be adjusted
          DO j=i+1,CurrentPath % NumberOfElements
            IF(.NOT. (BreakElement(j) .OR. DeleteElement(j))) &
                 CALL Fatal("ValidateCrevasseGroups","Programming error in maintaining aspect ratio")
            IF(DeleteElement(j)) CYCLE
            !This is the next 'break element' after i
            !Determine which nodes we keep

            IF((CurrentPath % NodeNumbers(j) /= &
                 Mesh % Elements(CurrentPath % ElementNumbers(j)) % NodeIndexes(1)) .OR. &
                 (CurrentPath % NodeNumbers(j+1) /= &
                 Mesh % Elements(CurrentPath % ElementNumbers(j)) % NodeIndexes(2))) THEN

              CALL Fatal("ValidateCrevassePaths", "Chain building error")
            END IF

            Mesh % Elements(CurrentPath % ElementNumbers(i)) % NodeIndexes(2) = &
                 Mesh % Elements(CurrentPath % ElementNumbers(j)) % NodeIndexes(2)

            !We now want to delete it, because we only keep one from each broken pair
            DeleteElement(j) = .TRUE.
            EXIT !we paired this one, move on
          END DO
        END DO

        !Delete the elements and nodes
        IF(COUNT(DeleteElement) > 0) THEN
          !elements
          ALLOCATE(WorkInt(COUNT(.NOT. DeleteElement)))
          WorkInt = PACK(CurrentPath % ElementNumbers,.NOT.DeleteElement)

          DEALLOCATE(CurrentPath % ElementNumbers)
          ALLOCATE(CurrentPath % ElementNumbers(SIZE(WorkInt)))

          CurrentPath % ElementNumbers = WorkInt
          CurrentPath % NumberOfElements = SIZE(WorkInt)
          DEALLOCATE(WorkInt)

          !nodes
          ALLOCATE(WorkInt(COUNT(.NOT. FarNode)))
          WorkInt = PACK(CurrentPath % NodeNumbers, .NOT.FarNode)

          DEALLOCATE(CurrentPath % NodeNumbers)
          ALLOCATE(CurrentPath % NodeNumbers(SIZE(WorkInt)))

          CurrentPath % NodeNumbers = WorkInt
          CurrentPath % NumberOfNodes = SIZE(WorkInt)
          DEALLOCATE(WorkInt)
        END IF

        DEALLOCATE(FarNode, Constriction, ConstrictDirection, BreakElement, DeleteElement)
        CurrentPath => CurrentPath % Next
      END DO
    END IF !Snakey

    !Update Left, Right & Extent
    CALL ComputePathExtent(CrevassePaths, Mesh % Nodes, .TRUE.)

    !-----------------------------------------------------
    ! Move nodes from crevassepaths which aren't projectable
    !-----------------------------------------------------
    !  1) Path elements are ordered as a chain
    !  2) Path % Element(i) has nodes i, i+1
    !
    !  Go through CrevassePath nodes, marking those
    !  which are 'shadowed' by further away elements.
    !-----------------------------------------------------

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))

      ALLOCATE(PathMoveNode(CurrentPath % NumberOfNodes))
      PathMoveNode = .FALSE.

      DO i=1,CurrentPath % NumberOfNodes
        n = CurrentPath % NodeNumbers(i)
        DO j=1,CurrentPath % NumberOfElements
          ElNo = CurrentPath % ElementNumbers(j)
          NodeNums = Mesh % Elements(ElNo) % NodeIndexes
          IF(ANY(NodeNums == n)) CYCLE !Node is in element, skip
          !Check if node lies between element nodes
          IF( (Mesh % Nodes % y(NodeNums(1)) > Mesh % Nodes % y(n)) .NEQV. &
               (Mesh % Nodes % y(NodeNums(2)) > Mesh % Nodes % y(n)) ) THEN
            !Check the node is in front of the element

            A = MINLOC(Mesh % Nodes % z(NodeNums),1)
            B = MAXLOC(Mesh % Nodes % z(NodeNums),1)
            CCW = ((Mesh % Nodes % y(n) - Mesh % Nodes % y(NodeNums(A))) * &
                 (Mesh % Nodes % z(NodeNums(B)) - Mesh % Nodes % z(NodeNums(A)))) > &
                 ((Mesh % Nodes % z(n) - Mesh % Nodes % z(NodeNums(A))) * &
                 (Mesh % Nodes % y(NodeNums(B)) - Mesh % Nodes % y(NodeNums(A))))

            ToLeft = Mesh % Nodes % y(NodeNums(A)) > Mesh % Nodes % y(NodeNums(B))

            IF(CCW .EQV. ToLeft) THEN
              !Node should be removed
              PathMoveNode(i) = .TRUE.
              EXIT
            END IF

          END IF
        END DO
      END DO

      IF(Debug) THEN
        PRINT *,'Path ',CurrentPath % ID,' has ',&
             COUNT(PathMoveNode),' nodes which need to be shifted.'

        DO i=1,CurrentPath % NumberOfNodes
          IF(.NOT. PathMoveNode(i)) CYCLE
          PRINT *,'Need to move node: ',i,' y: ',&
               Mesh % Nodes % y(CurrentPath % NodeNumbers(i)),&
               ' z: ',Mesh % Nodes % z(CurrentPath % NodeNumbers(i))

        END DO
      END IF

      !Now that nodes have been marked as shadowed, identify chains
      !and the location of the node to which these groups of nodes should be moved.
      Shifted = .TRUE.
      Start = 1
      DO WHILE(Shifted)
        Shifted = .FALSE.

        DO i=Start,CurrentPath % NumberOfNodes
          IF(PathMoveNode(i)) THEN
            IF(.NOT. Shifted) THEN
              Shifted = .TRUE.
              FirstIndex = i
            END IF
            LastIndex = i
          ELSE
            IF(Shifted) EXIT
          END IF
        END DO
        IF(.NOT. Shifted) EXIT

        !We have identified a chain from FirstIndex to LastIndex which need to be moved.
        !They should be moved to either FirstIndex-1 or LastIndex+1
        !(Whichever is further back)
        !Note special case at start and end of path
        IF(FirstIndex == 1) THEN
          ShiftToMe = CurrentPath % NodeNumbers(LastIndex+1)
        ELSE IF(LastIndex == CurrentPath % NumberOfNodes) THEN
          ShiftToMe = CurrentPath % NodeNumbers(FirstIndex-1)
        ELSE IF(Mesh % Nodes % z(CurrentPath % NodeNumbers(FirstIndex-1)) <&
             Mesh % Nodes % z(CurrentPath % NodeNumbers(LastIndex+1))) THEN
          ShiftToMe = CurrentPath % NodeNumbers(FirstIndex-1)
        ELSE
          ShiftToMe = CurrentPath % NodeNumbers(LastIndex+1)
        END IF

        Mesh % Nodes % y(CurrentPath % NodeNumbers(FirstIndex:LastIndex)) = &
             Mesh % Nodes % y(ShiftToMe)

        IF(Debug) PRINT *,'Shifting nodes ',FirstIndex,' to ',LastIndex,&
             ' to point: ',Mesh % Nodes % y(ShiftToMe)
        Start = LastIndex + 1
      END DO

      DEALLOCATE(PathMoveNode)
      CurrentPath => CurrentPath % Next
    END DO

    !NOTE: probably not really necessary here, Shifted nodes don't extend
    !the extent
    !Update Left, Right & Extent
    CALL ComputePathExtent(CrevassePaths, Mesh % Nodes, .TRUE.)

    !--------------------------------------------------------
    ! Remove crevassepaths which are contained within others.
    !--------------------------------------------------------
    !  1) All crevasse paths start and end on the calving front.
    !  2) Crevasse paths can't cross each other.
    !
    !  Thus, iff a crevasse path is surrounded laterally by
    !  another single crevasse path, we remove it, because
    !  it must be contained by the larger one.
    !--------------------------------------------------------

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))

       OtherPath => CrevassePaths
       DO WHILE(ASSOCIATED(OtherPath))
          IF(ASSOCIATED(OtherPath, CurrentPath)) THEN
             OtherPath => OtherPath % Next
             CYCLE
          END IF

          IF((CurrentPath % Left >= OtherPath % Left) .AND. &
               (CurrentPath % Right <= OtherPath % Right)) THEN!contained within
             CurrentPath % Valid = .FALSE.
             IF(Debug) PRINT *,'Debug, marked path ',CurrentPath % ID,' for deletion &
                  &because its contained within path ',OtherPath % ID
          END IF
          OtherPath => OtherPath % Next
       END DO

       CurrentPath => CurrentPath % Next
    END DO

    !Actually remove previous marked
    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))
       WorkPath => CurrentPath % Next

       IF(.NOT. CurrentPath % Valid) THEN
          IF(ASSOCIATED(CurrentPath,CrevassePaths)) CrevassePaths => WorkPath
          CALL RemoveCrevassePath(CurrentPath)
          IF(Debug) CALL Info("ValidateCrevassePaths","Removing a crevasse path")
       END IF
       CurrentPath => WorkPath
    END DO

    !-------------------------------------------------
    ! Check for paths partly obscuring each other
    !  (fully obscured are dealt with above)
    !-------------------------------------------------
    ! If paths partially overlap, the overlapping nodes
    ! of whichever path is seaward are moved.
    ! i.e. the larger calving event takes precedent
    !-------------------------------------------------

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))

      OtherPath => CrevassePaths
      DO WHILE(ASSOCIATED(OtherPath))
        IF(ASSOCIATED(OtherPath, CurrentPath)) THEN
          OtherPath => OtherPath % Next
          CYCLE
        END IF

        IF((CurrentPath % Left < OtherPath % Right) .EQV. &
             (OtherPath % Left < CurrentPath % Right)) THEN !overlap

          IF(Debug) PRINT *,'Debug, paths: ',CurrentPath % ID, OtherPath % ID,' partially overlap'

          !Is the other path to the right or left?
          OtherRight = CurrentPath % Right < OtherPath % Right

          !Check not fully contained - should have been dealt with above
          IF((CurrentPath % Right > OtherPath % Right) .NEQV. &
               (CurrentPath % Left > OtherPath % Left)) THEN
            CALL Warn("ValidateCrevassePaths","Encountered full overlap which &
                 &should already have been taken care of! OK if this is rare, &
                 &otherwise maybe programming error")
          END IF

          IF(OtherRight) THEN
            RightPath => OtherPath
            LeftPath => CurrentPath
          ELSE
            RightPath => CurrentPath
            LeftPath => OtherPath
          END IF

          !Find the left and rightmost nodes of the two paths
          DO i=1,LeftPath % NumberOfNodes
            IF(Debug) PRINT *,'Debug, node ',i,' of leftpath: ',&
                 Mesh % Nodes % y(LeftPath % NodeNumbers(i)), LeftPath % Right

            IF(Mesh % Nodes % y(LeftPath % NodeNumbers(i)) >= LeftPath % Right) LeftIdx = i
          END DO

          DO i=1,RightPath % NumberOfNodes
            IF(Debug) PRINT *,'Debug, node ',i,' of rightpath: ',&
                 Mesh % Nodes % y(RightPath % NodeNumbers(i)), RightPath % Left

            IF(Mesh % Nodes % y(RightPath % NodeNumbers(i)) <= RightPath % Left) RightIdx = i
          END DO

          !See which is further forward.
          ShiftRightPath = Mesh % Nodes % z(LeftPath % NodeNumbers(LeftIdx)) < &
               Mesh % Nodes % z(RightPath % NodeNumbers(RightIdx))

          IF(ShiftRightPath) THEN
            ShiftTo = Mesh % Nodes % y(LeftPath % NodeNumbers(LeftIdx))
            DO i=1,RightPath % NumberOfNodes
              IF(Mesh % Nodes % y(RightPath % NodeNumbers(i)) < ShiftTo) THEN
                IF(Debug) PRINT *,'Debug, overlap shifting right node ',i,' path '&
                     ,RightPath % ID,' from ', Mesh % Nodes % y(RightPath % NodeNumbers(i)),&
                     ' to ',ShiftTo
                Mesh % Nodes % y(RightPath % NodeNumbers(i)) = ShiftTo
              END IF
            END DO
            CALL ComputePathExtent(RightPath, Mesh % Nodes, .FALSE.)

          ELSE
            ShiftTo = Mesh % Nodes % y(RightPath % NodeNumbers(RightIdx))
            DO i=1,LeftPath % NumberOfNodes
              IF(Mesh % Nodes % y(LeftPath % NodeNumbers(i)) > ShiftTo) THEN
                IF(Debug) PRINT *,'Debug, overlap shifting left node ',i,' path ',&
                     LeftPath % ID,' from ',Mesh % Nodes % y(LeftPath % NodeNumbers(i)),&
                     ' to ',ShiftTo
                Mesh % Nodes % y(LeftPath % NodeNumbers(i)) = ShiftTo
              END IF
            END DO
            CALL ComputePathExtent(LeftPath, Mesh % Nodes, .FALSE.)

          END IF
        END IF

        OtherPath => OtherPath % Next
      END DO

      CurrentPath => CurrentPath % Next
    END DO

    !-----------------------------------------------------------------------
    ! Remove elements whose nodes are in a vertical line
    !     (to prevent potential issues in interp)
    !-----------------------------------------------------------------------
    ! This occurs due to the shifting which occurs above.
    ! NOTE: This breaks the assumption that element(i) has nodes (i) & (i+1)
    ! It also breaks the chain! Currently OK but don't rely on this below this
    ! point, or in Calving3D.F90
    !-----------------------------------------------------------------------

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))

      ALLOCATE(DeleteElement(CurrentPath % NumberOfElements),&
           DeleteNode(CurrentPath % NumberOfNodes))
      DeleteElement = .FALSE.
      DeleteNode = .FALSE.

      DO i=1,CurrentPath % NumberOfElements
        !Element i is composed of nodes i,i+1
        IF(Mesh % Nodes % y(CurrentPath % NodeNumbers(i)) == &
             Mesh % Nodes % y(CurrentPath % NodeNumbers(i+1))) THEN
          DeleteElement(i) = .TRUE.
          IF(Debug) PRINT *,'Debug, deleting element: ',i,' from path: ',&
               CurrentPath % ID,' because its a straight line'
        END IF
      END DO

      IF(DeleteElement(1)) DeleteNode(1) = .TRUE.
      IF(DeleteElement(SIZE(DeleteElement))) DeleteNode(SIZE(DeleteNode)) = .TRUE.

      DO i=2,CurrentPath % NumberOfNodes-1
        IF(DeleteElement(i-1) .AND. DeleteElement(i)) DeleteNode(i) = .TRUE.
      END DO

      !Delete them
      IF(COUNT(DeleteElement) > 0) THEN
        !elements
        ALLOCATE(WorkInt(COUNT(.NOT. DeleteElement)))
        WorkInt = PACK(CurrentPath % ElementNumbers,.NOT.DeleteElement)

        DEALLOCATE(CurrentPath % ElementNumbers)
        ALLOCATE(CurrentPath % ElementNumbers(SIZE(WorkInt)))

        CurrentPath % ElementNumbers = WorkInt
        CurrentPath % NumberOfElements = SIZE(WorkInt)
        DEALLOCATE(WorkInt)

        !nodes
        ALLOCATE(WorkInt(COUNT(.NOT. DeleteNode)))
        WorkInt = PACK(CurrentPath % NodeNumbers, .NOT.DeleteNode)

        DEALLOCATE(CurrentPath % NodeNumbers)
        ALLOCATE(CurrentPath % NodeNumbers(SIZE(WorkInt)))

        CurrentPath % NodeNumbers = WorkInt
        CurrentPath % NumberOfNodes = SIZE(WorkInt)
        DEALLOCATE(WorkInt)
      END IF

      DEALLOCATE(DeleteElement, DeleteNode)
      CurrentPath => CurrentPath % Next
    END DO

    !--------------------------------------------------------
    ! Put the mesh back
    !--------------------------------------------------------
    CALL RotateMesh(Mesh, UnRotationMatrix)

  END SUBROUTINE ValidateCrevassePaths

  !Calculates the left and rightmost extent, and the difference (width) of
  !Path, given the node locations in Nodes.
  SUBROUTINE ComputePathExtent(CrevassePaths, Nodes, DoAll)
    TYPE(CrevassePath_t), POINTER :: CrevassePaths
    TYPE(Nodes_t), POINTER :: Nodes
    LOGICAL :: DoAll
    !-----------------------------------------------
    TYPE(CrevassePath_t), POINTER :: CurrentPath
    INTEGER :: n

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))
       CurrentPath % Left = HUGE(1.0_dp)
       CurrentPath % Right = -1.0*HUGE(1.0_dp)

       n = CurrentPath % NumberOfNodes

       CurrentPath % Left = MINVAL(Nodes % y(CurrentPath % NodeNumbers))
       CurrentPath % Right = MAXVAL(Nodes % y(CurrentPath % NodeNumbers))

       CurrentPath % Extent = CurrentPath % Right - CurrentPath % Left

       CurrentPath => CurrentPath % Next

       IF(.NOT. DoAll) EXIT
    END DO

  END SUBROUTINE ComputePathExtent

  !-----------------------------------------------------------------------------
  ! Returns the Path ID of the CrevassePath_t which contains the given element
  !     0 if not found
  !-----------------------------------------------------------------------------
  FUNCTION ElementPathID(CrevassePaths, ElementNo) RESULT(ID)
    TYPE(CrevassePath_t), POINTER :: CrevassePaths
    INTEGER :: ElementNo, ID
    !----------------------------------------------
    TYPE(CrevassePath_t), POINTER :: CurrentPath

    ID = 0

    CurrentPath => CrevassePaths
    DO WHILE(ASSOCIATED(CurrentPath))
       IF(ASSOCIATED(CurrentPath % ElementNumbers)) THEN
          IF(ANY(CurrentPath % ElementNumbers == ElementNo)) THEN
             ID = CurrentPath % ID
             EXIT
          END IF
       END IF
       CurrentPath => CurrentPath % Next
    END DO
    
  END FUNCTION ElementPathID

  !-----------------------------------------------------------------------------
  ! Constructs groups of nodes which fall below a given threshold for some variable
  ! Used with the result of ProjectCalving, it groups nodes which have crevasse 
  ! penetration beyond the threshold.
  !-----------------------------------------------------------------------------
  SUBROUTINE FindCrevasseGroups(Mesh, Variable, Neighbours, Threshold, Groups)
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Variable
    INTEGER, POINTER :: Neighbours(:,:)
    TYPE(CrevasseGroup3D_t), POINTER :: Groups, CurrentGroup
    REAL(KIND=dp) :: Threshold
    !---------------------------------------
    INTEGER :: i, ID
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: VPerm(:)
    INTEGER, ALLOCATABLE :: WorkInt(:)
    LOGICAL, ALLOCATABLE :: Condition(:)
    LOGICAL :: First, Debug
    
    Debug = .FALSE.

    Values => Variable % Values
    VPerm => Variable % Perm

    ALLOCATE(Condition(Mesh % NumberOfNodes))
    DO i=1, Mesh % NumberOfNodes

       IF(VPerm(i) <= 0) THEN
          Condition(i) = .FALSE.
       ELSE IF(Values(VPerm(i)) < Threshold) THEN
          Condition(i) = .TRUE.
       ELSE
          Condition(i) = .FALSE.
       END IF

    END DO

    First = .TRUE.
    ID = 1
    DO i=1,Mesh % NumberOfNodes
       IF(.NOT. Condition(i)) CYCLE

       IF(Debug) PRINT *,'PE:', ParEnv % MyPE,' debug, new group'

       IF(First) THEN
          ALLOCATE(CurrentGroup)
          Groups => CurrentGroup
          First = .FALSE.
       ELSE
          ALLOCATE(CurrentGroup % Next)
          CurrentGroup % Next % Prev => CurrentGroup
          CurrentGroup => CurrentGroup % Next
       END IF
       
       CurrentGroup % ID = ID
       ID = ID + 1

       ALLOCATE(CurrentGroup % NodeNumbers(500))
       CurrentGroup % NumberOfNodes = 1

       !Add node to group and switch it off
       CurrentGroup % NodeNumbers(CurrentGroup % NumberOfNodes) = i
       Condition(i) = .FALSE.

       !Search neighbours
       CALL SearchNeighbours(i, Neighbours, CurrentGroup, Condition)

       ALLOCATE(WorkInt(CurrentGroup % NumberOfNodes))
       WorkInt = CurrentGroup % NodeNumbers(1:CurrentGroup % NumberOfNodes)
       DEALLOCATE(CurrentGroup % NodeNumbers)
       ALLOCATE(CurrentGroup % NodeNumbers(CurrentGroup % NumberOfNodes))
       CurrentGroup % NodeNumbers = WorkInt
       DEALLOCATE(WorkInt)

       CALL UpdateCGrpBB(CurrentGroup, Mesh)
    END DO

    IF(Debug) THEN
       CurrentGroup => Groups
       i=1
       DO WHILE(ASSOCIATED(CurrentGroup))
          PRINT *,'group: ',i,' has ', CurrentGroup % NumberOfNodes,' nodes.'
          i = i + 1
          CurrentGroup => CurrentGroup % Next
       END DO
    END IF
  END SUBROUTINE FindCrevasseGroups

  SUBROUTINE DeallocateCrevasseGroup(CGrp)
    TYPE(CrevasseGroup3D_t), POINTER :: CGrp
    
    IF(ASSOCIATED(CGrp % Next)) CGrp % Next % Prev => CGrp % Prev
    IF(ASSOCIATED(CGrp % Prev)) CGrp % Prev % Next => CGrp % Next
    
    IF(ASSOCIATED(CGrp % NodeNumbers)) DEALLOCATE(CGrp % NodeNumbers)
    IF(ASSOCIATED(CGrp % FrontNodes)) DEALLOCATE(CGrp % FrontNodes)
    IF(ASSOCIATED(CGrp % BoundaryNodes)) DEALLOCATE(CGrp % BoundaryNodes)

    DEALLOCATE(CGrp)

  END SUBROUTINE DeallocateCrevasseGroup

  !Update the Bounding Box of a CrevasseGroup
  SUBROUTINE UpdateCGrpBB(CGrp, Mesh)
    TYPE(CrevasseGroup3D_t), POINTER :: CGrp
    TYPE(Mesh_t), POINTER :: Mesh

    CGrp % BoundingBox(1) = MINVAL(Mesh % Nodes % x(CGrp % NodeNumbers))
    CGrp % BoundingBox(2) = MAXVAL(Mesh % Nodes % x(CGrp % NodeNumbers))
    CGrp % BoundingBox(3) = MINVAL(Mesh % Nodes % y(CGrp % NodeNumbers))
    CGrp % BoundingBox(4) = MAXVAL(Mesh % Nodes % y(CGrp % NodeNumbers))

  END SUBROUTINE UpdateCGrpBB

  !Add a list of points to a CrevasseGroup3D object
  !Don't need to pass the mesh because we're just adding 
  !point indices
  SUBROUTINE AddNodesToGroup(Group, Points, PointCount)
    TYPE(CrevasseGroup3D_t), POINTER :: Group
    INTEGER :: Points(:)
    INTEGER, POINTER :: NewNodeNumbers(:)
    INTEGER :: PointCount, NewNumberOfNodes
    
    NewNumberOfNodes = Group % NumberOfNodes + PointCount
    ALLOCATE(NewNodeNumbers(NewNumberOfNodes))

    NewNodeNumbers(1:Group % NumberOfNodes) = Group % NodeNumbers
    NewNodeNumbers(Group % NumberOfNodes+1:NewNumberOfNodes) = Points(1:PointCount) 

    !Update count
    Group % NumberOfNodes = NewNumberOfNodes
    
    !Point Group to new node list
    DEALLOCATE(Group % NodeNumbers)
    Group % NodeNumbers => NewNodeNumbers
    NULLIFY(NewNodeNumbers)
  END SUBROUTINE AddNodesToGroup

  !------------------------------------------------------------
  ! Routine to recursively search neighbours and put them 
  ! in the current group
  ! Adapted from 2D Calving
  !------------------------------------------------------------
  RECURSIVE SUBROUTINE SearchNeighbours(nodenum, Neighbours, Group, Condition)
    INTEGER :: nodenum
    INTEGER, POINTER :: Neighbours(:,:)
    TYPE(CrevasseGroup3D_t), POINTER :: Group
    LOGICAL, ALLOCATABLE :: Condition(:)
    !------------------------------------------------
    INTEGER :: i, neighbourindex, NoNeighbours

    NoNeighbours = COUNT(Neighbours(nodenum,:) > 0)
    DO i = 1,NoNeighbours
       neighbourindex = Neighbours(nodenum,i)
       IF(.NOT. Condition(neighbourindex)) CYCLE

       Group % NumberOfNodes = Group % NumberOfNodes + 1

       !check space
       IF(Group % NumberOfNodes > SIZE(Group % NodeNumbers)) THEN
          PRINT *, 'Debug, need more space, allocating: ', 2*SIZE(Group % NodeNumbers)
          CALL DoubleIntVectorSize(Group % NodeNumbers)
          PRINT *, 'Debug, new size: ', SIZE(Group % NodeNumbers)
       END IF
       
       Group % NodeNumbers(Group % NumberOfNodes) = neighbourindex

       !Switch it off so it doesn't get readded
       Condition(neighbourindex) = .FALSE.

       CALL SearchNeighbours(neighbourindex, Neighbours, Group, Condition)
    END DO

  END SUBROUTINE SearchNeighbours

  !Marks recursive neighbours with same int
  RECURSIVE SUBROUTINE MarkNeighbours(nodenum, Neighbours, Array, Mark)
    INTEGER :: nodenum
    INTEGER, POINTER :: Array(:)
    LOGICAL, POINTER :: Neighbours(:,:)
    !------------------------------------------------
    INTEGER :: i, Mark

    DO i = 1,SIZE(Neighbours,1)
       IF(.NOT. Neighbours(nodenum,i)) CYCLE
       IF(Array(i)==Mark) CYCLE !already got

       Array(i) = Mark
       CALL MarkNeighbours(i, Neighbours, Array, Mark)
    END DO

  END SUBROUTINE MarkNeighbours

  !-------------------------------------------------------------
  ! Given a CrevasseGroup3D object, finds and stores boundary nodes
  ! BoundaryMask is a logical array TRUE where node sits on a 
  ! mesh (not group) boundary
  ! Note: Not used
  !-------------------------------------------------------------
  SUBROUTINE GetGroupBoundaryNodes(Group, Neighbours, BoundaryMask)
    TYPE(CrevasseGroup3D_t), POINTER :: Group
    INTEGER, POINTER :: Neighbours(:,:)
    LOGICAL :: BoundaryMask(:)
    !-----------------------------------------
    INTEGER :: i, j, node, BNodes, NoNeighbours, neighbour
    INTEGER, ALLOCATABLE :: WorkInt(:)
    LOGICAL :: IsBoundaryNode

    IF(ASSOCIATED(Group % BoundaryNodes)) &
         DEALLOCATE(Group % BoundaryNodes)

    ALLOCATE(Group % BoundaryNodes(100))
    Group % BoundaryNodes = 0
    BNodes = 0

    DO i=1, Group % NumberOfNodes
       IsBoundaryNode = .FALSE.
       node = Group % NodeNumbers(i)

       IF(BoundaryMask(node)) THEN
          IsBoundaryNode = .TRUE.
       ELSE
          NoNeighbours = COUNT(Neighbours(node, :) > 0)
          DO j=1,NoNeighbours
             neighbour = Neighbours(node, j)
             IF(ANY(Group % NodeNumbers == neighbour)) CYCLE

             !Only get here if there's a node NOT in the group
             IsBoundaryNode = .TRUE.
             EXIT
          END DO
       END IF

       IF(IsBoundaryNode) THEN
          BNodes = BNodes + 1
          IF(BNodes > SIZE(Group % BoundaryNodes)) &
               CALL DoubleIntVectorSize(Group % BoundaryNodes)
          Group % BoundaryNodes(BNodes) = node
       END IF
    END DO

    ALLOCATE(WorkInt(BNodes))
    WorkInt = Group % BoundaryNodes(1:BNodes)
    DEALLOCATE(Group % BoundaryNodes)
    ALLOCATE(Group % BoundaryNodes(BNodes))
    Group % BoundaryNodes = WorkInt
    DEALLOCATE(WorkInt)

    !TODO: Order boundary nodes (clockwise?)
  END SUBROUTINE GetGroupBoundaryNodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function to detect if a given node lies within
  ! a 3D crevasse group (physically, not 'graph'ically
  ! Note: not used...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION NodeInCrevasseGroup(NodeNumber, Nodes, CrevasseGroup) RESULT(InGroup)
    INTEGER :: NodeNumber
    TYPE(Nodes_t) :: Nodes
    TYPE(CrevasseGroup3D_t) :: CrevasseGroup
!--------------------------------------------------
    LOGICAL :: InGroup
    REAL(KIND=dp) :: node_x,node_y, BB(4)

    IF(ANY(CrevasseGroup % NodeNumbers == NodeNumber)) & 
         CALL Fatal("NodeInCrevasseGroup", "Scanning for node which&
         &belongs to CrevasseGroup. This is not intended usage!")

    node_x = Nodes % x(NodeNumber)
    node_y = Nodes % y(NodeNumber)

    BB = CrevasseGroup % BoundingBox

    IF(node_x < BB(1) .OR. node_x > BB(2) .OR. &
         node_y < BB(3) .OR. node_y > BB(4)) THEN
       
       InGroup = .FALSE.
       RETURN

    END IF

    CALL Fatal("NodeInCrevasseGroup", "Haven't finished implementing this yet!")

    !Recursively look at node neighbours, stopping when we reach a group member,
    !until we reach freedom (a mesh boundary node) or give up (node is contained
    !within crevassegroup, and this tells us about the topology of the group)

    !RETURN should not just be a logical, this should be repurposed to inform
    !about which boundary it reached.
    
  END FUNCTION NodeInCrevasseGroup

  !Doubles the size of a pointer integer array
  !This version takes a Pointer argument, should
  !be used with care...
  SUBROUTINE DoubleIntVectorSizeP(Vec, fill)
    INTEGER, POINTER :: Vec(:)
    INTEGER, OPTIONAL :: fill
    !----------------------------------------
    INTEGER, ALLOCATABLE :: WorkVec(:)
    
    ALLOCATE(WorkVec(SIZE(Vec)))
    WorkVec = Vec
    
    DEALLOCATE(Vec)
    ALLOCATE(Vec(2*SIZE(WorkVec)))

    IF(PRESENT(fill)) THEN
       Vec = fill
    ELSE
       Vec = 0
    END IF

    Vec(1:SIZE(WorkVec)) = WorkVec

  END SUBROUTINE DoubleIntVectorSizeP

  !Doubles the size of a pointer integer array
  !Allocatable array version
  SUBROUTINE DoubleIntVectorSizeA(Vec, fill)
    INTEGER, ALLOCATABLE :: Vec(:)
    INTEGER, OPTIONAL :: fill
    !----------------------------------------
    INTEGER, ALLOCATABLE :: WorkVec(:)
    
    ALLOCATE(WorkVec(SIZE(Vec)))
    WorkVec = Vec
    
    DEALLOCATE(Vec)
    ALLOCATE(Vec(2*SIZE(WorkVec)))

    IF(PRESENT(fill)) THEN
       Vec = fill
    ELSE
       Vec = 0
    END IF

    Vec(1:SIZE(WorkVec)) = WorkVec

  END SUBROUTINE DoubleIntVectorSizeA


  !-----------------------------------------------------------------------------
  !Given a Nodes_t object, removes the nodes specified by RemoveLogical array
  !Optionally, user may provide a list of node numbers (NodeNums), from which
  !relevant nodes will also be removed
  SUBROUTINE RemoveNodes(InNodes, RemoveLogical, NodeNums)
    TYPE(Nodes_t) :: InNodes, WorkNodes
    LOGICAL, ALLOCATABLE :: RemoveLogical(:)
    INTEGER :: i,counter
    INTEGER, POINTER, OPTIONAL :: NodeNums(:)
    INTEGER, ALLOCATABLE :: WorkNodeNums(:)

    WorkNodes % NumberOfNodes = SIZE(InNodes % x) - COUNT(RemoveLogical)

    ALLOCATE(WorkNodes % x(WorkNodes % NumberOfNodes),&
         WorkNodes % y(WorkNodes % NumberOfNodes),&
         WorkNodes % z(WorkNodes % NumberOfNodes))
    IF(PRESENT(NodeNums)) ALLOCATE(WorkNodeNums(WorkNodes % NumberOfNodes))

    counter = 1
    DO i=1,InNodes % NumberOfNodes
       IF(.NOT. RemoveLogical(i)) THEN
          WorkNodes % x(counter) = InNodes % x(i)
          WorkNodes % y(counter) = InNodes % y(i)
          WorkNodes % z(counter) = InNodes % z(i)
          WorkNodeNums(counter) = NodeNums(i)

          counter = counter + 1
       END IF
    END DO

    DEALLOCATE(InNodes % x, InNodes % y, InNodes % z )
    ALLOCATE(InNodes % x(WorkNodes % NumberOfNodes), &
         InNodes % y(WorkNodes % NumberOfNodes), &
         InNodes % z(WorkNodes % NumberOfNodes))

    IF(PRESENT(NodeNums)) THEN
       DEALLOCATE(NodeNums)
       ALLOCATE(NodeNums(WorkNodes % NumberOfNodes))
    END IF

    InNodes % NumberOfNodes = WorkNodes % NumberOfNodes
    InNodes % x = WorkNodes % x
    InNodes % y = WorkNodes % y
    InNodes % z = WorkNodes % z
    NodeNums = WorkNodeNums

    DEALLOCATE(WorkNodes % x, WorkNodes % y, WorkNodes % z)
    IF(PRESENT(NodeNums)) DEALLOCATE(WorkNodeNums)
  END SUBROUTINE RemoveNodes

  !------------------------------------------------------------------------------
  !> Sort an index array, and change the order of an real array accordingly.
  !> Stolen from GeneralUtils, modified so as to leave the initial index array in order
  !------------------------------------------------------------------------------
  SUBROUTINE MySortF( n,c,b )
    !------------------------------------------------------------------------------
    INTEGER :: n,c(:)
    INTEGER, ALLOCATABLE :: a(:)
    REAL(KIND=dp) :: b(:)
    !------------------------------------------------------------------------------

    INTEGER :: i,j,l,ir,ra
    REAL(KIND=dp) :: rb
    !------------------------------------------------------------------------------

    ALLOCATE(a(SIZE(c)))
    a = c

    IF ( n <= 1 ) RETURN

    l = n / 2 + 1
    ir = n
    DO WHILE( .TRUE. )

       IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
       ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
             a(1) = ra
             b(1) = rb
             RETURN
          END IF
       END IF
       i = l
       j = l + l
       DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
             IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
             a(i) = a(j)
             b(i) = b(j)
             i = j
             j = j + i
          ELSE
             j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
    END DO

    DEALLOCATE(a)

    !------------------------------------------------------------------------------
  END SUBROUTINE MySortF
  !------------------------------------------------------------------------------



  !Returns an ordered list of nodenumbers which specify an edge of a domain,
  !where the edge is determined by the overlap between the two provided permutations
  !NOTE: Returned domain edge is valid only on boss partition (PE=0)
  SUBROUTINE GetDomainEdge(Model, Mesh, TopPerm, EdgeMaskName, OrderedNodes, OrderedNodeNums, &
       Parallel, Simplify, MinDist) 

    IMPLICIT NONE

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: TopPerm(:)
    CHARACTER(MAX_NAME_LEN) :: EdgeMaskName
    TYPE(Nodes_t) :: OrderedNodes, UnorderedNodes
    LOGICAL :: Parallel
    LOGICAL, OPTIONAL :: Simplify
    REAL(KIND=dp), OPTIONAL :: MinDist
    !----------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(NeighbourList_T), ALLOCATABLE :: PartNeighbourList(:)
    INTEGER :: i,j,k,m,n,prev,next,part_start,find_start,find_fin,find_stride,put_start,&
         put_fin, counter,NoNodes, NoNodesOnEdge, NoNeighbours, neigh, Segments, TotSegSplits, &
         direction, index, segnum, soff, foff, target_nodenum, next_nodenum, EdgeBCtag, GlobalNN
    INTEGER :: comm, ierr !MPI stuff
    INTEGER, POINTER :: UnorderedNodeNums(:)=>NULL(), OrderedNodeNums(:), &
         UOGlobalNodeNums(:)=>NULL(), OrderedGlobalNodeNums(:)=>NULL()
    INTEGER, ALLOCATABLE :: NeighbourPartsList(:), PartNodesOnEdge(:), &
         disps(:), nodenum_disps(:), PartOrder(:,:), MyCornerNodes(:), MyNeighbourParts(:), &
         NewSegStart(:), PartSegments(:), SegStarts_Gather(:), WorkInt(:), NodeNeighbours(:,:), &
         GlobalCorners(:), CornerParts(:), PCornerCounts(:)
    LOGICAL :: Debug, ActivePart, Boss, Simpl, NotThis, Found, ThisBC
    LOGICAL, ALLOCATABLE :: OnEdge(:), ActivePartList(:), RemoveNode(:), IsCornerNode(:)
    REAL(KIND=dp) :: prec
    CHARACTER(MAX_NAME_LEN) :: FuncName

    TYPE AllocIntList_t
       INTEGER, DIMENSION(:), POINTER :: Indices
    END TYPE AllocIntList_t
    TYPE(AllocIntList_t), ALLOCATABLE :: PartSegStarts(:)
    
    !------------------------------------------------
    ! Change in strategy:
    !
    ! Previously, used stiffness matrix to determine connectivity, but 
    ! this causes problems when multiple nodes on the boundary reside
    ! in the same top surface tri element:
    !
    !            *===*===*---*===*===*
    !    from this one---^\ /
    !                      *
    !                      ^-- we want this one
    !
    !  Various versions of this issue can occur...
    !
    !  SO, instead of using the stiffness matrix, we should
    !  check all the boundary elements on the relevant SIDE
    !  boundary (e.g. calving front, right sidewall...), 
    !  looking for elements containing nodes for which the 
    !  top mask is true.
    ! 
    !  Because of the extruded structure of the mesh, nodes
    !  within the same boundary quad will always be neighbours,
    !  and each node shall have no more than 2 neighbours.
    !----------------------------------------------------

    FuncName = "GetDomainEdge"
    Debug = .FALSE.
    ActivePart = .TRUE.

    NoNodes = SIZE(TopPerm) !total number of nodes in domain/partition

    IF(Parallel) THEN
       comm = ELMER_COMM_WORLD
       Boss = (ParEnv % MyPE == 0)
    ELSE
       Boss = .TRUE. !only one part in serial, so it's in charge of computation
    END IF

    IF(PRESENT(Simplify)) THEN
       Simpl = Simplify
    ELSE 
       Simpl = .FALSE.
    END IF

    ALLOCATE(OnEdge(NoNodes), NodeNeighbours(NoNodes,2))
    OnEdge = .FALSE.
    NodeNeighbours = -1

    !Find correct BC from logical
    DO i=1,Model % NumberOfBCs
       ThisBC = ListGetLogical(Model % BCs(i) % Values,EdgeMaskName,Found)
       IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
       EdgeBCtag =  Model % BCs(i) % Tag
       EXIT
    END DO

    !Cycle boundary elements, marking nodes on edge and finding neighbours
    DO i=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       IF(Element % BoundaryInfo % Constraint /= EdgeBCtag) CYCLE !elem not on lateral boundary
       IF(.NOT. ANY(TopPerm(Element % NodeIndexes) > 0)) CYCLE !elem contains no nodes on top

       IF(GetElementFamily(Element) == 1) &
            CALL Fatal(FuncName, "101 Elements are supposed to be a thing of the past!")

       !Cycle nodes in element
       DO j=1,Element % TYPE % NumberOfNodes
          IF(.NOT. TopPerm(Element % NodeIndexes(j)) > 0) CYCLE
          OnEdge(Element % NodeIndexes(j)) = .TRUE.

          !Cycle nodes in element
          DO k=1,Element % TYPE % NumberOfNodes
             IF(j==k) CYCLE
             IF(.NOT. TopPerm(Element % NodeIndexes(k))>0) CYCLE
             DO m=1,2 !fill NodeNeighbours
                IF(NodeNeighbours(Element % NodeIndexes(j),m) /= -1) CYCLE
                NodeNeighbours(Element % NodeIndexes(j),m) = Element % NodeIndexes(k)
                EXIT
             END DO
             IF(.NOT. ANY(NodeNeighbours(Element % NodeIndexes(j),:) == Element % NodeIndexes(k))) &
                  CALL Fatal(FuncName,'Identified more than two neighbours')
          END DO
       END DO

    END DO

    NoNodesOnEdge = COUNT(OnEdge)
    IF(NoNodesOnEdge == 1) THEN
      CALL Fatal(FuncName, "A single node identified on boundary, shouldnt be possible. &
           &Someone is messing around with 101 elements.")
    END IF

    ALLOCATE(UnorderedNodeNums(NoNodesOnEdge),&
         OrderedNodeNums(NoNodesOnEdge))
    OrderedNodeNums = -1 !initialize to invalid value

    j = 0
    DO i=1,NoNodes
       IF(.NOT. OnEdge(i)) CYCLE
       j = j + 1
       UnorderedNodeNums(j) = i
    END DO

    !Cycle nodes on edge, looking for one with only one neighbour (a corner)
    IF(NoNodesOnEdge > 1) THEN

       ALLOCATE(IsCornerNode(NoNodesOnEdge))
       IsCornerNode = .FALSE.

       DO i=1,NoNodesOnEdge
          IsCornerNode(i) = COUNT(NodeNeighbours(UnOrderedNodeNums(i),:) == -1) == 1
          IF(COUNT(NodeNeighbours(UnOrderedNodeNums(i),:) == -1) == 2) &
               CALL Fatal(FuncName, "Found an isolated node on edge")
       END DO

       IF(MOD(COUNT(IsCornerNode),2) /= 0) THEN
          WRITE(Message,'(A,i0)') "Found an odd number of&
               & corner nodes in partition: ",ParEnv % MyPE
          CALL Fatal(FuncName, Message)
       END IF

       Segments = COUNT(IsCornerNode) / 2
       IF(Debug .AND. Segments > 1) PRINT *, &
            'Partition ',ParEnv % MyPE, ' has ',Segments,' line segments on boundary.'

       ALLOCATE(NewSegStart(Segments-1))
       ALLOCATE(MyCornerNodes(COUNT(IsCornerNode)))

       counter = 1
       DO i=1,NoNodesOnEdge 
          IF(IsCornerNode(i)) THEN
             MyCornerNodes(counter) = i
             counter = counter + 1
          END IF
       END DO

       counter = 1
       DO k=1,Segments
          
          IF(k==1) THEN
             OrderedNodeNums(counter) = UnorderedNodeNums(MyCornerNodes(1))
          ELSE
             DO i=2, SIZE(MyCornerNodes)
                IF(ANY(OrderedNodeNums == UnorderedNodeNums(MyCornerNodes(i)))) THEN
                   CYCLE
                ELSE
                   OrderedNodeNums(counter) = UnorderedNodeNums(MyCornerNodes(i))
                   EXIT
                END IF
             END DO
          END IF
          counter = counter + 1

          !----------------------------------------------------
          !   Move along from corner, filling in order
          !----------------------------------------------------
          DO i=counter,NoNodesOnEdge
             Found = .FALSE.
             IF(OrderedNodeNums(i-1) == -1) CALL Abort()

             DO j=1,2
                IF(NodeNeighbours(OrderedNodeNums(i-1),j) == -1) CYCLE !First and last nodes, corner
                IF(ANY(OrderedNodeNums(1:i-1) == NodeNeighbours(OrderedNodeNums(i-1),j))) &
                     CYCLE !already in list

                OrderedNodeNums(i) = NodeNeighbours(OrderedNodeNums(i-1),j)
                Found = .TRUE.
             END DO

             IF(.NOT. Found) EXIT
          END DO

          counter = i

          IF(counter >= NoNodesOnEdge) EXIT !this should be redundant...
          NewSegStart(k) = counter
       END DO

    ELSE !Either 1 or 0 nodes found, not an active boundary partition
       !0 node case, obvious
       !1 node case, if THIS partition only has one node on the boundary,
       !this same node must be caught by two other partitions, so we aren't needed.
       ALLOCATE(NewSegStart(0), MyCornerNodes(0))
       ActivePart = .FALSE.
       Segments = 0
       NoNodesOnEdge = 0 !simplifies mpi comms later
       IF(.NOT.Parallel) CALL Fatal(FuncName,&
            "Found either 1 or 0 nodes in a serial run, this isn't a valid boundary edge!")
    END IF


    !gather corner count - replaces 101 element detection
    ALLOCATE(PCornerCounts(ParEnv % PEs),disps(ParEnv % PEs))

    CALL MPI_AllGather(SIZE(MyCornerNodes), 1, MPI_INTEGER, PCornerCounts, &
         1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)

    disps(1) = 0
    DO i=2, ParEnv % PEs
       disps(i) = disps(i-1) + PCornerCounts(i-1)
    END DO

    ALLOCATE(GlobalCorners(SUM(PCornerCounts)),&
         CornerParts(SUM(PCornerCounts)))

    !gather corner nodenums
    CALL MPI_AllGatherV(Mesh % ParallelInfo % GlobalDOFs(UnorderedNodeNums(MyCornerNodes)), &
         SIZE(MyCornerNodes), MPI_INTEGER, GlobalCorners, PCornerCounts, disps, &
         MPI_INTEGER, ELMER_COMM_WORLD, ierr)

    !note which partition sent each corner node
    counter = 1
    DO i=1,ParEnv % PEs
       IF(PCornerCounts(i) == 0) CYCLE
       CornerParts(counter:counter+PCornerCounts(i)-1) = i-1
       counter = counter + PCornerCounts(i)
    END DO

    !Quick check:
    DO i=1,SIZE(GlobalCorners)
      counter = COUNT(GlobalCorners == GlobalCorners(i))
      IF(counter > 2) CALL Fatal(FuncName,"Programming error in partition segment detection, node found too many times!")
    END DO

    !Now GlobalCorners and CornerParts tell us which partitions found corner nodes
    !(i.e. nodes which will join other segments)

    !Remember that, in parallel, we're using local rather than global node numbers
    IF(Parallel) THEN
       IF(ActivePart) THEN
          ALLOCATE(MyNeighbourParts(Segments*2))

          DO i=1,Segments*2 !Find neighbour partition numbers

             IF(i==1) THEN
                n = OrderedNodeNums(1)
             ELSE IF(i==Segments*2) THEN
                n = OrderedNodeNums(NoNodesOnEdge)
             ELSE IF(MOD(i,2)==0) THEN
                n = OrderedNodeNums(NewSegStart(i/2)-1)
             ELSE
                n = OrderedNodeNums(NewSegStart(i/2))
             END IF
             
             MyNeighbourParts(i) = -1 !default if not caught in loop below
             GlobalNN = Mesh % ParallelInfo % GlobalDOFs(n)
             DO j=1,SIZE(GlobalCorners)
               IF(GlobalCorners(j) /= GlobalNN) CYCLE
               IF(CornerParts(j) == ParEnv % MyPE) CYCLE
               MyNeighbourParts(i) = CornerParts(j)
               IF( .NOT. (ANY(Mesh % ParallelInfo % NeighbourList(n) % Neighbours &
                    == MyNeighbourParts(i)))) CALL Fatal(FuncName, &
                    "Failed sanity check on neighbour partition detection.")
             END DO
          END DO
       ELSE
          ALLOCATE(MyNeighbourParts(0))
       END IF

       IF(Boss) ALLOCATE(PartSegments(ParEnv % PEs))

       CALL MPI_GATHER(Segments, 1, MPI_INTEGER, PartSegments, &
            1, MPI_INTEGER,  0, comm, ierr)

       IF(Boss) THEN

          TotSegSplits = 0
          DO i=1,SIZE(PartSegments)
             TotSegSplits = TotSegSplits + MAX(PartSegments(i)-1,0)
          END DO

          ALLOCATE(nodenum_disps(ParEnv % PEs), &
               PartNodesOnEdge(ParEnv % PEs), &
               NeighbourPartsList(SUM(PartSegments)*2), &
               PartNeighbourList(ParEnv % PEs), &
               SegStarts_Gather(TotSegSplits))

          DO i=1,ParEnv % PEs
             ALLOCATE(PartNeighbourList(i) % Neighbours(PartSegments(i)*2))
          END DO

          disps(1) = 0
          DO i=2, ParEnv % PEs
             disps(i) = disps(i-1) + MAX(PartSegments(i-1)-1,0)
          END DO

       END IF

       !Get found count from each part to boss
       CALL MPI_GATHER(NoNodesOnEdge, 1, MPI_INTEGER, PartNodesOnEdge, &
            1, MPI_INTEGER, 0, comm ,ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       IF(Debug .AND. Boss) THEN
          PRINT *, 'boss size(SegStarts_Gather): ', SIZE(SegStarts_Gather)
          PRINT *, 'boss PartSegments: ', PartSegments
          PRINT *, 'boss disps:', disps
       END IF
       
       IF(Boss) THEN
          ALLOCATE(WorkInt(ParEnv % PEs))
          WorkInt = MAX(PartSegments-1,0)
       END IF

       CALL MPI_GATHERV(NewSegStart, MAX(Segments-1,0), MPI_INTEGER, SegStarts_Gather, &
            WorkInt, disps, MPI_INTEGER, 0, comm, ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       IF(Boss) THEN
          ALLOCATE(PartSegStarts(ParEnv % PEs))
          DO i=1,ParEnv % PEs
             j = PartSegments(i)
             ALLOCATE(  PartSegStarts(i) % Indices(MAX((j - 1),0)))
             IF(j > 1) THEN
                IF(Debug) PRINT *, 'debug disps(i),j', disps(i),j
                PartSegStarts(i) % Indices = SegStarts_Gather(1+disps(i) : (1+disps(i) + (j-1)-1) ) 
             END IF
             IF(Debug) PRINT *, i,' partsegstarts: ', PartSegStarts(i) % Indices
          END DO

          disps(1) = 0
          DO i=2, ParEnv % PEs
             disps(i) = disps(i-1) + PartSegments(i-1)*2
          END DO

          WorkInt = PartSegments*2
       END IF

       !Get neighbour part numbers from each part to boss
       CALL MPI_GATHERV(MyNeighbourParts, Segments*2, MPI_INTEGER, NeighbourPartsList, &
            WorkInt, disps, MPI_INTEGER, 0, comm ,ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       IF(Debug .AND. Boss) PRINT *, 'DEBUG, NewSegStart: ', NewSegStart

       IF(Boss) THEN
          ActivePartList = (PartNodesOnEdge > 0)

          !Here we account for shared nodes on partition boundaries
          OrderedNodes % NumberOfNodes = SUM(PartNodesOnEdge) - (SIZE(NeighbourPartsList)/2 - 1)
          UnorderedNodes % NumberOfNodes = SUM(PartNodesOnEdge) !but they are still present when gathered...

          ALLOCATE(PartOrder(SIZE(NeighbourPartsList)/2,2),&
               OrderedNodes % x(OrderedNodes % NumberOfNodes),&
               OrderedNodes % y(OrderedNodes % NumberOfNodes),&
               OrderedNodes % z(OrderedNodes % NumberOfNodes),&
               UnorderedNodes % x(UnorderedNodes % NumberOfNodes),&
               UnorderedNodes % y(UnorderedNodes % NumberOfNodes),&
               UnorderedNodes % z(UnorderedNodes % NumberOfNodes),&
               UOGlobalNodeNums(UnorderedNodes % NumberOfNodes),&
               OrderedGlobalNodeNums(OrderedNodes % NumberOfNodes))

          nodenum_disps(1) = 0
          DO i=2, ParEnv % PEs
             nodenum_disps(i) = nodenum_disps(i-1) + PartNodesOnEdge(i-1)
          END DO

          IF(Debug) THEN
             PRINT *, 'debug disps: ', disps
             PRINT *, 'debug nodenum_disps: ', nodenum_disps
             PRINT *, 'debug neighbourpartslist: ',NeighbourPartsList
             PRINT *, 'Partition Segments: ',PartSegments
          END IF
       END IF

       !-----------------------------------------------------------
       ! Gather node coords from all partitions
       ! Note, they're going into 'UnorderedNodes': though they are ordered
       ! within their partition, the partitions aren't ordered...
       !-----------------------------------------------------------

       !Global Node Numbers
       CALL MPI_GATHERV(Mesh % ParallelInfo % GlobalDOFs(OrderedNodeNums),&
            NoNodesOnEdge,MPI_INTEGER,&
            UOGlobalNodeNums,PartNodesOnEdge,&
            nodenum_disps,MPI_INTEGER,0,comm, ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       !X coords
       CALL MPI_GATHERV(Mesh % Nodes % x(OrderedNodeNums),&
            NoNodesOnEdge,MPI_DOUBLE_PRECISION,&
            UnorderedNodes % x,PartNodesOnEdge,&
            nodenum_disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       !Y coords
       CALL MPI_GATHERV(Mesh % Nodes % y(OrderedNodeNums),&
            NoNodesOnEdge,MPI_DOUBLE_PRECISION,&
            UnorderedNodes % y,PartNodesOnEdge,&
            nodenum_disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       !Z coords
       CALL MPI_GATHERV(Mesh % Nodes % z(OrderedNodeNums),&
            NoNodesOnEdge,MPI_DOUBLE_PRECISION,&
            UnorderedNodes % z,PartNodesOnEdge,&
            nodenum_disps,MPI_DOUBLE_PRECISION,0,comm, ierr)
       IF(ierr /= MPI_SUCCESS) CALL Fatal(FuncName,"MPI Error!")
       CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       !-----------------------------------------------------------
       ! Determine order of partitions by linking neighbours and
       ! checking globalnodenumbers where appropriate
       !-----------------------------------------------------------

       IF(Boss) THEN
          !Notes: NeighbourPartsList is zero indexed, like PEs
          !PartOrder is 1 indexed
          !disps is 1 indexed. So disps(NeighbourPartsList+1)

          PartOrder = 0 !init
          direction = 0
          prev = -1

          !First fill in PartNeighbourList % Neighbours
          DO i=1,ParEnv % PEs
             IF(PartSegments(i)==0) CYCLE
             PartNeighbourList(i) % Neighbours = &
                  NeighbourPartsList( (1+disps(i)) : (1+disps(i) + (PartSegments(i)*2) - 1) )
             !There is the possibility of missing an end (-1) due to partition
             !landing right on corner
             DO j=1,SIZE(PartNeighbourList(i) % Neighbours)
                IF(PartNeighbourList(i) % Neighbours(j) == -1) CYCLE
                IF(PartSegments(PartNeighbourList(i) % Neighbours(j)+1) < 1) THEN
                   IF(Debug) PRINT *, 'Neighbour ',PartNeighbourList(i) % Neighbours(j)+1,&
                        "isn't really on boundary, so changing to -1"
                   PartNeighbourList(i) % Neighbours(j) = -1
                END IF
             END DO

             IF(Debug) PRINT *, i, ': Neighbours: ', PartNeighbourList(i) % Neighbours
             !find a corner partition
             IF(ANY(PartNeighbourList(i) % Neighbours == prev)) next = i
          END DO

          IF(Debug) THEN
             PRINT *, 'Debug GetDomainEdge, globalno, unorderednodes % x: '
             DO i=1,SIZE(UOGlobalNodeNums)
                PRINT *, i, UOGlobalNodeNums(i), UnorderedNodes % x(i)
             END DO

             PRINT *, 'debug nodenum_disps: '
             DO i=1, SIZE(nodenum_disps)
                PRINT *, i,'  ',nodenum_disps(i)
             END DO
          END IF


          counter = 1

          DO WHILE(.TRUE.)
             IF((COUNT(PartNeighbourList(next) % Neighbours == prev) == 1) .OR. &
                (prev == -1)) THEN
                DO j=1,SIZE(PartNeighbourList(next) % Neighbours)
                   IF(PartNeighbourList(next) % Neighbours(j) == prev) THEN
                      index = j
                      EXIT
                   END IF
                END DO
             ELSE !Neighbours on both sides, so need to inspect globalnodenumbers
                IF(Debug) PRINT *, 'debug, two matches'
                DO j=1,SIZE(PartNeighbourList(next) % Neighbours)
                   IF(PartNeighbourList(next) % Neighbours(j) == prev) THEN

                      segnum = ((j-1)/2) + 1
                      direction = (2 * MOD(j, 2)) - 1
                      
                      IF(segnum == 1) THEN
                         soff = 0
                      ELSE
                         soff = PartSegStarts(next) % Indices(segnum - 1) - 1
                      END IF
                      IF(segnum == PartSegments(next)) THEN
                         foff = 0
                      ELSE
                         foff = -1 * (PartNodesOnEdge(next) - PartSegStarts(next) % Indices(segnum) + 1)
                      END IF

                      IF(direction > 0) THEN
                         next_nodenum = UOGlobalNodeNums(1 + nodenum_disps(next) + soff)
                      ELSE
                         !one node before (-1) the next partition's (+1) nodes
                         IF(next == ParEnv % PEs) THEN
                            k = SIZE(UOGlobalNodeNums)
                         ELSE
                            k = 1 + nodenum_disps(next+1) - 1
                         END IF
                         next_nodenum = UOGlobalNodeNums(k + foff)
                      END IF
                      IF(Debug) THEN
                         PRINT *, 'debug, next_nodenum: ', next_nodenum
                         PRINT *, 'debug, target_nodenum: ', target_nodenum
                      END IF
                      IF(next_nodenum == target_nodenum) THEN
                         index = j
                         EXIT
                      END IF
                   END IF
                END DO
             END IF

             segnum = ((index-1)/2) + 1
             direction = (2 * MOD(index, 2)) - 1
             PartOrder(counter,1) = next - 1
             PartOrder(counter,2) = direction * segnum
             counter = counter + 1

             IF(Debug) THEN
                PRINT *, 'index: ', index
                PRINT *, 'segnum: ', segnum
                PRINT *, 'direction: ',direction
                PRINT *, 'next: ', next
                PRINT *, 'prev: ', prev
             END IF

             prev = next - 1
             j = next
             next = PartNeighbourList(next) % Neighbours(index + direction)

             !In case of two matches, need a target node to find
             IF(segnum == 1) THEN
                soff = 0
             ELSE
                soff = PartSegStarts(j) % Indices(segnum - 1) - 1
             END IF
             IF(segnum == PartSegments(j)) THEN
                foff = 0
             ELSE
                foff = -1 * (PartNodesOnEdge(j) - PartSegStarts(j) % Indices(segnum) + 1)
             END IF

             IF(direction < 0) THEN
                target_nodenum = UOGlobalNodeNums(1 + nodenum_disps(prev+1) + soff)
             ELSE
                IF(prev + 1 == ParEnv % PEs) THEN
                   k = SIZE(UOGlobalNodeNums)
                ELSE
                   k = 1 + nodenum_disps(prev+1+1) - 1
                END IF
                !one node before (-1) the next partition's (+1) nodes
                target_nodenum = UOGlobalNodeNums(k + foff)
             END IF

             !wipe them out so we don't accidentally come back this way
             PartNeighbourList(j) % Neighbours(index:index+direction:direction) = 0

             IF(next == -1) EXIT
             next = next + 1
          END DO

          IF(Debug) PRINT *, 'Debug GetDomainEdge, part order:', PartOrder

       END IF

       !-----------------------------------------------------------
       ! Put nodes collected from partitions into order 
       !-----------------------------------------------------------

       IF(Boss) THEN
          put_start = 1

          DO i=1,SIZE(PartOrder,1)
             j = PartOrder(i,1) + 1
             segnum = PartOrder(i,2)

             IF(j==0) CALL Abort()

             foff = 0
             soff = 0
             IF(PartSegments(j) > 1) THEN
                IF(Debug) THEN
                   PRINT *, 'Debug GetDomainEdge, extracting nodes from segmented partition'
                   PRINT *, 'Debug GetDomainEdge, segnum: ', segnum
                   PRINT *, 'Debug GetDomainEdge, partnodes: ', PartNodesOnEdge(j)
                   PRINT *, 'Debug GetDomainEdge, PartSegStarts(j) % Indices: ',&
                        PartSegStarts(j) % Indices
                   PRINT *, 'Debug GetDomainEdge, nodenum_disps(j): ',nodenum_disps(j)
                END IF

                IF(ABS(segnum) == 1) THEN

                   soff = 0
                ELSE
                   soff = PartSegStarts(j) % Indices(ABS(segnum) - 1) - 1
                END IF
                IF(ABS(segnum) == PartSegments(j)) THEN
                   foff = 0
                ELSE
                   foff = -1 * (PartNodesOnEdge(j) - PartSegStarts(j) % Indices(ABS(segnum)) + 1)
                END IF
             END IF

             part_start = 1 + nodenum_disps(j) !where are this partitions nodes?
             IF(segnum > 0) THEN
                find_start = part_start + soff
                find_fin = part_start + PartNodesOnEdge(j) - 1 + foff
                find_stride = 1
             ELSE
                find_fin = part_start + soff
                find_start = part_start + PartNodesOnEdge(j) - 1 + foff
                find_stride = -1
             END IF

             put_fin = put_start + ABS(find_start - find_fin)
             IF(Debug) THEN
                PRINT *, 'Debug, find start, end: ',find_start, find_fin, find_stride
                PRINT *, 'Debug, put start, end: ',put_start, put_fin
                PRINT *, 'Total slots: ',SIZE(OrderedNodes % x)
             END IF

             OrderedNodes % x(put_start:put_fin) = &
                  UnorderedNodes % x(find_start:find_fin:find_stride)
             OrderedNodes % y(put_start:put_fin) = &
                  UnorderedNodes % y(find_start:find_fin:find_stride)
             OrderedNodes % z(put_start:put_fin) = &
                  UnorderedNodes % z(find_start:find_fin:find_stride)
             OrderedGlobalNodeNums(put_start:put_fin) = &
                  UOGlobalNodeNums(find_start:find_fin:find_stride)

             put_start = put_fin !1 node overlap
          END DO

          DEALLOCATE(OrderedNodeNums)
          ALLOCATE(OrderedNodeNums(OrderedNodes % NumberOfNodes))
          OrderedNodeNums = OrderedGlobalNodeNums

          IF(Debug) THEN
             PRINT *, 'Debug GetDomainEdge, globalno, orderednodes % x: '
             DO i=1,SIZE(OrderedNodes % x)
                PRINT *, OrderedNodeNums(i), OrderedNodes % x(i)
             END DO
          END IF
       END IF

    ELSE !serial
       OrderedNodes % NumberOfNodes = NoNodesOnEdge
       ALLOCATE(OrderedNodes % x(OrderedNodes % NumberOfNodes),&
            OrderedNodes % y(OrderedNodes % NumberOfNodes),&
            OrderedNodes % z(OrderedNodes % NumberOfNodes))

       OrderedNodes % x = Mesh % Nodes % x(OrderedNodeNums)
       OrderedNodes % y = Mesh % Nodes % y(OrderedNodeNums)
       OrderedNodes % z = Mesh % Nodes % z(OrderedNodeNums)

       !No action required on OrderedNodeNums...
    END IF

    !-------------------------------------------------------------
    ! Simplify geometry by removing interior nodes on any straight
    ! lines if requested
    !-------------------------------------------------------------
    IF(Simpl .AND. Boss) THEN
       ALLOCATE(RemoveNode(OrderedNodes % NumberOfNodes))
       RemoveNode = .FALSE.

       DO i=2,OrderedNodes % NumberOfNodes-1 !Test all interior nodes
          IF(Debug) THEN
             PRINT *, (NodesGradXY(OrderedNodes,i,i-1))
             PRINT *, (NodesGradXY(OrderedNodes,i+1,i))
             PRINT *, 'DIFF: ',ABS(NodesGradXY(OrderedNodes,i,i-1) -&
                  NodesGradXY(OrderedNodes,i+1,i))
             PRINT *, ''
          END IF

          !Need to determine numerical precision of input datapoints
          !i.e. after how many decimal places are values constant
          !e.g. 0.23000000... or 99999...
          prec = MAX(RealAeps(OrderedNodes % x(i)),RealAeps(OrderedNodes % y(i)))

          IF(ABS(NodesGradXY(OrderedNodes,i,i-1) - NodesGradXY(OrderedNodes,i+1,i)) < prec) THEN
             RemoveNode(i) = .TRUE.
          END IF
       END DO

       IF(COUNT(RemoveNode) > 0) THEN

          CALL RemoveNodes(OrderedNodes, RemoveNode, OrderedNodeNums)

          IF(Debug) THEN
             PRINT *, 'Debug GetDomainEdge, Simplify removing: ', COUNT(RemoveNode), ' nodes'
             DO i=1,OrderedNodes % NumberOfNodes
                PRINT *, 'Debug GetDomainEdge, node: ',i
                PRINT *, 'x: ',OrderedNodes % x(i),'y: ',OrderedNodes % y(i)
             END DO
          END IF !debug

       END IF !removing any nodes
       DEALLOCATE(RemoveNode)
    END IF !simplify       

    !-------------------------------------------------------------
    ! Remove any nodes which are closer together than MinDist, if
    ! this is specified.
    !-------------------------------------------------------------
    IF(PRESENT(MinDist) .AND. Boss) THEN
       !Cycle all nodes, remove any too close together
       !This won't guarantee that the new domain edge is *within* the old one
       !but could be adapted to do so
       ALLOCATE(RemoveNode(OrderedNodes % NumberOfNodes))
       RemoveNode = .FALSE.
       DO i=2,OrderedNodes % NumberOfNodes-1 !Test all interior nodes
          j = i - 1
          DO WHILE(RemoveNode(j))
             j = j-1
          END DO

          IF(NodeDist2D(OrderedNodes, i, j) < MinDist) THEN
             RemoveNode(i) = .TRUE.
             IF(Debug) THEN
                PRINT *, 'Debug GetDomainEdge, MinDist, removing node ',i,' too close to: ', j
                PRINT *, 'Debug GetDomainEdge, MinDist, dist: ',NodeDist2D(OrderedNodes, i, j)
             END IF
          END IF
       END DO

       IF(COUNT(RemoveNode) > 0) THEN

          CALL RemoveNodes(OrderedNodes, RemoveNode, OrderedNodeNums)

          IF(Debug) THEN
             PRINT *, 'Debug GetDomainEdge, MinDist removing: ', COUNT(RemoveNode), ' nodes'
             DO i=1,OrderedNodes % NumberOfNodes
                PRINT *, 'Debug GetDomainEdge, node: ',i
                PRINT *, 'x: ',OrderedNodes % x(i),'y: ',OrderedNodes % y(i)
             END DO
          END IF !debug

       END IF !removing any nodes
       DEALLOCATE(RemoveNode)
    END IF !MinDist

    !------------ DEALLOCATIONS ------------------

    DEALLOCATE(OnEdge, UnorderedNodeNums, GlobalCorners, CornerParts, PCornerCounts)

    IF(Boss .AND. Parallel) THEN !Deallocations
       DEALLOCATE(UnorderedNodes % x, &
            UnorderedNodes % y, &
            UnorderedNodes % z, &
            PartNodesOnEdge, &
            disps, nodenum_disps, &
            PartOrder, &
            UOGlobalNodeNums, &
            OrderedGlobalNodeNums)
    END IF

  END SUBROUTINE GetDomainEdge

  ! Copies over time variables and creates coordinate vars. Basically pinched
  ! from AddMeshCoordinatesAndTime() and Multigrid
  SUBROUTINE CopyIntrinsicVars(OldMesh, NewMesh)
    IMPLICIT NONE
    
    TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER :: WorkVar
    !----------------------------------------------------------
    NULLIFY( Solver )

    CALL VariableAdd( NewMesh % Variables, NewMesh,Solver, &
         'Coordinate 1',1,NewMesh % Nodes % x )

    CALL VariableAdd(NewMesh % Variables,NewMesh,Solver, &
         'Coordinate 2',1,NewMesh % Nodes % y )

    CALL VariableAdd(NewMesh % Variables,NewMesh,Solver, &
         'Coordinate 3',1,NewMesh % Nodes % z )
    
    WorkVar => VariableGet( OldMesh % Variables, 'Time', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Time', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Periodic Time', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Periodic Time', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Timestep', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Timestep', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Timestep size', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Timestep size', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Timestep interval', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Timestep interval', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Coupled iter', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Coupled iter', 1, WorkVar % Values )

    WorkVar => VariableGet( OldMesh % Variables, 'Nonlin iter', ThisOnly=.TRUE.)
    CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, 'Nonlin iter', 1, WorkVar % Values )

  END SUBROUTINE CopyIntrinsicVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Function to rotate a mesh by rotationmatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE RotateMesh(Mesh, RotationMatrix)

    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: RotationMatrix(3,3), NodeHolder(3)
    INTEGER :: i

    DO i=1,Mesh % NumberOfNodes
       NodeHolder(1) = Mesh % Nodes % x(i)
       NodeHolder(2) = Mesh % Nodes % y(i)
       NodeHolder(3) = Mesh % Nodes % z(i)

       NodeHolder = MATMUL(RotationMatrix,NodeHolder)

       Mesh % Nodes % x(i) = NodeHolder(1)
       Mesh % Nodes % y(i) = NodeHolder(2)
       Mesh % Nodes % z(i) = NodeHolder(3)
    END DO

  END SUBROUTINE RotateMesh
  
  SUBROUTINE DeallocateElement(Element)
   
    IMPLICIT NONE
    TYPE(Element_t) :: Element

    IF ( ASSOCIATED( Element % NodeIndexes ) ) &
         DEALLOCATE( Element % NodeIndexes )
    Element % NodeIndexes => NULL()

    IF ( ASSOCIATED( Element % EdgeIndexes ) ) &
         DEALLOCATE( Element % EdgeIndexes )
    Element % EdgeIndexes => NULL()

    IF ( ASSOCIATED( Element % FaceIndexes ) ) &
         DEALLOCATE( Element % FaceIndexes )
    Element % FaceIndexes => NULL()

    IF ( ASSOCIATED( Element % DGIndexes ) ) &
         DEALLOCATE( Element % DGIndexes )
    Element % DGIndexes => NULL()

    IF ( ASSOCIATED( Element % BubbleIndexes ) ) &
         DEALLOCATE( Element % BubbleIndexes )
    Element % BubbleIndexes => NULL()

    IF ( ASSOCIATED( Element % PDefs ) ) &
         DEALLOCATE( Element % PDefs )
    Element % PDefs => NULL()

  END SUBROUTINE DeallocateElement

  !Identify front elements connected to the bed, which are sufficiently horizontal
  !to warrant reclassification as basal elements.
  !Note, only does elements currently connected to the bed. i.e. one row per dt
  !Returns:
  !   NewBasalNode(:), LOGICAL true where frontal node becomes basal
  !   ExFrontalNode(:), LOGICAL true where a frontal node no longer 
  !      belongs to its front column (though it may still be on the front...)
  !
  ! NOTE, if an error in this subroutine, could be element
  ! which sits between 2 NewBasalElems

  SUBROUTINE ConvertFrontalToBasal(Model, Mesh, FrontMaskName, BotMaskName, &
       ZThresh, NewBasalNode, FoundSome)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: ZThresh
    LOGICAL :: FoundSome
    LOGICAL, POINTER :: NewBasalNode(:), ExFrontalNode(:), NewBasalElem(:)
    CHARACTER(MAX_NAME_LEN) :: FrontMaskName, BotMaskName
    !-------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    TYPE(Solver_t), POINTER :: NullSolver => NULL()
    TYPE(Element_t), POINTER :: Element, New303Elements(:,:), WorkElements(:)
    INTEGER :: i,j,k,n,dummyint, ierr, FrontBCtag, BasalBCtag, count303, &
         CountSharedExFrontal, CountSharedNewBasal, SharedExGlobal(2), &
         SharedNewGlobal(2), OldElemCount, NewElemCount
    INTEGER, POINTER :: NodeIndexes(:), AllSharedExGlobal(:)=>NULL(), &
         AllSharedNewGlobal(:)=>NULL(), FrontPerm(:), BotPerm(:)
    REAL(KIND=dp) :: Normal(3)
    LOGICAL :: ThisBC, Found, Debug
    CHARACTER(MAX_NAME_LEN) :: FuncName

    FoundSome = .FALSE.
    FuncName = "ConvertFrontalToBasal"
    Debug = .FALSE.

    n = Mesh % NumberOfNodes
    ALLOCATE(NewBasalNode(n),&
         ExFrontalNode(n),&
         FrontPerm(n),&
         BotPerm(n),&
         NewBasalElem(Mesh % NumberOfBulkElements+1: &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements))

    NewBasalNode = .FALSE.
    ExFrontalNode = .FALSE.
    NewBasalElem = .FALSE.

    CALL MakePermUsingMask( Model, NullSolver, Mesh, BotMaskName, &
         .FALSE., BotPerm, dummyint)
    CALL MakePermUsingMask( Model, NullSolver, Mesh, FrontMaskName, &
         .FALSE., FrontPerm, dummyint)

    !Find frontal BC from logical
    DO i=1,Model % NumberOfBCs
       ThisBC = ListGetLogical(Model % BCs(i) % Values,FrontMaskName,Found)
       IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
       FrontBCtag =  Model % BCs(i) % Tag
       EXIT
    END DO

    !Find basal BC from logical
    DO i=1,Model % NumberOfBCs
       ThisBC = ListGetLogical(Model % BCs(i) % Values,BotMaskName,Found)
       IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
       BasalBCtag =  Model % BCs(i) % Tag
       EXIT
    END DO

    CountSharedExFrontal = 0
    CountSharedNewBasal = 0
    SharedExGlobal = 0
    SharedNewGlobal = 0

    !---------------------------------------------------
    ! Find elements for conversion, and set node switches
    !---------------------------------------------------
    DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(i)
       IF(Element % BoundaryInfo % Constraint /= FrontBCtag) CYCLE !not on front
       IF(Element % TYPE % ElementCode == 101) CYCLE

       NodeIndexes => Element % NodeIndexes

       IF(.NOT. (ANY(BotPerm(NodeIndexes) > 0) )) CYCLE !not connected to bed

       n = Element % TYPE % NumberOfNodes

       ALLOCATE(Nodes % x(n), Nodes % y(n), Nodes % z(n))

       Nodes % x = Mesh % Nodes % x(NodeIndexes)
       Nodes % y = Mesh % Nodes % y(NodeIndexes)
       Nodes % z = Mesh % Nodes % z(NodeIndexes)

       Normal = NormalVector(Element, Nodes)

       !compare element normal to threshold
       IF(Normal(3) < ZThresh) THEN
          FoundSome = .TRUE.

          !Nodes currently on bed become 'ex frontal nodes'
          !Nodes not currently on bed become 'new basal nodes'
          DO j=1,SIZE(NodeIndexes)

             IF(BotPerm(NodeIndexes(j)) > 0) THEN
                IF(.NOT. ExFrontalNode(NodeIndexes(j))) THEN !maybe already got in another elem

                   ExFrontalNode(NodeIndexes(j)) = .TRUE.

                   !If node is in another partition, need to pass this info
                   IF(SIZE(Mesh % ParallelInfo % NeighbourList(NodeIndexes(j)) % Neighbours)>1) THEN
                      CountSharedExFrontal = CountSharedExFrontal + 1
                      IF(CountSharedExFrontal > 2) CALL Fatal(FuncName, &
                           "Found more than 2 ExFrontalNodes on partition boundary...")

                      SharedExGlobal(CountSharedExFrontal) = Mesh % ParallelInfo % GlobalDofs(NodeIndexes(j))
                   END IF
                END IF
             ELSE
                IF(.NOT. NewBasalNode(NodeIndexes(j))) THEN !maybe already got in another elem

                   NewBasalNode(NodeIndexes(j)) = .TRUE.

                   !If node is in another partition, need to pass this info
                   IF(SIZE(Mesh % ParallelInfo % NeighbourList(NodeIndexes(j)) % Neighbours)>1) THEN
                      CountSharedNewBasal = CountSharedNewBasal + 1
                      IF(CountSharedNewBasal > 2) CALL Fatal(FuncName, &
                           "Found more than 2 NewBasalNodes on partition boundary...")

                      SharedNewGlobal(CountSharedNewBasal) = &
                           Mesh % ParallelInfo % GlobalDofs(NodeIndexes(j))
                   END IF
                END IF


             END IF

          END DO

          NewBasalElem(i) = .TRUE.
          IF(Debug) PRINT *, ParEnv % MyPE, 'Debug, converting element: ',i,&
               ' with nodes: ', NodeIndexes
       END IF

       DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)
    END DO

    !Distribute information about shared frontal nodes
    !which are no longer on the front.
    !NOTE: we may also need to pass NewBasalNodes...
    IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, shared ex frontal nodes: ',SharedExGlobal
    IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, shared new basal nodes: ',SharedNewGlobal

    ALLOCATE(AllSharedExGlobal(2*ParEnv % PEs),&
         AllSharedNewGlobal(2*ParEnv % PEs))

    CALL MPI_ALLGATHER(SharedExGlobal,2,MPI_INTEGER,&
         AllSharedExGlobal,2,MPI_INTEGER, ELMER_COMM_WORLD, ierr)
    CALL MPI_ALLGATHER(SharedNewGlobal,2,MPI_INTEGER,&
         AllSharedNewGlobal,2,MPI_INTEGER, ELMER_COMM_WORLD, ierr)

    DO i=1,Mesh % NumberOfNodes
       IF(FrontPerm(i) <= 0) CYCLE
       IF(ANY(AllSharedExGlobal == Mesh % ParallelInfo % GlobalDOFs(i))) THEN
          ExFrontalNode(i) = .TRUE.
          FoundSome = .TRUE.
          IF(Debug) PRINT *, ParEnv % MyPE, ' Debug, received shared exfrontalnode: ',i
       END IF
       IF(ANY(AllSharedNewGlobal == Mesh % ParallelInfo % GlobalDOFs(i))) THEN
          NewBasalNode(i) = .TRUE.
          FoundSome = .TRUE.
          IF(Debug) PRINT *, ParEnv % MyPE, ' Debug, received shared newbasalnode: ',i
       END IF
    END DO

    !------------------------------------------------------------------------------
    ! Cycle front elements, looking for those to convert 404 -> 303 for front interp
    ! And, also, a rare case where one element is sandwiched between shared ExFrontalNodes
    ! In this case, cycle
    !------------------------------------------------------------------------------
    DO j=1,2

       count303 = 0
       DO i=Mesh % NumberOfBulkElements + 1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

          Element => Mesh % Elements(i)
          IF(Element % BoundaryInfo % Constraint /= FrontBCtag) CYCLE !not on front
          IF(Element % TYPE % ElementCode == 101) CYCLE
          IF(NewBasalElem(i)) CYCLE !element disappears from front entirely

          NodeIndexes => Element % NodeIndexes

          IF(.NOT. (ANY(BotPerm(NodeIndexes) > 0) )) CYCLE
          IF(.NOT. ANY(ExFrontalNode(NodeIndexes))) CYCLE !Not affected

          IF(j==2 .AND. Debug) PRINT *, ParEnv % MyPE, ' Debug, switching element: ',&
               i,' with nodeindexes ', NodeIndexes

          IF(COUNT(ExFrontalNode(NodeIndexes)) /= 1) CYCLE

          !iff only change one row of elements at at time, we only get here
          !through elements to the side which become 303
          count303 = count303 + 1

          !First time we just count and allocate...
          IF(j==2) THEN
             DO k=1,2
                New303Elements(count303,k) % TYPE => GetElementType( 303, .FALSE. )
                New303Elements(count303,k) % NDOFs = 3
                New303Elements(count303,k) % ElementIndex = i
                New303Elements(count303,k) % BodyID = Element % BodyID

                ALLOCATE(New303Elements(count303,k) % NodeIndexes(3))
             END DO

             !The temporary frontal element
             New303Elements(count303,1) % NodeIndexes = &
                  PACK(NodeIndexes, (.NOT. ExFrontalNode(NodeIndexes)))

             !The temporary basal element
             New303Elements(count303,2) % NodeIndexes = &
                  PACK(NodeIndexes, ( (BotPerm(NodeIndexes)>0) .OR. NewBasalNode(NodeIndexes) ) )

             DO k=1,2
                ALLOCATE(New303Elements(count303,k) % BoundaryInfo)
                New303Elements(count303,k) % BoundaryInfo % Left  => Element % BoundaryInfo % Left
                New303Elements(count303,k) % BoundaryInfo % Right => Element % BoundaryInfo % Right

                IF(k==1) THEN
                   n = FrontBCtag
                ELSE
                   n = BasalBCtag
                END IF

                New303Elements(count303,k) % BoundaryInfo % Constraint = n
             END DO

             IF(Debug) PRINT *, ParEnv % MyPE, ' debug, new frontal element ',i,' has nodes: ', &
                  New303Elements(count303,1) % NodeIndexes

             IF(Debug) PRINT *, ParEnv % MyPE, ' debug, new basal element ',i,' has nodes: ', &
                  New303Elements(count303,2) % NodeIndexes
          END IF
       END DO

       IF(j==1) THEN
          ALLOCATE(New303Elements(count303,2))
       END IF

    END DO

    !-------------------------------------------------------
    ! Now modify mesh % elements accordingly
    !-------------------------------------------------------
    IF(FoundSome) THEN

       OldElemCount = Mesh % NumberOfBulkElements + &
            Mesh % NumberOfBoundaryElements
       NewElemCount = OldElemCount + count303

       ALLOCATE(WorkElements(NewElemCount))
       WorkElements(1:OldElemCount) = Mesh % Elements(1:OldElemCount)

       DO i=1,count303
          n = New303Elements(i,1) % ElementIndex

          Element => WorkElements(n)

          CALL FreeElementStuff(Element)
 
          Element = New303Elements(i,1)
          Element => WorkElements(OldElemCount + i)

          Element = New303Elements(i,2)
          Element % ElementIndex = OldElemCount + i
       END DO

       ! Change constraint on NewBasalElem
       DO i=LBOUND(NewBasalElem,1),UBOUND(NewBasalElem,1)
          IF(.NOT. NewBasalElem(i)) CYCLE
          WorkElements(i) % BoundaryInfo % Constraint = BasalBCtag
       END DO

       DEALLOCATE(Mesh % Elements)
       Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements + count303
       Mesh % Elements => WorkElements
    END IF

    CALL SParIterAllReduceOR(FoundSome)

    NULLIFY(WorkElements)

    !TODO: Free New303Elements
    DEALLOCATE(AllSharedExGlobal, AllSharedNewGlobal, &
         NewBasalElem, FrontPerm, BotPerm, ExFrontalNode)

  END SUBROUTINE ConvertFrontalToBasal

  SUBROUTINE FreeElementStuff(Element)
    TYPE(Element_t), POINTER :: Element
    IF(ASSOCIATED(Element % NodeIndexes)) DEALLOCATE(Element % NodeIndexes)    
    IF(ASSOCIATED(Element % EdgeIndexes)) DEALLOCATE(Element % EdgeIndexes)
    IF(ASSOCIATED(Element % FaceIndexes)) DEALLOCATE(Element % FaceIndexes)
    IF(ASSOCIATED(Element % BubbleIndexes)) DEALLOCATE(Element % BubbleIndexes)
    IF(ASSOCIATED(Element % DGIndexes)) DEALLOCATE(Element % DGIndexes)
    IF(ASSOCIATED(Element % PDefs)) DEALLOCATE(Element % PDefs)
  END SUBROUTINE FreeElementStuff


  !Turns off (or back on) a specified solver, and adds a string "Save Exec When"
  ! to solver % values to allow it to be switched back on to the correct setting.
  SUBROUTINE SwitchSolverExec(Solver, Off)

    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    LOGICAL :: Off
    !-----------------------------------------
    CHARACTER(MAX_NAME_LEN) :: SaveExecWhen
    LOGICAL :: Found

    SaveExecWhen = ListGetString(Solver % Values, "Save Exec When", Found)
    IF(.NOT. Found) THEN
      SaveExecWhen = ListGetString(Solver % Values, 'Exec Solver', Found)
      IF(.NOT. Found) SaveExecWhen = 'always'
      CALL ListAddString(Solver % Values, 'Save Exec When', SaveExecWhen)
    END IF

    IF(Off) THEN

      !Turning the solver off
      Solver % SolverExecWhen = SOLVER_EXEC_NEVER
      CALL ListAddString(Solver % Values, 'Exec Solver', 'Never')

    ELSE

      CALL ListAddString(Solver % Values, 'Exec Solver', SaveExecWhen)

      SELECT CASE( SaveExecWhen )
      CASE( 'never' )
        Solver % SolverExecWhen = SOLVER_EXEC_NEVER
      CASE( 'always' )
        Solver % SolverExecWhen = SOLVER_EXEC_ALWAYS
      CASE( 'after simulation', 'after all' )
        Solver % SolverExecWhen = SOLVER_EXEC_AFTER_ALL
      CASE( 'before simulation', 'before all' )
        Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_ALL
      CASE( 'before timestep' )
        Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_TIME
      CASE( 'after timestep' )
        Solver % SolverExecWhen = SOLVER_EXEC_AFTER_TIME
      CASE( 'before saving' )
        Solver % SolverExecWhen = SOLVER_EXEC_AHEAD_SAVE
      CASE( 'after saving' )
        Solver % SolverExecWhen = SOLVER_EXEC_AFTER_SAVE
      CASE DEFAULT
        CALL Fatal("SwitchSolverExec","Programming error here...")
      END SELECT

    END IF

  END SUBROUTINE SwitchSolverExec

  SUBROUTINE PlanePointIntersection ( pp, pnorm, p1, p2, p_intersect, found_intersection )
    !Get the intersection point between a line and plane in 3D
    ! Plane defined by point "pp" and norm "pnorm", line defined by points "p1" and "p2"
    ! Intersection returned in p_intersect
    !found_intersection = .FALSE. if they happen to be parallel

    REAL(KIND=dp) :: pp(3), pnorm(3), p1(3), p2(3), p_intersect(3)
    LOGICAL :: found_intersection
    !----------------------------
    REAL(KIND=dp) :: pl(3), dist

    pl = p2 - p1

    IF(ABS(DOT_PRODUCT(pl,pnorm)) < EPSILON(1.0_dp)) THEN
      !Line and plane are parallel...
      found_intersection = .FALSE.
      RETURN
    END IF

    dist = DOT_PRODUCT((pp - p1), pnorm) / DOT_PRODUCT(pl,pnorm)

    p_intersect = p1 + dist*pl
    found_intersection = .TRUE.

  END SUBROUTINE PlanePointIntersection

  SUBROUTINE LineSegmentsIntersect ( a1, a2, b1, b2, intersect_point, does_intersect )
    ! Find if two 2D line segments intersect
    ! Line segment 'a' runs from point a1 => a2, same for b

    IMPLICIT NONE

    REAL(KIND=dp) :: a1(2), a2(2), b1(2), b2(2), intersect_point(2)
    LOGICAL :: does_intersect
    !-----------------------
    REAL(KIND=dp) :: r(2), s(2), rxs, bma(2), t, u


    does_intersect = .FALSE.
    intersect_point = 0.0_dp

    r = a2 - a1
    s = b2 - b1

    rxs = VecCross2D(r,s)

    IF(rxs == 0.0_dp) RETURN

    bma = b1 - a1

    t = VecCross2D(bma,s) / rxs
    u = VecCross2D(bma,r) / rxs

    IF(t < 0.0_dp .OR. t > 1.0_dp .OR. u < 0.0_dp .OR. u > 1.0_dp) RETURN

    intersect_point = a1 + (t * r)
    does_intersect = .TRUE.

  END SUBROUTINE LineSegmentsIntersect

  SUBROUTINE LinesIntersect ( a1, a2, b1, b2, intersect_point, does_intersect )
    ! Find where two 2D lines intersect
    ! Line 'a' explicitly defined by points a1, a2 which lie on line, same for b
    ! based on LineSegmentsIntersect above

    IMPLICIT NONE

    REAL(KIND=dp) :: a1(2), a2(2), b1(2), b2(2), intersect_point(2)
    LOGICAL :: does_intersect
    !-----------------------
    REAL(KIND=dp) :: r(2), s(2), rxs, bma(2), t, u


    does_intersect = .TRUE.

    intersect_point = 0.0_dp

    r = a2 - a1
    s = b2 - b1

    rxs = VecCross2D(r,s)

    IF(rxs == 0.0_dp) THEN
      does_intersect = .FALSE.
      RETURN
    ENDIF

    bma = b1 - a1

    t = VecCross2D(bma,s) / rxs
    u = VecCross2D(bma,r) / rxs

    intersect_point = a1 + (t * r)

  END SUBROUTINE LinesIntersect

  FUNCTION VecCross2D(a, b) RESULT (c)
    REAL(KIND=dp) :: a(2), b(2), c

    c = a(1)*b(2) - a(2)*b(1)

  END FUNCTION VecCross2D

  !This subroutine should identify discrete calving events for the
  !purposes of local remeshing. For now it returns 1
  SUBROUTINE CountCalvingEvents(Model, Mesh,CCount)
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),POINTER :: Mesh
    INTEGER :: CCount

    Ccount = 1
  END SUBROUTINE CountCalvingEvents
END MODULE CalvingGeometry

