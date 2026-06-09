!/******************************************************************************
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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh refinement and splitting: equal refinement, quad splitting, level-set splitting.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshSplit

    USE MeshBasics
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

  FUNCTION SplitMeshEqual(Mesh,h) RESULT( NewMesh )
!------------------------------------------------------------------------------
    REAL(KIND=dp), OPTIONAL :: h(:)
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:),xh(:)
    INTEGER :: i, j, k, n, NewElCnt, NodeCnt, EdgeCnt, FaceCnt, Node, ParentId, Diag, NodeIt
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Eparent,Face,Faces(:)
    INTEGER, POINTER :: Child(:,:)
    INTEGER :: n1,n2,n3,EoldNodes(4),FaceNodes(4),EdgeNodes(2) ! Only linears so far
    INTEGER :: FaceNumber,Edge1,Edge2,Edge3,Edge4,Node12,Node23,Node34,Node41,Node31
    REAL(KIND=dp) :: dxyz(3,3),Dist(3),r,s,t,h1,h2
    TYPE(PElementDefs_t), POINTER :: PDefs
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER, ALLOCATABLE :: FacePerm(:), BulkPerm(:)
    LOGICAL :: Parallel
    CHARACTER(*), PARAMETER :: Caller = 'SplitMeshEqual'
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN

    CALL Info( Caller, 'Mesh splitting works for first order elements 303, 404, 504, (706) and 808.', Level = 6 )

    DO i=1,Mesh % NumberOfBulkElements
      SELECT CASE(Mesh % Elements(i) % TYPE % ElementCode/100)
      CASE(6)
        CALL Fatal(Caller,'Pyramids not supported, sorry.')
      END SELECT
    END DO

    NewMesh => AllocateMesh()

    NewMesh % SingleMesh = Mesh % SingleMesh
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. NewMesh % SingleMesh )

    
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT.EdgesPresent) CALL FindMeshEdges( Mesh )

    CALL ResetTimer(Caller)

    CALL Info( Caller, '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( Caller, Message, Level=6 )
!
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = Mesh % NumberOfNodes + Mesh % NumberOfEdges
!
!   For quad faces add one node in the center:
!   ------------------------
    ALLOCATE(FacePerm(Mesh % NumberOfFaces)); FacePerm = 0
    FaceCnt = 0
    DO i = 1, Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       IF( Face % TYPE % NumberOfNodes == 4 ) THEN
         NodeCnt = NodeCnt+1
         FaceCnt = FaceCnt+1
         FacePerm(i) = NodeCnt
       END IF
    END DO    
    IF(FaceCnt>0) CALL Info( Caller,'Added '//I2S(FaceCnt)//' nodes in the center of faces',Level=10)

!
!   For quads and bricks, count centerpoints:
!   -----------------------------------------
    NodeIt = 0
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(4,8)
          NodeCnt = NodeCnt + 1
          NodeIt = NodeIt + 1
       END SELECT
    END DO    
    IF(NodeIt>0) CALL Info( Caller,'Added '//I2S(NodeIt)//' nodes in the center of bulks',Level=10)

!
!   new mesh nodecoordinate arrays:
!   -------------------------------
    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:Mesh % NumberOfNodes) = u(1:Mesh % NumberOfNodes)
    y(1:Mesh % NumberOfNodes) = v(1:Mesh % NumberOfNodes)
    z(1:Mesh % NumberOfNodes) = w(1:Mesh % NumberOfNodes)

! what is h? - pointer to nodal element size
    IF (PRESENT(h)) THEN
      ALLOCATE(xh(SIZE(x)))
      xh(1:SIZE(h)) = h
    END IF
!
!   add edge centers:
!   -----------------
    j =  Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
       j = j + 1
       Edge => Mesh % Edges(i)
       k = Edge % TYPE % NumberOfNodes
       IF (PRESENT(h)) THEN
         h1=h(Edge % NodeIndexes(1))
         h2=h(Edge % NodeIndexes(2))
         r=1._dp/(1+h1/h2)
         x(j) = r*u(Edge%NodeIndexes(1))+(1-r)*u(Edge%NodeIndexes(2))
         y(j) = r*v(Edge%NodeIndexes(1))+(1-r)*v(Edge%NodeIndexes(2))
         z(j) = r*w(Edge%NodeIndexes(1))+(1-r)*w(Edge%NodeIndexes(2))
         xh(j)=r*h1+(1-r)*h2
       ELSE
         x(j) = SUM(u(Edge % NodeIndexes))/k
         y(j) = SUM(v(Edge % NodeIndexes))/k
         z(j) = SUM(w(Edge % NodeIndexes))/k
       END IF
    END DO    
    CALL Info(Caller,'Added edge centers to the nodes list.', Level=15 )  

!   add quad face centers for bricks and prisms(wedges):
!   ----------------------------
    j = Mesh % NumberOfNodes + Mesh % NumberOfEdges
    DO i=1,Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       k = Face % TYPE % NumberOfNodes
       IF( k==4 ) THEN
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Face % EdgeIndexes(2))
            h2=xh(n+Face % EdgeIndexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Face % EdgeIndexes(3))
            h2=xh(n+Face % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Face,u(Face % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Face,v(Face % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Face,w(Face % NodeIndexes),r,s)
            xh(j) = InterpolateInElement2D(Face,h(Face % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Face % NodeIndexes))/k
            y(j) = SUM(v(Face % NodeIndexes))/k
            z(j) = SUM(w(Face % NodeIndexes))/k
          END IF
       END IF
    END DO    
    CALL Info(Caller,'Added face centers to the nodes list.', Level=15 )

!   add centerpoint for quads & bricks:
!   -----------------------------------
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       k = Eold % TYPE % NumberOfNodes
       SELECT CASE( Eold % TYPE % ElementCode / 100 )

       CASE(4)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Eold % Edgeindexes(2))
            h2=xh(n+Eold % Edgeindexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Eold % EdgeIndexes(3))
            h2=xh(n+Eold % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Eold,u(Eold % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Eold,v(Eold % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Eold,w(Eold % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       CASE(8)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes+Mesh % NumberOfEdges
            h1=xh(n+Eold % FaceIndexes(4))
            h2=xh(n+Eold % FaceIndexes(6))
            r=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(5))
            h2=xh(n+Eold % FaceIndexes(3))
            s=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(2))
            h2=xh(n+Eold % FaceIndexes(1))
            t=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement3D(Eold,u(Eold % NodeIndexes),r,s,t)
            y(j) = InterpolateInElement3D(Eold,v(Eold % NodeIndexes),r,s,t)
            z(j) = InterpolateInElement3D(Eold,w(Eold % NodeIndexes),r,s,t)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       END SELECT
    END DO
    CALL Info(Caller,'Added quad and brick centers to the nodes list.', Level=15 )

    
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt
!
!   Update bulk elements:
!   =====================
!
!   First count new elements:
!   -------------------------
    NewElCnt = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode/100 )

!      Each element will be divided into 2**Dim new elements:
!      ------------------------------------------------------
       CASE(1)
          NewElCnt = NewElCnt + 1 ! lines
       CASE(2)
          NewElCnt = NewElCnt + 2 ! lines
       CASE(3)
          NewElCnt = NewElCnt + 4 ! trias
       CASE(4)
          NewElCnt = NewElCnt + 4 ! quads
       CASE(5)
          NewElCnt = NewElCnt + 8 ! tetras
       CASE(7)
          NewElCnt = NewElCnt + 8 ! prisms (wedges)
       CASE(8)
          NewElCnt = NewElCnt + 8 ! hexas
       END SELECT
    END DO

    WRITE( Message, * ) 'Count of new elements : ', NewElCnt
    CALL Info( Caller, Message, Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info(Caller,'New mesh allocated.', Level=10 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 8 )
    CALL Info(Caller,'Array for bulk elements allocated.', Level=10 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes
    EdgeCnt = Mesh % NumberOfEdges

!
!   Index to old edge/quad/hexa centerpoint node in the new mesh nodal arrays:
!   ---------------------------------------------------------------------
    Node = NodeCnt + EdgeCnt + FaceCnt
!
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)

       SELECT CASE( Eold % TYPE % ElementCode )
       CASE(101)
!
!         Copy point element
!         ------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 1)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)

       CASE(202)
!
!         Split edge to two edges
!         ------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 2)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 2)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)

       CASE(303)
!
!         Split triangle to four triangles from
!         edge centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,3) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,4) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt

       CASE(404)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split quad to four new quads from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)


       CASE(504)
!
!         Split tetra to 8 new elements from
!         corners and edge centerpoints:
!         ----------------------------------
!
!         1st new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt

!         Then the annoying part; we still have to split the
!         remaining octahedron into four elements. This can
!         be done in three ways of which only one preserves
!         the minimum angle condition (Delaunay splitting):
!         --------------------------------------------------
          dxyz(1,1) = x(Eold % EdgeIndexes(4) + NodeCnt) &
                    - x(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(2,1) = y(Eold % EdgeIndexes(4) + NodeCnt) &
                    - y(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(3,1) = z(Eold % EdgeIndexes(4) + NodeCnt) &
                    - z(Eold % EdgeIndexes(2) + NodeCnt)

          dxyz(1,2) = x(Eold % EdgeIndexes(5) + NodeCnt) &
                    - x(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(2,2) = y(Eold % EdgeIndexes(5) + NodeCnt) &
                    - y(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(3,2) = z(Eold % EdgeIndexes(5) + NodeCnt) &
                    - z(Eold % EdgeIndexes(3) + NodeCnt)

          dxyz(1,3) = x(Eold % EdgeIndexes(6) + NodeCnt) &
                    - x(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(2,3) = y(Eold % EdgeIndexes(6) + NodeCnt) &
                    - y(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(3,3) = z(Eold % EdgeIndexes(6) + NodeCnt) &
                    - z(Eold % EdgeIndexes(1) + NodeCnt)

          Dist(1) = SQRT( dxyz(1,1)**2 + dxyz(2,1)**2 + dxyz(3,1)**2 )
          Dist(2) = SQRT( dxyz(1,2)**2 + dxyz(2,2)**2 + dxyz(3,2)**2 )
          Dist(3) = SQRT( dxyz(1,3)**2 + dxyz(2,3)**2 + dxyz(3,3)**2 )

          Diag = 1  ! The default diagonal for splitting is between edges 2-4
          IF (Dist(2) < Dist(1) .AND. Dist(2) < Dist(3)) Diag = 2 ! Edges 3-5
          IF (Dist(3) < Dist(1) .AND. Dist(3) < Dist(2)) Diag = 3 ! Edges 1-6

          SELECT CASE( Diag )
          CASE(1)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
          CASE(2)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
          CASE(3)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt

          END SELECT


       CASE(706)
!
!         Split prism to 8 new prism from edge
!         centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt 
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt 
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt 
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))

!
!         3rd new element (near node 3)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(9) + NodeCnt

!
!         4th new element (bottom center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt

!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt

!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
!
!         8th new element (top half, center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt



       CASE(808)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split brick to 8 new bricks from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 8)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(7) = Node
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(6))
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(8) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(6) = Node
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(12)+ NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(5) = Node
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(5))
!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(8) + NodeCnt
!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Node
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(2))
!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(12)+ NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(8) = Eold % NodeIndexes(8)
!
!         8th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(7) = Eold % NodeIndexes(7)
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(7) + NodeCnt

       CASE DEFAULT
          WRITE( Message,* ) 'Element type ', Eold % TYPE % ElementCode, &
              ' not supported by the multigrid solver.'
          CALL Fatal( Caller, Message )
       END SELECT
    END DO

!
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt

!
!   Update boundary elements:
!   NOTE: Internal boundaries not taken care of...:!!!!
!   ---------------------------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements

       j = i + Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(j)
!
!      get parent of the boundary element:
!      -----------------------------------
       Eparent => Eold % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED(Eparent) ) &
          eParent => Eold % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED( Eparent ) ) CYCLE

       ParentId = Eparent % ElementIndex

       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(1)
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 1 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
          DO j=NewElCnt-1,1,-1
            Eold => NewMesh % Elements(j)
            IF(Eold % Type % ElementCode/100==2) THEN
              IF(ANY(Eold % NodeIndexes==Enew % NodeIndexes(1))) THEN
                IF(.NOT.ASSOCIATED(Enew % BoundaryInfo % Left)) THEN
                  Enew % BoundaryInfo % Left => Eold
                ELSE
                  Enew % BoundaryInfo % Right => Eold
                  EXIT
                END IF
              END IF
            END IF
          END DO

       CASE(2)
!
!         Line segments:
!         ==============
!
!         which edge of the parent element are we ?
!         -----------------------------------------
          Found = .FALSE.
          DO Edge1=1,SIZE(Eparent % EdgeIndexes)
            Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
            Found = ANY(Eold % NodeIndexes(1:2) == Edge % NodeIndexes(1) ) .AND. &
                ANY(Eold % NodeIndexes(1:2) == Edge % NodeIndexes(2) )
            IF(Found) EXIT
          END DO
          IF(.NOT. Found) THEN
            CALL Fatal(Caller,'Could not find parent edge with nodes: '//&
                I2S(Eold % NodeIndexes(1))//' '//I2S(Eold % NodeIndexes(2)))
          END IF

!
!         index of the old edge centerpoint in the
!         new mesh nodal arrays:
!         ----------------------------------------
          Node = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------

          Found = .FALSE.

          n1 = 4 
          IF( Eparent % TYPE % ElementCode > 500 ) n1 = 8
          
          DO j=1,n1
            Eptr => NewMesh % Elements( Child(ParentId,j) )
            n = Eptr % TYPE % NumberOfNodes

            ! The parent is unique! Hence it is enough to find a parent with both matches.
            Found =  ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(1) ) .AND. &
                ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(2) )
            IF ( Found ) EXIT
          END DO


          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
                    
          DO j=1,n1
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found =  ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(1) ) .AND. &
                 ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(2) )
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr

       CASE(3)
!
!         Trias:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:3) = Eold % NodeIndexes(1:3)
          CALL sort( 3, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:3) = Face % NodeIndexes(1:3)
             CALL sort( 3, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) ) EXIT

          END DO
!
!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node31 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node31
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr

       CASE(4)
!
!         Quads:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:4) = Eold % NodeIndexes(1:4)
          CALL sort( 4, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:4) = Face % NodeIndexes(1:4)
             CALL sort( 4, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) .AND. &
                  EoldNodes(4) == FaceNodes(4) ) EXIT

          END DO

!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Fourth edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          DO Edge4 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge4) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node = FacePerm(Eparent % FaceIndexes(FaceNumber)) ! faces mid-point
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node34 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
          Node41 = Eparent % EdgeIndexes(Edge4) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Node41
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n) ) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 )  CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          Enew % NodeIndexes(4) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node41
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Node34
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Node34
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
       END SELECT
    END DO

!
!   Update new mesh boundary element counts:
!   ----------------------------------------
    NewMesh % NumberOfBoundaryElements = NewElCnt - &
            NewMesh % NumberOfBulkElements
    NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
    NewMesh % MaxElementNodes = Mesh % MaxElementNodes

    j = 0
    DO i=1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
      Enew => NewMesh % Elements(i)

      IF ( Enew % DGDOFs>0 ) THEN
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF

      IF (i<=NewMesh % NumberOfBulkElements) THEN
         PDefs => Enew % PDefs

         IF(ASSOCIATED(PDefs)) THEN
           CALL AllocatePDefinitions(Enew)
           Enew % PDefs = PDefs

           ! All elements in actual mesh are not edges
           Enew % PDefs % isEdge = .FALSE.

           ! If element is of type tetrahedron and is a p element,
           ! do the Ainsworth & Coyle trick
           IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
            CALL GetRefPElementNodes( Enew % Type,  Enew % Type % NodeU, &
                 Enew % Type % NodeV, Enew % Type % NodeW )
         END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO

    CALL Info( Caller, '******** New mesh ********', Level=6 )
    WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
    CALL Info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
    CALL Info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
    CALL Info( Caller, Message, Level=6 )


    ! Information of the new system size, also in parallel
    !----------------------------------------------------------------------
    ParTmp(1) = Mesh % NumberOfNodes
    ParTmp(2) = Mesh % NumberOfBulkElements
    ParTmp(3) = Mesh % NumberOfBoundaryElements
    ParTmp(4) = NewMesh % NumberOfNodes
    ParTmp(5) = NewMesh % NumberOfBulkElements
    ParTmp(6) = NewMesh % NumberOfBoundaryElements

    IF( .FALSE. .AND. Parallel ) THEN
      CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

      CALL Info(Caller,'Information on parallel mesh sizes')
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' bulk elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(3),' boundary elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(4),' nodes'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(5),' bulk elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(6),' boundary elements'
      CALL Info(Caller,Message)
    END IF


    CALL CheckTimer(Caller,Delete=.TRUE.)

!
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    IF( Parallel ) THEN
      CALL UpdateParallelMesh( Mesh, NewMesh )
    END IF
!
!
!   Finalize:
!   ---------
    DEALLOCATE( Child )
    IF(.NOT.EdgesPresent) THEN
      CALL Info(Caller,'Releasing edges from the old mesh as they are not needed!',Level=20)
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    ELSE
      CALL Info(Caller,'Generating edges in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshEdges( NewMesh )
    END IF

    ! Our boundary may be a circle, cylinder or sphere surface.
    ! Honor those shapes when splitting the mesh!
    CALL FollowCurvedBoundary( CurrentModel, NewMesh, .FALSE. ) 
    
    
!call writemeshtodisk( NewMesh, "." )
!stop
CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelMesh( Mesh, NewMesh )
!------------------------------------------------------------------------------
       TYPE(Mesh_t), TARGET :: Mesh, NewMesh
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Edge, Face, Element, BoundaryElement
       INTEGER :: i,j,k,l,m,n,p,q, istat
       INTEGER, POINTER :: IntCnts(:),IntArray(:),Reorder(:)
       INTEGER, ALLOCATABLE :: list1(:), list2(:)
       LOGICAL, ALLOCATABLE :: InterfaceTag(:)

       INTEGER :: jedges
       LOGICAL :: Found
!------------------------------------------------------------------------------


!      Update mesh interfaces for parallel execution.
!      ==============================================
!
!      Try to get an agreement about the  global numbering
!      of new mesh nodes among set of processes solving
!      this specific eq. Also allocate and generate
!      all other control information needed in parallel
!      execution:
!      ----------------------------------------------------
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) &
         CALL Fatal( 'UpdateParallelMesh', 'Allocate error.' )
       CALL AllocateVector( NewMesh % ParallelInfo % GInterface,n  )
       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )

       DO i=1,n
          NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % GInterface = .FALSE.
       NewMesh % ParallelInfo % GInterface(1:n) = Mesh % ParallelInfo % GInterface

       NewMesh % ParallelInfo % GlobalDOFs = 0
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = &
          Mesh % ParallelInfo % GlobalDOFs
!
!      My theory is, that a new node will be an
!      interface node only if all the edge or face
!      nodes which contribute to its existence are
!      interface nodes (the code immediately below
!      will only count sizes):
!      -------------------------------------------
!

       ! New version based on edges and faces (2. March 2007):
       !=====================================================
       SELECT CASE( CoordinateSystemDimension() )
          
       CASE(2)
          !
          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % GInterface(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'
          !
          ! Determine possible interface edges:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfEdges ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfEdges
             Edge => Mesh % Edges(i)
             IF( ASSOCIATED(Edge % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Edge % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % GInterface( Edge % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          !
          ! Eliminate false positives based on BoundaryElement -data:
          !----------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements( Mesh % NumberOfBulkElements + i )
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED( Element ) ) &
                  Element => BoundaryElement % BoundaryInfo % Right
             IF( .NOT.ASSOCIATED( Element ) ) CYCLE
             IF( .NOT.ASSOCIATED( Element % EdgeIndexes ) ) CYCLE
             
             ALLOCATE( list1( SIZE( BoundaryElement % NodeIndexes )))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort( SIZE(list1), list1 )
             
             DO j = 1,Element % TYPE % NumberOfEdges
                k = Element % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                IF( SIZE( Edge % NodeIndexes ) /= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Edge % NodeIndexes )))
                list2 = Edge % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO

                DEALLOCATE(list2)
                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfEdges
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Edge => Mesh % Edges(i)
             
             ! This is just for the edge count:
             !---------------------------------
             IF( NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + i) ) CYCLE
             
             ! Mark interface nodes and count edges:
             !--------------------------------------
             NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + i) = .TRUE.
             p = p+1

          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 2*p ! check
          
       CASE(3)

          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % GInterface(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'

          ! Determine possible interface faces:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfFaces ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( ASSOCIATED(Face % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Face % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % GInterface( Face % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          
          ! Eliminate false interface faces based on BoundaryElement -data:
          !----------------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED(Element) ) &
                Element => BoundaryElement % BoundaryInfo % Right
              IF( .NOT.ASSOCIATED(Element) ) CYCLE
              IF( .NOT.ASSOCIATED(Element % FaceIndexes) ) CYCLE
             
             ALLOCATE(list1(SIZE(BoundaryElement % NodeIndexes)))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort(SIZE(list1),list1)
             
             DO j = 1,Element % TYPE % NumberOfFaces
                k = Element % FaceIndexes(j)
                Face => Mesh % Faces(k)
                IF(SIZE(Face % NodeIndexes)/= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Face % NodeIndexes )))
                list2 = Face % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO
                
                DEALLOCATE(list2)

                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Count interface faces:
          !-----------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( InterfaceTag(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface faces'
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Face => Mesh % Faces(i)
             
             DO j = 1,SIZE( Face % EdgeIndexes )
                k = Face % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                
                ! This is just for the edge count:
                !---------------------------------
                IF( NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + k) ) CYCLE
                
                ! Mark interface nodes and count edges:
                !--------------------------------------
                NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + k) = .TRUE.
                p = p+1
             END DO
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 3*p ! check
          
       END SELECT

!======================================================================================================
       j = p
       jedges = p

!      For bricks, check also the faces:
!      ---------------------------------
       DO i = 1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i) 
          IF( Face % TYPE % NumberOfNodes == 4 ) THEN
             IF ( ALL( Mesh % ParallelInfo % GInterface( Face % NodeIndexes ) ) ) THEN
                NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes &
                     + Mesh % NumberOfEdges + i ) = .TRUE.
                j = j + 1
                k = k + Face % TYPE % NumberOfNodes
             END IF
          END IF
       END DO

!      CALL AllocateVector( IntCnts,  j )
!      CALL AllocateVector( IntArray, k )
!
!      Old mesh nodes were copied as is...
!
       IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % Neighbourlist ) ) THEN
         CALL Fatal('UpdateParallelMesh','Original mesh has no NeighbourList!')
       END IF
       
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
           CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, 1 )
           NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = ParEnv % MyPe
         ELSE
           CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, &
               SIZE( Mesh % ParallelInfo % Neighbourlist(i) % Neighbours) )
           NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
               Mesh % ParallelInfo % NeighbourList(i) % Neighbours
         END  IF
       END DO
!
!      Take care of the new mesh internal nodes.
!      Parallel global numbering will take care
!      of the interface nodes:
!      ----------------------------------------
       DO i=Mesh % NumberOfNodes+1, NewMesh % NumberOfNodes
          IF ( .NOT. NewMesh % ParallelInfo % GInterface(i) ) THEN
            CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours,1 )
            NewMesh % ParallelInfo % NeighbourList(i) %  Neighbours(1) = ParEnv % MyPE
          END IF
       END DO
!
!      Copy global indices of edge and/or face nodes
!      to temporary work arrays:
!      ---------------------------------------------
!
! check also this:
!      j = 0
!      k = 0
!      DO i = 1,Mesh % NumberOfEdges
!         Edge => Mesh % Edges(i)
!         
!         ! Added check for parent elements 25.2.2007:
!         Found = .NOT.( ASSOCIATED(edge % boundaryinfo % left) &
!              .AND.  ASSOCIATED(edge % boundaryinfo % right) )
!         
!         IF ( ALL(Mesh % ParallelInfo % GInterface(Edge % NodeIndexes)) .AND. Found ) THEN
!            j = j + 1
!            IntCnts(j) = Edge % TYPE % NumberOfNodes
!            IntArray( k+1:k+IntCnts(j) ) = &
!                 Mesh % Parallelinfo % GlobalDOFs(Edge % NodeIndexes)
!            CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!            k = k + IntCnts(j)
!         END IF
!      END DO
!      !
!      ! For bricks, check also the faces:
!      ! ---------------------------------
!      DO i = 1,Mesh % NumberOfFaces
!         Face => Mesh % Faces(i)
!         IF( Face % TYPE % NumberOfNodes == 4 ) THEN
!            IF ( ALL( Mesh % ParallelInfo % GInterface(Face % NodeIndexes) ) ) THEN
!               j = j + 1
!               IntCnts(j) = Face % TYPE % NumberOfNodes
!               IntArray(k+1:k+IntCnts(j)) = &
!                    Mesh % ParallelInfo % GlobalDOFs(Face % NodeIndexes)
!               CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!               k = k + IntCnts(j)
!            END IF
!         END IF
!      END DO
!
!      Finally the beef, do the exchange of new
!      interfaces. The parallel global numbering
!      subroutine will also do reordering of the
!      nodes, hence the reorder array:
!      -------------------------------------------
       CALL AllocateVector( Reorder, NewMesh % NumberOfNodes )
       Reorder = [ (i, i=1,NewMesh % NumberOfNodes) ]

       k = NewMesh % Nodes % NumberOfNodes - Mesh % Nodes % NumberOfNodes


       CALL ResetTimer('ParallelGlobalNumbering')
       CALL ParallelGlobalNumbering( NewMesh, Mesh, k, Reorder )
       CALL CheckTimer('ParallelGlobalNumbering',Level=7,Delete=.TRUE.)

       
!      Account for the reordering of the nodes:
!      ----------------------------------------
       DO i=1,NewMesh % NumberOfBulkElements + &
           NewMesh % NumberOfBoundaryElements
         NewMesh % Elements(i) % NodeIndexes = &
             Reorder( NewMesh % Elements(i) % NodeIndexes )
       END DO

!      DEALLOCATE( IntCnts, IntArray, Reorder )
       !      DEALLOCATE( Reorder )
       

!------------------------------------------------------------------------------
    END SUBROUTINE UpdateParallelMesh
  END FUNCTION SplitMeshEqual

  SUBROUTINE SplitMeshQuads(Mesh, Vlist) 
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    INTEGER :: i, j, k, k2, n, AddCnt, NewElCnt, nBulkElems, nBoundaryElems 
    LOGICAL :: Found, FacesPresent, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Parent
    TYPE(PElementDefs_t), POINTER :: PDefs
    LOGICAL :: Parallel, IsBulkElement, IsOddCut
    LOGICAL :: Is15, Is24, Is35, Is26, Is34, Is16
    INTEGER, ALLOCATABLE :: CutCorner(:), MinCorner(:), BulkElementOffset(:)
    INTEGER :: TypeCnt(0:8), LocalMap(4), LocalMin(3), LocalCut(3)
    INTEGER :: CutChanges, MeshDim, CutComb, iCut
    TYPE(Element_t), POINTER :: Face
    TYPE(Element_t), POINTER :: NewElements(:) => NULL()
    CHARACTER(*), PARAMETER :: Caller="SplitMeshQuads"

!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN
    IF(.NOT. ASSOCIATED(VList)) RETURN
    IF(Mesh % MeshDim < 2 ) RETURN
    
    IF( Mesh % MeshDim == 2 ) THEN
      IF(.NOT. ListGetLogical( Vlist,'Split Mesh Quads',Found ) ) RETURN              
    ELSE
      IF(.NOT. ListGetLogical( Vlist,'Split Mesh Prisms',Found ) ) RETURN              
    END IF
    CALL Info( Caller,'Splitting all quadrilaterals into triangles in '//I2S(Mesh % MeshDim)//'D',Level=5)
      
    TypeCnt = 0
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      j = Mesh % Elements(i) % TYPE % ElementCode/100
      TypeCnt(j) = TypeCnt(j) + 1
    END DO
    
    IF(TypeCnt(4) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(4))//' quad elements',Level=8)
    IF(TypeCnt(6) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(6))//' pyramid elements',Level=8)
    IF(TypeCnt(7) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(7))//' prism elements',Level=8)    
    IF(TypeCnt(8) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(8))//' hexahedron elements',Level=8)

    !DO i=0,8
    !  PRINT *,'TypeCount:',i,TypeCnt(i)
    !END DO
    
    IF(TypeCnt(6) + TypeCnt(8) > 0 ) THEN
      CALL Fatal(Caller,'Not implemented yet for pyramids and hexahedrons!')
    END IF
    
    IF(Mesh % MeshDim == 3 .AND. TypeCnt(7) == 0) THEN
      CALL Warn(Caller,'No wedges exist, doing nothing!')
      RETURN
    END IF
    
    CALL ResetTimer(Caller)

    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh )
    
    AddCnt = TypeCnt(4) + 2*TypeCnt(7)    
    CALL Info(Caller,'Number of elements added by splitting is '//I2S(AddCnt),Level=6)
    
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
    MeshDim = Mesh % MeshDim
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    FacesPresent = ASSOCIATED(Mesh % Faces)
        
    NewElCnt = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements + AddCnt
    CALL Info(Caller,'Count of new elements: '//I2S(NewElCnt),Level=7)
    CALL AllocateVector( NewElements, NewElCnt )
    CALL Info(Caller,'New elements allocated.', Level=10 )
    
    ! We do not need to all the complex stuff in 2D.
    IF( MeshDim < 3 ) GOTO 1

    ! First negotiate the split direction of the mesh.
    ! We start from assuming that cut direction is though off local
    ! indexes and the one owning the smallest global index calls the shots.    
    !----------------------------------------------------------------------
    CALL Info(Caller,'Trying to find consistent splitting direction in 3D',Level=12)
    IF(.NOT.FacesPresent) CALL FindMeshFaces3D( Mesh )
    ALLOCATE(CutCorner(Mesh % NumberOfFaces), MinCorner(Mesh % NumberOfFaces))
    
    ! Initialize with smallest index of the face
    DO i=1,Mesh % NumberOfFaces
      Face => Mesh % Faces(i)
      MinCorner(i) = MINVAL(Face % NodeIndexes)
    END DO
    CutCorner = MinCorner

    BLOCK
      INTEGER :: DoMax
      INTEGER, POINTER :: Inds(:)
      REAL(KIND=dp) :: d13,d24
    
      DoMax = 0
      IF(LIstGetLogical(Vlist,'Split Mesh Prisms Min',Found )) DoMax = 1
      IF(LIstGetLogical(Vlist,'Split Mesh Prisms Max',Found )) DoMax = -1

      ! Optionally cut the 3D meshes such that the shorter (longer) diagonal
      ! is used to cut the quad faces. 
      IF(DoMax /= 0) THEN
        DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF(Face % TYPE % ElementCode /= 404) CYCLE
          Inds => Face % NodeIndexes
          
          ! Compute |r1-r3|^2 
          d13 = (x(Inds(1))-x(Inds(3)))**2 + &
              (y(Inds(1))-y(Inds(3)))**2 + (z(Inds(1))-z(Inds(3)))**2 
          
          ! Compute |r2-r4|^2 
          d24 = (x(Inds(2))-x(Inds(4)))**2 + &
              (y(Inds(2))-y(Inds(4)))**2 + (z(Inds(2))-z(Inds(4)))**2 
          
          IF(DoMax * d13 < DoMax * d24) THEN
            CutCorner(i) = MIN(Inds(1),Inds(3))
          ELSE
            CutCorner(i) = MIN(Inds(2),Inds(4))
          END IF
        END DO
      END IF
    END BLOCK
          
    
    
    DO WHILE(.TRUE.)
      CutChanges = 0
      
      DO i=1,Mesh % NumberOfBulkElements         
        Eold => Mesh % Elements(i)        
        SELECT CASE( Eold % TYPE % ElementCode )
          
        CASE(706)
          ! Faces 1 & 2 are triangles
          ! Faces 3, 4 and 5 are married
          ! There are 8 ways to cut the faces but only 6 of those are legal. 
          ! WedgeFaceMap(3,:) = (/ 1,2,5,4 /)
          ! WedgeFaceMap(4,:) = (/ 2,3,6,5 /)
          ! WedgeFaceMap(5,:) = (/ 3,1,4,6 /)
                            
          LocalMin(1:3) = MinCorner(Eold % FaceIndexes(3:5))
          LocalCut(1:3) = CutCorner(Eold % FaceIndexes(3:5))

          ! How are we cutting the three faces? 
          Is24 = ANY( LocalCut(1) == Eold % NodeIndexes([2,4])) 
          Is26 = ANY( LocalCut(2) == Eold % NodeIndexes([2,6])) 
          Is16 = ANY( LocalCut(3) == Eold % NodeIndexes([1,4])) 
          Is15 = .NOT. Is24
          Is35 = .NOT. Is26
          Is34 = .NOT. Is16

          ! Code the cases to numbers [111,222]
          ! Cases 121 and 212 are not allowed!
          CutComb = 111
          IF(.NOT. Is24) CutComb = CutComb+100  ! Is15
          IF(.NOT. Is26) CutComb = CutComb+10   ! Is35
          IF(.NOT. Is16) CutComb = CutComb+1    ! Is34
          
          IF( CutComb == 121 .OR. CutComb == 212 ) THEN
            iCut = 0
            DO k=1,3
              IF(LocalCut(k) == MAXVAL(LocalCut(1:3))) EXIT
            END DO

            IF(k==1) THEN
              IF(CutComb == 121 ) iCut = 1 !-> 221 
              IF(CutComb == 212 ) iCut = 2 !-> 112
            ELSE IF(k==2) THEN
              IF(CutComb == 121 ) iCut = 2 !-> 111
              IF(CutComb == 212 ) iCut = 3 !-> 222
            ELSE IF(k==3) THEN
              IF(CutComb == 121 ) iCut = 3 !-> 122
              IF(CutComb == 212 ) iCut = 1 !-> 211
            ELSE
              CALL Fatal(Caller,'Invalid value for k!')
            END IF

            IF(icut>0) THEN
              CutCorner(Eold % FaceIndexes(2+k)) = Eold % NodeIndexes(iCut)
              CutChanges = CutChanges + 1                            
            ELSE
              CALL Fatal(Caller,'Could not fix invalid cut!')
            END IF
          END IF
        END SELECT
      END DO

      CALL Info(Caller,'Number of switches in cut direction: '//I2S(CutChanges))
      IF(CutChanges == 0) EXIT      
    END DO
    
    ! Jump directly here if we have 2D mesh. 
1   CONTINUE

    ! We need to register the offset coming from split.
    ALLOCATE(BulkElementOffset(Mesh % NumberOfBulkElements))
    BulkElementOffset = 0
    NewElCnt = 0
            
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Eold => Mesh % Elements(i)
      IsBulkElement = (i <= Mesh % NumberOfBulkElements ) 
      
      IF(IsBulkElement) THEN
        ! For bulk elements store the offset so we can remap the parents more easily. 
        BulkElementOffset(i) = NewElCnt 
      END IF

      SELECT CASE( Eold % TYPE % ElementCode )

      CASE(101,202,303,504)
        ! Copy elements without quad faces as is. 
        NewElCnt = NewElCnt + 1
        Enew => NewElements(NewElCnt)
        Enew = Eold
        Enew % ElementIndex = NewElCnt         
        CALL AllocateVector( ENew % NodeIndexes, Eold % TYPE % NumberOfNodes )
        Enew % NodeIndexes = Eold % NodeIndexes

        IF(.NOT. IsBulkElement ) THEN
          CALL UpdateParentElements()
        END IF

      CASE(404)
        IF( MeshDim == 3 ) THEN
          Parent => Eold % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED(Parent)) THEN
            Parent => Eold % BoundaryInfo % Right
          END IF          
          Face => Find_Face(Mesh, Eold, Parent ) 
          IsOddCut = ANY( CutCorner(Face % ElementIndex) == Face % NodeIndexes([1,3]) )
        ELSE
          Face => Eold
          IsOddCut = .TRUE.
        END IF
                    
        DO j=1,2         
          IF( IsOddCut ) THEN
            IF(j==1) THEN
              LocalMap(1:3) = [1,2,3]
            ELSE
              LocalMap(1:3) = [3,4,1]
            END IF
          ELSE
            IF(j==1) THEN
              LocalMap(1:3) = [4,1,2]
            ELSE
              LocalMap(1:3) = [2,3,4]
            END IF
          END IF

          NewElCnt = NewElCnt + 1
          Enew => NewElements(NewElCnt)
          Enew = Eold
          Enew % TYPE => GetElementType(303)
          Enew % ElementIndex = NewElCnt         
          CALL AllocateVector( ENew % NodeIndexes, 3 )
          Enew % NodeIndexes(1:3) = Face % NodeIndexes(LocalMap(1:3))
          Enew % EdgeIndexes => NULL()
          Enew % FaceIndexes => NULL()
          
          IF(.NOT. IsBulkElement ) THEN
            CALL UpdateParentElements()
          END IF
        END DO

      CASE(706)

        LocalCut(1:3) = CutCorner(Eold % FaceIndexes(3:5))

        ! How are we cutting the three faces? 
        Is24 = ANY( LocalCut(1) == Eold % NodeIndexes([2,4])) 
        Is26 = ANY( LocalCut(2) == Eold % NodeIndexes([2,6])) 
        Is16 = ANY( LocalCut(3) == Eold % NodeIndexes([1,4])) 
        Is15 = .NOT. Is24
        Is35 = .NOT. Is26
        Is34 = .NOT. Is16

        ! Code the cases to numbers [111,222]
        CutComb = 111
        IF(.NOT. Is24) CutComb = CutComb+100  ! Is15
        IF(.NOT. Is26) CutComb = CutComb+10   ! Is35
        IF(.NOT. Is16) CutComb = CutComb+1    ! Is34
              
        DO j=1,3                   
          SELECT CASE(CutComb)
          CASE(111)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,6]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,2,3,6]
            END IF

          CASE(112)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,6]
            ELSE 
              LocalMap(1:4) = [2,3,6,4]
            END IF

          CASE(122)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,3]
            ELSE 
              LocalMap(1:4) = [4,5,3,6]
            END IF

          CASE(221)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,3,5,6]
            END IF

          CASE(222)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,3]
            ELSE 
              LocalMap(1:4) = [4,5,3,6]
            END IF

          CASE(211)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,6]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,2,3,6]
            END IF

          CASE DEFAULT
            CALL Fatal(Caller,'Unknown case for split: '//I2S(CutComb))

          END SELECT

          NewElCnt = NewElCnt + 1
          Enew => NewElements(NewElCnt)
          Enew = Eold
          Enew % TYPE => GetElementType(504)
          Enew % ElementIndex = NewElCnt         
          CALL AllocateVector( ENew % NodeIndexes, 4 )
          Enew % NodeIndexes = Eold % NodeIndexes(LocalMap(1:4))          
          Enew % EdgeIndexes => NULL()
          Enew % FaceIndexes => NULL()
        END DO
      END SELECT
            
      IF(i==Mesh % NumberOfBulkElements) THEN
        nBulkElems = NewElCnt
      END IF
    END DO

    ! Release old elements and replace them with new elements and element counts
    CALL ReleaseMeshElements( Mesh ) 

    Mesh % Elements => NewElements
    Mesh % NumberOfBulkElements = nBulkElems
    Mesh % NumberOfBoundaryElements = NewElCnt - nBulkElems
            
    ! These are now conservative and could be updated
    ! NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs

    IF( MeshDim == 3 ) THEN
      Mesh % MaxElementNodes = 4
    ELSE
      Mesh % MaxElementNodes = 3
    END IF
    
#if 0
    j = 0
    DO i=1,Mesh % NumberOfBulkElements
      Enew => NewElements(i)        
      IF ( Enew % DGDOFs>0 ) THEN
        Enew % DGDofs == Enew % TYPE % NumberOfNodes
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF
    END DO
#endif

    DO i=1,NewElCnt 
      IF (i<=Mesh % NumberOfBulkElements) THEN
        Enew => Mesh % Elements(i)
        PDefs => Enew % PDefs
        IF(ASSOCIATED(PDefs)) THEN
          CALL AllocatePDefinitions(Enew)
          Enew % PDefs = PDefs
          
          ! All elements in actual mesh are not edges
          Enew % PDefs % isEdge = .FALSE.
          
          ! If element is of type tetrahedron and is a p element,
          ! do the Ainsworth & Coyle trick
          IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
          CALL GetRefPElementNodes( Enew % TYPE,  Enew % TYPE % NodeU, &
              Enew % TYPE % NodeV, Enew % TYPE % NodeW )
        END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO
    
    CALL CheckTimer(Caller,Delete=.TRUE.)
    
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    IF( Parallel ) THEN
      ! We don't have any updates for parallel mesh.
      ! The nodes stay the same.
    END IF

    ! Release old edges and faces since they don't match the new mesh.
    !-----------------------------------------------------------------
    CALL ReleaseMeshEdgeTables( Mesh )
    CALL ReleaseMeshFaceTables( Mesh )
    !Mesh % Faces => NULL()
    !Mesh % Edges => NULL()
    !Mesh % NumberOfFaces = 0 
    !Mesh % NumberOfEdges = 0 
    
    IF( FacesPresent ) THEN 
      CALL Info(Caller,'Generating faces in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshFaces3D( Mesh )
    END IF
    IF( EdgesPresent ) THEN 
      CALL Info(Caller,'Generating faces in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshEdges( Mesh )
    END IF

    !CALL CheckMeshInfo( Mesh )     
    !CALL writemeshtodisk( Mesh, "koe" )

  CONTAINS

    
    SUBROUTINE UpdateParentElements()

      INTEGER :: j,m,lcnt,l,nCands,ElemCode,BulkOffset
      TYPE(Element_t), POINTER :: CandParent, Parent
      LOGICAL :: hit
      
      ALLOCATE( Enew % BoundaryInfo )
      Enew % BoundaryInfo = Eold % BoundaryInfo 

      DO j=1,2
        IF(j==1) THEN
          Parent => Eold % BoundaryInfo % Left
        ELSE
          Parent => Eold % BoundaryInfo % Right
        END IF
        IF(.NOT. ASSOCIATED( Parent ) ) CYCLE

        ElemCode = Parent % TYPE % ElementCode

        ! Depending on the elementtype the original element has been split into several candidate elements. 
        SELECT CASE(ElemCode)
        CASE(202)
          nCands = 1
        CASE(303)
          nCands = 1
        CASE(404)
          nCands = 2
        CASE(504)
          nCands = 1
        CASE(706)
          nCands = 3
        END SELECT

        hit = .FALSE.

        BulkOffset = BulkElementOffset(Parent % ElementIndex)

        DO k=1,nCands
          CandParent => NewElements(BulkOffset+k)
          m = Enew % Type % NumberOfNodes
          lcnt = 0
          DO l=1,m
            IF(ANY(Enew % NodeIndexes(l) == CandParent % NodeIndexes)) lcnt = lcnt + 1
          END DO
          IF(lcnt == m) THEN
            IF(j==1) THEN
              Enew % BoundaryInfo % Left => CandParent
            ELSE
              Enew % BoundaryInfo % Right => CandParent
            END IF
            hit = .TRUE.
            EXIT
          END IF
        END DO
            
        IF(.NOT. Hit) THEN
          PRINT *,'Not Found:',j,k,lcnt,nCands,ElemCode,Parent % ElementIndex, BulkOffset
          PRINT *,'This:',Eold % NodeIndexes
          PRINT *,'Parent:',Parent % NodeIndexes
          DO k=1,nCands
            CandParent => NewElements(BulkOffset+k)
            PRINT *,'Cands:',CandParent % NodeIndexes
          END DO
          CALL Fatal('UpdateParentElements','Could not find parent for type '//I2S(ElemCode))
        END IF
        
      END DO

    END SUBROUTINE UpdateParentElements
    
     
   END SUBROUTINE SplitMeshQuads

  FUNCTION SplitMeshLevelset(Mesh,Vlist) RESULT( NewMesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist    
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: phi(:)
    INTEGER, ALLOCATABLE :: EdgeSplit(:)
    LOGICAL, ALLOCATABLE :: CutNode(:)
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: SplitReady    
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:)
    REAL(KIND=dp) :: Eps
    INTEGER, POINTER :: NodeIndexes(:), EdgeIndexes(:)    
    INTEGER :: i, j, j2, j3, k, k2, k3, l, l2, l3, m, n, &
        n_old, n_new, n_cut, n_split, n_neg, n_pos
    INTEGER :: NoHits, NewElCnt, BCCnt, prevl, &
        NodeCnt, FaceCnt, Node, ParentId 
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Parent 
    INTEGER, POINTER :: Child(:,:)
    REAL(KIND=dp) :: h1,h2,hprod,r,s1,s2 
    REAL(KIND=dp), POINTER :: stime(:)
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER :: BodyOffset, SgnNode, BodyCount, LevelsetBC
    LOGICAL :: PosOffset, BulkParent, Parallel
    CHARACTER(:), ALLOCATABLE :: str       
    CHARACTER(*), PARAMETER :: Caller = 'SplitMeshLevelset'
        
!------------------------------------------------------------------------------
    CALL Info( Caller, 'Splitting finite element mesh at zero levelset!', Level = 5 )

    IF ( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn(Caller,'Original mesh not associated!')
      RETURN
    END IF
        
    CALL ResetTimer(Caller)
    
    DO i=1,Mesh % NumberOfBulkElements
      n = Mesh % Elements(i) % TYPE % ElementCode 
      IF( n /= 303 .AND. n /= 504 ) THEN
        CALL Fatal(Caller,'Only linear triangles and tets can be split: '//I2S(n))
      END IF
    END DO
    
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh )
        
    CALL Info( Caller, '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( Caller, Message, Level=6 )

    ! At this stage the coordinates have not been added as variable.
    ! We cannot use the UDF's if these are not available. Also time
    ! is needed by default in some calls. 
    Var => VariableGet( Mesh % Variables,'time')
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )   
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )    
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )    
      ALLOCATE(stime(1)); stime(1) = 0.0_dp
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Time',DOFs=1, Values=sTime )
      CurrentModel % Variables => Mesh % Variables
    END IF
    
    ! Initialize the levelset function for all nodes
    n_old = Mesh % NumberOfNodes
    ALLOCATE( Phi(n_old) )
    
    str = ListGetString( Vlist,'Levelset Variable', Found)
    IF( Found ) THEN
      Var => VariableGet(Mesh % Variables, str)
      IF(.NOT. ASSOCIATED(Var) ) THEN
        CALL Fatal(Caller,'"Levelset Variable" requested, but not available: '//TRIM(str))
      END IF
      Phi = 1.0_dp
      ! We revert to nodal indexes since it will be easier in the future!
      DO i=1,n_old
        j = Var % Perm(i)
        IF(j>0) Phi(i) = Var % Values(j)
      END DO
    ELSE      
      DO i=1,n_old
        Phi(i) = ListGetRealAtNode(Vlist,'Levelset Function', i, Found)
        IF(.NOT. Found ) THEN
          CALL Fatal(Caller,'"Levelset Function" needed to enrich the mesh!')             
        END IF
      END DO
    END IF
    
    Eps = ListGetCReal( Vlist,'Levelset Epsilon',Found )
    IF(.NOT. Found ) Eps = 1.0e-3
    
    n_pos = COUNT( Phi > 0.0 )
    n_neg = COUNT( Phi < 0.0 ) 
        
    BodyOffset = ListGetInteger( Vlist,'Levelset Body Offset',Found ) 
    PosOffset = ListGetLogical( Vlist,'Levelset Offset Positive',Found ) 
    LevelsetBC = ListGetInteger( Vlist,'Levelset Boundary',Found )
    IF(.NOT. Found) LevelsetBC = CurrentModel % NumberOfBCs
    
    IF( Parallel ) THEN
      n_pos = ParallelReduction(n_pos) 
      n_neg = ParallelReduction(n_neg)
    END IF
    
    CALL Info(Caller,'Positive and negative values: '&
        //I2S(n_pos)//' vs. '//I2S(n_neg),Level=7)    
    
    IF( n_pos == 0 .OR. n_neg == 0 ) THEN
      CALL Warn(Caller,'Nothing to do, no zero levelset available!')
      RETURN
    END IF
       
    ! We need edges in order to do the splitting!
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT. EdgesPresent) CALL FindMeshEdges( Mesh )
        
    ALLOCATE( EdgeSplit(Mesh % NumberOfEdges), CutNode(n_old) )    
    EdgeSplit = 0
    CutNode = .FALSE.
        
    j = 0
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      h1 = Phi(NodeIndexes(1))
      h2 = Phi(NodeIndexes(2))
      hprod = h1*h2
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        IF( r <= Eps ) THEN
          CutNode(NodeIndexes(2)) = .TRUE.
        ELSE IF(1.0-r < Eps ) THEN
          CutNode(NodeIndexes(1)) = .TRUE.
        ELSE
          j = j+1 
          EdgeSplit(i) = j
        END IF
      ELSE IF( ABS(hprod) < 1.0d-20 ) THEN
        IF(ABS(h1) < 1.0e-20) CutNode(NodeIndexes(1)) = .TRUE. 
        IF(ABS(h2) < 1.0e-20) CutNode(NodeIndexes(2)) = .TRUE.
      END IF
    END DO
    
    n_new = j
    CALL Info(Caller,'Number of additional nodes: '//I2S(n_new),Level=6)

    j = COUNT( CutNode )
    CALL Info(Caller,'Number of cut nodes: '//I2S(j),Level=6)
    
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = n_old + n_new 

!   Create the new mesh
!   -------------------------------
    NewMesh => AllocateMesh()    
    NewMesh % SingleMesh = Mesh % SingleMesh
    NewMesh % Name = Mesh % Name   

    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:n_old) = u
    y(1:n_old) = v
    z(1:n_old) = w

!   add new nodes where edges are split:
!   ------------------------------------
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      j = EdgeSplit(i)
      IF( j > 0 ) THEN
        j = j + n_old
        h1 = Phi(NodeIndexes(1))
        h2 = Phi(NodeIndexes(2))
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        x(j) = r*u(NodeIndexes(1)) + (1-r)*u(NodeIndexes(2))
        y(j) = r*v(NodeIndexes(1)) + (1-r)*v(NodeIndexes(2))
        z(j) = r*w(NodeIndexes(1)) + (1-r)*w(NodeIndexes(2))
      END IF
    END DO

    CALL Info(Caller,'Added new nodes on the splitted edges.', Level=10 )  

    
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt

!   Update bulk elements:
!   =====================
!
!   First count maximum number of new elements:
!   -------------------------------------------
    NewElCnt = 0
    BodyCount = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Eold => Mesh % Elements(i)
      j = 1

      Found = .FALSE.
      IF( ASSOCIATED( Eold % EdgeIndexes ) ) THEN
        Found = ANY(EdgeSplit(Eold % EdgeIndexes) > 0 )
      ELSE
        CALL Fatal(Caller,'No edges for element: '//I2S(i))
      END IF
      
      IF( Found ) THEN
        SELECT CASE( Eold % TYPE % ElementCode/100 )                
        CASE(2)
          j = 2
        CASE(3)
          j = 3
        CASE(5)
          j = 6
        END SELECT
        ! There will also be additional BC elements on the cut!
        j = j + 1
      END IF
      NewElCnt = NewElCnt + j
    END DO
    
    CALL Info( Caller,'Maximum estimated count of new elements: '//I2S(NewElCnt), Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info(Caller,'New mesh allocated.', Level=20 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 6 )
    Child = 0
    CALL Info(Caller,'Array for bulk elements allocated.', Level=20 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes

!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)
       NodeIndexes => Eold % NodeIndexes       
       n = Eold % TYPE % NumberOfNodes                
       n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )

       ! We continue splitting until the element is exhausted
       SplitReady = .FALSE.

       ! Split elements to no more than 6 pieces
       DO m = 1,6
         NewElCnt = NewElCnt + 1
         Child(i,m) = NewElCnt
         Enew => NewMesh % Elements(NewElCnt)

         Enew = Eold
         Enew % TYPE => Eold % TYPE
         Enew % BodyId = Eold % BodyId
         Enew % PartIndex = Eold % PartIndex
         Enew % ElementIndex = NewElCnt
         Enew % NDOFs = Eold % NDOFs
         Enew % EdgeIndexes => NULL()
         Enew % FaceIndexes => NULL()
         Enew % BoundaryInfo => NULL()
         
         CALL AllocateVector( ENew % NodeIndexes, n)
         
         IF( n_split == 0 ) THEN
           Enew % NodeIndexes = NodeIndexes
           DO j=1,n
             IF(.NOT. CutNode(NodeIndexes(j)) ) THEN
               ! This is a representative node that is used to determine the sign of the
               ! new elements in order to decide whether to add offset for body or not. 
               SgnNode = j
               EXIT
             END IF
           END DO
           
           SplitReady = .TRUE.
         ELSE           
           n_cut = COUNT( CutNode(NodeIndexes) )
       
           IF ( Eold % TYPE % ElementCode == 303 ) THEN         
             ! Split triangle to four triangles split on one or two edges
             !-----------------------------------------------------------
             IF( n_split == 2 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               IF( m == 1 ) THEN
                 ! There are two ways to split the triangle.
                 ! Choose the one with shorter diameter.
                 s1 = (x(NodeIndexes(j)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (y(NodeIndexes(j)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (z(NodeIndexes(j)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2
                 s2 = (x(NodeIndexes(j2)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (y(NodeIndexes(j2)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (z(NodeIndexes(j2)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2
                 Enew % NodeIndexes(1) = NodeIndexes(j)
                 Enew % NodeIndexes(2) = NodeIndexes(j2)                 
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 ELSE
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
                 END IF
                 SgnNode = j
               ELSE IF(m==2) THEN
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(1) = NodeIndexes(j)
                   SgnNode = j
                 ELSE
                   Enew % NodeIndexes(1) = NodeIndexes(j2)                   
                   SgnNode = j2
                 END IF
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                
               ELSE IF(m==3) THEN
                 Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = NodeIndexes(j3)
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF

             ELSE IF( n_split == 1 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               ! One cut result to splitted elements only if the opposing node is cut through
               IF( .TRUE. .OR. CutNode(NodeIndexes(j3)) ) THEN
                 IF(m==1) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j2)
                   Enew % NodeIndexes(3) = NodeIndexes(j3)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j2
                   ELSE
                     SgnNode = j3
                   END IF
                 ELSE IF(m==2) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j3)
                   Enew % NodeIndexes(3) = NodeIndexes(j)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j
                   ELSE
                     SgnNode = j3
                   END IF 
                   SplitReady = .TRUE.
                 END IF
               ELSE
                 Enew % NodeIndexes = NodeIndexes
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF
             ELSE
               CALL Fatal(Caller,'Triangle can only deal with 1 and 2 splits!')
             END IF
           ELSE              
             CALL Fatal(Caller,'Element type '//I2S(Eold % TYPE % ElementCode)//&
                 ' not supported by the levelset splitter.')
           END IF
         END IF

         ! Set offset for inside/outside elements of the zero levelset.
         ! The SgnNode is a representative node the sign of which tells whether we are inside
         ! or outside. 
         IF( PosOffset ) THEN
           IF( Phi(NodeIndexes(SgnNode)) > 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset
             BodyCount = BodyCount + 1
           END IF
         ELSE
           IF( Phi(NodeIndexes(SgnNode)) < 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset            
             BodyCount = BodyCount + 1
           END IF
         END IF
         IF( SplitReady ) EXIT
       END DO
     END DO
     
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt
    
    CALL Info(Caller,'Number of elements inside: '//I2S(BodyCount),Level=7)
    
   
!   Update boundary elements:
!   ---------------------------------------------------

    BCCnt = 0
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IF( i == Mesh % NumberOfBulkElements + 1 ) THEN
        CALL Info(Caller,'Number of boundary elements from bulk cuts: '//I2S(BCCnt))           
        BCCnt = 0
      END IF
     
      Eold => Mesh % Elements(i)
      NodeIndexes => Eold % NodeIndexes             
      BulkParent = ( i <= Mesh % NumberOfBulkElements )
      n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )
      n_cut = COUNT( CutNode(NodeIndexes) )

      ! Elements created from bulk cuts require some splits or cuts.
      ! Existing boundary elements remain even without cuts.
      IF( BulkParent ) THEN
        IF( n_split + n_cut <= 1 ) CYCLE
      END IF
  
      SplitReady = .FALSE.
      
      ! Each existing boundary element may be cut to several pieces
      ! For triangles this is just max two!
      DO m=1,10          
        BCCnt = BCCnt + 1
        NewElCnt = NewElCnt + 1
        IF( NewElCnt > SIZE( NewMesh % Elements ) ) THEN
          CALL Fatal(Caller,'Too few elements allocated: '//I2S(NewElCnt))
        END IF
       
        Enew => NewMesh % Elements(NewElCnt)
        
        ALLOCATE(Enew % BoundaryInfo)         
        Enew % PartIndex = Eold % PartIndex
        Enew % ElementIndex = NewElCnt
        
        n = 2
        Enew % TYPE => GetElementType(202)
        CALL AllocateVector( ENew % NodeIndexes, n)
        Enew % NDOFs = n
        Enew % EdgeIndexes => NULL()
        Enew % FaceIndexes => NULL()
                
        IF( BulkParent ) THEN
          ! There are the new boundary elements that come from splitting the mesh
          ! at zero levelset. Give the boundary a new index. 
          Enew % BoundaryInfo % Constraint = LevelsetBC
          
          IF ( Eold % TYPE % ElementCode == 303 ) THEN         
            IF( n_split == 2 ) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
            ELSE IF( n_split == 1 .AND. n_cut == 1) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
              END DO
              !j2 = MODULO(j,3)+1
              !j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
              DO j2=1,3
                IF( CutNode(NodeIndexes(j2)) ) EXIT
              END DO
              Enew % NodeIndexes(2) = NodeIndexes(j2)
            ELSE IF( n_cut == 2) THEN
              DO j=1,3
                IF( .NOT. CutNode(NodeIndexes(j) ) ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = NodeIndexes(j2)
              Enew % NodeIndexes(2) = NodeIndexes(j3)
            ELSE
              CALL Fatal(Caller,'Can only deal with 2 or 1+1 splits!')
            END IF
          ELSE              
            CALL Fatal(Caller,'Element type '//I2S(Eold % TYPE % ElementCode)//&
                ' not supported by the levelset splitting.')
          END IF          
          SplitReady = .TRUE.
          
        ELSE
          ! Each existing boundary element may be cut to several pieces
          Enew % BoundaryInfo = Eold % BoundaryInfo
          
          IF( n_split == 0 ) THEN
            ! If no edge is split the element stays as is
            Enew % NodeIndexes = Eold % NodeIndexes
            SplitReady = .TRUE.
            
          ELSE IF( Eold % TYPE % ElementCode == 202 ) THEN
            IF(m==1) THEN
              Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
            ELSE IF(m==2) THEN
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
              Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
              SplitReady = .TRUE.
            END IF
          ELSE
            CALL Fatal(Caller,'Cannot do this element yet!')
          END IF
        END IF

         
        prevl = 0
        DO k=1,2
          ! Pointer to the found left/right bulk element
          Eptr => NULL()

          IF( BulkParent ) THEN
            ! If the boundary results from splitting existing elements then
            ! the parent is the existing bulk elements. 
            Parent => Mesh % Elements(i)
          ELSE
            ! If boundary results from existing boundary elements then the potential
            ! parents are the children of the old parents. 
            IF( k==1 ) THEN
              Parent => Eold % BoundaryInfo % Left
            ELSE            
              Parent => Eold % BoundaryInfo % Right
            END IF
            IF(.NOT. ASSOCIATED(Parent)) CYCLE
          END IF

          ! Find the correct parent among the splitted children of the
          ! initial bulk elements. There may be 1 or several children. 
          DO k2 = 1, 6            
            l = Child( Parent % ElementIndex, k2 )
            IF(l==0) CYCLE
            NoHits = 0
            
            IF( BulkParent ) THEN
              IF( k==2 .AND. l == prevl ) CYCLE
            END IF

            IF(l == 0 .OR. l > SIZE( NewMesh % Elements) ) THEN
            ! This is left for debugging...
#if 1
              PRINT *,'Size Elements:',l, SIZE(NewMesh % Elements)
              PRINT *,'Child:',m,n_split,n_cut,k2,BulkParent
              PRINT *,'Parent index:',i,Parent % ElementIndex
              PRINT *,'Parents children',Child( Parent % ElementIndex, :)
              PRINT *,'Parent indexes:',Parent % NodeIndexes
              PRINT *,'Enew:',Enew % NodeIndexes
              PRINT *,'Eold:',Eold % NodeIndexes
              PRINT *,'Eold edges:',Eold % EdgeIndexes
              PRINT *,'Eold edge node indexes:',Mesh % Edges(Eold % EdgeIndexes(1) ) % NodeIndexes
              DO l2=1,6
                IF(Child( Parent % ElementIndex, l2) == 0 ) EXIT
                PRINT *,'Parent:',l2,NewMesh % &
                    Elements(Child( Parent % ElementIndex,l2)) %  NodeIndexes
              END DO
              PRINT *,'old node indexes:',Mesh % Elements(Parent % ElementIndex) %  NodeIndexes
              PRINT *,'old edge indexes:',Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes
              PRINT *,'old cut indexes:',CutNode(Mesh % Elements(Parent % ElementIndex) %  NodeIndexes)
              PRINT *,'old split indexes:',EdgeSplit(Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes)
#endif
              EXIT
            END IF

            Eptr => NewMesh % Elements(l)

            DO l2 = 1,Enew % Type % NumberOfNodes 
              DO l3 = 1, Eptr % TYPE % NumberOfNodes
                IF( Enew % NodeIndexes(l2) == Eptr % NodeIndexes(l3) ) THEN
                  NoHits = NoHits + 1
                  EXIT
                END IF
              END DO
            END DO
            
            IF( NoHits == n ) EXIT
          END DO
          
          IF( NoHits == n ) THEN
            IF( k==1) THEN
              prevl = l
              Enew % BoundaryInfo % Left => Eptr
            ELSE
              Enew % BoundaryInfo % Right => Eptr
            END IF
          ELSE
            IF(k==1) CALL Warn(Caller,'Could not find even 1 parent!')
          END IF
            
        END DO
       
       ! When we have created all the new boundary elements resulting from splitting
       ! the master element then proceed to next element. 
       IF(SplitReady) EXIT
     END DO
   END DO


!   Update new mesh element counts:
!   -------------------------------
   CALL Info(Caller,'Number of total elements: '//I2S(NewElCnt),Level=7)
    
!   Update new mesh boundary element counts:
!   ----------------------------------------
   NewMesh % NumberOfBoundaryElements = NewElCnt - &
       NewMesh % NumberOfBulkElements
   NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
   NewMesh % MaxElementNodes = Mesh % MaxElementNodes
   
    
   CALL Info( Caller, '******** New mesh ********', Level=6 )
   WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
   CALL Info( Caller, Message, Level=6 )


   ! Information of the new system size, also in parallel
   !----------------------------------------------------------------------
   ParTmp(1) = Mesh % NumberOfNodes
   ParTmp(2) = Mesh % NumberOfBulkElements
   ParTmp(3) = Mesh % NumberOfBoundaryElements
   ParTmp(4) = NewMesh % NumberOfNodes
   ParTmp(5) = NewMesh % NumberOfBulkElements
   ParTmp(6) = NewMesh % NumberOfBoundaryElements
   
   IF( .FALSE. .AND. Parallel ) THEN
     CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
     
     CALL Info(Caller,'Information on parallel mesh sizes',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(1))//' nodes',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(2))//' bulk elements',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(3))//' boundary elements',Level=8)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(4))//' nodes',Level=5)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(5))//' bulk elements',Level=5)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(6))//' boundary elements',Level=5)
   END IF

    
   ! Update structures needed for parallel execution:
   !--------------------------------------------------
   IF( Parallel ) THEN
     CALL UpdateParallelInfo( Mesh, NewMesh )
   END IF

   ! Finalize:
   !-----------
   IF(.NOT.EdgesPresent) THEN
     CALL ReleaseMeshEdgeTables( Mesh )
     CALL ReleaseMeshFaceTables( Mesh )
   ELSE
     CALL FindMeshEdges( NewMesh )
   END IF

   CALL CheckTimer(Caller,Delete=.TRUE.)

   CALL Info(Caller,'Mesh was enriched with zero levelset',Level=8)

 CONTAINS
    
!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelInfo( Mesh, NewMesh )
!------------------------------------------------------------------------------
      TYPE(Mesh_t), TARGET :: Mesh, NewMesh
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Edge
      INTEGER :: i,j1,j2,n,n0,m,istat
      LOGICAL :: Found
!------------------------------------------------------------------------------
!
!      Update mesh interfaces for parallel execution.
!      ==============================================
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) CALL Fatal( Caller, 'Allocate error.' )
       DO i=1,n
         NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       CALL AllocateVector( NewMesh % ParallelInfo % GInterface,n  )       
       NewMesh % ParallelInfo % GInterface = .FALSE.

       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )
       NewMesh % ParallelInfo % GlobalDOFs = 0

       ! Inherit the old parallel data
       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % GInterface(1:n) = Mesh % ParallelInfo % GInterface
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = Mesh % ParallelInfo % GlobalDOFs
       DO i=1,n
         m = SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) 
         ALLOCATE( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(m) )
         NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO

       n0 = ParallelReduction(MAXVAL(Mesh % ParallelInfo % GlobalDofs),2)       
       CALL Info(Caller,'Offset for parallel numbering of new nodes: '//I2S(n0))

       ! We need global numbering for the edges that we use for the unique numbering of new nodes
       CALL SParEdgeNumbering(Mesh)
       
       DO i=1,Mesh % NumberOfEdges
         j = EdgeSplit(i)
         IF(j==0) CYCLE
         Edge => Mesh % Edges(j)

         ! Make a unique parallel number for the new nodes introduced at split edges
         NewMesh % ParallelInfo % GlobalDOFs(n+j) = n0 + Edge % GElementIndex         

         j1 = Edge % NodeIndexes(1)
         j2 = Edge % NodeIndexes(2)
         m = CountSameIntegers(Mesh % ParallelInfo % NeighbourList(j1) % Neighbours, &
             Mesh % ParallelInfo % NeighbourList(j2) % Neighbours, &
             NewMesh % ParallelInfo % NeighbourList(n+j) % Neighbours ) 
         NewMesh % ParallelInfo % GInterface(n+j) = (m>1)
       END DO
       
    END SUBROUTINE UpdateParallelInfo
    
  END FUNCTION SplitMeshLevelset


END MODULE MeshSplit
