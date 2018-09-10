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

!Remesh the calving model using MMG3D
SUBROUTINE CalvingRemeshMMG( Model, Solver, dt, Transient )

  USE MeshUtils
  USE CalvingGeometry
  USE MMG3DTools

  IMPLICIT NONE

#include "mmg/mmg3d/libmmgtypesf.h"

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !--------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh,NewMesh,NewMeshR
  TYPE(Element_t),POINTER :: Element, ParentElem
  TYPE(Element_t), ALLOCATABLE, TARGET :: PElements(:)
  TYPE(Nodes_t), TARGET :: Remesh_nodes
  INTEGER :: i,j, NNodes,NBulk, Nbdry, ierr, my_cboss,MyPE, PEs,CCount, counter, &
       my_calv_front,calv_front, ncalv_parts, group_calve, comm_calve, group_world,ecode, NElNodes,&
       NBulk_send, NBdry_send, NBulk_fixed, NBdry_fixed,test_min, GlNode_max,NNode_remesh,NEl_remesh,&
       adjList(4)
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, ALLOCATABLE :: Prnode_count(:), cgroup_membs(:),&
       PNVerts(:),PNTris(:), PNTetras(:),FixedElems(:), GDOFs_send(:),PGDOFs_send(:),&
       el_stream(:), pel_stream(:), pnbulk_send(:), pnbdry_send(:), pnbulk_fixed(:), &
       pnbdry_fixed(:),pcalv_front(:),PElem_parts(:),GtoLNN(:), LtoGNN(:),&
       disps(:),el_stream_sizes(:),elem_fixed(:),Pelem_fixed(:),pnode_parts(:),Nodeno_map(:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), test_lset(:), NodeCoords_send(:),PNodeCoords_send(:),&
       PLset_send(:)
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), fixed_node(:), &
       elem_send(:), bulk_fixed(:), bdry_fixed(:),bulk_send(:), bdry_send(:), RmElem(:), RmNode(:)
  LOGICAL :: ImBoss
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  SolverParams => GetSolverParams()
  SolverName = "CalvingRemeshMMG"

  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes
  NBulk = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements

  calv_front = 1
  MyPe = ParEnv % MyPE
  PEs = ParEnv % PEs

  !Check all elements are tetras or tris:
  DO i=1,NBulk+Nbdry
    Element => Mesh % Elements(i)
    ecode = Element % TYPE % ElementCode
    IF(ecode /= 101 .AND. ecode /= 202 .AND. ecode /= 303 .AND. ecode /= 504) THEN
      PRINT *,MyPE,' has unsupported element type: ',ecode
      CALL Fatal(SolverName, "Invalid Element Type!")
    END IF
  END DO

  !For testing - set a calving levelset function to mimick calving events
  !-------------------

  test_point = (/491600.0, -2290000.0,79.9166641235/)
  test_thresh = 1500.0
  remesh_thresh = 2000.0

  ALLOCATE(test_dist(NNodes),&
       test_lset(NNodes),&
       remeshed_node(NNodes),&
       calved_node(NNodes),&
       fixed_node(NNodes),&
       )

  test_dist = 0.0
  remeshed_node = .FALSE.
  calved_node = .FALSE.
  fixed_node = .FALSE.

  !Compute involved nodes
  DO i=1, NNodes
    test_dist(i) = (((Mesh % Nodes % x(i) - test_point(1))**2.0) + &
         ((Mesh % Nodes % y(i) - test_point(2))**2.0) + &
         ((Mesh % Nodes % z(i) - test_point(3))**2.0)) ** 0.5
  END DO

  calved_node = test_dist < test_thresh
  remeshed_node = test_dist < (test_thresh + remesh_thresh)
  test_lset = test_dist - test_thresh

  my_calv_front = 0
  IF(ANY(remeshed_node)) THEN
    my_calv_front = 1 !this is integer (not logical) so we can later move to multiple calving fronts
  END IF

  !TODO - could make this more efficent by cutting out some of the elements from the calved region.
  !This region will be remeshed too, but we discard it, so the closer we can get to the edge of the
  !calving event, the better.

  !Identify all elements which need to be sent (including those fixed in place)
  !Here we also find extra nodes which are just beyond the remeshing threshold
  !but which are in fixed elements, thus also need to be sent
  IF(my_calv_front > 0) THEN
    ALLOCATE(elem_fixed(nbulk+nbdry), elem_send(nbulk+nbdry))
    elem_fixed = 0 !integer rather than logical to facilitate sending in packed element structure
    elem_send = .FALSE.

    !Each partition identifies (based on nodes), elements which need to be transferred
    DO i=1,Nbulk+Nbdry
      Element => Mesh % Elements(i)
      NodeIndexes => Element % NodeIndexes
      NElNodes = Element % TYPE % NumberOfNodes
      IF(ANY(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
        elem_send(i) = .TRUE.
        IF(.NOT. ALL(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
          elem_fixed(i) = 1
          fixed_node(NodeIndexes(1:NElNodes)) = .TRUE.
        END IF
      END IF
    END DO
    
    remeshed_node = remeshed_node .OR. fixed_node

    DEALLOCATE(fixed_node) !reuse later
  END IF


  !This does nothing yet but it will be important - determine
  !the discrete calving zones, each of which will be seperately remeshed
  !by a nominated boss partition
  !  CALL CountCalvingEvents(Model, Mesh, ccount)
  ccount = 1

  !Negotiate local calving boss - just one front for now
  !-----------------------
  ALLOCATE(pcalv_front(PEs))
  pcalv_front = 0
  CALL MPI_AllGather(my_calv_front, 1, MPI_INTEGER, pcalv_front, &
       1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)

  PRINT *,mype,' debug calv_front: ',my_calv_front
  PRINT *,mype,' debug calv_fronts: ',pcalv_front

  !Create an MPI_COMM for each calving region, allow gathering instead of 
  !explicit send/receive
  ncalv_parts = COUNT(pcalv_front==calv_front) !only one front for now...
  ALLOCATE(cgroup_membs(ncalv_parts))
  counter = 0
  DO i=1,PEs
    IF(pcalv_front(i) == calv_front) THEN !one calve front
      counter = counter + 1
      cgroup_membs(counter) = i-1
    END IF
  END DO
  PRINT *,mype,' debug group: ',cgroup_membs

  PRINT *,'making group'
  CALL MPI_Comm_Group( ELMER_COMM_WORLD, group_world, ierr)
  CALL MPI_Group_Incl( group_world, ncalv_parts, cgroup_membs, group_calve, ierr)
  CALL MPI_Comm_create( ELMER_COMM_WORLD, group_calve, COMM_CALVE, ierr)
  PRINT *,' made group'

  !Partition is involved in calving remeshing
  IF(My_Calv_Front>0) THEN

    my_cboss = MINLOC(pcalv_front, 1, MASK=pcalv_front==calv_front)-1
    PRINT *,MyPe,' debug calving boss: ',my_cboss
    ImBoss = MyPE == my_cboss



    IF(ImBoss) THEN
      ALLOCATE(Prnode_count(ncalv_parts))
      Prnode_count = 0
    END IF

    !Gather info from partitions which are involved
    !-----------------------

    !root is always zero because cboss is always lowest member of new group (0)
    CALL MPI_Gather(COUNT(remeshed_node), 1, MPI_INTEGER, Prnode_count, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)

    IF(ImBoss) THEN
      ALLOCATE(PNverts(ncalv_parts),&
           PNtetras(ncalv_parts),&
           PNtris(ncalv_parts),&
           PGDOFs_send(SUM(Prnode_count)+(2*ncalv_parts)), &
           Plset_send(SUM(Prnode_count)), &
           PNodeCoords_send(SUM(Prnode_count)*3), &
           disps(ncalv_parts)&
      )
    ELSE
      ALLOCATE(disps(1), prnode_count(1))
    END IF

    CALL PackNodesToSend(Mesh, (remeshed_node), GDOFs_send, NodeCoords_send, 3)
    PRINT *,ParEnv % MyPE, ' debug count send: ',COUNT(remeshed_node),' size: ',SIZE(GDOFs_send)

    CALL MPI_BARRIER(COMM_CALVE)
    IF(ImBoss) THEN
      disps(1) = 0
      DO i=2,ncalv_parts
        disps(i) = disps(i-1) + prnode_count(i-1) + 2
      END DO
    END IF

    !Gather node indices
    CALL MPI_GatherV(GDOFs_send, COUNT(remeshed_node)+2, MPI_INTEGER, PGDOFs_send, Prnode_count+2, &
         disps, MPI_INTEGER, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    IF(ImBoss) THEN
      disps(1) = 0
      DO i=2,ncalv_parts
        disps(i) = disps(i-1) + prnode_count(i-1)
      END DO
    END IF

    !Gather node positions
    CALL MPI_GatherV(NodeCoords_send, COUNT(remeshed_node)*3, MPI_DOUBLE_PRECISION, PNodeCoords_send,&
         Prnode_count*3, disps*3, MPI_DOUBLE_PRECISION, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    !Gather the level set function
    CALL MPI_GatherV(PACK(test_lset,remeshed_node), COUNT(remeshed_node), MPI_DOUBLE_PRECISION, &
          Plset_send, Prnode_count, disps, MPI_DOUBLE_PRECISION, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    !Note - we don't need to send the fixed_node info because cboss can determine this from the
    !set of fixed elements.

    IF(ImBoss) CALL UnpackNodesSent(PGDOFs_send,PNodeCoords_send, Remesh_Nodes, 3, pnode_parts)
    NNode_remesh = SIZE(PGDOFs_send)

    IF(ImBoss) PRINT *,'debug got nodes: ',SIZE(Remesh_Nodes % x),SIZE(PGDOfs_send),MINVAL(PGDOfs_send),MAXVAL(PGDOfs_send)

    CALL PackElemsToSend(Mesh, elem_send, el_stream, elem_fixed)

    IF(ImBoss) ALLOCATE(el_stream_sizes(ncalv_parts))

    CALL MPI_Gather(SIZE(el_stream), 1, MPI_INTEGER, el_stream_sizes, 1,&
         MPI_INTEGER, 0,  COMM_CALVE, ierr)

    IF(ImBoss) THEN
      disps = 0
      DO i=2,ncalv_parts
        disps(i) = disps(i-1) + el_stream_sizes(i-1)
      END DO
      ALLOCATE(PEl_stream(SUM(el_stream_sizes)))
    END IF

    CALL MPI_GatherV(el_stream, SIZE(el_stream), MPI_INTEGER, PEl_stream,&
         el_stream_sizes, disps, MPI_INTEGER, 0, COMM_CALVE,ierr)

    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")


    !Now boss does the remeshing etc
    !--------------------------------
    IF(ImBoss) THEN
      CALL UnpackElemsSent(PEl_stream, PElements,  3, PElem_parts, pelem_fixed)

      NEl_remesh = SIZE(PElements)
      PRINT *,' Debug size PElements: ',NEl_remesh

      GlNode_max = 0
      DO i=1,NEl_remesh
        IF(.NOT. ASSOCIATED(PElements(i) % NodeIndexes)) &
             CALL Fatal(SolverName,"Programm error: unassociated elements")
        GlNode_max = MAX(MAXVAL(PElements(i) % NodeIndexes),GlNode_max)
      END DO
      PRINT *,' debug max elem nodeindices: ',GlNode_max
      !-------------------------

      !Need to construct local node numbers and perms, and fixed_node
      ALLOCATE(LtoGNN(NNode_remesh),GtoLNN(GLNode_max),fixed_node(NNode_remesh))
      DO i=1,NNode_remesh
        LtoGNN(i) = PGDofs_send(i)
        GtoLNN(PGDofs_send(i)) = i
      END DO

      DO i=1,NEl_remesh
        IF(PElem_fixed(i) == 1) fixed_node(GtoLNN(PElements(i) % NodeIndexes)) = .TRUE.
      END DO
      PRINT *,' debug fixed count, size: ',COUNT(fixed_node),SIZE(fixed_node), COUNT(pelem_fixed==1),SiZE(pelem_fixed)

      !Point temporary mesh to elements and nodes
      NewMesh => AllocateMesh()
      NewMesh % Elements => PElements
      NewMesh % Nodes => Remesh_Nodes

      NewMesh % NumberOfNodes = NNode_remesh

      NewMesh % NumberOfBulkElements = 0
      NewMesh % NumberOfBoundaryElements = 0
      DO i=1,NEl_remesh
        NElNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
        NewMesh % Elements(i) % NodeIndexes(1:NElNodes) = &
             GtoLNN(NewMesh % Elements(i) % NodeIndexes(1:NElNodes))
        IF(NewMesh % Elements(i) % TYPE % DIMENSION == 3) THEN
          NewMesh % NumberOfBulkElements = NewMesh % NumberOfBulkElements + 1
          ! PRINT *,'input elem ',i,' body id: ',NewMesh % Elements(i) % BodyID
        ELSE
          NewMesh % NumberOfBoundaryElements = NewMesh % NumberOfBoundaryElements + 1
          ! PRINT *,'input elem ',i,' bdry id: ',NewMesh % Elements(i) % BoundaryInfo % Constraint
          ! PRINT *,'input elem ',i,' bdry body id: ',NewMesh % Elements(i) %  BodyID
        END IF
      END DO

      !Initialise MMG datastructures
      mmgMesh = 0
      mmgSol  = 0

      CALL MMG3D_Init_mesh(MMG5_ARG_start, &
           MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
           MMG5_ARG_end)

      CALL Set_MMG3D_Mesh(NewMesh)

      !Request isosurface discretization
      CALL MMG3D_Set_iparameter(mmgMesh, mmgSol, MMGPARAM_iso, 1,ierr)

      !Set geometric parameters for remeshing
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmin,&
          25.0_dp,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmax,&
          400.0_dp,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hausd,&
           5.0_dp,ierr)

      !Feed in the level set distance
      CALL MMG3D_SET_SOLSIZE(mmgMesh, mmgSol, MMG5_Vertex, NNode_remesh ,MMG5_Scalar, ierr)
      DO i=1,NNode_remesh
        CALL MMG3D_Set_scalarSol(mmgSol,&
             Plset_send(i), &
             i,ier)
      END DO

      !Set required nodes and elements
      DO i=1,NNode_remesh
        IF(fixed_node(i)) CALL MMG3D_SET_REQUIREDVERTEX(mmgMesh,i,ierr)
        ! PRINT *,'Fixing node: ',i,NewMesh % Nodes % x(i), NewMesh % Nodes % y(i)
      END DO

      DO i=1,NEl_remesh
        IF(PElem_fixed(i)==1) THEN
          ! PRINT *,'Fixing element ',i,' type: ',NewMesh % Elements(i) % TYPE % ElementCode
          IF(NewMesh % Elements(i) % TYPE % DIMENSION == 3) THEN
            CALL MMG3D_SET_REQUIREDTETRAHEDRON(mmgMesh,i,ierr)
          ELSEIF(NewMesh % Elements(i) % TYPE % DIMENSION == 2) THEN
            CALL MMG3D_SET_REQUIREDTRIANGLE(mmgMesh,i,ierr)
          END IF
        END IF
      END DO

      !> 4) (not mandatory): check if the number of given entities match with mesh size
      CALL MMG3D_Chk_meshData(mmgMesh,mmgSol,ierr)
      IF ( ier /= 1 ) CALL EXIT(105)

      !> ------------------------------ STEP  II --------------------------
      !! remesh function
      CALL MMG3D_mmg3dls(mmgMesh,mmgSol,ierr)

      ! !New remeshing procedure--------------
      ! DO i=1,NNode_remesh
      !   CALL MMG3D_Set_scalarSol(mmgSol,&
      !        20.0_dp, &
      !        i,ierr)
      ! END DO
      ! CALL MMG3D_Set_iparameter(mmgMesh, mmgSol, MMG3D_IPARAM_iso, 0,ierr)
      ! CALL MMG3D_mmg3dlib(mmgMesh, mmgSol, ierr)
      ! !----------------------------

      CALL MMG3D_SaveMesh(mmgMesh,"test_out.mesh",LEN(TRIM("test_out.mesh")),ierr)

      CALL GET_MMG3D_MESH(NewMeshR)

      NNodes = NewMeshR % NumberOfNodes
      NBulk = NewMeshR % NumberOfBulkElements
      NBdry = NewMeshR % NumberOfBoundaryElements
      PRINT *, 'NewMeshR nonodes: ',NNodes, NBulk, NBdry

      !Clear out unneeded elements
      !Body elems with BodyID 3 (2) are the calving event (remaining domain)
      !Boundary elems seem to have tags 0,1,10...
      ! 0 seems to be the edge triangles, 1 the outer surface (all bcs) and 2 the new level set surf
      ALLOCATE(RmElem(NBulk+NBdry), RmNode(NNodes))
      RmElem = .FALSE.
      RmNode = .TRUE.

      DO i=1,NBulk
        Element => NewMeshR % Elements(i)
        NElNodes = Element % TYPE % NumberOfNodes
        ! PRINT *,'new element ',i,' body id: ',Element % BodyID
        IF(Element % TYPE % ElementCode /= 504) &
             PRINT *,'Programming error, bulk type: ',Element % TYPE % ElementCode
        IF(NElNodes /= 4) PRINT *,'Programming error, bulk nonodes: ',NElNodes

        !outbody - mark for deletion and move on
        IF(Element % BodyID == 3) THEN
          RmElem(i) = .TRUE.
          CYCLE
        END IF

        !Mark any nodes found in valid bulk elements 
        RmNode(Element % NodeIndexes(1:NElNodes)) = .FALSE.
      END DO

      !Same thru the boundary elements
      !If their parent elem is body 3, delete, *unless* they have constraint 10 (the calving face)
      !in which case use mmg3d function to get the adjacent tetras to find the valid parent
      DO i=NBulk+1, NBulk + NBdry
        Element => NewMeshR % Elements(i)
        NElNodes = Element % TYPE % NumberOfNodes

        IF(Element % TYPE % ElementCode /= 303) THEN
          ! PRINT *,'Deleting element code ',Element % TYPE % ElementCode
          RmElem(i) = .TRUE.
          CYCLE
        END IF

        ParentElem => Element % BoundaryInfo % Left
        PRINT *,'Parent elem type, body id, bcid: ',ParentElem % TYPE % ElementCode, &
             ParentElem % BodyID, Element % BoundaryInfo % Constraint

        !Not needed
        IF(ParentElem % BodyID == 3 .AND. Element % BoundaryInfo % Constraint /= 10) THEN
          ! PRINT *,'Triangle not needed: ',i,Element % BoundaryInfo % Constraint
          RmElem(i) = .TRUE.

        !Calving front - find the parent elem in the remaining domain:
        ELSEIF(ParentElem % BodyID == 3) THEN
!          Element % BoundaryInfo % Constraint = 1 !TODO <- deal with these properly
          CALL MMG3D_Get_AdjaTet(mmgMesh, ParentElem % ElementIndex, adjList,ier)
          DO j=1,4
            IF(adjlist(j) == 0) CYCLE
            IF(NewMeshR % Elements(adjlist(j)) % BodyID == 2) THEN
              Element % BoundaryInfo % Left => NewMeshR % Elements(adjList(j))
              EXIT
            END IF
          END DO
        END IF
      END DO

      !WORKING HERE - need an efficient way to remove the relevant element structures etc
      ! AND RENUMBER THE NODES to account for deletion.

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
        Element => NewMeshR % Elements(i)
        IF(RmElem(i)) CYCLE
        DO j=1,Element % TYPE % NumberOfNodes
          Element % NodeIndexes(j) = Nodeno_map(Element % NodeIndexes(j))
          IF(Element % NodeIndexes(j) == 0) CALL Fatal(SolverName, &
               "Programming error: mapped nodeno = 0")
        END DO
      END DO

      !QUICK TEST - THIS WILL LEAK MEMORY!!!!
      NewMeshR % Elements = PACK(NewMeshR % Elements, (.NOT. RmElem))
      NewMeshR % Nodes % x = PACK(NewMeshR % Nodes % x, (.NOT. RmNode))
      NewMeshR % Nodes % y = PACK(NewMeshR % Nodes % y, (.NOT. RmNode))
      NewMeshR % Nodes % z = PACK(NewMeshR % Nodes % z, (.NOT. RmNode))
      NewMeshR % NumberOfBulkElements = COUNT(.NOT. RmElem(1:nbulk))
      NewMeshR % NumberOfBoundaryElements = COUNT(.NOT. RmElem(nbulk+1:nbulk+nbdry))
      NewMeshR % NumberOfNodes = COUNT(.NOT. RmNode)

      ! DO i=NewMeshR % NumberOfBulkElements + 1, NewMeshR % NumberOfBulkElements + NewMeshR % NumberOfBoundaryElements
      !   PRINT *,'elem type: ',NewMeshR % Elements(i) % TYPE % ElementCode,ASSOCIATED(NewMesh % Elements(i) % BoundaryInfo)
      !   PRINT *,'constraint: ',NewMesh % Elements(i) % BoundaryInfo % Constraint
      ! END DO
      PRINT *,'NewMeshR counts ',SIZE(NewMeshR % Elements),NewMeshR % NumberOfBulkElements+&
           NewMeshR % NumberOfBoundaryElements, SIZE(NewMeshR % Nodes % X), NewMeshR % NumberOfNodes 
!      CALL WriteMeshToDisk2(Model, NewMeshR, "./outmesh")

    END IF


    CALL MPI_BARRIER(COMM_CALVE,ierr)

    IF(ImBoss) DEALLOCATE(Remesh_nodes % x, Remesh_nodes % y, Remesh_nodes % z)
   CALL MPI_COMM_FREE(COMM_CALVE,ierr)
   CALL MPI_GROUP_FREE(group_world,ierr)
  END IF


  CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)


  !Potentially useful functions
  !CALL SyncNeighbours(ParEnv)
  !CALL SParFaceNumbering, SParEdgeNumbering
  !FindMeshEdges? SetMedgeEdgeFaceDOFs, SetMeshMaxDofs
  !MeshStabParams!

  !CALL InspectMesh(NewMesh)


  !Info contained in an ElmerGrid parallel mesh file:
  ! Node positions
  ! global NN
  ! Shared node info (partitions)  eqv  % NeighbourList(n) % Neighbours
  ! Element NodeNums <- should pass global
  ! Element type
  ! Boundary element parents
  ! Boundary element physical IDs (BC info)
  ! Bulk element bodies 


END SUBROUTINE CalvingRemeshMMG
