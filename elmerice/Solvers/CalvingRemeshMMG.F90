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

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !--------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh,NewMesh
  TYPE(Element_t),POINTER :: Element
  TYPE(Element_t), ALLOCATABLE :: PElements(:)
  TYPE(Nodes_t) :: Remesh_nodes
  INTEGER :: i, NNodes,NBulks, Nbdry, ierr, my_cboss,MyPE, PEs,CCount, counter, &
       my_calv_front,calv_front, ncalv_parts, group_calve, comm_calve, group_world,ecode, NElNodes,&
       NBulk_send, NBdry_send, NBulk_fixed, NBdry_fixed,test_min, test_max
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, ALLOCATABLE :: Pcnode_count(:), Prnode_count(:), cgroup_membs(:),&
       PNVerts(:),PNTris(:), PNTetras(:),FixedElems(:), GDOFs_send(:),PGDOFs_send(:),&
       el_stream(:), pel_stream(:), pnbulk_send(:), pnbdry_send(:), pnbulk_fixed(:), &
       pnbdry_fixed(:),pcalv_front(:),PElem_parts(:),&
       disps(:),el_stream_sizes(:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), NodeCoords_send(:),PNodeCoords_send(:)
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), elem_fixed(:),&
       elem_send(:), bulk_fixed(:), bdry_fixed(:),bulk_send(:), bdry_send(:)
  LOGICAL :: ImBoss
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
  SolverParams => GetSolverParams()
  SolverName = "CalvingRemeshMMG"

  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes
  NBulks = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements

  calv_front = 1
  MyPe = ParEnv % MyPE
  PEs = ParEnv % PEs

  !Check all elements are tetras or tris:
  DO i=1,NBulks+Nbdry
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
  remesh_thresh = 10000.0

  ALLOCATE(test_dist(NNodes),&
       remeshed_node(NNodes),&
       calved_node(NNodes),&
       )

  test_dist = 0.0
  remeshed_node = .FALSE.
  calved_node = .FALSE.

  !Compute involved nodes
  DO i=1, NNodes
    test_dist(i) = (((Mesh % Nodes % x(i) - test_point(1))**2.0) + &
         ((Mesh % Nodes % y(i) - test_point(2))**2.0) + &
         ((Mesh % Nodes % z(i) - test_point(3))**2.0)) ** 0.5
  END DO

  !This does nothing yet but it will be important - determine
  !the discrete calving zones, each of which will be seperately remeshed
  !by a nominated boss partition
!  CALL CountCalvingEvents(Model, Mesh, ccount)
  ccount = 1

  calved_node = test_dist < test_thresh
  remeshed_node = test_dist < (test_thresh + remesh_thresh)

  my_calv_front = 0
  IF(ANY(remeshed_node)) THEN
    my_calv_front = 1 !this is integer (not logical) so we can later move to multiple calving fronts
  END IF

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
      ALLOCATE(Pcnode_count(ncalv_parts), Prnode_count(ncalv_parts))
      Pcnode_count = 0
      Prnode_count = 0
    END IF

    !Gather info from partitions which are involved
    !-----------------------

    !root is always zero because cboss is always lowest member of new group (0)
    CALL MPI_Gather(COUNT(calved_node), 1, MPI_INTEGER, Pcnode_count, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr) 

    CALL MPI_Gather(COUNT(remeshed_node), 1, MPI_INTEGER, Prnode_count, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)

    IF(ImBoss) THEN
      ALLOCATE(PNverts(ncalv_parts),&
           PNtetras(ncalv_parts),&
           PNtris(ncalv_parts),&
           PNBulk_fixed(ncalv_parts),&
           PNBulk_send(ncalv_parts),&
           PNBdry_fixed(ncalv_parts),&
           PNBdry_send(ncalv_parts),&
           PGDOFs_send(SUM(Prnode_count)), &
           PNodeCoords_send(SUM(Prnode_count)*3), &
           disps(ncalv_parts)&
      )
    ELSE
      ALLOCATE(disps(1), prnode_count(1))
    END IF

    CALL PackNodesToSend(Mesh, remeshed_node, GDOFs_send, NodeCoords_send, 3)

    IF(ImBoss) THEN
      disps(1) = 0
      DO i=2,ncalv_parts
        disps(i) = disps(i-1) + prnode_count(i-1)
      END DO
      PRINT *,mype,' size gdofs: ',SIZE(GDOFs_send),'nodes count: ',SIZE(nodecoords_send)
      PRINT *,mype,' disps*3: ',disps*3,' pnodecoord_send: ', SIZE(pnodecoords_send)
    END IF

    CALL MPI_GatherV(GDOFs_send, COUNT(remeshed_node), MPI_INTEGER, PGDOFs_send, Prnode_count, &
         disps, MPI_INTEGER, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    CALL MPI_GatherV(NodeCoords_send, COUNT(remeshed_node)*3, MPI_DOUBLE_PRECISION, PNodeCoords_send,&
         Prnode_count*3, disps*3, MPI_DOUBLE_PRECISION, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    PRINT *,mype,'gathered' 

    IF(ImBoss) CALL UnpackNodesSent(PGDOFs_send,PNodeCoords_send, Remesh_Nodes, DIM=3)

    IF(ImBoss) PRINT *,'debug got nodes: ',SIZE(Remesh_Nodes % x),SIZE(PGDOfs_send),MINVAL(PGDOfs_send),MAXVAL(PGDOfs_send)

    PRINT *,mype,' gathered 1'

    ALLOCATE(elem_fixed(nbulks+nbdry), elem_send(nbulks+nbdry))
    elem_fixed = .FALSE.
    elem_send = .FALSE.

    !Each partition identifies (based on nodes), elements which need to be transferred
    DO i=1,NBulks+Nbdry
      Element => Mesh % Elements(i)
      NodeIndexes => Element % NodeIndexes
      NElNodes = Element % TYPE % NumberOfNodes
      IF(ANY(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
        elem_send(i) = .TRUE.
        IF(ALL(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
          elem_fixed(i) = .TRUE.
        END IF
      END IF
    END DO
    nbulk_send = COUNT(elem_send(1:nbulks))
    nbdry_send = COUNT(elem_send(nbulks+1:nbulks+nbdry))
    nbulk_fixed = COUNT(elem_fixed(1:nbulks))
    nbdry_fixed = COUNT(elem_fixed(nbulks+1:nbulks+nbdry))

    PRINT *,mype,' about to pack'

    CALL PackElemsToSend(Mesh, elem_send, el_stream)

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

    IF(ImBoss) CALL UnpackElemsSent(PEl_stream, PElements,  3, PElem_parts)

    IF(ImBoss) THEN
      test_max = 0
      test_min = HUGE(test_min)
      PRINT *,' DEBUG SIZE PELEMENTS: ',SIZE(PElements)
      DO i=1,SIZE(PElements)
        IF(.NOT. ASSOCIATED(PElements(i) % NodeIndexes)) &
             CALL Fatal(SolverName,"Programm error: unassociated elements")
        test_max = MAX(MAXVAL(PElements(i) % NodeIndexes),test_max)
        test_min = MIN(MINVAL(PElements(i) % NodeIndexes),test_min)
      END DO
      PRINT *,' debug max min elem nodeindices: ',test_min, test_max
    END IF
    PRINT *,mype,' packed ',SIZE(elem_send)

    !and send them to boss
    CALL MPI_Gather(nbulk_send, 1, MPI_INTEGER, pnbulk_send, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)
    CALL MPI_Gather(nbdry_send, 1, MPI_INTEGER, pnbdry_send, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)
    CALL MPI_Gather(nbulk_fixed, 1, MPI_INTEGER, pnbulk_fixed, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)
    CALL MPI_Gather(nbdry_fixed, 1, MPI_INTEGER, pnbdry_fixed, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)

    PRINT *,mype,' gathered counts'

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
