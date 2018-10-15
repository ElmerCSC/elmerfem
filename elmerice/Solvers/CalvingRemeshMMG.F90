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

!Remesh the calving model using MMG3D - runs in parallel but remeshing is serial!
!Takes a level set which defines a calving event (or multiple calving events). Level
! set is negative inside a calving event, and positive in the remaining domain. This
! hasn't actually been implemented yet, we use a test function.

! Strategy:
!----------------

! - Use Mesh % Repartition and RedistributeMesh to send relevant (calving)
!   mesh regions to the nominated remeshing partition

! - Remesh with MMG - done serially on nominated partition (TODO - improve nomination of partition)

! - Global node/element number renegotiation

! - Zoltan to compute new partitioning

! - RedistributeMesh back to target processors

! - Interpolate variables etc

! - Continue simulation

SUBROUTINE CalvingRemeshMMG( Model, Solver, dt, Transient )

  USE MeshUtils
  USE CalvingGeometry
  USE MeshPartition
  USE MeshRemeshing

  IMPLICIT NONE

#include "mmg/mmg3d/libmmgtypesf.h"

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !--------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh,GatheredMesh,NewMeshR,FinalMesh
  TYPE(Element_t),POINTER :: Element, ParentElem
  INTEGER :: i,j,k,NNodes,GNBulk, GNBdry, GNNode, NBulk, Nbdry, ierr, &
       my_cboss,MyPE, PEs,CCount, counter, GlNode_min, GlNode_max,adjList(4),front_BC_ID, &
       my_calv_front,calv_front, ncalv_parts, group_calve, comm_calve, group_world,ecode, NElNodes
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, ALLOCATABLE :: Prnode_count(:), cgroup_membs(:),disps(:),elem_fixed(:), &
       PGDOFs_send(:),pcalv_front(:),GtoLNN(:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh, hmin, hmax, hgrad, hausd
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), test_lset(:), Ptest_lset(:), Gtest_lset(:)
  REAL(KIND=dp), POINTER :: PArray(:,:)=> NULL()
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), fixed_node(:), fixed_elem(:), &
       elem_send(:), RmElem(:), RmNode(:)
  LOGICAL :: ImBoss, Found, Debug=.FALSE.
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

  front_BC_id = 1

  hmin = ListGetConstReal(SolverParams, &
       "Mesh Hmin", Found)
  IF(.NOT. Found) hmin = 20.0

  hmax = ListGetConstReal(SolverParams, &
       "Mesh Hmax", Found)
  IF(.NOT. Found) hmax = 4000.0

  hgrad = ListGetConstReal(SolverParams, &
       "Mesh Hgrad", Found)
  IF(.NOT. Found) hgrad = 0.5

  hausd = ListGetConstReal(SolverParams, &
       "Mesh Hausd", Found)
  IF(.NOT. Found) hausd = 20.0

  PArray => ListGetConstRealArray(SolverParams, &
       "Test Point", Found)
  IF(Found) THEN
    test_point = PArray(1:3,1)
  ELSE 
    test_point = (/491600.0, -2290000.0,79.9166641235/)
  END IF

  test_thresh = ListGetConstReal(SolverParams, &
       "Test Calve Dist", Found)
  IF(.NOT. Found) test_thresh = 1500.0

  remesh_thresh = ListGetConstReal(SolverParams, &
       "Test Remesh Dist", Found)
  IF(.NOT. Found) remesh_thresh = 1000.0

  PRINT *,ParEnv % MyPE,' hmin: ',hmin
  PRINT *,ParEnv % MyPE,' hmax: ',hmax
  PRINT *,ParEnv % MyPE,' hausd: ',hausd
  PRINT *,ParEnv % MyPE,' test_point: ',test_point

  ALLOCATE(test_dist(NNodes),&
       test_lset(NNodes),&
       remeshed_node(NNodes),&
       calved_node(NNodes),&
       fixed_node(NNodes)&
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
  ALLOCATE(elem_fixed(nbulk+nbdry), elem_send(nbulk+nbdry))
  elem_fixed = 0 !integer rather than logical to facilitate sending in packed element structure
  elem_send = .FALSE.

  IF(my_calv_front > 0) THEN
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

  IF(Debug) THEN
    PRINT *,mype,' debug calv_front: ',my_calv_front
    PRINT *,mype,' debug calv_fronts: ',pcalv_front
  END IF

  ncalv_parts = COUNT(pcalv_front==calv_front) !only one front for now...
  ALLOCATE(cgroup_membs(ncalv_parts))
  counter = 0
  DO i=1,PEs
    IF(pcalv_front(i) == calv_front) THEN !one calve front
      counter = counter + 1
      cgroup_membs(counter) = i-1
    END IF
  END DO
  IF(Debug) PRINT *,mype,' debug group: ',cgroup_membs


  !Create an MPI_COMM for each calving region, allow gathering instead of 
  !explicit send/receive
  CALL MPI_Comm_Group( ELMER_COMM_WORLD, group_world, ierr)
  CALL MPI_Group_Incl( group_world, ncalv_parts, cgroup_membs, group_calve, ierr)
  CALL MPI_Comm_create( ELMER_COMM_WORLD, group_calve, COMM_CALVE, ierr)
  !--------------------------

  !Work out to whom I send mesh parts for calving
  !TODO - this is currently the lowest partition number which has some calving nodes,
  !but it'd be more efficient to take the partition with the *most* calved nodes
  IF(My_Calv_Front>0) THEN

    my_cboss = MINLOC(pcalv_front, 1, MASK=pcalv_front==calv_front)-1
    IF(Debug) PRINT *,MyPe,' debug calving boss: ',my_cboss
    ImBoss = MyPE == my_cboss

    IF(ImBoss) THEN
      ALLOCATE(Prnode_count(ncalv_parts))
      Prnode_count = 0
    END IF
  ELSE
    ImBoss = .FALSE.
  END IF


  !Redistribute the mesh (i.e. partitioning) so that cboss partitions
  !contain entire calving/remeshing regions.
  IF(.NOT. ASSOCIATED(Mesh % Repartition)) THEN
    ALLOCATE(Mesh % Repartition(NBulk+NBdry), STAT=ierr)
    IF(ierr /= 0) PRINT *,ParEnv % MyPE,' couldnt allocate Mesh % Repartition'
  END IF

  Mesh % Repartition = ParEnv % MyPE + 1
  DO i=1,NBulk+NBdry
    IF(elem_send(i)) Mesh % Repartition(i) = my_cboss+1
  END DO

  GatheredMesh => RedistributeMesh(Model, Mesh, .TRUE., .FALSE.)
  !Confirmed that boundary info is for Zoltan at this point

  IF(Debug) THEN
    PRINT *,ParEnv % MyPE,' gatheredmesh % nonodes: ',GatheredMesh % NumberOfNodes
    PRINT *,ParEnv % MyPE,' gatheredmesh % neelems: ',GatheredMesh % NumberOfBulkElements, &
         GatheredMesh % NumberOfBoundaryElements
  END IF

  !Now we have the gathered mesh, need to send:
  ! - Remeshed_node, test_lset, elem_fixed
  ! - we can convert this into an integer code on elements (0 = leave alone, 1 = remeshed, 2 = fixed)
  ! - or we can simply send test_lset and recompute on calving_boss

  IF(My_Calv_Front>0) THEN
    !root is always zero because cboss is always lowest member of new group (0)
    CALL MPI_Gather(COUNT(remeshed_node), 1, MPI_INTEGER, Prnode_count, 1, &
         MPI_INTEGER, 0, COMM_CALVE,ierr)

    IF(ImBoss) THEN

      IF(Debug) PRINT *,'boss debug prnode_count: ', Prnode_count
      GLNode_max = MAXVAL(GatheredMesh % ParallelInfo % GlobalDOFs)
      GLNode_min = MINVAL(GatheredMesh % ParallelInfo % GlobalDOFs)
      GNBulk = GatheredMesh % NumberOfBulkElements
      GNBdry = GatheredMesh % NumberOfBoundaryElements
      GNNode = GatheredMesh % NumberOfNodes

      ALLOCATE(PGDOFs_send(SUM(Prnode_count)), &
           Ptest_lset(SUM(Prnode_count)), &
           disps(ncalv_parts),GtoLNN(GLNode_min:GLNode_max),&
           gtest_lset(GNNode),&
           fixed_node(GNNode),&
           fixed_elem(GNBulk + GNBdry))

      fixed_node = .FALSE.
      fixed_elem = .FALSE.
      gtest_lset = remesh_thresh + 500.0 !Ensure any far (unshared) nodes are fixed

      !Compute the global to local map
      DO i=1,GNNode
        GtoLNN(GatheredMesh % ParallelInfo % GlobalDOFs(i)) = i
      END DO

    ELSE
      ALLOCATE(disps(1), prnode_count(1))
    END IF

    !Compute the offset in the gathered array from each part
    IF(ImBoss) THEN
      disps(1) = 0
      DO i=2,ncalv_parts
        disps(i) = disps(i-1) + prnode_count(i-1)
      END DO
    END IF

    !Gather the level set function
    CALL MPI_GatherV(PACK(test_lset,remeshed_node), COUNT(remeshed_node), MPI_DOUBLE_PRECISION, &
          Ptest_lset, Prnode_count, disps, MPI_DOUBLE_PRECISION, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    !Gather the GDOFs
    CALL MPI_GatherV(PACK(Mesh % ParallelInfo % GlobalDOFs,remeshed_node), &
         COUNT(remeshed_node), MPI_INTEGER, PGDOFs_send, Prnode_count, &
         disps, MPI_INTEGER, 0, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    !Nodes with levelset value greater than remesh_thresh are kept fixed
    !as are any elements with any fixed_nodes
    !This implies the levelset is a true distance function, everywhere in the domain.
    IF(ImBoss) THEN
      DO i=1,SUM(Prnode_count)
        k = GtoLNN(PGDOFs_send(i))
        IF(k==0) CALL Fatal(SolverName, "Programming error in Global to Local NNum map")
        Gtest_lset(k) = Ptest_lset(i)
      END DO

      fixed_node = Gtest_lset > remesh_thresh

      DO i=1, GNBulk + GNBdry
        Element => GatheredMesh % Elements(i)
        IF(ANY(fixed_node(Element % NodeIndexes))) THEN
          fixed_elem(i) = .TRUE.
        END IF
      END DO
    END IF

  END IF

  !Nominated parition does the remeshing
  IF(ImBoss) THEN

      !Initialise MMG datastructures
      mmgMesh = 0
      mmgSol  = 0

      CALL MMG3D_Init_mesh(MMG5_ARG_start, &
           MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
           MMG5_ARG_end)

      CALL Set_MMG3D_Mesh(GatheredMesh, .TRUE.)

      !Request isosurface discretization
      CALL MMG3D_Set_iparameter(mmgMesh, mmgSol, MMGPARAM_iso, 1,ierr)

      !Set geometric parameters for remeshing
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmin,&
           hmin,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmax,&
           hmax,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hausd,&
           hausd,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hgrad,&
           hgrad,ierr)

      !Feed in the level set distance
      CALL MMG3D_SET_SOLSIZE(mmgMesh, mmgSol, MMG5_Vertex, GNNode ,MMG5_Scalar, ierr)
      DO i=1,GNNode
        CALL MMG3D_Set_scalarSol(mmgSol,&
             Gtest_lset(i), &
             i,ierr)
      END DO

      IF(Debug) THEN
        PRINT *,' boss fixed node: ',COUNT(fixed_node),SIZE(fixed_node)
        PRINT *,' boss fixed elem: ',COUNT(fixed_elem),SIZE(fixed_elem)
      END IF

      !Set required nodes and elements
      DO i=1,GNNode
        IF(fixed_node(i)) THEN
          CALL MMG3D_SET_REQUIREDVERTEX(mmgMesh,i,ierr)
        END IF
      END DO

      DO i=1,GNBulk + GNBdry
        IF(fixed_elem(i)) THEN
          IF(GatheredMesh % Elements(i) % TYPE % DIMENSION == 3) THEN
            CALL MMG3D_SET_REQUIREDTETRAHEDRON(mmgMesh,i,ierr)
          ELSEIF(GatheredMesh % Elements(i) % TYPE % DIMENSION == 2) THEN
            CALL MMG3D_SET_REQUIREDTRIANGLE(mmgMesh,i-GNBulk,ierr)
          END IF
        END IF
      END DO

      !> 4) (not mandatory): check if the number of given entities match with mesh size
      CALL MMG3D_Chk_meshData(mmgMesh,mmgSol,ierr)
      IF ( ierr /= 1 ) CALL EXIT(105)

      !> ------------------------------ STEP  II --------------------------
      !! remesh function
      CALL MMG3D_mmg3dls(mmgMesh,mmgSol,ierr)

      CALL MMG3D_SaveMesh(mmgMesh,"test_out.mesh",LEN(TRIM("test_out.mesh")),ierr)

      CALL Get_MMG3D_Mesh(NewMeshR, .TRUE.)

      NewMeshR % Name = Mesh % Name

      NNodes = NewMeshR % NumberOfNodes
      NBulk = NewMeshR % NumberOfBulkElements
      NBdry = NewMeshR % NumberOfBoundaryElements
      IF(DEBUG) PRINT *, 'NewMeshR nonodes: ',NNodes, NBulk, NBdry
      
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

        IF(Element % TYPE % ElementCode /= 504) &
             PRINT *,'Programming error, bulk type: ',Element % TYPE % ElementCode
        IF(NElNodes /= 4) PRINT *,'Programming error, bulk nonodes: ',NElNodes

        !outbody - mark for deletion and move on
        IF(Element % BodyID == 3) THEN
          RmElem(i) = .TRUE.
          !Deal with an MMG3D eccentricity - returns erroneous GlobalDOF = 10
          !on some split calving elements
          NewMeshR % ParallelInfo % GlobalDOFs(Element % NodeIndexes) = 0
          CYCLE
        END IF

        CALL MMG3D_Get_AdjaTet(mmgMesh, i, adjList,ierr)
        IF(ALL(adjList == 0)) RmElem(i) = .TRUE.

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
          RmElem(i) = .TRUE.
          CYCLE
        END IF

        ParentElem => Element % BoundaryInfo % Left

        IF(ParentElem % BodyID == 3) THEN
          
          !Not needed
          IF(Element % BoundaryInfo % Constraint /= 10) THEN
            !TODO, test constraint == 10 for other BC numbers on front

            RmElem(i) = .TRUE.

          !Switch parent elem to elem in remaining domain
          ELSE
            CALL MMG3D_Get_AdjaTet(mmgMesh, ParentElem % ElementIndex, adjList,ierr)
            DO j=1,4
              IF(adjlist(j) == 0) CYCLE
              IF(NewMeshR % Elements(adjlist(j)) % BodyID == 2) THEN
                Element % BoundaryInfo % Left => NewMeshR % Elements(adjList(j))
                EXIT
              END IF
            END DO
          END IF

        !Edge case - unconnected bulk element
        ELSE IF(RmElem(ParentElem % ElementIndex)) THEN
          RmElem(i) = .TRUE.
        END IF

      END DO

      !Set constraint 10 (the newly formed calving front) to front_BC_id
      !and (temporarily (TODO)) set 0 to 4 so that it doesn't break WriteMeshToDisk2
      !NOTE: I think this will be 10 * the previous BC ID...
      DO i=NBulk+1, NBulk + NBdry
        Element => NewMeshR % Elements(i)
        IF(Element % BoundaryInfo % Constraint == 10) &
             Element % BoundaryInfo % Constraint = front_BC_id

        IF(Element % BoundaryInfo % Constraint == 0) &
             Element % BoundaryInfo % Constraint = 4
      END DO
     

      !TODO - issue - MMG3Dls will mark new nodes 'required' if they split a previous
      !boundary element. but only *front* boundary elements? Have confirmed it ONLY returns
      !required & 10 for calving front nodes, and it's not to do with manually setting required stuff.
      !But, the model does complain about required entities if the hole is internal, but no weird 
      !GIDs are returned.
      IF(Debug) THEN
        DO i=1,GatheredMesh % NumberOfNodes
          PRINT *, ParEnv % MyPE,' debug old ',i,&
               ' GDOF: ',GatheredMesh % ParallelInfo % GlobalDOFs(i),&
               ' xyz: ',&
               GatheredMesh % Nodes % x(i),&
               GatheredMesh % Nodes % y(i),&
               GatheredMesh % Nodes % z(i),fixed_node(i)
          IF(fixed_node(i)) THEN
            IF(.NOT. ANY(NewMeshR % ParallelInfo % GlobalDOFs == &
                 GatheredMesh % ParallelInfo % GlobalDOFs(i))) CALL Fatal(SolverName, &
                 "Missing required node in output!")
          END IF
        END DO
      END IF

      !Chop out the flagged elems and nodes
      CALL CutMesh(NewMeshR, RmNode, RmElem)

      IF(Debug) THEN
        DO i=1,NewMeshR % NumberOfNodes
          PRINT *, ParEnv % MyPE,' debug new ',i,&
               ' GDOF: ',NewMeshR % ParallelInfo % GlobalDOFs(i),&
               ' xyz: ',&
               NewMeshR % Nodes % x(i),&
               NewMeshR % Nodes % y(i),&
               NewMeshR % Nodes % z(i)

          IF(NewMeshR % ParallelInfo % GlobalDOFs(i) == 0) CYCLE
          IF(.NOT. ANY(GatheredMesh % ParallelInfo % GlobalDOFs == &
               NewMeshR % ParallelInfo % GlobalDOFs(i))) CALL Warn(SolverName, &
               "Unexpected GID")
        END DO
      END IF

      !Update parallel info from old mesh nodes (shared node neighbours)
      CALL MapNewParallelInfo(GatheredMesh, NewMeshR)

      CALL ReleaseMesh(GatheredMesh)
      GatheredMesh => NewMeshR
      NewMeshR => NULL()
   END IF

   !Wait for all partitions to finish
   IF(My_Calv_Front>0) THEN
     CALL MPI_BARRIER(COMM_CALVE,ierr)
     CALL MPI_COMM_FREE(COMM_CALVE,ierr)
     CALL MPI_GROUP_FREE(group_world,ierr)
   END IF
   CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)

   !Now each partition has GatheredMesh, we need to renegotiate globalDOFs
   CALL RenumberGDOFs(Mesh, GatheredMesh)

   !and global element numbers
   CALL RenumberGElems(GatheredMesh)

   !Some checks on the new mesh
   !----------------------------
   DO i=1,GatheredMesh % NumberOfNodes
     IF(GatheredMesh % ParallelInfo % INTERFACE(i)) THEN
       IF(.NOT. ASSOCIATED(GatheredMesh % ParallelInfo % Neighbourlist(i) % Neighbours)) &
            CALL Fatal(SolverName, "Neighbourlist not associated!")
       IF(SIZE(GatheredMesh % ParallelInfo % Neighbourlist(i) % Neighbours) < 2) &
            CALL Fatal(SolverName, "Neighbourlist too small!")
     END IF
   END DO

   DO i=1,GatheredMesh % NumberOfNodes
     IF(.NOT. ASSOCIATED(GatheredMesh % ParallelInfo % NeighbourList(i) % Neighbours)) &
          CALL Fatal(SolverName, 'Unassociated Neighbourlist % Neighbours')
     IF(GatheredMesh % ParallelInfo % GlobalDOFs(i) == 0) &
          CALL Fatal(SolverName, 'Bad GID 0')
   END DO

   IF(ANY(GatheredMesh % Elements % GElementIndex <= 0)) CALL Fatal(SolverName, 'Element has ID 0')


   !Call zoltan to determine redistribution of mesh
   ! then do the redistribution
   !-------------------------------
   CALL Zoltan_Interface( Model, GatheredMesh )

   FinalMesh => RedistributeMesh(Model, GatheredMesh, .TRUE., .FALSE.)

   Model % Mesh => FinalMesh
   Solver % Mesh => FinalMesh
   Model % Meshes => FinalMesh

   CALL ReleaseMesh(GatheredMesh)

   !TODO here - call something like SwitchMesh from CalvingRemesh.F90 to 
   !handle the interpolation of the variables, the reconstruction of the
   !solver matrices etc.

END SUBROUTINE CalvingRemeshMMG
