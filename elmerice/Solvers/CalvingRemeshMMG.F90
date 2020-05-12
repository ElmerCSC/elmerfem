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
  TYPE(Variable_t), POINTER :: CalvingVar
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh,GatheredMesh,NewMeshR,NewMeshRR,FinalMesh
  TYPE(Element_t),POINTER :: Element, ParentElem
  INTEGER :: i,j,k,NNodes,GNBulk, GNBdry, GNNode, NBulk, Nbdry, ierr, &
       my_cboss,MyPE, PEs,CCount, counter, GlNode_min, GlNode_max,adjList(4),front_BC_ID, &
       my_calv_front,calv_front, ncalv_parts, group_calve, comm_calve, group_world,ecode, NElNodes,&
       target_bodyid,gdofs(4)
  INTEGER, POINTER :: NodeIndexes(:), geom_id
  INTEGER, ALLOCATABLE :: Prnode_count(:), cgroup_membs(:),disps(:), &
       PGDOFs_send(:),pcalv_front(:),GtoLNN(:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh, hmin, hmax, hgrad, hausd
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), test_lset(:), Ptest_lset(:), Gtest_lset(:), &
       target_length(:,:)
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), fixed_node(:), fixed_elem(:), &
       elem_send(:), RmElem(:), RmNode(:),new_fixed_node(:), new_fixed_elem(:)
  LOGICAL :: ImBoss, Found, Isolated, Debug=.TRUE.
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

  !Get Body ID
  target_bodyid = Mesh % Elements(1) % BodyID
  IF(target_bodyid /= 1) CALL Warn(SolverName, "Body ID is not 1, this case might not be well handled")
  !For testing - set a calving levelset function to mimic calving events
  !-------------------

  !TODO - unhardcode (detect?) this
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

  remesh_thresh = ListGetConstReal(SolverParams, &
       "Test Remesh Dist", Found)
  IF(.NOT. Found) remesh_thresh = 1000.0

  PRINT *,ParEnv % MyPE,' hmin: ',hmin
  PRINT *,ParEnv % MyPE,' hmax: ',hmax
  PRINT *,ParEnv % MyPE,' hausd: ',hausd
  PRINT *,ParEnv % MyPE,' remesh_thresh: ',remesh_thresh


  !Get the calving levelset function (-ve inside calving event, +ve in intact ice)
  !-------------------
  CalvingVar => VariableGet(Mesh % Variables, "Calving Lset", .TRUE., UnfoundFatal=.TRUE.)

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

  test_lset = CalvingVar % Values(CalvingVar % Perm(:)) !TODO - quick&dirty, possibly zero perm?
  calved_node = test_lset < 0.0
  remeshed_node = test_lset < remesh_thresh

  my_calv_front = 0
  IF(ANY(remeshed_node)) THEN
    my_calv_front = 1 !this is integer (not logical) so we can later move to multiple calving fronts
  END IF

  !TODO - could make this more efficient by cutting out some of the elements from the calved region.
  !This region will be remeshed too, but we discard it, so the closer we can get to the edge of the
  !calving event, the better.

  !Identify all elements which need to be sent (including those fixed in place)
  !Here we also find extra nodes which are just beyond the remeshing threshold
  !but which are in fixed elements, thus also need to be sent
  ALLOCATE(elem_send(nbulk+nbdry))
  elem_send = .FALSE.

  IF(my_calv_front > 0) THEN
    !Each partition identifies (based on nodes), elements which need to be transferred
    DO i=1,Nbulk
      Element => Mesh % Elements(i)
      CALL Assert(Element % ElementIndex == i, SolverName,"Misnumbered bulk element.")
      NodeIndexes => Element % NodeIndexes
      NElNodes = Element % TYPE % NumberOfNodes
      IF(ANY(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
        elem_send(i) = .TRUE.
        IF(.NOT. ALL(remeshed_node(NodeIndexes(1:NElNodes)))) THEN
          fixed_node(NodeIndexes(1:NElNodes)) = .TRUE.
        END IF
      END IF
    END DO

    remeshed_node = remeshed_node .OR. fixed_node

    !Cycle boundary elements, checking parent elems
    !BC elements follow parents
    DO i=NBulk+1, NBulk+NBdry
      Element => Mesh % Elements(i)
      ParentElem => Element % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(ParentElem)) THEN
        ParentElem => Element % BoundaryInfo % Right
      END IF
      CALL Assert(ASSOCIATED(ParentElem),SolverName,"Boundary element has no parent!")
      IF(elem_send(ParentElem % ElementIndex)) elem_send(i) = .TRUE.
    END DO

    DEALLOCATE(fixed_node) !reuse later
  END IF


  !This does nothing yet but it will be important - determine
  !the discrete calving zones, each of which will be separately remeshed
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
    IF(ierr /= 0) PRINT *,ParEnv % MyPE,' could not allocate Mesh % Repartition'
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
  ! - Remeshed_node, test_lset
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

  !Nominated partition does the remeshing
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

      CALL Get_MMG3D_Mesh(NewMeshR, .TRUE., new_fixed_node, new_fixed_elem)

      NewMeshR % Name = Mesh % Name

      NNodes = NewMeshR % NumberOfNodes
      NBulk = NewMeshR % NumberOfBulkElements
      NBdry = NewMeshR % NumberOfBoundaryElements
      IF(DEBUG) PRINT *, 'NewMeshR nonodes, nbulk, nbdry: ',NNodes, NBulk, NBdry
      
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
        ELSE IF(Element % BodyID == 2) THEN
          Element % BodyID = target_bodyid
        ELSE
          PRINT *,'Erroneous body id: ',Element % BodyID
          CALL Fatal(SolverName, "Bad body id!")
        END IF

        !Get rid of any isolated elements (I think?)
        ! This may not be necessary, if the calving levelset is well defined
        CALL MMG3D_Get_AdjaTet(mmgMesh, i, adjList,ierr)
        IF(ALL(adjList == 0)) THEN

          !check if this is truly isolated or if it's at a partition boundary
          isolated = .TRUE.
          gdofs = NewMeshR % ParallelInfo % GlobalDOFs(Element % NodeIndexes(1:NELNodes))
          DO j=1,GatheredMesh % NumberOfNodes
            IF(ANY(gdofs == GatheredMesh % ParallelInfo % GlobalDOFs(j)) .AND. &
                 GatheredMesh % ParallelInfo % INTERFACE(j)) THEN
              isolated = .FALSE.
              EXIT
            END IF
          END DO

          IF(isolated) THEN
            RmElem(i) = .TRUE.
            CALL WARN(SolverName, 'Found an isolated body element.')
          END IF
        END IF

        !Mark all nodes in valid elements for keeping
        RmNode(Element % NodeIndexes(1:NElNodes)) = .FALSE.
      END DO

      !Same thru the boundary elements
      !If their parent elem is body 3, delete, *unless* they have constraint 10 (the calving face)
      !in which case use mmg3d function to get the adjacent tetras to find the valid parent
      DO i=NBulk+1, NBulk + NBdry
        Element => NewMeshR % Elements(i)
        NElNodes = Element % TYPE % NumberOfNodes

        !Get rid of non-triangles
        IF(Element % TYPE % ElementCode /= 303) THEN
          RmElem(i) = .TRUE.
          CYCLE
        END IF

        !Get rid of boundary elements without BCID (newly created internal)
        IF(Element % BoundaryInfo % Constraint == 0) THEN
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
              IF(NewMeshR % Elements(adjlist(j)) % BodyID == target_bodyid) THEN
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
      !NOTE: I think this will be 10 * the previous BC ID...
      !And set all BC BodyIDs based on constraint
      DO i=NBulk+1, NBulk + NBdry

        IF(RmElem(i)) CYCLE

        Element => NewMeshR % Elements(i)
        geom_id => Element % BoundaryInfo % Constraint
        IF(geom_id == 10) geom_id = front_BC_id

        CALL Assert((geom_id > 0) .AND. (geom_id <= Model % NumberOfBCs),&
             SolverName,"Unexpected BC element body id!")

        Element % BodyId  = ListGetInteger( &
             Model % BCs(geom_id) % Values, 'Body Id', Found, 1, Model % NumberOfBodies )
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

      new_fixed_elem = PACK(new_fixed_elem, .NOT. RmElem)
      new_fixed_node = PACK(new_fixed_node, .NOT. RmNode)

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

      !Anisotropic mesh improvement following calving cut:
      !TODO - unhardcode! Also this doesn't work very well yet.
      ALLOCATE(target_length(NewMeshR % NumberOfNodes,3))
      target_length(:,1) = 300.0
      target_length(:,2) = 300.0
      target_length(:,3) = 50.0

      CALL RemeshMMG3D(NewMeshR, target_length, NewMeshRR, new_fixed_node, new_fixed_elem)

      !Update parallel info from old mesh nodes (shared node neighbours)
      CALL MapNewParallelInfo(GatheredMesh, NewMeshRR)

      CALL ReleaseMesh(GatheredMesh)
      CALL ReleaseMesh(NewMeshR)
      GatheredMesh => NewMeshRR
      NewMeshR => NULL()
      NewMeshRR => NULL()
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

   !Call zoltan to determine redistribution of mesh
   ! then do the redistribution
   !-------------------------------

   CALL Zoltan_Interface( Model, GatheredMesh )

   FinalMesh => RedistributeMesh(Model, GatheredMesh, .TRUE., .FALSE.)

   CALL CheckMeshQuality(FinalMesh)

   FinalMesh % Name = TRIM(Mesh % Name)
   FinalMesh % OutputActive = .TRUE.
   FinalMesh % Changed = .TRUE. 

   !Actually switch the model's mesh
   CALL SwitchMesh(Model, Solver, Mesh, FinalMesh)

   !After SwitchMesh because we need GroundedMask
   CALL EnforceGroundedMask(Mesh)

   !Recompute mesh bubbles etc
   CALL MeshStabParams( Model % Mesh )

   !Release the old mesh
   CALL ReleaseMesh(GatheredMesh)

CONTAINS

  SUBROUTINE CheckMeshQuality(Mesh)

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t),POINTER :: Element, Parent
    REAL(KIND=dp) :: U,V,W,detJ,Basis(10), mean 
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,n,l,k, count
    LOGICAL :: stat
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="CheckMeshQuality"

    n = Mesh % MaxElementNodes
    ALLOCATE(ElementNodes % x(n),&
         ElementNodes % y(n),&
         ElementNodes % z(n))

    !Check all BC elements have parents
    DO i=Mesh % NumberOfBulkElements+1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)
      Parent => Element % BoundaryInfo % Left
      IF( .NOT. ASSOCIATED(Parent) ) THEN
        Parent => Element % BoundaryInfo % Right
      END IF
      IF( .NOT. ASSOCIATED( Parent ) ) THEN
        PRINT *,ParEnv % MyPE,i,' BC element without parent! constraint: ',Element % BoundaryInfo % constraint, &
             ' body id: ',Element % BodyID,' nodes: ',Element % NodeIndexes,&
             ' global: ',Mesh % ParallelInfo % GlobalDOFs(Element%NodeIndexes)
        CALL Fatal(SolverName, "BC Element without parent!")
      END IF
    END DO

    !check for duplicate element & node indices (locally only)
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      IF(Mesh % Elements(i) % GElementIndex <= 0) CALL Fatal(SolverName, 'Element has ID 0')
      DO j=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        IF(i==j) CYCLE
        IF(Mesh % Elements(i) % GElementIndex == Mesh % Elements(j) % GElementIndex) THEN
          PRINT *,ParEnv % MyPE,' elements ',i,j,' have same GElementIndex: ',&
               Mesh % Elements(j) % GElementIndex
          CALL Fatal(SolverName, "Duplicate GElementIndexes!")
        END IF
      END DO
    END DO

    DO i=1,Mesh % NumberOfNodes
      IF(Mesh % ParallelInfo % GlobalDOFs(i) <= 0) CALL Fatal(SolverName, 'Node has ID 0')
      DO j=1,Mesh % NumberOfNodes
        IF(i==j) CYCLE
        IF(Mesh % ParallelInfo % GlobalDOFs(i) == Mesh % ParallelInfo % GlobalDOFs(j)) THEN
          PRINT *,ParEnv % MyPE,' nodes ',i,j,' have same GlobalDOF: ',&
               Mesh % ParallelInfo % GlobalDOFs(j)
          CALL Fatal(SolverName, "Duplicate GlobalDOFs!")
        END IF
      END DO
    END DO

    !Check element detj etc
    DO j=1,2
      IF(j==1) mean = 0.0
      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        !Check element for duplicate node indexes
        DO k=1,n
          DO l=1,n
            IF(l==k) CYCLE
            IF(NodeIndexes(k) == NodeIndexes(l)) THEN
              WRITE(Message, '(A,i0,A)') "Mesh Element ",i," has duplicate node indexes!"
              CALL Fatal(FuncName,Message)
            END IF
          END DO
        END DO

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

        stat = ElementInfo( Element,ElementNodes,U,V,W,detJ, &
             Basis )
        !Check detj - warn if deviates from mean, fatal if <= 0
        IF(j==2) THEN
          IF(detj <= 0.0) THEN
            WRITE(Message, '(A,i0,A)') "Element ",j," has detj <= 0" 
            CALL Fatal(FuncName, Message)
          ELSE IF(detj < mean/10.0 .OR. detj > mean*10.0) THEN
            WRITE(Message, '(i0,A,i0,A,F10.2,A,F10.2,A)') ParEnv % MyPE,' element ',&
                 i,' detj (',detj,') deviates from mean (',mean,')'
            CALL Warn(FuncName, Message)
          END IF
        ELSE
          mean = mean + detj
        END IF
      END DO
      IF(j==1) mean = mean / Mesh % NumberOfBulkElements
    END DO




  END SUBROUTINE CheckMeshQuality

  !Takes a mesh with GroundedMask defined on the base, and
  !ensures that grounded nodes remain grounded
  !i.e. sets z = min zs bottom wherever GroundedMask>-0.5
  SUBROUTINE EnforceGroundedMask(Mesh)
    TYPE(Mesh_t), POINTER :: Mesh
    !-------------------------
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Variable_t), POINTER :: GMaskVar
    REAL(KIND=dp), POINTER :: GMask(:)
    REAL(KIND=dp) :: zb
    INTEGER :: i,n
    INTEGER, POINTER :: GMaskPerm(:)
    CHARACTER(MAX_NAME_LEN) :: FuncName="EnforceGroundedMask", GMaskVarName

    GMaskVarName = "GroundedMask"
    GMaskVar => VariableGet(Mesh % Variables, GMaskVarName, .TRUE.)
    IF(.NOT.ASSOCIATED(GMaskVar)) THEN
      CALL Info(FuncName, "Didn't find GroundedMask, so not enforcing bed height",Level=5)
      RETURN
    END IF

    Material => GetMaterial(Mesh % Elements(1)) !TODO, this is not generalised

    GMask => GMaskVar % Values
    GMaskPerm => GMaskVar % Perm

    DO i=1,Mesh % NumberOfNodes
      IF(GMaskPerm(i) == 0) CYCLE
      zb = ListGetRealAtNode(Material, "Min Zs Bottom",i,UnfoundFatal=.TRUE.)

      !Floating -> check no penetration
      !Grounded -> set to bedrock height
      IF(GMask(GMaskPerm(i)) < -0.5) THEN
        IF(Mesh % Nodes % z(i) < zb) Mesh % Nodes % z(i) = zb
      ELSE
        Mesh % Nodes % z(i) = zb
      END IF
    END DO

  END SUBROUTINE EnforceGroundedMask

END SUBROUTINE CalvingRemeshMMG
