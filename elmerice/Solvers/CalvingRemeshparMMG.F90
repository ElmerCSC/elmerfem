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

SUBROUTINE CalvingRemeshParMMG( Model, Solver, dt, Transient )

    USE MeshUtils
    USE CalvingGeometry
    USE MeshPartition
    USE MeshRemeshing

    IMPLICIT NONE

#include "parmmg/libparmmgtypesf.h"

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !--------------------------------------
  TYPE(Variable_t), POINTER :: CalvingVar,DistanceVar
  TYPE(ValueList_t), POINTER :: SolverParams, MeshParams
  TYPE(Mesh_t),POINTER :: Mesh,GatheredMesh,NewMeshR,NewMeshRR,FinalMesh, DistributedMesh
  TYPE(Element_t),POINTER :: Element, ParentElem
  INTEGER :: i,j,k,NNodes,GNBulk, GNBdry, GNNode, NBulk, Nbdry, ierr, &
       my_cboss,MyPE, PEs,CCount, counter, GlNode_min, GlNode_max,adjList(4),&
       front_BC_ID, front_body_id, my_calv_front,calv_front, ncalv_parts, &
       group_calve, comm_calve, group_world,ecode, NElNodes, target_bodyid,gdofs(4), &
       PairCount,RPairCount, NCalvNodes, croot, nonCalvBoss, RNNodes, RNElems
  INTEGER, POINTER :: NodeIndexes(:), geom_id
  INTEGER, ALLOCATABLE :: Prnode_count(:), cgroup_membs(:),disps(:), &
       PGDOFs_send(:),pcalv_front(:),GtoLNN(:),EdgePairs(:,:),REdgePairs(:,:), ElNodes(:),&
       Nodes(:), IslandCounts(:), pNCalvNodes(:,:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh, hmin, hmax, hgrad, hausd, newdist
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), test_lset(:), Ptest_lset(:), Gtest_lset(:), &
       target_length(:,:), Ptest_dist(:), Gtest_dist(:)
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), fixed_node(:), fixed_elem(:), &
       elem_send(:), RmElem(:), RmNode(:),new_fixed_node(:), new_fixed_elem(:), FoundNode(:,:), &
       UsedElem(:), NewNodes(:), RmIslandNode(:), RmIslandElem(:)
  LOGICAL :: ImBoss, Found, Isolated, Debug,DoAniso,NSFail,CalvingOccurs,&
       RemeshOccurs,CheckFlowConvergence, HasNeighbour, Distributed
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, CalvingVarName
  TYPE(Variable_t), POINTER :: TimeVar
  INTEGER :: Time, remeshtimestep, proc, idx, island, node
  REAL(KIND=dp) :: TimeReal, PreCalveVolume, PostCalveVolume, CalveVolume

  SolverParams => GetSolverParams()
  SolverName = "CalvingRemeshParMMG"
  Debug=.TRUE.
  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes
  NBulk = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements

  calv_front = 1
  MyPe = ParEnv % MyPE
  PEs = ParEnv % PEs

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  TimeReal = TimeVar % Values(1)
  Time = INT(TimeReal)

  !for first time step calculate mesh volume
  IF(Time == 1) THEN
   CALL MeshVolume(Mesh, .TRUE., PreCalveVolume)
   IF(MyPe == 0) PRINT*, 'First timestep precalve volume = ', PreCalveVolume
  END IF

  !Check all elements are tetras or tris:
  DO i=1,NBulk+Nbdry
    Element => Mesh % Elements(i)
    ecode = Element % TYPE % ElementCode
    IF(ecode /= 101 .AND. ecode /= 202 .AND. ecode /= 303 .AND. ecode /= 504) THEN
      PRINT *,MyPE,' has unsupported element type: ',ecode
      CALL Fatal(SolverName, "Invalid Element Type!")
    END IF
  END DO

  !Get main Body ID
  target_bodyid = Mesh % Elements(1) % BodyID
  IF(target_bodyid /= 1) CALL Warn(SolverName, "Body ID is not 1, this case might not be well handled")

  !TODO - unhardcode (detect?) this
  front_BC_id = 1
  front_body_id =  ListGetInteger( &
       Model % BCs(front_bc_id) % Values, 'Body Id', Found, 1, Model % NumberOfBodies )

  hmin = ListGetConstReal(SolverParams, "Mesh Hmin",  DefValue=20.0_dp)
  hmax = ListGetConstReal(SolverParams, "Mesh Hmax",  DefValue=4000.0_dp)
  hgrad = ListGetConstReal(SolverParams,"Mesh Hgrad", DefValue=0.5_dp)
  hausd = ListGetConstReal(SolverParams, "Mesh Hausd",DefValue=20.0_dp)
  remesh_thresh = ListGetConstReal(SolverParams,"Remeshing Distance", DefValue=1000.0_dp)
  CalvingVarName = ListGetString(SolverParams,"Calving Variable Name", DefValue="Calving Lset")
  remeshtimestep = ListGetInteger(SolverParams,"Remesh TimeStep", DefValue=2)

  IF(ParEnv % MyPE == 0) THEN
    PRINT *,ParEnv % MyPE,' hmin: ',hmin
    PRINT *,ParEnv % MyPE,' hmax: ',hmax
    PRINT *,ParEnv % MyPE,' hgrad: ',hgrad
    PRINT *,ParEnv % MyPE,' hausd: ',hausd
    PRINT *,ParEnv % MyPE,' remeshing distance: ',remesh_thresh
    PRINT *,ParEnv % MyPE,' remeshing every ', remeshtimestep
  END IF

  !If FlowSolver failed to converge (usually a result of weird mesh), large unphysical
  !calving events can be predicted. So, turn off CalvingOccurs, and ensure a remesh
  !Also undo this iterations mesh update
  NSFail = ListGetLogical(Model % Simulation, "Flow Solution Failed",CheckFlowConvergence)
  IF(CheckFlowConvergence) THEN
    IF(NSFail) THEN
      CalvingOccurs = .FALSE.
      RemeshOccurs = .TRUE.
      CALL Info(SolverName, "Remeshing but not calving because NS failed to converge.")
    END IF
  ELSE
     CalvingOccurs=.TRUE.
     RemeshOccurs=.TRUE.
  END IF

  ALLOCATE(remeshed_node(NNodes),&
       fixed_node(NNodes))
  remeshed_node = .FALSE.
  fixed_node = .FALSE.

  ! TO DO some other wah to define remeshed nodes
  DistanceVar  => VariableGet(Mesh % Variables, "Distance", .TRUE., UnfoundFatal=.TRUE.)
  ALLOCATE(test_dist(NNodes))
  test_dist = DistanceVar  % Values(DistanceVar % Perm(:))
  remeshed_node = test_dist < remesh_thresh

  !Get the calving levelset function (-ve inside calving event, +ve in intact ice)
  !-------------------
  IF (CalvingOccurs) THEN
     CalvingVar => VariableGet(Mesh % Variables, "Calving Lset", .TRUE., UnfoundFatal=.TRUE.)

     ALLOCATE(test_lset(NNodes),&
          calved_node(NNodes)&
          )

     calved_node = .FALSE.
     test_lset = CalvingVar % Values(CalvingVar % Perm(:)) !TODO - quick&dirty, possibly zero perm?
     calved_node = test_lset < 0.0
     ! calving front boundary nodes may have lset value greater than remesh dist
     DO i=1, NNodes
        newdist = MINVAL((/test_lset(i), test_dist(i)/))
        remeshed_node(i) = newdist < remesh_thresh
     END DO
  END IF

  my_calv_front = 0
  IF(ANY(remeshed_node)) THEN
    my_calv_front = 1 !this is integer (not logical) so we can later move to multiple calving fronts
  END IF

  NCalvNodes = COUNT(remeshed_node)

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

  ! loop to allow negotiation for multiple calving fronts
  ! This communication allows the determination of which part of the mesh
  ! has the most calvnodes and will become cboss
  ALLOCATE(pNCalvNodes(MAXVAL(pcalv_front),PEs))
  DO i=1, MAXVAL(pcalv_front)
    CALL MPI_ALLGATHER(NCalvNodes, 1, MPI_INTEGER, pNCalvNodes(i,:), &
    1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)
  END DO

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
  !Now cboss set to part with most calving nodes
  !NonCalvBoss will take any non calving nodes from the cboss
  ! this is still *TODO*
  IF(My_Calv_Front>0) THEN

    ! assume currently only one calving front
    !my_cboss = MINLOC(pNCalvNodes(1,:), 1, MASK=pNCalvNodes(1,:)==MAXVAL(pNCalvNodes(1,:)))-1
    my_cboss = 0
    nonCalvBoss = MINLOC(pNCalvNodes(1,:), 1, MASK=pNCalvNodes(1,:)==MINVAL(pNCalvNodes(1,:)))-1
    IF(Debug) PRINT *,MyPe,' debug calving boss: ',my_cboss
    ImBoss = MyPE == my_cboss

    IF(ImBoss) THEN
      ALLOCATE(Prnode_count(ncalv_parts))
      Prnode_count = 0
    END IF
  ELSE
    ! only used if cboss has non calving nodes
    nonCalvBoss = MINLOC(pNCalvNodes(1,:), 1, MASK=pNCalvNodes(1,:)==MINVAL(pNCalvNodes(1,:)))-1
    ImBoss = .FALSE.
  END IF

  !croot location is i when cboss = groupmember(i)
  DO i=1, ncalv_parts
    IF(my_cboss /= cgroup_membs(i)) CYCLE
    croot = i-1
  END DO

  !Redistribute the mesh (i.e. partitioning) so that cboss partitions
  !contain entire calving/remeshing regions.
  IF(.NOT. ASSOCIATED(Mesh % Repartition)) THEN
    ALLOCATE(Mesh % Repartition(NBulk+NBdry), STAT=ierr)
    IF(ierr /= 0) PRINT *,ParEnv % MyPE,' could not allocate Mesh % Repartition'
  END IF

  !TODO send non calving nodes on ImBoss to nocalvboss
  Mesh % Repartition = ParEnv % MyPE + 1
  DO i=1,NBulk+NBdry
    IF(elem_send(i)) Mesh % Repartition(i) = my_cboss+1
    IF(ImBoss .AND. .NOT. elem_send(i)) THEN
      WRITE(Message, '(A,i0,A)') "ImBoss sending Element ",i," to NonCalvBoss"
      CALL WARN(SolverName, Message)
      Mesh % Repartition(i) = NonCalvBoss+1
    END IF
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
         MPI_INTEGER, croot, COMM_CALVE,ierr)

    IF(ImBoss) THEN

      IF(Debug) PRINT *,'boss debug prnode_count: ', Prnode_count
      GLNode_max = MAXVAL(GatheredMesh % ParallelInfo % GlobalDOFs)
      GLNode_min = MINVAL(GatheredMesh % ParallelInfo % GlobalDOFs)
      GNBulk = GatheredMesh % NumberOfBulkElements
      GNBdry = GatheredMesh % NumberOfBoundaryElements
      GNNode = GatheredMesh % NumberOfNodes

      ALLOCATE(PGDOFs_send(SUM(Prnode_count)), &
           Ptest_dist(SUM(Prnode_count)), &
           disps(ncalv_parts),GtoLNN(GLNode_min:GLNode_max),&
           gtest_dist(GNNode),&
           fixed_node(GNNode),&
           fixed_elem(GNBulk + GNBdry))

      IF(CalvingOccurs) ALLOCATE(Ptest_lset(SUM(Prnode_count)), &
                                gtest_lset(GNNode))
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

    IF (CalvingOccurs) THEN
      !Gather the level set function
      CALL MPI_GatherV(PACK(test_lset,remeshed_node), COUNT(remeshed_node), MPI_DOUBLE_PRECISION, &
            Ptest_lset, Prnode_count, disps, MPI_DOUBLE_PRECISION, croot, COMM_CALVE,ierr)
      IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")
    END IF
    !Gather the distance to front, but let it output to Ptest_lset to avoid repetitive code
    CALL MPI_GatherV(PACK(test_dist,remeshed_node), COUNT(remeshed_node), MPI_DOUBLE_PRECISION, &
          Ptest_dist, Prnode_count, disps, MPI_DOUBLE_PRECISION, croot, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")
    !END IF
    !Gather the GDOFs
    CALL MPI_GatherV(PACK(Mesh % ParallelInfo % GlobalDOFs,remeshed_node), &
         COUNT(remeshed_node), MPI_INTEGER, PGDOFs_send, Prnode_count, &
         disps, MPI_INTEGER, croot, COMM_CALVE,ierr)
    IF(ierr /= MPI_SUCCESS) CALL Fatal(SolverName,"MPI Error!")

    !Nodes with levelset value greater than remesh_thresh are kept fixed
    !as are any elements with any fixed_nodes
    !This implies the levelset is a true distance function, everywhere in the domain.
    IF(ImBoss) THEN
      DO i=1,SUM(Prnode_count)
        k = GtoLNN(PGDOFs_send(i))
        IF(k==0) CALL Fatal(SolverName, "Programming error in Global to Local NNum map")
        IF(CalvingOccurs) Gtest_lset(k) = Ptest_lset(i)
        Gtest_dist(k) = Ptest_dist(i)
      END DO

      ! calving front boundary nodes may have lset value greater than remesh dist
      ! so use dist so front boundary nodes aren't fixed
      IF(CalvingOccurs) THEN
        DO i=1,GNNode
          newdist = MINVAL((/Gtest_lset(i), Gtest_dist(i)/))
          fixed_node(i) = newdist > remesh_thresh
        END DO
      ELSE
        fixed_node = Gtest_dist > remesh_thresh
      END IF

      DO i=1, GNBulk + GNBdry
        Element => GatheredMesh % Elements(i)
        IF(ANY(fixed_node(Element % NodeIndexes))) THEN
          fixed_elem(i) = .TRUE.
        END IF
      END DO
    END IF
  END IF ! My calving front > 1

  !Nominated partition does the remeshing
  IF(ImBoss) THEN
    IF (CalvingOccurs) THEN
      !Initialise MMG datastructures
      mmgMesh = 0
      mmgLs  = 0
      mmgMet  = 0

      CALL MeshVolume(GatheredMesh, .FALSE., PreCalveVolume)

      CALL MMG3D_Init_mesh(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppLs,mmgLs, &
          MMG5_ARG_end)

      CALL GetCalvingEdgeNodes(GatheredMesh, .FALSE., EdgePairs, PairCount)
      !  now Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)
      CALL Set_MMG3D_Mesh(GatheredMesh, .TRUE., EdgePairs, PairCount)

      !Request isosurface discretization
      CALL MMG3D_Set_iparameter(mmgMesh, mmgLs, MMGPARAM_iso, 1,ierr)

      !set angle detection on (1, default) and set threshold angle (85 degrees)
      !TODO - here & in MeshRemesh, need a better strategy than automatic detection
      !i.e. feed in edge elements.

      !I think these are defunct as they are called in RemeshMMG
      ! this is unless cutting lset and anisotrophic remeshing take place in one step 
      !print *, 'first call of angle detection $$$ - turned on '
      CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgLs,MMGPARAM_angle, &
           0,ierr)

      !!! angle detection changes calving in next time steps dramatically 
      !! if on automatically turns angle on
      !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_angleDetection,&
      !     85.0_dp,ierr)

      !Turn on debug (1)
      CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgLs,MMGPARAM_debug, &
           1,ierr)

      !Turn off automatic aniso remeshing (0)
      CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgLs,MMGPARAM_aniso, &
           0,ierr)

      !Set geometric parameters for remeshing
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgLs,MMGPARAM_hmin,&
           hmin,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgLs,MMGPARAM_hmax,&
           hmax,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgLs,MMGPARAM_hausd,&
           hausd,ierr)
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgLs,MMGPARAM_hgrad,&
           hgrad,ierr)

      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgLs,MMGPARAM_rmc,&
           0.0001_dp,ierr)

      !Feed in the level set distance
      CALL MMG3D_SET_SOLSIZE(mmgMesh, mmgLs, MMG5_Vertex, GNNode ,MMG5_Scalar, ierr)
      DO i=1,GNNode
        CALL MMG3D_Set_scalarSol(mmgLs,&
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
      CALL MMG3D_Chk_meshData(mmgMesh,mmgLs,ierr)
      IF ( ierr /= 1 ) CALL EXIT(105)

      CALL MMG3D_SaveMesh(mmgMesh,"test_prels.mesh",LEN(TRIM("test_prels.mesh")),ierr)
      CALL MMG3D_SaveSol(mmgMesh, mmgLs,"test_prels.sol",LEN(TRIM("test_prels.sol")),ierr)
      !> ------------------------------ STEP  II --------------------------
      !! remesh function
      ! mmg5.5 not using isosurface discretization. More robust to remesh seperately
      ! addtionally computationally lighter as iceberg are not finely remeshed
      CALL MMG3D_mmg3dls(mmgMesh,mmgLs,0_dp,ierr)

      IF ( ierr == MMG5_STRONGFAILURE ) THEN
        PRINT*,"BAD ENDING OF MMG3DLS: UNABLE TO SAVE MESH", Time
      ELSE IF ( ierr == MMG5_LOWFAILURE ) THEN
        PRINT*,"BAD ENDING OF MMG3DLS", time
      ENDIF

      CALL MMG3D_SaveMesh(mmgMesh,"test_ls.mesh",LEN(TRIM("test_ls.mesh")),ierr)
      CALL MMG3D_SaveSol(mmgMesh, mmgLs,"test_ls.sol",LEN(TRIM("test_ls.sol")),ierr)

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
                 GatheredMesh % ParallelInfo % GInterface(j)) THEN
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

      !! Release mmg mesh
      CALL MMG3D_Free_all(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppLs,mmgLs, &
          MMG5_ARG_end)

      !MMG3DLS returns constraint = 10 on newly formed boundary elements
      !(i.e. the new calving front). Here it is set to front_BC_id
      !And set all BC BodyIDs based on constraint
      DO i=NBulk+1, NBulk + NBdry

        IF(RmElem(i)) CYCLE

        Element => NewMeshR % Elements(i)
        geom_id => Element % BoundaryInfo % Constraint

        IF(geom_id == 10) THEN
          geom_id = front_BC_id

          Element % BodyId  = front_body_id
        END IF

        CALL Assert((geom_id > 0) .AND. (geom_id <= Model % NumberOfBCs),&
             SolverName,"Unexpected BC element body id!")

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

      NNodes = NewMeshR % NumberOfNodes
      NBulk = NewMeshR % NumberOfBulkElements
      NBdry = NewMeshR % NumberOfBoundaryElements

      !sometimes some icebergs are not removed fully from the MMG levelset
      !find all connected nodes in groups and remove in groups not in the
      !main icebody
      !limit here of 10 possible mesh 'islands'
      ALLOCATE(FoundNode(10, NNodes), ElNodes(4), &
              UsedElem(NBulk), NewNodes(NBulk))
      FoundNode = .FALSE.! count of nodes found
      UsedElem = .FALSE. !count of elems used
      NewNodes=.FALSE. !count of new elems in last sweep
      island=0 ! count of different mesh islands
      HasNeighbour=.FALSE. ! whether node has neighour
      DO WHILE(COUNT(FoundNode) < NNodes)
        ! if last sweep has found no new neighbour nodes then move onto new island
        IF(.NOT. HasNeighbour) THEN
          island = island + 1
          IF(Island > 10) CALL FATAL(SolverName, 'More than 10 icebergs left after levelset')
          !find first unfound node
          Inner: DO i=1, NNodes
            IF(ANY(FoundNode(:, i))) CYCLE
            NewNodes(i) = .TRUE.
            FoundNode(island, i) = .TRUE.
            EXIT Inner
          END DO Inner
        END IF
        HasNeighbour=.FALSE.
        nodes=PACK((/ (i,i=1,SIZE(NewNodes)) /), NewNodes .eqv. .TRUE.)
        NewNodes=.FALSE.
        DO i=1, Nbulk
          IF(UsedElem(i)) CYCLE !already used elem
          DO j=1, SIZE(Nodes)
            node=Nodes(j)
            Element => NewMeshR % Elements(i)
            ElNodes = Element % NodeIndexes
            IF(.NOT. ANY(ElNodes==node)) CYCLE ! doesn't contain node
            UsedElem(i) = .TRUE. ! got nodes from elem
            DO k=1,SIZE(ElNodes)
              idx = Element % NodeIndexes(k)
              IF(idx == Node) CYCLE ! this nodes
              IF(FoundNode(island,idx)) CYCLE ! already got
              FoundNode(island, idx) = .TRUE.
              counter = counter + 1
              HasNeighbour = .TRUE.
              NewNodes(idx) = .TRUE.
            END DO
          END DO
        END DO
      END DO

      ALLOCATE(IslandCounts(10))
      DO i=1,10
        IslandCounts(i) = COUNT(FoundNode(i, :))
      END DO
      Island = COUNT(IslandCounts > 0)
      DEALLOCATE(IslandCounts) ! reused

      ! if iceberg not removed mark nodes and elems
      IF(Island > 1) THEN
        CALL WARN(SolverName, 'Mesh island found after levelset- removing...')
        ALLOCATE(RmIslandNode(NNodes), RmIslandElem(Nbulk+NBdry),&
                IslandCounts(Island))
        RmIslandNode = .FALSE.
        RmIslandElem = .FALSE.
        counter=0
        DO i=1, 10
          IF(COUNT(FoundNode(i, :)) == 0) CYCLE
          counter=counter+1
          IslandCounts(counter) = COUNT(FoundNode(i, :))
        END DO
        DO i=1, Island
          IF(IslandCounts(i) < MAXVAL(IslandCounts)) THEN
            DO j=1, NNodes
              IF(.NOT. FoundNode(i,j)) CYCLE! if not part of island
              RmIslandNode(j) = .TRUE.
            END DO
          END IF
        END DO
        ! mark all elems with rm nodes as need removing
        nodes = PACK((/ (i,i=1,SIZE(RmIslandNode)) /), RmIslandNode .eqv. .TRUE.)
        DO i=1, Nbulk+Nbdry
          Element => NewMeshR % Elements(i)
          ElNodes = Element % NodeIndexes
          !any([(any(A(i) == B), i = 1, size(A))])
          IF( .NOT. ANY([(ANY(ElNodes(i)==Nodes), i = 1, SIZE(ElNodes))])) CYCLE
          RmIslandElem(i) = .TRUE.
        END DO

        !Chop out missed iceberg if they exist
        CALL CutMesh(NewMeshR, RmIslandNode, RmIslandElem)

        !modify rmelem and rmnode to include island removal
        counter=0
        DO i=1, SIZE(RmElem)
          IF(RmElem(i)) CYCLE
          counter=counter+1
          IF(RmIslandElem(Counter)) RmElem(i) =.TRUE.
        END DO
        counter=0
        DO i=1, SIZE(RmNode)
          IF(RmNode(i)) CYCLE
          counter=counter+1
          IF(RmIslandNode(Counter)) RmNode(i) =.TRUE.
        END DO
        CALL INFO(SolverName, 'Mesh island removed', level=10)
      END IF

    END IF ! CalvingOccurs

      DoAniso = .TRUE.
      IF(DoAniso .AND. CalvingOccurs) THEN

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

        !Parameters for anisotropic remeshing are set in the Materials section, or can &
        !optionally be passed as a valuelist here. They have prefix RemeshMMG3D
        !TODO - apparently there is beta-testing capability to do both levelset cut and anisotropic
        !remeshing in the same step:
        !https://forum.mmgtools.org/t/level-set-and-anisotropic-mesh-optimization/369/3
        ! GetCalvingEdgeNodes detects all shared boundary edges, to keep them sharp

        ! calculate calved volume
        CALL MeshVolume(NewMeshR, .FALSE., PostCalveVolume)

        CalveVolume = PreCalveVolume - PostCalveVolume

        PRINT*, 'CalveVolume', CalveVolume
      END IF
    END IF

    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    Distributed = .TRUE.
    IF(DoAniso .AND. .NOT. Distributed) THEN
      ! remeshing but no calving
      IF(ImBoss .AND. .NOT. CalvingOccurs) NewMeshR => GatheredMesh

      !remeshing for calvign and no calving
      IF(ImBoss) CALL GetCalvingEdgeNodes(NewMeshR, .FALSE., REdgePairs, RPairCount)
            !  now Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)
      MeshParams => GetMaterial(Model % Mesh % Elements(1))
      CALL SequentialRemeshParMMG(Model, NewMeshR, NewMeshRR, ImBoss, REdgePairs, &
                RPairCount,new_fixed_node,new_fixed_elem, MeshParams)
      IF(ImBoss) THEN
            CALL ReleaseMesh(NewMeshR)
            NewMeshR => NewMeshRR
            NewMeshRR => NULL()

          !Update parallel info from old mesh nodes (shared node neighbours)
          CALL MapNewParallelInfo(GatheredMesh, NewMeshR)
          CALL ReleaseMesh(GatheredMesh)
          ! CALL ReleaseMesh(NewMeshR)
          GatheredMesh => NewMeshR
          NewMeshR => NULL()
          NewMeshRR => NULL()
      END IF
    END IF

    IF(DoAniso .AND. Distributed) THEN
      ! remeshing but no calving
      IF(ImBoss .AND. .NOT. CalvingOccurs) THEN
        NewMeshR => GatheredMesh
      !ELSE IF (ImBoss .AND. CalvingOccurs) THEN
        !CALL MapNewParallelInfo(GatheredMesh, NewMeshR)
        !DO i=1, NewMeshR % NumberOfNodes
          !PRINT*, 'boss', i, NewMeshR % ParallelInfo % GlobalDOFs(i)
          !NewMeshR % ParallelInfo % GlobalDOFs
        !END DO
      END IF

      IF(.NOT. ImBoss) THEN
        NewMeshR => AllocateMesh( 0, 0, 0, InitParallel = .TRUE.)
        ALLOCATE( NewMeshR % ParallelInfo % GlobalDofs( 0 ))
        !NewMeshR % ParallelInfo % GlobalDofs = 0
      END IF

      !IF(.NOT. ImBoss) PRINT*, ParEnv % MyPE, 'meshcheck', NewMeshR % ParallelInfo % GlobalDOFs
      !Now each partition has GatheredMesh, we need to renegotiate globalDOFs
      !CALL RenumberGDOFs(Mesh, NewMeshR)

      !and global element numbers
      !CALL RenumberGElems(NewMeshR)


      ! share fixed nodes
      ! assumption here is globaldofs == nodenumber
      !

      IF (ImBoss) THEN
        RNNodes = NewMeshR % NumberOfNodes
        RNElems = NewMeshR % NumberOfBulkElements + NewMeshR % NumberOfBoundaryElements

        DO i=1, RNNodes
          NewMeshR % ParallelInfo % GlobalDOFs(i) = i
        END DO
        DO i=1, RNElems
          NewMeshR % Elements(i) % GElementIndex = i
        END DO
      END IF
      CALL MPI_BCAST(RNNodes, 1, MPI_INTEGER, my_cboss, ELMER_COMM_WORLD, ierr)
      CALL MPI_BCAST(RNElems, 1, MPI_INTEGER, my_cboss, ELMER_COMM_WORLD, ierr)
      IF(.NOT. ImBoss) ALLOCATE(new_fixed_node(RNNodes), new_fixed_elem(RNElems))
      CALL MPI_BCAST(new_fixed_node,RNNodes,MPI_LOGICAL, my_cboss, ELMER_COMM_WORLD, ierr)
      CALL MPI_BCAST(new_fixed_elem,RNElems, MPI_LOGICAL, my_cboss, ELMER_COMM_WORLD, ierr)

      CALL Zoltan_Interface( Model, NewMeshR, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )

      DistributedMesh => RedistributeMesh(Model, NewMeshR, .TRUE., .FALSE.)
      CALL ReleaseMesh(NewMeshR)

      !remeshing for calvign and no calving
      !not sure if this works in parallel
      CALL GetCalvingEdgeNodes(DistributedMesh, .FALSE., REdgePairs, RPairCount)
      !  now Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)
      MeshParams => GetMaterial(Model % Mesh % Elements(1))
      CALL DistributedRemeshParMMG(Model, DistributedMesh, NewMeshRR, REdgePairs, &
                RPairCount,new_fixed_node,new_fixed_elem, MeshParams)
      CALL ReleaseMesh(NewMeshR)
      NewMeshR => NewMeshRR
      NewMeshRR => NULL()

      !Update parallel info from old mesh nodes (shared node neighbours)
      CALL MapNewParallelInfo(GatheredMesh, NewMeshR)
      CALL ReleaseMesh(GatheredMesh)
      ! CALL ReleaseMesh(NewMeshR)
      GatheredMesh => NewMeshR
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
     IF(GatheredMesh % ParallelInfo % GInterface(i)) THEN
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

   ParEnv % IsNeighbour(:)  = .FALSE.
   DO i=1, Mesh % NumberOfNodes
     IF ( ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours) ) THEN
       DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
         proc = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)
         IF ( ParEnv % Active(proc+1).AND.proc/=ParEnv % MYpe) &
             ParEnv % IsNeighbour(proc+1) = .TRUE.
       END DO
     END IF
   END DO

   !Call zoltan to determine redistribution of mesh
   ! then do the redistribution
   !-------------------------------

   CALL Zoltan_Interface( Model, GatheredMesh, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )

   FinalMesh => RedistributeMesh(Model, GatheredMesh, .TRUE., .FALSE.)

   CALL CheckMeshQuality(FinalMesh)

   FinalMesh % Name = TRIM(Mesh % Name)
   FinalMesh % OutputActive = .TRUE.
   FinalMesh % Changed = .TRUE. 

   !Actually switch the model's mesh
   CALL SwitchMesh(Model, Solver, Mesh, FinalMesh)

   !After SwitchMesh because we need GroundedMask
   CALL EnforceGroundedMask(Model, Mesh)

   !Calculate mesh volume
   CALL MeshVolume(Model % Mesh, .TRUE., PostCalveVolume)
   IF(MyPe == 0) PRINT*, 'Post calve volume: ', PostCalveVolume, ' after timestep', Time

   !Recompute mesh bubbles etc
   CALL MeshStabParams( Model % Mesh )

   !Release the old mesh
   CALL ReleaseMesh(GatheredMesh)

   !remove mesh update
   CALL ResetMeshUpdate(Model, Solver)

END SUBROUTINE CalvingRemeshParMMG

! test routine for parallel remeshing of a distributed mesh
SUBROUTINE ParMMGRemeshing ( Model, Solver, dt, TransientSimulation )

  USE MeshUtils
  USE CalvingGeometry
  USE MeshPartition
  USE MeshRemeshing

  IMPLICIT NONE

!-----------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!-----------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, OutMesh, FinalMesh
  TYPE(ValueList_t), POINTER :: SolverParams
  REAL(kind=dp), POINTER :: WorkReal(:)
  INTEGER, POINTER :: WorkPerm(:)
  INTEGER, ALLOCATABLE :: EdgePairs(:,:)
  LOGICAL :: Parallel
  INTEGER :: i,j,k,n, PairCount
  !------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = "ParMMGRemeshing"

#ifdef HAVE_PARMMG

  IF(ParEnv % PEs > 1) THEN
      Parallel = .TRUE.
  ELSE
      CALL FATAL(Solvername, 'Needs to be run in parallel')
  END IF

  Mesh => Model % Mesh
  SolverParams => GetSolverParams()

  n = Mesh % NumberOfNodes
  ALLOCATE(WorkPerm(n), WorkReal(n))
  WorkReal = 0.0_dp
  WorkPerm = [(i,i=1,n)]
  CALL VariableAdd(Mesh % Variables, Mesh, Solver, "dummy", &
        1, WorkReal, WorkPerm)
  NULLIFY(WorkReal,WorkPerm)

  CALL GetCalvingEdgeNodes(Mesh, .FALSE., EdgePairs, PairCount)

  CALL DistributedRemeshParMMG(Model, Mesh,OutMesh, EdgePairs, PairCount,Params=SolverParams)

  !Now each partition has GatheredMesh, we need to renegotiate globalDOFs
  CALL RenumberGDOFs(Mesh, OutMesh)

  !and global element numbers
  CALL RenumberGElems(OutMesh)

  CALL Zoltan_Interface( Model, OutMesh, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )

  FinalMesh => RedistributeMesh(Model, OutMesh, .TRUE., .FALSE.)

  FinalMesh % Name = TRIM(Mesh % Name)
  FinalMesh % OutputActive = .TRUE.
  FinalMesh % Changed = .TRUE.

  !Actually switch the model's mesh
  ! also releases old mesh
  CALL SwitchMesh(Model, Solver, Mesh, FinalMesh)

  CALL ReleaseMesh(OutMesh)
  !CALL ReleaseMesh(Mesh)

#else
  CALL FATAL(SolverName, 'ParMMG needs to be installed to use this solver!')
#endif

END SUBROUTINE ParMMGRemeshing