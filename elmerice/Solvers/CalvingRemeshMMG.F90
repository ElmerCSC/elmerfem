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
! set is negative inside a calving event, and positive in the remaining domain.

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

#include "mmg/common/libmmgtypesf.h"
#ifndef MMGVERSION_H
#define MMG_VERSION_LT(MAJOR,MINOR) 1
#endif

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !--------------------------------------
  TYPE(Variable_t), POINTER :: CalvingVar,DistanceVar,mmgVar
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh,GatheredMesh,NewMeshR,NewMeshRR,FinalMesh,ParMetisMesh
  TYPE(Element_t),POINTER :: Element, ParentElem
  INTEGER :: i,j,k,NNodes,GNBulk, GNBdry, GNNode, NBulk, Nbdry, ierr, &
       my_cboss,MyPE, PEs,CCount, counter, GlNode_min, GlNode_max,adjList(4),&
       front_BC_ID, front_body_id, my_calv_front,calv_front, ncalv_parts, &
       group_calve, comm_calve, group_world,ecode, NElNodes, target_bodyid,gdofs(4), &
       PairCount,RPairCount, NCalvNodes, croot, nonCalvBoss,&
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges
  INTEGER, POINTER :: NodeIndexes(:), geom_id
  INTEGER, ALLOCATABLE :: Prnode_count(:), cgroup_membs(:),disps(:), &
       PGDOFs_send(:),pcalv_front(:),GtoLNN(:),EdgePairs(:,:),REdgePairs(:,:), ElNodes(:),&
       Nodes(:), IslandCounts(:), pNCalvNodes(:,:), TetraQuality(:)
  REAL(KIND=dp) :: test_thresh, test_point(3), remesh_thresh, hmin, hmax, hgrad, hausd, &
       newdist, Quality, PauseVolumeThresh, MaxBergVolume, RmcValue
  REAL(KIND=dp), ALLOCATABLE :: test_dist(:), test_lset(:), Ptest_lset(:), Gtest_lset(:), &
       target_length(:,:), Ptest_dist(:), Gtest_dist(:), hminarray(:), hausdarray(:)
  REAL(KIND=dp), POINTER :: WorkArray(:,:) => NULL()
  LOGICAL, ALLOCATABLE :: calved_node(:), remeshed_node(:), fixed_node(:), fixed_elem(:), &
       elem_send(:), RmElem(:), RmNode(:),new_fixed_node(:), new_fixed_elem(:), FoundNode(:,:), &
       UsedElem(:), NewNodes(:), RmIslandNode(:), RmIslandElem(:), PartGotNodes(:)
  LOGICAL :: ImBoss, Found, Isolated, Debug,DoAniso,NSFail,CalvingOccurs,&
       RemeshOccurs,CheckFlowConvergence, NoNewNodes, RSuccess, Success,&
       SaveMMGMeshes, SaveMMGSols, PauseSolvers, PauseAfterCalving, FixNodesOnRails, &
       SolversPaused, NewIceberg, GotNodes(4), CalvingFileCreated=.FALSE., SuppressCalv,&
       DistributedMesh,SaveTerminus
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, CalvingVarName, MeshName, SolName, &
       premmgls_meshfile, mmgls_meshfile, premmgls_solfile, mmgls_solfile,&
       RepartMethod, Filename
  TYPE(Variable_t), POINTER :: TimeVar
  INTEGER :: Time, remeshtimestep, proc, idx, island, node, MaxLSetIter, mmgloops
  REAL(KIND=dp) :: TimeReal, PreCalveVolume, PostCalveVolume, CalveVolume, LsetMinQuality

  SAVE :: WorkArray, CalvingFileCreated

  SolverParams => GetSolverParams()
  SolverName = "CalvingRemeshMMG"
  Debug=.FALSE.
  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes
  NBulk = Mesh % NumberOfBulkElements
  NBdry = Mesh % NumberOfBoundaryElements

  calv_front = 1
  MyPe = ParEnv % MyPE
  PEs = ParEnv % PEs

#if MMG_VERSION_LT(5,6)
  PRINT*, SolverName, ': Starting MMG'
#elif MMG_VERSION_LT(5,8)
  PRINT*, SolverName, ': Starting MMG'
#else
  CALL FATAL(SolverName, 'Calving code only works with MMG 5.5')
#endif

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  TimeReal = TimeVar % Values(1)
  Time = INT(TimeReal)

  ! create mmg variable as this is added to one proc in remeshing changes following merge Nov/23
  mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop', ThisOnly = .TRUE.)
  IF(.NOT. ASSOCIATED(mmgVar) ) THEN
    CALL VariableAddVector( Model % Mesh % Variables,Model % Mesh,&
        Name='MMG Loop',Global=.TRUE.)
    mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop' )
  END IF

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
  !For testing - set a calving levelset function to mimic calving events
  !-------------------

  !TODO - unhardcode (detect?) this
  front_BC_id = 1
  front_body_id =  ListGetInteger( &
       Model % BCs(front_bc_id) % Values, 'Body Id', Found, 1, Model % NumberOfBodies )

  WorkArray => ListGetConstRealArray(SolverParams, "Mesh Hmin", Found)
  IF(.NOT. Found) CALL FATAL(SolverName, 'Provide hmin input array to be iterated through: "Mesh Hmin"')
  MaxLsetIter= SIZE(WorkArray(:,1))
  ALLOCATE(hminarray(MaxLsetIter))
  hminarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  hmax = ListGetConstReal(SolverParams, "Mesh Hmax",  DefValue=4000.0_dp)
  hgrad = ListGetConstReal(SolverParams,"Mesh Hgrad", DefValue=0.5_dp)
  WorkArray => ListGetConstRealArray(SolverParams, "Mesh Hausd", Found)
  IF(.NOT. Found) CALL FATAL(SolverName, 'Provide hausd input array to be iterated through: "Mesh Hausd"')
  IF(MaxLsetIter /= SIZE(WorkArray(:,1))) CALL FATAL(SolverName, 'The number of hmin options &
          must equal the number of hausd options')
  ALLOCATE(hausdarray(MaxLsetIter))
  hausdarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  remesh_thresh = ListGetConstReal(SolverParams,"Remeshing Distance", DefValue=1000.0_dp)
  LsetMinQuality = ListGetConstReal(SolverParams,"Mesh Min Quality", DefValue=0.00001_dp)
  RmcValue = ListGetConstReal(SolverParams,"Mesh Rmc Value", DefValue=1e-15_dp)
  CalvingVarName = ListGetString(SolverParams,"Calving Variable Name", DefValue="Calving Lset")
  SaveMMGMeshes = ListGetLogical(SolverParams,"Save MMGLS Meshes", DefValue=.FALSE.)
  SaveMMGSols = ListGetLogical(SolverParams,"Save MMGLS Sols", DefValue=.FALSE.)
  IF(SaveMMGMeshes) THEN
    premmgls_meshfile = ListGetString(SolverParams, "Pre MMGLS Mesh Name", UnfoundFatal = .TRUE.)
    mmgls_meshfile = ListGetString(SolverParams, "MMGLS Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF
  IF(SaveMMGSols) THEN
    premmgls_solfile = ListGetString(SolverParams, "Pre MMGLS Sol Name", UnfoundFatal = .TRUE.)
    mmgls_solfile = ListGetString(SolverParams, "MMGLS Output Sol Name", UnfoundFatal = .TRUE.)
  END IF
  PauseAfterCalving = ListGetLogical(SolverParams, "Pause After Calving Event", Found)
  IF(.NOT. Found) THEN
     CALL Info(SolverName, "Can't find 'Pause After Calving Event' logical in Solver section, &
          & assuming True")
     PauseAfterCalving = .TRUE.
  END IF
  FixNodesOnRails = ListGetLogical(SolverParams,"Fix Nodes On Rails", DefValue=.TRUE.)
  SuppressCalv = ListGetLogical(SolverParams,"Suppress Calving", DefValue=.FALSE.)
  SaveTerminus = ListGetLogical(SolverParams,"Save Terminus", DefValue=.TRUE.)

  IF(ParEnv % MyPE == 0) THEN
    PRINT *,ParEnv % MyPE,' hmin: ',hminarray
    PRINT *,ParEnv % MyPE,' hmax: ',hmax
    PRINT *,ParEnv % MyPE,' hgrad: ',hgrad
    PRINT *,ParEnv % MyPE,' hausd: ',hausdarray
    PRINT *,ParEnv % MyPE,' remeshing distance: ',remesh_thresh
    PRINT *,ParEnv % MyPE,' Max Levelset Iterations ', MaxLsetIter
  END IF

  !If FlowSolver failed to converge (usually a result of weird mesh), large unphysical
  !calving events can be predicted. So, turn off CalvingOccurs, and ensure a remesh
  !Also undo this iterations mesh update
  NSFail = ListGetLogical(Model % Simulation, "Flow Solution Failed",CheckFlowConvergence)
  CalvingOccurs = ListGetLogical( Model % Simulation, 'CalvingOccurs', Found )
  RemeshOccurs=.TRUE.
  IF(CheckFlowConvergence) THEN
    IF(NSFail) THEN
      CalvingOccurs = .FALSE.
      RemeshOccurs = .TRUE.
      CALL Info(SolverName, "Remeshing but not calving because NS failed to converge.")
    END IF
  END IF
  IF(SuppressCalv) CalvingOccurs = .FALSE.

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
    my_cboss = MINLOC(pNCalvNodes(1,:), 1, MASK=pNCalvNodes(1,:)==MAXVAL(pNCalvNodes(1,:)))-1
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
    my_cboss = MINLOC(pNCalvNodes(1,:), 1, MASK=pNCalvNodes(1,:)==MAXVAL(pNCalvNodes(1,:)))-1
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

  IF(ASSOCIATED(Mesh % Repartition)) THEN
    DEALLOCATE(Mesh % Repartition)
    Mesh % Repartition => NULL()
  END IF

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
      IF(CalvingOccurs) gtest_lset = remesh_thresh + 500.0 !Ensure any far (unshared) nodes are fixed

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

      CALL MeshVolume(GatheredMesh, .FALSE., PreCalveVolume)

      CALL GetCalvingEdgeNodes(GatheredMesh, .FALSE., EdgePairs, PairCount)

      mmgloops = 0

10 CONTINUE

      mmgloops=mmgloops+1
      hmin = hminarray(mmgloops)
      hausd = hausdarray(mmgloops)

      WRITE(Message, '(A,F10.5,A,F10.5)') 'Applying levelset with Hmin ',Hmin, ' and Hausd ', Hausd
      CALL INFO(SolverName, Message)
      !Initialise MMG datastructures
      mmgMesh = 0
      mmgLs  = 0
      mmgMet  = 0

      CALL MMG3D_Init_mesh(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppLs,mmgLs, &
          MMG5_ARG_end)

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
           RmcValue,ierr)

      CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMGPARAM_nosizreq,&
           1,ierr)

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
      IF(SaveMMGMeshes) THEN
        WRITE(MeshName, '(A,i0,A)') TRIM(premmgls_meshfile), time, '.mesh'
        CALL MMG3D_SaveMesh(mmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
      END IF
      IF(SaveMMGSols) THEN
        WRITE(SolName, '(A,i0,A)') TRIM(premmgls_solfile), time, '.sol'
        CALL MMG3D_SaveSol(mmgMesh, mmgLs,SolName,LEN(TRIM(SolName)),ierr)
      END IF
      !> ------------------------------ STEP  II --------------------------
      !! remesh function
      ! mmg5.5 not using isosurface discretization. More robust to remesh seperately
      ! addtionally computationally lighter as iceberg are not finely remeshed
      CALL MMG3D_mmg3dls(mmgMesh,mmgLs,0_dp,ierr)

      IF ( ierr == MMG5_STRONGFAILURE ) THEN
        PRINT*,"BAD ENDING OF MMG3DLS: UNABLE TO SAVE MESH", Time
        !! Release mmg mesh
        CALL MMG3D_Free_all(MMG5_ARG_start, &
            MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
            MMG5_ARG_end)
        WRITE(Message, '(A,F10.5,A,F10.5)') 'Levelset failed with Hmin ',Hmin, ' and Hausd ', Hausd
        CALL WARN(SolverName, Message)
        IF(mmgloops==MaxLsetIter) THEN
          Success=.FALSE.
          CalvingOccurs = .FALSE.
          GO TO 20
        ELSE
          GO TO 10
        END IF
      ELSE IF ( ierr == MMG5_LOWFAILURE ) THEN
        PRINT*,"LOW ENDING OF MMG3DLS", time
      ENDIF
      IF(SaveMMGMeshes) THEN
        WRITE(MeshName, '(A,i0,A)') TRIM(mmgls_meshfile), time, '.mesh'
        CALL MMG3D_SaveMesh(mmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
      END IF
      IF(SaveMMGSols) THEN
        WRITE(SolName, '(A,i0,A)') TRIM(mmgls_solfile), time, '.sol'
        CALL MMG3D_SaveSol(mmgMesh, mmgLs,SolName,LEN(TRIM(SolName)),ierr)
      END IF

      CALL MMG3D_Get_meshSize(mmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

      counter=0
      Do i=1, NTetras
        CALL MMG3D_Get_TetrahedronQuality(mmgMesh, mmgSol, i, Quality)
        IF(Quality == 0) CALL WARN(SolverName, 'Levelset could not determine elem quality')
        IF(Quality <= LsetMinQuality) counter = counter+1
      END DO
      IF ( Counter > 0 ) THEN
        PRINT*,"Bad elements detected - reruning levelset"
        !! Release mmg mesh
        CALL MMG3D_Free_all(MMG5_ARG_start, &
            MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
            MMG5_ARG_end)
          WRITE(Message, '(A,F10.5,A,F10.5)') 'Levelset failed with Hmin ',Hmin, ' and Hausd ', Hausd
        CALL WARN(SolverName, Message)
        IF(mmgloops==MaxLsetIter) THEN
          Success=.FALSE.
          CalvingOccurs = .FALSE.
          GO TO 20
        ELSE
          GO TO 10
        END IF
        !STOP MMG5_STRONGFAILURE
      ENDIF

      CALL Get_MMG3D_Mesh(NewMeshR, .TRUE., new_fixed_node, new_fixed_elem, Calving=.TRUE.)

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

      !output calving stats to file and find maxbergvolume
      CALL CalvingStatsMMG(SolverParams, NewMeshR, RmNode, RmElem, &
          CalvingFileCreated, MaxBergVolume)

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
              UsedElem(NBulk))
      FoundNode = .FALSE.! count of nodes found
      UsedElem = .FALSE. !count of elems used
      island=0 ! count of different mesh islands
      NoNewNodes=.TRUE. ! whether node has neighour
      DO WHILE(COUNT(FoundNode) < NNodes)
        IF(NoNewNodes) THEN
          island = island + 1
          NewIceberg = .TRUE.
        END IF
        NoNewNodes = .TRUE.
        DO i=1, NBulk
          IF(UsedElem(i)) CYCLE
          Element => NewMeshR % Elements(i)
          ElNodes = Element % NodeIndexes
          ! if there are not any matching nodes and its not a new iceberg
          DO j=1, 4
            GotNodes(j) = ANY(FoundNode(:, ElNodes(j)))
          END DO
          IF(ALL(.NOT. GotNodes) .AND. .NOT. NewIceberg) CYCLE
          NewIceberg = .FALSE.
          UsedElem(i) = .TRUE.
          FoundNode(island, ElNodes) = .TRUE.
          NoNewNodes = .FALSE.
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

20    CONTINUE

      DEALLOCATE(hminarray, hausdarray)

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
        !ALLOCATE(target_length(NewMeshR % NumberOfNodes,3))
        !target_length(:,1) = 300.0
        !target_length(:,2) = 300.0
        !target_length(:,3) = 50.0

        !Parameters for anisotropic remeshing are set in the Materials section, or can &
        !optionally be passed as a valuelist here. They have prefix RemeshMMG3D
        !TODO - apparently there is beta-testing capability to do both levelset cut and anisotropic
        !remeshing in the same step:
        !https://forum.mmgtools.org/t/level-set-and-anisotropic-mesh-optimization/369/3
        ! GetCalvingEdgeNodes detects all shared boundary edges, to keep them sharp

        ! calculate calved volume
        CALL MeshVolume(NewMeshR, .FALSE., PostCalveVolume)

        CalveVolume = PreCalveVolume - PostCalveVolume

        PRINT*, 'CalveVolume', CalveVolume, 'MaxBergVolume', MaxBergVolume

        CALL GetCalvingEdgeNodes(NewMeshR, .FALSE., REdgePairs, RPairCount)
        !  now Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)
        CALL RemeshMMG3D(Model, NewMeshR, NewMeshRR,REdgePairs, RPairCount,&
            NodeFixed=new_fixed_node, ElemFixed=new_fixed_elem, Success=RSuccess)
        IF(.NOT. RSuccess) THEN
          CALL WARN(SolverName, 'Remeshing failed - no calving at this timestep')
          GO TO 30 ! remeshing failed so no calving at this timestep
        END IF
        CALL ReleaseMesh(NewMeshR)
        NewMeshR => NewMeshRR
        NewMeshRR => NULL()

        !Update parallel info from old mesh nodes (shared node neighbours)
        CALL MapNewParallelInfo(GatheredMesh, NewMeshR)

      ELSE IF (DoAniso) THEN
         ! remeshing but no calving
        CALL GetCalvingEdgeNodes(GatheredMesh, .FALSE., REdgePairs, RPairCount)
        !  now Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)
        CALL RemeshMMG3D(Model, GatheredMesh, NewMeshR,REdgePairs, RPairCount,&
            NodeFixed=fixed_node, ElemFixed=fixed_elem, Success=RSuccess)
                !Update parallel info from old mesh nodes (shared node neighbours)
        IF(.NOT. RSuccess) THEN
          CALL WARN(SolverName, 'Remeshing failed - no calving at this timestep')
          GO TO 30 ! remeshing failed so no calving at this timestep
        END IF
        CALL MapNewParallelInfo(GatheredMesh, NewMeshR)

      ELSE ! Not DoAniso

        !Update parallel info from old mesh nodes (shared node neighbours)
        IF (CalvingOccurs) CALL MapNewParallelInfo(GatheredMesh, NewMeshR)

      END IF

      CALL ReleaseMesh(GatheredMesh)
      ! CALL ReleaseMesh(NewMeshR)
      GatheredMesh => NewMeshR
      NewMeshR => NULL()
      NewMeshRR => NULL()
   END IF ! ImBoss

30 CONTINUE

   !Wait for all partitions to finish
   IF(My_Calv_Front>0) THEN
     CALL MPI_BARRIER(COMM_CALVE,ierr)
     CALL MPI_COMM_FREE(COMM_CALVE,ierr)
     CALL MPI_GROUP_FREE(group_world,ierr)
   END IF
   CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)

   CALL MPI_BCAST(CalvingFileCreated, 1, MPI_LOGICAL, my_cboss, ELMER_COMM_WORLD, ierr)
   CALL MPI_BCAST(RSuccess, 1, MPI_LOGICAL, my_cboss, ELMER_COMM_WORLD, ierr)
   IF(.NOT. RSuccess) THEN
      IF(ImBoss .AND. CalvingFileCreated) THEN
        Filename = ListGetString(SolverParams,"Calving Stats File Name", Found)
        IF(.NOT. Found) THEN
          CALL WARN('CalvingStat', 'Output file name not given so using CalvingStats.txt')
          Filename = "CalvingStats.txt"
        END IF
        OPEN( 36, FILE=filename, STATUS='UNKNOWN', ACCESS='APPEND')
        WRITE(36, '(A,i0)') 'Remeshing failed: ', GetTimestep()
        CLOSE(36)
      END IF

      CALL ReleaseMesh(GatheredMesh)
      ! Release GatheredMesh % Redistribution
      IF(ASSOCIATED(GatheredMesh % Repartition)) THEN
          DEALLOCATE(GatheredMesh % Repartition)
          GatheredMesh % Repartition => NULL()
      END IF
      SolversPaused = ListGetLogical(Model % Simulation, 'Calving Pause Solvers', DefValue=.FALSE.)
      !remove mesh update
      IF(.NOT. SolversPaused) THEN
        CALL ResetMeshUpdate(Model, Solver)
      ELSE ! solvers paused so mesh must have changed since last free surface
        Mesh % MeshTag = Mesh % MeshTag + 1
      END IF

      ! make sure solvers are unpaused so front advances
      CALL PauseCalvingSolvers(Model, SolverParams, .FALSE.)
      RETURN
   END IF

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

   IF(SaveTerminus) CALL SaveTerminusPosition(Model, Solver, GatheredMesh, ImBoss)

   !Call zoltan to determine redistribution of mesh
   ! then do the redistribution
   !-------------------------------

   RepartMethod = ListGetString(Model % Solver % Values,"Repartition Method", Found)
   SELECT CASE( RepartMethod )
   CASE( 'parmetis' ) ! bit of a hack here. Parmetis requires fully distributed mesh so ensure all parts have one element
      ALLOCATE(PartGotNodes(ParEnv % PEs))
      CALL MPI_ALLGATHER(GatheredMesh % NumberOfNodes > 0, 1, MPI_LOGICAL, PartGotNodes, &
          1, MPI_LOGICAL, ELMER_COMM_WORLD, ierr)

      DistributedMesh = ALL(PartGotNodes)

      IF(.NOT. ASSOCIATED(GatheredMesh % Repartition)) THEN
        ALLOCATE(GatheredMesh % Repartition(GatheredMesh % NumberOfBulkElements+&
                GatheredMesh % NumberOfBoundaryElements), STAT=ierr)
        IF(ierr /= 0) PRINT *,ParEnv % MyPE,' could not allocate Mesh % Repartition'
      END IF

      GatheredMesh % Repartition = ParEnv % MyPE + 1
      IF(ImBoss) THEN
        counter=0
        DO i=1,ParEnv % PEs
          IF(PartGotNodes(i)) CYCLE ! got nodes
          counter=counter+1
          IF(counter > GatheredMesh % NumberOfNodes) &
            CALL FATAL(SolverName, 'CalvBoss does not have enough elems to share for ParMetis repartitioning')
          GatheredMesh % Repartition(counter) = i
          DO j=GatheredMesh % NumberOfBulkElements+1, GatheredMesh % NumberOfBulkElements+&
            GatheredMesh % NumberOfBoundaryElements
            Element => GatheredMesh % Elements(j)
            ParentElem => Element % BoundaryInfo % Left
            IF(.NOT. ASSOCIATED(ParentElem)) THEN
              ParentElem => Element % BoundaryInfo % Right
            END IF
            CALL Assert(ASSOCIATED(ParentElem),SolverName,"Boundary element has no parent!")
            IF(ParentElem % ElementIndex == counter) THEN
              GatheredMesh % Repartition(j) = i
            END IF
          END DO
        END DO
      END IF

      ParMetisMesh => RedistributeMesh(Model, GatheredMesh, .TRUE., .FALSE.)

      IF(ASSOCIATED(ParMetisMesh % Repartition)) THEN
        DEALLOCATE(ParMetisMesh % Repartition)
        ParMetisMesh % Repartition => NULL()
      END IF

      CALL ReleaseMesh(GatheredMesh)
      GatheredMesh => ParMetisMesh
      ParMetisMesh => NULL()
   END SELECT

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

   CALL CheckBaseFreeSurface(Model, Mesh, 0.01_dp)

   IF(FixNodesOnRails) CALL EnforceLateralMargins(Model, Solver % Values)

   !Calculate mesh volume
   CALL MeshVolume(Model % Mesh, .TRUE., PostCalveVolume)
   IF(MyPe == 0) PRINT*, 'Post calve volume: ', PostCalveVolume, ' after timestep', Time

   !Recompute mesh bubbles etc
   CALL MeshStabParams( Model % Mesh )

   !Release the old mesh
   CALL ReleaseMesh(GatheredMesh)

   ! Release GatheredMesh % Redistribution
   IF(ASSOCIATED(GatheredMesh % Repartition)) THEN
      DEALLOCATE(GatheredMesh % Repartition)
      GatheredMesh % Repartition => NULL()
   END IF

   !remove mesh update
   CALL ResetMeshUpdate(Model, Solver)

   !pause solvers?
   IF(PauseAfterCalving) THEN
    IF(ImBoss) THEN
        PauseVolumeThresh = ListGetConstReal(SolverParams, "Pause Solvers Minimum Iceberg Volume", Found)
        IF(.NOT. Found) THEN
          CALL WARN(SolverName, "'Pause Solvers Minimum Iceberg Volume' not provided do assuming is 0")
          PauseVolumeThresh = 0.0_dp
        END IF

        PauseSolvers = MaxBergVolume > PauseVolumeThresh .AND. RSuccess
    END IF

    CALL MPI_BCAST(PauseSolvers, 1, MPI_LOGICAL, my_cboss, ELMER_COMM_WORLD, ierr)

    CALL PauseCalvingSolvers(Model, SolverParams, PauseSolvers)
  END IF

END SUBROUTINE CalvingRemeshMMG

SUBROUTINE CheckFlowConvergenceMMG( Model, Solver, dt, Transient )

  USE CalvingGeometry

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !-------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t) :: RemeshSolver
  TYPE(Variable_t), POINTER :: FlowVar, TimeVar
  TYPE(ValueList_t), POINTER :: Params, FuncParams
  LOGICAL :: Parallel, Found, CheckFlowDiverge=.TRUE., CheckFlowMax, FirstTime=.TRUE.,&
       NSDiverge, NSFail, NSTooFast, ChangeMesh, NewMesh, CalvingOccurs, GotSaveMetric,&
       CalvingPauseSolvers, PNSFail=.FALSE.
  REAL(KIND=dp) :: SaveNewtonTol, MaxNSDiverge, MaxNSValue, FirstMaxNSValue, FlowMax,&
       SaveFlowMax, Mag, NSChange, SaveDt, SaveRelax,SaveMeshHmin,SaveMeshHmax,&
       SaveMeshHgrad,SaveMeshHausd, SaveMetric, MeshChange, NewMetric, NewMeshDist,&
       SaveMeshDist, NewMeshHGrad
  REAL(KIND=dp), ALLOCATABLE :: SaveMeshHausdArray(:,:), SaveMeshHminArray(:,:), &
       NewMeshHausdArray(:,:), NewMeshHminArray(:,:)
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:), WorkArray(:,:)
  INTEGER :: i,j,SaveNewtonIter,Num, ierr, FailCount, TimeIntervals, SaveSSiter,&
       MaxRemeshIter,SaveNLIter,CurrentNLIter,NewMeshCount
  CHARACTER(MAX_NAME_LEN) :: FlowVarName, SolverName, EqName, RemeshEqName

  SAVE ::SaveNewtonTol, SaveNewtonIter, SaveFlowMax, FirstTime, FailCount,&
       SaveRelax,SaveMeshHminArray,SaveMeshHmax,SaveMeshHausdArray,SaveMeshHgrad, &
       SaveSSiter, MaxRemeshIter, SaveNLIter, NewMesh, NewMeshCount,&
       SaveMetric, GotSaveMetric, MeshChange, NewMeshHausdArray, NewMeshHminArray,&
       NewMetric, SaveMeshDist, NewMeshDist, NewMeshHGrad, PNSFail

  Mesh => Solver % Mesh
  SolverName = 'CheckFlowConvergenceMMG'
  Params => Solver % Values
  Parallel = (ParEnv % PEs > 1)
  FuncParams => GetMaterial(Mesh % Elements(1)) !TODO, this is not generalised
  FlowVarName = ListGetString(Params,'Flow Solver Name',Found)
  IF(.NOT. Found) FlowVarName = "Flow Solution"
  FlowVar => VariableGet(Mesh % Variables, FlowVarName, .TRUE., UnfoundFatal=.TRUE.)

  RemeshEqName = ListGetString(Params,'Remesh Equation Name',Found)
  IF(.NOT. Found) RemeshEqName = "remesh"

  !Get a handle to the remesh solver
  Found = .FALSE.
  DO j=1,Model % NumberOfSolvers
    IF(ListGetString(Model % Solvers(j) % Values, "Equation") == RemeshEqName) THEN
      RemeshSolver = Model % Solvers(j)
      Found = .TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT. Found) CALL Fatal(SolverName, "Failed to get handle to Remesh Solver.")

  IF(FirstTime) THEN

    FailCount = 0
    NewMeshCount = 0

    SaveNewtonTol = ListGetConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", Found)
    IF(.NOT. Found) SaveNewtonTol = 1.0E-3
    SaveNewtonIter = ListGetInteger(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", Found)
    IF(.NOT. Found) SaveNewtonIter = 15
    SaveNLIter = ListGetInteger(FlowVar % Solver % Values, &
         "Nonlinear System Max Iterations", Found)
    IF(.NOT. Found) SaveNewtonIter = 50

    SaveRelax = ListGetConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Relaxation Factor", Found)
    IF(.NOT. Found) SaveRelax = 1.0
    PRINT *, 'TO DO: replace with MMG stuff, hausdorff?'
    ! Use RemeshMMG3D Hmin in Material or Mesh Hmin in RemeshSolver
    ! should remesh without calving/cutting as well, so take RemeshMMG3D
    !SaveMeshHmin = ListGetConstReal(FuncParams, "MMG Hmin", Found, UnfoundFatal=.TRUE.)
    SaveMeshHmax = ListGetConstReal(FuncParams, "MMG Hmax", Found, UnfoundFatal=.TRUE.)
    !SaveMeshHausd = ListGetConstReal(FuncParams, "MMG Hausd", Found, UnfoundFatal=.TRUE.)
    SaveMeshHgrad = ListGetConstReal(FuncParams, "MMG Hgrad", Found, UnfoundFatal=.TRUE.)
    SaveMeshDist = ListGetConstReal(RemeshSolver % Values, "Remeshing Distance", Found, UnfoundFatal=.TRUE.)

    WorkArray => ListGetConstRealArray(FuncParams, "MMG Hmin", Found)
    IF(.NOT. Found) CALL FATAL(SolverName, 'Provide hmin input array to be iterated through: "Mesh Hmin"')
    MaxRemeshIter= SIZE(WorkArray(:,1))

    PRINT*, 'workarraysize', SIZE(WorkArray(1,:)), SIZE(WorkArray(:,1))

    ALLOCATE(SaveMeshHminArray(MaxRemeshIter, 1), NewMeshHminArray(MaxRemeshIter,1))
    SaveMeshHminArray = WorkArray
    NULLIFY(WorkArray)
    WorkArray => ListGetConstRealArray(FuncParams, "MMG Hausd", Found)
    IF(.NOT. Found) CALL FATAL(SolverName, 'Provide hmin input array to be iterated through: "Mesh Hausd"')
    IF(MaxRemeshIter /= SIZE(WorkArray(:,1))) CALL FATAL(SolverName, 'The number of hmin options &
            must equal the number of hausd options')
    ALLOCATE(SaveMeshHausdArray(MaxRemeshIter,1), NewMeshHausdArray(MaxRemeshIter,1))
    SaveMeshHausdArray = WorkArray
    NULLIFY(WorkArray)

    SaveMetric = ListGetConstReal( FuncParams, 'GlacierMeshMetric Min LC', GotSaveMetric)

    MeshChange = ListGetConstReal( Params, 'Mesh Change Percentage', Found)
    IF(.NOT. Found) THEN
      MeshChange = 10.0_dp
      CALL WARN(SolverName, "Not found Mesh Percentage Change so assuming 10%")
    END IF
    MeshChange = (MeshChange/100.0_dp)+1

    SaveSSiter = ListGetInteger(Model % Simulation, "Steady State Max Iterations", Found)
    IF(.NOT. Found) SaveSSiter = 1

    TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation, &
         'Timestep Sizes', Found, UnfoundFatal=.TRUE.)
    IF(SIZE(TimestepSizes,1) > 1) CALL Fatal(SolverName,&
         "Calving solver requires a single constant 'Timestep Sizes'")
    SaveDt = TimestepSizes(1,1)

    NewMeshHminArray = SaveMeshHminArray
    NewMeshHausdArray = SaveMeshHausdArray
    NewMetric = SaveMetric
    NewMeshDist = SaveMeshDist
    NewMeshHGrad = SaveMeshHgrad
  ELSE
    SaveDt = ListGetCReal( CurrentModel % Simulation,'Timestep Size' )
  END IF

  ! since "Calving solver requires a single constant 'Timestep Sizes'"
  TimeIntervals = ListGetInteger(Model % Simulation, "Timestep Intervals", UnfoundFatal = .TRUE.)

  !Get current simulation time
  TimeVar => VariableGet(Model % Variables, 'Time', .TRUE.)

  MaxNSDiverge = ListGetConstReal(Params, "Maximum Flow Solution Divergence", CheckFlowDiverge)
  MaxNSValue = ListGetConstReal(Params, "Maximum Velocity Magnitude", CheckFlowMax)
  FirstMaxNSValue = ListGetConstReal(Params, "First Time Max Expected Velocity", Found)
  IF(.NOT. Found .AND. CheckFlowDiverge) THEN
    CALL Info(SolverName, "'First Time Max Expected Velocity' not found, setting to 1.0E4")
    FirstMaxNSValue = 1.0E4
  END IF

  !====================================================!
  !---------------- DO THINGS! ------------------------!
  !====================================================!

  NSFail=.FALSE.
  NSDiverge=.FALSE.
  NSTooFast=.FALSE.

  !In addition to checking for absolute failure (% NonlinConverged = 0), we can check
  !for suspiciously large shift in the max variable value (this usually indicates a problem)
  !and we can also check for unphysically large velocity values
  IF(CheckFlowDiverge .OR. CheckFlowMax) THEN

    FlowMax = 0.0_dp
    DO i=1,Mesh % NumberOfNodes
      Mag = 0.0_dp

      DO j=1,FlowVar % DOFs-1
        Mag = Mag + (FlowVar % Values( (FlowVar % Perm(i)-1)*FlowVar % DOFs + j ) ** 2.0_dp)
      END DO
      Mag = Mag ** 0.5_dp
      FlowMax = MAX(FlowMax, Mag)
    END DO

    CALL MPI_AllReduce(MPI_IN_PLACE, FlowMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ELMER_COMM_WORLD, ierr)
  END IF

  IF(CheckFlowDiverge) THEN
    !First time, there's no previous velocity against which to check divergence.
    !This is somewhat messy because of the separate 'Maximum Velocity Magnitude'
    IF(FirstTime) SaveFlowMax = MIN(FlowMax,FirstMaxNSValue)

    NSChange = FlowMax / SaveFlowMax
    IF(ParEnv % MyPE == 0) PRINT *,'Debug, Flow Max (old/new): ',SaveFlowMax, FlowMax,' NSChange: ',NSChange

    IF(NSChange > MaxNSDiverge) THEN
      NSDiverge = .TRUE.
      CALL Info(SolverName,"Large change in maximum velocity suggests dodgy&
           &Navier-Stokes solution.")
    END IF
    IF(.NOT. NSDiverge) SaveFlowMax = FlowMax
  END IF

  IF(CheckFlowMax) THEN
    IF(FlowMax > MaxNSValue) THEN
      NSTooFast = .TRUE.
      CALL Info(SolverName,"Large maximum velocity suggests dodgy&
           &Navier-Stokes solution.")
    END IF
  END IF

  NSFail = FlowVar % NonlinConverged < 1 .OR. NSDiverge .OR. NSTooFast
  ! Joe note: I commented out Eef's testing here during merge:
  ! PRINT *, 'temporarily set NSFail=True for testing'
  ! NSFail=.TRUE.
  IF(ParEnv % MyPE == 0) PRINT*, 'Debug', FlowVar % NonlinConverged, NSDiverge, NSTooFast
  IF(NSFail) THEN
    CALL Info(SolverName, "Skipping solvers except Remesh because NS failed to converge.")

    FailCount = FailCount + 1
    IF(ParEnv % MyPE == 0) PRINT *, 'FailCount=',FailCount
    ! Joe note: I commented out Eef's testing here during merge:
    ! PRINT *, 'Temporarily set failcount to 2, to force remeshing!'
    ! FailCount=2
    IF(NewMeshCount > 3) THEN
       CALL Fatal(SolverName, "Don't seem to be able to recover from NS failure, giving up...")
    END IF

    !Set the clock back one second less than a timestep size.
    !This means next time we are effectively redoing the same timestep
    !but without any solvers which are dependent on (t > told) to reallocate
    TimeVar % Values(1) = TimeVar % Values(1) - SaveDt + (1.0/(365.25 * 24 * 60 * 60.0_dp))

    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", 0.0_dp)
    CALL ListAddInteger( FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", 10000)

    IF(failcount > 1) ChangeMesh = .TRUE.
    IF(.NOT. NSDiverge .AND. .NOT. NSTooFast) THEN
      !ie fails from nonconvergence but flow looks ok
      CurrentNLIter = ListGetInteger(FlowVar % Solver % Values, &
          "Nonlinear System Max Iterations", Found)
      IF(CurrentNLIter >= SaveNLIter*4) THEN !clearly not going to converge
        CALL ListAddInteger( FlowVar % Solver % Values, &
        "Nonlinear System Max Iterations", SaveNLIter)
        !ChangeMesh = .FALSE.
      ELSE
        CALL ListAddInteger( FlowVar % Solver % Values, &
            "Nonlinear System Max Iterations", CurrentNLIter*2)
        !ChangeMesh = .FALSE.
      END IF
    END IF

    !If this is the second failure in a row, fiddle with the mesh
    !PRINT *, 'TO DO, optimize MMG parameters, need to change remesh distance as well?'
    IF(ChangeMesh) THEN
       NewMeshHminArray = NewMeshHminArray/MeshChange
       NewMeshHausdArray = NewMeshHausdArray/MeshChange
       NewMetric = NewMetric/MeshChange
       NewMeshDist = NewMeshDist*MeshChange
       NewMeshHGrad = NewMeshHGrad/MeshChange

       CALL Info(SolverName,"NS failed twice, fiddling with the mesh... ")
       CALL Info(SolverName,"Temporarily slightly change MMG params ")
       CALL ListAddConstRealArray(FuncParams, "MMG Hmin", MaxRemeshIter, 1, NewMeshHminArray)
       !CALL ListAddConstReal(FuncParams, "MMG Hmax", SaveMeshHmax*1.1_dp)
       CALL ListAddConstReal(FuncParams, "MMG Hgrad", NewMeshHGrad) !default 1.3
       CALL ListAddConstReal(RemeshSolver % Values, "Remeshing Distance", NewMeshDist)
       CALL ListAddConstRealArray(FuncParams, "MMG Hausd", MaxRemeshiter, 1, NewMeshHausdArray)
       CALL ListAddInteger( FlowVar % Solver % Values, "Nonlinear System Max Iterations", SaveNLIter)
       IF(GotSaveMetric) CALL ListAddConstReal(FuncParams, "GlacierMeshMetric Min LC", NewMetric)
      NewMesh = .TRUE.
      NewMeshCount = NewMeshCount + 1
    ELSE
      NewMesh = .FALSE.
    END IF

    IF( .NOT. (NSTooFast .OR. NSDiverge)) THEN
      !---Not quite converging---!

      CALL ListAddConstReal(FlowVar % Solver % Values, &
           "Nonlinear System Relaxation Factor", 0.9_dp)

    ELSE
      !---Solution blowing up----!

      !Set var values to zero so it doesn't mess up viscosity next time
      FlowVar % Values = 0.0_dp

      !TODO: What else? different linear method? more relaxation?
    END IF

    !prevent calving on next time step
    CalvingOccurs = .FALSE.
    CALL SParIterAllReduceOR(CalvingOccurs)
    CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )
  ELSE
! set original values back
    FailCount = 0
    NewMeshCount = 0

    NewMeshHminArray = SaveMeshHminArray
    NewMeshHausdArray = SaveMeshHausdArray
    NewMetric = SaveMetric
    NewMeshDist = SaveMeshDist
    NewMeshHGrad = SaveMeshHgrad

    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Newton After Tolerance", SaveNewtonTol)
    CALL ListAddInteger( FlowVar % Solver % Values, &
         "Nonlinear System Newton After Iterations", SaveNewtonIter)
    CALL ListAddConstReal(FlowVar % Solver % Values, &
         "Nonlinear System Relaxation Factor", SaveRelax)
    CALL ListAddInteger( FlowVar % Solver % Values, "Nonlinear System Max Iterations", SaveNLIter)
    CALL ListAddConstRealArray(FuncParams, "MMG Hmin", MaxRemeshIter, 1, SaveMeshHminArray)
    !CALL ListAddConstReal(FuncParams, "MMG Hmax", SaveMeshHmax)
    CALL ListAddConstReal(RemeshSolver % Values, "Remeshing Distance", SaveMeshDist)
    CALL ListAddConstReal(FuncParams, "MMG Hgrad", SaveMeshHgrad)
    CALL ListAddConstRealArray(FuncParams, "MMG Hausd", MaxRemeshIter, 1, SaveMeshHausdArray)
    IF(GotSaveMetric) CALL ListAddConstReal(FuncParams, "GlacierMeshMetric Min LC", SaveMetric)

  END IF

  !Set a simulation switch to be picked up by Remesher
  CALL ListAddLogical( Model % Simulation, 'Flow Solution Failed', NSFail )

  CalvingPauseSolvers = ListGetLogical( Model % Simulation, 'Calving Pause Solvers', Found )
  IF(.NOT. Found) THEN
    CALL INFO(SolverName, 'Cannot find Calving Pause Solvers so assuming calving has not pasued time.')
    CalvingPauseSolvers = .TRUE.
  END IF
  IF(.NOT. CalvingPauseSolvers .OR. FirstTime .OR. NSFail .OR. PNSFail) THEN
    !Switch off solvers
    DO Num = 1,999
      WRITE(Message,'(A,I0)') 'Switch Off Equation ',Num
      EqName = ListGetString( Params, Message, Found)
      IF( .NOT. Found) EXIT

      Found = .FALSE.
      DO j=1,Model % NumberOfSolvers
        IF(ListGetString(Model % Solvers(j) % Values, "Equation") == EqName) THEN
          Found = .TRUE.
          !Turn off (or on) the solver
          !If NS failed to converge, (switch) off = .true.
          CALL SwitchSolverExec(Model % Solvers(j), NSFail)
          EXIT
        END IF
      END DO

      IF(.NOT. Found) THEN
        WRITE (Message,'(A,A,A)') "Failed to find Equation Name: ",EqName,&
            " to switch off after calving."
        CALL Fatal(SolverName,Message)
      END IF
    END DO
  END IF

  IF(NSFail) THEN
    !CALL ListAddConstReal( Model % Simulation, 'Timestep Size', PseudoSSdt)
    !CALL ListAddInteger( Model % Simulation, 'Steady State Max Iterations', 1)
    CALL ListAddInteger( Model % Simulation, 'Timestep Intervals', TimeIntervals + 1)
  ELSE
    !CALL ListAddConstReal( Model % Simulation, 'Timestep Size', SaveDt)
    CALL ListAddInteger( Model % Simulation, 'Steady State Max Iterations', SaveSSiter)
  END IF

  FirstTime = .FALSE.
  PNSFail = NSFail

END SUBROUTINE CheckFlowConvergenceMMG
