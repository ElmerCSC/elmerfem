MODULE MeshRemeshing

USE DefUtils
USE Types
USE MeshUtils

IMPLICIT NONE

#ifdef HAVE_MMG
#include "mmg/mmg3d/libmmg3df.h"
#endif

#ifdef HAVE_MMG
INTEGER :: MMGPARAM_hausd = MMG3D_DPARAM_hausd
INTEGER :: MMGPARAM_hmin = MMG3D_DPARAM_hmin
INTEGER :: MMGPARAM_hmax = MMG3D_DPARAM_hmax
INTEGER :: MMGPARAM_iso = MMG3D_IPARAM_iso
INTEGER :: MMGPARAM_hgrad = MMG3D_DPARAM_hgrad
INTEGER :: MMGPARAM_angle = MMG3D_IPARAM_angle
MMG5_DATA_PTR_T :: mmgMesh
MMG5_DATA_PTR_T :: mmgSol
#endif

CONTAINS

!============================================
!============================================
!            MMG3D SUBROUTINES
!============================================
!============================================


SUBROUTINE Set_MMG3D_Mesh(Mesh, Parallel)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel

#ifdef HAVE_MMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ref,ierr
  INTEGER, ALLOCATABLE :: NodeRefs(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Set_MMG3D_Mesh"
  IF(CoordinateSystemDimension() /= 3) CALL Fatal("MMG3D","Only works for 3D meshes!")

  ALLOCATE(NodeRefs(6))

  IF(Parallel) CALL Assert(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs), FuncName,&
       "Parallel sim but no ParallelInfo % GlobalDOFs")

  nverts = Mesh % NumberOfNodes
  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0

  nbulk = Mesh % NumberOfBulkElements
  nbdry = Mesh % NumberOfBoundaryElements

  !Cycle mesh elements gathering type count
  DO i=1,nbulk+nbdry
    Element => Mesh % Elements(i)
    SELECT CASE(Element % TYPE % ElementCode)
    CASE(101)
      Warn101 = .TRUE.
    CASE(202)
      Warn202 = .TRUE.
    CASE(303)
      NTris = NTris + 1
    CASE(404)
      NQuads = NQuads + 1
    CASE(504)
      NTetras = NTetras + 1
    CASE(605)
      CALL Fatal("MMG3D","can't handle pyramid elements (605)")
    CASE(706)
      NPrisms = NPrisms + 1
    CASE(808)
      CALL Fatal("MMG3D","can't handle brick/hexahedral elements (808)")
    CASE DEFAULT
      PRINT *,'Bad element type: ',Element % TYPE % ElementCode
      CALL Fatal("MMG3D","Unsupported element type")
    END SELECT
  END DO

  IF(Warn101) CALL Warn("MMG3D","101 elements detected - these won't be remeshed")
  IF(Warn202) CALL Warn("MMG3D","202 elements detected - these won't be remeshed")

  !args: mesh, nvertices, ntetra, nprisms, ntriangless, nquads, nedges
  CALL MMG3D_Set_meshSize(mmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ierr)
  IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
       'CALL TO MMG3D_Set_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_meshSize DONE'

  ref = 0
  DO i=1,NVerts
    IF(Parallel) ref = Mesh % ParallelInfo % GlobalDOFs(i)
    CALL MMG3D_Set_vertex(mmgMesh, Mesh%Nodes%x(i), &
         Mesh%Nodes%y(i),Mesh%Nodes%z(i), ref, i, ierr)
!         Mesh%Nodes%y(i),Mesh%Nodes%z(i), 0, Mesh % ParallelInfo % GlobalDOFs(i), ierr)
    !PRINT *,'debug: mesh point: ',Mesh%Nodes%x(i), &
    !     Mesh%Nodes%y(i),Mesh%Nodes%z(i), i

    IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_Set_vertex FAILED')
  END DO
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_vertex DONE'

  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0

  !Cycle mesh elements, sending them to MMG3D
  DO i=1,nbulk+nbdry
    Element => Mesh % Elements(i)
    NNodes = Element % TYPE % NumberOfNodes

    NodeIndexes => Element % NodeIndexes
    NodeRefs(1:NNodes) = NodeIndexes(1:NNodes)
!    NodeRefs(1:NNodes) = Mesh % ParallelInfo % GlobalDOFs(NodeIndexes(1:NNodes))

!    PRINT *,'debug, elem ',i,' noderefs: ',NodeRefs(1:NNodes)
    SELECT CASE(Element % TYPE % ElementCode)
    CASE(303)
      ntris = ntris + 1
      CALL MMG3D_Set_triangle(mmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           Element % BoundaryInfo % Constraint, ntris, ierr)
    CASE(404)
      nquads = nquads + 1
      CALL MMG3D_Set_quadrilateral(mmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  Element % BoundaryInfo % Constraint, nquads,ierr)
    CASE(504)
      ntetras = ntetras + 1
      CALL MMG3D_Set_tetrahedron(mmgMesh,NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  Element % BodyID, ntetras,ierr)
    CASE(706)
      nprisms = nprisms + 1
      CALL MMG3D_Set_Prism(mmgMesh, NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4), NodeRefs(5), NodeRefs(6), Element % BodyID, nprisms,ierr)
    CASE DEFAULT
    END SELECT
  END DO
  IF (DEBUG) PRINT *,'--**-- MMG3D - Set elements DONE'

#else
     CALL FATAL('Set_MMG3D_Mesh',&
        'Remeshing utility MMG3D has not been installed')
#endif

END SUBROUTINE Set_MMG3D_Mesh

SUBROUTINE Set_MMG3D_Parameters(SolverParams)

  TYPE(ValueList_t), POINTER :: SolverParams

#ifdef HAVE_MMG
  REAL(KIND=dp) :: Pval
  LOGICAL :: AngleDetect
  INTEGER :: verbosity,MeMIncrease,Bucket,GMshOption,ierr
  REAL(KIND=dp) :: Hmin, Hmax, HSiz
  LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
  LOGICAL :: Found,Debug=.FALSE.

  ! Minimal mesh size:  hmin
  Hmin = GetConstReal( SolverParams, 'hmin', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,&
         Hmin,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hmin> Failed')
  END IF
  ! Maximal mesh size - hmax
  Hmax = GetConstReal( SolverParams, 'hmax', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,&
         Hmax,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hmax> Failed')
  END IF

  hsiz = GetConstReal( SolverParams, 'hsiz', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hsiz,&
         hsiz,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hsiz> Failed')
  END IF

!!! PARAMS: generic options (debug, mem, verbosity)
  ! [val] Set the verbosity level to n
  Verbosity=GetInteger( SolverParams,'verbosity',Found)
  IF (Found) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, &
         Verbosity,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
  ! [val] Set the maximal memory size to n MBytes.
  MemIncrease=GetInteger(SolverParams,'Increase Memory',Found)
  IF (FOUND) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_mem,&
         MemIncrease,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
  ! [0/1] Turn on the debug mode.
  DebugMode=GetLogical(SolverParams,'Debug Mode',Found)
  IF (Found .AND. DebugMode) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_debug,1, &
         ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
!!! PARAMS
  ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
  ! MMG3D_DPARAM_hausd default est 0.01 semble bien trop petit;
  !  il semble qu'il faille mettre une taille > taille des elements.
  Pval=GetConstReal( SolverParams, 'hausd', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,&
         Pval,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hausd> Failed')
  END IF
  ! Control gradation
  Pval=GetConstReal( SolverParams, 'hgrad', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,&
         Pval,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hgrad> Failed')
  END IF

!!! OTHER PARAMETERS: NOT ALL TESTED
  Pval = GetConstReal( SolverParams, 'Angle detection',Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
         Pval,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <Angle detection> Failed')
  ENDIF
  ! !< [1/0], Avoid/allow surface modifications */ 
  AngleDetect=GetLogical(SolverParams,'No Angle detection',Found)
  IF (Found.AND.AngleDetect) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle, &
         1,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER <No Angle detection> Failed')
  END IF
  ! [1/0] Avoid/allow point insertion
  NoInsert=GetLogical(SolverParams,'No insert',Found)
  IF (Found .AND. NoInsert) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noinsert,&
         1,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER <No insert> Failed') 
  END IF
  ! [1/0] Avoid/allow edge or face flipping
  NoSwap=GetLogical(SolverParams,'No swap',Found)
  IF (Found .AND. NoSwap) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noswap,&
         1,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No swap>Failed')
  END IF
  ! [1/0] Avoid/allow point relocation
  NoMove=GetLogical(SolverParams,'No move',Found)
  IF (Found .AND. NoMove) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nomove,&
         1,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No move> Failed')
  END IF
  ! [1/0] Avoid/allow surface modifications
  NoSurf=GetLogical(SolverParams,'No surf',Found)
  IF (Found .AND. NoSurf) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nosurf,&
         0,ierr)
    IF ( ierr == 0 ) CALL FATAL('M/MGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No surf> Failed')
  END IF
!!!
#else
     CALL FATAL('Set_MMG3D_Parameters',&
        'Remeshing utility MMG3D has not been installed')
#endif
END SUBROUTINE Set_MMG3D_Parameters


SUBROUTINE Get_MMG3D_Mesh(NewMesh, Parallel, FixedNodes, FixedElems)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  LOGICAL :: Parallel
  LOGICAL, OPTIONAL, ALLOCATABLE :: FixedNodes(:), FixedElems(:)
  !------------------------------------------------------------------------------

#ifdef HAVE_MMG

  TYPE(Element_t),POINTER ::  Element
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent,ied
  INTEGER :: ii,kk
  LOGICAL :: Debug

  !> a) get the size of the mesh: vertices, tetra,prisms, triangles, quads,edges
  CALL MMG3D_Get_meshSize(mmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)
  IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
       'CALL TO MMGS_Get_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Get_meshSize DONE'    

  ! INITIALISE THE NEW MESH STRUCTURE
  !NPrisms, NQuads should be zero
  NewMesh => AllocateMesh( NTetras, NTris, NVerts, .TRUE.)

  NewMesh % Name = "MMG3D_Output"
  NewMesh%MaxElementNodes=4
  NewMesh%MaxElementDOFs=4
  NewMesh%MeshDim=3

  IF(PRESENT(FixedNodes)) THEN
    ALLOCATE(FixedNodes(NVerts))
    FixedNodes = .FALSE.
  END IF
  IF(PRESENT(FixedElems)) THEN
    ALLOCATE(FixedElems(NTetras+NTris))
    FixedNodes = .FALSE.
  END IF

  IF(NPrisms /= 0) CALL Fatal("MMG3D", "Programming Error: MMG3D returns prisms")
  IF(NQuads /= 0) CALL Fatal("MMG3D", "Programming Error: MMG3D returns quads")
  NewMesh % NumberOfNodes = NVerts
  NewMesh % NumberOfBulkElements = NTetras
  NewMesh % NumberOfBoundaryElements = NTris

  IF(Parallel) THEN
    NewMesh % ParallelInfo % Interface = .FALSE.
  END IF

  !! GET NEW VERTICES
  NewMesh % Nodes % z = 0._dp
  Do ii=1,NVerts
    CALL MMG3D_Get_vertex(mmgMesh,&
         NewMesh % Nodes % x(ii),&
         NewMesh % Nodes % y(ii),&
         NewMesh % Nodes % z(ii),&
         ref,corner,required,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_vertex FAILED')
    IF(Parallel) THEN
      IF(required > 0) THEN
        NewMesh % ParallelInfo % GlobalDOFs(ii) = ref
      ELSE
        !GlobalDOF undefined - need to negotiate w/ other parts
        NewMesh % ParallelInfo % GlobalDOFs(ii) = 0
      END IF
    END IF
    IF(PRESENT(FixedNodes)) FixedNodes(ii) = required > 0
  End do

  IF (DEBUG) PRINT *,'MMG3D_Get_vertex DONE'

  !! GET NEW TETRAS
  DO ii=1,NewMesh % NumberOfBulkElements
    Element => NewMesh % Elements(ii)
    Element % TYPE => GetElementType( 504 )
    Element % NDOFs = Element % TYPE % NumberOfNodes
    Element % ElementIndex = ii
    Element % PartIndex = ParEnv % myPE
    CALL AllocateVector(Element % NodeIndexes, 4)
    NodeIndexes => Element % NodeIndexes
    CALL MMG3D_Get_tetrahedron(mmgMesh, &
         NodeIndexes(1), &
         NodeIndexes(2), &
         NodeIndexes(3), &
         NodeIndexes(4), &
         Element % BodyId, & !TODO - many tetras end up with very high BodyIDs
         required,ierr)
    IF(PRESENT(FixedElems)) FixedElems(ii) = required > 0
  END DO
  IF (DEBUG) PRINT *,'MMG3D_Get_tets DONE'

  !! Get BC Elements
  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1

    Element => NewMesh % Elements(kk)

    Element % TYPE => GetElementType( 303 )
    Element % NDOFs = Element % TYPE % NumberOfNodes
    Element % ElementIndex = kk
    Element % PartIndex = ParEnv % myPE

    CALL AllocateVector(Element % NodeIndexes, 3)
    NodeIndexes => Element % NodeIndexes
    CALL MMG3D_Get_triangle(mmgMesh, &
         NodeIndexes(1), &
         NodeIndexes(2), &
         NodeIndexes(3), &
         ref   , &
         required,ierr)
    IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_triangle FAILED')

    Allocate(Element % BoundaryInfo)
    Element % BoundaryInfo % Constraint=ref

    IF(PRESENT(FixedElems)) FixedElems(kk) = required > 0
  END DO

  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1

    CALL MMG3D_GET_TetFromTria(mmgMesh, &
         ii,parent,ied,ierr)

    IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_TetFromTria FAILED')
    Element => NewMesh % Elements(kk)
    Element % BoundaryInfo % Left => NewMesh % Elements(parent) !TODO - parent ID offset?
  END DO
#else
     CALL FATAL('Get_MMG3D_Mesh',&
        'Remeshing utility MMG3D has not been installed')
#endif

END SUBROUTINE Get_MMG3D_Mesh

!Subroutine to negotiate new global node numbers between partitions
!as a necessary precursor to repartitioning the mesh.
!Expects to receive OldMesh with valid GlobalDOFs, and new mesh 
!in which all GlobalDOFs are either present in the OldMesh, or set to zero.
!We allow here that each partition may have any number of nodes (inc. zero)
!Assumes that NewMesh doesn't have any nodes which are both *shared* and *unmarked*
SUBROUTINE RenumberGDOFs(OldMesh,NewMesh)
  TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
  !----------------------
  INTEGER :: i,j,k,n,counter, OldNN, NewNN, nglobal_pool, nlocal_pool,&
       ierr, my_maxgdof, mingdof, maxgdof, request,unused,&
       Need, PNeed(ParEnv % PEs), pnlocal_pool(ParEnv % PEs), disps(ParEnv % PEs), &
       PNeed_tot
  INTEGER, ALLOCATABLE :: old_gdofs(:), new_gdofs(:), global_pool(:),&
       local_pool(:), work_int(:), GDOF_remap(:), send_to(:)
  INTEGER, POINTER :: ngdof_ptr(:), ogdof_ptr(:)
  LOGICAL :: Root, Debug = .FALSE.
  LOGICAL, ALLOCATABLE :: AvailGDOF(:), pool_duplicate(:), im_using(:), used(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="RenumberGDOFs"

  Root = ParEnv % MyPE == 0
  OldNN = OldMesh % NumberOfNodes
  NewNN = NewMesh % NumberOfNodes

  ngdof_ptr => NewMesh % ParallelInfo % GlobalDOFs
  ogdof_ptr => OldMesh % ParallelInfo % GlobalDOFs

  ALLOCATE(old_gdofs(OldNN), &
       new_gdofs(NewNN), &
       AvailGDOF(OldNN))

  AvailGDOF = .FALSE.

  old_gdofs = ogdof_ptr
  CALL Sort(OldNN,old_gdofs)
  new_gdofs = ngdof_ptr
  CALL Sort(NewNN,new_gdofs)

  my_maxgdof = old_gdofs(oldNN)
  CALL MPI_ALLREDUCE( my_maxgdof, maxgdof, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)

  IF(Debug) PRINT *,ParEnv % MyPE,' debug max gdof :', maxgdof

  !Locally check the consistency of the old mesh
  DO i=1,OldNN-1
    IF(old_gdofs(i) == old_gdofs(i+1)) &
         CALL Fatal(FuncName,"OldMesh has duplicate GlobalDOFs")
  END DO
  IF(ANY(old_gdofs <= 0)) CALL Fatal(FuncName,"OldMesh has at least 1 GlobalDOFs <= 0")

  !Determine locally which GlobalDOFs are available for assignment
  DO i=1, OldNN
    k = SearchI(NewNN, new_gdofs, old_gdofs(i))
    IF(k == 0) AvailGDOF(i) = .TRUE.
  END DO

  nlocal_pool = COUNT(AvailGDOF)
  ALLOCATE(local_pool(COUNT(AvailGDOF)))
  local_pool = PACK(old_gdofs,AvailGDOF)
  CALL Assert(ALL(local_pool > 0), "RenumberGDOFs", "Programming error: some local pool <= 0")

  !Gather pool of available global nodenums to root
  !--------------------------------------
  !Note: previously allowed partitions to fill from local pool at this stage
  !but this assumes that no previously shared nodes are destroyed. Not necessarily
  !the case.

  CALL MPI_GATHER(nlocal_pool, 1, MPI_INTEGER, pnlocal_pool, 1, &
       MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  IF(Root) THEN
    IF(Debug) PRINT *,'pnlocal_pool: ',pnlocal_pool
    nglobal_pool = SUM(pnlocal_pool)
    ALLOCATE(global_pool(nglobal_pool), pool_duplicate(nglobal_pool))
    pool_duplicate = .FALSE.
    disps(1) = 0
    DO i=2,ParEnv % PEs
      disps(i) = disps(i-1) + pnlocal_pool(i-1)
    END DO
  END IF

  CALL MPI_GatherV(local_pool, nlocal_pool, MPI_INTEGER, global_pool, &
       pnlocal_pool,disps, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  DEALLOCATE(local_pool)

  !Sort & remove duplicates from global pool
  IF(Root) THEN
    CALL Sort(nglobal_pool, global_pool)
    DO i=2,nglobal_pool
      IF(global_pool(i) == global_pool(i-1)) pool_duplicate(i) = .TRUE.
    END DO
    nglobal_pool = COUNT(.NOT. pool_duplicate)

    ALLOCATE(work_int(nglobal_pool))
    work_int = PACK(global_pool, .NOT. pool_duplicate)

    DEALLOCATE(global_pool)
    CALL MOVE_ALLOC(work_int, global_pool)

    IF(Debug) THEN
      PRINT *,'Removed ',COUNT(pool_duplicate),' duplicate global pool entries'
      PRINT *,'new global pool: ',global_pool
    END IF
  END IF

  !Now we have a pool of nodenums (global_pool) no longer
  !used by each partition, but we need to check for those 
  !nodes which were simply passed from one partition to another
  !in OldMesh -> NewMesh
  CALL MPI_BCAST(nglobal_pool,1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
  IF(.NOT. Root) ALLOCATE(global_pool(nglobal_pool))
  ALLOCATE(im_using(nglobal_pool))
  im_using = .FALSE.
  CALL MPI_BCAST(global_pool,nglobal_pool, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  DO i=1, nglobal_pool
    k = SearchI(NewNN, new_gdofs, global_pool(i))
    IF(k /= 0) im_using(i) = .TRUE.
    IF(k /= 0 .AND. Debug) PRINT *,ParEnv % MyPE,' using ',global_pool(i),' from global pool'
  END DO

  !Logical reduction to determine which global dofs are actually in use
  IF(Root) ALLOCATE(used(nglobal_pool))
  CALL MPI_REDUCE(im_using, used, nglobal_pool, MPI_LOGICAL, MPI_LOR, 0, ELMER_COMM_WORLD, ierr)

  IF(ROOT) THEN
    IF(Debug) PRINT *, COUNT(used),' are in use of ',nglobal_pool
    ALLOCATE(work_int(COUNT(.NOT. used)))
    work_int = PACK(global_pool, .NOT. used)
    DEALLOCATE(global_pool)
    CALL MOVE_ALLOC(work_int, global_pool)
    nglobal_pool = SIZE(global_pool)
  ELSE
    DEALLOCATE(global_pool)
  END IF

  !Gather how many global nodenums required by each partition
  Need = COUNT(ngdof_ptr == 0)
  CALL MPI_GATHER(Need, 1, MPI_INTEGER, PNeed, 1, MPI_INTEGER, &
       0, ELMER_COMM_WORLD, ierr)

  !Either: 
  ! SUM(PNeed) > nglobal_pool : generate additional globalDOFs
  ! SUM(PNeed) < nglobal_pool : delete excess globalDOFs
  ! SUM(PNeed) == nglobal_pool : all good (rare)

  IF(Root) THEN
    PNeed_tot = SUM(PNeed)

    IF(PNeed_tot > nglobal_pool) THEN

      ALLOCATE(work_int(PNeed_tot))
      work_int = 0
      work_int(1:nglobal_pool) = global_pool
      DO i=1, PNeed_tot - nglobal_pool
        work_int(nglobal_pool + i) = maxgdof + i
      END DO
      DEALLOCATE(global_pool)
      CALL MOVE_ALLOC(work_int, global_pool)
      nglobal_pool = PNeed_tot
    END IF

    !Mark global dofs to send to each part
    !-1 means excess (to be destroyed)
    ALLOCATE(send_to(nglobal_pool))
    send_to = -1 !default
    counter = 1
    DO i=1,ParEnv % PEs
      send_to(counter:counter+PNeed(i)-1) = i-1
      counter = counter + PNeed(i)
    END DO
  END IF

  !Set up iRECV if required
  IF(Need > 0) THEN
    ALLOCATE(local_pool(Need))
    CALL MPI_iRECV(local_pool,Need,MPI_INTEGER,0,1001,&
       ELMER_COMM_WORLD, request, ierr)
  END IF

  !Root sends out nodes from pool
  IF(Root) THEN
    IF(Debug) PRINT *,'Pneed: ',PNeed
    counter = 1
    DO i=1,ParEnv % PEs
      IF(PNeed(i) <= 0) CYCLE

      IF(ANY(send_to(counter:counter+PNeed(i)-1) /= i-1)) &
           CALL Fatal(FuncName, "Programming error in send pool")

      CALL MPI_SEND(global_pool(counter:counter+PNeed(i)-1),PNeed(i),&
           MPI_INTEGER, i-1, 1001, ELMER_COMM_WORLD, ierr)

      counter = counter + PNeed(i)
    END DO
  END IF

  IF(Need > 0) CALL MPI_Wait(request, MPI_STATUS_IGNORE, ierr)

  !Fill from pool as required
  IF(Need > 0) THEN
    counter = 0
    DO i=1,NewNN
      IF(ngdof_ptr(i) == 0) THEN
        counter = counter + 1
        ngdof_ptr(i) = local_pool(counter)
      END IF
    END DO
  END IF


  !Root packs and sends unused gdofs
  IF(Root) THEN
    unused = COUNT(send_to == -1)
    IF(Debug) PRINT *,ParEnv % MyPE,' unused: ',unused

    ALLOCATE(work_int(unused))
    work_int = PACK(global_pool,  send_to == -1)
    DEALLOCATE(global_pool)
    CALL MOVE_ALLOC(work_int, global_pool)

    IF(Debug) THEN
      PRINT *,' unused, size(global_pool): ', unused, SIZE(global_pool)
      PRINT *,'Final global pool: ',global_pool
    END IF
  END IF

  CALL MPI_BCAST(unused,1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  IF(.NOT. Root) THEN
    ALLOCATE(global_pool(unused))
  END IF

  CALL MPI_BCAST(global_pool,unused,MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)

  !Renumber GDOFs to ensure contiguity if required
  IF(unused > 0) THEN
    IF(Debug) PRINT *,ParEnv % MyPE,' final global pool: ',global_pool

    MaxGDOF = MAXVAL(ngdof_ptr)
    MinGDOF = MINVAL(ngdof_ptr)
    IF(MinGDOF <= 0) CALL Fatal(FuncName, "Programming error: at least one gdof == 0")

    new_gdofs = ngdof_ptr
    CALL Sort(NewNN, new_gdofs)

    ALLOCATE(GDOF_remap(MinGDOF:MaxGDOF))
    DO i=1,NewNN
        GDOF_remap(new_gdofs(i)) = new_gdofs(i) - SearchIntPosition(global_pool, new_gdofs(i))
        IF(Debug) PRINT *,ParEnv % MyPE,' debug gdof map: ',new_gdofs(i), &
             SearchIntPosition(global_pool, new_gdofs(i)), GDOF_remap(new_gdofs(i))
    END DO

    ngdof_ptr = GDOF_remap(ngdof_ptr)
  END IF

END SUBROUTINE RenumberGDOFs

!Because elements aren't shared, renumbering is easier
!Simply number contiguously in each partition
!NOTE: won't work for halo elems
SUBROUTINE RenumberGElems(Mesh)
  TYPE(Mesh_t), POINTER :: Mesh
  !---------------------------------
  INTEGER :: i,NElem,MyStart,ierr
  INTEGER, ALLOCATABLE :: PNElem(:)

  NElem = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
  ALLOCATE(PNElem(ParEnv % PEs))

  CALL MPI_ALLGATHER(NElem,1,MPI_INTEGER, PNElem, 1, MPI_INTEGER, ELMER_COMM_WORLD, ierr)

  MyStart = SUM(PNElem(1:ParEnv % MyPE))
  DO i=1,NElem
    Mesh % Elements(i) % GElementIndex = MyStart + i
  END DO

END SUBROUTINE RenumberGElems

!Based on a previous mesh with valid nodal parallelinfo (% Interface & % NeighbourList)
!map that info onto NewMesh which shares at least some GlobalDOFs. Intended use is
!to enable reconnection of a parallel mesh part which has been remeshed internally, but
!whose partition boundaries remain as they were. e.g. CalvingRemeshMMG.F90
SUBROUTINE MapNewParallelInfo(OldMesh, NewMesh)
  TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
  !---------------------------------
  INTEGER, ALLOCATABLE :: GtoNewLMap(:)
  INTEGER :: i,k,n,MaxNGDof, MinNGDof
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="MapNewBCInfo"
  
  MinNGDof = HUGE(MinNGDof)
  MaxNGDof = 0
  DO i=1,NewMesh % NumberOfNodes
    k = NewMesh % ParallelInfo % GlobalDOFs(i)
    IF(k == 0) CYCLE
    MinNGDof = MIN(k,MinNGDof)
    MaxNGDof = MAX(k, MaxNGDof)
  END DO

  ALLOCATE(GtoNewLMap(MinNGDof:MaxNGDof))
  GtoNewLMap = 0

  DO i=1,NewMesh % NumberOfNodes
    k = NewMesh % ParallelInfo % GlobalDOFs(i)
    IF(k == 0) CYCLE
    GtoNewLMap(k) = i
  END DO

  DO i=1,OldMesh % NumberOfNodes
    IF(OldMesh % ParallelInfo % INTERFACE(i)) THEN
      k = OldMesh % ParallelInfo % GlobalDOFs(i)
      IF(k < LBOUND(GToNewLMap,1) .OR. k > UBOUND(GToNewLMap,1)) THEN
        CALL Warn(FuncName, "Interface node from OldMesh missing in NewMesh")
        CYCLE
      ELSEIF(GToNewLMap(k) == 0) THEN
        CALL Warn(FuncName, "Interface node from OldMesh missing in NewMesh")
        CYCLE
      END IF
      k = GToNewLMap(k)
      NewMesh % ParallelInfo % Interface(k) = .TRUE.

      n = SIZE(OldMesh % ParallelInfo % Neighbourlist(i) % Neighbours)
      IF(ASSOCIATED(NewMesh % ParallelInfo % Neighbourlist(k) % Neighbours)) &
           DEALLOCATE(NewMesh % ParallelInfo % Neighbourlist(k) % Neighbours)
      ALLOCATE(NewMesh % ParallelInfo % Neighbourlist(k) % Neighbours(n))
      NewMesh % ParallelInfo % Neighbourlist(k) % Neighbours = &
           OldMesh % ParallelInfo % Neighbourlist(i) % Neighbours
    END IF
  END DO

  DO i=1,NewMesh % NumberOfNodes
    IF(.NOT. ASSOCIATED(NewMesh % ParallelInfo % Neighbourlist(i) % Neighbours)) THEN
      ALLOCATE(NewMesh % ParallelInfo % Neighbourlist(i) % Neighbours(1))
      NewMesh % ParallelInfo % Neighbourlist(i) % Neighbours(1) = ParEnv % MyPE
    END IF
  END DO
END SUBROUTINE MapNewParallelInfo


!A subroutine for a 3D mesh (in serial) with MMG3D
!Inputs:
!   InMesh - the initial mesh
!   Metric - 2D real array specifying target metric
!   NodeFixed, ElemFixed - Optional mask to specify 'required' entities
!Output:
!   OutMesh - the improved mesh
!
SUBROUTINE RemeshMMG3D(InMesh,TargetLength,OutMesh,NodeFixed,ElemFixed)

  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:)
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  !-----------
  REAL(KIND=dp), ALLOCATABLE :: Metric(:,:)
  REAL(KIND=dp) :: hsiz(3)
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType
  LOGICAL :: Debug, Parallel, TensorMet
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName

#ifdef HAVE_MMG

  Debug = .TRUE.
  Parallel = ParEnv % PEs > 1
  FuncName = "RemeshMMG3D"

  !Scalar, vector, tensor metric?
  MetricDim = SIZE(TargetLength,2)
  NNodes = InMesh % NumberOfNodes
  NBulk = InMesh % NumberOfBulkElements
  NBdry = InMesh % NumberOfBoundaryElements

  IF(MetricDim == 1) THEN
    TensorMet = .FALSE.
  ELSE IF(MetricDim == 3) THEN
    TensorMet = .TRUE.
  ELSE
    CALL Fatal(FuncName,"Unexpected Metric Dimension!")
  END IF


  !Convert target length to metric
  IF(.NOT. TensorMet) THEN

    SolType = MMG5_Scalar
    ALLOCATE(Metric(NNodes,1))
    !TODO - check this - in anistropic case, Metric = edge length?
    Metric(:,1) = TargetLength(:,1)

  ELSE

    SolType = MMG5_Tensor
    !Upper triangle of symmetric tensor: 11,12,13,22,23,33
    ALLOCATE(Metric(NNodes,6))
    Metric = 0.0_dp
    DO i=1,NNodes
      !Metric = 1.0/(edge_length**2)
      Metric(i,1) = 1.0 / (TargetLength(i,1)**2.0)
      Metric(i,4) = 1.0 / (TargetLength(i,2)**2.0)
      Metric(i,6) = 1.0 / (TargetLength(i,3)**2.0)
    END DO

  END IF

!  hsiz = (/2000.0,2000.0,150.0/)
  mmgMesh = 0
  mmgSol  = 0

  CALL MMG3D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  CALL SET_MMG3D_MESH(InMesh,Parallel)

  !Set the metric values at nodes
  CALL MMG3D_Set_SolSize(mmgMesh, mmgSol, MMG5_Vertex, NNodes, SolType,ierr)
  IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")

  DO i=1,NNodes
    IF(TensorMet) THEN
      IF(Debug) PRINT *,'debug sol at ',i,' is: ',Metric(i,:)
      CALL MMG3D_Set_TensorSol(mmgSol,Metric(i,1),Metric(i,2),Metric(i,3),&
           Metric(i,4),Metric(i,5),Metric(i,6),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set tensor solution at vertex")
    ELSE
      CALL MMG3D_Set_ScalarSol(mmgSol,Metric(i,1),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
    END IF
  END DO

  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,&
       20.0_dp,ierr)

  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,&
         1.4_dp,ierr)

  !Take care of fixed nodes/elements if requested
  IF(PRESENT(NodeFixed)) THEN
    DO i=1,NNodes
      IF(NodeFixed(i)) THEN
        CALL MMG3D_SET_REQUIREDVERTEX(mmgMesh,i,ierr)
      END IF
    END DO
  END IF

  IF(PRESENT(ElemFixed)) THEN
    DO i=1,NBulk + NBdry
      IF(ElemFixed(i)) THEN
        IF(i <= NBulk) THEN
          CALL MMG3D_SET_REQUIREDTETRAHEDRON(mmgMesh,i,ierr)
        ELSE
          CALL MMG3D_SET_REQUIREDTRIANGLE(mmgMesh,i-NBulk,ierr)
        END IF
      END IF
    END DO
  END IF


  IF (DEBUG) PRINT *,'--**-- SET MMG3D PARAMETERS '
  ! CALL SET_MMG3D_PARAMETERS(SolverParams)

  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ierr)
  IF ( ierr == MMG5_STRONGFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
    STOP MMG5_STRONGFAILURE
  ELSE IF ( ierr == MMG5_LOWFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB"
  ENDIF
  IF (DEBUG) PRINT *,'--**-- MMG3D_mmg3dlib DONE'

  CALL MMG3D_SaveMesh(mmgMesh,"test_MMG3D_remesh.mesh",LEN(TRIM("test_MMG3D_remesh.mesh")),ierr)

  !! GET THE NEW MESH
  CALL GET_MMG3D_MESH(OutMesh,Parallel)



#else
  CALL Fatal(FuncName, "Remeshing utility MMG3D has not been installed")
#endif

END SUBROUTINE RemeshMMG3D

END MODULE MeshRemeshing
