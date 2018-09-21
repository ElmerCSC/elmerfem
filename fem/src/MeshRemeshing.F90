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


SUBROUTINE Get_MMG3D_Mesh(NewMesh, Parallel)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  LOGICAL :: Parallel
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
  NewMesh => AllocateMesh( NTetras, NTris, NVerts)

  NewMesh % Name = "MMG3D_Output"
  NewMesh%MaxElementNodes=4
  NewMesh%MaxElementDOFs=4
  NewMesh%MeshDim=3

  IF(NPrisms /= 0) CALL Fatal("MMG3D", "Programming Error: MMG3D returns prisms")
  IF(NQuads /= 0) CALL Fatal("MMG3D", "Programming Error: MMG3D returns quads")
  NewMesh % NumberOfNodes = NVerts
  NewMesh % NumberOfBulkElements = NTetras
  NewMesh % NumberOfBoundaryElements = NTris

  CALL AllocateVector( NewMesh % Nodes % x, NewMesh % NumberOfNodes)
  CALL AllocateVector( NewMesh % Nodes % y, NewMesh % NumberOfNodes)
  CALL AllocateVector( NewMesh % Nodes % z, NewMesh % NumberOfNodes)
  CALL AllocateVector( NewMesh % Elements, NewMesh % NumberOfBulkElements + &
       NewMesh % NumberOfBoundaryElements )

  IF(Parallel) THEN
    ALLOCATE(NewMesh % ParallelInfo % GlobalDOFs(NVerts))
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
  End do

  IF (DEBUG) PRINT *,'MMG3D_Get_vertex DONE'

  !! GET NEW TETRAS
  DO ii=1,NewMesh % NumberOfBulkElements
    Element => NewMesh % Elements(ii)
    Element % TYPE => GetElementType( 504 )
    Element % NDOFs = Element % TYPE % NumberOfNodes
    Element % ElementIndex = ii !TODO
    Element % PartIndex = ParEnv % myPE !TODO
    CALL AllocateVector(Element % NodeIndexes, 4)
    NodeIndexes => Element % NodeIndexes
    CALL MMG3D_Get_tetrahedron(mmgMesh, &
         NodeIndexes(1), &
         NodeIndexes(2), &
         NodeIndexes(3), &
         NodeIndexes(4), &
         Element % BodyId, &
         required,ierr)
  END DO
  IF (DEBUG) PRINT *,'MMG3D_Get_tets DONE'

  !! Get BC Elements
  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1

    Element => NewMesh % Elements(kk)

    Element % TYPE => GetElementType( 303 )
    Element % NDOFs = Element % TYPE % NumberOfNodes
    Element % ElementIndex = kk !TODO 
    Element % PartIndex = ParEnv % myPE !TODOx

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
!Local numbering should already be handled by mesh gluing. We assume
!here that each partition may have any number of nodes (inc. zero)
SUBROUTINE RenumberGDOFs(OldMesh,NewMesh)
  TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
  !----------------------
  INTEGER :: i,j,k, OldNN, NewNN
  INTEGER, ALLOCATABLE :: work_int(:)
  LOGICAL, ALLOCATABLE :: NeedGDOF(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="RenumberNodes"

  OldNN = OldMesh % NumberOfNodes
  NewNN = NewMesh % NumberOfNodes

  !Serially check the consistency of the old mesh
  ALLOCATE(work_int(OldNN), NeedGDOF(NewNN))
  work_int = OldMesh % ParallelInfo % GlobalDOFs
  CALL Sort(OldNN,work_int)
  DO i=1,OldNN-1
    IF(work_int(i) == work_int(i+1)) &
         CALL Fatal(FuncName,"OldMesh has duplicate GlobalDOFs")
    IF(work_int(i) == 0) &
         CALL Fatal(FuncName,"OldMesh has at least 1 GlobalDOFs = 0")
  END DO

  !Get a total node count - requires ParallelInfo % NeighbourList


  ! CALL MPI_ALLREDUCE(NewNN,gnodes,1, &
  !      MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

  !Determine which nodes need to get GlobalDOFs
  NeedGDOF = NewMesh % ParallelInfo % GlobalDOFs == 0

  !See if we can fill our new node GlobalDOFs from nodes which no longer exist

END SUBROUTINE RenumberGDOFs

END MODULE MeshRemeshing
