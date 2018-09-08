MODULE MMG3DTools

  USE DefUtils
  IMPLICIT NONE

#include "mmg/mmg3d/libmmg3df.h"

  LOGICAL :: Debug = .TRUE.
  INTEGER :: ier

  MMG5_DATA_PTR_T :: mmgMesh
  MMG5_DATA_PTR_T :: mmgSol

CONTAINS

SUBROUTINE SET_MMG3D_MESH(Mesh)

  USE DefUtils
  IMPLICIT NONE

  TYPE(Mesh_t), POINTER :: Mesh

  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry
  INTEGER, ALLOCATABLE :: NodeRefs(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.

  IF(CoordinateSystemDimension() /= 3) CALL Fatal("MMG3D","Only works for 3D meshes!")

  ALLOCATE(NodeRefs(6))

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
  CALL MMG3D_Set_meshSize(mmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ier)
  IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
       'CALL TO MMG3D_Set_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_meshSize DONE'

  DO i=1,NVerts
    CALL MMG3D_Set_vertex(mmgMesh, Mesh%Nodes%x(i), &
         Mesh%Nodes%y(i),Mesh%Nodes%z(i), 0, i, ier)
!         Mesh%Nodes%y(i),Mesh%Nodes%z(i), 0, Mesh % ParallelInfo % GlobalDOFs(i), ier)
    !PRINT *,'debug: mesh point: ',Mesh%Nodes%x(i), &
    !     Mesh%Nodes%y(i),Mesh%Nodes%z(i), i

    IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
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
      CALL MMG3D_Set_triangle(mmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), 1, ntris,ier)
    CASE(404)
      nquads = nquads + 1
      CALL MMG3D_Set_quadrilateral(mmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  1, nquads,ier)
    CASE(504)
      ntetras = ntetras + 1
      CALL MMG3D_Set_tetrahedron(mmgMesh,NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  2, ntetras,ier)
    CASE(706)
      nprisms = nprisms + 1
      CALL MMG3D_Set_Prism(mmgMesh, NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4), NodeRefs(5), NodeRefs(6), 2, nprisms,ier)
    CASE DEFAULT
    END SELECT
  END DO
  IF (DEBUG) PRINT *,'--**-- MMG3D - Set elements DONE'


END SUBROUTINE SET_MMG3D_MESH

SUBROUTINE SET_MMG3D_PARAMETERS(SolverParams)
  USE DefUtils
  IMPLICIT NONE
  TYPE(ValueList_t), POINTER :: SolverParams

  REAL(KIND=dp) :: Pval
  LOGICAL :: AngleDetect
  INTEGER :: verbosity,MeMIncrease,Bucket,GMshOption
  REAL(KIND=dp) :: Hmin, Hmax, HSiz
  LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
  LOGICAL :: Found

  ! Minimal mesh size:  hmin
  Hmin = GetConstReal( SolverParams, 'hmin', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,&
         Hmin,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hmin> Failed')
  END IF
  ! Maximal mesh size - hmax
  Hmax = GetConstReal( SolverParams, 'hmax', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,&
         Hmax,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hmax> Failed')
  END IF

  hsiz = GetConstReal( SolverParams, 'hsiz', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hsiz,&
         hsiz,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hsiz> Failed')
  END IF

!!! PARAMS: generic options (debug, mem, verbosity)
  ! [val] Set the verbosity level to n
  Verbosity=GetInteger( SolverParams,'verbosity',Found)
  IF (Found) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, &
         Verbosity,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
  ! [val] Set the maximal memory size to n MBytes.
  MemIncrease=GetInteger(SolverParams,'Increase Memory',Found)
  IF (FOUND) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_mem,&
         MemIncrease,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
  ! [0/1] Turn on the debug mode.
  DebugMode=GetLogical(SolverParams,'Debug Mode',Found)
  IF (Found .AND. DebugMode) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_debug,1, &
         ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER Failed')
  END IF
!!! PARAMS
  ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
  ! MMG3D_DPARAM_hausd default est 0.01 semble bien trop petit;
  !  il semble qu'il faille mettre une taille > taille des elements.
  Pval=GetConstReal( SolverParams, 'hausd', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,&
         Pval,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hausd> Failed')
  END IF
  ! Control gradation
  Pval=GetConstReal( SolverParams, 'hgrad', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,&
         Pval,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <hgrad> Failed')
  END IF

!!! OTHER PARAMETERS: NOT ALL TESTED
  Pval = GetConstReal( SolverParams, 'Angle detection',Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
         Pval,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_DPARAMETER <Angle detection> Failed')
  ENDIF
  ! !< [1/0], Avoid/allow surface modifications */ 
  AngleDetect=GetLogical(SolverParams,'No Angle detection',Found)
  IF (Found.AND.AngleDetect) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle, &
         1,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER <No Angle detection> Failed')
  END IF
  ! [1/0] Avoid/allow point insertion
  NoInsert=GetLogical(SolverParams,'No insert',Found)
  IF (Found .AND. NoInsert) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noinsert,&
         1,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver', &
         'CALL TO MMG3D_SET_IPARAMETER <No insert> Failed') 
  END IF
  ! [1/0] Avoid/allow edge or face flipping
  NoSwap=GetLogical(SolverParams,'No swap',Found)
  IF (Found .AND. NoSwap) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noswap,&
         1,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No swap>Failed')
  END IF
  ! [1/0] Avoid/allow point relocation
  NoMove=GetLogical(SolverParams,'No move',Found)
  IF (Found .AND. NoMove) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nomove,&
         1,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No move> Failed')
  END IF
  ! [1/0] Avoid/allow surface modifications
  NoSurf=GetLogical(SolverParams,'No surf',Found)
  IF (Found .AND. NoSurf) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nosurf,&
         0,ier)
    IF ( ier == 0 ) CALL FATAL('M/MGSolver',&
         'CALL TO MMG3D_SET_IPARAMETER <No surf> Failed')
  END IF
!!!
END SUBROUTINE SET_MMG3D_PARAMETERS


SUBROUTINE GET_MMG3D_MESH(NewMesh)
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER ::  Element
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ier
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent,ied
  INTEGER :: ii,kk


  !> a) get the size of the mesh: vertices, tetra,prisms, triangles, quads,edges
  CALL MMG3D_Get_meshSize(mmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ier)
  IF ( ier == 0 ) CALL FATAL('MMGSolver',&
       'CALL TO MMGS_Get_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Get_meshSize DONE'    

  ! INITIALISE THE NEW MESH STRUCTURE
  NewMesh => AllocateMesh()
  ! IF (IncrementMeshNumber) THEN
  !   write(NewMesh % Name,'(A,A,I0)') TRIM(OutPutFileName),'_N',MeshNumber
  ! ELSE
  !   NewMesh % Name=TRIM(OutPutFileName)
  ! END IF
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

  !! GET NEW VERTICES
  NewMesh % Nodes % z = 0._dp
  Do ii=1,NVerts
    CALL MMG3D_Get_vertex(mmgMesh,&
         NewMesh % Nodes % x(ii),&
         NewMesh % Nodes % y(ii),&
         NewMesh % Nodes % z(ii),&
         ref,corner,required,ier)
    IF ( ier == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_vertex FAILED')
  End do
  IF (DEBUG) PRINT *,'MMG3D_Get_vertex DONE'

  !! GET NEW TRIANGLES
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
         required,ier)
  END DO
  IF (DEBUG) PRINT *,'MMG3D_Get_tets DONE'

  !! Get BC Elements
  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTetras
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
         required,ier)
    IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_triangle FAILED')

    Allocate(Element % BoundaryInfo)
    Element % BoundaryInfo % Constraint=ref

  END DO

  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1

    CALL MMG3D_GET_TetFromTria(mmgMesh, &
         ii,parent,ied,ier)

    IF ( ier /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_TetFromTria FAILED')
    Element => NewMesh % Elements(kk)
    Element % BoundaryInfo % Left => NewMesh % Elements(parent) !TODO - parent ID offset?
  END DO

END SUBROUTINE GET_MMG3D_MESH

END MODULE MMG3DTools
