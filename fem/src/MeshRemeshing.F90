!/*****************************************************************************
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
! ****************************************************************************/
!
! ******************************************************************************
! *
! *  Authors: Joe Todd
! *
! ****************************************************************************/

MODULE MeshRemeshing

USE DefUtils
USE Types
USE MeshUtils
USE MeshPartition
USE SparIterComm

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
INTEGER :: MMGPARAM_angleDetection = MMG3D_DPARAM_angleDetection
INTEGER :: MMGPARAM_debug = MMG3D_IPARAM_debug
INTEGER :: MMGPARAM_rmc = MMG3D_DPARAM_rmc
INTEGER :: MMGPARAM_nosurf = MMG3D_IPARAM_nosurf
INTEGER :: MMGPARAM_aniso = MMG3D_IPARAM_anisosize
MMG5_DATA_PTR_T :: mmgMesh
MMG5_DATA_PTR_T :: mmgSol
MMG5_DATA_PTR_T :: mmgMet
MMG5_DATA_PTR_T :: mmgLs
#endif

#ifdef HAVE_PARMMG
#include "parmmg/libparmmgtypesf.h"
#endif

#ifdef HAVE_PARMMG
!INTEGER :: PMMGPARAM_hausd = PMMG_DPARAM_hausd
!INTEGER :: PMMGPARAM_hmin = PMMG_DPARAM_hmin
!INTEGER :: PMMGPARAM_hmax = PMMG_DPARAM_hmax
!INTEGER :: PMMGPARAM_iso = PMMG_IPARAM_iso
!INTEGER :: PMMGPARAM_hgrad = PMMG_DPARAM_hgrad
!INTEGER :: PMMGPARAM_angle = PMMG_IPARAM_angle
!INTEGER :: PMMGPARAM_angleDetection = PMMG_DPARAM_angleDetection
!INTEGER :: PMMGPARAM_debug = PMMG_IPARAM_debug
!INTEGER :: PMMGPARAM_nosurf = PMMG_IPARAM_nosurf
!INTEGER :: PMMGPARAM_aniso = PMMG_IPARAM_anisosize
MMG5_DATA_PTR_T :: pmmgMesh
MMG5_DATA_PTR_T :: pmmgMet
#endif

CONTAINS

!============================================
!============================================
!            MMG3D SUBROUTINES
!============================================
!============================================


SUBROUTINE Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount

#ifdef HAVE_MMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ref,ierr
  INTEGER, ALLOCATABLE :: NodeRefs(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.FALSE.,Elem202
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
  IF(Present(PairCount)) NEdges= PairCount

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
    ! ref = GDOF + 10 to avoid an input ref of 10 being confused
    ! with mmg output ref of 10 which occurs on some new nodes
    ! GDOF = ref - 10
    IF(Parallel) ref = Mesh % ParallelInfo % GlobalDOFs(i) + 10
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

  !! use element pairs '202' elements
  Elem202 = (PRESENT(EdgePairs))
  IF (Elem202) THEN
    DO i=1, PairCount
      NEdges = NEdges + 1 
      CALL MMG3D_Set_edge(mmgMesh, EdgePairs(1,i), EdgePairs(2,i), 1, nedges, ierr)
      CALL MMG3D_Set_ridge(mmgMesh, nedges, ierr)
      !CALL MMG3D_Set_requiredEdge(mmgMesh, Nedges, ierr)
    END DO
  END IF

  

  IF (DEBUG) PRINT *, '--**-- MMG3D - Set edge elements DONE'

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
  INTEGER, ALLOCATABLE :: BC2BodyMap(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent,ied
  INTEGER :: i,ii,kk,NoBCs
  LOGICAL :: Found, Debug

  !Set up a map of BoundaryInfo % Constraint to % BodyID
  NoBCs = CurrentModel % NumberOfBCs
  ALLOCATE(BC2BodyMap(NoBCs))
  BC2BodyMap = 0
  DO i=1,NoBCs
    BC2BodyMap(i) = ListGetInteger( &
           CurrentModel % BCs(i) % Values, 'Body Id', Found)
    IF(.NOT. Found) BC2BodyMap(i) = 0
  END DO

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
        ! ref = GDOF + 10 to avoid an input ref of 10 being confused
        ! with mmg output ref of 10 which occurs on some new nodes
        ! GDOF = ref - 10
        NewMesh % ParallelInfo % GlobalDOFs(ii) = ref - 10
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

    IF(ref > 0 .AND. ref <= CurrentModel % NumberOfBCs) THEN
      Element % BodyId = BC2BodyMap(ref)
    END IF

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
SUBROUTINE RemeshMMG3D(Model, InMesh,OutMesh,EdgePairs,PairCount,NodeFixed,ElemFixed,Params,Success)

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  TYPE(ValueList_t), POINTER, OPTIONAL :: Params
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount
  LOGICAL :: Success
  !-----------
  TYPE(Mesh_t), POINTER :: WorkMesh
  TYPE(ValueList_t), POINTER :: FuncParams, Material
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:),hminarray(:),hausdarray(:)
  REAL(KIND=dp), POINTER :: WorkReal(:,:,:) => NULL(), WorkArray(:,:) => NULL()
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName, MeshName, SolName, &
       premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile

  SAVE :: WorkReal,WorkArray

#ifdef HAVE_MMG

  Debug = .TRUE.
  Parallel = ParEnv % PEs > 1
  FuncName = "RemeshMMG3D"

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  !Optionally pass valuelist, by default use the Simulation section
  IF(PRESENT(Params)) THEN
    FuncParams => Params
  ELSE
    FuncParams => GetMaterial(InMesh % Elements(1)) !TODO, this is not generalised
  END IF

  !Get parameters from valuelist
  !hausd, hmin, hmax, hgrad, anisoflag, the metric
  !Scalar, vector, tensor metric?

  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hmin", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hmin"')
  MaxRemeshIter= SIZE(WorkArray(:,1))
  ALLOCATE(hminarray(MaxRemeshIter))
  hminarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  hmax = ListGetConstReal(FuncParams, "RemeshMMG3D Hmax",  Default=4000.0_dp)
  hgrad = ListGetConstReal(FuncParams,"RemeshMMG3D Hgrad", Default=0.5_dp)
  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hausd", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hausd"')
  IF(MaxRemeshIter /= SIZE(WorkArray(:,1))) CALL FATAL(FuncName, 'The number of hmin options &
          must equal the number of hausd options')
  ALLOCATE(hausdarray(MaxRemeshIter))
  hausdarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  RemeshMinQuality = ListGetConstReal(FuncParams, "RemeshMMG3D Min Quality",Default=0.0001_dp)
  AnisoFlag = ListGetLogical(FuncParams, "RemeshMMG3D Anisotropic", Default=.TRUE.)
  SaveMMGMeshes = ListGetLogical(FuncParams,"Save RemeshMMG3D Meshes", Default=.FALSE.)
  SaveMMGSols = ListGetLogical(FuncParams,"Save RemeshMMG3D Sols", Default=.FALSE.)
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre RemeshMMG3D Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "RemeshMMG3D Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre RemeshMMG3D Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "RemeshMMG3D Output Sol Name", UnfoundFatal = .TRUE.)
  END IF

  NNodes = InMesh % NumberOfNodes
  NBulk = InMesh % NumberOfBulkElements
  NBdry = InMesh % NumberOfBoundaryElements

  IF(AnisoFlag) THEN

    WorkMesh => Model % Mesh
    Model % Mesh => InMesh

    SolType = MMG5_Tensor
    !Upper triangle of symmetric tensor: 11,12,13,22,23,33
    ALLOCATE(Metric(NNodes,6))
    Metric = 0.0
    DO i=1,NNodes
      NodeNum = i

      CALL ListGetRealArray(FuncParams,"RemeshMMG3D Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)

      !Metric = 1.0/(edge_length**2)
      Metric(i,1) = 1.0 / (WorkReal(1,1,1)**2.0)
      Metric(i,4) = 1.0 / (WorkReal(2,1,1)**2.0)
      Metric(i,6) = 1.0 / (WorkReal(3,1,1)**2.0)

    END DO

    Model % Mesh => WorkMesh
    WorkMesh => NULL()

  ELSE

    SolType = MMG5_Scalar
    ALLOCATE(Metric(NNodes, 1))
    DO i=1,NNodes
      NodeNum = i
      Metric(i,:) = ListGetReal(FuncParams,"RemeshMMG3D Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
    END DO

  END IF

  body_offset = CurrentModel % NumberOfBCs + CurrentModel % NumberOfBodies + 1
  nBCs = CurrentModel % NumberOfBCs
  DO i=1,InMesh % NumberOfBulkElements
    InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID + body_offset
  END DO

  mmgloops=0
  Success=.TRUE.

10 CONTINUE

  mmgloops = mmgloops+1
  hmin = hminarray(mmgloops)
  Hausd = hausdarray(mmgloops)

  WRITE(Message, '(A,F10.5,A,F10.5)') 'Applying levelset with Hmin ',Hmin, ' and Hausd ', Hausd
  CALL INFO(FuncName, Message)

  mmgMesh = 0
  mmgSol  = 0

  !---------------------------------
  ! Issue here: MMG3D will 'helpfully' add any
  ! missing boundary triangles, assigning them
  ! % BoundaryInfo % Constraint = Parent % BodyID
  !
  ! To deal with this, temporarily offset BodyID of
  ! all body elems by Model % NumberOfBodies +
  ! Model % NumberOfBCs + 1, then afterwards we
  ! delete all the extra BC elems & revert the bodyID
  !----------------------------------

  CALL MMG3D_Init_mesh(MMG5_ARG_start, &
       MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
       MMG5_ARG_end)

  IF (Present(PairCount)) THEN
    CALL SET_MMG3D_MESH(InMesh,Parallel,EdgePairs,PairCount)
  ELSE
    CALL SET_MMG3D_MESH(InMesh,Parallel)
  END IF

  !Set the metric values at nodes
  CALL MMG3D_Set_SolSize(mmgMesh, mmgSol, MMG5_Vertex, NNodes, SolType,ierr)
  IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")

  DO i=1,NNodes
    IF(AnisoFlag) THEN
      IF(Debug) PRINT *,'debug sol at ',i,' is: ',Metric(i,:)
      CALL MMG3D_Set_TensorSol(mmgSol,Metric(i,1),Metric(i,2),Metric(i,3),&
           Metric(i,4),Metric(i,5),Metric(i,6),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set tensor solution at vertex")
    ELSE
      CALL MMG3D_Set_ScalarSol(mmgSol,Metric(i,1),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
    END IF
  END DO

  !Turn on debug (1)
  CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMGPARAM_debug, &
  1,ierr)
   
  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmin,&
       hmin,ierr)
  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmax,&
       hmax,ierr)
  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,&
       hausd,ierr)
  CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,&
       hgrad,ierr)
  ! allow surface modifications
  CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMGPARAM_nosurf,&
       0,ierr)

  !Turn off sharp angle detection (0)   
  CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle, &
       0,ierr)
  !Option to set angle detection threshold:
  !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
  !      85.0_dp,ierr)


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

  IF(SaveMMGMeshes) THEN
    WRITE(MeshName, '(A,i0,A)') TRIM(premmg_meshfile), time, '.mesh'
    CALL MMG3D_SaveMesh(mmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
  END IF
  IF(SaveMMGSols) THEN
    WRITE(SolName, '(A,i0,A)') TRIM(premmg_solfile), time, '.sol'
    CALL MMG3D_SaveSol(mmgMesh, mmgSol,SolName,LEN(TRIM(SolName)),ierr)
  END IF
  IF (DEBUG) PRINT *,'--**-- SET MMG3D PARAMETERS '
  ! CALL SET_MMG3D_PARAMETERS(SolverParams)

  CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ierr)
  IF ( ierr == MMG5_STRONGFAILURE .OR. ierr == MMG5_LOWFAILURE ) THEN
    PRINT*,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH"
    !! Release mmg mesh
    CALL MMG3D_Free_all(MMG5_ARG_start, &
        MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
        MMG5_ARG_end)
    WRITE(Message, '(A,F10.5,A,F10.5)') 'Remesh failed with Hmin ',Hmin, ' and Hausd ', Hausd
    CALL WARN(FuncName, Message)
    IF(mmgloops==MaxRemeshIter) THEN
      Success=.FALSE.
      GO TO 20
    ELSE
      GO TO 10
    END IF
    !STOP MMG5_STRONGFAILURE
  ENDIF

  IF (DEBUG) PRINT *,'--**-- MMG3D_mmg3dlib DONE'

  IF(SaveMMGMeshes) THEN
    WRITE(MeshName, '(A,i0,A)') TRIM(mmg_meshfile), time, '.mesh'
    CALL MMG3D_SaveMesh(mmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
  END IF
  IF(SaveMMGSols) THEN
    WRITE(SolName, '(A,i0,A)') TRIM(mmg_solfile), time, '.sol'
    CALL MMG3D_SaveSol(mmgMesh, mmgSol,SolName,LEN(TRIM(SolName)),ierr)
  END IF

  CALL MMG3D_Get_meshSize(mmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

  counter=0
  Do i=1, NTetras
    CALL MMG3D_Get_TetrahedronQuality(mmgMesh, mmgSol, i, Quality)
    IF(Quality == 0) CALL WARN(FuncName, 'Remeshing could not determine elem quality')
    IF(Quality <= RemeshMinQuality) counter = counter+1
  END DO
  IF(Debug) PRINT*, 'Bad Element Count: ', Counter
  IF ( Counter > 0 ) THEN
    PRINT*,"Bad elements detected - reruning remeshing"
    !! Release mmg mesh
    CALL MMG3D_Free_all(MMG5_ARG_start, &
        MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
        MMG5_ARG_end)
    WRITE(Message, '(A,F10.5,A,F10.5)') 'Remesh failed with Hmin ',Hmin, ' and Hausd ', Hausd
    CALL WARN(FuncName, Message)
    IF(mmgloops==MaxRemeshIter) THEN
      Success=.FALSE.
      PRINT*, 'remesh failed on bad elements'
      GO TO 20
    ELSE
      GO TO 10
    END IF
    !STOP MMG5_STRONGFAILURE
  ENDIF

  !! GET THE NEW MESH
  CALL GET_MMG3D_MESH(OutMesh,Parallel)

  !! Release mmg mesh
  CALL MMG3D_Free_all(MMG5_ARG_start, &
  MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
  MMG5_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements

  !Reset the BodyIDs (see above)
  DO i=1,NBulk
    OutMesh % Elements(i) % BodyID = OutMesh % Elements(i) % BodyID - body_offset
  END DO

  !And delete the unneeded BC elems
  ALLOCATE(RmElement(NBulk+NBdry))
  RmElement = .FALSE.
  DO i=NBulk+1,NBulk+NBdry
    Element => OutMesh % Elements(i)
    IF(Element % BoundaryInfo % Constraint > nBCs) THEN
      RmElement(i) = .TRUE.
    END IF
  END DO
  CALL CutMesh(OutMesh, RmElem=RmElement)

20 CONTINUE

  ! if remeshing has failed need to reset body ids
  IF(.NOT. Success) THEN
    DO i=1,InMesh % NumberOfBulkElements
      InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID - body_offset
    END DO
  END IF

#else
  CALL Fatal(FuncName, "Remeshing utility MMG3D has not been installed")
#endif

END SUBROUTINE RemeshMMG3D

!A subroutine for a 3D mesh (in parallel) with ParMMG3D
!Inputs:
!   InMesh - the initial mesh
!   Metric - 2D real array specifying target metric
!   NodeFixed, ElemFixed - Optional mask to specify 'required' entities
!Output:
!   OutMesh - the improved mesh
!
SUBROUTINE SequentialRemeshParMMG(Model, InMesh,OutMesh,Boss,EdgePairs,PairCount,NodeFixed,ElemFixed,Params)

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Boss
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount
  LOGICAL :: Success
  !-----------
  TYPE(Mesh_t), POINTER :: WorkMesh
  TYPE(ValueList_t), POINTER :: FuncParams, Material
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:),hminarray(:),hausdarray(:)
  REAL(KIND=dp), POINTER :: WorkReal(:,:,:) => NULL(), WorkArray(:,:) => NULL()
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName, MeshName, SolName, &
       premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile

  SAVE :: WorkReal,WorkArray

#ifdef HAVE_PARMMG

  Debug = .TRUE.
  Parallel = ParEnv % PEs > 1
  FuncName = "RemeshSeqParMMG3D"

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  ! params must be provided as not a global mesh
  FuncParams => Params

  !Get parameters from valuelist
  !hausd, hmin, hmax, hgrad, anisoflag, the metric
  !Scalar, vector, tensor metric?

  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hmin", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hmin"')
  MaxRemeshIter= SIZE(WorkArray(:,1))
  ALLOCATE(hminarray(MaxRemeshIter))
  hminarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  hmax = ListGetConstReal(FuncParams, "RemeshMMG3D Hmax",  Default=4000.0_dp)
  hgrad = ListGetConstReal(FuncParams,"RemeshMMG3D Hgrad", Default=0.5_dp)
  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hausd", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hausd"')
  IF(MaxRemeshIter /= SIZE(WorkArray(:,1))) CALL FATAL(FuncName, 'The number of hmin options &
            must equal the number of hausd options')
  ALLOCATE(hausdarray(MaxRemeshIter))
  hausdarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  RemeshMinQuality = ListGetConstReal(FuncParams, "RemeshMMG3D Min Quality",Default=0.0001_dp)
  AnisoFlag = ListGetLogical(FuncParams, "RemeshMMG3D Anisotropic", Default=.TRUE.)
  SaveMMGMeshes = ListGetLogical(FuncParams,"Save RemeshMMG3D Meshes", Default=.FALSE.)
  SaveMMGSols = ListGetLogical(FuncParams,"Save RemeshMMG3D Sols", Default=.FALSE.)
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre RemeshMMG3D Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "RemeshMMG3D Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre RemeshMMG3D Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "RemeshMMG3D Output Sol Name", UnfoundFatal = .TRUE.)
  END IF


  IF(Boss) THEN
    NNodes = InMesh % NumberOfNodes
    NBulk = InMesh % NumberOfBulkElements
    NBdry = InMesh % NumberOfBoundaryElements

    IF(AnisoFlag) THEN

      WorkMesh => Model % Mesh
      Model % Mesh => InMesh

      SolType = MMG5_Tensor
      !Upper triangle of symmetric tensor: 11,12,13,22,23,33
      ALLOCATE(Metric(NNodes,6))
      Metric = 0.0
      DO i=1,NNodes
        NodeNum = i

        CALL ListGetRealArray(FuncParams,"RemeshMMG3D Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)

        !Metric = 1.0/(edge_length**2)
        Metric(i,1) = 1.0 / (WorkReal(1,1,1)**2.0)
        Metric(i,4) = 1.0 / (WorkReal(2,1,1)**2.0)
        Metric(i,6) = 1.0 / (WorkReal(3,1,1)**2.0)

      END DO

      Model % Mesh => WorkMesh
      WorkMesh => NULL()

    ELSE

      SolType = MMG5_Scalar
      ALLOCATE(Metric(NNodes, 1))
      DO i=1,NNodes
        NodeNum = i
        Metric(i,:) = ListGetReal(FuncParams,"RemeshMMG3D Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
      END DO

    END IF

    body_offset = CurrentModel % NumberOfBCs + CurrentModel % NumberOfBodies + 1
    nBCs = CurrentModel % NumberOfBCs
    DO i=1,InMesh % NumberOfBulkElements
      InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID + body_offset
    END DO


    mmgloops=0
    Success=.TRUE.

10 CONTINUE

    mmgloops = mmgloops+1
    hmin = hminarray(mmgloops)
    Hausd = hausdarray(mmgloops)

    WRITE(Message, '(A,F10.5,A,F10.5)') 'Applying levelset with Hmin ',Hmin, ' and Hausd ', Hausd
    CALL INFO(FuncName, Message)

  END IF

  pmmgMesh = 0
  pmmgMet  = 0

  !---------------------------------
  ! Issue here: MMG3D will 'helpfully' add any
  ! missing boundary triangles, assigning them
  ! % BoundaryInfo % Constraint = Parent % BodyID
  !
  ! To deal with this, temporarily offset BodyID of
  ! all body elems by Model % NumberOfBodies +
  ! Model % NumberOfBCs + 1, then afterwards we
  ! delete all the extra BC elems & revert the bodyID
  !----------------------------------
  PRINT*, 'Initiating parmmg mesh', 'My PE', ParEnv % MyPE

  CALL PMMG_Init_parMesh(PMMG_ARG_start, &
      PMMG_ARG_ppParMesh,pmmgMesh, PMMG_ARG_pMesh,PMMG_ARG_pMet, &
      PMMG_ARG_dim,%val(3),PMMG_ARG_MPIComm,%val(ELMER_COMM_WORLD), &
      PMMG_ARG_end)

  IF(Boss) THEN
    PRINT*, 'Setting mesh ....'
    IF (Present(PairCount)) THEN
      CALL SET_ParMMG_MESH(InMesh,Parallel,EdgePairs,PairCount)
    ELSE
      CALL SET_ParMMG_MESH(InMesh,Parallel)
    END IF
    PRINT*, 'Finished setting mesh ....', pmmgMesh

    PRINT*, 'Setting met size', pmmgMesh
    !Set the metric values at nodes
    CALL PMMG_Set_MetSize(pmmgMesh, MMG5_Vertex, NNodes, SolType,ierr)
    IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")

    PRINT*, 'Set met size', pmmgMesh

    DO i=1,NNodes
      IF(AnisoFlag) THEN
        IF(Debug) PRINT *,'debug sol at ',i,' is: ',Metric(i,:)
        CALL PMMG_Set_TensorMet(pmmgMesh,Metric(i,1),Metric(i,2),Metric(i,3),&
            Metric(i,4),Metric(i,5),Metric(i,6),i,ierr)
        IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set tensor solution at vertex")
      ELSE
        CALL PMMG_Set_ScalarMet(pmmgMesh,Metric(i,1),i,ierr)
        IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
      END IF
    END DO
  END IF

  !Turn on debug (1)
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_debug, &
       1,ierr)
   
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmin,&
       hmin,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmax,&
       hmax,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hausd,&
       hausd,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hgrad,&
       hgrad,ierr)
  ! allow surface modifications
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_nosurf,&
       0,ierr)

  !Turn off sharp angle detection (0)   
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle, &
       0,ierr)
  !Option to set angle detection threshold:
  !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
  !      85.0_dp,ierr)

  IF(Boss) THEN
    !Take care of fixed nodes/elements if requested
    IF(PRESENT(NodeFixed)) THEN
      DO i=1,NNodes
        IF(NodeFixed(i)) THEN
          CALL PMMG_SET_REQUIREDVERTEX(pmmgMesh,i,ierr)
        END IF
      END DO
    END IF

    IF(PRESENT(ElemFixed)) THEN
      DO i=1,NBulk + NBdry
        IF(ElemFixed(i)) THEN
          IF(i <= NBulk) THEN
            CALL PMMG_SET_REQUIREDTETRAHEDRON(pmmgMesh,i,ierr)
          ELSE
            CALL PMMG_SET_REQUIREDTRIANGLE(pmmgMesh,i-NBulk,ierr)
          END IF
        END IF
      END DO
    END IF

    IF(SaveMMGMeshes) THEN
      WRITE(MeshName, '(A,i0,A)') TRIM(premmg_meshfile), time, '.mesh'
      CALL PMMG_SaveMesh_Centralized(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
    END IF
    IF(SaveMMGSols) THEN
      WRITE(SolName, '(A,i0,A)') TRIM(premmg_solfile), time, '.sol'
      CALL PMMG_SaveMet_Centralized(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
    END IF
    IF (DEBUG) PRINT *,'--**-- SET PMMG3D PARAMETERS '
    ! CALL SET_MMG3D_PARAMETERS(SolverParams)
  END IF

  CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
  CALL PMMG_parmmglib_centralized(pmmgMesh,ierr)
  IF ( ierr == PMMG_STRONGFAILURE .OR. ierr == PMMG_LOWFAILURE ) THEN
    PRINT*,"BAD ENDING OF PMMGLIB: UNABLE TO SAVE MESH"
    !! Release mmg mesh
    CALL PMMG_Free_all ( PMMG_ARG_start,     &
        PMMG_ARG_ppParMesh,pmmgMesh,         &
        PMMG_ARG_end)
    WRITE(Message, '(A,F10.5,A,F10.5)') 'Remesh failed with Hmin ',Hmin, ' and Hausd ', Hausd
    CALL WARN(FuncName, Message)
    IF(mmgloops==MaxRemeshIter) THEN
      Success=.FALSE.
      GO TO 20
    ELSE
      GO TO 10
    END IF
    !STOP MMG5_STRONGFAILURE
  ENDIF

  IF (DEBUG) PRINT *,'--**-- PMMG_parmmglib_centralized DONE'

  IF(Boss) THEN
    IF(SaveMMGMeshes) THEN
      WRITE(MeshName, '(A,i0,A)') TRIM(mmg_meshfile), time, '.mesh'
      CALL PMMG_SaveMesh_Centralized(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
    END IF
    IF(SaveMMGSols) THEN
      WRITE(SolName, '(A,i0,A)') TRIM(mmg_solfile), time, '.sol'
      CALL PMMG_SaveMet_Centralized(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
    END IF

    CALL PMMG_Get_meshSize(pmmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

    !! GET THE NEW MESH
    PRINT*, 'Recovering mesh'
    CALL GET_ParMMG_MESH(OutMesh,Parallel)
  END IF

  PRINT*, 'Releasing mesh'
  !! Release mmg mesh
  CALL PMMG_Free_all ( PMMG_ARG_start,     &
       PMMG_ARG_ppParMesh,pmmgMesh,         &
       PMMG_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements

  IF(Boss) THEN
    !Reset the BodyIDs (see above)
    DO i=1,NBulk
      OutMesh % Elements(i) % BodyID = OutMesh % Elements(i) % BodyID - body_offset
    END DO

    !And delete the unneeded BC elems
    ALLOCATE(RmElement(NBulk+NBdry))
    RmElement = .FALSE.
    DO i=NBulk+1,NBulk+NBdry
      Element => OutMesh % Elements(i)
      IF(Element % BoundaryInfo % Constraint > nBCs) THEN
        RmElement(i) = .TRUE.
      END IF
    END DO
    CALL CutMesh(OutMesh, RmElem=RmElement)
  END IF

20 CONTINUE

  IF(Boss) THEN
    ! if remeshing has failed need to reset body ids
    IF(.NOT. Success) THEN
      DO i=1,InMesh % NumberOfBulkElements
        InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID - body_offset
      END DO
    END IF
  END IF

#else
  CALL Fatal(FuncName, "Remeshing utility PMMG has not been installed")
#endif

END SUBROUTINE SequentialRemeshParMMG

SUBROUTINE Set_ParMMG_Mesh(Mesh, Parallel, EdgePairs, PairCount)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount

#ifdef HAVE_PARMMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,j,k, NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ref,ierr,&
      NoNeighbours, counter
  INTEGER, ALLOCATABLE :: NodeRefs(:), NeighbourList(:), NSharedNodes(:), SharedNodes(:,:),&
      SharedNodesGlobal(:,:)
  LOGICAL, ALLOCATABLE :: IsNeighbour(:)
  INTEGER, POINTER :: Neighbours(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.TRUE.,Elem202
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Set_ParMMG_Mesh"
  IF(CoordinateSystemDimension() /= 3) CALL Fatal("ParMMG","Only works for 3D meshes!")

  ALLOCATE(NodeRefs(6))

  IF(Parallel) CALL Assert(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs), FuncName,&
       "Parallel sim but no ParallelInfo % GlobalDOFs")

  nverts = Mesh % NumberOfNodes
  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0
  IF(Present(PairCount)) NEdges= PairCount

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
  CALL PMMG_Set_meshSize(pmmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ierr)
  PRINT*, 'Set MEsh SIze', ParEnv % MyPE, nverts
  IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
       'CALL TO MMG3D_Set_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_meshSize DONE', ParEnv % Mype

  ref = 0
  DO i=1,NVerts
    ! ref = GDOF + 10 to avoid an input ref of 10 being confused
    ! with mmg output ref of 10 which occurs on some new nodes
    ! GDOF = ref - 10
    IF(Parallel) ref = Mesh % ParallelInfo % GlobalDOFs(i) + 10
    CALL PMMG_Set_vertex(pmmgMesh, Mesh%Nodes%x(i), &
         Mesh%Nodes%y(i),Mesh%Nodes%z(i), ref, i, ierr)
    IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO MMG3D_Set_vertex FAILED')
  END DO
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_vertex DONE', ParEnv % Mype

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

    SELECT CASE(Element % TYPE % ElementCode)
    CASE(303)
      ntris = ntris + 1
      CALL PMMG_Set_triangle(pmmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           Element % BoundaryInfo % Constraint, ntris, ierr)
    CASE(404)
      nquads = nquads + 1
      CALL PMMG_Set_quadrilateral(pmmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  Element % BoundaryInfo % Constraint, nquads,ierr)
    CASE(504)
      ntetras = ntetras + 1
      CALL PMMG_Set_tetrahedron(pmmgMesh,NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4),  Element % BodyID, ntetras,ierr)
    CASE(706)
      nprisms = nprisms + 1
      CALL PMMG_Set_Prism(pmmgMesh, NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           NodeRefs(4), NodeRefs(5), NodeRefs(6), Element % BodyID, nprisms,ierr)
    CASE DEFAULT
    END SELECT
  END DO
  IF (DEBUG) PRINT *,'--**-- MMG3D - Set elements DONE', ParEnv % Mype

  !! use element pairs '202' elements
  Elem202 = (PRESENT(EdgePairs))
  IF (Elem202) THEN
    DO i=1, PairCount
      NEdges = NEdges + 1
      CALL PMMG_Set_edge(pmmgMesh, EdgePairs(1,i), EdgePairs(2,i), 1, nedges, ierr)
      CALL PMMG_Set_ridge(pmmgMesh, nedges, ierr)
      !CALL MMG3D_Set_requiredEdge(mmgMesh, Nedges, ierr)
    END DO
  END IF

  IF (DEBUG) PRINT *, '--**-- MMG3D - Set edge elements DONE', ParEnv % Mype

  ! use nodes to set mpi comms
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_APImode, 1, ierr)

  ! set the number of neighbours
  ALLOCATE(IsNeighbour(ParEnv % PEs))
  IsNeighbour = .FALSE.
  DO i=1, Mesh % NumberOfNodes
    Neighbours => Mesh %  ParallelInfo % NeighbourList(i) % Neighbours
    DO j=1, SIZE(Neighbours)
      IF(Neighbours(j) == ParEnv % MyPE) CYCLE
      IsNeighbour(Neighbours(j)+1) = .TRUE.
    END DO
  END DO

  NoNeighbours = COUNT(IsNeighbour)
  NeighbourList = PACK( (/ (i, i=1, ParEnv % PEs) /), IsNeighbour)

  ALLOCATE(NSharedNodes(NoNeighbours))
  NSharedNodes = 0
  DO i=1, Mesh % NumberOfNodes
    Neighbours => Mesh %  ParallelInfo % NeighbourList(i) % Neighbours
    DO j=1, SIZE(Neighbours)
      DO k=1, NoNeighbours
        IF(Neighbours(j) == NeighbourList(k)-1) &
          NSharedNodes(k) = NSharedNodes(k) + 1
      END DO
    END DO
  END DO

  ALLOCATE(SharedNodes(NoNeighbours, MAXVAL(NSharedNodes)), &
          SharedNodesGlobal(NoNeighbours, MAXVAL(NSharedNodes)))
  SharedNodes = 0
  DO i=1, NoNeighbours
    counter = 0
    DO j=1, Mesh % NumberOfNodes
      Neighbours => Mesh %  ParallelInfo % NeighbourList(j) % Neighbours
      DO k=1, SIZE(Neighbours)
        IF(Neighbours(k) == NeighbourList(i)-1) THEN
          Counter = counter + 1
          SharedNodes(i,counter) = j
          SharedNodesGlobal(i, counter) = Mesh % ParallelInfo % GlobalDOFs(j)
        END IF
      END DO
    END DO
  END DO

  CALL PMMG_SET_NUMBEROFNODECOMMUNICATORS(pmmgMesh, NoNeighbours, ierr)

  ! which node has what neighbours
  DO i=1, NoNeighbours
    CALL PMMG_SET_ITHNODECOMMUNICATORSIZE(pmmgMesh, i-1, NeighbourList(i)-1, &
          NSharedNodes(i), ierr)
    ! 1 is unorder nodes, 0 assumed ordered same on both sides
    ! elmer mesh is usually ordered but not after remeshing
    CALL PMMG_Set_ithNodeCommunicator_nodes(pmmgMesh, i-1, SharedNodes(i,:),&
          SharedNodesGlobal(i,:), 1, ierr)
  END DO

#else
     CALL FATAL('Set_ParMMG_Mesh',&
        'Remeshing utility ParMMG has not been installed')
#endif

END SUBROUTINE Set_ParMMG_Mesh

SUBROUTINE Get_ParMMG_Mesh(NewMesh, Parallel, FixedNodes, FixedElems)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  LOGICAL :: Parallel
  LOGICAL, OPTIONAL, ALLOCATABLE :: FixedNodes(:), FixedElems(:)
  !------------------------------------------------------------------------------

#ifdef HAVE_PARMMG

  TYPE(Element_t),POINTER ::  Element
  INTEGER, POINTER :: NodeIndexes(:), Neighbours(:), NSharedNodes(:), SharedNodes(:,:),&
          GlobalNums(:), GlobalNums2(:)
  INTEGER, ALLOCATABLE :: BC2BodyMap(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr, NoNeighbours,&
          OutProc, counter, GlobalID, owner, unique, ntot
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent,ied
  INTEGER :: i,j,ii,kk,NoBCs
  LOGICAL :: Found, Debug

  Debug= .TRUE.

  !Set up a map of BoundaryInfo % Constraint to % BodyID
  NoBCs = CurrentModel % NumberOfBCs
  ALLOCATE(BC2BodyMap(NoBCs))
  BC2BodyMap = 0
  DO i=1,NoBCs
    BC2BodyMap(i) = ListGetInteger( &
           CurrentModel % BCs(i) % Values, 'Body Id', Found)
    IF(.NOT. Found) BC2BodyMap(i) = 0
  END DO

  !> a) get the size of the mesh: vertices, tetra,prisms, triangles, quads,edges
  CALL PMMG_Get_meshSize(pmmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)
  IF ( ierr == 0 ) CALL FATAL('ParMMGSolver',&
       'CALL TO PMMG_Get_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- PMMG_Get_meshSize DONE'    

  ! INITIALISE THE NEW MESH STRUCTURE
  !NPrisms, NQuads should be zero
  NewMesh => AllocateMesh( NTetras, NTris, NVerts, .TRUE.)

  NewMesh % Name = "ParMMG_Output"
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

  IF(NPrisms /= 0) CALL Fatal("ParMMG", "Programming Error: ParMMG returns prisms")
  IF(NQuads /= 0) CALL Fatal("ParMMG", "Programming Error: ParMMG returns quads")
  NewMesh % NumberOfNodes = NVerts
  NewMesh % NumberOfBulkElements = NTetras
  NewMesh % NumberOfBoundaryElements = NTris

  IF(Parallel) THEN
    NewMesh % ParallelInfo % Interface = .FALSE.
  END IF

  !! GET NEW VERTICES
  NewMesh % Nodes % z = 0._dp
  Do ii=1,NVerts
    CALL PMMG_Get_vertex(pmmgMesh,&
         NewMesh % Nodes % x(ii),&
         NewMesh % Nodes % y(ii),&
         NewMesh % Nodes % z(ii),&
         ref,corner,required,ierr)
    IF ( ierr == 0 ) CALL FATAL('ParMMGSolver',&
         'CALL TO  ParMMG_Get_vertex FAILED')
    IF(PRESENT(FixedNodes)) FixedNodes(ii) = required > 0 .AND. ref > 10
  End do

  IF (DEBUG) PRINT *,'ParMMG_Get_vertex DONE'

  !! GET NEW TETRAS
  DO ii=1,NewMesh % NumberOfBulkElements
    Element => NewMesh % Elements(ii)
    Element % TYPE => GetElementType( 504 )
    Element % NDOFs = Element % TYPE % NumberOfNodes
    Element % ElementIndex = ii
    Element % PartIndex = ParEnv % MyPE
    CALL AllocateVector(Element % NodeIndexes, 4)
    NodeIndexes => Element % NodeIndexes
    CALL PMMG_Get_tetrahedron(pmmgMesh, &
         NodeIndexes(1), &
         NodeIndexes(2), &
         NodeIndexes(3), &
         NodeIndexes(4), &
         Element % BodyId, & !TODO - many tetras end up with very high BodyIDs
         required,ierr)
    IF(PRESENT(FixedElems)) FixedElems(ii) = required > 0
  END DO
  IF (DEBUG) PRINT *,'ParMMG_Get_tets DONE'

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
    CALL PMMG_Get_triangle(pmmgMesh, &
         NodeIndexes(1), &
         NodeIndexes(2), &
         NodeIndexes(3), &
         ref   , &
         required,ierr)
    IF ( ierr /= 1 ) CALL FATAL('ParMMGSolver',&
         'CALL TO  PMMG_Get_triangle FAILED')

    Allocate(Element % BoundaryInfo)
    Element % BoundaryInfo % Constraint=ref

    IF(ref > 0 .AND. ref <= CurrentModel % NumberOfBCs) THEN
      Element % BodyId = BC2BodyMap(ref)
    END IF

    IF(PRESENT(FixedElems)) FixedElems(kk) = required > 0
  END DO

  IF (DEBUG) PRINT *,'ParMMG_Get_triangle DONE'

  ! currently segfaults as no parmmg routine to recover parent id
  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1

    CALL PMMG_GET_TetFromTria(pmmgMesh, &
         ii,parent,ied,ierr)

    IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_TetFromTria FAILED')
    Element => NewMesh % Elements(kk)
    Element % BoundaryInfo % Left => NewMesh % Elements(parent) !TODO - parent ID offset?
  END DO

  ! get parallel info back
  ! get number of neighbours
  CALL PMMG_GET_NUMBEROFNODECOMMUNICATORS(pmmgMesh,NoNeighbours,ierr)

  ALLOCATE(Neighbours(NoNeighbours), NSharedNodes(NoNeighbours))
  DO i=1, NoNeighbours
    OutProc = i-1
    CALL PMMG_Get_ithNodeCommunicatorSize(pmmgMesh, outProc, Neighbours(i), NSharedNodes(i), ierr)
  END DO

  ALLOCATE(SharedNodes(NoNeighbours, MAXVAL(NSharedNodes)))
  DO i=1, NoNeighbours
    OutProc = i-1
    CALL PMMG_Get_NodeCommunicator_nodesf(pmmgMesh, OutProc, SharedNodes(i, 1:NSharedNodes(i)), ierr)
  END DO

  IF(.NOT. ASSOCIATED(NewMesh % ParallelInfo % INTERFACE)) &
    ALLOCATE(NewMesh % ParallelInfo % Interface(NewMesh % NumberOfNodes))
  NewMesh % ParallelInfo % Interface = .FALSE.
  DO i=1, NewMesh % NumberOfNodes
    counter = 1
    DO j=1, NoNeighbours
      IF(ANY(SharedNodes(j,1:NSharedNodes(j)) == i)) counter=counter+1
    END DO

    ALLOCATE(NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(counter))
    counter = 1
    NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(Counter) = ParEnv % MyPE
    DO j=1, NoNeighbours
      IF(ANY(SharedNodes(j,1:NSharedNodes(j)) == i)) THEN
        counter = counter + 1
        NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(counter) = Neighbours(j)
        NewMesh % ParallelInfo % INTERFACE(i) = .TRUE.
      END IF
    END DO
  END DO

  DO ii=1, NewMesh % NumberOfNodes
    CALL PMMG_Get_VertexGloNum(pmmgMesh, GlobalID, owner, ierr)
    NewMesh % ParallelInfo % GlobalDOFs(ii) = GlobalID
  END DO

  CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

#else
     CALL FATAL('Get_ParMMG_Mesh',&
        'Remeshing utility ParMMG has not been installed')
#endif

END SUBROUTINE Get_ParMMG_Mesh

!A subroutine for a 3D mesh (in parallel) with ParMMG3D
!Inputs:
!   InMesh - the initial mesh
!   Metric - 2D real array specifying target metric
!   EdgePairs, PairCount - 202 edge elems so angle detection is not required
!   NodeFixed, ElemFixed - Optional mask to specify 'required' entities
!Output:
!   OutMesh - the improved mesh
!
SUBROUTINE DistributedRemeshParMMG(Model, InMesh,OutMesh,EdgePairs,PairCount,NodeFixed,ElemFixed,Params)

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount
  LOGICAL :: Success
  !-----------
  TYPE(Mesh_t), POINTER :: WorkMesh
  TYPE(ValueList_t), POINTER :: FuncParams, Material
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:),hminarray(:),hausdarray(:),&
      WorkReal(:,:,:)
  REAL(KIND=dp), POINTER :: WorkArray(:,:) => NULL()
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality, TargetX
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, ElemBodyID, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName, MeshName, SolName, &
       premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile

  SAVE :: WorkReal,WorkArray

#ifdef HAVE_PARMMG

  Debug = .TRUE.
  Parallel = ParEnv % PEs > 1
  FuncName = "RemeshDistParMMG3D"

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  ! params must be provided as not a global mesh
  FuncParams => Params

  !Get parameters from valuelist
  !hausd, hmin, hmax, hgrad, anisoflag, the metric
  !Scalar, vector, tensor metric?

  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hmin", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hmin"')
  MaxRemeshIter= SIZE(WorkArray(:,1))
  ALLOCATE(hminarray(MaxRemeshIter))
  hminarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  hmax = ListGetConstReal(FuncParams, "RemeshMMG3D Hmax",  Default=4000.0_dp)
  hgrad = ListGetConstReal(FuncParams,"RemeshMMG3D Hgrad", Default=0.5_dp)
  WorkArray => ListGetConstRealArray(FuncParams, "RemeshMMG3D Hausd", Found)
  IF(.NOT. Found) CALL FATAL(FuncName, 'Provide hmin input array to be iterated through: "Mesh Hausd"')
  IF(MaxRemeshIter /= SIZE(WorkArray(:,1))) CALL FATAL(FuncName, 'The number of hmin options &
            must equal the number of hausd options')
  ALLOCATE(hausdarray(MaxRemeshIter))
  hausdarray = WorkArray(:,1)
  NULLIFY(WorkArray)
  RemeshMinQuality = ListGetConstReal(FuncParams, "RemeshMMG3D Min Quality",Default=0.0001_dp)
  AnisoFlag = ListGetLogical(FuncParams, "RemeshMMG3D Anisotropic", Default=.TRUE.)
  SaveMMGMeshes = ListGetLogical(FuncParams,"Save RemeshMMG3D Meshes", Default=.FALSE.)
  SaveMMGSols = ListGetLogical(FuncParams,"Save RemeshMMG3D Sols", Default=.FALSE.)
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre RemeshMMG3D Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "RemeshMMG3D Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre RemeshMMG3D Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "RemeshMMG3D Output Sol Name", UnfoundFatal = .TRUE.)
  END IF

  NNodes = InMesh % NumberOfNodes
  NBulk = InMesh % NumberOfBulkElements
  NBdry = InMesh % NumberOfBoundaryElements

  IF(AnisoFlag) THEN

    WorkMesh => Model % Mesh
    Model % Mesh => InMesh

    SolType = MMG5_Tensor
    !Upper triangle of symmetric tensor: 11,12,13,22,23,33
    ALLOCATE(Metric(NNodes,6))
    ALLOCATE(WorkReal(3,1,1))
    Metric = 0.0
    DO i=1,NNodes
      NodeNum = i

      !CALL ListGetRealArray(FuncParams,"RemeshMMG3D Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)

      TargetX = TimeVar % Values(1) / 10
      IF(ABS(InMesh % Nodes % x(i) - TargetX) < 0.1) THEN
        WorkReal(1,1,1) = 0.03_dp
      ELSE IF(ABS(InMesh % Nodes % x(i) - TargetX) < 0.3) THEN
        WorkReal(1,1,1) = 0.1_dp - 0.07_dp * (ABS(ABS(InMesh % Nodes % x(i) - TargetX) - 0.3)/0.3)
      ELSE
        WorkReal(1,1,1) = 0.1_dp
      END IF

      WorkReal(2,1,1) = WorkReal(1,1,1)
      WorkReal(3,1,1) = WorkReal(1,1,1)

      Metric(i,1) = 1.0 / (WorkReal(1,1,1)**2.0)
      Metric(i,4) = 1.0 / (WorkReal(2,1,1)**2.0)
      Metric(i,6) = 1.0 / (WorkReal(3,1,1)**2.0)

    END DO
    DEALLOCATE(WorkReal)
    Model % Mesh => WorkMesh
    WorkMesh => NULL()

  ELSE

    SolType = MMG5_Scalar
    ALLOCATE(Metric(NNodes, 1))
    DO i=1,NNodes
      NodeNum = i
      Metric(i,:) = ListGetReal(FuncParams,"RemeshMMG3D Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
    END DO

  END IF

  body_offset = CurrentModel % NumberOfBCs + CurrentModel % NumberOfBodies + 1
  nBCs = CurrentModel % NumberOfBCs
  DO i=1,InMesh % NumberOfBulkElements
    InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID + body_offset
  END DO


  mmgloops=0
  Success=.TRUE.

10 CONTINUE

  mmgloops = mmgloops+1
  hmin = hminarray(mmgloops)
  Hausd = hausdarray(mmgloops)

  WRITE(Message, '(A,F10.5,A,F10.5)') 'Applying levelset with Hmin ',Hmin, ' and Hausd ', Hausd
  CALL INFO(FuncName, Message)

  pmmgMesh = 0
  pmmgMet  = 0

  !---------------------------------
  ! Issue here: MMG3D will 'helpfully' add any
  ! missing boundary triangles, assigning them
  ! % BoundaryInfo % Constraint = Parent % BodyID
  !
  ! To deal with this, temporarily offset BodyID of
  ! all body elems by Model % NumberOfBodies +
  ! Model % NumberOfBCs + 1, then afterwards we
  ! delete all the extra BC elems & revert the bodyID
  !----------------------------------

  CALL PMMG_Init_parMesh(PMMG_ARG_start, &
      PMMG_ARG_ppParMesh,pmmgMesh, PMMG_ARG_pMesh,PMMG_ARG_pMet, &
      PMMG_ARG_dim,%val(3),PMMG_ARG_MPIComm,%val(ELMER_COMM_WORLD), &
      PMMG_ARG_end)

  IF (Present(PairCount)) THEN
    CALL SET_ParMMG_MESH(InMesh,Parallel,EdgePairs,PairCount)
  ELSE
    CALL SET_ParMMG_MESH(InMesh,Parallel)
  END IF

  !Set the metric values at nodes
  CALL PMMG_Set_MetSize(pmmgMesh, MMG5_Vertex, NNodes, SolType,ierr)
  IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")

  DO i=1,NNodes
    IF(AnisoFlag) THEN
      !IF(Debug) PRINT *,'debug sol at ',i,' is: ',Metric(i,:)
      CALL PMMG_Set_TensorMet(pmmgMesh,Metric(i,1),Metric(i,2),Metric(i,3),&
          Metric(i,4),Metric(i,5),Metric(i,6),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set tensor solution at vertex")
    ELSE
      CALL PMMG_Set_ScalarMet(pmmgMesh,Metric(i,1),i,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
    END IF
  END DO

  !Turn on debug (1)
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_debug, &
       1,ierr)

  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmin,&
       hmin,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmax,&
       hmax,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hausd,&
       hausd,ierr)
  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hgrad,&
       hgrad,ierr)
  ! allow surface modifications
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_nosurf,&
       0,ierr)

  ! compute globaldofs
  CALL PMMG_SET_IPARAMETER(pmmgMesh, PMMGPARAM_globalnum,&
       1, ierr)

  !Turn off sharp angle detection (0)
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle, &
       0,ierr)
  !Option to set angle detection threshold:
  !CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_angleDetection,&
  !      85.0_dp,ierr)

  !Take care of fixed nodes/elements if requested
  IF(PRESENT(NodeFixed)) THEN
    DO i=1,NNodes
      IF(NodeFixed(InMesh % ParallelInfo % GlobalDOFs(i))) THEN
        CALL PMMG_SET_REQUIREDVERTEX(pmmgMesh,i,ierr)
      END IF
    END DO
  END IF

  IF(PRESENT(ElemFixed)) THEN
    DO i=1,NBulk + NBdry
      IF(ElemFixed(InMesh % Elements(i) % GElementIndex)) THEN
        IF(InMesh % Elements(i) % TYPE % NumberOfNodes == 4) THEN
          CALL PMMG_SET_REQUIREDTETRAHEDRON(pmmgMesh,i,ierr)
        ELSE
          CALL PMMG_SET_REQUIREDTRIANGLE(pmmgMesh,i-NBulk,ierr)
        END IF
      END IF
    END DO
  END IF

  !! need to set face communicators eg neighbour procs either nodes or faces

  IF(SaveMMGMeshes) THEN
    WRITE(MeshName, '(A,i0,A)') TRIM(premmg_meshfile), time, '.mesh'
    CALL PMMG_SaveMesh_Distributed(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
  END IF
  IF(SaveMMGSols) THEN
    WRITE(SolName, '(A,i0,A)') TRIM(premmg_solfile), time, '.sol'
    CALL PMMG_SaveMet_Distributed(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
  END IF
  IF (DEBUG) PRINT *,'--**-- SET PMMG3D PARAMETERS '
  ! CALL SET_MMG3D_PARAMETERS(SolverParams)

  CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
  CALL PMMG_parmmglib_distributed(pmmgMesh,ierr)
  IF(ierr == PMMG_LOWFAILURE) THEN
    PRINT*, ParEnv % MyPE, 'low failure'
  ELSE IF ( ierr == PMMG_STRONGFAILURE .OR. ierr == PMMG_LOWFAILURE ) THEN
    PRINT*,"BAD ENDING OF PMMGLIB: UNABLE TO SAVE MESH", ParEnv % MyPE
    !! Release mmg mesh
    CALL FATAL('bad', 'ending of remeshing')
    CALL PMMG_Free_all ( PMMG_ARG_start,     &
        PMMG_ARG_ppParMesh,pmmgMesh,         &
        PMMG_ARG_end)
    WRITE(Message, '(A,F10.5,A,F10.5)') 'Remesh failed with Hmin ',Hmin, ' and Hausd ', Hausd
    CALL WARN(FuncName, Message)
    IF(mmgloops==MaxRemeshIter) THEN
      Success=.FALSE.
      GO TO 20
    ELSE
      GO TO 10
    END IF
    !STOP MMG5_STRONGFAILURE
  ENDIF

  IF (DEBUG) PRINT *,'--**-- PMMG_parmmglib_centralized DONE'

  IF(SaveMMGMeshes) THEN
    WRITE(MeshName, '(A,i0,A)') TRIM(mmg_meshfile), time, '.mesh'
    CALL PMMG_SaveMesh_Distributed(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
  END IF
  IF(SaveMMGSols) THEN
    WRITE(SolName, '(A,i0,A)') TRIM(mmg_solfile), time, '.sol'
    CALL PMMG_SaveMet_Distributed(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
  END IF

  CALL PMMG_Get_meshSize(pmmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

  !! GET THE NEW MESH
  CALL GET_ParMMG_MESH(OutMesh,Parallel)

  !! Release mmg mesh
  CALL PMMG_Free_all ( PMMG_ARG_start,     &
       PMMG_ARG_ppParMesh,pmmgMesh,         &
       PMMG_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements

  !Reset the BodyIDs (see above)
  DO i=1,NBulk
    OutMesh % Elements(i) % BodyID = OutMesh % Elements(i) % BodyID - body_offset
  END DO

  !And delete the unneeded BC elems
  ALLOCATE(RmElement(NBulk+NBdry))
  RmElement = .FALSE.
  DO i=NBulk+1,NBulk+NBdry
    Element => OutMesh % Elements(i)
    IF(Element % BoundaryInfo % Constraint > nBCs) THEN
      RmElement(i) = .TRUE.
    END IF
    ! on parallel interface
    IF(Element % BoundaryInfo % Constraint == 0 .AND. i>NBulk) THEN
      RmElement(i) = .TRUE.
    END IF
  END DO
  CALL CutMesh(OutMesh, RmElem=RmElement)

20 CONTINUE

  ! if remeshing has failed need to reset body ids
  IF(.NOT. Success) THEN
    DO i=1,InMesh % NumberOfBulkElements
      InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID - body_offset
    END DO
  END IF

#else
  CALL Fatal(FuncName, "Remeshing utility PMMG has not been installed")
#endif

END SUBROUTINE DistributedRemeshParMMG

END MODULE MeshRemeshing
