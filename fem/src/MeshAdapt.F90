!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh adaptation utilities including interfaces to MMG & Zoltan
!------------------------------------------------------------------------------

MODULE MeshAdapt

!TODO - preprocessor commands

USE DefUtils
USE Types
USE MeshUtils

#ifdef HAVE_ZOLTAN
USE Zoltan
#endif

IMPLICIT NONE

#ifdef HAVE_MMG
#include "mmg/mmg3d/libmmg3df.h"
#endif

#ifdef HAVE_MMG
INTEGER :: MMGPARAM_hausd = MMG3D_DPARAM_hausd
INTEGER :: MMGPARAM_hmin = MMG3D_DPARAM_hmin
INTEGER :: MMGPARAM_hmax = MMG3D_DPARAM_hmax
INTEGER :: MMGPARAM_iso = MMG3D_IPARAM_iso

MMG5_DATA_PTR_T :: mmgMesh
MMG5_DATA_PTR_T :: mmgSol
#endif

CONTAINS

!============================================
!============================================
!            MMG3D SUBROUTINES
!============================================
!============================================


SUBROUTINE Set_MMG3D_Mesh(Mesh)

  TYPE(Mesh_t), POINTER :: Mesh

#ifdef HAVE_MMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr
  INTEGER, ALLOCATABLE :: NodeRefs(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.FALSE.

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
  CALL MMG3D_Set_meshSize(mmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ierr)
  IF ( ierr /= 1 ) CALL FATAL('MMGSolver',&
       'CALL TO MMG3D_Set_meshSize FAILED')
  IF (DEBUG) PRINT *,'--**-- MMG3D_Set_meshSize DONE'

  DO i=1,NVerts
    CALL MMG3D_Set_vertex(mmgMesh, Mesh%Nodes%x(i), &
         Mesh%Nodes%y(i),Mesh%Nodes%z(i), 0, i, ierr)
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
           Element % BoundaryInfo % Constraint, ntris,ierr)
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


SUBROUTINE Get_MMG3D_Mesh(NewMesh)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
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
         ref,corner,required,ierr)
    IF ( ierr == 0 ) CALL FATAL('MMGSolver',&
         'CALL TO  MMG3D_Get_vertex FAILED')
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


!============================================
!============================================
!          ZOLTAN SUBROUTINES
!============================================
!============================================

SUBROUTINE Zoltan_Interface( Model, Solver, dt, Transient )

  USE MeshUtils

#ifdef HAVE_ZOLTAN
  USE Zoltan
#endif

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------

#ifdef HAVE_ZOLTAN
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element,Element2,Face, MFacePtr(:)
  TYPE(Graph_t) :: LocalGraph
  REAL(KIND=dp) :: t1,t2
  INTEGER :: i,j,k,l,m,nn,ierr,NNodes,NBulk,NElnodes,counter,DIM,&
       NFaces,NIFFaces,max_elemno
  INTEGER, ALLOCATABLE :: ElemAdj(:), ElemStart(:),ParElemAdj(:), ParElemStart(:),&
       ParElemIdx(:),ParElemAdjProc(:),FaceIFIDX(:),FaceOrder(:),sharecount(:),&
       ParElemMap(:)
  TYPE(NeighbourList_t), POINTER :: MFaceIFList(:)
  LOGICAL, POINTER :: MFaceIF(:)
  
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Zoltan_Interface"

  !Zoltan things
  TYPE(Zoltan_Struct), POINTER :: zz_obj
  INTEGER(Zoltan_INT) :: zierr,numGIDEntries, numLidEntries, numImport, numExport
  INTEGER(Zoltan_INT),DIMENSION(:), POINTER :: importGlobalGids, importLocalGids, importProcs, &
       importToPart,exportGlobalGids,exportLocalGids, exportProcs, exportToPart
  REAL(Zoltan_FLOAT) :: version
  LOGICAL :: changes,Debug

  TYPE ElemTable_t
     INTEGER :: counter=0
     INTEGER, ALLOCATABLE :: Idx(:)
  END TYPE ElemTable_T
  TYPE(ElemTable_t),ALLOCATABLE :: NodeElems(:),ElemElems(:)


  Mesh => Model % Mesh
  NNodes = Mesh % NumberOfNodes
  NBulk = Mesh % NumberOfBulkElements
  DIM = CoordinateSystemDimension()

  zierr = Zoltan_Initialize(version)
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to initialize Zoltan partitioner")

  NULLIFY(zz_obj)

  zz_obj => Zoltan_Create(ELMER_COMM_WORLD)

  zierr = Zoltan_Set_Param(zz_obj, "LB_METHOD", "GRAPH")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan param LB_METHOD")
  zierr = Zoltan_Set_Param(zz_obj, "LB_APPROACH", "REFINE") !REPARTITION/REFINE <- faster
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: LB_APPROACH")
  zierr = Zoltan_Set_Param(zz_obj, "GRAPH_PACKAGE", "PHG")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: GRAPH_PACKAGE")
  zierr = Zoltan_Set_Param(zz_obj, "NUM_GID_ENTRIES", "1")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: NUM_GID_ENTRIES")
  zierr = Zoltan_Set_Param(zz_obj, "NUM_LID_ENTRIES", "1")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: NUM_LID_ENTRIES")
  zierr = Zoltan_Set_Param(zz_obj, "RETURN_LISTS", "ALL")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: RETURN_LISTS")
  zierr = Zoltan_Set_Param(zz_obj, "OBJ_WEIGHT_DIM", "0")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: OBJ_WEIGHT_DIM")
  zierr = Zoltan_Set_Param(zz_obj, "EDGE_WEIGHT_DIM", "0")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: EDGE_WEIGHT_DIM")
  zierr = Zoltan_Set_Param(zz_obj, "DEBUG_LEVEL", "0")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: DEBUG_LEVEL")
  zierr = Zoltan_Set_Param(zz_obj, "CHECK_GRAPH", "0")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: CHECK_GRAPH")
  zierr = Zoltan_Set_Param(zz_obj, "PHG_MULTILEVEL", "0")
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: PHG_MULTILEVEL")

  !Callback functions to query number of elements and the element data
  zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zoltNumObjs)
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: 12")
  zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
  IF(zierr /= 0) CALL Fatal(FuncName,"Unable to set Zoltan Parameter: 13")

  !Callback functions to query number of edges and the edge data
  zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_EDGES_FN_TYPE,zoltNumEdges)
  IF(zierr /= 0) CALL Fatal(FuncName,"meow14")
  zierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_EDGE_LIST_FN_TYPE,zoltGetEdgeList)
  IF(zierr /= 0) CALL Fatal(FuncName,"meow15")


  !Parameters to set to get graph repartitioning (without relying on libparmetis?!)
  !LB_APPROACH = REPARTITION (OR REFINE)
  !LB_METHOD = GRAPH

  !Need the following functions for repartitioning:

  ! ZOLTAN_NUM_OBJ_FN
  ! ZOLTAN_OBJ_LIST_FN
  ! ZOLTAN_NUM_EDGES_MULTI_FN or ZOLTAN_NUM_EDGES_FN 
  ! ZOLTAN_EDGE_LIST_MULTI_FN or ZOLTAN_EDGE_LIST_FN

  ! ZOLTAN_OBJ_SIZE_MULTI_FN or ZOLTAN_OBJ_SIZE_FN  - Optional for LB_APPROACH=Repartition.
  ! ZOLTAN_PART_MULTI_FN or ZOLTAN_PART_FN - Optional for LB_APPROACH=Repartition and for REMAP=1. 

  t1 = CPUTime()
  CALL ElmerMeshToDualGraph(Mesh, LocalGraph)
  t2 = CPUTime()

  ! t1 = CPUTime()
  ! CALL GetBulkElemAdjacency( Mesh,ElemAdj,ElemStart )
  ! t2 = CPUTime()

  IF(Debug) THEN
    PRINT *,'size graph: ',SIZE(LocalGraph % ptr), SIZE(LocalGraph % ind), LocalGraph % n
    PRINT *,'start graph locs: ',LocalGraph % ptr(1:100)
    PRINT *,'start graph neighs: ',LocalGraph % ind(1:100)
  END IF
  !CALL GetParallelElemAdjacency( Mesh, ParElemAdj, ParElemStart)
  CALL MeshParallelDualGraph( Mesh, ParElemAdj, ParElemStart, ParElemIdx, ParElemAdjProc )

  PRINT *,ParEnv % MyPE,' final sizes ',SIZE(ParElemAdj), SIZE(ParElemStart), SIZE(ParElemIdx)

  !Construct map of ParElemIdx
  ALLOCATE(ParElemMap(NBulk))
  ParElemMap = 0
  DO i=1,SIZE(ParElemIDX)
    ParElemMap(ParElemIDX(i)) = i
  END DO

  CALL MPI_ALLREDUCE(MAXVAL(Mesh % Elements % GElementIndex), max_elemno, 1, MPI_INTEGER, &
       MPI_MAX, ELMER_COMM_WORLD, ierr)
  PRINT *,parenv % mype,' max elemno: ',max_elemno

  !Combine local and parallel element idx, and convert to global
  DO i=1,NBulk
    k = LocalGraph % ptr(i+1) - LocalGraph % ptr(i)
    IF(ParElemMap(i) > 0) k = k + ParElemStart(ParElemMap(i)+1) - ParElemStart(ParElemMap(i))
  END DO


  numGidEntries = 1
  numLidEntries = 1
  
  zierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
       numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
       numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)
  IF(zierr /= 0) CALL Fatal("meow","meow16")

  PRINT *,ParEnv % MyPE,' zoltan import: ',numImport, ' export: ',numExport,' total: ',NBulk


  CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER FUNCTION zoltNumObjs(DATA, ierr)
    use zoltan 
    implicit none

    ! Local declarations
    INTEGER(Zoltan_INT), INTENT(in) :: DATA(*)
    INTEGER(ZOLTAN_INT), INTENT(out) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zoltNumObjs = NBulk
    ierr = ZOLTAN_OK
    PRINT *,ParEnv % MyPE,' nbulk: ',zoltNumObjs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END FUNCTION zoltNumObjs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGetObjs (data, num_gid_entries, num_lid_entries, global_ids, & 
                        local_ids, wgt_dim, obj_wgts, ierr)
    use zoltan
    implicit none

    INTEGER(ZOLTAN_INT), INTENT(in) :: DATA(*)
    !TYPE(Mesh_t), POINTER, INTENT(in) :: DATA
    ! TYPE(Zoltan_User_Data_1) :: DATA
    integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
    integer(ZOLTAN_INT), intent(out) :: global_ids(*)
    integer(ZOLTAN_INT), intent(out) :: local_ids(*)
    integer(ZOLTAN_INT), intent(in) :: wgt_dim 
    real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    ! local declarations
    integer :: i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i= 1, NBulk
      global_ids(i) = Mesh % Elements(i) % GElementIndex
      local_ids(i) = i
    end do
    
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION zoltNumEdges(DATA, num_gid_entries, num_lid_entries, global_id, local_id, ierr)
    INTEGER(ZOLTAN_INT), INTENT(in) :: DATA  
    INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries, num_lid_entries  
    INTEGER(Zoltan_INT), INTENT(IN) :: global_id  
    INTEGER(Zoltan_INT), INTENT(IN) :: local_id  
    INTEGER(Zoltan_INT), INTENT(OUT) :: ierr

    zoltNumEdges = (LocalGraph % ptr(local_id+1) - LocalGraph % ptr(local_id))

    IF(ParElemMap(local_id) /= 0) THEN
      zoltNumEdges = zoltNumEdges + &
           (ParElemStart(ParElemMap(local_id)+1) - ParElemStart(ParElemMap(local_id)))
    END IF

    ! PRINT *,parenv % mype,' debug ',local_id,ParElemMap(local_id),' num edges: ',zoltNumEdges
    ierr = 0
  END FUNCTION zoltNumEdges
  
  SUBROUTINE zoltGetEdgeList(DATA, num_gid_entries, num_lid_entries, global_id, local_id, &
       nbor_global_id, nbor_procs, wgt_dim, ewgts, ierr)

    !TYPE(Mesh_t), POINTER, INTENT(in) :: DATA
    ! TYPE(Zoltan_User_Data_1) :: DATA
    INTEGER(ZOLTAN_INT), INTENT(in) :: DATA  
    INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries, num_lid_entries  
    INTEGER(Zoltan_INT), INTENT(IN) :: global_id
    INTEGER(Zoltan_INT), INTENT(IN) :: local_id
    INTEGER(Zoltan_INT), INTENT(OUT) :: nbor_global_id(*)
    INTEGER(Zoltan_INT), INTENT(OUT) :: nbor_procs(*)
    INTEGER(Zoltan_INT), INTENT(IN) :: wgt_dim(*)
    REAL(Zoltan_FLOAT), INTENT(OUT) :: ewgts(*)
    INTEGER(Zoltan_INT), INTENT(OUT) :: ierr 
    !----------------
    INTEGER :: counter,i,k,nlocal,nother

    nlocal = (LocalGraph % ptr(local_id+1) - LocalGraph % ptr(local_id))

    !Something wrong here
    DO i=1,nlocal
      k = i + LocalGraph % ptr(local_id) - 1
      ! PRINT *,ParEnv % MyPE,' debug setting nbor_procs',i,' nlocal: ',nlocal,&
      !      Mesh % Elements(LocalGraph % ind(k)) % GElementIndex
      nbor_global_id(i) = Mesh % Elements(LocalGraph % ind(k)) % GElementIndex
      nbor_procs(i) = ParEnv % MyPE
    END DO

    nother = 0
    IF(ParElemMap(local_id) /= 0) THEN

      nother = ParElemStart(ParElemMap(local_id)+1) - ParElemStart(ParElemMap(local_id))
      DO i=1,nother
        j = ParElemStart(ParElemMap(local_id)) + i - 1
        k = i + nlocal
        nbor_global_id(k) = ParElemAdj(j)
        nbor_procs(k) = ParElemAdjProc(j)
      END DO
    END IF

    IF(ANY(nbor_global_id(1:nlocal+nother) < 1)) CALL FATAL(FuncName,"Bad gid 0")
    IF(ANY(nbor_global_id(1:nlocal+nother) > max_elemno)) CALL FATAL(FuncName,"Bad gid max")

    !Some checks on the data - slow!
    DO i=1,nlocal+nother-1
      IF(ANY(nbor_global_id(i+1:nlocal+nother) == nbor_global_id(i))) THEN
        PRINT *,ParEnv % MyPE,' duplicates! :',nbor_global_id(i),' with ',nbor_global_id(1:nlocal+nother)
      END IF
    END DO
    IF(ANY(nbor_procs(1:nlocal+nother) > 3)) CALL FATAL(FuncName,"Bad proc max")
    IF(ANY(nbor_procs(1:nlocal+nother) < 0)) CALL FATAL(FuncName,"Bad proc min")

    ! PRINT *,ParEnv % MyPE,' debug el: ',global_id,&
    !      'local ',local_id,' of ',NBulk,&
    !      ParElemMap(local_id),' global neighs neigh: ',&
    !      nbor_global_id(1:nlocal+nother)

         

  END SUBROUTINE ZoltGetEdgeList

#else
     CALL FATAL('Zoltan_Interface',&
        'Repartitioning utility Zoltan (Trilinos) has not been installed')
#endif
END SUBROUTINE Zoltan_Interface

END MODULE MeshAdapt
