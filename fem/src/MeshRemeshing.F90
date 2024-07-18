!/****************************************************************************
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

!> \ingroup ElmerLib
!> \{

MODULE MeshRemeshing

USE Types
USE Lists
USE Messages
USE MeshUtils
USE MeshPartition
USE SparIterComm

IMPLICIT NONE

#ifdef HAVE_MMG
#include "mmg/mmg3d/libmmg3df.h"
#include "mmg/mmg2d/libmmg2df.h"
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
INTEGER :: MMGPARAM_hgradreq = MMG3D_DPARAM_hgradreq
INTEGER :: MMGPARAM_nosizreq = MMG3D_IPARAM_nosizreq
MMG5_DATA_PTR_T :: mmgMesh
MMG5_DATA_PTR_T :: mmgSol
MMG5_DATA_PTR_T :: mmgMet
MMG5_DATA_PTR_T :: mmgLs

MMG5_DATA_PTR_T :: pmmgMesh
#endif

#ifdef HAVE_PARMMG
#include "parmmg/libparmmgtypesf.h"
#endif

INTEGER, PARAMETER :: ElmerBCOffset = 1000

CONTAINS

!============================================
!============================================
!            MMG3D SUBROUTINES
!============================================
!============================================

  
SUBROUTINE Set_MMG3D_Mesh(Mesh, Parallel, EdgePairs, PairCount, Solver)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount
  TYPE(Solver_t), POINTER, OPTIONAL :: Solver
  
#ifdef HAVE_MMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)
  
  INTEGER :: i,j,NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ref,ierr
  INTEGER, ALLOCATABLE :: NodeRefs(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.FALSE.,Elem202,Found
  INTEGER, POINTER :: Perm(:)
  CHARACTER(:), ALLOCATABLE :: EquationName                              
  CHARACTER(*), PARAMETER :: FuncName="Set_MMG3D_Mesh"

  IF(CoordinateSystemDimension() /= 3) CALL Fatal(FuncName,"Only works for 3D meshes!")

  ALLOCATE(NodeRefs(6))

  IF(Parallel) CALL Assert(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs), FuncName,&
      "Parallel sim but no ParallelInfo % GlobalDOFs")
  
  IF( PRESENT( Solver ) ) THEN
    CALL Info(FuncName,'Setting only the active part of mesh')
    EquationName = ListGetString( Solver % Values, 'Equation', Found)                 
    Perm => Solver % Variable % Perm
    IF(.NOT. ASSOCIATED(Perm)) THEN
      CALL Fatal(FuncName,'Perm vector not associated!?')
    END IF
    NVerts = COUNT( Perm(1:Mesh % NumberOfNodes) > 0 )
  ELSE
    NULLIFY( Perm )
    NVerts = Mesh % NumberOfNodes    
  END IF
  CALL Info(FuncName,'Set number of nodes: '//I2S(Nverts),Level=20)
  
  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0
  IF(PRESENT(PairCount)) NEdges= PairCount

  DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    Element => Mesh % Elements(i)
    IF( PRESENT( Solver ) ) THEN
      IF( i <= Mesh % NumberOfBulkElements ) THEN
        IF ( .NOT. CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
      ELSE
        NodeIndexes => Element % NodeIndexes
        IF(ANY(Perm(NodeIndexes) == 0 ) ) CYCLE
      END IF
    END IF
    
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
      CALL Fatal(FuncName,"can't handle pyramid elements (605)")
    CASE(706)
      NPrisms = NPrisms + 1
    CASE(808)
      CALL Fatal(FuncName,"can't handle brick/hexahedral elements (808)")
    CASE DEFAULT
      CALL Fatal(FuncName,"Unsupported element type: "//TRIM(I2S(Element % TYPE % ElementCode)))
    END SELECT
  END DO

  Nbulk = NTetras + NPrisms
  Nbdry = NTris + Nquads

  CALL Info(FuncName,'Set number of bulk elements: '//I2S(Nbulk),Level=20)
  CALL Info(FuncName,'Set number of boundary elements: '//I2S(Nbdry),Level=20)

  
  IF(Warn101) CALL Warn(FuncName,"101 elements detected - these won't be remeshed")
  IF(Warn202) CALL Warn(FuncName,"202 elements detected - these won't be remeshed")
  
  !args: mesh, nvertices, ntetra, nprisms, ntriangless, nquads, nedges
  CALL MMG3D_Set_meshSize(mmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ierr)
  IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to MMG3D_Set_meshSize failed!')
  CALL Info(FuncName,'MMG3D_Set_meshSize done',Level=20)
  
  ref = 0
  DO i=1,Mesh % NumberOfNodes
    j = i
    IF( ASSOCIATED(Perm)) THEN
      j = Perm(i)
      IF(j==0) CYCLE
    END IF
    ! ref = GDOF + 10 to avoid an input ref of 10 being confused
    ! with mmg output ref of 10 which occurs on some new nodes
    ! GDOF = ref - 10
    IF(Parallel) ref = Mesh % ParallelInfo % GlobalDOFs(i) + 10
    CALL MMG3D_Set_vertex(mmgMesh, Mesh % Nodes % x(i), &
        Mesh % Nodes % y(i),Mesh % Nodes % z(i), ref, j, ierr)
    IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to MMG3D_Set_vertex failed!')
  END DO
  CALL Info(FuncName,'MMG3D_Set_vertex done',Level=20)

  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0

  ! Cycle mesh elements, sending them to MMG3D
  DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    Element => Mesh % Elements(i)
    NNodes = Element % TYPE % NumberOfNodes
    NodeIndexes => Element % NodeIndexes

    IF(PRESENT(Solver)) THEN
      NodeRefs(1:NNodes) = Perm(NodeIndexes(1:NNodes))
      IF(ANY(NodeRefs(1:NNodes) == 0 ) ) CYCLE      
      IF(i<=Mesh % NumberOfBulkElements) THEN
        ! We could be unlucky and have a bulk elements which  is actually not active even though all Perm>0
        IF ( .NOT. CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
      END IF
    ELSE
      NodeRefs(1:NNodes) = NodeIndexes(1:NNodes)
    END IF
    
!    NodeRefs(1:NNodes) = Mesh % ParallelInfo % GlobalDOFs(NodeIndexes(1:NNodes))

!    PRINT *,'debug, elem ',i,' noderefs: ',NodeRefs(1:NNodes)
    SELECT CASE(Element % TYPE % ElementCode)
    CASE(303)
      ntris = ntris + 1
      CALL MMG3D_Set_triangle(mmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
          Element % BoundaryInfo % Constraint+ElmerBCOffset, ntris, ierr)
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
  CALL Info(FuncName,'Set elements done',Level=20)


  ! Freeze nodes at the interface of remeshed and omitted volume
  IF( PRESENT( Solver ) ) THEN
    BLOCK
      LOGICAL, ALLOCATABLE :: TagNode(:)
      ALLOCATE(TagNode(Mesh % NumberOfNodes))
      TagNode = .FALSE.

      DO i=1,Mesh % NumberOfBulkElements 
        Element => Mesh % Elements(i)
        IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
        NodeIndexes => Element % NodeIndexes
        TagNode(NodeIndexes) = .TRUE.
      END DO

      j = 0
      DO i=1,Mesh % NumberOfNodes
        IF(TagNode(i) .AND. Perm(i) > 0 ) THEN
          CALL MMG3D_SET_REQUIREDVERTEX(mmgMesh,Perm(i),ierr)
          j = j+1
        END IF
      END DO
      CALL Info(FuncName,'Fixed number of nodes at interface: '//I2S(j),Level=10)
    END BLOCK
  END IF
    
  
  ! use element pairs '202' elements
  Elem202 = (PRESENT(EdgePairs))
  IF (Elem202) THEN
    DO i=1, PairCount
      NEdges = NEdges + 1 
      CALL MMG3D_Set_edge(mmgMesh, EdgePairs(1,i), EdgePairs(2,i), 1, nedges, ierr)
      CALL MMG3D_Set_ridge(mmgMesh, nedges, ierr)
      !CALL MMG3D_Set_requiredEdge(mmgMesh, Nedges, ierr)
    END DO
  END IF

  CALL Info(FuncName,'Set edge elements done',Level=20)

#else
  CALL Fatal('Set_MMG3D_Mesh',&
        'Remeshing utility MMG3D has not been installed')
#endif

END SUBROUTINE Set_MMG3D_Mesh


SUBROUTINE Check_Parameters_Obsolite(SolverParams)

  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Checked = .FALSE.

  IF(Checked) RETURN
  
  IF( ListCheckPrefix( SolverParams,'RemeshMMG3D') ) THEN
    CALL Fatal('Check_Parameters_Obsolite','Use "MMG" as prefix instead of "RemeshMMG3D"')
  END IF
  CALL ListObsoliteFatal(SolverParams,'hmin','MMG hmin')
  CALL ListObsoliteFatal(SolverParams,'hmax','MMG hmax')
  CALL ListObsoliteFatal(SolverParams,'hsiz','MMG hsiz')
  CALL ListObsoliteFatal(SolverParams,'hausd','MMG hausd')
  CALL ListObsoliteFatal(SolverParams,'hgrad','MMG hgrad')
  CALL ListObsoliteFatal(SolverParams,'verbosity','MMG verbosity')
  CALL ListObsoliteFatal(SolverParams,'angle detection','MMG angle detection')
  CALL ListObsoliteFatal(SolverParams,'no angle detection','MMG no angle detection')
  CALL ListObsoliteFatal(SolverParams,'increase memory','MMG increase memory')
  CALL ListObsoliteFatal(SolverParams,'no insert','MMG noinsert')
  CALL ListObsoliteFatal(SolverParams,'no swap','MMG no swap')
  CALL ListObsoliteFatal(SolverParams,'no move','MMG no move')
  CALL ListObsoliteFatal(SolverParams,'no surf','MMG no surf')

  Checked = .TRUE.
  
END SUBROUTINE Check_Parameters_Obsolite
  
  

SUBROUTINE Set_MMG3D_Parameters(SolverParams, ReTrial)

  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL, OPTIONAL :: ReTrial
  
#ifdef HAVE_MMG
  REAL(KIND=dp) :: Pval
  LOGICAL :: NoAngleDetect
  INTEGER :: verbosity,MeMIncrease,Bucket,GMshOption,ierr
  REAL(KIND=dp) :: Hmin, Hmax, HSiz
  LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
  LOGICAL :: Found,Debug=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Set_MMG3D_Parameters"

  IF(.NOT. ASSOCIATED(SolverParams) ) THEN
    CALL Fatal(FuncName,'"SolverParams" is not associated!')
  END IF

  CALL Check_Parameters_Obsolite(SolverParams)
  
  ! Minimal mesh size:  hmin
  hmin = ListGetCReal( SolverParams,'adaptive min h', Found ) 
  IF(.NOT. Found) Hmin = ListGetCReal( SolverParams,'mmg hmin', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,Hmin,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMG3D_SET_DPARAMETER <hmin> Failed')
  END IF
  
  ! Maximal mesh size - hmax
  hmax = ListGetCReal( SolverParams,'adaptive max h', Found ) 
  IF(.NOT. Found) Hmax = ListGetCReal( SolverParams, 'mmg hmax', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,Hmax,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMG3D_SET_DPARAMETER <hmax> Failed')
  END IF

  hsiz = ListGetConstReal( SolverParams, 'mmg hsiz', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hsiz,hsiz,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMG3D_SET_DPARAMETER <hsiz> Failed')
  END IF

  ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
  ! adaptive hausd used in 3D calving 
  Pval = ListGetCReal( SolverParams,'adaptive hausd', Found ) 
  IF(.NOT. FOund) Pval = ListGetCReal( SolverParams, 'mmg hausd', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMG3D_SET_DPARAMETER <hausd> Failed')
  END IF

  ! Control gradation
  Pval = ListGetConstReal( SolverParams, 'mmg hgrad', Found)
  IF (Found) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMG3D_SET_DPARAMETER <hgrad> Failed')
  END IF

  ! If this is a ReTrial then we only change some real valued keywords!
  IF( PRESENT( ReTrial ) ) THEN
    IF( ReTrial ) RETURN
  END IF
    
!!! PARAMS: generic options (debug, mem, verbosity)
  ! [val] Set the verbosity level to n
  Verbosity = ListGetInteger( SolverParams,'mmg verbosity',Found)
  IF (Found) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_verbose,Verbosity,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
        'Call to MMG3D_SET_IPARAMETER <verbose> Failed')
  END IF

  ! [val] Set the maximal memory size to n MBytes.
  MemIncrease = ListGetInteger(SolverParams,'mmg Increase Memory',Found)
  IF (Found) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_mem,MemIncrease,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
        'Call to MMG3D_SET_IPARAMETER <mem> Failed')
  END IF

  ! [0/1] Turn on the debug mode.
  DebugMode = ListGetLogical(SolverParams,'mmg Debug Mode',Found)
  IF(.NOT. Found) DebugMode = ListGetLogical(SolverParams,'Debug Mode',Found)
  IF (DebugMode) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_debug,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
        'Call to MMG3D_SET_IPARAMETER <debug> Failed')
  END IF

  ! !< [1/0], Avoid/allow automatic angle detection */
  NoAngleDetect = ListGetLogical(SolverParams,'mmg No Angle detection',Found)
  IF (NoAngleDetect) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle,0,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
        'Call to MMG3D_SET_IPARAMETER <No Angle detection> Failed')
  ELSE
    !This is the default!
    !CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle,1,ierr)
    !IF ( ierr == 0 ) CALL Fatal(FuncName,&
    !    'Call to MMG3D_SET_IPARAMETER <No Angle detection> Failed')
    !!! mmg angle detection angle
    Pval = ListGetCReal( SolverParams, 'mmg Angle detection',Found)
    IF (Found) THEN
      CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
          Pval,ierr)
      IF ( ierr == 0 ) CALL Fatal(FuncName,&
          'Call to MMG3D_SET_DPARAMETER <Angle detection> Failed')
    ELSE
      CALL WARN(FuncName, "Using mmg default value for automatic angle detection")
    ENDIF
  END IF

  ! [1/0] Avoid/allow point insertion
  NoInsert = ListGetLogical(SolverParams,'MMG No insert',Found)
  IF (NoInsert) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noinsert,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_IPARAMETER <No insert> Failed') 
  END IF

  ! [1/0] Avoid/allow edge or face flipping
  NoSwap = ListGetLogical(SolverParams,'MMG No swap',Found)
  IF (NoSwap) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_noswap,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <No swap> Failed')
  END IF

  ! [1/0] Avoid/allow point relocation
  NoMove = ListGetLogical(SolverParams,'MMG No move',Found)
  IF (NoMove) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nomove,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <No move> Failed')
  END IF

  ! [1/0], Preserve triangles at interface of 2 domains with same reference
  IF( ListGetLogical(SolverParams,'MMG open body',Found) ) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_opnbdy,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <opnbdy> Failed')
  END IF

  NoSurf = ListGetLogical(SolverParams,'MMG No surf',Found)
  IF(.NOT. Found) NoSurf = .TRUE.
  IF (NoSurf) THEN
    CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_nosurf,0,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <No surf> Failed')
  END IF
  
  Pval = ListGetCReal( SolverParams, 'MMG HgradReq',Found)
  IF( Found ) THEN
    CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_HgradReq,&
        Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
        'Call to MMG3D_SET_DPARAMETER <Angle detection> Failed')
  END IF
    
!!!
#else
     CALL Fatal('Set_MMG3D_Parameters',&
        'Remeshing utility MMG3D has not been installed')
#endif
   END SUBROUTINE Set_MMG3D_Parameters

   
SUBROUTINE Set_PMMG_Parameters(SolverParams, ReTrial )

  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL, OPTIONAL :: ReTrial
  
#ifdef HAVE_PARMMG
  REAL(KIND=dp) :: Pval
  LOGICAL :: NoAngleDetect
  INTEGER :: verbosity,MeMIncrease,Bucket,GMshOption,ierr,niter
  REAL(KIND=dp) :: Hmin, Hmax, HSiz
  LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
  LOGICAL :: Found,Debug=.FALSE.,ParMMG
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Set_PMMG_Parameters"

  CALL Check_Parameters_Obsolite(SolverParams)
  
  ! Minimal mesh size:  hmin
  hmin = ListGetCReal( SolverParams,'adaptive min h', Found ) 
  IF(.NOT. Found) Hmin = ListGetCReal( SolverParams,'mmg hmin', Found)
  
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmin,Hmin,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_DPARAMETER <hmin> Failed')
  END IF
  
  ! Maximal mesh size - hmax
  hmax = ListGetCReal( SolverParams,'adaptive max h', Found ) 
  IF(.NOT. Found) Hmax = ListGetCReal( SolverParams, 'mmg hmax', Found)
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmax,Hmax,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
        'Call to MMG3D_SET_DPARAMETER <hmax> Failed')
  END IF

  !hsiz = ListGetConstReal( SolverParams, 'mmg hsiz', Found)
  !IF(.NOT. Found) hsiz = ListGetConstReal( SolverParams, 'hsiz', Found)
  !IF (Found) THEN
  !  CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hsiz,hsiz,ierr)
  !  IF ( ierr == 0 ) CALL Fatal('MMGSolver', &
  !       'Call to MMG3D_SET_DPARAMETER <hsiz> Failed')
  !END IF

  ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
  ! MMG3D_DPARAM_hausd default est 0.01 semble bien trop petit;
  !  il semble qu'il faille mettre une taille > taille des elements.
  Pval = ListGetCReal( SolverParams,'adaptive hausd', Found ) 
  IF(.NOT. FOund) Pval = ListGetCReal( SolverParams, 'mmg hausd', Found)
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hausd,Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_DPARAMETER <hausd> Failed')
  END IF

  ! Control gradation
  Pval = ListGetConstReal( SolverParams, 'mmg hgrad', Found)
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hgrad,Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_DPARAMETER <hgrad> Failed')
  END IF

  Pval = ListGetConstReal( SolverParams, 'mmg hgradreq', Found, minv=1.3_dp, maxv=4.0_dp)
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hgradreq,Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_DPARAMETER <hgradreq> Failed')
  END IF


  
  ! If this is a ReTrial then we only change some real valued keywords!
  IF( PRESENT( ReTrial ) ) THEN
    IF( ReTrial ) RETURN
  END IF
    
!!! PARAMS: generic options (debug, mem, verbosity)
  ! [val] Set the verbosity level to n
  !Verbosity = ListGetInteger( SolverParams,'mmg verbosity',Found)
  !IF(.NOT. Found) Verbosity = ListGetInteger( SolverParams,'verbosity',Found)
  !IF (Found) THEN
  !  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_verbose,Verbosity,ierr)
  !  IF ( ierr == 0 ) CALL Fatal('MMGSolver',&
  !       'Call to MMG3D_SET_IPARAMETER Failed')
  !END IF

  ! [val] Set the maximal memory size to n MBytes.
  !MemIncrease = ListGetInteger(SolverParams,'Increase Memory',Found)
  !IF (Found) THEN
  !  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_mem,MemIncrease,ierr)
  !  IF ( ierr == 0 ) CALL Fatal('MMGSolver', &
  !       'Call to MMG3D_SET_IPARAMETER Failed')
  !END IF

  ! [0/1] Turn on the debug mode.
  DebugMode = ListGetLogical(SolverParams,'mmg Debug Mode',Found)
  IF (DebugMode) THEN
    CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_debug,1,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER Failed')
  END IF

!!! OTHER PARAMETERS: NOT ALL TESTED
  Pval = ListGetCReal( SolverParams, 'mmg Angle detection',Found)
  IF (Found) THEN
    CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_angleDetection,&
         Pval,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_DPARAMETER <Angle detection> Failed')
  ENDIF

  ! !< [1/0], Avoid/allow surface modifications */ 
  NoAngleDetect = ListGetLogical(SolverParams,'mmg No Angle detection',Found)
  IF (NoAngleDetect) THEN
    CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle,0,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName, &
         'Call to MMG3D_SET_IPARAMETER <No Angle detection> Failed')
  END IF

  ! [1/0] Avoid/allow point insertion
  !NoInsert = ListGetLogical(SolverParams,'No insert',Found)
  !IF (NoInsert) THEN
  !  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_noinsert,1,ierr)
  !  IF ( ierr == 0 ) CALL Fatal('MMGSolver', &
  !       'Call to MMG3D_SET_IPARAMETER <No insert> Failed') 
  !END IF

  ! [1/0] Avoid/allow edge or face flipping
  !NoSwap = ListGetLogical(SolverParams,'No swap',Found)
  !IF (NoSwap) THEN
  !  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_noswap,1,ierr)
  !  IF ( ierr == 0 ) CALL Fatal('MMGSolver',&
  !       'Call to MMG3D_SET_IPARAMETER <No swap>Failed')
  !END IF

!  [1/0] Avoid/allow point relocation
! NoMove = ListGetLogical(SolverParams,'No move',Found)
! IF (NoMove) THEN
!   CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_nomove,1,ierr)
!   IF ( ierr == 0 ) CALL Fatal('MMGSolver',&
!        'Call to MMG3D_SET_IPARAMETER <No move> Failed')
! END IF

  ! [1/0] Avoid/allow surface modifications
  NoSurf = ListGetLogical(SolverParams,'mmg No surf',Found)
  IF (NoSurf) THEN
    CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_nosurf,0,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <No surf> Failed')
  END IF

  niter = ListGetInteger(SolverParams,'mmg niter',Found ) 
  IF( Found ) THEN
    CALL PMMG_SET_IPARAMETER(pmmgMesh, PMMGPARAM_niter,niter,ierr) 
    IF ( ierr == 0 ) CALL Fatal(FuncName,&
         'Call to MMG3D_SET_IPARAMETER <Niter> Failed')
  END IF

!!!
#else
  CALL Fatal('Set_PMMG_Parameters',&
      'Remeshing utility MMG3D has not been installed')
#endif
END SUBROUTINE Set_PMMG_Parameters
   

SUBROUTINE Get_MMG3D_Mesh(NewMesh, Parallel, FixedNodes, FixedElems, Calving)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  LOGICAL :: Parallel
  LOGICAL, OPTIONAL, ALLOCATABLE :: FixedNodes(:), FixedElems(:)
  LOGICAL, OPTIONAL :: Calving
  !------------------------------------------------------------------------------

#ifdef HAVE_MMG

  TYPE(Element_t),POINTER ::  Element, Element0
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER, ALLOCATABLE :: BC2BodyMap(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent,ied
  INTEGER :: i,ii,kk,t,NoBCs, MinIndex, MaxIndex, maxnodes
  LOGICAL :: Found, Debug, Combine
  INTEGER :: nt0, np0, na0, np1
  INTEGER, POINTER :: Perm(:)
  INTEGER, ALLOCATABLE :: Perm0(:),ElemInd0(:)
  TYPE(Mesh_t), POINTER :: OldMesh
  TYPE(Solver_t), POINTER :: Solver
  CHARACTER(:), ALLOCATABLE :: EquationName                              
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Get_MMG3D_Mesh"

  
  ! Set up a map of BoundaryInfo % Constraint to % BodyID
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
  IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to MMGS_Get_meshSize failed!')
  CALL Info(FuncName,'MMG3D_Get_meshSize done',Level=20)
  IF(NPrisms /= 0) CALL Fatal(FuncName, "Programming Error: MMG3D returns prisms")
  IF(NQuads /= 0) CALL Fatal(FuncName, "Programming Error: MMG3D returns quads")

  maxnodes = 4
  nt0 = 0; np0 = 0; na0 = 0
  Combine = ListGetLogical( CurrentModel % Solver % Values,'Keep unmeshed regions',Found)
  IF( Combine ) THEN
    Solver => CurrentModel % Solver
    OldMesh => Solver % Mesh
    
    CALL Info(FuncName,'Merging with unmeshed regions',Found )
    ! This still related to the old mesh
    Perm => Solver % Variable % Perm
    IF(.NOT. ASSOCIATED(Perm)) THEN
      CALL Fatal(FuncName,'Perm vector not assoicated!?')
    END IF

    EquationName = ListGetString( Solver % Values, 'Equation', Found)                
    nt0 = 0
    DO t=1,OldMesh % NumberOfBulkElements
      Element => OldMesh % Elements(t)
      IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
      nt0 = nt0 +1 
    END DO

    IF(nt0 == 0) THEN
      CALL Info('FuncName','No passive regions, skipping combining of meshes!')
      Combine = .FALSE.
    ELSE
      CALL Info('FuncName','Elements in passive regions: '//I2S(nt0))
    END IF
  END IF

  IF( Combine ) THEN        
    ALLOCATE(Perm0(OldMesh % NumberOfNodes))
    Perm0 = 0

    ! Count saved boundary elements
    na0 = 0
    DO t=1,OldMesh % NumberOfBoundaryElements
      Element => OldMesh % Elements( OldMesh % NumberOfBulkElements + t )
      IF(ALL(Perm(Element % NodeIndexes) > 0) ) CYCLE
      na0 = na0 + 1
      maxnodes = MAX(Element % TYPE % NumberOfNodes,maxnodes)
    END DO

    IF(na0 > 0) THEN
      ALLOCATE(ElemInd0(OldMesh % NumberOfBulkElements))
      ElemInd0 = 0
    END IF

    ! Count saved bulk elements
    nt0 = 0
    DO t=1,OldMesh % NumberOfBulkElements
      Element => OldMesh % Elements(t)
      IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
      nt0 = nt0 + 1
      Perm0(Element % NodeIndexes) = 1
      IF(na0 > 0) ElemInd0(t) = nt0
    END DO

    ! Count saved nodes
    np0 = 0
    np1 = 0
    DO i=1,OldMesh % NumberOfNodes
      IF( Perm0(i) > 0 ) THEN
        IF( Perm(i) == 0 ) THEN
          ! Internal saved node
          np0 = np0 + 1
          Perm0(i) = Nverts + np0
        ELSE
          ! Shared node
          np1 = np1 + 1
          Perm0(i) = -np1
        END IF
      END IF
    END DO
  END IF
    
  NewMesh => AllocateMesh( NTetras+nt0, NTris+na0, NVerts+np0, ParEnv % PEs>1)
  NewMesh % Name = "MMG3D_Output"
  NewMesh % MaxElementNodes = maxnodes
  NewMesh % MaxElementDOFs = maxnodes
  NewMesh % MeshDim = 3
  NewMesh % NumberOfNodes = NVerts+np0
  NewMesh % NumberOfBulkElements = NTetras+nt0
  NewMesh % NumberOfBoundaryElements = NTris+na0

  IF( Combine ) THEN
    ALLOCATE( NewMesh % InvPerm( NewMesh % NumberOfBulkElements ) )
    NewMesh % InvPerm = 0    
  END IF
  
  !PRINT *,'Mesh Counts',NewMesh % NumberOfNodes,NewMesh % NumberOfBulkElements,&
  !    NewMesh % NumberOfBoundaryElements 
  
  IF(PRESENT(FixedNodes)) THEN
    ALLOCATE(FixedNodes(NVerts+np0))
    FixedNodes = .FALSE.
  END IF
  IF(PRESENT(FixedElems)) THEN
    ALLOCATE(FixedElems(NTetras+NTris+na0+nt0))
    FixedNodes = .FALSE.
  END IF

  IF(Parallel) THEN
    NewMesh % ParallelInfo % GInterface = .FALSE.
  END IF

  !! GET NEW VERTICES
  DO ii=1,NVerts
    CALL MMG3D_Get_vertex(mmgMesh,&
         NewMesh % Nodes % x(ii),&
         NewMesh % Nodes % y(ii),&
         NewMesh % Nodes % z(ii),&
         ref,corner,required,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to  MMG3D_Get_vertex failed!')
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
  END DO

  IF( Combine ) THEN
    DO i=1,OldMesh % NumberOfNodes
      ii = Perm0(i)
      IF( ii > Nverts ) THEN
        NewMesh % Nodes % x(ii) = OldMesh % Nodes % x(i)
        NewMesh % Nodes % y(ii) = OldMesh % Nodes % y(i)
        NewMesh % Nodes % z(ii) = OldMesh % Nodes % z(i)
      END IF
    END DO
  END IF
  CALL Info(FuncName,'MMG3D_Get_vertex done',Level=20)

  IF( Combine ) THEN
    BLOCK
      REAL(KIND=dp) :: x0,y0,z0,s2,s2min,s2minmax,s
      INTEGER :: j,jmin
      s2minmax = 0.0_dp
      DO i=1,OldMesh % NumberOfNodes 
        ii = Perm0(i)
        IF( ii < 0 ) THEN
          x0 = OldMesh % Nodes % x(i)
          y0 = OldMesh % Nodes % y(i)
          z0 = OldMesh % Nodes % z(i)          
          s2min = HUGE(s2min)
          jmin = 0
          DO j=1,Nverts
            s2 = (Newmesh % Nodes % x(j) - x0)**2 + (Newmesh % Nodes % y(j) - y0)**2 + &
                (Newmesh % Nodes % z(j) - z0)**2
            IF( s2 < s2min ) THEN
              s2min = s2
              jmin = j
            END IF
          END DO
          s2minmax = MAX(s2minmax, s2min)
          Perm0(i) = jmin
        END IF
      END DO
      !PRINT *,'MaxMin dx:',SQRT(s2minmax)
      IF( SQRT(s2minmax) > 1.0e-8 ) THEN
        PRINT *,'s2minmax:',s2minmax
        CALL Fatal(FuncName,'Could not accurately find the boundary nodes in the new mesh!')
      END IF
    END BLOCK
  END IF
  
  ! Get new tetrahedral elements
  MinIndex = HUGE(MinIndex)
  MaxIndex = 0
  DO ii = 1,NewMesh % NumberOfBulkElements
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
    IF(PRESENT(FixedElems)) FixedElems(ii) = ( required > 0 ) 

    MinIndex = MIN( MinIndex, MINVAL( NodeIndexes(1:4) ) )
    MaxIndex = MAX( MaxIndex, MAXVAL( NodeIndexes(1:4) ) )
  END DO
  IF (DEBUG) PRINT *,'MMG3D_Get_tets DONE'

  IF( Combine ) THEN
    nt0 = 0
    DO t=1,OldMesh % NumberOfBulkElements        
      Element0 => OldMesh % Elements(t)
      IF ( CheckElementEquation( CurrentModel, Element0, EquationName ) ) CYCLE
      nt0 = nt0 +1 

      Element => NewMesh % Elements(Ntetras+nt0)
      Element % TYPE => Element0 % Type
      Element % NDOFs = Element0 % NDOFs
      Element % ElementIndex = Ntetras+nt0
      Element % GElementIndex = Ntetras+nt0
      Element % PartIndex = ParEnv % myPE
      CALL AllocateVector(Element % NodeIndexes, SIZE(Element0 % NodeIndexes))
      Element % NodeIndexes = Perm0(Element0 % NodeIndexes)

      ! This may be needed for efficient interpolation of remaining dofs
      NewMesh % InvPerm(Ntetras+nt0) = t 

      IF( ANY(Element % NodeIndexes < 0 ) ) THEN        
        PRINT *,'Perm0 small:',Perm0(Element0 % NodeIndexes), np0
      END IF
      IF( ANY(Element % NodeIndexes > Nverts + np0 ) ) THEN        
        PRINT *,'Perm0 big:',Perm0(Element0 % NodeIndexes), np0
      END IF
      Element % BodyId = Element0 % BodyId
    END DO
  END IF
  
  CALL Info(FuncName,'MMG3D_Get_tets done',Level=20)

  ! Get new boundary elements
  kk = NewMesh % NumberOfBulkElements
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
    IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to MMG3D_Get_triangle failed!')

    ! using mmgls new boundary (zero contour) is set to 10 by mmg
    IF(Calving .AND. ref == 10) ref = ref + ElmerBCOffset

    ref = ref - ElmerBCOffset
    IF(ref<0) ref = 0

    Allocate(Element % BoundaryInfo)
    Element % BoundaryInfo % Constraint = ref

    IF(ref > 0 .AND. ref <= CurrentModel % NumberOfBCs) THEN
      Element % BodyId = BC2BodyMap(ref)
    END IF

    IF(PRESENT(FixedElems)) FixedElems(kk) = ( required > 0 )
    
    CALL MMG3D_GET_TetFromTria(mmgMesh,ii,parent,ied,ierr)
    
    IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to  MMG3D_Get_TetFromTria failed!')
    Element % BoundaryInfo % Left => NewMesh % Elements(parent) !TODO - parent ID offset?
  END DO
  

  IF( na0 > 0 ) THEN
    kk = NewMesh % NumberOfBulkElements + Ntris
    na0 = 0
    DO t=1,OldMesh % NumberOfBoundaryElements
      Element0 => OldMesh % Elements( OldMesh % NumberOfBulkElements + t )
      IF(ALL(Perm(Element0 % NodeIndexes) > 0) ) CYCLE      

      kk = kk+1
      Element => NewMesh % Elements(kk)        
      Element % TYPE => Element0 % Type
      Element % NDOFs = Element0 % NDOFs
      Element % ElementIndex = kk
      Element % PartIndex = ParEnv % myPE

      CALL AllocateVector(Element % NodeIndexes, SIZE(Element0 % NodeIndexes))
      Element % NodeIndexes = Perm0(Element0 % NodeIndexes)

      ALLOCATE(Element % BoundaryInfo)
      Element % BoundaryInfo % Constraint = Element0 % BoundaryInfo % Constraint

      IF( ASSOCIATED(Element0 % BoundaryInfo % Left) ) THEN
        i = ElemInd0(Element0 % BoundaryInfo % Left % ElementIndex)          
        IF(i>0) Element % BoundaryInfo % Left => NewMesh % Elements(i+nt0)
      END IF
      IF( ASSOCIATED(Element0 % BoundaryInfo % Right) ) THEN
        i = ElemInd0(Element0 % BoundaryInfo % Right % ElementIndex)          
        IF(i>0) Element % BoundaryInfo % Right => NewMesh % Elements(i+nt0)
      END IF
    END DO
  END IF
             
! CALL SetMeshMaxDOFs(NewMesh)

  ! this throws up lots of warning in 3D calving
  ! since calving uses multiple inputs turning off during merge (15/11/23)
  IF(.NOT. Calving) CALL Finalize_MMG_Mesh(NewMesh)

  
#else
     CALL Fatal('Get_MMG3D_Mesh',&
        'Remeshing utility MMG3D has not been installed')
#endif

END SUBROUTINE Get_MMG3D_Mesh

! Subroutine to negotiate new global node numbers between partitions
! as a necessary precursor to repartitioning the mesh.
! Expects to receive OldMesh with valid GlobalDOFs, and new mesh 
! in which all GlobalDOFs are either present in the OldMesh, or set to zero.
! We allow here that each partition may have any number of nodes (inc. zero)
! Assumes that NewMesh doesn't have any nodes which are both *shared* and *unmarked*
!-----------------------------------------------------------------------------------
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
  CHARACTER(*), PARAMETER :: FuncName="RenumberGDOFs"

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

! Because elements aren't shared, renumbering is easier
! Simply number contiguously in each partition
! NOTE: won't work for halo elems
!---------------------------------------------------------
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


! Based on a previous mesh with valid nodal parallelinfo (% GInterface & % NeighbourList)
! map that info onto NewMesh which shares at least some GlobalDOFs. Intended use is
! to enable reconnection of a parallel mesh part which has been remeshed internally, but
! whose partition boundaries remain as they were. e.g. CalvingRemeshMMG.F90
!---------------------------------------------------------------------------------------------
SUBROUTINE MapNewParallelInfo(OldMesh, NewMesh)
  TYPE(Mesh_t), POINTER :: OldMesh, NewMesh
  !---------------------------------
  INTEGER, ALLOCATABLE :: GtoNewLMap(:)
  INTEGER :: i,k,n,MaxNGDof, MinNGDof
  CHARACTER(*), PARAMETER :: FuncName="MapNewBCInfo"
  
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
    IF(OldMesh % ParallelInfo % GInterface(i)) THEN
      k = OldMesh % ParallelInfo % GlobalDOFs(i)
      IF(k < LBOUND(GToNewLMap,1) .OR. k > UBOUND(GToNewLMap,1)) THEN
        CALL Warn(FuncName, "Interface node from OldMesh missing in NewMesh")
        CYCLE
      ELSEIF(GToNewLMap(k) == 0) THEN
        CALL Warn(FuncName, "Interface node from OldMesh missing in NewMesh")
        CYCLE
      END IF
      k = GToNewLMap(k)
      NewMesh % ParallelInfo % GInterface(k) = .TRUE.

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


! A subroutine for a 3D mesh (in serial) with MMG3D
! Inputs:
!   InMesh - the initial mesh
!   Metric - 2D real array specifying target metric
!   NodeFixed, ElemFixed - Optional mask to specify 'required' entities
! Output:
!   OutMesh - the improved mesh
!----------------------------------------------------------------------------------
SUBROUTINE RemeshMMG3D(Model, InMesh,OutMesh,EdgePairs,PairCount,&
    NodeFixed,ElemFixed,Params,Hvar,Solver,Success)

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  TYPE(ValueList_t), POINTER, OPTIONAL :: Params
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount
  TYPE(Variable_t), POINTER, OPTIONAL :: HVar
  TYPE(Solver_t), POINTER, OPTIONAL :: Solver
  LOGICAL :: Success
  !-----------
  TYPE(Mesh_t), POINTER :: WorkMesh
  TYPE(ValueList_t), POINTER :: FuncParams, Material
  TYPE(Variable_t), POINTER :: TimeVar, MMGVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:)
  REAL(KIND=dp), POINTER :: WorkReal(:,:,:) => NULL(), hminarray(:,:) => NULL(),&
       hausdarray(:,:) => NULL()
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols, &
      UseHvar, UseTargetLength, MultipleInputs
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(:), ALLOCATABLE :: FuncName, &
        premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile
  CHARACTER(MAX_NAME_LEN) :: MeshName, SolName
  SAVE :: WorkReal

#ifdef HAVE_MMG

  Debug = .TRUE.
  Parallel = ParEnv % PEs > 1
  FuncName = "RemeshMMG3D"

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop', ThisOnly = .TRUE.)
  IF(.NOT. ASSOCIATED(mmgVar) ) THEN        
    CALL VariableAddVector( Model % Mesh % Variables,Model % Mesh,&
        Name='MMG Loop',Global=.TRUE.)
    mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop' )   
  END IF
  mmgVar % Values(1) = 0.0_dp
  
  
  ! Optionally pass valuelist, by default use the Simulation section
  IF(PRESENT(Params)) THEN
    FuncParams => Params
  ELSE
    i = ListGetInteger( CurrentModel % Bodies(InMesh % Elements(1) % BodyId) % Values, &
        'Material')      
    FuncParams => CurrentModel % Materials(i) % Values  !TODO, this is not generalised
  END IF
  
  MaxRemeshIter = ListGetInteger( FuncParams,'MMG Remesh Max Iterations', Found )
  IF(.NOT. Found ) MaxRemeshIter = 10
  MultipleInputs = ListGetLogical( FuncParams,'MMG Multiple Inputs', Found )
  IF(.NOT. Found) MultipleInputs = .FALSE.
  IF(MultipleInputs) THEN
    hminarray => ListGetConstRealArray(FuncParams, "MMG Hmin", Found)
    hausdarray => ListGetConstRealArray(FuncParams, "MMG Hausd", Found)
  END IF
   
  RemeshMinQuality = ListGetConstReal(FuncParams, "MMG Min Quality", Found, DefValue=0.0001_dp)
  
  SaveMMGMeshes = ListGetLogical(FuncParams,"Save MMG Meshes", Found)
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre MMG Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "MMG Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF
  
  SaveMMGSols = ListGetLogical(FuncParams,"Save MMG Sols", Found)
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre MMG Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "MMG Output Sol Name", UnfoundFatal = .TRUE.)
  END IF

  NNodes = InMesh % NumberOfNodes
  NBulk = InMesh % NumberOfBulkElements
  NBdry = InMesh % NumberOfBoundaryElements

  UseHvar = PRESENT(Hvar)
  IF(UseHvar) UseHvar = ASSOCIATED(HVar)
  
  UseTargetLength = ListCheckPresent( FuncParams,'MMG Target Length' )

  
  IF( UseHvar ) THEN
    CALL Info(FuncName,'Using external field for mesh metric: '//TRIM(HVar % Name),Level=10)
    IF( HVar % Dofs == 1 ) THEN
      SolType = MMG5_Scalar
      AnisoFlag = .FALSE.
    ELSE
      CALL Fatal(FuncName,'Implemented so far only for 1 dofs!')
    END IF      
  ELSE IF( UseTargetLength ) THEN
    AnisoFlag = ListGetLogical(FuncParams, "MMG Anisotropic", Found ) 
    IF(.NOT. Found) AnisoFlag = .TRUE.
    
    IF(AnisoFlag) THEN
      CALL Info(FuncName,'Using anisotropic mesh metric',Level=15)
      WorkMesh => Model % Mesh
      Model % Mesh => InMesh
      
      SolType = MMG5_Tensor
      !Upper triangle of symmetric tensor: 11,12,13,22,23,33
      ALLOCATE(Metric(NNodes,6))
      Metric = 0.0
      DO i=1,NNodes
        NodeNum = i        
        CALL ListGetRealArray(FuncParams,"MMG Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)
        
        !Metric = 1.0/(edge_length**2)
        Metric(i,1) = 1.0 / (WorkReal(1,1,1)**2.0)
        Metric(i,4) = 1.0 / (WorkReal(2,1,1)**2.0)
        Metric(i,6) = 1.0 / (WorkReal(3,1,1)**2.0)
      END DO
            
      Model % Mesh => WorkMesh
      WorkMesh => NULL()
    ELSE
      CALL Info(FuncName,'Using isotropic mesh metric',Level=15)
      SolType = MMG5_Scalar
      ALLOCATE(Metric(NNodes, 1))
      DO i=1,NNodes
        NodeNum = i
        Metric(i,:) = ListGetReal(FuncParams,"MMG Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
      END DO
    END IF
  END IF
    
  nBCs = CurrentModel % NumberOfBCs
  body_offset = nBCs + CurrentModel % NumberOfBodies + 1

  ! The feature where some of the elements are conserved and some are remeshed
  ! does not really like the offset. Instead of fixing it we here only avoid the
  ! problem...
  IF( ListGetLogical( FuncParams,'Keep Unmeshed regions',Found ) ) THEN
    body_offset = 0
  END IF
  
  IF( body_offset > 0 ) THEN
    DO i=1,InMesh % NumberOfBulkElements
      InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID + body_offset
    END DO
  END IF
    
  DO mmgloops = 1, MaxRemeshIter 

    CALL Info(FuncName,'Applying remeshing trial: '//I2S(mmgloops),Level=5)

    Success = .TRUE.
    IF( mmgloops > 1 ) THEN
      !! Redoing adaptive mesh, release the previous mmg mesh
      CALL MMG3D_Free_all(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, MMG5_ARG_end)      
    END IF
        
    ! Enable external depende on "mmg loop"
    mmgVar % Values(1) = 1.0_dp * mmgloops 
    
    mmgMesh = 0
    mmgSol  = 0
    CALL MMG3D_Init_mesh(MMG5_ARG_start, &
        MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
        MMG5_ARG_end)

    IF(MultipleInputs) THEN
      CALL ListAddConstReal(FuncParams, 'adaptive min h', hminarray(mmgloops, 1))
      CALL ListAddConstReal(FuncParams, 'adaptive hausd', hausdarray(mmgloops, 1))
    END IF
           
    ! If this is retrial then get only selected parameters again that may depend on the variable "MMG Loop".
    CALL Set_MMG3D_Parameters(FuncParams, mmgloops > 1 )
            
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
    IF (PRESENT(PairCount)) THEN
      CALL Info(FuncName,'Using provided edge pairs: '//I2S(PairCount),Level=10)
      CALL SET_MMG3D_MESH(InMesh,Parallel,EdgePairs,PairCount,Solver=Solver)
    ELSE
      CALL SET_MMG3D_MESH(InMesh,Parallel,Solver=Solver)
    END IF
    
    ! Set the metric values at nodes
    IF( UseTargetLength .OR. UseHVar ) THEN
      CALL MMG3D_Set_SolSize(mmgMesh, mmgSol, MMG5_Vertex, NNodes, SolType,ierr)
      IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")

      DO i=1,NNodes
        IF(AnisoFlag) THEN
          CALL MMG3D_Set_TensorSol(mmgSol,Metric(i,1),Metric(i,2),Metric(i,3),&
              Metric(i,4),Metric(i,5),Metric(i,6),i,ierr)
          IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set tensor solution at vertex")
        ELSE 
          IF( UseHvar ) THEN
            CALL MMG3D_Set_ScalarSol(mmgSol,HVar % Values(i),i,ierr)
          ELSE
            CALL MMG3D_Set_ScalarSol(mmgSol,Metric(i,1),i,ierr)
          END IF
          IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
        END IF
      END DO
    END IF
      
    !Turn on debug (1)
    !CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMGPARAM_debug, 1,ierr)
    
    !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmin,hmin,ierr)
    !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMGPARAM_hmax,hmax,ierr)
    !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hausd,hausd,ierr)
    !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad,hgrad,ierr)
    ! allow surface modifications
    !CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMGPARAM_nosurf,0,ierr)
    
    !Turn off sharp angle detection (0)   
    !CALL MMG3D_SET_IPARAMETER(mmgMesh,mmgSol,MMG3D_IPARAM_angle,0,ierr)
    !Option to set angle detection threshold:
    !CALL MMG3D_SET_DPARAMETER(mmgMesh,mmgSol,MMG3D_DPARAM_angleDetection,&
    !      85.0_dp,ierr)
    
    !Take care of fixed nodes/elements if requested
    IF(PRESENT(NodeFixed)) THEN
      IF( InfoActive(20) ) THEN
        j = COUNT(NodeFixed)
        CALL Info(FuncName,'Set required vertexes: '//I2S(j))
      END IF
      DO i=1,NNodes
        IF(NodeFixed(i)) THEN
          CALL MMG3D_SET_REQUIREDVERTEX(mmgMesh,i,ierr)
        END IF
      END DO
    END IF

    IF(PRESENT(ElemFixed)) THEN
      IF( InfoActive(20) ) THEN
        j = COUNT(ElemFixed)
        CALL Info(FuncName,'Set reuired elements: '//I2S(j))
      END IF
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
    
    CALL MMG3D_mmg3dlib(mmgMesh,mmgSol,ierr)
    
    IF ( ierr == MMG5_STRONGFAILURE .OR. ierr == MMG5_LOWFAILURE ) THEN
      CALL Warn(FuncName,'MMG3DLib resulted to error, trying remeshing')
      Success = .FALSE.
      CYCLE
    END IF
    
    CALL Info(FuncName,'MMG3D_mmg3dlib done',Level=20)
    
    IF(SaveMMGMeshes) THEN
      WRITE(MeshName, '(A,I0,A)') TRIM(mmg_meshfile), time, '.mesh'
      CALL MMG3D_SaveMesh(mmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
    END IF
    IF(SaveMMGSols) THEN
      WRITE(SolName, '(A,I0,A)') TRIM(mmg_solfile), time, '.sol'
      CALL MMG3D_SaveSol(mmgMesh, mmgSol,SolName,LEN(TRIM(SolName)),ierr)
    END IF
    
    CALL MMG3D_Get_meshSize(mmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

    counter=0
    DO i=1, NTetras
      CALL MMG3D_Get_TetrahedronQuality(mmgMesh, mmgSol, i, Quality)
      IF(Quality == 0) CALL Warn(FuncName, 'Remeshing could not determine elem quality')
      IF(Quality <= RemeshMinQuality) counter = counter+1
    END DO
    
    IF ( Counter > 0 ) THEN
      CALL Info(FuncName,'Bad element count: '//TRIM(I2S(counter)),Level=20)
      CALL Warn(FuncName,'Bad elements detected - rerunning remeshing')
      Success = .FALSE.
      CYCLE
    END IF
    
    IF(Success) EXIT
  END DO
  

  ! Transfer the new mesh into Elmer mesh format
  CALL GET_MMG3D_MESH(OutMesh,Parallel, Calving=MultipleInputs)

  ! Release the mesh from mmg mesh format
  CALL MMG3D_Free_all(MMG5_ARG_start, &
  MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
  MMG5_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements

  ! Reset the BodyIDs (see above)
  IF( body_offset > 0 ) THEN 
    OutMesh % Elements(1:Nbulk) % BodyID = OutMesh % Elements(1:Nbulk) % BodyID - body_offset

    i=InMesh % NumberOfBulkElements
    InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID - body_offset
  END IF
    
  ! And delete the unneeded BC elems
  ! Important note for calving. mmg adds boundary nodes to upstream user defined boundary.
  ! these need to be removed
  ! this is a temp fix as calving algo has mutliple inputs
  IF(MultipleInputs) THEN
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

  RETURN

    
20 CONTINUE

  ! if remeshing has failed need to reset body ids
  IF(.NOT. Success .AND. body_offset > 0 ) THEN
    i=InMesh % NumberOfBulkElements
    InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID - body_offset
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
SUBROUTINE SequentialRemeshParMMG(Model, InMesh,OutMesh,Boss,EdgePairs,PairCount,&
    NodeFixed,ElemFixed,Params)

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
  TYPE(Variable_t), POINTER :: TimeVar, MMGVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:),hminarray(:),hausdarray(:)
  REAL(KIND=dp), POINTER :: WorkReal(:,:,:) => NULL() 
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(:), ALLOCATABLE :: MeshName, SolName, &
        premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile
  CHARACTER(*), PARAMETER :: FuncName = "SequentialRemeshParMMG3D"
  SAVE :: WorkReal 

#ifdef HAVE_PARMMG

  Debug = .TRUE.
  Parallel = ( ParEnv % PEs > 1 )

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop', ThisOnly = .TRUE.)
  IF(.NOT. ASSOCIATED(mmgVar) ) THEN        
    CALL VariableAddVector( Model % Mesh % Variables,Model % Mesh,&
        Name='MMG Loop',Global=.TRUE.)
    mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop' )   
  END IF
  mmgVar % Values(1) = 0.0_dp

  ! params must be provided as not a global mesh
  FuncParams => Params
  
  MaxRemeshIter = ListGetInteger( FuncParams,'MMG Remesh Max Iterations', Found )
  IF(.NOT. Found ) MaxRemeshIter = 10
    
  !Get parameters from valuelist
  !hausd, hmin, hmax, hgrad, anisoflag, the metric
  !Scalar, vector, tensor metric?


  RemeshMinQuality = ListGetConstReal(FuncParams, "MMG Min Quality",Found, DefValue=0.0001_dp)
  AnisoFlag = ListGetLogical(FuncParams, "MMG Anisotropic", DefValue=.TRUE.)

  SaveMMGMeshes = ListGetLogical(FuncParams,"Save RemeshMMG3D Meshes", DefValue=.FALSE.)
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre RemeshMMG3D Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "MMG Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF

  SaveMMGSols = ListGetLogical(FuncParams,"Save RemeshMMG3D Sols", DefValue=.FALSE.)
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre RemeshMMG3D Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "MMG Output Sol Name", UnfoundFatal = .TRUE.)
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
        CALL ListGetRealArray(FuncParams,"MMG Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)

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
        Metric(i,:) = ListGetReal(FuncParams,"MMG Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
      END DO

    END IF

    body_offset = CurrentModel % NumberOfBCs + CurrentModel % NumberOfBodies + 1
    body_offset = 0

    nBCs = CurrentModel % NumberOfBCs
    IF( body_offset > 0 ) THEN
      i=InMesh % NumberOfBulkElements
      InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID + body_offset
    END IF
      

    DO mmgloops = 1, MaxRemeshIter 
      
      Success = .TRUE.
      MMGVar % Values(1) = 1.0_dp * mmgloops
      
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
      CALL Info(FuncName,'Initiating parmmg mesh',Level=20)
      
      pmmgMesh = 0
      CALL PMMG_Init_parMesh(PMMG_ARG_start, &
          PMMG_ARG_ppParMesh,pmmgMesh, PMMG_ARG_pMesh,PMMG_ARG_pMet, &
          PMMG_ARG_dim,%val(3),PMMG_ARG_MPIComm,%val(ELMER_COMM_WORLD), &
          PMMG_ARG_end)
      
      ! If this is retrial then get only selected parameters again that may depend on the variable "MMG Loop".
      CALL Set_PMMG_Parameters(FuncParams, mmgloops > 1 )
        
      !hmin = hminarray(mmgloops)
      !Hausd = hausdarray(mmgloops)
      
      !WRITE(Message, '(A,F10.5,A,F10.5)') 'Applying levelset with Hmin ',Hmin, ' and Hausd ', Hausd
      !CALL INFO(FuncName, Message)
            
      IF(Boss) THEN
        CALL Info(FuncName,'Setting mesh ....',Level=20)
        IF (PRESENT(PairCount)) THEN
          CALL SET_ParMMG_MESH(InMesh,Parallel,EdgePairs,PairCount)
        ELSE
          CALL SET_ParMMG_MESH(InMesh,Parallel)
        END IF
        
        CALL Info(FuncName,'Setting met size',Level=20)
        !Set the metric values at nodes
        CALL PMMG_Set_MetSize(pmmgMesh, MMG5_Vertex, NNodes, SolType,ierr)
        IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set solution size.")
        
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
      CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_debug,1,ierr)   
      !CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmin,hmin,ierr)
      !CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hmax,hmax,ierr)
      !CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hausd,hausd,ierr)
      !CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_hgrad,hgrad,ierr)
      ! allow surface modifications
      CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_nosurf,0,ierr)
      
      !Turn off sharp angle detection (0)   
      CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle,0,ierr)
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
        IF( SaveMMGMeshes .OR. SaveMMGSols ) THEN
          CALL Info(FuncName,'Saving of PMMG files finished', Level=20)
        END IF
      END IF

      CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
      CALL PMMG_parmmglib_centralized(pmmgMesh,ierr)
      IF ( ierr == PMMG_STRONGFAILURE .OR. ierr == PMMG_LOWFAILURE ) THEN
        CALL Warn(FuncName,'BAD ENDING OF PMMGLIB: UNABLE TO SAVE MESH')
        !! Release mmg mesh
        CALL PMMG_Free_all ( PMMG_ARG_start,     &
            PMMG_ARG_ppParMesh,pmmgMesh,         &
            PMMG_ARG_end)
        Success = .FALSE.
      END IF

      IF(Success) EXIT
      
    END DO

  END IF
    
  CALL Info(FuncName,'PMMG_parmmglib_centralized done',Level=20)

  IF(Boss) THEN
    IF(SaveMMGMeshes) THEN
      MeshName = TRIM(mmg_meshfile)//I2S(time)//'.mesh'
      CALL PMMG_SaveMesh_Centralized(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
    END IF
    IF(SaveMMGSols) THEN
      SolName = TRIM(mmg_solfile)//I2S(time)//'.sol'
      CALL PMMG_SaveMet_Centralized(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
    END IF

    CALL PMMG_Get_meshSize(pmmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

    !! GET THE NEW MESH
    CALL Info(FuncName,'Recovering mesh',Level=20)
    CALL GET_ParMMG_MESH(OutMesh,Parallel,Calving=.FALSE.)
  END IF

  CALL Info(FuncName,'Releasing mesh',Level=20)
  !! Release mmg mesh
  CALL PMMG_Free_all ( PMMG_ARG_start,     &
       PMMG_ARG_ppParMesh,pmmgMesh,         &
       PMMG_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements

  IF(Boss) THEN
    !Reset the BodyIDs (see above)
    IF( body_offset > 0 ) THEN
      OutMesh % Elements(1:Nbulk) % BodyID = OutMesh % Elements(1:Nbulk) % BodyID - body_offset
    END IF

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
    IF(.NOT. Success .AND. body_offset > 0 ) THEN
      i=InMesh % NumberOfBulkElements
      InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID - body_offset
    END IF
  END IF

#else
  CALL Fatal(FuncName, "Remeshing utility PMMG has not been installed")
#endif

END SUBROUTINE SequentialRemeshParMMG


SUBROUTINE Set_ParMMG_Mesh(Mesh, Parallel, EdgePairs, PairCount, FreezeInternalArg)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel
  LOGICAL, OPTIONAL :: FreezeInternalArg
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount

#ifdef HAVE_PARMMG
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER :: i,j,k, NNodes,NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ref,ierr,&
      NoNeighbours, counter, MaxNeighbours
  INTEGER, ALLOCATABLE :: NodeRefs(:), NeighbourList(:), NSharedNodes(:), SharedNodes(:,:),&
      SharedNodesGlobal(:,:)
  LOGICAL, ALLOCATABLE :: IsNeighbour(:)
  INTEGER, POINTER :: Neighbours(:)
  LOGICAL :: Warn101=.FALSE., Warn202=.FALSE.,Debug=.FALSE.,Elem202, FreezeInternal
  CHARACTER(*), PARAMETER :: FuncName="Set_ParMMG_Mesh"
  IF(CoordinateSystemDimension() /= 3) CALL Fatal("ParMMG","Only works for 3D meshes!")

  ALLOCATE(NodeRefs(6))

  CALL Info(FuncName,'Changing Elmer mesh format into MMG mesh format in parallel',Level=10)
  
  IF(Parallel) CALL Assert(ASSOCIATED(Mesh % ParallelInfo % GlobalDOFs), FuncName,&
       "Parallel sim but no ParallelInfo % GlobalDOFs")

  FreezeInternal = .FALSE.
  IF(PRESENT(FreezeInternalArg)) FreezeInternal = FreezeInternalArg

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

    IF(.NOT.FreezeInternal) THEN
      IF(i>nbulk) THEN
         IF(Element % BoundaryInfo % Constraint<=0) CYCLE
      END IF
    END IF

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
      CALL Fatal(FuncName,"can't handle pyramid elements (605)")
    CASE(706)
      NPrisms = NPrisms + 1
    CASE(808)
      CALL Fatal(FuncName,"can't handle brick/hexahedral elements (808)")
    CASE DEFAULT
      CALL Fatal(FuncName,"Unsupported element type: "&
          //TRIM(I2S(Element % TYPE % ElementCode)))
    END SELECT
  END DO

  IF(Warn101) CALL Warn(FuncName,"101 elements detected - these won't be remeshed")
  IF(Warn202) CALL Warn(FuncName,"202 elements detected - these won't be remeshed")

  !args: mesh, nvertices, ntetra, nprisms, ntriangless, nquads, nedges
  CALL PMMG_Set_meshSize(pmmgMesh,nverts,ntetras,nprisms,ntris,nquads,nedges,ierr)

  IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to PMMG_Set_meshSize failed!')
  CALL Info(FuncName,'PMMG_Set_meshSize done',Level=20)

  ref = 0
  DO i=1,NVerts
    ! ref = GDOF + 10 to avoid an input ref of 10 being confused
    ! with mmg output ref of 10 which occurs on some new nodes
    ! GDOF = ref - 10
    IF(Parallel) ref = Mesh % ParallelInfo % GlobalDOFs(i) + 10
    CALL PMMG_Set_vertex(pmmgMesh, Mesh%Nodes%x(i), &
         Mesh%Nodes%y(i),Mesh%Nodes%z(i), ref, i, ierr)
    IF ( ierr /= 1 ) CALL Fatal(FuncName,'Call to PMMG_Set_vertex failed!')
  END DO
  CALL Info(FuncName,'PMMG_Set_vertex done',Level=20)

  ntetras = 0
  nprisms = 0
  ntris = 0
  nquads = 0
  nedges = 0

  !Cycle mesh elements, sending them to MMG3D
  DO i=1,nbulk+nbdry
    Element => Mesh % Elements(i)
    IF(.NOT.FreezeInternal) THEN
      IF(i>nbulk) THEN
         IF(Element % BoundaryInfo % Constraint<=0) CYCLE
      END IF
    END IF

    NNodes = Element % TYPE % NumberOfNodes

    NodeIndexes => Element % NodeIndexes
    NodeRefs(1:NNodes) = NodeIndexes(1:NNodes)

    SELECT CASE(Element % TYPE % ElementCode)
    CASE(303)
      ntris = ntris + 1
      CALL PMMG_Set_triangle(pmmgMesh,  NodeRefs(1), NodeRefs(2), NodeRefs(3), &
           Element % BoundaryInfo % Constraint+ElmerBCOffset, ntris, ierr)
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
  CALL Info(FuncName,'Set elements done',Level=20)

  !! use element pairs '202' elements
  Elem202 = (PRESENT(EdgePairs))
  IF (Elem202) THEN
    DO i=1, PairCount
      NEdges = NEdges + 1
      CALL PMMG_Set_edge(pmmgMesh, EdgePairs(1,i), EdgePairs(2,i), 1, nedges, ierr)
      CALL PMMG_Set_ridge(pmmgMesh, nedges, ierr)
      !CALL MMG3D_Set_requiredEdge(mmgMesh, Nedges, ierr)
    END DO
    CALL Info(FuncName,'ParMMG - Set edge elements done')
  END IF
 
  ! use nodes to set mpi comms
  CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_APImode, 1, ierr)

  ! set the number of neighbours
  ALLOCATE(IsNeighbour(ParEnv % PEs))
  IsNeighbour = .FALSE.
  MaxNeighbours = 0
  DO i=1, Mesh % NumberOfNodes
    Neighbours => Mesh %  ParallelInfo % NeighbourList(i) % Neighbours
    k = SIZE(Neighbours) 
    MaxNeighbours = MAX(MaxNeighbours,k) 
    DO j=1, k
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

  IF( InfoActive(30) ) THEN
    PRINT *,'NoNeighbours:',NoNeighbours, MaxNeighbours
    PRINT *,'NeigbourList:',ParEnv % MyPe, NeighbourList
    PRINT *,'SharedNodes:',ParEnv % MyPe, NSharedNodes
  END IF
    
  IF( SUM( NSharedNodes ) == 0 ) THEN
    CALL Fatal(FuncName,'No shared nodes?!')
  END IF
  
  
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
     CALL Fatal('Set_ParMMG_Mesh',&
        'Remeshing utility ParMMG has not been installed')
#endif

END SUBROUTINE Set_ParMMG_Mesh


SUBROUTINE Get_ParMMG_Mesh(NewMesh, Parallel, FixedNodes, FixedElems, Calving)

  !------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: NewMesh
  LOGICAL :: Parallel
  LOGICAL :: Calving
  LOGICAL, OPTIONAL, ALLOCATABLE :: FixedNodes(:), FixedElems(:)
  !------------------------------------------------------------------------------

#ifdef HAVE_PARMMG

  TYPE(Element_t),POINTER ::  Element
  INTEGER, POINTER :: NodeIndexes(:), Neighbours(:), NSharedNodes(:), SharedNodes(:),&
          GlobalNums(:), GlobalNums2(:), NodeNeigh(:), NodeNeigh0(:)
  INTEGER, ALLOCATABLE :: BC2BodyMap(:)
  INTEGER :: NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, nbulk, nbdry,ierr, NoNeighbours,&
          OutProc, counter, GlobalID, owner, unique, ntot
  INTEGER :: ref,corner,required,ridge
  INTEGER :: parent(2),ied(2)
  INTEGER :: i,j,k,ii,kk,NoBCs,MinIndex,MaxIndex,MinIndexBC,MaxIndexBC,imin
  LOGICAL :: Found, Debug
  LOGICAL, ALLOCATABLE :: UsedNode(:)
  CHARACTER(*), PARAMETER :: FuncName="Get_ParMMG_Mesh"

  CALL Info(FuncName,'Getting Elmer mesh format from MMG mesh format in parallel',Level=10)
  
  Debug= .FALSE.

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
  IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to PMMG_Get_meshSize failed!')
  CALL Info(FuncName,'PMMG_Get_meshSize done',Level=20)

  ! INITIALISE THE NEW MESH STRUCTURE
  !NPrisms, NQuads should be zero

  NewMesh => AllocateMesh( NTetras, NTris, NVerts, .TRUE.)

  NewMesh % Name = "ParMMG_Output"
  NewMesh % MaxElementNodes = 4
  NewMesh % MaxElementDOFs = 4
  NewMesh % MeshDim = 3

  IF(PRESENT(FixedNodes)) THEN
    ALLOCATE(FixedNodes(NVerts))
    FixedNodes = .FALSE.
  END IF
  IF(PRESENT(FixedElems)) THEN
    ALLOCATE(FixedElems(NTetras+NTris))
    FixedNodes = .FALSE.
  END IF

  IF(NPrisms /= 0) CALL Fatal(FuncName, "Programming Error: ParMMG returns prisms")
  IF(NQuads /= 0) CALL Fatal(FuncName, "Programming Error: ParMMG returns quads")
  NewMesh % NumberOfNodes = NVerts
  NewMesh % NumberOfBulkElements = NTetras
  NewMesh % NumberOfBoundaryElements = NTris

  IF(Parallel) THEN
    NewMesh % ParallelInfo % GInterface = .FALSE.
  END IF

  !! GET NEW VERTICES
  NewMesh % Nodes % z = 0._dp
  DO ii=1,NVerts
    CALL PMMG_Get_vertex(pmmgMesh,&
         NewMesh % Nodes % x(ii),&
         NewMesh % Nodes % y(ii),&
         NewMesh % Nodes % z(ii),&
         ref,corner,required,ierr)
    IF ( ierr == 0 ) CALL Fatal(FuncName,'Call to  ParMMG_Get_vertex failed!')
    IF(PRESENT(FixedNodes)) FixedNodes(ii) = required > 0 .AND. ref > 10
  END DO
  ALLOCATE( UsedNode(NVerts) )
  UsedNode = .FALSE.
  
  
  CALL Info(FuncName,'ParMMG_Get_vertex done',Level=20)

  !! GET NEW TETRAS
  MinIndex = HUGE(MinIndex)
  MaxIndex = 0
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

    MinIndex = MIN(MinIndex,MINVAL(NodeIndexes(1:4)))
    MaxIndex = MAX(MaxIndex,MAXVAL(NodeIndexes(1:4)))
    UsedNode(NodeIndexes(1:4)) = .TRUE.
  END DO
  CALL Info(FuncName,'ParMMG_Get_tets done',Level=20)

  !PRINT *,'UsedNodes:',ParEnv % MyPe, COUNT(UsedNode)
  
  !! Get BC Elements
  MinIndexBC = HUGE(MinIndexBC)
  MaxIndexBC = 0
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
    IF ( ierr /= 1 ) CALL Fatal('ParMMGSolver',&
         'Call to  PMMG_Get_triangle failed!')

    ref = ref - ElmerBCOffset
    if ( ref < 0 ) ref=0

    Allocate(Element % BoundaryInfo)
    Element % BoundaryInfo % Constraint=ref

    IF(ref > 0 .AND. ref <= CurrentModel % NumberOfBCs) THEN
      Element % BodyId = BC2BodyMap(ref)
    END IF

    IF(PRESENT(FixedElems)) FixedElems(kk) = required > 0

    MinIndexBC = MIN(MinIndexBC,MINVAL(NodeIndexes(1:3)))
    MaxIndexBC = MAX(MaxIndexBC,MAXVAL(NodeIndexes(1:3)))
    UsedNode(NodeIndexes(1:3)) = .TRUE.
  END DO

  CALL Info(FuncName,'ParMMG_Get_triangle done',Level=20)

  ! Currently segfaults as no parmmg routine to recover parent id
#if 1
  kk=NewMesh % NumberOfBulkElements
  DO ii=1,NTris
    kk = kk + 1
    CALL PMMG_GET_TetsFromTria(pmmgMesh,ii,Parent,ied,ierr)

    IF ( ierr /= 1 ) CALL Fatal('MMGSolver','Call to  MMG3D_Get_TetFromTria failed!')

    Element => NewMesh % Elements(kk)
    Element % BoundaryInfo % Left  => Null()
    Element % BoundaryInfo % Right => Null()
 
    IF(Parent(1)>0.AND.Parent(1)<=kk) THEN
      Element % BoundaryInfo % Left  => NewMesh % Elements(parent(1))
      IF(Parent(2)>0.AND.Parent(2)<=kk) THEN
        Element % BoundaryInfo % Right => NewMesh % Elements(parent(2))
      END IF
    ELSE IF(Parent(2)>0.AND.Parent(2)<=kk) THEN
      Element % BoundaryInfo % Left  => NewMesh % Elements(parent(2))
    END IF
  END DO
#else
  BLOCK
    INTEGER :: i,j,l,k,p

    PRINT*,ParEnv % MyPE, ' start '
    kk = NewMesh % NumberOfBulkElements
    DO ii=1,NTris
      kk =  kk + 1

      Element => NewMesh % Elements(kk)
      IF(.NOT.ASSOCIATED(Element % BoundaryInfo)) STOP

      Element % BoundaryInfo % Left  => Null()
      Element % BoundaryInfo % Right => Null()

      DO k=1,Newmesh % NumberOfBulkElements
        l = 0 
        DO j=1,4
          IF (Element % NodeIndexes(1)==NewMesh % Elements(k) % NodeIndexes(j)) l = l+1
        END DO
        IF(l<1) EXIT

        DO j=1,4
          IF (Element % NodeIndexes(2)==NewMesh % Elements(k) % NodeIndexes(j)) l = l+1
        END DO
        IF(l<2) EXIT

        DO j=1,4
          IF (Element % NodeIndexes(3)==NewMesh % Elements(k) % NodeIndexes(j)) l = l+1
        END DO

        IF (l==3) THEN
          IF ( ASSOCIATED(Element % BoundaryInfo % Left ) ) THEN
            Element % BoundaryInfo % Right => NewMesh % Elements(k)
            EXIT
          ELSE
            Element % BoundaryInfo % Left  => NewMesh % Elements(k)
!           EXIT
          END IF
        END IF
      END DO
    END DO
    PRINT*,ParEnv % MyPE, ' done '
  END BLOCK
#endif

  ! get parallel info back
  ! get number of neighbours
  CALL PMMG_GET_NUMBEROFNODECOMMUNICATORS(pmmgMesh,NoNeighbours,ierr)
     
  ALLOCATE(Neighbours(NoNeighbours), NSharedNodes(NoNeighbours))
  DO i=1, NoNeighbours
    OutProc = i-1
    CALL PMMG_Get_ithNodeCommunicatorSize(pmmgMesh, outProc, Neighbours(i), NSharedNodes(i), ierr)
  END DO

  IF( InfoActive(30) ) THEN
    PRINT *,'NoNeighbours:',ParEnv % MyPe, NoNeighbours
    PRINT *,'Neighbours:',ParEnv % MyPe, Neighbours
    PRINT *,'Neighbours:',ParEnv % MyPe, NSharedNodes    
  END IF

  IF(.NOT. ASSOCIATED(NewMesh % ParallelInfo % GInterface)) &
      ALLOCATE(NewMesh % ParallelInfo % GInterface(NewMesh % NumberOfNodes))
  NewMesh % ParallelInfo % GInterface = .FALSE.
  IF(.NOT. ASSOCIATED(NewMesh % ParallelInfo % NeighbourList)) &
      ALLOCATE(NewMesh % ParallelInfo % NeighbourList(NewMesh % NumberOfNodes))
  ALLOCATE(NodeNeigh0(ParEnv % PEs))
  NodeNeigh0 = 0
  
  ALLOCATE(SharedNodes(MAXVAL(NSharedNodes)))
  DO i=1, NoNeighbours
    OutProc = i-1
    CALL PMMG_Get_ithNodeCommunicator_nodes(pmmgMesh, OutProc, SharedNodes(1:NSharedNodes(i)), ierr)
    
    DO j=1,NSharedNodes(i)
      k = SharedNodes(j)
      counter = 0
      NULLIFY(NodeNeigh)
      NodeNeigh => NewMesh % ParallelInfo % NeighbourList(k) % Neighbours
      IF(ASSOCIATED(NodeNeigh)) THEN
        counter = SIZE(NodeNeigh)
        NodeNeigh0(1:counter) = NodeNeigh(1:counter)
        DEALLOCATE(NodeNeigh)
        ALLOCATE(NodeNeigh(counter+1))
        NodeNeigh(1:counter) = NodeNeigh0(1:counter)
      ELSE
        counter = 1 
        ALLOCATE(NodeNeigh(counter+1))
        NodeNeigh(counter) = ParEnv % MyPe
      END IF

      NodeNeigh(counter+1) = Neighbours(i)
      NewMesh % ParallelInfo % NeighbourList(k) % Neighbours => NodeNeigh       
      NewMesh % ParallelInfo % GInterface(k) = .TRUE.
    END DO
  END DO

  DO k=1,NewMesh % NumberOfNodes 
    NodeNeigh => NewMesh % ParallelInfo % NeighbourList(k) % Neighbours
    IF(ASSOCIATED(NodeNeigh)) THEN
      counter = SIZE(NodeNeigh)

      ! Choose the owner to be the smallest partition index
      imin = 1
      DO i=2,counter
        IF( NodeNeigh(i) < NodeNeigh(imin) ) THEN
          imin = i
        END IF
      END DO
      IF(imin /= 1) THEN
        j = NodeNeigh(1)
        NodeNeigh(1) = NodeNeigh(imin) 
        NodeNeigh(imin) = j
      END IF
    ELSE
      ALLOCATE(NodeNeigh(1))
      NodeNeigh(1) = ParEnv % MyPe
      NewMesh % ParallelInfo % NeighbourList(k) % Neighbours => NodeNeigh
    END IF
  END DO
  DEALLOCATE(NodeNeigh0)
  
  DO ii=1, NewMesh % NumberOfNodes
    CALL PMMG_Get_VertexGloNum(pmmgMesh, GlobalID, owner, ierr)
    NewMesh % ParallelInfo % GlobalDOFs(ii) = GlobalID
  END DO

  CALL Info(FuncName,'Before comm barrier',Level=20)
  
  CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

  CALL Finalize_MMG_Mesh(NewMesh)

#else
     CALL Fatal('Get_ParMMG_Mesh',&
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
SUBROUTINE DistributedRemeshParMMG(Model, InMesh,OutMesh,EdgePairs,PairCount,&
    NodeFixed,ElemFixed,Params,HVar,Angle)

  TYPE(Model_t) :: Model
  TYPE(Mesh_t), POINTER :: InMesh, OutMesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL, ALLOCATABLE, OPTIONAL :: NodeFixed(:), ElemFixed(:)
  INTEGER, ALLOCATABLE, OPTIONAL :: EdgePairs(:,:)
  INTEGER, OPTIONAL :: PairCount  
  REAL(KIND=dp), OPTIONAL :: Angle
  TYPE(Variable_t), POINTER, OPTIONAL :: Hvar
  LOGICAL :: Success
  !-----------
  TYPE(Mesh_t), POINTER :: WorkMesh
  TYPE(ValueList_t), POINTER :: FuncParams, Material
  TYPE(Variable_t), POINTER :: TimeVar, MMGVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: TargetLength(:,:), Metric(:,:),hminarray(:),hausdarray(:)
  REAL(KIND=dp), POINTER :: WorkReal(:,:,:) => NULL()
  REAL(KIND=dp) :: hsiz(3),hmin,hmax,hgrad,hausd,RemeshMinQuality,Quality, TargetX
  INTEGER :: i,j,MetricDim,NNodes,NBulk,NBdry,ierr,SolType,body_offset,&
       nBCs,NodeNum(1), MaxRemeshIter, mmgloops, ElemBodyID, &
       NVerts, NTetras, NPrisms, NTris, NQuads, NEdges, Counter, Time
  INTEGER, ALLOCATABLE :: TetraQuality(:)
  LOGICAL :: Debug, Parallel, AnisoFlag, Found, SaveMMGMeshes, SaveMMGSols, FreezeInternal
  LOGICAL, ALLOCATABLE :: RmElement(:)
  CHARACTER(:), ALLOCATABLE :: MeshName, SolName, &
        premmg_meshfile, mmg_meshfile, premmg_solfile, mmg_solfile
  CHARACTER(*), PARAMETER :: FuncName = "DistributedRemeshParMMG3D"
  SAVE :: WorkReal 

#ifdef HAVE_PARMMG

  Debug = .FALSE.
  Parallel = ParEnv % PEs > 1

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  Time = INT(TimeVar % Values(1))

  mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop', ThisOnly = .TRUE.)
  IF(.NOT. ASSOCIATED(mmgVar) ) THEN        
    CALL VariableAddVector( Model % Mesh % Variables,Model % Mesh,&
        Name='MMG Loop',Global=.TRUE.)
    mmgVar => VariableGet( Model % Mesh % Variables,'MMG Loop' )   
  END IF
  mmgVar % Values(1) = 0.0_dp

  ! params must be provided as not a global mesh
  IF(.NOT. ASSOCIATED( Params ) ) THEN
    CALL Fatal(FuncName,'"Params" not associated!')
  END IF
  FuncParams => Params
  
  MaxRemeshIter = ListGetInteger( FuncParams,'MMG Remesh Max Iterations', Found )
  IF(.NOT. Found ) MaxRemeshIter = 10
    
  RemeshMinQuality = ListGetConstReal(FuncParams, "MMG Min Quality",Found, DefValue=0.0001_dp)

  SaveMMGMeshes = ListGetLogical(FuncParams,"Save RemeshMMG3D Meshes", Found )
  IF(SaveMMGMeshes) THEN
    premmg_meshfile = ListGetString(FuncParams, "Pre RemeshMMG3D Mesh Name", UnfoundFatal = .TRUE.)
    mmg_meshfile = ListGetString(FuncParams, "MMG Output Mesh Name", UnfoundFatal = .TRUE.)
  END IF

  SaveMMGSols = ListGetLogical(FuncParams,"Save RemeshMMG3D Sols", Found ) 
  IF(SaveMMGSols) THEN
    premmg_solfile = ListGetString(FuncParams, "Pre RemeshMMG3D Sol Name", UnfoundFatal = .TRUE.)
    mmg_solfile = ListGetString(FuncParams, "MMG Output Sol Name", UnfoundFatal = .TRUE.)
  END IF

  NNodes = InMesh % NumberOfNodes
  NBulk = InMesh % NumberOfBulkElements
  NBdry = InMesh % NumberOfBoundaryElements

  IF( PRESENT( Hvar ) ) THEN
    IF( HVar % Dofs == 1 ) THEN
      SolType = MMG5_Scalar
      AnisoFlag = .FALSE.
    ELSE
      CALL Fatal(FuncName,'Implemented so far only for 1 dofs!')
    END IF      
  ELSE    
    AnisoFlag = ListGetLogical(FuncParams, "MMG Anisotropic", Found, DefValue=.TRUE.)  
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
        CALL ListGetRealArray(FuncParams,"MMG Target Length", WorkReal, 1, NodeNum, UnfoundFatal=.TRUE.)        
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
        Metric(i,:) = ListGetReal(FuncParams,"MMG Target Length", 1, NodeNum, UnfoundFatal=.TRUE.)
      END DO      
    END IF
  END IF
    
  nBCs = CurrentModel % NumberOfBCs
  body_offset = nBCs + CurrentModel % NumberOfBodies + 1
! body_offset = 0
  
  IF( body_offset > 0 ) THEN
    i=InMesh % NumberOfBulkElements
    InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID + body_offset
  END IF

  DO mmgloops = 1, MaxRemeshIter 
    CALL Info(FuncName,'Applying levelset trial MMG3D: '//TRIM(I2S(mmgloops)),Level=5)

    Success = .TRUE.

    ! Enable external depende on "mmg loop"
    mmgVar % Values(1) = 1.0_dp * mmgloops 

    pmmgMesh = 0
    
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

    ! If this is retrial then get only selected parameters again that may depend on the variable "MMG Loop".
    CALL Set_PMMG_Parameters(FuncParams, mmgloops > 1 )
        
    FreezeInternal = ListGetLogical( Params, 'MMG freeze internal boundaries', Found )
    IF (Present(PairCount)) THEN
      CALL SET_ParMMG_MESH(InMesh,Parallel,EdgePairs,PairCount,FreezeInternalArg=FreezeInternal)
    ELSE
      CALL SET_ParMMG_MESH(InMesh,Parallel,FreezeInternalArg=FreezeInternal)
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
        IF( PRESENT( HVar ) ) THEN
          CALL PMMG_Set_ScalarMet(pmmgMesh,HVar % Values(i),i,ierr)
        ELSE
          CALL PMMG_Set_ScalarMet(pmmgMesh,Metric(i,1),i,ierr)
        END IF
        IF(ierr /= 1) CALL Fatal(FuncName, "Failed to set scalar solution at vertex")
      END IF
    END DO

    ! compute globaldofs
    CALL PMMG_SET_IPARAMETER(pmmgMesh, PMMGPARAM_globalnum,1, ierr)

#if 0 
    IF(PRESENT(Angle)) THEN
      !Turn on sharp angle detection (1)
      CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle, &
          1,ierr)
      !Option to set angle detection threshold:
      CALL PMMG_SET_DPARAMETER(pmmgMesh,PMMGPARAM_angleDetection,&
          Angle,ierr)
    ELSE
      !Turn off sharp angle detection (0)
      CALL PMMG_SET_IPARAMETER(pmmgMesh,PMMGPARAM_angle, &
          0,ierr)
    END IF
#endif
    
    !Take care of fixed nodes/elements if requested
    IF(PRESENT(NodeFixed)) THEN
      DO i=1,NNodes
        IF(NodeFixed(InMesh % ParallelInfo % GlobalDOFs(i))) THEN
          CALL PMMG_SET_REQUIREDVERTEX(pmmgMesh,i,ierr)
        END IF
      END DO
    END IF

    IF(FreezeInternal) THEN
      DO i=nBulk+1,nBulk+nBdry
        IF ( inMesh % Elements(i) % Type % ElementCode /= 303 ) CYCLE
!       IF ( inMesh % Elements(i) % BoundaryInfo % Constraint<=0 ) CYCLE
 
        IF ( ASSOCIATED(inMesh % Elements(i) % BoundaryInfo % Right) ) THEN
          CALL PMMG_SET_REQUIREDTRIANGLE(pmmgMesh,i-nBulk,ierr)
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
    IF( SaveMMGMeshes .OR. SaveMMGSols ) THEN
      CALL Info(FuncName,'Saving of PMMG files finished', Level=20)
    END IF
      
    CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
    CALL PMMG_parmmglib_distributed(pmmgMesh,ierr)

    ! need to check that if one process returns failure all do
    ! based of serial remeshing failure routine so may need refinement
    IF ( ierr == PMMG_STRONGFAILURE .OR. ierr == PMMG_LOWFAILURE ) THEN
      CALL Warn(FuncName,'BAD ENDING OF PMMGLIB: UNABLE TO SAVE MESH')
      Success=.FALSE.
    ENDIF
    
    IF( Success ) EXIT
    
    IF( mmgloops > 1 ) THEN
      !! Redoing adaptive mesh, release the previous mmg mesh
      CALL MMG3D_Free_all(MMG5_ARG_start, &
          MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
          MMG5_ARG_end)      
    END IF
    
    IF( mmgloops == MaxRemeshIter ) GOTO 20
    
    CALL Info(FuncName,'PMMG_parmmglib_centralized done',Level=20)
  END DO

  IF(SaveMMGMeshes) THEN
    MeshName = TRIM(mmg_meshfile)//I2S(time)//'.mesh'
    CALL PMMG_SaveMesh_Distributed(pmmgMesh,MeshName,LEN(TRIM(MeshName)),ierr)
  END IF
  IF(SaveMMGSols) THEN
    SolName = TRIM(mmg_solfile)//I2S(time)//'.sol'
    CALL PMMG_SaveMet_Distributed(pmmgMesh,SolName,LEN(TRIM(SolName)),ierr)
  END IF

  CALL PMMG_Get_meshSize(pmmgMesh,NVerts,NTetras,NPrisms,NTris,NQuads,NEdges,ierr)

  !! GET THE NEW MESH
  CALL GET_ParMMG_MESH(OutMesh,Parallel,Calving=.FALSE.)

  !! Release mmg mesh
  CALL PMMG_Free_all ( PMMG_ARG_start,     &
       PMMG_ARG_ppParMesh,pmmgMesh,         &
       PMMG_ARG_end)

  NBulk = OutMesh % NumberOfBulkElements
  NBdry = OutMesh % NumberOfBoundaryElements
  
  !Reset the BodyIDs (see above)
  IF( body_offset > 0 ) THEN
    OutMesh % Elements(1:Nbulk) % BodyID = OutMesh % Elements(1:NBulk) % BodyID - body_offset
    i=InMesh % NumberOfBulkElements
    InMesh % Elements(1:i) % BodyID = InMesh % Elements(1:i) % BodyID - body_offset
  END IF

  !And delete the unneeded BC elems
  IF( ListGetLogical( FuncParams,'Adaptive Cat at Levelset', Found ) ) THEN
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
  END IF
    

  ! if remeshing has failed need to reset body ids
20 IF(.NOT. Success) THEN
    DO i=1,InMesh % NumberOfBulkElements
      InMesh % Elements(i) % BodyID = InMesh % Elements(i) % BodyID - body_offset
    END DO
  END IF

#else
  CALL Fatal(FuncName, "Remeshing utility PMMG has not been installed")
#endif

END SUBROUTINE DistributedRemeshParMMG


!------------------------------------------------------------------------------
! 2D remeshing routine copied from MMG2DSolver.F90 to library of Elmer.
! Here is the original author information:
! 
! ******************************************************************************
! *
! *  Author: F. Gillet-Chaulet (IGE)
! *  Email:  fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:    http://elmerice.elmerfem.org
! *
! *  Original Date: 13-07-2017, 
! *****************************************************************************
!------------------------------------------------------------------------------
  FUNCTION GET_MMG2D_MESH(MeshNumber,OutputFilename) RESULT(NewMesh)
    IMPLICIT NONE
!------------------------------------------------------------------------------
    INTEGER :: MeshNumber
    TYPE(Mesh_t), POINTER :: NewMesh
    CHARACTER(LEN=MAX_NAME_LEN) :: OutPutFileName
!------------------------------------------------------------------------------
#if HAVE_MMG
    TYPE(Element_t),POINTER ::  Element, Element0
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: np,nt,na,nq,ier,na0,np0,np1,nt0
    INTEGER :: ref,corner,required,ridge
    INTEGER :: parent,ied
    INTEGER :: tt, ii, jj, kk, ll, i, j, t, maxnodes
    LOGICAL :: Debug = .FALSE.
    LOGICAL :: Combine, Found
    INTEGER, POINTER :: Perm(:)
    INTEGER, ALLOCATABLE :: Perm0(:),ElemInd0(:)
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Mesh_t), POINTER :: OldMesh
    CHARACTER(:), ALLOCATABLE :: EquationName                              
    CHARACTER(*), PARAMETER :: FuncName="Get_MMG2D_Mesh"

    CALL Info(FuncName,'Copy Mesh to MMG format',Level=20)
    
    Solver => CurrentModel % Solver      
    OldMesh => Solver % Mesh
    
    !> a) get the size of the mesh: vertices,  triangles, edges
    
#if MMG_VERSION_LT(5,5) 
    CALL MMG2D_Get_meshSize(mmgMesh,np,nt,na,ier)
#else
    CALL MMG2D_Get_meshSize(mmgMesh,np,nt,nq,na,ier)
    IF (nq /= 0) CALL Fatal(FuncName,'Sorry no support for 404 elements')
#endif

    IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMGS_Get_meshSize failed!')
    CALL Info(FuncName,'MMG2D_Get_meshSize done',Level=30)
      

    nt0 = 0; np0 = 0; na0 = 0
    Combine = ListGetLogical( Solver % Values,'Keep unmeshed regions',Found )

    IF( Combine ) THEN
      CALL Info(FuncName,'Merging with unmeshed regions',Found )
      ! This still related to the old mesh
      Perm => Solver % Variable % Perm
      IF(.NOT. ASSOCIATED(Perm)) THEN
        CALL Fatal(FuncName,'Perm vector not assoicated!?')
      END IF
      
      EquationName = ListGetString( Solver % Values, 'Equation', Found)                
      nt0 = 0
      DO t=1,OldMesh % NumberOfBulkElements
        Element => OldMesh % Elements(t)
        IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
        nt0 = nt0 +1 
      END DO
      
      IF(nt0 == 0) THEN
        CALL Info('FuncName','No passive regions, skipping combining of meshes!')
        Combine = .FALSE.
      ELSE
        CALL Info('FuncName','Elements in passive regions: '//I2S(nt0))
      END IF
    END IF
      
    maxnodes = 3
    IF( Combine ) THEN        
      ALLOCATE(Perm0(OldMesh % NumberOfNodes))
      Perm0 = 0

      na0 = 0
      DO t=1,OldMesh % NumberOfBoundaryElements
        Element => OldMesh % Elements( OldMesh % NumberOfBulkElements + t )
        IF(ALL(Perm(Element % NodeIndexes) > 0) ) CYCLE
        na0 = na0 + 1
        maxnodes = MAX(Element % TYPE % NumberOfNodes,maxnodes)
      END DO

      IF(na0 > 0) THEN
        ALLOCATE(ElemInd0(OldMesh % NumberOfBulkElements))
        ElemInd0 = 0
      END IF

      nt0 = 0
      DO t=1,OldMesh % NumberOfBulkElements
        Element => OldMesh % Elements(t)
        IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
        nt0 = nt0 + 1
        Perm0(Element % NodeIndexes) = 1
        IF(na0 > 0) ElemInd0(t) = nt0
      END DO

      np0 = 0
      np1 = 0
      DO i=1,OldMesh % NumberOfNodes
        IF( Perm0(i) > 0 ) THEN
          IF( Perm(i) == 0 ) THEN
            ! Internal saved node
            np0 = np0 + 1
            Perm0(i) = np + np0
          ELSE
            ! Shared node
            np1 = np1 + 1
            Perm0(i) = -np1
          END IF
        END IF
      END DO
    END IF
    
    ! Initialize the new mesh stucture
    NewMesh => AllocateMesh(nt + nt0,na + na0,np + np0,ParEnv%PEs > 1)
    IF (MeshNumber > 0 ) THEN
      WRITE(NewMesh % Name,'(A,A,I0)') TRIM(OutPutFileName),'_N',MeshNumber
    ELSE
      NewMesh % Name = TRIM(OutPutFileName)
    END IF

    NewMesh % MaxElementNodes = maxnodes
    NewMesh % MeshDim = OldMesh % MeshDim    
    NewMesh % NumberOfNodes = np + np0 
    NewMesh % NumberOfBulkElements = nt + nt0 
    NewMesh % NumberOfBoundaryElements = na + na0
    
    CALL AllocateVector( NewMesh % Nodes % x, np + np0 )
    CALL AllocateVector( NewMesh % Nodes % y, np + np0 ) 
    CALL AllocateVector( NewMesh % Nodes % z, np + np0 ) 
    CALL AllocateVector( NewMesh % Elements, nt+nt0 + na+na0 )

    IF( Combine ) THEN
      ALLOCATE( NewMesh % InvPerm( NewMesh % NumberOfBulkElements ) )
      NewMesh % InvPerm = 0    
    END IF
    
    ! Get new nodes in 2d
    NewMesh % Nodes % z = 0._dp
    DO ii=1,np
      CALL MMG2D_Get_vertex(mmgMesh,&
          NewMesh % Nodes % x(ii), NewMesh % Nodes % y(ii),&
          ref,corner,required,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,'Call to  MMG2D_Get_vertex failed!')
    END DO

    IF( Combine ) THEN
      DO i=1,OldMesh % NumberOfNodes
        ii = Perm0(i)
        IF( ii > np ) THEN
          NewMesh % Nodes % x(ii) = OldMesh % Nodes % x(i)
          NewMesh % Nodes % y(ii) = OldMesh % Nodes % y(i)
        END IF
      END DO
    END IF
    CALL Info(FuncName,'MMG2D_Get_vertex done',Level=30)   

    IF( Combine ) THEN
      BLOCK
        REAL(KIND=dp) :: x0,y0,s2,s2min,s2minmax,s
        INTEGER :: j,jmin
        s2minmax = 0.0_dp
        DO i=1,OldMesh % NumberOfNodes 
          ii = Perm0(i)
          IF( ii < 0 ) THEN
            x0 = OldMesh % Nodes % x(i)
            y0 = OldMesh % Nodes % y(i)
            s2min = HUGE(s2min)
            jmin = 0
            DO j=1,np
              s2 = (Newmesh % Nodes % x(j) - x0)**2 + (Newmesh % Nodes % y(j) - y0)**2
              IF( s2 < s2min ) THEN
                s2min = s2
                jmin = j
              END IF
            END DO
            s2minmax = MAX(s2minmax, s2min)
            Perm0(i) = jmin
          END IF
        END DO        
        !PRINT *,'MaxMin dx:',SQRT(s2minmax)
        IF( SQRT(s2minmax) > 1.0e-8 ) THEN
          CALL Fatal(FuncName,'Could not accurately find the boundary nodes in the new mesh!')
        END IF
      END BLOCK
    END IF
    
    ! Get new triangle elements
    DO tt=1,NewMesh % NumberOfBulkElements
      Element => NewMesh % Elements(tt)
      Element % TYPE => GetElementType( 303 )
      Element % NDOFs = Element % TYPE % NumberOfNodes
      Element % ElementIndex = tt
      Element % GElementIndex = tt
      Element % PartIndex = ParEnv % myPE
      CALL AllocateVector(Element % NodeIndexes, 3)
      NodeIndexes => Element % NodeIndexes
      CALL MMG2D_Get_triangle(mmgMesh, &
          NodeIndexes(1), NodeIndexes(2), NodeIndexes(3), &
          Element % BodyId, required, ier)
    END DO

    IF( Combine ) THEN
      nt0 = 0
      DO t=1,OldMesh % NumberOfBulkElements        
        Element0 => OldMesh % Elements(t)
        IF ( CheckElementEquation( CurrentModel, Element0, EquationName ) ) CYCLE
        nt0 = nt0 +1 
               
        Element => NewMesh % Elements(nt+nt0)
        Element % TYPE => Element0 % Type
        Element % NDOFs = Element0 % NDOFs
        Element % ElementIndex = nt+nt0
        Element % GElementIndex = nt+nt0
        Element % PartIndex = ParEnv % myPE
        CALL AllocateVector(Element % NodeIndexes, SIZE(Element0 % NodeIndexes))
        Element % NodeIndexes = Perm0(Element0 % NodeIndexes)
        
        ! This may be needed for efficient interpolation of remaining dofs
        NewMesh % InvPerm(nt+nt0) = t 
                
        IF( ANY(Element % NodeIndexes < 0 ) ) THEN        
          PRINT *,'Perm0 small:',Perm0(Element0 % NodeIndexes), np0
        END IF
        IF( ANY(Element % NodeIndexes > np + np0 ) ) THEN        
          PRINT *,'Perm0 big:',Perm0(Element0 % NodeIndexes), np0
        END IF          
        Element % BodyId = Element0 % BodyId
      END DO
    END IF
       
    CALL Info(FuncName,'MMG2D_Get_triangle done',Level=30)

    
    ! Get new boundary elements
    kk = NewMesh % NumberOfBulkElements
    DO ii=1,na
      kk = kk + 1
      Element => NewMesh % Elements(kk)

      Element % TYPE => GetElementType( 202 )
      Element % NDOFs = Element % TYPE % NumberOfNodes
      Element % ElementIndex = kk
      Element % PartIndex = ParEnv % myPE

      CALL AllocateVector(Element % NodeIndexes, 2)
      NodeIndexes => Element % NodeIndexes
      CALL MMG2D_Get_edge(mmgMesh,NodeIndexes(1),NodeIndexes(2), &
          ref, ridge,required,ier)
      IF ( ier /= 1 ) CALL Fatal(FuncName,'Call to  MMG2D_Get_Edge failed!')

      ALLOCATE(Element % BoundaryInfo)
      Element % BoundaryInfo % Constraint = ref

      CALL MMG2D_GET_TRIFROMEDGE(mmgMesh,ii,parent,ied,ier)
      IF( parent < 1 .OR. parent > nt ) THEN
        PRINT *,'Edge:',ii,na,nt,nt0,parent,ied
        CALL Fatal(FuncName,'Parent out of range')
      END IF
      IF ( ier /= 1 ) CALL Fatal(FuncName,'Call to  MMG2D_Get_TRIFROMEDGE failed!')      
      Element % BoundaryInfo % Left => NewMesh % Elements(parent)
    END DO


    kk=NewMesh % NumberOfBulkElements+na
    IF( na0 > 0 ) THEN
      na0 = 0
      DO t=1,OldMesh % NumberOfBoundaryElements
        Element0 => OldMesh % Elements( OldMesh % NumberOfBulkElements + t )
        IF(ALL(Perm(Element0 % NodeIndexes) > 0) ) CYCLE      
        
        kk = kk+1
        Element => NewMesh % Elements(kk)        
        Element % TYPE => Element0 % Type
        Element % NDOFs = Element0 % NDOFs
        Element % ElementIndex = kk
        Element % PartIndex = ParEnv % myPE
        
        CALL AllocateVector(Element % NodeIndexes, SIZE(Element0 % NodeIndexes))
        Element % NodeIndexes = Perm0(Element0 % NodeIndexes)

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Constraint = Element0 % BoundaryInfo % Constraint

        IF( ASSOCIATED(Element0 % BoundaryInfo % Left) ) THEN
          i = ElemInd0(Element0 % BoundaryInfo % Left % ElementIndex)          
          IF(i>0) Element % BoundaryInfo % Left => NewMesh % Elements(i+nt0)
        END IF
        IF( ASSOCIATED(Element0 % BoundaryInfo % Right) ) THEN
          i = ElemInd0(Element0 % BoundaryInfo % Right % ElementIndex)          
          IF(i>0) Element % BoundaryInfo % Right => NewMesh % Elements(i+nt0)
        END IF
      END DO
    END IF
             
!   CALL SetMeshMaxDOFs(NewMesh)

    CALL Finalize_MMG_Mesh(NewMesh)

BLOCK

     CHARACTER(:), ALLOCATABLE :: ElementDef, ElementDef0
     INTEGER ::  i,Def_DOFs(10,6)
     TYPE(Solver_t), POINTER :: Solver
     LOGICAL :: stat, GotMesh = .FALSE.


      Solver => CurrentModel % Solver


      Def_Dofs = -1; Def_Dofs(:,1)=1

      ! Define what kind of element we are working with in this solver
      !-----------------------------------------------------------------
      ElementDef = ListGetString( Solver % Values, 'Element', stat )

      IF ( .NOT. stat ) THEN
        IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', stat ) ) THEN
           Solver % Def_Dofs(:,:,4) = 0  ! The final value is set when calling LoadMesh2
           IF ( .NOT. GotMesh ) Def_Dofs(:,4) = MAX(Def_Dofs(:,4),0 )
           i=i+1
           Solver % DG = .TRUE.
!          CYCLE
        ELSE
           ElementDef = "n:1"
        END IF
      END IF


      ElementDef0 = ElementDef
      DO WHILE(.TRUE.)
        j = INDEX( ElementDef0, '-' )
        IF (j>0) THEN
          !
          ! Read the element definition up to the next flag which specifies the
          ! target element set
          !
          ElementDef = ElementDef0(1:j-1)
        ELSE
          ElementDef = ElementDef0
        END IF
        !  Calling GetDefs fills Def_Dofs arrays:
        CALL GetDefs( ElementDef, Solver % Def_Dofs, Def_Dofs(:,:), .NOT. GotMesh )
        IF(j>0) THEN
          ElementDef0 = ElementDef0(j+1:)
        ELSE
          EXIT
        END IF
      END DO

      CALL PrepareMesh( CurrentModel, NewMesh, ParEnv % PEs>1 , Def_Dofs )
END BLOCK

CONTAINS



!------------------------------------------------------------------------------
!> This subroutine is used to fill Def_Dofs array of the solver structure.
!> Note that this subroutine makes no attempt to figure out the index of
!> the body, so all bodies are assigned with the same element definition.
!> A similar array of reduced dimension is also filled so as to figure out
!> the maximal-complexity definition over all solvers which use the same
!> global mesh.
!------------------------------------------------------------------------------
    SUBROUTINE GetDefs(ElementDef, Solver_Def_Dofs, Def_Dofs, Def_Dofs_Update)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN) :: ElementDef     !< an element definition string
      INTEGER, INTENT(OUT) :: Solver_Def_Dofs(:,:,:) !< Def_Dofs of the solver structure
      INTEGER, INTENT(INOUT) :: Def_Dofs(:,:)        !< holds the maximal-complexity definition on global mesh
      LOGICAL, INTENT(IN) :: Def_Dofs_Update         !< is .TRUE. when the definition refers to the global mesh
!------------------------------------------------------------------------------
      INTEGER, POINTER :: ind(:)
      INTEGER, TARGET :: Family(10)
      INTEGER :: i,j,l,n
      LOGICAL :: stat

      Family = [1,2,3,4,5,6,7,8,9,10]

      ! The default assumption is that the given element definition is applied 
      ! to all basic element families (note that the element sets 9 and 10 are
      ! not included since the explicit choice of the target family is 
      ! a part of the element definition string when the target index is
      ! deduced to be 9 or 10).
      !
      ind => Family(1:8)
      !
      ! If the element family is specified, change the target family 
      !
      IF (SEQL(ElementDef, 'point') )     ind => Family(1:1)
      IF (SEQL(ElementDef, 'line') )      ind => Family(2:2)
      IF (SEQL(ElementDef, 'tri') )       ind => Family(3:3)
      IF (SEQL(ElementDef, 'quad') )      ind => Family(4:4)
      IF (SEQL(ElementDef, 'tetra') )     ind => Family(5:5)
      IF (SEQL(ElementDef, 'pyramid') )   ind => Family(6:6)
      IF (SEQL(ElementDef, 'prism') )     ind => Family(7:7)
      IF (SEQL(ElementDef, 'brick') )     ind => Family(8:8)
      IF (SEQL(ElementDef, 'tri_face') )  ind => Family(9:9)
      IF (SEQL(ElementDef, 'quad_face') ) ind => Family(10:10)

      n = INDEX(ElementDef,'-')
      IF (n<=0) n=LEN_TRIM(ElementDef)
          
      j = INDEX( ElementDef(1:n), 'n:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,1) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,1) = MAX(Def_Dofs(ind,1), l)
      END IF
          
      j = INDEX( ElementDef(1:n), 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,2) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,2) = MAX(Def_Dofs(ind,2), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,3) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,3) = MAX(Def_Dofs(ind,3), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l

        ! Zero value triggers discontinuous approximation within LoadMesh2,
        ! substitute the default negative initialization value to avoid troubles:
        IF (l == 0) l = -1

        Solver_Def_Dofs(ind,:,4) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Solver_Def_Dofs(ind,:,4) = 0
          IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef(1:n), 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,5) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,5) = MAX(Def_Dofs(ind,5), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2)=='%' ) THEN
          ! Seeing a p-element definition starting as p:% means that a 
          ! a special keyword construct is used so that the degree of
          ! approximation can be evaluated by calling a MATC function.
          ! This special case is handled elsewhere and we now postpone
          ! setting the right value.
          Solver_Def_Dofs(ind,:,6) = 0
        ELSE
          READ( ElementDef(j+2:), * ) l
          Solver_Def_Dofs(ind,:,6) = l
          IF ( Def_Dofs_Update ) Def_Dofs(ind,6) = MAX(Def_Dofs(ind,6), l )
         END IF
      END IF

!------------------------------------------------------------------------------
    END SUBROUTINE GetDefs
!------------------------------------------------------------------------------

    

    
#else
  CALL Fatal('Get_MMG2D_Mesh', "Remeshing utility MMG has not been installed")
#endif
  
  END FUNCTION GET_MMG2D_MESH


  SUBROUTINE SET_MMG2D_SOL(Mesh,MeshSize,Scalar)
    IMPLICIT NONE
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: MeshSize
    LOGICAL :: Scalar
    REAL(KIND=dp) :: M11,M22,M12
    INTEGER :: NVert
    INTEGER :: ier
    INTEGER :: ii,jj
    LOGICAL :: Debug = .FALSE.
    CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="Set_MMG2D_Sol"

#ifdef HAVE_MMG
    IF(.NOT. ASSOCIATED(Mesh) ) THEN
      CALL Fatal(FuncName,'Mesh not associated!')
    END IF
    IF(.NOT. ASSOCIATED(MeshSize) ) THEN
      CALL Fatal(FuncName,'MeshSize variable not associated!')
    END IF

    NVert = Mesh % NumberOfNodes
    IF (Scalar) THEN
      CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NVert,MMG5_Scalar,ier)
    ELSE
      CALL MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NVert,MMG5_Tensor,ier)
    END IF
    IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMG2D_Set_SolSize failed!')
    CALL Info(FuncName,'MMG2D_Set_solSize done',Level=30)

    DO ii=1,NVert
      jj = ii
      IF( ASSOCIATED( MeshSize % Perm ) ) jj = MeshSize % Perm(ii)

      IF (Scalar) THEN
        CALL MMG2D_Set_scalarSol(mmgSol,&
            MeshSize%Values(jj),ii,ier)
      ELSE
        M11 = MeshSize % Values(3*jj-2)
        M22 = MeshSize % Values(3*jj-1)
        M12 = MeshSize % Values(3*jj-0)
        CALL MMG2D_Set_tensorSol(mmgSol,M11,M12,M22,ii,ier)
      ENDIF
      IF ( ier == 0 ) CALL Fatal(Funcname,'Call to MMG2D_Set_scalarSo failed!')
    END DO    
    CALL Info(FuncName,'MMG2D_Set_tensorSol done',Level=30)

#else
  CALL Fatal(FuncName,'Remeshing utility MMG has not been installed!')
#endif
        
  END SUBROUTINE Set_MMG2D_Sol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Set_MMG2D_Mesh(Mesh,Solver)
    IMPLICIT NONE
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER, OPTIONAL :: Solver
    
    CHARACTER(*), PARAMETER :: FuncName="Set_MMG2D_Mesh"
#ifdef HAVE_MMG
    TYPE(Element_t),POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: NVert,NEle,NEdge,Nquad
    INTEGER :: n
    INTEGER :: ier
    INTEGER :: ii,jj,tt
    INTEGER :: Ind
    LOGICAL :: Found, Debug = .FALSE.
    INTEGER, POINTER :: Perm(:)
    CHARACTER(:), ALLOCATABLE :: EquationName                              

    
    CALL Info(FuncName,'Setting 2D mesh using MMG',Level=20)
    
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal(FuncName,'Mesh not associated!')
    END IF

    IF( PRESENT( Solver ) ) THEN
      CALL Info(FuncName,'Setting only the active part of mesh')
      Perm => Solver % Variable % Perm
      IF(.NOT. ASSOCIATED(Perm)) THEN
        CALL Fatal(FuncName,'Perm vector not assoicated!?')
      END IF
      ii = Mesh % NumberOfNodes
      NVert = COUNT( Perm(1:ii) > 0 )

      EquationName = ListGetString( Solver % Values, 'Equation', Found)                
      NEle = 0
      DO tt=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(tt)
        IF ( .NOT. CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
        NEle = NEle +1 
      END DO
      
      NEdge = 0
      DO tt=1,Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements( Mesh % NumberOfBulkElements + tt )
        NodeIndexes => Element % NodeIndexes
        IF(ALL(Perm(NodeIndexes) > 0 ) ) NEdge = NEdge + 1 
      END DO
      IF( ParEnv % PEs > 1 ) THEN
        CALL Fatal(FuncName,'Cannot pick part of parallel mesh yet!')
      END IF
    ELSE
      NULLIFY( Perm )
      NVert = Mesh % NumberOfNodes
      NEle = Mesh % NumberOfBulkElements
      NEdge = Mesh % NumberOfBoundaryElements
    END IF
    ! support only 303 elements no 404
    Nquad=0

    !PRINT *,'Nele:',NVert, Nele, Nedge, Mesh % NumberOfBulkElements - Nele
    
    CALL Info(FuncName,'Setting mesh for MMG2D',Level=20)
    
#if MMG_VERSION_LT(5,5) 
    CALL MMG2D_Set_meshSize(mmgMesh,NVert,NEle,NEdge,ier)
#else
    CALL MMG2D_Set_meshSize(mmgMesh,NVert,NEle,Nquad,NEdge,ier)
#endif
    IF ( ier == 0 ) CALL Fatal(FuncName,&
        'Call to MMG2D_Set_meshSize failed!')
    CALL Info(FuncName,'MMG2D_Set_meshSize done',Level=30)

    DO ii=1,Mesh % NumberOfNodes
      jj = ii
      IF( ASSOCIATED( Perm ) ) THEN
        jj = Perm(ii)
        IF(jj==0) CYCLE
      END IF
      CALL MMG2D_Set_vertex(mmgMesh, &
          Mesh%Nodes%x(ii), Mesh%Nodes%y(ii), &
          0, jj , ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMG2D_Set_vertex failed!')
    END DO
    CALL Info(FuncName,'MMG2D_Set_vertex done',Level=30)

    jj = 0
    DO tt=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(tt) 

      IF( PRESENT( Solver ) ) THEN
        IF ( .NOT. CheckElementEquation( CurrentModel, Element, EquationName ) ) CYCLE
      END IF
      jj = jj+1
      
      IF( ParEnv % PEs > 1 ) THEN
        Ind = Element % GElementIndex
      ELSE
        ind = jj
      END IF
      
      NodeIndexes => Element % NodeIndexes
      
      IF (Element % TYPE % ElementCode /= 303) &
          CALL Fatal(FuncName,'Work only with 303 elements')
      n = Element % TYPE % NumberOfNodes

      IF( ASSOCIATED( Perm ) ) THEN
        CALL MMG2D_Set_triangle(mmgMesh, &
            Perm(NodeIndexes(1)), Perm(NodeIndexes(2)), Perm(NodeIndexes(3)), &
            Element % BodyId, jj, ier)
      ELSE
        CALL MMG2D_Set_triangle(mmgMesh, &
            NodeIndexes(1), NodeIndexes(2), NodeIndexes(3), &
            Element % BodyId, jj, ier)
      END IF        
      IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMG2D_Set_triangle failed!')
    END DO
    CALL Info(FuncName,'MMG2D_Set_triangle done',Level=30)

    jj = 0
    DO tt=1,Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements( Mesh % NumberOfBulkElements + tt )
      NodeIndexes => Element % NodeIndexes
      IF( ASSOCIATED(Perm) ) THEN
        IF(ANY(Perm(NodeIndexes) == 0) ) CYCLE
      END IF
      jj = jj + 1
      n = Element % TYPE % NumberOfNodes

      SELECT CASE( Element % TYPE % ElementCode )
      CASE(101)
        CYCLE
      CASE(202)
        CONTINUE
      CASE DEFAULT
        CALL Fatal(FuncName,'Work only with 202 boundary elements')
      END SELECT

      IF( ASSOCIATED(Perm) ) THEN
        CALL MMG2D_Set_Edge(mmgMesh, &
            Perm(NodeIndexes(1)), Perm(NodeIndexes(2)), &
            Element % BoundaryInfo % Constraint, jj, ier)
      ELSE
        CALL MMG2D_Set_Edge(mmgMesh, &
            NodeIndexes(1), NodeIndexes(2), &
            Element % BoundaryInfo % Constraint, jj, ier)
      END IF
      IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMG2D_Set_Edge failed!')
    END DO
    CALL Info(FuncName,'MMG2D_Set_Edge done',Level=30)
#else
    CALL Fatal(FuncName,'Remeshing utility MMG has not been installed')
#endif
    
  END SUBROUTINE SET_MMG2D_MESH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SET_MMG2D_Parameters(SolverParams)
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(*), PARAMETER :: FuncName="Set_MMG2D_Parameters"
#ifdef HAVE_MMG   
    REAL(KIND=dp) :: hsiz,Pval
    INTEGER :: ier
    LOGICAL :: NoAngleDetect
    INTEGER :: verbosity,MeMIncrease,Bucket,GMSHoption     
    LOGICAL :: DebugMode,NoInsert,NoSwap,NoMove,NoSurf
    LOGICAL :: Found, Stat
    INTEGER :: istat
    REAL(KIND=dp) :: hmin, hmax


    CALL Info(FuncName,'Setting MMG2D Parameters',Level=20)

    CALL Check_Parameters_Obsolite(SolverParams)
    
    ! Minimal mesh size:  hmin
    Hmin = ListGetConstReal( SolverParams, 'mmg hmin', Found)
    IF (Found) THEN
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hmin,Hmin,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <hmin> Failed')
    END IF

    ! Maximal mesh size - hmax
    Hmax = ListGetConstReal( SolverParams, 'mmg hmax', Found)
    IF (Found) THEN
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,Hmax,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <hmax> Failed')
    END IF

    hsiz = ListGetConstReal( SolverParams, 'mmg hsiz', Found)
    IF (Found) THEN
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hsiz,hsiz,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <hsiz> Failed')
    END IF
    
    ! Set the verbosity level to n
    Verbosity = ListGetInteger( SolverParams,'mmg verbosity',Found)
    IF (Found) THEN
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,Verbosity,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,&
          'Call to MMG2D_SET_IPARAMETER Failed')
    END IF
    
    ! Set the maximal memory size to n MBytes.
    MemIncrease = ListGetInteger(SolverParams,'mmg Increase Memory',Found)
    IF (FOUND) THEN
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_mem,MemIncrease,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_IPARAMETER Failed')
    END IF

    ! Turn on the debug mode (0/1).
    Stat = ListGetLogical(SolverParams,'mmg Debug Mode',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat=1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_debug,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,&
          'Call to MMG2D_SET_IPARAMETER Failed')
    END IF

    ! Control global Hausdorff distance (on all the boundary surfaces of the mesh)
    ! MMG2D_DPARAM_hausd default est 0.01 semble bien trop petit;
    !  il semble qu'il faille mettre une taille > taille des elements.
    Pval = ListGetConstReal( SolverParams, 'mmg hausd', Found)
    IF (Found) THEN
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hausd,Pval,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <hausd> Failed')
    END IF

    ! Control gradation
    Pval = ListGetConstReal( SolverParams, 'mmg hgrad', Found)
    IF (Found) THEN
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_hgrad,Pval,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <hgrad> Failed')
    END IF

    ! OTHER PARAMETERS: NOT ALL TESTED
    Pval = ListGetConstReal( SolverParams, 'mmg Angle detection',Found)
    IF (Found) THEN      
      CALL MMG2D_SET_DPARAMETER(mmgMesh,mmgSol,MMG2D_DPARAM_angleDetection, Pval,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_DPARAMETER <Angle detection> Failed')
    ENDIF

    ! [1/0], Avoid/allow surface modifications  
    Stat = ListGetLogical(SolverParams,'mmg No Angle detection',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat=1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_angle,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_IPARAMETER <No Angle detection> Failed')
    END IF

    ! [1/0] Avoid/allow point insertion
    Stat = ListGetLogical(SolverParams,'mmg No insert',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat=1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_noinsert,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName, &
          'Call to MMG2D_SET_IPARAMETER <No insert> Failed') 
    END IF

    ! [1/0] Avoid/allow edge or face flipping
    Stat = ListGetLogical(SolverParams,'mmg No swap',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat=1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_noswap,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,&
          'Call to MMG2D_SET_IPARAMETER <No swap>Failed')
    END IF

    ! [1/0] Avoid/allow point relocation
    Stat = ListGetLogical(SolverParams,'mmg no move',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat = 1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_nomove,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,&
          'Call to MMG2D_SET_IPARAMETER <No move> Failed')
    END IF
    
    ! [1/0] Avoid/allow surface modifications
    Stat = ListGetLogical(SolverParams,'mmg No surf',Found)
    IF (Found) THEN
      istat=0; IF(Stat) istat=1
      CALL MMG2D_SET_IPARAMETER(mmgMesh,mmgSol,MMG2D_IPARAM_nosurf,istat,ier)
      IF ( ier == 0 ) CALL Fatal(FuncName,&
          'Call to MMG2D_SET_IPARAMETER <No surf> Failed')
    END IF
#else
    CALL Fatal(FuncName, "Remeshing utility MMG has not been installed")
#endif
        
  END SUBROUTINE Set_MMG2D_Parameters



  FUNCTION MMG2D_ReMesh( RefMesh, Hvar, Solver) RESULT ( NewMesh ) 

    TYPE(Mesh_t), POINTER :: NewMesh, RefMesh
    TYPE(Variable_t), POINTER, OPTIONAL :: Hvar
    TYPE(Solver_t), POINTER, OPTIONAL :: Solver 
    CHARACTER(*), PARAMETER :: FuncName="MMG2D_ReMesh"
#ifdef HAVE_MMG    
    TYPE(ValueList_t), POINTER :: SolverParams
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: FileName
    LOGICAL :: Found, Numbering
    INTEGER :: ier, MeshNumber = 0

    CALL Info(FuncName,'Performing Remeshing using MMG library!')
        
    SolverParams => CurrentModel % Solver % Values
    Mesh => CurrentModel % Solver % Mesh
        
    CALL Info(FuncName,'Initialization of MMG',Level=20)
    mmgMesh = 0
    mmgSol  = 0
    CALL MMG2D_Init_mesh(MMG5_ARG_start, &
        MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
        MMG5_ARG_end)
    
    CALL SET_MMG2D_Parameters(SolverParams)

    IF( PRESENT( Solver ) ) THEN
      CALL SET_MMG2D_Mesh(RefMesh, Solver=Solver)
    ELSE
      CALL SET_MMG2D_Mesh(RefMesh)
    END IF
      
    IF(PRESENT(HVar)) THEN
      IF(ASSOCIATED(HVar)) THEN
        CALL Info(FuncName,'Set the local size field!',Level=20)      
        CALL SET_MMG2D_Sol(Mesh,HVar,Hvar % Dofs == 1 )
      END IF
    END IF
      
    CALL Info(FuncName,'Check the mesh data',Level=20)              
    CALL MMG2D_Chk_meshData(mmgMesh,mmgSol,ier)
    IF ( ier == 0 ) CALL Fatal(FuncName,'Call to MMG2D_Chk_meshData failed!')

    IF( ListGetLogical(SolverParams,'Save Initial MMG Mesh', Found ) ) THEN
      filename = "MMGini.mesh"      
      CALL MMG2D_SaveMesh(mmgMesh,TRIM(filename),len_TRIM(filename),ier)
      filename = "MMGini.sol"      
      CALL MMG2D_SaveSol(mmgMesh,mmgSol,TRIM(filename),len_TRIM(filename),ier)
    END IF

    CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
    IF ( ier == MMG5_STRONGFAILURE ) THEN
      CALL Fatal(FuncName,'Bad ending of MMG2DLIB: unable to save mesh!')
    ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
      CALL Warn(FuncName,'Bad ending of MMG2DLIB: trying to continue!')
    ENDIF
    
    IF( ListGetLogical(SolverParams,'Save Final MMG Mesh', Found ) ) THEN
      filename = "MMGfinal.mesh"      
      CALL MMG2D_SaveMesh(mmgMesh,TRIM(filename),len_TRIM(filename),ier)
    END IF
    
    Numbering = ListGetLogical(SolverParams,'Increment Mesh Number',Found)
    IF(.NOT. Found) Numbering = ListGetLogical(SolverParams,'Adaptive Mesh Numbering',Found)
    IF(.NOT. Found) Numbering = .TRUE.
    IF( Numbering ) MeshNumber = MeshNumber + 1
 
    filename = ListGetString( SolverParams, 'Adaptive Mesh Name', Found )
    IF(.NOT. Found) filename = 'RefinedMesh'

    NewMesh => Get_MMG2D_Mesh(MeshNumber,Filename)

    CALL Info(FuncName,'Remeshing by MMG2D library finished!',Level=12)
#else
    CALL Fatal(FuncName,'Remeshing utility MMG has not been installed!')
#endif
    
  END FUNCTION MMG2D_ReMesh



  SUBROUTINE Finalize_MMG_Mesh(Mesh)

     TYPE(Mesh_t), POINTER :: Mesh

     CHARACTER(:), ALLOCATABLE :: ElementDef, ElementDef0
     INTEGER ::  i,j,Def_DOFs(10,6)
     TYPE(Solver_t), POINTER :: Solver
     LOGICAL :: stat, GotMesh = .FALSE.


      Solver => CurrentModel % Solver


      Def_Dofs = -1; Def_Dofs(:,1)=1

      ! Define what kind of element we are working with in this solver
      !-----------------------------------------------------------------
      ElementDef = ListGetString( Solver % Values, 'Element', stat )

      IF ( .NOT. stat ) THEN
        IF ( ListGetLogical( Solver % Values, 'Discontinuous Galerkin', stat ) ) THEN
           Solver % Def_Dofs(:,:,4) = 0  ! The final value is set when calling LoadMesh2
           IF ( .NOT. GotMesh ) Def_Dofs(:,4) = MAX(Def_Dofs(:,4),0 )
           i=i+1
           Solver % DG = .TRUE.
!          CYCLE
        ELSE
           ElementDef = "n:1"
        END IF
      END IF

      ElementDef0 = ElementDef
      DO WHILE(.TRUE.)
        j = INDEX( ElementDef0, '-' )
        IF (j>0) THEN
          !
          ! Read the element definition up to the next flag which specifies the
          ! target element set
          !
          ElementDef = ElementDef0(1:j-1)
        ELSE
          ElementDef = ElementDef0
        END IF
        !  Calling GetDefs fills Def_Dofs arrays:
        CALL GetDefs( ElementDef, Solver % Def_Dofs, Def_Dofs(:,:), .NOT. GotMesh )
        IF(j>0) THEN
          ElementDef0 = ElementDef0(j+1:)
        ELSE
          EXIT
        END IF
      END DO

      CALL PrepareMesh( CurrentModel, Mesh, ParEnv % PEs>1 , Def_Dofs )

CONTAINS

!------------------------------------------------------------------------------
!> This subroutine is used to fill Def_Dofs array of the solver structure.
!> Note that this subroutine makes no attempt to figure out the index of
!> the body, so all bodies are assigned with the same element definition.
!> A similar array of reduced dimension is also filled so as to figure out
!> the maximal-complexity definition over all solvers which use the same
!> global mesh.
!------------------------------------------------------------------------------
    SUBROUTINE GetDefs(ElementDef, Solver_Def_Dofs, Def_Dofs, Def_Dofs_Update)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN) :: ElementDef     !< an element definition string
      INTEGER, INTENT(OUT) :: Solver_Def_Dofs(:,:,:) !< Def_Dofs of the solver structure
      INTEGER, INTENT(INOUT) :: Def_Dofs(:,:)        !< holds the maximal-complexity definition on global mesh
      LOGICAL, INTENT(IN) :: Def_Dofs_Update         !< is .TRUE. when the definition refers to the global mesh
!------------------------------------------------------------------------------
      INTEGER, POINTER :: ind(:)
      INTEGER, TARGET :: Family(10)
      INTEGER :: i,j,l,n
      LOGICAL :: stat

      Family = [1,2,3,4,5,6,7,8,9,10]

      ! The default assumption is that the given element definition is applied 
      ! to all basic element families (note that the element sets 9 and 10 are
      ! not included since the explicit choice of the target family is 
      ! a part of the element definition string when the target index is
      ! deduced to be 9 or 10).
      !
      ind => Family(1:8)
      !
      ! If the element family is specified, change the target family 
      !
      IF (SEQL(ElementDef, 'point') )     ind => Family(1:1)
      IF (SEQL(ElementDef, 'line') )      ind => Family(2:2)
      IF (SEQL(ElementDef, 'tri') )       ind => Family(3:3)
      IF (SEQL(ElementDef, 'quad') )      ind => Family(4:4)
      IF (SEQL(ElementDef, 'tetra') )     ind => Family(5:5)
      IF (SEQL(ElementDef, 'pyramid') )   ind => Family(6:6)
      IF (SEQL(ElementDef, 'prism') )     ind => Family(7:7)
      IF (SEQL(ElementDef, 'brick') )     ind => Family(8:8)
      IF (SEQL(ElementDef, 'tri_face') )  ind => Family(9:9)
      IF (SEQL(ElementDef, 'quad_face') ) ind => Family(10:10)

      n = INDEX(ElementDef,'-')
      IF (n<=0) n=LEN_TRIM(ElementDef)
          
      j = INDEX( ElementDef(1:n), 'n:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,1) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,1) = MAX(Def_Dofs(ind,1), l)
      END IF
          
      j = INDEX( ElementDef(1:n), 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,2) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,2) = MAX(Def_Dofs(ind,2), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,3) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,3) = MAX(Def_Dofs(ind,3), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l

        ! Zero value triggers discontinuous approximation within LoadMesh2,
        ! substitute the default negative initialization value to avoid troubles:
        IF (l == 0) l = -1

        Solver_Def_Dofs(ind,:,4) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Solver_Def_Dofs(ind,:,4) = 0
          IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef(1:n), 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Solver_Def_Dofs(ind,:,5) = l
        IF ( Def_Dofs_Update ) Def_Dofs(ind,5) = MAX(Def_Dofs(ind,5), l )
      END IF
          
      j = INDEX( ElementDef(1:n), 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2)=='%' ) THEN
          ! Seeing a p-element definition starting as p:% means that a 
          ! a special keyword construct is used so that the degree of
          ! approximation can be evaluated by calling a MATC function.
          ! This special case is handled elsewhere and we now postpone
          ! setting the right value.
          Solver_Def_Dofs(ind,:,6) = 0
        ELSE
          READ( ElementDef(j+2:), * ) l
          Solver_Def_Dofs(ind,:,6) = l
          IF ( Def_Dofs_Update ) Def_Dofs(ind,6) = MAX(Def_Dofs(ind,6), l )
         END IF
      END IF

!------------------------------------------------------------------------------
    END SUBROUTINE GetDefs
!------------------------------------------------------------------------------

END SUBROUTINE Finalize_MMG_Mesh
    
  
END MODULE MeshRemeshing
