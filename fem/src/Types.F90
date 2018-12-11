!/*****************************************************************************/
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
! *  Authors: Juha Ruokolainen, Jouni Malinen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> Type definitions for Elmer library.
!------------------------------------------------------------------------------

#include "../config.h"

MODULE Types
 
   USE Messages
   USE iso_varying_string
   USE, INTRINSIC :: ISO_C_BINDING
#ifdef _OPENMP
   USE omp_lib 
#endif 

#ifdef HAVE_LUA
   USE Lua
#endif

   INTEGER, PARAMETER :: MAX_NAME_LEN = 128, MAX_STRING_LEN=2048
   ! Parameter for internal blocking
   INTEGER, PARAMETER :: VECTOR_BLOCK_LENGTH = 128
   ! Parameter for internally avoiding calls to BLAS
   INTEGER, PARAMETER :: VECTOR_SMALL_THRESH = 9

#if defined(ARCH_32_BITS)
   INTEGER, PARAMETER :: AddrInt = SELECTED_INT_KIND(9)
#else
   INTEGER, PARAMETER :: AddrInt = SELECTED_INT_KIND(18)
#endif

   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)

   REAL(KIND=dp), PARAMETER :: AEPS = 10 * EPSILON(1.0_dp), &
         PI = 3.1415926535897932384626433832795_dp
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: MATRIX_CRS  = 1, &
                        MATRIX_BAND = 2, &
                        MATRIX_SBAND = 3, & 
                        MATRIX_LIST = 4
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: SOLVER_EXEC_NEVER      = -1, &
                        SOLVER_EXEC_ALWAYS     =  0, &
                        SOLVER_EXEC_AHEAD_ALL  =  1, &
                        SOLVER_EXEC_AHEAD_TIME =  2, &
                        SOLVER_EXEC_AFTER_ALL  =  3, &
                        SOLVER_EXEC_AFTER_TIME =  4, &
                        SOLVER_EXEC_AHEAD_SAVE =  5, &
                        SOLVER_EXEC_AFTER_SAVE =  6, &
                        SOLVER_EXEC_PREDCORR = 7

  INTEGER, PARAMETER :: SOLVER_MODE_DEFAULT = 0, &    ! normal pde
	                SOLVER_MODE_AUXILIARY = 1, &  ! no fem machinery (SaveData)
	                SOLVER_MODE_ASSEMBLY = 2, &   ! coupled solver with single block
	                SOLVER_MODE_COUPLED = 3, &    ! coupled solver with multiple blocks
	                SOLVER_MODE_BLOCK = 4, &      ! block solver
	                SOLVER_MODE_GLOBAL = 5, &     ! lumped variables (no mesh)
	                SOLVER_MODE_MATRIXFREE = 6, & ! normal field, no matrix
                        SOLVER_MODE_STEPS = 7         ! as the legacy but splitted to different steps

  INTEGER, PARAMETER :: PROJECTOR_TYPE_DEFAULT = 0, &  ! unspecified constraint matrix
                        PROJECTOR_TYPE_NODAL = 1, &    ! nodal projector
                        PROJECTOR_TYPE_GALERKIN = 2    ! Galerkin projector

  INTEGER, PARAMETER :: DIRECT_NORMAL = 0, & ! Normal direct method
                        DIRECT_PERMON = 1    ! Permon direct method
  
!------------------------------------------------------------------------------
  CHARACTER, PARAMETER :: Backslash = ACHAR(92)
!------------------------------------------------------------------------------



  
#ifndef USE_ISO_C_BINDINGS
INTERFACE
  SUBROUTINE Envir(a,b,len)
     USE, INTRINSIC :: ISO_C_BINDING
     INTEGER(C_INT) :: len
     CHARACTER(C_CHAR) :: a(*), b(*)
  END SUBROUTINE Envir

  SUBROUTINE SystemC(str)
     USE, INTRINSIC :: ISO_C_BINDING
     CHARACTER(C_CHAR) :: str(*)
  END SUBROUTINE SystemC

  SUBROUTINE MakeDirectory(str)
     USE, INTRINSIC :: ISO_C_BINDING
     CHARACTER(C_CHAR) :: str(*)
  END SUBROUTINE MakeDirectory

  SUBROUTINE Matc(cmd,VALUE,len)
     USE, INTRINSIC :: ISO_C_BINDING
     INTEGER(C_INT) :: len
     CHARACTER(C_CHAR) :: cmd(*), VALUE(*)
  END SUBROUTINE Matc
END INTERFACE
#endif

#ifdef HAVE_MUMPS
  INCLUDE 'dmumps_struc.h'
#endif


  TYPE BasicMatrix_t
    INTEGER :: NumberOfRows
    INTEGER, ALLOCATABLE :: Rows(:), Cols(:), Diag(:)
    INTEGER, ALLOCATABLE :: GRows(:), RowOwner(:)
    REAL(KIND=dp), ALLOCATABLE :: Values(:),MassValues(:), &
        DampValues(:),ILUValues(:),PrecValues(:)
  END TYPE BasicMatrix_t


  TYPE SubVector_t
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: rnorm, bnorm, xnorm
    REAL(KIND=dp), ALLOCATABLE :: rhs(:)
    REAL(KIND=dp), ALLOCATABLE :: DiagScaling(:)
    TYPE(Solver_t), POINTER :: Solver => NULL()
  END TYPE SubVector_t

  TYPE SubMatrix_t
    TYPE(Matrix_t), POINTER :: Mat
    TYPE(Matrix_t), POINTER :: PrecMat
  END TYPE SubMatrix_t

  TYPE BlockMatrix_t
    INTEGER :: NoVar = 0, MaxSize, TotSize
    INTEGER, POINTER :: Offset(:)
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: rnorm, bnorm, xnorm
    TYPE(SubMatrix_t), ALLOCATABLE :: SubMatrix(:,:)
    LOGICAL, ALLOCATABLE :: SubMatrixActive(:,:)
    TYPE(SubVector_t), POINTER :: SubVector(:) => NULL()
    INTEGER, POINTER :: BlockStruct(:)
    LOGICAL :: GotBlockStruct
    LOGICAL, ALLOCATABLE :: SubMatrixTranspose(:,:)
  END TYPE BlockMatrix_t

#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)                                 
  TYPE CPardiso_struct                                                          
    INTEGER :: n                                                                
    INTEGER :: mtype                                                            
    INTEGER :: msglvl                                                           
    INTEGER :: maxfct                                                           
    INTEGER :: mnum                                                             
    INTEGER :: nrhs                                                             
    INTEGER, POINTER CONTIG :: ia(:) => NULL(), ja(:) => NULL()
    REAL(kind=dp), POINTER CONTIG :: aa(:) => NULL(), rhs(:) => NULL(), &
          x(:) => NULL()
    INTEGER, POINTER CONTIG :: IParm(:) => NULL()    
    INTEGER(KIND=AddrInt), POINTER CONTIG :: ID(:) => NULL()                           
  END TYPE CPardiso_struct
#endif     

  TYPE Matrix_t
    TYPE(Matrix_t), POINTER :: Child => NULL(), Parent => NULL(), CircuitMatrix => Null(), &
        ConstraintMatrix=>NULL(), EMatrix=>NULL(), AddMatrix=>NULL(), CollectionMatrix=>NULL()

    INTEGER :: NumberOfRows, ExtraDOFs=0, ParallelDOFs=0

    TYPE(Solver_t), POINTER :: Solver => NULL()

    LOGICAL :: NoDirichlet = .FALSE.
    REAL(KIND=dp), ALLOCATABLE :: Dvalues(:)
    LOGICAL, ALLOCATABLE :: ConstrainedDOF(:)

    INTEGER :: Subband, FORMAT, SolveCount, Comm=-1
    LOGICAL :: Ordered, Lumped, Symmetric, COMPLEX, DGMatrix, Cholesky

    INTEGER :: ProjectorBC,ProjectorType

    TYPE(ListMatrix_t), POINTER :: ListMatrix(:) => NULL()

    INTEGER, POINTER :: Perm(:)=>NULL(),InvPerm(:)=>NULL(), Gorder(:)=>NULL(), EPerm(:)=>NULL()
    INTEGER, ALLOCATABLE :: GRows(:), RowOwner(:)
    INTEGER, POINTER CONTIG :: Rows(:)=>NULL(),Cols(:)=>NULL(), Diag(:)=>NULL()

    REAL(KIND=dp), POINTER CONTIG :: RHS(:)=>NULL(),BulkRHS(:)=>NULL(),RHS_im(:)=>NULL(),Force(:,:)=>NULL()
    REAL(KIND=dp), POINTER CONTIG :: BulkResidual(:)=>NULL()

    REAL(KIND=dp),  POINTER CONTIG :: Values(:)=>NULL(), ILUValues(:)=>NULL(), &
               DiagScaling(:) => NULL(), TValues(:) => Null()

    REAL(KIND=dp), ALLOCATABLE :: extraVals(:)
    REAL(KIND=dp) :: RhsScaling
    REAL(KIND=dp),  POINTER CONTIG :: MassValues(:)=>NULL(),DampValues(:)=>NULL(), &
        BulkValues(:)=>NULL(), BulkMassValues(:)=>NULL(), BulkDampValues(:)=>NULL(), &
        PrecValues(:)=>NULL()

#ifdef HAVE_FETI4I
    TYPE(C_PTR) :: PermonMatrix = C_NULL_PTR, PermonSolverInstance = C_NULL_PTR
#endif
#ifdef HAVE_MUMPS
    TYPE(dmumps_struc), POINTER :: MumpsID => NULL() ! Global distributed Mumps
    TYPE(dmumps_struc), POINTER :: MumpsIDL => NULL() ! Local domainwise Mumps
#endif
#if defined(HAVE_MKL) || defined(HAVE_PARDISO)
    INTEGER, POINTER :: PardisoParam(:) => NULL()
    INTEGER(KIND=AddrInt), POINTER :: PardisoID(:) => NULL()
#endif
#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)                                 
    TYPE(CPardiso_struct), POINTER :: CPardisoID => NULL()                      
#endif 
#ifdef HAVE_SUPERLU
    INTEGER(KIND=AddrInt) :: SuperLU_Factors=0
#endif
#ifdef HAVE_UMFPACK
    INTEGER(KIND=AddrInt) :: UMFPack_Numeric=0
#endif
#ifdef HAVE_CHOLMOD
    INTEGER(KIND=AddrInt) :: Cholmod=0
#endif
#ifdef HAVE_HYPRE
    INTEGER(KIND=C_INTPTR_T) :: Hypre=0
#endif
#ifdef HAVE_TRILINOS
    INTEGER(KIND=C_INTPTR_T) :: Trilinos=0
#endif
    INTEGER(KIND=AddrInt) :: SpMV=0

    INTEGER(KIND=AddrInt) :: MatVecSubr = 0

    INTEGER, POINTER CONTIG :: ILURows(:)=>NULL(),ILUCols(:)=>NULL(),ILUDiag(:)=>NULL()

!   For Complex systems, not used yet!:
!   -----------------------------------
    COMPLEX(KIND=dp), POINTER :: CRHS(:)=>NULL(),CForce(:,:)=>NULL()
    COMPLEX(KIND=dp),  POINTER :: CValues(:)=>NULL(),CILUValues(:)=>NULL()
    COMPLEX(KIND=dp),  POINTER :: CMassValues(:)=>NULL(),CDampValues(:)=>NULL()

! For Flux Corrected Transport 
    REAL(KIND=dp), POINTER :: FCT_D(:) => NULL()
    REAL(KIND=dp), POINTER :: MassValuesLumped(:) => NULL()

    TYPE(ParallelInfo_t), POINTER :: ParallelInfo=>NULL()
    TYPE(SParIterSolverGlobalD_t), POINTER :: ParMatrix=>NULL()
  END TYPE Matrix_t
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Typedefs for parallel solver 
!------------------------------------------------------------------------------

  TYPE ParEnv_t
     INTEGER                          :: PEs
     INTEGER                          :: MyPE
     LOGICAL                          :: Initialized
     INTEGER                          :: ActiveComm
     LOGICAL, DIMENSION(:), POINTER   :: Active
     LOGICAL, DIMENSION(:), POINTER   :: IsNeighbour
     LOGICAL, DIMENSION(:), POINTER   :: SendingNB
     INTEGER                          :: NumOfNeighbours
     INTEGER                          :: NumberOfThreads = 1
   END TYPE ParEnv_t


  TYPE GlueTableT
     INTEGER, DIMENSION(:), POINTER :: Rows=>NULL(), &
                Cols=>NULL(), Inds=>NULL(), RowOwner=>NULL()
  END TYPE GlueTableT


  TYPE VecIndicesT
     INTEGER, DIMENSION(:), POINTER :: RevInd=>NULL()
  END TYPE VecIndicesT


  TYPE IfVecT
     REAL(KIND=dp), DIMENSION(:), POINTER :: IfVec=>NULL()
  END TYPE IfVecT


  TYPE RHST
     REAL(KIND=dp), DIMENSION(:), POINTER :: RHSvec=>NULL()
     INTEGER, DIMENSION(:), POINTER :: RHSind=>NULL()
  END TYPE RHST


  TYPE DPBufferT
     REAL(KIND=dp), DIMENSION(:), POINTER :: DPBuf=>NULL()
  END TYPE DPBufferT


  TYPE ResBufferT
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ResVal
     INTEGER, DIMENSION(:), ALLOCATABLE :: ResInd
  END TYPE ResBufferT


  TYPE IfLColsT
     INTEGER, DIMENSION(:), POINTER :: IfVec=>NULL()
  END TYPE IfLColsT


  TYPE SplittedMatrixT
     TYPE (BasicMatrix_t), DIMENSION(:), POINTER :: IfMatrix=>NULL()
     TYPE (Matrix_t), POINTER :: InsideMatrix=>NULL()
     TYPE (BasicMatrix_t), DIMENSION(:), POINTER :: NbsIfMatrix=>NULL()

     TYPE (VecIndicesT), DIMENSION(:), POINTER :: VecIndices=>NULL()
     TYPE (IfVecT), DIMENSION(:), POINTER :: IfVecs=>NULL()
     TYPE (IfLColsT), DIMENSION(:), POINTER :: IfORows=>NULL()
     TYPE (IfLColsT), DIMENSION(:), POINTER :: IfLCols=>NULL()
     TYPE (GlueTableT), POINTER :: GlueTable=>NULL()
     TYPE (RHST), DIMENSION(:), POINTER :: RHS=>NULL()
     TYPE (ResBufferT), DIMENSION(:), POINTER :: ResBuf=>NULL()
     REAL(KIND=dp), POINTER CONTIG :: &
           Work(:,:)=>NULL(),TmpXVec(:)=>NULL(),TmpRVec(:)=>NULL()
  END TYPE SplittedMatrixT


  TYPE SParIterSolverGlobalD_t
     TYPE (SplittedMatrixT), POINTER :: SplittedMatrix=>NULL()
     TYPE (Matrix_t), POINTER :: Matrix=>NULL()
     TYPE (ParallelInfo_t), POINTER :: ParallelInfo=>NULL()
     TYPE(ParEnv_t) :: ParEnv
     INTEGER :: DOFs, RelaxIters
  END TYPE SParIterSolverGlobalD_t

  TYPE(SParIterSolverGlobalD_t), POINTER :: ParMatrix

!-------------------------------------------------------------------------------

   !
   ! Basis function type
   !
   TYPE BasisFunctions_t 
      INTEGER :: n
      INTEGER, POINTER :: p(:)=>NULL(),q(:)=>NULL(),r(:)=>NULL()
      REAL(KIND=dp), POINTER :: coeff(:)=>NULL()
   END TYPE BasisFunctions_t


   !
   ! Element type description 
   !
   TYPE ElementType_t
     TYPE(ElementType_t),POINTER :: NextElementType ! this is a list of types

     INTEGER :: ElementCode                         ! numeric code for element

     INTEGER :: BasisFunctionDegree, &              ! linear or quadratic
                NumberOfNodes, &                
                NumberOfEdges, &                
                NumberOfFaces, &                
                DIMENSION                           ! 1=line, 2=surface, 3=volume

     INTEGER :: GaussPoints,GaussPoints2, GaussPoints0 ! number of gauss points to use

     REAL(KIND=dp) :: StabilizationMK               ! stab.param. depending on
                                                    ! interpolation type

     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:)
     REAL(KIND=dp), DIMENSION(:), POINTER CONTIG :: NodeU, NodeV, NodeW
   END TYPE ElementType_t

!------------------------------------------------------------------------------
   TYPE ValueListEntry_t
     INTEGER :: TYPE
     TYPE(ValueListEntry_t), POINTER :: Next => Null()

     REAL(KIND=dp), POINTER :: TValues(:), Cumulative(:) => NULL()
     REAL(KIND=dp), POINTER :: FValues(:,:,:), CubicCoeff(:)=>NULL()
     INTEGER :: Fdim = 0 
     
     LOGICAL :: LValue
     INTEGER, POINTER :: IValues(:)

     INTEGER(KIND=AddrInt) :: PROCEDURE

     REAL(KIND=dp) :: Coeff = 1.0_dp
     CHARACTER(LEN=MAX_NAME_LEN) :: CValue

     INTEGER :: NameLen,DepNameLen = 0
     CHARACTER(LEN=MAX_NAME_LEN) :: Name,DependName

#ifdef DEVEL_LISTCOUNTER 
     INTEGER :: Counter = 0
#endif

#ifdef HAVE_LUA
     LOGICAL :: LuaFun = .FALSE.
     !CHARACTER(len=:), ALLOCATABLE :: LuaCmd
#endif
     
   END TYPE ValueListEntry_t

   TYPE ValueList_t
     TYPE(ValueListEntry_t), POINTER :: Head => Null()
   END TYPE ValueList_t

   
   TYPE VariableTable_t     
     TYPE(Variable_t), POINTER :: Variable
   END TYPE VariableTable_t

   
   ! This is a tentative data type to speed up the retrieval of parameters
   ! at elements.
   !----------------------------------------------------------------------
   TYPE ValueHandle_t
     INTEGER :: ValueType = -1
     INTEGER :: SectionType = -1
     INTEGER :: ListId = -1
     LOGICAL :: BulkElement
     TYPE(Element_t), POINTER :: Element => NULL()
     TYPE(ValueList_t), POINTER :: List => Null()
     TYPE(ValueList_t), POINTER :: Ptr  => Null()
     TYPE(Nodes_t), POINTER :: Nodes
     INTEGER, POINTER :: Indexes
     INTEGER :: n
     INTEGER :: nValuesVec = 0
     REAL(KIND=dp), POINTER :: ValuesVec(:) => NULL()
     REAL(KIND=dp), POINTER :: Values(:) => NULL()
     REAL(KIND=dp), POINTER :: ParValues(:,:) => NULL()
     INTEGER :: ParNo = 0
     INTEGER :: IValue, DefIValue = 0
     REAL(KIND=dp) :: RValue, DefRValue = 0.0_dp
     INTEGER :: Rdim = 0
     REAL(KIND=dp), POINTER :: RTensor(:,:) => NULL()
     REAL(KIND=dp), POINTER :: RTensorValues(:,:,:) => NULL()
     LOGICAL :: LValue, DefLValue = .FALSE.
     CHARACTER(LEN=MAX_NAME_LEN) :: CValue
     INTEGER :: CValueLen
     LOGICAL :: Found
     CHARACTER(LEN=MAX_NAME_LEN) :: Name
     LOGICAL :: Initialized = .FALSE.
     LOGICAL :: AllocationsDone = .FALSE.
     LOGICAL :: ConstantEverywhere = .FALSE.
     LOGICAL :: GlobalEverywhere = .FALSE.
     LOGICAL :: GlobalInList = .FALSE.
     LOGICAL :: EvaluateAtIP = .FALSE.
     LOGICAL :: SomeVarAtIp = .FALSE.
     LOGICAL :: SomewhereEvaluateAtIP = .FALSE.
     LOGICAL :: NotPresentAnywhere = .FALSE.
     LOGICAL :: UnfoundFatal = .FALSE.
     REAL(KIND=dp) :: minv, maxv
     LOGICAL :: GotMinv = .FALSE., GotMaxv = .FALSE.
     TYPE(VariableTable_t) :: VarTable(32)
     INTEGER :: VarCount

     TYPE(ValueHandle_t), POINTER :: HandleIm
   END TYPE ValueHandle_t


   TYPE VariableHandle_t     
     TYPE(Variable_t), POINTER :: Variable=>NULL()
     REAL(KIND=dp),POINTER :: Values(:)=>NULL()
     INTEGER,POINTER :: Perm(:)=>NULL()
     INTEGER :: tstep = 0
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: ActiveElement = .FALSE.
     REAL(KIND=dp) :: ElementValues(100)
     INTEGER :: n = 0
   END TYPE VariableHandle_t
   
   
!------------------------------------------------------------------------------

   TYPE MaterialArray_t
     TYPE(ValueList_t), POINTER :: Values => Null()
   END TYPE MaterialArray_t

!------------------------------------------------------------------------------

   TYPE BoundaryConditionArray_t
     INTEGER :: TYPE,Tag
     TYPE(Matrix_t), POINTER :: PMatrix => NULL()
     LOGICAL :: PMatrixGalerkin = .FALSE.
     TYPE(ValueList_t), POINTER :: Values => Null()
   END TYPE BoundaryConditionArray_t

!------------------------------------------------------------------------------

   TYPE InitialConditionArray_t
     INTEGER :: TYPE,Tag
     TYPE(ValueList_t), POINTER :: Values => Null()
   END TYPE InitialConditionArray_t

!------------------------------------------------------------------------------

    TYPE ComponentArray_t
      TYPE(ValueList_t), POINTER :: Values => Null()
    END TYPE ComponentArray_t

!------------------------------------------------------------------------------

    TYPE BodyForceArray_t
      TYPE(ValueList_t), POINTER :: Values => Null()
    END TYPE BodyForceArray_t

!------------------------------------------------------------------------------

    TYPE BoundaryArray_t
      TYPE(ValueList_t), POINTER :: Values => Null()
    END TYPE BoundaryArray_t

!------------------------------------------------------------------------------

    TYPE BodyArray_t
      TYPE(ValueList_t), POINTER :: Values => Null()
    END TYPE BodyArray_t

!------------------------------------------------------------------------------

    TYPE EquationArray_t
      TYPE(ValueList_t), POINTER :: Values => Null()
    END TYPE EquationArray_t

!------------------------------------------------------------------------------

!   TYPE SimulationInfo_t
!     TYPE(ValueList_t) :: Values => Null()
!   END TYPE SimulationInfo_t

!------------------------------------------------------------------------------
   INTEGER, PARAMETER :: Variable_on_nodes  = 0
   INTEGER, PARAMETER :: Variable_on_edges  = 1
   INTEGER, PARAMETER :: Variable_on_faces  = 2
   INTEGER, PARAMETER :: Variable_on_nodes_on_elements = 3
   INTEGER, PARAMETER :: Variable_on_gauss_points = 4
   INTEGER, PARAMETER :: Variable_on_elements = 5
   INTEGER, PARAMETER :: Variable_global = 6

   
   
   TYPE IntegrationPointsTable_t
     INTEGER :: IPCount = 0
     INTEGER, POINTER :: IPOffset(:)
     !TYPE(GaussIntegrationPoints_t), POINTER :: IPs
   END TYPE IntegrationPointsTable_t
      
   
!  TYPE Variable_Component_t
!     CHARACTER(LEN=MAX_NAME_LEN) :: Name
!     INTEGER :: DOFs, Type
!  END TYPE Variable_Component_t

   TYPE Variable_t
     TYPE(Variable_t), POINTER :: Next => NULL()
     TYPE(Variable_t), POINTER :: EVar => NULL() 
     INTEGER :: NameLen = 0
     CHARACTER(LEN=MAX_NAME_LEN) :: Name

     TYPE(Solver_t), POINTER :: Solver => NULL()
     LOGICAL :: Valid, Output
     TYPE(Mesh_t), POINTER :: PrimaryMesh => NULL()

     LOGICAL :: ValuesChanged = .FALSE.

! Some variables are created from pointers to the primary variables
     LOGICAL :: Secondary = .FALSE.

     INTEGER :: TYPE = Variable_on_nodes

     INTEGER :: DOFs = 0
     INTEGER, POINTER          :: Perm(:) => NULL()
     REAL(KIND=dp)             :: Norm=0, PrevNorm=0,NonlinChange=0, SteadyChange=0
     INTEGER :: NonlinConverged=-1, SteadyConverged=-1, NonlinIter=-1
     INTEGER :: LinConverged=-1
     COMPLEX(KIND=dp), POINTER :: EigenValues(:) => NULL(), &
          EigenVectors(:,:) => NULL()
     REAL(KIND=dp), POINTER :: ConstraintModes(:,:) => NULL()
     INTEGER, POINTER :: ConstraintModesIndeces(:) => NULL()
     INTEGER :: NumberOfConstraintModes = 0
     REAL(KIND=dp), POINTER :: Values(:) => NULL() ,&
          PrevValues(:,:) => NULL(), &
          PValues(:) => NULL(), NonlinValues(:) => NULL(), &
          SteadyValues(:) => NULL()
     LOGICAL, POINTER :: UpperLimitActive(:) => NULL(), LowerLimitActive(:) => NULL()
     COMPLEX(KIND=dp), POINTER :: CValues(:) => NULL()
     TYPE(IntegrationPointsTable_t), POINTER :: IPTable => NULL()
   END TYPE Variable_t

!------------------------------------------------------------------------------
   TYPE ListMatrixEntry_t
     INTEGER :: Index = -1
     REAL(KIND=dp) :: Value = 0.0
     TYPE(ListMatrixEntry_t), POINTER :: Next => NULL()
   END TYPE ListMatrixEntry_t

   TYPE ListMatrixEntryPool_t
      TYPE(ListMatrixEntry_t), ALLOCATABLE :: Entries(:)
      INTEGER :: NextIndex = 0
      TYPE(ListMatrixEntryPool_t), POINTER :: Next => NULL()
   END type ListMatrixEntryPool_t

   TYPE ListMatrixPool_t
     TYPE(ListMatrixEntryPool_t), POINTER :: EntryPool => NULL()
     TYPE(ListMatrixEntry_t), POINTER :: Deleted => NULL()
     INTEGER :: PoolSize = 0
   END TYPE ListMatrixPool_t
   
   TYPE ListMatrix_t
     INTEGER :: Degree, Level
     TYPE(ListMatrixEntry_t), POINTER :: Head
   END TYPE ListMatrix_t

   TYPE ListMatrixArray_t
     TYPE(ListMatrix_t), ALLOCATABLE :: Rows(:)
     TYPE(ListMatrixPool_t), ALLOCATABLE :: Pool(:)
#ifdef _OPENMP
     INTEGER(KIND=omp_lock_kind), ALLOCATABLE :: RowLocks(:)
#endif
   END TYPE ListMatrixArray_t
   
!------------------------------------------------------------------------------

   TYPE Factors_t 
     INTEGER :: NumberOfFactors = 0, NumberOfImplicitFactors = 0
     INTEGER, POINTER :: Elements(:) => NULL()
     REAL(KIND=dp), POINTER :: Factors(:) => NULL()
   END TYPE Factors_t

!-------------------------------------------------------------------------------

   TYPE BoundaryInfo_t
     TYPE(Factors_t), POINTER :: GebhardtFactors=>NULL()
     INTEGER :: Constraint = 0, OutBody = -1
     TYPE(Element_t), POINTER :: Left =>NULL(), Right=>NULL()
   END TYPE BoundaryInfo_t

!-------------------------------------------------------------------------------

   TYPE ElementData_t
     TYPE(ElementData_t), POINTER :: Next=>NULL()
     TYPE(varying_string) :: Name
     REAL(KIND=dp), POINTER :: Values(:)=>NULL()
   END TYPE ElementData_t

!-------------------------------------------------------------------------------

   TYPE Element_t
     TYPE(ElementType_t), POINTER :: TYPE => NULL()

     LOGICAL :: Copy = .FALSE.

     INTEGER :: BodyId=0, Splitted=0
     REAL(KIND=dp) :: StabilizationMK,hK

     TYPE(BoundaryInfo_t),  POINTER :: BoundaryInfo => NULL()

     INTEGER :: ElementIndex=-1, GElementIndex=-1, PartIndex=-1, NDOFs=0, BDOFs=0, DGDOFs=0
     INTEGER, DIMENSION(:), POINTER CONTIG :: &
         NodeIndexes => NULL(), EdgeIndexes   => NULL(), &
         FaceIndexes => NULL(), BubbleIndexes => NULL(), &
         DGIndexes   => NULL()
!DIR$ ATTRIBUTES ALIGN:64::NodeIndexes, EdgeIndexes, FaceIndexes, BubbleIndexes, DGIndexes

     TYPE(PElementDefs_t), POINTER :: PDefs=>NULL()
     TYPE(ElementData_t),  POINTER :: PropertyData=>NULL()
   END TYPE Element_t

!-------------------------------------------------------------------------------

   TYPE PElementDefs_t
      INTEGER :: P
      INTEGER :: TetraType       ! Type of p tetrahedron={0,1,2}
      LOGICAL :: isEdge          ! Is element an edge or face?
      INTEGER :: GaussPoints     ! Number of gauss points to use when using p elements
      LOGICAL :: pyramidQuadEdge ! Is element an edge of pyramid quad face?
      INTEGER :: localNumber     ! Local number of an edge or face for element on boundary
   END TYPE PElementDefs_t

!-------------------------------------------------------------------------------

   TYPE NeighbourList_t
     INTEGER, DIMENSION(:), POINTER :: Neighbours
   END TYPE NeighbourList_t

!------------------------------------------------------------------------------

   !
   ! Coordinate and vector type definition, coordinate arrays must be allocated
   ! prior to use of variables of this type.
   !
   TYPE Nodes_t
     INTEGER :: NumberOfNodes
     REAL(KIND=dp), ALLOCATABLE :: xyz(:,:)
     REAL(KIND=dp), POINTER CONTIG :: x(:)=>NULL()
     REAL(KIND=dp), POINTER CONTIG :: y(:)=>NULL()
     REAL(KIND=dp), POINTER CONTIG :: z(:)=>NULL()
!DIR$ ATTRIBUTES ALIGN:64::x,y,z,xyz
   END TYPE Nodes_t

!------------------------------------------------------------------------------

   TYPE QuadrantPointer_t
     TYPE(Quadrant_t), POINTER :: Quadrant
   END TYPE QuadrantPointer_t

!------------------------------------------------------------------------------

   TYPE Quadrant_t
     INTEGER, DIMENSION(:), POINTER :: Elements
     REAL(KIND=dp) :: SIZE, MinElementSize, BoundingBox(6)
     INTEGER :: NElemsInQuadrant
     TYPE(QuadrantPointer_t), DIMENSION(:), POINTER :: ChildQuadrants
   END TYPE Quadrant_t

!------------------------------------------------------------------------------

   TYPE Projector_t
     TYPE(Projector_t), POINTER :: Next
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Matrix_t), POINTER :: Matrix, TMatrix
   END TYPE Projector_t


!------------------------------------------------------------------------------

   TYPE ParallelInfo_t
     INTEGER :: NumberOfIfDOFs
     LOGICAL, POINTER               :: INTERFACE(:)
     INTEGER, POINTER               :: GlobalDOFs(:)
     TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
     INTEGER, POINTER               :: Gorder(:) => Null()

     LOGICAL, POINTER               :: FaceInterface(:)
     TYPE(NeighbourList_t),POINTER  :: FaceNeighbourList(:)

     LOGICAL, POINTER               :: EdgeInterface(:)
     TYPE(NeighbourList_t),POINTER  :: EdgeNeighbourList(:)
   END TYPE ParallelInfo_t

!------------------------------------------------------------------------------

   TYPE Mesh_t
     CHARACTER(MAX_NAME_LEN) :: Name
     TYPE(Mesh_t), POINTER   :: Next,Parent,Child

     TYPE(Projector_t), POINTER :: Projector
     TYPE(Quadrant_t), POINTER  :: RootQuadrant

     LOGICAL :: Changed, OutputActive, Stabilize
     INTEGER :: SavesDone, AdaptiveDepth

     TYPE(Factors_t), POINTER :: ViewFactors(:)

     TYPE(ParallelInfo_t) :: ParallelInfo
     TYPE(Variable_t), POINTER :: Variables

     TYPE(Nodes_t), POINTER :: Nodes
     TYPE(Element_t), DIMENSION(:), POINTER :: Elements, Edges, Faces
     TYPE(Nodes_t), POINTER :: NodesOrig
     TYPE(Nodes_t), POINTER :: NodesMapped

     LOGICAL :: DisContMesh 
     INTEGER, POINTER :: DisContPerm(:)
     INTEGER :: DisContNodes
     
     INTEGER, POINTER :: InvPerm(:)

     INTEGER :: NumberOfNodes, NumberOfBulkElements, NumberOfEdges, &
                NumberOfFaces, NumberOfBoundaryElements, MeshDim = 0, MaxDim = 0, PassBCcnt=0
     INTEGER :: MinEdgeDOFs, MinFaceDOFs
     INTEGER :: MaxElementNodes, MaxElementDOFs, MaxEdgeDOFs, MaxFaceDOFs, MaxBDOFs

     LOGICAL :: EntityWeightsComputed 
     REAL(KIND=dp), POINTER :: BCWeight(:), BodyForceWeight(:),&
         BodyWeight(:), MaterialWeight(:)

     INTEGER, POINTER :: RePartition(:) => NULL()
     
   END TYPE Mesh_t

   TYPE Graph_t
     INTEGER :: n
     INTEGER, ALLOCATABLE :: ptr(:), ind(:)
   END type Graph_t
   
   TYPE Graphcolour_t
     INTEGER :: nc
     INTEGER, ALLOCATABLE :: colours(:)
   END TYPE Graphcolour_t

   TYPE MortarBC_t 
     TYPE(Matrix_t), POINTER :: Projector => NULL()
     INTEGER, POINTER :: Perm(:) => NULL()
     REAL(KIND=dp), POINTER :: Rhs(:) => NULL()
     REAL(KIND=dp), POINTER :: Diag(:) => NULL()
     LOGICAL, POINTER :: Active(:) => NULL()
     REAL(KIND=dp) :: SlaveScale = 1.0_dp
     REAL(KIND=dp) :: MasterScale = 1.0_dp
     LOGICAL :: LumpedDiag = .TRUE.
   END TYPE MortarBC_t

!------------------------------------------------------------------------------

!   TYPE Constants_t
!     REAL(KIND=dp) :: Gravity(4)
!     REAL(KIND=dp) :: StefanBoltzmann
!   END TYPE Constants_t

!------------------------------------------------------------------------------

    TYPE Solver_t
      INTEGER :: SolverId = 0
      TYPE(ValueList_t), POINTER :: Values => NULL()

      INTEGER :: TimeOrder,DoneTime,Order,NOFEigenValues=0
      INTEGER :: TimesVisited = 0
      INTEGER(KIND=AddrInt) :: PROCEDURE, LinBeforeProc, LinAfterProc

      REAL(KIND=dp) :: Alpha,Beta,dt

      LOGICAL :: NewtonActive = .FALSE.
      
      INTEGER :: SolverExecWhen
      INTEGER :: SolverMode

      INTEGER :: MultiGridLevel,  MultiGridTotal, MultiGridSweep
      LOGICAL :: MultiGridSolver, MultiGridEqualSplit
      TYPE(Mesh_t), POINTER :: Mesh => NULL()

      INTEGER, POINTER :: ActiveElements(:) => NULL()
      INTEGER, POINTER :: InvActiveElements(:) => NULL()
      INTEGER :: NumberOfActiveElements
      INTEGER, ALLOCATABLE ::  Def_Dofs(:,:,:)

      TYPE(BlockMatrix_t), POINTER :: BlockMatrix => NULL()
      TYPE(Matrix_t),   POINTER :: Matrix => NULL()
      TYPE(Variable_t), POINTER :: Variable => NULL()

      TYPE(Matrix_t), POINTER :: ConstraintMatrix => NULL()
      TYPE(MortarBC_t), POINTER :: MortarBCs(:) => NULL()
      LOGICAL :: MortarBCsChanged = .FALSE., ConstraintMatrixVisited = .FALSE.
      INTEGER(KIND=AddrInt) :: MortarProc, &
          BoundaryElementProcedure=0, BulkElementProcedure=0

      TYPE(Graph_t), POINTER :: ColourIndexList => NULL(), BoundaryColourIndexList => NULL()
      INTEGER :: CurrentColour = 0, CurrentBoundaryColour = 0
      INTEGER :: DirectMethod = DIRECT_NORMAL
      LOGICAL :: GlobalBubbles = .FALSE., DG = .FALSE.
#ifdef USE_ISO_C_BINDINGS
      TYPE(C_PTR) :: CWrap = C_NULL_PTR
#endif
      TYPE(IntegrationPointsTable_t), POINTER :: IPTable => NULL()
    END TYPE Solver_t

!------------------------------------------------------------------------------


!-------------------Circuit stuff----------------------------------------------
  TYPE CircuitVariable_t
    LOGICAL :: isIvar, isVvar
    INTEGER :: BodyId, valueId, ImValueId, dofs, pdofs, Owner, ComponentId
    TYPE(Component_t), POINTER :: Component => Null()
    REAL(KIND=dp), ALLOCATABLE :: A(:), B(:)
    REAL(KIND=dp), ALLOCATABLE :: SourceRe(:), SourceIm(:), Mre(:), Mim(:)
    INTEGER, ALLOCATABLE :: EqVarIds(:)
  END TYPE CircuitVariable_t
  
  TYPE Component_t
    REAL(KIND=dp) :: Inductance=0._dp, Resistance=0._dp, Conductance = 0._dp, ElArea, &
                     N_j, coilthickness, i_multiplier_re, i_multiplier_im, nofturns
    INTEGER :: polord, nofcnts, BodyId, ComponentId
    INTEGER, POINTER :: ElBoundaries(:) => Null()
    INTEGER, POINTER :: BodyIds(:) => Null()
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    TYPE(CircuitVariable_t), POINTER :: ivar, vvar
    LOGICAL :: UseCoilResistance = .FALSE.
  END TYPE Component_t

  TYPE Circuit_t
    REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Mre(:,:), Mim(:,:), Area(:)
    INTEGER, ALLOCATABLE :: ComponentIds(:), Perm(:)
    LOGICAL :: UsePerm = .FALSE., Harmonic
    INTEGER :: n, m, n_comp,CvarDofs
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: names(:), source(:)
    TYPE(Component_t), POINTER :: Components(:)
    TYPE(CircuitVariable_t), POINTER :: CircuitVariables(:)
  END TYPE Circuit_t
!-------------------Circuit stuff----------------------------------------------

!------------------------------------------------------------------------------
    TYPE Model_t
!------------------------------------------------------------------------------
!
!     Coodrinate system dimension + type
!
      INTEGER :: DIMENSION, CoordinateSystem
!
!     Model dimensions
!
      INTEGER :: NumberOfBulkElements, &
                 NumberOfNodes,        &
                 NumberOfBoundaryElements
!
!     Simulation input data, that concern the model as a whole
!
      TYPE(ValueList_t), POINTER :: Simulation => Null()
!
!     Variables
!
      TYPE(Variable_t), POINTER  :: Variables => NULL()
!
!     Some physical constants, that will be read from the database or set by
!     other means: gravity direction/intensity and Stefan-Boltzmann constant)
!
      TYPE(ValueList_t), POINTER :: Constants => Null()
!
!     Types  of  equations (flow,heat,...) and  some  parameters (for example
!     laminar or turbulent flow or type of convection model for heat equation,
!     etc.)
!
      INTEGER :: NumberOfEquations = 0
      TYPE(EquationArray_t), POINTER :: Equations(:) => NULL()
!
!     Active electrical components
!
      INTEGER :: NumberOfComponents = 0
      TYPE(ComponentArray_t), POINTER :: Components(:) => NULL()
!
!     Active bodyforces: (bussinesq approx., heatsource, freele chosen
!     bodyforce...)
!
      INTEGER :: NumberOfBodyForces = 0
      TYPE(BodyForceArray_t), POINTER :: BodyForces(:) => NULL()
!
!     Initial conditions for field variables
!
      INTEGER :: NumberOfICs = 0
      TYPE(InitialConditionArray_t), POINTER :: ICs(:) => NULL()
!
!     Boundary conditions
!
      INTEGER :: NumberOfBCs = 0
      TYPE(BoundaryConditionArray_t), POINTER :: BCs(:) => NULL()
!
!     For free surfaces the curvatures...
!
      INTEGER, POINTER :: FreeSurfaceNodes(:) => NULL()
      REAL(KIND=dp), POINTER :: BoundaryCurvatures(:) => NULL()
!
!     Material parameters
!
      INTEGER :: NumberOfMaterials = 0
      TYPE(MaterialArray_t), POINTER :: Materials(:) => NULL()
!
!     Active bodies, every element has a pointer to a body, body has
!     material,ICs,bodyforces and equations
!
      INTEGER :: NumberOfBodies  = 0
      TYPE(BodyArray_t), POINTER :: Bodies(:) => NULL()
!
!      Boundary to boundary condition mapping
!
      INTEGER :: NumberOfBoundaries = 0
      INTEGER, POINTER :: BoundaryId(:) => NULL()
      TYPE(BoundaryArray_t), POINTER :: Boundaries(:) => NULL()
!
!     Linear equation solvers
!
      INTEGER :: NumberOfSolvers = 0
      TYPE(Solver_t), POINTER :: Solvers(:) => NULL()
!
!     Node coordinates + info for parallel computations
!
      TYPE(Nodes_t), POINTER :: Nodes => NULL()
!
!     Max number of nodes in any one element in this model
!
      INTEGER :: MaxElementNodes = 0
!
!     Elements
!
      TYPE(Element_t), POINTER :: Elements(:) => NULL()
!
!     For reference the current element in process   
!
      TYPE(Element_t), POINTER :: CurrentElement => NULL()
!
!     These are for internal use,   number of potentially nonzero elements
!     in stiffness and mass matrices (with one dof), and number of nonzero
!     elements in rows of the matrices.
!
      INTEGER :: TotalMatrixElements = 0
      INTEGER, POINTER :: RowNonzeros(:) => NULL()

      TYPE(Mesh_t), POINTER :: Meshes

      TYPE(Mesh_t),   POINTER :: Mesh   => NULL()
      TYPE(Solver_t), POINTER :: Solver => NULL()
      
      ! Circuits:
      INTEGER, POINTER :: n_Circuits=>Null(), Circuit_tot_n=>Null()
      TYPE(Matrix_t), POINTER :: CircuitMatrix => Null()
      TYPE(Circuit_t), POINTER :: Circuits(:) => Null()
      TYPE(Solver_t), POINTER :: ASolver    
      
      LOGICAL :: HarmonicCircuits

    END TYPE Model_t

    TYPE(Model_t),  POINTER :: CurrentModel
    TYPE(Matrix_t), POINTER :: GlobalMatrix

    INTEGER :: ELMER_COMM_WORLD = -1
!------------------------------------------------------------------------------
END MODULE Types
!------------------------------------------------------------------------------
!> \}
