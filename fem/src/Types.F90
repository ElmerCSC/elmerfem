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
 
!  USE Messages
   USE, INTRINSIC :: ISO_C_BINDING
#ifdef _OPENMP
   USE omp_lib 
#endif 

   USE Lua
   IMPLICIT NONE

   INTEGER, PARAMETER :: MAX_NAME_LEN = 128, MAX_STRING_LEN=2048, MAX_PATH_LEN=4096
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
                        SOLVER_EXEC_PREDCORR = 7,    &
                        SOLVER_EXEC_WHENCREATED = 8, &
                        SOLVER_EXEC_AFTER_CONTROL = 9
                        

  INTEGER, PARAMETER :: SOLVER_MODE_DEFAULT = 0, &    ! normal pde
                        SOLVER_MODE_AUXILIARY = 1, &  ! no fem machinery (SaveData)
                        SOLVER_MODE_ASSEMBLY = 2, &   ! coupled solver with single block
                        SOLVER_MODE_COUPLED = 3, &    ! coupled solver with multiple blocks
                        SOLVER_MODE_BLOCK = 4, &      ! block solver
                        SOLVER_MODE_GLOBAL = 5, &     ! lumped variables (no mesh)
                        SOLVER_MODE_MATRIXFREE = 6, & ! normal field, no matrix
                        SOLVER_MODE_STEPS = 7         ! as the legacy but split to different steps

  INTEGER, PARAMETER :: PROJECTOR_TYPE_DEFAULT = 0, &  ! unspecified constraint matrix
                        PROJECTOR_TYPE_NODAL = 1, &    ! nodal projector
                        PROJECTOR_TYPE_GALERKIN = 2, & ! Galerkin projector
                        PROJECTOR_TYPE_INTEGRAL = 3, & ! Integral type of constraint
                        PROJECTOR_TYPE_ROBIN = 4, &    ! Robin type of constraint
                        PROJECTOR_TYPE_NITSCHE = 5     ! Projector for Nitsche interface conditions
                                              
  INTEGER, PARAMETER :: DIRECT_NORMAL = 0, & ! Normal direct method
                        DIRECT_PERMON = 1    ! Permon direct method

  ! Operations used in ExchangeSourceVec
  INTEGER, PARAMETER :: OPER_SUM = 0, OPER_MIN = 1, OPER_MAX = 2, OPER_MEAN = 3
  
!------------------------------------------------------------------------------
  CHARACTER, PARAMETER :: Backslash = ACHAR(92)
!------------------------------------------------------------------------------


#ifdef HAVE_MUMPS
  INCLUDE 'smumps_struc.h'
  INCLUDE 'cmumps_struc.h'
  INCLUDE 'dmumps_struc.h'
  INCLUDE 'zmumps_struc.h'
#endif

  TYPE ArgStr_t
     CHARACTER(:), ALLOCATABLE :: astr
  END TYPE ArgStr_t


  TYPE IdxList_t
    INTEGER, ALLOCATABLE :: Ind(:)
  END TYPE IdxList_t


! Parallel collection state for the AMGX interface.
!
! AMGX wants whole owned rows on the owning rank, while Elmer keeps partial rows
! split by column ownership (see SplitMatrix in SParIterSolver), so AMGXSolver
! gathers the matrix itself into > Matrix % CollectionMatrix <. Everything that
! gather needs in order to repeat cheaply lives here, and lives with the matrix:
! it used to be held in SAVEd locals of AMGXSolver, which meant two solvers
! using AMGX on different matrices tore each other's state down on every call.
!
! The structure of the collected matrix is fixed once built, so the row and
! column indices are exchanged once and later solves ship only values:
!   LocalMap  entry of Matrix % Values -> slot in CollectionMatrix % Values,
!             zero for rows this partition does not own
!   SendIdx   per neighbour, the Matrix % Values entries to ship, in order
!   RecvIdx   per neighbour, the slots to add the arriving values into, in the
!             same order
  TYPE AMGXCollection_t
    INTEGER :: ng = 0                  !< global number of owned rows
    INTEGER :: nnzA = -1               !< structure of the source matrix the
    INTEGER :: nrowsA = -1             !<   cached pattern was built against
    LOGICAL :: PatternReady = .FALSE.
    INTEGER, ALLOCATABLE :: APerm(:), iLPerm(:), part_vec(:)
    INTEGER, ALLOCATABLE :: GlobalToLocal(:), SendTo(:)
    INTEGER, ALLOCATABLE :: LocalMap(:)
    TYPE(IdxList_t), ALLOCATABLE :: SendIdx(:), RecvIdx(:)
  END TYPE AMGXCollection_t


  TYPE BasicMatrix_t
    INTEGER :: NumberOfRows
    INTEGER, ALLOCATABLE :: Rows(:), Cols(:), Diag(:)
    INTEGER, ALLOCATABLE :: GRows(:), RowOwner(:)
    REAL(KIND=dp), ALLOCATABLE :: Values(:),MassValues(:), &
        DampValues(:),ILUValues(:),PrecValues(:)
!   Block view of a complex interface matrix, as for Matrix_t: NumberOfRows/2
!   block rows and one COMPLEX per 2x2 block. The local column index is folded
!   in from IfLCols, and rows this rank owns and entries with no local column
!   are dropped at build time, so the product needs no test at all.
    INTEGER, ALLOCATABLE :: BRows(:), BCols(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: CValues(:)
  END TYPE BasicMatrix_t


  TYPE SubVector_t
    TYPE(Variable_t), POINTER :: Var => NULL()
    REAL(KIND=dp) :: rnorm, bnorm, xnorm
    REAL(KIND=dp), ALLOCATABLE :: rhs(:)
    REAL(KIND=dp), ALLOCATABLE :: DiagScaling(:)
    TYPE(Solver_t), POINTER :: Solver => NULL()
    LOGICAL :: AddVector = .FALSE.
  END TYPE SubVector_t

  TYPE SubMatrix_t
    TYPE(Matrix_t), POINTER :: Mat => NULL()
    TYPE(Matrix_t), POINTER :: PrecMat => NULL()
    LOGICAL :: ParallelSquareMatrix = .TRUE.
    LOGICAL :: ParallelIsolatedMatrix = .FALSE.
    INTEGER, POINTER :: ParPerm(:) => NULL()
  END TYPE SubMatrix_t

  TYPE BlockMatrix_t
    INTEGER :: NoVar = 0, MaxSize, TotSize
    INTEGER, POINTER :: Offset(:) => NULL()
    INTEGER, POINTER :: ParOffset(:) => NULL()
    INTEGER, POINTER :: BlockPerm(:) => NULL()
    INTEGER, POINTER :: ParBlockPerm(:) => NULL()
    INTEGER, POINTER :: ParPerm(:) => NULL()
    TYPE(Solver_t), POINTER :: Solver => NULL()
    REAL(KIND=dp) :: rnorm, bnorm, xnorm
    TYPE(SubMatrix_t), ALLOCATABLE :: SubMatrix(:,:)
    LOGICAL, ALLOCATABLE :: SubMatrixActive(:,:)
    TYPE(SubVector_t), POINTER :: SubVector(:) => NULL()
    INTEGER, POINTER :: BlockStruct(:) => NULL()
    INTEGER, POINTER :: InvBlockStruct(:) => NULL()
    TYPE(Matrix_t), POINTER :: ParentMatrix => NULL()
    LOGICAL :: GotBlockStruct
    LOGICAL, ALLOCATABLE :: SubMatrixTranspose(:,:)
    INTEGER :: NoIters = 0
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


#ifdef HAVE_ROCALUTION
  TYPE Matrix_arr_t
     TYPE(Matrix_t), POINTER :: M
  END TYPE Matrix_arr_t

  TYPE RocParams_t
    TYPE(Matrix_t), POINTER :: Rmatrix => Null()
    TYPE(Matrix_arr_t), POINTER :: IMatrix(:) => Null()
    INTEGER, POINTER :: CntPerm(:)=> Null(), LocPerm(:) => Null(), gOffset(:) => Null()
  END TYPE RocParams_t
#endif


#ifdef HAVE_MUMPS
  !> Plan for turning MUMPS's distributed solution into Elmer's local rows.
  !> ICNTL(21)=1 leaves the solution spread over the ranks in a distribution
  !> MUMPS chooses, reported in ISOL_loc, which has nothing to do with the mesh
  !> partitioning. Redistributing it used to mean scattering into a vector of
  !> global length on every rank and allreducing that. It can instead be one
  !> MPI_ALLTOALLV, because ContinuousNumbering hands each rank a contiguous
  !> range of global indices for the rows it owns, so the owner of an index is
  !> a search in a table of range bounds rather than a lookup that would have
  !> to be communicated. Everything here is fixed by the factorization, so it
  !> is built once and reused by every solve that follows.
  TYPE MumpsSolPlan_t
    INTEGER :: nGlob = 0                     !< global rows, in plan units
    INTEGER :: gStart = 0                    !< our owned range is (gStart,gStart+nOwn]
    INTEGER :: nOwn = 0
    LOGICAL :: Built = .FALSE.               !< the per-solve plan below is ready
    INTEGER :: nSend = 0, nRecv = 0          !< slots we send / receive
    INTEGER, ALLOCATABLE :: gBound(:)        !< 0:PEs, owned index ranges
    INTEGER, ALLOCATABLE :: sCnt(:), sDsp(:) !< what we send, per rank
    INTEGER, ALLOCATABLE :: rCnt(:), rDsp(:) !< what we receive, per rank
    INTEGER, ALLOCATABLE :: sPerm(:)         !< sol_loc index of each send slot
    INTEGER, ALLOCATABLE :: rRow(:)          !< local row of each receive slot
    INTEGER, ALLOCATABLE :: isol(:)          !< the ISOL_loc it was built for
  END TYPE MumpsSolPlan_t
#endif

  TYPE Matrix_t
    TYPE(Matrix_t), POINTER :: Child => NULL(), Parent => NULL(), CircuitMatrix => NULL(), &
        ConstraintMatrix=>NULL(), EMatrix=>NULL(), AddMatrix=>NULL(), CollectionMatrix=>NULL()

    INTEGER :: NumberOfRows, ExtraDOFs=0, ParallelDOFs=0

    ! Number of degrees of freedom in sparse matrix such that there is always a (ndeg x ndeg) dense block
    INTEGER :: ndeg = -1
    
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

    REAL(KIND=dp), POINTER CONTIG :: RhsAdjoint(:)=>NULL()
    
    REAL(KIND=dp),  POINTER CONTIG :: Values(:)=>NULL(), ILUValues(:)=>NULL(), &
               DiagScaling(:) => NULL(), TValues(:) => NULL(), Values_im(:) => NULL()

    REAL(KIND=dp), ALLOCATABLE :: extraVals(:)
    REAL(KIND=dp) :: RhsScaling=1.0_dp, AveScaling=1.0_dp 
    INTEGER :: ScalingMethod = 0
    REAL(KIND=dp),  POINTER CONTIG :: MassValues(:)=>NULL(),DampValues(:)=>NULL(), &
        BulkValues(:)=>NULL(), BulkMassValues(:)=>NULL(), BulkDampValues(:)=>NULL(), &
        PrecValues(:)=>NULL(), HaloValues(:)=>NULL(), HaloMassValues(:)=>NULL()

#ifdef HAVE_FETI4I
    TYPE(C_PTR) :: PermonMatrix = C_NULL_PTR, PermonSolverInstance = C_NULL_PTR
#endif
#ifdef HAVE_MUMPS
    TYPE(dmumps_struc), POINTER :: MumpsID => NULL() ! Global distributed Mumps
    TYPE(dmumps_struc), POINTER :: MumpsIDL => NULL() ! Local domainwise Mumps

    TYPE(zmumps_struc), POINTER :: ZMumpsID => NULL() ! Global distributed Mumps
    TYPE(zmumps_struc), POINTER :: ZMumpsIDL => NULL() ! Local domainwise Mumps

    TYPE(smumps_struc), POINTER :: SMumpsID => NULL() ! Global distributed Mumps
    TYPE(smumps_struc), POINTER :: SMumpsIDL => NULL() ! Local domainwise Mumps

    TYPE(cmumps_struc), POINTER :: CMumpsID => NULL() ! Global distributed Mumps
    TYPE(cmumps_struc), POINTER :: CMumpsIDL => NULL() ! Local domainwise Mumps

    ! Redistribution plan for the distributed solution, valid as long as the
    ! factorization above is. Released by FreeMumpsFactorizations.
    TYPE(MumpsSolPlan_t), POINTER :: MumpsSolPlan => NULL()
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
#ifdef HAVE_ROCALUTION
    TYPE(RocParams_t) :: RocParams
#endif
    INTEGER(KIND=C_INTPTR_T) :: AMGX=0, AMGXMV=0
    TYPE(AMGXCollection_t), POINTER :: AMGXColl => NULL()
    INTEGER(KIND=AddrInt) :: SpMV=0

    TYPE(C_FUNPTR) :: MatVecSubr = C_NULL_FUNPTR

    INTEGER, POINTER CONTIG :: ILURows(:)=>NULL(),ILUCols(:)=>NULL(),ILUDiag(:)=>NULL()

!   For Complex systems, not used yet!:
!   -----------------------------------
    COMPLEX(KIND=dp), POINTER :: CRHS(:)=>NULL(),CForce(:,:)=>NULL()
    COMPLEX(KIND=dp),  POINTER :: CValues(:)=>NULL(),CILUValues(:)=>NULL()
    COMPLEX(KIND=dp),  POINTER :: CMassValues(:)=>NULL(),CDampValues(:)=>NULL()

!   Optional block CRS view of a complex matrix: NumberOfRows/2 block rows and
!   one COMPLEX coefficient per 2x2 [Re -Im; Im Re] block, in BRows/BCols and
!   the CValues slot above. Built on demand by CRS_BuildBlockCRS, used by
!   CRS_BlockComplexMatrixVectorMultiply. Not present unless asked for.
!
!   CPrecValues is the same view over PrecValues when a separate preconditioning
!   matrix exists. It shares BRows and BCols outright: DefaultUpdatePrecC glues
!   the preconditioner through the same structure as the matrix, so the only
!   thing that differs is the coefficients. The complex ILU factorizes
!   PrecValues when it is present, so without this the view is bypassed exactly
!   on the cases that have one -- which includes every VectorHelmholtz case
!   carrying a damping coefficient.
    INTEGER, POINTER :: BRows(:)=>NULL(), BCols(:)=>NULL(), BDiag(:)=>NULL()
    COMPLEX(KIND=dp), POINTER :: CPrecValues(:)=>NULL()

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
     INTEGER                          :: PEs  = 1
     INTEGER                          :: MyPE = 0
     LOGICAL                          :: Initialized = .FALSE.
     INTEGER                          :: ActiveComm  = 0
     LOGICAL, DIMENSION(:), POINTER   :: Active => Null()
     LOGICAL, DIMENSION(:), POINTER   :: IsNeighbour => Null()
     INTEGER                          :: NumOfNeighbours = 0
     INTEGER                          :: NumberOfThreads = 1
     LOGICAL                          :: ExternalInit = .FALSE.
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


  TYPE RealBuf_t
    REAL(KIND=dp), ALLOCATABLE :: rbuf(:)
  END TYPE RealBuf_t

  ! The gather that fills the send buffer of one neighbour in SParMatrixVector.
  ! Which interface row ends up in which slot of the buffer is fixed by the
  ! matrix structure, so it is worked out once instead of by rescanning every
  ! RowOwner array on every product. The slots are grouped into segments, one
  ! per contributing interface block, leaving one pointer dereference per
  ! segment and one indexed load per slot.
  TYPE MVPackT
    INTEGER :: nseg = 0
    INTEGER, ALLOCATABLE :: SegIf(:)     ! interface block a segment reads from
    INTEGER, ALLOCATABLE :: SegStart(:)  ! first slot of a segment, size nseg+1
    INTEGER, ALLOCATABLE :: Row(:)       ! row of that block, one per slot
  END TYPE MVPackT

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
     ! Persistent MPI communication buffers for SParMatrixVector (allocated once at setup)
     INTEGER, ALLOCATABLE :: MVNeigh(:)
     INTEGER, ALLOCATABLE :: MVSendSize(:)
     INTEGER, ALLOCATABLE :: MVRecvSize(:)
     TYPE(RealBuf_t), ALLOCATABLE :: MVSendBuf(:)
     TYPE(RealBuf_t), ALLOCATABLE :: MVRecvBuf(:)
     INTEGER, ALLOCATABLE :: MVRequests(:)
     INTEGER, ALLOCATABLE :: MVSendRequests(:)
     TYPE(MVPackT), ALLOCATABLE :: MVPack(:)
  END TYPE SplittedMatrixT


  TYPE SParIterSolverGlobalD_t
     TYPE (SplittedMatrixT), POINTER :: SplittedMatrix=>NULL()
     TYPE (Matrix_t), POINTER :: Matrix=>NULL()
     INTEGER :: DOFs, RelaxIters
     ! The active partitions and the neighbours are a property of the matrix,
     ! not of the solver: one solver may own several matrices of differing
     ! parallel structure. This is the owner of the > Active < and
     ! > IsNeighbour < arrays, > Solver % ParEnv < only mirrors it.
     TYPE(ParEnv_t) :: ParEnv
     TYPE (ParallelInfo_t), POINTER :: ParallelInfo=>NULL()
  END TYPE SParIterSolverGlobalD_t

  TYPE(SParIterSolverGlobalD_t), POINTER :: ParMatrix => NULL()

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
   INTEGER, PARAMETER :: ELEM_BASIS_CACHE_SIZE = 64

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

     TYPE(BasisFunctions_t), POINTER :: BasisFunctions(:) => NULL()
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: NodeU, NodeV, NodeW
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: P_NodeU, P_NodeV, P_NodeW
     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: N_NodeU, N_NodeV, N_NodeW
     ! Reference basis function cache — keyed by (u,v,w). This structure is
     ! SHARED by every element of this type across all OpenMP threads, so the
     ! cache is given a per-thread column (last index) to keep both lookup and
     ! fill lock-free and race-free during threaded assembly (see ElemInfo's
     ! ElementInfo). BasisCacheCount(tid)=0 means that thread's cache is empty.
     INTEGER, ALLOCATABLE :: BasisCacheCount(:)          ! (nthreads)
     REAL(KIND=dp), ALLOCATABLE :: BasisCacheU(:,:)      ! (ELEM_BASIS_CACHE_SIZE, nthreads)
     REAL(KIND=dp), ALLOCATABLE :: BasisCacheV(:,:)
     REAL(KIND=dp), ALLOCATABLE :: BasisCacheW(:,:)
     REAL(KIND=dp), ALLOCATABLE :: BasisCache(:,:,:)     ! (ELEM_BASIS_CACHE_SIZE, n_nodes, nthreads)
     REAL(KIND=dp), ALLOCATABLE :: dBasisCache(:,:,:,:)  ! (ELEM_BASIS_CACHE_SIZE, n_nodes, 3, nthreads)
   END TYPE ElementType_t

!------------------------------------------------------------------------------
   TYPE ValueListEntry_t
     INTEGER :: Type
     TYPE(ValueListEntry_t), POINTER :: Next => NULL()

     REAL(KIND=dp), POINTER :: TValues(:) => NULL()
     REAL(KIND=dp), POINTER :: Cumulative(:) => NULL()
     REAL(KIND=dp), POINTER :: FValues(:,:,:) => NULL()
     REAL(KIND=dp), POINTER :: CubicCoeff(:) => NULL()
     INTEGER :: Fdim = 0 
     
     LOGICAL :: LValue
     INTEGER, POINTER :: IValues(:) => NULL()

     TYPE(C_FUNPTR) :: PROCEDURE = C_NULL_FUNPTR

     REAL(KIND=dp) :: Coeff = 1.0_dp    
     CHARACTER(:), ALLOCATABLE :: CValue

     INTEGER :: NameLen,DepNameLen = 0
     CHARACTER(:), ALLOCATABLE :: Name,DependName

#ifdef DEVEL_LISTCOUNTER 
     INTEGER :: Counter = 0
#endif
#ifdef DEVEL_LISTUSAGE
     INTEGER :: Counter = 0
#endif

     LOGICAL :: LuaFun = .FALSE.
     INTEGER :: partag = 0
     LOGICAL :: disttag = .FALSE.
   END TYPE ValueListEntry_t

   TYPE ValueList_t
     TYPE(ValueListEntry_t), POINTER :: Head => NULL()
   END TYPE ValueList_t

   
   TYPE VariableTable_t     
     TYPE(Variable_t), POINTER :: Variable => NULL()
     TYPE(ValueListEntry_t), POINTER :: Keyword => NULL()
     REAL(KIND=dp) :: ParamValue 
     INTEGER :: tstep = 0
   END TYPE VariableTable_t

   
   ! This is a tentative data type to speed up the retrieval of parameters
   ! at elements.
   !----------------------------------------------------------------------
   TYPE ValueHandle_t
     INTEGER :: ValueType = -1
     INTEGER :: SectionType = -1
     INTEGER :: ListId = -9999
     LOGICAL :: BulkElement
     TYPE(Element_t), POINTER :: Element => NULL()
     TYPE(ValueList_t), POINTER :: List => NULL()
     TYPE(ValueList_t), POINTER :: Ptr  => NULL()
     TYPE(Nodes_t), POINTER :: Nodes => NULL()
     INTEGER, POINTER :: Indexes => NULL()
     INTEGER :: n
     INTEGER :: nValuesVec = 0
     REAL(KIND=dp), POINTER :: ValuesVec(:) => NULL()
     REAL(KIND=dp), POINTER :: Values(:) => NULL()
     REAL(KIND=dp), POINTER :: ParValues(:,:) => NULL()
     LOGICAL, POINTER :: ParUsed(:) => NULL()
     INTEGER :: ParNo = 0
     INTEGER :: IValue, DefIValue = 0
     REAL(KIND=dp) :: RValue, DefRValue = 0.0_dp
     INTEGER :: Rdim = 0
     REAL(KIND=dp), POINTER :: RTensor(:,:) => NULL()
     REAL(KIND=dp), POINTER :: RTensorValues(:,:,:) => NULL()
     LOGICAL :: LValue, DefLValue = .FALSE.
     CHARACTER(:), ALLOCATABLE :: CValue
     INTEGER :: CValueLen
     LOGICAL :: Found
     CHARACTER(:), ALLOCATABLE :: Name
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
     INTEGER :: VarCount = 0
     INTEGER :: IntVarCount = 0
     TYPE(ValueHandle_t), POINTER :: HandleIm => NULL()
     TYPE(ValueHandle_t), POINTER :: Handle2 => NULL()
     TYPE(ValueHandle_t), POINTER :: Handle3 => NULL()
   END TYPE ValueHandle_t


   TYPE VariableHandle_t     
     TYPE(Variable_t), POINTER :: Variable=>NULL()
     REAL(KIND=dp),POINTER :: Values(:)=>NULL()
     REAL(KIND=dp),POINTER :: ipValues(:)=>NULL()
     REAL(KIND=dp),POINTER :: ipValues3D(:,:)=>NULL()     
     INTEGER :: ipN = 0     
     INTEGER,POINTER :: Perm(:)=>NULL()
     INTEGER :: dofs
     INTEGER :: tstep = 0
     TYPE(Element_t), POINTER :: Element=>NULL()
     LOGICAL :: ActiveElement = .FALSE.
     LOGICAL :: Found
     INTEGER :: Indexes(100)     
     INTEGER :: n = 0
   END TYPE VariableHandle_t
   
   
!------------------------------------------------------------------------------

   TYPE MaterialArray_t
     TYPE(ValueList_t), POINTER :: Values => NULL()
   END TYPE MaterialArray_t

!------------------------------------------------------------------------------

   TYPE BoundaryConditionArray_t
     INTEGER :: Tag=0
     TYPE(Matrix_t), POINTER :: PMatrix => NULL()
     TYPE(ValueList_t), POINTER :: Values => NULL()
   END TYPE BoundaryConditionArray_t

!------------------------------------------------------------------------------

   TYPE InitialConditionArray_t
     INTEGER :: Tag=0
     TYPE(ValueList_t), POINTER :: Values => NULL()
   END TYPE InitialConditionArray_t

!------------------------------------------------------------------------------

    TYPE ComponentArray_t
      TYPE(ValueList_t), POINTER :: Values => NULL()
    END TYPE ComponentArray_t

!------------------------------------------------------------------------------

    TYPE BodyForceArray_t
      TYPE(ValueList_t), POINTER :: Values => NULL()
    END TYPE BodyForceArray_t

!------------------------------------------------------------------------------

    TYPE BoundaryArray_t
      TYPE(ValueList_t), POINTER :: Values => NULL()
    END TYPE BoundaryArray_t

!------------------------------------------------------------------------------

    TYPE BodyArray_t
      TYPE(ValueList_t), POINTER :: Values => NULL()
    END TYPE BodyArray_t

!------------------------------------------------------------------------------

    TYPE EquationArray_t
      TYPE(ValueList_t), POINTER :: Values => NULL()
    END TYPE EquationArray_t

!------------------------------------------------------------------------------

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
     INTEGER, POINTER :: IPOffset(:) => NULL()
     !TYPE(GaussIntegrationPoints_t), POINTER :: IPs
   END TYPE IntegrationPointsTable_t
      
   
   TYPE Variable_t
     TYPE(Variable_t), POINTER :: Next => NULL()
     TYPE(Variable_t), POINTER :: EVar => NULL() 
     INTEGER :: NameLen = 0
     CHARACTER(:), ALLOCATABLE :: Name

     TYPE(Solver_t), POINTER :: Solver => NULL()
     LOGICAL :: Valid = .TRUE.
     LOGICAL :: Output = .TRUE.
     TYPE(Mesh_t), POINTER :: PrimaryMesh => NULL()

     LOGICAL :: ValuesChanged = .TRUE.
     LOGICAL :: DgAveraged = .FALSE.
     
! Some variables are created from pointers to the primary variables
     LOGICAL :: Secondary = .FALSE.

     INTEGER :: TYPE = Variable_on_nodes

     INTEGER :: DOFs = 0
     INTEGER, POINTER          :: Perm(:) => NULL()
     LOGICAL :: PeriodicFlipActive = .FALSE.
     REAL(KIND=dp)             :: Norm=0.0_dp, PrevNorm=0.0_dp,&
         NonlinChange=0.0_dp, SteadyChange=0.0_dp
     INTEGER :: NonlinConverged=-1, SteadyConverged=-1, NonlinIter=-1
     INTEGER :: LinConverged=-1
     COMPLEX(KIND=dp), POINTER :: EigenValues(:) => NULL(), &
          EigenVectors(:,:) => NULL()
     REAL(KIND=dp), POINTER :: ConstraintModes(:,:) => NULL()
     INTEGER, POINTER :: ConstraintModesIndeces(:) => NULL()
     REAL(KIND=dp), POINTER :: ConstraintModesWeights(:) => NULL()
     INTEGER :: NumberOfConstraintModes = -1
     LOGICAL :: FrozenMode = .FALSE.
     REAL(KIND=dp), POINTER :: Values(:) => NULL() ,&
          PrevValues(:,:) => NULL(), &
          PValues(:) => NULL(), NonlinValues(:) => NULL(), &
          SteadyValues(:) => NULL(), DeltaValues(:) => NULL()
     LOGICAL, POINTER :: UpperLimitActive(:) => NULL(), LowerLimitActive(:) => NULL()
     REAL(KIND=dp), POINTER :: UpperLimit(:) => NULL(), LowerLimit(:) => NULL()
     COMPLEX(KIND=dp), POINTER :: CValues(:) => NULL()
     TYPE(IntegrationPointsTable_t), POINTER :: IPTable => NULL()
   END TYPE Variable_t

!------------------------------------------------------------------------------
   TYPE ListMatrixEntry_t
     INTEGER :: Index = -1
     REAL(KIND=dp) :: val = 0.0_dp
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
     TYPE(ListMatrixEntry_t), POINTER :: Head => NULL()
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
     INTEGER, ALLOCATABLE :: Elements(:)
     REAL(KIND=dp), ALLOCATABLE :: Factors(:)
   END TYPE Factors_t

!-------------------------------------------------------------------------------

   TYPE BoundaryInfo_t
     TYPE(Factors_t), POINTER :: RadiationFactors => NULL()
     INTEGER :: Constraint = 0, OutBody = -1
     REAL(KIND=dp), ALLOCATABLE :: Radiators(:)
     TYPE(Element_t), POINTER :: Left =>NULL(), Right=>NULL()
   END TYPE BoundaryInfo_t

!-------------------------------------------------------------------------------

   TYPE ElementData_t
     TYPE(ElementData_t), POINTER :: Next=>NULL()
     CHARACTER(:), ALLOCATABLE :: Name
     REAL(KIND=dp), POINTER :: Values(:)=>NULL()
   END TYPE ElementData_t

!-------------------------------------------------------------------------------

   TYPE Element_t
     TYPE(ElementType_t), POINTER :: TYPE => NULL()

     LOGICAL :: Copy = .FALSE.

     INTEGER :: BodyId=0, Splitted=0, Status=0
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
      LOGICAL :: Serendipity=.TRUE.     ! Is this a serendipity of complete basis element ?
      INTEGER :: localNumber     ! Local number of an edge or face for element on boundary
      TYPE(Element_t), POINTER :: localParent => NULL()     ! Local number of an edge or face for element on boundary
   END TYPE PElementDefs_t

!-------------------------------------------------------------------------------

   TYPE NeighbourList_t
     INTEGER, DIMENSION(:), POINTER :: Neighbours=>NULL()
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
     TYPE(Quadrant_t), POINTER :: Quadrant=>NULL()
   END TYPE QuadrantPointer_t

!------------------------------------------------------------------------------

   TYPE Quadrant_t
     INTEGER, DIMENSION(:), POINTER :: Elements
     REAL(KIND=dp) :: SIZE, MinElementSize, BoundingBox(6)
     INTEGER :: NElemsInQuadrant
     TYPE(QuadrantPointer_t), DIMENSION(:), POINTER :: ChildQuadrants=>NULL()
   END TYPE Quadrant_t

!------------------------------------------------------------------------------

   TYPE Projector_t
     TYPE(Projector_t), POINTER :: Next=>NULL()
     TYPE(Mesh_t), POINTER :: Mesh=>NULL()
     TYPE(Matrix_t), POINTER :: Matrix=>NULL(), TMatrix=>NULL()
   END TYPE Projector_t


!------------------------------------------------------------------------------

   TYPE ParallelInfo_t
     INTEGER :: NumberOfIfDOFs
     LOGICAL, POINTER               :: GInterface(:) => NULL()
     INTEGER, POINTER               :: GlobalDOFs(:) => NULL()
     TYPE(NeighbourList_t),POINTER  :: NeighbourList(:) => NULL()
     INTEGER, POINTER               :: Gorder(:) => NULL()

     LOGICAL, POINTER               :: FaceInterface(:) => NULL()
     TYPE(NeighbourList_t),POINTER  :: FaceNeighbourList(:) => NULL()

     LOGICAL, POINTER               :: EdgeInterface(:) => NULL()
     TYPE(NeighbourList_t),POINTER  :: EdgeNeighbourList(:) => NULL()
     LOGICAL                        :: NothingShared = .FALSE.
   END TYPE ParallelInfo_t

!------------------------------------------------------------------------------

   TYPE NormalTangential_t     
     CHARACTER(:), ALLOCATABLE :: NormalTangentialName
     INTEGER :: NormalTangentialNOFNodes = 0
     INTEGER, POINTER :: BoundaryReorder(:) => NULL()
     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:)  => NULL()
     REAL(KIND=dp), POINTER :: BoundaryTangent1(:,:) => NULL()
     REAL(KIND=dp), POINTER :: BoundaryTangent2(:,:) => NULL()
   END TYPE NormalTangential_t

   TYPE FactorsStore_t
     TYPE(Factors_t), POINTER :: VF(:) => NULL()
   END TYPE FactorsStore_t 

   TYPE Mesh_t
     CHARACTER(:), ALLOCATABLE :: Name
     TYPE(Mesh_t), POINTER   :: Next => NULL(), Parent => NULL(), Child => NULL()

     TYPE(Projector_t), POINTER :: Projector => NULL()
     TYPE(Quadrant_t), POINTER  :: RootQuadrant => NULL()

     LOGICAL :: Changed, OutputActive, Stabilize = .FALSE.
     INTEGER :: SavesDone, AdaptiveDepth, MeshTag = 1
     LOGICAL :: AdaptiveFinished = .FALSE.

     TYPE(Factors_t), POINTER :: ViewFactors(:)=>NULL()
     TYPE(FactorsStore_t), ALLOCATABLE :: VFStore(:)

     TYPE(ParallelInfo_t) :: ParallelInfo
     TYPE(Variable_t), POINTER :: Variables => NULL()

     TYPE(Nodes_t), POINTER :: Nodes => NULL()
     TYPE(Element_t), DIMENSION(:), POINTER :: Elements => NULL(), &
         Edges => NULL(), Faces => NULL()
     TYPE(Nodes_t), POINTER :: NodesOrig => NULL()
     TYPE(Nodes_t), POINTER :: NodesMapped => NULL()

     INTEGER :: SolverId = 0
     
     LOGICAL :: DisContMesh = .FALSE.  
     INTEGER, POINTER :: DisContPerm(:) => NULL()
     INTEGER :: DisContNodes = 0

     INTEGER, POINTER :: PeriodicPerm(:) => NULL()
     LOGICAL, POINTER :: PeriodicFlip(:) => NULL()
     
     INTEGER, POINTER :: InvPerm(:) => NULL()

     INTEGER :: NumberOfNodes, NumberOfBulkElements, NumberOfEdges, &
                NumberOfFaces, NumberOfBoundaryElements, MeshDim = 0, MaxDim = 0, PassBCcnt=0
     INTEGER :: MinEdgeDOFs, MinFaceDOFs
     INTEGER :: MaxElementNodes, MaxElementDOFs, MaxEdgeDOFs, MaxFaceDOFs, MaxBDOFs
     INTEGER :: MaxNDOFs ! The maximum of nodal DOFs per node (created with a flag "n:")

     LOGICAL :: EntityWeightsComputed 
     REAL(KIND=dp), POINTER :: BCWeight(:) => NULL(), BodyForceWeight(:) => NULL(),&
         BodyWeight(:) => NULL(), MaterialWeight(:) => NULL()

     INTEGER, POINTER :: RePartition(:) => NULL()
     TYPE(NeighbourList_t), POINTER :: Halo(:) => NULL()
     LOGICAL :: HaveHalo = .FALSE.

     LOGICAL :: SingleMesh = .FALSE.
   END TYPE Mesh_t

   TYPE Graph_t
     INTEGER :: n
     INTEGER, ALLOCATABLE :: ptr(:), ind(:)
   END type Graph_t
   
   TYPE Graphcolour_t
     INTEGER :: nc
     INTEGER, POINTER :: colours(:) => Null()
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
     INTEGER :: RowOffset = 0
   END TYPE MortarBC_t

   TYPE TabulatedBasisAtIp_t
     REAL(KIND=dp), POINTER :: Basis(:) => NULL()
     REAL(KIND=dp), POINTER :: dBasisdx(:,:) => NULL()
     REAL(KIND=dp) :: Weight = 0.0_dp
   END TYPE TabulatedBasisAtIp_t

   TYPE LumpedModel_t
     LOGICAL :: IsComplex = .FALSE.
     INTEGER :: CurrentRow = -1
     INTEGER :: NoModes = 0
     INTEGER :: CntModes = 0
     REAL(KIND=dp), POINTER :: CMatrix(:,:) => NULL()
     REAL(KIND=dp), POINTER :: CMatrixIm(:,:) => NULL()                
     REAL(KIND=dp), POINTER :: Crhs(:) => NULL()
     REAL(KIND=dp), POINTER :: CrhsIm(:) => NULL()
     REAL(KIND=dp), POINTER :: ImpRe(:) => NULL()
     REAL(KIND=dp), POINTER :: ImpIm(:) => NULL()     
   END TYPE LumpedModel_t


   TYPE LocalSystemStorage_t
     INTEGER :: n = -1                            ! Size of local matrix assembled
     REAL(KIND=dp), ALLOCATABLE ::  K(:,:), F(:)  ! Local stiffness matrix and force
     INTEGER :: eind = -1                         ! Pointer to the active element that holds this element matrix storage.
   END TYPE LocalSystemStorage_T
   
   
!------------------------------------------------------------------------------

    TYPE Solver_t
      INTEGER :: SolverId = 0
      TYPE(ValueList_t), POINTER :: Values => NULL()

      INTEGER :: TimeOrder=0,DoneTime=0,Order=0,NOFEigenValues=0
      INTEGER :: TimesVisited = 0
      TYPE(C_FUNPTR) :: PROCEDURE = C_NULL_FUNPTR, LinBeforeProc = C_NULL_FUNPTR, LinAfterProc = C_NULL_FUNPTR

      REAL(KIND=dp) :: Alpha,Beta,dt

      REAL(KIND=dp) :: AitkenRelax = 1.0_dp
      
      LOGICAL :: NewtonActive = .FALSE.
      LOGICAL :: PeriodicFlipActive = .FALSE.
      
      INTEGER :: SolverExecWhen=-1
      INTEGER :: SolverMode=-1

      INTEGER :: MultiGridLevel=-1, MultiGridTotal=-1, MultiGridSweep=-1
      LOGICAL :: MultiGridSolver=.FALSE., MultiGridEqualSplit=.FALSE.
      TYPE(Mesh_t), POINTER :: Mesh => NULL()
      INTEGER :: MeshTag = 1
      LOGICAL :: MeshChanged = .FALSE.
      
      INTEGER, POINTER :: ActiveElements(:) => NULL()
      INTEGER, POINTER :: InvActiveElements(:) => NULL()
      INTEGER :: NumberOfActiveElements=0
      INTEGER, ALLOCATABLE ::  Def_Dofs(:,:,:)

      TYPE(BlockMatrix_t), POINTER :: BlockMatrix => NULL()
      TYPE(Matrix_t),   POINTER :: Matrix => NULL()
      TYPE(Variable_t), POINTER :: Variable => NULL()

      TYPE(Matrix_t), POINTER :: ConstraintMatrix => NULL()
      TYPE(MortarBC_t), POINTER :: MortarBCs(:) => NULL()
      LOGICAL :: MortarBCsChanged = .FALSE., ConstraintMatrixVisited = .FALSE.
      TYPE(C_FUNPTR) :: BoundaryElementProcedure = C_NULL_FUNPTR, BulkElementProcedure = C_NULL_FUNPTR

      TYPE(Graph_t), POINTER :: ColourIndexList => NULL()
      TYPE(Graph_t), POINTER :: BoundaryColourIndexList => NULL()
      INTEGER :: CurrentColour = 0, CurrentBoundaryColour = 0
      INTEGER :: DirectMethod = DIRECT_NORMAL
      LOGICAL :: GlobalBubbles = .FALSE.
      LOGICAL :: DG = .FALSE.
      TYPE(C_PTR) :: CWrap = C_NULL_PTR
      TYPE(IntegrationPointsTable_t), POINTER :: IPTable => NULL()
      LOGICAL :: Parallel = .FALSE.

      TYPE(NormalTangential_t) :: NormalTangential

      INTEGER :: NumberOfConstraintModes = -1 
      TYPE(LumpedModel_t), POINTER :: Lumped => NULL()

      INTEGER :: LocalSystemMode = -1
      TYPE(LocalSystemStorage_t), POINTER :: LocalSystem(:) => NULL()

      TYPE(ParEnv_t) :: ParEnv
      REAL(KIND=dp), POINTER :: CutInterp(:) => NULL()
    END TYPE Solver_t

!------------------------------------------------------------------------------


!-------------------Circuit stuff----------------------------------------------
  TYPE CircuitVariable_t
    ! Default initialized for the same reason as Component_t below. ImValueId in
    ! particular is only assigned for a harmonic circuit, while
    ! GetComponentCurrent() tests it with "> 0" for any circuit - so in a
    ! transient run it used to decide on uninitialized memory, and when that
    ! happened to be a valid index it read the real current as the imaginary one.
    LOGICAL :: isIvar = .FALSE., isVvar = .FALSE.
    INTEGER :: BodyId = 0, valueId = 0, ImValueId = 0, dofs = 0, pdofs = 0, &
         Owner = 0, ComponentId = 0
    TYPE(Component_t), POINTER :: Component => NULL()
    REAL(KIND=dp), ALLOCATABLE :: A(:), B(:)
    REAL(KIND=dp), ALLOCATABLE :: SourceRe(:), SourceIm(:), Mre(:), Mim(:)
    INTEGER, ALLOCATABLE :: EqVarIds(:)
  END TYPE CircuitVariable_t
  
  TYPE Component_t
    ! Every field is default initialized on purpose. ReadComponents() only fills
    ! the ones that its coil type needs - ElArea and N_j for instance are set for
    ! 'stranded' and 'foil winding' but not for 'massive' or for a resistor - yet
    ! consumers such as GetComponentArea() and the 'Stranded Coil N_j' keyword
    ! read them for any component. Without the initializers those reads returned
    ! whatever was on the heap. Zero also doubles as a detectable "not set".
    REAL(KIND=dp) :: Inductance=0._dp, Resistance=0._dp, Conductance = 0._dp, &
         ElArea=0._dp, N_j=0._dp, coilthickness=1._dp, i_multiplier_re=0._dp, &
         i_multiplier_im=0._dp, nofturns=0._dp, &
         VoltageFactor=1._dp, SymmetryCoeff=1._dp
    INTEGER :: polord=0, nofcnts=0, BodyId=0, ComponentId=0
    INTEGER, POINTER :: ElBoundaries(:) => NULL()
    INTEGER, POINTER :: BodyIds(:) => NULL()
    ! Indices into the circuit solver's active and boundary element lists for the
    ! elements belonging to this component, built once by
    ! BuildComponentElementLists(). Ascending, so that walking them backwards
    ! reproduces the "DO q=GetNOFActive(),1,-1" order the assembly used when it
    ! rescanned every element for every component.
    INTEGER, ALLOCATABLE :: ElemIdx(:), BCElemIdx(:)
    CHARACTER(:), ALLOCATABLE :: CoilType, ComponentType
    TYPE(CircuitVariable_t), POINTER :: ivar=>NULL(), vvar=>NULL()
    LOGICAL :: UseCoilResistance = .FALSE.
  END TYPE Component_t

!-------------------Circuit stuff----------------------------------------------
  TYPE Circuit_t
    REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Mre(:,:), Mim(:,:), Area(:)
    INTEGER, ALLOCATABLE :: ComponentIds(:), Perm(:)
    LOGICAL :: UsePerm = .FALSE., Harmonic=.FALSE., Parallel=.FALSE.
    INTEGER :: n=0, m=0, n_comp=0,CvarDofs=0
    CHARACTER(MAX_NAME_LEN), ALLOCATABLE :: names(:), source(:)
    TYPE(Component_t), POINTER :: Components(:)=>NULL()
    TYPE(CircuitVariable_t), POINTER :: CircuitVariables(:)=>NULL()
    TYPE(Solver_t), POINTER :: ASolver => NULL()
  END TYPE Circuit_t

!> Everything one circuit solver owns.
!>
!> This used to be a set of slots directly on Model_t, plus a handful of module
!> level SAVE variables in CircuitUtils, which limited a run to one circuit
!> solver: a second one would overwrite the first one's circuits, matrix and
!> build record. Collecting it here makes the state per instance, so that
!> Model % CircuitModels can hold one entry per circuit solver.
!>
!> Instances are addressed through Model % CircuitModel, which the circuit
!> solver entry points set to their own container. See GetCircuitModel() and
!> SetCircuitModel() in CircuitUtils.
  TYPE CircuitModel_t
    ! Which solvers this instance belongs to. Solver is the circuit solver that
    ! owns the container, ASolver the FEM solver whose matrix the circuit
    ! equations are appended to. Both are borrowed pointers into Model % Solvers.
    INTEGER :: SolverId = 0
    TYPE(Solver_t), POINTER :: Solver => NULL()
    TYPE(Solver_t), POINTER :: ASolver => NULL()

    ! The circuit equations themselves.
    INTEGER :: n_Circuits = 0, Circuit_tot_n = 0
    TYPE(Circuit_t), POINTER :: Circuits(:) => NULL()
    TYPE(Matrix_t), POINTER :: CircuitMatrix => NULL()
    LOGICAL :: Harmonic = .FALSE.

    ! Cache invalidation. Generation is a ticket drawn from a module wide
    ! counter in CircuitUtils, so it is unique over all instances and over all
    ! rebuilds. Routines that cache something behind a "first time through" test
    ! compare a saved copy of it against CircuitsGeneration(), and therefore
    ! re-derive both when this instance is rebuilt and when the active instance
    ! changes under them. BuiltNm/BuiltTotN/BuiltMesh record what the structures
    ! were built against; see CircuitsCheckStale().
    INTEGER :: Generation = 0
    INTEGER :: BuiltNm = -1, BuiltTotN = -1
    TYPE(Mesh_t), POINTER :: BuiltMesh => NULL()

    ! Prefix for the MATC symbols this instance reads its definitions from
    ! ("Circuits", "C.<p>.*"). Empty is the historical global namespace.
    CHARACTER(:), ALLOCATABLE :: MatcPrefix

    ! Driver state that used to be SAVEd inside the circuit solver entry points,
    ! and so was shared by every instance. Crt holds the circuit variable values
    ! of the previous timestep, sized Circuit_tot_n; Tstep is the timestep it was
    ! last refreshed on; MultName is the name of the A solver's Lagrange
    ! multiplier variable the values are read from.
    REAL(KIND=dp), ALLOCATABLE :: Crt(:)
    INTEGER :: Tstep = -1
    LOGICAL :: Parallel = .FALSE.
    CHARACTER(:), ALLOCATABLE :: MultName
  END TYPE CircuitModel_t
!-------------------Circuit stuff----------------------------------------------

!------------------------------------------------------------------------------
    TYPE Model_t
!------------------------------------------------------------------------------
!
!     Coordinate system dimension + type
!
      INTEGER :: DIMENSION, CoordinateSystem
!
!     Model dimensions
!
      INTEGER :: NumberOfBulkElements=0, &
                 NumberOfNodes=0,        &
                 NumberOfBoundaryElements=0
!
!     Simulation input data, that concern the model as a whole
!
      TYPE(ValueList_t), POINTER :: Simulation => NULL()
!
!     Variables
!
      TYPE(Variable_t), POINTER  :: Variables => NULL()

!     External control of the simulation to sweep over parameter space.
      TYPE(ValueList_t), POINTER :: Control => NULL()
      
!     Some physical constants, that will be read from the database or set by
!     other means: gravity direction/intensity and Stefan-Boltzmann constant)
!
      TYPE(ValueList_t), POINTER :: Constants => NULL()
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

      TYPE(Mesh_t), POINTER :: Meshes => NULL()

      TYPE(Mesh_t),   POINTER :: Mesh   => NULL()
      TYPE(Solver_t), POINTER :: Solver => NULL()
      
      ! Circuits: one container per circuit solver, plus a pointer to the one
      ! whose equations are currently being handled. Everything in the circuit
      ! package reads its state through CircuitModel rather than from here.
      TYPE(CircuitModel_t), POINTER :: CircuitModels(:) => NULL()
      TYPE(CircuitModel_t), POINTER :: CircuitModel => NULL()

! Tag counts to speed things up
      INTEGER :: NumberOfDistTags=-1,NumberOfParTags=-1
      
    END TYPE Model_t

    TYPE(Model_t),  POINTER :: CurrentModel => NULL()
    TYPE(Matrix_t), POINTER :: GlobalMatrix => NULL()

    INTEGER :: ELMER_COMM_WORLD = -1

    CHARACTER(len=MAX_NAME_LEN) :: ExecID
!------------------------------------------------------------------------------
END MODULE Types
!------------------------------------------------------------------------------
!> \}
