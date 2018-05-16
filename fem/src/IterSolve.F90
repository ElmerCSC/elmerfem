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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! ****************************************************************************/

#include "huti_fdefs.h"

!> \ingroup ElmerLib 
!> \{




!------------------------------------------------------------------------------
!> Module containing a the control for iterative solvers for linear systems
!> that come with the Elmer suite.
!------------------------------------------------------------------------------
MODULE IterSolve

   USE Lists
   USE CRSMatrix
   USE BandMatrix
   USE IterativeMethods

#ifdef USE_ISO_C_BINDINGS
   USE huti_sfe
#endif

   IMPLICIT NONE

   !/*
   ! * Iterative method selection
   ! */
   INTEGER, PARAMETER, PRIVATE :: ITER_BiCGStab     =           320
   INTEGER, PARAMETER, PRIVATE :: ITER_TFQMR        =           330
   INTEGER, PARAMETER, PRIVATE :: ITER_CG           =           340
   INTEGER, PARAMETER, PRIVATE :: ITER_CGS          =           350
   INTEGER, PARAMETER, PRIVATE :: ITER_GMRES        =           360
   INTEGER, PARAMETER, PRIVATE :: ITER_BiCGStab2    =           370
   INTEGER, PARAMETER, PRIVATE :: ITER_SGS          =           380
   INTEGER, PARAMETER, PRIVATE :: ITER_JACOBI       =           390
   INTEGER, PARAMETER, PRIVATE :: ITER_RICHARDSON   =           391
   INTEGER, PARAMETER, PRIVATE :: ITER_BICGSTABL    =           400
   INTEGER, PARAMETER, PRIVATE :: ITER_GCR          =           410
   INTEGER, PARAMETER, PRIVATE :: ITER_IDRS         =           420

   !/*
   ! * Preconditioning type code
   ! */
   INTEGER, PARAMETER, PRIVATE :: PRECOND_NONE      =           500
   INTEGER, PARAMETER, PRIVATE :: PRECOND_DIAGONAL  =           510
   INTEGER, PARAMETER, PRIVATE :: PRECOND_ILUn      =           520
   INTEGER, PARAMETER, PRIVATE :: PRECOND_ILUT      =           530
   INTEGER, PARAMETER, PRIVATE :: PRECOND_MG        =           540
   INTEGER, PARAMETER, PRIVATE :: PRECOND_BILUn     =           550
   INTEGER, PARAMETER, PRIVATE :: PRECOND_Vanka     =           560
   INTEGER, PARAMETER, PRIVATE :: PRECOND_Circuit   =           570

   INTEGER, PARAMETER :: stack_max=64
   INTEGER :: stack_pos=0
   LOGICAL :: FirstCall(stack_max)

CONTAINS


#ifndef HUTI_MAXTOLERANCE
#define HUTI_MAXTOLERANCE dpar(2)
#endif
#ifndef HUTI_SGSPARAM
#define HUTI_SGSPARAM dpar(3)
#endif
#ifndef HUTI_BICGSTABL_L
#define HUTI_BICGSTABL_L ipar(16)
#endif
#ifndef HUTI_DIVERGENCE
#define HUTI_DIVERGENCE 3
#endif
#ifndef HUTI_GCR_RESTART
#define HUTI_GCR_RESTART ipar(17)
#endif
#ifndef HUTI_IDRS_S
#define HUTI_IDRS_S ipar(18)
#endif

!------------------------------------------------------------------------------
!> Dummy preconditioner, if linear system scaling is active this corresponds
!> to diagonal preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE pcond_dummy(u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM)
    INTEGER :: i
!------------------------------------------------------------------------------
    !$OMP PARALLEL DO
    DO i=1,HUTI_NDIM
       u(i) = v(i)
    END DO
    !$OMP END PARALLEL DO
!------------------------------------------------------------------------------
  END SUBROUTINE pcond_dummy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Complex dummy preconditioner, if linear system scaling is active this corresponds
!> to diagonal preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE pcond_dummy_cmplx(u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM)
!------------------------------------------------------------------------------
    u = v 
!------------------------------------------------------------------------------
  END SUBROUTINE pcond_dummy_cmplx
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The routine that decides which linear system solver to call, and calls it.
!> There are two main sources of iterations within Elmer.
!> 1) The old HUTiter C++ library that includes the most classic iterative Krylov
!>    methods.
!> 2) The internal MODULE IterativeMethods that includes some classic iterative
!>    methods and also some more recent Krylov methods. 
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE IterSolver( A,x,b,Solver,ndim,DotF, &
              NormF,MatvecF,PrecF,StopcF )
!------------------------------------------------------------------------------
#ifdef USE_ISO_C_BINDINGS
    USE huti_sfe
#endif
    USE ListMatrix
    USE SParIterGlobals
    IMPLICIT NONE

!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
    TYPE(Matrix_t), TARGET :: A
    INTEGER, OPTIONAL :: ndim
    INTEGER(KIND=AddrInt), OPTIONAL :: DotF, NormF, MatVecF, PrecF, StopcF
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: Adiag,CM,PrecMat,SaveGlobalM

    REAL(KIND=dp) :: dpar(HUTI_DPAR_DFLTSIZE),stopfun
!   external stopfun
    REAL(KIND=dp), ALLOCATABLE :: work(:,:)
    INTEGER :: i,j,k,N,ipar(HUTI_IPAR_DFLTSIZE),wsize,istat,IterType,PCondType,ILUn,Blocks
    LOGICAL :: Internal, NullEdges
    LOGICAL :: ComponentwiseStopC, NormwiseStopC, RowEquilibration
    LOGICAL :: Condition,GotIt, Refactorize,Found,GotDiagFactor,Robust
    LOGICAL :: ComplexSystem
    
    REAL(KIND=dp) :: ILUT_TOL, DiagFactor

    TYPE(ValueList_t), POINTER :: Params

    CHARACTER(LEN=MAX_NAME_LEN) :: str

    EXTERNAL MultigridPrec
    EXTERNAL NormwiseBackwardError, ComponentwiseBackwardError
    EXTERNAL NormwiseBackwardErrorGeneralized
#ifndef USE_ISO_C_BINDINGS
    INTEGER  :: HUTI_D_BICGSTAB, HUTI_D_BICGSTAB_2, HUTI_D_TFQMR, &
        HUTI_D_CG, HUTI_D_CGS, HUTI_D_GMRES
    EXTERNAL :: HUTI_D_BICGSTAB, HUTI_D_BICGSTAB_2, HUTI_D_TFQMR, &
        HUTI_D_CG, HUTI_D_CGS, HUTI_D_GMRES
    
    INTEGER  :: HUTI_Z_BICGSTAB, HUTI_Z_BICGSTAB_2, HUTI_Z_TFQMR, &
        HUTI_Z_CG, HUTI_Z_CGS, HUTI_Z_GMRES
    EXTERNAL :: HUTI_Z_BICGSTAB, HUTI_Z_BICGSTAB_2, HUTI_Z_TFQMR, &
        HUTI_Z_CG, HUTI_Z_CGS, HUTI_Z_GMRES

    REAL(KIND=dp) :: ddot, dnrm2, dznrm2
    EXTERNAL :: ddot, dnrm2, dznrm2
    
    COMPLEX(KIND=dp) :: zdotc
    EXTERNAL :: zdotc
#endif
    
    INTEGER(KIND=Addrint) :: dotProc, normProc, pcondProc, &
        pcondrProc, mvProc, iterProc, StopcProc
    INTEGER(KIND=Addrint) :: AddrFunc
#ifdef USE_ISO_C_BINDINGS
    INTEGER :: astat
    COMPLEX(KIND=dp), ALLOCATABLE :: xC(:), bC(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: workC(:,:)
    EXTERNAL :: AddrFunc    
#endif

    INTERFACE
      SUBROUTINE VankaCreate(A,Solver)
        USE Types
        TYPE(Matrix_t) :: A
        TYPE(Solver_t) :: Solver
      END SUBROUTINE VankaCreate
      
      SUBROUTINE VankaPrec(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        REAL(KIND=dp) :: u(*),v(*)
      END SUBROUTINE VankaPrec

      SUBROUTINE CircuitPrec(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        REAL(KIND=dp) :: u(*),v(*)
      END SUBROUTINE CircuitPrec

      SUBROUTINE CircuitPrecComplex(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        COMPLEX(KIND=dp) :: u(*),v(*)
      END SUBROUTINE CircuitPrecComplex
    END INTERFACE
!------------------------------------------------------------------------------
    N = A % NumberOfRows
    IF ( PRESENT(ndim) ) n=ndim
    
    ipar = 0
    dpar = 0.0D0
    pcondrProc = 0
!------------------------------------------------------------------------------
    Params => Solver % Values
    str = ListGetString( Params,'Linear System Iterative Method',Found )
    IF( .NOT. Found ) THEN
      CALL Warn('IterSolver','> Linear System Iterative Method < not found, using BiCGstab')
      str = 'bicgstab'      
    ELSE
      CALL Info('IterSolver','Using iterative method: '//TRIM(str),Level=9)
    END IF
    
    ComplexSystem = ListGetLogical( Params,'Linear System Complex',Found ) 
    IF( .NOT. Found ) ComplexSystem = A % COMPLEX 

    IF( ListGetLogical( Params,'Linear System Skip Complex',GotIt ) ) THEN
      CALL Info('IterSolver','This time skipping complex treatment',Level=20)
      A % COMPLEX = .FALSE.
      ComplexSystem = .FALSE.
    END IF
            
    IF( ComplexSystem ) THEN
      CALL Info('IterSolver','Matrix is complex valued',Level=10)
    ELSE
      CALL Info('IterSolver','Matrix is real valued',Level=12)
    END IF
   
    SELECT CASE(str)
    CASE('bicgstab2')
      IterType = ITER_BiCGStab2
    CASE('bicgstabl')
      IterType = ITER_BICGstabl
    CASE('bicgstab')
      IterType = ITER_BiCGStab
    CASE('tfqmr')
      IterType = ITER_TFQMR
    CASE('cgs')
      IterType = ITER_CGS
    CASE('cg')
      IterType = ITER_CG
    CASE('gmres')
      IterType = ITER_GMRES
    CASE('sgs')
      IterType = ITER_SGS
    CASE('jacobi')
      IterType = ITER_jacobi
    CASE('richardson')
      IterType = ITER_richardson
    CASE('gcr')
      IterType = ITER_GCR
    CASE('idrs')
      IterType = ITER_IDRS
    CASE DEFAULT
      IterType = ITER_BiCGStab
    END SELECT
    
!------------------------------------------------------------------------------

    HUTI_WRKDIM = 0
    Internal = .FALSE.
    
    SELECT CASE ( IterType )
      
      ! Solvers from HUTiter
      !-------------------------------------------------------       
    CASE (ITER_BiCGStab)
      HUTI_WRKDIM = HUTI_BICGSTAB_WORKSIZE
      
    CASE (ITER_BiCGStab2)
      HUTI_WRKDIM = HUTI_BICGSTAB_2_WORKSIZE
      
    CASE (ITER_TFQMR)
      HUTI_WRKDIM = HUTI_TFQMR_WORKSIZE
      
    CASE (ITER_CG)
      HUTI_WRKDIM = HUTI_CG_WORKSIZE
      
    CASE (ITER_CGS)
      HUTI_WRKDIM = HUTI_CGS_WORKSIZE
      
    CASE (ITER_GMRES)
      HUTI_GMRES_RESTART = ListGetInteger( Params,&
          'Linear System GMRES Restart',  GotIt ) 
      IF ( .NOT. GotIT ) HUTI_GMRES_RESTART = 10
      HUTI_WRKDIM = HUTI_GMRES_WORKSIZE + HUTI_GMRES_RESTART
      
      ! Solvers from IterativeMethods.src
      !-------------------------------------------------------       
    CASE (ITER_SGS)
      HUTI_WRKDIM = 1
      HUTI_SGSPARAM = ListGetConstReal( Params,'SGS Overrelaxation Factor',&
          GotIt,minv=0.0_dp,maxv=2.0_dp)
      IF(.NOT. GotIt) HUTI_SGSPARAM = 1.8
      Internal = .TRUE.
      
    CASE (ITER_Jacobi, ITER_Richardson)
      HUTI_WRKDIM = 1
      Internal = .TRUE.
      
    CASE (ITER_GCR)
      HUTI_WRKDIM = 1
      HUTI_GCR_RESTART = ListGetInteger( Params, &
          'Linear System GCR Restart',  GotIt ) 
      IF ( .NOT. GotIT ) HUTI_GCR_RESTART = ListGetInteger( Params, &
          'Linear System Max Iterations', minv=1 )
      Internal = .TRUE.
      
    CASE (ITER_BICGSTABL)
      HUTI_WRKDIM = 1
      HUTI_BICGSTABL_L = ListGetInteger( Params,'BiCGstabl polynomial degree',&
          GotIt,minv=2)
      IF(.NOT. GotIt) HUTI_BICGSTABL_L = 2
      Internal = .TRUE.

    CASE (ITER_IDRS)
      HUTI_WRKDIM = 1
      HUTI_IDRS_S = ListGetInteger( Params,'IDRS parameter',GotIt,minv=1)
      IF(.NOT. GotIt) HUTI_IDRS_S = 4
      Internal = .TRUE.
      
    END SELECT
!------------------------------------------------------------------------------
    
    wsize = HUTI_WRKDIM
    
    StopcProc = 0
    IF (PRESENT(StopcF)) THEN
       StopcProc = StopcF
       HUTI_STOPC = HUTI_USUPPLIED_STOPC
    ELSE
       ComponentwiseStopC = ListGetLogical(Params,'Linear System Componentwise Backward Error',GotIt)
       IF (ComponentwiseStopC) THEN
          StopcProc = AddrFunc(ComponentwiseBackwardError)
          HUTI_STOPC = HUTI_USUPPLIED_STOPC
       ELSE
          NormwiseStopC = ListGetLogical(Params,'Linear System Normwise Backward Error',GotIt)
          IF (NormwiseStopC) THEN
             RowEquilibration = ListGetLogical(Params,'Linear System Row Equilibration',GotIt)
             IF (RowEquilibration) THEN
                StopcProc = AddrFunc(NormwiseBackwardError)              
             ELSE
                StopcProc = AddrFunc(NormwiseBackwardErrorGeneralized)
             END IF
             HUTI_STOPC = HUTI_USUPPLIED_STOPC
          ELSE
             HUTI_STOPC = HUTI_TRESID_SCALED_BYB
          END IF
       END IF
    END IF
    HUTI_NDIM  = N
    
    HUTI_DBUGLVL  = ListGetInteger( Params, &
        'Linear System Residual Output', GotIt )
    IF ( .NOT.Gotit ) HUTI_DBUGLVL = 1
    
    IF ( Parenv % myPE /= 0 ) HUTI_DBUGLVL=0
    
    HUTI_MAXIT = ListGetInteger( Params, &
        'Linear System Max Iterations', minv=1 )
    
    HUTI_MINIT = ListGetInteger( Params, &
        'Linear System Min Iterations', GotIt )
    
#ifdef USE_ISO_C_BINDINGS
    IF( ComplexSystem ) THEN
        ALLOCATE(workC(N/2,wsize), stat=istat)
        IF ( istat /= 0 ) THEN
            CALL Fatal( 'IterSolve', 'Memory allocation failure.' )
        END IF
        workC = cmplx(0,0,dp)
    ELSE
        ALLOCATE(work(N,wsize), stat=istat)
        IF ( istat /= 0 ) THEN
            CALL Fatal( 'IterSolve', 'Memory allocation failure.' )
        END IF
        !$OMP PARALLEL PRIVATE(j)
        DO j=1,wsize
           !$OMP DO
           DO i=1,N
              work(i,j) = real(0,dp)
           END DO
           !$OMP END DO
        END DO
        !$OMP END PARALLEL
    END IF
#else
    ALLOCATE( work(N,wsize),stat=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'IterSolve', 'Memory allocation failure.' )
    END IF
    work=0._dp
#endif

    IF ( ALL(x == 0.0) ) x = 1.0d-8
    HUTI_INITIALX = HUTI_USERSUPPLIEDX
    
    HUTI_TOLERANCE = ListGetCReal( Params, &
        'Linear System Convergence Tolerance' )
    
    HUTI_MAXTOLERANCE = ListGetCReal( Params, &
        'Linear System Divergence Limit', GotIt)
    IF(.NOT. GotIt) HUTI_MAXTOLERANCE = 1.0d20
    
    IF( ListGetLogical( Params,'Linear System Robust',GotIt) ) THEN
      HUTI_ROBUST = 1
      HUTI_ROBUST_TOLERANCE = ListGetCReal( Params,'Linear System Robust Tolerance',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_TOLERANCE = HUTI_TOLERANCE**(2.0/3.0)
      HUTI_ROBUST_MAXTOLERANCE = ListGetCReal( Params,'Linear System Robust Limit',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_MAXTOLERANCE = SQRT( HUTI_TOLERANCE )      
      HUTI_ROBUST_STEPSIZE = ListGetCReal( Params,'Linear System Robust Margin',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_STEPSIZE = 1.1_dp
      HUTI_ROBUST_MAXBADIT = ListGetInteger( Params,'Linear System Robust Max Iterations',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_MAXBADIT = HUTI_MAXIT / 2
    ELSE
      HUTI_ROBUST = 0
    END IF


    IF( ListGetLogical( Params,'IDRS Smoothing',GotIt) ) THEN
      HUTI_SMOOTHING = 1
    ELSE
      HUTI_SMOOTHING = 0
    END IF
      
    
!------------------------------------------------------------------------------

    
    IF ( .NOT. PRESENT(PrecF) ) THEN
      str = ListGetString( Params, &
          'Linear System Preconditioning',gotit )
      IF ( .NOT.gotit ) str = 'none'
      
      A % Cholesky = ListGetLogical( Params, &
          'Linear System Symmetric ILU', Gotit )
      
      ILUn = -1
      IF ( str == 'none' ) THEN
        PCondType = PRECOND_NONE

      ELSE IF ( str == 'diagonal' ) THEN
        PCondType = PRECOND_DIAGONAL

      ELSE IF ( str == 'ilut' ) THEN
        ILUT_TOL = ListGetCReal( Params, &
            'Linear System ILUT Tolerance',GotIt )
        PCondType = PRECOND_ILUT

      ELSE IF ( SEQL(str, 'ilu') ) THEN
        ILUn = NINT(ListGetCReal( Params, &
            'Linear System ILU Order', gotit ))
        IF ( .NOT.gotit ) &
            ILUn = ICHAR(str(4:4)) - ICHAR('0')
        IF ( ILUn  < 0 .OR. ILUn > 9 ) ILUn = 0
        PCondType = PRECOND_ILUn

      ELSE IF ( SEQL(str, 'bilu') ) THEN
        ILUn = ICHAR(str(5:5)) - ICHAR('0')
        IF ( ILUn  < 0 .OR. ILUn > 9 ) ILUn = 0
        IF( Solver % Variable % Dofs == 1) THEN
          CALL Warn('IterSolver','BILU for one dofs is equal to ILU!')
          PCondType = PRECOND_ILUn
        ELSE
          PCondType = PRECOND_BILUn
        END IF

      ELSE IF ( str == 'multigrid' ) THEN
        PCondType = PRECOND_MG

      ELSE IF ( str == 'vanka' ) THEN
        PCondType = PRECOND_VANKA

      ELSE IF ( str == 'circuit' ) THEN
        ILUn = ListGetInteger( Params, 'Linear System ILU Order', gotit )
        IF(.NOT.Gotit ) ILUn=-1
        PCondType = PRECOND_Circuit

      ELSE
        PCondType = PRECOND_NONE
        CALL Warn( 'IterSolve', 'Unknown preconditioner type, feature disabled.' )
      END IF
      
      IF ( .NOT. ListGetLogical( Params, 'No Precondition Recompute',GotIt ) ) THEN
        n = ListGetInteger( Params, 'Linear System Precondition Recompute', GotIt )
        IF ( n <= 0 ) n = 1
        
        Refactorize = ListGetLogical( Params, 'Linear System Refactorize', Gotit )
        IF ( .NOT. Gotit ) Refactorize = .TRUE.
        
        IF (.NOT.(ASSOCIATED(A % ILUValues).OR.ASSOCIATED(A % CILUValues)).OR. &
                  (Refactorize.AND.MOD(A % SolveCount, n)==0) ) THEN


          IF ( A % FORMAT == MATRIX_CRS ) THEN

            ! Optionally one may emphasize the diagonal entries in the linear system
            ! to make the preconditioning more stable.
            !-------------------------------------------------------------------------          
            DiagFactor = ListGetCReal( Params,'Linear System ILU Factor',GotIt ) 
            GotDiagFactor = ( DiagFactor > EPSILON( DiagFactor ) ) 
            IF( GotDiagFactor ) THEN
              CALL Info('IterSolver','Applying diagonal relaxation for ILU', Level=8)
              DiagFactor = DiagFactor + 1.0_dp
              A % Values( A % Diag ) = DiagFactor * A % Values( A % Diag )      
            END IF

            IF ( ComplexSystem ) THEN
              IF ( PCondType == PRECOND_ILUn .OR. (PCondType==PRECOND_Circuit.AND.ILUn>=0) ) THEN
                NullEdges = ListGetLogical(Params, 'Edge Basis', GotIt)
                CM => A % ConstraintMatrix
                IF(NullEdges.OR.ASSOCIATED(CM)) THEN

                  IF(ASSOCIATED(A % ILURows)) DEALLOCATE(A % ILURows)
                  IF(ASSOCIATED(A % ILUCols)) DEALLOCATE(A % ILUCols)
                  IF(ASSOCIATED(A % ILUDiag)) DEALLOCATE(A % ILUDiag)
                  IF(ASSOCIATED(A % CILUValues)) DEALLOCATE(A % CILUValues)

                  PrecMat => AllocateMatrix()
                  PrecMat % FORMAT = MATRIX_LIST
                  PrecMat % CIluValues => NULL()

                  IF(ASSOCIATED(CM)) THEN
                    DO i=CM % NumberOfRows,1,-1
                      k = i + A % NumberOfRows
                      CALL List_AddMatrixIndex( PrecMat % ListMatrix,k,k)
                      IF(MOD(k,2)==0) THEN
                        CALL List_AddMatrixIndex(PrecMat % ListMatrix, k, k-1)
                      ELSE
                        CALL List_AddMatrixIndex(PrecMat % ListMatrix, k, k+1)
                      END IF

                      DO j=CM % Rows(i+1)-1,CM % Rows(i),-1
                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, CM % Cols(j), CM % Values(j))

                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             CM % Cols(j), i + A % NumberOfRows, CM % Values(j))
                      END DO
                    END DO
                  END IF

                  k = A % NumberOfRows - A % ExtraDOFs
                  DO i=A % NumberOfRows,1,-1
                    IF(i>k) THEN
                       CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i)
                       IF(MOD(i,2)==0) THEN
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i-1)
                       ELSE
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i+1)
                       END IF
                    ELSE IF (NullEdges) THEN
                       CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, i, 1._dp)
                       IF(MOD(i,2)==0) THEN
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i-1)
                       ELSE
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i+1)
                       END IF
                    END IF

                    DO j=A % Rows(i+1)-1,A % Rows(i),-1
                      IF (i>k .OR. A % Cols(j)>k .OR. .NOT.NullEdges) THEN
                        CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, A % Cols(j), A % Values(j))
                      END IF
                    END DO
                  END DO

                  CALL List_ToCRSMatrix(PrecMat)
                  Condition = CRS_ComplexIncompleteLU(PrecMat,ILUn)

                  A % ILURows => PrecMat % IluRows
                  A % ILUCols => PrecMat % IluCols
                  A % ILUDiag => PrecMat % IluDiag
                  A % CILUvalues => PrecMat % CIluValues

                  DEALLOCATE(PrecMat % Values)
                  IF(.NOT.ASSOCIATED(A % ILURows,PrecMat % Rows)) DEALLOCATE(PrecMat % Rows)
                  IF(.NOT.ASSOCIATED(A % ILUCols,PrecMat % Cols)) DEALLOCATE(PrecMat % Cols)
                  IF(.NOT.ASSOCIATED(A % ILUDiag,PrecMat % Diag)) DEALLOCATE(PrecMat % Diag)
                  DEALLOCATE(PrecMat)
                ELSE
                  Condition = CRS_ComplexIncompleteLU(A,ILUn)
                END IF
              ELSE IF ( PCondType == PRECOND_ILUT ) THEN
                Condition = CRS_ComplexILUT( A,ILUT_TOL )
              END IF
            ELSE IF (ILUn>=0) THEN
              SELECT CASE(PCondType)
              CASE(PRECOND_ILUn, PRECOND_Circuit)
                NullEdges = ListGetLogical(Params, 'Edge Basis', GotIt)
                CM => A % ConstraintMatrix
                IF(NullEdges.OR.ASSOCIATED(CM)) THEN

                  IF(ASSOCIATED(A % ILURows)) DEALLOCATE(A % ILURows)
                  IF(ASSOCIATED(A % ILUCols)) DEALLOCATE(A % ILUCols)
                  IF(ASSOCIATED(A % ILUDiag)) DEALLOCATE(A % ILUDiag)
                  IF(ASSOCIATED(A % ILUValues)) DEALLOCATE(A % ILUValues)

                  PrecMat => AllocateMatrix()
                  PrecMat % FORMAT = MATRIX_LIST

                  IF(ASSOCIATED(CM)) THEN
                    DO i=CM % NumberOfRows,1,-1
                      CALL List_AddMatrixIndex( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, i + A % NumberOFrows)

                      DO j=CM % Rows(i+1)-1,CM % Rows(i),-1
                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, CM % Cols(j), CM % Values(j))

                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             CM % Cols(j), i + A % NumberOfRows, CM % Values(j))
                      END DO
                    END DO
                  END IF

                  k = A % NumberOfRows - A % ExtraDOFs
                  DO i=A % NumberOfRows,1,-1
                    IF(i>k) THEN
                       CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i)
                    ELSE IF (NullEdges) THEN
                       CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, i, 1._dp)
                    END IF

                    DO j=A % Rows(i+1)-1,A % Rows(i),-1
                      IF (i>k .OR. A % Cols(j)>k .OR. .NOT.NullEdges) THEN
                        CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, A % Cols(j), A % Values(j))
                      END IF
                    END DO
                  END DO

                  CALL List_ToCRSMatrix(PrecMat)
                  Condition = CRS_IncompleteLU(PrecMat,ILUn)

                  A % ILURows => PrecMat % IluRows
                  A % ILUCols => PrecMat % IluCols
                  A % ILUDiag => PrecMat % IluDiag
                  A % ILUvalues => PrecMat % IluValues

                  DEALLOCATE(PrecMat % Values)
                  IF(.NOT.ASSOCIATED(A % ILURows,PrecMat % Rows)) DEALLOCATE(PrecMat % Rows)
                  IF(.NOT.ASSOCIATED(A % ILUCols,PrecMat % Cols)) DEALLOCATE(PrecMat % Cols)
                  IF(.NOT.ASSOCIATED(A % ILUDiag,PrecMat % Diag)) DEALLOCATE(PrecMat % Diag)
                  DEALLOCATE(PrecMat)
                ELSE
                  Condition = CRS_IncompleteLU(A,ILUn)
                END IF
              CASE(PRECOND_ILUT)
                Condition = CRS_ILUT( A,ILUT_TOL )
              CASE(PRECOND_BILUn)
                Blocks = Solver % Variable % Dofs
                IF ( Blocks <= 1 ) THEN
                  Condition = CRS_IncompleteLU(A,ILUn)
                ELSE
                  IF( .NOT. ASSOCIATED( A % ILUValues ) ) THEN
                    Adiag => AllocateMatrix()
                    CALL CRS_BlockDiagonal(A,Adiag,Blocks)
                    Condition = CRS_IncompleteLU(Adiag,ILUn)
                    A % ILURows   => Adiag % ILURows
                    A % ILUCols   => Adiag % ILUCols
                    A % ILUValues => Adiag % ILUValues
                    A % ILUDiag   => Adiag % ILUDiag                 
                    IF (ILUn > 0) THEN
                      DEALLOCATE(Adiag % Rows,Adiag % Cols, Adiag % Diag, Adiag % Values)
                    END IF
                    DEALLOCATE( Adiag )
                  ELSE
                    Condition = CRS_IncompleteLU(A,ILUn)
                  END IF
                END IF
              CASE(PRECOND_VANKA)
                !                  CALL VankaCreate( A, SolverParam )
              END SELECT
            END IF

            IF( GotDiagFactor ) THEN
              CALL Info('IterSolver','Reverting diagonal relaxation for ILU', Level=10)
               A % Values( A % Diag ) = A % Values( A % Diag ) / DiagFactor
            END IF

          ELSE
            IF ( PCondType == PRECOND_ILUn ) THEN
              CALL Warn( 'IterSolve', 'No ILU Preconditioner for Band Matrix format,' )
              CALL Warn( 'IterSolve', 'using Diagonal preconditioner instead...' )
              PCondType = PRECOND_DIAGONAL
            END IF
          END IF
        END IF
      END IF
    END IF
    
    A % SolveCount = A % SolveCount + 1
!------------------------------------------------------------------------------

    IF ( PRESENT(MatvecF) ) THEN
      mvProc = MatvecF
    ELSE
      IF ( .NOT. ComplexSystem ) THEN
        mvProc = AddrFunc( CRS_MatrixVectorProd )
      ELSE
        mvProc = AddrFunc( CRS_ComplexMatrixVectorProd )
      END IF
    END IF
    
    IF ( PRESENT(dotF) ) THEN
      dotProc = dotF
    ELSE
      dotProc = 0
    END IF
    
    IF ( PRESENT(normF) ) THEN
      normProc = normF
    ELSE
      normProc = 0
    END IF
    
    
    IF ( PRESENT(PrecF) ) THEN
      pcondProc = PrecF
    ELSE
      SELECT CASE( PCondType )
      CASE (PRECOND_NONE)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( pcond_dummy )
        ELSE
          pcondProc = AddrFunc( pcond_dummy_cmplx  )
        END IF
        
      CASE (PRECOND_DIAGONAL)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CRS_DiagPrecondition )
        ELSE
          pcondProc = AddrFunc( CRS_ComplexDiagPrecondition )
        END IF
        
      CASE (PRECOND_ILUn, PRECOND_ILUT, PRECOND_BILUn )
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CRS_LUPrecondition )
        ELSE
          pcondProc = AddrFunc( CRS_ComplexLUPrecondition )
        END IF

      CASE (PRECOND_MG)
        pcondProc = AddrFunc( MultiGridPrec )
        
      CASE (PRECOND_VANKA)
        pcondProc = AddrFunc( VankaPrec )

      CASE (PRECOND_Circuit)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CircuitPrec )
        ELSE
          pcondProc = AddrFunc( CircuitPrecComplex )
        END IF
        
      CASE DEFAULT
        pcondProc = 0
      END SELECT
    END IF
    

    IF ( .NOT. ComplexSystem ) THEN
      SELECT CASE ( IterType )

       ! Solvers from HUTiter library 
       !-------------------------------------------------------       
      CASE (ITER_BiCGStab)
        iterProc = AddrFunc( HUTI_D_BICGSTAB )
      CASE (ITER_BiCGStab2)
        iterProc = AddrFunc( HUTI_D_BICGSTAB_2 )
      CASE (ITER_TFQMR)
        iterProc = AddrFunc( HUTI_D_TFQMR )
      CASE (ITER_CG)
        iterProc = AddrFunc( HUTI_D_CG )
      CASE (ITER_CGS)
        iterProc = AddrFunc( HUTI_D_CGS )
      CASE (ITER_GMRES)
        iterProc = AddrFunc( HUTI_D_GMRES )
        
        ! Solvers from IterativeMethods.src 
        !-------------------------------------------------------
      CASE (ITER_SGS)
        iterProc = AddrFunc( itermethod_sgs )
      CASE (ITER_JACOBI)
        iterProc = AddrFunc( itermethod_jacobi )
      CASE (ITER_RICHARDSON)
        iterProc = AddrFunc( itermethod_richardson )
      CASE (ITER_GCR)
        iterProc = AddrFunc( itermethod_gcr )
      CASE (ITER_BICGSTABL)
        iterProc = AddrFunc( itermethod_bicgstabl )
      CASE (ITER_IDRS)
        iterProc = AddrFunc( itermethod_idrs )
        
      END SELECT
      
      IF( Internal ) THEN
        IF ( dotProc  == 0 ) dotProc = AddrFunc(ddot)
        IF ( normProc == 0 ) normproc = AddrFunc(dnrm2)
        IF( HUTI_DBUGLVL == 0) HUTI_DBUGLVL = HUGE( HUTI_DBUGLVL )        
      END IF
      
    ELSE
      HUTI_NDIM = HUTI_NDIM / 2
      SELECT CASE ( IterType )

        ! Solvers from HUTiter library 
        !-------------------------------------------------------       
      CASE (ITER_BiCGStab)
        iterProc = AddrFunc( HUTI_Z_BICGSTAB )
      CASE (ITER_BiCGStab2)
        iterProc = AddrFunc( HUTI_Z_BICGSTAB_2 )
      CASE (ITER_TFQMR)
        iterProc = AddrFunc( HUTI_Z_TFQMR )
      CASE (ITER_CG)
        iterProc = AddrFunc( HUTI_Z_CG )
      CASE (ITER_CGS)
        iterProc = AddrFunc( HUTI_Z_CGS )
      CASE (ITER_GMRES)
        iterProc = AddrFunc( HUTI_Z_GMRES )
        
        ! Solvers from IterativeMethods.src 
        !-------------------------------------------------------
      CASE (ITER_GCR)
        iterProc = AddrFunc( itermethod_z_gcr )
      CASE (ITER_BICGSTABL)
        iterProc = AddrFunc( itermethod_z_bicgstabl )
      CASE (ITER_IDRS)
        iterProc = AddrFunc( itermethod_z_idrs )

      END SELECT
      
      IF( Internal ) THEN
        IF ( dotProc  == 0 ) dotProc = AddrFunc(zdotc)
        IF ( normProc == 0 ) normproc = AddrFunc(dznrm2)
        IF( HUTI_DBUGLVL == 0) HUTI_DBUGLVL = HUGE( HUTI_DBUGLVL )
      END IF
      
    END IF
    
!------------------------------------------------------------------------------

    stack_pos = stack_pos+1
    IF(stack_pos>stack_max) THEN
      CALL Fatal('IterSolver', 'Recursion too deep ('//TRIM(I2S(stack_pos))//' vs '//TRIM(I2S(stack_max))//')')
    ELSE IF(stack_pos<=0) THEN
      CALL Fatal('IterSolver', 'eh')
    END IF
    FirstCall(stack_pos) = .TRUE.

    SaveGlobalM => GlobalMatrix
    GlobalMatrix => A
    
#ifdef USE_ISO_C_BINDINGS
    IF ( ComplexSystem ) THEN
      ! Associate xC and bC with complex variables
      ALLOCATE(xC(HUTI_NDIM), bC(HUTI_NDIM), STAT=astat)
      IF (astat /= 0) THEN
        CALL Fatal('IterSolve','Unable to allocate memory for complex arrays')
      END IF
      ! Initialize xC and bC
      DO i=1,HUTI_NDIM
        xC(i) = cmplx(x(2*i-1),x(2*i),dp)
      END DO
      DO i=1,HUTI_NDIM
        bC(i) = cmplx(b(2*i-1),b(2*i),dp)
      END DO

      CALL Info('IterSolver','Calling complex iterative solver',Level=32)
      CALL IterCall( iterProc, xC, bC, ipar, dpar, workC, &
          mvProc, pcondProc, pcondrProc, dotProc, normProc, stopcProc )

      ! Copy result back
      DO i=1,HUTI_NDIM
        x(2*i-1) = REAL(REAL(xC(i)),dp)
        x(2*i) = REAL(AIMAG(xC(i)),dp)
      END DO
      DEALLOCATE(bC,xC)
    ELSE
      CALL Info('IterSolver','Calling real valued iterative solver',Level=32)
      CALL IterCall( iterProc, x, b, ipar, dpar, work, &
          mvProc, pcondProc, pcondrProc, dotProc, normProc, stopcProc )
    ENDIF
#else
    CALL Info('IterSolver','Calling iterative solver',Level=32)   
    CALL IterCall( iterProc, x, b, ipar, dpar, work, &
        mvProc, pcondProc, pcondrProc, dotProc, normProc, stopcProc )
#endif
    GlobalMatrix => SaveGlobalM

    
    stack_pos=stack_pos-1
    
    IF ( ComplexSystem ) HUTI_NDIM = HUTI_NDIM * 2
!------------------------------------------------------------------------------
    IF ( HUTI_INFO /= HUTI_CONVERGENCE .AND. ParEnv % myPE==0 ) THEN
      CALL Info('IterSolve','Returned return code: '//TRIM(I2S(HUTI_INFO)),Level=15)
      IF( HUTI_INFO == HUTI_DIVERGENCE ) THEN
        CALL NumericalError( 'IterSolve', 'System diverged over maximum tolerance.')
      ELSE IF( HUTI_INFO == HUTI_MAXITER ) THEN
        CALL NumericalError( 'IterSolve', 'Too many iterations was needed.')
      END IF
    END IF
!------------------------------------------------------------------------------
#ifdef USE_ISO_C_BINDINGS
    IF ( ComplexSystem ) THEN
        DEALLOCATE( workC )
    ELSE
        DEALLOCATE( work )
    END IF
#else 
    DEALLOCATE( work )
#endif

!------------------------------------------------------------------------------
  END SUBROUTINE IterSolver
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> This routine may be used to either inform user or terminate following
!> convergence/numerical issues, based on a flag in the SIF. Default
!> behaviour terminates execution.
!-----------------------------------------------------------------------
   SUBROUTINE NumericalError( Caller, String, Fatal )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     LOGICAL, OPTIONAL :: Fatal
!-----------------------------------------------------------------------
     LOGICAL :: GlobalNumFatal, SolverNumFatal, IsFatal, Found
!-----------------------------------------------------------------------

     !Fatality logic:
     ! 1) Respect calling routine's wishes if present
     ! 2) Respect solver specific option if present
     ! 3) Respect global abort flag if present
     ! 4) Otherwise fatal (backwards compatibility)

     IF(PRESENT(Fatal)) THEN
       IsFatal = Fatal
     ELSE
       SolverNumFatal = ListGetLogical( CurrentModel % Solver % Values, &
            'Linear System Abort Not Converged', Found)
       IF(Found) THEN
         IsFatal = SolverNumFatal
       ELSE
         GlobalNumFatal = ListGetLogical(CurrentModel % Simulation,&
            'Global Abort Not Converged',Found)
         IF(Found) THEN
           IsFatal = GlobalNumFatal
         ELSE
           IsFatal = .TRUE.
         END IF
       END IF
     END IF

     IF ( OutputLevelMask(0) ) THEN
       IF(IsFatal) THEN
         WRITE( *, '(A,A,A,A)', ADVANCE='YES' ) &
              'NUMERICAL ERROR:: ', TRIM(Caller), ': ', TRIM(String)
       ELSE
         WRITE( *, '(A,A,A,A)', ADVANCE='YES' ) &
              'NUMERICAL WARNING:: ', TRIM(Caller), ': ', TRIM(String)
       END IF
       CALL FLUSH(6)
     END IF

     IF(IsFatal) STOP

!-----------------------------------------------------------------------
   END SUBROUTINE NumericalError
!-----------------------------------------------------------------------

END MODULE IterSolve

!> \} ElmerLib
