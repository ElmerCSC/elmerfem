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
! *  This routine is modified from the arpack examples driver dndrv3 to
! *  suit the needs of ELMER.
! *
! ******************************************************************************
! *
! * \Original Authors
! *     Richard Lehoucq
! *     Danny Sorensen
! *     Chao Yang
! *    Dept. of Computational &
! *    Applied Mathematics
! *    Rice University
! *    Houston, Texas
! *
! * FILE: ndrv3.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
! *
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
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!----------------------------------------------------------------------------
!> Parallel solver for eigenvalue problems.
!----------------------------------------------------------------------------


MODULE ParallelEigenSolve

   USE Messages
   USE Lists, ONLY : ListAddLogical, ListGetConstReal, ListGetInteger, &
                     ListGetLogical, ListGetString, ListRemove
   USE GeneralUtils, ONLY : seql, sortc, I2S
   USE SParIterGlobals
   USE CRSMatrix, ONLY : CRS_ComplexILUT, CRS_ComplexIncompleteLU, CRS_ILUT, &
                         CRS_IncompleteLU, CRS_LUSolve
   USE DirectSolve, ONLY : DirectSolver
   USE Multigrid, ONLY : MultiGridSolve
   USE EigenSolve, ONLY : ArpackSetWhich, EigenSystemSorting
   USE MatrixScaling, ONLY : ScaleLinearSystem, BackScaleLinearSystem
   USE ParallelUtils, ONLY : ParallelVector, ParallelInitSolve, ParallelIter, &
                             ParallelReduction, ParallelMatrixVector, ParallelNorm, &
                             ParallelDot, ParallelMatrix, PartitionVector, ParallelUpdateResult

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
     SUBROUTINE ParallelArpackEigenSolve( Solver,A,N,NEIG,EigValues,EigVectors )
!------------------------------------------------------------------------------

! to suit the needs of ELMER.
!
!  Oct 21 2000, Juha Ruokolainen
!
!\Original Authors
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
! FILE: ndrv3.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
!
!
!

#define PARALLEL_FOR_REAL
#ifdef PARALLEL_FOR_REAL
#  if defined(ELMER_HAVE_MPI_MODULE)
      USE mpi
#  endif
#endif

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: A
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef PARALLEL_FOR_REAL
#  if defined(ELMER_HAVE_MPIF_HEADER)
      INCLUDE "mpif.h"
#  endif
#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      TYPE(Matrix_t), POINTER :: PMatrix
      REAL(KIND=dp), TARGET :: WORKD(3*N), RESID(N)
      REAL(KIND=dp), POINTER CONTIG :: x(:), b(:)
      INTEGER :: IPARAM(11), IPNTR(14)
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      REAL(KIND=dp), ALLOCATABLE :: WORKL(:), D(:,:), WORKEV(:), V(:,:)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      CHARACTER ::     BMAT*1, Which*2, DirectMethod*32
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, k, l, p, pn, IERR, iter, &
                       NCONV, maxitr, ishfts, mode, istat, DOFs, LinIter
      LOGICAL   ::     First, Stat, Direct = .FALSE., &
                       Iterative = .FALSE., NewSystem
      LOGICAL :: Factorize, FreeFactorize,FoundFactorize,FoundFreeFactorize
      LOGICAL :: NoPrecRecompute, FoundNoPrecRecompute
      REAL(KIND=dp), TARGET :: Solution(n), Solution_im(n),ForceVector(n)
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL, s, Residual(n), LinConv, ILUTOL
!
      REAL(KIND=dp), POINTER CONTIG ::  SaveRhs(:)
      REAL(KIND=dp), POINTER :: SaveValues(:)
      CHARACTER(:), ALLOCATABLE :: str, Method

      INTERFACE
        SUBROUTINE BlockSolveExt(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
        END SUBROUTINE BlockSolveExt
      END INTERFACE

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

!
      Solution    = 0
      ForceVector = 0
      Residual    = 0
      IPNTR = 0

      DOFs = Solver % Variable % DOFs
      CALL ParallelInitSolve( A, Solution, ForceVector, Residual )

      PMatrix => ParallelMatrix(A)
      PN = PMatrix % NumberOfRows

!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix. A     |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G'.) NEV is the number of eigenvalues to be       |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                     N <= MAXN,                     |
!     |                   NEV <= MAXNEV,                   |
!     |               NEV + 1 <= NCV <= MAXNCV             |
!     %----------------------------------------------------%
!
      NCV = ListGetInteger( Solver % Values, 'Eigen System Lanczos Vectors', stat )
      IF ( .NOT. stat ) NCV = 3*NEIG + 1

      IF ( NCV <=  NEIG ) THEN
         CALL Fatal( 'ParallelEigenSolve', &
               'Number of Lanczos vectors must exceed the number of eigenvalues.' )
      END IF

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV,3), &
         WORKEV(3*NCV), V(PN,NCV), CHOOSE(NCV), STAT=istat )

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'ParallelEigenSolve', 'Memory allocation error.' )
      END IF
!
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      TOL = ListGetConstReal( Solver % Values, 'Eigen System Convergence Tolerance', stat )
      IF ( .NOT. stat ) THEN
         TOL = 100 * ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance' )
      END IF
      IDO   = 0
      kinfo = 0
      lWORKL = 3*NCV**2 + 6*NCV
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of DSAUPD is used     |
!     | (IPARAM(7) = 2).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      BMAT  = 'G'
      CALL ArpackSetWhich( Solver % Values, A % Lumped, Mode, Which )

      Maxitr = ListGetInteger( Solver % Values, 'Eigen System Max Iterations', stat )
      IF ( .NOT. stat ) Maxitr = 300
!
      IPARAM = 0
      IPARAM(1) = ishfts
      IPARAM(3) = maxitr
      IPARAM(7) = mode

      SigmaR = 0.0d0
      SigmaI = 0.0d0
      V = 0.0d0

!     Compute LU-factors for (A-\sigma M) (if consistent mass matrix)
!
      Factorize = ListGetLogical( Solver % Values, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Solver % Values, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Solver % Values, &
                'Linear System Free Factorization', FoundFreeFactorize )
      CALL ListAddLogical( Solver % Values,  &
                     'Linear System Free Factorization',.FALSE. )

      Method = ListGetString( Solver % Values, 'Linear System Solver', stat )

      Direct = Method == 'direct'

      DirectMethod = ListGetString( Solver % Values, &
           'Linear System Direct Method', stat )

      IF ( .NOT. A % Lumped ) THEN
         SigmaR = ListGetConstReal( Solver % Values,'Eigen System Shift', stat )
         IF ( SigmaR /= 0.0d0 ) THEN
            A % Values = A % Values - SigmaR * A % MassValues
         END IF

         IF ( Direct ) THEN
            SELECT CASE( DirectMethod )
            CASE( 'mumps' )
            CASE DEFAULT
               Stat = CRS_ILUT(A, 0.0d0)
            END SELECT
         END IF
      END IF

      IF ( .NOT. Direct .OR. DirectMethod /= 'mumps' )  THEN
        LinIter = ListGetInteger( Solver % Values, 'Linear System Max Iterations', stat )
        IF ( .NOT. Stat ) LinIter = 1000
        LinConv = ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance', stat )
        IF ( .NOT. Stat ) LinConv = 1.0D-9
!        Preconditioning:
!        ----------------
        str = ListGetString( Solver % Values, 'Linear System Preconditioning', stat )

        k = 0
        IF ( str == 'ilut' )  THEN
          ILUTOL = ListGetConstReal( Solver % Values, &
               'Linear System ILUT Tolerance' )

          stat = CRS_ILUT( PMatrix, ILUTOL )
        ELSE IF ( SEQL(str, 'ilu') ) THEN
           IF(LEN(str)>=4) k = ICHAR(str(4:4)) - ICHAR('0')
           IF ( k<0 .OR. k>9 ) k = 0
           PMatrix % Cholesky = ListGetLogical( Solver % Values, &
                'Linear System Symmetric ILU', stat )
           stat = CRS_IncompleteLU( PMatrix, k, Solver % Values )
        END IF
      END IF

      IF ( .NOT. ASSOCIATED( A % RHS ) ) THEN
         ALLOCATE( A % RHS(N) )
         A % RHS = 0._dp
      END IF
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      iter = 1
      NewSystem = .TRUE.

      Iterative = ListGetString( Solver % Values, &
        'Linear System Solver', stat ) == 'iterative'

      NoPrecRecompute = ListGetLogical( Solver % Values, 'No Precondition Recompute', FoundNoPrecRecompute )
      IF ( Iterative .AND. NoPrecRecompute ) &
        CALL ListAddLogical(Solver % Values, 'No Precondition Recompute', .FALSE.)

      DO WHILE( ido /= 99 )
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         IF ( A % Symmetric ) THEN
            CALL PDSAUPD ( A % Comm, ido, BMAT, PN, Which, NEIG, TOL, &
              RESID, NCV, V, PN, IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )
         ELSE
            CALL PDNAUPD ( A % Comm, ido, BMAT, PN, Which, NEIG, TOL, &
              RESID, NCV, V, PN, IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )
         END IF
!
         IF (ido == -1 .OR. ido == 1) THEN
            CALL Info( 'ParallelEigenSolve', '.', .TRUE., Level=5 )
            iter = iter + 1
!---------------------------------------------------------------------
!             Perform  y <--- OP*x = inv[M]*A*x   (lumped mass)
!                      ido =-1 inv(A-sigmaR*M)*M*x
!                      ido = 1 inv(A-sigmaR*M)*z
!---------------------------------------------------------------------
            IF ( ido == -1 .OR. A % Lumped ) THEN
               ! ido==-1: need M*x first; lumped always needs M*x regardless of ido
               x => WORKD(IPNTR(1):IPNTR(1)+PN-1)
               b => WORKD(IPNTR(2):IPNTR(2)+PN-1)
               CALL MGmv( A, x, b, .FALSE., .TRUE. )
               DO i=1,PN
                 x(i) = b(i)
               END DO
            END IF

            x => Solution
            b => ForceVector

            CALL PartitionVector(A, x, WORKD(IPNTR(2):IPNTR(2)+PN-1))
            IF ( ido == 1 .AND. .NOT. A % Lumped ) THEN
               CALL PartitionVector(A, b, WORKD(IPNTR(3):IPNTR(3)+PN-1))
            ELSE
               CALL PartitionVector(A, b, WORKD(IPNTR(1):IPNTR(1)+PN-1))
            END IF

            ! Some strategies (such as 'block') may depend on that these are set properly
            ! to reflect the linear problem under study.
            SaveRhs => A % rhs
            A % rhs => ForceVector

            SELECT CASE( Method )
            CASE('multigrid')
              CALL MultiGridSolve( A, x, b, &
                  DOFs, Solver, Solver % MultiGridLevel, NewSystem )
            CASE('iterative')
              CALL ParallelIter( A, A % ParallelInfo, DOFs, &
                     x, b, Solver, A % ParMatrix )
            CASE('block')
              CALL BlockSolveExt( A, x, b, Solver )
            CASE ('direct')
              CALL DirectSolver( A, x, b, Solver )
            CASE DEFAULT
              CALL Fatal('ParallelEigenSolve','Unknown linear system method: '//TRIM(Method))
            END SELECT
            CALL ParallelInitSolve( A, x, b, Residual )

            A % rhs => SaveRhs

            CALL ParallelVector(A, WORKD(IPNTR(2):IPNTR(2)+PN-1), x)
         ELSE IF (ido == 2) THEN
!
!           %-----------------------------------------%
!           |         Perform  y <--- M*x.            |
!           | Need the matrix vector multiplication   |
!           | routine here that takes WORKD(IPNTR(1)) |
!           | as the input and returns the result to  |
!           | WORKD(IPNTR(2)).                        |
!           %-----------------------------------------%
!
            IF ( A % Lumped ) THEN
               DO i=0,n-1
                  WORKD( IPNTR(2)+i ) = WORKD( IPNTR(1)+i ) * &
                    A % MassValues( A % Diag(i+1) )
               END DO
            ELSE
               x => WORKD(IPNTR(1):IPNTR(1)+PN-1)
               b => WORKD(IPNTR(2):IPNTR(2)+PN-1)
               CALL MGmv( A, x, b, .FALSE., .TRUE. )
            END IF
         END IF

         IF ( NewSystem .AND. ido /= 2 ) THEN
            IF ( Iterative ) THEN
               CALL ListAddLogical( Solver % Values,  'No Precondition Recompute', .TRUE. )
            ELSE
               CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
            END IF
            NewSystem = .FALSE.
         END IF
      END DO

      IF ( FoundFactorize ) THEN
        CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', Factorize )
      ELSE
        CALL ListRemove( Solver % Values, 'Linear System Refactorize' )
      END IF

      IF ( .NOT. FoundFreeFactorize ) THEN
        CALL ListRemove( Solver % Values, 'Linear System Free Factorization' )
      ELSE
        CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', FreeFactorize )
      END IF

      IF ( FoundNoPrecRecompute ) THEN
        CALL ListAddLogical( Solver % Values, 'No Precondition Recompute', NoPrecRecompute )
      ELSE
        CALL ListRemove( Solver % Values, 'No Precondition Recompute' )
      END IF

!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      IF ( kinfo /= 0 ) THEN
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD  |
!        %--------------------------%
!
         WRITE( Message, * ) 'Error with DNAUPD, info = ',kinfo
         CALL Fatal( 'ParallelEigenSolve', Message )
!
      END IF
!
!     %-------------------------------------------%
!     | No fatal errors occurred.                 |
!     | Post-Process using DSEUPD.                |
!     |                                           |
!     | Computed eigenvalues may be extracted.    |
!     |                                           |
!     | Eigenvectors may also be computed now if  |
!     | desired.  (indicated by rvec = .true.)    |
!     %-------------------------------------------%
!
      D = 0.0d0
      IF ( A % Symmetric ) THEN
         CALL pDSEUPD ( A % Comm, .TRUE., 'A', Choose, D, V, PN, SigmaR,  &
            BMAT, PN, Which, NEIG, TOL, RESID, NCV, V, PN, &
            IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
      ELSE
         CALL pDNEUPD ( A % Comm, .TRUE., 'A', Choose, D, D(1,2), &
            V, PN, SigmaR, SigmaI, WORKEV, BMAT, PN, &
            Which, NEIG, TOL, RESID, NCV, V, PN, &
            IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
      END IF

!     %----------------------------------------------%
!     | Eigenvalues are returned in the First column |
!     | of the two dimensional array D and the       |
!     | corresponding eigenvectors are returned in   |
!     | the First NEV columns of the two dimensional |
!     | array V if requested.  Otherwise, an         |
!     | orthogonal basis for the invariant subspace  |
!     | corresponding to the eigenvalues in D is     |
!     | returned in V.                               |
!     %----------------------------------------------%

      IF (IERR /= 0) THEN
!
!        %------------------------------------%
!        | Error condition:                   |
!        | Check the documentation of DNEUPD. |
!        %------------------------------------%
!
         WRITE( Message, * ) ' Error with DNEUPD, info = ', IERR
         CALL Fatal( 'ParallelEigenSolve', Message )
      END IF
!
!     %------------------------------------------%
!     | Print additional convergence information |
!     %------------------------------------------%
!
      IF ( kinfo == 1 ) THEN
         CALL Fatal( 'ParallelEigenSolve', 'Maximum number of iterations reached.' )
      ELSE IF ( kinfo == 3 ) THEN
         CALL Fatal( 'ParallelEigenSolve', &
            'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
      END IF
!
!     Sort the eigenvalues to ascending order:
!        ----------------------------------------
      ALLOCATE( Perm(NEIG) )
      Perm = (/ (i, i=1,NEIG) /)
      DO i=1,NEIG
         EigValues(i) = CMPLX( D(i,1), D(i,2),KIND=dp )
      END DO
      CALL SortC( NEIG, EigValues, Perm )
!
!     Extract the values to ELMER structures:
!     -----------------------------------------
      CALL Info( 'ParallelEigenSolve', ' ', Level=4 )
      CALL Info( 'ParallelEigenSolve', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
      CALL Info( 'ParallelEigenSolve', ' ', Level=4 )
      WRITE( Message, * ) 'The convergence criterion is ', TOL
      CALL Info( 'ParallelEigenSolve', Message, Level=4 )
      WRITE( Message, * ) ' The number of converged Ritz values is ', IPARAM(5)
      CALL Info( 'ParallelEigenSolve', Message, Level=4 )
      CALL Info( 'ParallelEigenSolve', ' ', Level=4 )
      CALL Info( 'ParallelEigenSolve', 'Computed Eigen Values: ', Level=3 )
      CALL Info( 'ParallelEigenSolve', '--------------------------------', Level=3 )

      ! Restore matrix values, if modified when using shift:
      ! ---------------------------------------------------
      IF ( SigmaR /= 0.0d0 ) THEN
         A % Values = A % Values + SigmaR * A % MassValues
      END IF

      k = 1
      DO i=1,NEIG
        p = Perm(i)
        WRITE( Message, * ) i,EigValues(i)
        CALL Info( 'ParallelEigenSolve', Message, Level=3 )

        k = 1
        DO j=1,p-1
           IF ( D(j,2) == 0 ) THEN
              k = k + 1
           ELSE
              k = k + 2
           END IF
        END DO

        Residual = 0.0d0
        Solution = 0.0d0; Solution_im = 0.0d0

        DO j=1,PN
          IF ( D(p,2) /= 0.0d0 ) THEN
            A % ParMatrix % SplittedMatrix % TmpXVec(j) = V(j,k)
          ELSE
            A % ParMatrix % SplittedMatrix % TmpXVec(j) = V(j,k)
          END IF
        END DO

        CALL ParallelUpdateResult( A, Solution, Residual )

        DO j=1,PN
          IF ( D(p,2) /= 0.0d0 ) THEN
            A % ParMatrix % SplittedMatrix % TmpXVec(j) = V(j,k+1)
          ELSE
            A % ParMatrix % SplittedMatrix % TmpXVec(j) = 0.0d0
          END IF
        END DO
        CALL ParallelUpdateResult( A, Solution_im, Residual )


        DO j=1,N
          EigVectors(i,j) = CMPLX( Solution(j), Solution_im(j), KIND=dp )
        END DO

      END DO
      CALL Info( 'ParallelEigenSolve', '--------------------------------',Level=3 )

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm )
#else
      CALL Fatal( 'ParallelEigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
#endif
!
!------------------------------------------------------------------------------
     END SUBROUTINE ParallelArpackEigenSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ParallelArpackEigenSolveComplex( Solver,Matrix,N,NEIG,EigValues,EigVectors )
!------------------------------------------------------------------------------
!> Parallel solution of eigenvalue problems using PARPACK library,
!> complex-valued version. The structure follows ArpackEigenSolveComplex.
!
!  Oct 21 2000, Juha Ruokolainen
!
!\Original Authors
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
! FILE: ndrv3.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
!------------------------------------------------------------------------------

#ifdef PARALLEL_FOR_REAL
#  if defined(ELMER_HAVE_MPI_MODULE)
      USE mpi
#  endif
#endif

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: Matrix
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef PARALLEL_FOR_REAL
#  if defined(ELMER_HAVE_MPIF_HEADER)
      INCLUDE "mpif.h"
#  endif
#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      TYPE(Matrix_t), POINTER :: PMatrix
      INTEGER :: IPARAM(11), IPNTR(14), Perm(NEIG)
      LOGICAL, ALLOCATABLE :: Choose(:)
      COMPLEX(KIND=dp) :: D(NEIG)
      REAL(KIND=dp), ALLOCATABLE :: RWORK(:)
      REAL(KIND=dp), ALLOCATABLE, TARGET :: Residual(:), Solution(:), ForceVector(:)

      COMPLEX(KIND=dp), POINTER CONTIG :: x(:), b(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET :: WORKL(:), WORKEV(:), V(:,:), WORKD(:), RESID(:), xx(:)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      CHARACTER ::     BMAT*1, Which*2
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, p, pn, IERR, iter, &
                       maxitr, mode, istat, DOFs, MinPn
      LOGICAL   ::     Stat, Iterative, NewSystem, ScaleSystem, NullifyGuess
      LOGICAL   ::     Factorize, FreeFactorize, FoundFactorize, FoundFreeFactorize, StandardForm
      LOGICAL   ::     NoPrecRecompute, FoundNoPrecRecompute

      COMPLEX(KIND=dp) :: Sigma
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL
      REAL(KIND=dp), POINTER CONTIG :: SaveRhs(:)

      CHARACTER(:), ALLOCATABLE :: Method

      INTEGER :: me
      TYPE(NeighbourList_t), POINTER :: OwnerList(:)

      TYPE(ValueList_t), POINTER :: Params

      CHARACTER(*), PARAMETER :: Caller = 'ParallelArpackEigenSolveComplex'

      INTERFACE
        SUBROUTINE BlockSolveExt(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
        END SUBROUTINE BlockSolveExt
      END INTERFACE

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

      CALL Info(Caller,'Parallel eigen system solution with complex matrices!',Level=6)

      Params => Solver % Values

      ALLOCATE( Residual(2*n), Solution(2*n), ForceVector(2*n), STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( Caller, 'Memory allocation error.' )
      END IF

      Solution    = 0
      ForceVector = 1
      Residual    = 0
      IPNTR = 0

      DOFs = Solver % Variable % DOFs
      CALL ParallelInitSolve( Matrix, Solution, ForceVector, Residual )

      PMatrix => ParallelMatrix(Matrix)
      pn = PMatrix % NumberOFRows / 2

!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix. A     |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G'). NEV is the number of eigenvalues to be       |
!     | approximated. The user can modify NEV, NCV, WHICH  |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum. However, the      |
!     | following conditions must be satisfied:            |
!     |               NEV + 2 <= NCV                       |
!     |                   NCV <= N_p, with N_p local       |
!     |                          size (see below)          |
!     %----------------------------------------------------%
!

      ! pzneupd checks NCV against the local size pn (ncv > n gives ierr = -3),
      ! so NCV is limited by the smallest partition. It seems that this check
      ! should be corrected on PARPACK side, but now we just inform what can
      ! go wrong:
      MinPn = ParallelReduction(pn,1)
      NCV = ListGetInteger( Params, 'Eigen System Lanczos Vectors', stat )
      IF ( stat ) THEN
        IF ( NCV > MinPn ) THEN
          CALL Fatal( Caller, 'Eigen System Lanczos Vectors ('//I2S(NCV)// &
              ') exceeds the size of the smallest partition ('//I2S(MinPn)// &
              '), which PARPACK (pzneupd) does not allow' )
        END IF
      ELSE
        NCV = MIN( 3*NEIG + 1, MinPn )
        IF ( NCV < 3*NEIG + 1 ) THEN
          CALL Warn( Caller, 'Default number of Lanczos vectors reduced from '// &
              I2S(3*NEIG+1)//' to '//I2S(NCV)//' by the size of the smallest partition' )
          CALL Warn( Caller, 'Convergence may suffer, consider fewer or better balanced partitions' )
        END IF
      END IF

      ! pzneupd also requires NCV >= NEIG+2
      IF ( NCV <= NEIG + 1 ) THEN
        IF ( .NOT. stat .AND. NCV == MinPn ) THEN
          CALL Fatal( Caller, 'Smallest partition ('//I2S(MinPn)// &
              ' rows) is too small for '//I2S(NEIG)//' eigenvalues, NCV >= NEIG+2 is needed' )
        ELSE
          CALL Fatal( Caller, 'Number of Lanczos vectors must be at least '// &
              'the number of eigenvalues plus two.' )
        END IF
      END IF

!
!     %--------------------------------------------------%
!     | The work array WORKL is used in PZNAUPD as       |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below. The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine       |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in PZNAUPD to start the Arnoldi        |
!     | iteration.                                       |
!     %--------------------------------------------------%
!
      lWORKL = 3*NCV**2 + 6*NCV

      TOL = ListGetConstReal( Params, 'Eigen System Convergence Tolerance', stat )
      IF ( .NOT. stat ) THEN
         TOL = 100 * ListGetConstReal( Params, 'Linear System Convergence Tolerance' )
      END IF
      IDO   = 0
      kinfo = 0

      ALLOCATE( WORKL(lWORKL), rwork(NCV), &
         WORKEV(2*NCV), V(pn,NCV), CHOOSE(NCV), WORKD(3*pn), RESID(pn), xx(pn), STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( Caller, 'Memory allocation error.' )
      END IF

      WORKL = 0
      D = 0
      WORKEV = 0
      V = 0.0d0
      CHOOSE=.FALSE.
      WORKD=0
      RESID = 0
      xx = 0

!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of PZNAUPD is used     |
!     | by default. All these options may be              |
!     | changed by the user. For details, see the         |
!     | documentation in PZNAUPD.                         |
!     %---------------------------------------------------%
!
      IF (Matrix % Lumped) THEN
        CALL Warn(Caller, 'No implementation for a lumped matrix in Mode 2')
      END IF

      CALL ArpackSetWhich( Params, .FALSE., Mode, Which )

      Maxitr = ListGetInteger( Params, 'Eigen System Max Iterations', stat, DefValue = 300 )

      IPARAM = 0
      IPARAM(1) = 1
      IPARAM(3) = maxitr

      ! The generalized mode of ARPACK uses the M-inner product and so needs M to be
      ! Hermitian positive semidefinite. If this cannot be guaranteed, use instead the
      ! Euclidean inner product and apply ARPACK to the standard eigenvalue problem
      ! for OP = (A-sigma*M)^{-1}*M. The eigenvalues are then recovered as
      ! lambda = sigma + 1/theta, with theta an eigenvalue of OP.
      StandardForm = ListGetLogical(Params, 'Eigen System Standard Form', stat)
      IF (StandardForm) THEN
        BMAT = 'I'
        Mode = 1
      ELSE
        BMAT  = 'G'
      END IF
      CALL Info(Caller, 'The ARPACK routine pznaupd will be called in mode '//I2S(mode), Level=12)
      IPARAM(7) = mode

!     Compute LU-factors for (A-\sigma M) (if consistent mass matrix):
!     ----------------------------------------------------------------
      Factorize = ListGetLogical( Params, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Params, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Params, &
                'Linear System Free Factorization', FoundFreeFactorize )
      CALL ListAddLogical( Params,  &
                     'Linear System Free Factorization',.FALSE. )

      IF (ListGetLogical(Params, 'Linear System Skip Scaling', stat)) THEN
        CALL Info(Caller, 'This time skipping scaling', Level=20)
        ScaleSystem = .FALSE.
      ELSE
        ScaleSystem = ListGetLogical(Params, 'Linear System Scaling', stat, DefValue = .TRUE.)
        ! Row equilibration wouldn't lead to a Hermitian matrix
        IF (ListGetLogical(Params, 'Linear System Row Equilibration', stat)) THEN
          CALL Fatal(Caller, 'Set Linear System Row Equilibration = False for the eigen solution')
        END IF
      END IF

      Method = ListGetString( Params, 'Linear System Solver', stat )
      IF (.NOT. stat) Method = 'direct'
      Iterative = Method == 'iterative'

      ! Mode 3:
      SigmaR = ListGetConstReal( Params,'Eigen System Shift', stat )
      SigmaI = ListGetConstReal( Params,'Eigen System Shift Im', stat )
      Sigma = CMPLX(SigmaR,SigmaI,KIND=dp)

      IF ( ABS(Sigma) > AEPS ) CALL ShiftComplexMatrix( Matrix, -Sigma )

      IF (ScaleSystem) CALL ScaleLinearSystem(Solver, Matrix, Shift = Sigma)

      ! The parallel matrix holds a copy of the values. Update it so that the shifted
      ! and scaled matrices are used in the matrix-vector products:
      CALL ParallelInitSolve( Matrix, Solution, ForceVector, Residual )

      IF ( .NOT. ASSOCIATED( Matrix % RHS ) ) THEN
         ALLOCATE(Matrix % RHS(2*n)); Matrix % RHS = 0._dp
      END IF
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      iter = 0
      NewSystem = .TRUE.

      NoPrecRecompute = ListGetLogical( Params, 'No Precondition Recompute', FoundNoPrecRecompute )
      IF ( Iterative .AND. NoPrecRecompute ) &
        CALL ListAddLogical(Params, 'No Precondition Recompute', .FALSE.)

      NullifyGuess = ListGetLogical(Params, 'Linear System Nullify Guess', stat)

      me = ParEnv % MyPe
      OwnerList =>  Matrix % ParallelInfo % NeighbourList

      DO WHILE( ido /= 99 )
!
!        %---------------------------------------------%
!        | Repeatedly call the routine PZNAUPD and     |
!        | take actions indicated by parameter IDO     |
!        | until either convergence is indicated or    |
!        | maxitr has been exceeded.                   |
!        %---------------------------------------------%
         CALL PZNAUPD( Matrix % Comm, ido, BMAT, pn, Which, NEIG, TOL, &
           RESID, NCV, V, pn, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, kinfo )

         IF (ido == -1 .OR. ido == 1) THEN
            iter = iter + 1
            IF(InfoActive(20)) THEN
              CALL Info(Caller,'Arpack reverse communication calls: '//I2S(iter))
            ELSE
              CALL Info(Caller, '.', .TRUE., Level=5)
            END IF
!---------------------------------------------------------------------
!           Perform  y = OP*x, with
!                    OP*x = inv(A-sigma*M)*M*x for mode 3 and ido =-1
!                    OP*x = inv(A-sigma*M)*z, with z returned by pznaupd, for mode 3 and ido = 1
!---------------------------------------------------------------------
            IF ( ido == 1 .AND. .NOT. StandardForm ) THEN
               b => WORKD(IPNTR(3):IPNTR(3)+pn-1)
            ELSE
               x => WORKD(IPNTR(1):IPNTR(1)+pn-1)
               b => xx
               CALL CMGmv( Matrix, x, b, .FALSE., .TRUE. )
            END IF
            x => WORKD(IPNTR(2):IPNTR(2)+pn-1)

            ! Map the owned values to the local real-valued vectors of the partition
            ! (the right-hand side vanishes at the dofs owned by other partitions):
            ForceVector = 0._dp
            Solution    = 0._dp
            j = 0
            DO i=0,n-1
              IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                j = j + 1
                IF ( .NOT. NullifyGuess ) THEN
                  Solution(2*i+1) = REAL(x(j))
                  Solution(2*i+2) = AIMAG(x(j))
                END IF
                ForceVector(2*i+1) = REAL(b(j))
                ForceVector(2*i+2) = AIMAG(b(j))
              END IF
            END DO

            ! Some strategies (such as 'block') may depend on that these are set properly
            ! to reflect the linear problem under study.
            SaveRhs => Matrix % rhs
            Matrix % rhs => ForceVector

            SELECT CASE( Method )
            CASE('multigrid')
              CALL MultiGridSolve( Matrix, Solution, ForceVector, &
                  DOFs, Solver, Solver % MultiGridLevel, NewSystem )
            CASE('iterative')
              CALL ParallelIter( Matrix, Matrix % ParallelInfo, DOFs, &
                  Solution, ForceVector, Solver, Matrix % ParMatrix )
            CASE('block')
              CALL BlockSolveExt( Matrix, Solution, ForceVector, Solver )
            CASE('direct')
              CALL DirectSolver( Matrix, Solution, ForceVector, Solver )
            CASE DEFAULT
              CALL Fatal(Caller,'Unknown linear system method: '//TRIM(Method))
            END SELECT

            Matrix % rhs => SaveRhs

            j = 0
            DO i=0,n-1
              IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                j = j + 1
                x(j) = CMPLX(Solution(2*i+1), Solution(2*i+2), KIND=dp)
              END IF
            END DO

         ELSE IF (ido == 2) THEN
!
!           %-----------------------------------------%
!           |         Perform  y <--- M*x.            |
!           | Need the matrix vector multiplication   |
!           | routine here that takes WORKD(IPNTR(1)) |
!           | as the input and returns the result to  |
!           | WORKD(IPNTR(2)).                        |
!           %-----------------------------------------%
!
            x => WORKD(IPNTR(1):IPNTR(1)+pn-1)
            b => WORKD(IPNTR(2):IPNTR(2)+pn-1)
            CALL CMGmv( Matrix, x, b, .FALSE., .TRUE. )
         END IF

         IF ( NewSystem .AND. ido /= 2 ) THEN
            IF ( Iterative ) THEN
               CALL ListAddLogical( Params,  'No Precondition Recompute', .TRUE. )
            ELSE
               CALL ListAddLogical( Params, 'Linear System Refactorize', .FALSE. )
            END IF
            NewSystem = .FALSE.
         END IF
      END DO

      IF ( FoundFactorize ) THEN
        CALL ListAddLogical( Params, 'Linear System Refactorize', Factorize )
      ELSE
        CALL ListRemove( Params, 'Linear System Refactorize' )
      END IF

      IF ( FoundFreeFactorize ) THEN
        CALL ListAddLogical( Params, 'Linear System Free Factorization', FreeFactorize )
      ELSE
        CALL ListRemove( Params, 'Linear System Free Factorization' )
      END IF

      IF ( FoundNoPrecRecompute ) THEN
        CALL ListAddLogical( Params, 'No Precondition Recompute', NoPrecRecompute )
      ELSE
        CALL ListRemove( Params, 'No Precondition Recompute' )
      END IF

!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      IF ( kinfo /= 0 ) THEN
!
!        %---------------------------%
!        | Error message, check the  |
!        | documentation in PZNAUPD  |
!        %---------------------------%
        SELECT CASE(kinfo)
        CASE(1)
          CALL Fatal( Caller, 'Maximum number of iterations reached.' )
        CASE(3)
          CALL Fatal( Caller, 'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
        CASE DEFAULT
          WRITE( Message, * ) 'Error with PZNAUPD, info = ',kinfo
          CALL Fatal( Caller, Message )
        END SELECT
      END IF
!
!     %-------------------------------------------%
!     | No fatal errors occurred.                 |
!     | Post-Process using PZNEUPD.               |
!     |                                           |
!     | Computed eigenvalues may be extracted.    |
!     |                                           |
!     | Eigenvectors may also be computed now if  |
!     | desired.  (indicated by rvec = .true.)    |
!     %-------------------------------------------%
!
      D = 0.0_dp
      CALL PZNEUPD( Matrix % Comm, .TRUE., 'A', Choose, D, V, pn, Sigma, WORKEV, BMAT, pn, &
         Which, NEIG, TOL, RESID, NCV, V, PN, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, IERR )

!     %----------------------------------------------%
!     | Eigenvalues are returned in the one-         |
!     | dimensional array D and the                  |
!     | corresponding eigenvectors are returned in   |
!     | the first NEV columns of the two dimensional |
!     | array V if requested.  Otherwise, an         |
!     | orthogonal basis for the invariant subspace  |
!     | corresponding to the eigenvalues in D is     |
!     | returned in V.                               |
!     %----------------------------------------------%

      IF (IERR /= 0) THEN
!
!        %------------------------------------%
!        | Error condition:                   |
!        | Check the documentation of PZNEUPD |
!        %------------------------------------%
!
         WRITE( Message, * ) ' Error with PZNEUPD, info = ', IERR
         CALL Fatal( Caller, Message )
      END IF

      ! Transform the eigenvalues of OP to those of the original problem:
      IF (StandardForm) THEN
        DO i=1,NEIG
          IF (ABS(D(i)) > TINY(1.0_dp)) THEN
            D(i) = Sigma + 1.0_dp/D(i)
          ELSE
            D(i) = CMPLX(HUGE(1.0_dp), 0.0_dp, KIND=dp)
          END IF
        END DO
      END IF

!     Sort the eigenvalues:
!     ---------------------
      DO i=1,NEIG
         EigValues(i) = D(i)
      END DO
      CALL EigenSystemSorting( Params, NEIG, Perm, EigValues )
!
!     Extract the values to ELMER structures:
!     -----------------------------------------
      CALL Info( Caller, ' ', Level=4 )
      CALL Info( Caller, 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
      CALL Info( Caller, ' ', Level=4 )
      WRITE( Message,'(A,ES12.3)') 'Convergence criterion is: ', TOL
      CALL Info( Caller, Message, Level=7 )
      CALL Info( Caller,'Number of eigensystem iterations is: '//I2S(iter),Level=4)
      CALL Info( Caller,'Number of converged Ritz values is: '//I2S(IPARAM(5)),Level=4)
      CALL Info( Caller, ' ', Level=7 )
      CALL Info( Caller, 'Computed Eigen Values: ', Level=4 )
      CALL Info( Caller, '--------------------------------', Level=7 )

      EigVectors = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      DO i=1,NEIG
        p = Perm(i)
        WRITE( Message,'(I0,A,2ES15.6)') i,': ',EigValues(i)
        CALL Info( Caller, Message, Level=4 )

        DO j=0,pn-1
          Matrix % ParMatrix % SplittedMatrix % TmpXVec(2*j+1) = REAL(V(j+1,p))
          Matrix % ParMatrix % SplittedMatrix % TmpXVec(2*j+2) = AIMAG(V(j+1,p))
        END DO

        CALL ParallelUpdateResult( Matrix, Solution, Residual )

        DO j=0,n-1
          EigVectors(i,j+1) = CMPLX( Solution(2*j+1), Solution(2*j+2), KIND=dp )
        END DO
      END DO

      ! Note that EigVectors need not be the same array as Solver % Variable % EigenVectors:
      IF (ScaleSystem) CALL BackScaleLinearSystem( Solver, Matrix, EigenVectors = EigVectors )

      ! Restore matrix values, if modified when using shift:
      ! ---------------------------------------------------
      IF ( ABS(Sigma) > AEPS ) CALL ShiftComplexMatrix( Matrix, Sigma )

      IF ( ListGetLogical( Params, 'Eigen System Compute Residuals', stat ) ) THEN
        CALL Warn( Caller, 'Eigen System Compute Residuals is not yet available in parallel' )
      END IF
      CALL Info( Caller, '--------------------------------',Level=4 )

      DEALLOCATE( WORKL, rwork, WORKEV, V, CHOOSE, WORKD, RESID, xx, &
          Residual, Solution, ForceVector )

      CALL Info(Caller,'Finished eigen system solution!',Level=8)
#else
      CALL Fatal( Caller, 'Arpack Eigen System Solver not available.' )
#endif
#endif
!
!------------------------------------------------------------------------------
     END SUBROUTINE ParallelArpackEigenSolveComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE MGmv( A, x, b, Update, UseMass )
!------------------------------------------------------------------------------
       REAL(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update, UseMass
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       LOGICAL :: mass, updt
!------------------------------------------------------------------------------
       mass = .FALSE.
       IF (PRESENT(UseMass)) Mass = UseMass
       updt = .FALSE.
       IF (PRESENT(Update)) updt = Update

       CALL ParallelMatrixVector( A,x,b,updt,mass )
!------------------------------------------------------------------------------
    END SUBROUTINE MGmv
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE CMGmv( A, x, b, Update, UseMass )
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
       LOGICAL, OPTIONAL :: Update, UseMass
       TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
       LOGICAL :: mass, updt
       INTEGER :: i
       REAL(KIND=dp), ALLOCATABLE :: xx(:), bb(:)
!------------------------------------------------------------------------------
       mass = .FALSE.
       IF (PRESENT(UseMass)) Mass = UseMass
       updt = .FALSE.
       IF (PRESENT(Update)) updt = Update

       ALLOCATE(xx(2*size(x)),bb(2*size(b)))
       DO i=0,size(x)-1
         xx(2*i+1) = REAL(x(i+1))
         xx(2*i+2) = AIMAG(x(i+1))
       END DO

       bb = 0._dp
       CALL ParallelMatrixVector( A,xx,bb,updt,mass )

       DO i=0,size(x)-1
         b(i+1) = CMPLX(bb(2*i+1),bb(2*i+2),KIND=dp)
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE CMGmv
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Compute A % Values = A % Values + s * A % MassValues for a complex-valued
!> matrix in the real-valued storage where the complex entry a = ar + i ai
!> is represented by the block [ar -ai; ai ar].
!------------------------------------------------------------------------------
    SUBROUTINE ShiftComplexMatrix( A, s )
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       COMPLEX(KIND=dp), INTENT(IN) :: s
!------------------------------------------------------------------------------
       COMPLEX(KIND=dp) :: c, m
       INTEGER :: i, j
!------------------------------------------------------------------------------
       DO i=1,A % NumberOfRows,2
         DO j=A % Rows(i),A % Rows(i+1)-1,2
           c = CMPLX( A % Values(j), -A % Values(j+1), KIND=dp )
           m = CMPLX( A % MassValues(j), -A % MassValues(j+1), KIND=dp )
           c = c + s * m
           A % Values(j) = REAL(c)
           A % Values(j+1) = -AIMAG(c)
         END DO
       END DO
       DO i=2,A % NumberOfRows,2
         DO j=A % Rows(i),A % Rows(i+1)-1,2
           c = CMPLX( A % Values(j+1), A % Values(j), KIND=dp )
           m = CMPLX( A % MassValues(j+1), A % MassValues(j), KIND=dp )
           c = c + s * m
           A % Values(j) = AIMAG(c)
           A % Values(j+1) = REAL(c)
         END DO
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE ShiftComplexMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE ParallelEigenSolve
!------------------------------------------------------------------------------


!> \}
