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

   USE CRSMatrix
   USE IterSolve
   USE Multigrid
   USE ParallelUtils

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
     SUBROUTINE ParallelArpackEigenSolve( Solver,A,N,NEIG,EigValues,EigVectors )
!------------------------------------------------------------------------------

! the suit the needs of ELMER.
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

      IMPLICIT NONE


      TYPE(Matrix_t), POINTER :: A
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#define PARALLEL_FOR_REAL
#ifdef PARALLEL_FOR_REAL
INCLUDE "mpif.h"
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
      REAL(KIND=dp), TARGET :: Solution(n), Solution_im(n),ForceVector(n)
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL, s, Residual(n), LinConv, ILUTOL
!
      REAL(KIND=dp), POINTER CONTIG ::  SaveRhs(:)
      REAL(KIND=dp), POINTER :: SaveValues(:)
      CHARACTER(LEN=MAX_NAME_LEN) :: str, Method

      INTEGER :: me
      TYPE(NeighbourList_t), POINTER :: OwnerList(:)

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

!
      Solution    = 0
      ForceVector = 0
      Residual    = 0

      DOFs = Solver % Variable % DOFs
      CALL ParallelInitSolve( A, Solution, ForceVector, Residual )

      PMatrix => ParallelMatrix(A) 
      PN = PMatrix % NumberOFRows

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
         CALL Fatal( 'EigenSolve', & 
               'Number of Lanczos vectors must exceed the number of eigenvalues.' )
      END IF

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV,3), &
         WORKEV(3*NCV), V(PN,NCV), CHOOSE(NCV), STAT=istat )

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'EigenSolve', 'Memory allocation error.' )
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
      IF ( A % Lumped ) THEN
         Mode  =  2
         SELECT CASE( ListGetString(Solver % Values,'Eigen System Select',stat) )
         CASE( 'smallest magnitude' )
              Which = 'SM'
         CASE( 'largest magnitude')
              Which = 'LM'
         CASE( 'smallest real part')
              Which = 'SR'
         CASE( 'largest real part')
              Which = 'LR'
         CASE( 'smallest imag part' )
              Which = 'SI'
         CASE( 'largest imag part' )
              Which = 'LI'
         CASE DEFAULT
              Which = 'SM'
         END SELECT
      ELSE
         Mode  = 3
         SELECT CASE( ListGetString(Solver % Values,'Eigen System Select',stat) )
         CASE( 'smallest magnitude' )
              Which = 'LM'
         CASE( 'largest magnitude')
              Which = 'SM'
         CASE( 'smallest real part')
              Which = 'LR'
         CASE( 'largest real part')
              Which = 'SR'
         CASE( 'smallest imag part' )
              Which = 'LI'
         CASE( 'largest imag part' )
              Which = 'SI'
         CASE DEFAULT
              Which = 'LM'
         END SELECT
      END IF

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
                'Linear System Refactorize', FoundFreeFactorize )
      CALL ListAddLogical( Solver % Values,  &
                     'Linear System Free Factorization',.FALSE. )

      Method = ListGetString( Solver % Values, &
           'Linear System Solver', stat )

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
!        Preconditiong:
!        --------------
        str = ListGetString( Solver % Values, 'Linear System Preconditioning', stat )

        IF ( str == 'ilut' )  THEN
          ILUTOL = ListGetConstReal( Solver % Values, &
               'Linear System ILUT Tolerance' )

          stat = CRS_ILUT( PMatrix, ILUTOL )
        ELSE IF ( SEQL(str, 'ilu') ) THEN
           k = ICHAR(str(4:4)) - ICHAR('0')
           IF ( k<0 .OR. k>9 ) k = 0
           PMatrix % Cholesky = ListGetLogical( Solver % Values, &
                'Linear System Symmetric ILU', stat )
           stat = CRS_IncompleteLU( PMatrix, k )
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

      stat = ListGetLogical(Solver % Values,  'No Precondition Recompute', stat  )
      IF ( Iterative .AND. Stat ) &
        CALL ListAddLogical(Solver % Values, 'No Precondition Recompute', .FALSE.)

      me = ParEnv % MyPe
      OwnerList =>  A % ParallelInfo % NeighbourList

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
!           WRITE( Message, * ) ' Arnoldi iteration: ', Iter
!           CALL Info( 'EigenSolve', Message, Level=5 )
            CALL Info( 'EigenSolve', '.', .TRUE., Level=5 )
            iter = iter + 1
!---------------------------------------------------------------------
!             Perform  y <--- OP*x = inv[M]*A*x   (lumped mass)
!                      ido =-1 inv(A-sigmaR*M)*M*x 
!                      ido = 1 inv(A-sigmaR*M)*z
!---------------------------------------------------------------------
            IF ( .NOT. A % Lumped .AND. ido == 1 ) THEN
               x => Solution
               b => ForceVector

               CALL PartitionVector(A, x, WORKD(IPNTR(2):IPNTR(2)+PN-1))
               CALL PartitionVector(A, b, WORKD(IPNTR(3):IPNTR(3)+PN-1))

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
                        x,b, Solver, A % ParMatrix )
               CASE('block')
                 CALL BlockSolveExt( A, x,b, Solver )
               CASE ('direct')
                 CALL DirectSolver( A, x,b, Solver )
               CASE DEFAULT
                 CALL Fatal('EigenSolve','Unknown linear system method: '//TRIM(Method))
               END SELECT
               CALL ParallelInitSolve( A, x, b, Residual )

               A % rhs => SaveRhs

               CALL ParallelVector(A, WORKD(IPNTR(2):IPNTR(2)+PN-1), x)
            ELSE
               x => WORKD(IPNTR(1):IPNTR(1)+PN-1)
               b => WORKD(IPNTR(2):IPNTR(2)+PN-1)
               CALL MGmv( A, x, b, .FALSE., .TRUE. )
               DO i=1,PN
                 x(i) = b(i)
               END DO

               x => Solution
               b => ForceVector

               CALL PartitionVector(A, x, WORKD(IPNTR(2):IPNTR(2)+PN-1))
               CALL PartitionVector(A, b, WORKD(IPNTR(1):IPNTR(1)+PN-1))

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
                 CALL Fatal('EigenSolve','Unknown linear system method: '//TRIM(Method))
               END SELECT
               CALL ParallelInitSolve( A, x, b, Residual )

               A % rhs => SaveRhs

               CALL ParallelVector(A, WORKD(IPNTR(2):IPNTR(2)+PN-1), x)
            END IF
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
         CALL Fatal( 'EigenSolve', Message )
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
         CALL pDSEUPD ( ELMER_COMM_WORLD, .TRUE., 'A', Choose, D, V, PN, SigmaR,  &
            BMAT, PN, Which, NEIG, TOL, RESID, NCV, V, PN, &
            IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
      ELSE
         CALL pDNEUPD ( ELMER_COMM_WORLD, .TRUE., 'A', Choose, D, D(1,2), &
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
         CALL Fatal( 'EigenSolve', Message )
      END IF
!
!     %------------------------------------------%
!     | Print additional convergence information |
!     %------------------------------------------%
!
      IF ( kinfo == 1 ) THEN
         CALL Fatal( 'EigenSolve', 'Maximum number of iterations reached.' )
      ELSE IF ( kinfo == 3 ) THEN
         CALL Fatal( 'EigenSolve', &
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
      CALL Info( 'EigenSolve', ' ', Level=4 )
      CALL Info( 'EigenSolve', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
      CALL Info( 'EigenSolve', ' ', Level=4 )
      WRITE( Message, * ) 'The convergence criterion is ', TOL
      CALL Info( 'EigenSolve', Message, Level=4 )
      WRITE( Message, * ) ' The number of converged Ritz values is ', IPARAM(5)
      CALL Info( 'EigenSolve', Message, Level=4 )
      CALL Info( 'EigenSolve', ' ', Level=4 )
      CALL Info( 'EigenSolve', 'Computed Eigen Values: ', Level=3 )
      CALL Info( 'EigenSolve', '--------------------------------', Level=3 )

      ! Restore matrix values, if modified when using shift:
      ! ---------------------------------------------------
      IF ( SigmaR /= 0.0d0 ) THEN
         A % Values = A % Values + SigmaR * A % MassValues
      END IF

      k = 1
      DO i=1,NEIG
        p = Perm(i)
        WRITE( Message, * ) i,EigValues(i)
        CALL Info( 'EigenSolve', Message, Level=3 )

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
      CALL Info( 'EigenSolve', '--------------------------------',Level=3 )

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm )
#else
      CALL Fatal( 'EigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
#endif
!
!------------------------------------------------------------------------------
     END SUBROUTINE ParallelArpackEigenSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ParallelArpackEigenSolveComplex( Solver,Matrix,N,NEIG,EigValues,EigVectors )
!------------------------------------------------------------------------------

! the suit the needs of ELMER.
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

      IMPLICIT NONE


      TYPE(Matrix_t), POINTER :: Matrix
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#define PARALLEL_FOR_REAL
#ifdef PARALLEL_FOR_REAL
INCLUDE "mpif.h"
#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      TYPE(Matrix_t), POINTER :: PMatrix
      INTEGER :: IPARAM(32), IPNTR(32)
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      REAL(KIND=dp) :: RWORK(2*n)

      COMPLEX(KIND=dp), POINTER CONTIG :: x(:), b(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET :: WORKL(:), D(:), WORKEV(:), V(:,:), WORKD(:), RESID(:), xx(:)

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

      COMPLEX(KIND=dp) :: Sigma, c, m

      LOGICAL :: Factorize, FreeFactorize,FoundFactorize,FoundFreeFactorize
      REAL(KIND=dp), TARGET :: SigmaR, SigmaI, TOL, s, Residual(2*n), Solution(2*n), &
              ForceVector(2*n), LinConv, ILUTOL
!
      REAL(KIND=dp), POINTER :: SaveValues(:)
      CHARACTER(LEN=MAX_NAME_LEN) :: str

      INTEGER :: me
      TYPE(NeighbourList_t), POINTER :: OwnerList(:)

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      Solution    = 0
      ForceVector = 1
      Residual    = 0

      DOFs = Solver % Variable % DOFs
      CALL ParallelInitSolve( Matrix, Solution, ForceVector, Residual )

      PMatrix => ParallelMatrix(Matrix)
      PN = PMatrix % NumberOFRows/2

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
      IF ( .NOT. stat ) NCV = MIN( 3*NEIG + 1, pn )
      NCV = NINT(ParallelReduction(1._dp*NCV,1))

      IF ( NCV <=  NEIG ) THEN
         CALL Fatal( 'EigenSolve', & 
               'Number of Lanczos vectors must exceed the number of eigenvalues.' )
      END IF

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV), &
         WORKEV(3*NCV), V(PN,NCV), CHOOSE(NCV), WORKD(3*pn), RESID(pn), xx(pn), STAT=istat )

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'EigenSolve', 'Memory allocation error.' )
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
      IF ( Matrix % Lumped ) THEN
         Mode  =  2
         SELECT CASE( ListGetString(Solver % Values,'Eigen System Select',stat) )
         CASE( 'smallest magnitude' )
              Which = 'SM'
         CASE( 'largest magnitude')
              Which = 'LM'
         CASE( 'smallest real part')
              Which = 'SR'
         CASE( 'largest real part')
              Which = 'LR'
         CASE( 'smallest imag part' )
              Which = 'SI'
         CASE( 'largest imag part' )
              Which = 'LI'
         CASE DEFAULT
              Which = 'SM'
         END SELECT
      ELSE
         Mode  = 3
         SELECT CASE( ListGetString(Solver % Values,'Eigen System Select',stat) )
         CASE( 'smallest magnitude' )
              Which = 'LM'
         CASE( 'largest magnitude')
              Which = 'SM'
         CASE( 'smallest real part')
              Which = 'LR'
         CASE( 'largest real part')
              Which = 'SR'
         CASE( 'smallest imag part' )
              Which = 'LI'
         CASE( 'largest imag part' )
              Which = 'SI'
         CASE DEFAULT
              Which = 'LM'
         END SELECT
      END IF

      Maxitr = ListGetInteger( Solver % Values, 'Eigen System Max Iterations', stat )
      IF ( .NOT. stat ) Maxitr = 300
!
      IPARAM = 0
      IPARAM(1) = ishfts
      IPARAM(3) = maxitr 
      IPARAM(7) = mode

      SigmaR = 0.0d0
      SigmaI = 0.0d0
      Sigma = 0._dp
      V = 0.0d0

!     Compute LU-factors for (A-\sigma M) (if consistent mass matrix)
!
      Factorize = ListGetLogical( Solver % Values, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Solver % Values, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Solver % Values, &
                'Linear System Refactorize', FoundFreeFactorize )
      CALL ListAddLogical( Solver % Values,  &
                     'Linear System Free Factorization',.FALSE. )

      Direct = ListGetString( Solver % Values, &
           'Linear System Solver', stat ) == 'direct'

      DirectMethod = ListGetString( Solver % Values, &
           'Linear System Direct Method', stat )

      IF ( .NOT. Matrix % Lumped ) THEN
         SigmaR = ListGetConstReal( Solver % Values,'Eigen System Shift', stat )
         SigmaI = ListGetConstReal( Solver % Values,'Eigen System Shift im', stat )
         Sigma = CMPLX(SigmaR,SigmaI,KIND=dp)
         IF ( Sigma /= 0._dp ) THEN
            DO i=1,Matrix % NumberOfRows,2
              DO j=Matrix % Rows(i),Matrix % Rows(i+1)-1,2
                 c = CMPLX( Matrix % Values(j), -Matrix % Values(j+1), KIND=dp)
                 m = CMPLX( Matrix % MassValues(j), -Matrix % MassValues(j+1), KIND=dp)
                 c = c - Sigma * m
                 Matrix % Values(j) = REAL(c)
                 Matrix % Values(j+1) = -AIMAG(c)
              END DO
            END DO
            DO i=2,Matrix % NumberOfRows,2
              DO j=Matrix % Rows(i),Matrix % Rows(i+1)-1,2
                 c = CMPLX( Matrix % Values(j+1), Matrix % Values(j), KIND=dp)
                 m = CMPLX( Matrix % MassValues(j+1), Matrix % MassValues(j), KIND=dp)
                 c = c - Sigma * m
                 Matrix % Values(j) = AIMAG(c)
                 Matrix % Values(j+1) = REAL(c)
              END DO
            END DO
         END IF

         IF ( Direct ) THEN
            SELECT CASE( DirectMethod )
            CASE( 'mumps' )
            CASE DEFAULT
               Stat = CRS_ComplexILUT(Matrix, 0._dp)
            END SELECT
         END IF
      END IF

      IF ( .NOT. Direct .OR. DirectMethod /= 'mumps' )  THEN
        LinIter = ListGetInteger( Solver % Values, 'Linear System Max Iterations', stat )
        IF ( .NOT. Stat ) LinIter = 1000
        LinConv = ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance', stat )
        IF ( .NOT. Stat ) LinConv = 1.0d-9
!        Preconditiong:
!        --------------
        str = ListGetString( Solver % Values, 'Linear System Preconditioning', stat )

        IF ( str == 'ilut' )  THEN
          ILUTOL = ListGetConstReal( Solver % Values, &
               'Linear System ILUT Tolerance' )

          stat = CRS_ComplexILUT( PMatrix, ILUTOL )
        ELSE IF ( SEQL(str, 'ilu') ) THEN
           k = ICHAR(str(4:4)) - ICHAR('0')
           IF ( k<0 .OR. k>9 ) k = 0
           stat = CRS_ComplexIncompleteLU( PMatrix, k )
        END IF
      END IF

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

      Iterative = ListGetString( Solver % Values, &
        'Linear System Solver', stat ) == 'iterative'

      stat = ListGetLogical(Solver % Values,  'No Precondition Recompute', stat  )
      IF ( Iterative .AND. Stat ) &
        CALL ListAddLogical(Solver % Values, 'No Precondition Recompute', .FALSE.)

      me = ParEnv % MyPe
      OwnerList =>  Matrix % ParallelInfo % NeighbourList

      DO WHILE( ido /= 99 )
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
         CALL PZNAUPD ( Matrix % Comm, ido, BMAT, pn, Which, NEIG, TOL, &
           RESID, NCV, V, pn, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, kinfo )

         IF (ido == -1 .OR. ido == 1) THEN
            iter = iter + 1
!           WRITE( Message, * ) ' Arnoldi iteration: ', Iter
!           CALL Info( 'EigenSolve', Message, Level=5 )
            CALL Info( 'EigenSolve', '.', .TRUE., Level=5 )
!---------------------------------------------------------------------
!             Perform  y <--- OP*x = inv[M]*A*x   (lumped mass)
!                      ido =-1 inv(A-sigmaR*M)*M*x 
!                      ido = 1 inv(A-sigmaR*M)*z
!---------------------------------------------------------------------
            IF ( .NOT. Matrix % Lumped .AND. ido == 1 ) THEN
               x => WORKD(IPNTR(2):IPNTR(2)+pn-1)
               b => WORKD(IPNTR(3):IPNTR(3)+pn-1)

               ForceVector = 0._dp
               Solution    = 0._dp
               IF ( Direct .AND. DirectMethod == 'mumps' ) THEN
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        ForceVector(2*i+1) = REAL(b(j))
                        ForceVector(2*i+2) = AIMAG(b(j))
                     END IF
                  END DO
                  CALL DirectSolver( Matrix,Solution,ForceVector,Solver )
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        x(j) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                     END IF
                  END DO
               ELSE IF ( Solver % MultiGridSolver ) THEN
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        Solution(2*i+1)    = REAL(x(j))
                        Solution(2*i+2)    = AIMAG(x(j))
                        ForceVector(2*i+1) = REAL(b(j))
                        ForceVector(2*i+2) = AIMAG(b(j))
                     END IF
                  END DO
                  CALL MultiGridSolve( Matrix, Solution, ForceVector,  DOFs, &
                        Solver, Solver % MultiGridLevel, NewSystem )
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        x(j) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                     END IF
                  END DO
               ELSE
                  DO i=0,pn-1
                    Solution(2*i+1)    = REAL(x(i+1))
                    Solution(2*i+2)    = AIMAG(x(i+1))
                    ForceVector(2*i+1) = REAL(b(i+1))
                    ForceVector(2*i+2) = AIMAG(b(i+1))
                  END DO
                  CALL BiCGParEigen( Matrix,Solution,ForceVector,Residual,LinIter,LinConv )
                  DO i=0,pn-1
                    x(i+1) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                  END DO
               END IF
               CALL ParallelInitSolve( Matrix, Solution, ForceVector, Residual )
            ELSE
               x => WORKD(IPNTR(1):IPNTR(1)+pn-1)
               b => xx
               CALL CMGmv( Matrix, x, b, .FALSE., .TRUE. )

               x => WORKD(IPNTR(2):IPNTR(2)+pn-1)

               ForceVector = 0._dp
               Solution    = 0._dp
               IF ( Direct .AND. DirectMethod == 'mumps' ) THEN
                  j = 0
                  DO i=0,n-1
                    IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                       j = j + 1
                       ForceVector(2*i+1) = REAL(b(j))
                       ForceVector(2*i+2) = AIMAG(b(j))
                    END IF
                  END DO
                  CALL DirectSolver( Matrix,Solution,ForceVector,Solver )
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        x(j) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                     END IF
                  END DO
               ELSE IF ( Solver % MultiGridSolver ) THEN
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        Solution(2*i+1)    = REAL(x(j))
                        Solution(2*i+2)    = AIMAG(x(j))
                        ForceVector(2*i+1) = REAL(b(j))
                        ForceVector(2*i+2) = AIMAG(b(j))
                     END IF
                  END DO
                  CALL MultiGridSolve( Matrix, Solution, ForceVector,  DOFs, &
                        Solver, Solver % MultiGridLevel, NewSystem )
                  j = 0
                  DO i=0,n-1
                     IF ( OwnerList(2*i+1) % Neighbours(1) == me ) THEN
                        j = j + 1
                        x(j) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                     END IF
                  END DO
               ELSE
                  DO i=0,pn-1
                    Solution(2*i+1)    = REAL(x(i+1))
                    Solution(2*i+2)    = AIMAG(x(i+1))
                    ForceVector(2*i+1) = REAL(b(i+1))
                    ForceVector(2*i+2) = AIMAG(b(i+1))
                  END DO
                  CALL BiCGParEigen( Matrix,Solution,ForceVector,Residual,LinIter,LinConv )
                  DO i=0,pn-1
                    x(i+1) = CMPLX(Solution(2*i+1),Solution(2*i+2),KIND=dp)
                  END DO
               END IF
               CALL ParallelInitSolve( Matrix, Solution, ForceVector, Residual )
            END IF
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
            IF ( Matrix % Lumped ) THEN
               DO i=0,n-1
                  WORKD( IPNTR(2)+i ) = WORKD( IPNTR(1)+i ) * &
                    Matrix % MassValues( Matrix % Diag(i+1) )
               END DO
            ELSE
               x => WORKD(IPNTR(1):IPNTR(1)+PN-1)
               b => WORKD(IPNTR(2):IPNTR(2)+PN-1)
               CALL CMGmv( Matrix, x, b, .FALSE., .TRUE. )
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

      IF ( FoundFreeFactorize ) THEN
        CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', FreeFactorize )
      ELSE
        CALL ListRemove( Solver % Values, 'Linear System Free Factorization' )
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
         CALL Fatal( 'EigenSolve', Message )
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

      d = 0.0_dp
      CALL pzNEUPD ( Matrix % Comm, .TRUE., 'A', Choose, D, V, pn, Sigma, WORKEV, BMAT, pn, &
         Which, NEIG, TOL, RESID, NCV, V, PN, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, IERR )


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
         CALL Fatal( 'EigenSolve', Message )
      END IF
!
!     %------------------------------------------%
!     | Print additional convergence information |
!     %------------------------------------------%
!
      IF ( kinfo == 1 ) THEN
         CALL Fatal( 'EigenSolve', 'Maximum number of iterations reached.' )
      ELSE IF ( kinfo == 3 ) THEN
         CALL Fatal( 'EigenSolve', &
            'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
      END IF      
!
!     Sort the eigenvalues to ascending order:
!        ----------------------------------------
      ALLOCATE( Perm(NEIG) )
      Perm = (/ (i, i=1,NEIG) /)
      DO i=1,NEIG
         EigValues(i) = D(i)
      END DO
      CALL SortC( NEIG, EigValues, Perm )
!
!     Extract the values to ELMER structures:
!     -----------------------------------------
      CALL Info( 'EigenSolve', ' ', Level=4 )
      CALL Info( 'EigenSolve', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
      CALL Info( 'EigenSolve', ' ', Level=4 )
      WRITE( Message, * ) 'The convergence criterion is ', TOL
      CALL Info( 'EigenSolve', Message, Level=4 )
      WRITE( Message, * ) ' The number of converged Ritz values is ', IPARAM(5)
      CALL Info( 'EigenSolve', Message, Level=4 )
      CALL Info( 'EigenSolve', ' ', Level=4 )
      CALL Info( 'EigenSolve', 'Computed Eigen Values: ', Level=3 )
      CALL Info( 'EigenSolve', '--------------------------------', Level=3 )

      ! Restore matrix values, if modified when using shift:
      ! ---------------------------------------------------
      IF ( Sigma /= 0._dp ) THEN
        DO i=1,Matrix % NumberOfRows,2
          DO j=Matrix % Rows(i),Matrix % Rows(i+1)-1,2
             c = CMPLX( Matrix % Values(j), -Matrix % Values(j+1), KIND=dp)
             m = CMPLX( Matrix % MassValues(j), -Matrix % MassValues(j+1), KIND=dp)
             c = c + Sigma * m
             Matrix % Values(j) = REAL(c)
             Matrix % Values(j+1) = -AIMAG(c)
          END DO
        END DO
        DO i=2,Matrix % NumberOfRows,2
          DO j=Matrix % Rows(i),Matrix % Rows(i+1)-1,2
             c = CMPLX( Matrix % Values(j+1), Matrix % Values(j), KIND=dp)
             m = CMPLX( Matrix % MassValues(j+1), Matrix % MassValues(j), KIND=dp)
             c = c + Sigma * m
             Matrix % Values(j) = AIMAG(c)
             Matrix % Values(j+1) = REAL(c)
          END DO
        END DO
      END IF

      k = 1
      DO i=1,NEIG
        p = Perm(i)
        WRITE( Message, * ) i,EigValues(i)
        CALL Info( 'EigenSolve', Message, Level=3 )

        DO j=0,pn-1
          Matrix % ParMatrix % SplittedMatrix % TmpXVec(2*j+1) = REAL(V(j+1,p))
          Matrix % ParMatrix % SplittedMatrix % TmpXVec(2*j+2) = AIMAG(V(j+1,p))
        END DO

        CALL ParallelUpdateResult( Matrix, Solution, Residual )

        DO j=0,n-1
          EigVectors(i,j+1) = CMPLX( Solution(2*j+1), Solution(2*j+2), KIND=dp )
        END DO

      END DO
      CALL Info( 'EigenSolve', '--------------------------------',Level=3 )

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm, WORKD, RESID, xx )
#else
      CALL Fatal( 'EigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
#endif
!
!------------------------------------------------------------------------------
     END SUBROUTINE ParallelArpackEigenSolveComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE Jacobi( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A, M
       INTEGER :: Rounds
       REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
       INTEGER :: i,n
!------------------------------------------------------------------------------
       DO i=1,Rounds
          CALL MGmv( A, x, r )
          r(1:n) = b(1:n) - r(1:n)

          r(1:n) = r(1:n) / M % Values(M % Diag)
          x(1:n) = x(1:n) + r(1:n)
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE Jacobi
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE CGParEigen( A, x, b, r, Rounds, Conv )
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A
       INTEGER :: Rounds
       REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
       REAL(KIND=dp) :: alpha,rho,oldrho, Conv
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: M
       INTEGER :: i,n
       REAL(KIND=dp), POINTER CONTIG :: Mx(:),Mb(:),Mr(:), Z(:), P(:), Q(:)
       REAL(KIND=dp) :: RNorm
!------------------------------------------------------------------------------
       M => ParallelMatrix( A, Mx, Mb, Mr )
       n = M % NumberOfRows
       M % RHS(1:n) = b(1:n)

       ALLOCATE( Z(n), P(n), Q(n) )

       CALL MGmv( A, Mx, Mr )
       Mr(1:n) = Mb(1:n) - Mr(1:n)

       DO i=1,Rounds
          Z(1:n) = Mr(1:n)
          CALL CRS_LUSolve( n, M, Z )
          rho = MGdot( n, Mr, Z )

          IF ( i == 1 ) THEN
             P(1:n) = Z(1:n)
          ELSE
             P(1:n) = Z(1:n) + rho * P(1:n) / oldrho
          END IF

          CALL MGmv( A, P, Q )
          alpha  = rho / MGdot( n, P, Q )
          oldrho = rho

          Mx(1:n) = Mx(1:n) + alpha * P(1:n)
          Mr(1:n) = Mr(1:n) - alpha * Q(1:n)

          RNorm = MGnorm( n, Mr ) / MGNorm( n, Mb )
          IF ( RNorm < Conv ) EXIT
       END DO

       WRITE( Message, * ) 'Iters: ', i, RNorm
       CALL Info( 'CGParEigen', Message, Level=4 ) 

       DEALLOCATE( Z, P, Q )

       x(1:n)= Mx(1:n)
       b(1:n)= Mb(1:n)
!------------------------------------------------------------------------------
    END SUBROUTINE CGParEigen
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE BiCGParEigen( A, x, b, r, Rounds, Conv )
!------------------------------------------------------------------------------
       TYPE(Matrix_t), POINTER :: A,M
       INTEGER :: Rounds
       REAL(KIND=dp) :: Conv
       REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
       INTEGER :: i,n
       REAL(KIND=dp) :: RNorm
       REAL(KIND=dp) :: alpha,beta,omega,rho,oldrho
       REAL(KIND=dp), ALLOCATABLE :: Ri(:),P(:),V(:),T(:),T1(:),T2(:),S(:)
       REAL(KIND=dp), POINTER CONTIG :: Mx(:),Mb(:),Mr(:)
!------------------------------------------------------------------------------
       M => ParallelMatrix( A, Mx, Mb, Mr )
       n = M % NumberOfRows
       M % RHS(1:n) = b(1:n)

       CALL MGmv( A, Mx, Mr )
       Mr(1:n) = Mb(1:n) - Mr(1:n)

       ALLOCATE( Ri(n),P(n),V(n),T(n),T1(n),T2(n),S(n) )

       Ri(1:n) = Mr(1:n)
       P(1:n) = 0
       V(1:n) = 0
       omega  = 1
       alpha  = 0
       oldrho = 1

       DO i=1,Rounds
          rho = MGdot( n, Mr, Ri )

          beta = alpha * rho / ( oldrho * omega )
          P(1:n) = Mr(1:n) + beta * (P(1:n) - omega*V(1:n))

          V(1:n) = P(1:n)
          CALL CRS_LUSolve( n, M, V )
          T1(1:n) = V(1:n)
          CALL MGmv( A, T1, V )

          alpha = rho / MGdot( n, Ri, V )

          S(1:n) = Mr(1:n) - alpha * V(1:n)

          T(1:n) = S(1:n)
          CALL CRS_LUSolve( n, M, T )
          T2(1:n) = T(1:n)
          CALL MGmv( A, T2, T )
          omega = MGdot( n,T,S ) / MGdot( n,T,T )

          oldrho = rho
          Mr(1:n) = S(1:n) - omega*T(1:n)
          Mx(1:n) = Mx(1:n) + alpha*T1(1:n) + omega*T2(1:n)

          RNorm = MGnorm( n, Mr ) / MGNorm( n, Mb )
          IF ( RNorm < Conv ) EXIT
       END DO

!      WRITE( Message, * ) 'Iters: ', i, RNorm
!      CALL Info( 'BiCGParEigen', Message, Level=4 ) 


       DEALLOCATE( Ri,P,V,T,T1,T2,S )

10     Continue

       x(1:n) = Mx(1:n)
       b(1:n) = Mb(1:n)
!------------------------------------------------------------------------------
    END SUBROUTINE BiCGParEigen
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
       s = ParallelNorm( n, x )
!------------------------------------------------------------------------------
    END FUNCTION MGnorm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGdot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       REAL(KIND=dp) :: s
       REAL(KIND=dp) CONTIG :: x(:),y(:)
!------------------------------------------------------------------------------
       s = ParallelDot( n, x, y )
!------------------------------------------------------------------------------
    END FUNCTION MGdot
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
END MODULE ParallelEigenSolve

!> \}
