!/*****************************************************************************/
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
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  This module provides various solution techniques for eigenvalue problems. 
! *  The module is based on the ARPACK example driver dndrv3 that has been modified
! *  to fit the needs of ElmerSolver.
! *
! *
! *  FILE: ndrv3.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
! *
! *  Example Authors: Richard Lehoucq, Danny Sorensen and Chao Yang
! *                   Dept. of Computational & Applied Mathematics
! *                   Rice University
! *                   Houston, Texas
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 21 Oct 2000
! *
! *****************************************************************************/


!----------------------------------------------------------------------------
!> Module containing ARPACK routines for the solution of Eigenvalue problems.
!----------------------------------------------------------------------------
!> \ingroup ElmerLib 
!> \{

MODULE EigenSolve

   IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Solution of Eigen value problems using ARPACK library. 
!------------------------------------------------------------------------------
     SUBROUTINE ArpackEigenSolve( Solver,Matrix,N,NEIG,EigValues,EigVectors )
!------------------------------------------------------------------------------
      USE CRSMatrix
      USE IterSolve
      USE Multigrid

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: Matrix, A
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG, DPERM(n)
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      REAL(KIND=dp), TARGET :: WORKD(3*N), RESID(N),bb(N),xx(N)
      REAL(KIND=dp), POINTER CONTIG :: x(:), b(:)
      INTEGER :: IPARAM(11), IPNTR(14)
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      REAL(KIND=dp), ALLOCATABLE :: WORKL(:), D(:,:), WORKEV(:), V(:,:)
      CHARACTER(LEN=MAX_NAME_LEN) :: Method

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      CHARACTER ::     BMAT*1, Which*2, DirectMethod*100
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, k, l, p, IERR, iter, &
                       NCONV, maxitr, ishfts, mode, istat, Dofs
      LOGICAL   ::     First, Stat, Direct = .FALSE., &
                       Iterative = .FALSE., NewSystem, Damped, Stability, &
                       NormalizeToUnity

      LOGICAL :: Factorize, FreeFactorize,FoundFactorize,FoundFreeFactorize
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL, r

      COMPLEX(KIND=dp) :: s
!
      REAL(KIND=dp), POINTER CONTIG :: SaveValues(:), SaveRhs(:)
      TYPE(ValueList_t), POINTER :: Params

!     %--------------------------------------%
!     | Check if system is damped and if so, |
!     | move to other subroutine             |
!     %--------------------------------------%

      Params => Solver % Values
      Damped = ListGetLogical( Params, 'Eigen System Damped', stat )
      IF ( .NOT. stat ) Damped = .FALSE.

      IF ( Damped ) THEN
         CALL ArpackDampedEigenSolve( Solver, Matrix, 2*N, 2*NEIG, &
              EigValues,EigVectors )
         RETURN
      END IF

!     %----------------------------------------%
!     | Check if stability analysis is defined |
!     | and if so move to other subroutine     |
!     %----------------------------------------%

      Stability = ListGetLogical( Params, 'stability analysis', stat )
      IF ( .NOT. stat ) Stability = .FALSE.

      IF ( Stability ) THEN
         CALL ArpackStabEigenSolve( Solver, Matrix, N, NEIG, EigValues,EigVectors )
         RETURN
      END IF

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
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
      NCV = ListGetInteger( Params, 'Eigen System Lanczos Vectors', stat )
      IF ( .NOT. stat ) NCV = 3*NEIG + 1

      IF ( NCV <=  NEIG ) THEN
         CALL Fatal( 'EigenSolve', & 
               'Number of Lanczos vectors must exceed the number of eigenvalues.' )
      END IF

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV,3), &
                WORKEV(3*NCV), V(n,NCV), CHOOSE(NCV), STAT=istat )

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
!

      TOL = ListGetConstReal( Params, 'Eigen System Convergence Tolerance', stat )
      IF ( .NOT. stat ) THEN
         TOL = 100 * ListGetConstReal( Params, 'Linear System Convergence Tolerance' )
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
         SELECT CASE( ListGetString(Params,'Eigen System Select',stat) )
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
         SELECT CASE( ListGetString(Params,'Eigen System Select',stat) )
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

      Maxitr = ListGetInteger( Params, 'Eigen System Max Iterations', stat )
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
      Factorize = ListGetLogical( Params, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Params, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Params, &
                'Linear System Free Fctorization', FoundFreeFactorize )
      CALL ListAddLogical( Params,  &
                     'Linear System Free Factorization',.FALSE. )

      IF ( .NOT. Matrix % Lumped ) THEN
        SigmaR = ListGetConstReal( Params,'Eigen System Shift', stat )
        IF ( SigmaR /= 0.0d0 ) THEN
          Matrix % Values = Matrix % Values - SigmaR * Matrix % MassValues
        END IF

        Method = ListGetString( Params,'Linear System Solver', stat )         
        IF ( Method == 'direct' ) THEN
          DirectMethod = ListGetString( Params, &
              'Linear System Direct Method', stat )
          
          SELECT CASE( DirectMethod )
          CASE('umfpack', 'big umfpack', 'mumps', 'superlu', 'pardiso')
          CASE DEFAULT
            Stat = CRS_ILUT(Matrix, 0.0d0)
          END SELECT
        END IF
      END IF

!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      iter = 1
      NewSystem = .TRUE.

      Iterative = ListGetString( Params, &
               'Linear System Solver', stat ) == 'iterative'

      stat = ListGetLogical( Params,  'No Precondition Recompute', stat  )
      IF ( Iterative .AND. Stat ) THEN
         CALL ListAddLogical( Params, 'No Precondition Recompute', .FALSE. )
      END IF

      DO WHILE( ido /= 99 )

!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         IF ( Matrix % Symmetric ) THEN
            CALL DSAUPD ( ido, BMAT, n, Which, NEIG, TOL, &
              RESID, NCV, V, n, IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )
         ELSE
            CALL DNAUPD ( ido, BMAT, n, Which, NEIG, TOL, &
              RESID, NCV, v, n, IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )
         END IF
 
         IF (ido == -1 .OR. ido == 1) THEN

           WRITE( Message,'(A,I0)') 'Arpack reverse communication calls: ',Iter
           CALL Info( 'EigenSolve', Message, Level=6 )
!            CALL Info( 'EigenSolve', '.', .TRUE., Level=5 )
            iter = iter + 1

!---------------------------------------------------------------------
!             Perform  y <--- OP*x = inv[M]*A*x   (lumped mass)
!                      ido =-1 inv(A-sigmaR*M)*M*x 
!                      ido = 1 inv(A-sigmaR*M)*z
!---------------------------------------------------------------------

            IF ( Matrix % Lumped ) THEN
              CALL CRS_MatrixVectorMultiply( Matrix, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
              DO i=0,n-1
                WORKD( IPNTR(1)+i ) = WORKD( IPNTR(2)+i )
              END DO
              DO i=0,n-1
                WORKD( IPNTR(2)+i ) = WORKD( IPNTR(1)+i ) / &
                    Matrix % MassValues( Matrix % Diag(i+1) )
              END DO
            ELSE              
             
              Dofs = Solver % Variable % Dofs 
              A => Matrix
              x => workd(ipntr(2):ipntr(2)+n-1)

              IF ( ido == -1 ) THEN
                SaveValues => A % Values
                A % Values => A % MassValues
                CALL CRS_MatrixVectorMultiply( A, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
                A % Values => SaveValues
                
                DO i=0,n-1
                  WORKD( IPNTR(1)+i ) = WORKD( IPNTR(2)+i )
                END DO
                b => workd(ipntr(1):ipntr(1)+n-1)
              ELSE
                b => workd(ipntr(3):ipntr(3)+n-1)                              
              END IF

              ! Some strategies (such as 'block') may depend on that these are set properly 
              ! to reflect the linear problem under study.            
              SaveRhs => A % rhs
              A % rhs => b

              SELECT CASE( Method ) 
              CASE('multigrid')
                CALL MultiGridSolve( A, x, b, &
                    DOFs, Solver, Solver % MultiGridLevel, NewSystem )
              CASE('iterative')
                CALL IterSolver( A, x, b, Solver )
              CASE('block')
                CALL BlockSolveExt( A, x, b, Solver )
              CASE ('direct')
                CALL DirectSolver( A, x, b, Solver )
              CASE DEFAULT
                CALL Fatal('EigenSolve','Unknown linear system method: '//TRIM(Method))
              END SELECT

              A % rhs => SaveRhs
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
               SaveValues => Matrix % Values
               Matrix % Values => Matrix % MassValues
               CALL CRS_MatrixVectorMultiply( Matrix, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
               Matrix % Values => SaveValues
            END IF
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

      IF ( .NOT. FoundFreeFactorize ) THEN
        CALL ListRemove( Params, 'Linear System Free Factorization' )
      ELSE
        CALL ListAddLogical( Params, 'Linear System Free Factorization', FreeFactorize )
      END IF

!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
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
      ELSE 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         D = 0.0d0
         IF ( Matrix % Symmetric ) THEN
            CALL DSEUPD ( .TRUE., 'A', Choose, D, V, N, SigmaR,  &
               BMAT, n, Which, NEIG, TOL, RESID, NCV, V, N, &
               IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
         ELSE
            CALL DNEUPD ( .TRUE., 'A', Choose, D, D(1,2), &
               V, N, SigmaR, SigmaI, WORKEV, BMAT, N, &
               Which, NEIG, TOL, RESID, NCV, V, N, &
               IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
         END IF
 
!        %----------------------------------------------%
!        | Eigenvalues are returned in the First column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the First NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         IF (IERR /= 0) THEN 
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
! 
            WRITE( Message, * ) ' Error with DNEUPD, info = ', IERR
            CALL Fatal( 'EigenSolve', Message )
         END IF
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         IF ( kinfo == 1 ) THEN
            CALL Fatal( 'EigenSolve', 'Maximum number of iterations reached.' )
         ELSE IF ( kinfo == 3 ) THEN
            CALL Fatal( 'EigenSolve', 'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
         END IF      
!
!        Sort the eigenvalues to ascending order:
!        ----------------------------------------
         ALLOCATE( Perm(NEIG) )
         Perm = (/ (i, i=1,NEIG) /)
         DO i=1,NEIG
            EigValues(i) = CMPLX( D(i,1), D(i,2),KIND=dp )
         END DO
         CALL SortC( NEIG, EigValues, Perm )
         IF( MINVAL( Perm ) < 1 .OR. MAXVAL( Perm ) > NEIG ) THEN
           CALL Fatal('EigenSolve','Reordering of EigenValues failed')
         END IF

!
!        Extract the values to Elmer structures:
!        -----------------------------------------
         CALL Info( 'EigenSolve', ' ', Level=4 )
         CALL Info( 'EigenSolve', 'Eigen system solution complete: ', Level=4 )
         CALL Info( 'EigenSolve', ' ', Level=4 )
         WRITE( Message,'(A,ES12.3)') 'Convergence criterion is: ', TOL
         CALL Info( 'EigenSolve', Message, Level=4 )
         WRITE( Message,'(A,I0)' ) 'Number of converged Ritz values is: ', IPARAM(5)
         CALL Info( 'EigenSolve', Message, Level=4 )
         WRITE( Message,'(A,I0)') 'Number of update iterations taken: ',  IPARAM(3)
         CALL Info( 'EigenSolve', Message, Level=4 )
         CALL Info( 'EigenSolve', ' ', Level=4 )
         WRITE( Message,'(A,I0,A)') 'Computed ',NEIG,' Eigen Values'
         CALL Info( 'EigenSolve',Message, Level=4 )
         CALL Info( 'EigenSolve', '--------------------------------', Level=4 )

         ! Restore matrix values, if modified when using shift:
         ! ---------------------------------------------------
         IF ( SigmaR /= 0.0d0 ) THEN
           Matrix % Values = Matrix % Values + SigmaR * Matrix % MassValues
         END IF

         ! extract vectors:
         ! ----------------
         k = 1
         DO i=1,NEIG
           p = Perm(i)
           WRITE( Message,'(I0,A,2ES15.6)') i,': ',EigValues(i)
           CALL Info( 'EigenSolve', Message, Level=4 )

           k = 1
           DO j=1,p-1
              IF ( D(j,2) == 0 ) THEN
                 k = k + 1
              ELSE
                 k = k + 2
              END IF
           END DO

           CALL Info('EigenSolve','Copying Eigenvectors to solution',Level=12)

           DO j=1,N
              IF ( D(p,2) /= 0.0d0 ) THEN
                 EigVectors(i,j) = CMPLX( V(j,k),V(j,k+1),KIND=dp )
              ELSE
                 EigVectors(i,j) = CMPLX( V(j,k),0.0d0,KIND=dp )
              END IF
           END DO

           ! Normalization moved to ScaleEigenVectors

         END DO
                    
         IF ( ListGetLogical( Params, 'Eigen System Compute Residuals', stat ) ) THEN
           CALL Info('EigenSolve','Computing eigen system residuals',Level=8)
           CALL CheckResiduals( Matrix, Neig, EigValues, EigVectors )
         END IF
         CALL Info( 'EigenSolve', '--------------------------------',Level=4 )
      END IF

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm )

#else
      CALL Fatal( 'EigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
!
!------------------------------------------------------------------------------
    END SUBROUTINE ArpackEigenSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Scale both real and complex valued eigenvectors. 
!------------------------------------------------------------------------------
    SUBROUTINE ScaleEigenVectors( Matrix, EigVectors, NoEigen, NormalizeToUnity)

      USE CRSMatrix
      USE IterSolve
      USE Multigrid

      IMPLICIT NONE

      TYPE(Matrix_t), TARGET :: Matrix
      COMPLEX(KIND=dp) :: EigVectors(:,:)
      INTEGER :: n, NoEigen
      LOGICAL :: NormalizeToUnity

      INTEGER :: i,j,k,l, mk, mj
      REAL(KIND=dp) :: r
      COMPLEX(KIND=dp) :: s, s1, mx

      CALL Info('ScaleEigenVectors','Scaling eigen vectors',Level=10)

      ! Real case: Normalize eigenvector  (x) so that x^T(M x) = 1
      ! Complex case: Normalize eigenvectors (x) so that x^H(M x) = 1
      ! (probably already done, but no harm in redoing!)
      ! Optionally normalize such that the maximum amplitude is set one.
      ! -----------------------------------------------------------------------------
      
      n = Matrix % NumberOfRows
      IF ( Matrix % Complex ) n = n / 2

      DO i = 1, NoEigen

        IF(  Matrix % COMPLEX ) THEN
          s = 0.0_dp
          IF( NormalizeToUnity ) THEN
            DO j=1,n
              s1 = EigVectors(i,j) * CONJG(EigVectors(i,j))
              IF( ABS( s1 ) > ABS( s ) ) s = s1
            END DO
          ELSE IF ( Matrix % Lumped ) THEN
            DO j=1,n
              s = s + ABS( EigVectors(i,j) )**2 * Matrix % MassValues( Matrix % Diag(2*j-1) )
            END DO
          ELSE
            DO j=1,Matrix % NumberOfRows,2
              DO l=Matrix % Rows(j),Matrix % Rows(j+1)-1,2
                mx = CMPLX( Matrix % MassValues(l), -Matrix % MassValues(l+1), KIND=DP )
                mj  = (j-1)/2 + 1
                mk  = (Matrix % Cols(l)-1)/2 + 1
                s = s + mx * CONJG( EigVectors(i,mj) ) * EigVectors(i,mk)
              END DO
            END DO
          END IF
          
          s = CMPLX( ParallelReduction( REAL(s) ), ParallelReduction( AIMAG(s) ), KIND=dp )
                  
          IF ( ABS(s) > 0 ) THEN
            s = SQRT(s) 
            WRITE(Message,'(A,2ES12.3)') 'Normalizing Eigenvector with: ',REAL(s),AIMAG(s)
            CALL Info('EigenSolve',Message,Level=12)
            EigVectors(i,1:n) = EigVectors(i,1:n) / s
          ELSE
            CALL Warn('EigenSolve','Eigenmode has zero amplitude!')
          END IF
        ELSE          
          r = 0.0_dp
          IF( NormalizeToUnity ) THEN
            DO j=1,n
              r = MAX( r, ABS(EigVectors(i,j))**2 )
            END DO
          ELSE IF ( Matrix % Lumped ) THEN
            DO j=1,n
              r= r + ABS(EigVectors(i,j))**2 * &
                  Matrix % MassValues(Matrix % Diag(j))
            END DO
          ELSE
            r = 0
            DO j=1,n
              DO l=Matrix % Rows(j), Matrix % Rows(j+1)-1
                r = r +  CONJG(EigVectors(i,j)) * Matrix % MassValues(l) * EigVectors(i,Matrix % Cols(l))
              END DO
            END DO
          END IF
          
          r = ParallelReduction(r) 

          IF( ABS(r - 1) < EPSILON( r ) ) THEN
            CALL Info('EigenSolve','Eigenmode already normalized!',Level=12)              
          ELSE IF ( ABS(r) > 0 ) THEN            
            r = SQRT( r ) 
            WRITE(Message,'(A,ES12.3)') 'Normalizing Eigenvector with: ',r
            CALL Info('EigenSolve',Message,Level=12)
            EigVectors(i,:) = EigVectors(i,:) / r
          ELSE
            CALL Warn('EigenSolve','Eigenmode has zero amplitude!')
          END IF
        END IF
          
      END DO

    END SUBROUTINE ScaleEigenVectors
!------------------------------------------------------------------------------

     
!------------------------------------------------------------------------------
    SUBROUTINE CheckResiduals( Matrix, n, Eigs, EigVectors )
!------------------------------------------------------------------------------
      USE CRSMatrix
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER :: i,n,sz
      COMPLEX(KIND=dp) :: Eigs(:), EigVectors(:,:)

      REAL(KIND=dp), ALLOCATABLE :: x(:), y(:)

      sz = Matrix % NumberOfRows
      ALLOCATE( x(sz), y(sz) )
      DO i=1,n
        Matrix % Values = Matrix % Values - REAL(Eigs(i)) * Matrix % MassValues
        x = REAL( EigVectors(i,1:sz) )
        CALL CRS_MatrixVectorMultiply( Matrix, x, y )
        Matrix % Values = Matrix % Values + REAL(Eigs(i)) * Matrix % MassValues

        WRITE( Message, * ) 'L^2 Norm of the residual: ', i, SQRT(SUM(y**2))
        CALL Info( 'CheckResiduals', Message, Level = 3 )
      END DO
      DEALLOCATE( x,y )
!------------------------------------------------------------------------------
END SUBROUTINE CheckResiduals
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solution of Eigen value problems using ARPACK library, stabilized version. 
!------------------------------------------------------------------------------
     SUBROUTINE ArpackStabEigenSolve( Solver, &
          Matrix, N, NEIG, EigValues, EigVectors )
!------------------------------------------------------------------------------
      USE CRSMatrix
      USE IterSolve
      USE Multigrid

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: Matrix
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG, DPERM(n)
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      REAL(KIND=dp), TARGET :: WORKD(3*N), RESID(N)
      REAL(KIND=dp), POINTER CONTIG :: x(:), b(:)
      INTEGER :: IPARAM(11), IPNTR(14)
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      REAL(KIND=dp), ALLOCATABLE :: WORKL(:), D(:,:), V(:,:)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      CHARACTER ::     BMAT*1, Which*2, DirectMethod*200
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, k, l, p, IERR, iter, &
                       NCONV, maxitr, ishfts, mode, istat
      LOGICAL   ::     First, Stat, Direct = .FALSE., FactorizeSetting, &
                       Iterative = .FALSE., NewSystem, &
                       rvec = .TRUE.
      LOGICAL :: Factorize, FreeFactorize,FoundFactorize,FoundFreeFactorize
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL

      COMPLEX(KIND=dp) :: s
!
      REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
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
      IF ( Matrix % Lumped ) THEN
         CALL Error( 'BucklingEigenSolve', &
              'Lumped matrices are not allowed in stability analysis.' )
      END IF

      NCV = 3 * NEIG + 1

      lWORKL = NCV*(NCV+8)

      ALLOCATE( WORKL(lWORKL), D(NCV,2), V(N,NCV), CHOOSE(NCV), STAT=istat )

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'StabEigenSolve', 'Memory allocation error.' )
      END IF

      TOL = ListGetConstReal( Solver % Values, 'Eigen System Convergence Tolerance', stat )
      IF ( .NOT. stat ) THEN
         TOL = 100 * ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance' )
      END IF
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
      IDO   = 0
      kinfo = 0
      ishfts = 1
      BMAT  = 'G'
      Mode = 2

      SELECT CASE( ListGetString( Solver % Values,'Eigen System Select', stat ) )
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
      CASE( 'smallest algebraic' )
         Which = 'LA'
      CASE( 'largest algebraic' )
         Which = 'SA'
      CASE DEFAULT
         Which = 'LM'
      END SELECT

      Maxitr = ListGetInteger( Solver % Values, 'Eigen System Max Iterations', stat )
      IF ( .NOT. stat ) Maxitr = 300

      IPARAM = 0
      IPARAM(1) = ishfts
      IPARAM(3) = maxitr 
      IPARAM(7) = mode

      SigmaR = 0.0d0
      SigmaI = 0.0d0
      V = 0.0d0
      D = 0.0d0

      Factorize = ListGetLogical( Solver % Values, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Solver % Values, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Solver % Values, &
                'Linear System Refactorize', FoundFreeFactorize )
      CALL ListAddLogical( Solver % Values,  &
                     'Linear System Free Factorization',.FALSE. )

      Direct = ListGetString( Solver % Values, &
           'Linear System Solver', stat ) == 'direct'
      IF ( Direct ) THEN
         DirectMethod = ListGetString( Solver % Values, &
           'Linear System Direct Method', stat )

         SELECT CASE( DirectMethod )
         CASE('umfpack', 'big umfpack','mumps', 'superlu', 'pardiso' )
         CASE DEFAULT
            Stat = CRS_ILUT(Matrix, 0.0d0)
         END SELECT
      END IF

      Iterative = ListGetString( Solver % Values, &
               'Linear System Solver', stat ) == 'iterative'

      stat = ListGetLogical( Solver % Values, 'No Precondition Recompute', stat  )

      IF ( Iterative .AND. Stat ) THEN
         CALL ListAddLogical( Solver % Values, 'No Precondition Recompute', .FALSE. )
      END IF
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
      iter = 1
      NewSystem = .TRUE.

      DO WHILE( ido /= 99 )

         CALL DSAUPD( ido, BMAT, n, Which, NEIG, TOL, &
              RESID, NCV, V, n, IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )

         IF( ido==-1 .OR. ido==1 ) THEN
!           WRITE( Message, * ) 'Arpack reverse communication calls: ', Iter
!           CALL Info( 'StabEigenSolve', Message, Level=5 )
            CALL Info( 'StabEigenSolve', '.', .TRUE., Level=5 )
            iter = iter + 1
         END IF

         SELECT CASE( ido )
         CASE( -1, 1 )
            SaveValues => Matrix % Values
            Matrix % Values => Matrix % MassValues
            CALL CRS_MatrixVectorMultiply( Matrix, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
            Matrix % Values => SaveValues
            
            DO i=0,n-1
               WORKD( IPNTR(1)+i ) = WORKD( IPNTR(2)+i )
            END DO
            
            IF ( Direct ) THEN
              CALL DirectSolver( Matrix,WORKD(IPNTR(2)),WORKD(IPNTR(1)), Solver )
            ELSE               
               x => workd(ipntr(2):ipntr(2)+n-1)
               b => workd(ipntr(1):ipntr(1)+n-1)
               
               IF ( Solver % MultiGridSolver ) THEN
                  CALL MultiGridSolve( Matrix, x, b, Solver % Variable % DOFs,  &
                       Solver, Solver % MultiGridLevel, NewSystem )
               ELSE
                  CALL IterSolver( Matrix, x, b, Solver )
               END IF
               
            END IF

         CASE( 2 )            
            CALL CRS_MatrixVectorMultiply( Matrix, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )

         END SELECT

         IF ( NewSystem .AND. ido /= 2 ) THEN
            IF ( Iterative ) THEN
               CALL ListAddLogical( Solver % Values,  'No Precondition Recompute', .TRUE. )
            ELSE
               CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
            END IF
            NewSystem = .FALSE.
         END IF

!-----------------------------------------------------------------------------------------      
      END DO  ! ido == 99

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
!
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      IF ( kinfo /= 0 ) THEN
         WRITE( Message, * ) 'Error with DSAUPD, info = ',kinfo
         CALL Fatal( 'StabEigenSolve', Message )
!
      ELSE 
         D = 0.0d0
         rvec = .TRUE.

         CALL DSEUPD ( rvec, 'A', Choose, D, V, N, SigmaR, &
              BMAT, n, Which, NEIG, TOL, RESID, NCV, V, N, &
              IPARAM, IPNTR, WORKD, WORKL, lWORKL, IERR )
            
!        %----------------------------------------------%
!        | Eigenvalues are returned in the First column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the First NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         IF (IERR /= 0) THEN
            WRITE( Message, * ) ' Error with DSEUPD, info = ', IERR
            CALL Fatal( 'StabEigenSolve', Message )
         END IF
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         IF ( kinfo == 1 ) THEN
            CALL Fatal( 'StabEigenSolve', 'Maximum number of iterations reached.' )
         ELSE IF ( kinfo == 3 ) THEN
            CALL Fatal( 'StabEigenSolve', 'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
         END IF      
!
!        Sort the eigenvalues to ascending order:
!        ----------------------------------------
         ALLOCATE( Perm(NEIG) )
         Perm = (/ (i, i=1,NEIG) /)
         DO i=1,NEIG
            EigValues(i) = CMPLX( 1.0d0 / D(i,1), D(i,2),KIND=dp )
         END DO

         CALL SortC( NEIG, EigValues, Perm )
!
!        Extract the values to ELMER structures:
!        -----------------------------------------
         CALL Info( 'StabEigenSolve', ' ', Level=4 )
         CALL Info( 'StabEigenSolve', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
         CALL Info( 'StabEigenSolve', ' ', Level=4 )
         WRITE( Message, * ) 'The convergence criterion is ', TOL
         CALL Info( 'StabEigenSolve', Message, Level=4 )
         WRITE( Message, * ) ' The number of converged Ritz values is ', IPARAM(5)
         CALL Info( 'StabEigenSolve', Message, Level=4 )
         CALL Info( 'StabEigenSolve', ' ', Level=4 )
         CALL Info( 'StabEigenSolve', 'Computed Eigen Values: ', Level=4 )
         CALL Info( 'StabEigenSolve', '--------------------------------', Level=4 )
         k = 1

         DO i=1,NEIG
           p = Perm(i)
           WRITE( Message, * ) i,EigValues(i)
           CALL Info( 'StabEigenSolve', Message, Level=4 )

           k = 1
           DO j=1,p-1
              IF ( D(j,2) == 0 ) THEN
                 k = k + 1
              ELSE
                 k = k + 2
              END IF
           END DO

           DO j=1,N
              IF ( D(p,2) /= 0.0d0 ) THEN
                 EigVectors(i,j) = CMPLX( V(j,k),V(j,k+1),KIND=dp )
              ELSE
                 EigVectors(i,j) = CMPLX( V(j,k),0.0d0,KIND=dp )
              END IF
           END DO
           
           ! Normalizatin moved to ScaleEigenVectors
        END DO

        CALL Info( 'StabEigenSolve', '--------------------------------',Level=4 )

      END IF

      DO i = 1,Neig
         IF( REAL(EigValues(i)) < 0.0d0 ) THEN
            EigVectors(i,:) = EigVectors(i,:) * CMPLX(0.0d0,1.0d0,KIND=dp)
         END IF
      END DO

      DEALLOCATE( WORKL, D, V, CHOOSE, Perm )
#else
      CALL Fatal( 'StabEigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
!
!------------------------------------------------------------------------------
    END SUBROUTINE ArpackStabEigenSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ArpackEigenSolveComplex( Solver,Matrix,N,NEIG, &
                      EigValues, EigVectors )
!------------------------------------------------------------------------------
!> Solution of Eigen value problems using ARPACK library, complex valued version. 
!------------------------------------------------------------------------------
      USE CRSMatrix
      USE IterSolve
      USE Multigrid

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: Matrix, A
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG, DPERM(n)
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef USE_ARPACK
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      COMPLEX(KIND=dp) :: WORKD(3*N), RESID(N)
      INTEGER :: IPARAM(11), IPNTR(14)
      INTEGER, ALLOCATABLE :: Perm(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: WORKL(:), D(:), WORKEV(:), V(:,:)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      CHARACTER ::     BMAT*1, Which*2
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, k, l, p, IERR, iter, &
                       NCONV, maxitr, ishfts, mode, istat, dofs
      LOGICAL   ::     First, Stat, Direct = .FALSE., FoundFactorize,&
                       Iterative = .FALSE., NewSystem, Factorize, FreeFactorize, FoundFreeFactorize

      CHARACTER(LEN=MAX_NAME_LEN) :: DirectMethod, Method
      COMPLEX(KIND=dp) :: Sigma = 0.0d0, s
      REAL(KIND=dp), TARGET :: x(2*n), b(2*n)
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL, RWORK(N)
!
      REAL(KIND=dp), POINTER CONTIG :: SaveValues(:), SaveRhs(:)

      TYPE(ValueList_t), POINTER :: Params
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
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
      Params => Solver % Values

      NCV   = 3*NEIG+1

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV), &
         WORKEV(3*NCV), V(n,NCV+1), CHOOSE(NCV), STAT=istat )

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'EigenSolveComplex', 'Memory allocation error.' )
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
!
      TOL = ListGetConstReal( Solver % Values, 'Eigen System Convergence Tolerance', stat )
      IF ( .NOT. stat ) THEN
         TOL = 100 * ListGetConstReal( Solver % Values, 'Linear System Convergence Tolerance' )
      END IF

      lWORKL = 3*NCV**2 + 6*NCV 
      IDO   = 0
      kinfo = 0
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
         SELECT CASE(ListGetString( Solver % Values, 'Eigen System Select',stat) )
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
         SELECT CASE(ListGetString( Solver % Values, 'Eigen System Select',stat) )
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
!
      Maxitr = ListGetInteger( Solver % Values, 'Eigen System Max Iterations', stat )
      IF ( .NOT. stat ) Maxitr = 300

      IPARAM = 0
      IPARAM(1) = ishfts
      IPARAM(3) = maxitr 
      IPARAM(7) = mode

      SigmaR = 0
      SigmaI = 0
      V = 0

!     Compute LU-factors for (A-\sigma M) (if consistent mass matrix):
!     ----------------------------------------------------------------
      Factorize = ListGetLogical( Params, &
            'Linear System Refactorize', FoundFactorize )
      CALL ListAddLogical( Params, 'Linear System Refactorize',.TRUE. )

      FreeFactorize = ListGetLogical( Params, &
                'Linear System Refactorize', FoundFreeFactorize )

      CALL ListAddLogical( Params,  &
                     'Linear System Free Factorization',.FALSE. )

      IF ( .NOT. Matrix % Lumped ) THEN
        SigmaR = ListGetConstReal( Params,'Eigen System Shift', stat )
        SigmaI = ListGetConstReal( Params,'Eigen System Shift Im', stat )
        Sigma = CMPLX(SigmaR,SigmaI)

        IF ( Sigma /= 0._dp ) THEN
          Matrix % Values = Matrix % Values - Sigma * Matrix % MassValues
        END IF

        Method = ListGetString( Params,'Linear System Solver', stat )         
        IF ( Method == 'direct' ) THEN
          DirectMethod = ListGetString( Params, &
              'Linear System Direct Method', stat )
          
          SELECT CASE( DirectMethod )
          CASE('umfpack', 'big umfpack', 'mumps', 'superlu', 'pardiso')
          CASE DEFAULT
            Stat = CRS_ComplexILUT(Matrix, 0._dp)
          END SELECT
        END IF
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

      stat = ListGetLogical( Solver % Values,  'No Precondition Recompute', stat  )

      IF ( Iterative .AND. Stat ) THEN
         CALL ListAddLogical( Solver % Values, 'No Precondition Recompute', .FALSE. )
      END IF

      A => Matrix

      DO WHILE( ido /= 99 )
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         CALL ZNAUPD ( ido, BMAT, n, Which, NEIG, TOL, &
           RESID, NCV, v, n, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, kinfo )

         IF (ido == -1 .OR. ido == 1) THEN
!           WRITE( Message, * ) 'Arpack reverse communication calls: ', Iter
!           CALL Info( 'EigenSolveComplex', Message, Level=5 )
            CALL Info( 'EigenSolveComplex', '.', .TRUE., Level=5 )
            Iter = Iter + 1
!---------------------------------------------------------------------
!           Perform  y <--- OP*x = inv[M]*A*x   (lumped mass)
!                    ido =-1 inv(A-sigmaR*M)*M*x 
!                    ido = 1 inv(A-sigmaR*M)*z
!---------------------------------------------------------------------

            ! Some strategies (such as 'block') may depend on that these are set properly 
            ! to reflect the linear problem under study.            
            SaveRhs => A % rhs
            A % rhs => b
            Dofs = Solver % Variable % Dofs

            IF ( ido == -1 ) THEN
              SaveValues => A % Values
              A % Values => A % MassValues
              CALL CRS_ComplexMatrixVectorMultiply( A, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
              A % Values => SaveValues
              
              DO i=0,n-1
                b(2*i+1) = REAL(  WORKD( IPNTR(2)+i ) )
                b(2*i+2) = AIMAG( WORKD( IPNTR(2)+i ) )
              END DO
            ELSE
              DO i=0,n-1
                b(2*i+1) = REAL(  WORKD( IPNTR(3)+i ) )
                b(2*i+2) = AIMAG( WORKD( IPNTR(3)+i ) )
              END DO
            END IF

            DO i=0,n-1
               x(2*i+1) = REAL(  WORKD( IPNTR(2)+i ) )
               x(2*i+2) = AIMAG( WORKD( IPNTR(2)+i ) )
            END DO

            SELECT CASE( Method ) 
            CASE('multigrid')
              CALL MultiGridSolve( A, x, b, &
                  DOFs, Solver, Solver % MultiGridLevel, NewSystem )
            CASE('iterative')
              CALL IterSolver( A, x, b, Solver )
            CASE('block')
              CALL BlockSolveExt( A, x, b, Solver )
            CASE ('direct')
              CALL DirectSolver( A, x, b, Solver )
            CASE DEFAULT
              CALL Fatal('EigenSolve','Unknown linear system method: '//TRIM(Method))
            END SELECT

            A % rhs => SaveRhs

            DO i=0,n-1
              WORKD( IPNTR(2)+i ) = CMPLX( x(2*i+1), x(2*i+2),KIND=dp )
            END DO
!
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
            SaveValues => A % Values
            A % Values => A % MassValues
            CALL CRS_ComplexMatrixVectorMultiply(A, WORKD(IPNTR(1)), WORKD(IPNTR(2)) )
            A % Values => SaveValues
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

       CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .TRUE. )
       CALL ListAddLogical( Solver % Values, &
                           'Linear System Free Factorization', .TRUE. )
!
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
         CALL Fatal( 'EigenSolveComplex', Message )
!
      ELSE 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
!           
         D = 0.0d0
         CALL ZNEUPD ( .TRUE., 'A', Choose, D, V, N, Sigma, WORKEV, BMAT, N, Which, NEIG, &
           TOL, RESID, NCV, V, N, IPARAM, IPNTR, WORKD, WORKL, lWORKL, RWORK, IERR )

!        %----------------------------------------------%
!        | Eigenvalues are returned in the First column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the First NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         IF (IERR /= 0) THEN 
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
! 
            WRITE( Message, * ) ' Error with DNEUPD, info = ', IERR
            CALL Fatal( 'EigenSolveComplex', Message )
         END IF
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         IF ( kinfo == 1 ) THEN
            CALL Fatal( 'EigenSolveComplex', 'Maximum number of iterations reached.' )
         ELSE IF ( kinfo == 3 ) THEN
            CALL Fatal( 'EigenSolveComplex', 'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
         END IF      
!
!        Sort the eigenvalues to ascending order:
!        ----------------------------------------
         ALLOCATE( Perm(NEIG) )
         Perm = (/ (i, i=1,NEIG) /)
         DO i=1,NEIG
            EigValues(i) = D(i)
         END DO
         CALL SortC( NEIG, EigValues, Perm )

!
!        Extract the values to ELMER structures:
!        -----------------------------------------
         CALL Info( 'EigenSolveComplex', ' ', Level=4 )
         CALL Info( 'EigenSolveComplex', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
         CALL Info( 'EigenSolveComplex', ' ', Level=4 )
         WRITE( Message, * ) 'The convergence criterion is ', TOL
         CALL Info( 'EigenSolveComplex', Message, Level=4 )
         WRITE( Message, * ) ' The number of converged Ritz values is ', IPARAM(5)
         CALL Info( 'EigenSolveComplex', Message, Level=4 )
         CALL Info( 'EigenSolveComplex', ' ', Level=4 )
         CALL Info( 'EigenSolveComplex', 'Computed Eigen Values: ', Level=4 )
         CALL Info( 'EigenSolveComplex', '--------------------------------', Level=4 )

         ! Restore matrix values, if modified when using shift:
         ! ---------------------------------------------------
         IF ( Sigma /= 0._dp ) THEN
           Matrix % Values = Matrix % Values + Sigma * Matrix % MassValues
         END IF

         k = 1
         DO i=1,NEIG
            p = Perm(i)
            WRITE( Message, * ) i,EigValues(i)
            CALL Info( 'EigenSolveComplex', Message, Level=4 )

            DO j=1,N
               EigVectors(i,j) = V(j,p)
            END DO

            ! Scaling moved to ScaleEigenVectors
         END DO

         IF ( ListGetLogical( Params, 'Eigen System Compute Residuals', stat ) ) THEN
           CALL Info('EigenSolve','Computing eigen system residuals',Level=8)
           CALL CheckResidualsComplex( Matrix, Neig, EigValues, EigVectors )
         END IF
         CALL Info( 'EigenSolveComplex', '--------------------------------',Level=4 )
!
      END IF

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm )
#else
      CALL Fatal( 'EigenSolveComplex', 'Arpack Eigen System Solver not available.' )
#endif
!
!------------------------------------------------------------------------------
     END SUBROUTINE ArpackEigenSolveComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE CheckResidualsComplex( Matrix, n, Eigs, EigVectors )
!------------------------------------------------------------------------------
      USE CRSMatrix
      TYPE(Matrix_t), POINTER :: Matrix
      INTEGER :: i,j,k,n,sz
      COMPLEX(KIND=dp) :: Eigs(:), EigVectors(:,:)

      REAL(KIND=dp), ALLOCATABLE, TARGET :: vals(:)
      REAL(KIND=dp), POINTER CONTIG :: svals(:)
      COMPLEX(KIND=dp) :: c,m
      COMPLEX(KIND=dp), ALLOCATABLE :: x(:), y(:)

      sz = Matrix % NumberOfRows/2
      ALLOCATE( x(sz), y(sz), vals(size(matrix % values)) ); vals=0
      DO i=1,n
        DO j=1,sz
          DO k=Matrix % Rows(2*j-1), Matrix % Rows(2*j)-1,2
            c = CMPLX(Matrix % Values(k), -Matrix % Values(k+1),KIND=dp)
            m = CMPLX(Matrix % MassValues(k), -Matrix% MassValues(k+1),KIND=dp)
            c = c - eigs(i) * m
            vals(k) = REAL(c)
            vals(k+1) = -AIMAG(c)
          END DO
        END DO

        x = EigVectors(i,:)
        svals => Matrix % Values
        Matrix % Values => vals
        CALL CRS_ComplexMatrixVectorMultiply( Matrix, x, y )
        Matrix % Values => svals

        WRITE( Message, * ) 'L^2 Norm of the residual: ', i, SQRT(SUM(ABS(y)**2))
        CALL Info( 'CheckResiduals', Message, Level = 3 )
      END DO
      DEALLOCATE( x,y )
!------------------------------------------------------------------------------
END SUBROUTINE CheckResidualsComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ArpackDampedEigenSolve( Solver, KMatrix, N, NEIG, EigValues, &
          EigVectors )
!------------------------------------------------------------------------------
!> Solution of Eigen value problems using ARPACK library, damped version. 
!------------------------------------------------------------------------------
      USE CRSMatrix
      USE IterSolve
      USE Multigrid
      USE ElementUtils

      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: KMatrix
      TYPE(Solver_t), TARGET :: Solver
      INTEGER :: N, NEIG, DPERM(n)
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:,:)

#ifdef USE_ARPACK

!     %--------------%
!     | Local Arrays |
!     %--------------%

      TYPE(Matrix_t), POINTER :: MMatrix, BMatrix
      REAL(KIND=dp), TARGET :: WORKD(3*N), RESID(N)
      REAL(KIND=dp), POINTER CONTIG :: x(:), b(:)
      INTEGER :: IPARAM(11), IPNTR(14)
      INTEGER, ALLOCATABLE :: Perm(:), kMap(:)
      LOGICAL, ALLOCATABLE :: Choose(:)
      REAL(KIND=dp), ALLOCATABLE :: WORKL(:), D(:,:), WORKEV(:), V(:,:)
      CHARACTER(LEN=MAX_NAME_LEN) :: str
      COMPLEX(KIND=dp) :: s, EigTemp(NEIG)

!     %---------------%
!     | Local Scalars |
!     %---------------%

      CHARACTER ::     BMAT*1, Which*2
      INTEGER   ::     IDO, NCV, lWORKL, kinfo, i, j, k, l, p, IERR, iter, &
                       NCONV, maxitr, ishfts, mode, istat, DampedMaxIter, ILU
      LOGICAL   ::     First, Stat, NewSystem, UseI = .FALSE.
      REAL(KIND=dp) :: SigmaR, SigmaI, TOL, DampedTOL, IScale

!     %-------------------------------------%
!     | So far only iterative solver        |
!     | and non-lumped matrixes are allowed |
!     %-------------------------------------%

      IF ( KMatrix % Lumped ) THEN
         CALL Error( 'DampedEigenSolve', 'Lumped matrixes are not allowed' )
      END IF

      IF (  ListGetString( Solver % Values, 'Linear System Solver', Stat ) &
           == 'direct' ) THEN
         CALL Error( 'DampedEigenSolve', 'Direct solver is not allowed' )
      END IF

      IF ( Solver % MultiGridSolver ) THEN
         CALL Error( 'DampedEigenSolve', 'MultiGrid solver is not allowed' )
      END IF

      Stat = ListGetLogical( Solver % Values, &
           'No Precondition Recompute', Stat  )

      IF ( Stat ) THEN
         CALL ListAddLogical( Solver % Values, &
              'No Precondition Recompute', .FALSE. )
      END IF


!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

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

      NCV = 3 * NEIG + 1

      ALLOCATE( WORKL(3*NCV**2 + 6*NCV), D(NCV,3), &
         WORKEV(3*NCV), V(n,NCV), CHOOSE(NCV), STAT=istat )

      CHOOSE = .FALSE.
      Workl = 0.0d0
      workev = 0.0d0
      v = 0.0d0
      d = 0.0d0

      IF ( istat /= 0 ) THEN
         CALL Fatal( 'DampedEigenSolve', 'Memory allocation error.' )
      END IF

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

      TOL = ListGetConstReal( Solver % Values, &
           'Eigen System Convergence Tolerance', Stat )

      IF ( .NOT. Stat ) THEN
         TOL = 100 * ListGetConstReal( Solver % Values, &
              'Linear System Convergence Tolerance' )
      END IF

      DampedMaxIter = ListGetInteger( Solver % Values, &
           'Linear System Max Iterations', Stat, 1 )

      IF ( .NOT. Stat ) DampedMaxIter = 100

      DampedTOL = ListGetConstReal( Solver % Values, &
           'Linear System Convergence Tolerance', Stat )

      IF ( .NOT. Stat ) DampedTOL = TOL / 100

      UseI = ListGetLogical( Solver % Values, &
                     'Eigen System Use Identity', Stat )

      IF ( .NOT. Stat ) UseI = .TRUE.
!      IF ( .NOT. Stat ) UseI = .FALSE.   Changed by Antti 2004-02-18

      IDO   = 0
      kinfo = 0
      lWORKL = 3*NCV**2 + 6*NCV 
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of DSAUPD is used     |
!     | (IPARAM(7) = 2).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
      
      ishfts = 1
      BMAT  = 'G'
      Mode  = 3
      
      SELECT CASE( ListGetString(Solver % Values, 'Eigen System Select',Stat) )
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

      Maxitr = ListGetInteger(Solver % Values,'Eigen System Max Iterations',Stat)
      IF ( .NOT. Stat ) Maxitr = 300

      IPARAM = 0
      IPARAM(1) = ishfts
      IPARAM(3) = maxitr 
      IPARAM(7) = mode

      SigmaR = 0.0d0
      SigmaI = 0.0d0
      V = 0.0d0
!------------------------------------------------------------------------------
! Create M and B matrixes
!------------------------------------------------------------------------------
      MMatrix => AllocateMatrix()
      MMatrix % Values => KMatrix % MassValues
      MMatrix % NumberOfRows = KMatrix % NumberOfRows
      MMatrix % Rows => KMatrix % Rows
      MMatrix % Cols => KMatrix % Cols
      MMatrix % Diag => KMatrix % Diag

      IScale = MAXVAL( ABS( MMatrix % Values ) )

      BMatrix => AllocateMatrix()
      BMatrix % NumberOfRows = KMatrix % NumberOfRows
      BMatrix % Rows => KMatrix % Rows
      BMatrix % Cols => KMatrix % Cols
      BMatrix % Diag => KMatrix % Diag
      BMatrix % Values => KMatrix % DampValues
!------------------------------------------------------------------------------
!     ILU Preconditioning
!------------------------------------------------------------------------------
      str = ListGetString( Solver % Values, 'Linear System Preconditioning', Stat )

      IF ( .NOT. Stat ) THEN
         CALL Warn( 'DampedEigenSolve', 'Using ILU0 preconditioning' )
         ILU = 0
      ELSE
         IF ( str == 'none' .OR. str == 'diagonal' .OR. &
              str == 'ilut' .OR. str == 'multigrid' ) THEN

           ILU = 0
           CALL Warn( 'DampedEigenSolve', 'Useing ILU0 preconditioning' )
         ELSE IF ( SEQL(str,'ilu') ) THEN
           ILU = ICHAR(str(4:4)) - ICHAR('0')
           IF ( ILU  < 0 .OR. ILU > 9 ) ILU = 0
         ELSE
           ILU = 0
           CALL Warn( 'DampedEigenSolve','Unknown preconditioner type, using ILU0' )
         END IF
      END IF

      KMatrix % Cholesky = ListGetLogical( Solver % Values,  &
              'Linear System Symmetric ILU', Stat )
      Stat = CRS_IncompleteLU( KMatrix, ILU )
      IF ( .NOT. UseI ) Stat = CRS_IncompleteLU( MMatrix, ILU )

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

      iter = 1
      DO WHILE( ido /= 99 )

!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
         CALL DNAUPD ( ido, BMAT, n, Which, NEIG, TOL, RESID, NCV, v, n, &
                 IPARAM, IPNTR, WORKD, WORKL, lWORKL, kinfo )

         SELECT CASE(ido)
         CASE(-1 )
            !
            ! ido =-1 inv(A)*M*x:
            !----------------------------
            x => workd( ipntr(1) : ipntr(1)+n-1 )
            b => workd( ipntr(2) : ipntr(2)+n-1 )

            CALL EigenMGmv2( n/2, MMatrix, x, b, UseI, IScale )
            DO i=1,n
              x(i) = b(i)
            END DO
            CALL EigenBiCG( n, KMatrix, MMatrix, BMatrix, &
                 b, x, DampedMaxIter, DampedTOL, UseI, IScale )

         CASE( 1 )
            ! 
            ! ido =-1 inv(A)*z:
            !--------------------------
!           WRITE( Message, * ) 'Arpack reverse communication calls: ', Iter
!           CALL Info( 'DampedEigenSolve', Message, Level=5 )
            CALL Info( 'DampedEigenSolve', '.', .TRUE., Level=5 )
            iter = iter + 1

            x => workd( ipntr(2) : ipntr(2)+n-1 )
            b => workd( ipntr(3) : ipntr(3)+n-1 )

            CALL EigenBiCG( n, KMatrix, MMatrix, BMatrix, &
                  x, b, DampedMaxIter, DampedTOL, UseI,IScale )

         CASE( 2 )
!           %-----------------------------------------%
!           |         Perform  y <--- M*x.            |
!           | Need the matrix vector multiplication   |
!           | routine here that takes WORKD(IPNTR(1)) |
!           | as the input and returns the result to  |
!           | WORKD(IPNTR(2)).                        |
!           %-----------------------------------------%
            x => workd( ipntr(1): ipntr(1)+n-1 )
            b => workd( ipntr(2): ipntr(2)+n-1 )
            CALL EigenMGmv2( N/2, MMatrix, x, b, UseI, IScale )

         CASE DEFAULT
         END SELECT

      END DO

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
         CALL Fatal( 'DampedEigenSolve', Message )
      ELSE 
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%

         D = 0.0d0
         CALL DNEUPD ( .TRUE., 'A', Choose, D, D(1,2), V, N, SigmaR, SigmaI, &
         WORKEV, BMAT, N, Which, NEIG, TOL, RESID, NCV, V, N, IPARAM, IPNTR, &
                        WORKD, WORKL, lWORKL, IERR )
!        %----------------------------------------------%
!        | Eigenvalues are returned in the First column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the First NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

         IF ( IERR /= 0 ) THEN 
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
            WRITE( Message, * ) ' Error with DNEUPD, info = ', IERR
            CALL Fatal( 'DampedEigenSolve', Message )
         END IF

!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%

         IF ( kinfo == 1 ) THEN
            CALL Fatal( 'DampedEigenSolve', 'Maximum number of iterations reached.' )
         ELSE IF ( kinfo == 3 ) THEN
            CALL Fatal( 'DampedEigenSolve', &
                 'No shifts could be applied during implicit Arnoldi update, try increasing NCV.' )
         END IF      

!        Sort the eigenvalues to ascending order:
!        ( and keep in mind the corresponding vector )
!        ---------------------------------------------
         DO i = 1, NEIG
            EigTemp(i) = CMPLX( D(i,1), D(i,2), KIND=dp )
         END DO

         ALLOCATE( kMap( NEIG ) )
         kMap(1) = 1
         DO i = 2, NEIG
            IF ( AIMAG( EigTemp(i-1) ) == 0 ) THEN
               kMap(i) = kMap(i-1) + 1
            ELSE IF ( EigTemp(i) == CONJG( EigTemp(i-1) ) ) THEN
               kMap(i) = kMap(i-1)
            ELSE
               kMap(i) = kMap(i-1) + 2
            END IF
         END DO

         ALLOCATE( Perm( NEIG ) )
         Perm = (/ (i, i=1,NEIG) /)
         CALL SortC( NEIG, EigTemp, Perm )
         
!        Extract the values to ELMER structures:
!        -----------------------------------------
         CALL Info( 'DampedEigenSolve', ' ', Level=4 )
         CALL Info( 'DampedEigenSolve', 'EIGEN SYSTEM SOLUTION COMPLETE: ', Level=4 )
         CALL Info( 'DampedEigenSolve', ' ', Level=4 )

         WRITE( Message, * ) 'The convergence criterion is ', TOL
         CALL Info( 'DampedEigenSolve', Message, Level=4 )

         WRITE( Message, * ) ' The number of converged Ritz values is ', &
              IPARAM(5)
         CALL Info( 'DampedEigenSolve', Message, Level=4 )

         CALL Info( 'DampedEigenSolve', ' ', Level=4 )
         CALL Info( 'DampedEigenSolve', 'Computed Eigen Values: ', Level=4 )
         CALL Info( 'DampedEigenSolve', '--------------------------------', Level=4 )

!------------------------------------------------------------------------------
! Extracting the right eigen values and vectors.
!------------------------------------------------------------------------------

! Take the first ones separately
!------------------------------------------------------------------------------
         EigValues(1) = EigTemp(1)         
         
         WRITE( Message, * ) 1,EigValues(1)
         CALL Info( 'DampedEigenSolve', Message, Level=4 )
         
         p = Perm(1)
         k = kMap(p)
         IF( AIMAG( EigValues(1) ) == 0 ) THEN
            DO j = 1, N/2
               EigVectors(1,j) = CMPLX( V(j,k), 0.0d0,KIND=dp )
            END DO
         ELSE
            DO j = 1, N/2
               EigVectors(1,j) = CMPLX( V(j,k), V(j,k+1),KIND=dp )
            END DO
         END IF

! Then take the rest of requested values
!------------------------------------------------------------------------------
         l = 2
         DO i = 2, NEIG/2
            IF ( AIMAG( EigValues(i-1) ) /= 0 .AND. &
                 ABS(AIMAG(EigTemp(l))) == ABS(AIMAG(EigValues(i-1)))) l=l+1
            
            EigValues(i) = EigTemp(l)
            IF ( AIMAG( EigValues(i) ) < 0 ) THEN
               EigValues(i) = CONJG( EigValues(i) )
            END IF
            
            WRITE( Message, * ) i,EigValues(i)
            CALL Info( 'DampedEigenSolve', Message, Level=4 )
            
            p = Perm(l)
            k = kMap(p)
            IF( AIMAG( EigValues(i) ) == 0 ) THEN
               DO j = 1, N/2
                  EigVectors(i,j) = CMPLX( V(j,k), 0.0d0,KIND=dp )
               END DO
            ELSE
               DO j = 1, N/2
                  EigVectors(i,j) = CMPLX( V(j,k), V(j,k+1),KIND=dp )
               END DO
            END IF

            l = l + 1
         END DO
               
         DO i = 1, NEIG/2
            s = 0.0d0
            DO j=1,N/2
               DO k = MMatrix % Rows(j), MMatrix % Rows(j+1)-1
                  s = s + MMatrix % Values(k) * &
                       CONJG( EigVectors(i,j) ) * EigVectors(i,MMatrix % Cols(k))
               END DO
            END DO
            IF ( ABS(s) > 0 ) EigVectors(i,:) = EigVectors(i,:) / SQRT(s)
         END DO

         CALL Warn('DampedEigenSolve','Check that the scaling is not done twice if you call this!')
         
         ! Standard scaling moved to ScaleEigenVectors

         
         CALL Info( 'DampedEigenSolve', '--------------------------------',Level=4 )
      END IF

      DEALLOCATE( WORKL, D, WORKEV, V, CHOOSE, Perm )

      NULLIFY( MMatrix % Rows, MMatrix % Cols, MMatrix % Diag, MMatrix % Values )
      CALL FreeMatrix( MMatrix )

      NULLIFY( BMatrix % Rows, BMatrix % Cols, BMatrix % Diag, BMatrix % Values )
      CALL FreeMatrix( BMatrix )
#else
      CALL Fatal( 'DampedEigenSolve', 'Arpack Eigen System Solver not available.' )
#endif
!------------------------------------------------------------------------------
     END SUBROUTINE ArpackDampedEigenSolve
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE EigenBiCG( n, KMatrix, MMatrix, BMatrix, x, b, Rounds, TOL, UseI, IScale )
!------------------------------------------------------------------------------
      USE CRSMatrix

      TYPE(Matrix_t), POINTER :: KMatrix, MMatrix, BMatrix
      INTEGER :: Rounds
      REAL(KIND=dp) ::  TOL, IScale
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      LOGICAL :: UseI
!------------------------------------------------------------------------------
      INTEGER :: i, n
      REAL(KIND=dp) :: alpha, beta, omega, rho, oldrho, bnorm
      REAL(KIND=dp) :: r(n), Ri(n), P(n), V(n), S(n), &
           T(n), T1(n), T2(n), Tmp(n/2)
!------------------------------------------------------------------------------
      CALL EigenMGmv1( n/2, KMatrix, MMatrix, BMatrix, x, r, UseI, IScale )
      r(1:n) = b(1:n) - r(1:n)

      Ri(1:n) = r(1:n)
      P(1:n) = 0
      V(1:n) = 0
      omega  = 1
      alpha  = 0
      oldrho = 1
      Tmp = 0.0d0

      bnorm = EigenMGdot( n,b,b )

      CALL Info( 'EigenBiCG', '--------------------' )
      CALL Info( 'EigenBiCG', 'Begin BiCG iteration' )
      CALL Info( 'EigenBiCG', '--------------------' )

      DO i=1,Rounds
         rho = EigenMGdot( n, r, Ri )
         
         beta = alpha * rho / ( oldrho * omega )
         P(1:n) = r(1:n) + beta * (P(1:n) - omega*V(1:n))
!------------------------------------------------------------------------------
         Tmp(1:n/2) = P(1:n/2)

         IF ( .NOT. UseI ) THEN
            CALL CRS_LUSolve( n/2, MMatrix, Tmp(1:n/2) )
         ELSE
            Tmp(1:n/2) = Tmp(1:n/2) / IScale
         END IF

         V(n/2+1:n) = Tmp(1:n/2)

         Tmp(1:n/2) = P(n/2+1:n)
         CALL CRS_LUSolve( n/2, KMatrix, Tmp(1:n/2) )
         V(1:n/2) = -1*Tmp(1:n/2)

         T1(1:n) = V(1:n)         
         CALL EigenMGmv1( n/2, KMatrix, MMatrix, BMatrix, T1, V, UseI, IScale )
!------------------------------------------------------------------------------
         alpha = rho / EigenMGdot( n, Ri, V )         
         S(1:n) = r(1:n) - alpha * V(1:n)         
!------------------------------------------------------------------------------
         Tmp(1:n/2) = S(1:n/2)

         IF ( .NOT. UseI ) THEN
            CALL CRS_LUSolve( n/2, MMatrix, Tmp(1:n/2) )
         ELSE
            Tmp(1:n/2) = Tmp(1:n/2) / IScale
         END IF

         T(n/2+1:n) = Tmp(1:n/2)

         Tmp(1:n/2) = S(n/2+1:n)
         CALL CRS_LUSolve( n/2, KMatrix, Tmp(1:n/2) )
         T(1:n/2) = -1*Tmp(1:n/2)

         T2(1:n) = T(1:n)
         CALL EigenMGmv1( n/2, KMatrix, MMatrix, BMatrix, T2, T, UseI, IScale )
!------------------------------------------------------------------------------
         omega = EigenMGdot( n,T,S ) / EigenMGdot( n,T,T )         
         oldrho = rho
         r(1:n) = S(1:n) - omega*T(1:n)
         x(1:n) = x(1:n) + alpha*T1(1:n) + omega*T2(1:n)
!------------------------------------------------------------------------------
         WRITE(*,*) i,EigenMGdot( n,r,r ) / bnorm

         IF ( EigenMGdot( n,r,r ) / bnorm < TOL ) THEN
            CALL EigenMGmv1( n/2, KMatrix, MMatrix, BMatrix, x, r, UseI, IScale )
            r(1:n) = b(1:n) - r(1:n)

            WRITE( Message,* ) 'Correct residual:', EigenMGdot( n,r,r ) / bnorm
            CALL Info( 'EigenBiCG', Message )

            IF ( EigenMGdot( n,r,r ) / bnorm < TOL ) EXIT
         END IF
      END DO

      IF ( EigenMGdot( n,r,r ) / bnorm >= TOL ) THEN
         CALL Fatal( 'EigenBiCG', 'Failed to converge' )
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE EigenBiCG
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION EigenMGdot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
      USE Types
      INTEGER :: n
      REAL(KIND=dp) :: s, x(:), y(:)
      
      s = DOT_PRODUCT( x(1:n), y(1:n) )
!------------------------------------------------------------------------------
    END FUNCTION EigenMGdot
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE EigenMGmv1( n, KMatrix, MMatrix, BMatrix, x, b, UseI, IScale )
!------------------------------------------------------------------------------
      USE CRSMatrix

      INTEGER :: n
      TYPE(Matrix_t), POINTER :: KMatrix, MMatrix, BMatrix
      REAL(KIND=dp) :: IScale
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      LOGICAL :: UseI

      REAL(KIND=dp) :: Tmp(n)

      Tmp = 0.0d0
      b = 0.0d0

      IF ( .NOT. UseI ) THEN
         CALL CRS_MatrixVectorMultiply( MMatrix, x(n+1:2*n), Tmp(1:n) )
         b(1:n) = b(1:n) + Tmp(1:n)
      ELSE
         b(1:n) = x(n+1:2*n) * IScale
      END IF

      CALL CRS_MatrixVectorMultiply( KMatrix, x(1:n), Tmp(1:n) )
      b(n+1:2*n) = b(n+1:2*n) - Tmp(1:n)

      CALL CRS_MatrixVectorMultiply( BMatrix, x(n+1:2*n), Tmp(1:n) )
      b(n+1:2*n) = b(n+1:2*n) - Tmp(1:n)
!------------------------------------------------------------------------------
    END SUBROUTINE EigenMGmv1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE EigenMGmv2( n, MMatrix, x, b, UseI, IScale )
!------------------------------------------------------------------------------
      USE CRSMatrix

      INTEGER :: n
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      REAL(KIND=dp) :: IScale
      TYPE(Matrix_t), POINTER :: MMatrix
      LOGICAL :: UseI

      IF ( .NOT. UseI ) THEN
         CALL CRS_MatrixVectorMultiply( MMatrix, x(1:n), b(1:n) )
      ELSE
         b(1:n) = x(1:n) * IScale
      END IF
      CALL CRS_MatrixVectorMultiply( MMatrix, x(n+1:2*n), b(n+1:2*n) )
!------------------------------------------------------------------------------
    END SUBROUTINE EigenMGmv2
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE EigenSolve
!------------------------------------------------------------------------------

!> \} ElmerLib
