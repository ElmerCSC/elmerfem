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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2001
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!> Module containing the smoothers used in multigrid solvers.
!-----------------------------------------------------------------------------
 
MODULE Smoothers

  USE Types
  USE CRSMatrix
  USE Lists
  USE ParallelUtils  

  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
    FUNCTION MGSmooth( Solver, A, Mesh, x, b, r, Level, DOFs, &
        PreSmooth, LowestSmooth, CF) RESULT(RNorm)
!------------------------------------------------------------------------------
      USE ParallelUtils
      TYPE(Solver_t), POINTER :: Solver
      TYPE(Matrix_t), POINTER :: A
      TYPE(Mesh_t) :: Mesh
      INTEGER :: Level, DOFs
      REAL(KIND=dp), TARGET CONTIG :: x(:),b(:),r(:)
      REAL(KIND=dp) :: RNorm, rphi=5.0_dp
      LOGICAL, OPTIONAL :: PreSmooth, LowestSmooth
      INTEGER, POINTER, OPTIONAL :: CF(:)
!------------------------------------------------------------------------------
      CHARACTER(LEN=MAX_NAME_LEN) :: IterMethod, im
      LOGICAL :: Parallel, Found, Lowest, Pre
      TYPE(Matrix_t), POINTER :: M
      INTEGER :: i, j, k, n, Rounds, InvLevel, me
      INTEGER, POINTER :: Iters(:)
      REAL(KIND=dp), POINTER CONTIG :: Mx(:),Mb(:),Mr(:)
      REAL(KIND=dp) :: Omega, Bnorm, TOL
      REAL(KIND=dp), POINTER :: TmpArray(:,:)
      REAL(KIND=dp), ALLOCATABLE :: Q(:), Z(:), Ri(:), T(:), &
             T1(:), T2(:), S(:), V(:), Pr(:), dx(:),diag(:)
!------------------------------------------------------------------------------
      TYPE( IfLColsT), POINTER :: IfL, IfO
      INTEGER :: row
      TYPE (BasicMatrix_t), POINTER :: CurrIf
!------------------------------------------------------------------------------
      
      SAVE Z, Pr, Q, Ri, T, T1, T2, S, V

      Parallel = ParEnv % PEs > 1


      IF ( .NOT. Parallel ) THEN
        M  => A
        Mx => x
        Mb => b
        Mr => r

        n = A % NumberOfRows
        ALLOCATE(Diag(n))
        Diag = A % Values(A % Diag)
      ELSE
        CALL ParallelUpdateSolve( A,x,r )
        M => ParallelMatrix( A, Mx, Mb, Mr )

        n = M % NumberOfRows
        ALLOCATE(Diag(n))
        Diag = M % Values(M % Diag)

#if 0
! obsolete (at least for now...)
        DO i = 1, ParEnv % PEs
          CurrIf => A % ParMatrix % SplittedMatrix % IfMatrix(i)
          IF ( CurrIf % NumberOfRows == 0 ) CYCLE

          IfL => A % ParMatrix % SplittedMatrix % IfLCols(i)
          IfO => A % ParMatrix % SplittedMatrix % IfORows(i)
          DO j = 1, CurrIf % NumberOfRows
            IF ( Currif % RowOwner(j) /= ParEnv % MyPE ) CYCLE

            row = IfO % IfVec(j)
            DO k = CurrIf % Rows(j), CurrIf % Rows(j+1) - 1
              IF ( IfL % IfVec(k) == row ) &
                Diag(row) = Diag(row) + CurrIf % Values(k)
            END DO
          END DO
        END DO
#endif
      END IF
      
      n = M % NumberOfRows
      InvLevel = 1 + Solver % MultiGridTotal - Level

      Lowest = .FALSE.
      IF( PRESENT( LowestSmooth ) ) Lowest = LowestSmooth

      Pre = .FALSE.
      IF( PRESENT( PreSmooth ) ) Pre = PreSmooth


!      Smoothing iterative method:
!      ---------------------------
      IF( Lowest )THEN
        IterMethod = ListGetString( Solver % Values, 'MG Lowest Smoother', Found )
      ELSE
        Found = .FALSE.
      END IF
      IF(.NOT. Found) THEN
        IterMethod = ListGetString( Solver % Values, 'MG Smoother', Found )
      END IF
      IF ( .NOT. Found ) THEN
        IterMethod = ListGetString( Solver % Values, &
            'Linear System Iterative Method', Found )
      END IF
      IF ( .NOT. Found ) THEN
        IF( DOFs == 1) THEN
          IterMethod = 'sgs'
        ELSE
          IterMethod = 'bsgs'
        END IF
      END IF


      Rounds = 0
      IF(Lowest) THEN
        Rounds = ListGetInteger( Solver % Values,'MG Lowest Smoothing Iterations',Found)
      ELSE IF( Pre ) THEN
        Iters => ListGetIntegerArray( Solver % Values,'MG Pre Smoothing Iterations',Found)
        IF(Found) THEN
          Rounds = Iters(MIN(InvLevel,SIZE(Iters)))
        ELSE        
          Rounds = 1
        END IF
      ELSE
        Iters => ListGetIntegerArray( Solver % Values,'MG Post Smoothing Iterations',Found)
        IF(Found) THEN
          Rounds = Iters(MIN(InvLevel,SIZE(Iters)))
        ELSE        
          Rounds = 1
        END IF
      END IF

      IF( Rounds == 0 ) THEN
        CALL Info('MGSmooth','Zero smoothing rounds given, doing nothing.')
        GOTO 10
      END IF

      IF( IterMethod == 'direct1d' .OR. IterMethod == 'psgs' ) THEN
        IF( .NOT. PRESENT( CF ) ) THEN
          CALL Fatal('MGSmooth','Smoother requires CF clustering info: '//TRIM(IterMethod))
        END IF
      END IF

      TOL = ListGetConstReal( Solver % Values, 'MG Smoother Reduction TOL',Found)

      RNorm = MGnorm( n, Mr ) 
      IF ( Rounds <= 0 ) RETURN

      SELECT CASE( IterMethod )
      CASE( 'cg' )
        ALLOCATE( Z(n), Pr(n), Q(n) )
        
      CASE( 'bicgstab' )
        ALLOCATE( Pr(n), Ri(n), T(n), T1(n), T2(n), S(n), V(n) )

      CASE( 'direct1d' ) 
        ALLOCATE( dx(n) )

      END SELECT

      TmpArray => ListGetConstRealArray(Solver % Values,'MG Smoother Relaxation Factor',Found)
      IF( ASSOCIATED(TmpArray)) THEN
        Omega = TmpArray(MIN(InvLevel,SIZE(TmpArray,1)),1)
      ELSE
        Omega = 1.0d0
      END IF

      IF( Pre ) THEN
        CALL Info('MGSmooth','Applying pre-smoother: '//TRIM(IterMethod), Level=10 )
      ELSE
        CALL Info('MGSmooth','Applying post-smoother: '//TRIM(IterMethod), Level=10 )
      END IF


      SELECT CASE( IterMethod )
      CASE( 'jacobi' ) 
        CALL Jacobi( n, A, M, Mx, Mb, Mr, Rounds )
        
      CASE( 'gs' )                         
        CALL GS( n, A, M, Mx, Mb, Mr, Rounds )

      CASE( 'bgs' )                         
        CALL BGS( n, A, M, Mx, Mb, Mr, DOFs,Rounds )
       
      CASE( 'sgs' )                                     
        CALL SGS( n, A, M, Mx, Mb, Mr, Rounds)

      CASE( 'isgs' )                                     
        CALL InternalSGS( n, A, M, x, b, r, Rounds)

      CASE( 'jacobi+isgs' )                                             
        CALL SmoothedJacobi( n, A, M, Mx, Mb, Mr, Omega, Rounds )
        IF(Parallel) CALL ParallelUpdateResult(A,x,r)

        CALL InternalSGS( n, A, M, x, b, r, Rounds)
        IF(Parallel) CALL ParallelUpdateSolve(A,x,r)

        CALL SmoothedJacobi( n, A, M, Mx, Mb, Mr, Omega, Rounds )

      CASE( 'bsgs' )                                     
        CALL BSGS( n, A, M, Mx, Mb, Mr, DOFs, Rounds)
       
      CASE( 'wjacobi' )                                     
        CALL SmoothedJacobi( n, A, M, Mx, Mb, Mr, Omega, Rounds )
        
      CASE( 'wgs' )                                   
        CALL SmoothedGS( n, A, M, Mx, Mb, Mr, Omega, Rounds )
        
      CASE( 'wsgs' )                                     
        CALL SmoothedSGS( n, A, M, Mx, Mb, Mr, Omega, Rounds)
        
      CASE( 'csgs' )                                     
        CALL CSGS( n, A, M, Mx, Mb, Mr, Rounds)
        
      CASE( 'cjacobi' )                                     
        CALL CJacobi( n, A, M, Mx, Mb, Mr, Rounds )
        
      CASE( 'psgs' )                                     
        CALL PostSGS( n, A, M, Mx, Mb, Mr, CF, Rounds)

      CASE( 'direct1d' )                                     
        CALL Direct1dSmoother( n, A, M, Mx, Mb, Mr, CF, Rounds)

      CASE( 'cg' )
        CALL CG( n, A, M, Mx, Mb, Mr, Rounds )

      CASE( 'ccg' )
        CALL CCG( n, A, M, Mx, Mb, Mr, Rounds )
       
      CASE( 'bicgstab' )
       CALL BiCG( n, A, M, Mx, Mb, Mr, Rounds )

      CASE( 'uzawa' )                                   
        CALL Uzawa( n, A, M, Mx, Mb, Mr, Rounds )

      CASE( 'vanka' )                                   
        CALL Vanka( n, A, M, Mx, Mb, Mr, Rounds )

      CASE( 'test gs' )                                   
        CALL TestGS( n, A, M, Mx, Mb, Mr, Rounds )

      CASE DEFAULT
        CALL Warn('MGSmooth','Unknown smoother - '//TRIM(IterMethod)//' using Jacobi')
        CALL Jacobi( n, A, M, Mx, Mb, Mr, Rounds )
      END SELECT


      SELECT CASE( Itermethod )
      CASE( 'cg' )
        DEALLOCATE( Z, Pr, Q)
        
      CASE( 'bicgstab' )
        DEALLOCATE( Pr, Ri, T, T1, T2, S, V )
      END SELECT

10    CONTINUE

      CALL MGmv( A, x, r, .TRUE. )
      r = b - r
      IF ( Parallel ) THEN
        DO i=1,SIZE(Mr)
          Mr(i) = Mb(i) - Mr(i)
        END DO
      END IF
      RNorm = MGnorm( n, Mr ) 

      CALL Info('MGSmooth','Smoothing finished',Level=12)

!------------------------------------------------------------------------------

    CONTAINS 

!------------------------------------------------------------------------------
      FUNCTION MGnorm( n, x ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp)  :: s
        REAL(KIND=dp) CONTIG :: x(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = SQRT( DOT_PRODUCT( x(1:n), x(1:n) ) )
        ELSE
          s = ParallelNorm( n, x )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION MGnorm
!------------------------------------------------------------------------------

      
!------------------------------------------------------------------------------
      FUNCTION MGdot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
        INTEGER :: n
        REAL(KIND=dp)  :: s
        REAL(KIND=dp) CONTIG :: x(:),y(:)
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = ParallelDot( n, x, y )
        END IF
!------------------------------------------------------------------------------
      END FUNCTION MGdot
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION MGCdot( n, x, y ) RESULT(s)
!------------------------------------------------------------------------------
       INTEGER :: n
       COMPLEX(KIND=dp)  :: s
       COMPLEX(KIND=dp) CONTIG :: x(:),y(:)
!------------------------------------------------------------------------------
       s = DOT_PRODUCT( x(1:n), y(1:n) )
!------------------------------------------------------------------------------
    END FUNCTION MGCdot
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE MGmv( A, x, b, Update )
!------------------------------------------------------------------------------
        REAL(KIND=dp) CONTIG :: x(:), b(:)
        TYPE(Matrix_t), POINTER :: A
        LOGICAL, OPTIONAL :: Update
!------------------------------------------------------------------------------
        IF ( .NOT. Parallel ) THEN
          CALL CRS_MatrixVectorMultiply( A, x, b )
        ELSE
          IF ( PRESENT( Update ) ) THEN
            CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
          ELSE
            CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
          END IF
        END IF
!------------------------------------------------------------------------------
      END SUBROUTINE MGmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE MGCmv( A, x, b, Update )
!------------------------------------------------------------------------------
        COMPLEX(KIND=dp) CONTIG :: x(:), b(:)
        LOGICAL, OPTIONAL :: Update
        TYPE(Matrix_t), POINTER :: A
!------------------------------------------------------------------------------
        CALL CRS_ComplexMatrixVectorMultiply( A, x, b )
!------------------------------------------------------------------------------
      END SUBROUTINE MGCmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE Jacobi( n, A, M, x, b, r, Rounds)
!-------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,n
!------------------------------------------------------------------------------
        DO i=1,Rounds
          CALL MGmv(A, x, r)
          DO j=1,n
            r(j) = b(j) - r(j)
            x(j) = x(j) + r(j) / Diag(j)
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE Jacobi
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE SmoothedJacobi( n, A, M, x, b, r, w, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) :: w
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,n
!------------------------------------------------------------------------------
        DO i=1,Rounds
          CALL MGmv( A, x, r )
          DO j=1,n
            r(j) = b(j) - r(j)
            x(j) = x(j) + w * r(j) / Diag(j)
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE SmoothedJacobi
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE CJacobi( n, A, M, rx, rb, rr, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: n,Rounds
        REAL(KIND=dp) CONTIG :: rx(:),rb(:),rr(:)
!------------------------------------------------------------------------------
        COMPLEX(KIND=dp) :: x(n/2),b(n/2),r(n/2)
        INTEGER :: i,j,diag
!------------------------------------------------------------------------------

        DO i=1,n/2
          r(i) = CMPLX( rr(2*i-1), rr(2*i),KIND=dp )
          x(i) = CMPLX( rx(2*i-1), rx(2*i),KIND=dp )
          b(i) = CMPLX( rb(2*i-1), rb(2*i),KIND=dp )
        END DO
        
        DO j=1,Rounds
          CALL MGCmv( A, x, r )
          r(1:n/2) = b(1:n/2) - r(1:n/2)
          
          DO i=1,n/2
            diag = M % diag(2*i-1)
            r(i) = r(i) / CMPLX( M % Values(diag), M % Values(diag+1),KIND=dp)
            x(i) = x(i) + r(i)
          END DO
        END DO
        
        DO i=1,n/2
          rr(2*i-1) =  REAL( r(i) )
          rr(2*i-0) =  AIMAG( r(i) )
          rx(2*i-1) =  REAL( x(i) )
          rx(2*i-0) =  AIMAG( x(i) )
          rb(2*i-1) =  REAL( b(i) )
          rb(2*i-0) =  AIMAG( b(i) )
        END DO

!------------------------------------------------------------------------------
      END SUBROUTINE CJacobi
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE GS( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG  :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------
     
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
                
        DO k=1,Rounds
          
          DO i=1,n
            s = 0.0_dp
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE GS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Richards iteration with preconditioning by lumped mass matrix.
!------------------------------------------------------------------------------
      SUBROUTINE Richards( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG  :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: t, s, q
        INTEGER, POINTER :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------
     
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
                
        DO k=1,Rounds
          
          DO i=1,n
            s = 0.0_dp
            q = 0.0_dp
            DO j=Rows(i),Rows(i+1)-1
              t = Values(j)
              s = s + t * x(Cols(j)) 
              q = q + t
            END DO
            
            r(i) = (b(i)-s) / q
            x(i) = x(i) + r(i)
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE Richards
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE BGS( n, A, M, x, b, r, DOFs, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: DOFs, Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,n,id,dof
        REAL(KIND=dp) :: s(DOFs)
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
!------------------------------------------------------------------------------
     
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
                
        DO k=1,Rounds          
          DO i=1,n / DOFs
            s = 0.0d0
            DO dof=1,DOFs
              id = (i-1)*DOFs + dof
              DO j=Rows(id),Rows(id+1)-1
                s(dof) = s(dof) + x(Cols(j)) * Values(j)
              END DO
            END DO
            DO dof=1,DOFs
              id = (i-1)*DOFs + dof             
              r(id) = (b(id)-s(dof)) / Diag(j)
              x(id) = x(id) + r(id)
            END DO
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE BGS
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
      SUBROUTINE SmoothedGS( n, A, M, x, b, r, w, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) :: w
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
!------------------------------------------------------------------------------
     
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values

       
        DO k=1,Rounds
          DO i=1,n
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + w * r(i)
          END DO
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE SmoothedGS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE SGS( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
        
        DO k=1,Rounds
          DO i=1,n
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
          
          DO i=n,1,-1
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
        END DO
      END SUBROUTINE SGS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Internal symmetric-gauss-seidel for parallel computations
!------------------------------------------------------------------------------
      SUBROUTINE InternalSGS( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        
        Rows   => A % Rows
        Cols   => A % Cols 
        Values => A % Values
        
        DO k=1,Rounds
          DO i=1,A % NumberOFRows
            ! Skip the interface elements as the gauss-seidel cannot be used to update them
            IF( Parallel ) THEN
              IF( A % ParallelInfo % Interface(i) ) CYCLE
            END IF

            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / A % Values(A % Diag(i))
            x(i) = x(i) + r(i)
          END DO
          
          DO i=A % NumberOfRows,1,-1
            IF( Parallel ) THEN
              IF( A % ParallelInfo % Interface(i) ) CYCLE
            END IF

            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / (A % Values(A % Diag(i)))
            x(i) = x(i) + r(i)
          END DO
        END DO
      END SUBROUTINE InternalSGS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Block Symmetric Gauss Seidel 
!------------------------------------------------------------------------------
      SUBROUTINE BSGS( n, A, M, x, b, r, DOFs, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: DOFs, Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,n,id,dof
        REAL(KIND=dp) :: s(DOFs)
        INTEGER, POINTER CONTIG  :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values

        DO k=1,Rounds
          DO i=1,n/DOFs
            s = 0.0d0
            DO dof = 1,DOFs
              id = (i-1)*DOFs + dof
              DO j=Rows(id),Rows(id+1)-1
                s(dof) = s(dof) + x(Cols(j)) * Values(j)
              END DO
            END DO
            DO dof = 1,DOFs
              id = (i-1)*DOFs + dof
              r(id) = (b(id)-s(dof)) / Diag(id)
              x(id) = x(id) + r(id)
            END DO
          END DO
          
          DO i=n/DOFs,1,-1
            s = 0.0d0
            DO dof = 1,DOFs
              id = (i-1)*DOFs + dof
              DO j=Rows(id),Rows(id+1)-1
                s(dof) = s(dof) + x(Cols(j)) * Values(j)
              END DO
            END DO
            DO dof = 1,DOFs
              id = (i-1)*DOFs + dof
              r(id) = (b(id)-s(dof)) / Diag(id)
              x(id) = x(id) + r(id)
            END DO
          END DO
        END DO
      END SUBROUTINE BSGS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE SmoothedSGS( n, A, M, x, b, r, w, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) :: w
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
        
        DO k=1,Rounds
          DO i=1,n
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + w * r(i)
          END DO
          
          DO i=n,1,-1
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + w * r(i)
          END DO
        END DO
      END SUBROUTINE SmoothedSGS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE CSGS( n, A, M, rx, rb, rr, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: rx(:),rb(:),rr(:)
        INTEGER :: i,j,k,n,l
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        COMPLEX(KIND=dp) :: r(n/2),b(n/2),x(n/2),s
        
        DO i=1,n/2
          r(i) = CMPLX( rr(2*i-1), rr(2*i), KIND=dp )
          x(i) = CMPLX( rx(2*i-1), rx(2*i), KIND=dp )
          b(i) = CMPLX( rb(2*i-1), rb(2*i), KIND=dp )
        END DO
        
        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values
        
        DO k=1,Rounds
          DO i=1,n/2
            s = 0.0d0
            
            DO j=Rows(2*i-1),Rows(2*i)-1,2             
              s = s + x((Cols(j)+1)/2) * CMPLX( Values(j), -Values(j+1),KIND=dp)
            END DO
            
            j = A % Diag(2*i-1)
            r(i) = (b(i)-s) / CMPLX( Values(j), -Values(j+1),KIND=dp )
            x(i) = x(i) + r(i)
          END DO
          
          DO i=n/2,1,-1
            s = 0.0d0
            
            DO j=Rows(2*i-1),Rows(2*i)-1,2             
              s = s + x((Cols(j)+1)/2) * CMPLX( Values(j), -Values(j+1),KIND=dp)
            END DO
            
            j = A % Diag(2*i-1)
            r(i) = (b(i)-s) / CMPLX( Values(j), -Values(j+1),KIND=dp )
            x(i) = x(i) + r(i)
          END DO
          
        END DO
        
        DO i=1,n/2
          rr(2*i-1) =  REAL( r(i) )
          rr(2*i-0) =  AIMAG( r(i) )
          rx(2*i-1) =  REAL( x(i) )
          rx(2*i-0) =  AIMAG( x(i) )
          rb(2*i-1) =  REAL( b(i) )
          rb(2*i-0) =  AIMAG( b(i) )
        END DO
        
      END SUBROUTINE CSGS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE PostSGS( n, A, M, x, b, r, f, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        INTEGER, POINTER :: f(:)
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,n
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
        
        n = M % NumberOfRows
        Rows   => M % Rows
        Cols   => M % Cols
        Values => M % Values
        
        DO k=1,Rounds
          
          DO i=1,n
            IF(f(i) /= 0) CYCLE
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
          DO i=1,n
            IF(f(i) == 0) CYCLE
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
          
          DO i=n,1,-1
            IF(f(i) /= 0) CYCLE
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
          DO i=n,1,-1
            IF(f(i) == 0) CYCLE
            s = 0.0d0
            DO j=Rows(i),Rows(i+1)-1
              s = s + x(Cols(j)) * Values(j)
            END DO
            r(i) = (b(i)-s) / Diag(i)
            x(i) = x(i) + r(i)
          END DO
          
        END DO
      END SUBROUTINE PostSGS
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> For some cases a smoother that only works with the local strong connections
!> might be ideal. Given the clustering "f" just picks the entries i and j 
!> such that f(i)=f(j) and use the reduced matrix in a direct solver to smooth
!> the system. Note that the method guarantees that the linear system actually
!> consists of a number of local problems that are fairly small in size. 
!------------------------------------------------------------------------------
      SUBROUTINE Direct1dSmoother( n, A, M, x, b, r, f, Rounds )
!------------------------------------------------------------------------------
        USE DirectSolve
        USE MeshUtils
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        INTEGER, POINTER :: f(:)
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
        INTEGER :: i,j,k,kb,n
        REAL(KIND=dp) :: s,rowsum,frac
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:),ColsB(:),RowsB(:),DiagB(:),&
            TopPointer(:),BotPointer(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:),ValuesB(:)
        TYPE(Matrix_t), POINTER :: Acluster => NULL()
        INTEGER, POINTER :: NodeLayer(:)
        INTEGER :: NoLayers,jc,kc,layer0,dlayer,klayer(-20:20),iter,NoBlocks,&
            ii,jj,mi,mj


        TYPE(Mesh_t), POINTER :: Mesh        

        SAVE :: Acluster, NodeLayer, NoLayers

        n = A % NumberOfRows
        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values

        Mesh => CurrentModel % Mesh
        NoBlocks = Dofs


        ! First time, compute the size of the reduced matrix
        !---------------------------------------------------
        IF(.NOT. ASSOCIATED( Acluster ) ) THEN
          ALLOCATE( Acluster ) 
          ALLOCATE( Acluster % Rows(n+1) )
          Acluster % NumberOfRows = n
          Acluster % Rows = 0
          RowsB => Acluster % Rows

          kb = 1
          RowsB(kb) = 1
          DO i=1,n
            IF( NoBlocks > 1 ) THEN
              ii = ( i - 1) / NoBlocks + 1
            ELSE
              ii = i
            END IF
            DO k=Rows(i),Rows(i+1)-1
              j = Cols(k)
              IF( NoBlocks > 1 ) THEN
                jj = ( j - 1) / NoBlocks + 1
              ELSE
                jj = j
              END IF
              IF( f(ii) == f(jj) ) THEN
                kb = kb + 1
              END IF
            END DO
            RowsB(i+1) = kb
          END DO
          
          ALLOCATE( Acluster % Cols(kb-1), Acluster % Values(kb-1), &
              Acluster % Diag(n) ) 
          
          WRITE(Message,'(A,F8.3,A)') '1D matrix size fraction: ',100.0_dp*kb/SIZE(Values),' %'      
          CALL Info('Direct1dSmoother',Message)

          ! PRINT *,'making matrix structure'          
          Acluster % Cols = 0
          ColsB => Acluster % Cols
          DiagB => Acluster % Diag
          ValuesB => Acluster % Values
          ValuesB = 0.0_dp

          kb = 1
          DO i=1,n
            IF( NoBlocks > 1 ) THEN
              ii = ( i - 1) / NoBlocks + 1
            ELSE
              ii = i
            END IF            
            DO k=Rows(i),Rows(i+1)-1
              j = Cols(k)
              IF( NoBlocks > 1 ) THEN
                jj = ( j - 1) / NoBlocks + 1
              ELSE
                jj = j
              END IF
              IF( i == j ) DiagB(i) = kb
              IF( f(ii) == f(jj) ) THEN
                ValuesB(kb) = ValuesB(kb) + Values(k)
                ColsB(kb) = j
                kb = kb + 1                
              END IF
            END DO
          END DO

          ! PRINT *,'finding top and bottom nodes'
          CALL DetectExtrudedStructure( Mesh, Solver, &
              NumberOfLayers = NoLayers, NodeLayer = NodeLayer )
        END IF


        ! Note: this assumes that the matrix is the same!!
        ! Make better if you start using it more
        !--------------------------------------------------
        ValuesB => Acluster % Values
        RowsB => Acluster % Rows
        ColsB => Acluster % Cols
        DiagB => Acluster % Diag
        ValuesB = 0.0_dp

        ! Now pick up the values for the reduced matrix.
        ! Initialization of ColsB could actually be done 
        ! as a preprocessing step...
        !--------------------------------------------------
        ValuesB = 0.0_dp
        kb = 1        
        klayer = 0

        DO i=1,n
          IF( NoBlocks > 1 ) THEN
            ii = ( i - 1) / NoBlocks + 1
            mi = MODULO( i-1, NoBlocks )
            layer0 = NoBlocks * NodeLayer(ii) + mi 
          ELSE
            ii = i
            layer0 = NodeLayer(ii)
          END IF
                    
          ! First find the dlayer->klayer mapping for the cluster
          !------------------------------------------------------
          DO k=Rows(i),Rows(i+1)-1
            j = Cols(k)

            IF( NoBlocks > 1 ) THEN
              jj = ( j - 1) / NoBlocks + 1
              mj = MODULO( j-1, NoBlocks )
            ELSE
              jj = j
            END IF
            
            IF( f(ii) == f(jj) ) THEN
              IF( NoBlocks > 1 ) THEN
                dlayer = NoBlocks * NodeLayer(jj) +  mj - layer0
              ELSE
                dlayer = NodeLayer(j) - layer0
              END IF
              IF( ABS( dlayer ) > 20 ) THEN
	        PRINT *,'dlayer',dlayer
                CALL Fatal('Direct1dSmoother','Offset in indeces too big!')
              END IF
              klayer( dlayer ) = kb
              kb = kb + 1
            END IF
          END DO
          
          ! Using the mapping map values in all columns to the cluster one
          !----------------------------------------------------------------
          DO k=Rows(i),Rows(i+1)-1
            j = Cols(k)  
            
            IF( NoBlocks > 1 ) THEN
              jj = ( j - 1) / NoBlocks + 1
              mj = MODULO( j-1, NoBlocks )
              dlayer = NoBlocks * NodeLayer(jj) + mj - layer0
            ELSE
              jj = j
              dlayer = NodeLayer(j) - layer0
            END IF

            IF( ABS( dlayer ) > 20 ) THEN
	      PRINT *,'dlayer:',dlayer
              CALL Fatal('Direct1dSmoother','Offset in indeces too big2!')
            END IF
            kc = klayer( dlayer )
            ValuesB(kc) = ValuesB(kc) + Values(k)
          END DO
        END DO
        
        ! Perform given number of rounds
        ! For this smoother one is probably a good value most often
        !------------------------------------------------------------
        DO iter=1,Rounds
          ! PRINT *,'X init range:',MINVAL(x),MAXVAL(x),SUM(x)/SIZE(x)

          CALL MGmv( A, x, r )
          r(1:n) = b(1:n) - r(1:n)
          
          ! Make the correction that is caused when the matrix values are lumped on the 
          ! vertical lines. The objective is that if x is solution of Ax=b then it will
          ! not be modified by this smoother. 
          !----------------------------------------------------------------------------
          klayer = 0
          DO i=1,n

            IF( NoBlocks > 1 ) THEN
              ii = ( i - 1) / NoBlocks + 1
              mi = MODULO( i-1, NoBlocks )
              layer0 = NoBlocks * NodeLayer(ii) + mi 
            ELSE
              ii = i
              layer0 = NodeLayer(i)
            END IF

            ! Get the numbers pointing to different layers [-1,0,1]
            DO k=Rows(i),Rows(i+1)-1
              j = Cols(k)

              IF( NoBlocks > 1 ) THEN
                jj = ( j - 1) / NoBlocks + 1
                mj = MODULO( j-1, NoBlocks )
              ELSE
                jj = j
              END IF
              
              IF( f(ii) == f(jj) ) THEN
                IF( NoBlocks > 1 ) THEN
                  dlayer = NoBlocks * NodeLayer(jj) + mj - layer0
                ELSE
                  dlayer = NodeLayer(j) - layer0
                END IF
                klayer( dlayer ) = j
              END IF
            END DO

            ! For the non-cluster entries perform the lumping to the 1d cluster
            DO k=Rows(i),Rows(i+1)-1
              j = Cols(k)            

              IF( NoBlocks > 1 ) THEN
                jj = ( j - 1) / NoBlocks + 1
                mj = MODULO( j-1, NoBlocks )
              ELSE
                jj = j
              END IF
               
              IF( f(ii) /= f(jj) ) THEN
                IF( NoBlocks > 1 ) THEN
                  dlayer = NoBlocks * NodeLayer(jj) + mj - layer0
                ELSE
                  dlayer = NodeLayer(j) - layer0
                END IF
                jc = klayer( dlayer )
                r(i) = r(i) + Values(k) * (x(j) - x(jc)) 
              END IF
            END DO
          END DO

          ! Solve the system using direct solver (e.g. umfpack) and update the solution
          CALL DirectSolver( Acluster, dx, r, Solver )
          ! This would be full direct solve, and of course does not make sense except for testing
          ! CALL DirectSolver( A, dx, r, Solver )

          x(1:n) = x(1:n) + dx(1:n)
        END DO
        ! PRINT *,'X fin range:',MINVAL(x),MAXVAL(x),SUM(x)/SIZE(x)

      END SUBROUTINE Direct1dSmoother
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Conjugate gradient as a smoother. 
!------------------------------------------------------------------------------
      SUBROUTINE CG( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A,M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,n
        REAL(KIND=dp) :: alpha,rho,oldrho,residual0,residual
!------------------------------------------------------------------------------
        CALL MGmv( A, x, r )
        r(1:n) = b(1:n) - r(1:n)
        residual0 = MGnorm(n, r)
        
        DO i=1,Rounds
          Z(1:n) = r(1:n)
          CALL CRS_LUSolve( n, M, Z )
          rho = MGdot( n, r, Z )
          
          IF ( i == 1 ) THEN
            Pr(1:n) = Z(1:n)
          ELSE
            Pr(1:n) = Z(1:n) + rho * Pr(1:n) / oldrho
          END IF
          
          CALL MGmv( A, Pr, Q )
          alpha  = rho / MGdot( n, Pr, Q )
          oldrho = rho
          
          x(1:n) = x(1:n) + alpha * Pr(1:n)
          r(1:n) = r(1:n) - alpha * Q(1:n)

          residual = MGnorm(n, r)
          IF(residual<TOL*residual0) EXIT
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE CG
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Complex valued conjugate gradient as a smoother.
!------------------------------------------------------------------------------
      SUBROUTINE CCG( n, A, M, rx, rb, rr, Rounds )
!------------------------------------------------------------------------------
        INTEGER :: i,n, Rounds
        TYPE(Matrix_t), POINTER :: A,M
        REAL(KIND=dp) CONTIG :: rx(:),rb(:),rr(:)
        COMPLEX(KIND=dp) :: alpha,rho,oldrho
        COMPLEX(KIND=dp) :: r(n/2),b(n/2),x(n/2)
        COMPLEX(KIND=dp) :: Z(n), Pc(n), Q(n)
!------------------------------------------------------------------------------
        DO i=1,n/2
          r(i) = CMPLX( rr(2*i-1), rr(2*i),KIND=dp )
          x(i) = CMPLX( rx(2*i-1), rx(2*i),KIND=dp )
          b(i) = CMPLX( rb(2*i-1), rb(2*i),KIND=dp )
        END DO
        
        CALL MGCmv( A, x, r )
        r(1:n/2) = b(1:n/2) - r(1:n/2)
        
        DO i=1,Rounds
          Z(1:n/2) = r(1:n/2)
          CALL CRS_ComplexLUSolve( n, M, Z )
          rho = MGCdot( n/2, r, Z )
          
          IF ( i == 1 ) THEN
            Pc(1:n/2) = Z(1:n/2)
          ELSE
            Pc(1:n/2) = Z(1:n/2) + rho * Pc(1:n/2) / oldrho
          END IF
          
          CALL MGCmv( A, Pc, Q )
          alpha  = rho / MGCdot( n/2, Pc, Q )
          oldrho = rho
          
          x(1:n/2) = x(1:n/2) + alpha * Pc(1:n/2)
          r(1:n/2) = r(1:n/2) - alpha * Q(1:n/2)
        END DO

        DO i=1,n/2
          rr(2*i-1) =  REAL( r(i) )
          rr(2*i-0) =  AIMAG( r(i) )
          rx(2*i-1) =  REAL( x(i) )
          rx(2*i-0) =  AIMAG( x(i) )
          rb(2*i-1) =  REAL( b(i) )
          rb(2*i-0) =  AIMAG( b(i) )
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE CCG
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Uzawa algorithm as a smoother.
!------------------------------------------------------------------------------
      SUBROUTINE Uzawa( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        USE LinearAlgebra

        TYPE(Matrix_t), POINTER :: A,M
        INTEGER :: Rounds, n
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,l,t,k0,k1,dofs,nn,elem,it
        INTEGER, POINTER :: ind(:)
        TYPE(Variable_t), POINTER :: Var
        REAL(KIND=dp) :: st_norm
        REAL(KIND=dp), ALLOCATABLE :: px(:), pb(:), pr(:)
!------------------------------------------------------------------------------

         Var => VariableGet( Mesh % Variables,  &
               CurrentModel % Solver % Variable % Name, ThisOnly=.TRUE. )

         dofs = Var % DOFs

         ALLOCATE( px(n), pr(n), pb(n) )

         CALL MGMv( A,x,r )
         r(1:n) = b(1:n) - r(1:n)
         st_norm = SQRT(SUM(r**2))

IF ( Rounds == 0 ) RETURN

DO it=1,200
         k = 0
         DO i=dofs,n,dofs
           k = k + 1
           px(k) = x(i) 
           pb(k) = b(i) 
           pr(k) = r(i) 
         END DO

         k = 0
         r = b
         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           DO j=A % Rows(i),A % Rows(i+1)-1
             IF ( MOD(A % Cols(j),dofs)==0 ) CYCLE
             r(i) = r(i) - A % Values(j)*x(A % Cols(j))
           END DO
         END DO

         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           r(i) = r(i)/(rphi*A % Values(A % Diag(i)))
         END DO


         k = 0
         DO i=dofs,n,dofs
           k = k+1
           DO j=A % Rows(i),A % Rows(i+1)-1
             IF ( MOD(A % Cols(j),dofs)==0 ) CYCLE
             pb(k) = pb(k) - A % Values(j)*x(A % Cols(j))
             pb(k) = pb(k) - A % Values(j)*r(A % Cols(j))
           END DO
         END DO

         CALL BiCGUzawa( k,A,M,px,pb,pr,Rounds, st_norm*1.d-2 )
         k = 0
         DO i=dofs,n,dofs
           k = k + 1
           x(i) = px(k)
         END DO

         px = x(1:n)
         k = 0
         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           k = k+1
           px(i) = px(i) + r(i)
         END DO

         k = 0
         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           k = k+1
           DO j=A % Rows(i)+dofs-1,A % Rows(i+1)-1,dofs
             px(i) = px(i) - A % Values(j)*x(A % Cols(j)) / &
                   (rphi*A % Values(A % Diag(i)))
           END DO
         END DO
         x(1:n) = px

         CALL MGMv( A,x,r )
         r(1:n) = b(1:n) - r(1:n)
         PRINT*,'AAAAAAAAAA: ', it, Rounds, st_norm*0.5_dp, SQRT(SUM(r**2))

         IF ( it > Rounds ) THEN
           IF ( SQRT(SUM(r**2)) < 0.5_dp*st_norm ) EXIT
         END IF
END DO
!------------------------------------------------------------------------------
      END SUBROUTINE Uzawa
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE UzawaMv( A, x, b )
!------------------------------------------------------------------------------
        TYPE(Matrix_t) :: A
        INTEGER :: i,j,k,n,q,dofs=3
        REAL(KIND=dp) CONTIG :: x(:),b(:)
        REAL(KIND=dp), ALLOCATABLE :: temp(:)
!------------------------------------------------------------------------------
         n = A % NumberOfRows
         ALLOCATE(temp(n))

         b = 0._dp
         temp = 0._dp

         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           DO j=A % Rows(i)+dofs-1,A % Rows(i+1)-1,dofs
             q = A % Cols(j)/dofs
             temp(i) = temp(i) + A % Values(j)*x(q)
           END DO
         END DO

         DO i=1,n
           IF ( MOD(i,dofs)==0 ) CYCLE
           temp(i) = temp(i)/(rphi*A % Values(A % Diag(i)))
         END DO

         k = 0
         DO i=dofs,n,dofs
           k = k + 1
           DO j=A % Rows(i),A % Rows(i+1)-1
             IF ( MOD(A % Cols(j),dofs)==0 ) CYCLE
             b(k) = b(k) - A % Values(j)*temp(A % Cols(j))
           END DO
         END DO

         k = 0
         DO i=dofs,n,dofs
           k = k + 1
           DO j=A % Rows(i)+dofs-1,A % Rows(i+1)-1,dofs
             q = A % Cols(j)/dofs
             b(k) = b(k) + A % Values(j)*x(q)
           END DO
         END DO
!------------------------------------------------------------------------------
     END SUBROUTINE UzawaMv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE UzawaPcond( A, b )
!------------------------------------------------------------------------------
        TYPE(Matrix_t) :: A
        REAL(KIND=dp) CONTIG :: b(:)
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: s
       INTEGER :: i,k,dofs=3,n

        n = A % NumberOfRows
        k = 0
        DO i=dofs,n,dofs
          k = k + 1
          s = A % Values(A % Diag(i))
          IF ( ABS(s) > 100*AEPS ) b(k) = b(k) / s
        END DO
!------------------------------------------------------------------------------
     END SUBROUTINE UzawaPcond
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE BiCGUzawa( n, A, M, x, b, r, Rounds, reps )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A,M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,n
        REAL(KIND=dp) :: alpha,beta,omega,rho,oldrho, reps
        REAL(KIND=dp), ALLOCATABLE :: Pr(:), Ri(:), V(:), S(:), T(:), T1(:), T2(:)
!------------------------------------------------------------------------------
        ALLOCATE( Pr(n), Ri(n), T(n), T1(n), T2(n), S(n), V(n) )
        PR = 0; Ri=0; T=0; T1=0; T2=0; S=0; V=0

        CALL UzawaMv( A, x, r )
        r(1:n) = b(1:n) - r(1:n)
        
        Ri(1:n) = r(1:n)
        Pr(1:n) = 0
        V(1:n) = 0
        omega  = 1
        alpha  = 0
        oldrho = 1
        
        DO i=1,200
          rho = MGdot( n, r, Ri )
          
          beta = alpha * rho / ( oldrho * omega )
          Pr(1:n) = r(1:n) + beta * (Pr(1:n) - omega*V(1:n))
          
          V(1:n) = Pr(1:n)
          CALL UzawaPcond( A,V )
          T1(1:n) = V(1:n)
          CALL UzawaMv( A, T1, V )

          alpha = rho / MGdot( n, Ri, V )
          S(1:n) = r(1:n) - alpha * V(1:n)
          
          T(1:n) = S(1:n)
          CALL UzawaPcond( A,T )
          T2(1:n) = T(1:n)
          CALL UzawaMv( A, T2, T )

          omega = MGdot( n,T,S ) / MGdot( n,T,T )
          oldrho = rho
          r(1:n) = S(1:n) - omega*T(1:n)
          x(1:n) = x(1:n) + alpha*T1(1:n) + omega*T2(1:n)

          CALL UzawaMv( A, x, t1 )
          t1(1:n) = b(1:n) - t1(1:n)
          IF ( SQRT(SUM(t1**2)) < reps ) EXIT
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE BiCGUzawa
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE BiCG( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A,M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,n
        REAL(KIND=dp) :: alpha,beta,omega,rho,oldrho
!------------------------------------------------------------------------------

        IF ( Rounds <= 0 ) RETURN

        CALL MGmv( A, x, r )
        r(1:n) = b(1:n) - r(1:n)
        
        Ri(1:n) = r(1:n)
        Pr(1:n) = 0
        V(1:n) = 0
        omega  = 1
        alpha  = 0
        oldrho = 1
        
        DO i=1,Rounds
          rho = MGdot( n, r, Ri )
          
          beta = alpha * rho / ( oldrho * omega )
          Pr(1:n) = r(1:n) + beta * (Pr(1:n) - omega*V(1:n))
          
          V(1:n) = Pr(1:n)
          CALL CRS_LUSolve( n, M, V )
          T1(1:n) = V(1:n)
          CALL MGmv( A, T1, V )

          alpha = rho / MGdot( n, Ri, V )
          
          S(1:n) = r(1:n) - alpha * V(1:n)
          
          T(1:n) = S(1:n)
          CALL CRS_LUSolve( n, M, T )
          T2(1:n) = T(1:n)
          CALL MGmv( A, T2, T )
          omega = MGdot( n,T,S ) / MGdot( n,T,T )

          oldrho = rho
          r(1:n) = S(1:n) - omega*T(1:n)
          x(1:n) = x(1:n) + alpha*T1(1:n) + omega*T2(1:n)
        END DO
!------------------------------------------------------------------------------
      END SUBROUTINE BiCG
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE Vanka( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
USE linearalgebra
        TYPE(Matrix_t), POINTER :: A,M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,l,n,t,nn,k1,k2,elem,nsize,it
        REAL(KIND=dp) :: s

        LOGICAL :: NS

        TYPE(Variable_t), POINTER :: Var
        TYPE(Element_t), POINTER :: element

        INTEGER, ALLOCATABLE :: ind(:)
        REAL(KIND=dp), ALLOCATABLE :: AL(:,:), h(:)
!------------------------------------------------------------------------------
        Var => VariableGet( Mesh % Variables, &
                 CurrentModel % Solver % Variable % Name, ThisOnly=.TRUE. )

        NS = ( GetVarName( Var ) == 'flow solution' ) 
        
        elem = Mesh % NumberOfBulkElements
        nsize = Mesh % MaxElementDOFs*DOFs
        ALLOCATE( AL(nsize,nsize), ind(nsize), h(nsize) )
        AL = 0._dp

        DO it=1,Rounds
          DO t=1,elem
            element => Mesh % Elements(t)
            IF ( ANY(Var % Perm(Element % NodeIndexes)<=0) ) CYCLE

            nn = Element % TYPE % NumberOfnodes
            nsize = nn*DOFs

            k = 0
            DO i=1,nn
              l = Element % Nodeindexes(i)
              DO j=1,DOFs
                k = k +1
                ind(k) = DOFs*(Var % Perm(l)-1)+j
              END DO
            END DO

            DO i=1,nsize
              j = ind(i)
              s = 0._dp
              DO k=A % Rows(j),A % Rows(j+1)-1
                s = s + A % Values(k)*x(A % Cols(k))
              END DO
              h(i) = b(j)-s
            END DO

            IF (NS) THEN
              DO i=1,nsize
                DO j=1,nsize
                  AL(i,j) = CRS_GetMatrixElement( A,ind(i),ind(j) )
                END DO
              END DO
              CALL SolveLinSys( nsize,SIZE(AL,1),AL,h )
            ELSE
              h(1:nsize)=h(1:nsize)/A % Values(A % Diag(ind(1:nsize)))
            END IF

            x(ind(1:nsize))=x(ind(1:nsize))+h(1:nsize)
          END DO
        END DO

        DEALLOCATE( AL, h, ind )
!------------------------------------------------------------------------------
      END SUBROUTINE Vanka
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE SolveLinSys( N,LDa,A,x )
!------------------------------------------------------------------------------
        INTEGER  N,IPIV(N),LDa,info
        DOUBLE PRECISION  A(LDa,*),x(n)

        IF ( N .LE. 0 ) RETURN
        CALL DGETRF( N,N,A,LDa,IPIV,INFO )
        CALL DGETRS( 'N',N,1,A,LDa,IPIV,X,N,INFO )
!------------------------------------------------------------------------------
      END SUBROUTINE SolveLinSys
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Gauss-Seidel smoother with plenty of comments for testing purposes.
!------------------------------------------------------------------------------
      SUBROUTINE TestGS( n, A, M, x, b, r, Rounds )
!------------------------------------------------------------------------------
        TYPE(Matrix_t), POINTER :: A, M
        INTEGER :: Rounds
        REAL(KIND=dp) CONTIG :: x(:),b(:),r(:)
!------------------------------------------------------------------------------
        INTEGER :: i,j,k,l,n,o,nsize
        REAL(KIND=dp) :: s
        INTEGER, POINTER CONTIG :: Cols(:),Rows(:)
        REAL(KIND=dp), POINTER CONTIG :: Values(:)
!------------------------------------------------------------------------------
     
        n = A % NumberOfRows
        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values
        
        PRINT *,'TestGS: Starting',&
            ASSOCIATED(Rows),ASSOCIATED(Cols),ASSOCIATED(Values),ASSOCIATED(M % diag)
        PRINT *,'TestGS: Sizes',&
            SIZE(Cols),SIZE(Values),SIZE(Rows)-1,SIZE(M % diag),SIZE(b),SIZE(r)
        PRINT *,'TestGS: MinVal',&
            MINVAL(Rows),MINVAL(Cols),MINVAL(Values),MINVAL(M % diag)
        PRINT *,'TestGS: MaxVal',&
            MAXVAL(Rows),MAXVAL(Cols),MAXVAL(Values),MAXVAL(M % diag)
        
        nsize = SIZE(Cols)
        
        
        DO k=1,Rounds
          
          DO i=1,n
            s = 0.0d0
            
            DO j=Rows(i),Rows(i+1)-1
              IF(j<1 .OR. j>nsize) THEN
                PRINT *,'TestGs A:',i,j
              END IF
              o = Cols(j)
              IF(o<1 .OR. o>n) THEN
                PRINT *,'TestGs B:',i,j,o
              END IF
              s = s + x(o) * Values(j)
            END DO
            
            l = M % diag(i)
            IF(l<1 .OR. l>nsize) THEN
              PRINT *,'TestGs C:',i,j,o,l
            END IF
            r(i) = (b(i)-s) / M % Values(l)
            x(i) = x(i) + r(i)
          END DO
        END DO
        
        PRINT *,'TestGS: Finished'
!------------------------------------------------------------------------------
      END SUBROUTINE TestGS
!------------------------------------------------------------------------------

    END FUNCTION MGSmooth

END MODULE Smoothers

!> \}
