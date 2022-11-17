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
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!> Linear Algebra: LU-decomposition & matrix inverse & nonsymmetric full 
!>  matrix eigenvalues  (don't use this for anything big, use for example 
!>  LAPACK routines instead...)
!------------------------------------------------------------------------------
MODULE LinearAlgebra

  USE Types
  IMPLICIT NONE

 CONTAINS

  SUBROUTINE InvertMatrix( A,n )

    INTEGER :: n 
    REAL(KIND=dp) :: A(:,:)

    REAL(KIND=dp) :: s
    INTEGER :: i,j,k
    INTEGER :: pivot(n)
    LOGICAL :: erroneous

   ! /*
   !  *  AP = LU
   !  */
   CALL LUDecomp( a,n,pivot,erroneous  )

   IF (erroneous) CALL Fatal('InvertMatrix', 'inversion needs successfull LU-decomposition')

   ! /*  
   !  *  INV(U)
   !  */
   DO i=N-1,1,-1
     DO j=N,i+1,-1
       s = -A(i,j)
       DO k=i+1,j-1
         s = s - A(i,k)*A(k,j)
       END DO
       A(i,j) = s
     END DO
   END DO

   ! /*
   !  * INV(L)
   !  */
   DO i=n-1,1,-1
     DO j=n,i+1,-1
       s = 0.D00
       DO k=i+1,j
         s = s - A(j,k)*A(k,i)
       END DO
       A(j,i) = A(i,i)*s
     END DO
   END DO
  
   ! /* 
   !  * A  = INV(AP)
   !  */
   DO i=1,n
     DO j=1,n
       s = 0.0D0
       DO k=MAX(i,j),n
         IF ( k /= i ) THEN
           s = s + A(i,k)*A(k,j)
         ELSE
           s = s + A(k,j)
         END IF
       END DO
       A(i,j) = s
     END DO
   END DO

   ! /*
   !  * A = INV(A) (at last)
   !  */
   DO i=n,1,-1
     IF ( pivot(i) /= i ) THEN
       DO j = 1,n
         s = A(i,j)
         A(i,j) = A(pivot(i),j)
         A(pivot(i),j) = s
       END DO
     END IF
   END DO

 END SUBROUTINE InvertMatrix


 SUBROUTINE LUSolve( n,A,x,pivot_in )
   REAL(KIND=dp) :: A(:,:)
   REAL(KIND=dp) :: x(n)
   INTEGER :: n 
   INTEGER, OPTIONAL :: pivot_in(:)

   REAL(KIND=dp) :: s
   INTEGER :: i,j
   INTEGER :: pivot(n)
   LOGICAL :: erroneous

   ! /*
   !  *  AP = LU
   !  */
   IF(PRESENT(pivot_in)) THEN
     Pivot(1:n) = Pivot_in(1:n)
   ELSE
     CALL LUDecomp( A,n,pivot,erroneous )
     IF (erroneous) CALL Fatal('LUSolve', 'LU-decomposition fails')
   END IF

   ! Forward substitute
   DO i=1,n
      s = x(i)
      DO j=1,i-1
         s = s - A(i,j) * x(j)
      END DO
      x(i) = A(i,i) * s 
   END DO

   !
   ! Backward substitute (solve x from Ux = z)
   DO i=n,1,-1
      s = x(i)
      DO j=i+1,n
         s = s - A(i,j) * x(j)
      END DO
      x(i) = s
   END DO

   DO i=n,1,-1
      IF ( pivot(i) /= i ) THEN
         s = x(i)
         x(i) = x(pivot(i))
         x(pivot(i)) = s
      END IF
   END DO

  END SUBROUTINE LUSolve


!----------------------------------------------------------------------
!> LU- decomposition by gaussian elimination. Row pivoting is used.
!  
!> result : AP = L'U ; L' = LD; pivot[i] is the swapped column number
!> for column i.
! 
!> Result is stored in place of original matrix.
!----------------------------------------------------------------------
  SUBROUTINE LUDecomp( a,n,pivot,erroneous )

    REAL(KIND=dp), DIMENSION (:,:) :: a
    INTEGER :: n
    INTEGER, DIMENSION (:) :: pivot
    LOGICAL, OPTIONAL :: erroneous

    INTEGER :: i,j,k,l
    REAL(KIND=dp) :: swap

    IF (PRESENT(erroneous)) erroneous = .FALSE.

    DO i=1,n
      j = i
      DO k=i+1,n
        IF ( ABS(A(i,k)) > ABS(A(i,j)) ) j = k
      END DO

      IF ( ABS(A(i,j)) == 0.0d0) THEN
        CALL Error( 'LUDecomp', 'Matrix is singular.' )
        IF (PRESENT(erroneous)) erroneous = .TRUE.
        RETURN
      END IF

      pivot(i) = j

      IF ( j /= i ) THEN
        DO k=1,i
          swap = A(k,j)
          A(k,j) = A(k,i)
          A(k,i) = swap
        END DO
      END IF

      DO k=i+1,n
        A(i,k) = A(i,k) / A(i,i)
      END DO

      DO k=i+1,n
        IF ( j /= i ) THEN
          swap = A(k,i)
          A(k,i) = A(k,j)
          A(k,j) = swap
        END IF

        DO  l=i+1,n
          A(k,l) = A(k,l) - A(k,i) * A(i,l)
        END DO
      END DO
    END DO

    pivot(n) = n
    IF ( ABS(A(n,n)) == 0.0d0 ) THEN
      CALL Error( 'LUDecomp',  'Matrix is (at least almost) singular.' )
    END IF

    DO i=1,n
      IF ( ABS(A(i,i)) == 0.0d0 ) THEN
        CALL Error( 'LUSolve', 'Matrix is singular.' )
        IF (PRESENT(erroneous)) erroneous = .TRUE.
        RETURN       
      END IF
      A(i,i) = 1.0_dp / A(i,i)
    END DO

  END SUBROUTINE LUDecomp



  SUBROUTINE ComplexInvertMatrix( A,n )

    COMPLEx(KIND=dp), DIMENSION(:,:) :: A
    INTEGER :: n 

    COMPLEX(KIND=dp) :: s
    INTEGER :: i,j,k
    INTEGER :: pivot(n)
    LOGICAL :: erroneous

   ! /*
   !  *  AP = LU
   !  */
   CALL ComplexLUDecomp( a,n,pivot,erroneous )

   IF (erroneous) CALL Fatal('ComplexInvertMatrix', 'inversion needs successfull LU-decomposition')

   DO i=1,n
     IF ( ABS(A(i,i))==0.0d0 ) THEN
       CALL Error( 'ComplexInvertMatrix', 'Matrix is singular.' )
       RETURN       
     END IF
     A(i,i) = 1.0D0/A(i,i)
   END DO

   ! /*  
   !  *  INV(U)
   !  */
   DO i=N-1,1,-1
     DO j=N,i+1,-1
       s = -A(i,j)
       DO k=i+1,j-1
         s = s - A(i,k)*A(k,j)
       END DO
       A(i,j) = s
     END DO
   END DO

   ! /*
   !  * INV(L)
   !  */
   DO i=n-1,1,-1
     DO j=n,i+1,-1
       s = 0.D00
       DO k=i+1,j
         s = s - A(j,k)*A(k,i)
       END DO
       A(j,i) = A(i,i)*s
     END DO
   END DO
  
   ! /* 
   !  * A  = INV(AP)
   !  */
   DO i=1,n
     DO j=1,n
       s = 0.0D0
       DO k=MAX(i,j),n
         IF ( k /= i ) THEN
           s = s + A(i,k)*A(k,j)
         ELSE
           s = s + A(k,j)
         END IF
       END DO
       A(i,j) = s
     END DO
   END DO

   ! /*
   !  * A = INV(A) (at last)
   !  */
   DO i=n,1,-1
     IF ( pivot(i) /= i ) THEN
       DO j = 1,n
         s = A(i,j)
         A(i,j) = A(pivot(i),j)
         A(pivot(i),j) = s
       END DO
     END IF
   END DO

 END SUBROUTINE ComplexInvertMatrix

!-------------------------------------------------------------------------
!> LU- decomposition by gaussian elimination for complex valued linear system.
!> Row pivoting is used.
! 
!> result : AP = L'U ; L' = LD; pivot[i] is the swapped column number
!> for column i.
!
!> Result is stored in place of original matrix.
!-------------------------------------------------------------------------

  SUBROUTINE ComplexLUDecomp( a,n,pivot,erroneous )

    COMPLEX(KIND=dp), DIMENSION (:,:) :: a
    INTEGER :: n
    INTEGER, DIMENSION (:) :: pivot
    LOGICAL, OPTIONAL :: erroneous

    INTEGER :: i,j,k,l
    COMPLEX(KIND=dp) :: swap

    IF (PRESENT(erroneous)) erroneous = .FALSE.

    DO i=1,n
      j = i
      DO k=i+1,n
        IF ( ABS(A(i,k)) > ABS(A(i,j)) ) j = k
      END DO

      IF ( ABS(A(i,j))==0.0d0 ) THEN
        CALL Error( 'ComplexLUDecomp', 'Matrix is singular.' )
        IF (PRESENT(erroneous)) erroneous = .TRUE.        
        RETURN
      END IF

      pivot(i) = j

      IF ( j /= i ) THEN
        DO k=1,i
          swap = A(k,j)
          A(k,j) = A(k,i)
          A(k,i) = swap
        END DO
      END IF

      DO k=i+1,n
        A(i,k) = A(i,k) / A(i,i)
      END DO

      DO k=i+1,n
        IF ( j /= i ) THEN
          swap = A(k,i)
          A(k,i) = A(k,j)
          A(k,j) = swap
        END IF

        DO  l=i+1,n
          A(k,l) = A(k,l) - A(k,i) * A(i,l)
        END DO
      END DO
    END DO

    pivot(n) = n
    IF ( ABS(A(n,n))==0.0d0 ) THEN
      CALL Error( 'ComplexLUDecomp', 'Matrix is (at least almost) singular.' )
    END IF
  END SUBROUTINE ComplexLUDecomp

  ! --------------------------------------------------
  !> Solve eigenvalues of a nonsymmetric matrix A(n,n)
  !> The matrix is modified in the process.
  ! --------------------------------------------------
  SUBROUTINE EigenValues( A, n, Vals )
    IMPLICIT NONE
    REAL(KIND=dp) :: A(:,:)
    INTEGER :: n
    COMPLEX(KIND=dp) :: Vals(:)

    INTEGER :: iter, i, j, k
    REAL(KIND=dp) :: b, s, t

    IF (n==1 ) THEN
      Vals(1) = A(1,1)
      RETURN
    END IF

    CALL Hesse(A, n)

    DO iter=1,1000
      DO i=1,n-1
        s = AEPS * ( ABS(A(i,i))+ABS(A(i+1,i+1)) )
        IF ( ABS(A(i+1,i)) < s ) A(i+1,i) = 0.0
      END DO
    
      i = 1
      DO WHILE( .TRUE. )
        DO j=i,n-1
          IF (A(j+1,j) /= 0 ) EXIT
        END DO

        DO k=j,n-1
          IF (A(k+1,k) == 0 ) EXIT
        END DO
        i = k;
        IF ( i >= n .OR. k-j+1 >= 3 ) EXIT
      END DO
    
      IF (k-j+1 < 3) EXIT
      CALL Francis(A(j:,j:), k-j+1);
    END DO

    j = 0
    i = 1
    DO WHILE( i<n ) 
      IF (A(i+1,i) == 0) THEN
        j = j + 1
        Vals(j) = A(i,i)
      ELSE
        b = ( A(i,i)+A(i+1,i+1) ); s=b*b
        t = A(i,i)*A(i+1,i+1) - A(i,i+1)*A(i+1,i)
        s = s - 4*t
        IF  (s >= 0) THEN
          s = SQRT(s)
          j = j + 1
          IF ( b>0 ) THEN
            Vals(j) = (b+s)/2
          ELSE
            Vals(j) = 2*t/(b-s)
          END IF
          j = j + 1
          IF ( b > 0 ) THEN
            Vals(j) = 2*t/(b+s)
          ELSE
            Vals(j) = (b-s)/2
          END IF
        ELSE
          j = j + 1
          Vals(j)  = CMPLX(b, SQRT(-s),KIND=dp)/2
          j = j + 1
          Vals(j)  = CMPLX(b,-SQRT(-s),KIND=dp)/2
        END IF
        i = i + 1
      END IF
      i = i + 1
    END DO
  
    IF (A(n, n-1) == 0) THEN
      j = j + 1
      Vals(j) = A(n,n)
    END IF

CONTAINS

    SUBROUTINE vbcalc( x,v,b,beg,end )
      IMPLICIT NONE
      REAL(KIND=dp) ::  x(:),v(:),b
      INTEGER ::  beg, end

      REAL(KIND=dp) :: alpha,m
      INTEGER :: i

      m = MAXVAL( ABS(x(beg:end)) )

      IF ( m == 0 ) THEN
        v(beg:end) = 0
      ELSE
        alpha = 0
        m = 1 / m
        DO i=beg,end
          v(i) = m*x(i)
          alpha = alpha+v(i)*v(i)
        END DO
        alpha = SQRT(alpha)
        b = 1/(alpha*(alpha+ABS(v(beg))))
        IF ( v(beg) < 0 ) THEN
          v(beg) = v(beg) - alpha
        ELSE
          v(beg) = v(beg) + alpha
        END IF
      END IF
    END SUBROUTINE vbcalc


    SUBROUTINE Hesse(H, dim)
      IMPLICIT NONE
      INTEGER ::  dim
      REAL(KIND=dp) :: H(:,:)

      REAL(KIND=dp) :: b, s
      INTEGER ::  i, j, k;
      REAL(KIND=dp), ALLOCATABLE ::  v(:), x(:)

      ALLOCATE( x(dim), v(dim) )

      DO i=1,dim-1
        DO j=i+1,dim
          x(j) = H(j,i)
        END DO

        CALL vbcalc(x, v, b, i+1, dim )

        IF (v(i+1) == 0) EXIT

        DO j=i+2, dim
          x(j) = v(j)/v(i+1)
          v(j) = b*v(i+1)*v(j)
        END DO
        v(i+1) = b*v(i+1)*v(i+1)

        DO j=1,dim
          s = 0.0
          DO k=i+1,dim
            s = s + H(j,k)*v(k)
          END DO
          H(j,i+1) = H(j,i+1) - s
          DO k=i+2,dim
            H(j,k) = H(j,k) - s*x(k)
          END DO
        END DO

        DO j=1,dim
          s = H(i+1,j)
          DO k=i+2,dim
            s = s + H(k,j)*x(k)
          END DO
          DO k=i+1,dim
            H(k,j) = H(k,j) - s*v(k)
          END DO
        END DO
        H(i+2:dim,i) = 0
      END DO

      DEALLOCATE( x, v )
    END SUBROUTINE Hesse


    SUBROUTINE Francis( H, dim )
      IMPLICIT NONE
      INTEGER :: dim
      REAL(KIND=dp) :: H(:,:)

      INTEGER ::  i, i1, j, k, m, n, end
      REAL(KIND=dp) :: x(3), v(3), b, s, t, v0i

      n = dim; m = n-1

      t = H(m,m)*H(n,n) - H(m,n)*H(n,m)
      s = H(m,m)+H(n,n)

      x(1) = H(1,1)*H(1,1)+H(1,2)*H(2,1)-s*H(1,1)+t
      x(2) = H(2,1)*(H(1,1)+H(2,2)-s)
      x(3) = H(2,1)*H(3,2);

      CALL vbcalc(x, v, b, 1, 3)

      IF (v(1) == 0) RETURN

      ! DO i=2,3
      !   x(i) = v(i)/v(1);
      !   v(i) = b * v(1) * v(i);
      ! END DO
      x(2) = v(2)/v(1)
      v(2) = b*v(1)*v(2)
      x(3) = v(3)/v(1)
      v(3) = b*v(1)*v(3)
      x(1) = 1; v(1) = b*v(1)*v(1)

      DO i=1,dim
        ! DO j=1,3
        !   s = s + H(i,j) * v(j)
        ! END DO
        s = H(i,1)*v(1) + H(i,2)*v(2) + H(i,3)*v(3)
   
        ! DO j=1,3
        !   H(i,j) = H(i,j) - s * x(j)
        ! END DO
        H(i,1) = H(i,1) - s
        H(i,2) = H(i,2) - s*x(2)
        H(i,3) = H(i,3) - s*x(3)
      END DO

      DO i=1,dim
        ! s = 0.0d0
        ! DO j=1,3
        !   s = s + H(j,i) * x(j)
        ! END DO
        s = h(1,i) + H(2,i)*x(2) + H(3,i)*x(3)

        ! DO j=1,3
        !   H(j,i) = H(j,i) - s * v(j)
        ! END DO
        H(1,i) = H(1,i) - s*v(1)
        H(2,i) = H(2,i) - s*v(2)
        H(3,i) = H(3,i) - s*v(3)
      END DO

      DO i=1,dim-1
        end = MIN(2, dim-i-1)
        DO j=1,end+1
          x(j) = H(i+j,i)
        END DO

        CALL vbcalc(x, v, b, 1, end+1)

        IF (v(1) == 0) EXIT

        DO j=2,end+1
          x(j) = v(j)/v(1);
          v(j) = b*v(1)*v(j)
        END DO
        x(1) = 1; v(1) = b*v(1)*v(1)

        DO j=1,dim
          s = 0.0
          DO k=1,end+1
            s = s + H(j,i+k)*v(k)
          END DO
          DO k=1,end+1
            H(j,i+k) = H(j,i+k) - s*x(k)
          END DO
        END DO

        DO j=1,dim
          s = H(i+1,j)
          DO k=2,end+1
            s = s + H(i+k,j)*x(k)
          END DO
          DO k=1,end+1
            H(i+k,j) = H(i+k,j) - s*v(k)
          END DO
        END DO
        H(i+2:dim,i) = 0
      END DO
    END SUBROUTINE Francis
  END SUBROUTINE EigenValues

END MODULE LinearAlgebra

!> \} ElmerLib
