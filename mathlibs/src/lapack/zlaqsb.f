      SUBROUTINE ZLAQSB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            KD, LDAB, N
      DOUBLE PRECISION   AMAX, SCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   S( * )
      COMPLEX*16         AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAQSB equilibrates a symmetric band matrix A using the scaling
*  factors in the vector S.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U'*U or A = L*L' of the band
*          matrix A, in the same storage format as A.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  S       (output) DOUBLE PRECISION array, dimension (N)
*          The scale factors for A.
*
*  SCOND   (input) DOUBLE PRECISION
*          Ratio of the smallest S(i) to the largest S(i).
*
*  AMAX    (input) DOUBLE PRECISION
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies whether or not equilibration was done.
*          = 'N':  No equilibration.
*          = 'Y':  Equilibration was done, i.e., A has been replaced by
*                  diag(S) * A * diag(S).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if scaling should be done
*  based on the ratio of the scaling factors.  If SCOND < THRESH,
*  scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if scaling should
*  be done based on the absolute size of the largest matrix element.
*  If AMAX > LARGE or AMAX < SMALL, scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
      ELSE
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A is stored in band format.
*
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = MAX( 1, J-KD ), J
                  AB( KD+1+I-J, J ) = CJ*S( I )*AB( KD+1+I-J, J )
   10          CONTINUE
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, MIN( N, J+KD )
                  AB( 1+I-J, J ) = CJ*S( I )*AB( 1+I-J, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of ZLAQSB
*
      END
