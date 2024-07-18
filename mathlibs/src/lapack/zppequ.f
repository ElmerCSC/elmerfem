      SUBROUTINE ZPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
      DOUBLE PRECISION   AMAX, SCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   S( * )
      COMPLEX*16         AP( * )
*     ..
*
*  Purpose
*  =======
*
*  ZPPEQU computes row and column scalings intended to equilibrate a
*  Hermitian positive definite matrix A in packed storage and reduce
*  its condition number (with respect to the two-norm).  S contains the
*  scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix
*  B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal.
*  This choice of S puts the condition number of B within a factor N of
*  the smallest possible condition number over all possible diagonal
*  scalings.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The upper or lower triangle of the Hermitian matrix A, packed
*          columnwise in a linear array.  The j-th column of A is stored
*          in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*
*  S       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, S contains the scale factors for A.
*
*  SCOND   (output) DOUBLE PRECISION
*          If INFO = 0, S contains the ratio of the smallest S(i) to
*          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
*          large nor too small, it is not worth scaling by S.
*
*  AMAX    (output) DOUBLE PRECISION
*          Absolute value of largest matrix element.  If AMAX is very
*          close to overflow or very close to underflow, the matrix
*          should be scaled.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the i-th diagonal element is nonpositive.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, JJ
      DOUBLE PRECISION   SMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPPEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Initialize SMIN and AMAX.
*
      S( 1 ) = DBLE( AP( 1 ) )
      SMIN = S( 1 )
      AMAX = S( 1 )
*
      IF( UPPER ) THEN
*
*        UPLO = 'U':  Upper triangle of A is stored.
*        Find the minimum and maximum diagonal elements.
*
         JJ = 1
         DO 10 I = 2, N
            JJ = JJ + I
            S( I ) = DBLE( AP( JJ ) )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   10    CONTINUE
*
      ELSE
*
*        UPLO = 'L':  Lower triangle of A is stored.
*        Find the minimum and maximum diagonal elements.
*
         JJ = 1
         DO 20 I = 2, N
            JJ = JJ + N - I + 2
            S( I ) = DBLE( AP( JJ ) )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   20    CONTINUE
      END IF
*
      IF( SMIN.LE.ZERO ) THEN
*
*        Find the first non-positive diagonal element and return.
*
         DO 30 I = 1, N
            IF( S( I ).LE.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   30    CONTINUE
      ELSE
*
*        Set the scale factors to the reciprocals
*        of the diagonal elements.
*
         DO 40 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   40    CONTINUE
*
*        Compute SCOND = min(S(I)) / max(S(I))
*
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF
      RETURN
*
*     End of ZPPEQU
*
      END
