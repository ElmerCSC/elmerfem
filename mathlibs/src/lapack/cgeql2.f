      SUBROUTINE CGEQL2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGEQL2 computes a QL factorization of a complex m by n matrix A:
*  A = Q * L.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, if m >= n, the lower triangle of the subarray
*          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L;
*          if m <= n, the elements on and below the (n-m)-th
*          superdiagonal contain the m by n lower trapezoidal matrix L;
*          the remaining elements, with the array TAU, represent the
*          unitary matrix Q as a product of elementary reflectors
*          (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k) . . . H(2) H(1), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
*  A(1:m-k+i-1,n-k+i), and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF, CLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEQL2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = K, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        A(1:m-k+i-1,n-k+i)
*
         ALPHA = A( M-K+I, N-K+I )
         CALL CLARFG( M-K+I, ALPHA, A( 1, N-K+I ), 1, TAU( I ) )
*
*        Apply H(i)' to A(1:m-k+i,1:n-k+i-1) from the left
*
         A( M-K+I, N-K+I ) = ONE
         CALL CLARF( 'Left', M-K+I, N-K+I-1, A( 1, N-K+I ), 1,
     $               CONJG( TAU( I ) ), A, LDA, WORK )
         A( M-K+I, N-K+I ) = ALPHA
   10 CONTINUE
      RETURN
*
*     End of CGEQL2
*
      END
