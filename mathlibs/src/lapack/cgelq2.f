      SUBROUTINE CGELQ2( M, N, A, LDA, TAU, WORK, INFO )
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
*  CGELQ2 computes an LQ factorization of a complex m by n matrix A:
*  A = L * Q.
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
*          On exit, the elements on and below the diagonal of the array
*          contain the m by min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX array, dimension (M)
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
*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
*  A(i,i+1:n), and tau in TAU(i).
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
      EXTERNAL           CLACGV, CLARF, CLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
         CALL XERBLA( 'CGELQ2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
*
         CALL CLACGV( N-I+1, A( I, I ), LDA )
         ALPHA = A( I, I )
         CALL CLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA,
     $                TAU( I ) )
         IF( I.LT.M ) THEN
*
*           Apply H(i) to A(i+1:m,i:n) from the right
*
            A( I, I ) = ONE
            CALL CLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ),
     $                  A( I+1, I ), LDA, WORK )
         END IF
         A( I, I ) = ALPHA
         CALL CLACGV( N-I+1, A( I, I ), LDA )
   10 CONTINUE
      RETURN
*
*     End of CGELQ2
*
      END
