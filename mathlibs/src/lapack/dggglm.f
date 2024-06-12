      SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), WORK( * ),
     $                   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGGLM solves a general Gauss-Markov linear model (GLM) problem:
*
*          minimize || y ||_2   subject to   d = A*x + B*y
*              x
*
*  where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
*  given N-vector. It is assumed that M <= N <= M+P, and
*
*             rank(A) = M    and    rank( A B ) = N.
*
*  Under these assumptions, the constrained equation is always
*  consistent, and there is a unique solution x and a minimal 2-norm
*  solution y, which is obtained using a generalized QR factorization
*  of A and B.
*
*  In particular, if matrix B is square nonsingular, then the problem
*  GLM is equivalent to the following weighted linear least squares
*  problem
*
*               minimize || inv(B)*(d-A*x) ||_2
*                   x
*
*  where inv(B) denotes the inverse of B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrices A and B.  N >= 0.
*
*  M       (input) INTEGER
*          The number of columns of the matrix A.  0 <= M <= N.
*
*  P       (input) INTEGER
*          The number of columns of the matrix B.  P >= N-M.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
*          On entry, the N-by-M matrix A.
*          On exit, A is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,P)
*          On entry, the N-by-P matrix B.
*          On exit, B is destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, D is the left hand side of the GLM equation.
*          On exit, D is destroyed.
*
*  X       (output) DOUBLE PRECISION array, dimension (M)
*  Y       (output) DOUBLE PRECISION array, dimension (P)
*          On exit, X and Y are the solutions of the GLM problem.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N+M+P).
*          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
*          where NB is an upper bound for the optimal blocksizes for
*          DGEQRF, SGERQF, DORMQR and SORMRQ.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  ===================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, LOPT, LWKOPT, NB, NB1, NB2, NB3, NB4, NP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DGGQRF, DORMQR, DORMRQ, DTRSV,
     $                   XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NP = MIN( N, P )
      NB1 = ILAENV( 1, 'DGEQRF', ' ', N, M, -1, -1 )
      NB2 = ILAENV( 1, 'DGERQF', ' ', N, M, -1, -1 )
      NB3 = ILAENV( 1, 'DORMQR', ' ', N, M, P, -1 )
      NB4 = ILAENV( 1, 'DORMRQ', ' ', N, M, P, -1 )
      NB = MAX( NB1, NB2, NB3, NB4 )
      LWKOPT = M + NP + MAX( N, P )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LWORK.LT.MAX( 1, N+M+P ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGGLM', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the GQR factorization of matrices A and B:
*
*            Q'*A = ( R11 ) M,    Q'*B*Z' = ( T11   T12 ) M
*                   (  0  ) N-M             (  0    T22 ) N-M
*                      M                     M+P-N  N-M
*
*     where R11 and T22 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL DGGQRF( N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ),
     $             WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = WORK( M+NP+1 )
*
*     Update left-hand-side vector d = Q'*d = ( d1 ) M
*                                             ( d2 ) N-M
*
      CALL DORMQR( 'Left', 'Transpose', N, 1, M, A, LDA, WORK, D,
     $             MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      LOPT = MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
*     Solve T22*y2 = d2 for y2
*
      CALL DTRSV( 'Upper', 'No transpose', 'Non unit', N-M,
     $            B( M+1, M+P-N+1 ), LDB, D( M+1 ), 1 )
      CALL DCOPY( N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 )
*
*     Set y1 = 0
*
      DO 10 I = 1, M + P - N
         Y( I ) = ZERO
   10 CONTINUE
*
*     Update d1 = d1 - T12*y2
*
      CALL DGEMV( 'No transpose', M, N-M, -ONE, B( 1, M+P-N+1 ), LDB,
     $            Y( M+P-N+1 ), 1, ONE, D, 1 )
*
*     Solve triangular system: R11*x = d1
*
      CALL DTRSV( 'Upper', 'No Transpose', 'Non unit', M, A, LDA, D, 1 )
*
*     Copy D to X
*
      CALL DCOPY( M, D, 1, X, 1 )
*
*     Backward transformation y = Z'*y
*
      CALL DORMRQ( 'Left', 'Transpose', P, 1, NP,
     $             B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y,
     $             MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO )
      WORK( 1 ) = M + NP + MAX( LOPT, INT( WORK( M+NP+1 ) ) )
*
      RETURN
*
*     End of DGGGLM
*
      END
