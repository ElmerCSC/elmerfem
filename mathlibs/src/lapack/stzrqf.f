      SUBROUTINE STZRQF( M, N, A, LDA, TAU, INFO )
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
      REAL               A( LDA, * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  This routine is deprecated and has been replaced by routine STZRZF.
*
*  STZRQF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
*  to upper triangular form by means of orthogonal transformations.
*
*  The upper trapezoidal matrix A is factored as
*
*     A = ( R  0 ) * Z,
*
*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
*  triangular matrix.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= M.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the leading M-by-N upper trapezoidal part of the
*          array A must contain the matrix to be factorized.
*          On exit, the leading M-by-M upper triangular part of A
*          contains the upper triangular matrix R, and elements M+1 to
*          N of the first M rows of A, with the array TAU, represent the
*          orthogonal matrix Z as a product of M elementary reflectors.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) REAL array, dimension (M)
*          The scalar factors of the elementary reflectors.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The factorization is obtained by Householder's method.  The kth
*  transformation matrix, Z( k ), which is used to introduce zeros into
*  the ( m - k + 1 )th row of A, is given in the form
*
*     Z( k ) = ( I     0   ),
*              ( 0  T( k ) )
*
*  where
*
*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
*                                                 (   0    )
*                                                 ( z( k ) )
*
*  tau is a scalar and z( k ) is an ( n - m ) element vector.
*  tau and z( k ) are chosen to annihilate the elements of the kth row
*  of X.
*
*  The scalar tau is returned in the kth element of TAU and the vector
*  u( k ) in the kth row of A, such that the elements of z( k ) are
*  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in
*  the upper triangular part of A.
*
*  Z is given by
*
*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K, M1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMV, SGER, SLARFG, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STZRQF', -INFO )
         RETURN
      END IF
*
*     Perform the factorization.
*
      IF( M.EQ.0 )
     $   RETURN
      IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
      ELSE
         M1 = MIN( M+1, N )
         DO 20 K = M, 1, -1
*
*           Use a Householder reflection to zero the kth row of A.
*           First set up the reflection.
*
            CALL SLARFG( N-M+1, A( K, K ), A( K, M1 ), LDA, TAU( K ) )
*
            IF( ( TAU( K ).NE.ZERO ) .AND. ( K.GT.1 ) ) THEN
*
*              We now perform the operation  A := A*P( k ).
*
*              Use the first ( k - 1 ) elements of TAU to store  a( k ),
*              where  a( k ) consists of the first ( k - 1 ) elements of
*              the  kth column  of  A.  Also  let  B  denote  the  first
*              ( k - 1 ) rows of the last ( n - m ) columns of A.
*
               CALL SCOPY( K-1, A( 1, K ), 1, TAU, 1 )
*
*              Form   w = a( k ) + B*z( k )  in TAU.
*
               CALL SGEMV( 'No transpose', K-1, N-M, ONE, A( 1, M1 ),
     $                     LDA, A( K, M1 ), LDA, ONE, TAU, 1 )
*
*              Now form  a( k ) := a( k ) - tau*w
*              and       B      := B      - tau*w*z( k )'.
*
               CALL SAXPY( K-1, -TAU( K ), TAU, 1, A( 1, K ), 1 )
               CALL SGER( K-1, N-M, -TAU( K ), TAU, 1, A( K, M1 ), LDA,
     $                    A( 1, M1 ), LDA )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of STZRQF
*
      END
