      SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGTTRS solves one of the systems of equations
*     A * X = B,  A**T * X = B,  or  A**H * X = B,
*  with a tridiagonal matrix A using the LU factorization computed
*  by ZGTTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) COMPLEX*16 array, dimension (N-1)
*          The (n-1) multipliers that define the matrix L from the
*          LU factorization of A.
*
*  D       (input) COMPLEX*16 array, dimension (N)
*          The n diagonal elements of the upper triangular matrix U from
*          the LU factorization of A.
*
*  DU      (input) COMPLEX*16 array, dimension (N-1)
*          The (n-1) elements of the first super-diagonal of U.
*
*  DU2     (input) COMPLEX*16 array, dimension (N-2)
*          The (n-2) elements of the second super-diagonal of U.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row IPIV(i).  IPIV(i) will always be either
*          i or i+1; IPIV(i) = i indicates a row interchange was not
*          required.
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the matrix of right hand side vectors B.
*          On exit, B is overwritten by the solution vectors X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            ITRANS, J, JB, NB
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGTTS2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOTRAN = ( TRANS.EQ.'N' .OR. TRANS.EQ.'n' )
      IF( .NOT.NOTRAN .AND. .NOT.( TRANS.EQ.'T' .OR. TRANS.EQ.
     $    't' ) .AND. .NOT.( TRANS.EQ.'C' .OR. TRANS.EQ.'c' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Decode TRANS
*
      IF( NOTRAN ) THEN
         ITRANS = 0
      ELSE IF( TRANS.EQ.'T' .OR. TRANS.EQ.'t' ) THEN
         ITRANS = 1
      ELSE
         ITRANS = 2
      END IF
*
*     Determine the number of right-hand sides to solve at a time.
*
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'ZGTTRS', TRANS, N, NRHS, -1, -1 ) )
      END IF
*
      IF( NB.GE.NRHS ) THEN
         CALL ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL ZGTTS2( ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ),
     $                   LDB )
   10    CONTINUE
      END IF
*
*     End of ZGTTRS
*
      END
