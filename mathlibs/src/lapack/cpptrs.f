      SUBROUTINE CPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX            AP( * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CPPTRS solves a system of linear equations A*X = B with a Hermitian
*  positive definite matrix A in packed storage using the Cholesky
*  factorization A = U**H*U or A = L*L**H computed by CPPTRF.
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
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**H*U or A = L*L**H, packed columnwise in a linear
*          array.  The j-th column of U or L is stored in the array AP
*          as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTPSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPPTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U'*U.
*
         DO 10 I = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL CTPSV( 'Upper', 'Conjugate transpose', 'Non-unit', N,
     $                  AP, B( 1, I ), 1 )
*
*           Solve U*X = B, overwriting B with X.
*
            CALL CTPSV( 'Upper', 'No transpose', 'Non-unit', N, AP,
     $                  B( 1, I ), 1 )
   10    CONTINUE
      ELSE
*
*        Solve A*X = B where A = L*L'.
*
         DO 20 I = 1, NRHS
*
*           Solve L*Y = B, overwriting B with X.
*
            CALL CTPSV( 'Lower', 'No transpose', 'Non-unit', N, AP,
     $                  B( 1, I ), 1 )
*
*           Solve L'*X = Y, overwriting B with X.
*
            CALL CTPSV( 'Lower', 'Conjugate transpose', 'Non-unit', N,
     $                  AP, B( 1, I ), 1 )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of CPPTRS
*
      END
