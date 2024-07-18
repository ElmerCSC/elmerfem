      SUBROUTINE SGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,
     $                   IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ),
     $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ),
     $                   FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  SGTRFS improves the computed solution to a system of linear
*  equations when the coefficient matrix is tridiagonal, and provides
*  error bounds and backward error estimates for the solution.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) REAL array, dimension (N-1)
*          The (n-1) subdiagonal elements of A.
*
*  D       (input) REAL array, dimension (N)
*          The diagonal elements of A.
*
*  DU      (input) REAL array, dimension (N-1)
*          The (n-1) superdiagonal elements of A.
*
*  DLF     (input) REAL array, dimension (N-1)
*          The (n-1) multipliers that define the matrix L from the
*          LU factorization of A as computed by SGTTRF.
*
*  DF      (input) REAL array, dimension (N)
*          The n diagonal elements of the upper triangular matrix U from
*          the LU factorization of A.
*
*  DUF     (input) REAL array, dimension (N-1)
*          The (n-1) elements of the first superdiagonal of U.
*
*  DU2     (input) REAL array, dimension (N-2)
*          The (n-2) elements of the second superdiagonal of U.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row IPIV(i).  IPIV(i) will always be either
*          i or i+1; IPIV(i) = i indicates a row interchange was not
*          required.
*
*  B       (input) REAL array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input/output) REAL array, dimension (LDX,NRHS)
*          On entry, the solution matrix X, as computed by SGTTRS.
*          On exit, the improved solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  FERR    (output) REAL array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) REAL array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) REAL array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Internal Parameters
*  ===================
*
*  ITMAX is the maximum number of steps of iterative refinement.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      REAL               THREE
      PARAMETER          ( THREE = 3.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      CHARACTER          TRANSN, TRANST
      INTEGER            COUNT, I, J, KASE, NZ
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGTTRS, SLACON, SLAGTM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGTRFS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
         TRANSN = 'N'
         TRANST = 'T'
      ELSE
         TRANSN = 'T'
         TRANST = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = 4
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 110 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL SCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL SLAGTM( TRANS, N, 1, -ONE, DL, D, DU, X( 1, J ), LDX, ONE,
     $                WORK( N+1 ), N )
*
*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward
*        error bound.
*
         IF( NOTRAN ) THEN
            IF( N.EQ.1 ) THEN
               WORK( 1 ) = ABS( B( 1, J ) ) + ABS( D( 1 )*X( 1, J ) )
            ELSE
               WORK( 1 ) = ABS( B( 1, J ) ) + ABS( D( 1 )*X( 1, J ) ) +
     $                     ABS( DU( 1 )*X( 2, J ) )
               DO 30 I = 2, N - 1
                  WORK( I ) = ABS( B( I, J ) ) +
     $                        ABS( DL( I-1 )*X( I-1, J ) ) +
     $                        ABS( D( I )*X( I, J ) ) +
     $                        ABS( DU( I )*X( I+1, J ) )
   30          CONTINUE
               WORK( N ) = ABS( B( N, J ) ) +
     $                     ABS( DL( N-1 )*X( N-1, J ) ) +
     $                     ABS( D( N )*X( N, J ) )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               WORK( 1 ) = ABS( B( 1, J ) ) + ABS( D( 1 )*X( 1, J ) )
            ELSE
               WORK( 1 ) = ABS( B( 1, J ) ) + ABS( D( 1 )*X( 1, J ) ) +
     $                     ABS( DL( 1 )*X( 2, J ) )
               DO 40 I = 2, N - 1
                  WORK( I ) = ABS( B( I, J ) ) +
     $                        ABS( DU( I-1 )*X( I-1, J ) ) +
     $                        ABS( D( I )*X( I, J ) ) +
     $                        ABS( DL( I )*X( I+1, J ) )
   40          CONTINUE
               WORK( N ) = ABS( B( N, J ) ) +
     $                     ABS( DU( N-1 )*X( N-1, J ) ) +
     $                     ABS( D( N )*X( N, J ) )
            END IF
         END IF
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         S = ZERO
         DO 50 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) /
     $             ( WORK( I )+SAFE1 ) )
            END IF
   50    CONTINUE
         BERR( J ) = S
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND.
     $       COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL SGTTRS( TRANS, N, 1, DLF, DF, DUF, DU2, IPIV,
     $                   WORK( N+1 ), N, INFO )
            CALL SAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use SLACON to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 60 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   60    CONTINUE
*
         KASE = 0
   70    CONTINUE
         CALL SLACON( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ),
     $                KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)**T).
*
               CALL SGTTRS( TRANST, N, 1, DLF, DF, DUF, DU2, IPIV,
     $                      WORK( N+1 ), N, INFO )
               DO 80 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
   80          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 90 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
   90          CONTINUE
               CALL SGTTRS( TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV,
     $                      WORK( N+1 ), N, INFO )
            END IF
            GO TO 70
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 100 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  100    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  110 CONTINUE
*
      RETURN
*
*     End of SGTRFS
*
      END
