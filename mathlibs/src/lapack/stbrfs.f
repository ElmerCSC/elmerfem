      SUBROUTINE STBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,
     $                   LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               AB( LDAB, * ), B( LDB, * ), BERR( * ),
     $                   FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  STBRFS provides error bounds and backward error estimates for the
*  solution to a system of linear equations with a triangular band
*  coefficient matrix.
*
*  The solution matrix X must be computed by STBTRS or some other
*  means before entering this routine.  STBRFS does not do iterative
*  refinement because doing so cannot improve the backward error.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals or subdiagonals of the
*          triangular band matrix A.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  AB      (input) REAL array, dimension (LDAB,N)
*          The upper or lower triangular band matrix A, stored in the
*          first kd+1 rows of the array. The j-th column of A is stored
*          in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          If DIAG = 'U', the diagonal elements of A are not referenced
*          and are assumed to be 1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input) REAL array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input) REAL array, dimension (LDX,NRHS)
*          The solution matrix X.
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
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      CHARACTER          TRANST
      INTEGER            I, J, K, KASE, NZ
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SLACON, STBMV, STBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
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
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STBRFS', -INFO )
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
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = KD + 2
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 250 J = 1, NRHS
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A or A', depending on TRANS.
*
         CALL SCOPY( N, X( 1, J ), 1, WORK( N+1 ), 1 )
         CALL STBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK( N+1 ),
     $               1 )
         CALL SAXPY( N, -ONE, B( 1, J ), 1, WORK( N+1 ), 1 )
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
         DO 20 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   20    CONTINUE
*
         IF( NOTRAN ) THEN
*
*           Compute abs(A)*abs(X) + abs(B).
*
            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 40 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 30 I = MAX( 1, K-KD ), K
                        WORK( I ) = WORK( I ) +
     $                              ABS( AB( KD+1+I-K, K ) )*XK
   30                CONTINUE
   40             CONTINUE
               ELSE
                  DO 60 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 50 I = MAX( 1, K-KD ), K - 1
                        WORK( I ) = WORK( I ) +
     $                              ABS( AB( KD+1+I-K, K ) )*XK
   50                CONTINUE
                     WORK( K ) = WORK( K ) + XK
   60             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 80 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 70 I = K, MIN( N, K+KD )
                        WORK( I ) = WORK( I ) + ABS( AB( 1+I-K, K ) )*XK
   70                CONTINUE
   80             CONTINUE
               ELSE
                  DO 100 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 90 I = K + 1, MIN( N, K+KD )
                        WORK( I ) = WORK( I ) + ABS( AB( 1+I-K, K ) )*XK
   90                CONTINUE
                     WORK( K ) = WORK( K ) + XK
  100             CONTINUE
               END IF
            END IF
         ELSE
*
*           Compute abs(A')*abs(X) + abs(B).
*
            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 120 K = 1, N
                     S = ZERO
                     DO 110 I = MAX( 1, K-KD ), K
                        S = S + ABS( AB( KD+1+I-K, K ) )*
     $                      ABS( X( I, J ) )
  110                CONTINUE
                     WORK( K ) = WORK( K ) + S
  120             CONTINUE
               ELSE
                  DO 140 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 130 I = MAX( 1, K-KD ), K - 1
                        S = S + ABS( AB( KD+1+I-K, K ) )*
     $                      ABS( X( I, J ) )
  130                CONTINUE
                     WORK( K ) = WORK( K ) + S
  140             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 160 K = 1, N
                     S = ZERO
                     DO 150 I = K, MIN( N, K+KD )
                        S = S + ABS( AB( 1+I-K, K ) )*ABS( X( I, J ) )
  150                CONTINUE
                     WORK( K ) = WORK( K ) + S
  160             CONTINUE
               ELSE
                  DO 180 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 170 I = K + 1, MIN( N, K+KD )
                        S = S + ABS( AB( 1+I-K, K ) )*ABS( X( I, J ) )
  170                CONTINUE
                     WORK( K ) = WORK( K ) + S
  180             CONTINUE
               END IF
            END IF
         END IF
         S = ZERO
         DO 190 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) /
     $             ( WORK( I )+SAFE1 ) )
            END IF
  190    CONTINUE
         BERR( J ) = S
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
         DO 200 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
  200    CONTINUE
*
         KASE = 0
  210    CONTINUE
         CALL SLACON( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ),
     $                KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)').
*
               CALL STBSV( UPLO, TRANST, DIAG, N, KD, AB, LDAB,
     $                     WORK( N+1 ), 1 )
               DO 220 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  220          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 230 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  230          CONTINUE
               CALL STBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB,
     $                     WORK( N+1 ), 1 )
            END IF
            GO TO 210
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 240 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  240    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  250 CONTINUE
*
      RETURN
*
*     End of STBRFS
*
      END
