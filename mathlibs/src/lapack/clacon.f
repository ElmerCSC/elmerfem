      SUBROUTINE CLACON( N, V, X, EST, KASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
*     ..
*     .. Array Arguments ..
      COMPLEX            V( N ), X( N )
*     ..
*
*  Purpose
*  =======
*
*  CLACON estimates the 1-norm of a square, complex matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) COMPLEX array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) COMPLEX array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         where A' is the conjugate transpose of A, and CLACON must be
*         re-called with all the other parameters unchanged.
*
*  EST    (output) REAL
*         An estimate (a lower bound) for norm(A).
*
*  KASE   (input/output) INTEGER
*         On the initial call to CLACON, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from CLACON, KASE will again be 0.
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named CONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  Last modified:  April, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E0, TWO = 2.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),
     $                   CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
*     ..
*     .. External Functions ..
      INTEGER            ICMAX1
      REAL               SCSUM1, SLAMCH
      EXTERNAL           ICMAX1, SCSUM1, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 90, 120 )JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 130
      END IF
      EST = SCSUM1( N, X, 1 )
*
      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,
     $               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
*
   40 CONTINUE
      J = ICMAX1( N, X, 1 )
      ITER = 2
*
*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN
*
*     ................ ENTRY   (JUMP = 3)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
   70 CONTINUE
      CALL CCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SCSUM1( N, V, 1 )
*
*     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )
     $   GO TO 100
*
      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,
     $               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
*
*     ................ ENTRY   (JUMP = 4)
*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
*
   90 CONTINUE
      JLAST = J
      J = ICMAX1( N, X, 1 )
      IF( ( ABS( X( JLAST ) ).NE.ABS( X( J ) ) ) .AND.
     $    ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
*
*     ITERATION COMPLETE.  FINAL STAGE.
*
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
*
*     ................ ENTRY   (JUMP = 5)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
  120 CONTINUE
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL CCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*
  130 CONTINUE
      KASE = 0
      RETURN
*
*     End of CLACON
*
      END
