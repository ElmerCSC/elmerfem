      SUBROUTINE SLARRB( N, D, L, LD, LLD, IFIRST, ILAST, SIGMA, RELTOL,
     $                   W, WGAP, WERR, WORK, IWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N
      REAL               RELTOL, SIGMA
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), L( * ), LD( * ), LLD( * ), W( * ),
     $                   WERR( * ), WGAP( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given the relatively robust representation(RRR) L D L^T, SLARRB
*  does ``limited'' bisection to locate the eigenvalues of L D L^T,
*  W( IFIRST ) thru' W( ILAST ), to more accuracy. Intervals
*  [left, right] are maintained by storing their mid-points and
*  semi-widths in the arrays W and WERR respectively.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) REAL array, dimension (N)
*          The n diagonal elements of the diagonal matrix D.
*
*  L       (input) REAL array, dimension (N-1)
*          The n-1 subdiagonal elements of the unit bidiagonal matrix L.
*
*  LD      (input) REAL array, dimension (N-1)
*          The n-1 elements L(i)*D(i).
*
*  LLD     (input) REAL array, dimension (N-1)
*          The n-1 elements L(i)*L(i)*D(i).
*
*  IFIRST  (input) INTEGER
*          The index of the first eigenvalue in the cluster.
*
*  ILAST   (input) INTEGER
*          The index of the last eigenvalue in the cluster.
*
*  SIGMA   (input) REAL
*          The shift used to form L D L^T (see SLARRF).
*
*  RELTOL  (input) REAL
*          The relative tolerance.
*
*  W       (input/output) REAL array, dimension (N)
*          On input, W( IFIRST ) thru' W( ILAST ) are estimates of the
*          corresponding eigenvalues of L D L^T.
*          On output, these estimates are ``refined''.
*
*  WGAP    (input/output) REAL array, dimension (N)
*          The gaps between the eigenvalues of L D L^T. Very small
*          gaps are changed on output.
*
*  WERR    (input/output) REAL array, dimension (N)
*          On input, WERR( IFIRST ) thru' WERR( ILAST ) are the errors
*          in the estimates W( IFIRST ) thru' W( ILAST ).
*          On output, these are the ``refined'' errors.
*
*****Reminder to Inder --- WORK is never used in this subroutine *****
*  WORK    (input) REAL array, dimension (???)
*          Workspace.
*
*  IWORK   (input) INTEGER array, dimension (2*N)
*          Workspace.
*
*****Reminder to Inder --- INFO is never set in this subroutine ******
*  INFO    (output) INTEGER
*          Error flag.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0E0, TWO = 2.0E0, HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CNT, I, I1, I2, INITI1, INITI2, J, K, NCNVRG,
     $                   NEIG, NINT, NRIGHT, OLNINT
      REAL               DELTA, EPS, GAP, LEFT, MID, PERT, RIGHT, S,
     $                   THRESH, TMP, WIDTH
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Precision' )
      I1 = IFIRST
      I2 = IFIRST
      NEIG = ILAST - IFIRST + 1
      NCNVRG = 0
      THRESH = RELTOL
      DO 10 I = IFIRST, ILAST
         IWORK( I ) = 0
         PERT = EPS*( ABS( SIGMA )+ABS( W( I ) ) )
         WERR( I ) = WERR( I ) + PERT
         IF( WGAP( I ).LT.PERT )
     $      WGAP( I ) = PERT
   10 CONTINUE
      DO 20 I = I1, ILAST
         IF( I.EQ.1 ) THEN
            GAP = WGAP( I )
         ELSE IF( I.EQ.N ) THEN
            GAP = WGAP( I-1 )
         ELSE
            GAP = MIN( WGAP( I-1 ), WGAP( I ) )
         END IF
         IF( WERR( I ).LT.THRESH*GAP ) THEN
            NCNVRG = NCNVRG + 1
            IWORK( I ) = 1
            IF( I1.EQ.I )
     $         I1 = I1 + 1
         ELSE
            I2 = I
         END IF
   20 CONTINUE
*
*     Initialize the unconverged intervals.
*
      I = I1
      NINT = 0
      RIGHT = ZERO
   30 CONTINUE
      IF( I.LE.I2 ) THEN
         IF( IWORK( I ).EQ.0 ) THEN
            DELTA = EPS
            LEFT = W( I ) - WERR( I )
*
*           Do while( CNT(LEFT).GT.I-1 )
*
   40       CONTINUE
            IF( I.GT.I1 .AND. LEFT.LE.RIGHT ) THEN
               LEFT = RIGHT
               CNT = I - 1
            ELSE
               S = -LEFT
               CNT = 0
               DO 50 J = 1, N - 1
                  TMP = D( J ) + S
                  S = S*( LD( J ) / TMP )*L( J ) - LEFT
                  IF( TMP.LT.ZERO )
     $               CNT = CNT + 1
   50          CONTINUE
               TMP = D( N ) + S
               IF( TMP.LT.ZERO )
     $            CNT = CNT + 1
               IF( CNT.GT.I-1 ) THEN
                  DELTA = TWO*DELTA
                  LEFT = LEFT - ( ABS( SIGMA )+ABS( LEFT ) )*DELTA
                  GO TO 40
               END IF
            END IF
            DELTA = EPS
            RIGHT = W( I ) + WERR( I )
*
*           Do while( CNT(RIGHT).LT.I )
*
   60       CONTINUE
            S = -RIGHT
            CNT = 0
            DO 70 J = 1, N - 1
               TMP = D( J ) + S
               S = S*( LD( J ) / TMP )*L( J ) - RIGHT
               IF( TMP.LT.ZERO )
     $            CNT = CNT + 1
   70       CONTINUE
            TMP = D( N ) + S
            IF( TMP.LT.ZERO )
     $         CNT = CNT + 1
            IF( CNT.LT.I ) THEN
               DELTA = TWO*DELTA
               RIGHT = RIGHT + ( ABS( SIGMA )+ABS( RIGHT ) )*DELTA
               GO TO 60
            END IF
            WERR( I ) = LEFT
            W( I ) = RIGHT
            IWORK( N+I ) = CNT
            NINT = NINT + 1
            I = CNT + 1
         ELSE
            I = I + 1
         END IF
         GO TO 30
      END IF
*
*     While( NCNVRG.LT.NEIG )
*
      INITI1 = I1
      INITI2 = I2
   80 CONTINUE
      IF( NCNVRG.LT.NEIG ) THEN
         OLNINT = NINT
         I = I1
         DO 100 K = 1, OLNINT
            NRIGHT = IWORK( N+I )
            IF( IWORK( I ).EQ.0 ) THEN
               MID = HALF*( WERR( I )+W( I ) )
               S = -MID
               CNT = 0
               DO 90 J = 1, N - 1
                  TMP = D( J ) + S
                  S = S*( LD( J ) / TMP )*L( J ) - MID
                  IF( TMP.LT.ZERO )
     $               CNT = CNT + 1
   90          CONTINUE
               TMP = D( N ) + S
               IF( TMP.LT.ZERO )
     $            CNT = CNT + 1
               CNT = MAX( I-1, MIN( NRIGHT, CNT ) )
               IF( I.EQ.NRIGHT ) THEN
                  IF( I.EQ.IFIRST ) THEN
                     GAP = WERR( I+1 ) - W( I )
                  ELSE IF( I.EQ.ILAST ) THEN
                     GAP = WERR( I ) - W( I-1 )
                  ELSE
                     GAP = MIN( WERR( I+1 )-W( I ), WERR( I )-W( I-1 ) )
                  END IF
                  WIDTH = W( I ) - MID
                  IF( WIDTH.LT.THRESH*GAP ) THEN
                     NCNVRG = NCNVRG + 1
                     IWORK( I ) = 1
                     IF( I1.EQ.I ) THEN
                        I1 = I1 + 1
                        NINT = NINT - 1
                     END IF
                  END IF
               END IF
               IF( IWORK( I ).EQ.0 )
     $            I2 = K
               IF( CNT.EQ.I-1 ) THEN
                  WERR( I ) = MID
               ELSE IF( CNT.EQ.NRIGHT ) THEN
                  W( I ) = MID
               ELSE
                  IWORK( N+I ) = CNT
                  NINT = NINT + 1
                  WERR( CNT+1 ) = MID
                  W( CNT+1 ) = W( I )
                  W( I ) = MID
                  I = CNT + 1
                  IWORK( N+I ) = NRIGHT
               END IF
            END IF
            I = NRIGHT + 1
  100    CONTINUE
         NINT = NINT - OLNINT + I2
         GO TO 80
      END IF
      DO 110 I = INITI1, INITI2
         W( I ) = HALF*( WERR( I )+W( I ) )
         WERR( I ) = W( I ) - WERR( I )
  110 CONTINUE
*
      RETURN
*
*     End of SLARRB
*
      END
