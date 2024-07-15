      SUBROUTINE ZLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
     $                   WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N, OFFSET
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      DOUBLE PRECISION   VN1( * ), VN2( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAQP2 computes a QR factorization with column pivoting of
*  the block A(OFFSET+1:M,1:N).
*  The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A. N >= 0.
*
*  OFFSET  (input) INTEGER
*          The number of rows of the matrix A that must be pivoted
*          but no factorized. OFFSET >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is
*          the triangular factor obtained; the elements in block
*          A(OFFSET+1:M,1:N) below the diagonal, together with the
*          array TAU, represent the orthogonal matrix Q as a product of
*          elementary reflectors. Block A(1:OFFSET,1:N) has been
*          accordingly pivoted, but no factorized.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  JPVT    (input/output) INTEGER array, dimension (N)
*          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(i) = 0,
*          the i-th column of A is a free column.
*          On exit, if JPVT(i) = k, then the i-th column of A*P
*          was the k-th column of A.
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
*
*  VN1     (input/output) DOUBLE PRECISION array, dimension (N)
*          The vector with the partial column norms.
*
*  VN2     (input/output) DOUBLE PRECISION array, dimension (N)
*          The vector with the exact column norms.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
*    X. Sun, Computer Science Dept., Duke University, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      COMPLEX*16         CONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MN, OFFPI, PVT
      DOUBLE PRECISION   TEMP, TEMP2
      COMPLEX*16         AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLARF, ZLARFG, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DZNRM2
      EXTERNAL           IDAMAX, DZNRM2
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M-OFFSET, N )
*
*     Compute factorization.
*
      DO 20 I = 1, MN
*
         OFFPI = OFFSET + I
*
*        Determine ith pivot column and swap if necessary.
*
         PVT = ( I-1 ) + IDAMAX( N-I+1, VN1( I ), 1 )
*
         IF( PVT.NE.I ) THEN
            CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            VN1( PVT ) = VN1( I )
            VN2( PVT ) = VN2( I )
         END IF
*
*        Generate elementary reflector H(i).
*
         IF( OFFPI.LT.M ) THEN
            CALL ZLARFG( M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ), 1,
     $                   TAU( I ) )
         ELSE
            CALL ZLARFG( 1, A( M, I ), A( M, I ), 1, TAU( I ) )
         END IF
*
         IF( I.LT.N ) THEN
*
*           Apply H(i)' to A(offset+i:m,i+1:n) from the left.
*
            AII = A( OFFPI, I )
            A( OFFPI, I ) = CONE
            CALL ZLARF( 'Left', M-OFFPI+1, N-I, A( OFFPI, I ), 1,
     $                  DCONJG( TAU( I ) ), A( OFFPI, I+1 ), LDA,
     $                  WORK( 1 ) )
            A( OFFPI, I ) = AII
         END IF
*
*        Update partial column norms.
*
         DO 10 J = I + 1, N
            IF( VN1( J ).NE.ZERO ) THEN
               TEMP = ONE - ( ABS( A( OFFPI, J ) ) / VN1( J ) )**2
               TEMP = MAX( TEMP, ZERO )
               TEMP2 = ONE + 0.05D0*TEMP*( VN1( J ) / VN2( J ) )**2
               IF( TEMP2.EQ.ONE ) THEN
                  IF( OFFPI.LT.M ) THEN
                     VN1( J ) = DZNRM2( M-OFFPI, A( OFFPI+1, J ), 1 )
                     VN2( J ) = VN1( J )
                  ELSE
                     VN1( J ) = ZERO
                     VN2( J ) = ZERO
                  END IF
               ELSE
                  VN1( J ) = VN1( J )*SQRT( TEMP )
               END IF
            END IF
   10    CONTINUE
*
   20 CONTINUE
*
      RETURN
*
*     End of ZLAQP2
*
      END
