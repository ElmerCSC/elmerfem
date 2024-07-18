      SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N
      REAL               SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      COMPLEX            A( LDA, * ), RHS( * )
*     ..
*
*  Purpose
*  =======
*
*  CGESC2 solves a system of linear equations
*
*            A * X = scale* RHS
*
*  with a general N-by-N matrix A using the LU factorization with
*  complete pivoting computed by CGETC2.
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input) COMPLEX array, dimension (LDA, N)
*          On entry, the  LU part of the factorization of the n-by-n
*          matrix A computed by CGETC2:  A = P * L * U * Q
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1, N).
*
*  RHS     (input/output) COMPLEX array, dimension N.
*          On entry, the right hand side vector b.
*          On exit, the solution vector X.
*
*  IPIV    (iput) INTEGER array, dimension (N).
*          The pivot indices; for 1 <= i <= N, row i of the
*          matrix has been interchanged with row IPIV(i).
*
*  JPIV    (iput) INTEGER array, dimension (N).
*          The pivot indices; for 1 <= j <= N, column j of the
*          matrix has been interchanged with column JPIV(j).
*
*  SCALE    (output) REAL
*           On exit, SCALE contains the scale factor. SCALE is chosen
*           0 <= SCALE <= 1 to prevent owerflow in the solution.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               BIGNUM, EPS, SMLNUM
      COMPLEX            TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASWP, CSCAL, SLABAD
*     ..
*     .. External Functions ..
      INTEGER            ICAMAX
      REAL               SLAMCH
      EXTERNAL           ICAMAX, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
*     Set constant to control overflow
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Apply permutations IPIV to RHS
*
      CALL CLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
*
*     Solve for L part
*
      DO 20 I = 1, N - 1
         DO 10 J = I + 1, N
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
   10    CONTINUE
   20 CONTINUE
*
*     Solve for U part
*
      SCALE = ONE
*
*     Check for scaling
*
      I = ICAMAX( N, RHS, 1 )
      IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
         TEMP = CMPLX( ONE / TWO, ZERO ) / ABS( RHS( I ) )
         CALL CSCAL( N, TEMP, RHS( 1 ), 1 )
         SCALE = SCALE*REAL( TEMP )
      END IF
      DO 40 I = N, 1, -1
         TEMP = CMPLX( ONE, ZERO ) / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         DO 30 J = I + 1, N
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
   30    CONTINUE
   40 CONTINUE
*
*     Apply permutations JPIV to the solution (RHS)
*
      CALL CLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      RETURN
*
*     End of CGESC2
*
      END
