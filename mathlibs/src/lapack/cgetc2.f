      SUBROUTINE CGETC2( N, A, LDA, IPIV, JPIV, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CGETC2 computes an LU factorization, using complete pivoting, of the
*  n-by-n matrix A. The factorization has the form A = P * L * U * Q,
*  where P and Q are permutation matrices, L is lower triangular with
*  unit diagonal elements and U is upper triangular.
*
*  This is a level 1 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the n-by-n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*          If U(k, k) appears to be less than SMIN, U(k, k) is given the
*          value of SMIN, giving a nonsingular perturbed system.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1, N).
*
*  IPIV    (output) INTEGER array, dimension (N).
*          The pivot indices; for 1 <= i <= N, row i of the
*          matrix has been interchanged with row IPIV(i).
*
*  JPIV    (output) INTEGER array, dimension (N).
*          The pivot indices; for 1 <= j <= N, column j of the
*          matrix has been interchanged with column JPIV(j).
*
*  INFO    (output) INTEGER
*           = 0: successful exit
*           > 0: if INFO = k, U(k, k) is likely to produce overflow if
*                one tries to solve for x in Ax = b. So U is perturbed
*                to avoid the overflow.
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
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IP, IPV, J, JP, JPV
      REAL               BIGNUM, EPS, SMIN, SMLNUM, XMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGERU, CSWAP, SLABAD
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Set constants to control overflow
*
      INFO = 0
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Factorize A using complete pivoting.
*     Set pivots less than SMIN to SMIN
*
      DO 40 I = 1, N - 1
*
*        Find max element in matrix A
*
         XMAX = ZERO
         DO 20 IP = I, N
            DO 10 JP = I, N
               IF( ABS( A( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( A( IP, JP ) )
                  IPV = IP
                  JPV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( I.EQ.1 )
     $      SMIN = MAX( EPS*XMAX, SMLNUM )
*
*        Swap rows
*
         IF( IPV.NE.I )
     $      CALL CSWAP( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA )
         IPIV( I ) = IPV
*
*        Swap columns
*
         IF( JPV.NE.I )
     $      CALL CSWAP( N, A( 1, JPV ), 1, A( 1, I ), 1 )
         JPIV( I ) = JPV
*
*        Check for singularity
*
         IF( ABS( A( I, I ) ).LT.SMIN ) THEN
            INFO = I
            A( I, I ) = CMPLX( SMIN, ZERO )
         END IF
         DO 30 J = I + 1, N
            A( J, I ) = A( J, I ) / A( I, I )
   30    CONTINUE
         CALL CGERU( N-I, N-I, -CMPLX( ONE ), A( I+1, I ), 1,
     $               A( I, I+1 ), LDA, A( I+1, I+1 ), LDA )
   40 CONTINUE
*
      IF( ABS( A( N, N ) ).LT.SMIN ) THEN
         INFO = N
         A( N, N ) = CMPLX( SMIN, ZERO )
      END IF
      RETURN
*
*     End of CGETC2
*
      END
