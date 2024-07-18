      SUBROUTINE CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, J1, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
*  in an upper triangular matrix pair (A, B) by an unitary equivalence
*  transformation.
*
*  (A, B) must be in generalized Schur canonical form, that is, A and
*  B are both upper triangular.
*
*  Optionally, the matrices Q and Z of generalized Schur vectors are
*  updated.
*
*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'
*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'
*
*
*  Arguments
*  =========
*
*  WANTQ   (input) LOGICAL
*          .TRUE. : update the left transformation matrix Q;
*          .FALSE.: do not update Q.
*
*  WANTZ   (input) LOGICAL
*          .TRUE. : update the right transformation matrix Z;
*          .FALSE.: do not update Z.
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  A       (input/output) COMPLEX arrays, dimensions (LDA,N)
*          On entry, the matrix A in the pair (A, B).
*          On exit, the updated matrix A.
*
*  LDA     (input)  INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) COMPLEX arrays, dimensions (LDB,N)
*          On entry, the matrix B in the pair (A, B).
*          On exit, the updated matrix B.
*
*  LDB     (input)  INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  Q       (input/output) COMPLEX array, dimension (LDZ,N)
*          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,
*          the updated matrix Q.
*          Not referenced if WANTQ = .FALSE..
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= 1;
*          If WANTQ = .TRUE., LDQ >= N.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ,N)
*          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,
*          the updated matrix Z.
*          Not referenced if WANTZ = .FALSE..
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= 1;
*          If WANTZ = .TRUE., LDZ >= N.
*
*  J1      (input) INTEGER
*          The index to the first block (A11, B11).
*
*  INFO    (output) INTEGER
*           =0:  Successful exit.
*           =1:  The transformed matrix pair (A, B) would be too far
*                from generalized Schur form; the problem is ill-
*                conditioned. (A, B) may have been partially reordered,
*                and ILST points to the first row of the current
*                position of the block being moved.
*
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  In the current code both weak and strong stability tests are
*  performed. The user can omit the strong stability test by changing
*  the internal logical parameter WANDS to .FALSE.. See ref. [2] for
*  details.
*
*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*
*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,
*      Department of Computing Science, Umea University, S-901 87 Umea,
*      Sweden, 1994. Also as LAPACK Working Note 87. To appear in
*      Numerical Algorithms, 1996.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               TEN
      PARAMETER          ( TEN = 10.0E+0 )
      INTEGER            LDST
      PARAMETER          ( LDST = 2 )
      LOGICAL            WANDS
      PARAMETER          ( WANDS = .TRUE. )
*     ..
*     .. Local Scalars ..
      LOGICAL            STRONG, WEAK
      INTEGER            I, M
      REAL               CQ, CZ, EPS, SA, SB, SCALE, SMLNUM, SS, SUM,
     $                   THRESH, WS
      COMPLEX            CDUM, F, G, SQ, SZ
*     ..
*     .. Local Arrays ..
      COMPLEX            S( LDST, LDST ), T( LDST, LDST ), WORK( 8 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CLARTG, CLASSQ, CROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      M = LDST
      WEAK = .FALSE.
      STRONG = .FALSE.
*
*     Make a local copy of selected block in (A, B)
*
      CALL CLACPY( 'Full', M, M, A( J1, J1 ), LDA, S, LDST )
      CALL CLACPY( 'Full', M, M, B( J1, J1 ), LDB, T, LDST )
*
*     Compute the threshold for testing the acceptance of swapping.
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      SCALE = REAL( CZERO )
      SUM = REAL( CONE )
      CALL CLACPY( 'Full', M, M, S, LDST, WORK, M )
      CALL CLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
      CALL CLASSQ( 2*M*M, WORK, 1, SCALE, SUM )
      SA = SCALE*SQRT( SUM )
      THRESH = MAX( TEN*EPS*SA, SMLNUM )
*
*     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
*     using Givens rotations and perform the swap tentatively.
*
      F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
      G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
      SA = ABS( S( 2, 2 ) )
      SB = ABS( T( 2, 2 ) )
      CALL CLARTG( G, F, CZ, SZ, CDUM )
      SZ = -SZ
      CALL CROT( 2, S( 1, 1 ), 1, S( 1, 2 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( 2, T( 1, 1 ), 1, T( 1, 2 ), 1, CZ, CONJG( SZ ) )
      IF( SA.GE.SB ) THEN
         CALL CLARTG( S( 1, 1 ), S( 2, 1 ), CQ, SQ, CDUM )
      ELSE
         CALL CLARTG( T( 1, 1 ), T( 2, 1 ), CQ, SQ, CDUM )
      END IF
      CALL CROT( 2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, CQ, SQ )
      CALL CROT( 2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, CQ, SQ )
*
*     Weak stability test: |S21| + |T21| <= O(EPS F-norm((S, T)))
*
      WS = ABS( S( 2, 1 ) ) + ABS( T( 2, 1 ) )
      WEAK = WS.LE.THRESH
      IF( .NOT.WEAK )
     $   GO TO 20
*
      IF( WANDS ) THEN
*
*        Strong stability test:
*           F-norm((A-QL'*S*QR, B-QL'*T*QR)) <= O(EPS*F-norm((A, B)))
*
         CALL CLACPY( 'Full', M, M, S, LDST, WORK, M )
         CALL CLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
         CALL CROT( 2, WORK, 1, WORK( 3 ), 1, CZ, -CONJG( SZ ) )
         CALL CROT( 2, WORK( 5 ), 1, WORK( 7 ), 1, CZ, -CONJG( SZ ) )
         CALL CROT( 2, WORK, 2, WORK( 2 ), 2, CQ, -SQ )
         CALL CROT( 2, WORK( 5 ), 2, WORK( 6 ), 2, CQ, -SQ )
         DO 10 I = 1, 2
            WORK( I ) = WORK( I ) - A( J1+I-1, J1 )
            WORK( I+2 ) = WORK( I+2 ) - A( J1+I-1, J1+1 )
            WORK( I+4 ) = WORK( I+4 ) - B( J1+I-1, J1 )
            WORK( I+6 ) = WORK( I+6 ) - B( J1+I-1, J1+1 )
   10    CONTINUE
         SCALE = REAL( CZERO )
         SUM = REAL( CONE )
         CALL CLASSQ( 2*M*M, WORK, 1, SCALE, SUM )
         SS = SCALE*SQRT( SUM )
         STRONG = SS.LE.THRESH
         IF( .NOT.STRONG )
     $      GO TO 20
      END IF
*
*     If the swap is accepted ("weakly" and "strongly"), apply the
*     equivalence transformations to the original matrix pair (A,B)
*
      CALL CROT( J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      CALL CROT( N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, CQ, SQ )
      CALL CROT( N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, CQ, SQ )
*
*     Set  N1 by N2 (2,1) blocks to 0
*
      A( J1+1, J1 ) = CZERO
      B( J1+1, J1 ) = CZERO
*
*     Accumulate transformations into Q and Z if requested.
*
      IF( WANTZ )
     $   CALL CROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, CZ, CONJG( SZ ) )
      IF( WANTQ )
     $   CALL CROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, CQ, CONJG( SQ ) )
*
*     Exit with INFO = 0 if swap was successfully performed.
*
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   20 CONTINUE
      INFO = 1
      RETURN
*
*     End of CTGEX2
*
      END
