      SUBROUTINE CGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGEHD2 reduces a complex general matrix A to upper Hessenberg form H
*  by a unitary similarity transformation:  Q' * A * Q = H .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to CGEBAL; otherwise they should be
*          set to 1 and N respectively. See Further Details.
*          1 <= ILO <= IHI <= max(1,N).
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the n by n general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the unitary matrix Q as a product of elementary
*          reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  TAU     (output) COMPLEX array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of (ihi-ilo) elementary
*  reflectors
*
*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*  exit in A(i+2:ihi,i), and tau in TAU(i).
*
*  The contents of A are illustrated by the following example, with
*  n = 7, ilo = 2 and ihi = 6:
*
*  on entry,                        on exit,
*
*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*  (                         a )    (                          a )
*
*  where a denotes an element of the original matrix A, h denotes a
*  modified element of the upper Hessenberg matrix H, and vi denotes an
*  element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX            ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF, CLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEHD2', -INFO )
         RETURN
      END IF
*
      DO 10 I = ILO, IHI - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         ALPHA = A( I+1, I )
         CALL CLARFG( IHI-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
         A( I+1, I ) = ONE
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL CLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),
     $               A( 1, I+1 ), LDA, WORK )
*
*        Apply H(i)' to A(i+1:ihi,i+1:n) from the left
*
         CALL CLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1,
     $               CONJG( TAU( I ) ), A( I+1, I+1 ), LDA, WORK )
*
         A( I+1, I ) = ALPHA
   10 CONTINUE
*
      RETURN
*
*     End of CGEHD2
*
      END
