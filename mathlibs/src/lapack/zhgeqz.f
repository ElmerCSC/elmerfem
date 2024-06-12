      SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,
     $                   ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
     $                   RWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZHGEQZ implements a single-shift version of the QZ
*  method for finding the generalized eigenvalues w(i)=ALPHA(i)/BETA(i)
*  of the equation
*
*       det( A - w(i) B ) = 0
*
*  If JOB='S', then the pair (A,B) is simultaneously
*  reduced to Schur form (i.e., A and B are both upper triangular) by
*  applying one unitary tranformation (usually called Q) on the left and
*  another (usually called Z) on the right.  The diagonal elements of
*  A are then ALPHA(1),...,ALPHA(N), and of B are BETA(1),...,BETA(N).
*
*  If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the unitary
*  transformations used to reduce (A,B) are accumulated into the arrays
*  Q and Z s.t.:
*
*       Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*
*       Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': compute only ALPHA and BETA.  A and B will not
*                 necessarily be put into generalized Schur form.
*          = 'S': put A and B into generalized Schur form, as well
*                 as computing ALPHA and BETA.
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not modify Q.
*          = 'V': multiply the array Q on the right by the conjugate
*                 transpose of the unitary tranformation that is
*                 applied to the left side of A and B to reduce them
*                 to Schur form.
*          = 'I': like COMPQ='V', except that Q will be initialized to
*                 the identity first.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not modify Z.
*          = 'V': multiply the array Z on the right by the unitary
*                 tranformation that is applied to the right side of
*                 A and B to reduce them to Schur form.
*          = 'I': like COMPZ='V', except that Z will be initialized to
*                 the identity first.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows and
*          columns 1:ILO-1 and IHI+1:N.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the N-by-N upper Hessenberg matrix A.  Elements
*          below the subdiagonal must be zero.
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to upper triangular form.
*          If JOB='E', then on exit A will have been destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max( 1, N ).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.  Elements
*          below the diagonal must be zero.
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to upper triangular form.
*          If JOB='E', then on exit B will have been destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max( 1, N ).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*          The diagonal elements of A when the pair (A,B) has been
*          reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N
*          are the generalized eigenvalues.
*
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          The diagonal elements of B when the pair (A,B) has been
*          reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N
*          are the generalized eigenvalues.  A and B are normalized
*          so that BETA(1),...,BETA(N) are non-negative real numbers.
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
*          If COMPQ='N', then Q will not be referenced.
*          If COMPQ='V' or 'I', then the conjugate transpose of the
*             unitary transformations which are applied to A and B on
*             the left will be applied to the array Q on the right.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          If COMPZ='N', then Z will not be referenced.
*          If COMPZ='V' or 'I', then the unitary transformations which
*             are applied to A and B on the right will be applied to the
*             array Z on the right.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (A,B) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (A,B) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO-N+1,...,N should be correct.
*          > 2*N:     various "impossible" errors.
*
*  Further Details
*  ===============
*
*  We assume that complex ABS works as long as its value is less than
*  overflow.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT
      DOUBLE PRECISION   ABSB, ANORM, ASCALE, ATOL, BNORM, BSCALE, BTOL,
     $                   C, SAFMIN, TEMP, TEMP2, TEMPR, ULP
      COMPLEX*16         ABI22, AD11, AD12, AD21, AD22, CTEMP, CTEMP2,
     $                   CTEMP3, ESHIFT, RTDISC, S, SHIFT, SIGNBC, T,
     $                   U12, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN,
     $                   SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.N ) THEN
         INFO = -8
      ELSE IF( LDB.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHGEQZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
*     WORK( 1 ) = CMPLX( 1 )
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DCMPLX( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = ZLANHS( 'F', IN, A( ILO, ILO ), LDA, RWORK )
      BNORM = ZLANHS( 'F', IN, B( ILO, ILO ), LDB, RWORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*
*
*     Set Eigenvalues IHI+1:N
*
      DO 10 J = IHI + 1, N
         ABSB = ABS( B( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( J, J ) / ABSB )
            B( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, B( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, A( 1, J ), 1 )
            ELSE
               A( J, J ) = A( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            B( J, J ) = CZERO
         END IF
         ALPHA( J ) = A( J, J )
         BETA( J ) = B( J, J )
   10 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 190
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever
*        Row operations modify columns whatever:ILASTM
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = CZERO
      MAXIT = 30*( IHI-ILO+1 )
*
      DO 170 JITER = 1, MAXIT
*
*        Check for too many iterations.
*
         IF( JITER.GT.MAXIT )
     $      GO TO 180
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: A(j,j-1)=0  or  j=ILO
*           2: B(j,j)=0
*
*        Special case: j=ILAST
*
         IF( ILAST.EQ.ILO ) THEN
            GO TO 60
         ELSE
            IF( ABS1( A( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               A( ILAST, ILAST-1 ) = CZERO
               GO TO 60
            END IF
         END IF
*
         IF( ABS( B( ILAST, ILAST ) ).LE.BTOL ) THEN
            B( ILAST, ILAST ) = CZERO
            GO TO 50
         END IF
*
*        General case: j<ILAST
*
         DO 40 J = ILAST - 1, ILO, -1
*
*           Test 1: for A(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS1( A( J, J-1 ) ).LE.ATOL ) THEN
                  A( J, J-1 ) = CZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for B(j,j)=0
*
            IF( ABS( B( J, J ) ).LT.BTOL ) THEN
               B( J, J ) = CZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  IF( ABS1( A( J, J-1 ) )*( ASCALE*ABS1( A( J+1,
     $                J ) ) ).LE.ABS1( A( J, J ) )*( ASCALE*ATOL ) )
     $                ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 20 JCH = J, ILAST - 1
                     CTEMP = A( JCH, JCH )
                     CALL ZLARTG( CTEMP, A( JCH+1, JCH ), C, S,
     $                            A( JCH, JCH ) )
                     A( JCH+1, JCH ) = CZERO
                     CALL ZROT( ILASTM-JCH, A( JCH, JCH+1 ), LDA,
     $                          A( JCH+1, JCH+1 ), LDA, C, S )
                     CALL ZROT( ILASTM-JCH, B( JCH, JCH+1 ), LDB,
     $                          B( JCH+1, JCH+1 ), LDB, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     IF( ILAZR2 )
     $                  A( JCH, JCH-1 ) = A( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS1( B( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 60
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 70
                        END IF
                     END IF
                     B( JCH+1, JCH+1 ) = CZERO
   20             CONTINUE
                  GO TO 50
               ELSE
*
*                 Only test 2 passed -- chase the zero to B(ILAST,ILAST)
*                 Then process as in the case B(ILAST,ILAST)=0
*
                  DO 30 JCH = J, ILAST - 1
                     CTEMP = B( JCH, JCH+1 )
                     CALL ZLARTG( CTEMP, B( JCH+1, JCH+1 ), C, S,
     $                            B( JCH, JCH+1 ) )
                     B( JCH+1, JCH+1 ) = CZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL ZROT( ILASTM-JCH-1, B( JCH, JCH+2 ), LDB,
     $                             B( JCH+1, JCH+2 ), LDB, C, S )
                     CALL ZROT( ILASTM-JCH+2, A( JCH, JCH-1 ), LDA,
     $                          A( JCH+1, JCH-1 ), LDA, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     CTEMP = A( JCH+1, JCH )
                     CALL ZLARTG( CTEMP, A( JCH+1, JCH-1 ), C, S,
     $                            A( JCH+1, JCH ) )
                     A( JCH+1, JCH-1 ) = CZERO
                     CALL ZROT( JCH+1-IFRSTM, A( IFRSTM, JCH ), 1,
     $                          A( IFRSTM, JCH-1 ), 1, C, S )
                     CALL ZROT( JCH-IFRSTM, B( IFRSTM, JCH ), 1,
     $                          B( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL ZROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   30             CONTINUE
                  GO TO 50
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 70
            END IF
*
*           Neither test passed -- try next J
*
   40    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = 2*N + 1
         GO TO 210
*
*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   50    CONTINUE
         CTEMP = A( ILAST, ILAST )
         CALL ZLARTG( CTEMP, A( ILAST, ILAST-1 ), C, S,
     $                A( ILAST, ILAST ) )
         A( ILAST, ILAST-1 ) = CZERO
         CALL ZROT( ILAST-IFRSTM, A( IFRSTM, ILAST ), 1,
     $              A( IFRSTM, ILAST-1 ), 1, C, S )
         CALL ZROT( ILAST-IFRSTM, B( IFRSTM, ILAST ), 1,
     $              B( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL ZROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
*
   60    CONTINUE
         ABSB = ABS( B( ILAST, ILAST ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( ILAST, ILAST ) / ABSB )
            B( ILAST, ILAST ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( ILAST-IFRSTM, SIGNBC, B( IFRSTM, ILAST ), 1 )
               CALL ZSCAL( ILAST+1-IFRSTM, SIGNBC, A( IFRSTM, ILAST ),
     $                     1 )
            ELSE
               A( ILAST, ILAST ) = A( ILAST, ILAST )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, ILAST ), 1 )
         ELSE
            B( ILAST, ILAST ) = CZERO
         END IF
         ALPHA( ILAST ) = A( ILAST, ILAST )
         BETA( ILAST ) = B( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 190
*
*        Reset counters
*
         IITER = 0
         ESHIFT = CZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 160
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST.  We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
   70    CONTINUE
         IITER = IITER + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute the Shift.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.NE.IITER ) THEN
*
*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
*           the bottom-right 2x2 block of A inv(B) which is nearest to
*           the bottom-right element.
*
*           We factor B as U*D, where U has unit diagonals, and
*           compute (A*inv(D))*inv(U).
*
            U12 = ( BSCALE*B( ILAST-1, ILAST ) ) /
     $            ( BSCALE*B( ILAST, ILAST ) )
            AD11 = ( ASCALE*A( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*A( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*A( ILAST-1, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            AD22 = ( ASCALE*A( ILAST, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            ABI22 = AD22 - U12*AD21
*
            T = HALF*( AD11+ABI22 )
            RTDISC = SQRT( T**2+AD12*AD21-AD11*AD22 )
            TEMP = DBLE( T-ABI22 )*DBLE( RTDISC ) +
     $             DIMAG( T-ABI22 )*DIMAG( RTDISC )
            IF( TEMP.LE.ZERO ) THEN
               SHIFT = T + RTDISC
            ELSE
               SHIFT = T - RTDISC
            END IF
         ELSE
*
*           Exceptional shift.  Chosen for no particularly good reason.
*
            ESHIFT = ESHIFT + DCONJG( ( ASCALE*A( ILAST-1, ILAST ) ) /
     $               ( BSCALE*B( ILAST-1, ILAST-1 ) ) )
            SHIFT = ESHIFT
         END IF
*
*        Now check for two consecutive small subdiagonals.
*
         DO 80 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            CTEMP = ASCALE*A( J, J ) - SHIFT*( BSCALE*B( J, J ) )
            TEMP = ABS1( CTEMP )
            TEMP2 = ASCALE*ABS1( A( J+1, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS1( A( J, J-1 ) )*TEMP2.LE.TEMP*ATOL )
     $         GO TO 90
   80    CONTINUE
*
         ISTART = IFIRST
         CTEMP = ASCALE*A( IFIRST, IFIRST ) -
     $           SHIFT*( BSCALE*B( IFIRST, IFIRST ) )
   90    CONTINUE
*
*        Do an implicit-shift QZ sweep.
*
*        Initial Q
*
         CTEMP2 = ASCALE*A( ISTART+1, ISTART )
         CALL ZLARTG( CTEMP, CTEMP2, C, S, CTEMP3 )
*
*        Sweep
*
         DO 150 J = ISTART, ILAST - 1
            IF( J.GT.ISTART ) THEN
               CTEMP = A( J, J-1 )
               CALL ZLARTG( CTEMP, A( J+1, J-1 ), C, S, A( J, J-1 ) )
               A( J+1, J-1 ) = CZERO
            END IF
*
            DO 100 JC = J, ILASTM
               CTEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -DCONJG( S )*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = CTEMP
               CTEMP2 = C*B( J, JC ) + S*B( J+1, JC )
               B( J+1, JC ) = -DCONJG( S )*B( J, JC ) + C*B( J+1, JC )
               B( J, JC ) = CTEMP2
  100       CONTINUE
            IF( ILQ ) THEN
               DO 110 JR = 1, N
                  CTEMP = C*Q( JR, J ) + DCONJG( S )*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = CTEMP
  110          CONTINUE
            END IF
*
            CTEMP = B( J+1, J+1 )
            CALL ZLARTG( CTEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
            B( J+1, J ) = CZERO
*
            DO 120 JR = IFRSTM, MIN( J+2, ILAST )
               CTEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -DCONJG( S )*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = CTEMP
  120       CONTINUE
            DO 130 JR = IFRSTM, J
               CTEMP = C*B( JR, J+1 ) + S*B( JR, J )
               B( JR, J ) = -DCONJG( S )*B( JR, J+1 ) + C*B( JR, J )
               B( JR, J+1 ) = CTEMP
  130       CONTINUE
            IF( ILZ ) THEN
               DO 140 JR = 1, N
                  CTEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -DCONJG( S )*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = CTEMP
  140          CONTINUE
            END IF
  150    CONTINUE
*
  160    CONTINUE
*
  170 CONTINUE
*
*     Drop-through = non-convergence
*
  180 CONTINUE
      INFO = ILAST
      GO TO 210
*
*     Successful completion of all QZ steps
*
  190 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 200 J = 1, ILO - 1
         ABSB = ABS( B( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( J, J ) / ABSB )
            B( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, B( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, A( 1, J ), 1 )
            ELSE
               A( J, J ) = A( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            B( J, J ) = CZERO
         END IF
         ALPHA( J ) = A( J, J )
         BETA( J ) = B( J, J )
  200 CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  210 CONTINUE
      WORK( 1 ) = DCMPLX( N )
      RETURN
*
*     End of ZHGEQZ
*
      END
