      SUBROUTINE CTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      REAL               RWORK( * ), S( * ), SEP( * )
      COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  CTRSNA estimates reciprocal condition numbers for specified
*  eigenvalues and/or right eigenvectors of a complex upper triangular
*  matrix T (or of any matrix Q*T*Q**H with Q unitary).
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for
*          eigenvalues (S) or eigenvectors (SEP):
*          = 'E': for eigenvalues only (S);
*          = 'V': for eigenvectors only (SEP);
*          = 'B': for both eigenvalues and eigenvectors (S and SEP).
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute condition numbers for all eigenpairs;
*          = 'S': compute condition numbers for selected eigenpairs
*                 specified by the array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
*          condition numbers are required. To select condition numbers
*          for the j-th eigenpair, SELECT(j) must be set to .TRUE..
*          If HOWMNY = 'A', SELECT is not referenced.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) COMPLEX array, dimension (LDT,N)
*          The upper triangular matrix T.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input) COMPLEX array, dimension (LDVL,M)
*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
*          (or of any Q*T*Q**H with Q unitary), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VL, as returned by
*          CHSEIN or CTREVC.
*          If JOB = 'V', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
*
*  VR      (input) COMPLEX array, dimension (LDVR,M)
*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
*          (or of any Q*T*Q**H with Q unitary), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VR, as returned by
*          CHSEIN or CTREVC.
*          If JOB = 'V', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
*
*  S       (output) REAL array, dimension (MM)
*          If JOB = 'E' or 'B', the reciprocal condition numbers of the
*          selected eigenvalues, stored in consecutive elements of the
*          array. Thus S(j), SEP(j), and the j-th columns of VL and VR
*          all correspond to the same eigenpair (but not in general the
*          j-th eigenpair, unless all eigenpairs are selected).
*          If JOB = 'V', S is not referenced.
*
*  SEP     (output) REAL array, dimension (MM)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the selected eigenvectors, stored in consecutive
*          elements of the array.
*          If JOB = 'E', SEP is not referenced.
*
*  MM      (input) INTEGER
*          The number of elements in the arrays S (if JOB = 'E' or 'B')
*           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
*
*  M       (output) INTEGER
*          The number of elements of the arrays S and/or SEP actually
*          used to store the estimated condition numbers.
*          If HOWMNY = 'A', M is set to N.
*
*  WORK    (workspace) COMPLEX array, dimension (LDWORK,N+1)
*          If JOB = 'E', WORK is not referenced.
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
*
*  RWORK   (workspace) REAL array, dimension (N)
*          If JOB = 'E', RWORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The reciprocal of the condition number of an eigenvalue lambda is
*  defined as
*
*          S(lambda) = |v'*u| / (norm(u)*norm(v))
*
*  where u and v are the right and left eigenvectors of T corresponding
*  to lambda; v' denotes the conjugate transpose of v, and norm(u)
*  denotes the Euclidean norm. These reciprocal condition numbers always
*  lie between zero (very badly conditioned) and one (very well
*  conditioned). If n = 1, S(lambda) is defined to be 1.
*
*  An approximate error bound for a computed eigenvalue W(i) is given by
*
*                      EPS * norm(T) / S(i)
*
*  where EPS is the machine precision.
*
*  The reciprocal of the condition number of the right eigenvector u
*  corresponding to lambda is defined as follows. Suppose
*
*              T = ( lambda  c  )
*                  (   0    T22 )
*
*  Then the reciprocal condition number is
*
*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
*
*  where sigma-min denotes the smallest singular value. We approximate
*  the smallest singular value by the reciprocal of an estimate of the
*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
*  defined to be abs(T(1,1)).
*
*  An approximate error bound for a computed right eigenvector VR(i)
*  is given by
*
*                      EPS * norm(T) / SEP(i)
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SOMCON, WANTBH, WANTS, WANTSP
      CHARACTER          NORMIN
      INTEGER            I, IERR, IX, J, K, KASE, KS
      REAL               BIGNUM, EPS, EST, LNRM, RNRM, SCALE, SMLNUM,
     $                   XNORM
      COMPLEX            CDUM, PROD
*     ..
*     .. Local Arrays ..
      COMPLEX            DUMMY( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICAMAX
      REAL               SCNRM2, SLAMCH
      COMPLEX            CDOTC
      EXTERNAL           LSAME, ICAMAX, SCNRM2, SLAMCH, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACON, CLACPY, CLATRS, CSRSCL, CTREXC, SLABAD,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
*     Set M to the number of eigenpairs for which condition numbers are
*     to be computed.
*
      IF( SOMCON ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
*
      INFO = 0
      IF( .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -13
      ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRSNA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) )
     $         RETURN
         END IF
         IF( WANTS )
     $      S( 1 ) = ONE
         IF( WANTSP )
     $      SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
      KS = 1
      DO 50 K = 1, N
*
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) )
     $         GO TO 50
         END IF
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            PROD = CDOTC( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
            RNRM = SCNRM2( N, VR( 1, KS ), 1 )
            LNRM = SCNRM2( N, VL( 1, KS ), 1 )
            S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
*
         END IF
*
         IF( WANTSP ) THEN
*
*           Estimate the reciprocal condition number of the k-th
*           eigenvector.
*
*           Copy the matrix T to the array WORK and swap the k-th
*           diagonal element to the (1,1) position.
*
            CALL CLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            CALL CTREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, K, 1, IERR )
*
*           Form  C = T22 - lambda*I in WORK(2:N,2:N).
*
            DO 20 I = 2, N
               WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20       CONTINUE
*
*           Estimate a lower bound for the 1-norm of inv(C'). The 1st
*           and (N+1)th columns of WORK are used to store work vectors.
*
            SEP( KS ) = ZERO
            EST = ZERO
            KASE = 0
            NORMIN = 'N'
   30       CONTINUE
            CALL CLACON( N-1, WORK( 1, N+1 ), WORK, EST, KASE )
*
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Solve C'*x = scale*b
*
                  CALL CLATRS( 'Upper', 'Conjugate transpose',
     $                         'Nonunit', NORMIN, N-1, WORK( 2, 2 ),
     $                         LDWORK, WORK, SCALE, RWORK, IERR )
               ELSE
*
*                 Solve C*x = scale*b
*
                  CALL CLATRS( 'Upper', 'No transpose', 'Nonunit',
     $                         NORMIN, N-1, WORK( 2, 2 ), LDWORK, WORK,
     $                         SCALE, RWORK, IERR )
               END IF
               NORMIN = 'Y'
               IF( SCALE.NE.ONE ) THEN
*
*                 Multiply by 1/SCALE if doing so will not cause
*                 overflow.
*
                  IX = ICAMAX( N-1, WORK, 1 )
                  XNORM = CABS1( WORK( IX, 1 ) )
                  IF( SCALE.LT.XNORM*SMLNUM .OR. SCALE.EQ.ZERO )
     $               GO TO 40
                  CALL CSRSCL( N, SCALE, WORK, 1 )
               END IF
               GO TO 30
            END IF
*
            SEP( KS ) = ONE / MAX( EST, SMLNUM )
         END IF
*
   40    CONTINUE
         KS = KS + 1
   50 CONTINUE
      RETURN
*
*     End of CTRSNA
*
      END
