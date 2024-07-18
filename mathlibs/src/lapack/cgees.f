      SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
     $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
*     ..
*
*  Purpose
*  =======
*
*  CGEES computes for an N-by-N complex nonsymmetric matrix A, the
*  eigenvalues, the Schur form T, and, optionally, the matrix of Schur
*  vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
*
*  Optionally, it also orders the eigenvalues on the diagonal of the
*  Schur form so that selected eigenvalues are at the top left.
*  The leading columns of Z then form an orthonormal basis for the
*  invariant subspace corresponding to the selected eigenvalues.

*  A complex matrix is in Schur form if it is upper triangular.
*
*  Arguments
*  =========
*
*  JOBVS   (input) CHARACTER*1
*          = 'N': Schur vectors are not computed;
*          = 'V': Schur vectors are computed.
*
*  SORT    (input) CHARACTER*1
*          Specifies whether or not to order the eigenvalues on the
*          diagonal of the Schur form.
*          = 'N': Eigenvalues are not ordered:
*          = 'S': Eigenvalues are ordered (see SELECT).
*
*  SELECT  (input) LOGICAL FUNCTION of one COMPLEX argument
*          SELECT must be declared EXTERNAL in the calling subroutine.
*          If SORT = 'S', SELECT is used to select eigenvalues to order
*          to the top left of the Schur form.
*          IF SORT = 'N', SELECT is not referenced.
*          The eigenvalue W(j) is selected if SELECT(W(j)) is true.
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten by its Schur form T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  SDIM    (output) INTEGER
*          If SORT = 'N', SDIM = 0.
*          If SORT = 'S', SDIM = number of eigenvalues for which
*                         SELECT is true.
*
*  W       (output) COMPLEX array, dimension (N)
*          W contains the computed eigenvalues, in the same order that
*          they appear on the diagonal of the output Schur form T.
*
*  VS      (output) COMPLEX array, dimension (LDVS,N)
*          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
*          vectors.
*          If JOBVS = 'N', VS is not referenced.
*
*  LDVS    (input) INTEGER
*          The leading dimension of the array VS.  LDVS >= 1; if
*          JOBVS = 'V', LDVS >= N.
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  BWORK   (workspace) LOGICAL array, dimension (N)
*          Not referenced if SORT = 'N'.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*          > 0: if INFO = i, and i is
*               <= N:  the QR algorithm failed to compute all the
*                      eigenvalues; elements 1:ILO-1 and i+1:N of W
*                      contain those eigenvalues which have converged;
*                      if JOBVS = 'V', VS contains the matrix which
*                      reduces A to its partially converged Schur form.
*               = N+1: the eigenvalues could not be reordered because
*                      some eigenvalues were too close to separate (the
*                      problem is very ill-conditioned);
*               = N+2: after reordering, roundoff changed values of
*                      some complex eigenvalues so that leading
*                      eigenvalues in the Schur form no longer satisfy
*                      SELECT = .TRUE..  This could also be caused by
*                      underflow due to scaling.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTST, WANTVS
      INTEGER            HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO,
     $                   ITAU, IWRK, K, MAXB, MAXWRK, MINWRK
      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
*     ..
*     .. Local Arrays ..
      REAL               DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY,
     $                   CLASCL, CTRSEN, CUNGHR, SLABAD, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
         INFO = -10
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by CHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      MINWRK = 1
      IF( INFO.EQ.0 .AND. ( LWORK.GE.1 .OR. LQUERY ) ) THEN
         MAXWRK = N + N*ILAENV( 1, 'CGEHRD', ' ', N, 1, N, 0 )
         MINWRK = MAX( 1, 2*N )
         IF( .NOT.WANTVS ) THEN
            MAXB = MAX( ILAENV( 8, 'CHSEQR', 'SN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'CHSEQR', 'SN', N, 1,
     $          N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, HSWORK, 1 )
         ELSE
            MAXWRK = MAX( MAXWRK, N+( N-1 )*
     $               ILAENV( 1, 'CUNGHR', ' ', N, 1, N, -1 ) )
            MAXB = MAX( ILAENV( 8, 'CHSEQR', 'EN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'CHSEQR', 'EN', N, 1,
     $          N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, HSWORK, 1 )
         END IF
         WORK( 1 ) = MAXWRK
      END IF
      IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL CLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
*     Permute the matrix to make it more nearly triangular
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      IBAL = 1
      CALL CGEBAL( 'P', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      ITAU = 1
      IWRK = N + ITAU
      CALL CGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $             LWORK-IWRK+1, IERR )
*
      IF( WANTVS ) THEN
*
*        Copy Householder vectors to VS
*
         CALL CLACPY( 'L', N, N, A, LDA, VS, LDVS )
*
*        Generate unitary matrix in VS
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL CUNGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
      END IF
*
      SDIM = 0
*
*     Perform QR iteration, accumulating Schur vectors in VS if desired
*     (CWorkspace: need 1, prefer HSWORK (see comments) )
*     (RWorkspace: none)
*
      IWRK = ITAU
      CALL CHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, W, VS, LDVS,
     $             WORK( IWRK ), LWORK-IWRK+1, IEVAL )
      IF( IEVAL.GT.0 )
     $   INFO = IEVAL
*
*     Sort eigenvalues if desired
*
      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA )
     $      CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, W, N, IERR )
         DO 10 I = 1, N
            BWORK( I ) = SELECT( W( I ) )
   10    CONTINUE
*
*        Reorder eigenvalues and transform Schur vectors
*        (CWorkspace: none)
*        (RWorkspace: none)
*
         CALL CTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, W, SDIM,
     $                S, SEP, WORK( IWRK ), LWORK-IWRK+1, ICOND )
      END IF
*
      IF( WANTVS ) THEN
*
*        Undo balancing
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL CGEBAK( 'P', 'R', N, ILO, IHI, RWORK( IBAL ), N, VS, LDVS,
     $                IERR )
      END IF
*
      IF( SCALEA ) THEN
*
*        Undo scaling for the Schur form of A
*
         CALL CLASCL( 'U', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL CCOPY( N, A, LDA+1, W, 1 )
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of CGEES
*
      END
