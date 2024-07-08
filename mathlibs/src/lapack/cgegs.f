      SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  This routine is deprecated and has been replaced by routine CGGES.
*
*  CGEGS computes for a pair of N-by-N complex nonsymmetric matrices A,
*  B:  the generalized eigenvalues (alpha, beta), the complex Schur
*  form (A, B), and optionally left and/or right Schur vectors
*  (VSL and VSR).
*
*  (If only the generalized eigenvalues are needed, use the driver CGEGV
*  instead.)
*
*  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
*  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
*  is singular.  It is usually represented as the pair (alpha,beta),
*  as there is a reasonable interpretation for beta=0, and even for
*  both being zero.  A good beginning reference is the book, "Matrix
*  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
*
*  The (generalized) Schur form of a pair of matrices is the result of
*  multiplying both matrices on the left by one unitary matrix and
*  both on the right by another unitary matrix, these two unitary
*  matrices being chosen so as to bring the pair of matrices into
*  upper triangular form with the diagonal elements of B being
*  non-negative real numbers (this is also called complex Schur form.)
*
*  The left and right Schur vectors are the columns of VSL and VSR,
*  respectively, where VSL and VSR are the unitary matrices
*  which reduce A and B to Schur form:
*
*  Schur form of (A,B) = ( (VSL)**H A (VSR), (VSL)**H B (VSR) )
*
*  Arguments
*  =========
*
*  JOBVSL   (input) CHARACTER*1
*          = 'N':  do not compute the left Schur vectors;
*          = 'V':  compute the left Schur vectors.
*
*  JOBVSR   (input) CHARACTER*1
*          = 'N':  do not compute the right Schur vectors;
*          = 'V':  compute the right Schur vectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the first of the pair of matrices whose generalized
*          eigenvalues and (optionally) Schur vectors are to be
*          computed.
*          On exit, the generalized Schur form of A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX array, dimension (LDB, N)
*          On entry, the second of the pair of matrices whose
*          generalized eigenvalues and (optionally) Schur vectors are
*          to be computed.
*          On exit, the generalized Schur form of B.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX array, dimension (N)
*  BETA    (output) COMPLEX array, dimension (N)
*          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
*          j=1,...,N  are the diagonals of the complex Schur form (A,B)
*          output by CGEGS.  The  BETA(j) will be non-negative real.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VSL     (output) COMPLEX array, dimension (LDVSL,N)
*          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*          (See "Purpose", above.)
*          Not referenced if JOBVSL = 'N'.
*
*  LDVSL   (input) INTEGER
*          The leading dimension of the matrix VSL. LDVSL >= 1, and
*          if JOBVSL = 'V', LDVSL >= N.
*
*  VSR     (output) COMPLEX array, dimension (LDVSR,N)
*          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*          (See "Purpose", above.)
*          Not referenced if JOBVSR = 'N'.
*
*  LDVSR   (input) INTEGER
*          The leading dimension of the matrix VSR. LDVSR >= 1, and
*          if JOBVSR = 'V', LDVSR >= N.
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*          To compute the optimal value of LWORK, call ILAENV to get
*          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:
*          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;
*          the optimal LWORK is N*(NB+1).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) REAL array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          =1,...,N:
*                The QZ iteration failed.  (A,B) are not in Schur
*                form, but ALPHA(j) and BETA(j) should be correct for
*                j=INFO+1,...,N.
*          > N:  errors that usually indicate LAPACK problems:
*                =N+1: error return from CGGBAL
*                =N+2: error return from CGEQRF
*                =N+3: error return from CUNMQR
*                =N+4: error return from CUNGQR
*                =N+5: error return from CGGHRD
*                =N+6: error return from CHGEQZ (other than failed
*                                               iteration)
*                =N+7: error return from CGGBAK (computing VSL)
*                =N+8: error return from CGGBAK (computing VSR)
*                =N+9: error return from CLASCL (various places)
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),
     $                   CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT,
     $                   ILO, IRIGHT, IROWS, IRWORK, ITAU, IWORK,
     $                   LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SAFMIN, SMLNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY,
     $                   CLASCL, CLASET, CUNGQR, CUNMQR, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           ILAENV, LSAME, CLANGE, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVSL, 'N' ) ) THEN
         IJOBVL = 1
         ILVSL = .FALSE.
      ELSE IF( LSAME( JOBVSL, 'V' ) ) THEN
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      END IF
*
      IF( LSAME( JOBVSR, 'N' ) ) THEN
         IJOBVR = 1
         ILVSR = .FALSE.
      ELSE IF( LSAME( JOBVSR, 'V' ) ) THEN
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      END IF
*
*     Test the input arguments
*
      LWKMIN = MAX( 2*N, 1 )
      LWKOPT = LWKMIN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB1 = ILAENV( 1, 'CGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'CUNMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'CUNGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = N*(NB+1)
         WORK( 1 ) = LOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEGS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = SLAMCH( 'E' )*SLAMCH( 'B' )
      SAFMIN = SLAMCH( 'S' )
      SMLNUM = N*SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
*
      IF( ILASCL ) THEN
         CALL CLASCL( 'G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
      END IF
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
*
      IF( ILBSCL ) THEN
         CALL CLASCL( 'G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWORK = IRIGHT + N
      IWORK = 1
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWORK ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 1
         GO TO 10
      END IF
*
*     Reduce B to triangular form, and initialize VSL and/or VSR
*
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 2
         GO TO 10
      END IF
*
      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ),
     $             LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 3
         GO TO 10
      END IF
*
      IF( ILVSL ) THEN
         CALL CLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                VSL( ILO+1, ILO ), LDVSL )
         CALL CUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL,
     $                WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK,
     $                IINFO )
         IF( IINFO.GE.0 )
     $      LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 4
            GO TO 10
         END IF
      END IF
*
      IF( ILVSR )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )
*
*     Reduce to generalized Hessenberg form
*
      CALL CGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL,
     $             LDVSL, VSR, LDVSR, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 5
         GO TO 10
      END IF
*
*     Perform QZ algorithm, computing Schur vectors if desired
*
      IWORK = ITAU
      CALL CHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ),
     $             LWORK+1-IWORK, RWORK( IRWORK ), IINFO )
      IF( IINFO.GE.0 )
     $   LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         IF( IINFO.GT.0 .AND. IINFO.LE.N ) THEN
            INFO = IINFO
         ELSE IF( IINFO.GT.N .AND. IINFO.LE.2*N ) THEN
            INFO = IINFO - N
         ELSE
            INFO = N + 6
         END IF
         GO TO 10
      END IF
*
*     Apply permutation to VSL and VSR
*
      IF( ILVSL ) THEN
         CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSL, LDVSL, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 7
            GO TO 10
         END IF
      END IF
      IF( ILVSR ) THEN
         CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSR, LDVSR, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 8
            GO TO 10
         END IF
      END IF
*
*     Undo scaling
*
      IF( ILASCL ) THEN
         CALL CLASCL( 'U', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
         CALL CLASCL( 'G', -1, -1, ANRMTO, ANRM, N, 1, ALPHA, N, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
      END IF
*
      IF( ILBSCL ) THEN
         CALL CLASCL( 'U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
         CALL CLASCL( 'G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 9
            RETURN
         END IF
      END IF
*
   10 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of CGEGS
*
      END
