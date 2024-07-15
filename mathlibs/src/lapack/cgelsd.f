      SUBROUTINE CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
     $                   WORK, LWORK, RWORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               RWORK( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CGELSD computes the minimum-norm solution to a real linear least
*  squares problem:
*      minimize 2-norm(| b - A*x |)
*  using the singular value decomposition (SVD) of A. A is an M-by-N
*  matrix which may be rank-deficient.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  The problem is solved in three steps:
*  (1) Reduce the coefficient matrix A to bidiagonal form with
*      Householder tranformations, reducing the original problem
*      into a "bidiagonal least squares problem" (BLS)
*  (2) Solve the BLS using a divide and conquer approach.
*  (3) Apply back all the Householder tranformations to solve
*      the original least squares problem.
*
*  The effective rank of A is determined by treating as zero those
*  singular values which are less than RCOND times the largest singular
*  value.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
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
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X. NRHS >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A has been destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the M-by-NRHS right hand side matrix B.
*          On exit, B is overwritten by the N-by-NRHS solution matrix X.
*          If m >= n and RANK = n, the residual sum-of-squares for
*          the solution in the i-th column is given by the sum of
*          squares of elements n+1:m in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M,N).
*
*  S       (output) REAL array, dimension (min(M,N))
*          The singular values of A in decreasing order.
*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
*
*  RCOND   (input) REAL
*          RCOND is used to determine the effective rank of A.
*          Singular values S(i) <= RCOND*S(1) are treated as zero.
*          If RCOND < 0, machine precision is used instead.
*
*  RANK    (output) INTEGER
*          The effective rank of A, i.e., the number of singular values
*          which are greater than RCOND*S(1).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK must be at least 1.
*          The exact minimum amount of workspace needed depends on M,
*          N and NRHS. As long as LWORK is at least
*              2 * N + N * NRHS
*          if M is greater than or equal to N or
*              2 * M + M * NRHS
*          if M is less than N, the code will execute correctly.
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*
*  RWORK   (workspace) REAL array, dimension at least
*             10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS +
*             (SMLSIZ+1)**2
*          if M is greater than or equal to N or
*             10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS +
*             (SMLSIZ+1)**2
*          if M is less than N, the code will execute correctly.
*          SMLSIZ is returned by ILAENV and is equal to the maximum
*          size of the subproblems at the bottom of the computation
*          tree (usually about 25), and
*             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
*
*  IWORK   (workspace) INTEGER array, dimension (LIWORK)
*          LIWORK >= 3 * MINMN * NLVL + 11 * MINMN,
*          where MINMN = MIN( M,N ).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*          > 0:  the algorithm for computing the SVD failed to converge;
*                if INFO = i, i off-diagonal elements of an intermediate
*                bidiagonal form did not converge to zero.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Ren-Cang Li, Computer Science Division, University of
*       California at Berkeley, USA
*     Osni Marques, LBNL/NERSC, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ,
     $                   LDWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM,
     $                   MNTHR, NRWORK, NWORK, SMLSIZ
      REAL               ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEBRD, CGELQF, CGEQRF, CLACPY,
     $                   CLALSD, CLASCL, CLASET, CUNMBR,
     $                   CUNMLQ, CUNMQR, SLABAD, SLASCL,
     $                   SLASET, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH, ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, 'CGELSD', ' ', M, N, NRHS, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, MAXMN ) ) THEN
         INFO = -7
      END IF
*
      SMLSIZ = ILAENV( 9, 'CGELSD', ' ', 0, 0, 0, 0 )
*
*     Compute workspace.
*     (Note: Comments in the code beginning "Workspace:" describe the
*     minimal amount of workspace needed at that point in the code,
*     as well as the preferred amount for good performance.
*     NB refers to the optimal block size for the immediately
*     following subroutine, as returned by ILAENV.)
*
      MINWRK = 1
      IF( INFO.EQ.0 ) THEN
         MAXWRK = 0
         MM = M
         IF( M.GE.N .AND. M.GE.MNTHR ) THEN
*
*           Path 1a - overdetermined, with many more rows than columns.
*
            MM = N
            MAXWRK = MAX( MAXWRK, N*ILAENV( 1, 'CGEQRF', ' ', M, N, -1,
     $               -1 ) )
            MAXWRK = MAX( MAXWRK, NRHS*ILAENV( 1, 'CUNMQR', 'LC', M,
     $               NRHS, N, -1 ) )
         END IF
         IF( M.GE.N ) THEN
*
*           Path 1 - overdetermined or exactly determined.
*
            MAXWRK = MAX( MAXWRK, 2*N+( MM+N )*
     $               ILAENV( 1, 'CGEBRD', ' ', MM, N, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, 2*N+NRHS*
     $               ILAENV( 1, 'CUNMBR', 'QLC', MM, NRHS, N, -1 ) )
            MAXWRK = MAX( MAXWRK, 2*N+( N-1 )*
     $               ILAENV( 1, 'CUNMBR', 'PLN', N, NRHS, N, -1 ) )
            MAXWRK = MAX( MAXWRK, 2*N+N*NRHS )
            MINWRK = MAX( 2*N+MM, 2*N+N*NRHS )
         END IF
         IF( N.GT.M ) THEN
            IF( N.GE.MNTHR ) THEN
*
*              Path 2a - underdetermined, with many more columns
*              than rows.
*
               MAXWRK = M + M*ILAENV( 1, 'CGELQF', ' ', M, N, -1, -1 )
               MAXWRK = MAX( MAXWRK, M*M+4*M+2*M*
     $                  ILAENV( 1, 'CGEBRD', ' ', M, M, -1, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS*
     $                  ILAENV( 1, 'CUNMBR', 'QLC', M, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )*
     $                  ILAENV( 1, 'CUNMLQ', 'LC', N, NRHS, M, -1 ) )
               IF( NRHS.GT.1 ) THEN
                  MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
               ELSE
                  MAXWRK = MAX( MAXWRK, M*M+2*M )
               END IF
               MAXWRK = MAX( MAXWRK, M*M+4*M+M*NRHS )
            ELSE
*
*              Path 2 - underdetermined.
*
               MAXWRK = 2*M + ( N+M )*ILAENV( 1, 'CGEBRD', ' ', M, N,
     $                  -1, -1 )
               MAXWRK = MAX( MAXWRK, 2*M+NRHS*
     $                  ILAENV( 1, 'CUNMBR', 'QLC', M, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, 2*M+M*
     $                  ILAENV( 1, 'CUNMBR', 'PLN', N, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, 2*M+M*NRHS )
            END IF
            MINWRK = MAX( 2*M+N, 2*M+M*NRHS )
         END IF
         MINWRK = MIN( MINWRK, MAXWRK )
         WORK( 1 ) = CMPLX( MAXWRK, 0 )
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGELSD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         GO TO 10
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
*
*     Get machine parameters.
*
      EPS = SLAMCH( 'P' )
      SFMIN = SLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Scale A if max entry outside range [SMLNUM,BIGNUM].
*
      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM.
*
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         CALL SLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 10
      END IF
*
*     Scale B if max entry outside range [SMLNUM,BIGNUM].
*
      BNRM = CLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM.
*
         CALL CLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM.
*
         CALL CLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
*     If M < N make sure B(M+1:N,:) = 0
*
      IF( M.LT.N )
     $   CALL CLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )
*
*     Overdetermined case.
*
      IF( M.GE.N ) THEN
*
*        Path 1 - overdetermined or exactly determined.
*
         MM = M
         IF( M.GE.MNTHR ) THEN
*
*           Path 1a - overdetermined, with many more rows than columns
*
            MM = N
            ITAU = 1
            NWORK = ITAU + N
*
*           Compute A=Q*R.
*           (RWorkspace: need N)
*           (CWorkspace: need N, prefer N*NB)
*
            CALL CGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),
     $                   LWORK-NWORK+1, INFO )
*
*           Multiply B by transpose(Q).
*           (RWorkspace: need N)
*           (CWorkspace: need NRHS, prefer NRHS*NB)
*
            CALL CUNMQR( 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B,
     $                   LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*           Zero out below R.
*
            IF( N.GT.1 ) THEN
               CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ),
     $                      LDA )
            END IF
         END IF
*
         ITAUQ = 1
         ITAUP = ITAUQ + N
         NWORK = ITAUP + N
         IE = 1
         NRWORK = IE + N
*
*        Bidiagonalize R in A.
*        (RWorkspace: need N)
*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
*
         CALL CGEBRD( MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,
     $                INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors of R.
*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ),
     $                B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'U', SMLSIZ, N, NRHS, S, RWORK( IE ), B, LDB,
     $                RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ),
     $                IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of R.
*
         CALL CUNMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ),
     $                B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+
     $         MAX( M, 2*M-4, NRHS, N-3*M ) ) THEN
*
*        Path 2a - underdetermined, with many more columns than rows
*        and sufficient workspace for an efficient algorithm.
*
         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ),
     $       M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         NWORK = M + 1
*
*        Compute A=L*Q.
*        (CWorkspace: need 2*M, prefer M+M*NB)
*
         CALL CGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ),
     $                LWORK-NWORK+1, INFO )
         IL = NWORK
*
*        Copy L to WORK(IL), zeroing out above its diagonal.
*
         CALL CLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL CLASET( 'U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ),
     $                LDWORK )
         ITAUQ = IL + LDWORK*M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M
*
*        Bidiagonalize L in WORK(IL).
*        (RWorkspace: need M)
*        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)
*
         CALL CGEBRD( M, M, WORK( IL ), LDWORK, S, RWORK( IE ),
     $                WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ),
     $                LWORK-NWORK+1, INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors of L.
*        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK,
     $                WORK( ITAUQ ), B, LDB, WORK( NWORK ),
     $                LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'U', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB,
     $                RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ),
     $                IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of L.
*
         CALL CUNMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK,
     $                WORK( ITAUP ), B, LDB, WORK( NWORK ),
     $                LWORK-NWORK+1, INFO )
*
*        Zero out below first M rows of B.
*
         CALL CLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )
         NWORK = ITAU + M
*
*        Multiply transpose(Q) by B.
*        (CWorkspace: need NRHS, prefer NRHS*NB)
*
         CALL CUNMLQ( 'L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B,
     $                LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      ELSE
*
*        Path 2 - remaining underdetermined cases.
*
         ITAUQ = 1
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M
*
*        Bidiagonalize A.
*        (RWorkspace: need M)
*        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
*
         CALL CGEBRD( M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1,
     $                INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors.
*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ),
     $                B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'L', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB,
     $                RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ),
     $                IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of A.
*
         CALL CUNMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ),
     $                B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      END IF
*
*     Undo scaling.
*
      IF( IASCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
*
   10 CONTINUE
      WORK( 1 ) = CMPLX( MAXWRK, 0 )
      RETURN
*
*     End of CGELSD
*
      END
