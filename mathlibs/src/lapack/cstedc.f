      SUBROUTINE CSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,
     $                   LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), E( * ), RWORK( * )
      COMPLEX            WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEDC computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See SLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original Hermitian matrix
*                  also.  On entry, Z contains the unitary matrix used
*                  to reduce the original matrix to tridiagonal form.
*
*  N       (input) INTEGER
*          The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ,N)
*          On entry, if COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If eigenvectors are desired, then LDZ >= max(1,N).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.
*          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) REAL array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of the array RWORK.
*          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.
*          If COMPZ = 'V' and N > 1, LRWORK must be at least
*                         1 + 3*N + 2*N*lg N + 3*N**2 ,
*                         where lg( N ) = smallest integer k such
*                         that 2**k >= N.
*          If COMPZ = 'I' and N > 1, LRWORK must be at least
*                         1 + 4*N + 2*N**2 .
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the RWORK array,
*          returns this value as the first entry of the RWORK array, and
*          no error message related to LRWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.
*          If COMPZ = 'V' or N > 1,  LIWORK must be at least
*                                    6 + 6*N + 5*N*lg N.
*          If COMPZ = 'I' or N > 1,  LIWORK must be at least
*                                    3 + 5*N .
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            END, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL,
     $                   LRWMIN, LWMIN, M, SMLSIZ, START
      REAL               EPS, ORGNRM, P, TINY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANST
      EXTERNAL           ILAENV, LSAME, SLAMCH, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CLACRM, CLAED0, CSTEQR, CSWAP, SLASCL,
     $                   SLASET, SSTEDC, SSTEQR, SSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( N.LE.1 .OR. ICOMPZ.LE.0 ) THEN
         LWMIN = 1
         LIWMIN = 1
         LRWMIN = 1
      ELSE
         LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( ICOMPZ.EQ.1 ) THEN
            LWMIN = N*N
            LRWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1
            LRWMIN = 1 + 4*N + 2*N**2
            LIWMIN = 3 + 5*N
         END IF
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -8
      ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN 
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN 
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 )
     $      Z( 1, 1 ) = ONE
         RETURN 
      END IF
*
      SMLSIZ = ILAENV( 9, 'CSTEDC', ' ', 0, 0, 0, 0 )
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures SSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default.
*
*     If COMPZ = 'N', use SSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SSTERF( N, D, E, INFO )
         RETURN 
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPZ.EQ.0 ) THEN
            CALL SSTERF( N, D, E, INFO )
            RETURN 
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            CALL CSTEQR( 'I', N, D, E, Z, LDZ, RWORK, INFO )
            RETURN 
         ELSE
            CALL CSTEQR( 'V', N, D, E, Z, LDZ, RWORK, INFO )
            RETURN 
         END IF
      END IF
*
*     If COMPZ = 'I', we simply call SSTEDC instead.
*
      IF( ICOMPZ.EQ.2 ) THEN
         CALL SLASET( 'Full', N, N, ZERO, ONE, RWORK, N )
         LL = N*N + 1
         CALL SSTEDC( 'I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1,
     $                IWORK, LIWORK, INFO )
         DO 20 J = 1, N
            DO 10 I = 1, N
               Z( I, J ) = RWORK( ( J-1 )*N+I )
   10       CONTINUE
   20    CONTINUE
         RETURN 
      END IF
*
*     From now on, only option left to be handled is COMPZ = 'V',
*     i.e. ICOMPZ = 1.
*
*     Scale.
*
      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO )
     $   RETURN
*
      EPS = SLAMCH( 'Epsilon' )
*
      START = 1
*
*     while ( START <= N )
*
   30 CONTINUE
      IF( START.LE.N ) THEN
*
*     Let END be the position of the next subdiagonal entry such that
*     E( END ) <= TINY or END = N if no such subdiagonal exists.  The
*     matrix identified by the elements between START and END
*     constitutes an independent sub-problem.
*
         END = START
   40    CONTINUE
         IF( END.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( END ) ) )*SQRT( ABS( D( END+1 ) ) )
            IF( ABS( E( END ) ).GT.TINY ) THEN
               END = END + 1
               GO TO 40
            END IF
         END IF
*
*        (Sub) Problem determined.  Compute its size and solve it.
*
         M = END - START + 1
         IF( M.GT.SMLSIZ ) THEN
            INFO = SMLSIZ
*
*           Scale.
*
            ORGNRM = SLANST( 'M', M, D( START ), E( START ) )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M,
     $                   INFO )
            CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ),
     $                   M-1, INFO )
*
            CALL CLAED0( N, M, D( START ), E( START ), Z( 1, START ),
     $                   LDZ, WORK, N, RWORK, IWORK, INFO )
            IF( INFO.GT.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +
     $                MOD( INFO, ( M+1 ) ) + START - 1
               RETURN 
            END IF
*
*           Scale back.
*
            CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M,
     $                   INFO )
*
         ELSE
            CALL SSTEQR( 'I', M, D( START ), E( START ), RWORK, M,
     $                   RWORK( M*M+1 ), INFO )
            CALL CLACRM( N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N,
     $                   RWORK( M*M+1 ) )
            CALL CLACPY( 'A', N, M, WORK, N, Z( 1, START ), LDZ )
            IF( INFO.GT.0 ) THEN
               INFO = START*( N+1 ) + END
               RETURN 
            END IF
         END IF
*
         START = END + 1
         GO TO 30
      END IF
*
*     endwhile
*
*     If the problem split any number of times, then the eigenvalues
*     will not be properly ordered.  Here we permute the eigenvalues
*     (and the associated eigenvectors) into ascending order.
*
      IF( M.NE.N ) THEN
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 60 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 50 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
   50       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
   60    CONTINUE
      END IF
*
      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CSTEDC
*
      END
