      SUBROUTINE SSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,
     $                   LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,
     $                   LDZ, WORK, IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,
     $                   N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL               AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),
     $                   W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SSBGVX computes selected eigenvalues, and optionally, eigenvectors
*  of a real generalized symmetric-definite banded eigenproblem, of
*  the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
*  and banded, and B is also positive definite.  Eigenvalues and
*  eigenvectors can be selected by specifying either all eigenvalues,
*  a range of values or a range of indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  KA      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
*
*  KB      (input) INTEGER
*          The number of superdiagonals of the matrix B if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first ka+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
*
*          On exit, the contents of AB are destroyed.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KA+1.
*
*  BB      (input/output) REAL array, dimension (LDBB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix B, stored in the first kb+1 rows of the array.  The
*          j-th column of B is stored in the j-th column of the array BB
*          as follows:
*          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
*
*          On exit, the factor S from the split Cholesky factorization
*          B = S**T*S, as returned by SPBSTF.
*
*  LDBB    (input) INTEGER
*          The leading dimension of the array BB.  LDBB >= KB+1.
*
*  Q       (output) REAL array, dimension (LDQ, N)
*          If JOBZ = 'V', the n-by-n matrix used in the reduction of
*          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,
*          and consequently C to tridiagonal form.
*          If JOBZ = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  If JOBZ = 'N',
*          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*SLAMCH('S').
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) REAL array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
*          eigenvectors, with the i-th column of Z holding the
*          eigenvector associated with W(i).  The eigenvectors are
*          normalized so Z**T*B*Z = I.
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) REAL array, dimension (7N)
*
*  IWORK   (workspace/output) INTEGER array, dimension (5N)
*
*  IFAIL   (input) INTEGER array, dimension (M)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvalues that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0 : successful exit
*          < 0 : if INFO = -i, the i-th argument had an illegal value
*          <= N: if INFO = i, then i eigenvectors failed to converge.
*                  Their indices are stored in IFAIL.
*          > N : SPBSTF returned an error code; i.e.,
*                if INFO = N + i, for 1 <= i <= N, then the leading
*                minor of order i of B is not positive definite.
*                The factorization of B could not be completed and
*                no eigenvalues or eigenvectors were computed.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, UPPER, VALEIG, WANTZ
      CHARACTER          ORDER, VECT
      INTEGER            I, IINFO, INDD, INDE, INDEE, INDIBL, INDISP,
     $                   INDIWO, INDWRK, ITMP1, J, JJ, NSPLIT
      REAL               TMP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV, SLACPY, SPBSTF, SSBGST, SSBTRD,
     $                   SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KA.LT.0 ) THEN
         INFO = -5
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -8
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 ) THEN
         INFO = -12
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -14
      ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
         INFO = -15
      ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
         INFO = -16
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSBGVX', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Form a split Cholesky factorization of B.
*
      CALL SPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem.
*
      CALL SSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ,
     $             WORK, IINFO )
*
*     Reduce symmetric band matrix to tridiagonal form.
*
      INDD = 1
      INDE = INDD + N
      INDWRK = INDE + N
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL SSBTRD( VECT, UPLO, N, KA, AB, LDAB, WORK( INDD ),
     $             WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO )
*
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then call SSTERF or SSTEQR.  If this fails for some
*     eigenvalue, then try SSTEBZ.
*
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.
     $    ( ABSTOL.LE.ZERO ) ) THEN
         CALL SCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
         IF( .NOT.WANTZ ) THEN
            CALL SSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL SLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
            CALL SSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,
     $                   WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call SSTEBZ and, if eigenvectors are desired,
*     call SSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL SSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL,
     $             WORK( INDD ), WORK( INDE ), M, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL SSTEIN( N, WORK( INDD ), WORK( INDE ), M, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
*
*        Apply transformation matrix used in reduction to tridiagonal
*        form to eigenvectors returned by SSTEIN.
*
         DO 20 J = 1, M
            CALL SCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL SGEMV( 'N', N, N, ONE, Q, LDQ, WORK, 1, ZERO,
     $                  Z( 1, J ), 1 )
   20    CONTINUE
      END IF
*
   30 CONTINUE
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   40       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of SSBGVX
*
      END
