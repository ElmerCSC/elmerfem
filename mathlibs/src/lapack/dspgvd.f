      SUBROUTINE DSPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
     $                   LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSPGVD computes all the eigenvalues, and optionally, the eigenvectors
*  of a real generalized symmetric-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
*  B are assumed to be symmetric, stored in packed format, and B is also
*  positive definite.
*  If eigenvectors are desired, it uses a divide and conquer algorithm.
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
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the symmetric matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*
*          On exit, the contents of AP are destroyed.
*
*  BP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the symmetric matrix
*          B, packed columnwise in a linear array.  The j-th column of B
*          is stored in the array BP as follows:
*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
*
*          On exit, the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T, in the same storage
*          format as B.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
*          eigenvectors.  The eigenvectors are normalized as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK >= 1.
*          If JOBZ = 'N' and N > 1, LWORK >= 2*N.
*          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.
*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  DPPTRF or DSPEVD returned an error code:
*             <= N:  if INFO = i, DSPEVD failed to converge;
*                    i off-diagonal elements of an intermediate
*                    tridiagonal form did not converge to zero;
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
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
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            J, LGN, LIWMIN, LWMIN, NEIG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DPPTRF, DSPEVD, DSPGST, DTPMV, DTPSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LGN = 0
         LIWMIN = 1
         LWMIN = 1
      ELSE
         LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( 2**LGN.LT.N )
     $      LGN = LGN + 1
         IF( WANTZ ) THEN
            LIWMIN = 3 + 5*N
            LWMIN = 1 + 5*N + 2*N*LGN + 2*N**2
         ELSE
            LIWMIN = 1
            LWMIN = 2*N
         END IF
      END IF
*
      IF( ITYPE.LT.0 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -11
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPGVD', -INFO )
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
*     Form a Cholesky factorization of BP.
*
      CALL DPPTRF( UPLO, N, BP, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem and solve.
*
      CALL DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
      CALL DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
      LWMIN = MAX( DBLE( LWMIN ), DBLE( WORK( 1 ) ) )
      LIWMIN = MAX( DBLE( LIWMIN ), DBLE( IWORK( 1 ) ) )
*
      IF( WANTZ ) THEN
*
*        Backtransform eigenvectors to the original problem.
*
         NEIG = N
         IF( INFO.GT.0 )
     $      NEIG = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
*
*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
*           backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
*
            DO 10 J = 1, NEIG
               CALL DTPSV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ),
     $                     1 )
   10       CONTINUE
*
         ELSE IF( ITYPE.EQ.3 ) THEN
*
*           For B*A*x=(lambda)*x;
*           backtransform eigenvectors: x = L*y or U'*y
*
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
*
            DO 20 J = 1, NEIG
               CALL DTPMV( UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ),
     $                     1 )
   20       CONTINUE
         END IF
      END IF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of DSPGVD
*
      END
