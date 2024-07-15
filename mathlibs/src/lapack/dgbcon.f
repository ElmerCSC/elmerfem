      SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, KL, KU, LDAB, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGBCON estimates the reciprocal of the condition number of a real
*  general band matrix A, in either the 1-norm or the infinity-norm,
*  using the LU factorization computed by DGBTRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as
*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  ANORM   (input) DOUBLE PRECISION
*          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*          If NORM = 'I', the infinity-norm of the original matrix A.
*
*  RCOND   (output) DOUBLE PRECISION
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, J, JP, KASE, KASE1, KD, LM
      DOUBLE PRECISION   AINVNM, SCALE, SMLNUM, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLACON, DLATBS, DRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.2*KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KD = KL + KU + 1
      LNOTI = KL.GT.0
      KASE = 0
   10 CONTINUE
      CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            IF( LNOTI ) THEN
               DO 20 J = 1, N - 1
                  LM = MIN( KL, N-J )
                  JP = IPIV( J )
                  T = WORK( JP )
                  IF( JP.NE.J ) THEN
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
                  CALL DAXPY( LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
   20          CONTINUE
            END IF
*
*           Multiply by inv(U).
*
            CALL DLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ),
     $                   INFO )
         ELSE
*
*           Multiply by inv(U').
*
            CALL DLATBS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N,
     $                   KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ),
     $                   INFO )
*
*           Multiply by inv(L').
*
            IF( LNOTI ) THEN
               DO 30 J = N - 1, 1, -1
                  LM = MIN( KL, N-J )
                  WORK( J ) = WORK( J ) - DDOT( LM, AB( KD+1, J ), 1,
     $                        WORK( J+1 ), 1 )
                  JP = IPIV( J )
                  IF( JP.NE.J ) THEN
                     T = WORK( JP )
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  END IF
   30          CONTINUE
            END IF
         END IF
*
*        Divide X by 1/SCALE if doing so will not cause overflow.
*
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 40
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   40 CONTINUE
      RETURN
*
*     End of DGBCON
*
      END
