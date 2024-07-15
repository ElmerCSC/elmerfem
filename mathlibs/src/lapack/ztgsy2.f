      SUBROUTINE ZTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,
     $                   INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N
      DOUBLE PRECISION   RDSCAL, RDSUM, SCALE
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), E( LDE, * ), F( LDF, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTGSY2 solves the generalized Sylvester equation
*
*              A * R - L * B = scale *   C               (1)
*              D * R - L * E = scale * F
*
*  using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
*  (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
*  N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
*  (i.e., (A,D) and (B,E) in generalized Schur form).
*
*  The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
*  scaling factor chosen to avoid overflow.
*
*  In matrix notation solving equation (1) corresponds to solve
*  Zx = scale * b, where Z is defined as
*
*         Z = [ kron(In, A)  -kron(B', Im) ]             (2)
*             [ kron(In, D)  -kron(E', Im) ],
*
*  Ik is the identity matrix of size k and X' is the transpose of X.
*  kron(X, Y) is the Kronecker product between the matrices X and Y.
*
*  If TRANS = 'C', y in the conjugate transposed system Z'y = scale*b
*  is solved for, which is equivalent to solve for R and L in
*
*              A' * R  + D' * L   = scale *  C           (3)
*              R  * B' + L  * E'  = scale * -F
*
*  This case is used to compute an estimate of Dif[(A, D), (B, E)] =
*  = sigma_min(Z) using reverse communicaton with ZLACON.
*
*  ZTGSY2 also (IJOB >= 1) contributes to the computation in ZTGSYL
*  of an upper bound on the separation between to matrix pairs. Then
*  the input (A, D), (B, E) are sub-pencils of two matrix pairs in
*  ZTGSYL.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          = 'N', solve the generalized Sylvester equation (1).
*          = 'T': solve the 'transposed' system (3).
*
*  IJOB    (input) INTEGER
*          Specifies what kind of functionality to be performed.
*          =0: solve (1) only.
*          =1: A contribution from this subsystem to a Frobenius
*              norm-based estimate of the separation between two matrix
*              pairs is computed. (look ahead strategy is used).
*          =2: A contribution from this subsystem to a Frobenius
*              norm-based estimate of the separation between two matrix
*              pairs is computed. (DGECON on sub-systems is used.)
*          Not referenced if TRANS = 'T'.
*
*  M       (input) INTEGER
*          On entry, M specifies the order of A and D, and the row
*          dimension of C, F, R and L.
*
*  N       (input) INTEGER
*          On entry, N specifies the order of B and E, and the column
*          dimension of C, F, R and L.
*
*  A       (input) COMPLEX*16 array, dimension (LDA, M)
*          On entry, A contains an upper triangular matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the matrix A. LDA >= max(1, M).
*
*  B       (input) COMPLEX*16 array, dimension (LDB, N)
*          On entry, B contains an upper triangular matrix.
*
*  LDB     (input) INTEGER
*          The leading dimension of the matrix B. LDB >= max(1, N).
*
*  C       (input/ output) COMPLEX*16 array, dimension (LDC, N)
*          On entry, C contains the right-hand-side of the first matrix
*          equation in (1).
*          On exit, if IJOB = 0, C has been overwritten by the solution
*          R.
*
*  LDC     (input) INTEGER
*          The leading dimension of the matrix C. LDC >= max(1, M).
*
*  D       (input) COMPLEX*16 array, dimension (LDD, M)
*          On entry, D contains an upper triangular matrix.
*
*  LDD     (input) INTEGER
*          The leading dimension of the matrix D. LDD >= max(1, M).
*
*  E       (input) COMPLEX*16 array, dimension (LDE, N)
*          On entry, E contains an upper triangular matrix.
*
*  LDE     (input) INTEGER
*          The leading dimension of the matrix E. LDE >= max(1, N).
*
*  F       (input/ output) COMPLEX*16 array, dimension (LDF, N)
*          On entry, F contains the right-hand-side of the second matrix
*          equation in (1).
*          On exit, if IJOB = 0, F has been overwritten by the solution
*          L.
*
*  LDF     (input) INTEGER
*          The leading dimension of the matrix F. LDF >= max(1, M).
*
*  SCALE   (output) DOUBLE PRECISION
*          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions
*          R and L (C and F on entry) will hold the solutions to a
*          slightly perturbed system but the input matrices A, B, D and
*          E have not been changed. If SCALE = 0, R and L will hold the
*          solutions to the homogeneous system with C = F = 0.
*          Normally, SCALE = 1.
*
*  RDSUM   (input/output) DOUBLE PRECISION
*          On entry, the sum of squares of computed contributions to
*          the Dif-estimate under computation by ZTGSYL, where the
*          scaling factor RDSCAL (see below) has been factored out.
*          On exit, the corresponding sum of squares updated with the
*          contributions from the current sub-system.
*          If TRANS = 'T' RDSUM is not touched.
*          NOTE: RDSUM only makes sense when ZTGSY2 is called by
*          ZTGSYL.
*
*  RDSCAL  (input/output) DOUBLE PRECISION
*          On entry, scaling factor used to prevent overflow in RDSUM.
*          On exit, RDSCAL is updated w.r.t. the current contributions
*          in RDSUM.
*          If TRANS = 'T', RDSCAL is not touched.
*          NOTE: RDSCAL only makes sense when ZTGSY2 is called by
*          ZTGSYL.
*
*  INFO    (output) INTEGER
*          On exit, if INFO is set to
*            =0: Successful exit
*            <0: If INFO = -i, input argument number i is illegal.
*            >0: The matrix pairs (A, D) and (B, E) have common or very
*                close eigenvalues.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      INTEGER            LDZ
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, LDZ = 2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            I, IERR, J, K
      DOUBLE PRECISION   SCALOC
      COMPLEX*16         ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            IPIV( LDZ ), JPIV( LDZ )
      COMPLEX*16         RHS( LDZ ), Z( LDZ, LDZ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZGESC2, ZGETC2, ZLATDF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input parameters
*
      INFO = 0
      IERR = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.2 ) ) THEN
         INFO = -2
      ELSE IF( M.LE.0 ) THEN
         INFO = -3
      ELSE IF( N.LE.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
         INFO = -12
      ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
         INFO = -14
      ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
         INFO = -16
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSY2', -INFO )
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
*
*        Solve (I, J) - system
*           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
*           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
*        for I = M, M - 1, ..., 1; J = 1, 2, ..., N
*
         SCALE = ONE
         SCALOC = ONE
         DO 30 J = 1, N
            DO 20 I = M, 1, -1
*
*              Build 2 by 2 system
*
               Z( 1, 1 ) = A( I, I )
               Z( 2, 1 ) = D( I, I )
               Z( 1, 2 ) = -B( J, J )
               Z( 2, 2 ) = -E( J, J )
*
*              Set up right hand side(s)
*
               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )
*
*              Solve Z * x = RHS
*
               CALL ZGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 )
     $            INFO = IERR
               IF( IJOB.EQ.0 ) THEN
                  CALL ZGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 K = 1, N
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              C( 1, K ), 1 )
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              F( 1, K ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
               ELSE
                  CALL ZLATDF( IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL,
     $                         IPIV, JPIV )
               END IF
*
*              Unpack solution vector(s)
*
               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )
*
*              Substitute R(I, J) and L(I, J) into remaining equation.
*
               IF( I.GT.1 ) THEN
                  ALPHA = -RHS( 1 )
                  CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, C( 1, J ), 1 )
                  CALL ZAXPY( I-1, ALPHA, D( 1, I ), 1, F( 1, J ), 1 )
               END IF
               IF( J.LT.N ) THEN
                  CALL ZAXPY( N-J, RHS( 2 ), B( J, J+1 ), LDB,
     $                        C( I, J+1 ), LDC )
                  CALL ZAXPY( N-J, RHS( 2 ), E( J, J+1 ), LDE,
     $                        F( I, J+1 ), LDF )
               END IF
*
   20       CONTINUE
   30    CONTINUE
      ELSE
*
*        Solve transposed (I, J) - system:
*           A(I, I)' * R(I, J) + D(I, I)' * L(J, J) = C(I, J)
*           R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
*        for I = 1, 2, ..., M, J = N, N - 1, ..., 1
*
         SCALE = ONE
         SCALOC = ONE
         DO 80 I = 1, M
            DO 70 J = N, 1, -1
*
*              Build 2 by 2 system Z'
*
               Z( 1, 1 ) = DCONJG( A( I, I ) )
               Z( 2, 1 ) = -DCONJG( B( J, J ) )
               Z( 1, 2 ) = DCONJG( D( I, I ) )
               Z( 2, 2 ) = -DCONJG( E( J, J ) )
*
*
*              Set up right hand side(s)
*
               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )
*
*              Solve Z' * x = RHS
*
               CALL ZGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 )
     $            INFO = IERR
               CALL ZGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 K = 1, N
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ),
     $                           1 )
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ),
     $                           1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
*
*              Unpack solution vector(s)
*
               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )
*
*              Substitute R(I, J) and L(I, J) into remaining equation.
*
               DO 50 K = 1, J - 1
                  F( I, K ) = F( I, K ) + RHS( 1 )*DCONJG( B( K, J ) ) +
     $                        RHS( 2 )*DCONJG( E( K, J ) )
   50          CONTINUE
               DO 60 K = I + 1, M
                  C( K, J ) = C( K, J ) - DCONJG( A( I, K ) )*RHS( 1 ) -
     $                        DCONJG( D( I, K ) )*RHS( 2 )
   60          CONTINUE
*
   70       CONTINUE
   80    CONTINUE
      END IF
      RETURN
*
*     End of ZTGSY2
*
      END
