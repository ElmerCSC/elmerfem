      SUBROUTINE ZLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASYF computes a partial factorization of a complex symmetric matrix
*  A using the Bunch-Kaufman diagonal pivoting method. The partial
*  factorization has the form:
*
*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
*        ( 0  U22 ) (  0   D  ) ( U12' U22' )
*
*  A  =  ( L11  0 ) ( D    0  ) ( L11' L21' )  if UPLO = 'L'
*        ( L21  I ) ( 0   A22 ) (  0    I   )
*
*  where the order of D is at most NB. The actual order is returned in
*  the argument KB, and is either NB or NB-1, or N if N <= NB.
*  Note that U' denotes the transpose of U.
*
*  ZLASYF is an auxiliary routine called by ZSYTRF. It uses blocked code
*  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
*  A22 (if UPLO = 'L').
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NB      (input) INTEGER
*          The maximum number of columns of the matrix A that should be
*          factored.  NB should be at least 2 to allow for 2-by-2 pivot
*          blocks.
*
*  KB      (output) INTEGER
*          The number of columns of A that were actually factored.
*          KB is either NB-1 or NB, or N if N <= NB.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, A contains details of the partial factorization.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If UPLO = 'U', only the last KB elements of IPIV are set;
*          if UPLO = 'L', only the first KB elements are set.
*
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  W       (workspace) COMPLEX*16 array, dimension (LDW,NB)
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W.  LDW >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*               has been completed, but the block diagonal matrix D is
*               exactly singular.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     $                   KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, ROWMAX
      COMPLEX*16         D11, D21, D22, R1, T, Z
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      EXTERNAL           LSAME, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZGEMV, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Initialize ALPHA for use in choosing pivot block size.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Factorize the trailing columns of A using the upper triangle
*        of A and working backwards, and compute the matrix W = U12*D
*        for use in updating A11
*
*        K is the main loop index, decreasing from N in steps of 1 or 2
*
*        KW is the column of W which corresponds to column K of A
*
         K = N
   10    CONTINUE
         KW = NB + K - N
*
*        Exit from loop
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
*
*        Copy column K of A to column KW of W and update it
*
         CALL ZCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL ZGEMV( 'No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA,
     $                  W( K, KW+1 ), LDW, CONE, W( 1, KW ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( W( K, KW ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.GT.1 ) THEN
            IMAX = IZAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = CABS1( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column KW-1 of W and update it
*
               CALL ZCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL ZCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                     W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N )
     $            CALL ZGEMV( 'No transpose', K, N-K, -CONE,
     $                        A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW,
     $                        CONE, W( 1, KW-1 ), 1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = IMAX + IZAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = CABS1( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IZAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, KW-1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column KW-1 of W to column KW
*
                  CALL ZCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
*
*                 interchange rows and columns K-1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            KKW = NB + KK - N
*
*           Updated column KP is already stored in column KKW of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL ZCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL ZCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
*
*              Interchange rows KK and KP in last KK columns of A and W
*
               CALL ZSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL ZSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column KW of W now holds
*
*              W(k) = U(k)*D(k)
*
*              where U(k) is the k-th column of U
*
*              Store U(k) in column k of A
*
               CALL ZCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = CONE / A( K, K )
               CALL ZSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-by-2 pivot block D(k): columns KW and KW-1 of W now
*              hold
*
*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
*
*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
*              of U
*
               IF( K.GT.2 ) THEN
*
*                 Store U(k) and U(k-1) in columns k and k-1 of A
*
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = CONE / ( D11*D22-CONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KSTEP
         GO TO 10
*
   30    CONTINUE
*
*        Update the upper triangle of A11 (= A(1:k,1:k)) as
*
*        A11 := A11 - U12*D*U12' = A11 - U12*W'
*
*        computing blocks of NB columns at a time
*
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
*
*           Update the upper triangle of the diagonal block
*
            DO 40 JJ = J, J + JB - 1
               CALL ZGEMV( 'No transpose', JJ-J+1, N-K, -CONE,
     $                     A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, CONE,
     $                     A( J, JJ ), 1 )
   40       CONTINUE
*
*           Update the rectangular superdiagonal block
*
            CALL ZGEMM( 'No transpose', 'Transpose', J-1, JB, N-K,
     $                  -CONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW,
     $                  CONE, A( 1, J ), LDA )
   50    CONTINUE
*
*        Put U12 in standard form by partially undoing the interchanges
*        in columns k+1:n
*
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF( JP.NE.JJ .AND. J.LE.N )
     $      CALL ZSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LE.N )
     $      GO TO 60
*
*        Set KB to the number of columns factorized
*
         KB = N - K
*
      ELSE
*
*        Factorize the leading columns of A using the lower triangle
*        of A and working forwards, and compute the matrix W = L21*D
*        for use in updating A22
*
*        K is the main loop index, increasing from 1 in steps of 1 or 2
*
         K = 1
   70    CONTINUE
*
*        Exit from loop
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
*
*        Copy column K of A to column K of W and update it
*
         CALL ZCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA,
     $               W( K, 1 ), LDW, CONE, W( K, K ), 1 )
*
         KSTEP = 1
*
*        Determine rows and columns to be interchanged and whether
*        a 1-by-1 or 2-by-2 pivot block will be used
*
         ABSAKK = CABS1( W( K, K ) )
*
*        IMAX is the row-index of the largest off-diagonal element in
*        column K, and COLMAX is its absolute value
*
         IF( K.LT.N ) THEN
            IMAX = K + IZAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = CABS1( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           Column K is zero: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              no interchange, use 1-by-1 pivot block
*
               KP = K
            ELSE
*
*              Copy column IMAX to column K+1 of W and update it
*
               CALL ZCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL ZCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ),
     $                     1 )
               CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ),
     $                     LDA, W( IMAX, 1 ), LDW, CONE, W( K, K+1 ),
     $                     1 )
*
*              JMAX is the column-index of the largest off-diagonal
*              element in row IMAX, and ROWMAX is its absolute value
*
               JMAX = K - 1 + IZAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = CABS1( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IZAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, K+1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 no interchange, use 1-by-1 pivot block
*
                  KP = K
               ELSE IF( CABS1( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 interchange rows and columns K and IMAX, use 1-by-1
*                 pivot block
*
                  KP = IMAX
*
*                 copy column K+1 of W to column K
*
                  CALL ZCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
*
*                 interchange rows and columns K+1 and IMAX, use 2-by-2
*                 pivot block
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
*
*           Updated column KP is already stored in column KK of W
*
            IF( KP.NE.KK ) THEN
*
*              Copy non-updated column KK to column KP
*
               A( KP, K ) = A( KK, K )
               CALL ZCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL ZCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
*
*              Interchange rows KK and KP in first KK columns of A and W
*
               CALL ZSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL ZSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-by-1 pivot block D(k): column k of W now holds
*
*              W(k) = L(k)*D(k)
*
*              where L(k) is the k-th column of L
*
*              Store L(k) in column k of A
*
               CALL ZCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = CONE / A( K, K )
                  CALL ZSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-by-2 pivot block D(k): columns k and k+1 of W now hold
*
*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
*
*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
*              of L
*
               IF( K.LT.N-1 ) THEN
*
*                 Store L(k) and L(k+1) in columns k and k+1 of A
*
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = CONE / ( D11*D22-CONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
*
*              Copy D(k) to A
*
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
         END IF
*
*        Store details of the interchanges in IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KSTEP
         GO TO 70
*
   90    CONTINUE
*
*        Update the lower triangle of A22 (= A(k:n,k:n)) as
*
*        A22 := A22 - L21*D*L21' = A22 - L21*W'
*
*        computing blocks of NB columns at a time
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           Update the lower triangle of the diagonal block
*
            DO 100 JJ = J, J + JB - 1
               CALL ZGEMV( 'No transpose', J+JB-JJ, K-1, -CONE,
     $                     A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, CONE,
     $                     A( JJ, JJ ), 1 )
  100       CONTINUE
*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                     K-1, -CONE, A( J+JB, 1 ), LDA, W( J, 1 ),
     $                     LDW, CONE, A( J+JB, J ), LDA )
  110    CONTINUE
*
*        Put L21 in standard form by partially undoing the interchanges
*        in columns 1:k-1
*
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF( JP.NE.JJ .AND. J.GE.1 )
     $      CALL ZSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GE.1 )
     $      GO TO 120
*
*        Set KB to the number of columns factorized
*
         KB = K - 1
*
      END IF
      RETURN
*
*     End of ZLASYF
*
      END
