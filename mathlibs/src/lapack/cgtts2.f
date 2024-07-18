      SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            ITRANS, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
*     ..
*
*  Purpose
*  =======
*
*  CGTTS2 solves one of the systems of equations
*     A * X = B,  A**T * X = B,  or  A**H * X = B,
*  with a tridiagonal matrix A using the LU factorization computed
*  by CGTTRF.
*
*  Arguments
*  =========
*
*  ITRANS  (input) INTEGER
*          Specifies the form of the system of equations.
*          = 0:  A * X = B     (No transpose)
*          = 1:  A**T * X = B  (Transpose)
*          = 2:  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) COMPLEX array, dimension (N-1)
*          The (n-1) multipliers that define the matrix L from the
*          LU factorization of A.
*
*  D       (input) COMPLEX array, dimension (N)
*          The n diagonal elements of the upper triangular matrix U from
*          the LU factorization of A.
*
*  DU      (input) COMPLEX array, dimension (N-1)
*          The (n-1) elements of the first super-diagonal of U.
*
*  DU2     (input) COMPLEX array, dimension (N-2)
*          The (n-2) elements of the second super-diagonal of U.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row IPIV(i).  IPIV(i) will always be either
*          i or i+1; IPIV(i) = i indicates a row interchange was not
*          required.
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the matrix of right hand side vectors B.
*          On exit, B is overwritten by the solution vectors X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX            TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( ITRANS.EQ.0 ) THEN
*
*        Solve A*X = B using the LU factorization of A,
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
   10       CONTINUE
*
*           Solve L*x = b.
*
            DO 20 I = 1, N - 1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               ELSE
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               END IF
   20       CONTINUE
*
*           Solve U*x = b.
*
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 )
     $         B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                       D( N-1 )
            DO 30 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )*
     $                     B( I+2, J ) ) / D( I )
   30       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 10
            END IF
         ELSE
            DO 60 J = 1, NRHS
*
*           Solve L*x = b.
*
               DO 40 I = 1, N - 1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
                  ELSE
                     TEMP = B( I, J )
                     B( I, J ) = B( I+1, J )
                     B( I+1, J ) = TEMP - DL( I )*B( I, J )
                  END IF
   40          CONTINUE
*
*           Solve U*x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               IF( N.GT.1 )
     $            B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                          D( N-1 )
               DO 50 I = N - 2, 1, -1
                  B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )*
     $                        B( I+2, J ) ) / D( I )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ITRANS.EQ.1 ) THEN
*
*        Solve A**T * X = B.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
   70       CONTINUE
*
*           Solve U**T * x = b.
*
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 )
     $         B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 80 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )*
     $                     B( I-2, J ) ) / D( I )
   80       CONTINUE
*
*           Solve L**T * x = b.
*
            DO 90 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               END IF
   90       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 70
            END IF
         ELSE
            DO 120 J = 1, NRHS
*
*           Solve U**T * x = b.
*
               B( 1, J ) = B( 1, J ) / D( 1 )
               IF( N.GT.1 )
     $            B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
               DO 100 I = 3, N
                  B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-
     $                        DU2( I-2 )*B( I-2, J ) ) / D( I )
  100          CONTINUE
*
*           Solve L**T * x = b.
*
               DO 110 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                     B( I, J ) = TEMP
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
*
*        Solve A**H * X = B.
*
         IF( NRHS.LE.1 ) THEN
            J = 1
  130       CONTINUE
*
*           Solve U**H * x = b.
*
            B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
            IF( N.GT.1 )
     $         B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) /
     $                     CONJG( D( 2 ) )
            DO 140 I = 3, N
               B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )*B( I-1, J )-
     $                     CONJG( DU2( I-2 ) )*B( I-2, J ) ) /
     $                     CONJG( D( I ) )
  140       CONTINUE
*
*           Solve L**H * x = b.
*
            DO 150 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - CONJG( DL( I ) )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                  B( I, J ) = TEMP
               END IF
  150       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 130
            END IF
         ELSE
            DO 180 J = 1, NRHS
*
*           Solve U**H * x = b.
*
               B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
               IF( N.GT.1 )
     $            B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) /
     $                        CONJG( D( 2 ) )
               DO 160 I = 3, N
                  B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )*
     $                        B( I-1, J )-CONJG( DU2( I-2 ) )*
     $                        B( I-2, J ) ) / CONJG( D( I ) )
  160          CONTINUE
*
*           Solve L**H * x = b.
*
               DO 170 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - CONJG( DL( I ) )*
     $                           B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                     B( I, J ) = TEMP
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
*
*     End of CGTTS2
*
      END
