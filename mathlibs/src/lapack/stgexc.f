      SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, IFST, ILST, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  STGEXC reorders the generalized real Schur decomposition of a real
*  matrix pair (A,B) using an orthogonal equivalence transformation
*
*                 (A, B) = Q * (A, B) * Z',
*
*  so that the diagonal block of (A, B) with row index IFST is moved
*  to row ILST.
*
*  (A, B) must be in generalized real Schur canonical form (as returned
*  by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2
*  diagonal blocks. B is upper triangular.
*
*  Optionally, the matrices Q and Z of generalized Schur vectors are
*  updated.
*
*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)'
*         Q(in) * B(in) * Z(in)' = Q(out) * B(out) * Z(out)'
*
*
*  Arguments
*  =========
*
*  WANTQ   (input) LOGICAL
*          .TRUE. : update the left transformation matrix Q;
*          .FALSE.: do not update Q.
*
*  WANTZ   (input) LOGICAL
*          .TRUE. : update the right transformation matrix Z;
*          .FALSE.: do not update Z.
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the matrix A in generalized real Schur canonical
*          form.
*          On exit, the updated matrix A, again in generalized
*          real Schur canonical form.
*
*  LDA     (input)  INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) REAL array, dimension (LDB,N)
*          On entry, the matrix B in generalized real Schur canonical
*          form (A,B).
*          On exit, the updated matrix B, again in generalized
*          real Schur canonical form (A,B).
*
*  LDB     (input)  INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  Q       (input/output) REAL array, dimension (LDZ,N)
*          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.
*          On exit, the updated matrix Q.
*          If WANTQ = .FALSE., Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= 1.
*          If WANTQ = .TRUE., LDQ >= N.
*
*  Z       (input/output) REAL array, dimension (LDZ,N)
*          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.
*          On exit, the updated matrix Z.
*          If WANTZ = .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= 1.
*          If WANTZ = .TRUE., LDZ >= N.
*
*  IFST    (input/output) INTEGER
*  ILST    (input/output) INTEGER
*          Specify the reordering of the diagonal blocks of (A, B).
*          The block with row index IFST is moved to row ILST, by a
*          sequence of swapping between adjacent blocks.
*          On exit, if IFST pointed on entry to the second row of
*          a 2-by-2 block, it is changed to point to the first row;
*          ILST always points to the first row of the block in its
*          final position (which may differ from its input value by
*          +1 or -1). 1 <= IFST, ILST <= N.
*
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 4*N + 16.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*           =0:  successful exit.
*           <0:  if INFO = -i, the i-th argument had an illegal value.
*           =1:  The transformed matrix pair (A, B) would be too far
*                from generalized Schur form; the problem is ill-
*                conditioned. (A, B) may have been partially reordered,
*                and ILST points to the first row of the current
*                position of the block being moved.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            HERE, LWMIN, NBF, NBL, NBNEXT
*     ..
*     .. External Subroutines ..
      EXTERNAL           STGEX2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input arguments.
*
      INFO = 0
      LWMIN = MAX( 1, 4*N+16 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQ.LT.1 .OR. WANTQ .AND. ( LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. ( LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -12
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STGEXC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the first row of the specified block and find out
*     if it is 1-by-1 or 2-by-2.
*
      IF( IFST.GT.1 ) THEN
         IF( A( IFST, IFST-1 ).NE.ZERO )
     $      IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( A( IFST+1, IFST ).NE.ZERO )
     $      NBF = 2
      END IF
*
*     Determine the first row of the final block
*     and find out if it is 1-by-1 or 2-by-2.
*
      IF( ILST.GT.1 ) THEN
         IF( A( ILST, ILST-1 ).NE.ZERO )
     $      ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( A( ILST+1, ILST ).NE.ZERO )
     $      NBL = 2
      END IF
      IF( IFST.EQ.ILST )
     $   RETURN
*
      IF( IFST.LT.ILST ) THEN
*
*        Update ILST.
*
         IF( NBF.EQ.2 .AND. NBL.EQ.1 )
     $      ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 )
     $      ILST = ILST + 1
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap with next one below.
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1-by-1 or 2-by-2.
*
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( A( HERE+NBF+1, HERE+NBF ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, HERE, NBF, NBNEXT, WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
*
*           Test if 2-by-2 block breaks into two 1-by-1 blocks.
*
            IF( NBF.EQ.2 ) THEN
               IF( A( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1-by-1 blocks, each of which
*           must be swapped individually.
*
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( A( HERE+3, HERE+2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, HERE+1, 1, NBNEXT, WORK, LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1-by-1 blocks.
*
               CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                      LDZ, HERE, 1, 1, WORK, LWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  ILST = HERE
                  RETURN
               END IF
               HERE = HERE + 1
*
            ELSE
*
*              Recompute NBNEXT in case of 2-by-2 split.
*
               IF( A( HERE+2, HERE+1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2-by-2 block did not split.
*
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE, 1, NBNEXT, WORK, LWORK,
     $                         INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
*
*                 2-by-2 block did split.
*
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 1
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 1
               END IF
*
            END IF
         END IF
         IF( HERE.LT.ILST )
     $      GO TO 10
      ELSE
         HERE = IFST
*
   20    CONTINUE
*
*        Swap with next one below.
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1-by-1 or 2-by-2.
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( A( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, HERE-NBNEXT, NBNEXT, NBF, WORK, LWORK,
     $                   INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
*
*           Test if 2-by-2 block breaks into two 1-by-1 blocks.
*
            IF( NBF.EQ.2 ) THEN
               IF( A( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1-by-1 blocks, each of which
*           must be swapped individually.
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( A( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, HERE-NBNEXT, NBNEXT, 1, WORK, LWORK,
     $                   INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1-by-1 blocks.
*
               CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                      LDZ, HERE, NBNEXT, 1, WORK, LWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  ILST = HERE
                  RETURN
               END IF
               HERE = HERE - 1
            ELSE
*
*             Recompute NBNEXT in case of 2-by-2 split.
*
               IF( A( HERE, HERE-1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2-by-2 block did not split.
*
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE-1, 2, 1, WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
*
*                 2-by-2 block did split.
*
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 1
                  CALL STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                         Z, LDZ, HERE, 1, 1, WORK, LWORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 1
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST )
     $      GO TO 20
      END IF
      ILST = HERE
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of STGEXC
*
      END
