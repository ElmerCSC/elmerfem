      SUBROUTINE SLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, L, LDC, M, N
      REAL               TAU
*     ..
*     .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLARZ applies a real elementary reflector H to a real M-by-N
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*
*  H is a product of k elementary reflectors as returned by STZRZF.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  L       (input) INTEGER
*          The number of entries of the vector V containing
*          the meaningful part of the Householder vectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  V       (input) REAL array, dimension (1+(L-1)*abs(INCV))
*          The vector v in the representation of H as returned by
*          STZRZF. V is not used if TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) REAL
*          The value tau in the representation of H.
*
*  C       (input/output) REAL array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) REAL array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMV, SGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w( 1:n ) = C( 1, 1:n )
*
            CALL SCOPY( N, C, LDC, WORK, 1 )
*
*           w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )' * v( 1:l )
*
            CALL SGEMV( 'Transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V,
     $                  INCV, ONE, WORK, 1 )
*
*           C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
*
            CALL SAXPY( N, -TAU, WORK, 1, C, LDC )
*
*           C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
*                               tau * v( 1:l ) * w( 1:n )'
*
            CALL SGER( L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ),
     $                 LDC )
         END IF
*
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w( 1:m ) = C( 1:m, 1 )
*
            CALL SCOPY( M, C, 1, WORK, 1 )
*
*           w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
*
            CALL SGEMV( 'No transpose', M, L, ONE, C( 1, N-L+1 ), LDC,
     $                  V, INCV, ONE, WORK, 1 )
*
*           C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
*
            CALL SAXPY( M, -TAU, WORK, 1, C, 1 )
*
*           C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
*                               tau * w( 1:m ) * v( 1:l )'
*
            CALL SGER( M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ),
     $                 LDC )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of SLARZ
*
      END
