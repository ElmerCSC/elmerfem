      REAL             FUNCTION SCSUM1( N, CX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX            CX( * )
*     ..
*
*  Purpose
*  =======
*
*  SCSUM1 takes the sum of the absolute values of a complex
*  vector and returns a single precision result.
*
*  Based on SCASUM from the Level 1 BLAS.
*  The change is to use the 'genuine' absolute value.
*
*  Contributed by Nick Higham for use with CLACON.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements in the vector CX.
*
*  CX      (input) COMPLEX array, dimension (N)
*          The vector whose elements will be summed.
*
*  INCX    (input) INTEGER
*          The spacing between successive values of CX.  INCX > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL               STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      SCSUM1 = 0.0E0
      STEMP = 0.0E0
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 )
     $   GO TO 20
*
*     CODE FOR INCREMENT NOT EQUAL TO 1
*
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      SCSUM1 = STEMP
      RETURN
*
*     CODE FOR INCREMENT EQUAL TO 1
*
   20 CONTINUE
      DO 30 I = 1, N
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      SCSUM1 = STEMP
      RETURN
*
*     End of SCSUM1
*
      END
