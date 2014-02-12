      DOUBLE PRECISION   FUNCTION PDLAMCH( ICTXT, CMACH )
      include "mpif.h"
*
*  -- ScaLAPACK auxilliary routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 1995
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
      INTEGER            ICTXT
*     ..
*
*  Purpose
*  =======
*
*  PDLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle in which the computation takes
*          place.
*
*  CMACH   (global input) CHARACTER*1
*          Specifies the value to be returned by PDLAMCH:
*          = 'E' or 'e',   PDLAMCH := eps
*          = 'S' or 's ,   PDLAMCH := sfmin
*          = 'B' or 'b',   PDLAMCH := base
*          = 'P' or 'p',   PDLAMCH := eps*base
*          = 'N' or 'n',   PDLAMCH := t
*          = 'R' or 'r',   PDLAMCH := rnd
*          = 'M' or 'm',   PDLAMCH := emin
*          = 'U' or 'u',   PDLAMCH := rmin
*          = 'L' or 'l',   PDLAMCH := emax
*          = 'O' or 'o',   PDLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IDUMM
      DOUBLE PRECISION   TEMP, TEMP1
*     ..
*     .. External Subroutines ..
*      EXTERNAL           DGAMN2D, DGAMX2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, LSAME
*     ..
*     .. Executable Statements ..
*
      TEMP1 = DLAMCH( CMACH )
*
      IF( LSAME( CMACH, 'E' ).OR.LSAME( CMACH, 'S' ).OR.
     $    LSAME( CMACH, 'M' ).OR.LSAME( CMACH, 'U' ) ) THEN
          CALL MPI_ALLREDUCE( TEMP1, TEMP, 1, MPI_DOUBLE_PRECISION,
     $                        MPI_MAX, ICTXT, IDUMM )
*         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, TEMP, 1, IDUMM,
*     $                 IDUMM, 1, -1, IDUMM )
      ELSE IF( LSAME( CMACH, 'L' ).OR.LSAME( CMACH, 'O' ) ) THEN
          CALL MPI_ALLREDUCE( TEMP1, TEMP, 1, MPI_DOUBLE_PRECISION,
     $                        MPI_MIN, ICTXT, IDUMM )
*         CALL DGAMN2D( ICTXT, 'All', ' ', 1, 1, TEMP, 1, IDUMM,
*     $                 IDUMM, 1, -1, IDUMM )
      ELSE
          TEMP = TEMP1
      END IF
*
      PDLAMCH = TEMP
*
*     End of PDLAMCH
*
      END
