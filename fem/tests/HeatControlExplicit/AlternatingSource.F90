! This is a user defined source that uses the center point
! to alternative the sign of the source. 
!---------------------------------------------------------------
FUNCTION AlternatingSource( Model, n, x ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x(2), f

  REAL(KIND=dp) :: tlim = 0.07, tcenter, time
  INTEGER :: sgn = 1
  
  SAVE sgn
  
  tcenter = x(1)
  time = x(2)

  ! Switch sgn of source term when we go beyond [-tlim,tlim]
  IF( sgn == 1 .AND. tcenter > tlim ) THEN
    IF(ParEnv % MyPe == 0) PRINT *,'Neg source:',time,tcenter
    sgn = -1
  END IF
  IF( sgn == -1 .AND. tcenter < -tlim ) THEN
    IF(ParEnv % MyPe == 0) PRINT *,'Pos source:',time,tcenter
    sgn = 1
  END IF
    
  f = sgn * 1.0_dp 
    
END FUNCTION AlternatingSource

