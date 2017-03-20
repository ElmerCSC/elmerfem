FUNCTION linkFS (Model, nodenumber, Array) RESULT(TopSurface)
  USE DefUtils
  IMPLICIT NONE
  ! in
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: TopSurface, Array(4)
  ! internal
  REAL (KIND=dp) :: Time, FS, X, Y, Z, R, L, H0, N, EPS, a, pr, hpr, lim
  INTEGER :: DIM
  LOGICAL :: FirstTime = .TRUE.

  SAVE FirstTime, DIM


  DIM = CoordinateSystemDimension()
  !PRINT *, "Node=", nodenumber, "A=", Array(1:4)
  Time = Array(1)
  IF (Time == 1) THEN 
     FirstTime = .TRUE.
  ELSE
     FirstTime = .FALSE.
  END IF
  FS = Array(2)
  X = Array(3)
  IF (DIM == 3) THEN
     Y = Array(4)
  ELSE
     Y = 0.0_dp
  END IF
  EPS = 100.0


  IF (FirstTime) THEN
     R = (X**2.0_dp + Y**2.0_dp)**0.5_dp
     L =750000.0_dp
     H0 = 3575.1_dp
     n = 3.0_dp
     IF (R < L) THEN
        TopSurface = H0 * (1.0_dp - (R/L)**((N+1)/N) )**(N/(2.0_dp*N+2.0_dp)) + EPS
        !        TopSurface = H0 * (1.0_dp - (R/l)**() ) + EPS
        !         TopSurface = H0 
     ELSE  
        TopSurface = EPS
     END IF
     ! PRINT *, "X=", X, "Y=", Y, "R=", R, "S=", TopSurface
     FirstTime = .FALSE.
  ELSE
     TopSurface = FS
  END IF
  RETURN
END FUNCTION linkFS
