PROGRAM runoptf90

  USE elmerparam

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)

  REAL (KIND=dp) :: x(2), res(2)
  INTEGER :: i,j, n

  n = 2
  x(1) = 1.0
  DO i = 1, 3
    x(2) = 1.0
    DO j = 1, 2
      res = elmer_param(n,x)
      x(2) = 2 * x(2) 
    END DO
    x(1) = 4 * x(1)
  END DO

END PROGRAM runoptf90
