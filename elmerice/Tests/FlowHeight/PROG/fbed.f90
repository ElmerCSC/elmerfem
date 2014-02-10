FUNCTION fbed(x,y) 
USE types
IMPLICIT NONE
REAL(KIND=dp) :: x, y, fbed
REAL(KIND=dp) :: Lxy, fsurf

Lxy = 2.0e3_dp

fbed =  -1000.0 + 100.0*SIN(4*2.0*pi*x/Lxy)

END FUNCTION fbed
