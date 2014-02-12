FUNCTION fbed(x,y) 
USE types
IMPLICIT NONE
REAL(KIND=dp) :: x, y, fbed
REAL(KIND=dp) :: Lxy, fsurf

Lxy = 5.0e3_dp

fbed = fsurf(x,y) -1000.0d0*(1.0 - 0.5d0*SIN(2.0*Pi*x/Lxy)*SIN(2.0*Pi*y/Lxy))

END FUNCTION fbed
