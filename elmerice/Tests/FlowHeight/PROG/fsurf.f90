FUNCTION fsurf(x,y) 
USE types
IMPLICIT NONE
REAL(KIND=dp) :: x, y, fsurf, Lxy

Lxy = 2.0e3

fsurf =  100.0*SIN(2.0*pi*x/Lxy)

END FUNCTION fsurf
