FUNCTION fsurf(x,y) 
USE types
IMPLICIT NONE
REAL(KIND=dp) :: x, y, fsurf

fsurf = - x * TAN(0.5d0*Pi/180.0d0)  

END FUNCTION fsurf
