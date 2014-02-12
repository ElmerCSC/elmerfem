FUNCTION fbed(x,y) 
USE types
IMPLICIT NONE
REAL(KIND=dp) :: x, y, fbed
REAL(KIND=dp) :: Lxy, fsurf

Lxy = 100.0    

fbed = 10.0*(COS(6.0*Pi*x/Lxy)*COS(6.0*Pi*y/Lxy) -1.0)

END FUNCTION fbed
