FUNCTION Youngs( Model, n, x ) RESULT( s )
  USE Types
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x,s,s1,s2,s3,xx,yy
  
  xx = Model % Nodes % x(n)
  yy = Model % Nodes % y(n)
  
  s =  1.0d0 / SQRT( (xx-11.0)**2 + (yy-4.9)**2 )
  
END FUNCTION Youngs


FUNCTION InFlow( Model, n, x ) RESULT( vin )
  USE Types
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: yy,x,vin,v0,vt
  
  xx = Model % Nodes % x(n)
  yy = Model % Nodes % y(n)
  
  v0 = (-(yy**2)+10*yy)/25
  
  IF(x < 8.0) THEN
    vt = x/8.0
  ELSE
    vt = 1.0
  END IF
  
  vin = v0*vt
END FUNCTION InFlow

