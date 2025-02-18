FUNCTION FunA( Model, n, t ) RESULT( s )
  USE DefUtils
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,s
  
  s = t
  
END FUNCTION FunA


FUNCTION FunB( Model, n, t ) RESULT( s )
  USE DefUtils
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,s
  REAL(KIND=dp) :: xx,yy
  
  xx = Model % Nodes % x(n)
  yy = Model % Nodes % y(n)

  s = xx + yy
  
END FUNCTION FunB


FUNCTION FunC( Model, n, x  )  RESULT ( s ) 
  USE DefUtils
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x(*), s

  s = x(1)**2 + x(2)**2

END FUNCTION FunC




