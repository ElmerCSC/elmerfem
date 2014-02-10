MODULE Blowers

  USE Types
  USE DefUtils
  IMPLICIT NONE

  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp) :: x1,y1,r1,x2,y2,r2,vx1,vy1,vx2,vy2
  
  SAVE :: Visited,x1,y1,r1,x2,y2,r2,vx1,vy1,vx2,vy2


CONTAINS

  FUNCTION InsideBlower(Model,n) RESULT ( BlowerInd ) 
    INTEGER :: n 
    INTEGER :: BlowerInd
    TYPE(Model_t) :: Model

    REAL(KIND=dp) :: x,y
    TYPE(Valuelist_t), POINTER :: List

    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
   
    IF(.NOT. Visited ) THEN
      List => Model % Simulation 
      x1 = GetCReal( List,'Blower 1 x')
      y1 = GetCReal( List,'Blower 1 y')
      r1 = GetCReal( List,'Blower 1 r')
      vx1 = GetCReal( List,'Blower 1 vx')
      vy1 = GetCReal( List,'Blower 1 vy')
      PRINT *,'Blower 1:',x1,y1,r1,vx1,vy1

      x2 = GetCReal( List,'Blower 2 x')
      y2 = GetCReal( List,'Blower 2 y')
      r2 = GetCReal( List,'Blower 2 r')
      vx2 = GetCReal( List,'Blower 2 vx')
      vy2 = GetCReal( List,'Blower 2 vy')
      PRINT *,'Blower 1:',x2,y2,r2,vx2,vy2
    
      Visited = .TRUE.
    END IF

    
    IF( (x-x1)**2 + (y-y1)**2 < r1**2 ) THEN
      BlowerInd = 1
    ELSE IF( (x-x2)**2 + (y-y2)**2 < r2**2 ) THEN
      BlowerInd = 2 
    ELSE
      BlowerInd = 0
    END IF

  END FUNCTION InsideBlower


END MODULE Blowers


FUNCTION BlowerVelox( Model, n, t ) RESULT(f)
  USE Types
  USE DefUtils
  USE Blowers
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f
  
  INTEGER :: ind
  
  ind = InsideBlower(Model, n)
  IF( ind == 1 ) THEN
    f = vx1
  ELSE IF( ind == 2 ) THEN
    f = vx2 
  ELSE
    f = 0.0_dp
  END IF
END FUNCTION BlowerVelox


FUNCTION BlowerVeloy( Model, n, t ) RESULT(f)
  USE Types
  USE DefUtils
  USE Blowers
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f
  
  INTEGER :: ind
  
  ind = InsideBlower(Model, n)
  IF( ind == 1 ) THEN
    f = vy1
  ELSE IF( ind == 2 ) THEN
    f = vy2 
  ELSE
    f = 0.0_dp
  END IF
END FUNCTION BlowerVeloy


FUNCTION BlowerVeloCond( Model, n, t ) RESULT(f)
  USE Types
  USE DefUtils
  USE Blowers
  
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t,f
  
  INTEGER :: ind
  
  ind = InsideBlower(Model, n)
  IF( ind > 0 ) THEN
    f = 1.0_dp
  ELSE
    f = -1.0_dp
  END IF
END FUNCTION BlowerVeloCond
