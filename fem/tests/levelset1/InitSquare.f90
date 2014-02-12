
  FUNCTION InitSquare( Model, n, t ) RESULT(f)
    USE Types
    USE SolverUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n,i,j
    REAL(KIND=dp) :: t,x,y,r,Lx,Ly,dist(4),f,a1,a2,b1,b2
    
    SAVE Visited, a1, a2, b1, b2
    
    IF(.NOT. Visited) THEN
      Lx = 1.0
      Ly = 1.0
      
      a1 = 0.25 * Lx 
      a2 = 0.75 * Lx
      b1 = 0.25 * Ly
      b2 = 0.75 * Ly
      Visited = .TRUE.
    END IF
    
    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
    
    ! distance from horizontal lines
    IF( x < a1 ) THEN
      dist(1) = SQRT( (x-a1)**2.0 + (y-b1)**2.0)
      dist(2) = SQRT( (x-a1)**2.0 + (y-b2)**2.0)
    ELSE IF( x > a2 ) THEN     
      dist(1) = SQRT( (x-a2)**2.0 + (y-b1)**2.0)
      dist(2) = SQRT( (x-a2)**2.0 + (y-b2)**2.0)
    ELSE     
      dist(1) = ABS(y - b1)
      dist(2) = ABS(y - b2)
    END IF
    
    ! distance from vertical lines
    IF( y < b1 ) THEN
      dist(3) = SQRT( (x-a1)**2.0 + (y-b1)**2.0)
      dist(4) = SQRT( (x-a2)**2.0 + (y-b1)**2.0)
    ELSE IF( y > b2 ) THEN     
      dist(3) = SQRT( (x-a1)**2.0 + (y-b2)**2.0)
      dist(4) = SQRT( (x-a2)**2.0 + (y-b2)**2.0)
    ELSE     
      dist(3) = ABS(x - a1)
      dist(4) = ABS(x - a2)
    END IF
    
    ! minimum distance to the sides of the rectangle
    f = MINVAL( ABS( dist ) )
    
    ! set negative sign inside the rectangle 
    IF( .NOT. (x > a1 .AND. x < a2 .AND. y > b1 .AND. y < b2) ) THEN
      f = -f
    END IF

  END FUNCTION InitSquare
