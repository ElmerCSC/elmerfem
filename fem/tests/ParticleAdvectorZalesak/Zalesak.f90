!------------------------------------------------------
! Pijl et al, 4.1.2. Zalesak's rotating disk
!------------------------------------------------------

  FUNCTION InitZalesak( Model, n, t ) RESULT(f)
    USE DefUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n,i,j
    REAL(KIND=dp) :: t,x,y,r,Lx,Ly,x0,y0,r0,w0,dist(4),f,a1,a2,b1,b2,f1,f2
    
    SAVE Visited, x0, y0, r0, w0, a1, a2, b1, b2
    
    IF(.NOT. Visited) THEN
      Lx = ListGetConstReal(Model % Simulation,'Zalesak Lx',GotIt)
      IF(.NOT. GotIt) THEN
        Lx = MAXVAL( Model % Nodes % x )
        Lx = ParallelReduction( Lx ) 
      END IF
      Ly = ListGetConstReal(Model % Simulation,'Zalesak Ly',GotIt)
      IF(.NOT. GotIt) THEN
        Ly = MAXVAL( Model % Nodes % y ) 
        Ly = ParallelReduction( Ly )
      END IF
      
      x0 = 0.5 * Lx
      y0 = 0.75 * Ly
      r0 = (3.0/20.0) * Lx
      w0 = r0 / 3.0d0

      a1 = x0 - w0/2
      a2 = x0 + w0/2
      b1 = -1.0
      b2 = y0 + r0 - w0
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
    f1 = MINVAL( ABS( dist ) )
    
    ! set negative sign inside the rectangle 
    IF( x > a1 .AND. x < a2 .AND. y > b1 .AND. y < b2) THEN
      f1 = -f1
    END IF
      
    ! inside of a circle
    f2 = r0 - SQRT( (x-x0)**2.0 + (y-y0)**2.0 )

    ! take a union of the rectangle and circle
    IF( f2 > 0.0 ) THEN
      IF( f1 < 0.0 ) THEN
        f = f1
      ELSE
        f = MIN(f1,f2)
      END IF
    ELSE 
      IF( f1 > 0.0 ) THEN
        f = f2
      ELSE 
        f = MIN(f1,f2)
      END IF
    END IF
     
  END FUNCTION InitZalesak



  FUNCTION RotateVeloX( Model, n, t ) RESULT(f)
    USE Types
    USE SolverUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n,i,j
    REAL(KIND=dp) :: t,x,y,r,Lx,Ly,x0,y0,f
    
    SAVE Visited, x0, y0
    
    IF(.NOT. Visited) THEN
      Lx = ListGetConstReal(Model % Simulation,'Zalesak Lx',GotIt)
      IF(.NOT. GotIt) THEN
        Lx = MAXVAL( Model % Nodes % x )
        Lx = ParallelReduction( Lx ) 
      END IF
      Ly = ListGetConstReal(Model % Simulation,'Zalesak Ly',GotIt)
      IF(.NOT. GotIt) THEN
        Ly = MAXVAL( Model % Nodes % y ) 
        Ly = ParallelReduction( Ly )
      END IF
      x0 = Lx / 2.0
      y0 = Ly / 2.0
      Visited = .TRUE.
    END IF
    
    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    

    f = -(y-y0)

  END FUNCTION RotateVeloX  


  FUNCTION RotateVeloY( Model, n, t ) RESULT(f)
    USE Types
    USE SolverUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n,i,j
    REAL(KIND=dp) :: t,x,y,r,Lx,Ly,x0,y0,f
    
    SAVE Visited, x0, y0
    
    IF(.NOT. Visited) THEN
      Lx = ListGetConstReal(Model % Simulation,'Zalesak Lx',GotIt)
      IF(.NOT. GotIt) THEN
        Lx = MAXVAL( Model % Nodes % x )
        Lx = ParallelReduction( Lx ) 
      END IF
      Ly = ListGetConstReal(Model % Simulation,'Zalesak Ly',GotIt)
      IF(.NOT. GotIt) THEN
        Ly = MAXVAL( Model % Nodes % y ) 
        Ly = ParallelReduction( Ly )
      END IF
      x0 = Lx / 2.0
      y0 = Ly / 2.0
      Visited = .TRUE.
    END IF
    
    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    

    f = x - x0

  END FUNCTION RotateVeloY


