
  FUNCTION InitDist( Model, n, t ) RESULT(f)
    USE Types
    USE SolverUtils
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    LOGICAL :: Visited = .FALSE., GotIt
    INTEGER :: n
    REAL(KIND=dp) :: t,x,y,a1,a2,f
    
    SAVE Visited, a1, a2

    IF(.NOT. Visited) THEN
      a1 = 5.0_dp
      a2 = 10.0_dp
      Visited = .TRUE.
    END IF
    
    f = 0.0_dp

    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    

    ! cut-off Gaussian
    ! f = MIN(1.0_dp, a1*EXP(-a2*(x**2+y**2)))

    ! sharp cylinder
    IF( (x-0.5)**2 + y**2 < 0.3**2 ) f = 1.0_dp

    ! sharp square
    ! IF( ABS(x-0.5) < 0.3 .AND. ABS(y) < 0.3 ) f = 1.0_dp


  END FUNCTION InitDist
