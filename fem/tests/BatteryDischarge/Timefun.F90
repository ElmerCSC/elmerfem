FUNCTION TimeFun(Model, n, tind) RESULT(dt)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: tind, dt

  LOGICAL :: Visited = .FALSE.,maxReached=.FALSE.
  REAL(KIND=dp) :: dtmax, dt0, q

  SAVE Visited, MaxReached, dtmax, dt0, q

  IF(.NOT. Visited ) THEN
    dtmax = ListGetConstReal( Model % Simulation,'Max Timestep')
    dt0 = ListGetConstReal( Model % Simulation,'First Timestep')
    q = ListGetConstReal( Model % Simulation,'Timestep Ratio')
  END IF

  ! Since q^(tind-1) minght be huge we use constant timestep after reaching the max
  IF( maxReached ) THEN
    dt = dtmax
    RETURN
  END IF

  ! Timestep growing exponetially
  dt = dt0 * q**(tind-1.0_dp)

  ! Until max timestep is reached
  IF( dt > dtmax ) THEN
    dt = dtmax
    maxReached = .TRUE.
  END IF
  
END FUNCTION TimeFun
