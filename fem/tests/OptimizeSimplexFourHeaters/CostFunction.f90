!------------------------------------------------------------------------------
! Compute the cost function value when the desired temperature for each heater
! is known to be 20, 30, 40 and 50 degs.
!------------------------------------------------------------------------------
FUNCTION CostFunction( Model, n, x ) RESULT( s )

  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE SolverUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x,s

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Variable_t), POINTER :: Var
  REAL (KIND=DP) :: Cost, Tmax
  REAL(KIND=dp) :: x0(4), y0(4), xx, yy
  INTEGER :: i0(4), i
  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp), POINTER :: Temp(:)
  INTEGER, POINTER :: TempPerm(:)
  

  SAVE i0, Visited

  Var => VariableGet( Model % Variables, 'Temperature', .TRUE.)
  IF ( .NOT. ASSOCIATED( Var ) )  THEN
    CALL FATAL('CostFunction','variable Temperature not found')
  ELSE
    Temp => Var % Values
    TempPerm => Var % Perm
  END IF

  IF(.NOT. Visited) THEN
    x0 = Model % Nodes % x(1)
    y0 = Model % Nodes % y(1)
    i0 = 1
    
    DO i=1,Model % NumberOfNodes
      xx = Model % Nodes % x(i)
      yy = Model % Nodes % y(i)
      IF((xx-2.0)**2.0+(yy-2.0)**2.0 < (x0(1)-2.0)**2.0 + (y0(1)-2.0)**2.0) THEN
        x0(1) = xx
        y0(1) = yy
        i0(1) = i
      END IF
      IF((xx-6.0)**2.0+(yy-2.0)**2.0 < (x0(2)-6.0)**2.0 + (y0(2)-2.0)**2.0) THEN
        x0(2) = xx
        y0(2) = yy
        i0(2) = i
      END IF
      IF((xx-2.0)**2.0+(yy-6.0)**2.0 < (x0(3)-2.0)**2.0 + (y0(3)-6.0)**2.0) THEN
        x0(3) = xx
        y0(3) = yy
        i0(3) = i
      END IF
      IF((xx-6.0)**2.0+(yy-6.0)**2.0 < (x0(4)-6.0)**2.0 + (y0(4)-6.0)**2.0) THEN
        x0(4) = xx
        y0(4) = yy
        i0(4) = i
      END IF
    END DO
    
    PRINT *,'Node indexes',i0
    PRINT *,'X-coords',x0
    PRINT *,'Y-coords',y0
    Visited = .TRUE.

    PRINT *,'T1',Temp(TempPerm(i0(1)))
    PRINT *,'T2',Temp(TempPerm(i0(2)))
    PRINT *,'T3',Temp(TempPerm(i0(3)))
    PRINT *,'T4',Temp(TempPerm(i0(4)))
  END IF

  Cost = ABS(Temp(TempPerm(i0(1)))-20.0) + ABS(Temp(TempPerm(i0(2)))-30.0) + &
      ABS(Temp(TempPerm(i0(3)))-40.0) + ABS(Temp(TempPerm(i0(4)))-50.0)


  s = Cost

!------------------------------------------------------------------------------
END FUNCTION CostFunction
!------------------------------------------------------------------------------


