FUNCTION CatenaryCurve( Model, n, x ) RESULT(y)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x, y

  
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Visited = .FALSE.,Found
  REAL(KIND=dp) :: x1,x2,dy,dx,xm,fa,dfa,a,da
  INTEGER :: i
  
  SAVE :: Mesh,xm,dx,a,dy
  

  y = 0._dp

  IF(.NOT. Visited ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation

    x1 = GetCReal(Params,'Catenary Start')
    x2 = GetCReal(Params,'Catenary End')    
    dy = GetCReal(Params,'Catenary Height')

    dx = (x2-x1) / 2
    xm = (x1+x2) / 2

    ! Newton's iteration
    a = dx    
    DO i=1,100
      fa = a*(COSH(dx/a)-1)-dy
      dfa = COSH(dx/a)-1-(dx/a)*SINH(dx/a)
      da = -fa/dfa
      a = a + da
      PRINT *,'a_i:',i,a,da,fa,dfa,fa+dy
      IF(ABS(da) < 1.0d-10 ) EXIT
    END DO
    
    Visited = .TRUE.
  END IF

  ! default displacement
  y = 0.0_dp  

  ! curve is symmetric around mid point
  x = Mesh % Nodes % x(n) - xm
  IF( ABS(x) > dx ) RETURN

  y = a*(COSH(x/a)-COSH(dx/a))
    
END FUNCTION CatenaryCurve

