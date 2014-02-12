!-------------------------------------------------------------
! Gives the upper limit as determined by a sphere in a fixed
! position: R^2=x^2+y^2+(z-H)^2.
!-------------------------------------------------------------
  FUNCTION SphereBottom( Model, n, t ) RESULT(f)
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t,f
!-------------------------------------------------------------
    REAL(KIND=dp) :: x,y,z,R,Disc,H

    x = Model % Nodes % x(n)
    y = Model % Nodes % y(n)    
    z = Model % Nodes % z(n)
    R = 3.0_dp
    
    R = ListGetCReal( Model % Simulation,'Sphere Radius')
    H = ListGetCReal( Model % Simulation,'Sphere Z')

    Disc = R**2 - x**2 - y**2
    IF( Disc < 0 ) THEN
      f = HUGE(f)
    ELSE
      f = H - SQRT(Disc)
!      PRINT *,'Contact:',x,y,z,H,f
    END IF

  END FUNCTION SphereBottom
!-------------------------------------------------------------

