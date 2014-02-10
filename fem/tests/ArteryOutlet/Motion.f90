FUNCTION Motion( Model,n,t ) RESULT(uy)
USE Types

IMPLICIT NONE

TYPE(Model_t) :: Model

INTEGER :: n
DOUBLE PRECISION :: t, uy, x, umax, tpulse

x = Model % Nodes % x(n)

! Total displacement
umax = 0.05

! Duration of the motion (umax)
tpulse = 0.6

x = Model % Nodes % x(n)

IF (t < tpulse/2) THEN
   uy = (umax/2)*(SIN((2*t-tpulse/2)/tpulse*pi)+1)
ELSE
   uy = umax
END IF

END FUNCTION Motion
