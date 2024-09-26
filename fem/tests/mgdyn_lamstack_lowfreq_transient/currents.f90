
module ToroidCurrents2
  contains

    function rho(r, r0, y)
      USE Types
      real(kind=dp) :: rho, r, r0, y

       rho = sqrt((r-r0)**2+y**2)
    end function rho

    function currentInToroidR(r, r0, y, turns, I) result (curr)
      USE Types
      real(kind=dp) :: r, r0, y, curr, rho1, rho2, turns, I

      rho1 = rho(8d-3, r0, y)
      rho2 = rho(10d-3, r0, y)

      curr = turns * I /(2*pi*(rho2-rho1)*r)*(-y/rho(r, r0, y))

    end function currentInToroidR

    function currentInToroidY(r, r0, y, turns, I) result (curr)
      USE Types
      real(kind=dp) :: r, r0, y, curr, rho1, rho2, turns, I

      rho1 = rho(8d-3, r0, y)
      rho2 = rho(10d-3, r0, y)

      curr = turns * I /(2*pi*(rho2-rho1)*r)*((r-r0)/rho(r, r0, y))

    end function currentInToroidY

End module ToroidCurrents2

FUNCTION currdens1( model, n, args) RESULT(curr)
  USE DefUtils
  Use ToroidCurrents2
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: x, y, z, args(4), curr, f, omega, &
                   theta, r0, r, turns, I, t

  x = args(1)
  y = args(2)
  z = args(3)
  t = args(4)

  r0 = 45d-3
  turns = 100d0

  f = 500d0
  omega = 2d0*pi*f

  I = 2d0 * sin(omega * t)

  r = sqrt(x**2 + z**2)

  theta = atan2(x,z)

  curr = currentInToroidR(r, r0, y, turns, I) * sin(theta)

END FUNCTION currdens1

FUNCTION currdens2( model, n, args) RESULT(curr)
  USE DefUtils
  Use ToroidCurrents2
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: x, y, z, args(4), curr, f, omega, &
                   theta, r0, r, turns, I, t

  x = args(1)
  y = args(2)
  z = args(3)
  t = args(4)

  r0 = 45d-3
  turns = 100d0

  f = 500d0
  omega = 2d0*pi*f

  I = 2d0 * sin(omega * t)

  r = sqrt(x**2 + z**2)

  theta = atan2(x,z)

  curr = currentInToroidY(r, r0, y, turns, I)

END FUNCTION currdens2

FUNCTION currdens3( model, n, args) RESULT(curr)
  USE DefUtils
  Use ToroidCurrents2
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: x, y, z, args(4), curr, f, omega, &
                   theta, r0, r, turns, I, t

  x = args(1)
  y = args(2)
  z = args(3)
  t = args(4)

  r0 = 45d-3
  turns = 100d0

  f = 500d0
  omega = 2d0*pi*f

  I = 2d0 * sin(omega * t)

  r = sqrt(x**2 + z**2)

  theta = atan2(x,z)

  curr = currentInToroidR(r, r0, y, turns, I) * cos(theta)

END FUNCTION currdens3
