module loop_source
  use DefUtils
  IMPLICIT none

  real(kind=dp) :: r_loop = 10_dp-7
  real(kind=dp) :: I_loop = 1._dp
  real(kind=dp) :: OX = 0.1_dp
  real(kind=dp) :: OY = 0.1_dp
  real(kind=dp) :: OZ = 0.1_dp
  real(kind=dp) :: reps = 1.d-7

  contains
    function a_phi(angle, r, x, z, I)
      real(kind=dp), intent(in) :: angle, r, x, z, I
      real(kind=dp) :: a_phi
      real(kind=dp), parameter :: mu_0 = 4.d-7*pi


      a_phi=mu_0/4._dp*(r**2._dp*I*x)/(r**2._dp+x**2._dp+z**2._dp)**(3._dp/2._dp) * &
                     (1._dp+(15._dp*r**2* x**2)*(8._dp*(r**2._dp+x**2._dp+z**2._dp)**2._dp))
    end function a_phi
end module

FUNCTION SourceFunX(Model, n, t) RESULT(f)
  USE DefUtils
  USE loop_source
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: xq, yq, zq, angle

  Mesh => Model % Mesh

  xq = Mesh % Nodes % x(n)-OX
  yq = Mesh % Nodes % y(n)-OY
  zq = Mesh % Nodes % z(n)-OZ

  angle = atan2(yq, xq)
  f = cos(angle) * a_phi(angle, r_loop, xq, zq, I_loop)
END FUNCTION SourceFunX

FUNCTION SourceFunY(Model, n, t) RESULT(f)
  USE DefUtils
  USE loop_source
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: xq, yq, zq, angle

  Mesh => Model % Mesh

  xq = Mesh % Nodes % x(n)-OX
  yq = Mesh % Nodes % y(n)-OY
  zq = Mesh % Nodes % z(n)-OZ

  angle = atan2(yq, xq)
  f = sin(angle) * a_phi(angle, r_loop, xq, zq, I_loop)
END FUNCTION SourceFunY
