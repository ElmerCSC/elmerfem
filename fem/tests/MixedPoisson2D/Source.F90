FUNCTION SourceFun(Model, n, t) RESULT(f)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: xq, yq, zq

  Mesh => Model % Mesh

  xq = Mesh % Nodes % x(n)
  yq = Mesh % Nodes % y(n)
  zq = Mesh % Nodes % z(n)

  !f = 32.0d0 * ( yq-yq**2 + xq-xq**2 )
  !f = 8.0d0 * xq
  f = 0.0d0

END FUNCTION SourceFun
