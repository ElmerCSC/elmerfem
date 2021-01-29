FUNCTION InletVelo(Model, n, t) RESULT(v)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, v

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: xq, yq, zq

  Mesh => Model % Mesh

  xq = Mesh % Nodes % x(n)
  yq = Mesh % Nodes % y(n)
  zq = Mesh % Nodes % z(n)

  v = 100.0*(1.0e-4-xq**2-yq**2)

END FUNCTION InletVelo
