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

  f = 128.0d0 * ( (yq**2 - 1.0d0/4.0d0)*(zq**2 - 1.0d0/4.0d0) + &
      (xq**2 - 1.0d0/4.0d0)*(zq**2 - 1.0d0/4.0d0) + &  
      (xq**2 - 1.0d0/4.0d0)*(yq**2 - 1.0d0/4.0d0) )

END FUNCTION SourceFun
