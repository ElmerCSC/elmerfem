FUNCTION SphereDistance( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f

  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp) :: x,y,z
  TYPE(Mesh_t), POINTER :: Mesh
  
  SAVE Visited, Mesh
  
  IF( .NOT. Visited ) THEN
    Mesh => GetMesh()
    Visited = .TRUE.
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  f = SQRT(x*x+y*y+z*z)
  
END FUNCTION SphereDistance

