SUBROUTINE HeatCondUniso( Model, n, t, f )
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t
  REAL(KIND=dp) :: f(:,:)  ! this size needs to be consistent with the sif file!

!  REAL(KIND=dp) :: x,y,z
!  TYPE(Mesh_t), POINTER :: Mesh
    
!  Mesh => GetMesh()
!  x = Mesh % Nodes % x(n)   
!  y = Mesh % Nodes % y(n)   
!  z = Mesh % Nodes % z(n)   

  PRINT *,'size',t,SIZE(f)

  f(1,1) = 1.0_dp + t / 5.0_dp
  f(1,2) = 0.0_dp 
  f(2,1) = 0.0_dp
  f(2,2) = 20.0_dp - 2.0_dp * t
  
END SUBROUTINE HeatCondUniso
