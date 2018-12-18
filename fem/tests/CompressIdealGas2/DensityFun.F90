! Density (\rho) and density sensitivity to pressure (d\rho/dp) for ideal gases
! Formulas defined around reference pressure. 

FUNCTION GasDensity( Model, n, dpres ) RESULT( rho)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: dpres, rho

  REAL(KIND=dp) :: Pres0, Rho0, Ratio, InvRatio
  LOGICAL :: Visited  = .FALSE.
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: pres, comp

  SAVE Visited, Pres0, Rho0, Ratio, InvRatio


  IF( .NOT. Visited ) THEN
    Mesh => GetMesh()

    Material => GetMaterial( Mesh % Elements(1) )

    Ratio = ListGetCReal( Material,'Specific Heat Ratio')
    InvRatio = 1.0_dp / Ratio
    Rho0 = ListGetCReal( Material,'Equilibrium Density')
    Pres0 = ListGetCReal( Material,'Equilibrium Pressure')

    Visited = .TRUE.
  END IF

  rho = Rho0 * ( dpres/pres0 + 1.0_dp )**InvRatio     

END FUNCTION GasDensity


FUNCTION GasDrhodp( Model, n, dpres ) RESULT( drhodp )
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: dpres, drhodp

  REAL(KIND=dp) :: Pres0, Rho0, Ratio, InvRatio
  LOGICAL :: Visited  = .FALSE.
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: coeff, pres, rho, comp

  SAVE Visited, Pres0, Rho0, Ratio, InvRatio, Coeff


  IF( .NOT. Visited ) THEN
    Material => GetMaterial()
    Mesh => GetMesh()
    Ratio = ListGetCReal( Material,'Specific Heat Ratio')
    InvRatio = 1.0_dp / Ratio
    Rho0 = ListGetCReal( Material,'Equilibrium Density')
    Pres0 = ListGetCReal( Material,'Equilibrium Pressure')

    Coeff = Rho0 / ( Ratio * Pres0**InvRatio )

    Visited = .TRUE.
  END IF

  drhodp = Coeff * (dpres+pres0)**(InvRatio-1)

END FUNCTION GasDrhodp


