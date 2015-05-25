FUNCTION gravity( model, n, temperature ) RESULT(Gforce)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, Gforce

  Gforce = -9.82*(1.205*(300/(273+temperature)))

END FUNCTION gravity

FUNCTION heating( model, n, JH ) RESULT(JHmass)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: JH, JHmass, density

  density = 2700._dp
  
  JHmass = JH/density

END FUNCTION heating

FUNCTION getAirViscosity( model, n, x ) RESULT(viscosity)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: viscosity, x

  if ( x > 0.40_dp ) then
    viscosity = 1e-3
  else
    viscosity = 1e-5
  end if

END FUNCTION getAirViscosity


FUNCTION getTcondition( model, n, args ) RESULT(Condition)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: args(2), Condition, x, y

  x=args(1)
  y=args(2)

  if ( x > 0.450_dp ) then
    Condition = 1d0
  else
    Condition = -1d0
  end if

END FUNCTION getTcondition

FUNCTION getAlWindingSigma( model, n, T ) RESULT(eff_sigma)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, T, rho, rho0, rho1, ft, it, ff
  
  TYPE(Valuelist_t), POINTER :: Material
  LOGICAL :: FOUND
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getAlWindingSigma', 'Material not found')
 
  ft = GetConstReal(Material, 'Foil Layer Thickness', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getAlWindingSigma', 'Foil Layer Thickness not found in Material section')

  it = GetConstReal(Material, 'Insulator Layer Thickness', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getAlWindingSigma', 'Insulator Layer Thickness not found in Material section')

  rho1 = 0.113518_dp
  rho0 = -0.678964_dp
  
  rho = (rho1 * (273.15_dp+T) + rho0)*1.d-9
  alsigma = 1._dp/rho
  
  ! filling factor
  ff = ft/(ft+it)
  ! corrected conductivity
  eff_sigma = ff * alsigma
  
END FUNCTION getAlWindingSigma

