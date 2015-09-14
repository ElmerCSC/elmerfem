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
  TYPE(Element_t), POINTER :: Element
  TYPE(Valuelist_t), POINTER :: Material
  Logical :: Found
  
  Element => GetCurrentElement()
  Material => GetMaterial(Element, Found)
  IF (.not. Found) CALL Fatal('heating', 'Material not found')
  
  density = GetConstReal(Material, 'density', Found)
  IF (.NOT. FOUND) CALL Fatal('heating', 'density not found in Material section')
  
  JHmass = JH/density

END FUNCTION heating

FUNCTION getAirViscosity( model, n, x ) RESULT(viscosity)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: viscosity, x, viscosityLowerLimit
  
  TYPE(Valuelist_t), POINTER :: Material
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getAirViscosity', 'Material not found')
  
  viscosityLowerLimit = GetConstReal(Material, 'Viscosity Lower Limit', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Viscosity Lower Limit not found in Material section')
  
  if (x > viscosityLowerLimit) then !( x > 0.40_dp ) then
    viscosity = 1e-3
  else
    viscosity = 1e-5
  end if

END FUNCTION getAirViscosity

FUNCTION getLinearMethod( model, n, x ) RESULT(temperature)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, x, x1, x2, t1, t2, k, b
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getLinearMethod', 'Body Parameters not found!')
  
  x1 = GetConstReal(BodyParams, 'Temperature Limit Position 1', Found)
  IF (.NOT. FOUND) CALL Fatal('getLinearMethod', 'Temperature Limit Position 1 not found in Body section!')
  
  x2 = GetConstReal(BodyParams, 'Temperature Limit Position 2', Found)
  IF (.NOT. FOUND) CALL Fatal('getLinearMethod', 'Temperature Limit Position 2 not found in Body section!')
  
  t1 = GetConstReal(BodyParams, 'Initial Air Temperature', Found)
  IF (.NOT. FOUND) CALL Fatal('getLinearMethod', 'Initial Air Temperature not found in Body section!')
  
  t2 = GetConstReal(BodyParams, 'Temperature', Found)
  IF (.NOT. FOUND) CALL Fatal('getLinearMethod', 'Temperature not found in Body section!')
  
  IF (x < x1) THEN
    temperature = t1
  ELSEIF (x1 < x .AND. x < x2) THEN
    k = (t2-t1)/(x2-x1)
    b = t1 - k*x1
    temperature = k*x + b
  ELSE
    temperature = t2
  END IF
  
END FUNCTION getLinearMethod


FUNCTION getTcondition( model, n, args ) RESULT(Condition)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(Valuelist_t), POINTER :: BF
  INTEGER :: n
  REAL(KIND=dp) :: args(2), Condition, x, y, limit
  LOGICAL :: FOUND

  x=args(1)
  y=args(2)
  
  BF => GetBodyForce() 
  IF (.not. ASSOCIATED(BF)) CALL Fatal('getTcondition', 'Body Force not found')
  
  limit = GetConstReal(BF, 'Temperature Condition Limit Position', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getTcondition', 'Temperature Condition Limit Position not found in Body Force section')
  
  if ( x > limit ) then
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

