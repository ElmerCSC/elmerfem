FUNCTION getWindingSigma( model, n, T ) RESULT(eff_sigma)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, T, rho, rho0, rho1, ft, it, ff
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: FOUND

  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getWindingSigma', 'Body Parameters not found!')
  
  ft = GetConstReal(BodyParams, 'Foil Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Foil Layer Thickness not found in Body section')

  it = GetConstReal(BodyParams, 'Insulator Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Insulator Layer Thickness not found in Body section')
  
  rho1 = GetConstReal(Material, 'rho1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho1 not found in Material section')
  
  rho0 = GetConstReal(Material, 'rho0', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho0 not found in Material section')
  
  rho = (rho1 * (273.15_dp+T) + rho0)*1.d-9
  alsigma = 1._dp/rho
  
  ! filling factor
  ff = ft/(ft+it)
  ! corrected conductivity
  eff_sigma = ff * alsigma
  
END FUNCTION getWindingSigma

SUBROUTINE getWindingk( model, n, dummyArgument,Conductivity )
  ! modules needed
  USE DefUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp),POINTER ::  Conductivity(:,:)
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  REAL(KIND=dp) :: ft, fc, it, ic

  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingk', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getWindingk', 'Body Parameters not found!')
  
  ft = GetConstReal(BodyParams, 'Foil Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Foil Layer Thickness not found in Body section')

  fc = GetConstReal(Material, 'Foil Layer Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Foil Layer Heat Condictivity not found in Material section')
 
  it = GetConstReal(BodyParams, 'Insulator Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Thickness not found in Body section')

  ic = GetConstReal(BodyParams, 'Insulator Layer Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Heat Conductivity not found in Body section')
  
  Conductivity(1,1) = ((ft + it) * fc * ic) / (it * fc + ft * ic)
  Conductivity(1,2) = (ft * fc + it * ic) / (ft + it)
  !WRITE(Message,*)  Conductivity(1,1), Conductivity(1,2)
  !CALL Info('getWindingk', Message, Level = 5)
 
END SUBROUTINE getWindingk
