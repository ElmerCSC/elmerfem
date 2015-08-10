FUNCTION getWindingSigmaConstantTemperature( model, n, T ) RESULT(eff_sigma)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, WindingT, T, rho, rho0, rho1, ft, it, ff, c
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: FOUND, FOUNDFT, FOUNDIT, FOUNDC

  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getWindingSigma', 'Body Parameters not found!')
  
  ft = GetConstReal(BodyParams, 'Foil Layer Thickness', FoundFT)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Foil Layer Thickness not found in Body section')

  it = GetConstReal(BodyParams, 'Insulator Layer Thickness', FoundIT)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Insulator Layer Thickness not found in Body section')
  
  c = GetConstReal(BodyParams, 'Insulator Portion', FoundC)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Insulator Layer Thickness not found in Body section')
  
  IF (((.NOT. FOUNDFT).OR.(.NOT. FOUNDIT)).AND.(.NOT. FOUNDC)) &
    CALL Fatal('getWindingSigma', 'Insulator Values not found in Body section')
  
  WindingT = GetConstReal(BodyParams, 'Winding Temperature', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Winding Temperature not found in Body section')
  
  rho1 = GetConstReal(Material, 'rho1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho1 not found in Material section')
  
  rho0 = GetConstReal(Material, 'rho0', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho0 not found in Material section')
  
  rho = (rho1 * (273.15_dp+WindingT) + rho0)*1.d-9
  alsigma = 1._dp/rho
  
  ! filling factor
  IF (.NOT. FoundC) THEN
    ff = ft/(ft+it)
  ELSE
    ff = 1-c
  END IF
  ! corrected conductivity
  eff_sigma = ff * alsigma
  
END FUNCTION getWindingSigmaConstantTemperature


FUNCTION getWindingSigma( model, n, T ) RESULT(eff_sigma)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, T, rho, rho0, rho1, ft, it, ff, c
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: FOUND, FOUNDFT, FOUNDIT, FOUNDC

  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getWindingSigma', 'Body Parameters not found!')
  
  ft = GetConstReal(BodyParams, 'Foil Layer Thickness', FoundFT)

  it = GetConstReal(BodyParams, 'Insulator Layer Thickness', FoundIT)

  c = GetConstReal(BodyParams, 'Insulator Portion', FoundC)
  
  IF (((.NOT. FOUNDFT).OR.(.NOT. FOUNDIT)).AND.(.NOT. FOUNDC)) &
    CALL Fatal('getWindingSigma', 'All needed values not found in Body section')
    
  rho1 = GetConstReal(Material, 'rho1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho1 not found in Material section')
  
  rho0 = GetConstReal(Material, 'rho0', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'rho0 not found in Material section')
  
  rho = (rho1 * (273.15_dp+T) + rho0)*1.d-9
  alsigma = 1._dp/rho
  
  ! filling factor
  IF (.NOT. FoundC) THEN
    ff = ft/(ft+it)
  ELSE
    ff = 1-c
  END IF
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
  REAL(KIND=dp) ::  Conductivity(:)
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  REAL(KIND=dp) :: ft, mc, it, ic, st, Ksi, Ktot

  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND, FOUNDFT, FOUNDST
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingk', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getWindingk', 'Body Parameters not found!')
  
  ft = GetConstReal(BodyParams, 'Foil Layer Thickness', FoundFT)
  
  st = GetConstReal(BodyParams, 'Strand Thickness', FoundST)
  
  mc = GetConstReal(Material, 'Material Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Material Heat Condictivity not found in Material section')
  
  it = GetConstReal(BodyParams, 'Insulator Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Thickness not found in Body section')

  ic = GetConstReal(BodyParams, 'Insulator Layer Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Heat Conductivity not found in Body section')

  IF ((.NOT. FOUNDFT) .AND. (.NOT. FOUNDST)) CALL Fatal('getWindingk', 'Material Thickness not found in Body section')
  
  IF (.NOT. FOUNDST) THEN
    Conductivity(1) = ((ft + it) * mc * ic) / (it * mc + ft * ic)
    Conductivity(2) = (ft * mc + it * ic) / (ft + it)
  ELSE
    Ksi = (st*mc)+(it*ic)/(st+it) !((st+it) * mc * ic)/((ic*st)+(mc+it))
    Ktot = ((st+it)*Ksi*ic)/((it*Ksi)+(st*ic))     !(Ksi*st + ic*it)/(st+it)
    Conductivity(1) = Ktot
    Conductivity(2) = Ktot
  END IF
  !WRITE(Message,*)  Conductivity(1,1), Conductivity(1,2)
  !CALL Info('getWindingk', Message, Level = 5)
END SUBROUTINE getWindingk
