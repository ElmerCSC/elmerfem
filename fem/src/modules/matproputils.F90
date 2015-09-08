MODULE trapezoid
USE DefUtils
CONTAINS
  SUBROUTINE trapezoid_for_data(data_x, size_x, data_y,size_y,  a,b,integral,average)
  !-----------------------------------------------------------------------------------------
  !*****************************************************************************************
  !
  ! Calculate the integral for data point with start and end point for discrete data points.
  ! Return value of integral and average value, i.e. integral/(end point - start point).
  !
  ! ARGUMENTS: 
  ! data_x   - x data for integration
  ! data_y   - y data for integration
  ! a   - Lower limit of integration
  ! b    - Upper limit of integration
  ! size_x - size of the data_x
  ! size_y - size of the data_y
  !
  ! RETURN:
  ! integral - Result of integration
  ! average  - Average value of the function within interwall
  !*****************************************************************************************

  IMPLICIT NONE
  REAL(KIND = dp) a, b, integral,average
  REAL(KIND = dp) h
  REAL(KIND = dp) point_x1,point_y1, point_x2,point_y2
  REAL(KIND = dp), POINTER :: fit_x(:), fit_y(:) 
  REAL(KIND = dp) :: data_x(size_x), data_y(size_y) 
  INTEGER nint, size_x, size_y, loc,loc_start,loc_end
  INTEGER i
  LOGICAL loc_extra

  nint = 0

  !PRINT *, "limits", a,b

  ! Find start point 
  CALL find_point_in_array(data_x,size_x,data_y,size_y,a,point_x1,point_y1,loc_start, loc_extra)
  IF (loc_extra) THEN
    nint = nint 
  ENDIF
  
  ! Find end point
  CALL find_point_in_array(data_x,size_x,data_y,size_y,b,point_x2,point_y2,loc_end, loc_extra)  
  nint =  nint +  loc_end-loc_start +1   
  IF (loc_extra) THEN
    nint = nint + 1
  ENDIF   
  
  integral = compute_integral(nint, data_x, data_y, loc_start, point_x1, point_y1, point_x2, point_y2)
  
  average = integral/(point_x2-point_x1)


  END SUBROUTINE trapezoid_for_data
  
  SUBROUTINE find_point_in_array(data_x,size_x,data_y,size_y,test_point, point_x,point_y,loc, loc_extra)
  !-----------------------------------------------------------------------------------------
  !*****************************************************************************************
  !
  !  Find the given point in array with x_data points.  
  !  Return the index of the location in array.
  !  If such point is not found from the list, interpolate linearly the location.
  !  If the given point is smaller (larger) than smallest (largest) data point, 
  !    use the smallest (largest) data y point. For point between datapoints, use simple 
  !    linear interpolation.
  !    
  !
  !  EXAMPLE DATA:
  !            ***     *-------
  !           *  **---*
  !          * 
  !  ----***
  !
  !  
  !  * =  datapoints
  !  - =  interpolated data
  !
  !  ARGUMENTS:
  !
  !  REAL(KIND = dp) :: data_x 
  !  REAL(KIND = dp) :: data_y
  !  INTEGER :: size_x
  !  INTEGER :: size_y
  !  REAL(KIND = dp) :: test_point
  !     INPUT: the x-coordinate to be tested. 
  ! 
  ! RETURN:
  !
  !  REAL(KIND = dp) :: point_x - the starting/ending x point for integration
  !  REAL(KIND = dp) :: point_y - the starting/ending y point for integration   
  ! INTEGER :: loc - the index of the point. If not found, use 0 for point smaller than data min, 
  !            otherwise use the     nearest point smaller than the test_point 
  ! LOGICAL :: loc_extra - .TRUE. if extra points are created
  !*****************************************************************************************

  INTEGER size_x, size_y, loc
  LOGICAL loc_extra
  REAL(KIND = dp) data_x(size_x), data_y(size_y)
  REAL(KIND = dp) test_point, point_x, point_y
  REAL(KIND = dp) line_a,line_b

  loc_extra = .TRUE.

  DO i = 1,size_x,1
    IF (test_point == data_x(i)) THEN
       point_x =  test_point     
       loc = i
       loc_extra = .FALSE.
       point_y = data_y(loc)
       EXIT
    ENDIF
  ENDDO

  !find the nearest point from array if the excact location is not found

  IF (loc_extra) THEN
     IF (test_point < data_x(1)) THEN
        point_x = test_point
        point_y = data_y(1)
        loc = 0
        loc_extra = .TRUE.
     ELSEIF (test_point > data_x(size_x)) THEN
        point_x = test_point
        point_y = data_y(size_y)
        loc = size_x
        loc_extra = .TRUE.
     ELSE
        DO i = 2,size_x,1
           IF (test_point > data_x(i-1) .AND. test_point < data_x(i)) THEN
              loc = i-1
              loc_extra = .TRUE.
              IF (data_x(i) /= data_x(i-1)) THEN
                 line_a = (data_y(i)-data_y(i-1))/(data_x(i)-data_x(i-1))
              ELSE
                 line_a = 0
              ENDIF
              line_b = data_y(i-1) - line_a*data_x(i-1)
              point_x =  test_point
              point_y = line_a * point_x + line_b
              EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDIF
  END  SUBROUTINE find_point_in_array
  
  FUNCTION compute_integral(nint, data_x, data_y, loc_start, point_x1, point_y1, point_x2, point_y2) RESULT (integral)
  !-----------------------------------------------------------------------------------------
  !*****************************************************************************************
  !
  !  Compute the integral according to data array. 
  !    
  !
  !  ARGUMENTS:
  !
  !  INTEGER :: nint - number of datapoints
  !  REAL(KIND = dp) :: data_x 
  !  REAL(KIND = dp) :: data_y
  !  INTEGER :: loc_start - start index for integration in data
  !  REAL(KIND = dp) :: point_x1 - starting x-point for integration 
  !  REAL(KIND = dp) :: point_x2 - ending x-point for integration 
  !  REAL(KIND = dp) :: point_y1 - starting y-point for integration 
  !  REAL(KIND = dp) :: point_y2 - ending y-point for integration 
  ! 
  ! RETURN:
  !
  !  REAL(KIND = dp) :: integral - teh value of the integral
  !
  !*****************************************************************************************
    INTEGER :: nint, i, loc_start
    REAL(KIND = dp) :: point_x1, point_y1, point_x2, point_y2
    REAL(KIND = dp) :: data_x(:), data_y(:)
    REAL(KIND = dp) :: fit_x(nint), fit_y(nint)
    REAL(KIND = dp) :: integral, h
    
    fit_x(1) = point_x1
    fit_y(1) = point_y1
    
    DO i = 2,nint-1,1
       fit_x(i) = data_x(loc_start+i-1)
       fit_y(i) = data_y(loc_start+i-1)
    END DO
    
    fit_x(nint) = point_x2
    fit_y(nint) = point_y2

    integral = 0.d0
    DO i = 1, nint-1,1
       h = abs(fit_x(i+1) - fit_x(i))
       integral = integral + 0.5*h*(fit_y(i)+fit_y(i+1))
    END DO

  END FUNCTION compute_integral

END MODULE trapezoid

MODULE CompUtils

CONTAINS
  !------------------------------------------------------------------------------
  FUNCTION GetComponentParams(Element) RESULT (ComponentParams)
  !------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
      
    INTEGER :: i
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: BodyParams, ComponentParams
    LOGICAL :: Found
      
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('GetCompParams', 'Component parameters not found')
      
    i = GetInteger(BodyParams, 'Component', Found)

    IF (.NOT. Found) CALL Fatal ('GetComponentParams', 'Body not associated to any Component!')

    ComponentParams => CurrentModel % Components(i) % Values
        
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('CircuitsAndDynamicsHarmonic', &
                                                           'Component parameters not found!')
      
  !------------------------------------------------------------------------------
  END FUNCTION GetComponentParams
  !------------------------------------------------------------------------------
END MODULE CompUtils


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !  
!             FUNCTIONS CALLABLE FROM ELMER SIF                             !
!                                                                           ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION getExpCoeffEpikotem( model, n, T ) RESULT(ExpCoeff_average)
!===================================================================
! Calculate the thermal expansion coefficient for Epikotem 
!===================================================================
USE DefUtils
USE trapezoid
IMPLICIT None
TYPE(Model_t) :: model
INTEGER :: n



REAL(KIND = dp) T,Tref,a,b,f, x,integral,ExpCoeff_average
REAL(KIND = dp), DIMENSION(63) :: data_y, data_x
INTEGER i, size_x,size_y
TYPE(Valuelist_t), POINTER :: Material
LOGICAL :: FOUND

Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')

Tref = GetConstReal(Material, 'Reference Temperature', FOUND)
IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Reference Temperature not found in Material section')

! Epikotem data
data_y = (/29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, &
           29.969, 29.969, 30.3096, 30.6502, 30.6502, 30.9907, 30.9907, 31.3313, &
           31.6718, 31.6718, 32.0124, 32.3529, 32.6935, 33.0341, 33.3746, 33.7152, &
           34.0557, 34.3963, 35.0774, 37.1207, 38.8235, 40.8669, 42.9102, 44.9536, &
           47.678, 50.743, 53.808, 57.2136, 60.2786, 63.3437, 66.4087, 69.8142, 72.8793, &
           75.9443, 79.0093, 81.7337, 84.1176, 86.8421, 89.226, 91.6099, 93.9938, 95.6966, &
           96.7183, 98.0805, 99.1022, 100.124, 100.124, 100.124, 100.124, 100.124, 100.124, &
           100.124, 100.124, 100.124/)
data_x = (/0.0, 0.236686, 1.89349, 4.26036, 6.62722, 8.99408, 11.3609, 13.7278, &
           16.0947, 18.4615, 20.8284, 22.9586, 25.3254, 27.6923, 30.0592, 32.426, &
           34.7929, 36.9231, 39.2899, 41.6568, 44.0237, 46.3905, 48.7574, 50.8876, &
           53.2544, 55.6213, 57.9882, 60.1183, 62.0118, 63.9053, 66.0355, 67.929, &
           69.8225, 70.7692, 71.716, 72.6627, 73.6095, 74.5562, 75.503, 76.213, &
           77.1598, 78.1065, 79.0533, 80.0, 81.4201, 83.0769, 84.7337, 86.1538, &
           87.8107, 89.4675, 91.3609, 93.7278, 95.858, 98.2249, 100.355, 102.722, &
           105.089, 107.456, 109.822, 112.189, 114.556, 116.923, 119.29/)

size_x = SIZE(data_x)
size_y = SIZE(data_y)

IF (T >= Tref) THEN
  a = Tref
  b = T
ELSE
  a = T
  b = Tref
ENDIF

CALL trapezoid_for_data(data_x, size_x, data_y,size_y, a,b,integral,ExpCoeff_average)

ExpCoeff_average=ExpCoeff_average*1.d-6

END FUNCTION getExpCoeffEpikotem

FUNCTION getExpCoeffEpikotemFgnetComposite( model, n, T ) RESULT(ExpCoeff_average)
!===================================================================
! Calculate the thermal expansion coefficient for Epikotem 
!===================================================================
USE DefUtils
USE trapezoid
IMPLICIT None
TYPE(Model_t) :: model
INTEGER :: n



REAL(KIND = dp) T,Tref,a,b,f, x,integral,ExpCoeff_average, ExpCoeff_fgnet, Vf_ek, Vf_fgnet
REAL(KIND = dp), DIMENSION(63) :: data_y, data_x
INTEGER i, size_x,size_y
TYPE(Valuelist_t), POINTER :: Material
LOGICAL :: FOUND

Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')

Tref = GetConstReal(Material, 'Reference Temperature', FOUND)
IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Reference Temperature not found in Material section')

Vf_fgnet = GetConstReal(Material, 'Volume Fraction Material 2', FOUND)
IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Volume Fraction Material 2 not found in Material section')

Vf_ek = 1.0 - Vf_fgnet

! Epikotem data
data_y = (/29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, 29.969, &
           29.969, 29.969, 30.3096, 30.6502, 30.6502, 30.9907, 30.9907, 31.3313, &
           31.6718, 31.6718, 32.0124, 32.3529, 32.6935, 33.0341, 33.3746, 33.7152, &
           34.0557, 34.3963, 35.0774, 37.1207, 38.8235, 40.8669, 42.9102, 44.9536, &
           47.678, 50.743, 53.808, 57.2136, 60.2786, 63.3437, 66.4087, 69.8142, 72.8793, &
           75.9443, 79.0093, 81.7337, 84.1176, 86.8421, 89.226, 91.6099, 93.9938, 95.6966, &
           96.7183, 98.0805, 99.1022, 100.124, 100.124, 100.124, 100.124, 100.124, 100.124, &
           100.124, 100.124, 100.124/)
data_x = (/0.0, 0.236686, 1.89349, 4.26036, 6.62722, 8.99408, 11.3609, 13.7278, &
           16.0947, 18.4615, 20.8284, 22.9586, 25.3254, 27.6923, 30.0592, 32.426, &
           34.7929, 36.9231, 39.2899, 41.6568, 44.0237, 46.3905, 48.7574, 50.8876, &
           53.2544, 55.6213, 57.9882, 60.1183, 62.0118, 63.9053, 66.0355, 67.929, &
           69.8225, 70.7692, 71.716, 72.6627, 73.6095, 74.5562, 75.503, 76.213, &
           77.1598, 78.1065, 79.0533, 80.0, 81.4201, 83.0769, 84.7337, 86.1538, &
           87.8107, 89.4675, 91.3609, 93.7278, 95.858, 98.2249, 100.355, 102.722, &
           105.089, 107.456, 109.822, 112.189, 114.556, 116.923, 119.29/)

size_x = SIZE(data_x)
size_y = SIZE(data_y)

IF (T >= Tref) THEN
  a = Tref
  b = T
ELSE
  a = T
  b = Tref
ENDIF

CALL trapezoid_for_data(data_x, size_x, data_y,size_y, a,b,integral,ExpCoeff_average)

ExpCoeff_fgnet=5.3
ExpCoeff_average=ExpCoeff_average*Vf_ek+ExpCoeff_fgnet*Vf_fgnet

ExpCoeff_average=ExpCoeff_average*1.d-6

END FUNCTION getExpCoeffEpikotemFgnetComposite


FUNCTION getWindingSigmaConstantTemperature( model, n, T ) RESULT(eff_sigma)
  USE DefUtils
  USE CompUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, WindingT, T, rho, rho0, rho1, ft, it, ff, c
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, ComponentParams
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: FOUND, FOUNDFT, FOUNDIT, FOUNDC

  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigmaConstant', 'Material not found')
  
  Element => GetCurrentElement()
  !CircuitVariableBody => Model % Bodies (Element % BodyId)
  !BodyParams => CircuitVariableBody % Values
  ComponentParams => GetComponentParams(Element)
  IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('getWindingSigmaConstant', 'Body Parameters not found!')
  
  ft = GetConstReal(ComponentParams, 'Foil Layer Thickness', FoundFT)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Foil Layer Thickness not found in Body section')

  it = GetConstReal(ComponentParams, 'Insulator Layer Thickness', FoundIT)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Insulator Layer Thickness not found in Body section')
  
  c = GetConstReal(ComponentParams, 'Insulator Portion', FoundC)
!  IF (.NOT. FOUND) CALL Fatal('getWindingSigma', 'Insulator Layer Thickness not found in Body section')
  
  IF (((.NOT. FOUNDFT).OR.(.NOT. FOUNDIT)).AND.(.NOT. FOUNDC)) &
    CALL Fatal('getWindingSigmaConstant', 'Insulator Values not found in Component section')
  
  WindingT = GetConstReal(ComponentParams, 'Winding Temperature', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigmaConstant', 'Winding Temperature not found in Component section')
  
  rho1 = GetConstReal(Material, 'rho1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigmaConstant', 'rho1 not found in Material section')
  
  rho0 = GetConstReal(Material, 'rho0', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getWindingSigmaConstant', 'rho0 not found in Material section')
  
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
  USE CompUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: eff_sigma, alsigma, T, rho, rho0, rho1, ft, it, ff, c
  
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  TYPE(Valuelist_t), POINTER :: Material, ComponentParams
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: FOUND, FOUNDFT, FOUNDIT, FOUNDC

  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingSigma', 'Material not found')
  
  Element => GetCurrentElement()
  ComponentParams => GetComponentParams(Element)
  IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('getWindingSigma', 'Component Parameters not found!')
  
  ft = GetConstReal(ComponentParams, 'Foil Layer Thickness', FoundFT)

  it = GetConstReal(ComponentParams, 'Insulator Layer Thickness', FoundIT)

  c = GetConstReal(ComponentParams, 'Insulator Portion', FoundC)
  
  IF (((.NOT. FOUNDFT).OR.(.NOT. FOUNDIT)).AND.(.NOT. FOUNDC)) &
    CALL Fatal('getWindingSigma', 'All needed values not found in Component section')
    
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE getElasticModulusComposite( model, n, dummyArgument,ElasticModulus )
  ! modules needed
  USE DefUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  ElasticModulus(:,:)
  REAL(KIND=dp) :: Em,Ef,Vf,E_perpendicular, E_parallel,G12,Gm,Gf
  REAL(KIND=dp) :: Ex, Ey, Ez, nu,nu_xy,nu_xz,nu_yx,nu_yz,nu_zx,nu_zy, D,GPA
  REAL(KIND=dp) :: nu_parallel, nu_perpendicular, nu_m, nu_f

  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND



  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingk', 'Material not found') 
  
  Em = GetConstReal(Material, 'Youngs Modulus Material 1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getElasticModulus', 'Youngs Modulus for Material 1 not found')
  Ef = GetConstReal(Material, 'Youngs Modulus Material 2', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getElasticModulus', 'Youngs Modulus for Material 2 not found')
  Vf = GetConstReal(Material, 'Volume Fraction Material 2', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getElasticModulus', 'Volume Fraction for Material 2 not found')

  nu_m = GetConstReal(Material, 'Poisson ratio Material 1', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getElasticModulus', 'Poisson ratio for Material 1 not found')
  nu_f = GetConstReal(Material, 'Poisson ratio Material 2', FOUND)
  IF (.NOT. FOUND) CALL Fatal('getElasticModulus', 'Poisson ratio for Material 2 not found')

  E_parallel = Em*(1-Vf) + Ef*Vf
  E_perpendicular = 1/(((1-Vf)/Em)+(Vf/Ef))

  nu_parallel = nu_m*(1-Vf) + nu_f*Vf
  !nu_perpendicular = 1/(((1-Vf)/nu_m)+(Vf/nu_f))
  
  nu_perpendicular = (E_perpendicular/E_parallel)*nu_parallel

  Gm = Em/(2*(1+nu_m))
  Gf = Em/(2*(1+nu_f))

  G12 = Gm*Gf/((1-Vf)*Gf+Vf*Gm)

  Ex = 1/(((1-Vf)/Em)+(Vf/Ef)) ! perpendicular
  Ey = Em*(1-Vf) + Ef*Vf ! parallel, direction of fibres
  Ez = Ex ! perpendicular, same as in x-direction

  !Ex nu_yx = Ey nu_xy

  nu_xy = nu_perpendicular 
  nu_xz = nu_perpendicular
  nu_yx = nu_parallel
  nu_yz = nu_parallel
  nu_zx = nu_perpendicular  
  nu_zy = nu_perpendicular 
  
  D = 1/(1-nu_xy*nu_yx-nu_yz*nu_zy-nu_zx*nu_xz-2*nu_xy*nu_yz*nu_zx)
 
  ElasticModulus(1,1) = (1-nu_yz*nu_zy)*D*Ex  
  ElasticModulus(1,2) = (nu_yx+nu_zx*nu_yz)*D*Ex
  ElasticModulus(1,3) = (nu_zx+nu_yx*nu_zy)*D*Ex
  ElasticModulus(1,4) = 0
 
  ElasticModulus(2,2) = (1-nu_zx*nu_xz)*D*Ey  
  ElasticModulus(2,1) = (nu_xy+nu_xz*nu_zy)*D*Ey
  ElasticModulus(2,3) = (nu_zy+nu_zx*nu_xy)*D*Ey
  ElasticModulus(2,4) = 0

  ElasticModulus(3,3) = (1-nu_xy*nu_yx)*D*Ez  
  ElasticModulus(3,1) = (nu_xz+nu_xy*nu_yz)*D*Ez
  ElasticModulus(3,2) = (nu_yz+nu_xz*nu_yx)*D*Ez
  ElasticModulus(3,4) = 0

  ElasticModulus(4,1) = 0 
  ElasticModulus(4,2) = 0 
  ElasticModulus(4,3) = 0 
  ElasticModulus(4,4) = G12 ! Ex/(2*(1+nu_yx)) 
 
END SUBROUTINE getElasticModulusComposite


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE getWindingk( model, n, dummyArgument,Conductivity )
  ! modules needed
  USE DefUtils
  USE CompUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  Conductivity(:)
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  REAL(KIND=dp) :: ft, mc, it, ic, st, Ksi, Ktot

  TYPE(Valuelist_t), POINTER :: Material, ComponentParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND, FOUNDFT, FOUNDST
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getWindingk', 'Material not found')
  
  Element => GetCurrentElement()
  ComponentParams => GetComponentParams(Element)
  IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('getWindingk', 'Component Parameters not found!')
  
  ft = GetConstReal(ComponentParams, 'Foil Layer Thickness', FoundFT)
  
  st = GetConstReal(ComponentParams, 'Strand Thickness', FoundST)
  
  mc = GetConstReal(Material, 'Material Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Material Heat Condictivity not found in Material section')
  
  it = GetConstReal(ComponentParams, 'Insulator Layer Thickness', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Thickness not found in Component section')

  ic = GetConstReal(ComponentParams, 'Insulator Layer Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getWindingk', 'Insulator Layer Heat Conductivity not found in Component section')

  IF ((.NOT. FOUNDFT) .AND. (.NOT. FOUNDST)) CALL Fatal('getWindingk', 'Material Thickness not found in Body section')
  
  IF (.NOT. FOUNDST) THEN
    Conductivity(1) = ((ft + it) * mc * ic) / (it * mc + ft * ic)
    Conductivity(2) = (ft * mc + it * ic) / (ft + it)
  ELSE
    Ksi = (st*mc)+(it*ic)/(st+it) 
    Ktot = ((st+it)*Ksi*ic)/((it*Ksi)+(st*ic))
    Conductivity(1) = Ktot
    Conductivity(2) = Ktot
  END IF
  !WRITE(Message,*)  Conductivity(1,1), Conductivity(1,2)
  !CALL Info('getWindingk', Message, Level = 5)
END SUBROUTINE getWindingk

SUBROUTINE getEpikotemFbnetk( model, n, dummyArgument,Conductivity )
  ! modules needed
  USE DefUtils
  USE CompUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  Conductivity(:)
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  REAL(KIND=dp) :: et, mc, it, ic, st, Ksi, Ktot

  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND, FOUNDFT, FOUNDST
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getEpikotemFbnetk', 'Material not found')
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getEpikotemFbnetk', 'Body Parameters not found!')
  
  et = GetConstReal(BodyParams, 'Epikotem Layer Thickness', FoundFT)
  
  mc = GetConstReal(Material, 'Material Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getEpikotemFbnetk', 'Material Heat Conductivity not found in Material section')
  
  it = GetConstReal(BodyParams, 'Net Size', Found)
  IF (.NOT. FOUND) CALL Fatal('getEpikotemFbnetk', 'Net Size not found in Body section')

  ic = GetConstReal(BodyParams, 'Cast Heat Conductivity', Found)
  IF (.NOT. FOUND) CALL Fatal('getEpikotemFbnetk', 'Cast Heat Conductivity not found in Body section')

  IF ((.NOT. FOUNDFT)) CALL Fatal('getEpikotemFbnetk', 'Epicotem Layer Thickness not found in Body section')
  
  Conductivity(1) = ((et + it) * mc * ic) / (it * mc + et * ic)   ! perp
  Conductivity(2) = (et * mc + et * ic) / (et + et)               ! par
  !WRITE(Message,*)  Conductivity(1,1), Conductivity(1,2)
  !CALL Info('getEpikotemFbnetk', Message, Level = 5)
END SUBROUTINE getEpikotemFbnetk

SUBROUTINE getFbnetVolumeFraction( model, n, dummyArgument, volumeFraction )
  ! modules needed
  USE DefUtils
  USE CompUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  volumeFraction
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody

  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Element => GetCurrentElement()
  CircuitVariableBody => Model % Bodies (Element % BodyId)
  BodyParams => CircuitVariableBody % Values
  IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('getFbnetVolumeFraction', 'Component Parameters not found!')
  
  volumeFraction = GetConstReal(BodyParams, 'Volume Fraction', Found)
  IF (.NOT. FOUND) CALL Fatal('getFbnetVolumeFraction', 'Fbnet Volume Fraction not found in Component section')
 
END SUBROUTINE getFbnetVolumeFraction


SUBROUTINE getInsulationVolumeFraction( model, n, dummyArgument, volumeFraction )
  ! modules needed
  USE DefUtils
  USE CompUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  volumeFraction
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody

  TYPE(Valuelist_t), POINTER :: Material, ComponentParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getInsulationVolumeFraction', 'Material not found')
  
  Element => GetCurrentElement()
  ComponentParams => GetComponentParams(Element)
  IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('getInsulationVolumeFraction', 'Component Parameters not found!')
  
  volumeFraction = GetConstReal(ComponentParams, 'Insulation Volume Fraction', Found)
  IF (.NOT. FOUND) CALL Fatal('getInsulationVolumeFraction', 'Insulation Volume Fraction not found in Component section')
 
END SUBROUTINE getInsulationVolumeFraction

SUBROUTINE getHeatExpCoeffComposite( model, n, dummyArgument, Coefficient )
  ! modules needed
  USE DefUtils
  USE CompUtils
 
  IMPLICIT None
 
  ! variables in function header
  TYPE(Model_t) :: model
  INTEGER :: n
  REAL(KIND=dp) :: dummyArgument
 
  ! variables needed inside function
  REAL(KIND=dp) ::  Coefficient(:)
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody
  REAL(KIND=dp) :: heatExpCoeff1, heatExpCoeff2, volumeFraction1, volumeFraction2
  REAL(KIND=dp) :: E1, E2, div

  TYPE(Valuelist_t), POINTER :: Material, ComponentParams
  TYPE(Element_t), POINTER :: Element
  Logical :: FOUND
  
  Material => GetMaterial()
  IF (.not. ASSOCIATED(Material)) CALL Fatal('getHeatExpCoeffComposite', 'Material not found')
  
  heatExpCoeff1 = GetConstReal(Material, 'Heat Expansion Coefficient Material 1', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Heat Expansion Coefficient not found in Material section')
  
  heatExpCoeff2 = GetConstReal(Material, 'Heat Expansion Coefficient Material 2', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Heat Expansion Coefficient not found in Material section')
  
  volumeFraction2 = GetConstReal(Material, 'Volume Fraction Material 2', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Volume Fraction Material 2 not found in Material section')
  
  volumeFraction1 = 1.0 - volumeFraction2
  
  E1 = GetConstReal(Material, 'Youngs Modulus Material 1', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Youngs Modulus Material 1 not found in Material section')
  
  E2 = GetConstReal(Material, 'Youngs Modulus Material 2', Found)
  IF (.NOT. FOUND) CALL Fatal('getHeatExpCoeffComposite', 'Youngs Modulus Material 2 not found in Material section')
  
  div = E1*volumeFraction1 + E2*volumeFraction2
  Coefficient(1) = (E1*heatExpCoeff1*volumeFraction1 + E2*heatExpCoeff2*volumeFraction2)/div
  Coefficient(2) = (volumeFraction1 * heatExpCoeff1) + (volumeFraction2 * heatExpCoeff2)
  
END SUBROUTINE getHeatExpCoeffComposite
