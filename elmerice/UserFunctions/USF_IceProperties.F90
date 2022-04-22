!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *
! *  Authors: Denis Cohen, Thomas Zwinger
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! Contains functions for material properties of ice:
! - IceConductivity_SI
! - IceConductivity_m_Mpa_a
! - IceCapacity_SI
! - IceCapacity_m_MPa_a
! - IcePressureMeltingPoint
!
MODULE IceProperties
  USE DefUtils
  USE Types
  
  IMPLICIT None

  CONTAINS

  !==============================================================================
  FUNCTION GetIceConductivity(temp) RESULT(cond)
  !==============================================================================


    IMPLICIT None
    REAL(KIND=dp) :: temp, cond

    
    
    cond = 9.828_dp*exp(-5.7d-03*temp)

  END FUNCTION GetIceConductivity
  
  !==============================================================================
  FUNCTION GetIceCapacity(temp) RESULT(capac)
  !==============================================================================

    USE DefUtils

    IMPLICIT None

    REAL(KIND=dp) :: temp, capac



    capac = 146.3_dp + (7.253_dp * temp)
  END FUNCTION GetIceCapacity
  
  !==============================================================================
  FUNCTION GetIcePressureMeltingPoint(ClausiusClapeyron, press) RESULT(Tpmp)
  !==============================================================================

    USE DefUtils
    
    IMPLICIT None

    REAL(KIND=dp) :: Tpmp, ClausiusClapeyron, press

    Tpmp = 273.15_dp - ClausiusClapeyron*MAX(press, 0.0_dp)

  END FUNCTION GetIcePressureMeltingPoint

END MODULE IceProperties

!==============================================================================
FUNCTION IceConductivity(Model, Node, temp) RESULT(cond)
!==============================================================================
  USE IceProperties
  
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond

  REAL(KIND=dp) :: scalingfactor
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER  :: Material
  LOGICAL :: Found
  
  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  scalingfactor = GetConstReal(Material,"Heat Conductivity Scaling Factor",Found)
  IF (.NOT.Found) scalingfactor = 1.0_dp
  cond = scalingfactor * GetIceConductivity(temp)
END FUNCTION IceConductivity

!==============================================================================
FUNCTION IceCapacity(Model, Node, temp) RESULT(capac)
!==============================================================================
  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  REAL(KIND=dp) :: scalingfactor
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER  :: Material
  LOGICAL :: Found
  
  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  scalingfactor = GetConstReal(Material,"Heat Capacity Scaling Factor",Found)
  IF (.NOT.Found) scalingfactor = 1.0_dp
  capac = scalingfactor * GetIceCapacity(temp)
END FUNCTION IceCapacity

!==============================================================================
FUNCTION IcePressureMeltingPoint(Model, Node, press) RESULT(Tpmp)
!==============================================================================

  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Tpmp, press

  INTEGER :: N
  REAL(KIND=dp) :: ClausiusClapeyron, tmpoffset
  TYPE(ValueList_t), POINTER :: Constants
  LOGICAL :: FirstTime = .TRUE., GotIt, InCelsius
  REAL(KIND=dp) :: scalingfactor
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER  :: Material
  
  SAVE FirstTime, ClausiusClapeyron

  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  scalingfactor = GetConstReal(Material,"Pressure Scaling Factor",GotIt)
  IF (.NOT.GotIt) scalingfactor = 1.0_dp
  InCelsius = GetLogical(Material, "Pressure Melting Point Celsius", GotIt)
  IF (InCelsius) THEN
    tmpoffset = 273.15_dp
  ELSE
    tmpoffset = 0.0_dp
  END IF
  IF (FirstTime) THEN
    FirstTime = .FALSE.
    Constants => GetConstants()
    IF (.NOT.ASSOCIATED(Constants)) CALL FATAL("IcePressureMeltingPoint","No Constants associated.")
    ClausiusClapeyron = GetConstReal( Constants, 'Clausius Clapeyron Constant', GotIt)
    IF (.NOT.GotIt) THEN
      ClausiusClapeyron = 9.8d-08
      CALL INFO("IcePressureMeltingPoint","No entry found for >Clausius Clapeyron Constant<.",Level=9)
      CALL INFO("IcePressureMeltingPoint","Setting to 9.8d-08 (SI units)",Level=9)
    END IF
  END IF
  Tpmp = GetIcePressureMeltingPoint(ClausiusClapeyron,scalingfactor*press) + tmpoffset

END FUNCTION IcePressureMeltingPoint

!==============================================================================
FUNCTION ArrheniusFactor(Model, Node, InputArray) RESULT(ArrhF)
!==============================================================================

  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: InputArray(2), Arrhf

  INTEGER :: N
  REAL(KIND=dp) ::  Temp, Pressure
  REAL(KIND=dp) :: Tempoffset, GlenExponent, GlenEnhancementFactor, &
       ClausiusClapeyron, GasConstant, LimitTemperature,&
        Tpmp, Trel, scalingfactor, RateFactor(2), ActivationEnergy(2)
  TYPE(ValueList_t), POINTER :: Constants
  LOGICAL ::  GotIt, InCelsius, ScaleRateFactor
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER  :: Material
  


  
  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  Constants => GetConstants()

  GasConstant = GetConstReal(Constants, "Gas Constant",GotIt)
  IF (.NOT.GotIt) GasConstant = 8.314_dp
  ClausiusClapeyron = GetConstReal( Constants, 'Clausius Clapeyron Constant', GotIt)
  IF (.NOT.GotIt) THEN
    ClausiusClapeyron = 9.8d-08
    CALL INFO("IcePressureMeltingPoint","No entry found for >Clausius Clapeyron Constant<.",Level=5)
    CALL INFO("IcePressureMeltingPoint","Setting to 9.8d-08 (SI units)",Level=5)
  END IF
  scalingfactor = GetConstReal(Material,"Pressure Scaling Factor",GotIt)
  IF (.NOT.GotIt) scalingfactor = 1.0_dp
  InCelsius = GetLogical(Material, "Temperature in Celsius", GotIt)
  IF (InCelsius) THEN
    Tempoffset = 0.0_dp
  ELSE
    Tempoffset = 273.15_dp
  END IF

  Temp = InputArray(1) - Tempoffset
  Pressure = InputArray(2)
  
  GlenExponent = GetConstReal(Material,"Glen Exponent",GotIt)
   IF (.NOT.GotIt) THEN
     CALL INFO("IceProperties(ArrheniusFactor)",&
          '"Glen Exponent" not found. Setting to 3.0',Level=5)
    GlenExponent = 3.0_dp
  END IF  
  LimitTemperature = GetConstReal(Material,"Limit Temperature",GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusFactor)",&
         '"Limit Temperature" not found. Setting to -10.0',Level=5)
    LimitTemperature = -10.0_dp
  END IF
  RateFactor(1) = GetConstReal(Material,"Rate Factor 1", GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusFactor)",&
         '"Rate Factor 1" not found. Setting to (SI value) 2.89165e-13',Level=5)
    RateFactor(1) = 2.89165d-13
  END IF
  RateFactor(2) = GetConstReal(Material,"Rate Factor 2", GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusFactor)",&
         '"Rate Factor 2" not found. Setting to (SI value) 2.42736e-02',Level=5)
    RateFactor(2) = 2.42736d-02
  END IF
  ScaleRateFactor = GetLogical(Material,"Scale Rate Factors",GotIt)
  IF (ScaleRateFactor) RateFactor = RateFactor*(scalingfactor**GlenExponent)
  
  ActivationEnergy(1) = GetConstReal(Material,"Activation Energy 1", GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusEnergy)",&
         '"Activation Energy 1" not found. Setting to (SI value) 60.0e03',Level=5)
    ActivationEnergy(1) = 60.0d03
  END IF
  ActivationEnergy(2) = GetConstReal(Material,"Activation Energy 2", GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusEnergy)",&
         '"Activation Energy 2" not found. Setting to (SI value) 115.0e3',Level=5)
    ActivationEnergy(2) = 115.d03
  END IF
 
  GlenEnhancementFactor = GetConstReal(Material,"Glen Enhancement Factor",GotIt)
  IF (.NOT.GotIt) THEN
    CALL INFO("IceProperties(ArrheniusEnergy)",&
         '"Glen Enhancement Factor" not found. Setting to 1.0', Level=5)
    GlenEnhancementFactor = 1.0_dp
  END IF

  ! need to scale pressure 
  Tpmp = GetIcePressureMeltingPoint(ClausiusClapeyron, scalingfactor * Pressure)
   
  Trel = MAX(MIN(Temp - Tpmp,0.0_dp),-60.0_dp)
  IF (Trel<-10) THEN
    Arrhf = RateFactor(1) * EXP( -ActivationEnergy(1)/( GasConstant* (273.15 + Trel)))
  ELSE
    Arrhf= RateFactor(2) * EXP( -ActivationEnergy(2)/( GasConstant* (273.15 + Trel)))
  END IF
END FUNCTION ArrheniusFactor
