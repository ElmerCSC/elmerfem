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

