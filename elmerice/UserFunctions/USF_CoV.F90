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
! *  Authors: F. Gillet-Chaulet (IGE-Grenoble-FR)
! *           R. Gladstone (Uni Lapland)
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *  Date modifications:
! *
! * 
! *****************************************************************************
!#######################################################################
!#
!#  A collection of USER FUNCTIONS to perform variable changes in
!inverse methods; i.e. beta=10^a or beta=a^2
!#
!#######################################################################
!# Compute VarOut=10^VarIn
FUNCTION TenPowerA(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut = 10._dp**(VarIn)
  
End FUNCTION TenPowerA
!# Compute DTenPowerA/DA=ln(10)*10^A
!# VarIn=A
FUNCTION TenPowerA_d(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut = (10.0**(VarIn))*log(10.0)
  
End FUNCTION TenPowerA_d
!# Compute DJDA from DJDB if B=10^A: DJDA=DJDB*ln(10)*10^A
!# DJDB=VarIn(1)
!# A=VarIn(2)
FUNCTION Derivative_TenPowerA(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut
  
  VarOut = VarIn(1)*(10.0**(VarIn(2)))*log(10.0)
  
End FUNCTION Derivative_TenPowerA
!# Compute VarOut=Log10(VarIn)
FUNCTION Log10A(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut=log10(VarIn)
  
End FUNCTION Log10A
!# Compute VarOut=VarIn*VarIn
FUNCTION Asquare(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut = VarIn*VarIn
END FUNCTION Asquare
!# Compute Compute dA^2/dA=2*A
!# VarIn=A
FUNCTION Asquare_d(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut = 2*VarIn
END FUNCTION Asquare_d
!# Compute DJDA from DJDB if B=A^2: DJDA=DJDB*2A
!# DJDB=VarIn(1)
!# A=VarIn(2)
FUNCTION Derivative_Asquare(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn(2),VarOut
  
  VarOut = 2.0*VarIn(1)*VarIn(2)
  
End FUNCTION Derivative_Asquare
!# Compute VarOut=sqrt(VarIn)
FUNCTION SQRTA(Model,nodenumber,VarIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: VarIn,VarOut
  
  VarOut = sqrt(VarIn)
END FUNCTION SQRTA

!# Compute VarOut=10^VarIn, masked
! Input arguments:
!  ArgIn(1)  VarIn, the variable upon which to operate
!  ArgIn(2)  mask, a mask variable, typically GroundedMask
FUNCTION TenPowerA_masked(Model,nodenumber,ArgIn) RESULT(VarOut)
  USE DefUtils
  IMPLICIT none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: ArgIn(2),VarOut

  REAL(kind=dp) :: VarIn,mask

  VarIn = ArgIn(1)
  mask  = ArgIn(2)

  ! Note that the GroundedMask default behaviour is that -1.0 indicates
  ! floating shelves, 0.0 indicates grounding line nodes and 1.0
  ! indicates grounded nodes.
  IF (mask.LE.0.0_dp) THEN
    VarOut = 0.0_dp
  ELSE
    VarOut = 10._dp**(VarIn)
  END IF
END FUNCTION TenPowerA_Masked

!# Compute DTenPowerA/DA=ln(10)*10^A, masked
! Input arguments:
!  ArgIn(1)  VarIn, the variable upon which to operate
!  ArgIn(2)  mask, a mask variable, typically GroundedMask
FUNCTION TenPowerA_d_Masked(Model,nodenumber,ArgIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: ArgIn(2),VarOut

  REAL(kind=dp) :: VarIn,mask

  VarIn = ArgIn(1)
  mask  = ArgIn(2)
  
  IF (mask.LE.0.0_dp) THEN
    VarOut = 0.0_dp
  ELSE
    VarOut = (10.0**(VarIn))*log(10.0)
  END IF
  
END FUNCTION TenPowerA_d_Masked

!# This function can be used, for example, when computing
!# viscosity as a function of enhancement factor and an
!# initial guess.
FUNCTION Asquare_Scaled(Model,nodenumber,ArgIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: ArgIn(2),VarOut

  REAL(kind=dp) :: VarIn
  REAL(kind=dp) :: VarScale ! a scaling variable
  
  VarIn    = ArgIn(1)
  VarScale = ArgIn(2)

  VarOut   = VarIn*VarIn*VarScale

END FUNCTION Asquare_Scaled

!# This function can be used, for example, for differentiating
!# by parts the viscosity, when computing viscosity as a
!# function of enhancement factor and an initial guess.
FUNCTION Asquare_d_Scaled(Model,nodenumber,ArgIn) RESULT(VarOut)
  USE DefUtils
  implicit none
  !-----------------
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(kind=dp) :: ArgIn(2),VarOut

  REAL(kind=dp) :: VarIn
  REAL(kind=dp) :: VarScale ! a scaling variable

  VarIn    = ArgIn(1)
  VarScale = ArgIn(2)
  
  VarOut = 2*VarIn*VarScale
  
END FUNCTION Asquare_d_Scaled

