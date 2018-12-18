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
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2014/02/20 - Start wrinting User functions needed for the coupling between
! * the hydrology model (GlaDS) and the ice flow (SSA).
! *****************************************************************************
!> USF_CouplingGlaDS_SSA.F90
!> 
!> HorizontalVelo returns the horizontal velocity (positive)
!> 
FUNCTION HorizontalVelo (Model, nodenumber, x) RESULT(ub)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: x, ub          
   INTEGER :: nodenumber
       
   TYPE(ValueList_t), POINTER :: Constants
   TYPE(Variable_t), POINTER :: FlowVariable, Gravity
   REAL(KIND=dp), POINTER :: FlowValues(:)
   INTEGER, POINTER :: FlowPerm(:)
   INTEGER :: DIM, i
   REAL (KIND=dp) :: velo(2)
   LOGICAL :: GotIt

   CHARACTER(LEN=MAX_NAME_LEN) :: USFName

   USFName = 'HorizontalVelo'
   DIM = CoordinateSystemDimension()

   ! Get the variables to compute ub : the fields velocity 1 and 2
   FlowVariable => VariableGet( Model % Variables, 'ssavelocity' )
   
   IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
   ELSE
      CALL FATAL(USFName, 'Need SSA Solver, variable >SSAVelocity< not associated !!!')
   END IF

   velo = 0.0_dp
   DO i=1, DIM
     velo(i) = FlowValues( (DIM)*(FlowPerm(Nodenumber)-1) + i ) 
   END DO
   ub = SQRT(velo(1)**2+velo(2)**2)
 END FUNCTION HorizontalVelo
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Function OverburdenPressure return the cryostatic pressure 
!> 
 FUNCTION OverburdenPressure (Model, nodenumber, x) RESULT(IcePress)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: x, IcePress          
   INTEGER :: nodenumber
   TYPE(ValueList_t), POINTER :: Constants
   TYPE(Variable_t), POINTER :: ZbVariable, ZsVariable
   REAL(KIND=dp), POINTER :: ZbValues(:), ZsValues(:)
   INTEGER, POINTER :: ZbPerm(:), ZsPerm(:)
   INTEGER :: DIM, i
   REAL (KIND=dp) :: zb, zs, gravity, rhoi
   LOGICAL :: GotIt

   CHARACTER(LEN=MAX_NAME_LEN) :: USFName
   
   USFName = 'OverburdenPressure'
   DIM = CoordinateSystemDimension()
   
   ! Get the constants needed to compute the overburden ice pressure
   Constants => GetConstants()
   gravity = GetConstReal( Constants, 'Gravity Norm', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(a)')'Keyword >Gravity Norm< not found in Constant section'
      CALL FATAL(USFName, Message)
   END IF
   rhoi = GetConstReal ( Constants, 'Ice Density', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(a)')'Keyword >Ice Density< not found in Constant section'
      CALL FATAL(USFName, Message)
   END IF

   ! Get the variables to compute IcePress
   ZbVariable => VariableGet( Model % Variables, 'zb' )
   IF ( ASSOCIATED( ZbVariable ) ) THEN
      ZbPerm    => ZbVariable % Perm
      ZbValues  => ZbVariable % Values
   ELSE
      CALL FATAL(USFName, 'Variable >Zb< not associated !!!')
   END IF
     
   ! Get the variables to compute IcePress 
   ZsVariable => VariableGet( Model % Variables, 'zs' )
   IF ( ASSOCIATED( ZsVariable ) ) THEN
      ZsPerm    => ZsVariable % Perm
      ZsValues  => ZsVariable % Values
   ELSE
      CALL FATAL(USFName, 'Variable >Zs< not associated !!!')
   END IF
   zb = 0.0_dp
   zs = 0.0_dp

   zb = ZbValues( ZbPerm(Nodenumber))
   zs = ZsValues( ZsPerm(Nodenumber)) 
   IcePress = rhoi*gravity*(zs-zb)
   IF (IcePress .LT. 0) THEN
     print*, "Ice pressure value :", IcePress
   END IF
   END FUNCTION OverburdenPressure
