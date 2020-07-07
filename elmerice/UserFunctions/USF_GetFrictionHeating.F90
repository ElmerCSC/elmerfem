!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Module containing a functions for friction heat
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Juha Ruokolainen,  Joe Todd
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Current date: 20 July 2012 (Martina)
! *
! *****************************************************************************/


! contains function for adding friction heat to geothermal heat flux at the basis and exporting it
! ATTENTION: consistent units for capacity and conductivy are needed.
! ATTENTION: The keyword "Friction Heat" does make reference to Strainheating, not Friction heat!
! usage:




FUNCTION getFrictionHeat(  Model, Node, DummyInput)RESULT(frictionheat)
  
  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyInput, frictionHeat
  !----------------------------------------------------------------------------
  
  INTEGER :: DIM, i, j, Ind(3,3)
  REAL(KIND=dp), POINTER :: FlowValues(:),NormalValues(:),StressValues(:)
  REAL(KIND=dp) :: normal(3), velo(3), un, ut, Sig(3,3), Sn(3), snn, snt
  INTEGER, POINTER :: FlowPerm(:),StressPerm(:), NormalPerm(:)
  LOGICAL :: FirstTime=.TRUE.,UnFoundFatal=.TRUE.
  TYPE(Variable_t), POINTER :: FlowVar,StressVariable, NormalVar
  CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
  
  SAVE FirstTime, DIM, FunctionName,Ind
  
  IF (FirstTime) THEN
     WRITE(FunctionName, '(A)') 'getFrictionHeat'
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     DO i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
  END IF
  
  ! Get the variable velocity
  !---------------------------
  FlowVar => VariableGet( Model % Variables, 'Flow Solution',UnFoundFatal=UnFoundFatal)
  FlowPerm    => FlowVar % Perm
  FlowValues  => FlowVar % Values
  
  ! Get the stress variable
  !------------------------
  StressVariable => VariableGet( Model % Variables, 'Stress',UnFoundFatal=UnFoundFatal)
  StressPerm    => StressVariable % Perm
  StressValues  => StressVariable % Values
  
  ! Get the variable for normal vector
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
  NormalPerm => NormalVar % Perm
  NormalValues => NormalVar % Values
  
  DO i=1, DIM
     normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
  END DO
  
  !Tangential velocity
  un = SUM(velo(1:DIM)*(normal(1:DIM))) 
  ut = SQRT(SUM( (velo(1:DIM))**2.0_dp ) - (un**2.0_dp) )

  !Tangential Stress
  DO i=1, DIM
     DO j= 1, DIM
        Sig(i,j) =  &
             StressValues( 2*DIM *(StressPerm(Node)-1) + Ind(i,j) )
     END DO
  END DO
  DO i=1, DIM
     Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
  END DO
  Snn = SUM( Sn(1:DIM) * normal(1:DIM) )
  Snt = SQRT( SUM(Sn(1:DIM)**2.0_dp) - (Snn**2.0_dp))

  frictionHeat =ut*Snt
END FUNCTION getFrictionHeat

!/******************************************************************************
! *
! *  Module containing a functions for friction heat based
! *       on residuals from (Navier-)Stokes solver
! *  This function should be preferably used to compute the heat production
! *  at the base, as it utilizes the natural way to couple the flow
! *  and temperature solution
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Juha Ruokolainen, Martina Schäfer, Olivier Gagliardini
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Current date:  28 August 2014 (Martina/Thomas)
! *
! *****************************************************************************/

FUNCTION getFrictionLoads(  Model, Node, DummyInput )RESULT(frictionLoad)
  
  USE DefUtils

  IMPLICIT NONE
  
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: DummyInput, frictionLoad
  !----------------------------------------------------------------------------

  INTEGER :: DIM, i, other_body_id
  REAL(KIND=dp), POINTER :: FlowValues(:),FlowLoadValues(:),NormalValues(:),MaskValues(:)
  REAL(KIND=dp) :: normal(3), velo(3), normalvelocity, flowload(3), tangvelocity(3)
  INTEGER, POINTER :: FlowPerm(:),FlowLoadPerm(:), NormalPerm(:), MaskPerm(:)
  LOGICAL :: FirstTime=.TRUE., GotIt,UnFoundFatal,UseMask = .FALSE., Warned=.FALSE.
  TYPE(Variable_t), POINTER :: FlowVar,FlowLoadVar, NormalVar, MaskVar
  TYPE(ValueList_t), POINTER :: Equation
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolutionName, FlowLoadsName, MaskName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='USF_GetFrictionHeating(getFrictionLoads)'
  TYPE(Element_t), POINTER ::  BoundaryElement, ParentElement
  
  SAVE FirstTime,  DIM, UseMask

  
  IF (FirstTime) THEN
    !WRITE(FunctionName,'(A)') 'USF_GetFrictionHeating(getFrictionLoads)'    

    DIM = CoordinateSystemDimension()
  END IF
   
  ! Get variable names from Equation section
  !-----------------------------------------
  BoundaryElement => Model % CurrentElement
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
    ParentElement => BoundaryElement % BoundaryInfo % Right
    IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
    ParentElement => BoundaryElement % BoundaryInfo % Right
    IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  Equation => GetEquation(Element=ParentElement,Found=GotIt)
  IF (.NOT.ASSOCIATED(Equation) .OR. .NOT.GotIt) THEN
    IF (FirstTime) THEN
      WRITE (Message,'(A,I3)') 'No "Equation" found. Using default values for variables'
      CALL WARN(FunctionName,Message)    
      Warned = .TRUE.
    END IF
    WRITE(FlowSolutionName,'(A)') 'Flow Solution'
    WRITE(FlowLoadsName,'(A)') TRIM(FlowSolutionName)//' Loads'
    UseMask = .FALSE.
  ELSE
    FlowSolutionName = GetString( Equation , 'Flow Solution Name', GotIt )
    IF (.NOT.GotIt) THEN
      WRITE(FlowSolutionName,'(A)') 'Flow Solution'
      IF (FirstTime) THEN
        WRITE(Message,'(A,A)') 'Using default name for flow solution: ', &
             FlowSolutionName
        CALL WARN(FunctionName,Message)
        Warned = .TRUE.
      END IF
    END IF
    FlowLoadsName = GetString( Equation , 'Flow Loads Name', GotIt )  
    IF (.NOT. GotIt) THEN
      WRITE(FlowLoadsName,'(A)') TRIM(FlowSolutionName)//' Loads'
      IF (FirstTime) THEN
        WRITE(Message,'(A,A)') 'Using default name for flow solution loads: ', &
             FlowLoadsName
        CALL WARN(FunctionName,Message)
        Warned = .TRUE.
      END IF
    END IF
    MaskName = GetString( Equation , 'Friction Load Mask', UseMask )    
    IF (UseMask) THEN
      WRITE (Message, '(A,A)') '>Friction Load Mask< found and set to ', &
           TRIM(MaskName)
      CALL INFO(FunctionName, Message, Level=1)
    END IF
  END IF
  IF (Warned .AND. FirstTime) &
       CALL WARN(FunctionName,"All Warnings will be further omitted")
  FirstTime = .FALSE.
  
  ! Get the variable velocity
  !---------------------------
  FlowVar => VariableGet( Model % Variables, TRIM(FlowSolutionName),UnFoundFatal=UnFoundFatal)
  FlowPerm    => FlowVar % Perm
  FlowValues  => FlowVar % Values

  
  ! Get the Stokes loads
  !---------------------------
  FlowLoadVar => VariableGet( Model % Variables, TRIM(FlowLoadsName),UnFoundFatal=UnFoundFatal)
  FlowLoadPerm    => FlowLoadVar % Perm
  FlowLoadValues  => FlowLoadVar % Values
  

  ! Get the variable for normal vector
  !-----------------------------------
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
  NormalPerm => NormalVar % Perm
  NormalValues => NormalVar % Values

  IF (UseMask) THEN
    MaskVar => VariableGet( Model % Variables, TRIM(MaskName),UnFoundFatal=UnFoundFatal)
    MaskPerm    => MaskVar % Perm
    MaskValues  => MaskVar % Values
  END IF

  IF (UseMask .AND. (MaskValues(MaskPerm(Node))  .LE. 0.0_dp)) THEN
    frictionLoad = 0.0_dp
  ELSE
    DO i=1, DIM
      normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
      velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
      flowload(i) = FlowLoadValues( (DIM+1)*(FlowLoadPerm(Node)-1) + i )
    END DO
 
    normalvelocity    = SUM( velo(1:DIM) * normal(1:DIM) )

    DO i=1, DIM
      tangvelocity(i) = velo(i) - normalvelocity * normal(i)
    END DO
    
    frictionLoad = &
         MAX( (-1.0_dp * (SUM(tangvelocity(1:DIM) * flowLoad(1:DIM)))), 0.0_dp)
  END IF
END FUNCTION getFrictionLoads
