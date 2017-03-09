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
! *  A Module for computing the sum of scanned fields.
! *
! *  Author:  Eelis Takala
! *  Email:   eelis.takala@gmail.com
! *  Web:     http://www.trafotek.fi
! *  Address: Trafotek
! *           Kaarinantie 700
! *           20540 Turku
! *
! *  Original Date: 9.3.2017
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
SUBROUTINE SumFieldSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: SumFieldVar, FieldVar, ScanVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: ASolver => NULL()
  INTEGER :: N, ScanInt
  LOGICAL :: Found, First=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN) :: FieldName
 
  TYPE(ValueList_t), POINTER :: SolverParam

  SAVE FieldVar, SumFieldVar
!------------------------------------------------------------------------------
 
  IF (First) THEN
    First = .FALSE.
 
    Mesh => Model % Mesh
    N = Mesh % MaxElementDOFs

    SumFieldVar => VariableGet( Mesh % Variables, Solver % Variable % Name )
    IF(.NOT. ASSOCIATED(SumFieldVar)) THEN
      CALL Fatal('FieldSumSolver',TRIM(Solver % Variable % Name)//' not associated!')
    END IF
    
    FieldName = ListGetString(GetSolverParams(), 'Field Name')
    FieldVar => VariableGet( Mesh % Variables, FieldName)
    IF(.NOT. ASSOCIATED(FieldVar)) THEN
      CALL Fatal('FieldSumSolver',TRIM(FieldName)//' not associated!')
    END IF

  END IF
 
  ScanVar => VariableGet(Solver % Mesh % Variables, 'scan')
  IF (.NOT. ASSOCIATED(ScanVar)) CALL Fatal('SumFieldSolver', 'Scan variable not found.')
  ScanInt = INT( ScanVar % Values(1) ) 

  IF (SIZE(SumFieldVar % Values) == SIZE(FieldVar % Values)) &
    CALL Fatal('SumFieldSolver','Summed fields are of different size than &
    the defined Sum Field Variable.')

  IF( ScanInt == 1 ) THEN
    SumFieldVar % Values = FieldVar % Values
  ELSE 
    SumFieldVar % Values = FieldVar % Values + SumFieldVar % Values
  END IF 

!------------------------------------------------------------------------------
END SUBROUTINE SumFieldSolver
!------------------------------------------------------------------------------
