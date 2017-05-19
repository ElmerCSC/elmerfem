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
! *  A Module for treating scanned field solutions.
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
SUBROUTINE ScannedFieldSolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  INTEGER :: i, FieldInt, ScanMax, ScanInt, ScanSolverInt, NofFields, istat
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: FieldName
  LOGICAL :: Found

  SolverParams => GetSolverParams()

  i=1
  DO WHILE(.TRUE.)
    IF ( .NOT.ListCheckPresent(SolverParams, &
          "Exported Variable "//TRIM(i2s(i))) ) EXIT
    i = i + 1
  END DO

  NofFields = 1
  DO WHILE(.TRUE.)
    IF ( .NOT.ListCheckPresent(SolverParams, &
          "Field Name "//TRIM(i2s(NofFields))) ) EXIT
    NofFields = NofFields + 1
  END DO
  NofFields = NofFields - 1

  CALL ListAddString( SolverParams, 'Variable', '-nooutput SF_dummy' )
  CALL ListAddInteger( SolverParams, 'Number of Fields', NofFields)

  ScanSolverInt = GetInteger(SolverParams,'Scan Solver',Found) 
  IF (.NOT. Found) CALL Fatal('ScannedFieldSolver_Init','Scan Solver not found.')

  ScanMax = GetInteger(Model % Solvers(ScanSolverInt) % Values, 'Scanning Loops', Found)
  IF (.NOT. Found) CALL Fatal('ScannedFieldSolver_Init','Scan Loops not found.')

  CALL ListAddInteger(SolverParams, 'Scan Max', ScanMax)

  DO FieldInt = 1, NofFields
    FieldName = GetString(SolverParams,'Field Name '//TRIM(i2s(FieldInt)),Found) 
    IF (.NOT. Found) CALL Fatal('ScannedFieldSolver_Init','Field Name '//TRIM(i2s(FieldInt))//' not found.')

    DO ScanInt = 1, ScanMax 
      CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
        'Scan '//TRIM(i2s(ScanInt))//" "//TRIM(FieldName))
      i=i+1
    END DO

    CALL ListAddString( SolverParams, "Exported Variable "//TRIM(i2s(i)), &
      '-nooutput '//TRIM(FieldName)//' Dummy')
    i=i+1
 
  END DO

!------------------------------------------------------------------------------
END SUBROUTINE ScannedFieldSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE ScannedFieldSolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: SumFieldVar, ScanVar, ScanFieldVar, FieldVar
  TYPE VariablePtr_t
    TYPE(Variable_t), POINTER :: Var
  END TYPE VariablePtr_t
  TYPE(VariablePtr_t), POINTER :: FieldVars(:)
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: ASolver => NULL()
  INTEGER :: N, ScanInt, istat, NofFields, FieldInt, ScanMax 
  LOGICAL :: Found, First=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: FieldNames(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: ScanFieldName, SumFieldName 
 
  TYPE(ValueList_t), POINTER :: SolverParam

  SAVE FieldVars, FieldNames, NofFields, Mesh
!------------------------------------------------------------------------------
 
  IF (First) THEN
    First = .FALSE.
 
    Mesh => Model % Mesh
    N = Mesh % MaxElementDOFs

    NofFields = ListGetInteger(GetSolverParams(), 'Number of Fields')

    ALLOCATE(FieldNames(NofFields), FieldVars(NofFields), STAT=istat)
    IF ( istat /= 0 ) CALL Fatal('ScannedFieldSolver','Memory allocation error')
    
    DO FieldInt = 1, NofFields
      FieldNames(FieldInt) = ListGetString(GetSolverParams(), 'Field Name '//TRIM(i2s(FieldInt)))
      FieldVars(FieldInt) % Var => VariableGet( Mesh % Variables, FieldNames(FieldInt))
      IF(.NOT. ASSOCIATED(FieldVars(FieldInt) % Var)) THEN
        CALL Fatal('ScannedFieldSolver',TRIM(FieldNames(FieldInt))//' not associated!')
      END IF
    END DO

  END IF
 
  CALL Info('ScannedFieldSolver','-------------------------------------------',Level=6)
  CALL Info('ScannedFieldSolver','Computing the scan field solutions. ',Level=5)
  CALL Info('ScannedFieldSolver','-------------------------------------------',Level=6)

  ScanVar => VariableGet(Solver % Mesh % Variables, 'scan')
  IF (.NOT. ASSOCIATED(ScanVar)) CALL Fatal('ScannedFieldSolver', 'Scan variable not found.')
  ScanInt = INT( ScanVar % Values(1) ) 

  ScanMax = GetInteger(GetSolverParams(), 'Scan Max', Found)
  IF (.NOT. Found) CALL Fatal('ScannedFieldSolver', 'The maximum number of scan loops not found.')

  DO FieldInt = 1, NofFields
    FieldVar => FieldVars(FieldInt) % Var 
    CALL Info('ScannedFieldSolver','Field '//TRIM(FieldNames(FieldInt))//'.',Level=5)
    ScanFieldName = 'scan '//TRIM(i2s(ScanInt))//" "//TRIM(FieldNames(FieldInt))
    ScanFieldVar => VariableGet( Mesh % Variables, ScanFieldName)
    IF(.NOT. ASSOCIATED(ScanFieldVar)) THEN
      CALL Fatal('ScannedFieldSolver', TRIM(ScanFieldName)//' not associated!')
    END IF

    CALL MakeVarSimilarToModelVar(ScanFieldVar, FieldVar) 
    IF (SIZE(ScanFieldVar % Values) .NE. SIZE(FieldVar % Values)) &
      CALL Fatal('ScannedFieldSolver','Scanned fields are of different size than &
      the defined Scan Field Variable.')

    SumFieldName = TRIM(FieldNames(FieldInt))//' Dummy'
    SumFieldVar => VariableGet( Mesh % Variables, SumFieldName)
    IF(.NOT. ASSOCIATED(ScanFieldVar)) THEN
      CALL Fatal('ScannedFieldSolver', TRIM(SumFieldName)//' not associated')
    END IF

    CALL MakeVarSimilarToModelVar(SumFieldVar, FieldVar) 
    IF (SIZE(ScanFieldVar % Values) .NE. SIZE(SumFieldVar % Values)) &
      CALL Fatal('ScannedFieldSolver','Summed fields are of different size than &
      the defined Sum Field Variable.')

    ScanFieldVar % Values = FieldVar % Values
    IF( ScanInt == 1 ) THEN
      SumFieldVar % Values = FieldVar % Values
    ELSE
      SumFieldVar % Values = FieldVar % Values + SumFieldVar % Values
    END IF 

    IF ( ScanInt == ScanMax ) THEN
      FieldVar % Values = SumFieldVar % Values
    END IF

  END DO

CONTAINS

  SUBROUTINE MakeVarSimilarToModelVar(Var, ModelVar) 
    TYPE(Variable_t), POINTER :: Var, ModelVar
    INTEGER :: istat

    IF (SIZE(Var % Values) .NE. SIZE(ModelVar % Values)) THEN
      DEALLOCATE(Var % Values)
      ALLOCATE(Var % Values(SIZE(ModelVar % Values)), STAT=istat)
      IF ( istat /= 0 ) CALL Fatal('MakeVarSimilarToModelVar','Memory allocation error')
    END IF

    IF (.NOT. ASSOCIATED(Var % Perm,  ModelVar % Perm)) Var % Perm => ModelVar % Perm
    IF ( Var % DOFs .NE. ModelVar % DOFs) Var % DOFs = ModelVar % DOFs
    IF ( Var % TYPE .NE. ModelVar % TYPE ) Var % TYPE = ModelVar % TYPE

  END SUBROUTINE MakeVarSimilarToModelVar
 
!------------------------------------------------------------------------------
END SUBROUTINE ScannedFieldSolver
!------------------------------------------------------------------------------
