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
! *  Authors: Antti Pursula
! *  Email:   Antti.Pursula@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2003
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for repeated reading of an existing solution file.
!> The intended use for this is in postprocessing if the user needs to perform some 
!> additional step that was not done initially.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ReloadSolution( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: RestartFile, OutputName
  LOGICAL :: GotIt, ContReading
  INTEGER :: FirstStep, LastStep, StartingStep, StepsBetween
  INTEGER :: k, Round = 0, Visit = 0
  
  SAVE Round, FirstStep, LastStep, Visit, StepsBetween, ContReading

!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------
  IF ( Visit == 0 ) THEN
    FirstStep = GetInteger( Solver % Values, 'Reload Range Minimum', GotIt )
    IF ( .NOT. GotIt )  FirstStep = 1

    LastStep = GetInteger( Solver % Values, 'Reload Range Maximum', GotIt )
    IF ( .NOT. GotIt )  LastStep = 100000

    StartingStep = GetInteger( Solver % Values, 'Reload Starting Position', GotIt )
    IF ( GotIt )  Round = Round + StartingStep - FirstStep

    StepsBetween = GetInteger( Solver % Values, 'Reload Reading Intervals', GotIt )
    IF ( .NOT. GotIt )  StepsBetween = 1

    ContReading = GetLogical( Solver % Values, 'Continuous Reading', GotIt )
    IF ( .NOT. GotIt )  ContReading = .TRUE.

    IF ( ContReading ) THEN
      IF ( StartingStep /= 1 .OR. FirstStep /= 1 ) THEN
        ContReading = .FALSE.
        CALL Info( 'ReloadSolution', 'Continuous reading disabled', Level=8 )
      END IF
    END IF

  END IF

!------------------------------------------------------------------------------

  IF ( MOD( Visit, StepsBetween ) == 0 ) THEN

    RestartFile = GetString( Solver % Values, 'Reload Solution File', GotIt )

    IF ( GotIt ) THEN

      IF ( Visit == 0 ) THEN
        WRITE( Message, * ) 'Reading old solution from file ', TRIM(RestartFile)
        CALL Info( 'ReloadSolution', Message, Level=4 )
      END IF

      Round = Round + 1
      k = MOD( Round, 1+LastStep-FirstStep )
      IF ( k == 0 ) THEN
         k = LastStep
      ELSE
         k = k + FirstStep - 1
      END IF

      Mesh => Solver % Mesh

      IF ( LEN_TRIM(Mesh % Name) > 0 ) THEN
        OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
        WRITE( Message, *) 'Opening file: ', TRIM(OutputName)
      ELSE
        OutputName = TRIM(RestartFile)
        WRITE( Message, *) 'Opening file: ', TRIM(OutputName)
      END IF
      CALL Info( 'ReloadSolution', Message, Level=4 )

      IF(ParEnv % PEs > 1) OutputName = TRIM(OutputName) // '.' // i2s(ParEnv % MyPe)

      WRITE( Message, *) 'Loading Timestep: ', k
      CALL Info( 'ReloadSolution', Message, Level=4 )
      IF ( ContReading ) THEN
        CALL LoadRestartFile( OutputName, k, Mesh, .TRUE. )
      ELSE
        CALL LoadRestartFile( OutputName, k, Mesh )
      END IF

    END IF

  END IF

  Visit = Visit + 1
!------------------------------------------------------------------------------
END SUBROUTINE ReloadSolution
!------------------------------------------------------------------------------
