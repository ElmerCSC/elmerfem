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
! *
! * A dummy solver to generate a pressure preconditioning matrix for block
! * preconditioning. Some default initializations for using this especially in 
! * connection with the ParStokes solver are introduced to simplify the writing 
! * of sif files. 
! *
! ******************************************************************************
! *
! *  Authors: Mika Malinen & Juha Ruokolainen
! *  Email:   mika.malinen@csc.fi & Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2009-06-24
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
SUBROUTINE PressurePrecond_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  CALL ListAddString(SolverParams, 'Exec Solver', 'before simulation' )
  CALL ListAddLogical(SolverParams, 'Variable Output', .FALSE.) 
  CALL ListAddNewLogical(SolverParams, 'Bubbles in Global System', .FALSE.)  
  CALL ListAddLogical(SolverParams, 'Skip Compute Nonlinear Change', .TRUE.) 
  CALL ListAddLogical(SolverParams, 'Back Rotate N-T Solution', .FALSE.) 

  CALL ListAddNewString(SolverParams, 'Variable', 'P')

  CALL ListAddNewString(SolverParams, 'Linear System Solver', 'Iterative')
  CALL ListAddNewString(SolverParams, 'Linear System Iterative Method', 'BiCGStab2') 
  CALL ListAddNewInteger(SolverParams, 'Linear System Max Iterations', 1000)
  CALL ListAddNewString(SolverParams, 'Linear System Preconditioning', 'Diagonal') 
  CALL ListAddNewConstReal(SolverParams, 'Linear System Convergence Tolerance', 1.0d-6)
  CALL ListAddNewLogical(SolverParams, 'Linear System Abort Not Converged', .FALSE.)

!------------------------------------------------------------------------------
END SUBROUTINE PressurePrecond_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE PressurePrecond( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------

!  print *, 'A dummy solver subroutine for generating pressure preconditioning matrix'

!------------------------------------------------------------------------------
END SUBROUTINE PressurePrecond
!------------------------------------------------------------------------------
