!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  THE FRONT.
! *
! *  This solver's assembly moved into HeatSolve, which now IS it: the
! *  vectorized/threaded implementation that used to live in this file, under
! *  this file's own name, is the default heat equation solver. What is left
! *  here is a front, kept only so that a sif naming "HeatSolveVec"
! *  "HeatSolver" -- as tests written against the old name still do -- keeps
! *  working unchanged.
! *
! *  The delegation goes through GetProcAddr and ExecSolver, which is how the
! *  core itself invokes a solver, so nothing links this file against
! *  HeatSolve.so or HeatSolveLegacy.so, and any of the three can be rebuilt
! *  alone.
! *
! *  "Legacy Assembly = Logical True" in the solver section reaches the
! *  original scalar-element solver (HeatSolveLegacy.F90) instead, exactly as
! *  it does when named directly through "HeatSolve".
! *
! *****************************************************************************/
MODULE HeatSolveVecFront
  USE DefUtils
  USE LoadMod, ONLY: ExecSolver
  IMPLICIT NONE IMPLICIT_EXTERNAL

CONTAINS

!------------------------------------------------------------------------------
!> Whether this sif asks for the original scalar-element assembly rather than
!> the front.
!------------------------------------------------------------------------------
  FUNCTION LegacyAssembly( Solver ) RESULT( Legacy )
    TYPE(Solver_t) :: Solver
    LOGICAL :: Legacy, Found

    Legacy = ListGetLogical( Solver % Values, 'Legacy Assembly', Found )
  END FUNCTION LegacyAssembly

!------------------------------------------------------------------------------
!> Call one entry point of the given target solver library with this solver.
!> The name is resolved at run time, as the core resolves any solver, so this
!> file stays independent of both delegation targets.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateTo( TargetFile, Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: TargetFile, Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( TRIM(TargetFile)//' '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'HeatSolver', &
        'This solver is a front for '//TRIM(TargetFile)//' and "'//TRIM(Entry)// &
        '" could not be found. Is '//TRIM(TargetFile)//'.so installed beside '// &
        'HeatSolveVec.so?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateTo

!------------------------------------------------------------------------------
!> Which file to delegate one call to: the original scalar-element solver when
!> "Legacy Assembly" is requested, the vectorized/threaded solver otherwise.
!------------------------------------------------------------------------------
  SUBROUTINE Delegate( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    IF ( LegacyAssembly( Solver ) ) THEN
      CALL DelegateTo( 'HeatSolveLegacy', 'HeatSolverLegacy'//Entry, &
          Model, Solver, dt, Transient )
    ELSE
      CALL DelegateTo( 'HeatSolve', 'HeatSolver'//Entry, &
          Model, Solver, dt, Transient )
    END IF
  END SUBROUTINE Delegate

END MODULE HeatSolveVecFront


!------------------------------------------------------------------------------
!> HeatSolve's own pre-pass that sets up a p-bubble/SUPG "Element" string for
!> convection. Legacy has no such hook, so under "Legacy Assembly" this is a
!> no-op, exactly as if the sif had named "HeatSolveLegacy" directly (which
!> the core would not find a "_Init0" entry point for at all).
!------------------------------------------------------------------------------
SUBROUTINE HeatSolver_Init0( Model,Solver,dt,Transient )
  USE HeatSolveVecFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  IF ( LegacyAssembly( Solver ) ) RETURN
  CALL DelegateTo( 'HeatSolve', 'HeatSolver_Init0', Model, Solver, dt, Transient )
END SUBROUTINE HeatSolver_Init0


!------------------------------------------------------------------------------
SUBROUTINE HeatSolver_init( Model,Solver,dt,Transient )
  USE HeatSolveVecFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL Delegate( '_Init', Model, Solver, dt, Transient )
END SUBROUTINE HeatSolver_init


!------------------------------------------------------------------------------
SUBROUTINE HeatSolver( Model,Solver,dt,Transient )
  USE HeatSolveVecFront
  IMPLICIT NONE IMPLICIT_EXTERNAL
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL Delegate( '', Model, Solver, dt, Transient )
END SUBROUTINE HeatSolver
