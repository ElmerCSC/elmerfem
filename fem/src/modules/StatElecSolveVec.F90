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
! *  This solver's assembly moved into StatElecSolve, which now IS it: the
! *  vectorized/threaded implementation that used to live in this file, under
! *  this file's own name, is the default static electric field solver. What
! *  is left here is a front, kept only so that a sif naming
! *  "StatElecSolveVec" "StatElecSolver" -- as tests written against the old
! *  name still do -- keeps working unchanged.
! *
! *  The delegation goes through GetProcAddr and ExecSolver, which is how the
! *  core itself invokes a solver, so nothing links this file against
! *  StatElecSolve.so or StatElecSolveLegacy.so, and any of the three can be
! *  rebuilt alone.
! *
! *  "Legacy Assembly = Logical True" in the solver section reaches the
! *  original scalar-element solver (StatElecSolveLegacy.F90) instead, exactly
! *  as it does when named directly through "StatElecSolve".
! *
! *  "Permittivity Of Vacuum" convention: StatElecSolve.F90 now defaults this
! *  to 1.0 (natural/relative units) when a sif does not set it. Sifs naming
! *  THIS file directly predate that and were tuned against the real physical
! *  default, so this front injects that value as a default (only if the sif
! *  has not already set its own) before delegating to the engine, on every
! *  path that reaches it -- keeping those sifs' behavior unchanged. A sif
! *  reached through "Legacy Assembly" does not need this: StatElecSolveLegacy
! *  has always defaulted to 1.0 on its own, unaffected by any of this.
! *
! *****************************************************************************/
MODULE StatElecSolveVecFront
  USE DefUtils
  USE LoadMod, ONLY: ExecSolver
  IMPLICIT NONE

  REAL(KIND=dp), PARAMETER :: PhysicalEps0 = 8.854187817e-12_dp

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
!> Default "Permittivity Of Vacuum" to the real physical value, only if the
!> sif has not already set its own -- see the file header.
!------------------------------------------------------------------------------
  SUBROUTINE InjectPhysicalEps0( Model )
    TYPE(Model_t) :: Model

    IF( .NOT. ASSOCIATED( Model % Constants ) ) THEN
      Model % Constants => ListAllocate()
    END IF
    CALL ListAddNewConstReal( Model % Constants,'Permittivity Of Vacuum', PhysicalEps0 )
  END SUBROUTINE InjectPhysicalEps0

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
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'StatElecSolver', &
        'This solver is a front for '//TRIM(TargetFile)//' and "'//TRIM(Entry)// &
        '" could not be found. Is '//TRIM(TargetFile)//'.so installed beside '// &
        'StatElecSolveVec.so?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateTo

!------------------------------------------------------------------------------
!> Which file to delegate one call to: the original scalar-element solver when
!> "Legacy Assembly" is requested, the vectorized/threaded solver otherwise --
!> injecting the physical default Eps0 for the latter, on every entry point,
!> since which one runs first is not this file's to assume (see StressSolve's
!> own note on "_Init0" not always running).
!------------------------------------------------------------------------------
  SUBROUTINE Delegate( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    IF ( LegacyAssembly( Solver ) ) THEN
      CALL DelegateTo( 'StatElecSolveLegacy', 'StatElecSolverLegacy'//Entry, &
          Model, Solver, dt, Transient )
    ELSE
      CALL InjectPhysicalEps0( Model )
      CALL DelegateTo( 'StatElecSolve', 'StatElecSolver'//Entry, &
          Model, Solver, dt, Transient )
    END IF
  END SUBROUTINE Delegate

END MODULE StatElecSolveVecFront


!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_Init0( Model,Solver,dt,Transient )
  USE StatElecSolveVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  ! Legacy has no complex/harmonic-mode setup at all, so there is nothing to
  ! delegate to for it here -- just skip it, same as StatElecSolve.F90 does.
  IF ( LegacyAssembly( Solver ) ) RETURN
  CALL Delegate( '_Init0', Model, Solver, dt, Transient )
END SUBROUTINE StatElecSolver_Init0


!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_init( Model,Solver,dt,Transient )
  USE StatElecSolveVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL Delegate( '_Init', Model, Solver, dt, Transient )
END SUBROUTINE StatElecSolver_init


!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver( Model,Solver,dt,Transient )
  USE StatElecSolveVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL Delegate( '', Model, Solver, dt, Transient )
END SUBROUTINE StatElecSolver


!------------------------------------------------------------------------------
!> The vectorized solver's optional postprocessing step. Only ever reached
!> when NOT running under "Legacy Assembly" (StatElecSolve's own _init only
!> sets "PostSolver Active" on that path), but delegates properly either way in
!> case that keyword gets set some other way.
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_post( Model,Solver,dt,Transient )
  USE StatElecSolveVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  IF ( LegacyAssembly( Solver ) ) RETURN
  CALL InjectPhysicalEps0( Model )
  CALL DelegateTo( 'StatElecSolve', 'StatElecSolver_post', Model, Solver, dt, Transient )
END SUBROUTINE StatElecSolver_post
