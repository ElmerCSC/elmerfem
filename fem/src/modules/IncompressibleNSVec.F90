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
! *  This solver's code moved into IncompressibleNS.F90, under that plain
! *  name -- there never was a non-vectorized "IncompressibleNS" to
! *  distinguish it from, so unlike StatCurrentSolve/StatElecSolve there is no
! *  "Legacy Assembly" here, only a straight rename. What is left in this file
! *  is a front, kept only so that a sif naming "IncompressibleNSVec"
! *  "IncompressibleNSSolver" -- as the many existing tests do -- keeps
! *  working unchanged.
! *
! *  The delegation goes through GetProcAddr and ExecSolver, which is how the
! *  core itself invokes a solver, so nothing links this file against
! *  IncompressibleNS.so, and either can be rebuilt alone.
! *
! *****************************************************************************/
MODULE IncompressibleNSVecFront
  USE DefUtils
  USE LoadMod, ONLY: ExecSolver
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Call one entry point of IncompressibleNS.so with this solver. The name is
!> resolved at run time, as the core resolves any solver, so this file stays
!> independent of it.
!------------------------------------------------------------------------------
  SUBROUTINE DelegateTo( Entry, Model, Solver, dt, Transient )
    CHARACTER(LEN=*) :: Entry
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: dt
    LOGICAL :: Transient

    TYPE(C_FUNPTR) :: Proc

    Proc = GetProcAddr( 'IncompressibleNS '//TRIM(Entry), abort = .FALSE. )
    IF ( .NOT. C_ASSOCIATED( Proc ) ) CALL Fatal( 'IncompressibleNSSolver', &
        'This solver is a front for IncompressibleNS and "'//TRIM(Entry)// &
        '" could not be found. Is IncompressibleNS.so installed beside '// &
        'IncompressibleNSVec.so?' )

    CALL ExecSolver( Proc, Model, Solver, dt, Transient )
  END SUBROUTINE DelegateTo

END MODULE IncompressibleNSVecFront


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver_Init0( Model,Solver,dt,Transient )
  USE IncompressibleNSVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL DelegateTo( 'IncompressibleNSSolver_Init0', Model, Solver, dt, Transient )
END SUBROUTINE IncompressibleNSSolver_Init0


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver_init( Model,Solver,dt,Transient )
  USE IncompressibleNSVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL DelegateTo( 'IncompressibleNSSolver_Init', Model, Solver, dt, Transient )
END SUBROUTINE IncompressibleNSSolver_init


!------------------------------------------------------------------------------
SUBROUTINE IncompressibleNSSolver( Model,Solver,dt,Transient )
  USE IncompressibleNSVecFront
  IMPLICIT NONE
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  CALL DelegateTo( 'IncompressibleNSSolver', Model, Solver, dt, Transient )
END SUBROUTINE IncompressibleNSSolver
