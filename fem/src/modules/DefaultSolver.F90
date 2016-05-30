SUBROUTINE DefaultSolver_prebulk( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  CALL DefaultInitialize()

END SUBROUTINE DefaultSolver_prebulk

SUBROUTINE DefaultSolver_postbulk( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  CALL DefaultFinishBulkAssembly()

END SUBROUTINE DefaultSolver_postbulk

SUBROUTINE DefaultSolver_postboundary( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------

  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

END SUBROUTINE DefaultSolver_postboundary

SUBROUTINE DefaultSolver_solve( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  REAL(KIND=dp) :: Norm

  Norm = DefaultSolve()

END SUBROUTINE DefaultSolver_solve

