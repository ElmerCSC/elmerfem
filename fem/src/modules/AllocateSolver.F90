!------------------------------------------------------------------------------
SUBROUTINE AllocateSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------

  CALL ListAddNewString(Solver % Values,'Exec Solver','never')

!------------------------------------------------------------------------------
END SUBROUTINE AllocateSolver_init
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!> This is a dummy solver with the sole purpose of allocating space for
!> other solvers.
!------------------------------------------------------------------------------
SUBROUTINE AllocateSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------

  CALL Info('AllocateSolver','This solver does nothing, it is just for allocating stuff')

!------------------------------------------------------------------------------
END SUBROUTINE AllocateSolver
!------------------------------------------------------------------------------
