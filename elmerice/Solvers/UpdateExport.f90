
      SUBROUTINE UpdateExport( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE SolverUtils
      IMPLICIT NONE

      TYPE(Model_t) :: Model
      TYPE(Solver_t):: Solver
      REAL(KIND=dp) :: dt
      LOGICAL :: TransientSimulation

      call UpdateExportedVariables(Solver)

      END SUBROUTINE UpdateExport
