
RECURSIVE SUBROUTINE DummySolver( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
 
  IMPLICIT NONE


  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: Var
  LOGICAL :: CalvingOccurs, Found

  PRINT *,"Calving 3D Test Solver"
  PRINT *,"***********************************"

  Var => Solver % Variable
  IF (.NOT. ASSOCIATED(Var)) THEN
     CALL FATAL('DummySolver','No Variable associated')
  END IF

  CalvingOccurs = ListGetLogical(Model % Simulation, 'CalvingOccurs', Found)
  IF(.NOT. Found) CalvingOccurs = .FALSE.

  IF(CalvingOccurs) THEN
    Var % Norm = 1
  ELSE
    Var % Norm = 0
  END IF
END SUBROUTINE DummySolver
