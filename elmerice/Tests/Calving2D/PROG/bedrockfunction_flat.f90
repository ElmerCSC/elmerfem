
RECURSIVE SUBROUTINE getBedrock( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
 
  IMPLICIT NONE
  INTERFACE
     FUNCTION initbedrock(Model, nodenumber, x) RESULT(elevation)
       USE DefUtils
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL (KIND=dp) :: x, elevation
     END FUNCTION initbedrock
  END INTERFACE

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
  TYPE(ValueList_t), Pointer :: BC
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: VarPerm(:)
  INTEGER :: VarDOFs, i, k
  REAL(KIND=dp) :: x
  REAL(KIND=dp), POINTER :: VarValues(:)
  LOGICAL :: GotIt

  CALL INFO("getBedrock","Computing bedrock distribution", Level=1)

  Var => Solver % Variable
  IF (ASSOCIATED(Var)) THEN
     VarPerm => Var % Perm
     VarDOFs =  Var % DOFs
     VarValues => Var % Values
  ELSE
     CALL FATAL('getBedrock','No Variable associated')
  END IF
  k=0
  DO i = 1,Model % NumberOfNodes
     IF (VarPerm(i) > 0) THEN
        x = Solver % Mesh % Nodes % x(i)
        VarValues(VarPerm(i)) = initbedrock(Model,i,x)
        PRINT *, "bed:", x, VarValues(VarPerm(i))
     END IF
  END DO
END SUBROUTINE getBedrock

FUNCTION initbedrock(Model, nodenumber, x) RESULT(elevation)
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: x, elevation,z0, slope
    slope = 1.0_dp/40.0
    z0 = -550.0_dp
    elevation = z0 - slope*x
END FUNCTION initbedrock

FUNCTION initsurface(Model, nodenumber, x) RESULT(elevation)
  USE ElementDescription
  USE DefUtils

  IMPLICIT NONE
  INTERFACE
     FUNCTION initbedrock(Model, nodenumber, x) RESULT(elevation)
       USE DefUtils
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL (KIND=dp) :: x, elevation
     END FUNCTION initbedrock
  END INTERFACE

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL (KIND=dp) :: x, elevation
  
  elevation = initbedrock(Model, nodenumber, x) + 680.0_dp
END FUNCTION initsurface
