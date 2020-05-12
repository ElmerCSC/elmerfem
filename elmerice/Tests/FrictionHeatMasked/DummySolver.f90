
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
  TYPE(ValueList_t), Pointer :: BC
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: VarPerm(:)
  INTEGER :: VarDOFs, i, j, k, N, t
  REAL(KIND=dp), POINTER :: VarValues(:)
  LOGICAL :: GotIt

  PRINT *,"DummySolver"
  PRINT *,"***********************************"

  Var => Solver % Variable
  IF (ASSOCIATED(Var)) THEN
     VarPerm => Var % Perm
     VarDOFs =  Var % DOFs
     VarValues => Var % Values
  ELSE
     CALL FATAL('DummySolver','No Variable associated')
  END IF
  k=0
 ! DO i = 1,Model % NumberOfNodes
 !    DO j=1,VarDOFs
 !       k = k + 1
 !       VarValues(VarDOFs*(VarPerm(i) - 1)+j) = k
 !    END DO
 ! END DO
 ! VarValues = 0.0_dp
  DO t=1, Solver % Mesh % NumberOfBoundaryElements
      ! get element information
      Element => GetBoundaryElement(t)
      IF ( .NOT.ActiveBoundaryElement() ) CYCLE
      BC => GetBC()
      n = GetElementNOFNodes()
      VarValues(VarPerm(Element % NodeIndexes)) = ListGetReal(BC,TRIM(Solver % Variable % Name),n,Element % NodeIndexes,GotIt)
   END DO
  
END SUBROUTINE DummySolver
