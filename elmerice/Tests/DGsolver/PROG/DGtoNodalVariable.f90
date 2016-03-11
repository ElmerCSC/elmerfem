! This Solver is to write a DG variable into a Nodal Variable 
! There are two solver
! the first write the DG exported variable into a 'true' Nodal variable
! the second write this 'true' nodal variable into the primary DG variable
RECURSIVE SUBROUTINE DGtoNodalVariable1(Model, Solver, Timestep, TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------ internal variables ---------------------------
  REAL(KIND=dp) :: z, D
  INTEGER :: i, j, k, t, n, dummyInt, Active, Indexes(128)
  LOGICAL :: GotIt,UnFoundFatal
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: InName, OutName, SolverName
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: InSol, OutSol

  SolverName = 'DGtoNodalVariable1'
  Mesh => GetMesh()
  n = Mesh % MaxElementNodes

  SolverParams => GetSolverParams()
  InName = GetString(SolverParams,'Input Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a)')'Keyword >Input Variable Name< not found in Solver section'
     CALL FATAL(SolverName, Message)
  END IF

  OutName = GetString(SolverParams,'Output Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a)')'Keyword >Output Variable Name< not found in Solver section'
     CALL FATAL(SolverName, Message)
  END IF

  InSol => VariableGet(Solver % Mesh % Variables,InName,UnFoundFatal=UnFoundFatal)

  OutSol => VariableGet(Solver % Mesh % Variables,OutName,UnFoundFatal=UnFoundFatal)

  Active = GetNOFActive()
  DO t = 1, Active  
     Element => GetActiveElement( t )
     n = GetElementNOfNodes( Element )
     dummyInt = GetElementDOFs( Indexes )
     DO i=1,n
        j = Indexes(i) 
        OutSol % Values(OutSol % Perm(j)) = InSol % Values(InSOl % Perm(j))
     END DO
  END DO 
END SUBROUTINE DGtoNodalVariable1 

RECURSIVE SUBROUTINE DGtoNodalVariable2(Model, Solver, Timestep, TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------ internal variables ---------------------------
  REAL(KIND=dp) :: z, D
  INTEGER :: i, j, k, t, n, dummyInt, Active, Indexes(128)
  LOGICAL :: GotIt,UnFoundFatal
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: InName, OutName, SolverName
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: OutSol
  REAL(KIND=dp), ALLOCATABLE :: InValues(:)

  SolverName = 'DGtoNodalVariable2'
  Mesh => GetMesh()
  n = Mesh % MaxElementNodes
  ALLOCATE(InValues(n))

  SolverParams => GetSolverParams()
  InName = GetString(SolverParams,'Input Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a)')'Keyword >Input Variable Name< not found in Solver section'
     CALL FATAL(SolverName, Message)
  END IF

  OutName = GetString(SolverParams,'Output Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a)')'Keyword >Output Variable Name< not found in Solver section'
     CALL FATAL(SolverName, Message)
  END IF

  OutSol => VariableGet(Solver % Mesh % Variables,OutName,UnFoundFatal=UnFoundFatal)

  Active = GetNOFActive()
  DO t = 1, Active  
     Element => GetActiveElement( t )
     n = GetElementNOfNodes( Element )
     dummyInt = GetElementDOFs( Indexes )
     CALL GetScalarLocalSolution(InValues,TRIM(InName))
     DO i=1,n
        j = Indexes(i) 
        OutSol % Values(OutSol % Perm(j)) = InValues(i)
     END DO
  END DO 
  DEALLOCATE(InValues)
END SUBROUTINE DGtoNodalVariable2 
