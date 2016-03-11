! This Solver initialize the two DG variable (have to be run to time) 
RECURSIVE SUBROUTINE InitializeDGVariable(Model, Solver, Timestep, TransientSimulation)
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
  TYPE(Nodes_t) :: ElementNodes
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, SolverName
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: VarSol
  REAL(KIND=dp), ALLOCATABLE :: Depth(:)

  SolverName = 'InitializeDGVariable'
  Mesh => GetMesh()
  n = Mesh % MaxElementNodes
  ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
           Depth(n))

  SolverParams => GetSolverParams()
  VariableName = GetString(SolverParams,'Initialized Variable Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a)')'Keyword >Initialized Variable Name< not found in Solver section'
     CALL FATAL(SolverName, Message)
  END IF

  VarSol => VariableGet(Solver % Mesh % Variables,VariableName,UnFoundFatal=UnFoundFatal)

  Active = GetNOFActive()
  DO t = 1, Active  
     Element => GetActiveElement( t )
     n = GetElementNOfNodes( Element )
     dummyInt = GetElementDOFs( Indexes )
     CALL GetElementNodes( ElementNodes )
     ! Non DG other variable should be read like this
     ! Don't use their perm and the j indexes as it will be wrong 
     ! when the solver is called for a DG variable
     CALL GetScalarLocalSolution(Depth,'Depth')
     DO i=1,n
        j = Indexes(i) 
        z = ElementNodes % z(i)
        D = Depth(i) 
        VarSol % Values(VarSol % Perm(j)) = 1.0_dp - 0.6_dp*z/100.0 
        write(*,*)t,i,j,k, z+D
        
     END DO
  END DO 

  DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z, Depth)
END SUBROUTINE InitializeDGVariable
