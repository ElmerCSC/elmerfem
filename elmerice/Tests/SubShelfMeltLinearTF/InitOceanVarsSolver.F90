! InitOceanVarsSolver: sets auxiliary variables needed by SubShelfMeltLinearTF
! from SIF constants, making them available before the melt solver runs.

SUBROUTINE InitOceanVarsSolver( Model, Solver, dt, TransientSimulation )
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t),  INTENT(INOUT) :: Model
  TYPE(Solver_t), TARGET, INTENT(INOUT) :: Solver
  REAL(KIND=dp),  INTENT(IN) :: dt
  LOGICAL,        INTENT(IN) :: TransientSimulation

  TYPE(Variable_t), POINTER :: Zb_var, GM_var, Toce_var, Soce_var, Bedrock_var
  REAL(KIND=dp)              :: Zb_val, GM_val, T_val, S_val, Bedrock_val
  LOGICAL                    :: Found, BedrockFound
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER                    :: ii

  SolverParams => GetSolverParams()

  Zb_val = GetConstReal(SolverParams, 'Zb Value',   Found)
  IF (.NOT. Found) CALL FATAL('InitOceanVarsSolver', 'Zb Value not found')
  GM_val = GetConstReal(SolverParams, 'GroundedMask Value', Found)
  IF (.NOT. Found) CALL FATAL('InitOceanVarsSolver', 'GroundedMask Value not found')
  T_val  = GetConstReal(SolverParams, 'temp_oce_post Value', Found)
  IF (.NOT. Found) CALL FATAL('InitOceanVarsSolver', 'temp_oce_post Value not found')
  S_val  = GetConstReal(SolverParams, 'sal_oce_post Value', Found)
  IF (.NOT. Found) CALL FATAL('InitOceanVarsSolver', 'sal_oce_post Value not found')

  Bedrock_val = GetConstReal(SolverParams, 'bedrock Value', BedrockFound)

  Zb_var    => VariableGet(Solver % Mesh % Variables, 'Zb')
  GM_var    => VariableGet(Solver % Mesh % Variables, 'GroundedMask')
  Toce_var  => VariableGet(Solver % Mesh % Variables, 'temp_oce_post')
  Soce_var  => VariableGet(Solver % Mesh % Variables, 'sal_oce_post')
  Bedrock_var => VariableGet(Solver % Mesh % Variables, 'bedrock')

  IF (.NOT. ASSOCIATED(Zb_var))   CALL FATAL('InitOceanVarsSolver', 'Zb not found')
  IF (.NOT. ASSOCIATED(GM_var))   CALL FATAL('InitOceanVarsSolver', 'GroundedMask not found')
  IF (.NOT. ASSOCIATED(Toce_var)) CALL FATAL('InitOceanVarsSolver', 'temp_oce_post not found')
  IF (.NOT. ASSOCIATED(Soce_var)) CALL FATAL('InitOceanVarsSolver', 'sal_oce_post not found')
  IF (BedrockFound .AND. (.NOT. ASSOCIATED(Bedrock_var))) CALL FATAL('InitOceanVarsSolver', 'bedrock not found')

  Zb_var   % Values = Zb_val
  GM_var   % Values = GM_val
  Toce_var % Values = T_val
  Soce_var % Values = S_val
  IF (BedrockFound) Bedrock_var % Values = Bedrock_val

  CALL INFO('InitOceanVarsSolver','Ocean variables initialised', Level=4)

END SUBROUTINE InitOceanVarsSolver
