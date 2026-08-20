SUBROUTINE CheckDiffusionCurrentMMS(Model, Solver, dt, TransientSimulation)
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  TYPE(Variable_t), POINTER :: Potential, Current
  REAL(KIND=dp), PARAMETER :: Faraday = 96485.33212_dp
  REAL(KIND=dp), PARAMETER :: Conductivity = 19.903500454680898_dp
  REAL(KIND=dp), PARAMETER :: ProfileCoefficient = -8.57142857143e-8_dp
  REAL(KIND=dp), PARAMETER :: PotentialTolerance = 2.0e-9_dp
  REAL(KIND=dp), PARAMETER :: CurrentTolerance = 2.0e-7_dp
  REAL(KIND=dp), PARAMETER :: NegativeControlFloor = 1.0e-5_dp
  REAL(KIND=dp) :: Xi, Shape, ExactPotential
  REAL(KIND=dp) :: PotentialError, CurrentMaximum, NegativeControlError
  INTEGER :: Node, Permutation, Component, Dimension

  Potential => VariableGet(Model % Mesh % Variables, 'Potential')
  Current => VariableGet(Model % Mesh % Variables, 'Volume Current')
  IF (.NOT. ASSOCIATED(Potential)) &
      CALL Fatal('CheckDiffusionCurrentMMS', 'Potential variable not found')
  IF (.NOT. ASSOCIATED(Current)) &
      CALL Fatal('CheckDiffusionCurrentMMS', 'Volume Current variable not found')

  Dimension = CoordinateSystemDimension()
  PotentialError = 0.0_dp
  CurrentMaximum = 0.0_dp
  NegativeControlError = 0.0_dp

  DO Node = 1, Model % NumberOfNodes
    Xi = Model % Nodes % x(Node) / 1.0e-3_dp
    Shape = Xi * (1.0_dp - Xi)
    ExactPotential = -Faraday * ProfileCoefficient * Shape / Conductivity

    Permutation = Potential % Perm(Node)
    IF (Permutation <= 0) &
        CALL Fatal('CheckDiffusionCurrentMMS', 'Potential undefined at a node')
    PotentialError = MAX(PotentialError, &
        ABS(Potential % Values(Permutation) - ExactPotential))

    ! Negative control: the feature-off (zero-potential) field must not pass
    ! the manufactured-potential gate.
    NegativeControlError = MAX(NegativeControlError, ABS(ExactPotential))

    Permutation = Current % Perm(Node)
    IF (Permutation <= 0) &
        CALL Fatal('CheckDiffusionCurrentMMS', 'Volume Current undefined at a node')
    DO Component = 1, Dimension
      CurrentMaximum = MAX(CurrentMaximum, &
          ABS(Current % Values(Current % DOFs * (Permutation - 1) + Component)))
    END DO
  END DO

  WRITE(*,'(A,ES16.8)') 'MMS max potential error       = ', PotentialError
  WRITE(*,'(A,ES16.8)') 'MMS max total current         = ', CurrentMaximum
  WRITE(*,'(A,ES16.8)') 'MMS negative-control error    = ', NegativeControlError

  IF (NegativeControlError <= NegativeControlFloor) CALL Fatal( &
      'CheckDiffusionCurrentMMS', 'Negative control did not fail the gate')
  IF (PotentialError > PotentialTolerance) CALL Fatal( &
      'CheckDiffusionCurrentMMS', 'Manufactured potential is outside tolerance')
  IF (CurrentMaximum > CurrentTolerance) CALL Fatal( &
      'CheckDiffusionCurrentMMS', 'Total current does not cancel')
END SUBROUTINE CheckDiffusionCurrentMMS
