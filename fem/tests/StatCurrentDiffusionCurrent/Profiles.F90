! Three-ion manufactured concentration profiles on 0 <= x <= 1 mm.
FUNCTION C1Profile(Model, Node, Time) RESULT(Value)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Time, Value, Xi

  Xi = Model % Nodes % x(Node) / 1.0e-3_dp
  Value = 1000.0_dp + 100.0_dp * Xi * (1.0_dp - Xi)
END FUNCTION C1Profile

FUNCTION C2Profile(Model, Node, Time) RESULT(Value)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Time, Value, Xi

  Xi = Model % Nodes % x(Node) / 1.0e-3_dp
  Value = 800.0_dp - 71.4285714286_dp * Xi * (1.0_dp - Xi)
END FUNCTION C2Profile

FUNCTION C3Profile(Model, Node, Time) RESULT(Value)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Time, Value, Xi

  Xi = Model % Nodes % x(Node) / 1.0e-3_dp
  Value = 1800.0_dp + 28.5714285714_dp * Xi * (1.0_dp - Xi)
END FUNCTION C3Profile
