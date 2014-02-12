SUBROUTINE TrescaSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: Stress, Tresca
  INTEGER, POINTER :: StressPerm(:), TrescaPerm(:)
  INTEGER :: i, j, nodes, StressDOFs, TrescaDOFs, ok
  REAL(KIND=dp), POINTER :: StressValues(:), TrescaValues(:)
  REAL(KIND=dp) :: StressComponent(6), Matrix(3,3), Eigs(3), Work(8)
  REAL(KIND=dp) :: SigmaMax, SigmaMin
  SAVE TrescaValues
!------------------------------------------------------------------------------

  ! Get stress tensor:
  !--------------------
  Stress => VariableGet( Solver % Mesh % Variables, 'Stress' )
  
  IF( ASSOCIATED(Stress) ) THEN
     StressPerm => Stress % Perm
     StressDOFs = Stress % DOFs
     StressValues => Stress % Values
  ELSE
     STOP 'Error: Stress was not found'
  END IF

  IF( StressDOFs /= 6 ) THEN
     STOP 'Error: Unexpected number of DOFs for stress'
  END IF

  nodes = SIZE( StressValues ) / StressDOFs

  ! Introduce Tresca stress:
  !--------------------------
  Tresca => VariableGet( Solver % Mesh % Variables, 'Tresca' )

  IF( .NOT.ASSOCIATED(Tresca) ) THEN
     TrescaPerm => StressPerm
     TrescaDOFs = 1
     ALLOCATE( TrescaValues( nodes ) )
     CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
          Solver, 'Tresca', TrescaDOFs, TrescaValues, TrescaPerm )
  END IF

  ! Compute Tresca stress in node points:
  !--------------------------------------
  DO i = 1, nodes
     DO j = 1, StressDOFs
        StressComponent(j) = StressValues(StressDOFs*(i-1)+j)
     END DO

     Matrix(1,1) = StressComponent(1)
     Matrix(2,2) = StressComponent(2)
     Matrix(3,3) = StressComponent(3)
     Matrix(1,2) = StressComponent(4)
     Matrix(2,3) = StressComponent(5)
     Matrix(1,3) = StressComponent(6)

     CALL DSYEV('N', 'U',  3, Matrix, 3, Eigs, Work, 8, ok)

     IF( ok /= 0 ) THEN
        STOP 'Error: DSYEV failed'
     END IF

     SigmaMax = MAXVAL( Eigs )
     SigmaMin = MINVAL( Eigs )

     TrescaValues(i) = SigmaMax - SigmaMin
  END DO

END SUBROUTINE TrescaSolver
!------------------------------------------------------------------------------
