SUBROUTINE bedrock( Model,Solver,dt,TransientSimulation )

!*************************************************************************
!
! creates a synthetic bedrock variable 
!
!*************************************************************************

  USE DefUtils
  IMPLICIT NONE

!-----------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t), POINTER :: SolverParams

  INTEGER :: ii, tt, nn, jj, DIM
  INTEGER, POINTER :: Permutation(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: x, y, z

  REAL(KIND=dp)  :: fbed


!-----------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  DIM = CoordinateSystemDimension()

  SolverParams => GetSolverParams()

  ! the bedrock is never changed, so it is filled in the firsttime
  DO tt = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(tt)
    nn = GetElementNOFNodes()

    CALL GetElementNodes( Nodes )
    DO ii = 1, nn
      IF ( Permutation(Element % NodeIndexes(ii)) == 0 ) CYCLE
      
      x = Model % Nodes % x(Element % NodeIndexes(ii))

      IF (DIM==3) THEN
            y = Model % Nodes % y(Element % NodeIndexes(ii))
            VariableValues(Permutation(Element % NodeIndexes(ii))) = fbed(x,y) 
      ELSE IF (DIM==2) THEN
            VariableValues(Permutation(Element % NodeIndexes(ii))) = fbed(x,0.0)
      END IF

    END DO

  END DO

END SUBROUTINE bedrock
