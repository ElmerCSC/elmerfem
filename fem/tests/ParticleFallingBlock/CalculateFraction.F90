!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation,
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!------------------------------------------------------------------------------
SUBROUTINE CalculateFraction( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: VarTot, VarB
  REAL(KIND=dp) :: MinW, MaxW, SumTot, SumB, SumF
  INTEGER :: i, zerocount
  !------------------------------------------------------------------------------

  VarTot => VariableGet( Solver % Mesh % Variables,'Particle Weight')
  VarB => VariableGet( Solver % Mesh % Variables,'Particle Weight B')

  MinW = MINVAL( VarTot % Values )
  MaxW = MAXVAL( VarTot % Values )
  SumTot = SUM( VarTot % Values )
  SumB = SUM( VarB % Values ) 

  PRINT *,'Calculate Fraction'
  PRINT *,'Min and max weight:',minw, maxw
  IF( minw > EPSILON( minw ) ) THEN
    PRINT *,'Ratio of max and min w:',maxw / minw
  END IF
  PRINT *,'Ratio of b and tot w:',SumB / SumTot

  zerocount = 0
  DO i = 1, SIZE( VarTot % Values )
    IF( VarTot % Values(i) < EPSILON( MinW ) ) THEN
      Solver % Variable % Values(i) = -1.0_dp
      zerocount = zerocount + 1
    ELSE  
      Solver % Variable % Values(i) = VarB % Values(i) / VarTot % Values(i)
    END IF
  END DO
  Solver % Variable % Norm = 1.0_dp

  SumF = SUM( Solver % Variable % Values )
  PRINT *,'Sum of fraction:',Sumf
  PRINT *,'Undetermined nodes:',zerocount
  
END SUBROUTINE CalculateFraction
