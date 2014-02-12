!
!
!  Function a(D) and b(D) from Gagliardini and Meyssonier, 1997.
!  modified to fulfill the condition 3x a(D) >= 2x b(D) for D > 0.1 
!  for n=3 only 

FUNCTION ParameterA ( Model, nodenumber, D ) RESULT(a)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: a, D, DD     

    IF (D >= 1.0_dp) THEN
      a = 1.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.1 ) DD = 0.1_dp
      a = EXP( 16.02784497_dp - 19.25_dp * DD )

    ELSE
      a =  1.0_dp  + 2.0/3.0 * (1.0_dp - D) 
      a = a / ( D**1.5 )

    End If
END FUNCTION ParameterA 


FUNCTION ParameterB ( Model, nodenumber, D ) RESULT(b)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: b, D, DD 

    IF (D >= 1.0_dp) THEN
      b = 0.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.1 ) DD = 0.1_dp
      b = EXP( 16.75024378_dp - 22.51_dp * DD )

    ELSE
      b = (1.0_dp - D)**(1.0/3.0) 
      b = 3.0/4.0 * ( b / (3.0 * (1 - b)) )**1.5

    End If
END FUNCTION ParameterB 
