FUNCTION MyWeertman( Model, n, v ) RESULT(wcoeff)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: v, wcoeff

  REAL(KIND=dp) :: wexp=1.0/3.0_dp, beta=0.02


  
  wcoeff = MIN(beta * v**(wexp-1.0_dp),1.0e20)

  !PRINT *,'w:',wcoeff,v,wexp,beta
  
END FUNCTION MyWeertman

