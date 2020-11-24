FUNCTION MySource( Model, n, tx ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: tx(2), f

  !PRINT *,'tx:',tx
  
  f = SIN(tx(1)+tx(2))
    
END FUNCTION MySource

