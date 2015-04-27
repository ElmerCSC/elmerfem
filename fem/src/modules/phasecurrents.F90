!!!!!!!!!!!!!!!!! PHASE 1
FUNCTION P1( model, n, time ) RESULT(current)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(ValueList_t), POINTER :: BF
  LOGICAL :: Found
  INTEGER :: n
  REAL(KIND=dp) :: time, current, ISC, freq, coeff
  
  BF => Model % BodyForces(1) % Values
  
  ISC = GetConstReal(BF, 'P1:Current', Found)
  IF (.NOT. FOUND) CALL WARN('P1', 'P1:Current not found in Body Force 1 section.')
  freq = GetConstReal(BF, 'P1:Frequency', Found)
  IF (.NOT. FOUND) CALL WARN('P1', 'P1:Frequency not found in Body Force 1 section.')
  !ISC = 23.09_dp
  !freq = 50._dp 
  !print *, "ISC = ", ISC
  !print *, "freq = ", freq
  !CALL FATAL('P1','Stop this.')

  coeff = 1._dp
  IF (2*pi*freq*time <= pi/2) coeff = 4._dp*freq*time 
 
  current = coeff * sqrt(2._dp)*ISC*sin(2._dp*pi*freq*time)
 
END FUNCTION P1
 
!!!!!!!!!!!!!!!!! PHASE 2
 
FUNCTION P2( model, n, time ) RESULT(current)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(ValueList_t), POINTER :: BF
  LOGICAL :: Found
  INTEGER :: n
  REAL(KIND=dp) :: time, current, ISC, freq, Phase, coeff
  
  BF => Model % BodyForces(1) % Values
  
  ISC = GetConstReal(BF, 'P2:Current', Found)
  IF (.NOT. FOUND) CALL WARN('P2', 'P2:Current not found in Body Force 1 section.')
  freq = GetConstReal(BF, 'P2:Frequency', Found)
  IF (.NOT. FOUND) CALL WARN('P2', 'P2:Frequency not found in Body Force 1 section.')
  Phase = 2._dp*pi/3._dp
  
  !ISC = 23.09_dp
  !freq = 50._dp
  !print *, "ISC = ", ISC
  !print *, "freq = ", freq
  !CALL FATAL('P2','Stop this.')
  
  coeff = 1._dp
  IF (2*pi*freq*time <= pi/2) coeff = 4._dp*freq*time 
 
  current = coeff * sqrt(2._dp)*ISC*sin(2._dp*pi*freq*time-Phase)
 
END FUNCTION P2
 
!!!!!!!!!!!!!!!!! PHASE 3
 
 
FUNCTION P3( model, n, time ) RESULT(current)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(ValueList_t), POINTER :: BF
  LOGICAL :: Found
  INTEGER :: n
  REAL(KIND=dp) :: time, current, ISC, freq, Phase, coeff
  
  BF => Model % BodyForces(1) % Values
  
  ISC = GetConstReal(BF, 'P3:Current', Found)
  IF (.NOT. FOUND) CALL WARN('P3', 'P3:Current not found in Body Force 1 section.')
  freq = GetConstReal(BF, 'P3:Frequency', Found)
  IF (.NOT. FOUND) CALL WARN('P3', 'P3:Frequency not found in Body Force 1 section.')
  Phase = 2*2._dp*pi/3._dp
  
  !ISC = 23.09_dp
  !freq = 50._dp
  !print *, "ISC = ", ISC
  !print *, "freq = ", freq
  !CALL FATAL('P3','Stop this.')
  
  coeff = 1._dp
  IF (2*pi*freq*time <= pi/2) coeff = 4._dp*freq*time 
 
  current = coeff * sqrt(2._dp)*ISC*sin(2._dp*pi*freq*time-Phase)
 
END FUNCTION P3
