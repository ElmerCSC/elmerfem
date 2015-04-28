!!!!!!!!!!!!!!!!! PHASE 1
FUNCTION source( model, n, time ) RESULT(current)
  USE DefUtils
  IMPLICIT None
  TYPE(Model_t) :: model
  TYPE(ValueList_t), POINTER :: BF
  LOGICAL :: Found
  INTEGER :: n
  REAL(KIND=dp) :: time, current, I, freq, coeff, phase
  
  BF => Model % BodyForces(1) % Values
  
  CALL ListPushNameSpace(TRIM(ListGetActiveName())//':')
  I = GetConstReal(BF, 'Current', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': I not found in Body Force 1 section.')
  
  freq = GetConstReal(BF, 'Frequency', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': Frequency not found in Body Force 1 section.')
  
  phase = GetConstReal(BF, 'Phase', Found)
  IF (.NOT. FOUND) CALL FATAL('source', ListGetActiveName()//': phase not found in Body Force 1 section.')
  
  WRITE(Message,*) TRIM(ListGetActiveName())//': I = ', I
  CALL Info('source', Message, Level=5 )
  WRITE(Message,*) TRIM(ListGetActiveName())//': frequency = ', freq
  CALL Info('source', Message, Level=5 )
  WRITE(Message,*) TRIM(ListGetActiveName())//': phase = ', phase
  CALL Info('source', Message, Level=5 )

  CALL ListPopNameSpace()
  
  coeff = 1._dp
  IF (2*pi*freq*time <= pi/2) coeff = 4._dp*freq*time 
 
  current = coeff * sqrt(2._dp)*I*sin(2._dp*pi*freq*time-phase)
 
END FUNCTION source
