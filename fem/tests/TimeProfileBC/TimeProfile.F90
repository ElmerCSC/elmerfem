! Returns a piecewise linear profile such that each control point may be
! a function of time or some other global variable.
!--------------------------------------------------------------------------
FUNCTION TimeProfile( Model, n, tx ) RESULT( f )
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: tx, f

  INTEGER :: tn, tn0 = -1, m, i
  LOGICAL :: Found
  TYPE(ValueListEntry_t), POINTER :: ptr
  TYPE(ValueList_t), POINTER :: BC
  CHARACTER(LEN=MAX_NAME_LEN) :: Name

  SAVE tn0, ptr, Name

  BC => GetBC()
  tn = GetTimestep()

  IF( tn /= tn0 ) THEN
    tn0 = tn

    Name = "Timeprofile"
    
    ptr => ListFind(BC,Name,Found)
    IF(.NOT. Found ) CALL Fatal('TimeProfile','Could not find item: '//TRIM(Name))

    IF( ptr % TYPE /= LIST_TYPE_VARIABLE_SCALAR ) THEN
      CALL Fatal('TimeProfile','Item should be variable scalar: '//TRIM(Name))        
    END IF
    
    IF ( ptr % PROCEDURE /= 0 ) THEN
      CALL Fatal('TimeProfile','Item should not be a function: '//TRIM(Name))             
    END IF
      
    m = SIZE( ptr % Fvalues(1,1,:) )
    
    DO i=1,m
      f = ListGetCReal(BC,TRIM(Name)//' point '//I2S(i),UnfoundFatal=.TRUE.)
      ptr % Fvalues(1,1,i) = f
    END DO

    PRINT *,'Updated temperature profile: ',ptr % Fvalues(1,1,:)
  END IF
    
  f = InterpolateCurve( ptr % TValues,ptr % FValues(1,1,:), &
      tx, ptr % CubicCoeff )

END FUNCTION TimeProfile
  
