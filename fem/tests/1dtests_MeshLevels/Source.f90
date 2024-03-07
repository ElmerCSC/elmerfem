FUNCTION Source( Model, n, f ) RESULT(h)
   USE Lists
   TYPE(Model_t) :: Model
   INTEGER :: n
   REAL(KIND=dp) :: f, h, x, s

   INTEGER :: i, Step, OldStep = -1
   TYPE(Variable_t), POINTER :: T, P

    T => VariableGet( Model % Variables, 'Time' )
    Step = NINT( T % Values(1) )

    P => VariableGet( Model % Variables, 'Potential' )

!   First check the previous result
!   -------------------------------
    IF ( Step>1 .AND. Step /= OldStep ) THEN
       DO i=1,Model % NumberOfNodes
          x = Model % Nodes % x(i)
          s = P % Values( P % Perm(i) )
          SELECT CASE(Step)
             CASE(2) ! can we solve a zero field ?
                IF ( s /= 0 ) &
                   CALL Fatal( 'Source', 'Not able to integrate a zero field ?')
             CASE(3) ! how about a  parabola ?
                IF ( ABS(s+(x**2-PI*x)/2) > 1.0d-8 ) &
                   CALL Fatal( 'Source', 'Not able to integrate parabola ?' )
             CASE(4) ! how about a third degree polynomial
                IF ( ABS(s+(x**3-PI**2*x)/6) > 1.0d-8 ) &
                   CALL Fatal( 'Source', 'Not able to integrate 3d deg polynomial?' )
             CASE(5) ! how about a sine
                IF ( ABS(s-2*sin(x)) > 1.0d-4 ) &
                   CALL Fatal( 'Source', 'Not able to integrate sine?' )
          END SELECT
       END DO
       OldStep = Step
    END IF

!   Then define new one
!   -------------------
    x = Model % Nodes % x(n)

    SELECT CASE(Step)
       CASE(1) ! can we solve a zero field ?
          h = 0
       CASE(2) ! how about a  parabola ?
          h = 1
       CASE(3) ! how about a third degree polynomial
          h = x
       CASE(4) ! how about a sine
          h = f
       CASE(5) ! how about a sine
          h = f
    END SELECT
END FUNCTION Source
