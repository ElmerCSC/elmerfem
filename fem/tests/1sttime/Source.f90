FUNCTION Source( Model, n, f ) RESULT(h)
   USE Lists
   TYPE(Model_t) :: Model
   INTEGER :: n
   REAL(KIND=dp) :: f, h, x, s, Time, OldTime = 0

   INTEGER :: i, Case = 1
   TYPE(Variable_t), POINTER :: T, P

    T => VariableGet( Model % Variables, 'Time' )
    Time = T % Values(1)
    P => VariableGet( Model % Variables, 'Potential' )

!   First check the previous result
!   -------------------------------
    IF ( Time > OldTime + 1 + 1.0d-8 ) THEN
       OldTime = FLOOR( Time )
       x = Model % Nodes % x(1)
       s = P % Values( P % Perm(1) )
       SELECT CASE(Case)
          CASE(1) ! can we solve straight line
              IF ( ABS(f-1.0d0) > 1.0d-10 ) &
                CALL Fatal( 'Source', 'Not able to integrate straight line ?' )
              Model % Solver % Matrix % Force = 0.0d0
              P % Values = 0.0d0
              P % PrevValues(:,:) = 0.0d0
          CASE(2) ! how about a  parabola ?
              IF ( ABS(s-0.5d0) > 1.0d-6 ) &
                CALL Fatal( 'Source', 'Not able to integrate parabola ?' )
              Model % Solver % Matrix % Force = 0.0d0
              P % Values = 0.0d0
              P % PrevValues(:,:) = 0.0d0
          CASE(3) ! how about a sine
             IF ( ABS(s-sin(1.0d0)) > 1.0d-3 ) &
                CALL Fatal( 'Source', 'Not able to integrate sine?' )
       END SELECT
       Case = Case + 1
    END IF

    Time = Time - FLOOR(Time)
    IF ( ABS(Time) < 1.0d-8 ) Time = 1

!   Then define new one
!   -------------------
    SELECT CASE(Case)
       CASE(1) ! how about a  parabola ?
          h = 1
       CASE(2) ! how about a third degree polynomial
          h = Time
       CASE(3) ! how about a sine
          h = cos(Time)
       CASE(4) ! how about a sine
          h = 0.0d0
    END SELECT
END FUNCTION Source
