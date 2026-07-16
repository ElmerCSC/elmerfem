
!!!!!! ISMIP Test A - u, v, w, p CI selon SIA

FUNCTION uSIAB ( Model, nodenumber, x) RESULT(u)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, y, z, L, u, S, H
   REAL(KIND=dp) ::  tanalpha, ub, beta 
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE L,  tanalpha

   IF (FirstTime) THEN
       FirstTime=.False.
       L = MAXVAL(Model % Nodes % x)
       tanalpha = TAN(0.1_dp*Pi/180.0_dp)
   End If
   
! term (rho.g a / eta_0)^3 /4 = 9.45699e-14 m-3a-1   

     x = Model%Nodes%x(nodenumber)
     y = Model%Nodes%y(nodenumber)
     z = Model%Nodes%z(nodenumber)

     S = -tanalpha*x
     H = 1000.0_dp 

!    ub = (rho.g H sin(alpha) ) / beta^2
     beta = MAX(1.0_dp+SIN(2.0_dp*Pi*x/L) , 0.1_dp) 
     ub = 14.760702_dp / beta 

     u = ub + 9.45699d-14 * ((S-y)**4._dp - H**4._dp)
   
1000 Format(a1)
1002 Format(a1,g14.8)
   
END FUNCTION uSIAB   

FUNCTION wSIAB ( Model, nodenumber, x) RESULT(v)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, i
   REAL(KIND=dp) :: x, y, z, L, v, S, B, H, dHx 
   REAL(KIND=dp) ::  tanalpha 
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE L, tanalpha

   IF (FirstTime) THEN
       FirstTime=.False.
       L = MAXVAL(Model % Nodes % x)
       tanalpha = TAN(0.1_dp*Pi/180.0_dp)
   End If
   
! term (rho.g a / eta_0)^3 /4 = 9.45699e-14 m-3a-1   

     x = Model%Nodes%x(nodenumber)
     y = Model%Nodes%y(nodenumber)
     z = Model%Nodes%z(nodenumber)

     S = -tanalpha*x
     H = 1000.0_dp 
     B = S-H
     dHx = 0.0_dp                                         

     v = -9.45699d-14 * (tanalpha * ( (S-z)**4.0_dp - H**4.0_dp) + 4.0_dp*dHx*(z-B)*H**3.0_dp) 
   
1000 Format(a1)
1002 Format(a1,g14.8)
   
END FUNCTION wSIAB   


FUNCTION pSIAB ( Model, nodenumber, x) RESULT(p)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, i
   REAL(KIND=dp) :: x, y, z, p, S, H 
   REAL(KIND=dp) :: tanalpha

   tanalpha = TAN(0.1_dp*Pi/180.0_dp)

   
! term rho.g  = 0.0089271   

     x = Model%Nodes%x(nodenumber)
     y = Model%Nodes%y(nodenumber)
     z = Model%Nodes%z(nodenumber)

     S = -tanalpha*x

     p = 0.0089271_dp * (S - z) 
   
END FUNCTION pSIAB   

FUNCTION Sliding ( Model, nodenumber, dummyr) RESULT(C)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, i
   REAL(KIND=dp) :: dummyr, x, y, C, L, omega
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime, L, omega

   IF (FirstTime) THEN
     FirstTime=.FALSE.
     L = MAXVAL(Model % Nodes % x)
     omega = 2.0_dp * Pi /L
     WRITE(Message,*) 'Lenght of domain:', L,' Corresponding wavefrequency:',omega 
     CALL INFO('ISMIP_C(Sliding)',Message,Level=1)
   END IF

   x = Model % Nodes % x(nodenumber)
   y = Model % Nodes % y(nodenumber)
   C  = 1.0d-03 *(1.0_dp + SIN(omega* x)*SIN(omega * y))

END FUNCTION Sliding


      Function ViscExp( model ) RESULT(n)

      USE types
      USE DefUtils

      Implicit None
      type(model_t) :: model
      real(kind=dp) :: n
      type(variable_t), pointer :: timevar

      timevar => VariableGet( Model % Variables,'Time')

      IF ( NINT(TimeVar % Values(1)) == 1 ) THEN
          n = 1.0_dp
      ELSE
          n = 1.0_dp/3.0_dp
      END IF
      END FUNCTION ViscExp
                                                     

