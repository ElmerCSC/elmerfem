      MODULE MacAyealContants
        USE DefUtils
        REAL(kind=dp),PARAMETER :: Lx=200.0d03
        REAL(kind=dp),PARAMETER :: Ly=50.0d03
        REAL(kind=dp),PARAMETER :: yearinsec = 31557600.0_dp !365.25*24*60*60
        REAL(kind=dp),PARAMETER :: MPa=1.0d06
      END MODULE MacAyealContants
      
      FUNCTION Zs(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx(2) !x,y
       REAL(kind=dp) :: VarOut

       VarOut=500.0-1.0e-03*tx(1)+20.0*(sin(3.0*pi*tx(1)/Lx)*sin(2.0*pi*tx(2)/Ly))

      END FUNCTION Zs

      FUNCTION H(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx(2) !x,y
       REAL(kind=dp) :: VarOut

        VarOut=1500.0-2.0e-3*tx(1)
      END FUNCTION H

      FUNCTION Zb(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx(2) !x,y
       REAL(kind=dp) :: VarOut
       REAL (KIND=dp) :: Zs,H

        VarOut=Zs(Model,nodenumber,tx)-H(Model,nodenumber,tx)

      END FUNCTION Zb

      FUNCTION betaSquare(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx(2) !x,y
       REAL(kind=dp) :: VarOut

       REAL(kind=dp) :: F1,F2,beta

       F1=sin(3.0*pi*tx(1)/Lx)*sin(pi*tx(2)/Ly)
       F2=sin(pi*tx(1)/(2.0*Lx))*cos(4.0*pi*tx(2)/Ly)
       beta=5.0e3*F1+5.0e03*F2
       VarOut=beta*beta/(MPa*yearinsec)

      END FUNCTION betaSquare

      FUNCTION beta(Model,nodenumber,tx) RESULT(VarOut)
        USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx(2) !x,y
       REAL(kind=dp) :: VarOut
       REAL(kind=dp) :: betaSquare

       VarOut=sqrt(betaSquare(Model,nodenumber,tx))

      END FUNCTION beta 

      FUNCTION INFLOW(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx !y
       REAL(kind=dp) :: VarOut

       VarOut=4.753e-6*yearinsec*(sin(2.0*pi*(Ly-tx)/Ly)+2.5*sin(pi*(Ly-tx)/Ly))
      END FUNCTION INFLOW
      FUNCTION OUTFLOW(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx !y
       REAL(kind=dp) :: VarOut

       VarOut=1.584e-5*yearinsec*(sin(2.0*pi*(Ly-tx)/Ly)+2.5*sin(pi*(Ly-tx)/Ly)+0.5*sin(3.0*pi*(Ly-tx)/Ly))
      END FUNCTION OUTFLOW

      FUNCTION Density(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx !dumy
       REAL(kind=dp) :: VarOut

       VarOut=917.0/(MPa*yearinsec**2)

      END FUNCTION Density
      FUNCTION Gravity(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx !dumy
       REAL(kind=dp) :: VarOut

       VarOut=-9.81*yearinsec**2

      END FUNCTION Gravity

      FUNCTION Viscosity(Model,nodenumber,tx) RESULT(VarOut)
       USE MacAyealContants
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: tx !dumy
       REAL(kind=dp) :: VarOut

       VarOut=1.8e8*1.0e-6*(2.0*yearinsec)**(-1.0/3.0)

      END FUNCTION Viscosity
