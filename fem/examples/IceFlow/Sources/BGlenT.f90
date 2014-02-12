!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Non-relative viscosities and Temperature dependancy      !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       Function BGlenT(Tc,W)
       
       USE defGrid
       
       Implicit None
       Real(kind=dp) :: BGlenT
       Real(kind=dp), Intent(in) :: Tc                     ! Temperature en d Celsius
       Real(kind=dp), Intent(in), Dimension(7) :: W        ! Glen law parameters
       ! W(1)=B0 Glen's fluidity parameter
       ! W(2)=n Glen's law exponent
       ! W(3)=Q1 Energy activation for Tc<Tl
       ! W(4)=Q2 Energy activation for Tc>Tl
       ! W(5)=T0 temperature reference for B0
       ! W(6)=Tl temperature for Q1->Q2
       Real(kind=dp), parameter :: Tzero=273.15
       Real(kind=dp) :: Q, DT
       
       
       Q= W(3)
       If (Tc.GT.W(6)) Q=W(4)
       Q=Q/W(7)
       DT=-1./(Tc+Tzero)+1./(W(5)+Tzero)

       BGlenT=W(1)*Exp(Q*DT)

       End
