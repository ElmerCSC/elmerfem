!!! calcul de a4 a partir de a2 par fermeture IBOF (Chung 2002)
!! a2 rentre dans l'ordre : 11, 22, 33, 12, 23 ,13
!!!a4  sort dans l'ordre : 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
!
!------------------------------------------------------------------------------
      subroutine IBOF(a2,a4)

      USE Types
       
       implicit none
       Real(dp),dimension(6),intent(in):: a2  
       Real(dp),dimension(9),intent(out):: a4  
       Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
       Real(dp):: b_11,b_22,b_12,b_13,b_23
       Real(dp):: aPlusa

       Real(dp),dimension(21) :: vec
       Real(dp),dimension(3,21) :: Mat
       Real(dp),dimension(6) :: beta
       Real(dp) :: Inv2,Inv3
       integer :: i,j      

       
       a_11=a2(1)
       a_22=a2(2)
       a_33=a2(3)
       a_12=a2(4)
       a_23=a2(5)
       a_13=a2(6)

!       a_11=MIN(MAX(a2(1),0._dp),1._dp)
!       a_22=MIN(MAX(a2(2),0._dp),1._dp)
!       aPlusa=a_11+a_22
!       If(aPlusa.GT.1._dp) then
!               a_11=a_11/aPlusa
!               a_22=a_22/aPlusa
!       EndIf
!       a_33=1._dp-a_11-a_22
!       If(a2(4).GT.0._dp) then
!               a_12=Min(a2(4),0.5_dp)
!       Else
!               a_12=Max(a2(4),-0.5_dp)
!       Endif
!       If(a2(5).GT.0._dp) then
!               a_23=Min(a2(5),0.5_dp)
!       Else
!               a_23=Max(a2(5),-0.5_dp)
!       Endif
!       If(a2(6).GT.0._dp) then
!               a_13=Min(a2(6),0.5_dp)
!       Else
!               a_13=Max(a2(6),-0.5_dp)
!       Endif
       
       

!       write(*,*) 'a2'
!       write(*,*) a2

       !Coefficients 

      Mat(1,1)=0.217774509809788e+02_dp
      Mat(1,2)=-.297570854171128e+03_dp
      Mat(1,3)=0.188686077307885e+04_dp
      Mat(1,4)=-.272941724578513e+03_dp
      Mat(1,5)=0.417148493642195e+03_dp
      Mat(1,6)=0.152038182241196e+04_dp
      Mat(1,7)=-.137643852992708e+04_dp
      Mat(1,8)=-.628895857556395e+03_dp
      Mat(1,9)=-.526081007711996e+04_dp
      Mat(1,10)=-.266096234984017e+03_dp
      Mat(1,11)=-.196278098216953e+04_dp
      Mat(1,12)=-.505266963449819e+03_dp
      Mat(1,13)=-.110483041928547e+03_dp
      Mat(1,14)=0.430488193758786e+04_dp
      Mat(1,15)=-.139197970442470e+02_dp
      Mat(1,16)=-.144351781922013e+04_dp
      Mat(1,17)=-.265701301773249e+03_dp
      Mat(1,18)=-.428821699139210e+02_dp
      Mat(1,19)=-.443236656693991e+01_dp
      Mat(1,20)=0.309742340203200e+04_dp
      Mat(1,21)=0.386473912295113e+00_dp
      Mat(2,1)=-.514850598717222e+00_dp
      Mat(2,2)=0.213316362570669e+02_dp
      Mat(2,3)=-.302865564916568e+03_dp
      Mat(2,4)=-.198569416607029e+02_dp
      Mat(2,5)=-.460306750911640e+02_dp
      Mat(2,6)=0.270825710321281e+01_dp
      Mat(2,7)=0.184510695601404e+03_dp
      Mat(2,8)=0.156537424620061e+03_dp
      Mat(2,9)=0.190613131168980e+04_dp
      Mat(2,10)=0.277006550460850e+03_dp
      Mat(2,11)=-.568117055198608e+02_dp
      Mat(2,12)=0.428921546783467e+03_dp
      Mat(2,13)=0.142494945404341e+03_dp
      Mat(2,14)=-.541945228489881e+04_dp
      Mat(2,15)=0.233351898912768e+02_dp
      Mat(2,16)=0.104183218654671e+04_dp
      Mat(2,17)=0.331489412844667e+03_dp
      Mat(2,18)=0.660002154209991e+02_dp
      Mat(2,19)=0.997500770521877e+01_dp
      Mat(2,20)=0.560508628472486e+04_dp
      Mat(2,21)=0.209909225990756e+01_dp
      Mat(3,1)=0.203814051719994e+02_dp
      Mat(3,2)=-.283958093739548e+03_dp
      Mat(3,3)=0.173908241235198e+04_dp
      Mat(3,4)=-.195566197110461e+03_dp
      Mat(3,5)=-.138012943339611e+03_dp
      Mat(3,6)=0.523629892715050e+03_dp
      Mat(3,7)=0.859266451736379e+03_dp
      Mat(3,8)=-.805606471979730e+02_dp
      Mat(3,9)=-.468711180560599e+04_dp
      Mat(3,10)=0.889580760829066e+01_dp
      Mat(3,11)=-.782994158054881e+02_dp
      Mat(3,12)=-.437214580089117e+02_dp
      Mat(3,13)=0.112996386047623e+01_dp
      Mat(3,14)=0.401746416262936e+04_dp
      Mat(3,15)=0.104927789918320e+01_dp
      Mat(3,16)=-.139340154288711e+03_dp
      Mat(3,17)=-.170995948015951e+02_dp
      Mat(3,18)=0.545784716783902e+00_dp
      Mat(3,19)=0.971126767581517e+00_dp
      Mat(3,20)=0.141909512967882e+04_dp
      Mat(3,21)=0.994142892628410e+00_dp

       
       ! calcul des invariants
       Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
       Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
       
     ! polynome complet de degre 5 des 2 invariants.
         vec(1)=1._dp
         vec(2)=Inv2
         vec(3)=vec(2)*vec(2)
         vec(4)=Inv3
         vec(5)=vec(4)*vec(4)
         vec(6)=vec(2)*vec(4)
         vec(7)=vec(3)*vec(4)
         vec(8)=vec(2)*vec(5)
         vec(9)=vec(2)*vec(3)
         vec(10)=vec(5)*vec(4)
         vec(11)=vec(9)*vec(4)
         vec(12)=vec(3)*vec(5)
         vec(13)=vec(2)*vec(10)
         vec(14)=vec(3)*vec(3)
         vec(15)=vec(5)*vec(5)
         vec(16)=vec(14)*vec(4)
         vec(17)=vec(12)*vec(2)
         vec(18)=vec(12)*vec(4)
         vec(19)=vec(2)*vec(15)
         vec(20)=vec(14)*vec(2)
         vec(21)=vec(15)*vec(4)

       ! calcul des beta_bar (cf annexe C Chung)
       ! attention beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
       !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

       ! calcul des trois beta en fonction du polynome
         beta(:)=0._dp
         Do i=1,3
          Do j=1,21
            beta(i)=beta(i)+Mat(i,j)*vec(j)
          End do
         End do
          
       ! calcul des 3 autres pour avoir la normalisation
         beta(4)=3._dp*(-1._dp/7._dp+beta(1)*(1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
                  beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
                  beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
                  16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

         beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
                  7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
                  beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
                  8._dp*Inv2*Inv2/5._dp))/7._dp

         beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
                   6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

        ! pour avoir les beta_bar
        Do i=1,6
         beta(i)=beta(i)/3._dp
        End do
         beta(2)=beta(2)/2._dp
         beta(5)=beta(5)/2._dp
         beta(6)=beta(6)/2._dp

        !! calcul des 5 b=a.a
        b_11=a_11*a_11+a_12*a_12+a_13*a_13
        b_22=a_22*a_22+a_12*a_12+a_23*a_23
        b_12=a_11*a_12+a_12*a_22+a_13*a_23
        b_13=a_11*a_13+a_12*a_23+a_13*a_33
        b_23=a_12*a_13+a_22*a_23+a_23*a_33

        !Calcul des 9 termes de a4

        a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
         6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
        a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
         6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

        a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
         beta(3)*(b_11*b_22+2._dp*b_12*b_12)


         a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
          beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
          (b_11*b_23+2._dp*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
          beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
          (b_22*b_13+2._dp*b_12*b_23)


         a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
          3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
         a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
         a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
          3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
         a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12


          !modif fab 01/09/04
          ! si a11=0 => on annule les termes de a4
          !IF(a_11.EQ.0_dp) then
                 ! write(*,*) 'a11 = 0'
           !       a4(1)=0._dp
           !       a4(3)=0._dp
           !       a4(4)=0._dp
           !       a4(5)=0._dp
           !       a4(6)=0._dp
           !       a4(7)=0._dp
           !       a4(9)=0._dp
           !Endif
                  
          

         End 
!**************************************************************************************
