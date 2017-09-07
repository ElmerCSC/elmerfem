!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> tools for utilizing the GOLF anisotropic law
      !
      ! Definition of the interpolation grid 
      !
  
      Module DefGrid

!     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)  ! If not using
                                                         ! with Elmer   
      USE Types    ! If using with Elmer 
      
      Real, Parameter :: kmin=0.002_dp ! valeur de ki mimum
      Integer, Parameter :: Ndiv=30    ! Ndiv+2 Number of points along ik1
      Integer, Parameter :: Ntot=813   ! Total number of points
      Integer, Parameter :: NetaI=4878 ! 6*4884 length of EtaI
      Integer, Parameter, Dimension(32) :: Nk2 = & 
                (/ -1,  46,  93, 139, 183, 226, 267, 307, 345, 382,&
                 417, 451, 483, 514, 543, 571, 597, 622, 645, 667,& 
                 687, 706, 723, 739, 753, 766, 777, 787, 795, 802,&
                 807, 811/) 
      
      End Module DefGrid
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Non-relative viscosities and Temperature dependency      !!!!!!
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!   Expression of the viscosity matrix in the reference frame    !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  from  the 6 eta6 viscoisties of the matrice law  
!
!  S = eta6(k) tr(M_kd)M_k^D + eta6(k+3)(M_k d + d M_k)^D 
!     
!   get the Voigt eta36  expressed in the general reference frame
!      
!     Si = Aij dj 
!     Aij = eta(i,j) i=1,6; j=1,6 non-symetric matrix 
!      where S=(S11,S22,S33,S12,S23,S31)      
!      and   d=(d11,d22,d33,2 d12,2 d23,2 d31)      
!      Voigt notation as in ELMER
!

       Subroutine ViscGene(eta6,Angle,eta36)
       
        use defgrid
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: Angle    ! Euler Angles        
       Real(kind=dp), Intent(in), Dimension(6) :: Eta6     ! 6 Viscosities of the
                                                  ! matrice law

       Real(kind=dp), Intent(out), Dimension(6,6) :: eta36  ! Viscosity matrix in 
                                                   ! the reference frame
        
       Real(kind=dp) :: p,t,o,st,ct
       Real(kind=dp), Dimension(3,3) :: Q
       Real(kind=dp) :: dQk 
       Real(kind=dp), Dimension(6), Parameter :: coef=(/1.,1.,1.,.0,.0,.0/)
       Integer, Dimension(6), Parameter :: ik=(/1,2,3,1,2,3/)
       Integer, Dimension(6), Parameter :: jk=(/1,2,3,2,3,1/)
       Integer :: k,m,n
       Integer :: i,j 
       

!  Angle = Phi, Theta, Omega
       
       p=Angle(1)
       t=Angle(2)
       o=Angle(3)

! terms of the rotation matrix from RG to RO

        ct = cos(t)
        Q(3,3) = ct
        Q(1,3) = sin(o)
        Q(2,3) = cos(o)
        Q(3,1) = sin(p)
        Q(3,2) = cos(p)

        Q(1,1) = Q(3,2)*Q(2,3) - Q(3,1)*Q(1,3)*ct 
        Q(1,2) = Q(3,1)*Q(2,3) + Q(3,2)*Q(1,3)*ct 
        Q(2,1) = - Q(3,2)*Q(1,3) - Q(3,1)*Q(2,3)*ct 
        Q(2,2) = - Q(3,1)*Q(1,3) + Q(3,2)*Q(2,3)*ct 

        st = sin(t)
        Q(1,3) = Q(1,3)*st
        Q(2,3) = Q(2,3)*st
        Q(3,1) = Q(3,1)*st
        Q(3,2) = -Q(3,2)*st

! 36  terms of the Voigt matrix in the reference frame 
        Do m=1,6
          Do n=1,6
            eta36(m,n) = 0.
          End Do
        End Do
        Do k=1,3
          dQk=Q(k,1)*Q(k,1)+Q(k,2)*Q(k,2)+Q(k,3)*Q(k,3)
          Do m=1,6
            Do n=1,6
              eta36(m,n)=eta36(m,n)+eta6(k)*Q(k,ik(n))*Q(k,jk(n))*&
               (Q(k,ik(m))*Q(k,jk(m))-1./3.*coef(m)*dQk)    
              eta36(m,n) = eta36(m,n) - 2./3.*eta6(k+3)*coef(m)* &
               Q(k,ik(n))*Q(k,jk(n))
            End Do
          End Do

          eta36(1,1)=eta36(1,1)+eta6(k+3)*Q(k,1)*Q(k,1)*2.
          eta36(1,4)=eta36(1,4)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(1,6)=eta36(1,6)+eta6(k+3)*Q(k,1)*Q(k,3)

          eta36(2,2)=eta36(2,2)+eta6(k+3)*Q(k,2)*Q(k,2)*2.
          eta36(2,4)=eta36(2,4)+eta6(k+3)*Q(k,2)*Q(k,1)
          eta36(2,5)=eta36(2,5)+eta6(k+3)*Q(k,2)*Q(k,3)
          
          eta36(3,3)=eta36(3,3)+eta6(k+3)*Q(k,3)*Q(k,3)*2.
          eta36(3,5)=eta36(3,5)+eta6(k+3)*Q(k,3)*Q(k,2)
          eta36(3,6)=eta36(3,6)+eta6(k+3)*Q(k,3)*Q(k,1)

          eta36(4,1)=eta36(4,1)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,2)=eta36(4,2)+eta6(k+3)*Q(k,1)*Q(k,2)
          eta36(4,4)=eta36(4,4)+eta6(k+3)*(dQk-Q(k,3)*Q(k,3))*0.5
          eta36(4,5)=eta36(4,5)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(4,6)=eta36(4,6)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5

          eta36(5,2)=eta36(5,2)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,3)=eta36(5,3)+eta6(k+3)*Q(k,2)*Q(k,3)
          eta36(5,4)=eta36(5,4)+eta6(k+3)*Q(k,1)*Q(k,3)*0.5
          eta36(5,5)=eta36(5,5)+eta6(k+3)*(dQk-Q(k,1)*Q(k,1))*0.5
          eta36(5,6)=eta36(5,6)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5

          eta36(6,1)=eta36(6,1)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,3)=eta36(6,3)+eta6(k+3)*Q(k,1)*Q(k,3)
          eta36(6,4)=eta36(6,4)+eta6(k+3)*Q(k,2)*Q(k,3)*0.5
          eta36(6,5)=eta36(6,5)+eta6(k+3)*Q(k,1)*Q(k,2)*0.5
          eta36(6,6)=eta36(6,6)+eta6(k+3)*(dQk-Q(k,2)*Q(k,2))*0.5
        End Do
         
         End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   sens = 1 : tri des ki pour avoir k1 < k2 < k3
!   sens =-1 : tri des visc apres calcul avec k1 < k2 < k3 
!
!
       Subroutine triki(ki0,ki,visc,ordre,sens)

       use defgrid
       
       implicit none
       
       Real(kind=dp), Dimension(3) ::  ki0,ki
       real(kind=dp), Dimension(6) :: visc,b
       real(kind=dp) :: a
       Integer :: sens,i,j
       Integer, Dimension(3) :: ordre 
       
       

!
!    Passage pour trier les ki
!
       If (sens.EQ.1) Then
         Do i=1,3
           ki(i)=ki0(i)
           ordre(i)=i
         End Do
         Do j=2,3
          a=ki(j)
          Do i=j-1,1,-1
          If (ki(i).LE.a) Goto 20
          ki(i+1)=ki(i)
          ordre(i+1)=ordre(i)
          End Do
  20      Continue
         ki(i+1)=a
         ordre(i+1)=j
         End Do
!
!   Passage pour remettre les viscosite dans le bon ordre
!
       ElseIf (sens.EQ.-1) Then

         Do i=1,6
         b(i)=visc(i)
         End Do
         Do i=1,3
          visc(ordre(i))=b(i)
          visc(ordre(i)+3)=b(i+3)
        End Do

       Else
         Write(*,*)'triki.f : sens <> 1 ou -1' 
       Stop
       End If
       End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Interpolation quadratique d'une quantite Q definie en x1,x2,x3
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      Function InterP(t,x,Q)

      use defgrid
!
      Implicit None
      Real(kind=dp), Dimension(3) :: x,Q
      Real(kind=dp) :: t,InterP,d12,d23
      Real(kind=dp) :: Ip

!
      d12=x(2)-x(1)
      d23=x(3)-x(2)
      Ip=Q(1)*(x(2)-t)*(x(3)-t)/((d12+d23)*d12)
      Ip=Ip+Q(2)*(t-x(1))*(x(3)-t)/(d12*d23)
      Ip=Ip+Q(3)*(t-x(1))*(t-x(2))/((d12+d23)*d23)

      InterP = Ip
      Return
      End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Interpolation quadratique d'une quantite Q definie en 9 points 
!
!         y3 ---     Q7 ---- Q8 ------ Q9
!                    |        |         | 
!                    |        |         |
!         y2 ---     Q4 ----- Q5 ----- Q6
!                    |        |         | 
!                    |        |         |
!         y1 ---     Q1 ----- Q2 ----- Q3
!
!                    |        |         |
!                   x1        x2       x3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      Function InterQ9(x,y,xi,yi,Q)

      use defgrid
!
      Implicit None
      Real(KIND=dp), Dimension(3) ::  xi,yi,a
      Real(kind=dp), Dimension(9) ::  Q
      Real(kind=dp) ::  InterQ9,InterP
      Real(kind=dp) ::  Ip,x,y
      Integer  i

!
         a(1)=InterP(x,xi,Q(1))
         a(2)=InterP(x,xi,Q(4))
         a(3)=InterP(x,xi,Q(7))
         Ip=InterP(y,yi,a)

         InterQ9 = Ip

      Return
      End
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!							              !!!!!!
!!!!!!       Orthotropic Polycrystalline Ice Law from LGGE            !!!!!!
!!!!!! 							              !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   give the Voigt viscosity matrix eta36 expressed 
!                                    in the general reference frame
!     Si = Aij dj 
!     Aij = eta36(i,j) i=1,6; j=1,6 non-symetric matrix 
!      where S=(S11,S22,S33,S12,S23,S31)      
!      and   d=(d11,d22,d33,2d12,2d23,2d31)      
!         as in Elmer
!
       
       Subroutine OPILGGE_ai_nl(ai,Angle,etaI,eta36)
       
       USE DEFGRID
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: ai       ! Texture parameters 
       Real(kind=dp), Dimension(3) :: ki       ! Texture parameters 
       Real(kind=dp), Intent(in), Dimension(3) :: Angle    ! Euler Angles               
       Real(kind=dp), Intent(in), Dimension(:) :: etaI ! Grid Viscosities 
       Real(kind=dp), Intent(out), Dimension(6,6) :: eta36
       Real(kind=dp), Dimension(6) :: eta6
       Real(kind=dp) :: aplusa
       Integer :: i
        
            Do i=1,3
              ki(i)=ai(i)
            End do

            Do i=1,3
             ki(i)=min(Max(ki(i),0._dp),1._dp)
            End do
            aplusa=ki(1)+ki(2)+ki(3)
            If (aplusa.GT.1._dp) then
                    do i=1,3
                     ki(i)=ki(i)/aplusa
                    end do
            endif
       
       ! Value of the 6 relative viscosities in the orthotropic frame

       Call ViscMat_ai(ki,eta6,etaI)

       ! Viscosities in the reference frame

       Call ViscGene(eta6,Angle,eta36)

       End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccc									cccc
!ccccc       subroutine de calcul des viscosites Analytiques            cccc
!ccccc       Les viscosites ont ete calculees sur une grille            cccc
!ccccc       avec le sous-prog makeEtaI.f                               cccc
!ccccc         Elles passent dans etaI(6*814)=EtaI(4884)                cccc
!ccccc 									cccc      
!ccccc      !!! En entree a1,a2,a3 les valeurs propres du tenseur       cccc
!ccccc                                          d'orientation           cccc
!ccccc 									cccc      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
       Subroutine ViscMat_ai(ai0,eta6,etaI)
       
       USE defGrid
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: ai0
       Real(kind=dp), Intent(out), Dimension(6) :: eta6
       Real(kind=dp), Intent(in), Dimension(NetaI) :: etaI
       Real(kind=dp), Dimension(3) :: a1i,a2i
       Real(kind=dp) :: InterQ9
       Real(kind=dp), Dimension(3) :: ai 
       Real(kind=dp), Dimension(9) :: etaN
       Real(kind=dp) :: Delta
       Real(kind=dp) ::  a1,a2
       Real(kind=dp), parameter ::  UnTier = 0.3333333333333333333333333333333333333333_dp  
       Integer, Dimension(3) :: ordre
       Integer :: i,j,n
       Integer :: ik1,ik2
       Integer :: N4,N5,N6
              
       
       Delta = UnTier / Ndiv
!
! tri des ki 
! pour avoir k1 < k2 < k3 
!
!
        Call triki(ai0,ai,eta6,ordre,1)
!
        a1 = ai(1)
        a2 = ai(2)
!
! Position de a1,a2 dans EtaI
! calcul des indices ik1,ik2      
! ik1,ik2 indice du noeud en haut a droite du pt (a1,a2) 
!
         ik1 = Int((a1 + Delta)/Delta) + 1
         ik2 = Int((a2 + Delta)/Delta) + 1

! Si ik1 + 2ik2 -3(Ndiv+1) >0 => on est sur la frontiere a2=a3
! ( a1+2a2=1)
!  => ik1=ik1-1 sauf si ik1=2 , ik2=ik2-1         
!         
         If ((ik1+2*ik2-3*(Ndiv+1)).GE.0) Then
          If ((ik1.NE.2).And.((ik1+2*ik2-3*(Ndiv+1)).NE.0).And.  &
          (abs((a1-Delta*(ik1-1))/a1).GT.1.0E-5))  ik1=ik1-1
          ik2=ik2-1 
         End If
         If (ik1.EQ.1) ik1=2
!
! Indice N4,N5 et N6 des points 4,5 et 6 dans EtaI
!  
         
         N4 = Nk2(ik1-1) + ik2 - ik1 + 3         
         N5 = Nk2(ik1) + ik2 - ik1 + 2         
         N6 = Nk2(ik1+1) + ik2 - ik1 + 1         

!
! Remplissage etaN(1 a 9)
!  7 - 8 - 9  
!  4 - 5 - 6 
!  1 - 2 - 3
!
       Do i=1,3
       a1i(i)=Delta*(ik1-3+i)
       a2i(i)=Delta*(ik2-3+i)
       End Do
!
       Do n=1,6
         Do i=1,3              
           etaN(1+(i-1)*3) = etaI(6*(N4-3+i)+n)
           etaN(2+(i-1)*3) = etaI(6*(N5-3+i)+n)
           etaN(3+(i-1)*3) = etaI(6*(N6-3+i)+n)
         End Do
!
! interpolation sur le Q9  
!
!        If ( (a1 < a1i(1)) .OR. &
!             (a1 > a1i(3)) .OR. &
!             (a2 < a2i(1)) .OR. &
!             (a2 > a2i(3)) ) Then
!           write(*,*)a1,a1i
!           write(*,*)a2,a2i
!           Stop
!         End If
            
         eta6(n)=InterQ9(a1,a2,a1i,a2i,etaN)

       End Do
!
! tri des eta
!
        Call triki(ai0,ai,eta6,ordre,-1)

       Return 
       End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


         End 
!**************************************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute Eigenvalues (ai) and euler angles (Euler)
!  of the second order orientation tensor (a2)
!   using lapack DGEEV lapack routine
!--------------------------------------------------------
      subroutine R2Ro(a2,dim,ai,Euler)

      use Types    ! types d'Elmer

      implicit none

      Real(dp),dimension(6),intent(in) :: a2
      Real(dp),dimension(3),intent(out) :: ai,Euler
      Real(dp),dimension(3,3) :: A,EigenVec
      Real(dp) :: Dumy(1,3),EI(3),Work(24)
      integer :: dim               ! dimension  (2D-3D)
      integer :: i,infor 

      Do i=1,3
         A(i,i)=a2(i)
      End do
         A(1,2)=a2(4)
         A(2,1)=A(1,2)
         A(2,3)=a2(5)
         A(3,2)=A(2,3)
         A(1,3)=a2(6)
         A(3,1)=A(1,3)
      

      CALL DGEEV('N','V',3,A,3,ai,EI,Dumy,1,EigenVec,3,Work,24,infor)

      if (infor.ne.0) &
      CALL FATAL('GolfLaw,R2R0', 'failed to compute fabric eignevalues')

     ! need a right handed orthonormal basis to compute euler angles;
     ! not guarantee by DGEEV.
      EigenVec(1,3)=EigenVec(2,1)*EigenVec(3,2)-EigenVec(3,1)*EigenVec(2,2)
      EigenVec(2,3)=EigenVec(3,1)*EigenVec(1,2)-EigenVec(1,1)*EigenVec(3,2)
      EigenVec(3,3)=EigenVec(1,1)*EigenVec(2,2)-EigenVec(2,1)*EigenVec(1,2)

      Euler(2)=Acos(EigenVec(3,3))
      if (abs(Euler(2)).gt.tiny(Euler(2))) then !3D euler angles 
        Euler(1)=ATAN2(EigenVec(1,3),-EigenVec(2,3))
        Euler(3)=ATAN2(EigenVec(3,1),EigenVec(3,2))
      else ! only one rotation of angle phi
        Euler(3)=0.0
        Euler(1)=ATAN2(EigenVec(2,1),EigenVec(1,1))
      end if

      RETURN
      END
       

