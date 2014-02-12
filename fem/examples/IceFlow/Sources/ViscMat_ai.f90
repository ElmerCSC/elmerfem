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
