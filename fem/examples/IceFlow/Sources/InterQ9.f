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
