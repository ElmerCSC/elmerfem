cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Interpolation quadratique d'une quantite Q definie en x1,x2,x3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      Function InterP(t,x,Q)

      use defgrid
c
      Implicit None
      Real(kind=dp), Dimension(3) :: x,Q
      Real(kind=dp) :: t,InterP,d12,d23
      Real(kind=dp) :: Ip

c
      d12=x(2)-x(1)
      d23=x(3)-x(2)
      Ip=Q(1)*(x(2)-t)*(x(3)-t)/((d12+d23)*d12)
      Ip=Ip+Q(2)*(t-x(1))*(x(3)-t)/(d12*d23)
      Ip=Ip+Q(3)*(t-x(1))*(t-x(2))/((d12+d23)*d23)

      InterP = Ip
      Return
      End
