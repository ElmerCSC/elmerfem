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
