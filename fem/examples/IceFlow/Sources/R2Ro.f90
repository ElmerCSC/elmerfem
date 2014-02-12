!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine de calcul des 3 : - valeurs propres de a2
!                              - angles d'Euler
!
!--------------------------------------------------------
      subroutine R2Ro(a2,dim,ai,Euler)

      use Types    ! types d'Elmer

      implicit none

      Real(dp),dimension(6),intent(in) :: a2
      integer :: dim               ! dimension  (2D-3D)
      Real(dp),dimension(3),intent(out) :: ai,Euler
      Real(dp) :: disc
      Real(dp), parameter :: Prec=1.e-06
      Real(dp) :: a,b,c,aplusa
      Complex :: Sol(3)
      Integer :: typeSol
      Integer :: i


      SELECT CASE(dim)

! ******** case 2D **********
      CASE(2)
        
     ! a2 in the orthotropic frame
     !   IF (abs(a2(4)).LT.Prec) Then
     !      ai(1:3)=a2(1:3)
     !      Euler=0._dp
     !      
     !   ELSE
        
           disc = (a2(1)+a2(2))*(a2(1)+a2(2)) - &
              4.*(a2(1)*a2(2)-a2(4)*a2(4))

           IF (disc.LT.0._dp) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)
        
            ai(2) = 0.5 * (a2(1)+a2(2)+disc)
            ai(1) = 0.5 * (a2(1)+a2(2)-disc)

            Euler=0._dp
            Euler(1) = atan2(a2(4),a2(1)-ai(2))


            Do i=1,2
             ai(i)=Max(ai(i),0._dp)
            End do

            
            aplusa=ai(1)+ai(2)
            If (aplusa.GT.1._dp) then
               !  write(*,*) 'depasse 1 dans R2R0',aplusa,ai(1),ai(2)
                    do i=1,2
                     ai(i)=ai(i)/aplusa
                    end do
            endif
            
            ai(3) = 1._dp-ai(1)-ai(2)
        
      !  END IF

! *********** case 3D **************
      CASE(3)
        ! a2 in the orthotropic frame.
        IF ((abs(a2(4)).LT.prec).and.(abs(a2(5)).LT.prec) &
           .and.(abs(a2(6)).LT.prec) ) THEN
              ai(1:3)=a2(1:3)
              Euler=0.
        ! plane (1,2) isotropic      
        ELSE IF ((abs(a2(5)).LT.prec).and.(abs(a2(6)).LT.prec)) THEN

           disc = (a2(1)+a2(2))*(a2(1)+a2(2)) - &
              4.*(a2(1)+a2(2)+a2(4)*a2(4))
              

           IF (disc.LT.0._dp) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)

            
            ! 2 eigenvalues
            ai(2) = 0.5 * (a2(1)+a2(2)+disc)
            ai(1) = 0.5 * (a2(1)+a2(2)-disc)

            ! eigenvalues beetwen 0 and 1

            Do i=1,2
             ai(i)=Max(ai(i),0._dp)
            End do
            aplusa=ai(1)+ai(2)
            If (aplusa.GT.1._dp) then
                    do i=1,2
                     ai(i)=ai(i)/aplusa
                    end do
            endif

            ai(3) = 1._dp-ai(1)-ai(2)

            Euler=0.
            Euler(1) = atan2(a2(4),a2(1)-ai(2))
                
         ! plane (2,3) isotropic
         ELSE IF ((abs(a2(4)).LT.prec).and.(abs(a2(6)).LT.prec)) THEN

           disc = (a2(2)+a2(3))*(a2(2)+a2(3)) - &
              4.*(a2(2)+a2(3)+a2(5)*a2(5))

           IF (disc.LT.0._dp) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)
            
            ai(2) = 0.5 * (a2(1)+a2(3)+disc)
            ai(3) = 0.5 * (a2(1)+a2(3)-disc)
        
            Do i=2,3
             ai(i)=Max(ai(i),0._dp)
            End do
            aplusa=ai(2)+ai(3)
            If (aplusa.GT.1._dp) then
                    do i=2,3
                     ai(i)=ai(i)/aplusa
                    end do
            endif

            ai(1) = 1._dp-ai(2)-ai(3)

            Euler=0.
            Euler(2) = atan2(a2(5),a2(2)-ai(3))

          ! plane (1,3) isotropic
         ELSE IF ((abs(a2(4)).LT.prec).and.(abs(a2(5)).LT.prec)) THEN
         
           disc = (a2(1)+a2(3))*(a2(1)+a2(3)) - &
              4.*(a2(1)+a2(3)+a2(6)*a2(6))

           IF (disc.LT.0._dp) THEN
              print *, 'Pb to calculate eigenvalues of a2'
              stop
           END IF

            disc = sqrt(disc)

            ai(3) = 0.5 * (a2(1)+a2(3)+disc)
            ai(1) = 0.5 * (a2(1)+a2(3)-disc)
        
            Do i=1,3,2
             ai(i)=Max(ai(i),0._dp)
            End do
            aplusa=ai(1)+ai(3)
            If (aplusa.GT.1._dp) then
                    do i=1,3,2
                     ai(i)=ai(i)/aplusa
                    end do
            endif

            ai(2) = 1._dp-ai(1)-ai(3)

            
            Euler(1) = -0.5*Pi
            Euler(2) = atan2(a2(6),a2(1)-ai(3))
            Euler(3)= 0.5*Pi

         ! general case
        ELSE
                
          ! Recherche des 3 valeurs propres ai(i)
          ! equation caracteristique de la forme l^3+al^2+bl+c=0
                
             a=-1.
             b=a2(1)*(a2(2)+a2(3))+a2(2)*a2(3) &
                -a2(5)*a2(5)-a2(4)*a2(4)-a2(6)*a2(6)
             c=-a2(1)*a2(2)*a2(3)-2.*a2(4)*a2(6)*a2(5) &
               +a2(1)*a2(5)*a2(5)+a2(2)*a2(6)*a2(6) &
               +a2(3)*a2(4)*a2(4)

             Call degre3(a,b,c,Sol,typeSol)
              If (typeSol.EQ.1) Then
                Do i=1,3
                 ai(i)=real(Sol(i))
                End Do
              Else
                Write(*,*)'Solution complex in R2Ro'
                stop
              END IF
              
            !! Euler Angles                                                           
                                              
             write(*,*) 'pas encore pret'
             stop
         
        END IF

      CASE DEFAULT
        write(*,*) 'dimens',dim
        print *, 'Eigenvalues of a2, dimension not available'
      
      END SELECT

!      write(*,*)'------------------------------------------------------'
!      write(*,'((a),6(e13.5,x))') 'a2',a2(:)
!      write(*,'((a),3(e13.5,x))') 'ai',ai(:)
!      write(*,'((a),3(e13.5,x))') 'Euler',Euler(:)
!      write(*,*) '-----------------------------------------------------'

      

      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solution equation du troisieme degre
! 
! type=1 3 racines reelles
! type=2 1 racine reelle + 2 complexes
!
!
      Subroutine degre3(a,b,c,sol,type)
      
      Use Types

      Implicit None
      Real(kind=dp) ::   a,b,c,R1,R
      Real(kind=dp) ::  p,q,phi
      Integer  :: i,type
      Complex  :: sol(3)


      p=b-a*a/3.
      q=2.*a*a*a/27.-a*b/3.+c

      R=q*q/4.+p*p*p/27.

      If (R.GT.0.) Then
       R1=(-q/2.+sqrt(R))**(1./3.)+(-q/2.-sqrt(R))**(1./3.)
       sol(1)=cmplx(R1,0.)
       sol(2)=cmplx(-R1/2.,sqrt(3.*R1*R1+4.*p)/2.)
       sol(3)=cmplx(-R1/2.,-sqrt(3.*R1*R1+4.*p)/2.)
       type=2
      Else If (R.LT.0.) Then
        phi=Acos(-q*sqrt(-27./(4.*p*p*p)))
        sol(1)=cmplx(2.*sqrt(-p/3.)*cos(phi/3.),0.)
        sol(2)=cmplx(2.*sqrt(-p/3)*cos((2.*Pi+phi)/3.),0.)
        sol(3)=cmplx(2.*sqrt(-p/3)*cos((4.*Pi+phi)/3.),0.)
        type=1
      Else
         sol(1)=cmplx(3.*q/p,0.)
         sol(2)=cmplx(-3.*q/(2.*p),0.)
         sol(3)=cmplx(-3.*q/(2.*p),0.)
         type=1
      End If
      Do i=1,3
        sol(i)=sol(i)- cmplx(a/3.,0.)
      End Do


      End
        
