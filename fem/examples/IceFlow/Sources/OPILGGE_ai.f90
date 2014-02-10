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
       
       Subroutine OPILGGE_ai(ki,Angle,Tc,W,etaI,eta36)
       
       USE DEFGRID
       
       Implicit None
       Real(kind=dp), Intent(in), Dimension(3) :: ki       ! Texture parameters 
       Real(kind=dp), Intent(in), Dimension(3) :: Angle    ! Euler Angles        
       Real(kind=dp), Intent(in) :: Tc                     ! Temperature en d Celsius
       Real(kind=dp), Intent(in), Dimension(7) :: W        ! Glen law parameters
       ! W(1)=B0 Glen's fluidity parameter
       ! W(2)=n Glen's law exponent
       ! W(3)=Q1 Energy activation for Tc<Tl
       ! W(4)=Q2 Energy activation for Tc>Tl
       ! W(5)=T0 temperature reference for B0
       ! W(6)=Tl temperature limite for Q1->Q2
       ! W(7)=R  Gas constant
       Real(kind=dp), Intent(in), Dimension(:) :: etaI ! Grid Viscosities 
       Real(kind=dp), Intent(out), Dimension(6,6) :: eta36
       Real(kind=dp), Dimension(6) :: eta6
       Real(kind=dp) :: Bg,BGlenT
       Integer :: i
        
       
       ! Value of the 6 relative viscosities in the orthotropic frame

       Call ViscMat_ai(ki,eta6,etaI)

       ! Non-relative viscosities and temperature dependency

       Bg=BGlenT(Tc,W)

       Do i=1,6
         eta6(i)=eta6(i)/Bg
       End Do


       ! Viscosities in the reference frame

       Call ViscGene(eta6,Angle,eta36)

       End
