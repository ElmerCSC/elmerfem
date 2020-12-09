       FUNCTION Thickness(Model,nodenumber,x) RESULT(H)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: x !x
       REAL(kind=dp) :: H
       LOGICAL,SAVE :: FirstTime=.TRUE.
       REAL(kind=dp),SAVE :: Hgl,dhdx,Vgl

       IF (FirstTime) THEN
         CALL RAMP_PARAMETERS(Model,Vgl,Hgl,dhdx)
         FirstTime=.FALSE.
       ENDIF

        H = Hgl + dhdx * x
       
       End FUNCTION Thickness


       FUNCTION Velocity(Model,nodenumber,x) RESULT(U)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: x !x
       REAL(kind=dp) :: U
       LOGICAL,SAVE :: FirstTime=.TRUE.
       REAL(kind=dp),SAVE :: Vgl,Hgl,dhdx
       REAL(kind=dp),SAVE :: A,n,rhoi,rhow,gravity
       REAL(kind=dp) :: astar,Hn
       REAL(kind=dp) :: A_star,H_n

       IF (FirstTime) THEN
        CALL PHYSICAL_PARAMETERS(Model,A,n,rhoi,rhow,gravity)
        CALL RAMP_PARAMETERS(Model,Vgl,Hgl,dhdx)
        FirstTime=.FALSE.
       END IF

       astar=A_star(A,n,rhoi,rhow,gravity)
       Hn=H_n(x,n,Hgl,dhdx)

       U=Vgl+astar*Hn

       END FUNCTION Velocity

       FUNCTION SMB(Model,nodenumber,x) RESULT(adot)
       USE DefUtils
       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp),INTENT(IN) :: x !x
       REAL(kind=dp) :: adot
       LOGICAL,SAVE :: FirstTime=.TRUE.
       REAL(kind=dp),SAVE :: Vgl,Hgl,dhdx
       REAL(kind=dp),SAVE :: A,n,rhoi,rhow,gravity
       REAL(kind=dp) :: astar,H,U
       REAL(kind=dp) :: A_star,Thickness,Velocity

       IF (FirstTime) THEN
        CALL PHYSICAL_PARAMETERS(Model,A,n,rhoi,rhow,gravity)
        CALL RAMP_PARAMETERS(Model,Vgl,Hgl,dhdx)
        FirstTime=.FALSE.
       END IF

       astar=A_star(A,n,rhoi,rhow,gravity)
       H=Thickness(Model,nodenumber,x)
       U=Velocity(Model,nodenumber,x)

       adot=astar*H**(n+1)+dhdx*U

       END FUNCTION SMB


       FUNCTION H_n(x,n,Hgl,dhdx) RESULT(Hn)
       USE DefUtils
       implicit none
       REAL(kind=dp),INTENT(IN) :: x,n,Hgl,dhdx
       REAL(kind=dp) :: Hn

         Hn=-(1._dp/(n+1))*(1._dp/dhdx)*(Hgl**(n+1))*(1._dp-(1._dp+dhdx*x/Hgl)**(n+1)) 

       END FUNCTION H_n
       
       FUNCTION A_star(A,n,rhoi,rhow,gravity) RESULT(Astar)
       USE DefUtils
       implicit none
       REAL(kind=dp),INTENT(IN) :: A,n,rhoi,rhow,gravity
       REAL(kind=dp) :: Astar
       REAL(kind=dp) :: rho_star

        rho_star=rhoi*(rhow-rhoi)/rhow
        Astar=A*(0.25*rho_star*abs(gravity))**n
      
       END FUNCTION A_star 
       


       SUBROUTINE RAMP_PARAMETERS(Model,Vgl,Hgl,dhdx)
       USE DefUtils
       implicit none
       TYPE(Model_t) :: Model
       REAL(kind=dp),INTENT(OUT) :: Vgl,Hgl,dhdx
        
         Hgl = ListGetConstReal( Model % Constants, 'RAMP Hgl', UnFoundFatal=.TRUE. )
         dhdx = ListGetConstReal( Model % Constants, 'RAMP dhdx', UnFoundFatal=.TRUE. )
         Vgl = ListGetConstReal( Model % Constants, 'RAMP Vgl', UnFoundFatal=.TRUE. )

       END SUBROUTINE RAMP_PARAMETERS

       SUBROUTINE PHYSICAL_PARAMETERS(Model,A,n,rhoi,rhow,gravity)
       USE DefUtils
       implicit none
       TYPE(Model_t) :: Model
       REAL(kind=dp),INTENT(OUT) :: A,n,rhoi,rhow,gravity
        A = ListGetConstReal( Model % Constants, 'RAMP RateFactor', UnFoundFatal=.TRUE. )
        n = ListGetConstReal( Model % Constants, 'RAMP Glen', UnFoundFatal=.TRUE. )
        rhoi = ListGetConstReal( Model % Constants, 'RAMP rhoi', UnFoundFatal=.TRUE. )
        rhow = ListGetConstReal( Model % Constants, 'RAMP rhow', UnFoundFatal=.TRUE. )
        gravity = ListGetConstReal( Model % Constants, 'RAMP gravity', UnFoundFatal=.TRUE. )
       END SUBROUTINE PHYSICAL_PARAMETERS
