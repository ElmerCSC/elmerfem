!/*****************************************************************************/
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
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
!/******************************************************************************
! *
! * Surface Boundary Condition for steady state thermal regime
! * Gilbert, A., Sinisalo, A., Gurung, T. R., Fujita, K., Maharjan, S. B., Sherpa, T. C., & Fukuda, T. (2020). 
! * The influence of water percolation through crevasses on the thermal regime of a Himalayan mountain glacier. 
! * The Cryosphere, 14(4), 1273–1288. https://doi.org/10.5194/tc-14-1273-2020
! *
! ******************************************************************************
! *
! *  Authors: Adrien Gilbert
! *
! *  Original Date: May 2020
! *
! *****************************************************************************/
SUBROUTINE SurfEnthBoundarySolver( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

IMPLICIT NONE
TYPE(Model_t) :: Model
TYPE(Variable_t), POINTER :: Accumulation,Rad_fact_var,SurfGrad1,Surfgrad2
TYPE(Variable_t), POINTER :: MB,Dens,Firn,SE,Depth,Melting,Refreeze,Raining,PotRad
TYPE(Solver_t), POINTER :: Solver
TYPE(Element_t),POINTER :: Element
TYPE(ValueList_t), POINTER :: SolverParams
INTEGER, POINTER :: NodeIndexes(:)

INTEGER :: n,i,j,cont,nb_surf,nb_vert,io,nb_year,nb_day,it
REAL(KIND=dp) :: f, z, deg_pos,accu_ref,a,accu,melt_local,accu_ice, temp_10m,rain,t_simu,deg_jour,T,dt
REAL(KIND=dp) :: z_precip,grad_accu,grad,z_temp,seuil_precip,seuil_fonte,Rad_fact_snow,Rad_fact_ice,Rad_fact
REAL(KIND=dp) :: Pfact,temp_correc,surimposed_ice_fact,firn_param,deg1,deg2,precip_correc,sigma,melt

REAL(KIND=dp) :: Sr,rho_w,rho_ice,L_heat,rho_surf,T0,g1,g2,g3,Mean_Temp_Air,x_output,y_output

REAL(KIND=dp) :: fx,fy,slop,asp,S0,dr,lat,L,term1,term2,term3,tau_r,tau_d,tau_b,srad,sinAlpha,R,M,Is
REAL(KIND=dp) :: Ir,Iday,I0,hsr,hs,cos_i,dS,Idiff,reflec,Norm

real (KIND=dp), dimension(:),allocatable :: zprof_1D,TempAir,Precip,PotRadNodes,DensNodes,FirnNodes
real (KIND=dp),dimension(365) :: TempAirMean, PrecipMean,TempAirMeanTry

character(LEN=MAX_NAME_LEN) :: filename,filename2

logical :: first_time=.true.,TransientSimulation, GotIt, PrecipData, node_output=.false.,GotNode=.false.,Output1D,SigmaOK
logical :: OutputFirn,OutputDens,OutputMB,OutputMelting,OutputAccumulation
logical :: OutputRefreeze,OutputRad_Fact_var,OutputRaining,OutputPotRad

save first_time,nb_surf,nb_vert,TempAir,TempAirMean,nb_day,nb_year,PrecipData,Precip,PrecipMean,PotRadNodes,DensNodes,FirnNodes


!===============================================================================
!Initialization=================================================================

if (first_time) then
	first_time=.false.

	filename= GetString(Model % Constants,'AirTemperatureFile', GotIt)
	
	IF (.NOT.GotIt) THEN
      CALL FATAL('Surface Boundary', 'No file for daily air temperature defined in Constants section (AirTemperatureFile=)')
    END IF
	
	filename2= GetString(Model % Constants,'PrecipFile', PrecipData)
	
	IF (.NOT.PrecipData) THEN
		IF (ParEnv % MyPE<1) then
			CALL WARN('Surface Boundary', 'No file for daily precipitation defined in Constants section (PrecipFile=)')
			print*,'Using constant precipitation rate defined by the "Precip" keyword'
		ENDIF
    END IF
	
	open(1,file=filename,status='old')
	nb_day = 0
	do
		read(1,*,iostat=io)
		nb_day = nb_day + 1
		if (io/=0) exit
	enddo
	close(1)
	
	nb_day=nb_day-1
	
	nb_year=floor(nb_day/365.25)
	
	IF (PrecipData) then
		open(1,file=filename2,status='old')
		cont = 0
		do
			read(1,*,iostat=io)
			cont = cont + 1
			if (io/=0) exit
		enddo
		close(1)
		cont=cont-1
	
		IF (cont.ne.nb_day) then
			print*, 'Lenght precip file = ',cont
			print*, 'Lenght temperature file = ', nb_day
			CALL FATAL('Surface Boundary', 'Precip and temperature files have to be same lenght')
		ENDIF
	
		allocate(Precip(nb_day))
	
		open(1,file=filename2,status='old')
		do i=1,nb_day
			read(1,*) Precip(i)
		enddo
		close(1)
	
	!Mean annual precipitation cycle for steady state thermal regime
		PrecipMean(:)=0.0
		cont=0
		do i=1,nb_year
			do j=1,365
				cont=cont+1
				PrecipMean(j)=PrecipMean(j)+Precip(cont)/nb_year
			enddo
		enddo
	
	ENDIF
	
	allocate(TempAir(nb_day))
	TempAirMean(:)=0.0
	TempAir(:)=0.0
		
	IF (ParEnv % MyPE < 1) THEN	
	
	write(*,*) '-------------------------------------------------------------------------------'	
	write(*,*) '...Compute Steady State mass balance and surface boundary for Enthalpy Solver...'
	write(*,*) '-------------------------------------------------------------------------------'
	
	ENDIF
		
	open(1,file=filename,status='old')
	do i=1,nb_day
		read(1,*) TempAir(i)
	enddo
	close(1)
	
	nb_year=floor(nb_day/365.25)
	
!Mean annual temperature cycle for steady state thermal regime
	TempAirMean(:)=0.0
	cont=0
	do i=1,nb_year
		do j=1,365
			cont=cont+1
			TempAirMean(j)=TempAirMean(j)+TempAir(cont)/nb_year
		enddo
	enddo
!------------------------------------------------------------

! Adding random noise that reproduce the correct amount of positive degree day in the mean annual cycle
!----------------------------------------------------------------------------------------------------------------------

! Positive degree day in the data:
cont=0
deg1=0.0
do i=1,nb_year
	do j=1,365
		cont=cont+1
		if (TempAir(cont)>0.0) then
			deg1=deg1+TempAir(cont)/nb_year
		endif
	enddo
enddo

	SigmaOK=.false.
	sigma=0.0
	It=0
	
DO WHILE (.not.SigmaOK)
!Gaussian random nose to take into account daily variability	
	TempAirMeanTry=TempAirMean
	sigma=sigma+0.001
	It=It+1

	call srand(860000)
	do i=1,365
		TempAirMeanTry(i)=TempAirMeanTry(i)+norm_rand(0.0d0,sigma)
	enddo
	
! Positive degree day in the mean annual cycle:
	deg2=0.0
	do i=1,365
		if (TempAirMeanTry(i)>0.0) then
			deg2=deg2+TempAirMeanTry(i)
		endif
	enddo

!Check if agreement between data and mean cylce, if not, increase sigma
	if ((deg2>deg1+0.5).or.(deg2<deg1-0.5)) then
		SigmaOK=.false.
	else
		SigmaOK=.true.
		TempAirMean=TempAirMeanTry
		IF (ParEnv % MyPE < 1) THEN	
			write(*,*) 'Standard Deviation for mean annual cycle =',sigma,' degreeC'
		ENDIF
	endif

	if (It>10000) then
		CALL FATAL('Surface Boundary', 'Mean annual cycle does not match the amount of positive degree day of the data')
	endif

ENDDO


!Get number of surface nodes-----------------------------------
!--------------------------------------------------------------

	Depth => VariableGet( Model % Variables, 'Depth')
	IF ( .not. ASSOCIATED( Depth ) ) THEN
		CALL FATAL('Surface Boundary','Need variable Depth from Flowdepth Solver')
	ENDIF
	

	cont=0

	DO n=1,model % NumberOfNodes
		if (Depth % Values (Depth % perm (n))==0.0) then
			cont=cont+1
		endif
	ENDDO

	nb_surf=cont
	nb_vert=model % NumberOfNodes/cont
	
	
	allocate(FirnNodes(model % NumberOfNodes))
	allocate(PotRadNodes(model % NumberOfNodes))
	allocate(DensNodes(model % NumberOfNodes))

!end first_time
endif

!===============================================================================
!Get parameter and allocate 1D variables========================================

allocate(zprof_1D(nb_vert))

Depth => VariableGet( Model % Variables, 'Depth')
	IF ( .not. ASSOCIATED( Depth ) ) THEN
		CALL FATAL('Surface Boundary','Need variable Depth from Flowdepth Solver')
	ENDIF
	
SurfGrad1 => VariableGet( Model % Variables, 'SurfGrad1')
	IF ( .not. ASSOCIATED( SurfGrad1 ) ) THEN
		CALL FATAL('Surface Boundary','Need variable SurfGrad1 from Flowdepth Solver')
	ENDIF

SurfGrad2 => VariableGet( Model % Variables, 'SurfGrad2')
	IF ( .not. ASSOCIATED( SurfGrad2 ) ) THEN
		CALL FATAL('Surface Boundary','Need variable SurfGrad2 from Flowdepth Solver')
	ENDIF


SE => VariableGet( Model % Variables, 'Surf Enth')

Firn => VariableGet( Model % Variables, 'Firn')
	IF ( .not. ASSOCIATED( Firn ) ) THEN
		OutputFirn=.false.
		ELSE
		OutputFirn=.true.
	ENDIF
	
Dens => VariableGet(Model % Mesh % Variables, 'Densi' )
	IF ( .not. ASSOCIATED( Dens ) ) THEN
		OutputDens=.false.
		ELSE
		OutputDens=.true.
	ENDIF

MB => VariableGet( Model % Variables, 'Mass Balance')
	IF ( .not. ASSOCIATED( MB ) ) THEN
		OutputMB=.false.
		ELSE
		OutputMB=.true.
	ENDIF
	
Melting => VariableGet( Model % Variables, 'Melting')
	IF ( .not. ASSOCIATED( Melting ) ) THEN
		OutputMelting=.false.
		ELSE
		OutputMelting=.true.
	ENDIF

Accumulation => VariableGet( Model % Variables, 'Accu')
	IF ( .not. ASSOCIATED( Accumulation ) ) THEN
		OutputAccumulation=.false.
		ELSE
		OutputAccumulation=.true.
	ENDIF
Refreeze => VariableGet( Model % Variables, 'Refreeze')
	IF ( .not. ASSOCIATED( Refreeze ) ) THEN
		OutputRefreeze=.false.
		ELSE
		OutputRefreeze=.true.
	ENDIF
	
Rad_Fact_var => VariableGet( Model % Variables, 'Rad_Fact')
	IF ( .not. ASSOCIATED( Rad_Fact_var ) ) THEN
		OutputRad_Fact_var=.false.
		ELSE
		OutputRad_Fact_var=.true.
	ENDIF
	
Raining => VariableGet( Model % Variables, 'Rain')
	IF ( .not. ASSOCIATED( Raining ) ) THEN
		OutputRaining=.false.
		ELSE
		OutputRaining=.true.
	ENDIF
	
PotRad => VariableGet( Model % Variables, 'PotRad')
	IF ( .not. ASSOCIATED( PotRad ) ) THEN
		OutputPotRad=.false.
		ELSE
		OutputPotRad=.true.
	ENDIF


z_temp=GetConstReal(Model % Constants, "z_temp")
z_precip=GetConstReal(Model % Constants, "z_precip")
seuil_precip=GetConstReal(Model % Constants, "seuil_precip")
seuil_fonte=GetConstReal(Model % Constants, "seuil_fonte")
surimposed_ice_fact=GetConstReal(Model % Constants, "super_ice")
grad=GetConstReal(Model % Constants, "GradTemp")
Rad_fact_ice=GetConstReal(Model % Constants, "RadFact_ice")
Rad_fact_snow=GetConstReal(Model % Constants, "RadFact_snow")
deg_jour=GetConstReal(Model % Constants, "Deg_jour")
firn_param=GetConstReal(Model % Constants, "firn_param")
grad_accu=GetConstReal(Model % Constants, "GradPrecip")
temp_correc=GetConstReal(Model % Constants, "TempCorrec",GotIt)
	IF (.not.GotIt) THEN
		temp_correc=0.0
	ENDIF


x_output=GetConstReal(Model % Constants, "X_output1D", Output1D)
y_output=GetConstReal(Model % Constants, "Y_output1D", Output1D)


IF (PrecipData) then
	precip_correc=GetConstReal(Model % Constants, "PrecipCorrec",GotIt)
	IF (.not.GotIt) THEN
		precip_correc=1.0
	ENDIF
ELSE
	Pfact=GetConstReal(Model % Constants, "Precip",GotIt)
	IF (.not.GotIt) THEN
		CALL FATAL('Surface Boundary','No precipition file, need to define mean precipition (Precip = )')
	ENDIF
ENDIF

rho_surf=GetConstReal(Model % Constants, "rho_surf")
rho_w=GetConstReal(Model % Constants, "rho_w")
rho_ice=GetConstReal(Model % Constants, "rho_ice")
Sr=GetConstReal(Model % Constants, "Sr")
L_heat=GetConstReal(Model % Constants, "L_heat")
T0=GetConstReal(Model % Constants, "T_ref_enthalpy")



DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
!==============================================================================
!Compute potential solar radiation=================================================
	
S0 = 1367         
dr= 0.0174532925
lat=GetConstReal(Model % Constants, "Latitude")
reflec=0.6

fx=SurfGrad1 % Values (SurfGrad1 % perm (n))
fy=SurfGrad2 % Values (SurfGrad2 % perm (n))

slop=atan(sqrt(fx**2+fy**2))
asp=atan2(fx,fy)*(-1)
L=lat*dr

term1 = sin(L)*cos(Slop) - cos(L)*sin(Slop)*cos(Asp)
term2 = cos(L)*cos(Slop) + sin(L)*sin(Slop)*cos(Asp)
term3 = sin(Slop)*sin(Asp)

srad=0.0
do i=1,365 
    
I0 = S0 * (1.0 + 0.0344*cos(360.0*dr*(194.0+i)/365.0))
dS = 23.45 * dr* sin(360.0*dr * ( (284.0+i)/365.0 ) )
hsr = acos(-tan(L)*tan(dS))
It=nint(12.0*(1.0+hsr/Pi)-12.0*(1.0-hsr/Pi))
Iday=0
     do j=1,It 
             
        hs=hsr-15.0*dr*j     
        sinAlpha = sin(L)*sin(dS)+cos(L)*cos(dS)*cos(hs)
        M=sqrt(1229.0+((614.0*sinAlpha))**2)-614.0*sinAlpha
        tau_b = 0.56 * (exp(-0.65*M) + exp(-0.095*M))
        tau_d = 0.271-0.294*tau_b
        tau_r = 0.271+0.706*tau_b
        cos_i = (sin(dS)*term1) + (cos(dS)*cos(hs)*term2) + (cos(dS)*term3*sin(hs))
        Is = I0 * tau_b
        R = Is * cos_i
		if (R<0.0) then
			R=0
		endif
        Idiff = I0 * tau_d * cos(Slop)*cos(Slop)/2.0 * sinAlpha
        Ir = I0 * reflec * tau_r * sin(Slop)*sin(Slop)/2.0 * sinAlpha
        R= R + Idiff + Ir
		if (R<0.0) then
			R=0
		endif
         Iday=Iday+R
      enddo
srad = srad + Iday/365.0/24.0
enddo  

PotRadNodes(n)=srad

IF (OutputPotRad) THEN	
PotRad % Values (PotRad % perm (n)) = srad
ENDIF

ENDIF
ENDDO


!===============================================================================
!Run Mass Balance model=========================================================

	DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
	
	
		z= model % nodes % z(n)
		accu_ref=0.0
		deg_pos=0.0
		rain=0.0

!===============================================================================
!Get mean steady climate forcing================================================		
		Mean_Temp_Air=0.0
		melt_local=0.0
		
		DO i=1,nb_day
			
			if (PrecipData) then
			Pfact=Precip(i)*365.25*precip_correc
			endif
		
			T=TempAir(i)+grad*(z_temp-z)+temp_correc
			Mean_Temp_Air=Mean_Temp_Air+T/nb_day
    
			melt=(T-seuil_fonte)*deg_jour+rad_fact_snow*PotRadNodes(n)
			if (melt>0) then
				melt_local=melt_local+melt
			endif
	
			if ((T<=seuil_precip).and.(Pfact>1.0e-6)) then
				accu_ref=accu_ref+Pfact/365.25*(1.0+(z-z_precip)*grad_accu)
			endif

			if ((T>seuil_precip).and.(Pfact>1.0e-6)) then
				rain=rain+Pfact/365.25*(1.0+(z-z_precip)*grad_accu)
			endif

		ENDDO


melt_local=melt_local/nb_year
accu_ref=accu_ref/nb_year
rain=rain/nb_year

!===============================================================================
!Compute MB (Gilbert et al., 2016)==============================================

accu=accu_ref
accu_ice = min(surimposed_ice_fact*accu,melt_local)

a=-melt_local-accu_ice*(1.0+1.0/((1.0-rho_surf/rho_ice)*rho_w/rho_surf))+accu

if (a>=0.0) then
	Rad_fact=Rad_fact_snow
else
	Rad_fact=Rad_fact_ice-(Rad_fact_ice-Rad_fact_snow)*&
	&(accu-accu_ice*(1.0+1.0/((1.0-rho_surf/rho_ice)*rho_w/rho_surf)))/melt_local
endif

melt_local=0.0
		
DO i=1,nb_day
	T=TempAir(i)+grad*(z_temp-z)+temp_correc
	melt=(T-seuil_fonte)*deg_jour+rad_fact*PotRadNodes(n)
	if (melt>0) then
			melt_local=melt_local+melt
	endif
ENDDO

melt_local=melt_local/nb_year

if (OutputMB) then
MB % values (MB % perm(n)) = (accu_ice+accu-melt_local)/(rho_ice/rho_w)
endif
if (OutputMelting) then
Melting % values (Melting % perm(n)) = melt_local
endif
if (OutputRaining) then
Raining % values (Raining % perm(n)) = rain
endif
if (OutputAccumulation) then
Accumulation % values (Accumulation % perm(n)) = accu
endif
if (OutputRad_Fact_var) then
Rad_fact_var % values (Rad_fact_var % perm(n)) = Rad_fact
endif

!===============================================================================
!Compute Firn Thickness=========================================================

FirnNodes(n) = a*firn_param !Firn % values (Firn % perm(n)) + a*dt - Firn % values (Firn % perm(n))*dt/firn_param

if (FirnNodes(n)<0.0) then
	FirnNodes(n)= 0.0
endif

if (OutputFirn) then
Firn % values (Firn % perm(n)) = FirnNodes(n)
endif

!===============================================================================
!Compute Density profile and get 1D vertical profiles===========================

do i=1,nb_vert
	cont=n-(i-1)*nb_surf

	if (FirnNodes(n)>1.0) then
		DensNodes(cont)=rho_surf+Depth % values (Depth % perm(cont))&
		&/(FirnNodes(n)*2.0*rho_w)*(rho_ice**2-rho_surf**2)
	else
		DensNodes(cont)=rho_ice
	endif

	if (DensNodes(cont)>rho_ice) then
		DensNodes(cont)=rho_ice
	endif

	zprof_1D(i)= model % nodes % z(cont)
	
	if (OutputDens) then
		Dens % Values(Dens % Perm(cont)) = DensNodes(cont)
	endif
enddo



!Check if export 1D output for this node ---------------------------------
!-------------------------------------------------------------------------
if ((.not.GotNode).and.(Output1D)) then
	if ((model % nodes % x(n)>x_output-100).and.(model % nodes % x(n)<x_output+100)) then
		if ((model % nodes % y(n)>y_output-100).and.(model % nodes % y(n)<y_output+100)) then
			node_output=.true.
			GotNode=.true.
		else
			node_output=.false.
		endif
	else
		node_output=.false.
	endif
else
	node_output=.false.
endif
!------------------------------------------------------------------------- 

!===============================================================================
!Get temperature bellow active layer for diricklet surface boundary of enthalpy
   

 call SolveTemp_1D(Model,Solver,Element,zprof_1D,FirnNodes(n),TempAirMean,temp_10m,n)

 SE % values (SE % perm(n)) = 3.626*temp_10m**2+146.3*temp_10m-T0*(146.3+3.626*T0)

 
endif ! If surface node
ENDDO ! On nodes

  CALL DefaultInitialize()
  
CONTAINS

!===============================================================================
!Solve 1D temperature evolution on vertical profile with CrankNicholson sheme===


subroutine SolveTemp_1D(Model,Solver,Element,z_prof_elmer,firn_thick,TempSurf,temp_10m,surf_node_nb)
  implicit none
  
  TYPE(Model_t) :: Model
  TYPE(Solver_t), POINTER :: Solver
  TYPE(Element_t),POINTER :: Element
  
  real(KIND=dp), dimension(:),allocatable       :: diag,lowerdiag,upperdiag,timeterm,beta,b,ldiag,refreezing,old_refreezing
  real(KIND=dp), dimension(:),allocatable       :: temp,dens,dens_ref,old_temp,water,old_water,vit
  real(KIND=dp), dimension(:),intent(in)        :: z_prof_elmer, TempSurf
  real(KIND=dp), intent(out)                    :: temp_10m
  real(KIND=dp), intent(in)                     :: firn_thick
  integer, intent(in)        					:: surf_node_nb
  
  real (KIND=dp)                               :: prevdxks,newdxks,sumbeta,beta1,beta2,ks,rhos,cont,tsurf,z,thick,z2,dz_ref,z3
  real(KIND=dp)								   :: neige,fonte,dz_ini,total_refreeze,old_temp10m,thick1D,Tair_mean
  real(KIND=dp)								   :: dz,dz_elmer,w1,dt,cpice,melt,accu,T,thick_ref,old_dz,old_zsurf,zsurf
  integer                                      :: n1D,i,j,n_elmer,ii,code,time,cont2,ii_ref,day
  logical									   :: converge

  character*100 :: filename
  
dt=3600*24
cpice=2050.0

n_elmer=size(z_prof_elmer)
thick = abs(z_prof_elmer(n_elmer)-z_prof_elmer(1))
dz_elmer=thick/(n_elmer-1)
dz_ini=0.06
thick1D=10.0

if (thick<thick1D) then
	n1D=floor(thick/dz_ini)+1
	dz=thick/(n1D-1)
else
	n1D=floor(thick1D/dz_ini)+1
    dz=thick1D/(n1D-1)
	thick=thick1D
endif


thick_ref=thick
dz_ref=dz
  
if (.not.allocated(diag)) then
     allocate(diag(n1D))
     allocate(upperdiag(n1D)) 
     allocate(lowerdiag(n1D)) 
     allocate(timeterm(n1D)) 
     allocate(beta(n1D))
	 allocate(b(n1D))
	 allocate(ldiag(n1D))
	 allocate(temp(n1D))
	 allocate(dens(n1D))
	 allocate(water(n1D))
	 allocate(old_water(n1D))
	 allocate(old_temp(n1D))
	 allocate(dens_ref(n1D))
	 allocate(refreezing(n1D))
	 allocate(old_refreezing(n1D))
	 allocate(vit(n1D))
endif
  

  
!===============================================================================
!Initilize density============================================================== 

z=0.0
do i=1,n1D
	if (firn_thick>0.0) then
		dens_ref(i)=rho_surf+z/(firn_thick*2.0*rho_w)*(rho_ice**2-rho_surf**2)
	else
		dens_ref(i)=rho_ice
	endif
	if (dens_ref(i)>rho_ice) then
		dens_ref(i)=rho_ice
	endif
	z=z+dz
enddo

dens=dens_ref
water=0.0
old_water=0.0

!===============================================================================
!Initilize firn/snow thickness and temp profile================================= 

  neige=0.0
  Tair_mean=0.0
  
do day=1,365
	
	    if (PrecipData) then
		Pfact=PrecipMean(day)*365.25*precip_correc
		endif
	
		T=TempSurf(day)+grad*(z_temp-z_prof_elmer(1))+temp_correc
        Tair_mean=Tair_mean+T/365.0
		
		melt=(T-seuil_fonte)*deg_jour+rad_fact_snow*PotRadNodes(surf_node_nb)
		if (melt<0) then
			melt=0.0
		endif
		if ((T<seuil_precip).and.(Pfact>1.0e-6)) then
			accu=Pfact/365.25*(1.0+(z_prof_elmer(1)-z_precip)*grad_accu)
		endif
		
		neige=neige+accu-melt
		
		if (neige<0.0) then
			neige=0.0
		endif
enddo

temp=Tair_mean+273.15

neige=neige*firn_param

temp_10m=0.0
converge=.false.

 if (node_output) then
 open(1,file='Output1D.dat')
 endif


!===============================================================================
!Time loop =====================================================================

  do while (.not.converge) 
  
  old_temp10m=temp_10m
  cont2=0
  temp_10m=0.0
  fonte=0.0
  refreezing=0.0
  total_refreeze=0.0

  
	do day=1,365
	
!===============================================================================
!Compute melt, surface temperature and snow thickness===========================
  
  
    T=TempSurf(day)+grad*(z_temp-z_prof_elmer(1))+temp_correc
	if (PrecipData) then
	Pfact=PrecipMean(day)*365.25*precip_correc
	endif
	
    melt=(T-seuil_fonte)*deg_jour+rad_fact_snow*PotRadNodes(surf_node_nb)
	
	if (melt<0) then
		melt=0.0
	endif


    if ((T<=seuil_precip).and.(Pfact>1.0e-6)) then
       	accu=Pfact/365.25*(1.0+(z_prof_elmer(1)-z_precip)*grad_accu)
    endif
	if ((T>seuil_precip).and.(Pfact>1.0e-6)) then
		rain=Pfact/365.25*(1.0+(z_prof_elmer(1)-z_precip)*grad_accu)
	else
		rain=0.0
	endif
	
	neige=neige+accu-melt-neige/firn_param/365.25
	
	if (neige<0.0) then
		neige=0.0
	endif
	tsurf=T+273.15
	if (tsurf>273.15) then
	tsurf=273.15
	endif

	fonte=fonte+melt+rain
	
!===============================================================================
!Thickness change for seasonal snow=============================================

 old_zsurf = z_prof_elmer(1)+(thick-thick_ref)
 old_dz = dz
 old_temp = temp
 old_water = water
 old_refreezing = refreezing
 
 thick=thick+(accu-melt-neige/firn_param/365.25)/(rho_surf/rho_w)
 
 if (thick<thick_ref) then
	thick=thick_ref
 endif
 
dz=thick/(n1D-1)
zsurf=z_prof_elmer(1)+(thick-thick_ref)

z=zsurf

do i=1,n1D

!Interpolation from old grid to new grid------------------------------------

	if (old_zsurf>=z) then	

		ii=floor((old_zsurf-z)/old_dz)+1
		z2=old_zsurf-(ii-1)*old_dz

		if ((z2-z)<0.0) then
			w1=1.0-(z-z2)/old_dz
		else
			w1=(z2-z)/old_dz
		endif
	
		if (ii>=n1D) then
			ii=n1D-1
		endif
	
		temp(i)=old_temp(ii)*(1.0-w1)+w1*old_temp(ii+1)
		water(i)=old_water(ii)*(1.0-w1)+w1*old_water(ii+1)
		refreezing(i)=old_refreezing(ii)*(1.0-w1)+w1*old_refreezing(ii+1)
	else
		temp(i)=tsurf
		water(i)=0.0
		refreezing(i)=0.0
	endif
	
!Interpolation from elmer to new grid-----------------------------------

	if (z_prof_elmer(1)>=z) then	
		
		ii_ref=floor((z_prof_elmer(1)-z)/dz_ref)+1
		z3=z_prof_elmer(1)-(ii_ref-1)*dz_ref
	
		if ((z3-z)>0.0) then
			w1=(z3-z)/dz_ref
		else
			w1=1.0-(z-z3)/dz_ref
		endif
		
		if (ii_ref>=n1D) then
			ii_ref=n1D-1
		endif
	
		dens(i)=dens_ref(ii_ref)*(1.0-w1)+w1*dens_ref(ii_ref+1)

	else
		dens(i)=rho_surf

	endif
	
    z=z-dz
	
enddo  

!---------------------------------------------------
!Output---------------------------------------------
 if (node_output) then
 if (mod(day,5)==1) then
    z=zsurf
    do i=1,n1D
    write(1,*) temp(i),z,dens(i),water(i),refreezing(i),n1D
    z=z-dz
    enddo
 endif
 endif
!--------------------------------------------------------
!--------------------------------------------------------

!===============================================================================
!Solve 1D diffusion/advection + percolation/refreezing==========================

!Desactivate advection:
vit=0.0

!---------------------------------------------------
!Diffusion/advection--------------------------------

     prevdxks=0.0
     z=0.0
     do i=1,n1D
        rhos=dens(i)
		ks=0.0000025*rhos**2-0.000123*rhos+0.024
		newdxks=dz/ks
        beta(i)=1/(prevdxks+newdxks)
        prevdxks=newdxks
        z=z+dz
     enddo
	 
     beta(1)=0.0
     z=dz
	
! intermediate layers-------
     do i=2,n1D-1
        timeterm(i)=dz*dens(i) * cpice / dt
        beta1=beta(i)
        beta2=beta(i+1)

        sumbeta= beta1+beta2
        diag(i)=timeterm(i) + sumbeta - dens(i)*cpice*vit(i)*(1.0-z/thick_ref)/(dens(i)/rho_w)
        upperdiag(i)=-beta2 + dens(i)*cpice*vit(i)*(1.0-z/thick_ref)/(dens(i)/rho_w)
        lowerdiag(i)=-beta1

        z=z+dz
     end do

! top layer-------
     timeterm(1)=dz * dens(1) * cpice / dt
     diag(1)=1.0
     upperdiag(1)=0.0
     lowerdiag(1)=0.0

! bottom layer-------
     ! condition de Neuman
     timeterm(n1D)=dz* dens(n1D) * cpice / dt
     diag(n1D)=timeterm(n1D) + beta(n1D)
     upperdiag(n1D)=0.0
     lowerdiag(n1D)=-beta(n1D)
 
!Force vector-------------
!-------------------------

!intermediate layers:
  z=dz
  do i=2,n1D-1
     b(i)=(timeterm(i) - beta(i)-beta(i+1)+dens(i)*cpice*vit(i)*(1.0-z/thick_ref)/(dens(i)/rho_w))*temp(i)&
& + beta(i)*temp(i-1) + (beta(i+1)-dens(i)*cpice*vit(i)*(1.0-z/thick_ref)/(dens(i)/rho_w))*temp(i+1)
     z=z+dz
  end do
  
! top layer:
  b(1)=tsurf
  ldiag=diag

! bottome layer :
! condition de Neuman
b(n1D)=(timeterm(n1D) - beta(n1D))*temp(n1D) + beta(n1D)*temp(n1D-1)

!Solve Matrix--------------
!--------------------------
   call tridag(lowerdiag,ldiag,upperdiag,b,temp,n1D,code)
 
 
 
!--------------------------------------------------------
! Percolation/Refreezing --------------------------------
  
 do i=1,n1D
	temp(i)=temp(i)+water(i)*L_heat/(cpice*dens(i))
	if (temp(i)>273.15) then
	
	    refreezing(i)=refreezing(i)+(water(i)-(temp(i)-273.15)*(cpice*dens(i))/L_heat)*dz/rho_w
		total_refreeze=total_refreeze+(water(i)-(temp(i)-273.15)*(cpice*dens(i))/L_heat)*dz/rho_w
		
		water(i)=(temp(i)-273.15)*(cpice*dens(i))/L_heat
		temp(i)=273.15
	else
	    refreezing(i)=refreezing(i)+water(i)*dz/rho_w
		total_refreeze=total_refreeze+water(i)*dz/rho_w
		
		water(i)=0.0
	endif
enddo
 
  cont=(melt+rain)*L_heat*rho_w
     j=1

     do while (cont.ne.0.0)
		if ((dens(j)>800.0)) then
			cont=0.0
		endif

		temp(j)=temp(j)+cont/(cpice*dens(j)*dz)

        if (temp(j)>273.15) then
		
		    refreezing(j)=refreezing(j)+(cont-(temp(j)-273.15)*(cpice*dens(j)*dz))/(L_heat*rho_w)
			total_refreeze=total_refreeze+(cont-(temp(j)-273.15)*(cpice*dens(j)*dz))/(L_heat*rho_w)
			
			cont=(temp(j)-273.15)*(cpice*dens(j)*dz)
							
			
		if (water(j)<Sr*rho_w*(1.0-dens(j)/rho_ice)) then
				
				
				water(j)=water(j)+cont/dz/L_heat
			
				if (water(j)>=Sr*rho_w*(1.0-dens(j)/rho_ice)) then
				
					cont=(water(j)-Sr*rho_w*(1.0-dens(j)/rho_ice))*dz*L_heat
					water(j)=Sr*rho_w*(1.0-dens(j)/rho_ice)
					
				else
					cont=0.0
					
				endif
		endif
						
            temp(j)=273.15
        else
            refreezing(j)=refreezing(j)+cont/(L_heat*rho_w)
			total_refreeze=total_refreeze+cont/(L_heat*rho_w)
			cont=0.0
			
        endif
     

         j=j+1
		 
		 if (j==n1D) then
			cont=0.0
		 endif
	 
     enddo

!---------------------------------------------------------------
!---------------------------------------------------------------


	temp_10m=temp(n1D)+temp_10m
	cont2=cont2+1
 
enddo

!Annual mean at bottom
!---------------------------------------------------------------
temp_10m=temp_10m/cont2

!Convergence criteria-------------------------------------------

if (abs(temp_10m-old_temp10m)<0.005) then
converge=.true.
endif

enddo

if (OutPutRefreeze) then
Refreeze % values (Refreeze % perm(n)) = total_refreeze
endif
  
 end subroutine SolveTemp_1D
 
 
 
!===============================================================================
!===============================================================================
  !*****************************************************************
  ! Solves for a vector U of length N the tridiagonal linear set
  ! M U = R, where A, B and C are the three main diagonals of matrix
  ! M(N,N), the other terms are 0. R is the right side vector.
  !*****************************************************************
  
SUBROUTINE TRIDAG(A,B,C,R,U,N,CODE)
  implicit none
  INTEGER       :: N,J
  REAL  (KIND=dp)       :: BET,GAM(N),A(N),B(N),C(N),R(N),U(N)
  INTEGER       :: CODE

  IF(B(1).EQ.0.D0) THEN
     CODE=1
     RETURN
  END IF

  BET=B(1)
  U(1)=R(1)/BET
  DO J=2,N                    !Decomposition and forward substitution
    GAM(J)=C(J-1)/BET
    BET=B(J)-A(J)*GAM(J)
    IF(BET.EQ.0.D0) THEN            !Algorithm fails
      CODE=2
      RETURN
    END IF
    U(J)=(R(J)-A(J)*U(J-1))/BET
  END DO

  DO J=N-1,1,-1                     !Back substitution
    U(J)=U(J)-GAM(J+1)*U(J+1)
  END DO
  
  CODE=0
  RETURN
END SUBROUTINE TRIDAG

!===============================================================================
!===============================================================================




function norm_rand(mean, std_dev)
    real(KIND=dp) :: norm_rand
    real(KIND=dp), intent(in) :: mean, std_dev
    real(KIND=dp) :: x, y, r
    real(KIND=dp), save :: spare
    logical, save :: has_spare
	

    ! use a spare saved from a previous run if one exists
    if (has_spare) then
        has_spare = .FALSE.
        norm_rand = mean + (std_dev * spare)
        return
    else
        r = 1.0
        do while ( r >= 1.0 )
            ! generate random number pair between 0 and 1
            x=rand()
            y=rand()
            ! normalise random numbers to be in square of side-length = R
            x = (x * 2.0) - 1.0
            y = (y * 2.0) - 1.0
            r = x*x + y*y
        end do

        ! calculate the co-efficient to multiply random numbers x and y
        ! by to achieve normal distribution
        r = sqrt((-2.0 * log(r)) / r)

        norm_rand = mean + (std_dev * x * r)
        spare = y * r
        has_spare = .TRUE.
        return
    end if
end function norm_rand



END SUBROUTINE SurfEnthBoundarySolver



