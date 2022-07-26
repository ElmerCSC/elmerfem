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
! *  Authors: Thomas Zwinger, Denis Cohen, Juha Hartikainen
! *  Email:  thomas Zwinger [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date:  January 2017  -               
! * 
! *****************************************************************************
!>  Solute transport solver
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> solute (salt ions) transport equation for enhanced permafrost model
!> \ingroup Solvers
SUBROUTINE PermafrostSoluteTransport( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material
  TYPE(Variable_t), POINTER :: PressureVar,SalinityVar,PorosityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor, MaxSalinity, MaxSalinityValue, AverageCorrectedMaxValue, &
       MinSalinity, MinSalinityValue, AverageCorrectedMinValue
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE., ElementWiseRockMaterial, ActiveMassMatrix = .TRUE., &
       InitializeSteadyState = .FALSE., CorrectValues=.FALSE., ExtForce=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostSoluteTransport'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, VarName, TemperatureName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName
  TYPE(ValueHandle_t) :: Temperature_h, Pressure_h, Salinity_h, Porosity_h, Load_h

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,&
       CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       ActiveMassMatrix, InitializeSteadyState, &
       CorrectValues, MinSalinity, MaxSalinity
  !------------------------------------------------------------------------------
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, 'Computing solute transport           ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()

  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
  ExtForce = GetLogical(Params,'Compute External Force fc', Found)

  Salinity => Solver % Variable % Values
  IF (.NOT.ASSOCIATED(Salinity)) THEN
    WRITE(Message,*) "Variable for solute fraction not associated"
    CALL FATAL(Solvername,Message)
  END IF
  SalinityPerm => Solver % Variable % Perm

  ! check, whether we assume steady state (despite transient run)
  ! this can come handy to produce a balance-pressure field at the
  ! start of the simulation
  !---------------------------------------------------------------  
  IF (FirstTime) THEN
    InitializeSteadyState = GetLogical(Params,'Initialize Steady State',Found)
  END IF

  IF (InitializeSteadyState) THEN
    IF (GetTimeStep() == 1) THEN
      CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=6)
      ActiveMassMatrix = .FALSE.
    ELSE 
      CALL INFO(SolverName,"Switching mass matrix to active after initializing with steady state",Level=6)
      ActiveMassMatrix = .TRUE.
      InitializeSteadyState = .FALSE.
    END IF
  END IF
  
  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  IF (FirstTime) THEN
    ! handles to local variables
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )
    ! Handle to salinity Source (possible description of heat source at elements/IP's) 
    CALL ListInitElementKeyword( Load_h,'Body Force','Salinity Source' ) 


    CorrectValues = GetLogical(Params, 'Correct Values', Found)
    IF (CorrectValues) THEN
      CALL INFO(SolverName,' "Correct Values" set - negative salinities will be reset to minimum value')
      MinSalinity = GetConstReal(Params, 'Minimum Salinity', Found)
      IF (.NOT.Found) THEN
        CALL WARN(SolverName,' "Minimum Salinity" no found, setting to 1.0E-12')
        MinSalinity = 1.0d-12
      ELSE
        WRITE(Message,*) ' "Minimum Salinity" set to ', MinSalinity
        CALL INFO(SolverName,Message,Level=5)
      END IF
      MaxSalinity = GetConstReal(Params, 'Maximum Salinity', Found)
      IF (.NOT.Found) THEN
        CALL WARN(SolverName,' "Maximum Salinity" no found, setting to 2.0')
        MaxSalinity = 2.0_dp
      ELSE
        WRITE(Message,*) ' "Maximum Salinity" set to ', MaxSalinity
        CALL INFO(SolverName,Message,Level=5)
      END IF
    END IF
  END IF


  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL INFO( SolverName,'Nonlinear iteration '//TRIM(I2S(iter))//&
        ' out of '//TRIM(I2S(maxiter)), Level=4)
    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      Material => GetMaterial()
      IF (FirstTime) THEN        

        ! check, whether we have globally or element-wise defined values of rock-material parameters
        ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
        IF (ElementWiseRockMaterial) THEN
          WRITE (Message,*) 'Found "Element Rock Material File"'
          CALL INFO(SolverName,Message,Level=5)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=5)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material)
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=5)
          FirstTime = .FALSE.
        END IF
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )        
      END IF

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        CALL INFO(SolverName,'Permafrost Phase Change Model" set to '//TRIM(PhaseChangeModel),Level=9)
      END IF

      CALL LocalMatrixSolute(  Element, Element % ElementIndex, Active, n, nd+nb,&
           CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial,ActiveMassMatrix, ExtForce)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()

    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCSolute(  Element, n, nd+nb,&
             CurrentSoluteMaterial, CurrentSolventMaterial,&
             NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)
        ! PRINT *,"Solute:",t,"of",Active,":",n, nb
      END IF
    END DO
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    ! correct values to positive (if requested)
    !------------------------------------------
    IF (CorrectValues) THEN
      J=0
      K=0
      AverageCorrectedMinValue=0.0_dp
      AverageCorrectedMaxValue=0.0_dp
      MinSalinityValue = MinSalinity
      MaxSalinityValue = MaxSalinity
      DO I = 1, Solver % Mesh % NumberOfNodes
        IF (SalinityPerm(I) < 1) CYCLE
        IF (Salinity(SalinityPerm(I)) < MinSalinity) THEN
          AverageCorrectedMinValue=AverageCorrectedMinValue+Salinity(SalinityPerm(I))
          MinSalinityValue = MIN(Salinity(SalinityPerm(I)), MinSalinityValue)
          Salinity(SalinityPerm(I))=MinSalinity          
          J=J+1
        END IF
        IF (Salinity(SalinityPerm(I)) > MaxSalinity) THEN
          AverageCorrectedMaxValue=AverageCorrectedMaxValue+Salinity(SalinityPerm(I))
          MaxSalinityValue = MAX(Salinity(SalinityPerm(I)), MaxSalinityValue)
          Salinity(SalinityPerm(I))=MaxSalinity          
          K=K+1
        END IF
      END DO
      WRITE(Message,*) 'Corrected ',J,' values that where smaller than ',MinSalinity,'.'
      CALL INFO(SolverName,Message,Level=5)
      WRITE(Message,*) 'Min. corrected salinity value:',  MinSalinityValue,&
           '. Average corrected:',AverageCorrectedMinValue/(1.0_dp*J)      
      CALL INFO(SolverName,Message,Level=5)
      WRITE(Message,*) 'Corrected ',K,' values that where lager than ',MaxSalinity,'.'
      CALL INFO(SolverName,Message,Level=5)
      WRITE(Message,*) 'Max. corrected salinity value:',  MaxSalinityValue,&
           '. Average corrected:',AverageCorrectedMaxValue/(1.0_dp*J)      
      CALL INFO(SolverName,Message,Level=5)
    END IF

    ! non-linear iteration converged?
    !--------------------------------
    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
  END DO

  CALL DefaultFinish()
CONTAINS

  !------------------------------------------------------------------------------
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSolute(  Element, ElementID, NoElements, n, nd,&
       CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, ActiveMassMatrix,&
       ExtForce)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
!!$    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
!!$         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial,ActiveMassMatrix, ExtForce !GivenGWflux, 
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: vstarAtIP(3)   ! needed in equation
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),TemperatureTimeDer,PressureTimeDer,JgwDAtIP(3),&
         KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,KgwppAtIP(3,3),fwAtIP,mugwAtIP!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP, & ! material properties at IP
         CcYcTAtIP, CcYcPAtIP, CcYcYcAtIP, rhocPAtIP, rhocYcAtIP, rhocTAtIP,& ! material properties at IP
         DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3), fcAtIP(3), extforceFlux(3),&
         DispersionCoefficient, MolecularDiffusionCoefficent ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n), XiBefore
    REAL(KIND=dp), POINTER :: gWork(:,:), XiAtIp(:)
    INTEGER :: i,j,t,p,q,DIM, RockMaterialID, IPPerm
    INTEGER, POINTER :: XiAtIPPerm(:)
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,&
         ConstantDispersion=.FALSE.,ConstantDiffusion=.FALSE.,CryogenicSuction=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixSolute)'
    TYPE(Variable_t), POINTER  :: XiAtIPVar
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    XiAtIPVar => VariableGet( Solver % Mesh % Variables, 'Xi')
    IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
      WRITE(Message,*) 'Variable Xi is not associated'
      CALL FATAL(SolverName,Message)
    END IF
    XiAtIPPerm => XiAtIPVar % Perm
    XiAtIp => XiAtIPVar % Values

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Solute Source', Found )

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    DispersionCoefficient = GetConstReal(Material,"Dispersion Coefficient", ConstantDispersion)
    MolecularDiffusionCoefficent = &
         GetConstReal(Material,"Molecular Diffusion Coefficent", ConstantDiffusion)
    
    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    !      PRINT *,"Here0"
    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL INFO(FunctionName,'Number of Gauss points for 1st element:'&
          //TRIM(I2S(IP % n)),Level=7)
      CALL Info(FunctionName,'Elemental n:'//TRIM(I2S(n))//' nd:'&
          //TRIM(I2S(nd))//' nd:'//TRIM(I2S(nb)),Level=7)
    END IF

    DO t=1,IP % n
      IPPerm = XiAtIPPerm(ElementID) + t

      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = ListGetElementReal( Temperature_h,  Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
      PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
      PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
      SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')
      gradpAtIP = ListGetElementRealGrad( Pressure_h,dBasisdx,Element,Found)
      gradTAtIP = ListGetElementRealGrad( Temperature_h,dBasisdx,Element,Found)

      vstarAtIP = 0.0_dp ! CHANGE to SUM(  Basis(1:N) * NodalRockVelocity(1:N) )

      !Materialproperties needed at IP

      ! water/ice densitities
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)
      !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        XiBefore =  XiAtIP(IPPerm)
        CALL  GetXiHartikainen(RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
        !PRINT *, "SoluteTransport", XiAtIP(IPPerm), IPPerm, XiBefore
      END SELECT

      ! solute and rock densities and derivatives
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      rhocPAtIP = rhocP(CurrentSoluteMaterial,rhocAtIP,ConstVal)
      rhocYcAtIP = rhocYc(CurrentSoluteMaterial,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
      rhocTAtIP = rhocT(CurrentSoluteMaterial,rhocAtIP,T0,TemperatureAtIP,ConstVal)
      rhosAtIP = rhos(RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!

      ! capacities of solutes      
      CcYcTAtIP = CcYcT(rhocTAtIP,PorosityAtIP,SalinityAtIP)
      CcYcPAtIP = CcYcP(rhocPAtIP,PorosityAtIP, SalinityAtIP)
      CcYcYcAtIP = CcYcYc(rhocAtIP,rhocYcAtIP,PorosityAtIP, SalinityAtIP)

      ! groundwater viscosity is pulled outtside flux computation, as needed anyway
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)

      !PRINT *, "Solute: Compute Flux"
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      KgwAtIP = GetKgw(RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP(IPPerm),MinKgw)
      !PRINT *, "Solute: Kgw", KgwAtIP(1,1)
      fwAtIP = fw(RockMaterialID,CurrentSolventMaterial,&
           Xi0tilde,rhowAtIP,XiAtIP(IPPerm),GasConstant,TemperatureAtIP)
      KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
      IF (CryogenicSuction) THEN
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
      ELSE
        KgwppAtIP = KgwAtIP
      END IF
      !PRINT *,"Solute: KgwppAtIP",KgwppAtIP
      rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP)
      !IF (SalinityAtIP > 0.2_dp) THEN
      !  PRINT *,"Solute: rhogw", rhogwAtIP, rhowAtIP,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP          
      !END IF
      ! gradT and gradP have been moved upwards, as needed elsewhere

      JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,&
           Gravity,rhogwAtIP,DIM,CryogenicSuction)
      !PRINT *, 'SoluteTransport: Jgw=', JgwDAtIP(1)*(365.25*3600*24), JgwDAtIP(2)*(365.25*3600*24)
      !PRINT *, KgwppAtIP(1,1),KgwAtIP(1,1),gradpAtIP,&
      !     Gravity,rhogwAtIP,DIM,CryogenicSuction

!!$      END IF

      ! parameters for diffusion-dispersion flow
      !r12AtIP = GetR(CurrentSoluteMaterial,GasConstant,rhocAtIP,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP)
      r12AtIP = &
           GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,rhowAtIP,rhocAtIP,&
           XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP)
      !IF (r12AtIP(2) > 1.2_dp) PRINT *,"Salinity: R2", r12AtIp(2)

      !DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)

      !PRINT *, "Solute: SalinityAtIP", SalinityAtIP
      !PRINT *, "Solute: Dm", DmAtIP,CurrentSoluteMaterial % rhoc0,CurrentSolventMaterial % muw0,TemperatureAtIP
      IF (ConstantDispersion) THEN
        KcAtIP = GetConstKC(DispersionCoefficient)
        !PRINT *,"DispersionCoefficient",KcAtIP
      ELSE
        IF (ConstantDiffusion) THEN
          DmAtIP = MolecularDiffusionCoefficent
        ELSE
          DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,CurrentSoluteMaterial % rhoc0,&
               CurrentSolventMaterial % muw0,TemperatureAtIP)
        END IF
        KcAtIP = GetKc(RockMaterialID,DmAtIP,XiAtIP(IPPerm),JgwDAtIP,PorosityAtIP)
        !PRINT *,"Solute: Kc", KcAtIP(1,1), DmAtIP,XiAtIP(IPPerm),JgwDAtIP(1:2),PorosityAtIP
      END IF
      KcYcYcAtIP = GetKcYcYc(KcAtIP,r12AtIP)
      !PRINT *,"Solute: KcYcYc", KcYcYcAtIP(1,1)
      fcAtIP = GetFc(rhocAtIP,rhowAtIP,Gravity,r12AtIP,XiTAtIP,XiPAtIP,XiAtIP(IPPerm),gradPAtIP,gradTAtIP) 

      Weight = IP % s(t) * DetJ
      !PRINT *,"Solute:",DIM,Weight
      !PRINT *, "Solute: rhoc, Porosity", rhocAtIP,PorosityAtIP

      !PRINT *, rhocAtIP
      DO p=1,nd
        DO q=1,nd
          ! diffusion term (Porosity rhoc KcYcYc.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + &
                   Weight * PorosityAtIP * rhocAtIP * KcYcYcAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
            END DO
          END DO
          !PRINT *, "Solute:  KcYcYcAtIP", KcYcYcAtIP(1,1)
          ! advection term (CcYcYc * (v* .grad(u)),v) ! V* not implemented yet and set to zero
          ! -----------------------------------
          !STIFF (p,q) = STIFF(p,q) &
          !     + Weight * CcYcYcAtIP * SUM(vstarAtIP(1:DIM)*dBasisdx(q,1:dim)) * Basis(p)

          ! left overs from partial integration of fluxes
          ! (rhoc/Xi) (u, Jgw.grad(v))
          STIFF (p,q) = STIFF(p,q) &
               - Weight * rhocAtIP * Basis(q) * SUM(JgwDAtIP(1:DIM) * dBasisdx(p,1:DIM))/XiAtIP(IPPerm)
          !PRINT *,'SoluteTransport:', JgwDAtIP(1:DIM), rhocAtIP
          ! porosity rhoc  (u,(Kc.fc).grad(v))
          
 

          ! time derivative (CcYcYc*du/dt,v):
          ! ------------------------------
          IF (ActiveMassMatrix) &
               MASS(p,q) = MASS(p,q) + Weight * CcYcYcAtIP  * Basis(q) * Basis(p)
        END DO
      END DO

      ! if we use ext. force
      IF (ExtForce) THEN
        DO p=1,nd
          DO q=1,nd
            DO i=1,DIM
              extforceFlux =  SUM(KcAtIP(i,1:DIM)*fcAtIP(1:DIM))
            END DO
            STIFF (p,q) = STIFF(p,q) &
                 - Weight * PorosityAtIP * rhocAtIP * Basis(q) * SUM(extforceFlux(1:DIM) * dBasisdx(p,1:DIM))
          END DO
        END DO
      END IF

      LoadAtIP = LoadAtIP !+ TemperatureTimeDer * CcYcTAtIP + PressureTimeDer * CcYcPAtIP 


      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)

    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSolute
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSolute(Element, n, nd,&
       CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    USE DefUtils
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial!GivenGWflux, 
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), ImposedSalinity(n), JgwDN(n),F,JgwDNAtIP,Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP, gradTAtIP(3),gradPAtIP(3)       
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    LOGICAL :: Stat,Found,ConstVal,FluxCondition,GWFluxCondition,WeakDirichletCond,ConstantsRead=.FALSE.
    INTEGER :: i,t,p,q,DIM,body_id, other_body_id, material_id, RockMaterialID
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition, ParentMaterial
    TYPE(Element_t), POINTER ::  ParentElement
    TYPE(Nodes_t) :: Nodes    
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostSoluteTransport (LocalMatrixBCSolute)'
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde, XiTAtIP, XiPAtIP, XiYcAtIP,XiEtaAtIP
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP,&
         rhowAtIP, rhoiAtIP, rhocAtIP !needed by XI
    SAVE Nodes,DIM, ConstantsRead, N0, DeltaT, T0, p0, eps, GasConstant, Gravity
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      !PRINT *, "BCSolute: (Constantsread) ", GasConstant, N0, DeltaT, T0, p0, eps, Gravity, ConstantsRead
    END IF


    ! inquire parent element and material
    other_body_id = Element % BoundaryInfo % outbody
    IF (other_body_id < 1) THEN ! only one body in calculation
      ParentElement => Element % BoundaryInfo % Right
      IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => Element % BoundaryInfo % Left
    ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
      ParentElement => Element % BoundaryInfo % Right
      IF (ParentElement % BodyId == other_body_id) ParentElement => Element % BoundaryInfo % Left
    END IF
    ! all the above was just so we can get the material properties of the parent element...
    body_id = ParentElement % BodyId
    material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', Found)
    IF (.NOT.Found) CALL FATAL(FunctionName,'Parent Material ID not found')

    ParentMaterial => Model % Materials(material_id) % Values
    IF (.NOT. ASSOCIATED(ParentMaterial)) THEN
      WRITE(Message,*)&
           'No material values found for body no ', body_id,&
           ' under material id ', material_id
      CALL FATAL(FunctionName,Message)
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(ParentElement)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ParentElement % ElementIndex ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(ParentMaterial,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(ParentMaterial,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    !Check, whether we have a prescribed solute flow
    Flux(1:n)  = GetReal( BoundaryCondition,'Solute Flow', FluxCondition )
    ! check, whether we have a prescribed groundwater flux
    JgwDN(1:n)  = GetReal( BoundaryCondition,'Groundwater Flux', GWFluxCondition )
    ! Check, whether we have a weakly imposed Dirichlet condition
    ImposedSalinity(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)

    ! if none of the above, we can call it a day
    IF (.NOT.(FluxCondition .OR. GWFluxCondition .OR. WeakDirichletCond)) RETURN

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL INFO(FunctionName,'Number of Gauss points for 1st element:'&
          //TRIM(I2S(IP % n)),Level=7)
    END IF
    
    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    !PRINT *,"BCSolute:",deltaInElement,eps,DeltaT,T0,GasConstant

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! we need XiAtIP and rhocAtIP only if we have a non-zero groundwater flux or solute flow
      IF (FluxCondition .OR. GWFluxCondition) THEN ! else spare us the computation

        ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
        TemperatureAtIP = ListGetElementRealParent( Temperature_h, Basis, Element, Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
        PorosityAtIP = ListGetElementRealParent( Porosity_h, Basis, Element, Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
        PressureAtIP = ListGetElementRealParent( Pressure_h, Basis, Element, Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
        SalinityAtIP = ListGetElementRealParent( Salinity_h, Basis, Element, Found)
        IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')

        !Materialproperties needed at IP

        ! water/ice densitities
        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal) !!
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)
        !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,GlobalRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          ! NB: XiTAtIP, XiPAtIP not needed
        CASE DEFAULT ! Hartikainen model
          CALL  GetXiHartikainen(RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .TRUE.,.FALSE., .FALSE., .FALSE.,.FALSE.) ! we need to compute, as IP's on boundary elements deviate from bulk
          ! NB: XiTAtIP, XiPAtIP, XiYcAtIP not needed
          !PRINT *, "SoluteTransportBC", XiAtIP
        END SELECT
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)

        ! Check, whether we have a weakly imposed Dirichlet condition
        IF (GivenGWFlux) THEN
          JgwDNAtIP = SUM(Basis(1:n)*JgwDN(1:n))
          ! contribution from partial integration of groundwater flux term (always on)
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) &
                   + Weight * rhocAtIP * Basis(q) * Basis(p) * JgwDNAtIP/XiAtIP
            END DO
          END DO
        END IF

        !PRINT *,"BCSolute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

        ! Given flux:
        ! -----------
        IF (Fluxcondition) THEN
          F = SUM(Basis(1:n)*Flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) - Weight * PorosityAtIP * rhocAtIP * F * Basis(1:nd)
          !PRINT *,"Salinity BC: Flux:", F, PorosityAtIP , rhocAtIP , Weight, Weight * PorosityAtIP * rhocAtIP * F * Basis(1:nd)
        END IF
      END IF

      ! Given salinity, weakly imposed
      !----------------------------------------------------------------------
      IF (WeakDirichletCond) THEN
        SalinityAtIP = SUM(Salinity(1:n)*Basis(1:n))
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * C * SalinityAtIP * Basis(1:nd)
      END IF
    END DO
    !PRINT *, "Salinity BC: Flux:",SUM(FORCE(1:nd))
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCSolute
  !------------------------------------------------------------------------------

END SUBROUTINE PermafrostSoluteTransport
