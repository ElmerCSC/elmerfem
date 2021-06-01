!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) asny later version.
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
!>  Solver for permafrost heat transfer problem including phase change
!---------------------------------------------------------------------------------------------  
!-----------------------------------------------------------------------------
!> heat transfer equation for enhanced permafrost model
!-----------------------------------------------------------------------------
SUBROUTINE PermafrostHeatTransfer_init( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: DT
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation_init'
  LOGICAL :: OutputXi = .FALSE. , OutputFlux=.FALSE., Found
  TYPE(Variable_t), POINTER :: XiAtIPVar
  INTEGER, POINTER :: XiAtIPPerm(:)
  REAL(KIND=dp), POINTER :: XiAtIp(:)
  INTEGER :: I
  !------------------------------------------------------------------------------
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, '  Initializing heat transfer         ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  SolverParams => GetSolverParams()
  OutputXi = GetLogical(SolverParams, 'Output Xi', Found)
  IF (.NOT.Found) OutputXi = .FALSE.
  !PRINT *,SolverName,OutputXi
  IF (OutputXi) THEN
    CALL INFO(SolverName,'Output of IP variable "Xi" ',Level=6)
    CALL ListAddString( SolverParams, NextFreeKeyword('Exported Variable',SolverParams),'-IP Xi' )
  ELSE
    CALL INFO(SolverName,'No output of IP variable "Xi" ',Level=6)
    CALL ListAddString( SolverParams, NextFreeKeyword('Exported Variable',SolverParams),'-IP -nooutput Xi' )
  END IF
  CALL INFO(SolverName,'Added variable Xi',Level=6)
  
  ! Add linear system defaults: BiCGStab+ILU0
  CALL ListAddNewString(SolverParams,'Linear System Solver','Iterative')
  CALL ListAddNewString(SolverParams,'Linear System Iterative Method','BiCGStab')
  CALL ListAddNewString(SolverParams,'Linear System Preconditioning','ILU0')
  CALL ListAddNewInteger(SolverParams,'Linear System Max Iterations',500)
  CALL ListAddNewInteger(SolverParams,'Linear System Residual Output',10)
  CALL ListAddNewConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-08_dp)
  ! Add Nonlinear system defaults
  CALL ListAddNewConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0e-05_dp)
  CALL ListAddNewInteger(SolverParams,'Nonlinear System Max Iterations',50) 

  CALL INFO( SolverName, ' Done Initializing      ',Level=6 )
END SUBROUTINE PermafrostHeatTransfer_init
!------------------------------------------------------------------------------
!> heat transfer equation for enhanced permafrost model
!> \ingroup Solvers
SUBROUTINE PermafrostHeatTransfer( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: DummyGWfluxVar
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  REAL(KIND=dp) :: Norm
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE., FluxOutput = .FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists=.FALSE.,&
       InitializeSteadyState=.FALSE.,ActiveMassMatrix=.TRUE.,&
       NoSalinity=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName, XiAtIPName
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h,&
       Vstar1_h, Vstar2_h, Vstar3_h

  SAVE DIM,FirstTime,AllocationsDone,FluxOutput,DepthName,XiAtIPName,&
       CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,ComputeDt,DepthExists,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h, &
       Vstar1_h, Vstar2_h, Vstar3_h
  !------------------------------------------------------------------------------
  CALL INFO( SolverName, '-------------------------------------',Level=4 )
  CALL INFO( SolverName, 'Computing heat transfer              ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------',Level=4 )

  IF (FirstTime) THEN
    DIM = CoordinateSystemDimension()
    ! Handles to other system variables
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )

    ! Handles to time derivatives of system variables
    CALL ListInitElementKeyword( PressureVelo_h, 'Material', 'Pressure Velocity Variable' )
    CALL ListInitElementKeyword( SalinityVelo_h, 'Material', 'Salinity Velocity Variable' )

    ! Handle to Heat Source (possible description of heat source at elements/IP's) 
    CALL ListInitElementKeyword( Load_h,'Body Force','Heat Source' )

    ! Handles to advection velocities
    CALL ListInitElementKeyword( Vstar1_h,'Material','Convection Velocity 1')
    CALL ListInitElementKeyword( Vstar2_h,'Material','Convection Velocity 2')
    IF (DIM > 2) &
         CALL ListInitElementKeyword( Vstar3_h,'Material','Convection Velocity 3')
    
    ! Handles to other variables
    CALL ListInitElementKeyword( Depth_h, 'Material', 'Depth Variable' )

  END IF
  
  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)

  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1
  
  ! check, whether we assume steady state (despite transient run)
  ! this can come handy to produce a balance-pressure field at the
  ! start of the simulation
  !---------------------------------------------------------------  
  IF (FirstTime) &
       InitializeSteadyState = GetLogical(Params,'Initialize Steady State',Found)
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
  
  ! check, whether an output variable for groundwater flux exists
  !--------------------------------------------------------------
  IF (FirstTime) THEN
    FluxOutput = GetLogical(Params,'Computed Groundwater Flux',Found)
    IF (.NOT.Found) FluxOutput = .FALSE.
    IF (FluxOutput) THEN
      DO I=1,DIM
        DummyGWfluxVar => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux '//TRIM(I2S(i)))
        FluxOutput = ASSOCIATED(DummyGWfluxVar)
        IF (.NOT.FluxOutput) EXIT
      END DO
      IF (FluxOutput) THEN
        WRITE (Message,*) 'Groundwater flow will be read from: Groundwater Flux {1..',DIM,'}'
        CALL INFO(SolverName,Message,Level=5)
      END IF
    END IF
  END IF  
  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL INFO( SolverName, "Nonlinear iteration "&
        //TRIM(I2S(iter))//" out of "//TRIM(I2S(maxiter)),Level=4)
    
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
          CALL INFO(SolverName,Message,Level=6)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=6)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (GlobalRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=6)
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
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF

      CALL LocalMatrixHTEQ(  Element, Element % ElementIndex, Active, n, nd+nb,&
           CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial,&
           ActiveMassMatrix,FluxOutput)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()

    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCHTEQ(  Element, n, nd+nb )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO
  
  CALL DefaultFinish()

CONTAINS

  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixHTEQ(  Element, ElementID, NoElements, n, nd,&
       CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial,&
       ActiveMassMatrix,FluxOutput)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial,ActiveMassMatrix, FluxOutput
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: DepthAtIP,RefDepth,CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,&
         KgwppAtIP(3,3),fwAtIP,mugwAtIP,DtdAtIP(3,3)!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         PressureVeloAtIP,SalinityVeloAtIP,&
         StiffPQ, meanfactor, vstarAtIP(3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,IPPerm,DIM, RockMaterialID, FluxDOFs
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,&
         CryogenicSuction=.FALSE.,HydroGeo=.FALSE.,ComputeFlux=.TRUE.,&
         NoSalinity=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
    TYPE(Variable_t), POINTER :: XiAtIPVar, GWfluxVar1, GWfluxVar2, GWfluxVar3
    INTEGER, POINTER :: XiAtIPPerm(:),GWfluxPerm(:)
    REAL(KIND=dp), POINTER :: XiAtIP(:), FluxAtElem(:)

    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, ConstVal,DIM, GasConstant, N0,DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    gradTAtIP = 0.0_dp
    gradPAtIP = 0.0_dp
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    XiAtIPVar => VariableGet( Solver % Mesh % Variables, 'Xi')
    IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
      WRITE(Message,*) 'Variable Xi is not associated'
      CALL FATAL(SolverName,Message)
    END IF
    XiAtIPPerm => XiAtIPVar % Perm
    XiAtIp => XiAtIPVar % Values

    IF (FluxOutput) THEN
      GWfluxVar1 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 1')
      GWfluxPerm => GWfluxVar1 % Perm
      IF (DIM > 1) THEN
        GWfluxVar2 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 2')
      END IF	
      IF (DIM == 3) THEN
        GWfluxVar3 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 3')
      END IF
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)

    NoSalinity = GetLogical(Material,'No Salinity',Found)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)

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

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

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
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      ! The heat source term
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
      !IF (LoadAtIP > 0.0_dp) PRINT *,"HTEQ:LoadAtIP", LoadAtIP      
      ! Contribution from Radiogenic Heat Production
      DepthAtIP = ListGetElementReal( Depth_h, Basis, Element, DepthExists, GaussPoint=t)
      IF (DepthExists) &
           RefDepth = GetConstReal(Material,'Radiogenic Reference Depth',DepthExists)
      IF (DepthExists)  &
        LoadAtIP = LoadAtIP  + RadiogenicHeatProduction(RockMaterialID,DepthAtIP,RefDepth)

      ! System variables (Temperature, Porosity, Pressure, Salinity) at IP
      PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
      PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
      IF (NoSalinity) THEN
        SalinityAtIP = 0.0_dp
      ELSE
        SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) THEN
          CALL INFO(SolverName,'Salinity not found - setting to zero',Level=7)
          NoSalinity=.TRUE.
        END IF
      END IF
      TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
      !IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')



      ! Time derivatives of system variables at IP
      PressureVeloAtIP = 0.0_dp
      SalinityVeloAtIP = 0.0_dp
      IF (ActiveMassMatrix) THEN
        PressureVeloAtIP = ListGetElementReal( PressureVelo_h, Basis, Element, Found, GaussPoint=t)
        SalinityVeloAtIP = ListGetElementReal( SalinityVelo_h, Basis, Element, Found, GaussPoint=t)
      END IF

      ! bedrock deformation velocity at IP
      vstarAtIP = 0.0_dp
      vstarAtIP(1) = ListGetElementReal( Vstar1_h, Basis, Element, Found, GaussPoint=t)
      vstarAtIP(2) = ListGetElementReal( Vstar2_h, Basis, Element, Found, GaussPoint=t)
      IF (DIM > 2) &
           vstarAtIP(3) = ListGetElementReal( Vstar3_h, Basis, Element, Found, GaussPoint=t)

      !Materialproperties needed for computing Xi at IP
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      IPPerm = XiAtIPPerm(ElementID) + t
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
        CALL  GetXiHartikainen (RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
      rhosAtIP = rhos(RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"HTEQ: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

      ! heat capacities
      csAtIP   = cs(RockMaterialID,&
           T0,TemperatureAtIP,ConstVal)
      cwAtIP   = cw(CurrentSolventMaterial,&
           T0,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP,ConstVal)
      ciAtIP   = ci(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      ccAtIP   = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)

      ! latent heat
      hiAtIP = hi(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      hwAtIP = hw(CurrentSolventMaterial,&
           T0,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP,ConstVal)

      ! heat conductivity at IP
      ksthAtIP = GetKalphath(GlobalRockMaterial % ks0th(RockMaterialID),&
           GlobalRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
      !IF (ksthAtIP > 3.01_dp) &
!     PRINT *,GlobalRockMaterial % ks0th(RockMaterialID), RockMaterialID, ElementID
      kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP(IPPerm),&
           SalinityATIP,PorosityAtIP,meanfactor)

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP(IPPerm),XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP(IPPerm),SalinityAtIP)

      ! compute groundwater flux for advection term
      CGTpAtIP = GetCGTp(rhoiAtIP,hiAtIP,hwAtIP,XiPAtIP,PorosityAtIP)
      CGTycAtIP = GetCGTyc(rhoiAtIP,hiAtIP,hwAtIP,XiYcAtIP,PorosityAtIP)

      ! groundwater flux
      !-----------------
      IF (FluxOutput) THEN
        JgwDAtIP(1) = GWfluxVar1 % Values(GWfluxPerm(ElementID) + t)
	IF (DIM > 1) THEN
          JgwDAtIP(2) = GWfluxVar2 % Values(GWfluxPerm(ElementID) + t)
          IF (DIM == 3) &
               JgwDAtIP(3) = GWfluxVar3 % Values(GWfluxPerm(ElementID) + t)
	END IF
      ELSE
        JgwDAtIP = 0.0_dp
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP(IPPerm),MinKgw)
        fwAtIP = fw(RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP(IPPerm),GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        !PRINT *,"HTEQ: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP)

        ! gradients at IP
        gradTAtIP = ListGetElementRealGrad( Temperature_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Temperature gradient')
        gradpAtIP = ListGetElementRealGrad( Pressure_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Pressure gradient')

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        !PRINT *,"HTEQ: JgwD=(",JgwDAtIP(1:DIM)*365.5*24.0*3600.0,")"        
      END IF
      
      ! add thermal dispersion in Hydro-Geological Mode
      !------------------------------------------------
      IF (HydroGeo) THEN
        DtdAtIP = GetDtd(RockMaterialID,XiAtIP(IPPerm),PorosityAtIP,JgwDAtIP)
        DO I=1,DIM
          DO J=1,DIM
            KGTTAtIP(I,J) = KGTTAtIP(I,J) + CGWTTAtIP * DtdAtIP(I,J)
          END DO
        END DO
      END IF

      Weight = IP % s(t) * DetJ
      !PRINT *, "Weight=", Weight
      !KGTTAtIP = 0
      !KGTTAtIP(1,1) = 3.0_dp
      !KGTTAtIP(2,2) = 3.0_dp
      DO p=1,nd
        DO q=1,nd          
          ! diffusion term (KGTTAtIP.grad(u),grad(v)):
          DO i=1,DIM
            DO j=1,DIM
              Stiff(p,q) = Stiff(p,q) + Weight * KGTTAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)
              !PRINT *,"cond", Weight," *", KGTTAtIP(i,j)," *", dBasisdx(p,j),"*", dBasisdx(q,i)
            END DO
          END DO
          ! advection term due to groundwater velocity (CgwTT * (Jgw.grad(u)),v)
          ! --------------------------------------------------------------------
          STIFF (p,q) = STIFF(p,q) +&
               Weight * CgwTTAtIP * SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p) !
          !PRINT *,"adv", Weight," *", CgwTTAtIP," *", SUM(JgwDAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)
          ! advection term due to bedrock velocity (CGTT * (vstar.grad(u)),v)
          ! --------------------------------------------------------------------
          STIFF (p,q) = STIFF(p,q) +&
               Weight * CGTTAtIP * SUM(vstarAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)
          ! time derivative (CGTT*du/dt,v):
          ! ------------------------------
          IF (ActiveMassMatrix) &
               MASS(p,q) = MASS(p,q) + Weight * (CGTTAtIP) * Basis(q) * Basis(p) !
          !PRINT *,"storage", CGTTAtIP, "*",Basis(q) * Basis(p) 
        END DO
      END DO
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
      ! temperature force due to pressure change (CGTp * (dp/dt + vstar grad p),u)
      FORCE(1:nd) = FORCE(1:nd) - &
           Weight * CGTpAtIP*(PressureVeloAtIP + SUM(vstarAtIP(1:DIM)*gradPAtIP(1:DIM))) * Basis(1:nd)
      ! temperature force due to salinity change (CGTy * (dyc/dt + vstar grad yc),u)
      FORCE(1:nd) = FORCE(1:nd) - &
           Weight * CGTycAtIP*(SalinityVeloAtIP + SUM(vstarAtIP(1:DIM)*gradPAtIP(1:DIM))) * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixHTEQ
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCHTEQ( Element, n, nd )
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Fluxcondition,Robincondition
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BoundaryCondition

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BoundaryCondition,'Heat Flux', FluxCondition )
    Coeff(1:n) = GetReal( BoundaryCondition,'Heat Transfer Coefficient', RobinCondition )
    Ext_t(1:n) = GetReal( BoundaryCondition,'External Temperature', RobinCondition )

    IF (FluxCondition .OR. RobinCondition)  THEN
      ! Numerical integration:
      !-----------------------
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ

        ! Evaluate terms at the integration point:
        !------------------------------------------

        ! Given flux:
        ! -----------
        F = SUM(Basis(1:n)*flux(1:n))

        ! Robin condition (C*(u-u_0)):
        ! ---------------------------
        C = SUM(Basis(1:n)*coeff(1:n))
        Ext = SUM(Basis(1:n)*ext_t(1:n))

        IF (Robincondition) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C*Ext * Basis(1:nd)
        ELSE IF (Fluxcondition) THEN
          !FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
          FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
          !PRINT *,"LocalMatrixBCHTEQ:",F
        END IF
      END DO
    END IF
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCHTEQ
  !------------------------------------------------------------------------------

END SUBROUTINE PermafrostHeatTransfer
!------------------------------------------------------------------------------
