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
!>  Solver for saturated gorundwater flux including solutes and phase change
!---------------------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlow_Init( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: DT
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundwaterFlow_init'
  TYPE(Variable_t), POINTER :: ReferenceVar
  
  LOGICAL :: OffsetDensity = .FALSE. , Found
  !------------------------------------------------------------------------------

  CALL INFO( SolverName, '-------------------------------------------',Level=4 )
  CALL INFO( SolverName, '  Initializing Permafrost Groundwater Flow ',Level=4 )
  CALL INFO( SolverName, '-------------------------------------------',Level=4 )
  SolverParams => GetSolverParams()
  
  IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    CALL ListAddString( SolverParams, 'Variable', 'GWPressure' )
    CALL WARN( SolverName, 'Variable not found. Adding default "GWPressure"')
  END IF

  ! Add linear system defaults: BCGStab+ILU0
  CALL ListAddNewString(SolverParams,'Linear System Solver','Iterative')
  CALL ListAddNewString(SolverParams,'Linear System Iterative Method','BiCGStab')
  CALL ListAddNewString(SolverParams,'Linear System Preconditioning','ILU0')
  CALL ListAddNewInteger(SolverParams,'Linear System Max Iterations',500)
  CALL ListAddNewInteger(SolverParams,'Linear System Residual Output',10)
  CALL ListAddNewConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-08_dp)
  ! Add Nonlinear system defaults
  CALL ListAddnewConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0e-05_dp)
  CALL ListAddNewInteger(SolverParams,'Nonlinear System Max Iterations',50) 
  
  CALL INFO( SolverName, '  Done Initialization',Level=4)
  CALL INFO( SolverName, '-------------------------------------------',Level=4 )
END SUBROUTINE PermafrostGroundwaterFlow_Init
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlow( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, DummyDtVar,SalinityDtVar,&
       DummyGWfluxVar,StressInvVar
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, Active,iter, maxiter, istat,StressInvDOFs
  INTEGER,PARAMETER :: io=22
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), DummyDtPerm(:),SalinityDtPerm(:),&
       StressInvPerm(:),DummyGWfluxPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), DummyDt(:),SalinityDt(:),&
       DummyGWflux(:),StressInv(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.,GivenGWFlux,ElementWiseRockMaterial, DummyLog=.FALSE.,&
       InitializeSteadyState=.FALSE., ActiveMassMatrix=.TRUE., ComputeDeformation=.FALSE.,DeformationExists=.FALSE.,&
       StressInvAllocationsDone=.FALSE.,StressInvDtAllocationsDone=.FALSE.,&
       HydroGeo=.FALSE.,ComputeDt=.FALSE.,FluxOutput=.FALSE.,&
       TemperatureTimeDerExists=.FALSE.,SalinityTimeDerExists=.FALSE., OffsetDensity=.FALSE., &
       ComputeFreshwaterHead=.FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  !CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, StressInvName, &
  CHARACTER(LEN=MAX_NAME_LEN) ::    VarName,PhaseChangeModel,ElementRockMaterialName,DeformationName
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       TemperatureBC_h, PressureBC_h, SalinityBC_h, PorosityBC_h, Recharge_h, GWFlux_h,&
       TemperatureDt_h, SalinityDt_h, StressInv_h,StressInvDt_h,Vstar1_h, Vstar2_h, Vstar3_h
  TYPE(VariableHandle_t) :: NormalVar_h
  

  ! Additional variables to output salinity at BCs with conditional dirichlet 
  INTEGER, POINTER :: BCFluxPerm(:)
  INTEGER, SAVE :: BCFluxNodes 
  TYPE(Variable_t), POINTER :: BCFluxVar
    
  
  SAVE DIM,FirstTime,AllocationsDone,CurrentSoluteMaterial,CurrentSolventMaterial,&
       ElementWiseRockMaterial, ComputeDeformation, FluxOutput,&
       StressInvAllocationsDone, StressInvDtAllocationsDone, OffsetDensity, &
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       TemperatureDt_h, SalinityDt_h, StressInv_h, StressInvDt_h, &
       Vstar1_h, Vstar2_h, Vstar3_h, NormalVar_h, Recharge_h, GWFlux_h,&
       TemperatureBC_h, PressureBC_h, SalinityBC_h, PorosityBC_h,&
       ActiveMassMatrix, InitializeSteadyState, HydroGeo, ComputeDt, &
       ComputeFreshwaterHead
  !------------------------------------------------------------------------------
  CALL DefaultStart()

  Params => GetSolverParams()

  ! Initial settings
  !---------------------------------------------------------------  
  IF (FirstTime) THEN
    ! check, whether we assume steady state (despite transient run)
    ! this can come handy to produce a balance-pressure field at the
    ! start of the simulation
    InitializeSteadyState = GetLogical(Params,'Initialize Steady State',Found)
    ! inquire whether to include time-derivative terms in force vector
    ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)
    IF (ComputeDt) THEN
      CALL INFO(SolverName,"Computing time derivatives in force vector",Level=4)
    ELSE
      CALL INFO(SolverName,"Omitting time derivatives in force vector",Level=4)
    END IF
    ! inquire whether to compute deformation force term
    ComputeDeformation = GetLogical(Params,'Compute Deformation',Found)
    IF (ComputeDeformation) THEN
      CALL INFO(SolverName,"Including stress invariant derivative in force vector",Level=4)
    ELSE
      CALL INFO(SolverName,"Omitting stress invariant derivative in force vector",Level=4)
    END IF
    ComputeFreshwaterHead = GetLogical(Params,'Compute Freshwater Head',Found)
      IF (ComputeFreshwaterHead) &
           CALL INFO(SolverName,'Computing freshwater head',Level=4)
  END IF
  
  IF (InitializeSteadyState) THEN
    IF (GetTimeStep() == 1) THEN
      CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=4)
      ActiveMassMatrix = .FALSE.
    ELSE 
      CALL INFO(SolverName,"Switching mass matrix to active after initializing with steady state",Level=4)
      ActiveMassMatrix = .TRUE.
      InitializeSteadyState = .FALSE.
    END IF
  END IF
  
  
  maxiter = ListGetInteger( Params, &
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1


  !FluxOutput = GetLogical(Params,'Groundwater Flux Output',Found)
  !StressInvName =  ListGetString(params,'Ground Stress Invariant Variable Name',ComputeDeformation)
  !DeformationName = ListGetString(params,'Ground Deformation Variable Name ',DeformationExists)
  
  ! solver variable
  Pressure => Solver % Variable % Values
  PressurePerm => Solver % Variable % Perm
  VarName = Solver % Variable % Name

  
  IF (FirstTime) THEN
    DIM = CoordinateSystemDimension()
    ! Handles to all variables    
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )
    ! Handles to all variables at BC's
    CALL ListInitElementKeyword( PressureBC_h, 'Boundary Condition', 'Pressure Variable' )
    CALL ListInitElementKeyword( SalinityBC_h, 'Boundary Condition', 'Salinity Variable' )
    CALL ListInitElementKeyword( PorosityBC_h, 'Boundary Condition', 'Porosity Variable' )
    CALL ListInitElementKeyword( TemperatureBC_h, 'Boundary Condition','Temperature Variable' )
    !Handles to values on boundaries
    CALL ListInitElementKeyword( Recharge_h, 'Boundary Condition', 'Freshwater Recharge' )
    CALL ListInitElementKeyword( GWFlux_h, 'Boundary Condition', 'Groundwater Flux' )
    !CALL ListInitElementKeyword( StressInv_h, 'Material', 'Stress Invariant Variable' )
    ! Handles to advection velocities
    CALL ListInitElementKeyword( Vstar1_h,'Material','Convection Velocity 1')
    CALL ListInitElementKeyword( Vstar2_h,'Material','Convection Velocity 2')
    IF (DIM > 2) &
         CALL ListInitElementKeyword( Vstar3_h,'Material','Convection Velocity 3')
    ! Handles to time-derivatives
    IF (ComputeDt) THEN
      CALL ListInitElementKeyword( TemperatureDt_h, 'Material', 'Temperature Velocity Variable' )
      CALL ListInitElementKeyword( SalinityDt_h, 'Material', 'Salinity Velocity Variable' )
    END IF
    IF (ComputeDeformation) THEN
      CALL ListInitElementKeyword( StressInvDt_h, 'Material', 'Stress Invariant Velocity Variable' )
    END IF
    ! Handle for surface normal
    CALL ListInitElementVariable( NormalVar_h, 'Normal Vector' )
  END IF
  


  IF (FirstTime) THEN
    OffsetDensity = GetLogical(Model % Constants,'Permafrost Offset Density', Found)
    IF (.NOT.Found) THEN
      OffsetDensity = .FALSE.
    END IF
  END IF
  
  ! check, whether an output variable for groundwater flux exists
  !--------------------------------------------------------------
  IF (FirstTime) THEN
    DO I=1,DIM
      DummyGWfluxVar => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux '//TRIM(I2S(i)))
      FluxOutput = ASSOCIATED(DummyGWfluxVar)
      IF (.NOT.FluxOutput) EXIT
    END DO
    IF (FluxOutput) THEN
      CALL INFO(SolverName,'Groundwater flow will be written to: Groundwater Flux {1..'&
          //TRIM(I2S(DIM))//'}',Level=4)
    END IF
  END IF
  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL DefaultInitialize()
    !------------------------------------------------------------------------------
    Active = Solver % NumberOfActiveElements
    DO t=1,Active
      Element => GetActiveElement(t)
      IF (.NOT.ASSOCIATED(Element)) CYCLE      
      ! cycle halo elements
      !-------------------
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      
      Material => GetMaterial(Element)
      IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I0)') 'No Material found for element no. ',Element % ElementIndex
        CALL FATAL(SolverName,Message)
      END IF

      
      ! inquire whether to use hydro-geo simplifications
      HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
      IF(HydroGeo) THEN
        ComputeDt = .FALSE.
        IF (FirstTime) THEN
          CALL INFO(SolverName,"Using hydro-geo simplifications.",Level=9)
          CALL INFO(SolverName,"Switching time derivatives in force vector off",Level=9)
        END IF
      END IF
      
      PhaseChangeModel = ListGetString(Material, &
           'Permafrost Phase Change Model', Found )
      IF (Found) THEN
        WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
        CALL INFO(SolverName,Message,Level=9)
      END IF

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
          NumberOfRockRecords = ReadPermafrostRockMaterial( Material )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          PRINT *, "NumberOfRockRecords=", NumberOfRockRecords
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=5)
          FirstTime = .FALSE.
        END IF
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        IF (.NOT.ASSOCIATED(CurrentSolventMaterial)) &
             CALL FATAL(Solvername,'Solvent Material not associated')
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        IF (.NOT.ASSOCIATED(CurrentSoluteMaterial)) &
             CALL FATAL(Solvername,'Solute Material not associated')
      END IF

      N  = GetElementNOFNodes()
      ND = GetElementNOFDOFs()
      NB = GetElementNOFBDOFs()

      ! compose element-wise contributions to matrix and R.H.S

      CALL LocalMatrixDarcy( Model, Element, Element % ElementIndex, N, ND+NB, Active,  &
           CurrentSoluteMaterial,CurrentSolventMaterial,&
           PhaseChangeModel,ElementWiseRockMaterial, ActiveMassMatrix, &
           ComputeDt,ComputeDeformation,HydroGeo, OffsetDensity,FluxOutput, &
           ComputeFreshwaterHead)
    END DO
    OffsetDensity = .FALSE. ! make sure to not read after first time
    CALL DefaultFinishBulkAssembly()
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCDarcy(  Element, t, n, nd+nb, CurrentSoluteMaterial)
      END IF
    END DO
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    !Solve the system:
    !--------------------
    Norm = DefaultSolve()
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT


  END DO

  CALL DefaultFinish()

       
  IF( ListGetLogical( Params,'Compute BC Flux',Found ) ) THEN
    CALL Info(Solvername,'Computing flux for Dirichlet BCs for salinity',Level=6)
    DummyGWfluxVar => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux')
    IF(.NOT. ASSOCIATED(DummyGWFluxVar) ) THEN
      CALL Fatal(SolverName,'Groundwater Flux not available!')
    END IF
    BCFluxVar => VariableGet( Solver % Mesh % Variables,'BC Flux')
    IF(.NOT. ASSOCIATED( BCFluxVar ) ) THEN
      CALL Info(SolverName,'Creating permutation for boundary flux',Level=12)
      ALLOCATE( BCFluxPerm( Solver % Mesh % NumberOfNodes ) )
      BCFluxPerm = 0
      CALL MakePermUsingMask( Model, Solver, Solver % Mesh,'Salinity',.FALSE.,&
          BCFluxPerm, BCFluxNodes )
      CALL Info(SolverName,'Creating variable for boundary flux of size: '&
          //TRIM(I2S(BCFluxNodes)),Level=12)
      CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver,&
          'BC Flux', DummyGWFluxVar % DOFs, Perm = BCFluxPerm )
    END IF
    BCFluxVar => VariableGet( Solver % Mesh % Variables,'BC Flux')
    CALL Ip2DgSwapper( Solver % Mesh, DummyGWFluxVar, BCFluxVar )
  END IF
  

CONTAINS
  ! PermafrostGroundWaterFlow : Assembly of the matrix entries arising from the bulk elements 

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixDarcy( Model, Element, ElementID, n, nd, NoElements,&
       CurrentSoluteMaterial,CurrentSolventMaterial,&
       PhaseChangeModel, ElementWiseRockMaterial, ActiveMassMatrix,&
       ComputeDt,ComputeDeformation,HydroGeo,OffsetDensity,FluxOutput,ComputeFreshwaterHead)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL :: ElementWiseRockMaterial,ActiveMassMatrix,ComputeDt,&
         ComputeDeformation,HydroGeo,OffsetDensity,FluxOutput,&
         ComputeFreshwaterHead
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CgwppAtIP,CgwpTAtIP,CgwpYcAtIP,CgwpI1AtIP,KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),&
         meanfactor,MinKgw,gradTAtIP(3),gradPAtIP(3),gradYcAtIP(3),fluxTAtIP(3),fluxgAtIP(3),vstarAtIP(3) ! needed in equation
    REAL(KIND=dp) :: JgwDAtIP(3),JcFAtIP(3), DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3),&
         fcAtIP(3), DispersionCoefficient, MolecularDiffusionCoefficent ! from salinity transport
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT  XiAtIP,
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP,bijAtIP(2,2),bijYctIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: ks0th,e1,bs,rhos0,cs0,Xi0,eta0,alphaL,alphaT,RadGen,acs(0:5),&
         as0,aas(0:5),ks0,cks(0:5)  ! stuff coming from RockMaterial
    INTEGER :: acsl,aasl,cksl       ! stuff coming from RockMaterial
    REAL(KIND=dp) :: EGAtIP,nuGAtIP,kappaGAtIP ! bedrock deformation
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps,Gravity(3)! real constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,&
         rhospAtIP,rhosTAtIP,&
         rhowTAtIP,rhowPAtIP,rhowYcAtIP,&
         rhoiPAtIP,rhoiTAtIP,&
         rhocPAtIP,rhocTAtIP,rhocYcAtIP,&
         rhogwPAtIP,rhogwTAtIP,rhogwYcAtIP, rhoGAtIP, rhow0
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ,elevationAtIp
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,KPorosityAtIP,SalinityAtIP,PressureAtIP
    REAL(KIND=dp) :: TemperatureDtAtIP,SalinityDtAtIP,PressureDtAtIP,StressInvDtAtIP 
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp) , POINTER :: gWork(:,:)
    !REAL(KIND=dp) , ALLOCATABLE :: CgwpI1AtNodes(:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID, FluxDOFs,IPPerm,IPPermRhogw,IPPermFreshwaterHead
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE., ConstantDispersion=.FALSE.,&
         ConstantDiffusion=.FALSE., CryogenicSuction=.FALSE., swaptensor=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow', &
         FunctionName='Permafrost (LocalMatrixDarcy)'
    TYPE(Variable_t), POINTER :: XiAtIPVar, RhoOffsetAtIPVar,GWfluxVar1,GWfluxVar2,GWfluxVar3,FreshwaterHeadAtIPVar
    INTEGER, POINTER :: XiAtIPPerm(:), RhoOffsetAtIPPerm(:),GWfluxPerm(:),FreshwaterHeadAtIPPerm(:)
    REAL(KIND=dp), POINTER :: XiAtIP(:), RhoOffsetAtIP(:),GWFluxVal(:),FreshwaterHeadAtIP(:)
    
    SAVE Nodes, ConstantsRead, ConstVal, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
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
    XiAtIP => XiAtIPVar % Values

    ! find needed variables
    IF (ComputeFreshwaterHead) THEN
      FreshwaterHeadAtIPVar => VariableGet( Model % Mesh % Variables, 'Freshwater Head')
      IF (.NOT.ASSOCIATED(FreshwaterHeadAtIPVar)) THEN
        WRITE(Message,*) '"Compute Freshwater Head" set but variable "Freshwater Head" is not associated'
      END IF
      FreshwaterHeadAtIPPerm => FreshwaterHeadAtIPVar % Perm  
      FreshwaterHeadAtIP => FreshwaterHeadAtIPVar % Values
    END IF
    IF (OffsetDensity) THEN
      RhoOffsetAtIPVar => VariableGet( Model % Mesh % Variables, 'Reference Offset Density')
      IF (.NOT.ASSOCIATED(RhoOffsetAtIPVar)) THEN
        WRITE(Message,*) '"Permafrost Offset Density" is set, but variable "Reference Offset Density" is not associated'
        CALL FATAL(SolverName,Message)
      END IF
      RhoOffsetAtIPPerm => RhoOffsetAtIPVar % Perm  
      RhoOffsetAtIP => RhoOffsetAtIPVar % Values
    END IF
    IF (FluxOutput) THEN
      GWfluxVar1 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 1')
      GWfluxPerm => GWfluxVar1 % Perm
      GWfluxVar2 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 2')
      IF (DIM == 3) THEN
        GWfluxVar3 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 3')
      END IF
    END IF
      
    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      LOAD(1:n) = GetReal( BodyForce,'Groundwater source', Found )   
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",Found)
    IF (.NOT.Found) THEN
      CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
      meanfactor = 1.0_dp
    END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    swaptensor = GetLogical(Material,'Swap Tensor',Found)
    
    NoSalinity = GetLogical(Material,'No Salinity',Found)
    
    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (ConstVal) &
        CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    DispersionCoefficient = GetConstReal(Material,"Dispersion Coefficient", ConstantDispersion)
    MolecularDiffusionCoefficent = GetConstReal(Material,"Molecular Diffusion Coefficent", ConstantDiffusion)
    CryogenicSuction = GetLogical(Material,"Compute Cryogenic Suction", Found)

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )   
    IF( Element % ElementIndex == 1 ) THEN
      CALL INFO(SolverName,'Number of Gauss points for 1st element:'&
          //TRIM(I2S(IP % n)),Level=31)
      CALL Info(SolverName,'Elemental n:'//TRIM(I2S(n))//' nd:'&
          //TRIM(I2S(nd))//' nd:'//TRIM(I2S(nb)),Level=31)
    END IF

    
    DO t=1,IP % n
      IPPerm = XiAtIPPerm(ElementID) + t
      
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:N) * LOAD(1:N) )

      ! from coordinate system
      Weight = IP % s(t) * DetJ

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
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
      IF (ComputeDeformation) THEN
        StressInvDtAtIP = ListGetElementReal( StressInvDt_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) &
             CALL WARN(SolverName,'"Stress Invariant Velocity" not found - setting to zero')
      END IF
             

      ! Variable gradients at IP
      gradTAtIP = ListGetElementRealGrad( Temperature_h,dBasisdx,Element,Found)
      IF (.NOT.Found) CALL FATAL(SolverName,'Unable to find Temperature gradient')
      gradYcAtIP = ListGetElementRealGrad( Salinity_h,dBasisdx,Element,Found)
      IF (.NOT.Found) CALL FATAL(SolverName,'Unable to find Salintiy gradient')
      gradpAtIP = ListGetElementRealGrad( Pressure_h,dBasisdx,Element,Found)
      IF (.NOT.Found) CALL FATAL(SolverName,'Unable to find Pressure gradient')

      ! Time derivatives of other variables
      IF (ComputeDt) THEN
        TemperatureDtAtIP = ListGetElementReal( TemperatureDt_h,Basis,Element, Found, GaussPoint=t)
        IF (.NOT.Found)  CALL FATAL(SolverName,'Temperature Velocity variable not found')
        SalinityDtAtIP = ListGetElementReal( SalinityDt_h,Basis,Element, Found, GaussPoint=t)
        IF (.NOT.Found)  CALL FATAL(SolverName,'Salinity Velocity variable not found')
      END IF


      ! Materialproperties (basically densities and derivatives of it)
      ! needed at IP for Xi computation (anything ELSE thereafter)
      rhosAtIP = rhos(RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      IF (ConstVal) THEN
        rhospAtIP = 0.0_dp
        rhosTAtIP = 0.0_dp
      ELSE
        rhospAtIP =rhosp(RockMaterialID,rhosAtIP,p0,PressureAtIP)
        rhosTAtIP =rhosT(RockMaterialID,rhosAtIP,T0,TemperatureAtIP)
      END IF
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      Xi0Tilde = GetXi0Tilde(RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson') ! classic simpified Anderson model
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
        CALL  GetXiHartikainen(RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.)
        IF (XiAtIP(IPPerm) .NE. XiAtIP(IPPerm)) THEN
          PRINT *, "Darcy: XiAtIP", B1AtIP,D1AtIP,Xi0Tilde
          PRINT *, "Darcy:  XiAtIP", deltaInElement, GlobalRockMaterial % e1(RockMaterialID), bijAtIP
          PRINT *, "Darcy:  XiAtIP", Xi0tilde,SalinityAtIP
          PRINT *, "Darcy: XiAtIP", B1AtIP*B1AtIP + D1AtIP !1.0/(1.0 + 0.5*B1AtIP + SQRT(B1AtIP*B1AtIP + D1AtIP)
          CALL FATAL(SolverName,"XiAtIP is NaN")
        END IF
      END SELECT

      ! on Xi (directly or indirectly) dependent material parameters (incl. updates) at IP
      rhowAtIP  = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal) ! update

      IF (ConstVal) THEN
        rhowPAtIP = 0.0_dp
        rhowTAtIP = 0.0_dp
        rhoiPAtIP = 0.0_dp
        rhoiTAtIP = 0.0_dp    
      ELSE
        rhowPAtIP = rhowP(CurrentSolventMaterial,rhowAtIP,p0,PressureAtIP) ! update with new rhowAtIP
        rhowTAtIP = rhowT(CurrentSolventMaterial,rhowAtIP,T0,TemperatureAtIP)
        rhoiPAtIP = rhoiP(CurrentSolventMaterial,rhoiAtIP,p0,PressureAtIP)
        rhoiTAtIP = rhoiT(CurrentSolventMaterial,rhoiAtIP,T0,TemperatureAtIP)
      END IF
      IF (NoSalinity) THEN
        rhocAtIP    = 0.0_dp
        rhocPAtIP   = 0.0_dp
        rhocTAtIP   = 0.0_dp
        rhocYcAtIP  = 0.0_dp
        rhowYcAtIP  = 0.0_dp
        rhogwYcAtIP = 0.0_dp      
      ELSE
        rhocAtIP    = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        rhocPAtIP   = rhocP(CurrentSoluteMaterial,rhocAtIP,ConstVal)
        rhocTAtIP   = rhocT(CurrentSoluteMaterial,rhocAtIP,T0,TemperatureAtIP,ConstVal)
        rhocYcAtIP  = rhocYc(CurrentSoluteMaterial,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
        rhowYcAtIP  = rhowYc(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP)
        rhogwYcAtIP = rhogwYc(rhowAtIP, rhocAtIP, rhowYcAtIP,rhocYcAtIP,XiAtIP(IPPerm),SalinityAtIP)
      END IF
      rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP(IPPerm),SalinityAtIP)
      IF (OffsetDensity) THEN

        IPPermRhogw = RhoOffsetAtIPPerm(ElementID) + t
        ! density of whole ground (rock + water + ice + solutes)
        rhoGAtIP = rhoG(rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP(IPPerm))
        IF (rhoGAtIP .NE. rhoGAtIP) THEN
          PRINT *,rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP
          CALL FATAL(FunctionName,'rhoGAtIP is NaN')
        END IF
        RhoOffsetAtIP(IPPermRhoGW) = rhoGAtIP - rhogwAtIP
      END IF


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ComputeFreshwaterHead) THEN
        IPPermFreshwaterHead = FreshwaterHeadAtIPPerm(ElementID) + t
        SELECT CASE(DIM)
        CASE(1)
          elevationAtIp = SUM( Basis(1:N)*Solver % Mesh % Nodes % x(1:N) )
        CASE(2)
          elevationAtIp = SUM( Basis(1:N)*Solver % Mesh % Nodes % y(1:N) )
        CASE DEFAULT
          elevationAtIp = SUM( Basis(1:N)*Solver % Mesh % Nodes % z(1:N) )
        END SELECT
        rhow0 = CurrentSolventMaterial % rhow0
        FreshwaterHeadAtIP(IPPermFreshwaterHead) = &
             (rhow0/rhogwAtIP)*(PressureAtIP/(rhow0*SQRT(SUM(gravity(1:DIM)*gravity(1:DIM)))) + elevationAtIp) &
             + ((rhogwAtIP - rhow0)/rhow0) * elevationAtIP
      END IF

      rhogwPAtIP = rhogwP(rhowPAtIP,rhocPAtIP,XiAtIP(IPPerm),SalinityAtIP)
      rhogwTAtIP = rhogwT(rhowTAtIP,rhocTAtIP,XiAtIP(IPPerm),SalinityAtIP)

      !IF ((rhogwAtIP .NE. rhogwAtIP) THEN ! sanity check
      !  PRINT *,"rhowgAtIP:",rhogwAtIP,XiAtIP(IPPerm),SalinityAtIP
      !  STOP
      !END IF

      ! conductivities at IP
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      KgwAtIP = GetKgw(RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP(IPPerm),MinKgw)
      KgwpTAtIP = 0.0_dp
      KgwppAtIP = 0.0_dp
      IF (CryogenicSuction) THEN
        fwAtIP = fw(RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP(IPPerm),GasConstant,TemperatureAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
      ELSE
        KgwppAtIP = KgwAtIP
        fwAtIP = 0.0_dp
        KgwpTAtIP = 0.0_dp
      END IF
 
      ! Elastic properties at IP
      EGAtIP = EG(CurrentSolventMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
      nuGAtIP = nuG(CurrentSolventMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
      kappaGAtIP = kappaG(EGAtIP,nuGAtIP)
      
      ! capacities at IP
      IF (HydroGeo) THEN   ! Simplifications: Xip=0 Xi=1 kappas=0
        CgwppAtIP = PorosityAtIP * rhogwPAtIP + rhogwAtIP * kappaGAtIP
      ELSE
        CgwppAtIP = GetCgwpp(rhogwAtIP,rhoiAtIP,rhosAtIP,rhogwPAtIP,rhoiPAtIP,rhosPAtIP,&
             kappaGAtIP,XiAtIP(IPPerm),XiPAtIP,&
             RockMaterialID,PorosityAtIP)
        !PRINT *,"CgwppAtIP:",CgwppAtIP," KgwppAtIP", KgwppAtIP
      END IF
      CgwpTAtIP = GetCgwpT(rhogwAtIP,rhoiAtIP,rhosAtIP,rhogwTAtIP,rhoiTAtIP,rhosTAtIP,XiAtIP(IPPerm),XiTAtIP,PorosityAtIP)
      IF (.NOT.NoSalinity) THEN
        CgwpYcAtIP = GetCgwpYc(rhogwAtIP,rhoiAtIP,rhogwYcAtIP,XiAtIP(IPPerm),XiYcAtIP,PorosityAtIP)     
      END IF
      IF (ComputeDeformation) THEN
        CgwpI1AtIP = GetCgwpI1(rhogwAtIP,rhoiAtIP,XiAtIP(IPPerm),&
             kappaGAtIP,RockMaterialID)
        IF (CgwpI1AtIP > 1.0d-03) THEN ! sanity check
          PRINT *,"CgwpI1AtIP", CgwpI1AtIP,&
               XiAtIP(IPPerm),rhogwAtIP,rhoiAtIP,kappaGAtIP, EGAtIP, nuGAtIP
          STOP
        END IF
      END IF

 

      
      !IF ( (.NOT.ConstantDispersion) .OR. FluxOutput) THEN
      IF (FluxOutput) THEN
         JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,&
             Gravity,rhogwAtIP,DIM,CryogenicSuction)
        !PRINT *, "JgwDAtIP", JgwDAtIP
        IF (FluxOutput) THEN
          GWfluxVar1 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(1)
          GWfluxVar2 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(2)
          IF (DIM == 3) GWfluxVar3 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(3)
        END IF
      END IF
      
      JcFAtIP = 0._dp
      IF (.NOT.NoSalinity) THEN  
        IF (ConstantDispersion) THEN
          KcAtIP = GetConstKC(DispersionCoefficient)
        ELSE
          IF (ConstantDiffusion) THEN
            DmAtIP = MolecularDiffusionCoefficent
          ELSE
            DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)
          END IF
          KcAtIP = GetKc(RockMaterialID,DmAtIP,XiAtIP(IPPerm),JgwDAtIP,PorosityAtIP)
        END IF
        ! parameters for diffusion-dispersion flow
        r12AtIP = GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,&
             rhowAtIP,rhocAtIP,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP)
        fcAtIP = GetFc(rhocAtIP,rhowAtIP,Gravity,r12AtIP,XiTAtIP,XiPAtIP,XiAtIP(IPPerm),gradPAtIP,gradTAtIP)
        KcYcYcAtIP = GetKcYcYc(KcAtIP,r12AtIP)
        JcFAtIP = GetJcF(KcYcYcAtIP,KcAtIP,fcAtIP,gradYcAtIP,SalinityAtIP)        
      END IF

      ! bedrock deformation velocity at IP
      vstarAtIP(1) = ListGetElementReal( Vstar1_h, Basis, Element, Found, GaussPoint=t)
      vstarAtIP(2) = ListGetElementReal( Vstar2_h, Basis, Element, Found, GaussPoint=t)
      IF (DIM > 2) &
           vstarAtIP(3) = ListGetElementReal( Vstar3_h, Basis, Element, Found, GaussPoint=t)
      
      ! fluxes other than pressure induced at IP
      DO i=1,DIM
        fluxTAtIP(i) =  SUM(KgwpTAtIP(i,1:DIM)*gradTAtIP(1:DIM))
        fluxgAtIP(i) = rhogwAtIP * SUM(KgwAtIP(i,1:DIM)*Gravity(1:DIM))   !!
        ! insert missing JcF here
        IF ((fluxgAtIP(i) .NE. fluxgAtIP(i)) .OR. (fluxTAtIP(i) .NE. fluxTAtIP(i))) THEN
          PRINT *, "NaN in r.h.s. of Darcy fluxes"
          PRINT *, "flux(",i,")= Jgwg",fluxgAtIP(i),"+ JgwpT", fluxTAtIP(i)
          PRINT *, "KgwAtIP=",KgwAtIP(i,1:DIM)
          PRINT *, "KgwpTAtIP=",KgwpTAtIP(i,1:DIM)
          PRINT *, "gradTAtIP(1:DIM):", gradTAtIP(1:DIM)
          IF (KgwpTAtIP(i,1) .NE. KgwpTAtIP(i,1)) PRINT *,CryogenicSuction,fwAtIP,XiPAtIP,KgwAtIP
          PRINT *, "rhowAtIP=",rhowAtIP," rhocAtIP=",rhocAtIP
          PRINT *, "XiAtIP(",IPPerm,")=",XiAtIP(IPPerm)
          STOP
        END IF
      END DO

      ! composition of the matrix:
      ! -----------------------------------
      ! this term can be switched off
      ! time derivative (Cgwpp*dp/dt,v):
      IF (ActiveMassMatrix) THEN
        DO p=1,nd
          DO q=1,nd
            MASS(p,q) = MASS(p,q) + Weight * CgwppAtIP * Basis(q) * Basis(p)
          END DO
        END DO
      END IF
      
      DO p=1,nd
        DO q=1,nd
          ! advection term due to bedrock velocity (Cgwpp * (vstar.grad(u)),v)
          ! --------------------------------------------------------------------
          STIFF (p,q) = STIFF(p,q) + Weight * &
             CgwppAtIP * SUM(vstarAtIP(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! diffusion term ( Kgwpp * grad(p),grad(v))
          ! div(J_gwp) = d(Kgwpp_i,j dp/dx_j)/dx_i
          ! -----------------------------------
          StiffPQ = 0.0
          DO i=1,DIM
            DO j=1,DIM
              StiffPQ = StiffPQ +  rhogwAtIP * KgwppAtIP(i,j) * dBasisdx(p,j)* dBasisdx(q,i)              
            END DO
          END DO
          STIFF(p,q) = STIFF(p,q) + Weight * StiffPQ
        END DO
      END DO
      
      ! body forces
      !--------------------------------------
      DO p=1,nd     
        FORCE(p) = FORCE(p) + Weight * rhogwAtIP *  SUM(fluxgAtIP(1:DIM)*dBasisdx(p,1:DIM))
        FORCE(p) = FORCE(p) - &
             Weight * Basis(p) * GlobalRockMaterial % etak(RockMaterialID) *&
             (rhocAtIP - rhowAtIP)* SUM(JcFAtIP(1:DIM)*dBasisdx(p,1:DIM))
        FORCE(p) = FORCE(p) + Weight * LoadAtIP * Basis(p)
        IF (CryogenicSuction) &
             FORCE(p) = FORCE(p) + Weight * rhogwAtIP * SUM(fluxTAtIP(1:DIM)*dBasisdx(p,1:DIM))
        IF (ComputeDt) THEN
          FORCE(p) = FORCE(p) + Weight * CgwpTAtIP * Basis(p)* TemperatureDtAtIP !dT/dt + v* grad T
          IF (.NOT. NoSalinity) &
               FORCE(p) = FORCE(p) + Weight * CgwpYcAtIP * Basis(p) * SalinityDtAtIP ! dyc/dt + v* grad yc
        END IF
        IF (ComputeDeformation)  THEN       
          FORCE(p) = FORCE(p) &
               - Weight*CgwpI1AtIP* Basis(p)* StressInvDtAtIP
          !IF (GetTimeStep() > 2) &
          !     PRINT *,"CgwpI1AtIP:",CgwpI1AtIP,"NodalStressInvDt(1:N):",NodalStressInvDt(1:N)
        END IF

      END DO
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixDarcy
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCDarcy( Element,  element_id, n, nd, CurrentSoluteMaterial)

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: n, nd, element_id
    TYPE(Element_t), POINTER :: Element
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    !------------------------------------------------------------------------------
    !REAL(KIND=dp) :: Flux(n), Coeff(n), Pressure(n), FluxAtIP, Weight
    REAL(KIND=dp) :: FluxAtIP, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,WeakPressure(nd), PressureCond(n)
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(nd)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: fluxgAtIP(3), rhogwAtIP, KgwAtIP(3,3), XiAtIP, mugwAtIP, MinKgw 

    REAL(KIND=dp) :: PressureAtIP, PorosityAtIP, SalinityAtIP, TemperatureAtIP, NormalAtIP(3)
    !REAL(KIND=dp), POINTER :: Nvector(:)
    LOGICAL :: Stat,Found,FluxCondition,WeakDirichletCond,ConstVal,ConstantsRead,Recharge
    INTEGER :: i,t,p,q,dim,body_id, other_body_id, material_id, RockMaterialID
    INTEGER, POINTER :: NodeIndexes(:)!, NPerm(:)
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition, ParentMaterial
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER ::  ParentElement
    !TYPE(Variable_t), POINTER :: NormalSolution
    
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlow (LocalMatrixBCDarcy)'

    SAVE Nodes,DIM,ConstantsRead,GasConstant, N0, DeltaT, T0, p0, eps, Gravity
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN
    
    FluxCondition = .FALSE.
    Recharge = .FALSE.

    ! just using WeakPressure as dummy destination to inquire whether we have a
    ! Dirichlet condition and can stop composing weakly imposed conditions
    WeakPressure(1:n) = GetReal( BoundaryCondition,TRIM(VarName), Found)
    IF (Found) THEN
      PressureCond(1:n) = 1.0_dp
      PressureCond(1:n) = GetReal(BoundaryCondition, TRIM(VarName)//" Condition", Found)
      IF (.NOT.(ANY(PressureCond <= 0.0_dp))) RETURN
    END IF
    
    WeakDirichletCond = .FALSE.

    WeakPressure(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)
   
    IF (WeakDirichletCond) THEN
      CALL INFO(FunctionName,'Setting weak condition (ignoring recharge and flux)',Level=12)
    ELSE
      FluxCondition = .TRUE.
    END IF
    
    IF(.NOT.ConstantsRead) THEN
      ConstantsRead = &
           ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      !PRINT *, "BCSolute: (Constantsread) ", GasConstant, N0, DeltaT, T0, p0, eps, Gravity, ConstantsRead
    END IF

    IF (FluxCondition) THEN
      body_id = GetInteger(BoundaryCondition,'Permafrost Target Body', Found)   
      ! inquire parent element and material
      IF (Found) THEN
        ParentElement => Element % BoundaryInfo % Right
        IF (body_id .NE. ParentElement % BodyId) THEN
          ParentElement => Element % BoundaryInfo % Left
        END IF                
      ELSE    
        other_body_id = Element % BoundaryInfo % outbody      
        IF (other_body_id < 1) THEN ! only one body in calculation
          ParentElement => Element % BoundaryInfo % Right
          IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => Element % BoundaryInfo % Left
        ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
          ParentElement => Element % BoundaryInfo % Right
          IF (ParentElement % BodyId == other_body_id) ParentElement => Element % BoundaryInfo % Left
        END IF
        body_id = ParentElement % BodyId
      END IF
        
      ! all the above was just so we can get the material properties of the parent element...      
      material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', Found)
      IF (.NOT.Found) CALL FATAL(FunctionName,'Parent Material ID in BC not found')
      !PRINT *,"Parent Material ID",  body_id,  material_id
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
        RockMaterialID = ParentElement % ElementIndex  ! each element has it's own set of parameters
      ELSE
        !RockMaterialID = ListGetInteger(ParentMaterial,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
        RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found)
        IF (.NOT.Found) THEN
          PRINT *,'ParentElement % ElementIndex',ParentElement % ElementIndex
          PRINT *,"Rock Material ID",RockMaterialID
        END IF
      END IF

      ConstVal = GetLogical(ParentMaterial,'Constant Permafrost Properties',Found)
      MinKgw = GetConstReal( Material, &
           'Hydraulic Conductivity Limit', Found)
      IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
           MinKgw = 1.0D-14
      
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

 
    

    ! Numerical integration:
    !-----------------------
    !IF (FluxCondition .OR. WeakDirichletCond) THEN ! spare us, if natural BC
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )
        Weight = IP % s(t) * DetJ
        !NodeIndexes => Element % NodeIndexes
        FluxAtIP = ListGetElementReal(Recharge_h, Basis, Element, Recharge)

        IF (Recharge) THEN
          FluxCondition = .TRUE.
          ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
          !TemperatureAtIP = ListGetElementReal( TemperatureBC_h, Basis=Basis, Element=Element, Found=Found, GaussPoint=t)
          !IF (FluxAtIP > 0.0) PRINT *, "Recharge=", FluxAtIP
          TemperatureAtIP = ListGetElementRealParent( Temperature_h, Basis=Basis, Element=Element, Found=Found )
          IF (.NOT.Found) CALL FATAL(SolverName,'Temperature in BC not found')
          PorosityAtIP = ListGetElementRealParent( Porosity_h, Basis=Basis, Element=Element, Found=Found)
          IF (.NOT.Found) CALL FATAL(SolverName,'Porosity in BC not found')
          PressureAtIP = ListGetElementRealParent( Pressure_h, Basis=Basis, Element=Element, Found=Found)
          !SalinityAtIP = ListGetElementReal( SalinityBC_h, Basis, Element, Found, GaussPoint=t)

          SalinityAtIP = 0.0_dp ! WE ASSUME FRESHWATER INFLOW!!!
          rhogwAtIP =  rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)! WE ASSUME FRESHWATER INFLOW!!!
          FluxAtIP = FluxAtIP*rhogwAtIP
        ELSE
          FluxCondition = .FALSE.
          FluxAtIP = ListGetElementReal(GWFlux_h, Basis, Element, FluxCondition)
          IF (FluxCondition) PRINT *, 'Groundwater Flux',  FluxAtIP
        END IF
        
       
        IF (Fluxcondition) THEN
          FORCE(1:nd) = FORCE(1:nd) + Weight * FluxAtIP * Basis(1:nd)
          !PRINT *,"FluxCondition", body_id, material_id
          ! Given pressure, weakly imposed
          !----------------------------------------------------------------------
        ELSE IF (WeakDirichletCond) THEN
          PressureAtIP = SUM(WeakPressure(1:n)*Basis(1:n))
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * PressureAtIP * Basis(1:nd)
        END IF
      END DO
      CALL DefaultUpdateEquations(STIFF,FORCE)
    !END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCDarcy
  !------------------------------------------------------------------------------

END SUBROUTINE PermafrostGroundwaterFlow
!------------------------------------------------------------------------------
