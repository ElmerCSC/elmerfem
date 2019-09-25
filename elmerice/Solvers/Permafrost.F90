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
!>  Solvers for enhanced permafrost problem 
!---------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!> Solver for groundwater flow of the enhanced permafrost model
!    (i.e. Darcy Flow representing saturated aquifer)
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

  CALL Info( SolverName, '-------------------------------------------',Level=1 )
  CALL Info( SolverName, '  Initializing Permafrost Groundwater Flow      ',Level=1 )
  CALL Info( SolverName, '-------------------------------------------',Level=1 )
  SolverParams => GetSolverParams()
  
  IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
    CALL ListAddString( SolverParams, 'Variable', 'GWPressure' )
    CALL WARN( SolverName, 'Variable not found. Adding default "GWPressure"')
  END IF

  ! Add linear system defaults: BCGStab+ILU0
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
       CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
       CALL ListAddString(SolverParams,'Linear System Iterative Method','BiCGStab')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
       CALL ListAddString(SolverParams,'Linear System Preconditioning','ILU0')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
       CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-08_dp)
  ! Add Nonlinear system defaults
  IF(.NOT. ListCheckPresent(SolverParams,'Nonlinear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0e-05_dp)
  IF(.NOT. ListCheckPresent(SolverParams,'Nonlinear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Nonlinear System Max Iterations',50) 
  
  CALL Info( SolverName, '  Done Initialization',Level=1)
  CALL Info( SolverName, '-------------------------------------------',Level=1 )
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
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
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
  !REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
  !     NodalPressure(:), DummyNodalGWflux(:,:), NodalStressInv(:),NodalDeformation(:),&
  !     NodalStressInvDt(:),NodalTemperatureDt(:), NodalSalinityDt(:),&
  !     NodalDummyDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.FALSE., NoSalinity=.FALSE.,GivenGWFlux,ElementWiseRockMaterial, DummyLog=.FALSE.,&
       InitializeSteadyState=.FALSE., ActiveMassMatrix=.TRUE., ComputeDeformation=.FALSE.,DeformationExists=.FALSE.,&
       StressInvAllocationsDone=.FALSE.,StressInvDtAllocationsDone=.FALSE.,&
       HydroGeo=.FALSE.,ComputeDt=.FALSE.,FluxOutput=.FALSE.,&
       TemperatureTimeDerExists=.FALSE.,SalinityTimeDerExists=.FALSE., OffsetDensity=.FALSE.
  CHARACTER :: DimensionString
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow'
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, StressInvName, &
       VarName,PhaseChangeModel,ElementRockMaterialName,DeformationName
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       TemperatureDt_h, SalinityDt_h, StressInv_h,StressInvDt_h,Vstar1_h, Vstar2_h, Vstar3_h
  
  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       ElementWiseRockMaterial, ComputeDeformation, FluxOutput,&
       StressInvAllocationsDone, StressInvDtAllocationsDone, OffsetDensity, &
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       TemperatureDt_h, SalinityDt_h, StressInv_h, StressInvDt_h, &
       Vstar1_h, Vstar2_h, Vstar3_h, &
       ActiveMassMatrix, InitializeSteadyState, HydroGeo, ComputeDt
       !  NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure,NodalStressInv, &
      ! NodalStressInvDt,NodalTemperatureDt,NodalDummyDt,NodalSalinityDt, DummyNodalGWflux, &
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
    IF (.NOT.Found) ComputeDt = .FALSE.
    IF (ComputeDt) THEN
      CALL INFO(SolverName,"Computing time derivatives in force vector",Level=1)
    ELSE
      CALL INFO(SolverName,"Ommitting time derivatives in force vector",Level=1)
    END IF
    ! inquire whether to compute deformation force term
    ComputeDeformation = GetLogical(Params,'Compute Deformation',Found)
    IF (ComputeDeformation) THEN
      CALL INFO(SolverName,"Including stress invariant derivative in force vector",Level=1)
    ELSE
      CALL INFO(SolverName,"Ommitting stress invariant derivative in force vector",Level=1)
    END IF
  END IF
  
  IF (InitializeSteadyState) THEN
    IF (GetTimeStep() == 1) THEN
      CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=1)
      ActiveMassMatrix = .FALSE.
    ELSE 
      CALL INFO(SolverName,"Switching mass matrix to active after initializing with steady state",Level=1)
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
      WRITE (DimensionString,'(I1)') I
      DummyGWfluxVar => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux '//TRIM(DimensionString))
      IF (.NOT.ASSOCIATED(DummyGWfluxVar)) THEN
        FluxOutput = .FALSE.
      ELSE
        FluxOutput = .TRUE.       
      END IF
      IF (.NOT.FluxOutput) EXIT
    END DO
    IF (FluxOutput) THEN
      WRITE (Message,*) 'Groundwater flow will be written to: Groundwater Flux {1..',DIM,'}'
      CALL INFO(SolverName,Message,Level=1)
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
      
      ! inquire whether to use hydro-geo simplifications
      HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
      IF (.NOT.Found) HydroGeo = .FALSE.
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
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          PRINT *, "NumberOfRockRecords=", NumberOfRockRecords
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        IF (.NOT.ASSOCIATED(CurrentSolventMaterial)) &
             CALL FATAL(Solvername,'Solvent Material not associated')
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
        IF (.NOT.ASSOCIATED(CurrentSoluteMaterial)) &
             CALL FATAL(Solvername,'Solute Material not associated')
      END IF
      IF (.NOT.ASSOCIATED(Material)) THEN
        WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
        CALL FATAL(SolverName,Message)
      END IF
      N  = GetElementNOFNodes()
      ND = GetElementNOFDOFs()
      NB = GetElementNOFBDOFs()

      ! compose element-wise contributions to matrix and R.H.S

      CALL LocalMatrixDarcy( Model, Element, t, N, ND+NB, Active,  &
           CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
           PhaseChangeModel,ElementWiseRockMaterial, ActiveMassMatrix, &
           ComputeDt,ComputeDeformation,HydroGeo, OffsetDensity,FluxOutput)
    END DO
    OffsetDensity = .FALSE. ! sure to not read after first time
    CALL DefaultFinishBulkAssembly()
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCDarcy(  Element, n, nd+nb )
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

CONTAINS
  ! PermafrostGroundWaterFlow : Assembly of the matrix entries arising from the bulk elements 

  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixDarcy( Model, Element, ElementID, n, nd, NoElements,&
       CurrentRockMaterial, CurrentSoluteMaterial,CurrentSolventMaterial,&
       PhaseChangeModel, ElementWiseRockMaterial, ActiveMassMatrix,&
       ComputeDt,ComputeDeformation,HydroGeo,OffsetDensity,FluxOutput)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    LOGICAL :: ElementWiseRockMaterial,ActiveMassMatrix,ComputeDt,&
         ComputeDeformation,HydroGeo,OffsetDensity,FluxOutput
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CgwppAtIP,CgwpTAtIP,CgwpYcAtIP,CgwpI1AtIP,KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),&
         meanfactor,MinKgw,gradTAtIP(3),gradPAtIP(3),gradYcAtIP(3),fluxTAtIP(3),fluxgAtIP(3),vstarAtIP(3) ! needed in equation
    REAL(KIND=dp) :: JgwDAtIP(3),JcFAtIP(3), DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3),&
         fcAtIP(3), DispersionCoefficient ! from salinity transport
    REAL(KIND=dp) :: Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT  XiAtIP,
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP,bijAtIP(2,2),bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: ks0th,e1,bs,rhos0,cs0,Xi0,eta0,Kgwh0(3,3),qexp,alphaL,alphaT,RadGen,acs(0:5),&
         as0,aas(0:5),ks0,cks(0:5)  ! stuff comming from RockMaterial
    INTEGER :: acsl,aasl,cksl       ! stuff comming from RockMaterial
    REAL(KIND=dp) :: EGAtIP,nuGAtIP,kappaGAtIP ! bedrock deformation
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0,eps,Gravity(3)! real constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,&
         rhospAtIP,rhosTAtIP,&
         rhowTAtIP,rhowPAtIP,rhowYcAtIP,&
         rhoiPAtIP,rhoiTAtIP,&
         rhocPAtIP,rhocTAtIP,rhocYcAtIP,&
         rhogwPAtIP,rhogwTAtIP,rhogwYcAtIP, rhoGAtIP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,StiffPQ
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,KPorosityAtIP,SalinityAtIP,PressureAtIP
    REAL(KIND=dp) :: TemperatureDtAtIP,SalinityDtAtIP,PressureDtAtIP,StressInvDtAtIP 
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp) , POINTER :: gWork(:,:)
    !REAL(KIND=dp) , ALLOCATABLE :: CgwpI1AtNodes(:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID, FluxDOFs,IPPerm,IPPermRhogw
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE., ConstVal=.FALSE., ConstantDispersion=.FALSE.,&
         CryogenicSuction=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostGroundWaterFlow', &
         FunctionName='Permafrost (LocalMatrixDarcy)'
    CHARACTER :: DimensionString
    TYPE(Variable_t), POINTER :: XiAtIPVar, RhoOffsetAtIPVar,GWfluxVar1,GWfluxVar2,GWfluxVar3
    INTEGER, POINTER :: XiAtIPPerm(:), RhoOffsetAtIPPerm(:),GWfluxPerm(:)
    REAL(KIND=dp), POINTER :: XiAtIP(:), RhoOffsetAtIP(:),GWFluxVal(:)
    
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


    IF (FluxOutput) THEN
      GWfluxVar1 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 1')
      GWfluxPerm => GWfluxVar1 % Perm
      GWfluxVar2 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 2')
      GWfluxVar2 % Perm = GWfluxPerm
      IF (DIM == 3) THEN
        GWfluxVar3 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 3')
        GWfluxVar3 % Perm = GWfluxPerm
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

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF
    DispersionCoefficient = GetConstReal(Material,"Dispersion Coefficient", ConstantDispersion)
    CryogenicSuction = GetLogical(Material,"Compute Cryogenic Suction", Found)
    IF (.NOT.Found) CryogenicSuction = .FALSE.


    ! check, whether we have globally or element-wise defined values of rock-material parameters
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
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
      SalinityAtIP = 0.0_dp
      SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')
      IF (ComputeDeformation) THEN
        StressInvDtAtIP = 0.0_dp
        StressInvDtAtIP =  ListGetElementReal( StressInvDt_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) &
             CALL WARN(SolverName,'"Stress Invariant Velocity" not found - setting to zero')
      END IF
             

      ! Variable gradients at IP
      gradTAtIP  = 0._dp
      gradYcAtIP = 0._dp
      gradpAtIP  = 0._dp
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
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      IF (ConstVal) THEN
        rhospAtIP = 0.0_dp
        rhosTAtIP = 0.0_dp
      ELSE
        rhospAtIP =rhosp(CurrentRockMaterial,RockMaterialID,rhosAtIP,p0,PressureAtIP)
        rhosTAtIP =rhosT(CurrentRockMaterial,RockMaterialID,rhosAtIP,T0,TemperatureAtIP)
      END IF
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)      
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson') ! classic simpified Anderson model
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)        
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.)
        IF (XiAtIP(IPPerm) .NE. XiAtIP(IPPerm)) THEN
          PRINT *, "Darcy: XiAtIP", B1AtIP,D1AtIP,Xi0Tilde
          PRINT *, "Darcy:  XiAtIP", deltaInElement, CurrentRockMaterial % e1(RockMaterialID), bijAtIP
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
        RhoOffsetAtIPVar => VariableGet( Model % Mesh % Variables, 'Reference Offset Density')
        IF (.NOT.ASSOCIATED(RhoOffsetAtIPVar)) THEN
          WRITE(Message,*) '"Permafrost Offset Density" is set, but variable "Reference Offset Density" is not associated'
          CALL FATAL(SolverName,Message)
        END IF
        RhoOffsetAtIPPerm => RhoOffsetAtIPVar % Perm  
        RhoOffsetAtIP => RhoOffsetAtIPVar % Values
        IPPermRhogw = RhoOffsetAtIPPerm(ElementID) + t
        rhoGAtIP = rhoG(rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP(IPPerm))
        IF (rhoGAtIP .NE. rhoGAtIP) THEN
          PRINT *,rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP
          CALL FATAL(FunctionName,'rhoGAtIP is NaN')
        END IF
        RhoOffsetAtIP(IPPermRhoGW) = rhoGAtIP - rhogwAtIP
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
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP(IPPerm),MinKgw)
      KgwpTAtIP = 0.0_dp
      KgwppAtIP = 0.0_dp
      IF (CryogenicSuction) THEN
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP(IPPerm),GasConstant,TemperatureAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
      ELSE
        KgwppAtIP = KgwAtIP
        fwAtIP = 0.0_dp
        KgwpTAtIP = 0.0_dp
      END IF
 

      ! Elastic properties at IP
      EGAtIP = EG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
      nuGAtIP = nuG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
      kappaGAtIP = kappaG(EGAtIP,nuGAtIP)
      
      ! capacities at IP
      IF (HydroGeo) THEN   ! Simplifications: Xip=0 Xi=1 kappas=0
        CgwppAtIP = PorosityAtIP * rhogwPAtIP + rhogwAtIP * kappaGAtIP
      ELSE
        CgwppAtIP = GetCgwpp(rhogwAtIP,rhoiAtIP,rhogwPAtIP,rhoiPAtIP,rhosPAtIP,&
             kappaGAtIP,XiAtIP(IPPerm),XiPAtIP,&
             CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        !PRINT *,"CgwppAtIP:",CgwppAtIP," KgwppAtIP", KgwppAtIP
      END IF
      CgwpTAtIP = GetCgwpT(rhogwAtIP,rhoiAtIP,rhogwTAtIP,rhoiTAtIP,rhosTAtIP,XiAtIP(IPPerm),XiTAtIP,PorosityAtIP)
      IF (.NOT.NoSalinity) THEN
        CgwpYcAtIP = GetCgwpYc(rhogwAtIP,rhoiAtIP,rhogwYcAtIP,XiAtIP(IPPerm),XiYcAtIP,PorosityAtIP)     
      END IF
      IF (ComputeDeformation) THEN
        CgwpI1AtIP = GetCgwpI1(rhogwAtIP,rhoiAtIP,XiAtIP(IPPerm),&
             kappaGAtIP,CurrentRockMaterial,RockMaterialID)
        IF (CgwpI1AtIP > 1.0d-03) THEN ! sanity check
          PRINT *,"CgwpI1AtIP", CgwpI1AtIP,&
               XiAtIP(IPPerm),rhogwAtIP,rhoiAtIP,kappaGAtIP, EGAtIP, nuGAtIP
          STOP
        END IF
      END IF

      ! parameters for diffusion-dispersion flow
      r12AtIP = GetR(CurrentSoluteMaterial,CurrentSolventMaterial,GasConstant,&
           rhowAtIP,rhocAtIP,XiAtIP(IPPerm),TemperatureAtIP,SalinityAtIP)
      DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,rhocAtIP,mugwAtIP,TemperatureAtIP)
      
      IF ( (.NOT.ConstantDispersion) .OR. FluxOutput) THEN
        JgwDAtIP = 0.0_dp
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,&
             Gravity,rhogwAtIP,DIM,CryogenicSuction)
        !PRINT *, "JgwDAtIP", JgwDAtIP
        IF (FluxOutput) THEN
          GWfluxVar1 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(1)
          GWfluxVar2 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(2)
          IF (DIM == 3) GWfluxVar2 % Values(GWfluxPerm(ElementID) + t) = JgwDAtIP(2)
        END IF
      END IF
      
      JcFAtIP = 0._dp
      IF (.NOT.NoSalinity) THEN  
        IF (ConstantDispersion) THEN
          KcAtIP = GetConstKC(DispersionCoefficient)
        ELSE
          KcAtIP = GetKc(CurrentRockMaterial,RockMaterialID,DmAtIP,XiAtIP(IPPerm),JgwDAtIP,PorosityAtIP)
        END IF  
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
             Weight * Basis(p) * CurrentRockMaterial % etak(RockMaterialID) *&
             (rhocAtIP - rhowAtIP)* SUM(JcFAtIP(1:DIM)*dBasisdx(p,1:DIM))
        FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
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
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixDarcy
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCDarcy( Element, n, nd)

    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Pressure(n), F,Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,PressureAtIP
    REAL(KIND=dp) :: MASS(nd,nd),STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), PARAMETER :: C=1000.0_dp
    LOGICAL :: Stat,Found,FluxCondition,WeakDirichletCond
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BoundaryCondition
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlow (LocalMatrixBCDarcy)'

    SAVE Nodes,DIM
    !------------------------------------------------------------------------------
    BoundaryCondition => GetBC()
    IF (.NOT.ASSOCIATED(BoundaryCondition) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BoundaryCondition,'Groundwater Flux', FluxCondition )
    ! Check, whether we have a weakly imposed Dirichlet condition
    Pressure(1:n) = GetReal( BoundaryCondition,'Imposed '// TRIM(VarName), WeakDirichletCond)

    ! Numerical integration:
    !-----------------------
    IF (FluxCondition .OR. WeakDirichletCond) THEN ! spare us, if natural BC
      IP = GaussPoints( Element )
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ
        ! Given flux:
        ! -----------
        IF (Fluxcondition) THEN

          F = SUM(Basis(1:n)*flux(1:n))
          FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
          ! Given pressure, weakly imposed
          !----------------------------------------------------------------------
        ELSE IF (WeakDirichletCond) THEN
          PressureAtIP = SUM(Pressure(1:n)*Basis(1:n))
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
            END DO
          END DO
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * PressureAtIP * Basis(1:nd)
        END IF
      END DO
      CALL DefaultUpdateEquations(STIFF,FORCE)
    END IF
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBCDarcy
  !------------------------------------------------------------------------------

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------  
END SUBROUTINE PermafrostGroundwaterFlow
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization for the primary solver, PermafrostGroundwaterFlux. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlux_Init( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: DT
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: dim
  CHARACTER(LEN=MAX_NAME_LEN) :: EqName, VarName, FluxName, GradName
  LOGICAL :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux_Init", &
       FluxVariableName
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  dim = CoordinateSystemDimension()
  CALL INFO('PermafrostGroundwaterFlux_init','---------------------------------------------------------',Level=1)
  CALL INFO('PermafrostGroundwaterFlux_init','Initializing computations for permafrost groundwater flow',Level=1)
    CALL INFO('PermafrostGroundwaterFlux_init','---------------------------------------------------------',Level=1)
  
  IF( dim < 2 .OR. dim > 3 ) THEN
    CALL Fatal('PermafrostGroundwaterFlux_init','Flux computation makes sense only in 2D and 3D')
  END IF


  VarName = TRIM('Permafrost GroundWater')
  !VarName = TRIM('GW')

  IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
    EqName = ListGetString( SolverParams,'Equation')
    CALL ListAddString( SolverParams, 'Variable','-nooutput '//TRIM(EqName)//'_temp' )
  END IF

  FluxName = TRIM(VarName)//' Flux'
  CALL Info('PermafrostGroundwaterFlux_init','Saving flux to: '//TRIM(FluxName), Level=1) 
  IF(dim == 2) THEN
    FluxVariableName=TRIM(FluxName)//'['//TRIM(FluxName)//':2]'
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':2]')
  ELSE IF(dim == 3) THEN
    FluxVariableName=TRIM(FluxName)//'['//TRIM(FluxName)//':3]'
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),&
         TRIM(FluxName)//'['//TRIM(FluxName)//':3]')
  ELSE
    CALL FATAL('PermafrostGroundwaterFlux_init','Wrong dimension of problem')
  END IF
  CALL ListAddString( SolverParams,&
       NextFreeKeyword('Exported Variable',SolverParams),&
       FluxVariableName)
  WRITE(Message,*) 'Added ',TRIM(FluxVariableName),' as variable'
  CALL INFO('PermafrostGroundwaterFlux_init',Message,Level=3)
  IF( GetLogical( SolverParams,'Calculate Flux Abs',GotIt) ) THEN
    FluxName = TRIM(VarName)//' Flux_abs'
    CALL Info('PermafrostGroundwaterFlux_init','Saving flux abs to: '//FluxName) 
    CALL ListAddString( SolverParams,&
         NextFreeKeyword('Exported Variable',SolverParams),TRIM(FluxName))
  END IF

  CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

  CALL ListAddLogical( SolverParams,'Skip Compute Nonlinear Change',.TRUE.)

  ! Add linear system defaults: cg+diagonal
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
       CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
       CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
       CALL ListAddString(SolverParams,'Linear System Preconditioning','diagonal')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
       CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-08_dp)

  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostGroundwaterFlux_Init
!> Solver for groundwater flow of the enhanced permafrost model
!------------------------------------------------------------------------------
SUBROUTINE PermafrostGroundwaterFlux( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

  USE CoordinateSystems
  USE DefUtils
  USE PermaFrostMaterials

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: Transient      !< Steady state or transient simulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName, PotName
  INTEGER :: i,j,k,dim,DOFs,firstmag
  LOGICAL :: GotIt
  REAL(KIND=dp) :: Unorm, Totnorm, val
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  TYPE(Variable_t), POINTER :: FluxSol
  TYPE FieldTable_t
     REAL(KIND=dp), POINTER :: Values(:)
     INTEGER, POINTER :: Perm(:)
  END TYPE FieldTable_t
  TYPE(FieldTable_t) :: Fields(3)
  TYPE(Variable_t), POINTER :: PressureVar,TemperatureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       DummyGWfluxVar
  INTEGER,POINTER :: PressurePerm(:), TemperaturePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       DummyGWfluxPerm(:)
  INTEGER :: NumberOfRockRecords
  REAL(KIND=dp),POINTER :: Pressure(:), Temperature(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       DummyGWflux(:)
!!$  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalTemperature(:), NodalSalinity(:),&
!!$       NodalPressure(:), NodalGWflux(:,:),NodalTemperatureDt(:),NodalPressureDt(:),&
!!$       NodalSalinityDt(:) ! all dummies
  LOGICAL :: FirstTime=.TRUE.,AllocationsDone,ConstantPorosity, NoSalinity,GivenGWFlux=.FALSE.,&
       UnfoundFatal=.TRUE.,ComputeDt=.FALSE.,DummyLog,   ComputeFluxAtIP = .FALSE., Found
  CHARACTER(LEN=MAX_NAME_LEN) :: TemperatureName, PorosityName, SalinityName, PressureName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName="PermafrostGroundwaterFlux"
  TYPE(ValueHandle_t) :: Temperature_h, Pressure_h, Salinity_h, Porosity_h
  ! -------------------------------------------------------------
  SAVE  SaveRHS,AllocationsDone,FirstTime,&
       Temperature_h, Pressure_h, Salinity_h, Porosity_h
  !NodalPorosity,NodalTemperature,NodalSalinity,NodalPressure, &
  !     NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,
  ! -------------------------------------------------------------

  IF (FirstTime) THEN
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )
    FirstTime=.FALSE.
  END IF


  CALL Info( SolverName, '-------------------------------------',Level=3 )
  CALL Info( SolverName, 'Computing the groundwater flux       ',Level=3 )
  CALL Info( SolverName, '-------------------------------------',Level=3 )

  dim = CoordinateSystemDimension()
  !------------------------------------------------------------------------------
  !  Check what needs to be computed
  !------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN

  SolverParams => GetSolverParams()
  Dofs = Dim

  ComputeFluxAtIP = GetLogical(SolverParams,'Compute Flux At IP',Found)
  IF (.NOT.Found) ComputeFluxAtIP = .FALSE.
  
  !ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)

  ! Read Variables
!!$  CALL AssignVars(Solver,Model,AllocationsDone,&
!!$       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
!!$       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
!!$       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
!!$       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
!!$       DummyGWfluxVar,DummyGWfluxVar,DummyGWfluxVar, &       
!!$       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
!!$       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
!!$       DummyGWfluxPerm, DummyGWfluxPerm,DummyGWfluxPerm, &
!!$       Temperature, Pressure, Porosity,Salinity,&
!!$       TemperatureDt, PressureDt, SalinityDt,&
!!$       DummyGWflux,DummyGWflux,DummyGWflux, &       
!!$       DummyLog, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt,SolverName)
  !-------------------------------------------------------------------------------
  ! If only one component is used use the scalar equation, otherwise use an
  ! auxiliary variable to store all the dimensions
  !-------------------------------------------------------------------------------
  Varname = TRIM('Permafrost Groundwater')

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 1',UnFoundFatal=UnFoundFatal )
  Fields(1) % Values => FluxSol % Values
  Fields(1) % Perm => FluxSol % Perm

  FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 2',UnFoundFatal=UnFoundFatal )
  Fields(2) % Values => FluxSol % Values
  Fields(2) % Perm => FluxSol % Perm

  IF( dim == 3 ) THEN
    FluxSol => VariableGet( Solver % Mesh % Variables, TRIM(VarName)//' Flux 3',UnFoundFatal=UnFoundFatal )
    Fields(3) % Values => FluxSol % Values
    Fields(3) % Perm => FluxSol % Perm
  END IF
  at0 = RealTime()

  IF (ComputeFluxAtIP)  THEN
    CALL FluxAtIP(Fields)
  ELSE
    CALL DefaultInitialize()

    ALLOCATE(ForceVector(SIZE(Solver % Matrix % RHS),DOFs))  
    ForceVector = 0.0_dp
    SaveRHS => Solver % Matrix % RHS

    CALL BulkAssembly()
    CALL DefaultFinishAssembly()

    at1 = RealTime()
    WRITE(Message,* ) 'Assembly Time: ',at1-at0
    CALL Info( SolverName, Message, Level=5 )
    !        
    !------------------------------------------------------------------------------     


    TotNorm = 0.0_dp
    DO i=1,Dofs
      WRITE(Message,'(A,I1,A,I1)') "Working on DOF ",i," out of ",Dofs
      CALL INFO(SolverName,Message,Level=3)
      Solver % Matrix % RHS => ForceVector(:,i)
      UNorm = DefaultSolve()
      WRITE( Message, * ) 'Norm of DOF: ',i,'=',UNorm ** 2.0_dp
      CALL INFO(SolverName,Message,Level=3)
      TotNorm = TotNorm + Unorm ** 2.0_dp
      Fields(i) % Values = Solver % Variable % Values
      !Fields(i) % Values = 1.0_dp * i
    END DO
    
    DEALLOCATE( ForceVector )  
    Solver % Matrix % RHS => SaveRHS
    TotNorm = SQRT(TotNorm)
    Solver % Variable % Norm = Totnorm


  !------------------------------------------------------------------------------     

    at2 = RealTime()
    WRITE(Message,* ) 'Solution Time: ',at2-at1
    CALL Info( SolverName, Message, Level=4 )
    
    WRITE( Message, * ) 'Result Norm: ',TotNorm
    CALL Info( SolverName, Message, Level=4 )
    
    CALL Info( SolverName, 'All done',Level=4 )
    CALL Info( SolverName, '-------------------------------------',Level=6 )
  END IF


CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE FluxAtIP(Fields)
    IMPLICIT NONE
    TYPE(FieldTable_t) :: Fields(3)
    !------------------------------------------------------------------------------
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: elem,t,i,j,k,p,q,n,nd, DIM,Rank, RockMaterialID, Active, IPPerm
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,detJ,GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    LOGICAL :: Found, ConstVal=.FALSE.,ConstantsRead=.FALSE., FirstTime=.TRUE., ElementWiseRockMaterial,&
         CryogenicSuction=.FALSE.
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: GasConstant, N0, meanfactor,DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         JgwDAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, &
         bijAtIP(2,2), bijYcAtIP(2,2),gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuf
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP  
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
!!$    REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:), NodalSalinity(:), NodalPressure(:),NodalPorosity(:)
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlux (BulkAssembly)'
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel,ElementRockMaterialName
       ! -------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, meanfactor,GasConstant, N0,DeltaT, eps,T0, p0,Gravity,&
         CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
         FirstTime,ElementWiseRockMaterial
    ! -------------------------------------------------------------
    
    n = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )

    !ALLOCATE( NodalPressure(N),NodalPorosity(N),NodalTemperature(N),NodalSalinity(N) )
    ALLOCATE( Basis(N), dBasisdx(N,3) )
    
    Active = Solver % NumberOFActiveElements
    DO elem = 1,Active

      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      Material => GetMaterial(Element)

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
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
          PRINT *, "NumberOfRockRecords", NumberOfRockRecords
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )        
      END IF

      IF(.NOT.ConstantsRead) THEN
        ConstantsRead = &
             ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      END IF

      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

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

      ! check, whether we have globally or element-wise defined values of rock-material parameters
      IF (ElementWiseRockMaterial) THEN
        RockMaterialID = elem  ! each element has it's own set of parameters
      ELSE
        RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
      END IF

      ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
      IF (.NOT.Found) THEN
        ConstVal = .FALSE.
      ELSE
        IF (ConstVal) &
             CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
      END IF

      deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
             IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

        ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
        TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
        PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
        PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
        SalinityAtIP = 0.0_dp
        SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')


        ! Variable gradients at IP
        gradTAtIP  = 0._dp
        gradPAtIP  = 0._dp
        gradTAtIP = ListGetElementRealGrad( Temperature_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Temperature gradient')
        gradpAtIP = ListGetElementRealGrad( Pressure_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Pressure gradient')

        !Materialproperties needed at IP
        rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        !        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        rhowAtIP =  rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiTAtIP = &
               XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiPAtIP   = &
               XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
        CASE DEFAULT ! Hartikainen model
          CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .TRUE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.)
        END SELECT
        rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = 0.0_dp
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = 0.0_dp
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        !PRINT *,"Flux:", rhogwAtIP,rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP
        Weight = IntegStuff % s(t) * detJ

 

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        DO i=1,DIM
          IPPerm = Fields(i)% Perm(elem) + t
          Fields(i)% Values(IPPerm) = JgwDAtIP(i)
        END DO
      END DO
    END DO
    DEALLOCATE(Basis, dBasisdx)
!!$    ,&
!!$         NodalPressure,NodalPorosity,NodalTemperature,NodalSalinity )

    !------------------------------------------------------------------------------
  END SUBROUTINE FluxAtIP
  !------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    INTEGER :: elem,t,i,j,k,p,q,n,nd, DIM,Rank, RockMaterialID, Active
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,detJ,GradAtIp(3)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    LOGICAL :: Found, ConstVal=.FALSE.,ConstantsRead=.FALSE., FirstTime=.TRUE., ElementWiseRockMaterial,&
         CryogenicSuction=.FALSE.
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: GasConstant, N0, meanfactor,DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: KgwAtIP(3,3),KgwppAtIP(3,3),KgwpTAtIP(3,3),MinKgw,gradTAtIP(3),gradPAtIP(3),&
         JgwDAtIP(3) ! needed in equation
    REAL(KIND=dp) :: XiAtIP,Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,ksthAtIP  ! function values needed for KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, &
         bijAtIP(2,2), bijYcAtIP(2,2),gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) :: fwAtIP, mugwAtIP !  JgwD stuf
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP  
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
!!$    REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:), NodalSalinity(:), NodalPressure(:),NodalPorosity(:)
    REAL(KIND=dp) :: TemperatureAtIP,PorosityAtIP,SalinityAtIP,PressureAtIP
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='PermafrostGroundwaterFlux (BulkAssembly)'
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel,ElementRockMaterialName
    ! -------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, meanfactor,GasConstant, N0,DeltaT, eps,T0, p0,Gravity,&
         CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
         FirstTime,ElementWiseRockMaterial
    ! -------------------------------------------------------------

    n = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )

    ALLOCATE( STIFF(n,n), FORCE(dofs,n) )
    ALLOCATE( Basis(n), dBasisdx(n,3) )
!!$    ALLOCATE( NodalPressure(N),NodalPorosity(N),NodalTemperature(N),NodalSalinity(N) )

    Active = Solver % NumberOFActiveElements
    DO elem = 1,Active

      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      Material => GetMaterial(Element)

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
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
          !PRINT *, "NumberOfRockRecords", NumberOfRockRecords
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF
        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
          FirstTime = .FALSE.
        END IF
        CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
        CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )        
      END IF

      IF(.NOT.ConstantsRead) THEN
        ConstantsRead = &
             ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
      END IF

      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()

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

      ! check, whether we have globally or element-wise defined values of rock-material parameters
      IF (ElementWiseRockMaterial) THEN
        RockMaterialID = elem  ! each element has it's own set of parameters
      ELSE
        RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
      END IF

      ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
      IF (.NOT.Found) THEN
        ConstVal = .FALSE.
      ELSE
        IF (ConstVal) &
             CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
      END IF

      deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
             IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )

        ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
        TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
        PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
        PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
        SalinityAtIP = 0.0_dp
        SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')

        ! Variable gradients at IP
        gradTAtIP  = 0._dp
        gradpAtIP  = 0._dp
        gradTAtIP = ListGetElementRealGrad( Temperature_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Temperature gradient')
        gradpAtIP = ListGetElementRealGrad( Pressure_h,dBasisdx,Element,Found)
        IF (.NOT.Found) CALL FATAL(SolverName,'Unable to compute Pressure gradient')

        !Materialproperties needed at IP
        rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        !        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        rhowAtIP =  rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiTAtIP = &
               XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          XiPAtIP   = &
               XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
        CASE DEFAULT ! Hartikainen model
          CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .TRUE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.)
        END SELECT
        rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
        rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0Tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = 0.0_dp
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = 0.0_dp
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        !PRINT *,"Flux:", rhogwAtIP,rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP
        Weight = IntegStuff % s(t) * detJ

        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
          END DO
        END DO

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        DO i=1,dim
          FORCE(i,1:nd) = FORCE(i,1:nd) + Weight *  JgwDAtIP(i) * Basis(1:nd)
        END DO
      END DO

      !------------------------------------------------------------------------------
      !      Update global matrices from local matrices 
      !------------------------------------------------------------------------------
      Solver % Matrix % Rhs => SaveRhs
      CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd) )
      !      END IF

      DO i=1,Dofs
        Solver % Matrix % RHS => ForceVector(:,i)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO

    END DO

    DEALLOCATE( STIFF, FORCE, Basis, dBasisdx)
!!$    ,&
!!$         NodalPressure,NodalPorosity,NodalTemperature,NodalSalinity )

    !------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostGroundwaterFlux
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!  Compute invariant of stress tensor and its time-derivative
!-----------------------------------------------------------------------------

SUBROUTINE PermafrostStressInvariant( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostStressInvariant'
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, StressVariableName, PressureName
  LOGICAL :: Found, FirstTime=.TRUE., NoPressure, UpdatePrev=.FALSE.,SteadyState=.FALSE.
  TYPE(Variable_t), POINTER :: InvariantVar, InvariantVeloVar, StressVariableVar, PressureVar
  INTEGER, POINTER :: InvariantPerm(:), InvariantVeloPerm(:), StressVariablePerm(:), PressurePerm(:)
  INTEGER :: I, DIM, StressVariableDOFs, CurrentTime, activenodes
  REAL (KIND=dp), POINTER :: Invariant(:), InvariantVelo(:), StressVariable(:),&
       InvariantPrev(:,:), Pressure(:), PrevPressure(:)
  REAL (KIND=dp) :: AverageInvariant
  !------------------------------------------------------------------------------
  SAVE FirstTime, CurrentTime, DIM, NoPressure,&
       PressureVar, Pressure, PressurePerm, PressureName,&
       StressVariableVar, StressVariable, StressVariablePerm, StressVariableDOFs, &
       StressVariableName, InvariantVeloVar, InvariantVeloPerm, InvariantVelo,SteadyState
 
  
  CALL Info( SolverName, '-----------------------------------------',Level=1 )
  CALL Info( SolverName, ' Computing  Permafrost Stress Invariant',Level=1 )
  CALL Info( SolverName, '-----------------------------------------',Level=1 )

  SolverParams => GetSolverParams()


  IF (FirstTime) THEN
    CurrentTime=GetTimeStep()
    UpdatePrev = .FALSE.
    DIM=CoordinateSystemDimension()
  ELSE
    IF (CurrentTime .NE. GetTimeStep()) THEN
      UpdatePrev = .TRUE.
      CurrentTime=GetTimeStep()
    ELSE
      UpdatePrev = .FALSE.
    END IF
  END IF

  IF (FirstTime .OR. Model % Mesh % Changed) THEN
    CALL INFO( SolverName, 'Initialization step:',Level=9 )
    SteadyState = GetLogical(SolverParams,'Steady State',Found)
    IF (SteadyState) THEN
      CALL INFO (SolverName,'Computing steady state only',Level=1)
      UpdatePrev = .FALSE.
    END IF
    PressureName = ListGetString(SolverParams, &
         'Pressure Variable', Found )
    IF (.NOT.Found) THEN
      CALL WARN(SolverName," 'Pressure Variable' not found. Using default 'Pressure' ")
      WRITE(PressureName,'(A)') 'Pressure'
    ELSE
      WRITE(Message,'(A,A)') "'Pressure Variable' found and set to: ", PressureName
      CALL INFO(SolverName,Message,Level=9)
    END IF
    PressureVar => VariableGet(Solver % Mesh % Variables,PressureName)
    IF (.NOT.ASSOCIATED(PressureVar)) THEN
      NULLIFY(Pressure)
      WRITE(Message,'(A,A,A)') 'Pressure Variable "', TRIM(PressureName), '" not associated'
      CALL FATAL(SolverName,Message)
    ELSE
      Pressure => PressureVar % Values
      IF (ASSOCIATED(Pressure)) THEN
        PressurePerm => PressureVar % Perm
        PrevPressure => PressureVar % PrevValues(:,1)
        NoPressure = .FALSE.
        WRITE(Message,'(A,A,A)') 'Pressure Variable "', TRIM(PressureName), '" associated'
        CALL INFO(SolverName,Message,Level=9)
      ELSE
        CALL FATAL(SolverName,'Pressure values not associated')
      END IF
    END IF
    CALL INFO( SolverName, 'Initialization completed',Level=9 )
  END IF
  
  StressVariableName = ListGetString(SolverParams,'Stress Variable Name',Found)
  IF (.NOT.Found) CALL FATAL(SolverName,' "Stress Variable Name" not found')
  StressVariableVar => VariableGet( Solver % Mesh % Variables, StressVariableName )
  IF ( ASSOCIATED(StressVariableVar)) THEN
    StressVariablePerm => StressVariableVar % Perm
    StressVariable => StressVariableVar % Values
    StressVariableDOFs = StressVariableVar % DOFs
  ELSE
    CALL FATAL(SolverName, TRIM(StressVariableName)//' not found')
  END IF


  ! END IF

    
  InvariantVar =>  Solver % Variable
  IF ( ASSOCIATED(InvariantVar)) THEN    
    VariableName = InvariantVar % Name
    Invariant => Solver % Variable % Values
    InvariantPerm => Solver % Variable % Perm
    InvariantPrev => Solver % Variable % PrevValues
    WRITE (Message,*) 'Solver variable ', TRIM(VariableName),', found'
    CALL INFO(SolverName,Message,Level=9)
  ELSE
    CALL FATAL(SolverName,'Solver variable not associated')
  END IF

  InvariantVeloVar => VariableGet( Solver % Mesh % Variables, TRIM(VariableName) // " Velocity" )
  IF ( ASSOCIATED(InvariantVeloVar)) THEN
    InvariantVeloPerm => InvariantVeloVar % Perm
    InvariantVelo =>InvariantVeloVar  % Values
     CALL INFO(SolverName, TRIM(VariableName) // " Velocity"//' found',Level=9)
  ELSE
    CALL FATAL(SolverName, TRIM(VariableName) // " Velocity"//' not found')
  END IF
  !!!!!!!!!!!!!
  !!!! VariableAdd for TRIM(VariableName) // " Velocity"
  
  activenodes = 0
  InvariantVelo = 0.0_dp  
  DO I = 1,Solver % Mesh % Nodes % NumberOfNodes
    IF (InvariantPerm(I) == 0) CYCLE
    IF (PressurePerm(I) == 0) THEN
      WRITE (Message,*) 'No entry for pressure variable',PressureName,'at point',I
      CALL FATAL(SolverName,Message)
    END IF
    IF ( UpdatePrev .AND. .NOT.(SteadyState) ) THEN
      InvariantPrev(InvariantPerm(I),1) =  Invariant(InvariantPerm(I))  
    END IF
    !PRINT *,"PermafrostStressInvariant:",Pressure(PressurePerm(I))
    !PRINT *,InvariantPerm(I)
    !PRINT *,Invariant(InvariantPerm(I))
    !PRINT *, (StressVariablePerm(I)-1)*StressVariableDOFs
    Invariant(InvariantPerm(I)) = &
         GetFirstInvariant(StressVariable,Pressure(PressurePerm(I)),(StressVariablePerm(I)-1)*StressVariableDOFs,DIM)
    IF (SteadyState) THEN
      InvariantVelo(InvariantVeloPerm(I)) = 0.0
    ELSE
      InvariantVelo(InvariantVeloPerm(I)) = &
           (Invariant(InvariantPerm(I))-InvariantPrev(InvariantPerm(I),1))/dt
      !PRINT *,"PermafrostStressInvariant:",InvariantVelo(InvariantVeloPerm(I))
    END IF
    !IF (InvariantVelo(InvariantVeloPerm(I)) /= 0.0_dp) THEN
    !  PRINT *, "InvariantVelo:", InvariantVelo(InvariantVeloPerm(I)), Invariant(InvariantPerm(I)), InvariantPrev(InvariantPerm(I),1)
    !END IF
    activenodes = activenodes + 1
    AverageInvariant = AverageInvariant + Invariant(InvariantPerm(I))
    IF (FirstTime .AND. .NOT.(SteadyState)) THEN
      InvariantPrev(InvariantPerm(I),1) = Invariant(InvariantPerm(I))
    END IF
  END DO
  AverageInvariant = AverageInvariant/DBLE(activenodes)
  WRITE(Message,*) 'Average invariant of ', activenodes,' out of ', Solver % Mesh % Nodes % NumberOfNodes,&
       ' active nodes:', AverageInvariant
  CALL INFO(SolverName,Message,Level=3)
  FirstTime = .FALSE.
  CONTAINS
    FUNCTION GetFirstInvariant(Stress,PressureAtPoint,Position,DIM) RESULT(FirstInvariant)
      REAL (KIND=dp) ::  FirstInvariant
      REAL (KIND=dp), POINTER :: Stress(:)
      REAL (KIND=dp) :: PressureAtPoint
      INTEGER :: DIM, Position
      !--------
      INTEGER :: I
      
      FirstInvariant = 0.0_dp
      ! tr(sigma - 1p)
      DO I=1,DIM
        FirstInvariant = FirstInvariant + Stress(Position+I) - PressureAtPoint
      END DO
      
    END FUNCTION GetFirstInvariant  
END SUBROUTINE PermafrostStressInvariant
  
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
  CHARACTER :: DimensionString
  LOGICAL :: OutputXi = .FALSE. , OutputFlux=.FALSE., Found
  TYPE(Variable_t), POINTER :: XiAtIPVar
  INTEGER, POINTER :: XiAtIPPerm(:)
  REAL(KIND=dp), POINTER :: XiAtIp(:)
  INTEGER :: I
  !------------------------------------------------------------------------------
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, '  Initializing heat transfer         ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  SolverParams => GetSolverParams()
  OutputXi = GetLogical(SolverParams, 'Output Xi', Found)
  IF (.NOT.Found) OutputXi = .FALSE.
  !PRINT *,SolverName,OutputXi
  IF (OutputXi) THEN
    !PRINT *,"Hello"
    CALL INFO(SolverName,'Output of IP variable "Xi" ',Level=1)
    CALL ListAddString( SolverParams, NextFreeKeyword('Exported Variable',SolverParams),'-IP -dofs 1 Xi' )
    CALL INFO(SolverName,'Added variable Xi',Level=1)
  ELSE
    CALL INFO(SolverName,'No output of IP variable "Xi" ',Level=1)
    CALL ListAddString( SolverParams, NextFreeKeyword('Exported Variable',SolverParams),'-IP -nooutput -dofs 1 Xi' )
    CALL INFO(SolverName,'Added variable Xi',Level=1)
  END IF
  
  ! Add linear system defaults: BiCGStab+ILU0
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
       CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
       CALL ListAddString(SolverParams,'Linear System Iterative Method','BiCGStab')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
       CALL ListAddString(SolverParams,'Linear System Preconditioning','ILU0')
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
       CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
  IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-08_dp)
  ! Add Nonlinear system defaults
  IF(.NOT. ListCheckPresent(SolverParams,'Nonlinear System Convergence Tolerance')) &
       CALL ListAddConstReal(SolverParams,'Nonlinear System Convergence Tolerance',1.0e-05_dp)
  IF(.NOT. ListCheckPresent(SolverParams,'Nonlinear System Max Iterations')) &
       CALL ListAddInteger(SolverParams,'Nonlinear System Max Iterations',50) 

  CALL Info( SolverName, ' Done Initializing      ',Level=1 )
END SUBROUTINE PermafrostHeatTransfer_init
!------------------------------------------------------------------------------
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
  TYPE(Variable_t), POINTER :: TemperatureVar,PressureVar,PorosityVar,SalinityVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       DepthVar, DummyGWfluxVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:), DepthPerm(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:), Depth(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., FluxOutput = .FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists=.FALSE.,&
       InitializeSteadyState=.FALSE.,ActiveMassMatrix=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostHeatEquation'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName, XiAtIPName
  CHARACTER :: DimensionString
  TYPE(ValueHandle_t) :: Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h,&
       Vstar1_h, Vstar2_h, Vstar3_h

  SAVE DIM,FirstTime,AllocationsDone,FluxOutput,DepthName,XiAtIPName,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,ComputeDt,DepthExists,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       PressureVelo_h, SalinityVelo_h, Depth_h, &
       Vstar1_h, Vstar2_h, Vstar3_h
  !------------------------------------------------------------------------------
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, 'Computing heat transfer              ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )

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
      CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=1)
      ActiveMassMatrix = .FALSE.
    ELSE 
      CALL INFO(SolverName,"Switching mass matrix to active after initializing with steady state",Level=1)
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
        WRITE (DimensionString,'(I1)') I
        DummyGWfluxVar => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux '//TRIM(DimensionString))
        IF (.NOT.ASSOCIATED(DummyGWfluxVar)) THEN
          FluxOutput = .FALSE.
        ELSE
          FluxOutput = .TRUE.       
        END IF
        IF (.NOT.FluxOutput) EXIT
      END DO
      IF (FluxOutput) THEN
        WRITE (Message,*) 'Groundwater flow will be read from: Groundwater Flux {1..',DIM,'}'
        CALL INFO(SolverName,Message,Level=1)
      END IF
    END IF
  END IF  
  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    WRITE(Message,*) "Nonlinear iteration ", iter, " out of ", maxiter
    CALL INFO( SolverName, Message, Level=3)
    
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
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
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

      CALL LocalMatrixHTEQ(  Element, t, Active, n, nd+nb,&
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
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
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial,&
       ActiveMassMatrix,FluxOutput)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
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
         CryogenicSuction=.FALSE.,HydroGeo=.FALSE.,ComputeFlux=.TRUE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER :: DimensionString
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
        GWfluxVar2 % Perm = GWfluxPerm
      END IF	
      IF (DIM == 3) THEN
        GWfluxVar3 => VariableGet( Solver % Mesh % Variables, 'Groundwater Flux 3')
        GWfluxVar3 % Perm = GWfluxPerm
      END IF
    END IF

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
    IF (.NOT.Found) HydroGeo = .FALSE.

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

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
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      LoadAtIP = 0.0_dp ! init      
      ! The heat soruce term
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
      !IF (LoadAtIP > 0.0_dp) PRINT *,"HTEQ:LoadAtIP", LoadAtIP      
      ! Contribution from Radiogenic Heat Production
      DepthAtIP = ListGetElementReal( Depth_h, Basis, Element, DepthExists, GaussPoint=t)
      IF (DepthExists) &
           RefDepth = GetConstReal(Material,'Radiogenic Reference Depth',DepthExists)
      IF (DepthExists)  &
        LoadAtIP = LoadAtIP  + RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,DepthAtIP,RefDepth)

      ! System variables (Temperature, Porosity, Pressure, Salinity) at IP
      PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
      PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
      IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
      SalinityAtIP = 0.0_dp
      SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
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
      vstarAtIP(1) = ListGetElementReal( Vstar1_h, Basis, Element, Found, GaussPoint=t)
      vstarAtIP(2) = ListGetElementReal( Vstar2_h, Basis, Element, Found, GaussPoint=t)
      IF (DIM > 2) &
           vstarAtIP(3) = ListGetElementReal( Vstar3_h, Basis, Element, Found, GaussPoint=t)

      !Materialproperties needed for computing Xi at IP
      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      IPPerm = XiAtIPPerm(ElementID) + t
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen (CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP(IPPerm),XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP(IPPerm),SalinityAtIP,ConstVal)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP(IPPerm),TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"HTEQ: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

      ! heat capacities
      csAtIP   = cs(CurrentRockMaterial,RockMaterialID,&
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
      ksthAtIP = GetKalphath(CurrentRockMaterial % ks0th(RockMaterialID),&
           CurrentRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
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
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP(IPPerm),MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
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
        DtdAtIP = GetDtd(CurrentRockMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP,JgwDAtIP)
        DO I=1,DIM
          DO J=1,DIM
            KGTTAtIP(I,J) = KGTTAtIP(I,J) + CGWTTAtIP * DtdAtIP(I,J)
          END DO
        END DO
      END IF

      Weight = IP % s(t) * DetJ
      !PRINT *, "Weight=", Weight
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
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
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

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostHeatTransfer
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!> solute (salt ions) transport equation for enhanced permafrost model
!-----------------------------------------------------------------------------
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
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE., ElementWiseRockMaterial, ActiveMassMatrix = .TRUE., &
       InitializeSteadyState = .FALSE.
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostSoluteTransport'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, VarName, TemperatureName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName
  TYPE(ValueHandle_t) :: Temperature_h, Pressure_h, Salinity_h, Porosity_h, Load_h
  
  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       ElementWiseRockMaterial,&
       Load_h, Temperature_h, Pressure_h, Salinity_h, Porosity_h,&
       ActiveMassMatrix, InitializeSteadyState
  !------------------------------------------------------------------------------
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, 'Computing solute transport           ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL DefaultStart()

  VarName = Solver % Variable % Name
  Params => GetSolverParams()

  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)

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
  IF (FirstTime) &
       InitializeSteadyState = GetLogical(Params,'Initialize Steady State',Found)
  IF (InitializeSteadyState) THEN
    IF (GetTimeStep() == 1) THEN
      CALL INFO(SolverName,"Initializing with steady state (no mass matrix)",Level=1)
      ActiveMassMatrix = .FALSE.
    ELSE 
      CALL INFO(SolverName,"Switching mass matrix to active after initializing with steady state",Level=1)
      ActiveMassMatrix = .TRUE.
      InitializeSteadyState = .FALSE.
    END IF
  END IF
  
  IF (FirstTime) THEN
    ! handles to local variables
    CALL ListInitElementKeyword( Temperature_h, 'Material', 'Temperature Variable' )
    CALL ListInitElementKeyword( Pressure_h, 'Material', 'Pressure Variable' )
    CALL ListInitElementKeyword( Salinity_h, 'Material', 'Salinity Variable' )
    CALL ListInitElementKeyword( Porosity_h, 'Material', 'Porosity Variable' )
      ! Handle to salinity Source (possible description of heat source at elements/IP's) 
    CALL ListInitElementKeyword( Load_h,'Body Force','Salinity Source' )
  END IF
  
  maxiter = ListGetInteger( Params,&
       'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    WRITE(Message,*) "Nonlinear iteration ", iter, " out of ", maxiter
    CALL INFO( SolverName, Message, Level=3)
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
          CALL INFO(SolverName,Message,Level=3)
          CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        END IF
        IF (ElementWiseRockMaterial) THEN
          ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
          NumberOfRockRecords = &
               ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        ELSE
          NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
        END IF

        IF (NumberOfRockRecords < 1) THEN
          CALL FATAL(SolverName,'No Rock Material specified')
        ELSE
          CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
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

      CALL LocalMatrixSolute(  Element, t, Active, n, nd+nb,&
           CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
           NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial,ActiveMassMatrix)
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()

    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBCSolute(  Element, t, Active, n, nd+nb,&
             CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
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
    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()
CONTAINS

  !------------------------------------------------------------------------------
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixSolute(  Element, ElementID, NoElements, n, nd,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, ActiveMassMatrix)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
!!$    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
!!$         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: ElementWiseRockMaterial,ActiveMassMatrix !GivenGWflux, 
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
         DmAtIP, r12AtIP(2), KcAtIP(3,3), KcYcYcAtIP(3,3), fcAtIP(3), extforceFlux(3), DispersionCoefficient ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor, DummyTensor(3,3)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n), XiBefore
    REAL(KIND=dp), POINTER :: gWork(:,:), XiAtIp(:)
    INTEGER :: i,j,t,p,q,DIM, RockMaterialID, IPPerm
    INTEGER, POINTER :: XiAtIPPerm(:)
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,ConstantDispersion,CryogenicSuction=.FALSE.
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
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

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

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)
    !      PRINT *,"Here0"
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
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
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
      !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP(IPPerm) = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp(IPPerm),0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        XiBefore =  XiAtIP(IPPerm)
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
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
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!

      ! capacities of solutes      
      CcYcTAtIP = CcYcT(rhocTAtIP,PorosityAtIP,SalinityAtIP)
      CcYcPAtIP = CcYcP(rhocPAtIP,PorosityAtIP, SalinityAtIP)
      CcYcYcAtIP = CcYcYc(rhocAtIP,rhocYcAtIP,PorosityAtIP, SalinityAtIP)

      ! groundwater viscosity is pulled outtside flux computation, as needed anyway
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)

      JgwDAtIP = 0.0_dp
      !PRINT *, "Solute: Compute Flux"
      mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
           XiAtIP(IPPerm),T0,SalinityAtIP,TemperatureAtIP,ConstVal)
      KgwAtIP = 0.0_dp
      KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
           mugwAtIP,XiAtIP(IPPerm),MinKgw)
      !PRINT *, "Solute: Kgw", KgwAtIP(1,1)
      fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
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
      DmAtIP = Dm(CurrentSoluteMaterial,N0,GasConstant,CurrentSoluteMaterial % rhoc0,&
           CurrentSolventMaterial % muw0,TemperatureAtIP)
      !PRINT *, "Solute: SalinityAtIP", SalinityAtIP
      !PRINT *, "Solute: Dm", DmAtIP,CurrentSoluteMaterial % rhoc0,CurrentSolventMaterial % muw0,TemperatureAtIP
      IF (ConstantDispersion) THEN
        KcAtIP = GetConstKC(DispersionCoefficient)
        !PRINT *,"DispersionCoefficient",KcAtIP
      ELSE
        DummyTensor = GetConstKC(3.565d-06)
        KcAtIP = GetKc(CurrentRockMaterial,RockMaterialID,DmAtIP,XiAtIP(IPPerm),JgwDAtIP,PorosityAtIP)
        !PRINT *,"Solute: Kc", KcAtIP(1,1), DummyTensor(1,1), DmAtIP,XiAtIP(IPPerm),JgwDAtIP(1:2),PorosityAtIP
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
          DO i=1,DIM
            extforceFlux =  SUM(KcAtIP(i,1:DIM)*fcAtIP(1:DIM))
          END DO
          STIFF (p,q) = STIFF(p,q) &
               - Weight * PorosityAtIP * rhocAtIP * Basis(q) * SUM(extforceFlux(1:DIM) * dBasisdx(p,1:DIM))

          ! time derivative (CcYcYc*du/dt,v):
          ! ------------------------------
          IF (ActiveMassMatrix) &
               MASS(p,q) = MASS(p,q) + Weight * CcYcYcAtIP  * Basis(q) * Basis(p)
        END DO
      END DO


      LoadAtIP = LoadAtIP !+ TemperatureTimeDer * CcYcTAtIP + PressureTimeDer * CcYcPAtIP 


      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixSolute
  !------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBCSolute(Element, ElementID, NoElements, n, nd,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    USE DefUtils
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
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
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(ParentMaterial,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(ParentMaterial,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

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
    IP = GaussPoints( Element )

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
        TemperatureAtIP = ListGetElementReal( Temperature_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Temperature not found')
        PorosityAtIP = ListGetElementReal( Porosity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Porosity not found')
        PressureAtIP = ListGetElementReal( Pressure_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL FATAL(SolverName,'Pressure not found')
        SalinityAtIP = 0.0_dp
        SalinityAtIP = ListGetElementReal( Salinity_h, Basis, Element, Found, GaussPoint=t)
        IF (.NOT.Found) CALL WARN(SolverName,'Salinity not found - setting to zero')

        !Materialproperties needed at IP

        ! water/ice densitities
        rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal) !!
        rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
        Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)
        !PRINT *,"Solute: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

        ! unfrozen pore-water content at IP
        SELECT CASE(PhaseChangeModel)
        CASE('anderson')
          XiAtIP = &
               GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
               CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
               T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
          ! NB: XiTAtIP, XiPAtIP not needed
        CASE DEFAULT ! Hartikainen model
          CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
               CurrentSoluteMaterial,CurrentSolventMaterial,&
               TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
               Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
               GasConstant,p0,T0,&
               XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
               .TRUE.,.FALSE., .FALSE., .FALSE.,.FALSE.) ! we need to compute, as IP's on boundary elements deviate from bulk
          ! NB: XiTAtIP, XiPAtIP, XiYcAtIP not needed
          PRINT *, "SoluteTransportBC", XiAtIP
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

  ! Perform static condensation in case bubble dofs are present
  !------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
    !------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    !------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
END SUBROUTINE PermafrostSoluteTransport

!==============================================================================
!>  initialization of IP variable to constant value
!==============================================================================
SUBROUTINE IPVariableInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: IPVar
  REAL(KIND=dp), POINTER :: IPVarValue(:)
  TYPE(ValueHandle_t) :: InitialIPVar_h
  REAL(KIND=dp) :: InitValue, detJ
  INTEGER :: IPVarDOFs, I, t, ICId, N, istat
  INTEGER, POINTER :: IPVarPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="IPVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: IPVariableName
  TYPE(GaussIntegrationPoints_t), TARGET :: IP
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  LOGICAL :: Visited = .FALSE., Found, ReadFromIC=.FALSE., stat
  REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:)
  SAVE Visited

  IF (Visited) RETURN

  SolverParams => GetSolverParams()

  IPVariableName = ListGetString(SolverParams, &
       'IP Variable', Found )
  IF (.NOT.Found) THEN
    CALL FATAL(SolverName, ' "IP Variable" not found - you have to provide one')
  ELSE
    WRITE (Message,*) ' "IP Variable ": ', TRIM(IPVariableName),' found' 
    CALL INFO(SolverName, Message,Level=1)
  END IF
  IPVar => VariableGet( Solver % Mesh % Variables, IPVariableName,Found,UnfoundFatal=.TRUE. )
  
  IF ( ASSOCIATED( IPVar ) ) THEN
    IPVarPerm    => IPVar % Perm
    IPVarValue  => IPVar % Values
    IPVarDOFs = IPVar % DOFs
    InitValue = GetConstReal(SolverParams,TRIM(IPVariableName),Found)
    ReadFromIC = .NOT.(Found)
  ELSE
    CALL FATAL(SolverName, 'Could not find "IP Variable"')
  END IF

  IF (ReadFromIC) THEN
    CALL ListInitElementKeyword( InitialIPVar_h,'Initial Condition',TRIM(IPVariableName) )
    WRITE(Message,*) IPVariableName, ' from corresponding initial condition'
    N = 2 * MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE(Basis(N),dBasisdx(N,3),stat=istat)
  ELSE
    WRITE(Message,*) IPVariableName, ' to constant',  InitValue
    IPVarValue = InitValue
  END IF
  
  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing ip variable           ', Level=1)
  CALL Info(SolverName, Message, Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  Visited = .TRUE.
  
  IF (ReadFromIC) THEN
    DO i = 1,  Solver % NumberOFActiveElements
      Element => GetActiveElement(i)
      IP = GaussPoints( Element )
      CALL GetElementNodes( Nodes )
      ICid = GetICId( Element, Found )
      IF (.NOT.Found) CALL FATAL(SolverName,'Corresponding "Initial Condition" not found')

      DO t=1,IP % n
        IF (ReadFromIC) THEN
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          InitValue = ListGetElementReal(InitialIPVar_h, Basis, Element, Found, GaussPoint=t)
          IF (.NOT.Found) CALL FATAL(SolverName,"Initial value not found in IC")
        END IF
        IPVarValue((IPVarPerm(i)*IPVarDOFs) + t*IPVarDOFs) = InitValue
        IPVarValue((IPVarPerm(i)*IPVarDOFs) + t*IPVarDOFs) = 1.0
      END DO
    END DO
    DEALLOCATE(Basis, dBasisdx)
  END IF
  CALL INFO(SolverName,"Itialisation Done",Level=1)
END SUBROUTINE IPVariableInit
!==============================================================================
!>  initialization of Porosity to given reference value in material
!==============================================================================
SUBROUTINE PorosityInit_old(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
  INTEGER, POINTER :: PorosityPerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:)
  REAL(KIND=dp), ALLOCATABLE :: NodalHits(:)
  INTEGER :: DIM, i, j, k, NumberOfRockRecords,RockMaterialID,CurrentNode,Active,totalunset,totalset
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PorosityInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,ElementRockMaterialName
  LOGICAL :: Visited = .False., Found, GotIt,ElementWiseRockMaterial

  SAVE Visited,ElementWiseRockMaterial,NodalHits
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------

  ! Execute solver only once at beginning
  IF (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing porosity to reference ', Level=1)
  CALL Info(SolverName, 'levels in material file            ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  ! Get variables
  DIM = CoordinateSystemDimension()

  ! Get info
  SolverParams => GetSolverParams()

  PorosityName = ListGetString(SolverParams, &
       'Porosity Variable', GotIt )
  IF (.NOT.GotIt) THEN
    PorosityName = "Porosity"
    CALL WARN(SolverName, ' "Porosity Variable" not found - trying default "Porosity"')
  END IF
  PorosityVariable => VariableGet( Solver % Mesh % Variables, PorosityName,GotIt )

  IF ( ASSOCIATED( PorosityVariable ) ) THEN
    PorosityPerm    => PorosityVariable % Perm
    PorosityValues  => PorosityVariable % Values
  ELSE
    CALL FATAL(SolverName, 'Could not find "Porosity Variable"')
  END IF

  ! Loop over elements
  Active = Solver % NumberOFActiveElements
  IF (.NOT.Visited) &
       ALLOCATE(NodalHits(Solver % Mesh % NumberOfNodes))

  NodalHits = 0.0_dp
  PorosityValues = 0.0_dp

  DO i = 1, Active
    CurrentElement => GetActiveElement(i)

    IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
    
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial(CurrentElement)
    
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (.NOT.Visited) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = &
           ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      !PRINT *,"PorosityInit:",TRIM(ElementRockMaterialName),ElementWiseRockMaterial
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
      END IF
      dim = CoordinateSystemDimension()
      Visited=.True.
    END IF

    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF

    ! Loop over nodes of element
    DO k = 1, GetElementNOFNodes(CurrentElement)
      CurrentNode = CurrentElement % NodeIndexes(k)!PorosityPerm(CurrentElement % NodeIndexes(k))
      NodalHits(CurrentNode) = NodalHits(CurrentNode) + 1.0_dp
      PorosityValues(CurrentNode) = &
           PorosityValues(CurrentNode) + CurrentRockMaterial % eta0(RockMaterialID)
      !PRINT *,CurrentNode,CurrentRockMaterial % eta0(RockMaterialID)
    END DO
  END DO

  totalset = 0
  totalunset = 0
  ! norm the result
  DO i = 1, Solver % Mesh % NumberOfNodes
    
    CurrentNode = PorosityPerm(i)
    IF ((NodalHits(i) > 0) .AND. (CurrentNode > 0)) THEN
      PorosityValues(CurrentNode) =  PorosityValues(CurrentNode)/(NodalHits(i))
      totalset = totalset + 1
    ELSE 
      PorosityValues(CurrentNode) =  0.0_dp
      totalunset = totalunset + 1
      IF (CurrentNode > 0) THEN
        WRITE(Message,*) 'Porosity value for active node ',CurrentNode,' has not been initiated' 
        CALL WARN(SolverName,Message)
      END IF
    END IF
  END DO

  WRITE(Message,*) 'Active elements:',Active,'. Initiated:',totalset,' of total meshpoints',  Solver % Mesh % NumberOfNodes

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Done Initializing porosity         ', Level=1)
  CALL Info(SolverName, Message, Level=1)
  IF (totalunset > 0) THEN
    WRITE(Message,*) totalunset, ' points in set'
    CALL WARN(SolverName,Message)
  END IF
  CALL Info(SolverName, '-----------------------------------', Level=1)
  !==============================================================================
END SUBROUTINE PorosityInit_Old

SUBROUTINE PorosityInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
  INTEGER, POINTER :: PorosityPerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:)
  REAL(KIND=dp), ALLOCATABLE :: NodalHits(:)
  INTEGER :: DIM, i, j, k, NumberOfRockRecords,RockMaterialID,CurrentNode,Active,totalunset,totalset
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PorosityInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,ElementRockMaterialName
  LOGICAL :: Visited = .False., Found, GotIt,ElementWiseRockMaterial

  SAVE Visited,ElementWiseRockMaterial,NodalHits
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------

  ! Execute solver only once at beginning
  IF (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing porosity to reference ', Level=1)
  CALL Info(SolverName, 'levels in material file            ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  ! Get variables
  DIM = CoordinateSystemDimension()

  ! Get info
  SolverParams => GetSolverParams()

  PorosityName = ListGetString(SolverParams, &
       'Porosity Variable', GotIt )
  IF (.NOT.GotIt) THEN
    PorosityName = "Porosity"
    CALL WARN(SolverName, ' "Porosity Variable" not found - trying default "Porosity"')
  END IF
  PorosityVariable => VariableGet( Solver % Mesh % Variables, PorosityName,GotIt )

  IF ( ASSOCIATED( PorosityVariable ) ) THEN
    PorosityPerm    => PorosityVariable % Perm
    PorosityValues  => PorosityVariable % Values
  ELSE
    CALL FATAL(SolverName, 'Could not find "Porosity Variable"')
  END IF
  !==============================================================================
  ! Loop over elements
  Active = Solver % NumberOFActiveElements
  PorosityValues = 0.0000001_dp ! some not complete insane default

  DO i = 1, Active
    CurrentElement => GetActiveElement(i)

    IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
    
    !NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial(CurrentElement)
    
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (.NOT.Visited) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = &
           ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      !PRINT *,"PorosityInit:",TRIM(ElementRockMaterialName),ElementWiseRockMaterial
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
      END IF
      dim = CoordinateSystemDimension()
      Visited=.True.
    END IF

    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF

    CurrentNode = PorosityPerm(i)
    !PRINT *,CurrentNode, i
    PorosityValues(CurrentNode) = CurrentRockMaterial % eta0(RockMaterialID)
  END DO
  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing porosity to reference ', Level=1)
  CALL Info(SolverName, 'done                               ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)
END SUBROUTINE PorosityInit
!==============================================================================
!>  Evolution of Porosity 
!==============================================================================
SUBROUTINE PermafrostPorosityEvolution( Model, Solver, Timestep, TransientSimulation )
 USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: PorosityVariable, TemperatureVar, PressureVar, StrainVar
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
  INTEGER, POINTER :: PorosityPerm(:), StrainPerm(:), TemperaturePerm(:), PressurePerm(:), NodeIndexes(:)
  REAL(KIND=dp), POINTER :: PorosityValues(:), Strain(:), Temperature(:), Pressure(:),&
       NodalStrain(:), NodalTemperature(:), NodalPressure(:),&
       PrevNodalTemperature(:), PrevNodalPressure(:),&
       PrevTemperature(:), PrevPressure(:)
  REAL(KIND=dp), ALLOCATABLE :: PrevStrainInvariant(:)
  REAL(KIND=dp) :: aux, Nodalrhos, PrevNodalrhos, StrainInvariant,&
            GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3)
  INTEGER :: DIM, i, j, k, N, NumberOfRockRecords,RockMaterialID,CurrentNode,Active,&
       StrainDOFs,TemperatureDOFS,PressureDOFs,totalunset,totalset,istat
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="PermafrostPorosityEvolution"
  CHARACTER(LEN=MAX_NAME_LEN) :: PorosityName,PressureName,TemperatureName,StrainVarName,ElementRockMaterialName
  LOGICAL :: FirstTime=.TRUE.,FirstVisit=.TRUE.,Found,GotIt,ElementWiseRockMaterial,&
       StrainVarExists, TemperatureVarExists, PressureVarExists,ConstVal
  !------------------------------
  SAVE FirstTime,FirstVisit,ElementWiseRockMaterial,&
       NodalStrain, NodalTemperature, NodalPressure,&
       PrevNodalTemperature, PrevNodalPressure,&
       StrainVar, TemperatureVar, PressureVar,&
       Strain, PrevStrainInvariant, Temperature, Pressure,&
       StrainDOFs,TemperatureDOFs,PressureDOFs,&
       PrevTemperature, PrevPressure, &
       StrainPerm, TemperaturePerm, PressurePerm,&
       StrainVarExists, TemperatureVarExists, PressureVarExists, &
       DIM, N0, GasConstant, DeltaT, T0, p0, eps, Gravity,&
       CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------
  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, ' computing evolution of porosity   ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)

  ! Get info and solver variable
  SolverParams => GetSolverParams()
  PorosityVariable => Solver % Variable
  IF (.NOT.ASSOCIATED(PorosityVariable)) CALL FATAL(SolverName,'No variable for "Porosity" associated')
  PorosityName = TRIM( Solver % Variable % Name )
  IF (TRIM(PorosityName) .NE. "porosity") THEN
    WRITE (Message,*) TRIM(PorosityName),' is not the expected "Porosity" - hopefully on purpose'
    CALL WARN(SolverName, Message)
  END IF

  PorosityPerm    => PorosityVariable % Perm
  PorosityValues  => PorosityVariable % Values
  IF (.NOT.ASSOCIATED( PorosityVariable % PrevValues )) THEN
    ALLOCATE(PorosityVariable % PrevValues(SIZE(PorosityValues), 1))
  END IF
  PorosityVariable % PrevValues(1:SIZE(PorosityValues),1) =&
       PorosityVariable % Values(1:SIZE(PorosityValues))
  
  IF(FirstTime) THEN
    IF(.NOT.(&
         ReadPermafrostConstants(Model, SolverName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)))&
         CALL FATAL(SolverName,"Errors in reading constants")
  END IF
  
  ! assign needed variables
  !(NB: we rather skip AssignVar routine, as ONLY Temperature and Pressure are needed)
  IF (FirstTime .OR. (Solver % Mesh % Changed ) ) THEN
    CALL INFO(SolverName,"Initialisation",Level=1)
    StrainVarName = GetString(SolverParams,'Strain Variable',Found)
    IF (.NOT.Found) THEN
      WRITE(StrainVarName,*) 'Strain'
      WRITE(Message,*) '"Strain Variable" not found - assuming default value: ',TRIM(StrainVarName)
      CALL WARN(SolverName,Message)
    ELSE
      WRITE(Message,*) '"Strain Variable" found and set to: ',TRIM(StrainVarName)
      CALL INFO(SolverName,Message,Level=1)
    END IF
    CALL AssignSingleVar(Solver,Model,NodalStrain,StrainVar,StrainPerm, Strain, &
         StrainVarName,StrainDOFs,StrainVarExists)
    IF (ASSOCIATED(StrainVar)) THEN
      IF (.NOT.FirstTime) DEALLOCATE(PrevStrainInvariant)
      ALLOCATE(PrevStrainInvariant(SIZE(StrainPerm)),stat=istat)
    ELSE
      CALL FATAL(SolverName,'No "Strain Varaible" associated')
    END IF
    TemperatureName = GetString(SolverParams,'Temperature Variable',Found)
    IF (.NOT.Found) THEN
      CALL WARN(SolverName,' "Temperature Variable" not found - assuming default value "Temperature" ')
      WRITE(TemperatureName,*) 'Temperature'
          ELSE
      WRITE(Message,*) '"Temperature Variable" found and set to: ',TRIM(TemperatureName)
      CALL INFO(SolverName,Message,Level=1)
    END IF
    CALL AssignSingleVar(Solver,Model,NodalTemperature,TemperatureVar,&
         TemperaturePerm, Temperature, &
         TemperatureName,TemperatureDOFS,TemperatureVarExists,&
         PrevNodalVariable=PrevNodalTemperature, PrevVariable=PrevTemperature)

    
    PressureName = GetString(SolverParams,'Pressure Variable',Found)
    IF (.NOT.Found) THEN
      CALL WARN(SolverName,' "Pressure Variable" not found - assuming default value "Pressure" ')
      WRITE(PressureName,*) 'Pressure'
    ELSE
      WRITE(Message,*) ' "Pressure Variable" found and set to: ',PressureName
      CALL INFO(SolverName,Message,Level=1)
    END IF
    CALL AssignSingleVar(Solver,Model,NodalPressure,PressureVar,&
         PressurePerm, Pressure, &
         PressureName,PressureDOFs,PressureVarExists,&
         PrevNodalVariable=PrevNodalPressure, PrevVariable=PrevPressure)

    FirstVisit = GetLogical(SolverParams,'Initialize Time Derivatives',Found)
    IF (.NOT.Found) FirstVisit = .FALSE.
  END IF
  ! some sanity check
  IF (.NOT.ASSOCIATED(Temperature) .OR. .NOT.ASSOCIATED(TemperaturePerm))&
       CALL FATAL(SolverName,'Values of temperature variable not found')
  IF (.NOT.ASSOCIATED(PrevTemperature))&
       CALL FATAL(SolverName,'Previous values of temperature variable not found')
  IF (.NOT.ASSOCIATED(Pressure) .OR. .NOT.ASSOCIATED(PressurePerm))&
       CALL FATAL(SolverName,'Values of pressure variable not found')
  IF (.NOT.ASSOCIATED(PrevPressure))&
       CALL FATAL(SolverName,'Previous values of pressure variable not found')
  ! Loop over elements
  Active = Solver % NumberOFActiveElements

  DO i = 1, Active
    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes
    Material => GetMaterial(CurrentElement)
    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(SolverName,'"Constant Permafrost Properties" set to true',Level=15)
    END IF
    IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material pointer found')
    IF (FirstTime) THEN
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = ListGetString(Material,"Element Rock Material File",ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
      END IF
      dim = CoordinateSystemDimension()
      FirstTime = .FALSE.
    END IF

    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = i
    ELSE      
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', GotIt,UnfoundFatal=.TRUE.)
      IF (.NOT.GotIt) CALL FATAL(SolverName,"Rock Material ID not found")
    END IF
    N = GetElementNOFNodes(CurrentElement)
    CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,NodalTemperature,Temperature,TemperatureDOFs)
    CALL ReadSingleVar(N,CurrentElement,PressurePerm,NodalPressure,Pressure,PressureDOFs)
    IF (FirstVisit) THEN ! write current values if starting
      CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,PrevNodalTemperature,Temperature,TemperatureDOFs)
      CALL ReadSingleVar(N,CurrentElement,PressurePerm,PrevNodalPressure,Pressure,PressureDOFs)
    ELSE
      CALL ReadSingleVar(N,CurrentElement,TemperaturePerm,PrevNodalTemperature,PrevTemperature,TemperatureDOFs)
      CALL ReadSingleVar(N,CurrentElement,PressurePerm,PrevNodalPressure,PrevPressure,PressureDOFs)
    END IF

    ! Loop over nodes of element
    DO k = 1, N
      CurrentNode = CurrentElement % NodeIndexes(k)
      PrevNodalrhos = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,&
           PrevNodalTemperature(k),&
           PrevNodalPressure(k),ConstVal)
      IF (PrevNodalrhos .NE. PrevNodalrhos) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for PrevNodalrhos"
        PRINT *,"PermafrostPorosityEvolution: ",&
           PrevNodalTemperature(k),&
           PrevNodalPressure(k),&
           RockMaterialID,T0,p0,ConstVal,&
           CurrentElement % NodeIndexes(k)
        CALL FATAL(SolverName,'Exiting')
      END IF
      Nodalrhos = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,&
           NodalTemperature(k),&
           NodalPressure(k),ConstVal)
      IF (Nodalrhos .NE. Nodalrhos) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for Nodalrhos"
        PRINT *,"PermafrostPorosityEvolution: ",&
           NodalTemperature(k),&
           NodalPressure(k),&
           RockMaterialID,T0,p0,ConstVal,&
           CurrentElement % NodeIndexes(k)
        CALL FATAL(SolverName,'Exiting')
      END IF
      StrainInvariant = 0.0
      ! first 1..DIM elements of StrainRate variable are th ediagonal entries
      DO J=1,DIM
        PrevStrainInvariant(StrainPerm(CurrentNode)) = StrainInvariant
        StrainInvariant = StrainInvariant &
             + Strain((StrainPerm(CurrentNode)-1)*StrainDOFs + J)
      END DO
      aux = 1.0_dp - (PrevNodalrhos/Nodalrhos)*&
           (1.0_dp + PrevStrainInvariant(StrainPerm(CurrentNode)))/(1.0_dp + StrainInvariant)
      PorosityValues(PorosityPerm(CurrentNode)) = &
           PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1)*(1.0_dp - aux) + aux
      IF (PorosityValues(PorosityPerm(CurrentNode)) <= 0.0) THEN
        WRITE(Message,*) "Reset Invalid value: ",&
             PorosityValues(PorosityPerm(CurrentNode)),&
             "=", PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1),&
             "*", (1.0_dp - aux), "+", aux
        !CALL FATAL(SolverName, Message)
        CALL WARN(SolverName, Message)
        PorosityValues(PorosityPerm(CurrentNode)) = 1.0d-08 !!!!! REPLACE 
      END IF
      IF (PorosityValues(PorosityPerm(CurrentNode)) .NE. PorosityValues(PorosityPerm(CurrentNode))) THEN
        PRINT *,"PermafrostPorosityEvolution: ","Found weird number for PorosityValues"
        PRINT *,"PermafrostPorosityEvolution: PorosityValues(",PorosityPerm(CurrentNode),")=",&
             PorosityValues(PorosityPerm(CurrentNode))
        PRINT *,"PermafrostPorosityEvolution:  Prev=",PorosityVariable % PrevValues(PorosityPerm(CurrentNode),1)
        PRINT *,"PermafrostPorosityEvolution: ",PrevNodalrhos,Nodalrhos,StrainInvariant
        CALL FATAL(SolverName,'Exiting')
      END IF
    END DO
  END DO
  FirstVisit = .FALSE.
END SUBROUTINE PermafrostPorosityEvolution

  
!==============================================================================
!>  initialization of arbitrary scalar given by nodal file
!==============================================================================
SUBROUTINE NodalVariableInit(Model, Solver, Timestep, TransientSimulation )
  !==============================================================================

  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: Timestep
  LOGICAL :: TransientSimulation

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: NodalVariable
  TYPE(ValueList_t), POINTER :: SolverParams,Material
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER, POINTER :: NodalVariablePerm(:)
  INTEGER,PARAMETER :: io=26
  INTEGER, ALLOCATABLE :: GlobalToLocalPerm(:)
  REAL(KIND=dp), POINTER :: NodalVariableValues(:)
  REAL(KIND=dp) :: InputField, InitValue, ValueOffset
  INTEGER :: DIM, i, j, CurrentNode, NumberOfNodes, MaxNumberOfGNodes, MinNumberOfGNodes,&
       OK,  counter, localGlobalRange
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName="NodalVariableInit"
  CHARACTER(LEN=MAX_NAME_LEN) :: NodalVariableName,NodalVariableFileName
  LOGICAL :: Visited = .FALSE., Found, Parallel, GotIt, FromFile=.FALSE.

  !SAVE Visited
  !,DIM,CurrentRockMaterial,NumberOfRockRecords

  !------------------------------------------------------------------------------

  ! Execute solver only once at beginning
  !if (Visited) RETURN

  CALL Info(SolverName, '-----------------------------------', Level=1)
  CALL Info(SolverName, 'Initializing variable to reference ', Level=1)
  CALL Info(SolverName, 'levels (either file or IC)         ', Level=1)
  CALL Info(SolverName, '-----------------------------------', Level=1)


  DIM = CoordinateSystemDimension()
  Parallel = (ParEnv % PEs > 1)
  Mesh => GetMesh()

  ! Get variable to fill in
  SolverParams => GetSolverParams()

  NodalVariableName = ListGetString(SolverParams, &
       'Nodal Variable', GotIt )
  IF (.NOT.GotIt) THEN
    CALL FATAL(SolverName, ' "Nodal Variable" not found')
  END IF
  NodalVariable => VariableGet( Mesh % Variables, NodalVariableName,GotIt )
  IF (.NOT.GotIt) CALL FATAL(SolverName,"Variable not found")

  IF ( ASSOCIATED( NodalVariable ) ) THEN
    NodalVariablePerm    => NodalVariable % Perm
    NodalVariableValues  => NodalVariable % Values
    WRITE (Message,*) 'Reading variable ',TRIM(NodalVariableName)
    CALL INFO(SolverName,Message,Level=3)
  ELSE
    WRITE (Message,*) 'Could not find ',TRIM(NodalVariableName)
    CALL FATAL(SolverName, Message)
  END IF
  NodalVariableValues = 0.0_dp

  NodalVariableFileName = ListGetString(SolverParams, &
       'Nodal Variable File', FromFile )

  ValueOffset = GetConstReal(SolverParams,'Variable Offset',GotIt)
  IF (.NOT.GotIt) THEN
    ValueOffset = 0.0_dp
  ELSE
    WRITE (Message,*) ' "Variable Offset" found and set to: ', ValueOffset
    CALL INFO(SolverName,Message,Level=3)
  END IF

  IF (.NOT.FromFile) THEN
    InitValue = GetConstReal(SolverParams,TRIM(NodalVariableName),Found)
    IF (.NOT.Found) THEN
      WRITE(Message,*) 'No entry for ',TRIM(NodalVariableName),&
           ' found in Solver section (IC version not implemented)'
      CALL FATAL(SolverName,Message)
    END IF
    NodalVariableValues = InitValue + ValueOffset
  ELSE
    NumberOfNodes = Mesh % NumberOfNodes

    IF (Parallel) THEN
      MaxNumberOfGNodes = MAXVAL(Mesh % ParallelInfo % GlobalDOFs)
      MinNumberOfGNodes = MINVAL(Mesh % ParallelInfo % GlobalDOFs)
      !localGlobalRange = MaxNumberOfGNodes - MinNumberOfGNodes
      IF (MaxNumberOfGNodes <= MinNumberOfGNodes) CALL FATAL(SolverName,"No nodes in parallel domain")
      ALLOCATE(GlobalToLocalPerm(MinNumberOfGNodes:MaxNumberOfGNodes), STAT=OK)
      IF (OK /= 0) CALL FATAL(SolverName,"Allocation error of GlobalToLocalPerm")
      GlobalToLocalPerm = 0
      DO I=1,NumberOfNodes
        GlobalToLocalPerm(Mesh % ParallelInfo % GlobalDOFs(I)) = I
      END DO
      PRINT *, TRIM(SolverName),": ParENV:",ParEnv % MyPE,".  Global Nodal Numbers from",&
           MinNumberOfGNodes,"to",MaxNumberOfGNodes
    ELSE
      MinNumberOfGNodes = 1
      MaxNumberOfGNodes = NumberOfNodes   
    END IF

    OPEN(unit = io, file = TRIM(NodalVariableFileName), status = 'old',action='read',iostat = ok)
    IF (ok /= 0) THEN
      WRITE(Message,'(A,A)') 'Unable to open file ',TRIM(NodalVariableFileName)
      CALL FATAL(TRIM(SolverName),TRIM(message))
    ELSE
      !------------------------------------------------------------------------------
      ! Read in the number of records ordered in global node-numbering
      ! in file (first line integer)
      !------------------------------------------------------------------------------
      DO J=1,MaxNumberOfGNodes ! all or in parallel up to max global index
        READ (io, *, END=70, IOSTAT=OK, ERR=80) counter, InputField
        IF (counter .NE. J) CALL FATAL(SolverName,'No concecutive numbering in file')
        IF (J < MinNumberOfGNodes) CYCLE

        IF (Parallel) THEN
          I = GlobalToLocalPerm(J)
          IF (I == 0) CYCLE ! point in range, but not in partition        
        ELSE
          I=J
        END IF
        !IF ((NodalVariablePerm(I)<1) .OR. (NodalVariablePerm(I)>NumberOfNodes)) THEN
        !  PRINT *, "NodalVariableInit:", ParEnv % myPE, "NodalVariablePerm(",I,")=",&
        !       NodalVariablePerm(I),">",NumberOfNodes
        !  CALL FATAL(SolverName,'No corresponding entry of target variable')
        !END IF
        NodalVariableValues(NodalVariablePerm(I)) = InputField + ValueOffset
        ! PRINT *,i,counter
      END DO
      !PRINT *, "END", i,counter
70    IF (J-1 .NE. MaxNumberOfGNodes) THEN
        WRITE (Message,*) 'Number of records ',i,' in file ',&
             TRIM(NodalVariableFileName),' does not match number of nodes ',&
             NumberOfNodes, ' in mesh'
        CALL FATAL(SolverName,Message)
      END IF
      CLOSE (io)
      IF (Parallel) &
           DEALLOCATE(GlobalToLocalPerm)
      RETURN
80    CALL FATAL(SolverName,"I/O error")
    END IF
  END IF
END SUBROUTINE NodalVariableInit

!==============================================================================
!> output of material parameter
!==============================================================================
SUBROUTINE PermafrostMaterialOutput( Model,Solver,dt,TransientSimulation )
  !==============================================================================
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
  TYPE(Variable_t), POINTER :: PressureVar,PorosityVar,SalinityVar,TemperatureVar,&
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       GWFluxVar1,GWFluxVar2,GWFluxVar3
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,TensorComponent(2)
  INTEGER,PARAMETER :: io=24
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),PorosityPerm(:),SalinityPerm(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       OutputPropertyPerm(:),GWfluxPerm1(:),GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       OutputProperty(:),GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:), NodalGwFlux(:,:),&
       NodalTemperatureDt(:),NodalPressureDt(:),NodalSalinityDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE., &
       ComputeDt=.FALSE.,ComputeXiT=.FALSE.,ElementWiseRockMaterial,GivenGWflux
  !CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostMaterialOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, TemperatureName,&
       PhaseChangeModel, VariableName,ElementRockMaterialName

  SAVE DIM,FirstTime,AllocationsDone,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       NumberOfRockRecords,NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGwFlux,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       ElementWiseRockMaterial
  !------------------------------------------------------------------------------
  Params => GetSolverParams()
  VariableName = ListGetString(Params,"Output Property",Found)
  ComputeDt = GetLogical(Params,'Compute Time Derivatives',Found)

  IF (.NOT.Found) &
       CALL INFO(SolverName, ' "Output Property" not found')

  TensorComponent(1:2) = ListGetInteger(Params,"Output Property Component",Found)
  IF (.NOT.Found) TensorComponent = 1

  CALL DefaultInitialize()

  ! Assign output variables
  OutputProperty => Solver % Variable % Values
  OutputPropertyPerm => Solver % Variable % Perm

  ! Read Variables
  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar, &
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux,DIM,ComputeDt,SolverName)

  Active = GetNOFActive()


  DO t=1,Active
    Element => GetActiveElement(t)      
    n  = GetElementNOFNodes(Element)
    Material => GetMaterial(Element)
    PhaseChangeModel = ListGetString(Material, &
         'Permafrost Phase Change Model', Found )
    IF (Found) THEN
      WRITE (Message,'(A,A)') '"Permafrost Phase Change Model" set to ', TRIM(PhaseChangeModel)
      CALL INFO(SolverName,Message,Level=9)
    END IF

    IF (FirstTime) THEN
      dim = CoordinateSystemDimension()
      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
        PRINT *, "NumberOfRockRecords", NumberOfRockRecords
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF
      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
        FirstTime = .FALSE.
      END IF
      CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
      CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )        
    END IF

    CALL ReadVars(N,Element,Model,Material,&
         NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
         Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
         TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
         GWfluxPerm1,GWfluxPerm2,GWfluxPerm3,&
         NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
         PorosityName,SolverName,DIM)

    CALL LocalMatrixMaterialOutput( Element, t, Active, n, nd,&
         NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
         NodalGWflux, GivenGWflux,&
         CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
         NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, &
         VariableName,TensorComponent)
  END DO

  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()
  IF (TRIM(VariableName) == 'kgw') THEN
    DO I=1,Solver % Mesh % NumberOfNodes
      OutputProperty(OutputPropertyPerm(i)) = EXP(OutputProperty(OutputPropertyPerm(i)))
    END DO
  END IF
  CALL INFO("SolverName","Computation of " // TRIM(VariableName) // "done",Level=1)

CONTAINS
  ! Assembly of the matrix entries arising from the bulk elements
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixMaterialOutput( Element, ElementID, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, GivenGWflux,&
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial, &
       VariableName,TensorComponent)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords, TensorComponent(2)
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel

    !    INTEGER :: n,TensorComponent(2)
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
         ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,hiAtIP,hwAtIP  ! function values needed for C's and KGTT
    REAL(KIND=dp) :: B1AtIP,B2AtIP,DeltaGAtIP, bijAtIP(2,2), bijYcAtIP(2,2),&
         gwaAtIP,giaAtIP,gwaTAtIP,giaTAtIP,gwapAtIP,giapAtIP !needed by XI
    REAL(KIND=dp) ::  gradTAtIP(3),gradPAtIP(3),JgwDAtIP(3),KgwAtIP(3,3),KgwpTAtIP(3,3),MinKgw,KgwppAtIP(3,3),fwAtIP,mugwAtIP!  JgwD stuff
    REAL(KIND=dp) :: deltaInElement,D1AtIP,D2AtIP
    REAL(KIND=dp) :: GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3) ! constants read only once
    REAL(KIND=dp) :: rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,rhogwAtIP,csAtIP,cwAtIP,ciAtIP,ccAtIP ! material properties at IP
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,Weight,LoadAtIP,&
         TemperatureAtIP,PorosityAtIP,PressureAtIP,SalinityAtIP,&
         StiffPQ, meanfactor
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    REAL(KIND=dp), POINTER :: gWork(:,:)
    REAL(KIND=DP) :: PropertyAtIP
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,ConstantMeanFactor=.FALSE.,&
         CryogenicSuction=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName, VariableName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixMaterialOutput)'
    !------------------------------------------------------------------------------
    SAVE Nodes, ConstantsRead, DIM, GasConstant, N0,DeltaT, T0, p0,eps,Gravity
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

    ! Get stuff from SIF BodyForce section
    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
         LOAD(1:n) = GetReal( BodyForce,'Heat Source', Found )

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

    meanfactor = GetConstReal(Material,"Conductivity Arithmetic Mean Weight",ConstantMeanFactor)
    
    !IF (.NOT.Found) THEN
    !  CALL INFO(FunctionName,'"Conductivity Arithmetic Mean Weight" not found. Using default unity value.',Level=9)
    !  meanfactor = 1.0_dp
    !END IF
    MinKgw = GetConstReal( Material, &
         'Hydraulic Conductivity Limit', Found)
    IF (.NOT.Found .OR. (MinKgw <= 0.0_dp))  &
         MinKgw = 1.0D-14

    deltaInElement = delta(CurrentSolventMaterial,eps,DeltaT,T0,GasConstant)

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))
      
      
      !Materialproperties needed at IP

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal) !!
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen(CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.)
      END SELECT
      
      IF (.NOT.ConstantMeanFactor) & 
           meanfactor = 1.0_dp - PorosityAtIP * XiAtIp
      
      !Materialproperties needed at IP:
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)

      ! heat capacities
      csAtIP   = cs(CurrentRockMaterial,RockMaterialID,&
           T0,TemperatureAtIP,ConstVal)
      cwAtIP   = cw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
      ciAtIP   = ci(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      ccAtIP   = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)

      ! latent heat
      hiAtIP = hi(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      hwAtIP = hw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)

      ! heat conductivity at IP
      ksthAtIP = GetKalphath(CurrentRockMaterial % ks0th(RockMaterialID),&
           CurrentRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
      kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)      
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP,&
           SalinityATIP,PorosityAtIP,meanfactor)

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
      CGTpAtIP = GetCGTp(rhoiAtIP,hiAtIP,hwAtIP,XiPAtIP,PorosityAtIP) !NEW
      CGTycAtIP = GetCGTyc(rhoiAtIP,hiAtIP,hwAtIP,XiYcAtIP,PorosityAtIP) !NEW

      ! groundwater flux at IP
      JgwDAtIP = 0.0_dp
      IF (GivenGWFlux) THEN
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N)) 
        END DO
      ELSE        
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO
        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
      END IF

      ! select parameter name for output
      SELECT CASE(VariableName)
      CASE('kgw')
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP =  GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        IF (KgwAtIP(TensorComponent(1),TensorComponent(2)) <= 0.0_dp) STOP
        PropertyAtIP = KgwAtIP(TensorComponent(1),TensorComponent(2))
        !IF (PropertyAtIP < 0.0_dp) THEN
        !  PRINT *, "KgwAtIP < 0 =", KgwAtIP, muw0,muw0,XiAtIP,rhow0,qexp,Kgwh0,MinKgw!
        ! STOP
        ! ELSE
        !   PRINT *, "KgwAtIP=", PropertyAtIP
        ! END IF
      CASE('rhogw')       
        PropertyAtIP = rhogwAtIP
      CASE('gwa')
        PropertyAtIP = gwaAtIP
      CASE('mugw')
        PropertyAtIP = mugwAtIP
      CASE('kgwpt')
        PropertyAtIP = KgwpTAtIP(TensorComponent(1),TensorComponent(2))
      CASE('kgwpp')
        PropertyAtIP = KgwppAtIP(TensorComponent(1),TensorComponent(2))
      CASE('CgwTT')
        PropertyAtIP = CGTTAtIP
      CASE('KGTT')
        PropertyAtIP = KGTTAtIP(TensorComponent(1),TensorComponent(2))
      CASE DEFAULT
        WRITE(Message,*) ' Variable "', TRIM(VariableName), '" not implemented.'
        CALL FATAL(SolverName,Message)
      END SELECT

      Weight = IP % s(t) * DetJ

      DO p=1,n
        DO q=1,n
          Stiff(p,q) = Stiff(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO
      FORCE(1:n) = FORCE(1:n) + Weight * PropertyAtIP * Basis(1:n)

    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE)
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixMaterialOutput
END SUBROUTINE PermafrostMaterialOutput


!==============================================================================
!> output of material parameter at IP
!==============================================================================
SUBROUTINE PermafrostIPOutput( Model,Solver,dt,TransientSimulation )
  !==============================================================================
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
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWfluxVar1,GWfluxVar2,GWfluxVar3,DepthVar
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: i,j,k,l,n,nb, nd,t, DIM, ok, NumberOfRockRecords, active,iter, maxiter, istat,DepthDOFs
  INTEGER,PARAMETER :: io=23
  INTEGER,POINTER :: TemperaturePerm(:), PressurePerm(:),&
       PorosityPerm(:),SalinityPerm(:),GWfluxPerm1(:),&
       TemperatureDtPerm(:), PressureDtPerm(:), SalinityDtPerm(:),&
       GWfluxPerm2(:),GWfluxPerm3(:)
  REAL(KIND=dp) :: Norm, meanfactor
  REAL(KIND=dp),POINTER :: Temperature(:), Pressure(:), Porosity(:), Salinity(:),&
       TemperatureDt(:), PressureDt(:), SalinityDt(:),&
       GWflux1(:),GWflux2(:),GWflux3(:)
  REAL(KIND=dp),POINTER :: NodalPorosity(:), NodalPressure(:), NodalSalinity(:),&
       NodalTemperature(:),NodalGWflux(:,:),NodalDepth(:),&
       NodalTemperatureDt(:), NodalSalinityDt(:),&
       NodalPressureDt(:)
  LOGICAL :: Found, FirstTime=.TRUE., AllocationsDone=.FALSE.,&
       ConstantPorosity=.TRUE., NoSalinity=.TRUE., NoPressure=.TRUE.,GivenGWFlux=.FALSE.,&
       ComputeDt=.FALSE.,ElementWiseRockMaterial, DepthExists
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: VariableBaseName(:)
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: SolverName='PermafrostIPOutput'
  CHARACTER(LEN=MAX_NAME_LEN) :: PressureName, PorosityName, SalinityName, GWfluxName, PhaseChangeModel,&
       ElementRockMaterialName,VarName, DepthName
  TYPE(ValueHandle_t) :: Load_h

  SAVE DIM,FirstTime,AllocationsDone,GivenGWFlux,DepthName,&
       CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,NumberOfRockRecords,&
       NodalPorosity,NodalPressure,NodalSalinity,NodalTemperature,NodalGWflux,NodalDepth,&
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt,&
       ElementWiseRockMaterial,ComputeDt

  CALL Info( SolverName, '-------------------------------------',Level=1 )
  CALL Info( SolverName, ' Assignment of IP variables          ',Level=1 )
  CALL Info( SolverName, '-------------------------------------',Level=1 )

  CALL AssignVars(Solver,Model,AllocationsDone,&
       NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux, &
       NodalTemperatureDt,NodalPressureDt,NodalSalinityDt, &
       TemperatureVar, PressureVar, PorosityVar,SalinityVar, &
       TemperatureDtVar, PressureDtVar, SalinityDtVar,&
       GWFluxVar1,GWFluxVar2,GWFluxVar3, &
       TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm, &       
       TemperatureDtPerm, PressureDtPerm, SalinityDtPerm, &
       GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
       Temperature, Pressure, Porosity,Salinity,&
       TemperatureDt, PressureDt, SalinityDt,&
       GWFlux1,GWFlux2,GWFlux3, &
       NoPressure, NoSalinity,ConstantPorosity,GivenGWFlux, DIM, ComputeDt, SolverName)

  Active = GetNOFActive()

  DO t=1,Active
    Element => GetActiveElement(t)
    Material => GetMaterial()
    IF (FirstTime) THEN

      ! check, whether we have globally or element-wise defined values of rock-material parameters
      ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
      IF (ElementWiseRockMaterial) THEN
        WRITE (Message,*) 'Found "Element Rock Material File"'
        CALL INFO(SolverName,Message,Level=3)
        CALL INFO(SolverName,'Using element-wise rock material definition',Level=3)
      END IF
      IF (ElementWiseRockMaterial) THEN
        ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
        NumberOfRockRecords = &
             ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,Solver,DIM)
      ELSE
        NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
      END IF

      IF (NumberOfRockRecords < 1) THEN
        CALL FATAL(SolverName,'No Rock Material specified')
      ELSE
        CALL INFO(SolverName,'Permafrost Rock Material read',Level=3)
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

    CALL ReadVars(N,Element,Model,Material,&
         NodalTemperature,NodalPressure,NodalPorosity,NodalSalinity,NodalGWflux,&
         Temperature, Pressure, Porosity,Salinity,GWFlux1,GWFlux2,GWFlux3,&
         TemperaturePerm, PressurePerm, PorosityPerm,SalinityPerm,&
         GWfluxPerm1, GWfluxPerm2,GWfluxPerm3, &
         NoSalinity,NoPressure,ConstantPorosity,GivenGWFlux,&
         PorosityName,SolverName,DIM)
    PRINT *,"SetIPValues"
    CALL SetIPValues(  Element, t, Active, n, nd+nb,&
         NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&           
         NodalGWflux, NodalDepth, GivenGWflux, DepthExists, &
         CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
         NumberOfRockRecords, PhaseChangeModel,ElementWiseRockMaterial)

  END DO

CONTAINS
  SUBROUTINE SetIPValues(   Element, ElementID, NoElements, n, nd,&
       NodalTemperature, NodalPressure, NodalPorosity, NodalSalinity,&
       NodalGWflux, NodalDepth,GivenGWflux,DepthExists, &
       CurrentRockMaterial, CurrentSoluteMaterial, CurrentSolventMaterial,&
       NumberOfRockRecords, PhaseChangeModel, ElementWiseRockMaterial)
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n, nd, ElementID, NoElements, NumberOfRockRecords
    TYPE(Element_t), POINTER :: Element
    TYPE(RockMaterial_t),POINTER :: CurrentRockMaterial
    TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
    TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
    REAL(KIND=dp) :: NodalTemperature(:), NodalSalinity(:),&
         NodalGWflux(:,:), NodalPorosity(:), NodalPressure(:),NodalDepth(:)
    LOGICAL, INTENT(IN) :: GivenGWflux, ElementWiseRockMaterial,DepthExists
    CHARACTER(LEN=MAX_NAME_LEN) :: PhaseChangeModel
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: RefDepth,CGTTAtIP, CgwTTAtIP, CGTpAtIP, CGTycAtIP,KGTTAtIP(3,3)   ! needed in equation
    REAL(KIND=dp) :: XiAtIP, Xi0Tilde,XiTAtIP,XiPAtIP,XiYcAtIP,XiEtaAtIP,&
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
         StiffPQ, meanfactor
    REAL(KIND=dp), POINTER :: gWork(:,:)
    INTEGER :: i,t,p,q,DIM, RockMaterialID
    LOGICAL :: Stat,Found, ConstantsRead=.FALSE.,ConstVal=.FALSE.,CryogenicSuction=.FALSE.,HydroGeo=.FALSE.
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    CHARACTER(LEN=MAX_NAME_LEN) :: MaterialFileName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='Permafrost(LocalMatrixHTEQ)'
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

    ! Get stuff from SIF Material section
    Material => GetMaterial(Element)
    IF (ElementWiseRockMaterial) THEN
      RockMaterialID = ElementID  ! each element has it's own set of parameters
    ELSE
      RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
    END IF

    IF (DepthExists) THEN
      RefDepth = GetConstReal(Material,'Radiogenic Reference Depth',Found)
      IF (Found) THEN
        DO I=1,N

          !               RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth)
          !PRINT *,"HTEQ: RGEN",RadiogenicHeatProduction(CurrentRockMaterial,RockMaterialID,NodalDepth(I),RefDepth), NodalDepth(I)
        END DO

      END IF
    END IF

    HydroGeo = GetLogical(Material,'Hydrogeological Model',Found)
    IF (.NOT.Found) HydroGeo = .FALSE.

    ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
    IF (.NOT.Found) THEN
      ConstVal = .FALSE.
    ELSE
      IF (ConstVal) &
           CALL INFO(FunctionName,'"Constant Permafrost Properties" set to true',Level=9)
    END IF

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

    ! Loop all Gauss-points
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      ! The heat soruce term
      !LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found, GaussPoint=t)
      !IF (LoadAtIP > 0.0_dp) PRINT *,"HTEQ:LoadAtIP", LoadAtIP
      ! Contribution from other heat source
      !LoadAtIP = LoadAtIP + SUM( Basis(1:n) * LOAD(1:n) )

      ! Variables (Temperature, Porosity, Pressure, Salinity) at IP
      TemperatureAtIP = SUM( Basis(1:N) * NodalTemperature(1:N) )
      PorosityAtIP = SUM( Basis(1:N) * NodalPorosity(1:N))
      PressureAtIP = SUM( Basis(1:N) * NodalPressure(1:N))      
      SalinityAtIP = SUM( Basis(1:N) * NodalSalinity(1:N))

      
      !Materialproperties needed for computing Xi at IP

      rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
      rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!      
      Xi0Tilde = GetXi0Tilde(CurrentRockMaterial,RockMaterialID,PorosityAtIP)

      ! unfrozen pore-water content at IP
      SELECT CASE(PhaseChangeModel)
      CASE('anderson')
        XiAtIP = &
             GetXiAnderson(0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiTAtIP = &
             XiAndersonT(XiAtIP,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)
        XiPAtIP   = &
             XiAndersonP(XiAtIp,0.011_dp,-0.66_dp,9.8d-08,&
             CurrentSolventMaterial % rhow0,CurrentRockMaterial % rhos0(RockMaterialID),&
             T0,TemperatureAtIP,PressureAtIP,PorosityAtIP)       
      CASE DEFAULT ! Hartikainen model
        CALL  GetXiHartikainen (CurrentRockMaterial,RockMaterialID,&
             CurrentSoluteMaterial,CurrentSolventMaterial,&
             TemperatureAtIP,PressureAtIP,SalinityAtIP,PorosityAtIP,&
             Xi0tilde,deltaInElement,rhowAtIP,rhoiAtIP,&
             GasConstant,p0,T0,&
             XiAtIP,XiTAtIP,XiYcAtIP,XiPAtIP,XiEtaAtIP,&
             .FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.)
      END SELECT

      !Materialproperties needed at IP:
      rhowAtIP = rhowupdate(CurrentSolventMaterial,rhowAtIP,XiAtIP,SalinityAtIP,ConstVal)
      rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)!!
      rhocAtIP = rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"SetIPValues: rhowAtIP, rhoiAtIP, rhosAtIP", rhowAtIP, rhoiAtIP, rhosAtIP

      ! heat capacities
      csAtIP   = cs(CurrentRockMaterial,RockMaterialID,&
           T0,TemperatureAtIP,ConstVal)
      cwAtIP   = cw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"SetIPValues:cw",T0,TemperatureAtIP,SalinityAtIP,cw0,&
      !     acw,bcw,acwl,bcwl
      !PRINT *, "SetIPValues:cwAtIP", cwAtIP, "cw0",cw0,"acw",acw,"bcw",bcw,"T0",T0,SalinityAtIP,TemperatureAtIP,PressureAtIP
      ciAtIP   = ci(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      !ci(ci0,aci,T0,TemperatureAtIP,PressureAtIP)  !!
      ccAtIP   = cc(CurrentSoluteMaterial,&
           T0,TemperatureAtIP,SalinityAtIP,ConstVal)
      !PRINT *,"SetIPValues: cw,ci,cs,cc",cwAtIP,ciAtIP,csAtIP,ccAtIP

      ! latent heat
      hiAtIP = hi(CurrentSolventMaterial,&
           T0,TemperatureAtIP,ConstVal)
      hwAtIP = hw(CurrentSolventMaterial,&
           T0,XiAtIP,TemperatureAtIP,SalinityAtIP,ConstVal)
 

      ! heat conductivity at IP
      ksthAtIP = GetKalphath(CurrentRockMaterial % ks0th(RockMaterialID),&
           CurrentRockMaterial % bs(RockMaterialID),T0,TemperatureAtIP)
      kwthAtIP = GetKalphath(CurrentSolventMaterial % kw0th,CurrentSolventMaterial % bw,T0,TemperatureAtIP)
      kithAtIP = GetKalphath(CurrentSolventMaterial % ki0th,CurrentSolventMaterial % bi,T0,TemperatureAtIP)
      kcthAtIP = GetKalphath(CurrentSoluteMaterial % kc0th,CurrentSoluteMaterial % bc,T0,TemperatureAtIP)
      KGTTAtIP = GetKGTT(ksthAtIP,kwthAtIP,kithAtIP,kcthAtIP,XiAtIP,&
           SalinityATIP,PorosityAtIP,meanfactor)
 

      ! heat capacities at IP
      CGTTAtIP = &
           GetCGTT(XiAtIP,XiTAtIP,rhosAtIP,rhowAtIP,rhoiAtIP,rhocAtIP,&
           cwAtIP,ciAtIP,csAtIP,ccAtIP,hiAtIP,hwAtIP,&
           PorosityAtIP,SalinityAtIP)
 
      CgwTTAtIP = GetCgwTT(rhowAtIP,rhocAtIP,cwAtIP,ccAtIP,XiAtIP,SalinityAtIP)
 
      ! compute groundwater flux for advection term

      CGTpAtIP = GetCGTp(rhoiAtIP,hiAtIP,hwAtIP,XiPAtIP,PorosityAtIP) !NEW
      CGTycAtIP = GetCGTyc(rhoiAtIP,hiAtIP,hwAtIP,XiYcAtIP,PorosityAtIP) !NEW
      ! groundwater flux
      !PRINT *, "SetIPValues:KGTTAtIP", KGTTAtIP,"CgwTTAtIP",CgwTTAtIP
      JgwDAtIP = 0.0_dp
      IF (GivenGWFlux) THEN
        !PRINT *, "SetIPValues: Interpolate Flux"
        DO I=1,DIM
          JgwDAtIP(I) = SUM( Basis(1:N) * NodalGWflux(I,1:N)) 
        END DO
      ELSE        
        !PRINT *, "SetIPValues: Compute Flux"
        mugwAtIP = mugw(CurrentSolventMaterial,CurrentSoluteMaterial,&
             XiAtIP,T0,SalinityAtIP,TemperatureAtIP,ConstVal)
        KgwAtIP = 0.0_dp
        KgwAtIP = GetKgw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             mugwAtIP,XiAtIP,MinKgw)
        fwAtIP = fw(CurrentRockMaterial,RockMaterialID,CurrentSolventMaterial,&
             Xi0tilde,rhowAtIP,XiAtIP,GasConstant,TemperatureAtIP)
        KgwpTAtIP = GetKgwpT(fwAtIP,XiTAtIP,KgwAtIP)
        IF (CryogenicSuction) THEN
          KgwppAtIP = GetKgwpp(fwAtIP,XiPAtIP,KgwAtIP)
        ELSE
          KgwppAtIP = KgwAtIP
        END IF
        !PRINT *,"SetIPValues: KgwppAtIP",KgwppAtIP
        rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)

        DO i=1,DIM
          gradTAtIP(i) =  SUM(NodalTemperature(1:N)*dBasisdx(1:N,i))
          gradPAtIP(i) =  SUM(NodalPressure(1:N) * dBasisdx(1:N,i))
        END DO

        JgwDAtIP = GetJgwD(KgwppAtIP,KgwpTAtIP,KgwAtIP,gradpAtIP,gradTAtIP,Gravity,rhogwAtIP,DIM,CryogenicSuction)
        !PRINT *,"SetIPValues: JgwD=(",JgwDAtIP(1:DIM)*365.5*24.0*3600.0,")"
      END IF

      ! add thermal dispersion in Hydro-Geological Mode
      IF (HydroGeo) THEN
        DtdAtIP = GetDtd(CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP,JgwDAtIP)
        DO I=1,DIM
          DO J=1,DIM
            KGTTAtIP(I,J) = KGTTAtIP(I,J) + CGWTTAtIP * DtdAtIP(I,J)
          END DO
        END DO
      END IF


    END DO
  END SUBROUTINE SetIPValues
END SUBROUTINE PermafrostIPOutput
!---------------------------------------------------------------------------------------------
! Functions needed for permafrost model (might be shifted to USF-directory)
!---------------------------------------------------------------------------------------------
FUNCTION GetKGuu(Model,IPNo,PorosityAtIP) RESULT(KGuuAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: PorosityAtIP, KGuuAtIP(6,6)
  !-----------
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: EGAtIP, nuGAtIP
  TYPE(Variable_t), POINTER :: XiAtIPVar
  INTEGER, POINTER :: XiAtIPPerm(:)
  REAL(KIND=dp), POINTER :: XiAtIP(:)
  LOGICAL :: FirstTime = .TRUE., ElementWiseRockMaterial, Found
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetKGuu)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentSolventMaterial,CurrentRockMaterial,DIM,ElementWiseRockMaterial


  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  XiAtIPVar => VariableGet( Model % Mesh % Variables, 'Xi')
  IF (.NOT.ASSOCIATED(XiAtIPVar)) THEN
    WRITE(Message,*) 'Variable Xi is not associated'
    CALL FATAL(FunctionName,Message)
  END IF
  XiAtIPPerm => XiAtIPVar % Perm
  XiAtIp => XiAtIPVar % Values
  IPPerm = XiAtIPPerm(t) + IPNo

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=3)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=3)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=3)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found, UnfoundFatal=.TRUE.)
  END IF

  EGAtIP = EG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
  nuGAtIP = nuG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP(IPPerm),PorosityAtIP)
  KGuuAtIP = KGuu(EGAtIP,nuGAtIP,DIM)
END FUNCTION GetKGuu
!---------------------------------------------------------------------------------------------
FUNCTION GetBetaG(Model,IPNo,ArgumentsAtIP) RESULT(betaGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), betaGAtIP
  !--------------
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: EGAtIP, nuGAtIP, PorosityAtIP, XiAtIP
  LOGICAL :: FirstTime = .TRUE., ElementWiseRockMaterial, Found
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetBetaG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentRockMaterial,DIM,ElementWiseRockMaterial

  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=3)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=3)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=3)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found, UnfoundFatal=.TRUE.)
  END IF
  PorosityAtIP = ArgumentsAtIP(1)
  XiAtIP = ArgumentsAtIP(2)
  betaGAtIP = betaG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
END FUNCTION GetBetaG
  !---------------------------------------------------------------------------------------------
FUNCTION GetNuG(Model,IPNo,ArgumentsAtIP) RESULT(nuGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), nuGAtIP
  
  !-----
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: PorosityAtIP, XiAtIP
  LOGICAL :: Found, FirstTime = .TRUE., ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetNuG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentRockMaterial,CurrentSolventMaterial,DIM,ElementWiseRockMaterial
  
  IF (FirstTime) CALL INFO("Permafrost(GetNuG)","Initializing",Level=1)
  PorosityAtIP=ArgumentsAtIP(1)
  XiAtIP=ArgumentsAtIP(2)
  !XiAtIP=1.0_dp 
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  t = Element % ElementIndex
  Material => GetMaterial(Element)

  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=3)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=3)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=3)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
  END IF

  nuGAtIP = nuG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
  !PRINT *,"getNuG:", nuGAtIP, XiAtIp(IPPerm),PorosityAtIP
END FUNCTION GetNuG
!---------------------------------------------------------------------------------------------
FUNCTION GetEG(Model,DummyIPNo,ArgumentsAtIP) RESULT(EGAtIP)
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: DummyIPNo
  REAL(KIND=dp) :: ArgumentsAtIP(2), EGAtIP
  !-----
  TYPE(Solver_t) :: DummySolver
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t),POINTER :: Element
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  INTEGER :: RockMaterialID, NumberOfRockRecords, DIM, t, i, IPPerm
  REAL(KIND=dp) :: PorosityAtIP, XiAtIP
  LOGICAL :: Found,FirstTime = .TRUE., ElementWiseRockMaterial
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'PermafrostMaterials (GetNuG)'
  !-----------
  SAVE FirstTime,NumberOfRockRecords,CurrentRockMaterial,CurrentSolventMaterial,DIM,ElementWiseRockMaterial

  PorosityAtIP=ArgumentsAtIP(1)
  XiAtIP=ArgumentsAtIP(2)
  !PRINT *, "GetEG:", PorosityAtIP, XiAtIP
  
  Element => Model % CurrentElement
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  !t = Element % ElementIndex
  Material => GetMaterial(Element)
  
  IF (FirstTime .OR. (Model % Mesh % Changed)) THEN
    DIM =  CoordinateSystemDimension()

    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=3)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=3)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
    END IF

    IF (NumberOfRockRecords < 1) THEN
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=3)
      FirstTime = .FALSE.
    END IF
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t  ! each element has it's own set of parameters
  ELSE
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found,UnfoundFatal=.TRUE.)
  END IF
  EGAtIP = EG(CurrentSolventMaterial,CurrentRockMaterial,RockMaterialID,XiAtIP,PorosityAtIP)
  !IF ((EGAtIP < 1.0d05) .OR. (EGAtIP > 1.0d07)) PRINT *,"GetEG",EGAtIP,XiAtIP,PorosityAtIP
END FUNCTION GetEG
!---------------------------------------------------------------------------------------------
FUNCTION GetElasticityForce(Model,IPNo,ArgumentsAtIP) RESULT(EforceAtIP) ! needs arguments Temperature, Pressure, Porosity, Salinity
  USE DefUtils
  USE PermaFrostMaterials
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER, INTENT(IN) :: IPNo
  REAL(KIND=dp) :: ArgumentsAtIP(5), EforceAtIP
  !--------------
  REAL(KIND=dp) :: TemperatureAtIP, PressureAtIP, PorosityAtIP, SalinityAtIP,XiAtIP,&
       rhogwAtIP, rhosAtIP, rhowAtIP,rhocAtIP, rhoiAtIP,rhoGAtIP,&
       GasConstant, N0, DeltaT, T0, p0, eps, Gravity(3),InitialOffsetRhoAtIP
  TYPE(Variable_t), POINTER :: RhoOffsetAtIPVar
  INTEGER, POINTER :: RhoOffsetAtIPPerm(:)
  REAL(KIND=dp), POINTER :: RhoOffsetAtIP(:)
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  INTEGER ::  DIM, t,NumberOfRockRecords, RockMaterialID, IPPerm
  TYPE(SolventMaterial_t), POINTER :: CurrentSolventMaterial
  TYPE(SoluteMaterial_t), POINTER :: CurrentSoluteMaterial
  TYPE(RockMaterial_t), POINTER :: CurrentRockMaterial
  TYPE(Solver_t) :: DummySolver
  CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName = 'Permafrost (GetElasticityForce)'
  CHARACTER(LEN=MAX_NAME_LEN) :: ElementRockMaterialName
  LOGICAL :: Found,FirstTime=.TRUE.,ConstVal=.FALSE., ConstantsRead = .FALSE.,&
       ElementWiseRockMaterial, OffsetDensity=.FALSE.

  SAVE ConstantsRead,ElementWiseRockMaterial,GasConstant, DIM, N0, DeltaT, T0, p0, eps, Gravity,&
       NumberOfRockRecords,FirstTime,CurrentRockMaterial,CurrentSoluteMaterial,CurrentSolventMaterial,&
       OffsetDensity

  
  IF (.NOT.ConstantsRead) THEN
    ConstantsRead = &
         ReadPermafrostConstants(Model, FunctionName, DIM, GasConstant, N0, DeltaT, T0, p0, eps, Gravity)
    OffsetDensity = GetLogical(Model %Constants,'Permafrost Offset Density', Found)
    IF (.NOT.Found) THEN
      OffsetDensity = .FALSE.
      CALL WARN(FunctionName,'No offset for groundwater pressure included - might lead to artifial high compression')
    ELSE
      CALL INFO(FunctionName,'Offset groundwater pressure is activated',Level=1)
    END IF
  END IF
  Element => Model % CurrentElement
  t = Element % ElementIndex
  IF (.NOT.ASSOCIATED(Element)) CALL FATAL(FunctionName,'Element not associated')
  IF (OffsetDensity) THEN
    RhoOffsetAtIPVar => VariableGet( Model % Mesh % Variables, 'Reference Offset Density')
    IF (.NOT.ASSOCIATED(RhoOffsetAtIPVar)) THEN
      WRITE(Message,*) '"Permafrost Offset Density" is set, but variable "Reference Offset Density" is not associated'
      CALL FATAL(FunctionName,Message)
    END IF
    RhoOffsetAtIPPerm => RhoOffsetAtIPVar % Perm
    IPPerm = RhoOffsetAtIPPerm(t) + IPNo  
    RhoOffsetAtIP => RhoOffsetAtIPVar % Values
    InitialOffsetRhoAtIP = RhoOffsetAtIP(IPPerm)
  ELSE
    InitialOffsetRhoAtIP = 0.0_dp
  END IF

  Material => GetMaterial(Element)
  ConstVal = GetLogical(Material,'Constant Permafrost Properties',Found)
  IF (FirstTime) THEN
    ! check, whether we have globally or element-wise defined values of rock-material parameters
    ElementRockMaterialName = GetString(Material,'Element Rock Material File',ElementWiseRockMaterial)
    IF (ElementWiseRockMaterial) THEN
      WRITE (Message,*) 'Found "Element Rock Material File"'
      CALL INFO(FunctionName,Message,Level=3)
      CALL INFO(FunctionName,'Using element-wise rock material definition',Level=3)
    END IF
    IF (ElementWiseRockMaterial) THEN
      ! read element-wise material parameter (CurrentRockMaterial will have one entry each element)
      NumberOfRockRecords = &
           ReadPermafrostElementRockMaterial(CurrentRockMaterial,ElementRockMaterialName,DummySolver,DIM,SkipInit=.TRUE.)
    ELSE
      NumberOfRockRecords =  ReadPermafrostRockMaterial( Material,Model % Constants,CurrentRockMaterial )
    END IF
    IF (NumberOfRockRecords < 1) THEN
      PRINT *, "NumberOfRockRecords=", NumberOfRockRecords
      CALL FATAL(FunctionName,'No Rock Material specified')
    ELSE
      CALL INFO(FunctionName,'Permafrost Rock Material read',Level=3)
      FirstTime = .FALSE.
    END IF
    CALL ReadPermafrostSoluteMaterial( Material,Model % Constants,CurrentSoluteMaterial )
    CALL SetPermafrostSolventMaterial( CurrentSolventMaterial )
  END IF

  IF (ElementWiseRockMaterial) THEN
    RockMaterialID = t
  ELSE      
    RockMaterialID = ListGetInteger(Material,'Rock Material ID', Found)
    IF (.NOT.Found) CALL FATAL(FunctionName,"Rock Material ID not found")
  END IF

  TemperatureAtIP = ArgumentsAtIP(1)
  PressureAtIP    = ArgumentsAtIP(2)
  PorosityAtIP    = ArgumentsAtIP(3)
  SalinityAtIP    = ArgumentsAtIP(4)
  XiAtIP          = ArgumentsAtIP(5)
  
  rhocAtIP =  rhoc(CurrentSoluteMaterial,T0,p0,XiAtIP,TemperatureAtIP,PressureAtIP,SalinityAtIP,ConstVal)
  !IF (rhocAtIP .NE. rhocAtIP) CALL FATAL(FunctionName,'rhocAtIP is NaN')  
  rhowAtIP = rhow(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  !IF (rhowAtIP .NE. rhocAtIP) CALL FATAL(FunctionName,'rhowAtIP is NaN') 
  rhoiAtIP = rhoi(CurrentSolventMaterial,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  rhosAtIP = rhos(CurrentRockMaterial,RockMaterialID,T0,p0,TemperatureAtIP,PressureAtIP,ConstVal)
  rhogwAtIP = rhogw(rhowAtIP,rhocAtIP,XiAtIP,SalinityAtIP)
  !!!! IF (Time == 0) store that one at IP

  
  rhoGAtIP = rhoG(rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP)
  IF (rhoGAtIP .NE. rhoGAtIP) THEN
    PRINT *,rhosAtIP,rhogwAtIP,rhoiAtIP,PorosityAtIP,SalinityAtIP,XiAtIP
    CALL FATAL(FunctionName,'rhoGAtIP is NaN')
  END IF
  
  !EforceAtIP = -rhoGAtIP * SQRT(SUM(Gravity(1:3)*Gravity(1:3)))
  EforceAtIP = (rhoGAtIP - InitialOffsetRhoAtIP)* Gravity(DIM)
  !IF (t==100 )PRINT *, "GetElasticityForce: EforceAtIP=", EforceAtIP, "=(",rhoGAtIP, "-", InitialOffsetRhoAtIP,")*", Gravity(DIM)
  !     PorosityAtIP,SalinityAtIP,XiAtIP, "   *********"  
END FUNCTION GetElasticityForce







