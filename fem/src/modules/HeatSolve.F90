!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Module containing a solver for heat equation
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE HeatSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'HeatSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  
  Params => GetSolverParams()

  ! Default variable name
  CALL ListAddNewString( Params,'Variable','Temperature')

  ! These use one flag to call library features to compute automatically
  ! a conductivity matrix.
  IF( ListGetLogical(Params,'Calculate Conductance Matrix',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Matrix Symmetric',.TRUE.)
    CALL ListAddNewString( Params,'Constraint Modes Matrix Filename',&
        'ThermalConductanceMatrix.dat',.FALSE.)
    CALL ListRenameAllBC( Model,'Conductivity Body','Constraint Mode Temperature')
  END IF

  ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVIVTY
  CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
  
END SUBROUTINE HeatSolver_Init


!------------------------------------------------------------------------------
!> Subroutine for solving the energy a.k.a. heat equation in various coordinate systems.
!> \ingroup Solvers
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE HeatSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DiffuseConvective
     USE DiffuseConvectiveGeneral
     USE Differentials
     USE Radiation
     USE MaterialModels
     USE Adaptive
     USE DefUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     INTEGER, PARAMETER :: PHASE_SPATIAL_1 = 1
     INTEGER, PARAMETER :: PHASE_SPATIAL_2 = 2
     INTEGER, PARAMETER :: PHASE_TEMPORAL  = 3
    
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver

     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix

     INTEGER :: i,j,k,l,m,n,nd,t,tt,iter,k1,k2,body_id,eq_id,istat,LocalNodes,bf_id

     TYPE(Nodes_t)   :: ElementNodes
     TYPE(Element_t),POINTER :: Element,Parent,RadiationElement

     REAL(KIND=dp) :: RelativeChange, &
           Norm,PrevNorm,Text,S,C,C1,Emissivity,StefanBoltzmann, &
           ReferencePressure=0.0d0, SpecificHeatRatio

     CHARACTER(LEN=MAX_NAME_LEN) :: RadiationFlag,ConvectionFlag

     INTEGER :: PhaseChangeModel
     CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel, StabilizeFlag, VarName

     INTEGER, POINTER :: NodeIndexes(:)
     LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles,NewtonLinearization = .FALSE., &
         Found, GotIt, HeatFluxBC, HeatGapBC, GotMeltPoint, IsRadiation, IsRadiosity, InfBC
! Which compressibility model is used
     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, ConvectionField
     INTEGER :: CompressibilityModel

     LOGICAL :: AllocationsDone = .FALSE.,PhaseSpatial=.FALSE., &
        PhaseChange=.FALSE., CheckLatentHeatRelease=.FALSE., FirstTime, &
        SmartHeaterControl, IntegralHeaterControl, HeaterControlLocal, SmartTolReached=.FALSE., &
        TransientHeaterControl, SmartHeaterAverage, ConstantBulk, SaveBulk, &
	TransientAssembly, Converged, AnyMultiply, NeedFlowSol
     LOGICAL, POINTER :: SmartHeaters(:), IntegralHeaters(:)

     TYPE(Variable_t), POINTER :: TempSol,FlowSol,HeatSol,CurrentSol, MeshSol, DensitySol
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER, POINTER :: TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:)

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter,MDOFs, &
         SmartHeaterBC, SmartHeaterNode, DoneTime=0, bc_elem, nb, NOFactive
     REAL(KIND=dp) :: NonlinearTol,NewtonTol,SmartTol,Relax, &
            SaveRelax,dt,dt0,CumulativeTime, VisibleFraction, PowerScaling=1.0, PrevPowerScaling=1.0, &
            PowerRelax, PowerTimeScale, PowerSensitivity, xave, yave, Normal(3), &
	    dist, mindist, ControlPoint(3), HeatTransferMultiplier

     REAL(KIND=dp), POINTER :: Temperature(:),PrevTemperature(:),FlowSolution(:), &
       ElectricCurrent(:), PhaseChangeIntervals(:,:),ForceVector(:), &
       PrevSolution(:), HC(:), Hwrk(:,:,:),MeshVelocity(:), XX(:), YY(:),ForceHeater(:),&
       RealWork(:,:)

     REAL(KIND=dp), ALLOCATABLE :: vals(:)
     REAL(KIND=dp) :: Jx,Jy,Jz,JAbs, Power, MeltPoint, IntHeatSource

     INTEGER, ALLOCATABLE, SAVE :: Indexes(:), SaveIndexes(:)

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
       FORCE(:), U(:), V(:), W(:), MU(:,:),TimeForce(:), &
       Density(:), LatentHeat(:), HeatTransferCoeff(:), &
       HeatCapacity(:), Enthalpy(:), EnthalpyFraction(:), Viscosity(:), LocalTemperature(:), &
       NodalVal(:), ElectricConductivity(:), Permeability(:), Work(:), C0(:), &
       Pressure(:), dPressuredt(:), GasConstant(:),AText(:), HeaterArea(:), &
       HeaterTarget(:), HeaterScaling(:), HeaterDensity(:), HeaterSource(:), &
       HeatExpansionCoeff(:), ReferenceTemperature(:), PressureCoeff(:), &
       PhaseVelocity(:,:), HeatConductivityIso(:), &
       PerfusionRate(:), PerfusionDensity(:), PerfusionHeatCapacity(:), PerfusionRefTemperature(:)

     REAL(KIND=dp), ALLOCATABLE :: Areas(:), Emiss(:), Reflect(:)
     LOGICAL :: Spectral, Radiosity
     
     SAVE U, V, W, MU, MASS, STIFF, LOAD, PressureCoeff, &
       FORCE, ElementNodes, HeatConductivity, HeatCapacity, HeatTransferCoeff, &
       Enthalpy, EnthalpyFraction, Density, LatentHeat, PhaseVelocity, AllocationsDone, Viscosity, TimeForce, &
       LocalNodes, LocalTemperature, Work, ElectricConductivity, &
       NodalVal, Permeability, C0, dPressuredt, Pressure, &
       GasConstant,AText,Hwrk, XX, YY, ForceHeater, Power, HeaterArea, HeaterTarget, &
       HeaterScaling, HeaterDensity, HeaterSource, SmartHeaters, IntegralHeaters, SmartTolReached,    &
       ReferenceTemperature, HeatExpansionCoeff, PrevPowerScaling, PowerScaling, &
       MeltPoint, DoneTime, SmartHeaterNode, SmartHeaterBC, SmartHeaterAverage, &
       HeatConductivityIso, &
       PerfusionRate, PerfusionDensity, PerfusionHeatCapacity, PerfusionRefTemperature


     INTERFACE
        FUNCTION HeatSolver_Boundary_Residual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatSolver_Boundary_Residual

        FUNCTION HeatSolver_Edge_Residual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION HeatSolver_Edge_Residual

        FUNCTION HeatSolver_Inside_Residual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatSolver_Inside_Residual
     END INTERFACE

     REAL(KIND=dp) :: at,at0,totat,st,totst,t1


     CALL Info('HeatSolver','-------------------------------------------',Level=6)
     CALL Info('HeatSolver','Solving the energy equation for temperature',Level=5)

     IF( ListCheckPresentAnyBC( Model,'Heat Gap') ) THEN
       CALL Warn('HeatSolver','The old way of dealing with HeatGap is obsolite!!')
     END IF

     
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     SolverParams => GetSolverParams()

     Radiosity = GetLogical( SolverParams,'Radiosity Model',Found ) 
     Spectral = GetLogical( SolverParams,'Spectral Model',Found )
     IF( Spectral ) Radiosity = .TRUE.

     IF(.NOT. Radiosity ) THEN
       IF( ListCheckPresent( SolverParams,'Radiator Coordinates' ) .OR. &
           ListCheckPresentAnyBodyForce(Model,'Radiator Coordinates') ) THEN
         CALL Fatal('HeatSolve','For radiative point sources use HeatSolveVec instead!')
       END IF
     END IF

!------------------------------------------------------------------------------
!    The View and Gebhart factors may change. If this is necessary, this is 
!    done within this subroutine. The routine is called in the
!    start as it may affect the matrix topology.
!    Newton lineariarization option is needed only when there is radiation.
!------------------------------------------------------------------------------
     IsRadiation = ListCheckPresentAnyBC( Model,'Radiation')
     
     IF( IsRadiation .AND. .NOT. Radiosity ) THEN
       CALL RadiationFactors( Solver, .FALSE., .FALSE.)       
     END IF

     ! The solver matrix, permutation etc. may change because radiation factors are recomputed
     StiffMatrix => GetMatrix()
     ForceVector => Solver % Matrix % RHS
     TempSol => Solver % Variable
     TempPerm    => TempSol % Perm
     Temperature => TempSol % Values
     VarName = GetVarName( TempSol ) 
  
     LocalNodes = COUNT( TempPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN
     IF(SIZE(Temperature) < LocalNodes) LocalNodes = SIZE(Temperature)
     
     
     NeedFlowSol = .FALSE.
     DO i=1,Model % NumberOfEquations
       ConvectionFlag = GetString( Model % Equations(i) % Values, 'Convection', Found )
       IF ( ConvectionFlag == 'computed' ) THEN
         NeedFlowSol = .TRUE.
         EXIT
       END IF
     END DO

     FlowSol => NULL()
     IF( NeedFlowSol ) THEN
       ConvectionField = GetString( SolverParams, 'Temperature Convection Field', Found )     
       IF ( Found ) THEN       
         FlowSol => VariableGet( Solver % Mesh % Variables, ConvectionField )
       ELSE
         FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
       END IF
      
       IF ( ASSOCIATED( FlowSol ) ) THEN
         FlowPerm     => FlowSol % Perm
         NSDOFs       =  FlowSol % DOFs
         FlowSolution => FlowSol % Values
       ELSE
         CALL Fatal('HeatSolver','Flow is "computed" but not flow field available!')
       END IF
     END IF
       
     DensitySol => VariableGet( Solver % Mesh % Variables, 'Density' )

     ! Check whether we have some heater controls. This will affect initialization stuff. 
     SmartHeaterControl = ListCheckPresentAnyBodyForce( Model,'Smart Heater Control')
     IntegralHeaterControl = ListCheckPresentAnyBodyForce( Model,'Integral Heat Source')
   
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN
        N = Solver % Mesh % MaxElementDOFs

        IF ( AllocationsDone ) THEN
          DEALLOCATE(  &
                 U, V, W, MU,           &
                 Pressure,              &
                 dPressureDt,           &
                 PressureCoeff,        &
                 ElementNodes % x,      &
                 ElementNodes % y,      &
                 ElementNodes % z,      &
                 Density,Work,          &
                 LatentHeat,            &
                 PhaseVelocity,         &
                 ElectricConductivity,  &
                 Permeability,          &
                 Viscosity,C0,          &
                 HeatTransferCoeff,     &
                 HeatExpansionCoeff,    &
                 ReferenceTemperature,  &
                 MASS,       &
                 LocalTemperature,      &
                 HeatCapacity,Enthalpy, &
                 EnthalpyFraction,      &
                 NodalVal,       &
                 GasConstant, AText,    &
                 HeatConductivity,      &
                 STIFF,LOAD,            &
                 Indexes, SaveIndexes,  &
                 FORCE, TimeForce,      &
                 HeatConductivityIso,   &
                 PerfusionRate,         &
                 PerfusionDensity,      &
                 PerfusionHeatCapacity, &
                 PerfusionRefTemperature )
       END IF

       ALLOCATE( &
                 Indexes(N), SaveIndexes(N),           &
                 U( N ),   V( N ),  W( N ),            &
                 MU( 3,N ),                            &
                 Pressure( N ),                        &
                 dPressureDt( N ),                     &
                 PressureCoeff( N ),                   &
                 ElementNodes % x( N ),                &
                 ElementNodes % y( N ),                &
                 ElementNodes % z( N ),                &
                 Density( N ),Work( N ),               &
                 LatentHeat( N ),                      &
                 PhaseVelocity(3, N),                  &
                 ElectricConductivity(N),              &
                 Permeability(N),                      &
                 Viscosity(N),C0(N),                   &
                 HeatTransferCoeff( N ),               &
                 HeatExpansionCoeff( N ),              &
                 ReferenceTemperature( N ),            &
                 MASS(  2*N,2*N ),                     &
                 LocalTemperature( N ),                &
                 HeatCapacity( N ),Enthalpy( N ),      &
                 EnthalpyFraction( N ),                &
                 NodalVal( N ),                 &
                 GasConstant( N ),AText( N ),          &
                 HeatConductivity( 3,3,N ),            &
                 HeatConductivityIso( N ),             &
                 STIFF( 2*N,2*N ),LOAD( N ), &
                 FORCE( 2*N ), TimeForce(2*N), &
                 PerfusionRate( N ),     &
                 PerfusionDensity( N ),      &
                 PerfusionHeatCapacity( N ), &
                 PerfusionRefTemperature( N ), &
                 STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'HeatSolve', 'Memory allocation error' )
       END IF


       IF ( SmartHeaterControl .OR. IntegralHeaterControl) THEN          
          n = Model % NumberOfBodyForces
          IF ( AllocationsDone ) DEALLOCATE( HeaterArea, HeaterDensity, HeaterSource, &
               HeaterScaling, HeaterTarget, SmartHeaters, IntegralHeaters )
          ALLOCATE( HeaterArea(n), HeaterDensity(n), HeaterSource(n), &
               HeaterScaling(n), HeaterTarget(n), SmartHeaters(n), &
               IntegralHeaters(n) )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          SmartHeaters = .FALSE.
          IntegralHeaters = .FALSE.
       END IF
       
       IF( SmartHeaterControl ) THEN
          IF ( AllocationsDone ) DEALLOCATE( XX, YY, ForceHeater  )
          n = SIZE( Temperature )
          ALLOCATE( XX( n ), YY(n), ForceHeater( n ), STAT=istat )
          IF ( istat /= 0 ) THEN
             CALL Fatal( 'HeatSolve', 'Memory allocation error' )
          END IF
          XX = 0.0d0 
          YY = 0.0d0
          ForceHeater = 0.0d0
       END IF
       
       NULLIFY( Hwrk )
       AllocationsDone = .TRUE.
    END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
     dt = Timestep
     Constants => GetConstants()
     IF( IsRadiation ) THEN
       StefanBoltzmann = ListGetConstReal( Model % Constants, &
                     'Stefan Boltzmann',UnfoundFatal=.TRUE.)
     END IF

!------------------------------------------------------------------------------
     Stabilize = GetLogical( SolverParams,'Stabilize',Found )

     UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
     IF ( .NOT.Found ) UseBubbles = .TRUE.

     StabilizeFlag = GetString( SolverParams, &
          'Stabilization Method',Found )

     SELECT CASE(StabilizeFlag)
     CASE('vms')
       Stabilize = .FALSE.
       UseBubbles= .FALSE.
     CASE('stabilized')
       Stabilize = .TRUE.
       UseBubbles = .FALSE.
     CASE('bubbles')
       Stabilize = .FALSE.
       UseBubbles = .TRUE.
     END SELECT

     NonlinearIter = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 1

     NonlinearTol  = GetConstReal( SolverParams, &
                     'Nonlinear System Convergence Tolerance',    Found )

     IF( IsRadiation ) THEN
       NewtonTol     = GetConstReal( SolverParams, &
                      'Nonlinear System Newton After Tolerance',  Found )
       NewtonIter    = GetInteger(   SolverParams, &
                      'Nonlinear System Newton After Iterations', Found )
     ELSE
       NewtonTol = 1.0_dp
       NewtonIter =  0
     END IF
     IF ( NewtonIter == 0) NewtonLinearization = .TRUE.

     Relax = GetCReal( SolverParams,'Nonlinear System Relaxation Factor',Found )
     IF ( .NOT.Found ) Relax = 1

     TransientAssembly = TransientSimulation
     dt0 = ListGetConstReal(SolverParams,'Steady State Transition Timestep',Found)
     IF(.NOT. Found) dt0 = ListGetConstReal(SolverParams,'Smart Heater Time Scale',Found)

     IF(Found .AND. dt > dt0) TransientAssembly = .FALSE.

     
     AnyMultiply = ListCheckPresentAnyMaterial( Model, 'Heat Transfer Multiplier' ) 

!------------------------------------------------------------------------------

     TransientHeaterControl = .FALSE.
     IF(SmartHeaterControl) THEN

       ! Mark the smart heaters 
       SmartHeaters = .FALSE.
       bf_id = 0
       DO i = 1,Model % NumberOfBodyForces
         IF( ListGetLogical( Model % BodyForces(i) % Values, &
             'Smart Heater Control', Found ) ) THEN
           SmartHeaters(i) = .TRUE.	     
           bf_id = i
         END IF
       END DO

       ! Find the BC that controls the heater 
       ! If not found assume that smart heater is related to phase change 
       MeltPoint = GetCReal( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Temperature',GotMeltPoint)           
              
       SmartHeaterAverage = .FALSE.
       SmartHeaterNode = ListGetInteger( Model % BodyForces(bf_id) % Values,&
           'Smart Heater Control Node',GotIt) 
       IF(.NOT. GotIt) THEN
         RealWork => ListGetConstRealArray( Model % BodyForces(bf_id) % Values,&
             'Smart Heater Control Point',GotIt) 
         IF( GotIt ) THEN
           ControlPoint(1:3) = RealWork(1:3,1)
           
           mindist = HUGE( mindist )
           DO l=1,Model % NumberOfNodes
             IF( TempPerm(l) == 0 ) CYCLE
             
             jx = Model % Mesh % Nodes % x(l)
             jy = Model % Mesh % Nodes % y(l)
             jz = Model % Mesh % Nodes % z(l)
             
             dist = (ControlPoint(1)-jx)**2 + (ControlPoint(2)-jy)**2 + (ControlPoint(3)-jz)**2
             IF( dist < mindist ) THEN
               mindist = dist
               SmartHeaterNode = l
             END IF
           END DO
         END IF

         WRITE(Message,*) 'Found Control Point at distance:',SQRT(mindist)
         CALL Info('HeatSolve',Message)
         WRITE(Message,*) 'Control Point Index:',SmartHeaterNode
         CALL Info('HeatSolve',Message)        
       END IF
       
       IF( .NOT. GotMeltPoint .OR. SmartHeaterNode == 0) THEN
         GotIt = .FALSE.
         Found = .FALSE.
         SmartHeaterBC = 0
         
         DO i=1,Model % NumberOfBCs
           GotIt = ListGetLogical( Model % BCs(i) % Values,'Smart Heater Boundary', Found ) 
           IF(GotIt) THEN
             SmartHeaterBC = i
             EXIT
           END IF
         END DO
         IF(.NOT. GotIt) THEN
           DO i=1,Model % NumberOfBCs
             GotIt = ListGetLogical( Model % BCs(i) % Values,'Phase Change', Found ) 
             IF(GotIt) THEN
               SmartHeaterBC = i
               EXIT
             END IF
           END DO
         END IF
         IF(SmartHeaterBC == 0) THEN
           CALL Fatal('HeatSolve','Smart Heater Boundary / Phase Change is undefined')
         END IF
         
         MeltPoint = GetCReal( Model % BCs(SmartHeaterBC) % Values,&
             'Smart Heater Temperature',Found)
         IF(.NOT. Found) THEN
           DO k=1, Model % NumberOfMaterials
             MeltPoint = GetCReal( Model % Materials(k) % Values, &
                 'Melting Point', Found )
             IF(Found) EXIT
           END DO
           IF(.NOT. Found) THEN
             CALL Fatal('HeatSolver','Smart Heater Temperature / Melting Point is undefined')
           END IF
         END IF
         
         ! Find the node related to temperature control 
         SmartHeaterAverage = ListGetLogical(Solver % Values,'Smart Heater Average', Found)
         IF(.NOT. SmartHeaterAverage) THEN
           jx = -HUGE(jx)
           DO k = Model % Mesh % NumberOfBulkElements + 1, &
               Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements
             
             Element => Model % Mesh % Elements(k)
             
             IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
               DO l=1,Element % TYPE % NumberOfNodes
                 IF ( Model % Mesh % Nodes % x(Element % NodeIndexes(l)) >= jx ) THEN
                   j = Element % NodeIndexes(l) 
                   jx = Model % Mesh % Nodes % x(Element % NodeIndexes(l))
                 END IF
               END DO
             END IF
           END DO
           SmartHeaterNode = j
         END IF
       END IF

        SmartTol  = GetConstReal( SolverParams, &
             'Smart Heater Control After Tolerance',  Found )
        IF(.NOT. Found) THEN
          SmartTolReached = .TRUE.
          SmartTol = 1.0
        END IF   
     
        PowerTimeScale = ListGetConstReal(Solver % Values, &
             'Smart Heater Time Scale',Found)

        IF(TransientSimulation .AND. dt < PowerTimeScale) THEN
           TransientHeaterControl = .TRUE.
           CALL Info( 'HeatSolve', 'Using Transient Heater Control')
        ELSE
           TransientHeaterControl = .FALSE.
           CALL Info( 'HeatSolve', 'Using Steady-state Heater Control')
        END IF
        
        IF(Solver % DoneTime /= DoneTime) THEN
           PrevPowerScaling = PowerScaling
           DoneTime = Solver % DoneTime
        END IF
     END IF

     IF( IntegralHeaterControl) THEN
        CALL Info( 'HeatSolve', 'Using Integral Heater Control')       
        IntegralHeaters = .FALSE.
        DO i = 1,Model % NumberOfBodyForces
           IntegralHeaters(i) = ListCheckPresent( Model % BodyForces(i) % Values, &
                'Integral Heat Source')
        END DO
     END IF

!------------------------------------------------------------------------------

     ConstantBulk = GetLogical( SolverParams, 'Constant Bulk System', Found )
     SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Save Bulk System', Found )
     SaveBulk = ConstantBulk .OR. GetLogical( SolverParams, 'Calculate Loads', Found )

!------------------------------------------------------------------------------

     SaveRelax = Relax
     CumulativeTime = 0.0d0
     HeaterControlLocal = .FALSE.

     IF(isRadiation) THEN
       nb = Solver % Mesh % NumberOfBoundaryElements
       ALLOCATE(Areas(nb), Emiss(nb), Reflect(nb) )
       Areas=0; Emiss=0; Reflect=0
     END IF
!------------------------------------------------------------------------------
     FirstTime = .TRUE.

     ALLOCATE(PrevSolution(LocalNodes))
     
     DO WHILE( CumulativeTime < Timestep-1.0d-12 .OR. .NOT. TransientSimulation )
!------------------------------------------------------------------------------
!    The first time around this has been done by the caller...
!------------------------------------------------------------------------------
     IF ( TransientSimulation .AND. .NOT.FirstTime ) &
       CALL InitializeTimestep(Solver)
     FirstTime = .FALSE.
!------------------------------------------------------------------------------
!    Save current solution
!------------------------------------------------------------------------------
     PrevSolution = Temperature(1:LocalNodes)
     IF ( TransientSimulation ) THEN
       PrevTemperature => Solver % Variable % PrevValues(:,1)
     END IF
!------------------------------------------------------------------------------
     
     totat = 0.0d0
     totst = 0.0d0

     Norm = Solver % Variable % Norm

     CALL DefaultStart()
     
     
     DO iter=1,NonlinearIter
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       WRITE( Message,* ) 'TEMPERATURE ITERATION', iter
       CALL Info( 'HeatSolve', Message, Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', 'Starting Assembly...', Level=4 )

500    IF ( ConstantBulk .AND. ASSOCIATED(Solver % Matrix % BulkValues) ) THEN
         Solver % Matrix % Values = Solver % Matrix % BulkValues
         Solver % Matrix % RHS = Solver % Matrix % BulkRHS
         GOTO 1000
       END IF
            
       IF(Radiosity) THEN
         CALL RadiationFactors( Solver, .FALSE., NewtonLinearization)         
         CALL TabulateBoundaryAverages(Solver % Mesh, Emiss, Reflect)         
       ELSE
         CALL TabulateBoundaryAverages(Solver % Mesh, Emiss) 
       END IF
         
!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
 
       IF ( SmartHeaterControl .OR. IntegralHeaterControl ) THEN
          IF( SmartHeaterControl) ForceHeater = 0.0d0
          HeaterArea = 0.0d0
          HeaterSource = 0.0d0
          HeaterScaling = 1.0d0
          HeaterDensity = 0.0d0
          HeaterTarget = 0.0d0
          HeaterControlLocal = .FALSE.

          DO t=1,Solver % NumberOfActiveElements             
             Element => GetActiveElement(t)             
             BodyForce => GetBodyForce()
             
             IF ( .NOT. ASSOCIATED( BodyForce ) ) CYCLE
             bf_id = GetBodyForceId()
             
             IF( .NOT. (SmartHeaters(bf_id) .OR. IntegralHeaters(bf_id) ) ) CYCLE

             n = GetElementNOFNodes()

             Material => GetMaterial()

             Density(1:n) = GetReal( Material, 'Density' )
             Load(1:n) = GetReal( BodyForce, 'Heat Source' )

             s = ElementArea( Solver % Mesh, Element, n )

             IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
                  CurrentCoordinateSystem() == CylindricSymmetric ) s = 2 * PI * s

             HeaterSource(bf_id) = HeaterSource(bf_id) + s * SUM(Density(1:n) * Load(1:n)) / n
             HeaterArea(bf_id) = HeaterArea(bf_id) + s
             HeaterDensity(bf_id) = HeaterDensity(bf_id) + s * SUM( Density(1:n) ) / n
          END DO

          DO i = 1,Model % NumberOfBodyForces
             IF( IntegralHeaters(i) .OR. SmartHeaters(i) ) THEN
                HeaterDensity(i) = HeaterDensity(i) / HeaterArea(i)
             END IF
             IF(IntegralHeaters(i)) THEN
                HeaterTarget(i) = GetCReal(  Model % BodyForces(i) % Values, &
                     'Integral Heat Source', Found )
                HeaterScaling(i) = HeaterTarget(i) / HeaterSource(i)
             END IF
          END DO
       END IF

!------------------------------------------------------------------------------
       body_id = -1
       NULLIFY(Material)
!------------------------------------------------------------------------------
!      Bulk elements
!------------------------------------------------------------------------------
       CALL StartAdvanceOutput( 'HeatSolve', 'Assembly:' )
       NofActive = GetNOFActive()

       DO t=1,NofActive
         
         CALL AdvanceOutput(t,NofActive)
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where temperature 
!        should be calculated
!------------------------------------------------------------------------------
         Element => GetActiveElement(t)

!------------------------------------------------------------------------------
         IF ( Element % BodyId /= body_id ) THEN
!------------------------------------------------------------------------------
           Equation => GetEquation()
           ConvectionFlag = GetString( Equation, 'Convection', Found )

           Material => GetMaterial()
!------------------------------------------------------------------------------
           CompressibilityFlag = GetString( Material, &
                 'Compressibility Model', Found)
           IF ( .NOT.Found ) CompressibilityModel = Incompressible

           SELECT CASE( CompressibilityFlag )

             CASE( 'incompressible' )
               CompressibilityModel = Incompressible

             CASE( 'user defined' )
               CompressibilityModel = UserDefined1

             CASE( 'perfect gas', 'perfect gas equation 1' )
               CompressibilityModel = PerfectGas1

             CASE( 'thermal' )
               CompressibilityModel = Thermal

             CASE DEFAULT
               CompressibilityModel = Incompressible
           END SELECT
!------------------------------------------------------------------------------

           PhaseModel = GetString( Equation, 'Phase Change Model',Found )
           IF(.NOT. Found) PhaseModel = GetString( Material, 'Phase Change Model',Found )

           PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')
           IF ( PhaseChange ) THEN
              CheckLatentHeatRelease = GetLogical( Equation, &
                   'Check Latent Heat Release',Found )
           END IF
         END IF
!------------------------------------------------------------------------------

         n = GetElementNOFNodes()
         CALL GetElementNodes( ElementNodes )

         CALL GetScalarLocalSolution( LocalTemperature )
!------------------------------------------------------------------------------
!        Get element material parameters
!------------------------------------------------------------------------------
         HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity', Found )

         CALL ListGetRealArray( Material,'Heat Conductivity',Hwrk,n, &
                      Element % NodeIndexes )
         HeatConductivity = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
             HeatConductivity( i,i,1:n ) = Hwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
             HeatConductivity(i,i,1:n) = Hwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
               HeatConductivity( i,j,1:n ) = Hwrk(i,j,1:n)
             END DO
           END DO
         END IF
!------------------------------------------------------------------------------

         IF ( CompressibilityModel == PerfectGas1 ) THEN

           ! Read Specific Heat Ratio:
           !--------------------------
           SpecificHeatRatio = GetConstReal( Material, &
               'Specific Heat Ratio', Found )
           IF ( .NOT.Found ) SpecificHeatRatio = 5.d0/3.d0

           ! For an ideal gas, \gamma, c_p and R are really a constant
           ! GasConstant is an array only since HeatCapacity formally is:
           !-------------------------------------------------------------
           GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) * &
               HeatCapacity(1:n) / SpecificHeatRatio

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 1.0_dp
         ELSE IF ( CompressibilityModel == Thermal ) THEN
           ReferenceTemperature(1:n) = GetReal( Material, 'Reference Temperature' )
           HeatExpansionCoeff(1:n) = GetReal( Material, 'Heat Expansion Coefficient' )

           Density(1:n) = GetReal( Material,'Density' )
           Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                ( LocalTemperature(1:n) - ReferenceTemperature(1:n) ) )

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) &
             PressureCoeff(1:n) = LocalTemperature(1:n) * HeatExpansionCoeff(1:n) / &
                   ( 1-HeatExpansionCoeff(1:n)*( &
                               LocalTemperature(1:n)-ReferenceTemperature(1:n)) )
         ELSE IF ( CompressibilityModel == UserDefined1 ) THEN
           IF ( ASSOCIATED( DensitySol ) ) THEN
             CALL GetScalarLocalSolution( Density, 'Density' ) 
           ELSE
             Density(1:n) = GetReal( Material,'Density' )
           END IF
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
         ELSE
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
           Density(1:n) = GetReal( Material, 'Density' )
         END IF

!------------------------------------------------------------------------------
! Take pressure deviation p_d as the dependent variable, p = p_0 + p_d
! for PerfectGas, read p_0
!------------------------------------------------------------------------------
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', Found)
           IF ( .NOT.Found ) ReferencePressure = 0.0d0
         END IF
!------------------------------------------------------------------------------

         HeaterControlLocal = .FALSE.
         Load = 0.0D0
         Pressure = 0.0d0
         dPressuredt = 0.0d0
!------------------------------------------------------------------------------
!        Check for convection model
!------------------------------------------------------------------------------
         C1 = 1.0D0
         U = 0._dp
         V = 0._dp
         W = 0._dp

         MU = 0.0d0
         CALL GetVectorLocalSolution( MU, 'Mesh Velocity' )

         IF ( ConvectionFlag == 'constant' ) THEN

           U(1:n) = GetReal( Material, 'Convection Velocity 1', Found )
           IF ( .NOT. Found ) &
              U(1:n) = GetReal( Equation, 'Convection Velocity 1', Found )
           V(1:n) = GetReal( Material, 'Convection Velocity 2', Found )
           IF ( .NOT. Found ) &
             V(1:n) = GetReal( Equation, 'Convection Velocity 2', Found )
           W(1:n) = GetReal( Material, 'Convection Velocity 3', Found )
           IF ( .NOT. Found ) &
             W(1:n) = GetReal( Equation, 'Convection Velocity 3', Found )

         ELSE IF ( ConvectionFlag == 'computed' ) THEN
           DO i=1,n
             k = FlowPerm(Element % NodeIndexes(i))
             IF ( k > 0 ) THEN
!------------------------------------------------------------------------------
               Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
               SELECT CASE( CompressibilityModel )
                 CASE( PerfectGas1 )
                   Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalTemperature(i) )
               END SELECT
               IF ( TransientSimulation ) THEN
                 dPressureDt(i) = ( FlowSolution(NSDOFs*k) - &
                     FlowSol % PrevValues(NSDOFs*k,1) ) / dt
               END IF
!------------------------------------------------------------------------------

               SELECT CASE( NSDOFs )
               CASE(3)
                 U(i) = FlowSolution( NSDOFs*k-2 )
                 V(i) = FlowSolution( NSDOFs*k-1 )
                 W(i) = 0.0D0

               CASE(4)
                 U(i) = FlowSolution( NSDOFs*k-3 )
                 V(i) = FlowSolution( NSDOFs*k-2 )
                 W(i) = FlowSolution( NSDOFs*k-1 )
               END SELECT
             ELSE
               U(i) = 0.0d0
               V(i) = 0.0d0
               W(i) = 0.0d0
             END IF
           END DO
         ELSE
           IF ( ALL(MU==0) ) C1 = 0.0D0 
         END IF

         HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
 
!------------------------------------------------------------------------------
!        Check if modelling Phase Change with Eulerian approach 
!------------------------------------------------------------------------------
         PhaseSpatial = .FALSE.
         IF (  PhaseChange ) THEN
           CALL EffectiveHeatCapacity()
         END IF

         Viscosity = 0.0d0
!------------------------------------------------------------------------------
!        Add body forces, if any
!------------------------------------------------------------------------------
         BodyForce => GetBodyForce()
         IF ( ASSOCIATED( BodyForce ) ) THEN
           bf_id = GetBodyForceId()
!------------------------------------------------------------------------------
!          Frictional viscous heating
!------------------------------------------------------------------------------
           IF ( GetLogical( BodyForce, 'Friction Heat',Found) ) THEN
              Viscosity(1:n) = GetReal( Material,'Viscosity' )
           END IF
!------------------------------------------------------------------------------
!          Given heat source
!------------------------------------------------------------------------------
           Load(1:n) = GetReal( BodyForce, 'Volumetric Heat Source', Found )
           IF(.NOT. Found ) THEN
             Load(1:n) = Density(1:n) *  GetReal( BodyForce, 'Heat Source', Found )
           END IF
             
           IF ( SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
              IF(  SmartHeaters(bf_id) ) THEN
               HeaterControlLocal = .TRUE.
               IF( TransientHeaterControl ) THEN
                 Load(1:n) = PrevPowerScaling * Load(1:n)
                 HeaterScaling(bf_id) = PrevPowerScaling
               END IF
             END IF
           END IF

           IF ( IntegralHeaterControl ) THEN
              IF( IntegralHeaters(bf_id) ) THEN
                 Load(1:n) = Load(1:n) * HeaterScaling(bf_id) 
              END IF
           END IF
         END IF
	
         C0 = 0.0_dp
!------------------------------------------------------------------------------
! Note at this point HeatCapacity = \rho * c_p OR \rho * (c_p - R)
! and C1 = 0 (diffusion) or 1 (convection)
!------------------------------------------------------------------------------
	 
!------------------------------------------------------------------------------
!          Perfusion (added as suggested by Matthias Zenker)
!------------------------------------------------------------------------------
         IF( ASSOCIATED(BodyForce) ) THEN
           PerfusionRate(1:n) = GetReal( BodyForce, 'Perfusion Rate', Found )
         
           IF ( Found ) THEN
             PerfusionRefTemperature(1:n) = GetReal( BodyForce, 'Perfusion Reference Temperature' )
             PerfusionDensity(1:n) = GetReal( BodyForce, 'Perfusion Density' )
             PerfusionHeatCapacity(1:n) = GetReal( BodyForce, 'Perfusion Heat Capacity' )
             C0(1:n) = PerfusionHeatCapacity(1:n) * PerfusionRate(1:n) * PerfusionDensity(1:n) 
             Load(1:n) = Load(1:n) + C0(1:n) * PerfusionRefTemperature(1:n)           
           END IF
         END IF

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, UseBubbles, Element, n, ElementNodes )

!------------------------------------------------------------------------------
         ELSE
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalTemperature, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, Element, n, ElementNodes )
!------------------------------------------------------------------------------
         END IF
!------------------------------------------------------------------------------

         ! The heat equation may have lower dimensional elements active also.
         ! For example, heat transfer through a pipe could be expressed by 1d elements.
         ! Then the multiplier should be the area of the pipe when included in 3D mesh.
         IF( AnyMultiply ) THEN
           HeatTransferMultiplier = GetCReal( Material, 'Heat Transfer Multiplier', Found )
           IF( Found ) THEN
             MASS = HeatTransferMultiplier * MASS
             STIFF = HeatTransferMultiplier * STIFF
             FORCE = HeatTransferMultiplier * FORCE
           END IF
         END IF


         IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN

           IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
             CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF

           CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
               ForceHeater, FORCE, n, 1, TempPerm(Element % NodeIndexes) )
         ELSE
            Bubbles = UseBubbles .AND. .NOT.Stabilize .AND. &
            ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )
            
!------------------------------------------------------------------------------
!           If time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
            TimeForce  = 0.0_dp
            IF ( TransientAssembly ) THEN
               IF ( ConstantBulk ) THEN
                 CALL DefaultUpdateMass( MASS )
               ELSE
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
               END IF
            ELSE IF ( Solver % NOFEigenValues>0 ) THEN
              CALL DefaultUpdateDamp(MASS)
            END IF
!------------------------------------------------------------------------------
!           Update global matrices from local matrices
!------------------------------------------------------------------------------
            IF (  Bubbles ) THEN
               CALL Condensate( N, STIFF, FORCE, TimeForce )
            END IF

            CALL DefaultUpdateEquations( STIFF, FORCE )
         END IF
!------------------------------------------------------------------------------
      END DO     !  Bulk elements
!------------------------------------------------------------------------------
      
      CALL DefaultFinishBulkAssembly()


1000  CONTINUE

     

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      DO bc_elem = 1, Solver % Mesh % NumberOfBoundaryElements
        
        Element => GetBoundaryElement(bc_elem)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

        ! This checks whether there are any Dirichlet conditions on the 
        ! smart heater boundary. If there are the r.h.s. must be zero as 
        ! there can possibly not be any effect on temperature.
        !-----------------------------------------------------------------
	IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN
          IF( ListCheckPresent(BC, Varname) ) THEN
             nd = GetElementDOFs(Indexes)
             ForceHeater(TempPerm(Indexes(1:nd))) = 0.0_dp
          END IF
        END IF

        HeatFluxBC = GetLogical( BC, 'Heat Flux BC', Found )
        IF ( Found .AND. .NOT. HeatFluxBC ) CYCLE

        HeatGapBC = ListGetLogical( BC, 'Heat Gap', Found )
        CALL AddHeatFluxBC()

        IF ( HeatGapBC ) THEN
          CALL FindGapIndexes( Element, Indexes, n )
          SaveIndexes(1:n) = Element % NodeIndexes
          Element % NodeIndexes = Indexes(1:n)
          CALL AddHeatFluxBC()
          Element % NodeIndexes = SaveIndexes(1:n)
        END IF

      END DO   ! Neumann & Newton BCs
!------------------------------------------------------------------------------


      IF ( TransientSimulation .AND. ConstantBulk ) CALL AddGlobalTime()

      CALL DefaultFinishBoundaryAssembly()
      CALL DefaultFinishAssembly()
      CALL Info( 'HeatSolve', 'Assembly done', Level=4 )

      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      at = CPUTime() - at
      st = CPUTime()

      PrevNorm = Norm

      IF(SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
      
        IF(.NOT. TransientHeaterControl) THEN
          ! These are control loops. Do use them to check convergence or advance the
          ! nonlinear iteration flag.
          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.TRUE.)
          CALL ListAddLogical(SolverParams,'Skip Advance Nonlinear Iter',.TRUE.)

          Relax = GetCReal( SolverParams,'Nonlinear System Relaxation Factor', Found )
          
          IF ( Found .AND. Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,&
                'Nonlinear System Relaxation Factor', 1.0d0 )
          ELSE
            Relax = 1.0d0
          END IF          

          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              ForceHeater, XX, Norm, 1, Solver )         
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, YY, Norm, 1, Solver )

          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.FALSE.)
          CALL ListAddLogical(SolverParams,'Skip Advance Nonlinear Iter',.FALSE.)
        ELSE                    
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, Temperature, Norm, 1, Solver )
          YY = Temperature
        END IF

        IF(.NOT. SmartHeaterAverage) THEN
          xave = XX(TempPerm(SmartHeaterNode))
          yave = YY(TempPerm(SmartHeaterNode))
        ELSE          
          xave = 0.0d0
          yave = 0.0d0
          j = 0
          
          DO k = Model % Mesh % NumberOfBulkElements + 1, &
              Model % Mesh % NumberOfBulkElements + Model % Mesh % NumberOfBoundaryElements            

            Element => Model % Mesh % Elements(k)            
            IF ( Element % BoundaryInfo % Constraint == SmartHeaterBC ) THEN
              l = Element % TYPE % NumberOfNodes
              j = j + l
              xave = xave + SUM( XX(TempPerm(Element % NodeIndexes)) )
              yave = yave + SUM( YY(TempPerm(Element % NodeIndexes)) )
            END IF
          END DO
          xave = xave / j
          yave = yave / j 
          CALL ListAddConstReal(Model % Simulation,'res: Smart Heater Temperature',yave)
        END IF

        IF(.NOT. TransientHeaterControl) THEN
          IF ( ASSOCIATED(Solver % Variable % NonlinValues) ) THEN
            Solver % Variable % NonlinValues = Temperature
          END IF

          PowerScaling = (MeltPoint - yave) / xave 
          Temperature = YY + PowerScaling * XX

          ! The change is computed separately for the controlled temperature field
          !-----------------------------------------------------------------------
          CALL ComputeChange(Solver,.FALSE.,LocalNodes,Temperature)
          Norm = Solver % Variable % Norm

        END IF

        IF(dt > PowerTimeScale) THEN
          IF ( Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', Relax )
          END IF
        END IF
      ELSE
!------------------------------------------------------------------------------
!     Check stepsize for nonlinear iteration
!------------------------------------------------------------------------------
        IF( DefaultLinesearch( Converged ) ) GOTO 500
        IF( Converged ) EXIT

        Norm = DefaultSolve()
      END IF


      IF( SmartHeaterControl .OR. IntegralHeaterControl) THEN
         CALL Info( 'HeatSolve', 'Heater Control Information', Level=4 )
         DO i=1,Model % NumberOfBodyForces
            IF( .NOT. (SmartHeaters(i) .OR. IntegralHeaters(i))) CYCLE
            IF( SmartHeaters(i) )  HeaterScaling(i) = PowerScaling

            WRITE( Message, '(A,T35,I15)' ) 'Heater for body: ', i
            CALL Info( 'HeatSolve', Message, Level=4 )
            IF(SmartHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Smart heater'
            IF(IntegralHeaters(i)) WRITE( Message, '(A,T35,A)' ) 'Heater type:','Integral heater'
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Volume (m^3): ', HeaterArea(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            s = HeaterSource(i) * HeaterScaling(i)
            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater Power (W): ', s
            CALL Info( 'HeatSolve', Message, Level=4 )

            WRITE( Message,'(A,T35,ES15.4)' ) 'Heater scaling: ', HeaterScaling(i)
            CALL Info( 'HeatSolve', Message, Level=4 )
            WRITE( Message, '(A,T35,ES15.4)' ) 'Heater Power Density (W/kg): ', s/(HeaterDensity(i) * HeaterArea(i))
            CALL Info( 'HeatSolve', Message, Level=4 )

            CALL ListAddConstReal(Model % Simulation,'res: Heater Power Scaling '//I2S(i),HeaterScaling(i))
            CALL ListAddConstReal(Model % Simulation,'res: Heater Power Density '//I2S(i),&
                 s/(HeaterDensity(i) * HeaterArea(i)))
         END DO
      END IF


      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'HeatSolve', Message, Level=4 )
!------------------------------------------------------------------------------
!     If modelling phase change (and if requested by the user), check if any
!     node has jumped over the phase change interval, and if so, reduce
!     timestep and or relaxation and recompute.
!------------------------------------------------------------------------------
      IF (PhaseChange .AND. CheckLatentHeatRelease ) THEN
!------------------------------------------------------------------------------
        IF ( CheckLatentHeat() ) THEN
          Temperature(1:LocalNodes) = PrevSolution
          Norm = PrevNorm

          IF ( TransientSimulation ) THEN
            dt = dt / 2
            Solver % dt = dt
            WRITE( Message, * ) &
                  'Latent heat release check: reducing timestep to: ',dt
            CALL Info( 'HeatSolve', Message, Level=4 )
          ELSE
            Relax = Relax / 2
            CALL  ListAddConstReal( Solver % Values,  &
                 'Nonlinear System Relaxation Factor', Relax )
            WRITE( Message, * ) &
                 'Latent heat release check: reducing relaxation to: ',Relax
            CALL Info( 'HeatSolve', Message, Level=4 )
          END IF

          CYCLE
        END IF
        IF ( .NOT.TransientSimulation ) PrevSolution=Temperature(1:LocalNodes)
      END IF
!------------------------------------------------------------------------------
     
      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'HeatSolve', Message, Level=4 )

      IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
               NewtonLinearization = .TRUE.
      Converged =  ( Solver % Variable % NonlinConverged > 0 ) .AND. &
          ( .NOT. SmartHeaterControl .OR. SmartTolReached )
      IF( Converged ) EXIT

      IF(SmartHeaterControl) THEN
        IF ( RelativeChange < SmartTol ) THEN
          SmartTolReached = .TRUE.
          YY = Temperature
        END IF
      END IF
      
!------------------------------------------------------------------------------
    END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------
    IF(TransientHeaterControl) THEN
      PowerRelax = GetCReal(Solver % Values,'Smart Heater Relaxation Factor', GotIt)
      IF(.NOT. GotIt) PowerRelax = 1.0_dp
      PowerSensitivity = ListGetConstReal(Solver % Values,'Smart Heater Power Sensivity',GotIt)
      IF(.NOT. GotIt) PowerSensitivity = 4.0_dp
      PowerScaling = PowerScaling * (1 + PowerSensitivity * PowerRelax * (MeltPoint/yave - 1.0d0) ) 

      IF( ListGetLogical( Solver % Values,'Smart Heater Transient Speedup',GotIt ) ) THEN
        Temperature = Temperature * (1 + PowerRelax * (MeltPoint/yave - 1.0d0)   )     
      END IF
      YY = Temperature
    END IF

!------------------------------------------------------------------------------
!   Compute cumulative time done by now and time remaining
!------------------------------------------------------------------------------
    IF ( .NOT. TransientSimulation ) EXIT
    CumulativeTime = CumulativeTime + dt
    dt = Timestep - CumulativeTime

   END DO ! time interval
   Solver % dt = Timestep

   CALL DefaultFinish()
   
!------------------------------------------------------------------------------
   CALL  ListAddConstReal( Solver % Values,  &
        'Nonlinear System Relaxation Factor', SaveRelax )
!------------------------------------------------------------------------------

   DEALLOCATE( PrevSolution )

   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) THEN
     IF(.NOT. ListGetLogical( Solver % Values,'Library Adaptivity',Found )) THEN
       CALL RefineMesh( Model,Solver,Temperature,TempPerm, &
           HeatSolver_Inside_Residual, HeatSolver_Edge_Residual, &
           HeatSolver_Boundary_Residual )
     END IF
   END IF
     
CONTAINS


!------------------------------------------------------------------------------
! To save some time tabulate data needed for the diffuse gray radiation. 
!------------------------------------------------------------------------------
  SUBROUTINE TabulateBoundaryAverages( Mesh, Emiss, Reflect )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), ALLOCATABLE :: Emiss(:)
    REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: Reflect(:)
 !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    INTEGER :: bindex, nb, n, j, noactive
    INTEGER :: ElemInds(12)
    REAL(KIND=dp) :: NodalVal(12), Absorp
    
     nb = Mesh % NumberOfBoundaryElements
     NoActive = 0
     
     DO j=1,nb
       bindex = j + Mesh % NumberOfBulkElements
       Element => Mesh % Elements(bindex)

       BC => GetBC(Element)
       IF(.NOT. ASSOCIATED( BC ) ) CYCLE

       IF( .NOT. ListCheckPresent( BC,'Radiation' ) ) CYCLE

       IF(.NOT. ALLOCATED( Emiss ) ) THEN
         ALLOCATE( Emiss(nb) )
         Emiss = 0.0_dp
         IF( PRESENT( Reflect ) ) THEN
           ALLOCATE( Reflect(nb) )
           Reflect = 0.0_dp
         END IF
       END IF

       n = GetElementNOFNodes(Element)

       NodalVal(1:n) = GetReal(BC,'Emissivity',Found)
       IF (Found) THEN
         Emiss(j) = SUM(NodalVal(1:n)) / n
         IF( PRESENT(Reflect)) THEN
           NodalVal(1:n) = GetReal(BC,'Absorptivity',Found)
           IF(Found) THEN
             Absorp = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp = Emiss(j)
           END IF
           NodalVal(1:n) = GetReal(BC,'Reflectivity',Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp
           END IF
         END IF
       ELSE
         NodalVal(1:n) = GetParentMatProp('Emissivity',Element)
         Emiss(j) = SUM(NodalVal(1:n)) / n
         IF( PRESENT( Reflect ) ) THEN
           NodalVal(1:n) = GetParentMatProp('Absorptivity',Element, Found)
           IF(Found) THEN
             Absorp = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp = Emiss(j)
           END IF
           NodalVal(1:n) = GetParentMatProp('Reflectivity',Element, Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp
           END IF
         END IF
       END IF
     END DO

   END SUBROUTINE TabulateBoundaryAverages
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE AddHeatFluxBC()
!------------------------------------------------------------------------------
      CALL GetElementNodes( ElementNodes )

      HeatTransferCoeff = 0.0D0
      LOAD  = 0.0D0
!------------------------------------------------------------------------------
!     BC: -k@T/@n = \epsilon\sigma(T^4 - Text^4)
!------------------------------------------------------------------------------
      RadiationFlag = GetString( BC, 'Radiation', Found )

      IF ( Found .AND. RadiationFlag /= 'none' ) THEN
        Emissivity = Emiss(bc_elem)

!------------------------------------------------------------------------------
        IsRadiosity = .FALSE.
        IF (  RadiationFlag == 'idealized' ) THEN
          AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
          IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
        ELSE
          IF( Radiosity ) THEN
            CALL RadiosityRadiation( Model, Solver, Element, & 
                n, Temperature, TempPerm, ForceVector )
            IsRadiosity = .TRUE.
          ELSE          
            CALL DiffuseGrayRadiation( Model, Solver, Element, & 
                Temperature, TempPerm, ForceVector, VisibleFraction, Text)
            
            IF( GetLogical( BC, 'Radiation Boundary Open', Found) ) THEN
              AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
              IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )
              IF( VisibleFraction >= 1.0_dp ) THEN
                Atext(1:n) = Text
              ELSE
                Atext(1:n) = ( (1 - VisibleFraction) * Atext(1:n)**4 + &
                    VisibleFraction * Text**4 ) ** 0.25_dp
              END IF
            ELSE
              AText(1:n) = Text
            END IF
          END IF
        END IF
!------------------------------------------------------------------------------
!       Add our own contribution to surface temperature (and external
!       if using linear type iteration or idealized radiation)
!------------------------------------------------------------------------------
        IF(.NOT. IsRadiosity ) THEN
          DO j=1,n
            k = TempPerm(Element % NodeIndexes(j))
            Text = AText(j)
            IF ( .NOT. HeatGapBC .AND. NewtonLinearization ) THEN
              HeatTransferCoeff(j) = Emissivity * 4*Temperature(k)**3 * &
                  StefanBoltzmann
              LOAD(j) = Emissivity*(3*Temperature(k)**4+Text**4) * &
                  StefanBoltzmann
            ELSE
              HeatTransferCoeff(j) = Emissivity * (Temperature(k)**3 +   &
                  Temperature(k)**2*Text+Temperature(k)*Text**2 + Text**3) * &
                  StefanBoltzmann 
              LOAD(j) = HeatTransferCoeff(j) * Text
            END IF
          END DO
        END IF
      END IF  ! of radition
!------------------------------------------------------------------------------

      Work(1:n)  = GetReal( BC, 'Heat Transfer Coefficient',Found )
      IF ( Found ) THEN
       AText(1:n) = GetReal( BC, 'External Temperature',Found )
       DO j=1,n
!------------------------------------------------------------------------------
!         BC: -k@T/@n = \alpha(T - Text)
!------------------------------------------------------------------------------
          k = TempPerm(Element % NodeIndexes(j))
          LOAD(j) = LOAD(j) + Work(j) * AText(j)
          HeatTransferCoeff(j) = HeatTransferCoeff(j) + Work(j)
        END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = (rho*L)*v.n 
!     Heating related to pulling is possible only in ss cases where pull velocity
!     is desrcibed.
!------------------------------------------------------------------------------

      IF( GetLogical( BC, 'Phase Change',Found ) ) THEN
         PhaseVelocity(1,1:n) = GetReal( BC,'Phase Velocity 1', Found  )
         PhaseVelocity(2,1:n) = GetReal( BC,'Phase Velocity 2', Found  )
         PhaseVelocity(3,1:n) = GetReal( BC,'Phase Velocity 3', Found  )
  
         ! Ensure that the latent heat and density come from the same side
         LatentHeat(1:n) = GetParentMatProp( 'Latent Heat', &
              UElement = Element, UParent = Parent )
         IF(.NOT. ASSOCIATED(Parent) ) THEN
           CALL Warn('HeatSolve','Parent not associated')
         ELSE
           k = GetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
           Density(1:n) = GetReal( Model % Materials(k) % Values, 'Density' )
         END IF

         ! This could be rather put as a new type of BC into the assembly routine and 
         ! then the Normal could be taken at the proper Gaussian integration points. 
         Normal = NormalVector( Element, ElementNodes, 0.0_dp, 0.0_dp, .TRUE. )

         DO i=1,n
            LOAD(i) = LOAD(i) + &
                 LatentHeat(i) * Density(i) * SUM( Normal(1:3) * PhaseVelocity(1:3,i))
         END DO
      END IF

!------------------------------------------------------------------------------
!     BC: -k@T/@n = g
!------------------------------------------------------------------------------
      LOAD(1:n) = LOAD(1:n) +  GetReal( BC, 'Heat Flux', Found )

      InfBC = ListGetLogical( BC,'Infinity BC '//TRIM(VarName),GotIt)
      IF( InfBC ) THEN
        AText(1:n) = GetReal( BC,'Infinity BC '//TRIM(VarName)//' Offset',GotIt)
        ! currently only isotropic heat conductivity supported
        HeatConductivityIso(1:n) = GetParentMatProp('Heat Conductivity',Element,GotIt)
        IF(.NOT. GotIt) THEN
          CALL Fatal( 'HeatSolver','Could not find > Heat Conductivity < for parent!' )           
        END IF
      END IF


!------------------------------------------------------------------------------
!     Get element matrix and rhs due to boundary conditions ...
!------------------------------------------------------------------------------
      IF ( CurrentCoordinateSystem() == Cartesian ) THEN
        CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
            LOAD,HeatTransferCoeff,InfBC,HeatConductivityIso,AText(1:n),&
            Element,n,ElementNodes )
      ELSE
        IF( InfBC ) THEN
          CALL Fatal('HeatSolver','Infinity BC not implemented only for cartersian case!')
        END IF
        CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
            LOAD,HeatTransferCoeff,Element,n,ElementNodes ) 
      END IF

!------------------------------------------------------------------------------
!     Update global matrices from local matrices
!------------------------------------------------------------------------------
      IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
        MASS = 0.d0
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF

      IF ( HeatGapBC ) &
        CALL AddHeatGap( Solver, Element, STIFF, TempPerm)

      CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatFluxBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddGlobalTime()
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,n
      REAL(KIND=dp) :: FORCE(1)
      REAL(KIND=dp), POINTER :: SaveValues(:) => NULL()
      SAVE STIFF, MASS, X
      REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),MASS(:,:), X(:,:)

      IF ( .NOT.ASSOCIATED(Solver % Variable % Values, SaveValues) ) THEN
         IF ( ALLOCATED(STIFF) ) DEALLOCATE( STIFF,MASS,X )
         n = 0
         DO i=1,Solver % Matrix % NumberOfRows
           n = MAX( n,Solver % Matrix % Rows(i+1)-Solver % Matrix % Rows(i) )
         END DO
         k = SIZE(Solver % Variable % PrevValues,2)
         ALLOCATE( STIFF(1,n),MASS(1,n),X(n,k) )
 
         SaveValues => Solver % Variable % Values
      END IF

      DO i=1,Solver % Matrix % NumberOFRows
        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
          n=n+1
          STIFF(1,n) = Solver % Matrix % Values(j)
          MASS(1,n)  = Solver % Matrix % MassValues(j)
          X(n,:) = Solver % Variable % PrevValues(Solver % Matrix % Cols(j),:)
        END DO
        FORCE(1) = Solver % Matrix % RHS(i)
        Solver % Matrix % Force(i,1) = FORCE(1)
        k = MIN( Solver % DoneTime, Solver % Order )
        CALL BDFLocal( n, dt, MASS, STIFF, FORCE, X, k )

        n = 0
        DO j=Solver % Matrix % Rows(i),Solver % Matrix % Rows(i+1)-1
           n=n+1
          Solver % Matrix % Values(j) = STIFF(1,n)
        END DO
        Solver % Matrix % RHS(i) = FORCE(1)
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddGlobalTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE DiffuseGrayRadiation( Model, Solver, Element,  &
      Temperature, TempPerm, ForceVector,AngleFraction, Text)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: TempPerm(:)
      REAL(KIND=dp) :: Temperature(:), ForceVector(:)
      REAL(KIND=dp) :: AngleFraction, Text
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Area, Asum, gEmissivity, Base(12), Load, TransCoeff
      REAL(KIND=dp), POINTER :: Fact(:)
      INTEGER :: i,j,k,l,m,ImplicitFactors, nf,nr, bindex, nb
      INTEGER, POINTER :: ElementList(:)
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------

      Asum = 0.0_dp
      IF ( .NOT. NewtonLinearization ) THEN
        Text = ComputeRadiationLoad( Model, Solver % Mesh, Element, &
           Temperature, TempPerm, Emissivity, AngleFraction, Areas, Emiss )

      ELSE   !  Full Newton-Raphson solver
!------------------------------------------------------------------------------
!       Go through surfaces (j) this surface (i) is getting
!       radiated from.
!------------------------------------------------------------------------------

        ElementList => Element % BoundaryInfo % RadiationFactors % Elements

        nf = Element % BoundaryInfo % RadiationFactors % NumberOfFactors

        bindex = Element % ElementIndex - Solver % Mesh % NumberOfBulkElements
        Area  = Areas(bindex)
        CALL GetBase( Base, Element, n, ElementNodes )

        Fact => Element % BoundaryInfo % RadiationFactors % Factors

        DO j=1,nf
          RadiationElement => Solver % Mesh % Elements(ElementList(j))
          Text = Fact(j)
          Asum = Asum + Text

!------------------------------------------------------------------------------
!         Gebhart factors are given elementwise at the center
!         of the element, so take average of nodal temperatures
!         (or integrate over surface j)
!------------------------------------------------------------------------------

          k = RadiationElement % TYPE % NumberOfNodes
          ImplicitFactors = Element % BoundaryInfo % RadiationFactors % NumberOfImplicitFactors
          IF(ImplicitFactors == 0) &
              ImplicitFactors = Element % BoundaryInfo % RadiationFactors % NumberOfFactors

          IF(j <= ImplicitFactors) THEN
            
            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )**(1._dp/4._dp)
!------------------------------------------------------------------------------
!          Linearization of the G_jiT^4_j term
!------------------------------------------------------------------------------
           LOAD = -3 * Text * S**4 * StefanBoltzmann
           TransCoeff = -4 * Text * S**3 * StefanBoltzmann
!------------------------------------------------------------------------------
!          Integrate the contribution of surface j over surface i
!          and add to global matrix
!------------------------------------------------------------------------------
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              DO l=1,k
                k2 = TempPerm( RadiationElement % NodeIndexes(l) )
                CALL AddToMatrixElement( StiffMatrix,k1,k2,TransCoeff*Base(m)/k )
              END DO
              ForceVector(k1) = ForceVector(k1) + Load*Base(m)
            END DO

          ELSE
            S = (SUM( Temperature( TempPerm( RadiationElement % &
                NodeIndexes))**4 )/k )
            
            LOAD = Text * S * StefanBoltzmann
            
            DO m=1,n
              k1 = TempPerm( Element % NodeIndexes(m) )
              ForceVector(k1) = ForceVector(k1) + LOAD*Base(m)
            END DO            
          END IF 

        END DO

!------------------------------------------------------------------------------
!       We have already added all external temperature contributions
!       to the matrix for the Newton type iteration
!------------------------------------------------------------------------------
        AngleFraction = Asum / Emissivity
        Text = 0.0

      END IF  !  of newton-raphson

    END SUBROUTINE DiffuseGrayRadiation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE RadiosityRadiation( Model, Solver, Element,  &
      n, Temperature, TempPerm, ForceVector)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: n
      INTEGER :: TempPerm(:)
      REAL(KIND=dp) :: Temperature(:), ForceVector(:)      
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Emis1, Refl1, RadCoeffAtIp, RadLoadAtIp, TempAtIp, s, x, y, z, DetJ
      REAL(KIND=dp) :: Basis(27)
      REAL(KIND=dp), POINTER :: Fact(:)
      INTEGER :: t,p,q,k1,k2
      INTEGER, POINTER :: pIndexes(:)
      TYPE(GaussIntegrationPoints_t), TARGET :: IP
      LOGICAL :: stat
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------
      
      pIndexes => Element % NodeIndexes     
      
      Emis1 = Emissivity            
      IF(.NOT.Spectral ) THEN
        Refl1 = Reflect(bc_elem)
        Emis1 = Emis1 / Refl1
      END IF
      Fact => Element % BoundaryInfo % RadiationFactors % Factors
      
      TempAtIp = SUM(Temperature(TempPerm(pIndexes(1:n))))/n

      IF(NewtonLinearization) THEN
        RadLoadAtIp =  3 * Emis1 * TempAtIp**4 * StefanBoltzmann + &
             (Fact(1) - Fact(2)*TempAtIp)
        RadCoeffAtIp = 4 * Emis1 * TempAtIp**3 * StefanBoltzmann - Fact(2)
      ELSE
        RadCoeffAtIp = Emis1 * StefanBoltzmann * TempAtIp**3
        RadLoadAtIp = Fact(1)
      END IF

      IP = GaussPoints( Element )
      
      DO t=1,IP % n
        stat = ElementInfo( Element,ElementNodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
        s = detJ * IP % s(t)
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( ElementNodes % x(1:n)*Basis )
          y = SUM( ElementNodes % y(1:n)*Basis )
          z = SUM( ElementNodes % z(1:n)*Basis )
          s = s * CoordinateSqrtMetric( x,y,z )
        END IF
        DO p=1,n
          k1 = TempPerm( pIndexes(p) )
          DO q=1,n
            k2 = TempPerm( pIndexes(q) )
            CALL AddToMatrixElement( Solver % Matrix,k1,k2,s*Basis(p)*Basis(q)*RadCoeffAtIP)
          END DO
          ForceVector(k1) = ForceVector(k1) + s *Basis(p) * RadLoadAtIp
        END DO
      END DO
      
    END SUBROUTINE RadiosityRadiation
!------------------------------------------------------------------------------


    
!------------------------------------------------------------------------------
    SUBROUTINE EffectiveHeatCapacity()
      LOGICAL :: Found, Specific, GotFraction
      REAL(KIND=dp), ALLOCATABLE :: dT(:)
      REAL(KIND=dp) :: dT0

!------------------------------------------------------------------------------
!     See if temperature gradient indside the element is large enough 
!     to use  the c_p = SQRT( (dH/dx)^2 / (dT/dx)^2 ), otherwise
!     use c_p = dH/dT, or if in time dependent simulation, use
!     c_p = (dH/dt) / (dT/dt), if requested. 
!------------------------------------------------------------------------------

      SELECT CASE(PhaseModel)
!------------------------------------------------------------------------------
        CASE( 'spatial 1' )
          PhaseChangeModel = PHASE_SPATIAL_1
!------------------------------------------------------------------------------

        CASE( 'spatial 2' )
!------------------------------------------------------------------------------
! Check if local variation of temperature is large enough to actually use the
! Spatial 2 model. Should perhaps be scaled to element size (or actually
! compute the gradient, but this will do for now...).
!------------------------------------------------------------------------------
          s = MAXVAL(LocalTemperature(1:n))-MINVAL(LocalTemperature(1:n))
          IF ( s < AEPS ) THEN
            PhaseChangeModel = PHASE_SPATIAL_1
          ELSE
            PhaseChangeModel = PHASE_SPATIAL_2
          END IF

!------------------------------------------------------------------------------
! Note that here HeatCapacity is miused for saving dT.
!------------------------------------------------------------------------------
        CASE('temporal')
          IF ( TransientSimulation )  THEN
            ALLOCATE( dT(n) )
            dT(1:n) = Temperature(TempPerm(Element % NodeIndexes)) - &
                     PrevTemperature(TempPerm(Element % NodeIndexes))

            IF ( ANY(ABS(dT(1:n)) < AEPS) ) THEN
              PhaseChangeModel = PHASE_SPATIAL_1
            ELSE
              PhaseChangeModel = PHASE_TEMPORAL
            END IF
          ELSE
             PhaseChangeModel = PHASE_SPATIAL_1
          END IF

!------------------------------------------------------------------------------
        CASE DEFAULT
          PhaseChangeModel = PHASE_SPATIAL_1

      END SELECT
!------------------------------------------------------------------------------

      PhaseSpatial = ( PhaseChangeModel == PHASE_SPATIAL_2 )
      Specific = ListCheckPresent( Material,'Specific Enthalpy')

      EnthalpyFraction(1:n) = ListGetReal(Material,'Enthalpy Fraction',&
          n,Element % NodeIndexes,GotFraction)          
 

!-----------------------------------------------------------------------------
      SELECT CASE( PhaseChangeModel )

!------------------------------------------------------------------------------
! This phase change model is available only for some type of real entries 
! that have an implemented analytical derivation rule.
!-----------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_1 )

        Work(1:n) = ListGetReal( Material, &
            'Effective Heat Capacity', n,Element % NodeIndexes, Found )
        IF ( .NOT. Found ) THEN
          dT0 = ListGetCReal( Material,'Enthalpy Temperature Differential',Found )
          IF(.NOT. Found) dT0 = 1.0d-3
          IF( Specific ) THEN
            Work(1:n) = ListGetDerivValue( Material, &
                'Specific Enthalpy', n,Element % NodeIndexes, dT0 )
            Work(1:n) = Density(1:n) * Work(1:n)
          ELSE
            Work(1:n) = ListGetDerivValue( Material, &
                'Enthalpy', n,Element % NodeIndexes, dT0 )
          END IF
        END IF

        IF( GotFraction ) THEN
          HeatCapacity(1:n) = HeatCapacity(1:n) + EnthalpyFraction(1:n) * Work(1:n) 
        ELSE
          HeatCapacity(1:n) = HeatCapacity(1:n) + Work(1:n) 
        END IF
          
!---------------------------------------------------------------------------------------
! Note that for the 'spatial 2' model the evaluation of c_p is done in each integration
! point and thus Enthalphy and PhaseSpatial flag are used instead of HeatCapacity directly.
!-----------------------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_2 )
        IF( Specific ) THEN
          Enthalpy(1:n) = ListGetReal(Material,'Specific Enthalpy',n,Element % NodeIndexes)
          Enthalpy(1:n) = Density(1:n) * Enthalpy(1:n)
        ELSE
          Enthalpy(1:n) = ListGetReal(Material,'Enthalpy',n,Element % NodeIndexes)          
        END IF
        
        IF( GotFraction ) THEN
          CALL Warn('EffectiveHeatCapacity',&
              '> Enthalpy Fraction < not treated yet by spatial 2 phase change')
        END IF

          
!------------------------------------------------------------------------------
      CASE( PHASE_TEMPORAL )
        ! When retrieving the value of enthalphy on the previous timestep 
        ! the relevant entries of the Temperature solution in the global vector
        ! are tampered in order to make the ListGetReal command work as wanted. 
        ! 1) Values at current temperature     
        !------------------------------------------------------------------------
        IF( Specific ) THEN
          Work(1:n) = ListGetReal( Material,'Specific Enthalpy',n,Element % NodeIndexes )
        ELSE
          Work(1:n) = ListGetReal( Material,'Enthalpy',n,Element % NodeIndexes )
        END IF

        ! 2) Values at previous temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) 

        IF( Specific ) THEN
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Specific Enthalpy', &
              n,Element % NodeIndexes )          
          Work(1:n) = Density(1:n) * Work(1:n) / dT(1:n)
        ELSE
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Enthalpy', &
              n,Element % NodeIndexes )
          Work(1:n) = Work(1:n) / dT(1:n)
        END IF

        IF( GotFraction ) THEN
          HeatCapacity(1:n) = HeatCapacity(1:n) + EnthalpyFraction(1:n) * Work(1:n) 
        ELSE
          HeatCapacity(1:n) = HeatCapacity(1:n) + Work(1:n) 
        END IF


        ! Revert to current temperature
        Temperature(TempPerm(Element % NodeIndexes)) = & 
            PrevTemperature(TempPerm(Element % NodeIndexes)) + dT(1:n)        

!------------------------------------------------------------------------------
      END SELECT

!------------------------------------------------------------------------------
    END SUBROUTINE EffectiveHeatCapacity
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    FUNCTION CheckLatentHeat() RESULT(Failure)
!------------------------------------------------------------------------------
      LOGICAL :: Failure, PhaseChange, CheckLatentHeatRelease
      INTEGER :: t, eq_id, body_id
      CHARACTER(LEN=MAX_NAME_LEN) :: PhaseModel
!------------------------------------------------------------------------------

      Failure = .FALSE.
!------------------------------------------------------------------------------
      DO t=1,Solver % Mesh % NumberOfBulkElements
!------------------------------------------------------------------------------
!       Check if this element belongs to a body where temperature 
!       has been calculated
!------------------------------------------------------------------------------
        Element => Solver % Mesh % Elements(t)

        NodeIndexes => Element % NodeIndexes
        IF ( ANY( TempPerm( NodeIndexes ) <= 0 ) ) CYCLE

        body_id = Element % Bodyid
        eq_id = ListGetInteger( Model % Bodies(body_id) % Values, &
            'Equation', minv=1, maxv=Model % NumberOfEquations )

        PhaseModel = ListGetString( Model % Equations(eq_id) % Values, &
                          'Phase Change Model',Found )

        PhaseChange = Found .AND. (PhaseModel(1:4) /= 'none')

        IF ( PhaseChange ) THEN
          CheckLatentHeatRelease = ListGetLogical(Model % Equations(eq_id) % &
                    Values, 'Check Latent Heat Release',Found )
        END IF
        IF ( .NOT. ( PhaseChange .AND. CheckLatentHeatRelease ) ) CYCLE

        n = Element % TYPE % NumberOfNodes
!------------------------------------------------------------------------------
!       Set the current element pointer in the model structure to
!       reflect the element being processed
!------------------------------------------------------------------------------
        Model % CurrentElement => Element
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!       Get element material parameters
!------------------------------------------------------------------------------
        k = ListGetInteger( Model % Bodies(body_id) % Values,'Material', &
                minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(k) % Values

        PhaseChangeIntervals => ListGetConstRealArray( Material, &
                        'Phase Change Intervals' )

        DO k=1,n
          i = TempPerm( NodeIndexes(k) )
          DO j=1,SIZE(PhaseChangeIntervals,2)
            IF ( ( Temperature(i)  < PhaseChangeIntervals(1,j) .AND. &
                   PrevSolution(i) > PhaseChangeIntervals(2,j) ).OR. &
                 ( Temperature(i)  > PhaseChangeIntervals(2,j) .AND. &
                   PrevSolution(i) < PhaseChangeIntervals(1,j) )  ) THEN
              Failure = .TRUE.
              EXIT
            END IF
          END DO
          IF ( Failure ) EXIT
        END DO
        IF ( Failure ) EXIT
      END DO
!------------------------------------------------------------------------------
    END FUNCTION CheckLatentHeat
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
   SUBROUTINE GetBase( Base, Element, n, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Base(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n,  m

     REAL(KIND=dp) :: Basis(n), DetJ

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     Base = 0._dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
!------------------------------------------------------------------------------
       DO p=1,N
         Base(p) = Base(p) + s * Basis(p)
       END DO
     END DO
   END SUBROUTINE GetBase
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE IntegOverA( BoundaryMatrix, BoundaryVector, &
     LOAD, NodalAlpha, Element, n, m, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), LOAD,NodalAlpha

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n,  m

     REAL(KIND=dp) :: Basis(n), DetJ

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0_dp
     BoundaryMatrix = 0.0_dp
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     Force = LOAD
     Alpha = NodalAlpha / m

     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
       IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
         x = SUM( Nodes % x(1:n)*Basis )
         y = SUM( Nodes % y(1:n)*Basis )
         z = SUM( Nodes % z(1:n)*Basis )
         s = s * CoordinateSqrtMetric( x,y,z )
       END IF
!------------------------------------------------------------------------------
!      Force = SUM( LOAD(1:n) * Basis )
!      Alpha = SUM( NodalAlpha(1:n) * Basis )

       DO p=1,N
         DO q=1,M
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + s * Alpha * Basis(p)
         END DO
       END DO

       DO p=1,N
         BoundaryVector(p) = BoundaryVector(p) + s * Force * Basis(p)
       END DO
     END DO
   END SUBROUTINE IntegOverA
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindGapIndexes( Element, Indexes, n )
!------------------------------------------------------------------------------
      TYPE(Element_t) :: Element
      INTEGER :: n,Indexes(:)
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      Left  => Element % BoundaryInfo % Left
      Right => Element % BoundaryInfo % Right

      IF ( .NOT.ASSOCIATED(Left) .OR. .NOT.ASSOCIATED(Right) ) RETURN

      l = 0
      DO i=1,n
        Parent => Left
        k = Element % NodeIndexes(i)

        IF ( ANY( Parent % NodeIndexes == k ) ) &
          Parent => Right

        x0 = ElementNodes % x(i)
        y0 = ElementNodes % y(i)
        z0 = ElementNodes % z(i)
        DO j=1,Parent % TYPE % NumberOfNodes
          k = Parent % NodeIndexes(j)
          x = Solver % Mesh % Nodes % x(k) - x0
          y = Solver % Mesh % Nodes % y(k) - y0
          z = Solver % Mesh % Nodes % z(k) - z0
          IF ( x**2 + y**2 + z**2 < AEPS ) EXIT
        END DO
        Indexes(i) = k
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE FindGapIndexes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE AddHeatGap( Solver, Element, STIFF, TempPerm )
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: TempPerm(:)
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l, Ind(n)
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      CALL FindGapIndexes( Element, Ind, n )
      DO i=1,n
        DO j=1,n
          k = TempPerm( Element % NodeIndexes(i) )
          l = TempPerm( Ind(j) )
          IF ( k > 0 .AND. l > 0 ) THEN
            CALL AddToMatrixElement( Solver % Matrix,k,l,-STIFF(i,j) )
          END IF
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatGap
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE HeatSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION HeatSolver_Boundary_Residual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Radiation

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,n,l,t,dim,Pn,En,nd
     LOGICAL :: stat, Found
     INTEGER, ALLOCATABLE :: Indexes(:)

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), ExtTemperature(:), &
       TransferCoeff(:), EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
       dBasisdx(:,:), Temperature(:), Flux(:), NodalEmissivity(:)

     REAL(KIND=dp) :: Conductivity, Emissivity, StefanBoltzmann

     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, gx, gy, gz

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE., Dirichlet
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Indicator = 0.0d0
     Gnorm     = 0.0d0

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT
!    
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     en = Edge % TYPE % NumberOfNodes
     pn = Element % TYPE % NumberOfNodes

     ALLOCATE(EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     nd = GetElementNOFDOFs(Element)
     ALLOCATE( Temperature(nd), Basis(nd), ExtTemperature(nd), &
        TransferCoeff(en), x(en), y(en), z(en), EdgeBasis(nd), &
        dBasisdx(nd,3), NodalConductivity(nd), Flux(nd), &
        NodalEmissivity(nd), Indexes(nd) ) 

     nd = GetElementDOFs(Indexes,Element)

     ALLOCATE(Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )
     Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))


     DO l = 1,en
       DO k = 1,pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO
!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------

     Indicator    = 0.0d0
     EdgeLength   = 0.0d0
     ResidualNorm = 0.0d0

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                 'Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        Dirichlet = ListCheckPresent( Model % BCs(j) % Values,'Temperature')       
        
!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:en) = ListGetReal( Model % BCs(j) % Values, &
          'Heat Flux', en, Edge % NodeIndexes, Found )

!       ...convective heat transfer:
!       ----------------------------
        TransferCoeff(1:en) =  ListGetReal( Model % BCs(j) % Values, &
          'Heat Transfer Coefficient', en, Edge % NodeIndexes, Found )

        ExtTemperature(1:en) = ListGetReal( Model % BCs(j) % Values, &
          'External Temperature', en, Edge % NodeIndexes, Found )

!       ...black body radiation:
!       ------------------------
        Emissivity      = 0.0d0
        StefanBoltzmann = 0.0d0

        SELECT CASE(ListGetString(Model % BCs(j) % Values,'Radiation',Found))
           !------------------
           CASE( 'idealized' )
           !------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                  ListGetConstReal( Model % Constants,'Stefan Boltzmann',UnfoundFatal=.TRUE. )

           !---------------------
           CASE( 'diffuse gray' )
           !---------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

              ExtTemperature(1:en) =  ComputeRadiationLoad( Model, &
                      Mesh, Edge, Quant, Perm, Emissivity )
        END SELECT

!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Heat Conductivity', Hwrk, en, Edge % NodeIndexes )
        NodalConductivity(1:en) = Hwrk(1,1,1:en)

!       elementwise nodal solution:
!       ---------------------------
        nd = GetElementDOFs(Indexes,Element)
        Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0

        IntegStuff = GaussPoints(Edge)

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
                    EdgeBasis, dBasisdx )
           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ
           ELSE
              gx = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
              gy = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
              gz = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, gx, gy, gz )

              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:en) * x(1:en) )
           v = SUM( EdgeBasis(1:en) * y(1:en) )
           w = SUM( EdgeBasis(1:en) * z(1:en) )
           stat = ElementInfo(Element,Nodes, u, v, w, detJ,Basis,dBasisdx )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Conductivity = SUM( NodalConductivity(1:en) * EdgeBasis(1:en) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:en) * EdgeBasis(1:en) )

!          convective ...:
!          ----------------
           Residual = Residual + SUM(TransferCoeff(1:en) * EdgeBasis(1:en)) * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) )

!          black body radiation...:
!          -------------------------
           Residual = Residual + &
                Emissivity * StefanBoltzmann * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) ** 4 - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) ** 4 )

!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,dim
                 Residual = Residual + Conductivity  * &
                    SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(k)

                 Gnorm = Gnorm + s * (Conductivity * &
                       SUM(dBasisdx(1:nd,k) * Temperature(1:nd)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,dim
                 DO l=1,dim
                    Residual = Residual + Metric(k,l) * Conductivity  * &
                       SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
                      SUM(dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l))**2
                 END DO
              END DO
           END IF

           EdgeLength   = EdgeLength + s
           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * Residual ** 2
           END IF
        END DO
        EXIT
     END DO

     IF ( CoordinateSystemDimension() == 3 ) EdgeLength = SQRT(EdgeLength)

!    Gnorm = EdgeLength * Gnorm
     Indicator = EdgeLength * ResidualNorm
!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Boundary_Residual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION HeatSolver_Edge_Residual(Model,Edge,Mesh,Quant,Perm) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,dim,En,Pn,nd
     INTEGER, ALLOCATABLE :: Indexes(:)
     LOGICAL :: stat, Found
     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), &
            EdgeBasis(:), Basis(:), dBasisdx(:,:), Temperature(:)

     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump, Conductivity

     REAL(KIND=dp) :: u, v, w, s, detJ

     REAL(KIND=dp) :: Residual, ResidualNorm, Area

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

     !    Initialize:
     !    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
     END DO

     Grad = 0.0d0
!
!    ---------------------------------------------

     n = Mesh % MaxElementDOFs
     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     en = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( NodalConductivity(en), EdgeBasis(en), Basis(n), &
        dBasisdx(n,3), x(en), y(en), z(en), Temperature(n), Indexes(n) )

!    Integrate square of jump over edge:
!    -----------------------------------
     ResidualNorm = 0.0d0
     EdgeLength   = 0.0d0
     Indicator    = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
           v = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
           w = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! 
        ! Compute flux over the edge as seen by elements
        ! on both sides of the edge:
        ! ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           pn = Element % TYPE % NumberOfNodes

           DO j = 1,en
              DO k = 1,pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM(EdgeBasis(1:en) * x(1:en))
           v = SUM(EdgeBasis(1:en) * y(1:en))
           w = SUM(EdgeBasis(1:en) * z(1:en))
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           nd = GetElementDOFs(Indexes,Element)
           Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
           Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
           Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))

           stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx)
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                   'Heat Conductivity', Hwrk,en, Edge % NodeIndexes )

           NodalConductivity(1:en) = Hwrk( 1,1,1:en )
           Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
!
!          Temperature at element nodal points:
!          ------------------------------------
           Temperature(1:nd) = Quant( Perm(Indexes(1:nd)) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,dim
              Grad(j,i) = Conductivity * SUM( dBasisdx(1:nd,j) * Temperature(1:nd) )
           END DO
        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0d0
        DO k=1,dim
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,dim
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
     END DO

     IF (dim==3) EdgeLength = SQRT(EdgeLength)
     Indicator = EdgeLength * ResidualNorm

!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Edge_Residual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION HeatSolver_Inside_Residual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,n,t,dim,nd
     INTEGER, ALLOCATABLE :: Indexes(:)

     LOGICAL :: stat, Found, Compressible, VolSource
     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalCapacity(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
     REAL(KIND=dp), ALLOCATABLE :: Velo(:,:), Pressure(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Temperature(:), PrevTemp(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:)

     REAL(KIND=dp) :: u, v, w, s, detJ, Density, Capacity

     REAL(KIND=dp) :: SpecificHeatRatio, ReferencePressure, dt
     REAL(KIND=dp) :: Source, Residual, ResidualNorm, Area, Conductivity

     TYPE( ValueList_t ), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Fnorm     = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
     IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN

     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

!    Alllocate local arrays
!    ----------------------
     nd = GetElementNOFDOFs(Element)
     n = GetElementNOFNodes(Element)
     ALLOCATE( NodalDensity(nd), NodalCapacity(nd), NodalConductivity(nd),      &
         Velo(3,nd), Pressure(nd), NodalSource(nd), Temperature(nd), PrevTemp(nd), &
         Basis(nd), dBasisdx(nd,3), ddBasisddx(nd,3,3), Indexes(nd) )
!
!    Element nodal points:
!    ---------------------
     ALLOCATE( Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )

     nd = GetElementDOFs(Indexes,Element)
     Nodes % x = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z = Mesh % Nodes % z(Indexes(1:nd))
!
!    Elementwise nodal solution:
!    ---------------------------
     Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))
!
!    Check for time dep.
!    -------------------
     PrevTemp(1:nd) = Temperature(1:nd)
     dt = Model % Solver % dt
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Temperature', .TRUE. )
        PrevTemp(1:nd) = Var % PrevValues(Var % Perm(Indexes(1:nd)),1)
     END IF
!
!    Material parameters: conductivity, heat capacity and density
!    -------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL ListGetRealArray(Material,'Heat Conductivity',Hwrk,n,Element % NodeIndexes)
     NodalConductivity(1:n) = Hwrk( 1,1,1:n)
     NodalDensity(1:n) = GetReal(Material, 'Density', Found )
     NodalCapacity(1:n) = GetReal(Material, 'Heat Capacity', Found )
!
!    Check for compressible flow equations:
!    --------------------------------------
     Compressible = .FALSE.

     IF (  ListGetString( Material, 'Compressibility Model', Found ) == &
                 'perfect gas equation 1' ) THEN

        Compressible = .TRUE.

        Pressure = 0.0d0
        Var => VariableGet( Mesh % Variables, 'Pressure', .TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Pressure(1:n) = Var % Values(Var % Perm(Indexes(1:n)))
        END IF

        ReferencePressure = GetConstReal( Material, 'Reference Pressure' )
        SpecificHeatRatio = GetConstReal( Material, 'Specific Heat Ratio' )
        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * NodalCapacity(1:n) * Temperature(1:n) )
     END IF
!
!    Get (possible) convection velocity at the nodes of the element:
!    ----------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                minv=1, maxv=Model % NumberOFEquations )

     Velo = 0.0d0
     SELECT CASE( ListGetString( Model % Equations(k) % Values, &
                         'Convection', Found ) )

        !-----------------
        CASE( 'constant' )
        !-----------------

           Velo(1,1:n) = GetReal( Material, 'Convection Velocity 1',Found )
           Velo(2,1:n) = GetReal( Material, 'Convection Velocity 2',Found )
           Velo(3,1:n) = GetReal( Material, 'Convection Velocity 3',Found )

        !-----------------
        CASE( 'computed' )
        !-----------------

           Var => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              IF ( ALL(Var % Perm(Indexes(1:n)) > 0 ) ) THEN
                 Velo(1,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
                 Var => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(2,1:n) = Var % Values(Var % Perm(Indexes(1:n)) )
                 Var => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(3,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
              END IF
           END IF

     END SELECT

!
!    Heat source:
!    ------------
!
     k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',Found, &
                 1, Model % NumberOFBodyForces)

     NodalSource = 0.0d0
     IF( k > 0 ) THEN
       NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
           'Volumetric Heat Source',VolSource ) 
       IF( .NOT. VolSource ) THEN
         NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
             'Heat Source',  Found )
       END IF
     END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0d0
     Area = 0.0d0

     IntegStuff = GaussPoints( Element )
     ddBasisddx = 0

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM(Basis(1:nd) * Nodes % x(1:nd))
           v = SUM(Basis(1:nd) * Nodes % y(1:nd))
           w = SUM(Basis(1:nd) * Nodes % z(1:nd))

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Capacity     = SUM(NodalCapacity(1:n) * Basis(1:n))
        Density      = SUM(NodalDensity(1:n) * Basis(1:n))
        Conductivity = SUM(NodalConductivity(1:n) * Basis(1:n))
!
!       Residual of the convection-diffusion (heat) equation:
!        R = \rho * c_p * (@T/@t + u.grad(T)) - &
!            div(C grad(T)) + p div(u) - h,
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = \rho * c_p * (@T/@t + u^j T_{,j}) - &
!          g^{jk} (C T_{,j}}_{,k} + p div(u) - h
!       ---------------------------------------------------
!

        IF( VolSource ) THEN
          Residual = -SUM( NodalSource(1:n) * Basis(1:n) )
        ELSE
          Residual = -Density * SUM( NodalSource(1:n) * Basis(1:n) )
        END IF
          
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,dim
!
!             - grad(C).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                 SUM( NodalConductivity(1:n) * dBasisdx(1:n,j) )

!
!             - C div(grad(T)):
!             -------------------
!
              Residual = Residual - Conductivity * &
                 SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,j) )
           END DO
        ELSE
           DO j=1,dim
              DO k=1,dim
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Conductivity * &
                    SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,dim
                    Residual = Residual + Metric(j,k) * Conductivity * &
                      Symb(j,k,l) * SUM( Temperature(1:nd) * dBasisdx(1:nd,l) )
                 END DO
              END DO
           END DO
        END IF

!       + \rho * c_p * (@T/@t + u.grad(T)):
!       -----------------------------------
        Residual = Residual + Density * Capacity *  &
           SUM((Temperature(1:nd)-PrevTemp(1:nd))*Basis(1:nd)) / dt

        DO j=1,dim
           Residual = Residual + &
              Density * Capacity * SUM(Velo(j,1:n) * Basis(1:n)) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) )
        END DO


        IF ( Compressible ) THEN
!
!          + p div(u) or p u^j_{,j}:
!          -------------------------
!
           DO j=1,dim
              Residual = Residual + &
                 SUM( Pressure(1:n) * Basis(1:n) ) * &
                      SUM( Velo(j,1:n) * dBasisdx(1:n,j) )

              IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO k=1,dim
                    Residual = Residual + &
                       SUM( Pressure(1:n) * Basis(1:n) ) * &
                           Symb(j,k,j) * SUM( Velo(k,1:n) * Basis(1:n) )
                 END DO
              END IF
           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
        DO i=1,dim
           Fnorm = Fnorm + s * (Density *SUM(NodalSource(1:n)*Basis(1:n)))**2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
     END DO

!    Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm
!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Inside_Residual
!------------------------------------------------------------------------------
