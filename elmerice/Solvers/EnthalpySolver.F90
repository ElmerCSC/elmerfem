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
!
!/******************************************************************************
! *
! *  Module containing a solver for enthalpy equation in glaciology
! *
! ******************************************************************************
! *
! *  Authors: Adrien Gilbert, 
! *
! *  Original Date: Jun 2014
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
   SUBROUTINE EnthalpySolver( Model,Solver,Timestep,TransientSimulation )
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
         Found, GotIt, HeatFluxBC, HeatGapBC, GotMeltPoint, IsRadiation, InfBC, UnFoundFatal=.TRUE.
! Which compressibility model is used
     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, ConvectionField
     INTEGER :: CompressibilityModel

     LOGICAL :: AllocationsDone = .FALSE.,PhaseSpatial=.FALSE., &
        PhaseChange=.FALSE., CheckLatentHeatRelease=.FALSE., FirstTime, &
        SmartHeaterControl, IntegralHeaterControl, HeaterControlLocal, SmartTolReached=.FALSE., &
        TransientHeaterControl, SmartHeaterAverage, ConstantBulk, SaveBulk, &
	TransientAssembly
     LOGICAL, POINTER :: SmartHeaters(:), IntegralHeaters(:)

     TYPE(Variable_t), POINTER :: TempSol,FlowSol,HeatSol,CurrentSol, MeshSol, DensitySol
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER, POINTER :: EnthalpyPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:)

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter,MDOFs, &
         SmartHeaterBC, SmartHeaterNode, DoneTime=0
     REAL(KIND=dp) :: NonlinearTol,NewtonTol,SmartTol,Relax, &
            SaveRelax,dt,dt0,CumulativeTime, VisibleFraction, PowerScaling=1.0, PrevPowerScaling=1.0, &
            PowerRelax, PowerTimeScale, PowerSensitivity, xave, yave, Normal(3), &
	    dist, mindist, ControlPoint(3)

     REAL(KIND=dp), POINTER :: Enthalpy_h(:),PrevEnthalpy_h(:),FlowSolution(:), &
       ElectricCurrent(:), PhaseChangeIntervals(:,:),ForceVector(:), &
       PrevSolution(:), HC(:), Hwrk(:,:,:),MeshVelocity(:), XX(:), YY(:),ForceHeater(:),&
       RealWork(:,:)

     REAL(KIND=dp), ALLOCATABLE :: vals(:)
     REAL(KIND=dp) :: Jx,Jy,Jz,JAbs, Power, MeltPoint, IntHeatSource

     REAL(KIND=dp) :: Tref,Tm,P,hm,hi,L_heat,A_cap,B_cap,Ptriple,Psurf,beta

     CHARACTER(LEN=MAX_NAME_LEN) :: PressureName,WaterName,PhaseEnthName,TempName

     TYPE(Variable_t), POINTER :: PressureVariable,PhaseChangeEnthVar,WaterVar,TemphomoVar
     REAL(KIND=dp), POINTER :: PressureValues(:),PhaseChangeEnthValues(:)
     INTEGER, POINTER :: PressurePerm(:),PhaseChangeEnthPerm(:)

      INTEGER, ALLOCATABLE, SAVE :: Indexes(:), SaveIndexes(:)

     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
       FORCE(:), U(:), V(:), W(:), MU(:,:),TimeForce(:), &
       Density(:), LatentHeat(:), HeatTransferCoeff(:), &
       HeatCapacity(:),WaterDiffusivity(:), Enthalpy(:), Viscosity(:), LocalEnthalpy_h(:), &
       NodalEmissivity(:), ElectricConductivity(:), Permeability(:), Work(:), C0(:), &
       Pressure(:), dPressuredt(:), GasConstant(:),AText(:), HeaterArea(:), &
       HeaterTarget(:), HeaterScaling(:), HeaterDensity(:), HeaterSource(:), &
       HeatExpansionCoeff(:), ReferenceEnthalpy_h(:), PressureCoeff(:), &
       PhaseVelocity(:,:), HeatConductivityIso(:), &
       PerfusionRate(:), PerfusionDensity(:), PerfusionHeatCapacity(:), PerfusionRefEnthalpy_h(:)

     SAVE U, V, W, MU, MASS, STIFF, LOAD, PressureCoeff, &
       FORCE, ElementNodes, HeatConductivity, HeatCapacity,WaterDiffusivity, HeatTransferCoeff, &
       Enthalpy, Density, LatentHeat, PhaseVelocity, AllocationsDone, Viscosity, TimeForce, &
       LocalNodes, LocalEnthalpy_h, Work, ElectricConductivity, &
       NodalEmissivity, Permeability, C0, dPressuredt, Pressure, &
       GasConstant,AText,Hwrk, XX, YY, ForceHeater, Power, HeaterArea, HeaterTarget, &
       HeaterScaling, HeaterDensity, HeaterSource, SmartHeaters, IntegralHeaters, SmartTolReached,    &
       ReferenceEnthalpy_h, HeatExpansionCoeff, PrevPowerScaling, PowerScaling, &
       MeltPoint, DoneTime, SmartHeaterNode, SmartHeaterBC, SmartHeaterAverage, &
       HeatConductivityIso, &
       PerfusionRate, PerfusionDensity, PerfusionHeatCapacity, PerfusionRefEnthalpy_h


     INTERFACE
        FUNCTION HeatBoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatBoundaryResidual

        FUNCTION HeatEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION HeatEdgeResidual

        FUNCTION HeatInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION HeatInsideResidual
     END INTERFACE

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1
#else
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
#endif

!------------------------------------------------------------------------------
!    The View and Gebhardt factors may change. If this is necessary, this is 
!    done within this subroutine. The routine is called in the
!    start as it may affect the matrix topology.
!    Newton lineariarization option is needed only when there is radiation.
!------------------------------------------------------------------------------
     IsRadiation = ListCheckPresentAnyBC( Model,'Radiation')
     IF( IsRadiation ) THEN
       CALL RadiationFactors( Solver, .FALSE.)
     END IF

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     StiffMatrix => GetMatrix()
     ForceVector => Solver % Matrix % RHS

     TempSol => Solver % Variable
     EnthalpyPerm    => TempSol % Perm
     Enthalpy_h => TempSol % Values
     VarName = GetVarName( TempSol ) 
  
     LocalNodes = COUNT( EnthalpyPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     SolverParams => GetSolverParams()
     ConvectionField = GetString( SolverParams, 'Enthalpy_h Convection Field', Found )

     IF ( Found ) THEN
       FlowSol => VariableGet( Solver % Mesh % Variables, ConvectionField )
     ELSE
       FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
     END IF

     IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm     => FlowSol % Perm
       NSDOFs       =  FlowSol % DOFs
       FlowSolution => FlowSol % Values
     END IF

     DensitySol => VariableGet( Solver % Mesh % Variables, 'Enthalpy Density' )

     ! Check whether we have some heater controls. This will affect initialization stuff. 
     SmartHeaterControl = ListCheckPresentAnyBodyForce( Model,'Smart Heater Control')
     IntegralHeaterControl = ListCheckPresentAnyBodyForce( Model,'Integral Heat Source')
   
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
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
                 ReferenceEnthalpy_h,  &
                 MASS,       &
                 LocalEnthalpy_h,      &
                 HeatCapacity,Enthalpy,WaterDiffusivity, &
                 NodalEmissivity,       &
                 GasConstant, AText,    &
                 HeatConductivity,      &
                 STIFF,LOAD,            &
                 Indexes, SaveIndexes,  &
                 FORCE, TimeForce,      &
                 HeatConductivityIso,   &
                 PerfusionRate,         &
                 PerfusionDensity,      &
                 PerfusionHeatCapacity, &
                 PerfusionRefEnthalpy_h )
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
                 ReferenceEnthalpy_h( N ),            &
                 MASS(  2*N,2*N ),                     &
                 LocalEnthalpy_h( N ),                &
                 HeatCapacity( N ),Enthalpy( N ), WaterDiffusivity( N ),     &
                 NodalEmissivity( N ),                 &
                 GasConstant( N ),AText( N ),          &
                 HeatConductivity( 3,3,N ),            &
                 HeatConductivityIso( N ),             &
                 STIFF( 2*N,2*N ),LOAD( N ), &
                 FORCE( 2*N ), TimeForce(2*N), &
                 PerfusionRate( N ),     &
                 PerfusionDensity( N ),      &
                 PerfusionHeatCapacity( N ), &
                 PerfusionRefEnthalpy_h( N ), &
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
          n = SIZE( Enthalpy_h )
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
                     'Stefan Boltzmann' )
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

     Relax = GetCReal( SolverParams, &
               'Nonlinear System Relaxation Factor',Found )
     IF ( .NOT.Found ) Relax = 1

     TransientAssembly = TransientSimulation
     dt0 = ListGetConstReal(SolverParams,'Steady State Transition Timestep',Found)
     IF(.NOT. Found) dt0 = ListGetConstReal(SolverParams,'Smart Heater Time Scale',Found)

     IF(Found .AND. dt > dt0) TransientAssembly = .FALSE.


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
           'Smart Heater Enthalpy_h',GotMeltPoint)           
              
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
             IF( EnthalpyPerm(l) == 0 ) CYCLE
             
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
             'Smart Heater Enthalpy_h',Found)
         IF(.NOT. Found) THEN
           DO k=1, Model % NumberOfMaterials
             MeltPoint = GetCReal( Model % Materials(k) % Values, &
                 'Melting Point', Found )
             IF(Found) EXIT
           END DO
           IF(.NOT. Found) THEN
             CALL Fatal('HeatSolver','Smart Heater Enthalpy_h / Melting Point is undefined')
           END IF
         END IF
         
         ! Find the node related to Enthalpy_h control 
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
     PrevSolution = Enthalpy_h(1:LocalNodes)
     IF ( TransientSimulation ) THEN
       PrevEnthalpy_h => Solver % Variable % PrevValues(:,1)
     END IF
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     Norm = Solver % Variable % Norm


!Calculate Phase change enthalpy ==================================================================
!==================================================================================================

	  A_cap = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >Enthalpy Heat Capacity A< &
         & defined in model constants. Using >7.253< as default (Paterson, 1994).')
         A_cap = 7.253
      END IF
	  
	  B_cap = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >Enthalpy Heat Capacity B< & 
         & defined in model constants. Using >146.3< as default (Paterson, 1994).')
         B_cap = 146.3
      END IF

  PressureName = GetString(Constants,'Pressure Variable', GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >Pressure Variable< defined. Using >Pressure< as default.')
         WRITE(PressureName,'(A)') 'Pressure'
      ELSE
         WRITE(Message,'(a,a)') 'Variable Name for pressure: ', PressureName
         CALL INFO('EnthalpySolve',Message,Level=12)
      END IF

      Tref = GetConstReal(Model % Constants, "T_ref_enthalpy",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >Reference Enthalpy< defined. Using >200K< as default.')
         Tref = 200.0
      END IF
	  
	  beta = GetConstReal(Model % Constants, "beta_clapeyron",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >beta_clapeyron< defined. Using >9.74 10-2 K Mpa-1< as default.')
         beta = 0.0974
      END IF
	  
	  Psurf = GetConstReal(Model % Constants, "P_surf",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >Psurf< defined. Using >1.013 10-1 MPa < as default.')
         Psurf = 0.1013
      END IF
	  
	  Ptriple = GetConstReal(Model % Constants,"P_triple",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >P_triple< defined. Using >0.061173 MPa < as default.')
         Ptriple = 0.061173
      END IF

      PressureVariable => VariableGet(Solver % Mesh %Variables,PressureName)
      IF ( ASSOCIATED( PressureVariable ) ) THEN
       PressurePerm    => PressureVariable % Perm
       PressureValues => PressureVariable % Values
      END IF

 PhaseEnthName = GetString(SolverParams , 'Exported Variable 1', Found )

     IF (.NOT.Found) &
        CALL FATAL('EnthalpySolver','No value > Exported Variable 1 (enthalpy phase change) < found in Solver')

      WaterVar => VariableGet(Solver % Mesh % Variables, WaterName)
      TemphomoVar => VariableGet(Solver % Mesh % Variables,TempName)

      PhaseChangeEnthVar => &
               VariableGet(Solver % Mesh % Variables,PhaseEnthName)
      IF ( ASSOCIATED( PhaseChangeEnthVar ) ) THEN
         PhaseChangeEnthPerm => PhaseChangeEnthVar % Perm    
         PhaseChangeEnthValues => PhaseChangeEnthVar % Values  
      END IF

do i=1,model % NumberOfNodes
  P=PressureValues (PressurePerm(i))+Psurf
  Tm=273.16-beta*(P-Ptriple)! Clapeyron
  PhaseChangeEnthValues (PhaseChangeEnthPerm(i)) = A_cap*0.5*Tm*Tm+B_cap*Tm-Tref*(A_cap*0.5*Tref+B_cap) ! use CP(T)=B_cap+A_cap*T
enddo

!=====================================================================================================================
!=====================================================================================================================


     DO iter=1,NonlinearIter
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       WRITE( Message,* ) 'Enthalpy_h ITERATION', iter
       CALL Info( 'HeatSolve', Message, Level=4 )
       CALL Info( 'HeatSolve', '-------------------------------------',Level=4 )
       CALL Info( 'HeatSolve', ' ', Level=4 )
       CALL Info( 'HeatSolve', 'Starting Assembly...', Level=4 )

       IF ( ConstantBulk .AND. ASSOCIATED(Solver % Matrix % BulkValues) ) THEN
         Solver % Matrix % Values = Solver % Matrix % BulkValues
         Solver % Matrix % RHS = Solver % Matrix % BulkRHS
         GOTO 1000
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

             Density(1:n) = GetReal( Material, 'Enthalpy Density' )
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
       DO t=1,Solver % NumberOfActiveElements

         CALL AdvanceOutput(t,GetNOFActive())
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where Enthalpy_h 
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

         CALL GetScalarLocalSolution( LocalEnthalpy_h )
!------------------------------------------------------------------------------
!        Get element material parameters
!------------------------------------------------------------------------------
         HeatCapacity(1:n) = 1.0 
 	 WaterDiffusivity(1:n) = GetReal( Material, 'Enthalpy Water Diffusivity', Found )

IF (.NOT.Found) THEN
         CALL FATAL('EnthalpySolver', 'No value for >Enthalpy Water Diffusivity< defined in material.')
ENDIF
	

         CALL ListGetRealArray( Material,'Enthalpy Heat Diffusivity',Hwrk,n, &
                      Element % NodeIndexes )
					 
         HeatConductivity = 0.0d0
         IF ( SIZE(Hwrk,1) == 1 ) THEN
           DO i=1,3
		DO j=1,n
			IF (PhaseChangeEnthValues(PhaseChangeEnthPerm(Element % NodeIndexes(j)))<&
				&Enthalpy_h(EnthalpyPerm(Element % NodeIndexes(j)))) then
             			HeatConductivity( i,i,j ) = WaterDiffusivity(j)
			ELSE
				HeatConductivity( i,i,j ) = Hwrk( 1,1,j )
			ENDIF
		ENDDO
           END DO
         ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Hwrk,1))
		DO j=1,n
			IF (PhaseChangeEnthValues(PhaseChangeEnthPerm(Element % NodeIndexes(j)))<&
				&Enthalpy_h(EnthalpyPerm(Element % NodeIndexes(j)))) then
				HeatConductivity(i,i,j) = WaterDiffusivity(j)
			ELSE
			        HeatConductivity(i,i,j) = Hwrk(i,1,j)
			ENDIF
		ENDDO
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Hwrk,1))
             DO j=1,MIN(3,SIZE(Hwrk,2))
		DO k=1,n
			IF (PhaseChangeEnthValues(PhaseChangeEnthPerm(Element % NodeIndexes(k)))<&
				&Enthalpy_h(EnthalpyPerm(Element % NodeIndexes(k)))) then
				HeatConductivity( i,j,k ) = WaterDiffusivity(k)
			ELSE
				HeatConductivity( i,j,k ) = Hwrk(i,j,k)
			ENDIF

		ENDDO
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
           ReferenceEnthalpy_h(1:n) = GetReal( Material, 'Reference Enthalpy_h' )
           HeatExpansionCoeff(1:n) = GetReal( Material, 'Heat Expansion Coefficient' )

           Density(1:n) = GetReal( Material,'Enthalpy Density' )
           Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                ( LocalEnthalpy_h(1:n) - ReferenceEnthalpy_h(1:n) ) )

           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) &
             PressureCoeff(1:n) = LocalEnthalpy_h(1:n) * HeatExpansionCoeff(1:n) / &
                   ( 1-HeatExpansionCoeff(1:n)*( &
                               LocalEnthalpy_h(1:n)-ReferenceEnthalpy_h(1:n)) )
         ELSE IF ( CompressibilityModel == UserDefined1 ) THEN
           IF ( ASSOCIATED( DensitySol ) ) THEN
             CALL GetScalarLocalSolution( Density, 'Enthalpy Density' ) 
           ELSE
             Density(1:n) = GetReal( Material,'Enthalpy Density' )
           END IF
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
         ELSE
           PressureCoeff(1:n) = GetReal( Material, 'Pressure Coefficient', Found )
           IF ( .NOT. Found ) PressureCoeff(1:n) = 0.0_dp
           Density(1:n) = GetReal( Material, 'Enthalpy Density' )
         END IF

!------------------------------------------------------------------------------
! Take pressure deviation p_d as the dependent variable, p = p_0 + p_d
! for PerfectGas, read p_0
!------------------------------------------------------------------------------
         IF ( CompressibilityModel /= Incompressible ) THEN
           ReferencePressure = ListGetConstReal( Material, &
               'Reference Pressure', Found,UnFoundFatal=UnFoundFatal)
            !Previous default value: ReferencePressure = 0.0d0
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

         ELSE IF ( ConvectionFlag == 'computed' .AND. &
              ASSOCIATED(FlowSolution) ) THEN
           DO i=1,n
             k = FlowPerm(Element % NodeIndexes(i))
             IF ( k > 0 ) THEN
!------------------------------------------------------------------------------
               Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
               SELECT CASE( CompressibilityModel )
                 CASE( PerfectGas1 )
                   Density(i)  = Pressure(i) / &
                       ( GasConstant(i) * LocalEnthalpy_h(i) )
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
         ELSE IF (ConvectionFlag=='computed' ) THEN
           CALL Warn( 'HeatSolver', 'Convection model specified ' //  &
                 'but no associated flow field present?' )
         ELSE
           IF ( ALL(MU==0) ) C1 = 0.0D0 
         END IF
!------------------------------------------------------------------------------
!        Check if modelling Phase Change with Eulerian approach 
!------------------------------------------------------------------------------
         PhaseSpatial = .FALSE.
         IF (  PhaseChange ) THEN
           CALL EffectiveHeatCapacity()
         ELSE
           HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
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
           Load(1:n) = Density(1:n) *  GetReal( BodyForce, 'Heat Source', Found )
           
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
         PerfusionRate(1:n) = GetReal( BodyForce, 'Perfusion Rate', Found )
         
         IF ( Found ) THEN
           PerfusionRefEnthalpy_h(1:n) = GetReal( BodyForce, 'Perfusion Reference Enthalpy_h' )
           PerfusionDensity(1:n) = GetReal( BodyForce, 'Perfusion Density' )
           PerfusionHeatCapacity(1:n) = GetReal( BodyForce, 'Perfusion Heat Capacity' )
           C0(1:n) = PerfusionHeatCapacity(1:n) * PerfusionRate(1:n) * PerfusionDensity(1:n) 
           Load(1:n) = Load(1:n) + C0(1:n) * PerfusionRefEnthalpy_h(1:n)           
         END IF

!------------------------------------------------------------------------------
!        Get element local matrices, and RHS vectors
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalEnthalpy_h, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, UseBubbles, Element, n, ElementNodes )

!------------------------------------------------------------------------------
         ELSE
!------------------------------------------------------------------------------
           CALL DiffuseConvectiveGenCompose( &
               MASS, STIFF, FORCE, LOAD, &
               HeatCapacity, C0, C1*HeatCapacity(1:n), HeatConductivity, &
               PhaseSpatial, LocalEnthalpy_h, Enthalpy, U, V, W, &
               MU(1,1:n),MU(2,1:n),MU(3,1:n), Viscosity, Density, Pressure, &
               dPressureDt, PressureCoeff, CompressibilityModel /= Incompressible, &
               Stabilize, Element, n, ElementNodes )
!------------------------------------------------------------------------------
         END IF
!------------------------------------------------------------------------------

         IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN

           IF ( TransientAssembly .AND. .NOT. ConstantBulk ) THEN
             CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF

           CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
               ForceHeater, FORCE, n, 1, EnthalpyPerm(Element % NodeIndexes) )
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
      DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        n = GetElementNOFNodes()

        ! Check that the dimension of element is suitable for fluxes
        IF( .NOT. PossibleFluxElement(Element) ) CYCLE

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

        ! This checks whether there are any Dirichlet conditions on the 
        ! smart heater boundary. If there are the r.h.s. must be zero as 
        ! there can possibly not be any effect on Enthalpy_h.
        !-----------------------------------------------------------------
	IF ( HeaterControlLocal .AND. .NOT. TransientHeaterControl) THEN
          IF( ListCheckPresent(BC, Varname) ) THEN
             nd = GetElementDOFs(Indexes)
             ForceHeater(EnthalpyPerm(Indexes(1:nd))) = 0.0_dp
          END IF
        END IF

        HeatFluxBC = GetLogical( BC, 'Enthalpy Heat Flux BC', Found )
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

      CALL DefaultFinishAssembly()
      CALL Info( 'HeatSolve', 'Assembly done', Level=4 )

      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      at = CPUTime() - at
      st = CPUTime()

      PrevNorm = Norm

      IF(SmartHeaterControl .AND. NewtonLinearization .AND. SmartTolReached) THEN
      
        IF(.NOT. TransientHeaterControl) THEN
 
          CALL ListAddLogical(SolverParams, &
              'Skip Compute Nonlinear Change',.TRUE.)

          Relax = GetCReal( SolverParams, &
              'Nonlinear System Relaxation Factor', Found )
          
          IF ( Found .AND. Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values, &
                'Nonlinear System Relaxation Factor', 1.0d0 )
          ELSE
            Relax = 1.0d0
          END IF          

          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              ForceHeater, XX, Norm, 1, Solver )
         
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, YY, Norm, 1, Solver )

          CALL ListAddLogical(SolverParams,'Skip Compute Nonlinear Change',.FALSE.)
        ELSE          
          CALL SolveSystem( Solver % Matrix, ParMatrix, &
              Solver % Matrix % RHS, Enthalpy_h, Norm, 1, Solver )
          YY = Enthalpy_h
        END IF

        IF(.NOT. SmartHeaterAverage) THEN
          xave = XX(EnthalpyPerm(SmartHeaterNode))
          yave = YY(EnthalpyPerm(SmartHeaterNode))
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
              xave = xave + SUM( XX(EnthalpyPerm(Element % NodeIndexes)) )
              yave = yave + SUM( YY(EnthalpyPerm(Element % NodeIndexes)) )
            END IF
          END DO
          xave = xave / j
          yave = yave / j 
          CALL ListAddConstReal(Model % Simulation,'res: Smart Heater Enthalpy_h',yave)
        END IF

        IF(.NOT. TransientHeaterControl) THEN
          IF ( ASSOCIATED(Solver % Variable % NonlinValues) ) THEN
            Solver % Variable % NonlinValues = Enthalpy_h
          END IF

          PowerScaling = (MeltPoint - yave) / xave 
          Enthalpy_h = YY + PowerScaling * XX

          ! The change is computed separately for the controlled Enthalpy_h field
          !-----------------------------------------------------------------------
          CALL ComputeChange(Solver,.FALSE.,LocalNodes,Enthalpy_h)
          Norm = Solver % Variable % Norm

        END IF

        IF(dt > PowerTimeScale) THEN
          IF ( Relax /= 1.0d0 ) THEN
            CALL ListAddConstReal( Solver % Values,  &
                'Nonlinear System Relaxation Factor', Relax )
          END IF
        END IF
      ELSE
        Norm = DefaultSolve()
      END IF


      IF( SmartHeaterControl .OR. IntegralHeaterControl) THEN
         
         CALL ListAddConstReal(Model % Simulation,'res: Heater Power Scaling',PowerScaling)
         
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
            
            IF( SmartHeaters(i)) CALL ListAddConstReal(Model % Simulation,'res: Heater Power Density',&
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
          Enthalpy_h(1:LocalNodes) = PrevSolution
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
        IF ( .NOT.TransientSimulation ) PrevSolution=Enthalpy_h(1:LocalNodes)
      END IF
!------------------------------------------------------------------------------
     
      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm   : ',Norm
      CALL Info( 'HeatSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'HeatSolve', Message, Level=4 )

      IF ( RelativeChange < NewtonTol .OR. iter >= NewtonIter ) &
               NewtonLinearization = .TRUE.
      IF ( RelativeChange < NonlinearTol .AND. &
          (.NOT. SmartHeaterControl .OR. SmartTolReached)) EXIT

      IF(SmartHeaterControl) THEN
        IF ( RelativeChange < SmartTol ) THEN
          SmartTolReached = .TRUE.
          YY = Enthalpy_h
        END IF
      END IF
      
!------------------------------------------------------------------------------
    END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------
    IF(TransientHeaterControl) THEN
      PowerRelax = GetCReal(Solver % Values,'Smart Heater Relaxation Factor', GotIt)
      IF(.NOT. GotIt) PowerRelax = 1.0_dp
      PowerSensitivity = ListGetConstReal(Solver % Values,'Smart Heater Power Sensivity',GotIt,UnFoundFatal=UnFoundFatal)
      !Previous default value: PowerSensitivity = 4.0_dp
      PowerScaling = PowerScaling * (1 + PowerSensitivity * PowerRelax * (MeltPoint/yave - 1.0d0) ) 

      IF( ListGetLogical( Solver % Values,'Smart Heater Transient Speedup',GotIt ) ) THEN
        Enthalpy_h = Enthalpy_h * (1 + PowerRelax * (MeltPoint/yave - 1.0d0)   )     
      END IF
      YY = Enthalpy_h
    END IF


!------------------------------------------------------------------------------
!   Compute temperature homologous and water content
!------------------------------------------------------------------------------

 WaterName = GetString(SolverParams , 'Exported Variable 2', Found )

     IF (.NOT.Found) &
        CALL FATAL('EnthalpySolver','No value > Exported Variable 2 (water content) < found in Solver')

 TempName = GetString(SolverParams , 'Exported Variable 3', Found )

     IF (.NOT.Found) &
        CALL FATAL('EnthalpySolver','No value > Exported Variable 3 (temperature) < found in Solver')

      WaterVar => VariableGet(Solver % Mesh % Variables, WaterName)
      TemphomoVar => VariableGet(Solver % Mesh % Variables,TempName)

      L_heat = GetConstReal(Model % Constants, "L_heat",GotIt)
      IF (.NOT.GotIt) THEN
         CALL WARN('EnthalpySolver', 'No Keyword >L_heat< defined in model constants. Using >334000Jkg-1< as default.')
         L_heat = 334000.0
      END IF
	  
do i=1,Model % NumberOfNodes

hi = Enthalpy_h (EnthalpyPerm (i) )
hm = PhaseChangeEnthValues (PhaseChangeEnthPerm(i))

  if (hi<hm) then
    WaterVar % values ( WaterVar % perm (i) ) = 0.0
    TemphomoVar % values ( TemphomoVar % perm (i) ) = &
    & (-B_cap+(B_cap**2+A_cap*2*(A_cap*0.5*Tref**2+B_cap*Tref+hi))**0.5 ) / A_cap - 273.16 ! Use CP(T)=B_cap+A_cap*T
  else
    TemphomoVar % values ( TemphomoVar % perm (i) ) =  &
    & (-B_cap+(B_cap**2+A_cap*2*(A_cap*0.5*Tref**2+B_cap*Tref+hm))**0.5 ) / A_cap - 273.16 ! Use CP(T)=B_cap+A_cap*T
    WaterVar % values ( WaterVar % perm (i) ) = (hi-hm)/L_heat
  endif



enddo

!------------------------------------------------------------------------------
!   Compute cumulative time done by now and time remaining
!------------------------------------------------------------------------------
    IF ( .NOT. TransientSimulation ) EXIT
    CumulativeTime = CumulativeTime + dt
    dt = Timestep - CumulativeTime

   END DO ! time interval
   Solver % dt = Timestep

!------------------------------------------------------------------------------
   CALL  ListAddConstReal( Solver % Values,  &
        'Nonlinear System Relaxation Factor', SaveRelax )
!------------------------------------------------------------------------------

   DEALLOCATE( PrevSolution )

   IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) &
      CALL RefineMesh( Model,Solver,Enthalpy_h,EnthalpyPerm, &
            HeatInsideResidual, HeatEdgeResidual, HeatBoundaryResidual )

CONTAINS



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

      IF ( Found .AND. RadiationFlag(1:4) /= 'none' ) THEN
        NodalEmissivity(1:n) = GetReal(BC, 'Emissivity', Found)
        IF(.NOT. Found) THEN
           NodalEmissivity(1:n) = GetParentMatProp( 'Emissivity' )
        END IF
        Emissivity = SUM( NodalEmissivity(1:n) ) / n

!------------------------------------------------------------------------------
        IF (  RadiationFlag(1:9) == 'idealized' ) THEN
          AText(1:n) = GetReal( BC, 'Radiation External Enthalpy_h',Found )
          IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Enthalpy_h' )
        ELSE
          CALL DiffuseGrayRadiation( Model, Solver, Element, & 
              Enthalpy_h, EnthalpyPerm, ForceVector, VisibleFraction, Text)

          IF( GetLogical( BC, 'Radiation Boundary Open', Found) ) THEN
            AText(1:n) = GetReal( BC, 'Radiation External Enthalpy_h',Found )
            IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Enthalpy_h' )
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
!------------------------------------------------------------------------------
!       Add our own contribution to surface Enthalpy_h (and external
!       if using linear type iteration or idealized radiation)
!------------------------------------------------------------------------------
        DO j=1,n
          k = EnthalpyPerm(Element % NodeIndexes(j))
          Text = AText(j)

          IF ( .NOT. HeatGapBC .AND. NewtonLinearization ) THEN
             HeatTransferCoeff(j) = Emissivity * 4*Enthalpy_h(k)**3 * &
                               StefanBoltzmann
             LOAD(j) = Emissivity*(3*Enthalpy_h(k)**4+Text**4) * &
                               StefanBoltzmann
          ELSE
             HeatTransferCoeff(j) = Emissivity * (Enthalpy_h(k)**3 + &
             Enthalpy_h(k)**2*Text+Enthalpy_h(k)*Text**2 + Text**3) * &
                               StefanBoltzmann 
             LOAD(j) = HeatTransferCoeff(j) * Text
          END IF
        END DO
      END IF  ! of radition
!------------------------------------------------------------------------------

      Work(1:n)  = GetReal( BC, 'Heat Transfer Coefficient',Found )
      IF ( Found ) THEN
       AText(1:n) = GetReal( BC, 'External Enthalpy_h',Found )
       DO j=1,n
!------------------------------------------------------------------------------
!         BC: -k@T/@n = \alpha(T - Text)
!------------------------------------------------------------------------------
          k = EnthalpyPerm(Element % NodeIndexes(j))
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
           Density(1:n) = GetReal( Model % Materials(k) % Values, 'Enthalpy Density' )
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
      LOAD(1:n) = LOAD(1:n) +  GetReal( BC, 'Enthalpy Heat Flux', Found )
      WaterDiffusivity(1:n) = GetReal( Material, 'Enthalpy Water Diffusivity', Found )

      InfBC = ListGetLogical( BC,'Infinity BC '//TRIM(VarName),GotIt)
      IF( InfBC ) THEN
        AText(1:n) = GetReal( BC,'Infinity BC '//TRIM(VarName)//' Offset',GotIt)
        ! currently only isotropic heat conductivity supported
        HeatConductivityIso(1:n) = GetParentMatProp('Enthalpy Heat Diffusivity',Element,GotIt)

		DO k=1,n
			IF (PhaseChangeEnthValues(PhaseChangeEnthPerm(Element % NodeIndexes(k)))<&
				&Enthalpy_h(EnthalpyPerm(Element % NodeIndexes(k)))) then
				HeatConductivityIso(k) = WaterDiffusivity(k)
			ENDIF
		ENDDO


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
        CALL AddHeatGap( Solver, Element, STIFF, EnthalpyPerm)

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
      Enthalpy_h, EnthalpyPerm, ForceVector,AngleFraction, Text)
!------------------------------------------------------------------------------
      TYPE(Model_t)  :: Model
      TYPE(Solver_t) :: Solver
      TYPE(Element_t), POINTER :: Element
      INTEGER :: EnthalpyPerm(:)
      REAL(KIND=dp) :: Enthalpy_h(:), ForceVector(:)
      REAL(KIND=dp) :: AngleFraction, Text
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Area, Asum
      INTEGER :: i,j,k,l,m,ImplicitFactors
      INTEGER, POINTER :: ElementList(:)
!------------------------------------------------------------------------------
!     If linear iteration compute radiation load
!------------------------------------------------------------------------------


      Asum = 0.0d0
      IF ( .NOT. NewtonLinearization ) THEN
        Text = ComputeRadiationLoad( Model, Solver % Mesh, Element, &
                 Enthalpy_h, EnthalpyPerm, Emissivity, AngleFraction)
      ELSE   !  Full Newton-Raphson solver
!------------------------------------------------------------------------------
!       Go trough surfaces (j) this surface (i) is getting
!       radiated from.
!------------------------------------------------------------------------------

        Area  = ElementArea( Solver % Mesh, Element, n )
        ElementList => Element % BoundaryInfo % GebhardtFactors % Elements

        DO j=1,Element % BoundaryInfo % GebhardtFactors % NumberOfFactors

          RadiationElement => Solver % Mesh % Elements( ElementList(j) )

          Text = ComputeRadiationCoeff(Model,Solver % Mesh,Element,j) / ( Area )
          Asum = Asum + Text
!------------------------------------------------------------------------------
!         Gebhardt factors are given elementwise at the center
!         of the element, so take avarage of nodal Enthalpy_hs
!         (or integrate over surface j)
!------------------------------------------------------------------------------

          k = RadiationElement % TYPE % NumberOfNodes
          ImplicitFactors = Element % BoundaryInfo % GebhardtFactors % NumberOfImplicitFactors
          IF(ImplicitFactors == 0) &
              ImplicitFactors = Element % BoundaryInfo % GebhardtFactors % NumberOfFactors

          IF(j <= ImplicitFactors) THEN
            
            S = (SUM( Enthalpy_h( EnthalpyPerm( RadiationElement % &
                NodeIndexes))**4 )/k )**(1.0d0/4.0d0)
!------------------------------------------------------------------------------
!         Linearization of the G_jiT^4_j term
!------------------------------------------------------------------------------
            HeatTransferCoeff(1:n) = -4 * Text * S**3 * StefanBoltzmann
            LOAD(1:n) = -3 * Text * S**4 * StefanBoltzmann
!------------------------------------------------------------------------------
!         Integrate the contribution of surface j over surface i
!         and add to global matrix
!------------------------------------------------------------------------------
            CALL IntegOverA( STIFF, FORCE, LOAD, &
                HeatTransferCoeff, Element, n, k, ElementNodes ) 
            
            IF ( TransientAssembly ) THEN
              MASS = 0.d0
              CALL Add1stOrderTime( MASS, STIFF, &
                  FORCE,dt,n,1,EnthalpyPerm(Element % NodeIndexes),Solver )
            END IF
            
            DO m=1,n
              k1 = EnthalpyPerm( Element % NodeIndexes(m) )
              DO l=1,k
                k2 = EnthalpyPerm( RadiationElement % NodeIndexes(l) )
                CALL AddToMatrixElement( StiffMatrix,k1, &
                    k2,STIFF(m,l) )
              END DO
              ForceVector(k1) = ForceVector(k1) + FORCE(m)
            END DO

          ELSE

            S = (SUM( Enthalpy_h( EnthalpyPerm( RadiationElement % &
                NodeIndexes))**4 )/k )
            
            HeatTransferCoeff(1:n) = 0.0d0
            LOAD(1:n) = Text * S * StefanBoltzmann
            
            CALL IntegOverA( STIFF, FORCE, LOAD, &
                HeatTransferCoeff, Element, n, k, ElementNodes ) 
            
            DO m=1,n
              k1 = EnthalpyPerm( Element % NodeIndexes(m) )
              ForceVector(k1) = ForceVector(k1) + FORCE(m)
            END DO
            
          END IF 

        END DO

!------------------------------------------------------------------------------
!       We have already added all external Enthalpy_h contributions
!       to the matrix for the Newton type iteration
!------------------------------------------------------------------------------
        AngleFraction = Asum / Emissivity
        Text = 0.0

      END IF  !  of newton-raphson

    END SUBROUTINE DiffuseGrayRadiation
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE EffectiveHeatCapacity()
      LOGICAL :: Found, Specific
      REAL(KIND=dp), ALLOCATABLE :: dT(:)

!------------------------------------------------------------------------------
!     See if Enthalpy_h gradient indside the element is large enough 
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
! Check if local variation of Enthalpy_h is large enough to actually use the
! Spatial 2 model. Should perhaps be scaled to element size (or actually
! compute the gradient, but this will do for now...).
!------------------------------------------------------------------------------
          s = MAXVAL(LocalEnthalpy_h(1:n))-MINVAL(LocalEnthalpy_h(1:n))
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
            dT(1:n) = Enthalpy_h(EnthalpyPerm(Element % NodeIndexes)) - &
                     PrevEnthalpy_h(EnthalpyPerm(Element % NodeIndexes))

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

!-----------------------------------------------------------------------------
      SELECT CASE( PhaseChangeModel )

!------------------------------------------------------------------------------
! This phase change model is available only for some type of real entries 
! that have an implemented analytical derivation rule.
!-----------------------------------------------------------------------------
      CASE( PHASE_SPATIAL_1 )
        HeatCapacity(1:n) = ListGetReal( Material, &
             'Effective Heat Capacity', n,Element % NodeIndexes, Found )
        IF ( .NOT. Found ) THEN
          IF( Specific ) THEN
            HeatCapacity(1:n) = ListGetDerivValue( Material, &
                'Specific Enthalpy', n,Element % NodeIndexes )
            HeatCapacity(1:n) = Density(1:n) * HeatCapacity(1:n)
          ELSE
            HeatCapacity(1:n) = ListGetDerivValue( Material, &
                'Enthalpy', n,Element % NodeIndexes )            
          END IF
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
          
!------------------------------------------------------------------------------
      CASE( PHASE_TEMPORAL )
        ! When retrieving the value of enthalphy on the previous timestep 
        ! the relevant entries of the Enthalpy_h solution in the global vector
        ! are tampered in order to make the ListGetReal command work as wanted. 
        ! 1) Values at current Enthalpy_h     
        !------------------------------------------------------------------------
        IF( Specific ) THEN
          Work(1:n) = ListGetReal( Material,'Specific Enthalpy',n,Element % NodeIndexes )
        ELSE
          Work(1:n) = ListGetReal( Material,'Enthalpy',n,Element % NodeIndexes )
        END IF

        ! 2) Values at previous Enthalpy_h
        Enthalpy_h(EnthalpyPerm(Element % NodeIndexes)) = & 
            PrevEnthalpy_h(EnthalpyPerm(Element % NodeIndexes)) 

        IF( Specific ) THEN
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Specific Enthalpy', &
              n,Element % NodeIndexes )          
          HeatCapacity(1:n) = Density(1:n) * Work(1:n) / dT(1:n)
       ELSE
          Work(1:n) = Work(1:n) - ListGetReal( Material,'Enthalpy', &
              n,Element % NodeIndexes )
          HeatCapacity(1:n) = Work(1:n) / dT(1:n)
        END IF

        ! Revert to current Enthalpy_h
        Enthalpy_h(EnthalpyPerm(Element % NodeIndexes)) = & 
            PrevEnthalpy_h(EnthalpyPerm(Element % NodeIndexes)) + dT(1:n)

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
!       Check if this element belongs to a body where Enthalpy_h 
!       has been calculated
!------------------------------------------------------------------------------
        Element => Solver % Mesh % Elements(t)

        NodeIndexes => Element % NodeIndexes
        IF ( ANY( EnthalpyPerm( NodeIndexes ) <= 0 ) ) CYCLE

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
          i = EnthalpyPerm( NodeIndexes(k) )
          DO j=1,SIZE(PhaseChangeIntervals,2)
            IF ( ( Enthalpy_h(i)  < PhaseChangeIntervals(1,j) .AND. &
                   PrevSolution(i) > PhaseChangeIntervals(2,j) ).OR. &
                 ( Enthalpy_h(i)  > PhaseChangeIntervals(2,j) .AND. &
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
   SUBROUTINE IntegOverA( BoundaryMatrix, BoundaryVector, &
     LOAD, NodalAlpha, Element, n, m, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LOAD(:),NodalAlpha(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n,  m

     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
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
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                  Basis,dBasisdx )

       s = SqrtElementMetric * S_Integ(t)
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
       Force = SUM( LOAD(1:n) * Basis )
       Alpha = SUM( NodalAlpha(1:n) * Basis )

       DO p=1,N
         DO q=1,M
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
                  s * Alpha * Basis(p) / m
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
    SUBROUTINE AddHeatGap( Solver, Element, STIFF, EnthalpyPerm )
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      REAL(KIND=dp) :: STIFF(:,:)
      INTEGER :: EnthalpyPerm(:)
      TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Parent,Left,Right
      INTEGER :: i,j,k,l, Ind(n)
      REAL(KIND=dp) :: x0,y0,z0,x,y,z
!------------------------------------------------------------------------------
      CALL FindGapIndexes( Element, Ind, n )
      DO i=1,n
        DO j=1,n
          k = EnthalpyPerm( Element % NodeIndexes(i) )
          l = EnthalpyPerm( Ind(j) )
          IF ( k > 0 .AND. l > 0 ) THEN
            CALL AddToMatrixElement( Solver % Matrix,k,l,-STIFF(i,j) )
          END IF
        END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE AddHeatGap
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  END SUBROUTINE EnthalpySolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION HeatBoundaryResidual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Radiation

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,n,l,t,DIM,Pn,En
     LOGICAL :: stat, Found

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), ExtEnthalpy_h(:), &
       TransferCoeff(:), EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
       dBasisdx(:,:), Enthalpy_h(:), Flux(:), NodalEmissivity(:)

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

     En = Edge % TYPE % NumberOfNodes
     Pn = Element % TYPE % NumberOfNodes

     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( Nodes % x(Pn), Nodes % y(Pn), Nodes % z(Pn) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( Enthalpy_h(Pn), Basis(Pn), ExtEnthalpy_h(En), &
        TransferCoeff(En), x(En), y(En), z(En), EdgeBasis(En), &
        dBasisdx(Pn,3), NodalConductivity(En), Flux(En), &
        NodalEmissivity(En) ) 

     DO l = 1,En
       DO k = 1,Pn
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
!                 'Enthalpy Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        s = ListGetConstReal( Model % BCs(j) % Values,'Enthalpy_h',Dirichlet )

!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Enthalpy Heat Flux', En, Edge % NodeIndexes, Found )

!       ...convective heat transfer:
!       ----------------------------
        TransferCoeff(1:En) =  ListGetReal( Model % BCs(j) % Values, &
          'Heat Transfer Coefficient', En, Edge % NodeIndexes, Found )

        ExtEnthalpy_h(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'External Enthalpy_h', En, Edge % NodeIndexes, Found )

!       ...black body radiation:
!       ------------------------
        Emissivity      = 0.0d0
        StefanBoltzmann = 0.0d0

        SELECT CASE(ListGetString(Model % BCs(j) % Values,'Radiation',Found))
           !------------------
           CASE( 'idealized' )
           !------------------

              NodalEmissivity(1:En) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', En, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:En) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:En)) / En

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

           !---------------------
           CASE( 'diffuse gray' )
           !---------------------

              NodalEmissivity(1:En) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', En, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:En) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:En)) / En

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

              ExtEnthalpy_h(1:En) =  ComputeRadiationLoad( Model, &
                      Mesh, Edge, Quant, Perm, Emissivity )
        END SELECT

!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Enthalpy Heat Diffusivity', Hwrk, En, Edge % NodeIndexes )

        NodalConductivity( 1:En ) = Hwrk( 1,1,1:En )

!       elementwise nodal solution:
!       ---------------------------
        Enthalpy_h(1:Pn) = Quant( Perm(Element % NodeIndexes) )

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0

        IntegStuff = GaussPoints( Edge )

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
              gx = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
              gy = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
              gz = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
      
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, gx, gy, gz )

              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
                 Basis, dBasisdx )
!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Conductivity = SUM( NodalConductivity(1:En) * EdgeBasis(1:En) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:En) * EdgeBasis(1:En) )

!          convective ...:
!          ----------------
           Residual = Residual + SUM(TransferCoeff(1:En) * EdgeBasis(1:En)) * &
                     ( SUM( Enthalpy_h(1:Pn) * Basis(1:Pn) ) - &
                       SUM( ExtEnthalpy_h(1:En) * EdgeBasis(1:En) ) )

!          black body radiation...:
!          -------------------------
           Residual = Residual + &
                Emissivity * StefanBoltzmann * &
                     ( SUM( Enthalpy_h(1:Pn) * Basis(1:Pn) ) ** 4 - &
                       SUM( ExtEnthalpy_h(1:En) * EdgeBasis(1:En) ) ** 4 )

!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,DIM
                 Residual = Residual + Conductivity  * &
                    SUM( dBasisdx(1:Pn,k) * Enthalpy_h(1:Pn) ) * Normal(k)

                 Gnorm = Gnorm + s * (Conductivity * &
                       SUM(dBasisdx(1:Pn,k) * Enthalpy_h(1:Pn)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,DIM
                 DO l=1,DIM
                    Residual = Residual + Metric(k,l) * Conductivity  * &
                       SUM( dBasisdx(1:Pn,k) * Enthalpy_h(1:Pn) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
                      SUM(dBasisdx(1:Pn,k) * Enthalpy_h(1:Pn) ) * Normal(l))**2
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

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF

!    Gnorm = EdgeLength * Gnorm
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( Enthalpy_h, Basis, ExtEnthalpy_h, TransferCoeff,  &
        x, y, z, EdgeBasis, dBasisdx, NodalConductivity, Flux, &
        NodalEmissivity ) 
!------------------------------------------------------------------------------
  END FUNCTION HeatBoundaryResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION HeatEdgeResidual( Model, Edge, Mesh, Quant, Perm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Edge
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element, Bndry

     INTEGER :: i,j,k,l,n,t,DIM,En,Pn
     LOGICAL :: stat, Found
     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), &
            EdgeBasis(:), Basis(:), dBasisdx(:,:), Enthalpy_h(:)

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
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
     END DO

     Grad = 0.0d0
!
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left
     n = Element % TYPE % NumberOfNodes

     Element => Edge % BoundaryInfo % Right
     n = MAX( n, Element % TYPE % NumberOfNodes )

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     En = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( NodalConductivity(En), EdgeBasis(En), Basis(n), &
        dBasisdx(n,3), x(En), y(En), z(En), Enthalpy_h(n) )

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
           u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
           v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
           w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )

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
           Pn = Element % TYPE % NumberOfNodes

           DO j = 1,En
              DO k = 1,Pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
           Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
           Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
             Basis, dBasisdx )
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                   'Enthalpy Heat Diffusivity', Hwrk,En, Edge % NodeIndexes )

           NodalConductivity( 1:En ) = Hwrk( 1,1,1:En )
           Conductivity = SUM( NodalConductivity(1:En) * EdgeBasis(1:En) )
!
!          Enthalpy_h at element nodal points:
!          ------------------------------------
           Enthalpy_h(1:Pn) = Quant( Perm(Element % NodeIndexes) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,DIM
              Grad(j,i) = Conductivity * SUM( dBasisdx(1:Pn,j) * Enthalpy_h(1:Pn) )
           END DO
        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0d0
        DO k=1,DIM
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,DIM
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
     END DO

     IF ( CoordinateSystemDimension() == 3 ) THEN
        EdgeLength = SQRT(EdgeLength)
     END IF
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, x, y, z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
     DEALLOCATE( NodalConductivity, EdgeBasis, Basis, dBasisdx, Enthalpy_h)

!------------------------------------------------------------------------------
  END FUNCTION HeatEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION HeatInsideResidual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,n,t,DIM

     LOGICAL :: stat, Found, Compressible
     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalCapacity(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
     REAL(KIND=dp), ALLOCATABLE :: Velo(:,:), Pressure(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Enthalpy_h(:), PrevTemp(:)
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
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

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
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( NodalDensity(n), NodalCapacity(n), NodalConductivity(n),       &
         Velo(3,n), Pressure(n), NodalSource(n), Enthalpy_h(n), PrevTemp(n), &
         Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3) )
!
!    Elementwise nodal solution:
!    ---------------------------
     Enthalpy_h(1:n) = Quant( Perm(Element % NodeIndexes) )
!
!    Check for time dep.
!    -------------------
     PrevTemp(1:n) = Enthalpy_h(1:n)
     dt = Model % Solver % dt
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Enthalpy_h', .TRUE. )
        PrevTemp(1:n) = Var % PrevValues(Var % Perm(Element % NodeIndexes),1)
     END IF
!
!    Material parameters: conductivity, heat capacity and density
!    -------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL ListGetRealArray( Material, &
                  'Enthalpy Heat Diffusivity', Hwrk,n, Element % NodeIndexes )

     NodalConductivity( 1:n ) = Hwrk( 1,1,1:n )

     NodalDensity(1:n) = ListGetReal( Material, &
            'Enthalpy Density', n, Element % NodeIndexes, Found )

     NodalCapacity(1:n) = 1.0 
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
           Pressure(1:n) = &
               Var % Values( Var % Perm(Element % NodeIndexes) )
        END IF

        ReferencePressure = ListGetConstReal( Material, &
                   'Reference Pressure' )

        SpecificHeatRatio = ListGetConstReal( Material, &
                   'Specific Heat Ratio' )

        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * NodalCapacity(1:n) * Enthalpy_h(1:n) )
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

           Velo(1,1:n) = ListGetReal( Material, &
              'Convection Velocity 1', n, Element % NodeIndexes, Found )

           Velo(2,1:n) = ListGetReal( Material, &
              'Convection Velocity 2', n, Element % NodeIndexes, Found )

           Velo(3,1:n) = ListGetReal( Material, &
              'Convection Velocity 3', n, Element % NodeIndexes, Found )

        !-----------------
        CASE( 'computed' )
        !-----------------

           Var => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              IF ( ALL( Var % Perm( Element % NodeIndexes ) > 0 ) ) THEN
                 Velo(1,1:n) = Var % Values(Var % Perm(Element % NodeIndexes))
   
                 Var => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(2,1:n) = Var % Values( &
                              Var % Perm(Element % NodeIndexes ) )
   
                 Var => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(3,1:n) = Var % Values( &
                             Var % Perm( Element % NodeIndexes ) )
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
     IF ( Found .AND. k > 0  ) THEN
        NodalSource(1:n) = ListGetReal( Model % BodyForces(k) % Values, &
               'Heat Source', n, Element % NodeIndexes, Found )
     END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0d0
     Area = 0.0d0

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Capacity     = SUM( NodalCapacity(1:n) * Basis(1:n) )
        Density      = SUM( NodalDensity(1:n) * Basis(1:n) )
        Conductivity = SUM( NodalConductivity(1:n) * Basis(1:n) )
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
        Residual = -Density * SUM( NodalSource(1:n) * Basis(1:n) )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,DIM
!
!             - grad(C).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Enthalpy_h(1:n) * dBasisdx(1:n,j) ) * &
                 SUM( NodalConductivity(1:n) * dBasisdx(1:n,j) )

!
!             - C div(grad(T)):
!             -------------------
!
              Residual = Residual - Conductivity * &
                 SUM( Enthalpy_h(1:n) * ddBasisddx(1:n,j,j) )
           END DO
        ELSE
           DO j=1,DIM
              DO k=1,DIM
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Enthalpy_h(1:n) * dBasisdx(1:n,j) ) * &
                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Conductivity * &
                    SUM( Enthalpy_h(1:n) * ddBasisddx(1:n,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,DIM
                    Residual = Residual + Metric(j,k) * Conductivity * &
                      Symb(j,k,l) * SUM( Enthalpy_h(1:n) * dBasisdx(1:n,l) )
                 END DO
              END DO
           END DO
        END IF

!       + \rho * c_p * (@T/@t + u.grad(T)):
!       -----------------------------------
        Residual = Residual + Density * Capacity *  &
           SUM((Enthalpy_h(1:n)-PrevTemp(1:n))*Basis(1:n)) / dt

        DO j=1,DIM
           Residual = Residual + &
              Density * Capacity * SUM( Velo(j,1:n) * Basis(1:n) ) * &
                    SUM( Enthalpy_h(1:n) * dBasisdx(1:n,j) )
        END DO


        IF ( Compressible ) THEN
!
!          + p div(u) or p u^j_{,j}:
!          -------------------------
!
           DO j=1,DIM
              Residual = Residual + &
                 SUM( Pressure(1:n) * Basis(1:n) ) * &
                      SUM( Velo(j,1:n) * dBasisdx(1:n,j) )

              IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO k=1,DIM
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
        DO i=1,DIM
           Fnorm = Fnorm + s * ( Density * &
             SUM( NodalSource(1:n) * Basis(1:n) ) ) ** 2
        END DO

        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
     END DO

!    Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm

     DEALLOCATE( NodalDensity, NodalCapacity, NodalConductivity,    &
         Velo, Pressure, NodalSource, Enthalpy_h, PrevTemp, Basis, &
         dBasisdx, ddBasisddx )
!------------------------------------------------------------------------------
  END FUNCTION HeatInsideResidual
!------------------------------------------------------------------------------

