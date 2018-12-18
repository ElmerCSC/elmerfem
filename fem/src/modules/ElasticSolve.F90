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
! *  Authors: Juha Ruokolainen, Mikko Lyly, Mika Malinen
! *  Email:   Juha.Ruokolainen@csc.fi, Mika.Malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Initializations for the primary solver: ElasticSolver 
!------------------------------------------------------------------------------
SUBROUTINE ElasticSolver_Init0( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: MixedFormulation, Found
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  MixedFormulation = GetLogical(SolverParams, 'Mixed Formulation', Found) .AND. &
      GetLogical(SolverParams, 'Neo-Hookean Material', Found)

  IF( MixedFormulation ) THEN
    CALL ListAddNewString( SolverParams, "Element", "p:2" )
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE ElasticSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE ElasticSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: dim, i, DOFs
  LOGICAL :: Found, AxialSymmetry, MixedFormulation
  LOGICAL :: CalculateStrains, CalculateStresses
  LOGICAL :: CalcPrincipalAngle, CalcPrincipal
  LOGICAL :: CalcPrincipalStress, CalcPrincipalStrain
  LOGICAL :: UseUMAT, OutputStateVars
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  AxialSymmetry = CurrentCoordinateSystem() == AxisSymmetric
  MixedFormulation = GetLogical(SolverParams, 'Mixed Formulation', Found) .AND. &
      GetLogical(SolverParams, 'Neo-Hookean Material', Found)

  IF ( .NOT. ListCheckPresent( SolverParams, 'Variable') ) THEN
    dim = CoordinateSystemDimension()
    IF (MixedFormulation) THEN
      DOFs = dim + 1
      SELECT CASE(dim)
      CASE(2)
        CALL ListAddString( SolverParams, 'Variable', 'MixedSol[Disp:2 Pres:1]' )
      CASE(3)
        CALL ListAddString( SolverParams, 'Variable', 'MixedSol[Disp:3 Pres:1]' )
      END SELECT
    ELSE
      DOFs = dim
      CALL ListAddString( SolverParams, 'Variable', 'Displacement' )
    END IF

    CALL ListAddInteger( SolverParams, 'Variable DOFs', DOFs )
  END IF

  CALL ListAddInteger( SolverParams,'Time derivative order', 2 )
  CALL ListAddNewLogical( SolverParams,'Bubbles in Global System',.TRUE.)
  CALL ListAddNewLogical( SolverParams,'Displace Mesh At Init',.TRUE.)

  CalculateStrains = GetLogical(SolverParams, 'Calculate Strains', Found)
  CalculateStresses = GetLogical(SolverParams, 'Calculate Stresses', Found)

  !-------------------------------------------------------------------------------
  ! If stress computation is requested somewhere, then enforce it:
  !--------------------------------------------------------------------------------
  IF( .NOT. CalculateStresses ) THEN
     CalculateStresses = ListGetLogicalAnyEquation( Model,'Calculate Stresses')
     IF ( CalculateStresses ) CALL ListAddLogical( SolverParams,'Calculate Stresses',.TRUE.)
  END IF

  CalcPrincipal = GetLogical(SolverParams, 'Calculate Principal', Found)
  CalcPrincipalAngle = GetLogical(SolverParams, 'Calculate PAngle', Found)
  IF (CalcPrincipalAngle) CalcPrincipal = .TRUE. ! Principal angle computation enforces component calculation

  !----------------------------------------------------------------------------------------------------
  CalcPrincipalStress = CalculateStresses .AND. CalcPrincipal
  CalcPrincipalStrain = CalculateStrains .AND. CalcPrincipal


  IF ( CalculateStresses ) THEN
     IF (AxialSymmetry) THEN
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Stress[Stress_xx:1 Stress_zz:1 Stress_yy:1 Stress_xy:1]' )
     ELSE
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Stress[Stress_xx:1 Stress_yy:1 Stress_zz:1 Stress_xy:1 Stress_yz:1 Stress_xz:1]' )
     END IF

     CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), 'vonMises' )

     IF (CalcPrincipalStress) THEN
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Principal Stress[Principal Stress:3]' )
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Tresca' )

        IF (CalcPrincipalAngle) THEN
           CALL ListAddString( SolverParams,&
                NextFreeKeyword('Exported Variable ',SolverParams), &
                '-dofs 9 Principal Angle' )
        END IF
     END IF
  END IF

  IF (CalculateStrains) THEN
     IF (AxialSymmetry) THEN
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Strain[Strain_xx:1 Strain_zz:1 Strain_yy:1 Strain_xy:1]' )
     ELSE
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Strain[Strain_xx:1 Strain_yy:1 Strain_zz:1 Strain_xy:1 Strain_yz:1 Strain_xz:1]' )
     END IF

     IF (CalcPrincipalStrain) THEN
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Principal Strain[Principal Strain:3]' )
             
     END IF
  END IF

  UseUMAT = ListGetLogical(SolverParams, 'Use UMAT', Found)
  IF (UseUMAT) THEN
    OutputStateVars = GetLogical(SolverParams, 'Output State Variables', Found)
    IF (OutputStateVars) THEN
      CALL ListAddString(SolverParams, NextFreeKeyword('Exported Variable ', SolverParams), &
          '-dofs 3 -ip StateVar[D1:1 D2:1 D3:1]' )
      CALL ListAddString(SolverParams, NextFreeKeyword('Exported Variable ', SolverParams), &
          '-dofs 9 -ip StateDir[PDir1_x:1 PDir1_y:1 PDir1_z:1 PDir2_x:1 PDir2_y:1 PDir2_z:1 PDir3_x:1 PDir3_y:1 PDir3_z:1]')      
    END IF
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE ElasticSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solver for the general non-linear elasticity equations.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ElasticSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------

  USE Adaptive
  USE DefUtils
  USE MaterialModels
  USE StressLocal
  
  IMPLICIT NONE

!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: StiffMatrix, PMatrix
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Variable_t), POINTER :: StressSol, TempSol, FlowSol, Var, StateSol, StateDir
  TYPE(ValueList_t), POINTER :: SolverParams, Material, BC, Equation, BodyForce
  TYPE(Nodes_t) :: ElementNodes, ParentNodes, FlowNodes
  TYPE(Element_t), POINTER :: CurrentElement, ParentElement, FlowElement
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  

  LOGICAL :: GotForceBC, GotFSIBC, GotIt, NewtonLinearization = .FALSE., Isotropic = .TRUE., &
       RotateModuli, LinearModel = .FALSE., MeshDisplacementActive, NeoHookeanMaterial = .FALSE., &
       AxialSymmetry
  LOGICAL :: UseUMAT, InitializeStateVars, OutputStateVars, HenckyStrain
  LOGICAL :: LargeDeflection
  LOGICAL :: MixedFormulation
  LOGICAL :: PseudoTraction, GlobalPseudoTraction
  LOGICAL :: PlaneStress, CalculateStrains, CalculateStresses
  LOGICAL :: CalcPrincipalAngle, CalcPrincipal
  LOGICAL :: CalcPrincipalStress, CalcPrincipalStrain
  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: CompressibilityDefined = .FALSE.
  LOGICAL :: NormalSpring, NormalTangential
  LOGICAL :: Converged, NoExternalLoads
  LOGICAL :: Parallel, Scanning

  INTEGER, POINTER :: TempPerm(:),StressPerm(:),PressPerm(:),NodeIndexes(:), &
          Indices(:), FlowPerm(:), AdjacentNodes(:)

  INTEGER :: dim,i,j,k,l,m,n,nd,nb,ntot,t,iter,NDeg,k1,k2,STDOFs,LocalNodes,istat
  INTEGER :: NewtonIter, NonlinearIter, MinNonlinearIter, FlowNOFNodes
  INTEGER :: CoordinateSystem
  INTEGER :: NPROPS, NSTATEV, MaxIntegrationPoints = 1, N_Gauss


  REAL(KIND=dp), POINTER :: Temperature(:),Pressure(:),Displacement(:), UWrk(:,:), &
       Work(:,:), ForceVector(:), Velocity(:,:), FlowSolution(:), SaveValues(:), &
       NodalStrain(:), NodalStress(:), VonMises(:), &
       PrincipalStress(:), PrincipalStrain(:), Tresca(:), PrincipalAngle(:)
  REAL(KIND=dp), POINTER :: PointwiseStateV(:,:), PointwiseStateV0(:,:), MaterialConstants(:,:)
  REAL(KIND=dp), POINTER :: TotalSol(:) => NULL()
  REAL(KIND=dp), POINTER CONTIG :: ValuesSaved(:) => NULL()
  REAL(KIND=dp), POINTER :: Pb(:)

  REAL(KIND=dp), ALLOCATABLE :: LocalMassMatrix(:,:),LocalStiffMatrix(:,:),&
       LocalDampMatrix(:,:),LoadVector(:,:),InertialLoad(:,:), Viscosity(:), LocalForce(:), &
       LocalTemperature(:),ElasticModulus(:,:,:),PoissonRatio(:), Density(:), &
       Damping(:), HeatExpansionCoeff(:,:,:),Alpha(:,:),Beta(:), &
       ReferenceTemperature(:),BoundaryDispl(:),LocalDisplacement(:,:), PrevSOL(:), &
       PrevLocalDisplacement(:,:), SpringCoeff(:,:,:), LocalExternalForce(:)
         
  REAL(KIND=dp) :: UNorm, TransformMatrix(3,3), Tdiff, Normal(3), s, UnitNorm, DragCoeff
  REAL(KIND=dp) :: Norm, NonlinTol, NonlinRes0, NonlinRes 

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0
#else
  REAL(KIND=dp) :: at,at0,CPUTime,RealTime
#endif

  CHARACTER(LEN=MAX_NAME_LEN) :: str, CompressibilityFlag
  CHARACTER(LEN=MAX_NAME_LEN) :: EquationName
!------------------------------------------------------------------------------
  SAVE LocalMassMatrix,LocalStiffMatrix,LocalDampMatrix,LoadVector,InertialLoad, Viscosity, &
       LocalForce,ElementNodes,ParentNodes,FlowNodes,Alpha,Beta, &
       LocalTemperature,AllocationsDone,ReferenceTemperature,BoundaryDispl, &
       ElasticModulus, PoissonRatio,Density,Damping,HeatExpansionCoeff, &
       LocalDisplacement, Velocity, Pressure, PrevSOL, CalculateStrains, CalculateStresses, &
       OutputStateVars, NodalStrain, NodalStress, VonMises, PrincipalStress, PrincipalStrain, &
       Tresca, PrincipalAngle, CalcPrincipalAngle, CalcPrincipal, &
       PrevLocalDisplacement, SpringCoeff, Indices
  SAVE NPROPS, NSTATEV, MaxIntegrationPoints, PointwiseStateV, PointwiseStateV0, &
      InitializeStateVars, TotalSol, LocalExternalForce
!-----------------------------------------------------------------------------------------------------
  INTERFACE
    FUNCTION ElastBoundaryResidual( Model,Edge,Mesh,Quant,Perm, Gnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END FUNCTION ElastBoundaryResidual

    FUNCTION ElastEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END FUNCTION ElastEdgeResidual

    FUNCTION ElastInsideResidual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END FUNCTION ElastInsideResidual
  END INTERFACE


  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  CALL Info( 'ElasticSolve', 'Starting Solver', Level=10 )
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

  SolverParams => GetSolverParams()
  Mesh => GetMesh()
  dim = CoordinateSystemDimension()
  CoordinateSystem = CurrentCoordinateSystem()
  AxialSymmetry = CoordinateSystem == AxisSymmetric
  Parallel = ParEnv % PEs > 1
  Scanning = ListGetString(Model % Simulation, 'Simulation Type', GotIt) == 'scanning'

  StressSol => Solver % Variable
  StressPerm     => StressSol % Perm
  STDOFs         =  StressSol % DOFs
  Displacement   => StressSol % Values
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  LocalNodes = COUNT( StressPerm > 0 )
  IF ( LocalNodes <= 0 ) RETURN

  TempSol => VariableGet( Mesh % Variables, 'Temperature' )
  IF ( ASSOCIATED( TempSol) ) THEN
     TempPerm    => TempSol % Perm
     Temperature => TempSol % Values
  END IF

  FlowSol => VariableGet( Mesh % Variables, 'Flow Solution' )
  IF ( ASSOCIATED( FlowSol) ) THEN
     FlowPerm => FlowSol % Perm
     k = SIZE( FlowSol % Values )
     FlowSolution => FlowSol % Values
  END IF

  MeshDisplacementActive = ListGetLogical( SolverParams, &
       'Displace Mesh', GotIt )
  IF ( .NOT. GotIt ) MeshDisplacementActive = .TRUE.

  IF ( AllocationsDone .AND. MeshDisplacementActive ) THEN
     CALL DisplaceMesh( Mesh, Displacement, -1, StressPerm, STDOFs, UpdateDirs=dim )
  END IF

  !-------------------------------------------------------------------------
  !    Check how material behaviour is defined: 
  !-------------------------------------------------------------------------
  UseUMAT = ListGetLogical( SolverParams, 'Use UMAT', GotIt )
  IF (UseUMAT .AND. TransientSimulation) CALL Fatal('ElasticSolve', &
      'UMAT version does not yet support transient simulation')

  NeoHookeanMaterial = ListGetLogical( SolverParams, 'Neo-Hookean Material', GotIt )
  IF (NeoHookeanMaterial) Isotropic = .TRUE.
  MixedFormulation = NeoHookeanMaterial .AND. &
      ListGetLogical( SolverParams, 'Mixed Formulation', GotIt )
  IF (MixedFormulation .AND. (STDOFs /= (dim + 1))) CALL Fatal('ElasticSolve', &
      'With mixed formulation variable DOFs should equal to space dimensions + 1')

  !------------------------------------------------------------------------------
  !     Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Mesh % Changed ) THEN
     N = Mesh % MaxElementDOFs

     IF ( AllocationsDone ) THEN
        DEALLOCATE( &
             BoundaryDispl, &
             ReferenceTemperature, &
             HeatExpansionCoeff, &
             LocalTemperature, &
             Pressure, Velocity, &
             ElasticModulus, PoissonRatio, &
             Density, Damping, &
             LocalForce, LocalExternalForce, Viscosity, &
             LocalMassMatrix,  &
             LocalStiffMatrix,  &
             LocalDampMatrix,  &
             LoadVector, InertialLoad, Alpha, Beta, &
             LocalDisplacement, &
             PrevLocalDisplacement, &
             SpringCoeff, &
             Indices)
     END IF

     ALLOCATE( &
          BoundaryDispl( N ), &
          ReferenceTemperature( N ), &
          HeatExpansionCoeff( 3,3,N ), &
          LocalTemperature( N ), &
          Pressure( N ), Velocity( 3,N ), &
          ElasticModulus( 6,6,N ), PoissonRatio( N ), &
          Density( N ), Damping( N ), &
          LocalForce( STDOFs*N ), LocalExternalForce( STDOFs*N ), Viscosity( N ), &
          LocalMassMatrix(  STDOFs*N,STDOFs*N ),  &
          LocalStiffMatrix( STDOFs*N,STDOFs*N ),  &
          LocalDampMatrix( STDOFs*N,STDOFs*N ),  &
          LoadVector( 4,N ), InertialLoad(3,N), Alpha( 3,N ), Beta( N ), &
          LocalDisplacement( 4,N ), &
          PrevLocalDisplacement( 4,N ), &
          SpringCoeff( N,3,3 ), &
          Indices(N), &
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'ElasticSolve',  'Memory allocation error.' )
     END IF

     IF (UseUMAT .AND. (.NOT. AllocationsDone) ) THEN
        !------------------------------------------------------------------------------
        ! The UMAT description of the material behavior may depend on a list of state
        ! variables that depend on the material point. This description is produced
        ! at each integration point. The following is used for finding the right size
        ! of an array holding this information and allocating the array.
        ! TO DO: Take into account the possibility Mesh % Changed = .TRUE. 
        !------------------------------------------------------------------------------
        M = GetNOFActive()
        DO t=1,M
           CurrentElement => GetActiveElement(t)
           Material => GetMaterial()

           NPROPS = GetInteger( Material, 'Number of Material Constants', GotIt)
           IF (.NOT. GotIt) CALL Fatal('ElasticSolve', &
                'Number of Material Constants for UMAT must be specified')
           NSTATEV = GetInteger( Material, 'Number of State Variables', GotIt)
           IF (.NOT. GotIt) CALL Fatal('ElasticSolve', &
                'Number of Material Constants for UMAT must be specified')

           IntegStuff = GaussPoints( CurrentElement )

           MaxIntegrationPoints = MAX( IntegStuff % n, MaxIntegrationPoints )     
        END DO

        ! ---------------------------------------------------------------------
        ! PointwiseStateV is now allocated for keeping the state variables as
        ! they evolve during the nonlinear iteration to obtain the solution
        ! at the new time level m+1. Allocation is done for the given number 
        ! of state variables + 9 additional variables which are three energy 
        ! variables and six stress components:
        ! ---------------------------------------------------------------------
        CALL AllocateArray(PointwiseStateV, M * MaxIntegrationPoints, NSTATEV+9)
        PointwiseStateV = 0.0d0
        ! ----------------------------------------------------------------------
        ! We also create a similar variable PointwiseStateV0 which keeps the
        ! the state variables that describe the material state corresponding to 
        ! the converged solution at the previous time level m. The right values
        ! of the state variables corresponding to the initial state can be found
        ! by making an extra UMAT call. Whether this call is needed is indicated
        ! by the last extra entry: a value < 0 means that the state variables have
        ! not yet been initiated by the extra call.
        ! ----------------------------------------------------------------------
        CALL AllocateArray(PointwiseStateV0, M * MaxIntegrationPoints, NSTATEV+10)
        PointwiseStateV0 = 0.0d0
        PointwiseStateV0(:,NSTATEV+10) = -1.0d0
        InitializeStateVars = GetLogical(SolverParams, 'Initialize State Variables', &
            GotIt)
     END IF

     !----------------------------------------------------------------
     ! Check whether strains and stresses are computed...
     !----------------------------------------------------------------
     CalculateStrains = GetLogical(SolverParams, 'Calculate Strains', GotIt )    
     CalculateStresses = GetLogical(SolverParams, 'Calculate Stresses', GotIt ) 

     IF (UseUMAT) THEN
        OutputStateVars = GetLogical(SolverParams, 'Output State Variables', GotIt)
        ! Principal tensors are not yet available:
        CalcPrincipal = .FALSE.
        CalcPrincipalAngle = .FALSE.
     ELSE
        OutputStateVars = .FALSE.
        CalcPrincipal = GetLogical(SolverParams, 'Calculate Principal', GotIt )     
        CalcPrincipalAngle = GetLogical(SolverParams, 'Calculate PAngle', GotIt )
        ! Principal angle computation enforces component calculation:
        IF (CalcPrincipalAngle) CalcPrincipal = .TRUE. 
     END IF
     AllocationsDone = .TRUE.
  END IF

  !---------------------------------------------------------------------------------------------------
  !    Set pointers to the variables containing the stress and strain fields:
  !--------------------------------------------------------------------------------------------------
  CalcPrincipalStress = CalculateStresses .AND. CalcPrincipal
  CalcPrincipalStrain = CalculateStrains .AND. CalcPrincipal
  IF ( CalculateStresses ) THEN
     Var => VariableGet( Mesh % Variables, 'Stress', .TRUE. )
     IF ( ASSOCIATED( Var ) ) THEN
        StressPerm  => Var % Perm
        NodalStress => Var % Values
     ELSE  
        CALL Fatal('ElasticSolver','Variable > Stress < does not exits!')
     END IF

     Var => VariableGet( Mesh % Variables, 'VonMises',.TRUE. )
     IF ( ASSOCIATED( Var ) ) THEN
        VonMises => Var % Values
     ELSE
        CALL Fatal('ElasticSolver','Variable > vonMises < does not exits!')
     END IF

     IF (CalcPrincipalStress) THEN
        Var => VariableGet( Mesh % Variables, 'Principal Stress',.TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           PrincipalStress => Var % Values
        ELSE                 
           CALL Fatal('ElasticSolver','Variable > Principal Stress < does not exits!')
        END IF

        Var => VariableGet( Mesh % Variables, 'Tresca',.TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Tresca => Var % Values
        ELSE
           CALL Fatal('ElasticSolver','Variable > Tresca < does not exits!')
        END IF

        IF (CalcPrincipalAngle) THEN
           Var => VariableGet( Mesh % Variables, 'Principal Angle' )                 
           IF ( ASSOCIATED( Var ) ) THEN
              PrincipalAngle => Var % Values
           ELSE
              CALL Fatal('ElasticSolver','Variable > Principal Angle < does not exits!')
           END IF
        END IF
     END IF
  END IF

  IF (CalculateStrains) THEN
     Var => VariableGet( Mesh % Variables, 'Strain' )
     IF ( ASSOCIATED( Var ) ) THEN
        NodalStrain => Var % Values
     ELSE
        CALL Fatal('ElasticSolver','Variable > Strain < does not exits!')
     END IF
     IF (CalcPrincipalStrain) THEN
        Var => VariableGet( Mesh % Variables, 'Principal Strain' )
        IF ( ASSOCIATED( Var ) ) THEN
           PrincipalStrain => Var % Values
        ELSE
           CALL Fatal('ElasticSolver','Variable > Principal Strain < does not exits!')
        END IF
     END IF
  END IF

  IF (UseUMAT .AND. OutputStateVars) THEN
    StateSol => VariableGet(Mesh % Variables, 'StateVar')
    IF ( .NOT. ASSOCIATED(StateSol) ) THEN
      CALL Fatal('ElasticSolver','Variable > StateVar < does not exits!')
    END IF
    StateDir => VariableGet(Mesh % Variables, 'StateDir')
    IF ( .NOT. ASSOCIATED(StateDir) ) THEN
      CALL Fatal('ElasticSolver','Variable > StateDir < does not exits!')
    END IF
  END IF

  ALLOCATE( PrevSOL(SIZE(Displacement)) )
  PrevSOL = Displacement
  IF (UseUMAT) THEN
    IF (.NOT. ASSOCIATED(StiffMatrix % BulkRHS)) &
        ALLOCATE(StiffMatrix % BulkRHS(SIZE(StiffMatrix % RHS)))
    StiffMatrix % BulkRHS = 0.0d0
    
    IF (.NOT. ASSOCIATED(TotalSol)) ALLOCATE( TotalSol(SIZE(Displacement)) )

    IF (Scanning .AND. .NOT. ASSOCIATED( StressSol % PrevValues )) THEN
      ALLOCATE(StressSol % PrevValues(SIZE(Displacement), 1))
      StressSol % PrevValues = 0.0d0
    END IF

    CALL ListAddLogical(SolverParams, 'Skip Compute Nonlinear Change', .TRUE.)
    HenckyStrain = ListGetLogical(SolverParams, 'Use Hencky strain', GotIt)
  END IF

  !------------------------------------------------------------------------------
  !    Do some additional initialization, and go for it
  !------------------------------------------------------------------------------
  NonlinearIter = ListGetInteger( SolverParams, &
       'Nonlinear System Max Iterations', GotIt )
  IF ( .NOT. GotIt ) NonlinearIter = 1
  NonlinTol = GetConstReal(SolverParams, 'Nonlinear System Convergence Tolerance')
  MinNonlinearIter = ListGetInteger( SolverParams, &
       'Nonlinear System Min Iterations', GotIt )

  LinearModel = ListGetLogical( SolverParams, &
       'Elasticity Solver Linear', GotIt )
  LargeDeflection = ListGetLogical(SolverParams, 'Large Deflection', GotIt)
  IF (.NOT. GotIt) LargeDeflection = .TRUE.

  IF (.NOT. LargeDeflection) HenckyStrain = .FALSE.

  GlobalPseudoTraction = GetLogical( SolverParams, 'Pseudo-Traction', GotIt)

  
  CALL DefaultStart()

  DO iter=1,NonlinearIter

     at  = CPUTime()
     at0 = RealTime()

     CALL Info( 'ElasticSolve', ' ', Level=4 )
     CALL Info( 'ElasticSolve', ' ', Level=4 )
     CALL Info( 'ElasticSolve', &
          '-------------------------------------', Level=4 )
     WRITE( Message, * ) 'ELASTICITY ITERATION   ', iter
     CALL Info( 'ElasticSolve', Message, Level=4 )
     CALL Info( 'ElasticSolve', &
          '-------------------------------------', Level=4 )
     CALL Info( 'ElasticSolve', ' ', Level=4 )
     CALL Info( 'ElasticSolve', 'Starting assembly...', Level=4 )

     IF (UseUMAT) TotalSol(:) = Displacement(:)

     !------------------------------------------------------------------------------
100  CALL DefaultInitialize()
     !------------------------------------------------------------------------------
     DO t=1,GetNOFActive()

        IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                (Solver % NumberOfActiveElements-t) / &
                (1.0*Solver % NumberOfActiveElements)), ' % done'

           CALL Info( 'ElasticSolve', Message, Level=5 )
           at0 = RealTime()
        END IF

        CurrentElement => GetActiveElement(t)
        CALL GetElementNodes(ElementNodes, CurrentElement)
        NodeIndexes => CurrentElement % NodeIndexes

        n = GetElementNOFNodes()
        nd = GetElementDOFs( Indices )
        nb = GetElementNOFBDOFs()
        ntot = nd + nb

        !-----------------------------------------------------------------------------------
        !        Get the material parameters relating to the constitutive law:
        !------------------------------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()

        PlaneStress = GetLogical( Equation, 'Plane Stress', GotIt )
        PoissonRatio = 0.0d0

        IF (UseUMAT) THEN
           CALL GetConstRealArray( Material, MaterialConstants, 'Material Constants', GotIt)
           IF ( SIZE(MaterialConstants,1) /= NPROPS) &
                CALL Fatal('ElasticSolve','Check the size of Material Constants array')
        ELSE
           IF (NeoHookeanMaterial) THEN
              ElasticModulus(1,1,1:n) = ListGetReal( Material, &
                   'Youngs Modulus', n, NodeIndexes, GotIt )
           ELSE
              CALL InputTensor( ElasticModulus, Isotropic, &
                   'Youngs Modulus', Material, n, NodeIndexes )
              !------------------------------------------------------------------------------
              ! Check whether the rotation transformation of elastic modulus is necessary...
              !------------------------------------------------------------------------------
              RotateModuli = GetLogical( Material, 'Rotate Elasticity Tensor', GotIt )
              IF ( RotateModuli ) THEN
                 DO i=1,3
                    RotateModuli = .FALSE.
                    IF( i == 1 ) THEN
                       CALL GetConstRealArray( Material, UWrk, &
                            'Material Coordinates Unit Vector 1', GotIt, CurrentElement )
                    ELSE IF( i == 2 ) THEN
                       CALL GetConstRealArray( Material, UWrk, &
                            'Material Coordinates Unit Vector 2', GotIt, CurrentElement )
                    ELSE                
                       CALL GetConstRealArray( Material, UWrk, &
                            'Material Coordinates Unit Vector 3', GotIt, CurrentElement )
                    END IF

                    IF( GotIt ) THEN
                       UnitNorm = SQRT( SUM( Uwrk(1:3,1)**2 ) )
                       IF( UnitNorm < EPSILON( UnitNorm ) ) THEN
                          CALL Fatal('ElasticSolve','Given > Material Coordinate Unit Vector < too short!')
                       END IF
                       TransformMatrix(i,1:3) = Uwrk(1:3,1) / UnitNorm  
                       RotateModuli = .TRUE.
                    END IF
                    IF( .NOT. RotateModuli  ) CALL Fatal( 'ElasticSolve', &
                         'No unit vectors found but > Rotate Elasticity Tensor < set True?' )
                 END DO
              END IF
           END IF
           IF (Isotropic) PoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio' )
        END IF

        HeatExpansionCoeff = 0.0D0
        DO i=1,3
           HeatExpansionCoeff(i,i,1:n) = GetReal( Material,'Heat Expansion Coefficient', GotIt )
        END DO
        ReferenceTemperature(1:n) = GetReal( Material, 'Reference Temperature', GotIt )

        Density(1:n) = GetReal( Material, 'Density', GotIt )
        Damping(1:n) = GetReal( Material, 'Damping' ,GotIt )
        !------------------------------------------------------------------------------
        !        Set body forces
        !------------------------------------------------------------------------------
        BodyForce => GetBodyForce()

        LoadVector = 0.0D0
        InertialLoad = 0.0D0

        IF ( ASSOCIATED(BodyForce) ) THEN
           LoadVector(1,1:n) = GetReal( BodyForce, 'Stress Bodyforce 1', GotIt )
           LoadVector(2,1:n) = GetReal( BodyForce, 'Stress Bodyforce 2', GotIt )
           IF ( dim > 2 ) THEN
              LoadVector(3,1:n) = GetReal( BodyForce, 'Stress Bodyforce 3', GotIt )
           END IF

           InertialLoad(1,1:n) = GetReal( BodyForce, 'Inertial Bodyforce 1', GotIt )
           InertialLoad(2,1:n) = GetReal( BodyForce, 'Inertial Bodyforce 2', GotIt )

           IF ( dim > 2 ) THEN
              InertialLoad(3,1:n) = GetReal(  BodyForce, 'Inertial Bodyforce 3', GotIt )
           END IF
        END IF

        !------------------------------------------------------------------------------
        !        Get values of field variables:
        !------------------------------------------------------------------------------
        LocalTemperature = 0.0D0
        IF (UseUMAT) THEN
           IF ( ASSOCIATED(TempSol) ) THEN
              DO i=1,n
                 k = TempPerm(NodeIndexes(i))
                 LocalTemperature(i) = Temperature(k)
              END DO
           ELSE
              DO i=1,n
                 LocalTemperature(i) = ReferenceTemperature(i)
              END DO
           END IF
        ELSE
           IF ( ASSOCIATED(TempSol) ) THEN
              DO i=1,n
                 k = TempPerm(NodeIndexes(i))
                 LocalTemperature(i) = Temperature(k) - ReferenceTemperature(i)
              END DO
           END IF
        END IF

        LocalDisplacement = 0.0D0
        DO i=1,nd
           k = StressPerm(Indices(i))
           DO j=1,STDOFs
              LocalDisplacement(j,i) = Displacement(STDOFs*(k-1)+j)
           END DO
        END DO

        ! ----------------------------------------------------------------
        ! Some material models may need the displacement field at the
        ! previous time/load step
        ! ----------------------------------------------------------------
        PrevLocalDisplacement = 0.0D0
        IF (UseUMAT) THEN
          IF (TransientSimulation) THEN
            DO i=1,nd
              k = StressPerm(Indices(i))
              DO j=1,STDOFs
                PrevLocalDisplacement(j,i) = Solver % Variable % PrevValues(STDOFs*(k-1)+j,3)
              END DO
            END DO
          ELSE
            IF (Scanning) THEN
              DO i=1,nd
                k = StressPerm(Indices(i))
                DO j=1,STDOFs
                  PrevLocalDisplacement(j,i) = Solver % Variable % PrevValues(STDOFs*(k-1)+j,1)
                END DO
              END DO
            END IF
          END IF
        END IF

        IF ( LinearModel ) LocalDisplacement = 0.0d0

        !-------------------------------------------------------------------------------------------
        !        Select subroutine to integrate the element matrix and vector
        !-------------------------------------------------------------------------------------------
        IF ( CoordinateSystem == Cartesian .OR. AxialSymmetry) THEN
           IF (UseUMAT) THEN
              ! ------------------------------------------------------------------------------
              ! This branch assumes that the material behavior is defined 
              ! via an umat subroutine. The umat routine should specify
              ! a material response function which gives the Cauchy stress
              ! as a function of the strain tensor and state variables.
              !-------------------------------------------------------------------------------
              CALL LocalMatrixWithUMAT(LocalMassMatrix, LocalDampMatrix, &
                   LocalStiffMatrix, LocalForce, LocalExternalForce, dt, LoadVector, InertialLoad, &
                   MaterialConstants, NPROPS, PointwiseStateV, PointwiseStateV0, NSTATEV, &
                   MaxIntegrationPoints, InitializeStateVars, Density, Damping, AxialSymmetry, &
                   PlaneStress, LargeDeflection, HenckyStrain, CurrentElement, n, nd, ntot, STDOFs, &
                   ElementNodes, LocalDisplacement, PrevLocalDisplacement, LocalTemperature, &
                   t, Iter)

              !CALL Fatal( 'ElasticSolve', 'This version does not offer an umat interface' )

              ! ---------------------------------------------------------------------------
              ! Create a RHS vector which contains just the contribution of external loads
              ! for the purpose of nonlinear error estimation:
              ! ---------------------------------------------------------------------------
              IF (Iter == 1) THEN
                ValuesSaved => StiffMatrix % RHS
                StiffMatrix % RHS => StiffMatrix % BulkRHS
                CALL DefaultUpdateForce(LocalExternalForce)
                Solver % Matrix % RHS => ValuesSaved
              END IF

           ELSE
              !-------------------------------------------------------
              ! The following are used for handling cases where
              ! the material response function gives the second
              ! Piola-Kirchhoff stress
              !--------------------------------------------------------
              IF (NeoHookeanMaterial) THEN
                 CALL NeoHookeanLocalMatrix( LocalMassMatrix, LocalDampMatrix, &
                      LocalStiffMatrix, LocalForce, LoadVector, InertialLoad, ElasticModulus, &
                      PoissonRatio,Density,Damping,AxialSymmetry,PlaneStress,HeatExpansionCoeff, &
                      LocalTemperature,CurrentElement,n,ntot,ElementNodes,LocalDisplacement, &
                      MixedFormulation)
              ELSE
                 CALL LocalMatrix( LocalMassMatrix, LocalDampMatrix, &
                      LocalStiffMatrix,LocalForce, LoadVector, InertialLoad, ElasticModulus, &
                      PoissonRatio,Density,Damping,AxialSymmetry,PlaneStress,HeatExpansionCoeff, &
                      LocalTemperature,CurrentElement,n,ntot,ElementNodes,LocalDisplacement, &
                      Isotropic, RotateModuli, TransformMatrix)
              END IF
           END IF
        ELSE
           CALL Fatal('ElasticSolve', 'Unsupported coordinate system')
        END IF
        !------------------------------------------------------------------------------
        !        If time dependent simulation, add mass matrix to global 
        !        matrix and global RHS vector
        !------------------------------------------------------------------------------
        IF ( TransientSimulation ) THEN
           CALL Default2ndOrderTime( LocalMassMatrix, LocalDampMatrix, &
                LocalStiffMatrix, LocalForce )
        END IF
        !------------------------------------------------------------------------------
        !        Update global matrices from local matrices 
        !------------------------------------------------------------------------------
        CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
        !------------------------------------------------------------------------------
     END DO

     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     !     Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     DO t = 1,GetNOFBoundaryElements()
        CurrentElement =>  GetBoundaryElement(t)
        IF ( CurrentElement % TYPE % ElementCode == 101 ) CYCLE
        IF (.NOT. ActiveBoundaryElement()) CYCLE

        n  = GetElementNOFNodes()
        ntot = GetElementNOFDOFs()

        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
           LoadVector = 0.0D0
           Alpha      = 0.0D0
           Beta       = 0.0D0
           SpringCoeff = 0.0d0
           !------------------------------------------------------------------------------
           ! The components of surface forces 
           !------------------------------------------------------------------------------
           GotForceBC = .FALSE.
           LoadVector(1,1:n) = GetReal( BC, 'Surface Traction 1', GotIt )
           IF (.NOT. GotIt) LoadVector(1,1:n) = GetReal( BC, 'Force 1', GotIt )
           GotForceBC = GotForceBC .OR. GotIt

           LoadVector(2,1:n) = GetReal( BC, 'Surface Traction 2', GotIt )
           IF (.NOT. GotIt) LoadVector(2,1:n) = GetReal( BC, 'Force 2', GotIt )
           GotForceBC = GotForceBC .OR. GotIt

           LoadVector(3,1:n) = GetReal( BC, 'Surface Traction 3', GotIt )
           IF (.NOT. GotIt) LoadVector(3,1:n) = GetReal( BC, 'Force 3', GotIt )
           GotForceBC = GotForceBC .OR. GotIt

           Beta(1:n) = GetReal( BC, 'Normal Surface Traction', GotIt )
           IF (.NOT. GotIt) Beta(1:n) = GetReal( BC, 'Normal Force', gotIt )
           GotForceBC = GotForceBC .OR. GotIt

           GotForceBC = GotForceBC .OR. GetLogical( BC, 'Force BC', GotIt )

           SpringCoeff(1:n,1,1) =  GetReal( BC, 'Spring', NormalSpring )
           IF ( .NOT. NormalSpring ) THEN
              DO i=1,dim
                 SpringCoeff(1:n,i,i) = GetReal( BC, ComponentName('Spring',i), GotIt)
              END DO

              DO i=1,dim
                 DO j=1,dim
                    IF (ListCheckPresent(BC,'Spring '//TRIM(i2s(i))//i2s(j) )) &
                         SpringCoeff(1:n,i,j)=GetReal( BC, 'Spring '//TRIM(i2s(i))//i2s(j), GotIt)
                 END DO
              END DO
           END IF

           GotFSIBC = GetLogical( BC, 'FSI BC', GotIt )
           IF(.NOT. GotIt ) GotFSIBc = ASSOCIATED( FlowSol  ) 

           IF ( .NOT. GotForceBC .AND. .NOT. GotFSIBC .AND. ALL( SpringCoeff==0.0d0 ) ) CYCLE

           PseudoTraction = GetLogical( BC, 'Pseudo-Traction', GotIt)
           IF(.NOT. GotIt ) PseudoTraction = GlobalPseudoTraction
           !------------------------------------------------------------------------------

           ParentElement => CurrentElement % BoundaryInfo % Left

           IF ( .NOT. ASSOCIATED( ParentElement ) ) THEN
              ParentElement => CurrentElement % BoundaryInfo % Right
           ELSE
              IF ( ANY(StressPerm(ParentElement % NodeIndexes)==0 )) &
                   ParentElement => CurrentElement % BoundaryInfo % Right
           END IF

           nd = GetElementDOFs(Indices, ParentElement)
           CALL GetElementNodes( ParentNodes, ParentElement )

           LocalDisplacement = 0.0D0
           IF( .NOT. LinearModel ) THEN
              DO l=1,nd
                 k = StressPerm(Indices(l))
                 DO j=1,STDOFs
                    LocalDisplacement(j,l) = Displacement(STDOFs*(k-1)+j)
                 END DO
              END DO
           END IF

           NULLIFY( FlowElement )
           FlowNOFNodes = 1

           ! Note: Here the flow solution is not interpolated using the full p-basis
           IF ( GotFSIBC ) THEN
              FlowElement => CurrentElement % BoundaryInfo % Left

              IF ( .NOT. ASSOCIATED(FlowElement) ) THEN
                 FlowElement => CurrentElement % BoundaryInfo % Right
              ELSE
                 IF ( ANY(FlowPerm(FlowElement % NodeIndexes)==0 )) THEN
                    FlowElement => CurrentElement % BoundaryInfo % Right
                 END IF
              END IF

              IF ( ASSOCIATED(FlowElement) ) THEN
                 FlowNOFNodes = 0
                 FlowNOFNodes = FlowElement % TYPE % NumberOfNodes
                 AdjacentNodes => FlowElement % NodeIndexes

                 CALL GetElementNodes( FlowNodes, FlowElement )

                 DO j=1,FlowNOFNodes
                    k = StressPerm(AdjacentNodes(j))
                    IF ( k /= 0 ) THEN
                       k = STDOFs*(k-1)
                       FlowNodes % x(j) = FlowNodes % x(j) + PrevSOL( k+1 )

                       IF ( STDOFs > 1 ) &
                            FlowNodes % y(j) = FlowNodes % y(j) + PrevSOL( k+2 )

                       IF ( STDOFs > 2 ) &
                            FlowNodes % z(j) = FlowNodes % z(j) + PrevSOL( k+3 )
                    END IF
                 END DO

                 Velocity = 0.0D0
                 DO l=1,FlowNOFNodes
                    k = FlowPerm(AdjacentNodes(l))
                    DO j=1,FlowSol % DOFs-1
                       Velocity(j,l) = FlowSolution(FlowSol % DOFs*(k-1)+j)
                    END DO
                    Pressure(l) = FlowSolution(FlowSol % DOFs*k)
                 END DO

                 j = ListGetInteger( Model % Bodies(FlowElement % BodyId) &
                      % Values,'Material', minv=1, maxv=Model % NumberOFMaterials )
                 Material => Model % Materials(j) % Values
                 
                 Viscosity(1:FlowNOFNodes) = ListGetReal( &
                     Material,'Viscosity',FlowNOFNodes,AdjacentNodes,gotIt )
                 
                 CompressibilityFlag = ListGetString( Material, &
                     'Compressibility Model', GotIt )
                 
                 CompressibilityDefined = .FALSE.
                 IF ( GotIt ) THEN
                   IF ( CompressibilityFlag /= 'incompressible' ) THEN 
!.AND. &
!                       CompressibilityFlag /= 'artificial compressible') THEN
                     CompressibilityDefined = .TRUE.
                   END IF
                 END IF
                 
                 DragCoeff = ListGetCReal( BC,'FSI Drag Multiplier',GotIt)
                 IF(GotIt) THEN
                   Viscosity(1:FlowNOFNodes) = DragCoeff * Viscosity(1:FlowNOFNodes) 
                 END IF

              END IF
           END IF

           NormalTangential = GetLogical( BC, 'Normal-Tangential ' // & 
                GetVarName(Solver % Variable), GotIt )

           IF ( CoordinateSystem == Cartesian .OR. AxialSymmetry) THEN
              CALL LocalBoundaryMatrix( LocalStiffMatrix, LocalForce, &
                   LoadVector, SpringCoeff, NormalSpring, Alpha, Beta, LocalDisplacement, &
                   CurrentElement, n, ntot, ParentElement, ParentElement % TYPE % NumberOfNodes, &
                   nd, ParentNodes, FlowElement, FlowNOFNodes, FlowNodes, Velocity,  &
                   Pressure, Viscosity, Density, CompressibilityDefined, AxialSymmetry, &
                   NormalTangential, PseudoTraction, MixedFormulation, LargeDeflection)

              IF (UseUmat .AND. Iter == 1) THEN
                ! ---------------------------------------------------------------------------
                ! Update the RHS vector which contains just the contribution of external loads
                ! for the purpose of nonlinear error estimation:
                ! ---------------------------------------------------------------------------
                ValuesSaved => StiffMatrix % RHS
                StiffMatrix % RHS => StiffMatrix % BulkRHS
                CALL DefaultUpdateForce(LocalForce)
                Solver % Matrix % RHS => ValuesSaved
              END IF              

           ELSE
              CALL Fatal('ElasticSolve', 'Unsupported coordinate system')
           END IF

           !------------------------------------------------------------------------------
           !           Update global matrices from local matrices (will also affect
           !           LocalStiffMatrix and LocalForce if transient simulation is on).
           !------------------------------------------------------------------------------

           IF ( TransientSimulation ) THEN
              LocalDampMatrix = 0._dp
              LocalMassMatrix = 0._dp

              CALL Default2ndOrderTime( LocalMassMatrix, LocalDampMatrix, &
                   LocalStiffMatrix, LocalForce )
           END IF

           CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )              
        END IF
     END DO
     !------------------------------------------------------------------------------
     CALL DefaultFinishBoundaryAssembly()
     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     IF (UseUMAT) THEN
       ! ---------------------------------------------------------------------------------
       ! Now the solution variable is the solution increment while the sif-file specifies
       ! the Dirichlet BCs for the complete field. Modify BCs so that the right BC
       ! is obtained for the solution increment.
       ! ---------------------------------------------------------------------------------
       IF (ALLOCATED(StiffMatrix % ConstrainedDOF)) THEN
         DO i=1,StiffMatrix % NumberOfRows
           IF (StiffMatrix % ConstrainedDOF(i)) THEN
             StiffMatrix % DValues(i) = StiffMatrix % DValues(i) - Displacement(i)
           END IF
         END DO
         CALL EnforceDirichletConditions(Solver, StiffMatrix, ForceVector)
       END IF

       ! The initial guess for the displacement increment:
       Displacement = 0.0d0

       ! ---------------------------------------------------------------------------------
       ! Check whether the nonlinear iteration can be terminated:
       ! ---------------------------------------------------------------------------------
       IF (Iter == 1) THEN

         IF (Parallel) THEN
           IF (.NOT. ASSOCIATED(StiffMatrix % ParMatrix)) &
               CALL ParallelInitMatrix(Solver, StiffMatrix)

           PMatrix => StiffMatrix % ParMatrix % SplittedMatrix % InsideMatrix
           IF (.NOT. ASSOCIATED(PMatrix % RHS)) &
               ALLOCATE(PMatrix % RHS(PMatrix % NumberOfRows))

           ! Temporarily set the parallel rhs vector to be the plain source vector:
           CALL ParallelUpdateRHS(StiffMatrix, StiffMatrix % BulkRHS)
           Pb => PMatrix % RHS
           Norm = MAXVAL(ABS(Pb))
           Norm = ParallelReduction(Norm,2)
         ELSE
           Norm = MAXVAL(ABS(StiffMatrix % BulkRHS(:)))
         END IF

         NoExternalLoads = Norm < AEPS       
         IF (NoExternalLoads) THEN
           ! This appears to be a purely BC-loaded case, switch to using a different criterion
           ! (use absolute norm, this can be hard ...):
           CALL Info('ElasticSolver', 'No pressure external loads ... ', Level=4)
           CALL Info('ElasticSolver', &
               'Switch to using absolute norm in the nonlinear error estimation',  Level=4)
           CALL Info('ElasticSolver', &
               'This may give a hard stopping criterion',  Level=4)
           NonlinRes0 = 1.0d0
         ELSE
           ! Compute the 2-norm of the external load vector
           IF (Parallel)  THEN
             Norm = 0.0d0
             DO i=1,PMatrix % NumberOfRows
               Norm = Norm + Pb(i)**2
             END DO
             NonlinRes0 = SQRT(ParallelReduction(Norm))
           ELSE
             NonlinRes0 = SQRT(SUM(Solver % Matrix % BulkRHS(:)**2))
           END IF
         END IF
       END IF


       IF (Parallel) THEN
         ! Employ BulkRHS vector to estimate the size of the current residual (RHS):
         Solver % Matrix % BulkRHS = Solver % Matrix % RHS
         CALL ParallelUpdateRHS(StiffMatrix, Solver % Matrix % BulkRHS)
         Norm = 0.0d0
         DO i=1,PMatrix % NumberOfRows
           Norm = Norm + Pb(i)**2
         END DO
         NonlinRes = SQRT(ParallelReduction(Norm)) / NonlinRes0
       ELSE
         NonlinRes = SQRT(SUM(StiffMatrix % RHS(:)**2)) / NonlinRes0
       END IF
       WRITE(Message,'(a,I4,ES12.3)') 'Residual for nonlinear iterate', &
           Iter-1, NonLinRes
       CALL Info('ElasticitySolver', Message, Level=3)        
       IF (NonlinRes < NonlinTol .AND. (iter-1) >= MinNonlinearIter) THEN
         WRITE(Message,'(a)') 'Nonlinear iteration is terminated succesfully'
         CALL Info('ElasticitySolver', Message, Level=3)
         ! Save the state variables corresponding to the converged nonlinear
         ! solution to the array holding the previous solution state:
         PointwiseStateV0(:,1:NStatev+9) = PointwiseStateV(:,1:NStateV+9)
         Displacement(:) = TotalSol(:)
         IF (Scanning) StressSol % PrevValues(:,1) = Displacement(:)
         EXIT
       END IF

     ELSE
       IF ( DefaultLinesearch( Converged ) ) GOTO 100
       IF( iter >= MinNonlinearIter .AND. Converged ) EXIT
     END IF

     !------------------------------------------------------------------------------
     !     Solve the system and check for convergence
     !------------------------------------------------------------------------------
     UNorm = DefaultSolve()

     IF (UseUmat) THEN
       Displacement(:) = TotalSol(:) + Displacement(:)
       IF (iter==NonlinearIter) THEN
         WRITE(Message,'(a)') 'The maximum of nonlinear iterations reached: Terminating...'
         CALL Info('ElasticitySolver', Message, Level=3)        
         ! Save the state variables corresponding to the converged nonlinear
         ! solution to the array holding the previous solution state:
         PointwiseStateV0(:,1:NStatev+9) = PointwiseStateV(:,1:NStateV+9)
         IF (Scanning) StressSol % PrevValues(:,1) = Displacement(:)
         EXIT
       END IF
     ELSE
       !----------------------------------------------------------------------------------
       IF ( ( Solver % Variable % NonlinConverged == 1 .OR. iter==NonlinearIter ) .AND. &
           ( iter >= MinNonlinearIter ) ) EXIT
     END IF

  !------------------------------------------------------------------------------
  END DO ! of nonlinear iter
  !------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !   Perform strain and stress computation...
  !-----------------------------------------------------------------------------
  IF (CalculateStrains .OR. CalculateStresses) THEN
     CALL Info('ElasticSolve','Computing postprocessing fields')
     IF (UseUMAT) THEN
        CALL GenerateStressVariable(PointwiseStateV, NodalStress, StressPerm, MaxIntegrationPoints, &
             NStateV, CalculateStresses, AxialSymmetry)

        CALL GenerateStrainVariable(Displacement, NodalStrain, StressPerm, CalculateStrains, &
            AxialSymmetry, LargeDeflection)
     ELSE
        CALL ComputeStressAndStrain( Displacement, NodalStrain, NodalStress, VonMises, StressPerm, &
             PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle, AxialSymmetry, NeoHookeanMaterial, &
             CalculateStrains, CalculateStresses, CalcPrincipal, CalcPrincipalAngle, MixedFormulation)
     END IF
  END IF

  !-----------------------------------------------------------------------------
  !   Write the state variable solution ...
  !-----------------------------------------------------------------------------
  IF (UseUMAT .AND. OutputStateVars) THEN
    CALL GenerateStateVariable(PointwiseStateV, MaxIntegrationPoints, NStateV, StateSol, StateDir)
  END IF


  IF ( ListGetLogical(SolverParams, 'Adaptive Mesh Refinement', GotIt) ) THEN
     IF (UseUmat .OR. NeoHookeanMaterial) THEN
        CALL Info('ElasticSolve','Adaptive Mesh Refinement is not available') 
     ELSE
        CALL RefineMesh( Model, Solver, Displacement, StressPerm, &
             ElastInsideResidual, ElastEdgeResidual, ElastBoundaryResidual )

        IF ( MeshDisplacementActive ) THEN
           StressSol => Solver % Variable
           IF ( .NOT.ASSOCIATED( Mesh, Model % Mesh ) ) &
                CALL DisplaceMesh( Mesh, StressSol % Values, 1, &
                StressSol % Perm, StressSol % DOFs, .FALSE. )
        END IF
     END IF
  END IF

  IF ( MeshDisplacementActive ) THEN
     CALL Info('ElasticSolve','Displacing the mesh with computed displacement field')
     CALL DisplaceMesh( Mesh, Displacement, 1, StressPerm, STDOFs, .FALSE., dim )
  END IF

  DEALLOCATE( PrevSOL )

  CALL DefaultFinish()
  
  CALL Info('ElasticSolver','All done',Level=4)
  CALL Info('ElasticSolver','------------------------------------------',Level=4)

!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
! This subroutine uses the subroutine umat (Abaqus software convention for 
! defining a user-supplied material model) to get the material model.
! This subroutine assumes that a stress response function for the Cauchy
! stress is supplied (originally Elmer has employed Piola-Kirchhoff stresses).
! This is still a development version.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixWithUMAT(MassMatrix, DampMatrix, StiffMatrix, ForceVector, &
       ExternalForceVector, dt, LoadVector, InertialLoad, MaterialConstants, &
       NrInProps, PointwiseStateV, PointwiseStateV0, NStateV, MaxMaterialPoints, &
       InitializeStateVars, NodalDensity, NodalDamping, AxialSymmetry, PlaneStress, &
       LargeDeflection, HenckyStrain, Element, n, nd, ntot, dofs, Nodes, NodalDisplacement, &
       PrevNodalDisplacement, NodalTemperature, ElementIndex, IterationIndex)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MassMatrix(:,:), DampMatrix(:,:), StiffMatrix(:,:)
    REAL(KIND=dp) :: ExternalForceVector(:)
    REAL(KIND=dp) :: ForceVector(:), dt, LoadVector(:,:), InertialLoad(:,:)
    REAL(KIND=dp), POINTER :: MaterialConstants(:,:)
    INTEGER :: NrInProps
    REAL(KIND=dp), POINTER :: PointwiseStateV(:,:), PointwiseStateV0(:,:)
    INTEGER :: NStateV, MaxMaterialPoints
    LOGICAL :: InitializeStateVars
    REAL(KIND=dp) :: NodalDensity(:), NodalDamping(:)
    LOGICAL :: AxialSymmetry, PlaneStress, LargeDeflection
    LOGICAL :: HenckyStrain
    TYPE(Element_t) :: Element
    INTEGER :: n, nd, ntot, dofs
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: NodalDisplacement(:,:), PrevNodalDisplacement(:,:)
    REAL(KIND=dp) :: NodalTemperature(:)
    INTEGER :: ElementIndex
    INTEGER :: IterationIndex     ! The iteration index to resolve the nonlinearity
    !------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    REAL(KIND=dp) :: SymBasis1(3,3), SymBasis2(3,3), SymBasis3(3,3)
    REAL(KIND=dp) :: SymBasis4(3,3), SymBasis5(3,3), SymBasis6(3,3)
    REAL(KIND=dp) :: Identity(3,3)

    REAL(KIND=dp) :: Basis(ntot), dBasis(ntot,3), SqrtElementMetric
    REAL(KIND=dp) :: B(6,ntot*dofs)
    REAL(KIND=dp) :: Force(3), InertialForce(3)
    REAL(KIND=dp) :: Density, Damping, Temperature
    REAL(KIND=dp) :: Strain(3,3), Strain0(3,3)
    REAL(KIND=dp) :: Stress(3,3), Stress0(3,3), Stress1(3,3), Stress2(3,3)
    REAL(KIND=dp) :: dStress1(3,3)
    REAL(KIND=dp) :: Grad(3,3), Grad0(3,3)
    REAL(KIND=dp) :: InvDefG(3,3), DetDefG
    REAL(KIND=dp) :: InvDefG0(3,3), DetDefG0
    REAL(KIND=dp) :: C(3,3), InvC(3,3)
!    REAL(KIND=dp) :: EigenC(3,3)

    REAL(KIND=dp) :: WorkTensor1(3,3), WorkTensor2(3,3), WorkTensor3(3,3)
    REAL(KIND=dp) :: WorkVec1(1:6,1), WorkVec2(1:6,1)
    REAL(KIND=dp) :: s, u, v, w, r

    INTEGER :: i, j, k, l, p, q, t, dim, cdim, totdofs

    LOGICAL :: stat

    ! -----------------------------------------------------------------------------
    ! Variables for calling an UMAT subroutine that defines the material behaviour.
    ! The commenting with !* means that the variable is supported by Elmer. All
    ! variables have sufficient sizes for 3-D cases. For the use of these variables
    ! see also the definition of the subroutine umat.
    ! -----------------------------------------------------------------------------
    DOUBLE PRECISION :: StressVec(6)          !*
    DOUBLE PRECISION :: StateV(NStateV)       !* 
    DOUBLE PRECISION :: StressDer(6,6)        !*
    DOUBLE PRECISION :: EnergyElast           !* 
    DOUBLE PRECISION :: EnergyPlast           !* 
    DOUBLE PRECISION :: EnergyVisc            !* 
    DOUBLE PRECISION :: rpl
    DOUBLE PRECISION :: ddsddt(6)
    DOUBLE PRECISION :: drplde(6)
    DOUBLE PRECISION :: drpldt
    DOUBLE PRECISION :: stran(6)              !*
    DOUBLE PRECISION :: dstran(6)             !*
    DOUBLE PRECISION :: TimeAtStep(2)         !*
    DOUBLE PRECISION :: dtime                 !*
    DOUBLE PRECISION :: Temp                  !*
    DOUBLE PRECISION :: dTemp = 0.0d0         !  Zero for isothermal conditions
    DOUBLE PRECISION :: predef(1) = 0.0d0
    DOUBLE PRECISION :: dpred(1) = 0.0d0
    character(len=80) :: cmname               !*
    INTEGER :: ndi                            !*
    INTEGER :: nshr                           !*
    INTEGER :: ntens                          !*
    !    INTEGER :: NStateV                   !* Specified in the subroutine call
    DOUBLE PRECISION :: InProps(NrInProps)    !*
    !    INTEGER :: NrInProps                 !* Specified in the subroutine call
    DOUBLE PRECISION :: coords(3) = 0.0d0     !  TO DO: use this to provide the current coordinates
    DOUBLE PRECISION :: drot(3,3)     
    DOUBLE PRECISION :: pnewdt = 3.0d0
    DOUBLE PRECISION :: celent = 1.0d0        !* TO DO: use this to provide the element size 
    DOUBLE PRECISION :: DefG0(3,3)            !*
    DOUBLE PRECISION :: DefG(3,3)             !*
    !    INTEGER :: ElementIndex              !* Specified in the subroutine call
    INTEGER :: npt                            !*
    INTEGER :: layer = 1                     
    INTEGER :: kspt = 1
    INTEGER :: kstep = 1
    INTEGER :: kinc = 1

    !------------------------------------------------------------------------------
    ! Variables for computing eigenvector basis via calling the lapack dsyev
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: QWork(3,3), EigenVals(3), PriWork(102)
    INTEGER :: PriLWork=102, PriInfo=0
    !------------------------------------------------------------------------------

    IF (PlaneStress) CALL Fatal('LocalMatrixWithUMAT', 'Cannot yet handle plane stress case')
    IF (ntot > nd) CALL Fatal('LocalMatrixWithUMAT', 'Static condensation of bubbles is missing')

    ! ---------------------------------------------------------------------------------
    ! Six basis vectors for expressing symmetric tensors: The components of the 
    ! engineering strain vector [E_11 E_22 E_33 2E_12 2E_13 2E_23] are thus the 
    ! components of the strain tensor with respect to this basis
    ! ---------------------------------------------------------------------------------
    SymBasis1(1:3,1:3) = RESHAPE((/ 1,0,0,0,0,0,0,0,0 /),(/ 3,3 /))
    SymBasis2(1:3,1:3) = RESHAPE((/ 0,0,0,0,1,0,0,0,0 /),(/ 3,3 /))
    SymBasis3(1:3,1:3) = RESHAPE((/ 0,0,0,0,0,0,0,0,1 /),(/ 3,3 /)) 
    SymBasis4(1:3,1:3) = RESHAPE((/ 0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis5(1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0 /),(/ 3,3 /))
    SymBasis6(1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0 /),(/ 3,3 /))

    Identity(1:3,1:3) = RESHAPE((/ 1,0,0,0,1,0,0,0,1 /),(/ 3,3 /))

    ! The dimensionality (dim) of the state of stress:
    ! ---------------------------------------------------------------------------------
    cdim = CoordinateSystemDimension()
    IF (PlaneStress) THEN
      dim = cdim
    ELSE
      ! Axial symmetry, plane strain and 3-D:
      dim = 3
    END IF

    ! Define the array size for some umat variables: 
    ! ---------------------------------------------------------------------------------
    SELECT CASE(cdim)
    CASE(2)
      ! In plane stress case the third normal stress component is zero, but these size 
      ! definitions allow the umat subroutine to return the third normal strain which 
      ! cannot be reproduced from the 2-D displacement solution. 
      ! TO DO: indicate the plane stress condition via the material model name?
      ndi = 3
      nshr = 1
    CASE(3)
      ndi = 3
      nshr = 3
    END SELECT
    ntens = ndi + nshr

    totdofs = dofs * ntot

    ! --------------------------------------------------------------------------
    ! Specify some UMAT variables ...
    ! --------------------------------------------------------------------------
    DO i = 1,NrInProps
       InProps(i) = MaterialConstants(i,1)  
    END DO

    ! The material model name is used as a switch between some examples of implementation:
    IF (HenckyStrain) THEN
      cmname = 'hencky-stvenant-kirchhoff'//CHAR(0)
    ELSE
      IF (LargeDeflection) THEN
        cmname = 'stvenant-kirchhoff'//CHAR(0)
      ELSE
        cmname = 'linear isotropic'//CHAR(0)
      END IF
    END IF

    dtime = dt
    TimeAtStep(1:2) = GetTime() - dt
    DRot = Identity

    ! ------------------------------------
    ! Integration stuff
    ! ------------------------------------   
    IntegStuff = GaussPoints( Element )

    ForceVector = 0.0D0
    ExternalForceVector = 0.0D0
    StiffMatrix = 0.0D0
    MassMatrix  = 0.0D0
    DampMatrix  = 0.0d0

    DO t=1,IntegStuff % n
      B = 0.0d0
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      !------------------------------------------------------------------------------
      ! Basis function values & derivatives at the integration point
      !------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, SqrtElementMetric, &
          Basis, dBasis )

      s = SqrtElementMetric * IntegStuff % s(t)
      IF (AxialSymmetry) THEN
        r = SUM( Basis(1:n) * Nodes % x(1:n) )
        s = s * r
      END IF

      !---------------------------------------------------------------------------
      ! Force at integration point
      !----------------------------------------------------------------------------   
      DO i=1,cdim
        Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
        InertialForce(i) = SUM( InertialLoad(i,1:n)*Basis(1:n) )
      END DO

      !---------------------------------------------------------------------------
      ! Temperature, density and damping
      !---------------------------------------------------------------------------
      Temp = SUM( NodalTemperature(1:n)*Basis(1:n) )
      Density = SUM( NodalDensity(1:n)*Basis(1:n) )
      Damping = SUM( NodalDamping(1:n)*Basis(1:n) )

      !--------------------------------------------------------------------
      ! Compute the formulation variables (the displacement and deformation
      ! gradients):
      !--------------------------------------------------------------------
      Grad = 0.0d0
      Grad0 = 0.0d0
      IF (AxialSymmetry) THEN
        Grad(1,1) = SUM( NodalDisplacement(1,1:nd) * dBasis(1:nd,1) )
        Grad(1,3) = SUM( NodalDisplacement(1,1:nd) * dBasis(1:nd,2) ) 
        Grad(2,2) = 1.0d0/r * SUM( NodalDisplacement(1,1:nd) * Basis(1:nd) )
        Grad(3,1) = SUM( NodalDisplacement(2,1:nd) * dBasis(1:nd,1) )
        Grad(3,3) = SUM( NodalDisplacement(2,1:nd) * dBasis(1:nd,2) )

        Grad0(1,1) = SUM( PrevNodalDisplacement(1,1:nd) * dBasis(1:nd,1) )
        Grad0(1,3) = SUM( PrevNodalDisplacement(1,1:nd) * dBasis(1:nd,2) ) 
        Grad0(2,2) = 1.0d0/r * SUM( PrevNodalDisplacement(1,1:nd) * Basis(1:nd) )
        Grad0(3,1) = SUM( PrevNodalDisplacement(2,1:nd) * dBasis(1:nd,1) )
        Grad0(3,3) = SUM( PrevNodalDisplacement(2,1:nd) * dBasis(1:nd,2) )          
      ELSE
        ! Note that in the plane stress case we don't have means to create the fully
        ! consistent displacement gradient in the third direction:
        Grad(1:cdim,1:cdim) = MATMUL(NodalDisplacement(1:cdim,1:nd),dBasis(1:nd,1:cdim))
        Grad0(1:cdim,1:cdim) = MATMUL(PrevNodalDisplacement(1:cdim,1:nd),dBasis(1:nd,1:cdim))
      END IF
      DefG = Identity + Grad
      DefG0 = Identity + Grad0

      IF (LargeDeflection) THEN
        SELECT CASE( dim )
        CASE( 1 )
          DetDefG = DefG(1,1)
        CASE( 2 )
          DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
        CASE( 3 )
          DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
              DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
              DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
        END SELECT

        !-------------------------------------------------------------
        !  InvDefG will be the inverse of the deformation gradient
        !-------------------------------------------------------------
        InvDefG = DefG
        CALL InvertMatrix( InvDefG, dim )       
      END IF

      SELECT_STRAIN_MEASURE: IF (HenckyStrain) THEN
        ! TO DO:
        ! - warn if a Hencky umat has not been implemented
        ! --------------------------------------------
        ! The right Cauchy-Green deformation tensor 
        ! --------------------------------------------
        C = MATMUL( TRANSPOSE(DefG0), DefG0 )

        ! -----------------------------------------------------------
        ! Compute the spectral decomposition of C
        ! -----------------------------------------------------------
        DO i=1,3
          k = i
          DO j=k,3
            QWork(i,j) = C(i,j)
          END DO
        END DO
        CALL DSYEV('V', 'U', 3, QWork, 3, EigenVals, PriWork, PriLWork, PriInfo)
        IF (PriInfo /= 0) THEN
          CALL Fatal( 'ElasticSolve', 'DSYEV cannot generate eigen basis')          
        END IF

        Strain0 = 0.0d0
!        EigenC = MATMUL( TRANSPOSE(QWork), MATMUL(C,QWork) )       
        Strain0(1,1) = LOG(SQRT(EigenVals(1)))
        Strain0(2,2) = LOG(SQRT(EigenVals(2)))       
        Strain0(3,3) = LOG(SQRT(EigenVals(3)))
        ! Transform back to the original coordinates:
        Strain0 = MATMUL(QWork, MATMUL(Strain0,TRANSPOSE(QWork)))

        ! -----------------------------------------------------------
        ! Repeat for the current right Cauchy-Green deformation tensor:
        ! -----------------------------------------------------------
        C = MATMUL( TRANSPOSE(DefG), DefG )

        DO i=1,3
          k = i
          DO j=k,3
            QWork(i,j) = C(i,j)
          END DO
        END DO
        CALL DSYEV('V', 'U', 3, QWork, 3, EigenVals, PriWork, PriLWork, PriInfo)
        IF (PriInfo /= 0) THEN
          CALL Fatal( 'ElasticSolve', 'DSYEV cannot generate eigen basis')          
        END IF

        Strain = 0.0d0
        Strain(1,1) = LOG(SQRT(EigenVals(1)))
        Strain(2,2) = LOG(SQRT(EigenVals(2)))       
        Strain(3,3) = LOG(SQRT(EigenVals(3)))
        Strain = MATMUL(QWork, MATMUL(Strain,TRANSPOSE(QWork)))

        ! NOTE: The differentiation of the Hencky strain is done via a truncated 
        ! series expansion which may become inaccurate for large strains. 
        ! However, this inaccuracy should not break the consistency
        ! of the solution method: if nonlinear iterations converge, we shoud have
        ! a solution. That is, inaccuracy of the strain expansion has the effect
        ! that the Newton iteration is replaced by inexact Newton iteration.
        ! Currently no warnings are given for the possibility that the strain
        ! expansion may not be accurate:
        !
        !IF ( ANY(EigenVals(:) >= 2.0d0) .OR. ANY(Eigenvals(:) <= 0.5d0) ) &
        !     CALL Fatal( 'ElasticSolve', 'Series expansion for Hencky strain too short!')
      ELSE
        ! ---------------------------------------------------------------------------
        ! If the Hencky strain is not used, we use the standard material strain tensor 
        ! or its linearization
        ! ---------------------------------------------------------------------------
        Strain = 0.0d0
        Strain0 = 0.0d0
        Strain(1:dim,1:dim) = 0.5d0 * (Grad(1:dim,1:dim) + TRANSPOSE(Grad(1:dim,1:dim)))
        Strain0(1:dim,1:dim) = 0.5d0 * (Grad0(1:dim,1:dim) + TRANSPOSE(Grad0(1:dim,1:dim)))

        IF (LargeDeflection) THEN
          Strain(1:dim,1:dim) = Strain(1:dim,1:dim) + 0.5d0 * &
              MATMUL(TRANSPOSE(Grad(1:dim,1:dim)),Grad(1:dim,1:dim))
          Strain0(1:dim,1:dim) = Strain0(1:dim,1:dim) + 0.5d0 * &
              MATMUL(TRANSPOSE(Grad0(1:dim,1:dim)),Grad0(1:dim,1:dim))
        END IF

      END IF SELECT_STRAIN_MEASURE

      ! The umat (engineering) strain variable giving the strain before the increment:
      Stran(1) = Strain0(1,1)
      Stran(2) = Strain0(2,2)
      Stran(3) = Strain0(3,3)
      Stran(4) = 2.0d0 * Strain0(1,2)
      Stran(5) = 2.0d0 * Strain0(1,3)
      Stran(6) = 2.0d0 * Strain0(2,3)

      ! The umat variable giving the candidate for the strain increment:
      dStran(1) = Strain(1,1) - Strain0(1,1)
      dStran(2) = Strain(2,2) - Strain0(2,2)
      dStran(3) = Strain(3,3) - Strain0(3,3)
      dStran(4) = 2.0d0 * (Strain(1,2) - Strain0(1,2))
      dStran(5) = 2.0d0 * (Strain(1,3) - Strain0(1,3))
      dStran(6) = 2.0d0 * (Strain(2,3) - Strain0(2,3))
        
      ! ----------------------------------------------------------------------------
      ! Obtain the Cauchy stress and the stress response function derivative 
      ! via UMAT interface. If the state variables have not been initiated to correspond
      ! the initial state (stress-free initial condition is supposed), we first make
      ! an extra UMAT call to obtain the state variables if requested in the sif file.
      ! ----------------------------------------------------------------------------
      INITIALIZE_STATE_VARIABLES: IF ( InitializeStateVars .AND. &
          PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+10) < 0.0d0 ) THEN

        DO i=1,NStateV
          StateV(i) = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,i)
        END DO

        EnergyElast = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+1)
        EnergyPlast = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+2)
        EnergyVisc = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+3)
        StressVec(1:3+nshr) = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t, &
            NStateV+4:NStateV+6+nshr)
        
        ! We insert the identity tensor as the deformation gradient so the initial solution 
        ! should be the zero-displacement solution:
        stran = 0.0d0
        dstran = 0.0d0
        CALL UMAT(StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
            EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
            stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
            predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
            drot, pnewdt, celent, Identity, Identity, ElementIndex, t, layer, kspt, kstep, kinc)

        I = (ElementIndex-1)*MaxMaterialPoints+t
        IF ( ANY(StressVec(1:3+nshr) /= PointwiseStateV0(I,NStateV+4:NStateV+6+nshr)) ) THEN
          CALL Fatal('ElasticSolve','State variables initialization is changing stress')
        END IF

        ! Update the state variables storage (energy variables are not updated):
        DO i=1,NStateV
          PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,i) = StateV(i)
        END DO

        i = (ElementIndex-1)*MaxMaterialPoints+t
        PointwiseStateV0(i,NStateV+10) = 1.0d0
      END IF INITIALIZE_STATE_VARIABLES

      ! -----------------------------------------------------------------------------
      ! Perform the actual UMAT call. Before this, get the state variables and 
      ! the stress as specified at the previous time/load level for converged solution:
      ! -----------------------------------------------------------------------------
      DO i=1,NStateV
        StateV(i) = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,i)
      END DO
      EnergyElast = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+1)
      EnergyPlast = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+2)
      EnergyVisc = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t,NStateV+3)

      StressVec(1:3+nshr) = PointwiseStateV0((ElementIndex-1)*MaxMaterialPoints+t, &
          NStateV+4:NStateV+6+nshr)

      CALL UMAT(StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)

      ! ---------------------------------------------------------------------------
      ! Update data which gives the state variables corresponding to the current 
      ! nonlinear iterate.
      ! ---------------------------------------------------------------------------
      DO i=1,NStateV
        PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,i) = StateV(i)
      END DO
      PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,NStateV+1) = EnergyElast
      PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,NStateV+2) = EnergyPlast
      PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,NStateV+3) = EnergyVisc

      PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,NStateV+4:NStateV+6) = &
          StressVec(1:3)
      PointwiseStateV((ElementIndex-1)*MaxMaterialPoints+t,NStateV+7:NStateV+6+nshr) = &
          StressVec(4:3+nshr)

      STIFFMATRIX_FOR_CHOSEN_STRAIN: IF (.NOT. LargeDeflection) THEN
        ! ----------------------------------------
        ! Create the strain-displacement matrix B:
        ! ----------------------------------------
        DO p=1,ntot
          DO i=1,cdim
            B(i,(p-1)*dofs+i) = dBasis(p,i)
          END DO
          DO i=1,nshr
            SELECT CASE(i)
            CASE(1)
              B(ndi+i,(p-1)*dofs+1) = dBasis(p,2)
              B(ndi+i,(p-1)*dofs+2) = dBasis(p,1)
            CASE(2)
              B(ndi+i,(p-1)*dofs+1) = dBasis(p,3)
              B(ndi+i,(p-1)*dofs+3) = dBasis(p,1)
            CASE(3)
              B(ndi+i,(p-1)*dofs+2) = dBasis(p,3)
              B(ndi+i,(p-1)*dofs+3) = dBasis(p,2)
            END SELECT
          END DO
        END DO
       
        CALL StrainEnergyDensity(StiffMatrix, StressDer, B, ntens, totdofs, s)

        ! Internal force terms for the residual vector:
        ForceVector(1:totdofs) = ForceVector(1:totdofs) - MATMUL( TRANSPOSE(B(1:ntens,1:totdofs)), &
            StressVec(1:ntens) ) * s

        ! External forces:
        DO p=1,ntot
          DO i=1,cdim
            ExternalForceVector(dofs*(p-1)+i) = ExternalForceVector(dofs*(p-1)+i) +  ( &
                Basis(p) * Force(i) + Basis(p) * InertialForce(i) * Density) * s
            ForceVector(dofs*(p-1)+i) = ForceVector(dofs*(p-1)+i) +  ( &
                Basis(p) * Force(i) + Basis(p) * InertialForce(i) * Density) * s
          END DO
        END DO

      ELSE
        ! -------------------------------------------------------------------------
        ! THIS BRANCH CONTAINS AN IMPLEMENTATION FOR THE COMBINATION OF A NONLINEAR
        ! STRAIN AND CAUCHY STRESS. 
        !
        ! Now utilize the UMAT output to obtain the Newton linearization. First form
        ! the current Cauchy stress sigma_{n+1}^{(k)} (with k = IterationIndex-1 so 
        ! that IterationIndex = k+1 is associated with the variables to be solved) as 
        ! a symmetric tensor:
        !---------------------------------------------------------------------------
        Stress = StressVec(1)*SymBasis1 + StressVec(2)*SymBasis2 + &
            StressVec(3)*SymBasis3 + 2.0d0*StressVec(4)*SymBasis4 
        IF (nshr == 3) Stress = Stress + 2.0d0*StressVec(5)*SymBasis5 + &
            2.0d0*StressVec(6)*SymBasis6

        !--------------------------------------------------
        ! The first Piola-Kirchhoff stress
        !--------------------------------------------------
        Stress1 = DetDefG * MATMUL(Stress,TRANSPOSE(InvDefG))

        DO p = 1,ntot
          DO i = 1,cdim
            !------------------------------------------------------------------------
            ! Grad will now be the displacement gradient corresponding to 
            ! the displacement test function
            ! -----------------------------------------------------------------------
            Grad = 0.0d0
            IF (AxialSymmetry) THEN
              SELECT CASE(i)
              CASE (1)
                Grad(1,1) = dBasis(p,1)
                Grad(1,3) = dBasis(p,2)
                Grad(2,2) = 1.0d0/r * Basis(p)
              CASE (2)
                Grad(3,1) = dBasis(p,1)
                Grad(3,3) = dBasis(p,2)                   
              END SELECT
            ELSE
              Grad(i,:) = dBasis(p,:)
            END IF
            !--------------------------------------------------------------------------------------
            ! The following is for handling the part of DS(F)[U], with S the first Piola-Kirchhoff 
            ! stress and U the increment of the deformation gradient. We manipulate the innerproduct
            ! <DS(F)[U],Grad> such that <DS(F)[U],Grad> = <U,W> with W the result of the manipulation. 
            ! First the part of W that do not depend on the response function derivative:
            !--------------------------------------------------------------------------------------
            WorkTensor2 = detDefG * TRACE( MATMUL(Stress,MATMUL(Grad,InvDefG)), dim) * &
                TRANSPOSE(InvDefG) - detDefG * MATMUL(TRANSPOSE(InvDefG), MATMUL(TRANSPOSE(Grad), &
                MATMUL(Stress, TRANSPOSE(InvDefG))))  

            !---------------------------------------------------------------------------------------
            ! The rest of W originating from the response function derivative: 
            !---------------------------------------------------------------------------------------
            WorkTensor1 = detDefG * MATMUL(Grad,InvDefG)
            WorkTensor1 = 0.5d0 * (WorkTensor1 + TRANSPOSE(WorkTensor1))
            WorkVec1(1,1) = WorkTensor1(1,1)
            WorkVec1(2,1) = WorkTensor1(2,2)
            WorkVec1(3,1) = WorkTensor1(3,3)
            WorkVec1(4,1) = WorkTensor1(1,2) +  WorkTensor1(2,1)
            IF (nshr == 3) THEN
              WorkVec1(5,1) = WorkTensor1(1,3) +  WorkTensor1(3,1)
              WorkVec1(6,1) = WorkTensor1(2,3) +  WorkTensor1(3,2)
            END IF
            WorkVec2(1:ntens,1) = MATMUL(TRANSPOSE(StressDer(1:ntens,1:ntens)),WorkVec1(1:ntens,1))
            ! The following is based on the assumed symmetry of StressDer(:,:):
            WorkTensor3 = WorkVec2(1,1)*SymBasis1 + WorkVec2(2,1)*SymBasis2 + &
                WorkVec2(3,1)*SymBasis3 + 2.0d0*WorkVec2(4,1)*SymBasis4
            IF (nshr == 3) THEN
              WorkTensor3 = WorkTensor3 + 2.0d0*WorkVec2(5,1)*SymBasis5 + 2.0d0*WorkVec2(6,1)*SymBasis6
            END IF            
            
            ! -------------------------------------------------------------------------------------
            ! The computation of the differential of the Hencky strain function is based on
            ! its truncated series expansion. 
            ! TO DO: The following involves the differential of the Hencky strain function.
            ! For some reason it doesn't appear to give convergence. Therefore we still omit this and
            ! replace the Hencky strain differential by the differential of the Lagrangian
            ! strain. This is expected to work for reasonably small straining. Find a remedy!
            ! -------------------------------------------------------------------------------------
            !IF (HenckyStrain) THEN
            IF (.FALSE.) THEN
              WorkTensor1 = WorkTensor3
              ! Compute the derivative C'= Dg(WorkTensor1) with g the matrix
              ! square root function:
              WorkTensor3 = MATMUL( TRANSPOSE(QWork), MATMUL(WorkTensor1,QWork) ) 
              WorkVec1(1,1) = 1.0d0/(2.0d0*sqrt(EigenVals(1))) * WorkTensor3(1,1)    
              WorkVec1(2,1) = 1.0d0/(2.0d0*sqrt(EigenVals(2))) * WorkTensor3(2,2)
              WorkVec1(3,1) = 1.0d0/(2.0d0*sqrt(EigenVals(3))) * WorkTensor3(3,3)
              WorkVec1(4,1) = 1.0d0/(sqrt(EigenVals(1)) + sqrt(EigenVals(2))) * &
                  (WorkTensor3(1,2) + WorkTensor3(2,1))
              IF (nshr == 3) THEN
                WorkVec1(5,1) = 1.0d0/(sqrt(EigenVals(1)) + sqrt(EigenVals(3))) * &
                    (WorkTensor3(1,3) + WorkTensor3(3,1))
                WorkVec1(6,1) = 1.0d0/(sqrt(EigenVals(2)) + sqrt(EigenVals(3))) * &
                    (WorkTensor3(2,3) + WorkTensor3(3,2))
              END IF
              WorkTensor3 = WorkVec1(1,1)*SymBasis1 + WorkVec1(2,1)*SymBasis2 + &
                  WorkVec1(3,1)*SymBasis3 + WorkVec1(4,1)*SymBasis4
              IF (nshr == 3) THEN
                WorkTensor3 = WorkTensor3 + WorkVec1(5,1)*SymBasis5 + &
                    WorkVec1(6,1)*SymBasis6
              END IF
              WorkTensor3 = MATMUL( QWork, MATMUL(WorkTensor3,TRANSPOSE(QWork)) )

              WorkTensor3 = 4.0d0 * WorkTensor3 - WorkTensor1
            END IF

            !--------------------------------------------------------------------------
            ! dStress1 is for evaluating the contribution <DS(F^k)[U],Grad> 
            ! in the form <U,W>. Compute W from its splitting:
            !---------------------------------------------------------------------------
            dStress1 = WorkTensor2 + MATMUL(DefG,WorkTensor3)

            IF (AxialSymmetry) THEN
              DO q = 1,ntot
                DO j = 1,cdim
                  SELECT CASE(j)
                  CASE(1)
                    StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        + (dBasis(q,1)*dStress1(1,1) + dBasis(q,2)*dStress1(1,3) &
                        + 1.0d0/r*Basis(q)*dStress1(2,2))*s
                  CASE(2)
                    StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        + (dBasis(q,1)*dStress1(3,1) + dBasis(q,2)*dStress1(3,3) ) * s
                  END SELECT
                END DO
              END DO

              ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                  +(Basis(p)*Force(i)*DetDefG &
                  +Basis(p)*InertialForce(i)*Density)*s

              ExternalForceVector(cdim*(p-1)+i) = ExternalForceVector(cdim*(p-1)+i) &
                  +(Basis(p)*Force(i)*DetDefG &
                  +Basis(p)*InertialForce(i)*Density)*s

              SELECT CASE(i)
              CASE(1)
                ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                    -(dBasis(p,1) * Stress1(1,1) + dBasis(p,2) * Stress1(1,3) &
                    + 1.0d0/r * Basis(p) * Stress1(2,2)) * s 
              CASE(2)
                ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                    -(dBasis(p,1) * Stress1(3,1) + dBasis(p,2) * Stress1(3,3)) * s 
              END SELECT
            ELSE
              DO q = 1,ntot
                DO j = 1,cdim
                  StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                      = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                      + DOT_PRODUCT(dBasis(q,:),dStress1(j,:))*s
                END DO
              END DO

              ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                  +(Basis(p)*Force(i)*DetDefG &
                  +Basis(p)*InertialForce(i)*Density &
                  -DOT_PRODUCT(dBasis(p,:),Stress1(i,:)))*s

              ExternalForceVector(cdim*(p-1)+i) = ExternalForceVector(cdim*(p-1)+i) &
                  +(Basis(p)*Force(i)*DetDefG &
                  +Basis(p)*InertialForce(i)*Density) * s

            END IF

          END DO
        END DO

      END IF STIFFMATRIX_FOR_CHOSEN_STRAIN

      ! Integrate mass matrix:
      ! ----------------------
      DO p = 1,ntot
        DO q = 1,ntot
          DO i = 1,cdim
            MassMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                = MassMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                + Basis(p)*Basis(q)*Density*s
          END DO
        END DO
      END DO

      ! Utilize the Rayleigh damping:
      ! -----------------------------
      DampMatrix = Damping * MassMatrix

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixWithUMAT
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Perform the operation
!
!    A = A + C' * B * C * s
!
! with
!
!    Size( A ) = n x n
!    Size( B ) = m x m
!    Size( C ) = m x n
!------------------------------------------------------------------------------
  SUBROUTINE StrainEnergyDensity(A, B, C, m, n, s)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: A(:,:)
    REAL(KIND=dp), INTENT(IN) :: B(:,:), C(:,:)
    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: s
!------------------------------------------------------------------------------
    A(1:n,1:n) = A(1:n,1:n) + s * MATMUL(TRANSPOSE(C(1:m,1:n)),MATMUL(B(1:m,1:m),C(1:m,1:n))) 
!------------------------------------------------------------------------------
  END SUBROUTINE StrainEnergyDensity
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
! The template for including a material model definition written in the form of
! an Abaqus user subroutine (UMAT). The arguments which can be supposed to be 
! supported by Elmer are capitalized:
!------------------------------------------------------------------------------
  SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
       rpl, ddsddt, drplde, drpldt, STRAN, DSTRAN, TIME, DTIME, TEMP, dTemp, &
       predef, dpred, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, &
       coords, drot, pnewdt, celent, DFRGRD0, DFRGRD1, NOEL, NPT, layer, kspt, &
       kstep, kinc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: STRESS(NTENS)
    ! Requirement for Elmer: At the time of calling the Cauchy stress T_n before
    ! the time/load increment is given
    ! Requirement for umat:  The stress T_{n+1}^{(k)} corresponding to the 
    ! current approximation of the strain increment (DSTRAN) must be returned. 
    ! If the strain increment is defined to be zero in the beginning of the
    ! nonlinear iteration, Elmer will generate a candidate for the strain increment
    ! by assuming purely elastic increment characterized by DDSDDE.

    REAL(KIND=dp), INTENT(INOUT) :: STATEV(NSTATEV)
    ! Requirement for Elmer: The state variables Q_n as specified at the 
    ! previous time/load level for converged solution are given.
    ! Requirement for umat:  The state variables Q_{n+1}^{(k)} corresponding to 
    ! the current approximation of the strain increment must be returned. If 
    ! convergence is attained, these values will be saved and associated with the 
    ! converged solution (cf. the input values)

    REAL(KIND=dp), INTENT(OUT) :: DDSDDE(NTENS,NTENS)
    ! The derivative of (Cauchy) stress response function with respect to the 
    ! strain evaluated for the current approximation must be returned

    REAL(KIND=dp), INTENT(INOUT) :: SSE, SPD, SCD
    ! Requirement for Elmer: Provide specific strain energy (sse), plastic 
    ! dissipation (spd) and creep dissipation (scd) at the previous time/load 
    ! level (these are supposed to be declared to be state variables)
    ! Requirement for umat:  The values of the energy variables corresponding to 
    ! the current approximation may be returned

    REAL(KIND=dp), INTENT(OUT) :: rpl
    ! The mechanical heating power (volumetric)

    REAL(KIND=dp), INTENT(OUT) :: ddsddt(NTENS), drplde(NTENS), drpldt

    REAL(KIND=dp), INTENT(IN) :: STRAN(NTENS)
    ! This gives the strains before the time/load increment.
    ! The strain can be computed from the deformation gradient, so this
    ! argument can be considered to be redundant. Elmer provides
    ! this information anyway. Abaqus assumes that the logarithmic strain 
    ! is used, but Elmer may also use other strain measures.

    REAL(KIND=dp), INTENT(IN) :: DSTRAN(NTENS)
    ! The current candidate for the strain increment to obtain the current 
    ! candidate for the stress. In principle this could be computed from the 
    ! deformation gradient; cf. the variable stran.

    REAL(KIND=dp), INTENT(IN) :: TIME(2)
    ! Both entries give time before the time/load increment (the time for the last
    ! converged solution

    REAL(KIND=dp), INTENT(IN) :: DTIME
    ! The time increment

    REAL(KIND=dp), INTENT(IN) :: TEMP
    ! Temperature before the time/load increment

    REAL(KIND=dp), INTENT(IN) :: dtemp
    ! Temperature increment associated wíth the time/load increment. Currently
    ! Elmer assumes isothermal conditions during the load increment.

    REAL(KIND=dp), INTENT(IN) :: predef(1), dpred(1)
    ! These are just dummy variables for Elmer

    CHARACTER(len=80), INTENT(IN) :: CMNAME
    ! The material model name

    INTEGER, INTENT(IN) :: NDI
    ! The number of direct stress components

    INTEGER, INTENT(IN) :: NSHR
    ! The number of the engineering shear strain components

    INTEGER, INTENT(IN) :: NTENS 
    ! The size of the array containing the stress or strain components

    INTEGER, INTENT(IN) :: NSTATEV
    ! The number of state variables associated with the material model

    REAL(KIND=dp), INTENT(IN) :: PROPS(NPROPS)
    ! An array of material constants

    INTEGER, INTENT(IN) :: NPROPS
    ! The number of the material constants

    REAL(KIND=dp), INTENT(IN) :: coords(3)
    ! The coordinates of the current point could be specified

    REAL(KIND=dp), INTENT(IN) :: drot(3,3)
    ! No support for keeping track of rigid body rotations 
    ! (the variable is initialized to the identity)

    REAL(KIND=dp), INTENT(INOUT) :: pnewdt
    ! Currently, suggesting a new size of time increment does not make any impact

    REAL(KIND=dp), INTENT(IN) :: celent
    ! The element size is not yet provided by Elmer

    REAL(KIND=dp), INTENT(IN) :: DFRGRD0(3,3)
    ! The deformation gradient before the time/load increment (at the previous 
    ! time/load level for converged solution)

    REAL(KIND=dp), INTENT(IN) :: DFRGRD1(3,3)
    ! The deformation gradient corresponding to the current approximation
    ! (cf. the return value of STRESS variable) 

    INTEGER, INTENT(IN) :: NOEL
    ! The element number

    INTEGER, INTENT(IN) :: NPT
    ! The integration point number

    INTEGER, INTENT(IN) :: layer, kspt, kstep, kinc
    ! kstep and kinc could be provided to give information on the incrementation
    ! procedure
!------------------------------------------------------------------------------
    ! Local variables:
    REAL(KIND=dp) :: SymBasis(6,3,3)
    REAL(KIND=dp) :: Identity(3,3), B(3,3), C(3,3), Strain(3,3), S(3,3), Sigma(3,3)
    REAL(KIND=dp) :: WorkMat(3,3)
    REAL(KIND=dp) :: StrainVec(ntens), Stress2(ntens)
    REAL(KIND=dp) :: DetDefG
    REAL(KIND=dp) :: nu, E, LambdaLame, MuLame
    REAL(KIND=dp) :: EigenVals(3), PriWork(102)

    INTEGER :: i, j, k
    INTEGER :: PriLWork=102, PriInfo

    LOGICAL :: LargeDeflection, HenckyStrain
!------------------------------------------------------------------------------

    ! This example creates the basic isotropic model, so there is no state variables
    ! to be handled. 

    ! Use the material model name cmname as a switch:
    HenckyStrain = cmname(1:6)=='hencky'
    IF (cmname(1:2)=='st' .OR. HenckyStrain) THEN
      LargeDeflection = .TRUE.

      SymBasis(1,1:3,1:3) = RESHAPE((/ 1,0,0,0,0,0,0,0,0 /),(/ 3,3 /))
      SymBasis(2,1:3,1:3) = RESHAPE((/ 0,0,0,0,1,0,0,0,0 /),(/ 3,3 /))
      SymBasis(3,1:3,1:3) = RESHAPE((/ 0,0,0,0,0,0,0,0,1 /),(/ 3,3 /)) 
      SymBasis(4,1:3,1:3) = RESHAPE((/ 0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0 /),(/ 3,3 /))
      SymBasis(5,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0 /),(/ 3,3 /))
      SymBasis(6,1:3,1:3) = RESHAPE((/ 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0 /),(/ 3,3 /))
      Identity(1:3,1:3) = RESHAPE((/ 1,0,0,0,1,0,0,0,1 /),(/ 3,3 /))

      B = MATMUL(dfrgrd1, TRANSPOSE(dfrgrd1))
      C = MATMUL(TRANSPOSE(dfrgrd1), dfrgrd1)
      IF (HenckyStrain) THEN
        ! -----------------------------------------------------------
        ! Compute the spectral decomposition of C
        ! -----------------------------------------------------------
        DO i=1,3
          k = i
          DO j=k,3
            WorkMat(i,j) = C(i,j)
          END DO
        END DO
        CALL DSYEV('V', 'U', 3, WorkMat, 3, EigenVals, PriWork, PriLWork, PriInfo)
        IF (PriInfo /= 0) THEN
          CALL Fatal( 'UMAT', 'DSYEV cannot generate eigen basis')          
        END IF

        Strain = 0.0d0
        Strain(1,1) = LOG(SQRT(EigenVals(1)))
        Strain(2,2) = LOG(SQRT(EigenVals(2)))       
        Strain(3,3) = LOG(SQRT(EigenVals(3)))
        ! Transform back to the original coordinates:
        Strain = MATMUL(WorkMat, MATMUL(Strain,TRANSPOSE(WorkMat)))
      ELSE
        ! This example uses the Lagrangian (Green-St Venant) strain tensor:
        Strain = 0.5d0 * (C - Identity)
      END IF
      
      DO i=1,ndi
        StrainVec(i) = Strain(i,i)
      END DO
      DO i=1,nshr
        SELECT CASE(i)
        CASE(1)
          StrainVec(ndi+i) = Strain(1,2)+Strain(2,1)
        CASE(2)
          StrainVec(ndi+i) = Strain(1,3)+Strain(3,1)
        CASE(3)
          StrainVec(ndi+i) = Strain(2,3)+Strain(3,2)
        END SELECT
      END DO

      DetDefG = Dfrgrd1(1,1) * ( Dfrgrd1(2,2)*Dfrgrd1(3,3) - Dfrgrd1(2,3)*Dfrgrd1(3,2) ) + &
              Dfrgrd1(1,2) * ( Dfrgrd1(2,3)*Dfrgrd1(3,1) - Dfrgrd1(2,1)*Dfrgrd1(3,3) ) + &
              Dfrgrd1(1,3) * ( Dfrgrd1(2,1)*Dfrgrd1(3,2) - Dfrgrd1(2,2)*Dfrgrd1(3,1) )
    ELSE
      LargeDeflection = .FALSE.
    END IF
 
    ! Get Young's modulus and the Poisson ratio:
    E = Props(2)
    nu = Props(3)
    
    LambdaLame = E * nu / ( (1.0d0+nu) * (1.0d0-2.0d0*nu) )
    MuLame = E / (2.0d0 * (1.0d0 + nu))
    
    IF (LargeDeflection) THEN
      ! --------------------------------------------------------------------------------
      ! Currently the only nonlinear constitutive law in the umat form is the St. Venant-
      ! Kirchhoff model. In this case we compute the current stress directly by using the 
      ! supplied deformation gradient, so that the strain increment is not used. 
      ! In addition, since it seems that the exact differentiation of the response function 
      ! for the Cauchy stress cannot be done in a straightforward manner, we now make only
      ! a partial approximation. The Cauchy stress is given by
      ! 
      !     sigma(F) = 1/det(F) F S(E(F)) F^T
      !
      ! We however consider only the depedence on the strain as
      !
      !     sigma(.) = 1/det(F) F S(.) F^T
      !
      ! This simplification makes the nonlinear iteration to be an inexact Newton 
      ! method whose performance may deteriorate for large strains. If the convergence is 
      ! attained, the solution nevertheless obeys the St. Venant-Kirchhoff law since
      ! there are no approximations in the computation of the residual.
      ! --------------------------------------------------------------------------------
      ! The constitutive matrix relating the second Piola-Kirchhoff stress and
      ! the strain tensor:
      ddsdde = 0.0d0
      ddsdde(1:ndi,1:ndi) = LambdaLame
      DO i=1,ntens
        ddsdde(i,i) = ddsdde(i,i) + MuLame
      END DO
      DO i=1,ndi
        ddsdde(i,i) = ddsdde(i,i) + MuLame
      END DO
      Stress2 = MATMUL(ddsdde,StrainVec)

      ! The second Piola-Kirchhoff stress in the tensor form:
      S = 0.0d0
      DO i=1,ndi
        S = S + Stress2(i)*SymBasis(i,:,:)
      END DO
      DO i=1,nshr
        S = S + 2.0d0 * Stress2(ndi+i) * SymBasis(ndi+i,:,:)
      END DO

      ! The Cauchy stress tensor:
      Sigma = 1.0d0/DetDefG * MATMUL(dfrgrd1, MATMUL(S,TRANSPOSE(dfrgrd1)))

      DO i=1,ndi
        Stress(i) = Sigma(i,i)
      END DO
      DO i=1,nshr
        SELECT CASE(i)
        CASE(1)
          Stress(ndi+i) = Sigma(1,2)
        CASE(2)
          Stress(ndi+i) = Sigma(1,3)
        CASE(3)
          Stress(ndi+i) = Sigma(2,3)
        END SELECT
      END DO

      ! The derivative: The part corresponding to lambda * tr(E) I
      ddsdde = 0.0d0
      WorkMat = LambdaLame * 1/DetDefG * B
      DO i=1,ndi
        DO j=1,ndi
          ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
        END DO
        DO j=1,nshr
          SELECT CASE(j)
          CASE(1)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
          CASE(2)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
          CASE(3)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
          END SELECT
        END DO
      END DO

      ! The rest corresponding to  2 * mu * E
      DO i=1,ndi
        WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(i,:,:), TRANSPOSE(dfrgrd1)))
        DO j=1,ndi
          ddsdde(j,i) = ddsdde(j,i) + WorkMat(j,j)
        END DO
        DO j=1,nshr
          SELECT CASE(j)
          CASE(1)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,2)
          CASE(2)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(1,3)
          CASE(3)
            ddsdde(ndi+j,i) = ddsdde(ndi+j,i) + WorkMat(2,3)
          END SELECT
        END DO
      END DO

      DO i=1,nshr
        WorkMat = 2.0d0 * MuLame * 1/DetDefG * MATMUL(dfrgrd1, MATMUL(SymBasis(ndi+i,:,:), TRANSPOSE(dfrgrd1)))
        DO j=1,ndi
          ddsdde(j,ndi+i) = ddsdde(j,ndi+i) + 1.0d0 * WorkMat(j,j)
        END DO
        DO j=1,nshr
          SELECT CASE(j)
          CASE(1)
            ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,2)
          CASE(2)
            ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(1,3)
          CASE(3)
            ddsdde(ndi+j,ndi+i) = ddsdde(ndi+j,ndi+i) + 1.0d0 * WorkMat(2,3)
          END SELECT
        END DO
      END DO

    ELSE
      ddsdde = 0.0d0
      ddsdde(1:ndi,1:ndi) = LambdaLame
      DO i=1,ntens
        ddsdde(i,i) = ddsdde(i,i) + MuLame
      END DO
      DO i=1,ndi
        ddsdde(i,i) = ddsdde(i,i) + MuLame
      END DO
      ! We have a linear response function, so the following update is precise
      ! (no higher-order terms related to the notion of differentiability).
      ! Note that we could also define
      !
      !        stress = stress_response_function(stran + dstran)
      ! or
      !        stress = stress_response_function(dfrgrd1)
      !
      ! which may be the precise definition of the functionality required. 
      stress = stress + MATMUL(ddsdde,dstran)
      ! So, for this model, the other way to return the stress:
      !stress = MATMUL(ddsdde,stran+dstran)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE UMAT
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( MassMatrix,DampMatrix,StiffMatrix,ForceVector, &
       LoadVector, InertialLoad, ElasticModulus, NodalPoisson, NodalDensity, NodalDamping, &
       AxialSymmetry,PlaneStress,NodalHeatExpansion, NodalTemperature, Element, n, ntot, &
       Nodes, LocalDisplacement, Isotropic, RotateModuli, TransformMatrix )
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: StiffMatrix(:,:),MassMatrix(:,:),DampMatrix(:,:), &
         NodalHeatExpansion(:,:,:), ElasticModulus(:,:,:)
    REAL(KIND=dp) :: NodalTemperature(:),NodalDensity(:), &
         NodalDamping(:),LoadVector(:,:), InertialLoad(:,:)
    REAL(KIND=dp) :: LocalDisplacement(:,:), TransformMatrix(3,3)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalPoisson

    LOGICAL :: AxialSymmetry,PlaneStress, Isotropic, RotateModuli

    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes

    INTEGER :: n, ntot, N_gauss
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: Basis(ntot)
    REAL(KIND=dp) :: dBasisdx(ntot,3),SqrtElementMetric

    REAL(KIND=dp) :: Force(3), InertialForce(3), NodalLame1(n),NodalLame2(n),Density, &
         Damping,Lame1,Lame2
    REAL(KIND=dp) :: Grad(3,3),Identity(3,3),DetDefG,CofG(3,3),TrueForce(3), G(6,6)
    REAL(KIND=dp) ::  DefG(3,3), Strain(3,3), Stress2(3,3), Stress1(3,3)

    REAL(KIND=dp) :: dDefG(3,3),dStrain(3,3),dStress2(3,3),dStress1(3,3)
    REAL(KIND=dp) :: dDefGU(3,3),dStrainU(3,3),dStress2U(3,3),dStress1U(3,3)

    REAL(KIND=dp) :: Load(3),Temperature, GradBasis(3,3)
    REAL(KIND=dp), DIMENSION(3,3) :: HeatExpansion

    INTEGER :: i,j,k,l,p,q,t,dim,cdim

    REAL(KIND=dp) :: s,u,v,w,r

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    INTEGER :: N_Integ

    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

    LOGICAL :: stat
 !------------------------------------------------------------------------------
    cdim = CoordinateSystemDimension()

    !--------------------------------------------
    ! The dimensionality of the state of stress:
    !---------------------------------------------
    IF (AxialSymmetry) THEN
       dim = 3
    ELSE
       dim = cdim
    END IF

    IF (Isotropic) THEN 
       IF ( PlaneStress ) THEN
          NodalLame1(1:n) = ElasticModulus(1,1,1:n) * NodalPoisson(1:n) /  &
               ( (1.0d0 - NodalPoisson(1:n)**2) )
       ELSE
          NodalLame1(1:n) = ElasticModulus(1,1,1:n) * NodalPoisson(1:n) /  &
               (  (1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)) )
       END IF

       NodalLame2(1:n) = ElasticModulus(1,1,1:n)  / ( 2* (1.0d0 + NodalPoisson(1:n)) )
    END IF


    ForceVector = 0.0D0
    StiffMatrix = 0.0D0
    MassMatrix  = 0.0D0
    DampMatrix  = 0.0d0

    Identity = 0.0D0
    DO i = 1,dim
       Identity(i,i) = 1.0D0
    END DO

    !-------------------------------------------------------
    !    Integration stuff
    !-------------------------------------------------------    
    IntegStuff = GaussPoints( element )

    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

    DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

       !------------------------------------------------------------------------------
       !       Basis function values & derivatives at the integration point
       !------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
            Basis,dBasisdx )

       s = SqrtElementMetric * S_Integ(t)
       IF (AxialSymmetry) THEN
          r = SUM( Basis(1:n) * Nodes % x(1:n) )
          s = s * r
       END IF
       !------------------------------------------------------------------------------
       !       Force at integration point
       !-----------------------------------------------------------------------------   
       Force = 0.0D0
       DO i=1,cdim
          Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
          InertialForce(i) = SUM( InertialLoad(i,1:n)*Basis(1:n) )
       END DO

       IF (Isotropic) THEN
          !-------------------------------------------------
          ! Lame parameters at the integration point
          !------------------------------------------------
          Lame1 = SUM( NodalLame1(1:n)*Basis(1:n) )
          Lame2 = SUM( NodalLame2(1:n)*Basis(1:n) )
          Density = SUM( NodalDensity(1:n)*Basis(1:n) )
          Damping = SUM( NodalDamping(1:n)*Basis(1:n) )

          !------------------------------------------------------------------
          ! Deformation gradient etc. evaluated using the current solution:
          !------------------------------------------------------------------
          Grad = 0.0d0
          IF (AxialSymmetry) THEN
             Grad(1,1) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,1) )
             Grad(1,3) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,2) ) 
             Grad(2,2) = 1.0d0/r * SUM( LocalDisplacement(1,1:ntot) * Basis(1:ntot) )
             Grad(3,1) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,1) )
             Grad(3,3) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,2) )
          ELSE           
             Grad(1:dim,1:dim) = MATMUL(LocalDisplacement(1:dim,1:ntot),dBasisdx(1:ntot,1:dim))
          END IF
          DefG = Identity + Grad
          Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
          Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
          Stress1 = MATMUL(DefG,Stress2)

          SELECT CASE( dim )
          CASE( 1 )
             DetDefG = DefG(1,1)
          CASE( 2 )
             DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
          CASE( 3 )
             DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
                  DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
                  DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
          END SELECT

          !-----------------------------------------------------------------------
          !  Gateaux derivatives of the solution with respect to the displacement:
          !  ---------------------------------------------------------------------
          dDefGU = Grad
          dStrainU = (MATMUL(TRANSPOSE(DefG),dDefGU) &
               + MATMUL(TRANSPOSE(dDefGU),DefG))/2.0D0
          dStress2U = 2.0D0*Lame2*dStrainU + Lame1*TRACE(dStrainU,dim)*Identity
          dStress1U = MATMUL(dDefGU,Stress2) + MATMUL(DefG,dStress2U)

          !----------------------------------------------------------------------------
          ! Loop over the test functions (stiffness matrix for Newton linearization):
          ! ---------------------------------------------------------------------------
          DO p = 1,ntot
             DO i = 1,cdim
                !------------------------------------------------------------------------
                !  Gateaux derivatives of the solution with respect to the test functions:
                ! -----------------------------------------------------------------------
                dDefG = 0.0D0
                IF (AxialSymmetry) THEN
                   SELECT CASE(i)
                   CASE (1)
                      dDefG(1,1) = dBasisdx(p,1)
                      dDefG(1,3) = dBasisdx(p,2)
                      dDefG(2,2) = 1.0d0/r * Basis(p)
                   CASE (2)
                      dDefG(3,1) = dBasisdx(p,1)
                      dDefG(3,3) = dBasisdx(p,2)                   
                   END SELECT
                ELSE                 
                   dDefG(i,:) = dBasisdx(p,:)
                END IF

                dStrain = (MATMUL(TRANSPOSE(DefG),dDefG) &
                     + MATMUL(TRANSPOSE(dDefG),DefG))/2.0D0
                dStress2 = 2.0D0*Lame2*dStrain + Lame1*TRACE(dStrain,dim)*Identity
                dStress1 = MATMUL(dDefG,Stress2) + MATMUL(DefG,dStress2)

                IF (AxialSymmetry) THEN

                   ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                        +(Basis(p)*Force(i)*DetDefG &
                        +Basis(p)*InertialForce(i)*Density &
                        -DDOT_PRODUCT(dDefG,Stress1,dim) &
                        +DDOT_PRODUCT(dDefG,dStress1U,dim))*s

                   DO q = 1,ntot
                      DO j = 1,cdim
                         SELECT CASE(j)
                         CASE(1)
                            StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                                 = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                                 + (dBasisdx(q,1)*dStress1(1,1) + dBasisdx(q,2)*dStress1(1,3) &
                                 + 1.0d0/r*Basis(q)*dStress1(2,2))*s
                         CASE(2)
                            StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                                 = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                                 + (dBasisdx(q,1)*dStress1(3,1) + dBasisdx(q,2)*dStress1(3,3) ) * s
                         END SELECT
                      END DO
                   END DO

                ELSE

                   ForceVector(dim*(p-1)+i) = ForceVector(dim*(p-1)+i) &
                        +(Basis(p)*Force(i)*DetDefG &
                        +Basis(p)*InertialForce(i)*Density &
                        -DOT_PRODUCT(dBasisdx(p,:),Stress1(i,:)) &
                        +DOT_PRODUCT(dBasisdx(p,:),dStress1U(i,:)))*s

                   DO q = 1,ntot
                      DO j = 1,dim
                         StiffMatrix(dim*(p-1)+i,dim*(q-1)+j) &
                              = StiffMatrix(dim*(p-1)+i,dim*(q-1)+j) &
                              + DOT_PRODUCT(dBasisdx(q,:),dStress1(j,:))*s
                      END DO
                   END DO
                END IF
             END DO
          END DO

       ELSE
          ! print *, 'anisotropy active...'
          !--------------------------------------------------------------------------
          ! Anisotropic material is handled in this branch. 
          !-------------------------------------------------------------------------
          IF (dim /= 3 ) &
               CALL Fatal( 'ElasticSolve',  'Material anistropy implemented only for 3-d' )
          IF (AxialSymmetry) &
               CALL Fatal('ElasticSolve', 'Axially symmetric option is not supported for anisotropic materials')

          G = 0.0d0
          DO i=1,SIZE(ElasticModulus,1)
             DO j=1,SIZE(ElasticModulus,2)
                G(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
             END DO
          END DO

          IF ( RotateModuli ) THEN
             CALL RotateElasticityMatrix( G, TransformMatrix, dim )
          END IF

          !-------------------------------------------------------------------------
          ! Compute the formulation variables for the current solution iterate
          !--------------------------------------------------------------------
          Grad = MATMUL(LocalDisplacement(:,1:ntot),dBasisdx)
          DefG = Identity + Grad
          Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0

          SELECT CASE( dim )
          CASE( 1 )
             DetDefG = DefG(1,1)
          CASE( 2 )
             DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
          CASE( 3 )
             DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
                  DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
                  DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
          END SELECT

          !-------------------------------------------------------------
          ! The second Piola-Kirchhoff stress for the current iterate
          !--------------------------------------------------------------
          CALL Strain2Stress(Stress2, Strain, G, dim, .FALSE.)         
          !--------------------------------------------------
          ! The first Piola-Kirchhoff stress
          !--------------------------------------------------
          Stress1 = MATMUL(DefG,Stress2)

          !-----------------------------------------------------------------
          ! dStress2U will be the derivative term Dg(F_k)[grad u_k] with
          ! g the response function giving the second Piola-Kirchhoff stress
          ! in terms of the deformation gradient F
          !------------------------------------------------------------------
          dDefGU = Grad
          dStrainU = (MATMUL(TRANSPOSE(DefG),Grad) &
               + MATMUL(TRANSPOSE(Grad),DefG))/2.0D0
          CALL Strain2Stress(dStress2U, dStrainU, G, dim, .FALSE.)    
          !-------------------------------------------------------------
          ! dStress1U presents the derivative term DS(F_k)[grad u_k] with
          ! S the first  Piola-Kirchhoff stress
          !-------------------------------------------------------------
          dStress1U = MATMUL(Grad,Stress2) + MATMUL(DefG,dStress2U)

          !---------------------------------------------------------
          ! Newton iteration:
          !------------------------------------------------
          DO p = 1,ntot
             DO i = 1,dim
                !------------------------------------------------------------------------
                ! Grad will now be the velocity gradient corresponding to the velocity
                ! test function
                ! -----------------------------------------------------------------------
                Grad = 0.0d0
                Grad(i,:) = dBasisdx(p,:)

                !---------------------------------------------------------------------
                ! dStress2 will correspond to the term (G*)dStrainU, with G* the adjoint
                ! of the elasticity tensor and the strain field dStrainU defined as 
                ! follows: 
                !------------------------------------------------------------------
                dStrainU = (MATMUL(TRANSPOSE(DefG),Grad) &
                     + MATMUL(TRANSPOSE(Grad),DefG))/2.0D0
                CALL Strain2Stress(dStress2, dStrainU, TRANSPOSE(G), dim, .FALSE.)                  

                !-------------------------------------------------------------
                ! Then dStress1 relates to having an equivalent expression for
                ! the derivative DS(F_k)[grad u_{k+1}] with S the first  
                ! Piola-Kirchhoff stress.
                !-------------------------------------------------------------
                dStress1 = MATMUL(Grad,Stress2) + MATMUL(DefG,dStress2)

                ForceVector(dim*(p-1)+i) = ForceVector(dim*(p-1)+i) &
                     +(Basis(p)*Force(i)*DetDefG &
                     +Basis(p)*InertialForce(i)*Density &
                     -DOT_PRODUCT(dBasisdx(p,:),Stress1(i,:)) &
                     +DOT_PRODUCT(dBasisdx(p,:),dStress1U(i,:)))*s

                DO q = 1,ntot
                   DO j = 1,dim
                      StiffMatrix(dim*(p-1)+i,dim*(q-1)+j) &
                           = StiffMatrix(dim*(p-1)+i,dim*(q-1)+j) &
                           + DOT_PRODUCT(dBasisdx(q,:),dStress1(j,:))*s
                   END DO
                END DO
             END DO
          END DO
       END IF


       !      Integrate mass matrix:
       !      ----------------------
       DO p = 1,ntot
          DO q = 1,ntot
             DO i = 1,cdim

                MassMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                     = MassMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                     + Basis(p)*Basis(q)*Density*s

             END DO
          END DO
       END DO

       !      Utilize the Rayleigh damping:
       !      -----------------------------
       DampMatrix = Damping * MassMatrix

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE NeoHookeanLocalMatrix( MassMatrix,DampMatrix,StiffMatrix,ForceVector, &
       LoadVector, InertialLoad, NodalYoung, NodalPoisson, NodalDensity, NodalDamping, &
       AxialSymmetry, PlaneStress, NodalHeatExpansion, NodalTemperature, Element, n, ntot, &
       Nodes, LocalDisplacement, MixedFormulation )
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: StiffMatrix(:,:),MassMatrix(:,:),DampMatrix(:,:), &
         NodalHeatExpansion(:,:,:), NodalYoung(:,:,:)
    REAL(KIND=dp) :: NodalTemperature(:),NodalDensity(:), &
         NodalDamping(:),LoadVector(:,:), InertialLoad(:,:)
    REAL(KIND=dp) :: LocalDisplacement(:,:)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector,NodalPoisson

    LOGICAL :: AxialSymmetry, PlaneStress, MixedFormulation

    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes

    INTEGER :: n, ntot, N_gauss
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    REAL(KIND=dp) :: Basis(ntot)
    REAL(KIND=dp) :: dBasisdx(ntot,3),SqrtElementMetric

    REAL(KIND=dp) :: Force(3), InertialForce(3), NodalLame1(n),NodalLame2(n),Density, &
         Damping,Lame1,Lame2,NodalPressure(ntot),Pressure,NodalPressurePar(n),PressurePar
    REAL(KIND=dp) :: Grad(3,3),InvC(3,3),Identity(3,3),DetDefG,CofG(3,3),TrueForce(3)
    REAL(KIND=dp) :: DefG(3,3), InvDefG(3,3),Strain(3,3), Stress2(3,3), Stress1(3,3)
    REAL(KIND=dp) :: dDefG(3,3),dStrain(3,3),dStress2(3,3),dStress1(3,3)
    REAL(KIND=dp) :: dDefGU(3,3),dStrainU(3,3),dStress2U(3,3),dStress1U(3,3)

    REAL(KIND=dp) :: Load(3),Temperature, GradBasis(3,3)
    REAL(KIND=dp), DIMENSION(3,3) :: HeatExpansion
    REAL(KIND=dp) :: s,u,v,w,r

    INTEGER :: i,j,k,l,p,q,t,dim,cdim,DOFs

    INTEGER :: N_Integ

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    LOGICAL :: stat
!------------------------------------------------------------------------------

    cdim = CoordinateSystemDimension()
    !--------------------------------------------
    ! The dimensionality of the state of stress:
    !---------------------------------------------
    IF (AxialSymmetry) THEN
       dim = 3
    ELSE
       dim = cdim
    END IF

    !------------------------------------------------------------------------------
    ! If the mixed formulation is employed, the auxiliary variable is used to
    ! handle the terms that would grow without a limit as the Poisson ratio approaches
    ! the value 1/2. The code used in the standard case is reused by redefining 
    ! the Lame parameter mu and adding remaining terms afterwards.
    !------------------------------------------------------------------------------
    IF( MixedFormulation ) THEN

      IF (PlaneStress) CALL Warn( 'ElasticSolve',  &
          'Mixed formulation does not support plane stress: plane strain assumed instead' )

      DOFs = cdim + 1
      ! To reuse the code, set the lambda parameter to zero and instead
      ! introduce the epsilon parameter = 1/lambda:
      NodalLame1(1:n) = 0.0d0
      IF ( ALL( ABS(NodalPoisson(1:n)) < AEPS ) ) THEN
        CALL Fatal( 'ElasticSolve',  &
            'Mixed formulation with the zero Poisson ratio is not allowed' )
      ELSE
        NodalPressurePar(1:n) = (1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)) / &
            ( NodalYoung(1,1,1:n) * NodalPoisson(1:n)  )
      END IF
       
    ELSE
      DOFs = cdim

      IF ( PlaneStress ) THEN
        NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) /  &
            ( (1.0d0 - NodalPoisson(1:n)**2) )
      ELSE
        NodalLame1(1:n) = NodalYoung(1,1,1:n) * NodalPoisson(1:n) /  &
            (  (1.0d0 + NodalPoisson(1:n)) * (1.0d0 - 2.0d0*NodalPoisson(1:n)) )
      END IF
    END IF

    NodalLame2(1:n) = NodalYoung(1,1,1:n)  / ( 2* (1.0d0 + NodalPoisson(1:n)) )

    ForceVector = 0.0D0
    StiffMatrix = 0.0D0
    MassMatrix  = 0.0D0
    DampMatrix  = 0.0d0

    Identity = 0.0D0
    DO i = 1,dim
       Identity(i,i) = 1.0D0
    END DO

    IntegStuff = GaussPoints( element )

    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

    DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

       !------------------------------------------------------------------------------
       !     Basis function values & derivatives at the integration point
       !------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric,Basis,dBasisdx )

       s = SqrtElementMetric * S_Integ(t)
       IF (AxialSymmetry) THEN
          r = SUM( Basis(1:n) * Nodes % x(1:n) )
          s = s * r
       END IF

       !------------------------------------------------------------------------
       !     Force at integration point
       !------------------------------------------------------------------------   
       Force = 0.0D0
       DO i=1,cdim
          Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
          InertialForce(i) = SUM( InertialLoad(i,1:n)*Basis(1:n) )
       END DO
       !-----------------------------------------------------------------------
       !     Material properties at the integration point
       !-----------------------------------------------------------------------
       Lame1 = SUM( NodalLame1(1:n)*Basis(1:n) )
       Lame2 = SUM( NodalLame2(1:n)*Basis(1:n) )
       IF (MixedFormulation) THEN
         Pressure = SUM( LocalDisplacement(DOFs,1:n) * Basis(1:n) )   
         Lame2 = Lame2 + Pressure
       END IF

       Density = SUM( NodalDensity(1:n)*Basis(1:n) )
       Damping = SUM( NodalDamping(1:n)*Basis(1:n) )

       !--------------------------------------------------------------------
       ! Compute the formulation variables for the current solution iterate
       !--------------------------------------------------------------------
       Grad = 0.0d0
       IF (AxialSymmetry) THEN
          Grad(1,1) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(1,3) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,2) ) 
          Grad(2,2) = 1.0d0/r * SUM( LocalDisplacement(1,1:ntot) * Basis(1:ntot) )
          Grad(3,1) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(3,3) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,2) )
       ELSE           
          Grad(1:dim,1:dim) = MATMUL(LocalDisplacement(1:dim,1:ntot),dBasisdx(1:ntot,1:dim))
       END IF
       DefG = Identity + Grad

       SELECT CASE( dim )
       CASE( 1 )
          DetDefG = DefG(1,1)
       CASE( 2 )
          DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
       CASE( 3 )
          DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
               DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
               DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
       END SELECT

       InvC = MATMUL( TRANSPOSE(DefG), DefG )
       InvDefG = DefG
       !-------------------------------------------------------------
       !  InvC will now be the inverse of the right Cauchy-Green tensor
       !-------------------------------------------------------------
       CALL InvertMatrix( InvC, dim )
       CALL InvertMatrix( InvDefG, dim )       
       !-------------------------------------------------------------
       ! The second Piola-Kirchhoff stress for the current iterate
       !--------------------------------------------------------------
       Stress2 = Lame1/2.0d0 * (DetDefG - 1.0d0) * (DetDefG + 1.0d0) * InvC + &
            Lame2 * (Identity - InvC)
       !--------------------------------------------------
       ! The first Piola-Kirchhoff stress
       !--------------------------------------------------
       Stress1 = MATMUL(DefG,Stress2)

       !-----------------------------------------------------------------
       ! dStress2U gives the derivative term DG(F_k)[grad u_k] with
       ! G the response function giving the second Piola-Kirchhoff stress
       ! in terms of the deformation gradient F
       !------------------------------------------------------------------
       dStress2U =  Lame1 * DetDefG**2 * TRACE( MATMUL(Grad,InvDefG), dim ) * InvC - &
            Lame1/2.0d0 * (DetDefG - 1.0d0) * (DetDefG + 1.0d0) * &
            MATMUL( InvC, & 
            MATMUL( MATMUL(TRANSPOSE(DefG),Grad) + MATMUL(TRANSPOSE(Grad),DefG), InvC) ) + & 
            Lame2 * MATMUL( InvC, & 
            MATMUL( MATMUL(TRANSPOSE(DefG),Grad) + MATMUL(TRANSPOSE(Grad),DefG), InvC) )   

       !-------------------------------------------------------------
       ! dStress1U presents the derivative term DS(F_k)[grad u_k] with
       ! S the first  Piola-Kirchhoff stress
       !-------------------------------------------------------------
       dStress1U = MATMUL(Grad,Stress2) + MATMUL(DefG,dStress2U)

       !---------------------------------------------------------
       ! Newton iteration:
       !------------------------------------------------
       DO p = 1,ntot
          DO i = 1,cdim
             !------------------------------------------------------------------------
             ! Grad will now be the velocity gradient corresponding to the velocity
             ! test function
             ! -----------------------------------------------------------------------
             Grad = 0.0d0
             IF (AxialSymmetry) THEN
                SELECT CASE(i)
                CASE (1)
                   Grad(1,1) = dBasisdx(p,1)
                   Grad(1,3) = dBasisdx(p,2)
                   Grad(2,2) = 1.0d0/r * Basis(p)
                CASE (2)
                   Grad(3,1) = dBasisdx(p,1)
                   Grad(3,3) = dBasisdx(p,2)                   
                END SELECT
             ELSE
                Grad(i,:) = dBasisdx(p,:)
             END IF

             !-----------------------------------------------------------------
             ! dStress2 gives the derivative term DG(F_k)[grad v] with
             ! G the response function giving the second Piola-Kirchhoff stress
             ! in terms of the deformation gradient F and v the test function
             !------------------------------------------------------------------
             dStress2 = Lame1 * DetDefG**2 * TRACE( MATMUL(Grad,InvDefG), dim ) * InvC - &
                  Lame1/2.0d0 * (DetDefG - 1.0d0) * (DetDefG + 1.0d0) * &
                  MATMUL( InvC, & 
                  MATMUL( MATMUL(TRANSPOSE(DefG),Grad) + MATMUL(TRANSPOSE(Grad),DefG), InvC) ) + & 
                  Lame2 * MATMUL( InvC, & 
                  MATMUL( MATMUL(TRANSPOSE(DefG),Grad) + MATMUL(TRANSPOSE(Grad),DefG), InvC) )  

             !-------------------------------------------------------------
             ! dStress1 is the derivative DS(F_k)[grad v] with
             ! S the first  Piola-Kirchhoff stress      
             !-------------------------------------------------------------
             dStress1 = MATMUL(Grad,Stress2) + MATMUL(DefG,dStress2)

             IF (AxialSymmetry) THEN
                ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                     +(Basis(p)*Force(i)*DetDefG &
                     +Basis(p)*InertialForce(i)*Density &
                     -DDOT_PRODUCT(Grad,Stress1,dim) &
                     +DDOT_PRODUCT(Grad,dStress1U,dim))*s
                
                DO q = 1,ntot
                   DO j = 1,cdim
                      SELECT CASE(j)
                      CASE(1)
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(1,1) + dBasisdx(q,2)*dStress1(1,3) &
                              + 1.0d0/r*Basis(q)*dStress1(2,2))*s
                      CASE(2)
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(3,1) + dBasisdx(q,2)*dStress1(3,3) ) * s
                      END SELECT
                   END DO
                END DO
             ELSE               
                ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                     +(Basis(p)*Force(i)*DetDefG &
                     +Basis(p)*InertialForce(i)*Density &
                     -DOT_PRODUCT(dBasisdx(p,:),Stress1(i,:)) &
                     +DOT_PRODUCT(dBasisdx(p,:),dStress1U(i,:)))*s

                DO q = 1,ntot
                   DO j = 1,dim
                      StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                           = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                           + DOT_PRODUCT(dBasisdx(q,:),dStress1(j,:))*s
                   END DO
                END DO
             END IF
          END DO
       END DO

       !--------------------------------------------
       !      Integrate mass matrix:
       !-------------------------------------------
       DO p = 1,ntot
          DO q = 1,ntot
             DO i = 1,cdim
                MassMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                     = MassMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                     + Basis(p)*Basis(q)*Density*s
             END DO
          END DO
       END DO

       !-------------------------------------
       !      Utilize the Rayleigh damping:
       !-------------------------------------
       DampMatrix = Damping * MassMatrix
       
       !-------------------------------------------------------------------------------
       ! Add remaining terms which relate to having the pressure variable as an unknown: 
       !-------------------------------------------------------------------------------
       IF (MixedFormulation) THEN 

         PressurePar = SUM( NodalPressurePar(1:n)*Basis(1:n) )
         Grad = DefG - Identity

         !-------------------------------------------------------------
         ! The constraint equation to determine the pressure variable:
         !-------------------------------------------------------------
         DO p = 1,n
           ! Use Newton's method:
           ForceVector(DOFs*p) = ForceVector(DOFs*p) - 0.5d0 * (DetDefG**2 - 1.0d0) * Basis(p) * s + &
               DetDefG**2 * TRACE(MATMUL(Grad,InvDefG),dim) * Basis(p) * s

           DO q = 1,ntot
             DO i = 1,cdim
               IF (AxialSymmetry) THEN
                 SELECT CASE(i)
                 CASE(1)
                   StiffMatrix(DOFs*p,DOFs*(q-1)+i) = StiffMatrix(DOFs*p,DOFs*(q-1)+i) + &
                       DetDefG**2 * ( dBasisdx(q,1) * InvDefG(1,1) + dBasisdx(q,2) * InvDefG(3,1) + &
                       Basis(q)/r * InvDefG(2,2) ) *  Basis(p) * s
                 CASE(2)
                   StiffMatrix(DOFs*p,DOFs*(q-1)+i) = StiffMatrix(DOFs*p,DOFs*(q-1)+i) + &
                       DetDefG**2 * ( dBasisdx(q,1) * InvDefG(1,3) + dBasisdx(q,2) * InvDefG(3,3) ) * &
                       Basis(p) * s                  
                 END SELECT
               ELSE
                 ! Use Newton's method:
                 StiffMatrix(DOFs*p,DOFs*(q-1)+i) = StiffMatrix(DOFs*p,DOFs*(q-1)+i) + DetDefG**2 * &
                     SUM( dBasisdx(q,1:cdim) * InvDefG(1:cdim,i) ) *  Basis(p) * s
               END IF
             END DO

             IF (q > n) CYCLE
             StiffMatrix(DOFs*p,DOFs*q) = StiffMatrix(DOFs*p,DOFs*q) + PressurePar * Basis(p) * Basis(q) * s

           END DO
         END DO

         !-------------------------------------------------------------
         ! Modify rows corresponding to the displacements:
         !-------------------------------------------------------------
         DO p = 1,ntot
           DO i = 1,cdim
             !------------------------------------------------------------------------
             ! Grad will now be the velocity gradient corresponding to the velocity
             ! test function
             ! -----------------------------------------------------------------------
             Grad = 0.0d0

             IF (AxialSymmetry) THEN
               SELECT CASE(i)
               CASE (1)
                 Grad(1,1) = dBasisdx(p,1)
                 Grad(1,3) = dBasisdx(p,2)
                 Grad(2,2) = 1.0d0/r * Basis(p)
               CASE (2)
                 Grad(3,1) = dBasisdx(p,1)
                 Grad(3,3) = dBasisdx(p,2)                   
               END SELECT

             ELSE
               Grad(i,:) = dBasisdx(p,:)
             END IF

             ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                 + Pressure * dBasisdx(p,i) * s &
                 - Pressure * DDOT_PRODUCT(TRANSPOSE(InvDefG),Grad,dim) * s

             IF ( AxialSymmetry .AND. (i==1) ) ForceVector(DOFs*(p-1)+i) = &
                 ForceVector(DOFs*(p-1)+i) + Pressure * Basis(p)/r * s

             DO q = 1,ntot
               DO j = 1,cdim
                 IF (AxialSymmetry) THEN
                   SELECT CASE(j)
                   CASE(1)
                     StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         - Pressure * ( dBasisdx(q,1) * Grad(1,1) &
                         + dBasisdx(q,2) * Grad(1,3) &
                         + Basis(q)/r * Grad(2,2) ) * s 
                   CASE(2)
                     StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         - Pressure * ( dBasisdx(q,1) * Grad(3,1) &
                         + dBasisdx(q,2) * Grad(3,3) ) * s 
                   END SELECT
                 ELSE
                   StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                       = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                       - Pressure * DOT_PRODUCT(dBasisdx(q,:),Grad(j,:))*s
                 END IF
               END DO

               IF (q <= n) THEN
                 StiffMatrix(DOFs*(p-1)+i,DOFs*q) &
                     = StiffMatrix(DOFs*(p-1)+i,DOFs*q) - Basis(q) * &
                     DDOT_PRODUCT(TRANSPOSE(InvDefG),Grad,dim) * s 
               END IF

             END DO
           END DO
         END DO
       END IF
     END DO

     IF( MixedFormulation) THEN 
        ! Use just the lowest-order basis for the pressure variable:
        DO p = n+1,ntot
           i = DOFs * p
           ForceVector(i)   = 0.0d0
           StiffMatrix(i,:) = 0.0d0
           StiffMatrix(:,i) = 0.0d0       
           StiffMatrix(i,i) = 1.0d0
        END DO
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE NeoHookeanLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix(BoundaryMatrix,BoundaryVector, &
       LoadVector, NodalSpringCoeff,NormalSpring,NodalAlpha,NodalBeta, &
       LocalDisplacement,Element,n,ntot,Parent,pn,pntot,ParentNodes,Flow,fn, &
       FlowNodes,Velocity,Pressure,NodalViscosity, NodalDensity, &
       CompressibilityDefined, AxialSymmetry, NormalTangential, &
       PseudoTraction, MixedFormulation, LargeDeflection)
!------------------------------------------------------------------------------
    USE Integration
    USE LinearAlgebra
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:)
    REAL(KIND=dp) :: LoadVector(:,:), NodalSpringCoeff(:,:,:)
    LOGICAL :: NormalSpring
    REAL(KIND=dp) :: NodalAlpha(:,:), NodalBeta(:)
    REAL(KIND=dp) :: LocalDisplacement(:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, ntot
    TYPE(Element_t), POINTER :: Parent
    INTEGER :: pn, pntot
    TYPE(Nodes_t) :: ParentNodes
    TYPE(Element_t), POINTER :: Flow
    INTEGER :: fn
    TYPE(Nodes_t) :: FlowNodes
    REAL(KIND=dp) :: Velocity(:,:), Pressure(:), NodalViscosity(:), NodalDensity(:)
    LOGICAL :: CompressibilityDefined, AxialSymmetry, NormalTangential
    LOGICAL :: PseudoTraction, MixedFormulation, LargeDeflection
!----------------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot)
    REAL(KIND=dp) :: dBasisdx(ntot,3),SqrtElementMetric, MetricTerm
    REAL(KIND=dp) :: x(n),y(n),z(n), fx(n), fy(n), fz(n), Density, tm(3)

    REAL(KIND=dp) :: PBasis(pntot)
    REAL(KIND=dp) :: PdBasisdx(pntot,3),PSqrtElementMetric

    REAL(KIND=dp) :: FBasis(fn)
    REAL(KIND=dp) :: FdBasisdx(fn,3),FSqrtElementMetric

    REAL(KIND=dp) :: u,v,w,s,r,ParentU,ParentV,ParentW
    REAL(KIND=dp) :: FlowStress(3,3),Viscosity
    REAL(KIND=dp) :: Force(3),Alpha(3),Beta,Normal(3),RefNormal(3),Identity(3,3)
    REAL(KIND=dp) :: Grad(3,3),DefG(3,3),DetDefG,CofDefG(3,3),ScaleFactor
    REAL(KIND=dp) :: SpringCoeff(3,3)
    REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

    INTEGER :: i,j,t,q,p,dim,N_Integ,DOFs

    LOGICAL :: stat,pstat

    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!----------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    dim = CoordinateSystemDimension()
    IF (MixedFormulation) THEN
      DOFs = dim + 1
    ELSE
      DOFs = dim
    END IF

    Identity = 0.0D0
    DO i = 1,dim
       Identity(i,i) = 1.0D0
    END DO

    BoundaryVector = 0.0D0
    BoundaryMatrix = 0.0D0

    DO i = 1,n
       DO j = 1,pn
          IF( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
             x(i) = Parent % TYPE % NodeU(j)
             y(i) = Parent % TYPE % NodeV(j)
             z(i) = Parent % TYPE % NodeW(j)
             EXIT
          END IF
       END DO
    END DO

    IF ( ASSOCIATED( Flow ) ) THEN
       DO i = 1,n
          DO j = 1,fn
             IF ( Element % NodeIndexes(i) == Flow % NodeIndexes(j) ) THEN
                fx(i) = Flow % TYPE % NodeU(j)
                fy(i) = Flow % TYPE % NodeV(j)
                fz(i) = Flow % TYPE % NodeW(j)
                EXIT
             END IF
          END DO
       END DO
    END IF

    IP = GaussPoints( Element )

    DO t=1,IP % n
       u = IP % U(t)
       v = IP % V(t)
       w = IP % W(t)

       stat = ElementInfo( Element, Nodes, u, v, w, SqrtElementMetric, Basis, dBasisdx )
       
       s = SqrtElementMetric * IP % s(t)
       IF (AxialSymmetry) THEN
          r = SUM( Basis(1:n) * Nodes % x(1:n) )
          s = s * r
       END IF

       !     Calculate the basis functions for the parent element:
       !     -----------------------------------------------------
       ParentU = SUM( Basis(1:n)*x(1:n) )
       ParentV = SUM( Basis(1:n)*y(1:n) )
       ParentW = SUM( Basis(1:n)*z(1:n) )

       Pstat= ElementInfo( Parent,ParentNodes,ParentU,ParentV,ParentW, &
            PSqrtElementMetric,PBasis,PdBasisdx )

       ! Compute the cofactor matrix of the deformation gradient from the previous step:
       ! --------------------------------------------------------------------------------
       IF (LargeDeflection) THEN
         Grad = 0.0d0
         DefG = 0.0d0
         Grad(1:dim,1:dim) = MATMUL(LocalDisplacement(1:dim,1:pntot),PdBasisdx(1:pntot,1:dim))
         DefG = Identity + Grad

         SELECT CASE( dim )
         CASE(1)
           DetDefG = DefG(1,1)
         CASE(2)
           DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
         CASE(3)
           DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
               DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
               DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
         END SELECT
         IF (AxialSymmetry) THEN
           DetDefG = (1.0d0 + SUM(PBasis(1:pntot)*LocalDisplacement(1,1:pntot))/r) * DetDefG
         END IF
         CALL InvertMatrix( DefG, dim )        ! Inverse of the deformation gradient
         CofDefG = 0.0d0
         CofDefG(1:dim,1:dim) = DetDefG*TRANSPOSE( DefG(1:dim,1:dim) )   ! Cofactor of the deformation gradient
       ELSE
         CofDefG = Identity
       END IF


       ! Calculate traction from the flow solution:
       ! ------------------------------------------
       IF ( ASSOCIATED( Flow ) ) THEN
          ParentU = SUM( Basis(1:n)*fx(1:n) )
          ParentV = SUM( Basis(1:n)*fy(1:n) )
          ParentW = SUM( Basis(1:n)*fz(1:n) )

          Pstat = ElementInfo( Flow,FlowNodes,ParentU,ParentV,ParentW, &
               FSqrtElementMetric,FBasis,FdBasisdx )

          Grad = MATMUL( Velocity(:,1:fn),FdBasisdx )
          Density    = SUM( NodalDensity(1:fn) * FBasis )
          Viscosity  = SUM( NodalViscosity(1:fn) * FBasis )

          Viscosity = EffectiveViscosity( Viscosity,Density,Velocity(1,:),Velocity(2,:), &
               Velocity(3,:),FlowElement,FlowNodes,fn,fn,ParentU,ParentV,ParentW,LocalIP=t)
          Viscosity  = SUM( NodalViscosity(1:fn) * FBasis )

          FlowStress = Viscosity * ( Grad + TRANSPOSE(Grad) )

          DO i=1,dim
             FlowStress(i,i) = FlowStress(i,i) - SUM( Pressure(1:fn)*FBasis )
             IF( CompressibilityDefined ) THEN
                FlowStress(i,i) = FlowStress(i,i) - Viscosity * (2.0d0/3.0d0)*TRACE(Grad,dim)
             END IF
          END DO
       END IF

       ! -------------------------------------------------------
       ! Normal vector and its transformation:
       ! -------------------------------------------------------
       RefNormal = NormalVector( Element,Nodes,u,v,.TRUE. )
       Normal = MATMUL(CofDefG,RefNormal)
       ! ----------------------------------------------------------------------------------
       ! The metric term that relates the surface area elements in the deformed and
       ! the reference configuration (this is unrelated to the finite element mapping):
       !  -----------------------------------------------------------------------------
       MetricTerm = SQRT( SUM( Normal(1:dim)*Normal(1:dim) ) ) 
       ! -----------------------------------------------------------------------------------
       ! Note that basically all traction BCs yield nonlinear contributions. Here all
       ! dependencies on the solution are estimated simply by using the previous iterate.
       ! First the surface force that is normal to the deformed surface, i.e.
       ! T(x,t)m(x) = beta * m(x) yielding the condition s = beta * cof(F)n for the
       ! pseudo-traction s corresponding to the first Piola-Kirchhoff stress:
       ! -----------------------------------------------------------------------------------
       Force = SUM( NodalBeta(1:n)*Basis(1:n) ) * Normal
       IF ( ASSOCIATED( Flow ) ) THEN
          Force = Force + MATMUL( FlowStress, Normal )
       END IF
       DO q=1,ntot
          DO i=1,dim
             BoundaryVector((q-1)*DOFs+i) = BoundaryVector((q-1)*DOFs+i) + &
                  s * Basis(q) * Force(i)
          END DO
       END DO

       DO i=1,dim
          Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
       END DO

       IF (PseudoTraction) THEN
          ! ----------------------------------------------------------------------------
          ! The pseudo-traction which corresponds to the first Piola-Kirchhoff stress
          ! and measures the true force per unit undeformed area:
          ! ----------------------------------------------------------------------------
          DO q=1,ntot
             DO i=1,dim
                BoundaryVector((q-1)*DOFs+i) = BoundaryVector((q-1)*DOFs+i) + &
                     s * Basis(q) * Force(i)
             END DO
          END DO
       ELSE
          ! ------------------------------------------------------------------------------------------
          ! The true surface force the material description of which is given componentwise with 
          ! respect to the frame of reference (the metric term arises here as the pseudo-traction
          ! vector corresponding to the first Piola-Kirchhoff stress expresses surface force per unit
          ! area in the reference configuration):
          ! ---------------------------------------------------------------------------------------
          DO q=1,ntot
             DO i=1,dim
                BoundaryVector((q-1)*DOFs+i) = BoundaryVector((q-1)*DOFs+i) + &
                     s * Basis(q) * Force(i) * MetricTerm
             END DO
          END DO
       END IF

       ! ---------------------------------------------------------------------------------------------
       ! Spring terms on the boundary: These contributions are defined with respect the undeformed 
       ! configuration.
       ! -------------------------------------------------------------------------------------------
       IF (NormalSpring) THEN
          SpringCoeff(1,1) = SUM(Basis(1:n)*NodalSpringCoeff(1:n,1,1))
          DO p=1,ntot
             DO i=1,dim 
                DO q=1,ntot
                   DO j=1,dim 
                      BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+j) = BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+j) + &
                           SpringCoeff(1,1) * Basis(q) * RefNormal(j) * Basis(p) * RefNormal(i) * s
                   END DO
                END DO
             END DO
          END DO
       ELSE
          DO i=1,dim
             DO j=1,dim
                SpringCoeff(i,j) = SUM(Basis(1:n)*NodalSpringCoeff(1:n,i,j))
             END DO
          END DO
          ! TO DO: More general spring conditions should be treated here
          DO p=1,ntot
             DO i=1,dim 
                DO q=1,ntot
                   DO j=1,dim 
                      BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+j) = BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+j) + &
                           SpringCoeff(i,j) * Basis(q) * Basis(p) * s
                   END DO
                END DO
             END DO
          END DO
       END IF


       !     NOTE: Currently Alpha parameter is set to be zero so the following has no effect:
       !     ---------------------------------------------------------------------------------
       IF (.FALSE.) THEN
          DO i=1,dim
             Alpha(i) = SUM( NodalAlpha(i,1:n)*Basis(1:n) )
          END DO

          DO p=1,ntot
             DO q=1,ntot
                DO i=1,dim
                   BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+i) =  &
                        BoundaryMatrix((p-1)*DOFs+i,(q-1)*DOFs+i) + &
                        s * Alpha(i) * Basis(q) * Basis(p)
                END DO
             END DO
          END DO
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalBoundaryMatrix
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
  SUBROUTINE GenerateStrainVariable(Displacement, NodalStrain, Perm, CalculateStrains, &
      AxialSymmetry, LargeDeflection)
!---------------------------------------------------------------------------------------------------
!  This subroutine creates the strain field. In the case of plane stress this routine does not
!  however generate the strain component E_33. The principal strains are not computed.
!---------------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: Displacement(:), NodalStrain(:)
    INTEGER, POINTER :: Perm(:)
    LOGICAL :: CalculateStrains, AxialSymmetry, LargeDeflection
!--------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: StSolver
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(ValueList_t), POINTER :: Equation

    LOGICAL :: FirstTime = .TRUE., Found, OptimizeBW, GlobalBubbles, Stat, UseMask   
    LOGICAL :: Factorize,  FoundFactorize, FreeFactorize, FoundFreeFactorize
    INTEGER, POINTER :: Permutation(:), Indices(:)
    INTEGER :: dim, elem, n, nd, i, k, l, p, q, N_Gauss, Ind(9), StrainDim

    REAL(KIND=dp), POINTER :: StrainTemp(:)
    REAL(KIND=dp), ALLOCATABLE :: SForceG(:), LocalDisplacement(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Mass(:,:), Force(:), SForce(:), Basis(:), dBasisdx(:,:)

    REAL(KIND=dp) :: Identity(3,3), Strain(3,3), Grad(3,3)
    REAL(KIND=dp) :: u, v, w, Weight, detJ, r, res

    CHARACTER(LEN=MAX_NAME_LEN) :: eqname

    SAVE FirstTime, StSolver, Permutation, Force, SForceG, StrainTemp, Eqname, Nodes, UseMask
    SAVE StrainDim
 !--------------------------------------------------------------------------------
    IF (.NOT. CalculateStrains) RETURN

    IF (AxialSymmetry) THEN
      dim = CoordinateSystemDimension()
    ELSE
      dim = 3
    END IF

    n = Solver % Mesh % MaxElementDOFs
    ALLOCATE( Indices(n), &
         LocalDisplacement(3,n), &
         Mass(n,n), &
         Force(n), &
         SForce(6*n), &
         Basis(n), &
         dBasisdx(n,3) )

    IF (FirstTime) THEN
       ALLOCATE( StSolver )
       StSolver = Solver

       ALLOCATE( Permutation( SIZE(Solver % Variable % Perm) ) )

       CALL ListSetNameSpace('strain:')

       OptimizeBW = GetLogical( StSolver % Values, 'Optimize Bandwidth', Found )
       IF ( .NOT. Found ) OptimizeBW = .TRUE.
       GlobalBubbles = GetLogical( StSolver % Values, 'Bubbles in Global System', Found )
       IF ( .NOT. Found ) GlobalBubbles = .TRUE.

       IF( ListGetLogicalAnyEquation( Model,'Calculate Strains' ) ) THEN
          UseMask = .TRUE.
          eqname = 'Calculate Strains'
       ELSE
          UseMask = .FALSE.
          eqname = TRIM( ListGetString( StSolver % Values,'Equation') )
       END IF
       StSolver % Matrix => CreateMatrix( Model, Solver, Solver % Mesh, Permutation, &
            1, MATRIX_CRS, OptimizeBW, eqname, GlobalBubbles=GlobalBubbles )

       ALLOCATE( StSolver % Matrix % RHS(StSolver % Matrix % NumberOfRows) )
       StSolver % Matrix % Comm = Solver % Matrix % Comm      

       IF (AxialSymmetry) THEN
          StrainDim = 4
       ELSE
          StrainDim = 6
       END IF
       ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StrainDim) )

       ALLOCATE( StrainTemp(StSolver % Matrix % NumberOfRows) )
       StrainTemp = 0.0d0

       CALL VariableAdd( StSolver % Mesh % Variables, StSolver % Mesh, StSolver, &
            'StrainTemp', 1, StrainTemp, Permutation, Output=.FALSE. )
       StSolver % Variable => VariableGet( StSolver % Mesh % Variables, 'StrainTemp' )

       FirstTime = .FALSE.
    ELSE
       CALL ListSetNameSpace('strain:')
    END IF

    Model % Solver => StSolver
    NodalStrain = 0.0d0
    SForceG = 0.0d0

    IF (AxialSymmetry) THEN
       Ind = (/ 1, 4, 4, 3, 0, 0, 0, 0, 0 /)
    ELSE
       Ind = (/ 1, 4, 6, 4, 2, 5, 6, 5, 3 /)
    END IF

    CALL DefaultInitialize()
    !------------------------------------------------------------------------
    ! Assembly loop 
    !------------------------------------------------------------------------
    DO elem = 1, Solver % NumberOfActiveElements
       Element => GetActiveElement(elem, Solver)
       n  = GetElementNOFNodes()
       nd = GetElementDOFs( Indices )

       CALL GetElementNodes( Nodes )
       CALL GetVectorLocalSolution( LocalDisplacement, USolver=Solver )

       Equation => GetEquation()
       !---------------------------------------
       ! Check if strains wanted for this body:
       ! ---------------------------------------
       IF( UseMask ) THEN
          IF(.NOT. GetLogical( Equation, eqname, Found )) CYCLE
       END IF

       IntegStuff = GaussPoints( element )

       Mass = 0.0d0
       Force = 0.0d0
       SForce = 0.0d0        
       Strain = 0.0d0

       DO t=1,IntegStuff % n
          u = IntegStuff % u(t)
          v = IntegStuff % v(t)
          w = IntegStuff % w(t)
          Weight = IntegStuff % s(t)

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx ) 
          Weight = Weight * detJ

          Grad = MATMUL( LocalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
          IF (AxialSymmetry) THEN
             r = SUM(Basis(1:n) * Nodes % x(1:n))
             Grad(3,3) = 1.0d0/r * SUM(LocalDisplacement(1,1:nd) * Basis(1:nd))
          END IF
          
          Strain = (TRANSPOSE(Grad)+Grad)/2.0D0
          IF (LargeDeflection) Strain = Strain + MATMUL(TRANSPOSE(Grad),Grad)/2.0D0

          DO p=1,nd
             DO q=1,nd
                Mass(p,q) = Mass(p,q) + Weight * Basis(q) * Basis(p)
             END DO

             DO i=1,dim
                DO j=i,dim
                   k = Ind( dim*(i-1)+j )
                   SForce(StrainDim*(p-1)+k) = SForce(StrainDim*(p-1)+k) + Weight * Strain(i,j) * Basis(p)
                END DO
             END DO
             IF (AxialSymmetry) &
                  SForce(StrainDim*(p-1)+2) = SForce(StrainDim*(p-1)+2) + Weight * Strain(3,3) * Basis(p)
          END DO
       END DO

       CALL DefaultUpdateEquations( Mass, Force ) 

       !--------------------------------
       ! Assemble global RHS vectors:
       !--------------------------------
       DO p=1,nd
          l = Permutation(Indices(p))
          DO i=1,StrainDim
             SForceG(StrainDim*(l-1)+i) = SForceG(StrainDim*(l-1)+i) + SForce(StrainDim*(p-1)+i)
          END DO
       END DO
    END DO


    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    CALL Info('ElasticSolve','Calculating strain components',Level=7)

    Factorize = GetLogical( SolverParams, 'Linear System Refactorize', FoundFactorize )
    FreeFactorize = GetLogical( SolverParams, &
         'Linear System Free Factorization', FoundFreeFactorize )

    CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .FALSE. )
    CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', .FALSE. )   

    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .TRUE.)
    n = SIZE(StSolver % Variable % Values)

    DO i=1,StrainDim
       IF (AxialSymmetry) THEN
          SELECT CASE(i)
          CASE(1)
             CALL Info('ElasticSolve','Strain Component 11',Level=3)
          CASE(2)
             CALL Info('ElasticSolve','Strain Component 33',Level=3)
          CASE(3)
             CALL Info('ElasticSolve','Strain Component 22',Level=3)                
          CASE(4)
             CALL Info('ElasticSolve','Strain Component 12',Level=3)              
          END SELECT
       ELSE
          SELECT CASE(i)
          CASE(1)
             CALL Info('ElasticSolve','Strain Component 11',Level=3)
          CASE(2)
             CALL Info('ElasticSolve','Strain Component 22',Level=3)
          CASE(3)
             CALL Info('ElasticSolve','Strain Component 33',Level=3)                
          CASE(4)
             CALL Info('ElasticSolve','Strain Component 12',Level=3)
          CASE(5)
             CALL Info('ElasticSolve','Strain Component 23',Level=3)                
          CASE(6)
             CALL Info('ElasticSolve','Strain Component 13',Level=3)
          END SELECT
       END IF

       StSolver % Matrix % RHS = SForceG(i::StrainDim)
       StSolver % Variable % Values = 0.0d0

       res = DefaultSolve()
       WRITE( Message, '(a,g15.8)') 'Solution Norm:', ComputeNorm(StSolver,n)
       CALL Info( 'GenerateStrainVariable', Message, Level=3 )

       DO l=1,SIZE( Permutation )
          IF ( Permutation(l) <= 0 ) CYCLE
          NodalStrain(StrainDim*(Perm(l)-1)+i) = StSolver % Variable % Values(Permutation(l))
       END DO

    END DO

    IF ( FoundFactorize ) THEN
       CALL ListAddLogical( SolverParams, 'Linear System Refactorize', Factorize )
    ELSE
       CALL ListRemove( SolverParams, 'Linear System Refactorize' )
    END IF

    IF ( .NOT. FoundFreeFactorize ) THEN
       CALL ListRemove( SolverParams, 'Linear System Free Factorization' )
    ELSE
       CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', FreeFactorize )
    END IF

    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .FALSE.)

    DEALLOCATE( Indices, &
         LocalDisplacement, &
         MASS, &
         Force, &
         SForce, &
         Basis, &
         dBasisdx )

    Model % Solver => Solver
    CALL ListSetNameSpace('')

    CALL Info('ElasticSolve','Finished strain postprocessing',Level=7)
!--------------------------------------------------------------------------------
  END SUBROUTINE GenerateStrainVariable
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
  SUBROUTINE GenerateStressVariable( PointwiseStateV, NodalStress, Perm, &
       MaxIntegrationPoints, NStateV, CalculateStress, AxialSymmetry)
!--------------------------------------------------------------------------------
!   This subroutine generates the stress field for material models which
!   depend on a list of state variables
!--------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: PointwiseStateV(:,:)    
    REAL(KIND=dp), POINTER :: NodalStress(:)
    INTEGER, POINTER :: Perm(:) 
    INTEGER :: MaxIntegrationPoints, NStateV
    LOGICAL :: CalculateStress, AxialSymmetry
 !---------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: StSolver
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(ValueList_t), POINTER :: Equation

    LOGICAL :: FirstTime = .TRUE., Found, OptimizeBW, GlobalBubbles, Stat, UseMask   
    LOGICAL :: Factorize,  FoundFactorize, FreeFactorize, FoundFreeFactorize

    INTEGER, POINTER :: Permutation(:), Indices(:)
    INTEGER :: dim, elem, n, nd, i, k, l, p, q, N_Gauss, Ind(6), StressDim

    REAL(KIND=dp), POINTER :: StressTemp(:)
    REAL(KIND=dp), ALLOCATABLE :: SForceG(:)
    REAL(KIND=dp), ALLOCATABLE :: Mass(:,:), Force(:), SForce(:), Basis(:), dBasisdx(:,:)

    REAL(KIND=dp) :: u, v, w, Weight, detJ, res

    CHARACTER(LEN=MAX_NAME_LEN) :: eqname

    SAVE FirstTime, StSolver, Permutation, Force, SForceG, StressTemp, Eqname, Nodes, UseMask
    SAVE StressDim
 !--------------------------------------------------------------------------------------------
    IF (.NOT. CalculateStress) RETURN

    dim = CoordinateSystemDimension()

    n = Solver % Mesh % MaxElementDOFs
    ALLOCATE( Indices(n), &
         Mass(n,n), &
         Force(n), &
         SForce(6*n), &
         Basis(n), &
         dBasisdx(n,3) )

    IF (FirstTime) THEN
       ALLOCATE( StSolver )
       StSolver = Solver

       ALLOCATE( Permutation( SIZE(Solver % Variable % Perm) ) )

       CALL ListSetNameSpace('stress:')

       OptimizeBW = GetLogical( StSolver % Values, 'Optimize Bandwidth', Found )
       IF ( .NOT. Found ) OptimizeBW = .TRUE.
       GlobalBubbles = GetLogical( StSolver % Values, 'Bubbles in Global System', Found )
       IF ( .NOT. Found ) GlobalBubbles = .TRUE.

       IF( ListGetLogicalAnyEquation( Model,'Calculate Stresses' ) ) THEN
          UseMask = .TRUE.
          eqname = 'Calculate Stresses'
       ELSE
          UseMask = .FALSE.
          eqname = TRIM( ListGetString( StSolver % Values,'Equation') )
       END IF
       StSolver % Matrix => CreateMatrix( Model, Solver, Solver % Mesh, Permutation, &
            1, MATRIX_CRS, OptimizeBW, eqname, GlobalBubbles=GlobalBubbles )

       ALLOCATE( StSolver % Matrix % RHS(StSolver % Matrix % NumberOfRows) )
       StSolver % Matrix % Comm = Solver % Matrix % Comm      

       IF (AxialSymmetry) THEN
          StressDim = 4
       ELSE
          StressDim = 6
       END IF
       ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StressDim) )

       ALLOCATE( StressTemp(StSolver % Matrix % NumberOfRows) )
       StressTemp = 0.0d0

       CALL VariableAdd( StSolver % Mesh % Variables, StSolver % Mesh, StSolver, &
            'StressTemp', 1, StressTemp, Permutation, Output=.FALSE. )
       StSolver % Variable => VariableGet( StSolver % Mesh % Variables, 'StressTemp' )

       FirstTime = .FALSE.
    ELSE
       CALL ListSetNameSpace('stress:')
    END IF

    Model % Solver => StSolver
    NodalStress = 0.0d0
    SForceG = 0.0d0

    IF (AxialSymmetry) THEN
       Ind = (/ 4, 5, 6, 8, 7, 9 /)
    ELSE
       Ind = (/ 4, 5, 6, 7, 9, 8 /)
    END IF

    CALL DefaultInitialize()
    !------------------------------------------------------------------------
    ! Assembly loop 
    !------------------------------------------------------------------------
    DO elem = 1, Solver % NumberOfActiveElements
       Element => GetActiveElement(elem, Solver)
       n  = GetElementNOFNodes()
       nd = GetElementDOFs( Indices )
       CALL GetElementNodes( Nodes )

       k = (elem-1)*MaxIntegrationPoints

       Equation => GetEquation()
       !---------------------------------------
       ! Check if stresses wanted for this body:
       ! ---------------------------------------
       IF( UseMask ) THEN
          IF(.NOT. GetLogical( Equation, eqname, Found )) CYCLE
       END IF

       IntegStuff = GaussPoints( element )

       Mass = 0.0d0
       Force = 0.0d0
       SForce = 0.0d0        

       DO t=1,IntegStuff % n
          u = IntegStuff % u(t)
          v = IntegStuff % v(t)
          w = IntegStuff % w(t)
          Weight = IntegStuff % s(t)

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx ) 
          Weight = Weight * detJ

          DO p=1,nd
             DO q=1,nd
                Mass(p,q) = Mass(p,q) + Weight * Basis(q) * Basis(p)
             END DO

             DO i=1,StressDim
                SForce(StressDim*(p-1)+i) = SForce(StressDim*(p-1)+i) + Weight * &
                     PointwiseStateV(k+t,NStateV+Ind(i)) * Basis(p)
             END DO
          END DO
       END DO

       CALL DefaultUpdateEquations( Mass, Force ) 

       !--------------------------------
       ! Assemble global RHS vectors:
       !--------------------------------
       DO p=1,nd
          l = Permutation(Indices(p))
          DO i=1,StressDim
             SForceG(StressDim*(l-1)+i) = SForceG(StressDim*(l-1)+i) + SForce(StressDim*(p-1)+i)
          END DO
       END DO
    END DO

    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    CALL Info('ElasticSolve','Calculating stress components',Level=7)

    Factorize = GetLogical( SolverParams, 'Linear System Refactorize', FoundFactorize )
    FreeFactorize = GetLogical( SolverParams, &
         'Linear System Free Factorization', FoundFreeFactorize )

    CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .FALSE. )
    CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', .FALSE. )   

    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .TRUE.)

    n = SIZE(StSolver % Variable % Values)
    DO i=1,StressDim
       IF (AxialSymmetry) THEN
          SELECT CASE(i)
          CASE(1)
             CALL Info('ElasticSolve','Stress Component 11',Level=3)
          CASE(2)
             CALL Info('ElasticSolve','Stress Component 33',Level=3)
          CASE(3)
             CALL Info('ElasticSolve','Stress Component 22',Level=3)                
          CASE(4)
             CALL Info('ElasticSolve','Stress Component 12',Level=3)              
          END SELECT
       ELSE
          SELECT CASE(i)
          CASE(1)
             CALL Info('ElasticSolve','Stress Component 11',Level=3)
          CASE(2)
             CALL Info('ElasticSolve','Stress Component 22',Level=3)
          CASE(3)
             CALL Info('ElasticSolve','Stress Component 33',Level=3)                
          CASE(4)
             CALL Info('ElasticSolve','Stress Component 12',Level=3)
          CASE(5)
             CALL Info('ElasticSolve','Stress Component 23',Level=3)                
          CASE(6)
             CALL Info('ElasticSolve','Stress Component 13',Level=3)
          END SELECT
       END IF

       StSolver % Matrix % RHS = SForceG(i::StressDim)
       StSolver % Variable % Values = 0.0d0

       res = DefaultSolve()
       WRITE( Message, '(a,g15.8)') 'Solution Norm:', ComputeNorm(StSolver,n)
       CALL Info( 'GenerateStressVariable', Message, Level=3 )

       DO l=1,SIZE( Permutation )
          IF ( Permutation(l) <= 0 ) CYCLE
          NodalStress(StressDim*(Perm(l)-1)+i) = StSolver % Variable % Values(Permutation(l))
       END DO
    END DO

    IF ( FoundFactorize ) THEN
       CALL ListAddLogical( SolverParams, 'Linear System Refactorize', Factorize )
    ELSE
       CALL ListRemove( SolverParams, 'Linear System Refactorize' )
    END IF

    IF ( .NOT. FoundFreeFactorize ) THEN
       CALL ListRemove( SolverParams, 'Linear System Free Factorization' )
    ELSE
       CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', FreeFactorize )
    END IF

    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .FALSE.)

    DEALLOCATE( Indices, &
         MASS, &
         Force, &
         SForce, &
         Basis, &
         dBasisdx )

    Model % Solver => Solver
    CALL ListSetNameSpace('')

    CALL Info('ElasticSolve','Finished stress postprocessing',Level=7)
!----------------------------------------------------------------------------------
  END SUBROUTINE GenerateStressVariable
!----------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------
  SUBROUTINE ComputeStressAndStrain( Displacement, NodalStrain, NodalStress, VonMises, Perm, &
       PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle, AxialSymmetry, &
       NeoHookeanMaterial, CalculateStrains, CalculateStresses, CalcPrincipal, &
       CalcPrincipalAngle, MixedFormulation)
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: Displacement(:), NodalStrain(:), NodalStress(:), VonMises(:), &
         PrincipalStress(:), PrincipalStrain(:), Tresca(:), PrincipalAngle(:) 
    INTEGER, POINTER :: Perm(:)
    LOGICAL :: CalculateStrains, CalculateStresses, CalcPrincipal, CalcPrincipalAngle, &
         NeoHookeanMaterial, AxialSymmetry, MixedFormulation
!--------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: StSolver
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(ValueList_t), POINTER :: Equation

    INTEGER, POINTER :: Permutation(:), Indices(:)

    INTEGER :: dim, cdim, n, nd, elem, i, j, k, l, p, q, t, Ind(9), StrainDim, DOFs

    REAL(KIND=dp), POINTER :: StressTemp(:)
    REAL(KIND=dp), ALLOCATABLE :: ForceG(:), SForceG(:), LocalDisplacement(:,:), &
         Mass(:,:), Force(:), SForce(:), Basis(:), dBasisdx(:,:), &
         NodalLame1(:), NodalLame2(:)

    REAL(KIND=dp) :: Strain(3,3), Stress(3,3), Stress2(3,3), Grad(3,3), DefG(3,3), Identity(3,3), &
         InvC(3,3), InvDefG(3,3), u, v, w, Weight, detJ, res, Lame1, Lame2, nu, DetDefG, G(6,6), r, &
         ScalePar, Pres

    LOGICAL :: FirstTime = .TRUE., Found, OptimizeBW, GlobalBubbles, Stat, &
         Factorize,  FoundFactorize, FreeFactorize, FoundFreeFactorize, PlaneStress, &
         Isotropic, UseMask, LimiterOn, ContactOn, ResidualOn

    CHARACTER(LEN=MAX_NAME_LEN) :: eqname

    SAVE StSolver, Permutation, FirstTime, ForceG, SForceG, Nodes, StressTemp, Eqname, StrainDim, UseMask
    !---------------------------------------------------------------------------------------------
    ! These variables are needed for Principal stress calculation;
    ! they are quite small and allocated even if principal stress calculation
    ! is not requested
    !------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: PriCache(3,3), PriTmp, PriW(3),PriWork(102)
    INTEGER       :: PriN=3, PriLWork=102, PriInfo=0
    REAL(KIND=dp) :: PriAngT1=0, PriAngT2=0, PriAngV(3)=0
 !----------------------------------------------------------------------------------------------
    cdim = CoordinateSystemDimension()

    ! The dimensionality of the stress/strain state:
    IF (AxialSymmetry) THEN
       dim = 3
    ELSE
       dim = cdim
    END IF

    IF (MixedFormulation) THEN
      DOFs = cdim + 1 
    ELSE
      DOFs = cdim
    END IF

    n = Solver % Mesh % MaxElementDOFs
    ALLOCATE( Indices(n), &
         LocalDisplacement(4,n), &
         Mass(n,n), &
         Force(6*n), &
         SForce(6*n), &
         Basis(n), &
         dBasisdx(n,3), &
         NodalLame1(n), &
         NodalLame2(n) )   

    IF (FirstTime) THEN
       ALLOCATE( StSolver )
       StSolver = Solver

       ALLOCATE( Permutation( SIZE(Solver % Variable % Perm) ) )
       ! Permutation = Perm

       CALL ListSetNameSpace('stress:')

       OptimizeBW = GetLogical( StSolver % Values, 'Optimize Bandwidth', Found )
       IF ( .NOT. Found ) OptimizeBW = .TRUE.
       GlobalBubbles = GetLogical( StSolver % Values, 'Bubbles in Global System', Found )
       IF ( .NOT. Found ) GlobalBubbles = .TRUE.

       IF( ListGetLogicalAnyEquation( Model,'Calculate Stresses' ) ) THEN
          UseMask = .TRUE.
          eqname = 'Calculate Stresses'
       ELSE
          UseMask = .FALSE.
          eqname = TRIM( ListGetString( StSolver % Values,'Equation') )
       END IF
       StSolver % Matrix => CreateMatrix( Model, Solver, Solver % Mesh, Permutation, &
            1, MATRIX_CRS, OptimizeBW,eqname, GlobalBubbles=GlobalBubbles )

       ALLOCATE( StSolver % Matrix % RHS(StSolver % Matrix % NumberOfRows) )
       StSolver % Matrix % Comm = Solver % Matrix % Comm      

       IF (AxialSymmetry) THEN
          StrainDim = 4
       ELSE
          StrainDim = 6
       END IF

       IF (CalculateStrains) ALLOCATE( ForceG(StSolver % Matrix % NumberOfRows*StrainDim) )
       IF (CalculateStresses) ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StrainDim) )

       ALLOCATE( StressTemp(StSolver % Matrix % NumberOfRows) )
       StressTemp   = 0.0d0

       CALL VariableAdd( StSolver % Mesh % Variables, StSolver % Mesh, StSolver, &
            'StressTemp', 1, StressTemp, Perm, Output=.FALSE. )
       StSolver % Variable => VariableGet( StSolver % Mesh % Variables, 'StressTemp' )

       FirstTime = .FALSE.
    ELSE
       CALL ListSetNameSpace('stress:')
    END IF

    LimiterOn = ListGetLogical( StSolver % Values,'Apply Limiter', Found ) 
    IF( LimiterOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Apply Limiter',.FALSE.)
    END IF
    ContactOn = ListGetLogical( StSolver % Values,'Apply Contact BCs', Found ) 
    IF( ContactOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Apply Contact BCs',.FALSE.)
    END IF
    ResidualOn = ListGetLogical( StSolver % Values,'Linear System Residual Mode', Found ) 
    IF( ResidualOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Linear System Residual Mode',.FALSE.)
    END IF
    

    Model % Solver => StSolver
    IF (AxialSymmetry) THEN
       Ind = (/ 1, 4, 4, 3, 0, 0, 0, 0, 0 /)
    ELSE
       Ind = (/ 1, 4, 6, 4, 2, 5, 6, 5, 3 /)
    END IF
    IF (CalculateStrains) THEN
       NodalStrain = 0.0d0
       ForceG      = 0.0d0
    END IF
    IF (CalculateStresses) THEN
       NodalStress = 0.0d0
       SForceG      = 0.0d0
    END IF
    CALL DefaultInitialize()

    !------------------------------------------------------------------------
    ! Assembly loop 
    !------------------------------------------------------------------------
    DO elem = 1, Solver % NumberOfActiveElements
       Element => GetActiveElement(elem, Solver)
       n  = GetElementNOFNodes()
       nd = GetElementDOFs( Indices )

       CALL GetElementNodes( Nodes )
       CALL GetVectorLocalSolution( LocalDisplacement, USolver=Solver )

       !-------------------------------------------------------------------
       ! Find material parameters
       !--------------------------------------------------------------------
       Equation => GetEquation()
       Material => GetMaterial()

       ! Check if stresses wanted for this body:
       ! ---------------------------------------
       IF( UseMask ) THEN
          IF(.NOT. GetLogical( Equation, eqname, Found )) THEN
             PRINT *,'not active:',TRIM(eqname)
             CYCLE
          END IF
       END IF

       IF (NeoHookeanMaterial) THEN
          Isotropic = .TRUE.
          ElasticModulus(1,1,1:n) = ListGetReal( Material, &
               'Youngs Modulus', n, Indices, Found )
       ELSE
          CALL InputTensor( ElasticModulus, Isotropic, &
               'Youngs Modulus', Material, n, Indices )        
       END IF

       !------------------------------------------------------------------------------
       ! Check whether the rotation transformation of elastic moduli is necessary...
       !------------------------------------------------------------------------------
       RotateModuli = GetLogical( Material, 'Rotate Elasticity Tensor', Found )
       IF ( RotateModuli .AND. (.NOT. Isotropic) ) THEN
          RotateModuli = .FALSE.
          DO i=1,3
             IF( i == 1 ) THEN
                CALL GetConstRealArray( Material, UWrk, &
                     'Material Coordinates Unit Vector 1', Found, Element )
             ELSE IF( i == 2 ) THEN
                CALL GetConstRealArray( Material, UWrk, &
                     'Material Coordinates Unit Vector 2', Found, Element )
             ELSE                
                CALL GetConstRealArray( Material, UWrk, &
                     'Material Coordinates Unit Vector 3', Found, Element )
             END IF

             IF( Found ) THEN
                UnitNorm = SQRT( SUM( Uwrk(1:3,1)**2 ) )
                IF( UnitNorm < EPSILON( UnitNorm ) ) THEN
                   CALL Fatal('ElasticSolve','Given > Materia Coordinate Unit Vector < too short!')
                END IF
                TransformMatrix(i,1:3) = Uwrk(1:3,1) / UnitNorm  
                RotateModuli = .TRUE.
             ELSE 
                TransformMatrix(i,1:3) = 0.0_dp
                TransformMatrix(i,i) = 1.0_dp
             END IF
          END DO

          IF( .NOT. RotateModuli  ) THEN
             CALL Fatal( 'ElasticSolve', &
                  'No unit vectors found but > Rotate Elasticity Tensor < set True?' )
          END IF
       END IF

       PoissonRatio = 0.0d0
       IF (Isotropic) THEN
          PoissonRatio(1:n) = ListGetReal( Material, 'Poisson Ratio', n, Indices )
          IF (MixedFormulation) THEN
            NodalLame1(1:n) = 0.0d0
            PlaneStress = .FALSE.
          ELSE
            !-----------------------------------------------------------------------------------
            ! In the case of plane stress alter the definition of the Lame (lambda) parameter
            ! so that the plane stress components are directly obtained in terms of the
            ! plane strain components. The strain E_33 can then be expressed as
            ! E_33 = -nu/(1-nu)(E_11 + E_22).
            !-----------------------------------------------------------------------------------
            PlaneStress = GetLogical( Equation, 'Plane Stress', Found )
            IF ( PlaneStress ) THEN
              NodalLame1(1:n) = ElasticModulus(1,1,1:n) * PoissonRatio(1:n) /  &
                  ( (1.0d0 - PoissonRatio(1:n)**2) )
            ELSE
              NodalLame1(1:n) = ElasticModulus(1,1,1:n) * PoissonRatio(1:n) /  &
                  (  (1.0d0 + PoissonRatio(1:n)) * (1.0d0 - 2.0d0*PoissonRatio(1:n)) )
            END IF
          END IF
          NodalLame2(1:n) = ElasticModulus(1,1,1:n)  / ( 2* (1.0d0 + PoissonRatio(1:n)) )
       END IF


       Identity = 0.0D0
       DO i = 1,cdim
          Identity(i,i) = 1.0D0
       END DO
       IF (AxialSymmetry .OR. (Isotropic .AND. (.NOT. PlaneStress))) Identity(3,3) = 1.0D0

       IntegStuff = GaussPoints( Element )      
       Strain = 0.0d0
       Stress = 0.0d0
       Mass = 0.0d0
       Force = 0.0d0      
       SForce = 0.0d0        

       DO t=1,IntegStuff % n
          u = IntegStuff % u(t)
          v = IntegStuff % v(t)
          w = IntegStuff % w(t)
          Weight = IntegStuff % s(t)

          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx ) 
          Weight = Weight * detJ

          IF (Isotropic) THEN
             Lame1 = SUM( NodalLame1(1:n)*Basis(1:n) )
             Lame2 = SUM( NodalLame2(1:n)*Basis(1:n) )
             nu = SUM( PoissonRatio(1:n)*Basis(1:n) )
          ELSE
             G = 0.0d0
             DO i=1,SIZE(ElasticModulus,1)
                DO j=1,SIZE(ElasticModulus,2)
                   G(i,j) = SUM( Basis(1:n) * ElasticModulus(i,j,1:n) )
                END DO
             END DO

             IF ( RotateModuli ) THEN
                CALL RotateElasticityMatrix( G, TransformMatrix, dim )
             END IF
          END IF

          Grad = 0.0d0
          Grad(1:cdim,1:cdim) = MATMUL( LocalDisplacement(1:cdim,1:nd), dBasisdx(1:nd,1:cdim) )
          IF (AxialSymmetry) THEN
             r = SUM(Basis(1:n) * Nodes % x(1:n))
             Grad(3,3) = 1.0d0/r * SUM(LocalDisplacement(1,1:nd) * Basis(1:nd))
          END IF
          DefG = Identity + Grad

          SELECT CASE( dim )
          CASE( 1 )
             DetDefG = DefG(1,1)
          CASE( 2 )
             DetDefG = DefG(1,1)*DefG(2,2) - DefG(1,2)*DefG(2,1)
          CASE( 3 )
             DetDefG = DefG(1,1) * ( DefG(2,2)*DefG(3,3) - DefG(2,3)*DefG(3,2) ) + &
                  DefG(1,2) * ( DefG(2,3)*DefG(3,1) - DefG(2,1)*DefG(3,3) ) + &
                  DefG(1,3) * ( DefG(2,1)*DefG(3,2) - DefG(2,2)*DefG(3,1) )
          END SELECT

          Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
          IF (Isotropic .AND. PlaneStress) &
               Strain(3,3) = -nu/(1.0d0-nu)*(Strain(1,1)+Strain(2,2))

          IF (NeoHookeanMaterial) THEN
             IF (MixedFormulation) THEN
               Pres = -SUM(LocalDisplacement(DOFs,1:n) * Basis(1:n))
             ELSE
               Pres = Lame1/2.0d0 * (DetDefG - 1.0d0) * (DetDefG + 1.0d0)
             END IF
             InvC = MATMUL( TRANSPOSE(DefG), DefG )
             InvDefG = DefG
             !-------------------------------------------------------------
             !  InvC will now be the inverse of the right Cauchy-Green tensor
             !-------------------------------------------------------------
             CALL InvertMatrix( InvC, dim )
             CALL InvertMatrix( InvDefG, dim )       
             !-------------------------------------------------------------
             ! The second Piola-Kirchhoff stress for the current iterate
             !--------------------------------------------------------------
             Stress2 =  Pres * InvC + Lame2 * (Identity - InvC)
          ELSE
             IF (.NOT. Isotropic) THEN
                CALL Strain2Stress(Stress2, Strain, G, dim, .FALSE.) 
             ELSE
                Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
             END IF
          END IF
          Stress =  1.0d0/DetDefG * MATMUL( MATMUL(DefG,Stress2), TRANSPOSE(DefG) )

          DO p=1,nd
             DO q=1,nd
                Mass(p,q) = Mass(p,q) + Weight * Basis(q) * Basis(p)
             END DO

             IF (AxialSymmetry) THEN
                IF (CalculateStrains) THEN
                   DO i=1,2
                      DO j=i,2
                         k = Ind( 2*(i-1)+j )
                         Force(4*(p-1)+k) = Force(4*(p-1)+k) + Weight * Strain(i,j) * Basis(p)
                      END DO
                   END DO
                   Force(4*(p-1)+2) = Force(4*(p-1)+2) + Weight * Strain(3,3) * Basis(p)
                END IF
                IF (CalculateStresses) THEN
                   DO i=1,2
                      DO j=i,2
                         k = Ind( 2*(i-1)+j )
                         SForce(4*(p-1)+k) = SForce(4*(p-1)+k) + Weight * Stress(i,j) * Basis(p)
                      END DO
                   END DO
                   SForce(4*(p-1)+2) = SForce(4*(p-1)+2) + Weight * Stress(3,3) * Basis(p)
                END IF
             ELSE
                IF (CalculateStrains) THEN
                   DO i=1,3
                      DO j=i,3
                         k = Ind( 3*(i-1)+j )
                         Force(6*(p-1)+k) = Force(6*(p-1)+k) + Weight * Strain(i,j) * Basis(p)
                      END DO
                   END DO
                END IF
                IF (CalculateStresses) THEN
                   DO i=1,3
                      DO j=i,3
                         k = Ind( 3*(i-1)+j )
                         SForce(6*(p-1)+k) = SForce(6*(p-1)+k) + Weight * Stress(i,j) * Basis(p)
                      END DO
                   END DO
                END IF
             END IF
          END DO
       END DO

       CALL DefaultUpdateEquations( Mass, Force )

       !--------------------------------
       ! Assemble global RHS vectors:
       !--------------------------------   
       IF (CalculateStrains) THEN
          DO p=1,nd
             l = Permutation(Indices(p))
             DO i=1,StrainDim
                ForceG(StrainDim*(l-1)+i) = ForceG(StrainDim*(l-1)+i) + Force(StrainDim*(p-1)+i)
             END DO
          END DO
       END IF

       IF (CalculateStresses) THEN
          DO p=1,nd
             l = Permutation(Indices(p))
             DO i=1,StrainDim
                SForceG(StrainDim*(l-1)+i) = SForceG(StrainDim*(l-1)+i) + SForce(StrainDim*(p-1)+i)
             END DO
          END DO
       END IF

    END DO

    Factorize = GetLogical( SolverParams, 'Linear System Refactorize', FoundFactorize )
    FreeFactorize = GetLogical( SolverParams, &
         'Linear System Free Factorization', FoundFreeFactorize )

    CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .FALSE. )
    CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', .FALSE. )   
    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .TRUE.)

    n = SIZE(StSolver % Variable % Values)
    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    IF (CalculateStrains) THEN
       CALL Info('ElasticSolve','Calculating strain components',Level=7)
       DO i=1,StrainDim
          IF (AxialSymmetry) THEN
             SELECT CASE(i)
             CASE(1)
                CALL Info('ElasticSolve','Strain Component 11',Level=3)
             CASE(2)
                CALL Info('ElasticSolve','Strain Component 33',Level=3)
             CASE(3)
                CALL Info('ElasticSolve','Strain Component 22',Level=3)                
             CASE(4)
                CALL Info('ElasticSolve','Strain Component 12',Level=3)              
             END SELECT
          ELSE
             SELECT CASE(i)
             CASE(1)
                CALL Info('ElasticSolve','Strain Component 11',Level=3)
             CASE(2)
                CALL Info('ElasticSolve','Strain Component 22',Level=3)
             CASE(3)
                CALL Info('ElasticSolve','Strain Component 33',Level=3)                
             CASE(4)
                CALL Info('ElasticSolve','Strain Component 12',Level=3)
             CASE(5)
                CALL Info('ElasticSolve','Strain Component 23',Level=3)                
             CASE(6)
                CALL Info('ElasticSolve','Strain Component 13',Level=3)
             END SELECT
          END IF

          StSolver % Matrix % RHS = ForceG(i::StrainDim)
          StSolver % Variable % Values = 0.0d0

          res = DefaultSolve()
          WRITE( Message, '(a,g15.8)') 'Solution Norm:', ComputeNorm(StSolver,n)
          CALL Info( 'ComputeStressAndStrain', Message, Level=3 )

          DO l=1,SIZE( Permutation )
             IF ( Permutation(l) <= 0 ) CYCLE
             NodalStrain(StrainDim*(Perm(l)-1)+i) = StSolver % Variable % Values(Permutation(l))
          END DO
       END DO
    END IF

    IF (CalculateStresses) THEN
       CALL Info('ElasticSolve','Calculating stress components',Level=7)
       DO i=1,StrainDim
          IF (AxialSymmetry) THEN
             SELECT CASE(i)
             CASE(1)
                CALL Info('ElasticSolve','Stress Component 11',Level=3)
             CASE(2)
                CALL Info('ElasticSolve','Stress Component 33',Level=3)
             CASE(3)
                CALL Info('ElasticSolve','Stress Component 22',Level=3)                
             CASE(4)
                CALL Info('ElasticSolve','Stress Component 12',Level=3)              
             END SELECT
          ELSE
             SELECT CASE(i)
             CASE(1)
                CALL Info('ElasticSolve','Stress Component 11',Level=3)
             CASE(2)
                CALL Info('ElasticSolve','Stress Component 22',Level=3)
             CASE(3)
                CALL Info('ElasticSolve','Stress Component 33',Level=3)                
             CASE(4)
                CALL Info('ElasticSolve','Stress Component 12',Level=3)
             CASE(5)
                CALL Info('ElasticSolve','Stress Component 23',Level=3)                
             CASE(6)
                CALL Info('ElasticSolve','Stress Component 13',Level=3)
             END SELECT
          END IF

          StSolver % Matrix % RHS = SForceG(i::StrainDim)
          StSolver % Variable % Values = 0.0d0

          res = DefaultSolve()
          WRITE( Message, '(a,g15.8)') 'Solution Norm:', ComputeNorm(StSolver,n)
          CALL Info( 'ComputeStressAndStrain', Message, Level=3 )

          DO l=1,SIZE( Permutation )
             IF ( Permutation(l) <= 0 ) CYCLE
             NodalStress(StrainDim*(Perm(l)-1)+i) = StSolver % Variable % Values(Permutation(l))
          END DO
       END DO


       ! Von Mises stress from the component nodal values:
       ! -------------------------------------------------
       VonMises = 0
       IF (Identity(3,3) < 1.0d0) Identity(3,3) = 1.0d0
       Stress = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE

          IF (AxialSymmetry) THEN
             p = 0
             DO j=1,2
                DO k=1,2
                   p = p + 1
                   q = 4 * (Perm(i)-1) + IND(p)
                   Stress(j,k) = NodalStress(q)
                END DO
             END DO
             q = 4 * (Perm(i)-1) + 2
             Stress(3,3) = NodalStress(q)
          ELSE
             p = 0
             DO j=1,3
                DO k=1,3
                   p = p + 1
                   q = 6 * (Perm(i)-1) + IND(p)
                   Stress(j,k) = NodalStress(q)
                END DO
             END DO
          END IF

          Stress(:,:) = Stress(:,:) - TRACE(Stress(:,:),3) * Identity/3
          DO j=1,3
             DO k=1,3
                VonMises(Perm(i)) = VonMises(Perm(i)) + Stress(j,k)**2
             END DO
          END DO
       END DO
       VonMises = SQRT( 3.0d0 * VonMises / 2.0d0 )
    END IF


    IF ( FoundFactorize ) THEN
       CALL ListAddLogical( SolverParams, 'Linear System Refactorize', Factorize )
    ELSE
       CALL ListRemove( SolverParams, 'Linear System Refactorize' )
    END IF

    IF ( .NOT. FoundFreeFactorize ) THEN
       CALL ListRemove( SolverParams, 'Linear System Free Factorization' )
    ELSE
       CALL ListAddLogical( SolverParams, 'Linear System Free Factorization', FreeFactorize )
    END IF
    CALL ListAddLogical(StSolver % Values, 'Skip Compute Nonlinear Change', .FALSE.)


    !----------------------------------------------
    ! The principal and Tresca stresses:
    !--------------------------------------------------
    IF (CalcPrincipal .AND. CalculateStresses) THEN
       CALL Info('ElasticSolve','Calculating principal stresses',Level=7)
       PriCache = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE       
          IF (AxialSymmetry) THEN
             DO j=1,4
                q = 4 * (Perm(i)-1) + j
                IF (j==4) THEN
                   PriCache(1,3) = NodalStress(q)
                ELSE
                   PriCache(j,j) = NodalStress(q)
                END IF
             END DO
          ELSE          
             DO j=1,6
                q = 6 * (Perm(i)-1) + j
                IF (j>3) THEN
                   SELECT CASE(j)
                   CASE(4)
                      PriCache(1,2) = NodalStress(q)
                   CASE(5)
                      PriCache(2,3) = NodalStress(q)
                   CASE(6)
                      PriCache(1,3) = NodalStress(q)
                   END SELECT
                ELSE
                   PriCache(j,j) = NodalStress(q)
                END IF
             END DO
          END IF

          !-----------------------------------------------------------------------------
          ! Use lapack to solve the eigenvalues (i.e. the principal stresses)
          !-----------------------------------------------------------------------------
          CALL DSYEV( 'V', 'U', 3, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
          IF (PriInfo /= 0) THEN 
             CALL Fatal( 'ElasticSolve', 'DSYEV cannot generate eigen basis')
          END IF

          DO l=1,3
             ! The eigenvalues are returned in the opposite order: 
             PrincipalStress(3 * (Perm(i)-1 )+l) = PriW(4-l)                        
          END DO

          IF (CalcPrincipalAngle) THEN
             DO k=1,3
                PrincipalAngle(9 * (Perm(i)-1) + 3*(k-1) + 1) = ACOS(PriCache(1,4-k))
                PrincipalAngle(9 * (Perm(i)-1) + 3*(k-1) + 2) = ACOS(PriCache(2,4-k))
                PrincipalAngle(9 * (Perm(i)-1) + 3*(k-1) + 3) = ACOS(PriCache(3,4-k))
             END DO
          END IF

          ! Tresca:                        
          Tresca(Perm(i)) = (PrincipalStress(3*(Perm(i)-1) +1) - &
               PrincipalStress(3*(Perm(i)-1) +2))/2
          PriTmp = (PrincipalStress(3*(Perm(i)-1) +2) - &
               PrincipalStress(3*(Perm(i)-1) +3))/2
          IF (PriTmp > Tresca(Perm(i)) ) Tresca(Perm(i)) = PriTmp

          PriTmp = (PrincipalStress(3*(Perm(i)-1) +1) - &
               PrincipalStress(3*(Perm(i)-1) +3))/2
          IF (PriTmp > Tresca(Perm(i)) ) Tresca(Perm(i)) = PriTmp
          
       END DO
    END IF

    IF (CalcPrincipal .AND. CalculateStrains) THEN
       CALL Info('ElasticSolve','Calculating principal strains',Level=7)
       PriCache = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE
          IF (AxialSymmetry) THEN
             DO j=1,4
                q = 4 * (Perm(i)-1) + j
                IF (j==4) THEN
                   PriCache(1,3) = NodalStrain(q)
                ELSE
                   PriCache(j,j) = NodalStrain(q)
                END IF
             END DO
          ELSE
             DO j=1,6
                q = 6 * (Perm(i)-1) + j
                IF (j>3) THEN
                   SELECT CASE(j)
                   CASE(4)
                      PriCache(1,2) = NodalStrain(q)
                   CASE(5)
                      PriCache(2,3) = NodalStrain(q)
                   CASE(6)
                      PriCache(1,3) = NodalStrain(q)
                   END SELECT
                ELSE
                   PriCache(j,j) = NodalStrain(q)
                END IF
             END DO
          END IF

          ! Use lapack to solve eigenvalues:
          CALL DSYEV( 'N', 'U', 3, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
          IF (PriInfo /= 0) THEN 
             CALL Fatal( 'ElasticSolve', 'DSYEV cannot generate eigen basis')
          END IF

          DO l=1,3
             PrincipalStrain(3 * (Perm(i)-1 )+l) = PriW(4-l)
          END DO
       END DO
    END IF

    DEALLOCATE( Indices, &
         LocalDisplacement, &
         MASS, &
         FORCE, &
         SForce, &
         Basis, &
         dBasisdx,&
         NodalLame1, &
         NodalLame2 )  

    Model % Solver => Solver

    CALL ListSetNameSpace('')

    IF( LimiterOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Apply Limiter',.TRUE.)
    END IF
    IF( ContactOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Apply Contact BCs',.TRUE.)
    END IF
    IF( ResidualOn ) THEN
      CALL ListAddLogical( StSolver % Values,'Linear System Residual Mode',.TRUE.)
    END IF

    CALL Info('ElasticSolve','Finished postprocessing',Level=7)
!--------------------------------------------------------------------------------
  END SUBROUTINE ComputeStressAndStrain
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
  SUBROUTINE GenerateStateVariable(PointwiseStateV, MaxIntegrationPoints, NStateV, &
      StateSol, StateDir)
!--------------------------------------------------------------------------------
!   This subroutine generates the field for visualizing the state variables of
!   the umat material model. This assumes that the first (six) state variables 
!   define a tensor variable whose principal components are here solved. Whether
!   calling this subroutine is feasible depends on case.   
!--------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: PointwiseStateV(:,:) 
    INTEGER :: MaxIntegrationPoints, NStateV
    TYPE(Variable_t), POINTER :: StateSol, StateDir
!---------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER :: elem, dofs, dirdofs, i, j, k, t
    INTEGER :: PriLWork=102, PriInfo=0
    REAL(KIND=dp) :: Work(3,3), EigenVals(3), PriWork(102)
!---------------------------------------------------------------------------------
    IF (NStateV < 6) CALL Fatal('GenerateStateVariable', &
        'At least 6 state variables should exist')

    dofs = StateSol % DOFs
    dirdofs = StateDir % DOFs
    DO elem = 1, Solver % NumberOfActiveElements
      Element => GetActiveElement(elem, Solver)
      IntegStuff = GaussPoints(Element)
      k = (elem-1)*MaxIntegrationPoints

      DO t=1,IntegStuff % n
        Work(1,1) = PointwiseStateV(k+t,1)
        Work(2,2) = PointwiseStateV(k+t,2)
        Work(3,3) = PointwiseStateV(k+t,3)
        Work(1,2) = PointwiseStateV(k+t,4)
        Work(1,3) = PointwiseStateV(k+t,5)
        Work(2,3) = PointwiseStateV(k+t,6)

        CALL DSYEV('V', 'U', 3, Work, 3, EigenVals, PriWork, PriLWork, PriInfo)
        IF (PriInfo /= 0) THEN
          CALL Fatal( 'GenerateStateVariable', 'DSYEV cannot generate eigen basis')          
        END IF

        j = StateSol % Perm(Element % ElementIndex) + t
        StateSol % Values(dofs*(j-1)+1) = EigenVals(1)
        StateSol % Values(dofs*(j-1)+2) = EigenVals(2)
        StateSol % Values(dofs*(j-1)+3) = EigenVals(3)

        StateDir % Values(dirdofs*(j-1)+1:dirdofs*(j-1)+3) = Work(1:3,1)
        StateDir % Values(dirdofs*(j-1)+4:dirdofs*(j-1)+6) = Work(1:3,2)
        StateDir % Values(dirdofs*(j-1)+7:dirdofs*(j-1)+9) = Work(1:3,3)

        ! Write stress 11 as the first component:
        !StateSol % Values(dofs*(j-1)+1) = PointwiseStateV(k+t,NStateV+4)
        ! Write stress 22 as the second component:
        !StateSol % Values(dofs*(j-1)+2) = PointwiseStateV(k+t,NStateV+5)
      END DO
    END DO
!--------------------------------------------------------------------------------
  END SUBROUTINE GenerateStateVariable
!--------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION TRACE(A,N) RESULT(B)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I
!------------------------------------------------------------------------------
    B = 0.0D0
    DO i = 1,N
       B = B + A(i,i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION TRACE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION DDOT_PRODUCT(A,B,N) RESULT(C)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B(:,:),C
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I,J
!------------------------------------------------------------------------------
    C = 0.0D0
    DO I = 1,N
       DO J = 1,N
          C = C + A(I,J)*B(I,J)
       END DO
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DDOT_PRODUCT
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ElasticSolver
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
   FUNCTION ElastBoundaryResidual( Model, Edge, Mesh, Quant, Perm, Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
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

     INTEGER :: i,j,k,n,l,t,dim,DOFs,Pn,En
     LOGICAL :: stat, GotIt

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Normal(3), EdgeLength
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area
     REAL(KIND=dp) :: Dir(3)

     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Density
     REAL(KIND=dp) :: Temperature
     REAL(KIND=dp) :: Lame1
     REAL(KIND=dp) :: Lame2
     REAL(KIND=dp) :: Damping
     REAL(KIND=dp) :: HeatExpansionCoeff
     REAL(KIND=dp) :: ReferenceTemperature
     REAL(KIND=dp) :: Identity(3,3), YoungsAverage
     REAL(KIND=dp) :: Grad(3,3), DefG(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:), ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), dEdgeBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:), ExtPressure(:)
     REAL(KIND=dp), ALLOCATABLE :: Force(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalYoungsModulus(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame1(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame2(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDamping(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalHeatExpansionCoeff(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalReferenceTemperature(:)

     LOGICAL :: PlaneStress
     INTEGER :: eq_id
     TYPE(ValueList_t), POINTER :: Material
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Gnorm = 0.0d0

     Identity = 0.0d0
     DO i=1,3
        Identity(i,i) = 1.0d0
     END DO

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

     DOFs = dim
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
!    
!    --------------------------------------------------
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

     ALLOCATE( x(En), y(En), z(En), EdgeBasis(En), dEdgeBasisdx(En,3),         &
       Basis(Pn), dBasisdx(Pn,3), Force(3,En), ExtPressure(En),                &
       NodalDisplacement(3,Pn), NodalYoungsModulus(En), Nodaldensity(En),      &
       NodalTemperature(Pn), NodalLame1(En), NodalLame2(En), NodalDamping(Pn), &
       NodalPoissonRatio(En), NodalHeatExpansionCOeff(En),                     &
       NodalReferenceTemperature(En) )

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
     Indicator     = 0.0d0
     EdgeLength    = 0.0d0
     YoungsAverage = 0.0d0
     ResidualNorm  = 0.0d0

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!        IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                  'Flow Force BC', gotIt ) ) CYCLE

!
!       Logical parameters:
!       -------------------
        eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
             minv=1, maxv=Model % NumberOfEquations )
        PlaneStress = ListGetLogical( Model % Equations(eq_id) % Values,'Plane Stress',GotIt )
!
!       Material parameters:
!       --------------------
        k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                minv=1, maxv=Model % NumberOFMaterials )
        Material => Model % Materials(k) % Values
        NodalYoungsModulus(1:En) = ListGetReal( Material,'Youngs Modulus', &
             En, Edge % NodeIndexes, GotIt )
        NodalPoissonRatio(1:En) = ListGetReal( Material, 'Poisson Ratio', &
             En, Edge % NodeIndexes, GotIt )
        NodalTemperature(1:En) = ListGetReal( Material,'Temperature', &
             En, Edge % NodeIndexes, GotIt )
        NodalReferenceTemperature(1:En) = ListGetReal( Material,'Reference Temperature', &
             En, Edge % NodeIndexes, GotIt )
        NodalDensity(1:En) = ListGetReal( Material,'Density',En,Edge % NodeIndexes, GotIt )
        NodalDamping(1:En) = ListGetReal( Material,'Damping',En,Edge % NodeIndexes, GotIt )
        HeatExpansionCoeff   = 0.0D0
        
        IF ( PlaneStress ) THEN
           NodalLame1(1:En) = NodalYoungsModulus(1:En) * NodalPoissonRatio(1:En) /  &
                ( (1.0d0 - NodalPoissonRatio(1:En)**2) )
        ELSE
           NodalLame1(1:En) = NodalYoungsModulus(1:En) * NodalPoissonRatio(1:En) /  &
                (  (1.0d0 + NodalPoissonRatio(1:En)) * ( 1.0d0 - 2.0d0*NodalPoissonRatio(1:En) ) )
        END IF

        NodalLame2(1:En) = NodalYoungsModulus(1:En)  / ( 2.0d0*(1.0d0 + NodalPoissonRatio(1:En) ) )
!
!       Given traction:
!       ---------------
        Force = 0.0d0

        Force(1,1:En) = ListGetReal( Model % BCs(j) % Values, &
            'Force 1', En, Edge % NodeIndexes, GotIt )

        Force(2,1:En) = ListGetReal( Model % BCs(j) % Values, &
            'Force 2', En, Edge % NodeIndexes, GotIt )

        Force(3,1:En) = ListGetReal( Model % BCs(j) % Values, &
            'Force 3', En, Edge % NodeIndexes, GotIt )

!       Force in normal direction:
!       ---------------------------
        ExtPressure(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'Normal Force', En, Edge % NodeIndexes, GotIt )

!       If dirichlet BC for displacement in any direction given,
!       nullify force in that directon:
!       ------------------------------------------------------------------
        Dir = 1.0d0
        s = ListGetConstReal( Model % BCs(j) % Values, 'Displacement 1', GotIt )
        IF ( GotIt ) Dir(1) = 0

        s = ListGetConstReal( Model % BCs(j) % Values, 'Displacement 2', GotIt )
        IF ( GotIt ) Dir(2) = 0

        s = ListGetConstReal( Model % BCs(j) % Values, 'Displacement 3', GotIt )
        IF ( GotIt ) Dir(3) = 0
!
!       Elementwise nodal solution:
!       ---------------------------
        NodalDisplacement = 0.0d0
        DO k=1,DOFs
           NodalDisplacement(k,1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
        END DO
!
!       Integration:
!       ------------
        EdgeLength    = 0.0d0
        YoungsAverage = 0.0d0
        ResidualNorm  = 0.0d0

        IntegStuff = GaussPoints( Edge )

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
               EdgeBasis, dEdgeBasisdx )

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

           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           u = SUM( EdgeBasis(1:En) * x(1:En) )
           v = SUM( EdgeBasis(1:En) * y(1:En) )
           w = SUM( EdgeBasis(1:En) * z(1:En) )

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
              Basis, dBasisdx )

           Lame1 = SUM( NodalLame1(1:En) * EdgeBasis(1:En) )
           Lame2 = SUM( NodalLame2(1:En) * EdgeBasis(1:En) )
!
!          Stress tensor on the edge:
!          --------------------------
           Grad = MATMUL( NodalDisplacement(:,1:Pn),dBasisdx(1:Pn,:) )
           DefG = Identity + Grad
           Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
           Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
           Stress1 = MATMUL(DefG,Stress2)
!
!          Given force at the integration point:
!          -------------------------------------
           Residual = 0.0d0
           Residual = MATMUL( Force(:,1:En), EdgeBasis(1:En) ) - &
                 SUM( ExtPressure(1:En) * EdgeBasis(1:En) ) * Normal

           Residual = Residual - MATMUL( Stress1, Normal ) * Dir

           EdgeLength   = EdgeLength + s
           ResidualNorm = ResidualNorm + s * SUM( Residual(1:dim) ** 2 )
           YoungsAverage = YoungsAverage + &
                   s * SUM( NodalYoungsModulus(1:En) * EdgeBasis(1:En) )
        END DO
        EXIT
     END DO

     IF ( YoungsAverage > AEPS ) THEN
        YoungsAverage = YoungsAverage / EdgeLength
        Indicator = EdgeLength * ResidualNorm / YoungsAverage
     END IF

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( x, y, z, EdgeBasis, dEdgeBasisdx, Basis, dBasisdx,  &
      Force, ExtPressure, NodalDisplacement, NodalYoungsModulus,     &
      Nodaldensity, NodalTemperature, NodalLame1, NodalLame2, NodalDamping, &
      NodalPoissonRatio, NodalHeatExpansionCOeff, NodalReferenceTemperature )

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION TRACE(A,N) RESULT(B)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I
!------------------------------------------------------------------------------
    B = 0.0D0
    DO i = 1,N
       B = B + A(i,i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION TRACE
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   END FUNCTION ElastBoundaryResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION ElastEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT( Indicator )
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

     INTEGER :: i,j,k,l,n,t,dim,DOFs,En,Pn
     LOGICAL :: stat, GotIt

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Stress(3,3,2), Jump(3), Identity(3,3)
     REAL(KIND=dp) :: Normal(3)
     REAL(KIND=dp) :: Displacement(3)
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Density
     REAL(KIND=dp) :: Temperature
     REAL(KIND=dp) :: Lame1
     REAL(KIND=dp) :: Lame2
     REAL(KIND=dp) :: Damping
     REAL(KIND=dp) :: HeatExpansionCoeff
     REAL(KIND=dp) :: ReferenceTemperature
     REAL(KIND=dp) :: YoungsAverage
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual, ResidualNorm, EdgeLength
     REAL(KIND=dp) :: Grad(3,3), DefG(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)

     LOGICAL :: PlaneStress
     INTEGER :: eq_id
     TYPE(ValueList_t), POINTER :: Material

     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:), ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), Basis(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalYoungsModulus(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame1(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame2(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDamping(:)
     REAL(KIND=dp), ALLOCATABLE :: x(:), y(:), z(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalHeatExpansionCoeff(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalReferenceTemperature(:)

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

     DOFs = dim
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs - 1

     Metric = 0.0d0
     Identity = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
        Identity(i,i) = 1.0d0
     END DO
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

     ALLOCATE( Basis(n), EdgeBasis(En), dBasisdx(n,3), x(En), y(En), z(En),   &
       NodalDisplacement(3,n), NodalYoungsModulus(En), NodalPoissonRatio(En), &
       NodalDensity(en), NodalTemperature(n), NodalLame1(En), NodalLame2(En), &
       NodalDamping(En), NodalHeatExpansionCoeff(En), NodalReferenceTemperature(En) )


!    Integrate square of jump over edge:
!    ------------------------------------
     ResidualNorm  = 0.0d0
     EdgeLength    = 0.0d0
     Indicator     = 0.0d0
     Grad          = 0.0d0
     YoungsAverage = 0.0d0

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

        Stress = 0.0d0
        DO i = 1,2
           IF ( i==1 ) THEN
              Element => Edge % BoundaryInfo % Left
           ELSE
              Element => Edge % BoundaryInfo % Right
           END IF

           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE

           Pn = Element % TYPE % NumberOfNodes
           Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
           Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
           Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)

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

           stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
               Basis, dBasisdx )
!
!          Logical parameters:
!          -------------------
           eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                  minv=1, maxv=Model % NumberOFEquations )

           PlaneStress = ListGetLogical( Model % Equations(eq_id) % Values,'Plane Stress',GotIt )
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOfMaterials )

           Material => Model % Materials(k) % Values

           NodalYoungsModulus(1:En) = ListGetReal( Material,'Youngs Modulus', &
                En, Edge % NodeIndexes, GotIt )
           YoungsModulus = SUM( NodalYoungsModulus(1:En) * EdgeBasis(1:En) )

           NodalPoissonRatio(1:En) = ListGetReal( Material, 'Poisson Ratio', &
                En, Edge % NodeIndexes, GotIt )
           PoissonRatio = SUM( NodalPoissonRatio(1:En) * EdgeBasis(1:En) )

           NodalTemperature(1:En) = ListGetReal( Material,'Temperature', &
                En, Edge % NodeIndexes, GotIt )
           Temperature = SUM( NodalTemperature(1:En) * EdgeBasis(1:En) )

           NodalReferenceTemperature(1:En) = ListGetReal( Material,'Reference Temperature', &
                En, Edge % NodeIndexes, GotIt )
           ReferenceTemperature = SUM( NodalReferenceTemperature(1:En) * EdgeBasis(1:En) )

           NodalDensity(1:En) = ListGetReal( Material,'Density',En,Edge % NodeIndexes, GotIt )
           Density = SUM( NodalDensity(1:En) * EdgeBasis(1:En) )

           NodalDamping(1:En) = ListGetReal( Material,'Damping',En,Edge % NodeIndexes, GotIt )
           Damping = SUM( NodalDamping(1:En) * EdgeBasis(1:En) )

           HeatExpansionCoeff   = 0.0D0

           IF ( PlaneStress ) THEN
              NodalLame1(1:En) = NodalYoungsModulus(1:En) * NodalPoissonRatio(1:En) /  &
                   ( (1.0d0 - NodalPoissonRatio(1:En)**2) )
           ELSE
              NodalLame1(1:En) = NodalYoungsModulus(1:En) * NodalPoissonRatio(1:En) /  &
                   (  (1.0d0 + NodalPoissonRatio(1:En)) * ( 1.0d0 - 2.0d0*NodalPoissonRatio(1:En) ) )
           END IF

           NodalLame2(1:En) = NodalYoungsModulus(1:En)  / ( 2.0d0*(1.0d0 + NodalPoissonRatio(1:En) ) )

           Lame1 = SUM( NodalLame1(1:En) * EdgeBasis(1:En) )
           Lame2 = SUM( NodalLame2(1:En) * EdgeBasis(1:En) )
!
!          Elementwise nodal solution:
!          ---------------------------
           NodalDisplacement = 0.0d0
           DO k=1,DOFs
              NodalDisplacement(k,1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
           END DO
!
!          Stress tensor on the edge:
!          --------------------------
           Grad = MATMUL(NodalDisplacement(:,1:Pn),dBasisdx(1:Pn,:) )
           DefG = Identity + Grad
           Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
           Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
           Stress1 = MATMUL(DefG,Stress2)
           Stress(:,:,i) = Stress1

        END DO

        EdgeLength  = EdgeLength + s
        Jump = MATMUL( ( Stress(:,:,1) - Stress(:,:,2)), Normal )
        ResidualNorm = ResidualNorm + s * SUM( Jump(1:dim) ** 2 )

        YoungsAverage = YoungsAverage + s * YoungsModulus

     END DO

     YoungsAverage = YoungsAverage / EdgeLength
     Indicator = EdgeLength * ResidualNorm / YoungsAverage

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( Basis, EdgeBasis, dBasisdx, x, y, z,   &
       NodalDisplacement, NodalYoungsModulus, NodalPoissonRatio,  &
       NodalDensity, NodalTemperature, NodalLame1, NodalLame2,    &
       NodalDamping, NodalHeatExpansionCoeff, NodalReferenceTemperature )

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION TRACE(A,N) RESULT(B)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I
!------------------------------------------------------------------------------
    B = 0.0D0
    DO i = 1,N
       B = B + A(i,i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION TRACE
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   END FUNCTION ElastEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION ElastInsideResidual( Model, Element,  &
                      Mesh, Quant, Perm, Fnorm ) RESULT( Indicator )
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

     INTEGER :: i,j,k,l,m,n,t,dim,DOFs

     LOGICAL :: stat, GotIt

     TYPE( Variable_t ), POINTER :: Var


     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp) :: Density
     REAL(KIND=dp) :: YoungsModulus
     REAL(KIND=dp) :: PoissonRatio
     REAL(KIND=dp) :: Lame1
     REAL(KIND=dp) :: Lame2
     REAL(KIND=dp) :: Damping
     REAL(KIND=dp) :: HeatExpansionCoeff
     REAL(KIND=dp) :: ReferenceTemperature
     REAL(KIND=dp) :: Displacement(3),Identity(3,3)
     REAL(KIND=dp) :: Grad(3,3), DefG(3,3), Strain(3,3), Stress1(3,3), Stress2(3,3)
     REAL(KIND=dp) :: YoungsAverage, Energy
     REAL(KIND=dp) :: Temperature

     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalYoungsModulus(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalPoissonRatio(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame1(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalLame2(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDamping(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalDisplacement(:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalHeatExpansionCoeff(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalReferenceTemperature(:)
     REAL(KIND=dp), ALLOCATABLE :: Stress(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: NodalTemperature(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalForce(:,:), Veloc(:,:), Accel(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)

     INTEGER :: eq_id

     LOGICAL :: PlaneStress, Transient

     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp), POINTER :: Gravity(:,:)
     REAL(KIND=dp) :: Residual(3), ResidualNorm, Area

     TYPE(ValueList_t), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Fnorm = 0.0d0
     Indicator = 0.0d0

     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

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

     DOFs = dim 
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), NodalDensity(n), &
      NodalYoungsModulus(n), NodalPoissonRatio(n), NodalLame1(n), NodalLame2(n), &
      NodalDamping(n), NodalDisplacement(3,n), NodalHeatExpansionCoeff(n), &
      NodalReferenceTemperature(n), Stress(3,3,n), NodalTemperature(n),    &
      NodalForce(3,n), Veloc(3,n), Accel(3,n), Basis(n), dBasisdx(n,3) )

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)
!
!    Logical parameters:
!    -------------------
     eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
              minv=1, maxv=Model % NumberOfEquations )

     PlaneStress = ListGetLogical( Model % Equations(eq_id) % Values, &
          'Plane Stress',GotIt )
!
!    Material parameters:
!    --------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
             minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     NodalYoungsModulus(1:n) = ListGetReal( Material,'Youngs Modulus', &
          n, Element % NodeIndexes, GotIt )

     NodalPoissonRatio(1:n) = ListGetReal( Material, 'Poisson Ratio', &
          n, Element % NodeIndexes, GotIt )

     NodalTemperature(1:n) = ListGetReal( Material,'Temperature', &
          n, Element % NodeIndexes, GotIt )

     NodalReferenceTemperature(1:n) = ListGetReal( Material,'Reference Temperature', &
          n, Element % NodeIndexes, GotIt )

!
!    Check for time dep.
!    -------------------

     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Transient = .TRUE.
        Var => VariableGet( Model % Variables, 'Displacement', .TRUE. )
        DO i=1,DOFs
           Veloc(i,1:n) = Var % PrevValues(DOFs*(Var % Perm(Element % NodeIndexes)-1)+i,1)
           Accel(i,1:n) = Var % PrevValues(DOFs*(Var % Perm(Element % NodeIndexes)-1)+i,2)
        END DO

        NodalDensity(1:n) = ListGetReal( Material,'Density', &
               n, Element % NodeIndexes, GotIt )

        NodalDamping(1:n) = ListGetReal( Material,'Damping', &
               n, Element % NodeIndexes, GotIt )
     ELSE
        Transient = .FALSE.
     END IF

     HeatExpansionCoeff   = 0.0D0

     IF ( PlaneStress ) THEN
        NodalLame1(1:n) = NodalYoungsModulus(1:n) * NodalPoissonRatio(1:n) /  &
             ( (1.0d0 - NodalPoissonRatio(1:n)**2) )
     ELSE
        NodalLame1(1:n) = NodalYoungsModulus(1:n) * NodalPoissonRatio(1:n) /  &
             (  (1.0d0 + NodalPoissonRatio(1:n)) * ( 1.0d0 - 2.0d0*NodalPoissonRatio(1:n) ) )
     END IF

     NodalLame2(1:n) = NodalYoungsModulus(1:n)  / ( 2.0d0*(1.0d0 + NodalPoissonRatio(1:n) ) )
!
!    Elementwise nodal solution:
!    ---------------------------
     NodalDisplacement = 0.0d0
     DO k=1,DOFs
        NodalDisplacement(k,1:n) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
     END DO
!
!    Body Forces:
!    ------------
     k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt, &
                    1, Model % NumberOfBodyForces )

     NodalForce = 0.0d0

     IF ( GotIt .AND. k > 0  ) THEN

        NodalForce(1,1:n) = NodalForce(1,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Stress BodyForce 1', &
             n, Element % NodeIndexes, GotIt )
        
        NodalForce(2,1:n) = NodalForce(2,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Stress BodyForce 2', &
             n, Element % NodeIndexes, GotIt )
        
        NodalForce(3,1:n) = NodalForce(3,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Stress BodyForce 3', &
             n, Element % NodeIndexes, GotIt )

     END IF

     Identity = 0.0D0
     DO i = 1,DIM
        Identity(i,i) = 1.0D0
     END DO
!
!    Values of the stress tensor at node points:
!    -------------------------------------------
     Grad = 0.0d0
     DO i = 1,n
        u = Element % TYPE % NodeU(i)
        v = Element % TYPE % NodeV(i)
        w = Element % TYPE % NodeW(i)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx )

        Lame1 = NodalLame1(i)
        Lame2 = NodalLame2(i)

        Grad = 0.0d0
        Grad = MATMUL(NodalDisplacement(:,1:N),dBasisdx(1:N,:) )
        DefG = Identity + Grad
        Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
        Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
        Stress1 = MATMUL(DefG,Stress2)
        Stress(:,:,i) = Stress1

     END DO
!
!    Integrate square of residual over element:
!    ------------------------------------------
     ResidualNorm = 0.0d0
     Fnorm = 0.0d0
     Area = 0.0d0
     Energy = 0.0d0
     YoungsAverage = 0.0d0

     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF
!
!       Residual of the diff.equation:
!       ------------------------------
        Residual = 0.0d0
        DO i = 1,Dim
           Residual(i) = SUM( NodalForce(i,1:n) * Basis(1:n) )

           IF ( Transient ) THEN
              Residual(i) = Residual(i) + SUM( NodalDensity(1:n) * Basis(1:n) ) * &
                                 SUM( Accel(i,1:n) * Basis(1:n) )

              Residual(i) = Residual(i) + SUM( NodalDamping(1:n) * Basis(1:n) ) * &
                                 SUM( Veloc(i,1:n) * Basis(1:n) )
           END IF

           DO j = 1,Dim
              DO k = 1,n
                 Residual(i) = Residual(i) + Stress(i,j,k)*dBasisdx(k,j)
              END DO
           END DO
        END DO
!
!       Dual norm of the load:
!       ----------------------
        DO i = 1,Dim
           Fnorm = Fnorm + s * SUM( NodalForce(i,1:n) * Basis(1:n) ) ** 2
        END DO

        YoungsAverage = YoungsAverage + s * SUM( NodalYoungsModulus(1:n) * Basis(1:n) )

!       Energy:
!       -------
        Grad = 0.0d0
        Grad = MATMUL(NodalDisplacement(:,1:N),dBasisdx(1:N,:) )
        DefG = Identity + Grad
        Strain = (TRANSPOSE(Grad)+Grad+MATMUL(TRANSPOSE(Grad),Grad))/2.0D0
        Stress2 = 2.0D0*Lame2*Strain + Lame1*TRACE(Strain,dim)*Identity
        Stress1 = MATMUL(DefG,Stress2)
        Energy = Energy + s*DDOTPROD(Strain,Stress1,Dim)/2.0d0

        Area = Area + s
        ResidualNorm = ResidualNorm + s * SUM( Residual(1:dim) ** 2 )

     END DO

     YoungsAverage = YoungsAverage / Area
     Fnorm = Energy
     Indicator = Area * ResidualNorm / YoungsAverage

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, NodalDensity,      &
      NodalYoungsModulus, NodalPoissonRatio, NodalLame1, NodalLame2, &
      NodalDamping, NodalDisplacement, NodalHeatExpansionCoeff,      &
      NodalReferenceTemperature, Stress, NodalTemperature,           &
      NodalForce, Veloc, Accel, Basis, dBasisdx )

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION TRACE(A,N) RESULT(B)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I
!------------------------------------------------------------------------------
    B = 0.0D0
    DO i = 1,N
       B = B + A(i,i)
    END DO
!------------------------------------------------------------------------------
  END FUNCTION TRACE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION DDOTPROD(A,B,N) RESULT(C)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    DOUBLE PRECISION :: A(:,:),B(:,:),C
    INTEGER :: N
!------------------------------------------------------------------------------
    INTEGER :: I,J
!------------------------------------------------------------------------------
    C = 0.0D0
    DO i = 1,N
       DO j = 1,N
          C = C + A(i,j)*B(i,j)
       END DO
    END DO
!------------------------------------------------------------------------------
  END FUNCTION DDOTPROD
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END FUNCTION ElastInsideResidual
!------------------------------------------------------------------------------
