!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
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
  
  CALL ListAddLogical( SolverParams,'Solid Solver',.TRUE.)
  
!------------------------------------------------------------------------------
END SUBROUTINE ElasticSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE ElasticSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE StressLocal, ONLY: StressFieldDefinition, SymTensorComponents
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: dim, i, DOFs
  LOGICAL :: Found, AxialSymmetry, MixedFormulation
  LOGICAL :: MaxwellMaterial
  LOGICAL :: CalculateStrains, CalculateStresses
  LOGICAL :: CalcPrincipalAngle, CalcPrincipal
  LOGICAL :: CalcPrincipalStress, CalcPrincipalStrain
  LOGICAL :: OutputStateVars
  INTEGER :: NState
  TYPE(ValueList_t), POINTER :: Material
  CHARACTER(LEN=MAX_NAME_LEN) :: str
  CHARACTER(*), PARAMETER :: Caller = 'ElasticSolver_init'
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  ! Cylindric Symmetric means axisymmetric with a swirling component, which this
  ! solver has no degree of freedom for -- there is no u_theta in a 2D run. It is
  ! accepted as a synonym because sifs use it as one: fem/tests/CylComAxi is a pure
  ! axisymmetric benchmark, by its own header, that declares Cylindric Symmetric.
  AxialSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  MixedFormulation = GetLogical(SolverParams, 'Mixed Formulation', Found) .AND. &
      GetLogical(SolverParams, 'Neo-Hookean Material', Found)

  dim = CoordinateSystemDimension()

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

  !----------------------------------------------------------------------------
  ! MAXWELL VISCOELASTICITY, set up exactly as StressSolve sets it up, because a
  ! sif written for there has to mean the same thing here. Four of these are not
  ! defaults but impositions, and they are what makes the scheme the scheme:
  !
  ! - the law relates a stress RATE to a strain rate, so the time integration is
  !   first order and single step: BDF of order one. This is also why the keyword
  !   above had to start honouring "Time derivative order" before this could work;
  ! - the lag stress lives at the integration points, in an -ip variable whose
  !   components are the independent ones of a symmetric tensor. Under axial
  !   symmetry that is what makes room for the hoop component at all;
  ! - at least two nonlinear iterations, since the lag stress is updated during the
  !   assembly and the first pass therefore solves with a stress from the previous
  !   step rather than this one;
  ! - and the axisymmetric incompressible combination is refused, as it is there:
  !   the pressure would reach the hoop equation both through the lag stress and
  !   through the mixed formulation.
  !
  ! "Maxwell material" is accepted in the solver section as well as in a material,
  ! and copied onto every material when it is given there -- again as StressSolve
  ! does, since the assembly reads it per material.
  !----------------------------------------------------------------------------
  MaxwellMaterial = ListGetLogicalAnyMaterial( Model, 'Maxwell material' )
  IF ( .NOT. MaxwellMaterial ) THEN
    MaxwellMaterial = GetLogical( SolverParams, 'Maxwell material', Found )
    IF ( MaxwellMaterial ) THEN
      DO i=1,Model % NumberOfMaterials
        CALL ListAddLogical( Model % Materials(i) % Values, 'Maxwell material', .TRUE. )
      END DO
    END IF
  END IF

  IF ( MaxwellMaterial ) THEN
    CALL ListAddString( SolverParams, 'Timestepping Method', 'BDF' )
    CALL ListAddInteger( SolverParams, 'BDF Order', 1 )
    CALL ListAddInteger( SolverParams, 'Time derivative Order', 1 )

    IF ( AxialSymmetry .AND. GetLogical( SolverParams, 'Incompressible', Found ) ) &
        CALL Fatal( Caller, 'Maxwell material with "Incompressible" is not '// &
        'supported in axisymmetric coordinates' )

    CALL ListAddString( SolverParams, &
        NextFreeKeyword('Exported Variable ',SolverParams), &
        '-dofs '//I2S(SymTensorComponents(dim,AxialSymmetry))//' -ip ve_stress' )

    i = GetInteger( SolverParams, 'Nonlinear System Min Iterations', Found )
    CALL ListAddInteger( SolverParams, 'Nonlinear System Min Iterations', MAX(i,2) )
    i = GetInteger( SolverParams, 'Nonlinear System Max Iterations', Found )
    CALL ListAddInteger( SolverParams, 'Nonlinear System Max Iterations', MAX(i,2) )
  END IF

  ! Second order in time is this solver's default, but let a sif ask for the first
  ! order transient StressSolve has always offered: it is the same keyword, and
  ! ListAddInteger would overwrite the value the user gave. The core reads it too
  ! (Solver % TimeOrder), so the assembly and the time history have to agree.
  ! ListAddNew, so the order one that a Maxwell material just imposed survives.
  CALL ListAddNewInteger( SolverParams,'Time derivative order', 2 )
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
     ! One layout for every 2D case, axisymmetric or not, and the same one
     ! StressSolve writes: (11,22,33,12) with the 23 and 13 shears dropped, those
     ! being identically zero in two dimensions. The out-of-plane 33 is kept, and
     ! is the hoop component under axial symmetry.
     CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          TRIM(StressFieldDefinition('Stress',dim)) )

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
     CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          TRIM(StressFieldDefinition('Strain',dim)) )

     IF (CalcPrincipalStrain) THEN
        CALL ListAddString( SolverParams,&
             NextFreeKeyword('Exported Variable ',SolverParams), &
             'Principal Strain[Principal Strain:3]' )
             
     END IF
  END IF

  IF (.NOT. ListCheckPresentAnyMaterial(Model, 'UMAT Subroutine') ) RETURN

  
  ! Following definitions only apply to UMAT

  OutputStateVars = GetLogical(SolverParams, 'Output State Variables', Found)

  IF ( dim == 3 ) THEN
    IF (OutputStateVars) THEN
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          '-ip UmatStress[UmatStress_xx:1 UmatStress_yy:1 UmatStress_zz:1 UmatStress_xy:1 UmatStress_yz:1 UmatStress_xz:1]' )
    ELSE
      str = 'UmatStress[UmatStress_xx:1 UmatStress_yy:1 UmatStress_zz:1 UmatStress_xy:1 UmatStress_yz:1 UmatStress_xz:1]'
      str = '-nooutput -ip '//TRIM(str)
      CALL ListAddString( SolverParams, NextFreeKeyword('Exported Variable ',SolverParams), str)
    END IF
  ELSE
    IF (OutputStateVars) THEN
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          '-ip UmatStress[UmatStress_xx:1 UmatStress_zz:1 UmatStress_yy:1 UmatStress_xy:1]' )
    ELSE
      CALL ListAddString( SolverParams,&
          NextFreeKeyword('Exported Variable ',SolverParams), &
          '-nooutput -ip UmatStress[UmatStress_xx:1 UmatStress_zz:1 UmatStress_yy:1 UmatStress_xy:1]' )
    END IF
  END IF

  IF (OutputStateVars) THEN
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable ',SolverParams), &
        '-dofs 3 -ip UmatEnergy' )
  ELSE
    CALL ListAddString( SolverParams,&
        NextFreeKeyword('Exported Variable ',SolverParams), &
        '-nooutput -dofs 3 -ip UmatEnergy' )
  END IF

  Nstate = 0
  DO i=1,Model % NumberOfMaterials
    Material => Model % Materials(i) % Values
    IF( ListCheckPresent( Material,'UMAT Subroutine') ) THEN
      Nstate = MAX(Nstate, GetInteger( Material, 'Number of State Variables', Found))
      IF (.NOT. Found) CALL Fatal(Caller, &
          'Number of Material Constants for UMAT must be specified')
    END IF
  END DO
  
  CALL Info(Caller,'Maximum number of state variables in UMAT: '//I2S(Nstate),Level=7)
  
  ! Create variables for some state variables of a user-defined material model (UMAT):
  ! Note that Elmer does not like length of zero for the variables.
  IF( NState > 0 ) THEN
    IF (OutputStateVars) THEN
      str = '-dofs '//I2S(NState)//' -ip UmatState'
    ELSE
      str = '-nooutput -dofs '//I2S(NState)//' -ip UmatState'
    END IF
    CALL ListAddString(SolverParams, NextFreeKeyword('Exported Variable ', SolverParams), str )
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
  USE Constitutive
  USE ModelLumping
  USE MainUtils, ONLY : SetGlobalBubblesFlag
  
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
  TYPE(Variable_t), POINTER :: StressSol, TempSol, FlowSol, Var
  TYPE(ValueList_t), POINTER :: SolverParams, Material, PrevMaterial, BC, Equation, BodyForce
  TYPE(Nodes_t) :: ElementNodes, ParentNodes, FlowNodes
  TYPE(Element_t), POINTER :: CurrentElement, ParentElement, FlowElement
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  

  LOGICAL :: GotForceBC, GotFSIBC, GotSpring, GotIt, NewtonLinearization = .FALSE., &
      Isotropic = .TRUE., RotateModuli, LinearModel = .FALSE., MeshDisplacementActive, &
      NeoHookeanMaterial = .FALSE., AxialSymmetry
  ! InputTensor reports whether the heat expansion coefficient was given as a
  ! single scalar. Nothing here needs to know -- the coefficient is expanded onto
  ! the diagonal either way -- but the argument is not optional.
  LOGICAL :: IsotropicHeatExpansion
  LOGICAL :: UseUMAT, InitializeStateVars, HenckyStrain
  LOGICAL :: LargeDeflection
  LOGICAL :: MixedFormulation
  ! "Incompressible": the LINEAR mixed formulation, a pressure unknown with the
  ! div u = 0 constraint. Read per call and not SAVEd, like the rest of the
  ! per-call state -- this solver can be entered while it is already running.
  LOGICAL :: LinearIncompressible
  ! Either mixed formulation: all that the boundary assembly and the stress
  ! postprocessing need to know is that the local stride carries a pressure.
  LOGICAL :: PressureUnknown
  ! Maxwell viscoelasticity: whether any material in the model asks for it, the
  ! integration point lag stress it keeps its history in, and how many independent
  ! components that tensor has in this configuration.
  LOGICAL :: MaxwellMaterial, MaxwellHere
  ! The prestressed eigen analyses: which of the two keywords is active, whether the
  ! term is live on THIS pass (it is not on the static first one), and the sif's own
  ! "Eigen Analysis" setting, which the two passes toggle between.
  LOGICAL :: StabilityAnalysis, GeometricStiffness, GeometricActive, OrigEigenAnalysis
  ! "Quasi Stationary": the transient run without its inertial term.
  LOGICAL :: QuasiStationary
  ! "Gravitational Prestress Advection": whether this element's body force asks for
  ! the term, and the rho*g it is scaled by.
  LOGICAL :: GotGPA
  ! "Stress Pressure", the pressure-like body load.
  LOGICAL :: GotPressureLoad
  ! Handles for material data and loads read AT the integration point, shared by the
  ! assembly and the postprocessing. They are keyword descriptors and nothing else --
  ! the same for every instance of this solver -- so SAVEing them is safe where SAVEing
  ! state would not be.
  TYPE(ValueHandle_t), SAVE :: YoungIP_h, PoissonIP_h, BetaIP_h, LoadIP_h(4)
  LOGICAL, SAVE :: IPHandlesDone = .FALSE.
  ! Material data and loads given AT THE INTEGRATION POINTS rather than at the nodes:
  ! "Youngs Modulus at IP", "Poisson Ratio at IP", "Heat Expansion Coefficient IP" and
  ! "Stress Bodyforce at IP". What they change is only WHERE the value comes from.
  LOGICAL :: EvalYoungIP, EvalPoissonIP, EvalBetaIP, EvalLoadIP
  TYPE(Variable_t), POINTER :: VeStress => NULL()
  INTEGER :: nve
  LOGICAL :: PseudoTraction, GlobalPseudoTraction
  LOGICAL :: PlaneStress, CalculateStrains, CalculateStresses
  LOGICAL :: CalcPrincipalAngle, CalcPrincipal
  LOGICAL :: CalcPrincipalStress, CalcPrincipalStrain
  LOGICAL :: AllocationsDone = .FALSE., HarmonicAnalysis
  ! Whether the SAVEd per-element storage has to be (re)allocated for THIS call.
  LOGICAL :: Realloc
  LOGICAL :: ConstantBulkMatrix, ConstantBulkSystem, ConstantSystem
  LOGICAL :: ConstantBulkMatrixInUse
  LOGICAL :: CompressibilityDefined = .FALSE.
  LOGICAL :: NormalSpring, NormalTangential
  LOGICAL :: Converged, NoExternalLoads
  LOGICAL :: Parallel, Scanning

  INTEGER, POINTER :: TempPerm(:),StressPerm(:),PressPerm(:),NodeIndexes(:), &
          Indices(:), FlowPerm(:), AdjacentNodes(:)

  INTEGER :: dim,i,j,k,l,m,n,nd,nb,ntot,t,iter,NDeg,STDOFs,LocalNodes,istat
  INTEGER :: NonlinearIter, MinNonlinearIter, FlowNOFNodes, previ
  INTEGER :: EigenModes, Passes
  INTEGER :: RelIntegOrder
  INTEGER :: CoordinateSystem
  INTEGER :: NPROPS, NSTATEV, MAXSTATEV

  REAL(KIND=dp), POINTER :: Temperature(:),Pressure(:),Displacement(:), UWrk(:,:), &
       Work(:,:,:) => NULL(), ForceVector(:), Velocity(:,:), FlowSolution(:), SaveValues(:), &
       NodalStrain(:), NodalStress(:), VonMises(:), &
       PrincipalStress(:), PrincipalStrain(:), Tresca(:), PrincipalAngle(:)
  REAL(KIND=dp), POINTER :: MaterialConstants(:,:)
  REAL(KIND=dp), POINTER :: TotalSol(:) => NULL()
  REAL(KIND=dp), POINTER CONTIG :: ValuesSaved(:) => NULL()
  REAL(KIND=dp), POINTER :: Pb(:)

  REAL(KIND=dp), ALLOCATABLE :: LocalMassMatrix(:,:),LocalStiffMatrix(:,:),&
       LocalDampMatrix(:,:),LoadVector(:,:),InertialLoad(:,:), Viscosity(:), LocalForce(:), &
       MaxwellViscosity(:), NodalGPA(:), NodalPressureLoad(:), &
       NodalStressLoad(:,:), NodalStrainLoad(:,:), &
       LocalTemperature(:),ElasticModulus(:,:,:),PoissonRatio(:), Density(:), &
       Damping(:), HeatExpansionCoeff(:,:,:),Alpha(:,:),Beta(:), &
       ReferenceTemperature(:),BoundaryDispl(:),LocalDisplacement(:,:), PrevSOL(:), &
       PrevLocalDisplacement(:,:), SpringCoeff(:,:,:), LocalExternalForce(:), &
       DisplacementRot(:), LocalForceSaved(:)
         
  REAL(KIND=dp) :: UNorm, TransformMatrix(3,3), Tdiff, Normal(3), s, UnitNorm, DragCoeff
  REAL(KIND=dp) :: Norm, NonlinTol, NonlinRes0, NonlinRes, time 
  REAL(KIND=dp) :: at,at0

  CHARACTER(LEN=MAX_NAME_LEN) :: str, CompressibilityFlag
  CHARACTER(LEN=MAX_NAME_LEN) :: UMATName 
  CHARACTER(LEN=80) :: UmatModel
  TYPE(C_FUNPTR) :: UMATSubrtn
  
  TYPE(Variable_t), POINTER :: UmatEnergyVar, UmatStressVar, UmatStateVar
  REAL(KIND=dp), POINTER :: UmatEnergy(:), UmatStress(:), UmatState(:)
  REAL(KIND=dp), POINTER :: UmatEnergy0(:),UmatStress0(:), UmatState0(:)
  LOGICAL, ALLOCATABLE :: UmatInitDone(:)
  LOGICAL :: AnyDamping, GotDamping, GotRayleighAlpha, GotRayleighBeta, NeedMass
  REAL(KIND=dp) :: RayleighAlpha, RayleighBeta

  ! "Time derivative order": whether the inertial term is integrated at all. Read
  ! once per call rather than per element as StressSolve does, and deliberately not
  ! SAVEd -- this solver can be entered while it is already running.
  LOGICAL :: SecondOrderTime
  
  ! Model lumping: six load cases whose reactions become one 6x6 spring matrix for
  ! the boundary. State of the run, deliberately NOT in any SAVE list -- it has to
  ! live across the six cases of THIS call and no longer, which is exactly the
  ! mistake PrevSOL used to make.
  LOGICAL :: UseModelLumping
  TYPE(ModelLumping_t) :: Lump

  ! Gates for the refusal of StressSolve-only keywords: true when one is set
  ! somewhere in the model, so that the exact per-element test is paid for only
  ! then. See where they are assigned.
  LOGICAL :: StressOnlyKeywords, StressLoadInBC, ImagLoadInBC

  ! Staged construction, "Update Reference Displacement": which bodies ask for it,
  ! and which of them are to have the solution copied into the reference at the end.
  ! Allocated per call rather than SAVEd, so a nested instance has its own.
  LOGICAL :: UpdateReference
  LOGICAL, ALLOCATABLE :: UpdatePresent(:), UpdateActive(:)
  TYPE(Variable_t), POINTER :: ReferenceSol

  CHARACTER(*), PARAMETER :: Caller = 'ElasticSolver'

  
!------------------------------------------------------------------------------
  SAVE LocalMassMatrix,LocalStiffMatrix,LocalDampMatrix,LoadVector,InertialLoad, Viscosity, &
       MaxwellViscosity, NodalGPA, NodalPressureLoad, &
       NodalStressLoad, NodalStrainLoad, Work, &
       LocalForce,ElementNodes,ParentNodes,FlowNodes,Alpha,Beta, &
       LocalTemperature,AllocationsDone,ReferenceTemperature,BoundaryDispl, &
       ElasticModulus, PoissonRatio,Density,Damping,HeatExpansionCoeff, &
       LocalDisplacement, Velocity, Pressure, CalculateStrains, CalculateStresses, &
       NodalStrain, NodalStress, VonMises, PrincipalStress, PrincipalStrain, &
       Tresca, PrincipalAngle, CalcPrincipalAngle, CalcPrincipal, &
       PrevLocalDisplacement, SpringCoeff, Indices
  SAVE MAXSTATEV, InitializeStateVars, TotalSol, LocalExternalForce
  SAVE UmatEnergyVar, UmatStressVar, UmatStateVar, UmatEnergy, UmatStress, UmatState, &
      UmatEnergy0, UmatStress0, UmatState0, UmatInitDone
!-----------------------------------------------------------------------------------------------------
  INTERFACE
    SUBROUTINE ElasticSolver_Boundary_Residual( Model,Edge,Mesh,Quant,Perm, Gnorm,Indicator)
      USE Types
      TYPE(Element_t) :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t) :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END SUBROUTINE ElasticSolver_Boundary_Residual

    SUBROUTINE ElasticSolver_Edge_Residual( Model,Edge,Mesh,Quant,Perm,Indicator)
      USE Types
      TYPE(Element_t) :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t) :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END SUBROUTINE ElasticSolver_Edge_Residual

    SUBROUTINE ElasticSolver_Inside_Residual( Model,Element,Mesh,Quant,Perm, Fnorm,Indicator)
      USE Types
      TYPE(Element_t) :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t) :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END SUBROUTINE ElasticSolver_Inside_Residual
  END INTERFACE


  !------------------------------------------------------------------------------
  !    Get variables needed for solution
  !------------------------------------------------------------------------------
  CALL Info( Caller, '----------------------------------',Level=5)
  CALL Info( Caller, 'Starting Elasticity Solver', Level=5 )
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  
  SolverParams => GetSolverParams()
  Mesh => GetMesh()
  dim = CoordinateSystemDimension()
  CoordinateSystem = CurrentCoordinateSystem()
  AxialSymmetry = CoordinateSystem == AxisSymmetric .OR. &
      CoordinateSystem == CylindricSymmetric
  
  IF ( .NOT. ( CoordinateSystem == Cartesian .OR. AxialSymmetry) ) THEN
    CALL Fatal(Caller, 'Unsupported coordinate system')
  END IF

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
    IF( .NOT. ListGetLogicalAnyBC( Model,'FSI BC' ) ) THEN
      CALL Warn(Caller,'Note that "FSI BC" is not activated automatically any more!')
    END IF
  ELSE
    IF( ListGetLogicalAnyBC( Model,'FSI BC' ) ) THEN
      CALL Warn(Caller,'FSI BC requires flow field that is not available')
    END IF
  END IF

  MeshDisplacementActive = ListGetLogical( SolverParams, &
       'Displace Mesh', GotIt )

  ! An eigen or harmonic analysis has no displacement field to speak of -- the
  ! solution is a set of modes, and what the variable happens to hold afterwards
  ! is not a deformation of the mesh. So do not displace by it unless asked to,
  ! which is what StressSolve has always done.
  IF ( .NOT. GotIt ) MeshDisplacementActive = .NOT. EigenOrHarmonicAnalysis()

  HarmonicAnalysis = getLogical( SolverParams, 'Harmonic Analysis', GotIt ) .OR. &
      getLogical( SolverParams,'Harmonic Mode',GotIt )

  ! Sometimes we might want to use this solver to provide also eigenmode or harmonic analysis.
  ! Then we need to add also the mass even though the system is not transient.
  IF( TransientSimulation ) THEN
    NeedMass = .FALSE.
  ELSE
    NeedMass = EigenOrHarmonicAnalysis() .OR. HarmonicAnalysis
  END IF

  ! Anything other than two means the first order transient: the mass matrix is
  ! dropped and the damping matrix carries the time derivative, as in StressSolve.
  ! Absent counts as two, so that a sif that never mentions the keyword keeps the
  ! behaviour this solver has always had even if _Init did not run.
  SecondOrderTime = .TRUE.
  IF( TransientSimulation ) THEN
    i = ListGetInteger( SolverParams,'Time derivative order', GotIt )
    IF( GotIt ) SecondOrderTime = ( i == 2 )
  END IF
    
  
  IF ( AllocationsDone .AND. MeshDisplacementActive ) THEN
     CALL DisplaceMesh( Mesh, Displacement, -1, StressPerm, STDOFs, UpdateDirs=dim )
  END IF

  !-------------------------------------------------------------------------
  !    Check how material behaviour is defined: 
  !-------------------------------------------------------------------------
  !
  ! The only way to make the umat 
  ! version active is to have "UMAT Subroutine" as specified.
  !
  UseUMAT = ListCheckPresentAnyMaterial(Model, 'UMAT Subroutine')
  IF (UseUMAT .AND. TransientSimulation) THEN
    CALL Fatal(Caller, 'UMAT version does not yet support transient simulation')
  END IF

  PrevMaterial => NULL()   
  NeoHookeanMaterial = ListGetLogical( SolverParams, 'Neo-Hookean Material', GotIt )
  IF (NeoHookeanMaterial) Isotropic = .TRUE.
  MixedFormulation = NeoHookeanMaterial .AND. &
      ListGetLogical( SolverParams, 'Mixed Formulation', GotIt )
  IF (MixedFormulation .AND. (STDOFs /= (dim + 1))) CALL Fatal(Caller, &
      'With mixed formulation variable DOFs should equal to space dimensions + 1')

  !----------------------------------------------------------------------------
  ! "Incompressible" -- StressSolve's mixed formulation for a LINEAR material
  ! that is truly incompressible. It is the small strain specialisation of what
  ! "Mixed Formulation" does for the neo-Hookean law: the same pressure unknown
  ! and the same lowest-order pressure basis, but the constraint is div u = 0
  ! rather than a linearised det F = 1, so it belongs in LocalMatrix and not in
  ! NeoHookeanLocalMatrix. The two keywords are therefore exclusive.
  !
  ! Note that this ignores "Poisson Ratio" and takes mu = E/3, which is what
  ! StressSolve does: the ratio is not read as a value near 1/2, it is 1/2 by
  ! construction and the volumetric part is carried by the pressure instead.
  !
  ! Only the tests whose data is to hand are made here. "Large Deflection" is read
  ! further down, and plane stress and isotropy are per element, so they are
  ! refused in the assembly loop where THIS element's equation and material are
  ! known -- a migration sif with bodies through both solvers must not be refused
  ! for what the other solver's bodies ask.
  !----------------------------------------------------------------------------
  LinearIncompressible = ListGetLogical( SolverParams, 'Incompressible', GotIt )
  IF ( LinearIncompressible ) THEN
    IF ( MixedFormulation .OR. NeoHookeanMaterial ) CALL Fatal(Caller, &
        '"Incompressible" is the linear mixed formulation and cannot be combined '// &
        'with the neo-Hookean one')
    IF ( UseUMAT ) CALL Fatal(Caller, &
        '"Incompressible" is implemented by LocalMatrix and not by the UMAT assembly')
    IF ( STDOFs /= dim + 1 ) CALL Fatal(Caller, &
        '"Incompressible" needs one more variable DOF than there are dimensions, '// &
        'for the pressure: '//I2S(STDOFs)//' given for dim '//I2S(dim))
    IF ( AxialSymmetry ) CALL Fatal(Caller, &
        '"Incompressible" is not implemented under axial symmetry: the hoop term '// &
        'reaches the constraint as well as the divergence, and StressSolve refuses '// &
        'the same combination')
    ! Accepted, not refused -- see the assembly loop for why -- but not in silence.
    ! Once, rather than once per timestep: the Maxwell sifs this is here for run
    ! twenty steps, and a warning repeated twenty times reads as noise.
    IF ( Solver % DoneTime <= 1 .AND. &
        ListGetLogicalAnyEquation( Model, 'Plane Stress' ) ) CALL Warn(Caller, &
        '"Incompressible" overrides "Plane Stress": div u = 0 in the plane is the '// &
        'plane strain constraint. StressSolve ignores the keyword here as well, '// &
        'without saying so.')
  END IF
  PressureUnknown = MixedFormulation .OR. LinearIncompressible

  !----------------------------------------------------------------------------
  ! MAXWELL VISCOELASTICITY. The law is the linear one with its stiffness scaled
  ! by xPhi = 1/(1 + mu/eta dt) and an additive stress offset carrying the lag
  ! stress from the previous step -- so it is a RHEOLOGY around the existing
  ! constitutive model rather than a new model, which is what §3's scouting said
  ! and what makes it cheap. What it does need is per-integration-point history,
  ! which is the "ve_stress" -ip variable that _Init exports, read as v0 and
  ! written as v through MaterialState_t.
  !
  ! Read per material and not per solver, because that is where StressSolve reads
  ! it and a body may be viscoelastic while its neighbour is not -- the model-wide
  ! test here only decides whether any of the machinery is needed at all.
  !----------------------------------------------------------------------------
  MaxwellMaterial = ListGetLogicalAnyMaterial( Model, 'Maxwell material' )
  IF ( MaxwellMaterial ) THEN
    IF ( .NOT. TransientSimulation ) CALL Fatal( Caller, &
        '"Maxwell material" is a rate law and needs a transient simulation' )
    IF ( UseUMAT .OR. NeoHookeanMaterial ) CALL Fatal( Caller, &
        '"Maxwell material" is implemented by LocalMatrix, so not with UMAT or the '// &
        'neo-Hookean model -- a user routine or a finite strain law owns its own '// &
        'rate behaviour' )

    ! The lag stress at the integration points. Its component count is the number
    ! of independent components of a symmetric tensor in this configuration, and
    ! _Init sizes the variable by the same rule; if the two ever disagree the
    ! indexing below would run into the neighbouring point, so say so instead.
    VeStress => VariableGet( Mesh % Variables, 've_stress' )
    IF ( .NOT. ASSOCIATED( VeStress ) ) CALL Fatal( Caller, &
        '"Maxwell material" is set but the "ve_stress" variable is missing' )

    nve = SymTensorComponents( dim, AxialSymmetry )
    IF ( VeStress % DOFs /= nve ) CALL Fatal( Caller, &
        'Variable "ve_stress" has '//I2S(VeStress % DOFs)//' components per point, '// &
        'expected '//I2S(nve) )

    IF ( .NOT. ASSOCIATED( VeStress % PrevValues ) ) THEN
      ALLOCATE( VeStress % PrevValues( SIZE(VeStress % Values), 1 ) )
      VeStress % PrevValues = 0.0_dp
    END IF
  END IF

  !----------------------------------------------------------------------------
  ! The affine offset channel -- "Stress Load" and "Strain Load" -- is
  ! implemented by LocalMatrix and by neither of the other two assembly
  ! routines. Refused here rather than read and dropped: a keyword that is
  ! parsed, interpolated and then silently ignored is the one failure mode this
  ! solver has produced repeatedly, and it looks exactly like a converged answer.
  !----------------------------------------------------------------------------
  IF ( UseUMAT .OR. NeoHookeanMaterial ) THEN
    IF ( ListCheckPrefixAnyBodyForce( Model, 'Stress Load' ) .OR. &
         ListCheckPrefixAnyBodyForce( Model, 'Strain Load' ) ) THEN
      CALL Fatal( Caller, '"Stress Load" and "Strain Load" are available for the '// &
          'linear elastic material models only, not with UMAT or the neo-Hookean model' )
    END IF
  END IF

  ! Thermal expansion travels the same channel, so the same applies -- but to the
  ! neo-Hookean model only. A UMAT owns its constitutive law entirely and is handed
  ! the temperature to do with as it likes, so the keyword is its business there
  ! rather than something this solver has dropped.
  IF ( NeoHookeanMaterial ) THEN
    IF ( ListCheckPresentAnyMaterial( Model, 'Heat Expansion Coefficient' ) ) THEN
      CALL Fatal( Caller, '"Heat Expansion Coefficient" is available for the linear '// &
          'elastic material models only, not with the neo-Hookean model' )
    END IF
  END IF

  AnyDamping = ListCheckPresentAnyMaterial( Model,"Damping" ) .OR. &
      ListCheckPrefixAnyMaterial( Model,"Rayleigh" )
  GotDamping = .FALSE.
  GotRayleighAlpha = .FALSE.
  GotRayleighBeta = .FALSE.
  
  !------------------------------------------------------------------------------
  !     Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  N = Mesh % MaxElementDOFs

  !----------------------------------------------------------------------------
  ! The storage below is SAVEd, so it is shared by every ElasticSolve instance in
  ! the sif and was sized by whichever of them ran first. A second instance with
  ! more degrees of freedom per node -- one under a mixed formulation beside a
  ! plain one -- would index arrays sized for the other. So the guard asks whether
  ! what is there is big enough for THIS instance, and not merely whether somebody
  ! allocated something once.
  !
  ! Grown and never shrunk: an instance that finds the arrays larger than it asked
  ! for leaves them alone, which keeps the one that allocated them safe. LocalForce
  ! carries both factors, STDOFs and the element's degrees of freedom, so it is the
  ! one to measure.
  !
  ! Two statements and not one expression, because Fortran does not promise to stop
  ! evaluating an .OR. once it is decided, and SIZE of an unallocated array is not
  ! a question with an answer.
  !----------------------------------------------------------------------------
  Realloc = .NOT. AllocationsDone .OR. Solver % MeshChanged
  IF ( .NOT. Realloc ) Realloc = SIZE( LocalForce ) < STDOFs * N

  IF ( Realloc ) THEN
     IF ( AllocationsDone ) THEN
        DEALLOCATE( &
             BoundaryDispl, &
             ReferenceTemperature, &
             HeatExpansionCoeff, &
             LocalTemperature, &
             Pressure, Velocity, &
             ElasticModulus, PoissonRatio, &
             Density, Damping, &
             LocalForce, LocalExternalForce, Viscosity, MaxwellViscosity, NodalGPA, &
             NodalPressureLoad, &
             LocalMassMatrix,  &
             LocalStiffMatrix,  &
             LocalDampMatrix,  &
             LoadVector, InertialLoad, Alpha, Beta, &
             NodalStressLoad, NodalStrainLoad, &
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
          MaxwellViscosity( N ), NodalGPA( N ), NodalPressureLoad( N ), &
          LocalMassMatrix(  STDOFs*N,STDOFs*N ),  &
          LocalStiffMatrix( STDOFs*N,STDOFs*N ),  &
          LocalDampMatrix( STDOFs*N,STDOFs*N ),  &
          LoadVector( 4,N ), InertialLoad(3,N), Alpha( 3,N ), Beta( N ), &
          NodalStressLoad( 6,N ), NodalStrainLoad( 6,N ), &
          LocalDisplacement( 4,N ), &
          PrevLocalDisplacement( 4,N ), &
          SpringCoeff( N,3,3 ), &
          Indices(N), &
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( Caller,  'Memory allocation error.' )
     END IF

     IF (UseUMAT .AND. (.NOT. AllocationsDone) ) THEN
       ! ---------------------------------------------------------------------
       ! Stress and energy variables are always created when a UMAT subroutine
       ! is used. Get pointers to these variables:
       ! ---------------------------------------------------------------------
              
       UmatEnergyVar => VariableGet( Mesh % Variables, 'UmatEnergy')
       IF (.NOT. ASSOCIATED( UmatEnergyVar ) ) THEN
         CALL Fatal(Caller,'Could not find variable "UmatEnergy"')
       END IF

       UmatStressVar => VariableGet( Mesh % Variables, 'UmatStress')
       IF (.NOT. ASSOCIATED( UmatStressVar ) ) THEN
         CALL Fatal(Caller,'Could not find variable "UmatStress"')
       END IF

       UmatEnergy => UmatEnergyVar % Values
       UmatStress => UmatStressVar % Values                

       UmatEnergy = 0.0_dp
       UmatStress = 0.0_dp

       ! ----------------------------------------------------------------------
       ! We also create similar variables with suffix "0" to keep the variable
       ! values corresponding to the converged solution at the previous time 
       ! level m. In addition to the stress and energy variables, we need to 
       ! save the state variables as they evolve during the nonlinear iteration 
       ! to obtain the solution at the new time level m+1. The right values
       ! of the state variables corresponding to the initial state can be found
       ! by making an extra UMAT call. Check whether this call is needed.
       ! ----------------------------------------------------------------------

       ALLOCATE( UmatEnergy0( SIZE( UmatEnergy ) ) ) 
       ALLOCATE( UmatStress0( SIZE( UmatStress ) ) ) 
              
       UmatEnergy0 = 0.0_dp
       UmatStress0 = 0.0_dp
       
       UmatStateVar => VariableGet( Mesh % Variables, 'UmatState')
       IF( ASSOCIATED( UmatStateVar ) ) THEN
         MaxStateV = UmatStateVar % Dofs
         CALL Info(Caller,'Maximum number of state variables in UMAT: '&
             //I2S(MaxStateV),Level=7)
         UmatState => UmatStateVar % Values         
         ALLOCATE( UmatState0( SIZE( UmatState ) ) )          
         UmatState = 0.0_dp         
         UmatState0 = 0.0_dp
       ELSE
         CALL Info(Caller,'Could not find variable "UmatState", assuming no state variable!')
         MaxStateV = 0
         UmatState => NULL()
         UmatState0 => NULL()
       END IF
      
       ALLOCATE( UmatInitDone( SIZE( UmatEnergy ) / 3 ) )
       UmatInitDone = .FALSE.
       
       InitializeStateVars = GetLogical(SolverParams, 'Initialize State Variables',GotIt)
     END IF

     AllocationsDone = .TRUE.
  END IF

  !----------------------------------------------------------------------------
  ! The damping, zeroed on every call as StressSolve zeroes it before its own
  ! assembly. The array is only written where a material names a damping, so with
  ! none named it would be read uninitialised -- and it IS read even then: the
  ! inertial term is dropped under "Quasi Stationary" only where the element
  ! carries no damping, and garbage there puts the mass silently back. It is SAVEd
  ! and shared between instances besides, so a second solver with no damping would
  ! otherwise inherit the first one's values.
  !----------------------------------------------------------------------------
  Damping = 0.0_dp

  !----------------------------------------------------------------------------
  ! What this solver is asked to postprocess. Read on EVERY call, and outside the
  ! allocation guard above, which is where these four used to sit -- so a second
  ! ElasticSolve instance in the same sif never read its own keywords at all and
  ! silently inherited the first one's, these variables being SAVEd as well.
  !
  ! The symptom was the familiar one: a first solver that does not ask for stresses
  ! leaves CalculateStresses false, the second one asks and is not heard, and its
  ! projected stress comes back a converged ZERO. Note that the size test on the
  ! guard does not reach this -- two instances with the same STDOFs need no
  ! reallocation, so the guard stays shut and the keywords would still go unread.
  ! They are two defects in one place, and this is the half that was measured.
  !----------------------------------------------------------------------------
  CalculateStrains = GetLogical(SolverParams, 'Calculate Strains', GotIt )
  CalculateStresses = GetLogical(SolverParams, 'Calculate Stresses', GotIt )

  IF (UseUMAT) THEN
     ! Principal tensors are not yet available:
     CalcPrincipal = .FALSE.
     CalcPrincipalAngle = .FALSE.
  ELSE
     CalcPrincipal = GetLogical(SolverParams, 'Calculate Principal', GotIt )
     CalcPrincipalAngle = GetLogical(SolverParams, 'Calculate PAngle', GotIt )
     ! Principal angle computation enforces component calculation:
     IF (CalcPrincipalAngle) CalcPrincipal = .TRUE.
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

  ! Allocated on every entry and deallocated before every return, and therefore
  ! deliberately NOT in the SAVE list above -- which it used to be, and that was a
  ! crash rather than an inefficiency.
  !
  ! THIS SOLVER CAN BE ENTERED WHILE IT IS ALREADY RUNNING. A block system names
  ! its assembly slaves with "Pre Solvers", DefaultStart activates them, and if a
  ! slave is another elasticity solver then the second instance reaches this line
  ! with the first instance's array still allocated: "Fortran runtime error:
  ! Attempting to allocate already allocated variable 'prevsol'". Two solid bodies
  ! coupled through a block system is exactly that arrangement, and it is what
  ! fem/tests/ElasticBeamSolidCoupling now covers.
  !
  ! The neighbours here were already correct: DisplacementRot and LocalForceSaved
  ! are in no SAVE list, so each invocation gets its own.
  ALLOCATE( PrevSOL(SIZE(Displacement)) )
  IF (UseUMAT) THEN
    ALLOCATE(DisplacementRot(SIZE(Displacement)))
    ALLOCATE(LocalForceSaved(SIZE(LocalForce)))
  END IF

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
  IF( NonlinearIter > 1 ) THEN
    NonlinTol = GetConstReal(SolverParams, 'Nonlinear System Convergence Tolerance')
  END IF
  MinNonlinearIter = ListGetInteger( SolverParams, &
       'Nonlinear System Min Iterations', GotIt )

  LinearModel = ListGetLogical( SolverParams, &
       'Elasticity Solver Linear', GotIt )
  LargeDeflection = ListGetLogical(SolverParams, 'Large Deflection', GotIt)
  IF (.NOT. GotIt) LargeDeflection = .TRUE.

  IF (.NOT. LargeDeflection) HenckyStrain = .FALSE.

  ! The linear mixed formulation is a small strain one: its constraint is div u = 0
  ! and the residual it rides assumes that one Newton step is exact, which holds
  ! only while the kinematics are linear.
  IF ( LinearIncompressible .AND. LargeDeflection ) CALL Fatal(Caller, &
      '"Incompressible" is a small strain formulation: set "Large Deflection = False", '// &
      'or use the neo-Hookean mixed formulation for a finite strain one')

  ! Maxwell likewise: StressSolve has no finite strain kinematics at all, so there is
  ! nothing to reproduce and the lag stress would be pushed forward by a deformation
  ! gradient the law knows nothing about.
  IF ( MaxwellMaterial .AND. LargeDeflection ) CALL Fatal(Caller, &
      '"Maxwell material" is implemented for small strain only, as in StressSolve: '// &
      'set "Large Deflection = False"')

  !-----------------------------------------------------------------------------
  ! MODEL LUMPING: reduce the body to a 6x6 spring matrix for one boundary, by
  ! solving six load cases -- three translations and three rotations, imposed
  ! either as displacements or as pure forces and moments -- and reading the
  ! reactions off the lumping boundary.
  !
  ! The whole of it is in ModelLumping.F90, which was extracted from StressSolve
  ! as shared code precisely so that this solver could call it, and which this
  ! solver then did not call at all: "Model Lumping" was silently inert here.
  ! Four calls do it, and the module needs nothing else from the driver.
  !
  ! THE SIX CASES ARE LOAD CASES, NOT NEWTON STEPS. They ride the nonlinear
  ! iteration counter, so both bounds are forced to six and the convergence tests
  ! never get a chance to stop early. That is only sound for a law whose residual
  ! is independent of the iterate -- which is what "Constant Bulk System" below
  ! already demands, and its Fatal covers large deflection, UMAT and neo-Hookean.
  !
  ! "Constant Bulk System" is ADDED to the list rather than merely assumed,
  ! because the same keyword does two things and lumping needs both: it makes the
  ! six cases share one assembled stiffness, and it is what makes
  ! DefaultFinishBulkAssembly save BulkValues -- which ModelLumpingSprings
  ! multiplies the solution by to recover the reactions. StressSolve sets its own
  ! internal flag here instead, so a lumping sif that omits the keyword gets the
  ! reuse without the saved bulk matrix; adding it to the list keeps the two in
  ! agreement by construction.
  !-----------------------------------------------------------------------------
  UseModelLumping = ListGetLogical( SolverParams, 'Model Lumping', GotIt )
  IF ( UseModelLumping ) THEN
    IF ( dim /= 3 ) CALL Fatal( Caller, 'Model lumping is implemented in 3D only' )
    NonlinearIter = 6
    MinNonlinearIter = 6
    CALL ListAddLogical( SolverParams, 'Constant Bulk System', .TRUE. )
  END IF

  !-----------------------------------------------------------------------------
  ! GEOMETRIC STIFFNESS, and STABILITY ANALYSIS, which is the same term put in a
  ! different matrix. Both are prestressed EIGENVALUE analyses and neither is a
  ! constitutive matter: the term is
  !
  !   INT sigma_ij d(phi_p)/dx_i d(phi_q)/dx_j , diagonal in the components,
  !
  ! with sigma the stress of the state the structure is already in. So it takes
  ! TWO passes and that is the whole of the orchestration: the first solves the
  ! static problem with the eigen analysis switched OFF, and the second builds the
  ! term from what the first found and solves the eigenproblem with it. Exactly as
  ! StressSolve does it, including forcing both iteration bounds to two so the
  ! convergence tests cannot stop after the static pass.
  !
  ! Where they differ is which matrix the term lands in, and that is what makes one
  ! a vibration problem and the other a buckling problem:
  !
  ! - "Geometric Stiffness" adds it to the STIFFNESS, so the eigenvalues are
  !   frequencies of a prestressed structure against the ordinary mass matrix;
  ! - "Stability Analysis" puts MINUS it in the MASS slot INSTEAD of the density
  !   mass, so the eigenvalues are the load multipliers at which the stiffness
  !   loses definiteness -- the critical loads. This is why the density mass matrix
  !   must not be assembled at all in that case.
  !-----------------------------------------------------------------------------
  StabilityAnalysis = ListGetLogical( SolverParams, 'Stability Analysis', GotIt )
  GeometricStiffness = ListGetLogical( SolverParams, 'Geometric Stiffness', GotIt )

  IF ( StabilityAnalysis .OR. GeometricStiffness ) THEN
    IF ( StabilityAnalysis .AND. GeometricStiffness ) CALL Fatal( Caller, &
        '"Stability Analysis" and "Geometric Stiffness" are the same term in '// &
        'different matrices and cannot both be active' )
    IF ( CoordinateSystem /= Cartesian ) CALL Fatal( Caller, &
        'The geometric stiffness is implemented in cartesian coordinates only, as '// &
        'in StressSolve' )
    IF ( LargeDeflection ) CALL Fatal( Caller, &
        'The geometric stiffness is already part of the tangent under '// &
        '"Large Deflection", where it comes from the current iterate rather than '// &
        'from a prestressed state. Set "Large Deflection = False" to ask for the '// &
        'prestressed eigen analysis these keywords mean' )
    IF ( UseUMAT .OR. NeoHookeanMaterial ) CALL Fatal( Caller, &
        'The geometric stiffness is implemented by LocalMatrix, so not with UMAT '// &
        'or the neo-Hookean model' )

    NonlinearIter = 2
    MinNonlinearIter = 2
    OrigEigenAnalysis = ListGetLogical( SolverParams, 'Eigen Analysis', GotIt )
  END IF

  !-----------------------------------------------------------------------------
  ! STRESSSOLVE KEYWORDS THIS SOLVER DOES NOT IMPLEMENT, refused rather than
  ! ignored. The two solvers are being merged, so sifs will be pointed here that
  ! were written for there, and a keyword that is read and dropped looks exactly
  ! like a converged answer. This solver has already produced that failure twice
  ! over: thermal expansion was parsed, interpolated, passed into the assembly and
  ! never read, and "Model Lumping" was inert while its module sat there shared.
  ! Both returned a plausible number. Neither said anything.
  !
  ! The list IS the remaining retirement backlog, which is why it is in one place.
  !
  ! Two kinds of check, and the split is about cost rather than taste. The solver's
  ! own section is read once here, so those tests are free and exact. The material
  ! and body force keywords cannot be resolved exactly without knowing which lists
  ! this solver's bodies reach, and testing them per element would cost a string
  ! comparison per element per keyword -- of the order of a tenth of the assembly,
  ! which is the same budget the batched constitutive entry point was written to
  ! recover. So a model-wide test gates them: nothing set anywhere costs one
  ! logical per element, and only a sif that does set one pays for the precise
  ! per-element test in the assembly loop, where it Fatals on the first element.
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! "Quasi Stationary": drop the inertial term from a transient run, keeping the
  ! density for whatever else reads it -- "no mass-matrix, despite finite densities",
  ! as the sif that asks for it puts it.
  !
  ! StressSolve conflates the mass with the damping in this one flag: it sets
  ! NeedMass from the keyword, and then any damping present turns NeedMass back on,
  ! which builds BOTH matrices again. That is reproduced rather than tidied, and the
  ! test is made on the element's own nodal values as it is made there -- see where
  ! the mass is integrated.
  !-----------------------------------------------------------------------------
  QuasiStationary = ListGetLogical( SolverParams, 'Quasi Stationary', GotIt )
  IF ( QuasiStationary .AND. .NOT. TransientSimulation ) CALL Warn( Caller, &
      '"Quasi Stationary" only has meaning in a transient simulation' )

  ! StressSolve lets an Equation section name the temperature field it couples to.
  ! This solver reads the variable called "Temperature" and nothing else, so the
  ! keyword would silently point at a field that is never fetched. Implementing it
  ! means a per-BODY lookup, since the Equation section is per body while this
  ! solver's temperature pointer is taken once per call.
  IF ( ListCheckPresentAnyEquation( Model, 'Temperature Name' ) ) CALL Fatal( Caller, &
      '"Temperature Name" is not implemented here: this solver couples to the '// &
      'variable named "Temperature"' )

  IF ( ListCheckPresent( SolverParams, 'Update Transient System' ) ) CALL Fatal( Caller, &
      '"Update Transient System" is not implemented here: this solver''s reuse of an '// &
      'assembled system is governed by the "Constant Bulk Matrix", "Constant Bulk '// &
      'System" and "Constant System" keywords instead' )

  ! Set anywhere in the model, in any material or body force? Then the assembly
  ! loop will test the lists this element actually uses. See the note above.
  StressOnlyKeywords = &
      ListCheckPrefixAnyMaterial( Model, 'Pre Stress' ) .OR. &
      ListCheckPrefixAnyMaterial( Model, 'Pre Strain' ) .OR. &
      ListCheckPrefixAnyMaterial( Model, 'Mesh Velocity' ) .OR. &
      ListCheckPresentAnyBodyForce( Model, 'Stress Pressure im' ) .OR. &
      ListCheckPresentAnyBodyForce( Model, 'Stress Bodyforce 1 im' ) .OR. &
      ListCheckPresentAnyBodyForce( Model, 'Stress Bodyforce 2 im' ) .OR. &
      ListCheckPresentAnyBodyForce( Model, 'Stress Bodyforce 3 im' )

  ! "Stress Load" is implemented as a body force here and not yet as a boundary
  ! condition, which StressSolve also reads it as. Gated the same way.
  StressLoadInBC = ListCheckPrefixAnyBC( Model, 'Stress Load' )

  ! The imaginary half of a boundary load, gated the same way. Real system only here.
  ImagLoadInBC = ListCheckPresentAnyBC( Model, 'Force 1 im' ) .OR. &
      ListCheckPresentAnyBC( Model, 'Force 2 im' ) .OR. &
      ListCheckPresentAnyBC( Model, 'Force 3 im' ) .OR. &
      ListCheckPresentAnyBC( Model, 'Normal Force im' )

  !-----------------------------------------------------------------------------
  ! STAGED CONSTRUCTION -- "Update Reference Displacement" in a body force. The
  ! keyword is a switch on a sign rather than a value: where its nodal values are
  ! mostly non-negative the body's displacement is to BECOME the new reference at
  ! the end of the step, and where they are mostly negative the body is solved
  ! measured FROM the reference already stored, by adding K u_ref to its load. A
  ! sif drives a construction sequence by flipping that sign over time, which is
  ! what fem/tests/staged_sim does.
  !
  ! The reference itself lives in a variable named "Reference Displacement", which
  ! the sif exports. StressSolve tests for the keyword only where that variable
  ! exists, so a sif that forgets to export it gets the keyword ignored without a
  ! word; here that is the one thing refused outright, since the keyword is then
  ! doing nothing and the answer is the unstaged one.
  !-----------------------------------------------------------------------------
  ReferenceSol => VariableGet( Mesh % Variables, 'Reference Displacement' )
  UpdateReference = .FALSE.

  IF ( ListCheckPresentAnyBodyForce( Model, 'Update Reference Displacement' ) ) THEN
    IF ( .NOT. ASSOCIATED( ReferenceSol ) ) CALL Fatal( Caller, &
        '"Update Reference Displacement" needs a variable named "Reference '// &
        'Displacement" to store the reference in: export one from this solver' )

    IF ( LargeDeflection ) CALL Fatal( Caller, &
        '"Update Reference Displacement" needs "Large Deflection = False": the '// &
        'reference is subtracted through the linear stiffness, which is not the '// &
        'tangent of a geometrically nonlinear step' )

    ALLOCATE( UpdatePresent( Model % NumberOfBodies ), &
              UpdateActive( Model % NumberOfBodies ) )
    UpdatePresent = .FALSE.
    UpdateActive = .FALSE.

    DO i = 1, Model % NumberOfBodies
      j = GetInteger( Model % Bodies(i) % Values, 'Body Force', GotIt )
      IF ( .NOT. GotIt ) CYCLE
      UpdatePresent(i) = ListCheckPresent( Model % BodyForces(j) % Values, &
          'Update Reference Displacement' )
    END DO
    UpdateReference = ANY( UpdatePresent )
  END IF

  !-----------------------------------------------------------------------------
  ! "Local Matrix Storage" lets the assembly build one element's local matrix and
  ! reuse it for every element the core has marked identical to it -- by
  ! "Local Matrix Identical" for the whole set, or "... Identical Bodies" per
  ! body. The bookkeeping is all in DefaultStart and UseLocalMatrixCopy; a solver
  ! opts in with the one test in the assembly loop below.
  !
  ! It is only sound while the local matrix depends on nothing element-specific,
  ! which rules out more here than it does in StressSolve. Only the stiffness and
  ! the force are stored, so mass and damping cannot be carried; and a tangent
  ! that reads the current solution differs element by element however identical
  ! the geometry, so the nonlinear paths are out too. Refuse rather than quietly
  ! assemble the wrong matrix -- being wrong here looks like a converged answer.
  !-----------------------------------------------------------------------------
  IF( ListGetLogical( SolverParams, 'Local Matrix Storage', GotIt ) ) THEN
    IF( NeedMass .OR. TransientSimulation ) CALL Fatal( Caller, &
        '"Local Matrix Storage" is applicable to steady cases only' )
    IF( UseUMAT .OR. NeoHookeanMaterial .OR. LargeDeflection ) CALL Fatal( Caller, &
        '"Local Matrix Storage" needs a linear material: set "Large Deflection = False"'//&
        ' and use neither a UMAT nor a neo-Hookean material' )
  END IF

  !-----------------------------------------------------------------------------
  ! Reusing what a previous pass assembled. The three keywords differ in how much
  ! is taken to be unchanged; see the restore block in the assembly below.
  !
  ! The same soundness conditions as the local matrix cache apply, for the same
  ! reason: a stiffness held over from the previous pass is only the right
  ! stiffness while it does not depend on the solution, and only the bulk matrix
  ! and its load are saved, so mass and damping cannot be carried over. StressSolve
  ! supports the transient case by integrating the restored matrix in time
  ! afterwards (its AddGlobalTime); that is not done here, so refuse it rather
  ! than appear to honour the keyword.
  !-----------------------------------------------------------------------------
  ! Relative order of the integration rule, as StressSolve has always taken it.
  ! Zero when the keyword is absent, which is the rule GaussPoints would have
  ! chosen anyway. Matters most for p-elements, where the default rule is the one
  ! the element declares; fem/tests/ElastPelem2dPmultg* are StressSolve cases that
  ! turn it down to keep a p-refined solve affordable.
  RelIntegOrder = ListGetInteger( SolverParams,'Relative Integration Order', GotIt )

  ConstantSystem     = ListGetLogical( SolverParams, 'Constant System', GotIt )
  ConstantBulkSystem = ListGetLogical( SolverParams, 'Constant Bulk System', GotIt )
  ConstantBulkMatrix = ListGetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )

  IF( ConstantSystem .OR. ConstantBulkSystem .OR. ConstantBulkMatrix ) THEN
    IF( NeedMass .OR. TransientSimulation ) CALL Fatal( Caller, &
        'The "Constant ..." keywords are applicable to steady cases only here' )
    IF( UseUMAT .OR. NeoHookeanMaterial .OR. LargeDeflection ) CALL Fatal( Caller, &
        'The "Constant ..." keywords need a linear material: set'//&
        ' "Large Deflection = False" and use neither a UMAT nor a neo-Hookean material' )
  END IF

  ! The tests below are in order of increasing boldness, so the weaker keyword
  ! wins if both are given. StressSolve orders them the same way.
  IF( ConstantSystem .AND. ( ConstantBulkSystem .OR. ConstantBulkMatrix ) ) THEN
    CALL Warn( Caller,'"Constant System" is superseded by the narrower '//&
        '"Constant Bulk System"/"Constant Bulk Matrix" given beside it' )
  END IF

  GlobalPseudoTraction = GetLogical( SolverParams, 'Pseudo-Traction', GotIt)


  ! If we need the previous timestep for UMAT, what is the step we need?
  previ = 0
  IF (UseUMAT) THEN
    IF( TransientSimulation ) THEN
      previ = 3
    ELSE IF( Scanning ) THEN
      previ = 1
    END IF
  END IF
  ! Maxwell needs it too, for the pressure of the previous step under the mixed
  ! formulation and for the velocity a non-Newtonian viscosity is a function of.
  ! Column one and not three: the second order transient keeps the displacement at
  ! t-dt in the third, but a Maxwell material is integrated first order, where it
  ! is the first -- which is also what GetVectorLocalSolution(tStep=-1) returns,
  ! the call StressSolve makes for the same data.
  IF ( MaxwellMaterial .AND. TransientSimulation ) previ = 1
  IF(previ > 0) THEN
    CALL Info('ElasticSolver','Taking previous displacement from PrevValues(:,'//I2S(previ)//')',Level=30)
  END IF
  
  
  time = GetTime()

  ! The geometry of the lumping boundary -- its area, centre and second moments --
  ! and the rigid body mass. Computed once, before the load cases, since every case
  ! reads it and none of it moves.
  IF ( UseModelLumping ) CALL ModelLumpingInit( Lump, Solver, Model )

  CALL DefaultStart()

  ! The keyword handles for integration point data, bound once for the whole run.
  ! Bound whether or not any sif asks for them, so that no binding happens inside an
  ! element loop.
  IF ( .NOT. IPHandlesDone ) THEN
    CALL ListInitElementKeyword( YoungIP_h, 'Material', 'Youngs Modulus' )
    CALL ListInitElementKeyword( PoissonIP_h, 'Material', 'Poisson Ratio' )
    CALL ListInitElementKeyword( BetaIP_h, 'Material', 'Heat Expansion Coefficient' )
    DO i=1,3
      CALL ListInitElementKeyword( LoadIP_h(i), 'Body Force', 'Stress Bodyforce '//I2S(i) )
    END DO
    CALL ListInitElementKeyword( LoadIP_h(4), 'Body Force', 'Stress Pressure' )
    IPHandlesDone = .TRUE.
  END IF

  DO iter=1,NonlinearIter

     ! The two passes of a prestressed eigen analysis. The first is a static solve
     ! whose only purpose is the stress state, so the eigen analysis is switched off
     ! for it and the geometric term is not yet there to add; the second has both.
     ! The keyword is rewritten rather than a flag kept, because it is the core that
     ! reads it when the system is solved.
     GeometricActive = .FALSE.
     IF ( StabilityAnalysis .OR. GeometricStiffness ) THEN
       GeometricActive = ( iter > 1 )
       CALL ListAddLogical( SolverParams, 'Eigen Analysis', &
           GeometricActive .AND. OrigEigenAnalysis )
     END IF

     at  = CPUTime()
     at0 = RealTime()

     CALL Info( Caller, ' ', Level=7 )
     CALL Info( Caller,'-------------------------------------', Level=5 )
     CALL Info( Caller,'ELASTICITY ITERATION '//I2S(iter),Level=4)
     CALL Info( Caller,'-------------------------------------', Level=5 )
     CALL Info( Caller, ' ', Level=7 )
     CALL Info( Caller, 'Starting assembly...', Level=7 )

     IF (UseUMAT) TotalSol(:) = Displacement(:)

     !------------------------------------------------------------------------------
100  CALL DefaultInitialize()

     !---------------------------------------------------------------------------
     ! Reuse what the previous pass assembled, at one of three depths. The bulk
     ! matrix and right hand side are saved into BulkValues/BulkRHS by
     ! DefaultFinishBulkAssembly below whenever any of these keywords is set, so
     ! there is nothing to restore on the first pass and the tests below simply
     ! fall through to a full assembly.
     !
     !   Constant Bulk Matrix -- the matrix is reused, the load is reassembled.
     !     The element loop still runs, but only the force is glued (see the
     !     ConstantBulkMatrixInUse test at its end).
     !   Constant Bulk System -- matrix and bulk load both reused; the element
     !     loop is skipped entirely and only the boundary conditions run.
     !   Constant System -- the boundary conditions are constant too, so nothing
     !     is assembled at all.
     !
     ! The motivating case is model lumping, whose six load cases share one
     ! stiffness and differ only in what is imposed on the lumping boundary.
     !---------------------------------------------------------------------------
     ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
         ASSOCIATED( Solver % Matrix % BulkValues )

     IF ( ASSOCIATED( Solver % Matrix % BulkValues ) ) THEN
       IF ( ConstantBulkMatrix .OR. ConstantBulkSystem .OR. ConstantSystem ) THEN
         Solver % Matrix % Values = Solver % Matrix % BulkValues
       END IF

       IF ( ConstantBulkSystem .OR. ConstantSystem ) THEN
         Solver % Matrix % RHS = Solver % Matrix % BulkRHS
       ELSE IF ( ConstantBulkMatrix ) THEN
         Solver % Matrix % RHS = 0.0_dp
       END IF

       IF ( ConstantBulkSystem ) GO TO 2000
       IF ( ConstantSystem )     GO TO 3000
     END IF
     !------------------------------------------------------------------------------
     DO t=1,GetNOFActive()
      
        IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                (Solver % NumberOfActiveElements-t) / &
                (1.0*Solver % NumberOfActiveElements)), ' % done'

           CALL Info( Caller, Message, Level=7 )
           at0 = RealTime()
        END IF

        CurrentElement => GetActiveElement(t)
        CALL GetElementNodes(ElementNodes, CurrentElement)
        NodeIndexes => CurrentElement % NodeIndexes

        n = GetElementNOFNodes()
        nd = GetElementDOFs( Indices )
        nb = GetElementNOFBDOFs()
        ntot = nd + nb

        ! The core has marked this element as identical to one already assembled,
        ! so its local matrix is in store and DefaultUpdateEquations will fetch it
        ! rather than read what is passed. Skip building it. See the guard on
        ! "Local Matrix Storage" above for when this is sound.
        IF( .NOT. ConstantBulkMatrixInUse ) THEN
          IF( UseLocalMatrixCopy( Solver, activeind = t ) ) GOTO 200
        END IF

        !-----------------------------------------------------------------------------------
        !        Get the material parameters relating to the constitutive law:
        !------------------------------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()

        ! One logical when nothing anywhere in the model asks for a StressSolve-only
        ! capability, which is the ordinary case; the exact test only when something
        ! does, and then it stops on this element. See where the gate is set.
        IF ( StressOnlyKeywords ) CALL RefuseStressSolveKeywords( Material )

        IF ( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
          IF ( UseUMAT ) THEN
            UMATName = ListGetString(Material, 'UMAT Subroutine', UnfoundFatal=.TRUE.)
            UMATSubrtn = GetProcAddr( UMATName )
            NPROPS = GetInteger(Material,'Number of Material Constants', GotIt)
            NSTATEV = GetInteger(Material,'Number of State Variables', GotIt)
          END IF
          PrevMaterial => Material
        END IF
        
        PlaneStress = GetLogical( Equation, 'Plane Stress', GotIt )
        PoissonRatio = 0.0d0

        ! Data given at the integration points, per material as StressSolve reads it.
        ! An IP-evaluated quantity is one whose keyword is a procedure or a table
        ! wanting the position rather than a nodal value -- the permafrost case passes
        ! functions of ground temperature and water content this way.
        !
        ! Read HERE, before anything reads the material: these decide whether the nodal
        ! reads below may happen at all, and ListGetReal refuses such a keyword outright.
        ! Read late they would carry the previous element's answer, or nothing at all on
        ! the first -- which is exactly how this was found.
        EvalYoungIP   = GetLogical( Material, 'Youngs Modulus at IP', GotIt )
        EvalPoissonIP = GetLogical( Material, 'Poisson Ratio at IP', GotIt )
        EvalBetaIP    = GetLogical( Material, 'Heat Expansion Coefficient IP', GotIt )

        ! The constraint div u = 0 in two dimensions is the plane STRAIN statement, so
        ! it contradicts plane stress -- where it is the out-of-plane strain that
        ! preserves the volume. StressSolve resolves that by ignoring the keyword
        ! outright: its mixed branch never consults it, and a sif setting both gets
        ! the same answer as one setting neither, measured. Every Maxwell sif in the
        ! tree sets both, so the combination has to be accepted rather than refused --
        ! but said out loud, once, where StressSolve says nothing. Warned before the
        ! loop; here the flag is only cleared, so the rest of the assembly and the
        ! postprocessing agree about what is being solved.
        IF ( LinearIncompressible ) PlaneStress = .FALSE.

        IF (UseUMAT) THEN
           CALL GetConstRealArray( Material, MaterialConstants, 'Material Constants', GotIt)
           IF ( SIZE(MaterialConstants,1) < NPROPS) &
                CALL Fatal(Caller,'Check the size of Material Constants array')
           UMATModel = ListGetString(Material, 'Name', GotIt)
        ELSE
           IF (NeoHookeanMaterial) THEN
              ElasticModulus(1,1,1:n) = ListGetReal( Material, &
                   'Youngs Modulus', n, NodeIndexes, GotIt )
            ELSE
              ! Asked for at the integration points instead, so it must not be read
              ! at the nodes at all: ListGetReal refuses a keyword whose value is a
              ! function of quantities that live on the points. Isotropy is assumed
              ! with it, which is the same assumption StressSolve makes -- a matrix
              ! valued modulus at the points is not offered by either.
              IF ( EvalYoungIP ) THEN
                 ElasticModulus = 0.0_dp
                 Isotropic = .TRUE.
              ELSE
                 CALL InputTensor( ElasticModulus, Isotropic, &
                      'Youngs Modulus', Material, n, NodeIndexes )
              END IF

              ! Isotropy is known only now, from this element's own material, which
              ! is where the other half of the "Incompressible" refusal belongs.
              IF ( LinearIncompressible .AND. .NOT. Isotropic ) CALL Fatal(Caller, &
                  '"Incompressible" is implemented for an isotropic material only, as '// &
                  'in StressSolve, which takes mu = E/3 rather than a matrix of moduli')
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
                          CALL Fatal(Caller,'Given > Material Coordinate Unit Vector < too short!')
                       END IF
                       TransformMatrix(i,1:3) = Uwrk(1:3,1) / UnitNorm  
                       RotateModuli = .TRUE.
                    END IF
                    IF( .NOT. RotateModuli  ) CALL Fatal( Caller, &
                         'No unit vectors found but > Rotate Elasticity Tensor < set True?' )
                 END DO
              END IF
           END IF
           ! Not read under "Incompressible": the ratio IS 1/2 there, lambda is
           ! dropped and mu = E/3, so whatever the material says about it is
           ! overridden before it is used. StressSolve guards the read the same
           ! way; without the guard every Maxwell sif -- none of which names a
           ! ratio -- warns once per element per timestep about a keyword whose
           ! value it would then discard.
           IF ( Isotropic .AND. .NOT. EvalPoissonIP .AND. .NOT. LinearIncompressible ) &
               PoissonRatio(1:n) = GetReal( Material, 'Poisson Ratio' )
        END IF
        
        ! Scalar, one value per direction, or a full tensor -- InputTensor decides
        ! from the shape of what the sif gave, and fills the diagonal in the first
        ! two cases. StressSolve reads the keyword through this same routine, so the
        ! spellings a thermal sif may use are the same on both solvers.
        !
        ! Only the DIAGONAL is used, here as in StressSolve: the thermal eigenstrain
        ! is diag(alpha) * dT, so an off-diagonal expansion coefficient is read and
        ! then ignored by both. That is pre-existing and deliberately left alone.
        IF ( EvalBetaIP ) THEN
          HeatExpansionCoeff = 0.0_dp
          IsotropicHeatExpansion = .TRUE.
          GotIt = .TRUE.
        ELSE
          CALL InputTensor( HeatExpansionCoeff, IsotropicHeatExpansion, &
              'Heat Expansion Coefficient', Material, n, NodeIndexes, GotIt )
        END IF
        ReferenceTemperature(1:n) = GetReal( Material, 'Reference Temperature', GotIt )
        
        Density(1:n) = GetReal( Material, 'Density', GotIt )

        ! Viscoelasticity is a property of THIS element's material: a body may be
        ! Maxwell while its neighbour is elastic, which is what the earth tests do.
        ! MaxwellMaterial above only says whether anything in the model asks for it.
        MaxwellHere = .FALSE.
        IF( MaxwellMaterial ) THEN
          MaxwellHere = GetLogical( Material, 'Maxwell material', GotIt )
          IF( MaxwellHere ) THEN
            MaxwellViscosity(1:n) = GetReal( Material, 'Viscosity', GotIt )
            IF( .NOT. GotIt ) CALL Fatal( Caller, &
                '"Maxwell material" needs a "Viscosity": without one the relaxation '// &
                'time is zero and the stiffness vanishes' )
            IF( .NOT. Isotropic ) CALL Fatal( Caller, &
                '"Maxwell material" is implemented for an isotropic material only: '// &
                'the relaxation needs a shear modulus, and an elasticity matrix does '// &
                'not offer one' )
          END IF
        END IF

        IF( AnyDamping ) THEN
          Damping(1:n) = GetReal( Material, 'Damping' ,GotDamping )
          RayleighAlpha = GetCReal( Material, 'Rayleigh Damping alpha',GotRayleighAlpha )
          RayleighBeta = GetCReal( Material, 'Rayleigh Damping beta', GotRayleighBeta )
        END IF
                
        !------------------------------------------------------------------------------
        !        Set body forces
        !------------------------------------------------------------------------------
        BodyForce => GetBodyForce()
        
        LoadVector = 0.0D0
        InertialLoad = 0.0D0
        NodalStressLoad = 0.0D0
        NodalStrainLoad = 0.0D0
        GotGPA = .FALSE.
        NodalGPA = 0.0D0
        GotPressureLoad = .FALSE.
        NodalPressureLoad = 0.0D0
        EvalLoadIP = .FALSE.

        IF ( ASSOCIATED(BodyForce) ) THEN
          ! Read at the points below instead when the keyword says so, and then not
          ! here: the same restriction as on the material data above.
          EvalLoadIP = GetLogical( BodyForce, 'Stress Bodyforce at IP', GotIt )

          IF( ListCheckPrefix(BodyForce,'Stress Bodyforce') .AND. .NOT. EvalLoadIP ) THEN
            LoadVector(1,1:n) = GetReal( BodyForce, 'Stress Bodyforce 1', GotIt )
            LoadVector(2,1:n) = GetReal( BodyForce, 'Stress Bodyforce 2', GotIt )
            IF ( dim > 2 ) THEN
              LoadVector(3,1:n) = GetReal( BodyForce, 'Stress Bodyforce 3', GotIt )
            END IF
          END IF

          IF( ListCheckPrefix(BodyForce,'Inertial Bodyforce') ) THEN
            InertialLoad(1,1:n) = GetReal( BodyForce, 'Inertial Bodyforce 1', GotIt )
            InertialLoad(2,1:n) = GetReal( BodyForce, 'Inertial Bodyforce 2', GotIt )
            IF ( dim > 2 ) THEN
              InertialLoad(3,1:n) = GetReal(  BodyForce, 'Inertial Bodyforce 3', GotIt )
            END IF
          END IF

          IF( STDOFS > dim ) THEN
            LoadVector(STDOFs,1:n) = GetReal( BodyForce, 'Stress Volume Source', GotIt )
          END IF

          ! "Stress Pressure": a pressure-like body load, which enters the force as
          ! p times the divergence of the test function rather than as a force
          ! density. Poroelasticity is what asks for it -- the ElmerIce permafrost
          ! case passes a groundwater pressure through it. NOT the same keyword as
          ! "Stress Volume Source" above, which feeds the constraint row of the
          ! mixed formulation.
          IF ( .NOT. EvalLoadIP ) NodalPressureLoad(1:n) = &
              GetReal( BodyForce, 'Stress Pressure', GotPressureLoad )

          !----------------------------------------------------------------------
          ! "Gravitational Prestress Advection": the restoring force that arises
          ! when a body under its own weight is displaced vertically, so that a
          ! column of material of density rho is advected through a gravity field g.
          ! One term in the stiffness, coupling every momentum row to the VERTICAL
          ! displacement, with rho*g given as "GPA Coeff". Read here and applied in
          ! LocalMatrix; it is not a constitutive property and does not go near the
          ! material model.
          !----------------------------------------------------------------------
          GotGPA = GetLogical( BodyForce, 'Gravitational Prestress Advection', GotIt )
          IF( GotGPA ) THEN
            NodalGPA(1:n) = GetReal( BodyForce, 'GPA Coeff', GotIt )
            IF( .NOT. GotIt ) CALL Warn( Caller, &
                '"Gravitational Prestress Advection" is set with no "GPA Coeff": '// &
                'the term is then identically zero' )
          END IF

          !----------------------------------------------------------------------
          ! An additive stress and an eigenstrain, both in Voigt form. Together
          ! they make the response affine rather than merely linear:
          !
          !   sigma = C : ( eps - "Strain Load" ) + "Stress Load"
          !
          ! which is StressSolve's reading of the same two keywords, and what
          ! LocalMatrix's offset channel implements below.
          !
          ! Zeroed before the test and not inside it: these arrays are SAVEd, so a
          ! body force that sets neither keyword would otherwise inherit whatever
          ! the previous element's body force did -- the trap LoadVector above is
          ! zeroed against for the same reason.
          !----------------------------------------------------------------------
          IF ( ListCheckPrefix( BodyForce, 'Stress Load' ) ) &
              CALL GetVoigtLoad( BodyForce, 'Stress Load', NodalStressLoad, n )
          IF ( ListCheckPrefix( BodyForce, 'Strain Load' ) ) &
              CALL GetVoigtLoad( BodyForce, 'Strain Load', NodalStrainLoad, n )
        END IF
                
        !------------------------------------------------------------------------------
        !        Get values of field variables:
        !------------------------------------------------------------------------------
        IF (UseUMAT) THEN
          LocalTemperature(1:n) = ReferenceTemperature(1:n) 
          IF ( ASSOCIATED(TempSol) ) THEN
            WHERE( TempPerm( NodeIndexes(1:n) ) > 0 )
              LocalTemperature(1:n) = Temperature(TempPerm(NodeIndexes(1:n)))
            END WHERE
          END IF
        ELSE
          LocalTemperature = 0.0D0
          IF ( ASSOCIATED(TempSol) ) THEN
            WHERE( TempPerm( NodeIndexes(1:n) ) > 0 )
              LocalTemperature(1:n) = Temperature(TempPerm(NodeIndexes(1:n))) - &
                  ReferenceTemperature(1:n)
            END WHERE
          END IF
        END IF

        LocalDisplacement = 0.0D0
        IF( .NOT. LinearModel ) THEN
          DO i=1,nd
            k = StressPerm(Indices(i))
            DO j=1,STDOFs
              LocalDisplacement(j,i) = Displacement(STDOFs*(k-1)+j)
            END DO
          END DO
        END IF
        
        ! ----------------------------------------------------------------
        ! Some material models may need the displacement field at the
        ! previous time/load step
        ! ----------------------------------------------------------------
        PrevLocalDisplacement = 0.0D0          
        IF( previ > 0 ) THEN
          DO i=1,nd
            k = StressPerm(Indices(i))
            DO j=1,STDOFs
              PrevLocalDisplacement(j,i) = Solver % Variable % PrevValues(STDOFs*(k-1)+j,previ)
            END DO
          END DO
        END IF
        
        !-------------------------------------------------------------------------------------------
        !        Select subroutine to integrate the element matrix and vector
        !-------------------------------------------------------------------------------------------
        IF (UseUMAT) THEN
          ! ------------------------------------------------------------------------------
          ! This branch assumes that the material behavior is defined 
          ! via an umat subroutine. The umat routine should specify
          ! a material response function which gives the Cauchy stress
          ! as a function of the strain tensor and state variables.
          !-------------------------------------------------------------------------------
          CALL LocalMatrixWithUMAT(LocalMassMatrix, LocalDampMatrix, &
              LocalStiffMatrix, LocalForce, LocalExternalForce, time, dt, LoadVector, InertialLoad, &
              MaterialConstants, NPROPS, NSTATEV, &
              InitializeStateVars, Density, Damping, AxialSymmetry, &
              PlaneStress, LargeDeflection, HenckyStrain, CurrentElement, n, nd, ntot, STDOFs, &
              ElementNodes, LocalDisplacement, PrevLocalDisplacement, LocalTemperature, &
              t, Iter, UMATModel)

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
                Isotropic, RotateModuli, TransformMatrix, LargeDeflection, &
                NodalStressLoad, NodalStrainLoad, LinearIncompressible, &
                MaxwellHere, MaxwellViscosity, PrevLocalDisplacement, GotGPA, NodalGPA, &
                GeometricActive, StabilityAnalysis, QuasiStationary, &
                GotPressureLoad, NodalPressureLoad, &
                EvalYoungIP, EvalPoissonIP, EvalBetaIP, EvalLoadIP )
          END IF
        END IF

        ! Staged construction: this body's displacement is measured from a stored
        ! reference rather than from the mesh. Between the local matrix and the
        ! update, since it needs the one to build a contribution to the other.
        IF ( UpdateReference ) CALL AddReferenceDisplacement( CurrentElement, n )

        IF( GotRayleighAlpha ) THEN
          LocalDampMatrix = LocalDampMatrix + RayleighAlpha * LocalMassMatrix
        END IF
        IF( GotRayleighBeta ) THEN
          LocalDampMatrix = LocalDampMatrix + RayleighBeta * LocalStiffMatrix
        END IF
        
        !------------------------------------------------------------------------------
        !        If time dependent simulation, add mass matrix to global 
        !        matrix and global RHS vector
        !------------------------------------------------------------------------------
        IF ( TransientSimulation ) THEN
           IF( SecondOrderTime ) THEN
             CALL Default2ndOrderTime( LocalMassMatrix, LocalDampMatrix, &
                  LocalStiffMatrix, LocalForce )
           ELSE
             CALL Default1stOrderTime( LocalDampMatrix, LocalStiffMatrix, LocalForce )
           END IF
        END IF
        !------------------------------------------------------------------------------
        !        Update global matrices from local matrices
        !------------------------------------------------------------------------------
200     IF ( ConstantBulkMatrixInUse ) THEN
          ! The matrix was restored wholesale, so only the load is wanted here.
          CALL DefaultUpdateForce( LocalForce )
        ELSE
          CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )
        END IF
        !------------------------------------------------------------------------------

        IF( NeedMass ) THEN
          CALL DefaultUpdateMass( LocalMassMatrix )
          IF( AnyDamping ) CALL DefaultUpdateDamp( LocalDampMatrix )
        END IF
        

     END DO

     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     !     Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
2000 DO t = 1,GetNOFBoundaryElements()
        CurrentElement =>  GetBoundaryElement(t)
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
           ! We assume that consistently either keyword type is used.
           !------------------------------------------------------------------------------
           GotForceBC = .FALSE.
           IF( ListCheckPrefix( BC,'Surface Traction' ) ) THEN
             GotForceBC = .TRUE.
             LoadVector(1,1:n) = GetReal( BC, 'Surface Traction 1', GotIt )
             LoadVector(2,1:n) = GetReal( BC, 'Surface Traction 2', GotIt )
             LoadVector(3,1:n) = GetReal( BC, 'Surface Traction 3', GotIt )
           ELSE IF( ListCheckPrefix( BC,'Force' ) ) THEN
             GotForceBC = .TRUE.
             LoadVector(1,1:n) = GetReal( BC, 'Force 1', GotIt )
             LoadVector(2,1:n) = GetReal( BC, 'Force 2', GotIt )
             LoadVector(3,1:n) = GetReal( BC, 'Force 3', GotIt )
           END IF
             
           Beta(1:n) = GetReal( BC, 'Normal Surface Traction', GotIt )
           IF (.NOT. GotIt) Beta(1:n) = GetReal( BC, 'Normal Force', gotIt )
           GotForceBC = GotForceBC .OR. GotIt

           ! The other way to drive a lumping load case: a pure force or a pure
           ! moment over the lumping boundary, in place of a prescribed
           ! displacement. It OVERWRITES the surface load read above, which is the
           ! routine's own doing -- it zeroes the array it is given -- and is what
           ! StressSolve does too. Nothing but a "Model Lumping Boundary" reaches
           ! this, so a sif that does not lump cannot see it.
           IF ( UseModelLumping .AND. .NOT. Lump % FixDisplacement ) THEN
             IF ( GetLogical( BC, 'Model Lumping Boundary', GotIt ) ) THEN
               ! ElementNodes still holds the last BULK element here, and the moment
               ! cases need the coordinates of THIS boundary element. Harmless to
               ! refresh: LocalBoundaryMatrix fetches its own nodes.
               CALL GetElementNodes( ElementNodes, CurrentElement )
               CALL ModelLumpingLoads( Lump, iter, ElementNodes, n, LoadVector )
               GotForceBC = .TRUE.
             END IF
           END IF

           GotSpring = ListCheckPrefix( BC,'Spring' )
           IF( GotSpring ) THEN           
             SpringCoeff(1:n,1,1) = GetReal( BC, 'Spring', NormalSpring )           
             IF ( .NOT. NormalSpring ) THEN
               DO i=1,dim
                 SpringCoeff(1:n,i,i) = GetReal( BC, ComponentName('Spring',i), GotIt)
               END DO
               DO i=1,dim
                 DO j=1,dim
                   IF (ListCheckPresent(BC,'Spring '//i2s(i)//i2s(j) )) &
                       SpringCoeff(1:n,i,j)=GetReal( BC, 'Spring '//i2s(i)//i2s(j), GotIt)
                 END DO
               END DO
             END IF
           END IF
             
           ! "Stress Load" as a BOUNDARY condition, which StressSolve reads and this
           ! solver implements only as a body force. Gated model-wide, so a sif
           ! without it anywhere pays one logical per boundary element.
           IF ( StressLoadInBC ) THEN
             IF ( ListCheckPrefix( BC, 'Stress Load' ) ) CALL Fatal( Caller, &
                 '"Stress Load" is implemented here as a body force and not yet as '// &
                 'a boundary condition' )
           END IF

           ! On the VALUE, not on the presence -- see the body force half for why.
           IF ( ImagLoadInBC ) THEN
             DO i = 1,4
               IF ( i <= 3 ) THEN
                 str = 'Force '//I2S(i)//' im'
               ELSE
                 str = 'Normal Force im'
               END IF
               IF ( .NOT. ListCheckPresent( BC, str ) ) CYCLE
               Work(1,1,1:n) = GetReal( BC, str, GotIt )
               IF ( ANY( Work(1,1,1:n) /= 0.0_dp ) ) CALL Fatal( Caller, &
                   'The imaginary part of a boundary load is not implemented here: '// &
                   'this solver assembles the real system only, so "'//TRIM(str)//'" '// &
                   'would be read and dropped' )
             END DO
           END IF

           GotFSIBC = GetLogical( BC, 'FSI BC', GotIt )

           IF ( .NOT. ( GotForceBC .OR. GotFSIBC .OR. GotSpring ) ) CYCLE

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
                   CompressibilityDefined = ( CompressibilityFlag /= 'incompressible' )  &
                       .OR. ( CompressibilityFlag /= 'artificial compressible') 
                 END IF
                 
                 DragCoeff = ListGetCReal( BC,'FSI Drag Multiplier',GotIt)
                 IF(GotIt) THEN
                   Viscosity(1:FlowNOFNodes) = DragCoeff * Viscosity(1:FlowNOFNodes) 
                 END IF

              END IF
           END IF

           NormalTangential = GetLogical( BC, 'Normal-Tangential ' // & 
                GetVarName(Solver % Variable), GotIt )

           CALL LocalBoundaryMatrix( LocalStiffMatrix, LocalForce, &
               LoadVector, SpringCoeff, GotSpring, NormalSpring, Alpha, Beta, LocalDisplacement, &
               CurrentElement, n, ntot, ParentElement, ParentElement % TYPE % NumberOfNodes, &
               nd, ParentNodes, FlowElement, FlowNOFNodes, FlowNodes, Velocity,  &
               Pressure, Viscosity, Density, CompressibilityDefined, AxialSymmetry, &
               NormalTangential, PseudoTraction, PressureUnknown, LargeDeflection)

           IF (UseUmat .AND. Iter == 1) THEN
             ! ---------------------------------------------------------------------------
             ! Update the RHS vector which contains just the contribution of external loads
             ! for the purpose of nonlinear error estimation:
             ! ---------------------------------------------------------------------------
             ValuesSaved => StiffMatrix % RHS
             StiffMatrix % RHS => StiffMatrix % BulkRHS
             LocalForceSaved = LocalForce
             CALL DefaultUpdateForce(LocalForce)
             LocalForce = LocalForceSaved
             Solver % Matrix % RHS => ValuesSaved
           END IF

           !------------------------------------------------------------------------------
           !           Update global matrices from local matrices (will also affect
           !           LocalStiffMatrix and LocalForce if transient simulation is on).
           !------------------------------------------------------------------------------

           IF ( TransientSimulation ) THEN
              LocalDampMatrix = 0._dp
              LocalMassMatrix = 0._dp

              IF( SecondOrderTime ) THEN
                CALL Default2ndOrderTime( LocalMassMatrix, LocalDampMatrix, &
                     LocalStiffMatrix, LocalForce )
              ELSE
                CALL Default1stOrderTime( LocalDampMatrix, LocalStiffMatrix, LocalForce )
              END IF
           END IF

           CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )              
        END IF
     END DO
     !------------------------------------------------------------------------------
     CALL DefaultFinishBoundaryAssembly()

     ! This is a matrix level routine for setting friction such that tangential
     ! traction is the normal traction multiplied by a coefficient.
3000 CALL SetImplicitFriction(Model, Solver,'Implicit Friction Coefficient',&
         'Friction Direction')
     
     CALL DefaultFinishAssembly()

     ! One load case imposed as a prescribed displacement of the lumping boundary --
     ! a pure translation or a pure rotation. Between the assembly and the Dirichlet
     ! conditions, which is where StressSolve puts it and what the routine's own
     ! comment requires: it sets boundary values directly into the matrix, and the
     ! sif's own conditions must be applied after so that they win where both speak.
     IF ( UseModelLumping .AND. Lump % FixDisplacement ) THEN
       CALL ModelLumpingDisplacements( Lump, Solver, Model, iter )
     END IF

     CALL DefaultDirichletBCs()

     IF (UseUMAT) THEN
       ! ---------------------------------------------------------------------------------
       ! Now the solution variable is the solution increment while the sif-file specifies
       ! the Dirichlet BCs for the complete field. Modify BCs so that the right BC
       ! is obtained for the solution increment.
       ! ---------------------------------------------------------------------------------
       DisplacementRot = Displacement
       CALL RotateNTSystemAll(DisplacementRot, StressPerm, STDOFs)
       
       IF (ALLOCATED(StiffMatrix % ConstrainedDOF)) THEN
         DO i=1,StiffMatrix % NumberOfRows
           IF (StiffMatrix % ConstrainedDOF(i)) THEN
             StiffMatrix % DValues(i) = StiffMatrix % DValues(i) - DisplacementRot(i)
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
       WRITE(Message,'(A,ES12.3)') 'Residual for nonlinear iterate '&
           //I2S(Iter-1)//': ',NonLinRes
       CALL Info('ElasticitySolver', Message, Level=5)        

       IF (NonlinRes < NonlinTol .AND. (iter-1) >= MinNonlinearIter) THEN
         CALL Info('ElasticitySolver','Nonlinear iteration is terminated succesfully',Level=5)

         ! Save the state variables corresponding to the converged nonlinear
         ! solution to the array holding the previous solution state:
         UmatEnergy0 = UmatEnergy
         UmatStress0 = UmatStress
         IF(ASSOCIATED(UmatState)) UmatState0 = UmatState
         
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

     ! One row or column of the lumped 6x6 stiffness, from the reactions this load
     ! case produced on the lumping boundary; after the sixth it is inverted and
     ! written out. Before the convergence tests below, not after: the sixth case
     ! leaves through one of those EXITs, and a row missed there is a matrix never
     ! written. The reactions come from BulkValues times the solution, which is why
     ! "Constant Bulk System" had to be added to the list and not merely assumed.
     IF ( UseModelLumping ) CALL ModelLumpingSprings( Lump, Solver, Model, iter )

     IF (UseUmat) THEN
       Displacement(:) = TotalSol(:) + Displacement(:)
       IF (iter==NonlinearIter) THEN
         CALL Info('ElasticitySolver', &
             'The maximum of nonlinear iterations reached: Terminating...', Level=5)        

         ! Save the state variables corresponding to the converged nonlinear
         ! solution to the array holding the previous solution state:
         UmatEnergy0 = UmatEnergy
         UmatStress0 = UmatStress
         IF(ASSOCIATED(UmatState)) UmatState0 = UmatState
         
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
  ! Staged construction: the bodies whose "Update Reference Displacement" came out
  ! non-negative take this step's solution as their new reference, so that the next
  ! stage is measured from where this one left off. Body by body, since a
  ! construction sequence has some bodies advancing while others are held.
  !-----------------------------------------------------------------------------
  IF ( UpdateReference ) THEN
    DO t = 1, Solver % NumberOfActiveElements
      CurrentElement => GetActiveElement(t)
      IF ( .NOT. UpdateActive( CurrentElement % BodyId ) ) CYCLE
      n = GetElementNOFNodes( CurrentElement )
      DO j = 1, n
        k = ReferenceSol % Perm( CurrentElement % NodeIndexes(j) )
        i = StressPerm( CurrentElement % NodeIndexes(j) )
        IF ( k == 0 .OR. i == 0 ) CYCLE
        ReferenceSol % Values( ReferenceSol % DOFs*(k-1)+1 : ReferenceSol % DOFs*k ) = &
            Displacement( ReferenceSol % DOFs*(i-1)+1 : ReferenceSol % DOFs*i )
      END DO
    END DO
  END IF

  !-----------------------------------------------------------------------------
  !   Perform strain and stress computation...
  !-----------------------------------------------------------------------------
  IF (CalculateStrains .OR. CalculateStresses) THEN
     CALL Info(Caller,'Computing postprocessing fields')

     !--------------------------------------------------------------------------
     ! An eigen or harmonic analysis has a displacement per mode, hence stresses
     ! per mode, so the postprocessing runs once for each. The nodal fields can
     ! only ever hold the last of them, so each mode's result is kept with the
     ! mode itself -- see ElasticityStoreEigenmode. A harmonic mode is complex and
     ! takes two passes, the real part and then the imaginary.
     !
     ! NOFEigenValues is zero in an ordinary solve, which is the single pass of
     ! the loop below with none of this entered.
     !
     ! The displacement is left holding the last mode rather than restored, which
     ! is what StressSolve does. The reported norm of a case like
     ! fem/tests/StrainCalculation03 is that mode's, so the two solvers have to
     ! agree here for its reference norm to survive them being merged. Displacing
     ! the mesh by a mode shape is the part that would be wrong, and the default
     ! set for "Displace Mesh" above stops that.
     !--------------------------------------------------------------------------
     EigenModes = Solver % NOFEigenValues
     Passes = 1
     IF ( EigenModes > 0 .AND. HarmonicAnalysis ) Passes = 2

     DO i=1,MAX( EigenModes, 1 )
        DO l=1,Passes
           IF ( EigenModes > 0 ) THEN
              CALL Info(Caller,'Computing stresses for eigenmode: '//I2S(i),Level=5)
              IF ( l == 1 ) THEN
                 Displacement = REAL( Solver % Variable % EigenVectors(i,:) )
              ELSE
                 Displacement = AIMAG( Solver % Variable % EigenVectors(i,:) )
              END IF
           END IF

           IF (UseUMAT) THEN
              CALL GenerateStressVariable(NodalStress, StressPerm, &
                  CalculateStresses, AxialSymmetry)

              CALL GenerateStrainVariable(Displacement, NodalStrain, StressPerm, CalculateStrains, &
                  AxialSymmetry, LargeDeflection)
           ELSE
              CALL ComputeStressAndStrain( Displacement, NodalStrain, NodalStress, VonMises, StressPerm, &
                   PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle, AxialSymmetry, NeoHookeanMaterial, &
                   CalculateStrains, CalculateStresses, CalcPrincipal, CalcPrincipalAngle, MixedFormulation, &
                   LargeDeflection, LinearIncompressible)
           END IF

           IF ( EigenModes > 0 ) CALL ElasticityStoreEigenmode( Solver, Mesh, i, l == 2 )
        END DO
     END DO
  END IF

  IF ( ListGetLogical(SolverParams, 'Adaptive Mesh Refinement', GotIt) ) THEN
     IF (UseUmat .OR. NeoHookeanMaterial) THEN
        CALL Info(Caller,'Adaptive Mesh Refinement is not available') 
     ELSE IF(.NOT.ListGetLogical(SolverParams, 'Library Adaptivity', GotIt ) ) THEN
        CALL RefineMesh( Model, Solver, Displacement, StressPerm, &
             ElasticSolver_Inside_Residual, ElasticSolver_Edge_Residual, ElasticSolver_Boundary_Residual )

        IF ( MeshDisplacementActive ) THEN
           StressSol => Solver % Variable
           IF ( .NOT.ASSOCIATED( Mesh, Model % Mesh ) ) &
                CALL DisplaceMesh( Mesh, StressSol % Values, 1, &
                StressSol % Perm, StressSol % DOFs, .FALSE. )
        END IF
     END IF
  END IF

  IF ( MeshDisplacementActive ) THEN
     CALL Info(Caller,'Displacing the mesh with computed displacement field')
     CALL DisplaceMesh( Mesh, Displacement, 1, StressPerm, STDOFs, .FALSE., dim )
  END IF

  DEALLOCATE( PrevSOL )
  IF (UseUmat) THEN
    DEALLOCATE(DisplacementRot)
    DEALLOCATE(LocalForceSaved)
  END IF

  CALL DefaultFinish()
  
  CALL Info('ElasticSolver','All done',Level=4)
  CALL Info('ElasticSolver','------------------------------------------',Level=4)

!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
!> Stop on a StressSolve keyword this solver does not implement, for the material
!> and body force of the element being assembled.
!>
!> Reached only when the model-wide gate says one of them is set somewhere, so the
!> string comparisons here are not on the ordinary path; see where that gate is
!> assigned. The lists tested are THIS element's, so a mixed sif -- some bodies
!> through StressSolve, some through here, which is what a migration looks like --
!> is not refused for what the other solver's bodies ask.
!>
!> Every entry is an item of the remaining unification backlog, and the message
!> says what is missing rather than merely that something is.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!> Staged construction, one element: decide from the sign of "Update Reference
!> Displacement" whether this body is to become the new reference or to be solved
!> measured from the existing one, and in the second case add K u_ref to the load.
!>
!> The sign convention is StressSolve's and it is a majority vote over the nodes,
!> not a test of one value -- a keyword given as a table over time crosses zero
!> between stages, and the nodes of one element can straddle the crossing.
!>
!> Only the NODAL block of the local matrix is contracted, which is what
!> StressSolve does: the reference is a nodal field and has nothing to say about
!> bubble degrees of freedom.
!------------------------------------------------------------------------------
  SUBROUTINE AddReferenceDisplacement( Element, n )
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BF
    INTEGER :: i, j, k, body, RefDofs, nref
    LOGICAL :: Found
    REAL(KIND=dp) :: UpdateRef(n), NodalRefD(STDOFs*n)
!------------------------------------------------------------------------------
    body = Element % BodyId
    IF ( body < 1 .OR. body > SIZE(UpdatePresent) ) RETURN
    IF ( .NOT. UpdatePresent(body) ) RETURN

    BF => GetBodyForce()
    IF ( .NOT. ASSOCIATED( BF ) ) RETURN

    UpdateRef(1:n) = GetReal( BF, 'Update Reference Displacement', Found )
    IF ( .NOT. Found ) RETURN

    ! Non-negative in the majority: this body's solution becomes the reference,
    ! which happens after the solve rather than here.
    IF ( COUNT( UpdateRef(1:n) < 0.0_dp ) <= COUNT( UpdateRef(1:n) >= 0.0_dp ) ) THEN
      UpdateActive(body) = .TRUE.
      RETURN
    END IF

    UpdateActive(body) = .FALSE.

    RefDofs = ReferenceSol % DOFs
    nref = RefDofs * n
    DO i = 1, n
      k = ReferenceSol % Perm( Element % NodeIndexes(i) )
      IF ( k == 0 ) THEN
        NodalRefD( RefDofs*(i-1)+1 : RefDofs*i ) = 0.0_dp
      ELSE
        NodalRefD( RefDofs*(i-1)+1 : RefDofs*i ) = &
            ReferenceSol % Values( RefDofs*(k-1)+1 : RefDofs*k )
      END IF
    END DO

    DO i = 1, nref
      LocalForce(i) = LocalForce(i) + &
          SUM( LocalStiffMatrix(i,1:nref) * NodalRefD(1:nref) )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddReferenceDisplacement
!------------------------------------------------------------------------------


  SUBROUTINE RefuseStressSolveKeywords( Material )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Material
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BF
    LOGICAL :: Found
    ! For the value test on the imaginary load below. N is the host's
    ! Mesh % MaxElementDOFs, so this covers any element.
    REAL(KIND=dp) :: Imag(N)
    INTEGER :: ic
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Material ) ) THEN
      IF ( ListCheckPrefix( Material, 'Pre Stress' ) .OR. &
           ListCheckPrefix( Material, 'Pre Strain' ) ) CALL Fatal( Caller, &
          '"Pre Stress" / "Pre Strain" are not implemented here. They enter '// &
          'StressSolve as a GEOMETRIC STIFFNESS and not as the additive stress '// &
          'that "Stress Load" / "Strain Load" give, which this solver does have' )

      ! The mesh velocity of a moving-mesh advection term, which StressSolve carries
      ! in its damping matrix. Group B in the plan document; no sif in the tree uses it.
      IF ( ListCheckPrefix( Material, 'Mesh Velocity' ) ) CALL Fatal( Caller, &
          '"Mesh Velocity" is not implemented here: StressSolve carries it as an '// &
          'advection term in the damping matrix, which this assembly does not build' )

    END IF

    BF => GetBodyForce()
    IF ( ASSOCIATED( BF ) ) THEN
      ! The imaginary half of a harmonic load. This solver assembles the real system
      ! only, so an imaginary body force would be read and dropped.
      !
      ! Refused on its VALUE and not on its presence, which matters here rather than
      ! being fastidious: HelmholtzStructure2 and 3 both declare
      ! "Stress Bodyforce 1 im = Real 0.0", and a zero imaginary load is no load at
      ! all -- nothing is dropped, so there is nothing to refuse. A refusal belongs on
      ! what would be lost, not on what was mentioned.
      DO ic = 1,3
        IF ( .NOT. ListCheckPresent( BF, 'Stress Bodyforce '//I2S(ic)//' im' ) ) CYCLE
        Imag(1:n) = GetReal( BF, 'Stress Bodyforce '//I2S(ic)//' im', Found )
        IF ( ANY( Imag(1:n) /= 0.0_dp ) ) CALL Fatal( Caller, &
            'The imaginary part of a body force is not implemented here: this solver '// &
            'assembles the real system only, so "Stress Bodyforce '//I2S(ic)//' im" '// &
            'would be read and dropped' )
      END DO

      ! "Stress Pressure" itself is implemented; only its imaginary half is not, and
      ! that on its value like the rest of the harmonic channel.
      IF ( ListCheckPresent( BF, 'Stress Pressure im' ) ) THEN
        Imag(1:n) = GetReal( BF, 'Stress Pressure im', Found )
        IF ( ANY( Imag(1:n) /= 0.0_dp ) ) CALL Fatal( Caller, &
            'The imaginary part of a body force is not implemented here: this solver '// &
            'assembles the real system only, so "Stress Pressure im" would be read '// &
            'and dropped' )
      END IF

    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE RefuseStressSolveKeywords
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Read a symmetric tensor given in Voigt form, nodewise, from a keyword list.
!>
!> Two spellings, because StressSolve accepts two and a sif written for it has to
!> keep working: the whole tensor as one array valued keyword,
!>
!>   Stress Load(6) = Real 1.0e5 0.0 0.0 0.0 0.0 0.0
!>
!> or one component at a time, "Stress Load 1" and so on. The array form wins
!> where both are given, which is the order StressSolve tries them in.
!>
!> ONE DELIBERATE DIFFERENCE from StressSolve: it accepts the componentwise form
!> for "Stress Load" only, and silently ignores "Strain Load 1". Here both
!> keywords take both spellings. That is a superset, so nothing written for
!> StressSolve changes meaning; the reverse direction is worth knowing about
!> before comparing the two solvers on a componentwise "Strain Load".
!------------------------------------------------------------------------------
  SUBROUTINE GetVoigtLoad( List, Name, Nodal, n )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: List
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp) :: Nodal(:,:)
    INTEGER :: n
!------------------------------------------------------------------------------
    LOGICAL :: Found
    INTEGER :: i, k
!------------------------------------------------------------------------------
    CALL GetRealArray( List, Work, Name, Found )
    IF ( Found ) THEN
      ! The first column of what may be given as a matrix, and at most the six
      ! independent components -- a longer keyword is the user's error, not a
      ! reason to write past the caller's array.
      k = MIN( SIZE(Work,1), 6 )
      Nodal(1:k,1:n) = Work(1:k,1,1:n)
      RETURN
    END IF

    DO i=1,6
      Nodal(i,1:n) = GetReal( List, TRIM(Name)//' '//I2S(i), Found )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetVoigtLoad
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! This subroutine uses the subroutine umat (Abaqus software convention for
! defining a user-supplied material model) to get the material model.
! This subroutine assumes that a stress response function for the Cauchy
! stress is supplied (originally Elmer has employed Piola-Kirchhoff stresses).
! A template subroutine UMAT_template located in the file 
!
!    .../fem/src/modules/UMATLib.F90) 
!
! provides a starting point for writing new user-supplied material models.
! An additional file which contains new UMAT material models can be named freely
! and it may contain several freely named subroutines that has the same arguments
! as the template subroutine UMAT_template. The Elmer solver keyword 
! "UMAT Subroutine" can be chosen to specify the file (that has been compiled
! with an elmerf90 command before simulation) and pick the subroutine desired. 
! NOTE: This is still a development version. For some examples see also
!       the directories .../fem/tests/UMAT_*
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixWithUMAT(MassMatrix, DampMatrix, StiffMatrix, ForceVector, &
       ExternalForceVector, time, dt, LoadVector, InertialLoad, MaterialConstants, &
       NrInProps, NStateV, &
       InitializeStateVars, NodalDensity, NodalDamping, AxialSymmetry, PlaneStress, &
       LargeDeflection, HenckyStrain, Element, n, nd, ntot, dofs, Nodes, NodalDisplacement, &
       PrevNodalDisplacement, NodalTemperature, ElementIndex, IterationIndex, &
       UMATModel)
    
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MassMatrix(:,:), DampMatrix(:,:), StiffMatrix(:,:)
    REAL(KIND=dp) :: ForceVector(:), ExternalForceVector(:)
    REAL(KIND=dp) :: time, dt
    REAL(KIND=dp) :: LoadVector(:,:), InertialLoad(:,:)
    REAL(KIND=dp), POINTER :: MaterialConstants(:,:)
    INTEGER :: NrInProps
    INTEGER :: NStateV
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
    CHARACTER(len=80) :: UMATModel
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
    INTEGER :: ipindex
    
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

    cmname = UMATModel

    dtime = dt
    TimeAtStep(1:2) = time - dt
    DRot = Identity

    ! ------------------------------------
    ! Integration stuff
    ! ------------------------------------   
    IntegStuff = GaussPoints( Element, RelOrder = RelIntegOrder )

    ForceVector = 0.0D0
    ExternalForceVector = 0.0D0
    StiffMatrix = 0.0D0
    MassMatrix  = 0.0D0
    DampMatrix  = 0.0d0

    
    DO t=1,IntegStuff % n
      ipindex = GetIpIndex( t, usolver=solver, element=element, ipvar = UmatEnergyVar )   
           
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
        ! Ordered (r, z, phi), hoop at 3, which is the UMAT interface's own
        ! (rr, zz, theta-theta, rz) and Elmer's convention everywhere outside
        ! this solver's assemblies. So the in-plane block is the same cdim
        ! product as the Cartesian branch below, with only the hoop appended.
        Grad(1:cdim,1:cdim) = MATMUL(NodalDisplacement(1:cdim,1:nd),dBasis(1:nd,1:cdim))
        Grad(3,3) = 1.0d0/r * SUM( NodalDisplacement(1,1:nd) * Basis(1:nd) )

        Grad0(1:cdim,1:cdim) = MATMUL(PrevNodalDisplacement(1:cdim,1:nd),dBasis(1:nd,1:cdim))
        Grad0(3,3) = 1.0d0/r * SUM( PrevNodalDisplacement(1,1:nd) * Basis(1:nd) )
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
          CALL Fatal( Caller, 'DSYEV cannot generate eigen basis')          
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
          CALL Fatal( Caller, 'DSYEV cannot generate eigen basis')          
        END IF

        Strain = 0.0d0
        Strain(1,1) = LOG(SQRT(EigenVals(1)))
        Strain(2,2) = LOG(SQRT(EigenVals(2)))       
        Strain(3,3) = LOG(SQRT(EigenVals(3)))
        Strain = MATMUL(QWork, MATMUL(Strain,TRANSPOSE(QWork)))

        ! NOTE: The differentiation of the Hencky strain is done via a truncated 
        ! series expansion which may become inaccurate for large strains. 
        ! However, this inaccuracy should not break the consistency
        ! of the solution method: if nonlinear iterations converge, we should have
        ! a solution. That is, inaccuracy of the strain expansion has the effect
        ! that the Newton iteration is replaced by inexact Newton iteration.
        ! Currently no warnings are given for the possibility that the strain
        ! expansion may not be accurate:
        !
        !IF ( ANY(EigenVals(:) >= 2.0d0) .OR. ANY(Eigenvals(:) <= 0.5d0) ) &
        !     CALL Fatal( Caller, 'Series expansion for Hencky strain too short!')
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
      ! With the axisymmetric ordering (r, z, phi) these slots are UMAT's own
      ! (11, 22, 33, 12) in both coordinate systems -- rr, zz, hoop, rz under
      ! axial symmetry and xx, yy, zz, xy in the plane -- so there is no
      ! axisymmetric special case left to write. Slots 5 and 6 are set beyond
      ! ntens = 4 there and simply not passed.
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

      ! -----------------------------------------------------------------------------
      ! Get the state variables and 
      ! the stress as specified at the previous time/load level for converged solution:
      ! -----------------------------------------------------------------------------
      EnergyElast = UmatEnergy0(3*(Ipindex-1)+1)
      EnergyPlast = UmatEnergy0(3*(Ipindex-1)+2)
      EnergyVisc = UmatEnergy0(3*(Ipindex-1)+3)

      StressVec(1:ntens) = UmatStress0(ntens*(Ipindex-1)+1:ntens*IpIndex)
      IF( NStateV > 0 ) THEN
        StateV(1:NstateV) = UmatState0(MaxStateV*(Ipindex-1)+1:MaxstateV*(IpIndex-1)+NStateV)
      END IF
        
      ! ----------------------------------------------------------------------------
      ! Obtain the Cauchy stress and the stress response function derivative 
      ! via UMAT interface. If the state variables have not been initiated to correspond
      ! the initial state (stress-free initial condition is supposed), we first make
      ! an extra UMAT call to obtain the state variables if requested in the sif file.
      ! ----------------------------------------------------------------------------
      INITIALIZE_STATE_VARIABLES: IF ( InitializeStateVars .AND. .NOT. UmatInitDone(ipIndex ) ) THEN
                
        ! We insert the identity tensor as the deformation gradient so the initial solution 
        ! should be the zero-displacement solution:
        stran = 0.0d0
        dstran = 0.0d0
                
        CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
            EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
            stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
            predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
            drot, pnewdt, celent, Identity, Identity, ElementIndex, t, layer, kspt, kstep, kinc)

        IF ( ANY(StressVec(1:ntens) /= UmatStress0(ntens*(Ipindex-1)+1:ntens*IpIndex) ) ) THEN
          CALL Fatal(Caller,'State variables initialization is changing stress')
        END IF

        ! Update the state variables storage (energy variables are not updated):
        IF( NStateV > 0 ) THEN
          UmatState0(MaxStateV*(Ipindex-1)+1:MaxstateV*(IpIndex-1)+NStateV) = StateV(1:NstateV)        
        END IF
          
        UmatInitDone(ipindex) = .TRUE.
      END IF INITIALIZE_STATE_VARIABLES

      ! -----------------------------------------------------------------------------
      ! Perform the actual UMAT call.
      ! -----------------------------------------------------------------------------      
      CALL UMATusersubrtn(UMATSubrtn, StressVec(1:ntens), StateV, StressDer(1:ntens,1:ntens), EnergyElast, &
          EnergyPlast, EnergyVisc, rpl, ddsddt(1:ntens), drplde(1:ntens), drpldt, &
          stran(1:ntens), dstran(1:ntens), TimeAtStep, dtime, Temp, dTemp, &
          predef, dpred, cmname, ndi, nshr, ntens, NStateV, InProps, NrInProps, coords, &
          drot, pnewdt, celent, DefG0, DefG, ElementIndex, t, layer, kspt, kstep, kinc)
        
      ! ---------------------------------------------------------------------------
      ! Update data which gives the state variables corresponding to the current 
      ! nonlinear iterate.
      ! ---------------------------------------------------------------------------
      UmatEnergy(3*(Ipindex-1)+1) = EnergyElast
      UmatEnergy(3*(Ipindex-1)+2) = EnergyPlast
      UmatEnergy(3*(Ipindex-1)+3) = EnergyVisc
      
      UmatStress(ntens*(Ipindex-1)+1:ntens*IpIndex) = StressVec(1:ntens)
      IF( NStateV > 0 ) THEN
        UmatState(MaxStateV*(Ipindex-1)+1:MaxstateV*(IpIndex-1)+NStateV) = StateV(1:NStateV) 
      END IF
        
      STIFFMATRIX_FOR_CHOSEN_STRAIN: IF (.NOT. LargeDeflection) THEN
        ! ----------------------------------------
        ! Create the strain-displacement matrix B:
        ! ----------------------------------------
        IF (AxialSymmetry) THEN
          ! Rows in UMAT's (rr, zz, theta-theta, rz) order. The hoop is row 3,
          ! not row 2: it used to sit in row 2 with zz in row 3, self-consistently
          ! with the old (r, phi, z) tensor ordering but swapped against the
          ! convention the user's own umat subroutine is written to.
          DO p=1,ntot
            B(1,(p-1)*dofs+1) = dBasis(p,1)
            B(2,(p-1)*dofs+2) = dBasis(p,2)
            B(3,(p-1)*dofs+1) = 1.0d0/r * Basis(p)
            B(4,(p-1)*dofs+1) = dBasis(p,2)
            B(4,(p-1)*dofs+2) = dBasis(p,1)
          END DO
        ELSE
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
        END IF

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
            StressVec(3)*SymBasis3
        SELECT CASE(nshr)
        CASE(1)
          ! SymBasis4 is the symmetric (1,2) basis, and slot 4 is the in-plane
          ! shear in both coordinate systems now that axial symmetry orders the
          ! components (r, z, phi): rz there, xy in the plane. The axisymmetric
          ! branch that mapped slot 4 onto SymBasis5, the (1,3) basis, went with
          ! the old ordering and is gone with it.
          Stress = Stress + 2.0d0*StressVec(4)*SymBasis4
        CASE(3)
          Stress = Stress + 2.0d0*StressVec(4)*SymBasis4 + &
              2.0d0*StressVec(5)*SymBasis5 + 2.0d0*StressVec(6)*SymBasis6
        END SELECT

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
                Grad(1,2) = dBasis(p,2)
                Grad(3,3) = 1.0d0/r * Basis(p)
              CASE (2)
                Grad(2,1) = dBasis(p,1)
                Grad(2,2) = dBasis(p,2)
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
            SELECT CASE(nshr)
            CASE(1)
              ! Slot 4 is the in-plane shear in both coordinate systems: see the
              ! note on the stress unpacking above.
              WorkVec1(4,1) = WorkTensor1(1,2) +  WorkTensor1(2,1)
            CASE(3)
              WorkVec1(4,1) = WorkTensor1(1,2) +  WorkTensor1(2,1)
              WorkVec1(5,1) = WorkTensor1(1,3) +  WorkTensor1(3,1)
              WorkVec1(6,1) = WorkTensor1(2,3) +  WorkTensor1(3,2)
            END SELECT

            WorkVec2(1:ntens,1) = MATMUL(TRANSPOSE(StressDer(1:ntens,1:ntens)),WorkVec1(1:ntens,1))
            ! The following is based on the assumed symmetry of StressDer(:,:):
            WorkTensor3 = WorkVec2(1,1)*SymBasis1 + WorkVec2(2,1)*SymBasis2 + &
                WorkVec2(3,1)*SymBasis3
            SELECT CASE(nshr)
            CASE(1)
              WorkTensor3 = WorkTensor3 + 2.0d0*WorkVec2(4,1)*SymBasis4
            CASE(3)
              WorkTensor3 = WorkTensor3 + 2.0d0*WorkVec2(4,1)*SymBasis4 + &
                  2.0d0*WorkVec2(5,1)*SymBasis5 + 2.0d0*WorkVec2(6,1)*SymBasis6
            END SELECT
            
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
                        + (dBasis(q,1)*dStress1(1,1) + dBasis(q,2)*dStress1(1,2) &
                        + 1.0d0/r*Basis(q)*dStress1(3,3))*s
                  CASE(2)
                    StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        = StiffMatrix(cdim*(p-1)+i,cdim*(q-1)+j) &
                        + (dBasis(q,1)*dStress1(2,1) + dBasis(q,2)*dStress1(2,2) ) * s
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
                    -(dBasis(p,1) * Stress1(1,1) + dBasis(p,2) * Stress1(1,2) &
                    + 1.0d0/r * Basis(p) * Stress1(3,3)) * s
              CASE(2)
                ForceVector(cdim*(p-1)+i) = ForceVector(cdim*(p-1)+i) &
                    -(dBasis(p,1) * Stress1(2,1) + dBasis(p,2) * Stress1(2,2)) * s
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
      IF( GotDamping ) THEN
        DO p = 1,ntot
          DO q = 1,ntot
            DO i = 1,cdim
              DampMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                  = DampMatrix(cdim*(p-1)+i,cdim*(q-1)+i) &
                  + Basis(p)*Basis(q)*Damping*s
            END DO
          END DO
        END DO
      END IF

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
  SUBROUTINE LocalMatrix( MassMatrix,DampMatrix,StiffMatrix,ForceVector, &
       LoadVector, InertialLoad, ElasticModulus, NodalPoisson, NodalDensity, NodalDamping, &
       AxialSymmetry,PlaneStress,NodalHeatExpansion, NodalTemperature, Element, n, ntot, &
       Nodes, LocalDisplacement, Isotropic, RotateModuli, TransformMatrix, &
       LargeDeflection, NodalStressLoad, NodalStrainLoad, LinearIncompressible, &
       MaxwellHere, NodalViscosity, PrevLocalDispl, GotGPA, NodalGPA, &
       GeometricActive, StabilityAnalysis, QuasiStationary, &
       GotPressureLoad, NodalPressureLoad, &
       EvalYoungIP, EvalPoissonIP, EvalBetaIP, EvalLoadIP )
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: StiffMatrix(:,:),MassMatrix(:,:),DampMatrix(:,:), &
         NodalHeatExpansion(:,:,:), ElasticModulus(:,:,:)
    REAL(KIND=dp) :: NodalTemperature(:),NodalDensity(:), &
         NodalDamping(:),LoadVector(:,:), InertialLoad(:,:)
    REAL(KIND=dp) :: NodalStressLoad(:,:), NodalStrainLoad(:,:)
    REAL(KIND=dp) :: LocalDisplacement(:,:), TransformMatrix(3,3)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalPoisson

    LOGICAL :: AxialSymmetry,PlaneStress, Isotropic, RotateModuli, LargeDeflection
    LOGICAL :: LinearIncompressible
    LOGICAL :: MaxwellHere, GotGPA
    LOGICAL :: GeometricActive, StabilityAnalysis, QuasiStationary, GotPressureLoad
    LOGICAL :: EvalYoungIP, EvalPoissonIP, EvalBetaIP, EvalLoadIP
    ! For the integration point reads, whose Found flag nothing here acts on: a
    ! keyword asked for at the points and absent gives zero, as it does at the nodes.
    LOGICAL :: Found
    ! Whether the inertial term is actually dropped at this element -- see the mass.
    LOGICAL :: QuasiActive
    REAL(KIND=dp) :: NodalPressureLoad(:)
    REAL(KIND=dp) :: NodalViscosity(:), PrevLocalDispl(:,:), NodalGPA(:)

    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes

    INTEGER :: n, ntot
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: Basis(ntot)
    REAL(KIND=dp) :: dBasisdx(ntot,3),SqrtElementMetric

    REAL(KIND=dp) :: Force(3), InertialForce(3), NodalLame1(n),NodalLame2(n),Density, &
         Damping,Lame1,Lame2
    REAL(KIND=dp) :: Grad(3,3),Identity(3,3),DetDefG,G(6,6)
    ! Required by the condensation, and unused here: the out-of-plane component it
    ! describes does no virtual work, so it concerns the postprocessor alone.
    REAL(KIND=dp) :: EzzC(3)
    REAL(KIND=dp) ::  DefG(3,3), Strain(3,3), Stress2(3,3), Stress1(3,3)

    ! The affine part of the response: a stress at zero strain and an eigenstrain
    ! the elastic law does not see. StressOffset is what the two collapse into.
    REAL(KIND=dp) :: StressLoad(6), StrainLoad(6), StressOffset(3,3), EigenStrain(3,3)
    LOGICAL :: GotOffset, GotVoigtLoad, NeedHeat

    ! Maxwell viscoelasticity at this integration point. xPhi is the relaxation
    ! factor 1/(1 + mu/eta dt) that scales the whole stiffness; LagStress is the
    ! history read from the previous step and NewLagStress what replaces it;
    ! ElasticStress is the UNRELAXED response, which the update needs and the
    ! assembly does not. Pres/Pres0 are the mixed formulation's pressure now and at
    ! the previous step, zero without it.
    REAL(KIND=dp) :: xPhi, ShearModulus, MaxVisc, MuDer
    REAL(KIND=dp) :: LagStress(3,3), NewLagStress(3,3), ElasticStress(3,3)
    REAL(KIND=dp) :: LagVec(6), Pres, Pres0
    REAL(KIND=dp) :: NodalVelo(3,ntot)
    INTEGER :: VeBase, VeIdx, nIP
    ! rho*g at this point, for the gravitational prestress advection term.
    REAL(KIND=dp) :: GPAatIP, PressureLoadAtIP
    REAL(KIND=dp) :: YoungAtIP, PoissonAtIP
    ! The stress the geometric stiffness is built from -- the constitutive response
    ! of the current iterate, before any relaxation or affine offset, which is the
    ! stress StressSolve's own separate LocalStress call for this term returns.
    REAL(KIND=dp) :: GeomStress(3,3), InnerProd

    REAL(KIND=dp) :: dDefG(3,3),dStress1(3,3)
    REAL(KIND=dp) :: dDefGU(3,3),dStrainU(3,3),dStress2U(3,3),dStress1U(3,3)

    ! The test function directions, one per (test function, component), and the
    ! constitutive response to all of them. Sized 3*ntot rather than cdim*ntot
    ! because an automatic array is dimensioned on entry to the routine, before
    ! cdim has been assigned; cdim is at most three, so the bound is safe and
    ! only the tail goes unused.
    REAL(KIND=dp) :: dDefGs(3,3,3*ntot), dStrains(3,3,3*ntot), dStresses(3,3,3*ntot)

    REAL(KIND=dp) :: Temperature, HeatExpansion(3,3)

    INTEGER :: i,j,k,l,p,q,t,dim,cdim,DOFs

    REAL(KIND=dp) :: s,u,v,w,r

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    INTEGER :: N_Integ

    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

    LOGICAL :: stat

    TYPE(MaterialPoint_t) :: MatPoint
    TYPE(MaterialResponse_t) :: MatResponse
    TYPE(MaterialState_t) :: MatState
    TYPE(MaterialModel_t) :: MatModel
    !> Sized for the largest model this routine selects, and passed as
    !> MatProps(1:nProps) so each model sees exactly its own layout.
    REAL(KIND=dp) :: MatProps(ANISOLIN_NPROPS), MatPropsT(ANISOLIN_NPROPS)
    INTEGER :: nProps
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

    !-------------------------------------------------------------------------
    ! The linear mixed formulation. The Poisson ratio is not read as a number
    ! near 1/2: it IS 1/2, so lambda is dropped altogether and the pressure
    ! unknown carries the volumetric part, while mu = E / (2(1+1/2)) = E/3.
    ! That is what StressSolve's own "Incompressible" branch assembles, term for
    ! term, and it is why the ratio in the material may be absent or absurd
    ! without changing the answer.
    !
    ! DOFs is the stride of the local matrix, one more than the displacement
    ! components, with the pressure last at each node.
    !-------------------------------------------------------------------------
    DOFs = cdim
    IF ( LinearIncompressible ) THEN
       DOFs = cdim + 1
       NodalLame1(1:n) = 0.0d0
       NodalLame2(1:n) = ElasticModulus(1,1,1:n) / 3.0d0
    END IF


    ForceVector = 0.0D0
    StiffMatrix = 0.0D0
    MassMatrix  = 0.0D0
    DampMatrix  = 0.0d0

    ! THE KINEMATIC identity, dim-diagonal, and deliberately not the one the
    ! constitutive models build for themselves. It does double duty here -- it is
    ! also DefG when the deflection is small -- so it belongs to the deformation
    ! and not to the stress.
    !
    ! The models' ReducedIdentity is CDim-diagonal plus the out-of-plane entry
    ! away from plane stress. The two agree in three dimensions, under plane
    ! stress and under axial symmetry, and differ in ONE case, 2D plane strain,
    ! where this one is diag(1,1,0) and theirs diag(1,1,1). Neither is wrong:
    ! the assembly wants only the components that do virtual work, since a plane
    ! weak form cannot see sigma_33, while the postprocessor wants everything
    ! reportable.
    !
    ! What the difference costs here is nothing, and that is checked rather than
    ! hoped for. It makes Stress2(3,3) nonzero in plane strain where it used to
    ! be zero, and that entry reaches the force vector and the stiffness only
    ! through the third row and column of DefG and dDefG -- both identically zero
    ! in the plane, this identity being what puts DefG(3,3) at zero. Bit-identity
    ! across fourteen probes confirms it.
    Identity = 0.0D0
    DO i = 1,dim
       Identity(i,i) = 1.0D0
    END DO

    !--------------------------------------------------------------------------
    ! ONE assembly, both materials. The two branches this replaced differed only
    ! in how the stress follows from the strain: Grad, DefG, DetDefG, the Gateaux
    ! terms, the Newton loop and the stiffness assembly were the same code
    ! written twice. That duplication had already cost two uninitialised reads, a
    ! non-conforming MATMUL that flang and gfortran resolved differently, and the
    ! whole axisymmetric anisotropic case -- the second copy simply never grew
    ! the hoop kinematics. Selecting the model here removes the copy rather than
    ! adding a third.
    !--------------------------------------------------------------------------
    IF (Isotropic) THEN
       MatModel = IsotropicLinearModel()
    ELSE
       MatModel = AnisotropicLinearModel()
    END IF
    nProps = MatModel % nProps

    ! The Gateaux shortcut below applies the law to a strain INCREMENT and takes
    ! the result for a directional derivative. That is a licence the model
    ! grants, not a property of the assembly, so it is asserted rather than
    ! assumed: a law that is not linear in its strain measure belongs in
    ! NeoHookeanLocalMatrix or LocalMatrixWithUMAT, which exist for that reason.
    IF ( .NOT. MatModel % StrainLinear ) CALL Fatal( Caller, 'Material model "'// &
        TRIM(MatModel % Name)//'" is not linear in the strain and needs its own assembly' )

    ! Axial symmetry with an anisotropic material used to be refused here, and the
    ! reason was the axis ORDERING: this assembly ordered the components
    ! (r, phi, z) with the hoop at index 2, while a matrix valued "Youngs Modulus",
    ! this solver's own postprocessor and StressSolve all use (r, z, phi) with the
    ! hoop at index 3. An isotropic law cannot see the difference, the trace and
    ! the identity being permutation blind; for an anisotropic C it is exactly
    ! C(2,2) against C(3,3).
    !
    ! The assemblies now order the components (r, z, phi) too, so the elasticity
    ! matrix is read in the same Voigt order it is written in -- (rr, zz, hoop,
    ! rz, z-hoop, r-hoop) -- and the anticipated permutation adapter turned out
    ! not to be needed at all. Aligning the conventions WAS the whole of the work.

    MatPoint % Dim = dim
    MatPoint % CDim = cdim
    MatPoint % PlaneStress = PlaneStress
    MatPoint % AxiSymmetric = AxialSymmetry
    MatPoint % Kinematics = MERGE( KINEMATICS_LARGE_DEFLECTION, &
        KINEMATICS_SMALL_STRAIN, LargeDeflection )

    ! Whether the affine channel is live at all, tested once per element rather
    ! than per integration point. Nothing set means not merely a zero offset but
    ! no offset arithmetic, so an ordinary run pays for three ANYs.
    !
    ! NodalTemperature is the difference from the reference temperature, so an
    ! isothermal body is all zeros here whatever its reference temperature was,
    ! and a material with no expansion coefficient contributes nothing either.
    GotVoigtLoad = ANY( NodalStressLoad(1:6,1:n) /= 0.0_dp ) .OR. &
                   ANY( NodalStrainLoad(1:6,1:n) /= 0.0_dp )
    NeedHeat = ANY( NodalTemperature(1:ntot) /= 0.0_dp ) .AND. &
               ANY( NodalHeatExpansion(1:3,1:3,1:n) /= 0.0_dp )
    GotOffset = GotVoigtLoad .OR. NeedHeat .OR. MaxwellHere

    !-------------------------------------------------------
    !    Integration stuff
    !-------------------------------------------------------
    IntegStuff = GaussPoints( element, RelOrder = RelIntegOrder )

    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

    !--------------------------------------------------------------------------
    ! MAXWELL: the element's slice of the integration point lag stress, and the
    ! history this step is measured from.
    !
    ! The lag stress that CONVERGED at the end of the previous timestep is the
    ! reference for the whole of this one, while the assembly overwrites the current
    ! values as it goes. So it is saved once here, at the first nonlinear iteration
    ! of the first coupled iteration, rather than tested for at every integration
    ! point -- which is where StressSolve saves it too.
    !
    ! The velocity that a non-Newtonian viscosity is a function of is the
    ! displacement increment over the step. Formed for the whole element here, as
    ! StressSolve forms it, so EffectiveViscosity sees the same nodal field.
    !--------------------------------------------------------------------------
    IF ( MaxwellHere ) THEN
       VeBase = VeStress % Perm( Element % ElementIndex )
       nIP = VeStress % Perm( Element % ElementIndex + 1 ) - VeBase
       IF ( nIP /= N_Integ ) CALL Fatal( Caller, 'Element '// &
            I2S(Element % ElementIndex)//' integrates over '//I2S(N_Integ)// &
            ' points but "ve_stress" was allocated '//I2S(nIP)//' of them' )

       IF ( GetNonlinIter() == 1 .AND. GetCoupledIter() == 1 ) &
            VeStress % PrevValues( nve*VeBase+1 : nve*(VeBase+nIP), 1 ) = &
            VeStress % Values( nve*VeBase+1 : nve*(VeBase+nIP) )

       NodalVelo = 0.0_dp
       DO i=1,cdim
          NodalVelo(i,1:ntot) = ( LocalDisplacement(i,1:ntot) - &
               PrevLocalDispl(i,1:ntot) ) / dt
       END DO
    END IF

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

       ! "Stress Bodyforce at IP": the same load, asked for at this point rather than
       ! interpolated from the nodes. The inertial load is not part of the keyword and
       ! stays nodal, as it is there.
       IF ( EvalLoadIP ) THEN
          DO i=1,cdim
             Force(i) = ListGetElementReal( LoadIP_h(i), Basis, GetCurrentElement(), Found, GaussPoint=t )
          END DO
       END IF

       ! Density and damping are properties of the material and not of whether it
       ! happens to be isotropic, and the mass and damping matrices below are
       ! assembled outside that branch. Interpolated here, in the common part, for
       ! that reason: assigned only inside the isotropic branch, as they were, they
       ! were READ UNINITIALISED for every anisotropic material -- by the mass matrix
       ! loop at the foot of this integration point, and by the inertial term of the
       ! anisotropic branch's own force vector.
       !
       ! The symptom was an intermittent NaN, which is what an uninitialised read
       ! looks like when the value is multiplied by something that is usually zero:
       ! stack garbage that happens to carry a NaN pattern propagates, garbage that
       ! does not is silently multiplied away. Valgrind put it at the gluing of the
       ! local matrix; three runs of the same case NaNed twice.
       Density = SUM( NodalDensity(1:n)*Basis(1:n) )
       Damping = SUM( NodalDamping(1:n)*Basis(1:n) )

       ! rho*g for the gravitational prestress advection, zero unless asked for.
       GPAatIP = 0.0_dp
       IF ( GotGPA ) GPAatIP = SUM( NodalGPA(1:n)*Basis(1:n) )

       ! The pressure-like body load, likewise -- and it travels with the body force
       ! when that is asked for at the integration points, which is where StressSolve
       ! reads it too: one keyword governs the four components together.
       PressureLoadAtIP = 0.0_dp
       IF ( EvalLoadIP ) THEN
          PressureLoadAtIP = ListGetElementReal( LoadIP_h(4), Basis, GetCurrentElement(), Found, &
              GaussPoint=t )
       ELSE IF ( GotPressureLoad ) THEN
          PressureLoadAtIP = SUM( NodalPressureLoad(1:n)*Basis(1:n) )
       END IF

       ! The thermal state at this point. NodalTemperature is already the
       ! DIFFERENCE from the reference temperature, taken where it is gathered.
       IF ( NeedHeat ) THEN
          Temperature = SUM( NodalTemperature(1:ntot)*Basis(1:ntot) )
          DO i=1,3
             DO j=1,3
                HeatExpansion(i,j) = SUM( NodalHeatExpansion(i,j,1:n)*Basis(1:n) )
             END DO
          END DO

          ! At the integration point instead, and on the DIAGONAL only -- which is all
          ! either solver reads of this coefficient, and all StressSolve fills from the
          ! same handle.
          IF ( EvalBetaIP ) THEN
             HeatExpansion = 0.0_dp
             DO i=1,3
                HeatExpansion(i,i) = ListGetElementReal( BetaIP_h, Basis, GetCurrentElement(), &
                    Found, GaussPoint=t )
             END DO
          END IF
       END IF

       !------------------------------------------------------------------------
       ! The material data the model reads, pre-evaluated here. A model left to
       ! look its own keywords up per integration point would spend a string
       ! comparison per call, one to two orders of magnitude more than the entire
       ! interface costs.
       !
       ! MatPropsT carries the ADJOINT of the same data, and it is not decoration.
       ! The tangent term below contracts the test function gradient with the
       ! adjoint of the elasticity tensor while the residual terms use the tensor
       ! itself, and for a C that is not symmetric those are different matrices.
       ! The anisotropic branch this replaces wrote it as TRANSPOSE(G) at its one
       ! call site, which sat INSIDE the p,i loop; here the transpose is taken
       ! once per integration point instead of 208 times per hex. A self-adjoint
       ! law simply hands over the same numbers twice.
       !------------------------------------------------------------------------
       IF (Isotropic) THEN
          !-------------------------------------------------
          ! Lame parameters at the integration point
          !------------------------------------------------
          Lame1 = SUM( NodalLame1(1:n)*Basis(1:n) )
          Lame2 = SUM( NodalLame2(1:n)*Basis(1:n) )

          !--------------------------------------------------------------------
          ! Material data given AT this integration point. Only the SOURCE of the
          ! two numbers changes: the Lame parameters follow from Young and Poisson
          ! by the same formulas as the nodal path above, plane stress and the
          ! incompressible override included, which is why they are rebuilt here
          ! rather than a second convention being introduced.
          !
          ! Either keyword may be given alone, so whichever is not asked for at the
          ! point is interpolated from the nodes as before.
          !--------------------------------------------------------------------
          IF ( EvalYoungIP .OR. EvalPoissonIP ) THEN
             IF ( EvalYoungIP ) THEN
                YoungAtIP = ListGetElementReal( YoungIP_h, Basis, GetCurrentElement(), Found, GaussPoint=t )
             ELSE
                YoungAtIP = SUM( ElasticModulus(1,1,1:n)*Basis(1:n) )
             END IF

             IF ( EvalPoissonIP ) THEN
                PoissonAtIP = ListGetElementReal( PoissonIP_h, Basis, GetCurrentElement(), Found, GaussPoint=t )
             ELSE
                PoissonAtIP = SUM( NodalPoisson(1:n)*Basis(1:n) )
             END IF

             IF ( PlaneStress ) THEN
                Lame1 = YoungAtIP * PoissonAtIP / ( 1.0d0 - PoissonAtIP**2 )
             ELSE
                Lame1 = YoungAtIP * PoissonAtIP / &
                    ( (1.0d0 + PoissonAtIP) * (1.0d0 - 2.0d0*PoissonAtIP) )
             END IF
             Lame2 = YoungAtIP / ( 2.0d0 * (1.0d0 + PoissonAtIP) )

             IF ( LinearIncompressible ) THEN
                Lame1 = 0.0d0
                Lame2 = YoungAtIP / 3.0d0
             END IF
          END IF

          MatProps(ISOLIN_LAME1) = Lame1
          MatProps(ISOLIN_LAME2) = Lame2
          MatPropsT(1:ISOLIN_NPROPS) = MatProps(1:ISOLIN_NPROPS)
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

          ! Plane strain confines the out-of-plane thermal expansion, and for an
          ! anisotropic material the stress that confinement produces has in-plane
          ! normal components. Folded into the in-plane coefficients here, which is
          ! where StressSolve folds it too -- and it has to happen BEFORE the
          ! condensation below, since it reads the out-of-plane couplings the
          ! condensation clears. Under plane stress the out-of-plane direction is
          ! free and there is nothing to fold; for an isotropic material neither
          ! solver folds anything, which is its own small inconsistency and left as
          ! it stands.
          IF ( NeedHeat .AND. dim == 2 .AND. .NOT. PlaneStress ) THEN
             HeatExpansion(1,1) = HeatExpansion(1,1) + HeatExpansion(3,3) * &
                 ( G(2,2)*G(1,3)-G(1,2)*G(2,3) ) / ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )
             HeatExpansion(2,2) = HeatExpansion(2,2) + HeatExpansion(3,3) * &
                 ( G(1,1)*G(2,3)-G(1,2)*G(1,3) ) / ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )
          END IF

          ! Reduced to the plane packing that a two-dimensional contraction
          ! expects. Handed the raw 6x6, a plane assembly reads C(3,3) -- the 33
          ! modulus -- as the shear modulus: 1346 in place of 385 on the test
          ! material, and wrong in the stiffness rather than only in the output.
          IF ( dim == 2 ) &
              CALL CondensePlaneElasticityMatrix( G, PlaneStress, EzzC )

          ! Flattened in Fortran's own column major order, so the model indexes
          ! back with 6*(j-1)+i and no packing convention has to be agreed twice.
          MatProps(ANISOLIN_C:ANISOLIN_C+ANISOLIN_NPROPS-1) = &
              RESHAPE( G, [ ANISOLIN_NPROPS ] )
          MatPropsT(ANISOLIN_C:ANISOLIN_C+ANISOLIN_NPROPS-1) = &
              RESHAPE( TRANSPOSE(G), [ ANISOLIN_NPROPS ] )
       END IF

       !------------------------------------------------------------------------
       ! The affine part of the response, gathered into ONE stress offset:
       !
       !   sigma = C : ( eps - eps0 ) + sigma0
       !         = C : eps + ( sigma0 - C : eps0 )
       !
       ! where eps0 is "Strain Load" and sigma0 is "Stress Load". Both are
       ! properties of this point, not of any test function, so the parenthesis is
       ! evaluated once per integration point.
       !
       ! THERMAL EXPANSION IS THE SAME CHANNEL, and it is here rather than in a
       ! term of its own for that reason: alpha * dT is an eigenstrain, so it adds
       ! to eps0 and the arithmetic below carries it. StressSolve keeps them apart
       ! -- its thermal term contracts G*C with the expansion coefficient in the
       ! force vector, its "Strain Load" multiplies C by hand first -- but they are
       ! one capability, and one of them was already implemented here as nothing at
       ! all: this routine took NodalHeatExpansion and NodalTemperature as
       ! arguments and read neither. Thermal strain was silently dropped, and no
       ! test in the tree ran ElasticSolve with a heat expansion coefficient.
       !
       ! THE EIGENSTRAIN GOES THROUGH THE MODEL. C : eps0 could be formed from the
       ! elasticity matrix directly, and for the isotropic law that is two lines;
       ! asking the model contracts it with the same code that contracts every
       ! other strain here, so an anisotropic C, a condensed plane C and whatever
       ! law comes next need no second path and no second convention.
       !
       ! WHY THIS NEEDS NO NEW ASSEMBLY ROUTINE. An affine law is not linear in
       ! the strain, and this assembly's Gateaux shortcut requires linearity --
       ! MatModel % StrainLinear is asserted above for exactly that reason. The
       ! offset survives it because of WHERE it is added: to Stress2 alone, below,
       ! while every tangent term applies the law to a strain INCREMENT and so
       ! excludes the offset by construction. Adding it inside the model instead
       ! would put it into the derivative too, and the stiffness would be wrong.
       !------------------------------------------------------------------------
       IF ( GotOffset ) THEN
          StressLoad = 0.0_dp
          StrainLoad = 0.0_dp
          IF ( GotVoigtLoad ) THEN
             DO i=1,6
                StressLoad(i) = SUM( NodalStressLoad(i,1:n)*Basis(1:n) )
                StrainLoad(i) = SUM( NodalStrainLoad(i,1:n)*Basis(1:n) )
             END DO
          END IF

          ! Voigt to tensor in the packing this configuration uses: the full six
          ! in 3D, the reduced (11,22,12) in the plane -- where slot three is the
          ! shear and not the out-of-plane normal -- and (rr,zz,hoop,rz) under
          ! axial symmetry. StressLocal's own converter, so the packing is agreed
          ! in one place and CDim is what selects it, the same argument
          ! StressSolve passes.
          CALL Vector62Tensor( StressLoad, StressOffset, cdim, AxialSymmetry )
          CALL Vector62Tensor( StrainLoad, EigenStrain, cdim, AxialSymmetry )

          ! Halved off the diagonal. A Voigt STRAIN vector carries engineering
          ! shear, twice the tensor component, and the models double it back when
          ! they pack for C -- so the round trip is exact, two being a power of
          ! two. It is also invisible: a factor of two lost here changes no
          ! isotropic diagonal case and no test that does not shear.
          DO i=1,3
             DO j=1,3
                IF ( i /= j ) EigenStrain(i,j) = EigenStrain(i,j) / 2.0_dp
             END DO
          END DO

          ! The thermal eigenstrain, alpha * dT, diagonal and unsheared -- so the
          ! halving above does not concern it and it is added after. Over Dim and
          ! not over three: Dim is already the dimension of the state of STRESS, so
          ! it is three under axial symmetry, where the hoop expands and carries
          ! stress, and two in the plane, where the out-of-plane expansion is not
          ! part of the plane system. The confined plane strain case is the
          ! exception and it has been dealt with above, by folding the out-of-plane
          ! coefficient into the in-plane ones.
          IF ( NeedHeat ) THEN
             DO i=1,dim
                EigenStrain(i,i) = EigenStrain(i,i) + HeatExpansion(i,i) * Temperature
             END DO
          END IF

          MatPoint % Strain = EigenStrain
          CALL MatModel % Stress( MatPoint, MatProps(1:nProps), MatState, MatResponse )
          StressOffset = StressOffset - MatResponse % Stress

          !-------------------------------------------------------------------
          ! MAXWELL: the lag stress carried over from the previous step is an
          ! additive stress too, so it goes through this same channel. The
          ! relaxation factor is formed here because the offset needs it as well
          ! as the stiffness does.
          !
          ! THE SIGN IS WHERE THE TWO SOLVERS PART, and it is not arbitrary:
          ! StressSolve keeps the NEGATIVE of an additive offset in its own
          ! StressLoad -- "StressLoad = MATMUL(C,StrainLoad) - StressLoad" is what
          ! establishes that convention -- and adds it to the force vector, where
          ! this assembly adds the offset to the stress and the residual then
          ! subtracts it. So what goes in here is minus the tensor
          ! ViscoElasticLoad hands back.
          !-------------------------------------------------------------------
          IF ( MaxwellHere ) THEN
             MaxVisc = SUM( NodalViscosity(1:n) * Basis(1:n) )
             MaxVisc = EffectiveViscosity( MaxVisc, Density, NodalVelo(1,:), &
                  NodalVelo(2,:), NodalVelo(3,:), GetCurrentElement(), Nodes, &
                  n, ntot, u, v, w, MuDer, LocalIP=t )

             ! The shear modulus IS the second Lame parameter, in both
             ! formulations: E/(2(1+nu)) where the ratio is read, and E/3 where the
             ! incompressible branch set it above. StressSolve writes the two cases
             ! out as separate expressions and they come to the same number.
             ShearModulus = Lame2
             xPhi = 1.0_dp / ( 1.0_dp + ShearModulus / MaxVisc * dt )

             ! The mixed formulation's pressure now and at the previous step. Over
             ! the corner nodes only, the pressure living on the lowest-order basis.
             Pres = 0.0_dp
             Pres0 = 0.0_dp
             IF ( LinearIncompressible ) THEN
                Pres  = SUM( Basis(1:n) * LocalDisplacement(DOFs,1:n) )
                Pres0 = SUM( Basis(1:n) * PrevLocalDispl(DOFs,1:n) )
             END IF

             VeIdx = nve * ( VeBase + t - 1 )
             CALL Vector62Tensor( VeStress % PrevValues(VeIdx+1:VeIdx+nve,1), &
                  LagStress, cdim, AxialSymmetry )

             StressOffset = StressOffset - xPhi * ( LagStress - Pres0*Identity )
          END IF
       END IF

       !------------------------------------------------------------------
       ! Deformation gradient etc. evaluated using the current solution:
       !------------------------------------------------------------------
       Grad = 0.0d0
       IF (AxialSymmetry) THEN
          ! Ordered (r, z, phi), so the hoop is index 3. This is Elmer's
          ! axisymmetric convention everywhere else -- ComputeStressAndStrain
          ! below, StressSolve's LocalStress, and the UMAT interface's own
          ! (rr, zz, theta-theta, rz) -- and this assembly used to be the one
          ! place ordering them (r, phi, z) with the hoop at index 2. Invisible
          ! for an isotropic law, since the trace and the identity are blind to
          ! which axis is which; for an anisotropic C it was the difference
          ! between C(2,2) and C(3,3), and it is why anisotropy could not come
          ! near axial symmetry until the two orderings were reconciled.
          Grad(1,1) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(1,2) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,2) )
          Grad(3,3) = 1.0d0/r * SUM( LocalDisplacement(1,1:ntot) * Basis(1:ntot) )
          Grad(2,1) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(2,2) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,2) )
       ELSE
          Grad(1:dim,1:dim) = MATMUL(LocalDisplacement(1:dim,1:ntot),dBasisdx(1:ntot,1:dim))
       END IF
       ! Small strain keeps the reference and current configurations
       ! coincident, which also makes DetDefG come out as one below.
       IF (LargeDeflection) THEN
          DefG = Identity + Grad
       ELSE
          DefG = Identity
       END IF
       Strain = (TRANSPOSE(Grad)+Grad)/2.0D0
       IF (LargeDeflection) Strain = Strain + MATMUL(TRANSPOSE(Grad),Grad)/2.0D0

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
       ! Response % Stress is INTENT(OUT) with a default initialiser, so the
       ! whole tensor is written on every call. That subsumes the zeroing this
       ! branch used to need: Strain2Stress with dim 2 wrote only the four
       ! in-plane entries and left row and column three carrying stack garbage,
       ! which reached the force vector through a dBasisdx(p,3) that is zero in
       ! the plane -- silent unless the garbage happened to be a NaN.
       MatPoint % Strain = Strain
       CALL MatModel % Stress( MatPoint, MatProps(1:nProps), MatState, MatResponse )
       Stress2 = MatResponse % Stress

       !----------------------------------------------------------------------
       ! MAXWELL: the whole response is relaxed by xPhi. Scaling the STRESS rather
       ! than the moduli the model was handed is what keeps the unrelaxed response
       ! available -- the history update below needs C : eps itself, and dividing it
       ! back out of a scaled stress would be both wasteful and inexact.
       !----------------------------------------------------------------------
       ! The stress the geometric stiffness is built from: the constitutive response
       ! at this iterate, taken before the relaxation below and before the affine
       ! offset, which is what StressSolve's separate LocalStress call for this term
       ! gives it.
       GeomStress = Stress2

       IF ( MaxwellHere ) THEN
          ElasticStress = Stress2
          Stress2 = xPhi * Stress2
       END IF

       ! The affine offset, and the ONE place it may be added: the residual stress
       ! and nothing that a derivative is taken of. See where it is formed above.
       IF ( GotOffset ) Stress2 = Stress2 + StressOffset

       !----------------------------------------------------------------------
       ! MAXWELL: advance the lag stress at this point. Written into the variable's
       ! current values, from which the next timestep's first iteration will save
       ! the converged result as its own history -- so the last nonlinear iteration
       ! to run is the one that counts, which is why the scheme asks for at least
       ! two of them.
       !----------------------------------------------------------------------
       IF ( MaxwellHere ) THEN
          NewLagStress = ( 1.0_dp - xPhi ) * ElasticStress + &
               xPhi * ( LagStress - Pres0*Identity ) + Pres*Identity
          CALL Tensor26Vector( NewLagStress, LagVec, cdim, AxialSymmetry )
          VeStress % Values(VeIdx+1:VeIdx+nve) = LagVec(1:nve)
       END IF

       !--------------------------------------------------
       ! The first Piola-Kirchhoff stress
       !--------------------------------------------------
       Stress1 = MATMUL(DefG,Stress2)

       !-----------------------------------------------------------------------
       !  Gateaux derivatives of the solution with respect to the displacement:
       !  ---------------------------------------------------------------------
       dDefGU = Grad
       dStrainU = (MATMUL(TRANSPOSE(DefG),dDefGU) &
            + MATMUL(TRANSPOSE(dDefGU),DefG))/2.0D0
       MatPoint % Strain = dStrainU
       CALL MatModel % Stress( MatPoint, MatProps(1:nProps), MatState, MatResponse )
       dStress2U = MatResponse % Stress
       ! Relaxed with the same factor as the stress it is the derivative of.
       IF ( MaxwellHere ) dStress2U = xPhi * dStress2U
       dStress1U = MATMUL(DefG,dStress2U)
       IF (LargeDeflection) dStress1U = dStress1U + MATMUL(dDefGU,Stress2)

       !------------------------------------------------------------------------
       !  Gateaux derivatives of the solution with respect to the test functions,
       !  for every (test function, component) at this integration point.
       !
       !  Kept apart from the assembly below so that the constitutive law is
       !  applied to all of them in ONE call. Those cdim*ntot evaluations -- 24
       !  per integration point on a trilinear hexahedron, 192 of the element's
       !  208 -- share this point, these Props and this State; only the strain
       !  differs, which is exactly the shape the batched entry point is for.
       !  One call per point rather than cdim*ntot of them is worth 21 of the 25
       !  percent the interface had cost this assembly, and what it saves is the
       !  mechanism rather than the arithmetic: see ConstitutiveStressBatch_i.
       !------------------------------------------------------------------------
       DO p = 1,ntot
          DO i = 1,cdim
             k = cdim*(p-1)+i
             dDefGs(:,:,k) = 0.0D0
             IF (AxialSymmetry) THEN
                ! (r, z, phi), hoop at 3, as Grad above.
                SELECT CASE(i)
                CASE (1)
                   dDefGs(1,1,k) = dBasisdx(p,1)
                   dDefGs(1,2,k) = dBasisdx(p,2)
                   dDefGs(3,3,k) = 1.0d0/r * Basis(p)
                CASE (2)
                   dDefGs(2,1,k) = dBasisdx(p,1)
                   dDefGs(2,2,k) = dBasisdx(p,2)
                END SELECT
             ELSE
                dDefGs(i,:,k) = dBasisdx(p,:)
             END IF

             dStrains(:,:,k) = (MATMUL(TRANSPOSE(DefG),dDefGs(:,:,k)) &
                  + MATMUL(TRANSPOSE(dDefGs(:,:,k)),DefG))/2.0D0
          END DO
       END DO

       ! The adjoint of the tensor, not the tensor: see MatPropsT above.
       CALL ConstitutiveStresses( MatModel, MatPoint, MatPropsT(1:nProps), &
            MatState, cdim*ntot, dStrains, dStresses )

       ! MAXWELL: the stiffness is relaxed by the same factor as the stress, which
       ! is the whole of what the rheology does to the tangent -- the lag stress
       ! enters the residual and not this.
       IF ( MaxwellHere ) dStresses(:,:,1:cdim*ntot) = xPhi * dStresses(:,:,1:cdim*ntot)

       !----------------------------------------------------------------------------
       ! Loop over the test functions (stiffness matrix for Newton linearization):
       ! ---------------------------------------------------------------------------
       DO p = 1,ntot
          DO i = 1,cdim
             k = cdim*(p-1)+i
             dDefG = dDefGs(:,:,k)
             dStress1 = MATMUL(DefG,dStresses(:,:,k))
             IF (LargeDeflection) dStress1 = dStress1 + MATMUL(dDefG,Stress2)

             IF (AxialSymmetry) THEN

                ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                     +(Basis(p)*Force(i)*DetDefG &
                     +Basis(p)*InertialForce(i)*Density &
                     +PressureLoadAtIP*dBasisdx(p,i) &
                     -DDOTPROD(dDefG,Stress1,dim) &
                     +DDOTPROD(dDefG,dStress1U,dim))*s

                DO q = 1,ntot
                   DO j = 1,cdim
                      SELECT CASE(j)
                      CASE(1)
                         ! (r, z, phi): the radial row, the r-z shear and the hoop.
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(1,1) + dBasisdx(q,2)*dStress1(1,2) &
                              + 1.0d0/r*Basis(q)*dStress1(3,3))*s
                      CASE(2)
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(2,1) + dBasisdx(q,2)*dStress1(2,2) ) * s
                      END SELECT
                   END DO
                END DO

             ELSE

                ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                     +(Basis(p)*Force(i)*DetDefG &
                     +Basis(p)*InertialForce(i)*Density &
                     +PressureLoadAtIP*dBasisdx(p,i) &
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

             !--------------------------------------------------------------------
             ! The pressure column of this momentum row: -p div v. The pressure is
             ! NOT put into the stress here, the way the affine offset of "Stress
             ! Load" is, because it is an UNKNOWN of the same system -- it has to
             ! reach the matrix, not the residual, or the monolithic solve becomes
             ! a Picard iteration on the pressure and the sifs that ask for one
             ! nonlinear iteration would return a half-solved system.
             !--------------------------------------------------------------------
             IF ( LinearIncompressible ) THEN
                DO q = 1,ntot
                   StiffMatrix(DOFs*(p-1)+i,DOFs*q) = StiffMatrix(DOFs*(p-1)+i,DOFs*q) &
                        - Basis(q)*dBasisdx(p,i)*s
                END DO
             END IF

             !--------------------------------------------------------------------
             ! The gravitational prestress advection term: rho*g times the gradient
             ! of the trial function against the test function's value, in the
             ! column of the VERTICAL displacement -- the last coordinate component.
             ! A stiffness and nothing else, so it needs no counterpart in the
             ! residual: the two residual terms of this assembly cancel for a linear
             ! law, leaving the load, and an addition to the matrix alone is
             ! therefore exactly an addition to the operator.
             !
             ! ADDED AS ITS SYMMETRIC PART, which is not a modelling choice made
             ! here but a reproduction of one made there. StressSolve symmetrises
             ! its whole local stiffness wholesale before gluing it in --
             ! "STIFF = (STIFF + TRANSPOSE(STIFF))/2" -- and every other term it
             ! assembles is already symmetric, so the operator that reaches its
             ! linear solver contains only the symmetric half of this one. Assembling
             ! the term as written instead moves the earth cases 2.6%, and their
             ! published reference norms are the symmetric-half ones.
             !
             ! THIS TERM IS GENUINELY NOT SYMMETRIC, and whether the symmetric half
             ! is the intended model is an OPEN QUESTION, not something this comment
             ! settles: it couples the gradient of the VERTICAL displacement to every
             ! test component, and its transpose is a different operator. Nothing in
             ! the sifs that use it forces the choice either -- they solve with GCR
             ! and do not declare "Linear System Symmetric", so they could carry a
             ! non-symmetric operator as they stand. So this reproduces StressSolve
             ! because the merger requires it to, and if the answer is that the full
             ! term belongs in the system, both solvers change together and the earth
             ! reference norms move with them.
             !--------------------------------------------------------------------
             IF ( GotGPA ) THEN
                DO q = 1,ntot
                   StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+cdim) = &
                        StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+cdim) &
                        + 0.5_dp*GPAatIP*dBasisdx(q,i)*Basis(p)*s
                   StiffMatrix(DOFs*(q-1)+cdim,DOFs*(p-1)+i) = &
                        StiffMatrix(DOFs*(q-1)+cdim,DOFs*(p-1)+i) &
                        + 0.5_dp*GPAatIP*dBasisdx(q,i)*Basis(p)*s
                END DO
             END IF
          END DO

          !-----------------------------------------------------------------------
          ! The constraint equation, once per test function rather than per
          ! component: -div u tested with the pressure basis. Its residual is
          ! identically zero, so there is nothing to add to ForceVector.
          !-----------------------------------------------------------------------
          IF ( LinearIncompressible ) THEN
             DO q = 1,ntot
                DO j = 1,cdim
                   StiffMatrix(DOFs*p,DOFs*(q-1)+j) = StiffMatrix(DOFs*p,DOFs*(q-1)+j) &
                        - dBasisdx(q,j)*Basis(p)*s
                END DO
             END DO
          END IF
       END DO


       !      Integrate mass matrix:
       !      ----------------------
       ! Not under stability analysis: there the mass slot carries the geometric
       ! stiffness instead, and a density mass added to it would turn the buckling
       ! eigenproblem into something that is neither buckling nor vibration.
       !
       ! Nor under "Quasi Stationary", which asks for the transient run without its
       ! inertial term. StressSolve conflates the two matrices in that one flag --
       ! it clears NeedMass from the keyword and then lets any damping present set it
       ! again, which restores BOTH -- so the test here is the one it makes, on this
       ! element's own nodal damping rather than on the keyword being mentioned.
       QuasiActive = QuasiStationary .AND. .NOT. ( ANY( NodalDamping(1:n) /= 0.0_dp ) &
           .OR. GotRayleighAlpha .OR. GotRayleighBeta )

       IF ( .NOT. ( StabilityAnalysis .OR. QuasiActive ) ) THEN
         DO p = 1,ntot
           DO q = 1,ntot
             DO i = 1,cdim
               MassMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   = MassMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   + Basis(p)*Basis(q)*Density*s
             END DO
           END DO
         END DO
       END IF

       !------------------------------------------------------------------------
       ! THE GEOMETRIC STIFFNESS of the state this iterate is in:
       !
       !   INT sigma_ij d(phi_p)/dx_i d(phi_q)/dx_j
       !
       ! diagonal in the displacement components, since it is the work done by the
       ! existing stress against the rotation of a fibre and that is the same for
       ! each component. Live on the second pass only -- see the loop that switches
       ! it on -- and into the stiffness or, negated, into the mass, which is what
       ! separates a prestressed vibration problem from a buckling one.
       !------------------------------------------------------------------------
       IF ( GeometricActive ) THEN
         DO p = 1,ntot
           DO q = 1,ntot
             InnerProd = 0.0d0
             DO i = 1,dim
               DO j = 1,dim
                 InnerProd = InnerProd + dBasisdx(p,i)*dBasisdx(q,j)*GeomStress(i,j)
               END DO
             END DO

             DO k = 1,cdim
               IF ( StabilityAnalysis ) THEN
                 MassMatrix(DOFs*(p-1)+k,DOFs*(q-1)+k) &
                     = MassMatrix(DOFs*(p-1)+k,DOFs*(q-1)+k) - s * InnerProd
               ELSE
                 StiffMatrix(DOFs*(p-1)+k,DOFs*(q-1)+k) &
                     = StiffMatrix(DOFs*(p-1)+k,DOFs*(q-1)+k) + s * InnerProd
               END IF
             END DO
           END DO
         END DO
       END IF

       !      Utilize the nodal damping:
       !      -----------------------------
       IF( GotDamping .AND. .NOT. QuasiActive ) THEN
         DO p = 1,ntot
           DO q = 1,ntot
             DO i = 1,cdim
               DampMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   = DampMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   + Basis(p)*Basis(q)*Damping*s
             END DO
           END DO
         END DO
       END IF

    END DO

    !--------------------------------------------------------------------------
    ! The pressure lives on the lowest-order basis only, so the bubble and other
    ! higher-order pressure degrees of freedom are eliminated with a unit
    ! diagonal -- the MINI element, and the same elimination StressSolve writes.
    ! It is done after the integration loop because it clears whole rows and
    ! columns, which an integration point may not do.
    !--------------------------------------------------------------------------
    IF ( LinearIncompressible ) THEN
       DO p = n+1,ntot
          i = DOFs*p
          ForceVector(i)   = 0.0d0
          StiffMatrix(i,:) = 0.0d0
          StiffMatrix(:,i) = 0.0d0
          StiffMatrix(i,i) = 1.0d0
       END DO
    END IF
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

    INTEGER :: n, ntot
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    REAL(KIND=dp) :: Basis(ntot)
    REAL(KIND=dp) :: dBasisdx(ntot,3),SqrtElementMetric

    REAL(KIND=dp) :: Force(4), InertialForce(3), NodalLame1(n),NodalLame2(n),Density, &
         Damping,Lame1,Lame2,NodalPressure(ntot),Pressure,NodalPressurePar(n),PressurePar
    REAL(KIND=dp) :: Grad(3,3),InvC(3,3),Identity(3,3),DetDefG
    REAL(KIND=dp) :: DefG(3,3), InvDefG(3,3),Strain(3,3), Stress2(3,3), Stress1(3,3)
    REAL(KIND=dp) :: dDefG(3,3),dStrain(3,3),dStress2(3,3),dStress1(3,3)
    REAL(KIND=dp) :: dDefGU(3,3),dStrainU(3,3),dStress2U(3,3),dStress1U(3,3)

    REAL(KIND=dp) :: Temperature, HeatExpansion(3,3)
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

      IF (PlaneStress) CALL Fatal( Caller,  &
          'Mixed formulation does not support plane stress')

      DOFs = cdim + 1
      ! To reuse the code, set the lambda parameter to zero and instead
      ! introduce the epsilon parameter = 1/lambda:
      NodalLame1(1:n) = 0.0d0
      IF ( ALL( ABS(NodalPoisson(1:n)) < AEPS ) ) THEN
        CALL Fatal( Caller,  &
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

    IntegStuff = GaussPoints( element, RelOrder = RelIntegOrder )

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
       ! We could have an entry for loss of volume
       DO i=1,dofs         
         Force(i) = SUM( LoadVector(i,1:n)*Basis(1:n) )
       END DO
       DO i=1,cdim
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
          ! (r, z, phi), hoop at 3 -- see the note in LocalMatrix.
          Grad(1,1) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(1,2) = SUM( LocalDisplacement(1,1:ntot) * dBasisdx(1:ntot,2) )
          Grad(3,3) = 1.0d0/r * SUM( LocalDisplacement(1,1:ntot) * Basis(1:ntot) )
          Grad(2,1) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,1) )
          Grad(2,2) = SUM( LocalDisplacement(2,1:ntot) * dBasisdx(1:ntot,2) )
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
                   Grad(1,2) = dBasisdx(p,2)
                   Grad(3,3) = 1.0d0/r * Basis(p)
                CASE (2)
                   Grad(2,1) = dBasisdx(p,1)
                   Grad(2,2) = dBasisdx(p,2)
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
                     -DDOTPROD(Grad,Stress1,dim) &
                     +DDOTPROD(Grad,dStress1U,dim))*s
                
                DO q = 1,ntot
                   DO j = 1,cdim
                      SELECT CASE(j)
                      CASE(1)
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(1,1) + dBasisdx(q,2)*dStress1(1,2) &
                              + 1.0d0/r*Basis(q)*dStress1(3,3))*s
                      CASE(2)
                         StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                              + (dBasisdx(q,1)*dStress1(2,1) + dBasisdx(q,2)*dStress1(2,2) ) * s
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
       !  Utilize the nodal damping:
       !-------------------------------------
       IF( GotDamping ) THEN
         DO p = 1,ntot
           DO q = 1,ntot
             DO i = 1,cdim
               DampMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   = DampMatrix(DOFs*(p-1)+i,DOFs*(q-1)+i) &
                   + Basis(p)*Basis(q)*Damping*s
             END DO
           END DO
         END DO
       END IF
       
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
                 ! (r, z, phi), and note InvDefG is indexed TRANSPOSED here --
                 ! the contraction is InvDefG(b,a) * dGrad(a,b), so the shear
                 ! partner of dGrad(1,2) is InvDefG(2,1) and not InvDefG(1,2).
                 CASE(1)
                   StiffMatrix(DOFs*p,DOFs*(q-1)+i) = StiffMatrix(DOFs*p,DOFs*(q-1)+i) + &
                       DetDefG**2 * ( dBasisdx(q,1) * InvDefG(1,1) + dBasisdx(q,2) * InvDefG(2,1) + &
                       Basis(q)/r * InvDefG(3,3) ) *  Basis(p) * s
                 CASE(2)
                   StiffMatrix(DOFs*p,DOFs*(q-1)+i) = StiffMatrix(DOFs*p,DOFs*(q-1)+i) + &
                       DetDefG**2 * ( dBasisdx(q,1) * InvDefG(1,2) + dBasisdx(q,2) * InvDefG(2,2) ) * &
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
                 Grad(1,2) = dBasisdx(p,2)
                 Grad(3,3) = 1.0d0/r * Basis(p)
               CASE (2)
                 Grad(2,1) = dBasisdx(p,1)
                 Grad(2,2) = dBasisdx(p,2)
               END SELECT

             ELSE
               Grad(i,:) = dBasisdx(p,:)
             END IF

             ForceVector(DOFs*(p-1)+i) = ForceVector(DOFs*(p-1)+i) &
                 + Pressure * dBasisdx(p,i) * s &
                 - Pressure * DDOTPROD(TRANSPOSE(InvDefG),Grad,dim) * s

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
                         + dBasisdx(q,2) * Grad(1,2) &
                         + Basis(q)/r * Grad(3,3) ) * s
                   CASE(2)
                     StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         = StiffMatrix(DOFs*(p-1)+i,DOFs*(q-1)+j) &
                         - Pressure * ( dBasisdx(q,1) * Grad(2,1) &
                         + dBasisdx(q,2) * Grad(2,2) ) * s
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
                     DDOTPROD(TRANSPOSE(InvDefG),Grad,dim) * s 
               END IF

             END DO
           END DO

           ! Source/drain for volume
           ForceVector(DOFs*p) = ForceVector(DOFs*p) &
               + Basis(p)*Force(dofs)*s  ! DetDefG - to multiply with this or not?                       
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
       LoadVector, NodalSpringCoeff,GotSpring,NormalSpring,NodalAlpha,NodalBeta, &
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
    LOGICAL :: GotSpring, NormalSpring
    REAL(KIND=dp) :: NodalAlpha(:,:), NodalBeta(:)
    REAL(KIND=dp) :: LocalDisplacement(:,:)
    TYPE(Element_t), TARGET :: Element
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
    REAL(KIND=dp) :: x(n),y(n),z(n), fx(n), fy(n), fz(n), Density
    REAL(KIND=dp) :: PBasis(pntot)
    REAL(KIND=dp) :: PdBasisdx(pntot,3),PSqrtElementMetric
    REAL(KIND=dp) :: FBasis(fn)
    REAL(KIND=dp) :: FdBasisdx(fn,3),FSqrtElementMetric
    REAL(KIND=dp) :: u,v,w,s,r,ParentU,ParentV,ParentW
    REAL(KIND=dp) :: FlowStress(3,3),Viscosity
    REAL(KIND=dp) :: Force(3),Alpha(3),Beta,Normal(3),RefNormal(3),FluidNormal(3),Identity(3,3)
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

    IP = GaussPoints( Element, RelOrder = RelIntegOrder )

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
          !
          ! In the case of an internal boundary (the model uses parent elements
          ! on both sides), the following command should create the normal vector
          ! pointing outwards from the structural body:
          !
          FluidNormal = NormalVector(Element,Nodes,u,v,Parent=Parent) 
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
         ! In the case of an internal boundary, the direction of RefNormal may not
         ! point outwards from the structural body, so we need to test to obtain
         ! the right sign:
         IF ( SUM( RefNormal * FluidNormal ) > 0 ) THEN
           Force = Force + MATMUL( FlowStress, Normal )
         ELSE
           Force = Force - MATMUL( FlowStress, Normal )
         END IF           
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

       IF( GotSpring ) THEN
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
    LOGICAL :: GlobalBubbles, Rebuilt
    INTEGER :: StrainDim
    REAL(KIND=dp), ALLOCATABLE :: SForceG(:)

    TYPE(NodalProjector_t), SAVE :: Proj
    PROCEDURE(ProjectedTensors_i) :: ElasticStrainAtIP

    SAVE SForceG, StrainDim
 !--------------------------------------------------------------------------------
    IF (.NOT. CalculateStrains) RETURN

    ! The saved matrix and permutation describe the mesh we last ran on. After a
    ! refinement they are stale, and gluing element contributions into them
    ! corrupts memory, so rebuild. As in ComputeStress the auxiliary solver object
    ! is kept and reassigned rather than freed: the variable added below records it
    ! as its owner, and VariableGet dereferences that owner
    ! (ListGetString(PVar % Solver % Values,'Equation')), so releasing it would
    ! leave the previous mesh's variable list pointing at freed memory.
    ! Resolved here rather than inside the projector, which sits below MainUtils.
    GlobalBubbles = SetGlobalBubblesFlag( Solver )

    CALL NodalProjectorSetup( Proj, Solver, 'strain:', 'Calculate Strains', &
        'StrainTemp', GlobalBubbles, Rebuilt )

    StSolver => Proj % PSolver

    IF ( Rebuilt ) THEN
       IF ( ALLOCATED(SForceG) ) DEALLOCATE( SForceG )
       StrainDim = SymTensorOutputComponents( CoordinateSystemDimension() )
       ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StrainDim) )
    END IF

    ! Limiters, contact conditions, residual mode, eigen/harmonic settings and the
    ! relaxation factor belong to the primary solve, not to an L2 fit; put aside
    ! until NodalProjectorEnd.
    CALL NodalProjectorBegin( Proj, Solver )
    NodalStrain = 0.0d0

    CALL NodalProjectorAssemble( Proj, StrainDim, AxialSymmetry, ElasticStrainAtIP, SForceG )


    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    CALL Info(Caller,'Calculating strain components',Level=7)

    CALL NodalProjectorSolve( Proj, 'Strain', StrainDim, SForceG, NodalStrain, Perm )

    CALL NodalProjectorEnd( Proj, Solver )

    CALL Info(Caller,'Finished strain postprocessing',Level=7)
!--------------------------------------------------------------------------------
  END SUBROUTINE GenerateStrainVariable
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
  SUBROUTINE GenerateStressVariable( NodalStress, Perm, &
       CalculateStress, AxialSymmetry)
!--------------------------------------------------------------------------------
!   This subroutine generates the stress field for material models which
!   depend on a list of state variables
!--------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: NodalStress(:)
    INTEGER, POINTER :: Perm(:) 
    LOGICAL :: CalculateStress, AxialSymmetry
 !---------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: StSolver
    LOGICAL :: GlobalBubbles, Rebuilt
    INTEGER :: StressDim
    REAL(KIND=dp), ALLOCATABLE :: SForceG(:)

    TYPE(NodalProjector_t), SAVE :: Proj
    PROCEDURE(ProjectedTensors_i) :: ElasticUmatStressAtIP

    SAVE SForceG, StressDim
 !--------------------------------------------------------------------------------------------
    IF (.NOT. CalculateStress) RETURN

    ! Rebuilt on mesh change -- see the note in GenerateStrainVariable.
    ! Resolved here rather than inside the projector, which sits below MainUtils.
    GlobalBubbles = SetGlobalBubblesFlag( Solver )

    CALL NodalProjectorSetup( Proj, Solver, 'stress:', 'Calculate Stresses', &
        'StressTemp', GlobalBubbles, Rebuilt )

    StSolver => Proj % PSolver

    IF ( Rebuilt ) THEN
       IF ( ALLOCATED(SForceG) ) DEALLOCATE( SForceG )

       ! One count now, where there used to be a StressDim of 4 in any 2D case
       ! written into a StressComponents-wide variable that was 6 unless
       ! axisymmetric -- so a plane case left two slots permanently unwritten.
       ! That mismatch was this routine's own TO DO and it is what is gone.
       StressDim = SymTensorOutputComponents( CoordinateSystemDimension() )

       ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StressDim) )
    END IF

    ! Limiters, contact conditions, residual mode, eigen/harmonic settings and the
    ! relaxation factor belong to the primary solve, not to an L2 fit; put aside
    ! until NodalProjectorEnd.
    CALL NodalProjectorBegin( Proj, Solver )
    NodalStress = 0.0d0

    CALL NodalProjectorAssemble( Proj, StressDim, AxialSymmetry, &
        ElasticUmatStressAtIP, SForceG )

    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    CALL Info(Caller,'Calculating stress components',Level=7)

    CALL NodalProjectorSolve( Proj, 'Stress', StressDim, SForceG, NodalStress, Perm )

    CALL NodalProjectorEnd( Proj, Solver )

    CALL Info(Caller,'Finished stress postprocessing',Level=7)
!----------------------------------------------------------------------------------
  END SUBROUTINE GenerateStressVariable
!----------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------
  SUBROUTINE ComputeStressAndStrain( Displacement, NodalStrain, NodalStress, VonMises, Perm, &
       PrincipalStress, PrincipalStrain, Tresca, PrincipalAngle, AxialSymmetry, &
       NeoHookeanMaterial, CalculateStrains, CalculateStresses, CalcPrincipal, &
       CalcPrincipalAngle, MixedFormulation, LargeDeflection, LinearIncompressible)
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: Displacement(:), NodalStrain(:), NodalStress(:), VonMises(:), &
         PrincipalStress(:), PrincipalStrain(:), Tresca(:), PrincipalAngle(:) 
    INTEGER, POINTER :: Perm(:)
    LOGICAL :: CalculateStrains, CalculateStresses, CalcPrincipal, CalcPrincipalAngle, &
         NeoHookeanMaterial, AxialSymmetry, MixedFormulation, LargeDeflection
    LOGICAL :: LinearIncompressible
!--------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: StSolver
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(ValueList_t), POINTER :: Equation

    INTEGER, POINTER :: Permutation(:), Indices(:)

    INTEGER :: dim, cdim, n, nd, elem, i, j, k, l, p, q, t, StrainDim, DOFs

    REAL(KIND=dp), POINTER :: StressTemp(:)
    REAL(KIND=dp), ALLOCATABLE :: ForceG(:), SForceG(:), LocalDisplacement(:,:), &
         Mass(:,:), Force(:), SForce(:), Basis(:), dBasisdx(:,:), &
         NodalLame1(:), NodalLame2(:)

    TYPE(MaterialPoint_t) :: MatPoint
    TYPE(MaterialResponse_t) :: MatResponse
    TYPE(MaterialState_t) :: MatState
    TYPE(MaterialModel_t) :: MatModel
    ! Sized by the largest layout any model selected here asks for, which is the
    ! anisotropic one's flattened 6x6.
    REAL(KIND=dp) :: MatProps(ANISOLIN_NPROPS)
    ! InvC, InvDefG and Pres went with the neo-Hookean block that moved into the
    ! constitutive interface. InvDefG was already dead before the move: it was
    ! assigned and inverted at every integration point and never read, so the
    ! postprocessing loop was paying for one matrix inversion per point for
    ! nothing.
    REAL(KIND=dp) :: Strain(3,3), Stress(3,3), Stress2(3,3), Grad(3,3), DefG(3,3), Identity(3,3), &
         u, v, w, Weight, detJ, res, Lame1, Lame2, nu, DetDefG, G(6,6), r
    ! The temperature difference at the integration point. Not called Temperature:
    ! that name belongs to the host's temperature FIELD, which this routine reads.
    REAL(KIND=dp) :: TempAtIp
    ! The pressure unknown of the linear mixed formulation at the integration point.
    REAL(KIND=dp) :: Pressure
    ! Material data read at the integration point rather than interpolated.
    LOGICAL :: EvalYoungIP, EvalPoissonIP, EvalBetaIP
    REAL(KIND=dp) :: YoungAtIP, PoissonAtIP
    LOGICAL :: NeedHeat
    ! The plane stress out-of-plane strain coefficients, filled by the condensation
    ! and meaningful only under plane stress. See CondensePlaneElasticityMatrix.
    REAL(KIND=dp) :: EzzC(3)

    LOGICAL :: FirstTime = .TRUE., Found, OptimizeBW, GlobalBubbles, Stat, &
         PlaneStress, &
         Isotropic, UseMask, LimiterOn, ContactOn, ResidualOn

    CHARACTER(LEN=MAX_NAME_LEN) :: eqname

    TYPE(NodalProjector_t), SAVE :: Proj
    LOGICAL :: Rebuilt

    SAVE ForceG, SForceG, Nodes, StrainDim
    !---------------------------------------------------------------------------------------------
    ! These variables are needed for Principal stress calculation;
    ! they are quite small and allocated even if principal stress calculation
    ! is not requested
    !------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: PriCache(3,3), PriTmp, PriW(3),PriWork(102)
    INTEGER       :: PriN=3, PriLWork=102, PriInfo=0
 !----------------------------------------------------------------------------------------------
    cdim = CoordinateSystemDimension()

    ! The dimensionality of the stress/strain state:
    IF (AxialSymmetry) THEN
       dim = 3
    ELSE
       dim = cdim
    END IF

    ! Either mixed formulation puts the pressure last at each node, so the stride
    ! is one more than the displacement components in both.
    IF (MixedFormulation .OR. LinearIncompressible) THEN
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

    ! Rebuilt on mesh change -- see the note in GenerateStrainVariable.
    ! Resolved here rather than inside the projector, which sits below MainUtils.
    GlobalBubbles = SetGlobalBubblesFlag( Solver )

    ! Unlike the other two callers this one registers the hidden variable against
    ! the field permutation rather than the projection's own, so it is passed
    ! explicitly. Whether that difference is deliberate is a separate question.
    CALL NodalProjectorSetup( Proj, Solver, 'stress:', 'Calculate Stresses', &
        'StressTemp', GlobalBubbles, Rebuilt, VarPerm = Perm )

    StSolver => Proj % PSolver
    Permutation => Proj % Perm
    UseMask = Proj % UseMask
    eqname = Proj % EqName

    IF ( Rebuilt ) THEN
       IF ( ALLOCATED(ForceG) ) DEALLOCATE( ForceG )
       IF ( ALLOCATED(SForceG) ) DEALLOCATE( SForceG )

       ! The mesh dimension, not the local "dim" above, which is the stress state's
       ! and is forced to 3 on the axis.
       StrainDim = SymTensorOutputComponents( cdim )

       IF (CalculateStrains) ALLOCATE( ForceG(StSolver % Matrix % NumberOfRows*StrainDim) )
       IF (CalculateStresses) ALLOCATE( SForceG(StSolver % Matrix % NumberOfRows*StrainDim) )
    END IF

    

    ! Limiters, contact conditions, residual mode, eigen/harmonic settings and the
    ! relaxation factor belong to the primary solve, not to an L2 fit; put aside
    ! until NodalProjectorEnd.
    CALL NodalProjectorBegin( Proj, Solver )
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
       ELSE IF ( GetLogical( Material, 'Youngs Modulus at IP', Found ) ) THEN
          ! At the integration points, so not readable at the nodes at all -- and
          ! isotropy assumed with it, as in the assembly and as in StressSolve.
          ElasticModulus = 0.0_dp
          Isotropic = .TRUE.
       ELSE
          CALL InputTensor( ElasticModulus, Isotropic, &
               'Youngs Modulus', Material, n, Indices )
       END IF
       EvalYoungIP = GetLogical( Material, 'Youngs Modulus at IP', Found )
       EvalPoissonIP = GetLogical( Material, 'Poisson Ratio at IP', Found )
       EvalBetaIP = GetLogical( Material, 'Heat Expansion Coefficient IP', Found )

       ! Selected per element rather than per integration point, since the choice
       ! turns on which keywords the material and solver gave and not on position.
       ! Neo-Hookean is tested first because it sets Isotropic itself.
       IF ( NeoHookeanMaterial ) THEN
          MatModel = NeoHookeanModel()
       ELSE IF ( Isotropic ) THEN
          MatModel = IsotropicLinearModel()
       ELSE
          MatModel = AnisotropicLinearModel()
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
                   CALL Fatal(Caller,'Given > Material Coordinate Unit Vector < too short!')
                END IF
                TransformMatrix(i,1:3) = Uwrk(1:3,1) / UnitNorm  
                RotateModuli = .TRUE.
             ELSE 
                TransformMatrix(i,1:3) = 0.0_dp
                TransformMatrix(i,i) = 1.0_dp
             END IF
          END DO

          IF( .NOT. RotateModuli  ) THEN
             CALL Fatal( Caller, &
                  'No unit vectors found but > Rotate Elasticity Tensor < set True?' )
          END IF
       END IF

       PoissonRatio = 0.0d0
       ! Only the isotropic branch below assigns this, and every reading of it is
       ! guarded by Isotropic -- but Fortran does not promise to stop evaluating an
       ! .AND. once it is decided, so an anisotropic material would have it read
       ! while undefined. Defined here instead, which changes no outcome.
       PlaneStress = .FALSE.
       IF (Isotropic) THEN
          ! As in the assembly, unread under "Incompressible" -- the branch at the
          ! end of this block overwrites both Lame parameters from Young alone.
          IF ( .NOT. EvalPoissonIP .AND. .NOT. LinearIncompressible ) &
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

          ! The linear mixed formulation, as in the assembly: the ratio is 1/2 by
          ! construction rather than read, so lambda goes and mu is E/3. Set after
          ! the branch above so that whatever the material said about the ratio is
          ! overridden here too, and the reported stress follows the same law the
          ! system was assembled from.
          IF ( LinearIncompressible ) THEN
            NodalLame1(1:n) = 0.0d0
            NodalLame2(1:n) = ElasticModulus(1,1,1:n) / 3.0d0
            PlaneStress = .FALSE.
          END IF
       ELSE IF ( dim == 2 ) THEN
          ! An anisotropic material in the plane needs the assumption too, since it
          ! decides which out-of-plane component the condensation leaves to be
          ! recovered. Read only when dim is 2, so that a three dimensional
          ! anisotropic run keeps the .FALSE. above and behaves exactly as before.
          PlaneStress = GetLogical( Equation, 'Plane Stress', Found )
       END IF


       ! The thermal state of this element, gathered exactly as the assembly loop
       ! gathers it: the temperature DIFFERENCE from the material's reference
       ! temperature, and the expansion coefficient in whatever shape the sif gave.
       ! What is reported below is then the ELASTIC strain, the thermal part taken
       ! out, and the stress that follows from it -- which is what StressSolve
       ! reports and the only reading consistent with the stress.
       NeedHeat = .FALSE.
       IF ( ASSOCIATED( TempSol ) .AND. .NOT. NeoHookeanMaterial ) THEN
          CALL InputTensor( HeatExpansionCoeff, IsotropicHeatExpansion, &
               'Heat Expansion Coefficient', Material, n, Indices, Found )
          IF ( Found ) THEN
             ReferenceTemperature(1:n) = ListGetReal( Material, 'Reference Temperature', &
                  n, Indices, Found )
             IF ( .NOT. Found ) ReferenceTemperature(1:n) = 0.0_dp
             LocalTemperature(1:n) = 0.0_dp
             WHERE( TempPerm( Element % NodeIndexes(1:n) ) > 0 )
                LocalTemperature(1:n) = Temperature(TempPerm(Element % NodeIndexes(1:n))) - &
                     ReferenceTemperature(1:n)
             END WHERE
             NeedHeat = ANY( LocalTemperature(1:n) /= 0.0_dp )
          END IF
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
           ! The projection is an L2 fit, so in axisymmetric coordinates it has to
           ! be weighted by the radius like any other volume integral. Without this
           ! the fit is made in the Cartesian metric and every spatially varying
           ! component comes out different from the correctly weighted one.
           IF (AxialSymmetry) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )

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

             ! In the plane the constitutive model wants the reduced three
             ! component packing, not the raw 6x6 -- see the failure mode spelled
             ! out in AnisotropicLinearStress. This is the single step that
             ! separated a solver with anisotropy in the plane from one without.
             IF ( dim == 2 ) &
                 CALL CondensePlaneElasticityMatrix( G, PlaneStress, EzzC )
          END IF

          Grad = 0.0d0
          Grad(1:cdim,1:cdim) = MATMUL( LocalDisplacement(1:cdim,1:nd), dBasisdx(1:nd,1:cdim) )
          IF (AxialSymmetry) THEN
             r = SUM(Basis(1:n) * Nodes % x(1:n))
             Grad(3,3) = 1.0d0/r * SUM(LocalDisplacement(1,1:nd) * Basis(1:nd))
          END IF
          IF (LargeDeflection) THEN
             DefG = Identity + Grad
          ELSE
             DefG = Identity
          END IF

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

          Strain = (TRANSPOSE(Grad)+Grad)/2.0D0
          IF (LargeDeflection) Strain = Strain + MATMUL(TRANSPOSE(Grad),Grad)/2.0D0

          ! The thermal part taken out, leaving the elastic strain -- and taken out
          ! HERE, before the out-of-plane recovery below, because the recovery is a
          ! statement about the elastic strain and StressSolve orders the two the
          ! same way. Over Dim, so the hoop is included under axial symmetry and the
          ! out-of-plane direction is not in the plane.
          IF ( NeedHeat ) THEN
             TempAtIp = SUM( LocalTemperature(1:n)*Basis(1:n) )
             DO i=1,dim
                Strain(i,i) = Strain(i,i) - &
                     SUM( HeatExpansionCoeff(i,i,1:n)*Basis(1:n) ) * TempAtIp
             END DO
          END IF

          ! Under plane stress the out-of-plane strain is determined by the in-plane
          ! ones and is not carried by the system, so it is recovered for output.
          ! The anisotropic form needs the coefficients the condensation set aside,
          ! and has a shear term the isotropic one does not.
          IF ( PlaneStress ) THEN
             IF ( Isotropic ) THEN
                Strain(3,3) = -nu/(1.0d0-nu)*(Strain(1,1)+Strain(2,2))
             ELSE
                Strain(3,3) = PlaneStressStrainZZ( EzzC, Strain )
             END IF
          END IF

          ! Every law goes through the constitutive interface now, so the choice
          ! here is only what to put in front of it. In the isotropic cases the
          ! identity the model builds for itself is the same modified one assembled
          ! above -- CDim-diagonal, out-of-plane entry only away from plane stress
          ! -- which is shared code in the interface rather than a coincidence.
          MatPoint % Strain = Strain
          MatPoint % Dim = dim
          MatPoint % CDim = cdim
          MatPoint % PlaneStress = PlaneStress
          MatPoint % AxiSymmetric = AxialSymmetry
          MatPoint % Kinematics = MERGE( KINEMATICS_LARGE_DEFLECTION, &
              KINEMATICS_SMALL_STRAIN, LargeDeflection )
          IF (NeoHookeanMaterial) THEN
             ! The deformation gradient, which only a large-deflection law reads,
             ! and the pressure if the system carries one as an unknown. That
             ! pressure is the whole reason this branch could not be a model until
             ! MaterialPoint_t had a place for solution data; the volumetric
             ! alternative to it is a constitutive relation and has moved into the
             ! model, which is why no dilatation formula is left here.
             MatPoint % F = DefG
             MatPoint % DetF = DetDefG
             MatPoint % PressureSupplied = MixedFormulation
             IF (MixedFormulation) &
                 MatPoint % Pressure = -SUM(LocalDisplacement(DOFs,1:n) * Basis(1:n))
             MatProps(NEOHOOKE_LAME1) = Lame1
             MatProps(NEOHOOKE_LAME2) = Lame2
          ELSE IF ( Isotropic ) THEN
             ! The same rebuilding as the assembly does when the data lives at the
             ! points: only the source of the two numbers differs, so the formulas are
             ! the ones above rather than a second convention.
             IF ( EvalYoungIP .OR. EvalPoissonIP ) THEN
                IF ( EvalYoungIP ) THEN
                   YoungAtIP = ListGetElementReal( YoungIP_h, Basis, GetCurrentElement(), &
                          Found, GaussPoint=t )
                ELSE
                   YoungAtIP = SUM( ElasticModulus(1,1,1:n)*Basis(1:n) )
                END IF
                IF ( EvalPoissonIP ) THEN
                   PoissonAtIP = ListGetElementReal( PoissonIP_h, Basis, GetCurrentElement(), &
                          Found, GaussPoint=t )
                ELSE
                   PoissonAtIP = SUM( PoissonRatio(1:n)*Basis(1:n) )
                END IF

                IF ( PlaneStress ) THEN
                   Lame1 = YoungAtIP * PoissonAtIP / ( 1.0d0 - PoissonAtIP**2 )
                ELSE
                   Lame1 = YoungAtIP * PoissonAtIP / &
                          ( (1.0d0 + PoissonAtIP) * (1.0d0 - 2.0d0*PoissonAtIP) )
                END IF
                Lame2 = YoungAtIP / ( 2.0d0 * (1.0d0 + PoissonAtIP) )
                nu = PoissonAtIP

                IF ( LinearIncompressible ) THEN
                   Lame1 = 0.0d0
                   Lame2 = YoungAtIP / 3.0d0
                END IF
             END IF
             MatProps(ISOLIN_LAME1) = Lame1
             MatProps(ISOLIN_LAME2) = Lame2
          ELSE
             ! The whole interpolated elasticity matrix, flattened -- the per-point
             ! material data question the interface was left holding.
             MatProps(ANISOLIN_C:ANISOLIN_C+ANISOLIN_NPROPS-1) = &
                 RESHAPE( G, [ ANISOLIN_NPROPS ] )
          END IF
          CALL MatModel % Stress( MatPoint, MatProps, MatState, MatResponse )
          Stress2 = MatResponse % Stress

          ! The linear mixed formulation's stress is sigma = 2 mu eps - p I, and the
          ! pressure part is added here rather than inside the model for the same
          ! reason the affine offset of "Stress Load" is: the assembly's tangent
          ! calls the law on strain increments through the same MaterialPoint_t, so a
          ! pressure carried in the model would enter the derivative as well. It is
          ! the identity of the state of STRESS, so it is CDim-diagonal here just as
          ! it is in the constraint the pressure came from.
          IF ( LinearIncompressible ) THEN
             Pressure = SUM( LocalDisplacement(DOFs,1:n) * Basis(1:n) )
             DO i = 1,cdim
                Stress2(i,i) = Stress2(i,i) - Pressure
             END DO
          END IF

          Stress =  1.0d0/DetDefG * MATMUL( MATMUL(DefG,Stress2), TRANSPOSE(DefG) )

          ! The plane strain counterpart of the recovery above: here it is the
          ! out-of-plane STRESS that the plane system determines without carrying.
          ! After the push-forward, since the modified identity leaves DefG(3,3)
          ! zero in the plane and the term would otherwise be annihilated. The
          ! isotropic case is not handled here because its Lame parameters already
          ! put the right value in Stress2(3,3).
          IF ( dim == 2 .AND. .NOT. Isotropic .AND. .NOT. PlaneStress ) &
              Stress(3,3) = Stress(3,3) + PlaneStrainStressZZ( G, Strain )

          CALL NodalProjectorMass( Mass, Basis, nd, Weight )
          IF (CalculateStrains) &
              CALL NodalProjectorTensor( Force, Basis, nd, StrainDim, Weight, Strain )
          IF (CalculateStresses) &
              CALL NodalProjectorTensor( SForce, Basis, nd, StrainDim, Weight, Stress )
       END DO

       CALL DefaultUpdateEquations( Mass, Force )

       !--------------------------------
       ! Assemble global RHS vectors:
       !--------------------------------   
       IF (CalculateStrains) &
           CALL NodalProjectorGlue( ForceG, Force, Permutation, Indices, nd, StrainDim )
       IF (CalculateStresses) &
           CALL NodalProjectorGlue( SForceG, SForce, Permutation, Indices, nd, StrainDim )

    END DO

    n = SIZE(StSolver % Variable % Values)
    !----------------------------------------------------------------------
    ! Linear solves componentwise...
    !-----------------------------------------------------------------------
    IF (CalculateStrains) THEN
       CALL Info(Caller,'Calculating strain components',Level=7)
       CALL NodalProjectorSolve( Proj, 'Strain', StrainDim, ForceG, NodalStrain, Perm )
    END IF

    IF (CalculateStresses) THEN
       CALL Info(Caller,'Calculating stress components',Level=7)
       CALL NodalProjectorSolve( Proj, 'Stress', StrainDim, SForceG, NodalStress, Perm )

       ! Von Mises stress from the component nodal values:
       ! -------------------------------------------------
       VonMises = 0
       IF (Identity(3,3) < 1.0d0) Identity(3,3) = 1.0d0
       Stress = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE

          q = StrainDim * (Perm(i)-1)
          CALL OutputVector2Tensor( NodalStress(q+1:q+StrainDim), StrainDim, Stress )

          Stress(:,:) = Stress(:,:) - TRACE(Stress(:,:),3) * Identity/3
          DO j=1,3
             DO k=1,3
                VonMises(Perm(i)) = VonMises(Perm(i)) + Stress(j,k)**2
             END DO
          END DO
       END DO
       VonMises = SQRT( 3.0d0 * VonMises / 2.0d0 )
    END IF


    !----------------------------------------------
    ! The principal and Tresca stresses:
    !--------------------------------------------------
    IF (CalcPrincipal .AND. CalculateStresses) THEN
       CALL Info(Caller,'Calculating principal stresses',Level=7)
       PriCache = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE
          q = StrainDim * (Perm(i)-1)
          CALL OutputVector2Tensor( NodalStress(q+1:q+StrainDim), StrainDim, PriCache )

          !-----------------------------------------------------------------------------
          ! Use lapack to solve the eigenvalues (i.e. the principal stresses)
          !-----------------------------------------------------------------------------
          CALL DSYEV( 'V', 'U', 3, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
          IF (PriInfo /= 0) THEN 
             CALL Fatal( Caller, 'DSYEV cannot generate eigen basis')
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
       CALL Info(Caller,'Calculating principal strains',Level=7)
       PriCache = 0.0d0
       DO i=1,SIZE( Perm )
          IF ( Perm(i) <= 0 ) CYCLE
          q = StrainDim * (Perm(i)-1)
          CALL OutputVector2Tensor( NodalStrain(q+1:q+StrainDim), StrainDim, PriCache )

          ! Use lapack to solve eigenvalues:
          CALL DSYEV( 'N', 'U', 3, PriCache, 3, PriW, PriWork, PriLWork, PriInfo )
          IF (PriInfo /= 0) THEN 
             CALL Fatal( Caller, 'DSYEV cannot generate eigen basis')
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

    CALL NodalProjectorEnd( Proj, Solver )


    CALL Info(Caller,'Finished postprocessing',Level=7)
!--------------------------------------------------------------------------------
  END SUBROUTINE ComputeStressAndStrain
!--------------------------------------------------------------------------------


!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ElasticSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The strain at one integration point, for GenerateStrainVariable's projection.
!>
!> File scope on purpose, and it has to be: an internal procedure passed as a
!> callback makes gfortran emit a stack trampoline, which marks the module as
!> needing an executable stack and stops it loading at all. See
!> ProjectedTensors_i. The consequence is that nothing here arrives by host
!> association -- the solver comes through Proj, and the rest is re-derived.
!------------------------------------------------------------------------------
SUBROUTINE ElasticStrainAtIP( Proj, Element, Nodes, n, nd, t, Basis, dBasisdx, T1, T2 )
!------------------------------------------------------------------------------
  USE StressLocal
  IMPLICIT NONE

  TYPE(NodalProjector_t) :: Proj
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n, nd, t
  REAL(KIND=dp) :: Basis(:), dBasisdx(:,:), T1(3,3), T2(3,3)
!------------------------------------------------------------------------------
  REAL(KIND=dp), ALLOCATABLE, SAVE :: LocalDisplacement(:,:)
  REAL(KIND=dp) :: Grad(3,3), r
  LOGICAL, SAVE :: LargeDeflection, AxialSymmetry
  LOGICAL :: Found
!------------------------------------------------------------------------------
  IF ( t == 1 ) THEN
    IF ( ALLOCATED(LocalDisplacement) ) THEN
      IF ( SIZE(LocalDisplacement,2) < Proj % Solver % Mesh % MaxElementDOFs ) &
          DEALLOCATE( LocalDisplacement )
    END IF
    IF ( .NOT. ALLOCATED(LocalDisplacement) ) &
        ALLOCATE( LocalDisplacement(3, Proj % Solver % Mesh % MaxElementDOFs) )

    ! Read rather than inherited, there being no host to inherit from. Both are
    ! cheap next to the element assembly, and reading them per element rather than
    ! once keeps this correct if either changes between calls.
    LargeDeflection = ListGetLogical( Proj % Solver % Values, 'Large Deflection', Found )
    IF ( .NOT. Found ) LargeDeflection = .TRUE.
    AxialSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric )

    CALL GetVectorLocalSolution( LocalDisplacement, USolver = Proj % Solver )
  END IF

  Grad = MATMUL( LocalDisplacement(:,1:nd), dBasisdx(1:nd,:) )
  IF ( AxialSymmetry ) THEN
    r = SUM( Basis(1:n) * Nodes % x(1:n) )
    Grad(3,3) = SUM( LocalDisplacement(1,1:nd) * Basis(1:nd) ) / r
  END IF

  T1 = ( TRANSPOSE(Grad) + Grad ) / 2.0_dp
  IF ( LargeDeflection ) T1 = T1 + MATMUL( TRANSPOSE(Grad), Grad ) / 2.0_dp
!------------------------------------------------------------------------------
END SUBROUTINE ElasticStrainAtIP
!------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
!> The UMAT's stress at one integration point, for GenerateStressVariable.
!>
!> File scope, like ElasticStrainAtIP and for the same reason. Nothing arrives by
!> host association, so the integration point variable and the component
!> permutation are found again here.
!>
!> The UMAT hands its stress back as a component vector in its own order rather
!> than as a tensor, so it is expanded into one for the assembly to pack again.
!> The round trip is exact: both directions are the same permutation, with zeros
!> for whatever the layout does not carry.
!--------------------------------------------------------------------------------
SUBROUTINE ElasticUmatStressAtIP( Proj, Element, Nodes, n, nd, t, Basis, dBasisdx, T1, T2 )
!--------------------------------------------------------------------------------
  USE StressLocal
  IMPLICIT NONE

  TYPE(NodalProjector_t) :: Proj
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  INTEGER :: n, nd, t
  REAL(KIND=dp) :: Basis(:), dBasisdx(:,:), T1(3,3), T2(3,3)
!--------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER, SAVE :: UmatStressVar => NULL()
  INTEGER, SAVE :: StressDofs, StressDim, Ind(6)
  INTEGER :: ipindex, i
  REAL(KIND=dp) :: Comp(6)
!--------------------------------------------------------------------------------
  IF ( t == 1 ) THEN
    UmatStressVar => VariableGet( Proj % Solver % Mesh % Variables, 'UmatStress' )
    StressDofs = UmatStressVar % Dofs
    StressDim = SymTensorOutputComponents( CoordinateSystemDimension() )

    ! Output slot -> the UmatStress slot it is read from. The UMAT keeps its own
    ! order, (rr,hoop,axial,rz) under axial symmetry and (11,22,33,12,13,23)
    ! otherwise, while the output is (11,22,33,12,23,13) with 33 out of plane. So
    ! axisymmetry swaps the hoop and axial slots, and 3D the last two shears.
    IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
         CurrentCoordinateSystem() == CylindricSymmetric ) THEN
      Ind = (/ 1, 3, 2, 4, 5, 6 /)
    ELSE
      Ind = (/ 1, 2, 3, 4, 6, 5 /)
    END IF
  END IF

  ipindex = GetIpIndex( t, USolver = Proj % Solver, Element = Element, &
      IpVar = UmatStressVar )

  Comp = 0.0_dp
  DO i=1,StressDim
    Comp(i) = UmatStressVar % Values( StressDofs*(ipindex-1) + Ind(i) )
  END DO

  CALL OutputVector2Tensor( Comp(1:StressDim), StressDim, T1 )
!--------------------------------------------------------------------------------
END SUBROUTINE ElasticUmatStressAtIP
!--------------------------------------------------------------------------------





!------------------------------------------------------------------------------
   SUBROUTINE ElasticSolver_Boundary_Residual( Model, Edge, Mesh, Quant, Perm, Gnorm, Indicator )
     USE StressLocal
     USE DefUtils
     IMPLICIT NONE
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     TYPE( Mesh_t )    :: Mesh
     TYPE( Element_t ) :: Edge
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
!------------------------------------------------------------------------------
     LOGICAL :: LargeDeflection, Found
!------------------------------------------------------------------------------
     ! These estimators used to be unconditionally geometrically nonlinear, which
     ! disagreed with the assembly from the moment that learned to honour the
     ! keyword: a case asking for small strain got a small-strain solution graded
     ! by a large-deflection error estimate. Same default as ElasticSolver's own
     ! reading, so a sif that does not mention it is unaffected.
     LargeDeflection = .TRUE.
     IF ( ASSOCIATED( Model % Solver ) ) THEN
       LargeDeflection = ListGetLogical( Model % Solver % Values, 'Large Deflection', Found )
       IF ( .NOT. Found ) LargeDeflection = .TRUE.
     END IF
     CALL ElasticityBoundaryResidual( Model, Edge, Mesh, Quant, Perm, Gnorm, Indicator, &
         LargeDeflection )
   END SUBROUTINE ElasticSolver_Boundary_Residual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ElasticSolver_Edge_Residual( Model,Edge,Mesh,Quant,Perm, Indicator )
     USE StressLocal
     USE DefUtils
     IMPLICIT NONE
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     TYPE(Mesh_t) :: Mesh
     TYPE(Element_t) :: Edge
     REAL(KIND=dp) :: Quant(:), Indicator(2)
!------------------------------------------------------------------------------
     LOGICAL :: LargeDeflection, Found
!------------------------------------------------------------------------------
     ! These estimators used to be unconditionally geometrically nonlinear, which
     ! disagreed with the assembly from the moment that learned to honour the
     ! keyword: a case asking for small strain got a small-strain solution graded
     ! by a large-deflection error estimate. Same default as ElasticSolver's own
     ! reading, so a sif that does not mention it is unaffected.
     LargeDeflection = .TRUE.
     IF ( ASSOCIATED( Model % Solver ) ) THEN
       LargeDeflection = ListGetLogical( Model % Solver % Values, 'Large Deflection', Found )
       IF ( .NOT. Found ) LargeDeflection = .TRUE.
     END IF
     CALL ElasticityEdgeResidual( Model, Edge, Mesh, Quant, Perm, Indicator, &
         LargeDeflection )
   END SUBROUTINE ElasticSolver_Edge_Residual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE ElasticSolver_Inside_Residual( Model, Element,  &
        Mesh, Quant, Perm, Fnorm, Indicator )
     USE StressLocal
     USE DefUtils
     IMPLICIT NONE
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     TYPE( Mesh_t )    :: Mesh
     TYPE( Element_t ) :: Element
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
!------------------------------------------------------------------------------
     LOGICAL :: LargeDeflection, Found
!------------------------------------------------------------------------------
     ! These estimators used to be unconditionally geometrically nonlinear, which
     ! disagreed with the assembly from the moment that learned to honour the
     ! keyword: a case asking for small strain got a small-strain solution graded
     ! by a large-deflection error estimate. Same default as ElasticSolver's own
     ! reading, so a sif that does not mention it is unaffected.
     LargeDeflection = .TRUE.
     IF ( ASSOCIATED( Model % Solver ) ) THEN
       LargeDeflection = ListGetLogical( Model % Solver % Values, 'Large Deflection', Found )
       IF ( .NOT. Found ) LargeDeflection = .TRUE.
     END IF
     CALL ElasticityInsideResidual( Model, Element, Mesh, Quant, Perm, Fnorm, Indicator, &
         LargeDeflection )
   END SUBROUTINE ElasticSolver_Inside_Residual
!------------------------------------------------------------------------------
