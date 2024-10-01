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

!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondFamily, LagrangeGauge, StaticConductivity, &
             ElectroDynamics
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(ValueListEntry_t), POINTER :: VariablePtr
  INTEGER, PARAMETER :: b_empty = 0, b_Piola = 1, &
      b_Secondorder = 2, b_Gauge = 4, b_Transient = 8, b_StaticCond = 16, &
      b_SecondFamily = 32
  INTEGER :: Paramlist
  CHARACTER(:), ALLOCATABLE :: ElemType

  TYPE(Solver_t), POINTER :: Solvers(:)
  INTEGER :: i,j,k,n
  CHARACTER(:), ALLOCATABLE :: eq
  INTEGER, POINTER :: ActiveSolvers(:)


  Paramlist = 0

  SolverParams => GetSolverParams()

  ElectroDynamics = ListGetLogical( SolverParams, 'Electrodynamics model', Found )
  IF ( ElectroDynamics) THEN
    CALL ListAddInteger( SolverParams, 'Time derivative order', 2)
  END IF

  StaticConductivity = ListGetLogical( SolverParams,'Static Conductivity',Found )
  IF( .NOT. Found ) THEN
    IF( ListCheckPrefixAnyBodyForce(Model, "Angular Velocity") .OR. &
        ListCheckPrefixAnyBodyForce(Model, "Lorentz Velocity") ) THEN
      CALL Info("WhitneyAVSolver_Init0", "Moving material triggers the use of scalar potential",Level=10)
      StaticConductivity = .TRUE.
    END IF

    IF( ListCheckPrefixAnyBC(Model, "Electric Current Density") ) THEN
      CALL Info("WhitneyAVSolver_Init0", &
          "> Electric Current Density < triggers the use of scalar potential",Level=10)    
      StaticConductivity = .TRUE.
    END IF
  END IF
  IF (.NOT. Transient .AND. StaticConductivity) THEN
    CALL Info("WhitneyAVSolver_Init0",'Including scalar potential in AV equation!',Level=6)
  END IF

  LagrangeGauge = GetLogical(SolverParams, 'Use Lagrange Gauge', Found)
  
  IF ( .NOT.ListCheckPresent(SolverParams, "Element") ) THEN
    ! We use one place where all the edge element keywords are defined and checked.
    CALL EdgeElementStyle(SolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE. )

    IF (PiolaVersion) Paramlist = Paramlist + b_Piola
    IF (SecondOrder) Paramlist = Paramlist + b_Secondorder
    IF (LagrangeGauge) Paramlist = Paramlist + b_Gauge
    IF (StaticConductivity) Paramlist = Paramlist + b_StaticCond
    IF (Transient .OR. ElectroDynamics) Paramlist = Paramlist + b_Transient
    IF (SecondFamily) ParamList = ParamList + b_SecondFamily
    
    SELECT CASE (Paramlist)
    CASE (b_Piola + b_Transient + b_Secondorder, &
         b_Piola + b_Gauge + b_Secondorder, &
         b_Piola + b_Transient + b_Secondorder + b_StaticCond, &
         b_Piola + b_Secondorder + b_StaticCond)
      ElemType = "n:1 e:2 -brick b:6 -prism b:2 -pyramid b:3 -quad_face b:4 -tri_face b:2"

      ! Note: only tets available for 2nd kind basis!
    CASE (b_Piola + b_Transient + b_SecondFamily, &
         b_Piola + b_Gauge + b_SecondFamily, &
         b_Piola + b_Transient + b_SecondFamily + b_StaticCond, &
         b_Piola + b_SecondFamily + b_StaticCond)
      ElemType = "n:1 e:2"

    CASE (b_Piola + b_Transient, &
         b_Piola + b_Transient + b_StaticCond, &
         b_Piola + b_Transient + b_Gauge)
      ElemType = "n:1 e:1 -brick b:3 -quad_face b:2" 

    CASE (b_Piola + b_Gauge)
      ElemType = "n:1 e:1 -brick b:3 -quad_face b:2" 

    CASE (b_Piola + b_Secondorder)
      ElemType = "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2" 

    CASE (b_Piola + b_Secondfamily)
      ElemType = "n:0 e:2"

    CASE (b_Piola)
      ElemType = "n:0 e:1 -brick b:3 -quad_face b:2"

    CASE (b_Piola + b_StaticCond )
      ElemType = "n:1 e:1 -brick b:3 -quad_face b:2" 

    CASE (b_Transient, &
         b_Transient + b_StaticCond, &
         b_StaticCond, &
         b_Gauge + b_Transient, &
         b_Gauge)
      ElemType = "n:1 e:1" 

    CASE (b_empty)
      ElemType = "n:0 e:1" 

    CASE default
      WRITE (Message,*) 'Unsupported degree-gauge-transient combination', Paramlist
      CALL Fatal('WhitneyAVSolver_Init0', Message)

    END SELECT

    CALL Info('WhitneyAVSolver_Init0','Setting element type to: "'//TRIM(ElemType)//'"',Level=6)
    CALL ListAddString( SolverParams,'Element',ElemType ) 

    
    IF( GetString(SolverParams,'Linear System Solver',Found) == 'block' ) THEN
      IF ( PiolaVersion ) THEN
        CALL Fatal('WhitneyAVSolver_Init0','Block strategy not applicable to piola version!')
      ELSE
        CALL ListAddLogical( SolverParams, "Optimize Bandwidth", .FALSE.)
      END IF
    END IF
  END IF

  IF (.NOT. Transient .AND. .NOT. ( StaticConductivity .OR. LagrangeGauge ) ) THEN
    CALL ListAddNewLogical( SolverParams,'Variable Output',.FALSE.)
  END IF
    
  CALL ListAddLogical( SolverParams,'Use Global Mass Matrix',.TRUE.) 

  ! This is for internal communication with the saving routines
  CALL ListAddLogical( SolverParams,'Hcurl Basis',.TRUE.)

  CALL ListAddNewString( SolverParams,'Variable','AV')
 
  IF (LagrangeGauge .AND. Transient .AND. &
      ListCheckPrefixAnyBC( Model, "Mortar BC" ) ) THEN
    CALL Info("WhitneyAVSolver_Init0", "Gauge field is not projected across mortar boundaries.") 
  END IF  
  
  ! THIS ENFORCES THE NEW STRATEGY !!!!
  CALL ListAddLogical( SolverParams,'Generic Source Fixing',.TRUE.)

  ! Add solvers, if "Helmholtz projection" requested
  IF (ListGetLogical(SolverParams, 'Helmholtz Projection', Found)) THEN
    Solvers => Model % Solvers
    n = Model % NumberOfSolvers
    Model % NumberOfSolvers = n+2

    ALLOCATE(Model % Solvers(Model % NumberOfSolvers))
    Model % Solvers(1:n) = Solvers

    DO i=n+1,n+2
      Model % Solvers(i) % PROCEDURE = 0
      Model % Solvers(i) % Matrix => Null()
      Model % Solvers(i) % Mesh => Null()
      Model % Solvers(i) % Variable => Null()
      Model % Solvers(i) % ActiveElements => Null()
      Model % Solvers(i) % NumberOfActiveElements = 0
    END DO

    DO i=1,Model % NumberOfSolvers
      IF(.NOT.ASSOCIATED(Model % Solvers(i) % Values)) &
         Model % Solvers(i) % Values => ListAllocate()
    END DO

    Eq =  ListGetString( SolverParams, 'Equation' )

    CALL ListAddIntegerArray( SolverParams, 'Post Solvers', 2, [n+1,n+2] )

    CALL ListAddString( Model % Solvers(n+1) % Values, 'Procedure', &
                'MagnetoDynamics HelmholtzProjectorT', CaseConversion=.FALSE. )
    CALL ListAddString( Model % Solvers(n+1) % Values, 'Equation', 'HP' )
    CALL ListAddString( Model % Solvers(n+1) % Values, 'Exec Solver', 'Never' )

    CALL ListAddString( Model % Solvers(n+2) % Values, 'Procedure', &
              'MagnetoDynamics RemoveKernelComponentT',CaseConversion=.FALSE. )
    CALL ListAddString( Model % Solvers(n+2) % Values, 'Equation', 'RMC' )
    CALL ListAddString( Model % Solvers(n+2) % Values, 'Exec Solver', 'Never' )

    DO i=1,Model % NumberOFEquations
      IF ( ListGetLogical( Model % Equations(i) % Values, Eq, Found ) ) THEN
        CALL ListAddLogical( Model % Equations(i) % Values, 'HP', .TRUE. )
        CALL ListAddLogical( Model % Equations(i) % Values, 'RMC' , .TRUE.)
      ELSE
        ActiveSolvers => ListGetIntegerArray( CurrentModel % Equations(i) % Values, &
                              'Active Solvers', Found )
        IF ( Found ) THEN
          DO k=1,SIZE(ActiveSolvers)
            IF ( ActiveSolvers(k) == Solver % SolverId ) THEN
              CALL ListAddLogical( Model % Equations(i) % Values, 'HP', .TRUE. )
              CALL ListAddLogical( Model % Equations(i) % Values, 'RMC', .TRUE. )
              EXIT
            END IF
          END DO
        END IF
      END IF
    END DO
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver_Init(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Found
  
  Mesh => GetMesh()
  IF( Mesh % MeshDim /= 3 ) THEN
    CALL Fatal('WhitneyAVSolver_Init','Solver requires 3D mesh!')
  END IF
  
  IF( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric ) THEN
    CALL Fatal('WhitneyAVSolver_Init','Solver not applicable to axially axisymmetric cases!')
  END IF

  ! Historically a real array could be used for H-B Curve.
  ! This dirty piece of code makes things backward compatible.
  BLOCK
    INTEGER :: i
    LOGICAL :: Cubic
    TYPE(ValueList_t), POINTER :: Material
    DO i=1,Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      IF( ListCheckPresent( Material, 'H-B Curve') ) THEN
        Cubic = GetLogical( Material, 'Cubic spline for H-B curve',Found)
        CALL ListRealArrayToDepReal(Material,'H-B Curve','dummy',&
            CubicTable=Cubic) !Monotone=.TRUE.)         
      END IF
    END DO
  END BLOCK
!------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVSolver_Init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>  Solve a vector potential A and scalar potential V from
! 
!>  sigma @A/@t + rot (1/mu) rot A + sigma grad(V) = J^s + curl(M^s) - sigma grad(V^s)
!>   -div(sigma*(@A/@t+grad(V))) = 0 .
!
!>  by using edge elements for A + nodal basis for V.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE MagnetoDynamicsUtils
  USE CircuitUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element, Edge

  REAL(KIND=dp) :: Norm, PrevDT=-1, RelChange
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, &
      BC, BodyParams, PrevMaterial

  INTEGER :: n,nb,nd,t,istat,i,j,k,l,nNodes,Active,FluxCount=0, &
          NoIterationsMax,NoIterationsMin, nsize

  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), POINTER :: VecPot(:)
  REAL(KIND=dp), POINTER :: Cwrk(:,:,:), Acoef_t(:,:,:) => NULL()
  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), Tcoef(:,:,:), &
                                GapLength(:), AirGapMu(:), LamThick(:), &
                                LamCond(:), Wbase(:), RotM(:,:,:), &
                                ThinLineCrossect(:),ThinLineCond(:)

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), DAMP(:,:), FORCE(:), &
            JFixFORCE(:), JFixVec(:,:),PrevSol(:)

  CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType
  LOGICAL :: LaminateStack, CoilBody, HasHBCurve, HasReluctivityFunction, &
      NewMaterial

  INTEGER, POINTER :: Perm(:)
  INTEGER, ALLOCATABLE :: FluxMap(:)
  LOGICAL, ALLOCATABLE, SAVE :: TreeEdges(:)
  LOGICAL :: Stat, TG, DoneAssembly=.FALSE., &
         SkipAssembly, ConstantSystem, ConstantBulk, JFix, JFixSolve, FoundRelax, &
         PiolaVersion, SecondOrder, LFact, LFactFound, EdgeBasis, &
         HasTensorReluctivity
  LOGICAL :: SteadyGauge, TransientGauge, TransientGaugeCollected=.FALSE., &
       HasStabC, RegularizeWithMass

  REAL(KIND=dp) :: Relax, gauge_penalize_c, gauge_penalize_m, gauge_coeff, &
      mass_reg_epsilon, newton_eps

  REAL(KIND=dp) :: NewtonTol
  INTEGER :: NewtonIter
  LOGICAL :: Newton
  
  TYPE(Variable_t), POINTER :: Var, JFixVar, CoordVar
  TYPE(Matrix_t), POINTER :: A
  TYPE(ListMatrix_t), POINTER, SAVE :: BasicCycles(:)
  TYPE(ValueList_t), POINTER :: CompParams
  TYPE(Matrix_t), POINTER :: CM=>NULL()
  
  INTEGER :: n_n, n_e
  INTEGER, POINTER :: Vperm(:), Aperm(:)
  REAL(KIND=dp), POINTER :: Avals(:), Vvals(:)

  CHARACTER(LEN=MAX_NAME_LEN):: CoilCurrentName
  TYPE(Variable_t), POINTER :: CoilCurrentVar
  REAL(KIND=dp) :: CurrAmp
  LOGICAL :: UseCoilCurrent, ElemCurrent, ElectroDynamics, EigenSystem

  TYPE(ValueHandle_t), SAVE :: mu_h 
  TYPE(Solver_t), POINTER :: pSolver

  REAL(KIND=dp) :: TorqueTol
  LOGICAL :: UseTorqueTol

  CHARACTER(*), PARAMETER :: Caller = 'WhitneyAVSolver'
  
  SAVE STIFF, LOAD, MASS, DAMP, FORCE, JFixFORCE, JFixVec, Tcoef, GapLength, AirGapMu, &
       Acoef, Cwrk, LamThick, LamCond, Wbase, RotM, AllocationsDone, &
       Acoef_t, ThinLineCrossect, ThinLineCond
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN	

  CALL Info(Caller,'',Level=6 )
  CALL Info(Caller,'-------------------------------------------------',Level=6 )
  IF( Transient ) THEN
    CALL Info(Caller,'Solving transient AV equations with edge elements',Level=5 )
  ELSE
    CALL Info(Caller,'Solving steady-state AV equations with edge elements',Level=5 )
  END IF
    
  SolverParams => GetSolverParams()
  pSolver => Solver

  ElectroDynamics = ListGetLogical( SolverParams, 'Electrodynamics model', Found )
  EigenSystem = ListGetLogical( SolverParams, 'Eigen Analysis', Found )
  
  CALL EdgeElementStyle(SolverParams, PiolaVersion, SecondOrder )
  IF (PiolaVersion) THEN
    CALL Info(Caller, &
        'Using Piola Transformed element basis functions',Level=4)
  END IF

  SteadyGauge = GetLogical(SolverParams, 'Use Lagrange Gauge', Found) .AND. .NOT. Transient
  TransientGauge = GetLogical(SolverParams, 'Use Lagrange Gauge', Found) .AND. Transient

  CoilCurrentName = GetString( SolverParams,'Current Density Name',UseCoilCurrent ) 
  IF(.NOT. UseCoilCurrent ) THEN
    UseCoilCurrent = GetLogical(SolverParams,'Use Nodal CoilCurrent',Found )
    IF(UseCoilCurrent) THEN
      CoilCurrentName = 'CoilCurrent'
    ELSE
      UseCoilCurrent = GetLogical(SolverParams,'Use Elemental CoilCurrent',Found )
      IF(UseCoilCurrent) CoilCurrentName = 'CoilCurrent e'
    END IF
  END IF
  ElemCurrent = .FALSE.
  IF( UseCoilCurrent ) THEN
    CoilCurrentVar => VariableGet(Solver % Mesh % Variables, CoilCurrentName )
    IF( ASSOCIATED( CoilCurrentVar ) ) THEN
      CALL Info(Caller,'Using precomputed field for current density: '//TRIM(CoilCurrentName),Level=5)
      IF( CoilCurrentVar % TYPE == Variable_on_nodes_on_elements ) THEN
        ElemCurrent = .TRUE.
      ELSE
        CALL Warn(Caller,'Precomputed CoilCurrent is not an elemental field!')
      END IF
    ELSE
      CALL Fatal(Caller,'Elemental current requested but not found:'//TRIM(CoilCurrentName))
    END IF
  END IF
  
  IF (SteadyGauge) THEN
    CALL Info("WhitneyAVSolver", "Utilizing Lagrange multipliers for gauge condition in steady state computation")
    IF(.not. ListCheckPresent( SolverParams, 'Linear System Refactorize') ) THEN
      CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
    END IF
  END IF

  IF(TransientGauge) THEN
    CALL Info("WhitneyAVSolver", "Utilizing Lagrange multipliers for gauge condition in transient computation")
    IF (.NOT. ListCheckPresent( SolverParams, "enforce exact dirichlet bcs" ) ) THEN
      CALL ListAddLogical(SolverParams,"enforce exact dirichlet bcs",.FALSE.)
      CALL Info("WhitneyAVSolver", "Setting 'enforce exact dirichlet bcs = Logical False'")
    END IF

    IF (.NOT. ListCheckPresent( SolverParams, "optimize bandwidth" ) ) THEN
      CALL ListAddLogical(SolverParams,"optimize bandwidth",.FALSE.)
      CALL Info("WhitneyAVSolver", "Setting 'Optimize Bandwidth = Logical False'") 

    ELSEIF (ListGetLogical(SolverParams, "Optimize bandwidth")) THEN
      CALL Warn("WhitneyAVSolver", &
           "Optimize bandwidth and use lagrange gauge in transient is known not to work. ")
    END IF

    IF(.not. ListCheckPresent( SolverParams, 'Linear System Refactorize') ) THEN
      CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
    END IF
    ! TODO: Check if there is mortar boundaries and report the above in that case only.
  END IF

  Newton = .FALSE.
  newton_eps = GetCReal(SolverParams, 'Newton epsilon', Found ) 
  IF(.NOT. Found) newton_eps = 1.0e-3
  
  mass_reg_epsilon = GetCReal(SolverParams, 'Mass regularize epsilon', RegularizeWithMass)
  IF (RegularizeWithMass .AND. mass_reg_epsilon == 0.0_dp) THEN
    RegularizeWithMass = .FALSE.
  END IF
  IF (RegularizeWithMass) THEN
    WRITE (Message, *) 'Mass regularization epsilon', mass_reg_epsilon
    CALL Info("WhitneyAVSolver", Message)
  END IF

  gauge_coeff = GetCReal(SolverParams, 'Lagrange Gauge coefficient',Found )
  IF(.NOT. Found ) gauge_coeff = 1.0_dp
  gauge_penalize_c = GetCReal(SolverParams, 'Lagrange Gauge Penalization coefficient', HasStabC)
  gauge_penalize_m = GetCReal(SolverParams, 'Lagrange Gauge Penalization coefficient mass', Found)
  HasStabC = HasStabC .OR. Found

  IF ( HasStabC ) THEN
    WRITE (Message, *) 'Lagrange Gauge penalization coefficient', gauge_penalize_c
    CALL Info(Caller, message)
    WRITE (Message, *) 'Lagrange Gauge penalization coefficient mass', gauge_penalize_m
    CALL Info(Caller, message)
  END IF

  ! Gauge tree, if requested or using direct solver:
  ! ------------------------------------------------
  TG = GetLogical(SolverParams,'Use tree gauge', Found)
  IF (.NOT. Found) THEN
    IF (.NOT. (SteadyGauge .OR. TransientGauge)) THEN
      IF( GetString(SolverParams,'Linear System Solver',Found)=='direct') THEN
        IF( PiolaVersion ) THEN
          CALL Fatal(Caller,'Direct solver (with tree gauge) is only possible with the lowest order edge basis!')
        ELSE
          CALL Info(Caller,'Defaulting to tree gauge when using direct solver')
          TG = .TRUE.
        END IF
      END IF
    END IF
  END IF

  IF( PiolaVersion .AND. TG ) THEN
    CALL Fatal(Caller,'Tree Gauge cannot be used in conjunction with Piola transformation')
  END IF

  TorqueTol = GetCReal(SolverParams,'Nonlinear System Torque Tolerance',UseTorqueTol)
  IF(UseTorqueTol) CALL Info(Caller,'Using additional nonlinear tolerance for torque',Level=10)

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()
  nNodes = Mesh % NumberOfNodes
  Perm => Solver % Variable % Perm
  Vecpot => Solver % Variable % Values

  IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough

     IF(ALLOCATED(FORCE)) THEN
       DEALLOCATE(FORCE, JFixFORCE, JFixVec, LOAD, STIFF, MASS, DAMP, TCoef, GapLength, AirGapMu, &
             Acoef, LamThick, LamCond, WBase, RotM, ThinLineCrossect, ThinLineCond )
     END IF

     ALLOCATE( FORCE(N), JFixFORCE(n), JFixVec(3,n), LOAD(7,N), STIFF(N,N), &
          MASS(N,N), DAMP(N,N), Tcoef(3,3,N), GapLength(N), &
          AirGapMu(N), Acoef(N), LamThick(N), &
          LamCond(N), Wbase(N), RotM(3,3,N),  &
          ThinLineCrossect(N), ThinLineCond(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( Caller, 'Memory allocation error.' )
     END IF

     IF(GetString(SolverParams,'Linear System Solver',Found)=='block') THEN
       n = Mesh % NumberOfNodes
       n_n = COUNT(Perm(1:n)>0)
       n_e = COUNT(Perm(n+1:)>0)
       ALLOCATE( Avals(n_e), Vvals(n_n) )
       Vvals = Vecpot(1:n)
       Avals = Vecpot(n+1:)

       ALLOCATE(Aperm(SIZE(Perm)),Vperm(SIZE(Perm)))
       Aperm = 0; Vperm = 0
       DO i=1,SIZE(AVals)
         Aperm(i+n) = i
       END DO
       DO i=1,SIZE(Vvals)
         Vperm(i)=i
       END DO

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, & 
                 GetVarName(Solver % Variable)//' 1',1,Vvals,Vperm)

       CALL VariableAdd(Mesh % Variables,Mesh,Solver, & 
                 GetVarName(Solver % Variable)//' 2',1,Avals,Aperm)
     END IF

     NULLIFY( Cwrk )

     IF (TransientGauge) THEN
       A => GetMatrix()
       IF (.NOT. TransientGaugeCollected) &
           CM => AddConstraintFromBulk(A, Solver % Matrix % ConstraintMatrix)
     END IF

     AllocationsDone = .TRUE.
  END IF

  ConstantSystem = GetLogical( SolverParams, 'Constant System', Found )
  ConstantBulk   = GetLogical( SolverParams, 'Constant Bulk System', Found )

  SkipAssembly = DoneAssembly.AND.(ConstantBulk.OR.ConstantSystem)

  JFix = GetLogical(SolverParams,'Fix input Current Density', Found)
  IF (.NOT. ( Found .OR. Transient ) ) THEN
    ! Only fix the current density if there is one
    JFix = ListCheckPrefixAnyBodyForce(Model, 'Current Density') 
  END IF
  JFixSolve = JFix

  IF (JFix) THEN
    JFixPhase = 1
    CALL JFixPotentialSolver(Model,Solver,dt,Transient)
    JFixVar => VariableGet(Mesh % Variables, 'JFix')    
    IF(.NOT. ASSOCIATED( JFixRhs ) ) THEN
      CALL Fatal(Caller,'JFixRhs should be associated!')
    END IF
    IF(.NOT. ASSOCIATED( JFixSurfacePerm ) ) THEN
      CALL Fatal(Caller,'JFixSurfacePerm should be associated!')
    END IF
    IF(.NOT. ALLOCATED( JFixSurfaceVec ) ) THEN
      CALL Fatal(Caller,'JFixSurfaceVec should be associated!')
    END IF   
  END IF

  ! 
  ! Use vec.pot. dofs only for convergence:
  ! ----------------------------------------
  CALL ListAddInteger(SolverParams,'Norm Permutation',nNodes+1)


  ! Resolve internal non.linearities, if requested:
  ! ----------------------------------------------
  NoIterationsMax = GetInteger( SolverParams, 'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found) NoIterationsMax = 1

  NoIterationsMin = GetInteger( SolverParams, 'Nonlinear System Min Iterations',Found)
  IF(.NOT. Found) NoIterationsMin = 1

  ! Use also these keyword for compatibility with ElmerGUI and old practices
  NewtonIter = GetInteger( SolverParams,&
      'Nonlinear System Newton After Iterations',Found ) 
  IF(.NOT. Found ) NewtonIter = NoIterationsMax

  NewtonTol = GetCReal( SolverParams,&
      'Nonlinear System Newton After Tolerance',Found )


! Not refactorizing seems to break things with gauges
  ! IF (SteadyGauge) THEN
  !   IF(.not. ListCheckPresent( SolverParams, 'Linear System Refactorize') ) THEN
  !     CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
  !   END IF
  ! END IF

  LFact = GetLogical( SolverParams,'Linear System Refactorize', LFactFound )
  IF ( dt /= PrevDT .AND. LFactFound .AND. .NOT. LFact ) THEN
    CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
  END IF
  EdgeBasis = .NOT.LFactFound .AND. GetLogical( SolverParams, 'Edge Basis', Found )

  CALL DefaultStart()
  
  DO i=1,NoIterationsMax
    Newton = GetLogical( SolverParams,'Newton-Raphson iteration',Found)
    IF(.NOT. Found ) Newton = ( i > NewtonIter .OR. Solver % Variable % NonlinChange < NewtonTol )
    IF( NoIterationsMax > 1 ) THEN
      CALL Info(Caller,'Nonlinear iteration: '//I2S(i),Level=8 )
    END IF

    IF( DoSolve(i) ) THEN
      IF(i>=NoIterationsMin) EXIT
    END IF
    IF( EdgeBasis ) CALL ListAddLogical(SolverParams,'Linear System Refactorize',.FALSE.)
    
    ! Currently assume that the source terms are constant over the nonlinear iteration
    JFixSolve = .FALSE.
  END DO
  IF ( EdgeBasis ) CALL ListRemove( SolverParams, 'Linear System Refactorize' )

  IF ( dt /= PrevDT .AND. LFactFound .AND. .NOT. LFact ) THEN
    CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .FALSE. )
  END IF
  PrevDT = dt

  IF(.NOT. UseTorqueTol) CALL CalculateLumpedParameters()
  
  CoordVar => VariableGet(Mesh % Variables,'Coordinates')
  IF(ASSOCIATED(CoordVar)) THEN
    DO i=1,Mesh % NumberOfNodes
      j = 3*(CoordVar % Perm(i)-1)
      CoordVar % Values(j+1) = Mesh % Nodes % x(i)
      CoordVar % Values(j+2) = Mesh % Nodes % y(i)
      CoordVar % Values(j+3) = Mesh % Nodes % z(i)
    END DO
  END IF

  CALL DefaultFinish()

  CALL Info(Caller,'All done',Level=8 )
  CALL Info(Caller,'-------------------------------------------',Level=8 )


CONTAINS

!------------------------------------------------------------------------------
  LOGICAL FUNCTION DoSolve(IterNo) RESULT(Converged)
!------------------------------------------------------------------------------
   IMPLICIT NONE
   CHARACTER(LEN=MAX_NAME_LEN) :: potname
   INTEGER :: i,j,k,t,n,nd,nb,IterNo

   REAL(KIND=dp)::TOL,Norm,PrevNorm, NonLinError, LinTol, RelTol, BaseTol
   LOGICAL :: Found, FoundMagnetization, CalculateNonlinearResidual, LFactFound
   LOGICAL :: AdaptiveTols, FoundAny, ConstraintActive, GotCoil

   TYPE(Matrix_t), POINTER :: MMatrix
   REAL(KIND=dp), POINTER :: Mx(:), Mb(:), Mr(:)
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpRVec, TmpRHSVec   
   CHARACTER(LEN=MAX_NAME_LEN) :: ConvergenceType
   REAL(KIND=dp),  POINTER CONTIG :: SaveValues(:), SaveRHS(:), ConstraintValues(:)
   REAL(KIND=dp) :: PrevTorque, TorqueErr, Torque
   
   
   SAVE TmpRHSVec, TmpRVec
  !-----------------
  !System assembly:
  !-----------------

  A => GetMatrix()
  IF (TransientGauge) ConstraintValues => CM % Values

300 CONTINUE

  IF ( SkipAssembly) THEN
    DO i=1,SIZE(A % RHS)
      A % RHS(i) = A % BulkRHS(i)
    END  DO      
    DO i=1,SIZE(A % Values)
      A % Values(i) = A % BulkValues(i)
    END DO
    IF ( ConstantBulk    ) GOTO 100
    IF ( ConstantSystem  ) GOTO 200
  END IF

  ! Timing
  CALL ResetTimer('MGDynAssembly')
  CALL DefaultInitialize()
  Active = GetNOFActive()
  
  IF( ListCheckPresentAnyMaterial(Model,'Reluctivity Function') ) THEN
    CALL ListInitElementKeyword( mu_h,'Material','Reluctivity Function',&
        EvaluateAtIp=.TRUE.,DummyCount=3)
  END IF
  
  PrevMaterial => NULL()
  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes() ! kulmat
     nd = GetElementNOFDOFs()  ! vapausasteet
     nb = GetElementNOFBDOFs()  ! sisäiset vapausasteet
             
     IF (SIZE(Tcoef,3) /= n) THEN
       DEALLOCATE(Tcoef)
       ALLOCATE(Tcoef(3,3,n), STAT=istat)
       IF ( istat /= 0 ) THEN
         CALL Fatal( Caller, 'Memory allocation error.' )
       END IF
     END IF
     
     LOAD = 0.0d0

     ! This way we don't have to inquire the list for all three components separately.
     ! Also writing of the sif file becomes more economical.
     
     BodyForce => GetBodyForce()
     FoundMagnetization = .FALSE.

     ! If the coil current field is elemental it is discontinuous and need not be limited
     ! to the body force. For nodal ones we don't have the same luxury.
     GotCoil = .FALSE.
     IF( UseCoilCurrent ) THEN
       IF( ElemCurrent .OR. ASSOCIATED(BodyForce) ) THEN
         CALL GetVectorLocalSolution( Load,UElement=Element,UVariable=CoilCurrentVar,Found=GotCoil)       
       END IF
     END IF
       

     IF ( ASSOCIATED(BodyForce) ) THEN       
       ! If not already given by CoilCurrent, request for current density
       IF( .NOT. GotCoil) THEN
         CALL GetRealVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
       END IF

       CurrAmp = ListGetCReal( BodyForce,'Current Density Multiplier',Found ) 
       IF(Found) Load(1:3,1:n) = CurrAmp * Load(1:3,1:n)
       
       CALL GetRealVector( BodyForce, Load(4:6,1:n), &
                'Magnetization', FoundMagnetization )
       Load(7,1:n) = GetReal( BodyForce, 'Electric Potential', Found )
     END IF

     Material => GetMaterial( Element )
     NewMaterial = .NOT. ASSOCIATED(Material, PrevMaterial)
     IF (NewMaterial) THEN              
       HasHBCurve = ListCheckPresent(Material, 'H-B Curve')
       HasReluctivityFunction = ListCheckPresent(Material,'Reluctivity Function')
       PrevMaterial => Material
     END IF
     
     IF(ASSOCIATED(Material).AND..NOT.FoundMagnetization) THEN
       CALL GetRealVector( Material, Load(4:6,1:n), 'Magnetization', FoundMagnetization )
     END IF

     CoilBody = .FALSE.
     CompParams => GetComponentParams( Element )

     CoilType = ''
     RotM = 0._dp
     ConstraintActive = .TRUE.
     IF (ASSOCIATED(CompParams)) THEN
       CoilType = GetString(CompParams, 'Coil Type', Found)
       IF (Found) THEN
         SELECT CASE (CoilType)
         CASE ('stranded')
            CoilBody = .TRUE.
         CASE ('massive')
            CoilBody = .TRUE.
         CASE ('foil winding')
            CoilBody = .TRUE.
            CALL GetElementRotM(Element, RotM, n)
         CASE DEFAULT
            CALL Fatal (Caller, 'Non existent Coil Type Chosen!')
         END SELECT
         ConstraintActive = GetLogical(CompParams, 'Activate Constraint', Found )
       END IF
     END IF

     LaminateStack = .FALSE.
     LaminateStackModel = ''
     HasTensorReluctivity = .FALSE.
     Acoef = 0.0d0
     Tcoef = 0.0d0
     IF ( ASSOCIATED(Material) ) THEN
       IF ( .NOT. ( HasHBCurve .OR. HasReluctivityFunction ) ) THEN
         CALL GetReluctivity(Material,Acoef_t,n,HasTensorReluctivity)
         IF (HasTensorReluctivity) THEN
           IF (SIZE(Acoef_t,1)==1 .AND. SIZE(Acoef_t,2)==1) THEN
             i = MIN(SIZE(Acoef), SIZE(Acoef_t,3))
             Acoef(1:i) = Acoef_t(1,1,1:i) 
             HasTensorReluctivity = .FALSE.
           END IF
         ELSE
           CALL GetReluctivity(Material,Acoef,n)
         END IF
         IF (HasTensorReluctivity) THEN
           IF (size(Acoef_t,1)/=3) CALL Fatal(Caller, &
               'Reluctivity tensor should be of size 3x3')
         END IF
       END IF
!------------------------------------------------------------------------------
!      Read conductivity values (might be a tensor)
!------------------------------------------------------------------------------
       Tcoef = GetElectricConductivityTensor(Element,n,'re',CoilBody,CoilType)       
       LaminateStackModel = GetString( Material, 'Laminate Stack Model', LaminateStack )
     END IF


     LamThick = 0.0_dp
     LamCond  = 0.0_dp
     IF (LaminateStack) THEN
       SELECT CASE(LaminateStackModel)
       CASE('low-frequency model')
         LamThick(1:n) = GetReal( Material, 'Laminate Thickness', Found )
         IF (.NOT. Found) CALL Fatal(Caller, 'Laminate Thickness not found!')

         LamCond(1:n) = GetReal( Material, 'Laminate Stack Conductivity', Found )
         IF (.NOT. Found) CALL Fatal(Caller, 'Laminate Stack Conductivity not found!')

       CASE DEFAULT
         CALL WARN(Caller, 'Nonexistent Laminate Stack Model chosen!')
       END SELECT
     END IF

     !Get element local matrix and rhs vector:
     !----------------------------------------
       CALL LocalMatrix( MASS, DAMP, STIFF, FORCE, JFixFORCE, JFixVec, LOAD, &
         Tcoef, Acoef, LaminateStack, LaminateStackModel, &
         LamThick, LamCond, CoilBody, CoilType, RotM, ConstraintActive, &
         Element, n, nd+nb, PiolaVersion, SecondOrder)
       
     ! Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     IF (Transient .OR. EigenSystem) THEN
       CALL DefaultUpdateMass(MASS)
       IF ( ElectroDynamics ) CALL DefaultUpdateDamp(DAMP)
     END IF

     ! Collect weak divergence constraint.
     !-----------------------------------------------------------------
     IF (Transient .AND. TransientGauge .AND. .NOT. TransientGaugeCollected) THEN
       CALL LocalConstraintMatrix(Element, n, nd+nb, PiolaVersion, SecondOrder)
     END IF

     CALL DefaultUpdateEquations(STIFF,FORCE)

     ! Memorize stuff for the fixing potential
     ! 1) Divergence of the source term
     ! 2) The source terms at the surface to determine the direction
     !-------------------------------------------------------------------
     IF( JFixSolve ) THEN
       JFixRhs(JFixVar % Perm(Element % NodeIndexes)) = &
           JFixRhs(JFixVar % Perm(Element % NodeIndexes)) + JFixFORCE(1:n)
       DO i=1,n
         j = JFixSurfacePerm(Element % NodeIndexes(i) )         
         IF( j > 0 ) JFixSurfaceVec(3*j-2:3*j) = &
             JFixSurfaceVec(3*j-2:3*j) + JFixVec(1:3,i)
       END DO
     END IF

   END DO

  CALL DefaultFinishBulkAssembly(BulkUpdate=ConstantBulk)


  ! If we are solving the fixing potential for this nonlinear iteration then
  ! add its contribution to the AV equation.
  IF( JFixSolve ) THEN    
    CALL Info(Caller,'Solving the fixing potential',Level=7)
    JFixPhase = 2
    
    IF( ListGetLogical( SolverParams,'Precomputed Fixing Term',Found ) ) THEN
      CALL Info(Caller,'Adding precomputed source term: g fix')
      Var => VariableGet( Mesh % Variables,'g fix')
      JFixRhs = JfixRhs + Var % Values
    END IF

    CALL JFixPotentialSolver(Model,Solver,dt,Transient)
   
    
    CALL Info(Caller,'Adding the fixing potential to the r.h.s. of AV equation',Level=10)   
    DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      nd = GetElementNOFDOFs()  
      nb = GetElementNOFBDOFs() 

      CALL LocalFixMatrix( FORCE, Element, n, nd+nb, PiolaVersion, SecondOrder)      
    END DO    
    CALL Info(Caller,'Finished adding the fixing potential',Level=10)   
  END IF


  ! This adds a precomputed source term to r.h.s. of the equation.
  ! Note that this is assumed to be already mapped to nodes. 
  IF( ListGetLogical( SolverParams,'Precomputed Source Term',Found ) ) THEN
    CALL Info(Caller,'Adding precomputed source term: g')
    Var => VariableGet( Mesh % Variables,'g')
    Solver % Matrix % Rhs = Solver % Matrix % Rhs + Var % Values
  END IF
  
  
100 CONTINUE

  !
  ! Robin type of BC in terms of H:
  !--------------------------------
  Active = GetNOFBoundaryElements()
  DO t=1,Active
     Element => GetBoundaryElement(t)
     BC=>GetBC()
     IF (.NOT. ASSOCIATED(BC) ) CYCLE

     SELECT CASE(GetElementFamily())
     CASE(1)
       CYCLE
     CASE(2)
       k = GetBoundaryEdgeIndex(Element,1); Element => Mesh % Edges(k)
     CASE(3,4)
       k = GetBoundaryFaceIndex(Element)  ; Element => Mesh % Faces(k)
     END SELECT
     IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

     Model % CurrentElement => Element
     nd = GetElementNOFDOFs(Element)
     n  = GetElementNOFNodes(Element)

     CALL GetRealVector( BC, Load(1:3,1:n), 'Magnetic Field Strength', Found )
     FoundAny = Found
     
     Acoef(1:n) = GetReal( BC, 'Magnetic Transfer Coefficient', Found ) 
     FoundAny = FoundAny .OR. Found
     
     Load(4,1:n) = GetReal( BC, 'Electric Current Density', Found )
     FoundAny = FoundAny .OR. Found 

     Load(5,1:n) = GetReal( BC, 'Electric Transfer Coefficient', Found )
     FoundAny = FoundAny .OR. Found

     Found = ListCheckPresent( BC, 'Perfect Electric Conductor' )
     FoundAny = FoundAny .OR. Found
 
     Found = ListCheckPresent( BC, 'Outflow BC' )
     FoundAny = FoundAny .OR. Found
     
     ThinLineCrossect = GetReal( BC, 'Thin Line Crossection Area', Found)

     IF (Found) THEN
       CALL Info("WhitneyAVSolver", "Found a Thin Line Element", level=10)
       ThinLineCond = GetReal(BC, 'Thin Line Conductivity', Found)
       IF (.NOT. Found) CALL Fatal('DoSolve','Thin Line Conductivity not found!')

       ! The basis function evaluation is not implemented for all basis types
       ! within LocalMatrixThinLine, check for the consistency:
       IF (.NOT. PiolaVersion) CALL Warn(Caller, &
           'The implementation of thin line element may need Use Piola Transform = True')

       CALL LocalMatrixThinLine(MASS,STIFF,FORCE,LOAD,ThinLineCrossect,ThinLineCond,Element,n,nd, &
           SecondOrder)
       CALL DefaultUpdateEquations(STIFF,FORCE,Element)
       IF (Transient .OR. EigenSystem) THEN
         IF(ElectroDynamics) THEN
           CALL DefaultUpdateDamp(MASS,Element)
         ELSE
           CALL DefaultUpdateMass(MASS,Element)
         END IF
       END IF
       CYCLE
     END IF
 
     !If air gap length keyword is detected, use air gap boundary condition
     GapLength=GetConstReal( BC, 'Air Gap Length', Found)
     IF (Found) THEN
       AirGapMu=GetConstReal( BC, 'Air Gap Relative Permeability', Found)
       IF (.NOT. Found) AirGapMu=1d0 ! if not found default to "air" property
       MASS = 0.0_dp
       CALL LocalMatrixAirGapBC(STIFF,FORCE,LOAD,GapLength,AirGapMu,Element,n,nd )
     ELSE IF(  FoundAny ) THEN
       CALL LocalMatrixBC(STIFF,MASS, FORCE,LOAD,Acoef,Element,n,nd )
     ELSE
       CYCLE
     END IF

     CALL DefaultUpdateEquations(STIFF,FORCE,Element)

     IF(ElectroDynamics) THEN
       CALL DefaultUpdateDamp(MASS,Element)
     ELSE
       CALL DefaultUpdateMass(MASS,Element)
     END IF
  END DO
  
  CALL DefaultFinishBoundaryAssembly(BulkUpdate=ConstantSystem)

  DoneAssembly = .TRUE.

  ! Check the timer
  CALL CheckTimer('MGDynAssembly', Delete=.TRUE.)
  
  
200 CONTINUE

  ! This is now automatically invoked as the time integration is set global in the Solver_init
  ! IF ( Transient.OR.EigenSystem) CALL Default1stOrderTimeGlobal()
  CALL DefaultFinishAssembly()

  ! Dirichlet BCs in terms of vector potential A:
  ! ---------------------------------------------
  IF ( TG ) THEN
    ! temporary fix to some scaling problem (to be resolved)...
    CALL ListAddLogical( SolverParams, 'Linear System Dirichlet Scaling', .FALSE.) 
  END IF


BLOCK
! Automatic BC for massive,foil coils outer boundaries, when "Activate Constraint" on!!

    TYPE(Element_t), POINTER :: Parent
    LOGICAL :: AutomaticBC
    INTEGER, POINTER :: Electrodes(:)

    A => GetMatrix()

    IF (.NOT.ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF

    IF(.NOT.ALLOCATED(A % DValues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF

    Active = GetNOFBoundaryElements()
    DO t = 1, Active
      Element => GetBoundaryElement(t)
      IF(.NOT.ActiveBoundaryElement()) CYCLE

      Parent => Element % BoundaryInfo % Right
      IF(ASSOCIATED(Parent)) CYCLE

      IF(ParEnv % PEs>1) THEN
        ! Assuming here that this is an internal boundary, if all elements nodes are
        ! interface nodes. Not foolproof i guess, but quite safe (?)
        IF (ALL(Solver % Mesh % ParallelInfo % GInterface(Element % NodeIndexes))) CYCLE
      END IF
 
      Parent => Element % BoundaryInfo % Left
      IF(.NOT.ASSOCIATED(Parent)) CYCLE

      CompParams => GetComponentParams(Parent)
      IF (.NOT. ASSOCIATED(CompParams)) CYCLE

      CoilType = GetString(CompParams, 'Coil Type', Found)
      IF(CoilType/='massive' .AND. CoilType/='foil winding') CYCLE

      ConstraintActive = GetLogical(CompParams,'Activate Constraint',Found )
      IF( .NOT. ConstraintActive ) CYCLE

      AutomaticBC = GetLogical( CompParams, 'Automatic electrode BC', Found )
      IF(.NOT.Found) AutomaticBC = .TRUE.

      IF(.NOT.AutomaticBC) CYCLE

      Electrodes =>  ListGetIntegerArray( CompParams, &
             'Electrode Boundaries', Found )

      IF(ASSOCIATED(Electrodes)) THEN
        IF(ALL(Electrodes/=Element % BoundaryInfo % Constraint)) CYCLE
      END IF

      DO i=1,Element % Type % NumberOfNodes
        j = Solver % Variable % Perm(Element % NodeIndexes(i))
        A % ConstrainedDOF(j) = .TRUE.
      END DO
    END DO
END BLOCK

  CALL DefaultDirichletBCs()

  ! Apply dirichlet BCs associated with weak divergence dofs
  IF (TransientGauge .AND. .NOT. TransientGaugeCollected) THEN
    CALL Info("WhitneyAVSolver", "Handling weak div boundary.", level=10)
    Active = GetNOFBoundaryElements()
    ELEMENT_LOOP: DO t = 1, Active
      Element => GetBoundaryElement(t)
      IF (.NOT. ActiveBoundaryElement(UElement=Element)) THEN
        CYCLE ELEMENT_LOOP
      END IF

      BC => GetBC()
      IF (.NOT. ASSOCIATED(BC)) THEN
        CYCLE ELEMENT_LOOP
      END IF

      IF (GetLogical ( BC, 'weak div dirichlet boundary', Found)) THEN
        n = GetElementNOFNodes(Element)
        DO i = 1, n
          j = Solver % Variable % Perm(Element % NodeIndexes(i))
          DO k=CM % Rows(j),CM % Rows(j+1)-1
            CM % Values(k) = 0.0_dp
            IF (CM % Cols(k) == j + A % NumberOfRows) CM % Values(k) = 1.0_dp
          END DO

          CM % RHS(j) = 0.0_dp

          k = CM % Diag(j)
          IF (k>0) THEN
            CM % Values(k) = 1.0_dp
          ELSE
            CALL Warn("WhitneyAVSolver", "there is no diagonal entry in constraint matrix.. Why is that?")
          END IF
        END DO
      END IF
    END DO ELEMENT_LOOP
    CALL Info("WhitneyAVSolver", "Done setting weak div dirichlet boundary.", level=10)

    CALL PackEdgeRows(CM, Model)
    CALL Info("WhitneyAVSolver", "Done removing unused constraints.", level=10)
    TransientGaugeCollected = .TRUE.
  END IF

  ! Dirichlet BCs in terms of magnetic flux density B:
  ! --------------------------------------------------
  CALL DirichletAfromB()
  CALL ConstrainUnused(A)

 
  IF (TG) THEN
    IF ( .NOT.ALLOCATED(TreeEdges) ) &
        CALL GaugeTree(Solver,Mesh,TreeEdges,FluxCount,FluxMap,Transient)
    CALL Info(Caller, 'Volume tree edges: '//i2s(COUNT(TreeEdges))// &
        ' of total: '//I2S(Mesh % NumberOfEdges),Level=5)
    
    DO i=1,SIZE(TreeEdges)
      IF(TreeEdges(i)) CALL SetDOFToValue(Solver,i,0._dp)
    END DO
  END IF

  IF (DefaultLineSearch(Converged)) RETURN
  IF ( Converged ) GOTO 10

  ! The following gives the user an option to adapt the linear system convergence tolerance
  ! adaptively:
  !-------------------------------------------------------------------------------------------
  AdaptiveTols = ListGetLogical(SolverParams, 'Linear System Adaptive Tolerance', Found)
  IF (AdaptiveTols) THEN
     IF (iterno == 1) THEN
        BaseTol = ListGetConstReal(SolverParams, 'Linear System Base Tolerance', Found)
        IF (Found) CALL ListAddConstReal(SolverParams, &
             'Linear System Convergence Tolerance', BaseTol)
     ELSE
        IF (.NOT. ListCheckPresent(SolverParams, 'Linear System Relative Tolerance')) THEN
           RelTol = 1.0d-2
        ELSE
           RelTol = ListGetConstReal(SolverParams, 'Linear System Relative Tolerance')
        END IF
        IF (.NOT. ListCheckPresent(SolverParams, 'Linear System Base Tolerance')) THEN
           BaseTol = 1.0d-3
        ELSE
           BaseTol = ListGetConstReal(SolverParams, 'Linear System Base Tolerance')
        END IF
        LinTol = MIN(BaseTol, RelTol * Solver % Variable % NonlinChange)
        CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance', &
             LinTol)
     END IF
  END IF

  norm = DefaultSolve()

  Converged = ( Solver % Variable % NonlinConverged == 1 )

  IF( UseTorqueTol ) THEN
    PrevTorque = Torque 
    CALL CalculateLumpedParameters(Torque)
    IF( iterNo >= 2 ) THEN
      TorqueErr = 2 * ABS(PrevTorque-Torque) / (ABS(PrevTorque)+ABS(Torque))
      IF( TorqueErr > TorqueTol ) THEN
        WRITE(Message,'(A,ES12.3)') 'Torque error at iteration '//I2S(iterNo)//':',TorqueErr
        CALL Info(Caller,Message,Level=6)
      END IF
    END IF
  END IF
    
  CALL Info(Caller,'Convergence status: '//I2S(Solver % Variable % NonlinConverged),Level=12)
  IF( Converged ) THEN
    CALL Info(Caller,'System has converged to tolerances after '//I2S(iterNo)//' iterations!',Level=12)
    IF( UseTorqueTol ) THEN
      IF( TorqueErr > TorqueTol ) THEN
        CALL Info(Caller,'Nonlinear system tolerance ok after '&
            //I2S(iterNo)//' but torque still wobbly!',Level=7)
        Converged = .FALSE.
      END IF
    END IF
  END IF
  
  IF( ListGetLogical( SolverParams,'Calculate Magnetic Norm',Found ) .OR. &
      ListGetLogical( SolverParams,'Use Magnetic Norm', Found ) ) THEN
    BLOCK
      REAL(KIND=dp) :: binteg, bmin, bmax, vinteg, bnorm

      binteg = 0.0_dp
      vinteg = 0.0_dp
      bmin = HUGE( bmin )
      bmax = -HUGE( bmax ) 

      Active = GetNOFActive()
      DO t=1,active
        Element => GetActiveElement(t)

        IF( ParEnv % PEs > 1 ) THEN
          IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
        END IF

        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs() 
        nb = GetElementNOFBDOFs()

        CALL AddLocalBNorm( Element, n, nd+nb, PiolaVersion, SecondOrder, binteg, vinteg, bmin, bmax)
      END DO

      IF( ParEnv % PEs > 1 ) THEN
        binteg = ParallelReduction( binteg )
        vinteg = ParallelReduction( vinteg ) 
        bmin = ParallelReduction( bmin,1 )
        bmax = ParallelReduction( bmax,2 )      
      END IF

      ! We used square to avoid taking root every time. 
      bmin = SQRT(bmin)
      bmax = SQRT(bmax)
      bnorm = SQRT(binteg/vinteg)
            
      WRITE( Message,'(A,ES15.6)') 'Magnetic field norm:',bnorm
      CALL Info(Caller, Message, Level=4) 
      
      CALL ListAddConstReal( Model % Simulation,'res: magnetic norm',bnorm )
      
      WRITE( Message,'(A,ES15.6)') 'Magnetic field minimum value:',bmin
      CALL Info(Caller, Message, Level=8 ) 

      WRITE( Message,'(A,ES15.6)') 'Magnetic field maximum value:',bmax
      CALL Info(Caller, Message, Level=8 ) 

      IF( ListGetLogical( SolverParams,'Use Magnetic Norm',Found ) ) THEN
        CALL Info(Caller,'Setting solver norm to magnetic norm!')
        Solver % Variable % Norm = bnorm
      END IF
      
    END BLOCK
  END IF
    

  
10 CONTINUE

  IF ( ALLOCATED(FluxMap) ) DEALLOCATE(FluxMap)
! IF ( ALLOCATED(TreeEdges) ) DEALLOCATE(TreeEdges)

! CALL WriteResults  ! debugging helper


!------------------------------------------------------------------------------
 END FUNCTION DoSolve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE ConstrainUnused(A)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(Matrix_t) :: A
!------------------------------------------------------------------------------
  INTEGER :: i,j,n

  REAL(kind=dp), ALLOCATABLE :: dDiag(:)
!------------------------------------------------------------------------------

  n = A % NumberOfRows
  ALLOCATE(dDiag(n)); dDiag = 0._dp

  DO i=1,n
    j = A % Diag(i)
    IF (j>0) dDiag(i) = A % Values(j)
  END DO
  IF (ParEnv % PEs>1) CALL ParallelSumVector(A, dDiag)

  n = SIZE(Solver % Variable % Perm)
  DO i=1,n
    J = Solver % Variable % Perm(i)
    IF(j>0) THEN
      IF ( dDiag(j)==0._dp ) THEN
        CALL ZeroRow(A,j)
        CALL SetMatrixElement(A,j,j,1._dp)
        A % RHS(j)=0._dp
        IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(j)=.TRUE.
      END IF
    END IF
  END DO
!------------------------------------------------------------------------------
 END SUBROUTINE ConstrainUnused
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumpedSurface()
!------------------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------------------
   REAL(KIND=dp) :: zzforce
   INTEGER :: i,n,nd,EdgeBasisDegree
   LOGICAL :: Found
   TYPE(ELement_t), POINTER :: Element, Parent
!------------------------------------------------------------------------------

   EdgeBasisDegree = 1
   IF (SecondOrder) EdgeBasisDegree = 2

   zzforce = 0
   DO i=1,Mesh % NumberOFBoundaryElements
     Element => GetBoundaryElement(i)
     IF (.NOT.GetLogical(GetBC(), 'Calculate Axial Force', Found ) ) CYCLE

     Parent => Element % BoundaryInfo % Left
     n  = GetELementNofNodes(Parent)
     nd = GetELementNofDOFs(Parent)
     CALL AxialForceSurf(zzforce,Element,n,nd,EdgeBasisDegree)
   END DO

   zzforce = ParallelReduction(zzforce)
   CALL ListAddConstReal(Model % Simulation,'res: Axial force(surf)', zzforce )

!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumpedSurface
!------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
! This is monolithic lumping routine copied from the 2D solver in the hope that
! we could use exactly same keywords and logic for torque computation etc.
!------------------------------------------------------------------------------
 SUBROUTINE CalculateLumpedParameters(Torque)
!------------------------------------------------------------------------------
   REAL(KIND=dp), OPTIONAL :: Torque

   REAL(KIND=dp) :: torq,TorqArea,IMoment,IA, &
       rinner,router,rmean,rdiff,ctorq,detJ,Weight,&
       Bp,Br,Bx,By,Bz,x,y,z,r,rho,zmin,zmax,AxialForce
   REAL(KIND=dp), ALLOCATABLE :: a(:),u(:),POT(:),dPOT(:), &
       pPot(:),Density(:)
   REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:), RotWBasis(:,:), WBasis(:,:)
   LOGICAL, ALLOCATABLE :: TorqueElem(:)
   INTEGER :: i,bfid,n,nd,nbf,PrevComm,NoSlices,NoTimes,np
   LOGICAL :: Found, Stat
   TYPE(ValueList_t),POINTER::Params
   TYPE(GaussIntegrationPoints_t) :: IP
   TYPE(Nodes_t) :: Nodes
   LOGICAL :: CalcTorque, CalcPot, CalcInert
   LOGICAL :: ThisTorque, ThisPot, ThisInert, ThisAxial, Parallel, HaveRange
   LOGICAL :: Visited = .FALSE.
   INTEGER :: EdgeBasisDegree, dim
   LOGICAL :: SliceAverage, WbaseFound
   REAL(KIND=dp) :: W(3)
   REAL(KIND=dp), POINTER :: Wpot(:)
   
   SAVE Visited, Nodes, Basis, dBasisdx, a, u, POT, dPOT, pPot, &
       Density, Ctorq, TorqueElem, RotWBasis, WBasis, zmin, zmax, &
       Wpot
   
!------------------------------------------------------------------------------

   CALL Info(Caller,'Calculating lumped parameters',Level=8)

   ! By construction
   dim = 3

   IF( dim == 3 ) THEN
     NoSlices = 1
     SliceAverage = .FALSE.
   ELSE    
     NoSlices = MAX(1,ListGetInteger( Model % Simulation,'Number Of Slices', SliceAverage ) )
   END IF

   NoTimes = ListGetInteger( Model % Simulation,'Number Of Times', Found )
   IF( NoTimes > 1 ) THEN
     PrevComm = ParEnv % ActiveComm
     ParEnv % ActiveComm = ParallelSlicesComm() 
   END IF
      
   ! Define whether we have something to compute
   ! Old way is to give radius range.
   ! New way uses rotor radius + neighboring element layer. 
   rinner = ListGetCRealAnyBody( Model,'r inner',CalcTorque )
   IF( CalcTorque ) THEN
     router = ListGetCRealAnyBody( Model,'r outer')
     rmean = (rinner+router)/2
     rdiff = (router-rinner)
     HaveRange = .TRUE.     
   ELSE
     rmean = ListGetConstReal( CurrentModel % Simulation,'Rotor Radius',CalcTorque)
     rdiff = ListGetConstReal( CurrentModel % Simulation,'Rotor Air Gap Width',Found)
     IF(.NOT. Found ) rdiff = 1.0e-3 * rmean
     HaveRange = .FALSE.
   END IF
   
   CalcPot = ListGetLogicalAnyBodyForce( Model,'Calculate Potential' )
   CalcInert = CalcTorque .AND. .NOT. Visited 
   
   
   IF( PRESENT(Torque) .AND. .NOT. CalcTorque ) THEN
     CALL Fatal(Caller,'Torque tolerance requested, but torque not computed!')
   END IF
        
   IF(.NOT. (CalcTorque .OR. CalcPot .OR. CalcInert) ) RETURN

   Parallel = ( ParEnv % PEs > 1 )
   
   nbf = Model % NumberOfBodyForces
   IF(.NOT. Visited ) THEN
     n = Model % Mesh % MaxElementDofs
     ALLOCATE( a(nbf), u(nbf), POT(n), dPOT(n), pPot(n), Density(n), &
         Basis(n), dBasisdx(n,3), RotWBasis(n,3),WBasis(n,3),Wpot(n))
   END IF

   IF( CalcTorque ) THEN
     torq = 0._dp
     TorqArea = 0._dp
     AxialForce = 0.0_dp
   END IF
   IF( CalcInert ) THEN
     IMoment = 0._dp
     IA = 0.0_dp
   END IF
   IF( CalcPot ) THEN   
     U=0._dp
     a=0._dp
   END IF

   IF (SecondOrder) THEN
     EdgeBasisDegree = 2
   ELSE
     EdgeBasisDegree = 1
   END IF
   
   IF(.NOT. Visited .AND. CalcTorque ) THEN
     zmin = HUGE(zmin)
     zmax = -HUGE(zmax)
     ALLOCATE( TorqueElem( GetNOFActive() ) )
     TorqueElem = .FALSE.
     
     DO i=1,GetNOFActive()
       Element => GetActiveElement(i)
     
       ThisTorque = .FALSE.
       
       n  = GetElementNOFNodes(Element)     
       CALL GetElementNodes( Nodes, Element )
       
       IF( HaveRange ) THEN       
         ! We are given range in classical Arkkio style. 
         ! Check how the center lies with respect to the range.
         x = SUM(Nodes % x(1:n))/n
         y = SUM(Nodes % y(1:n))/n
         r = SQRT(x**2+y**2)
         IF (r >= rinner .AND. r <= router) THEN
           TorqueElem(i) = .TRUE.
         END IF
       ELSE
         ! We are not given a range. Just take any element
         ! which has even one node at the given radius. 
         DO j=1,n
           x = Nodes % x(j)
           y = Nodes % y(j)
           r = SQRT(x**2+y**2)
           IF( ABS(r-rmean) < rdiff / 2 ) THEN
             TorqueElem(i) = .TRUE.
             EXIT
           END IF
         END DO
       END IF
       
       IF(TorqueElem(i)) THEN
         zmin = MIN(MINVAL(Nodes % z(1:n)), zmin)
         zmax = MAX(MAXVAL(Nodes % z(1:n)), zmin)
       END IF
     END DO
              
     i = COUNT( TorqueElem )
     i = ParallelReduction(i) 
     CALL Info(Caller,'Number of elements to compute torque: '//I2S(i))

     IF( dim == 3 ) THEN
       zmin = ParallelReduction(zmin,1)
       zmax = ParallelReduction(zmax,2)     
       WRITE(Message,'(A,2ES12.3)') 'Torque compute range in z-direction',zmin,zmax
       CALL Info(Caller,Message)
     END IF
   END IF

   
   DO i=1,GetNOFActive()
     Element => GetActiveElement(i)
     n = GetElementNOFNodes(Element)     
     
     ThisTorque = .FALSE.
     ThisPot = .FALSE.
     ThisInert = .FALSE.
     
     IF( CalcPot ) THEN
       Params => GetBodyForce(Element)
       IF(ASSOCIATED(Params)) THEN         
         ThisPot = GetLogical(Params,'Calculate Potential',Found)
         IF( ThisPot ) THEN
           bfid = GetBodyForceId(Element)           
           CALL GetLocalSolution(POT, UElement=Element)
           CALL GetLocalSolution(pPOT,tstep=-1,UElement=Element)
           IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
             dPot = (POT - pPOT)/dt
           ELSE
             dPot = 1.5_dp*POT - 2*pPOT
             CALL GetLocalSolution(pPOT,tstep=-2,UElement=Element)
             dPot = (dPOT + 0.5_dp*pPOT)/dt
           END IF
           IF( dim == 3 ) THEN
             CALL GetLocalSolution(Wpot,'W',UElement=Element)
             W = [0._dp, 0._dp, 1._dp]
             WbaseFound = ANY(ABS(Wpot(1:n)) > 1.0e-20 )
           END IF
         END IF
       END IF
     END IF
     
     IF( CalcInert ) THEN
       Params=>GetBodyParams(Element)
       IF(ASSOCIATED(Params)) THEN
         ThisInert = GetLogical(Params,'Calculate Inertial Moment',Found)
       END IF
       Density(1:n) = GetReal(GetMaterial(),'Density',Found,Element)
     END IF
     
     IF( CalcTorque ) THEN
       ThisTorque = TorqueElem(i)
       IF(ThisTorque .AND. .NOT. ThisPot ) THEN
         CALL GetLocalSolution(POT, UElement=Element)
       END IF
     END IF

     ! Only treat the element if we have something to compute
     IF( .NOT. (ThisPot .OR. ThisInert .OR. ThisTorque ) ) CYCLE
       
     nd = GetElementNOFDOFs(Element)
     IF(dim==3) THEN
       np = n*pSolver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
     END IF
     CALL GetElementNodes( Nodes, Element )

     ! Numerical integration:
     !-----------------------
     IF( dim == 3 ) THEN
       IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
           EdgeBasisDegree=EdgeBasisDegree)       
     ELSE
       IP = GaussPoints(Element)
     END IF

     DO t=1,IP % n       
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       IF( dim == 3 ) THEN
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
             RotBasis = RotWBasis, USolver = pSolver )
       ELSE
         stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
             IP % W(t), detJ, Basis, dBasisdx )
       END IF
       weight = IP % s(t) * detJ

       ! Coordinates of the integration point
       x = SUM(Nodes % x(1:n)*Basis(1:n))
       y = SUM(Nodes % y(1:n)*Basis(1:n))
       r = SQRT(x**2+y**2)
       
       IF(ThisPot ) THEN
         A(bfid) = A(bfid) + Weight
         IF( dim == 3 ) THEN
           IF( WbaseFound ) W = MATMUL(Wpot(1:n),dBasisdx(1:n,:))
           U(bfid) = U(bfid) + Weight * SUM(dPot(np+1:nd)*MATMUL(WBasis(1:nd-np,:),w))
         ELSE
           U(bfid) = U(bfid) + Weight * SUM(dPot(1:nd)*Basis(1:nd))
         END IF
       END IF
       
       IF( ThisTorque ) THEN                      
         IF( dim == 3 ) THEN
           Bx =  SUM(POT(np+1:nd) * RotWBasis(1:nd-np,1))
           By =  SUM(POT(np+1:nd) * RotWBasis(1:nd-np,2))
         ELSE           
           Bx =  SUM(POT(1:nd)*dBasisdx(1:nd,2))
           By = -SUM(POT(1:nd)*dBasisdx(1:nd,1))
         END IF

         Br =  x/r*Bx + y/r*By
         Bp = -y/r*Bx + x/r*By

         Torq = Torq + Weight * r*Br*Bp / (PI*4.0d-7*rdiff)
         TorqArea = TorqArea + Weight

         IF( dim == 3 ) THEN
           Bz =  SUM(POT(np+1:nd) * RotWBasis(1:nd-np,3))
           AxialForce = AxialForce + Weight*1*(Bx*Bz*x/r+By*Bz*y/r)/(PI*4.0d-7*rdiff)
         END IF
       END IF
       
       IF( ThisInert ) THEN
         IF( r < rmean ) THEN
           rho = SUM( density(1:n) * Basis(1:n) ) 
           IF( rho > EPSILON( rho ) ) THEN
             IA = IA + Weight
             U = U + Weight * r * rho
           END IF
         END IF
       END IF
     END DO
   END DO
     
   ! Finally perform parallel reduction if needed, and
   ! store the results for saving by SaveScalars.
   !-------------------------------------------------------------------------   
   IF( CalcPot ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       DO i=1,nbf
         a(i) = ParallelReduction(a(i)) / NoSlices
         u(i) = ParallelReduction(u(i)) / NoSlices
       END DO
     END IF     
     DO i=1,nbf
       IF(a(i)>0) THEN
         CALL ListAddConstReal(Model % Simulation,'res: Potential / bodyforce ' &
             //i2s(i),u(i)/a(i))
         CALL ListAddConstReal(Model % Simulation,'res: area / bodyforce ' &
             //i2s(i),a(i))
       END IF
     END DO
   END IF
   
   IF( CalcTorque ) THEN   
     ! Arkkios formula assumes that rinner and router are nicely aligned with elements.
     ! This may not the case, so the 1st time we make a geomeric correction. 
     IF(.NOT. Visited ) THEN
       WRITE(Message,'(A,ES15.4)') 'Air gap initial torque:', Torq
       CALL Info(Caller,Message,Level=6)

       TorqArea = ParallelReduction(TorqArea) / NoSlices       
       IF (TorqArea > EPSILON(TorqArea) ) THEN
         IF( dim == 3 ) THEN
           Ctorq = 2 * PI * rmean * rdiff * (zmax-zmin) / TorqArea
         ELSE
           Ctorq = 2 * PI * rmean * rdiff / TorqArea
         END IF
         
         WRITE(Message,'(A,F8.4)') 'Air gap correction initial:', cTorq
         CALL Info(Caller,Message,Level=4)
         
         ! The correction factor also corrects for the number of periods.
         ! We don't want that - so let us take back that and the torque
         ! can be compared to inertial moment of the sector still. 
         i = ListGetInteger( CurrentModel % Simulation,'Rotor Periods',Found )         
         IF( Parallel ) i = ParallelReduction( i, 2 ) 
         IF( i > 1 ) THEN
           WRITE(Message,'(A,I0)') 'Air gap correction rotor periods: ',i
           CALL Info(Caller,Message,Level=4)
           Ctorq = Ctorq / i 
         END IF         
       ELSE
         Ctorq = 1.0_dp
       END IF
       
       WRITE(Message,'(A,F8.4)') 'Air gap correction:', cTorq
       CALL Info(Caller,Message,Level=4)
       !CALL ListAddConstReal(Model % Simulation,'res: air gap correction', cTorq)
     END IF
       
     Torq = Ctorq * Torq
     
     IF( SliceAverage ) THEN
       ! Save slice torque even for one slice since then the output for scalars is the same
       ! for any number of slices.       
       WRITE(Message,'(A,ES15.4)') 'Air gap torque for slice'//I2S(ParEnv % MyPe)//':', Torq
       CALL Info(Caller,Message,Level=5)
       CALL ListAddConstReal(Model % Simulation,'res: air gap torque for slice', Torq)
     END IF
       
     ! But the averaging makes sense only for more than one slice
     Torq = ParallelReduction(Torq) / NoSlices
     WRITE(Message,'(A,ES15.4)') 'Air gap torque:', Torq
     CALL Info(Caller,Message,Level=5)
     CALL ListAddConstReal(Model % Simulation,'res: air gap torque', Torq)

     IF( dim == 3) THEN
       AxialForce = Ctorq * AxialForce 
       AxialForce = ParallelReduction(AxialForce)
       CALL ListAddConstReal(Model % Simulation,'res: Axial Force(vol)', AxialForce )
     END IF
           
     IF(PRESENT(Torque)) Torque = Torq

     ! This is placed here so that the legacy order of "res:" stuff is honored. 
     IF(ListGetLogicalAnyBC(Model,'Calculate Axial Force')) THEN
       CALL CalculateLumpedSurface()
     END IF
   END IF

     
   IF( CalcInert ) THEN
     IF( Parallel ) THEN
       IMoment = ParallelReduction(IMoment) / NoSlices
       IA = ParallelReduction(IA) / NoSlices
     END IF

     WRITE(Message,'(A,ES15.4)') 'Inertial volume:', IA
     CALL Info(Caller,Message,Level=7)

     WRITE(Message,'(A,ES15.4)') 'Inertial moment:', Imoment
     CALL Info(Caller,Message,Level=7)
       
     CALL ListAddConstReal(Model % Simulation,'res: inertial volume', IA)
     CALL ListAddConstReal(Model % Simulation,'res: inertial moment', IMoment)
   END IF
     
   Visited = .TRUE.
   
   ! Revert the communicatior back to original
   IF( NoTimes > 1 ) ParEnv % ActiveComm = PrevComm
   
!------------------------------------------------------------------------------
 END SUBROUTINE CalculateLumpedParameters
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE AxialForceSurf(U,Element,n,nd,EdgeBasisDegree)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n,nd,EdgeBasisDegree
    REAL(KIND=dp)::U
    TYPE(Element_t)::Element
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent
    REAL(KIND=dp) :: dBasisdx(nd,3),Basis(nd), DetJ, Pdetj, uu,v,w, &
             POT(nd),x,y,r,r0,r1,Wbasis(nd,3),RotWBasis(nd,3)
    REAL(KIND=dp) :: B(3,nd), Bx, By, Bz
    INTEGER :: t
    LOGICAL :: stat, Found
    TYPE(Nodes_t), SAVE :: Nodes, PNodes
    TYPE(GaussIntegrationPoints_t) :: IP
	!$OMP THREADPRIVATE(Nodes)

    CALL GetElementNodes( Nodes, Element )
    Parent => Element % BoundaryInfo % Left
    CALL GetElementNodes( PNodes, Parent )

    CALL GetLocalSolution(POT, UElement=Parent )
  
    !Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detJ, Basis, dBasisdx )
      CALL GetParentUVW(Element,GetElementNOFNodes(Element),Parent,n,uu,v,w,Basis)
      
      IF (PiolaVersion) THEN
        stat = EdgeElementInfo( Parent, PNodes, uu, v, w, &
              DetF = PDetJ, Basis = Basis, EdgeBasis = WBasis, RotBasis = RotWBasis, &
              BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
      ELSE
        stat = ElementInfo( Parent, PNodes, uu,v,w, pdetJ, Basis, dBasisdx )
        CALL GetEdgeBasis(Parent,WBasis,RotWBasis,Basis,dBasisdx)
      END IF

      x = SUM(Basis(1:n) * PNodes % x(1:n))
      y = SUM(Basis(1:n) * PNodes % y(1:n))
      r = SQRT(x**2 + y**2)
      x=x/r; y=y/r

      Bx =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,1)) 
      By =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,2))
      Bz =  SUM(POT(n+1:nd) * RotWBasis(1:nd-n,3))
      U = U + IP % s(t) * detJ * (Bx*Bz + By*Bz) /(PI*4.0d-7) !/ 2
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AxialForceSurf
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE LocalConstraintMatrix( Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
  IMPLICIT NONE

  TYPE(Element_t), POINTER :: Element
  INTEGER :: n, nd
  LOGICAL :: PiolaVersion, SecondOrder

  ! FE-Basis stuff
  REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
  REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
  REAL(KIND=dp) :: Permittivity(n), P_ip, MASS(nd,nd), DAMP(nd,nd), STIFF(nd,nd), FORCE(nd)

  INTEGER :: t, i, j, p, q, np, EdgeBasisDegree, r, s, Indexes(1:nd)
  TYPE(GaussIntegrationPoints_t) :: IP

  LOGICAL :: ElectroDynamics

  TYPE(Matrix_t), POINTER :: A
  REAL(KIND=dp), POINTER :: SaveRHS(:), SaveValues(:)

  TYPE(Nodes_t), SAVE :: Nodes
  !------------------------------------------------------------------------------
  IF (SecondOrder) THEN
    EdgeBasisDegree = 2
  ELSE
    EdgeBasisDegree = 1
  END IF

  CALL GetElementNodes( Nodes )

  STIFF = 0
  MASS  = 0
  DAMP  = 0
  FORCE = 0

  ElectroDynamics = ListGetLogical( SolverParams, 'Electrodynamics model', Found )
  IF ( ElectroDynamics ) THEN
    CALL GetPermittivity(GetMaterial(), Permittivity, n)
  END IF

  !Numerical integration:
  !----------------------
  IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
    EdgeBasisDegree=EdgeBasisDegree )

  np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
  DO t=1,IP % n

    stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
        RotBasis = RotWBasis, USolver = pSolver )

    IF ( TransientGauge ) THEN
      IF ( ElectroDynamics ) THEN
        P_ip = SUM( Basis(1:n) * Permittivity(1:n) )
        DO j = 1, np
          r = j
          DO i = 1,nd-np
            s = i+np
            DAMP(r,s) = DAMP(r,s) + P_ip*SUM(dBasisdx(j,:)*WBasis(i,:))*detJ*IP % s(t)
          END DO
          DO i = 1, np
            s = i
            STIFF(r,s) = STIFF(r,s) + P_ip*SUM(dBasisdx(j,:)*dBasisdx(i,:))*detJ*IP%s(t)
          END DO
        END DO
      ELSE
        DO j = 1, np
          r = j
          DO i = 1,nd-np
            s = i+np
            STIFF(r,s) = STIFF(r,s) + SUM(dBasisdx(j,:)*WBasis(i,:))*detJ*IP % s(t)
          END DO
          DO i = 1, np
            s = i
            STIFF(r,s) = STIFF(r,s) + gauge_penalize_c*SUM(dBasisdx(j,:)*dBasisdx(i,:))*detJ*IP%s(t) &
                 + gauge_penalize_m*Basis(j)*Basis(i)*detJ*IP%s(t)
          END DO
        END DO
      END IF
    END IF
  END DO

  IF(ElectroDynamics) CALL Default2ndOrderTimeR( MASS, DAMP, STIFF, FORCE )

  A => GetMatrix()

  SaveRHS  => A % RHS
  SaveValues => A % Values

  A % RHS => CM % RHS
  A % Values => CM % Values

  CALL DefaultUpdateEquations(STIFF, FORCE)

  A % RHS => SaveRHS
  A % Values => SaveValues

  p = GetElementDOFs(Indexes)
  DO j = 1, nd-np
    s = j+np
    CM % ConstrainedDOF(Solver % Variable % Perm(indexes(s))) = .TRUE.
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE LocalConstraintMatrix
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( MASS, DAMP, STIFF, FORCE, JFixFORCE, JFixVec, LOAD, &
            Tcoef, Acoef, LaminateStack, LaminateStackModel, &
            LamThick, LamCond, CoilBody, CoilType, RotM, ConstraintActive, &
            Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:), DAMP(:,:), JFixFORCE(:), JFixVec(:,:)
    REAL(KIND=dp) :: LOAD(:,:), Tcoef(:,:,:), Acoef(:), &
                     LamThick(:), LamCond(:)
    LOGICAL :: LaminateStack, CoilBody, ConstraintActive
    CHARACTER(LEN=MAX_NAME_LEN):: LaminateStackModel, CoilType
    REAL(KIND=dp) :: RotM(3,3,n)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Aloc(nd), JAC(nd,nd), mu, muder, B_ip(3), Babs
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), C(3,3), &
                     RotMLoc(3,3), velo(3), omega(3), omega_velo(3,n), &
                     lorentz_velo(3,n), VeloCrossW(3), RotWJ(3), CVelo(3), &
                     A_t(3,3), A_t_der(3,3), eps=1.0e-3, Permittivity(nd), P_ip
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ, L(3), G(3), M(3), JFixPot(nd)
    REAL(KIND=dp) :: LocalLamThick, LocalLamCond, CVeloSum
    REAL(KIND=dp), POINTER :: MuTensor(:,:)
    LOGICAL :: Stat, Found, HasVelocity, HasLorentzVelocity, HasAngularVelocity, LocalGauge
    INTEGER :: t, i, j, k, p, q, np, EdgeBasisDegree, mudim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp
    DAMP  = 0.0_dp
    JAC = 0.0_dp

    IF( JFix ) THEN
      ! If we are solving for the JFix field we cannot yet use it!
      ! This happens on the first iteration
      IF ( JFixSolve ) THEN
        JFixVec   = 0.0_dp
        JFixFORCE = 0.0_dp
      ELSE        
        JFixPot(1:n) = JFixVar % Values(JFixVar % Perm(Element % NodeIndexes))
      END IF
    END IF
      
    HasVelocity = .FALSE.
    LocalGauge = .FALSE.
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, omega_velo, 'Angular velocity', HasAngularVelocity)
      CALL GetRealVector( BodyForce, lorentz_velo, 'Lorentz velocity', HasLorentzVelocity)
      LocalGauge = GetLogical( BodyForce,'Local Lagrange Gauge', Found ) 
      HasVelocity = HasAngularVelocity .OR. HasLorentzVelocity
    END IF

    CALL GetPermittivity(GetMaterial(), Permittivity, n)

    IF ( HasHBCurve .OR. HasReluctivityFunction ) THEN
      CALL GetScalarLocalSolution(Aloc)
    END IF
    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
      
    !Numerical integration:
    !----------------------
    IP = GaussPointsAdapt(Element, Solver, EdgeBasis=.TRUE. )
!    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
!        EdgeBasisDegree=EdgeBasisDegree )

    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )

       IF ( HasHBCurve ) THEN
         B_ip = MATMUL( Aloc(np+1:nd), RotWBasis(1:nd-np,:) )
         babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )

         ! h-b
         IF( Newton ) THEN
           mu = ListGetFun( Material,'h-b curve',babs,dFdx=muder) / Babs
           muder = (muder-mu)/babs
         ELSE
           mu = ListGetFun( Material,'h-b curve',babs) / Babs
         END IF
       ELSE IF( HasReluctivityFunction ) THEN
         B_ip = MATMUL( Aloc(np+1:nd), RotWBasis(1:nd-np,:) )        
         babs = MAX( SQRT(SUM(B_ip**2)), 1.d-8 )
         mu = ListGetElementReal( mu_h, Basis, Element, &
             GaussPoint = t, Rdim=mudim, Rtensor=MuTensor, DummyVals = B_ip )
         IF (mudim < 2) CALL Fatal(Caller, &
             'Specify Reluctivity Function as a full (3x3)-tensor')
         A_t(1:3,1:3) = muTensor(1:3,1:3)         
         
         IF( Newton ) THEN
           ! Use central differencing 
           mu = ListGetElementReal( mu_h, Basis, Element, &
               GaussPoint = t, Rdim=mudim, Rtensor=MuTensor, DummyVals = (1+newton_eps)*B_ip )
           A_t_der(1:3,1:3) = muTensor(1:3,1:3)
           mu = ListGetElementReal( mu_h, Basis, Element, &
               GaussPoint = t, Rdim=mudim, Rtensor=MuTensor, DummyVals = (1-newton_eps)*B_ip )
           A_t_der(1:3,1:3) = ( A_t_der(1:3,1:3) - muTensor(1:3,1:3) ) / ( 2*newton_eps*babs)
         END IF
         
       ELSE IF( HasTensorReluctivity ) THEN
         IF (SIZE(Acoef_t,2) == 1) THEN
           A_t = 0.0d0
           DO i = 1,3
             A_t(i,i) = SUM(Basis(1:n)*Acoef_t(i,1,1:n))
           END DO
         ELSE
           DO i = 1,3
             DO j = 1,3
               A_t(i,j) = SUM(Basis(1:n)*Acoef_t(i,j,1:n))
             END DO
           END DO
         END IF
       ELSE
         mu = SUM( Basis(1:n) * Acoef(1:n) )
       END IF

       ! Compute convection type term coming from a rigid motion:
       ! --------------------------------------------------------
       IF (HasVelocity) THEN
         velo = 0.0_dp
         IF( HasAngularVelocity ) THEN
           omega(1) = SUM(basis(1:n)*omega_velo(1,1:n))
           omega(2) = SUM(basis(1:n)*omega_velo(2,1:n))
           omega(3) = SUM(basis(1:n)*omega_velo(3,1:n))
           DO i=1,n
             velo(1:3) = velo(1:3) + CrossProduct(omega_velo(1:3,i), [ &
                 basis(i) * Nodes % x(i), &
                 basis(i) * Nodes % y(i), &
                 basis(i) * Nodes % z(i)])
           END DO
         END IF
         IF( HasLorentzVelocity ) THEN
           velo(1:3) = velo(1:3) + [ &
               SUM(basis(1:n)*lorentz_velo(1,1:n)), &
               SUM(basis(1:n)*lorentz_velo(2,1:n)), &
               SUM(basis(1:n)*lorentz_velo(3,1:n))]
         END IF
       END IF

       ! Compute the conductivity tensor
       ! -------------------------------
       DO i=1,3
         DO j=1,3
           C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       ! Transform the conductivity tensor (in case of a foil winding):
       ! --------------------------------------------------------------
       IF (CoilBody .AND. CoilType /= 'massive') THEN
         DO i=1,3
           DO j=1,3
             RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )
           END DO
         END DO
         C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))
       END IF

       P_ip = SUM( Permittivity(1:n) * Basis(1:n) )

       M = MATMUL( LOAD(4:6,1:n), Basis(1:n) )
       L = MATMUL( LOAD(1:3,1:n), Basis(1:n) )
         
       LocalLamThick = SUM( Basis(1:n) * LamThick(1:n) )
       LocalLamCond = SUM( Basis(1:n) * LamCond(1:n) )

       ! Add -C * grad(V^s), where C is a tensor
       ! -----------------------------------------
       L = L-MATMUL(C, MATMUL(LOAD(7,1:n), dBasisdx(1:n,:)))

       IF( JFix ) THEN
         IF( JFixSolve ) THEN
           ! If we haven't solved for the disbalance of source terms assemble it here
           DO i = 1,n
             p = i
             JFixFORCE(p) = JFixFORCE(p) + SUM(L * dBasisdx(i,:)) * detJ * IP%s(t) 
             JFixVec(:,p) = JFixVec(:,p) + L * Basis(i) * detJ * IP%s(t)
           END DO
         ELSE         
           ! If we have already solved for the JFix potential use it here
           L = L - MATMUL(JFixPot(1:n), dBasisdx(1:n,:))
         END IF
       END IF
             
       ! ------------------------------------------------------------------
       ! Compute element stiffness matrix and force vector.
       ! If we calculate a coil, the nodal degrees of freedom are not used.
       ! ------------------------------------------------------------------
       IF ( ConstraintActive ) THEN
         ! --------------------------------------------------------
         ! The constraint equation involving the scalar potential:
         !     -div(C*(dA/dt+grad(V)))=0
         ! All terms that are added here depend on the electrical conductivity or permittivity,
         ! so they have an effect on a conductor only.
         ! --------------------------------------------------------
         CONDUCTOR: IF ( SUM(ABS(C)) > AEPS .OR. ElectroDynamics ) THEN
           IF ( Transient.OR.EigenSystem ) THEN
             DO p=1,np
               DO q=1,np

                 ! Compute the conductivity term <C grad V,grad v> for stiffness 
                 ! matrix (anisotropy taken into account)
                 ! -------------------------------------------

                 STIFF(p,q) = STIFF(p,q) + SUM(MATMUL(C, dBasisdx(q,:)) * dBasisdx(p,:))*detJ*IP % s(t)
                 IF(ElectroDynamics) THEN
                   DAMP(p,q) = DAMP(p,q) + P_ip*SUM( dBasisdx(q,:)*dBasisdx(p,:) )*detJ*IP % s(t)
                 END IF
               END DO
               DO j=1,nd-np
                 q = j+np

                 ! Compute the conductivity term <C A,grad v> for 
                 ! mass matrix (anisotropy taken into account)
                 ! -------------------------------------------
                 IF(ElectroDynamics) THEN
                   MASS(p,q) = MASS(p,q) + P_ip*SUM( dBasisdx(p,:)*WBasis(j,:) )*detJ*IP % s(t)
                   DAMP(p,q) = DAMP(p,q) + SUM(MATMUL(C, Wbasis(j,:))*dBasisdx(p,:))*detJ*IP % s(t)
                 ELSE
                   MASS(p,q) = MASS(p,q) + SUM(MATMUL(C, Wbasis(j,:))*dBasisdx(p,:))*detJ*IP % s(t)
                 END IF

                 ! Compute the conductivity term <C grad V, eta> for 
                 ! stiffness matrix (anisotropy taken into account)
                 ! ------------------------------------------------
                 STIFF(q,p) = STIFF(q,p) + SUM(MATMUL(C, dBasisdx(p,:))*WBasis(j,:))*detJ*IP % s(t)
                 IF(ElectroDynamics) THEN
                   DAMP(q,p) = DAMP(q,p) + P_ip*SUM( dBasisdx(p,:)*WBasis(j,:) )*detJ*IP % s(t)
                 END IF
               END DO
             END DO

           ELSE
             ! ---------------------------------------------------------------
             ! This is the steady state branch. 
             ! ------------------------------------------------------------------
             IF (.NOT. LaminateStack ) THEN
               DO p=1,np
                 DO q=1,np

                   ! Compute the conductivity term <C grad V,grad v> for stiffness 
                   ! matrix (anisotropy taken into account)
                   ! -------------------------------------------
                   STIFF(p,q) = STIFF(p,q) + SUM(MATMUL(C, dBasisdx(q,:)) * dBasisdx(p,:))*detJ*IP % s(t)
                 END DO

                 DO j=1,nd-np
                   q = j+np
                   ! The equation for the vector potential:
                   ! Compute the conductivity term <C grad V, eta> for 
                   ! stiffness matrix (anisotropy taken into account)
                   ! ------------------------------------------------
                   STIFF(q,p) = STIFF(q,p) + SUM(MATMUL(C, dBasisdx(p,:))*WBasis(j,:))*detJ*IP % s(t)
                 END DO
               END DO
             END IF
           END IF
         END IF CONDUCTOR
       END IF ! (.NOT. CoilBody)

       LORENTZ_EFFECT: IF ( HasVelocity .AND. .NOT. Transient) THEN
         !
         ! All terms that are added here depend on the electrical conductivity,
         ! so they have an effect on a conductor only.
         !
         A_CONDUCTOR: IF ( SUM(ABS(C)) > AEPS ) THEN
           !
           ! In the case of steady state model add the effect of v x curl A to 
           ! the electromagnetic field: 
           !
           DO p=1,np
             DO j=1,nd-np
               q = j+np
#ifndef __INTEL_COMPILER
               STIFF(p,q) = STIFF(p,q) - &
                   SUM(MATMUL(C,CrossProduct(velo, RotWBasis(j,:)))*dBasisdx(p,:))*detJ*IP % s(t)
#else
               ! Ifort workaround
               RotWJ(1:3) = RotWBasis(j,1:3)
               ! VeloCrossW(1:3) = CrossProduct(velo(1:3), RotWJ(1:3))
               ! CVelo(1:3)=MATMUL(C(1:3,1:3),VeloCrossW(1:3))
               CVelo(1:3) = C(1:3,1)*(velo(2)*RotWJ(3) - velo(3)*RotWJ(2))
               CVelo(1:3) = CVelo(1:3) + C(1:3,2)*(-velo(1)*RotWJ(3) + velo(3)*RotWJ(1))
               CVelo(1:3) = CVelo(1:3) + C(1:3,3)*(velo(1)*RotWJ(2) - velo(2)*RotWJ(1))
               CVeloSum = REAL(0,dp)
               DO k=1,3
                 CVeloSum = CVeloSum + CVelo(k)*dBasisdx(p,k)
               END DO
               STIFF(p,q) = STIFF(p,q) - CVeloSum*detJ*IP % s(t)
#endif
             END DO
           END DO

           DO i = 1,nd-np
             p = i+np
             DO j = 1,nd-np
               q = j+np          
               STIFF(p,q) = STIFF(p,q) - &
                   SUM(WBasis(i,:)*MATMUL(C,CrossProduct(velo, RotWBasis(j,:))))*detJ*IP%s(t)
             END DO
           END DO

         END IF A_CONDUCTOR
       END IF LORENTZ_EFFECT

       !-----------------------------------------------------------------
       ! The equations for the H(curl)-conforming part, i.e. the equation 
       ! for the vector potential
       !    C*dA/dt + curl(nu*curl(A)) + C*grad(V) =  
       !    J^s + curl(M^s) - C*grad(V^s),
       ! with the term C*grad(V) already handled above.
       ! -----------------------------------------------------------------
       DO i = 1,nd-np
         p = i+np
         FORCE(p) = FORCE(p) + (SUM(L*WBasis(i,:)) + &
            SUM(M*RotWBasis(i,:)))*detJ*IP%s(t) 
         DO j = 1,nd-np
           q = j+np

           IF (HasTensorReluctivity .OR. HasReluctivityFunction ) THEN
             STIFF(p,q) = STIFF(p,q) &
                 + SUM(RotWBasis(i,:) * MATMUL(A_t, RotWBasis(j,:)))*detJ*IP%s(t)
           ELSE
             STIFF(p,q) = STIFF(p,q) + mu * SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t)              
           END IF
           IF ( Newton ) THEN
             IF ( HasHBCurve ) THEN
               JAC(p,q) = JAC(p,q) + muder * SUM(B_ip(:)*RotWBasis(j,:)) * &
                   SUM(B_ip(:)*RotWBasis(i,:))*detJ*IP % s(t)/Babs
             ELSE IF( HasReluctivityFunction ) THEN
               ! Note: check this though for anisotropic cases!!!
               DO k=1,3
                 JAC(p,q) = JAC(p,q) + SUM(A_t_der(k,:) *  B_ip(k) * RotWBasis(j,:)) * &
                     SUM(B_ip(:)*RotWBasis(i,:))*detJ*IP % s(t)/Babs
               END DO
             END IF
           END IF

           ! Compute the conductivity term <C A,eta> for 
           ! mass matrix (anisotropy taken into account)
           ! This is not used in the case of stranded coil:
           ! ----------------------------------------------
           IF(ElectroDynamics) THEN
             MASS(p,q) = MASS(p,q) + P_ip*SUM(WBasis(j,:)*WBasis(i,:) )*detJ*IP % s(t)
             IF (CoilType /= 'stranded') THEN
               DAMP(p,q) = DAMP(p,q) + SUM(MATMUL(C, WBasis(j,:))*WBasis(i,:) )*detJ*IP % s(t)
             END IF
           ELSE
             IF (CoilType /= 'stranded') THEN
               MASS(p,q) = MASS(p,q) + SUM(MATMUL(C, WBasis(j,:))*WBasis(i,:) )*detJ*IP % s(t)
             END IF
           END IF

           ! Compute the low frequency eddy term for laminate stack model.
           ! Note that the conductivity term <C A, eta> above can be used to 
           ! introduce the anisotropic effect in the laminate stack. However, 
           ! in classical approach of the Low-Frequency model it is set 
           ! to zero (this is left to the user to decide).
           ! -------------------------------------------------------------------
           IF (LaminateStackModel=='low-frequency model') THEN
               MASS(p,q) = MASS(p,q) + LocalLamCond * LocalLamthick**2/12d0 * & 
                           SUM(RotWBasis(i,:)*RotWBasis(j,:)) * detJ*IP % s(t)
           END IF

         END DO
       END DO

       ! In steady state we can utilize the scalar variable
       ! for gauging the vector potential.
       IF ( SteadyGauge .OR. LocalGauge ) THEN
         DO j = 1, np
           q = j
           DO i = 1,nd-np
             p = i+np
             STIFF(q,p) = STIFF(q,p) - gauge_coeff * SUM(dBasisdx(j,:)*WBasis(i,:))*detJ*IP % s(t)
             STIFF(p,q) = STIFF(p,q) - gauge_coeff * SUM(dBasisdx(j,:)*WBasis(i,:))*detJ*IP % s(t)
           END DO
         END DO

         IF ( HasStabC ) THEN
           DO q = 1, np
             DO p = 1, np
               STIFF(p,q) = STIFF(p,q) + gauge_penalize_c*SUM(dBasisdx(q,:)*dBasisdx(p,:))*detJ*IP % s(t) &
                     + gauge_penalize_m*Basis(q)*Basis(p)*detJ*IP%s(t)
             END DO
           END DO
         END IF

       END IF
       
       ! Add mass-type term to the stiffness matrix for regularization if requested
       ! This turns the steady state system to "curl nu curl u + epsilon u = J"
       ! For discussion on regularization see, e.g., "Cassagrande, Hiptmair, Ostrowski, 
       ! An a priori error estimate for interior penalty discretizations of the Curl-Curl 
       ! operator on non-conforming meshes" section 6.
       
       IF ( RegularizeWithMass ) THEN
        DO j = 1, nd-np
          q = j + np
          DO i = 1, nd-np
            p = i + np
            STIFF(p,q) = STIFF(p,q) + mass_reg_epsilon*SUM(WBasis(i,:)*WBasis(j,:))*detJ*IP%s(t)
          END DO
        END DO
       END IF

    END DO

    IF ( Newton ) THEN
      IF( HasHBCurve .OR. HasReluctivityFunction ) THEN
        STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + JAC
        FORCE(1:nd) = FORCE(1:nd) + MATMUL(JAC,Aloc)
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE ConstraintMatrix( MASS, DAMP, STIFF, FORCE, &
         Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:), DAMP(:,:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), Permittivity(nd), P_ip
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, k, p, q, np, EdgeBasisDegree
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp
    DAMP  = 0.0_dp

    CALL GetPermittivity(GetMaterial(), Permittivity, n)

    !Numerical integration:
    !----------------------
    IP = GaussPointsAdapt(Element, Solver, EdgeBasis=.TRUE. )

    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )

       P_ip = SUM( Permittivity(1:n) * Basis(1:n) )
       DO p=1,np
         DO q=1,np
           STIFF(p,q) = STIFF(p,q) + P_ip*SUM( dBasisdx(q,:)*dBasisdx(p,:) )*detJ*IP % s(t)
         END DO
         DO j=1,nd-np
           q = j+np
           DAMP(p,q) = DAMP(p,q) + P_ip*SUM( Wbasis(j,:)*dBasisdx(p,:) )*detJ*IP % s(t)
         END DO
       END DO
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ConstraintMatrix
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE LocalFixMatrix( FORCE, Element, n, nd, PiolaVersion, SecondOrder )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3)
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ, L(3), JFixPot(nd)
    LOGICAL :: Stat 
    INTEGER :: t, i, p, np, EdgeBasisDegree
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes

!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
      EdgeBasisDegree = 2
    ELSE
      EdgeBasisDegree = 1
    END IF
    
    CALL GetElementNodes( Nodes )

    FORCE = 0.0_dp
    !CALL GetScalarLocalSolution( JFixPot, 'JFix')

    JFixPot(1:n) = JFixVar % Values( JFixVar % Perm( Element % NodeIndexes ) )
    
    IF( SUM(ABS(JFixPot(1:n))) < TINY(DetJ) ) RETURN
    
    ! Numerical integration:
    !----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree )

    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    DO t=1,IP % n

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )            

      L = MATMUL(-JFixPot(1:n), dBasisdx(1:n,:))
      DO i = 1,nd-np
        p = i+np
        FORCE(p) = FORCE(p) + SUM(L*WBasis(i,:)) * detJ * IP%s(t) 
      END DO
    END DO

    CALL DefaultUpdateForce(FORCE, Element )
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalFixMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  MASS, STIFF, FORCE, LOAD, Bcoef, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: LOAD(:,:), Bcoef(:)
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,L(3), t1(3), t2(3)
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), B, F, TC, c
    LOGICAL :: Stat, Found, PEC, OutFlow
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, k, ii,jj, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes, Element )

    PEC = GetLogical( GetBC(), 'Perfect Electric Conductor', Found)
    Outflow = GetLogical( GetBC(), 'Outflow BC', Found)
    IF(OutFlow) THEN
      c = GetCReal( GetBC(), 'Wave speed' )
    END IF

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)
    
    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))

    DO t=1,IP % n

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )
      
       B  = SUM(Basis(1:n) * Bcoef(1:n))
       L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))

       F  = SUM(LOAD(4,1:n)*Basis(1:n))
       TC = SUM(LOAD(5,1:n)*Basis(1:n))

       IF(PEC) THEN
         DO i = 1,nd-np
           p = i+np
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + 1d5*SUM(WBasis(i,:)*dBasisdx(q,:))*detJ*IP % s(t)
           END DO

           DO j = 1,nd-np
             q = j+np
             MASS(p,q) = MASS(p,q) + 1d5*SUM(WBasis(i,:)*Wbasis(j,:))*detJ*IP % s(t)
           END DO
         END DO
       END IF

       IF(Outflow) THEN
         DO i = 1,nd-np
           p = i+np
           DO j = 1,nd-np
             q = j+np
             MASS(p,q) = MASS(p,q) + SUM(WBasis(i,:)*Wbasis(j,:))*detJ*IP % s(t) / c
           END DO
         END DO
       END IF

       ! Compute element stiffness matrix and force vector:
       !---------------------------------------------------
       DO p=1,np
         FORCE(p) = FORCE(p) + F*Basis(p)*detJ*IP % s(t)
         DO q=1,np
           STIFF(p,q) = STIFF(p,q) + TC * Basis(p)*Basis(q)*detJ*IP % s(T)
         END DO
       END DO

       DO i = 1,nd-np
         p = i+np
         FORCE(p) = FORCE(p) - SUM(L*WBasis(i,:))*detJ*IP % s(t)
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + B * SUM(WBasis(i,:)*Wbasis(j,:))*detJ*IP % s(t)
         END DO
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  FUNCTION LocalFluxBC( LOAD, Element, n, nd ) RESULT(Bn)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: LOAD(:,:), Bn
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Edge, Parent
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,L(3),Ln
    REAL(KIND=dp) :: Normal(3)
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes,  Element )
    !
    ! Integrate (B,n) over boundary face:
    ! -----------------------------------
    IP = GaussPoints(Element)
    Bn = 0._dp
    DO t=1,IP % n
      stat = ElementInfo( Element,Nodes,IP % U(t),IP % V(t), &
                 IP % W(t),detJ,Basis,dBasisdx )

      Normal=NormalVector(Element,Nodes,IP % u(t),ip % v(t),.TRUE.)
      Ln = SUM(LOAD(4,1:n)*Basis(1:n))
      L  = MATMUL(LOAD(1:3,1:n), Basis(1:n))
      Bn = Bn + Detj * IP % S(t) * (Ln+SUM(L*Normal))
    END DO
!------------------------------------------------------------------------------
  END FUNCTION LocalFluxBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixAirGapBC(  STIFF, FORCE, LOAD, GapLength, AirGapMu, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: LOAD(:,:), GapLength(:), AirGapMu(:)
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3), localGapLength, muAir, muVacuum
    LOGICAL :: Stat
    INTEGER, POINTER :: EdgeMap(:,:)
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    EdgeBasisDegree = 1
    IF (SecondOrder) EdgeBasisDegree = 2

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    MASS  = 0.0_dp

    muVacuum = 4 * PI * 1d-7

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = RotWBasis, USolver = pSolver )

       localGapLength  = SUM(Basis(1:n) * GapLength(1:n))
       muAir  = SUM(Basis(1:n) * AirGapMu(1:n))
 
       DO i = 1,nd-np
         p = i+np
         DO j = 1,nd-np
           q = j+np
           STIFF(p,q) = STIFF(p,q) + localGapLength / (muAir*muVacuum) * &
              SUM(RotWBasis(i,:)*RotWBasis(j,:))*detJ*IP%s(t)
         END DO
       END DO  
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixAirGapBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixThinLine( MASS,STIFF, FORCE, LOAD, CrossectArea, Conductivity, Element, &
      n, nd, SecondOrder)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: LOAD(:,:), CrossectArea(:), Conductivity(:)
    REAL(KIND=dp) :: MASS(:,:),STIFF(:,:), FORCE(:)
    INTEGER :: n, nd
    LOGICAL :: SecondOrder
    TYPE(Element_t), POINTER :: Element, Parent, Edge
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3),Basis(n),dBasisdx(n,3),DetJ
    REAL(KIND=dp) :: C, Area
    LOGICAL :: Stat
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: t, i, j, np, p, q, EdgeBasisDegree

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, Element )

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2
    ELSE
      EdgeBasisDegree = 1
    END IF

    MASS  = 0.0_dp
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE., &
        EdgeBasisDegree=EdgeBasisDegree)

    np = n*MAXVAL(Solver % Def_Dofs(GetElementFamily(Element),:,1))
    DO t=1,IP % n

      stat = EdgeElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
          dBasisdx = dBasisdx, BasisDegree = EdgeBasisDegree, &
          ApplyPiolaTransform = .TRUE.)

       C = SUM(Basis(1:n) * Conductivity(1:n))
       Area = SUM(Basis(1:n) * CrossectArea(1:n))

       CONDUCTOR: IF ( ABS(C) > AEPS ) THEN
         DO p=1,np
           DO q=1,np

             ! Compute the conductivity term <C grad V x n,grad v x n> for stiffness 
             ! matrix (without anisotropy taken into account)
             ! -------------------------------------------

             STIFF(p,q) = STIFF(p,q) + Area * C * SUM(dBasisdx(q,:) * dBasisdx(p,:))*detJ*IP % s(t)

           END DO
           DO j=1,nd-np
             q = j+np

             ! Compute the conductivity term <C A x n,grad v x n> for 
             ! mass matrix (without anisotropy taken into account)
             ! -------------------------------------------
             MASS(p,q) = MASS(p,q) + Area * C * SUM(WBasis(j,:)*dBasisdx(p,:))*detJ*IP % s(t)

             ! Compute the conductivity term <C grad V x n, eta x n> for 
             ! stiffness matrix (without anisotropy taken into account)
             ! ------------------------------------------------
             STIFF(q,p) = STIFF(q,p) + Area * C * SUM(dBasisdx(p,:)*WBasis(j,:))*detJ*IP % s(t)
           END DO
         END DO

         DO i=1,nd-np
           p = i+np
           DO j=1,nd-np
             q = j+np

             ! Compute the conductivity term <C A x n, eta x n> for 
             ! mass matrix (without anisotropy taken into account)
             ! -------------------------------------------
             MASS(p,q) = MASS(p,q) + Area * C * SUM(WBasis(i,:)*Wbasis(j,:))*detJ*IP % s(t)
           END DO
         END DO

       END IF CONDUCTOR
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixThinLine
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
  SUBROUTINE DirichletAfromB()
!------------------------------------------------------------------------------
    USE ElementDescription, ONLY: GetEdgeMap

    IMPLICIT NONE
    REAL(KIND=dp) :: s,p(3),q(3),cx(3),r,xmin,ymin,zmin,xmax,ymax,zmax
    TYPE(ListMatrixEntry_t), POINTER :: Ltmp
    TYPE(Matrix_t), POINTER :: Smat
    TYPE(Nodes_t),SAVE :: Nodes
    TYPE(ValueList_t), POINTER :: BC

    LOGICAL :: Found, Found1,Found2,Found3,L1,L2,L3
    INTEGER :: i,j,k,l,m,t,ii,Faces,n,nd,Active,je1,je2,pe1,pe2

    TYPE(Element_t), POINTER :: Element, Edge, Edge1
    REAL(KIND=dp), ALLOCATABLE :: Bn(:)
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER, ALLOCATABLE :: dMap(:),FaceMap(:)
    LOGICAL, ALLOCATABLE :: FluxBoundaryEdge(:), CycleEdges(:), UsedFaces(:)
!------------------------------------------------------------------------------
    ALLOCATE(FluxBoundaryEdge(Mesh % NumberOFEdges)); FluxBoundaryEdge=.FALSE.

    Active = GetNOFBoundaryElements()

    DO t=1,Active
       Element => GetBoundaryElement(t)

       IF ( GetElementFamily()==1 ) CYCLE
       BC=>GetBC()
       IF (.NOT. ASSOCIATED(BC) ) CYCLE

       Found = ListCheckPrefix(BC,'Magnetic Flux Density')

       IF ( Found ) THEN
         SELECT CASE(GetElementFamily())
         CASE(2)
           CYCLE !what would it mean in 2D,at least with only B_z solved?
         CASE(3,4)
           k = GetBoundaryFaceIndex(Element); Element => Mesh % Faces(k)
         END SELECT
         IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
         FluxBoundaryEdge(Element % EdgeIndexes)=.TRUE.
       END IF
    END DO

    FluxCount = COUNT(FluxBoundaryEdge)
    IF ( FluxCount==0 ) THEN
      DEALLOCATE(FluxBoundaryEdge); RETURN
    END IF
    
    IF (.NOT.ALLOCATED(FluxMap) ) ALLOCATE(FluxMap(FluxCount))
    FluxCount = 0
    FluxMap   = 0
    DO i=1,Mesh % NumberOfEdges
      IF ( FluxBoundaryEdge(i) ) THEN
        FluxCount = FluxCount+1
        FluxMap(FluxCount) = i
      END IF
    END DO
    DEALLOCATE(FluxBoundaryEdge)
    
    DO i=1,FluxCount
      Edge => Mesh % Edges(FluxMap(i))
      Edge % BoundaryInfo % Left => NULL()
      Edge % BoundaryInfo % Right => NULL()
    END DO

    ALLOCATE(FaceMap(Mesh % NumberOfFaces)); FaceMap=0
    Faces = 0
   DO t=1,Active
      Element => GetBoundaryElement(t)

      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      Found = ListCheckPrefix(BC,'Magnetic Flux Density')
      IF ( .NOT. Found ) CYCLE

      k = GetBoundaryFaceIndex(Element); Element=>Mesh % Faces(k)
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
      Faces = Faces+1
      FaceMap(k) = Faces

      DO i=1,Element % TYPE % NumberOfNodes
        Edge => Mesh % Edges(Element % EdgeIndexes(i))
        IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Left)) THEN
           Edge % BoundaryInfo % Left => Element
        ELSE IF (.NOT.ASSOCIATED(Edge % BoundaryInfo % Right)) THEN
           Edge % BoundaryInfo % Right => Element
        END IF
      END DO
    END DO

    ! Make gauge tree for the boundary:
    ! ---------------------------------
    CALL GaugeTreeFluxBC(Solver,Mesh,TreeEdges,BasicCycles,FluxCount,FluxMap)

    CALL Info(Caller, 'Boundary tree edges: '//i2s(COUNT(TreeEdges(FluxMap))) // &
        ' of total: '//I2S(FluxCount),Level=5)
    
    ! Get (B,n) for BC faces:
    ! -----------------------
    ALLOCATE(Bn(Faces))
    DO t=1,Active
      Element => GetBoundaryElement(t)

      IF ( GetElementFamily()==1 ) CYCLE
      BC=>GetBC()
      IF (.NOT. ASSOCIATED(BC) ) CYCLE

      n  = GetElementNOFNodes(Element)
      CALL GetRealVector(BC,Load(1:3,1:n),'Magnetic Flux Density',Found1)
      LOAD(4,1:n) = GetReal(BC,'Magnetic Flux Density {n}',Found)

      IF (Found.OR.Found1) THEN
        k = GetBoundaryFaceIndex(Element)
        Element => Mesh % Faces(k)
        IF (.NOT.ActiveBoundaryElement(Element)) CYCLE        
        nd = GetElementNOFDOFs(Element)
        Bn(FaceMap(k))=LocalFluxBC(LOAD,Element,n,nd)
      END IF
    END DO

    !
    ! Calculate value for free edges using the fundamental loop basis
    ! generated by GaugeTreeFluxBC():
    ! ---------------------------------------------------------------
    ALLOCATE(CycleEdges(Mesh % NumberOFEdges), UsedFaces(Faces))
    CycleEdges = .FALSE.
    
    i = MAXVAL(BasicCycles(1:FluxCount) % Degree)
    ALLOCATE(dMap(i))

    Smat => GetMatrix()
    DO i=1,SIZE(BasicCycles)
      IF (BasicCycles(i) % Degree<=0 ) CYCLE

      ! 
      ! Extract loop edge indices: 
      ! --------------------------
      j = 0
      Ltmp => BasicCycles(i) % Head
      DO WHILE(ASSOCIATED(Ltmp))
        j = j + 1
        dMap(j) = Ltmp % Index; Ltmp => Ltmp % Next
      END DO
      IF ( j<= 0 ) CYCLE

      !
      ! Orient edges to form a polygonal path:
      ! --------------------------------------
      Edge  => Mesh % Edges(dMap(j))
      Edge1 => Mesh % Edges(dMap(j-1))
      IF ( ANY(Edge % NodeIndexes(1)==Edge1 % NodeIndexes) ) THEN
        l = Edge % NodeIndexes(1)
        Edge % NodeIndexes(1) = Edge % NodeIndexes(2)
        Edge % NodeIndexes(2) = l
      END IF
 
      DO k=j-1,1,-1
        Edge1 => Mesh % Edges(dMap(k))
        IF (Edge % NodeIndexes(2)==Edge1 % NodeIndexes(2)) THEN
          l = Edge1 % NodeIndexes(1)
          Edge1 % NodeIndexes(1) = Edge1 % NodeIndexes(2)
          Edge1 % NodeIndexes(2) = l
        END IF
        Edge => Edge1
      END DO

      !
      ! Try to find which way is inside...
      ! ----------------------------------
      Edge => Mesh % Edges(dMap(j))
      Element => Edge % BoundaryInfo % Left

      IF ( j==3 ) THEN
        m = 0
        DO k=1,3
          DO l=1,3
            IF (dMap(l)==Element % EdgeIndexes(k)) m=m+1
          END DO
        END DO
        L1 = m==3
        IF ( .NOT. L1 ) Element=>Edge % BoundaryInfo % Right
        S = Bn(FaceMap(Element % ElementIndex))
      ELSE
        ! If not a triangle, try a (planar) polygonal test. This
        ! will fail for general 3D paths. We'll spot the failure
        ! later by trial and error...Might be preferable to skip
        ! this altogether? Dunno....
        ! ------------------------------------------------------
        xmin=HUGE(xmin); xmax=-HUGE(xmax);
        ymin=HUGE(ymin); ymax=-HUGE(ymax);
        zmin=HUGE(zmin); zmax=-HUGE(zmax);
        DO k=1,j
          Edge1 => Mesh % Edges(dMap(k))
          DO l=1,2
            m = Edge1 % NodeIndexes(l)
            xmin = MIN(xmin,Mesh % Nodes % x(m))
            ymin = MIN(ymin,Mesh % Nodes % y(m))
            zmin = MIN(zmin,Mesh % Nodes % z(m))

            xmax = MAX(xmax,Mesh % Nodes % x(m))
            ymax = MAX(ymax,Mesh % Nodes % y(m))
            zmax = MAX(zmax,Mesh % Nodes % z(m))
          END DO
        END DO
        L1 = xmax-xmin > ymax-ymin
        L2 = xmax-xmin > zmax-zmin
        L3 = ymax-ymin > zmax-zmin
        IF ( l1 ) THEN
          l=1
          IF ( l3 ) THEN
            m=2; n=3
          ELSE
            m=3; n=2
          END IF
        ELSE
          IF ( l2 ) THEN
            l=1; m=2; n=3
          ELSE
            l=3; m=1; n=2
          END IF
        END IF
        cx(l) = SUM(Mesh % Nodes % x(Element % NodeIndexes))/3._dp
        cx(m) = SUM(Mesh % Nodes % y(Element % NodeIndexes))/3._dp
        cx(n) = SUM(Mesh % Nodes % z(Element % NodeIndexes))/3._dp

        L1 = .FALSE.
        DO k=j,1,-1
          Edge1 => Mesh % Edges(dMap(k))
          je1 = Edge1 % NodeIndexes(1)
          je2 = Edge1 % NodeIndexes(2)
          p(l) = Mesh % Nodes % x(je1)
          p(m) = Mesh % Nodes % y(je1)
          p(n) = Mesh % Nodes % z(je1)

          q(l) = Mesh % Nodes % x(je2)
          q(m) = Mesh % Nodes % y(je2)
          q(n) = Mesh % Nodes % z(je2)

          IF ((q(2)>cx(2)).NEQV.(p(2)>cx(2))) THEN
            IF (cx(1)<(p(1)-q(1))*(cx(2)-q(2))/(p(2)-q(2))+q(1)) L1=.NOT.L1
          END IF
        END DO
        IF (.NOT.L1) THEN
          IF (ASSOCIATED(Edge % BoundaryInfo % Right)) &
            Element=>Edge % BoundaryInfo % Right
        END IF

        ! Compute integral of (B,n) inside the cycle path
        ! -----------------------------------------------
        CycleEdges(dMap(1:j))=.TRUE.
        DO m=1,2
          S=0; UsedFaces = .FALSE.;
          IF( FloodFill(Element,CycleEdges,FaceMap,UsedFaces,Bn,S,0) )EXIT

          ! the in/out guess was wrong, try the other way:
          ! ----------------------------------------------
          IF (ASSOCIATED(Edge % BoundaryInfo % Right,Element)) THEN
            Element => Edge % BoundaryInfo % Left
          ELSE
            Element => Edge % BoundaryInfo % Right
          END IF

          IF(.NOT.ASSOCIATED(Element)) CALL Fatal('DirichletAfromB', 'Floodfill failing.')
        END DO

        CycleEdges(dMap(1:j))=.FALSE.
      END IF

      !
      ! Orient edge to parent triangle...
      ! ---------------------------------
      je1 = Edge % NodeIndexes(1)
      je2 = Edge % NodeIndexes(2)
      EdgeMap => GetEdgeMap(GetElementFamily(Element))
      DO t=1,Element % TYPE % NumberOfEdges
        pe1 = Element % NodeIndexes(EdgeMap(t,1))
        pe2 = Element % NodeIndexes(EdgeMap(t,2))
        IF (pe1==je1.AND.pe2==je2 .OR. pe1==je2.AND.pe2==je1) EXIT
      END DO
      IF ( pe1/=je1 ) S=-S

      !
      ! ...because we now know how to orient against outward normal:
      ! ------------------------------------------------------------
      CALL GetElementNodes(Nodes,Element)
      p = NormalVector(Element,Nodes,0._dp,0._dp)
      q = NormalVector(Element,Nodes,0._dp,0._dp,.TRUE.)
      IF ( SUM(p*q)<0 ) S=-S

      !
      ! Check whether some edges in the path have nonzero values,
      ! if so, substract from integral:
      ! ---------------------------------------------------------
      DO k=j-1,1,-1
        l = Perm(dMap(k)+nNodes)
        R = Smat % RHS(l)/Smat % Values(Smat % Diag(l))
        IF ( R==0 ) CYCLE

        Edge1 => Mesh % Edges(dMap(k))
        pe1=Edge1 % NodeIndexes(1)
        pe2=Edge1 % NodeIndexes(2)
        IF ( pe2<pe1 ) R=-R; S=S-R
      END DO

      !
      ! ...and finally we should have the edge value:
      ! ---------------------------------------------
      IF ( je2<je1 ) S=-S
      CALL SetDOFtoValue(Solver,dMap(j),S)
    END DO
    DEALLOCATE(dMap, CycleEdges, FaceMap, UsedFaces, Bn)
    !CALL List_FreeMatrix(SIZE(BasicCycles), BasicCycles)
!------------------------------------------------------------------------------
  END SUBROUTINE DirichletAfromB 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE FUNCTION FloodFill(Element,CycleEdges, &
          FaceMap,UsedFaces,Bn,CycleSum, level) RESULT(Found)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: e, Element
    REAL(KIND=dp) :: CycleSum, Bn(:)
    INTEGER :: i,j,n, FaceMap(:), level
    LOGICAL :: CycleEdges(:), UsedFaces(:), Found, L

    Found=.FALSE.
    IF (.NOT.ASSOCIATED(Element)) RETURN

    n=FaceMap(Element % ElementIndex)
    IF (UsedFaces(n)) THEN
      Found=.TRUE.; RETURN
    END IF
    UsedFaces(n)=.TRUE.
    CycleSum = CycleSum+Bn(n)

    DO i=1,Element % TYPE % NumberOfEdges
      j = Element % EdgeIndexes(i)

      IF ( CycleEdges(j) ) CYCLE

      e => Mesh % Edges(j) % BoundaryInfo % Right
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum,level+1)) RETURN
!     L=FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum,level+1)

      e => Mesh % Edges(j) % BoundaryInfo % Left
      IF(.NOT.FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum,level+1)) RETURN
!     L=FloodFill(e,CycleEdges,FaceMap,UsedFaces,Bn,CycleSum,level+1)
    END DO
    Found=.TRUE.; RETURN
!------------------------------------------------------------------------------
  END FUNCTION FloodFill
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  SUBROUTINE AddLocalBNorm( Element, n, nd, PiolaVersion, SecondOrder, &
      BInteg, vinteg, BMin, BMax )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: PiolaVersion, SecondOrder
    REAL(KIND=dp) :: BInteg, vinteg, BMin, BMax
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Aloc(nd), B_ip(3), B2
    REAL(KIND=dp) :: WBasis(nd,3), RotWBasis(nd,3),Basis(n),dBasisdx(n,3),DetJ
    LOGICAL :: Stat
    INTEGER :: t, i, j, np, EdgeBasisDegree
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (SecondOrder) THEN
       EdgeBasisDegree = 2
    ELSE
       EdgeBasisDegree = 1
    END IF

    CALL GetElementNodes( Nodes )
    CALL GetScalarLocalSolution(Aloc)
        
    ! Numerical integration:
    !------------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=EdgeBasisDegree )

    np = n*Solver % Def_Dofs(GetElementFamily(Element),Element % BodyId,1)
    DO t=1,IP % n            
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, &
          EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = pSolver ) 

      B_ip = MATMUL( Aloc(np+1:nd), RotWBasis(1:nd-np,:) )
      b2 = SUM(B_ip**2)
      
      BMin = MIN( BMin, b2 )
      BMax = MAX( BMax, b2 ) 
      BInteg = BInteg + b2 * detJ * IP % s(t)
      vinteg = vinteg + detJ * IP % s(t)
    END DO
  END SUBROUTINE AddLocalBNorm
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 END SUBROUTINE WhitneyAVSolver
!------------------------------------------------------------------------------




!/*****************************************************************************/
! *
! *  Utilities written as solvers to compute the Helmholtz projection P(A)
! *  of a curl-conforming vector field A. The projection can be obtained as 
! *  P(A) = A - W where  W is the curl-conforming field fitted to represent 
! *  grad Phi, with Phi being a H1-regular scalar field.
! * 
! *  This file contains the time-domain version of the transformation and also applies the
! *  correction to the V field within conducting regions.
! *
! *
! *  Authors: Mika Malinen, Juha Ruokolainen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 20, 2020
! *  Last Modified: June 18, 2021, Juha
! *
!******************************************************************************


!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjectorT_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i
  TYPE(ValueList_t), POINTER :: SolverParams, SParams
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  DO i=1,Model % NumberOfSolvers
    IF(ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO

  IF (i<=Model % NumberOfSolvers) THEN
    SParams => Model % Solvers(i) % Values
    CALL ListAddNewInteger( SolverParams,'Mortar BC Master Solver',i)
#if 0
    IF( GetLogical( SParams, 'Apply Mortar BCs', Found) ) THEN
      CALL ListAddLogical( SolverParams, 'Apply Mortar BCs', .TRUE. )
    END IF
    IF( GetLogical( SParams, 'Apply Conforming BCs', Found) ) THEN
      CALL ListAddLogical( SolverParams, 'Apply Conforming BCs', .TRUE. )
    END IF
    IF( GetLogical( SParams, 'Mortar BCs Additive', Found) ) THEN
      CALL ListAddLogical( SolverParams, 'Mortar BCs Additive', .TRUE. )
    END IF
    CALL ListAddLogical( SolverParams, 'Projector Skip Edges', .TRUE. )
#endif

    CALL ListCopyPrefixedKeywords(Model % Solvers(i) % Values, SolverParams, 'HelmholtzProjector:')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjectorT_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjectorT_Init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found
  INTEGER :: i,j
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
! CALL ListAddNewLogical(SolverParams, 'Linear System Refactorize', .FALSE.)

  CALL ListAddString( SolverParams, 'Variable', 'P' )
! CALL ListAddLogical( SolverParams, 'Variable Output',.FALSE. )

  DO i=1,Model % NumberOfSolvers
    IF(ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO
  CALL ListAddString( SolverParams, 'Potential Variable', GetVarName(Model % Solvers(i) % Variable))

  CALL ListAddLogical( SolverParams, 'Linear System Symmetric', .TRUE. )
  CALL ListAddString(  SolverParams, 'Linear System Solver', 'Iterative' )
  CALL ListAddString(  SolverParams, 'Linear System Preconditioning', 'ILU' )
  CALL ListAddInteger( SolverParams, 'Linear System Residual Output', 25 )
  CALL ListAddInteger( SolverParams, 'Linear System Max Iterations', 2000 )
  CALL ListAddString(  SolverParams, 'Linear System Iterative Method', 'CG' )
  CALL ListAddConstReal(   SolverParams, 'Linear System Convergence Tolerance', 1.0d-9 )

  DO j=1,Model % NumberOfBCs
    IF ( ListCheckPrefix( Model % BCs(j) % Values, &
               TRIM(GetVarName(Model % Solvers(i) % Variable)) // ' {e}' ) ) THEN
      CALL ListAddConstReal( Model % BCs(j) % Values, 'P', 0._dp )
    END IF
  END DO
!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjectorT_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute a H1-regular scalar field to obtain the Helmholtz projection P(A)
!> of a curl-conforming vector field A. Given the solution field Phi of this 
!> solver, the projection can be evaluated as P(A) = A - grad Phi.
!------------------------------------------------------------------------------
SUBROUTINE HelmholtzProjectorT(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix

  INTEGER :: i, j,k,l,n, n_pot, nd_pot, t
  INTEGER :: dim, PotDOFs
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:), PotSol(:)
  REAL(KIND=dp) :: Norm, Omega
  REAL(KIND=dp), POINTER :: SaveRHS(:), SOL(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName

  SAVE Stiff, Force, PotSol, AllocationsDone

!------------------------------------------------------------------------------
  CALL Info( 'HelmholtzProjector', '--------------------------------------------------',Level=12 )
  CALL Info( 'HelmholtzProjector', 'Computing fixing potential for the vector potential',Level=12 )

  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()
  
  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs

    ALLOCATE( FORCE(n), STIFF(n,n), PotSol(n), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'HelmholtzProjector', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (Found ) THEN
    CALL Info('HelmholtzProjector', 'Solver inherits potential '&
         //TRIM(PotName)//' from solver: '//I2S(i),Level=7)
  ELSE
    CALL Fatal('HelmholtzProjector', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  !
  ! Find some parameters to inherit the vector FE basis as defined in 
  ! the primary solver:
  !
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )
  IF (PiolaVersion) CALL Info('HelmholtzProjector', &
      'Using Piola-transformed finite elements', Level=5)
  
  !-----------------------
  ! System assembly:
  !----------------------
  CALL DefaultInitialize(Solver)

  active = GetNOFActive()
  DO t=1,active
    Element => GetActiveElement(t)
    !
    ! This solver relies on getting basis functions by calling a routine
    ! which returns a curl-conforming basis. It is thus assumed that
    ! the background mesh defines the number of Lagrange basis functions.
    !
    n = GetElementNOFNodes()
   
    ! The DOF counts for the potential (target) variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    CALL GetLocalSolution(PotSol, PotName, USolver=SolverPtr)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
        SecondOrder, n_pot, nd_pot, PotSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(STIFF, FORCE)
  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()

  !
  ! Finally, redefine the potential variable:
  ! -----------------------------------------
  IF( TransientSimulation ) THEN
    DO i=1,Solver % Mesh % NumberOfNodes
      IF( ASSOCIATED( Solver % Mesh % PeriodicPerm ) ) THEN
        IF( Solver % Mesh % PeriodicPerm(i) > 0 ) CYCLE
      END IF

      j = Solver % Variable % Perm(i)
      IF(j==0) CYCLE

      k = SolverPtr % Variable % Perm(i)        
      IF (k == 0) THEN
        CALL Fatal('HelmholtzProjector', &
          'The variable and potential permutations are nonmatching?')
      END IF

      SolverPtr % Variable % Values(k) = SolverPtr % Variable % Values(k) + &
        (Solver % Variable % Values(j) - Solver % Variable % PrevValues(j,1))/ dt
    END DO
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, dim, PiolaVersion, &
               SecondOrder, n_pot, nd_pot, PotSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: dim
    LOGICAL :: PiolaVersion, SecondOrder
    INTEGER :: n_pot, nd_pot      ! The size parameters of target field
    REAL(KIND=dp) :: PotSol(:)  ! The values of target field DOFS
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(3)
    REAL(KIND=dp) :: s, DetJ
    REAL(KIND=dp) :: WBasis(nd_pot-n_pot,3), CurlWBasis(nd_pot-n_pot,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2  
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
          EdgeBasisDegree=EdgeBasisDegree)
    ELSE
      EdgeBasisDegree = 1
      IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)
    END IF
          
    IF( dim == 2 .AND. .NOT. PiolaVersion ) THEN
      CALL Fatal('HelmholtzProjector', 'Use Piola Transform = True needed in 2D')
    END IF
    
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = CurlWBasis, USolver = SolverPtr ) 
      s = detJ * IP % s(t)

      A = MATMUL(PotSol(n_pot+1:nd_pot), WBasis(1:nd_pot-n_pot,:))

      DO p=1,n
        DO q=1,n
          STIFF(p,q) = STIFF(p,q) + SUM(dBasisdx(q,1:dim) * dBasisdx(p,1:dim)) * s
        END DO
      END DO

      DO p=1,n
        FORCE(p) = FORCE(p) + SUM(A * dBasisdx(p,:)) * s
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzProjectorT
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponentT_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: i,j
  CHARACTER(LEN=MAX_NAME_LEN) :: Avname
  LOGICAL :: Found, PiolaVersion, SecondOrder, SecondKind
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
! CALL ListAddLogical(SolverParams, 'Linear System Refactorize', .FALSE.)

! Kernel Variable = String "P"
! Potential Variable = String "AV"
! Linear System Symmetric = True
! Linear System Solver = "Iterative"
! Linear System Preconditioning = None
! Linear System Residual Output = 25
! Linear System Max Iterations = 2000
! Linear System Iterative Method = CG
! Linear System Convergence Tolerance = 1.0e-9

  CALL ListAddString( SolverParams, 'Variable', 'avm' )
  CALL ListAddLogical( SolverParams, 'Variable Output',.FALSE. )
  
  DO i=1,Model % NumberOfSolvers
    IF(ListGetLogical( Model % Solvers(i) % Values, 'Helmholtz Projection', Found)) EXIT
  END DO
  IF(i<=Model % NumberOfSolvers ) THEN
    CALL ListAddNewInteger( SolverParams,'Mortar BC Master Solver',i)
  END IF
    
  AVname = ListGetString( Model % Solvers(i) % Values, 'Variable' )
  
  j = index(AVname, '[')
  IF(j>0) AVname = AVname(1:j-1)
  CALL ListAddString( SolverParams, 'Potential Variable', AVName )

  IF (.NOT. ListCheckPresent(SolverParams, "Element")) THEN
    CALL EdgeElementStyle(Model % Solvers(i) % Values, PiolaVersion, SecondKind, SecondOrder )
    
    IF(PiolaVersion) CALL ListAddLogical(SolverParams, 'Use Piola Transform', PiolaVersion )
    IF(SecondOrder) CALL ListAddLogical(SolverParams, 'Quadratic Approximation', SecondOrder )
    IF(SecondKind) CALL ListAddLogical(SolverParams, 'Second Kind Basis', SecondKind )

    IF (SecondOrder) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:2 -brick b:6 -pyramid b:3 -prism b:2 -quad_face b:4 -tri_face b:2")
    ELSE IF( SecondKind ) THEN
      CALL ListAddString( SolverParams, "Element", "n:0 e:2")
    ELSE IF (PiolaVersion) THEN
      CALL ListAddString(SolverParams, "Element", &
          "n:0 e:1 -brick b:3 -quad_face b:2")
    ELSE
      CALL ListAddString( SolverParams, "Element", "n:0 e:1")
    END IF
  END IF

  CALL ListAddString( SolverParams, 'Kernel Variable', 'P' )

  CALL ListAddLogical( SolverParams, 'Linear System Symmetric', .TRUE. )
  CALL ListAddString(  SolverParams, 'Linear System Solver', 'Iterative' )
  CALL ListAddString(  SolverParams, 'Linear System Preconditioning', 'ILU' )
  CALL ListAddInteger( SolverParams, 'Linear System Residual Output', 25 )
  CALL ListAddInteger( SolverParams, 'Linear System Max Iterations', 2000 )
  CALL ListAddString(  SolverParams, 'Linear System Iterative Method', 'CG' )
  CALL ListAddConstReal( SolverParams, 'Linear System Convergence Tolerance', 1.0d-9 )
  CALL ListAddLogical( SolverParams,"Hcurl Basis",.TRUE.)

  CALL ListCopyPrefixedKeywords(Model % Solvers(i) % Values, SolverParams, 'RemoveKernelComponent:')
  
!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponentT_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Apply the Helmholtz projection on a curl-conforming vector field A
!>  when the kernel component grad phi of A (with respect to the curl operator)
!>  has been computed by using the subroutine HelmholtzProjector. This solver
!>  generates the representation W of grad phi in terms of the curl-conforming
!>  basis and finally redefines A := A - W, with W = grad phi. 
!------------------------------------------------------------------------------
SUBROUTINE RemoveKernelComponentT(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: SolverPtr, KerSolverPtr
  TYPE(Element_t), POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  LOGICAL :: PiolaVersion, SecondOrder
  LOGICAL :: ConstantBulkMatrix

  INTEGER :: dim, PotDOFs
  INTEGER :: i, j, k, n, nd, n_pot, nd_pot, t
  INTEGER :: istat, active

  REAL(KIND=dp), ALLOCATABLE, TARGET :: Stiff(:,:), Force(:), PhiSol(:), &
                    SOL(:), F(:)
  REAL(KIND=dp) :: Norm
  CHARACTER(LEN=MAX_NAME_LEN) :: PotName, Name
  REAL(KIND=dp), POINTER :: SaveRHS(:)
  TYPE(Variable_t), POINTER :: v

  SAVE Stiff, Force, PhiSol, AllocationsDone

  CALL Info( 'RemoveKernelComponent', '--------------------------------------------------',Level=12 )
  CALL Info( 'RemoveKernelComponent', 'Making the vector potential to be divergence free!',Level=12 )
  
  CALL DefaultStart()

  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF (.NOT. AllocationsDone) THEN
    n = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(n),STIFF(n,n),PhiSol(n),STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'RemoveKernelComponent', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !
  ! Find the variable which is projected:
  !
  PotName = GetString(SolverParams, 'Potential Variable', Found)
  IF (.NOT. Found ) PotName = 'av'

  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (PotName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with potential variable > '&
        //TRIM(PotName)//' < not found!')
  END IF

  !
  ! Find the variable which defines the kernel component:
  !
  Name = GetString(SolverParams, 'Kernel Variable', Found)
  IF (.NOT. Found ) Name = 'phi'
  V => VariableGet( Mesh % Variables, Name )

  Found = ASSOCIATED(v)
   
  IF (.NOT. Found ) THEN
    CALL Fatal('RemoveKernelComponent', 'Solver associated with kernel variable > '&
        //TRIM(Name)//' < not found!')
  END IF
  
  !
  ! Find some parameters to inherit the vector FE basis as defined in the primary solver:
  !
  CALL EdgeElementStyle(SolverPtr % Values, PiolaVersion, QuadraticApproximation = SecondOrder )
  IF (PiolaVersion) CALL Info('RemoveKernelComponent', &
      'Using Piola-transformed finite elements', Level=5)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize(Solver)

  DO t=1,active
    Element => GetActiveElement(t)

    n = GetElementNOFNodes()
    nd = GetElementNOFDOFs()
   
    ! The DOF counts for the potential variable: 
    n_pot = n*SolverPtr % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)
    nd_pot = GetElementNOFDOFs(USolver=SolverPtr)

    IF (nd /= nd_pot-n_pot) CALL Fatal('RemoveKernelComponent', &
     'Potential variable DOFs count /= the solver DOFs count')

    CALL GetLocalSolution(PhiSol, Name)

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix( STIFF, FORCE, Element, n, nd, dim, PiolaVersion, &
                SecondOrder, PhiSol )
    
    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(STIFF, FORCE)
  END DO

  CALL DefaultDirichletBCs()
  Norm = DefaultSolve()  

  !
  ! Finally, redefine the potential variable:
  !
  n = SIZE(Solver % Variable % Perm(:))
  IF (n ==  SIZE(SolverPtr % Variable % Perm(:))) THEN
    DO i=Solver % Mesh % NumberOfNodes+1,n
      IF( ASSOCIATED( Solver % Mesh % PeriodicPerm ) ) THEN
        IF( Solver % Mesh % PeriodicPerm(i) > 0 ) CYCLE
      END IF
      
      j = Solver % Variable % Perm(i)
      IF (j<=0) CYCLE

      k = SolverPtr % Variable % Perm(i)
      IF (k<=0) THEN
        CALL Fatal('RemoveKernelComponent', &
          'The variable and potential permutations are nonmatching?')
      END IF

      SolverPtr % Variable % Values(k) = SolverPtr % Variable % Values(k) - &
          Solver % Variable % Values(j)
    END DO
  ELSE
    CALL Fatal('RemoveKernelComponent', 'The variable and potential permutations differ')  
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Force, Element, n, nd, dim, PiolaVersion, &
              SecondOrder, PhiSol )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    REAL(KIND=dp) :: PhiSol(:)
    LOGICAL :: PiolaVersion, SecondOrder
!------------------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    LOGICAL :: Stat

    INTEGER :: i, j, p, q, t, EdgeBasisDegree 

    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), A(3)
    REAL(KIND=dp) :: s, DetJ
    REAL(KIND=dp) :: WBasis(nd,3), CurlWBasis(nd,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes)

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    IF (SecondOrder) THEN
      EdgeBasisDegree = 2  
    ELSE
      EdgeBasisDegree = 1
    END IF
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
        EdgeBasisDegree=EdgeBasisDegree)

    IF( dim == 2 .AND. .NOT. PiolaVersion ) THEN
      CALL Fatal('RemoveKernelComponent', '"Use Piola Transform = True" needed in 2D')
    END IF

    
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
          RotBasis = CurlWBasis, USolver = SolverPtr )       

      s = detJ * IP % s(t)
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + s * SUM(WBasis(q,:) * WBasis(p,:))
        END DO
      END DO

      A = MATMUL( PhiSol(1:n), dBasisdx(1:n,:) )
      DO q=1,nd
        FORCE(q) = FORCE(q) + s * SUM(A * WBasis(q,:))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE RemoveKernelComponentT
!------------------------------------------------------------------------------

