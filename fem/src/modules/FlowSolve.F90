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
! ****************************************************************************/

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
!> Initialization of the main solver: AdvectionDiffusionSolver
!------------------------------------------------------------------------------
   SUBROUTINE FlowSolver_init( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params

     Params => GetSolverParams()
     CALL ListAddInteger( Params,'Time Derivative Order',1 )

   END SUBROUTINE FlowSolver_init


!------------------------------------------------------------------------------
!> Solver for the Navier-Stokes equation in various different coordinate systems.
!------------------------------------------------------------------------------
   SUBROUTINE FlowSolver( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
    USE NavierStokes
    USE NavierStokesGeneral
    USE NavierStokesCylindrical
    USE Adaptive
    USE DefUtils
    USE FreeSurface
    USE ElementDescription, ONLY: GetEdgeMap
!------------------------------------------------------------------------------
    IMPLICIT NONE

     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER :: StiffMatrix
     
     INTEGER :: i,j,k,n,nb,nd,t,iter,LocalNodes,istat,q,m

     TYPE(ValueList_t),POINTER :: Material, BC, BodyForce, Equation
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,UNorm,Gravity(3),AngularVelocity(3), &
       Tdiff,s,Relaxation,NewtonTol,NewtonUBound,NonlinearTol, &
       ReferencePressure=0.0, SpecificHeatRatio, &
       PseudoCompressibilityScale=1.0, FreeSTol, res, MaxNorm

     INTEGER :: NSDOFs,NewtonIter,NewtonMaxIter,NonlinearIter,FreeSIter

     TYPE(Variable_t), POINTER :: DensitySol, TimeVar
     TYPE(Variable_t), POINTER :: FlowSol, TempSol, MeshSol

     INTEGER, POINTER :: FlowPerm(:),TempPerm(:), MeshPerm(:)
     REAL(KIND=dp), POINTER :: FlowSolution(:), Temperature(:), &
        gWork(:,:), ForceVector(:), LayerThickness(:), &
           SurfaceRoughness(:),MeshVelocity(:)

     REAL(KIND=dp), POINTER :: TempPrev(:)
     REAL(KIND=DP), POINTER :: Pwrk(:,:,:)

     LOGICAL :: Stabilize,NewtonLinearization = .FALSE., GotForceBC, GotIt, &
                  MBFlag, Convect  = .TRUE., NormalTangential, RelaxBefore, &
                  divDiscretization, GradPDiscretization, ComputeFree=.FALSE., &
                  Transient, Rotating, AnyRotating, OutOfPlaneFlow=.FALSE.,&
                  RecheckNewton=.FALSE., ImplicitFrictionDirection=.FALSE., &
                  LegacyBubbles=.FALSE.

! Which compressibility model is used
     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, StabilizeFlag, VarName
     CHARACTER(LEN=MAX_NAME_LEN) :: LocalCoords, FlowModel, DirectionName
     INTEGER :: CompressibilityModel, ModelCoords, ModelDim, NoActive
     INTEGER :: body_id,bf_id,eq_id,DIM
     INTEGER :: MidEdgeNodes(12), BrickFaceMap(6,4)
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     INTEGER, POINTER :: EdgeMap(:,:)


     INTEGER, SAVE :: Timestep, SaveTimestep=-1
     REAL(KIND=dp), ALLOCATABLE, SAVE :: pDensity0(:), pDensity1(:)
!
     LOGICAL :: AllocationsDone = .FALSE., FreeSurfaceFlag, &
         PseudoPressureExists, PseudoCompressible, Bubbles, P2P1, &
         Porous =.FALSE., PotentialForce=.FALSE., Hydrostatic=.FALSE., &
         MagneticForce =.FALSE., UseLocalCoords, PseudoPressureUpdate, &
         AllIncompressible

     REAL(KIND=dp),ALLOCATABLE :: MASS(:,:),STIFF(:,:), LoadVector(:,:), &
       Viscosity(:),FORCE(:), TimeForce(:), PrevDensity(:),Density(:),   &
       U(:),V(:),W(:),MU(:),MV(:),MW(:), Pressure(:),Alpha(:),Beta(:),   &
       ExtPressure(:),PrevPressure(:), HeatExpansionCoeff(:),            &
       ReferenceTemperature(:), Permeability(:),Mx(:),My(:),Mz(:),       &
       LocalTemperature(:), GasConstant(:), HeatCapacity(:),             &
       LocalTempPrev(:),SlipCoeff(:,:), PseudoCompressibility(:),        &
       PseudoPressure(:), Drag(:,:), PotentialField(:),    &
       PotentialCoefficient(:)

     SAVE U,V,W,MASS,STIFF,LoadVector,Viscosity, TimeForce,FORCE,ElementNodes,  &
       Alpha,Beta,ExtPressure,Pressure,PrevPressure, PrevDensity,Density,       &
       AllocationsDone,LocalNodes, HeatExpansionCoeff,ReferenceTemperature,     &
       Permeability,Mx,My,Mz,LayerThickness, SlipCoeff, SurfaceRoughness,       &
       LocalTemperature, GasConstant, HeatCapacity, LocalTempPrev,MU,MV,MW,     &
       PseudoCompressibilityScale, PseudoCompressibility, PseudoPressure,       &
       PseudoPressureExists, Drag, PotentialField, PotentialCoefficient, &
       ComputeFree, Indexes

      REAL(KIND=dp) :: at,at0,at1,totat,st,totst
!------------------------------------------------------------------------------

     INTERFACE
        FUNCTION FlowBoundaryResidual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
          INTEGER :: Perm(:)
        END FUNCTION FlowBoundaryResidual

        FUNCTION FlowEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Edge
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2)
          INTEGER :: Perm(:)
        END FUNCTION FlowEdgeResidual

        FUNCTION FlowInsideResidual( Model,Element,Mesh,Quant,Perm,Fnorm ) RESULT(Indicator)
          USE Types
          TYPE(Element_t), POINTER :: Element
          TYPE(Model_t) :: Model
          TYPE(Mesh_t), POINTER :: Mesh
          REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
          INTEGER :: Perm(:)
        END FUNCTION FlowInsideResidual
     END INTERFACE
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solving the system
!------------------------------------------------------------------------------
     CALL Info('FlowSolver','Solving the Navier-Stokes equations',Level=6)


     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
          
     CALL DefaultStart()
     
!    Check for local coordinate system

     LocalCoords = GetString( Solver % Values, 'Solver Coordinate System', &
          UseLocalCoords )

     IF ( UseLocalCoords ) THEN
        ModelCoords = Coordinates
        ModelDim = Model % DIMENSION
        SELECT CASE ( LocalCoords )
           CASE( 'cartesian 2d' )
              Coordinates = 1
              Model % DIMENSION = 2
              CALL Info( 'FlowSolve', 'Solver Coordinate System is cartesian 2d', LEVEL=31 )
           CASE( 'cartesian 3d' )
              Coordinates = 1
              Model % DIMENSION = 3
              CALL Info( 'FlowSolve', 'Solver Coordinate System is cartesian 3d', LEVEL=31 )
           CASE( 'axi symmetric' )
              Coordinates = 4
              Model % DIMENSION = 2
              CALL Info( 'FlowSolve', 'Solver Coordinate System is axi symmetric', LEVEL=31 )
           CASE( 'cylindric symmetric' )
              Coordinates = 3
              Model % DIMENSION = 3
              CALL Info( 'FlowSolve', 'Solver Coordinate System is cylindric symmetric', LEVEL=31 )
           CASE DEFAULT
              CALL Warn( 'FlowSolve', 'Solver coordinate system not recognised, using original' )
        END SELECT
     END IF

     ! check for Flow model, one of 'full', 'no convection', 'stokes':
     ! ---------------------------------------------------------------
     Transient = TransientSimulation
     Convect = .TRUE.
     FlowModel = GetString( GetSolverParams(), 'Flow Model', Gotit )
     
     SELECT CASE(FlowModel)
     CASE('no convection')
       Convect = .FALSE.
     CASE('stokes')
       Convect = .FALSE.
       Transient = .FALSE.
     CASE DEFAULT
       FlowModel = 'full'
     END SELECT

     DIM = CoordinateSystemDimension()

     FlowSol => Solver % Variable
     NSDOFs         =  FlowSol % DOFs
     FlowPerm       => FlowSol % Perm
     FlowSolution   => FlowSol % Values

     VarName = GetVarName(FlowSol)

     LocalNodes = COUNT( FlowPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     TempSol => VariableGet( Solver % Mesh % Variables, 'Temperature' )
     IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm     => TempSol % Perm
       Temperature  => TempSol % Values
       IF( Transient ) THEN
         IF ( ASSOCIATED(TempSol % PrevValues) ) TempPrev => TempSol % PrevValues(:,1)
       END IF
     END IF

     MeshSol => VariableGet( Solver % Mesh % Variables, 'Mesh Velocity')
     NULLIFY( MeshVelocity )
     IF ( ASSOCIATED( MeshSol ) ) THEN
       MeshPerm     => MeshSol % Perm
       MeshVelocity => MeshSol % Values
     END IF

     DensitySol => VariableGet( Solver % Mesh % Variables, 'Density' )

!------------------------------------------------------------------------------

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     UNorm = Solver % Variable % Norm


     ! Enable keyword used also in IncompressibleNSVec and HydrostaticNSVec
     IF(ListGetLogical(Solver % Values,'Constant-Viscosity Start',GotIt )) THEN
       IF(.NOT. AllocationsDone) THEN
         CALL ListAddConstReal( Solver % Values,'Newtonian Viscosity Condition',1.0_dp)
       END IF
     END IF
     
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------     
     IF ( .NOT.AllocationsDone .OR. Solver % MeshChanged ) THEN

       N = Solver % Mesh % MaxElementDOFs
       
       IF( AllocationsDone ) THEN
          DEALLOCATE(                               &
               U,  V,  W,                           &
               MU, MV, MW,                          &
               Indexes,                             &
               Pressure,                            &
               PrevPressure,                        &
               PseudoCompressibility,               &
               PrevDensity,Density,                 &
               LayerThickness,                      &
               SurfaceRoughness,                    &
               Permeability,                        &
               Mx,My,Mz,                            &
               SlipCoeff, Drag,                     &
               TimeForce,FORCE, Viscosity,          &
               MASS,  STIFF,                        &
               HeatExpansionCoeff,                  &
               GasConstant, HeatCapacity,           &
               ReferenceTemperature,                & 
               LocalTempPrev, LocalTemperature,     &
               PotentialField, PotentialCoefficient, &
               LoadVector, Alpha, Beta, &
               ExtPressure, STAT=istat )
       END IF

       ALLOCATE( U(N),  V(N),  W(N),                     &
                 MU(N), MV(N), MW(N),                    &
                 Indexes( N ),                           &
                 Pressure( N ),                          &
                 PrevPressure( N ),                      &
                 PseudoCompressibility( N ),             &
                 PrevDensity(N),Density( N ),            &
                 LayerThickness(N),                      &
                 SurfaceRoughness(N),                    &
                 Permeability(N),                        &
                 Mx(N),My(N),Mz(N),                      &
                 SlipCoeff(3,N), Drag(3,N),              &
                 TimeForce( 2*NSDOFs*N ),                &
                 FORCE( 2*NSDOFs*N ), Viscosity( N ), &
                 MASS(  2*NSDOFs*N,2*NSDOFs*N ),&
                 STIFF( 2*NSDOFs*N,2*NSDOFs*N ),&
                 HeatExpansionCoeff(N),                  &
                 GasConstant( N ), HeatCapacity( N ),    &
                 ReferenceTemperature(N),                & 
                 LocalTempPrev(N), LocalTemperature(N),  &
                 PotentialField( N ), PotentialCoefficient( N ), &
                 LoadVector( 4,N ), Alpha( N ), Beta( N ), &
                 ExtPressure( N ), STAT=istat )

       Drag = 0.0d0
       NULLIFY(Pwrk) 

       PseudoPressureExists = .FALSE.
       AllIncompressible = .TRUE.

       DO k=1,Model % NumberOfMaterials
         Material => Model % Materials(k) % Values
         CompressibilityFlag = ListGetString( Material, &
             'Compressibility Model', GotIt)
         IF (.NOT. gotIt ) CYCLE
         IF( CompressibilityFlag == 'artificial compressible') THEN
           PseudoPressureExists = .TRUE.
           AllIncompressible = .FALSE.
         ELSE IF ( CompressibilityFlag == 'incompressible' ) THEN
           CONTINUE
         ELSE
           AllIncompressible = .FALSE.
         END IF
       END DO

       IF ( PseudoPressureExists ) THEN
          IF ( AllocationsDone ) THEN
             DEALLOCATE( PseudoPressure )
          END IF
          n = SIZE( FlowSolution ) / NSDOFs
          ALLOCATE( PseudoPressure(n),STAT=istat ) 
       END IF

!------------------------------------------------------------------------------
!     This hack is needed  cause of the fluctuating pressure levels
!------------------------------------------------------------------------------
       IF( AllIncompressible ) THEN
         CALL Info('FlowSolve','Enforcing relative pressure relaxation',Level=8)
         CALL ListAddNewLogical( Solver % Values,'Relative Pressure Relaxation',.TRUE.)
       END IF
       
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'FlowSolve','Memory allocation error, Aborting.' )
       END IF
       
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------


     TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
     Timestep = NINT(Timevar % Values(1))
     IF ( SaveTimestep /= Timestep ) THEN
       IF ( ALLOCATED(pDensity0) ) pDensity0 = pDensity1
       SaveTimestep=Timestep
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
     IF ( GotIt ) THEN
       Gravity = gWork(1:3,1)*gWork(4,1)
     ELSE
       Gravity    =  0.00D0
       Gravity(2) = -9.81D0
     END IF

     AnyRotating = ListCheckPresentAnyBodyForce(Model,'Angular Velocity') .OR. &
         ListCheckPresentAnyBodyForce(Model,'Angular Velocity 1') .OR. &
         ListCheckPresentAnyBodyForce(Model,'Angular Velocity 2') .OR. &
         ListCheckPresentAnyBodyForce(Model,'Angular Velocity 3') 
         
!------------------------------------------------------------------------------

     ! Different options are:
     ! 1) stabilized
     ! 2) legacy pubbles
     ! 3) p-bubbles (or bubbles given in element definion)
     ! 4) p2p1
     ! 5) vms

     P2P1 = .FALSE.
     Bubbles   = ListGetLogical( Solver % Values,'Bubbles',GotIt )
     Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )
     LegacyBubbles = .FALSE.
     
     StabilizeFlag = ListGetString( Solver % Values, &
           'Stabilization Method', GotIt )
     IF ( .NOT. GotIt ) THEN
       IF ( Stabilize ) THEN
         StabilizeFlag = 'stabilized'
       ELSE IF ( Bubbles  ) THEN
         StabilizeFlag = 'bubbles'
       ELSE IF(ListCheckPresent(Solver % Values,'Element')) THEN
         StabilizeFlag = 'bubbles'
         Bubbles = .TRUE.
       ELSE
         CALL Info('FlowSolver','Defaulting to "stabilized" method')
         StabilizeFlag = 'stabilized'
         Stabilize = .TRUE.
       END IF
     ELSE
       IF (StabilizeFlag == 'p2/p1' .OR. StabilizeFlag == 'p2p1') THEN
         P2P1 = .TRUE.
       ELSE IF( StabilizeFlag == 'bubbles' ) THEN
         LegacyBubbles = .NOT. ListCheckPresent(Solver % Values,'Element') 
         Bubbles = .TRUE.
       ELSE IF(StabilizeFlag == 'stabilized' ) THEN
         Stabilize = .TRUE.
       ELSE IF(StabilizeFlag == 'pbubbles' ) THEN
         Bubbles = .TRUE.
       ELSE IF(StabilizeFlag == 'vms' ) THEN
         CONTINUE
       ELSE
         CALL Fatal('FlowSolver','Unknown "stabilization method": '//TRIM(StabilizeFlag))
       END IF
     END IF
     
     IF( Stabilize .AND. Bubbles ) THEN
       CALL Fatal('FlowSolver','You cant have stabilization and bubbles both!')
     END IF     
     IF( LegacyBubbles ) THEN
       CALL Info('FlowSolver','Using legacy bubbles (as opposed to elemental ones!)',Level=8)
     END IF
       
     DivDiscretization = ListGetLogical( Solver % Values, &
              'Div Discretization', GotIt )

     GradPDiscretization = ListGetLogical( Solver % Values, &
              'Gradp Discretization', GotIt )

     NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance',minv=0.0d0 )

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance', minv=0.0d0 )

     !Option to switch back to picard if convergence exceeds certain tolerance
     NewtonUBound = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton Max Tolerance', GotIt )
     IF(GotIt) RecheckNewton = .TRUE.

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations', minv=0 )
     IF ( NewtonIter == 0 ) NewtonLinearization = .TRUE.

     !Option to switch back to picard after NewtonMaxIter iterations
     NewtonMaxIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton Max Iterations', GotIt )
     RecheckNewton = RecheckNewton .OR. GotIt

     IF (GetLogical( GetSolverParams(), &
         'Nonlinear System Reset Newton',  GotIt)) NewtonLinearization=.FALSE.

     NonlinearIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Max Iterations', minv=0 )

     IF ( .NOT. ListCheckPresent( Solver % Values, &
        'Nonlinear System Norm Dofs' ) ) THEN
       CALL ListAddInteger( Solver % Values, 'Nonlinear System Norm DOFs', NSDOFs-1 )
     END IF

     FreeSTol = ListGetConstReal( Solver % Values, &
        'Free Surface After Tolerance', GotIt, minv=0.0d0 )
     IF ( .NOT. GotIt ) FreeSTol = HUGE(1.0d0)

     FreeSIter = ListGetInteger( Solver % Values, &
        'Free Surface After Iterations', GotIt, minv=0 )
     IF ( .NOT. GotIt ) FreeSIter = 0

     MaxNorm = ListGetConstReal( Solver % Values, &
        'Nonlinear System Max Norm Return', GotIt, minv=0.0d0 )
     IF ( .NOT. GotIt ) MaxNorm = HUGE(1.0_dp)

     DirectionName = ListGetString(Solver %Values, 'Implicit Friction Direction Vector', ImplicitFrictionDirection)
     IF (ImplicitFrictionDirection) THEN
       CALL Info('FlowSolver','"Implicit Friction Direction Vector" set to: '//TRIM(DirectionName),Level=10)
     END IF

!------------------------------------------------------------------------------
!    Check if free surfaces present
!------------------------------------------------------------------------------
     FreeSurfaceFlag = .FALSE.
     DO i=1,Model % NumberOfBCs
       FreeSurfaceFlag = FreeSurfaceFlag .OR. ListGetLogical( &
          Model % BCs(i) % Values,'Free Surface', GotIt )
       IF ( FreeSurfaceFlag ) EXIT
     END DO

     CALL CheckCircleBoundary()
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

 
     ! Initialize the pressure to be used in artificial compressibility 
     IF(PseudoPressureExists) THEN
       PseudoPressure = FlowSolution(NSDOFs:SIZE(FlowSolution):NSDOFs)

       WRITE(Message,'(A,T25,E15.4)') 'PseudoPressure mean: ',&
           SUM(PseudoPressure)/SIZE(PseudoPressure)
       CALL Info('FlowSolve',Message,Level=5)

       PseudoCompressibilityScale = ListGetConstReal( Model % Simulation, &
           'Artificial Compressibility Scaling',GotIt)      

       IF(.NOT.GotIt) PseudoCompressibilityScale = 1.0
       IF(Transient) THEN
         PseudoCompressibilityScale = PseudoCompressibilityScale / dt
       END IF
       PseudoPressureUpdate = ListGetLogical( Model % Simulation, &
           'Pseudo Pressure Update',GotIt)
       IF (.NOT.GotIt) PseudoPressureUpdate = .FALSE.
     END IF

     DO iter=1,NonlinearIter

       IF (PseudoPressureExists .AND. PseudoPressureUpdate) &
          PseudoPressure = FlowSolution(NSDOFs:SIZE(FlowSolution):NSDOFs)

       at  = CPUTime()
       at0 = RealTime()
       at1 = RealTime()

       CALL Info( 'FlowSolve', ' ', Level=4 )
       CALL Info( 'FlowSolve', ' ', Level=4 )
       CALL Info( 'FlowSolve', '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'NAVIER-STOKES ITERATION', iter 
       CALL Info( 'FlowSolve',Message, Level=4 )
       CALL Info( 'FlowSolve','-------------------------------------', Level=4 )
       CALL Info( 'FlowSolve', ' ', Level=4 )
       CALL Info( 'FlowSolve','Starting Assembly...', Level=4 )

!------------------------------------------------------------------------------
       !CALL InitializeToZero( StiffMatrix, ForceVector )
       CALL DefaultInitialize() 
!------------------------------------------------------------------------------

       bf_id   = -1
       body_id = -1

       CALL StartAdvanceOutput( 'FlowSolve', 'Assembly: ' )
       NoActive = GetNOFActive()
       
       DO t = 1,NoActive
         
         CALL AdvanceOutput( t, NoActive )

         Element => GetActiveElement(t)
         NodeIndexes => Element % NodeIndexes

!------------------------------------------------------------------------------

         IF ( Element % BodyId /= body_id ) THEN
           body_id = Element % BodyId

           eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation', &
                   minv=1, maxv=Model % NumberOfEquations )
           Equation => Model % Equations(eq_id) % Values

           bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )
           IF( bf_id > 0 ) THEN
             BodyForce => Model % BodyForces(bf_id) % Values
           END IF

           IF ( FlowModel == 'full' ) THEN
             Convect = ListGetLogical( Equation,'NS Convect', GotIt )
             IF ( .NOT. GotIt ) Convect = .TRUE.
           ENDIF

           k = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOfMaterials )
           Material => Model % Materials(k) % Values

!------------------------------------------------------------------------------
           CompressibilityFlag = ListGetString( Material, &
               'Compressibility Model', GotIt)
           IF ( .NOT.GotIt ) CompressibilityModel = Incompressible
           PseudoCompressible = .FALSE.
!------------------------------------------------------------------------------
           SELECT CASE( CompressibilityFlag )
!------------------------------------------------------------------------------
             CASE( 'incompressible' )
               CompressibilityModel = Incompressible

             CASE( 'perfect gas', 'perfect gas equation 1' )
               CompressibilityModel = PerfectGas1

             CASE( 'thermal' )
               CompressibilityModel = Thermal

             CASE( 'user defined' )
               CompressibilityModel = UserDefined1

             CASE( 'pressure dependent' )
               CompressibilityModel = UserDefined2

             CASE( 'artificial compressible' )
               CompressibilityModel = Incompressible 
               PseudoCompressible = .TRUE.

             CASE DEFAULT
               CompressibilityModel = Incompressible
!------------------------------------------------------------------------------
           END SELECT
!------------------------------------------------------------------------------

           Gotit = .FALSE.
           IF ( bf_id > 0 ) THEN
             MagneticForce = ListGetLogical( BodyForce,'Lorentz Force', gotIt )
             Hydrostatic = ListGetLogical( BodyForce,'Hydrostatic Pressure',gotIt )
           END IF
           IF ( .NOT. GotIt ) THEN
             Hydrostatic = ListGetLogical( Equation,'Hydrostatic Pressure',gotIt )
           END IF
!------------------------------------------------------------------------------
           
           Rotating = .FALSE.
           IF( bf_id > 0 ) THEN
             gWork => ListGetConstRealArray( BodyForce,'Angular Velocity',GotIt)
             IF ( GotIt ) THEN
               IF( Coordinates == Cartesian ) THEN
                 AngularVelocity = gWork(1:3,1)
                 Rotating = .TRUE.
               ELSE
                 CALL Fatal('FlowSolve','Rotating coordinate implemented only for cartesian coordinates')
               END IF
             ELSE
               AngularVelocity = 0.0_dp
             END IF
           END IF
         END IF
!------------------------------------------------------------------------------

         n = GetElementNOFNodes()
         IF( Stabilize ) THEN
           nb = 0
         ELSE IF( LegacyBubbles ) THEN
           nb = n
         ELSE
           nb = GetElementNOFBDOFs()
         END IF
         nd = GetElementDOFs( Indexes )

         CALL GetElementNodes( ElementNodes )

         SELECT CASE( NSDOFs )
           CASE(3)
             U(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd))-2)
             V(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd))-1)
             W(1:nd) = 0.0_dp
             IF (bf_id > 0 ) THEN
               W(1:n)  = ListGetReal(BodyForce,'Out Of Plane Velocity',&
                    n, NodeIndexes(1:n),OutOfPlaneFlow)
               IF (.NOT.OutOfPlaneFlow) &
                    W(1:n) = 0.0_dp
             END IF
           CASE(4)
             U(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd))-3)
             V(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd))-2)
             W(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd))-1)
         END SELECT

         MU(1:nd) = 0.0d0
         MV(1:nd) = 0.0d0
         MW(1:nd) = 0.0d0
         IF ( ASSOCIATED( MeshVelocity ) ) THEN
            SELECT CASE( MeshSol % DOFs )
            CASE(2)
               IF ( ALL( MeshPerm( Indexes(1:nd) ) > 0 ) ) THEN
                  MU(1:nd) = MeshVelocity(2*MeshPerm(Indexes(1:nd))-1)
                  MV(1:nd) = MeshVelocity(2*MeshPerm(Indexes(1:nd))-0)
               END IF

            CASE(3)
               IF ( ALL( MeshPerm( NodeIndexes ) > 0 ) ) THEN
                  MU(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-2)
                  MV(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-1)
                  MW(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-0)
               END IF
            END SELECT
         END IF
         
         LocalTemperature = 0.0d0
         LocalTempPrev    = 0.0d0
         IF ( ASSOCIATED( TempSol ) ) THEN
            IF ( ALL( TempPerm(NodeIndexes) > 0 ) ) THEN
               LocalTemperature(1:nd) = Temperature( TempPerm(Indexes(1:nd)) )
               IF ( Transient .AND. CompressibilityModel /= Incompressible) THEN
                 LocalTempPrev(1:nd) = TempPrev( TempPerm(Indexes(1:nd)) )
               END IF
            END IF
         END IF
         ReferencePressure = 0.0d0

         PrevDensity = 0.0d0
         Density = 0.0d0
!------------------------------------------------------------------------------
         SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
           CASE( Incompressible )
!------------------------------------------------------------------------------
             Pressure(1:nd) = FlowSolution( NSDOFs*FlowPerm(Indexes(1:nd)) )
             Density(1:n) = GetReal( Material, 'Density' )

             IF(PseudoCompressible) THEN
               Pressure(1:n) = GetReal( Material,'Artificial Pressure', GotIt )
               IF(.NOT. GotIt) THEN
                 Pressure(1:nd) = PseudoPressure(FlowPerm(Indexes(1:nd))) 
               ELSE
                 Pressure(n+1:nd) = 0.0d0
               END IF
               PseudoCompressibility(1:n) = PseudoCompressibilityScale * &
                   GetReal(Material,'Artificial Compressibility', GotIt )
               IF(.NOT. gotIt) PseudoCompressibility(1:n) = 0.0d0
             END IF

!------------------------------------------------------------------------------
           CASE( PerfectGas1 )

              ! Use  ReferenceTemperature in .sif file for fixed temperature
              ! field. At the moment can not have both fixed T ideal gas and
              ! Boussinesq force:
              !-------------------------------------------------------------
              IF ( .NOT. ASSOCIATED( TempSol ) ) THEN
                 LocalTemperature(1:n) = GetReal( Material, &
                   'Reference Temperature' )
                 LocalTempPrev = LocalTemperature
              END IF

              HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity', GotIt )


              ! Read Specific Heat Ratio:
              !--------------------------
              SpecificHeatRatio = ListGetConstReal( Material, &
                     'Specific Heat Ratio', GotIt )
              IF ( .NOT.GotIt ) SpecificHeatRatio = 5.d0/3.d0


              ! For an ideal gas, \gamma, c_p and R are really a constant
              ! GasConstant is an array only since HeatCapacity formally is
              !------------------------------------------------------------
              GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) *  &
                   HeatCapacity(1:n) / SpecificHeatRatio


              ! For ideal gases take pressure deviation p_d as the
              ! dependent variable: p = p_0 + p_d
              ! Read p_0
              !---------------------------------------------------
              ReferencePressure = ListGetConstReal( Material, &
                      'Reference Pressure', GotIt )
              IF ( .NOT.GotIt ) ReferencePressure = 0.0d0

              Pressure(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd)))
              IF ( Transient ) THEN
                PrevPressure(1:nd) = Solver % Variable % PrevValues( &
                          NSDOFs*FlowPerm(Indexes(1:nd)),1 )
              END IF
              Density(1:n) = ( Pressure(1:n) + ReferencePressure ) / &
                 ( GasConstant(1:n) * LocalTemperature(1:n) )

           CASE( UserDefined1 )
             Pressure(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd)) )
             IF ( ASSOCIATED( DensitySol ) ) THEN
               Density(1:nd) = DensitySol % Values( DensitySol % Perm(Indexes(1:nd)) )
               IF ( Transient ) THEN
                  PrevDensity(1:nd) = DensitySol % PrevValues( &
                       DensitySol % Perm(Indexes(1:nd)),1)
                END IF
             ELSE
               Density(1:n) = GetReal( Material,'Density' )
               IF ( Transient ) THEN
                 IF (.NOT.ALLOCATED(pDensity0))  THEN
                   ALLOCATE(pDensity0(LocalNodes), &
                            pDensity1(LocalNodes))
                 END IF

                 IF ( Timestep==1 ) &
                     pDensity0(Indexes(1:n)) = Density(1:n)
                 pDensity1(Indexes(1:n)) = Density(1:n)
                 PrevDensity(1:n) = pDensity0(Indexes(1:n))
               END IF
             END IF

           CASE( UserDefined2 )
             Density(1:n) = GetReal( Material,'Density' )
             Pressure(1:nd) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd)))

           CASE( Thermal )
             Pressure(1:n) = FlowSolution(NSDOFs*FlowPerm(Indexes(1:nd)))

             HeatExpansionCoeff(1:n) = GetReal( Material, &
               'Heat Expansion Coefficient' )

             ReferenceTemperature(1:n) = GetReal( Material, &
               'Reference Temperature' )

             Density(1:n) = GetReal( Material,'Density' )

             IF( Transient ) THEN
               PrevDensity(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                  ( LocalTempPrev(1:n) - ReferenceTemperature(1:n) ) )
             END IF
             Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                  ( LocalTemperature(1:n) - ReferenceTemperature(1:n) ) )
!------------------------------------------------------------------------------
         END SELECT
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!        Read in porous media defs
!------------------------------------------------------------------------------
         Porous = ListGetLogical( Material,'Porous Media', GotIt)
         IF(Porous) THEN
           CALL GetRealArray( Material,  Pwrk,'Porous Resistivity',GotIt)
           
           IF( .NOT. GotIt ) THEN
             Drag( 1,1:n) = GetReal( Material,'Porous Resistivity 1',GotIt )
	     Drag( 2,1:n) = GetReal( Material,'Porous Resistivity 2',GotIt ) 
             IF( NSDOFs -1 > 2 ) THEN
   	       Drag( 3,1:n) = GetReal( Material,'Porous Resistivity 3',GotIt ) 
             END IF
           ELSE IF ( SIZE(Pwrk,1) == 1 ) THEN
             DO i=1,NSDOFs-1
               Drag( i,1:n ) = Pwrk( 1,1,1:n )
             END DO
           ELSE 
             DO i=1,MIN(NSDOFs,SIZE(Pwrk,1))
               Drag(i,1:n) = Pwrk(i,1,1:n)
             END DO
           END IF
         END IF

!------------------------------------------------------------------------------
!        Viscosity = Laminar viscosity
!------------------------------------------------------------------------------
         Viscosity(1:n) = GetReal( Material,'Viscosity' )

!------------------------------------------------------------------------------
!        Set body forces, if any
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         IF ( bf_id > 0 ) THEN
           HeatExpansionCoeff   = 0.0D0
           ReferenceTemperature = 0.0D0

!------------------------------------------------------------------------------
!          Boussinesq body force & gravity
!------------------------------------------------------------------------------
           IF ( ListGetLogical( BodyForce,'Boussinesq',gotIt) ) THEN

             HeatExpansionCoeff(1:n) = GetReal( Material, &
                 'Heat Expansion Coefficient' )
             
             ReferenceTemperature(1:n) = GetReal( Material, &
                 'Reference Temperature' )

             DO i=1,n
               k = TempPerm(NodeIndexes(i))
               IF ( k > 0 ) THEN
                 IF ( Hydrostatic ) THEN
                   Tdiff = 1 - HeatExpansionCoeff(i) * &
                      (Temperature(k) - ReferenceTemperature(i))

                   IF ( Tdiff <= 0.0D0 ) THEN
                      CALL Warn( 'FlowSolve','Zero or negative density.' )
                   END IF
                 ELSE
                   Tdiff = -HeatExpansionCoeff(i) * &
                               (Temperature(k) - ReferenceTemperature(i))
                 END IF
  
                 LoadVector(1,i)   = Gravity(1) * Tdiff
                 LoadVector(2,i)   = Gravity(2) * Tdiff
                 IF ( NSDOFs > 3 ) THEN
                   LoadVector(3,i) = Gravity(3) * Tdiff
                 END IF
               END IF
             END DO
           ELSE IF ( Hydrostatic ) THEN
             LoadVector(1,1:n)   = Gravity(1)
             LoadVector(2,1:n)   = Gravity(2)
             IF ( NSDOFs > 3 ) LoadVector(3,1:n) = Gravity(3)
           END IF
!------------------------------------------------------------------------------
           LoadVector(1,1:n) = LoadVector(1,1:n) + ListGetReal( BodyForce, &
               'Flow Bodyforce 1',n,NodeIndexes,gotIt )
           
           LoadVector(2,1:n) = LoadVector(2,1:n) + ListGetReal( BodyForce, &
               'Flow Bodyforce 2',n,NodeIndexes,gotIt )
           
           IF ( NSDOFs > 3 ) THEN
             LoadVector(3,1:n) = LoadVector(3,1:n) + ListGetReal( BodyForce, &
                 'Flow Bodyforce 3',n,NodeIndexes,gotIt )
           END IF

!------------------------------------------------------------------------------
           
           PotentialForce = ListGetLogical( BodyForce,'Potential Force',gotIt) 
           IF(PotentialForce) THEN
             PotentialField(1:n) = ListGetReal( BodyForce, &
                 'Potential Field',n,NodeIndexes)             
             PotentialCoefficient(1:n) = ListGetReal( BodyForce, &
                 'Potential Coefficient',n,NodeIndexes)
           END IF


        
!------------------------------------------------------------------------------
         END IF ! of body forces

!------------------------------------------------------------------------------
!
! NOTE: LoadVector is multiplied by density inside *Navier* routines
!
         IF ( Transient ) THEN
           SELECT CASE( CompressibilityModel )
           CASE( PerfectGas1 )
             IF ( ASSOCIATED( TempSol ) ) THEN
               DO i=1,n
                 k = TempPerm(NodeIndexes(i))
                 IF ( k > 0 ) THEN
                    LoadVector(NSDOFs,i) = LoadVector(NSDOFs,i) + &
                      ( Temperature(k) - TempPrev(k) ) / dt
                 END IF
               END DO
             END IF
           CASE( UserDefined1, Thermal )
              DO i=1,n
                LoadVector(NSDOFs,i) = LoadVector(NSDOFs,i) - &
                  ( Density(i) - PrevDensity(i) ) / (Density(i)*dt)
              END DO
           END SELECT
         END IF

!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         SELECT CASE(Coordinates)
         CASE( Cartesian )
!------------------------------------------------------------------------------
           SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
             CASE( Incompressible,PerfectGas1,UserDefined1,UserDefined2,Thermal)
!------------------------------------------------------------------------------
! Density needed for steady-state, also pressure for transient
!------------------------------------------------------------------------------
               CALL NavierStokesCompose( MASS,STIFF,FORCE, LoadVector, &
                   Viscosity,Density,U,V,W,MU,MV,MW,ReferencePressure+Pressure(1:n), &
                   LocalTemperature, Convect, StabilizeFlag, CompressibilityModel, &
                   PseudoCompressible, PseudoCompressibility, GasConstant, Porous, &
                   Drag, PotentialForce, PotentialField, PotentialCoefficient,  &
                   MagneticForce, Rotating, AngularVelocity, DivDiscretization, &
                   GradPDiscretization, NewtonLinearization, Element,n,ElementNodes)
!------------------------------------------------------------------------------
           END SELECT
!------------------------------------------------------------------------------

         CASE( Cylindric,CylindricSymmetric,AxisSymmetric )
! Same comments as Cartesian
!------------------------------------------------------------------------------
           SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
             CASE( Incompressible,PerfectGas1)
!------------------------------------------------------------------------------
               CALL NavierStokesCylindricalCompose( &
                   MASS,STIFF,FORCE, &
                   LoadVector, Viscosity,Density,U,V,W,MU,MV,MW, &
                   ReferencePressure+Pressure(1:n),LocalTemperature,&
                   Convect, StabilizeFlag, CompressibilityModel /= Incompressible, &
                   PseudoCompressible, PseudoCompressibility, GasConstant, Porous, Drag, &
                   PotentialForce, PotentialField, PotentialCoefficient, &
                   MagneticForce, divDiscretization, gradpDiscretization, NewtonLinearization, &
                   Element,n,ElementNodes )
!------------------------------------------------------------------------------
             CASE DEFAULT
               CALL Fatal('FlowSolver','Missing compressibility model in cylindrical coordinates')
               

           END SELECT
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
         CASE DEFAULT
!------------------------------------------------------------------------------

           SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
           CASE( Incompressible,PerfectGas1)
             
             CALL NavierStokesGeneralCompose( &
                 MASS,STIFF,FORCE, &
                 LoadVector, Viscosity,Density,U,V,W,MU,MV,MW,Stabilize, &
                 NewtonLinearization,Element,n,ElementNodes )
             
           CASE DEFAULT
             CALL Fatal('FlowSolver','Missing compressibility model in general coordinates')
             
           END SELECT
           
!------------------------------------------------------------------------------
         END SELECT

         ! We do not have stabilized formulation for compressible fluids. 
         IF ( CompressibilityModel /= Incompressible .AND. &
                 StabilizeFlag == 'stabilized' ) THEN
            Bubbles = .TRUE.
            StabilizeFlag = 'bubbles'
         END IF

         ! Internally P2P1 is dealt as special case of bubbles. 
         IF ( Element % TYPE % BasisFunctionDegree <= 1 .AND. P2P1 ) THEN
            Bubbles = .TRUE.
            StabilizeFlag = 'bubbles'
         END IF
         
         ! If bubbles are requested, but not in element formulation.
         IF ( nb==0 .AND. Bubbles ) nb = n
           
!------------------------------------------------------------------------------
!        If time dependent simulation, add mass matrix to global 
!        matrix and global RHS vector
!------------------------------------------------------------------------------
         TimeForce = 0.0_dp
         IF ( Transient ) THEN
!------------------------------------------------------------------------------
!          NOTE: the following will replace STIFF and FORCE
!          with the combined information
!------------------------------------------------------------------------------
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
         END IF
         
         IF ( nb > 0 ) THEN
           CALL NSCondensate( nd, nb, NSDOFs-1, STIFF, FORCE, TimeForce )
         END IF
         
!------------------------------------------------------------------------------
!        Add local stiffness matrix and force vector to global matrix & vector
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------

      CALL DefaultFinishBulkAssembly()

      CALL Info( 'FlowSolve', 'Assembly done', Level=4 )

!------------------------------------------------------------------------------
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------
      NoActive = GetNOFBoundaryElements()

      DO t = 1,NoActive

        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        IF( dim - GetElementDim(Element) > 1 ) CYCLE
        
        n = GetElementNOFNodes()

        CALL GetElementNodes( ElementNodes )
        NodeIndexes => Element % NodeIndexes

        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE

!------------------------------------------------------------------------------
        GotForceBC = GetLogical( BC, 'Flow Force BC',gotIt )
        IF ( .NOT. gotIt ) GotForceBC = .TRUE.

        IF ( GotForceBC ) THEN
          LoadVector  = 0.0d0
          Alpha       = 0.0d0
          ExtPressure = 0.0d0
          Beta        = 0.0d0
          SlipCoeff   = 0.0d0
          STIFF = 0.0d0
          FORCE = 0.0d0

!------------------------------------------------------------------------------
!         (at the moment the following is done...)
!         BC: \tau \cdot n = \alpha n +  @\beta/@t + R_k u_k + F
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!         normal force BC: \tau\cdot n = \alpha n
!------------------------------------------------------------------------------
          !IF ( GetLogical( BC, 'Free Surface',gotIt) ) THEN
            Alpha(1:n) = GetReal( BC,'Surface Tension Coefficient', gotIt )
          !END IF

          ExtPressure(1:n) = GetReal( BC, 'External Pressure', GotForceBC )
          IF(.NOT. GotForceBC) ExtPressure(1:n) = GetReal( BC, 'Normal Pressure', GotForceBC )
!------------------------------------------------------------------------------
!         tangential force BC:
!         \tau\cdot n = @\beta/@t (tangential derivative of something)
!------------------------------------------------------------------------------
              
          IF ( ASSOCIATED( TempSol ) ) THEN
            Beta(1:n) = GetReal( BC, &
                'Surface Tension Expansion Coefficient',gotIt )

            IF ( gotIt ) THEN
              DO j=1,n
                k = TempPerm( NodeIndexes(j) )
                IF ( k>0 ) Beta(j) = 1.0_dp - Beta(j) * Temperature(k)
              END DO
              Beta(1:n) = Beta(1:n) * GetReal(BC, 'Surface Tension Coefficient' )
            ELSE
              Beta(1:n) = GetReal( BC,'Surface Tension Coefficient', gotIt ) 
            END IF
          END IF

!------------------------------------------------------------------------------
!         force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------

          LoadVector(1,1:n) =  GetReal( BC, 'Pressure 1', GotIt )
          LoadVector(2,1:n) =  GetReal( BC, 'Pressure 2', GotIt )
          LoadVector(3,1:n) =  GetReal( BC, 'Pressure 3', GotIt )
          LoadVector(4,1:n) =  GetReal( BC, 'Mass Flux', GotIt )

!------------------------------------------------------------------------------
!         slip boundary condition BC: \tau\cdot n = R_k u_k
!------------------------------------------------------------------------------

          SlipCoeff = 0.0d0
          SlipCoeff(1,1:n) =  GetReal( BC, 'Slip Coefficient 1',GotIt )
          SlipCoeff(2,1:n) =  GetReal( BC, 'Slip Coefficient 2',GotIt )
          SlipCoeff(3,1:n) =  GetReal( BC, 'Slip Coefficient 3',GotIt )

          NormalTangential = GetLogical( BC, &
                 'Normal-Tangential Velocity', GotIt )

          IF (.NOT.GotIt) THEN
            NormalTangential = GetLogical( BC, &
                   'Normal-Tangential '//GetVarName(Solver % Variable), GotIt )
          END IF
!------------------------------------------------------------------------------
          SELECT CASE( CurrentCoordinateSystem() )
            
          CASE( Cartesian )            
            CALL NavierStokesBoundary(  STIFF, FORCE, &
                LoadVector, Alpha, Beta, ExtPressure, SlipCoeff, NormalTangential,   &
                Element, n, ElementNodes )
            
          CASE( Cylindric, CylindricSymmetric,  AxisSymmetric )            
            CALL NavierStokesCylindricalBoundary( STIFF, &
                FORCE, LoadVector, Alpha, Beta, ExtPressure, SlipCoeff, &
                NormalTangential, Element, n, ElementNodes)
            
          CASE DEFAULT            
            CALL NavierStokesGeneralBoundary( STIFF, &
                FORCE, LoadVector, Alpha, Beta, ExtPressure, SlipCoeff, &
                Element, n, ElementNodes)
            
         END SELECT

!------------------------------------------------------------------------------

          IF ( GetLogical( BC, 'Wall Law',GotIt ) ) THEN
            !/*
            ! * TODO: note that the following is not really valid, the
            ! * pointer to the Material structure is from the remains
            ! * of the last of the bulk elements.
            ! */
            Density(1:n)   = GetParentMatProp( 'Density' )
            Viscosity(1:n) = GetParentMatProp( 'Viscosity' )

            LayerThickness(1:n) = GetReal( BC, 'Boundary Layer Thickness' )
            SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness',GotIt )

            SELECT CASE( NSDOFs )
              CASE(3)
                U(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-2 )
                V(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-1 )
                W(1:n) = 0.0d0
            
              CASE(4)
                U(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-3 )
                V(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-2 )
                W(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-1 )
            END SELECT

            CALL NavierStokesWallLaw( STIFF,FORCE,     &
              LayerThickness,SurfaceRoughness,Viscosity,Density,U,V,W, &
                     Element,n, ElementNodes )
          ELSE IF ( GetLogical( BC, 'VMS Wall', GotIt ) ) THEN
            Density(1:n)   = GetParentMatProp( 'Density' )
            Viscosity(1:n) = GetParentMatProp( 'Viscosity' )

            LayerThickness(1:n) = GetReal( BC, 'Boundary Layer Thickness', GotIt )
            SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness',GotIt )

            SELECT CASE( NSDOFs )
              CASE(3)
                U(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-2 )
                V(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-1 )
                W(1:n) = 0.0d0
            
              CASE(4)
                U(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-3 )
                V(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-2 )
                W(1:n) = FlowSolution( NSDOFs*FlowPerm(NodeIndexes)-1 )
            END SELECT
            CALL VMSWalls( STIFF,FORCE,     &
              LayerThickness,SurfaceRoughness,Viscosity,Density,U,V,W, &
                     Element,n, ElementNodes )

          END IF
!------------------------------------------------------------------------------

          IF ( Transient ) THEN
            MASS = 0.0d0
            CALL Default1stOrderTime( MASS, STIFF, FORCE )
          END IF

!------------------------------------------------------------------------------
!         Add local stiffness matrix and force vector to
!         global matrix & vector
!--------------------------------------------------------------------------
          CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
        END IF
      END DO
!------------------------------------------------------------------------------
      !
      ! IMPLEMENT NOSLIP WALL BC CODE:
      ! ------------------------------
      DO i=1,Model % NumberOFBCs
        BC => Model % BCs(i) % Values
        IF ( GetLogical(  BC, 'Noslip wall BC', gotit ) ) THEN
          IF ( VarName  == 'flow solution' ) THEN
            CALL ListAddConstReal( BC, 'Velocity 1', 0.0_dp )
            CALL ListAddConstReal( BC, 'Velocity 2', 0.0_dp )
            IF ( NSDOFs>3 ) CALL ListAddConstReal( BC, 'Velocity 3', 0.0_dp )
          ELSE
            DO j=1,NSDOFs-1
              CALL ListAddConstReal( BC, ComponentName( &
                 Solver % Variable % Name, j), 0.0_dp )
            END DO
          END IF
        END IF
      END DO

      CALL DefaultFinishBoundaryAssembly()
     
      !------------------------------------------------------------------------------
      !     Implicit Friction Boundaries
      !------------------------------------------------------------------------------     
      IF (ImplicitFrictionDirection) THEN
        ! This is a matrix level routine for setting friction such that tangential
        ! traction is the normal traction multiplied by a coefficient.
        CALL SetImplicitFriction(Model, Solver,'Implicit Friction Coefficient', DirectionName )
      ELSE
        CALL SetImplicitFriction(Model, Solver,'Implicit Friction Coefficient' )
      END IF
      
      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()
      CALL Info( 'FlowSolve', 'Dirichlet conditions done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      at = CPUTime() - at
      st = CPUTime()

      Unorm = DefaultSolve()

      ! If we have constant viscosity start then remove it already after first solution.
      IF(ListGetLogical(Solver % Values,'Constant-Viscosity Start',GotIt )) THEN
        CALL ListRemove( Solver % Values,'Newtonian Viscosity Condition')
      END IF
      
      st = CPUTIme()-st
      totat = totat + at
      totst = totst + st
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
      CALL Info( 'FlowSolve', Message, Level=4 )
      WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
      CALL Info( 'FlowSolve', Message, Level=4 )

      n = NSDOFs * LocalNodes


!------------------------------------------------------------------------------

      WRITE( Message, * ) 'Result Norm     : ',Solver % Variable % Norm
      CALL Info( 'FlowSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',Solver % Variable % NonlinChange
      CALL Info( 'FlowSolve', Message, Level=4 )
      
      RelativeChange = Solver % Variable % NonlinChange
      IF ( RelativeChange < NewtonTol .OR. &
             iter > NewtonIter ) NewtonLinearization = .TRUE.

      ! Check whether we should revert back to the more robust Picard iteration
      IF ( RecheckNewton .AND. NewtonLinearization ) THEN
        IF(RelativeChange > NewtonUBound) THEN
          NewtonLinearization = .FALSE.
          CALL Info('FlowSolve', 'Newton tolerance exceeded, switching back to picard', Level=6)
        END IF
        IF ( iter >= NewtonMaxIter) THEN
          NewtonLinearization = .FALSE.
          CALL Info('FlowSolve', 'Newton iteration limit exceeded, switching back to picard', Level=6)
        END IF
        IF ( RelativeChange > NewtonUBound) THEN
          NewtonLinearization = .FALSE.
          CALL Info('FlowSolve', 'Newton tolerance exceeded, switching back to picard', Level=6)
        END IF
      END IF
        
      IF ( RelativeChange < NonLinearTol .AND. Iter<NonlinearIter ) EXIT
      IF ( Solver % Variable % Norm > MaxNorm) THEN
         CALL Warn('FlowSolve', 'Exiting as nonlinear norm is above allowed maximum!')
         EXIT
      END IF

!------------------------------------------------------------------------------
!     If free surfaces in model, this will move the nodal points
!------------------------------------------------------------------------------
      IF ( FreeSurfaceFlag ) THEN
        IF ( RelativeChange < FreeSTol .OR. iter > FreeSIter ) ComputeFree = .TRUE.

        IF ( ComputeFree ) THEN
          MBFlag = GetLogical(GetSolverParams(), 'Internal Move Boundary', GotIt)
          IF ( MBFlag .OR. .NOT. GotIt ) THEN
            Relaxation = GetCReal( Solver % Values, &
                'Free Surface Relaxation Factor', GotIt )            
            IF ( .NOT.GotIt ) Relaxation = 1.0_dp
            CALL MoveBoundary( Model, Relaxation )
          END IF
        END IF
      END IF
!------------------------------------------------------------------------------
    END DO  ! of nonlinear iteration

    IF ( P2P1 ) THEN
      !----------------------------------------------------------------------------------------
      ! Replace the zero pressure solution at the nodes which are not needed in the linear
      ! pressure approximation by the interpolated values for right visualization:
      !----------------------------------------------------------------------------------------
      NoActive = GetNOFActive()

      DO t=1,NoActive
        ! First the midedge nodes:
        Element => GetActiveElement(t)
        IF ( Element % TYPE % BasisFunctionDegree <= 1 ) CYCLE

        nd = GetElementDOFs( Indexes )
        k = GetElementFamily()
        EdgeMap => GetEdgeMap(k)
        SELECT CASE(k)
        CASE (3)
          MidEdgeNodes(1:3) = (/ 4, 5, 6 /)
        CASE (4)
          MidEdgeNodes(1:4) = (/ 5, 6, 7, 8 /)
        CASE (5)
          MidEdgeNodes(1:6) = (/ 5, 6, 7, 8, 9, 10 /) 
        CASE (6)
          MidEdgeNodes(1:8) = (/ 6, 7, 8, 9, 10, 11, 12, 13 /)
        CASE (7)
          MidEdgeNodes(1:9) = (/ 7, 8, 9, 10, 11, 12, 13, 14, 15 /)
        CASE (8)
          MidEdgeNodes(1:12) = (/ 9, 10, 11, 12, 17, 18, 19, 20, 13, 14, 15, 16 /)
        END SELECT

        DO q=1,SIZE(EdgeMap,1)
          m = (dim+1) * Solver % Variable % Perm(Indexes(MidEdgeNodes(q)))
          i = (dim+1) * Solver % Variable % Perm(Indexes(EdgeMap(q,1)))
          j = (dim+1) * Solver % Variable % Perm(Indexes(EdgeMap(q,2)))
          Solver % Variable % Values(m) = 0.5d0 * ( Solver % Variable % Values(i) + &
              Solver % Variable % Values(j) )   
        END DO

        ! The pressure at the midface nodes for 409 elements: 
        IF (k==4 .AND. nd==9) THEN
          res = 0.0d0
          DO q=1,4
            i = (dim+1) * Solver % Variable % Perm(Indexes(q))
            res = res + Solver % Variable % Values(i)
          END DO
          m = (dim+1) * Solver % Variable % Perm(Indexes(9))
          Solver % Variable % Values(m) = 0.25d0 * res
        END IF

        ! The pressure at the midpoint and at the midface nodes for 827 elements:
        IF (k==8 .AND. nd==27) THEN
          BrickFaceMap(1,:) = (/ 1,2,6,5 /)         
          BrickFaceMap(2,:) = (/ 2,3,7,6 /)
          BrickFaceMap(3,:) = (/ 4,3,7,8 /)
          BrickFaceMap(4,:) = (/ 1,4,8,5 /)
          BrickFaceMap(5,:) = (/ 1,2,3,4 /)          
          BrickFaceMap(6,:) = (/ 5,6,7,8 /)
          DO j=1,6
            res = 0.0d0
            DO q=1,4
              i = (dim+1) * Solver % Variable % Perm(Indexes(BrickFaceMap(j,q)))
              res = res + Solver % Variable % Values(i)
            END DO
            m = (dim+1) * Solver % Variable % Perm(Indexes(20+j))
            Solver % Variable % Values(m) = 0.25d0 * res
          END DO

          res = 0.0d0
          DO q=1,8
            i = (dim+1) * Solver % Variable % Perm(Indexes(q))
            res = res + Solver % Variable % Values(i)
          END DO
          m = (dim+1) * Solver % Variable % Perm(Indexes(27))      
          Solver % Variable % Values(m) = 0.125d0 * res
        END IF

      END DO
    END IF

    IF (ListGetLogical(Solver % Values,'Adaptive Mesh Refinement',GotIt)) &
      CALL RefineMesh( Model,Solver,FlowSolution,FlowPerm, &
         FlowInsideResidual, FlowEdgeResidual, FlowBoundaryResidual ) 

!------------------------------------------------------------------------------
    CALL CheckCircleBoundary()
!------------------------------------------------------------------------------

    IF ( UseLocalCoords ) THEN
       Coordinates = ModelCoords
       Model % DIMENSION = ModelDim
    END IF

CONTAINS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CheckCircleBoundary()
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: x,y,phi,x0,y0,r
      LOGICAL :: GotIt
      INTEGER :: i,j,k,l
!------------------------------------------------------------------------------

      l = 0
      DO i=1,Model % NumberOfBCs
         IF ( .NOT.ListgetLogical( Model % BCs(i) % Values, &
                  'Circle Boundary', GotIt ) ) CYCLE

         x0 = ListGetConstReal( Model % BCs(i) % Values, 'Circle X', GotIt )
         IF ( .NOT. GotIt ) x0 = 0.0d0

         y0 = ListGetConstReal( Model % BCs(i) % Values, 'Circle Y', GotIt )
         IF ( .NOT. GotIt ) y0 = 0.0d0

         R  = ListGetConstReal( Model % BCs(i) % Values, 'Circle R', GotIt )
         IF ( .NOT. GotIt ) R = 1.0d0

         DO j=Solver % Mesh % NumberOfBulkElements+1, &
            Solver % Mesh % NumberOfBulkElements+ &
               Solver % Mesh % NumberOfBoundaryElements
            Element => Solver % Mesh % Elements(j)
            IF ( Element % BoundaryInfo % Constraint &
                 /= Model % BCs(i) % Tag ) CYCLE

            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes
            DO k=1,n
               x = Solver % Mesh % Nodes % x(NodeIndexes(k)) - x0
               y = Solver % Mesh % Nodes % y(NodeIndexes(k)) - y0

               phi = ATAN2( y,x )
               x = R * COS( phi ) 
               y = R * SIN( phi ) 

               Solver % Mesh % Nodes % x(NodeIndexes(k)) = x + x0
               Solver % Mesh % Nodes % y(NodeIndexes(k)) = y + y0
            END DO
            l = l + 1
        END DO
     END DO

     IF ( l > 0 ) THEN
        WRITE( Message, * ) 'Elements on Circle', l
        CALL Info( 'FlowSolve', Message, Level=6 )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CheckCircleBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE FlowSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Compute the residual of the Navier-Stokes equation for the boundary elements.
!------------------------------------------------------------------------------
  FUNCTION FlowBoundaryResidual( Model, Edge, Mesh, &
        Quant, Perm, Gnorm ) RESULT( Indicator )
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
     TYPE(Element_t), POINTER :: Element

     INTEGER :: i,k,n,l,t,bc,DIM,DOFs,Pn,En
     LOGICAL :: stat, GotIt, Compressible

     REAL(KIND=dp) :: Grad(3,3), Grad1(3,3), Stress(3,3), Normal(3), ForceSolved(3), &
                      EdgeLength,u, v, w, s, detJ
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), dEdgeBasisdx(:,:), Basis(:),dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE ::  x(:), y(:), z(:), ExtPressure(:)
     REAL(KIND=dp), ALLOCATABLE :: Temperature(:), Tension(:), SlipCoeff(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Velocity(:,:), Pressure(:), Force(:,:), NodalViscosity(:)

     REAL(KIND=dp) :: Residual(3), ResidualNorm, Viscosity, Slip, Dir(3)

     TYPE(Variable_t), POINTER  :: TempSol
     TYPE(ValueList_t), POINTER :: Material
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Gnorm     = 0.0d0

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

     DOFs = DIM + 1
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

     ALLOCATE( EdgeBasis(En), dEdgeBasisdx(En,3), Basis(Pn), dBasisdx(Pn,3), &
      x(En), y(En), z(En), ExtPressure(En), Temperature(Pn), Tension(En),    &
      SlipCoeff(3,En), Velocity(3,Pn), Pressure(Pn), Force(3,En), NodalViscosity(En) )

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

     DO bc=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(bc) % Values, &
!                 'Flow Force BC', gotIt ) ) CYCLE
!
!       Get material parameters:
!       ------------------------

        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )
        Material => Model % Materials(k) % Values

        NodalViscosity(1:En) = ListGetReal( Material, &
                 'Viscosity', En, Edge % NodeIndexes, GotIt )

        Compressible = .FALSE.
        IF ( ListGetString( Material, 'Compressibility Model', GotIt ) == &
               'perfect gas equation 1' ) Compressible = .TRUE.


!       Given traction:
!       ---------------
        Force = 0.0d0

        Force(1,1:En) = ListGetReal( Model % BCs(bc) % Values, &
            'Pressure 1', En, Edge % NodeIndexes, GotIt )

        Force(2,1:En) = ListGetReal( Model % BCs(bc) % Values, &
            'Pressure 2', En, Edge % NodeIndexes, GotIt )

        Force(3,1:En) = ListGetReal( Model % BCs(bc) % Values, &
            'Pressure 3', En, Edge % NodeIndexes, GotIt )

!
!       Force in normal direction:
!       ---------------------------
        ExtPressure(1:En) = ListGetReal( Model % BCs(bc) % Values, &
          'External Pressure', En, Edge % NodeIndexes, GotIt )

!
!       Slip BC condition:
!       ------------------
        SlipCoeff = 0.0d0
        SlipCoeff(1,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
             'Slip Coefficient 1',En,Edge % NodeIndexes,GotIt )

        SlipCoeff(2,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
             'Slip Coefficient 2',En,Edge % NodeIndexes,GotIt )

        SlipCoeff(3,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
             'Slip Coefficient 3',En,Edge % NodeIndexes,GotIt )

!
!       Surface tension induced by temperature gradient (or otherwise):
!       ---------------------------------------------------------------
        TempSol => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )

        IF ( ASSOCIATED( TempSol ) ) THEN
          Tension(1:En) = ListGetReal( Model % BCs(bc) % Values, &
           'Surface Tension Expansion Coefficient',En,Edge % NodeIndexes,gotIt )

           IF ( gotIt ) THEN
              DO n=1,En
                 k = TempSol % Perm( Edge % NodeIndexes(n) )
                 IF (k>0) Tension(n) = 1.0d0 - Tension(n) * TempSol % Values(k)
              END DO

              Tension(1:En) = Tension(1:En) * ListGetReal( &
                 Model % BCs(bc) % Values,'Surface Tension Coefficient', &
                               En, Edge % NodeIndexes ) 
           ELSE
              Tension(1:En) = ListGetReal( &
                  Model % BCs(bc) % Values,'Surface Tension Coefficient', &
                         En, Edge % NodeIndexes,gotIt ) 
           END IF
        ELSE
           Tension(1:En) = ListGetReal( &
               Model % BCs(bc) % Values,'Surface Tension Coefficient', &
                      En, Edge % NodeIndexes,gotIt ) 
        END IF

!
!       If dirichlet BC for velocity in any direction given,
!       nullify force in that direction:
!       ------------------------------------------------------------------
        Dir = 1
        s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 1', GotIt )
        IF ( GotIt ) Dir(1) = 0

        s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 2', GotIt )
        IF ( GotIt ) Dir(2) = 0

        s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 3', GotIt )
        IF ( GotIt ) Dir(3) = 0

!
!       Elementwise nodal solution:
!       ---------------------------
        Velocity = 0.0d0
        DO k=1,DOFs-1
           Velocity(k,1:Pn) = Quant(DOFs*Perm(Element % NodeIndexes)-DOFs + k)
        END DO
        Pressure(1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes) )

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

           Viscosity = SUM( NodalViscosity(1:En) * EdgeBasis(1:En) )

           Residual = 0.0d0
!
!          Given force at the integration point:
!          -------------------------------------
           Residual = Residual + MATMUL( Force(:,1:En), EdgeBasis(1:En) ) - &
                 SUM( ExtPressure(1:En) * EdgeBasis(1:En) ) * Normal

!
!          Slip velocity BC:
!          -----------------
           DO i=1,DIM
              Slip = SUM( SlipCoeff(i,1:En) * EdgeBasis(1:En) )
              Residual(i) = Residual(i) - &
                   Slip * SUM( Velocity(i,1:Pn) * Basis(1:Pn) )
           END DO

!
!          Tangential tension force:
!          -------------------------
           DO i=1,DIM
              Residual(i) = Residual(i) + &
                   SUM( dEdgeBasisdx(1:En,i) * Tension(1:En) )
           END DO

!
!          Force given by the computed solution:
!          -------------------------------------
!
!          Stress tensor on the boundary:
!          ------------------------------
           Grad = MATMUL( Velocity(:,1:Pn), dBasisdx(1:Pn,:) )

           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
              Grad1 = Grad
              DO i=1,DIM
                 DO k=1,DIM
                    DO l=1,DIM
                       Grad1(i,k) = Grad1(i,k) - &
                          Symb(k,l,i) * SUM ( Velocity(l,1:Pn) * Basis(1:Pn) )
                    END DO
                 END DO
              END DO

              Grad = 0.0d0
              DO i=1,DIM
                 DO k=1,DIM
                    DO l=1,DIM
                       Grad(i,k) = Grad(i,k) + Metric(k,l) * Grad1(i,l)
                    END DO
                 END DO
              END DO
           END IF

           Stress = Viscosity * ( Grad + TRANSPOSE(Grad) )
           Stress = Stress - Metric * SUM( Pressure(1:Pn) * Basis(1:Pn) )

           IF ( Compressible ) THEN
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 DO i=1,DIM
                    DO k=1,DIM
                       Stress(i,i) = Stress(i,i) - &
                           (2.0d0/3.0d0) * Viscosity * Grad(k,k)
                    END DO
                 END DO
              ELSE
                 DO i=1,DIM
                    DO k=1,DIM
                       DO l=1,DIM
                          Stress(i,k) = Stress(i,k) - &
                             Metric(i,k) * (2.0d0/3.0d0) * Viscosity * Grad(l,l)
                       END DO
                    END DO
                 END DO
              END IF
           END IF

           ForceSolved = MATMUL(Stress,Normal)
           Residual = Residual - ForceSolved * Dir

           EdgeLength = EdgeLength + s

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Gnorm = Gnorm + s * SUM( ForceSolved**2 )
              ResidualNorm = ResidualNorm + s * SUM( Residual(1:DIM) ** 2 )
           ELSE
              CALL InvertMatrix( Metric,3 )
              DO i=1,DIM
                 DO k=1,DIM
                    ResidualNorm = ResidualNorm + &
                            s * Metric(i,k) * Residual(i) * Residual(k)
                    Gnorm = GNorm + s * Metric(i,k) * &
                                        ForceSolved(i) * ForceSolved(k)
                 END DO
              END DO
           END IF
        END DO
        EXIT
     END DO

     IF ( CoordinateSystemDimension() == 3 ) EdgeLength = SQRT(EdgeLength)
     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( EdgeBasis, dEdgeBasisdx, Basis, dBasisdx, x, y, z,   &
      ExtPressure, Temperature, Tension,SlipCoeff, Velocity, Pressure, &
      Force, NodalViscosity )
!------------------------------------------------------------------------------
  END FUNCTION FlowBoundaryResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the residual of the Navier-Stokes equation for the edge elements.
!------------------------------------------------------------------------------
  FUNCTION FlowEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT( Indicator )
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
     TYPE(Element_t), POINTER :: Element

     INTEGER :: i,j,k,l,n,t,DIM,DOFs,En,Pn
     LOGICAL :: stat, GotIt

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Stress(3,3,2), Jump(3), Viscosity

     REAL(KIND=dp) :: Grad(3,3), Grad1(3,3), Normal(3)

     REAL(KIND=dp), ALLOCATABLE :: NodalViscosity(:), x(:), y(:), z(:), &
               EdgeBasis(:), Basis(:), dBasisdx(:,:)
     REAL(KIND=dp), ALLOCATABLE :: Velocity(:,:), Pressure(:)

     REAL(KIND=dp) :: ResidualNorm, EdgeLength, u, v, w, s, detJ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT

     DOFs = DIM + 1
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs - 1

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

     ALLOCATE( NodalViscosity(En), x(En), y(En), z(En), EdgeBasis(En), &
           Basis(n), dBasisdx(n,3), Velocity(3,n), Pressure(n) )

!    Integrate square of jump over edge:
!    ------------------------------------
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

           k = ListGetInteger( Model % Bodies( Element % BodyId) % Values, 'Material', &
                            minv=1, maxv=Model % NumberOfMaterials )

           NodalViscosity(1:En) = ListGetReal( &
               Model % Materials(k) % Values, 'Viscosity', &
                    En, Edge % NodeIndexes, GotIt )

           Viscosity = SUM( NodalViscosity(1:En) * EdgeBasis(1:En) )
!
!          Elementwise nodal solution:
!          ---------------------------
           Velocity = 0.0d0
           DO k=1,DOFs-1
              Velocity(k,1:Pn) = Quant(DOFs*Perm(Element % NodeIndexes)-DOFs+k)
           END DO
           Pressure(1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes) )
!
!          Stress tensor on the edge:
!          --------------------------
           Grad = MATMUL( Velocity(:,1:Pn), dBasisdx(1:Pn,:) )

           IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
              Grad1 = Grad
              DO j=1,DIM
                 DO k=1,DIM
                    DO l=1,DIM
                       Grad1(j,k) = Grad1(j,k) - &
                          Symb(k,l,j) * SUM ( Velocity(l,1:Pn) * Basis(1:Pn) )
                    END DO
                 END DO
              END DO

              Grad = 0.0d0
              DO j=1,DIM
                 DO k=1,DIM
                    DO l=1,DIM
                       Grad(j,k) = Grad(j,k) + Metric(k,l) * Grad1(j,l)
                    END DO
                 END DO
              END DO
           END IF

           Stress(:,:,i) = Viscosity * ( Grad + TRANSPOSE(Grad) )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO j=1,DIM
                 Stress(j,j,i) = Stress(j,j,i) - SUM( Pressure(1:Pn) * Basis(1:Pn))
                 DO k=1,DIM
                    Stress(j,j,i) = Stress(j,j,i) - (2.0d0/3.0d0)*Viscosity*Grad(k,k)
                 END DO
              END DO
           ELSE
              DO j=1,DIM
                 DO k=1,DIM
                    Stress(j,k,i) = Stress(j,k,i) - &
                           Metric(j,k) * SUM( Pressure(1:Pn) * Basis(1:Pn) )

                    DO l=1,DIM
                       Stress(j,k,i) = Stress(j,k,i) - &
                           Metric(j,k) * (2.0d0/3.0d0) * Viscosity * Grad(l,l)
                    END DO
                 END DO
              END DO
           END IF

        END DO

        EdgeLength = EdgeLength + s

        Jump = MATMUL( ( Stress(:,:,1) - Stress(:,:,2)), Normal )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           ResidualNorm = ResidualNorm + s * SUM( Jump(1:DIM) ** 2 )
        ELSE
           CALL InvertMatrix( Metric,3 )
           DO i=1,DIM
              DO j=1,DIM
                 ResidualNorm = ResidualNorm + s*Metric(i,j)*Jump(i)*Jump(j)
              END DO
           END DO
        END IF
     END DO

     Indicator = EdgeLength * ResidualNorm

     DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
     DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)

     DEALLOCATE( NodalViscosity, x, y, z, EdgeBasis, &
           Basis, dBasisdx, Velocity, Pressure )
!------------------------------------------------------------------------------
  END FUNCTION FlowEdgeResidual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute the residual of the Navier-Stokes equation for the bulk elements.
!------------------------------------------------------------------------------
   FUNCTION FlowInsideResidual( Model, Element,  &
          Mesh, Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), FNorm
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,m,n,t,DIM,DOFs

     LOGICAL :: stat, GotIt, Compressible, Convect

     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Density, Viscosity,u, v, w, s, detJ
     REAL(KIND=dp) :: Residual(4), ResidualNorm, Area, ReferencePressure, dt

     REAL(KIND=dp), ALLOCATABLE :: NodalViscosity(:), NodalDensity(:), &
            Basis(:),  dBasisdx(:,:), ddBasisddx(:,:,:)
     REAL(KIND=dp),ALLOCATABLE :: Velocity(:,:), Pressure(:)
     REAL(KIND=dp),ALLOCATABLE :: PrevVelo(:,:), PrevPres(:)
     REAL(KIND=dp),ALLOCATABLE :: Temperature(:), NodalForce(:,:)
     REAL(KIND=dp),ALLOCATABLE :: HeatCapacity(:), ReferenceTemperature(:), &
                      HeatExpansionCoeff(:)

     REAL(KIND=dp) :: SpecificHeatRatio

     REAL(KIND=dp), POINTER :: Gravity(:,:)

     TYPE(ValueList_t), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     FNorm = 0.0d0

     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

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

     DOFs = DIM + 1
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( NodalViscosity(n), NodalDensity(n), Basis(n), dBasisdx(n,3), &
        ddBasisddx(n,3,3), Velocity(3,n), Pressure(n), PrevVelo(3,n),  &
        PrevPres(n), Temperature(n), NodalForce(4,n), HeatCapacity(n), &
        ReferenceTemperature(n), HeatExpansionCoeff(n) )

!
!    Material parameters: density, viscosity, etc.
!    ----------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     NodalDensity(1:n) = ListGetReal( &
         Material, 'Density', n, Element % NodeIndexes, GotIt )

     NodalViscosity(1:n) = ListGetReal( &
         Material, 'Viscosity', n, Element % NodeIndexes, GotIt )

     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
                      minv=1, maxv=Model % NumberOfEquations   )

     Convect = ListGetLogical( Model % Equations(k) % Values, &
                   'NS Convect', GotIt )
     IF ( .NOT. GotIt ) Convect = .TRUE.
!
!    Elementwise nodal solution:
!    ---------------------------
     Velocity = 0.0d0
     DO k=1,DOFs-1
        Velocity(k,1:n) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
     END DO
     Pressure(1:n) = Quant( DOFs*Perm(Element % NodeIndexes) )

!
!    Check for time dep.
!    -------------------
     PrevPres(1:n)     = Pressure(1:n)
     PrevVelo(1:3,1:n) = Velocity(1:3,1:n)

     dt = Model % Solver % dt

     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Flow Solution', .TRUE. )

        PrevVelo = 0.0d0
        DO k=1,DOFs-1
           PrevVelo(k,1:n) = &
              Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes)-DOFs+k,1)
        END DO
        PrevPres(1:n)=Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes),1)
     END IF


!
!    Check for compressible flow equations:
!    --------------------------------------
     Compressible = .FALSE.

     IF (  ListGetString( Material, 'Compressibility Model', GotIt ) == &
                      'perfect gas equation 1' ) THEN

        Compressible = .TRUE.

        Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Temperature(1:n) = &
               Var % Values( Var % Perm(Element % NodeIndexes) )
        ELSE
           Temperature(1:n) = ListGetReal( Material, &
               'Reference Temperature',n,Element % NodeIndexes )
        END IF

        SpecificHeatRatio = ListGetConstReal( Material, &
                  'Specific Heat Ratio' )

        ReferencePressure = ListGetConstReal( Material, &
                   'Reference Pressure' )

        HeatCapacity(1:n) = ListGetReal( Material, &
                      'Heat Capacity',n,Element % NodeIndexes )

        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * HeatCapacity(1:n) * Temperature(1:n) )
     END IF
!
!    Body Forces:
!    ------------
!
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, &
       'Body Force', GotIt, 1, Model % NumberOfBodyForces )

     NodalForce = 0.0d0

     IF ( GotIt .AND. k > 0  ) THEN
!
!       Boussinesq approximation of heat expansion for
!       incompressible flow equations:
!
!       Density for the force term equals to
!
!       \rho = rho_0 (1-\beta(T-T_0)),
!
!       where \beta is the  heat expansion  coefficient,
!       T temperature and \rho_0 and T_0 correspond to
!       stress free state. Otherwise density is assumed
!       constant.
!       ----------------------------------------------
        IF (ListGetLogical(Model % BodyForces(k) % Values,'Boussinesq',GotIt)) THEN

           Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              Temperature(1:n) = &
                  Var % Values( Var % Perm(Element % NodeIndexes) )

              HeatExpansionCoeff(1:n) = ListGetReal( Material, &
                 'Heat Expansion Coefficient',n,Element % NodeIndexes )

              ReferenceTemperature(1:n) = ListGetReal( Material, &
                 'Reference Temperature',n,Element % NodeIndexes )

              Gravity => ListGetConstRealArray( Model % Constants, &
                             'Gravity' )

              k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
                        minv=1, maxv=Model % NumberOfEquations )

              IF ( ListGetLogical( Model % Equations(k) % Values, &
                            'Hydrostatic Pressure', GotIt) ) THEN
                 DO i=1,DIM
                    NodalForce(i,1:n) = ( 1 - HeatExpansionCoeff(1:n) * &
                       ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
                            Gravity(i,1) * Gravity(4,1)
                 END DO
              ELSE
                 DO i=1,DIM
                    NodalForce(i,1:n) = ( -HeatExpansionCoeff(1:n) * &
                       ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
                            Gravity(i,1) * Gravity(4,1)
                 END DO
              END IF
           END IF
        END IF

!
!       Given external force:
!       ---------------------
        NodalForce(1,1:n) = NodalForce(1,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 1', &
                  n, Element % NodeIndexes, GotIt )

        NodalForce(2,1:n) = NodalForce(2,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 2', &
                  n, Element % NodeIndexes, GotIt )

        NodalForce(3,1:n) = NodalForce(3,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 3', &
                  n, Element % NodeIndexes, GotIt )
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
            Basis, dBasisdx, ddBasisddx, .TRUE. )

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

        Density   = SUM( NodalDensity(1:n)   * Basis(1:n) )
        Viscosity = SUM( NodalViscosity(1:n) * Basis(1:n) )
!
!       Residual of the navier-stokes equations:
!
!       or more generally:
!
!       ----------------------------------------------------------
!
        Residual = 0.0d0
        DO i=1,DIM
!
!          given force:
!          -------------
           Residual(i) = -Density * SUM( NodalForce(i,1:n) * Basis(1:n) )
 
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!             + grad(p):
!             ----------
              Residual(i) = Residual(i) + SUM( Pressure(1:n) * dBasisdx(1:n,i) )

              DO j=1,DIM
!
!                - 2 ( \mu \epsilon^{ij} )_{,j}:
!                -------------------------------
                 Residual(i) = Residual(i) - Viscosity * &
                     SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,j) )

                 Residual(i) = Residual(i) - &
                      SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &
                          SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )

                  Residual(i) = Residual(i) - Viscosity * &
                      SUM( Velocity(j,1:n) * ddBasisddx(1:n,i,j) )

                  Residual(i) = Residual(i) - &
                      SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &
                          SUM( Velocity(j,1:n) * dBasisdx(1:n,i) )

                  IF ( Compressible ) THEN
!
!                    + (2/3) grad(\mu div(u)):
!                    -------------------------
                     Residual(i) = Residual(i) + &
                        Viscosity * ( 2.0d0 / 3.0d0 ) * &
                           SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,i) )

                     Residual(i) = Residual(i) + &
                         SUM( NodalViscosity(1:n) * dBasisdx(1:n,i) ) * &
                             SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )

                  END IF
              END DO

              IF ( Convect ) THEN
!
!                + \rho * (@u/@t + u.grad(u)):
!                -----------------------------
                 Residual(i) = Residual(i) + Density *  &
                     SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt

                 DO j=1,DIM
                    Residual(i) = Residual(i) + &
                        Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
                            SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )
                 END DO
              END IF
           ELSE
!             + g^{ij}p_{,j}:
!             ---------------
              DO j=1,DIM
                 Residual(i) = Residual(i) + Metric(i,j) * &
                      SUM( Pressure(1:n) * dBasisdx(1:n,i) )
              END DO

!             - g^{jk} (\mu u^i_{,k})_{,j}):
!             ------------------------------
              DO j=1,DIM
                 DO k=1,DIM
                    Residual(i) = Residual(i) -   &
                         Metric(j,k) * Viscosity * &
                         SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,k) )

                    DO l=1,DIM
                       Residual(i) = Residual(i) +  &
                            Metric(j,k) * Viscosity * Symb(j,k,l) * &
                            SUM( Velocity(i,1:n) * dBasisdx(1:n,l) )

                       Residual(i) = Residual(i) -  &
                            Metric(j,k) * Viscosity * Symb(l,j,i) * &
                            SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )

                       Residual(i) = Residual(i) -  &
                            Metric(j,k) * Viscosity * Symb(l,k,i) * &
                            SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )

                       Residual(i) = Residual(i) -  &
                            Metric(j,k) * Viscosity * dSymb(l,j,i,k) * &
                            SUM( Velocity(l,1:n) * Basis(1:n) )

                       DO m=1,DIM
                          Residual(i) = Residual(i) - Metric(j,k) * Viscosity *&
                                  Symb(m,k,i) * Symb(l,j,m) * &
                                        SUM( Velocity(l,1:n) * Basis(1:n) )

                          Residual(i) = Residual(i) + Metric(j,k) * Viscosity *&
                                  Symb(j,k,m) * Symb(l,m,i) * &
                                        SUM( Velocity(l,1:n) * Basis(1:n) )
                       END DO
                    END DO
                 END DO
              END DO

!             - g^{ik} (\mu u^j_{,k})_{,j}):
!             ------------------------------
              DO j=1,DIM
                 DO k=1,DIM
                    Residual(i) = Residual(i) -   &
                         Metric(i,k) * Viscosity * &
                         SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,k) )

                    DO l=1,DIM
                       Residual(i) = Residual(i) +  &
                            Metric(i,k) * Viscosity * Symb(j,k,l) * &
                            SUM( Velocity(j,1:n) * dBasisdx(1:n,l) )

                       Residual(i) = Residual(i) -  &
                            Metric(i,k) * Viscosity * Symb(l,j,j) * &
                            SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )

                       Residual(i) = Residual(i) -  &
                            Metric(i,k) * Viscosity * Symb(l,k,j) * &
                            SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )

                       Residual(i) = Residual(i) -  &
                            Metric(i,k) * Viscosity * dSymb(l,j,j,k) * &
                            SUM( Velocity(l,1:n) * Basis(1:n) )

                       DO m=1,DIM
                          Residual(i) = Residual(i) - Metric(i,k) * Viscosity *&
                                  Symb(m,k,j) * Symb(l,j,m) * &
                                        SUM( Velocity(l,1:n) * Basis(1:n) )

                          Residual(i) = Residual(i) + Metric(i,k) * Viscosity *&
                                  Symb(j,k,m) * Symb(l,m,j) * &
                                        SUM( Velocity(l,1:n) * Basis(1:n) )
                       END DO
                    END DO
                 END DO
              END DO

              IF ( Convect ) THEN
!
!                + \rho * (@u/@t + u^j u^i_{,j}):
!                --------------------------------
                 Residual(i) = Residual(i) + Density *  &
                     SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt

                 DO j=1,DIM
                    Residual(i) = Residual(i) + &
                         Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
                         SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )

                    DO k=1,DIM
                       Residual(i) = Residual(i) + &
                            Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
                            Symb(j,k,i) * SUM( Velocity(k,1:n) * Basis(1:n) )
                    END DO
                 END DO
              END IF
           END IF
        END DO

!
!       Continuity equation:
!       --------------------
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!
!          + \rho * div(u):
!          ----------------
           DO j=1,DIM
              Residual(DIM+1) = Residual(DIM+1) + &
                   Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
           END DO

           IF ( Compressible ) THEN
!
!             + u.grad(\rho):
!             ----------------
              DO j=1,DIM
                 Residual(DIM+1) = Residual(DIM+1) + &
                      SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
                           SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
              END DO
           END IF
        ELSE
!
!          + \rho * u^j_{,j}:
!          ------------------
           DO j=1,DIM
              Residual(DIM+1) = Residual(DIM+1) + &
                   Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )

              DO k=1,DIM
                 Residual(DIM+1) = Residual(DIM+1) + Density * &
                      Symb(k,j,j) * SUM( Velocity(k,1:n) * Basis(1:n) )
              END DO
           END DO

           IF ( Compressible ) THEN
!
!             + u^j \rho_{,j}:
!             ----------------
              DO j=1,DIM
                 Residual(DIM+1) = Residual(DIM+1) + &
                      SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
                      SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
              END DO
           END IF
        END IF

        DO i=1,DIM
           FNorm = FNorm + s * (Density * SUM(NodalForce(i,1:n)*Basis(1:n))**2)
        END DO 
        Area = Area + s

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           ResidualNorm = ResidualNorm + &
             s * (Element % hK**2 * SUM(Residual(1:dim)**2) + Residual(dim+1)**2 )
        ELSE
           CALL InvertMatrix( Metric,3 )
           DO i=1,dim
              DO j=1,dim
                 ResidualNorm = ResidualNorm + &
                    s * Element % hK **2 * Metric(i,j) * Residual(i) * Residual(j)
              END DO
           END DO
           ResidualNorm = ResidualNorm + s * Residual(dim+1)**2
        END IF
     END DO

!    FNorm = Area * FNorm
     Indicator = ResidualNorm

     DEALLOCATE( NodalViscosity, NodalDensity, Basis, dBasisdx,           &
        ddBasisddx, Velocity, Pressure, PrevVelo, PrevPres, Temperature,   &
        NodalForce, HeatCapacity, ReferenceTemperature, HeatExpansionCoeff,&
        Nodes % x, Nodes % y, Nodes % z )
!------------------------------------------------------------------------------
  END FUNCTION FlowInsideResidual
!------------------------------------------------------------------------------

!> \} 
