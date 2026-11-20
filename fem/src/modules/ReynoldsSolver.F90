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
!/******************************************************************************
! *
! *  Authors: Peter Råback
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 23.10.2007
! *  Modified by: Peter Råback
! *  Modification date: 1.9.2008
! *
! *****************************************************************************/



!------------------------------------------------------------------------------
!> Solves the transient/steady state Reynolds Equation that is a dimensinally
!> reduced form of the Stokes equation in the case of narrow channels.
!> The possible uses are in modeling of lubrication and squeezed-film damping, for example.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ReynoldsSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE SolverBasics

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element, Parent
  TYPE(ValueList_t), POINTER :: Params, Material, Equation, BC, BodyForce

  INTEGER, PARAMETER :: Compressibility_None = 1, Compressibility_Weak = 2, &
      Compressibility_GasIsothermal = 3, Compressibility_GasAdiabatic = 4, &
      Compressibility_Artificial = 5
  INTEGER, PARAMETER :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2

  INTEGER :: iter, i, j, k, l, n, nd, t, istat, mat_id, ent_id, body_id, mat_idold, &
      NoIterations, ViscosityType, CompressibilityType, itime, SensTimestep = -1
  INTEGER, POINTER :: NodeIndexes(:), PressurePerm(:)

  LOGICAL :: GotIt, GotIt2, GotIt3, stat, AllocationsDone = .FALSE., SubroutineVisited = .FALSE., &
      UseVelocity, Bubbles, ApplyLimiter, LinearModel, ManningModel, GotMinGap, &
      OpenSide,GotExt,GotFlux, GotVelo, AnyBC, GotPseudoPressure, Converged
  REAL(KIND=dp), POINTER :: Pressure(:)
  REAL(KIND=dp) :: Norm, ReferencePressure, HeatRatio, BulkModulus, &
      mfp0, Pres, Dens, ManningCoeff, GravityCoeff, MinGap, MinGradPres, &
      ACScale, FSIMult
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:), TimeForce(:), &
      Viscosity(:), GapHeight(:), NormalVelocity(:), Velocity(:,:), &
      Admittance(:), Impedance(:), ElemPressure(:), PrevElemPressure(:),  ExtPressure(:), &
      ElemDensity(:),ElemArtif(:),ExtPres(:),FluxPres(:),VeloPres(:),CoeffPres(:), &
      ElemPseudoPressure(:), PseudoPressure(:)
  TYPE(Variable_t), POINTER :: SensVar, PresVar

  CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityModel, CompressibilityModel, varname
  CHARACTER(*), PARAMETER :: Caller = 'ReynoldsSolver'

  SAVE ElementNodes, Viscosity, GapHeight, ElemArtif, ElemDensity, Velocity, NormalVelocity, &
      Admittance, FORCE, STIFF, MASS, TimeForce, ElemPressure, PrevElemPressure, ExtPressure, &
      AllocationsDone, ExtPres, FluxPres, VeloPres, CoeffPres, PseudoPressure, GotPseudoPressure, &
      ElemPseudoPressure, SensTimestep


  CALL Info(Caller,'---------------------------------------',Level=5)
  IF(TransientSimulation) THEN
    CALL Info(Caller,'Solving the transient Reynolds equation',Level=5)
  ELSE
    CALL Info(Caller,'Solving the steady-state Reynolds equation',Level=5)
  END IF
  CALL Info(Caller,'---------------------------------------',Level=5)

  CALL DefaultStart()

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------

  Params => GetSolverParams()
  Bubbles = GetLogical( Params, 'Bubbles', GotIt )

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF(Solver % Variable % Dofs /= 1) THEN
    CALL Fatal(Caller,'Impossible number of dofs! (should be 1)')
  END IF
  PresVar => Solver % Variable
  Pressure => PresVar % Values
  PressurePerm => PresVar % Perm
  Varname = TRIM(PresVar % Name)
  IF( COUNT( PressurePerm > 0 ) <= 0) RETURN

!------------------------------------------------------------------------------
! Do some initial stuff
!------------------------------------------------------------------------------

  ManningModel = GetLogical( Params,'Manning Model',GotIt)
  IF( ManningModel ) THEN
    GravityCoeff = GetCReal( CurrentModel % Constants,'Gravity Coefficient',GotIt)
    IF(.NOT. GotIt) GravityCoeff = 9.81_dp
  END IF

  AnyBC = ListGetLogicalAnyBC( Model,'Open Side') .OR. &
      ListCheckPresentAnyBC( Model,'Filmpressure Flux') .OR. &
      ListCheckPresentAnyBC( Model,'Filmpressure Velocity') .OR. &
      ListCheckPresentAnyBC( Model,'Filmpressure Transfer Coefficient')

  MinGap = ListGetCReal( Params,'Min Gap Height',GotMinGap)

  MinGradPres = ListGetCReal( Params,'Initial Pressure Gradient',GotIt)
  IF(.NOT. GotIt .OR. AllocationsDone ) THEN
    MinGradPres = EPSILON( MinGradPres )
  END IF

  NoIterations = GetInteger( Params,'Nonlinear System Max Iterations',GotIt)
  IF(.NOT. GotIt) NoIterations = 1


!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------


  IF ( .NOT. AllocationsDone  ) THEN
    n = Solver % Mesh % MaxElementNodes
    ! Local matrices: nodal + legacy bubble dofs, or all dofs of p-elements
    nd = MAX( 2*n, Solver % Mesh % MaxElementDOFs )

    ALLOCATE(ElementNodes % x(n),  &
        ElementNodes % y(n),       &
        ElementNodes % z(n),       &
        Viscosity(n),              &
        GapHeight(n),          &
        ExtPres(n), &
        FluxPres(n), &
        VeloPres(n), &
        CoeffPres(n), &
        ElemArtif(n), &
        ElemDensity(n), &
        Velocity(3,N),         &
        NormalVelocity(n),     &
        Admittance(n),         &
        FORCE( nd ),           &
        STIFF( nd, nd ), &
        MASS( nd, nd ), &
        TimeForce( nd ), &
        ElemPressure(nd), &
        PrevElemPressure(nd), &
        ExtPressure(n), &
        ElemPseudoPressure(n), &
        STAT=istat )
    IF ( istat /= 0 ) CALL FATAL(Caller,'Memory allocation error')

    GotPseudoPressure = .FALSE.
    DO k=1,Model % NumberOfMaterials
      Material => Model % Materials(k) % Values
      CompressibilityModel = ListGetString( Material,'Compressibility Model', GotIt)
      IF (.NOT. GotIt ) CYCLE
      IF( CompressibilityModel == 'artificial compressible') THEN
        GotPseudoPressure = .TRUE.
        ALLOCATE( PseudoPressure(SIZE(Pressure)),STAT=istat)
        EXIT
      END IF
    END DO

    AllocationsDone = .TRUE.

    IF(ListGetLogical( Params,'Skip First Solution', GotIt ) ) THEN
      CALL Info(Caller,'Skipping first solution phase completely!',Level=5)
      RETURN
    END IF
  END IF

  IF(GotPseudoPressure) THEN
    PseudoPressure = Pressure
    ACScale = ListGetConstReal( Model % Simulation, &
        'Artificial Compressibility Scaling',GotIt)
    IF(.NOT.GotIt) ACScale = 1.0
    !IF(Transient) ACScale = ACScale / dt
  END IF

  FsiMult = ListGetCReal( Params,'Fsi Velocity Multiplier',GotIt )
  IF(.NOT. GotIt) FsiMult = 1.0_dp
  itime = GetTimestep()
  IF(itime == 1) THEN
    IF( ListGetLogical( Params,'Steady State Fsi Start', GotIt )  ) THEN
      CALL Info(Caller,'Ignoring FSI velocity for the 1st timestep',Level=7)
      FsiMult = 0.0_dp
    END IF
  END IF


!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------

  mat_idold = 0

  CALL Info(Caller,'-------------------------------------------------',Level=5)

200 DO iter = 1,NoIterations

    LinearModel = ( iter == 1 ) .AND. ListGetLogical( Params,'Linear First Iteration',GotIt)

    WRITE(Message,'(A,T35,I5)') 'Reynolds iteration:',iter
    CALL Info(Caller,Message,Level=5)


100 CONTINUE
    CALL DefaultInitialize()

!    Do the bulk assembly:
!    ---------------------
    CALL GlobalBulkAssembly()
    CALL DefaultFinishBulkAssembly( )

!------------------------------------------------------------------------------
!    Neumann & Newton BCs:
!------------------------------------------------------------------------------
    IF(AnyBC) THEN
      CALL GlobalBoundaryAssemby()
    END IF
    CALL DefaultFinishBoundaryAssembly( )

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()


    ! Check stepsize for nonlinear iteration
    !------------------------------------------------------------------------------
    IF( DefaultLinesearch( Converged ) ) GOTO 100
    IF( Converged ) EXIT

!    Solve the system and we are done:
!    ---------------------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
  END DO

  IF(DefaultSensitivity()) THEN
    NoIterations = 1
    GOTO 200
  END IF


  IF( ListGetLogical( Params,'Gap Sensitivity', GotIt ) ) THEN
    CALL Info(Caller,'Computing FilmPressure sentivity to gap height',Level=5)

    SensVar => VariableGet( Model % Variables,TRIM(Varname)//' Gap Sensitivity')
    IF( .NOT. ASSOCIATED( SensVar ) ) THEN
      CALL Fatal(Caller,'> '//TRIM(Varname)//' gap sensitivity < should exist!')
    END IF

    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.TRUE.)
    ApplyLimiter = ListGetLogical( Params,'Apply Limiter', GotIt )
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .FALSE. )

    IF( ManningModel ) THEN
      CALL Fatal(Caller,'"Gap Sensitivity" is not available with the Manning model!')
    END IF

    ! The time derivative of the sensitivity needs its own history: the
    ! initial pressure does not depend on the gap, hence zero initially.
    IF( TransientSimulation ) THEN
      IF( .NOT. ASSOCIATED( SensVar % PrevValues ) ) THEN
        IF( .NOT. ASSOCIATED( PresVar % PrevValues ) ) THEN
          CALL Fatal(Caller,'No previous values of '//TRIM(Varname)//' for its gap sensitivity!')
        END IF
        ALLOCATE( SensVar % PrevValues( SIZE(SensVar % Values), SIZE(PresVar % PrevValues,2) ) )
        SensVar % PrevValues = 0.0_dp
        SensTimestep = GetTimestep()
      END IF
      IF( GetTimestep() /= SensTimestep ) THEN
        DO k = SIZE(SensVar % PrevValues,2), 2, -1
          SensVar % PrevValues(:,k) = SensVar % PrevValues(:,k-1)
        END DO
        SensVar % PrevValues(:,1) = SensVar % Values
        SensTimestep = GetTimestep()
      END IF
    END IF

    Solver % Variable => SensVar
    SensVar % Values = 0.0_dp

    CALL DefaultInitialize()
    CALL GlobalBulkAssembly( 1 )
    CALL DefaultFinishBulkAssembly( )
    IF(AnyBC) THEN
      CALL GlobalBoundaryAssemby( 1 )
    END IF
    CALL DefaultFinishBoundaryAssembly( )

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs( Ux = SensVar )

!    Solve the system and we are done:
!    ---------------------------------
    Norm = DefaultSolve()

    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.FALSE.)
    Solver % Variable => PresVar
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .TRUE. )
  END IF

  CALL DefaultFinish()

  CALL Info(Caller,'-------------------------------------------------',Level=5)

CONTAINS



  ! Cycle over the bulk elements and build the global linear system
  ! Optionally performs sensitisity analysis.
  !----------------------------------------------------------------------
  SUBROUTINE GlobalBulkAssembly( SensitivityMode )

    INTEGER, OPTIONAL :: SensitivityMode
    INTEGER :: SensMode

    SensMode = 0
    IF( PRESENT(SensitivityMode) ) SensMode = SensitivityMode
    mat_idold = -1


    DO t=1,Solver % NumberOfActiveElements

      Element => GetActiveElement(t)

      IF( Element % TYPE % ElementCode > 500 ) THEN
        CALL Fatal(Caller,'This is a reduced dimensional solver for 1D and 2D only!')
      END IF

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()

      ! The bubbles are condensed as n dofs after the n nodal ones
      IF( Bubbles .AND. nd /= n ) THEN
        CALL Fatal(Caller,'"Bubbles" requires plain nodal elements, got '&
            //I2S(nd)//' dofs for '//I2S(n)//' nodes')
      END IF

      CALL GetElementNodes( ElementNodes )
      CALL GetScalarLocalSolution( ElemPressure, UVariable = PresVar)

      IF( SensMode > 0 ) THEN
        IF( TransientSimulation ) THEN
          CALL GetScalarLocalSolution( PrevElemPressure, tstep = -1, UVariable = PresVar )
        END IF
      END IF

      IF( GotPseudoPressure ) THEN
        ElemPseudoPressure(1:n) = PseudoPressure( PressurePerm(Element % NodeIndexes))
      END IF


      body_id =  Element % Bodyid

      ent_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
      Equation => Model % Equations(ent_id) % Values

      mat_id = GetInteger( Model % Bodies( body_id ) % Values, 'Material')
      Material => Model % Materials(mat_id) % Values

      ent_id = GetInteger( Model % Bodies( body_id ) % Values, 'Body Force',GotIt)
      IF(ent_id>0) THEN
        BodyForce => Model % BodyForces(ent_id) % Values
      ELSE
        BodyForce => NULL()
      END IF

!------------------------------------------------------------------------------
!       Get velocities
!------------------------------------------------------------------------------

      Velocity = 0.0_dp
      UseVelocity = .FALSE.
      GotIt = .FALSE.; GotIt2 = .FALSE.; GotIt3 = .FALSE.
      IF( ListCheckPrefix( BodyForce,'Surface Velocity') ) THEN
        Velocity(1,1:n) = GetReal(BodyForce,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(BodyForce,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(BodyForce,'Surface Velocity 3',GotIt3)
        UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
      END IF
      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( Material,'Surface Velocity') ) THEN
          Velocity(1,1:n) = GetReal(Material,'Surface Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(Material,'Surface Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Material,'Surface Velocity 3',GotIt3)
          UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
        END IF
      END IF

      IF(.NOT. UseVelocity) THEN
        IF( ListCheckPrefix( BodyForce,'Tangent Velocity') ) THEN
          Velocity(1,1:n) = GetReal(BodyForce,'Tangent Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(BodyForce,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(BodyForce,'Tangent Velocity 3',GotIt3)
        END IF
        IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3)) THEN
          IF( ListCheckPrefix( Material,'Tangent Velocity') ) THEN
            Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt)
            Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
            Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
          END IF
        END IF
      END IF

      NormalVelocity(1:n) = GetReal(BodyForce,'Normal Velocity',GotIt)
      IF(.NOT. GotIt) NormalVelocity(1:n) = GetReal(Material,'Normal Velocity',GotIt)

      NormalVelocity(1:n) = NormalVelocity(1:n) + &
          FsiMult * GetReal(BodyForce,'Fsi Velocity',GotIt)
      IF(.NOT. GotIt) NormalVelocity(1:n) = NormalVelocity(1:n) + &
          FsiMult * GetReal(Material,'Fsi Velocity',GotIt)

      IF( ManningModel ) THEN
        ElemDensity(1:n) = GetReal( Material,'Density')
      END IF

!------------------------------------------------------------------------------
!       Get material parameters
!------------------------------------------------------------------------------

      GapHeight(1:n) = GetReal( Material,'Gap Height')
      IF(GotMinGap) GapHeight(1:n) = MAX(GapHeight(1:n),MinGap)

      Admittance(1:n) = GetReal( Material, 'Flow Admittance', GotIt)
      IF(ASSOCIATED(BodyForce)) THEN
        ExtPressure(1:n) = GetReal( BodyForce, 'External FilmPressure', GotIt)
      ELSE
        ExtPressure(1:n) = 0.0_dp
      END IF

      Viscosity(1:n) = GetReal( Material, 'Viscosity')

      CALL GetMaterialModels( Material, mat_id )

      IF( CompressibilityType == Compressibility_Artificial ) THEN
        ElemArtif(1:n) = GetReal( Material,'Artificial Compressibility',GotIt)
      END IF

      STIFF = 0.0_dp
      MASS = 0.0_dp
      FORCE = 0.0_dp

      CALL LocalBulkMatrix( MASS, STIFF, FORCE, Element, n, nd, ElementNodes, SensMode )

!------------------------------------------------------------------------------
!  In time dependent simulation add mass matrix to stiff matrix
!------------------------------------------------------------------------------
      TimeForce  = 0.0_dp
      IF ( TransientSimulation ) THEN
        CALL Default1stOrderTime( MASS,STIFF,FORCE)
      END IF

!------------------------------------------------------------------------------
!  Update global matrices from local matrices
!------------------------------------------------------------------------------
      IF (  Bubbles ) THEN
        CALL Condensate( N, STIFF, FORCE, TimeForce )
      END IF

      CALL DefaultUpdateEquations( STIFF, FORCE )

!------------------------------------------------------------------------------
    END DO

  END SUBROUTINE GlobalBulkAssembly


!------------------------------------------------------------------------------
!> Read the material models and their constant parameters. Cached on the
!> material id; set mat_idold = -1 to force a re-read. Used by both the bulk
!> and the boundary assembly so that a boundary element sees its parent's
!> material, not whichever material the last bulk element happened to have.
!------------------------------------------------------------------------------
  SUBROUTINE GetMaterialModels( Material, mat_id )
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: mat_id

    IF( mat_id == mat_idold ) RETURN
    mat_idold = mat_id

    IF( ManningModel ) THEN
      ManningCoeff = GetCReal(Material,'Manning Coefficient')
    END IF

    ReferencePressure = GetCReal( Material,'Reference Pressure', GotIt )
    ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
    IF(GotIt) THEN
      IF( ViscosityModel == 'newtonian') THEN
        ViscosityType = Viscosity_Newtonian
      ELSE IF( ViscosityModel == 'rarefied') THEN
        ViscosityType = Viscosity_Rarefied
        mfp0 = GetCReal(Material,'Mean Free Path')
      ELSE
        CALL Fatal(Caller,'Unknown viscosity model: '//TRIM(ViscosityModel))
      END IF
    ELSE
      ViscosityType = Viscosity_Newtonian
    END IF

    CompressibilityType = Compressibility_None
    IF( .NOT. LinearModel ) THEN
      CompressibilityModel = GetString(Material,'Compressibility Model',GotIt)
      IF(GotIt) THEN
        IF(CompressibilityModel == 'incompressible') THEN
          CompressibilityType = Compressibility_None
        ELSE IF(CompressibilityModel == 'weakly compressible') THEN
          CompressibilityType = Compressibility_Weak
          BulkModulus = GetCReal( Material, 'Bulk Modulus')
        ELSE IF(CompressibilityModel == 'isothermal ideal gas') THEN
          CompressibilityType = Compressibility_GasIsothermal
        ELSE IF(CompressibilityModel == 'adiabatic ideal gas') THEN
          CompressibilityType = Compressibility_GasAdiabatic
          HeatRatio = GetCReal( Material, 'Specific Heat Ratio')
        ELSE IF( CompressibilityModel == 'artificial compressible') THEN
          CompressibilityType = Compressibility_Artificial
        ELSE
          CALL Fatal(Caller,'Unknown compressibility model: '//TRIM(CompressibilityModel))
        END IF
      END IF
    END IF

  END SUBROUTINE GetMaterialModels
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(MassMatrix, StiffMatrix, ForceVector, &
      Element, n, nd, Nodes, SensMode )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MassMatrix(:,:), StiffMatrix(:,:), ForceVector(:)
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
    INTEGER :: SensMode
!------------------------------------------------------------------------------
    ! ElementInfo evaluates the bubbles only if Basis has room for 2*n
    REAL(KIND=dp) :: Basis(MAX(2*n,nd)),dBasisdx(MAX(2*n,nd),3), detJ
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S, MS, MM, MA, L, A, B, HR, SL(3), SLR, SLL(3), F, AS
    REAL(KIND=dp) :: Normal(3), Velo(3), NormalVelo, TangentVelo(3), Damp, Pres, PrevPres, &
        TotPres, GradPres(3), AbsGradPres, dPdt, Gap, Visc, mfp, Kn, Density, DensityDer, &
        PseudoPres, ExtPres, DensityDer2, RarFac, JC, dMSdp, dMSdh, dSLLdp(3), FS, FG, &
        DivHU, DivU, GradP_Test
    LOGICAL :: Stat, GotAC
    INTEGER :: i,p,q,t,DIM, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()
    GotAC = .FALSE.

    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = n
    IF ( Bubbles ) THEN
      NBasis = 2*n
      IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
    ELSE
      NBasis = nd
      IntegStuff = GaussPoints( Element )
    END IF

!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, u, v, w, DetJ, &
               Basis, dBasisdx, Bubbles = Bubbles)

      s = s * DetJ
      IF ( CoordSys /= Cartesian ) THEN
        x = SUM( Nodes % x(1:n) * Basis(1:n) )
        y = SUM( Nodes % y(1:n) * Basis(1:n) )
        z = SUM( Nodes % z(1:n) * Basis(1:n) )
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF

!------------------------------------------------------------------------------
!      Parameters at integration point
!------------------------------------------------------------------------------

      NormalVelo = SUM( Basis(1:n) * NormalVelocity(1:n) )
      IF(UseVelocity) THEN
        IF( ASSOCIATED( Element % BoundaryInfo)) THEN
          Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
        ELSE
          Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
        END IF
        DO i=1,3
          Velo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
        NormalVelo = NormalVelo + SUM( Normal * Velo)
        TangentVelo = Velo - NormalVelo * Normal
      ELSE
        DO i=1,3
          TangentVelo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
      END IF

      Damp = SUM(Basis(1:n) * Admittance(1:n))
      ExtPres = SUM(Basis(1:n) * ExtPressure(1:n))
      Pres = SUM(Basis(1:nd) * ElemPressure(1:nd))
      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      TotPres = ReferencePressure + Pres

      ! If we compute sensitivity of solution we need various derivatives of pressure
      !--------------------------------------------------------------------------------
      IF( SensMode > 0 ) THEN
        IF( TransientSimulation ) THEN
          PrevPres = SUM( Basis(1:nd) * PrevElemPressure(1:nd) )
          dPdt = ( Pres - PrevPres ) / dt
        ELSE
          dPdt = 0.0_dp
        END IF
        DO i = 1,3
          GradPres(i) = SUM( dBasisdx(1:nd,i) * ElemPressure(1:nd) )
        END DO
      END IF


!------------------------------------------------------------------------------
!  Different material models. The "density" is used only as a functional form,
!  not as absolute value.
!------------------------------------------------------------------------------

      ! RarFac = -dlog(Visc)/dlog(Kn), Kn being proportional to 1/(TotPres*Gap)
      RarFac = 0.0_dp
      SELECT CASE (ViscosityType)

      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))

      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)
        RarFac = 1.159_dp * 9.638_dp*Kn**1.159_dp / (1+9.638_dp*Kn**1.159_dp)

      END SELECT

!------------------------------------------------------------------------------

      DensityDer2 = 0.0_dp
      SELECT CASE (CompressibilityType)

      CASE (Compressibility_None)
        Density = 1.0d0
        DensityDer = 0.0d0

      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
        DensityDer = Density / BulkModulus
        DensityDer2 = DensityDer / BulkModulus

      CASE (Compressibility_GasIsothermal)
        Density = TotPres
        DensityDer = 1.0d0

      CASE(Compressibility_GasAdiabatic)
        Density = TotPres ** (1.0_dp/HeatRatio)
        DensityDer = (1/HeatRatio) * TotPres ** (1.0_dp/HeatRatio - 1.0_dp)
        DensityDer2 = (1/HeatRatio) * (1/HeatRatio - 1.0_dp) * &
            TotPres ** (1.0_dp/HeatRatio - 2.0_dp)

      CASE (Compressibility_Artificial )
        Density = 1.0d0
        DensityDer = 0.0d0
        GotAC = .TRUE.

      END SELECT

!------------------------------------------------------------------------------
!  Coefficients of the differential equation at integration point
!------------------------------------------------------------------------------

      ! Multipliers of p: Stiffness matrix
      IF( ManningModel ) THEN
        Density = Density * SUM( Basis(1:n) * ElemDensity(1:n) )
        DensityDer = DensityDer * SUM( Basis(1:n) * ElemDensity(1:n) )
        DO i = 1,3
          GradPres(i) = SUM( dBasisdx(1:nd,i) * ElemPressure(1:nd) )
        END DO
        AbsGradPres = SQRT( SUM( GradPres**2 ) )
        AbsGradPres = MAX( AbsGradPres, MinGradPres )
        MS = -SQRT(Density/(GravityCoeff*AbsGradPres)) * Gap**(5.0/3)  / (2**(2.0/3) * ManningCoeff)
      ELSE
        MS = -Density * Gap**3 / (12 * Visc)
      END IF

      IF(CompressibilityType == Compressibility_GasIsothermal) THEN
        HR = -Damp * ( Density + ExtPres ) / 2
      ELSE IF(CompressibilityType == Compressibility_GasAdiabatic) THEN
        HR = -Damp * ( Density  + ExtPres ** (1.0_dp/HeatRatio) ) / 2
      ELSE
        HR = -Damp * Density
      END IF

      ! Multipliers of dp/dt: Mass matrix
      MM = -DensityDer * Gap

      ! Multiplier of dp/dt in terms of artificial copressibility
      ! This is pseudotime, not real time...
      MA = 0.0_dp
      IF( GotAC ) THEN
        PseudoPres = SUM(Basis(1:n) * ElemPseudoPressure(1:n) )
        MA = ( -Density / dt ) * SUM( ElemArtif(1:n) * Basis(1:n) )
        MA = ACScale * MA
      END IF

      ! Normal velocity: right-hand-side force vector
      L = Density * NormalVelo

      ! Tangential velocity: Both rhs and matrix contribution
      SLR = 0.0_dp
      SLL = 0.0_dp

      DO i=1,dim
        ! The plane element automatically omits the derivative in normal direction
        SLR = SLR + 0.5_dp * Density * SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n))
      END DO
      ! Implicit part: coefficient of the pressure gradient
      SLL = -0.5_dp* DensityDer * Gap * TangentVelo

      ! Gap sensitivity S = dp/dh for a uniform change h of the gap: with R(p,h)
      ! the residual assembled for SensMode == 0, solve (dR/dp) S = -dR/dh.
      ! Here dR/dp adds to the Picard matrix the derivatives of the
      ! coefficients (JC: of S*phi_q*phi_p, dMSdp: of S*grad(p).grad(phi_p),
      ! dSLLdp: of S*grad(p)*phi_p). The r.h.s. is FS*phi_p + FG*grad(p).grad(phi_p).
      ! The leakage (Flow Admittance) and the normal velocity do not depend on h.
      !------------------------------------------------------------------------------
      IF( SensMode == 1 ) THEN
        DivHU = 0.0_dp
        DivU = 0.0_dp
        DO i=1,dim
          DivHU = DivHU + SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n) )
          DivU = DivU + SUM( dBasisdx(1:n,i) * Velocity(i,1:n) )
        END DO

        ! MS = -Density Gap**3 / (12 Visc), Visc = Visc(Kn)
        dMSdp = -DensityDer * Gap**3 / (12 * Visc) - MS * RarFac / TotPres
        dMSdh = 3 * MS / Gap - MS * RarFac / Gap

        ! HR*(p-ExtPres): HR depends on p through the density
        IF( CompressibilityType == Compressibility_GasIsothermal .OR. &
            CompressibilityType == Compressibility_GasAdiabatic ) THEN
          JC = -Damp * DensityDer / 2
        ELSE
          JC = -Damp * DensityDer
        END IF
        JC = JC * ( Pres - ExtPres )

        ! r.h.s. L = Density*NormalVelo and SLR = Density*DivHU/2
        JC = JC - DensityDer * NormalVelo - 0.5_dp * DensityDer * DivHU
        dSLLdp = -0.5_dp * DensityDer2 * Gap * TangentVelo

        ! Mass term MM*dp/dt, MM = -DensityDer*Gap
        JC = JC - DensityDer2 * Gap * dPdt

        FG = -dMSdh
        FS = 0.5_dp * DensityDer * SUM( TangentVelo * GradPres ) &
            + 0.5_dp * Density * DivU + DensityDer * dPdt
      END IF

!------------------------------------------------------------------------------
!      The Reynolds equation
!------------------------------------------------------------------------------
      DO p=1,NBasis
        F = 0.0_dp

        DO q=1,NBasis
          A = (MA + HR) * Basis(q) * Basis(p)
          DO i=1,DIM
            DO j=1,DIM
              A = A + MS * Metric(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
              A = A + SLL(j) * Metric(i,j) * dBasisdx(q,i) * Basis(p)
            END DO
          END DO

          IF( SensMode == 1 ) THEN
            A = A + JC * Basis(q) * Basis(p)
            DO i=1,DIM
              DO j=1,DIM
                A = A + dMSdp * Metric(i,j) * GradPres(i) * dBasisdx(p,j) * Basis(q)
                A = A + dSLLdp(j) * Metric(i,j) * GradPres(i) * Basis(q) * Basis(p)
              END DO
            END DO
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A

          IF( TransientSimulation ) THEN
            B = MM * Basis(q) * Basis(p)
            MassMatrix(p,q)  = MassMatrix(p,q)  + s * B
          END IF
        END DO

        IF( SensMode == 0 ) THEN
          F = L + SLR + HR * ExtPres
          IF(GotAC) F = F + MA * PseudoPres
        ELSE IF( SensMode == 1 ) THEN
          F = FS
          GradP_Test = 0.0_dp
          DO i=1,DIM
            DO j=1,DIM
              GradP_Test = GradP_Test + Metric(i,j) * GradPres(i) * dBasisdx(p,j)
            END DO
          END DO
          ForceVector(p) = ForceVector(p) + s * FG * GradP_Test
        END IF

        ForceVector(p) = ForceVector(p) + s * Basis(p) * F

      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalBulkMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Cycle over boundary elements and add the flux BCs.
!> Currently only such BC is the condition for open side.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalBoundaryAssemby( SensitivityMode )

    INTEGER, OPTIONAL :: SensitivityMode
    INTEGER :: SensMode

    SensMode = 0
    IF( PRESENT(SensitivityMode) ) SensMode = SensitivityMode
    mat_idold = -1
    DO t=1, Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE

      BC => GetBC()
      IF ( .NOT. ASSOCIATED( BC ) ) CYCLE

      n  = GetElementNOFNodes()

      OpenSide = GetLogical(BC,'Open Side',gotIt)
      FluxPres(1:n) = GetReal(BC,'Filmpressure Flux',GotFlux)
      VeloPres(1:n) = GetReal(BC,'Filmpressure Velocity',GotVelo)
      CoeffPres(1:n) = GetReal(BC,'Filmpressure Transfer Coefficient',GotIt)
      ExtPres(1:n) = GetReal(BC,'External FilmPressure',GotExt)
      IF(GotExt .NEQV. GotIt) THEN
        CALL Fatal(Caller,'Give neither or both keywords for Robin BC!')
      END IF

      IF(.NOT. (OpenSide .OR. GotExt .OR. GotFlux .OR. GotVelo ) ) CYCLE

!------------------------------------------------------------------------------
      nd = GetElementNOFDOFs()
      IF ( GetElementFamily() == 1 ) CYCLE
      NodeIndexes => Element % NodeIndexes

      IF ( ANY( PressurePerm(NodeIndexes(1:n)) == 0 ) ) CYCLE

      Parent => Element % BoundaryInfo % Left
      stat = ASSOCIATED( Parent )
      IF ( stat ) stat = ALL(PressurePerm(Parent % NodeIndexes) > 0)

      IF(.NOT. stat) THEN
        Parent => ELement % BoundaryInfo % Right
        stat = ASSOCIATED( Parent )
        IF ( stat ) stat = ALL(PressurePerm(Parent % NodeIndexes) > 0)
        IF ( .NOT. stat )  CALL Fatal( Caller, &
            'No proper parent element available for specified boundary' )
      END IF

      CALL GetElementNodes( ElementNodes )
      CALL GetScalarLocalSolution( ElemPressure, UVariable = PresVar )

      mat_id = GetInteger( Model % Bodies(Parent % BodyId) % Values,'Material')
      Material => Model % Materials(mat_id) % Values
      CALL GetMaterialModels( Material, mat_id )

      GapHeight(1:n) = GetReal(Material,'Gap Height')
      IF(GotMinGap) GapHeight(1:n) = MAX( GapHeight(1:n), MinGap )

      Viscosity(1:n) = GetReal( Material, 'Viscosity')

      IF( ManningModel ) THEN
        ElemDensity(1:n) = GetReal( Material,'Density')
      END IF

      STIFF = 0.0d0
      MASS = 0.0d0
      FORCE = 0.0d0

!------------------------------------------------------------------------------
!             Get element local matrix and rhs vector
!------------------------------------------------------------------------------
      CALL LocalBoundaryMatrix( MASS, STIFF, FORCE, Element, n, nd, ElementNodes, SensMode )

!------------------------------------------------------------------------------
!             Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
      IF ( TransientSimulation ) THEN
        MASS = 0.d0
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF

      CALL DefaultUpdateEquations( STIFF, FORCE )

 !------------------------------------------------------------------------------
    END DO


  END SUBROUTINE GlobalBoundaryAssemby
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix(MASS, STIFF, FORCE, Element, n, nd, Nodes, SensMode)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MASS(:,:), STIFF(:,:), FORCE(:)
    INTEGER :: n, nd, SensMode
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: DetJ,U,V,W,S
    REAL(KIND=dp) :: Basis(nd)
    REAL(KIND=dp) :: Visc, dl, mfp, Kn, Damp, TotPres, Pres, Density, Gap, A, B, &
        DensityDer, KnFac, AOpen, dAdp, dAdh
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,DIM,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n
      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, DetJ, Basis )

      s = s * DetJ

      Gap = SUM( GapHeight(1:n) * Basis(1:n) )
      Pres = SUM( ElemPressure(1:nd) * Basis(1:nd))
      TotPres = ReferencePressure + Pres

      ! KnFac = dlog(1/(Visc*dl))/dlog(Kn), Kn being proportional to 1/(TotPres*Gap)
      KnFac = 0.0_dp
      SELECT CASE (ViscosityType)

      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        dl = 0.8488_dp

      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)
        dl = 0.8488_dp*(1+2.676_dp*Kn**0.659_dp)
        KnFac = 1.159_dp * 9.638_dp*Kn**1.159_dp / (1+9.638_dp*Kn**1.159_dp) &
            - 0.659_dp * 2.676_dp*Kn**0.659_dp / (1+2.676_dp*Kn**0.659_dp)
      END SELECT

      Damp = Gap / (12 * dl * Visc)

!------------------------------------------------------------------------------

      DensityDer = 0.0_dp
      SELECT CASE (CompressibilityType)

      CASE (Compressibility_None)
        Density = 1.0d0

      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
        DensityDer = Density / BulkModulus

      CASE (Compressibility_GasIsothermal)
        Density = TotPres
        DensityDer = 1.0_dp

      CASE(Compressibility_GasAdiabatic)
        Density = TotPres ** (1.0_dp/HeatRatio)
        DensityDer = (1/HeatRatio) * TotPres ** (1.0_dp/HeatRatio - 1.0_dp)

      CASE (Compressibility_Artificial)
        Density = 1.0_dp

      END SELECT

      IF(ManningModel) THEN
        Density = Density * SUM( Basis(1:n) * ElemDensity(1:n) )
      END IF


!------------------------------------------------------------------------------
      A = 0.0_dp
      B = 0.0_dp
      IF( GotExt ) THEN
        A = -SUM(Basis(1:n) * CoeffPres(1:n) )
        B = A * SUM( Basis(1:n) * ExtPres(1:n) )
      END IF
      IF( GotFlux ) THEN
        B = B - SUM( Basis(1:n) * FluxPres(1:n) )
      END IF
      IF( GotVelo ) THEN
        B = B - Gap * SUM( Basis(1:n) * VeloPres(1:n) )
      END IF
      IF(OpenSide) THEN
        AOpen = -Damp * Gap * Density
        A = A + AOpen
      END IF

      ! Gap sensitivity: dR/dp S = -dR/dh, R = A(p,h) p - B(h) as above. The
      ! Robin and flux terms do not depend on the gap.
      IF( SensMode == 1 ) THEN
        B = 0.0_dp
        IF( GotVelo ) THEN
          B = B - SUM( Basis(1:n) * VeloPres(1:n) )
        END IF
        IF( OpenSide ) THEN
          ! AOpen = -Density * Gap**2 / (12 dl Visc), dl and Visc depending on Kn
          dAdp = AOpen * ( DensityDer / Density - KnFac / TotPres )
          dAdh = AOpen * ( 2.0_dp - KnFac ) / Gap
          A = A + dAdp * Pres
          B = B - dAdh * Pres
        END IF
      END IF

!------------------------------------------------------------------------------
      DO p=1,nd
        FORCE(p) = FORCE(p) + s * Basis(p) * B
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + s * Basis(q) * Basis(p) * A
        END DO
      END DO
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalBoundaryMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE ReynoldsSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver, i.e. ReynoldsSolver.
!------------------------------------------------------------------------------

SUBROUTINE ReynoldsSolver_init( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: Params, BC
  CHARACTER(LEN=MAX_NAME_LEN) :: Varname
  INTEGER :: i

  Params => GetSolverParams()

  VarName = ListGetString( Params,'Variable',Found)
  IF(.NOT. Found ) THEN
    Varname = 'FilmPressure'
    CALL ListAddString( Params, 'Variable', VarName )
  END IF

! The new way with generic limiters is a library functionality.
! The Poisson equation is assembled using different sign that the
! typical convention. Hence the load sign for limiters is opposite
!----------------------------------------------------------------
  CALL ListAddLogical( Params,'Limiter Load Sign Negative',.TRUE.)

  IF( ListGetLogical( Params,'Gap Sensitivity', Found ) ) THEN
    CALL ListAddStrinG( Params,NextFreeKeyword('Exported Variable',Params),&
        TRIM(VarName)//' Gap Sensitivity')

    ! Where the pressure is given it does not depend on the gap
    DO i=1,Model % NumberOfBCs
      BC => Model % BCs(i) % Values
      IF( .NOT. ListCheckPresent( BC, VarName ) ) CYCLE
      IF( ListCheckPresent( BC, TRIM(VarName)//' Gap Sensitivity' ) ) CYCLE
      CALL ListAddConstReal( BC, TRIM(VarName)//' Gap Sensitivity', 0.0_dp )
    END DO
  END IF

  IF( ListGetLogical( Params,'Calculate Sensitivity', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Filmpressure Sensitivity' )
    CALL ListAddString(Params,'Sensitivity Variable','Filmpressure Sensitivity')
  END IF

END SUBROUTINE ReynoldsSolver_init




!------------------------------------------------------------------------------
!> Solver various postprocessing fields from the solution of the Reynolds
!> equation. These include heating, forces, and flux.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ReynoldsPostprocess( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription
  USE SolverBasics

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: PressureVar, VarResult, SolverVar, SensVar
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params, Material, Equation, BodyForce

  INTEGER, PARAMETER :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2

  INTEGER :: iter, i, j, k, l, n, nd, t, istat, ent_id, body_id, mat_id, mat_idold,&
      ViscosityType, dim
  INTEGER :: Component, Components, Mode
  INTEGER, POINTER :: NodeIndexes(:), PressurePerm(:)
  REAL(KIND=dp), POINTER :: mWork(:,:)

  LOGICAL :: GotIt, GotIt2, GotIt3, stat, UseVelocity, AllocationsDone = .FALSE., &
      OpposingWall, CalculateMoment, ManningModel, GotMinGap

  REAL(KIND=dp), POINTER :: Pressure(:), gWork(:,:)
  REAL(KIND=dp) :: Norm, ReferencePressure, mfp0, HeatSlide, HeatPres, HeatTotal, &
      Pforce(3), Vforce(3), Sforce(3), TotForce, Moment(3), MomentAbout(3), AmbientPres, &
      ManningCoeff, GravityCoeff, MinGap
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), Viscosity(:), GapHeight(:), &
      Velocity(:,:), ElemPressure(:), SensPressure(:), BotHeight(:), ElemDensity(:)
  ! One r.h.s. for each component of each field, all with the same mass matrix
  INTEGER :: nrhs, np, RhsMode(16), RhsComp(16)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: ForceVector(:,:)
  REAL(KIND=dp), ALLOCATABLE :: LocalForce(:,:), SavedValues(:)
  REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
  LOGICAL :: GotFact, Refactorize, WarnNodal
  TYPE(Variable_t), POINTER :: CorrVar
  TYPE(Solver_t), POINTER :: PresSolver
  CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityModel, PressureName
  CHARACTER(*), PARAMETER :: Caller = 'ReynoldsPostprocess'


  SAVE ElementNodes, Viscosity, Velocity, GapHeight, ElemDensity, BotHeight, &
      FORCE, STIFF, ElemPressure, SensPressure, AllocationsDone


!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------

  CALL Info(Caller,'--------------------------------------------------',Level=5)
  CALL Info(Caller,'Computing postprocessing fields from film pressure',Level=5)
  CALL Info(Caller,'--------------------------------------------------',Level=5)

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------

  Params => GetSolverParams()

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF(Solver % Variable % Dofs /= 1) THEN
    CALL Fatal(Caller,'Impossible number of dofs! (should be 1)')
  END IF

  ManningModel = GetLogical(Params,'Manning Model',GotIt)

  PressureName = GetString(Params,'Reynolds Pressure Variable Name',GotIt)
  IF(.NOT. GotIt) PressureName = 'FilmPressure'

  PressureVar => VariableGet( Solver % Mesh % Variables, PressureName)
  IF(.NOT. ASSOCIATED(PressureVar)) THEN
    CALL Info(Caller,'Give pressure variable name with: "Reynolds Pressure Variable Name"',Level=3)
    CALL Fatal(Caller,'Could not find primary variable: '//TRIM(PressureName))
  END IF

  SensVar => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Gap Sensitivity')
  IF(ASSOCIATED(SensVar)) THEN
    CALL Info(Caller,'We got pressure sensitivity available!')
  END IF

  IF( ManningModel ) THEN
    gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
    IF ( GotIt ) THEN
      GravityCoeff = gWork(4,1)
    ELSE
      GravityCoeff = GetCReal( CurrentModel % Constants,'Gravity Coefficient',GotIt)
    END IF
    IF(.NOT. GotIt) GravityCoeff = 9.81_dp
  END IF

  MinGap = ListGetCReal( Params,'Min Gap Height',GotMinGap)

  AmbientPres = ListGetCReal( Params,'Ambient Pressure',GotIt)

  Pressure => PressureVar % Values
  PressurePerm => PressureVar % Perm
  ! No early return here even if this partition has no pressure dofs:
  ! this solver may still be active here, and the DefaultSolve and
  ! ParallelReduction calls below are collective over the active
  ! partitions. Elements without pressure are skipped in the loop.
  IF( COUNT( PressurePerm > 0 ) <= 0) THEN
    CALL Info(Caller,'No '//TRIM(PressureName)//' dofs in this partition',Level=8)
  END IF

  DIM = CoordinateSystemDimension()
  ReferencePressure = 0.0_dp

!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone  ) THEN
    N = Solver % Mesh % MaxElementNodes
    nd = MAX( n, Solver % Mesh % MaxElementDOFs )

    ALLOCATE(ElementNodes % x(n),  &
        ElementNodes % y(n),       &
        ElementNodes % z(n),       &
        Viscosity(n),              &
        GapHeight(n),          &
        ElemDensity(n), &
        BotHeight(n), &
        Velocity(3,n),         &
        FORCE(nd),           &
        STIFF(nd,nd), &
        ElemPressure(nd), &
        SensPressure(nd), &
        STAT=istat )

    IF ( istat /= 0 ) CALL FATAL(Caller,'Memory allocation error')

    AllocationsDone = .TRUE.
  END IF

  MomentAbout = 0.0_dp
  mWork => ListGetConstRealArray( Params,'Moment About',CalculateMoment)
  IF( CalculateMoment ) THEN
    i = MIN( dim, SIZE(mWork,1) )
    MomentAbout(1:i) = mWork(1:i,1)
  ELSE
    MomentAbout(1) = ListGetCReal( Params,'Moment About 1', CalculateMoment )
    MomentAbout(2) = ListGetCReal( Params,'Moment About 2', stat )
    CalculateMoment = stat .OR. CalculateMoment
    MomentAbout(3) = ListGetCReal( Params,'Moment About 3', stat )
    CalculateMoment = stat .OR. CalculateMoment
  END IF


!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------

  mat_idold = 0
  HeatSlide = 0.0_dp
  HeatPres = 0.0_dp
  Pforce = 0.0_dp
  Sforce = 0.0_dp
  Vforce = 0.0_dp
  Moment = 0.0_dp

  SolverVar => Solver % Variable
  IF( .NOT. ASSOCIATED( SolverVar ) ) THEN
    CALL Fatal(Caller,'Solver Variable not associated')
  END IF

  CALL Info(Caller,'Primary variable name: '//TRIM( Solver % variable % Name) )

  CorrVar => NULL()
  IF( ManningModel ) THEN
    CorrVar => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Corrected')
  END IF

  ! All the fields are L2 projections with the same mass matrix: collect
  ! one r.h.s. vector for each component of each requested field. The
  ! heating, if requested, is the solver variable itself and comes last.
  !-------------------------------------------------------------------
  nrhs = 0
  DO Mode = 1, 4
    VarResult => ModeVariable( Mode )
    IF(.NOT. ASSOCIATED(VarResult)) CYCLE
    DO Component = 1, VarResult % Dofs
      nrhs = nrhs + 1
      IF( nrhs > SIZE(RhsMode) ) CALL Fatal(Caller,'Too many fields to compute!')
      RhsMode(nrhs) = Mode
      RhsComp(nrhs) = Component
    END DO
  END DO

  IF( nrhs > 0 ) THEN
    SaveRHS => Solver % Matrix % RHS
    ALLOCATE( ForceVector( SIZE(SaveRHS), nrhs ), LocalForce( SIZE(FORCE), nrhs ) )
    ForceVector = 0.0_dp
    CALL DefaultInitialize()
  END IF

  ! The pressure solver may use different element dofs than this one:
  ! then only the nodal values (the vertex functions) are used.
  PresSolver => PressureVar % Solver
  IF(.NOT. ASSOCIATED(PresSolver)) PresSolver => Solver
  WarnNodal = .TRUE.

  DO t=1,Solver % NumberOfActiveElements

    Element => GetActiveElement(t)
    n  = GetElementNOFNodes()
    nd = GetElementNOFDOFs()

    NodeIndexes => Element % NodeIndexes
    ! The postprocessor may be active where the pressure is not
    IF( ANY( PressurePerm(NodeIndexes(1:n)) <= 0 ) ) CYCLE

    CALL GetElementNodes( ElementNodes )

    np = GetElementNOFDOFs( Element, PresSolver )
    IF( np /= nd ) THEN
      IF( WarnNodal .AND. np > n ) THEN
        CALL Warn(Caller,'Element dofs differ from those of '//TRIM(PressureName)//&
            ', using its nodal values only')
        WarnNodal = .FALSE.
      END IF
      np = n
    END IF
    CALL GetScalarLocalSolution( ElemPressure, UVariable = PressureVar )
    IF(ASSOCIATED(SensVar)) THEN
      CALL GetScalarLocalSolution( SensPressure, UVariable = SensVar )
    END IF

    body_id =  Element % Bodyid
    IF( body_id <= 0 ) CYCLE

    ent_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
    Equation => Model % Equations(ent_id) % Values

    mat_id = GetInteger( Model % Bodies( body_id ) % Values, 'Material')
    Material => Model % Materials(mat_id) % Values

    ent_id = GetInteger( Model % Bodies( body_id ) % Values, 'Body Force',GotIt)
    IF(ent_id>0) THEN
      BodyForce => Model % BodyForces(ent_id) % Values
    ELSE
      BodyForce => NULL()
    END IF

    !------------------------------------------------------------------------------
    !       Get velocities
    !------------------------------------------------------------------------------
    Velocity = 0.0_dp
    UseVelocity = .FALSE.
    GotIt = .FALSE.; GotIt2 = .FALSE.; GotIt3 = .FALSE.

    IF( ListCheckPrefix( BodyForce,'Surface Velocity') ) THEN
      Velocity(1,1:n) = GetReal(BodyForce,'Surface Velocity 1',GotIt)
      Velocity(2,1:n) = GetReal(BodyForce,'Surface Velocity 2',GotIt2)
      Velocity(3,1:n) = GetReal(BodyForce,'Surface Velocity 3',GotIt3)
      UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
    END IF
    IF(.NOT. UseVelocity) THEN
      IF( ListCheckPrefix( Material,'Surface Velocity') ) THEN
        Velocity(1,1:n) = GetReal(Material,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(Material,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(Material,'Surface Velocity 3',GotIt3)
        UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
      END IF
    END IF

    IF(.NOT. UseVelocity) THEN
      IF( ListCheckPrefix( BodyForce,'Tangent Velocity') ) THEN
        Velocity(1,1:n) = GetReal(BodyForce,'Tangent Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(BodyForce,'Tangent Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(BodyForce,'Tangent Velocity 3',GotIt3)
      END IF
      IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3 )) THEN
        IF( ListCheckPrefix( Material,'Tangent Velocity') ) THEN
          Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
        END IF
      END IF
    END IF

    !------------------------------------------------------------------------------
    !       Get material parameters
    !------------------------------------------------------------------------------
    GapHeight(1:n) = GetReal( Material,'Gap Height')
    IF(GotMinGap) GapHeight(1:n) = MAX(GapHeight(1:n), MinGap)

    Viscosity(1:n) = GetReal( Material, 'Viscosity')

    IF(mat_id /= mat_idold) THEN

      mat_idold = mat_id

      ReferencePressure = GetCReal( Material,'Reference Pressure', GotIt )
      ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
      IF(GotIt) THEN
        IF( ViscosityModel == 'newtonian') THEN
          ViscosityType = Viscosity_Newtonian
        ELSE IF( ViscosityModel == 'rarefied') THEN
          ViscosityType = Viscosity_Rarefied
          mfp0 = GetCReal(Material,'Mean Free Path')
        ELSE
          CALL Fatal(Caller,'Unknown viscosity model: '//TRIM(ViscosityModel))
        END IF
      ELSE
        ViscosityType = Viscosity_Newtonian
      END IF

      OpposingWall = ListGetLogical( Material,'Opposing Wall',GotIt)
    END IF

    IF( ManningModel ) THEN
      ManningCoeff = ListGetCReal( Material,'Manning Coefficient')
    END IF

    ! If we are solving the equation with the Manning's model we are actually solving for hydraulic pressure
    ! and the physical pressure is obtained as a postprocessing step. Now FEM equation needed it is just
    ! simple subtraction.
    IF( ASSOCIATED( CorrVar ) ) THEN
      ElemDensity(1:n) = GetReal( Material,'Density')
      BotHeight(1:n) = GetReal( Material,'Bedrock Elevation')
      CorrVar % Values(CorrVar % Perm(NodeIndexes)) = Pressure(PressurePerm(NodeIndexes)) &
          - GravityCoeff * ElemDensity(1:n) * ( BotHeight(1:n) + GapHeight(1:n) )
    END IF

    IF( nrhs == 0 ) CYCLE

    STIFF = 0.0_dp
    LocalForce = 0.0_dp

    CALL LocalMatrix( STIFF, LocalForce, Element, n, nd, np, ElementNodes )

    Solver % Matrix % RHS => SaveRHS
    CALL DefaultUpdateEquations( STIFF, LocalForce(:,1) )
    DO k=1,nrhs
      Solver % Matrix % RHS => ForceVector(:,k)
      CALL DefaultUpdateForce( LocalForce(:,k) )
    END DO
  END DO

  ! Solve the projections: the matrix (and with it the factorization or
  ! the preconditioner) is the same for all, only the r.h.s. differs.
  !-------------------------------------------------------------------
  IF( nrhs > 0 ) THEN
    ALLOCATE( SavedValues( SIZE(Solver % Matrix % Values) ) )
    SavedValues = Solver % Matrix % Values

    GotFact = ListCheckPresent( Params,'Linear System Refactorize' )
    IF( GotFact ) Refactorize = ListGetLogical( Params,'Linear System Refactorize',GotIt )

    DO k=1,nrhs
      Mode = RhsMode(k)
      Component = RhsComp(k)

      IF( k > 1 ) THEN
        ! What DefaultInitialize and the assembly would give again
        Solver % Matrix % Values = SavedValues
        IF( ALLOCATED(Solver % Matrix % ConstrainedDOF) ) Solver % Matrix % ConstrainedDOF = .FALSE.
        IF( ALLOCATED(Solver % Matrix % Dvalues) ) Solver % Matrix % Dvalues = 0.0_dp
        CALL ListAddLogical( Params,'Linear System Refactorize',.FALSE. )
      END IF
      Solver % Matrix % RHS => ForceVector(:,k)

      CALL DefaultFinishAssembly()

      ! There could be some periodic BCs hence the BCs
      !-----------------------------------------------
      CALL DefaultDirichletBCs()

      Norm = DefaultSolve()

      ! All but heating are computed one component at a time
      IF( Mode /= 4 ) THEN
        VarResult => ModeVariable( Mode )
        VarResult % Values(Component::VarResult % Dofs) = Solver % Variable % Values
      END IF
    END DO

    IF( GotFact ) THEN
      CALL ListAddLogical( Params,'Linear System Refactorize',Refactorize )
    ELSE
      CALL ListRemove( Params,'Linear System Refactorize' )
    END IF

    Solver % Matrix % RHS => SaveRHS
    DEALLOCATE( ForceVector, LocalForce, SavedValues )
  END IF

  ! Each partition integrated over its own elements only
  IF( ASSOCIATED( ModeVariable(1) ) ) THEN
    DO i=1,3
      Pforce(i) = ParallelReduction( Pforce(i) )
      Vforce(i) = ParallelReduction( Vforce(i) )
      Sforce(i) = ParallelReduction( Sforce(i) )
      Moment(i) = ParallelReduction( Moment(i) )
    END DO
    DO i=1,3
      WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Pressure force ',i,' (N):',Pforce(i)
      CALL Info(Caller,Message,Level=5)
      CALL ListAddConstReal( Model % Simulation,'res: Pressure force '&
          //I2S(i),Pforce(i))
    END DO
    DO i=1,3
      WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Sliding force ',i,' (N):',Vforce(i)
      CALL Info(Caller,Message,Level=5)
      CALL ListAddConstReal( Model % Simulation,'res: Sliding force '&
          //I2S(i),Vforce(i))
    END DO
    TotForce = SQRT( SUM((Pforce + Vforce)**2) )
    WRITE(Message,'(A,T35,ES15.4)') 'Reynolds force (N): ',TotForce
    CALL Info(Caller,Message,Level=5)
    CALL ListAddConstReal( Model % Simulation,'res: Reynolds force',TotForce)

    IF(ASSOCIATED(SensVar)) THEN
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Pressure spring ',i,' (N/m):',Sforce(i)
        CALL Info(Caller,Message,Level=5)
        CALL ListAddConstReal( Model % Simulation,'res: Pressure spring '&
            //I2S(i),Sforce(i))
      END DO
    END IF

    IF( CalculateMoment ) THEN
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Reynolds moment ',i,' (Nm):',Moment(i)
        CALL Info(Caller,Message,Level=5)
        CALL ListAddConstReal( Model % Simulation,'res: Reynolds moment '&
            //I2S(i),Moment(i))
      END DO
    END IF

  END IF

  IF( ASSOCIATED( ModeVariable(4) ) ) THEN
    HeatPres = ParallelReduction( HeatPres )
    HeatSlide = ParallelReduction( HeatSlide )
    HeatTotal = HeatPres + HeatSlide
    WRITE(Message,'(A,T35,ES15.4)') 'Pressure heating (W): ',HeatPres
    CALL Info(Caller,Message,Level=5)
    WRITE(Message,'(A,T35,ES15.4)') 'Sliding heating (W): ',HeatSlide
    CALL Info(Caller,Message,Level=5)
    WRITE(Message,'(A,T35,ES15.4)') 'Reynolds heating (W): ',HeatTotal
    CALL Info(Caller,Message,Level=5)
    CALL ListAddConstReal( Model % Simulation,'res: Reynolds heating',HeatTotal)
  END IF


  CALL Info(Caller,'Finished computing postprocessing fields',Level=8)
  CALL Info(Caller,'--------------------------------------------------',Level=8)


CONTAINS



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(STIFF, FORCE, Element, n, nd, np, Nodes)
!------------------------------------------------------------------------------
! The mass matrix and one r.h.s. column for each (RhsMode,RhsComp) pair.
! np is the number of pressure dofs to use (nd, or n for nodal values only).
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:,:)
    INTEGER :: n, nd, np
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), TARGET :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3), DetJ
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S
    REAL(KIND=dp) :: Spres, Sslide
    REAL(KIND=dp) :: Normal(3), Velo(3), TangentVelo(3), Pres, Gap, GradPres(3), &
        Visc, mfp, Kn, TotPres, source, Radius(3), SensPres
    LOGICAL :: Stat
    INTEGER :: i,k,p,t,Mode,Component, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )

!------------------------------------------------------------------------------
    DO t=1,IntegStuff % n

      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo(Element, Nodes, u, v, w, DetJ, Basis, dBasisdx)

      s = s * DetJ
      IF ( CoordSys /= Cartesian .OR. CalculateMoment ) THEN
        x = SUM( Nodes % x(1:n) * Basis(1:n) )
        y = SUM( Nodes % y(1:n) * Basis(1:n) )
        z = SUM( Nodes % z(1:n) * Basis(1:n) )
      END IF

      IF( CoordSys /= Cartesian ) THEN
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF

      IF( CalculateMoment ) THEN
        Radius(1) = X - MomentAbout(1)
        Radius(2) = Y - MomentAbout(2)
        Radius(3) = Z - MomentAbout(3)
      END IF

!------------------------------------------------------------------------------
!      Parameters at integration point
!------------------------------------------------------------------------------

      IF( ASSOCIATED( Element % BoundaryInfo)) THEN
        Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
      ELSE
        Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
      END IF

      IF( UseVelocity ) THEN
        DO i=1,3
          Velo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
        TangentVelo = Velo - SUM(Normal * Velo) * Normal
      ELSE
        DO i=1,3
          TangentVelo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
      END IF

      Pres = SUM(Basis(1:np) * ElemPressure(1:np))
      TotPres = Pres + ReferencePressure

      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      DO i=1,3
        GradPres(i) = SUM(dBasisdx(1:np,i) * ElemPressure(1:np))
      END DO

      SensPres = 0.0_dp
      IF( ASSOCIATED(SensVar) ) SensPres = SUM(Basis(1:np) * SensPressure(1:np))

!------------------------------------------------------------------------------
!  Different viscosity models. Should be consistent with the main solver.
!------------------------------------------------------------------------------

      SELECT CASE (ViscosityType)

      CASE (Viscosity_Newtonian)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))

      CASE (Viscosity_Rarefied)
        Visc = SUM(Basis(1:n) * Viscosity(1:n))
        mfp = mfp0 * ReferencePressure / TotPres
        Kn = mfp / ABS(Gap)
        Visc = Visc / (1+9.638_dp*Kn**1.159_dp)

      END SELECT

!------------------------------------------------------------------------------
!  Coefficients of the differential equation at integration point
!------------------------------------------------------------------------------

      TotPres = TotPres - AmbientPres

      DO p=1,nd
        STIFF(p,1:nd) = STIFF(p,1:nd) + s * Basis(p) * Basis(1:nd)
      END DO

      DO k=1,nrhs
        Mode = RhsMode(k)
        Component = RhsComp(k)

        SELECT CASE( Mode )

        CASE( 1 )
          ! Forces resulting from pressure and shear
          Spres = -TotPres * Normal( Component )
          IF( OpposingWall ) Spres = -Spres
          Spres = Spres - Gap * GradPres( Component ) / 2

          Sslide = Visc * TangentVelo(Component ) / Gap

          source = Spres + Sslide

          Pforce(Component) = Pforce(Component) + s * Spres
          Vforce(Component) = Vforce(Component) + s * Sslide

          IF( ASSOCIATED(SensVar) ) THEN
            Sforce(Component) = Sforce(Component) + s * SensPres * Normal( Component )
          END IF

          IF( CalculateMoment ) THEN
            IF( Component == 1 ) THEN
              Moment(2) = Moment(2) + Radius(3) * s * source
              Moment(3) = Moment(3) - Radius(2) * s * source
            ELSE IF ( Component == 2 ) THEN
              Moment(1) = Moment(1) - Radius(3) * s * source
              Moment(3) = Moment(3) + Radius(1) * s * source
            ELSE IF ( Component == 3 ) THEN
              Moment(1) = Moment(1) + Radius(2) * s * source
              Moment(2) = Moment(2) - Radius(1) * s * source
            END IF
          END IF

        CASE( 2 )
          ! Flux resulting from pressure gradient and sliding
          Spres = - (Gap**3 / (12 * Visc) ) * GradPres(Component)
          Sslide = Gap * TangentVelo(Component) / 2
          source = Spres + Sslide

        CASE( 3 )
          ! Mean velocity resulting from pressure gradient and sliding.
          ! Compared to above mainly missing one Gap.
          Spres = - (Gap**2 / (12 * Visc) ) * GradPres(Component)
          Sslide = TangentVelo(Component) / 2
          source = Spres + Sslide

        CASE( 4 )
          ! heating effect of pressure gradient and sliding
          Spres = (Gap**3 / (12 * Visc) ) * SUM(GradPres *GradPres )
          Sslide = (Visc / Gap) * SUM(TangentVelo * TangentVelo )
          source = Spres + Sslide

          HeatPres = HeatPres + s * Spres
          HeatSlide = HeatSlide + s * Sslide

        CASE DEFAULT
          CALL Fatal(Caller,'Unknow mode: '//I2S(Mode))

        END SELECT

        FORCE(1:nd,k) = FORCE(1:nd,k) + s * source * Basis(1:nd)
      END DO
    END DO

  END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
! The result field of each mode, if requested.
!------------------------------------------------------------------------------
  FUNCTION ModeVariable( Mode ) RESULT ( Var )
    INTEGER :: Mode
    TYPE(Variable_t), POINTER :: Var

    SELECT CASE( Mode )
    CASE( 1 )
      Var => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Force')
    CASE( 2 )
      Var => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Flux')
    CASE( 3 )
      Var => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Mean Velocity')
    CASE( 4 )
      Var => VariableGet( Solver % Mesh % Variables,TRIM(PressureName)//' Heating')
    CASE DEFAULT
      Var => NULL()
    END SELECT
  END FUNCTION ModeVariable

!------------------------------------------------------------------------------
  END SUBROUTINE ReynoldsPostprocess
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. ReynoldsPostprocess
!> \ingroup Solvers
!------------------------------------------------------------------------------
  SUBROUTINE ReynoldsPostprocess_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, Calculate, ManningModel
    INTEGER :: Dim,GivenDim,dofs
    CHARACTER(LEN=MAX_NAME_LEN) :: PressureName
    CHARACTER(*), PARAMETER :: Caller = 'ReynoldsPostprocess_init'

!------------------------------------------------------------------------------
    Params => GetSolverParams()
    Dim = CoordinateSystemDimension()

    PressureName = GetString(Params,'Reynolds Pressure Variable Name',Found)
    IF(.NOT. Found) PressureName = 'FilmPressure'

    ! If the heating is not computed use a temp variable for the scalar equations
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Heating',Found)
    IF( Calculate ) THEN
      CALL ListAddString( Params,'Variable', &
          TRIM(PressureName)//' Heating' )
    ELSE IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
      CALL Info(Caller,'Defaulting field name to: ReynoldsPost')
      CALL ListAddString( Params,'Variable', &
          '-nooutput ReynoldsPost' )
    END IF

    ManningModel = ListGetLogical( Params,'Manning Model',Found )
    IF( ManningModel ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),&
          TRIM(PressureName)//' corrected' )
    END IF

    ! The dofs of force is fixed by default to 3 since there is a normal component
    ! as well as tangential components of force.
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Force',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Force Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(Caller,'Creating "'//TRIM(PressureName)//'" Force with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Force' )
    END IF

    ! The dofs of flux is fixed by default 3 since there can be leakage
    ! due to perforation even in 2d case.
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Flux',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Flux Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(Caller,'Creating "'//TRIM(PressureName)//' Flux" with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Flux' )
    END IF

    Calculate = ListGetLogical(Params,'Calculate Mean Velocity',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Mean Velocity Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        dofs = 2
      ELSE
        dofs = 3
      END IF
      CALL Info(Caller,'Creating "'//TRIM(PressureName)//' mean velocity" with '&
          //I2S(dofs)//' components',Level=12)
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
          '-dofs '//I2S(dofs)//' '//TRIM(PressureName)//' Mean Velocity' )
    END IF

    CALL ListAddInteger( Params, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+ILU0
    CALL ListAddNewString(Params,'Linear System Solver','Iterative')
    CALL ListAddNewString(Params,'Linear System Iterative Method','cg')
    CALL ListAddNewString(Params,'Linear System Preconditioning','ILU0')
    CALL ListAddNewInteger(Params,'Linear System Max Iterations',500)
    CALL ListAddNewInteger(Params,'Linear System Residual Output',10)
    CALL ListAddNewConstReal(Params,'Linear System Convergence Tolerance',1.0e-10_dp)

!------------------------------------------------------------------------------
  END SUBROUTINE ReynoldsPostprocess_Init
!------------------------------------------------------------------------------
