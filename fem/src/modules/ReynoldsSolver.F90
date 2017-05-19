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
  USE SolverUtils

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
  TYPE(ValueList_t), POINTER :: Params, Material, Equation, BC 

  INTEGER, PARAMETER :: Compressibility_None = 1, Compressibility_Weak = 2, &
      Compressibility_GasIsothermal = 3, Compressibility_GasAdiabatic = 4
  INTEGER, PARAMETER :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2

  INTEGER :: iter, i, j, k, l, n, nd, t, istat, mat_id, eq_id, body_id, mat_idold, &
      NoIterations, ViscosityType, CompressibilityType
  INTEGER, POINTER :: NodeIndexes(:), PressurePerm(:)

  LOGICAL :: GotIt, GotIt2, GotIt3, stat, AllocationsDone = .FALSE., SubroutineVisited = .FALSE., &
      UseVelocity, SideCorrection, Bubbles, ApplyLimiter
  REAL(KIND=dp), POINTER :: Pressure(:)
  REAL(KIND=dp) :: Norm, ReferencePressure, HeatRatio, BulkModulus, &
      mfp0, Pres, Dens
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:), TimeForce(:), &
      Viscosity(:), GapHeight(:), NormalVelocity(:), Velocity(:,:), &
      Admittance(:), Impedance(:), ElemPressure(:), PrevElemPressure(:)
  TYPE(Variable_t), POINTER :: SensVar, SaveVar

  CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityModel, CompressibilityModel

  SAVE ElementNodes, Viscosity, GapHeight, Velocity, NormalVelocity, &
      Admittance, FORCE, STIFF, MASS, TimeForce, ElemPressure, PrevElemPressure, &
      AllocationsDone


  CALL Info('ReynoldsSolver','---------------------------------------',Level=5)
  IF(TransientSimulation) THEN
    CALL Info('ReynoldsSolver','Solving the transient Reynolds equation',Level=5)
  ELSE
    CALL Info('ReynoldsSolver','Solving the steady-state Reynolds equation',Level=5)    
  END IF
  CALL Info('ReynoldsSolver','---------------------------------------',Level=5)

  CALL DefaultStart()
  
!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  
  Params => GetSolverParams()
  Bubbles = GetLogical( Params, 'Bubbles', GotIt )

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF(Solver % Variable % Dofs /= 1) THEN
    CALL Fatal('ReynoldsSolver','Impossible number of dofs! (should be 1)')    
  END IF
  Pressure     => Solver % Variable % Values
  PressurePerm => Solver % Variable % Perm
  IF( COUNT( PressurePerm > 0 ) <= 0) RETURN

!------------------------------------------------------------------------------
! Do some initial stuff
!------------------------------------------------------------------------------

  SideCorrection = ListGetLogicalAnyBC( Model,'Open Side')

  NoIterations = GetInteger( Params,'Nonlinear System Max Iterations',GotIt)
  IF(.NOT. GotIt) NoIterations = 1
  
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------


  IF ( .NOT. AllocationsDone  ) THEN
    N = Solver % Mesh % MaxElementNodes

    ALLOCATE(ElementNodes % x( N ),  &
        ElementNodes % y( N ),       &
        ElementNodes % z( N ),       &
        Viscosity( N ),              &
        GapHeight(N),          &
        Velocity(3,N),         &
        NormalVelocity(N),     &
        Admittance(N),         &
        FORCE( 2*N ),           &
        STIFF( 2*N, 2*N ), &
        MASS( 2*N, 2*N ), &
        TimeForce( 2*N ), &
        ElemPressure(N), &
        PrevElemPressure(N), &
        STAT=istat )
    IF ( istat /= 0 ) CALL FATAL('ReynoldsSolver','Memory allocation error')

    AllocationsDone = .TRUE.
  END IF


!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
  
  mat_idold = 0

  CALL Info('ReynoldsSolver','-------------------------------------------------',Level=5)

  DO iter = 1,NoIterations

    WRITE(Message,'(A,T35,I5)') 'Reynolds iteration:',iter
    CALL Info('ReynoldsSolver',Message,Level=5)

    CALL DefaultInitialize()

!    Do the bulk assembly:
!    ---------------------
    CALL GlobalBulkAssembly()
    CALL DefaultFinishBulkAssembly( )

!------------------------------------------------------------------------------
!    Neumann & Newton BCs:
!------------------------------------------------------------------------------
    IF(SideCorrection) THEN
      CALL GlobalBoundaryAssemby()
    END IF
    
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

!    Solve the system and we are done:
!    ---------------------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT
  END DO
  
  IF( ListGetLogical( Params,'Gap Sensitivity', GotIt ) ) THEN
    CALL Info('ReynoldsSolver','Computing FilmPressure sentivity to gap height',Level=5)

    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.TRUE.)
    ApplyLimiter = ListGetLogical( Params,'Apply Limiter', GotIt )
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .FALSE. ) 
    
    SensVar => VariableGet( Model % Variables,'FilmPressure Gap Sensitivity')
    IF( .NOT. ASSOCIATED( SensVar ) ) THEN
      CALL Fatal('ReynoldsSolver','> Filmpressure gap sensitivity < should exist!')
    END IF
    SaveVar => Solver % Variable 
    Solver % Variable => SensVar

    CALL DefaultInitialize()
    CALL GlobalBulkAssembly( 1 )
    CALL DefaultFinishBulkAssembly( )
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs( Ux = SensVar )

!    Solve the system and we are done:
!    ---------------------------------
    Norm = DefaultSolve()

    CALL ListAddLogical(Params,'Skip Compute Nonlinear Change',.FALSE.)
    Solver % Variable => SaveVar
    IF( ApplyLimiter ) CALL ListAddLogical( Params,'Apply Limiter', .TRUE. ) 
  END IF

  CALL DefaultFinish() 

  CALL Info('ReynoldsSolver','-------------------------------------------------',Level=5)
  
CONTAINS  



  ! Cycle over the bulk elements and build the global linear system
  ! Optionally performs sensitisity analysis. 
  !----------------------------------------------------------------------
  SUBROUTINE GlobalBulkAssembly( SensitivityMode )

    INTEGER, OPTIONAL :: SensitivityMode
    INTEGER :: SensMode 
    
    SensMode = 0
    IF( PRESENT(SensitivityMode) ) SensMode = SensitivityMode


    DO t=1,Solver % NumberOfActiveElements

      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      
      CALL GetElementNodes( ElementNodes )
      CALL GetScalarLocalSolution( ElemPressure )

      IF( SensMode > 0 ) THEN
        IF( TransientSimulation ) THEN
          CALL GetScalarLocalSolution( PrevElemPressure, tstep = -1 )
        END IF
      END IF


      body_id =  Element % Bodyid

      eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
      Equation => Model % Equations(eq_id) % Values

      mat_id = GetInteger( Model % Bodies( body_id ) % Values, 'Material')
      Material => Model % Materials(mat_id) % Values

!------------------------------------------------------------------------------
!       Get velocities
!------------------------------------------------------------------------------        

      Velocity = 0.0_dp
      UseVelocity = .FALSE.
      IF( ListCheckPrefix( Equation,'Surface Velocity') ) THEN
        Velocity(1,1:n) = GetReal(Equation,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(Equation,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(Equation,'Surface Velocity 3',GotIt3)
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
        IF( ListCheckPrefix( Equation,'Tangent Velocity') ) THEN
          Velocity(1,1:n) = GetReal(Equation,'Tangent Velocity 1',GotIt) 
          Velocity(2,1:n) = GetReal(Equation,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Equation,'Tangent Velocity 3',GotIt3)
        END IF
        IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3)) THEN
          IF( ListCheckPrefix( Material,'Tangent Velocity') ) THEN
            Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt) 
            Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
            Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
          END IF
        END IF
        NormalVelocity(1:n) = GetReal(Equation,'Normal Velocity',GotIt)
        IF(.NOT. GotIt) NormalVelocity(1:n) = &
          GetReal(Material,'Normal Velocity',GotIt)
      END IF     

!------------------------------------------------------------------------------
!       Get material parameters
!------------------------------------------------------------------------------        

      GapHeight(1:n) = GetReal( Material,'Gap Height')
      Admittance(1:n) = GetReal( Material, 'Flow Admittance', GotIt)
      Viscosity(1:n) = GetReal( Material, 'Viscosity')

      IF(mat_id /= mat_idold) THEN                  

        mat_idold = mat_id

        ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
        IF(GotIt) THEN
          IF( ViscosityModel == 'newtonian') THEN
            ViscosityType = Viscosity_Newtonian
          ELSE IF( ViscosityModel == 'rarefied') THEN
            ViscosityType = Viscosity_Rarefied
            mfp0 = GetCReal(Material,'Mean Free Path')            
            ReferencePressure = GetCReal( Material,'Reference Pressure')           
          ELSE
            CALL Warn('ReynoldsSolver','Unknown viscosity model')
          END IF
        ELSE
          ViscosityType = Viscosity_Newtonian          
        END IF

        CompressibilityModel = GetString(Material,'Compressibility Model',GotIt)        
        IF(GotIt) THEN
          IF(CompressibilityModel == 'incompressible') THEN
            CompressibilityType = Compressibility_None
          ELSE IF(CompressibilityModel == 'weakly compressible') THEN
            CompressibilityType = Compressibility_Weak
            ReferencePressure = GetCReal( Material,'Reference Pressure',GotIt)           
            BulkModulus = GetCReal( Material, 'Bulk Modulus')
          ELSE IF(CompressibilityModel == 'isothermal ideal gas') THEN
            CompressibilityType = Compressibility_GasIsothermal
            ReferencePressure = GetCReal( Material,'Reference Pressure')           
          ELSE IF(CompressibilityModel == 'adiabatic ideal gas') THEN
            CompressibilityType = Compressibility_GasAdiabatic
            HeatRatio = GetCReal( Material, 'Specific Heat Ratio')
            ReferencePressure = GetCReal( Material,'Reference Pressure')                      
          ELSE
            CALL Warn('ReynoldsSolver','Unknown compressibility model')
          END IF
        ELSE          
          CompressibilityType = Compressibility_None
        END IF
      END IF

      STIFF = 0.0d0
      MASS = 0.0d0
      FORCE = 0.0d0
      
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
  SUBROUTINE LocalBulkMatrix(MassMatrix, StiffMatrix, ForceVector, &
      Element, n, nd, Nodes, SensMode )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MassMatrix(:,:), StiffMatrix(:,:), ForceVector(:)
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    INTEGER :: SensMode
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), SqrtElementMetric
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S, MS, MM, L, A, B, HR, SL(3), SLR, SLL(3), F
    REAL(KIND=dp) :: Normal(3), Velo(3), NormalVelo, TangentVelo(3), Damp, Pres, PrevPres, &
        TotPres, GradPres(3), dPdt, Gap, Visc, mfp, Kn, Density, DensityDer
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,DIM, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

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

      U = IntegStuff % u(t)
      V = IntegStuff % v(t)
      W = IntegStuff % w(t)
      S = IntegStuff % s(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx, Bubbles = Bubbles)
    
      s = s * SqrtElementMetric
      IF ( CoordSys /= Cartesian ) THEN
        X = SUM( Nodes % X(1:n) * Basis(1:n) )
        Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
        Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
        s = s * SqrtMetric
      END IF

!------------------------------------------------------------------------------
!      Parameters at integration point
!------------------------------------------------------------------------------

      IF(UseVelocity) THEN
        IF( ASSOCIATED( Element % BoundaryInfo)) THEN
          Normal = NormalVector( Element,Nodes,u,v,.TRUE. )
        ELSE
          Normal = NormalVector( Element,Nodes,u,v,.FALSE. )
        END IF
        DO i=1,3
          Velo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
        NormalVelo = SUM( Normal * Velo)
        TangentVelo = Velo - NormalVelo * Normal
      ELSE
        NormalVelo = SUM( Basis(1:n) * NormalVelocity(1:n) )
        DO i=1,3
          TangentVelo(i) = SUM( Basis(1:n) * Velocity(i,1:n) )
        END DO
      END IF

      Damp = SUM(Basis(1:n) * Admittance(1:n))
      Pres = SUM(Basis(1:n) * ElemPressure(1:n))
      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      TotPres = ReferencePressure + Pres

      ! If we compute sensitivity of solution we need various derivatives of pressure
      !--------------------------------------------------------------------------------
      IF( SensMode > 0 ) THEN
        IF( TransientSimulation ) THEN
          PrevPres = SUM( Basis(1:n) * PrevElemPressure(1:n) ) 
          dPdt = ( Pres - PrevPres ) / dt 
        ELSE 
          dPdt = 0.0_dp
        END IF        
        DO i = 1,3
          GradPres(i) = SUM( dBasisdx(1:n,i) * ElemPressure(1:n) )
        END DO
      END IF


!------------------------------------------------------------------------------
!  Different material models. The "density" is used only as a functional form,
!  not as absolute value.
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

      SELECT CASE (CompressibilityType) 

      CASE (Compressibility_None)
        Density = 1.0d0
        DensityDer = 0.0d0
        
      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
        DensityDer = Density / BulkModulus
        
      CASE (Compressibility_GasIsothermal)         
        Density = TotPres 
        DensityDer = 1.0d0
        
      CASE(Compressibility_GasAdiabatic) 
        Density = TotPres ** (1.0_dp/HeatRatio)
        DensityDer = (1/HeatRatio) * TotPres ** (1.0_dp/HeatRatio - 1.0_dp)
        
      END SELECT
      
!------------------------------------------------------------------------------
!  Coefficients of the differential equation at integration point
!------------------------------------------------------------------------------

      ! Multipliers of p: Stiffness matrix 
      MS = -Density * Gap**3 / (12 * Visc)
      HR = -Damp * Density

      ! Multipliers of dp/dt: Mass matrix 
      MM = -DensityDer * Gap
      
      ! Normal velocity: right-hand-side force vector
      L = Density * NormalVelo

      ! Tangential velocity: Both rhs and matix contribution
      SLR = 0.0d0
      SLL = 0.0d0

      DO i=1,dim
        ! The plane element automatically omits the derivative in normal direction
        SLR = SLR + 0.5_dp * Density * SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n)) 
      END DO
      ! Implicit part: coefficient of the pressure gradient
      SLL = -0.5_dp* DensityDer * Gap * Velo 

!------------------------------------------------------------------------------
!      The Reynolds equation
!------------------------------------------------------------------------------
      DO p=1,NBasis
        DO q=1,NBasis
          A = HR * Basis(q) * Basis(p)           
          DO i=1,DIM
            DO j=1,DIM
              A = A + MS * Metric(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
              A = A + SLL(j) * Metric(i,j) * dBasisdx(q,i) * Basis(p)
            END DO
          END DO          
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A 

          IF( TransientSimulation ) THEN
            B = MM * Basis(q) * Basis(p)
            MassMatrix(p,q)  = MassMatrix(p,q)  + s * B
          END IF
        END DO

        F = 0.0_dp
        IF( SensMode == 0 ) THEN
          F = L + SLR
        ELSE IF( SensMode == 1 ) THEN
          IF( TransientSimulation ) THEN
            F = -2.0 * DensityDer * dPdt
          END IF
          F = F - DensityDer * SUM( TangentVelo * GradPres )          
          DO i = 1,dim
            ! The plane element automatically omits the derivative in normal direction
            F = F - 1.5_dp * ( Density / Gap ) * SUM( dBasisdx(1:n,i) * Velocity(i,1:n) * GapHeight(1:n) ) 
          END DO
          F = F - 3.0_dp * Damp * Density * Pres / Gap 
          F = F - 3 * Density * NormalVelo / Gap
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
   SUBROUTINE GlobalBoundaryAssemby()

    DO t=1, Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE
      
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      IF ( GetElementFamily() == 1 ) CYCLE
      
      BC => GetBC()
      IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
      
      stat = GetLogical(BC,'Open Side',gotIt) 
      IF(.NOT. stat) CYCLE
!------------------------------------------------------------------------------
      NodeIndexes => Element % NodeIndexes
      
      IF ( ANY( PressurePerm(NodeIndexes(1:n)) == 0 ) ) CYCLE
      
      Parent => Element % BoundaryInfo % Left
      stat = ASSOCIATED( Parent )
      IF ( stat ) stat = stat .AND. ALL(PressurePerm(Parent % NodeIndexes) > 0)
      
      IF(.NOT. stat) THEN
        Parent => ELement % BoundaryInfo % Right            
        stat = ASSOCIATED( Parent )
        IF ( stat ) stat = stat .AND. ALL(PressurePerm(Parent % NodeIndexes) > 0)
        IF ( .NOT. stat )  CALL Fatal( 'ReynoldsSolver', &
            'No proper parent element available for specified boundary' )
      END IF
      
      Model % CurrentElement => Parent
      CALL GetElementNodes( ElementNodes )
      
      mat_id = GetInteger( Model % Bodies(Parent % BodyId) % Values,'Material')
      Material => Model % Materials(mat_id) % Values
      
      GapHeight(1:n) = GetReal(Material,'Gap Height')
      
      Viscosity(1:n) = GetReal( Material, 'Viscosity')
      
      STIFF = 0.0d0
      MASS = 0.0d0
      FORCE = 0.0d0
      
!------------------------------------------------------------------------------
!             Get element local matrix and rhs vector
!------------------------------------------------------------------------------
      
      CALL LocalBoundaryMatrix( MASS, STIFF, FORCE, Element, n, ElementNodes )
                    
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
  SUBROUTINE LocalBoundaryMatrix(MassMatrix, StiffMatrix, ForceVector, Element, n, Nodes)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MassMatrix(:,:), StiffMatrix(:,:), ForceVector(:)
    INTEGER :: n
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: Impedance 
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),X,Y,Z
    REAL(KIND=dp) :: Visc, dl, mfp, Kn, Damp, TotPres, Pres, Density, Gap, A
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
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx )
      
      s = s * SqrtElementMetric
      
      Gap = SUM( GapHeight(1:n) * Basis(1:n) )
      Pres = SUM( ElemPressure(1:n) * Basis(1:n))
      TotPres = ReferencePressure + Pres
      
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
      END SELECT
      
      Damp = Gap / (12 * dl * Visc)
      
!------------------------------------------------------------------------------
       
      SELECT CASE (CompressibilityType) 
        
      CASE (Compressibility_None)
        Density = 1.0d0
        
      CASE (Compressibility_Weak)
        Density = EXP(TotPres/BulkModulus)
        
      CASE (Compressibility_GasIsothermal)         
        Density = TotPres 
        
      CASE(Compressibility_GasAdiabatic) 
        Density = TotPres ** (1.0/HeatRatio)
        
      END SELECT

!------------------------------------------------------------------------------       
       A = -Damp * Gap * Density
!------------------------------------------------------------------------------
       DO p=1,n
         DO q=1,n
           StiffMatrix(p,q) = StiffMatrix(p,q) + s * Basis(q) * Basis(p) * A
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
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  TYPE(ValueList_t), POINTER :: Params 


  Params => GetSolverParams()

! The new way with generic limiters is a library functionality.
! The Poisson equation is assembled using different sign that the 
! typical convention. Hence the load sign for limiters is opposite
!----------------------------------------------------------------
  CALL ListAddLogical( Params,'Limiter Load Sign Negative',.TRUE.)

  IF( ListGetLogical( Params,'Gap Sensitivity', Found ) ) THEN
    CALL ListAddStrinG( Params,NextFreeKeyword('Exported Variable',Params),&
        'FilmPressure Gap Sensitivity')
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
  USE SolverUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: PressureVar, VarResult, SolverVar
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: Element, Parent
  TYPE(ValueList_t), POINTER :: Params, Material, Equation

  INTEGER, PARAMETER :: Viscosity_Newtonian = 1, Viscosity_Rarefied = 2

  INTEGER :: iter, i, j, k, l, n, nd, t, istat, eq_id, body_id, mat_id, mat_idold,&
      ViscosityType, dim
  INTEGER :: Component, Components, Mode
  INTEGER, POINTER :: NodeIndexes(:), PressurePerm(:)
  REAL(KIND=dp), POINTER :: mWork(:,:)	

  LOGICAL :: GotIt, GotIt2, GotIt3, stat, UseVelocity, AllocationsDone = .FALSE., &
      OpposingWall, CalculateMoment

  REAL(KIND=dp), POINTER :: Pressure(:)
  REAL(KIND=dp) :: Norm, ReferencePressure, mfp0, HeatSlide, HeatPres, HeatTotal, &
      Pforce(3), Vforce(3), TotForce, Moment(3), MomentAbout(3)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), Viscosity(:), GapHeight(:), &
      Velocity(:,:), ElemPressure(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: ViscosityModel, PressureName


  SAVE ElementNodes, Viscosity, Velocity, &
      GapHeight, FORCE, STIFF, ElemPressure, AllocationsDone

 
!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------

  CALL Info('ReynoldsPostprocess','--------------------------------------------------',Level=5)
  CALL Info('ReynoldsPostprocess','Computing postprocessing fields from film pressure',Level=5)
  CALL Info('ReynoldsPostprocess','--------------------------------------------------',Level=5)

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------

  Params => GetSolverParams()

  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF(Solver % Variable % Dofs /= 1) THEN
    CALL Fatal('ReynoldsPostprocess','Impossible number of dofs! (should be 1)')    
  END IF

  PressureName = GetString(Params,'Reynolds Pressure Variable Name',GotIt)
  IF(.NOT. GotIt) PressureName = 'FilmPressure'

  PressureVar => VariableGet( Solver % Mesh % Variables, PressureName)
  IF(.NOT. ASSOCIATED(PressureVar)) THEN
    CALL Warn('ReynoldsPostprocess','Could not get variable: '//TRIM(PressureName))
    RETURN
  END IF

  Pressure => PressureVar % Values
  PressurePerm => PressureVar % Perm
  IF( COUNT( PressurePerm > 0 ) <= 0) RETURN

  DIM = CoordinateSystemDimension()

!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone  ) THEN
    N = Solver % Mesh % MaxElementNodes

    ALLOCATE(ElementNodes % x( N ),  &
        ElementNodes % y( N ),       &
        ElementNodes % z( N ),       &
        Viscosity( N ),              &
        GapHeight(N),          &
        Velocity(3,N),         &
        FORCE( N ),           &
        STIFF( N, N ), &
        ElemPressure(N), &
        STAT=istat )

    IF ( istat /= 0 ) CALL FATAL('ReynoldsPostprocess','Memory allocation error')    

    AllocationsDone = .TRUE.
  END IF 

  mWork => ListGetConstRealArray( Params,'Moment About',CalculateMoment)
  IF( CalculateMoment ) THEN
    MomentAbout(1:dim) = mWork(1:dim,1)
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
  Vforce = 0.0_dp
  Moment = 0.0_dp

  SolverVar => Solver % Variable
  IF( .NOT. ASSOCIATED( SolverVar ) ) THEN
    CALL Fatal('ReynoldsPostprocess','Solver Variable not associated')
  END IF


  CALL Info('ReynoldsPostprocess','Primary variable name: '//TRIM( Solver % variable % Name) )

   
  DO Mode = 1, 3  

    IF( Mode == 1 ) THEN
      VarResult => VariableGet( Solver % Mesh % Variables,'FilmPressure Force')
      IF(.NOT. ASSOCIATED(VarResult)) CYCLE
    ELSE IF( Mode == 2 ) THEN
      VarResult => VariableGet( Solver % Mesh % Variables,'FilmPressure Flux')
      IF(.NOT. ASSOCIATED(VarResult)) CYCLE
    ELSE     
      VarResult => VariableGet( Solver % Mesh % Variables,'FilmPressure Heating')
      IF(.NOT. ASSOCIATED(VarResult)) CYCLE
    END IF
    Components = VarResult % Dofs

    DO Component = 1, Components

      CALL DefaultInitialize()

      !    Do the bulk assembly:
      !    ---------------------
      
      DO t=1,Solver % NumberOfActiveElements
        
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        
        CALL GetElementNodes( ElementNodes )
        
        NodeIndexes => Element % NodeIndexes         
        ElemPressure(1:n) = Pressure(PressurePerm(NodeIndexes(1:n)))
        
        body_id =  Element % Bodyid
	IF( body_id <= 0 ) CYCLE        

        eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation')
        Equation => Model % Equations(eq_id) % Values
        
        mat_id = GetInteger( Model % Bodies( body_id ) % Values, 'Material')
        Material => Model % Materials(mat_id) % Values

        !------------------------------------------------------------------------------
        !       Get velocities
        !------------------------------------------------------------------------------                
        Velocity(1,1:n) = GetReal(Equation,'Surface Velocity 1',GotIt)
        Velocity(2,1:n) = GetReal(Equation,'Surface Velocity 2',GotIt2)
        Velocity(3,1:n) = GetReal(Equation,'Surface Velocity 3',GotIt3)
        UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
        IF(.NOT. UseVelocity) THEN
          Velocity(1,1:n) = GetReal(Material,'Surface Velocity 1',GotIt)
          Velocity(2,1:n) = GetReal(Material,'Surface Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Material,'Surface Velocity 3',GotIt3)
          UseVelocity = GotIt .OR. GotIt2 .OR. GotIt3
        END IF
        
        IF(.NOT. UseVelocity) THEN
          Velocity(1,1:n) = GetReal(Equation,'Tangent Velocity 1',GotIt) 
          Velocity(2,1:n) = GetReal(Equation,'Tangent Velocity 2',GotIt2)
          Velocity(3,1:n) = GetReal(Equation,'Tangent Velocity 3',GotIt3)
          IF(.NOT. (GotIt .OR. GotIt2 .OR. GotIt3)) THEN
            Velocity(1,1:n) = GetReal(Material,'Tangent Velocity 1',GotIt) 
            Velocity(2,1:n) = GetReal(Material,'Tangent Velocity 2',GotIt2)
            Velocity(3,1:n) = GetReal(Material,'Tangent Velocity 3',GotIt3)
          END IF
        END IF
        

        !------------------------------------------------------------------------------
        !       Get material parameters
        !------------------------------------------------------------------------------                
        GapHeight(1:n) = GetReal(Material,'Gap Height')
        Viscosity(1:n) = GetReal( Material, 'Viscosity')
        
        IF(mat_id /= mat_idold) THEN                  
          
          mat_idold = mat_id
          
          ViscosityModel = GetString(Material,'Viscosity Model',GotIt)
          IF(GotIt) THEN
            IF( ViscosityModel == 'newtonian') THEN
              ViscosityType = Viscosity_Newtonian
            ELSE IF( ViscosityModel == 'rarefied') THEN
              ViscosityType = Viscosity_Rarefied
              mfp0 = GetCReal(Material,'Mean Free Path')            
              ReferencePressure = GetCReal( Material,'Reference Pressure')           
            ELSE
              CALL Warn('ReynoldsPostprocess','Unknown viscosity model')
            END IF
          ELSE
            ViscosityType = Viscosity_Newtonian          
          END IF

          OpposingWall = ListGetLogical( Material,'Opposing Wall',GotIt)
        END IF
        
        STIFF = 0.0d0
        FORCE = 0.0d0
        
        CALL LocalMatrix( STIFF, FORCE, Element, n, nd, ElementNodes) 
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
      
      !------------------------------------------------------------------------------
      
      CALL DefaultFinishAssembly()


      ! There could be some beriodic BCs hence the BCs
      !-----------------------------------------------
      CALL DefaultDirichletBCs()
      CALL Info( 'ReynoldsPostprocess', 'Dirichlet conditions done', Level=4 )
      
      ! Solve the system and we are done:
      !------------------------------------
      Norm = DefaultSolve()
      
      IF( Mode /= 3 ) THEN
        VarResult % Values(Component::Components) = Solver % Variable % Values
      END IF
    END DO


    IF( Mode == 1 ) THEN
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Pressure force ',i,' (N):',Pforce(i)
        CALL Info('ReynoldsPostprocess',Message,Level=5)
        CALL ListAddConstReal( Model % Simulation,'res: Pressure force '&
            //TRIM(I2S(i)),Pforce(i))
      END DO
      DO i=1,3
        WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Sliding force ',i,' (N):',Vforce(i)
        CALL Info('ReynoldsPostprocess',Message,Level=5)
        CALL ListAddConstReal( Model % Simulation,'res: Sliding force '&
            //TRIM(I2S(i)),Vforce(i))
      END DO
      TotForce = SQRT( SUM((Pforce + Vforce)**2) )
      WRITE(Message,'(A,T35,ES15.4)') 'Reynolds force (N): ',TotForce
      CALL Info('ReynoldsPostprocess',Message,Level=5)
      CALL ListAddConstReal( Model % Simulation,'res: Reynolds force',TotForce)

      IF( CalculateMoment ) THEN
        DO i=1,3
          WRITE(Message,'(A,I1,A,T35,ES15.4)') 'Reynolds moment ',i,' (Nm):',Moment(i)
          CALL Info('ReynoldsPostprocess',Message,Level=5)
          CALL ListAddConstReal( Model % Simulation,'res: Reynolds moment '&
              //TRIM(I2S(i)),Moment(i))
        END DO
      END IF


    ELSE IF( Mode == 2 ) THEN
      

    ELSE
      HeatTotal = HeatPres + HeatSlide
      WRITE(Message,'(A,T35,ES15.4)') 'Pressure heating (W): ',HeatPres
      CALL Info('ReynoldsPostprocess',Message,Level=5)
      WRITE(Message,'(A,T35,ES15.4)') 'Sliding heating (W): ',HeatSlide
      CALL Info('ReynoldsPostprocess',Message,Level=5)
      WRITE(Message,'(A,T35,ES15.4)') 'Reynolds heating (W): ',HeatTotal
      CALL Info('ReynoldsPostprocess',Message,Level=5)
      CALL ListAddConstReal( Model % Simulation,'res: Reynolds heating',HeatTotal)
    END IF

  END DO
  

CONTAINS



!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(StiffMatrix, ForceVector, Element, n, nd, Nodes)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: StiffMatrix(:,:), ForceVector(:)
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3), SqrtElementMetric
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: U, V, W, S, A
    REAL(KIND=dp) :: Vpres(3), Vslide(3), Spres, Sslide
    REAL(KIND=dp) :: Normal(3), Velo(3), TangentVelo(3), Pres, Gap, GradPres(3), &
        Visc, mfp, Kn, TotPres, source, Radius(3)
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0d0
    Metric(1,1) = 1.0d0
    Metric(2,2) = 1.0d0
    Metric(3,3) = 1.0d0

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------

    NBasis = nd
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
      stat = ElementInfo( Element, Nodes, U, V, W, SqrtElementMetric, &
          Basis, dBasisdx)
      
      s = s * SqrtElementMetric
      IF ( CoordSys /= Cartesian .OR. CalculateMoment ) THEN
        X = SUM( Nodes % X(1:n) * Basis(1:n) )
        Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
        Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
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
      
      Pres = SUM(Basis(1:n) * ElemPressure(1:n))
      TotPres = ReferencePressure + Pres
      Gap = SUM(Basis(1:n) * GapHeight(1:n))
      DO i=1,dim
        GradPres(i) = SUM(dBasisdx(1:n,i) * ElemPressure(1:n))
      END DO

!------------------------------------------------------------------------------
!  Different viscosity models. Should be consistant with the main solver.
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

      IF( Mode == 1 ) THEN
        ! Forces resulting from pressure and shear
        Spres = -TotPres * Normal( Component ) 
        IF( OpposingWall ) Spres = -Spres
        Spres = Spres - Gap * GradPres( Component ) / 2
        
        Sslide = Visc * TangentVelo(Component ) / Gap
        
        source = Spres + Sslide
        
        Pforce(Component) = Pforce(Component) + s * Spres
        Vforce(Component) = Vforce(Component) + s * Sslide

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
        
      ELSE IF( Mode == 2 ) THEN
        ! Flux resulting from pressure gradient and sliding 

        Spres = - (Gap**3 / (12 * Visc) ) * GradPres(Component)
        Sslide = Gap * TangentVelo(Component) / 2        
        ! add contribution of leaking
        
        source = Spres + Sslide 
      ELSE      
        ! heating effect of pressure gradient and sliding
        Spres = (Gap**3 / (12 * Visc) ) * SUM(GradPres *GradPres )
        Sslide = (Visc / Gap) * SUM(TangentVelo * TangentVelo )
        
        source = Spres + Sslide
        
        HeatPres = HeatPres + s * Spres
        HeatSlide = HeatSlide + s * Sslide
      END IF
      
      DO p=1,NBasis
        DO q=1,NBasis
          A = Basis(q) * Basis(p)           
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A 
        END DO
        ForceVector(p) = ForceVector(p) + s * Basis(p) * source
      END DO
    END DO

  END SUBROUTINE LocalMatrix

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
    LOGICAL :: Found, Calculate
    INTEGER :: Dim,GivenDim
!------------------------------------------------------------------------------
    Params => GetSolverParams()
    Dim = CoordinateSystemDimension()

    ! If the heating is not computed use a temp variable for the scalar equations
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Heating',Found)
    IF( Calculate ) THEN
      CALL ListAddString( Params,'Variable', &
          'FilmPressure Heating' )
    ELSE IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
      CALL ListAddString( Params,'Variable', &
          '-nooutput ReynoldsPost' )      
    END IF

    ! The dofs of force is fixed by default to 3 since there is a normal component
    ! as well as tangential components of force.
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Force',Found)
  
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Force Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
            '-dofs 2 FilmPressure Force' )
      ELSE
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
            '-dofs 3 FilmPressure Force' )
      END IF
    END IF

    ! The dofs of flux is fixed by default 3 since there can be leakage 
    ! due to perforation even in 2d case.
    !-------------------------------------------------------------------
    Calculate = ListGetLogical(Params,'Calculate Flux',Found)
    IF( Calculate ) THEN
      GivenDim = ListGetInteger(Params,'Calculate Flux Dim',Found)
      IF( Dim == 1 .OR. GivenDim == 2 ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            '-dofs 2 FilmPressure Flux' )
      ELSE
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            '-dofs 3 FilmPressure Flux' )
      END IF
    END IF

    CALL ListAddInteger( Params, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+diagonal
    IF(.NOT. ListCheckPresent(Params,'Linear System Solver')) &
      CALL ListAddString(Params,'Linear System Solver','Iterative')
    IF(.NOT. ListCheckPresent(Params,'Linear System Iterative Method')) &
      CALL ListAddString(Params,'Linear System Iterative Method','cg')
    IF(.NOT. ListCheckPresent(Params,'Linear System Preconditioning')) &
      CALL ListAddString(Params,'Linear System Preconditioning','diagonal')
    IF(.NOT. ListCheckPresent(Params,'Linear System Max Iterations')) &
      CALL ListAddInteger(Params,'Linear System Max Iterations',500)
    IF(.NOT. ListCheckPresent(Params,'Linear System Residual Output')) &
      CALL ListAddInteger(Params,'Linear System Residual Output',10)
    IF(.NOT. ListCheckPresent(Params,'Linear System Convergence Tolerance')) &
      CALL ListAddConstReal(Params,'Linear System Convergence Tolerance',1.0e-10_dp)

!------------------------------------------------------------------------------
  END SUBROUTINE ReynoldsPostprocess_Init
!------------------------------------------------------------------------------
