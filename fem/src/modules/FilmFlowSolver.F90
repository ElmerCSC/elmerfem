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
! * Module for the solution of reduced dimensional Navier-Stokes equation.
! * Here we assume that the average velocity of the film defines the full velocity
! * profile and hence the depth direction may be eliminated from the flow solution.
! * The intended use is film/channel flow but also cylindrical pipe flow is implemented.
! * This fills the gap between full Navier-Stokes solver and reduced dimensional
! * Reynolds solver. 
! *
! * The module is compatible with p-bubbles and/or p2/p1 elements, e.g.
! * 303b1    - triangle with one bubble
! * 303e1b1  - p2/p1 triangle!
! * 404b4    - quad with four bubbles
! * 404e1b1  - q2/q1 quad
! *
! * This module has been defived from a historical Navier-Stokes solver and
! * updated much later for problems involving challel flows.
! *
! *  Authors: Juha Ruokolainen, Peter Råback, Thomas Zwinger, Tómas Jóhannesson
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 27.10.2022
! *
!/*****************************************************************************/


!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver_init0( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  LOGICAL :: Found, Serendipity
  TYPE(ValueList_t), POINTER :: Params 
  Params => GetSolverParams()

  Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found)
  IF(.NOT.Found) Serendipity = .TRUE.

  CALL ListAddNewInteger(Params, 'Time derivative Order', 1)
  
  IF(Serendipity) THEN
    CALL ListAddNewString(Params,'Element','p:1 -line b:1 -tri b:1 -quad b:3')
  ELSE
    CALL ListAddNewString(Params,'Element','p:1 -line b:1 -tri b:1 -quad b:4')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE FilmFlowSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver_init(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params 
  LOGICAL :: Found, Found2
  INTEGER :: mdim
  CHARACTER(*), PARAMETER :: Caller = 'FilmFlowSolver_init'
!------------------------------------------------------------------------------ 
  Params => GetSolverParams() 
  
  mdim = ListGetInteger(Params,'Model Dimension',Found)
  IF(.NOT. Found ) THEN
    CALL Fatal(Caller,'Give "Model Dimension" i.e. the dimension of N-S equation!')    
  END IF
   
  IF ( mdim == 2 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow[FilmVelocity:2 FilmPressure:1]')
  ELSE IF( mdim == 1 ) THEN
    CALL ListAddNewString(Params, 'Variable', &
        'Flow[FilmSpeed:1 FilmPressure:1]')
  ELSE
    CALL Fatal(Caller,'This module does not make sense in dim: '//I2S(mdim))    
  END IF

  ! Study only velocity components in linear system
  CALL ListAddNewInteger(Params, 'Nonlinear System Norm DOFs', mdim )

  ! This should be true to incompressible flows where pressure level is not uniquely determined
  CALL ListAddNewLogical(Params, 'Relative Pressure Relaxation', .TRUE. )

  ! It makes sense to eliminate the bubbles to save memory and time
  CALL ListAddNewLogical(Params, 'Bubbles in Global System', .FALSE.)

  ! Use global mass matrix in time integration
  CALL ListAddNewLogical(Params, 'Global Mass Matrix', .TRUE.)

  IF (ListGetLogical(Params,'Calculate Friction Heating',Found)) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'FilmFriction HeatFlux' )
  END IF
  IF (ListGetLogical(Params,'Calculate Pressure Heating',Found2)) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'FilmPressure HeatFlux' )
  END IF
  IF( Found .OR. Found2 ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Heating Energy' )
  END IF
    
  
  
!------------------------------------------------------------------------------ 
END SUBROUTINE FilmFlowSolver_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE FilmFlowSolver( Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found, Convect, CSymmetry
  TYPE(Element_t),POINTER :: Element
  INTEGER :: i,n, nb, nd, t, istat, dim, mdim, BDOFs=1,Active,iter,maxiter,CoupledIter
  REAL(KIND=dp) :: Norm = 0, mingap, Grav, time, time0=-1.0_dp, MeltHeat, Cp, Ct
  TYPE(ValueList_t), POINTER :: Params, BodyForce, Material, BC
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), &
      FORCE(:), rho(:), gap(:), gap0(:), mu(:), ac(:), height(:), AcPres(:), &
      Velocity(:,:), MASS(:,:), AcPrevPressure(:), FsiRhs(:,:), PrevGap(:)
  REAL(KIND=dp), POINTER :: gWork(:,:), HeatingEnergy(:), FrictionHeatFlux(:), &
      PressureHeatFlux(:), HeatingW(:)
  LOGICAL :: GradP, LateralStrain, GotAc, SurfAc, UsePrevGap, GotGrav, GotHeight, &
      CalcFrictionHeating, CalcPressureHeating, CalcHeating, UseHeating, DoneWeight = .FALSE.
  TYPE(Variable_t), POINTER :: pVar, thisVar, hVar
  INTEGER :: GapDirection, FrictionModel 
  REAL(KIND=dp) :: GapFactor, Nm, TotHeating
  CHARACTER(:), ALLOCATABLE :: str
  CHARACTER(*), PARAMETER :: Caller = 'FilmFlowSolver'
  LOGICAL :: Debug, FirstRound=.TRUE.
  
  SAVE STIFF, MASS, LOAD, FORCE, rho, ac, gap, gap0, mu, height, AcPres, Velocity, &
      AcPrevPressure, AllocationsDone, pVar, GotAc, SurfAC, FsiRhs, PrevGap, &
      FrictionModel, time0, hVar, HeatingEnergy, FrictionHeatFlux, PressureHeatFlux, &
      HeatingW, FirstRound
!------------------------------------------------------------------------------

  CALL Info(Caller,'Computing reduced dimensional Navier-Stokes equations!')

  Mesh => GetMesh()
  Element => GetActiveElement(1)

  dim = CoordinateSystemDimension()
    
  Params => GetSolverParams()
  thisVar => Solver % Variable
  
  mdim = ListGetInteger( Params,'Model Dimension',UnFoundFatal=.TRUE.)
  Convect = GetLogical( Params, 'Convect', Found )
  GradP = GetLogical( Params, 'GradP Discretization', Found ) 
  LateralStrain = GetLogical( Params,'Lateral Strain',Found )
  mingap = ListGetCReal( Params,'Min Gap Height',Found )
  IF(.NOT. Found) mingap = 1.0e-20
  GotAC = ListCheckPresentAnyMaterial( Model,'Artificial Compressibility')

  UsePrevGap = ListGetLogical( Params,'Use Gap Average',Found )
  
  CoupledIter = GetCoupledIter()
  
  GapDirection = 0
  GapFactor = ListGetCReal( Params,'Gap Addition Factor',Found )
  IF( Found ) THEN  
    GapDirection = mdim+1
    IF( ABS( GapFactor ) > 1.0_dp ) THEN
      CALL Warn(Caller,'"Gap Addition Factor" greater to unity does not make sense!')
    END IF
  END IF

  CSymmetry = ListGetLogical( Params,'Axi Symmetric',Found )
  IF(.NOT. Found ) THEN 
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) 
  END IF

  FrictionModel = 0
  str = ListGetString( Params,'Friction model',Found )
  IF(Found ) THEN
    IF(str == 'laminar' ) THEN
      FrictionModel = 0
    ELSE IF( str == 'darcy' ) THEN
      FrictionModel = 1
    ELSE IF( str == 'darcy2' ) THEN
      FrictionModel = 2
    ELSE IF( str == 'manning' ) THEN
      FrictionModel = 3
    ELSE IF (str == 'manning2' ) THEN
      FrictionModel = 4
    ELSE
      CALL Fatal(Caller,'Uknown friction model: '//TRIM(str))
    END IF
    CALL Info(Caller,'Using friction model: '//TRIM(str),Level=7)
  END IF
  
  
  grav = 0.0_dp
  gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',GotGrav)
  IF(GotGrav) THEN
    grav = ABS(gWork(SIZE(gWork,1),1))
  END IF

  IF( FrictionModel == 2 ) THEN
    IF(.NOT. GotGrav) CALL Fatal(Caller,'Manning equation not possible without gravity!')
    IF(CSymmetry) CALL Fatal(Caller,'Manning equation not applicable to axial symmetry!')
  END IF

  CalcFrictionHeating = ListGetLogical(Params,'Calculate Friction Heating',Found)
  CalcPressureHeating = ListGetLogical(Params,'Calculate Pressure Heating',Found)
  CalcHeating = CalcFrictionHeating .OR. CalcPressureHeating

  UseHeating = ListGetLogical(Params,'Use Heating Source',Found)
    
  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    CALL Info(Caller,'Dimension of Navier-Stokes equation: '//I2S(mdim))
    CALL Info(Caller,'Dimension of coordinate system: '//I2S(dim))

    n = (mdim+1)*(Mesh % MaxElementDOFs+4*BDOFs)  ! just big enough for elemental arrays
    ALLOCATE( FORCE(n), LOAD(mdim+4,n), STIFF(n,n), MASS(n,n), &
        rho(n), ac(n), gap(n), gap0(n), height(n), mu(n), AcPres(n), Velocity(mdim+1,n), STAT=istat )
    Velocity = 0.0_dp
    IF ( istat /= 0 ) THEN
      CALL Fatal( Caller, 'Memory allocation error.' )
    END IF
    IF( GradP ) THEN
      CALL Info(Caller,'"Gradp Discretization" is set True',Level=10)
    END IF     

    pVar => VariableGet( Mesh % Variables,'FilmPressure')
    IF ( .NOT. ASSOCIATED(pVar) ) THEN
      CALL Fatal( Caller, 'Could not find required field "FilmPressure"!')
    END IF
    n = SIZE(pVar % Values) 
    IF( GotAC ) THEN
      CALL Info(Caller,'Using artificial compressibility for FSI emulation!') 
      ALLOCATE(AcPrevPressure(n),FsiRhs(2,n))
      AcPrevPressure = 0.0_dp
      FsiRhs = 0.0_dp
      SurfAc = ListGetLogical( Params,'Surface Compressibility',Found )
    END IF
    IF(UsePrevGap) THEN
      ALLOCATE(PrevGap(n))
      PrevGap = 0.0_dp
    END IF

    IF(CalcFrictionHeating) THEN
      CALL Info(Caller,'Registering friction heat flux',Level=7)
      hVar => VariableGet( Mesh % Variables,'FilmFriction HeatFlux')
      FrictionHeatFlux => hVar % Values
    END IF
    IF(CalcPressureHeating) THEN
      CALL Info(Caller,'Registering convected heat flux',Level=7)
      hVar => VariableGet( Mesh % Variables,'FilmPressure HeatFlux')
      PressureHeatFlux => hVar % Values
    END IF
    IF(CalcHeating) THEN    
      CALL Info(Caller,'Registering total energy production',Level=7)
      hVar => VariableGet( Mesh % Variables,'Heating Energy')
      HeatingEnergy => hVar % Values
      ALLOCATE(HeatingW(SIZE(HeatingEnergy)))
    END IF
          
    AllocationsDone = .TRUE.
  END IF

  IF( CalcHeating) THEN
    ! If we are visiting the same timestep several times only compute the nodal heat flux once.
    time = GetTime()
    IF(ABS(time-time0) < TINY(time)) THEN
      IF (CalcFrictionHeating .AND. CalcPressureHeating) THEN        
        HeatingEnergy = HeatingEnergy - dt * MAX(FrictionHeatFlux  + PressureHeatFlux, 0.0_dp)
      ELSE
        IF(CalcFrictionHeating) THEN
          HeatingEnergy = HeatingEnergy - dt * FrictionHeatFlux 
        END IF
        IF( CalcPressureHeating ) THEN
          HeatingEnergy = HeatingEnergy - dt * PressureHeatFlux
        END IF
      END IF
    END IF
    IF( CalcFrictionHeating ) FrictionHeatFlux = 0.0_dp
    IF( CalcPressureHeating ) PressureHeatFlux = 0.0_dp
    IF(.NOT. DoneWeight) HeatingW = 0.0_dp
    time0 = time
  END IF
      
  IF(GotAc) THEN  
    ! When we do more than one nonlinear iteration the pressure used for FSI iteration
    ! differs from the current pressure. Hence we memorize the pressure at the start.
    AcPrevPressure = pVar % Values
    FsiRhs = 0.0_dp
  END IF
   
  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Because the heating uses previous values of velocity we need at least two iterations!
  IF (CalcHeating) maxiter = MAX(2,maxiter)

  
  DO iter=1,maxiter    
    !Initialize the system and do the assembly:
    !----------------
    CALL DefaultInitialize()

    Newton = GetNewtonActive()
    
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      ! Volume forces:
      !---------------
      BodyForce => GetBodyForce()
      LOAD = 0.0d0
      IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1,1:n) = GetReal( BodyForce, 'FilmFlow Bodyforce 1', Found )
        IF(mdim>1) Load(2,1:n) = GetReal( BodyForce, 'FilmFlow Bodyforce 2', Found )
        Load(mdim+1,1:n) = GetReal( BodyForce, 'Normal Velocity', Found )
        Load(mdim+2,1:n) = GetReal( BodyForce, 'Fsi Velocity', Found )

        ! We are slightly misusing "Load" here to store these quantities. 
        Load(mdim+3,1:n) = GetReal( BodyForce, 'Flow Admittance', Found)
        Load(mdim+4,1:n) = GetReal( BodyForce, 'External FilmPressure', Found)
      END IF

      ! Material parameters:
      !---------------------
      Material => GetMaterial()
      rho(1:n) = GetReal( Material, 'Density' )
      mu(1:n)  = GetReal( Material, 'Viscosity' )
      gap(1:n) = GetReal( Material, 'Gap Height' )

      height(1:n) = GetReal( Material,'Bedrock Height',GotHeight) 
      
      IF(FrictionModel == 1 .OR. FrictionModel == 2 ) THEN
        nm = ListGetCReal( Material,'Darcy Roughness',UnfoundFatal=.TRUE.)
      ELSE IF(FrictionModel == 3 ) THEN
        nm = ListGetCReal( Material,'Manning coefficient',UnfoundFatal=.TRUE.)
      END IF

      IF(UseHeating) THEN
        MeltHeat = ListGetConstReal( Material, 'Latent Heat', UnFoundFatal=.TRUE.)         
      END IF

      IF(CalcPressureHeating) THEN
        Cp = ListGetConstReal( Material,'Water Heat Capacity', UnfoundFatal=.TRUE.) 
        Ct = ListGetConstReal( Material,'Pressure Melting Coefficient', UnfoundFatal=.TRUE.) 
      END IF
      
      WHERE(gap(1:n) < mingap )
        gap(1:n) = mingap
      END WHERE

      IF(UsePrevGap) THEN
        IF(CoupledIter == 1) THEN
          PrevGap(pVar % Perm(Element % NodeIndexes)) = gap(1:n)
        END IF
        gap0(1:n) = PrevGap(pVar % Perm(Element % NodeIndexes))
      END IF
              
      IF( GotAC ) THEN
        ac(1:n) = GetReal( Material,'Artificial Compressibility',Found )
        ! This is not the latest pressure but a pressure from the previous FSI iteration.
        AcPres(1:n) = AcPrevPressure(pVar % Perm(Element % NodeIndexes))
      END IF

      ! Get previous elementwise velocity iterate:
      ! Note: pressure is the dim+1 component here!
      !-------------------------------------------
      CALL GetVectorLocalSolution( Velocity )
        
      ! Get element local matrix and rhs vector:
      !-----------------------------------------
      CALL LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, rho, gap, gap0, height, &
          mu, ac, Velocity, AcPres, Element, n, nd, nd+nb, &
          dim, mdim, FirstRound )
      
      IF ( nb>0 ) THEN
        CALL NSCondensate( nd, nb, mdim, STIFF, FORCE )
      END IF

      IF ( Transient ) THEN
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
      END IF
      
      ! Update global matrix and rhs vector from local matrix & vector:
      !----------------------------------------------------------------
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    CALL DefaultFinishBulkAssembly()


    IF( GotAC ) THEN
      BLOCK
        REAL(KIND=dp) :: sorig, sfsi, coeff
        sorig = SUM(FsiRhs(1,:))
        sfsi = SUM(FsiRhs(2,:))
        coeff = 1.0_dp
        IF(sfsi /= 0.0) coeff = sorig / sfsi            
        IF(sfsi < sorig) coeff = 1.0_dp

        ! Just report incoming and outgoing total fluxes
        IF(InfoActive(20)) THEN
          PRINT *,'RHSComp:',sorig,sfsi,coeff
        END IF
      END BLOCK
    END IF
      
    DO t=1, Solver % Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE

      n = GetElementNOFNodes()      
      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE

      rho(1:n) = GetParentMatProp( 'Density', Element, Found )
      mu(1:n)  = GetParentMatProp( 'Viscosity', Element, Found )
      gap(1:n) = GetParentMatProp( 'Gap Height', Element, Found )

      WHERE(gap(1:n) < mingap )
        gap(1:n) = mingap
      END WHERE

      IF(UsePrevGap) THEN
        gap0(1:n) = PrevGap(pVar % Perm(Element % NodeIndexes))
      END IF
      
      DO i=1,mdim
        Load(i,1:n) = GetReal( BC, 'FilmPressure '//I2S(i), Found ) 
      END DO
      Load(mdim+1,1:n) = GetReal( BC, 'Mass Flux', Found )
      
      CALL LocalBoundaryMatrix(  MASS, STIFF, FORCE, Load, rho, gap, mu, &
          Element, n, dim, mdim )

      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    
    Norm = DefaultSolve()
    
    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
    FirstRound = .FALSE.
  END DO

  CALL DefaultFinish()

  IF( CalcHeating ) THEN   
    BLOCK
      REAL(KIND=dp) :: TotFlux, Area, cfix

      TotFlux = 0.0_dp
      IF(CalcFrictionHeating) &
          TotFlux = ParallelReduction(SUM(FrictionHeatFlux))
      IF(CalcPressureHeating) &
          TotFlux = TotFlux + ParallelReduction(SUM(FrictionHeatFlux))
      IF ( ParEnv % PEs > 1) THEN
        ! In parallel we need to sum up the terms at shared nodes. 
        IF( CalcFrictionHeating ) &
            CALL ParallelSumNodalVector( Mesh, FrictionHeatFlux, HVar % Perm )
        IF( CalcPressureHeating ) &
            CALL ParallelSumNodalVector( Mesh, PressureHeatFlux, HVar % Perm )
        IF(.NOT. DoneWeight ) &
          CALL ParallelSumNodalVector( Mesh, HeatingW, HVar % Perm )
      END IF
      DoneWeight = .TRUE.

      IF(CalcFrictionHeating) THEN
        WHERE(HeatingW > 0.0_dp )
          FrictionHeatFlux = FrictionHeatFlux / HeatingW
        END WHERE
      END IF

      IF(CalcPressureHeating) THEN
        WHERE(HeatingW > 0.0_dp )
          PressureHeatFlux = PressureHeatFlux / HeatingW
        END WHERE
      END IF
        
      ! This is just a tentative feature that would allow finding of steady state
      ! solutions. Not usuful generally. 
      IF(ListGetLogical(Params,'Enforce Zero Heating', Found ) ) THEN
        Area = ParallelReduction( SUM( HeatingW ) )
        cfix = -TotFlux / Area
        FrictionHeatFlux = FrictionHeatFlux + cfix
      END IF

      ! Total dissipated heat over time.
      IF (CalcFrictionHeating .AND. CalcPressureHeating) THEN        
        HeatingEnergy = HeatingEnergy + dt * MAX(FrictionHeatFlux  + PressureHeatFlux, 0.0_dp)
      ELSE
        IF(CalcFrictionHeating) THEN
          HeatingEnergy = HeatingEnergy + dt * FrictionHeatFlux 
        END IF
        IF( CalcPressureHeating ) THEN
          HeatingEnergy = HeatingEnergy + dt * PressureHeatFlux
        END IF
      END IF
    
      
        
      WRITE(Message,'(A,ES12.5)') 'Total heating power: ',TotFlux
      CALL Info(Caller, Message, Level=7)
    END BLOCK
  END IF
  

  BLOCK
    REAL(KIND=dp), POINTER :: Comp(:)

    IF( InfoActive(12) ) THEN
      n = SIZE(pVar % Values)
      DO i=1, Solver % Variable % dofs
        Comp => Solver % Variable % Values(i::Solver % Variable % Dofs)      
        CALL VectorValuesRange(Comp,n,'Velocity '//I2S(i))       
      END DO
      CALL VectorValuesRange(pVar % Values,n,'Pressure')       
      IF(GotAc) THEN
        CALL VectorValuesRange(AcPrevPressure,n,'Pressure0')       
        CALL VectorValuesRange(pVar % Values - AcPrevPressure,n,'PressureDiff '//I2S(CoupledIter))       
      END IF
    END IF      
  END BLOCK
  
  CALL Info(Caller,'All done',Level=12)

  
CONTAINS

  ! Eq (29) in:
  ! P. Praks and D. Brkić, Review of new flow friction equations:
  ! Constructing Colebrook’s explicit correlations accurately, Rev. int. métodos numér. cálc. diseño ing. (2020).
  ! Vol. 36, (3), 41 URL https://www.scipedia.com/public/Praks_Brkic_2020a
  ! https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae
  !--------------------------------------------------------------------------------------
  FUNCTION FrictionLawPraks(v,rho,nu,D,eps) RESULT (f)
    REAL(KIND=dp) :: v, rho, nu, D, eps, f
    REAL(KIND=dp) :: Re, A, B, C, x
    REAL(KIND=dp) :: MinRe
    LOGICAL :: Visited = .FALSE.

    SAVE Visited, MinRe

    IF(.NOT. Visited) THEN
      MinRe = ListGetCReal( Params,'Min Reynolds Number',Found )
      IF(.NOT. Found) MinRe = 4000.0_dp
      Visited = .TRUE.
    END IF

    ! The division by 2 fixes the inconsistancy between two scientific communities.
    Re = v*(D/2)*rho/nu

    IF(Re < MinRe) THEN
      Re = MinRe
      ! Enforce also the speed to be compatible with the min Re number!
      ! Note: this has effect also outside this routine!
      v = Re * nu / ((D/2) * rho) 
    END IF      
    
    A = Re * eps / 8.0897
    B = LOG(Re) - 0.779626
    x = A+B
    C = LOG(x)

    f = (0.8685972*(B-C+C/(x-0.5588*C+1.2079)))**(-0.5_dp)
    
  END FUNCTION FrictionLawPraks
    
  
!------------------------------------------------------------------------------
  SUBROUTINE LocalBulkMatrix(  MASS, STIFF, FORCE, LOAD, Nodalrho, NodalGap, &
      NodalGap0, NodalH, Nodalmu, NodalAC, NodalVelo, NodalAcPres, &
      Element, n, nd, ntot, dim, mdim , FirstRound)
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), NodalAC(:), Nodalrho(:), &
        NodalGap(:), NodalGap0(:), NodalH(:), NodalAcPres(:), NodalVelo(:,:)
    INTEGER :: dim, mdim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: FirstRound
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3)
    REAL(KIND=dp) :: DetJ,LoadAtIP(mdim+4),Velo(mdim), VeloGrad(mdim,mdim), gapGrad(mdim), &
        hGrad(mdim), presGrad(mdim)
    REAL(KIND=dp) :: NodalPres(n), NodalPrevPres(n)
    REAL(KIND=dp), POINTER :: A(:,:),F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, AcPres, gap, gap0, gap2, gapi, &
        ac, s, s0, s1, MinPres, MuCoeff, MinSpeed, Speed, h, q_p, q_f, Pres, &
        PrevPres, FlowAdm
    LOGICAL :: Visited = .FALSE.
    
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes, Visited, MinPres, MinSpeed
!------------------------------------------------------------------------------

    
    CALL GetElementNodes( Nodes )

    IF( GapDirection > 0 ) THEN
      SELECT CASE( GapDirection )
      CASE(1)
        Nodes % x(1:n) = Nodes % x(1:n) + GapFactor * NodalGap(1:n)
      CASE(2)
        Nodes % y(1:n) = Nodes % y(1:n) + GapFactor * NodalGap(1:n)
      CASE(3)
        Nodes % z(1:n) = Nodes % z(1:n) + GapFactor * NodalGap(1:n)
      END SELECT
      ! Does this have an effect?
      !PRINT *,'GapFactor:',GapFactor * NodalGap(1:n)
    END IF

    CALL GetLocalSolution( NodalPres,UElement=Element,UVariable=pVar)
    CALL GetLocalSolution( NodalPrevPres,UElement=Element,UVariable=pVar,tStep=-1)
    
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0
    gapGrad = 0.0_dp
    hGrad = 0.0_dp
    presGrad = 0.0_dp
    q_f = 0.0_dp
    q_p = 0.0_dp
    
    ! To my understanding we want to include the gap height to weight
    IF( Csymmetry ) THEN
      geomc = 2
    ELSE
      geomc = 1
    END IF
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element, PReferenceElement = .TRUE. )

    !IP = GaussPoints( Element )

    IF(.NOT. Visited) THEN
      CALL Info(Caller,'Number of integration points: '//I2S(IP % n))
      MinPres = ListGetConstReal( Params,'Min FilmPressure',Found )
      IF(.NOT. Found) MinPres = -HUGE(MinPres)
      MinSpeed = ListGetConstReal( Params,'Min Speed',Found )
      IF(.NOT. Found) MinSpeed = 1.0e-6
      Visited = .TRUE.
    END IF
    
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
         IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       s1 = s 
       s0 = s 
              
       ! Material parameters at the integration point:
       !----------------------------------------------      
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       gap = SUM( Basis(1:n) * NodalGap(1:n) ) 
       gap0 = SUM( Basis(1:n) * NodalGap0(1:n) ) 
       
       AcPres = SUM( NodalAcPres(1:n) * Basis(1:n) )
       AcPres = MAX(MinPres,AcPres)

       Pres = SUM(NodalPres(1:n) * Basis(1:n) )
       
       DO i=1,mdim
         gapGrad(i) = SUM( NodalGap(1:nd) * dBasisdx(1:nd,i) )
         presGrad(i) = SUM( NodalPres(1:nd) * dBasisdx(1:nd,i) )
       END DO

       IF( GotHeight ) THEN
         h = SUM( NodalH(1:nd) * Basis(1:nd) )
         DO i=1,mdim
           hGrad(i) = SUM( NodalH(1:nd) * dBasisdx(1:nd,i) ) 
         END DO
       ELSE
         IF(mdim == 1) THEN
           h = SUM( Nodes % y(1:n) * Basis(1:n) ) 
           hGrad(1) = SUM( Nodes % y(1:n) * dBasisdx(1:n,1) ) 
         ELSE
           h = SUM( Nodes % z(1:n) * Basis(1:n) ) 
           hGrad(1) = SUM( Nodes % z(1:n) * dBasisdx(1:n,1) ) 
           hGrad(2) = SUM( Nodes % z(1:n) * dBasisdx(1:n,2) ) 
         END IF
       END IF
         
       ! Previous velocity at the integration point:
       !--------------------------------------------
       Velo = MATMUL( NodalVelo(1:mdim,1:nd), Basis(1:nd) )
       VeloGrad = MATMUL( NodalVelo(1:mdim,1:nd), dBasisdx(1:nd,1:mdim) )
       Speed = SQRT(SUM(Velo(1:mdim)**2))
       
       IF( GotAC ) THEN
         ac = SUM( NodalAC(1:n) * Basis(1:n) ) / dt
         !IF(.NOT. SurfAc) ac = ac * gap
       END IF
       
       ! The source term at the integration point:
       !------------------------------------------
       DO i=1,mdim+4
         LoadAtIP(i) = SUM( Basis(1:n) * LOAD(i,1:n) )
       END DO

       IF ( Convect .AND. Newton ) THEN
         LoadAtIp(1:mdim) = LoadAtIp(1:mdim) + rho * MATMUL(VeloGrad(1:mdim,1:mdim),Velo(1:mdim))
       END IF

       ! Fsi velocity
       LoadAtIp(mdim+1:mdim+2) = geomc * LoadAtIp(mdim+1:mdim+2) 
             
       ! This is the Poisseille flow resistance
       IF(UsePrevGap) THEN
         ! This takes the analytical average when going from 1/d_0^2 to 1/d^2. 
         gap2 = gap*gap0
         gapi = SQRT(gap2)
       ELSE
         gap2 = gap**2
         gapi = gap
       END IF

       SELECT CASE( FrictionModel )
       CASE( 1, 2 ) 
         BLOCK
           REAL(KIND=dp) :: D, R, fd, GradZphi2
           Speed = MAX(MinSpeed,Speed)
           ! for a cross-section that is uniform along the tube or channel length, the wetted diameter is defined as Dh=4*A/P
           ! where A is the cross-sectional area of the flow, and P is the wetted perimeter of the cross-section, see
           ! https://en.wikipedia.org/wiki/Hydraulic_diameter
           ! this leads to consistent friction from the Colebrook–White equation whether Dh or Rh are used, see
           ! https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae             
           ! Note that for csummetry "gapi" is radius as the same formula works as well. 
           D = 2 * gapi
           fd = FrictionLawPraks(Speed,rho,mu,D,nm)
           IF( FrictionModel == 1 ) THEN
             MuCoeff = fd * rho * Speed / (2*D)
           ELSE
             GradZphi2 =  MAX(SUM((hGrad(1:mdim) + presGrad(1:mdim)/(rho*Grav))**2), 1.0E-09)
             MuCoeff = rho * SQRT(fd*Grav) * (2*gap)**(-1.0/2.0) * GradZphi2**(1.0/4.0)
           END IF
         END BLOCK

       CASE( 3 ) 
         BLOCK
           REAL(KIND=dp) :: GradZphi2
           GradZphi2 = MAX(SUM((hGrad(1:mdim) + presGrad(1:mdim)/(rho*Grav))**2), 1.0E-09)
           MuCoeff = nm * rho * Grav * (gapi/2)**(-2.0/3) * GradZphi2**(1.0/4.0)
         END BLOCK

       CASE( 4)
         BLOCK
           MuCoeff = rho * Grav * nm**2 * (gapi/2)**(-4.0/3) * Speed 
       END BLOCK
                    
       CASE DEFAULT 
         IF( CSymmetry ) THEN
           ! Note: gap is here the radius!
           MuCoeff = 8 * mu / gap2 
         ELSE
           MuCoeff = 12 * mu / gap2
         END IF

       END SELECT

       IF( CalcFrictionHeating ) THEN
         q_f = MuCoeff * gapi * Speed**2
       END IF
         
       IF( CalcPressureHeating ) THEN
         PrevPres = SUM(NodalPrevPres(1:n) * Basis(1:n))
         q_p = rho * gapi * Cp * Ct * ((Pres-PrevPres)/dt + SUM(Velo(1:mdim)*presGrad(1:mdim)))
         
       END IF
       IF (FirstRound) THEN
         q_f = 0.0_dp
         q_p = 0.0_dp
       END IF
       
       ! Finally, the elemental matrix & vector:
       !----------------------------------------       
       DO p=1,ntot
         DO q=1,ntot
           i = (mdim+1) * (p-1) + 1
           j = (mdim+1) * (q-1) + 1
           A => STIFF(i:i+mdim,j:j+mdim)
           M => MASS(i:i+mdim,j:j+mdim)

           DO i=1,mdim
             IF( Transient ) THEN
               M(i,i) = M(i,i) + s * rho * Basis(q) * Basis(p)
             END IF

             A(i,i) = A(i,i) + s * MuCoeff * Basis(q) * Basis(p)              

             DO j = 1,mdim
               IF( LateralStrain ) THEN 
                 A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
                 A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
               END IF
                 
               IF ( Convect ) THEN
                 A(i,i) = A(i,i) + s * rho * Velo(j) * dBasisdx(q,j) * Basis(p)
                 IF ( Newton ) THEN
                    A(i,j) = A(i,j) + s * rho * VeloGrad(i,j) * Basis(q) * Basis(p)
                 END IF
               END IF
             END DO
             
             ! Note that here the gap height must be included in the continuity equation
             IF( GradP ) THEN
               A(i,mdim+1) = A(i,mdim+1) + s * dBasisdx(q,i) * Basis(p)
               A(mdim+1,i) = A(mdim+1,i) - s * gap * rho * Basis(q) * dBasisdx(p,i)               
             ELSE
               A(i,mdim+1) = A(i,mdim+1) - s * Basis(q) * dBasisdx(p,i)
               A(mdim+1,i) = A(mdim+1,i) + s * gap * rho * dBasisdx(q,i) * Basis(p) & 
                   + geomc * s * rho * Basis(q) * gapGrad(i) * Basis(p)
             END IF
           END DO
             
           ! This is the implicit term in artificial compressibility for FSI coupling
           ! applied to thin film flow. 
           ! Div(u) + (c/dt)*p^(m) = (c/dt)*p^(m-1)
           ! See Raback et al., CFD Eccomas 2001.
           ! "FLUID-STRUCTURE INTERACTION BOUNDARY CONDITIONS BY ARTIFICIAL COMPRESSIBILITY".
           IF(GotAC) A(mdim+1,mdim+1) = A(mdim+1,mdim+1) + ac * s * rho * Basis(q) * Basis(p)              

           ! The implicit term for weakly enforce incoming flux. 
           A(mdim+1,mdim+1) = A(mdim+1,mdim+1) + s * rho * LoadAtIP(mdim+3) * Basis(q) * Basis(p)              
         END DO
         
         i = (mdim+1) * (p-1) + 1
         F => FORCE(i:i+mdim)
         
         ! This is the explit term in artificial compressibility for FSI coupling
         IF( GotAC ) F(mdim+1) = F(mdim+1) + ac * s * rho * Basis(p) * AcPres         

         ! Body force for velocity components and pressure
         F(1:mdim+1) = F(1:mdim+1) + s * rho * Basis(p) * LoadAtIp(1:mdim+1)

         ! Gravity for the slope
         IF(GotGrav) THEN
           F(1:mdim) = F(1:mdim) - s * rho * Grav * Basis(p) * hGrad(1:mdim)
         END IF
         
         ! Additional body force from FSI velocity
         F(mdim+1) = F(mdim+1) - s * rho * Basis(p) * LoadAtIp(mdim+2) 


         ! Robin condition for incoming flow in terms of (Flow Admittance) * (p - p_ext)
         F(mdim+1) = F(mdim+1) + s * rho * Basis(p) * LoadAtIp(mdim+3) * LoadAtIP(mdim+4) 

         
         IF(UseHeating) THEN
           ! Additional source term from friction melting.
           ! Continuity equation is weighted by gap so the BC term need not be divided by it
           ! as is the case for momentum equation.
           F(1:mdim+1) = F(1:mdim+1) + s * Basis(p) * MAX(q_f + q_p,0.0_dp) / (MeltHeat * rho)
         END IF
       END DO
         
       IF( CalcFrictionHeating ) THEN
         FrictionHeatFlux(hVar % Perm(Element % NodeIndexes)) = &
             FrictionHeatFlux(hVar % Perm(Element % NodeIndexes)) + s * Basis(1:n) * q_f
       END IF
       IF( CalcPressureHeating ) THEN
         PressureHeatFlux(hVar % Perm(Element % NodeIndexes)) = &
             PressureHeatFlux(hVar % Perm(Element % NodeIndexes)) + s * Basis(1:n) * q_p
       END IF
       IF(CalcHeating .AND. .NOT. DoneWeight) THEN
         ! We compute the normalization weight.
         HeatingW(hVar % Perm(Element % NodeIndexes)) = &
             HeatingW(hVar % Perm(Element % NodeIndexes)) + s * Basis(1:n)
       END IF
         
       ! These are just recorded in order to study the total forced
       ! and induced (by FSI coupling) fluxes. 
       IF(GotAC) THEN
         FsiRhs(1,ThisVar % Perm(Element % NodeIndexes)) = &
             FsiRhs(1,ThisVar % Perm(Element % NodeIndexes))  + &
             s * rho * LoadAtIp(mdim+1) * Basis(1:n)             
         
         FsiRhs(2,ThisVar % Perm(Element % NodeIndexes)) = &
             FsiRhs(2,ThisVar % Perm(Element % NodeIndexes))  + &
             s * rho * LoadAtIp(mdim+2) * Basis(1:n)             
       END IF
     END DO
     
   ! for p2/p1 elements set Dirichlet constraint for unused dofs,
   ! EliminateDirichlet will get rid of these:
   !-------------------------------------------------------------
    DO p = n+1,ntot
      i = (mdim+1) * p
      FORCE(i)   = 0.0d0
      MASS(:,i)  = 0.0d0
      MASS(i,:)  = 0.0d0
      STIFF(i,:) = 0.0d0
      STIFF(:,i) = 0.0d0
      STIFF(i,i) = 1.0d0
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalBulkMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalBoundaryMatrix(  MASS, STIFF, FORCE, NodalLoad, Nodalrho, NodalGap, &
      Nodalmu, Element, n, dim, mdim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), NodalLoad(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalGap(:)
    INTEGER :: dim, mdim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,Load(mdim+1)
    REAL(KIND=dp), POINTER :: F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, p, q, geomc
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, gap, ac, s
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0
    
    ! To my understanding we want to include the gap height to weight
    IF( Csymmetry ) THEN
      geomc = 2
    ELSE
      geomc = 1
    END IF
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
         IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ
              
       ! Material parameters at the integration point:
       !----------------------------------------------      
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )
       gap = SUM( Basis(1:n) * NodalGap(1:n) ) 

       DO i=1,mdim+1
         Load(i) = SUM( Basis(1:n) * NodalLoad(i,1:n) ) 
       END DO
         
       ! Finally, the elemental matrix & vector:
       !----------------------------------------       
       DO p=1,n
         i = (mdim+1) * (p-1) + 1
         F => FORCE(i:i+mdim)
         F = F + s * Basis(p) * Load
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalBoundaryMatrix
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
END SUBROUTINE FilmFlowSolver
!------------------------------------------------------------------------------
