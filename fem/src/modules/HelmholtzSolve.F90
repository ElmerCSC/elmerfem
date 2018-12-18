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
! *  Original Date: 04 Oct 2000
! *
! ****************************************************************************/

SUBROUTINE HelmholtzSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()
  CALL ListAddNewLogical( Params,'Linear System Complex',.TRUE.)

END SUBROUTINE HelmholtzSolver_init
  
  

!------------------------------------------------------------------------------
!> Solver for Helmholtz equation accounting also for variable density and 
!> convection field. Also includes a built-in interface for coupling with harmonic
!> velocity or displacement fields at the boundary. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE HelmholtzSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
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

  LOGICAL :: AllocationsDone = .FALSE., Bubbles, Found, UseConvection, &
      UseDensity, PlaneWave

  INTEGER :: iter, i, j, k, l, n, nd, t, istat, eq, LocalNodes
  REAL(KIND=dp) :: Norm, AngularFrequency, s


  TYPE(ValueList_t), POINTER :: Equation, Material, BodyForce, &
             BC, SolverParams, Simulation

  INTEGER :: NonlinearIter

  REAL(KIND=dp), ALLOCATABLE :: Load(:,:), Work(:), &
       SoundSpeed(:), Density(:), Damping(:), Impedance(:,:), ConvVelo(:,:)

  COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)

  SAVE STIFF, Work, Load, FORCE, &
       SoundSpeed, Density, Damping, Impedance, AllocationsDone, ConvVelo

#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1
#else
  REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
#endif

!-----------------------------------------------------------------------------
! Local variables for performing analyses with harmonic interfaces
!-----------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: FlowSol, DispSol
  LOGICAL :: FlowInterface, StructureInterface, stat 
  LOGICAL :: AnyFlowInterface, AnyStructureInterface, GotFrequency
  REAL(KIND=dp), POINTER :: Flow(:), Disp(:)
  COMPLEX(KIND=dp), POINTER :: DispEigen(:)
  INTEGER, POINTER ::  FlowPerm(:), DispPerm(:), PresPerm(:)
  INTEGER :: dim, FlowDofs, DispDofs, NoEigen
  TYPE(Element_t),POINTER :: Parent
  COMPLEX(KIND=dp), ALLOCATABLE :: WallVelocity(:,:)
  COMPLEX(KIND=dp) :: ImUnit
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName

  SAVE WallVelocity

!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     N = Solver % Mesh % MaxElementDOFs

     IF ( AllocationsDone ) THEN
        DEALLOCATE(            &
             Impedance,        &
             Work,             &
             FORCE,  STIFF,    &
             SoundSpeed, Density, ConvVelo, Damping, Load, &
             WallVelocity )
     END IF

     ALLOCATE( &
          Impedance( 2,N ),    &
          Work( N ),           &
          FORCE( 2*N ),        &
          STIFF( 2*N,2*N ),    &
          SoundSpeed( N ), Density( N ), ConvVelo(3,N), Damping( N ), Load( 2,N ), &
          WallVelocity(3,N), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'HelmholzSolve', 'Memory allocation error.' )
     END IF

     AllocationsDone = .TRUE.
  END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Do some additional initialization, and go for it
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  NonlinearIter = GetInteger( SolverParams, &
       'Nonlinear System Max Iterations', Found )

  IF ( .NOT.Found ) NonlinearIter = 1
  Bubbles = GetLogical( SolverParams, 'Bubbles', Found )

! Initially density was not used in the Helmholtz equation. However, if there 
! are several different densities it must be used and hence it was added later.
! For backward compatibility also the version where density is not provided should work.
!---------------------------------------------------------------------------------------
  UseDensity = .TRUE.
  DO i=1,Model % NumberOfMaterials
    IF( .NOT. ListCheckPresent( Model % Materials(i) % Values,'Density') ) THEN
      UseDensity = .FALSE.
      IF( Model % NumberOfMaterials > 1) THEN
        CALL Warn('HelmholtzSolver','There are multiple materials and density is not used, might result to errors!')
      END IF
      EXIT
    END IF
  END DO

  ! This flag could be needed in FSI iterations, for example
  CALL ListAddLogical( SolverParams,'Use Density', UseDensity )
  
  n = GetElementNOFNodes()
  Simulation => GetSimulation()
  dim = CoordinateSystemDimension()     
  GotFrequency = .FALSE.

  ! Check for flow or strcuture interface
  !--------------------------------------------------------
  WallVelocity = 0.0_dp
  ImUnit = CMPLX(0.0d0,1.0d0,KIND=dp) 
  PresPerm => Solver % Variable % Perm

  AnyFlowInterface = ListCheckPresentAnyBC(Model,'Flow Interface')
  AnyStructureInterface = ListCheckPresentAnyBC(Model,'Structure Interface')

  IF( AnyFlowInterface ) THEN
    VarName = GetString( SolverParams,'Velocity Variable Name',Found)
    IF(.NOT. Found ) VarName = 'flow'
    FlowSol => VariableGet( Solver % Mesh % Variables,VarName )
    IF ( ASSOCIATED(FlowSol) ) THEN
      Flow => FlowSol % Values
      FlowPerm => FlowSol % Perm
      FlowDofs = FlowSol % Dofs
    ELSE
      CALL Fatal('HelmholtzSolver','No flow variable associated:'//TRIM(VarName))
    END IF
    IF( FlowDofs /= 2*dim ) THEN
      CALL Fatal('HelmholtzSolver','Harmonic flow field should have 2*dim components')
    END IF
  END IF

  IF( AnyStructureInterface ) THEN
    VarName = GetString( SolverParams,'Displacement Variable Name',Found)
    IF(.NOT. Found ) VarName = 'displacement'
    DispSol => VariableGet( Solver % Mesh % Variables,VarName )
    IF ( ASSOCIATED(DispSol) ) THEN
      Disp => DispSol % Values
      DispPerm => DispSol % Perm
      DispDofs = DispSol % Dofs
    ELSE
      CALL Fatal('HelmholtzSolver','No displacement variable associated:'//TRIM(VarName))
    END IF

    NoEigen = GetInteger( SolverParams,'Displacement Variable EigenMode',Found)
    IF( NoEigen > 0 ) THEN
      IF( NoEigen > SIZE(DispSol % EigenValues)) THEN
        CALL Fatal('HelmholtzSolver','Requested eigenmode does not exist')
      END IF
      DispEigen => DispSol % EigenVectors(NoEigen,:)
      IF( GetLogical( SolverParams,'Displacement Variable Frequency',Found) ) THEN
        AngularFrequency = SQRT( DispSol % EigenValues(NoEigen))
	GotFrequency = .TRUE.
      END IF
      IF( DispDofs < dim ) THEN
        CALL Fatal('HelmholtzSolver','Eigenmode displacement field should have at least 1*dim components')
      END IF
    ELSE
      IF( DispDofs /= 2*dim ) THEN
        CALL Fatal('HelmholtzSolver','Harmonic displacement field should have 2*dim components')
      END IF
    END IF
  END IF


  ! Figure out angular frequency and save it for SaveScalars:
  !----------------------------------------------------------
  IF(.NOT. GotFrequency ) THEN 
    AngularFrequency = GetAngularFrequency(Found = GotFrequency )
  END IF

  IF(.NOT. GotFrequency ) THEN
    CALL Fatal('HelmholtzSolver','Could not fugure out Frquency!')
  END IF


  CALL ListAddConstReal( Model % Simulation, 'res: Frequency', AngularFrequency /(2*PI) )


  ! Check whether the equation lives on a convection field
  !-------------------------------------------------------
  UseConvection = .FALSE.
  ConvVelo = 0.0_dp

  DO i = 1, Model % NumberOfMaterials
    Material => Model % Materials(i) % Values

    UseConvection = ListCheckPresent( Material,'Convection Velocity 1')

    Found = ListCheckPresent( Material,'Convection Velocity 2')
    UseConvection = UseConvection .OR. Found

    IF( dim == 3 ) THEN
      Found = ListCheckPresent( Material,'Convection Velocity 3')
      UseConvection = UseConvection .OR. Found
    END IF
    IF( UseConvection ) EXIT
  END DO

!------------------------------------------------------------------------------
! Iterate over any nonlinearity of material or source
!------------------------------------------------------------------------------
  Norm = Solver % Variable % Norm
  totat = 0.0d0
  totst = 0.0d0

  CALL DefaultStart()
  
  DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
     at  = CPUTime()
     at0 = RealTime()

     CALL Info( 'HelmholtzSolve', ' ', Level=4 )
     CALL Info( 'HelmholtzSolve', '-------------------------------------', Level=4 )
     WRITE( Message, * ) 'Helmholtz iteration', iter
     CALL Info( 'HelmholtzSolve', Message, Level=4 )
     WRITE( Message, * ) 'Frequency (Hz): ', AngularFrequency/(2*PI)
     CALL Info( 'HelmholtzSolve', Message, Level=4 )
     CALL Info( 'HelmholtzSolve', '-------------------------------------', Level=4 )
     CALL Info( 'HelmholtzSolve', ' ', Level=4 )
     CALL Info( 'HelmholtzSolve', 'Starting Assembly', Level=4 )

     CALL DefaultInitialize()
!
!    Do the bulk assembly:
!    ---------------------

     CALL StartAdvanceOutput('HelmholtzSolve', 'Assembly:' )
!------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
!------------------------------------------------------------------------------

        CALL AdvanceOutput( t,Solver % NumberOFActiveElements )

!------------------------------------------------------------------------------
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()

!       Get equation & material parameters:
!       -----------------------------------
        Equation => GetEquation()
        Material => GetMaterial()

        Damping(1:n)    = GetReal( Material, 'Sound Damping', Found )
        SoundSpeed(1:n) = GetReal( Material, 'Sound Speed', Found )

        IF( UseDensity) Density(1:n) = GetReal( Material, 'Density', Found )

        IF( UseConvection ) THEN
          ConvVelo(1,1:n) = GetReal( Material, 'Convection Velocity 1', Found )
          ConvVelo(2,1:n) = GetReal( Material, 'Convection Velocity 2', Found )
          IF( dim == 3 ) THEN
            ConvVelo(3,1:n) = GetReal( Material, 'Convection Velocity 3', Found )
          END IF
        END IF

!       The source term on nodes:
!       -------------------------
        BodyForce => GetBodyForce()
        Load(1,1:n) = GetReal( BodyForce, 'Pressure Source 1', Found )
        Load(2,1:n) = GetReal( BodyForce, 'Pressure Source 2', Found )

!       Get element local matrix and rhs vector:
!       ----------------------------------------
        CALL LocalMatrix(  STIFF, FORCE, AngularFrequency, &
           SoundSpeed, ConvVelo, Damping, Load, Bubbles, Element, n, nd )

!       Update global matrix and rhs vector from local matrix & vector:
!       ---------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END DO

     CALL DefaultFinishBulkAssembly()

!---------------------------------------------------------------------------
! Standard Neumann & Newton BCs:
!---------------------------------------------------------------------------

     DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT.ActiveBoundaryElement() ) CYCLE

        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()

        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
          Load(1,1:n) = GetReal( BC, 'Wave Flux 1', Found )
          Load(2,1:n) = GetReal( BC, 'Wave Flux 2', Found )
          Impedance(1,1:n) = GetReal( BC, 'Wave Impedance 1', Found )
          Impedance(2,1:n) = GetReal( BC, 'Wave Impedance 2', Found )
                 
          IF( UseDensity ) THEN
            Density(1:n) = GetParentMatProp( 'Density', Element )
          END IF

          PlaneWave = GetLogical( BC,'Plane Wave BC',Found )
          IF( PlaneWave ) THEN 
            Impedance(1,1:n) = GetParentMatProp( 'Sound Speed', Element )
            Impedance(2,1:n) = 0.0_dp
          END IF

	  IF( UseConvection ) THEN
            ConvVelo(1,1:n) = GetParentMatProp( 'Convection Velocity 1', Element, Found )
            ConvVelo(2,1:n) = GetParentMatProp( 'Convection Velocity 2', Element, Found )
            IF( dim == 3 ) THEN
              ConvVelo(3,1:n) = GetParentMatProp( 'Convection Velocity 3', Element, Found )
            END IF
          END IF
          

          CALL LocalMatrixBoundary(  STIFF, FORCE, AngularFrequency, &
              Impedance, Load, Element, n, nd, ConvVelo )

          CALL DefaultUpdateEquations( STIFF, FORCE )
        END IF
     END DO

!-----------------------------------------------------------------------------
! Boundary conditions on harmonic flow or structure interfaces
!-----------------------------------------------------------------------------
     
     IF( AnyStructureInterface .OR. AnyFlowInterface ) THEN
       DO t=1, Solver % Mesh % NumberOfBoundaryElements
         Element => GetBoundaryElement(t)
         IF ( .NOT.ActiveBoundaryElement() ) CYCLE
         n = GetElementNOFNodes()
         IF ( GetElementFamily() == 1 ) CYCLE
         BC => GetBC()
         IF( .NOT. ASSOCIATED( BC ) ) CYCLE
         
         FlowInterface = ListGetLogical( BC, 'Flow Interface', Found )
         StructureInterface = ListGetLogical( BC, 'Structure Interface', Found )
         
         IF ( .NOT. (FlowInterface .OR. StructureInterface) ) CYCLE
         
         IF( FlowInterface ) THEN
           IF ( ANY( FlowPerm( Element % NodeIndexes(1:n) ) == 0 ) ) THEN
             CALL Fatal( 'HelmholtzSolve', 'Flow solution is not available on boundary')
           END IF           
           DO j=1,n
             k = FlowPerm( Element % NodeIndexes(j) ) 
             DO l=1,dim
               WallVelocity(l,j) = Flow( (k-1)*FlowDofs + 2*l-1 ) + &
                   ImUnit * Flow( (k-1)*FlowDofs + 2*l ) 
             END DO
           END DO
         ELSE IF( StructureInterface ) THEN
           IF ( ANY( DispPerm( Element % NodeIndexes(1:n) ) == 0 ) ) THEN
             CALL Fatal( 'HelmholtzSolve', 'Displacement solution is not available on boundary')
           END IF           
           IF( NoEigen > 0 ) THEN
             DO j=1,n
               k = DispPerm( Element % NodeIndexes(j) ) 
               DO l=1,dim
                 WallVelocity(l,j) = DispEigen( (k-1)*DispDofs + l )
               END DO
             END DO          
           ELSE
             DO j=1,n
               k = DispPerm( Element % NodeIndexes(j) ) 
               DO l=1,dim
                 WallVelocity(l,j) = Disp( (k-1)*DispDofs + 2*l-1 ) + &
                     ImUnit * Disp( (k-1)*DispDofs + 2*l )
               END DO
             END DO
           END IF
           WallVelocity = ImUnit * AngularFrequency * WallVelocity
         END IF
         
         ! Find the Helmholtz parent to determine the reference density.
         ! If density is used everywhere, then it is actually eliminated in
         ! this BC due to the scaling and hence unity is used instead.
         !----------------------------------------------------------------
         IF( UseDensity ) THEN
           Density = 1.0_dp
         ELSE
           Parent => Element % BoundaryInfo % Left
           stat = ASSOCIATED( Parent )
           IF (stat) stat = ALL(PresPerm(Parent % NodeIndexes(1:n)) > 0)
           IF ( .NOT. stat) THEN
             Parent => Element % BoundaryInfo % Right
             stat = ASSOCIATED( Parent )
             IF (stat) stat = ALL(PresPerm(Parent % NodeIndexes(1:n)) > 0)
             IF ( .NOT. stat )  CALL Fatal( 'HelmholtzSolve', &
                 'No parent element can be found for given boundary element' )
           END IF
           
           k = ListGetInteger( Model % Bodies(Parent % Bodyid) % Values, &
               'Material' )
           Material => Model % Materials(k) % Values
           Density(1:n) = ListGetReal( Material, 'Density', &
               n, Element % NodeIndexes(1:n) )
         END IF

         CALL LocalInterfaceMatrix(  STIFF, FORCE, AngularFrequency, Density, &
             Element, n, WallVelocity )
         
         CALL DefaultUpdateEquations( STIFF, FORCE )       
       END DO
     END IF


!------------------------------------------------------------------------------

     CALL DefaultFinishAssembly()
     CALL Info( 'HelmholtzSolve', 'Assembly done', Level=4 )

     CALL DefaultDirichletBCs()
!
!    Solve the system and we are done:
!    ---------------------------------
     at = CPUTime() - at
     st = CPUTime()
     Norm = DefaultSolve()

     st = CPUTIme()-st
     totat = totat + at
     totst = totst + st
     WRITE( Message, '(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
     CALL Info( 'HelmholtzSolve', Message, Level=4 )
     WRITE( Message, '(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
     CALL Info( 'HelmholtzSolve', Message, Level=4 )

!------------------------------------------------------------------------------

     IF( Solver % Variable % NonlinConverged == 1 ) EXIT

!------------------------------------------------------------------------------
  END DO ! of nonlinear iteration
!------------------------------------------------------------------------------

  CALL DefaultFinish()
  

CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, AngularFrequency, SoundSpeed, &
       ConvVelo, Damping, Load, Bubbles, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) ::  AngularFrequency, &
         SoundSpeed(:), Damping(:), Load(:,:), ConvVelo(:,:)
    COMPLEX(KIND=dp) :: STIFF(:,:), FORCE(:)
    LOGICAL :: Bubbles
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(2*nd),dBasisdx(2*nd,3)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,WaveNumber,M,D,L1,L2
    REAL(KIND=dp) :: DiffCoef(3,3), Velo(3), Rho
    COMPLEX(KIND=dp) :: A, B, ConvCoef
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim, NBasis, CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    REAL(KIND=dp) :: X,Y,Z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    Metric = 0.0_dp
    Metric(1,1) = 1.0_dp
    Metric(2,2) = 1.0_dp
    Metric(3,3) = 1.0_dp

    STIFF = 0.0_dp
    FORCE = 0.0_dp

    DiffCoef = 0.0_dp

!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    IF ( Bubbles ) THEN
       IntegStuff = GaussPoints( Element, Element % TYPE % GaussPoints2 )
       NBasis = 2*n
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
            Basis, dBasisdx, Bubbles=Bubbles )

       s = s * SqrtElementMetric
       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % X(1:n) * Basis(1:n) )
          Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
       END IF

!------------------------------------------------------------------------------
!      The source term and the coefficient of the time derivative and 
!      diffusion terms at the integration point
!------------------------------------------------------------------------------
       WaveNumber = AngularFrequency / SUM( SoundSpeed(1:n) * Basis(1:n) )

       D  =  WaveNumber * SUM( Damping(1:n) * Basis(1:n) )
       M  = -WaveNumber**2

       L1 = SUM( Load(1,1:n) * Basis(1:n) )
       L2 = SUM( Load(2,1:n) * Basis(1:n) )

       IF( UseDensity ) THEN
         Rho = SUM( Density(1:n) * Basis(1:n) ) 
       END IF

       DO i = 1,dim
         DiffCoef(i,i) = 1.0_dp
       END DO

       IF( UseConvection ) THEN 
         ConvCoef = 2.0_dp * SQRT((-1.0_dp,0.0_dp)) * WaveNumber

!        Scaled convection velocity
!        --------------------------
         Velo(1) = SUM( ConvVelo(1,1:n) * Basis(1:n) )
         Velo(2) = SUM( ConvVelo(2,1:n) * Basis(1:n) )
         Velo(3) = SUM( ConvVelo(3,1:n) * Basis(1:n) )
         Velo = Velo / SUM( SoundSpeed(1:n) * Basis(1:n) )

!        Diffusion and convection coefficients
!        -------------------------------------
      
         DO i = 1,dim
           DO j = 1,dim
             DiffCoef(i,j) = DiffCoef(i,j) - Velo(i)*Velo(j)
           END DO
         END DO
       END IF



!      Stiffness matrix and load vector
!      --------------------------------
       DO p=1,NBasis
          DO q=1,NBasis
             A = CMPLX( M, D,KIND=dp ) * Basis(q) * Basis(p)

             DO i=1,dim
               IF( UseConvection ) THEN
                 A = A + ConvCoef * Velo(i) * dBasisdx(q,i) * Basis(p)
               END IF
               DO j=1,dim
                 DO k = 1,dim
                   A = A + Metric(i,j) * DiffCoef(i,k) * dBasisdx(q,k) * dBasisdx(p,j)
                 END DO
               END DO
             END DO

             IF( UseDensity ) A = A / Rho
             STIFF(p,q) = STIFF(p,q) + s * A

          END DO
          
          B = Basis(p) * CMPLX( L1,L2,KIND=dp )
          IF( UseDensity ) B = B / Rho
          FORCE(p) = FORCE(p) + s * B
       END DO
    END DO
!------------------------------------------------------------------------------

    IF ( Bubbles ) THEN
       CALL CondensateP( n, n, STIFF, FORCE )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBoundary(  STIFF, FORCE, AngularFrequency, &
              Impedance, Load, Element, n, nd, ConvVelo )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: Impedance(:,:),Load(:,:)
    REAL(KIND=dp) :: AngularFrequency, ConvVelo(:,:)
    INTEGER :: n,nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Impedance1,Impedance2,L1,L2
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),X,Y,Z
    REAL(KIND=dp) :: Normal(3), Velo(3), NormVelo, TangVelo(3), Rho
    REAL(KIND=dp) :: Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    COMPLEX(KIND=dp) :: A, B, Admittance
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    STIFF = 0.0d0
    FORCE = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
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

       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % X(1:n) * Basis(1:n) )
          Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
       END IF

       Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

       Impedance1 = SUM( Impedance(1,1:n) * Basis(1:n) )
       Impedance2 = SUM( Impedance(2,1:n) * Basis(1:n) ) 
       IF ( ABS(Impedance1) < AEPS .AND. ABS(Impedance2) < AEPS) THEN
         Admittance = CMPLX(0.0d0,0.0d0,KIND=dp)
       ELSE         
         Admittance = CMPLX(0.0d0,1.0d0,KIND=dp) * AngularFrequency / CMPLX(Impedance1, Impedance2,KIND=dp)
       END IF

       IF( UseDensity ) THEN
         Rho = SUM( Density(1:n) * Basis(1:n) ) 
       END IF


       IF( UseConvection ) THEN
!        Scaled convection velocity
!        --------------------------
         Velo(1) = SUM( ConvVelo(1,1:n) * Basis(1:n) )
         Velo(2) = SUM( ConvVelo(2,1:n) * Basis(1:n) )
         Velo(3) = SUM( ConvVelo(3,1:n) * Basis(1:n) )
         Velo = Velo / SUM( SoundSpeed(1:n) * Basis(1:n) )
         NormVelo = SUM( Normal(1:dim) * Velo(1:dim) )
         TangVelo = Velo - Normal * NormVelo
      END IF

!------------------------------------------------------------------------------
       L1 = SUM( Load(1,1:n) * Basis(1:n) )
       L2 = SUM( Load(2,1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
       DO p=1,nd
          DO q=1,nd
             A = Admittance * Basis(q)*Basis(p)
             IF( UseConvection ) THEN
               A = A - NormVelo * Admittance * Basis(q)*Basis(p)
               A = A + SUM( TangVelo(1:dim) * dBasisdx(q,1:dim) ) * Basis(p)
             END IF
             IF( UseDensity ) A = A / Rho
             STIFF(p,q) = STIFF(p,q) + s * A
          END DO
          B = Basis(p) * CMPLX(L1,L2,KIND=dp)
          IF( UseDensity ) B = B / Rho
          FORCE(p) = FORCE(p) + s * B
       END DO
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalInterfaceMatrix(  STIFF, FORCE, AngularFrequency, Density, &
              Element, n, WallVelo )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: AngularFrequency, Density(:)
    COMPLEX(KIND=dp) :: WallVelo(:,:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Impedance1,Impedance2,L1,L2
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),X,Y,Z
    REAL(KIND=dp) :: Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: Normal(3), rho
    COMPLEX(KIND=dp) :: A, Admittance, NormVelo
    LOGICAL :: Stat
    INTEGER :: i,p,q,t,dim,CoordSys
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    CoordSys = CurrentCoordinateSystem()

    STIFF = 0.0d0
    FORCE = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
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
       
       IF ( CoordSys /= Cartesian ) THEN
          X = SUM( Nodes % X(1:n) * Basis(1:n) )
          Y = SUM( Nodes % Y(1:n) * Basis(1:n) )
          Z = SUM( Nodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
       END IF

       Normal = Normalvector(Element, Nodes, U, V, .TRUE.)

       NormVelo = 0.0_dp
       DO i=1,dim
         NormVelo = NormVelo + Normal(i) * SUM( WallVelo(i,1:n) * Basis(1:n) )
       END DO

       rho = SUM( Density(1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
       DO p=1,n
         FORCE(p) = FORCE(p) - s * Basis(p) * &
             ImUnit * rho * AngularFrequency * NormVelo
       END DO
!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalInterfaceMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE HelmholtzSolver
!------------------------------------------------------------------------------
