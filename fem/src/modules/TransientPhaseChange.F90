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
! *  Modified by: Peter Råback
! *  Modification date: 14.1.2010
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!>  Lagrangian & transient phase change solver for the liquid/solid interface.
!>  The equation is solved in two phases. First the velocity of the interface is
!>  computed and then the displacement. This enables independent smoothing of
!>  either velocity or displacement. There are also two different options for the 
!>  flux computations: internally or externally using the loads from the matrix
!>  residual. The latter one gives the more accurate approximation. The solver 
!>  includes some tailored pull control features that are originally intended
!>  for the modeling of Cz crystal growth but may also find other uses.
!>  Solve the free surface in the phase change problem using a transient algorithm. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE TransientPhaseChange( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: CurrentElement, Parent, Element
  TYPE(Variable_t), POINTER :: SurfSol, TempSol, HelpSol, LoadsSol
  TYPE(Nodes_t) :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff  
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Solver_t), POINTER :: PSolver 
  TYPE(ValueList_t), POINTER :: Params

  REAL(KIND=dp) :: Normal(3), u, v, w, UPull, PrevUpull, &
      Update, MaxUpdate, VelocityRelax, DispRelax, &
      surf, xx, yy, r, detJ, Temp, NonlinearTol, &
      d, HeatCond, s, CoordMin(3), CoordMax(3), RelativeChange, area, &
      Norm, maxds, maxds0, ds, MaxLoad, MinLoad, &
      pos0=0.0, prevpos0, Coeff, &
      MeanSurface, SpeedUp, LoadsRelax
  REAL(KIND=dp), POINTER :: Surface(:), PrevSurface(:), Temperature(:), ForceVector(:),  &
      x(:), y(:), z(:), Basis(:), dBasisdx(:,:), NodalTemp(:), &
      Conductivity(:), LatentHeat(:), Density(:), &
      Normals(:), Weights(:), SurfaceVelo(:), PrevSurfaceVelo(:), &
      CurrentLoads(:), PrevLoads(:)

  REAL(KIND=dp), ALLOCATABLE :: &          
      LocalStiffMatrix(:,:), LocalForceVector(:), LocalMassMatrix(:,:)  
  INTEGER :: i,j,k,t,n,nn,nd,pn,DIM,kl,kr,l, bc, Trip_node, axis_node, NonlinearIter, istat, &
       NElems,ElementCode,Next,Vertex,ii,imin,Node, iter, LiquidInd, Visited = -1, &
       SubroutineVisited = 0, NormalDirection, CoordMini(3), CoordMaxi(3), CoupledIter, &
       TimeStep, LoadsOrder
  INTEGER, POINTER :: NodeIndexes(:),TempPerm(:),SurfPerm(:),NormalsPerm(:)

  LOGICAL :: Stat, FirstTime = .TRUE., Debug, DoVelocityRelax, &
      PullControl, PullVelocitySet = .FALSE., IsoSurfAllocated, AllocationsDone = .FALSE., &
      UseLoads, AverageNormal, SurfaceVelocitySet = .FALSE., TriplePointFixed, &
      UseAverageLoads, UseFirstLoads
  LOGICAL, POINTER :: IsBoundaryNode(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, TemperatureName, str
  TYPE(ValueList_t), POINTER :: ValueList

  SAVE FirstTime, Trip_node, Axis_node,SubroutineVisited, prevpos0, &
      PrevSurfaceVelo, NormalDirection, ForceVector, PullControl, &
      Visited, Nodes, NodalTemp, Conductivity, LatentHeat, Density, &
      AllocationsDone, LocalStiffMatrix, LocalForceVector, LocalMassMatrix, &
      x, y, z, Basis, dBasisdx, norm, PullVelocitySet, &
      Normals, Weights, NormalsPerm, AverageNormal, SurfaceVelo, &
      SurfaceVelocitySet, CoordMax, CoordMin, CoordMaxi, CoordMini, UPull, &
      IsBoundaryNode, DoVelocityRelax, CurrentLoads, PrevLoads
  

  !------------------------------------------------------------------------------

  CALL Info('TransientPhaseChange','--------------------------------------------')
  CALL Info('TransientPhaseChange','Using transient algorithm for surface update')          
  CALL Info('TransientPhaseChange','--------------------------------------------')

  IF(.NOT. TransientSimulation ) THEN
    CALL Fatal('TransientPhaseChange','This only makes sense in a transient setting')
  END IF


  SubroutineVisited = SubroutineVisited + 1
  DIM = CoordinateSystemDimension()

  !------------------------------------------------------------------------------
  ! The variables needed for solution
  !------------------------------------------------------------------------------
  Params => GetSolverParams()

  PSolver => Solver
  SurfSol  => Solver % Variable
  Surface  => SurfSol % Values
  PrevSurface => SurfSol % PrevValues(:,1)
  SurfPerm => SurfSol % Perm
  IF(.NOT. ASSOCIATED (Surface) .OR. ALL(SurfPerm <= 0) ) THEN
    CALL Fatal('TransientPhaseChange','Surface field needed for Phase Change')
  END IF

  TemperatureName = ListGetString( Params, 'Phase Change Variable', Stat )
  IF(.NOT. Stat) TemperatureName = 'Temperature'

  TempSol => VariableGet( Solver % Mesh % Variables, TRIM(TemperatureName) )
  TempPerm    => TempSol % Perm
  Temperature => TempSol % Values
  IF(.NOT. ASSOCIATED (Temperature) .OR. ALL(TempPerm <= 0) ) THEN
    CALL Fatal('TransientPhaseChange','Temperature field needed for Phase Change')
  END IF

  NonlinearIter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations', stat )
  IF ( .NOT. stat ) NonlinearIter = 1    
  NonlinearTol  = ListGetConstReal( Params, &
      'Nonlinear System Convergence Tolerance', stat )

  PullControl = ListGetLogical( Params,'Pull Rate Control',stat)
  IF( PullControl ) THEN
    IF( dim == 3 ) THEN
      CALL Fatal('TransientPhaseChange','Pull rate control not implemented in 3D')
    END IF
    CALL Info('TransientPhaseChange','Using pull control for the phase change',Level=7)
  END IF

  TriplePointFixed = ListGetLogical( Params,'Triple Point Fixed',stat)
  IF( TriplePointFixed ) THEN
    IF( dim == 3 ) THEN
      CALL Fatal('TransientPhaseChange','Fixed triple point not implemented in 3D')
    END IF
    CALL Info('TransientPhaseChange','Fixing triple point position',Level=7)
  END IF

  HelpSol => VariableGet( Solver % Mesh % Variables, 'coupled iter')
  CoupledIter = HelpSol % Values(1)

  HelpSol => VariableGet( Solver % Mesh % Variables, 'timestep')
  TimeStep = HelpSol % Values(1)

!---------------------------------------------------------------------------------
! The first time the main axis of the free surface is determined
! and some permanent vectors related to the surface are allocated.
!---------------------------------------------------------------------------------

  IF(FirstTime) THEN
    CALL Info('TransientPhaseChange','Doing some first time initializations',Level=7)
    Trip_node = 0
    Axis_node = 0
    UPull = 0.0

    CoordMax = -HUGE(CoordMax)
    CoordMin = HUGE(CoordMin)
    
    DO k=1, Model % Mesh % NumberOfNodes
      IF( SurfPerm(k) <= 0) CYCLE
      
      DO j=1,DIM
        IF(j==1) xx = Model % Mesh % Nodes % x(k)
        IF(j==2) xx = Model % Mesh % Nodes % y(k)
        IF(j==3) xx = Model % Mesh % Nodes % z(k)
        IF(xx > CoordMax(j)) THEN
          CoordMax(j) = xx
          CoordMaxi(j) = k
        END IF
        IF(xx < CoordMin(j)) THEN
          CoordMin(j) = xx
          CoordMini(j) = k
        END IF
      END DO
    END DO
    
    ! Direction of minimum change is the normal direction if not given
    NormalDirection = ListGetInteger( Params,'Normal Direction',Stat)
    IF(.NOT. Stat) THEN
      j = 1
      DO i=1,DIM
        IF(CoordMax(i)-CoordMin(i) < CoordMax(j)-CoordMin(j)) THEN
          j = i
        END IF
      END DO
      NormalDirection = j
      CALL Info('TransientPhaseChange','Normal coordinate set to: '//TRIM(I2S(j)),Level=7)
    END IF
    
    ALLOCATE(IsBoundaryNode(SIZE(Surface)))
    IsBoundaryNode = .FALSE.

    ! In 2D the extremum points must be edge points
    IF( DIM == 2 ) THEN
      Trip_node = CoordMaxi(3 - NormalDirection)
      Axis_node = CoordMini(3 - NormalDirection) 
      IsBoundaryNode( SurfPerm(Trip_node) ) = .TRUE.
      IsBoundaryNode( SurfPerm(Axis_node) ) = .TRUE.
    END IF

    ! Otherwise select the points on edge using a flag
    DO t = 1, Solver % Mesh % NumberOfBoundaryElements      
      CurrentElement => GetBoundaryElement(t) 
      n  = GetElementNOFNodes()
      NodeIndexes => CurrentElement % NodeIndexes
      
      ValueList => GetBC()
      IF( .NOT. ListGetLogical(ValueList,'Phase Change Side',Stat)) CYCLE
      
      DO i=1,n
        j = SurfPerm(NodeIndexes(i))
        IF( j > 0 ) THEN
          IsBoundaryNode(j) = .TRUE.
        END IF
      END DO
    END DO

    n = COUNT( IsBoundaryNode )
    CALL Info('TransientPhaseChange','Number of boundary nodes: '//TRIM(I2S(n)),Level=7)

    n = Solver % Mesh % MaxElementNodes  
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        x(n), y(n), z(n), Basis(n), dBasisdx(n,3), NodalTemp(n), &
        Conductivity(n), LatentHeat(n), Density(n), &
        LocalStiffMatrix(n,n), LocalForceVector(n), LocalMassMatrix(n,n), &
        STAT=istat)
    IF ( istat /= 0 ) CALL Fatal( 'TransientPhaseChange', 'Memory allocation error 1.' )     
    

    ! Check whether normals are computed by an auxialiary solver
    !---------------------------------------------------------------------------------
    VariableName = ListGetString( Params, 'Normal Variable', Stat )
    IF(Stat) THEN
      HelpSol => VariableGet( Solver % Mesh % Variables, TRIM(VariableName), ThisOnly=.TRUE. )
    ELSE
      HelpSol => VariableGet( Solver % Mesh % Variables, 'Normals',ThisOnly=.TRUE. )
    END IF
    IF(ASSOCIATED(HelpSol)) THEN
      Normals => HelpSol % Values
      NormalsPerm => HelpSol % Perm
      AverageNormal = .TRUE.
      IF( NonlinearIter > 1 ) THEN
        CALL Warn('TransientPhaseChange','With external normal field there is no nonlinearity to iterate!')
      END IF
    ELSE
      AverageNormal = .FALSE.
    END IF

    ! The field is computed in two stages, 1st the velocity and then the displacement
    ! Ensure that also the velocity is allocated for
    !--------------------------------------------------------------------------------
    HelpSol => VariableGet( Solver % Mesh % Variables, &
         TRIM(ComponentName(Solver % Variable))//' Velo')
    IF(.NOT. ASSOCIATED(HelpSol)) THEN
      CALL Fatal('TransientPhaseChange','Surface Velo field should exist!')
    END IF
    SurfaceVelo => HelpSol % Values

    VelocityRelax = GetCReal(Params,'Velocity Relaxation Factor',DoVelocityRelax)
    IF(DoVelocityRelax) THEN
      ALLOCATE( PrevSurfaceVelo(SIZE(Surface)), STAT=istat)
      IF ( istat /= 0 ) CALL Fatal( 'TransientPhaseChange', 'Memory allocation error 2.' )     
      PrevSurfaceVelo = 0.0_dp
    END IF

    AllocationsDone = .TRUE.    
  END IF


  i =  ListGetInteger( Params,'Passive Steps',Stat)
  IF( i >= SubroutineVisited) GOTO 200

  ! The first pull velocity should always be set
  IF ( .NOT. PullVelocitySet ) THEN
    CurrentElement => GetActiveElement(1)
    k = GetMaterialId()
    Material => Model % Materials(k) % Values
    WRITE (str,'(A,I2)') 'Convection Velocity',NormalDirection
    UPull = ListGetConstReal( Material, str, Stat )
    PullVelocitySet = .TRUE.
  END IF

  ! Loads may be provided externally by using the 'Calculate Loads' flag in heat eq.
  !----------------------------------------------------------------------------------
  UseLoads = ListGetLogical( Params,'Use Nodal Loads',Stat) 

  IF(UseLoads) THEN
    LoadsSol => VariableGet( Solver % Mesh % Variables,TRIM(TemperatureName)//' Loads')
    IF(.NOT. ASSOCIATED(LoadsSol)) THEN
      CALL Fatal('TransientPhaseChange','Loads are requested to be used but missing!')
    ELSE
      CALL Info('TransientPhaseChange','Using nodal loads to determine phase change',Level=6)
    END IF

    IF(.NOT. ASSOCIATED( PrevLoads ) ) THEN
      CALL Info('TransientPhaseChange','Allocating for previous loads',Level=12)
      ALLOCATE( PrevLoads( SIZE( Surface ) ) )
      PrevLoads = 0.0_dp
    END IF
    IF(.NOT. ASSOCIATED( CurrentLoads ) ) THEN
      CALL Info('TransientPhaseChange','Allocating for current loads',Level=12)
      ALLOCATE( CurrentLoads( SIZE( Surface ) ) ) 
      CurrentLoads = 0.0_dp
    END IF

    UseFirstLoads = ListGetLogical( Params,'Use First Loads',Stat) 
    IF( .NOT. UseFirstLoads ) THEN
      IF( TimeStep > 1 .AND. CoupledIter == 1 ) THEN
        PrevLoads = CurrentLoads
      END IF
    END IF

    DO t=1,Solver % Mesh % NumberOfNodes
      i = SurfPerm(t)
      IF( i == 0) CYCLE
      j = LoadsSol % Perm(t)
      IF( j == 0 ) CYCLE

      CurrentLoads(i) = LoadsSol % Values(j)
    END DO

    IF( UseFirstLoads ) THEN
      IF( CoupledIter == 1 ) THEN
        PrevLoads = CurrentLoads 
      END IF
    END IF

    WRITE(Message,'(A,ES12.5)') 'Minimum nodal load at interface: ',MINVAL(CurrentLoads)
    CALL Info('TransientPhaseChange',Message ) 
    WRITE(Message,'(A,ES12.5)') 'Maximum nodal load at interface: ',MAXVAL(CurrentLoads)
    CALL Info('TransientPhaseChange',Message ) 
  END IF

  PrevUpull = Upull
  IF( DoVelocityRelax ) PrevSurfaceVelo = SurfaceVelo

  StiffMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS


  LoadsRelax = GetCReal(Params,'Loads Relaxation Factor',UseAverageLoads) 

  LoadsOrder = 0
  IF( UseAverageLoads .AND. ( .NOT. UseFirstLoads .OR. CoupledIter > 1 ) ) THEN
    LoadsOrder = MIN( 2, TimeStep - 1 ) 
  END IF


  DO iter = 1, NonlinearIter

    ! First solve the velocity field
    CALL DefaultInitialize()
    
    IF(UseLoads) THEN
      DO t=1,Solver % Mesh % NumberOfNodes
        i = SurfPerm(t)
        IF( i <= 0) CYCLE

        IF( LoadsOrder == 2 ) THEN
          ForceVector(i) = LoadsRelax * CurrentLoads(i) + &
              ( 1.0-LoadsRelax ) * PrevLoads(i)
        ELSE IF( LoadsOrder == 1 ) THEN
          ForceVector(i) = LoadsRelax * CurrentLoads(i) + &
              ( 1.0-LoadsRelax ) * PrevLoads(i)
        ELSE
          ForceVector(i) = CurrentLoads(i)
        END IF
      END DO
    END IF

    DO t = 1, Solver % NumberOfActiveElements         
      CurrentElement => GetActiveElement(t)

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      k = GetMaterialId()
      Material => Model % Materials(k) % Values
      NodeIndexes => CurrentElement % NodeIndexes
      
      LatentHeat(1:n) = ListGetReal( Material, 'Latent Heat', n, NodeIndexes )
      Density(1:n) = ListGetReal( Material, 'Density', n, NodeIndexes )
      Conductivity(1:n) = ListGetReal( Material,'Heat Conductivity', n, NodeIndexes )                
      
      CALL GetElementNodes( Nodes )
      CALL VelocityLocalMatrix( LocalStiffMatrix, LocalMassMatrix, LocalForceVector,&
          CurrentElement, n, nd, Nodes )
      
      CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForceVector )
    END DO

    CALL DefaultFinishBulkAssembly()
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()

    ! The Dirichlet conditions here are intended to enable periodic BCs only
    ! Actively setting the phase change to a specific value may be unphysical
    ! Feedback routines for pull velocity should be used instead.
    !-----------------------------------------------------------------------------
    CALL DefaultDirichletBCs()
    
    CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, SurfaceVelo, Norm, 1, Solver )        
    RelativeChange = Solver % Variable % NonlinChange
        
    IF ( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO
  
  IF( DoVelocityRelax ) THEN
    VelocityRelax = GetCReal(Params,'Velocity Relaxation Factor',Stat)
    SurfaceVelo = VelocityRelax * SurfaceVelo + (1-VelocityRelax) * PrevSurfaceVelo
  END IF

  WRITE( Message,'(A,ES12.5)') 'Minimum surface velocity: ',MINVAL(SurfaceVelo)
  CALL Info('TransientPhaseChange',Message,Level=12)
  WRITE( Message,'(A,ES12.5)') 'Maximum surface velocity: ',MAXVAL(SurfaceVelo)
  CALL Info('TransientPhaseChange',Message,Level=12)

  IF(PullControl .OR. TriplePointFixed) THEN      
    Upull = -SurfaceVelo(SurfPerm(Trip_node))
    IF(PullControl) THEN
      WRITE(Message,'(A,ES12.5)') 'Pull velocity: ', Upull
      CALL Info('TransientPhaseChange',Message) 
    END IF
  END IF
  

  ! Then solve the corresponding update in displacement field
  SpeedUp = ListGetConstReal( solver % Values,'Transient SpeedUp',Stat)
  IF(.NOT. Stat) SpeedUp = 1.0d0
  !-----------------------------------------------------------------------------------------
  ! If no smoothing is applied directly to the free surface then its update may be computed
  ! directly from the velocity field.
  !-----------------------------------------------------------------------------------------
  DispRelax = ListGetCReal(Params,'Surface Smoothing Factor',Stat) 

  ! Currently the PDE based version is disabled since without this one can used 
  ! higher order time discreatization
  IF( Solver % Order <= 1 .AND. ABS( DispRelax ) < EPSILON( DispRelax) ) THEN
    CALL Info('TransientPhaseChange','Updating surface by simple update',Level=7)
    Surface = PrevSurface + SpeedUp * dt * ( SurfaceVelo + Upull )
  ELSE
    CALL Info('TransientPhaseChange','Updating surfface by solving a PDE',Level=7)
    CALL DefaultInitialize()
    
    DO t = 1, Solver % NumberOfActiveElements       
      CurrentElement => GetActiveElement(t)

      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      NodeIndexes => CurrentElement % NodeIndexes
      
      CALL GetElementNodes( Nodes )
      CALL SurfaceLocalMatrix( LocalStiffMatrix, LocalMassMatrix, LocalForceVector,&
          CurrentElement, n, Nodes )
      
      CALL Default1stOrderTime( LocalMassMatrix,LocalStiffMatrix,&
          LocalForceVector)
      
      CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForceVector )
    END DO
    
    CALL DefaultFinishBulkAssembly()
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()

    ! The Dirichlet conditions here are intended to enable periodic BCs only
    CALL DefaultDirichletBCs()

    CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, Surface, Norm, 1, Solver )        
  END IF
  

  ! These integrate the pulled distance
  IF(PullControl .OR. TriplePointFixed) THEN
    IF(Visited /= Solver % DoneTime) THEN
      IF(Solver % DoneTime == 1) THEN
        prevpos0 = 0.0
      ELSE
        prevpos0 = pos0
      END IF
      Visited = Solver % DoneTime
    END IF
    
    pos0 = prevpos0 + UPull * dt
    
    ! This sets the maximum position of the crystal side that is still
    ! aligned with the pull direction. 
    IF(PullControl) CALL FindPullBoundary()
  END IF

  
200 CONTINUE
  IF(PullControl .OR. TriplePointFixed) THEN
    CALL ListAddConstReal(Model % Simulation,'res: Pull Position',pos0)       
    CALL ListAddConstReal( Model % Simulation,'res: Pull Velocity',UPull)
  END IF

  FirstTime = .FALSE.
  
!------------------------------------------------------------------------------
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE VelocityLocalMatrix( StiffMatrix, MassMatrix, ForceVector,&
      Element, n, nd, Nodes )
        
    ! external variables:
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)      
    INTEGER :: n, nd
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    
    ! internal variables:
    TYPE(Nodes_t) :: PNodes
    TYPE(Element_t), POINTER :: Parent      
    REAL(KIND=dp) :: Basis(3*n),dBasisdx(3*n,3), &
        X,Y,Z,U,V,W,S,detJ, TGrad(3,3),Flux,pu,pv,pw,pull, LocalHeat, LocalDens, &
        NodalTemp(3*n), xx(10),yy(10),zz(10),  NodalSurf(n), NodalNormal(3,n)
    REAL(KIND=dp) :: Velo, Ny, Normal(3), StabFactor, StabCoeff, DerSurf, xcoord
    LOGICAL :: Stat
    INTEGER :: i,j,k,l,t,p,q, pn, NBasis
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
 
    StabFactor = ListGetConstReal(Params,'Velocity Smoothing Factor',Stat )

    NodalSurf(1:n) = Surface( SurfPerm(NodeIndexes) )
   
    IF(AverageNormal) THEN
      IF( DIM == 2 ) THEN
        NodalNormal(1,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n))-1)
        NodalNormal(2,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n)))
        NodalNormal(3,1:n) = 0.0_dp
      ELSE
        NodalNormal(1,1:n) = Normals(3*NormalsPerm(NodeIndexes(1:n))-2)
        NodalNormal(2,1:n) = Normals(3*NormalsPerm(NodeIndexes(1:n))-1)
        NodalNormal(3,1:n) = Normals(3*NormalsPerm(NodeIndexes(1:n)))
      END IF
    END IF

    ALLOCATE( PNodes % x(10), PNodes % y(10), PNodes % z(10) )

!      Numerical integration:
!      ----------------------

    NBasis = nd
    IntegStuff = GaussPoints( Element )

    DO t = 1,IntegStuff % n
      
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)

      
!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,detJ,Basis,dBasisdx )
      s = s * detJ
      xcoord = SUM( Nodes % x(1:n) * Basis(1:n) )


! This is not really a physical equation and hence weighing with the radius is not necessary
!      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
!        s = s * xcoord
!      END IF
      
      IF(AverageNormal) THEN
        Normal(1) = SUM( Basis(1:n) * NodalNormal(1,1:n))
        Normal(2) = SUM( Basis(1:n) * NodalNormal(2,1:n))
        Normal(3) = SUM( Basis(1:n) * NodalNormal(3,1:n))
      ELSE
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )         
      END IF

      LocalHeat = SUM(Basis(1:n) * LatentHeat(1:n))
      LocalDens = SUM(Basis(1:n) * Density(1:n))
      StabCoeff = StabFactor * LocalHeat * LocalDens
      
      IF(UseLoads) THEN
        ! do nothing, loads already inserted

      ELSE 
        ! Compute the flux from normal derivaties
        TGrad = 0.0d0          
        l = 0
        DO i=1,2
          
          IF( i == 1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          
          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
          IF (ListGetLogical(Model % Materials(k) % Values,'Solid',stat)) THEN
            IF(l == 2) CALL Fatal('TransientPhaseChange','Both materials cannot be solid!')
            l = 2
          ELSE 
            IF(l == 1) CALL Fatal('TransientPhaseChange','Both materials cannot be liquid!')
            l = 1
          END IF
          
          pn = Parent % TYPE % NumberOfNodes
          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
          
          Conductivity(1:pn) = ListGetReal( Model % Materials(k) % Values, &
              'Heat Conductivity', pn, Parent % NodeIndexes )                     
          NodalTemp(1:pn) = Temperature( TempPerm(Parent % NodeIndexes) )
          
          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
          
          ! Calculate the basis functions for the parent element:
          !-----------------------------------------------------
          DO j = 1,n
            DO k = 1,pn
              IF( NodeIndexes(j) == Parent % NodeIndexes(k) ) THEN
                xx(j) = Parent % TYPE % NodeU(k)
                yy(j) = Parent % TYPE % NodeV(k)
                zz(j) = Parent % TYPE % NodeW(k)
                EXIT
              END IF
            END DO
          END DO
          
          pu = SUM( Basis(1:n) * xx(1:n) )
          pv = SUM( Basis(1:n) * yy(1:n) )
          pw = SUM( Basis(1:n) * zz(1:n) )
          
          PNodes % x(1:pn) = Solver % Mesh % Nodes % x(Parent % NodeIndexes)
          PNodes % y(1:pn) = Solver % Mesh % Nodes % y(Parent % NodeIndexes)
          PNodes % z(1:pn) = Solver % Mesh % Nodes % z(Parent % NodeIndexes)
          
          stat = ElementInfo( Parent, PNodes, pu, pv, pw, detJ, Basis, dBasisdx )
          
          DO j=1,DIM
            TGrad(l,j) = SUM( Conductivity(1:pn) * Basis(1:pn) ) * &
                SUM( dBasisdx(1:pn,j) * NodalTemp(1:pn) )
          END DO
        END DO
        
        Flux = SUM( (TGrad(1,:) - TGrad(2,:)) *  Normal)       
        stat = ElementInfo( Element,Nodes,U, V, W, detJ, Basis, dBasisdx )
      END IF

      ! Assembly the matrix
      DO p=1,nBasis        
        
        DO q=1,nBasis
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * Basis(p) * Basis(q) * &
              Normal(NormalDirection) * LocalDens * LocalHeat            
          
          ! Jump over the special nodes (when defined)
          IF( IsBoundaryNode( SurfPerm(NodeIndexes(p))) ) CYCLE

          DO i=1,DIM
            IF( i == NormalDirection ) CYCLE
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,i) * dBasisdx(p,i)
          END DO
        END DO
                
        ! transient part of heat flux
        IF(.NOT. UseLoads) THEN
          ForceVector(p) = ForceVector(p) - s * Basis(p) * Flux 
        END IF
        
      END DO

    END DO
    
    DEALLOCATE( PNodes % x, PNodes % y, PNodes % z )

!------------------------------------------------------------------------------
  END SUBROUTINE VelocityLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SurfaceLocalMatrix( StiffMatrix, MassMatrix, ForceVector,&
      Element, nCoord, Nodes )
        
    ! external variables:
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)      
    INTEGER :: nCoord
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    
    ! internal variables:
    REAL(KIND=dp) :: Basis(3*nCoord),dBasisdx(3*nCoord,3), &
        X,Y,Z,U,V,W,S,detJ,pull, NodalVelo(nCoord), xcoord, NodalNormal(3,nCoord)
    REAL(KIND=dp) :: Velo, StabFactor, StabCoeff
    LOGICAL :: Stat
    INTEGER :: i,j,k,l,t,p,q, n, pn
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Params,'Surface Smoothing Factor')
    StabCoeff = StabFactor / dt

    NodalVelo(1:n) = SurfaceVelo( SurfPerm(NodeIndexes) )


!      Numerical integration:
!      ----------------------

    IntegStuff = GaussPoints( Element )
    
    DO t = 1,IntegStuff % n
      
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
      s = IntegStuff % s(t)
      
!        Basis function values & derivatives at the integration point:
!        -------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,U,V,W,detJ,Basis,dBasisdx)
      s = s * detJ
      xcoord = SUM( Nodes % x(1:nCoord) * Basis(1:nCoord) )
            
      Velo = SpeedUp * ( SUM( Basis(1:n) * NodalVelo(1:n)) + Upull )

      DO p=1,n        
        DO q=1,n
          MassMatrix(p,q) = MassMatrix(p,q) + s * Basis(p) * Basis(q)  

          ! Jump over the special nodes (when defined)
          IF( IsBoundaryNode( SurfPerm(NodeIndexes(p))) ) CYCLE

          DO i=1,DIM
            IF( NormalDirection == i ) CYCLE
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,i) * dBasisdx(p,i)            
          END DO
        END DO        
        ForceVector(p) = ForceVector(p) + s * Basis(p) * Velo         
      END DO

    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE SurfaceLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    SUBROUTINE FindPullBoundary()

      REAL (KIND=dp) :: Ybot, Ymax, Ytop, Xtrip, x, y
      INTEGER :: t,k,i,n
      LOGICAL :: PullBoundary

      IF(NormalDirection /= 2) CALL Fatal('TransientPhaseChange',&
          'FindPullBoundary implemented only for lateral boundaries!')

      PullBoundary = .FALSE.
      Ybot = Solver % Mesh % Nodes % y(Trip_node)
      Ymax = MAXVAL(Solver % Mesh % Nodes % y)
      Ytop = Ymax
      Xtrip = Solver % Mesh % Nodes % x(Trip_node)
      
      DO t = Solver % Mesh % NumberOfBulkElements + 1, &
          Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements        
        
        CurrentElement => Solver % Mesh % Elements(t)
        Model % CurrentElement => CurrentElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        DO k=1, Model % NumberOfBCs
          IF ( Model % BCs(k) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
          IF( ListGetLogical(Model % BCs(k) % Values,'Pull Boundary',stat ) ) THEN
            PullBoundary = .TRUE.
            DO i = 1,n
              x = Solver % Mesh % Nodes % x(NodeIndexes(i))
              y = Solver % Mesh % Nodes % y(NodeIndexes(i))
              IF(y > Ybot .AND. ABS(x-Xtrip) > 1.0e-6 * (Ymax-Ybot)) Ytop = y
            END DO
          END IF
        END DO
      END DO
      
      IF (PullBoundary) THEN
        CALL ListAddConstReal(Model % Simulation,'res: Triple point position',Ybot)
        CALL ListAddConstReal(Model % Simulation,'res: Full pull position',Ytop)
      END IF

    END SUBROUTINE FindPullBoundary

!------------------------------------------------------------------------------
  END SUBROUTINE TransientPhaseChange
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: TransientPhaseChange.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE TransientPhaseChange_Init( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL ::  TransientSimulation
    REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName

    Params => GetSolverParams()

    VariableName = GetString(Params,'Variable')
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          '-nooutput '//TRIM(ComponentName(VariableName))//' Velo' )

    IF(.NOT. ListCheckPresent( Params,'Time Derivative Order') ) &
        CALL ListAddInteger( Params,'Time Derivative Order',1)


!------------------------------------------------------------------------------
END SUBROUTINE TransientPhaseChange_Init
!------------------------------------------------------------------------------
