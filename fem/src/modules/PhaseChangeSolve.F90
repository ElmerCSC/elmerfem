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
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Modified by: Jussi Heikonen, Ville Savolainen, Peter Raback
! *  Modification date: 2.9.2006
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization for the primary solver: PhaseChangeSolve
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PhaseChangeSolve_Init( Model,Solver,dt,TransientSimulation)
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
          '-nooutput '//TRIM(ComponentName(VariableName))//'Move' )

    IF (ListGetLogical(Params,'Use Average Velocity',Found)) &
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          '-nooutput '//TRIM(ComponentName(VariableName))//'MoveAve' )
!------------------------------------------------------------------------------
END SUBROUTINE PhaseChangeSolve_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Solve the free surface in the phase change problem using Lagrangian techniques. 
!> \deprecated This had been replaced by separate versions for transient and steady state phase change.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE PhaseChangeSolve( Model,Solver,dt,TransientSimulation )
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
  TYPE(Variable_t), POINTER :: SurfSol, TempSol, HelpSol, HelpSol2
  TYPE(Nodes_t) :: Nodes, PNodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff  
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Solver_t), POINTER :: PSolver 

  REAL(KIND=dp) :: Normal(3), u, v, w, UPull(3), PrevUpull(3), &
      Density, Update, MaxUpdate, MaxTempDiff, Relax, LocalRelax, AverageRelax, &
      AverageRelax2, surf, xx, yy, r, detJ, Temp, MeltPoint, NonlinearTol, &
      Temp1,Temp2,Temppi,dxmin,dymin, xmin, xmax, d, HeatCond, s, tave, prevtave, cvol, clim, ccum=1, &
      tabs, prevtabs, volabs, prevvolabs, CoordMin(3), CoordMax(3), RelativeChange, &
      dTdz, ttemp, NewtonAfterTol, area, volume, prevvolume, Eps, &
      Norm, maxds, maxds0, ds, FluxCorrect = 1.0, dpos, &
      pos0=0.0, prevpos0, ThermalInertiaFactor, Coeff, OrigNorm, &
      SteadyDt, AveragingDt, Trip_Temp, MeanSurface, PullRelax, &
      SpeedUp, PrevNorm
  REAL(KIND=dp), POINTER :: Surface(:), Temperature(:), ForceVector(:),  &
      x(:), y(:), z(:), Basis(:), dBasisdx(:,:), NodalTemp(:), &
      Conductivity(:), LatentHeat(:), TempDiff(:), Taverage(:), Taverage2(:), &
      Normals(:), Weights(:), SurfaceMove(:), SurfaceMoveAve(:)
  REAL (KIND=dp), ALLOCATABLE :: PrevTemp(:), IsoSurf(:,:)

  REAL(KIND=dp), ALLOCATABLE :: &          
      LocalStiffMatrix(:,:), LocalForceVector(:), LocalMassMatrix(:,:)
  
  INTEGER :: i,j,k,t,n,nn,pn,DIM,kl,kr,l, bc, Trip_node, NoBNodes, NonlinearIter, istat, &
       NElems,ElementCode,Next,Vertex,ii,imin,NewtonAfterIter,Node, iter, LiquidInd, Visited = -1, &
       SubroutineVisited = 0, NormalDir, TangentDirection, CoordMini(3), CoordMaxi(3), &
       Axis_node
  INTEGER, POINTER :: NodeIndexes(:),TempPerm(:),SurfPerm(:),NormalsPerm(:)

  LOGICAL :: Stat, FirstTime = .TRUE., Newton = .FALSE., UseTaverage, Debug, &
      PullControl, PullVelocitySet = .FALSE., IsoSurfAllocated, AllocationsDone = .FALSE., &
      TransientAlgo = .FALSE., SteadyAlgo = .FALSE., UseLoads, &
      AverageNormal, SurfaceVelocitySet = .FALSE., TriplePointFixed
  LOGICAL, POINTER :: NodeDone(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, str

  SAVE FirstTime, Trip_node, NoBNodes, SubroutineVisited, prevpos0, &
      PrevTemp, Newton, ccum, NormalDir, TangentDirection, MeltPoint, &
      ForceVector, FluxCorrect, NodeDone, Eps, PullControl, &
      Visited, Nodes, PNodes, NodalTemp, Conductivity, LatentHeat, &
      TempDiff, AllocationsDone, LocalStiffMatrix, LocalForceVector, LocalMassMatrix, &
      x, y, z, Basis, dBasisdx, norm, Axis_node, PullVelocitySet, &
      Normals, Weights, NormalsPerm, AverageNormal, SurfaceMove, SurfaceMoveAve, &
      SurfaceVelocitySet, CoordMax, CoordMin, CoordMaxi, CoordMini, UPull
  
  !------------------------------------------------------------------------------
  ! Decide which kind of algorithm to use for the current timestep size
  !------------------------------------------------------------------------------

  TransientAlgo = TransientSimulation
  IF(TransientAlgo) THEN
    SteadyDt = ListGetConstReal(Solver % Values,'Steady Transition Timestep',Stat)
    IF(Stat) TransientAlgo = (dt < SteadyDt)
  END IF
  SteadyAlgo = .NOT. TransientAlgo    

  CALL Info('PhaseChangeSolve',                  '--------------------------------------------')
  IF(SteadyAlgo) CALL Info('PhaseChangeSolve',   'Using steady algorithm to find the isotherm')      
  IF(TransientAlgo) CALL Info('PhaseChangeSolve','Using transient algorithm for surface update')          
  CALL Info('PhaseChangeSolve',                  '--------------------------------------------')

  SubroutineVisited = SubroutineVisited + 1

  DIM = CoordinateSystemDimension()
  IF(DIM /= 2) THEN
    CALL Fatal('PhaseChangeSolve','Implemented only in 2D')
  END IF

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------

  PSolver => Solver
  SurfSol  => Solver % Variable
  Surface  => SurfSol % Values
  SurfPerm => SurfSol % Perm
  IF(.NOT. ASSOCIATED (Surface) .OR. ALL(SurfPerm <= 0) ) THEN
    CALL Fatal('PhaseChangeSolve','Surface field needed for Phase Change')
  END IF

  VariableName = ListGetString( Solver % Values, 'Phase Change Variable', Stat )
  IF(Stat) THEN
    TempSol => VariableGet( Solver % Mesh % Variables, TRIM(VariableName) )
  ELSE
    TempSol => VariableGet( Solver % Mesh % Variables, 'Temperature' )
  END IF
  TempPerm    => TempSol % Perm
  Temperature => TempSol % Values
  IF(.NOT. ASSOCIATED (Temperature) .OR. ALL(TempPerm <= 0) ) THEN
    CALL Fatal('PhaseChangeSolve','Temperature field needed for Phase Change')
  END IF

  Relax = GetCReal( Solver % Values,  & 
      'Nonlinear System Relaxation Factor', stat )
  IF ( .NOT. stat ) Relax = 1.0d0
  NonlinearIter = ListGetInteger( Solver % Values, &
      'Nonlinear System Max Iterations', stat )
  IF ( .NOT. stat ) NonlinearIter = 1    
  NonlinearTol  = ListGetConstReal( Solver % Values, &
      'Nonlinear System Convergence Tolerance', stat )

  PullControl = ListGetLogical( Solver % Values,'Pull Rate Control',stat)
  TriplePointFixed = ListGetLogical( Solver % Values,'Triple Point Fixed',stat)

!---------------------------------------------------------------------------------
! The first time the main axis of the free surface is determined
! and some permanent vectors related to the surface are allocated.
!---------------------------------------------------------------------------------
  
  IF(FirstTime) THEN
    UPull = 0.0
    NoBNodes = 0    
    CoordMax = -HUGE(CoordMax)
    CoordMin = HUGE(CoordMin)
    
    DO k=1, Model % Mesh % NumberOfNodes
      IF( SurfPerm(k) <= 0) CYCLE
      NoBnodes = NoBnodes + 1

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
    
    ! Direction of minimum change
    j = 1
    DO i=1,DIM
      IF(CoordMax(i)-CoordMin(i) < CoordMax(j)-CoordMin(j)) THEN
        j = i
      END IF
    END DO
    NormalDir = j
    
    ! Direction of maximum change
    j = 1
    DO i=1,DIM
      IF(CoordMax(i)-CoordMin(i) > CoordMax(j)-CoordMin(j)) THEN
        j = i
      END IF
    END DO
    TangentDirection = j

    ! Find triple point:
    !-------------------
    Trip_node = CoordMaxi(TangentDirection)
    Axis_node = CoordMini(TangentDirection)

    n = Solver % Mesh % MaxElementNodes  

    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        PNodes % x(n), PNodes % y(n), PNodes % z(n), &
        x(n), y(n), z(n), Basis(n), dBasisdx(n,3), NodalTemp(n), &
        Conductivity(n), LatentHeat(n), TempDiff(n), &
        LocalStiffMatrix(n,n), LocalForceVector(n), LocalMassMatrix(n,n), &
        PrevTemp(NobNodes), &
        STAT=istat)
    IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error 1.' )     
    
    Nodes % x = 0.0d0
    Nodes % y = 0.0d0
    Nodes % z = 0.0d0
    PrevTemp = 0.0d0

    IF( SteadyAlgo ) THEN
      ALLOCATE( NodeDone( NobNodes ), STAT=istat)
      IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error 3.' )          
    END IF  

    VariableName = ListGetString( Solver % Values, 'Normal Variable', Stat )
    IF(Stat) THEN
      HelpSol => VariableGet( Solver % Mesh % Variables, TRIM(VariableName), ThisOnly=.TRUE. )
    ELSE
      HelpSol => VariableGet( Solver % Mesh % Variables, 'Normals',ThisOnly=.TRUE. )
    END IF
    IF(ASSOCIATED(HelpSol)) THEN
      Normals => HelpSol % Values
      NormalsPerm => HelpSol % Perm
      AverageNormal = .TRUE.
    ELSE
      AverageNormal = .FALSE.
    END IF

    HelpSol => VariableGet( Solver % Mesh % Variables, &
         TRIM(ComponentName(Solver % Variable))//'Move', stat)
    IF(ASSOCIATED(HelpSol)) THEN
      SurfaceMove => HelpSol % Values
    ELSE
      ALLOCATE( SurfaceMove(SIZE(Surface)), STAT=istat)
      IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error.' )           
      SurfaceMove = 0.0d0
      CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
          PSolver,TRIM(ComponentName(Solver % Variable))//'Move',1, &
              SurfaceMove,SurfPerm,Output=.FALSE.)
    END IF

    IF(ListGetLogical(Solver % Values,'Use Average Velocity',Stat)) THEN
      HelpSol => VariableGet( Solver % Mesh % Variables, &
             TRIM(ComponentName(Solver % Variable))//'MoveAve', stat)
      IF(ASSOCIATED(HelpSol)) THEN
        SurfaceMoveAve => HelpSol % Values
      ELSE
        ALLOCATE( SurfaceMoveAve(SIZE(Surface)), STAT=istat)
        IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error.' )           
        SurfaceMoveAve = 0.0d0
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
            PSolver,TRIM(ComponentName(Solver % Variable))//'MoveAve', &
               1,SurfaceMoveAve,SurfPerm,Output=.FALSE.)
      END IF
    END IF

    AllocationsDone = .TRUE.    
  END IF

  ! Find triple point temperature
  !-------------------------------

  Trip_Temp =  Temperature( TempPerm(Trip_node) )    


  i =  ListGetInteger( Solver % Values,'Passive Steps',Stat)
  IF( i >= SubroutineVisited) GOTO 200

  ! Find the constant material parameters
  !--------------------------------------

  CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements (1))
  n = CurrentElement % TYPE % NumberOfNodes
  NodeIndexes => CurrentElement % NodeIndexes
  k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values,'Material')
  Material => Model % Materials(k) % Values   
  MeltPoint = ListGetConstReal( Material,'Melting Point' )
  Density = ListGetConstReal( Material, 'Density' )

  ! The first velocity should always be set
  IF ( .NOT. PullVelocitySet ) THEN
    Upull(1) = ListGetConstReal( Material, 'Convection Velocity 1', Stat )
    Upull(2) = ListGetConstReal( Material, 'Convection Velocity 2', Stat )
    Upull(3) = 0.0
    PullVelocitySet = .TRUE.
  END IF

!--------------------------------------------------------------------
! The transient algorithm 
! In the transient algorithm a heat flux over the interface is computed and 
! it is assumed to be used solely in the melting of the solid into liquid. 
! This melting speed gives an estimate for the melting speed that may be 
! improved by iteration.
!--------------------------------------------------------------------

  IF ( TransientAlgo ) THEN
    StiffMatrix => Solver % Matrix
    ForceVector => Solver % Matrix % RHS

    UseLoads = ListGetLogical( Solver % Values,'Use Heat Load',Stat) 
    IF(UseLoads) THEN
      HelpSol => VariableGet( Solver % Mesh % Variables, 'Nodal Heat Load' )
    END IF
    
    PrevUpull = Upull
    PrevTemp = SurfaceMove
    
    ! nonlinear iteration could only be associated if normal is computed from the solution
    DO iter = 1, NonlinearIter
      
      ! First solve the velocity field
      CALL InitializeToZero( StiffMatrix, ForceVector )
      
      IF(UseLoads) THEN
        DO t=1,Solver % Mesh % NumberOfNodes
          i = SurfPerm(t)
          IF( i <= 0) CYCLE
          ForceVector(i) = HelpSol % Values(HelpSol % Perm(t))
        END DO
      END IF
      
      DO t = 1, Solver % NumberOfActiveElements         
        CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements (t))
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values,'Material')
        Material => Model % Materials(k) % Values
        
        LatentHeat(1:n) = ListGetReal( Material, 'Latent Heat', n, NodeIndexes )
        Conductivity(1:n) = ListGetReal( Material,'Heat Conductivity', n, NodeIndexes )                
        TempDiff(1:n) = Temperature( TempPerm(NodeIndexes) ) - MeltPoint
        
        CALL VelocityLocalMatrix( LocalStiffMatrix, LocalMassMatrix, LocalForceVector,&
            CurrentElement, n, Nodes )
        
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
            ForceVector, LocalForceVector, n, 1, SurfPerm(NodeIndexes) )
      END DO
      
      CALL FinishAssembly( Solver, ForceVector )

      ! No Dirihtlet conditions here since       
      ! One should not really try to force the phase change at some point, 
      ! rather use feedback to tune the pull velocity
      
      CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, SurfaceMove, Norm, 1, Solver )        
      RelativeChange = Solver % Variable % NonlinChange

      WRITE( Message, * ) 'Result Norm     : ',Norm
      CALL Info( 'PhaseChangeSolve', Message, Level=4 )

      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'PhaseChangeSolve', Message, Level=4 )

      IF ( RelativeChange < NonLinearTol ) EXIT
    END DO
            
    IF(ListGetLogical(Solver % Values,'Use Average Velocity',Stat)) THEN
      IF(SurfaceVelocitySet) THEN            
        LocalRelax = GetCReal(Solver % Values,'Velocity Averaging Factor')
        SurfaceMoveAve = LocalRelax * SurfaceMove + (1-LocalRelax) * SurfaceMoveAve               
      ELSE
        SurfaceMoveAve = SurfaceMove
        SurfaceVelocitySet = .TRUE.
      END IF
      LocalRelax = GetCReal(Solver % Values,'Velocity Relaxation Factor',Stat)
      IF(.NOT. Stat) LocalRelax = 1.0d0
      SurfaceMove = LocalRelax * SurfaceMove + (1-LocalRelax) * SurfaceMoveAve
    ELSE
      LocalRelax = GetCReal(Solver % Values,'Velocity Relaxation Factor',Stat)
      IF(.NOT. Stat) LocalRelax = 1.0d0
      SurfaceMove = LocalRelax * SurfaceMove + (1-LocalRelax) * PrevTemp
    END IF

    IF(PullControl .OR. TriplePointFixed) THEN      
      Upull(NormalDir) = -SurfaceMove(SurfPerm(Trip_node))
      IF(PullControl) THEN
        WRITE(Message,*) 'Pull velocity: ', Upull(NormalDir)
        CALL Info('PhaseChangeSolve',Message) 
      END IF
    END IF

      
    ! Then solve the corresponding update in displacement field
    SpeedUp = ListGetConstReal( solver % Values,'Transient SpeedUp',Stat)
    IF(.NOT. Stat) SpeedUp = 1.0d0
    !-----------------------------------------------------------------------------------------
    ! If no smoothing is applied directly to the free surface then its update may be computed
    ! directly from the velocity field.
    !-----------------------------------------------------------------------------------------
    IF( ListGetConstReal(Solver % Values,'Surface Smoothing Factor') < 1.0d-20) THEN
      Surface = Surface + SpeedUp * dt * ( SurfaceMove + Upull(NormalDir))
    ELSE
      CALL InitializeToZero( StiffMatrix, ForceVector )
      
      DO t = 1, Solver % NumberOfActiveElements       
        CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements (t))
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        CALL SurfaceLocalMatrix( LocalStiffMatrix, LocalMassMatrix, LocalForceVector,&
            CurrentElement, n, Nodes )
        
        CALL Add1stOrderTime( LocalMassMatrix, LocalStiffMatrix, &
            LocalForceVector, dt, n, 1, SurfPerm(NodeIndexes), Solver )
        
        CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
            ForceVector, LocalForceVector, n, 1, SurfPerm(NodeIndexes) )
      END DO
      
      CALL FinishAssembly( Solver, ForceVector )    
      CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, Surface, Norm, 1, Solver )        
    END IF
 
    IF(PullControl .OR. TriplePointFixed) THEN
      IF(Visited /= Solver % DoneTime) THEN
        IF(Solver % DoneTime == 1) THEN
          prevpos0 = 0.0
        ELSE
          prevpos0 = pos0
        END IF
        Visited = Solver % DoneTime
      END IF
      
      pos0 = prevpos0 + UPull(NormalDir) * dt
      
      IF(PullControl) THEN
        ! This sets the maximum position of the crystal side that is still
        ! aligned with the pull direction. 
        CALL FindPullBoundary()
      END IF

    END IF
  END IF

!--------------------------------------------------------------------
! Transient algorithm end
!--------------------------------------------------------------------


!----------------------------------------------------------------------------
! The Steady State the simulation is based on a geometric determination of the 
! isotherm. The solution may be accelerated using local or global Newton 
! type of iteration. It is activated only after the solution is quite accurate.
!-----------------------------------------------------------------------------

  IF( SteadyAlgo ) THEN

    UseTaverage = ListGetLogical( Solver % Values,'Average Temperature', stat )
    AveragingDt = ListGetConstReal(Solver % Values,'Steady Averaging Timestep',Stat)

    IF(UseTaverage) THEN
      HelpSol => VariableGet( Solver % Mesh % Variables, 'Taverage' )
      IF(.NOT. ASSOCIATED (HelpSol)) THEN
        ALLOCATE( Taverage( Model % NumberOfNodes ), STAT=istat )
        IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error.' )     
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
            PSolver, 'Taverage', 1, Taverage, TempPerm)
        HelpSol => VariableGet( Solver % Mesh % Variables, 'Taverage' )
      END IF

      NULLIFY(HelpSol2)
      AverageRelax2 = GetCReal(Solver % Values,'Temperature Averaging Factor',stat)
      IF(stat) THEN
        HelpSol2 => VariableGet( Solver % Mesh % Variables, 'Taverage Slow' )
        IF(.NOT. ASSOCIATED (HelpSol2)) THEN
          ALLOCATE( Taverage2( Model % NumberOfNodes ), STAT=istat )
          IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error.' )     
          CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
              PSolver, 'Taverage Slow', 1, Taverage2, TempPerm)
          HelpSol2 => VariableGet( Solver % Mesh % Variables, 'Taverage Slow' )
        END IF
      END IF

      IF(dt >= AveragingDt) THEN
        HelpSol % Values = TempSol % Values
        IF(ASSOCIATED(HelpSol2)) HelpSol2 % Values = TempSol % Values        
      ELSE      
        AverageRelax = GetCReal( Solver % Values,'Temperature Relaxation Factor')
        IF(.NOT. ASSOCIATED(HelpSol2)) THEN
          HelpSol % Values = (1.0d0-AverageRelax) * HelpSol % Values + AverageRelax * TempSol % Values
        ELSE
          HelpSol2 % Values = (1.0d0-AverageRelax2) * HelpSol2 % Values + AverageRelax2 * TempSol % Values
          HelpSol % Values = (1.0d0-AverageRelax) * HelpSol2 % Values + AverageRelax * TempSol % Values
        END IF
        PRINT *,'temperature set to average temperature'
        Temperature => HelpSol % Values
      END IF
    END IF

    NewtonAfterIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations', stat )    
    IF ( stat .AND. SubroutineVisited > NewtonAfterIter ) Newton = .TRUE.

    NewtonAfterTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance', stat )

    IF(Newton) THEN
      CALL Info( 'PhaseChangeSolve','Steady state newton formulation', Level=4 )
    ELSE
      CALL Info( 'PhaseChangeSolve','Steady state isoterm formulation', Level=4 )
    END IF

! Find the melting point:
!-----------------------

    IF( TriplePointFixed ) THEN
      MeltPoint = Temperature( TempPerm(Trip_node) )        
      WRITE(Message,*) 'Melting point set to triple point temperature',MeltPoint,Trip_node
      CALL Info('PhaseChangeSolve',Message)        
    END IF


! Create the isotherm for T=T_0:
!------------------------------

    IF(.NOT. Newton) THEN
      
      IsoSurfAllocated = .FALSE.
      xmin = HUGE(xmin)
      xmax = -HUGE(xmax)

100   NElems = 0
      
      DO t=1,Solver % Mesh % NumberOfBulkElements 
        
        CurrentElement => Solver % Mesh % Elements(t)
        
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes 
        ElementCode = CurrentElement % TYPE % ElementCode
        
        k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values,'Material')
        IF (.NOT. (ListGetLogical(Model % Materials(k) % Values, 'Solid', stat) .OR. &
            ListGetLogical(Model % Materials(k) % Values, 'Liquid', stat) )) CYCLE

        IF(ElementCode < 300) CYCLE

        IF ( ANY( TempPerm( NodeIndexes(1:n) ) <= 0 ) )  CYCLE
        TempDiff(1:n) = Temperature(TempPerm(NodeIndexes(1:n))) - MeltPoint

        IF( ALL ( TempDiff(1:n) < 0.0 ) ) CYCLE
        IF( ALL ( TempDiff(1:n) > 0.0 ) ) CYCLE       

        Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)

        Vertex = ElementCode/100
        n=0
        DO nn=1,Vertex
          next  = MODULO(nn,Vertex) + 1
          
          temp1 = TempDiff(nn)
          temp2 = TempDiff(next)
                    
          IF ( ( (temp1 < 0.0) .AND. (0.0 <= temp2) ) .OR. &
              ( (temp2 <= 0.0) .AND. (0.0 < temp1) ) ) THEN
            
            n = n + 1
            
            IF ( n <= 2 ) THEN
              NElems = NElems + 1              
              IF(IsoSurfAllocated) THEN
                IsoSurf(NElems,1) = Nodes % x(nn) + &
                    temp1 * ((Nodes % x(next) - Nodes % x(nn)) / (temp1-temp2))              
                IsoSurf(NElems,2) = Nodes % y(nn) + &
                    temp1 * ((Nodes % y(next) - Nodes % y(nn)) / (temp1-temp2))

                xmin = MIN( IsoSurf(Nelems,1), xmin ) 
                xmax = MAX( IsoSurf(Nelems,1), xmax )                 
              END IF
            ELSE
              CALL Warn('PhaseChangeSolve','Wiggly Isotherm')
              WRITE(Message,*) 'nodeindexes',nodeindexes(1:Vertex)
              CALL Warn('PhaseChangeSolve',Message)
              WRITE(Message,*) 'temperature',Temperature(TempPerm(NodeIndexes))
              CALL Warn('PhaseChangeSolve',Message)
              PRINT *,'Nodes % x',Nodes % x(n)
              PRINT *,'Nodes % y',Nodes % y(n)
            END IF
          END IF          
        END DO
        
        IF ( n == 1 ) THEN
          NElems = NElems - 1
          CYCLE
        END IF
        
        IF (IsoSurfAllocated .AND. n == 2) THEN
          IF ( IsoSurf(Nelems-1,TangentDirection) > IsoSurf(Nelems,TangentDirection) ) THEN
            Temppi = IsoSurf(Nelems-1,1)
            IsoSurf(Nelems-1,1) = IsoSurf(Nelems,1)
            IsoSurf(Nelems,1) = Temppi
            Temppi = IsoSurf(Nelems-1,2)
            IsoSurf(Nelems-1,2) = IsoSurf(Nelems,2)
            IsoSurf(Nelems,2) = Temppi
          END IF
        END IF
        
      END DO 

      IF(Nelems == 0) CALL Fatal('PhaseChangeSolve','Isotherm is empty thus cannot map phase change surface') 

      IF(.NOT. IsoSurfAllocated) THEN
        ALLOCATE( IsoSurf(Nelems+1,2))
        IsoSurfAllocated = .TRUE.
        WRITE(Message,*) 'Isotherm created with number of segments',Nelems
        CALL Info('PhaseChangeSolve',Message)
        GOTO 100
      END IF      
    END IF

    area = 0.0
    volume = 0.0
    tave = 0.0
    tabs = 0.0
    volabs = 0.0
    NodeDone = .FALSE.
    MaxTempDiff = 0.0


    DO t = 1, Solver % NumberOfActiveElements 

      CurrentElement => Solver % Mesh % Elements(Solver % ActiveElements(t))
      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes

      ElementCode = CurrentElement % TYPE % ElementCode
      IF(ElementCode < 200) CYCLE
      IF(ElementCode > 203) THEN
        CALL Fatal('PhaseChangeSolve','Implemented only for elements 202 and 203!')
        CYCLE
      END IF

      Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

      TempDiff(1:n) = Temperature( TempPerm(NodeIndexes(1:n)) ) - MeltPoint
      MaxTempDiff = MAX(MaxTempDiff, MAXVAL(ABS(TempDiff(1:n)))) 

      DO nn=1,n
 
        k = SurfPerm(NodeIndexes(nn))
        IF ( NodeDone(k) ) CYCLE
        NodeDone(k) = .TRUE.

        IF( TriplePointFixed .AND. NodeIndexes(nn) == trip_node) THEN
          SurfaceMove(k) = 0.0d0
          PRINT *,'triple point fixed by construction'
          CYCLE 
        END IF

        IF ( nn == 3 ) THEN
          Surface(k)=0.5*SUM(Surface(SurfPerm(NodeIndexes(1:2))))
          Update = 0.0
        END IF

        IF(TangentDirection == 1) THEN
          xx = Nodes % x(nn)
          yy = Nodes % y(nn)
        ELSE
          xx = Nodes % y(nn)
          yy = Nodes % x(nn)
        END IF


        IF ( .NOT. Newton ) THEN          
          
          ! Find the contour element that has the x-coordinate in closest to that of the
          ! free surface

          Eps = 1.0d-6 * ( xmax - xmin )

          dxmin = HUGE(dxmin)          
          dymin = HUGE(dymin)
          stat = .FALSE.

          DO i=1,Nelems-1,2

            IF ( (xx > IsoSurf(i,TangentDirection) - Eps) .AND. (xx < IsoSurf(i+1,TangentDirection) + Eps)) THEN
              dxmin = 0.0
              d = MIN( ABS(yy - IsoSurf(i,NormalDir)), ABS(yy - IsoSurf(i+1,NormalDir)) )              

              ! Punish for overlapping the boundaries
              d = d + MAX(0.0d0, IsoSurf(i,TangentDirection) - xx)
              d = d + MAX(0.0d0, xx - IsoSurf(i+1,TangentDirection) )

              IF(d <= dymin) THEN
                stat = .TRUE.
                dymin = d
                imin = i
              END IF
            END IF

            IF(.NOT. stat) THEN
              d = MIN( ABS(xx - IsoSurf(i,TangentDirection)), ABS(xx - IsoSurf(i+1,TangentDirection)) )
              IF (d <= dxmin) THEN
                dxmin = d
                imin = i
              END IF
            END IF
          END DO

          i = imin
          
          ! There may be a problem if the boundary cannot be mapped on an isotherm
          IF (.NOT. stat) THEN
            IF(dxmin > 1.0d-2* ABS(IsoSurf(i,TangentDirection)- IsoSurf(i+1,TangentDirection))) THEN
              CALL Warn('PhaseChangeSolve','Isotherm error?')
              WRITE(Message,*) 'Nodeindexes',NodeIndexes(nn)
              CALL Warn('PhaseChangeSolve',Message)
              WRITE(Message,*) 'x:',xx,' y:',yy
              CALL Warn('PhaseChangeSolve',Message)
              WRITE(Message,*) 'dxmin:',dxmin,' dymin:',dymin
              CALL Warn('PhaseChangeSolve',Message)
            END IF
          END IF

          IF ( ABS( IsoSurf(i+1,TangentDirection) - IsoSurf(i,TangentDirection) ) > AEPS ) THEN
            Update = IsoSurf(i,NormalDir) + ( xx - IsoSurf(i,TangentDirection) ) * &
                ( IsoSurf(i+1,NormalDir) - IsoSurf(i,NormalDir) ) / &
                ( IsoSurf(i+1,TangentDirection) - IsoSurf(i,TangentDirection) ) - yy
          ELSE
            Update = 0.5d0 * ( IsoSurf(i,NormalDir) + IsoSurf(i+1,NormalDir) ) - yy
          END IF

        END IF
              
        TTemp = Temperature( TempPerm(NodeIndexes(nn)) )        
       
        IF ( Newton ) THEN
          dTdz = TTemp - PrevTemp(k) 
          IF ( ABS(dTdz) < AEPS ) THEN
            CALL Warn( 'PhaseChangeSolve', 'Very small temperature update.' )
            dTdz = 1
          END IF
          Update = SurfaceMove(k) * ( MeltPoint - TTemp ) / dTdz
        END IF
         
        ! This enforcing is rather than by setting meltpoint to triple point temperature
        ! IF ( NodeIndexes(nn) == Trip_node ) Update = 0
     
        PrevTemp(k) = TTemp
        SurfaceMove(k) = Update
      END DO

      IntegStuff = GaussPoints( CurrentElement )      
      DO i=1,IntegStuff % n        
        
        u = IntegStuff % u(i)
        v = IntegStuff % v(i)
        w = IntegStuff % w(i)
        
        stat = ElementInfo( CurrentElement, Nodes, u, v, w, detJ, Basis, dBasisdx )
        
        s = IntegStuff % s(i) * detJ
        
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          s = s * SUM(Basis(1:n) * Nodes % x(1:n)) * 2.0 * PI
        END IF
        
        area = area + S
        volume = volume + S * SUM(Basis(1:n) * SurfaceMove(SurfPerm(NodeIndexes(1:n))))
        tave = tave + S * SUM(Basis(1:n) * TempDiff(1:n) )
        volabs = volabs + S * SUM(Basis(1:n) * ABS(SurfaceMove(SurfPerm(NodeIndexes(1:n)))))
        tabs = tabs + S * SUM(Basis(1:n) * ABS(TempDiff(1:n)) )
      END DO
      
    END DO


!    IF(.NOT. UseTAverage .AND. dt < SteadyDt) THEN
!      LocaRelax = Relax * GetCReal( Solver % Values,'Temperature Relaxation Factor')
!    ELSE
      LocalRelax = Relax 
!    END IF

    ! There are several different acceleration methods which are mainly inactive
    tave = tave / area
    tabs = tabs / area
    volume = volume / area
    volabs = volabs / area
    
    i = ListGetInteger(Solver % Values,'Lumped Newton After Iterations', Stat)
    IF(Stat .AND. SubroutineVisited > i) THEN
      
      j = ListGetInteger(Solver % Values,'Lumped Newton Mode', Stat)
      SELECT CASE( j ) 
      CASE( 1 )
        cvol = 0.5*(prevtave+tave)/(prevtave-tave)
        
      CASE( 2 )
        cvol = 0.5*(prevvolabs+volabs)/(prevvolabs-volabs)
        
      CASE( 3 )
        cvol = 0.5*(prevtabs+tabs)/(prevtabs-tabs)
        
      CASE DEFAULT
        cvol = 0.5*(prevvolume+volume)/(prevvolume-volume)

      END SELECT
      
      IF(cvol < 0.0) THEN
        cvol = 1.0
        ccum = 1.0
      END IF
      
      clim = ListGetConstReal(Solver % Values,'Lumped Newton Limiter', Stat)
      IF(.NOT. Stat) clim = 100.0
      cvol = MIN(clim,cvol)
      cvol = MAX(1.0/clim,cvol)
      
      ccum = ccum * cvol
      
      WRITE(Message,*) 'Lumped Newton relaxation: ', ccum
      CALL Info('PhaseChangeSolve',Message)
      
      LocalRelax = LocalRelax * ccum
    END IF

    SurfaceMove = LocalRelax * SurfaceMove
    Surface = Surface + SurfaceMove
    
    dpos = SurfaceMove(SurfPerm(Trip_node))

    MaxUpdate = MAXVAL(ABS(SurfaceMove)) / MAXVAL(ABS(Surface))

    IF ( ABS(MaxUpdate) < NewtonAfterTol ) Newton = .TRUE.    
    WRITE(Message,*) 'Maximum surface update: ', MaxUpdate
    CALL Info('PhaseChangeSolve',Message)

    WRITE(Message,*) 'Maximum temperature difference: ', MaxTempDiff
    CALL Info('PhaseChangeSolve',Message)

    Norm = SQRT( SUM( Surface**2 ) / SIZE( Surface ) )
    WRITE( Message, * ) 'Result Norm     : ',Norm
    CALL Info( 'PhaseChangeSolve', Message, Level=4 )
    Solver % Variable % Norm = Norm
    
    prevvolume = volume
    prevtave = tave
    prevtabs = tabs
    prevvolabs = volabs
    
    IF(IsoSurfAllocated) THEN
      IsoSurfAllocated = .FALSE.
      DEALLOCATE(IsoSurf)
    END IF
  END IF

!--------------------------------------------------------------------
! Steady algorithm end
!--------------------------------------------------------------------

  
200 CALL ListAddConstReal(Model % Simulation,'res: Triple point temperature',Trip_temp)
  CALL ListAddConstReal( Model % Simulation,'res: triple point movement',dpos)
  CALL ListAddConstReal(Model % Simulation,'res: Pull Position',pos0)       

  IF(PullControl) THEN
    IF(NormalDir == 1) THEN
      CALL ListAddConstReal( Model % Simulation,'res: Pull Velocity 1',UPull(NormalDir))
    ELSE         
      CALL ListAddConstReal( Model % Simulation,'res: Pull Velocity 2',UPull(NormalDir))        
    END IF
  END IF

  FirstTime = .FALSE.
  
!------------------------------------------------------------------------------
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE VelocityLocalMatrix( StiffMatrix, MassMatrix, ForceVector,&
      Element, nCoord, Nodes )
        
    ! external variables:
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:), ForceVector(:)      
    INTEGER :: nCoord
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    
    ! internal variables:
    TYPE(Nodes_t) :: PNodes
    TYPE(Element_t), POINTER :: Parent      
    REAL(KIND=dp) :: Basis(3*nCoord),dBasisdx(3*nCoord,3), &
        X,Y,Z,U,V,W,S,detJ, TGrad(3,3),Flux,pu,pv,pw,pull, LocalHeat, &
        NodalTemp(3*nCoord), xx(10),yy(10),zz(10),  NodalSurf(nCoord), NodalNormal(2,nCoord)
    REAL(KIND=dp) :: Velo, Ny, Normal(3), StabFactor, StabCoeff, DerSurf, xcoord
    LOGICAL :: Stat
    INTEGER :: i,j,k,l,t,p,q, n, pn
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Solver % Values,'Velocity Smoothing Factor')

    Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    NodalSurf(1:n) = Surface( SurfPerm(NodeIndexes) )
   
    IF(AverageNormal) THEN
      NodalNormal(1,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n))-1)
      NodalNormal(2,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n)))
    END IF

    IF(TransientAlgo) THEN
      ALLOCATE( PNodes % x(10), PNodes % y(10), PNodes % z(10) )
    END IF

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
     
!      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
!        s = s * xcoord
!      END IF
      
      IF(AverageNormal) THEN
        Normal(1) = SUM( Basis(1:n) * NodalNormal(1,1:n))
        Normal(2) = SUM( Basis(1:n) * NodalNormal(2,1:n))
        Normal(3) = 0.0d0
      ELSE
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )         
      END IF

      LocalHeat = SUM(Basis(1:n) * LatentHeat(1:n))
      StabCoeff = StabFactor * LocalHeat * Density 
      
 
      IF(UseLoads) THEN
        ! do nothing, loads already computed

      ELSE 
        ! Compute the flux from normal derivatives
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
            IF(l == 2) CALL Fatal('PhaseChangeSolve','Both materials cannot be solid!')
            l = 2
          ELSE 
            IF(l == 1) CALL Fatal('PhaseChangeSolve','Both materials cannot be liquid!')
            l = 1
          END IF
          
          pn = Parent % TYPE % NumberOfNodes
          k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
          
          Conductivity(1:pn) = ListGetReal( Model % Materials(k) % Values, &
              'Heat Conductivity', pn, Parent % NodeIndexes )                     
          NodalTemp(1:pn) = Temperature( TempPerm(Parent % NodeIndexes) )
          
          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
          
          !           Calculate the basis functions for the parent element:
          !           -----------------------------------------------------
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
      DO p=1,n        
        DO q=1,n
          StiffMatrix(p,q) = StiffMatrix(p,q) + s * Basis(p) * Basis(q) * &
              Normal(NormalDir) * Density * LocalHeat            
          IF(NodeIndexes(p) /= Trip_node) THEN
            StiffMatrix(p,q) = StiffMatrix(p,q) + &
                s * StabCoeff * dBasisdx(q,TangentDirection) * dBasisdx(p,TangentDirection)            
          END IF

          ! BC for tipple node 
!          IF(NodeIndexes(p) == Trip_node) THEN
!            StiffMatrix(p,q) = StiffMatrix(p,q) - &
!                StabCoeff * Basis(p) * dBasisdx(q,TangentDirection)            
!          END IF
        END DO
                
        ! transient part of heat flux
        IF(.NOT. (UseLoads)) THEN
          ForceVector(p) = ForceVector(p) - s * Basis(p) * Flux 
        END IF
        
      END DO

    END DO
    
    IF(TransientAlgo) DEALLOCATE( PNodes % x, PNodes % y, PNodes % z )

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
        X,Y,Z,U,V,W,S,detJ,pull, NodalVelo(nCoord), xcoord, NodalNormal(2,nCoord)
    REAL(KIND=dp) :: Velo, StabFactor, StabCoeff
    LOGICAL :: Stat
    INTEGER :: i,j,k,l,t,p,q, n, pn
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

    ForceVector = 0.0d0
    StiffMatrix = 0.0d0
    MassMatrix  = 0.0d0
    n = nCoord

    StabFactor = ListGetConstReal(Solver % Values,'Surface Smoothing Factor')
    StabCoeff = StabFactor / dt

    Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    NodalVelo(1:n) = SurfaceMove( SurfPerm(NodeIndexes) )

    IF(AverageNormal) THEN
      NodalNormal(1,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n))-1)
      NodalNormal(2,1:n) = Normals(2*NormalsPerm(NodeIndexes(1:n)))
    END IF


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
      
      IF(AverageNormal) THEN
        Normal(1) = SUM( Basis(1:n) * NodalNormal(1,1:n) )
        Normal(2) = SUM( Basis(1:n) * NodalNormal(2,1:n) ) 
        Normal(3) = 0.0d0
      ELSE
        Normal = NormalVector( Element, Nodes, u, v, .TRUE. )         
      END IF
      
      Velo = SUM( Basis(1:n) * NodalVelo(1:n)) + Upull(NormalDir)        
      Velo = SpeedUp * Velo 

      DO p=1,n        
        DO q=1,n
          MassMatrix(p,q) = MassMatrix(p,q) + s * Basis(p) * Basis(q)  
          StiffMatrix(p,q) = StiffMatrix(p,q) + &
              s * StabCoeff * dBasisdx(q,TangentDirection) * dBasisdx(p,TangentDirection)            
          
          ! BC for tipple node 
          IF(NodeIndexes(p) == Trip_node) THEN
            StiffMatrix(p,q) = StiffMatrix(p,q) - &
                StabCoeff * Basis(p) * dBasisdx(q,TangentDirection)            
          END IF
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

      IF(TangentDirection /= 1) CALL Fatal('PhaseChangeSolve',&
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
  END SUBROUTINE PhaseChangeSolve
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  FUNCTION MeltingHeat(Model, Node, t) RESULT(Flux)
!-------------------------------------------------------------------------------
! This subroutine computes the heat flux resulting from solidification 
! This 
!-------------------------------------------------------------------------------
  USE Types
  USE Lists
  USE ElementDescription
  IMPLICIT NONE
!-------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER:: Node
  REAL (KIND=dp):: Flux,t
!-------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: NormalSol
  INTEGER:: k,n,i
  INTEGER, POINTER :: NodeIndexes(:), NormalsPerm(:) 
  REAL (KIND=dp):: NodeLatentHeat, Density, NormalPull, u, v
  REAL (KIND=dp):: UPull(3) = (/ 0,0,0 /), Normal(3), ElemLatentHeat(4)
  REAL (KIND=dp), POINTER :: Normals(:)
  LOGICAL:: stat, NormalExist = .FALSE., Visited = .FALSE.
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t), POINTER :: CurrentElement, Parent
  
!------------------------------------------------------------------------------

  SAVE NormalExist, Normals, NormalsPerm, Nodes


  IF(.NOT. Visited) THEN
    NormalSol  => VariableGet( Model % Variables, 'Normals',ThisOnly=.TRUE. )
    IF(ASSOCIATED(NormalSol)) THEN
      NormalsPerm => NormalSol % Perm
      Normals => NormalSol % Values
      NormalExist = .TRUE.
    ELSE
      n = Model % Mesh % MaxElementNodes  
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
    END IF
    Visited = .TRUE.
  END IF

  CurrentElement => Model % CurrentElement
  NodeIndexes => CurrentElement % NodeIndexes
  n = CurrentElement % TYPE % NumberOfNodes
  
  k = ListGetInteger(Model % Bodies(CurrentElement % BodyId) % Values,'Material')
  ElemLatentHeat(1:n) = ListGetReal( Model % Materials(k) % Values, 'Latent Heat', n, NodeIndexes )
  
  DO i=1,n
    IF(NodeIndexes(i) == Node) EXIT
  END DO
  IF(NodeIndexes(i) /= Node) CALL Fatal('PhaseProcs','Node not found')
  NodeLatentHeat = ElemLatentHeat(i)
  
  UPull = 0.0
  UPull(1) = ListGetConstReal(Model % Simulation,'res: Pull Velocity 1',stat) 
  IF(.NOT. stat) UPull(1) = ListGetConstReal(Model % Materials(k) % Values,'Convection Velocity 1',stat) 

  UPull(2) = ListGetConstReal(Model % Simulation,'res: Pull Velocity 2',stat) 
  IF(.NOT. stat) UPull(2) = ListGetConstReal(Model % Materials(k) % Values,'Convection Velocity 2',stat) 
 
  Parent => CurrentElement % BoundaryInfo % Left
  k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
  IF (ListGetLogical(Model % Materials(k) % Values, 'Solid', stat)) THEN
    Density = ListGetConstReal( Model % Materials(k) % Values, 'Density' )
  ELSE   
    Parent => CurrentElement % BoundaryInfo % Right
    k = ListGetInteger(Model % Bodies(Parent % BodyId) % Values,'Material')
    Density = ListGetConstReal( Model % Materials(k) % Values, 'Density' )
  END IF

  IF(.NOT. ASSOCIATED(Parent)) CALL Fatal('PhaseProcs','Parent not found')
!------------------------------------------------------------------------------

  IF( NormalExist ) THEN
    Normal(1) = Normals( 2*NormalsPerm( Node ) - 1 )
    Normal(2) = Normals( 2*NormalsPerm( Node ))
    Normal(3) = 0.0d0
  ELSE
    Nodes % x(1:n) = Model % Nodes % x(NodeIndexes)
    Nodes % y(1:n) = Model % Nodes % y(NodeIndexes)
    Nodes % z(1:n) = Model % Nodes % z(NodeIndexes)
    
    ! For line segments the normal in the center is usually sufficient
    u = 0.0d0
    v = 0.0d0
    
    ! If inner boundary, Normal Target Body should be defined for the boundary 
    ! (if not, material density will be used to determine then normal direction 
    ! and should be defined for bodies on both sides):

    Normal = NormalVector( CurrentElement, Nodes, u, v, .TRUE. )
  END IF

  NormalPull = SUM( Normal * UPull )

  Flux = NodeLatentHeat * Density * NormalPull

!------------------------------------------------------------------------------
END FUNCTION MeltingHeat
!------------------------------------------------------------------------------
