/*****************************************************************************/
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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science 
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Date:    5.12.2008
! *  Edited:  16.12.2008 
! *
! *  Copyright 2008, CSC - IT Center for Science Ltd.
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!>  Lagrangian steady-state phase change solver for the liquid/solid interface.
!>  This phase change solver is based on the finding of a isotherm and mapping the 
!>  interface to the new isotherm. Also the previous temperature gradient may be 
!>  used to estimate the distance from the interface.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE SteadyPhaseChange( Model,Solver,dt,TransientSimulation )
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
  TYPE(Element_t), POINTER :: Element, Parent, Parent2
  TYPE(Variable_t), POINTER :: SurfSol, TempSol, HelpSol
  TYPE(Nodes_t) :: Nodes
  TYPE(GaussIntegrationPoints_t) :: IntegStuff  
  TYPE(ValueList_t), POINTER :: Material, BC
  TYPE(Solver_t), POINTER :: PSolver 
  TYPE(ValueList_t), POINTER :: Params
  
  REAL(KIND=dp) :: Normal(3), u, v, w, &
       Density, Update, MaxUpdate, MaxTempDiff, Relax, LocalRelax, &
       surf, xx, yy, r, detJ, Temp, MeltPoint, &
       Temp1,Temp2,Temppi,dxmin,dymin, xmin, xmax, d, HeatCond, s, tave, prevtave, cvol, clim, ccum=1, &
       tabs, prevtabs, volabs, prevvolabs, CoordMin(3), CoordMax(3), RelativeChange, &
       dTdz, ttemp, NewtonAfterTol, area, volume, prevvolume, Eps, &
       Norm, maxds, maxds0, ds, FluxCorrect = 1.0, dpos, Width, &
       pos0=0.0, prevpos0, Coeff, OrigNorm, Trip_Temp, SpeedUp, PrevNorm, &
       MaxAngle, MaxSurfaceMove, MaxSurface, x1, x2, y1, y2
  REAL(KIND=dp), POINTER :: Surface(:), SurfaceMove(:), Temperature(:), ForceVector(:),  &
       x(:), y(:), z(:), Basis(:), NodalTemp(:), &
       TempDiff(:),Weights(:),NewY(:)
  REAL (KIND=dp), ALLOCATABLE :: PrevTemp(:), IsoSurf(:,:)
  
  INTEGER :: i,j,k,t,n,nn,pn,DIM,kl,kr,l, Trip_node, NoBNodes, istat, &
       NElems,ElementCode,Next,Vertex,ii,imin,NewtonAfterIter,Node, iter, LiquidInd, Visited = -1, &
       SubroutineVisited = 0, NormalDirection, TangentDirection, CoordMini(3), CoordMaxi(3), &
       Axis_node, LiquidBody, SolidBody, NoPhaseElements
  INTEGER, POINTER :: Indexes(:),TempPerm(:),SurfPerm(:),PhaseElements(:)

  LOGICAL :: Stat, FirstTime = .TRUE., Newton = .FALSE., Debug, &
       IsoSurfAllocated, AllocationsDone = .FALSE., &
      UseLoads, SurfaceVelocitySet = .FALSE., TriplePointFixed, &
      GotSolid, GotLiquid, Active, GotIt
  LOGICAL, POINTER :: NodeDone(:), BoundaryMarker(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, str

  SAVE FirstTime, Trip_node, NoBNodes, SubroutineVisited, prevpos0, &
      PrevTemp, Newton, ccum, NormalDirection, TangentDirection, MeltPoint, &
      ForceVector, FluxCorrect, NodeDone, Eps,  &
      Visited, Nodes, NodalTemp, Width, &
      TempDiff, AllocationsDone, &
      x, y, z, Basis, norm, Axis_node, &
      Weights, SurfaceMove,  GotLiquid, GotSolid, LiquidBody, SolidBody, &
      SurfaceVelocitySet, CoordMax, CoordMin, CoordMaxi, CoordMini, &
      BoundaryMarker, IsoSurfAllocated, PhaseElements, NoPhaseElements, &
      prevtave, prevtabs, prevvolabs, prevvolume
  
!------------------------------------------------------------------------------

  CALL Info('SteadyPhaseChange',   '--------------------------------------------')
  CALL Info('SteadyPhaseChange',   'Using steady algorithm to find the isotherm')      
  CALL Info('SteadyPhaseChange',   '--------------------------------------------')

  SubroutineVisited = SubroutineVisited + 1

  DIM = CoordinateSystemDimension()
  IF(DIM /= 2) THEN
    CALL Fatal('SteadyPhaseChange','Implemented only in 2D')
  END IF

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  Params => GetSolverParams()

  PSolver => Solver
  NULLIFY(SurfSol)
  NULLIFY(Surface)
  SurfSol  => Solver % Variable
  IF(ASSOCIATED(SurfSol)) THEN
     Surface  => SurfSol % Values
     SurfPerm => SurfSol % Perm
     IF(.NOT. ASSOCIATED (Surface) .OR. ALL(SurfPerm <= 0) ) THEN
        CALL Fatal('SteadyPhaseChange','Surface field needed for Phase Change')
     END IF
     VariableName = ComponentName(Solver % Variable)
  ELSE
    VariableName = ListGetString( Params, 'Surface Variable', Stat )
    IF(.NOT. Stat) VariableName = 'Surface'
    SurfSol => VariableGet( Solver % Mesh % Variables, VariableName )
    Surface  => SurfSol % Values
    SurfPerm => SurfSol % Perm
    IF(.NOT. ASSOCIATED (Surface) .OR. ALL(SurfPerm <= 0)) THEN
       CALL Fatal('SteadyPhaseChange','Surface field needed for Phase Change')
    END IF
  END IF

! The change in surface variable
  HelpSol => VariableGet( Solver % Mesh % Variables, &
       TRIM(VariableName )//'Diff' )
  IF(ASSOCIATED(HelpSol)) THEN
    SurfaceMove => HelpSol % Values
  ELSE
    ALLOCATE( SurfaceMove(SIZE(Surface)), STAT=istat)
    IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error 1.' )           
    SurfaceMove = 0.0d0
    CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
        PSolver,TRIM(VariableName )//'Diff',1,SurfaceMove,SurfPerm)
  END IF

! variable of which the isotherm is used
  VariableName = ListGetString( Params, 'Phase Change Variable', Stat )
  IF( .NOT. Stat ) VariableName = 'Temperature'

  TempSol => VariableGet( Solver % Mesh % Variables, VariableName )
  TempPerm => TempSol % Perm
  Temperature => TempSol % Values
  IF(.NOT. ASSOCIATED (Temperature) ) THEN
     CALL Fatal('SteadyPhaseChange','Temperature field needed for Phase Change')
  END IF
  IF(ALL(TempPerm <= 0) ) THEN
     CALL Fatal('SteadyPhaseChange','Temperature field has zero Perm vector!')
  END IF

!---------------------------------------------------------------------------------
  Relax = GetCReal( Params,  & 
       'Nonlinear System Relaxation Factor', stat )
  IF ( .NOT. stat ) Relax = 1.0d0

!---------------------------------------------------------------------------------
! The first time the main axis of the free surface is determined
! and some permanent vectors related to the surface are allocated.
!---------------------------------------------------------------------------------
  
  IF(FirstTime) THEN
     ALLOCATE(BoundaryMarker(SIZE(SurfPerm)), STAT=istat)
     IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error 2.' )           

     BoundaryMarker = .FALSE.
     NoPhaseElements = 0
     DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        n = GetElementNOFNodes()
        IF ( GetElementFamily() == 1 ) CYCLE
        
        BC => GetBC()
        IF( .NOT. GetLogical( BC, 'Phase Change', GotIt ) ) CYCLE

        BoundaryMarker( Element % NodeIndexes ) = .TRUE.
        NoPhaseElements = NoPhaseElements + 1
     END DO    
     NoBNodes = COUNT ( BoundaryMarker )

     WRITE(Message,'(A,T35,I12)') 'Number of interface nodes:',NoBNodes
     CALL Info('SteadyPhaseChange',Message)

     IF( NoBNodes == 0 ) THEN
        CALL Warn('SteadyPhaseChange','There is no Phase Change boundary')
     END IF

     WRITE(Message,'(A,T35,I12)') 'Number of interface elements:',NoPhaseElements
     CALL Info('SteadyPhaseChange',Message)
    
     ALLOCATE(PhaseElements(NoPhaseElements), STAT=istat) 
     IF ( istat /= 0 ) CALL Fatal( 'PhaseChangeSolver', 'Memory allocation error 3.' )           
     PhaseElements = 0

     NoPhaseElements = 0
     DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        n = GetElementNOFNodes()
        IF ( GetElementFamily() == 1 ) CYCLE
        
        BC => GetBC()
        IF( .NOT. GetLogical( BC, 'Phase Change', GotIt ) ) CYCLE

        NoPhaseElements = NoPhaseElements + 1
        PhaseElements(NoPhaseElements) = t
     END DO    

     CoordMax = -HUGE(CoordMax)
     CoordMin = HUGE(CoordMin)

     DO k=1, Model % Mesh % NumberOfNodes
        IF( .NOT. BoundaryMarker(k) ) CYCLE

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
     NormalDirection = j

     ! Direction of maximum change
     j = 1
     DO i=1,DIM
        IF(CoordMax(i)-CoordMin(i) > CoordMax(j)-CoordMin(j)) THEN
           j = i
        END IF
     END DO
     Width = CoordMax(j) - CoordMin(j)
     TangentDirection = j


     WRITE(Message,'(A,T35,ES12.4)') 'Width of the interface',Width
     CALL Info('SteadyPhaseChange',Message)

     ! Find triple point:
     !-------------------
     Trip_node = CoordMaxi(TangentDirection)
     Axis_node = CoordMini(TangentDirection)
     
     WRITE(Message,'(A,T35,I12)') 'Index of the triple point: ',Trip_node
     CALL Info('SteadyPhaseChange',Message)
     
     WRITE(Message,'(A,T35,I12)') 'Index of the axis point: ',Axis_node
     CALL Info('SteadyPhaseChange',Message)


     ! Look if Liquid and Solid are determined in bodies
     !--------------------------------------------------
     GotLiquid = .FALSE.
     GotSolid = .FALSE.
     DO i=1,Model % NumberOfMaterials
        GotLiquid = GotLiquid .AND. &
             ListGetLogical( Model % Bodies(i) % Values,'Liquid', GotIt ) 
        GotSolid = GotSolid .AND. &
             ListGetLogical( Model % Bodies(i) % Values,'Solid', GotIt ) 
     END DO

     ! If not determine them by parenthood so that lower is liquid
     !------------------------------------------------------------
     SolidBody = 0
     LiquidBody = 0

     IF( .NOT. (GotLiquid .AND. GotSolid) ) THEN
        DO t=1, Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(t)
           
           n = GetElementNOFNodes()
           IF ( GetElementFamily() == 1 ) CYCLE
           
           BC => GetBC()
           IF( .NOT. GetLogical( BC, 'Phase Change', GotIt ) ) CYCLE
           
           Parent => Element % BoundaryInfo % Left
           Parent2 => Element % BoundaryInfo % Right
           IF( MAXVAL(Model % Mesh % Nodes % y(Parent % NodeIndexes)) > &
                MAXVAL(Model % Mesh % Nodes % y(Parent2 % NodeIndexes)) ) THEN
              IF(.NOT. GotSolid) SolidBody = Parent % BodyId
              IF(.NOT. GotLiquid) LiquidBody = Parent2 % BodyId
           ELSE
              IF(.NOT. GotSolid) SolidBody = Parent2 % BodyId
              IF(.NOT. GotLiquid) LiquidBody = Parent % BodyId
           END IF
           EXIT
        END DO

        WRITE(Message,'(A,T35,I12)') 'Body index for solid: ',SolidBody
        CALL Info('SteadyPhaseChange',Message)

        WRITE(Message,'(A,T35,I12)') 'Body index for liquid: ',LiquidBody
        CALL Info('SteadyPhaseChange',Message)
     END IF

     n = Solver % Mesh % MaxElementNodes  

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
          x(n), y(n), z(n), Basis(n), NodalTemp(n), &
          TempDiff(n), &
          PrevTemp(SIZE(Surface)), &
          STAT=istat)
     IF ( istat /= 0 ) CALL Fatal( 'SteadyPhaseChange', 'Memory allocation error 4.' )     

     Nodes % x = 0.0d0
     Nodes % y = 0.0d0
     Nodes % z = 0.0d0
     PrevTemp = 0.0d0
     
     ALLOCATE( NodeDone( SIZE(Surface) ), STAT=istat)
     IF (istat /= 0 ) CALL Fatal( 'SteadyPhaseChange', 'Memory allocation error 5.' )          
     
     AllocationsDone = .TRUE.    
  END IF

  Trip_Temp =  Temperature( TempPerm(Trip_node) )    

  i =  ListGetInteger( Params,'Passive Steps',Stat)
  IF( i >= SubroutineVisited) GOTO 200

!----------------------------------------------------------------------------

  NewtonAfterIter = ListGetInteger( Params, &
       'Nonlinear System Newton After Iterations', stat )    
  IF ( stat .AND. SubroutineVisited > NewtonAfterIter ) Newton = .TRUE.
  
  NewtonAfterTol = ListGetConstReal( Params, &
       'Nonlinear System Newton After Tolerance', stat )
  
  IF(Newton) THEN
    CALL Info( 'SteadyPhaseChange','Steady state newton formulation', Level=4 )
  ELSE
    CALL Info( 'SteadyPhaseChange','Steady state isoterm formulation', Level=4 )
  END IF

! Find the melting point:
!-----------------------

  TriplePointFixed = ListGetLogical( Params,'Triple Point Fixed',stat)
  IF(.NOT. Stat) TriplePointFixed = .FALSE.

  IF( TriplePointFixed ) THEN
    MeltPoint = Trip_Temp
    WRITE(Message,'(A,T35,ES12.4)') 'Melting point set: ',MeltPoint
    CALL Info('SteadyPhaseChange',Message)        
  ELSE
    DO k=1, Model % NumberOfMaterials
      MeltPoint = GetCReal( Model % Materials(k) % Values, &
          'Melting Point', GotIt )
      IF(GotIt) EXIT
    END DO
    IF( GotIt ) THEN
      WRITE(Message,'(A,T35,ES12.4)') 'Melting point found: ',MeltPoint
      CALL Info('SteadyPhaseChange',Message)        
    ELSE
      CALL Info('SteadyPhaseChange','Could not find melting point in any material!')        
    END IF
  END IF


! Create the isotherm for T=T_0:
!------------------------------

  IF(.NOT. Newton) THEN
    CALL CreateIsotherm()
  END IF


  area = 0.0
  volume = 0.0
  tave = 0.0
  tabs = 0.0
  volabs = 0.0
  NodeDone = .FALSE.
  MaxTempDiff = 0.0
  
  DO t=1, Solver % Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    
    n = GetElementNOFNodes()
    IF ( GetElementFamily() == 1 ) CYCLE
    
    BC => GetBC()
    IF( .NOT. GetLogical( BC, 'Phase Change', GotIt ) ) CYCLE
    
    Indexes => Element % NodeIndexes
    
    ElementCode = Element % TYPE % ElementCode
    IF(ElementCode < 200 .OR. ElementCode > 203) THEN
      CALL Fatal('PhaseChangeSolve','Implemented only for elements 202 and 203!')
      CYCLE
    END IF
    
    Nodes % x(1:n) = Solver % Mesh % Nodes % x(Indexes)
    Nodes % y(1:n) = Solver % Mesh % Nodes % y(Indexes)
    Nodes % z(1:n) = Solver % Mesh % Nodes % z(Indexes)
    
    TempDiff(1:n) = Temperature( TempPerm(Indexes(1:n)) ) - MeltPoint
    MaxTempDiff = MAX(MaxTempDiff, MAXVAL(ABS(TempDiff(1:n)))) 
    
    DO nn=1,n
      
      k = SurfPerm(Indexes(nn))
      IF ( NodeDone(k) ) CYCLE
      NodeDone(k) = .TRUE.
      
      IF( TriplePointFixed .AND. Indexes(nn) == trip_node) THEN
        SurfaceMove(k) = 0.0d0
        CALL Info('SteadyPhaseChange','Triple point position fixed')
        CYCLE 
      END IF
      
      ! For 2nd order set the middle node to be the mean 
      IF ( nn == 3 ) THEN
        SurfaceMove(k) = 0.5*SUM(SurfaceMove(SurfPerm(Indexes(1:2))))
        Update = 0.0
      END IF
      
      
      TTemp = Temperature( TempPerm(t) )
      
      IF ( Newton ) THEN
        dTdz = TTemp - PrevTemp(k) 
        IF ( ABS(dTdz) < AEPS ) THEN
          CALL Warn( 'SteadyPhaseChange', 'Very small temperature update.' )
          dTdz = 1
        END IF
        Update = SurfaceMove(k) * ( MeltPoint - TTemp ) / dTdz
      ELSE           
        ! Find the the contour element that has the x-coordinate in closest to the that of the
        ! free surface
        
        Eps = 1.0d-6 * ( xmax - xmin )
        
        IF(TangentDirection == 1) THEN
          xx = Nodes % x(nn)
          yy = Nodes % y(nn)
        ELSE
          xx = Nodes % y(nn)
          yy = Nodes % x(nn)
        END IF
        
        dxmin = HUGE(dxmin)          
        dymin = HUGE(dymin)
        stat = .FALSE.
        
        DO i=1,Nelems-1,2
          
          x1 = IsoSurf(i,TangentDirection)
          x2 = IsoSurf(i+1,TangentDirection)
          y1 = IsoSurf(i,NormalDirection)
          y2 = IsoSurf(i+1,NormalDirection)
          
          ! If node is in interval take the closest isotherm
          IF ( (xx > x1 - Eps) .AND. (xx < x2 + Eps)) THEN
            dxmin = 0.0
            d = MIN( ABS(yy - y1), ABS(yy - y2) )              
            
            ! Punish for overlapping the boundaries
            d = d + MAX(0.0d0, x1 - xx)
            d = d + MAX(0.0d0, xx - x2 )
            
            IF(d <= dymin) THEN
              stat = .TRUE.
              dymin = d
              imin = i
            END IF
          END IF
          
          ! If point not yet found in line segments check for close visinity              
          IF(.NOT. stat) THEN
            d = MIN( ABS(xx - x1), ABS(xx - x2) )
            IF (d <= dxmin) THEN
              dxmin = d
              imin = i
            END IF
          END IF
        END DO
        
        i = imin
        x1 = IsoSurf(i,TangentDirection)
        x2 = IsoSurf(i+1,TangentDirection)
        y1 = IsoSurf(i,NormalDirection)
        y2 = IsoSurf(i+1,NormalDirection)
        
        ! There may be a problem if the boundary cannot be mapped on an isotherm
        IF (.NOT. stat) THEN
          IF(dxmin > 1.0d-2* ABS(x1 - x2)) THEN
            CALL Warn('SteadyPhaseChange','Isotherm error?')
            WRITE(Message,*) 'Nodeindexes',Indexes(nn)
            CALL Warn('SteadyPhaseChange',Message)
            WRITE(Message,*) 'x:',xx,' y:',yy
            CALL Warn('SteadyPhaseChange',Message)
            WRITE(Message,*) 'dxmin:',dxmin,' dymin:',dymin
            CALL Warn('SteadyPhaseChange',Message)
          END IF
        END IF
        
        IF ( ABS( x2 - x1 ) > AEPS ) THEN
          Update = ( y1 - yy ) + &
              ( xx - x1 ) * ( y2 - y1 ) / ( x2 - x1 ) 
        ELSE
          Update = 0.5_dp * ( y1 + y2 ) - yy
        END IF
        
      END IF
      
            
      ! This enforcing is rather than by setting meltpoint to triple point temperature
      ! IF ( Indexes(nn) == Trip_node ) Update = 0
      
      PrevTemp(k) = TTemp
      SurfaceMove(k) = Update
    END DO

        
    IntegStuff = GaussPoints( Element )      
    DO i=1,IntegStuff % n        
      
      u = IntegStuff % u(i)
      v = IntegStuff % v(i)
      w = IntegStuff % w(i)
      
      stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
      
      s = IntegStuff % s(i) * detJ
      
      IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
        s = s * SUM(Basis(1:n) * Nodes % x(1:n)) * 2.0 * PI
      END IF
      
      area = area + S
      volume = volume + S * SUM(Basis(1:n) * SurfaceMove(SurfPerm(Indexes(1:n))))
      tave = tave + S * SUM(Basis(1:n) * TempDiff(1:n) )
      volabs = volabs + S * SUM(Basis(1:n) * ABS(SurfaceMove(SurfPerm(Indexes(1:n)))))
      tabs = tabs + S * SUM(Basis(1:n) * ABS(TempDiff(1:n)) )
    END DO
  END DO
  
  
  LocalRelax = Relax 
  
  ! There are several different acceleration methods which are mainly inactive
  tave = tave / area
  tabs = tabs / area
  volume = volume / area
  volabs = volabs / area
  
  i = ListGetInteger(Params,'Lumped Acceleration After Iterations', Stat)
  
  IF(Stat .AND. SubroutineVisited > i) THEN
    
     j = ListGetInteger(Params,'Lumped Acceleration Mode', Stat)
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
     
     clim = ListGetConstReal(Params,'Lumped Acceleration Limiter', Stat)
     IF(.NOT. Stat) clim = 100.0
     cvol = MIN(clim,cvol)
     cvol = MAX(1.0/clim,cvol)
     
     ccum = ccum * cvol
     
     WRITE(Message,'(A,T35,ES12.4)') 'Lumped Acceleration relaxation: ', ccum
     CALL Info('SteadyPhaseChange',Message)
     
     LocalRelax = LocalRelax * ccum
  END IF
  
  SurfaceMove = LocalRelax * SurfaceMove
  Surface = Surface + SurfaceMove
  
  dpos = SurfaceMove(SurfPerm(Trip_node))

  Norm = 0.0d0
  RelativeChange = 0.0d0
  DO k=1, Model % Mesh % NumberOfNodes
     IF( .NOT. BoundaryMarker(k) ) CYCLE
     Norm = Norm + Surface( SurfPerm(k) ) ** 2.0
     RelativeChange = RelativeChange + &
        SurfaceMove( SurfPerm(k) ) ** 2.0  
  END DO
  RelativeChange = SQRT(RelativeChange / Norm )
  Norm = SQRT(Norm / NoBNodes )
  
  MaxSurfaceMove = MAXVAL(ABS(SurfaceMove))
  MaxSurface = MAXVAL(ABS(Surface))
  MaxAngle = MaxSurface / Width
  
  IF ( ABS(RelativeChange) < NewtonAfterTol ) Newton = .TRUE.    
  
  
  WRITE(Message,'(A,T35,ES12.4)') 'Result Norm:', Norm
  CALL Info('SteadyPhaseChange',Message)

  WRITE(Message,'(A,T35,ES12.4)') 'Relative Change: ', RelativeChange
  CALL Info('SteadyPhaseChange',Message)

  WRITE(Message,'(A,T35,ES12.4)') 'Max Height:', MaxSurface
  CALL Info('SteadyPhaseChange',Message)

  WRITE(Message,'(A,T35,ES12.4)') 'Max Change in Height:', MaxSurfaceMove
  CALL Info('SteadyPhaseChange',Message)

  WRITE(Message,'(A,T35,ES12.4)') 'Max Angle:       ', MaxAngle
  CALL Info('SteadyPhaseChange',Message)

  WRITE(Message,'(A,T35,ES12.4)') 'Maximum temperature difference: ', MaxTempDiff
  CALL Info('SteadyPhaseChange',Message)
  
  Solver % Variable % Norm = Norm
  
  prevvolume = volume
  prevtave = tave
  prevtabs = tabs
  prevvolabs = volabs
 
  IF(IsoSurfAllocated) THEN
     DEALLOCATE(IsoSurf)
     IsoSurfAllocated = .FALSE.
  END IF
   
  IF( ListGetLogical(Params,'Internal Mesh Movement',GotIt)) THEN
     CALL BoxMoveMesh()
  END IF

  
200 CALL ListAddConstReal(Model % Simulation,'res: Triple point temperature',Trip_temp)
  CALL ListAddConstReal( Model % Simulation,'res: triple point movement',dpos)

  ! Add the coordinate y to the list of variables to save
  HelpSol => VariableGet( Solver % Mesh % Variables,'newy', ThisOnly=.TRUE.)
  IF(.NOT. ASSOCIATED(HelpSol)) THEN
     newy => Solver % Mesh % Nodes % y
     CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, &
          PSolver,'NewY',1,Newy)
  END IF

  FirstTime = .FALSE.
  
  
CONTAINS 



!-------------------------------------------------------------------------------------------
!> Subroutine creates isotherm where the temperature coincides with the melting temperature.
!-------------------------------------------------------------------------------------------
  SUBROUTINE CreateIsotherm()

    IsoSurfAllocated = .FALSE.
    xmin = HUGE(xmin)
    xmax = -HUGE(xmax)
    
100 NElems = 0
    
    DO t=1,Solver % Mesh % NumberOfBulkElements 
       
       Element => Solver % Mesh % Elements(t)
       ElementCode = Element % TYPE % ElementCode
       IF(ElementCode < 300) CYCLE
       
       k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material')
       Active = .FALSE.
       
        IF( GotLiquid ) THEN           
           IF( ListGetLogical(Model % Materials(k) % Values, 'Liquid', stat) ) &
                Active = .TRUE.
        ELSE
           IF( Element % BodyId == LiquidBody) Active = .TRUE.
        END IF
        
        IF( GotSolid ) THEN           
           IF( ListGetLogical(Model % Materials(k) % Values, 'Solid', stat) ) &
                Active = .TRUE.
        ELSE
           IF( Element % BodyId == SolidBody) Active = .TRUE.
        END IF
        
        IF(.NOT. Active) CYCLE
        
        n = Element % TYPE % NumberOfNodes
        Indexes => Element % NodeIndexes 
        IF ( ANY( TempPerm( Indexes(1:n) ) <= 0 ) )  CYCLE
        TempDiff(1:n) = Temperature(TempPerm(Indexes(1:n))) - MeltPoint
        
        IF( ALL ( TempDiff(1:n) < 0.0 ) ) CYCLE
        IF( ALL ( TempDiff(1:n) > 0.0 ) ) CYCLE       
        
        Nodes % x(1:n) = Solver % Mesh % Nodes % x(Indexes)
        Nodes % y(1:n) = Solver % Mesh % Nodes % y(Indexes)
        
        Vertex = ElementCode / 100
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
                 CALL Warn('SteadyPhaseChange','Wiggly Isotherm')
                 WRITE(Message,*) 'nodeindexes',Indexes(1:Vertex)
                 CALL Warn('SteadyPhaseChange',Message)
                 WRITE(Message,*) 'temperature',Temperature(TempPerm(Indexes))
                 CALL Warn('SteadyPhaseChange',Message)
                 PRINT *,'Nodes % x',Nodes % x(n)
                 PRINT *,'Nodes % y',Nodes % y(n)
              END IF
           END IF
        END DO
        
        IF ( n == 1 ) THEN
           IF( IsoSurfAllocated ) THEN
              IsoSurf(NElems,1) = 0.0_dp
              IsoSurf(NElems,2) = 0.0_dp
           END IF
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
     
     IF(Nelems == 0) CALL Fatal('SteadyPhaseChange','Isotherm is empty thus cannot map phase change surface') 
     
     IF(.NOT. IsoSurfAllocated) THEN
        ALLOCATE( IsoSurf(Nelems+1,2))
        IsoSurfAllocated = .TRUE.
        WRITE(Message,'(A,T35,I12)') 'Number of isotherm segments:',Nelems
        CALL Info('SteadyPhaseChange',Message)
        GOTO 100
     END IF

     ! The last one is just one extra node for safety
     IF(ListGetLogical(Params,'Save Isotherm',stat)) THEN
        CALL Info('SteadyPhaseChange','Isotherm saved to file isotherm.dat')
        OPEN (10,FILE='isotherm.dat')
        DO i=1,Nelems
           WRITE(10,*) i,IsoSurf(i,1),IsoSurf(i,2)
        END DO
        CLOSE(10)
     END IF


  END SUBROUTINE CreateIsotherm



!-------------------------------------------------------------------------------------------  
!> Internal mesh update strategy suitable for some simple geometries.
!> Assumes that the deformation is gradually decaying to the rectangle edges.
!-------------------------------------------------------------------------------------------  
  SUBROUTINE BoxMoveMesh()
    
    REAL(KIND=dp) :: x0, y0, ytop, ybot, coeff, &
         yc, dy, q, r, dx, dxmin, s1, s2
    REAL(KIND=dp), POINTER :: newy(:)
    INTEGER :: tmin
    LOGICAL :: hit
    
    coeff = 1.0_dp
    x0 = Solver % Mesh % Nodes % x(trip_node)
    y0 = Solver % Mesh % Nodes % y(trip_node)
    ytop = y0 + coeff * x0
    ybot = y0 - coeff * x0
    
    
    DO k = 1, Solver % Mesh % NumberOfNodes
       
       xx = Solver % Mesh % Nodes % x(k)
       IF(xx > x0 ) CYCLE
       
       yy = Solver % Mesh % Nodes % y(k)
       IF( yy > ytop ) CYCLE
       IF( yy < ybot ) CYCLE
       
       ! Skip nodes that are on the BC and do them last
       !--------------------------------------------------
       IF( BoundaryMarker(k) ) CYCLE
              
       ! Check if the node is fully within element
       !------------------------------------------
       hit = .FALSE.
       DO t = 1, NoPhaseElements
          Element => GetBoundaryElement(PhaseElements(t))
          
          n = GetElementNOFNodes()
          Indexes => Element % NodeIndexes
          
          ! So far only linear elements
          x1 = Solver % Mesh % Nodes % x( Indexes(1) )
          x2 = Solver % Mesh % Nodes % x( Indexes(2) )
          
          IF ( (xx - x1 ) * ( x2 - xx) >= 0.0_dp ) THEN
             hit = .TRUE. 
             EXIT
          END IF
       END DO
       
       ! If not, use the element with a closest node
       !--------------------------------------------
       IF(.NOT. hit) THEN
          tmin = 0
          dxmin = HUGE(dxmin)
          
          DO t = 1, NoPhaseElements
             Element => GetBoundaryElement(PhaseElements(t))
             
             n = GetElementNOFNodes()
             Indexes => Element % NodeIndexes
             
             DO i=1,2
                x1 = Solver % Mesh % Nodes % x( Indexes(i) )
                dx = ABS( x1 - xx )
                IF( dx < dxmin ) THEN
                   dxmin = dx
                   tmin = t
                END IF
             END DO
          END DO
          
          t = tmin
          Element => GetBoundaryElement(PhaseElements(t))
          
          n = GetElementNOFNodes()
          Indexes => Element % NodeIndexes
          
          x1 = Solver % Mesh % Nodes % x( Indexes(1) )
          x2 = Solver % Mesh % Nodes % x( Indexes(2) )
       END IF
 

       y1 = Solver % Mesh % Nodes % y( Indexes(1) )
       y2 = Solver % Mesh % Nodes % y( Indexes(2) )
       
       ! ratio at where the node is an line segment
       q = (xx - x1) / (x2 - x1)

       ! the corresponding y-coordinate
       yc = (1-q) * y1 + q * y2

       IF(yy > yc ) THEN
          r = 1.0_dp - (yy - yc) / ( ytop - yc )
       ELSE
          r = 1.0_dp - (yy - yc) / ( ybot - yc )
       END IF
       
       s1 = SurfaceMove( SurfPerm(Indexes(1)) )
       s2 = SurfaceMove( SurfPerm(Indexes(2)) )

       dy = (1-q) * s1 + q * s2
       Solver % Mesh % Nodes % y(k) = Solver % Mesh % Nodes % y(k) + r * dy
    END DO

    ! Finally do the nodes within the solver
    !-----------------------------------------
    DO k=1, Solver % Mesh % NumberOfNodes
       IF( .NOT. BoundaryMarker(k) ) CYCLE
       dy = SurfaceMove( SurfPerm(k) )
       Solver % Mesh % Nodes % y(k) =  Solver % Mesh % Nodes % y(k) + dy
    END DO
    
END SUBROUTINE BoxMoveMesh



!------------------------------------------------------------------------------
  END SUBROUTINE SteadyPhaseChange
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> This subroutine computes the heat flux resulting from solidification.
!> It may be computed when the velocity of the solidification front is known
!> a priori as is the case for various steady state pulling techniques.
!> \ingroup UDF
!-------------------------------------------------------------------------------
  FUNCTION MeltingHeat(Model, Node, t) RESULT(Flux)
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
  IF(NodeIndexes(i) /= Node) CALL Fatal('MeltingHeat','Node not found')
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

  IF(.NOT. ASSOCIATED(Parent)) CALL Fatal('MeltingHeat','Parent not found')
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



