!/*****************************************************************************
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
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.11.2005
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
!> Determines a timestep based on the maximum local Courant number. 
!> \ingroup UDF
!------------------------------------------------------------------------------
FUNCTION LevelSetTimestep( Model ) RESULT( dt )
  USE Types
  USE Lists
  USE Integration
  USE ElementDescription

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt

  TYPE(Solver_t), POINTER :: Solver
  TYPE(Variable_t), POINTER ::  TimeVariable, SurfSol
  REAL(KIND=dp), POINTER :: Surface(:), Surf(:)
  INTEGER, POINTER :: SurfPerm(:), NodeIndexes(:)
  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:)
 
  INTEGER :: i,j,k,n,t,elem,N_Integ,body_id,mat_id, dim, TimeIntervals
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), NodalVelo(:,:)
  REAL(KIND=dp) :: Val,Grad(3),Velo(3),NormalVelo,AbsVelo,GradAbs,detJ,&
      MaxNormVelo,MaxAbsVelo
  REAL(KIND=dp) :: s,u,v,w,dt0,prevdt,cumtime,dsmax
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
  LOGICAL :: stat, GotIt, AllocationsDone = .FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: LevelSetVariableName

  SAVE AllocationsDone, Basis, dBasisdx, NodalVelo, Surf, ElementNodes, prevdt


  TimestepSizes => ListGetConstRealArray( CurrentModel % Simulation,&
      'Timestep Sizes', GotIt )
  TimeIntervals = SIZE(TimestepSizes)
 
  IF(TimeIntervals > 1) THEN
    CALL Warn('LevelSetTimestep','Implemented only for one Time Interval')
  END IF

  dt0 = TimestepSizes(1,1)

  TimeVariable => VariableGet(Model % Variables, 'Time')
  cumtime = TimeVariable % Values(1)

  dim = CoordinateSystemDimension()

  ! The variable that should be reinitialized
  LevelSetVariableName = 'Surface'
  SurfSol => VariableGet( Model % Variables, TRIM(LevelSetVariableName) )
  IF(ASSOCIATED(SurfSol)) THEN
    Surface  => SurfSol % Values
    SurfPerm => SurfSol % Perm
  ELSE
    CALL Warn('LevelSetTimeStep','SurfSol does not exist: '//TRIM(LevelSetVariableName))
    RETURN
  END IF

  Solver => SurfSol % Solver
  DsMax = ListGetConstReal(Model % Simulation,'LevelSet Courant Number',GotIt)
  IF(.NOT. GotIt) DsMax = 1.0d0

  IF(.NOT. AllocationsDone) THEN
    N = Solver % Mesh % MaxElementNodes
    ALLOCATE( Basis(n), dBasisdx(n,3), NodalVelo(3,n), Surf(n), &
        ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    prevdt = dt0
    AllocationsDone = .TRUE. 
  END IF

  MaxNormVelo = 0.0d0
  MaxAbsVelo = 0.0d0

  DO elem=1,Solver % Mesh % NumberOfBulkElements
    
    CurrentElement => Solver % Mesh % Elements(elem)
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
    
    IF(ANY(SurfPerm(NodeIndexes) == 0)) CYCLE
    
    Surf(1:n) = Surface( SurfPerm(NodeIndexes) )
    IF(ALL(Surf(1:n) < 0.0d0) .OR. ALL(Surf(1:n) > 0.0d0) ) CYCLE

    ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
    ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
    ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
    
    Model % CurrentElement => CurrentElement
    body_id = CurrentElement % Bodyid    
    mat_id = ListGetInteger( Model % Bodies( body_id ) % Values, 'Material' )
    Material => Model % Materials(mat_id) % Values

!------------------------------------------------------------------------------
!         Computed or given velocity field
!------------------------------------------------------------------------------
       
    NodalVelo(1,1:n) = ListGetReal( Material,'Levelset Velocity 1',n,NodeIndexes,GotIt)
    NodalVelo(2,1:n) = ListGetReal( Material,'Levelset Velocity 2',n,NodeIndexes,GotIt)
    IF(dim == 3) NodalVelo(3,1:n) = ListGetReal( Material,'Levelset Velocity 3',n,NodeIndexes,GotIt)

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( CurrentElement )
    
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Maximum at any integration point
!------------------------------------------------------------------------------
    
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
         
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

      stat = ElementInfo( CurrentElement,ElementNodes,u,v,w,detJ, Basis,dBasisdx)
      
      DO i=1,dim
        Grad(i) = SUM( dBasisdx(1:n,i) * Surf(1:n) )
        Velo(i) = SUM( Basis(1:n) * NodalVelo(i,1:n) )
      END DO
      
      GradAbs = SQRT( SUM( Grad(1:dim) * Grad(1:dim) ) )         
      NormalVelo = SUM( Grad(1:dim) * Velo(1:dim) ) / GradAbs
      NormalVelo = NormalVelo / SQRT(detJ)
      MaxNormVelo = MAX(MaxNormVelo, ABS(NormalVelo))
      
      AbsVelo = SQRT(SUM (Velo(1:dim) * Velo(1:dim)) )
      AbsVelo = AbsVelo / SQRT(detJ)
      MaxAbsVelo = MAX(MaxAbsVelo, ABS(AbsVelo))
    END DO
  END DO

  dt = dt0
  IF( ListGetLogical(Model % Simulation,'LevelSet Timestep Directional',GotIt) ) THEN
    IF( MaxNormVelo * dt0 > dsMax) dt = dsMax / MaxNormVelo
  ELSE
    IF( MaxAbsVelo * dt0 > dsMax) dt = dsMax / MaxAbsVelo
  END IF
  
  IF( dt < dt0) THEN
    IF(dt > prevdt) dt = 0.5d0 * (dt + prevdt)
  END IF
  prevdt = dt
  
  WRITE(Message,'(a,ES12.3)') 'Levelset timestep',dt
  CALL Info( 'LevelSetTimestep',Message, Level=4 )
  CALL ListAddConstReal(Model % Simulation,'res: Levelset timestep',dt)

END FUNCTION LevelSetTimestep

