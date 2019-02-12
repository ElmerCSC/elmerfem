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
! *  Module for computing cost function of transient measurement data.
! *
! *  Authors: Peter Råback
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 12.02.2019
! *
! *****************************************************************************/

SUBROUTINE TransientCost_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  LOGICAL :: Found
  
  Name = ListGetString( Solver % Values, 'Equation',Found)
  IF(.NOT. Found ) Name = "CostFunction"
  CALL ListAddNewString( Solver % Values,'Variable',&
      '-nooutput -global '//TRIM(Name)//'_var')
  
END SUBROUTINE TransientCost_init
  

!------------------------------------------------------------------------------
!> Calculate transient cost function using tabulated data at measurement points.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE TransientCost( Model,Solver,dt,TransientSimulation )
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
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Element_t),POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(Variable_t), POINTER :: Var
  REAL(KIND=DP), POINTER :: DataCoordinates(:,:), DataWeights(:,:), DerDataWeights(:,:)
  REAL(KIND=dP), ALLOCATABLE :: Basis(:), LocalVals(:)
  INTEGER :: NoPoints, NoDims, Point, NormType, istat, i, j, n
  INTEGER, POINTER :: NodeIndexes(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: str    
  CHARACTER(*), PARAMETER :: Caller = 'TransientCost'
  LOGICAL :: Found, Stat, Hit, Visited = .FALSE.
  REAL(KIND=dp) :: TotCost, dCost, t0, t1, time, Coords(3), LocalCoords(3), &
      SimVal, RefVal, detJ

  TYPE ControlPoint_t
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: Basis(:)
    REAL(KIND=dp) :: Cost
  END TYPE ControlPoint_t
  
  TYPE(ControlPoint_t), ALLOCATABLE, TARGET :: ControlPoints(:)
  TYPE(ControlPoint_t), POINTER :: ControlPoint
 
  
  SAVE Visited, DataCoordinates, DataWeights, DerDataWeights, NoPoints, NoDims, &
      ControlPoints, Basis, LocalVals, t0, t1, NormType, Var
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------------')
  CALL Info(Caller,'Computing cost function for transient measurement data')
 
  Mesh => GetMesh()
  Params => GetSolverParams()
  
  IF( .NOT. Visited ) THEN
    DataCoordinates => ListGetConstRealArray(Params,'Data Coordinates',Found)
    IF(.NOT. Found) CALL Fatal(Caller,'Solver expects > Data Coordinates < ')
    
    NoPoints = SIZE(DataCoordinates,1)
    CALL Info(Caller,'Number of data points: '//TRIM(I2S(NoPoints)))

    NoDims = SIZE(DataCoordinates,2)
    CALL Info(Caller,'Dimension of data points: '//TRIM(I2S(NoDims)))
        
    DataWeights => ListGetConstRealArray(Params,&
        'Data Weights',Found)
    DerDataWeights => ListGetConstRealArray(Params,&
        'Derivative Data Weights',Found)
    IF( Found ) THEN
      CALL Warn(Caller,'>Derivative Data Weights< not dealt with yet!')           
    END IF

    t0 = ListGetConstReal(Params,'Start Time',Found )
    t1 = ListGetConstReal(Params,'Stop Time',Found )
    IF(.NOT. Found ) t1 = HUGE(t1)

    NormType = ListGetInteger( Params,'Norm Type',Found )
    IF( .NOT. Found ) NormType = 2
    
    str = GetString( Params,'Target Variable',Found )
    IF(.NOT. Found ) THEN
      str = 'Temperature'
      CALL Info(Caller,'>Target Variable< defaulted to: Temperature')
    END IF
    Var => VariableGet( Mesh % Variables, str )
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal(Caller,'Variable not present: '//TRIM(str))
    END IF
    IF( Var % Dofs > 1 ) THEN
      CALL Fatal(Caller,'Currently implemented only for scalar variables: '//TRIM(str))   
    END IF

    n = Mesh % MaxElementNodes
    
    CALL Info(Caller,'Searching for elements containing data points',Level=8)
    ALLOCATE(ControlPoints(NoPoints), Basis(n), LocalVals(n), STAT=istat)
    IF( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error')         
    
    DO i=1,NoPoints
      ControlPoint => ControlPoints(i)
      ControlPoint % Element => NULL()
      ControlPoint % Basis => NULL()
      ControlPoint % Cost = 0.0_dp
    END DO
      
    DO Point =1,NoPoints
      Coords(1:NoDims) = DataCoordinates(Point,1:NoDims)
      IF(NoDims < 3 ) Coords(NoDims+1:3) = 0.0_dp

      j = ClosestElementInMesh( Mesh, Coords )
      IF( j == 0 ) CYCLE

      Element => Mesh % Elements(j)
      CALL GetElementNodes( Nodes, Element )

      Hit = PointInElement( Element, Nodes, &
          Coords, LocalCoords, GlobalEps = 1.0_dp, LocalEps=0.1_dp )	          
      
      stat = ElementInfo( Element, Nodes, LocalCoords(1), LocalCoords(2), LocalCoords(3), &
          detJ, Basis )
      ControlPoints(Point) % Element => Element

      n = Element % TYPE % NumberOfNodes    
      ALLOCATE( ControlPoints(Point) % Basis(n) )
      ControlPoints(Point) % Basis(1:n) = Basis(1:n)
    END DO

    ! Add a starting value so that SaveScalars format stays the same
    CALL ListAddConstReal( Params,'res: transient cost',0.0_dp)

    Visited = .TRUE.
  END IF
  
  time = GetTime()

  IF( time < t0 ) THEN
    CALL Info(Caller,'Integration not started yet, exiting',Level=6)
    RETURN
  END IF

  IF( time > t1 ) THEN
    CALL Info(Caller,'Integration finished already, exiting',Level=6)
    RETURN
  END IF

  
  DO Point = 1, NoPoints
    ControlPoint => ControlPoints(Point)
    
    ! This is True only for one partition in parallel runs
    Element => ControlPoint % Element

    IF(.NOT. ASSOCIATED( Element ) ) CYCLE
    Model % CurrentElement => Element
    
    ! These are expected to be tables such that they depend on time
    str = 'Data '//TRIM(I2S(Point))
    RefVal = ListGetCReal( Params, str, Found )
    IF( .NOT. Found ) THEN
      CALL Fatal( Caller,'Measurement data not found: '//TRIM(str))
    END IF
        
    NodeIndexes => Element % NodeIndexes
    n = Element % Type % NumberOfNodes 
               
    IF( ANY(Var % Perm(NodeIndexes(1:n)) == 0)) THEN            
      CALL Fatal(Caller,'Cannot evaluate field with zero permutation')
    END IF

    LocalVals(1:n) = Var % Values(Var % Perm(NodeIndexes(1:n)))
    SimVal = SUM( ControlPoint % Basis(1:n) * LocalVals(1:n) )

    IF( NormType == 1 ) THEN
      dCost = dt * ABS(RefVal-SimVal)
    ELSE IF( NormType == 2 ) THEN
      dCost = dt * (RefVal-SimVal)**2
    ELSE
      CALL Fatal(Caller,'Invalid norm type: '//TRIM(I2S(NormType)))
    END IF
    
    IF( ASSOCIATED( DataWeights ) ) THEN
      dCost = dCost * DataWeights(Point,1)
    END IF

    ControlPoint % Cost = ControlPoint % Cost + dCost
  END DO

  TotCost = SUM( ControlPoints(1:NoPoints) % Cost )

  ! In parallel runs each partition only deals with the local hit and only the
  ! total cost is summed up. 
  TotCost = ParallelReduction( TotCost ) 

  Solver % Variable % Values = TotCost
  Solver % Variable % Norm = TotCost
  
  CALL ListAddConstReal( Params,'res: transient cost',TotCost)


END SUBROUTINE TransientCost
  
