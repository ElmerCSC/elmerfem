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
!
!/******************************************************************************
! *
! *  Generic utilities related to components
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 28 Sep 1998
! *
! *****************************************************************************/

!> Generic utilities related to components
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE LumpingUtils

   USE ElementUtils
   USE ModelDescription
   IMPLICIT NONE

 CONTAINS


!------------------------------------------------------------------------------
!> Compute reduction operators for a given component with given nodal vector field.
!> Force is simple sum of nodal forces
!> Moment is moment of nodal forces about a given center point
!> Torque is moment of nodal forces about a given rotational axis
!> If the given field is a elemental (DG) field it may be reduced by 
!> optional SetPerm reordering for minimal discontinuous set. 
!------------------------------------------------------------------------------
   SUBROUTINE ComponentNodalForceReduction(Model, Mesh, CompParams, NF, &
       Force, Moment, Torque, SetPerm ) 
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: CompParams
     TYPE(Variable_t), POINTER :: NF
     REAL(KIND=dp), OPTIONAL :: Moment(3), Force(3), Torque
     INTEGER, POINTER, OPTIONAL :: SetPerm(:)
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     TYPE(Nodes_t) :: Nodes
     LOGICAL, ALLOCATABLE :: VisitedNode(:)
     REAL(KIND=dp) :: Origin(3), Axis(3), P(3), F(3), v1(3), v2(3), &
         RotorRadius, rad, minrad, maxrad, eps
     REAL(KIND=dp), POINTER :: Pwrk(:,:)
     INTEGER :: t, i, j, k, n, dofs, globalnode, AirBody
     LOGICAL :: ElementalVar, Found, NeedLocation
     INTEGER, POINTER :: MasterEntities(:),NodeIndexes(:),DofIndexes(:)
     LOGICAL :: VisitNodeOnlyOnce     
     INTEGER :: FirstElem, LastElem
     LOGICAL :: BcMode, BulkMode, RotorMode, isParallel
     CHARACTER(*), PARAMETER :: Caller = 'ComponentNodalForceReduction'
    
     CALL Info(Caller,'Performing reduction for component: '&
         //TRIM(ListGetString(CompParams,'Name')),Level=10)

     IF(.NOT. (PRESENT(Torque) .OR. PRESENT(Moment) .OR. PRESENT(Force) ) ) THEN
       CALL Warn(Caller,'Nothing to compute!')
       RETURN
     END IF

     IF( PRESENT(Torque)) Torque = 0.0_dp
     IF( PRESENT(Moment)) Moment = 0.0_dp
     IF( PRESENT(Force)) Force = 0.0_dp

     isParallel = CurrentModel % Solver % Parallel


     eps = 1.0e-6
     BcMode = .FALSE.
     BulkMode = .FALSE.
     RotorMode = ListGetLogical( CompParams,'Rotor Mode',Found )
     IF( RotorMode ) THEN
       RotorRadius = ListGetConstReal( CurrentModel % Simulation,'Rotor Radius',Found )
       IF(.NOT. Found ) THEN
         CALL Fatal(Caller,'"Rotor Mode" requires "Rotor Radius"')
       END IF
     ELSE
       MasterEntities => ListGetIntegerArray( CompParams,'Master Bodies',BulkMode )     
       IF( .NOT. BulkMode ) THEN
         MasterEntities => ListGetIntegerArray( CompParams,'Master Boundaries', BCMode) 
       END IF                    
       IF(.NOT. (BulkMode .OR. BCMode ) ) THEN
         CALL Warn(Caller,'> Master Bodies < or > Master Boundaries < not given')
         RETURN
       END IF
     END IF
       
     NeedLocation = PRESENT( Moment ) .OR. PRESENT( Torque )

     ! User may specific origin and axis for torque computation
     ! By default (0,0,0) is the origin, and (0,0,1) the axis. 
     Pwrk => ListGetConstRealArray( CompParams,'Torque Origin',Found )
     IF( Found ) THEN
       IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
         CALL Fatal(Caller,'Size of > Torque Origin < should be 3!')
       END IF
       Origin = Pwrk(1:3,1)
     ELSE
       Origin = 0.0_dp
     END IF
     Pwrk => ListGetConstRealArray( CompParams,'Torque Axis',Found )
     IF( Found ) THEN
       IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
         CALL Fatal(Caller,'Size of > Torque Axis < should be 3!')
       END IF
       Axis = Pwrk(1:3,1)
       ! Normalize axis is it should just be used for the direction
       Axis = Axis / SQRT( SUM( Axis*Axis ) )
     ELSE
       Axis = 0.0_dp    
       Axis(3) = 1.0_dp  
     END IF

     ElementalVar = ( NF % TYPE == Variable_on_nodes_on_elements )
     IF( PRESENT( SetPerm ) .AND. .NOT. ElementalVar ) THEN
       CALL Fatal(Caller,'SetPerm is usable only for elemental fields')
     END IF

     dofs = NF % Dofs
     IF( dofs == 2 ) F(3) = 0.0_dp

     ! For nodal field compute only once each node
     ! For DG field each node is visited only once by construction
     VisitNodeOnlyOnce = .NOT. ElementalVar .OR. PRESENT(SetPerm)
     IF( VisitNodeOnlyOnce ) THEN
       IF( PRESENT( SetPerm ) ) THEN
         n = MAXVAL( SetPerm ) 
       ELSE
         n = Mesh % NumberOfNodes
       END IF
       ALLOCATE(VisitedNode( n ) )
       VisitedNode = .FALSE.
     END IF

     IF( BcMode ) THEN
       FirstElem = Mesh % NumberOfBulkElements + 1
       LastElem = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     ELSE       
       FirstElem = 1
       LastElem = Mesh % NumberOfBulkElements
     END IF

     ! This is a special reduction that only applies to rotors that are surrounded by airgap.
     ! Very special operator for electrical machines that removes the bookkeeping of the bodies
     ! that constitute the rotor. 
     AirBody = 0
     IF(RotorMode ) THEN       
       AirBody = ListGetInteger( CompParams,'Air Body',Found ) 
       IF(AirBody == 0) THEN
         DO t=FirstElem,LastElem
           Element => Mesh % Elements(t)
           n = Element % TYPE % NumberOfNodes
           CALL CopyElementNodesFromMesh( Nodes, Mesh, n, Element % NodeIndexes)         
           DO i=1,n
             rad = SQRT(Nodes % x(i)**2 + Nodes % y(i)**2)
             IF(i==1) THEN
               minrad = rad
               maxrad = rad
             ELSE
               minrad = MIN(minrad,rad)
               maxrad = MAX(maxrad,rad)
             END IF
           END DO

           ! The body is defined by an element that is at and inside the rotor radius. 
           IF(ABS(maxrad-RotorRadius) < eps .AND. minrad < RotorRadius*(1-eps) ) THEN
             AirBody = Element % BodyId
             EXIT
           END IF
         END DO
         AirBody = ParallelReduction(AirBody,2)         
         CALL Info(Caller,'Airgap inner body determined to be: '//I2S(AirBody),Level=12)           
         IF(AirBody==0) THEN
           CALL Fatal(Caller,'Could not define airgap inner body!')
         ELSE
           CALL ListAddInteger(CompParams,'Air Body',AirBody)
         END IF
       END IF
     END IF

    
     DO t=FirstElem,LastElem
       Element => Mesh % Elements(t)

       IF( BcMode ) THEN
         IF( ALL( MasterEntities /= Element % BoundaryInfo % Constraint ) ) CYCLE
       ELSE IF( BulkMode ) THEN
         IF( ALL( MasterEntities /= Element % BodyId ) ) CYCLE
       ELSE IF( RotorMode ) THEN
         IF( Element % BodyId == AirBody ) CYCLE
         CALL CopyElementNodesFromMesh( Nodes, Mesh, n, Element % NodeIndexes)         
         rad = SQRT((SUM(Nodes % x(1:n))/n)**2 + (SUM(Nodes % y(1:n))/n)**2)
         IF(rad > RotorRadius ) CYCLE         
       END IF

       n = Element % TYPE % NumberOfNodes
       NodeIndexes => Element % NodeIndexes 
       IF( ElementalVar ) THEN
         DofIndexes => Element % DGIndexes
       ELSE
         DofIndexes => NodeIndexes
       END IF

       DO i=1,n
         j = DofIndexes(i)        
         k = NF % Perm(j)
         IF( k == 0 ) CYCLE

         IF( VisitNodeOnlyOnce ) THEN
           IF( PRESENT( SetPerm ) ) j = SetPerm(j)
           IF( VisitedNode(j) ) CYCLE
           VisitedNode(j) = .TRUE.
         END IF

         globalnode = NodeIndexes(i)

         ! Only compute the parallel reduction once
         IF( isParallel ) THEN
           IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE

! This is (probably) not correct, the "nodal forces"-array is partial and should be summed --> comment out.
!          IF( VisitNodeOnlyOnce ) THEN           
!            IF( Mesh % ParallelInfo % NeighbourList(globalnode) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
!          END IF
         END IF
           
         F(1) = NF % Values( dofs*(k-1) + 1)
         F(2) = NF % Values( dofs*(k-1) + 2)
         IF( dofs == 3 ) THEN 
           F(3) = NF % Values( dofs*(k-1) + 3)
         END IF

         IF( PRESENT( Force ) ) THEN
           ! calculate simple sum
           Force = Force + F
         END IF

         IF( NeedLocation ) THEN
           P(1) = Mesh % Nodes % x(globalnode)
           P(2) = Mesh % Nodes % y(globalnode)
           P(3) = Mesh % Nodes % z(globalnode)

           v1 = P - Origin

           ! Calculate moment 
           IF( PRESENT( Moment ) ) THEN
             Moment = Moment + CrossProduct(v1,F)
           END IF

           ! Calculate torque around an axis
           IF( PRESENT( Torque ) ) THEN
             v1 = v1 - SUM(Axis*v1)*Axis
             v2 = CrossProduct(v1,F)
             Torque = Torque + SUM(Axis*v2)        
           END IF
         END IF

       END DO
     END DO

     IF( isParallel ) THEN
       IF( PRESENT( Force ) ) THEN
         DO i=1,3
           Force(i) = ParallelReduction(Force(i))
         END DO
       END IF
       
       IF( PRESENT( Moment ) ) THEN
         DO i=1,3
           Moment(i) = ParallelReduction(Moment(i))
         END DO
       END IF
       
       IF( PRESENT( Torque ) ) THEN
         Torque = ParallelReduction(Torque)
       END IF
     END IF
       
!------------------------------------------------------------------------------
   END SUBROUTINE ComponentNodalForceReduction
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Perform reduction from distributed data to components.
!> These are similar operations as the stastistical operations in SaveScalars.
!------------------------------------------------------------------------------
   FUNCTION ComponentNodalReduction(Model, Mesh, CompParams, Var, OperName ) RESULT ( OperX )
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: CompParams
     TYPE(Variable_t), POINTER :: Var
     REAL(KIND=dp) :: OperX
     CHARACTER(LEN=*) :: OperName
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     LOGICAL, ALLOCATABLE :: VisitedNode(:)
     INTEGER :: t, i, j, k, l, n, NoDofs, globalnode, sumi
     REAL(KIND=dp) :: X, Minimum, Maximum, AbsMinimum, AbsMaximum, SumX, SumXX, SumAbsX
     LOGICAL :: ElementalVar, Found
     INTEGER, POINTER :: MasterEntities(:),NodeIndexes(:),DofIndexes(:)
     LOGICAL :: VisitNodeOnlyOnce, Initialized
     INTEGER :: FirstElem, LastElem
     LOGICAL :: BcMode 


     CALL Info('ComponentNodalReduction','Performing reduction for component: '&
         //TRIM(ListGetString(CompParams,'Name')),Level=10)

     OperX = 0.0_dp

     BcMode = .FALSE.
     MasterEntities => ListGetIntegerArray( CompParams,'Master Bodies',Found ) 
     IF( .NOT. Found ) THEN
       MasterEntities => ListGetIntegerArray( CompParams,'Master Boundaries',Found ) 
       BcMode = .TRUE.
     END IF

     IF(.NOT. Found ) THEN
       CALL Warn('ComponentNodalReduction',&
           '> Master Bodies < or > Master Boundaries < not given')
       RETURN
     END IF

     IF( BcMode ) THEN
       FirstElem = Mesh % NumberOfBulkElements + 1
       LastElem = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     ELSE
       FirstElem = 1 
       LastElem = Mesh % NumberOfBulkElements
     END IF

     ElementalVar = ( Var % TYPE == Variable_on_nodes_on_elements )
     NoDofs = Var % Dofs
     Initialized = .FALSE.

     ! For nodal field compute only once each node
     ! For DG field each node is visited only once by construction
     VisitNodeOnlyOnce = .NOT. ElementalVar 
     IF( VisitNodeOnlyOnce ) THEN
       n = Mesh % NumberOfNodes
       ALLOCATE(VisitedNode( n ) )
       VisitedNode = .FALSE.
     END IF

     sumi = 0
     sumx = 0.0_dp
     sumxx = 0.0_dp
     sumabsx = 0.0_dp
     Maximum = 0.0_dp
     Minimum = 0.0_dp
     AbsMaximum = 0.0_dp
     AbsMinimum = 0.0_dp

     DO t=FirstElem,LastElem
       Element => Mesh % Elements(t)

       IF( BcMode ) THEN
         IF( ALL( MasterEntities /= Element % BoundaryInfo % Constraint ) ) CYCLE
       ELSE
         IF( ALL( MasterEntities /= Element % BodyId ) ) CYCLE
       END IF

       n = Element % TYPE % NumberOfNodes
       NodeIndexes => Element % NodeIndexes 
       IF( ElementalVar ) THEN
         DofIndexes => Element % DGIndexes
       ELSE
         DofIndexes => NodeIndexes
       END IF

       DO i=1,n
         j = DofIndexes(i)        
         IF( ASSOCIATED( Var % Perm ) ) THEN
           k = Var % Perm(j)
           IF( k == 0 ) CYCLE
         ELSE
           k = j
         END IF

         IF( VisitNodeOnlyOnce ) THEN
           IF( VisitedNode(j) ) CYCLE
           VisitedNode(j) = .TRUE.
         END IF

         globalnode = NodeIndexes(i)

         ! Only compute the parallel reduction once
         IF( ParEnv % PEs > 1 ) THEN
           IF( Mesh % ParallelInfo % NeighbourList(globalnode) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
         END IF

         IF(NoDofs <= 1) THEN
           x = Var % Values(k)
         ELSE
          x = 0.0d0
          DO l=1,NoDofs
            x = x + Var % Values(NoDofs*(k-1)+l) ** 2
          END DO
          x = SQRT(x)
        END IF

        IF(.NOT. Initialized) THEN
          Initialized = .TRUE.
          Maximum = x
          Minimum = x
          AbsMaximum = x
          AbsMinimum = x
        END IF

        sumi = sumi + 1
        sumx = sumx + x
        sumxx = sumxx + x*x
        sumabsx = sumabsx + ABS( x )
        Maximum = MAX(x,Maximum)
        Minimum = MIN(x,Minimum)
        IF(ABS(x) > ABS(AbsMaximum) ) AbsMaximum = x
        IF(ABS(x) < ABS(AbsMinimum) ) AbsMinimum = x
      END DO
    END DO
    
    
    sumi = ParallelReduction(sumi) 
    IF( sumi == 0 ) THEN
      CALL Warn('ComponentNodalReduction','No active nodes to reduced!')
      RETURN
    END IF


    SELECT CASE(OperName)
      
    CASE ('sum')
      sumx = ParallelReduction(sumx)
      operx = sumx

    CASE ('sum abs')
      sumx = ParallelReduction(sumabsx)
      operx = sumabsx
      
    CASE ('min')
      minimum = ParallelReduction(minimum,1)
      operx = Minimum

    CASE ('max')
      maximum = ParallelReduction(maximum,2)
      operx = Maximum
      
    CASE ('min abs')
      Absminimum = ParallelReduction(AbsMinimum,1)          
      operx = AbsMinimum
      
    CASE ('max abs')
      Absmaximum = ParallelReduction(AbsMaximum,2)     
      operx = AbsMaximum

    CASE ('range')
      minimum = ParallelReduction(minimum,1)     
      maximum = ParallelReduction(maximum,2)
      operx = Maximum - Minimum
      
    CASE ('mean')
      sumx = ParallelReduction(sumx)
      operx = sumx / sumi 
      
    CASE ('mean abs')
      sumx = ParallelReduction(sumabsx)
      operx = sumabsx / sumi 

    CASE ('variance')
      sumx = ParallelReduction(sumx)
      sumxx = ParallelReduction(sumxx)
      Operx = SQRT( sumxx/sumi-(sumx*sumx)/(sumi*sumi) )

    CASE DEFAULT 
      CALL Warn('ComponentNodalReduction','Unknown statistical operator!')

    END SELECT
      
    CALL Info('ComponentNodalReduction','Reduction operator finished',Level=12)

!------------------------------------------------------------------------------
  END FUNCTION ComponentNodalReduction
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Perform reduction from distributed data to components.
!> These are similar operations as the stastistical operations in SaveScalars.
!------------------------------------------------------------------------------
   FUNCTION ComponentIntegralReduction(Model, Mesh, CompParams, Var, &
       OperName, CoeffName, GotCoeff ) RESULT ( OperX )
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model 
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: CompParams 
     TYPE(Variable_t), POINTER :: Var 
     CHARACTER(LEN=*) :: OperName, CoeffName
     LOGICAL :: GotCoeff
     REAL(KIND=dp) :: OperX
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: t, i, j, k, n, NoDofs
     INTEGER, POINTER :: NodeIndexes(:), DofIndexes(:)
     REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,x,Grad(3)
     REAL(KIND=dp) :: func, CoeffAtIp, integral, vol
     LOGICAL :: ElementalVar, Found, Stat
     INTEGER :: PermIndexes(Model % MaxElementNodes)
     REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
     REAL(KIND=dp) :: Coeff(Model % MaxElementNodes)
     TYPE(GaussIntegrationPoints_t) :: IntegStuff
     INTEGER, POINTER :: MasterEntities(:)
     TYPE(Nodes_t), SAVE :: ElementNodes
     LOGICAL, SAVE :: AllocationsDone = .FALSE.
     INTEGER :: FirstElem, LastElem
     LOGICAL :: BcMode 


     CALL Info('ComponentIntegralReduction','Performing reduction for component: '&
         //TRIM(ListGetString(CompParams,'Name')),Level=10)

     OperX = 0.0_dp
     vol = 0.0_dp
     integral = 0.0_dp
     
     BcMode = .FALSE.
     MasterEntities => ListGetIntegerArray( CompParams,'Master Bodies',Found ) 
     IF( .NOT. Found ) THEN
       MasterEntities => ListGetIntegerArray( CompParams,'Master Boundaries',Found ) 
       BcMode = .TRUE.
     END IF

     IF(.NOT. Found ) THEN
       CALL Warn('ComponentIntegralReduction',&
           '> Master Bodies < or > Master Boundaries < not given')
       RETURN
     END IF

     IF( BcMode ) THEN
       FirstElem = Mesh % NumberOfBulkElements + 1
       LastElem = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     ELSE
       FirstElem = 1 
       LastElem = Mesh % NumberOfBulkElements
     END IF

     IF(.NOT. AllocationsDone ) THEN
       n = Model % MaxElementNodes
       ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )      
       AllocationsDone = .TRUE.
     END IF

     ElementalVar = ( Var % TYPE == Variable_on_nodes_on_elements )
     NoDofs = Var % Dofs
     CoeffAtIp = 1.0_dp

     DO t=FirstElem,LastElem
       Element => Mesh % Elements(t)
       IF( BcMode ) THEN
         IF( ALL( MasterEntities /= Element % BoundaryInfo % Constraint ) ) CYCLE
       ELSE
         IF( ALL( MasterEntities /= Element % BodyId ) ) CYCLE
       END IF

       n = Element % TYPE % NumberOfNodes
       NodeIndexes => Element % NodeIndexes 
       IF( ElementalVar ) THEN
         DofIndexes => Element % DGIndexes
       ELSE
         DofIndexes => NodeIndexes
       END IF

       IF( ASSOCIATED( Var % Perm ) ) THEN
         PermIndexes = Var % Perm( DofIndexes )
       ELSE
         PermIndexes = DofIndexes
       END IF

       IF ( ANY(PermIndexes == 0 ) ) CYCLE      

       ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
       ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
       ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
       
       IF(GotCoeff) THEN
         k = ListGetInteger( Model % Bodies( Element % BodyId ) % Values, &
             'Material', Found )
         Coeff(1:n) = ListGetReal( Model % Materials(k) % Values, &
             CoeffName, n, NodeIndexes(1:n) )
       END IF
       
!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
       
       DO i=1,IntegStuff % n
         U = IntegStuff % u(i)
         V = IntegStuff % v(i)
         W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
             Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
         s = SqrtElementMetric * IntegStuff % s(i)
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           x = 2 * PI * SUM( ElementNodes % x(1:n)*Basis(1:n) )
         END IF
         
         IF( GotCoeff ) THEN
           CoeffAtIp = SUM( Coeff(1:n) * Basis(1:n) )
         END IF
         vol = vol + S

         SELECT CASE(OperName)
           
         CASE ('volume')
           integral = integral + coeffAtIp * S
           
         CASE ('int','int mean')
           func = SUM( Var % Values(PermIndexes(1:n)) * Basis(1:n) )
           integral = integral + S * coeffAtIp * func 
           
         CASE ('int abs','int abs mean')
           func = ABS( SUM( Var % Values(PermIndexes(1:n)) * Basis(1:n) ) )
           integral = integral + S * coeffAtIp * func
           
         CASE ('diffusive energy')
           DO j = 1, 3
             Grad(j) = SUM( dBasisdx(1:n,j) *  Var % Values(PermIndexes(1:n) ) )
           END DO           
           integral = integral + s * CoeffAtIp * SUM( Grad * Grad )
           
         CASE ('convective energy')
           func = SUM( Var % Values(PermIndexes(1:n)) * Basis(1:n) )
           
           IF(NoDofs == 1) THEN
             func = SUM( Var % Values(PermIndexes(1:n)) * Basis(1:n) )
             integral = integral + s * coeffAtIp * func**2
           ELSE
             func = 0.0d0
             DO j=1,NoDofs
              func = SUM( Var % Values(NoDofs*(PermIndexes(1:n)-1)+j) * Basis(1:n) )
              integral = integral + s * coeffAtIp * func**2
            END DO
          END IF
          
        CASE ('potential energy')          
          func = SUM( Var % Values(PermIndexes(1:n)) * Basis(1:n) )
          integral = integral + s * coeffAtIp * func
          
        CASE DEFAULT 
          CALL Warn('ComponentIntegralReduction','Unknown operator')

        END SELECT

      END DO

    END DO

    integral = ParallelReduction( integral ) 

    SELECT CASE(OperName)
      
    CASE ('volume')        
      operx = integral
      
    CASE ('int')
      operx = integral
      
    CASE ('int abs')
      operx = integral
      
    CASE ('int mean')
      vol = ParallelReduction( vol ) 
      operx = integral / vol        
      
    CASE ('int abs mean')
      vol = ParallelReduction( vol ) 
      operx = integral / vol        
      
    CASE ('diffusive energy')
      operx = 0.5d0 * integral
      
    CASE ('convective energy')
      operx = 0.5d0 * integral
      
    CASE ('potential energy')
      operx = integral
      
    END SELECT
      
    CALL Info('ComponentIntegralReduction','Reduction operator finished',Level=12)

!------------------------------------------------------------------------------
  END FUNCTION ComponentIntegralReduction
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Each solver may include a list of dependent components that are updated
!> after the solver (or the nonlinear iteration related to it) has been executed.
!------------------------------------------------------------------------------
  SUBROUTINE UpdateDependentComponents( ComponentList )
    INTEGER, POINTER :: ComponentList(:)

    INTEGER :: i,j,NoVar
    CHARACTER(:), ALLOCATABLE :: OperName, VarName, CoeffName, TmpOper
    LOGICAL :: GotVar, GotOper, GotCoeff, VectorResult
    TYPE(ValueList_t), POINTER :: CompParams
    REAL(KIND=dp) :: ScalarVal, VectorVal(3), Power, Voltage
    TYPE(Variable_t), POINTER :: Var

    CALL Info('UpdateDepedentComponents','Updating Components to reflect new solution',Level=6)

    DO j=1,CurrentModel % NumberOfComponents

      IF( ALL( ComponentList /= j ) ) CYCLE

      CALL Info('UpdateDepedentComponents','Updating component: '//I2S(j))
      CompParams => CurrentModel % Components(j) % Values


      NoVar = 0
      DO WHILE( .TRUE. )
        NoVar = NoVar + 1
        OperName = ListGetString( CompParams,'Operator '//I2S(NoVar), GotOper)
        VarName = ListGetString( CompParams,'Variable '//I2S(NoVar), GotVar)
        CoeffName = ListGetString( CompParams,'Coefficient '//I2S(NoVar), GotCoeff)
        
        IF(.NOT. GotVar .AND. GotOper .AND. OperName == 'electric resistance') THEN
          VarName = 'Potential'
          GotVar = .TRUE.
          CALL Info('UpdateDependentComponents',&
              'Defaulting field to > Potential < for operator: '//TRIM(OperName),Level=8)
        END IF
        
        IF(.NOT. (GotVar .AND. GotOper ) ) EXIT

        Var => VariableGet( CurrentModel % Mesh % Variables, VarName ) 
        IF( .NOT. ASSOCIATED( Var ) ) THEN
          CALL Info('UpdateDependentComponents','Variable not available: '//TRIM(VarName))
          CYCLE
        END IF
        VectorResult = .FALSE.

        SELECT CASE( OperName ) 

        CASE('electric resistance')
          IF(.NOT. GotCoeff ) THEN
            CoeffName = 'electric conductivity'
            GotCoeff = .TRUE.
          END IF
          TmpOper = 'diffusive energy'
          Power = ComponentIntegralReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              TmpOper, CoeffName, GotCoeff )
          TmpOper = 'range'
          Voltage = ComponentNodalReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              TmpOper )
          ScalarVal = Voltage**2 / Power 
          CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName),ScalarVal )
 
        CASE ('sum','sum abs','min','max','min abs','max abs','range','mean','mean abs','variance')
          ScalarVal = ComponentNodalReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              OperName )
          CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName)//' '//TRIM(VarName),ScalarVal )

        CASE ('volume','int','int abs','int mean','int abs mean','diffusive energy',&
            'convective energy','potential energy')
          ScalarVal = ComponentIntegralReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              OperName, CoeffName, GotCoeff )
          CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName)//' '//TRIM(VarName),ScalarVal )

        CASE('force')
          CALL ComponentNodalForceReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              Force = VectorVal )
          VectorResult = .TRUE.

        CASE('moment')
          CALL ComponentNodalForceReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              Moment = VectorVal )
          VectorResult = .TRUE.

        CASE('torque')          
          CALL ComponentNodalForceReduction(CurrentModel, CurrentModel % Mesh, CompParams, Var, &
              Torque = ScalarVal )

        CASE DEFAULT
          CALL Fatal('UpdateDependentComponents','Uknown operator: '//TRIM(OperName))

        END SELECT

        IF( VectorResult ) THEN
          DO i=1,3
            WRITE( Message,'(A,ES15.6)') TRIM(OperName)//': '//TRIM(VarName)//' '&
                //I2S(i)//': ',ScalarVal
            CALL Info('UpdateDependentComponents',Message,Level=5)
            CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName)//': '&
                //TRIM(VarName)//' '//I2S(i),VectorVal(i) )                        
          END DO
        ELSE          
          WRITE( Message,'(A,ES15.6)') &
              'comp '//I2S(j)//': '//TRIM(OperName)//': '//TRIM(VarName)//': ',ScalarVal
          CALL Info('UpdateDependentComponents',Message,Level=5)
          CALL ListAddConstReal( CurrentModel % Simulation, &
              'res: comp '//I2S(j)//': '//TRIM(OperName)//' '//TRIM(VarName),ScalarVal )           
        END IF

      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateDependentComponents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Given a vector field compute line integral of Stokes theorem using
!> some geometric heuristics. 
!------------------------------------------------------------------------------
  FUNCTION ComponentStokesTheorem(Model, Mesh, Vlist, AVar, Surf ) RESULT ( FL ) 
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist
    TYPE(Variable_t), POINTER :: aVar
    LOGICAL :: Surf
    REAL(KIND=dp) :: FL
         
    TYPE(Matrix_t), POINTER :: NodeGraph
    REAL(KIND=dp), POINTER :: HelperArray(:,:)
    REAL(KIND=dp) :: Center(3), Coord(3), Coord0(3), Coord1(3), Coord2(3), &
        Normal(3), Tangent1(3), Tangent2(3), x1, y1, x2, y2, phi, phi1, phi2, &
        phisum, g, ssum, h, hmin, hmax, htol, hgoal, hrel, r, rmin, rmax, rtol
    COMPLEX :: gradv(3), Circ
    REAL(KIND=dp) :: EdgeVector(3), ds, r2, r2min, dsmax, dphi, dphimax
    INTEGER :: nsteps, sgn, t, i, j, k, i1, i2, j1, j2, l, &
        lp, n0, r2ind, imax, kmax, n
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, POINTER :: TargetBodies(:)
    LOGICAL :: Found, SaveLoop, Debug, GotHtol, GotRtol 
    TYPE(Element_t), POINTER :: Edge, Element
    LOGICAL, ALLOCATABLE :: NodeActive(:),Inside(:),Outside(:)       
    CHARACTER(*), PARAMETER :: Caller = 'ComponentStokesTheorem'
   
    CALL Info(Caller,'Calculating line integral for:'//TRIM(avar % name))
    
    TargetBodies => ListGetIntegerArray( VList,'Master Bodies',Found )
    IF( .NOT. Found ) TargetBodies => ListGetIntegerArray( VList,'Body',Found )
    IF( .NOT. Found ) CALL Fatal(Caller,'Stokes theorem requires > Master Bodies <') 

    HelperArray => ListGetConstRealArray( Vlist, 'Coil Center', UnfoundFatal = .TRUE. ) 
    Center(1:3) = HelperArray(1:3,1)

    HelperArray => ListGetConstRealArray( Vlist, 'Coil normal', UnfoundFatal = .TRUE. )
    Normal(1:3) = HelperArray(1:3,1)
    CALL TangentDirections(Normal, Tangent1, Tangent2)
    
    n = Mesh % NumberOfNodes 
    ALLOCATE(NodeActive(n))

    CALL SetActiveNodeSet()

    IF(Surf) THEN
      CALL ComputeCylinderIntegral()
    ELSE
      CALL ComputeLineIntegral()
    END IF

    FL = REAL(Circ)
    
  CONTAINS

    SUBROUTINE SetActiveNodeSet()

      ALLOCATE(Inside(n), Outside(n))
      Inside = .FALSE.
      Outside = .FALSE.    

      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF( ANY(TargetBodies == Element % BodyId) ) THEN
          Inside(Element % NodeIndexes) = .TRUE.
        ELSE
          Outside(Element % NodeIndexes) = .TRUE.
        END IF
      END DO

      !PRINT *,'Inside nodes:',COUNT(Inside)
      !PRINT *,'Outside nodes:',COUNT(Outside)

      ! We make the stokes theorem on the interface
      NodeActive = Inside .AND. Outside
      DEALLOCATE(Inside,Outside)

      n = COUNT(NodeActive)
      CALL Info(Caller,'Active nodes for edge path candidates: '//I2S(n))

      ! We may limit the active set of nodes through which path may be drawn.
      ! The idea could be to choose only the upper or lower surface of the coil. 
      htol = ListGetConstReal( Vlist,'Flux linkage height tolerance',GotHtol )
      rtol = ListGetConstReal( Vlist,'Flux linkage radius tolerance',GotRtol )
      IF(Found ) THEN      
        hgoal = ListGetConstReal( Vlist,'Flux linkage relative height',Found ) 
        hmin = HUGE(hmin)
        hmax = -HUGE(hmax)
        rmin = HUGE(rmin)
        rmax = -HUGE(rmax)
        DO k=1,2
          DO i=1,Mesh % NumberOfNodes
            IF(.NOT. NodeActive(i)) CYCLE
            Coord(1) = Mesh % Nodes % x(i)
            Coord(2) = Mesh % Nodes % y(i)
            Coord(3) = Mesh % Nodes % z(i)           
            h = SUM(Coord*Normal)
            r = SQRT(SUM(Coord*Tangent1)**2 + SUM(Coord*Tangent2)**2)

            IF(k==1) THEN
              hmin = MIN(h,hmin)
              hmax = MAX(h,hmax)
              rmin = MIN(r,rmin)
              rmax = MAX(r,rmax)
            ELSE
              IF(GotHtol) THEN
                hrel = (h-hmin)/(hmax-hmin)
                IF( ABS(hrel-hgoal) > htol ) NodeActive(i) = .FALSE.
              END IF
              IF(GotRtol) THEN
                IF( ABS(r-rmin)/rmin > rtol) NodeActive(i) = .FALSE.
              END IF
            END IF
          END DO

          IF(k==1) THEN
            IF( ParEnv % PEs > 1 ) THEN
              hmin = ParallelReduction(hmin,1)
              hmax = ParallelReduction(hmax,2)
              rmin = ParallelReduction(rmin,1)
              rmax = ParallelReduction(rmax,2)
            END IF
            !PRINT *,'Height interval in n-t system: ',hmin,hmax
            !PRINT *,'Radius interval in n-t system: ',hmin,hmax
          END IF
        END DO
        n = COUNT(NodeActive)
        CALL Info(Caller,'Active nodes after tolerance check: '//I2S(n))
      END IF
     
    END SUBROUTINE SetActiveNodeSet
      

    SUBROUTINE ComputeLineIntegral()

      ! Create a graph for node-to-edge connectivity
      !----------------------------------------------
      NodeGraph => AllocateMatrix()
      NodeGraph % FORMAT = MATRIX_LIST         
      DO i = Mesh % NumberOfEdges, 1, -1
        Edge => Mesh % Edges(i)    
        IF(ALL(NodeActive(Edge % NodeIndexes))) THEN
          DO j=1, Edge % TYPE % NumberOfNodes 
            CALL List_AddToMatrixElement( NodeGraph % ListMatrix,Edge % NodeIndexes(j),i,1.0_dp )
          END DO
        END IF
      END DO
      CALL List_ToCRSMatrix(NodeGraph)
      WRITE(Message,*) 'Nonzeros per row NodeGraph:',1.0_dp * SIZE(NodeGraph % Values) / NodeGraph % NumberOfRows
      CALL Info(Caller, Message)

      n0 = Mesh % NumberOfNodes

      SaveLoop = .FALSE.
      IF( SaveLoop ) THEN
        OPEN (10, FILE='Loop.dat' )
      END IF
      Debug = .FALSE.


      ! Find the node closest to the origin
      r2min = HUGE(r2min)
      r2ind = 0
      DO i=1,Mesh % NumberOfNodes
        IF(.NOT. NodeActive(i)) CYCLE
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        r2 = SUM((Coord-Center)**2)
        IF(r2 < r2min) THEN
          r2min = r2
          r2ind = i
          Coord0 = Coord
        END IF
      END DO

      IF( Debug ) THEN
        PRINT *,'Minimum distance:',SQRT(r2min),r2ind,NodeActive(r2ind)
        PRINT *,'Normal:',Normal
        PRINT *,'Tangent1:',Tangent1    
        PRINT *,'Tangent2:',Tangent2
        PRINT *,'Coord0:',Coord0
      END IF

      phisum = 0.0_dp
      nsteps = 0
      Circ = CMPLX(0.0_dp, 0.0_dp)
      ssum = 0.0_dp

      Coord1 = Coord0    
      i1 = r2ind 
      x1 = SUM(Coord1*Tangent1)
      y1 = SUM(Coord1*Tangent2) 
      phi1 = (180.0_dp/PI) * ATAN2(y1,x1)

      IF(SaveLoop) WRITE(10,*) nsteps, Coord1, phi1, phisum, ssum, REAL(Circ) 

      DO WHILE(.TRUE.)
        dsmax = 0.0_dp !-HUGE(dsmax) !0.0_dp
        kmax = 0

        ! Among the edges related to node "i1" find the one that has makes us further in
        ! minimizing the distance.
        DO j = NodeGraph % Rows(i1),NodeGraph % Rows(i1+1)-1
          k = NodeGraph % Cols(j)
          Edge => Mesh % Edges(k)
          NodeIndexes => Edge % NodeIndexes

          IF(NodeIndexes(1) == i1) THEN
            i2 = NodeIndexes(2)
          ELSE IF(NodeIndexes(2) == i1) THEN
            i2 = NodeIndexes(1)
          ELSE
            CALL Fatal(Caller,'Either node index in edge should be i1!')
          END IF
          IF(.NOT. NodeActive(i2)) CYCLE

          Coord(1) = Mesh % Nodes % x(i2)
          Coord(2) = Mesh % Nodes % y(i2)
          Coord(3) = Mesh % Nodes % z(i2)

          x2 = SUM(Coord*Tangent1)
          y2 = SUM(Coord*Tangent2)           
          phi = (180.0_dp/PI) * ATAN2(y2,x2)

          dphi = phi1-phi

          ! Deal with the discontinuity of angle around +/-180 deg
          IF(dphi < -180.0_dp) dphi = dphi+360.0_dp
          IF(dphi > 180.0_dp) dphi = dphi-360.0_dp

          IF(phisum < 315.0_dp ) THEN        
            ! This measure takes the shortest route at least in some case more robustly than some others.
            ds = dphi / SUM((Coord1-Coord)**2)
          ELSE
            ! After we only have ~45 degs left we find the node which approaches the starting node
            ! as well as possible. 
            ds = SUM( (Coord1-Coord0)**2 - (Coord-Coord0)**2 ) / SUM((Coord1-Coord)**2) 
          END IF

          IF( Debug ) THEN
            PRINT *,'kmax:',ds > dsmax, k,i2,ds,dphi,phi,Coord
          END IF

          ! Have we found a better edge candidate? Memorize that!
          IF( ds > dsmax ) THEN
            kmax = k
            imax = i2
            dsmax = ds
            Coord2 = Coord
            dphimax = dphi
            phi2 = phi
          END IF
        END DO

        ! When no way to get further we are node!
        IF( kmax == 0 ) THEN
          PRINT *,'Cands:',i1,NodeGraph % Rows(i1),NodeGraph % Rows(i1+1)-1

          CALL Fatal(Caller,'We had to stop because not route was found!')
          EXIT
        END IF

        nsteps = nsteps + 1

        ! Edge that goes to the minimum value
        k = kmax
        Edge => Mesh % Edges(k)    
        i2 = imax

        EdgeVector = Coord2 - Coord1

        ! Integration length to check optimality of the path
        ds = SQRT(SUM(EdgeVector**2))
        ssum = ssum + ds

        IF( MODULO(avar % dofs,3) == 0 ) THEN
          ! Integral over nodal field using the mean value
          j1 = avar % Perm(i1)
          j2 = avar % Perm(i2)          
          IF(j1==0 .OR. j2==0) CALL Fatal(Caller,'Nodal field missing on path!')
          DO k=1,3
            IF( avar % dofs == 3 ) THEN
              gradv(k) = ( avar % Values(3*(j1-1)+k) + avar % Values(3*(j2-1)+k) ) / 2
            ELSE
              gradv(k) = CMPLX(avar % Values(6*(j1-1)+k) + avar % Values(6*(j2-1)+k),&
                  avar % Values(6*(j1-1)+3+k) + avar % Values(6*(j2-1)+3+k) ) / 2
            END IF
          END DO
          Circ = Circ + SUM(gradv*EdgeVector)
        ELSE                
          ! Integral over edge field.

          ! Check the sign if the direction based on global edge direction rules
          ! If we do the path integral in the wrong direction compared to definition of edge switch the sign
          sgn = 1
          IF( ParEnv % PEs > 1 ) THEN                            
            i1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
            i2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
          END IF
          !IF(XOR(Edge % NodeIndexes(1) /= i1, i1 < i2) ) sgn = -sgn         
          IF(i1 > i2) sgn = -1

          j = avar % Perm(n0 + k)
          IF( j==0) CALL Fatal(Caller,'Edge field missing on path!')

          IF( avar % dofs == 1 ) THEN
            Circ = Circ + sgn * avar % Values(j)
          ELSE
            Circ = Circ + sgn * CMPLX( avar % Values(2*j-1),avar % Values(2*j) )
          END IF
        END IF

        ! Do not visit this point a 2nd time!
        !NodeActive(i1) = .FALSE.

        ! Continue from the end point
        i1 = imax
        Coord1 = Coord2
        phi1 = phi2

        phisum = phisum + dphimax      

        IF(Debug) PRINT *,'phisum:',nsteps,Coord1,phisum,ssum,phi1,REAL(Circ)


        IF( ABS(phisum) > 720.0_dp ) THEN
          CALL Fatal(Caller,'We circled twice around!?')
        END IF

        IF(SaveLoop) WRITE(10,*) nsteps, Coord1, phi1, phisum, ssum, REAL(Circ)

        ! We have come home to roost!
        IF(i1 == r2ind) THEN
          CALL Info(Caller,'Integration over full loop finished!')
          EXIT
        END IF
      END DO

      !PRINT *,'Path integral:',avar % dofs, targetbodies, nsteps, phisum, ssum, Circ

      CALL FreeMatrix(NodeGraph)
      IF(SaveLoop) CLOSE(10)
      
    END SUBROUTINE ComputeLineIntegral           


    SUBROUTINE ComputeCylinderIntegral()

      TYPE(Nodes_t), SAVE :: ElementNodes
      LOGICAL :: AllocationsDone = .FALSE.
      TYPE(GaussIntegrationPoints_t) :: IP
      REAL(KIND=dp) :: detJ, Area, s, TestVec(3)
      TYPE(ValueList_t), POINTER :: Params
      LOGICAL :: Stat, PiolaVersion, EdgeBasis 
      INTEGER :: np, nd, EdgeBasisDegree
      INTEGER, POINTER, SAVE :: Indexes(:)
      REAL(KIND=dp), POINTER, SAVE :: Basis(:), SOL(:,:), WBasis(:,:), dBasisdx(:,:), RotWBasis(:,:)

      CALL Info(Caller,'Estimating line integral from surface integral!')
      
      IF(.NOT. AllocationsDone ) THEN
        n = 2*Model % MaxElementNodes
        ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
            Basis(n), dBasisdx(n,3), WBasis(n,3), RotWBasis(n,3), SOL(6,n), Indexes(n) )      
        AllocationsDone = .TRUE.
      END IF

      Area = 0.0_dp
      Circ = 0.0_dp
      
      EdgeBasis = .FALSE.
      IF(avar % dofs <= 2) THEN
        EdgeBasis = .TRUE.
        Params => avar % Solver % Values

        CALL EdgeElementStyle(avar % Solver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
      END IF      
      
      DO t=1, Mesh % NumberOfFaces
        Element => Mesh % Faces(t)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes 

        IF(.NOT. ALL(NodeActive(NodeIndexes))) CYCLE
                
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

        dsmax = -HUGE(dsmax) 
        DO i=1,n
          Coord1(1) = ElementNodes % x(i)
          Coord1(2) = ElementNodes % y(i)
          Coord1(3) = ElementNodes % z(i)
          Coord1 = Coord1 - Center
          Coord1 = Coord1 - SUM(Normal*Coord1)*Normal
          DO j=i+1,n
            Coord2(1) = ElementNodes % x(j)
            Coord2(2) = ElementNodes % y(j)
            Coord2(3) = ElementNodes % z(j)
            Coord2 = Coord2 - Center
            Coord2 = Coord2 - SUM(Normal*Coord2)*Normal
            ds = SUM((Coord1-Coord2)**2)
            IF(ds > dsmax ) THEN
              dsmax = ds
              EdgeVector = Coord1-Coord2
            END IF
          END DO
        END DO
        EdgeVector = EdgeVector / SQRT(SUM(EdgeVector**2))

        TestVec = CrossProduct(Normal,Coord1)
        IF(SUM(TestVec*EdgeVector) > 0.0_dp) EdgeVector = -EdgeVector

        IF( EdgeBasis ) THEN
          nd = mGetElementDofs( Indexes, Uelement = Element, USolver = avar % Solver ) 
          np = COUNT(Indexes(1:nd) <= Mesh % NumberOfNodes)
          IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
              EdgeBasisDegree=EdgeBasisDegree)
        ELSE
          Indexes(1:n) = NodeIndexes(1:n)
          nd = n
          IP = GaussPoints( Element )          
        END IF
        
        DO i=1,avar % dofs 
          SOL(i,1:nd) = avar % values(avar % dofs*(avar % Perm(Indexes(1:nd))-1)+i)
        END DO
        
!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------

        DO l=1,IP % n
          IF(.NOT. EdgeBasis) THEN
            stat = ElementInfo( Element,ElementNodes,&
                IP % U(l),IP % V(l),IP % W(l), DetJ, Basis )             
          ELSE 
            stat = ElementInfo( Element, ElementNodes, IP % U(l), IP % V(l), &
                IP % W(l), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
                RotBasis = RotWBasis, USolver = avar % Solver )
          END IF

          !   stat = EdgeElementInfo(Element, ElementNodes, IP % U(l), IP % V(l), IP % W(l), &
          !       DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, dBasisdx = dBasisdx, &
          !       BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
          ! ELSE
          !   stat = ElementInfo(Element, ElementNodes, IP % U(l), IP % V(l), IP % W(l), &
          !       detJ, Basis, dBasisdx)           
          !   CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
          ! END IF
          
          s = DetJ * IP % s(l)            
          Area = Area + S

          SELECT CASE( avar % dofs )
          CASE( 1 )
            gradv = MATMUL(SOL(1,np+1:nd), WBasis(1:nd-np,:))
          CASE( 2 ) 
            gradv = CMPLX( MATMUL(SOL(1,np+1:nd), WBasis(1:nd-np,:)), &
                MATMUL(SOL(2,np+1:nd), WBasis(1:nd-np,:)))
          CASE( 3 )
            gradv = MATMUL(SOL(1:3,1:n),Basis(1:n))
          CASE( 6 )
            DO i=1,3
              gradv(i) = CMPLX( SUM(SOL(2*i-1,1:n)*Basis(1:n)), &
                  SUM(SOL(2*i,1:n)*Basis(1:n)) )
            END DO
          END SELECT
          
          Circ = Circ + s * SUM(gradv*EdgeVector)
        END DO
      END DO

      Circ = Circ / (hmax-hmin)       
      !PRINT *,'Path integral cyl:',avar % dofs, targetbodies, area, &
      !    area/((hmax-hmin)*2*PI), Circ
      
    END SUBROUTINE ComputeCylinderIntegral
          
  END FUNCTION ComponentStokesTheorem


!------------------------------------------------------------------------------
!> Given vector potential and current density compute the energy in the coil.
!> This is actually not energy, but twice the energy, since the values are
!> used to computed inductance matrix. 
!------------------------------------------------------------------------------
  FUNCTION ComponentCoilEnergy(Model, Mesh, MasterEntities, AVar, CVar, BCMode ) RESULT ( AIint ) 
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model    
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: MasterEntities(:) 
    TYPE(Variable_t), POINTER :: AVar, CVar
    LOGICAL, OPTIONAL :: BCMode 
    REAL(KIND=dp) :: AIint
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: t, i, j, k, l, n, np, nd, EdgeBasisDegree, t1, t2
    REAL(KIND=dp) :: volume
    LOGICAL :: Found
    LOGICAL :: Stat, PiolaVersion, EdgeBasis, DoBCs 
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    CHARACTER(*), PARAMETER :: Caller = 'ComponentCoilEnergy'
        
    IF(.NOT. ASSOCIATED(MasterEntities)) THEN
      CALL Fatal(Caller,'"MasterEntities" not associated!')
    END IF

    DoBCs = .FALSE.
    IF(PRESENT(BcMode)) DoBCs = BCMode 

    str = I2S(MasterEntities(1))
    DO i=2,SIZE(MasterEntities)
      str = TRIM(str)//' '//I2S(MasterEntities(i))
    END DO
    IF( DoBCs ) THEN
      CALL Info(Caller,'Performing reduction for BCs: '//TRIM(str),Level=10)
    ELSE
      CALL Info(Caller,'Performing reduction for bodies: '//TRIM(str),Level=10)
    END IF

    IF( DoBCs ) THEN
      t1 = Mesh % NumberOfBulkElements + 1
      t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      t1 = 1
      t2 = Mesh % NumberOfBulkElements
    END IF

    EdgeBasis = .FALSE.
    IF(avar % dofs <= 2) THEN
      EdgeBasis = .TRUE.      
      CALL EdgeElementStyle(avar % Solver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
    END IF

    AIint = 0.0_dp
    volume = 0.0_dp

    IF( CVar % Dofs /= 3 ) THEN
      CALL Fatal(Caller,'Expecting 3 components for current density!')
    END IF
    IF( ALL([1,3] /= AVar % Dofs) ) THEN
      CALL Fatal(Caller,'Expecting 1 or 3 components for vector potential!')
    END IF

    DO t=t1, t2
      Element => Mesh % Elements(t)
      IF( DoBCs ) THEN
        IF( ALL( MasterEntities /= Element % BoundaryInfo % Constraint ) ) CYCLE
      ELSE
        IF( ALL( MasterEntities /= Element % BodyId ) ) CYCLE
      END IF
      CALL LocalIntegElem()
    END DO

    !AIint = ParallelReduction( AIint ) 
    !Volume = ParallelReduction( volume ) 

    CALL Info(Caller,'Reduction operator finished',Level=12)

  CONTAINS
     
    SUBROUTINE LocalIntegElem()

      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t), SAVE :: ElementNodes
      LOGICAL, SAVE :: AllocationsDone = .FALSE.
      INTEGER, POINTER, SAVE :: EdgeIndexes(:)
      INTEGER, POINTER :: NodeIndexes(:), pIndexes(:)
      REAL(KIND=dp) :: DetJ,S,Cip(3),Aip(3)
      REAL(KIND=dp), POINTER, SAVE :: Basis(:), WBasis(:,:), dBasisdx(:,:), RotWBasis(:,:), &
          Aelem(:,:), Celem(:,:)    

      IF(.NOT. AllocationsDone ) THEN
        n = 2*Model % MaxElementNodes
        ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
            Basis(n), dBasisdx(n,3), WBasis(n,3), RotWBasis(n,3), Aelem(6,n), Celem(3,n), EdgeIndexes(n) )      
        AllocationsDone = .TRUE.
      END IF


      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes 

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

      IF( EdgeBasis ) THEN
        nd = mGetElementDofs( EdgeIndexes, Uelement = Element, USolver = avar % Solver ) 
        np = COUNT(EdgeIndexes(1:nd) <= Mesh % NumberOfNodes)
        IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
            EdgeBasisDegree=EdgeBasisDegree)
        pIndexes => EdgeIndexes
      ELSE
        nd = n
        IP = GaussPoints( Element )          
        pIndexes => NodeIndexes
      END IF

      ! Edge or nodal values of vector potential
      DO i=1,avar % dofs 
        Aelem(i,1:nd) = avar % values(avar % dofs*(avar % Perm(pIndexes(1:nd))-1)+i)
      END DO

      ! Nodal values of current density
      IF( cvar % TYPE == Variable_on_nodes_on_elements ) THEN
        pIndexes => Element % DGIndexes
      ELSE
        pIndexes => NodeIndexes
      END IF
      DO i=1,cvar % dofs
        Celem(i,1:n) = cvar % values(cvar % dofs*(cvar % Perm(pIndexes(1:n))-1)+i)
      END DO

      DO l=1,IP % n
        IF(.NOT. EdgeBasis) THEN
          stat = ElementInfo( Element,ElementNodes,&
              IP % U(l),IP % V(l),IP % W(l), DetJ, Basis )             
        ELSE 
          stat = ElementInfo( Element, ElementNodes, IP % U(l), IP % V(l), &
              IP % W(l), detJ, Basis, dBasisdx, EdgeBasis = WBasis, &
              RotBasis = RotWBasis, USolver = avar % Solver )
        END IF

        !stat = EdgeElementInfo(Element, ElementNodes, IP % U(l), IP % V(l), IP % W(l), &
        !      DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, dBasisdx = dBasisdx, &
        !      BasisDegree = EdgeBasisDegree, ApplyPiolaTransform = .TRUE.)
        !ELSE
        !  stat = ElementInfo(Element, ElementNodes, IP % U(l), IP % V(l), IP % W(l), &
        !      detJ, Basis, dBasisdx)           
        !  CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
        !END IF

        s = DetJ * IP % s(l)            

        ! Vector potential at IP
        SELECT CASE( avar % dofs )
        CASE( 1 )
          Aip = MATMUL(Aelem(1,np+1:nd), WBasis(1:nd-np,:))
          !CASE( 2 ) 
          !  Aip = CMPLX( MATMUL(Aelem(1,np+1:nd), WBasis(1:nd-np,:)), &
          !      MATMUL(Aelem(2,np+1:nd), WBasis(1:nd-np,:)))
        CASE( 3 )
          Aip = MATMUL(Aelem(1:3,1:n),Basis(1:n))
          !CASE( 6 ) 
          !  Aip = CMPLX( MATMUL(Aelem(1:2:5,1:n),Basis(1:n)), &
          !      MATMUL(Aelem(2:2:6,1:n),Basis(1:n)) )
        CASE DEFAULT
          CALL Fatal(Caller,'Invalid number of components for vector potential!')
        END SELECT

        ! Current density at IP
        Cip = MATMUL(Celem(1:3,1:n),Basis(1:n))

        AIint = AIint + s * SUM(Aip*Cip) 
        Volume = Volume + s
      END DO

    END SUBROUTINE LocalIntegElem
           
!------------------------------------------------------------------------------
  END FUNCTION ComponentCoilEnergy
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute integrals for determining the S parameters.
!------------------------------------------------------------------------------
  FUNCTION BoundaryWaveFlux(Model, Mesh, MasterEntities, Avar, InFlux, PortImp, PortBC ) &
      RESULT ( OutFlux ) 
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model    
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: MasterEntities(:) 
    TYPE(Variable_t), POINTER :: Avar
    COMPLEX(KIND=dp) :: OutFlux
    COMPLEX(KIND=dp) :: InFlux
    COMPLEX(KIND=dp) :: PortImp
    LOGICAL :: PortBC
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: t, i, j, k, l, n, np, nd, EdgeBasisDegree, t1, t2, bc_id
    LOGICAL :: Found, InitHandles
    LOGICAL :: Stat, PiolaVersion, EdgeBasis, UseGaussLaw
    TYPE(ValueList_t), POINTER :: BC
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    REAL(KIND=dp) :: area, omega
    COMPLEX(KIND=dp) :: int_norm, int_el, vol, curr, port_curr, trans, Zimp
    CHARACTER(*), PARAMETER :: Caller = 'BoundaryWaveFlux'

    area = 0.0_dp

    int_norm = 0.0_dp
    int_el = 0.0_dp

    vol = 0.0_dp
    curr = 0.0_dp
    port_curr = 0.0_dp
    trans = 0.0_dp
    
    IF(.NOT. ASSOCIATED(MasterEntities)) THEN
      CALL Fatal(Caller,'"MasterEntities" not associated!')
    END IF
    
    str = I2S(MasterEntities(1))
    DO i=2,SIZE(MasterEntities)
      str = TRIM(str)//' '//I2S(MasterEntities(i))
    END DO
    CALL Info(Caller,'Performing reduction for BCs: '//TRIM(str),Level=10)
    
    Omega = ListGetAngularFrequency(Found=Found)
    IF(.NOT. Found) CALL Fatal(Caller,'We need angular frequency!')
       
    t1 = Mesh % NumberOfBulkElements + 1
    t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    
    EdgeBasis = .FALSE.
    UseGaussLaw = .FALSE.

    IF(avar % dofs <= 2) THEN
      EdgeBasis = .TRUE.
      CALL EdgeElementStyle(avar % Solver % Values, PiolaVersion, BasisDegree = EdgeBasisDegree ) 
      UseGaussLaw = ListGetLogical(avar % solver % values, 'Use Gauss Law', Found)
    END IF
    
    OutFlux = 0.0_dp
    area = 0.0_dp

    IF( ALL([1,2,3,6] /= Avar % Dofs) ) THEN
      CALL Fatal(Caller,'Invalid number of components for vector potential!')
    END IF
    
    InitHandles = .TRUE.
    DO t=t1, t2
      Element => Mesh % Elements(t)
      bc_id = Element % BoundaryInfo % Constraint
      IF( ALL( MasterEntities /= bc_id ) ) CYCLE
      BC => Model % BCs(bc_id) % Values
      
      IF(EdgeBasis .AND. Element % Type % ElementCode > 300 ) THEN
        k = FindBoundaryFaceIndex(Mesh,Element)
        Element => Mesh % Faces(k)        
      END IF      
      IF(UseGaussLaw) THEN
        CALL LocalIntegBC_AV(BC, Element, InitHandles)
      ELSE
        CALL LocalIntegBC_E(BC, Element, InitHandles)
      END IF
    END DO

    Area = ParallelReduction(Area)
    trans = ParallelReduction(trans)
    Zimp = 1.0_dp / trans

    PortImp = Zimp
    
    IF( UseGaussLaw ) THEN
      vol = ParallelReduction(Vol)
      curr = ParallelReduction(curr)
      port_curr = ParallelReduction(port_curr)

      ! Still testing these!
      curr = -curr
      port_curr = -port_curr
      
      vol = vol / area

      PRINT *,'LumpedCurr av:',vol,curr,port_curr,Zimp *curr, CONJG(Zimp) * port_curr
      
      OutFlux = (vol + Zimp * curr) / (2*SQRT(REAL(Zimp)))
      InFlux = (vol - CONJG(Zimp) * port_curr  ) / (2*SQRT(REAL(Zimp)))      

      ! For now use just the average voltages
      OutFlux = vol !curr 
      InFlux = 1.0 !port_curr 
      PortImp = Zimp
    ELSE
      int_el = ParallelReduction(int_el)
      int_norm = ParallelReduction(int_norm)      
      OutFlux = int_el
      InFlux = int_norm 
      
      PRINT *,'LumpedCurr e:',int_el,int_norm,area,trans,Zimp

    END IF
    
    
    CALL Info(Caller,'Reduction operator finished',Level=12)
    
  CONTAINS

!-----------------------------------------------------------------------------
    SUBROUTINE LocalIntegBC_E( BC, Element, InitHandles )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: BC
      TYPE(Element_t), POINTER :: Element
      LOGICAL :: InitHandles
!------------------------------------------------------------------------------
      COMPLEX(KIND=dp) :: B, Zs, L(3), muinv, MagLoad(3), TemGrad(3), eps, &
          e_ip(3), e_ip_norm, e_ip_tan(3), f_ip_tan(3), imu, phi, eps0, mu0inv, epsr, mur
      REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),WBasis(:,:),RotWBasis(:,:), e_local(:,:)
      REAL(KIND=dp) :: weight, DetJ, Normal(3), cond, u, v, w, x, y, z, rob0
      TYPE(Nodes_t), SAVE :: ElementNodes, ParentNodes
      INTEGER, POINTER :: NodeIndexes(:), pIndexes(:), ParentIndexes(:)
      INTEGER, POINTER, SAVE :: EdgeIndexes(:)
      LOGICAL :: Stat, Found
      TYPE(GaussIntegrationPoints_t) :: IP
      INTEGER :: t, i, j, m, np, p, q, ndofs, n, nd
      LOGICAL :: AllocationsDone = .FALSE.
      TYPE(Element_t), POINTER :: Parent
      TYPE(ValueHandle_t), SAVE :: MagLoad_h, ElRobin_h, MuCoeff_h, Absorb_h, TemRe_h, TemIm_h
      TYPE(ValueHandle_t), SAVE :: CondCoeff_h, CurrDens_h, EpsCoeff_h
      INTEGER :: nactive
      
      SAVE AllocationsDone, WBasis, RotWBasis, Basis, dBasisdx, e_local, mu0inv, eps0
      
      ndofs = avar % dofs
      IF(.NOT. AllocationsDone ) THEN
        m = Mesh % MaxElementDOFs
        ALLOCATE( ElementNodes % x(m), ElementNodes % y(m), ElementNodes % z(m), EdgeIndexes(m), &
            ParentNodes % x(m), ParentNodes % y(m), ParentNodes % z(m), &
            WBasis(m,3), RotWBasis(m,3), Basis(m), dBasisdx(m,3), e_local(ndofs,m) )      
        AllocationsDone = .TRUE.
      END IF
      
      IF( InitHandles ) THEN
        CALL ListInitElementKeyword( ElRobin_h,'Boundary Condition','Electric Robin Coefficient',InitIm=.TRUE.)
        CALL ListInitElementKeyword( MagLoad_h,'Boundary Condition','Magnetic Boundary Load', InitIm=.TRUE.,InitVec3D=.TRUE.)
        CALL ListInitElementKeyword( Absorb_h,'Boundary Condition','Absorbing BC')
        CALL ListInitElementKeyword( TemRe_h,'Boundary Condition','TEM Potential')
        CALL ListInitElementKeyword( TemIm_h,'Boundary Condition','TEM Potential Im')

        CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
        CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
        CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
        Found = .FALSE.
        IF( ASSOCIATED( Model % Constants ) ) THEN
          mu0inv = ListGetConstReal( Model % Constants,'Permeability of Vacuum', Found )
          IF(mu0inv/=0) mu0inv = 1/mu0inv;
        END IF
        IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
        Found = .FALSE.
        IF( ASSOCIATED( Model % Constants ) ) THEN
          eps0 = ListGetConstReal ( Model % Constants,'Permittivity of Vacuum', Found )
        END IF
        IF(.NOT. Found ) eps0 = 8.854187817d-12           
        InitHandles = .FALSE.
      END IF

      imu = CMPLX(0.0_dp, 1.0_dp)
      rob0 = Omega * SQRT( eps0 / mu0inv )
      
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes 

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

      Parent => Element % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(Parent)) Parent => Element % BoundaryInfo % Right
      IF(.NOT. ASSOCIATED( Parent ) ) THEN
        CALL Fatal(Caller,'Model lumping requires parent element!')
      END IF

      IF( EdgeBasis ) THEN
        np = Parent % TYPE % NumberOfNodes
        ParentIndexes => Parent % NodeIndexes
        ParentNodes % x(1:np) = Mesh % Nodes % x(ParentIndexes(1:np))
        ParentNodes % y(1:np) = Mesh % Nodes % y(ParentIndexes(1:np))
        ParentNodes % z(1:np) = Mesh % Nodes % z(ParentIndexes(1:np))

        nd = mGetElementDofs( EdgeIndexes, Uelement = Parent, USolver = avar % Solver ) 
        np = COUNT(EdgeIndexes(1:nd) <= Mesh % NumberOfNodes)
        pIndexes => EdgeIndexes

        IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
            EdgeBasisDegree=EdgeBasisDegree)
      ELSE
        nd = n
        IP = GaussPoints( Element )          
        IF( avar % TYPE == Variable_on_nodes_on_elements ) THEN
          pIndexes => Element % DGIndexes
        ELSE         
          pIndexes => NodeIndexes
        END IF
      END IF

      ! Edge or nodal values of vector potential
      DO i=1,avar % dofs 
        e_local(i,1:nd) = avar % values(avar % dofs*(avar % Perm(pIndexes(1:nd))-1)+i)
      END DO
      
      Normal = NormalVector(Element, ElementNodes, Check=.TRUE.)

      ! Numerical integration:
      !-----------------------      
      DO t=1,IP % n  
        
        stat = ElementInfo( Element, ElementNodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )              
        weight = IP % s(t) * detJ

        ! Get material properties from parent element.
        !----------------------------------------------
        mur = ListGetElementComplex( MuCoeff_h, Basis, Parent, Found, GaussPoint = t )      
        IF( .NOT. Found ) mur = 1.0_dp
        muinv = mur * mu0inv

        epsr = ListGetElementComplex( EpsCoeff_h, Basis, Parent, Found, GaussPoint = t )      
        IF( .NOT. Found ) epsr = 1.0_dp
        eps = epsr * eps0
        
        Cond = ListGetElementReal( CondCoeff_h, Basis, Parent, Found, GaussPoint = t )
        
        IF( ListGetElementLogical( Absorb_h, Element, Found ) ) THEN
          B = imu * rob0 * SQRT( epsr / mur ) 
        ELSE        
          B = ListGetElementComplex( ElRobin_h, Basis, Element, Found, GaussPoint = t )
        END IF
                  
        Zs = 1.0_dp / (SQRT(REAL(muinv*eps)))

        MagLoad = ListGetElementComplex3D( MagLoad_h, Basis, Element, Found, GaussPoint = t )
        TemGrad = CMPLX( ListGetElementRealGrad( TemRe_h,dBasisdx,Element,Found), &
            ListGetElementRealGrad( TemIm_h,dBasisdx,Element,Found) )
        L = ( MagLoad + TemGrad ) / ( 2*B) 
                
        IF( EdgeBasis ) THEN
          ! In order to get the normal component of the electric field we must operate on the
          ! parent element. The surface element only has tangential components. 
          CALL FindParentUVW( Element, n, Parent, Parent % TYPE % NumberOfNodes, U, V, W, Basis ) 
          stat = ElementInfo( Parent, ParentNodes, u, v, w, detJ, Basis, dBasisdx, &
              EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = avar % Solver )
          e_ip(1:3) = CMPLX(MATMUL(e_local(1,np+1:nd),WBasis(1:nd-np,1:3)), MATMUL(e_local(2,np+1:nd),WBasis(1:nd-np,1:3)))       
        ELSE
          DO i=1,3
            e_ip(i) = CMPLX( SUM( Basis(1:n) * e_local(i,1:n) ), SUM( Basis(1:n) * e_local(i+3,1:n) ) )
          END DO
        END IF
        
        e_ip_norm = SUM(e_ip*Normal)
        e_ip_tan = e_ip - e_ip_norm * Normal

        ! Integral over electric field: This gives the phase
        int_el = int_el + weight * SUM(e_ip_tan * CONJG(L) )         

        ! Norm of electric field used for normalization
        int_norm = int_norm + weight * ABS( SUM( L * CONJG(L) ) ) 

        trans = trans + B * weight / Omega                
        area = area + weight        
      END DO
      
!------------------------------------------------------------------------------
    END SUBROUTINE LocalIntegBC_E
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
    SUBROUTINE LocalIntegBC_AV( BC, Element, InitHandles )
!------------------------------------------------------------------------------
      TYPE(ValueList_t), POINTER :: BC
      TYPE(Element_t), POINTER :: Element
      LOGICAL :: InitHandles
!------------------------------------------------------------------------------
      COMPLEX(KIND=dp) :: tc_ip, cd_ip, v_ip, ep_ip, eps0, eps, mu0inv, muinv, mur, epsr, &
          cond_ip, imu
      REAL(KIND=dp), ALLOCATABLE :: Basis(:),dBasisdx(:,:),v_local(:,:)
      REAL(KIND=dp) :: weight, DetJ 
      TYPE(Nodes_t), SAVE :: ElementNodes 
      INTEGER, POINTER :: NodeIndexes(:), pIndexes(:), ParentIndexes(:)
      LOGICAL :: Found
      TYPE(GaussIntegrationPoints_t) :: IP
      INTEGER :: t, i, j, m, np, p, q, ndofs, n, nd
      LOGICAL :: AllocationsDone = .FALSE.
      TYPE(Element_t), POINTER :: Parent, MatElement
      TYPE(ValueHandle_t), SAVE :: MuCoeff_h, EpsCoeff_h, CondCoeff_h, ExtPot_h
      TYPE(ValueHandle_t), SAVE :: TransferCoeff_h, ElCurrent_h, BCMat_h
      
      SAVE AllocationsDone, Basis, dBasisdx, v_local, mu0inv, eps0
      
      ndofs = avar % dofs
      IF(.NOT. AllocationsDone ) THEN
        m = Mesh % MaxElementDOFs
        ALLOCATE( ElementNodes % x(m), ElementNodes % y(m), ElementNodes % z(m), &
            Basis(m), dBasisdx(m,3), v_local(ndofs,m) )
        AllocationsDone = .TRUE.
      END IF

      ! BC given with these:
      ! Electric Transfer Coefficient      
      ! Electric Current Density / Incident Voltage
      IF( InitHandles ) THEN
        CALL ListInitElementKeyword( MuCoeff_h,'Material','Relative Reluctivity',InitIm=.TRUE.)      
        CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=.TRUE.)
        CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
        
        CALL ListInitElementKeyword( TransferCoeff_h,'Boundary Condition','Electric Transfer Coefficient',InitIm=.TRUE.)
        CALL ListInitElementKeyword( ElCurrent_h,'Boundary Condition','Electric Current Density',InitIm=.TRUE.)
        CALL ListInitElementKeyword( ExtPot_h,'Boundary Condition','Incident Voltage',InitIm=.TRUE.)

        CALL ListInitElementKeyword( BCMat_h,'Boundary Condition','Material')

        Found = .FALSE.
        IF( ASSOCIATED( Model % Constants ) ) THEN
          mu0inv = 1.0_dp / ListGetConstReal( Model % Constants,'Permeability of Vacuum', Found )
        END IF
        IF(.NOT. Found ) mu0inv = 1.0_dp / ( PI * 4.0d-7 )
        Found = .FALSE.
        IF( ASSOCIATED( Model % Constants ) ) THEN
          eps0 = ListGetConstReal ( Model % Constants,'Permittivity of Vacuum', Found )
        END IF
        IF(.NOT. Found ) eps0 = 8.854187817d-12           
        InitHandles = .FALSE.
      END IF

      imu = CMPLX(0.0_dp, 1.0_dp)
      
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes 

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

      Parent => Element % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(Parent)) Parent => Element % BoundaryInfo % Right
      IF(.NOT. ASSOCIATED( Parent ) ) THEN
        CALL Fatal(Caller,'Model lumping requires parent element!')
      END IF

      ! If we have material also define in the BC section then use it
      i = ListGetElementInteger( BCMat_h, Element, Found )
      IF( i > 0 ) THEN
        MatElement => Element
      ELSE
        MatElement => Parent
      END IF
              
      nd = n
      IP = GaussPoints( Element )          
      pIndexes => NodeIndexes

      ! Nodal values of scalar potential
      DO i=1,avar % dofs 
        v_local(i,1:nd) = avar % values(avar % dofs*(avar % Perm(pIndexes(1:nd))-1)+i)
      END DO
      
      ! Numerical integration:
      !-----------------------      
      DO t=1,IP % n  
        
        stat = ElementInfo( Element, ElementNodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )              
        weight = IP % s(t) * detJ

        ! Get material properties from parent element.
        !----------------------------------------------
        mur = ListGetElementComplex( MuCoeff_h, Basis, MatElement, Found, GaussPoint = t )      
        IF( .NOT. Found ) mur = 1.0_dp
        muinv = mur * mu0inv

        epsr = ListGetElementComplex( EpsCoeff_h, Basis, MatElement, Found, GaussPoint = t )      
        IF( .NOT. Found ) epsr = 1.0_dp
        eps = epsr * eps0

        cond_ip = ListGetElementReal( CondCoeff_h, Basis, MatElement, Found, GaussPoint = t )        
        cd_ip = ListGetElementComplex( ElCurrent_h, Basis, Element, Found, GaussPoint = t )

        tc_ip = ListGetElementComplex( TransferCoeff_h, Basis, Element, Found, GaussPoint = t )
        IF(Found) THEN
          ep_ip = ListGetElementComplex( ExtPot_h, Basis, Element, Found, GaussPoint = t )
          IF(Found) cd_ip = cd_ip + 2 * tc_ip * ep_ip
        END IF
        v_ip = CMPLX( SUM( Basis(1:n) * v_local(1,1:n) ), SUM( Basis(1:n) * v_local(2,1:n) ) )
                
        area = area + weight

        curr = curr - tc_ip * v_ip * weight
        port_curr = port_curr + cd_ip * weight
        trans = trans + tc_ip * weight         
        vol = vol + v_ip * weight          
      END DO
      
!------------------------------------------------------------------------------
    END SUBROUTINE LocalIntegBC_AV
!------------------------------------------------------------------------------
    
  END FUNCTION BoundaryWaveFlux

    
END MODULE LumpingUtils
!------------------------------------------------------------------------------


!> \}


