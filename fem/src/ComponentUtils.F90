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


MODULE ComponentUtils

   USE ElementUtils
   USE ModelDescription

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
     LOGICAL, ALLOCATABLE :: VisitedNode(:)
     REAL(KIND=dp) :: Origin(3), Axis(3), P(3), F(3), v1(3), v2(3)
     REAL(KIND=dp), POINTER :: Pwrk(:,:)
     INTEGER :: t, i, j, k, dofs, globalnode
     LOGICAL :: ElementalVar, Found, NeedLocation
     INTEGER, POINTER :: MasterEntities(:),NodeIndexes(:),DofIndexes(:)
     LOGICAL :: VisitNodeOnlyOnce     
     INTEGER :: FirstElem, LastElem
     LOGICAL :: BcMode 

     
     CALL Info('ComponentNodalForceReduction','Performing reduction for component: '&
         //TRIM(ListGetString(CompParams,'Name')),Level=10)

     IF(.NOT. (PRESENT(Torque) .OR. PRESENT(Moment) .OR. PRESENT(Force) ) ) THEN
       CALL Warn('ComponentNodalForceReduction','Nothing to compute!')
       RETURN
     END IF

     IF( PRESENT(Torque)) Torque = 0.0_dp
     IF( PRESENT(Moment)) Moment = 0.0_dp
     IF( PRESENT(Force)) Force = 0.0_dp

     BcMode = .FALSE.
     MasterEntities => ListGetIntegerArray( CompParams,'Master Bodies',Found )     
     IF( .NOT. Found ) THEN
       MasterEntities => ListGetIntegerArray( CompParams,'Master Boundaries',Found ) 
       BcMode = .TRUE.
     END IF

     IF(.NOT. Found ) THEN
       CALL Warn('ComponentNodalForceReduction',&
           '> Master Bodies < or > Master Boundaries < not given')
       RETURN
     END IF

     NeedLocation = PRESENT( Moment ) .OR. PRESENT( Torque )

     ! User may specific origin and axis for torque computation
     ! By default (0,0,0) is the origin, and (0,0,1) the axis. 
     Pwrk => ListGetConstRealArray( CompParams,'Torque Origin',Found )
     IF( Found ) THEN
       IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
         CALL Fatal('ComponentNodalForceReduction','Size of > Torque Origin < should be 3!')
       END IF
       Origin = Pwrk(1:3,1)
     ELSE
       Origin = 0.0_dp
     END IF
     Pwrk => ListGetConstRealArray( CompParams,'Torque Axis',Found )
     IF( Found ) THEN
       IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
         CALL Fatal('ComponentNodalForceReduction','Size of > Torque Axis < should be 3!')
       END IF
       Axis = Pwrk(1:3,1)
     ELSE
       Axis = 0.0_dp    
       Axis(3) = 1.0_dp  
     END IF

     ElementalVar = ( NF % TYPE == Variable_on_nodes_on_elements )
     IF( PRESENT( SetPerm ) .AND. .NOT. ElementalVar ) THEN
       CALL Fatal('ComponentNodalForceReduction','SetPerm is usable only for elemental fields')
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
         k = NF % Perm(j)
         IF( k == 0 ) CYCLE

         IF( VisitNodeOnlyOnce ) THEN
           IF( PRESENT( SetPerm ) ) j = SetPerm(j)
           IF( VisitedNode(j) ) CYCLE
           VisitedNode(j) = .TRUE.
         END IF

         globalnode = NodeIndexes(i)

         ! Only compute the parallel reduction once
         IF( ParEnv % PEs > 1 ) THEN
           IF( Mesh % ParallelInfo % NeighbourList(globalnode) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
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
             v1 = (1.0_dp - SUM(Axis*v1) ) * v1
             v2 = CrossProduct(v1,F)
             Torque = Torque + SUM(Axis*v2)        
           END IF
         END IF

       END DO
     END DO

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
     CHARACTER(LEN=MAX_NAME_LEN) :: OperName
     REAL(KIND=dp) :: OperX
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     LOGICAL, ALLOCATABLE :: VisitedNode(:)
     INTEGER :: t, i, j, k, NoDofs, globalnode, sumi
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
    
    
    sumi = NINT( ParallelReduction(1.0_dp * sumi) )
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
     CHARACTER(LEN=MAX_NAME_LEN) :: OperName, CoeffName
     LOGICAL :: GotCoeff
     REAL(KIND=dp) :: OperX
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: t, i, j, k, NoDofs
     INTEGER, POINTER :: NodeIndexes(:), DofIndexes(:)
     REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,Grad(3)
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
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName, VarName, CoeffName
    LOGICAL :: GotVar, GotOper, GotCoeff, VectorResult
    TYPE(ValueList_t), POINTER :: CompParams
    REAL(KIND=dp) :: ScalarVal, VectorVal(3)
    TYPE(Variable_t), POINTER :: Var

    CALL Info('UpdateDepedentComponents','Updating Components to reflect new solution',Level=6)

    DO j=1,CurrentModel % NumberOfComponents

      IF( ALL( ComponentList /= j ) ) CYCLE

      CALL Info('UpdateDepedentComponents','Updating component: '//TRIM(I2S(j)))
      CompParams => CurrentModel % Components(j) % Values


      NoVar = 0
      DO WHILE( .TRUE. )
        NoVar = NoVar + 1
        OperName = ListGetString( CompParams,'Operator '//TRIM(I2S(NoVar)), GotOper)
        VarName = ListGetString( CompParams,'Variable '//TRIM(I2S(NoVar)), GotVar)
        CoeffName = ListGetString( CompParams,'Coeffcient '//TRIM(I2S(NoVar)), GotCoeff)
        IF(.NOT. (GotVar .AND. GotOper ) ) EXIT

        Var => VariableGet( CurrentModel % Mesh % Variables, VarName ) 
        IF( .NOT. ASSOCIATED( Var ) ) THEN
          CALL Info('UpdateDependentComponents','Variable not available: '//TRIM(VarName))
          CYCLE
        END IF
        VectorResult = .FALSE.

        SELECT CASE( OperName ) 

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
                //TRIM(I2S(i))//': ',ScalarVal
            CALL Info('UpdateDependentComponents',Message,Level=5)
            CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName)//': '&
                //TRIM(VarName)//' '//TRIM(I2S(i)),VectorVal(i) )                        
          END DO
        ELSE          
          WRITE( Message,'(A,ES15.6)') TRIM(OperName)//': '//TRIM(VarName)//': ',ScalarVal
          CALL Info('UpdateDependentComponents',Message,Level=5)
          CALL ListAddConstReal( CompParams,'res: '//TRIM(OperName)//' '//TRIM(VarName),ScalarVal )           
        END IF

      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateDependentComponents
!------------------------------------------------------------------------------


 END MODULE ComponentUtils
!------------------------------------------------------------------------------


!> \}


