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
! *  Authors: Juha Ruokolainen, Antti Pursula
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Aug 2002
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. StatCurrentSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver_Init( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL ::  TransientSimulation
    REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
    LOGICAL :: Found, Calculate
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    INTEGER :: dim

    Params => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF (ListGetLogical(Params,'Calculate Joule Heating',Found)) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Joule Heating' )

    IF (ListGetLogical(Params,'Calculate Nodal Heating',Found)) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Nodal Joule Heating' )
    
    Calculate = ListGetLogical(Params,'Calculate Volume Current',Found)
    IF( Calculate ) THEN
      IF( Dim == 2 ) THEN
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'Volume Current[Volume Current:2]' )
      ELSE
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
            'Volume Current[Volume Current:3]' )
      END IF
    END IF

    ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVIVTY
    CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
        
!------------------------------------------------------------------------------
END SUBROUTINE StatCurrentSolver_Init
!------------------------------------------------------------------------------
    
!------------------------------------------------------------------------------
!>  Solve the Poisson equation for the electric potential and compute the 
!>  volume current and Joule heating
!------------------------------------------------------------------------------
  SUBROUTINE StatCurrentSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Differentials
     USE Adaptive
     IMPLICIT NONE
!------------------------------------------------------------------------------ 
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET:: Solver
     REAL (KIND=DP) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER  :: StiffMatrix
     TYPE(Element_t), POINTER :: CurrentElement
     TYPE(Nodes_t) :: ElementNodes

     REAL (KIND=DP), POINTER :: ForceVector(:), Potential(:)
     REAL (KIND=DP), POINTER :: ElField(:), VolCurrent(:)
     REAL (KIND=DP), POINTER :: Heating(:), NodalHeating(:)
     REAL (KIND=DP), POINTER :: EleC(:)
     REAL (KIND=DP), POINTER :: Cwrk(:,:,:)
     REAL (KIND=DP), ALLOCATABLE ::  Conductivity(:,:,:), &
       LocalStiffMatrix(:,:), Load(:), LocalForce(:)

     REAL (KIND=DP) :: Norm, HeatingTot, VolTot, CurrentTot, ControlTarget, ControlScaling = 1.0
     REAL (KIND=DP) :: Resistance, PotDiff
     REAL (KIND=DP) :: at, st, at0

     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: PotentialPerm(:)
     INTEGER :: i, j, k, n, t, istat, bf_id, LocalNodes, Dim, &
         iter, NonlinearIter
 
     LOGICAL :: AllocationsDone = .FALSE., gotIt, FluxBC
     LOGICAL :: CalculateField = .FALSE., ConstantWeights
     LOGICAL :: CalculateCurrent, CalculateHeating, CalculateNodalHeating
     LOGICAL :: ControlPower, ControlCurrent, Control

     TYPE(ValueList_t), POINTER :: Params
     TYPE(Variable_t), POINTER :: Var

     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName

     LOGICAL :: GetCondAtIp
     TYPE(ValueHandle_t) :: CondAtIp_h
     REAL(KIND=dp) :: CondAtIp
     
     SAVE LocalStiffMatrix, Load, LocalForce, &
          ElementNodes, CalculateCurrent, CalculateHeating, &
          AllocationsDone, VolCurrent, Heating, Conductivity, &
          CalculateField, ConstantWeights, &
          Cwrk, ControlScaling, CalculateNodalHeating

     INTERFACE
       FUNCTION StatCurrentSolver_Boundary_Residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
         USE Types
         TYPE(Element_t), POINTER :: Edge
         TYPE(Model_t) :: Model
         TYPE(Mesh_t), POINTER :: Mesh
         REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
         INTEGER :: Perm(:)
       END FUNCTION StatCurrentSolver_Boundary_Residual

       FUNCTION StatCurrentSolver_Edge_Residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
         USE Types
         TYPE(Element_t), POINTER :: Edge
         TYPE(Model_t) :: Model
         TYPE(Mesh_t), POINTER :: Mesh
         REAL(KIND=dp) :: Quant(:), Indicator(2)
         INTEGER :: Perm(:)
       END FUNCTION StatCurrentSolver_Edge_Residual

       FUNCTION StatCurrentSolver_Inside_Residual(Model, Element, Mesh, Quant, Perm, Fnorm) RESULT(Indicator)
         USE Types
         TYPE(Element_t), POINTER :: Element
         TYPE(Model_t) :: Model
         TYPE(Mesh_t), POINTER :: Mesh
         REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
         INTEGER :: Perm(:)
       END FUNCTION StatCurrentSolver_Inside_Residual
     END INTERFACE

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF(.NOT.ASSOCIATED(Solver % Matrix)) RETURN

     Potential     => Solver % Variable % Values
     PotentialPerm => Solver % Variable % Perm
     Params => GetSolverParams()

     LocalNodes = Model % NumberOfNodes
     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

     Norm = Solver % Variable % Norm
     DIM = CoordinateSystemDimension()

     ControlTarget = GetCReal( Params,'Power Control',ControlPower)
     IF(ControlPower) THEN
       ControlCurrent = .FALSE.
     ELSE
       ControlTarget = GetCReal( Params,'Current Control',ControlCurrent)
     END IF
     Control = ControlPower .OR. ControlCurrent

     ! To obtain convergence rescale the potential to the original BCs
     IF( Control ) THEN
       Potential = Potential / ControlScaling
       Solver % Variable % Norm = Solver % Variable % Norm / ControlScaling
     END IF

     NonlinearIter = ListGetInteger( Params, &
         'Nonlinear System Max Iterations', GotIt )
     IF ( .NOT. GotIt ) NonlinearIter = 1

     GetCondAtIp = ListGetLogical( Params,'Conductivity At Ip',GotIt )
     
!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged ) THEN
       N = Model % MaxElementNodes
 
       IF(AllocationsDone) THEN
         DEALLOCATE( ElementNodes % x, &
                   ElementNodes % y,   &
                   ElementNodes % z,   &
                   Conductivity,       &
                   LocalForce,         &
                   LocalStiffMatrix,   &
                   Load )
       END IF

       ALLOCATE( ElementNodes % x(N),   &
                 ElementNodes % y(N),   &
                 ElementNodes % z(N),   &
                 Conductivity(3,3,N),   &
                 LocalForce(N),         &
                 LocalStiffMatrix(N,N), &
                 Load(N),               &
                 STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'StatCurrentSolve', 'Memory allocation error.' )
       END IF

       NULLIFY( Cwrk )
 
       CalculateCurrent = ListGetLogical( Params, &
           'Calculate Volume Current', GotIt )
       IF ( CalculateCurrent ) THEN
         Var => VariableGet( Solver % Mesh % Variables,'Volume Current')
         IF( ASSOCIATED( Var) ) THEN
           VolCurrent => Var % Values
         ELSE
           CALL Fatal('StatCurrentSolver','Volume Current does not exist')
         END IF
       END IF
        
       CalculateHeating = ListGetLogicalAnyEquation( &
           Model,'Calculate Joule heating')
       IF ( .NOT. CalculateHeating )  &
           CalculateHeating = ListGetLogical( Params, &
           'Calculate Joule Heating', GotIt )
       IF ( CalculateHeating ) THEN
         Var => VariableGet( Solver % Mesh % Variables,'Joule Heating')
         IF( ASSOCIATED( Var) ) THEN
           Heating => Var % Values
         ELSE
           CALL Fatal('StatCurrentSolver','Joule Heating does not exist')
         END IF
       END IF

       CalculateNodalHeating = ListGetLogical( Params, &
           'Calculate Nodal Heating', GotIt )
       IF ( CalculateNodalHeating ) THEN
         Var => VariableGet( Solver % Mesh % Variables,'Nodal Joule Heating')
         IF( ASSOCIATED( Var) ) THEN
           NodalHeating => Var % Values
         ELSE
           CALL Fatal('StatCurrentSolver','Nodal Joule Heating does not exist')
         END IF
       END IF


       
       ConstantWeights = ListGetLogical( Params, &
           'Constant Weights', GotIt )

!------------------------------------------------------------------------------

       IF ( .NOT.ASSOCIATED( StiffMatrix % MassValues ) ) THEN
         ALLOCATE( StiffMatrix % Massvalues( LocalNodes ) )
         StiffMatrix % MassValues = 0.0d0
       END IF

!------------------------------------------------------------------------------
!      Add electric field to the variable list (disabled)
!------------------------------------------------------------------------------
       IF ( CalculateField ) THEN         
          CALL VariableAddVector( Solver % Mesh % Variables, Solver % Mesh, &
               Solver, 'Electric Field', dim, ElField, PotentialPerm)
       END IF
          
       AllocationsDone = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     EquationName = ListGetString( Params, 'Equation' )

     CALL Info( 'StatCurrentSolve', '-------------------------------------',Level=4 )
     CALL Info( 'StatCurrentSolve', 'STAT CURRENT SOLVER:  ', Level=4 )
     CALL Info( 'StatCurrentSolve', '-------------------------------------',Level=4 )

     CALL DefaultStart()
     
     DO iter = 1, NonlinearIter
       at  = CPUTime()
       at0 = RealTime()

       IF ( NonlinearIter > 1 ) THEN
         WRITE( Message, '(a,I0)' ) 'Static current iteration: ', iter
         CALL Info( 'StatCurrentSolve', Message, LEVEL=4 )
       END IF
       CALL Info( 'StatElecSolve', 'Starting Assembly...', Level=6 )

       CALL DefaultInitialize()

       !------------------------------------------------------------------------------

       !------------------------------------------------------------------------------
       !    Do the assembly
       !------------------------------------------------------------------------------

       IF( GetCondAtIp ) THEN
         CALL ListInitElementKeyword( CondAtIp_h,'Material','Electric Conductivity')
       END IF
         
       
       DO t = 1, Solver % NumberOfActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
               (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'

           CALL Info( 'StatCurrentSolve', Message, Level=5 )

           at0 = RealTime()
         END IF

         !------------------------------------------------------------------------------
         !        Check if this element belongs to a body where potential
         !        should be calculated
         !------------------------------------------------------------------------------
         CurrentElement => GetActiveElement(t)
         NodeIndexes => CurrentElement % NodeIndexes

         n = GetElementNOFNodes()

         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
         !------------------------------------------------------------------------------

         bf_id = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
             Values, 'Body Force', gotIt, minv=1, maxv=Model % NumberOfBodyForces )

         Load  = 0.0d0
         IF ( gotIt ) THEN
           Load(1:n) = ListGetReal( Model % BodyForces(bf_id) % Values, &
               'Current Source',n,NodeIndexes, Gotit )
         END IF

         IF( .NOT. GetCondAtIp ) THEN

           k = ListGetInteger( Model % Bodies(CurrentElement % BodyId) % &
               Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

           !------------------------------------------------------------------------------
           !      Read conductivity values (might be a tensor)
           !------------------------------------------------------------------------------
           
           CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Electric Conductivity', Cwrk, n, NodeIndexes )

           Conductivity = 0.0d0
           IF ( SIZE(Cwrk,1) == 1 ) THEN
             DO i=1,3
               Conductivity( i,i,1:n ) = Cwrk( 1,1,1:n )
             END DO
           ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
             DO i=1,MIN(3,SIZE(Cwrk,1))
               Conductivity(i,i,1:n) = Cwrk(i,1,1:n)
             END DO
           ELSE
             DO i=1,MIN(3,SIZE(Cwrk,1))
               DO j=1,MIN(3,SIZE(Cwrk,2))
                 Conductivity( i,j,1:n ) = Cwrk(i,j,1:n)
               END DO
             END DO
           END IF
         END IF
           
         !------------------------------------------------------------------------------
         !      Get element local matrix, and rhs vector
         !------------------------------------------------------------------------------
         CALL StatCurrentCompose( LocalStiffMatrix,LocalForce, &
             Conductivity,Load,CurrentElement,n,ElementNodes )
         !------------------------------------------------------------------------------
         !      Update global matrix and rhs vector from local matrix & vector
         !------------------------------------------------------------------------------

         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

         !------------------------------------------------------------------------------
       END DO
       CALL DefaultFinishBulkAssembly()

       !------------------------------------------------------------------------------
       !     Neumann boundary conditions
       !------------------------------------------------------------------------------
       DO t=Solver % Mesh % NumberOfBulkElements + 1, &
           Solver % Mesh % NumberOfBulkElements + &
           Solver % Mesh % NumberOfBoundaryElements

         CurrentElement => Solver % Mesh % Elements(t)

         DO i=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN

             !------------------------------------------------------------------------------
             !             Set the current element pointer in the model structure to
             !             reflect the element being processed
             !------------------------------------------------------------------------------
             Model % CurrentElement => CurrentElement
             !------------------------------------------------------------------------------
             n = CurrentElement % TYPE % NumberOfNodes
             NodeIndexes => CurrentElement % NodeIndexes
             IF ( ANY( PotentialPerm(NodeIndexes) <= 0 ) ) CYCLE

             FluxBC = ListGetLogical(Model % BCs(i) % Values, &
                 'Current Density BC',gotIt) 
             IF(GotIt .AND. .NOT. FluxBC) CYCLE

             !------------------------------------------------------------------------------
             !             BC: cond@Phi/@n = g
             !------------------------------------------------------------------------------
             Load = 0.0d0
             Load(1:n) = ListGetReal( Model % BCs(i) % Values,'Current Density', &
                 n,NodeIndexes,gotIt )
             IF(.NOT. GotIt) CYCLE

             ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
             ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
             ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

             !------------------------------------------------------------------------------
             !             Get element matrix and rhs due to boundary conditions ...
             !------------------------------------------------------------------------------
             CALL StatCurrentBoundary( LocalStiffMatrix, LocalForce,  &
                 Load, CurrentElement, n, ElementNodes )
             !------------------------------------------------------------------------------
             !             Update global matrices from local matrices
             !------------------------------------------------------------------------------

	     CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

      !------------------------------------------------------------------------------
           END IF ! of currentelement bc == bcs(i)
         END DO ! of i=1,model bcs
       END DO   ! Neumann BCs
       !------------------------------------------------------------------------------

       !------------------------------------------------------------------------------
       !    FinishAssembly must be called after all other assembly steps, but before
       !    Dirichlet boundary settings. Actually no need to call it except for
       !    transient simulations.
       !------------------------------------------------------------------------------
       CALL DefaultFinishAssembly()

       !------------------------------------------------------------------------------
       !    Dirichlet boundary conditions
       !------------------------------------------------------------------------------
       CALL DefaultDirichletBCs()

       at = CPUTime() - at
       WRITE( Message, * ) 'Assembly (s)          :',at
       CALL Info( 'StatCurrentSolve', Message, Level=5 )
       !------------------------------------------------------------------------------
       !    Solve the system and we are done.
       !------------------------------------------------------------------------------
       st = CPUTime()
       Norm = DefaultSolve()

       st = CPUTime() - st
       WRITE( Message, * ) 'Solve (s)             :',st
       CALL Info( 'StatCurrentSolve', Message, Level=5 )


!------------------------------------------------------------------------------
!    Compute the electric field from the potential: E = -grad Phi
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Compute the volume current: J = cond (-grad Phi)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Compute the Joule heating: H,tot = Integral (E . D)dV
!------------------------------------------------------------------------------
       
       IF ( Control .OR. CalculateCurrent .OR. CalculateHeating .OR. &
           CalculateNodalHeating ) THEN 
         CALL GeneralCurrent( Model, Potential, PotentialPerm )

         WRITE( Message, * ) 'Total Heating Power   :', Heatingtot
         CALL Info( 'StatCurrentSolve', Message, Level=4 )
         CALL ListAddConstReal( Model % Simulation, &
             'RES: Total Joule Heating', Heatingtot )
         
         PotDiff = DirichletDofsRange( Solver )
         
         IF( PotDiff > 0 ) THEN
           Resistance = PotDiff**2 / HeatingTot
           WRITE( Message, * ) 'Effective Resistance  :', Resistance
           CALL Info( 'StatCurrentSolve', Message, Level=4 )
           CALL ListAddConstReal( Model % Simulation, &
               'RES: Effective Resistance', Resistance )
         END IF
       END IF

       IF(Control ) THEN
         WRITE( Message, * ) 'Total Volume          :', VolTot
         CALL Info( 'StatCurrentSolve', Message, Level=4 )

         ControlScaling = 1.0_dp
         IF( ControlPower ) THEN
           ControlScaling = SQRT( ControlTarget / HeatingTot )
         ELSE IF( ControlCurrent ) THEN
           IF( PotDiff > 0.0d0 ) THEN
             CurrentTot = HeatingTot / PotDiff
             ControlScaling = ControlTarget / CurrentTot
             WRITE( Message, * ) 'Total Current         :', CurrentTot
             CALL Info( 'StatCurrentSolve', Message, Level=4 )
             CALL ListAddConstReal( Model % Simulation, &
                 'RES: TotalCurrent', CurrentTot )
           ELSE
             CALL Warn('StatCurrentSolver','Current cannot be determined without pot. difference')
           END IF
         END IF

         WRITE( Message, * ) 'Control Scaling       :', ControlScaling
         CALL Info( 'StatCurrentSolve', Message, Level=4 )
         CALL ListAddConstReal( Model % Simulation, &
             'RES: CurrentSolver Scaling', ControlScaling )
         Potential = ControlScaling * Potential
!         Solver % Variable % Norm = ControlScaling * Solver % Variable % Norm

         IF ( CalculateHeating )  Heating = ControlScaling**2 * Heating
         IF ( CalculateNodalHeating)  &
             NodalHeating = ControlScaling**2 * NodalHeating
         IF ( CalculateCurrent )  VolCurrent = ControlScaling * VolCurrent
       END IF

       IF( Solver % Variable % NonlinConverged > 0 ) EXIT

     END DO

     IF (ListGetLogical(Params, 'Adaptive Mesh Refinement', GotIt)) THEN
       IF(.NOT. ListGetLogical(Params,'Library Adaptivity',GotIt)) THEN
         CALL RefineMesh(Model, Solver, Solver % Variable % Values, Solver % Variable % Perm, &
             StatCurrentSolver_Inside_Residual, StatCurrentSolver_Edge_Residual, &
             StatCurrentSolver_Boundary_Residual)
       END IF
     END IF

     
    CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Potential')
   
    IF ( CalculateCurrent ) THEN
      CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Volume Current')
    END IF
    
    IF ( CalculateHeating ) THEN
      CALL InvalidateVariable( Model % Meshes, Solver % Mesh, 'Joule Heating')
    END IF

    IF ( CalculateNodalHeating ) THEN
      CALL InvalidateVariable( Model % Meshes, Solver % Mesh, &
          'Nodal Joule Heating')
    END IF

    CALL DefaultFinish()
    

!------------------------------------------------------------------------------
 
   CONTAINS

!------------------------------------------------------------------------------
!> Compute the Current and Joule Heating at model nodes.
!------------------------------------------------------------------------------
  SUBROUTINE GeneralCurrent( Model, Potential, Reorder )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: Potential(:)
    INTEGER :: Reorder(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes 
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    REAL(KIND=dp), POINTER :: U_Integ(:), V_Integ(:), W_Integ(:), S_Integ(:)
    REAL(KIND=dp), ALLOCATABLE :: SumOfWeights(:), tmp(:)
    REAL(KIND=dp) :: Conductivity(3,3,Model % MaxElementNodes)
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes)
    REAL(KIND=dp) :: dBasisdx(Model % MaxElementNodes,3)
    REAL(KIND=DP) :: SqrtElementMetric, ElemVol
    REAL(KIND=dp) :: ElementPot(Model % MaxElementNodes)
    REAL(KIND=dp) :: Current(3)
    REAL(KIND=dp) :: s, ug, vg, wg, Grad(3), EpsGrad(3)
    REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
    REAL(KIND=dp) :: HeatingDensity, x, y, z
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: N_Integ, t, tg, i, j, k
    LOGICAL :: Stat

!------------------------------------------------------------------------------

    ALLOCATE( Nodes % x( Model % MaxElementNodes ) )
    ALLOCATE( Nodes % y( Model % MaxElementNodes ) )
    ALLOCATE( Nodes % z( Model % MaxElementNodes ) )

    IF( CalculateHeating .OR. CalculateCurrent ) THEN
      ALLOCATE( SumOfWeights( Model % NumberOfNodes ) )
      SumOfWeights = 0.0d0
    END IF

    HeatingTot = 0.0d0
    VolTot = 0.0d0
    IF ( CalculateHeating )  Heating = 0.0d0
    IF ( CalculateNodalHeating)  NodalHeating = 0.0d0
    IF ( CalculateCurrent )  VolCurrent = 0.0d0

    IF( GetCondAtIp ) THEN
      CALL ListInitElementKeyword( CondAtIp_h,'Material','Electric Conductivity')
    END IF
     
!------------------------------------------------------------------------------
!   Go through model elements, we will compute on average of elementwise
!   fluxes to nodes of the model
!------------------------------------------------------------------------------
    DO t = 1,Solver % NumberOfActiveElements
!------------------------------------------------------------------------------
!        Check if this element belongs to a body where electrostatics
!        should be calculated
!------------------------------------------------------------------------------
       Element => Solver % Mesh % Elements( Solver % ActiveElements( t ) )
       Model % CurrentElement => Element
       NodeIndexes => Element % NodeIndexes

       IF ( Element % PartIndex /= ParEnv % MyPE ) CYCLE

       n = Element % TYPE % NumberOfNodes

       IF ( ANY(Reorder(NodeIndexes) == 0) ) CYCLE

       ElementPot(1:n) = Potential( Reorder( NodeIndexes(1:n) ) )
       
       Nodes % x(1:n) = Model % Nodes % x( NodeIndexes )
       Nodes % y(1:n) = Model % Nodes % y( NodeIndexes )
       Nodes % z(1:n) = Model % Nodes % z( NodeIndexes )

!------------------------------------------------------------------------------
!    Gauss integration stuff
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
       U_Integ => IntegStuff % u
       V_Integ => IntegStuff % v
       W_Integ => IntegStuff % w
       S_Integ => IntegStuff % s
       N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------

       IF( .NOT. GetCondAtIp ) THEN
         k = ListGetInteger( Model % Bodies( Element % BodyId ) % &
             Values, 'Material', minv=1, maxv=Model % NumberOfMaterials )

         CALL ListGetRealArray( Model % Materials(k) % Values, &
             'Electric Conductivity', Cwrk, n, NodeIndexes, gotIt )

         Conductivity = 0.0d0
         IF ( SIZE(Cwrk,1) == 1 ) THEN
           DO i=1,3
             Conductivity( i,i,1:n ) = Cwrk( 1,1,1:n )
           END DO
         ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
           DO i=1,MIN(3,SIZE(Cwrk,1))
             Conductivity(i,i,1:n) = Cwrk(i,1,1:n)
           END DO
         ELSE
           DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
               Conductivity( i,j,1:n ) = Cwrk(i,j,1:n)
             END DO
           END DO
         END IF
       END IF
         
!------------------------------------------------------------------------------
! Loop over Gauss integration points
!------------------------------------------------------------------------------

       HeatingDensity = 0.0d0
       Current = 0.0d0
       ElemVol = 0.0d0


       DO tg=1,N_Integ

          ug = U_Integ(tg)
          vg = V_Integ(tg)
          wg = W_Integ(tg)

!------------------------------------------------------------------------------
! Need SqrtElementMetric and Basis at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element, Nodes,ug,vg,wg, &
               SqrtElementMetric,Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          s = SqrtElementMetric * S_Integ(tg)

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
            x = SUM( Nodes % x(1:n)*Basis(1:n) )
            y = SUM( Nodes % y(1:n)*Basis(1:n) )
            z = SUM( Nodes % z(1:n)*Basis(1:n) )
            
            CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
            s = s * SqrtMetric * 2 * PI
          END IF

!------------------------------------------------------------------------------

          EpsGrad = 0.0d0
          IF( GetCondAtIp ) THEN
            CondAtIp = ListGetElementReal( CondAtIp_h, Basis, Element, Stat, GaussPoint = tg )
            DO j = 1, DIM
              Grad(j) = SUM( dBasisdx(1:n,j) * ElementPot(1:n) )
            END DO
            EpsGrad(1:dim) = CondAtIp * Grad(1:dim)
          ELSE
            DO j = 1, DIM
              Grad(j) = SUM( dBasisdx(1:n,j) * ElementPot(1:n) )
              DO i = 1, DIM
                EpsGrad(j) = EpsGrad(j) + SUM( Conductivity(j,i,1:n) * &
                    Basis(1:n) ) * SUM( dBasisdx(1:n,i) * ElementPot(1:n) )
              END DO
            END DO
          END IF

            
          VolTot = VolTot + s

          HeatingTot = HeatingTot + &
               s * SUM( Grad(1:DIM) * EpsGrad(1:DIM) )

          IF( CalculateHeating .OR. CalculateCurrent .OR. CalculateNodalHeating ) THEN
            HeatingDensity = HeatingDensity + &
                s * SUM( Grad(1:DIM) * EpsGrad(1:DIM) ) 
            DO j = 1,DIM
              Current(j) = Current(j) - EpsGrad(j) * s
            END DO
            
            ElemVol = ElemVol + s
          END IF

       END DO! of the Gauss integration points

!------------------------------------------------------------------------------
!   Weight with element area if required
!------------------------------------------------------------------------------

       IF( CalculateHeating .OR. CalculateCurrent ) THEN
         IF ( ConstantWeights ) THEN
           HeatingDensity = HeatingDensity / ElemVol
           Current(1:Dim) = Current(1:Dim) / ElemVol
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + 1
         ELSE
           SumOfWeights( Reorder( NodeIndexes(1:n) ) ) = &
               SumOfWeights( Reorder( NodeIndexes(1:n) ) ) + ElemVol
         END IF
       END IF
         
       IF ( CalculateHeating ) THEN
         Heating( Reorder(NodeIndexes(1:n)) ) = &
             Heating( Reorder(NodeIndexes(1:n)) ) + HeatingDensity
       END IF
       
       IF ( CalculateNodalHeating ) THEN
         NodalHeating( Reorder(NodeIndexes(1:n)) ) = &
             NodalHeating( Reorder(NodeIndexes(1:n)) ) + HeatingDensity
       END IF
         
       IF ( CalculateCurrent ) THEN
         DO j=1,DIM 
           VolCurrent(DIM*(Reorder(NodeIndexes(1:n))-1)+j) = &
               VolCurrent(DIM*(Reorder(NodeIndexes(1:n))-1)+j) + &
               Current(j)
         END DO
       END IF

    END DO! of the bulk elements

    IF ( CalculateHeating .OR. CalculateCurrent) THEN
      IF ( ParEnv % PEs > 1) THEN
        VolTot     = ParallelReduction(VolTot)
        HeatingTot = ParallelReduction(HeatingTot)
        
        IF ( CalculateCurrent) THEN
          ALLOCATE(tmp(SIZE(VolCurrent)/dim))
          DO i=1,dim
            tmp = VolCurrent(i::dim)
            CALL ParallelSumVector(Solver % Matrix, tmp)
            Volcurrent(i::dim) = tmp
          END DO
        END IF
        IF (CalculateHeating ) CALL ParallelSumVector(Solver % Matrix, Heating)
        CALL ParallelSumVector(Solver % Matrix, SumOfWeights)
      END IF
      
!------------------------------------------------------------------------------
!   Finally, compute average of the fluxes at nodes
!------------------------------------------------------------------------------
      DO i = 1, Model % NumberOfNodes
        IF ( ABS( SumOfWeights(i) ) > 0.0D0 ) THEN
          IF ( CalculateHeating )  Heating(i) = Heating(i) / SumOfWeights(i)
          DO j = 1, DIM
            IF ( CalculateCurrent )  VolCurrent(DIM*(i-1)+j) = &
                VolCurrent(DIM*(i-1)+j) /  SumOfWeights(i)
          END DO
        END IF
      END DO
      DEALLOCATE( SumOfWeights ) 
    END IF
      
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

!------------------------------------------------------------------------------
   END SUBROUTINE GeneralCurrent
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
     SUBROUTINE StatCurrentCompose( StiffMatrix,Force,Conductivity, &
                            Load,Element,n,Nodes )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: StiffMatrix(:,:),Force(:),Load(:), Conductivity(:,:,:)
       INTEGER :: n
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
 
       REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,C(3,3),x,y,z
       LOGICAL :: Stat

       INTEGER :: i,p,q,t,DIM
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------
       DIM = CoordinateSystemDimension()

       Force = 0.0d0
       StiffMatrix = 0.0d0
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                    Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
           y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
           z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
         END IF

         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
 
         S = S * SqrtElementMetric * SqrtMetric

         L = SUM( Load(1:n) * Basis )

         IF( GetCondAtIp ) THEN
           CondAtIp = ListGetElementReal( CondAtIp_h, Basis, Element, Stat, GaussPoint = t )
           C(1:dim,1:dim) = 0.0_dp
           DO i=1,dim
             C(i,i) = CondAtIp
           END DO
         ELSE
           DO i=1,DIM
             DO j=1,DIM
               C(i,j) = SUM( Conductivity(i,j,1:n) * Basis(1:n) )
             END DO
           END DO
         END IF
         
!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------
         DO p=1,N
           DO q=1,N
             A = 0.d0
             DO i=1,DIM
               DO J=1,DIM
                 A = A + C(i,j) * dBasisdx(p,i) * dBasisdx(q,j)
               END DO
             END DO
             StiffMatrix(p,q) = StiffMatrix(p,q) + S*A
           END DO
           Force(p) = Force(p) + S*L*Basis(p)
         END DO
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE StatCurrentCompose
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for boundary conditions
!>  of the electrostatic equation. 
!------------------------------------------------------------------------------
   SUBROUTINE StatCurrentBoundary( BoundaryMatrix, BoundaryVector, &
        LoadVector, Element, n, Nodes )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: BoundaryMatrix(:,:), BoundaryVector(:), LoadVector(:)
     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),SqrtElementMetric
     REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)

     REAL(KIND=dp) :: u,v,w,s,x,y,z
     REAL(KIND=dp) :: Force
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: t,q,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0d0
     BoundaryMatrix = 0.0d0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivates at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                 Basis,dBasisdx )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
         IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
           x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
           y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
           z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
         END IF

         CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
 
         s = S_Integ(t) * SqrtElementMetric * SqrtMetric

!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis )

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE StatCurrentBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE StatCurrentSolver
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  FUNCTION StatCurrentSolver_boundary_residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
  !------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
  !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Edge
  !------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes, EdgeNodes
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i, j, k, n, l, t, dim, Pn, En, nd
    LOGICAL :: stat, Found
    INTEGER, ALLOCATABLE :: Indexes(:)
    REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
    REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), &
        EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
        dBasisdx(:, :), Potential(:), Flux(:)
    REAL(KIND=dp) :: Normal(3), EdgeLength, gx, gy, gz, Conductivity
    REAL(KIND=dp) :: u, v, w, s, detJ
    REAL(KIND=dp) :: Residual, ResidualNorm
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    LOGICAL :: First = .TRUE., Dirichlet
    SAVE Hwrk, First
  !------------------------------------------------------------------------------
  
  !    Initialize:
  !    -----------
    IF (First) THEN
      First = .FALSE.
      NULLIFY (Hwrk)
    END IF
  
    Indicator = 0.0d0
    Gnorm = 0.0d0
  
    Metric = 0.0d0
    DO i = 1, 3
      Metric(i, i) = 1.0d0
    END DO
  
    SELECT CASE (CurrentCoordinateSystem())
    CASE (AxisSymmetric, CylindricSymmetric)
      dim = 3
    CASE DEFAULT
      dim = CoordinateSystemDimension()
    END SELECT
  !
  !    ---------------------------------------------
  
    Element => Edge % BoundaryInfo % Left
  
    IF (.NOT. ASSOCIATED(Element)) THEN
      Element => Edge % BoundaryInfo % Right
    ELSE IF (ANY(Perm(Element % NodeIndexes) <= 0)) THEN
      Element => Edge % BoundaryInfo % Right
    END IF
  
    IF (.NOT. ASSOCIATED(Element)) RETURN
    IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN
  
    en = Edge % TYPE % NumberOfNodes
    pn = Element % TYPE % NumberOfNodes
  
    ALLOCATE (EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en))
  
    EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
    EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
    EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)
  
    nd = GetElementNOFDOFs(Element)
    ALLOCATE (Potential(nd), Basis(nd), &
              x(en), y(en), z(en), EdgeBasis(nd), &
              dBasisdx(nd, 3), NodalConductivity(nd), Flux(nd), &
              Indexes(nd))
  
    nd = GetElementDOFs(Indexes, Element)
  
    ALLOCATE (Nodes % x(nd), Nodes % y(nd), Nodes % z(nd))
    Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
    Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
    Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))
  
    DO l = 1, en
      DO k = 1, pn
        IF (Edge % NodeIndexes(l) == Element % NodeIndexes(k)) THEN
          x(l) = Element % TYPE % NodeU(k)
          y(l) = Element % TYPE % NodeV(k)
          z(l) = Element % TYPE % NodeW(k)
          EXIT
        END IF
      END DO
    END DO
  !
  !    Integrate square of residual over boundary element:
  !    ---------------------------------------------------
  
    Indicator = 0.0d0
    EdgeLength = 0.0d0
    ResidualNorm = 0.0d0
  
    DO j = 1, Model % NumberOfBCs
      IF (Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag) CYCLE
  
  !
  !       Check if dirichlet BC given:
  !       ----------------------------
      Dirichlet = ListCheckPresent(Model % BCs(j) % Values, &
                                    ComponentName(Model % Solver % Variable))
      IF (.NOT. Dirichlet) THEN
        Dirichlet = ListCheckPrefix(Model % BCs(j) % Values, &
                                    'Constraint Mode')
      END IF
      ! TODO s = ListGetConstReal( Model % BCs(j) % Values,'Potential',Dirichlet )
  
  !       Get various flux bc options:
  !       ----------------------------
  
  !       ...given flux:
  !       --------------
      Flux(1:en) = ListGetReal(Model % BCs(j) % Values, &
                                'Electric Flux', en, Edge % NodeIndexes, Found)
  
  !       get material parameters:
  !       ------------------------
      k = ListGetInteger(Model % Bodies(Element % BodyId) % Values, 'Material', &
                          minv=1, maxv=Model % NumberOFMaterials)
  
      CALL ListGetRealArray(Model % Materials(k) % Values, &
                            'Electric Conductivity', Hwrk, en, Edge % NodeIndexes, stat)
      IF (.NOT. stat) THEN
        CALL Fatal('StatCurrentSolver_boundary_residual:','Electric Conductivity not found')
      END IF
      NodalConductivity(1:en) = Hwrk(1, 1, 1:en)
  
  !       elementwise nodal solution:
  !       ---------------------------
      nd = GetElementDOFs(Indexes, Element)
      Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
  
  !       do the integration:
  !       -------------------
      EdgeLength = 0.0d0
      ResidualNorm = 0.0d0
  
      IntegStuff = GaussPoints(Edge)
  
      DO t = 1, IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)
  
        stat = ElementInfo(Edge, EdgeNodes, u, v, w, detJ, &
                            EdgeBasis, dBasisdx)
        Normal = NormalVector(Edge, EdgeNodes, u, v, .TRUE.)
  
        IF (CurrentCoordinateSystem() == Cartesian) THEN
          s = IntegStuff % s(t) * detJ
        ELSE
          gx = SUM(EdgeBasis(1:en) * EdgeNodes % x(1:en))
          gy = SUM(EdgeBasis(1:en) * EdgeNodes % y(1:en))
          gz = SUM(EdgeBasis(1:en) * EdgeNodes % z(1:en))
          CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                    Symb, dSymb, gx, gy, gz)
          s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF
  
  !
  !          Integration point in parent element local
  !          coordinates:
  !          -----------------------------------------
        u = SUM(EdgeBasis(1:en) * x(1:en))
        v = SUM(EdgeBasis(1:en) * y(1:en))
        w = SUM(EdgeBasis(1:en) * z(1:en))
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
  
  !
  !          Conductivity at the integration point:
  !          --------------------------------------
        Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
  !
  !          given flux at integration point:
  !          --------------------------------
        Residual = -SUM(Flux(1:en) * EdgeBasis(1:en))
  
  !          flux given by the computed solution, and
  !          force norm for scaling the residual:
  !          -----------------------------------------
        IF (CurrentCoordinateSystem() == Cartesian) THEN
          DO k = 1, dim
            Residual = Residual + Conductivity * &
                        SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k)
  
            Gnorm = Gnorm + s * (Conductivity * &
                                  SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k))**2
          END DO
        ELSE
          DO k = 1, dim
            DO l = 1, dim
              Residual = Residual + Metric(k, l) * Conductivity * &
                          SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(l)
  
              Gnorm = Gnorm + s * (Metric(k, l) * Conductivity * &
                                    SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(l))**2
            END DO
          END DO
        END IF
  
        EdgeLength = EdgeLength + s
        IF (.NOT. Dirichlet) THEN
          ResidualNorm = ResidualNorm + s * Residual**2
        END IF
      END DO
      EXIT
    END DO
  
    IF (CoordinateSystemDimension() == 3) EdgeLength = SQRT(EdgeLength)
  
  !    Gnorm = EdgeLength * Gnorm
    Indicator = EdgeLength * ResidualNorm
  !------------------------------------------------------------------------------
    END FUNCTION StatCurrentSolver_boundary_residual
    !------------------------------------------------------------------------------
      
  
  FUNCTION StatCurrentSolver_edge_residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
  !------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
  
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Edge
  !------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes, EdgeNodes
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i, j, k, l, n, t, dim, En, Pn, nd
    INTEGER, ALLOCATABLE :: Indexes(:)
    LOGICAL :: stat
    REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
    REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
    REAL(KIND=dp) :: Conductivity
    REAL(KIND=dp) :: u, v, w, s, detJ
    REAL(KIND=dp) :: Grad(3, 3), Normal(3), EdgeLength, Jump
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), EdgeBasis(:), &
        Basis(:), dBasisdx(:, :), Potential(:)
    REAL(KIND=dp) :: ResidualNorm
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(ValueList_t), POINTER :: Material
  
    LOGICAL :: First = .TRUE.
    SAVE Hwrk, First
  !------------------------------------------------------------------------------
  
    !    Initialize:
    !    -----------
    IF (First) THEN
      First = .FALSE.
      NULLIFY (Hwrk)
    END IF
  
    SELECT CASE (CurrentCoordinateSystem())
    CASE (AxisSymmetric, CylindricSymmetric)
      dim = 3
    CASE DEFAULT
      dim = CoordinateSystemDimension()
    END SELECT
  
    Metric = 0.0d0
    DO i = 1, 3
      Metric(i, i) = 1.0d0
    END DO
    Grad = 0.0d0
  !
  !    ---------------------------------------------
  
    n = Mesh % MaxElementDOFs
    ALLOCATE (Nodes % x(n), Nodes % y(n), Nodes % z(n))
  
    en = Edge % TYPE % NumberOfNodes
    ALLOCATE (EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en))
  
    EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
    EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
    EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)
  
    ALLOCATE (NodalConductivity(en), EdgeBasis(en), Basis(n), &
              dBasisdx(n, 3), x(en), y(en), z(en), Potential(n), Indexes(n))
  
  !    Integrate square of jump over edge:
  !    -----------------------------------
    ResidualNorm = 0.0d0
    EdgeLength = 0.0d0
    Indicator = 0.0d0
  
    IntegStuff = GaussPoints(Edge)
  
    DO t = 1, IntegStuff % n
  
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
  
      stat = ElementInfo(Edge, EdgeNodes, u, v, w, detJ, &
                          EdgeBasis, dBasisdx)
  
      Normal = NormalVector(Edge, EdgeNodes, u, v, .FALSE.)
  
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        s = IntegStuff % s(t) * detJ
      ELSE
        u = SUM(EdgeBasis(1:en) * EdgeNodes % x(1:en))
        v = SUM(EdgeBasis(1:en) * EdgeNodes % y(1:en))
        w = SUM(EdgeBasis(1:en) * EdgeNodes % z(1:en))
  
        CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                  Symb, dSymb, u, v, w)
        s = IntegStuff % s(t) * detJ * SqrtMetric
      END IF
  
      !
      ! Compute flux over the edge as seen by elements
      ! on both sides of the edge:
      ! ----------------------------------------------
      DO i = 1, 2
        SELECT CASE (i)
        CASE (1)
          Element => Edge % BoundaryInfo % Left
        CASE (2)
          Element => Edge % BoundaryInfo % Right
        END SELECT
  !
  !          Can this really happen (maybe it can...)  ?
  !          -------------------------------------------
        IF (.NOT. ASSOCIATED(Element)) CYCLE
        IF (ANY(Perm(Element % NodeIndexes) <= 0)) CYCLE
  !
  !          Next, get the integration point in parent
  !          local coordinates:
  !          -----------------------------------------
        pn = Element % TYPE % NumberOfNodes
  
        DO j = 1, en
          DO k = 1, pn
            IF (Edge % NodeIndexes(j) == Element % NodeIndexes(k)) THEN
              x(j) = Element % TYPE % NodeU(k)
              y(j) = Element % TYPE % NodeV(k)
              z(j) = Element % TYPE % NodeW(k)
              EXIT
            END IF
          END DO
        END DO
  
        u = SUM(EdgeBasis(1:en) * x(1:en))
        v = SUM(EdgeBasis(1:en) * y(1:en))
        w = SUM(EdgeBasis(1:en) * z(1:en))
  !
  !          Get parent element basis & derivatives at the integration point:
  !          -----------------------------------------------------------------
        nd = GetElementDOFs(Indexes, Element)
        Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
        Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
        Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))
  
        stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
  !
  !          Material parameters:
  !          --------------------
        k = ListGetInteger(Model % Bodies( &
                            Element % BodyId) % Values, 'Material', &
                            minv=1, maxv=Model % NumberOFMaterials)
  
        Material => Model % Materials(k) % Values
        CALL ListGetRealArray(Material, &
                              'Electric Conductivity', Hwrk, en, Edge % NodeIndexes, stat)
        IF (.NOT. stat) THEN
          CALL Fatal('StatCurrentSolver_edge_residual:', 'Electric Conductivity not found')
        END IF
  
        NodalConductivity(1:en) = Hwrk(1, 1, 1:en)
        Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
  !
  !          Potential at element nodal points:
  !          ------------------------------------
        Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
  !
  !          Finally, the flux:
  !          ------------------
        DO j = 1, dim
          Grad(j, i) = Conductivity * SUM(dBasisdx(1:nd, j) * Potential(1:nd))
        END DO
      END DO
  
  !       Compute square of the flux jump:
  !       -------------------------------
      EdgeLength = EdgeLength + s
      Jump = 0.0d0
      DO k = 1, dim
        IF (CurrentCoordinateSystem() == Cartesian) THEN
          Jump = Jump + (Grad(k, 1) - Grad(k, 2)) * Normal(k)
        ELSE
          DO l = 1, dim
            Jump = Jump + &
                    Metric(k, l) * (Grad(k, 1) - Grad(k, 2)) * Normal(l)
          END DO
        END IF
      END DO
      ResidualNorm = ResidualNorm + s * Jump**2
    END DO
  
    IF (dim == 3) EdgeLength = SQRT(EdgeLength)
    Indicator = EdgeLength * ResidualNorm
  
    DEALLOCATE (Nodes % x, Nodes % y, Nodes % z)
  
    DEALLOCATE (EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
    DEALLOCATE (x, y, z, NodalConductivity, EdgeBasis, &
                Basis, dBasisdx, Potential)
  !------------------------------------------------------------------------------
    END FUNCTION StatCurrentSolver_edge_residual
    !------------------------------------------------------------------------------
      
  !------------------------------------------------------------------------------
  FUNCTION StatCurrentSolver_inside_residual(Model, Element, Mesh, Quant, Perm, Fnorm) RESULT(Indicator)
  !------------------------------------------------------------------------------
    USE DefUtils
  !------------------------------------------------------------------------------
    IMPLICIT NONE
  !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Perm(:)
    REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
  !------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    INTEGER :: i, j, k, l, n, t, dim, nd
    INTEGER, ALLOCATABLE :: Indexes(:)
    LOGICAL :: stat, Found
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
    REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
    REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Potential(:), PrevPot(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:, :), ddBasisddx(:, :, :)
    REAL(KIND=dp) :: u, v, w, s, detJ
    REAL(KIND=dp) :: Conductivity, dt
    REAL(KIND=dp) :: Residual, ResidualNorm, Area
    TYPE(ValueList_t), POINTER :: Material
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  
    LOGICAL :: First = .TRUE.
    SAVE Hwrk, First
  !------------------------------------------------------------------------------
  
  !    Initialize:
  !    -----------
    Indicator = 0.0d0
    Fnorm = 0.0d0
  !
  !    Check if this eq. computed in this element:
  !    -------------------------------------------
    IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN
  
    IF (First) THEN
      First = .FALSE.
      NULLIFY (Hwrk)
    END IF
  
    Metric = 0.0d0
    DO i = 1, 3
      Metric(i, i) = 1.0d0
    END DO
  
    SELECT CASE (CurrentCoordinateSystem())
    CASE (AxisSymmetric, CylindricSymmetric)
      dim = 3
    CASE DEFAULT
      dim = CoordinateSystemDimension()
    END SELECT
  
  !    Allocate local arrays
  !    ----------------------
    nd = GetElementNOFDOFs(Element)
    n = GetElementNOFNodes(Element)
    ALLOCATE (NodalConductivity(nd), &
              PrevPot(nd), NodalSource(nd), Potential(nd), &
              Basis(nd), dBasisdx(nd, 3), ddBasisddx(nd, 3, 3), Indexes(nd))
  !
  !    Element nodal points:
  !    ---------------------
    ALLOCATE (Nodes % x(nd), Nodes % y(nd), Nodes % z(nd))
  
    nd = GetElementDOFs(Indexes, Element)
    Nodes % x = Mesh % Nodes % x(Indexes(1:nd))
    Nodes % y = Mesh % Nodes % y(Indexes(1:nd))
    Nodes % z = Mesh % Nodes % z(Indexes(1:nd))
  !
  !    Elementwise nodal solution:
  !    ---------------------------
    Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
  !
  !    Check for time dep.
  !    -------------------
    PrevPot(1:nd) = Potential(1:nd)
    dt = Model % Solver % dt
    IF (ListGetString(Model % Simulation, 'Simulation Type') == 'transient') THEN
      Var => VariableGet(Model % Variables, 'Potential', .TRUE.)
      PrevPot(1:nd) = Var % PrevValues(Var % Perm(Indexes(1:nd)), 1)
    END IF
  !
  !    Material parameters: conductivity
  !    ---------------------------------
    k = ListGetInteger(Model % Bodies(Element % BodyId) % Values, 'Material', &
                        minv=1, maxv=Model % NumberOfMaterials)
  
    Material => Model % Materials(k) % Values
  
    CALL ListGetRealArray(Material, 'Electric Conductivity', Hwrk, n, Element % NodeIndexes, stat)
    IF (.NOT. stat) THEN
      CALL Fatal('StatCurrentSolver_inside_residual:', 'Electric Conductivity not found')
    END IF
    NodalConductivity(1:n) = Hwrk(1, 1, 1:n)
  
  !
  !    Current source density (source):
  !    --------------------------------
  !
    k = ListGetInteger( &
        Model % Bodies(Element % BodyId) % Values, 'Body Force', Found, &
        1, Model % NumberOFBodyForces)
  
    NodalSource = 0.0d0
    IF (Found .AND. k > 0) THEN
      NodalSource(1:n) = ListGetReal(Model % BodyForces(k) % Values, &
          'Current Source', n, Element % NodeIndexes, stat)
    END IF
  
  !
  !    Integrate square of residual over element:
  !    ------------------------------------------
  
    ResidualNorm = 0.0d0
    Area = 0.0d0
  
    IntegStuff = GaussPoints(Element)
    ddBasisddx = 0
  
    DO t = 1, IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)
  
      stat = ElementInfo(Element, Nodes, u, v, w, detJ, &
                          Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE.)
  
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        s = IntegStuff % s(t) * detJ
      ELSE
        u = SUM(Basis(1:nd) * Nodes % x(1:nd))
        v = SUM(Basis(1:nd) * Nodes % y(1:nd))
        w = SUM(Basis(1:nd) * Nodes % z(1:nd))
  
        CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                  Symb, dSymb, u, v, w)
        s = IntegStuff % s(t) * detJ * SqrtMetric
      END IF
  
      Conductivity = SUM(NodalConductivity(1:n) * Basis(1:n))
  !
  !       Residual of the current conservation equation:
  !
  !        R = -div(σ grad(u)) - s
  !       ---------------------------------------------------
  !
  !       or more generally:
  !
  !        R = -g^{jk} (σ u_{,j}}_{,k}) - s
  !       ---------------------------------------------------
  !
      Residual = -SUM(NodalSource(1:n) * Basis(1:n))
  
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        DO j = 1, dim
  !
  !             - grad(σ).grad(u):
  !             -------------------
          Residual = Residual - &
                      SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                      SUM(NodalConductivity(1:n) * dBasisdx(1:n, j))
  
  !
  !             - σ div(grad(u)):
  !             ------------------
          Residual = Residual - Conductivity * &
                      SUM(Potential(1:nd) * ddBasisddx(1:nd, j, j))
        END DO
      ELSE
        DO j = 1, dim
          DO k = 1, dim
  !
  !                - g^{jk} σ_{,k} u_{,j}:
  !                ------------------------
            Residual = Residual - Metric(j, k) * &
                        SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                        SUM(NodalConductivity(1:n) * dBasisdx(1:n, k))
  
  !
  !                - g^{jk} σ u_{,jk}:
  !                --------------------
            Residual = Residual - Metric(j, k) * Conductivity * &
                        SUM(Potential(1:nd) * ddBasisddx(1:nd, j, k))
  !
  !                + g^{jk} σ Γ_{jk}^l u_{,l}:
  !                ----------------------------
            DO l = 1, dim
              Residual = Residual + Metric(j, k) * Conductivity * &
                          Symb(j, k, l) * SUM(Potential(1:nd) * dBasisdx(1:nd, l))
            END DO
          END DO
        END DO
      END IF
  
  !
  !       Compute also force norm for scaling the residual:
  !       -------------------------------------------------
      DO i = 1, dim
        Fnorm = Fnorm + s * (SUM(NodalSource(1:n) * Basis(1:n)))**2
      END DO
      Area = Area + s
      ResidualNorm = ResidualNorm + s * Residual**2
    END DO
  
  !    Fnorm = Element % hk**2 * Fnorm
    Indicator = Element % hK**2 * ResidualNorm
  !------------------------------------------------------------------------------
  END FUNCTION StatCurrentSolver_inside_residual
!------------------------------------------------------------------------------
 
