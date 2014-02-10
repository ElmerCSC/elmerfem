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
! *  Authors: Mika Malinen, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Nov 2002
! *
! *****************************************************************************/
 
!------------------------------------------------------------------------------
!>  This is the wave equation solver for the convective transport equation. 
!>  As initial data one has to specify the field subject to the convection
!>  operator. This field should be declared using the sif-file flag Advection 
!>  Variable.  The rate of change of the field subject to the convection
!>  operator which is needed as initial condition at the beginning of time step
!>  is solved using RateOfChangeSolver and should be available to the solver
!>  (the name of this variable should be declared using the flag Rate Of Change
!>  Equation Variable).
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE TransportEquationSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Types
     USE Lists
     USE Integration
     USE ElementDescription
     USE SolverUtils
     USE ElementUtils

     IMPLICIT NONE

     TYPE VariablePtr_t
       TYPE(Variable_t), POINTER :: Var
     END TYPE VariablePtr_t

!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET:: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER  :: StiffMatrix
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: CurrentElement, Parent
     TYPE(Variable_t), POINTER :: FlowSol, Udot0Var 
     TYPE(VariablePtr_t), POINTER :: U0Var(:)
     TYPE(ValueList_t), POINTER :: Material

     INTEGER :: i,j,k,p,n,t,body_id,bf_id,istat,LocalNodes,&
         AdvectionVariableComponents, VelocityComponents
     INTEGER, SAVE :: ThisSolverCalls=0
     INTEGER, POINTER :: NodeIndexes(:), UPerm(:), U0Perm(:), &
         VelocityPerm(:)

     LOGICAL :: GotIt, stat, AllocationsDone = .FALSE.
     CHARACTER(LEN=MAX_NAME_LEN) :: AdvectionFlag, AdvectionVariable, &
         VariableName, EquationName
 
     REAL(KIND=dp), POINTER :: U(:), Udot0(:), ForceVector(:), Flow(:), &
         Velocity(:)
     REAL(KIND=dp) :: at,st,CPUTime, Norm, PrevNorm

     REAL(KIND=dp), ALLOCATABLE :: &
          LocalStiffMatrix(:,:),LocalForce(:), &
          LocalMassMatrix(:,:), LocalDampMatrix(:,:), &
          V1(:), V2(:), V3(:)

     SAVE LocalStiffMatrix, LocalForce, ElementNodes, & 
          AllocationsDone, LocalMassMatrix, LocalDampMatrix, &
          U0Var, V1, V2, V3      

!------------------------------------------------------------------------------
     ThisSolverCalls = ThisSolverCalls + 1

     IF (MOD(ThisSolverCalls,2)==1) THEN
       CALL Info('TransportEquationSolver', '')
       CALL Info('TransportEquationSolver', 'Starting initialization...')
       CALL Info('TransportEquationSolver', '')

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
       U     => Solver % Variable % Values
       UPerm => Solver % Variable % Perm
     
       LocalNodes = Model % NumberOfNodes
       StiffMatrix => Solver % Matrix
       ForceVector => StiffMatrix % RHS

       Norm = Solver % Variable % Norm

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
       AdvectionVariableComponents = Solver % Variable % DOFs

       IF ( .NOT. AllocationsDone ) THEN
         N = Model % MaxElementNodes
         
         ALLOCATE( ElementNodes % x( N ),                       &
             ElementNodes % y( N ),                             &
             ElementNodes % z( N ),                             &
             U0Var(AdvectionVariableComponents),                &
             LocalForce(N*AdvectionVariableComponents),         &
             LocalStiffMatrix(N*AdvectionVariableComponents,    &
             N*AdvectionVariableComponents),                    &
             LocalMassMatrix(N*AdvectionVariableComponents,     &
             N*AdvectionVariableComponents ),                   &
             LocalDampMatrix(N*AdvectionVariableComponents,     &
             N*AdvectionVariableComponents),                    &
             V1(N), V2(N), V3(N),                               &
             STAT=istat )
         
         IF ( istat /= 0 ) THEN
           CALL Fatal('TransportEquationSolver', &
               'Memory allocation error, aborting.')
         END IF

         AllocationsDone = .TRUE.
       END IF

!-----------------------------------------------------------------------------
!   Get initial conditions...
!-----------------------------------------------------------------------------
       AdvectionVariable = ListGetString(Solver % Values, &
           'Advection Variable', GotIt )

       IF (AdvectionVariableComponents == 1) THEN
         U0Var(1) % Var => VariableGet(Solver % Mesh % Variables, &
             AdvectionVariable)
         IF (ASSOCIATED( U0Var(1) % Var )) THEN
           IF (AdvectionVariableComponents /= U0Var(1) % Var % DOFs) &
               CALL Fatal('TransportEquationSolver', &
               'Solver Variable and Advection Variable Dofs do not equal')
         ELSE
           CALL Fatal('TransportEquationSolver', &
               'Specified Advection Variable does not exist')
         END IF
       ELSE
         DO i = 1, AdvectionVariableComponents
           WRITE(VariableName,'(A,A,I1)') TRIM(AdvectionVariable),' ', i
           U0Var(i) % Var => VariableGet(Solver % Mesh % Variables, &
               VariableName)
           IF (.NOT. ASSOCIATED( U0Var(i) % Var )) &
               CALL Fatal('TransportEquationSolver', &
               'Specified Advection Variable component does not exist')
         END DO
       END IF
       
       VariableName = ListGetString( Solver % Values, &
           'Rate Of Change Equation Variable' )
       Udot0Var => VariableGet(Solver % Mesh % Variables, VariableName)
       
!------------------------------------------------------------------------------
!    Get the type of Advection Velocity field and read the velocity field if 
!    computed by the Navier-Stokes equations solver
!------------------------------------------------------------------------------
       AdvectionFlag = ListGetString(Solver % Values, 'Advection')
       IF (AdvectionFlag == 'computed') THEN
         FlowSol => VariableGet( Solver % Mesh % Variables, &
             'Flow Solution' )
         IF ( ASSOCIATED( FlowSol ) ) THEN
           VelocityPerm => FlowSol % Perm
           Flow => FlowSol % Values  
           VelocityComponents = FlowSol % DOFs - 1 
           IF (VelocityComponents /= CoordinateSystemDimension()) &
               CALL Warn('TransportEquationSolver', &
               'Coordinate system and Advection Velocity dimensions unequal')
         ELSE
           CALL Fatal('TransportEquationSolver', &
               'Flow Solution does not exist')
         END IF
       ELSE
         IF (AdvectionFlag /= 'constant') CALL Fatal(     &
             'TransportEquationSolver', &
             'Advection flag should be either "computed" or "constant"') 
         VelocityComponents = CoordinateSystemDimension()         
       END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
       Solver % dt = 2.0d0 * dt
!------------------------------------------------------------------------------
!    Set initial conditions by updating (Solver % Variable)-fields
!------------------------------------------------------------------------------

       DO i=1,Model % NumberOfNodes 
         j = UPerm(i)
         k = U0Var(1) % Var % Perm(i)
         p = Udot0Var % Perm(i)
         IF(j>0 .AND. k>0 .AND. p>0) THEN
           DO t=1, AdvectionVariableComponents
             Solver % Variable % Values( &
                 (j-1)*AdvectionVariableComponents + t) = &
                 U0Var(t) % Var % Values(k)
             Solver % Variable % PrevValues( &
                 (j-1)*AdvectionVariableComponents + t,1) = &
                 Udot0Var % Values((p-1)*AdvectionVariableComponents + t)
             Solver % Variable % PrevValues( &
                 (i-1)*AdvectionVariableComponents + t,2) = 0.0d0
           END DO
         ELSE
           CALL Fatal('TransportEquationSolver',&
               'Nonmatching variable permutations') 
         END IF
       END DO

       at = CPUTime()

       CALL InitializeToZero( StiffMatrix, ForceVector )
       CALL InitializeTimestep( Solver )  ! Do substitution U(n-1) <- U(n)
                                          ! and read 'Bossak Alpha'

       EquationName = ListGetString( Solver % Values, 'Equation' )

!------------------------------------------------------------------------------
!    Do the assembly
!------------------------------------------------------------------------------
       DO t=1,Model % NumberOfBulkElements
!------------------------------------------------------------------------------
         CurrentElement => Solver % Mesh % Elements(t)
         IF ( .NOT. CheckElementEquation( Model, &
             CurrentElement, EquationName ) ) CYCLE

         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

         body_id = CurrentElement % Bodyid
         k = ListGetInteger( Model % Bodies(body_id) % Values,'Material' )
         Material => Model % Materials(k) % Values

         IF ( AdvectionFlag == 'constant' ) THEN
           SELECT CASE (VelocityComponents)
             CASE(1)
             V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                 NodeIndexes) 
             V2 = 0.0d0
             V3 = 0.0d0            
             CASE(2)
             V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                 NodeIndexes)
             V2(1:n) = ListGetReal( Material,'Advection Velocity 2',n, &
                 NodeIndexes)
             V3 = 0.0d0  
             CASE(3)
             V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                 NodeIndexes)
             V2(1:n) = ListGetReal( Material,'Advection Velocity 2',n, &
                 NodeIndexes)
             V3(1:n) = ListGetReal( Material,'Advection Velocity 3',n, &
                 NodeIndexes)
           END SELECT
         ELSE ! Use computed velocity field   
           DO i=1,n
             j = VelocityPerm( NodeIndexes(i) )
             SELECT CASE (VelocityComponents)
               CASE(1)
               V1(i) = Flow((j-1)*2+1)
               V2(i) = 0.0d0
               V3(i) = 0.0d0
               CASE (2)
               V1(i) = Flow((j-1)*3+1)
               V2(i) = Flow((j-1)*3+2)
               V3(i) = 0.0d0
               CASE (3)
               V1(i) = Flow((j-1)*4+1)
               V2(i) = Flow((j-1)*4+2)
               V3(i) = Flow((j-1)*4+3)
             END SELECT
           END DO
         END IF

!------------------------------------------------------------------------------
!      Get element local matrix, and rhs vector
!------------------------------------------------------------------------------
         CALL LocalMatrix(  LocalStiffMatrix, LocalDampMatrix, &
             LocalMassMatrix, LocalForce, CurrentElement, n, &
             AdvectionVariableComponents, ElementNodes, V1, V2, V3 )

         CALL Add2ndOrderTime2( LocalMassMatrix, LocalDampMatrix, &
             LocalStiffMatrix, LocalForce, 2.0d0*dt, N, &
             AdvectionVariableComponents, UPerm( NodeIndexes ), Solver )
       
!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
         CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
             ForceVector, LocalForce, n, AdvectionVariableComponents, &
             UPerm(NodeIndexes) ) 
       END DO

!------------------------------------------------------------------------------
!      Compute the matrix corresponding to the integral over boundaries
!------------------------------------------------------------------------------
       DO t = Solver % Mesh % NumberOfBulkElements + 1,  &
           Solver % Mesh % NumberOfBulkElements +  &
           Solver % Mesh % NumberOfBoundaryElements
         
         CurrentElement => Solver % Mesh % Elements(t)
         NodeIndexes => CurrentElement % NodeIndexes

         IF ( SIZE(NodeIndexes) > 1 ) THEN         

           Parent => CurrentELement % BoundaryInfo % Left
           stat = ASSOCIATED( Parent )
           IF (stat) stat = ALL(UPerm(Parent % NodeIndexes) > 0)
           IF ( .NOT. stat) THEN
             Parent => CurrentELement % BoundaryInfo % Right
             stat = ASSOCIATED( Parent ) 
             IF (stat) stat = ALL(UPerm(Parent % NodeIndexes) > 0)            
             IF ( .NOT. stat )  CALL Fatal( 'TransportEquationSolver', &
                 'No parent element can be found for given boundary element' )
           END IF
           IF ( .NOT. CheckElementEquation( Model, &
               Parent, EquationName ) ) CYCLE

           body_id = Parent % Bodyid
           k = ListGetInteger( Model % Bodies(body_id) % Values,'Material' )
           Material => Model % Materials(k) % Values

           n = CurrentElement % TYPE % NumberOfNodes
           ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

           IF ( AdvectionFlag == 'constant' ) THEN
             SELECT CASE (VelocityComponents)
               CASE(1)
               V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                   NodeIndexes) 
               V2 = 0.0d0
               V3 = 0.0d0            
               CASE(2)
               V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                   NodeIndexes)
               V2(1:n) = ListGetReal( Material,'Advection Velocity 2',n, &
                   NodeIndexes)
               V3 = 0.0d0  
               CASE(3)
               V1(1:n) = ListGetReal( Material,'Advection Velocity 1',n, &
                   NodeIndexes)
               V2(1:n) = ListGetReal( Material,'Advection Velocity 2',n, &
                   NodeIndexes)
               V3(1:n) = ListGetReal( Material,'Advection Velocity 3',n, &
                   NodeIndexes)
             END SELECT
           ELSE ! Use computed velocity field   
             DO i=1,n
               j = VelocityPerm( NodeIndexes(i) )
               SELECT CASE (VelocityComponents)
                 CASE(1)
                 V1(i) = Flow((j-1)*2+1)
                 V2(i) = 0.0d0
                 V3(i) = 0.0d0
                 CASE (2)
                 V1(i) = Flow((j-1)*3+1)
                 V2(i) = Flow((j-1)*3+2)
                 V3(i) = 0.0d0
                 CASE (3)
                 V1(i) = Flow((j-1)*4+1)
                 V2(i) = Flow((j-1)*4+2)
                 V3(i) = Flow((j-1)*4+3)
               END SELECT
             END DO
           END IF

           CALL LocalMatrixBoundary(  LocalStiffMatrix, LocalDampMatrix, &
               LocalMassMatrix, LocalForce, CurrentElement, n, &
               AdvectionVariableComponents, ElementNodes, V1, V2, V3 )

           CALL Add2ndOrderTime2( LocalMassMatrix, LocalDampMatrix, &
               LocalStiffMatrix, LocalForce, 2.0d0*dt, N, &
               AdvectionVariableComponents, UPerm( NodeIndexes ), Solver )

           CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
               ForceVector, LocalForce, n, AdvectionVariableComponents, &
               UPerm(NodeIndexes) ) 
         END IF
       END DO

!------------------------------------------------------------------------------
!    FinishAssembly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
      CALL FinishAssembly( Solver, ForceVector )

!------------------------------------------------------------------------------
!    Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------
      at = CPUTime() - at
!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
      st = CPUTime()
!------------------------------------------------------------------------------
!    Solve the primary unknown and update velocity and acceleration 
!    vectors
!------------------------------------------------------------------------------
      Norm = DefaultSolve()

      st = CPUTime() - st
      WRITE(Message,'(a,F8.2)') ' Assembly: (s)', at
      CALL Info( 'TransportEquationSolver', Message )
      WRITE(Message,'(a,F8.2)') ' Solve:    (s)', st
      CALL Info( 'TransportEquationSolver', Message )

      Solver % dt = dt

!------------------------------------------------------------------------------
!    Retrieve correct boundary conditions on the outflow boundary and 
!    update Advection Variable for the next time step. 
!------------------------------------------------------------------------------
      SELECT CASE(AdvectionVariableComponents)
        CASE(1)
        CALL ModifyBoundaryValues( Model, AdvectionVariable, 1, &
            AdvectionVariableComponents, U, UPerm)
      CASE DEFAULT
        DO i = 1, AdvectionVariableComponents
          WRITE(VariableName,'(A,A,I1)') TRIM(AdvectionVariable),' ', i
          CALL ModifyBoundaryValues( Model, VariableName, &
              i, AdvectionVariableComponents, U, UPerm)
        END DO
      END SELECT

      DO i=1,Model % NumberOfNodes 
        j = UPerm(i)
        k = U0Var(1) % Var % Perm(i)
        IF(j>0 .AND. k>0) THEN
          DO t = 1,AdvectionVariableComponents
            U0Var(t) % Var % Values(k) = U((j-1)*AdvectionVariableComponents+t)
          END DO
        ELSE
          CALL Fatal('TransportEquationSolver',&
              'Nonmatching variable permutations') 
        END IF
      END DO
        
!-----------------------------------------------------------------------------
    END IF         ! IF (MOD(ThisSolverCalls,2)==1)
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------- 
  CONTAINS
!------------------------------------------------------------------------------
    SUBROUTINE LocalMatrix( StiffMatrix, DampMatrix, MassMatrix, &
        Force, Element, n, VariableComponents, Nodes, V1, V2, V3 )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), &
          MassMatrix(:,:), DampMatrix(:,:),           &
          V1(:), V2(:), V3(:)
      INTEGER :: n, VariableComponents
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
      REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,temp1,temp2
      REAL(KIND=dp), DIMENSION(3) :: vel
      LOGICAL :: Stat
      INTEGER :: i,p,q,t,DIM
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

!------------------------------------------------------------------------------
       DIM = CoordinateSystemDimension()
       Force = 0.0d0
       StiffMatrix = 0.0d0
       DampMatrix = 0.0d0
       MassMatrix = 0.0d0
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
       DO t = 1, IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes, U, V, W, SqrtElementMetric, &
                    Basis, dBasisdx )
 
         S = S * SqrtElementMetric
         vel(1) = SUM( V1(1:n) * Basis(1:n) )
         vel(2) = SUM( V2(1:n) * Basis(1:n) )
         vel(3) = SUM( V3(1:n) * Basis(1:n) )

         DO i = 1, VariableComponents
           DO p = 1,N
             DO q = 1,N
             MassMatrix( (p-1)*VariableComponents+i,      &
                 (q-1)*VariableComponents+i) =            &
                 MassMatrix( (p-1)*VariableComponents+i,  &
                 (q-1)*VariableComponents+i)              &
                 + Basis(p) * Basis(q) * s

             temp1 = SUM( vel(1:DIM) * dBasisdx(q,1:DIM) )
             temp2 = SUM( vel(1:DIM) * dBasisdx(p,1:DIM) )

             StiffMatrix( (p-1)*VariableComponents+i,      &
                 (q-1)*VariableComponents+i) =             &
                 StiffMatrix( (p-1)*VariableComponents+i,  &
                 (q-1)*VariableComponents+i)               &
                 + temp1 * temp2 * s
             END DO
           END DO
         END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE LocalMatrixBoundary( StiffMatrix, DampMatrix, MassMatrix, &
         Force, Element, n, VariableComponents, Nodes, V1, V2, V3 )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), &
           MassMatrix(:,:), DampMatrix(:,:),           &
           V1(:), V2(:), V3(:)
       INTEGER :: n, VariableComponents
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,NormVel
       REAL(KIND=dp), DIMENSION(3) :: Vel
       REAL(KIND=dp) :: Normal(3)
       LOGICAL :: Stat
       INTEGER :: i,p,q,t,DIM,k
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
!------------------------------------------------------------------------------
       DIM = CoordinateSystemDimension()
       Force = 0.0d0
       StiffMatrix = 0.0d0
       DampMatrix = 0.0d0
       MassMatrix = 0.0d0
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
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
         vel(1) = SUM( V1(1:n) * Basis(1:n) )
         vel(2) = SUM( V2(1:n) * Basis(1:n) )
         vel(3) = SUM( V3(1:n) * Basis(1:n) )
         Normal = Normalvector(Element, Nodes, U, V, .TRUE.)
         NormVel = SUM( Normal(1:DIM) * Vel(1:DIM) )

         DO i = 1,VariableComponents
           DO p = 1,n
             DO q = 1,N
               DampMatrix( (p-1)*VariableComponents+i,      &
                   (q-1)*VariableComponents+i) =            &
                   DampMatrix( (p-1)*VariableComponents+i,  &
                   (q-1)*VariableComponents+i)              &
                   + Basis(p) * Basis(q) * NormVel* s
             END DO
           END DO
         END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE ModifyBoundaryValues( Model, Name, Component, &
         VariableComponents, U, Perm)
!------------------------------------------------------------------------------
       TYPE(Model_t) :: Model
       CHARACTER(LEN=*) :: Name 
       INTEGER :: Component, VariableComponents, Perm(:)
       REAL(KIND=dp), POINTER :: U(:) 
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: CurrentElement
       INTEGER, POINTER :: NodeIndexes(:)
       INTEGER :: i,j,k,n,t,k1,k2
       LOGICAL :: GotIt
       REAL(KIND=dp) :: Values(Model % MaxElementNodes),s
!------------------------------------------------------------------------------

       DO t = Model % NumberOfBulkElements + 1, &
           Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

         CurrentElement => Model % Elements(t)
         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes

         DO i=1,Model % NumberOfBCs
           IF ( CurrentElement % BoundaryInfo % Constraint == &
               Model % BCs(i) % Tag ) THEN
             Values(1:n) = ListGetReal( Model % BCs(i) % Values, &
                 Name, n, NodeIndexes, gotIt )
             IF ( gotIt ) THEN
               DO j=1,n
                 k = Perm(NodeIndexes(j))
                 IF ( k > 0 ) THEN
                   k = VariableComponents * (k-1) + Component
                   U(k) = Values(j)
                 END IF
               END DO
             END IF
           END IF
         END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE ModifyBoundaryValues
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the matrix containing other coefficients.
!------------------------------------------------------------------------------
     SUBROUTINE Add2ndOrderTime2( MassMatrix, DampMatrix, StiffMatrix,  &
         Force, dt, n, DOFs, NodeIndexes, Solver )
!------------------------------------------------------------------------------
!
! REAL(KIND=dp) :: MassMatrix(:,:)
!   INPUT:
!
! REAL(KIND=dp) :: DampMatrix(:,:)
!   INPUT:
!
! REAL(KIND=dp) :: StiffMatrix(:,:)
!   INOUT:
!   
! REAL(KIND=dp) :: Force(:)
!   INOUT:
!   
! REAL(KIND=dp) :: dt
!   INPUT: Simulation timestep size
!
! INTEGER :: n
!   INPUT: number of element nodes
!
! INTEGER :: DOFs
!   INPUT: variable degrees of freedom
!
! TYPE(Solver_t) :: Solver
!   INPUT: solver parameter list (used to get some options for time integration)
! 
!------------------------------------------------------------------------------
       TYPE(Solver_t) :: Solver

       REAL(KIND=dp) :: MassMatrix(:,:),DampMatrix(:,:), &
           StiffMatrix(:,:),Force(:),dt
       INTEGER :: n,DOFs
       INTEGER :: NodeIndexes(:)
!------------------------------------------------------------------------------
       LOGICAL :: GotIt
       INTEGER :: i,j,k,l
       CHARACTER(LEN=MAX_NAME_LEN) :: Method
       REAL(KIND=dp) :: s,t
       REAL(KIND=dp) :: X(DOFs*n),V(DOFs*N),A(DOFs*N)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
       IF ( Solver % Matrix % Lumped ) THEN
!------------------------------------------------------------------------------
         CALL Fatal('TransportEquationSolver',&
             '"Lumped" option is not available') 
!------------------------------------------------------------------------------
       END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get previous solution vectors and update current force
!-----------------------------------------------------------------------------
       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           L = DOFs * (NodeIndexes(i)-1) + j
           SELECT CASE(Method)
           CASE DEFAULT
             X(K) = Solver % Variable % PrevValues(L,3)
             V(K) = Solver % Variable % PrevValues(L,4)
             A(K) = Solver % Variable % PrevValues(L,5)
           END SELECT
           Solver % Matrix % Force(L,1) = Solver % Matrix % Force(L,1) + &
               Force(K)
         END DO
       END DO
!------------------------------------------------------------------------------
       CALL AverageAccelerationMethod( n*DOFs, dt, MassMatrix, DampMatrix, &
           StiffMatrix, Force, X, V, A )
!------------------------------------------------------------------------------
     END SUBROUTINE Add2ndOrderTime2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
     SUBROUTINE AverageAccelerationMethod( N, dt, MassMatrix, DampMatrix, &
         StiffMatrix, Force, X, V, A )
!------------------------------------------------------------------------------
       INTEGER :: N
       REAL(KIND=dp) :: Force(:),X(:),V(:),A(:),dt
       REAL(KIND=dp) :: Alpha,Beta, Gamma
       REAL(KIND=dp) :: MassMatrix(:,:),DampMatrix(:,:),StiffMatrix(:,:)
!------------------------------------------------------------------------------
       INTEGER :: i,j
       REAL(KIND=dp) :: s
!------------------------------------------------------------------------------

       Alpha = 0.0d0
       Gamma = 0.5d0 - Alpha
       Beta = (1.0d0 - Alpha)**2 / 4.0d0
       DO i=1,N
         s = 0.0d0
         DO j=1,N
           s = s + ( (1.0d0 - Alpha) / (Beta*dt**2) ) * MassMatrix(i,j) * X(j)
           s = s + ( (1.0d0 - Alpha) / (Beta*dt)) * MassMatrix(i,j) * V(j)
           
           s = s + ( Gamma / (Beta*dt) ) * DampMatrix(i,j) * X(j)
           
           s = s - StiffMatrix(i,j) * X(j)

           StiffMatrix(i,j) = StiffMatrix(i,j) +  &
               ( (1.0d0 - Alpha) / (Beta*dt**2) ) * MassMatrix(i,j) + &
               (Gamma / (Beta*dt)) * DampMatrix(i,j)
         END DO
         Force(i) = s
       END DO
!-----------------------------------------------------------------------------
     END SUBROUTINE AverageAccelerationMethod
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE TransportEquationSolver
!------------------------------------------------------------------------------





