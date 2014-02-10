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
!>  This subroutine solves the Galerkin approximation to the rate of change 
!>  of the field subject to the convection operator. This approximation 
!>  is used by the wave equation solver for the convective transport 
!>  equation (TransportEquationSolver). The field subject to the convection 
!>  operator is declared using sif-file flag Advection Variable and must be 
!>  available to the solver. 
!> \ingroup Solvers
!------------------------------------------------------------------------------
   SUBROUTINE RateOfChangeSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
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
     INTEGER :: i, j, k, p, n, t, body_id, istat, &
         LocalNodes, AdvectionVariableComponents, VelocityComponents
     INTEGER, SAVE :: ThisSolverCalls=0
     INTEGER, POINTER :: NodeIndexes(:), Udot0Perm(:), VelocityPerm(:)
 
     TYPE(Matrix_t), POINTER :: StiffMatrix
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: CurrentElement
     TYPE(Variable_t), POINTER :: VelocitySol, FlowSol
     TYPE(VariablePtr_t), POINTER :: U0Var(:)
     TYPE(ValueList_t), POINTER :: Material


     LOGICAL :: GotIt, AllocationsDone = .FALSE.
     CHARACTER(LEN=MAX_NAME_LEN) :: AdvectionFlag, AdvectionVariable, &
         VariableName, EquationName

     REAL(KIND=dp), POINTER :: Udot0(:), ForceVector(:), Flow(:), &
         Velocity(:)
     REAL(KIND=dp), ALLOCATABLE :: &
         LocalStiffMatrix(:,:), LocalForce(:), LocalU0(:), &
         V1(:), V2(:), V3(:)
     REAL(KIND=dp) :: at, st, CPUTime, Norm

     SAVE LocalStiffMatrix, LocalForce, ElementNodes, & 
          AllocationsDone, LocalU0, U0Var, V1, V2, V3

!------------------------------------------------------------------------------
     ThisSolverCalls = ThisSolverCalls + 1
     
     IF (MOD(ThisSolverCalls,2)==1) THEN
       CALL Info('RateOfChangeSolver', '') 
       CALL Info('RateOfChangeSolver', 'Starting initialization...') 
       CALL Info('RateOfChangeSolver', '') 

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
       Udot0     => Solver % Variable % Values
       Udot0Perm => Solver % Variable % Perm

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

         ALLOCATE( ElementNodes % x( N ),                            &
             ElementNodes % y( N ),                                  &
             ElementNodes % z( N ),                                  &
             LocalU0(N*AdvectionVariableComponents),                 &
             V1(N),                                                  &
             V2(N),                                                  &
             V3(N),                                                  &
             U0Var(AdvectionVariableComponents),                     &
             LocalForce(N*AdvectionVariableComponents),              &
             LocalStiffMatrix( N*AdvectionVariableComponents,        &
             N*AdvectionVariableComponents ),                        &
             STAT=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal('RateOfChangeSolver', &
               'Memory allocation error, aborting.')
         END IF

         AllocationsDone = .TRUE.
       END IF

!-----------------------------------------------------------------------------
!    Get the field for which the rate of change is computed here
!-----------------------------------------------------------------------------
       AdvectionVariable = ListGetString(Solver % Values, &
           'Advection Variable')

       IF (AdvectionVariableComponents == 1) THEN
         U0Var(1) % Var => VariableGet(Solver % Mesh % Variables, &
             AdvectionVariable)
         IF (ASSOCIATED( U0Var(1) % Var )) THEN
           IF (AdvectionVariableComponents /= U0Var(1) % Var % DOFs) &
               CALL Fatal('RateOfChangeSolver', &
               'Solver Variable and Advection Variable Dofs do not equal')
         ELSE
           CALL Fatal('RateOfChangeSolver', &
               'Specified Advection Variable does not exist')
         END IF
       ELSE
         DO i = 1, AdvectionVariableComponents
           WRITE(VariableName,'(A,A,I1)') TRIM(AdvectionVariable),' ', i
           U0Var(i) % Var => VariableGet(Solver % Mesh % Variables, &
               VariableName)
           IF (.NOT. ASSOCIATED( U0Var(i) % Var )) &
               CALL Fatal('RateOfChangeSolver', &
               'Specified Advection Variable component does not exist')
         END DO
       END IF

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
               CALL Warn('RateOfChangeSolver', &
               'Coordinate system and Advection Velocity dimensions unequal')
         ELSE
           CALL Fatal('RateOfChangeSolver', &
               'Flow Solution does not exist')
         END IF
       ELSE
         IF (AdvectionFlag /= 'constant') CALL Fatal(     &
             'RateOfChangeSolver', &
             'Advection flag should be either "computed" or "constant"')
         VelocityComponents = CoordinateSystemDimension() 
       END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------
       at = CPUTime()
       CALL InitializeToZero( StiffMatrix, ForceVector )
       EquationName = ListGetString( Solver % Values, 'Equation' )

!------------------------------------------------------------------------------
!    Do the assembly
!------------------------------------------------------------------------------
       DO t=1,Model % NumberOfBulkElements

         CurrentElement => Solver % Mesh % Elements(t)
         IF ( .NOT. CheckElementEquation( Model, &
             CurrentElement, EquationName ) ) CYCLE

         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes
         ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
         ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
         ElementNodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)

         DO p=1, AdvectionVariableComponents
           DO i=1,n
             j = Udot0Perm( NodeIndexes(i) ) 
             k = U0Var(p) % Var % Perm(NodeIndexes(i))
             IF(j>0 .AND. k>0) THEN
               LocalU0( (i-1)*AdvectionVariableComponents + p ) = &
                   U0Var(p) % Var % Values(k)
             ELSE
               CALL Fatal('RateOfChangeSolver', &
                   'Nonmatching variable permutations')
             END IF
           END DO
         END DO

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
         CALL LocalMatrix(  LocalStiffMatrix, LocalForce, &
             CurrentElement, n, AdvectionVariableComponents, ElementNodes, &
             LocalU0, V1, V2, V3 )

!------------------------------------------------------------------------------
!      Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
         CALL UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
             ForceVector, LocalForce, n, AdvectionVariableComponents, &
             Udot0Perm(NodeIndexes) ) 
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
!    FinishAssembly must be called after all other assembly steps, but before
!    Dirichlet boundary settings. Actually no need to call it except for
!    transient simulations.
!------------------------------------------------------------------------------
!    
!       CALL FinishAssembly( Solver,ForceVector )
!
!------------------------------------------------------------------------------
!    Dirichlet boundary conditions
!------------------------------------------------------------------------------
       DO i = 1, AdvectionVariableComponents
         VariableName = ComponentName(Solver % Variable,i)
         CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
             VariableName, i, AdvectionVariableComponents, Udot0Perm )
       END DO

!------------------------------------------------------------------------------
       at = CPUTime() - at
!------------------------------------------------------------------------------
!    Solve the system and we are done.
!------------------------------------------------------------------------------
       st = CPUTime()

       CALL SolveSystem( StiffMatrix, ParMatrix, ForceVector, &
           Udot0, Norm, AdvectionVariableComponents, Solver )  

       st = CPUTime() - st

       WRITE(Message,'(a,F8.2)') ' Assembly: (s)', at
       CALL Info( 'RateOfChangeSolver', Message )
       WRITE(Message,'(a,F8.2)') ' Solve:    (s)', st
       CALL Info( 'RateOfChangeSolver', Message )

     END IF

!------------------------------------------------------------------------------
   CONTAINS
!------------------------------------------------------------------------------
     SUBROUTINE LocalMatrix( StiffMatrix, Force, Element, n, &
         VariableComponents, Nodes, LocalU0, V1, V2, V3 )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: StiffMatrix(:,:), Force(:), LocalU0(:), &
            V1(:), V2(:), V3(:)
       REAL(KIND=dp), DIMENSION(3) :: vel
       INTEGER :: n, VariableComponents
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,temp
       LOGICAL :: Stat

       INTEGER :: i,p,q,t,DIM
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------

       DIM = CoordinateSystemDimension()
       Force = 0.0d0
       StiffMatrix = 0.0d0
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

         DO i = 1,VariableComponents
           DO p = 1,N
             DO q = 1,N
               StiffMatrix( (p-1)*VariableComponents+i, &
                   (q-1)*VariableComponents+i ) = &
                   StiffMatrix( (p-1)*VariableComponents+i, &
                   (q-1)*VariableComponents+i) &
                   + Basis(q) * Basis(p) * s

               temp = SUM( vel(1:DIM) * dBasisdx(q,1:DIM) )
               Force( (p-1)*VariableComponents + i) = &
                   Force( (p-1)*VariableComponents + i) - &
                   temp * LocalU0((q-1)*VariableComponents + i) * Basis(p) * s 
             END DO
           END DO
         END DO
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE RateOfChangeSolver
!------------------------------------------------------------------------------





