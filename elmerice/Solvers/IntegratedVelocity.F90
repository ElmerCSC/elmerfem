!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini
! *  Email:   gagliar@lgge.obs.ujf-grenoble.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 12 July 2010
! * 
! *****************************************************************************
!>   Solver to compute the vertically integrated horizontal velocity 
!>   Exported Variables IntVelocity dof = DIM-1
!>   Needs a velocity field from a NS type solver 
SUBROUTINE IntegratedVelocity( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!   Compute the depth integrated horizontal velocities = sum_zb^zs u dz
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
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
  TYPE(Element_t),POINTER :: CurrentElement, Element, ParentElement, &
                             BoundaryElement
  TYPE(Matrix_t),POINTER  :: StiffMatrix
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material
  TYPE(Variable_t), POINTER :: PointerToVariable, IntVeloSol, FlowVariable, HeightSol

  LOGICAL :: AllocationsDone = .FALSE., Found, OnSurface = .TRUE.,UnFoundFatal=.TRUE.

  INTEGER :: i, n, m, t, istat, DIM, COMP, other_body_id   
  INTEGER, POINTER :: Permutation(:), NodeIndexes(:), IntVeloPerm(:),&
                      FlowPerm(:), HeightPerm(:)
       
  REAL(KIND=dp), POINTER :: ForceVector(:)
  REAL(KIND=dp), POINTER :: VariableValues(:), IntVelo(:), FlowValues(:), Height(:)
  REAL(KIND=dp) :: Norm, cn, dd 

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
           NodalVar(:) 

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FlowSolverName

  SAVE STIFF, LOAD, FORCE, AllocationsDone, DIM, SolverName
  SAVE NodalVar 
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  WRITE(SolverName, '(A)') 'IntegratedVelocity'

  IntVeloSol => VariableGet( Solver % Mesh % Variables, 'Integrated Velocity',UnFoundFatal=UnFoundFatal)
  IntVelo => IntVeloSol % Values
  IntVeloPerm => IntVeloSol % Perm
  
  
  OnSurface = GetLogical( Solver % Values, 'On Surface', Found )
  IF (.Not.Found) OnSurface = .True.  !Default Value is True
  
  IF (OnSurface) THEN
     HeightSol => VariableGet( Solver % Mesh % Variables, 'Height',UnFoundFatal=UnFoundFatal)
     Height => HeightSol % Values
     HeightPerm => HeightSol % Perm
  ELSE
     HeightSol => VariableGet( Solver % Mesh % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
     Height => HeightSol % Values
     HeightPerm => HeightSol % Perm
  END IF
  
!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS, AIFlow, Porous, ...)
!------------------------------------------------------------------------------
      
    FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', Found )    
    IF (.NOT.Found) FlowSolverName = 'Flow Solution'
    FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
    IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
    ELSE
       CALL Info(SolverName, &
                      & 'No variable for velocity associated.', Level=4)
    END IF 
  
!--------------------------------------------------------------
!Allocate some permanent storage, this is done first time only:
!--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF, NodalVar) 

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), NodalVar(N), &
                          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( SolverName, 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
  END IF

     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS

DIM = CoordinateSystemDimension()

! Loop for ux and uy (in 3D)
DO COMP=1, DIM-1
! No non-linear iteration, no time dependency  
  VariableValues = 0.0d0
  Norm = Solver % Variable % Norm

  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
     n = GetElementNOFNodes()

     NodeIndexes => Element % NodeIndexes

     NodalVar(1:n) = FlowValues((DIM+1)*(FlowPerm(NodeIndexes(1:n))-1) + COMP)
     IF (.Not.OnSurface) NodalVar = -NodalVar

     CALL LocalMatrix (  STIFF, FORCE, Element, n , NodalVar )
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  

  ! Neumann conditions 
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     BoundaryElement => GetBoundaryElement(t)
     IF ( GetElementFamily() == 1 ) CYCLE
     NodeIndexes => BoundaryElement % NodeIndexes
     IF (ParEnv % myPe .NE. BoundaryElement % partIndex) CYCLE
     n = GetElementNOFNodes()

     
! Read the nodal velocity 
     NodalVar(1:n) = FlowValues((DIM+1)*(FlowPerm(NodeIndexes(1:n))-1) + COMP)
     IF (.Not.OnSurface) NodalVar = -NodalVar
          

     STIFF = 0.0D00
     FORCE = 0.0D00
     CALL LocalMatrixBC(  STIFF, FORCE, BoundaryElement, n, NodalVar)
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishAssembly()
  
! Dirichlet 
  CALL SetDirichletBoundaries( Model, StiffMatrix, ForceVector, &
          ComponentName('Integrated Velocity',COMP), 1,1, Permutation)

  CALL EnforceDirichletConditions( Solver, Solver % Matrix , Solver % Matrix % RHS)
  
  Norm = DefaultSolve()

  ! Save the solution on the right variable
  
  DO i = 1, Model % Mesh % NumberOfNodes
        IF (IntVeloPerm(i)>0) THEN
            IntVelo((DIM-1)*(IntVeloPerm(i)-1) + COMP) = VariableValues(Permutation(i)) 
            IF (Height(HeightPerm(i))>0.0_dp) IntVelo((DIM-1)*(IntVeloPerm(i)-1) + COMP) = &
                    IntVelo((DIM-1)*(IntVeloPerm(i)-1) + COMP) / Height(HeightPerm(i))
        END IF
  END DO
  
  
END DO !COMP


CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n, var)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:) , var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3), detJ,grad
    LOGICAL :: Stat
    INTEGER :: t, p,q ,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       grad  = SUM( var(1:n) * dBasisdx(1:n,dim) )
       FORCE(1:n) = FORCE(1:n) - grad * IP % s(t) * DetJ  * Basis(1:n)

       DO p=1,n
         DO q=1,n
           STIFF(p,q) = STIFF(p,q) + IP % S(t) * detJ * dBasisdx(q,dim)*dBasisdx(p,dim)
         END DO
       END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, Element, n, var ) 
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), var(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3), eta, grad 
    LOGICAL :: Stat
    INTEGER :: t, dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = CoordinateSystemDimension()

    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       grad  = SUM( var(1:n) * Basis(1:n) )

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + grad * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------
END SUBROUTINE IntegratedVelocity
!------------------------------------------------------------------------------


