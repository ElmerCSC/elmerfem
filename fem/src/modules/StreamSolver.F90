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
!-------------------------------------------------------------------------------
!> Solve the StreamFunction of a two-dimensional steady state flow field.
!> The intended usage of the solver is mainly for postprocessing purposes. 
!> \ingroup Solvers
!-------------------------------------------------------------------------------
SUBROUTINE StreamSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver  !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  TYPE(Variable_t), POINTER :: FlowSol
  TYPE(Element_t), POINTER :: CurrentElement

  CHARACTER(LEN=MAX_NAME_LEN) :: FlowVariableName

  REAL(KIND=dp), POINTER :: StreamFunction(:), FlowSolValues(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)
  REAL(KIND=dp) :: Norm, coeff, val, Anchor

  INTEGER, POINTER :: NodeIndexes(:), Perm(:), FlowSolPerm(:)
  INTEGER :: k, t, i, n, istat, NSDOFs, FirstNode

  TYPE(ValueList_t), POINTER :: SolverParams

  LOGICAL :: AllocationsDone = .FALSE., Found, Shifting, Scaling, StokesStream, &
      DirichletPoint, IsAxis

  SAVE STIFF, FORCE, LOAD, AllocationsDone

  CALL Info('StreamSolver',' ')
  CALL Info('StreamSolver','----------------------------')
  CALL Info('StreamSolver','STREAMSOLVER')
  CALL Info('StreamSolver','----------------------------')
  CALL Info('StreamSolver',' ')

!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT.ASSOCIATED( Solver % Matrix ) ) RETURN

  StreamFunction => Solver % Variable % Values
  Perm => Solver % Variable % Perm

  StiffMatrix => Solver % Matrix
!------------------------------------------------------------------------------
! Get initial values ( Default is FlowSolvers 'Flow Solution' )
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()
  FlowVariableName = GetString( SolverParams, &
     'Stream Function Velocity Variable', Found )
  IF ( .NOT. Found ) THEN
     CALL Info( 'StreamSolver', 'Stream Function Velocity Variable set to Flow Solution' )
     FlowVariableName = 'Flow Solution'
  END IF

  FlowSol => VariableGet( Solver % Mesh % Variables, FlowVariableName )
  IF ( ASSOCIATED( FlowSol ) ) THEN
     FlowSolPerm => FlowSol % Perm
     FlowSolValues => FlowSol % Values
     NSDOFs = FlowSol % DOFs
  ELSE
     CALL Warn( 'StreamSolver', 'No variable for velocity associated.' )
     CALL Warn( 'StreamSolver', 'Quitting execution of StreamSolver.' ) 
     RETURN
  END IF

!------------------------------------------------------------------------------
! Get keyword values
!------------------------------------------------------------------------------
  n = Solver % Mesh % NumberOfNodes
  FirstNode = GetInteger( SolverParams, 'Stream Function First Node', Found )
  IF ( Found ) THEN
    IF ( FirstNode > n ) THEN
      CALL Warn( 'StreamSolver', 'Given > Stream Function First Node < is too big.' )
      WRITE( Message, *) '> Stream Function First Node < set to ', n
      CALL Info( 'StreamSolver', Message )
      FirstNode = n
    ELSE IF ( FirstNode < 1 ) THEN
      CALL Warn( 'StreamSolver', 'Given > Stream Function First Node < is non-positive.' )
      CALL Info( 'StreamSolver', '> Stream Function First Node < set to 1.' )
      FirstNode = 1
    END IF
  ELSE
    CALL Info( 'StreamSolver', 'Keyword > Stream Function First Node < not given, use Dirichlet conditions!' )
  END IF

  Shifting = GetLogical( SolverParams, 'Stream Function Shifting', Found )
  IF ( .NOT. Found ) THEN
     CALL Info( 'StreamSolver', '> Stream Function Shifting < set to .TRUE.' )
     Shifting = .TRUE.
  END IF

  Scaling = GetLogical( SolverParams, 'Stream Function Scaling', Found )
  IF ( .NOT. Found ) THEN
     CALL Info( 'StreamSolver', '> Stream Function Scaling < set to .FALSE.' )
     Scaling = .FALSE.
  END IF

  DirichletPoint = GetLogical( SolverParams, 'Dirichlet Point', Found )
  Coeff = GetConstReal( SolverParams, 'Dirichlet Weight', Found )
  IF(.NOT. Found) Coeff = 1.0_dp

  Anchor = GetConstReal( SolverParams, 'Stream Function Penalty', Found )

  StokesStream = GetLogical( SolverParams, 'Stokes Stream Function', Found )
  IsAxis = CurrentCoordinateSystem() == AxisSymmetric .OR. &
         CurrentCoordinateSystem() == CylindricSymmetric
  IF ( Found ) THEN
    IF( IsAxis .AND. .NOT. StokesStream ) THEN
      CALL Warn( 'StreamSolver', 'Using normal stream function in axis symmetric case.' )
    ELSE IF ( .NOT. IsAxis .AND. StokesStream ) THEN
      CALL Warn( 'StreamSolver', 'Using Stokes stream function in cartesian case.' )
    END IF
  ELSE
    StokesStream = IsAxis
  END IF
  
!------------------------------------------------------------------------------
! Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays

     ALLOCATE( FORCE( N ),   LOAD( 2*N ), STIFF(N,N), STAT=istat ) 

     IF ( istat /= 0 ) CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     AllocationsDone = .TRUE.
  END IF
!------------------------------------------------------------------------------
! Initialize the system and do the assembly
!------------------------------------------------------------------------------
  CALL DefaultInitialize()
  
  DO t=1,Solver % NumberOfActiveElements
     CurrentElement => GetActiveElement(t)
     n = GetElementNOFNodes()
     NodeIndexes => CurrentElement % NodeIndexes
!------------------------------------------------------------------------------
     LOAD = 0.0d0
     DO i = 1,n
        k = FlowSolPerm( NodeIndexes(i) )
        k = (k-1) * NSDOFs
        LOAD( (i-1)*2+1 ) =  FlowSolValues( k+2 )
        LOAD( (i-1)*2+2 ) = -FlowSolValues( k+1 )
     END DO
!------------------------------------------------------------------------------
!     Get element local matrix and rhs vector
!------------------------------------------------------------------------------
     CALL LocalMatrix(  STIFF, FORCE, LOAD, CurrentElement, n, StokesStream )
!------------------------------------------------------------------------------
!     Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO! <- elements
   CALL DefaultFinishBulkAssembly()

   DO t=1,Solver % Mesh % NumberOfBoundaryElements
     CurrentElement => GetBoundaryElement(t)
     IF ( .NOT. ActiveBoundaryElement() ) CYCLE
     n = GetElementNOFNodes()
     NodeIndexes => CurrentElement % NodeIndexes
!------------------------------------------------------------------------------
     LOAD = 0.0d0
     DO i = 1,n
        k = FlowSolPerm( NodeIndexes(i) )
        k = (k-1) * NSDOFs
        LOAD( (i-1)*2 +1 ) = -FlowSolValues( k+2 )
        LOAD( (i-1)*2 +2 ) =  FlowSolValues( k+1 )
     END DO
!------------------------------------------------------------------------------
!    Get element local matrix and rhs vector
!------------------------------------------------------------------------------
     CALL LocalMatrixBC(  STIFF, FORCE, LOAD, CurrentElement, n, StokesStream )
!------------------------------------------------------------------------------
!     Update global matrix and rhs vector from local matrix & vector
!------------------------------------------------------------------------------
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO! <- elements
!------------------------------------------------------------------------------
   CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
!  Zero the row corresponding to the 'FirstNode':
!------------------------------------------------------------------------------
   IF( FirstNode > 0 ) THEN
     k = Solver % Variable % Perm( FirstNode )
     CALL UpdateDirichletDof( StiffMatrix, k, 0.0_dp )

     !val =  MAXVAL( ABS( StiffMatrix % Values) ) 
     !IF( DirichletPoint ) THEN
       

     !CALL ZeroRow( StiffMatrix, k )
     !StiffMatrix % RHS(k) = 0.0_dp
     !END IF
     !CALL AddToMatrixElement( StiffMatrix,k,k,val*coeff )
   END IF
   CALL DefaultDirichletBCs()
   !END IF
!------------------------------------------------------------------------------
!  Solve the system:
!------------------------------------------------------------------------------
   Norm = DefaultSolve()
!------------------------------------------------------------------------------
! Do Shifting and Scaling if needed
!------------------------------------------------------------------------------

  IF ( Shifting ) THEN
     StreamFunction = StreamFunction - MINVAL( StreamFunction )
  END IF

  IF ( Scaling ) THEN
    coeff = MAXVAL( StreamFunction ) - MINVAL( StreamFunction ) 
    IF ( coeff < AEPS ) THEN
      CALL Warn( 'StreamSolver', &
          'Maximum absolute value smaller than machine epsilon; cannot scale.' )
    ELSE
      StreamFunction = StreamFunction / coeff
    END IF
  END IF
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n,StokesStream )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: StokesStream
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), LoadAtIp(2)
    REAL(KIND=dp) :: DetJ, U, V, W, S, Radius
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim, i,j,k
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    STIFF  = 0.0d0
    FORCE  = 0.0d0
    CALL GetElementNodes( Nodes )
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )
    
    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
       
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, DetJ, Basis, dBasisdx )
       s = s * DetJ
       
!------------------------------------------------------------------------------
!      Load at the integration point
!------------------------------------------------------------------------------
       IF ( StokesStream ) THEN
         Radius = SUM(Nodes % x(1:n) * Basis(1:n))
         DO i = 1,dim
           LoadAtIP(i) = Radius * SUM( dBasisdx(1:n,i) * LOAD(i:2*n:2) )
         END DO
         LoadAtIP(1) = LoadAtIp(1) + SUM( Basis(1:n) * LOAD(1:2*n:2) )
       ELSE
         DO i = 1,dim
           LoadAtIP(i) = SUM( dBasisdx(1:n,i) * LOAD(i:2*n:2) )
         END DO
       END IF
!------------------------------------------------------------------------------
!      Finally, the elemental matrix & vector
!------------------------------------------------------------------------------       
       STIFF(1:n,1:n) = STIFF(1:n,1:n) &
           + s * MATMUL( dBasisdx, TRANSPOSE(dBasisdx) )

       DO i=1,n
         DO j=1,n
           STIFF(i,j) = STIFF(i,j) + s * Anchor * Basis(i) * Basis(j)
         END DO
       END DO

       FORCE(1:n) = FORCE(1:n) + s * SUM(LoadAtIp(1:dim)) * Basis(1:n)
       
    END DO! <- t eli integraatiopisteet

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( STIFF, FORCE, LOAD, Element, n,StokesStream )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: StokesStream
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), LoadAtIp(2)
    REAL(KIND=dp) :: DetJ, U, V, W, S, L, Normal(3)
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim, k
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    STIFF  = 0.0d0
    FORCE  = 0.0d0
    CALL GetElementNodes( Nodes )
!------------------------------------------------------------------------------
!   Numerical integration
!------------------------------------------------------------------------------
    IntegStuff = GaussPoints( Element )

    DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)
       
!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, DetJ, &
            Basis, dBasisdx )
       s = s * DetJ
       IF ( StokesStream ) s = s * SUM(Nodes % x(1:n)*Basis(1:n))
       
!------------------------------------------------------------------------------
!      Load at the integration point
!------------------------------------------------------------------------------
       LoadAtIP(1) = SUM( Basis(1:n) * LOAD(1:2*n:2) )
       LoadAtIP(2) = SUM( Basis(1:n) * LOAD(2:2*n:2) )
!------------------------------------------------------------------------------
!      Finally, the elemental matrix & vector
!------------------------------------------------------------------------------       
       Normal = NormalVector( Element, Nodes, u, v, .TRUE. )
       L = SUM( LoadAtIp(1:2) * Normal(1:2) )

       FORCE(1:n) = FORCE(1:n) + s * L * Basis(1:n)
    END DO! <- t eli integraatiopisteet

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE StreamSolver
!------------------------------------------------------------------------------
