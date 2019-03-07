SUBROUTINE EdgeElementSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the best approximation of the vector field U = (-y^3/3,x^3/3)
!  with respect to the energy norm using H(curl)-conforming basis functions. 
!  Additionally, compute the relative error of the solution using
!  the energy norm. Here the energy norm corresponds to the operator 
!  I + curl curl. This solver can thus be used for checking that
!  the convergence rate is quadratic as expected:
!
!             Error 
!  h=1/3   0.29277174E-02
!  h=1/6   0.68241500E-03
!  h=1/12  0.16795856E-03
!  h=1/24  0.41886157E-04

!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt
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
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: Var

  REAL(KIND=dp) :: Norm, u, v, w, Err, EK, SolNorm

  INTEGER :: n, ne, nf, nb, np, nd, t, istat, i, j, k, l, active, dim

  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A

  LOGICAL :: stat, PiolaVersion, ErrorEstimation, UseTabulatedBasis
  LOGICAL :: UseEnergyNorm

  INTEGER, ALLOCATABLE :: Indices(:)

  SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices
!------------------------------------------------------------------------------
  PiolaVersion = .TRUE.
  dim = CoordinateSystemDimension()

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(N), LOAD(6,N), STIFF(N,N), &
          Acoef(N), Indices(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'EdgeElementSolver', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  Solver % Matrix % COMPLEX = .FALSE.
  A => GetMatrix()

  ErrorEstimation = GetLogical( GetSolverParams(), 'Error Computation', Found)
  UseEnergyNorm = GetLogical( GetSolverParams(), 'Use Energy Norm', Found)
  UseTabulatedBasis = GetLogical( GetSolverParams(), 'Tabulate Basis', Found)

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize()

  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes() ! The number of nodes corresponding to the background mesh
     nd = GetElementNOFDOFs()  
     nb = GetElementNOFBDOFs()

     LOAD = 0.0d0
     IF (.FALSE.) THEN     
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) THEN
           Load(1,1:n) = GetReal( BodyForce, 'Body Source 1', Found )
           Load(2,1:n) = GetReal( BodyForce, 'Body Source 2', Found )
           Load(3,1:n) = GetReal( BodyForce, 'Body Source 3', Found )
        END IF
     END IF

     Acoef(1:n) = 1.0d0
     IF (.FALSE.) THEN
        Material => GetMaterial( Element )
        IF ( ASSOCIATED(Material) ) THEN
           Acoef(1:n) = GetReal( Material, 'Conductivity', Found )
           IF (.NOT. Found) CALL Fatal( 'NedelecSolve', 'Conductivity must be specified' )
        END IF
     END IF

     ! Perform an additional check that DOF counts are right:
     !-------------------------------------------------------------------
     IF (GetElementFamily() == 5 .AND. nd /= 6) THEN
        WRITE(Message,'(I2,A)') nd, 'DOFs Found'
        CALL Fatal('EdgeElementSolver','Indices for a tetrahedron erratic') 
     END IF
     IF (GetElementFamily() == 6 .AND. nd /= 10) THEN
        WRITE(Message,'(I2,A)') nd, 'DOFs Found'
        CALL Fatal('EdgeElementSolver','Indices for a pyramid erratic')
     END IF
     IF (GetElementFamily() == 7 .AND. nd /= 15) THEN
        WRITE(Message,'(I2,A)') nd, 'DOFs Found'
        CALL Fatal('EdgeElementSolver','Indices for a prism erratic')
     END IF
     IF (GetElementFamily() == 8 .AND. nd /= 27) THEN
       WRITE(Message,'(I2,A)') nd, 'DOFs Found'
       CALL Fatal('EdgeElementSolver','Indices for a brick erratic')
     END IF


     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix( STIFF, FORCE, LOAD, Acoef, Element, n, nd+nb, dim, UseTabulatedBasis)

     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()  

  !-------------------------------------------------------------------
  ! Compute the error norm for the model problem considered. Note that
  ! this part has not been modified to support the option
  ! Tabulate Basis = Logical True.
  !--------------------------------------------------------------------
  IF (ErrorEstimation) THEN
     Err = 0.0d0
     SolNorm = 0.0d0
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nd = GetElementDOFs( Indices )

        Load(1,1:nd) = Solver % Variable % Values( Solver % Variable % &
             Perm(Indices(1:nd)) )

        CALL MyComputeError(Load, Element, n, nd, dim, Err, SolNorm, UseEnergyNorm)
     END DO

     WRITE (*, '(A,E16.8)') 'Error Norm = ', SQRT(ParallelReduction(Err))/SQRT(ParallelReduction(SolNorm))
  END IF

CONTAINS




!---------------------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Acoef, Element, n, nd, dim, UseTabulatedBasis)
!---------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:,:), Acoef(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    LOGICAL :: UseTabulatedBasis 
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: nu
    REAL(KIND=dp) :: EBasis(nd,3), CurlEBasis(nd,3), F(3,3), G(3,3)
    REAL(KIND=dp) :: dBasisdx(n,3)
    REAL(KIND=dp) :: Basis(n), DetJ, xq, yq, zq, uq, vq, wq, sq
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np

    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    ! ----------------------------------------------------------------------------
    INTEGER :: PermVec(nd)
    REAL(KIND=dp) :: SignVec(nd)
    REAL(KIND=dp) :: ReadyBasis(nd,3), ReadyCurlBasis(nd,3)

    LOGICAL, SAVE :: BricksVisited = .FALSE., PrismsVisited = .FALSE.
    LOGICAL, SAVE :: PyramidsVisited = .FALSE., TetraVisited = .FALSE.

    REAL(KIND=dp), TARGET :: TetraBasis(6,12), TetraCurlBasis(6,12)
    REAL(KIND=dp), TARGET :: PyramidBasis(10,36), PyramidCurlBasis(10,36)
    REAL(KIND=dp), TARGET :: PrismBasis(15,54), PrismCurlBasis(15,54)
    REAL(KIND=dp), TARGET :: BrickBasis(27,81), BrickCurlBasis(27,81)
    REAL(KIND=dp), POINTER :: BasisTable(:,:), CurlBasisTable(:,:)

    SAVE PrismBasis, PrismCurlBasis, BrickBasis, BrickCurlBasis
    SAVE TetraBasis, TetraCurlBasis, PyramidBasis, PyramidCurlBasis
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0

    !-------------------------------------
    ! Numerical integration over element:
    !-------------------------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
         EdgeBasisDegree=2)    

    IF (UseTabulatedBasis .AND. PiolaVersion) THEN
       !----------------------------------------------------------------------------
       ! Obtain the data for permuting basis function positions and applying 
       ! sign reversions. This data is the same for all integration points.
       ! If elements of this type has not yet been visited, tabulate basis
       ! function values to avoid the recomputation.
       !----------------------------------------------------------------------------
       IF ( (GetElementFamily() == 5 .AND. TetraVisited) .OR. &
            (GetElementFamily() == 6 .AND. PyramidsVisited) .OR. &
            (GetElementFamily() == 7 .AND. PrismsVisited) .OR. &            
            (GetElementFamily() == 8 .AND. BricksVisited) ) THEN
          CALL ReorderingAndSignReversionsData(Element,Nodes,PermVec,SignVec)

       ELSE
          !---------------------------------------------------------------------------
          ! The first visit for this element type, tabulate the basis function values.
          !---------------------------------------------------------------------------
          CALL ReorderingAndSignReversionsData(Element,Nodes,PermVec,SignVec)
          DO t=1,IP % n
             stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), DetF=detJ, Basis=Basis, EdgeBasis=EBasis, RotBasis=CurlEBasis, &
                  ApplyPiolaTransform = .FALSE.)
             !------------------------------------------------------------------------
             ! Revert order and sign changes to the reference element default
             ! and tabulate the values for later usage
             !------------------------------------------------------------------------
             DO i=1,nd
                EBasis(i,1:3) = SignVec(i) * EBasis(i,1:3)
                CurlEBasis(i,1:3) = SignVec(i) * CurlEBasis(i,1:3)
             END DO

             SELECT CASE( GetElementFamily() )
             CASE(5)
                DO i=1,nd
                   j = PermVec(i)
                   TetraBasis(i,(t-1)*3+1:(t-1)*3+3) = EBasis(j,1:3)
                   TetraCurlBasis(i,(t-1)*3+1:(t-1)*3+3) = CurlEBasis(j,1:3)               
                END DO
             CASE(6)
                DO i=1,nd
                   j = PermVec(i)
                   PyramidBasis(i,(t-1)*3+1:(t-1)*3+3) = EBasis(j,1:3)
                   PyramidCurlBasis(i,(t-1)*3+1:(t-1)*3+3) = CurlEBasis(j,1:3)               
                END DO
             CASE(7)
                DO i=1,nd
                   j = PermVec(i)
                   PrismBasis(i,(t-1)*3+1:(t-1)*3+3) = EBasis(j,1:3)
                   PrismCurlBasis(i,(t-1)*3+1:(t-1)*3+3) = CurlEBasis(j,1:3)               
                END DO
             CASE(8)
                DO i=1,nd
                   j = PermVec(i)
                   BrickBasis(i,(t-1)*3+1:(t-1)*3+3) = EBasis(j,1:3)
                   BrickCurlBasis(i,(t-1)*3+1:(t-1)*3+3) = CurlEBasis(j,1:3)               
                END DO
             END SELECT
                
          END DO

          SELECT CASE( GetElementFamily() )
          CASE(5)
             TetraVisited = .true.
             CALL Info( 'EdgeElementSolver', 'TETRAHEDRAL BASIS FUNCTIONS WERE TABULATED', Level=10)
             WRITE(Message,'(A,I2)') 'Integration points ', IP % n
             CALL Info('EdgeElementSolver', Message, Level=10) 
          CASE(6)
             PyramidsVisited = .TRUE.
             CALL Info( 'EdgeElementSolve', 'PYRAMIDICAL BASIS FUNCTIONS WERE TABULATED', Level=10)
             WRITE(Message,'(A,I2)') 'Integration points ', IP % n
             CALL Info('EdgeElementSolver', Message, Level=10) 
          CASE(7)
             PrismsVisited = .TRUE.
             CALL Info( 'EdgeElementSolve', 'PRISM BASIS FUNCTIONS WERE TABULATED', Level=10)
             WRITE(Message,'(A,I2)') 'Integration points ', IP % n
             CALL Info('EdgeElementSolver', Message, Level=10) 
          CASE(8)
             BricksVisited = .TRUE.
             CALL Info( 'EdgeElementSolve', 'BRICK BASIS FUNCTIONS WERE TABULATED', Level=10)
             WRITE(Message,'(A,I2)') 'Integration points ', IP % n
             CALL Info('EdgeElementSolver', Message, Level=10) 
          END SELECT
       END IF
    END IF
    
    np = 0  ! Set np = n, if nodal dofs are employed; otherwise set np = 0

    DO t=1,IP % n

       IF (PiolaVersion) THEN
          IF (UseTabulatedBasis) THEN
             SELECT CASE(Element % TYPE % ElementCode / 100)
             CASE(5)
                BasisTable => TetraBasis(1:6,(t-1)*3+1:(t-1)*3+3)
                CurlBasisTable => TetraCurlBasis(1:6,(t-1)*3+1:(t-1)*3+3)
             CASE(6)
                BasisTable => PyramidBasis(1:10,(t-1)*3+1:(t-1)*3+3)
                CurlBasisTable => PyramidCurlBasis(1:10,(t-1)*3+1:(t-1)*3+3)                
             CASE(7)
                BasisTable => PrismBasis(1:15,(t-1)*3+1:(t-1)*3+3)
                CurlBasisTable => PrismCurlBasis(1:15,(t-1)*3+1:(t-1)*3+3) 
             CASE(8)
                BasisTable => BrickBasis(1:27,(t-1)*3+1:(t-1)*3+3)
                CurlBasisTable => BrickCurlBasis(1:27,(t-1)*3+1:(t-1)*3+3)                   
             CASE DEFAULT
                CALL Fatal( 'EdgeElementSolver', 'THE BASIS FUNCTIONS FOR THIS ELEMENT TYPE NONTABULATED' )
             END SELECT
             !----------------------------------------------------------------
             ! Permute, apply sign reversions and apply the Piola transform via
             ! calling EdgeElementInfo:
             !----------------------------------------------------------------
             DO i=1,nd
                ReadyBasis(PermVec(i),1:3) = BasisTable(i,1:3)
                ReadyCurlBasis(PermVec(i),1:3) = CurlBasisTable(i,1:3)
             END DO
             DO i=1,nd 
                ReadyBasis(i,1:3) = SignVec(i) * ReadyBasis(i,1:3)
                ReadyCurlBasis(i,1:3) = SignVec(i) * ReadyCurlBasis(i,1:3)
             END DO
             stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), DetF=detJ, Basis=Basis, EdgeBasis=EBasis, &
                  RotBasis=CurlEBasis, ApplyPiolaTransform = .TRUE., &
                  ReadyEdgeBasis=ReadyBasis, ReadyRotBasis = ReadyCurlBasis) 

          ELSE
             stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), DetF=detJ, Basis=Basis, EdgeBasis=EBasis, &
                  RotBasis=CurlEBasis, ApplyPiolaTransform = .TRUE., &
                  BasisDegree = 2)
          END IF

       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element, EBasis, CurlEBasis, Basis, dBasisdx)
       END IF

       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
       zq = SUM( Nodes % z(1:n) * Basis(1:n) )

       !nu = SUM( Basis(1:n) * Acoef(1:n) )
 
       !----------------------------------------------------------------
       ! The following branch could be used to produce the 
       ! Galerkin projection of a solution component for visualization.
       !------------------------------------------------------------------
       IF (np > 0) THEN
          DO p = 1,n
             DO q = 1,n       
                STIFF(p,q) = STIFF(p,q) + Basis(p) * Basis(q) * detJ * IP % s(t)    
             END DO

             DO q = 1,nd-np
                j = np + q
                ! The following is for plotting the x-component of the solution
                STIFF(p,j) = STIFF(p,j) - EBasis(q,1) * Basis(p) * detJ * IP % s(t)
             END DO
          END DO
       END IF

       !--------------------------------------------------------------
       ! The equations for H(curl)-conforming part...
       !---------------------------------------------------------------
       DO p = 1,nd-np
          !----------------------------
          ! The inner product (u,v)_E
          !----------------------------
          i = np + p
          DO q = 1,nd-np
             j = np + q
             STIFF(i,j) = STIFF(i,j) + 1.0d0 * &
                  SUM( EBasis(q,1:dim) * EBasis(p,1:dim) ) * detJ * IP % s(t) + &
                  CurlEBasis(q,3) * CurlEBasis(p,3) * detJ * IP % s(t)
          END DO
          
          !----------------------------------------
          ! RHS corresponding to the exact solution 
          !----------------------------------------
          FORCE(i) = FORCE(i) - yq**3/3.0d0 * EBasis(p,1) * detJ * IP % s(t) + &
               xq**3/3.0d0  * EBasis(p,2) * detJ * IP % s(t) + &
               (xq**2 + yq**2) * CurlEBasis(p,3) * detJ * IP % s(t)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!----------------------------------------------------------------------------------
  SUBROUTINE MyComputeError(LOAD, Element, n, nd, dim, EK, SolNorm, UseEnergyNorm)
!----------------------------------------------------------------------------------
    REAL(KIND=dp) :: Load(:,:), EK, SolNorm
    TYPE(Element_t), POINTER :: Element    
    INTEGER :: n, nd, dim
    LOGICAL :: UseEnergyNorm
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: EBasis(nd,3), CurlEBasis(nd,3)
    REAL(KIND=dp) :: Basis(n), DetJ, xq, yq, zq, uq, vq, wq, sq, &
         u(3), rotu(3), sol(3), rotsol(3), e(3), rote(3), F(3,3), G(3,3)
    REAL(KIND=dp) :: dBasisdx(n,3)
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    !-------------------------------------
    ! Numerical integration over element:
    !-------------------------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, EdgeBasisDegree=2) 

    np = 0  ! Set np = n, if nodal dofs are employed; otherwise set np = 0

    DO t=1,IP % n

       IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), F, G, detJ, Basis, EBasis, CurlEBasis, ApplyPiolaTransform = .TRUE., &
               BasisDegree = 2)
       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element, EBasis, CurlEBasis, Basis, dBasisdx)
       END IF

       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
!       zq = SUM( Nodes % z(1:n) * Basis(1:n) )

       u = 0.0d0
       u(1) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,1) )
       u(2) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,2) )
!       u(3) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,3) )

       rotu = 0.0d0
!       rotu(1) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,1) )
!       rotu(2) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,2) )
       rotu(3) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,3) )       

       ! Compute the square of the energy norm of the solution and error:
       sol = 0.0d0
       sol(1) = -yq**3/3.0d0
       sol(2) = xq**3/3.0d0

       rotsol = 0.0d0
       rotsol(3) = xq**2 + yq**2
       e(1:2) = sol(1:2) - u(1:2)
       e(3) = 0.0d0
       rote = 0.0d0
       rote(3) = rotsol(3) - rotu(3)

       IF (UseEnergyNorm) THEN
          ! Energy norm 
          SolNorm = SolNorm + (SUM( Sol(1:3) * Sol(1:3) ) + SUM( rotsol(1:3) * rotsol(1:3) )) * detJ
          EK = EK + (SUM( e(1:3) * e(1:3) ) + SUM( rote(1:3) * rote(1:3) )) * detJ
       ELSE
          ! L2 norm
          SolNorm = SolNorm + SUM( Sol(1:3) * Sol(1:3) )* detJ
          EK = EK + SUM( e(1:3) * e(1:3) )* detJ
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE MyComputeError
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE EdgeElementSolver
!------------------------------------------------------------------------------
