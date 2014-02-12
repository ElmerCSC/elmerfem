SUBROUTINE FlowdepthSolver( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: FlowdepthSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Flowdepth equation!
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
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable
  TYPE(Solver_t), POINTER :: PointerToSolver

  LOGICAL :: AllocationsDone = .FALSE., Found, CalcFree = .FALSE.

  INTEGER :: i, n, m, t, istat
  INTEGER, POINTER :: Permutation(:), NumberOfVisits(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), Surface(:), GradSurface1(:),GradSurface2(:)
  REAL(KIND=dp) :: Norm, Gradient,GradSurface(3)

  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)
       

  SAVE STIFF, LOAD, FORCE, Surface, GradSurface, AllocationsDone
!------------------------------------------------------------------------------
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  ! get Gradient (change of property
  ! with respect to unit-length of
  ! vertical direction
  !----------------------------------
  SolverParams => GetSolverParams()
  Gradient = GetConstReal( SolverParams, &
                      'Gradient',  Found )
  IF (.NOT. Found) THEN
     CALL WARN('FlowdepthSolve', 'No keyword >Gradient< found in section Solver')
     CALL WARN('FlowdepthSolve', 'Assuming value of -1')
     Gradient = -1.0D00
  ELSE
     WRITE(Message,'(A e12.4,A)') 'Gradient of ',Gradient,' applied'
     CALL INFO('FlowdepthSolve', Message,Level=1)
  END IF

  CalcFree = GetLogical(SolverParams, 'Calc Free Surface', Found)
  IF (.NOT. Found) THEN
     CalcFree = .FALSE.
  ELSE
     CALL INFO('FlowdepthSolve', 'Free surface variable will be calculated', Level=1)
  END IF

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed  ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     M = Model % Mesh % NumberOfNodes
     IF (AllocationsDone) DEALLOCATE(FORCE, LOAD, STIFF,&
          Surface, GradSurface1, GradSurface2, NumberOfVisits)

     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N),&
          Surface(M), GradSurface1(M), GradSurface2(M), NumberOfVisits(M),&
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'FlowdepthSolve', 'Memory allocation error.' )
     END IF
     
     IF (CalcFree) THEN
     ! Assign Variable for Residual (i.e., heat flux at boundaries)
     !-------------------------------------------------------------
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
             'FreeSurf', 1, Surface, Permutation )
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
             'FreeSurfGrad1', 1, GradSurface1, Permutation )
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
             'FreeSurfGrad2', 1, GradSurface2, Permutation )
     END IF

     AllocationsDone = .TRUE.
  END IF

  Surface = 1.0D00
  GradSurface1 = 2.0D00
  GradSurface2 = 3.0D00


  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()

  ! bulk assembly
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL LocalMatrix(  STIFF, FORCE, Element, n)
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  ! vN. conditions
  DO t=1,Solver % Mesh % NUmberOfBoundaryElements
     Element => GetBoundaryElement(t)
     n = GetElementNOFNodes()
     IF ( GetElementFamily() /= 1 ) THEN
        CALL LocalMatrixBC(  STIFF, FORCE, Element, n, Gradient )
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END IF
  END DO
  CALL DefaultFinishAssembly()

  ! Dirichlet 
  CALL DefaultDirichletBCs()
  !Solve the system
  Norm = DefaultSolve()

  IF (Calcfree) THEN   
     Surface = 0.0D00
     GradSurface1 = 0.0D00
     GradSurface2 = 0.0D00
     NumberOfVisits = 0
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        CALL GetSurfaceValue(Model, Surface, GradSurface1, GradSurface2,&
             VariableValues, Permutation, NumberOfVisits, Element, n)
     END DO
     DO i=1,Model % Mesh % NumberOfNodes
        GradSurface1(Permutation(i)) = GradSurface1(Permutation(i))/NumberOfVisits(i)
        GradSurface2(Permutation(i)) = GradSurface2(Permutation(i))/NumberOfVisits(i)
     END DO
  END IF
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE GetSurfaceValue(Model, Surface, GradSurface1, GradSurface2,&
       VariableValues, Permutation, NumberOfVisits, Element, n)
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: LocalSurf, GradSurface(3), Depth
    INTEGER :: n
    INTEGER, POINTER :: Permutation(:), NumberOfVisits(:)
    REAL(KIND=dp), POINTER :: Surface(:), GradSurface1(:), GradSurface2(:), VariableValues(:)
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),DetJ, z, U, V, W, SqrtElementMetric
    LOGICAL :: Stat
    INTEGER :: i,j,k,dim
    LOGICAL :: FirstTime
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    dim = CoordinateSystemDimension()
    ! loop over all nodes in Element
    DO i=1,N

       j = Element % NodeIndexes(i) ! get number of node in element in physical space
       NumberOfVisits(j) = NumberOfVisits(j) + 1

       ! get local coordinates of the point i inside the element
       U = Element % Type % NodeU(i)
       V = Element % Type % NodeV(i)
       W = Element % Type % NodeW(i)

       ! get local information on test-functions and derivatives of the point i
       stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
            Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE. )  
       IF (DIM == 2) THEN 
          z = Model % Nodes % y(j)
       ELSE IF (DIM == 3) THEN
          z = Model % Nodes % z(j)
       ELSE
          CALL FATAL('FlowdepthSolve', 'Flow depth for one-dimensional problem not defined!')
       END IF

       IF (NumberOfVisits(j) == 1) &
            Surface(Permutation(j)) = z + VariableValues(Permutation(j))  
!       DO k=1,n
!          GradSurface1(Permutation(j)) = GradSurface1(Permutation(j)) + dBasisdx(k,1)*VariableValues(Permutation(Element % NodeIndexes(k)))
!       END DO
       GradSurface1(Permutation(j)) = GradSurface1(Permutation(j)) +&
            SUM(dBasisdx(1:N,1)*VariableValues(Permutation(Element % NodeIndexes(1:N))))
       IF (DIM > 2) &
            GradSurface2(Permutation(j)) = GradSurface2(Permutation(j)) +&
            SUM(dBasisdx(1:N,2)*VariableValues(Permutation(Element % NodeIndexes(1:N))))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetSurfaceValue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),DetJ
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
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, Element, n, Gradient)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), Gradient
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), &
                      DetJ,Normal(3)
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

      Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE.)
      FORCE(1:n) = FORCE(1:n) + Gradient * IP % s(t) * DetJ * Normal(dim) * Basis(1:n)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------
END SUBROUTINE FlowdepthSolver
!------------------------------------------------------------------------------

SUBROUTINE getDistance( Model,Solver,dt,TransientSimulation )

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
  TYPE(Element_t),POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE., Found, Converged

  INTEGER :: n, t, istat, other_body_id
  REAL(KIND=dp) :: Norm, NormalGradient

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), FORCE(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType

  SAVE MASS, STIFF, LOAD, FORCE, AllocationsDone
!------------------------------------------------------------------------------


  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays

     ALLOCATE( FORCE(N), LOAD(N), MASS(N,N), STIFF(N,N), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'getDistance',&
             'Memory allocation error for matrix/vectors.' )
     END IF

     AllocationsDone = .TRUE.
  END IF

  !Read in solver parameters
  !-------------------------
  SolverParams => GetSolverParams()
  IF (.NOT. ASSOCIATED(SolverParams))&
       CALL FATAL('getDistance','No Solver section found')
  
  !Initialize the system and do the assembly:
  !------------------------------------------
  CALL DefaultInitialize()
  
  ! Assembly for the domain
  !------------------------
  DO t=1,Solver % NumberOfActiveElements
     
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     ! get material parameters
     Material => GetMaterial()
     IF (.NOT. ASSOCIATED(Material)) THEN
        WRITE(Message,'(A,I5,A)') 'No material for bulk element no. ',t,' found.'
        CALL FATAL('getDistance',Message)
     END IF
     
     ! get load for force vector
     LOAD = 1.0d00
     
     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, Element, n, TransientSimulation)
     
     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )
      
     !------------------------------------------------------------------------------
  END DO
!------------------------------------------------------------------------------
!     assembly of von Neumann boundary conditions
!------------------------------------------------------------------------------
  DO t=1, Solver % Mesh % NumberOfBoundaryElements
     
     Element => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement() ) CYCLE
     n = GetElementNOFNodes()
     IF ( GetElementFamily() == 1 ) CYCLE ! no von Neumann BC's on points
     BC => GetBC()
     
     FORCE = 0.0d00
     MASS = 0.0d00
     STIFF = 0.0d00
     LOAD = 0.0d00
     
     NormalGradient = GetConstReal(BC,'Normal Gradient', Found) 
     ! natural boundary condition
     IF (.NOT. Found) CYCLE
     ! natural boundary condition
     LOAD(1:n) = LOAD(1:n) + NormalGradient            
     ! do the assembly of the force vector
     CALL BoundaryCondition(LOAD, FORCE, Element, n)
     IF ( TransientSimulation ) THEN
        MASS = 0.0d00
        CALL Default1stOrderTime( MASS, STIFF, FORCE )
     END IF
     
     CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
  END DO ! end of assembly of von Neumann boundary conditions
!------------------------------------------------------------------------------
  CALL DefaultFinishAssembly()
!
!    Dirichlet BCs:
!    --------------
  CALL DefaultDirichletBCs()

  ! Solve the system
  ! ----------------
  Norm = DefaultSolve()
     
  WRITE( Message, * ) 'Result Norm   : ',Norm
  CALL Info( 'getDistance', Message, Level=4 )

!------------------------------------------------------------------------------


CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Element, n, TransientSimulation)
    IMPLICIT NONE
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:,:) :: MASS, STIFF
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
!    REAL(KIND=dp) :: HeatCapacity(:), HeatConductivity(:), Density(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: detJ, LoadAtIP
    LOGICAL :: Stat, getSecondDerivatives
    INTEGER :: t,i,j,DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
!-----------------------------------------------------------------
!   Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
    DO t=1,IP % n
       !Basis function values & derivatives at the integration point:
       !-------------------------------------------------------------
       getSecondDerivatives = .FALSE.
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)
       !The source term at the integration point:
       !-----------------------------------------
       LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
!-----------------------------------------------------------------    
!      Loop over j-components of matrices and over force vector
!-----------------------------------------------------------------    
       DO j=1,n
          FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
!-----------------------------------------------------------------    
!         Loop over i-components of matrices
!-----------------------------------------------------------------    
          DO i=1,n
             !The mass matrix, if needed
             !--------------------------
             IF (TransientSimulation) THEN
                MASS(i,j) = MASS(i,j)+ IP % s(t) * DetJ * &
                     Basis(i)*Basis(j)
             END IF

             !Finally, the stiffnes matrix:
             !------------------------------
             STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ * &
                  SUM(dBasisdx(i,1:DIM) * dBasisdx(j,1:DIM))
!----------------------------------------------------------------- 
          END DO ! end Loop over i-components of matrices
!----------------------------------------------------------------- 
       END DO ! end Loop over j-components of matrices and vector
!----------------------------------------------------------------- 
!-----------------------------------------------------------------
    END DO ! end Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE BoundaryCondition(LOAD, FORCE, Element, n)
    IMPLICIT NONE
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
         LocalHeatCapacity, LocalDensity
    LOGICAL :: stat, getSecondDerivatives
    INTEGER :: t,j
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    CALL GetElementNodes( Nodes )

    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
!-----------------------------------------------------------------
!   Loop over Gauss-points (boundary element Integration)
!-----------------------------------------------------------------
    DO t=1,IP % n
       !Basis function values & derivatives at the integration point:
       !-------------------------------------------------------------
       getSecondDerivatives = .FALSE.
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)
       !The source term at the integration point:
       !-----------------------------------------
       LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
!-----------------------------------------------------------------    
!      Loop over j-components of matrices and over force vector
!-----------------------------------------------------------------    
       DO j=1,n
          FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE BoundaryCondition
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE getDistance


SUBROUTINE signedDistance1( Model,Solver,dt,TransientSimulation )

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
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: DistanceSol

  LOGICAL :: AllocationsDone = .FALSE., Found, Converged

  INTEGER ::i,  n, t, istat, other_body_id, iter, NonlinearIter
  INTEGER, POINTER :: DistancePerm(:)
  REAL(KIND=dp) :: Norm, PrevNorm=0.0d00, NonlinearTol, RelativeChange, NonlinearRelax
  REAL(KIND=dp), POINTER :: Distance(:), PrevDistance(:)

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), FORCE(:), PrevDistanceElem(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType

  SAVE MASS, STIFF, LOAD, FORCE,&
       AllocationsDone, PrevNorm, PrevDistanceElem, PrevDistance
!------------------------------------------------------------------------------


  DistanceSol => Solver % Variable
  DistancePerm  => DistanceSol % Perm
  Distance => DistanceSol % Values
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     

     ALLOCATE( FORCE(N), LOAD(N), MASS(N,N), STIFF(N,N),&
          PrevDistanceElem(N), PrevDistance(SIZE(Solver % Variable % Values)),&
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'signedDistance',&
             'Memory allocation error for matrix/vectors.' )
     END IF

     PrevDistance = 1.0D-01

     AllocationsDone = .TRUE.
  END IF

  !Read in solver parameters
  !-------------------------
  SolverParams => GetSolverParams()
  IF (.NOT. ASSOCIATED(SolverParams))&
       CALL FATAL('signedDistance','No Solver section found')
  NonlinearIter = GetInteger(SolverParams, &
                     'Nonlinear System Max Iterations', Found)
  IF ( .NOT.Found ) NonlinearIter = 1
  NonlinearTol  = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) NonlinearTol = 1.0D-03
  NonlinearRelax  = GetConstReal( SolverParams, &
       'Nonlinear System Relaxation Factor',    Found )
  IF ( .NOT.Found ) &
       NonlinearRelax  = 1.0D00
     

  !loop as long non-linear system is converged or
  ! max iterations have been exceeded
  !-----------------------------------------------
  DO iter=1,NonlinearIter
     Converged = .FALSE.
     WRITE(Message,'(A,I5,A,I5)') 'Nonlinear iteration no.',iter,&
          ' of max. ', NonlinearIter
     CALL INFO('signedDistance',Message,level=1)

     !Initialize the system and do the assembly:
     !------------------------------------------
     CALL DefaultInitialize()

     ! Assembly for the domain
     !------------------------
     DO t=1,Solver % NumberOfActiveElements
     
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()

        PrevDistanceElem = PrevDistance(DistancePerm(Element % NodeIndexes))
        ! get load for force vector
        LOAD = 0.0d0
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) &
             LOAD(1:n) = GetReal( BodyForce, 'Distance Gradient', Found )      
     
        !Get element local matrix and rhs vector:
        !----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, PrevDistanceElem, Element, n, TransientSimulation)
     
        !Update global matrix and rhs vector from local matrix & vector:
        !---------------------------------------------------------------
!        IF ( TransientSimulation ) THEN
!           CALL Default1stOrderTime( MASS,STIFF,FORCE )
!        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
      
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------
!     assembly of von Neumann boundary conditions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     CALL DefaultFinishAssembly()
!
!    Dirichlet BCs:
!    --------------
     CALL DefaultDirichletBCs()

     ! Solve the system
     ! ----------------
     Norm = DefaultSolve()
     
     ! compute relative change of norm
     ! -------------------------------
     IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( 'signedDistance', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( 'signedDistance', Message, Level=4 )


     ! do we have do another round?
     ! ----------------------------
     IF ( RelativeChange < NonlinearTol ) THEN  ! no
        Converged = .TRUE.
        EXIT
     ELSE ! yes
        PrevNorm = Norm        
     END IF

     ! store actual values
     !--------------------
     DO i=1,SIZE(Solver % Variable % Values)
        PrevDistance(i) = NonlinearRelax * Distance(i)&
             + (1.0D00 - NonlinearRelax) *  PrevDistance(i)
     END DO
!------------------------------------------------------------------------------
  END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

  ! has non-linear solution converged?
  ! ----------------------------------
  IF ((.NOT.Converged) .AND. (NonlinearIter > 1)) THEN 
     WRITE( Message, * ) 'Nonlinear solution has not converged.',&
          ' Relative Change=',RelativeChange,'>',NonlinearTol
     CALL Warn('signedDistance', Message)
  ELSE
     WRITE( Message, * ) 'Nonlinear solution has converged after ',&
          iter,' steps.'
     CALL Info('signedDistance',Message,Level=1)
  END IF
     

!------------------------------------------------------------------------------

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, PrevDistanceElem, Element, n, TransientSimulation)
    IMPLICIT NONE
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:,:) :: MASS, STIFF
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD, PrevDistanceElem
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), PrevDistanceGradientAtIP(3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
         LocalHeatCapacity, LocalHeatConductivity, LocalDensity
    LOGICAL :: Stat, getSecondDerivatives
    INTEGER :: t,i,j,DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
!-----------------------------------------------------------------
!   Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
    DO t=1,IP % n
       !Basis function values & derivatives at the integration point:
       !-------------------------------------------------------------
       getSecondDerivatives = .FALSE.
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)
       !The source term at the integration point:
       !-----------------------------------------
       LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
       DO i=1,DIM
          PrevDistanceGradientAtIP(i) = SUM(dBasisdx(1:n,i)*PrevDistanceElem(1:n))
       END DO
!-----------------------------------------------------------------    
!      Loop over j-components of matrices and over force vector
!-----------------------------------------------------------------    
       DO j=1,n
          FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
!-----------------------------------------------------------------    
!         Loop over i-components of matrices
!-----------------------------------------------------------------    
          DO i=1,n
             !The mass matrix, if needed
             !--------------------------
!             IF (TransientSimulation) THEN
!                MASS(i,j) = MASS(i,j)+ IP % s(t) * DetJ * &
!                     Basis(i)*Basis(j)
!             END IF

             !Finally, the stiffnes matrix:
             !------------------------------
             STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ *&
                  Basis(j) * SUM(dBasisdx(i,1:DIM) * PrevDistanceGradientAtIP(1:DIM))
 !----------------------------------------------------------------- 
          END DO ! end Loop over i-components of matrices
!----------------------------------------------------------------- 
       END DO ! end Loop over j-components of matrices and vector
!----------------------------------------------------------------- 
!-----------------------------------------------------------------
    END DO ! end Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE signedDistance1
!------------------------------------------------------------------------------

!***********************
!***********************
SUBROUTINE signedDistance( Model,Solver,dt,TransientSimulation )

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
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: DistanceSol

  LOGICAL :: AllocationsDone = .FALSE., Found, Converged

  INTEGER ::i,  n, t, istat, other_body_id, iter, NonlinearIter
  INTEGER, POINTER :: DistancePerm(:)
  REAL(KIND=dp) :: Norm, PrevNorm=0.0d00, NonlinearTol, RelativeChange, NonlinearRelax
  REAL(KIND=dp), POINTER :: Distance(:), PrevDistance(:)

  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC, SolverParams
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), STIFF(:,:), LOAD(:), FORCE(:), PrevDistanceElem(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryType

  SAVE MASS, STIFF, LOAD, FORCE,&
       AllocationsDone, PrevNorm, PrevDistanceElem, PrevDistance
!------------------------------------------------------------------------------


  DistanceSol => Solver % Variable
  DistancePerm  => DistanceSol % Perm
  Distance => DistanceSol % Values
  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementNodes ! just big enough for elemental arrays
     

     ALLOCATE( FORCE(N), LOAD(N), MASS(N,N), STIFF(N,N),&
          PrevDistanceElem(N), PrevDistance(SIZE(Solver % Variable % Values)),&
          STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'signedDistance',&
             'Memory allocation error for matrix/vectors.' )
     END IF

     PrevDistance = 1.0D-01

     AllocationsDone = .TRUE.
  END IF

  !Read in solver parameters
  !-------------------------
  SolverParams => GetSolverParams()
  IF (.NOT. ASSOCIATED(SolverParams))&
       CALL FATAL('signedDistance','No Solver section found')
  NonlinearIter = GetInteger(SolverParams, &
                     'Nonlinear System Max Iterations', Found)
  IF ( .NOT.Found ) NonlinearIter = 1
  NonlinearTol  = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) NonlinearTol = 1.0D-03
  NonlinearRelax  = GetConstReal( SolverParams, &
       'Nonlinear System Relaxation Factor',    Found )
  IF ( .NOT.Found ) &
       NonlinearRelax  = 1.0D00
     

  !loop as long non-linear system is converged or
  ! max iterations have been exceeded
  !-----------------------------------------------
  DO iter=1,NonlinearIter
     Converged = .FALSE.
     WRITE(Message,'(A,I5,A,I5)') 'Nonlinear iteration no.',iter,&
          ' of max. ', NonlinearIter
     CALL INFO('signedDistance',Message,level=1)

     !Initialize the system and do the assembly:
     !------------------------------------------
     CALL DefaultInitialize()

     ! Assembly for the domain
     !------------------------
     DO t=1,Solver % NumberOfActiveElements
     
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()

        PrevDistanceElem = PrevDistance(DistancePerm(Element % NodeIndexes))
        ! get load for force vector
        LOAD = 0.0d0
        BodyForce => GetBodyForce()
        IF ( ASSOCIATED(BodyForce) ) &
             LOAD(1:n) = GetReal( BodyForce, 'Distance Gradient', Found )      
     
        !Get element local matrix and rhs vector:
        !----------------------------------------
        CALL LocalMatrix( MASS, STIFF, FORCE, LOAD, PrevDistanceElem, Element, n, TransientSimulation)
     
        !Update global matrix and rhs vector from local matrix & vector:
        !---------------------------------------------------------------
!        IF ( TransientSimulation ) THEN
!           CALL Default1stOrderTime( MASS,STIFF,FORCE )
!        END IF
        CALL DefaultUpdateEquations( STIFF, FORCE )
      
!------------------------------------------------------------------------------
     END DO
!------------------------------------------------------------------------------
!     assembly of von Neumann boundary conditions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     CALL DefaultFinishAssembly()
!
!    Dirichlet BCs:
!    --------------
     CALL DefaultDirichletBCs()

     ! Solve the system
     ! ----------------
     Norm = DefaultSolve()
     
     ! compute relative change of norm
     ! -------------------------------
     IF ( PrevNorm + Norm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     WRITE( Message, * ) 'Result Norm   : ',Norm
     CALL Info( 'signedDistance', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( 'signedDistance', Message, Level=4 )


     ! do we have do another round?
     ! ----------------------------
     IF ( RelativeChange < NonlinearTol ) THEN  ! no
        Converged = .TRUE.
        EXIT
     ELSE ! yes
        PrevNorm = Norm        
     END IF

     ! store actual values
     !--------------------
     DO i=1,SIZE(Solver % Variable % Values)
        PrevDistance(i) = NonlinearRelax * Distance(i)&
             + (1.0D00 - NonlinearRelax) *  PrevDistance(i)
     END DO
!------------------------------------------------------------------------------
  END DO ! of the nonlinear iteration
!------------------------------------------------------------------------------

  ! has non-linear solution converged?
  ! ----------------------------------
  IF ((.NOT.Converged) .AND. (NonlinearIter > 1)) THEN 
     WRITE( Message, * ) 'Nonlinear solution has not converged.',&
          ' Relative Change=',RelativeChange,'>',NonlinearTol
     CALL Warn('signedDistance', Message)
  ELSE
     WRITE( Message, * ) 'Nonlinear solution has converged after ',&
          iter,' steps.'
     CALL Info('signedDistance',Message,Level=1)
  END IF
     

!------------------------------------------------------------------------------

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, PrevDistanceElem, Element, n, TransientSimulation)
    IMPLICIT NONE
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:,:) :: MASS, STIFF
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD, PrevDistanceElem
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3), PrevDistanceGradientAtIP(3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
         LocalHeatCapacity, LocalHeatConductivity, LocalDensity
    LOGICAL :: Stat, getSecondDerivatives
    INTEGER :: t,i,j,DIM
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    DIM = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
!-----------------------------------------------------------------
!   Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
    DO t=1,IP % n
       !Basis function values & derivatives at the integration point:
       !-------------------------------------------------------------
       getSecondDerivatives = .FALSE.
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx, ddBasisddx, getSecondDerivatives)
       !The source term at the integration point:
       !-----------------------------------------
       LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
       DO i=1,DIM
          PrevDistanceGradientAtIP(i) = SUM(dBasisdx(1:n,i)*PrevDistanceElem(1:n))
       END DO
!-----------------------------------------------------------------    
!      Loop over j-components of matrices and over force vector
!-----------------------------------------------------------------    
       DO j=1,n
          FORCE(j) = FORCE(j) + IP % s(t) * DetJ * LoadAtIP * Basis(j)
!-----------------------------------------------------------------    
!         Loop over i-components of matrices
!-----------------------------------------------------------------    
          DO i=1,n
             !The mass matrix, if needed
             !--------------------------
!             IF (TransientSimulation) THEN
!                MASS(i,j) = MASS(i,j)+ IP % s(t) * DetJ * &
!                     Basis(i)*Basis(j)
!             END IF

             !Finally, the stiffnes matrix:
             !------------------------------
!             STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ *&
 !                 Basis(j) * SUM(dBasisdx(i,1:DIM) * PrevDistanceGradientAtIP(1:DIM))
             STIFF(i,j) = STIFF(i,j) + IP % s(t) * DetJ *&
                  Basis(j) * SQRT(ABS(SUM(dBasisdx(i,1:DIM) * PrevDistanceGradientAtIP(1:DIM))))
 !----------------------------------------------------------------- 
          END DO ! end Loop over i-components of matrices
!----------------------------------------------------------------- 
       END DO ! end Loop over j-components of matrices and vector
!----------------------------------------------------------------- 
!-----------------------------------------------------------------
    END DO ! end Loop over Gauss-points (element Integration)
!-----------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE signedDistance
!------------------------------------------------------------------------------




