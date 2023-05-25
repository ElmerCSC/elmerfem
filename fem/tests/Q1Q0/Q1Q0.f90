SUBROUTINE StokesSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
! A Navier-Stokes solver for Q1/Q0 element
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
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found, Convect
  TYPE(Element_t),POINTER :: Element

  INTEGER :: i,j,k,l,n, nb, nd, t, istat, dim, BDOFs=1,Active, dofs
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, Pressure

  TYPE(ValueList_t), POINTER :: BodyForce, Material
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), LOAD(:,:), &
    FORCE(:), rho(:), mu(:), Velocity(:,:), MASS(:,:)

  TYPE(Variable_t), POINTER :: P

  INTEGER, ALLOCATABLE, SAVE :: Indexes(:), pCount(:)
!------------------------------------------------------------------------------

   dim = CoordinateSystemDimension()
   Mesh => GetMesh()
   DOFs = Solver % Variable % DOFs

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     n = (dim+1)*(Mesh % MaxElementDOFs+BDOFs)  ! just big enough for elemental arrays
     ALLOCATE( FORCE(n), LOAD(n,4), STIFF(n,n), MASS(n,n), &
         rho(n), mu(n), Velocity(dim+1,n), Indexes(n),STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'StokesSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  Convect = GetLogical( GetSolverParams(), 'Convect', Found )
  IF ( .NOT. Found ) Convect = .TRUE.

  Velocity = 0.0d0

  P => VariableGet( Mesh % Variables, 'aPressure' )
  IF ( .NOT. ASSOCIATED(P) .OR. DOFs>dim  ) P=>NULL()

  IF ( ASSOCIATED(P) ) THEN
    ALLOCATE( pCount(SIZE(P % Values)) )
    pCount = 0
    P % Values = 0.0_dp
  END IF

  !Initialize the system and do the assembly:
  !------------------------------------------
  Active = GetNOFActive()
  CALL DefaultInitialize()
  DO t=1,Active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes()
     nb = GetElementNOFBDOFs()
     nd = GetElementDOFs(Indexes)

     ! Volume forces:
     !---------------
     BodyForce => GetBodyForce()
     LOAD = 0.0d0
     IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1:n,1) = GetReal( BodyForce, 'Source x', Found )
        Load(1:n,2) = GetReal( BodyForce, 'Source y', Found )
        Load(1:n,3) = GetReal( BodyForce, 'Source z', Found )
     END IF

     ! Material parameters:
     !---------------------
     Material => GetMaterial()
     rho(1:n) = GetReal( Material, 'Density' )
     mu(1:n)  = GetReal( Material, 'Viscosity' )

     ! Get previous elementwise velocity iterate:
     !-------------------------------------------
     CALL GetVectorLocalSolution( Velocity )

     ! Get element local matrix and rhs vector:
     !-----------------------------------------
     CALL LocalMatrix(  MASS, STIFF, FORCE, LOAD, rho, mu, &
         Velocity, Element, n, nd, nd+nb, dim )

     IF ( TransientSimulation ) &
       CALL Default1stOrderTime( MASS, STIFF, FORCE )

     IF (n==nd) THEN
       IF ( ASSOCIATED(P) ) THEN
         k = n*dim
         Pressure = 0._dp
         DO j=1,DOFs
           DO i=1,n
             l = DOFs*(i-1)+j
             Pressure = Pressure - Velocity(j,i)*STIFF(k+1,l)
           END DO
         END DO
         Pressure = Pressure / STIFF(k+1,k+1)
         P % Values(P % Perm(Indexes(1:n))) = &
           P % Values(P % Perm(Indexes(1:n))) + Pressure
         pCount(P % Perm(Indexes(1:n))) = pCount(P % Perm(Indexes(1:n)))+1
       END IF
       CALL LCondensateP( n*DOFs,STIFF )
     END IF

     ! Update global matrix and rhs vector from local matrix & vector:
     !----------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  PrevNorm = Norm
  Norm = DefaultSolve()

  IF ( ASSOCIATED(P) ) THEN
    IF ( nd>n ) THEN
      DO t=1,GetNOFActive()
        Element => GetActiveElement(t)
        n  = GetElementNOFNodes()
        nb = GetElementNOFBDOFs()
        nd = GetElementDOFs(Indexes)

        k=DOFs*(Solver % Variable % Perm(Indexes(nd))-1)+1
        Pressure = Solver % Variable % Values(k)

        P % Values(P % Perm(Indexes(1:n))) = &
          P % Values(P % Perm(Indexes(1:n))) + Pressure
        pCount(P % Perm(Indexes(1:n))) = pCount(P % Perm(Indexes(1:n)))+1
     END DO
    END IF

    WHERE( pCount>0 ) 
      P % Values = P % Values / PCount
    END WHERE
    DEALLOCATE( pCount )
  END IF

  RELC = 2.0d0 * ABS(PrevNorm - Norm) / (PrevNorm + Norm)
  IF ( RELC < 1.0d-2 ) Newton = .TRUE.

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  MASS, STIFF, FORCE, LOAD, Nodalrho, &
     Nodalmu, NodalVelo, Element, n, nd, ntot, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: MASS(:,:), STIFF(:,:), FORCE(:), LOAD(:,:)
    REAL(KIND=dp) :: Nodalmu(:), Nodalrho(:), NodalVelo(:,:)
    INTEGER :: dim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3), &
               DetJ,LoadAtIP(dim+1),Velo(dim), Grad(dim,dim)
    REAL(KIND=dp), POINTER :: A(:,:), F(:),M(:,:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: mu = 1.0d0, rho = 1.0d0, s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    MASS  = 0.0d0
    FORCE = 0.0d0

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       ! Material parameters at the integration point:
       !----------------------------------------------
       mu  = SUM( Basis(1:n) * Nodalmu(1:n) )
       rho = SUM( Basis(1:n) * Nodalrho(1:n) )

       ! Previous velocity at the integration point:
       !--------------------------------------------
       Velo = MATMUL( NodalVelo(1:dim,1:n), Basis(1:n) )
       Grad = MATMUL( NodalVelo(1:dim,1:n), dBasisdx(1:n,1:dim) )

       ! The source term at the integration point:
       !------------------------------------------
       LoadAtIP(1:dim+1) = MATMUL( Basis(1:n), LOAD(1:n,1:dim+1) )
       IF ( Convect .AND. Newton ) THEN
         LoadAtIp(1:dim) = LoadAtIp(1:dim) + rho * MATMUL(Grad,Velo)
       END IF

       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       l = dim*n + 1
       DO p=1,n
         DO q=1,n
           i = dim * (p-1)
           j = dim * (q-1)
           A => STIFF(i+1:i+dim,j+1:j+dim)
           M => MASS(i+1:i+dim,j+1:j+dim)
           DO i=1,dim
             M(i,i) = M(i,i) + s * rho * Basis(q) * Basis(p)
             DO j = 1,dim
               A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)
               A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)
               IF ( Convect ) THEN
                 A(i,i) = A(i,i) + s * rho * Velo(j) * dBasisdx(q,j) * Basis(p)
                 IF ( Newton ) THEN
                    A(i,j) = A(i,j) + s * rho * Grad(i,j) * Basis(q) * Basis(p)
                 END IF
               END IF
             END DO
           END DO
         END DO

         j = dim*(p-1)
         DO i=1,dim
           STIFF(j+i,l) = STIFF(j+i,l) - s * dBasisdx(p,i)
           STIFF(l,j+i) = STIFF(l,j+i) + s * rho * dBasisdx(p,i)
         END DO
         i = dim*(p-1)
         F => FORCE(i+1:i+dim)
         F = F + s * LoadAtIP(1:dim) * Basis(p) 
       END DO
       STIFF(l,l) = STIFF(l,l) + s * 1.d-8
    END DO

    DO p = n+1,ntot
      DO j=2,DOFs
        i = dim*(p-1)+j
        FORCE(i)   = 0.0d0
        MASS(:,i)  = 0.0d0
        MASS(i,:)  = 0.0d0
        STIFF(i,:) = 0.0d0
        STIFF(:,i) = 0.0d0
        STIFF(i,i) = 1.0d0
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LCondensateP( n,K )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  DESCRIPTION:
!     Subroutine for condensation of pressure dof,
!     Modifies given stiffness matrix.
!
!  ARGUMENTS:
!    INTEGER :: n
!      INPUT: Sum of nodal, edge and face degrees of freedom 
!
!    REAL(Kind=dp) :: K(:,:)
!      INOUT: Local stiffness matrix
!
!******************************************************************************
!------------------------------------------------------------------------------

    INTEGER :: n
    REAL(KIND=dp) :: K(:,:), Kbl(1,n), Klb(n,1), Kbb

    Kbb = K(n+1,n+1)
    Kbl(1,1:n) = K(n+1,1:n)
    Klb(1:n,1) = K(1:n,n+1)
    K(1:n,1:n) = K(1:n,1:n) - MATMUL(Klb,Kbl)/Kbb
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensateP
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE StokesSolver
!------------------------------------------------------------------------------
