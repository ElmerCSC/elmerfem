SUBROUTINE EdgeElementSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the best approximation of the vector field U = (1+z-y,1-z+x,1-x+y)
!  with respect to the L2 norm using H(curl)-conforming basis functions. 
!  Additionally, compute the relative error of the solution using
!  the energy norm corresponding to the operator I - curl curl.
!  For the basis functions obtained via the function EdgeElementInfo the exact 
!  solution should be in the FE solution space. This solver thus offers
!  a consistency check for creating discretizations based on the basis functions 
!  provided by the function EdgeElementInfo. This test can be performed basically
!  on any 3-D mesh if the element definition of the sif file is adjusted accordingly. 
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

  LOGICAL :: stat, PiolaVersion

  INTEGER, ALLOCATABLE :: Indeces(:)

  SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indeces
!------------------------------------------------------------------------------
  PiolaVersion = .TRUE.
  dim = CoordinateSystemDimension()

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(N), LOAD(6,N), STIFF(N,N), &
          Acoef(N), Indeces(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'NedelecSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  Solver % Matrix % COMPLEX = .FALSE.
  A => GetMatrix()

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

     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix( STIFF, FORCE, LOAD, Acoef, Element, n, nd+nb, dim)

     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO

  CALL DefaultDirichletBCs(PiolaCurlTransform=.TRUE.)

  Norm = DefaultSolve()  

  !-------------------------------------------------------------------
  ! Compute the error norm for the model problem considered:
  !--------------------------------------------------------------------
  Err = 0.0d0
  SolNorm = 0.0d0
  DO t=1,Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes()
     nd = GetElementDOFs( Indeces )

     Load(1,1:nd) = Solver % Variable % Values( Solver % Variable % &
          Perm(Indeces(1:nd)) )

     CALL MyComputeError(Load, Element, n, nd, dim, Err, SolNorm)
  END DO

  PRINT *, 'Error Norm= ', sqrt(Err)/sqrt(SolNorm)

CONTAINS




!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Acoef, Element, n, nd, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:,:), Acoef(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: nu
    REAL(KIND=dp) :: EBasis(nd,3), CurlEBasis(nd,3), F(3,3), G(3,3)
    REAL(KIND=dp) :: dBasisdx(n,3)
    REAL(KIND=dp) :: Basis(n), DetJ, xq, yq, zq, uq, vq, wq, sq
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np

    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0

    !-------------------------------------
    ! Numerical integration over element:
    !-------------------------------------
    IF (PiolaVersion) THEN
       IP = EdgeElementGaussPoints(Element % TYPE % ElementCode / 100)
    ELSE
       IP = GaussPoints(Element,RelOrder=1)
    END IF

    np = 0  ! Set np = n, if nodal dofs are employed; otherwise set np = 0

    DO t=1,IP % n

       IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), F, G, detJ, Basis, EBasis, CurlEBasis, ApplyPiolaTransform = .TRUE.)
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
          ! The inner product (u,v)_K
          !----------------------------
          i = np + p
          DO q = 1,nd-np
             j = np + q
             STIFF(i,j) = STIFF(i,j) + 1.0d0 * &
                  SUM( EBasis(q,1:dim) * EBasis(p,1:dim) ) * detJ * IP % s(t)
          END DO
          
          !----------------------------------------
          ! RHS corresponding to the exact solution 
          !----------------------------------------
          FORCE(i) = FORCE(i) +  1.0d0 * EBasis(p,1) * detJ * IP % s(t) + &
               1.0d0 * EBasis(p,2) * detJ * IP % s(t) + &
               1.0d0 * EBasis(p,3) * detJ * IP % s(t) + &
               (zq * EBasis(p,1) - xq *  EBasis(p,3)) * detJ * IP % s(t) + &
               (yq * EBasis(p,3) - zq *  EBasis(p,2)) * detJ * IP % s(t) + &
               (-yq * EBasis(p,1) + xq *  EBasis(p,2)) * detJ * IP % s(t) 
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE MyComputeError(LOAD, Element, n, nd, dim, EK, SolNorm)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Load(:,:), EK, SolNorm
    TYPE(Element_t), POINTER :: Element    
    INTEGER :: n, nd, dim
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
    IF (PiolaVersion) THEN
       IP = EdgeElementGaussPoints(Element % TYPE % ElementCode / 100)
    ELSE
       IP = GaussPoints(Element,RelOrder=1)
    END IF

    np = 0  ! Set np = n, if nodal dofs are employed; otherwise set np = 0

    DO t=1,IP % n

       IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), F, G, detJ, Basis, EBasis, CurlEBasis, ApplyPiolaTransform = .TRUE.)
       ELSE
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element, EBasis, CurlEBasis, Basis, dBasisdx)
       END IF

       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
       zq = SUM( Nodes % z(1:n) * Basis(1:n) )

       u(1) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,1) )
       u(2) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,2) )
       u(3) = SUM( Load(1,np+1:nd) * EBasis(1:nd-np,3) )

       rotu(1) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,1) )
       rotu(2) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,2) )
       rotu(3) = SUM( Load(1,np+1:nd) * CurlEBasis(1:nd-np,3) )       

       ! Compute the square of the energy norm of the solution and error:
       sol = 0.0d0
       sol(1) = 1.0d0 + zq - yq
       sol(2) = 1.0d0 - zq + xq
       sol(3) = 1.0d0 - xq + yq
       rotsol(1:3) = 2.0d0
       e(:) = sol(:) - u(:)  
       rote(:) = rotsol(:) - rotu(:)

       IF (.TRUE.) THEN
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
