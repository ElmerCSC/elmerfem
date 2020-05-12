!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation.
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!> Version supporting multithreading and SIMD friendly ElmerSolver
!> kernels. 
!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr
  LOGICAL :: Found, VecAsm
!------------------------------------------------------------------------------
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  nthr = 1
  !$ nthr = omp_get_max_threads()

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()

    totelem = 0

    nColours = GetNOFColours(Solver)
    VecAsm = (nColours > 1) .OR. (nthr == 1)

    CALL ResetTimer('ThreadedAssembly')

    !$OMP PARALLEL &
    !$OMP SHARED(Solver, Active, nColours, &
    !$OMP        VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb,col) &
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    DO col=1,nColours

       !$OMP SINGLE
       CALL Info('ModelPDEthreaded','Assembly of colour: '//TRIM(I2S(col)),Level=10)
       Active = GetNOFActive(Solver)
       !$OMP END SINGLE

       !$OMP DO
       DO t=1,Active
          Element => GetActiveElement(t)
          totelem = totelem + 1
          n  = GetElementNOFNodes(Element)
          nd = GetElementNOFDOFs(Element)
          nb = GetElementNOFBDOFs(Element)
          CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm )
       END DO
       !$OMP END DO
    END DO
    !$OMP END PARALLEL 

    CALL CheckTimer('ThreadedAssembly',Delete=.TRUE.)
    totelem = 0

    CALL DefaultFinishBulkAssembly()

    nColours = GetNOFBoundaryColours(Solver)
    VecAsm = (nColours > 1) .OR. (nthr == 1)

    CALL ResetTimer('ThreadedBCAssembly')
    
    !$OMP PARALLEL &
    !$OMP SHARED(Active, Solver, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb, col) & 
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    DO col=1,nColours
       !$OMP SINGLE
       CALL Info('ModelPDEthreaded','Assembly of boundary colour: '//TRIM(I2S(col)),Level=10)
       Active = GetNOFBoundaryActive(Solver)
       !$OMP END SINGLE

       !$OMP DO
       DO t=1,Active
          Element => GetBoundaryElement(t)
          ! WRITE (*,*) Element % ElementIndex
          totelem = totelem + 1
          IF(ActiveBoundaryElement(Element)) THEN
             n  = GetElementNOFNodes(Element)
             nd = GetElementNOFDOFs(Element)
             nb = GetElementNOFBDOFs(Element)
             CALL LocalMatrixBC(  Element, n, nd+nb, nb, VecAsm )
          END IF
       END DO
       !$OMP END DO
    END DO
    !$OMP END PARALLEL

    CALL CheckTimer('ThreadedBCAssembly',Delete=.TRUE.)
        
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged == 1 ) EXIT

  END DO

CONTAINS

! Assembly of the matrix entries arising from the bulk elements. SIMD version.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(IN) :: VecAsm
!------------------------------------------------------------------------------

    REAL(KIND=dp), ALLOCATABLE, TARGET, SAVE :: Coeffs(:,:), CoeffsInBasis(:,:)
    INTEGER, PARAMETER :: DIFF_IND = 1, CONV_IND = 2, REACT_IND= 3, TIME_IND = 4, &
         LOAD_IND = 5, VELO_IND = 6
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJ(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER CONTIG :: diff_coeff(:), conv_coeff(:), react_coeff(:), &
         time_coeff(:), Load(:), Velo(:,:), LoadAtIp(:), &
         DiffInBasis(:), ReactInBasis(:), ConvInBasis(:), VeloInBasis(:,:), rho(:)
    REAL(KIND=dp) :: Weight
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t), SAVE :: Nodes
    !$OMP THREADPRIVATE(Coeffs, CoeffsInBasis, Basis, dBasisdx, DetJ, &
    !$OMP               MASS, STIFF, FORCE, Nodes)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJ, Coeffs, CoeffsInBasis
!DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()
    IP = GaussPoints( Element )
    ngp = IP % n

    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
       IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
            DEALLOCATE(Basis,dBasisdx, DetJ, Coeffs, CoeffsInBasis, &
            MASS, STIFF, FORCE)
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
       ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJ(ngp), &
            Coeffs(n,VELO_IND+2), CoeffsInBasis(ngp,VELO_IND+2), &
            MASS(nd,nd), STIFF(nd,nd), FORCE(nd), STAT=allocstat)
       IF (allocstat /= 0) THEN
          CALL Fatal('ModelPDEthreaded::LocalMatrix','Local storage allocation failed')
       END IF
    END IF

    CALL GetElementNodesVec( Nodes, UElement=Element )

    ! Set up pointers to coefficients and to their projections in basis
    Load => Coeffs(1:n,LOAD_IND)
    diff_coeff => Coeffs(1:n,DIFF_IND)
    react_coeff => Coeffs(1:n,REACT_IND)
    conv_coeff => Coeffs(1:n,CONV_IND)
    time_coeff => Coeffs(1:n,TIME_IND)
    Velo => Coeffs(:,VELO_IND:VELO_IND+2)
    
    LoadAtIP => CoeffsInBasis(1:ngp,LOAD_IND)
    DiffInBasis => CoeffsInBasis(1:ngp,DIFF_IND)
    ReactInBasis => CoeffsInBasis(1:ngp,REACT_IND)
    ConvInBasis => CoeffsInBasis(1:ngp,CONV_IND)
    rho => CoeffsInBasis(1:ngp, TIME_IND)
    VeloInBasis => CoeffsInBasis(:,VELO_IND:VELO_IND+2)
    
    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    Load = 0._dp

    BodyForce => GetBodyForce(Element)
    IF ( ASSOCIATED(BodyForce) ) THEN
      CALL GetRealValues( BodyForce, 'field source', Load, Found, UElement=Element )
    END IF

    Material => GetMaterial(Element)

    CALL GetRealValues(Material, 'diffusion coefficient', diff_coeff, &
         Found, UElement=Element)
    IF (.NOT. Found) diff_coeff = 0._dp

    CALL GetRealValues(Material, 'reaction coefficient', react_coeff, &
         Found, UElement=Element)
    IF (.NOT. Found) react_coeff = 0._dp

    CALL GetRealValues(Material, 'convection coefficient', conv_coeff, &
         Found, UElement=Element)
    IF (.NOT. Found) conv_coeff = 0._dp

    CALL GetRealValues(Material, 'time derivative coefficient', time_coeff, &
         Found, UElement=Element)
    IF (.NOT. Found) time_coeff = 0._dp
    
    Velo = 0._dp
    DO i=1,dim
      ! Velo(1:n,i)=GetReal(Material,'convection velocity '//TRIM(I2S(i)),Found)
      CALL GetRealValues(Material, 'convection velocity '//TRIM(I2S(i)), Coeffs(1:n,VELO_IND+i-1), &
           Found, UElement=Element)
      IF (.NOT. Found) Coeffs(1:n,VELO_IND+i-1) = 0._dp
    END DO

    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJ, &
         SIZE(Basis,2), Basis, dBasisdx )

    ! Project all coefficients to the basis. On exit, CoeffsInBasis contains
    ! projections as (D, C, R, rho, LoadAtIP, a)
    CALL LinearForms_ProjectToU(ngp, n, Basis, Coeffs, CoeffsInBasis)
    
    ! The above call is equivalent to performing the loop
    ! DO t=1,IP % n
    !    D(t) = SUM(Basis(t,1:n)*diff_coeff(1:n))
    !    C(t) = SUM(Basis(t,1:n)*conv_coeff(1:n))
    !    R(t) = SUM(Basis(,t1:n)*react_coeff(1:n))   
    !    rho(t) = SUM(Basis(t,1:n)*time_coeff(1:n))
    !    LoadAtIP(t) = SUM( Basis(t,1:n) * LOAD(1:n) )
    !    DO i=1,dim
    !       a(t,i) = SUM(Velo(1:n,i)*Basis(t,1:n))
    !    END DO
    ! END DO
    
    ! Compute actual integration weights (recycle the memory space of DetJ)
    DO t=1,ngp
       DetJ(t) = IP % s(t)*Detj(t)
    END DO

    ! diffusion term 
    ! STIFF=STIFF+(D*grad(u),grad(v))
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJ, STIFF, diff_coeff)
    ! reaction term 
    ! STIFF=STIFF+(R*u,v)
    CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, STIFF, ReactInBasis)
    ! advection term
    ! STIFF=STIFF+(C*grad(u),v)
    CALL LinearForms_GradUdotU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, Basis, DetJ, STIFF, &
         ConvInBasis, VeloInBasis)
    ! time derivative
    ! MASS=MASS+(rho*du/dt,v)
    CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, MASS, rho)

    ! FORCE=FORCE+(u,f)
    CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, LoadAtIP, FORCE)

    ! The above calls are equivalent to 
    ! DO t=1,IP % n
       ! Weight = IP % s(t) * DetJ(t)
       ! Weight = DetJ(t)

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      ! STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
      !        D(t) * MATMUL( dBasisdx(t,1:nd,1:dim), TRANSPOSE( dBasisdx(t,1:nd,1:dim ) ))
      ! DO p=1,nd
      !   DO q=1,nd
          ! advection term (C*grad(u),v)
          ! -----------------------------------
          ! STIFF (p,q) = STIFF(p,q) + Weight * &
          !    C(t) * SUM(a(t,1:dim)*dBasisdx(t,q,1:dim)) * Basis(t,p)

          ! DO d=1,dim
          !   STIFF(p,q) = STIFF(p,q) + Weight * &
          !         C(t) * a(t,d)*dBasisdx(t,q,d) * Basis(t,p)
          ! END DO
          
          ! reaction term (R*u,v)
          ! -----------------------------------
          ! STIFF(p,q) = STIFF(p,q) + Weight * R(t)*Basis(t,q) * Basis(t,p)

          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          ! MASS(p,q) = MASS(p,q) + Weight * rho(t) * Basis(t,q) * Basis(t,p)
      !  END DO
      ! END DO

      ! FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP(t) * Basis(t,1:nd)
    ! END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------

! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, VecAsm )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: VecAsm
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), RobinCoeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    CALL GetRealValues(BC,'field flux', Flux, Found, UElement=Element )
    IF (.NOT. Found) Flux = 0._dp

    CALL GetRealValues(BC,'robin coefficient', RobinCoeff, Found, UElement=Element )
    IF (.NOT. Found) RobinCoeff = 0._dp

    CALL GetRealValues(BC,'external field', Ext_t, Found, UElement=Element )
    IF (.NOT. Found) Ext_t = 0._dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = SUM(Basis(1:n)*flux(1:n))

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = SUM(Basis(1:n)*RobinCoeff(1:n))
      Ext = SUM(Basis(1:n)*Ext_t(1:n))

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

! Perform static condensation in case bubble dofs are present
!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------
