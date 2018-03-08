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
  LOGICAL :: Found, VecAsm, InitHandles
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
    !$OMP SHARED(Solver, Active, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb,col, InitHandles) &
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
   
    DO col=1,nColours
      
      !$OMP SINGLE
      CALL Info('ModelPDEthreaded','Assembly of colour: '//TRIM(I2S(col)),Level=10)
      Active = GetNOFActive(Solver)
      !$OMP END SINGLE

      InitHandles = .TRUE.
      !$OMP DO
      DO t=1,Active
        Element => GetActiveElement(t)
        totelem = totelem + 1
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        nb = GetElementNOFBDOFs(Element)
        CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm, InitHandles )
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
    !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    DO col=1,nColours
       !$OMP SINGLE
       CALL Info('ModelPDEthreaded','Assembly of boundary colour: '//TRIM(I2S(col)),Level=10)
       Active = GetNOFBoundaryActive(Solver)
       !$OMP END SINGLE

       InitHandles = .TRUE. 
       !$OMP DO
       DO t=1,Active
          Element => GetBoundaryElement(t)
          ! WRITE (*,*) Element % ElementIndex
          totelem = totelem + 1
          IF(ActiveBoundaryElement(Element)) THEN
             n  = GetElementNOFNodes(Element)
             nd = GetElementNOFDOFs(Element)
             nb = GetElementNOFBDOFs(Element)
             CALL LocalMatrixBC(  Element, n, nd+nb, nb, VecAsm, InitHandles )
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
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm, InitHandles )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(IN) :: VecAsm
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJ(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, POINTER CONTIG :: DiffCoeff(:), ConvCoeff(:), ReactCoeff(:), &
        TimeCoeff(:), SourceCoeff(:), Velo1Coeff(:), Velo2Coeff(:), Velo3Coeff(:)
    REAL(KIND=dp), SAVE, POINTER CONTIG :: VeloCoeff(:,:)
    REAL(KIND=dp) :: Weight
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, DiffCoeff_h, ReactCoeff_h, TimeCoeff_h, ConvCoeff_h, &
        Velo1Coeff_h, Velo2Coeff_h, Velo3Coeff_h
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, DetJ, &
    !$OMP               MASS, STIFF, FORCE, Nodes, &
    !$OMP               SourceCoeff_h, DiffCoeff_h, ReactCoeff_h, TimeCoeff_h, &
    !$OMP               ConvCoeff_h, Velo1Coeff_h, Velo2Coeff_h, Velo3Coeff_h, &
    !$OMP               SourceCoeff, DiffCoeff, ReactCoeff, TimeCoeff, &
    !$OMP               ConvCoeff, Velo1Coeff, Velo2Coeff, Velo3Coeff, VeloCoeff )
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJ
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Field Source')
      CALL ListInitElementKeyword( DiffCoeff_h,'Material','Diffusion Coefficient')
      CALL ListInitElementKeyword( ReactCoeff_h,'Material','Reaction Coefficient')
      CALL ListInitElementKeyword( TimeCoeff_h,'Material','Time Derivative Coefficient')
      CALL ListInitElementKeyword( ConvCoeff_h,'Material','Convection Coefficient')
      CALL ListInitElementKeyword( Velo1Coeff_h,'Material','Convection Velocity 1')
      CALL ListInitElementKeyword( Velo2Coeff_h,'Material','Convection Velocity 2')
      CALL ListInitElementKeyword( Velo3Coeff_h,'Material','Convection Velocity 3')
      InitHandles = .FALSE.
    END IF

    
    dim = CoordinateSystemDimension()
    IP = GaussPoints( Element )
    ngp = IP % n

    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
      IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
            DEALLOCATE(Basis,dBasisdx, DetJ, MASS, STIFF, FORCE, VeloCoeff )
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJ(ngp), &
          MASS(nd,nd), STIFF(nd,nd), FORCE(nd), VeloCoeff(ngp,3), STAT=allocstat)
      
      IF (allocstat /= 0) THEN
        CALL Fatal('ModelPDEthreaded::LocalMatrix','Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodesVec( Nodes, UElement=Element )

    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp

    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJ, &
         SIZE(Basis,2), Basis, dBasisdx )

    ! Compute actual integration weights (recycle the memory space of DetJ)
    DO t=1,ngp
      DetJ(t) = IP % s(t) * Detj(t)
    END DO

    ! diffusion term     
    ! STIFF=STIFF+(D*grad(u),grad(v))
    DiffCoeff => ListGetElementRealVec( DiffCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJ, STIFF, DiffCoeff )
    END IF

    ! reaction term 
    ! STIFF=STIFF+(R*u,v)
    ReactCoeff => ListGetElementRealVec( ReactCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN    
      CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, STIFF, ReactCoeff )
    END IF
    
    ! time derivative
    ! MASS=MASS+(rho*du/dt,v)
    TimeCoeff => ListGetElementRealVec( TimeCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, MASS, TimeCoeff )
    END IF
      
    ! FORCE=FORCE+(u,f)
    SourceCoeff => ListGetElementRealVec( SourceCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, SourceCoeff, FORCE)
    END IF

    
    ! advection term
    ! STIFF=STIFF+(C*grad(u),v)
    ConvCoeff => ListGetElementRealVec( ConvCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN    
      Velo1Coeff => ListGetElementRealVec( Velo1Coeff_h, ngp, Basis, Element, Found ) 
      IF( Found ) VeloCoeff(1:ngp,1) = Velo1Coeff(1:ngp)
      ! This does not work
      !VeloCoeff(:,1) => Velo1Coeff(:) etc.

      Velo2Coeff => ListGetElementRealVec( Velo2Coeff_h, ngp, Basis, Element, Found ) 
      IF( Found ) VeloCoeff(1:ngp,2) = Velo2Coeff(1:ngp)

      IF( dim == 3 ) THEN
        Velo3Coeff => ListGetElementRealVec( Velo3Coeff_h, ngp, Basis, Element, Found ) 
        IF( Found ) VeloCoeff(1:ngp,3) = Velo3Coeff(1:ngp)
      END IF
      
      CALL LinearForms_GradUdotU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, Basis, DetJ, STIFF, &
          ConvCoeff, VeloCoeff )
    END IF          

    
    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------


  ! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, VecAsm, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: VecAsm
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), RobinCoeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC
       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Robin_h, Ext_h

    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Field Flux')
      CALL ListInitElementKeyword( Robin_h,'Boundary Condition','Robin Coefficient')
      CALL ListInitElementKeyword( Ext_h,'Boundary Condition','External Field')
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
           

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
      F = ListGetElementReal( Flux_h, Basis, Element, Found )
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
      END IF

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = ListGetElementReal( Robin_h, Basis, Element, Found )

      IF( Found ) THEN
        Ext = ListGetElementReal( Ext_h, Basis, Element, Found )
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * C * Ext * Basis(1:nd)
      END IF
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
