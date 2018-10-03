SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation!
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
  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, istat, active, dim, NonlinIter, iter, NonLinMinIter
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)
  TYPE(Variable_t), POINTER :: mat_ip_var
  type(valuehandle_t) :: grad_phi_h, u_mat_h, alpha_h
  REAL :: trad

  SAVE STIFF, LOAD, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  ! dim = CoordinateSystemDimension()
  dim = 2

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  mat_ip_var => VariableGet(Mesh % Variables, 'mat_ip')

  call ListInitElementKeyword(grad_phi_h, 'Material', 'grad_phi', defrvalue=1.0_dp )

  call ListInitElementKeyword(u_mat_h, 'Material', 'u_mat', DefRValue=1.0_dp)
  call ListInitElementKeyword(alpha_h, 'Material', 'alpha', DefRValue=1.0_dp)

  NonlinIter = GetInteger(GetSolverParams(), 'Nonlinear System Max Iterations',Found)
  IF(.NOT.Found) NonlinIter = 1

  NonlinMinIter = GetInteger(GetSolverParams(), 'Nonlinear System Min Iterations',Found)
  IF(.NOT.Found) NonlinMinIter = 1

  trad = ListGetCReal(GetSolverParams(), 'Traditional', Found)
  IF (.NOT. found) trad = 0.0_dp


  ! Nonlinear iteration
  DO iter = 1,NonlinIter
    Active = GetNOFActive()
    CALL DefaultInitialize()

    !System assembly:
    !----------------
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      LOAD = 0.0d0
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) &
          Load(1:n) = GetReal( BodyForce, 'Source', Found )

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd+nb, u_mat_h, alpha_h, trad)
      CALL LCondensate( nd, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
    END DO
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()


    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()
    IF( Solver % Variable % NonlinConverged == 1 .and. iter >= NonLinMinIter ) EXIT

  end do

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd, umat_h, alpha_h, trad)
!------------------------------------------------------------------------------
    use ISO_C_BINDING
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    REAL :: trad
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP, phi(nd), dut(nd,dim)
    real(kind=dp), TARGET :: u_trad(2)
    real(kind=dp) :: sigma
    REAL(kind=dp), POINTER :: u(:), du(:,:)
    real(kind=dp) :: rvalue, alpha
    real(kind=dp), pointer :: udu(:,:) !use only first 6 for 2d and all 12 for 3d in udu (forget 3d for now...)
    real(kind=dp), allocatable :: sigma_t(:,:)
    LOGICAL :: Stat
    INTEGER :: i,t, j, rdim
    TYPE(GaussIntegrationPoints_t) :: IP
    Type(ValueList_t), POINTER :: Material
    type(ValueListEntry_t), POINTER :: u_ptr

    TYPE(Nodes_t) :: Nodes
    procedure(), pointer :: u_mat
    type(c_funptr) :: c_u_ptr
    type(ValueHandle_t) :: umat_h, alpha_h

    SAVE Nodes, sigma_t
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    Material => GetMaterial()
    CALL GetScalarLocalSolution(phi)

    DO t=1,IP % n

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      rvalue = ListGetElementReal(umat_h, basis=basis, Found=found, &
          GaussPoint = t, rdim=rdim, rtensor=udu, IP=IP)
      alpha = ListGetElementReal(alpha_h, basis=basis, Found=Found, &
          GaussPoint = t)
      
      if(trad < 0.5_dp) then
        u(1:2) => udu(1:2,1)
        du(1:2,1:2) => udu(3:6,1)
      else
        u_trad(1:2) = alpha * matmul(phi(1:nd), dbasisdx(1:nd, 1:2))
        u(1:2) => u_trad(1:2)
      end if

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      if(trad < 0.5_dp) then
        do i=1,nd
          do j=1,dim
            dut(i,j) = sum(du(j,1:dim)*dbasisdx(i,1:dim))
          end do
        end do
      else
        dut(1:nd,1:dim) = alpha * dbasisdx(1:nd, 1:dim)
      end if

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             matmul(dut, transpose(dbasisdx(:,1:dim)))
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * & 
          ( LoadAtIP * Basis(1:nd) - (u(1)*dbasisdx(1:nd,1) + u(2) * dbasisdx(1:nd,2)) )
    END DO
    FORCE = FORCE + MATMUL(STIFF, phi)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

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
    K(1:n,1:n) = &
         K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine PoissonSolver_Derived_gradient(Variable, Element, IP, ipind, grad_phi)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Variable_t) :: Variable
  TYPE(Element_t), INTENT(IN) :: Element
  TYPE(GaussIntegrationPoints_t), INTENT(IN) :: IP
  INTEGER, INTENT(IN) :: ipind
  REAL(KIND=dp), intent(out) :: grad_phi(:)

  integer :: nd
!-------------------------------------------------------------------------------
  integer :: t, k, ngp, dim
  integer :: elemind  = -1
  logical :: estat, recache
  REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:), dBasisdx(:,:,:), DetJ(:)
  REAL(KIND=dp), ALLOCATABLE, SAVE :: local_phi(:), grad_phi_cache(:,:)
  TYPE(Nodes_t) :: Nodes
!-------------------------------------------------------------------------------
  !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, DetJ, local_phi, elemind, grad_phi_cache)

  SAVE Nodes
  recache = .false.
  if (elemind /= Element % ElementIndex) recache = .true.
  elemind = Element % ElementIndex

  IF (recache) THEN
    call GetElementNodesVec(Nodes, Uelement=Element, &
        USolver = Variable % Solver, UMesh = Variable % Solver % Mesh)

    ngp = IP % n
    nd = GetElementNOFDOFs(Element, USolver=Variable % Solver)
    dim = CoordinateSystemDimension()
    if (dim /= size(grad_phi,1)) then
      call fatal('PoissonSolver_Derived_gradient', &
          'Wrong output size.  Expecting '// trim(i2s(dim)) // ' got ' // trim(i2s(size(grad_phi,1))) // '.')
    end if

    IF (.not. allocated(grad_phi_cache) .or. &
        size(grad_phi_cache, 2) < dim .or. &
        size(grad_phi_cache, 1) < ngp) THEN
      ALLOCATE(grad_phi_cache(ngp, dim))
    END IF

    IF (.not. ALLOCATED(Basis) .or. size(basis,1) < nd .or. size(basis,2) < ngp) THEN
      ALLOCATE(basis(ngp, nd), dbasisdx(ngp, nd, 3), detj(ngp))
    END IF

    IF(.not. ALLOCATED(local_phi) .or. size(local_phi,1) < nd) THEN
      ALLOCATE(local_phi(nd))
    END IF

    CALL GetScalarLocalSolution(local_phi, UVariable=Variable)
    grad_phi_cache = 0.0_dp
    estat = ElementInfoVec(Element, Nodes, ngp, IP % U, ip%v, ip%w,&
        detJ, size(basis,2), basis, dbasisdx)
    DO t = 1, nd
      DO k = 1, ip % n
        grad_phi_cache(k,1) = grad_phi_cache(k,1) + dbasisdx(k,t,1) * local_phi(t)
      END DO
      DO k = 1, ip % n
        grad_phi_cache(k,2) = grad_phi_cache(k,2) + dbasisdx(k,t,2) * local_phi(t)
      END DO
    END DO
  END IF

  grad_phi(:) = grad_phi_cache(ipind, :)

!-------------------------------------------------------------------------------
END subroutine
!-------------------------------------------------------------------------------

! UDFs

!-------------------------------------------------------------------------------
FUNCTION source(Model, n, t) RESULT(s)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t
  REAL(KIND=dp) :: s

  s = 0.0_dp
  IF (t<0.5) s = 1.0_dp
!-------------------------------------------------------------------------------
END FUNCTION source
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine u_mat_new(Model, n, tau_v, u_du)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: tau_v(4) ! 2 dims for tau and 2 dims for v (which corresponds to grad phi)
  REAL(KIND=dp) :: u_du(:,:) ! 2 + 2x2 
!-------------------------------------------------------------------------------

  u_du = 0.0_dp
  u_du(1:2,1) = tau_v(1)*tau_v(3:4)
  u_du(1,1) = u_du(1,1)! + tau_v(2)
  u_du(3,1) = tau_v(1)
  u_du(6,1) = tau_v(1)
!-------------------------------------------------------------------------------
END SUBROUTINE u_mat_new
!------------------------------------------------------------------------------- 


!-------------------------------------------------------------------------------
FUNCTION alpha(Model, n, X) RESULT(Y)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X(2)
  REAL(KIND=dp) :: Y
!-------------------------------------------------------------------------------

  Y = X(1)
!-------------------------------------------------------------------------------
END FUNCTION alpha
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
FUNCTION trad(Model, n, X) RESULT(Y)
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: X
  REAL(KIND=dp) :: Y
!-------------------------------------------------------------------------------
  if (x < 1.5) then
    print *, "Traditional"
    y = 1.0_dp
  else
    print *, "Non-traditional"
    y = 0.0_dp
  end if

!-------------------------------------------------------------------------------
END FUNCTION trad
!-------------------------------------------------------------------------------
