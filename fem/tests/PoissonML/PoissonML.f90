SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation without constructing a global matrix!
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
  LOGICAL :: Found
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm, d_val, diag, rt, ct, eps, f
  INTEGER :: sz, i, j, k, l, m, n, nelem, nb, nd, t, Active, bActive, maxi

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, BC

  REAL(KIND=dp), POINTER :: x(:)
  REAL(KIND=dp), ALLOCATABLE :: b(:)

  INTEGER, ALLOCATABLE :: inds(:)

  ! structure listing degrees of freedom
  TYPE np_t
    REAL(KIND=dp) :: diag = 0._dp
    INTEGER :: eCount = 0
    INTEGER, ALLOCATABLE :: Elems(:)
  END TYPE np_t

  ! structure listing elements
  TYPE ed_t 
    INTEGER :: nd
    INTEGER, ALLOCATABLE :: dofIndeces(:)
    REAL(KIND=dp), ALLOCATABLE :: stiff(:,:), force(:)
  END TYPE ed_t

  TYPE(np_t), ALLOCATABLE :: np(:)
  INTEGER, ALLOCATABLE :: tp(:)
  TYPE(ed_t), ALLOCATABLE, TARGET :: ed(:)
!------------------------------------------------------------------------------

   Mesh => GetMesh()

   Active = GetNOFActive()
   bActive = GetNOFBoundaryElements()
   n = SIZE(Solver % Variable % Values)
   ALLOCATE(ed(Active+bActive), np(n))
   DO i=1,n
     ALLOCATE(np(i) % Elems(4))
   END DO

   !System assembly:
   !----------------
!$omp parallel do private(t,i,n,nb,nd,element,bodyforce,found)
   DO t=1,Active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes(Element)
     nd = GetElementNOFDOFs(Element)
     nb = GetElementNOFBDOFs(Element)

     BLOCK
       ! Per-iteration automatic arrays: STIFF/FORCE/LOAD used to be SAVEd
       ! allocatables shared via an OMP "private" clause, but privatizing an
       ! already-allocated allocatable gives each thread an UNALLOCATED copy
       ! (the allocation is not inherited) -- LocalMatrix's assumed-shape
       ! dummies then received an invalid actual argument, causing the
       ! intermittent SIGSEGV. Automatic arrays sized here are genuinely
       ! thread-local and correctly shaped for every call, no allocation or
       ! private-clause bookkeeping required.
       REAL(KIND=dp) :: STIFF(nd+nb,nd+nb), FORCE(nd+nb), LOAD(n)

       Load = 0.0_dp
       BodyForce => GetBodyForce(Element)
       IF ( ASSOCIATED(BodyForce) ) THEN
         Load(1:n) = GetReal( BodyForce, 'Source', Found )
       END IF

       !Get element local matrix and rhs vector:
       !----------------------------------------
       CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd+nb )
       CALL LCondensate( nd, nb, STIFF, FORCE )

       ed(t) % nd = nd
       ed(t) % Stiff = STIFF(1:nd,1:nd)
       ed(t) % Force = FORCE(1:nd)
     END BLOCK
     ed(t) % dofIndeces = [(i,i=1,nd)]
     nd =  GetElementDOFs(ed(t) % dofIndeces, Element)
     ed(t) % dofIndeces = Solver % Variable % Perm(ed(t) % dofIndeces)
   END DO
!$omp end parallel do


   ! BC's
!$omp parallel do private(t,i,j,n,nd,bc,found,element)
   DO t=1,GetNOFBoundaryElements()
     Element => GetBoundaryElement(t)

     n  = GetElementNOFNodes(Element)
     nd = GetElementNofDOFs(Element)
     j = t + Active

     BLOCK
       REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)

       BC => GetBC(Element)
       IF ( ASSOCIATED(BC)) THEN
         Load(1:n) = GetReal( BC, 'Flux', Found )
         CALL LocalMatrixBC(STIFF, FORCE, LOAD, Element, n, nd)
       ELSE
         FORCE(1:nd) = 0.0_dp
         STIFF(1:nd,1:nd) = 0.0_dp
       END IF

       ed(j) % nd = nd
       ed(j) % Stiff = STIFF(1:nd,1:nd)
       ed(j) % Force = FORCE(1:nd)
     END BLOCK
     ed(j) % dofIndeces = [(i,i=1,nd)]
     nd =  GetElementDOFs(ed(j) % dofIndeces, Element)
     ed(j) % dofIndeces = Solver % Variable % Perm(ed(j) % dofIndeces)
   END DO
!$omp end parallel do


   ! update dof structures, notably the dof element list
   DO t=1,SIZE(ed)
     IF ( t<=Active ) THEN
       Element => GetActiveElement(t)
     ELSE
       Element => GetBoundaryElement(t-Active)
     END IF
     nd = GetElementNOFDOFs(Element)
     DO i=1,nd
       j = ed(t) % dofIndeces(i)
       sz = np(j) % eCount
       IF ( sz>=SIZE(np(j) % Elems) ) THEN
         tp = np(j) % Elems(1:sz)
         DEALLOCATE(np(j) % Elems)
         ALLOCATE(np(j) % Elems(sz*2))
         np(j) % Elems(1:sz) = tp
       END IF
       np(j) % Elems(sz+1)=t
       np(j) % eCount = sz+1
       np(j) % diag = np(j) % diag + ed(t) % Stiff(i,i)
     END DO
   END DO


   ! Dirichlet conditions
   DO t=1,GetNOFBoundaryElements()
     Element => GetBoundaryElement(t)

     BC => GetBC(Element)
     IF (.NOT.ASSOCIATED(BC)) CYCLE

     d_val = GetCReal( BC, Solver % Variable % Name, Found )
     IF (.NOT. Found ) CYCLE

     n  = GetElementNOFNodes(Element)
     nd = GetElementNofDOFs(Element)

     inds = [(i,i=1,nd)]
     nd = GetElementDOFs(inds, Element)
     inds = Solver % Variable % Perm(inds)

     DO j=1,nd
       nelem = np(inds(j)) % eCount
       diag = np(inds(j)) % diag / nelem
       DO k=1,nelem
         m = np(inds(j)) % elems(k)
         DO l=1,SIZE(ed(m) % dofIndeces)
           IF(ed(m) % dofIndeces(l)==inds(j)) EXIT
         END DO
         IF(l>SIZE(ed(m) % dofIndeces)) STOP 'l'

         IF (j<=n) THEN
           ed(m) % force = ed(m) % force - ed(m) % stiff(:,l)*d_val
           ed(m) % force(l) = d_val * diag
         ELSE
           ed(m) % force(l) = 0.0_dp
         ENDIF
         ed(m) % stiff(:,l) = 0.0_dp
         ed(m) % stiff(l,:) = 0.0_dp
         ed(m) % stiff(l,l) = diag
       END DO
     END DO
   END DO

   x => Solver % Variable % Values
   n = SIZE(x)
   ALLOCATE(b(n))

!$omp parallel do private(i)
   DO i=1,n
     b(i) = 0
   END DO
!$omp end parallel do

   ! Assemble the global RHS from the per-element force vectors as an OpenMP
   ! reduction on the ALLOCATABLE array b. gfortran gives each thread a private
   ! allocatable copy (allocate/zero/combine/deallocate); that per-thread
   ! malloc/free churn used to intermittently surface as a SIGSEGV inside
   ! libgfortran's per-I/O setlocale during a later formatted WRITE (the
   ! post-solve ComputeChange norm print -- crash "after the last Linear System
   ! Timing line", all worker threads already parked). The root cause was NOT
   ! this reduction but libgfortran's non-thread-safe per-I/O setlocale on
   ! MinGW/UCRT (it merely amplified how often that fired); it is now
   ! neutralized once at startup by elmer_fix_numeric_locale()
   ! (fem/src/GFortranLocaleFix.c). With that fix in place the reduction is
   ! safe -- verified 2100/2100 clean at 4 threads, p:9 -- so we keep this
   ! clearer form rather than a hand-rolled atomic scatter.
!$omp parallel do private(i,j,k,nd) reduction(+:b)
   DO i=1,SIZE(ed)
     nd = SIZE(ed(i) % dofIndeces)
     DO j=1,nd
       k = ed(i) % dofIndeces(j)
       b(k) = b(k) + ed(i) % force(j)
     END DO
   END DO
!$omp end parallel do

   ct = Cputime(); rt=RealTime();

   eps  = GetCReal( GetSolverParams(), 'Linear System Convergence Tolerance', Found )
   maxi = GetInteger( GetSolverParams(), 'Linear System Max Iterations', Found )
   CALL PCG( n, x, b, eps, maxi )

   PRINT*,'Linear System Timing: ', cputime()-ct,realtime()-rt

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    ! Plain per-call local: a SAVEd/threadprivate Nodes here relied on
    ! !$omp threadprivate applying cleanly to a variable local to a
    ! CONTAINS-nested internal procedure, re-used across both threads and
    ! calls. Simpler and safe to just fetch it fresh every call instead.
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP, dx,dy,dz
    LOGICAL :: stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    ! Plain per-call local: a SAVEd/threadprivate Nodes here relied on
    ! !$omp threadprivate applying cleanly to a variable local to a
    ! CONTAINS-nested internal procedure, re-used across both threads and
    ! calls. Simpler and safe to just fetch it fresh every call instead.
    TYPE(Nodes_t) :: Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0_dp
    FORCE = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------


! Tailored local CG algo for speed testing (somewhat faster than any of the 
! library routines but not so much...)
!------------------------------------------------------------------------------
  SUBROUTINE  PCG( n, x, b, eps, maxiter )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: n, maxiter
    REAL(KIND=dp) :: x(n),b(n), eps
    REAL(KIND=dp) :: r(n), p(n), q(n), z(n)

    REAL(KIND=dp):: alpha, beta, rho, oldrho, residual, eps2, bnrm
    INTEGER :: iter
    REAL(KIND=dp), PARAMETER :: one = 1

    bnrm = rnrm(n,b)  ! bnrm = SUM(b*b)
    eps2 = eps**2*bnrm**2

    call rmv(n,x,r)         ! r = Ax
    CALL rsumb(n,r,-one,b)  ! r = b - r
    residual = rnrm(n,r)    ! residual = SUM(r*r)
    IF(residual<eps2) RETURN

    DO iter=1,maxiter
      call rprec(n,r,z)     ! z = r / diag(A)
      rho = rdot(n,r,z)     ! rho = SUM(r*z)
      
      IF (iter==1) THEN
        call rset(n,p,z)    ! p = z
      ELSE
        beta = rho/oldrho
        call rsumb(n,p,beta,z)    ! p = z + beta*p
      END IF
      oldrho = rho

      call rmv(n,p,q)             ! q = Ap
      alpha = rho/rdot(n,p,q)     ! alpha = rho/SUM(p*q)

      CALL rsuma(n,x,alpha,p)     ! x = x + alpha*p
      CALL rsuma(n,r,-alpha,q)    ! r = r - alpha*q

      residual = rnrm(n,r)        ! residual = SUM(r*r)
      IF(MOD(iter,25)==0) WRITE (*, '(I8, E11.4)') iter, SQRT(residual/bnrm**2)
      IF (residual < eps2) EXIT
    END DO

    call rmv(n,x,r)               ! r = Ax
    CALL rsumb(n,r,-one,b)        ! r = b - r
    residual = SQRT(rnrm(n,r))    ! residual = SUM(r*r)
    WRITE (*, *) iter, residual/bnrm
!------------------------------------------------------------------------------
  END SUBROUTINE PCG
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION rdot(n,p,q) RESULT(r)
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     REAL(KIND=dp) :: r
     REAL(KIND=dp), INTENT(in) :: p(:),q(:)

     INTEGER :: i

     r = 0
!$omp parallel do private(i) reduction(+:r)
     DO i=1,n
       r = r + p(i)*q(i)
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END FUNCTION rdot
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION rnrm(n,p) RESULT(r)
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     REAL(KIND=dp) :: r
     REAL(KIND=dp), INTENT(in) :: p(:)

     INTEGER :: i

     r = 0
!$omp parallel do private(i) reduction(+:r)
     DO i=1,n
       r = r + p(i)*p(i)
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END FUNCTION rnrm
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE rset(n,p,q)
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(in) :: q(:)
     INTEGER, INTENT(in) :: n
     REAL(KIND=dp), INTENT(out) :: p(:)

     INTEGER :: i

!$omp parallel do private(i)
     DO i=1,n
       p(i) = q(i)
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE rset
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE rsuma(n,p,alpha,q)
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(inout) :: p(:)
     INTEGER, INTENT(in) :: n
     REAL(KIND=dp), INTENT(in) :: q(:), alpha

     INTEGER :: i

!$omp parallel do private(i)
     DO i=1,n
       p(i) = p(i) + alpha*q(i)
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE rsuma
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE rsumb(n,p,beta,q)
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(inout) :: p(:)
     INTEGER, INTENT(in) :: n
     REAL(KIND=dp), INTENT(in) :: q(:), beta

     INTEGER :: i

!$omp parallel do private(i)
     DO i=1,n
       p(i) = q(i) + beta*p(i)
     END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE rsumb
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE rmv(n,u,v)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(in) :: u(n)
    INTEGER, INTENT(in) :: n
    REAL(KIND=dp), INTENT(inout) :: v(n)

    REAL(KIND=dp), ALLOCATABLE :: x(:)
    INTEGER :: i,nd
    INTEGER, POINTER :: inds(:)

!$omp parallel
!$omp do private(i)
     DO i=1,n
       v(i) = 0
     END DO
!$omp end do

!$omp do private(i,inds,x)
     DO i=1,SIZE(ed)
       inds => ed(i) % dofIndeces
       x = MATMUL(ed(i) % stiff, u(inds))
!$omp critical
       v(inds) = v(inds) + x
!$omp end critical
     END DO
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
  END SUBROUTINE rmv
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE rprec(n,u,v)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(in) :: u(n)
    INTEGER, INTENT(in) :: n
    REAL(KIND=dp), INTENT(out) :: v(n)

    INTEGER :: i
    REAL(KIND=dp), ALLOCATABLE :: x(:)

!!omp parallel do private(i,j,inds) shared(ed,n,u,v)
!      DO i=1,SIZE(ed)
!        inds = ed(i) % dofIndeces
!        DO j=1,size(inds)
!          v(inds(j)) = u(inds(j)) / ed(i) % stiff(j,j)
!        END DO
!      END DO
!!omp end parallel do

!$omp parallel do private(i)
    DO i=1,n
      v(i) = u(i) / np(i) % diag
    END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE rprec
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
