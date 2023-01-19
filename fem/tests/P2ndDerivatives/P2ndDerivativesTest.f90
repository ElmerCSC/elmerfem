SUBROUTINE TestSolver( Model,Solver,dt,TransientSimulation )
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
  TYPE(Element_t),POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, istat, active, family, qp
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

   !System assembly:
   !----------------
   Active = GetNofActive()
   qp = MAXVAL(Solver % Def_Dofs(:,1,6))
   PRINT*,'Testing 2nd derivatives with p(',i2s(qp), ')...'
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()
      CALL LocalMatrix(  Element, n, nd+nb, qp )
      WRITE(*,'(A)',ADVANCE='NO') '('//i2s(Element % Type % ElementCode)//') PASSED...'
   END DO
   WRITE(*,*) ''
   WRITE(*,*) ''
   Solver % Variable % Norm = 1

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, qp )
!------------------------------------------------------------------------------
    INTEGER :: n, nd, ef, qp
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),ddBasisddx(nd,3,3),DetJ,f(nd),diff(3)
    REAL(KIND=dp) :: dNodal(nd,nd,3), dx(nd,nd,3), ddxFromNodaldx(nd,3,3),ddiff(3,3)
    REAL(KIND=dp) :: x, y, z
    LOGICAL :: Stat
    INTEGER :: i,j,k,t,p,q, edim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    edim = Element % Type % Dimension

    IP = GaussPoints( Element )
    dNodal = 0
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
               IP % W(t), detJ, Basis, dBasisdx )

      x = sum(nodes % x(1:nd)*basis(1:nd))
      y = sum(nodes % y(1:nd)*basis(1:nd))
      z = sum(nodes % z(1:nd)*basis(1:nd))

      do p=1,nd
        do q=1,nd
          stiff(p,q) = stiff(p,q) + ip % s(t) * detJ* basis(p) * basis(q)
          dnodal(p,q,:) = dnodal(p,q,:) + ip % s(t) * detj * dBasisdx(p,:) * basis(q)
        end do
      
        SELECT CASE(edim)
        CASE(1)
          do i=0,qp
            Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * (x**i)
          end do
        CASE(2)
          do i=0,qp
            do j=0,qp
              if ( i+j>qp ) cycle
              Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * (x**i*y**j)
            end do
          end do
        CASE(3)
          do i=0,qp
            do j=0,qp
              do k=0,qp
                if ( i+j+k>qp ) cycle
                Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * (x**i*y**j*z**k)
              end do
            end do
          end do
        END SELECT
      end do
    END DO

    CALL InvertMatrix(STIFF,nd)
    f = MATMUL(STIFF,FORCE)

    DO i=1,nd
      DO j=1,3
        dx(i,:,j) = MATMUL(STIFF,dNodal(i,:,j))
      END DO
    END DO

    !Numerical integration:
    !----------------------
    STIFF = 0.0d0
    FORCE = 0.0d0
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t), detJ, Basis, dBasisdx, ddBasisddx, .TRUE. )

      x = SUM(Basis(1:nd)*Nodes % x(1:nd))
      y = SUM(Basis(1:nd)*Nodes % y(1:nd))
      z = SUM(Basis(1:nd)*Nodes % z(1:nd))

      ddxFromNodaldx = 0
      do i=1,nd
        do j=1,3
          do k=1,3
            ddxFromNodaldx(i,j,k) = sum(dx(i,1:nd,j)*dBasisdx(1:nd,k))
          end do
        end do
      end do

      IF (MAXVAL(ABS(ddxFromNodaldx-ddBasisddx))>1.d-7) STOP "ddx's don't match"
      diff=0
      ddiff=0
      SELECT CASE(edim)
      CASE(1)
        do i=0,qp
          if ( i>0 ) diff(1) = diff(1) + i*x**(i-1)
          if ( i>1 ) ddiff(1,1) = ddiff(1,1) + i*(i-1)*x**(i-2)
        end do
      CASE(2)
        do i=0,qp
          do j=0,qp
            if ( i+j>qp ) cycle
            if ( i>0 ) diff(1) = diff(1) + i*x**(i-1)*y**j
            if ( j>0 ) diff(2) = diff(2) + x**i*j*y**(j-1)

            if ( i>1 ) ddiff(1,1) = ddiff(1,1) + i*(i-1)*x**(i-2)*y**j
            if ( i>0 .and. j>0 ) ddiff(1,2) = ddiff(1,2) + i*x**(i-1)*j*y**(j-1)
            if ( j>1 ) ddiff(2,2) = ddiff(2,2) + x**i*j*(j-1)*y**(j-2)
          end do
        end do
      CASE(3)
        do i=0,qp
          do j=0,qp
            do k=0,qp
              if ( i+j+k>qp ) cycle
              if ( i>0 ) diff(1) = diff(1) + i*x**(i-1)*y**j*z**k
              if ( j>0 ) diff(2) = diff(2) + x**i*j*y**(j-1)*z**k
              if ( k>0 ) diff(3) = diff(3) + x**i*y**j*k*z**(k-1)

              if ( i>1 ) ddiff(1,1) = ddiff(1,1) + i*(i-1)*x**(i-2)*y**j*z**k
              if ( i>0 .and. j>0 ) ddiff(1,2) = ddiff(1,2) + i*x**(i-1)*j*y**(j-1)*z**k
              if ( i>0 .and. k>0 ) ddiff(1,3) = ddiff(1,3) + i*x**(i-1)*y**j*k*z**(k-1)
              if ( j>1 ) ddiff(2,2) = ddiff(2,2) + x**i*j*(j-1)*y**(j-2)*z**k
              if ( j>0 .and. k>0 ) ddiff(2,3) = ddiff(2,3) + x**i*j*y**(j-1)*k*z**(k-1)
            end do
          end do
        end do
      END SELECT

!     print*, sum(dbasisdx(1:nd,1)*f), diff(1)
!     print*, sum(dbasisdx(1:nd,2)*f), diff(2)
!     print*, sum(dbasisdx(1:nd,3)*f),  diff(3)
!     print*,'*'

!     print*,sum(ddbasisddx(1:nd,1,1)*f), ddiff(1,1), sum(ddxFromNodaldx(1:nd,1,1)*f)
!     print*,sum(ddbasisddx(1:nd,1,2)*f), ddiff(1,2), sum(ddxFromNodaldx(1:nd,1,2)*f)
!     print*,sum(ddbasisddx(1:nd,1,3)*f), ddiff(1,3), sum(ddxFromNodaldx(1:nd,1,3)*f)
!     print*,sum(ddbasisddx(1:nd,2,2)*f), ddiff(2,2), sum(ddxFromNodaldx(1:nd,2,2)*f)
!     print*,sum(ddbasisddx(1:nd,2,3)*f), ddiff(2,3), sum(ddxFromNodaldx(1:nd,2,3)*f)
!     print*,'-'

      IF(ABS( diff(1) - SUM(dBasisdx(1:nd,1)*f) )>1.d-7) STOP 'dx 1'
      IF(ABS( diff(2) - SUM(dBasisdx(1:nd,2)*f) )>1.d-7) STOP 'dx 2'
      IF(ABS( diff(3) - SUM(dBasisdx(1:nd,3)*f) )>1.d-7) STOP 'dx 3'

      IF(ABS(ddiff(1,1)-SUM(ddBasisddx(1:nd,1,1)*f))>1.d-6) STOP 'ddx 1'
      IF(ABS(ddiff(1,2)-SUM(ddBasisddx(1:nd,1,2)*f))>1.d-6) STOP 'ddx 2'
      IF(ABS(ddiff(1,3)-SUM(ddBasisddx(1:nd,1,3)*f))>1.d-6) STOP 'ddx 3'
      IF(ABS(ddiff(2,2)-SUM(ddBasisddx(1:nd,2,2)*f))>1.d-6) STOP 'ddx 4'
      IF(ABS(ddiff(2,3)-SUM(ddBasisddx(1:nd,2,3)*f))>1.d-6) STOP 'ddx 5'

      IF(ABS(ddiff(1,1)-SUM(ddxFromNodaldx(1:nd,1,1)*f))>1.d-6) STOP 'dfx 1'
      IF(ABS(ddiff(1,2)-SUM(ddxFromNodaldx(1:nd,1,2)*f))>1.d-6) STOP 'dfx 2'
      IF(ABS(ddiff(1,3)-SUM(ddxFromNodaldx(1:nd,1,3)*f))>1.d-6) STOP 'dfx 3'
      IF(ABS(ddiff(2,2)-SUM(ddxFromNodaldx(1:nd,2,2)*f))>1.d-6) STOP 'dfx 4'
      IF(ABS(ddiff(2,3)-SUM(ddxFromNodaldx(1:nd,2,3)*f))>1.d-6) STOP 'dfx 5'
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE TestSolver
!------------------------------------------------------------------------------
