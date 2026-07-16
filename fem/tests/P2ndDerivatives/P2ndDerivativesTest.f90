SUBROUTINE TestSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Test for p-element 2nd derivatives.
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
  INTEGER :: n, nb, nd, t, istat, active, qp
!------------------------------------------------------------------------------

  Mesh => GetMesh()

   Active = GetNofActive()
   qp = MAXVAL(Solver % Def_Dofs(:,1,6))

   PRINT*,'Testing 2nd derivatives with p(',i2s(qp), ')...'
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      WRITE(*,'(A)',ADVANCE='NO') '('//i2s(Element % Type % ElementCode)//')...'
      CALL LocalMatrix(  Element, n, nd, qp )
      WRITE(*,'(A)',ADVANCE='NO') 'PASSED '
   END DO
   WRITE(*,*) ''
   WRITE(*,*) ''
   Solver % Variable % Norm = 1

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, qp )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, ef, qp
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),ddBasisddx(nd,3,3),DetJ,f(nd),diff(3)
    REAL(KIND=dp) :: dNodal(nd,nd,3),dx(nd,nd,3),ddxFromNodaldx(nd,3,3),ddiff(3,3)
    REAL(KIND=dp) :: x, y, z, fx, s
    LOGICAL :: Stat, sq
    INTEGER :: i,j,k,t,p,q, edim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    edim = Element % Type % Dimension

    Nodes % x(1:n) = Element % Type % NodeU
    IF(edim>1) Nodes % y(1:n) = Element % Type % NodeV
    IF(edim>2) Nodes % z(1:n) = Element % Type % NodeW

    i = Element % Type % ElementCode
    sq = i==202 .or. i==404 .or. i==808

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
              if ( .not.sq .and.  i+j>qp ) cycle
              Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * (x**i*y**j)
            end do
          end do
        CASE(3)
          do i=0,qp
            do j=0,qp
              do k=0,qp
                if ( .not.sq .and. i+j+k>qp ) cycle
                Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * (x**i*y**j*z**k)
              end do
            end do
          end do
!               Force(p) = Force(p) + ip % s(t) * detJ * basis(p) * x * x
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

      IF(GetElementFamily()/=6 .AND. GetElementFamily()/=7) THEN
        DO i=1,nd
          IF (MAXVAL(ABS(ddxFromNodaldx(i,:,:)-ddBasisddx(i,:,:)))>1.d-8) THEN
            print*,i
            print*,ddxFromNodaldx(i,1,:), ddbasisddx(i,1,:)
            print*,ddxFromNodaldx(i,2,:), ddbasisddx(i,2,:)
            print*,ddxFromNodaldx(i,3,:), ddbasisddx(i,3,:)
            STOP "ddx's don't match"
          END IF
        END DO
      END IF
      diff=0
      fx = 0
      ddiff=0
      SELECT CASE(edim)
      CASE(1)
        do i=0,qp
          fx = fx + x**i
          if ( i>0 ) diff(1) = diff(1) + i*x**(i-1)
          if ( i>1 ) ddiff(1,1) = ddiff(1,1) + i*(i-1)*x**(i-2)
        end do
      CASE(2)
        do i=0,qp
          do j=0,qp
            if ( .not.sq .and. i+j>qp ) cycle
            fx = fx + x**i*y**j
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
              if ( .not.sq .and. i+j+k>qp ) cycle
              fx = fx + x**i*y**j*z**k
              if ( i>0 ) diff(1) = diff(1) + i*x**(i-1)*y**j*z**k
              if ( j>0 ) diff(2) = diff(2) + x**i*j*y**(j-1)*z**k
              if ( k>0 ) diff(3) = diff(3) + x**i*y**j*k*z**(k-1)

              if ( i>1 ) ddiff(1,1) = ddiff(1,1) + i*(i-1)*x**(i-2)*y**j*z**k
              if ( i>0 .and. j>0 ) ddiff(1,2) = ddiff(1,2) + i*x**(i-1)*j*y**(j-1)*z**k
              if ( i>0 .and. k>0 ) ddiff(1,3) = ddiff(1,3) + i*x**(i-1)*y**j*k*z**(k-1)
              if ( j>1 ) ddiff(2,2) = ddiff(2,2) + x**i*j*(j-1)*y**(j-2)*z**k
              if ( j>0 .and. k>0 ) ddiff(2,3) = ddiff(2,3) + x**i*j*y**(j-1)*k*z**(k-1)
              if ( k>1 ) ddiff(3,3) = ddiff(3,3) + x**i*y**j*k*(k-1)*z**(k-2)
            end do
          end do
        end do
!diff=0
!ddiff=0
!fx = x*x
!diff(1) = 2*x
!ddiff(1,1) = 2
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

      CALL CheckValue(fx,SUM(Basis(1:nd)*f), 1.0d-7, 'fx')

      CALL CheckValue(diff(1),SUM(dBasisdx(1:nd,1)*f), 2.5d-7, 'dx')
      CALL CheckValue(diff(2),SUM(dBasisdx(1:nd,2)*f), 2.5d-7, 'dy')
      CALL CheckValue(diff(3),SUM(dBasisdx(1:nd,3)*f), 2.5d-7, 'dz')

      CALL CheckValue(ddiff(1,1),SUM(ddBasisddx(1:nd,1,1)*f), 2.0d-5, 'dxx')
      CALL CheckValue(ddiff(1,2),SUM(ddBasisddx(1:nd,1,2)*f), 2.0d-5, 'dxy')
      CALL CheckValue(ddiff(1,3),SUM(ddBasisddx(1:nd,1,3)*f), 2.0d-5, 'dxz')
      CALL CheckValue(ddiff(2,2),SUM(ddBasisddx(1:nd,2,2)*f), 2.0d-5, 'dyy')
      CALL CheckValue(ddiff(2,3),SUM(ddBasisddx(1:nd,2,3)*f), 2.0d-5, 'dyz')
      CALL CheckValue(ddiff(3,3),SUM(ddBasisddx(1:nd,3,3)*f), 2.0d-5, 'dzz')
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CheckValue(f1,f2,eps,str)
!------------------------------------------------------------------------------
     REAL(KIND=dp)  :: f1,f2,eps
     REAL :: scal
     CHARACTER(*) :: str

     IF( ABS(f1) < 1.d-10 ) RETURN
     scal = MAX(ABS(f1),ABS(f2))
     IF(ABS(f1-f2)>scal*eps) THEN
       PRINT*,f1,f2,ABS(f1-f2), '>', scal*eps, str
       STOP "Test failed."
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CheckValue
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE TestSolver
!------------------------------------------------------------------------------
