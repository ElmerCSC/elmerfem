#ifdef USE_OMP4
#define VECT $OMP SIMD
#define ENDVECT $OMP END SIMD
#else
#define VECT DIR$ SIMD
#define ENDVECT
#endif

MODULE ElementBasisFunctions
    USE LocalTypes

    ! Module contains vectorized version of FE basis
    ! functions for selected elements



CONTAINS

    PURE SUBROUTINE LineNodalBasisVec(nvec, u, value)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        REAL(Kind=dp) :: c
        INTEGER :: j

        c=1D0/2D0
       ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = c*(1-u(j))
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = c*(1+u(j))
        END DO
!ENDVECT
    END SUBROUTINE LineNodalBasisVec

    PURE SUBROUTINE dLineNodalBasisVec(nvec, u, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u
        ! Variables
        REAL(Kind=dp), DIMENSION(:, :, :), INTENT(OUT) :: grad
        REAL(Kind=dp) :: c
        INTEGER :: j

        c=1D0/2D0
        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -c
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,1) = c
        END DO
!ENDVECT
    END SUBROUTINE dLineNodalBasisVec

    ! WARNING: this is not really for regular P triangle
    PURE SUBROUTINE TriangleNodalBasisVec(nvec, u, v, value)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        INTEGER :: j

        ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = 1-u(j)-v(j)
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = u(j)
        END DO
!ENDVECT
        ! Node 3
!VECT
        DO j=1,nvec
            value(j,3) = v(j)
        END DO
!ENDVECT
    END SUBROUTINE TriangleNodalBasisVec

    ! WARNING: this is not really for regular P triangle
    PURE SUBROUTINE dTriangleNodalBasisVec(nvec, u, v, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:,:), INTENT(OUT) :: grad
        INTEGER :: j

        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,1) = 1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,1) = 0D0
        END DO
!ENDVECT

        ! Second coordinate (eta)
!VECT
        DO j=1,nvec
            grad(j,1,2) = -1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,2) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,2) = 1D0
        END DO
!ENDVECT
    END SUBROUTINE dTriangleNodalBasisVec

    PURE SUBROUTINE QuadNodalBasisVec(nvec, u, v, value)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        REAL(Kind=dp) :: c
        INTEGER :: j

        c=1D0/4D0
        ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = c*(1-u(j))*(1-v(j))
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = c*(1+u(j))*(1-v(j))
        END DO
!ENDVECT
        ! Node 3
!VECT
        DO j=1,nvec
            value(j,3) = c*(1+u(j))*(1+v(j))
        END DO
!ENDVECT
        ! Node 4
!VECT
        DO j=1,nvec
            value(j,4) = c*(1-u(j))*(1+v(j))
        END DO
!ENDVECT
    END SUBROUTINE QuadNodalBasisVec

    PURE SUBROUTINE dQuadNodalBasisVec(nvec, u, v, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:,:), INTENT(OUT) :: grad
        REAL(Kind=dp) :: c
        INTEGER :: j

        c=1D0/4D0
        ! REAL(Kind=dp), INTENT(OUT) :: grad(nvec, 4, 2)
        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -c*(1-v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,1) = -grad(j,1,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,1) = c*(1+v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,1) = -grad(j,3,1)
        END DO
!ENDVECT

        ! Second coordinate (eta)
!VECT
        DO j=1,nvec
            grad(j,1,2) = -c*(1-u(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,2) = -c*(1+u(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,2) = -grad(j,2,2)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,2) = -grad(j,1,2)
        END DO
!ENDVECT
    END SUBROUTINE dQuadNodalBasisVec

    ! WARNING: this is not really for regular P tetra
    PURE SUBROUTINE TetraNodalBasisVec(nvec, u, v, w, value)
    IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        INTEGER :: j

        ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = 1-u(j)-v(j)-w(j)
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = u(j)
        END DO
!ENDVECT
        ! Node 3
!VECT
        DO j=1,nvec
            value(j,3) = v(j)
        END DO
!ENDVECT
        ! Node 4
!VECT
        DO j=1,nvec
            value(j,4) = w(j)
        END DO
!ENDVECT
    END SUBROUTINE TetraNodalBasisVec

    ! WARNING: this is not really for regular P tetra
    PURE SUBROUTINE dTetraNodalBasisVec(nvec, u, v, w, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:,:), INTENT(OUT) :: grad
        INTEGER :: j

        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,1) = 1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,1) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,1) = 0D0
        END DO
!ENDVECT

        ! Second coordinate (eta)
!VECT
        DO j=1,nvec
            grad(j,1,2) = -1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,2) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,2) = 1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,2) = 0D0
        END DO
!ENDVECT

        ! Third coordinate (zeta)
!VECT
        DO j=1,nvec
            grad(j,1,3) = -1D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,3) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,3) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,3) = 1D0
        END DO
!ENDVECT
    END SUBROUTINE dTetraNodalBasisVec

    ! WARNING: this is not really for regular P wedge
    PURE SUBROUTINE WedgeNodalBasisVec(nvec, u, v, w, value)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        REAL(Kind=dp) :: c
        INTEGER :: j

        c=1D0/2D0
        ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = c*(1-u(j)-v(j))*(1-w(j))
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = c*u(j)*(1-w(j))
        END DO
!ENDVECT
        ! Node 3
!VECT
        DO j=1,nvec
            value(j,3) = c*v(j)*(1-w(j))
        END DO
!ENDVECT
        ! Node 4
!VECT
        DO j=1,nvec
            value(j,4) = c*(1-u(j)-v(j))*(1+w(j))
        END DO
!ENDVECT
        ! Node 5
!VECT
        DO j=1,nvec
            value(j,5) = c*u(j)*(1+w(j))
        END DO
!ENDVECT
        ! Node 6
!VECT
        DO j=1,nvec
            value(j,6) = c*v(j)*(1+w(j))
        END DO
!ENDVECT
    END SUBROUTINE WedgeNodalBasisVec

    ! WARNING: this is not really for regular P wedge
    PURE SUBROUTINE dWedgeNodalBasisVec(nvec, u, v, w, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:,:), INTENT(OUT) :: grad
        REAL(Kind=dp) :: c
        INTEGER :: j

        c = 1D0/2D0
        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -c*(1-w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,1) = -grad(j,1,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,1) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,1) = -c*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,1) = -grad(j,4,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,1) = 0D0
        END DO
!ENDVECT

        ! Second coordinate (eta)
!VECT
        DO j=1,nvec
            grad(j,1,2) = -c*(1-w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,2) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,2) = -grad(j,1,2)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,2) = -c*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,2) = 0D0
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,2) = -grad(j,4,2)
        END DO
!ENDVECT

        ! Third coordinate (zeta)
!VECT
        DO j=1,nvec
            grad(j,1,3) = -c*(1-u(j)-v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,3) = -c*u(j)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,3) = -c*v(j)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,3) = -grad(j,1,3)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,3) = -grad(j,2,3)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,3) = -grad(j,3,3)
        END DO
!ENDVECT
    END SUBROUTINE dWedgeNodalBasisVec

    PURE SUBROUTINE BrickNodalBasisVec(nvec, u, v, w, value)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:), INTENT(OUT) :: value
        REAL(Kind=dp) :: c
        INTEGER :: j

        c = 1D0/8D0
        ! Node 1
!VECT
        DO j=1,nvec
            value(j,1) = c*(1-u(j))*(1-v(j))*(1-w(j))
        END DO
!ENDVECT
        ! Node 2
!VECT
        DO j=1,nvec
            value(j,2) = c*(1+u(j))*(1-v(j))*(1-w(j))
        END DO
!ENDVECT
        ! Node 3
!VECT
        DO j=1,nvec
            value(j,3) = c*(1+u(j))*(1+v(j))*(1-w(j))
        END DO
!ENDVECT
        ! Node 4
!VECT
        DO j=1,nvec
            value(j,4) = c*(1-u(j))*(1+v(j))*(1-w(j))
        END DO
!ENDVECT
        ! Node 5
!VECT
        DO j=1,nvec
            value(j,5) = c*(1-u(j))*(1-v(j))*(1+w(j))
        END DO
!ENDVECT
        ! Node 6
!VECT
        DO j=1,nvec
            value(j,6) = c*(1+u(j))*(1-v(j))*(1+w(j))
        END DO
!ENDVECT
        ! Node 7
!VECT
        DO j=1,nvec
            value(j,7) = c*(1+u(j))*(1+v(j))*(1+w(j))
        END DO
!ENDVECT
        ! Node 8
!VECT
        DO j=1,nvec
            value(j,8) = c*(1-u(j))*(1+v(j))*(1+w(j))
        END DO
!ENDVECT
    END SUBROUTINE BrickNodalBasisVec

    PURE SUBROUTINE dBrickNodalBasisVec(nvec, u, v, w, grad)
        IMPLICIT NONE

        ! Parameters
        INTEGER, INTENT(IN) :: nvec
        REAL(Kind=dp), DIMENSION(nvec), INTENT(IN) :: u,v,w
        ! Variables
        REAL(Kind=dp), DIMENSION(:,:,:), INTENT(OUT) :: grad
        REAL(Kind=dp) :: c
        INTEGER :: j

        c = 1D0/8D0
        ! First coordinate (xi)
!VECT
        DO j=1,nvec
            grad(j,1,1) = -c*(1-v(j))*(1-w(j))
        END DO
!ENDVECT

!VECT
        DO j=1,nvec
            grad(j,2,1) = -grad(j,1,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,1) = c*(1+v(j))*(1-w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,1) = -grad(j,3,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,1) = -c*(1-v(j))*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,1) = -grad(j,5,1)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,7,1) = c*(1+v(j))*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,8,1) = -grad(j,7,1)
        END DO
!ENDVECT

        ! Second coordinate (eta)
!VECT
        DO j=1,nvec
            grad(j,1,2) = -c*(1-u(j))*(1-w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,2) = -c*(1+u(j))*(1-w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,2) = -grad(j,2,2)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,2) = -grad(j,1,2)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,2) = -c*(1-u(j))*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,2) = -c*(1+u(j))*(1+w(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,7,2) = -grad(j,6,2)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,8,2) = -grad(j,5,2)
        END DO
!ENDVECT

        ! Third coordinate (zeta)
!VECT
        DO j=1,nvec
            grad(j,1,3) = -c*(1-u(j))*(1-v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,2,3) = -c*(1+u(j))*(1-v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,3,3) = -c*(1+u(j))*(1+v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,4,3) = -c*(1-u(j))*(1+v(j))
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,5,3) = -grad(j,1,3)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,6,3) = -grad(j,2,3)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,7,3) = -grad(j,3,3)
        END DO
!ENDVECT
!VECT
        DO j=1,nvec
            grad(j,8,3) = -grad(j,4,3)
        END DO
!ENDVECT
    END SUBROUTINE dBrickNodalBasisVec

END MODULE ElementBasisFunctions
