MODULE LinearForms
  USE Types, ONLY: dp
  USE Messages
  IMPLICIT NONE
  PRIVATE
  
  REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:), SAVE :: wrk
  !$OMP THREADPRIVATE(wrk)

  INTERFACE LinearForms_ProjectToU
    MODULE PROCEDURE LinearForms_ProjectToU_rank1, LinearForms_ProjectToU_rankn
  END INTERFACE LinearForms_ProjectToU

  PUBLIC LinearForms_GradUdotGradU, LinearForms_UdotF, LinearForms_ProjectToU
CONTAINS
  
  SUBROUTINE LinearForms_InitWorkSpace(n1, n2)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n1, n2
    INTEGER :: allocstat

    ! Reserve thread local workspace for routines within the module 
    IF (.NOT. ALLOCATED(wrk)) THEN
      ALLOCATE(wrk(n1,n2), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LinearForms::InitWorkSpace','Storage allocation failed')
      END IF
    ELSE IF (SIZE(wrk,1) < n1 .OR. SIZE(wrk,2) < n2) THEN
      DEALLOCATE(wrk)
      ALLOCATE(wrk(n1,n2), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LinearForms::InitWorkSpace','Storage allocation failed')
      END IF
    END IF
  END SUBROUTINE LinearForms_InitWorkSpace
  
  ! Compute bilinear form G=G+(alpha grad u, grad u) = grad u .dot. (alpha grad u) 
  SUBROUTINE LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: m, n, dim
    REAL(KIND=dp), INTENT(IN) :: GradU(:,:,:), weight(:)
    REAL(KIND=dp), INTENT(INOUT) :: G(:,:)
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: alpha(:)

    INTEGER :: i, j, k, ldbasis, ldwrk, ldk
    LOGICAL :: noAlphaWeight

    ! Set up workspace
    CALL LinearForms_InitWorkSpace(m, n)

    ldbasis = SIZE(GradU,1)
    ldwrk = SIZE(wrk,1)
    ldk = SIZE(G,1)

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO k=1, dim
      IF (noAlphaWeight) THEN
        DO j=1,n
          DO i=1,m
            wrk(i,j)=weight(i)*GradU(i,j,k)
          END DO
        END DO
      ELSE
        DO j=1,n
          DO i=1,m
            wrk(i,j)=weight(i)*alpha(i)*GradU(i,j,k)
          END DO
        END DO
      END IF

      ! Compute matrix \grad u \dot \grad u for dim=k
      CALL DGEMM('T', 'N', n, n, m, 1D0, GradU(1,1,k), ldbasis, &
              wrk, ldwrk, 1D0, G, ldk)
    END DO
  END SUBROUTINE LinearForms_GradUdotGradU

  SUBROUTINE LinearForms_ProjectToU_rank1(m, n, U, F, ProjectToU)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: U(:,:), F(:)
    REAL(KIND=dp), INTENT(INOUT) :: ProjectToU(:)

    CALL DGEMV('N', m, n, 1D0, U, SIZE(U,1), F, 1, 0D0, ProjectToU, 1)
  END SUBROUTINE LinearForms_ProjectToU_rank1

  SUBROUTINE LinearForms_ProjectToU_rankn(m, n, U, F, ProjectToU)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: U(:,:), F(:,:)
    REAL(KIND=dp), INTENT(INOUT) :: ProjectToU(:,:)
    
    CALL DGEMM('N', m, SIZE(F,2), n, 1D0, U, SIZE(U,1), F, SIZE(F,1), &
            0D0, ProjectToU, SIZE(ProjectToU,1))
  END SUBROUTINE LinearForms_ProjectToU_rankn

  ! Compute linear form UdotF=UdotF+(u,f) 
  SUBROUTINE LinearForms_UdotF(m, n, U, weight, F, UdotF)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp), INTENT(IN) :: U(:,:), F(:), weight(:)
    REAL(KIND=dp), INTENT(INOUT) :: UdotF(:)

    INTEGER :: i
    
    ! Set up workspace
    CALL LinearForms_InitWorkSpace(m, 1)

    !$OMP SIMD
    DO i=1,m
      wrk(i,1) = weight(i)*F(i)
    END DO
    !$OMP END SIMD

    CALL DGEMV('T', m, n, 1D0, U, SIZE(U,1), wrk(1,1), 1, 1D0, UdotF, 1)
  END SUBROUTINE LinearForms_UdotF

END MODULE LinearForms

