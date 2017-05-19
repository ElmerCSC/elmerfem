MODULE LinearForms
  USE Types, ONLY: dp, VECTOR_BLOCK_LENGTH
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

  SUBROUTINE LinearForms_FreeWorkSpace()
    IMPLICIT NONE

    IF (ALLOCATED(wrk)) DEALLOCATE(wrk)
  END SUBROUTINE LinearForms_FreeWorkSpace
  
  ! Compute bilinear form G=G+(alpha grad u, grad u) = grad u .dot. (alpha grad u) 
  SUBROUTINE LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n, dim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: GradU(:,:,:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: G(:,:)
    REAL(KIND=dp) CONTIG, INTENT(IN), OPTIONAL :: alpha(:)

    INTEGER :: i, ii, iin, j, k, kk, ldbasis, ldwrk, ldk, blklen
    LOGICAL :: noAlphaWeight

    ! Set up workspace
    CALL LinearForms_InitWorkSpace(VECTOR_BLOCK_LENGTH, n)

    ldbasis = SIZE(GradU,1)
    ldwrk = SIZE(wrk,1)
    ldk = SIZE(G,1)

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin=MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen=iin-ii+1
      DO k=1, dim
        IF (noAlphaWeight) THEN
          DO j=1,n
            !$OMP SIMD
            DO i=ii,iin
              wrk(i-ii+1,j)=weight(i)*GradU(i,j,k)
            END DO
          END DO
        ELSE
          DO j=1,n
            !$OMP SIMD
            DO i=ii,iin
              wrk(i-ii+1,j)=weight(i)*alpha(i)*GradU(i,j,k)
            END DO
          END DO
        END IF

        ! Compute matrix \grad u \dot \grad u for dim=k
        CALL DGEMM('T', 'N', n, n, blklen, &
                1D0, GradU(ii,1,k), ldbasis, &
                wrk, ldwrk, 1D0, G, ldk)
      END DO
    END DO
  END SUBROUTINE LinearForms_GradUdotGradU

  SUBROUTINE LinearForms_ProjectToU_rank1(m, n, U, F, ProjectToU)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp) CONTIG, INTENT(IN) :: U(:,:), F(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: ProjectToU(:)

    CALL DGEMV('N', m, n, 1D0, U, SIZE(U,1), F, 1, 0D0, ProjectToU, 1)
  END SUBROUTINE LinearForms_ProjectToU_rank1

  SUBROUTINE LinearForms_ProjectToU_rankn(m, n, U, F, ProjectToU)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp) CONTIG, INTENT(IN) :: U(:,:), F(:,:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: ProjectToU(:,:)
    
    CALL DGEMM('N', m, SIZE(F,2), n, 1D0, U, SIZE(U,1), F, SIZE(F,1), &
            0D0, ProjectToU, SIZE(ProjectToU,1))
  END SUBROUTINE LinearForms_ProjectToU_rankn

  ! Compute linear form UdotF=UdotF+(u,f) 
  SUBROUTINE LinearForms_UdotF(m, n, U, weight, F, UdotF, alpha)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp) CONTIG, INTENT(IN) :: U(:,:), F(:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: UdotF(:)
    REAL(KIND=dp) CONTIG, INTENT(IN), OPTIONAL :: alpha(:)

    INTEGER :: i, ii, iin, j, blklen
    LOGICAL :: noAlphaWeight

    ! Set up workspace
    CALL LinearForms_InitWorkSpace(VECTOR_BLOCK_LENGTH, 2)

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin = MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen= iin-ii+1
      ! Project local F to global basis

      IF (noAlphaWeight) THEN
        !$OMP SIMD
        DO i=ii,iin
          wrk(i-ii+1,2) = weight(i)*F(i)
        END DO
      ELSE
        !$OMP SIMD 
        DO i=ii,iin
          wrk(i-ii+1,2) = weight(i)*F(i)*alpha(i)
        END DO
      END IF

      CALL DGEMV('T', blklen, n, &
              1D0, U(ii,1), SIZE(U,1), wrk(1,2), 1, 1D0, UdotF, 1)
    END DO
  END SUBROUTINE LinearForms_UdotF

END MODULE LinearForms
