!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Mikko Byckling
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 31 May 2017
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!>  Module defining vectorized operations on common (bi)linear forms.
!-----------------------------------------------------------------------------


#include "../config.h"

MODULE LinearForms
  USE Types, ONLY: dp, VECTOR_BLOCK_LENGTH, VECTOR_SMALL_THRESH
  USE Messages
  IMPLICIT NONE
  PRIVATE

  INTERFACE LinearForms_ProjectToU
    MODULE PROCEDURE LinearForms_ProjectToU_rank1, LinearForms_ProjectToU_rankn
  END INTERFACE LinearForms_ProjectToU

  PUBLIC LinearForms_GradUdotGradU, LinearForms_UdotF, LinearForms_ProjectToU
CONTAINS

  ! Compute bilinear form G=G+(alpha grad u, grad u) = grad u .dot. (alpha grad u) 
  SUBROUTINE LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: m, n, dim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: GradU(:,:,:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: G(:,:)
    REAL(KIND=dp) CONTIG, INTENT(IN), OPTIONAL :: alpha(:)

    REAL(KIND=dp) :: wrk(VECTOR_BLOCK_LENGTH,n)
    INTEGER :: i, ii, iin, j, l, k, kk, ldbasis, ldwrk, ldk, blklen
    LOGICAL :: noAlphaWeight
!DIR$ ATTRIBUTES ALIGN:64::wrk

    ldbasis = SIZE(GradU,1)
    ldwrk = SIZE(wrk,1)
    ldk = SIZE(G,1)

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin=MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen=iin-ii+1
      
      IF (blklen < VECTOR_SMALL_THRESH) THEN
        ! Do not attempt to call BLAS for small cases to avoid preprocessing overhead
        IF (noAlphaWeight) THEN
          DO j=1,n
            !$OMP SIMD PRIVATE(l,k)
            DO i=1,n
!DIR$ LOOP COUNT MAX=3
!DIR$ UNROLL
              DO k=1,dim
                DO l=ii,iin
                  G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO j=1,n
            !$OMP SIMD PRIVATE(l,k)
            DO i=1,n
!DIR$ LOOP COUNT MAX=3
!DIR$ UNROLL
              DO k=1,dim
                DO l=ii,iin
                  G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)*alpha(l)
                END DO
              END DO
            END DO
          END DO
        END IF
      ELSE
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
      END IF
    END DO ! Vector blocks
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

    REAL(KIND=dp) :: wrk(VECTOR_BLOCK_LENGTH)
    INTEGER :: i, ii, iin, j, blklen, l
    LOGICAL :: noAlphaWeight
!DIR$ ATTRIBUTES ALIGN:64::wrk

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin = MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen= iin-ii+1
      ! Project local F to global basis

      IF (blklen < VECTOR_SMALL_THRESH) THEN
        IF (noAlphaWeight) THEN
          !$OMP SIMD PRIVATE(l)
          DO i=1,n
            DO l=ii,iin
              UdotF(i) = UdotF(i) + U(l,i)*F(l)*weight(l)
            END DO
          END DO
        ELSE
          !$OMP SIMD PRIVATE(l)
          DO i=1,n
            DO l=ii,iin
              UdotF(i) = UdotF(i) + U(l,i)*F(l)*weight(l)*alpha(l)
            END DO
          END DO
        END IF
      ELSE
        IF (noAlphaWeight) THEN
          !$OMP SIMD
          DO i=ii,iin
            wrk(i-ii+1) = weight(i)*F(i)
          END DO
        ELSE
          !$OMP SIMD 
          DO i=ii,iin
            wrk(i-ii+1) = weight(i)*F(i)*alpha(i)
          END DO
        END IF

        CALL DGEMV('T', blklen, n, &
              1D0, U(ii,1), SIZE(U,1), wrk, 1, 1D0, UdotF, 1)
      END IF
    END DO
  END SUBROUTINE LinearForms_UdotF

END MODULE LinearForms
