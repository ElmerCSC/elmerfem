!/*****************************************************************************/
! *
! * PlatiseFerroModel.F90
! *
! * Provides implementation of the inverse B-H Platise Ferromagnetic model
! * without hysteresis. Model functions are smooth and monotonic.
! *
! * References:
! *
! *  1. U. Platise, High Precision Wide Bandwidth Isolated Current Measurement
! *     in Networked Devices, Doctoral Dissertation, Jozef Stefan Institute,
! *     November 2021
! *  2. TODO: publish article regarding fast single term and n-term methods
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Author:  Uros Platise, ISOTEL Research
! *  Email:   uros@isotel.org
! *  Web:     http://isotel.org
! *
! *  Original Date:      23. October 2022
! *  Date modifications:
! *
! *****************************************************************************

FUNCTION PlatiseFerroModel( model, n, B ) RESULT( H )
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t)   :: model
  INTEGER         :: n, k
  REAL(KIND=dp)   :: B, H

  TYPE(ValueList_t), POINTER :: material
  LOGICAL         :: gotIt
  LOGICAL         :: firstVisited = .true.
  INTEGER         :: max_iter

  REAL(KIND=dp), POINTER  :: Hm_sqr_ptr(:,:), Bs_ptr(:,:)
  REAL(KIND=dp), ALLOCATABLE :: Hm_sqr(:), HHm_sqr(:), HHm_sqrt(:), Bs(:)
  REAL(KIND=dp)   :: B_ii, B_ii_, tol, tol_abs, u0

  SAVE Hm_sqr, HHm_sqr, HHm_sqrt, Bs, tol, max_iter, u0
  SAVE firstVisited

  IF (firstVisited) THEN
    firstVisited = .false.

    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
      CALL Fatal('PlatiseFerroModel', 'No material found')
    END IF

    u0 = GetConstReal( model % Constants, 'Permeability Of Vacuum', gotIt);
    IF (.NOT. gotIt) THEN
      u0 = 1.2566370614359173e-06
    END IF

    Hm_sqr_ptr => ListGetConstRealArray(material, 'PFM Dipoles Field Strength', gotIt);
    IF (.NOT. gotIt) THEN
      CALL Fatal('PlatiseFerroModel', 'PFM Dipoles Field Strength(n): undefined')
    END IF

    Bs_ptr => ListGetConstRealArray(material, 'PFM Dipoles Flux Density', gotIt);
    IF (.NOT. gotIt) THEN
      CALL Fatal('PlatiseFerroModel', 'PFM Dipoles Flux Density(n): undefined')
    END IF

    tol = GetConstReal(material, 'PFM Relative Tolerance', gotIt);
    IF (.NOT. gotIt) THEN
      tol = 1e-6
    END IF

    max_iter = GetConstReal(material, 'PFM Max Iterations', gotIt);
    IF (.NOT. gotIt) THEN
      max_iter = 1000 ! Largely depends on tol, typ it needs 5..8 iter on slopes, and <50 in saturated areas
    END IF

    IF (SIZE( Hm_sqr_ptr, 1 ) /= SIZE( Bs_ptr, 1 )) THEN
      CALL Fatal('PlatiseFerroModel', 'Array Sizes of PFM Dipoles Field Strength and Flux Density missmatch')
    END IF

    ALLOCATE( Hm_sqr(SIZE(Hm_sqr_ptr,1)), HHm_sqr(SIZE(Hm_sqr_ptr,1)), HHm_sqrt(SIZE(Hm_sqr_ptr,1)), Bs(SIZE(Bs_ptr,1)) )
    Hm_sqr   = Hm_sqr_ptr(:,1)**2
    Bs       = Bs_ptr(:,1)
  END IF

  H = 0

  IF (B /= 0) THEN
    tol_abs = ABS(tol * B)

    DO k = 1, max_iter
      HHm_sqr  = H**2 + Hm_sqr
      HHm_sqrt = SQRT(HHm_sqr)
      B_ii     = H*( SUM( Bs/HHm_sqrt ) + u0) - B
      B_ii_    = SUM( Bs*Hm_sqr*HHm_sqrt/(HHm_sqr**2) ) + u0
      H        = H - B_ii / B_ii_

      IF (IsNaN(H)) THEN
        CALL Fatal('PlatiseFerroModel', 'Inverse of Double Term Method diverged')
      END IF
      IF (ABS(B_ii) < tol_abs) THEN
        EXIT
      END IF
    END DO

    IF (k > max_iter) THEN
      CALL Fatal('PlatiseFerroModel', 'Inverse did not converge, may need to increase: PFM Max Iterations')
    END IF
  END IF

END FUNCTION PlatiseFerroModel

