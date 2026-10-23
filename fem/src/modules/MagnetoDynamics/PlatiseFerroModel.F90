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
! *  Date modifications: 18. October 2026 Added: Relative Permeability and
! *                                       Support for Multiple Materials
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
  INTEGER         :: mat_id, ndip, max_iter, i

  ! Parameters are per material: one slot per material id, filled on init
  TYPE PFM_Data_t
     REAL(KIND=dp), POINTER :: Hm_sqr(:) => NULL(), Bs(:) => NULL()
     REAL(KIND=dp)          :: tol, u
     INTEGER                :: max_iter
     LOGICAL                :: ready = .FALSE.
  END TYPE PFM_Data_t
  TYPE(PFM_Data_t), ALLOCATABLE, SAVE :: PFM(:)

  REAL(KIND=dp), POINTER           :: Hm_ptr(:,:), Bs_ptr(:,:)
  REAL(KIND=dp), POINTER           :: Hm_sqr(:), Bs(:)
  REAL(KIND=dp)                    :: B_ii, B_ii_, tol, tol_abs, u_r, u, Hm
  REAL(KIND=dp), ALLOCATABLE       :: HHm_sqr(:), HHm_sqrt(:)
  CHARACTER(LEN=MAX_NAME_LEN)      :: matname


  mat_id = GetMaterialId()
  IF (mat_id < 1) THEN
    CALL Fatal('PlatiseFerroModel', 'No material associated with the current element')
  END IF

  IF (.NOT. ALLOCATED(PFM)) THEN
    ALLOCATE(PFM(model % NumberOfMaterials))
  END IF
  IF (mat_id > SIZE(PFM)) THEN
    CALL Fatal('PlatiseFerroModel', 'Material index exceeds the number of materials')
  END IF

  IF (.NOT. PFM(mat_id) % ready) THEN
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
      CALL Fatal('PlatiseFerroModel', 'No material found')
    END IF
    matname = ListGetString(material,'Name',gotIt)

    u = GetConstReal( model % Constants, 'Permeability Of Vacuum', gotIt)
    IF (.NOT. gotIt) u = 1.2566370614359173d-06

    u_r = GetConstReal(material, 'Relative Permeability', gotIt)
    IF (gotIt) u = u_r * u

    Hm_ptr => ListGetConstRealArray(material, 'PFM Dipoles Field Strength', gotIt)
    IF (.NOT. gotIt) THEN
      CALL Fatal('PlatiseFerroModel', 'PFM Dipoles Field Strength(n): undefined')
    END IF

    Bs_ptr => ListGetConstRealArray(material, 'PFM Dipoles Flux Density', gotIt)
    IF (.NOT. gotIt) THEN
      CALL Fatal('PlatiseFerroModel', 'PFM Dipoles Flux Density(n): undefined')
    END IF

    tol = GetConstReal(material, 'PFM Relative Tolerance', gotIt)
    IF (.NOT. gotIt) tol = 1d-6

    max_iter = NINT(GetConstReal(material, 'PFM Max Iterations', gotIt))
    IF (.NOT. gotIt) max_iter = 1000 ! typ. 5..8 iter on slopes, <50 in saturated areas

    IF (SIZE( Hm_ptr, 1 ) /= SIZE( Bs_ptr, 1 )) THEN
      CALL Fatal('PlatiseFerroModel', 'Array Sizes of PFM Dipoles Field Strength and Flux Density missmatch')
    END IF

    ndip = SIZE(Hm_ptr,1)

    IF (ANY(Hm_ptr(:,1) <= 0) .OR. ANY(Bs_ptr(:,1) < 0)) THEN
      CALL Fatal('PlatiseFerroModel', &
          'PFM dipoles require Hm > 0 and Bs >= 0 in material "'//TRIM(matname)//'"')
    END IF

    ALLOCATE(PFM(mat_id) % Hm_sqr(ndip), PFM(mat_id) % Bs(ndip))
    PFM(mat_id) % Hm_sqr   = Hm_ptr(:,1)**2
    PFM(mat_id) % Bs       = Bs_ptr(:,1)
    PFM(mat_id) % u        = u
    PFM(mat_id) % tol      = tol
    PFM(mat_id) % max_iter = max_iter
    PFM(mat_id) % ready    = .TRUE.

    WRITE(Message,'(A,I0,A,A,A,I0,A,ES12.4)') 'Material ', mat_id, &
        ' "', TRIM(matname), '":'//ACHAR(9)//'dipoles=', ndip, ', u=', PFM(mat_id) % u
    CALL Info('PlatiseFerroModel',Message,Level=5)

    DO i = 1, ndip
      Hm = SQRT(PFM(mat_id) % Hm_sqr(i))
      WRITE(Message,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') '  dipole ', i, &
          ': Hm=', Hm, ', Bs=', PFM(mat_id) % Bs(i), &
          ', Bs/Hm=', PFM(mat_id) % Bs(i) / Hm
      CALL Info('PlatiseFerroModel', Message, Level=5)
    END DO
  END IF

  Hm_sqr   => PFM(mat_id) % Hm_sqr
  Bs       => PFM(mat_id) % Bs
  u        =  PFM(mat_id) % u
  tol      =  PFM(mat_id) % tol
  max_iter =  PFM(mat_id) % max_iter

  H = 0

  IF (B /= 0) THEN
    tol_abs = ABS(tol * B)

    ! Inverse of the n-term Platise model, solved by Newton on B:
    !   B(H)  = H * ( SUM(Bs_i/SQRT(Hm_i**2 + H**2)) + u )
    !   B'(H) = SUM(Bs_i*Hm_i**2/(Hm_i**2 + H**2)**1.5) + u
    !
    ! with Hm_sqr(i) = Hm_i**2, Bs(i) = Bs_i and u = mu_r * mu_0

    DO k = 1, max_iter
      HHm_sqr  = H**2 + Hm_sqr
      HHm_sqrt = SQRT(HHm_sqr)
      B_ii     = H*( SUM( Bs/HHm_sqrt ) + u) - B
      B_ii_    = SUM( Bs*Hm_sqr*HHm_sqrt/(HHm_sqr**2) ) + u
      H        = H - B_ii / B_ii_

      IF (IsNaN(H)) THEN
        CALL Fatal('PlatiseFerroModel', 'Inverse method diverged')
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
