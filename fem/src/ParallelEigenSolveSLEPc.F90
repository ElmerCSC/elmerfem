!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/

! ******************************************************************************
! *
! *  Authors: Sami Ilvonen
! *  Email:   sami.ilvonen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 24 Aug 2026
! *
! *****************************************************************************/

#include <slepc/finclude/slepceps.h>
MODULE ParallelEigenSolveSLEPc

   USE Messages
   IMPLICIT NONE

CONTAINS

   SUBROUTINE SLEPcEigenSolve(Matrix, A, N, NEIG, EigValues, EigVectors)
      USE Types
      USE slepceps
      IMPLICIT NONE

      TYPE(Matrix_t), POINTER :: Matrix, A
      INTEGER :: N, NEIG
      INTEGER, ALLOCATABLE :: dperm(:)
      COMPLEX(KIND=dp) :: EigValues(:), EigVectors(:, :)
      PetscErrorCode :: ierr

      PetscCallA(SLEPcInitialize(PETSC_NULL_CHARACTER, ierr))

      PetscCallA(SlepcFinalize(ierr))

   END SUBROUTINE SLEPcEigenSolve

END MODULE ParallelEigenSolveSLEPc
