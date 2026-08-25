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
!
!/******************************************************************************
! *
! *  Authors: Juhani Kataja
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

MODULE ADIOS2Utils
USE Types
use GeneralUtils
USE SParIterGlobals
USE SParIterComm
USE ADIOS2
IMPLICIT NONE

INTEGER, PARAMETER :: ADIOS2_ARRAY_GLOBAL = 1
INTEGER, PARAMETER :: ADIOS2_ARRAY_LOCAL = 2
INTEGER, PARAMETER :: ADIOS2_MAX_VARNAME_LEN = 512

! Global array support added
! NOTE: Global arrays are catenated in first dimension across ranks. 
! NOTE: Local arrays are named part_#/<varname> where # is MPI rank


TYPE :: AdiosWriter_t
  TYPE(adios2_adios), PRIVATE :: adios
  TYPE(adios2_io), PRIVATE :: io
  TYPE(adios2_engine), PRIVATE :: engine
  INTEGER(kind=4) :: array_kind

  INTEGER(kind=8), PRIVATE :: step_num

  LOGICAL, PRIVATE :: instep
  LOGICAL, PRIVATE :: Finalized = .false., initialized = .false.
  LOGICAL, PRIVATE :: write_offsets = .true.
  CONTAINS

  PROCEDURE, PUBLIC :: init => init_adios_t
  PROCEDURE, PUBLIC :: finalize => finalize_adios_t

  PROCEDURE, PRIVATE :: writer_real_t, writer_integer_t, writer_real_t_2
  PROCEDURE, PRIVATE :: get_adios_shape_n
  PROCEDURE, PRIVATE :: make_varname
  PROCEDURE, PRIVATE :: define_attribute_c, define_attribute_i8
  PROCEDURE :: begin_step, end_step


  GENERIC, PUBLIC :: define_attribute => define_attribute_c, define_attribute_i8
  GENERIC, PUBLIC :: write_data => writer_integer_t, writer_real_t, writer_real_t_2
  FINAL :: finalize_sub
END TYPE AdiosWriter_t

TYPE, extends(AdiosWriter_t) :: AdiosReader_t

  CONTAINS
  PROCEDURE, PRIVATE :: read_real_t, init_adios_r_t, finalize_adios_r_t

  PROCEDURE, PUBLIC :: init => init_adios_r_t
  PROCEDURE, PUBLIC :: finalize => finalize_adios_r_t

  GENERIC, PUBLIC :: read_data => read_real_t
  FINAL :: finalize_sub_r

END TYPE AdiosReader_t

CONTAINS

!----------------------------------------
! Reader specific subroutines

SUBROUTINE finalize_sub_r(this)
  IMPLICIT NONE
  TYPE(AdiosReader_t) :: this
  integer :: ierr
  ierr = this % finalize()
END SUBROUTINE

FUNCTION finalize_adios_r_t(this) result(ierr)
  IMPLICIT NONE
  CLASS(AdiosReader_t) :: this
  INTEGER :: ierr

  IF ((.NOT. this % finalized) .AND. this % initialized) THEN
    CALL adios2_close(this%engine, ierr)
    CALL adios2_finalize(this%adios, ierr)
  END IF
  this % finalized = .true.
  this % initialized = .false.

END FUNCTION finalize_adios_r_t

FUNCTION init_adios_r_t(this, fname, array_kind, write_offsets) result(ierr)
  IMPLICIT NONE
  CLASS(AdiosReader_t) :: this
  INTEGER :: ierr
  INTEGER, OPTIONAL :: array_kind
  LOGICAL, OPTIONAL :: write_offsets
  CHARACTER(*), INTENT(IN) :: fname

  if (present(array_kind)) then
    this % array_kind = array_kind
  else
    this % array_kind = ADIOS2_ARRAY_GLOBAL
  end if

  this % instep = .false.

  call adios2_init(this % adios, parenv % activecomm, ierr)
  call adios2_declare_io(this % io, this % adios, "ioReader", ierr)
  this % finalized = .false.
  call adios2_open(this % engine, this % io, fname, adios2_mode_read, ierr)
  IF (ierr .ne. 0) THEN
      CALL Warn("init_adios_r_t", "Failed to open stream " // fname)
      this % initialized =.false.
      RETURN
  END IF
  this % initialized = .true.
END FUNCTION init_adios_r_t

SUBROUTINE read_real_t(this, varname, X)
  IMPLICIT NONE
!---------
  CLASS(AdiosReader_t) :: this
  CHARACTER(*), INTENT(IN) :: varname
  REAL(KIND=dp), INTENT(OUT) :: X(:)
!---------
  TYPE(adios2_variable) :: var
  INTEGER(KIND=8), dimension(1) :: start_dim, count_dim
  INTEGER :: ierr, local_count, global_offset
  CHARACTER(ADIOS2_MAX_VARNAME_LEN) :: adios_varname

  adios_varname = trim(varname)

  IF (.NOT. this % initialized) THEN
    Call Warn("AdiosReader_t % read_real_t", "Reader not initialized")
    RETURN
  END IF

  call this % begin_step()
  call adios2_inquire_variable(var, this % io, adios_varname, ierr)

  IF(.NOT. var % valid) THEN
    Call Warn("AdiosReader_t % read_real_t", "Variable not found: "//trim(adios_varname))
    RETURN
  END IF

  SELECT CASE (this % array_kind)
    CASE (ADIOS2_ARRAY_GLOBAL)
      ! X must be pre-allocated to the correct local size same as done in write step
      ! Per rank offset is computed to match one computed by the writer
      local_count = size(X,1)
      global_offset = 0
      CALL MPI_Exscan(local_count, global_offset, 1, MPI_INTEGER, MPI_SUM, parenv%activecomm, ierr)
      start_dim(1) = global_offset
      count_dim(1) = local_count
      CALL adios2_set_selection(var, 1, start_dim, count_dim, ierr)
      CALL adios2_get(this%engine, var, X, adios2_mode_sync, ierr)

    CASE (ADIOS2_ARRAY_LOCAL)
    ! Not implemented yet
    CASE DEFAULT
  END SELECT

  call this % end_step()

END SUBROUTINE read_real_t

!----------------------------------------
! Writer specific subroutines

SUBROUTINE begin_step(this)

  class(AdiosWriter_t) :: this
  integer :: ierr

  if ((.not. this % instep) .and. this % initialized) then
    call adios2_begin_step(this%engine, ierr)
    this % instep = .true.
  end if

END SUBROUTINE

SUBROUTINE end_step(this)

  class(AdiosWriter_t) :: this
  integer :: ierr

  if (this % instep .and. this % initialized) then

    call adios2_end_step(this%engine, ierr)
    this % instep = .false.
end if

END SUBROUTINE

SUBROUTINE define_attribute_i8(this, attr_name, data)
  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), INTENT(IN) :: attr_name
  integer(kind=8) :: data
  type(adios2_attribute) :: attribute
  INTEGER :: ierr
  call adios2_define_attribute(attribute, this%io, attr_name, data, ierr)
END SUBROUTINE

SUBROUTINE define_attribute_c(this, attr_name, data)
  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), INTENT(IN) :: attr_name
  CHARACTER(*), intent(in) :: data
  type(adios2_attribute) :: attribute
  INTEGER :: ierr
  call adios2_define_attribute(attribute, this%io, attr_name, data, ierr)
END SUBROUTINE

SUBROUTINE get_adios_shape_n(this, array_dims, shape_dims, start_dims, count_dims, varname)
  IMPLICIT NONE
  CLASS(AdiosWriter_t) :: this
  INTEGER(kind=4), dimension(:), intent(in) :: array_dims
  INTEGER(KIND=8), dimension(:), intent(out) :: shape_dims, start_dims, count_dims
  integer(kind=8), allocatable :: sum_dims(:)
  integer :: ierr
  CHARACTER(ADIOS2_MAX_VARNAME_LEN), INTENT(IN) :: varname

  IF (this % array_kind .eq. ADIOS2_ARRAY_LOCAL) THEN
    shape_dims(1) = array_dims(1)
    start_dims(1) = 0
    count_dims(1) = array_dims(1)
  ELSEIF(this%array_kind .eq. ADIOS2_ARRAY_GLOBAL) THEN
    allocate(sum_dims(size(array_dims)))
    sum_dims(:) = array_dims(:)
    CALL MPI_AllReduce(array_dims, sum_dims, 1, MPI_INTEGER, MPI_SUM, parenv % activecomm, ierr)
    shape_dims(:) = sum_dims(:)
    sum_dims(:) = 0
    call mpi_exscan(array_dims, sum_dims, 1, MPI_INTEGER, MPI_SUM, parenv % activecomm, ierr)
    start_dims(:) = sum_dims(:)
    count_dims(:) = array_dims(:)

    IF (this % write_offsets) THEN
      BLOCK
        integer(kind=8), dimension(1) :: block_shape_dims, block_start_dims, block_count_dims
        type(adios2_variable) :: var
        block_shape_dims(1) = parenv%PEs
        block_start_dims(1) = parenv%mype
        block_count_dims(1) = 1
        CALL adios2_define_variable(var, this%io, trim(varname)//"_offsets", &
          adios2_type_integer8, 1, &
          block_shape_dims, block_start_dims, block_count_dims, &
          adios2_constant_dims, ierr)
        CALL adios2_put(this%engine, var, start_dims, ierr)
      END BLOCK
    END IF

  ELSE
    CALL Fatal('AdiosWriter_t', 'Unknown array_kind')
  END IF

END SUBROUTINE

SUBROUTINE make_varname(this, varname, adios_varname)

  IMPLICIT NONE

  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), intent(in) :: varname
  CHARACTER(ADIOS2_MAX_VARNAME_LEN), intent(out) :: adios_varname

  if (this % array_kind .eq. ADIOS2_ARRAY_LOCAL) THEN
    adios_varname = "part_" // i2s(ParEnv % MyPE) // "/" // trim(varname)
  else
    adios_varname = trim(varname)
  end if
END SUBROUTINE

FUNCTION init_adios_t(this, fname, array_kind, write_offsets) result(ierr)
  IMPLICIT NONE
  CLASS(AdiosWriter_t) :: this
  INTEGER :: ierr
  CHARACTER(*), intent(in) :: fname
  INTEGER, OPTIONAL :: array_kind
  LOGICAL, OPTIONAL :: write_offsets
  INTEGER :: array_kind_


  if (present(array_kind)) then
    this % array_kind = array_kind
  else
    this % array_kind = ADIOS2_ARRAY_GLOBAL
  end if

  this % step_num = 1
  this % instep = .false.

  if (present(write_offsets)) this % write_offsets = write_offsets

  CALL adios2_init(this % adios, parenv % activecomm, ierr)
  CALL adios2_declare_io(this % io, this % adios, "ioWriter", ierr)
  CALL adios2_open(this % engine, this % io, fname, adios2_mode_write, ierr)
  this % finalized = .false.
  this % initialized = .true.

END FUNCTION init_adios_t

FUNCTION finalize_adios_t(this) result(ierr)
  IMPLICIT NONE
  CLASS(AdiosWriter_t) :: this
  INTEGER :: ierr
  IF ((.NOT. this % finalized) .AND. this % initialized) THEN
    call adios2_flush_all(this%adios, ierr)
    CALL adios2_close(this%engine, ierr)
    CALL adios2_finalize(this%adios, ierr)
  END IF
  this % finalized = .true.
  this % initialized = .false.

END FUNCTION finalize_adios_t

SUBROUTINE finalize_sub(this) 
  IMPLICIT NONE
  TYPE(AdiosWriter_t) :: this
  INTEGER :: ierr
  ierr = this%finalize()
END SUBROUTINE finalize_sub


SUBROUTINE writer_integer_t(this, varname, x)

  IMPLICIT NONE

  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), intent(in) :: varname
  INTEGER(KIND=4), intent(in), dimension(:) :: x

  INTEGER(KIND=8), dimension(1) :: shape_dims, start_dims, count_dims
  INTEGER :: ierr
  CHARACTER(ADIOS2_MAX_VARNAME_LEN) :: adios_varname ! TODO: declare parameter max_adios_varname or use automatic allocation here
  TYPE(adios2_variable) :: var

  call this % make_varname(varname, adios_varname)

  call adios2_inquire_variable(var, this%io, adios_varname, ierr)
  if(.not. var % valid) then
    CALL this % get_adios_shape_n(shape(x), shape_dims, start_dims, count_dims, adios_varname)
    CALL adios2_define_variable(var, this%io, adios_varname, adios2_type_integer4, 1, &
      shape_dims, start_dims, count_dims, adios2_constant_dims, ierr)
  end if

  CALL adios2_put(this%engine, var, x, ierr)

END SUBROUTINE writer_integer_t

SUBROUTINE writer_real_t_2(this, varname, x)

  IMPLICIT NONE

  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), intent(in) :: varname
  REAL(KIND=dp), intent(in), dimension(:,:) :: x

  INTEGER(KIND=8), dimension(2) :: shape_dims, start_dims, count_dims
  INTEGER :: ierr
  CHARACTER(ADIOS2_MAX_VARNAME_LEN) :: adios_varname
  TYPE(adios2_variable) :: var

  call this % make_varname(varname, adios_varname)

  call adios2_inquire_variable(var, this%io, adios_varname, ierr)
  if(.not. var % valid) then
    CALL this % get_adios_shape_n(shape(x), shape_dims, start_dims, count_dims, adios_varname)
    CALL adios2_define_variable(var, this%io, adios_varname, adios2_type_double_precision, 2, &
      shape_dims, start_dims, count_dims, &
      adios2_constant_dims, ierr)
  end if
  CALL adios2_put(this%engine, var, x, ierr)

END SUBROUTINE writer_real_t_2

SUBROUTINE writer_real_t(this, varname, x)

  IMPLICIT NONE

  CLASS(AdiosWriter_t) :: this
  CHARACTER(*), intent(in) :: varname
  REAL(KIND=dp), intent(in), dimension(:) :: x

  INTEGER(KIND=8), dimension(1) :: shape_dims, start_dims, count_dims
  INTEGER :: ierr
  CHARACTER(ADIOS2_MAX_VARNAME_LEN) :: adios_varname
  TYPE(adios2_variable) :: var

  call this % make_varname(varname, adios_varname)

  call adios2_inquire_variable(var, this%io, adios_varname, ierr)
  if(.not. var % valid) then
    CALL this % get_adios_shape_n(shape(x), shape_dims, start_dims, count_dims, adios_varname)

    CALL adios2_define_variable(var, this%io, adios_varname, adios2_type_double_precision, 1, &
      shape_dims, start_dims, count_dims, adios2_constant_dims, ierr)
  end if
  CALL adios2_put(this%engine, var, x, ierr)

END SUBROUTINE writer_real_t

END MODULE
