
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
! *  Authors: Moritz Hanke (DKRZ), adapted by Thomas Zwinger, Clemens Schannwell
!    original code provided by Moritz Hanke (DKRZ) in the frame of the TerraDT project
! *  Email:   thomas.zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!> definitions for Elmer YAC coupling.
!------------------------------------------------------------------------------
#include "yac_config.h"

#define YAC_VERSION_GREATER_EQUAL(major, minor, patch) ( \
    (YAC_VERSION_MAJOR > major) || \
    (YAC_VERSION_MAJOR == major && YAC_VERSION_MINOR > minor) || \
    (YAC_VERSION_MAJOR == major && YAC_VERSION_MINOR == minor && YAC_VERSION_PATCH >= patch) \
)

!> Helper module for YAC utility functions used by Elmer coupling modules
!> Outside of MODULE elmer_coupling to avoid circular dependencies
MODULE elmer_coupling_utils

  USE yac, ONLY: YAC_ACTION_COUPLING, YAC_ACTION_PUT_FOR_RESTART, &
    YAC_ACTION_GET_FOR_RESTART, YAC_ACTION_REDUCTION, YAC_ACTION_NONE, &
    YAC_ACTION_OUT_OF_BOUND

  IMPLICIT NONE

  PUBLIC :: yac_action_to_string

CONTAINS

  !> Convert YAC action code to human-readable string
  !> @param action_code YAC action code (e.g., YAC_ACTION_COUPLING)
  !> @return String representation of the action code
  FUNCTION yac_action_to_string(action_code) RESULT(action_str)
    INTEGER, INTENT(IN) :: action_code
    CHARACTER(LEN=30) :: action_str

    SELECT CASE (action_code)
      CASE (YAC_ACTION_COUPLING)
        action_str = "YAC_ACTION_COUPLING"
      CASE (YAC_ACTION_PUT_FOR_RESTART)
        action_str = "YAC_ACTION_PUT_FOR_RESTART"
      CASE (YAC_ACTION_GET_FOR_RESTART)
        action_str = "YAC_ACTION_GET_FOR_RESTART"
      CASE (YAC_ACTION_REDUCTION)
        action_str = "YAC_ACTION_REDUCTION"
      CASE (YAC_ACTION_NONE)
        action_str = "YAC_ACTION_NONE"
      CASE (YAC_ACTION_OUT_OF_BOUND)
        action_str = "YAC_ACTION_OUT_OF_BOUND"
      CASE DEFAULT
        action_str = "UNKNOWN_ACTION"
    END SELECT
  END FUNCTION yac_action_to_string

END MODULE elmer_coupling_utils

MODULE elmer_ebfm_coupling

  USE yac, ONLY: yac_fdef_field, yac_fget_role_from_field_id, &
    yac_fget_field_datetime, yac_fget_field_role, yac_fget_field_timestep, &
    yac_ffield_has_metadata, yac_fget_field_metadata, yac_fget_points_size, &
    yac_fget_field_source, yac_fget, yac_fput, yac_fupdate, yac_fget_action, &
    YAC_TIME_UNIT_ISO_FORMAT, &
    YAC_ACTION_COUPLING, YAC_ACTION_GET_FOR_RESTART, &
    YAC_ACTION_PUT_FOR_RESTART, YAC_ACTION_REDUCTION, YAC_ACTION_NONE, &
    YAC_EXCHANGE_TYPE_SOURCE, YAC_EXCHANGE_TYPE_TARGET

  USE elmer_coupling_utils, ONLY: yac_action_to_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_elmer_ebfm_coupling
  PUBLIC :: construct_elmer_ebfm_coupling_post_sync
  PUBLIC :: destruct_elmer_ebfm_coupling
  PUBLIC :: elmer_ebfm_interface

  ! Fields Elmer receives from EBFM

  INTEGER :: t_ice_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: t_ice_field_name = "T_ice"
  INTEGER :: t_ice_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: t_ice_field(:,:)

  INTEGER :: smb_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: smb_field_name = "smb"
  INTEGER :: smb_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: smb_field(:,:)

  INTEGER :: runoff_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: runoff_field_name = "runoff"
  INTEGER :: runoff_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: runoff_field(:,:)

  ! Fields Elmer sends to EBFM

  INTEGER :: surface_height_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: surface_height_field_name = "h"
  INTEGER :: surface_height_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: surface_height_field(:,:)

CONTAINS

  SUBROUTINE construct_elmer_ebfm_coupling( &
        comp_id, corner_point_id, iso8601_timestep, cell_point_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: corner_point_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: iso8601_timestep

    INTEGER :: nbr_vertices, nbr_cells

    nbr_vertices = yac_fget_points_size(corner_point_id)
    nbr_cells = yac_fget_points_size(cell_point_id)

    ! register T_ice field in YAC
    CALL yac_fdef_field( &
      t_ice_field_name, comp_id, (/corner_point_id/), 1, t_ice_collection_size, &
      iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, t_ice_field_id);

    ! allocate and initialise T_ice field buffer
    ALLOCATE(t_ice_field(nbr_vertices, t_ice_collection_size))
    t_ice_field = 0.0

    ! register smb field in YAC
    CALL yac_fdef_field( &
      smb_field_name, comp_id, (/cell_point_id/), 1, smb_collection_size, &
      iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, smb_field_id);

    ! allocate and initialise smb field buffer
    ALLOCATE(smb_field(nbr_cells, smb_collection_size))
    smb_field = 0.0

    ! register runoff field in YAC
    CALL yac_fdef_field( &
      runoff_field_name, comp_id, (/cell_point_id/), 1, runoff_collection_size, &
      iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, runoff_field_id);

    ! allocate and initialise runoff field buffer
    ALLOCATE(runoff_field(nbr_cells, runoff_collection_size))
    runoff_field = 0.0

    ! register surface_height field in YAC
    CALL yac_fdef_field( &
      surface_height_field_name, comp_id, (/corner_point_id/), 1, surface_height_collection_size, &
      iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, surface_height_field_id);

    ! allocate and initialise surface_height field buffer
    ALLOCATE(surface_height_field(nbr_vertices, surface_height_collection_size))
    surface_height_field = 42.0

  END SUBROUTINE construct_elmer_ebfm_coupling

  SUBROUTINE construct_elmer_ebfm_coupling_post_sync( &
    is_root_rank, elmer_comp_name, elmer_grid_name)

    LOGICAL, INTENT(IN) :: is_root_rank
    CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
    CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name

    ! after synchronisation or the end of the definition phase YAC can be
    ! queried about various information

    IF (.NOT. is_root_rank) RETURN

    CALL print_field_info(elmer_comp_name, elmer_grid_name, t_ice_field_name)
    CALL print_field_info(elmer_comp_name, elmer_grid_name, smb_field_name)
    CALL print_field_info(elmer_comp_name, elmer_grid_name, runoff_field_name)
    CALL print_field_info(elmer_comp_name, elmer_grid_name, surface_height_field_name)

  CONTAINS

    SUBROUTINE print_field_info(elmer_comp_name, elmer_grid_name, field_name)

      CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
      CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name
      CHARACTER(LEN=*), INTENT(IN) :: field_name

      CHARACTER(LEN=:), ALLOCATABLE :: src_comp_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_grid_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_timestep
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_metadata

      IF (yac_fget_field_role( &
            elmer_comp_name, elmer_grid_name, field_name) == &
            YAC_EXCHANGE_TYPE_TARGET) THEN

#if (YAC_VERSION_MAJOR >= 3) && (YAC_VERSION_MINOR >= 6)
        CALL yac_fget_field_source( &
          elmer_comp_name, elmer_grid_name, field_name, &
          src_comp_name, src_grid_name, src_field_name);
#else
        src_comp_name = "ebfm"
        src_grid_name = "ebfm_grid"
        src_field_name = field_name
#endif

        IF (yac_fget_field_role( &
              src_comp_name, src_grid_name, src_field_name) == &
              YAC_EXCHANGE_TYPE_SOURCE) THEN

          src_field_timestep = &
            yac_fget_field_timestep( &
              src_comp_name, src_grid_name, src_field_name)

          IF (yac_ffield_has_metadata( &
                src_comp_name, src_grid_name, src_field_name)) THEN
            src_field_metadata = &
              yac_fget_field_metadata( &
                src_comp_name, src_grid_name, src_field_name)
          ELSE
            src_field_metadata = "N/A"
          END IF

          PRINT *, "ELMER: field ", field_name, ":"
          PRINT *, "ELMER:  - source:"
          PRINT *, "ELMER:    - component: ", src_comp_name
          PRINT *, "ELMER:    - grid:      ", src_grid_name
          PRINT *, "ELMER:    - timestep:  ", src_field_timestep
          PRINT *, "ELMER:    - metadata:  ", src_field_metadata

        END IF

      END IF

      ! TODO Add something analogously for fields where Elmer is the source?

    END SUBROUTINE print_field_info

  END SUBROUTINE construct_elmer_ebfm_coupling_post_sync

  SUBROUTINE elmer_ebfm_interface(is_root_rank)

    LOGICAL, INTENT(IN) :: is_root_rank

    INTEGER :: info, err

    ! checks whether the T_ice field is defined as a target in a couple
    IF (yac_fget_role_from_field_id(t_ice_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the T_ice field and print out some information
        CALL yac_fget_action(t_ice_field_id, info)
        PRINT *, "ELMER: call get for field: ", TRIM(t_ice_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(t_ice_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for T_ice field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, T_ice field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        t_ice_field_id, SIZE(t_ice_field, 1), SIZE(t_ice_field, 2), t_ice_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal T_ice field

      END IF
    END IF

    ! checks whether the smb field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(smb_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the smb field and print out some information
        CALL yac_fget_action(smb_field_id, info)
        PRINT *, "ELMER: call get for field: ", TRIM(smb_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(smb_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for smb field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, smb field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        smb_field_id, SIZE(smb_field, 1), SIZE(smb_field, 2), smb_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal smb field

      END IF
    END IF

    ! checks whether the runoff field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(runoff_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the runoff field and print out some information
        CALL yac_fget_action(runoff_field_id, info)
        PRINT *, "ELMER: call get for field: ", TRIM(runoff_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(runoff_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for runoff field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, runoff field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        runoff_field_id, SIZE(runoff_field, 1), SIZE(runoff_field, 2), runoff_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal runoff field

      END IF
    END IF

    ! checks whether the surface height field is defined as a source
    ! in a couple
    IF (yac_fget_role_from_field_id(surface_height_field_id) == &
        YAC_EXCHANGE_TYPE_SOURCE) THEN

      CALL yac_fget_action(surface_height_field_id, info)

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next put operation called for
        ! the surface_height field and print out some information
        PRINT *, "ELMER: call put for field: ", TRIM(surface_height_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(surface_height_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_PUT_FOR_RESTART) .OR. &
          (info == YAC_ACTION_REDUCTION)) THEN

        ! get data to be sent from elmer

        ! execute put operation for surface_height field
        ! * if this is a coupling timestep, this will block until the data has
        !   been received
        ! * if this is not a coupling timestep, surface_height field buffer
        !   is left untouched and routine will return immediately
        CALL yac_fput( &
          surface_height_field_id, SIZE(surface_height_field, 1), SIZE(surface_height_field, 2), surface_height_field, &
          info, err)
      ELSE IF (info == YAC_ACTION_NONE) THEN
        CALL yac_fupdate(surface_height_field_id)
      END IF
    END IF

  END SUBROUTINE elmer_ebfm_interface

  SUBROUTINE destruct_elmer_ebfm_coupling()

    ! clean up
    DEALLOCATE(t_ice_field, smb_field, runoff_field, surface_height_field)

  END SUBROUTINE destruct_elmer_ebfm_coupling

END MODULE elmer_ebfm_coupling


MODULE elmer_icon_coupling

  USE yac, ONLY: yac_fdef_field, yac_fdef_field_mask, yac_fget_role_from_field_id, &
    yac_fget_field_datetime, yac_fget_field_role, yac_fget_field_timestep, &
    yac_ffield_has_metadata, yac_fget_field_metadata, yac_fget_points_size, &
    yac_fget_field_source, yac_fget, yac_fput, yac_fupdate, yac_fget_action, &
    YAC_TIME_UNIT_ISO_FORMAT, &
    YAC_ACTION_COUPLING, YAC_ACTION_GET_FOR_RESTART, &
    YAC_ACTION_PUT_FOR_RESTART, YAC_ACTION_REDUCTION, YAC_ACTION_NONE, &
    YAC_EXCHANGE_TYPE_SOURCE, YAC_EXCHANGE_TYPE_TARGET

  USE elmer_coupling_utils, ONLY: yac_action_to_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_elmer_icon_coupling
  PUBLIC :: construct_elmer_icon_coupling_post_sync
  PUBLIC :: destruct_elmer_icon_coupling
  PUBLIC :: elmer_icon_interface

  INTEGER :: t_oce_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: t_oce_field_name = "temp_oce"
  INTEGER :: t_oce_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: t_oce_field(:,:)

  INTEGER :: sal_oce_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: sal_oce_field_name = "sal_oce"
  INTEGER :: sal_oce_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: sal_oce_field(:,:)

CONTAINS

  SUBROUTINE construct_elmer_icon_coupling( &
        comp_id, corner_point_id, iso8601_timestep, cell_point_id, boundary_corner_mask_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: corner_point_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: iso8601_timestep
    INTEGER, INTENT(IN) :: boundary_corner_mask_id

    INTEGER :: nbr_vertices

    nbr_vertices = yac_fget_points_size(corner_point_id)

    ! register ocean temperature field in YAC
    CALL yac_fdef_field_mask( &
      t_oce_field_name, comp_id, (/corner_point_id/), (/boundary_corner_mask_id/), 1, &
      t_oce_collection_size, iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, t_oce_field_id)

    ! allocate and initialise ocean temperature field buffer
    ALLOCATE(t_oce_field(nbr_vertices, t_oce_collection_size))

    ! register ocean salinity field in YAC
    CALL yac_fdef_field_mask( &
      sal_oce_field_name, comp_id, (/corner_point_id/), (/boundary_corner_mask_id/), 1, &
      sal_oce_collection_size, iso8601_timestep, YAC_TIME_UNIT_ISO_FORMAT, sal_oce_field_id)

    ! allocate and initialise ocean salinity field buffer
    ALLOCATE(sal_oce_field(nbr_vertices, sal_oce_collection_size))

  END SUBROUTINE construct_elmer_icon_coupling

  SUBROUTINE construct_elmer_icon_coupling_post_sync( &
    is_root_rank, elmer_comp_name, elmer_grid_name)

    LOGICAL, INTENT(IN) :: is_root_rank
    CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
    CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name

    ! after synchronisation or the end of the definition phase YAC can be
    ! queried about various information

    IF (.NOT. is_root_rank) RETURN

    CALL print_field_info(elmer_comp_name, elmer_grid_name, t_oce_field_name)
    CALL print_field_info(elmer_comp_name, elmer_grid_name, sal_oce_field_name)

  CONTAINS

    SUBROUTINE print_field_info(elmer_comp_name, elmer_grid_name, field_name)

      CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
      CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name
      CHARACTER(LEN=*), INTENT(IN) :: field_name

      CHARACTER(LEN=:), ALLOCATABLE :: src_comp_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_grid_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_name
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_timestep
      CHARACTER(LEN=:), ALLOCATABLE :: src_field_metadata

      IF (yac_fget_field_role( &
            elmer_comp_name, elmer_grid_name, field_name) == &
            YAC_EXCHANGE_TYPE_TARGET) THEN

#if YAC_VERSION_GREATER_EQUAL(3, 6, 0)
        CALL yac_fget_field_source( &
          elmer_comp_name, elmer_grid_name, field_name, &
          src_comp_name, src_grid_name, src_field_name);
#else
        src_comp_name = "icon"
        src_grid_name = "icon_grid"
        src_field_name = field_name
#endif

        IF (yac_fget_field_role( &
              elmer_comp_name, elmer_grid_name, field_name) == &
              YAC_EXCHANGE_TYPE_TARGET) THEN

          src_field_timestep = &
            yac_fget_field_timestep( &
              elmer_comp_name, elmer_grid_name, src_field_name)

          IF (yac_ffield_has_metadata( &
                src_comp_name, src_grid_name, src_field_name)) THEN
            src_field_metadata = &
              yac_fget_field_metadata( &
                src_comp_name, src_grid_name, src_field_name)
          ELSE
            src_field_metadata = "N/A"
          END IF

          PRINT *, "ELMER: field ", field_name, ":"
          PRINT *, "ELMER:  - source:"
          PRINT *, "ELMER:    - component: ", src_comp_name
          PRINT *, "ELMER:    - grid:      ", src_grid_name
          PRINT *, "ELMER:    - timestep:  ", src_field_timestep
          PRINT *, "ELMER:    - metadata:  ", src_field_metadata

        END IF

      END IF

    END SUBROUTINE print_field_info

  END SUBROUTINE construct_elmer_icon_coupling_post_sync

  SUBROUTINE elmer_icon_interface(is_root_rank)

    LOGICAL, INTENT(IN) :: is_root_rank

    INTEGER :: info, err

    ! checks whether the ocean temperature field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(t_oce_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the total cloud cover field and print out some information
        CALL yac_fget_action(t_oce_field_id, info)
        PRINT *, "call get for field: ", TRIM(t_oce_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(t_oce_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for ocean temperature field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, ocean temperature field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        t_oce_field_id, SIZE(t_oce_field, 1), SIZE(t_oce_field, 2), t_oce_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal ocean temperature field

      END IF
    END IF

    ! checks whether the ocean salinity field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(sal_oce_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (is_root_rank) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the precipitation flux field and print out some information
        CALL yac_fget_action(sal_oce_field_id, info)
        PRINT *, "call get for field: ", TRIM(sal_oce_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(sal_oce_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for ocean salinity field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, ocean salinity field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        sal_oce_field_id, SIZE(sal_oce_field, 1), SIZE(sal_oce_field, 2), sal_oce_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal ocean salinity field

      END IF
    END IF

  END SUBROUTINE elmer_icon_interface

  SUBROUTINE destruct_elmer_icon_coupling()

    ! clean up
    DEALLOCATE(t_oce_field, sal_oce_field)

  END SUBROUTINE destruct_elmer_icon_coupling

END MODULE elmer_icon_coupling

MODULE elmer_coupling

  USE mpi, ONLY: MPI_Comm_rank, MPI_Comm_size
  USE yac, ONLY: yac_fmpi_handshake, yac_fget_mpi_handshake_group_name, &
    yac_finit_comm, yac_fread_config_yaml, yac_fdef_comp, yac_fdef_grid, &
    yac_fset_global_index, yac_fdef_points, yac_fdef_mask_named, &
    yac_fsync_def, yac_fenddef, &
    yac_ffinalize, YAC_LOCATION_CELL, YAC_LOCATION_CORNER, &
    yac_fdef_calendar, YAC_PROLEPTIC_GREGORIAN, YAC_YEAR_OF_360_DAYS, &
    YAC_YEAR_OF_365_DAYS
  USE elmer_ebfm_coupling, ONLY: construct_elmer_ebfm_coupling, &
    construct_elmer_ebfm_coupling_post_sync, destruct_elmer_ebfm_coupling
  USE elmer_icon_coupling, ONLY: construct_elmer_icon_coupling, &
    construct_elmer_icon_coupling_post_sync, destruct_elmer_icon_coupling

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coupling_init, coupling_finalize
  PUBLIC :: coupling_setup
  PUBLIC :: coupler_get_code_id
  PUBLIC :: elmer_coupling_mpi_handshake

  INTEGER, PARAMETER, PRIVATE :: MAX_CHARLEN = 132
  INTEGER, PARAMETER, PUBLIC :: elmer_coupling_MAX_GROUPNAME_LEN = MAX_CHARLEN

  ! ROOT_RANK defines rank of this component taking care of logging.
  INTEGER, PARAMETER, PRIVATE :: ROOT_RANK = 0

  ! TODO: Allow to set component name from outside
  ! CHARACTER(LEN=MAX_CHARLEN) :: ELMER_COMP_NAME
  ! to make sure to have a single YAML file in case of multiple Elmer/Ice domains
  CHARACTER(LEN=MAX_CHARLEN), PARAMETER :: ELMER_COMP_NAME = "elmerice"
  CHARACTER(LEN=MAX_CHARLEN), PARAMETER :: ELMER_GRID_NAME = "elmer_grid"

  INTEGER :: comp_id

  ! True if this is the ROOT_RANK of this component.
  LOGICAL, PUBLIC :: is_root_rank

  LOGICAL :: couple_to_ebfm = .FALSE.
  LOGICAL :: couple_to_icon = .FALSE.

CONTAINS

  ! Wrapper for yac_fmpi_handshake
  SUBROUTINE elmer_coupling_mpi_handshake(comm, group_names, group_comms)
    INTEGER, INTENT(IN) :: comm
    CHARACTER(len=*), INTENT(IN) :: group_names(:)
    INTEGER, INTENT(OUT) :: group_comms(:)
    CALL yac_fmpi_handshake(comm, group_names, group_comms)
  END SUBROUTINE elmer_coupling_mpi_handshake

  SUBROUTINE coupler_get_code_id(label)
    CHARACTER(len=elmer_coupling_MAX_GROUPNAME_LEN), INTENT(OUT) :: label

#if YAC_VERSION_GREATER_EQUAL(3, 9, 0)
    ! YAC >= 3.9 provides a function to get the label from the YAC API
    label = yac_fget_mpi_handshake_group_name()
#else
    ! YAC < 3.9 does not provide this function, so we have to set it manually
    label = "yac"
#endif
  END SUBROUTINE coupler_get_code_id

  ! TODO: Refactor to also accept comp_name here
  ! SUBROUTINE coupling_init(elmer_rank, yac_comm, comp_name)
  SUBROUTINE coupling_init(elmer_rank, yac_comm)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: elmer_rank
    INTEGER, INTENT(IN) :: yac_comm

    ! CHARACTER(LEN=elmer_coupling_MAX_GROUPNAME_LEN), INTENT(IN) :: comp_name

    INTEGER :: ierror

    is_root_rank = (elmer_rank == ROOT_RANK)

    ! initialise YAC
    ! * is collective operation on yac_comm
    ! * should be called as early as possible
    ! * in case not all Elmer processes want to take part in the coupling
    !   this has to be modified
    !   (see:
    !     https://dkrz-sw.gitlab-pages.dkrz.de/yac/d4/d40/init_yac_detail.html)
    ! * will call MPI_Init, if not yet called by the user
    CALL yac_finit_comm (yac_comm)

    ! define component
    ! * is collective operation for all processes that initialised YAC
    ! TODO: Allow to set component name from outside; currently hard-coded
    ! ELMER_COMP_NAME = comp_name
    CALL yac_fdef_comp(ELMER_COMP_NAME, comp_id)

  END SUBROUTINE coupling_init

  !> Set up the coupling configuration for YAC
  !>
  !> @param lon_vertices Longitude coordinates of vertices (radians)
  !> @param lat_vertices Latitude coordinates of vertices (radians)
  !> @param lon_cells Longitude coordinates of cell centers (radians)
  !> @param lat_cells Latitude coordinates of cell centers (radians)
  !> @param cell_to_vertex Connectivity array (which vertices belong to which cell)
  !> @param num_vertices_per_cell Number of vertices for each cell
  !> @param cell_ids Global cell IDs
  !> @param vertex_ids Global vertex IDs
  !> @param grid_crs Coordinate reference system (CRS) of the grid (e.g., "EPSG:3413")
  !> @param coupling_config_file YAC coupling YAML configuration file
  !> @param calendar Calendar type (e.g., "proleptic_gregorian")
  !> @param iso8601_start_time Start time in ISO 8601 format (e.g., "2024-01-01T00:00:00Z")
  !> @param iso8601_end_time End time in ISO 8601 format (e.g., "2024-12-31T23:59:59Z")
  !> @param iso8601_timestep Timestep configuration string for YAC (e.g., "PT1H" for 1 hour)
  !> @param couple_to_ebfm_in Enable coupling to EBFM
  !> @param couple_to_icon_in Enable coupling to ICON
  !> @param boundary_corner_mask Logical mask indicating boundary corners
  SUBROUTINE coupling_setup(lon_vertices, lat_vertices, lon_cells, lat_cells, &
                            cell_to_vertex, num_vertices_per_cell, &
                            cell_ids, vertex_ids, &
                            grid_crs, &
                            coupling_config_file, &
                            calendar, &
                            iso8601_start_time, &
                            iso8601_end_time, &
                            iso8601_timestep, &
                            couple_to_ebfm_in, couple_to_icon_in, &
                            boundary_corner_mask)

    USE, INTRINSIC :: iso_c_binding, ONLY: C_INT, C_DOUBLE, C_CHAR

    IMPLICIT NONE

    ! Input arrays (already converted to lon/lat in radians)
    REAL(KIND=C_DOUBLE), INTENT(IN) :: lon_vertices(:)
    REAL(KIND=C_DOUBLE), INTENT(IN) :: lat_vertices(:)
    REAL(KIND=C_DOUBLE), INTENT(IN) :: lon_cells(:)
    REAL(KIND=C_DOUBLE), INTENT(IN) :: lat_cells(:)
    INTEGER, INTENT(IN) :: cell_to_vertex(:)
    INTEGER, INTENT(IN) :: num_vertices_per_cell(:)
    INTEGER, INTENT(IN) :: cell_ids(:)
    INTEGER, INTENT(IN) :: vertex_ids(:)
    CHARACTER(LEN=*), INTENT(IN) :: grid_crs
    CHARACTER(LEN=*), INTENT(IN) :: coupling_config_file

    ! Parameters for definition of time frame (must be consistent with other components)
    CHARACTER(LEN=*), INTENT(IN) :: calendar
    CHARACTER(LEN=*), INTENT(IN) :: iso8601_start_time
    CHARACTER(LEN=*), INTENT(IN) :: iso8601_end_time
    CHARACTER(LEN=*), INTENT(IN) :: iso8601_timestep

    LOGICAL, INTENT(IN) :: couple_to_ebfm_in, couple_to_icon_in
    LOGICAL, INTENT(IN) :: boundary_corner_mask(:)

    ! Local variables
    INTEGER :: grid_id, corner_point_id, cell_point_id, boundary_corner_mask_id
    INTEGER(KIND=C_INT) :: nbr_vertices, nbr_cells

    ! Store coupling flags in module variables for later use
    couple_to_ebfm = couple_to_ebfm_in
    couple_to_icon = couple_to_icon_in

    ! Infer sizes from input arrays
    nbr_vertices = SIZE(lon_vertices)
    nbr_cells = SIZE(lon_cells)

    ! register Elmer grid in YAC
    ! * is defined as an unstructured grid
    CALL yac_fdef_grid( &
      ELMER_GRID_NAME, nbr_vertices, nbr_cells, SUM(num_vertices_per_cell), &
      num_vertices_per_cell, lon_vertices, lat_vertices, cell_to_vertex, grid_id)

    ! define global ids for all cells and vertices
    ! * this information is used to identify cells/vertice between processes
    ! * if this information is not provided, YAC will generate it (resulting
    !   IDs then depend on the number of processes)
    CALL yac_fset_global_index(cell_ids, YAC_LOCATION_CELL, grid_id)
    CALL yac_fset_global_index(vertex_ids, YAC_LOCATION_CORNER, grid_id)

    ! define location at which field data will be provided
    ! * an arbitrary number of locations can be defined (e.g. at vertices,
    !   cell centers, edge middle points)
    CALL yac_fdef_points( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, lon_vertices, lat_vertices, &
      corner_point_id)
    CALL yac_fdef_points( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, lon_cells, lat_cells, cell_point_id)

    ! register boundary corner mask in YAC
    CALL yac_fdef_mask_named( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, boundary_corner_mask, &
      "boundary_corner_mask", boundary_corner_mask_id)

    ! construct coupling between Elmer/Ice and ICON
    IF (couple_to_icon) THEN
        CALL construct_elmer_icon_coupling( &
          comp_id, corner_point_id, iso8601_timestep, cell_point_id, &
          boundary_corner_mask_id)
    END IF
    IF (couple_to_ebfm) THEN
        CALL construct_elmer_ebfm_coupling( &
          comp_id, corner_point_id, iso8601_timestep, cell_point_id)
    END IF
    ! sychronizes all definitions between all components
    ! * afterwards the exchange information can be queried
    ! * this is optional
    CALL yac_fsync_def()

    ! construct coupling between Elmer/Ice and ICON (using sychronized
    ! information from all components)
    IF (couple_to_icon) THEN
        CALL construct_elmer_icon_coupling_post_sync( &
          is_root_rank, ELMER_COMP_NAME, ELMER_GRID_NAME)
    END IF
    IF (couple_to_ebfm) THEN
        CALL construct_elmer_ebfm_coupling_post_sync( &
          is_root_rank, ELMER_COMP_NAME, ELMER_GRID_NAME)
    END IF

    ! end of definition phase
    ! * collective operation for all processes that initialised YAC
    ! * exchanges information between all processes (various query routines are
    !   available to access this data; also from other components)
    ! * computes weights required for all defined couples
    CALL yac_fenddef()

  END SUBROUTINE coupling_setup

  SUBROUTINE coupling_finalize()

    IMPLICIT NONE

    IF (couple_to_ebfm) THEN
      CALL destruct_elmer_ebfm_coupling()
    END IF

    IF (couple_to_icon) THEN
      CALL destruct_elmer_icon_coupling()
    END IF

    ! finalise YAC
    ! * if user has called MPI_Init, he also has to call MPI_Finalize afterwards
    CALL yac_ffinalize()

  END SUBROUTINE coupling_finalize

END MODULE elmer_coupling







