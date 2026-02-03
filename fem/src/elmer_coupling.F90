
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

  USE yac

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

  USE yac
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
        comp_id, corner_point_id, timestepstring, cell_point_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: corner_point_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    INTEGER :: nbr_vertices, nbr_cells
    INTEGER :: i

    nbr_vertices = yac_fget_points_size(corner_point_id)
    nbr_cells = yac_fget_points_size(cell_point_id)

    ! register T_ice field in YAC
    CALL yac_fdef_field( &
      t_ice_field_name, comp_id, (/corner_point_id/), 1, t_ice_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, t_ice_field_id);

    ! allocate and initialise T_ice field buffer
    ALLOCATE(t_ice_field(nbr_vertices, t_ice_collection_size))
    t_ice_field = 0.0

    ! register smb field in YAC
    CALL yac_fdef_field( &
      smb_field_name, comp_id, (/cell_point_id/), 1, smb_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, smb_field_id);

    ! allocate and initialise smb field buffer
    ALLOCATE(smb_field(nbr_cells, smb_collection_size))
    smb_field = 0.0

    ! register runoff field in YAC
    CALL yac_fdef_field( &
      runoff_field_name, comp_id, (/corner_point_id/), 1, runoff_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, runoff_field_id);

    ! allocate and initialise runoff field buffer
    ALLOCATE(runoff_field(nbr_vertices, runoff_collection_size))
    runoff_field = 0.0

    ! register surface_height field in YAC
    CALL yac_fdef_field( &
      surface_height_field_name, comp_id, (/corner_point_id/), 1, surface_height_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, surface_height_field_id);

    ! allocate and initialise surface_height field buffer
    ALLOCATE(surface_height_field(nbr_vertices, surface_height_collection_size))
    surface_height_field = 42.0

  END SUBROUTINE construct_elmer_ebfm_coupling

  SUBROUTINE construct_elmer_ebfm_coupling_post_sync( &
    comm_rank, elmer_comp_name, elmer_grid_name)

    INTEGER, INTENT(IN) :: comm_rank
    CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
    CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name

    ! after synchronisation or the end of the definition phase YAC can be
    ! queried about various information

    IF (comm_rank /= 0) RETURN

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
        src_field_timestep = &
          yac_fget_field_timestep(src_comp_name, src_grid_name, src_field_name)
        src_field_metadata = &
          yac_fget_field_metadata(src_comp_name, src_grid_name, src_field_name)

        PRINT *, "ELMER: field ", field_name, ":"
        PRINT *, "ELMER:  - source:"
        PRINT *, "ELMER:    - component: ", src_comp_name
        PRINT *, "ELMER:    - grid:      ", src_grid_name
        PRINT *, "ELMER:    - timestep:  ", src_field_timestep
        PRINT *, "ELMER:    - metadata:  ", src_field_metadata

      END IF

      ! TODO Add something analogously for fields where Elmer is the source?

    END SUBROUTINE print_field_info

  END SUBROUTINE construct_elmer_ebfm_coupling_post_sync

  SUBROUTINE elmer_ebfm_interface(comm_rank)

    INTEGER, INTENT(IN) :: comm_rank

    INTEGER :: info, err

    PRINT *, "IN EBFM_INTERFACE" , comm_rank
    ! checks whether the T_ice field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(t_ice_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (comm_rank == 0) THEN

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
      PRINT *, "CALLING FGET for TICE" , comm_rank
      CALL yac_fget( &
        t_ice_field_id, SIZE(t_ice_field, 1), SIZE(t_ice_field, 2), t_ice_field, &
        info, err)
      PRINT *, "AFTER FGET for TICE" , comm_rank

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

      IF (comm_rank == 0) THEN

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

      IF (comm_rank == 0) THEN

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
    PRINT *, "BEFORE ICE_SHEET_HEIGHT" , comm_rank
    IF (yac_fget_role_from_field_id(surface_height_field_id) == &
        YAC_EXCHANGE_TYPE_SOURCE) THEN

      CALL yac_fget_action(surface_height_field_id, info)

      IF (comm_rank == 0) THEN

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
        PRINT *, "BEFORE FPUT for ICE_SHEET_HEIGHT" , comm_rank
        CALL yac_fput( &
          surface_height_field_id, SIZE(surface_height_field, 1), SIZE(surface_height_field, 2), surface_height_field, &
          info, err)
        PRINT *, "AFTER FPUT for ICE_SHEET_HEIGHT" , comm_rank
      ELSE IF (info == YAC_ACTION_NONE) THEN
        PRINT *, "BEFORE FUPDATE for ICE_SHEET_HEIGHT" , comm_rank
        CALL yac_fupdate(surface_height_field_id)
        PRINT *, "AFTER FUPDATE for ICE_SHEET_HEIGHT" , comm_rank
      END IF
    END IF

  END SUBROUTINE elmer_ebfm_interface

  SUBROUTINE destruct_elmer_ebfm_coupling()

    ! clean up
    DEALLOCATE(t_ice_field, smb_field, runoff_field, surface_height_field)

  END SUBROUTINE destruct_elmer_ebfm_coupling

END MODULE elmer_ebfm_coupling


MODULE elmer_icon_coupling

  USE yac
  USE elmer_coupling_utils, ONLY: yac_action_to_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_elmer_icon_coupling
  PUBLIC :: construct_elmer_icon_coupling_post_sync
  PUBLIC :: destruct_elmer_icon_coupling
  PUBLIC :: elmer_icon_interface

  INTEGER :: clt_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: clt_field_name = "tas"
  INTEGER :: clt_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: clt_field(:,:)

  INTEGER :: pr_field_id = -1
  CHARACTER(LEN=*), PARAMETER :: pr_field_name = "pr_snow"
  INTEGER :: pr_collection_size = 1
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE :: pr_field(:,:)


CONTAINS

  SUBROUTINE construct_elmer_icon_coupling( &
        comp_id, corner_point_id, timestepstring, cell_point_id)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: corner_point_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    INTEGER :: nbr_vertices, nbr_cells
    INTEGER :: i

    nbr_vertices = yac_fget_points_size(corner_point_id)
    nbr_cells = yac_fget_points_size(cell_point_id)

    ! register total cloud cover field in YAC
    CALL yac_fdef_field( &
      clt_field_name, comp_id, (/corner_point_id/), 1, clt_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, clt_field_id);

    ! allocate and initialise total cloud cover field buffer
    ALLOCATE(clt_field(nbr_vertices, clt_collection_size))
    clt_field = 0.0

    ! register precipitation flux field in YAC
    CALL yac_fdef_field( &
      pr_field_name, comp_id, (/cell_point_id/), 1, pr_collection_size, &
      timestepstring, YAC_TIME_UNIT_HOUR, pr_field_id)

    ! allocate and initialise precipitation flux field buffer
    ALLOCATE(pr_field(nbr_cells, pr_collection_size))
    pr_field = 0.0

  END SUBROUTINE construct_elmer_icon_coupling

  SUBROUTINE construct_elmer_icon_coupling_post_sync( &
    comm_rank, elmer_comp_name, elmer_grid_name)

    INTEGER, INTENT(IN) :: comm_rank
    CHARACTER(LEN=*), INTENT(IN) :: elmer_comp_name
    CHARACTER(LEN=*), INTENT(IN) :: elmer_grid_name

    ! after synchronisation or the end of the definition phase YAC can be
    ! queried about various information

    IF (comm_rank /= 0) RETURN

    CALL print_field_info(elmer_comp_name, elmer_grid_name, pr_field_name)
    CALL print_field_info(elmer_comp_name, elmer_grid_name, clt_field_name)

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
        src_field_timestep = &
          yac_fget_field_timestep(src_comp_name, src_grid_name, src_field_name)
        src_field_metadata = &
          yac_fget_field_metadata(src_comp_name, src_grid_name, src_field_name)

        

        
        PRINT *, "field ", field_name, ":"
        PRINT *, " - source:"
        PRINT *, "   - component: ", src_comp_name
        PRINT *, "   - grid:      ", src_grid_name
        PRINT *, "   - timestep:  ", src_field_timestep
        PRINT *, "   - metadata:  ", src_field_metadata

      END IF

    END SUBROUTINE print_field_info

  END SUBROUTINE construct_elmer_icon_coupling_post_sync

  SUBROUTINE elmer_icon_interface(comm_rank)

    INTEGER, INTENT(IN) :: comm_rank

    INTEGER :: info, err

    ! checks whether the total cloud cover field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(clt_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (comm_rank == 0) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the total cloud cover field and print out some information
        CALL yac_fget_action(clt_field_id, info)
        PRINT *, "call get for field: ", TRIM(clt_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(clt_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for total cloud cover field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, total cloud cover field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        clt_field_id, SIZE(clt_field, 1), SIZE(clt_field, 2), clt_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal total cloud cover field

      END IF
    END IF

    ! checks whether the precipitation flux field is defined as a target
    ! in a couple
    IF (yac_fget_role_from_field_id(pr_field_id) == &
        YAC_EXCHANGE_TYPE_TARGET) THEN

      IF (comm_rank == 0) THEN

        ! get the action executed by YAC in the next get operation called for
        ! the precipitation flux field and print out some information
        CALL yac_fget_action(pr_field_id, info)
        PRINT *, "call get for field: ", TRIM(pr_field_name), &
                 " datatime: ", TRIM(yac_fget_field_datetime(pr_field_id)), &
                 " action: ", TRIM(yac_action_to_string(info))
      END IF

      ! execute get operation for precipitation flux field
      ! * if this is a coupling timestep, this will block until the data has
      !   been received
      ! * if this is not a coupling timestep, precipitation flux field buffer
      !   is left untouched and routine will return immediately
      CALL yac_fget( &
        pr_field_id, SIZE(pr_field, 1), SIZE(pr_field, 2), pr_field, &
        info, err)

      ! if this was a coupling timestep
      IF ((info == YAC_ACTION_COUPLING) .OR. &
          (info == YAC_ACTION_GET_FOR_RESTART)) THEN

        ! prepare received data for elmer

        ! update elmer internal precipitation flux field

      END IF
    END IF

  END SUBROUTINE elmer_icon_interface

  SUBROUTINE destruct_elmer_icon_coupling()

    ! clean up
    DEALLOCATE(pr_field, clt_field)

  END SUBROUTINE destruct_elmer_icon_coupling

END MODULE elmer_icon_coupling

MODULE elmer_coupling

  USE mpi
  USE yac

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: coupling_init, coupling_finalize
  PUBLIC :: coupling_setup
  PUBLIC :: coupler_get_code_id
  PUBLIC :: elmer_coupling_mpi_handshake

  INTEGER, PARAMETER, PRIVATE :: MAX_CHARLEN = 132
  INTEGER, PARAMETER, PUBLIC :: elmer_coupling_MAX_GROUPNAME_LEN = MAX_CHARLEN

  ! TODO: Allow to set component name from outside
  ! CHARACTER(LEN=MAX_CHARLEN) :: ELMER_COMP_NAME
  ! to make sure to have a single YAML file in case of multiple Elmer/Ice domains
  CHARACTER(LEN=MAX_CHARLEN), PARAMETER :: ELMER_COMP_NAME = "elmerice"
  CHARACTER(LEN=MAX_CHARLEN), PARAMETER :: ELMER_GRID_NAME = "elmer_grid"

  INTEGER :: comp_id
  INTEGER :: comm_rank, comm_size

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
  ! SUBROUTINE coupling_init(coupling_config_file, elmer_comm, yac_comm, comp_name)
  SUBROUTINE coupling_init(coupling_config_file, elmer_comm, yac_comm)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: coupling_config_file
    INTEGER, INTENT(IN) :: elmer_comm
    INTEGER, INTENT(IN) :: yac_comm

    ! CHARACTER(LEN=elmer_coupling_MAX_GROUPNAME_LEN), INTENT(IN) :: comp_name

    INTEGER :: ierror


    ! initialise YAC
    ! * is collective operation on yac_comm
    ! * should be called as early as possible
    ! * in case not all Elmer processes want to take part in the coupling
    !   this has to be modified
    !   (see:
    !     https://dkrz-sw.gitlab-pages.dkrz.de/yac/d4/d40/init_yac_detail.html)
    ! * will call MPI_Init, if not yet called by the user
    PRINT *, "HELLO FROM COUPLING INIT", ELMER_COMP_NAME
    CALL yac_finit_comm (yac_comm)

    ! read configuration file
    ! * contains calendar, start- and end-date
    ! * defines couplings
    CALL yac_fread_config_yaml(TRIM(coupling_config_file))

    ! define component
    ! * is collective operation for all processes that initialised YAC
    ! TODO: Allow to set component name from outside; currently hard-coded
    ! ELMER_COMP_NAME = comp_name
    CALL yac_fdef_comp(ELMER_COMP_NAME, comp_id)

    ! get number of ranks for the elmer component
    ! (required for reading in the grid data)
    CALL MPI_Comm_rank(elmer_comm, comm_rank, ierror)
    CALL MPI_Comm_size(elmer_comm, comm_size, ierror)

  END SUBROUTINE coupling_init

  SUBROUTINE coupling_setup(grid_dir, num_parts, timestepstring)

    USE :: elmer_ebfm_coupling
    USE :: elmer_icon_coupling
    USE, INTRINSIC :: iso_c_binding, ONLY: C_INT, C_DOUBLE, C_PTR, C_F_POINTER, C_NULL_CHAR

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: grid_dir
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    INTEGER, INTENT(IN) :: num_parts

    INTEGER :: grid_id, corner_point_id, cell_point_id

    INTEGER(KIND=C_INT) :: nbr_vertices
    INTEGER(KIND=C_INT) :: nbr_cells
    TYPE(C_PTR)         :: num_vertices_per_cell_c_ptr ! int **
    TYPE(C_PTR)         :: x_vertices_c_ptr ! double **
    TYPE(C_PTR)         :: y_vertices_c_ptr ! double **
    TYPE(C_PTR)         :: x_cells_c_ptr ! double **
    TYPE(C_PTR)         :: y_cells_c_ptr ! double **
    TYPE(C_PTR)         :: cell_ids_c_ptr ! int **
    TYPE(C_PTR)         :: vertex_ids_c_ptr ! int **
    TYPE(C_PTR)         :: cell_to_vertex_c_ptr ! int **

    INTEGER(KIND=C_INT), POINTER :: num_vertices_per_cell_c(:)
    REAL(KIND=C_DOUBLE), POINTER :: x_vertices_c(:)
    REAL(KIND=C_DOUBLE), POINTER :: y_vertices_c(:)
    REAL(KIND=C_DOUBLE), POINTER :: x_cells_c(:)
    REAL(KIND=C_DOUBLE), POINTER :: y_cells_c(:)
    INTEGER(KIND=C_INT), POINTER :: cell_ids_c(:)
    INTEGER(KIND=C_INT), POINTER :: vertex_ids_c(:)
    INTEGER(KIND=C_INT), POINTER :: cell_to_vertex_c(:)

    INTEGER, ALLOCATABLE          :: num_vertices_per_cell(:)
    DOUBLE PRECISION, ALLOCATABLE :: x_vertices(:)
    DOUBLE PRECISION, ALLOCATABLE :: y_vertices(:)
    DOUBLE PRECISION, ALLOCATABLE :: x_cells(:)
    DOUBLE PRECISION, ALLOCATABLE :: y_cells(:)
    INTEGER, ALLOCATABLE          :: cell_ids(:)
    INTEGER, ALLOCATABLE          :: vertex_ids(:)
    INTEGER, ALLOCATABLE          :: cell_to_vertex(:)

    INTERFACE

      SUBROUTINE read_grid_c(grid_dir, rank, size, num_parts, &
                             nbr_vertices, nbr_cells, num_vertices_per_cell, &
                             x_vertices, y_vertices, x_cells, y_cells, &
                             cell_ids, vertex_ids, cell_to_vertex) &
        bind ( C, name='read_grid' )

        USE, INTRINSIC :: iso_c_binding

        CHARACTER(KIND=C_CHAR) :: grid_dir(*)
        INTEGER(KIND=C_INT), VALUE :: rank
        INTEGER(KIND=C_INT), VALUE :: size
        INTEGER(KIND=C_INT), VALUE :: num_parts
        INTEGER(KIND=C_INT)        :: nbr_vertices
        INTEGER(KIND=C_INT)        :: nbr_cells
        TYPE(C_PTR)                :: num_vertices_per_cell ! int **
        TYPE(C_PTR)                :: x_vertices ! double **
        TYPE(C_PTR)                :: y_vertices ! double **
        TYPE(C_PTR)                :: x_cells ! double **
        TYPE(C_PTR)                :: y_cells ! double **
        TYPE(C_PTR)                :: cell_ids ! int **
        TYPE(C_PTR)                :: vertex_ids ! int **
        TYPE(C_PTR)                :: cell_to_vertex ! int **

      END SUBROUTINE read_grid_c

      SUBROUTINE free_c ( ptr ) bind ( c, NAME='free' )

       USE, INTRINSIC :: iso_c_binding, ONLY: C_PTR

       TYPE(C_PTR), VALUE :: ptr

      END SUBROUTINE free_c

    END INTERFACE

    PRINT *, "READ GRID FROM FILE"
    ! get grid data from elmer component
    ! in the case of the dummy, we have to read it from file

    ! read grid data from file
    ! * each process only reads in its local part of the grid
    ! * in Elmer this information probably already available and does not have
    !   to be read from file
    PRINT *, "COMM_RANK", comm_rank
    PRINT *, "COMM_SIZE", comm_size
    PRINT *, "NUMPARTS", num_parts
    CALL read_grid_c( &
      TRIM(grid_dir) // c_null_char, comm_rank, comm_size, num_parts, &
      nbr_vertices, nbr_cells, num_vertices_per_cell_c_ptr, &
      x_vertices_c_ptr, y_vertices_c_ptr, x_cells_c_ptr, y_cells_c_ptr, &
      cell_ids_c_ptr, vertex_ids_c_ptr, cell_to_vertex_c_ptr)

    PRINT *, "After READINF GRID FROM FILE"
    CALL C_F_POINTER( &
      num_vertices_per_cell_c_ptr, num_vertices_per_cell_c, shape=[nbr_cells])
    CALL C_F_POINTER(x_vertices_c_ptr, x_vertices_c, shape=[nbr_vertices])
    CALL C_F_POINTER(y_vertices_c_ptr, y_vertices_c, shape=[nbr_vertices])
    CALL C_F_POINTER(x_cells_c_ptr, x_cells_c, shape=[nbr_cells])
    CALL C_F_POINTER(y_cells_c_ptr, y_cells_c, shape=[nbr_cells])
    CALL C_F_POINTER(cell_ids_c_ptr, cell_ids_c, shape=[nbr_cells])
    CALL C_F_POINTER(vertex_ids_c_ptr, vertex_ids_c, shape=[nbr_vertices])
    CALL C_F_POINTER( &
      cell_to_vertex_c_ptr, cell_to_vertex_c, shape=[SUM(num_vertices_per_cell_c)])

    num_vertices_per_cell = num_vertices_per_cell_c
    x_vertices = x_vertices_c
    y_vertices = y_vertices_c
    x_cells = x_cells_c
    y_cells = y_cells_c
    cell_ids = cell_ids_c
    vertex_ids = vertex_ids_c
    cell_to_vertex = cell_to_vertex_c + 1

    CALL free_c(num_vertices_per_cell_c_ptr)
    CALL free_c(x_vertices_c_ptr)
    CALL free_c(y_vertices_c_ptr)
    CALL free_c(x_cells_c_ptr)
    CALL free_c(y_cells_c_ptr)
    CALL free_c(cell_ids_c_ptr)
    CALL free_c(vertex_ids_c_ptr)
    CALL free_c(cell_to_vertex_c_ptr)


    ! register Elmer grid in YAC
    ! * is defined as an unstructured grid
    PRINT *, "BEFORE GRID DEF"
    CALL yac_fdef_grid( &
      ELMER_GRID_NAME, nbr_vertices, nbr_cells, SUM(num_vertices_per_cell), &
      num_vertices_per_cell, x_vertices, y_vertices, cell_to_vertex, grid_id)

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
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices, &
      corner_point_id)
    CALL yac_fdef_points( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, cell_point_id)

    PRINT *, "PRECIP TIMESTEP in HOURS", timestepstring
    ! construct coupling between Elmer/Ice and ICON
    !CALL construct_elmer_icon_coupling(comp_id, corner_point_id, timestepstring, cell_point_id)
    CALL construct_elmer_ebfm_coupling(comp_id, corner_point_id, timestepstring, cell_point_id)

    PRINT *, "AFTER constructing_elmer_ebfm_coupling", timestepstring
    ! sychronizes all definitions between all components
    ! * afterwards the exchange information can be queried
    ! * this is optional
    CALL yac_fsync_def()
    PRINT *, "AFTER synchronisation", timestepstring

    ! construct coupling between Elmer/Ice and ICON (using sychronized
    ! information from all components)
    !CALL construct_elmer_icon_coupling_post_sync( &
      !comm_rank, ELMER_COMP_NAME, ELMER_GRID_NAME)
    CALL construct_elmer_ebfm_coupling_post_sync( &
      comm_rank, ELMER_COMP_NAME, ELMER_GRID_NAME)

    ! end of definition phase
    ! * collective operation for all processes that initialised YAC
    ! * exchanges information between all processes (various query routines are
    !   available to access this data; also from other components)
    ! * computes weights required for all defined couples
    CALL yac_fenddef()
    PRINT *, "FINISHED COUPLING SETUP", timestepstring

  END SUBROUTINE coupling_setup

  SUBROUTINE coupling_finalize()

    USE elmer_ebfm_coupling
    USE elmer_icon_coupling

    IMPLICIT NONE

    PRINT *, "DESTRCUTING ELMER_ICON_COUPLING"
    !CALL destruct_elmer_icon_coupling()
    PRINT *, "DESTRCUTING ELMER_EBFM_COUPLING"
    CALL destruct_elmer_ebfm_coupling()

    ! finalise YAC
    ! * if user has called MPI_Init, he also has to call MPI_Finalize afterwards
    CALL yac_ffinalize()

  END SUBROUTINE coupling_finalize

END MODULE elmer_coupling







