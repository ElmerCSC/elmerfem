!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! *****************************************************************************
! * COUPLING WITH PDAF (Parallel Data Assimilation Framework)
! *  https://pdaf.awi.de/trac/wiki/WikiStart
! *  Code adapted from PDAF templates and examples
! *****************************************************************************
!*  Module for file operations for synthetic observations
MODULE mod_synobs
   USE DefUtils
   USE PDAFUtils, ONLY : mype_filter
   USE netcdf
   USE mod_io

CONTAINS
! *****************************************************************************
! Initialize observation file
!   File dimensions:
!    - dim_obs : process local obs dim (must be constant)
!    - steps : number of time steps (unlimited)
!   Variables:
!    - timestep(steps) : timesteps
!    - obs(dim_obs,steps) : observations
! *****************************************************************************
  SUBROUTINE init_file_syn_obs(dim_obs_p, file_obs)
    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_obs_p       !< process local obs dim
    character(len=*), intent(in) :: file_obs !< Name of observation file


! *** Local variables ***
    CHARACTER(*), PARAMETER :: Caller = "init_file_syn_obs"
    integer :: dimid_iter, dimid_obs              ! IDs for dimensions
    integer :: id_time, id_obs, id_var            ! IDs for variables
    integer :: dimids(2)                          ! Array of IDs
    integer :: fileid                             ! netcdf file id
    character(len=150) :: attstr                  ! String for attribute

  ! Print name file observation file
  IF (mype_filter == 0)  THEN
    WRITE(Message,'(a,a)') 'Initialize file for synthetic observations: ',trim(file_obs)
    CALL LOCAL_INFO(Caller,Message,Level=1)
  END IF

! *** Initialize file
  call nfcheck( NF90_CREATE(file_obs, 0, fileid))

  attstr  = 'Synthetic observations'
  call nfcheck( NF90_PUT_ATT(fileid, NF90_GLOBAL, 'title', trim(attstr))) 

  ! Define dimensions
  call nfcheck( NF90_DEF_DIM(fileid, 'dim_obs', dim_obs_p, dimid_obs))
  call nfcheck( NF90_DEF_DIM(fileid, 'steps',  NF90_UNLIMITED, dimid_iter))

  ! Define variables
  call nfcheck( NF90_DEF_VAR(fileid, 'timestep', NF90_INT, dimid_iter, id_time)) 

  dimids(1) = dimid_obs
  dimids(2) = dimid_iter
  call nfcheck( NF90_DEF_VAR(fileid, 'obs', NF90_DOUBLE, dimids, id_obs))
  call nfcheck( NF90_ENDDEF(fileid))

  ! Close file
  call nfcheck( NF90_CLOSE(fileid))

END SUBROUTINE init_file_syn_obs
! *****************************************************************************
! Write vector of synthetic observations to file
! *****************************************************************************
SUBROUTINE write_syn_obs(step, file_obs, dim_obs_p, observation_p)
  implicit none

! *** Arguments ***
  integer, intent(in) :: step        !< Currrent time step
  integer, intent(in) :: dim_obs_p   !< Dimension of full observation vector
  real(kind=dp), intent(in)    :: observation_p(dim_obs_p) !< Full observation vector
  character(len=*)    :: file_obs    !< Name of observation file

! *** Local variables ***
  CHARACTER(*), PARAMETER :: Caller = "write_syn_obs"
  integer :: iter                 ! file iteration counter
  integer :: fileid                       ! netcdf-ID of file
  integer :: id_time, id_obs  ! IDs for dimensions
  integer :: id_dimobs,id_dimsteps                 ! ID for dimension variable
  integer :: startv(2)         ! index arrays for writing

  integer :: dim_obs_store ! size of id_obs
  integer :: dimsteps      ! size of id_time

! ***********************************************
! *** Write observational information to file ***
! ***********************************************
  !* Open file
  call nfcheck( NF90_OPEN(trim(file_obs), NF90_WRITE, fileid))

  !* Read dimension dim_obs
  call nfcheck( NF90_INQ_DIMID(fileid, 'dim_obs', id_dimobs))
  call nfcheck( NF90_Inquire_dimension(fileid, id_dimobs, len=dim_obs_store))

  !* Consistency check
  IF (dim_obs_p /= dim_obs_store) then
    CALL PDAF_FATAL(Caller, &
          "Observation dimensions do not match :"//I2S(dim_obs_p)//' vs '//I2S(dim_obs_store))
  END IF

  !* get number of timesteps
  call nfcheck( NF90_INQ_DIMID(fileid, 'steps', id_dimsteps))
  call nfcheck( nf90_inquire_dimension(fileid, id_dimsteps, len=dimsteps))
  iter = dimsteps + 1

  IF (mype_filter == 0) THEN
    WRITE(Message,'(a)') 'Write synthetic observations at index '//I2S(iter)&
            //' to file: '//trim(file_obs)
    CALL LOCAL_INFO(Caller,Message,Level=1)
  END IF

  !* Write timestep
  call nfcheck( NF90_INQ_VARID(fileid, 'timestep', id_time))

  startv(1) = iter
  call nfcheck( NF90_PUT_VAR(fileid, id_time,step, start=startv(1:1)))

  !* Write vector of observations
  call nfcheck( NF90_INQ_VARID(fileid, 'obs', id_obs))

  startv(1) = 1
  startv(2) = iter
  call nfcheck( NF90_PUT_VAR(fileid, id_obs, observation_p, start=startv))

  !* Close file
  call nfcheck( NF90_CLOSE(fileid))
  
END SUBROUTINE write_syn_obs

! *****************************************************************************
! Read vector of synthetic observations from file
! *****************************************************************************
SUBROUTINE read_syn_obs(file_obs, dim_obs_p, observation_p, step)
  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_obs_p   !< Dimension of full observation vector
  real(kind=dp), intent(out)   :: observation_p(dim_obs_p) !< Full observation vector
  character(len=*)    :: file_obs    !< Name of observation file
  integer, intent(in) :: step        !< File time step to read

! *** Local variables ***
  CHARACTER(*), PARAMETER :: Caller = "read_syn_obs"
  integer :: iter                 ! file iteration counter
  integer :: fileid                       ! netcdf-ID of file
  integer :: id_obs                       ! IDs for dimensions
  integer :: id_time 
  integer :: id_dimobs                 ! ID for dimension variable
  integer :: id_dimsteps
  integer :: countv(2), startv(2)         ! index arrays for writing

  integer,allocatable :: obs_ts(:)
  integer :: dim_obs_store ! size of observation vector stored
  integer :: dimsteps
  integer :: i

! ***********************************************
! *** Read observations from file            ****          
! ***********************************************

  !* Open file
  call nfcheck( NF90_OPEN(trim(file_obs), NF90_NOWRITE, fileid))

  !* Read dimension dim_obs
  call nfcheck( NF90_INQ_DIMID(fileid, 'dim_obs', id_dimobs))
  call nfcheck( NF90_Inquire_dimension(fileid, id_dimobs, len=dim_obs_store))

  if (dim_obs_p /= dim_obs_store) then
      CALL PDAF_FATAL(Caller,'prolem with dim_obs: '//I2S(dim_obs_p)//' vs '//I2S(dim_obs_store))
  end if

  !* Read timesteps
  call nfcheck( NF90_INQ_DIMID(fileid, 'steps', id_dimsteps))
  call nfcheck( nf90_inquire_dimension(fileid, id_dimsteps, len=dimsteps))
  ALLOCATE(obs_ts(dimsteps))

  call nfcheck( NF90_INQ_VARID(fileid, 'timestep', id_time))
  call nfcheck( NF90_GET_VAR(fileid, id_time, obs_ts ))

  !* look for timestep corresponding to the current step
  iter=-1
  Do i=1,dimsteps
    IF (obs_ts(i) == step) THEN
       iter=i
       EXIT
    END IF
  End do
  IF (iter < 0 ) THEN
      CALL PDAF_FATAL(Caller,"Current timestep not found "//I2S(step))
  END IF

  IF (mype_filter == 0) THEN
    WRITE(Message,'(a)') 'Read synthetic observations at index '//I2S(iter)&
            //' from file: '//trim(file_obs)
    CALL LOCAL_INFO(Caller,Message,Level=1)
  END IF

  !* Read obs
  call nfcheck( NF90_INQ_VARID(fileid, 'obs', id_obs))
  startv(2) = iter
  countv(2) = 1
  startv(1) = 1
  countv(1) = dim_obs_store
  call nfcheck( NF90_GET_VAR(fileid, id_obs, observation_p, start=startv, count=countv))

  ! Close files
  call nfcheck( NF90_CLOSE(fileid))

  DEALLOCATE(obs_ts)
  
END SUBROUTINE read_syn_obs

END MODULE mod_synobs
