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
! * Called by filter PDAF_DA_GENOBS (used for synthetic observations experiments
! *   in twin experiments)
! * Called by every filter_pe
! * Called at each "assimilation step" when synthetic observations are generated
! * The Observations are saved in netcdf file - 1 file/observation type and filter_pe
! * see mod_synobs.F90 
SUBROUTINE get_obs_f_pdaf(step, dim_obs, observation)
  USE DefUtils
  USE PDAFUtils, &
       ONLY: mype_filter
  use mod_synobs, &
       only: write_syn_obs
  USE mod_observations, &
       ONLY: obs_size,obsinfo
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_obs              ! Dimension of obs. vector
  REAL(kind=dp), INTENT(in)    :: observation(dim_obs) ! Observation vector
  CHARACTER(LEN=MAX_PATH_LEN) :: file_synobs

  INTEGER :: i
  INTEGER :: id 
  INTEGER :: odim
  INTEGER :: count,st,en

  count=1
  DO i=1,obs_size
     id=obsinfo(i)%id
     file_synobs = 'syn_obs_'//TRIM(obsinfo(i)%name)//'_'//I2S(mype_filter)//'.nc'
     odim=obsinfo(i)%dim_obs_p
     st=count
     en=count+odim-1
     call write_syn_obs(step, file_synobs,odim,observation(st:en))
     count = count + odim
  END DO

END SUBROUTINE get_obs_f_pdaf

