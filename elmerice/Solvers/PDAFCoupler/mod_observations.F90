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


MODULE mod_observations

  USE DefUtils
  USE PDAF,ONLY: obs_f,obs_l
  USE PDAFomi, ONLY: PDAFomi_set_debug_flag
  USE mod_assimilation, ONLY: omi_debug_flag

  IMPLICIT NONE

  REAL(kind=dp), PARAMETER :: DEFAULT_RMS=-999._dp  ! default obs rms
  !* Derived type for observations informations
  TYPE obsinfo_t
     INTEGER :: id
     CHARACTER(len=MAX_NAME_LEN) :: name    ! Name of observed variable
     INTEGER :: dim_obs_p                   ! PE-local number of observations
     INTEGER :: Obs_type                    ! Obs type see below
     REAL(kind=dp) :: rms=DEFAULT_RMS
     CHARACTER(len=MAX_NAME_LEN) :: input_file 
     INTEGER :: sid=-1
     INTEGER :: dim
  END TYPE obsinfo_t

  TYPE(obsinfo_t), ALLOCATABLE :: obsinfo(:) ! variable containing all obs info (size=obs_size)
  INTEGER :: obs_size   ! number of different observations
  INTEGER :: delt_obs   ! constant timestep interval between assimilations
  INTEGER :: t_assim    ! next assimilation timestep

  INTEGER, PARAMETER :: STATE_ON_NODES=1, &  ! State variable observed at nodes
                        OBS_ON_NODES=2       ! Model variable observed at nodes


  !* instances of observation data types used in the observation operator
  !* see callback_obs_pdafomi
  TYPE(obs_f), TARGET, ALLOCATABLE :: Allobs(:)      ! full observation
  TYPE(obs_l), TARGET, ALLOCATABLE :: Allobs_l(:)    ! local observation
!$OMP THREADPRIVATE(Allobs_l)

END MODULE mod_observations
