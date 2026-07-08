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
! * Generic callback routines for observations
! * Called by all filter pes
! * We rely on specific routines for different observations types (only observations at nodes for now)
! *****************************************************************************
! *****************************************************************************
! * Intialise the observations
! *****************************************************************************
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)
   USE mod_observations
!  USE  PDAFUtils, ONLY : task_id
!  USE PDAF, ONLY : PDAF_get_memberid,PDAF_get_obsmemberid
  ! Include functions for different observations
  USE obs_AtNodes_pdafomi, ONLY: init_dim_obs_AtNodes
  USE mod_io
  USE mod_statevector_pdaf, ONLY: sfields,ofields

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_type ! Observation dimensions
  CHARACTER(*), PARAMETER :: Caller = "init_dim_obs_pdafomi"
  INTEGER :: i
  LOGICAL,SAVE :: FirstTime=.TRUE.

  !* set debug flag (will print many informations if activated)
  CALL PDAFomi_set_debug_flag(omi_debug_flag)

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************
  !!! the size might change if we don't assimilate the same number o obs...?
  IF (ALLOCATED(Allobs)) DEALLOCATE(Allobs)
  allocate(Allobs(obs_size))

  ! Initialize number of observations
  dim_obs_type = 0
  dim_obs = 0

  ! Call observation-specific routines
  DO i=1,obs_size

     SELECT CASE(obsinfo(i)%Obs_type)

      CASE(STATE_ON_NODES)
        call init_dim_obs_AtNodes(step, dim_obs_type,Allobs(i),obsinfo(i),sfields,FirstTime)

      CASE(OBS_ON_NODES)
        call init_dim_obs_AtNodes(step, dim_obs_type,Allobs(i),obsinfo(i),ofields,FirstTime)

      CASE DEFAULT
        CALL PDAF_FATAL(Caller,'Obs_type not available '//I2S(obsinfo(i)%Obs_type))
        
     END SELECT

     dim_obs = dim_obs + dim_obs_type
  END DO

  FirstTime=.FALSE.

END SUBROUTINE init_dim_obs_pdafomi

! *****************************************************************************
! * Observation operator
! *****************************************************************************
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)
  USE mod_observations
  USE PDAFUtils
  ! Include functions for different observations
  USE obs_AtNodes_pdafomi, ONLY: obs_op_state_AtNodes,obs_op_obs_AtNodes
  USE PDAF
  USE mod_io

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL(KIND=dp), INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL(KIND=dp), INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state

  CHARACTER(*), PARAMETER :: Caller = "obs_op_pdafomi"
  INTEGER :: i

  !* set debug flag (will print many informations if activated)
  CALL PDAFomi_set_debug_flag(omi_debug_flag)

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************
  DO i=1,obs_size

     SELECT CASE(obsinfo(i)%Obs_type)

      CASE(STATE_ON_NODES)
        CALL obs_op_state_AtNodes(dim_p, dim_obs, state_p, ostate,Allobs(i))

      CASE(OBS_ON_NODES)
        CALL obs_op_obs_AtNodes(dim_p, dim_obs, state_p, ostate,Allobs(i))

      CASE DEFAULT
        CALL PDAF_FATAL(Caller,'Obs_type not available '//I2S(obsinfo(i)%Obs_type))

     END SELECT
  END DO


END SUBROUTINE obs_op_pdafomi

! *****************************************************************************
! * Intialise the local observations
! *****************************************************************************
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)
  USE mod_observations
  ! Include functions for different observations
  USE obs_AtNodes_pdafomi, ONLY: init_dim_obs_l_AtNodes
  USE mod_io

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

  CHARACTER(*), PARAMETER :: Caller = "init_dim_obs_l_pdafomi"
  INTEGER :: i

  !* set debug flag (will print many informations if activated)
  CALL PDAFomi_set_debug_flag(omi_debug_flag)
! **********************************************
! *** Initialize local observation dimension ***
! **********************************************
  IF (ALLOCATED(Allobs_l)) DEALLOCATE(Allobs_l)
  ALLOCATE(Allobs_l(obs_size))

   DO i=1,obs_size

     SELECT CASE(obsinfo(i)%Obs_type)

      CASE(STATE_ON_NODES)
       CALL init_dim_obs_l_AtNodes(domain_p, step, dim_obs, dim_obs_l,Allobs(i),Allobs_l(i))

      CASE(OBS_ON_NODES)
       CALL init_dim_obs_l_AtNodes(domain_p, step, dim_obs, dim_obs_l,Allobs(i),Allobs_l(i))     

      CASE DEFAULT
        CALL PDAF_FATAL(Caller,'Obs_type not available '//I2S(obsinfo(i)%Obs_type))

     END SELECT
  END DO

END SUBROUTINE init_dim_obs_l_pdafomi
