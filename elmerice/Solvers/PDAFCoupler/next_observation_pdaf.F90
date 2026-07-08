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
!>  - Initialize information for next forecast phase
!>  - Called by every model integration task
!>     before at the initialization and then for
!>     each subsequent forecasts phase
!>  
!>  - Assume a constant ts interval delt_obs between assimilation
!>  - delt_obs is set at initialisation from the .sif
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)
  USE DefUtils
  USE PDAFUtils
  USE mod_observations, ONLY : delt_obs,t_assim
  USE mod_io

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: stepnow  !< Number of the current time step
  INTEGER, INTENT(out) :: nsteps   !< Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   !< Whether to exit forecasting (1 for exit)
  REAL(kind=dp), INTENT(out)    :: time     !< Current model (physical) time
! *** Local Variables ***  
  CHARACTER(*), PARAMETER :: Caller = "PDAF_ELMER :: next_observation_pdaf"
  INTEGER :: TimeStep

  
! ***********************************************************************
! *** Sanity check TimeStep should correspond to the provided stepnow ***
! ***********************************************************************
  TimeStep=GetTimeStep()
  IF (TimeStep /= stepnow) CALL PDAF_FATAL(Caller,'Time steps error '//I2S(TimeStep)//&
                                                  ' vs '//I2S(stepnow))
! *************************************************************
! *** Determine number of time steps until next observation ***
! *** delt_obs set in init_obs see INIT_PDAF
! *************************************************************
  nsteps = delt_obs

! **********************************************************************
! *** t_assim is used to collect MASK and ObsEns before assimilation ***
! ***   see PDAF_ASSIMILATE                                          ***
! **********************************************************************  
  t_assim = stepnow + nsteps

! *********************************
! *** Set current physical time ***
! *********************************
   time = GetTime()

! *********************
! *** Set exit flag ***
! *********************
   doexit = 0

! ***************************
! *** OutPut informations ***
! ***************************   
 IF (mype_ens == 0) THEN
  WRITE (message, '(A)') 'Current time step : '//I2S(TimeStep)
  CALL LOCAL_INFO(Caller,message,level=1)
  WRITE (message, '(A)') 'Next observation at time step : '//I2S(t_assim)
  CALL LOCAL_INFO(Caller,message,level=1)
 ENDIF

END SUBROUTINE next_observation_pdaf
