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
! * Ensemble initialisation
! * The routine is called by all filter PEs (corresponding to task_id=1 in the standard configuration)
! * Only called by PDAF3_init
! * Here we do nothing and assume that each member will do its own (different) initialisation
! *  => The state is not distributed to elmer if TimeStep == 0 (see distribute_state_pdaf)
! *****************************************************************************
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, ens_p, flag)
  USE DefUtils
  USE PDAFUtils, ONLY : mype_filter
  USE mod_io
  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype              !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 !< Size of ensemble
  REAL(kind=dp), INTENT(out)   :: state_p(dim_p)          !< PE-local model state
  REAL(kind=dp), INTENT(out)   :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  REAL(kind=dp), INTENT(out)   :: ens_p(dim_p, dim_ens)   !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 !< PDAF status flag

! *** local variables ***
  REAL(kind=dp),parameter :: default=-1.0e12
  CHARACTER(*), PARAMETER :: Caller = "init_ens_pdaf"

  IF (mype_filter == 0 ) THEN
     CALL LOCAL_INFO(Caller,"---------------------------------------------------------------",level=1)     
     CALL LOCAL_INFO(Caller,"NOTHING DONE; assume independent initialisation",level=1)
  END IF

  ens_p=default
  state_p=default

END SUBROUTINE init_ens_pdaf
