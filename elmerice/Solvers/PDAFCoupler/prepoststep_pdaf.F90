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
! * - Pre/Poststep routine for PDAF
! * - Used in all ensemble filters before and after the analysis
! * - Called by all filter pes
! * - Can be used for pre/post analysis of the forecast and analysis ensemble
! *
! * - Save ensemble mean and variance as elmer variables with suffix _meanf[f-a]
! *    and _var[f-a] for the forecast and analysys ensembles
! *    Only available on filter pes, i.e. also first model task (member)
! * - Compute observations diagnotics if do_omi_obsstats = TRUE
! *   Infos only printed to screeen for now .... 
! *
SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  USE DefUtils
  USE PDAFUtils, &     ! Parallelization
       ONLY: mype_filter, npes_filter, COMM_filter
  USE mod_assimilation, &      ! Assimilation variables
       ONLY: dim_state, do_omi_obsstats
  USE PDAF, &                  ! PDAF and PDAF-OMI diagnostic routines
       ONLY: PDAF_diag_variance,PDAF_diag_ensmean, PDAFomi_diag_obs_rmsd, PDAFomi_diag_diffstats
  USE mod_io
  USE mod_statevector_pdaf, ONLY : n_fields,sfields

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     !< Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   !< PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   !< PE-local dimension of observation vector
  REAL(kind=dp), INTENT(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  REAL(kind=dp), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  REAL(kind=dp), INTENT(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  INTEGER, INTENT(in) :: flag        !< PDAF status flag


! *** local variables ***
  INTEGER :: pdaf_status              ! status flag
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL(kind=dp) :: ens_stddev                  ! ensemble STDDEV = estimated RMS error
  REAL(kind=dp) :: mean_p(dim_p),variance_p(dim_p)
  ! Variables for parallelization - global fields
  ! Variables for observation diagnostics
  INTEGER :: nobs                     ! Number of active observation types
  REAL(kind=dp), POINTER :: obsrmsd(:)         ! Pointer to array of observation RMSDs
  REAL(kind=dp), POINTER :: obsstats(:,:)      ! Pointer to array  of observation statistics
  CHARACTER(*), PARAMETER :: Caller = "prepoststep_pdaf"
  INTEGER :: n,nf,node


! **********************
! *** INITIALIZATION ***
! **********************
  IF (mype_filter == 0) THEN
     CALL LOCAL_INFO(Caller,"------------------------------------------------------------",Level=1)
     IF (firsttime) THEN
        WRITE (Message, '(A)') 'Analyze initial state ensemble'
     ELSE
        IF (step<0) THEN
           WRITE (Message, '(Aa)') 'Analyze and write forecasted state ensemble'
        ELSE
           WRITE (Message, '(A)') 'Analyze and write assimilated state ensemble'
        END IF
     END IF
     CALL LOCAL_INFO(Caller,Message,Level=1)
  END IF


! ************************************************************
! *** Compute ensemble mean and standard deviation         ***
! *** (=RMS errors according to sampled covar matrix)      ***
! ************************************************************
! https://pdaf.awi.de/trac/wiki/DataAssimilationDiagnostics
  IF (dim_ens > 1) THEN
     CALL PDAF_diag_ensmean(dim_p, dim_ens, mean_p, ens_p, pdaf_status)
     !Rq. at initial time step>0
     !  but variables that are not created before first call to vtu output will not be saved??
     !   trick create the required variables now
     IF (firsttime) THEN
       CALL DISTRIBUTE(mean_p,dim_p,"_meanf")
       CALL DISTRIBUTE(mean_p,dim_p,"_meana")    
     ELSE        
       IF (step<0) THEN
         CALL DISTRIBUTE(mean_p,dim_p,"_meanf")
       ELSE
         CALL DISTRIBUTE(mean_p,dim_p,"_meana")      
       ENDIF
     END IF

     CALL PDAF_diag_variance(dim_p, dim_ens, mean_p, ens_p, variance_p, &
        ens_stddev, 0,1, COMM_filter, pdaf_status)
     IF (firsttime) THEN
       CALL DISTRIBUTE(variance_p,dim_p,"_varf")
       CALL DISTRIBUTE(variance_p,dim_p,"_vara")    
     ELSE
       IF (step<0) THEN
         CALL DISTRIBUTE(variance_p,dim_p,"_varf")
       ELSE
         CALL DISTRIBUTE(variance_p,dim_p,"_vara")      
       END IF
     ENDIF

     IF (mype_filter == 0) THEN
        WRITE (Message, '(A, es12.4)') &
             'ensemble standard deviation : ', ens_stddev
        CALL LOCAL_INFO(Caller,Message,Level=1)
     END IF
  END IF

! ***************************************
! *** Compute observation diagnostics ***
! ***************************************
! see https://pdaf.awi.de/trac/wiki/OMI_observation_diagnostics_PDAF3#PDAFomi_diag_stats
  IF (do_omi_obsstats) THEN
     ! Compute RMS deviation between observation and observed ensemble mean
     CALL PDAFomi_diag_obs_rmsd(nobs, obsrmsd, 1/(mype_filter+1))

     ! Compute statistics on deviation between observation and observed ensemble
     CALL PDAFomi_diag_diffstats(nobs, obsstats, 1/(mype_filter+1))
  END IF

  firsttime = .FALSE.

  IF (mype_filter == 0) THEN
     CALL LOCAL_INFO(Caller,"------------------------------------------------------------",Level=1)
  END IF

  CONTAINS

  !* Distribute ensemble means and variance to Elmer;
  !* As this is done by filter_pes corresponding to the first model task
  !* these variables will be available on the first member  
  SUBROUTINE DISTRIBUTE(state_p,dim_p,suffix)
  USE mod_statevector_pdaf, ONLY : n_fields,sfields
  IMPLICIT NONE
  TYPE(Variable_t), POINTER :: Var,RefVar
  TYPE(Mesh_t),POINTER :: Mesh
  REAL(kind=dp) :: state_p(dim_p)
  INTEGER :: dim_p
  CHARACTER(LEN=*) :: suffix
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  INTEGER, POINTER :: Perm(:)
  INTEGER :: nf,n,node
  
  Mesh => GetMesh()

  DO nf=1,n_fields
     VarName=TRIM(sfields(nf)%name)//TRIM(suffix)
     Var => VariableGet( Mesh % Variables,TRIM(VarName))
     IF (.NOT.ASSOCIATED(Var)) THEN
        RefVar => VariableGet( Mesh % Variables,TRIM(sfields(nf)%name),UnFoundFatal=.True.)
        CALL DefaultVariableAdd( TRIM(VarName), Perm = RefVar % Perm , Var = Var , Output=.True.)
     END IF

     Perm => Var % Perm
     DO n=1,sfields(nf)%dim
       node = sfields(nf)%StateToNode(n)
       IF (Perm(node) == 0) THEN
            CALL PDAF_FATAL(Caller,'Problem with permutation')
       END IF
        Var % Values ( Perm(node) ) = state_p( sfields(nf)%off + n )
     END DO

  END DO
  END SUBROUTINE DISTRIBUTE

END SUBROUTINE prepoststep_pdaf
