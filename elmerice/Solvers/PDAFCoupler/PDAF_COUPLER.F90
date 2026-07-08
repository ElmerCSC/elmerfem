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
! ******************************************************************************
! *
! *  Authors: Fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 4 Jun. 2026
! * 
! *****************************************************************************
! * COUPLING WITH PDAF (Parallel Data Assimilation Framework)
! *  https://pdaf.awi.de/trac/wiki/WikiStart
! * 
! * SOLVERS FOR THE ONLINE COUPLING :
! *  - INIT_PDAF : perform PDAF initilisation; executed <Before Simulation>; i.e. 
! *      before model time integration
! *  - PDAF_ASSIMILATE : Assimilation; executed <After TimeStep>; i.e. 
! *      at the end of each timestep 
! * OTHER SUBROUTINES:
! *  - init_parameters : pdaf parameters initialisation
! *  - init_obs :        observation parameters initialisation
! *****************************************************************************/

! *****************************************************************************
! * INIT_PDAF
! * to see : dim_lag,model_error, model_err_amp not set
! *****************************************************************************/
! *****************************************************************************
!> Adapted from init_pdaf in PDAF templates and examples
! *****************************************************************************
SUBROUTINE INIT_PDAF_Init0( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !--------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  ! ************************************************
  ! ***   Has to be executed before simulation   ***
  ! ************************************************
  CALL ListAddString(SolverParams, 'Exec Solver', 'before simulation' )

END SUBROUTINE INIT_PDAF_Init0

SUBROUTINE INIT_PDAF( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE PDAFUtils           ! Parallelization variables - Elmer Module
  USE PDAF                ! PDAF interface definitions
  USE mod_assimilation
  USE mod_statevector_pdaf, & ! Variables for statevector
        ONLY: setup_statevector,distribute_state_pdaf
  USE mod_io
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !--------------------------------------------------------------------------

  ! *** Local variables ***
  CHARACTER(*), PARAMETER :: Caller = "init_pdaf"
  INTEGER :: filter_param_i(2)       ! Integer parameter array for filter
  REAL(kind=dp) :: filter_param_r(1) ! Real parameter array for filter
  INTEGER :: status_pdaf             ! PDAF status flag

  ! *** External subroutines ***
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Time to next observation
              prepoststep_pdaf         ! Pre/poststep routine
  
! ***************************
! ***   Initialize PDAF   ***
! ***************************
  IF (mype_ens == 0) THEN
     WRITE (Message,'(A)')  '+++++++++++++++++++++++++++++++++++++++++++'     
     CALL LOCAL_INFO(Caller,Message,Level=1)
     WRITE (Message,'(A)')  '++  INITIALIZE PDAF - ONLINE MODE   +++++++'
     CALL LOCAL_INFO(Caller,Message,Level=1)
     WRITE (Message,'(A)')  '+++++++++++++++++++++++++++++++++++++++++++'     
     CALL LOCAL_INFO(Caller,Message,Level=1)
  END IF

! ***************************
! *** Ensemble settings   ***
! ***************************
  dim_ens = n_modeltasks ! full parallelisation dim_ens = n_modeltasks 
                         !  initialised during mpi initialisation 

! ***************************
! *** Options for DA methods
! ***************************
  call init_parameters()

  IF (filtertype==PDAF_DA_3DVAR) THEN
    WRITE (Message,'(A)') 'ERROR 3DVAR not available'
    CALL PDAF_FATAL(Caller,Message)
  END IF

! *** Possibly activate PDAF-OMI observation statistics ***
  IF (do_omi_obsstats) CALL PDAFomi_set_obs_diag(1)

! ***************************
! *** Define state vector ***
! ***************************
  call setup_statevector(dim_state, dim_state_p, screen)

! *******************************
! *** Observations parameters ***
! *******************************
   call init_obs()

! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** For all filters, PDAF_init is first called    ***
! *** specifying only the required parameters.      ***
! *** Further settings are done afterwards using    ***
! *** calls to PDAF_set_iparam & PDAF_set_rparam.   ***
! *****************************************************
! *** see https://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF
! ***   for all filter options

  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_r(1) = forget      ! Forgetting factor

  CALL PDAF3_init(filtertype, subtype, 0, &
       filter_param_i, 2,&
       filter_param_r, 1, &
       init_ens_pdaf, screen, status_pdaf)

  ! *** Additional parameter specifications ***
  ! *** -- These are all optional --        ***
  IF (filtertype==PDAF_DA_GENOBS) THEN
    CALL PDAF_set_iparam(3, seedset,status_pdaf)
  ELSE  
    ! Generic settings
    CALL PDAF_set_iparam(5, type_forget, status_pdaf)      ! Type of forgetting factor
    CALL PDAF_set_iparam(6, type_trans, status_pdaf)       ! Type of ensemble transformation
    CALL PDAF_set_iparam(7, type_sqrt, status_pdaf)        ! Type of transform square-root (SEIK-sub4/ESTKF)
    CALL PDAF_set_iparam(8, observe_ens, status_pdaf)      ! Whether to apply observation operator to ensemble mean
    CALL PDAF_set_iparam(9, type_obs_init, status_pdaf)    ! Initialize observation before or after call to prepoststep

    ! Setting for EnKF and LEnKF
    IF (filtertype==PDAF_DA_ENKF .OR. filtertype==PDAF_DA_LENKF) THEN
      CALL PDAF_set_iparam(4, rank_ana_enkf, status_pdaf) ! Rank of EVD in EnKF (0 for direct inversion)
    END IF

    ! Settings for NETF and LNETF
    IF (filtertype==PDAF_DA_NETF .OR. filtertype==PDAF_DA_LNETF) THEN
      CALL PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
      CALL PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
      CALL PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
      CALL PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
    END IF

    ! Settings for particle filter PF
    IF (filtertype==PDAF_DA_PF) THEN
      CALL PDAF_set_iparam(4, pf_noise_type, status_pdaf) ! Perturbation type
      CALL PDAF_set_iparam(6, pf_res_type, status_pdaf)   ! Resampling type
      CALL PDAF_set_iparam(7, type_winf, status_pdaf)     ! Type of weights inflation
      CALL PDAF_set_rparam(2, limit_winf, status_pdaf)    ! Limit for weights inflation
      CALL PDAF_set_rparam(3, pf_noise_amp, status_pdaf)  ! Noise amplitude
    END IF

    ! Settings for hybrid filter LKNETF
    IF (filtertype==PDAF_DA_LKNETF) THEN
       CALL PDAF_set_iparam(4, type_hyb, status_pdaf)      ! Choice of hybrid rule
       CALL PDAF_set_rparam(2, hyb_gamma, status_pdaf)     ! Hybrid filter weight for state
       CALL PDAF_set_rparam(3, hyb_kappa, status_pdaf)     ! Normalization factor for hybrid weight 
    END IF

  END IF

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (Message,'(A)') &
          'ERROR '//I2S(status_pdaf)// &
          ' in initialization of PDAF - stopping! (PE '//I2S(mype_ens)//')'
     CALL PDAF_FATAL(Caller,Message)
  END IF

! **********************
! *** Initialize IAU ***
! **********************
  CALL PDAF_iau_init(type_iau, steps_iau, status_pdaf)

! **********************************
! *** Prepare ensemble forecasts ***
! **********************************
  CALL PDAF_init_forecast(next_observation_pdaf, distribute_state_pdaf, &
                          prepoststep_pdaf, status_pdaf)

END SUBROUTINE INIT_PDAF

! *****************************************************************************
! * PDAF_ASSIMILATE
! *****************************************************************************
!>  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the routine of PDAF (PDAF3_assimilate), which checks
!! whether the forecast phase is completed. If so, the analysis step
!! is computed inside PDAF.
! *****************************************************************************
!> Adapted from assimilate_pdaf in PDAF templates and examples
! *****************************************************************************
SUBROUTINE PDAF_ASSIMILATE_Init0( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !--------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
  !------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  ! ************************************************
  ! ***   Has to be executed After TimeStep      ***
  ! ************************************************
  CALL ListAddString(SolverParams, 'Exec Solver', 'After TimeStep' )

END SUBROUTINE PDAF_ASSIMILATE_Init0

SUBROUTINE PDAF_ASSIMILATE( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE PDAFUtils                   ! Elmer Module: Parallelization variables
  USE PDAF, &                     ! PDAF  Module: interface definitions
       ONLY: PDAF3_assimilate, PDAF3_generate_obs, &
             PDAF_abort, PDAF_DA_GENOBS
  USE mod_assimilation            ! Local Module: DA Variables
  USE mod_statevector_pdaf, &     ! Local Module: State collection/distribution
       ONLY:  collect_state_pdaf,distribute_state_pdaf,collect_elmer_state, &
              n_ofields,ofields,n_mfields,mfields,&
              COLLECT_ObsEns,COLLECT_MASK

  USE mod_io                      ! Local Module for io
  USE mod_observations, ONLY : t_assim
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !--------------------------------------------------------------------------
! *** Local variables ***
  INTEGER :: status_pdaf          ! PDAF status flag
  CHARACTER(*), PARAMETER :: Caller = "assimilate_pdaf"
  INTEGER :: ts
  INTEGER :: j


! *** External subroutines ***
   EXTERNAL :: next_observation_pdaf, &       ! Provide time step of next observation
               prepoststep_pdaf               ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &          ! Provide number of local analysis domains
              init_dim_l_pdaf                 ! Initialize state dimension for local analysis domain

  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: init_dim_obs_pdafomi, &        ! Get dimension of full obs. vector for PE-local domain
              obs_op_pdafomi, &              ! Obs. operator for full obs. vector for PE-local domain
              init_dim_obs_l_pdafomi         ! Get dimension of obs. vector for local analysis domain

  ! Subroutine used for generating observations
  EXTERNAL :: get_obs_f_pdaf                 ! Get vector of synthetic observations from PDAF

  INTEGER :: dim_o,dim_m
  REAL(Kind=dp),DIMENSION(:),ALLOCATABLE :: tmp_p


  ts = GetTimeStep()

  ! *****************************************************************
  ! *** collect ensemble of "observations" variables and masks   ****
  ! *** in the filter pes                                        ****
  ! *** t_assim is set in next_observation_pdaf                  ****
  ! *****************************************************************
  IF (ts == t_assim) THEN
     IF (n_ofields > 0) THEN     
        dim_o=0     
        DO j=1,n_ofields
           dim_o = dim_o + ofields(j)%dim
        END DO    
        ALLOCATE(tmp_p(dim_o))
        CALL collect_elmer_state(dim_o,tmp_p,n_ofields,ofields)     
        CALL COLLECT_ObsEns(tmp_p,dim_o)
        DEALLOCATE(tmp_p)
     END IF
     IF (n_mfields > 0 ) THEN
        dim_m = 0
        DO j=1,n_mfields
           dim_m = dim_m + mfields(j)%dim
        END DO
        ALLOCATE(tmp_p(dim_m))
        CALL collect_elmer_state(dim_m,tmp_p,n_mfields,mfields)
        CALL COLLECT_MASK(tmp_p,dim_m)
        DEALLOCATE(tmp_p)
     END IF
  END IF
! *********************************
! *** Call assimilation routine ***
! *********************************
  IF (filtertype /= PDAF_DA_GENOBS) THEN
     ! Call universal PDAF3 interface routine
     CALL PDAF3_assimilate(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, &
          init_n_domains_pdaf, init_dim_l_pdaf, init_dim_obs_l_pdafomi, &
          prepoststep_pdaf, next_observation_pdaf, status_pdaf)
  ELSE
     ! Observation generation has its own OMI interface routine
     CALL PDAF3_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_f_pdaf, &
          prepoststep_pdaf, next_observation_pdaf, status_pdaf)
  END IF
! *************************
! *** Check status flag ***
! *************************
  IF (status_pdaf /= 0) THEN
     WRITE (Message,'(A)') &
          'ERROR '//I2S(status_pdaf)// &
          ' in PDAF3_assimilate - stopping! (PE '//I2S(mype_ens)//')'
     CALL PDAF_FATAL(Caller,Message)
  END IF


END SUBROUTINE PDAF_ASSIMILATE

! *****************************************************************************
! * init_parameters
! *****************************************************************************
!> Read parameters from .sif files in the PDAF_INIT solver section
! *****************************************************************************
!> Adapted from pdaf_parse in PDAF templates and examples
! *****************************************************************************
SUBROUTINE init_parameters()
  USE DefUtils
  USE mod_assimilation
   
  IMPLICIT NONE

  INTERFACE GetParam 
        procedure :: GetIParam,GetRParam,GetLParam
  END INTERFACE GetParam

  TYPE(ValueList_t), POINTER :: SolverParams
  INTEGER :: iparam
  REAL(kind=dp) :: rparam
  CHARACTER(len=MAX_NAME_LEN) :: handle

  SolverParams => GetSolverParams()

! **********************************
! *** Get options from .sif
! **********************************
! 
  handle = 'omi_debug_flag'     ! omi_debug_flag
  CALL GetParam(omi_debug_flag,handle)

  handle = 'twin_experiment'     ! twin experiment
  CALL GetParam(twin_experiment,handle)

  handle = 'seedset'             ! seedset for GenObs - 1->20
  CALL GetParam(seedset,handle)

  ! Observation settings
  handle = 'observe_ens'             ! (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  CALL GetParam(observe_ens,handle)

  handle = 'type_obs_init'           ! init obs. (0) before or (1) after call to prepostsstep
  CALL GetParam(type_obs_init,handle)

  handle = 'do_omi_obsstats'         ! Whether to let PDAF-OMI compute observation statistics
  CALL GetParam(do_omi_obsstats,handle)

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL GetParam(model_error,handle)

  handle = 'model_err_amp'           ! Amplitude of model error
  CALL GetParam(model_err_amp,handle)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL GetParam(screen,handle)

  ! Filters
  handle = 'filtertype'              ! Choose filter algorithm
  CALL GetParam(filtertype,handle)
  handle = 'subtype'                 ! Set subtype of filter
  CALL GetParam(subtype,handle)

  ! Control IAU
  handle = 'type_iau'                ! Set whether to use incremental updating
  CALL GetParam(type_iau,handle)
  handle = 'steps_iau'               ! Number of time steps over which IAU is applied
  CALL GetParam(steps_iau,handle)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  CALL GetParam(dim_lag,handle)

  ! Filter-specific settings
  handle = 'forget'                  ! Set forgetting factor
  CALL GetParam(forget,handle)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL GetParam(type_forget,handle)

  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  CALL GetParam(type_trans,handle)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL GetParam(type_sqrt,handle)
  handle = 'rank_ana_enkf'           ! Set rank for pseudo inverse in EnKF
  CALL GetParam(rank_ana_enkf,handle)

  ! Settings for localization in LSEIK/LETKF
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL GetParam(cradius,handle)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL GetParam(locweight,handle)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or radius for 1/e in exponential weighting
  CALL GetParam(sradius,handle)

  ! Settings for nonlinear filters
  handle = 'pf_res_type'             ! Resampling type for particle filter
  CALL GetParam(pf_res_type,handle)        
  handle = 'pf_noise_type'           ! Type of perturbing noise in PF
  CALL GetParam(pf_noise_type,handle)        
  handle = 'pf_noise_amp'            ! Amplitude of perturbing noise in PF
  CALL GetParam(pf_noise_amp,handle)        
  handle = 'type_winf'               ! Set type of weights inflation in NETF/LNETF
  CALL GetParam(type_winf,handle)
  handle = 'limit_winf'              ! Set limit for weights inflation
  CALL GetParam(limit_winf,handle)

  ! Hybrid weights for LKNETF
  handle = 'type_hyb'                ! Set type of hybrid weight
  CALL GetParam(type_hyb,handle)
  handle = 'hyb_gamma'               ! Set hybrid filter weight for state (1.0 LETKF, 0.0 LNETF)
  CALL GetParam(hyb_gamma,handle)
  handle = 'hyb_kappa'               ! Set hybrid norm (>1.0)
  CALL GetParam(hyb_kappa,handle)

  CONTAINS 
    !*******************************************************
    !*** Generic routines to get sif parameters
    !*******************************************************
    SUBROUTINE GetIParam(param,handle)
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: param
    character(len=MAX_NAME_LEN),INTENT(IN) :: handle
    INTEGER :: iparam
    LOGICAL :: Found
      iparam = ListGetInteger(SolverParams,handle,Found)
      IF (Found) param=iparam
    END SUBROUTINE GetIParam

    SUBROUTINE GetRParam(param,handle)
    IMPLICIT NONE
    REAL(kind=dp),INTENT(INOUT) :: param
    character(len=MAX_NAME_LEN),INTENT(IN) :: handle
    REAL(kind=dp) :: rparam
    LOGICAL :: Found
      
      rparam = ListGetConstReal(SolverParams,handle,Found)
      IF (Found) param=rparam

    END SUBROUTINE GetRParam

    SUBROUTINE GetLParam(param,handle)
    IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: param
    character(len=MAX_NAME_LEN),INTENT(IN) :: handle
    LOGICAL :: lparam
    LOGICAL :: Found
      
      lparam = ListGetLogical(SolverParams,handle,Found)
      IF (Found) param=lparam

    END SUBROUTINE GetLParam

END SUBROUTINE init_parameters

! ***************************************************************************************
! * init_obs
! ***************************************************************************************
!> Initialise the obsinfo object that contains informations about assimilated observations
SUBROUTINE init_obs()
  USE mod_observations
  USE mod_statevector_pdaf, ONLY : n_fields,sfields
  USE mod_io
  IMPLICIT NONE
  CHARACTER(*), PARAMETER :: Caller = "init_pdaf"
  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  REAL(kind=dp) :: rparam
  INTEGER :: i,j
  LOGICAL :: Found
  INTEGER :: obs_state_size,obs_var_size
  INTEGER :: id

  !*********************************************************************
  !* SolverParams (section INIT_PDAF Solver)
  !*********************************************************************
  SolverParams => GetSolverParams()

  !*********************************************************************
  !* (Constant) TimeStep interval between assimilations
  !*********************************************************************
  delt_obs = ListGetInteger(SolverParams,'delt_obs',UnFoundfatal=.True.)

  !*********************************************************************
  !* Observed variables corresponding to states variables
  !*********************************************************************
  DO j=1,1000
    VarName = ListGetString( SolverParams ,'Observed State '//I2S(j), Found )
    IF( .NOT. Found ) EXIT
  END DO
  obs_state_size=j-1

  !*********************************************************************
  !* Observed variables not corresponding to  states variables
  !*  and use in the observation operator e.g. velocity
  !*********************************************************************
  DO j=1,1000
    VarName = ListGetString( SolverParams ,'Observed Variable '//I2S(j), Found )
    IF( .NOT. Found ) EXIT
  END DO
  obs_var_size=j-1

  obs_size = obs_state_size + obs_var_size
  ALLOCATE(obsinfo(obs_size))

  !*********************************************************************
  !* Get required informations about observations
  !*********************************************************************
  DO j=1,obs_state_size
     obsinfo(j) % id = j
     !* Variable name
     obsinfo(j) % name = ListGetString( SolverParams ,'Observed State '//I2S(j), UnFoundFatal=.True. )
     !* coordinates dimension 1D,2D,3D
     obsinfo(j) % dim =  ListGetInteger( SolverParams ,'Observed State dimension '//I2S(j), UnFoundFatal=.True. )
     !* input file
     obsinfo(j) % input_file = ListGetString( SolverParams ,'Observed State input file '//I2S(j), UnFoundFatal=.True. )
     !* assume observation is given at nodes
     obsinfo(j) % obs_type = STATE_ON_NODES
     !* assume uniform rms for now
     rparam = ListGetConstReal(SolverParams, 'Observed State rms '//I2S(j),UnFoundFatal=.True.)
     obsinfo(j) % rms = rparam
     !* get the corresponding state variable id (id in sfields)
     DO i=1,n_fields
       IF (TRIM(obsinfo(j) % name) == TRIM(sfields(i) % name)) &
             obsinfo(j) % sid = i
     END DO
     IF (obsinfo(j) % sid < 1) THEN
             WRITE(Message,'(A)') 'No corresponding state variable for obs variable :'//TRIM(obsinfo(j) % name) 
             CALL PDAF_FATAL(Caller,Message)
     END IF
  END DO

  !* same for Observed Variables
  DO j=1,obs_var_size
     id = j + obs_state_size
     obsinfo(id) % id = id
     !* Variable name
     obsinfo(id) % name = ListGetString( SolverParams ,'Observed Variable '//I2S(j), UnFoundFatal=.True. )
     !* coordinates dimension 1D,2D,3D
     obsinfo(id) % dim =  ListGetInteger( SolverParams ,'Observed Variable dimension '//I2S(j), UnFoundFatal=.True. )
     !* input file
     obsinfo(id) % input_file = ListGetString( SolverParams ,'Observed Variable input file '//I2S(j), UnFoundFatal=.True. )
     !* assume observation is given at nodes
     obsinfo(id) % obs_type = OBS_ON_NODES
     !* assume uniform rms for now
     rparam = ListGetConstReal(SolverParams, 'Observed Variable rms '//I2S(j),UnFoundFatal=.True.)
     obsinfo(id) % rms = rparam
     !* get the corresponding observed variable id (id in ofields)
     obsinfo(id) % sid = j
  END DO
END SUBROUTINE init_obs

