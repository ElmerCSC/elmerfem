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
! *  Module containing  assimilation parameters
! *   All parameters are initialized with sensible default values
! *   All parameters can be changed via the .sif in the INIT_PDAF solver section
! *   Not everything tested...
MODULE mod_assimilation
  USE Types
  IMPLICIT NONE
  SAVE

!****************************************************************** 
! Variables **NOT** available as command line options!

  REAL, ALLOCATABLE :: coords_p(:,:)    !< Coordinates of process-local state vector entries
                                        !< needed to intiialize localization for LEnKF/ENSRF
  
  REAL(kind=dp) :: coords_l(3)                   !< Coordinates of local analysis domain
                                                 !< will be compared to thisobs % ncccord
!$OMP THREADPRIVATE(coords_l)                    


!******************************************************************
! **** PARAMETERS
!******************************************************************

!****************************************************************** 
! *** Activate omi_debug inofrmation
  INTEGER :: omi_debug_flag=0

!****************************************************************** 
! *** Perform twin experiment                                   ***
  LOGICAL :: twin_experiment=.FALSE.

! -----------------------------------------------------------------
! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! Settings for state vector size
  INTEGER :: dim_ens           !< Size of ensemble
  INTEGER :: dim_state         !< Global model state dimension
  INTEGER :: dim_state_p       !< Model state dimension for PE-local domain

! seedset for GenObs
  INTEGER :: seedset=1         ! integer set for GenOBS (1-20)

! Settings for model_error  
  LOGICAL :: model_error=.false.      !< Control application of model error
  REAL(kind=dp)  :: model_err_amp     !< Amplitude for model error

! Settings for observations - available as command line options
  INTEGER :: observe_ens=1     !< (0) apply H also to ensemble mean; (1) apply H only to ensemble states
  INTEGER :: type_obs_init=0   !< init obs. (0) before or (1) after call to prepostsstep
  LOGICAL :: do_omi_obsstats=.True. !< Whether to let OMI compute observation statistics

! General control of PDAF - available as command line options
  INTEGER :: screen=0     !< Control verbosity of PDAF
                          !< * (0) no outputs
                          !< * (1) progress info
                          !< * (2) add timings
                          !< * (3) debugging output
  INTEGER :: filtertype=4 !< Select filter algorithm:
                          !<   * SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !<   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
                          !<   LKNETF (11), PF (12), ENSRF/EAKF (13), GENOBS (100), 3DVAR (200)
  INTEGER :: subtype=0    !< Subtype of filter algorithm
                          !<   * SEIK:
                          !<     (0) ensemble forecast; new formulation
                          !<     (1) ensemble forecast; old formulation
                          !<     (2) SEIK with ensemble transformation
                          !<     (10) fixed error space basis
                          !<     (11) fixed state covariance matrix
                          !<   * LSEIK:
                          !<     (0) ensemble forecast;
                          !<     (2) LSEIK with ensemble transformation
                          !<     (10) fixed error space basis
                          !<     (11) fixed state covariance matrix
                          !<   * ETKF:
                          !<     (0) ETKF using T-matrix like SEIK
                          !<     (1) ETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LETKF:
                          !<     (0) LETKF using T-matrix like SEIK
                          !<     (1) LETKF following Hunt et al. (2007)
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * EnKF:
                          !<     (0) analysis for large observation dimension
                          !<     (1) analysis for small observation dimension
                          !<   * LEnKF:
                          !<     (0) standard analysis
                          !<   * ESTKF:
                          !<     (0) standard ESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to SEIK subtypes 2/3
                          !<   * LESTKF:
                          !<     (0) standard LESTKF 
                          !<       There are no fixed basis/covariance cases, as
                          !<       these are equivalent to LSEIK subtypes 2/3
                          !<   * NETF:
                          !<     (0) standard NETF 
                          !<   * LNETF:
                          !<     (0) standard LNETF 
                          !<   * LKNETF:
                          !<     (0) HNK: 2-step LKNETF with NETF before LETKF
                          !<     (1) HKN: 2-step LKNETF with LETKF before NETF
                          !<     (2) HSync: LKNETF synchronous
                          !<   * PF:
                          !<     (0) standard PF 
                          !<   * ENSRF/EAKF:
                          !<     (0) ENSRF with serial observation processing
                          !<     (1) EAKF with loca least square regression
                          !<   * 3D-Var:
                          !<     (0) parameterized 3D-Var
                          !<     (1) 3D Ensemble Var using LESTKF for ensemble update
                          !<     (2) 3D Ensemble Var using ESTKF for ensemble update
                          !<     (3) hybrid 3D-Var using LESTKF for ensemble update
                          !<     (4) hybrid 3D-Var using ESTKF for ensemble update
  INTEGER :: type_iau=0     !< Type of incremental updating:
                          !<     (0) no IAU
                          !<     (1) constant IAU weight
                          !<     (2) linear increase/decrease with maimum in middle of period
                          !<     (3) Null IAU: initialize increments arrays, but do not add increment
  INTEGER :: steps_iau    !< Number of time steps over which IAU is applied
  INTEGER :: dim_lag      !< Number of time instances for smoother

! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget=0 !< Type of forgetting factor
                           !<  SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                           !<   (0) fixed
                           !<   (1) global adaptive
                           !<   (2) local adaptive for LSEIK/LETKF/LESTKF
                           !<  NETF/LNETF/PF
                           !<   (0) apply inflation on forecast ensemble
                           !<   (2) apply inflation on analysis ensemble
  REAL(kind=dp) :: forget=1.0   !< Forgetting factor for filter analysis
  INTEGER :: dim_bias     !< dimension of bias vector

!    ! All localized filters
  REAL(kind=dp) :: cradius=-HUGE(1._dp)   !< Cut-off radius for local observation domain
  INTEGER :: locweight=0     !< * Type of localizing weighting of observations
                             !<   (0) constant weight of 1
                             !<   (1) exponentially decreasing with SRADIUS
                             !<   (2) use 5th-order polynomial
                             !<   (3) regulated localization of R with mean error variance
                             !<   (4) regulated localization of R with single-point error variance
  REAL(kind=dp) :: sradius=-HUGE(1._dp)    !< Support radius for 5th order polynomial
                           !<   or radius for 1/e for exponential weighting
!    ! ENKF
  INTEGER :: rank_ana_enkf=0 !< Rank to be considered for inversion of HPH in analysis of EnKF
                           !<  (0) for analysis w/o eigendecomposition
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF/NETF/LNETF/LKNETF
  INTEGER :: type_trans=0  !< Type of ensemble transformation 
                           !< * SEIK/LSEIK: 
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ETKF/LETKF with subtype=4: 
                           !< (0) use deterministic symmetric transformation
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T 
                           !< * ESTKF/LESTKF:
                           !< (0) use deterministic omega
                           !< (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           !< (2) use product of (0) with random orthonormal matrix with
                           !<     eigenvector (1,...,1)^T
                           !< * NETF/LNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
                           !< * LKNETF:
                           !< (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           !< (1) use identity transformation
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt=0   !< * Type of the transform matrix square-root 
                           !<   (0) symmetric square root
                           !<   (1) Cholesky decomposition
!    ! NETF/LNETF/PF
  INTEGER :: type_winf=0   !< Set weights inflation: 
                           !<   (0) no weights inflation
                           !<   (1) use N_eff/N>limit_winf
  REAL(kind=dp) :: limit_winf=HUGE(1._dp)    !< Limit for weights inflation: N_eff/N>limit_winf
!    ! hybrid LKNETF
  INTEGER :: type_hyb=0      !< * Type of hybrid weight:
                           !<   (0) use fixed hybrid weight hyb_gamma
                           !<   (1) use gamma_lin: (1 - N_eff/N_e)*hyb_gamma
                           !<   (2) use gamma_alpha: hybrid weight from N_eff/N>=hyb_gamma
                           !<   (3) use gamma_ska: 1 - min(s,k)/sqrt(hyb_kappa) with N_eff/N>=hyb_gamma
                           !<   (4) use gamma_sklin: 1 - min(s,k)/sqrt(hyb_kappa) >= 1-N_eff/N>=hyb_gamma
  REAL(kind=dp) :: hyb_gamma=1.0 !< Hybrid filter weight for state (1.0: LETKF, 0.0 LNETF)
  REAL(kind=dp) :: hyb_kappa=HUGE(1._dp)    !< Hybrid norm for using skewness and kurtosis

!    ! Particle filter
  INTEGER :: pf_res_type=1 !< * Resampling type for PF
                           !<   (1) probabilistic resampling
                           !<   (2) stochastic universal resampling
                           !<   (3) residual resampling        
  INTEGER :: pf_noise_type=0  !< * Resampling type for PF
                              !<   (0) no perturbations, (1) constant stddev, 
                              !<   (2) amplitude of stddev relative of ensemble variance
  REAL(kind=dp) :: pf_noise_amp=HUGE(1._dp)    !< Noise amplitude (>=0.0, only used if pf_noise_type>0)

END MODULE mod_assimilation
