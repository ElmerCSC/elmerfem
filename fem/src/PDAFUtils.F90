!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Fabien Gillet-Chaulet
! *  Email:   fabien.gillet-chaulet@univ-grenoble-alpes.fr
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 4 Jun 2026
! *
! ******************************************************************************/
!--------------------------------------------------------------------------------
!>  Module containing utilityes for the coupling with
!>  PDAF - Parallel Data Assimilation Framework -
!>  https://github.com/PDAF/PDAF
!>  Contains variables and codes for the parallel initilisation and finilisation
!>  Adapted from the online parallel model templates and test cases
!--------------------------------------------------------------------------------
MODULE PDAFUtils
  IMPLICIT NONE

  ! Variables for each model task
  INTEGER :: COMM_model=0               !< MPI communicator for model tasks
  INTEGER :: mype_model=0               !< PE rank in COMM_model
  INTEGER :: npes_model=1               !< Number of PEs in COMM_model

  ! Variables describing all processes involved in model integrations
  integer :: COMM_ens=0           !< Jont Communicator for entire ensemble
  integer :: mype_ens=0          !< Rank in COMM_ens
  integer :: npes_ens=1          !< Size of COMM_ens

  ! Variables describing the processes involved in the analysis step
  integer :: COMM_filter=0              !< MPI communicator processes in analysis step
  integer :: mype_filter=1              !< Process rank in COMM_da
  integer :: npes_filter=0              !< Number of processes in COMM_da

  ! Additional variables for use with PDAF
  INTEGER :: n_modeltasks=1             !< Number of parallel model tasks
  INTEGER :: task_id=1                  !< Index of my model task (1,...,n_modeltasks)

  ! Logical : TRUE if PDAF is used
  LOGICAL :: PDAF_INITIALISED=.FALSE.

  CONTAINS

  ! ** Initialise parallel communicators:
  ! **  fully parallel implementation case; i.e. number of model tasks is equal to the ensemble size
  ! **  https://pdaf.awi.de/trac/wiki/OnlineAdaptParallelization_PDAF3
  SUBROUTINE init_parallel_pdaf(screen, model_comm, model_comm_rank, model_comm_size)

    USE PDAF, ONLY: PDAF_parse, PDAF3_init_parallel
    IMPLICIT NONE

    ! *** Arguments ***
    INTEGER, INTENT(in)    :: screen           !< Whether screen information is shown

    ! Model parallelization variables 
    INTEGER, INTENT(inout) :: model_comm       !< Model MPI communicator for model tasks
    INTEGER, INTENT(inout) :: model_comm_size  !< Number of processes in model_comm
    INTEGER, INTENT(inout) :: model_comm_rank  !< Process rank in model_comm

    ! *** Local variables ***
    INTEGER :: dim_ens                         ! Ensemble size
    CHARACTER(len=32) :: handle                ! Handle for command line parser

    !** Parse ensemble size using PDAF_parse
    dim_ens=0
    handle = 'dim_ens'
    CALL PDAF_parse(handle, dim_ens)

    !** IF dim_ens = 0 do not start PDAF
    IF (dim_ens == 0) RETURN

    ! Set number of model tasks for fully-parallel mode
    n_modeltasks = dim_ens

    ! Store the parallelization variables provided by the model
    ! At this point, they describe all processes doing model integrations
    COMM_ens   = model_comm
    mype_ens   = model_comm_rank
    npes_ens   = model_comm_size

    ! Initialize ensemble parallelization
    CALL PDAF3_init_parallel(screen, 0, 1, dim_ens, n_modeltasks, &
       model_comm, model_comm_rank, model_comm_size, &
       COMM_filter, mype_filter, npes_filter, &
       task_id)

   ! Initialize parallelization variables for parallel_pdaf_mod
   COMM_model = model_comm
   mype_model = model_comm_rank
   npes_model = model_comm_size

   ! PDAF has been initilised and used
   PDAF_INITIALISED=.TRUE.

  END SUBROUTINE init_parallel_pdaf

  SUBROUTINE finalize_pdaf()
    USE PDAF, ONLY: PDAF_print_info, PDAF_finalize
    IMPLICIT NONE

  ! *** Show allocated memory for PDAF ***
    IF (mype_ens==0) CALL PDAF_print_info(10)

  ! *** Print PDAF timings onto screen ***
    IF (mype_ens==0) CALL PDAF_print_info(3)

  ! *** Deallocate PDAF arrays ***
    CALL PDAF_finalize()

  END SUBROUTINE finalize_pdaf

END MODULE PDAFUtils



