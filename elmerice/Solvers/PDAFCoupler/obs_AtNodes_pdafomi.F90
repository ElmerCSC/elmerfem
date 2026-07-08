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
MODULE obs_AtNodes_pdafomi
  USE DefUtils
  USE PDAFUtils, &
       ONLY: mype_filter,npes_filter    ! Rank of filter process
  USE PDAF, &
       ONLY: obs_f, obs_l,PDAF_DA_GENOBS   ! Declaration of observation data types
  USE mod_observations
  USE mod_io, ONLY: nfcheck
  USE netcdf
  IMPLICIT NONE

CONTAINS
!!
  SUBROUTINE init_dim_obs_AtNodes(step,dim_obs,thisobs,obsinfo,sfields,firsttime)

    USE PDAF, &                ! Include PDAF and PDAF-OMI routines
         ONLY: PDAF_local_type, PDAFomi_gather_obs, PDAFomi_set_localize_covar,PDAF_abort

    USE mod_assimilation, &    
         ONLY: filtertype, cradius, sradius, locweight,twin_experiment !, coords_p

    USE mod_observations, ONLY : delt_obs

    USE mod_synobs, &               ! Routines for synthetic observations
         ONLY: init_file_syn_obs, read_syn_obs
    USE mod_io
    USE mod_statevector_pdaf , ONLY : state_field

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector
    TYPE(obs_f), TARGET :: thisobs
    TYPE(obsinfo_t),INTENT(inout) :: obsinfo
    TYPE(state_field), INTENT(in) :: sfields(:) 
    LOGICAL,INTENT(in) :: firsttime

! *** Local variables ***
     CHARACTER(*), PARAMETER :: Caller = "init_dim_obs_AtNodes"
    CHARACTER(LEN=MAX_PATH_LEN) :: file_synobs
    INTEGER,SAVE :: iter=0
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL(kind=dp), ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL(kind=dp), ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL(kind=dp), ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    TYPE(Mesh_t),POINTER :: Mesh
    INTEGER :: node,NIndex
    REAL(kind=dp),parameter :: reps=1.0d-6
    INTEGER :: cnt
    INTEGER :: sid,sn
    LOGICAL :: Parallel
    ! netcdf variables
    INTEGER :: ncid
    INTEGER :: VarId
    INTEGER :: ndims
    INTEGER :: dimids(2)
    INTEGER :: nvals,ntime
    REAL(KIND=dp), ALLOCATABLE :: Values(:)
    REAL(KIND=dp) :: fillv
    INTEGER :: TimeIndex


    Parallel = (npes_filter > 1)

    Mesh => GetMesh()
! *********************************************
! *** Initialize full observation dimension ***
! *********************************************
    IF (mype_filter==0) THEN
       WRITE (Message,'(A)') 'Assimilate observations '//TRIM(obsinfo % name)
       CALL LOCAL_INFO(Caller,Message,Level=1)
    END IF

    ! Store whether to assimilate this observation type (used in routines below)
    thisobs%doassim = 1

    ! Specify whether to (1) use global observations for local filters,
    ! or (0) restrict the full observations to those relevant for a process domain
    thisobs%use_global_obs = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    thisobs%ncoord = obsinfo % dim
    IF ((obsinfo % dim < 1).OR.(obsinfo % dim > 3)) &
       CALL PDAF_FATAL(Caller,'Wrong dimension for observation '//TRIM(obsinfo % name))


    ! Open file
    call nfcheck( NF90_OPEN(trim(obsinfo % input_file), NF90_NOWRITE,ncid))
    call nfcheck( NF90_INQ_VARID(ncid,trim(obsinfo % name),VarId))
    call nfcheck( nf90_get_att(ncid, VarId,"_FillValue",fillv))
    call nfcheck( NF90_inquire_variable(ncid, VarId, ndims=ndims) )

    call nfcheck( nf90_inquire_variable(ncid, VarId, dimids = dimids(1:ndims)))

    ! here we assume that the first dim is the size of the variable
    ! and that the second should be time
    call nfcheck( nf90_inquire_dimension(ncid,dimids(1),len=nvals))
    ALLOCATE(Values(nvals))

    SELECT CASE(ndims)

      CASE(1)
        call nfcheck( nf90_get_var(ncid,VarId,Values))      

      CASE(2)        

        call nfcheck( nf90_inquire_dimension(ncid,dimids(2),len=ntime))
        TimeIndex = step / delt_obs
        IF ((TimeIndex < 1).OR.(TimeIndex > ntime)) THEN
          WRITE (*,*) 'init_dim_obs_AtNodes: Problem with Time Index '//I2S(TimeIndex)
          CALL PDAF_abort(1)
        END IF
        call nfcheck( nf90_get_var(ncid,VarId,Values,start=(/1,TimeIndex/),count=(/nvals/)))

      CASE DEFAULT

    END SELECT

    call nfcheck( NF90_CLOSE(ncid))

    !!!! TO SEE NOT OPTIMAL IF WORKING ON A BOUNDARY OF A WHOLE MESH
    !!!!  AS THE INPUT FILE MUST PROVIDE ALL MESH NODES....
    IF (Parallel) THEN
      ! count obs in domain
      cnt=0
      DO node=1,Mesh%NumberOfNodes
        ! node not owned by current partition
        IF (Mesh % ParallelInfo % NeighbourList(node) % Neighbours(1).NE.ParEnv % MyPE) CYCLE
        NIndex=Mesh % ParallelInfo % GlobalDOFs(node)     
        IF (NIndex > nvals) CALL PDAF_FATAL(Caller,'NIndex larger than nvals '//I2S(NIndex))
        IF (Values(NIndex) /= fillv) cnt=cnt+1        
      END DO
    ELSE 
      cnt = COUNT(Values /= fillv)
    END IF

    dim_obs_p = cnt
    obsinfo % dim_obs_p = dim_obs_p
       
      
    haveobs: if (dim_obs_p > 0) then
! +++++ Its only purpose is to let the program run with segmentation fault

    ! Allocate process-local observation arrays
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(thisobs%ncoord, dim_obs_p))

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    ALLOCATE(thisobs%id_obs_p(1 , dim_obs_p))

   cnt=0
   DO node=1,Mesh%NumberOfNodes

     IF (Parallel) THEN
        ! node not owned by current partition
        IF (Mesh % ParallelInfo % NeighbourList(node) % Neighbours(1).NE.ParEnv % MyPE) CYCLE
        NIndex=Mesh % ParallelInfo % GlobalDOFs(node)     
     ELSE
        NIndex=node     
     END IF

     IF (NIndex > nvals) CALL PDAF_FATAL(Caller,'NIndex larger than nvals '//I2S(NIndex))

     IF (Values(NIndex) /= fillv) THEN
        cnt=cnt+1        

        obs_p(cnt) = Values(NIndex)

        ! corresponding state id
        sid = obsinfo % sid
        if (sid < 1) CALL PDAF_FATAL(Caller,&
                'Observation '//TRIM(obsinfo % name)//' has no valid corresponding state')

        sn=sfields(sid)%NodeToState(node)
        if (sn < 1) CALL PDAF_FATAL(Caller,&
                'Observation '//TRIM(obsinfo % name)//' problem in NodeToState at node '//I2S(node))
        thisobs%id_obs_p(1,cnt) = sfields(sid)%off + sn

        ocoord_p(1,cnt)=Mesh%Nodes%x(node)
        IF (thisobs%ncoord > 1) &
            ocoord_p(2,cnt)=Mesh%Nodes%y(node)
        IF (thisobs%ncoord == 3) &
            ocoord_p(2,cnt)=Mesh%Nodes%z(node)    
     END IF
   END DO
   IF (cnt /= dim_obs_p) THEN
      WRITE (*,*) 'init_dim_obs_AtNodes: pb with obs size'
      CALL PDAF_abort(1)
   END IF

! **********************************
! *** Read PE-local observations ***
! **********************************
   IF ((abs(obsinfo%rms)-DEFAULT_RMS) > reps) THEN
      ivar_obs_p(:) =  1.0 / (obsinfo%rms * obsinfo%rms)
   ELSE
     WRITE (*,*) 'init_dim_obs_AtNodes: non uniform rms not supported yet...'
     CALL PDAF_abort(1)
   ENDIF



   else haveobs

       ! *** For dim_obs_p=0 allocate arrays with minimum size

       allocate(obs_p(1))
       allocate(ivar_obs_p(1))
       allocate(ocoord_p(2, 1))
       allocate(thisobs%id_obs_p(1, 1))

    end if haveobs

! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *** and overwrite real observation values.            ***
! *********************************************************
    IF (twin_experiment .AND. filtertype/=PDAF_DA_GENOBS) THEN
       ! Set file name (separate files for each process)
       file_synobs = 'syn_obs_'//TRIM(obsinfo%name)//'_'//I2S(mype_filter)//'.nc'

       ! Read synthetic observations
       CALL read_syn_obs(file_synobs, dim_obs_p, obs_p, step)
    END IF

! ****************************************
! *** Gather global observation arrays ***
! ****************************************
    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, cradius, dim_obs)

! ***************************************************************
! *** Provide localization information for LEnKF, EnsRF, EAKF ***
!    TO SEE
! ***************************************************************
!    IF (PDAF_local_type() > 1) THEN
!       CALL PDAFomi_set_localize_covar(thisobs, dim_state_p, n_dim, coords_p, &
!            locweight, cradius, sradius)
!    END IF

! ***************************************************************
! *** When generating synthetic observations: Initialize file ***
! ***************************************************************
    if (filtertype==PDAF_DA_GENOBS .and. firsttime) then
       ! Set file name (separate files for each process)
       file_synobs = 'syn_obs_'//TRIM(obsinfo%name)//'_'//I2S(mype_filter)//'.nc'

       ! Initialize file
       call init_file_syn_obs(dim_obs_p, file_synobs)
    end if

! ********************
! *** Finishing up ***
! ********************
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)
    DEALLOCATE(Values)

  END SUBROUTINE init_dim_obs_AtNodes



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_state_AtNodes(dim_p, dim_obs, state_p, ostate,thisobs)
    USE PDAF, &                ! Include PDAF-OMI routine
         ONLY: PDAFomi_obs_op_gridpoint
    USE PDAFUtils

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL(kind=dp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL(kind=dp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state
    TYPE(obs_f), TARGET :: thisobs

    CHARACTER(*), PARAMETER :: Caller = "obs_op_state_AtNodes"

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************
    ! Example: Observation operator for observed grid point values
    CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    !WRITE (*,*) TRIM(Caller),' 1 ',task_id,thisobs%id_obs_p(1,:)
    !WRITE (*,*) TRIM(Caller),' 2 ',task_id,ostate

  END SUBROUTINE obs_op_state_AtNodes

   SUBROUTINE obs_op_obs_AtNodes(dim_p, dim_obs, state_p, ostate,thisobs)
    USE PDAFUtils
    USE PDAF
    USE mod_statevector_pdaf, ONLY : ObsEns
    USE mod_io

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL(kind=dp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL(kind=dp), INTENT(inout) :: ostate(dim_obs)       !< Full observed state
    TYPE(obs_f) :: thisobs

    ! *** Local variables ***
    CHARACTER(*), PARAMETER :: Caller = "obs_op_obs_AtNodes"
    INTEGER :: i                           ! Counter
    REAL(kind=dp), ALLOCATABLE :: ostate_p(:)       ! local observed part of state vector
    INTEGER :: obsmemberid
    !REAL(KIND=dp) :: DIFF(dim_p)


    doassim: IF (thisobs%doassim == 1) THEN
       ! Consistency check
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: PDAFomi_obs_op_gridpoint - thisobs%id_obs_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       CALL PDAF_get_obsmemberid(obsmemberid)
       IF (obsmemberid == 0) &
           CALL PDAF_FATAL(Caller,'obsmemberid == 0 detected ')

       ! Rq. pour les filtresglobaux, c'est appelé au debut et apres l'analyse si omi_stat=True
       !WRITE(Message,'(A)') " : task_id="//I2S(task_id)//" obsmid="//I2S(obsmemberid)
       !call Info(Caller,message,level=3)

       !DIFF=state_p(:)-ObsEns(:,obsmemberid)
       !WRITE(*,*) TRIM(Message),SIZE(ObsEns,1),SIZE(ObsEns,2),SUM(DIFF)
       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = ObsEns(thisobs%id_obs_p(1, i),obsmemberid)
       ENDDO

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)

    !WRITE (*,*) TRIM(Caller),' 1 ',task_id,thisobs%id_obs_p(1,:)
    !WRITE (*,*) TRIM(Caller),' 2 ',task_id,ostate
       ! *** Clean up
       DEALLOCATE(ostate_p)

    END IF doassim

  END SUBROUTINE obs_op_obs_AtNodes
!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_AtNodes(domain_p, step, dim_obs, dim_obs_l,thisobs,thisobs_l)

    USE PDAF, &                ! Include PDAF-OMI routine
         ONLY: PDAFomi_init_dim_obs_l, PDAF_abort

    ! Include localization radius and local coordinates
    ! one can also set observation-specific values for the localization.
    USE mod_assimilation, &   
         ONLY: coords_l, cradius, locweight, sradius

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector
    TYPE(obs_f), TARGET :: thisobs
    TYPE(obs_l), TARGET :: thisobs_l

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************
    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.

    ! For cradius and sradius:
    ! If these are defined as scalar values, isotropic localization is used.
    ! If these are vectors, nonisotropic localization is used
    !   (their length has to be equal to thisobs%ncoord)

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, cradius, sradius, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_AtNodes

END MODULE obs_AtNodes_pdafomi
