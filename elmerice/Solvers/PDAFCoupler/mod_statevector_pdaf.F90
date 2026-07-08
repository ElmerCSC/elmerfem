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
! * Provide variable and routine to build the state vector and communicate between
! * Elmer model tasks and the filter tasks
MODULE mod_statevector_pdaf
  USE DefUtils
  USE PDAFUtils

  IMPLICIT NONE
  SAVE

  !< Fortran type storing size and offset of each model field in the state vector
  TYPE state_field
     INTEGER :: dim               ! size of field in state vector
     INTEGER :: off               ! offset of field in state vector
     CHARACTER(len=MAX_NAME_LEN) :: name    ! Name of field variable
     INTEGER :: SolverId=-1       ! Solver Id for state variables
     INTEGER, ALLOCATABLE :: StateToNode(:) ! State to node permutation table
     INTEGER, ALLOCATABLE :: NodeToState(:) ! Node to position in this state 
     INTEGER :: mask=-1           ! mask variable 
  END TYPE state_field

  !< number of fields in state vector
  INTEGER :: n_fields                   
  !< Vector of type variable holding dimension and offset of each field
  TYPE(state_field), ALLOCATABLE :: sfields(:)

  !< number of "observed" fields; i.e. variable that are not "state" variables
  !< but used in the observations operators (e.g. ice velocity) 
  INTEGER :: n_ofields
  TYPE(state_field), ALLOCATABLE :: ofields(:)

  !< mask variables; e.g. groundedmask where ice is floating and thus bed properties 
  !< are not seen but the model and must be excluded of the analysys if they are part
  !< of the state vector
  INTEGER :: n_mfields
  TYPE(state_field), ALLOCATABLE :: mfields(:)

  !< Perumation table from local analysys domain to mesh nodes
  INTEGER, DIMENSION(:), ALLOCATABLE :: DomainToNode

  !< pe_local full "observed" ensemble
   REAL(KIND=dp),ALLOCATABLE :: ObsEns(:,:)
  !< pe_local ReducedMask (gathered with MPI_MAX)
   REAL(KIND=dp),ALLOCATABLE :: ReducedMask(:)

CONTAINS
! ***********************************************************************
!* Initialize state vector from model fields
!* - Executed by each process participating in the model integration 
! ***********************************************************************
SUBROUTINE collect_state_pdaf(dim_p, state_p)
  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(kind=dp), INTENT(inout) :: state_p(dim_p)  !< local state vector

! *** local variables ***
  CHARACTER(*), PARAMETER :: Caller = "PDAF_ELMER: collect_state_pdaf"
  INTEGER :: TimeStep

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************
  TimeStep = GetTimeStep()
  WRITE (message, '(A)') 'Current TimeStep : '//I2S(TimeStep)
  call Info(Caller,message,level=3)

  CALL collect_elmer_state(dim_p,state_p,n_fields,sfields)

END SUBROUTINE collect_state_pdaf

! ***********************************************************************
!* Transfer Elmer variable to the "state" variable
! ***********************************************************************
SUBROUTINE collect_elmer_state(dim_p,state_p,n_f,sf)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: dim_p  !pe local dim os the "state"
  REAL(kind=dp), INTENT(inout) :: state_p(dim_p) ! pe local "state"
  INTEGER, INTENT(in) :: n_f    ! number of fields in sf
  TYPE(state_field) :: sf(n_f)  ! a state_field variable 


  CHARACTER(*), PARAMETER :: Caller = "PDAF_ELMER: collect_elmer_state"
  TYPE(Mesh_t),POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: Perm(:)
  INTEGER :: node
  INTEGER :: n,nf

  Mesh => GetMesh()

  DO nf=1,n_f
     VarName=TRIM(sf(nf)%name)
     Var => VariableGet( Mesh % Variables,VarName,UnFoundFatal=.TRUE.)
     Perm => Var % Perm

     DO n=1,sf(nf)%dim
       node = sf(nf)%StateToNode(n)
       IF (Perm(node) == 0) THEN
            CALL FATAL(Caller,'Problem with permutation')    
       END IF
       state_p( sf(nf)%off + n ) = Var % Values ( Perm(node) )
     END DO
  END DO
END SUBROUTINE collect_elmer_state
        
! ***********************************************************************
! *  Distribute state vector to Elmer variables
! *   Do nothing if TimeStep == 0 as we rely on Elmer initialisation for
! *     the members and nothing is done in init_ens_pdaf.F90
! *  Rq. We reset DoneTime to 0 for all solvers as the analysys introduce a 
! *    discontinuity that might me problematic if using integration scheme with order > 1
! ***********************************************************************
SUBROUTINE distribute_state_pdaf(dim_p, state_p)
  USE mod_io
  IMPLICIT NONE
  
! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p           !< PE-local state dimension
  REAL(kind=dp), INTENT(inout) :: state_p(dim_p)  !< PE-local state vector

! *** local variables ***
  CHARACTER(*), PARAMETER :: Caller = "distribute_state_pdaf"
  TYPE(Mesh_t),POINTER :: Mesh         !< Current Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  TYPE(Variable_t), POINTER :: Var
  INTEGER, POINTER :: Perm(:)
  INTEGER :: TimeStep
  INTEGER :: node
  INTEGER :: nf,n
  INTEGER :: i

! *************************************************
! *** Initialize model fields from state vector ***
! *** for process-local model domain            ***
!**************************************************
  Mesh => GetMesh()

  TimeStep = GetTimeStep()
  WRITE (message, '(A)') 'Current TimeStep : '//I2S(TimeStep)
  call Info(Caller,message,level=3)

  !* initialisation do nothing 
  IF (TimeStep == 0)  THEN
     IF (mype_ens == 0) THEN
       CALL LOCAL_INFO(Caller,"TimeStep == 0; Initial ensemble not distributed to Elmer",Level=1)
     END IF
     RETURN
  END IF

  DO nf=1,n_fields
     VarName=TRIM(sfields(nf)%name)
     Var => VariableGet( Mesh % Variables,VarName,UnFoundFatal=.TRUE.)
     Perm => Var % Perm

      DO n=1,sfields(nf)%dim
       node = sfields(nf)%StateToNode(n)
       IF (Perm(node) == 0) THEN
            CALL FATAL(Caller,'Problem with permutation')
       END IF
        Var % Values ( Perm(node) ) = state_p( sfields(nf)%off + n )
     END DO

  END DO

  !* Reset DoneTime to 0
  DO i = 1,CurrentModel % NumberOfSolvers
     CurrentModel % Solvers(i) % DoneTime = 0
  END DO

END SUBROUTINE distribute_state_pdaf

! ***********************************************************************************
! * Initialise state_field variables
! ************************************************************************************
  SUBROUTINE init_sfields(sf,n_f,Prefix)
    IMPLICIT NONE

    TYPE(state_field), ALLOCATABLE :: sf(:)
    INTEGER :: n_f
    CHARACTER(LEN=*) :: Prefix

! *** Local Variables ***
    CHARACTER(*), PARAMETER :: Caller = "PDAF_ELMER: init_sfields"
    TYPE(Mesh_t),POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    TYPE(Variable_t), POINTER :: Var

    INTEGER :: mask
    INTEGER :: j,cnt,node
    LOGICAL :: Found

    Mesh => GetMesh()
    Solver => CurrentModel % Solver

    DO j=1,100
      VarName = ListGetString( Solver%Values ,TRIM(Prefix)//' '//I2S(j), Found )
      IF( .NOT. Found ) EXIT
    END DO
    n_f = j - 1
    IF ( n_f == 0 ) RETURN

    ALLOCATE(sf(n_f))

    DO j=1,n_f
      VarName = ListGetString( Solver%Values ,TRIM(Prefix)//' '//I2S(j), UnFoundFatal=.TRUE.)

      Var => VariableGet( Mesh % Variables,VarName,UnFoundFatal=.TRUE.)
      IF (.NOT.ASSOCIATED(Var % Perm)) CALL FATAL(Caller,'Var % Perm not associated')
      IF (Var % Dofs > 1) CALL FATAL(Caller,'Var % Dofs > 1 not supported')
      IF (Var % Type /= Variable_on_nodes) CALL FATAL(Caller,'Var % Type not supported')

      IF (ASSOCIATED(Var % Solver)) THEN 
        !This is the main solver variable
        IF (ASSOCIATED(Var,Var % Solver % Variable)) &
                sf(j)%SolverId=Var % Solver % SolverId
      END IF

      ! name of the variable
      sf(j)%name = TRIM(VarName)
      ! dimension
      sf(j)%dim=COUNT(Var % Perm(:) > 0)

      mask = ListGetInteger( Solver%Values ,TRIM(Prefix)//' '//I2S(j)//' use mask', Found)
      IF (Found) sf(j)%mask=mask

      ALLOCATE(sf(j)%StateToNode(sfields(j)%dim))
      ALLOCATE(sf(j)%NodeToState(Mesh%NumberOfNodes))
      sf(j)%StateToNode = -1
      sf(j)%NodeToState = -1

      cnt=0
      DO node=1,Mesh%NumberOfNodes
       IF (Var % Perm(node) == 0) CYCLE
       cnt=cnt+1
       sf(j)%StateToNode(cnt)=node
       sf(j)%NodeToState(node)=cnt
     END DO
     IF (cnt /= sf(j)%dim) &
        CALL FATAL(Caller,'Error in field dimension for state variable '//I2S(j))
      
    END DO

! **************************************
! ***   Set offsets                  ***
! **************************************
    ! Define field offsets in state vector
    sf(1)%off = 0
    DO j = 2, n_f
       sf(j)%off = sf(j-1)%off + sf(j-1)%dim
    END DO

  END SUBROUTINE init_sfields

! ***********************************************************************************
! * Initialise the state vector
! ************************************************************************************
  SUBROUTINE setup_statevector(dim_state, dim_state_p, screen)
    USE PDAFUtils, &
         ONLY: COMM_ens,COMM_model,mype_model, npes_model, task_id

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(out) :: dim_state    !< Global dimension of state vector
    INTEGER, INTENT(out) :: dim_state_p  !< Local dimension of state vector
    INTEGER, INTENT(in)  :: screen       !< Verbosity flag

! *** Local variables ***
    INTEGER :: i                         ! Counters
    INTEGER :: MPIerr                    ! Error flag for MPI
    TYPE(Mesh_t),POINTER :: Mesh


! ***********************************
! *** Initialize the state vector ***
! ***********************************
    CALL init_sfields(sfields,n_fields,'State Variable')

! *** Set state vector dimension ***
    dim_state_p = SUM(sfields(:)%dim)

! *****************************************
! *** Initialize the observed variables ***
! *****************************************
    CALL init_sfields(ofields,n_ofields,'Observed Variable')

! *****************************************
! *** Initialize the observed variables ***
! *****************************************
    CALL init_sfields(mfields,n_mfields,'Mask Variable')

! **************************************
! ***   Allocate DomainToNode        ***
! **************************************
    Mesh => GetMesh()
    ALLOCATE(DomainToNode(Mesh%NumberOfNodes))

! *** Get global state dimension ***
    CALL MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)

! *** Write information about the state vector ***

    IF (task_id==1) THEN
       IF (mype_model==0) THEN
          WRITE (*,'(/a,2x,a)') 'ELMER-PDAF','+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'     
          WRITE (*,'(a,2x,a)') 'ELMER-PDAF', '*** Setup of state vector ***'
          WRITE (*,'(a,3x,a,i5)') 'ELMER-PDAF', '--- Number of fields in state vector:', n_fields
          WRITE (*,'(a,a7,3x,a2,4x,a8,5x,a9,6x,a6,6x,a8)') &
               'ELMER-PDAF','proc.','ID', 'variable', 'dimension', 'offset', 'SolverId'
       END IF

       IF ((mype_model==0 .AND. screen<=2) .OR. screen>2) THEN
          DO i = 1, n_fields
             WRITE (*,'(a, i6,2x,i4,4x,a10,2x,i10,2x,i10,i10)') &
                  'ELMER-PDAF', mype_model, i, TRIM(sfields(i)%name), sfields(i)%dim, sfields(i)%off,sfields(i)%SolverId
          END DO
       END IF

       IF (npes_model>1) THEN
          IF (screen>2 .OR. mype_model==0) WRITE (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'ELMER-PDAF', 'PE', mype_model, 'process-local full state dimension: ',dim_state_p
       END IF
       IF (mype_model==0) &
            WRITE (*,'(a,2x,a,1x,i10)') 'ELMER-PDAF', 'Global state dimension: ',dim_state

       IF (n_ofields > 0 ) THEN
        IF (mype_model==0) THEN       
          WRITE (* , *)      
          WRITE (*,'(a,3x,a,i5)') 'ELMER-PDAF', '--- Number of fields in observed vector:', n_ofields
            WRITE (*,'(a,a7,3x,a2,4x,a8,5x,a9,6x,a6,6x,a8)') &
                 'ELMER-PDAF','proc.','ID', 'variable', 'dimension', 'offset', 'SolverId'
         END IF
         IF ((mype_model==0 .AND. screen<=2) .OR. screen>2) THEN
          DO i = 1, n_ofields
             WRITE (*,'(a, i6,2x,i4,4x,a10,2x,i10,2x,i10,i10)') &
                  'ELMER-PDAF', mype_model, i, TRIM(ofields(i)%name), ofields(i)%dim, ofields(i)%off,ofields(i)%SolverId
          END DO
       END IF
       END IF

       IF (mype_model==0) &
          WRITE (*,'(a,2x,a)') 'ELMER-PDAF','+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'     
    END IF
    call MPI_Barrier(comm_ens, MPIerr)

  END SUBROUTINE setup_statevector

! ***********************************************************************************
! * Collect mask variable it filter_pes accross the ensemble using MPI_MAX
! ************************************************************************************
  SUBROUTINE COLLECT_MASK(mask_p,dim_m)
  USE PDAFUtils
  ! get the Couple communicator from th emain PDAF parallel module as this is not exposed 
  !  to Elmer otherwise
  USE PDAF_mod_parallel, ONLY : filterpe,modelpe,npes_couple,mype_couple,COMM_couple
  ! *** Arguments ***
  INTEGER, INTENT(in) :: dim_m                   !< PE-local mask dimension
  REAL(kind=dp), INTENT(inout) :: mask_p(dim_m)  !< local mask vector  
  LOGICAL,SAVE :: FirstTime=.TRUE.
  CHARACTER(*), PARAMETER :: Caller="COLLECT_MASK"
   INTEGER :: MPIerr

  IF (FirstTime) THEN
    ALLOCATE(ReducedMask(dim_m))
    FirstTime=.FALSE.
  ENDIF

  CALL MPI_REDUCE(mask_p,ReducedMask,dim_m,MPI_DOUBLE,MPI_MAX,0,COMM_couple, MPIerr)

  END SUBROUTINE COLLECT_MASK        

! ***************************************************************************************************
! * Collect te full observed ensemble in filters_pes
! *  Rq. this is specific for the full parallelisation where each model task integrate only 1 member
! ***************************************************************************************************
  SUBROUTINE COLLECT_ObsEns(state_p,dim_p)
  USE PDAFUtils        
  USE PDAF, ONLY : PDAF_get_memberid
  USE PDAF_mod_parallel, ONLY : filterpe,modelpe,npes_couple,mype_couple,COMM_couple
  USE mod_assimilation, ONLY: dim_ens

  ! *** Arguments ***
  INTEGER, INTENT(in) :: dim_p                    !< PE-local state dimension
  REAL(kind=dp), INTENT(inout) :: state_p(dim_p)  !< local state vector
  
  LOGICAL,SAVE :: FirstTime=.TRUE.
   CHARACTER(*), PARAMETER :: Caller="TEST_COLLECT_STATE"
   INTEGER :: memberid
   INTEGER :: MPIerr,MPIstatus(MPI_STATUS_SIZE)
   INTEGER :: pe_rank
   INTEGER :: col_frst

  !Rq. will always be 1, if only 1 member by model task => add consistency check??
  ! otherwise get memberid from model task.....
  CALL PDAF_get_memberid(memberid)

  IF (FirstTime) THEN
    IF (.NOT.filterpe) THEN
            !* non filter pes send their local state
            ALLOCATE(ObsEns(dim_p,1))
    ELSE
            !* filterpe collect the whole ensemble
            ALLOCATE(ObsEns(dim_p,dim_ens))
    END IF
    FirstTime=.FALSE.        
  ENDIF
  ObsEns(:,memberid)=state_p(:)

  IF (.NOT.filterpe .AND. npes_couple > 1) THEN
        ! Send sub-ensembles to couple PEs with rank 0
       CALL MPI_SEND(ObsEns, dim_p, MPI_DOUBLE_PRECISION, 0, mype_couple, &
            COMM_couple, MPIerr)
  ELSE
      DO pe_rank = 1, npes_couple - 1    
        col_frst = pe_rank + 1
        call MPI_recv(ObsEns(1, col_frst), &
                  dim_p , MPI_DOUBLE_PRECISION, &
                  pe_rank, pe_rank, COMM_couple, MPIstatus, MPIerr)
      END DO
  ENDIF
  END SUBROUTINE COLLECT_ObsEns

END MODULE mod_statevector_pdaf
