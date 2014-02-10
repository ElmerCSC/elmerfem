!------------------------------------------------------------------------------
! Vili Forsell
! Created: 13.6.2011
! Last Modified: 4.8.2011
!------------------------------------------------------------------------------
! This module contains functions for
! - getting dimensions sizes and NetCDF identifiers; GetAllDimensions()
! - getting data from NetCDF files; GetFromNetCDF()
! - handling NetCDF status errors; G_Error()
!------------------------------------------------------------------------------
MODULE NetCDFGeneralUtils
  USE DefUtils, ONLY: dp, MAX_NAME_LEN
  USE NetCDF
  USE Messages
  IMPLICIT NONE
  LOGICAL, PARAMETER :: DEBUG_UTILS = .FALSE.

  !--- A type for time dimension values
  TYPE TimeType_t
    LOGICAL :: is_defined
    REAL(KIND=dp) :: val ! Exact time value
    INTEGER :: id, & ! Dimension NetCDF id
              len, & ! Dimension length/size
              low, & ! Nearest lower index
              high ! Nearest higher index
    LOGICAL :: doInterpolation ! True, if the exact value is not an integer and, hence, requires interpolation
  END TYPE TimeType_t

  CONTAINS
  
    !------------------ CloseNetCDF() -----------------------
    !--- Closes the given NetCDF file
    SUBROUTINE CloseNetCDF( NCID )
      USE NetCDF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCID
      INTEGER :: status
      
      status = NF90_CLOSE(NCID)
      IF ( G_ERROR(status,'Failed to close NetCDF file.') ) THEN
        CALL abort()
      END IF
    END SUBROUTINE CloseNetCDF
  
    !------------------ GetDimension() ----------------------
    !--- Takes the NetCDF file identifier and dimension name (as in NetCDF) and gets the id and length of the dimension, or abort otherwise
    SUBROUTINE GetDimension( NCID, DIM_NAME, dim_id, dim_len )
    !--------------------------------------------------------
      USE Messages
      IMPLICIT NONE
      !--- Arguments
      INTEGER, INTENT(IN) :: NCID
      CHARACTER (*), INTENT(IN) :: DIM_NAME
      INTEGER, INTENT(OUT) :: dim_id, dim_len 
      
      !--- Variables
      CHARACTER :: tmp_name ! Temporary name
      INTEGER :: status ! Results and status information from NetCDF
      INTEGER, PARAMETER :: sentinel = -1 ! Default intial value in case of error

      !--- Initializations      
      dim_id = sentinel
      dim_len = sentinel
      
      !--- Get dimension information and check success
      status = NF90_INQ_DIMID(NCID,DIM_NAME,dim_id)
      IF ( .NOT. G_Error(status, 'Dimension identifier could not be found.') ) THEN
        status = NF90_INQUIRE_DIMENSION(NCID,dim_id,tmp_name,dim_len)
        IF ( G_Error(status, 'Dimension could not be inquired.') ) THEN
          dim_id = sentinel
          dim_len = sentinel
          CALL abort()
        END IF
      ELSE
        dim_id = sentinel
        CALL abort()
      END IF
      status = NF90_INQ_VARID(NCID,DIM_NAME,dim_id)
      WRITE(Message,'(A,A,A)') 'No variable corresponding to the dimension ', DIM_NAME ,' could be found.'
      IF ( G_Error(status, Message) ) THEN
        dim_id = sentinel
        CALL abort()
      END IF
    
    IF ( DEBUG_UTILS ) THEN ! Debug printouts
      WRITE (Message,'(A,A10,A,I5,A,I10,A)') 'Dimension: ', DIM_NAME,' with id ', dim_id, ' and size ', dim_len, ' read correctly'
      CALL Info('GridDataMapper',Message)
    END IF
    END SUBROUTINE GetDimension
   
  
    !------------------ GetAllDimensions() ----------------------
    !--- Takes the NetCDF file and name identifiers (as in NetCDF) and gets the ids and lengths of the dimensions, or abort otherwise
    SUBROUTINE GetAllDimensions( NCID, NAMES, dim_ids, dim_lens )
    !------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER (len = MAX_NAME_LEN), INTENT(IN) :: NAMES(:)
      INTEGER, INTENT(IN) :: NCID
      INTEGER, ALLOCATABLE, INTENT(OUT) :: dim_ids(:),dim_lens(:)
      INTEGER :: alloc_stat, nm

      IF ( (size(NAMES,1) .NE. size(dim_ids)) .AND. (size(dim_ids) .NE. size(dim_lens)) ) THEN
        CALL Fatal('GridDataMapper','GetAllDimensions() input dimensions do not agree!')
      END IF

      ! Allocates the result vectors
      ALLOCATE ( dim_ids(size(NAMES,1)), dim_lens(size(NAMES,1)), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF

      ! Collects the data for each name in order
      DO nm = 1,size(NAMES,1),1
        CALL GetDimension( NCID,NAMES(nm),dim_ids(nm),dim_lens(nm) )
      END DO
!      WRITE(*,*) 'End'
    END SUBROUTINE GetAllDimensions
  
   
    !----------------- GetFromNetCDF() --------------------
    !--- Reads the given variable name and returns the value from NetCDF grid
    !--- Does not require time, loc contains indexing
    FUNCTION GetFromNetCDF( NCID, VAR_ID, LOC, LOC_TIME, TIME, DIM_LENS, accessed, OUT_SIZE ) RESULT( success )
    !------------------------------------------------------
      USE NetCDF
      IMPLICIT NONE

      !------------------------------------------------------------------------------
      ! ARGUMENTS
      !------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: DIM_LENS(:), VAR_ID
      INTEGER, INTENT(IN) :: NCID, LOC(:), LOC_TIME
      TYPE(TimeType_t), INTENT(IN) :: TIME
      INTEGER, INTENT(IN) ::  OUT_SIZE(:) ! The sizes of each dimension for the return type
      LOGICAL :: success ! Output: TRUE if all ok
      REAL (KIND=dp), ALLOCATABLE, INTENT(INOUT) :: accessed(:) ! Later reshaped to proper dimensions

      !------------------------------------------------------------------------------
      ! VARIABLES
      !------------------------------------------------------------------------------
      INTEGER :: DIM_COUNT,TOTAL_SIZE
      INTEGER :: alloc_stat ! alloc_stat for allocation status
      INTEGER, ALLOCATABLE :: COUNT_VECTOR(:)
      ! COUNT_VECTOR is the amount of nodes taken
      ! starting from corresponding index vector locations (slabs of data)
      INTEGER, ALLOCATABLE :: locs(:,:) ! First column is left limit, second column is right limit
      INTEGER, ALLOCATABLE :: rev_perm(:) ! Reverse permuation for indices; necessary to compensate for the
      ! NetCDF Fortran access indexing being opposite to the one shown in network Common Data form Language (CDL)
      INTEGER :: loop, status,timeBias ! timeBias is the change to array indices caused by taking time dimension into account
      CHARACTER(len=64) :: tmpFormat
      
      !------------------------------------------------------------------------------
      ! INITIALIZATIONS
      !------------------------------------------------------------------------------

      ! Checks input size
      IF ( (size(OUT_SIZE) .NE. size(LOC)) .OR. (size(LOC) .NE. size(DIM_LENS)) ) THEN
        WRITE(Message, '(A,I3,A,I3,A,I3)') 'Number of dimensions differs between input coordinates: Output size ',&
                           size(OUT_SIZE), ', # of locations ', size(LOC), ', # of dimensions ', size(DIM_LENS)
        CALL Fatal('GridDataMapper',Message)
      END IF

      ! Checks if time dimension is taken into account (picks always just one point)
      timeBias = 0
      IF (TIME % IS_DEFINED) timeBias = 1

      DIM_COUNT = size(DIM_LENS,1) + timeBias ! May take one more dimension for time

      ! The same calculations would be done in any case; uses a little memory to save here
      ! The product of all sizes is the size of the one dimensional version of the NetCDF return value
      TOTAL_SIZE = 1
      DO loop = 1,size(OUT_SIZE,1),1
        TOTAL_SIZE = TOTAL_SIZE*OUT_SIZE(loop)
      END DO

      success = .FALSE. ! For checking if all went ok (allows later error recuperation)
      
      ALLOCATE ( accessed(TOTAL_SIZE), COUNT_VECTOR(DIM_COUNT),locs(DIM_COUNT,2),rev_perm(DIM_COUNT), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF
  
      accessed = 0

      ! In network Common Data form Language (CDL) notation infinite dimensions are first (f.ex. time)
      ! However, with Fortran access functions the order is required to be reversed, so it will be last
      rev_perm = (/ (DIM_COUNT - loop + 1, loop = 1, DIM_COUNT) /)

      ! If has time, then the first dimension is time and is set in count vector and locs
      ! When using via the mask rev_perm, the time index will then be last
      IF ( TIME % IS_DEFINED ) THEN
        COUNT_VECTOR(1) = 1
        locs(1,1) = LOC_TIME
      END IF
      COUNT_VECTOR(1+timeBias:size(OUT_SIZE)+timeBias) = OUT_SIZE(:)
      locs(1+timeBias:size(LOC)+timeBias,1) = LOC(:)
      
      locs(:,2) = locs(:,1) + COUNT_VECTOR(:) - 1 ! Covers the stencil area starting from left indices
      
      ! Checks each dimension range (and, hence, access attempt)
      IF ( TIME % IS_DEFINED ) THEN
        IF ( (locs(1,1) .LT. 1) .OR. (TIME % LEN .LT. locs(1,2)) ) THEN
          WRITE(tmpFormat,'(A,I3,A,I3,A)') '(A,/,', size(locs,2),'(I10),/,',size(locs,2),'(I10))'
          WRITE (*,tmpFormat) 'Locs: ', TRANSPOSE(locs)
          WRITE(tmpFormat,'(A,I3,A,I3,A)') '(A,/,',size(DIM_LENS),'(I10))'
          WRITE (*,tmpFormat) 'Dims: ', DIM_LENS
          CALL Fatal('GridDataMapper','Indexing time out of bounds.')
        END IF
      END IF

      DO loop = 1+timeBias,size(locs,1),1
        IF ( (locs(loop,1) .LT. 1) .OR. (DIM_LENS(loop-timeBias) .LT. locs(loop,2)) ) THEN
          WRITE(tmpFormat,'(A,I3,A,I3,A)') '(A,/,', size(locs,2),'(I10),/,',size(locs,2),'(I10))'
          WRITE (*,tmpFormat) 'Locs: ', TRANSPOSE(locs)
          WRITE(tmpFormat,'(A,I3,A,I3,A)') '(A,/,',size(DIM_LENS),'(I10))'
          WRITE (*,tmpFormat) 'Dims: ', DIM_LENS
          CALL Fatal('GridDataMapper','Indexing parameter(s) out of bounds.')
        END IF
      END DO

      !--- The dimensions and the locations have been read and checked; NetCDF accessing info is a-ok
 
      ! Access variable and take the values
      status = NF90_GET_VAR(NCID,var_id,accessed,locs(rev_perm(:),1),COUNT_VECTOR(rev_perm(:)))
      WRITE(Message,*) 'NetCDF variable access failed. Variable ID: ',var_id,' and taking ',&
         size(accessed,1), ' elements with true size ', TOTAL_SIZE,'Locs: ',locs,'Dims: ',DIM_LENS
      IF ( G_ERROR(status,Message) ) THEN
        accessed = 0
        CALL abort()
      END IF
      
      success = .TRUE. ! Successful
  
    END FUNCTION GetFromNetCDF
   


    !----------------- TimeValueToIndex() ---------------
    !--- Takes a NetCDF time value and converts it into an index
    FUNCTION TimeValueToIndex(NCID,TIME_NAME,DIM_ID,DIM_LEN,t_val,t_eps) RESULT(t_ind)
      IMPLICIT NONE
  
      !--- Arguments
      CHARACTER(len = MAX_NAME_LEN), INTENT(IN) :: TIME_NAME
      REAL(KIND=dp), INTENT(IN) :: t_val
      REAL(KIND=dp), INTENT(IN) :: t_eps ! The rounding for time
      INTEGER, INTENT(IN) :: NCID, DIM_ID, DIM_LEN
      REAL(KIND=dp) :: t_ind ! Output
  
       !--- Variables
      REAL(KIND=dp) :: t_min, t_max, t_tmp1(1), t_tmp2(2), t_diff
      INTEGER :: time_id, status
      INTEGER :: index_scalar(1), count_scalar(1)
  
      t_ind = -1.0_dp ! Initialization to out of bounds
      index_scalar = 1 ! Initialized to min value
      count_scalar = 2
      time_id = DIM_ID ! Last dimension is time
    
      ! 1) Inquire time variable's id
      status = NF90_INQ_VARID(NCID,TIME_NAME,time_id)
      IF ( G_Error(status,'NetCDF time variable name not found.') ) THEN
        RETURN
      END IF
  
      ! 2) Get the time range from NetCDF
      status = NF90_GET_VAR(NCID,time_id,t_tmp2,index_scalar,count_scalar)
      IF ( G_Error(status,'First NetCDF time value not found') ) THEN
        RETURN
      END IF
      t_min = t_tmp2(1)
      t_diff = t_tmp2(2) - t_tmp2(1)
  
      count_scalar = 1
      index_scalar = DIM_LEN ! Pick the last max value
   
      status = NF90_GET_VAR(NCID,time_id,t_tmp1,index_scalar,count_scalar)
      IF ( G_Error(status,'Last NetCDF time value not found') ) THEN
        RETURN
      END IF
      t_max = t_tmp1(1)

      ! 3) Use the time range to find the index for the time value (NetCDF variables uniform)
      IF ( t_val < t_min .OR. t_val > t_max ) THEN
        ! Sets to the nearest value if within tolerance, else error, to better compare real values
        IF ( t_val < t_min .AND. t_val + (t_eps*t_diff) >= t_min ) THEN
          t_ind = 1
        ELSE IF ( t_val > t_max .AND. t_val - (t_eps*t_diff) <= t_max ) THEN
          t_ind = DIM_LEN
        ELSE ! If no rounding possible; a real error
          WRITE (Message,'(A,F7.2,A,F7.2,A,F7.2,A,F7.2)') 'Input value ', t_val, &
                  ' is not within range [',t_min,', ',t_max,'] with step', t_diff
          CALL Fatal('GridDataMapper', Message)
        END IF
      ELSE  
        t_ind = ((t_val - t_min)/t_diff) + 1 ! Uniform grid: just remove the bias and normalize the difference out
      END IF
      ! No rounding for it is interpolated later on
      WRITE (Message, '(A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2)') 'Time index for given value ', &
                          t_val, ' is ', t_ind, ' over range [', t_min,',',t_max,'] with step ', t_diff
      CALL Info('GridDataMapper', Message)
  
    END FUNCTION TimeValueToIndex
 
    !-------------------- G_Error() ------------------------
    !----- Checks the status and if failure, prints the error message and returns .TRUE.
    !-------------------------------------------------------
    FUNCTION G_Error( status, msg ) RESULT(erred)
      IMPLICIT NONE
      
      !----- Declarations
      INTEGER, INTENT(IN) :: status ! Status value
      CHARACTER (len = *), INTENT(IN) :: msg ! Error details
      LOGICAL :: erred ! True, if all ok; False otherwise
      
      !----- Checks for errors
      erred = .FALSE.
      IF ( status .NE. NF90_NOERR ) THEN ! Error encountered
        erred = .TRUE.
        CALL Fatal( 'GridDataMapper', msg )
      END IF
    
    END FUNCTION G_Error

END MODULE NetCDFGeneralUtils

