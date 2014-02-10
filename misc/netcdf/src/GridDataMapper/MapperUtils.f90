!------------------------------------------------------------------------------
! Vili Forsell
! Created: 7.7.2011
! Last Modified: 13.7.2011
!------------------------------------------------------------------------------
! This module contains general functions for GridDataMapper
! - GetElmerNodeValue() ; Simplified access to Elmer Node's coordinates
! - GetElmerMinMax() ; Returns min, or max, value for given Elmer dimension
! - IntWidth() ; Returns the numeric width for an integer value
!     o Mainly used for determining widths for format fields to avoid prefix spaces
! - ListGetStrings() ; Collects an array of strings from SIF defined with names suffixed by running numbers
!     o Used to generalize string input
!------------------------------------------------------------------------------
MODULE MapperUtils

CONTAINS

  !---------------- GetElmerNodeValue() ---------------
  !--- Gets the value of the chosen node from the given dimension (1 = x, 2 = y, 3 = z)
  FUNCTION GetElmerNodeValue( Solver, node, dimE ) RESULT( node_val )
  !----------------------------------------------------
    USE DefUtils, ONLY: Solver_t, dp
    USE Messages, ONLY: Fatal
    IMPLICIT NONE

    TYPE(Solver_t), INTENT(IN) :: Solver
    INTEGER, INTENT(IN) :: node, dimE
    REAL(KIND=dp) :: node_val ! The output
    SELECT CASE (dimE)
      CASE (1)
        node_val = Solver % Mesh % Nodes % x(node)
      CASE (2)
        node_val = Solver % Mesh % Nodes % y(node)
      CASE (3)
        node_val = Solver % Mesh % Nodes % z(node)
      CASE DEFAULT
        CALL Fatal('GridDataMapper','GetElmerNodeValue(): Elmer dimension not found')
        node_val = 0
    END SELECT
     
  END FUNCTION GetElmerNodeValue

  !---------------- GetElmerMinMax() ---------------
  !--- Gets the minimum/maximum (chosen) value of the given dimension (1 = x, 2 = y, 3 = z)
  FUNCTION GetElmerMinMax( Solver, dimE, GET_MIN ) RESULT( node_val )
  !----------------------------------------------------
    USE DefUtils, ONLY: Solver_t, CoordinateSystemDimension, dp
    USE Messages, ONLY: Fatal
    IMPLICIT NONE

    TYPE(Solver_t), INTENT(IN) :: Solver
    INTEGER, INTENT(IN) :: dimE
    LOGICAL, INTENT(IN) :: GET_MIN
    REAL(KIND=dp) :: node_val ! The output

    SELECT CASE (dimE)
      CASE (1)
        IF ( GET_MIN ) THEN 
          node_val = MINVAL(Solver % Mesh % Nodes % x)
        ELSE 
          node_val = MAXVAL(Solver % Mesh % Nodes % x)
        END IF
      CASE (2)
        IF ( GET_MIN ) THEN
          node_val = MINVAL(Solver % Mesh % Nodes % y)
        ELSE 
          node_val = MAXVAL(Solver % Mesh % Nodes % y)
        END IF
      CASE (3)
        IF ( GET_MIN ) THEN
          node_val = MINVAL(Solver % Mesh % Nodes % z)
        ELSE 
          node_val = MAXVAL(Solver % Mesh % Nodes % z)
        END IF
      CASE DEFAULT
        CALL Fatal('GridDataMapper','GetAllElmerNodeValues(): Elmer dimension not found')
        node_val = 0
    END SELECT
     
  END FUNCTION GetElmerMinMax


  !----------------------- IntWidth() --------------
  !--- Finds the width of an integer; ignores sign
  FUNCTION IntWidth( NR ) RESULT( width )
  !-------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NR
    INTEGER :: val
    INTEGER :: width
    INTEGER :: comp

    width = 1
    comp = 10
    val = ABS(NR) + 1 ! Ignores sign; val >= 1

    ! Note that:
    ! 10^0 - 1 = 0 <= 0..9 <= 10 = 10^1 - 1,
    ! 10^1 - 1 = 0 <= 10..99 <= 10 = 10^2 - 1,
    ! 10^2 - 1 = 0 <= 100..999 <= 10 = 10^3 - 1, 
    ! and so forth,
    ! where 10 corresponds to width of the number with 10 based numbers (n == width).
    ! We get: 10^(n-1) <= NR + 1 <= 10^(n)

    ! Width usually small, so logarithmic search is not usually necessary.

    ! P: 10^0 = 1 <= val+1 .AND. width = 1 => For every j; 0 <= j <= i: 10^i <= val+1 .AND. i = 0 .AND. width = i + 1 = 0 + 1 = 1
    DO WHILE (comp < val)
      width = width + 1 ! width = i+1
      comp = 10*comp ! 10*10^(i) = 10^(i+1)
      ! I: For every j; 0 <= j <= i: 10^j <= val+1 .AND. Exists n; i < n: val+1 <= 10^n .AND. width = i+1
    END DO
    ! Q: (For every i; 0 <= i <= n-1: 10^i <= val+1) .AND. val+1 <= 10^n .AND. width = n
    !    => 10^(n-1) <= val+1 .AND. val+1 <= 10^n .AND. width = n

  END FUNCTION IntWidth

  !----------------- GetNetCDFAccessParameters() ------
  !--- Gets all access strings and their associated accessing order
  SUBROUTINE GetNetCDFAccessParameters( List,Variables,Constants,Permutation,Found )
  !----------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    !--- Arguments
    TYPE(ValueList_t), INTENT(IN), POINTER :: List ! A pointer to a list of values
    LOGICAL, OPTIONAL, INTENT(INOUT) :: Found(2) ! True, if found
    CHARACTER(LEN=MAX_NAME_LEN), INTENT(OUT), ALLOCATABLE :: Variables(:),Constants(:) ! For returned array data
    INTEGER, INTENT(OUT), ALLOCATABLE :: Permutation(:) ! NetCDF Access permutation; maps variables and constants to their proper locations
    ! Indexing order: Variables, then Constants

    !--- Variables
    TYPE(ValueList_t), POINTER :: ptr
    CHARACTER(LEN=MAX_NAME_LEN) :: tmpStr
    CHARACTER(LEN=MAX_NAME_LEN) :: Names(2) ! Names of the SIF variables
    LOGICAL, ALLOCATABLE :: IsVariable(:) ! Valid until total_size
    LOGICAL :: GotIt, PrevGotIt, OtherWasDefined
    INTEGER :: TOTAL, var_i, const_i, const_size, coord_size
    INTEGER :: loop, loop2, total_size, alloc_stat
    CHARACTER(LEN=10) :: tmpFormat

    !--- Initializations
    NULLIFY(ptr)
    IF ( PRESENT(Found) ) Found = .FALSE.

    TOTAL = GetInteger(List,'NetCDF Max Parameters', GotIt)
    IF ( .NOT. GotIt ) THEN
      CALL Warn('GridDataMapper',&
     'Please specify the maximum amount of NetCDF parameters with &
variable "NetCDF Max Parameters". Assumed to be 10 by default.')
      TOTAL = 10 ! Default 10, affects efficiency and halting
    ELSE IF ( TOTAL .LT. 1 ) THEN
      CALL Warn('GridDataMapper','Minimum of one parameter expected for "NetCDF Max Parameters". It is now set to 1.')
      TOTAL = 1
    END IF

    ALLOCATE ( IsVariable(TOTAL), STAT = alloc_stat )
    IF ( alloc_stat .NE. 0 ) THEN
      CALL Fatal('GridDataMapper','Memory ran out')
    END IF
    IsVariable = .FALSE. ! Defaulted to constant; even when over range

    ! The names of the SIF variables
    Names(1) = 'Coordinate Name'
    Names(2) = 'NetCDF Constant'

    ! Count the amount of defined strings for allocation
    total_size = 0
    coord_size = 0 ! Amount of existing coordinates/variables
    const_size = 0 ! Amount of existing constants
    PrevGotIt = .TRUE. ! If previous one was false, then the next one must not be true (otherwise omitted numbers)
    OtherWasDefined = .FALSE. ! True if both names defined for the same location
    DO loop = 1,TOTAL,1
      OtherWasDefined = .FALSE. ! (NO DUPLICATES)
      DO loop2 = 1,size(Names),1
        ! Find all Coordinates and Constants

        !--- Checks for existence of the Name (EXISTS CHECK)
        ! Tries to ensure that the given integer doesn't have extra spaces before it in char format
        ! Scales until a number with a width of 9 (limited by tmpFormat)
        WRITE(tmpFormat,'(A,I1,A)') '(A,A,I', IntWidth(loop),')'
        WRITE(tmpStr,tmpFormat) TRIM(Names(loop2)),' ',loop
        ptr => ListFind( List,tmpStr,GotIt ) ! NOTE: This will probably be a slow operation (linear/logarithmic)
 
        !--- Checks for duplicates and updates the sizes
        IF (GotIt) THEN
          IF (.NOT. ASSOCIATED(ptr)) THEN
            GotIt = .FALSE.
          ELSE IF ( OtherWasDefined ) THEN
            !--- The other Name was defined for the same index; conflict! (NO DUPLICATES)
            WRITE(Message,'(A,I3)') 'Both a Coordinate Name and a NetCDF Constant &
were defined for NetCDF variable access location ', loop
            CALL Fatal('GridDataMapper', Message)
          ELSE
            !--- Found and works (COORD_SIZE RIGHT, CONST_SIZE RIGHT, TOTAL_SIZE CORRECT)
            ! < Exists "Names 'loop2'" >
            total_size = loop ! Chooses the largest that exists
            IF ( loop2 .EQ. 1 ) THEN ! "Names" = "Coordinate Name"
              coord_size = coord_size + 1
              IsVariable(loop) = .TRUE. ! To find correct places later on
            ELSE ! "Names" = "NetCDF Constant"
              const_size = const_size + 1
              ! IsVariable defaulted to .FALSE.
            END IF
            ! Postcondition has been updated for all until loop and loop2 (NO DUPLICATES)
            OtherWasDefined = .TRUE.
          END IF
        END IF

        !--- Checks that there are no gaps with the numbering (TOTAL_SIZE CORRECT)
        IF ( GotIt .AND. (.NOT. PrevGotIt) ) THEN
          WRITE(Message,'(A,I3,A)') 'NetCDF access parameter(s) before parameter ', loop  ,' are missing'
          CALL Fatal('GridDataMapper',Message)
        END IF
      END DO

      !--- Both of the Names have been handled, so PrevGotIt can be updated for the next round (TOTAL_SIZE CORRECT)
      ! o If total_size < loop at the end(!) of loop, then there have been omissions; holds also for first round
      ! o Then, there must be no subsequent found parameters until the limit TOTAL; else there is a gap, which aborts immediately
      ! o NOTE: GotIt is not reliable enough for this check! F.ex. first name found, second not, would imply that nothing has been found.
      IF ( total_size .EQ. loop ) THEN 
        PrevGotIt = .TRUE.
      ELSE
        PrevGotIt = .FALSE.
      END IF
    END DO
    ! Postcondition for the loops:
    !< COORD_SIZE RIGHT: coord_size = |{ x | Exists "Coordinate Name 'x'" }| >
    !< CONST_SIZE RIGHT: const_size = |{ x | Exists "NetCDF Constant 'x'" }| >
    !< TOTAL_SIZE CORRECT: ( total_size = t .AND. (t = 0 .OR. (For Every i; 1 =< i <= t: Exists "Names 'i'")) .AND. ( .NOT. Exist i; t < i <= TOTAL: Exists "Names 'i'") >
    !< NO DUPLICATES: ( .NOT. Exists i; 1 =< i <= TOTAL: Exists "NetCDF Constant 'i'" .AND. Exists "Coordinate name 'i'" ) >

    IF ( total_size .LE. 0 ) RETURN ! No strings found

    !--- Input checked and sizes counted; now allocation and data retrieval

    ! Allocation
    ALLOCATE ( Variables(coord_size),Constants(const_size),Permutation(total_size), STAT = alloc_stat )
    IF ( alloc_stat .NE. 0 ) THEN
      CALL Fatal('GridDataMapper','Memory ran out')
    END IF
    Permutation = 0

    ! Getting the strings and the permutation data
    var_i = 0
    const_i = 0
    DO loop = 1,total_size,1
      !--- Pinpoint the correct type (correct until total_size)
      IF ( IsVariable(loop) ) THEN
        loop2 = 1 ! Variable/Coordinate
      ELSE
        loop2 = 2 ! Constant
      END IF

      WRITE(tmpFormat,'(A,I1,A)') '(A,A,I', IntWidth(loop) ,')'
      WRITE(tmpStr,tmpFormat) TRIM(Names(loop2)),' ',loop

      !--- Puts the data of the Names(loop2) in the right place
      IF ( IsVariable(loop) ) THEN 
        var_i = var_i + 1
        Variables(var_i) = GetString( List,tmpStr,GotIt )
        Permutation(var_i) = loop
      ELSE
        const_i = const_i + 1
        Constants(const_i) = GetString( List,tmpStr,GotIt )
        Permutation(coord_size + const_i) = loop
      END IF

      IF ( .NOT. GotIt ) THEN
        CALL Fatal('GridDataMapper','Obtained string did not exist after all')
      END IF
    END DO
    !--- Permutation's contents:
    ! Range 1 <= i <= coord_size: variable i's location during the loop
    ! Range coord_size + 1 <= j <= total_size: constant j's location during the loop
    ! In other words, indexing Permutation gives the right NetCDF access location for the access parameter

    IF ( coord_size .GT. 0 ) Found(1) = .TRUE.
    IF ( const_size .GT. 0 ) Found(2) = .TRUE.

  END SUBROUTINE GetNetCDFAccessParameters 

  !----------------- ListGetStrings() -----------------
  !--- Gets all strings defined with prefix "Name" and ending with " NR", where NR is a number in an array
  SUBROUTINE ListGetStrings( List,Name,Found,CValues )
  !----------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    !--- Arguments
    TYPE(ValueList_t), INTENT(IN), POINTER :: List ! A pointer to a list of values
    CHARACTER(LEN=*), INTENT(IN) :: Name ! Name of the SIF variable
    LOGICAL, OPTIONAL, INTENT(INOUT) :: Found ! True, if found
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: CValues(:) ! For returned array data

    !--- Variables
    TYPE(ValueList_t), POINTER :: ptr
    CHARACTER(LEN=MAX_NAME_LEN) :: tmpStr
    LOGICAL :: GotIt
    INTEGER :: loop, amount, alloc_stat
    CHARACTER(LEN=10) :: tmpFormat

    !--- Initializations
    NULLIFY(ptr)
    Found = .FALSE.

    ! Count the amount of defined strings for allocation
    amount = 0
    loop = 1
    GotIt = .TRUE.
    DO WHILE (GotIt)
      ! Tries to ensure that the given integer doesn't have extra spaces before it in char format
      ! Scales until a number with a width of 9 (limited by tmpFormat)
      WRITE(tmpFormat,'(A,I1,A)') '(A,A,I', IntWidth(loop),')'
      WRITE(tmpStr,tmpFormat) TRIM(Name),' ',loop
      ptr => ListFind( List,tmpStr,GotIt )

!      WRITE(*,*) 'TEMP: ', tmpStr

      ! Continues until first name is not found
      IF (GotIt) THEN
        IF (.NOT. ASSOCIATED(ptr)) THEN
          GotIt = .FALSE.
          RETURN
        ELSE
          amount = amount + 1
          loop = loop + 1
        END IF
      END IF
    END DO

    IF ( amount .LE. 0 ) RETURN ! No strings found

    ! Allocation
    ALLOCATE ( CValues(amount), STAT = alloc_stat )
    IF ( alloc_stat .NE. 0 ) THEN
      CALL Fatal('GridDataMapper','Memory ran out')
    END IF

    ! Getting the strings
    DO loop = 1,amount,1
      WRITE(tmpFormat,'(A,I1,A)') '(A,A,I', IntWidth(loop) ,')'
      WRITE(tmpStr,tmpFormat) TRIM(Name),' ',loop
      CValues(loop) = GetString( List,tmpStr,GotIt )
      IF (.NOT. GotIt) THEN
        CALL Fatal('GridDataMapper','Obtained string did not exist after all')
      END IF
    END DO
   
    Found = .TRUE.

  END SUBROUTINE ListGetStrings 

END MODULE MapperUtils
