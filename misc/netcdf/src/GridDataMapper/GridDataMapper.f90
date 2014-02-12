!------------------------------------------------------------------------------
! Peter Råback, Vili Forsell
! Created: 7.6.2011
! Last Modified: 4.8.2011
!------------------------------------------------------------------------------
SUBROUTINE GridDataMapper( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  This subroutine saves scalar values to a matrix.
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************

!******************************************************************************
! Notes:	o Performs a modification from measurements with uniform grid into Elmer with an unstructured grid
!		o Remember to see if incremental processing might be possible
!		o Some terminology: "nodes" and "edges" used for grid points and the lines connecting the for grid points and the lines connecting them
!******************************************************************************

  USE DefUtils, ONLY: dp, Solver_t, Model_t, Mesh_t,GetInteger, CoordinateSystemDimension, GetSolverParams, &
                      GetLogical, MAX_NAME_LEN, GetCReal
  USE Messages, ONLY: Info, Warn, Fatal, Message
  USE NetCDF
  USE NetCDFGridUtils, ONLY: UniformGrid_t, PrintGrid
  USE NetCDFInterpolate, ONLY: Interpolate, LinearInterpolation
  USE NetCDFGeneralUtils, ONLY: CloseNetCDF, TimeType_t
  USE MapperUtils, ONLY: GetElmerMinMax, GetElmerNodeValue
  USE CustomTimeInterpolation, ONLY: ChooseTimeInterpolation

  IMPLICIT NONE

  !------------------------------------------------------------------------------
  LOGICAL, PARAMETER :: DEBUG = .FALSE. ! Shows the basic debug info on grids and dimensions
  LOGICAL, PARAMETER :: DEBUG_MORE = .FALSE. ! Shows also debug printouts for each iteration
  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(TimeType_t) :: Time
  TYPE(UniformGrid_t) :: Grids(2)
  INTEGER :: k, node, DIM, MAX_STEPS
  INTEGER, POINTER :: FieldPerm(:)
  REAL(KIND=dp), POINTER :: Field(:)
  REAL(KIND=dp), ALLOCATABLE :: x(:)
  REAL(KIND=dp) :: interp_val, interp_val2, u1(2), u2(2)
  INTEGER, ALLOCATABLE :: dim_lens(:) ! Lengths for all dimensions
  INTEGER :: NCID, loop, Var_ID, alloc_stat
  CHARACTER (len = MAX_NAME_LEN) :: Coord_System, TimeInterpolationMethod
  REAL(KIND=dp) :: InterpMultiplier, InterpBias
  LOGICAL :: output_on

  !------------------------------------------------------------------------------
  ! General initializations
  !------------------------------------------------------------------------------

  CALL Info('GridDataMapper','-----------------------------------------', Level=4 )
  CALL Info('GridDataMapper','Getting field from grid data',Level=4) 

  Time % val = -1.0_dp
  Time % id = -1
  Time % len = -1
  Time % low = -1
  Time % high = -1
  Time % doInterpolation = .FALSE.
  Time % is_defined = .FALSE.
  
!  WRITE (*,*) 'val ', Time % val, ' id ', Time % id, ' len ', Time % len,&
! ' low ', Time % low, ' high ', Time % high, ' interp ', Time % doInterpolation

  IF ( TransientSimulation ) THEN
    MAX_STEPS = GetInteger( Model % Simulation, 'TimeStep Intervals' )
  ELSE
    MAX_STEPS = 1 ! Steady state
  END IF

  !-- Pointer declarations
  Mesh => Solver % Mesh
  Field => Solver % Variable % Values  ! This vector will get the field values now
  FieldPerm => Solver % Variable % Perm

  !------------------------------------------------------------------------------
  ! Initializing NetCDF values
  !------------------------------------------------------------------------------
  CALL InitNetCDF(Solver, NCID, Var_ID, dim_lens, Grids, Time, TransientSimulation, dt, MAX_STEPS, Coord_System)
  DIM = Grids(1) % COORD_COUNT ! Abbreviation [# of Elmer coordinates .EQ. # of NetCDF coordinates]

  !------------------------------------------------------------------------------
  ! Initializing Elmer mesh vectors, scaling and interpolation
  !------------------------------------------------------------------------------
  IF ( DIM .GT. 0 ) THEN
    ALLOCATE ( x(DIM), STAT = alloc_stat ) 
    IF ( alloc_stat .NE. 0 ) THEN
      CALL Fatal('GridDataMapper','Memory ran out')
    END IF

    !--- Scales the Grids, if applicable
    CALL InitScaling(Solver, Grids)
  END IF
 
  !------ Debug printouts -------------------------
  IF (DEBUG) THEN
    IF ( DIM .GT. 0 ) THEN
      PRINT *,'Initial Elmer Grid Bounding Box:'
      DO loop = 1,DIM,1
        PRINT *,'Coordinate ', loop,':', GetElmerMinMax(Solver,loop,.TRUE.), GetElmerMinMax(Solver,loop,.FALSE.)
      END DO
    END IF

    DO loop = 1,size(Grids,1),1
      CALL PrintGrid(Grids(loop),loop)
    END DO    
  END IF
  !------------------------------------------------

  !--- Initializes the interpolation variables
  CALL InitInterpolation( Solver, DIM, Grids, TimeInterpolationMethod, InterpMultiplier, InterpBias )

  !------------------------------------------------------------------------------
  ! INTERPOLATION LOOP
  !------------------------------------------------------------------------------
  ! NOTE: Interpolate() effectively ignores time, if "Time % is_defined" is .FALSE.
  !------------------------------------------------------------------------------

  output_on = .TRUE. ! If true, outputs some iteration information
  ! Go through the active nodes and perform interpolation
  DO node=1, Mesh % NumberOfNodes
    k = FieldPerm(node) 
    IF( k == 0 ) CYCLE

    IF ( DIM .GT. 0 ) THEN

      ! The point of interest
      DO loop = 1,DIM,1
        x(loop) = GetElmerNodeValue(Solver,node,Grids(1) % Elmer_perm(loop))
      END DO
 
      !--- Perform time interpolation, if needed 
      IF ( Time % doInterpolation ) THEN ! Two time values (otherwise not interpolated)
        IF ( .NOT. (Interpolate(Solver,NCID,Var_ID,dim_lens,Grids(1),& 
          Time, Time % low,interp_val,Coord_System,x) .AND. Interpolate(Solver,NCID,&
          Var_ID,dim_lens,Grids(2),Time,Time % high,interp_val2,Coord_System,x)) ) THEN
          CYCLE
        ELSE
          ! Time interpolation on already interpolated space values; save result in interp_val, use original time to weigh
          u1(1) = Time % low
          u1(2) = interp_val
          u2(1) = Time % high
          u2(2) = interp_val2
          interp_val = ChooseTimeInterpolation(Time % val,u1,u2,TimeInterpolationMethod,output_on) ! Chooses the time interpolation method
          ! See: CustomTimeInterpolation.f90
          output_on = .FALSE.
        END IF
      ELSE ! Holds: Time % low = Time % high = Time % val, also: Time % IS_DEFINED = .FALSE. works!
        IF (.NOT. Interpolate(Solver,NCID,Var_ID,dim_lens,Grids(1),&
                  Time,Time % low,interp_val,Coord_System,x) ) THEN
          CYCLE ! Ignore values for incompatible interpolation
        END IF
      END IF
    ELSE

      !------------------------------------------------------------------------------
      ! No coordinate dimensions; only time and constant dimensions.
      !------------------------------------------------------------------------------
      IF ( Time % doInterpolation ) THEN
        !--- Time interpolation in effect, perform as before
        IF (.NOT. (Interpolate(Solver,NCID,Var_ID,dim_lens,Grids(1),&
                     Time,Time % low,interp_val,Coord_System) .AND. &
                  Interpolate(Solver,NCID,Var_ID,dim_lens,Grids(2),&
                     Time,Time % high,interp_val2,Coord_System)) ) THEN
          CYCLE
        ELSE
          u1(1) = Time % low
          u1(2) = interp_val
          u2(1) = Time % high
          u2(2) = interp_val2
          interp_val = ChooseTimeInterpolation(Time % val,u1,u2,TimeInterpolationMethod,output_on)
          ! See: CustomTimeInterpolation.f90
          output_on = .FALSE.
        END IF
      ELSE
        !--- No time interpolation
        IF (.NOT. Interpolate(Solver,NCID,Var_ID,dim_lens,Grids(1),&
                     Time,Time % low,interp_val,Coord_System) ) THEN
          CYCLE
        END IF
      END IF
  
    END IF

    !------ Debug printouts -------------------------
    IF (DEBUG_MORE) THEN
       PRINT *,'Interpolation result: ', interp_val
    END IF
    !------------------------------------------------

    Field(k) = InterpMultiplier*interp_val + InterpBias ! Doesn't modify the result by default
  END DO

  !------------------------------------------------------------------------------
  ! Close the NetCDF file
  !------------------------------------------------------------------------------
  CALL CloseNetCDF(NCID)

  !------------------------------------------------------------------------------

  CALL Info('GridDataMapper','All done',Level=4)
  CALL Info('GridDataMapper', '-----------------------------------------', Level=4 )

CONTAINS

  !----------------- InitScaling() --------------------
  !--- Adds scaling to NetCDF Grids (by default does nothing)
  SUBROUTINE InitScaling(Solver, Grids)
  !----------------------------------------------------
    USE DefUtils, ONLY: Solver_t
    USE NetCDFGridUtils, ONLY: UniformGrid_t
    USE MapperUtils, ONLY: GetElmerMinMax
    IMPLICIT NONE
    
    !------------------------------------------------------------------------------
    ! ARGUMENTS
    !------------------------------------------------------------------------------
    TYPE(Solver_t), INTENT(IN) :: Solver
    TYPE(UniformGrid_t), INTENT(INOUT) :: Grids(:)

    !------------------------------------------------------------------------------
    ! VARIABLES
    !------------------------------------------------------------------------------
    LOGICAL :: Found, ENABLE_SCALING
    INTEGER :: DIM, alloc_stat
    REAL(KIND=dp), ALLOCATABLE :: x0e(:), x1e(:)

    !------------------------------------------------------------------------------
    ! Adds scaling, if it's set on and possible
    !------------------------------------------------------------------------------

    ENABLE_SCALING = GetLogical(GetSolverParams(Solver), "Enable Scaling", Found)
    DIM = Grids(1) % COORD_COUNT

    IF ( Found .AND. ENABLE_SCALING .AND. DIM .GT. 0 ) THEN

      CALL Warn('GridDataMapper','Elmer grid is scaled to match the NetCDF grid')

      ALLOCATE (x0e(DIM),x1e(DIM), STAT=alloc_stat)
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF

      !--- Collects the range of the Elmer mesh bounding box for scaling
      DO loop = 1,DIM,1
        x0e(loop) = GetElmerMinMax(Solver,loop,.TRUE.)
        x1e(loop) = GetElmerMinMax(Solver,loop,.FALSE.)
      END DO

      !--- Scales
      DO loop = 1,size(Grids,1),1
        ! First the scaling to same size (Eq. a( X1E(1)-X0E(1) ) = (X1(1)-X0(1)) ; ranges over a dimension are same. Solved for a, 1 if equal)
        Grids(loop) % scale(:) = (Grids(loop) % X1(1:DIM) - Grids(loop) % X0(1:DIM))/(X1E(1:DIM)-X0E(1:DIM)) ! Note: "/" and "*" elementwise operations for arrays in Fortran
        ! Second the vector to reach X0 from the scaled X0E (wherever it is)
        Grids(loop) % move(:) = Grids(loop) % X0(1:DIM) - Grids(loop) % scale(:)*X0E(1:DIM) ! zero, if equal
      END DO
    END IF
    !------------------------------------------------------------------------------

  END SUBROUTINE InitScaling

  !----------------- InitNetCDF() ---------------------
  !--- Gathers and initializes all the necessary NetCDF information for picking variables
  SUBROUTINE InitNetCDF(Solver, NCID, Var_ID, dim_lens, Grids, Time,&
                                   IS_TRANSIENT, STEP_SIZE, MAX_STEPS, Coord_System )
  !--------------------------------------------------

    USE NetCDFGridUtils, ONLY: PrintGrid, InitGrid, GetNetCDFGridParameters, Focus2DNetCDFGrid 
    USE NetCDFGeneralUtils, ONLY: GetAllDimensions, G_Error, TimeValueToIndex
    USE MapperUtils, ONLY: IntWidth, GetNetCDFAccessParameters
    USE NetCDF
    USE DefUtils, ONLY: dp, MAX_NAME_LEN, GetSolverParams, GetString, GetConstReal, ListGetString
    USE Messages, ONLY: Fatal, Message
    IMPLICIT NONE

    !------------------------------------------------------------------------------
    ! ARGUMENTS
    !------------------------------------------------------------------------------
    TYPE(Solver_t), INTENT(IN) :: Solver ! Elmer solver from SIF
    TYPE(TimeType_t), INTENT(INOUT) :: Time ! Time information
    TYPE(UniformGrid_t), INTENT(INOUT) :: Grids(:) ! NetCDF grids
    LOGICAL, INTENT(IN) :: IS_TRANSIENT ! For using Elmer time instead of a variable
    REAL(KIND=dp), INTENT(IN) :: STEP_SIZE ! Time step for Elmer transient system
    INTEGER, INTENT(IN) :: MAX_STEPS ! Max steps in Elmer transient system
    INTEGER, INTENT(OUT) :: NCID, Var_ID  !  NCID is the ID of the opened file, Var_ID the accessed variable id
    INTEGER, INTENT(INOUT), ALLOCATABLE :: dim_lens(:) ! Lengths for all dimensions
    CHARACTER (len = MAX_NAME_LEN), INTENT(OUT) :: Coord_System ! Coordinate system string

    !------------------------------------------------------------------------------
    ! NetCDF VARIABLES
    !------------------------------------------------------------------------------
    
    INTEGER :: DIM ! Used Elmer dimensions for data accessing
    INTEGER, ALLOCATABLE :: dim_ids(:) ! Ids for all dimensions
    LOGICAL :: Found(7) ! True if SIF definitions found
    CHARACTER (len = MAX_NAME_LEN) :: FileName ! File name for reading the data (of .nc format)
    CHARACTER (len = MAX_NAME_LEN) :: Var_Name, T_Name, Mask_Name
    CHARACTER (len = MAX_NAME_LEN), ALLOCATABLE :: Dim_Names(:) ! Contains the opened NetCDF dimension variables
    REAL(KIND=dp) :: Mask_Limit ! Value limit where masking starts to take effect
    INTEGER :: loop, alloc_stat,dim_count,status ! Status tells whether operations succeed
    INTEGER :: size_coord, size_const ! Temps for amounts of NetCDF coordinate and constant dimensions
    LOGICAL :: tmpBool ,IsTimeDependent ! True, if time is used when accessing NetCDF
    CHARACTER (len = MAX_NAME_LEN), ALLOCATABLE :: Coords(:), Constants(:) ! Contains the names for NetCDF coordinate and constant dimensions
    INTEGER, ALLOCATABLE :: Permutation(:) ! Contains the NetCDF Access permutation for Coords and Constants
    CHARACTER (len = 10) :: tmpFormat
  
    !------------------------------------------------------------------------------
    ! NetCDF INITIALIZATIONS
    !------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! 1) Gathers basic data from SIF
    !-------------------------------------------------------------------------------

    !--- Collects the input information from Solver Input File
    FileName = GetString( GetSolverParams(Solver), "File Name", Found(1) )
    Var_Name = GetString( GetSolverParams(Solver), "Var Name", Found(2) )

    !--- Collects the NetCDF accessing information (coordinates, constant dimensions, time)
    CALL GetNetCDFAccessParameters( GetSolverParams(Solver), Coords, Constants, Permutation, Found(3:4) )
    T_Name = GetString( GetSolverParams(Solver), "Time Name", IsTimeDependent ) ! If given, time is the last dimension
  
    ! Following parameters are needed for masking and coordinate system
    Mask_Name = GetString( GetSolverParams(Solver), "Mask Variable", Found(5) )
    Mask_Limit = GetConstReal( Solver % Values, "Mask Limit", Found(6) )
    Coord_System = GetString( GetSolverParams(Solver), "Coordinate System", Found(7) ) ! Any input is ok; only valid gives a conversion and error is given otherwise
    IF ( .NOT. Found(7) ) Coord_System = '' ! Initializes Coord_System, if not given

    !-------------------------------------------------------------------------------
    ! 2) Opening the NetCDF file and finding the variable
    !-------------------------------------------------------------------------------

    DO loop = 1,2,1
      IF ( .NOT. Found(loop) ) THEN ! Checks that the constants have been found successfully
        CALL Fatal('GridDataMapper', &
      "Unable to find a compulsory NetCDF Name Constant (the name of file or variable)")
      END IF
    END DO

    status = NF90_OPEN(FileName,NF90_NOWRITE,NCID) ! Read-only
    IF ( G_Error(status, "NetCDF file could not be opened") ) THEN
      CALL abort() ! End execution
    END IF

    ! Find variable to be accessed
    status = NF90_INQ_VARID(NCID,Var_Name,Var_ID)
    IF ( G_Error(status,'NetCDF variable name not found.') ) THEN
      CALL abort()
    END IF


    !-------------------------------------------------------------------------------
    ! 3) Finds the amounts for used coordinate and constant dimensions
    !-------------------------------------------------------------------------------

    !--- Defining coordinate, constant dimension and total dimension sizes
    size_coord = 0
    size_const = 0
    IF ( Found(3) ) size_coord = size(Coords,1)
    IF ( Found(4) ) size_const = size(Constants,1)
    dim_count = size_coord + size_const

    IF ( (dim_count .EQ. 0) .AND. (.NOT. IsTimeDependent) ) THEN
      CALL Fatal('GridDataMapper','Expected at least one indexing variable: time, coordinate, or constant.')
    END IF
    !> dim_count in range ( 0, (size_coord + size_const) ), and at least one possible NetCDF accessing parameter given

    !------------------------------------------------------------------------------
    ! 4) Obtains the used Elmer coordinate dimensionality (DIM)
    !------------------------------------------------------------------------------

    DIM = size_coord

    !-------------------------------------------------------------------------------
    ! 5) Forms an array of dimension names (constant and coordinate) for accessing
    !-------------------------------------------------------------------------------
    WRITE (Message,'(A,I3,A)') 'Using ',dim_count,' input dimensions.'
    CALL Info('GridDataMapper',Message)

    IF ( dim_count .GT. 0 ) THEN
   
      ! Form an array of wanted names for defining the information 
      ALLOCATE (Dim_Names(dim_count), dim_ids(dim_count), dim_lens(dim_count), STAT = alloc_stat)
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF
  
      ! Initialization in case of errors
      Dim_Names = 'none'
      dim_ids = -1
      dim_lens = -1
  
      ! Adds the coordinate NetCDF dimension names, if available
      IF ( size_coord .GT. 0 ) Dim_Names(1:size(Coords)) = Coords(:)
  
      ! Adds the constant NetCDF dimension names, if available
      IF ( size_const .GT. 0 ) THEN
        IF ( size_coord .GT. 0 ) THEN
          Dim_Names((size(Coords)+1):size(Dim_Names)) = Constants(:)
        ELSE
          Dim_Names(:) = Constants(:)
        END IF
      END IF
  
    END IF
    !> If there are dimensions, their names are in Dim_Names = [Coords, Constants]

    !-------------------------------------------------------------------------------
    ! 6) Gets dimensions on basis of the given names
    !-------------------------------------------------------------------------------
    IF ( dim_count .GT. 0 ) CALL GetAllDimensions(NCID,Dim_Names,dim_ids,dim_lens)
    !> dim_ids & dim_lens collected with the Dim_Names 

    !-------------------------------------------------------------------------------
    ! 7) Initializes the Grids (1/2)
    !------------------------------------------------------------------------------

    !--- Default initializations for all cases
    DO loop = 1,size(Grids,1),1
      CALL InitGrid(Grids(loop), dim_count, size_coord )
      Grids(loop) % access_perm = Permutation
    END DO

    !-------------------------------------------------------------------------------
    ! Connects the Elmer Mesh dimensions to NetCDF access parameters
    !-------------------------------------------------------------------------------


    !--- Initializations for the case of having some real content
    IF ( Grids(1) % IS_DEF ) THEN

      !--- Obtains the preliminary definining parameters for the NetCDF grid
      CALL GetNetCDFGridParameters(NCID,Grids(1),dim_ids,dim_lens) ! Normal grid parameters don't depend on time

      !--- Connects every NetCDF access coordinate to a Elmer coordinate
      DO loop = 1,Grids(1) % COORD_COUNT,1

        WRITE(tmpFormat,'(A,I1,A)') '(A,I',IntWidth(loop),',A)'
        WRITE(Message,tmpFormat) 'Coordinate ', Grids(1) % ACCESS_PERM( loop ), ' To Elmer Dimension'

        Grids(1) % Elmer_perm( loop ) = &
               GetInteger(GetSolverParams(Solver),Message,tmpBool)
        IF ( .NOT. tmpBool ) THEN
          WRITE(Message,'(A,I3,A)') 'Declared Coordinate Name ', Grids(1) % ACCESS_PERM( loop ),&
          ' is not connected to a Elmer dimension. Check all "Coordinate X To Elmer Dimension" definitions.'
          CALL Fatal('GridDataMapper',Message)
        END IF
      END DO 

      !--- Finds the index values for the found constant dimensions
      DO loop = 1,size(Grids(1) % const_vals),1
  
        WRITE(tmpFormat,'(A,I1,A)') '(A,I',IntWidth(loop),')'
        WRITE(Message,tmpFormat) 'NetCDF Constant Value ', Grids(1) % ACCESS_PERM( Grids(1) % COORD_COUNT + loop)

        Grids(1) % const_vals(loop) = GetInteger(GetSolverParams(Solver),Message,tmpBool)
        IF ( .NOT. tmpBool ) THEN
          WRITE(Message,'(A,I3,A)') 'Declared NetCDF Constant ', Grids(1) % ACCESS_PERM( Grids(1) % COORD_COUNT + loop),&
          ' has no corresponding index value. Check all NetCDF Constants.'
          CALL Fatal('GridDataMapper',Message)
        END IF
      END DO
      Grids(2) % const_vals = Grids(1) % const_vals
    END IF
    !> Grids initialized with constant values,x0,dx,nmax, and default for others; or Grid % IS_DEF is .FALSE.

    !-------------------------------------------------------------------------------
    ! 8) Initializes Time (if it is defined)
    !-------------------------------------------------------------------------------

    IF (IsTimeDependent) THEN

      ! Inform about the fate of time
      CALL Info('GridDataMapper','Time dimension taken into account.')

      Time % is_defined = .TRUE.
      CALL InitTime( Solver, NCID, T_Name, IS_TRANSIENT, STEP_SIZE, MAX_STEPS, Time )

      ! If there'll be interpolation, initialize both of the grids usable
      IF ( Grids(1) % IS_DEF ) THEN
        IF ( Time % doInterpolation ) THEN
          Grids(2) % x0(:) = Grids(1) % x0(:)
          Grids(2) % dx(:) = Grids(1) % dx(:)
          Grids(2) % nmax(:) = Grids(1) % nmax(:)
        END IF

        IF ( size_coord .EQ. 2 ) THEN
          IF ( Found(5) .AND. Found(6) ) THEN
            CALL Info('GridDataMapper','Two dimensional NetCDF grid focusing on basis of the given mask is in effect.')
            CALL Focus2DNetCDFGrid(NCID,Mask_Name,Mask_Limit,Grids(1),Time % low,dim_lens)
            IF ( Time % doInterpolation ) THEN ! Need to interpolate, return two different grid parameters
              
              Grids(2) % x0(:) = Grids(1) % x0(:) ! Copy the same results obtained for all grids
              Grids(2) % dx(:) = Grids(1) % dx(:)
              Grids(2) % nmax(:) = Grids(1) % nmax(:)
              CALL Focus2DNetCDFGrid(NCID,Mask_Name,Mask_Limit,Grids(2),Time % high,dim_lens)
            END IF
          ELSE
            ! No mask used
          END IF
        END IF
      END IF
    END IF
    !> Time initialized and Grids adjusted, if time-dependent focusing used

    !-------------------------------------------------------------------------------
    ! 9) Initializes the Grids (2/2) ( depends on Focus2DNetCDFGrid() )
    !-------------------------------------------------------------------------------
    IF ( Grids(1) % IS_DEF ) THEN
      DO loop = 1,size(Grids,1),1
        Grids(loop) % x1(:) = Grids(loop) % x0(:) +&
         (Grids(loop) % nmax(:)-1) * Grids(loop) % dx(:) ! In 3D case opposite points of the cube; if only one dimension, will be 0
      END DO

    END IF
    !> Grids initialized with proper x1's
 
    !-------------------------------------------------------------------------------
    !> Phase 1: Coord_System
    !> Phase 2: NCID, Var_ID
    !> Phase 4: DIM
    !> Phase 6: dim_ids, dim_lens
    !> Phase 8: Time
    !> Phase 9: Grids (scale(:), move(:) and eps(:) not touched)
    IF ( DEBUG ) THEN
      CALL PrintGrid(Grids(1),1)
      CALL PrintGrid(Grids(2),2)
    END IF
 
  END SUBROUTINE InitNetCDF


  !---------------- InitTime() ------------------------
  !--- Initializes the time values
  SUBROUTINE InitTime( Solver, NCID, T_Name, IS_TRANSIENT, STEP_SIZE, MAX_STEPS, TimeResult )
  !----------------------------------------------------
    USE NetCDFGeneralUtils, ONLY: TimeValueToIndex, GetDimension
    USE DefUtils, ONLY: MAX_NAME_LEN, GetSolverParams, GetLogical, GetCReal, GetTime, GetConstReal
    USE Messages, ONLY: Message, Info, Fatal
    IMPLICIT NONE
    !-------------------------------------------------------------------------------
    ! ARGUMENTS
    !-------------------------------------------------------------------------------
    !--- A) Input
    TYPE(Solver_t), INTENT(IN) :: Solver
    INTEGER, INTENT(IN) :: NCID
    LOGICAL, INTENT(IN) :: IS_TRANSIENT
    CHARACTER(len=MAX_NAME_LEN), INTENT(IN) :: T_Name
    REAL(KIND=dp), INTENT(IN) :: STEP_SIZE
    INTEGER, INTENT(IN) :: MAX_STEPS
    !--- B) Output
    TYPE(TimeType_t), INTENT(OUT) :: TimeResult

    !-------------------------------------------------------------------------------
    ! VARIABLES
    !-------------------------------------------------------------------------------
    LOGICAL :: IsTimeIndex, IsUserDefined, Found(5) ! True if SIF definitions found
    REAL(KIND=dp) :: TimeBias, Time, TimeEpsilon ! Biasing for every used Elmer time value/index

    !-------------------------------------------------------------------------------
    ! 1) Gathers data from SIF, applies default values
    !-------------------------------------------------------------------------------
    IsUserDefined = GetLogical( GetSolverParams(Solver), "User Defines Time", Found(1) ) ! Set to true, if old definitions are used
    IsTimeIndex = GetLogical( GetSolverParams(Solver), "Is Time Index", Found(2) ) ! If true, then the given time value is an index (Default: value)
    TimeBias = GetCReal( GetSolverParams(Solver), "NetCDF Starting Time", Found(3) ) ! Index, if IsTimeIndex is true; value otherwise
    TimeEpsilon = GetConstReal( GetSolverParams(Solver), "Epsilon Time", Found(5) ) ! The tolerance for the time value
    
    !-------------------------------------------------------------------------------
    ! 2) Collects dimension specific data
    !-------------------------------------------------------------------------------
    CALL GetDimension(NCID,T_Name, TimeResult % id, TimeResult % len)
    !> id and len set for TimeResult

    !-------------------------------------------------------------------------------
    ! 3) Sets default values, informs on effect
    !-------------------------------------------------------------------------------
    IF ( .NOT. Found(1) ) IsUserDefined = .FALSE.
    IF ( .NOT. Found(2) ) IsTimeIndex = .FALSE. ! Defaulted to False for values are more natural, otherwise set as given
    IF ( .NOT. Found(3) ) THEN
      IF ( IsTimeIndex ) THEN
        TimeBias = 1.0_dp ! Starts, by default, from index 1
      ELSE
        TimeBias = 0.0_dp ! Starts, by default, without value adjustment
      END IF
    ELSE
      ! Tell the user what is happening
      IF ( IsTimeIndex ) THEN
        WRITE (Message,'(A,F6.2)') 'Input time indices are adjusted (summed) by ', TimeBias
      ELSE
        WRITE (Message,'(A,F6.2)') 'Input time values are adjusted (summed) by ', TimeBias
      END IF
      CALL Info('GridDataMapper',Message)
    END IF
    IF ( .NOT. Found(5) ) THEN
      TimeEpsilon = 0.0_dp
    ELSE 
      WRITE (Message,'(A,F6.3)') 'Received Time Epsilon with value ', TimeEpsilon
      CALL Info('GridDataMapper',Message)
    END IF
    !> Default values set

    !-------------------------------------------------------------------------------
    ! 4) Handles transient cases (uses Elmer time, if not defined otherwise)
    !-------------------------------------------------------------------------------
    IF ( IS_TRANSIENT .AND. (.NOT. IsUserDefined ) ) THEN

      !-----------------------------------------------------------------------------
      ! A) Transient with Elmer time points
      !-----------------------------------------------------------------------------
      WRITE(Message, '(A,A)') 'Simulation is transient and user time input is ignored.',&
      ' (If own time scaling wanted, set SIF variable "User Defines Time" true and use user defined functions)'
      CALL Info('GridDataMapper', Message)

      !--- Get the time from Elmer
      Time = GetTime() 

      IF ( IsTimeIndex ) THEN  
        !---------------------------------------------------------------------------
        ! Bias time indices
        !---------------------------------------------------------------------------
        ! Indexing starts with step size in Elmer, bias it to first index

        !--- Check the bias doesn't cause over-indexing the NetCDF time range
        IF ( TimeBias .LT. 1.0_dp .OR. TimeBias .GT. TimeResult % LEN ) THEN
          WRITE (Message, '(A,F6.2,A,I5,A)') 'NetCDF Starting Time index ', TimeBias, &
                              ' does not fit within the NetCDF time index range (1,', TimeResult % LEN,')'
          CALL Fatal('GridDataMapper', Message)
        END IF

        !--- Checks the bias doesn't indirectly over-index due to iterations
        IF (MAX_STEPS .GT. ((TimeResult % LEN - (TimeBias - STEP_SIZE))/STEP_SIZE) ) THEN
          WRITE (Message, '(A,I5,A,F10.2,A,F6.2,A)') 'Defined amount of timestep intervals ', MAX_STEPS, &
                          ' is more than ', ((TimeResult % LEN - (TimeBias - STEP_SIZE))/STEP_SIZE),&
                          ', which is the maximum number of allowed size ', STEP_SIZE ,' steps on the NetCDF grid.'
          CALL Fatal('GridDataMapper',Message)
        END IF

        !--- Perform biasing (the bias is valid)
        Time = Time + (TimeBias - STEP_SIZE)
      ELSE
        !---------------------------------------------------------------------------
        ! Bias time values
        !---------------------------------------------------------------------------
        Time = Time + TimeBias
        ! NOTE: Time value is checked when it's converted to a time index with TimeValueToIndex()
        ! Does not check for eventual over-indexing!
      END IF

    ELSE 

      !-----------------------------------------------------------------------------
      ! B) For user defined and steady state
      !-----------------------------------------------------------------------------
      ! NOTE: User-defined/steady state never biased!
      Time = GetCReal( Solver % Values, "Time Point", Found(4) )

      IF ( .NOT. Found(4) ) THEN
        WRITE(Message,'(A,I3,A)') 'No time point given; specify it in the Solver Input File with name "Time Point"&
 or enable transient simulation!'
        CALL Fatal('GridDataMapper',Message)
      END IF

    END IF
    !> Time automatic and transient => Time biased
    !> Time user-defined/steady state => Time obtained directly from SIF

    !--- Converts the given time value into an appropriate time index
    IF (.NOT. IsTimeIndex) Time = TimeValueToIndex(NCID,T_Name,TimeResult % id,TimeResult % LEN,Time, TimeEpsilon)

    !--- Final check before letting through
    IF ( Time .LT. 1 .OR. Time .GT. TimeResult % LEN ) THEN
      WRITE (Message, '(A,F6.2,A,I5,A)') 'Time value ', Time, ' is out of range (1,', TimeResult % LEN, ')'
      CALL Fatal('GridDataMapper',Message)
    END IF
    !> Time values are converted to indices and checked that they are within range; eventual over-indexing not always checked
    !> I.e. Time is obtained and checked to be within range

    !-------------------------------------------------------------------------------
    ! 5) Finalizes TimeResult
    !-------------------------------------------------------------------------------

    ! Rounds to nearest integer index, if it is within epsilon range, to accomodate for insignificant differences between real values
    IF ( Found(5) .AND. (FLOOR( Time) .NE. CEILING( Time )) ) THEN
      WRITE (Message,'(A,F8.5,A,F8.5)') 'The time index with epsilon range: ', Time - FLOOR(Time), ' (+/-) ', TimeEpsilon
      CALL Info('GridDataMapper',Message)
      ! FLOOR(Time) < Time < CEILING(TIME) <=> 0 < Time - FLOOR(Time) < CEILING(TIME) - FLOOR(TIME) = 1
      ! If Epsilon makes the value exceed either limit, then the limit is where it is rounded; otherwise intact
      IF ( Time - FLOOR( Time ) .LE. TimeEpsilon ) THEN
        Time = FLOOR( Time )
      ELSE IF ( Time - FLOOR( Time ) .GE. 1.0_dp - TimeEpsilon ) THEN
        Time = CEILING( Time )
      END IF
    END IF
    TimeResult % val = Time
    TimeResult % low = FLOOR( Time )
    TimeResult % high = CEILING( Time )
    TimeResult % doInterpolation = (TimeResult % low .NE. TimeResult % high)
    IF ( TimeResult % len .EQ. 1 ) TimeResult % doInterpolation = .FALSE. ! No need to interpolate unit time
    IF ( TimeResult % low .LT. 1 ) TimeResult % low = 1
    IF ( TimeResult % high .GT. TimeResult % len ) TimeResult % high = TimeResult % len

    IF ( TimeResult % doInterpolation ) THEN
      WRITE (Message,'(A,F5.2,A)') 'Given time value ', TimeResult % val , ' is not an integer. Using time interpolation.'
      CALL Info('GridDataMapper',Message)
    ELSE
      WRITE (Message,'(A,F5.1,A)') 'Given time value ', TimeResult % val, ' is an integer. No time interpolation used.'
      CALL Info('GridDataMapper',Message)
    END IF

    !-------------------------------------------------------------------------------
  
  END SUBROUTINE InitTime


  !---------------- InitInterpolation() ---------------
  !--- Initializes the information needed for interpolation
  SUBROUTINE InitInterpolation( Solver, DIM, Grids, TimeInterpolationMethod, InterpMultiplier, InterpBias )
    USE DefUtils, ONLY: GetSolverParams, MAX_NAME_LEN, GetConstReal, GetString, GetCReal, CoordinateSystemDimension
    USE Messages
    IMPLICIT NONE

    !-------------------------------------------------------------------------------
    ! ARGUMENTS
    !-------------------------------------------------------------------------------
    TYPE(Solver_t), INTENT(IN) :: Solver
    INTEGER, INTENT(IN) :: DIM
    TYPE(UniformGrid_t), INTENT(INOUT) :: Grids(:)
    REAL(KIND=dp), INTENT(OUT) :: InterpMultiplier, InterpBias
    CHARACTER(len=MAX_NAME_LEN), INTENT(OUT) :: TimeInterpolationMethod

    !-------------------------------------------------------------------------------
    ! VARIABLES
    !-------------------------------------------------------------------------------
    LOGICAL :: tmpBool
    REAL(KIND=dp), ALLOCATABLE :: eps(:)
    INTEGER :: loop, alloc_stat

    !-------------------------------------------------------------------------------
    ! 1) Adds Epsilon values to Grids (if applicable)
    !-------------------------------------------------------------------------------
    IF ( Grids(1) % COORD_COUNT .GT. 0 ) THEN

      !--- Allocation
      ALLOCATE ( eps( Grids(1) % COORD_COUNT ), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF
      eps = 0

      ! Epsilons are the relative tolerances for the amount 
      ! the Elmer grid point misses the bounds of the NetCDF bounding box
      DO loop = 1,size(eps),1

        ! NOTE: There won't be more than 9 epsilons, because there can't be enough Elmer coordinates for that!
        WRITE(Message, '(A,I1)') 'Epsilon ', loop
        eps(loop) = GetConstReal(GetSolverParams(Solver), Message, tmpBool)
        IF ( .NOT. tmpBool ) THEN
           eps(loop) = 0
           WRITE(Message,'(A,I1,A,A)') 'Variable "Epsilon ', loop ,'" not given in Solver Input File. ',&
                                   'Using zero tolerance for NetCDF grid mismatches with Elmer mesh.'
           CALL Warn('GridDataMapper', Message)
        END IF
      END DO

      ! Sets the appropriate values to the Grids 
      Grids(1) % Eps(:) = eps(:) * Grids(1) % dx(:)
      Grids(2) % Eps(:) = eps(:) * Grids(2) % dx(:)
    END IF

    !-------------------------------------------------------------------------------
    ! 2) Sets the time interpolation method 
    !-------------------------------------------------------------------------------
    TimeInterpolationMethod = GetString( GetSolverParams(Solver), "Time Interpolation Method", tmpBool )
    IF ( .NOT. tmpBool ) THEN
      CALL Warn('GridDataMapper', 'SIF variable "Time Interpolation Method" not specified, using default settings.')
      TimeInterpolationMethod = ''
    ELSE
      WRITE (Message,'(A,A)' ) 'Received Time Interpolation Method SIF variable value: ', TimeInterpolationMethod
      CALL Info('GridDataMapper', Message)
    END IF
 
    !-------------------------------------------------------------------------------
    ! 3) Sets the final adjustments for the interpolation results (to allow relative/absolute time, f.ex.)
    !-------------------------------------------------------------------------------
    
    ! If absolute time is relative, the value of the time function is added to the interpolated location; otherwise it is multiplied with it
    InterpMultiplier = GetCReal( GetSolverParams(Solver), "Interpolation Multiplier", tmpBool ) ! Multiplies the final interpolation result by given number (relative time)
    IF ( .NOT. tmpBool ) InterpMultiplier = 1.0_dp ! Defaulted to 1, so it doesn't modify the result
  
    InterpBias = GetCReal( GetSolverParams(Solver), "Interpolation Bias", tmpBool ) ! Adds the bias to the final interpolation result (absolute time)
    IF ( .NOT. tmpBool ) InterpBias = 0.0_dp ! Defaulted to 0, so there is no bias on time
  
    !-------------------------------------------------------------------------------

  END SUBROUTINE InitInterpolation

!------------------------------------------------------------------------------
END SUBROUTINE GridDataMapper
!------------------------------------------------------------------------------

