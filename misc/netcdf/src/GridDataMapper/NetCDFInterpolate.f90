!------------------------------------------------------------------------------
! Peter Råback, Vili Forsell
! Created: 13.6.2011
! Last Modified: 13.7.2011
!------------------------------------------------------------------------------
! This module contains functions for
! - interpolating NetCDF data for an Elmer grid point (incl. coordinate transformation); Interpolate()
! - coordinate transformations for an Elmer grid point
!    o none by default (x,y,z) -> (x,y,z)
!    o cs2cs interface with parameters defined in Solver Input File
!    o cartesian-to-cylindrical transformation (x,y,z) -> (phi,r,z)
! - scaling Elmer mesh points to fit NetCDF data
!------------------------------------------------------------------------------
MODULE NetCDFInterpolate
  USE DefUtils, ONLY: dp, MAX_NAME_LEN
  USE NetCDFGeneralUtils, ONLY: GetFromNetCDF
  USE Messages
  IMPLICIT NONE

  INTERFACE
    !--- For connecting the C code, which accesses the cs2cs
    SUBROUTINE cs2cs_transform( coord, hasZ, isRad, elmer_proj, netcdf_proj, res ) BIND(c)
      USE iso_c_binding

      !--- Input parameters
      REAL(C_DOUBLE) :: coord(3)
      INTEGER(C_INT), VALUE :: hasZ, isRad
      CHARACTER(KIND=C_CHAR) :: elmer_proj(*), netcdf_proj(*)

      !--- Output parameters
      REAL(C_DOUBLE) :: res(3)
      
    END SUBROUTINE cs2cs_transform
  END INTERFACE

  LOGICAL :: DEBUG_INTERP = .FALSE.
  PRIVATE :: GetSolutionInStencil, CoordinateTransformation, ScaleMeshPoint, ChooseInterpolation, GetScalar

  CONTAINS

    !------------------ LinearInterpolation() ---------------------
    !--- Performs linear interpolation
    !--------------------------------------------------------------
    FUNCTION LinearInterpolation(x,u1,u2) RESULT(y)
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: u1(2),u2(2),x
      REAL(KIND=dp) :: y

      y = (((u2(2) - u1(2))/(u2(1) - u1(1)))*(x-u1(1)))+u1(2)
    END FUNCTION LinearInterpolation


    !------------------- BilinearInterpolation() -------------------
    !--- Performs bilinear interpolation on a stencil (2x2 matrix of corner values)
    !---  with given weights (a 2 dimensional vector)
    !---------------------------------------------------------------
    FUNCTION BiLinearInterpolation(stencil,weights) RESULT(y)
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: stencil(2,2), weights(2)
      REAL(KIND=dp) :: y

      y = stencil(1,1)*(1-weights(1))*(1-weights(2)) + &
          stencil(2,1)*weights(1)*(1-weights(2)) + &
          stencil(1,2)*(1-weights(1))*weights(2) + &
          stencil(2,2)*weights(1)*weights(2)

    END FUNCTION BiLinearInterpolation

    !------------------- TrilinearInterpolation() -------------------
    !--- Performs trilinear interpolation on a stencil (2x2x2 matrix of corner values)
    !---  with given weights (a 3 dimensional vector)
    !---------------------------------------------------------------
    FUNCTION TriLinearInterpolation(stencil,weights) RESULT(y)
      USE DefUtils
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: stencil(2,2,2), weights(3)
      REAL(KIND=dp) :: val1,val2,y

      val1 = BiLinearInterpolation(stencil(:,:,1),weights(1:2))
      val2 = BiLinearInterpolation(stencil(:,:,2),weights(1:2))
      y = val1*(1-weights(3)) + val2*weights(3)

    END FUNCTION TriLinearInterpolation

    !------------------ Interpolate() -----------------------------
    !--- Takes and interpolates one Elmer grid point to match NetCDF data; includes coordinate transformation
    !--- ASSUMES INPUT DIMENSIONS AGREE
    !--------------------------------------------------------------
    FUNCTION Interpolate(SOLVER,NCID,VAR_ID,DIM_LENS,GRID,&
            TIME,TIME_IND,interp_val,COORD_SYSTEM,X) RESULT( success )
      USE DefUtils, ONLY: Solver_t
      USE NetCDFGridUtils, ONLY: UniformGrid_t
      USE NetCDFGeneralUtils, ONLY: TimeType_t
      IMPLICIT NONE

      !------------------------------------------------------------------------------
      ! ARGUMENTS
      !------------------------------------------------------------------------------
      TYPE(Solver_t), INTENT(IN) :: SOLVER
      TYPE(UniformGrid_t), INTENT(IN) :: GRID
      TYPE(TimeType_t), INTENT(IN) :: TIME
      INTEGER, INTENT(IN) :: NCID,DIM_LENS(:),TIME_IND,VAR_ID
      CHARACTER(len = *), INTENT(IN) :: COORD_SYSTEM
      REAL(KIND=dp), INTENT(INOUT) :: interp_val ! Final Elmer point and interpolated value 
      REAL(KIND=dp), OPTIONAL, INTENT(IN) :: X(:)
      LOGICAL :: success ! Return value
      
      !------------------------------------------------------------------------------
      ! VARIABLES
      !------------------------------------------------------------------------------
      INTEGER :: alloc_stat, i
      INTEGER, ALLOCATABLE :: ind(:)
      REAL(KIND=dp), ALLOCATABLE :: Xf(:)

      IF ( PRESENT(X) .AND. GRID % COORD_COUNT .GT. 0 ) THEN
        !------------------------------------------------------------------------------
        ! Uses Elmer coordinates
        !------------------------------------------------------------------------------
        ALLOCATE (ind(GRID % DIMS), Xf(size(X)), STAT = alloc_stat)
        IF ( alloc_stat .NE. 0 ) THEN
          CALL Fatal('GridDataMapper','Interpolation vectors memory allocation failed')
        END IF
    
        !--- 1) Coordinate mapping from Elmer (x,y) to the one used by NetCDF
        Xf = CoordinateTransformation( SOLVER, X, COORD_SYSTEM )
  
        !--- 2) Scaling, if applicable
        ! NOTE! By default the SCALE consists of 1's, and MOVE 0's; hence, nothing happens
        !       without user specifically specifying so
        Xf = GRID % SCALE(:)*Xf + GRID % MOVE(:) ! Scales the mesh point within the NetCDF grid
  
        !--- 3) Index estimation
        ! Find the (i,j) indices [1,...,max] 
        ! Calculates the normalized difference vector; 
        ! i.e. the distance/indices to Elmer grid point x from the leftmost points of the NetCDF bounding box
        ind(1:GRID % COORD_COUNT) = CEILING( ( Xf(:) - GRID % X0(:) ) / GRID % DX(:) ) 
        ind((GRID % COORD_COUNT+1):GRID % DIMS) = GRID % CONST_VALS
 
        ! This could be done better, one could apply extrapolation 
        ! with a narrow layer.
        DO i = 1,size(Xf,1),1 ! NOTE: Does not modify the constant dimensions
    
          ! Checks that the estimated index is within the bounding box
          IF( ind(i) < 1 .OR. ind(i) >= GRID % NMAX(i) ) THEN
    
            ! If it's smaller than the leftmost index, but within tolerance (Eps), set it to lower bound; and vice versa
            IF( Xf(i) <= GRID % X0(i) .AND. Xf(i) >= GRID % X0(i) - GRID % EPS(i) ) THEN
              ind(i) = 1
            ELSE IF( Xf(i) >= GRID % X1(i) .AND. Xf(i) <= GRID % X1(i) + GRID % EPS(i) ) THEN
              ind(i) = GRID % NMAX(i)
            ELSE ! The index is too far to be salvaged
              WRITE (Message, '(A,F14.3,A,(F14.3),A,F14.3,A,I3,A,F14.3,A,F14.6,A)') &
                'Adjusted Elmer value is out of NetCDF bounds: ',&
               GRID % X0(i), ' <= ',Xf(i), ' <= ', GRID % X1(i), &
              ' over dimension ',i,' and originating from ', X(i),' despite error tolerance epsilon: ', GRID % EPS(i), '.'
              CALL Warn( 'GridDataMapper',Message)
              success = .FALSE.
              RETURN
            END IF
          END IF
        END DO

        !--- 4) Choose and perform interpolation 
        interp_val = ChooseInterpolation(NCID,VAR_ID,GRID,TIME,Xf,IND,TIME_IND,DIM_LENS)

      ELSE
        !------------------------------------------------------------------------------
        ! No Elmer coordinates (only constants/time)
        !------------------------------------------------------------------------------
        CALL GetScalar(NCID,VAR_ID,GRID,TIME,TIME_IND,DIM_LENS,interp_val)

      END IF  
 
      success = .TRUE.
  
    END FUNCTION Interpolate

    !----------------- ChooseInterpolation() ----------------------
    !--- Chooses the appropriate weighting, stencil and interpolation method
    FUNCTION ChooseInterpolation(NCID,VAR_ID,GRID,TIME,X_VAL,X_IND,TIME_IND,DIM_LENS) RESULT(interp_val)
    !--------------------------------------------------------------
      USE NetCDFGridUtils, ONLY: UniformGrid_t
      USE NetCDFGeneralUtils, ONLY: TimeType_t
      IMPLICIT NONE

      !------------------------------------------------------------------------------
      ! ARGUMENTS
      !------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NCID, VAR_ID, X_IND(:), DIM_LENS(:), TIME_IND
      TYPE(UniformGrid_t), INTENT(IN) :: GRID
      TYPE(TimeType_t), INTENT(IN) :: TIME
      REAL(KIND=dp), INTENT(IN) :: X_VAL(:)
      REAL(KIND=dp) :: interp_val ! The output

      !------------------------------------------------------------------------------
      ! VARIABLES
      !------------------------------------------------------------------------------
      INTEGER :: alloc_stat, loop
      INTEGER :: actual_size, & ! The actual amount of stencil values
                 actual_dims ! How many dimensions remain in use
      INTEGER :: last_loc ! Latest location on the array
      LOGICAL :: changed ! True, if the stencil has been changed
      INTEGER, ALLOCATABLE :: sizes(:) ! Possibly differing sizes of the stencils
      REAL(KIND=dp), ALLOCATABLE :: xi(:), weights(:), stencil_data(:)
      REAL(KIND=dp) :: u1(2),u2(2) ! For 1D (linear) interpolation
      REAL(KIND=dp), ALLOCATABLE :: stencilLine(:),stencilSqr(:,:),stencilCube(:,:,:) ! Adjusted if seems to over-index

      !------------------------------------------------------------------------------
      ! INITIALIZATIONS
      !------------------------------------------------------------------------------
      changed = .FALSE.
      actual_size = 1
      actual_dims = GRID % DIMS

      !--- Finds the locations which are on the edge of over-indexing, and limits them if necessary
      ALLOCATE ( sizes(GRID % DIMS), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out!')
      END IF
      sizes = 2

      DO loop = 1,size(X_IND,1)
        IF ( X_IND(loop) .EQ. DIM_LENS(loop) ) THEN ! Ignore last dimensions
          sizes(loop) = 1
          actual_dims = actual_dims - 1
          changed = .TRUE.
        END IF
        actual_size = actual_size * sizes(loop) ! Calculates the actual size of the stencil
      END DO
      
      !--- The sizes will change: need to allocate a temporary array for reshaping into a size suitable for interpolation
      IF ( changed ) THEN
        ALLOCATE ( stencil_data(actual_size), STAT = alloc_stat )
        IF ( alloc_stat .NE. 0 ) THEN
          CALL Fatal('GridDataMapper','Memory ran out!')
        END IF
      END IF

      !--- Allocates the weights and such
      ALLOCATE ( xi(actual_dims), weights(actual_dims), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out!')
      END IF

      !------------------------------------------------------------------------------
      ! A) Calculating the weights for each used dimension
      !------------------------------------------------------------------------------
      last_loc = 1 ! The latest handled location of the modified vectors
      DO loop = 1,size(sizes,1),1
        ! If the dimension is used, its size is larger than one (otherwise doesn't affect actual_size)
        IF ( sizes(loop) .GT. 1 ) THEN
          ! The value of the estimated NetCDF grid point
          xi(last_loc) = GRID % X0(loop) + (X_IND(loop)-1) * GRID % DX(loop)
  
          ! Interpolation weights, which are the normalized differences of the estimation from the adjacent grid point
          ! Can be negative if ceil for indices brings the value of xi higher than x
          !----------- Assume xi > x ------
          !  x0 + (ind-1)dx > x, where dx > 0
          ! <=> (ind-1)dx > x-x0
          ! <=> ind-1 > (x-x0)/dx
          ! <=> ceil((x-x0)/dx) > ((x-x0)/dx) + 1
          ! o Known  (x-x0)/dx  <= ceil((x-x0)/dx) < (x-x0)/dx+1 
          ! => Contradicts; Ergo, range ok.
          !--------------------------------
          ! p values should be within [0,1]
          ! 0 exactly when x = xi, 1 when (x-x0)/dx = ceil((x-x0)/dx) = ind
          weights(last_loc) = (X_VAL(loop)-xi(last_loc))/GRID % DX(loop)
          last_loc = last_loc + 1
        END IF
      END DO
  
      !------------------------------------------------------------------------------
      ! B) Obtaining the stencil values 
      !------------------------------------------------------------------------------
      !--- (Must be in original dimensions for the data acquisition to, f.ex., match the contents of X_IND)
      SELECT CASE ( GRID % DIMS )
        CASE (1) !-- 1D 
          ALLOCATE ( stencilLine(sizes(1)), STAT = alloc_stat )
          IF ( alloc_stat .NE. 0 ) THEN
            CALL Fatal('GridDataMapper','Memory ran out!')
          END IF

          CALL GetSolutionInStencil(NCID,VAR_ID,X_IND,TIME_IND,TIME,DIM_LENS,GRID % ACCESS_PERM,stencilLine = stencilLine)
          IF ( changed ) stencil_data = RESHAPE(stencilLine,SHAPE(stencil_data))

        CASE (2) !-- 2D
          ALLOCATE ( stencilSqr(sizes(1),sizes(2)), STAT = alloc_stat )
          IF ( alloc_stat .NE. 0 ) THEN
            CALL Fatal('GridDataMapper','Memory ran out!')
          END IF

          ! get data on stencil size(stencil)=(2,2), ind -vector describes the lower left corner
          CALL GetSolutionInStencil(NCID,VAR_ID,X_IND,TIME_IND,TIME,DIM_LENS,GRID % ACCESS_PERM,stencilSqr = stencilSqr)
          IF ( changed ) stencil_data = RESHAPE(stencilSqr,SHAPE(stencil_data))

        CASE (3) !-- 3D

          ALLOCATE ( stencilCube(sizes(1),sizes(2),sizes(3)), STAT = alloc_stat )
          IF ( alloc_stat .NE. 0 ) THEN
            CALL Fatal('GridDataMapper','Memory ran out!')
          END IF

          CALL GetSolutionInStencil(NCID,VAR_ID,X_IND,TIME_IND,TIME,DIM_LENS,GRID % ACCESS_PERM,stencilCube = stencilCube)
          IF ( changed ) stencil_data = RESHAPE(stencilCube,SHAPE(stencil_data))

        CASE DEFAULT !-- Error
          CALL Fatal('GridDataMapper','Cannot handle more than three variable dimensions!')
      END SELECT

      !------------------------------------------------------------------------------
      ! C) Interpolation
      !------------------------------------------------------------------------------
      !--- Allocates the interpolation arrays, if necessary
      SELECT CASE ( actual_size )
        CASE (1) ! Scalar

!          PRINT *, 'Scalar: ', stencil_data
          !-- stencil_data contains it all
          interp_val = stencil_data(1) ! NOTE: Returning a scalar follows only if changed is .TRUE.

        CASE (2) ! Line

          !--- Adjusts data/weights to proper size
          IF ( changed ) THEN
!            PRINT *, 'Line: ', stencil_data
            ALLOCATE ( stencilLine(2), STAT = alloc_stat )
            IF ( alloc_stat .NE. 0 ) THEN
              CALL Fatal('GridDataMapper','Memory ran out!')
            END IF
            stencilLine = stencil_data
          END IF

          !--- Linear interpolation
          !--- Note: X_IND is the point in the lower left corner; so, in this case the leftmost point of the linear line
          u1(1) = X_IND(1) ! Linear location from scalar x coord (integer),
          u1(2) = stencilLine(1) ! interpolation for the corresponding value
          u2(1) = X_IND(1) + 1 ! Same for the adjacent point before interpolation
          u2(2) = stencilLine(2)
          interp_val = LinearInterpolation((X_IND(1) + weights(1)),u1,u2)

        CASE (4) ! Square

          !--- Adjusts data/weights to proper size
          IF ( changed ) THEN
!            PRINT *, 'Square: ', stencil_data
            ALLOCATE ( stencilSqr(2,2), STAT = alloc_stat )
            IF ( alloc_stat .NE. 0 ) THEN
              CALL Fatal('GridDataMapper','Memory ran out!')
            END IF
            stencilSqr = RESHAPE(stencil_data,SHAPE(stencilSqr))
          END IF

          !--- Bilinear interpolation
          interp_val = BiLinearInterpolation(stencilSqr,weights)

        CASE (8) ! Cube

          !--- Adjusts data/weights to proper size
          IF ( changed ) THEN
!            PRINT *, 'Cube: ', stencil_data
            ALLOCATE ( stencilCube(2,2,2), STAT = alloc_stat )
            IF ( alloc_stat .NE. 0 ) THEN
              CALL Fatal('GridDataMapper','Memory ran out!')
            END IF
            stencilCube = RESHAPE(stencil_data,SHAPE(stencilCube))
          END IF

          !--- Trilinear interpolation
          interp_val = TriLinearInterpolation(stencilCube,weights)

        CASE DEFAULT ! Error
          WRITE(Message,'(A,I5,A)') 'Cannot interpolate a stencil of size ', actual_size, '.'
          CALL Fatal('GridDataMapper',Message)

      END SELECT

    END FUNCTION ChooseInterpolation

    !----------------- CoordinateTransformation() -----------------
    !--- Transforms input coordinates into the given coordinate system
    FUNCTION CoordinateTransformation( SOLVER, INPUT, COORD_SYSTEM ) RESULT( output )
    !--------------------------------------------------------------
      USE DefUtils, ONLY: dp, MAX_NAME_LEN, Solver_t, GetSolverParams, GetString, GetLogical
      USE iso_c_binding, ONLY: C_NULL_CHAR
      USE Messages
      IMPLICIT NONE
 
      !--- Input arguments
      TYPE(Solver_t), INTENT(IN) :: SOLVER
      CHARACTER(*), INTENT(IN) :: COORD_SYSTEM ! Some coordinate
      REAL(KIND=dp), INTENT(IN) :: INPUT(:) ! The input coordinates

      !--- Return value
      REAL(KIND=dp), ALLOCATABLE :: output(:) ! The output coordinates

      !--- Others
      REAL(KIND=dp) :: coord(3), res(3)
      INTEGER :: alloc_stat, hasZcoord
      LOGICAL :: found
      CHARACTER(len=MAX_NAME_LEN) :: elmer, netcdf

      !--- DEBUG printout
!      WRITE (*,*) 'Input ', input

      !--- Initializations
      coord = 0
      res = 0
      ALLOCATE ( output(size(INPUT)), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Coordinate transformation memory allocation failed')
      END IF

      !--- Selects the coordinate system transformation
      SELECT CASE (COORD_SYSTEM)

        !--- CS2CS Coordinate transformation
        CASE ('cs2cs')
!          CALL Info('GridDataMapper', 'Applies cs2cs coordinate transformation between the Elmer and NetCDF values!')

          !-- Gathers the coordinate information
          hasZcoord = 0
          IF ( size(INPUT,1) .GE. 3 ) THEN
            hasZcoord = 1
            coord(1:3) = INPUT(1:3)
          ELSE
            coord(1:2) = INPUT(1:2)
          END IF

          !--- DEBUG printout
!          WRITE (*,*) 'Coordinates ', coord

          !--- Picks up the constant data from the Elmer Solver parameters
          elmer = GetString(GetSolverParams(SOLVER),"CS2CS Elmer Projection", found)
          IF ( .NOT. found ) THEN
            CALL Fatal('GridDataMapper',&
  'CS2CS Transformation did not find Elmer projection information "CS2CS Elmer Projection" from the Solver Input File')
          END IF

          netcdf = GetString(GetSolverParams(SOLVER), "CS2CS NetCDF Projection", found)
          IF ( .NOT. found ) THEN
            CALL Fatal('GridDataMapper',&
  'CS2CS Transformation did not find NetCDF projection information "CS2CS NetCDF Projection" from the Solver Input File')
          END IF

          ! //C_NULL_CHAR's terminate the given strings with nulls to enable C compatibility
          IF ( GetLogical(GetSolverParams(SOLVER), "CS2CS Is Input Radians", found) .AND. found ) THEN
            CALL cs2cs_transform( coord, hasZcoord, 1, elmer//C_NULL_CHAR, netcdf//C_NULL_CHAR, res) ! True; is in radians
          ELSE
            CALL cs2cs_transform( coord, hasZcoord, 0, elmer//C_NULL_CHAR, netcdf//C_NULL_CHAR, res) ! False; is in degrees by default
          END IF

          !--- Sends the result as output
          IF ( size(output,1) .GE. 3 ) THEN
            output(1:3) = res(1:3)
          ELSE
            output(1:2) = res(1:2)
          END IF

          !--- DEBUG printout
!          WRITE (*,*) 'Result ', res, ' to out ', output

        CASE ('cylindrical')
!          CALL Info('GridDataMapper','Applies cylindrical coordinate transformation to cartesian coordinates!')

          !--- Transforms 3D Elmer grid points to cylindrical coordinates (phi,r,z)
          IF ( size(INPUT,1) .GE. 3) THEN
            output(1) = atan2( INPUT(2), INPUT(1) ) ! phi angle value
            output(2) = sqrt(INPUT(1)**2 + INPUT(2)**2) ! radius from the center of cylinder
            output(3) = INPUT(3) ! Height from zero level
          ELSE
            CALL Fatal('GridDataMapper','Cylindrical coordinate transformation requires at least three dimensional Elmer grid.')
          END IF

        !--- In default case, no coordinate transformation is applied
        CASE DEFAULT
!          WRITE (Message,'(A,A15,A)') 'No coordinate transformation applied: Unknown coordinate system "',&
!                coord_system, '". Check Solver Input File and the variable "Coordinate System"'
!          CALL Warn('GridDataMapper', Message)
          output(:) = INPUT(:)
      END SELECT

    END FUNCTION
  
    !------------------ ScaleMeshPoint() ------------------------
    !--- Takes an Elmer mesh point and moves and scales it within the NetCDF grid
    !--- Assumed that Elmer mesh and NetCDF grid should be 1:1, but aren't still completely matched
    !--- NOTE: Can be optimized by calculating move(:) and scales(:) before interpolation (constant over a mesh/grid combo)
    FUNCTION ScaleMeshPoint(X,X0,X1,X0E,X1E) RESULT( Xf )
    !------------------------------------------------------------
      USE Messages
      USE DefUtils, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp), INTENT(IN) :: X(:), &! The input Elmer point
                           X0(:), X1(:), & ! The limiting values (points) of NetCDF grid
                           X0E(:), X1E(:) ! The limiting values (points) of Elmer bounding box
      REAL(KIND=dp), ALLOCATABLE :: Xf(:), & ! Scaled value; the output
                       move(:), & ! Moves the Elmer min value to the NetCDF min value
                       scales(:) ! Scales the grids to same value range (NetCDF constant, Elmer varies)
      INTEGER :: alloc_stat ! For allocation

      !--- Initial checks and allocations

      ! All sizes are the same
      IF ( .NOT. ( (size(X) .EQ. size(X0))   .AND. (size(X0) .EQ. size(X1)) .AND. &
                   (size(X1) .EQ. size(X0E)) .AND. (size(X0E) .EQ. size(X1E)) ) ) THEN
        CALL Fatal( 'GridDataMapper', 'Scaling input point sizes do not match!')
      END IF
      ALLOCATE ( Xf(size(X)), move(size(X)), scales(size(X)), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out during scaling')
      END IF

      Xf = 0
      move = 0
      scales = 0
      !--- Calculates the modifications

      ! First the scaling to same size (Eq. a( X1E(1)-X0E(1) ) = (X1(1)-X0(1)) ; ranges over a dimension are same. Solved for a, 1 if equal)
      scales(:) = (X1(:)-X0(:))/(X1E(:)-X0E(:)) ! Note: "/" and "*" elementwise operations for arrays in Fortran

      ! Second the vector to reach X0 from the scaled X0E (wherever it is)
      move(:) = X0(:) - scales(:)*X0E(:) ! zero, if equal

      !--- Applies the modification
      Xf(:) = scales(:)*X(:) + move(:)

    END FUNCTION ScaleMeshPoint

    !------------------ GetScalar() -------------------------------
    !--- Gets a scalar value from NetCDF
    SUBROUTINE GetScalar( NCID, VAR_ID, GRID, TIME, TIME_IND, DIM_LENS, scalar )
    !--------------------------------------------------------------
      USE NetCDFGeneralUtils, ONLY: TimeType_t, GetFromNetCDF
      USE NetCDFGridUtils, ONLY: UniformGrid_t
      USE DefUtils, ONLY: dp
      IMPLICIT NONE

      !------------------------------------------------------------------------------
      ! ARGUMENTS
      !------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NCID, VAR_ID, TIME_IND, DIM_LENS(:)
      TYPE(UniformGrid_t), INTENT(IN) :: GRID
      TYPE(TimeType_t), INTENT(IN) :: TIME
      REAL(KIND=dp), INTENT(OUT) :: scalar

      !------------------------------------------------------------------------------
      ! VARIABLES
      !------------------------------------------------------------------------------
      REAL(KIND=dp), ALLOCATABLE :: data(:)
      INTEGER, ALLOCATABLE :: singleton(:)
      INTEGER :: alloc_stat

      !------------------------------------------------------------------------------
      ! Basic checks and data access with the constant values defined in the Grid
      !------------------------------------------------------------------------------
      ALLOCATE ( singleton(size(DIM_LENS)), STAT = alloc_stat )
      IF ( alloc_stat .NE. 0 ) THEN
        CALL Fatal('GridDataMapper','Memory ran out')
      END IF
      singleton = 1

      IF ( .NOT. GRID % IS_DEF ) CALL Fatal('GridDataMapper',&
                     'GetScalar requires a defined grid')
      IF ( GRID % DIMS .LE. GRID % COORD_COUNT ) CALL Fatal('GridDataMapper',&
                     'GetScalar requires constant values for accessing NetCDF')

      IF ( GetFromNetCDF(NCID,VAR_ID,GRID % CONST_VALS(GRID % ACCESS_PERM(:)),&
               TIME % low,TIME,DIM_LENS(GRID % ACCESS_PERM),data,singleton) ) THEN
        scalar = data(1)
      END IF

    END SUBROUTINE GetScalar
 
    !------------------ GetSolutionStencil() ----------------------
    !--- Gets a square matrix starting from the lower left index, the size is defined by input matrix stencil 
    SUBROUTINE GetSolutionInStencil( NCID,VAR_ID,X_IND,TIME_IND,TIME,DIM_LENS,ACC_PERM,stencilLine,stencilSqr,stencilCube )
    !--------------------------------------------------------------
      USE NetCDFGeneralUtils, ONLY: TimeType_t, GetFromNetCDF
      IMPLICIT NONE

      !--- Arguments
      INTEGER, INTENT(IN) :: NCID,X_IND(:),TIME_IND,DIM_LENS(:),VAR_ID,ACC_PERM(:)
      TYPE(TimeType_t), INTENT(IN) :: TIME
      REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: stencilLine(:),stencilSqr(:,:),stencilCube(:,:,:)

      !--- Variables
      INTEGER :: i,j,alloc_stat
      REAL(KIND=dp), ALLOCATABLE :: data(:)
      CHARACTER(len = 50) :: answ_format

!      WRITE (*,*) 'Stencil ', stencil(:,1), ' ; ', stencil(:,2) ,' X: ', X,' Y: ', Y   

      !--- Checks that the input exists and is unique
      IF ( PRESENT(stencilLine) .AND. (.NOT. (PRESENT(stencilSqr) .OR. PRESENT(stencilCube))) ) THEN !--- 1D

        ! Queries the stencil from NetCDF with associated error checks
        IF ( GetFromNetCDF(NCID,VAR_ID,X_IND(ACC_PERM(:)),&
            TIME_IND,TIME,DIM_LENS(ACC_PERM(:)),data,SHAPE(stencilLine)) ) THEN
   
          stencilLine = RESHAPE(data,SHAPE(stencilLine))
          IF ( DEBUG_INTERP ) THEN
            !------ Debug printouts -------------------------
            WRITE (*,*) 'STENCIL LINE:'
            WRITE (answ_format, *) '(', size(stencilLine,1),'(F10.4))'
            WRITE (*,answ_format) stencilLine(:)
            !------------------------------------------------
          END IF
        END IF

      ELSE IF ( PRESENT(stencilSqr) .AND. (.NOT. (PRESENT(stencilLine) .OR. PRESENT(stencilCube))) ) THEN !--- 2D

        ! Queries the stencil from NetCDF with associated error checks
        IF ( GetFromNetCDF(NCID,VAR_ID,X_IND(ACC_PERM(:)),&
TIME_IND,TIME,DIM_LENS(ACC_PERM(:)),data,SHAPE(stencilSqr)) ) THEN
  
          stencilSqr = RESHAPE(data,SHAPE(stencilSqr))
          IF ( DEBUG_INTERP ) THEN
            !------ Debug printouts -------------------------
            WRITE (*,*) 'STENCIL SQUARE:'
            DO i = 1,size(stencilSqr,2)
              WRITE (answ_format, *) '(', size(stencilSqr,1),'(F10.4))'
              WRITE (*,answ_format) stencilSqr(:,i)
            END DO
            !------------------------------------------------
          END IF
        END IF

      ELSE IF ( PRESENT(stencilCube) .AND. (.NOT. (PRESENT(stencilSqr) .OR. PRESENT(stencilLine))) ) THEN !--- 3D

        ! Queries the stencil from NetCDF with associated error checks
        IF ( GetFromNetCDF(NCID,VAR_ID,X_IND(ACC_PERM(:)),&
TIME_IND,TIME,DIM_LENS(ACC_PERM(:)),data,SHAPE(stencilCube)) ) THEN
  
         stencilCube = RESHAPE(data,SHAPE(stencilCube))
          IF ( DEBUG_INTERP ) THEN
            !------ Debug printouts -------------------------
            WRITE (*,*) 'STENCIL CUBE:'
            DO j = 1,size(stencilCube,3)
              DO i = 1,size(stencilCube,2)
                WRITE (answ_format, *) '(', size(stencilCube,1),'(F10.4))'
                WRITE (*,answ_format) stencilCube(:,i,j)
              END DO
              WRITE(*,'(A)') '----'
            END DO
            !------------------------------------------------
          END IF
        END IF
      ELSE
        CALL Fatal('GridDataMapper','Multiple, or no, stencils given for GetSolutionInStencil()!')
      END IF
      
    END SUBROUTINE GetSolutionInStencil


END MODULE NetCDFInterpolate
