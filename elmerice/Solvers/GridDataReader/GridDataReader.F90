!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
! *****************************************************************************
!
! ******************************************************************************
! *
! *  Authors: Peter R�back, Vili Forsell, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 7.6.2011
! *  Modification Date: 2.12.2011
! *
! *  Modified by Rupert Gladstone on 17.12.2012 to allow for 2D and non-
! *  uniform coordinate variables in the netcdf file (netcdf cells should be
! *  rectangular)
! *
! *  18.3.2013 / Peter R�back
! *  Improved possibility to deal with missing or errorness values.
! *
! ****************************************************************************


MODULE NetCDFInterface
  USE DefUtils
  USE NetCDF
  IMPLICIT NONE

  ! These global variables hide the NetCDF specific Id:s from the calling level.
  !----------------------------------------------------------------------------
  INTEGER :: FileId = -1
  INTEGER :: DimIds(4),CoordVarIds(4)
  INTEGER :: VarId
  INTEGER :: NetCDFStatus
  LOGICAL :: Debug = .FALSE.

  SAVE FileId, DimIds, CoordVarIds, VarId, Debug
  PRIVATE Fileid, DimIds, CoordVarIds, VarId, Debug, NetCDFStatus

  INTERFACE  NetCDFCoordVar
     MODULE PROCEDURE NetCDFCoordVar1D
     MODULE PROCEDURE NetCDFCoordVar2D
     MODULE PROCEDURE NetCDFCoordVar3D
  END INTERFACE

  CONTAINS

    !------------------------------------------------------------------------------------
    !> Gathers and initializes all the necessary NetCDF information for picking variables.
    !------------------------------------------------------------------------------------
    SUBROUTINE NetCDFInit(Params, NetDim, DimSize, CoordVarNDims, TimeSize, x0, dx, t0, dt, UniformCoords )
      !--------------------------------------------------
      TYPE(ValueList_t), POINTER, INTENT(IN) :: Params
      INTEGER, INTENT(OUT) :: NetDim                  ! Dimension of the NetCDF file
      INTEGER, INTENT(OUT) :: DimSize(:)              ! Lengths for all space dimensions
      INTEGER, INTENT(OUT) :: CoordVarNDims(:)        ! Number of dimensions of coordinate variables (time should always be one dimensional)
                                                      ! (assumptions will be made about which dimensions apply to which coord vars)
      INTEGER, INTENT(OUT) :: TimeSize                ! Length for time dimension
      REAL(KIND=dp), INTENT(OUT) :: x0(:)             ! Minimum coordinate values for the NetCDF grid
      REAL(KIND=dp), INTENT(OUT) :: dx(:)             ! Grid resolution for the NetCDF grid
      REAL(KIND=dp), INTENT(OUT) :: t0                ! Minimum time value
      REAL(KIND=dp), INTENT(OUT) :: dt                ! Time resolution
      LOGICAL, INTENT(OUT)       :: UniformCoords     ! Will only be true if all coordinate variables appear to be uniform
      !------------------------------------------------------------------------------
      LOGICAL :: Found
      CHARACTER (len=MAX_NAME_LEN) :: DimName         ! Name of the space or time dimension
      CHARACTER (len=MAX_NAME_LEN) :: CoordName       ! Name of the space or time coordinate variable (may be same as dimension)
      CHARACTER (len=MAX_NAME_LEN) :: FileName        ! File name for reading the data (of .nc format)
      CHARACTER (len=MAX_NAME_LEN) :: str
      INTEGER :: i, j, k, dimid, size, VarDim, IndVec(1)
      REAL(KIND=dp) :: FirstTwo(2), LastTwo(2), dx0, dx1, dx2

      ! Opening the NetCDF file
      !------------------------------------------------------------------------------
      FileName = GetString( Params, "Filename", Found )
      IF ( .NOT. Found ) CALL Fatal('GridDataReader','No > Filename < specified')
      NetCDFstatus = NF90_OPEN(FileName,NF90_NOWRITE,FileId)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal( 'GridDataReader', 'NetCDF file could not be opened: '//TRIM(FileName))
      END IF

      UniformCoords = .TRUE.
      x0 = 0.0_dp
      dx = 0.0_dp
      t0 = 0.0_dp
      dt = 0.0_dp
      DimIds = 0
      DimSize = 0
      CoordVarNDims = 0
      TimeSize = 0
      NetDim = 0

      DO i=1,4

        IF( i == 1 ) THEN
           DimName = GetString( Params, "X Dim Name", Found )
           IF ( .NOT. found ) THEN
              DimName = GetString( Params, "X Name", Found )
           END IF
        ELSE IF( i== 2 ) THEN
           DimName = GetString( Params, "Y Dim Name", Found )
           IF ( .NOT. found ) THEN
              DimName = GetString( Params, "Y Name", Found )
           END IF
        ELSE IF( i == 3 ) THEN
           DimName = GetString( Params, "Z Dim Name", Found )
           IF ( .NOT. found ) THEN
              DimName = GetString( Params, "Z Name", Found )
           END IF
        ELSE IF( i == 4 ) THEN
           DimName = GetString( Params, "Time Dim Name", Found )
           IF ( .NOT. found ) THEN
              DimName = GetString( Params, "Time Name", Found )
           END IF
        END IF

        IF(.NOT. Found ) THEN
          IF( i > 2 ) THEN
            CYCLE
          ELSE
            CALL Fatal('GridDataReader',"Unable to find compulsory coordinate name:"//TRIM(DimName))
          END IF
        END IF
        IF( i <= 3 ) NetDim = i

        ! Get the dimension id
        NetCDFstatus = NF90_INQ_DIMID(FileId,DimName,dimid)
        IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal('GridDataReader','Dimension identifier could not be found:'//TRIM(DimName))
        END IF

        IF( i == 1 ) THEN
           CoordName = GetString( Params, "X Var Name", Found )
        ELSE IF( i== 2 ) THEN
           CoordName = GetString( Params, "Y Var Name", Found )
        ELSE IF( i == 3 ) THEN
           CoordName = GetString( Params, "Z Var Name", Found )
        ELSE IF( i == 4 ) THEN
           CoordName = GetString( Params, "Time Var Name", Found )
        END IF

        IF ( .NOT. Found ) THEN
           CoordName = DimName
        END IF

        ! Get the variable id and check whether it is 1D or 2D
        NetCDFstatus = NF90_INQ_VARID(FileId,CoordName,varid)
        IF ( NetCDFstatus /= NF90_NOERR ) THEN
          WRITE(Message,'(A,I0)') 'Variable identifier could not be found: '//TRIM(CoordName)
          CALL Fatal('GridDataReader',Message)
        END IF

        CoordVarIds(i)=varid

        NetCDFstatus = NF90_INQUIRE_VARIABLE(FileId,varid,ndims=vardim)
        SELECT CASE ( VarDim )
        CASE (1)
           WRITE(Message,'(A,I0)') 'Found 1 dimensional coordinate variable > '&
                //TRIM(CoordName)//' < : ',VarDim
           CALL Info('GridDataReader',message,level=7)
        CASE(2)
           WRITE(Message,'(A,I0)') 'Found 2 dimensional coordinate variable (assuming non-uniform) > '&
                //TRIM(CoordName)//' < : ',VarDim
           CALL Info('GridDataReader',message,level=7)
           UniformCoords = .FALSE.
        CASE DEFAULT
           WRITE(Message,'(A,I0)') 'Invalid dimensions for coordinate variable > '&
                //TRIM(CoordName)//' < : ',VarDim
           CALL Fatal('GridDataReader',Message)
        END SELECT
        if (i <= 3) CoordVarNDims(i) = VarDim

        ! Get the size of the coordinate dimension
        NetCDFstatus = NF90_INQUIRE_DIMENSION(FileId,dimid,str,size)
        IF ( NetCDFstatus /= NF90_NOERR ) THEN
          CALL Fatal('GridDataReader','Dimension could not be inquired.')
        END IF
        IF (size <= 1) THEN
           IF (i == 4) THEN
              TimeSize = 1
              dt = 0
           ELSE
              CALL Fatal('GridDataReader','Scalar dimension encountered; No obtainable difference: '//TRIM(str))
           END IF
        END IF

        WRITE(Message,'(A,I0,A,I0)') 'Found dimension > '&
            //TRIM(CoordName)//' < with id ',dimid,' and size ',size
        CALL Info('GridDataReader',Message, Level=6 )

        ! A simple check whether the grid is uniform.  If it is, take the first two
        ! values for each dimension and compute the grid resolution from the data.
        !---------------------------------------------------------------------------

        IF ( (i == 4) .AND. (TimeSize == 1) ) THEN
           ! can't find first 2 if time is scalar, so treat separately
           FirstTwo = 0.0_dp
           IndVec = 1
           NetCDFstatus = NF90_GET_VAR(FileId,varid,FirstTwo(1),IndVec)
           IF ( NetCDFstatus /= NF90_NOERR ) THEN
              CALL Fatal('GridDataReader','NetCDF time values access failed.')
           END IF

        ELSE

           FirstTwo = 0.0_dp
           IndVec = 1
           NetCDFstatus = NF90_GET_VAR(FileId,varid,FirstTwo,IndVec)
           IF ( NetCDFstatus /= NF90_NOERR ) THEN
              CALL Fatal('GridDataReader','NetCDF dimension values access failed.')
           END IF

           IndVec = size - 1
           NetCDFstatus = NF90_GET_VAR(FileId,varid,LastTwo,IndVec)
           IF ( NetCDFstatus /= NF90_NOERR ) THEN
              CALL Fatal('GridDataReader','NetCDF dimension values access failed.')
           END IF

           DimIds(i) = dimid

           dx0 = FirstTwo(2)-FirstTwo(1)
           dx1 = LastTwo(2)-LastTwo(1)
           dx2 = (LastTwo(2)-FirstTwo(1))/(size-1)

           IF( ABS(dx1-dx0) > 1.0d-3 * ABS(dx0) ) THEN
              UniformCoords = .FALSE.
           END IF
           IF( ABS(dx2-dx0) > 1.0d-3 * ABS(dx0) ) THEN
              UniformCoords = .FALSE.
           END IF

           IF (UniformCoords) THEN
              WRITE(Message,'(A,ES12.3)') 'Grid parameter of dimension > '&
                   //TRIM(CoordName)//' < is ',dx0
              CALL Info('GridDataReader',Message, Level=6 )

              WRITE(Message,'(A,2ES12.3,A)') 'Range of dimension > '&
                   //TRIM(CoordName)//' < is [',FirstTwo(1),LastTwo(2),']'
              CALL Info('GridDataReader',Message, Level=6 )
           END IF
        END IF

        IF( i <= 3 ) THEN
           DimSize(i) = size
           x0(i) = FirstTwo(1)
           dx(i) = dx0
        ELSE
           t0 = FirstTwo(1)
           IF (TimeSize /= 1) THEN
              TimeSize = size
              dt = dx0
           END IF
        END IF

      END DO

    END SUBROUTINE NetCDFInit


    !-------------------------------------------------------------------------------
    !> Get data array for given netcdf coordinate variable
    !-------------------------------------------------------------------------------
    SUBROUTINE NetCDFCoordVar1D( i, values )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i ! index of coord var identifier

      REAL (KIND=dp), INTENT(OUT) :: values(:)

      NetCDFstatus = NF90_GET_VAR(FileId,CoordVarIds(i),values)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal('GridDataReader','NetCDF coord variable 1D read failed ')
      END IF

    END SUBROUTINE NetCDFCoordVar1D

    !-------------------------------------------------------------------------------
    SUBROUTINE NetCDFCoordVar2D( i, values )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i ! index of coord var identifier

      REAL (KIND=dp), INTENT(OUT) :: values(:,:)

      NetCDFstatus = NF90_GET_VAR(FileId,CoordVarIds(i),values)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal('GridDataReader','NetCDF coord variable 2D read failed ')
      END IF

    END SUBROUTINE NetCDFCoordVar2D

    !-------------------------------------------------------------------------------
    SUBROUTINE NetCDFCoordVar3D( i, values )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i ! index of coord var identifier

      REAL (KIND=dp), INTENT(OUT) :: values(:,:,:)

      NetCDFstatus = NF90_GET_VAR(FileId,CoordVarIds(i),values)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal('GridDataReader','NetCDF coord variable 3D read failed ')
      END IF

    END SUBROUTINE NetCDFCoordVar3D


    !-------------------------------------------------------------------------------
    !> Set NetCDF index of the variable to be mapped
    !-------------------------------------------------------------------------------
    SUBROUTINE NetCDFVariableInit( VarName )

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: VarName
      !-----------------------------------------------------------------------------

      VarId = 0
      NetCDFstatus = NF90_INQ_VARID(FileId,VarName,VarId)
      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        CALL Fatal('GridDataReader','NetCDF variable name not found: '//TRIM(VarName))
      END IF

      IF(Debug) PRINT *,'NetCDF variable index: ',TRIM(VarName), VarId

    END SUBROUTINE NetCDFVariableInit


    !----------------------------------------------------------------------------------------------
    !> Reads the given variable and returns data on the desired cell on it from NetCDF grid.
    !> There are two versions since ideally the dimensions of the vectors and are different.
    !> This is the 2D version.
    !----------------------------------------------------------------------------------------------
    SUBROUTINE NetCDFDataCell2D( outcome, DimIndex, TimeIndex )
      !------------------------------------------------------
      REAL (KIND=dp), INTENT(OUT) :: outcome(:,:,:)
      INTEGER, INTENT(IN) :: DimIndex(:)
      INTEGER, INTENT(IN) :: TimeIndex
      !------------------------------------------------------
      REAL (KIND=dp) :: stencil3D(2,2,1),stencil2D(2,2)
      INTEGER :: i,j,IndexVector3D(3), CountVector3D(3),&
          IndexVector2D(2),CountVector2D(2)

      ! Access variable and take the values
      !---------------------------------------------------------------------------------------------
      IF( TimeIndex == 0 ) THEN
        CountVector2D = (/ 2, 2 /)
        IndexVector2D = (/ DimIndex(1), DimIndex(2) /)
        NetCDFstatus = NF90_GET_VAR(FileId,VarId,stencil2D,IndexVector2D,CountVector2D)
        outcome(:,:,1) = stencil2D(:,:)
      ELSE
       CountVector3D = (/ 2, 2, 1 /)
        IndexVector3D = (/ DimIndex(1), DimIndex(2), TimeIndex /)
        NetCDFstatus = NF90_GET_VAR(FileId,VarId,stencil3D,IndexVector3D,CountVector3D)
        outcome(:,:,1) = stencil3D(:,:,1)
      END IF

      IF ( NetCDFstatus /= NF90_NOERR ) THEN
        PRINT *,'FileId:',FileId
        PRINT *,'VarId:',VarId
        IF( TimeIndex == 0 ) THEN
          PRINT *,'IndexVector:',IndexVector2D
          PRINT *,'CountVector:',CountVector2D
        ELSE
          PRINT *,'IndexVector:',IndexVector3D
          PRINT *,'CountVector:',CountVector3D
        END IF
        CALL Fatal('GridDataReader','NetCDF variable access failed in 2D.')
      END IF

    END SUBROUTINE NetCDFDataCell2D


    !----------------------------------------------------------------------------------------------
    !> The 3D version of the previous routine.
    !----------------------------------------------------------------------------------------------
   SUBROUTINE NetCDFDataCell3D( outcome, DimIndex, TimeIndex )
      !------------------------------------------------------
      REAL (KIND=dp), INTENT(OUT) :: outcome(:,:,:)
      INTEGER, INTENT(IN) :: DimIndex(:)
      INTEGER, INTENT(IN) :: TimeIndex
      !------------------------------------------------------
      REAL (KIND=dp) :: stencil4D(2,2,2,1),stencil3d(2,2,2)
      INTEGER :: i,j
      INTEGER :: IndexVector3D(3), CountVector3D(3), &
          IndexVector4D(4), CountVector4D(4)

      ! Access variable and take the values
      !---------------------------------------------------------------------------------------------
      IF( TimeIndex == 0 ) THEN
        CountVector3D = (/ 2, 2, 2 /)
        IndexVector3D = (/ DimIndex(1), DimIndex(2), DimIndex(3) /)
        NetCDFstatus = NF90_GET_VAR(FileId,VarId,stencil3D,IndexVector3D,CountVector3D)
        outcome(:,:,:) = stencil3D(:,:,:)
      ELSE
        CountVector4D = (/ 2, 2, 2, 1 /)
        IndexVector4D = (/ DimIndex(1), DimIndex(2), DimIndex(3), TimeIndex /)
        NetCDFstatus = NF90_GET_VAR(FileId,VarId,stencil4D,IndexVector4D,CountVector4D)
        outcome(:,:,:) = stencil4D(:,:,:,1)
      END IF

       IF ( NetCDFstatus /= NF90_NOERR ) THEN
        PRINT *,'FileId:',FileId
        PRINT *,'VarId:',VarId
        IF( TimeIndex == 0 ) THEN
          PRINT *,'IndexVector:',IndexVector3D
          PRINT *,'CountVector:',CountVector3D
        ELSE
          PRINT *,'IndexVector:',IndexVector4D
          PRINT *,'CountVector:',CountVector4D
        END IF
        CALL Fatal('GridDataReader','NetCDF variable access failed in 3D.')
      END IF

    END SUBROUTINE NetCDFDataCell3D


    !----------------------------------------------------------------------------------
    ! Closes the active NetCDF file.
    !----------------------------------------------------------------------------------
    SUBROUTINE NetCDFClose( )
      USE NetCDF
      IMPLICIT NONE
      INTEGER :: status

      status = NF90_CLOSE(FileId)
      IF ( status /= NF90_NOERR ) THEN ! Error encountered
        CALL Fatal( 'GridDataReader', 'Failed to close NetCDF file' )
      END IF

    END SUBROUTINE NetCDFClose

  END MODULE NetCDFInterface
  !----------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Solver for mapping data from uniform grids into Elmer mesh.
!> Currently netcdf interface has been implemented.
!------------------------------------------------------------------------------
SUBROUTINE GridDataReader( Model,Solver,dtime,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils
  USE MeshUtils
  USE ElementUtils
  USE NetCDFInterface

  IMPLICIT NONE

  TYPE array3D
     REAL(KIND=dp), POINTER :: values(:,:,:)
  END TYPE array3D

  !------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dtime
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  LOGICAL :: Debug = .FALSE.
  TYPE(Variable_t), POINTER :: FieldVar, PrevFieldVar
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(ValueList_t), POINTER :: Params
  TYPE(array3D) :: coordVar(3)
  INTEGER :: n,k,node,MeshDim, NetDim,iTime,nTime
  INTEGER, POINTER :: FieldPerm(:)
  REAL(KIND=dp), POINTER :: Field(:), FieldOldValues(:)=>NULL()
  REAL(KIND=dp) :: x(3),dx(3),x0(3),x1(3),u1(3),u2(3),dt,t0,val,&
                   Eps(3),Time,x0e(3),x1e(3),pTime,EpsTime,q,r
  INTEGER :: DimSize(3), CoordVarNDims(3), i
  INTEGER :: TimeSize, IntTimeIndex,tnmax, NoVar, InterpStatus
  INTEGER :: status, time_begin,time_end,MaskNodes
  INTEGER :: StatusCount(6)
  CHARACTER (len = MAX_NAME_LEN) :: str, VarName, TargetName, MaskName, &
      CoordSystem, TimeInterpolationMethod
  REAL(KIND=dp) :: Coeff, InterpMultiplier, InterpBias, TimeIndex, acc, MinMaxVals(2), &
      MissingVal
  LOGICAL :: Found, IsTime, DoCoordinateTransformation, DoCoordMapping, &
      DoScaling, DoBoundingBox, DoPeriodic, UniformCoords, HaveMinMax, DoNuInterpolation,&
      KeepOld
  INTEGER, POINTER :: CoordMapping(:), PeriodicDir(:)

  ! General initializations
  !------------------------------------------------------------------------------

  CALL Info('GridDataReader','-----------------------------------------', Level=5 )
  CALL Info('GridDataReader','Obtaining field(s) from grid NetCDF format',Level=5 )
  CALL ResetTimer('GridDataReader')

  !-- Pointer declarations
  PSolver => Solver
  Params => Solver % Values
  Mesh => Solver % Mesh
  MeshDim = Mesh % MeshDim

  ! Elmer resolution and coordinate system
  !------------------------------------------------------------------------------
  CALL InitEpsilon( Params, Eps, EpsTime )

  CoordSystem = GetString( Params,'Coordinate Transformation',&
      DoCoordinateTransformation)

  CoordMapping => ListGetIntegerArray( Params,'Coordinate Mapping',&
      DoCoordMapping )
  IF( DoCoordMapping ) THEN
    IF ( SIZE(CoordMapping) /= 3 ) THEN
      WRITE( Message, * ) 'Invalid size of coordinate mapping: ', SIZE(CoordMapping)
      CALL Fatal( 'GridDataReader', Message )
    END IF
    DO i=1,3
      IF( .NOT. ANY( CoordMapping == i ) ) THEN
        WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
        CALL Fatal( 'GridDataReader', Message )
      END IF
    END DO
  END IF

  PeriodicDir => ListGetIntegerArray( Params,'Periodic Directions',&
      DoPeriodic )

  KeepOld = ListGetLogical( Params, 'Out Of Bounds Retain Previous Value', Found)
  IF(.NOT. Found) KeepOld = .FALSE.

  !------------------------------------------------------------------------------
  ! Initialize NetCDF data locations, sizes and resolution
  !------------------------------------------------------------------------------
  CALL NetCDFInit(Params, NetDim, DimSize, CoordVarNDims, TimeSize, x0, dx, t0, dt, UniformCoords )

  IF( NetDim < 1 .OR. NetDim > 3 ) THEN
    CALL Fatal('GridDataReader','NetCDF dimensions must be either 2 or 3!')
  END IF
  IF (MeshDim > NetDim ) THEN
    CALL Info('GridDataReader','Omitting the extra dimensions of Elmer mesh',Level=6)
  END IF
  x1 = x0 + (DimSize-1) * dx

  !------------------------------------------------------------------------------
  ! If any coordinate variables have more than one dimension and or are non uniform
  ! then scrap assumptions about coordinates and read in the full coordinate
  ! variables for use in interpolation.
  !------------------------------------------------------------------------------
  DoNuInterpolation = (.NOT. UniformCoords) .OR. (MAXVAL(CoordVarNDims) > 1 )

  IF( DoNuInterpolation ) THEN
    WRITE(Message, '(A)') 'Coordinate variables are non-uniform and or more than &
        &one dimensional.  CELLS ASSUMED RECTANGULAR IN THIS CODE.'
    CALL Warn('GridDataReader',Message)

     DO i=1,3 ! spatial dimensions, time is not considered in this loop.
        SELECT CASE (CoordVarNDims(i))
        CASE (0)
           ALLOCATE(coordVar(i)%values(1,1,1))
        CASE (1)
           ALLOCATE(coordVar(i)%values(DimSize(i),1,1))
           CALL NetCDFCoordVar( i, coordVar(i)%values(:,1,1) )
        CASE (2)
           IF (i /= 3) THEN
              ! if first two coord variables are 2D assume their 2 dimensions are
              ! the first 2 dims, i.e. assume they both refer to the 2 horizontal
              ! dims.
              ALLOCATE(coordVar(i)%values(DimSize(1),DimSize(2),1))
              CALL NetCDFCoordVar( i, coordVar(i)%values(:,:,1) )
           ELSE
              CALL Fatal('GridDataReader',"Third coordinate variable has 2 dimensions, not expected, code some more...")
           END IF
        CASE(3)
              ALLOCATE(coordVar(i)%values(DimSize(1),DimSize(2),DimSize(3)))
              CALL NetCDFCoordVar( i, coordVar(i)%values(:,:,:) )
              CALL Fatal('GridDataReader',"Coordinate variable has 3 dimensions, not expected, code some more...")
        CASE DEFAULT
              CALL Fatal('GridDataReader',"Coordinate variable has more than 3 dimensions, not expected, code some more...")
        END SELECT
     END DO
  END IF


  !------------------------------------------------------------------------------
  ! Optionally map the Elmer mesh so that it coincides with the NetCDF mesh
  ! This is intended mainly for testing purposes etc.
  !
  !------------------------------------------------------------------------------
  DoScaling = ListGetLogical( Params,'Enable Scaling',Found )

  DoBoundingBox = ListGetLogical( Params,'Check Bounding Box',Found )
  IF( DoScaling ) THEN
    IF( DoCoordinateTransformation .OR. DoCoordMapping ) THEN
      CALL Fatal('GridDataReader','Cannot do scaling and mapping together!')
    ELSE
      DoBoundingBox = .TRUE.
    END IF
  END IF

  IF( DoBoundingBox ) THEN
    x = 0.0_dp
    x0e = HUGE( x0e )
    x1e = -HUGE( x1e )

    DO node=1, Mesh % NumberOfNodes
      IF( ASSOCIATED( FieldPerm ) ) THEN
        k = FieldPerm(node)
        IF( k == 0 ) CYCLE
      ELSE
        k = node
      END IF

      x(1) = Mesh % Nodes % x(node)
      x(2) = Mesh % Nodes % y(node)
      IF( NetDim == 3 ) x(3) = Mesh % Nodes % z(node)

      IF( DoCoordinateTransformation ) THEN
        x = LocalCoordinateTransformation( x, CoordSystem )
      END IF

      IF( DoCoordMapping ) THEN
        x = x( CoordMapping(1:NetDim) )
      END IF

      x0e = MIN( x0e, x )
      x1e = MAX( x1e, x )
    END DO

    DO i=1,NetDim
      WRITE(Message,'(A,I0,A,2ES12.3,A)') 'Initial Elmer coordinate ',i,' range is [',x0e(i),x1e(i),']'
      CALL Info('GridDataReader',Message, Level=6 )
    END DO

    IF( DoScaling ) THEN
      Found = .FALSE.
      DO i=1,NetDim
        q = (x1(i) - x0(i) ) / ( x1e(i) - x0e(i) )
        r = x0(i) - q * x0e(i)

        acc = MAX( ABS(r), ABS( q- 1.0) )

        ! If the bounding box of the two meshes is the same do nothing
        IF( acc < 1.0d-8 ) CYCLE

        Found = .TRUE.
        IF( i == 1 ) THEN
          Mesh % Nodes % x = r + q * Mesh % Nodes % x
        ELSE IF( i == 2 ) THEN
          Mesh % Nodes % y = r + q * Mesh % Nodes % y
        ELSE
          Mesh % Nodes % z = r + q * Mesh % Nodes % z
        END IF
      END DO

      IF( Found ) THEN
        x0e(1) = MINVAL(Mesh % Nodes % x)
        x0e(2) = MINVAL(Mesh % Nodes % y)
        x0e(3) = MINVAL(Mesh % Nodes % z)
        x1e(1) = MAXVAL(Mesh % Nodes % x)
        x1e(2) = MAXVAL(Mesh % Nodes % y)
        x1e(3) = MAXVAL(Mesh % Nodes % z)

        DO i=1,NetDim
          WRITE(Message,'(A,I0,A,2ES12.3,A)') 'Modified Elmer coordinate ',i,' range is [',x0e(i),x1e(i),']'
          CALL Info('GridDataReader',Message, Level=6 )
        END DO
      END IF
    ELSE
      IF( ALL( x0e >= x0 ) .AND. ALL( x1e <= x1 ) ) THEN
        CALL Info('GridDataReader','Elmer bounding box is within the NetCDF one!',Level=6)
      END IF

    END IF
  END IF



  !--------------------------------------------------------------------------------------
  ! Get the timestep at which interpolation is desired
  ! If the time does not coincide with a timestep in the file, two timesteps are needed.
  !--------------------------------------------------------------------------------------
  IF( TimeSize == 0 ) THEN
    CALL Info('GridDataReader','No time given, using 1st step',Level=6)
    IntTimeIndex = 0
    nTime = 1
    pTime = 1.0_dp
 ELSE
    CALL GetTimePoint(Params, t0, dt, TimeIndex )

    IntTimeIndex = NINT( TimeIndex )

    IF( ABS( TimeIndex - IntTimeIndex ) > EpsTime ) THEN
      IntTimeIndex = FLOOR( TimeIndex )
      nTime = 2
      pTime = 1 + IntTimeIndex - TimeIndex
      WRITE (Message,'(A,ES10.3,A)') 'Given time value ', TimeIndex , ' using time interpolation.'
    ELSE
      nTime = 1
      pTime = 1.0_dp
      WRITE (Message,'(A,ES10.3,A)') 'Given time value ', TimeIndex, ', no time interpolation used.'
    END IF
    CALL Info('GridDataReader',Message,Level=6)
  END IF

  IF( Debug ) THEN
    PRINT *,'nTime B',nTime,TimeIndex,pTime
  END IF

  !-------------------------------------------------------------------------------
  ! Loop over variables to be mapped.
  ! If no target variable is given it will be the same name as the primary variable.
  !-------------------------------------------------------------------------------
  NULLIFY( PrevFieldVar )

  NoVar = 0
  DO WHILE (.TRUE.)
    ! Get NetCDF variable
    !-----------------------------------------------------------------------------
    NoVar = NoVar + 1
    WRITE( str,'(A,I0)') 'Variable ',NoVar
    VarName = GetString( Params,str, Found )
    IF(.NOT. Found ) THEN
      IF( NoVar == 1 ) THEN
        CALL Fatal('GridDataReader','Calling subroutine without > Variable 1 < defined!')
      END IF
      EXIT
    END IF

    CALL Info('GridDataReader','Performing interpolation for variable: '//TRIM(VarName) )
    CALL NetCDFVariableInit( VarName )

    ! Get Elmer variable, if not present create it.
    !-------------------------------------------------------------------------------
    WRITE( str,'(A,I0)') 'Target Variable ',NoVar
    TargetName = GetString( Params,str,Found )
    IF( .NOT. Found ) TargetName = VarName
    FieldVar => VariableGet( Mesh % Variables,TargetName )
    IF( .NOT. ASSOCIATED( FieldVar ) ) THEN
      WRITE( str,'(A,I0)') 'Mask Name ',NoVar
      MaskName = GetString( Params,str, Found )

      NULLIFY(FieldPerm)
      ALLOCATE( FieldPerm( Mesh % NumberOfNodes ) )

      IF( Found ) THEN
        CALL MakePermUsingMask( Model, Solver, Mesh, MaskName,.FALSE.,FieldPerm,&
            MaskNodes,RequireLogical=.TRUE.)
        IF( MaskNodes == 0 ) THEN
          DEALLOCATE( FieldPerm )
          CALL Fatal('GridDataReader','No active nodes for mask: '//TRIM(MaskName))
        END IF

        WRITE (Message, '(A,I0,A,I0)') 'The mask > '//TRIM(MaskName)//' < resulted to '&
            ,MaskNodes,' nodes out of ',Mesh % NumberOfNodes
        CALL Info('GridDataReader',Message,Level=6)

        CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,TargetName,1,Perm=FieldPerm)
        FieldVar => VariableGet( Mesh % Variables,TargetName )
        NULLIFY(FieldPerm)
      ELSE
         FieldPerm = [(i,i=1,Mesh % NumberOfNodes)]
         CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,TargetName,1, Perm=FieldPerm)
        FieldVar => VariableGet( Mesh % Variables,TargetName )
      END IF
    END IF
    Field => FieldVar % Values
    FieldPerm => FieldVar % Perm

    !In case the user wants unfound values to retain their previous values
    !(i.e. advection out of domain)
    IF(KeepOld) THEN
       IF(ASSOCIATED(FieldOldValues)) DEALLOCATE(FieldOldValues)
       ALLOCATE(FieldOldValues(SIZE(Field)))
       FieldOldValues = Field
    END IF

    ! Set a constant background to the field. This can be done a priori
    ! since it does not interfere with the interpolation.
    !----------------------------------------------------------------------
    str = 'Interpolation Offset'
    InterpBias = GetCReal( Params,str,Found )
    IF( .NOT. Found ) THEN
      WRITE( str,'(A,I0)') TRIM(str)//' ',NoVar
      InterpBias = GetCReal( Params,str,Found )
      IF( .NOT. Found ) InterpBias = 0.0_dp
    END IF

    ! If the target variable is the same, then obviously it is so intentionally
    ! to combine values from two different data sets.
    !--------------------------------------------------------------------------
    IF( .NOT. ASSOCIATED( PrevFieldVar, FieldVar ) ) THEN
      Field = InterpBias
    END IF
    PrevFieldVar => FieldVar


    ! Multiply the field with a constant
    ! Here just the constant is obtained, multiplication is done later.
    !----------------------------------------------------------------------
    str = 'Interpolation Multiplier'
    InterpMultiplier = GetCReal( Params,str,Found )
    IF( .NOT. Found ) THEN
      WRITE( str,'(A,I0)') TRIM(str)//' ',NoVar
      InterpMultiplier = GetCReal( Params,str,Found )
      IF(.NOT. Found ) InterpMultiplier = 1.0_dp
    END IF

    ! Set the acceptable range of values that can be used in the interpolation
    !-------------------------------------------------------------------------
    str = 'Valid Min Value'
    MinMaxVals(1) = GetCReal( Params,str,Found )
    IF( .NOT. Found ) THEN
      WRITE( str,'(A,I0)') TRIM(str)//' ',NoVar
      MinMaxVals(1) = GetCReal( Params,str,Found )
      IF(.NOT. Found ) MinMaxVals(1) = -HUGE(MinMaxVals(1))
    END IF
    HaveMinMax = Found

    str = 'Valid Max Value'
    MinMaxVals(2) = GetCReal( Params,str,Found )
    IF( .NOT. Found ) THEN
      WRITE( str,'(A,I0)') TRIM(str)//' ',NoVar
      MinMaxVals(2) = GetCReal( Params,str,Found )
      IF(.NOT. Found ) MinMaxVals(2) = HUGE(MinMaxVals(2))
    END IF
    HaveMinMax = HaveMinMax .OR. Found

    str = 'Default Value'
    MissingVal = GetCReal( Params,str,Found )
    IF( .NOT. Found ) THEN
      WRITE( str,'(A,I0)') TRIM(str)//' ',NoVar
      MissingVal = GetCReal( Params,str,Found )
    END IF


    ! Loop over 1 or 2 timesteps
    !----------------------------------------------------------------------
    StatusCount = 0
    DO iTime = 1, nTime

      IF( iTime == 2 ) THEN
        IntTimeIndex = IntTimeIndex + 1
        pTime = 1.0_dp - pTime
      END IF
      Coeff = InterpMultiplier * pTime

      IF ( IntTimeIndex == 0 ) THEN
        CONTINUE
      ELSE IF ( IntTimeIndex < 1 .OR. IntTimeIndex > TimeSize ) THEN
        WRITE (Message, '(A,I0,A,I0,A)') 'Time value ', IntTimeIndex, ' is out of range (1,',TimeSize, ')'
        CALL Warn('GridDataReader',Message)
      END IF

      ! Go through the active nodes and perform interpolation
      !---------------------------------------------------------------------------
      x = 0.0_dp
      DO node=1, Mesh % NumberOfNodes
        IF( ASSOCIATED( FieldPerm ) ) THEN
           k = FieldPerm(node)
           IF( k == 0 ) CYCLE
        ELSE
           k = node
        END IF

        ! The Elmer point of interest
        ! Use the leading dimension of NetCDF data - not of Elmer.
        !-------------------------------------------------------------------------
        x(1) = Mesh % Nodes % x(node)
        x(2) = Mesh % Nodes % y(node)
        IF( NetDim == 3 ) x(3) = Mesh % Nodes % z(node)

        ! Coordinate mapping from Elmer (x,y) to the one used by NetCDF.
        !-------------------------------------------------------------------------
        IF( DoCoordinateTransformation ) THEN
          x = LocalCoordinateTransformation( x, CoordSystem )
        END IF

        IF( DoCoordMapping ) THEN
          x = x( CoordMapping(1:NetDim) )
        END IF

        IF( DoPeriodic ) THEN
          DO i=1,NetDim
            IF( PeriodicDir(i) > 0 ) THEN
              IF( x(i) > x1(i) ) x(i) = x(i) - (x1(i) - x0(i) )
              IF( x(i) < x0(i) ) x(i) = x(i) + (x1(i) - x0(i) )
            END IF
          END DO

        END IF

        IF ( DoNuInterpolation ) THEN
           InterpStatus = NUInterpolation(NetDim,x,DimSize,Eps,IntTimeIndex,val,coordVar,coordVarNDims)
        ELSE
          IF( HaveMinMax ) THEN
            InterpStatus = FDInterpolation(NetDim,x,DimSize,x0,dx,x1,Eps,IntTimeIndex,val,MinMaxVals)
          ELSE
            InterpStatus = FDInterpolation(NetDim,x,DimSize,x0,dx,x1,Eps,IntTimeIndex,val)
          END IF
        END IF

        IF( InterpStatus == 1 ) THEN
          CONTINUE
        ELSE IF( InterpStatus == 2 ) THEN
          val = MissingVal
        ELSE IF( InterpStatus == 3 .OR. InterpStatus == 4) THEN
          val = MinMaxVals(1)
        ELSE IF( InterpStatus == 5 .OR. InterpStatus == 6) THEN
          val = MinMaxVals(2)
        ELSE
          CALL Fatal('GridDataReader','Unknown InterpStatus!')
        END IF

        Field(k) = Field(k) + Coeff * val

        IF(KeepOld .AND. (InterpStatus == 2)) Field(k) = FieldOldValues(k)

        StatusCount(InterpStatus) = StatusCount(InterpStatus) + 1

      END DO
    END DO

    IF(KeepOld) DEALLOCATE(FieldOldValues)

    IF( StatusCount(1) > 0) THEN
      WRITE( Message,'(A,I0)')'Number of proper mappings  : ',StatusCount(1)
      CALL Info('GridDataReader',Message)
    END IF
    IF( StatusCount(2) > 0) THEN
      WRITE( Message,'(A,I0)')'Number of missing mappings : ',StatusCount(2)
      CALL Warn('GridDataReader',Message)
    END IF
    IF( StatusCount(3) > 0 ) THEN
      WRITE( Message,'(A,I0)')'Number of some too small values : ',StatusCount(3)
      CALL Warn('GridDataReader',Message)
    END IF
    IF( StatusCount(4) > 0 ) THEN
      WRITE( Message,'(A,I0)')'Number of all too small values : ',StatusCount(4)
      CALL Warn('GridDataReader',Message)
    END IF
    IF( StatusCount(5) > 0) THEN
      WRITE( Message,'(A,I0)')'Number of some too large values : ',StatusCount(5)
      CALL Warn('GridDataReader',Message)
    END IF
    IF( StatusCount(6) > 0) THEN
      WRITE( Message,'(A,I0)')'Number of all too large values : ',StatusCount(6)
      CALL Warn('GridDataReader',Message)
    END IF
  END DO

  CALL NetCDFClose()

  DO i = 1,3
     IF (ASSOCIATED(coordVar(i)%values)) NULLIFY(coordVar(i)%values)
  END DO

  CALL CheckTimer('GridDataReader',Delete=.TRUE.)
  CALL Info('GridDataReader','All done',Level=5)
  CALL Info('GridDataReader', '-----------------------------------------', Level=5 )

!  WRITE(Message,'(A,ES20.10)') 'NCTESTSTAT: ',SUM(field)
!  CALL Info('GridDataReader',Message,Level=2)

CONTAINS

  !------------------------------------------------------------------------------
  ! Initializes the resolution used in interpolation.
  !------------------------------------------------------------------------------
  SUBROUTINE InitEpsilon( Params, Eps, EpsTime )
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Params
    REAL(KIND=dp), INTENT(INOUT) :: Eps(:), EpsTime
    !-----------------------------------------------------------------
    LOGICAL :: Found
    REAL(KIND=dp) :: eX, eY

    ! Epsilons are the relative tolerances for the amount
    ! the Elmer grid point misses the bounds of the NetCDF bounding box
    !-------------------------------------------------------------------
    Eps = 0.0_dp
    Eps(1) = GetConstReal(Params, "X Epsilon", Found )
    IF ( .NOT. Found ) THEN
      CALL Warn('GridDataReader', 'Keyword > X Epsilon < not given, setting to default eps')
      Eps(1) = EPSILON( Eps(1) )
    END IF
    Eps(2) = GetConstReal(Params, "Y Epsilon", Found )
    IF ( .NOT. Found ) THEN
      CALL Info('GridDataReader', 'Keyword > Y Epsilon < not given, setting equal to > X Epsilon <',Level=6)
      Eps(2) = Eps(1)
    END IF
    IF( NetDim == 3 ) THEN
      Eps(3) = GetConstReal(Params, "Z Epsilon", Found )
      IF ( .NOT. Found ) THEN
        CALL Info('GridDataReader', 'Keyword > Z Epsilon < not given, setting equal to > X Epsilon <',Level=6)
        Eps(3) = Eps(1)
      END IF
    END IF

    EpsTime = GetConstReal(Params, "Time Epsilon", Found )
    IF(.NOT. Found) EpsTime = EPSILON( EpsTime )

    IF( Debug ) THEN
      PRINT *,'Eps',Eps,'EpsTime',EpsTime
    END IF

  END SUBROUTINE InitEpsilon


  !--------------------------------------------------------------------
  ! Initializes the time values
  !--------------------------------------------------------------------
  SUBROUTINE GetTimePoint( Params, t0, dt, TimeIndex )
  !----------------------------------------------------
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Params
    REAL(KIND=dp), INTENT(IN) :: t0, dt
    REAL(KIND=dp), INTENT(OUT) :: TimeIndex
    !-----------------------------------------------------------------------
    LOGICAL :: IsTimeIndex, Found
    REAL(KIND=dp) :: Time, Coeff
    INTEGER :: VisitedTimes =  0

    SAVE VisitedTimes

    VisitedTimes = VisitedTimes + 1
    IF( GetLogical( Params, "Is Time Counter", Found ) ) THEN
      TimeIndex = VisitedTimes
      RETURN
    END IF

    ! Get user-specified time or true physical time
    Time = GetCReal( Params, "Time Point", Found )
    IF( .NOT. Found ) THEN
      Time = GetTime()

      ! Add possible offset in time
      Coeff = GetCReal( Params, "Time Offset", Found )
      IF( Found ) Time = Time + Coeff

      ! Add possible multiplicator in time
      Coeff = GetCReal( Params, "Time Multiplier", Found )
      IF( Found ) Time = Coeff * Time
    END IF

    ! Check if time is assumed to be index variable, or true time
    IsTimeIndex = GetLogical( Params, "Is Time Index", Found )
    IF ( IsTimeIndex) THEN
      TimeIndex = Time
    ELSE
      TimeIndex =  1.0 + ( Time - t0 ) / dt
    END IF

  END SUBROUTINE GetTimePoint


  !--------------------------------------------------------------
  !> Performs linear interpolation
  !--------------------------------------------------------------
  FUNCTION LinearInterpolation(x,u1,u2) RESULT(y)
    REAL(KIND=dp), INTENT(IN) :: u1(2),u2(2),x
    REAL(KIND=dp) :: y

    y = (((u2(2) - u1(2))/(u2(1) - u1(1)))*(x-u1(1)))+u1(2)
  END FUNCTION LinearInterpolation

  !----------------------------------------------------------------------------
  !> Performs bilinear interpolation on a stencil (2x2 matrix of corner values)
  !----------------------------------------------------------------------------
  FUNCTION BiLinearInterpolation(stencil,weights) RESULT(val)
    REAL(KIND=dp), INTENT(IN) :: stencil(:,:), weights(:)
    REAL(KIND=dp) :: val,w(4)

    val = stencil(1,1)*(1-weights(1))*(1-weights(2)) + &
        stencil(2,1)*weights(1)*(1-weights(2)) + &
        stencil(1,2)*(1-weights(1))*weights(2) + &
        stencil(2,2)*weights(1)*weights(2)

  END FUNCTION BiLinearInterpolation

  !------------------------------------------------------------------------------
  !> Performs trilinear interpolation on a stencil (2x2x2 matrix of corner values)
  !------------------------------------------------------------------------------
  FUNCTION TriLinearInterpolation(stencil,weights) RESULT(val)
    REAL(KIND=dp), INTENT(IN) :: stencil(:,:,:), weights(:)
    REAL(KIND=dp) :: val, val1, val2

    val1 = BiLinearInterpolation(stencil(:,:,1),weights(1:2))
    val2 = BiLinearInterpolation(stencil(:,:,2),weights(1:2))
    val = val1*(1-weights(3)) + val2*weights(3)
  END FUNCTION TriLinearInterpolation

  !-------------------------------------------------------------------------------
  !> Non-Uniform interpolation.  Interpolation where the full coordinate variables
  !  are used because they are non-uniform or have more than one dimension.  Uses
  !  (bi)linear interpolation.
  !-------------------------------------------------------------------------------
  FUNCTION NUInterpolation(NetDim,x,DimSize,Eps,TimeIndex,val,coordVar,coordVarNDims) &
       RESULT( InterpStatus )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NetDim, DimSize(:), CoordVarNDims(:)
    INTEGER, INTENT(IN) :: TimeIndex
    REAL(KIND=dp), INTENT(IN) :: x(:),Eps(:)
    TYPE(array3D), INTENT(IN) :: coordVar(:)
    REAL(KIND=dp), INTENT(OUT) :: val ! Final Elmer point and interpolated value
    INTEGER :: InterpStatus

    INTEGER       :: Ind(NetDim) ! lower index of netcdf cell containing Elmer node
    REAL(KIND=dp) :: stencil(2,2,2)
    REAL(KIND=dp) :: Weights(netDim)
    LOGICAL :: success

    ! check coord variable dimensions
    IF (MAXVAL(CoordVarNDims) == 1) THEN
       ! all coord variable dimensions are 1D so search each separately
       DO i = 1,NetDim
          success = findCell1D(coordVar(i),x(i),Ind(i),weights(i))
          WRITE(Message, '(A)') 'Not yet tested this combination of netcdf &
               &coordinate variable dimensions &
               &(non-uniform, all single), pls remove this comment if it works...'
          CALL Warn('GridDataReader',Message)
       END DO

    ELSEIF ((CoordVarNDims(1) == 2) .AND. (CoordVarNDims(2) == 2)) THEN
       ! First 2 coordinate variables are 2D, handle them together
       success = findCell2D(coordVar(1:2),x(1:2),Ind(1:2),weights(1:2))

       IF (CoordVarNDims(3) == 1) THEN
          success = findCell1D(coordVar(3),x(3),Ind(3),weights(3))
          WRITE(Message, '(A)') 'Not yet tested this combination of netcdf &
               &coordinate variable dimensions &
               &(2,2,1), pls remove this comment if it works...'
          CALL Warn('GridDataReader',Message)
       END IF

    ELSE
       WRITE(Message, '(A)') 'Cannot handle this combination of netcdf &
            &coordinate variable dimensions, code some more...'
       CALL Fatal('GridDataReader',message)
    END IF

    SELECT CASE(netDim)
    CASE(2)
       CALL NetCDFDataCell2D( stencil, Ind, TimeIndex  )
       val = BiLinearInterpolation(stencil(:,:,1),weights)
    CASE(3)
       CALL NetCDFDataCell3D( stencil, Ind, TimeIndex  )
       val = TriLinearInterpolation(stencil,weights)
    CASE DEFAULT
       CALL FATAL('GridDataReader','expected 2 or 3 dimensional netcdf data')
    END SELECT

    InterpStatus = 1

  END FUNCTION NUInterpolation

  !-------------------------------------------------------------------------------
  FUNCTION findCell2D(coordVar,xe,ind,weights) &
       RESULT( success )

    IMPLICIT NONE

    TYPE(array3D), INTENT(IN)  :: coordVar(:)
    REAL(KIND=dp), INTENT(IN)  :: xe(:)

    INTEGER, INTENT(OUT)       :: ind(:)
    REAL(KIND=dp), INTENT(OUT) :: Weights(:)

    LOGICAL :: success

    REAL(KIND=dp),POINTER :: dist(:,:),distCum(:,:)
    REAL(KIND=dp)         :: xu(2),yu(2)   ! unit vectors in x and y coord directions
    REAL(KIND=dp)         :: xer(2)        ! Elmer node relative to netcdf ll cell corner
    INTEGER               :: nx,ny

    nx = SIZE(coordVar(1)%values(:,1,1))
    ny = SIZE(coordVar(1)%values(1,:,1))

    ALLOCATE( dist    ( nx,  ny   ) )
    ALLOCATE( distCum ( nx-1,ny-1 ) )

    ! distance of each netcdf point from the current Elmer node
    dist = SQRT( (xe(1)-coordVar(1)%values(1:nx,1:ny,1))**2 + &
         (xe(2)-coordVar(2)%values(1:nx,1:ny,1) )**2 )

    ! cumlative distance of four netcdf points (comprising a cell) from current ELmer node
    distCum = dist(1:nx-1,1:ny-1) + dist(2:nx,1:ny-1)+ &
         dist(2:nx,2:ny) + dist(1:nx-1,2:ny)

    ind = MINLOC(distCum)

    WRITE(Message,'(A,ES12.3)') 'Found cell with min dist ',MINVAL(dist)
    CALL Info('GridDataReader',message,level=8)

    IF(( xe(2) < MINVAL(coordVar(2)%values(ind(1):ind(1)+1,ind(2):ind(2)+1,1)) ) .OR. &
       ( xe(2) > MAXVAL(coordVar(2)%values(ind(1):ind(1)+1,ind(2):ind(2)+1,1)) ) .OR. &
       ( xe(1) < MINVAL(coordVar(1)%values(ind(1):ind(1)+1,ind(2):ind(2)+1,1)) ) .OR. &
       ( xe(1) > MAXVAL(coordVar(1)%values(ind(1):ind(1)+1,ind(2):ind(2)+1,1)) ) ) THEN
       WRITE(Message, '(A)') 'Elmer node not contained in netcdf cell, need to implement &
            &proper tolerance checks (Epsilon not currently used in this case)'
       CALL Warn('GridDataReader',Message)
    END IF

!    print*,shape(xer),shape(xe)
!    print*,11shape(coordVar(1)%values(ind(1),ind(2)))

    xer(1) = xe(1) - coordVar(1)%values(ind(1),ind(2),1)
    xer(2) = xe(2) - coordVar(2)%values(ind(1),ind(2),1)

    xu(1)  = coordVar(1)%values(ind(1)+1,ind(2),1) - coordVar(1)%values(ind(1),ind(2),1)
    xu(2)  = coordVar(2)%values(ind(1)+1,ind(2),1) - coordVar(2)%values(ind(1),ind(2),1)

    yu(1)  = coordVar(1)%values(ind(1),ind(2)+1,1) - coordVar(1)%values(ind(1),ind(2),1)
    yu(2)  = coordVar(2)%values(ind(1),ind(2)+1,1) - coordVar(2)%values(ind(1),ind(2),1)

    ! normalise xer and unit vectors
    xer(1) = xer(1) / SQRT(xu(1)**2+xu(2)**2)
    xer(2) = xer(2) / SQRT(yu(1)**2+yu(2)**2)
    xu     = xu / SQRT(xu(1)**2+xu(2)**2)
    yu     = yu / SQRT(yu(1)**2+yu(2)**2)

    weights(1) = SUM(xer * xu)
    weights(2) = SUM(xer * yu)

    DEALLOCATE(dist)
    DEALLOCATE(distCum)

    success = .TRUE.

  END FUNCTION findCell2D


  FUNCTION findCell1D(coordVar,xe,ind,weights) &
       RESULT( success )

    IMPLICIT NONE

    TYPE(array3D), INTENT(IN)  :: coordVar
    REAL(KIND=dp), INTENT(IN)  :: xe

    INTEGER, INTENT(OUT)       :: ind
    REAL(KIND=dp), INTENT(OUT) :: Weights

    LOGICAL :: success

    ind     = 0
    Weights = 0.0_dp

    success = .FALSE.

  END FUNCTION findCell1D

  !-------------------------------------------------------------------------------
  !> Interpolates one grid point given data on a finite difference stencil.
  !-------------------------------------------------------------------------------
  FUNCTION FDInterpolation(NetDim,x,DimSize,x0,dx,x1,Eps,TimeIndex,val,MinMax) &
       RESULT( InterpStatus )

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NetDim,DimSize(:)
    INTEGER, INTENT(IN) :: TimeIndex
    REAL(KIND=dp), INTENT(IN) :: x(:),x0(:),dx(:),x1(:),Eps(:)
    REAL(KIND=dp), OPTIONAL :: MinMax(:)
    REAL(KIND=dp), INTENT(INOUT) :: val ! Final Elmer point and interpolated value
    INTEGER :: InterpStatus
    !----------------------------------------------------------------------------------------
    INTEGER :: i, ind(3)
    REAL(KIND=dp) :: stencil(2,2,2)
    REAL(KIND=dp) :: Weights(3)
    REAL(KIND=dp) :: xi(3), xf(3)

!      WRITE (*,*) 'X: ', X
!      WRITE (*,*) 'X0: ', X0
!      WRITE (*,*) 'DX: ', DX
!      WRITE (*,*) 'NMAX: ', NMAX
!      WRITE (*,*) 'EPS: ', EPS
!      WRITE (*,*) 'TIME: ', TIMEIndex

    xf = x
    val = 0.0_dp

    DO i = 1,NetDim

      ! Find the (i,j) indices [1,...,max]
      ! Calculates the normalized difference vector;
      ! i.e. the distance/indices to a point x from the leftmost points of the FD box.
      !-------------------------------------------------------------------------------
      ind(i) = CEILING( ( xf(i) - x0(i) ) / dx(i) )

!      PRINT *,'Ind(i): ', i,ind(i),Xf(i),X0(i),Dx(i)

      ! Checks that the estimated index is within the bounding box
      IF( ind(i) < 1 .OR. ind(i) >= DimSize(i) ) THEN

        ! If it's smaller than the leftmost index, but within tolerance (Eps), set it to lower bound; and vice versa
        IF( xf(i) <= x0(i) .AND. xf(i) >= x0(i) - Eps(i) ) THEN
          ind(i) = 1
        ELSE IF( xf(i) >= x1(i) .AND. xf(i) <= x1(i) + Eps(i) ) THEN
          ind(i) = DimSize(i) - 1
        ELSE ! The index is too far to be salvaged
          WRITE (Message, '(A,I0,A,I0,A,F14.3,A)') 'ind(',i,') = ', ind(i), ' from Elmer coordinate ',&
              Xf(i), ' Not in bounding box'
          CALL Info( 'GridDataReader',Message,Level=10)

          InterpStatus = 2
          RETURN
        END IF
      END IF
    END DO

    ! The value of the estimated NetCDF grid point
    xi(1:NetDim) = x0(1:NetDim) + (ind(1:NetDim)-1) * dx(1:NetDim)

    ! Interpolation weights, which are the normalized differences of the estimation
    ! from lower left corner values.
    !------------------------------------------------------------------------------
    weights(1:NetDim) = (xf(1:NetDim)-xi(1:NetDim)) / dx(1:NetDim)

    IF( NetDim == 2 ) THEN
      CALL NetCDFDataCell2D( stencil, Ind, TimeIndex  )

      IF( PRESENT( MinMax ) ) THEN
        ! Too small value
        IF( MINVAL( stencil(1:2,1:2,1) ) < MinMax(1) ) THEN
          InterpStatus = 3
          IF( ALL( stencil(1:2,1:2,1) < MinMax(1) ) ) InterpStatus = 4
          RETURN
        END IF
        ! Too large value
        IF( MAXVAL( stencil(1:2,1:2,1) ) > MinMax(2) ) THEN
          InterpStatus = 5
          IF( ALL( stencil(1:2,1:2,1) > MinMax(2) ) ) InterpStatus = 6
          RETURN
        END IF
      END IF

      val = BiLinearInterpolation(stencil(:,:,1),weights)
    ELSE
      CALL NetCDFDataCell3D( stencil, Ind, TimeIndex  )

      IF( PRESENT( MinMax ) ) THEN
        IF( MINVAL( stencil(1:2,1:2,1:2) ) < MinMax(1) ) THEN
          InterpStatus = 3
          IF( ALL( stencil(1:2,1:2,1:2) < MinMax(1) ) ) InterpStatus = 4
          RETURN
        END IF
        IF( MAXVAL( stencil(1:2,1:2,1:2) ) > MinMax(2) ) THEN
          InterpStatus = 5
          IF( ALL( stencil(1:2,1:2,1:2) > MinMax(2) ) ) InterpStatus = 6
          RETURN
        END IF
      END IF

      val = TriLinearInterpolation(stencil,weights)
    END IF

    ! This is success
    InterpStatus = 1

  END FUNCTION FDInterpolation


  !------------------------------------------------------------------------------------
  ! Transforms Elmer input coordinates into the given coordinate system of the netCDF file.
  !------------------------------------------------------------------------------------
  FUNCTION LocalCoordinateTransformation( vec0, CoordSystem ) RESULT( vec1 )
    !--------------------------------------------------------------
    USE DefUtils, ONLY: dp
    USE Messages
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: CoordSystem  ! Some coordinate transformation
    REAL(KIND=dp), INTENT(IN) :: vec0(3)     ! input coordinates
    REAL(KIND=dp) :: vec1(3)                 ! output coordinates
    REAL(KIND=dp) :: vec2(3)
    REAL(KIND=dp), PARAMETER :: RAD_TO_DEG = 180.0 / PI
    REAL(KIND=dp), PARAMETER :: DEG_TO_RAD = PI / 180.0
    INTEGER :: i

    SELECT CASE ( CoordSystem )

    CASE ('lat-long')
      CALL Warn('GridDataReader','Applies latitude-longitude coordinate transformation; TODO!')
      vec1 = vec0

    CASE ('cylindrical')
      vec1(1) = SQRT( vec0(1)**2 + vec0(2)**2 )
      vec1(2) = RAD_TO_DEG * ATAN2( vec0(2), vec0(1) )
      vec1(3) = vec0(3)

    CASE ('none')
      vec1 = vec0

    CASE DEFAULT
      WRITE (Message,'(A)') 'No coordinate transformation applied: Unknown > Coordinate system < :"'// &
          TRIM(CoordSystem)
      CALL Warn('GridDataReader', Message)
      vec1 = vec0
    END SELECT

  END FUNCTION LocalCoordinateTransformation


!------------------------------------------------------------------------------
END SUBROUTINE GridDataReader
!------------------------------------------------------------------------------
