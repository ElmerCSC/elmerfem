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
! *  Authors: Thomas Zwinger, Martina Schäfer
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Interpolate scattered 2D data read from an ASCII file (x y Value) to the 
!  Elmer mesh.  This code has been checked in because it is thought to be 
!  widely used and therefore safer under version control rather than because it 
!  is thought to be the most suitable code for the purpose.  It may be that the 
!  Scattered2DDataInterpolator (also to be found under the elmerice section of 
!  the repository) is more robust and provides similar functionality, in which 
!  case "pointwise" should be removed from the repository once users have 
!  switched over to Scattered2DDataInterpolator.
!
!  Some revision history follows.
!
!!! version Martina, April 2013
!!! recent bugs fixed
!!! 30.11.2012: bull's eye  theInterpolatedValue = InData(closestPoints(1),DIM+1) (was  theInterpolatedValue = closestPoints(1),DIM+1)
!!! 14.1.2013: variables NoDim and VariableDirections need to be added in the SAVE command
!!! 19.3.2013: corrections of formatting strings so that they compile on taito (A A -> A,A and AI -> A,I)

!!! this version contains more flexibility for looping over the filename over time
!!! 3 possible cases for the filename have been implemented determined by the value of DataI
!!! if you don't want to use that, you have to skip all modifications made in relation to VariableDataName   
!!! usage:
! Solver 1
!  Exec Solver = "Before TimeStep"
!  Equation = "Pointwise Data"
!  Variable =  -nooutput "Dummy"
!  Variable DOFs = 1
!  Procedure = "pointwise-simplified-time" "InterpolatePointValue"
!  Nonlinear System Max Iterations = 1
!  Variable 1 = String "mb" !Variablename in Elmer to be read in

!  Variable DataI 1 = Integer 2
!0-> take only Data1 as filename
!1-> "Data1"+"ceiling(Timestep*Timestepsize)+Data Offset"+Data End
!2-> "Data1"+"ceiling(Timestep*Timestepsize)+Data Offset"+"-"+ceiling(Timestep*Timestepsize)+Data Offset+1" +Data End

!  Variable Data 1 = String "../cmb/cmb_xyz_"  
!  Variable Data Offset 1 = Integer 2005
!  Variable Data End 1 = String ".dat"

!  Variable 1 Supporting Points = Integer 3 !minimum of points to be used for interpolation
!  Variable 1 Dimensions = Integer 2 !dimension of variable, here needs a file with two columns, NOTE numbers have to be written as 0.04 and not 4e-2
!  Variable 1 Exponent = Real 3.0 ! inverse distance weighting exponent (has to be positive) - default 2.0
!  Variable 1 Area Scaling Factor = Real 2.0 ! increases the max search distance by factor 2 of the maximium distance of dataset poitns
!  Variable 1 Directions(2) = Integer 1 2 !which dimensions are these? Here direction 1 and 2 (x and y)
!  Exported Variable 1 = mb !Variablename in Elmer
!  Exported Variable 1 DOFS = Integer 1 !degrees of freedom
!End
!


RECURSIVE SUBROUTINE InterpolatePointValue( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE


  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  INTEGER :: istat, NoVariables, VariableNo, LocalNodes, DataChannel,&
       timearound=1, AllocationIncrement, DIM, dataread, elementnumber,&
       i, j, k, SupportingPoints(99), NoDim(99),VariableDirections(99,3), &
       SimulTime, VariableDataIName(99), VariableDataIName2(99), &
       skippednodes, allnodes
  INTEGER, POINTER :: Permutation(:),VarPerm(:),MaskVarPerm(:),NodeIndexes(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: InData(:,:), DummyIn(:)
  REAL(KIND=dp), POINTER :: Field(:), VarVal(:), MaskVarVal(:)
  LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., Found, GotVar,&
       VariablesExist, reinitiate,Found2,IsMasked(99)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, Name, DataName,DataNameI, &
       DataNameI2, temp,& 
       DataNameEnd,  VariableName(99) , VariableDataName(99), &
       VariableDataNameShort(99), VariableDataNameEnd(99), MaskName(99)
  TYPE(Variable_t), POINTER :: Var, MaskVar
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(ValueList_t), POINTER :: BC, Equation
  TYPE(Element_t), POINTER :: CurrentElement
  REAL (KIND=dp):: theInterpolatedValue, IDExponent(99), ScalingFactor(99)


  SAVE AllocationsDone, FirstTime, SolverName,&
       DIM,NoVariables, VariablesExist, LocalNodes,Permutation,&
       SupportingPoints, VariableName, MaskName,&
       VariableDataNameShort, SimulTime, NoDim, VariableDirections,&
       VariableDataIName, VariableDataNameEnd, VariableDataIName2, &
       IsMasked, IDExponent,ScalingFactor

  IF (FirstTime) WRITE(SolverName,'(A)') 'InterpolatePointValue'
  PointerToSolver => Solver
  IF ( .NOT. ASSOCIATED( PointerToSolver ) ) THEN
     CALL FATAL(SolverName, ' No Solver Pointer associated')
  END IF

  !----------------------------------------
  ! Do these things for the first time only
  !------------------------------------------------------------------------------------------------------------------
  IF (FirstTime) THEN
     LocalNodes = Model % NumberOfNodes
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     CALL INFO(SolverName,'(Re-)Initialization started.',Level=3)

     !-------------------------------------------------------------
     ! Find out how many variables should be read and interpolated
     !-------------------------------------------------------------
     NoVariables = 0
     GotVar = .TRUE.
     DO WHILE(GotVar) 
        NoVariables = NoVariables + 1
        IF (NoVariables > 99) &
             CALL FATAL(SolverName,'Number of parameters cannot exceed 99')
        IF(NoVariables < 10) THEN
           WRITE (Name,'(A,I2)') 'Variable',NoVariables
           WRITE (DataName,'(A,I2)') 'Variable Data',NoVariables
           WRITE (DataNameEnd,'(A,I2)') 'Variable Data End',NoVariables
           WRITE (DataNameI,'(A,I2)') 'Variable DataI',NoVariables
           WRITE (DataNameI2,'(A,I2)') 'Variable Data Offset',NoVariables
        ELSE
           WRITE (Name,'(A,I3)') 'Variable',NoVariables
           WRITE (DataName,'(A,I3)') 'Variable Data',NoVariables
           WRITE (DataNameEnd,'(A,I3)') 'Variable Data End',NoVariables
           WRITE (DataNameI,'(A,I3)') 'Variable DataI',NoVariables
           WRITE (DataNameI2,'(A,I3)') 'Variable Data Offset',NoVariables
        END IF
        VariableName(NoVariables) = ListGetString( Solver % Values, TRIM(Name), GotVar)


        IF(GotVar) THEN
           WRITE(Message,'(A,A,A)') TRIM(Name),': ', VariableName(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
        ELSE
           EXIT
        END IF

        VariableDataNameShort(NoVariables) = ListGetString( Solver % Values, TRIM(DataName), Found)
        VariableDataIName(NoVariables) = ListGetInteger( Solver % Values, TRIM(DataNameI),Found2)
        IF (Found2) THEN
           VariableDataNameEnd(NoVariables) = ListGetString( Solver % Values, TRIM(DataNameEnd), Found)
           VariableDataIName2(NoVariables) = ListGetInteger( Solver % Values, TRIM(DataNameI2),Found)
        ENDIF
        IF (Found) THEN
           WRITE(Message,'(A,A,A,A,A,I2,I6,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Data: ',&
                TRIM(VariableDataNameShort(NoVariables)), &
                VariableDataIName(NoVariables), VariableDataIName2(NoVariables),&
                VariableDataNameEnd(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
        ELSE
           WRITE(Message,'(A,A,A,I5,I5,A,A)') TRIM(Name),' Data file: ',&
                VariableDataNameShort(NoVariables), &
                VariableDataIName(NoVariables),&
                VariableDataIName2(NoVariables), VariableDataNameEnd(NoVariables),' not found.'
           CALL FATAL(SolverName,Message)
        END IF
        ! supporting points used for this variable
        SupportingPoints(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Supporting Points', Found)
        IF (.NOT.Found) THEN
           WRITE(Message,'(A,A)') TRIM(Name),&
                ' Number of supporting points not found - setting to 2'
           CALL WARN(SolverName,Message)
           SupportingPoints(NoVariables) = 2
        ELSE
           IF (SupportingPoints(NoVariables) > 99) THEN
              WRITE(Message,'(A,A,I6,A)') TRIM(Name),&
                   ' Number of supporting points: ', SupportingPoints(NoVariables), &
                   ' exceeds maximum value of 99 - resetting to 99'
              SupportingPoints(NoVariables) = 99
              CALL WARN(SolverName, Message)
           ELSEIF (SupportingPoints(NoVariables) < 2) THEN
              WRITE(Message,'(A,A,I6,A)') TRIM(Name),&
                   ' Number of supporting points: ', SupportingPoints(NoVariables), &
                   ' less than 2 - resetting to 2'
              SupportingPoints(NoVariables) = 2
              CALL WARN(SolverName, Message)
           ELSE
              WRITE(Message,'(A,A,I6)') TRIM(Name),&
                   ' Number of supporting points: ', SupportingPoints(NoVariables)
              CALL INFO(SolverName,Message,Level=4)
           END IF
        END IF
        ! Dimension of coordinates of input data for this variable
        NoDim(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Dimensions', Found)
        VariableDirections(NoVariables,1:3) = 0
        IF (.NOT.Found) THEN
           NoDIM(NoVariables) = DIM
           DO i=1,DIM
              VariableDirections(NoVariables,i) = i
           END DO
        ELSE
           VariableDirections(NoVariables,1:NoDim(NoVariables)) = &
                ListGetIntegerArray( Solver % Values,TRIM(Name) // ' Directions',Found)
           IF (Found) THEN
              WRITE(Message,'(A,A,A,I2,I2,I2,A)')&
                   'Directions for Variable ', TRIM(Name), ': ',VariableDirections(NoVariables,1:3)
              CALL INFO(SolverName,Message,Level=4)
           ELSE
              WRITE(Message,'(A,A,A)') &
                   TRIM(Name) // ' Dimensions', ' found, but no keyword ', TRIM(Name) // ' Directions'
              CALL FATAL(SolverName,Message)
           END IF
        END IF
        ! inverse distance exponent
        IDExponent(NoVariables) = &
             GetConstReal( Solver % Values, TRIM(Name) // ' Exponent', Found)
        IF (.NOT.Found) THEN
          IDExponent(NoVariables) = 2.0_dp
          WRITE(Message,'(A,A,A)')&
                   'Inverse distance >', TRIM(Name),' Exponent< not found - setting to 2.0'
          CALL WARN(SolverName,Message)
        ELSE
           IF (IDExponent(NoVariables) < 0.0_dp) THEN
              WRITE(Message,'(A,A,A,F8.2,A)')&
                   'Inverse distance >', TRIM(Name),' Exponent< value = ', IDExponent(NoVariables),' < 0.0 - setting to 1'
              CALL WARN(SolverName,Message)
              IDExponent(NoVariables) = 1.0_dp
           ELSE
              WRITE(Message,'(A,A,A,F8.2)')&
                   'Exponent for ', TRIM(Name), '= ', IDExponent(NoVariables)
              CALL INFO(SolverName,Message,Level=3)
           END IF
        END IF
        ScalingFactor(NoVariables) = &
             GetConstReal( Solver % Values, TRIM(Name) // ' Area Scaling Factor', Found)
        IF (.NOT.Found) THEN
           ScalingFactor(NoVariables) = 1.0_dp
           WRITE(Message,'(A,A,A)')&
                '>', TRIM(Name),' Area Scaling Factor< not found - setting to 1.0'
           CALL WARN(SolverName,Message)
        ELSE
           WRITE(Message,'(A,A,A,F8.2)')&
               'Area Scaling Factor for ', TRIM(Name), '= ', ScalingFactor(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
        END IF

        MaskName(NoVariables) = &
             ListGetString( Solver % Values, TRIM(Name) // ' Mask', Found)
        
        IF (.NOT.Found) THEN
           
           WRITE(Message,'(A,A,A)')&
                '>', TRIM(Name),' Mask< not found - interpolating whole area'
           CALL INFO(SolverName,Message,Level=3)
           IsMasked(NoVariables)=.FALSE.
        ELSE
           WRITE(Message,'(A,A,A,A)')&
               'Area Scaling Factor for ', TRIM(Name), ' is ', MaskName(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
           IsMasked(NoVariables)=.TRUE.
        END IF

     END DO
     NoVariables = NoVariables-1

     ! --------------------------------------------
     ! Allocate space for new variables to be added
     ! and add them
     ! --------------------------------------------
     IF(NoVariables > 0) VariablesExist = .TRUE.
     IF (VariablesExist) THEN 
        ALLOCATE(Permutation(LocalNodes))
        IF (.NOT.ASSOCIATED(Permutation)) &
             CALL FATAL(Solvername,'Failed to allocate permutation vector\n')
        DO i=1,LocalNodes
           Permutation(i) = i
        END DO
        DO VariableNo = 1, NoVariables
           Var => VariableGet( Model % Variables, TRIM(VariableName(VariableNo)), .TRUE.)     
           IF(.NOT. ASSOCIATED( Var ) ) THEN               
              ALLOCATE(Field(LocalNodes))
              Field(1:LocalNodes) = 1._dp * Permutation(1:LocalNodes) * VariableNo
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                   TRIM(VariableName(VariableNo)), 1, Field, Permutation )          
              WRITE(Message,'(A,I2,A,A,A)') 'Variable no. ',&
                   VariableNo, ' (',TRIM(VariableName(VariableNo)),') added.'
              CALL INFO(SolverName,Message,Level=3)
              NULLIFY( Field ) 
           END IF              
        END DO      
     ELSE
        CALL FATAL(SolverName, 'No valid variables found.')
     END IF
     CALL INFO(Solvername,'(Re-)Initialization done',Level=1)
  END IF ! FirstTime---------------------------------------------------------------------------------------------


  ! -------------------
  ! loop all variables
  ! that need to be
  ! interpolated
  ! -------------------
  IF (VariablesExist) THEN
     SimulTime=ceiling(GetTimeStep()*GetTimeStepSize())
     DO VariableNo = 1,NoVariables
        WRITE(Message,'(A,I3,A,I3)') 'Processing variable ',VariableNo,'/',NoVariables
        CALL INFO(SolverName,Message,Level=3)
        NULLIFY(Var,VarVal,VarPerm)
        Var => VariableGet( Model % Variables, TRIM(VariableName(VariableNo)), .TRUE.) 
        IF (.NOT.ASSOCIATED(Var)) THEN
           WRITE(Message,'(A,A)') 'Variable ', TRIM(VariableName(VariableNo)), ' not associated'
           CALL FATAL(SolverName,Message)
        END IF
        VarPerm  => Var % Perm
        VarVal => Var % Values
        IF (IsMasked(VariableNo)) THEN
           MaskVar => VariableGet( Model % Variables, TRIM(MaskName(VariableNo)), .TRUE.)
           IF (.NOT.ASSOCIATED(MaskVar)) THEN
              WRITE(Message,'(A,A,A)') MaskName(VariableNo),' given but no variable assigned - interpolating over whole area'
              CALL WARN(SolverName,Message)
           END IF
           MaskVarPerm  => MaskVar % Perm
           MaskVarVal => MaskVar % Values
        END IF
        !----------------------
        ! Read data and 
        ! to inquire the size
        !----------------------
        IF (ALLOCATED(DummyIn)) THEN 
           DEALLOCATE(DummyIn)
        END IF
        ALLOCATE(DummyIn(NoDim(VariableNo)))


        IF (VariableDataIName(VariableNo)==0) THEN
           VariableDataName(VariableNo)= VariableDataNameShort(VariableNo)
        ELSEIF (VariableDataIName(VariableNo)==2) THEN
           WRITE(temp,*) VariableDataIName2(VariableNo)+SimulTime
           VariableDataName(VariableNo) = TRIM(VariableDataNameShort(VariableNo))&
                //TRIM(ADJUSTL(temp))//"-"
           WRITE(temp,*) VariableDataIName2(VariableNo)+SimulTime+1
           VariableDataName(VariableNo)=TRIM(VariableDataName(VariableNo))&
                //TRIM(ADJUSTL(temp))//TRIM(VariableDataNameEnd(VariableNo))

        ELSEIF (VariableDataIName(VariableNo)==1) THEN
           WRITE(temp,*) VariableDataIName2(VariableNo)+SimulTime
           VariableDataName(VariableNo) = TRIM(VariableDataNameShort(VariableNo))&
                //TRIM(ADJUSTL(temp))//TRIM(VariableDataNameEnd(VariableNo))

        ELSE
           WRITE(Message,'(A,I2)') "Varialbe DataI has to be 0,1 or 2, but is" , VariableDataIName(VariableNo)
           CALL FATAL(SolverName,Message)
        ENDIF

        OPEN (15, FILE=VariableDataName(VariableNo), STATUS="UNKNOWN", IOSTAT=Istat)
        IF (Istat /= 0) THEN 
           WRITE(Message,'(A,A)') 'Error in opening file ', VariableDataName(VariableNo)
           CALL FATAL(SolverName,Message)
        END IF
        dataread = 0


        DO 
           dataread = dataread + 1     
           READ (15, *, END=10, IOSTAT=Istat, ERR=30) DummyIn(1:NoDim(VariableNo))
        END DO
10      CLOSE(15)       
        dataread = dataread - 1
        WRITE(Message,'(A,I10,A,A)') 'Found ', dataread, ' datasets in ', VariableDataName(VariableNo)

        CALL INFO(SolverName,Message,Level=3)

        IF (ALLOCATED(InData)) &
             DEALLOCATE(InData)

        ALLOCATE(InData(dataread,DIM+1),STAT=Istat)

        IF (istat /= 0) &
             CALL FATAL(SolverName, 'Allocation Error of input data array')
        !----------------------
        ! Read in data and 
        ! interpolate variables
        !----------------------

        OPEN (15, FILE=VariableDataName(VariableNo), STATUS="UNKNOWN", IOSTAT=Istat)

        IF (Istat /= 0) THEN 
           WRITE(Message,'(A,A)') 'Error in opening file ', VariableDataName(VariableNo)
           CALL FATAL(SolverName,Message)
        END IF

        DO i=1,dataread
           READ (15, *, END=20, IOSTAT=Istat, ERR=30) InData(i,1:NoDim(VariableNo)+1)
        END DO
20      CLOSE(15)
        WRITE(Message, '(A,I9,A,A,A,I3,A,A,A)') &
             'Data read  (',dataread,' datasets) from file ', TRIM(VariableDataName(VariableNo)),&
             ' for variable no. ', VariableNo,' (',&
             TRIM(VariableName(VariableNo)),')'
        CALL INFO(SolverName,Message,Level=1)  


        !------------------------------------------------------
        ! interpolate the values
        !------------------------------------------------------
        reinitiate = .TRUE.
        !------------------------------------------------------
        ! Loop all active elements of solver
        !------------------------------------------------------
        skippednodes = 0
        DO elementNumber=1,Solver % NumberOFActiveElements    
           CurrentElement => GetActiveElement(elementNumber)
           IF (.NOT.ASSOCIATED(CurrentElement)) CALL FATAL(SolverName,'Element pointer not associated')
           !------------------------------------------------------
           ! Loop all points of element
           !------------------------------------------------------
           DO i=1,GetElementNOFNodes(CurrentElement)
              allnodes = allnodes + 1
              IF ( (IsMasked(VariableNo)) .AND. (MaskVarVal(MaskVarPerm(CurrentElement % NodeIndexes( i ))) == 0.0_dp) ) THEN
                 skippednodes = skippednodes + 1
                 CYCLE
              END IF
              theInterpolatedValue = GetRadiallyInterpolatedValue(InData,&
                   Model % Nodes % x( CurrentElement % NodeIndexes( i ) ),&
                   Model % Nodes % y( CurrentElement % NodeIndexes( i ) ),&
                   Model % Nodes % z( CurrentElement % NodeIndexes( i ) ),&
                   dataread,&
                   NoDim(VariableNo),&
                   VariableDirections(VariableNo,1:3), &
                   IDExponent(VariableNo), &
                   ScalingFactor(VariableNo), &
                   reinitiate,&
                   SupportingPoints(VariableNo),&
                   SolverName) 
              reinitiate = .FALSE.
              VarVal(VarPerm(CurrentElement % NodeIndexes(i))) = theInterpolatedValue
           END DO
        END DO ! DO elementNumber
        WRITE (Message,'(A,A,A,F5.3)') &
             'Percentage of skipped nodes for variable', TRIM(VariableName(VariableNo)), ': ', FLOAT(skippednodes)/FLOAT(allnodes)
     END DO !DO VariableNo
  END IF
  RETURN

30 CLOSE(15)
  WRITE(Message,'(A,A)') 'Error in reading file ', VariableDataName(VariableNo)
  CALL FATAL(SolverName,Message)

CONTAINS
  !-----------------------------------------------------------------------------------------------------------
  FUNCTION GetRadiallyInterpolatedValue(InData,XI,YI,ZI,&
       NoInData,DIM,Directions,exponent,scalingfactor,reinitiate,numberOfSupportingPoints,SolverName)&
       RESULT(theInterpolatedValue)
    USE DefUtils
    REAL(KIND=dp) :: theInterpolatedValue
    REAL(KIND=dp)  :: InData(:,:), XI,YI,ZI,exponent,scalingfactor
    INTEGER  :: DIM, NoInData, numberOfSupportingPoints, Directions(3)
    LOGICAL :: reinitiate
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName


    REAL(KIND=dp) :: maxdistance, difference, actualdifference, minmaxxy(2,3), &
         radius, weightsum, weight, datasum, X(3)
    REAL(KIND=dp) :: distanceToPoint(99)
    INTEGER :: i,j,k,usedSupportingPoints,actualpoint, istat
    INTEGER :: closestPoints(99)
    LOGICAL :: isSmallerThanAny

    SAVE  maxdistance

    IF (numberOfSupportingPoints < 2) &
         CALL FATAL(SolverName // TRIM('(GetRadiallyInterpolatedValue)'),'The number of supporting points must be at least 2')
    IF (exponent < 0.0) &
         CALL  FATAL(SolverName // TRIM('(GetRadiallyInterpolatedValue)'),'The exponent should be larger or equal to unity')

    X(1:3) = 0.0_dp
    DO i=1,DIM
       SELECT CASE (Directions(i))
       CASE(1)
          X(i) = XI
       CASE(2)
          X(i) = YI
       CASE(3)
          X(i) = ZI
       CASE DEFAULT
          X(i) = 0.0_dp
       END SELECT
    END DO
    !--------------------
    ! (re)initialization
    !--------------------
    IF (reinitiate) THEN
       !-----------------
       ! get bounding box
       !-----------------
       minmaxxy(1,1:DIM) = InData(1,1:DIM)
       minmaxxy(2,1:DIM) = InData(1,1:DIM)
       DO i=1,NoInData
          DO j=1,DIM
             IF (InData(i,j) <  minmaxxy(1,j)) minmaxxy(1,j) = InData(i,j)
             IF (InData(i,j) >  minmaxxy(2,j)) minmaxxy(2,j) = InData(i,j)
          END DO
       END DO
       maxdistance = 0.0_dp
       DO j=1,DIM
          maxdistance = maxdistance + (minmaxxy(2,j) - minmaxxy(1,j))**2.0_dp
       END DO
       maxdistance = scalingfactor * sqrt(maxdistance)     
    END IF
    !-------------------
    ! get the supporting
    ! points for the 
    ! interpolation
    !-------------------
    distanceToPoint = maxdistance
    closestPoints = 0  
    DO i=1,NoInData       
       !-------------------------------
       ! get radius to input data point
       !-------------------------------
       radius = 0.0_dp
       DO j=1,DIM
          radius = radius + (X(j) -  InData(i,j))**2.0_dp
       END DO
       !--------------------------
       ! do we have a bull's eye?
       !--------------------------
       IF (radius > 1.0D-06) THEN 
          radius = SQRT(radius)
       ELSE
          CALL INFO(SolverName, 'Grid and data point coordinates match',Level=12)
          theInterpolatedValue = InData(i,DIM+1)
          RETURN
       END IF
       !----------------------------------------
       ! check whether one and if so which one 
       ! of the current supporting
       ! points has a longer distance
       !----------------------------------------
       actualdifference = distanceToPoint(numberOfSupportingPoints)
       actualpoint = 0
       isSmallerThanAny = .FALSE.
       DO j =1,numberOfSupportingPoints        
          difference = distanceToPoint(j) - radius
          IF (difference > 0.0_dp) THEN
             isSmallerThanAny = .TRUE.
             IF (difference < actualdifference) THEN
                actualpoint = j
                actualdifference = difference
             END IF
          END IF
       END DO
       !-------------
       ! if so, swap
       ! and reorder
       !-------------
       IF (isSmallerThanAny) THEN
          DO k=numberOfSupportingPoints-1,actualpoint,-1
             distanceToPoint(k+1) =  distanceToPoint(k)
             closestPoints(k+1) = closestPoints(k)
          END DO
          distanceToPoint(actualpoint) = radius
          closestPoints(actualpoint) = i
       END IF
    END DO
    !-----------
    ! interpolate
    !-----------
    weightsum = 0.0_dp
    theInterpolatedValue = 0.0_dp
    usedSupportingPoints = 0
    DO k=1,numberOfSupportingPoints
       IF (closestPoints(k) /= 0) THEN
          usedSupportingPoints = usedSupportingPoints + 1
          weight = (distanceToPoint(k))**(-exponent)
          theInterpolatedValue = theInterpolatedValue + weight * InData(closestPoints(k),DIM+1)
          weightsum = weightsum + weight
       END IF
    END DO
    IF (usedSupportingPoints < numberOfSupportingPoints) THEN
       WRITE(Message,'(A,F10.3,F15.3,F15.3,A,I3,A,I3)')&
            'Number of supporting points used for point (',&
            XI,YI,ZI,') =', usedSupportingPoints,&
            ' smaller than requested ', numberOfSupportingPoints
       CALL WARN(TRIM(Solvername) // '(GetRadiallyInterpolatedValue)',&
            Message)
    END IF
    IF (usedSupportingPoints == 0) THEN
       WRITE(Message,'(A,F10.3,F15.3,F15.3,A)') &
            'No supporting point for point (',&
            XI,YI,ZI,') found'
       CALL FATAL(TRIM(Solvername) // '(GetRadiallyInterpolatedValue)',&
            Message)
    END IF
    theInterpolatedValue = theInterpolatedValue/weightsum
    RETURN
  END FUNCTION GetRadiallyInterpolatedValue

END SUBROUTINE InterpolatePointValue


