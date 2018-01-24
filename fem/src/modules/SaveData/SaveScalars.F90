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
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Subroutine for saving scalar data to files
! *
! ******************************************************************************
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20 Nov 2001
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{

SUBROUTINE SaveScalars_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  INTEGER :: NormInd, LineInd, i
  LOGICAL :: GotIt, MarkFailed, AvoidFailed
  CHARACTER(LEN=MAX_NAME_LEN) :: Name

  
  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  NormInd = ListGetInteger( Solver % Values,'Show Norm Index',GotIt)
  IF( NormInd > 0 ) THEN
    Name = ListGetString( Solver % Values, 'Equation',GotIt)
    IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable',&
          '-nooutput -global '//TRIM(Name)//'_var')
    END IF
  END IF

  IF( ParEnv % MyPe == 0 ) THEN
    MarkFailed = ListGetLogical( Solver % Values,'Mark Failed Strategy',GotIt)
    AvoidFailed = ListGetLogical( Solver % Values,'Avoid Failed Strategy',GotIt)
    IF(.NOT. GotIt) AvoidFailed = MarkFailed
    
    IF( MarkFailed .OR. AvoidFailed ) THEN
      LineInd = ListGetInteger( Solver % Values,'Line Marker',GotIt)
      IF(.NOT. GotIt) THEN
        CALL Fatal('SaveScalars_init','Failed strategy marked requires > Line Marker <')
      END IF
      Name = 'FINISHED_MARKER_'//TRIM(I2S(LineInd))
    END IF

    IF( AvoidFailed ) THEN
      INQUIRE(FILE=TRIM(Name),EXIST=GotIt)
      IF( GotIt ) THEN     
        OPEN (10, FILE=Name)
        READ(10,*) i
        IF( i == 0 ) THEN
          CALL Fatal('SaveScalars_init','Strategy already failed before!')
        ELSE
          CLOSE(10)
        END IF
      END IF
    END IF
    
    ! Save a negative status during the execution such that if the
    ! program terminates the negative status will prevail
    IF( MarkFailed ) THEN
      CALL Info('SaveScalars_init','Saving False marker at start')
      i = 0
      OPEN(10,FILE=Name,STATUS='Unknown')
      WRITE(10,'(I0)') i
    END IF
  END IF


  
END SUBROUTINE SaveScalars_init



!------------------------------------------------------------------------------
!> This subroutine saves scalar values in ascii format to an external file.
!> This is a dynamically loaded solver with a standard interface.
!------------------------------------------------------------------------------
SUBROUTINE SaveScalars( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Interpolation

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: ParSolver
  TYPE(ValueList_t), POINTER :: Params
  TYPE(ValueListEntry_t), POINTER :: Lst
  TYPE(Variable_t), POINTER :: Var, OldVar, Var2, Var3
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Nodes_t) :: ElementNodes
  LOGICAL :: MovingMesh, GotCoeff, &
      GotIt, GotOper, GotParOper, GotVar, GotOldVar, ExactCoordinates, VariablesExist, &
      ComplexEigenVectors, ComplexEigenValues, IsParallel, ParallelWrite, LiveGraph, &
      FileAppend, SaveEigenValue, SaveEigenFreq, IsInteger, ParallelReduce, WriteCore, &
      Hit, SaveToFile, EchoValues, GotAny, BodyOper, BodyForceOper, &
      MaterialOper, MaskOper, GotMaskName, GotOldOper, ElementalVar, ComponentVar, &
      NodalOper, GotNodalOper
  LOGICAL, POINTER :: ValuesInteger(:)
  LOGICAL, ALLOCATABLE :: ActiveBC(:)

  REAL (KIND=DP) :: Minimum, Maximum, AbsMinimum, AbsMaximum, &
      Mean, Variance, MinDist, x, y, z, Vol, Intmean, intvar, &
      KineticEnergy, PotentialEnergy, &
      Coords(3), LocalCoords(3), TempCoordinates(3), Val, Val2, &
      Change = 0._dp, Norm = 0.0_dp, PrevNorm, ParallelHits, ParallelCands
  REAL (KIND=DP), ALLOCATABLE :: Values(:), &
      CoordinateBasis(:), ElementValues(:), BoundaryFluxes(:),BoundaryAreas(:)
  REAL (KIND=DP), POINTER :: PointCoordinates(:,:), LineCoordinates(:,:), WrkPntr(:)
  INTEGER, ALLOCATABLE :: BoundaryHits(:)
  INTEGER, POINTER :: PointIndex(:), NodeIndexes(:), CoordinateIndex(:), CoordinatesElemNo(:)
  INTEGER, ALLOCATABLE, TARGET :: ClosestIndex(:)
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: ValueNames(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: ScalarsFile, ScalarNamesFile, DateStr, &
      VariableName, OldVariableName, ResultPrefix, Suffix, Oper, Oper0, OldOper0, ParOper, Name, &
      CoefficientName, ScalarParFile, OutputDirectory, MinOper, MaxOper, &
      MaskName, SaveName
  INTEGER :: i,j,k,l,q,n,ierr,No,NoPoints,NoCoordinates,NoLines,NumberOfVars,&
      NoDims, NoDofs, NoOper, NoElements, NoVar, NoValues, PrevNoValues, DIM, &
      MaxVars, NoEigenValues, Ind, EigenDofs, LineInd, NormInd, CostInd, istat, nlen      
  INTEGER :: IntVal, FirstInd, LastInd 
  LOGICAL, ALLOCATABLE :: NodeMask(:)
  REAL (KIND=DP) :: CT, RT
#ifndef USE_ISO_C_BINDINGS
  REAL (KIND=DP) :: CPUTime, RealTime, CPUMemory
#endif


!------------------------------------------------------------------------------

  CALL Info('SaveScalars', '-----------------------------------------', Level=4 )
  CALL Info('SaveScalars','Saving scalar values of various kinds',Level=4)


  Mesh => GetMesh()
  DIM = CoordinateSystemDimension()
  Params => GetSolverParams()	

  MovingMesh = ListGetLogical(Params,'Moving Mesh',GotIt )

 
  ScalarsFile = ListGetString(Params,'Filename',SaveToFile )
  IF( SaveToFile ) THEN    
    ! Optionally number files by the number of partitions
    ! This makes the benchmarking more convenient since each case 
    ! may use the same command file
    IF(ListGetLogical(Params,'Partition Numbering',GotIt)) THEN
      i = INDEX( ScalarsFile,'.',.TRUE. )
      j = LEN_TRIM(ScalarsFile)
      IF(i > 0) THEN
        Suffix = ScalarsFile(i:j)
        WRITE( ScalarsFile,'(A,I0,A)') &
            ScalarsFile(1:i-1)//'_',ParEnv % PEs,Suffix(1:j-i+1)
      ELSE
        WRITE( ScalarsFile,'(A,I0)') &
            ScalarsFile(1:j)//'_',ParEnv % PEs 
      END IF
    END IF

    IF ( .NOT. FileNameQualified(ScalarsFile) ) THEN
      OutputDirectory = GetString( Params,'Output Directory',GotIt) 
      IF(.NOT. GotIt) OutputDirectory = GetString( Model % Simulation,&
          'Output Directory',GotIt) 
      IF( GotIt .AND. LEN_TRIM(OutputDirectory) > 0 ) THEN
        ScalarsFile = TRIM(OutputDirectory)// '/' //TRIM(ScalarsFile)
        CALL MakeDirectory( TRIM(OutputDirectory) // CHAR(0) )
      ELSE IF( LEN_TRIM(OutputPath ) > 0 ) THEN
        ScalarsFile = TRIM(OutputPath)// '/' //TRIM(ScalarsFile)
      END IF
    END IF

    IF(ListGetLogical(Params,'Filename Numbering',GotIt)) THEN
      ScalarsFile = NextFreeFilename( ScalarsFile ) 
    END IF

    ScalarNamesFile = TRIM(ScalarsFile) // '.' // TRIM("names")
    LiveGraph = ListGetLogical(Params,'Live Graph',GotIt) 
  END IF

  EchoValues = ListGetLogical( Params,'Echo Values',GotIt)
  IF(.NOT. GotIt) EchoValues = .NOT. SaveToFile

  ResultPrefix = ListGetString(Params,'Scalars Prefix',GotIt )
  IF(.NOT. gotIt) ResultPrefix = 'res:'


  IsParallel = .FALSE.
  WriteCore = .TRUE.
  ParallelWrite = .FALSE.
  ParallelReduce = .FALSE.

  IF( ParEnv % PEs > 1 ) THEN
    IsParallel = .TRUE.
    ParallelReduce = GetLogical( Params,'Parallel Reduce',GotIt)
    ParallelWrite = .NOT. ParallelReduce
    IF( ParEnv % MyPe > 0 .AND. ParallelReduce ) THEN
      EchoValues = .FALSE.
      WriteCore = .FALSE. 
    END IF
    OutputPE = ParEnv % MYPe
  END IF

  FileAppend = ListGetLogical( Params,'File Append',GotIt)

  NoLines = 0
  LineCoordinates => ListGetConstRealArray(Params,'Polyline Coordinates',gotIt)
  IF(gotIt) THEN
    NoLines = SIZE(LineCoordinates,1) / 2
    NoDims = SIZE(LineCoordinates,2)
  END IF

  NoPoints = 0
  PointIndex => ListGetIntegerArray( Params,'Save Points',GotIt)
  IF ( gotIt ) NoPoints = SIZE(PointIndex)
    
  NoCoordinates = 0
  NoElements = 0
  PointCoordinates => ListGetConstRealArray(Params,'Save Coordinates',gotIt)
  IF(gotIt) THEN
    NoDims = SIZE(PointCoordinates,2)
    ExactCoordinates = ListGetLogical(Params,'Exact Coordinates',GotIt )      

    IF( ParallelReduce .AND. .NOT. ExactCoordinates) THEN
      CALL Warn('SaveScalars','Only Exact Save Coordinates works in parallel, enforcing...')
      ExactCoordinates = .TRUE.
    END IF

    IF(ExactCoordinates) THEN            
      ! Look for the value at the given coordinate point really.
      NoElements = SIZE(PointCoordinates,1)
      GotIt = .FALSE.
      IF( .NOT. MovingMesh ) THEN
        CoordinatesElemNo => ListGetIntegerArray( Params,'Save Coordinate Elements',GotIt )
      END IF
      IF(.NOT. GotIt ) THEN
        CALL Info('SaveScalars','Searching for elements containing save coordinates',Level=8)
        ALLOCATE(ClosestIndex(NoElements), STAT=istat)
        IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error for CoordinateElemNo')         
        DO j=1,NoElements
          Coords(1:NoDims) = PointCoordinates(j,1:NoDims)
          IF(NoDims < 3 ) Coords(NoDims+1:3) = 0.0_dp
          ClosestIndex(j) = ClosestElementInMesh( Mesh, Coords )
        END DO
        CoordinatesElemNo => ClosestIndex
        IF( .NOT. MovingMesh ) THEN
          CALL ListAddIntegerArray( Params,'Save Coordinate Elements',&
              NoElements,ClosestIndex )
        END IF
      END IF
    ELSE      
      ! Find the indexes of minimum distances
      NoCoordinates = SIZE(PointCoordinates,1)
      GotIt = .FALSE.
      IF( .NOT. MovingMesh ) THEN
        CoordinateIndex => ListGetIntegerArray( Params,'Save Coordinate Indexes',GotIt )
      END IF
      IF( .NOT. GotIt ) THEN
        CALL Info('SaveScalars','Searching for closest nodes to coordinates',Level=8)
        ALLOCATE(ClosestIndex(NoCoordinates), STAT=istat)
        IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error for CoordinateIndex') 
        DO j=1,NoCoordinates 
          Coords(1:NoDims) = PointCoordinates(j,1:NoDims)
          IF(NoDims < 3 ) Coords(NoDims+1:3) = 0.0_dp
          ClosestIndex(j) = ClosestNodeInMesh( Mesh, Coords )
        END DO
        CoordinateIndex => ClosestIndex
        IF( .NOT. MovingMesh ) THEN
          CALL ListAddIntegerArray( Params,'Save Coordinate Indexes',&
              NoCoordinates,ClosestIndex )
        END IF
      END IF
    END IF
  END IF

!------------------------------------------------------------------------------

  n = Mesh % MaxElementNodes 
  ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), &
      ElementValues( n ), CoordinateBasis(n), STAT=istat)
  IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 1') 	

  n = MAX( Model % NumberOfBodies, MAX(Model % NumberOfBCs, NoLines))
  ALLOCATE( BoundaryFluxes(n), BoundaryAreas(n), BoundaryHits(n), STAT=istat )
  IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 2') 	
  
  ALLOCATE( ActiveBC( Model % NumberOfBCs ), STAT=istat )
  IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 3') 	


  ComplexEigenVectors = ListGetLogical(Params,'Complex Eigen Vectors',GotIt)
  
    
   
  !------------------------------------------------------------------------------
  ! Go through the variables and compute the desired statistical data
  !------------------------------------------------------------------------------
  NoValues = 0
  GotVar  = .TRUE.
  GotOper = .FALSE.
  GotOldOper = .FALSE.
  NULLIFY(OldVar)
  NoVar = 0
  MinOper = 'min'
  MaxOper = 'max'
  GotNodalOper = .FALSE.

  DO WHILE(GotVar .OR. GotOper)
    
    GotOper = .FALSE.
    NULLIFY(Var)
    
    NoVar = NoVar + 1
    WRITE (Name,'(A,I0)') 'Variable ',NoVar
  
    VariableName = ListGetString( Params, TRIM(Name), GotVar )


    IF(TRIM(VariableName) == 'cpu time' .OR. TRIM(VariableName) == 'cpu memory') THEN
      CALL Warn('SaveScalars','This variable should now be invoked as an operator: '//TRIM(VariableName))
      CYCLE
    END IF
    
    GotOldVar = .FALSE.

    IF(GotVar) THEN
      Var => VariableGet( Model % Variables, TRIM(VariableName) )
      IF ( .NOT. ASSOCIATED( Var ) )  THEN
        Var => VariableGet( Model % Variables, TRIM(VariableName)//' 1' )
        IF( ASSOCIATED( Var ) ) THEN
          CALL Info('SaveScalars','Treating a component variable: '//TRIM(VariableName),Level=8)
          ComponentVar = .TRUE.
          Var2 => VariableGet( Model % Variables, TRIM(VariableName)//' 2' )
          Var3 => VariableGet( Model % Variables, TRIM(VariableName)//' 3' )          
        ELSE
          CALL Fatal('SaveScalars','Requested variable does not exist: '//TRIM(VariableName))
        END IF
      ELSE
        ComponentVar = .FALSE.
      END IF    
      OldVar => Var
      OldVariableName = VariableName

      ! A 0D variable cannot really be much operated, hence save it as is
      !-------------------------------------------------------------------
      IF(SIZE(Var % Values) == Var % Dofs) THEN
        IsInteger = .FALSE.
        IF( VariableName == 'timestep' ) IsInteger = .TRUE.
        IF( VariableName == 'nonlin iter' ) IsInteger = .TRUE.
        IF( VariableName == 'coupled iter' ) IsInteger = .TRUE.

        IF( Var % Dofs == 1 ) THEN
          CALL AddToSaveList('value: '//TRIM(VariableName)//' scalar variable', &
                              Var % Values(1), IsInteger )
        ELSE
          DO j=1,Var % DOfs
            CALL AddToSaveList('value: '//ComponentName(VariableName,j)//' scalar variable', Var % Values(j))          
          END DO
        END IF
        CYCLE
      END IF

      WRITE (Name,'(A,I0)') 'Nodal Variable ',NoVar
      NodalOper = ListGetLogical(Params,TRIM(Name),GotNodalOper)   
    ELSE
      IF(ASSOCIATED(OldVar)) THEN
        Var => OldVar
        VariableName = OldVariableName
	GotOldVar = .TRUE.
      END IF
    END IF

    NoOper = NoVar     
    MaskOper = .FALSE.
    WRITE (Name,'(A,I0)') 'Operator ',NoOper
    Oper0 = ListGetString(Params,TRIM(Name),GotOper)
    IF(.NOT. GotOper .AND. GotOldOper ) Oper0 = OldOper0

    IF(.NOT. (GotOper .OR. GotVar ) ) CYCLE


    IF( ASSOCIATED( Var ) ) THEN
      CALL Info('SaveScalars','Treating variable: '//TRIM(VariableName),Level=12)
      ElementalVar = ( Var % TYPE == Variable_on_nodes_on_elements ) 
    END IF

    IF( GotOper ) THEN
      CALL Info('SaveScalars','Treating operator: '//TRIM(Oper0),Level=12)
      OldOper0 = Oper0
      GotOldOper = .TRUE.
    ELSE IF( GotOldOper ) THEN
      Oper0 = OldOper0
      GotOper = GotOldOper
    ELSE
      CALL Info('SaveScalars','No operator given for variable: '//TRIM(VariableName))
      CYCLE
    END IF

    
    BodyOper = .FALSE.
    BodyForceOper = .FALSE.      
    MaterialOper = .FALSE.
    nlen = LEN_TRIM(Oper0) 
    IF( Oper0(1:11) == 'body force ') THEN
      BodyForceOper = .TRUE.
      Oper = Oper0(12:nlen)
    ELSE IF( Oper0(1:5) == 'body ') THEN
      BodyOper = .TRUE.
      Oper = Oper0(6:nlen)
    ELSE IF( Oper0(1:9) == 'material ') THEN
      MaterialOper = .TRUE.
      Oper = Oper0(10:nlen)
    ELSE
      Oper = Oper0
    END IF
    MaskOper = ( BodyForceOper .OR. BodyOper .OR. MaterialOper )
    IF( MaskOper ) THEN
      CALL Info('SaveScalars','Operator to be masked: '//TRIM(Oper),Level=12)
    END IF


    WRITE (Name,'(A,I0)') 'Coefficient ',NoOper
    CoefficientName = ListGetString(Params,TRIM(Name),GotCoeff )

    WRITE (Name,'(A,I0)') 'Parallel Operator ',NoOper
    ParOper = ListGetString(Params,TRIM(Name),GotParOper)
    IF(.NOT. GotParOper) ParOper = OperToParOperMap(Oper)

    WRITE (Name,'(A,I0)') 'Mask Name ',NoOper
    MaskName = ListGetString(Params,TRIM(Name),GotMaskName)
    IF(.NOT. GotMaskName) THEN
      MaskName = 'save scalars'
    END IF

    IF( MaskOper ) THEN
      GotIt = .FALSE.
      IF( BodyOper ) THEN
        GotIt = ListGetLogicalAnyBody( Model, MaskName )
      ELSE IF( BodyForceOper ) THEN
        GotIt = ListGetLogicalAnyBodyForce( Model, MaskName )
      ELSE IF( MaterialOper ) THEN
        GotIt = ListGetLogicalAnyMaterial( Model, MaskName )
      END IF
      IF(.NOT. GotIt ) THEN
        CALL Warn('SaveScalars','Masked operators require mask: '//TRIM(MaskName))
      END IF
    END IF
    
    ActiveBC = .FALSE.
    DO j=1,Model % NumberOfBCs
      ActiveBC(j) =  &
          ListGetLogical(Model % BCs(j) % Values,'Flux Integrate',gotIt) .OR. &
          ListGetLogical(Model % BCs(j) % Values, MaskName, gotIt)
    END DO
    
    IF ( GotOper ) THEN
      
      SELECT CASE( Oper ) 
        
      CASE ('partitions')
      CASE ('partition checksum')
      CASE ('partition neighbours checksum')
      CASE ('cpu time')
      CASE ('wall time')
      CASE ('cpu memory') 
      CASE ('nodes')
      CASE ('elements')
      CASE ('bounding box')
      CASE ('area')
      CASE ('volume')
      CASE ('threads')

      CASE DEFAULT
        IF( .NOT. (GotVar .OR. GotOldVar ) ) THEN
          CALL Error('SaveScalars','Operator > '//TRIM(Oper)//' < requires variable')
          CYCLE
        END IF
      END SELECT


      ! Set default name for saving 
      IF(GotVar .OR. GotOldVar ) THEN
        SaveName = TRIM(Oper0)//': '//TRIM(VariableName)
        IF( GotMaskName ) THEN
          SaveName = TRIM(SaveName)//' mask '//TRIM(MaskName)
        END IF
        IF( GotNodalOper ) THEN
          IF( NodalOper ) THEN
            SaveName = TRIM(SaveName)//' nodal'
          ELSE
            SaveName = TRIM(SaveName)//' non-nodal'            
          END IF
        END IF
      END IF
        
      
      SELECT CASE(Oper)

      CASE ('partitions')
        Val = 1.0_dp * ParEnv % PEs 
        SaveName = 'value: number of partitions'
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)

      CASE ('threads')
        Val = 1.0_dp * ParEnv % NumberOfThreads 
        SaveName = 'value: number of threads'
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)

      CASE ('partition checksum')
        Val = 0.0_dp
        IF( ParEnv % PEs > 1 ) THEN
          Val = 1.0_dp * SUM( 1.0_dp * Mesh % ParallelInfo % GlobalDOFS ) 
          ! Give different partition different weight to create something like a checksum
          Val = ( ParEnv % MyPe + 1 ) * Val
        END IF
        SaveName = 'value: partition checksum'
        ! Don't use integer as type because it can exceed the bounds, long should be used
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)

      CASE ('partition neighbours checksum')
        Val = 0.0_dp
        IF( ParEnv % PEs > 1 ) THEN
          DO j=1,Mesh % NumberOfNodes 
            Val = Val + 1.0_dp * SUM( Mesh % ParallelInfo % NeighbourList(j) % Neighbours )
          END DO
          Val = ( ParEnv % MyPe + 1 ) * Val
        END IF
        SaveName = 'value: partition neighbours checksum'
        ! Don't use integer as type because it can exceed the bounds, long should be used
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)

      CASE ('cpu time')
        Val = CPUTime()
        SaveName = 'value: cpu time (s)'
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
      
      CASE ('wall time')
        Val = RealTime()
        SaveName = 'value: real time (s)'
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
      
      CASE ('cpu memory') 
        Val = CPUMemory()
        SaveName = 'value: maximum memory usage (kb)'
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)

      CASE ('nodes')
        Val = 1.0_dp * Solver % Mesh % NumberOfNodes
        SaveName = TRIM(Oper)//': '//TRIM(VariableName)
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)

      CASE ('elements')
        Val = 1.0_dp * Solver % Mesh % NumberOfBulkElements
        SaveName = TRIM(Oper)//': '//TRIM(VariableName)
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)

      CASE ('bounding box')
        Val = MINVAL( Solver % Mesh % Nodes % x ) 
        CALL AddToSaveList(TRIM(Oper)//' min x',Val,ParallelOperator=MinOper)
        Val = MAXVAL( Solver % Mesh % Nodes % x ) 
        CALL AddToSaveList(TRIM(Oper)//' max x',Val,ParallelOperator=MaxOper)
        Val = MINVAL( Solver % Mesh % Nodes % y ) 
        CALL AddToSaveList(TRIM(Oper)//' min y',Val,ParallelOperator=MinOper)
        Val = MAXVAL( Solver % Mesh % Nodes % y ) 
        CALL AddToSaveList(TRIM(Oper)//' max y',Val,ParallelOperator=MaxOper)
        IF( Solver % Mesh % MeshDim > 2 ) THEN
          Val = MINVAL( Solver % Mesh % Nodes % z ) 
          CALL AddToSaveList(TRIM(Oper)//' min z',Val,ParallelOperator=MinOper)
          Val = MAXVAL( Solver % Mesh % Nodes % z )
          CALL AddToSaveList(TRIM(Oper)//' max z',Val,ParallelOperator=MaxOper)
        END IF

      ! Operators that require a variable (except for 'area' and 'volume')
      !-----------------------------------------------------------------
      CASE ('norm')
        Val = Var % Norm
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
        
      CASE ('nonlin change')
        Val = Var % NonlinChange
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
        
      CASE ('steady state change')
        Val = Var % SteadyChange
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
        
      CASE ('nonlin iter')
        Val = Var % NonlinIter
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)
        
      CASE ('nonlin converged')
        Val = 1.0_dp * Var % NonlinConverged
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)
        
      CASE ('steady converged')
        Val = 1.0_dp * Var % SteadyConverged
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)
        
      CASE ('dofs')
        Val = 1.0_dp * SIZE(Var % Values)
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)

      CASE ('nans')  ! number of NaN:s in vector
        GotIt = .FALSE.
        j = 0
        DO i=1,SIZE(Var % Values)
          IF( Var % Values(i) /= Var % Values(i) ) THEN
            j = j + 1
          END IF
        END DO
        Val = 1.0_dp * j
        CALL AddToSaveList(SaveName,Val,.TRUE.,ParOper)
       
      CASE ('sum','sum abs','mean abs','max','max abs','min','min abs','mean','variance','range', &
          'sum square','mean square')
        IF( MaskOper ) CALL CreateNodeMask()
        IF( GotNodalOper ) THEN
          Val = VectorStatistics(Var,Oper,NodalOper)
        ELSE
          Val = VectorStatistics(Var,Oper)
        END IF
        CALL AddToSaveList(SaveName,Val,.FALSE.,ParOper)
        
      CASE ('deviation')
        IF( MaskOper ) CALL CreateNodeMask()
        Val = VectorMeanDeviation(Var,Oper)
        CALL AddToSaveList(SaveName, Val,.FALSE.,ParOper)
        
      CASE ('int','int mean','int abs','int abs mean','int variance','volume',&
          'potential energy','diffusive energy','convective energy')
        
        IF( MaskOper ) CALL CreateNodeMask()
        Val = BulkIntegrals(Var, Oper, GotCoeff, CoefficientName)
        IF(GotCoeff) THEN
          SaveName = TRIM(SaveName)//' with '//TRIM(CoefficientName)
        END IF
        CALL AddToSaveList(SaveName, Val,.FALSE.,ParOper)
        
      CASE('boundary sum','boundary dofs','boundary max','boundary max abs','boundary min',&
          'boundary min abs','boundary mean')
        
        IF( .NOT. ANY( ActiveBC ) ) THEN
          CALL Error('SaveScalars','No flag > '//TRIM(MaskName)// &
              ' < active for operator: '// TRIM(Oper))
          CYCLE
        END IF
 
        BoundaryHits = 0
        BoundaryFluxes = 0.0_dp         
        CALL BoundaryStatistics(Var, Oper, GotCoeff, &
            CoefficientName, BoundaryFluxes, BoundaryHits)
        
        GotAny = .FALSE.
        DO j=1,Model % NumberOfBCs
          IF( ActiveBC(j) ) THEN
            IF( TRIM(Oper) == 'boundary mean' ) THEN
              IF( ParEnv % PEs > 1 .AND. ParallelReduce ) THEN
                CALL Warn('SaveScalars','Operator > boundary mean < not implemented in parallel!')
              ELSE IF( BoundaryHits(j) > 0 ) THEN
                BoundaryFluxes(j) = BoundaryFluxes(j) / BoundaryHits(j)
              END IF
            END IF
            WRITE (Name,'(A,A,A,A,I0)') TRIM(Oper),': ',TRIM(VariableName),' over bc ',j
            CALL AddToSaveList( TRIM(Name), BoundaryFluxes(j),.FALSE.,ParOper)
          END IF
        END DO
        
      CASE ('boundary int','boundary int mean','area','diffusive flux','convective flux')
         
        IF( .NOT. ANY( ActiveBC ) ) THEN
          CALL Error('SaveScalars','No flag > '//TRIM(MaskName)// &
              '< active for operator: '// TRIM(Oper))
        ELSE
          BoundaryHits = 0
          BoundaryFluxes = 0.0_dp
          BoundaryAreas = 0.0_dp      
          Minimum = HUGE(Minimum) 
          Maximum = -HUGE(Maximum)
          
          CALL BoundaryIntegrals(Var, Oper, GotCoeff, CoefficientName,&
              BoundaryFluxes,BoundaryAreas,BoundaryHits)
          
          DO j=1,Model % NumberOfBCs
            IF( ActiveBC(j) ) THEN
              IF( TRIM(Oper) == 'boundary int mean' ) THEN
                IF( ParEnv % PEs > 1 .AND. ParallelReduce ) THEN
                  CALL Warn('SaveScalars','Operator > boundary int mean < not implemented in parallel!')
                ELSE IF( BoundaryAreas(j) > 0.0 ) THEN
                  BoundaryFluxes(j) = BoundaryFluxes(j) / BoundaryAreas(j)
                END IF
              END IF
              WRITE (Name,'(A,A,A,A,I0)') TRIM(Oper),': ',TRIM(VariableName),' over bc ',j
              CALL AddToSaveList( TRIM(Name), BoundaryFluxes(j),.FALSE.,ParOper)
            END IF
          END DO
          
          IF(TRIM(Oper) == 'diffusive flux' .OR. TRIM(Oper) == 'convective flux') THEN
            ParOper = 'min'
            WRITE (Name,'(A,A,A,A)') 'min ',TRIM(Oper),': ',TRIM(VariableName)
            CALL AddToSaveList( TRIM(Name), Minimum,.FALSE.,ParOper)
            
            ParOper = 'max'
            WRITE (Name,'(A,A,A,A)') 'max ',TRIM(Oper),': ',TRIM(VariableName)
            CALL AddToSaveList( TRIM(Name), Maximum,.FALSE.,ParOper)
          END IF
        END IF


        IF(NoLines > 0) THEN
          BoundaryHits = 0
          BoundaryFluxes = 0.0_dp
          BoundaryAreas = 0.0_dp
          
          CALL PolylineIntegrals(Var, Oper, GotCoeff, CoefficientName,&
              BoundaryFluxes,BoundaryAreas, BoundaryHits)
          
          DO j=1,NoLines
            IF( TRIM(Oper) == 'boundary int mean' ) THEN
              IF( BoundaryHits(j) > 0 ) THEN 
                BoundaryFluxes(j) = BoundaryFluxes(j) / BoundaryAreas(j)
              END IF
            END IF
            WRITE (Name,'(A,A,A,A,I0)') TRIM(Oper),': ',TRIM(VariableName),' over polyline ',j
            CALL AddToSaveList( TRIM(Name), BoundaryFluxes(j),.FALSE.,ParOper)
          END DO
        END IF
        
      CASE DEFAULT 
        
        WRITE (Message,'(A,A)') 'Unknown operator: ',TRIM(Oper)
        CALL WARN('SaveScalars',Message)
        
      END SELECT

    END IF
  END DO

  IF( NoVar > 0 ) THEN
    WRITE (Message,'(A,I0,A)') 'Performed ',NoVar-1,' reduction operations'
    CALL Info('SaveScalars',Message,Level=7)
  END IF


  !------------------------------------------------------------------------------
  !   Add eigenvalues on the list if told to
  !------------------------------------------------------------------------------

  SaveEigenValue = ListGetLogical( Params, 'Save Eigenvalues', GotIt )
  IF(.NOT. GotIt) & 
      SaveEigenValue = ListGetLogical( Params, 'Save Eigen values', GotIt )

  SaveEigenFreq = ListGetLogical( Params, 'Save Eigenfrequencies', GotIt )
  IF(.NOT. GotIt) &
      SaveEigenFreq = ListGetLogical( Params, 'Save Eigen Frequencies', GotIt ) 

  IF ( SaveEigenValue .OR. SaveEigenFreq ) THEN
    ComplexEigenValues = ListGetLogical(Params,'Complex Eigen Values',GotIt)
    IF(.NOT. GotIt) &
        ComplexEigenValues = ListGetLogical(Params,'Complex EigenValues',GotIt)
    
    l = 0
    DO i = 1, Model % NumberOfSolvers       
      IF ( Model % Solvers(i) % NOFEigenValues > 0 ) THEN
        DO k = 1, Model % Solvers(i) % NOFEigenValues
          
          Val = REAL( Model % Solvers(i) % Variable % EigenValues(k) )
          l = l + 1

          IF ( SaveEigenValue ) THEN
            WRITE( Name, '("eigen: Re value ", I0)' ) k
            CALL AddToSaveList(TRIM(Name), Val)
            IF( ComplexEigenValues ) THEN
              Val2 = AIMAG( Model % Solvers(i) % Variable % EigenValues(k) )
              WRITE( Name, '("eigen: Im value ", I0)' ) k
              CALL AddToSaveList(TRIM(Name), Val2)
            END IF
          END IF 

          IF ( SaveEigenFreq ) THEN
            WRITE( Name, '("eigen: frequency ", I0, " [Hz]")' ) k
            IF ( Val >= 0.0 ) THEN
              Val = SQRT(Val) / (2*PI)
            ELSE
              ! If the eigenvalue is negative take take take a square root of its absolute value and
              ! return a negative frequency. 
              Val = -SQRT(-Val) / (2*PI)
            END IF
            CALL AddToSaveList(TRIM(Name), Val)
          END IF

        END DO
      END IF
    END DO
    WRITE (Message,'(A,I0,A)') 'Found ',l,' Eigenvalues'
    CALL Info('SaveScalars',Message)
  END IF
    
  !------------------------------------------------------------------------------
  ! Get the info at given node points
  !------------------------------------------------------------------------------
  DO k=1,NoPoints+NoCoordinates

    IF(k <= NoPoints) THEN
      l = PointIndex(k)
    ELSE
      l = CoordinateIndex(k-NoPoints)
    END IF
    
    Var => Model % Variables

    DO WHILE( ASSOCIATED( Var ) )
      
      IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == Var % DOFs ) THEN
        
        CONTINUE

      ELSE IF (ASSOCIATED (Var % EigenVectors)) THEN

        NoEigenValues = SIZE(Var % EigenValues) 
        EigenDofs = SIZE( Var % EigenVectors(1,:) ) / SIZE( Var % Perm )

        IF(EigenDofs == Var % DOFs) THEN
          DO j=1,NoEigenValues
            DO i=1,Var % DOFs              
              Ind = Var % Perm(l)              
              IF( Ind > 0 ) THEN

                Val = REAL( Var % EigenVectors(j, Var%Dofs*(Ind-1)+i) )
                IF(Var % DOFs == 1) THEN
                  WRITE(Name,'("value: Re Eigen ",I0," ",A," at node ",I7)') j,TRIM(Var % Name),l
                ELSE
                  WRITE(Name,'("value: Re Eigen ",I0," ",A,I2," at node ",I7)') j,TRIM(Var % Name),i,l
                END IF
                CALL AddToSaveList( TRIM(Name), Val)
                
                IF(ComplexEigenVectors) THEN
                  Val2 = AIMAG( Var % EigenVectors(j, Var%Dofs*(Ind-1)+i) )
                  IF(Var % DOFs == 1) THEN
                    WRITE(Name,'("value: Im Eigen ",I0," ",A," at node ",I7)') j,TRIM(Var % Name),l
                  ELSE
                    WRITE(Name,'("value: Im Eigen ",I0," ",A,I2," at node ",I7)') j,TRIM(Var % Name),i,l
                  END IF
                  CALL AddToSaveList( TRIM(Name), Val2)
                END IF

              END IF
            END DO
          END DO
        END IF

      ELSE IF( Var % Dofs == 1) THEN          
        ! The variables exist always also as scalars, therefore omit vectors. 

        Ind = l
        IF(ASSOCIATED(Var % Perm)) Ind = Var % Perm(l)

        IF(Ind > 0) THEN
          WRITE(Name,'("value: ",A," at node ",I7)') TRIM( Var % Name ), l
          CALL AddToSaveList( TRIM(Name), Var % Values(Ind))        
        END IF

      END IF

      Var => Var % Next      
    END DO
  END DO

  IF( NoPoints + NoCoordinates > 0 ) THEN
    WRITE (Message,'(A,I0,A)') 'Tabulated all field values at ',NoPoints+NoCoordinates,' points'
    CALL Info('SaveScalars',Message)
  END IF
  
  !------------------------------------------------------------------------------
  ! Get the info at exact coordinates within elements
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! For parallel cases set the default to -HUGE and take the max of the values
  ParOper = 'max'

  DO k=1,NoElements        
    l = CoordinatesElemNo(k)
    IF( l > 0 ) THEN
      CurrentElement => Mesh % Elements(l)
      n = CurrentElement % TYPE % NumberOfNodes

      NodeIndexes => CurrentElement % NodeIndexes

      Coords(1:NoDims) = PointCoordinates(k,1:NoDims)
      IF(NoDims < 3 ) Coords(NoDims+1:3) = 0.0_dp

      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
               
      Hit = PointInElement( CurrentElement, ElementNodes, &
          Coords, LocalCoords, GlobalEps = 1.0_dp, LocalEps=0.1_dp )	          

      ElementValues(1:n) = 0.0d0          
      CoordinateBasis = 0.0_dp
      DO q=1,N
        ElementValues(q) = 1.0d0
        CoordinateBasis(q) = InterpolateInElement( CurrentElement, ElementValues, &
            LocalCoords(1), LocalCoords(2), LocalCoords(3) )
        ElementValues(q) = 0.0d0
      END DO
    ELSE
      IF( .NOT. IsParallel ) CYCLE
    END IF

    Var => Model % Variables
    DO WHILE( ASSOCIATED( Var ) )

      ElementalVar = ( Var % TYPE == Variable_on_nodes_on_elements ) 

      IF ( .NOT. Var % Output .OR. SIZE(Var % Values) == Var % DOFs) THEN
        CONTINUE 
      
      ELSE IF (ASSOCIATED (Var % EigenVectors)) THEN
        NoEigenValues = SIZE(Var % EigenValues) 
        EigenDofs = SIZE( Var % EigenVectors(1,:) ) / SIZE( Var % Perm )

        IF(EigenDofs == Var % DOFs) THEN
          DO j=1,NoEigenValues
            DO i=1,Var % DOFs

              Val = -HUGE(Val)
              Val2 = -HUGE(Val)
              GotIt = .FALSE.

              IF( ElementalVar ) THEN
                NodeIndexes => CurrentElement % DgIndexes
              ELSE
                NodeIndexes => CurrentElement % NodeIndexes 
              END IF

              IF( l == 0 ) THEN
                
              ELSE IF( ALL(Var % Perm(NodeIndexes(1:n)) > 0)) THEN
                ElementValues(1:n) = REAL( Var % EigenVectors(j,Var%Dofs*(Var % Perm(NodeIndexes(1:n))-1)+i) )
                Val = SUM( CoordinateBasis(1:n) * ElementValues(1:n) )             
                
                IF(ComplexEigenVectors) THEN
                  ElementValues(1:n) = AIMAG( Var % EigenVectors(j,Var%Dofs*(Var % Perm(NodeIndexes(1:n))-1)+i) )
                  Val2 = SUM( CoordinateBasis(1:n) * ElementValues(1:n) )             
                END IF

                GotIt = .TRUE.
              END IF
              
              IF( GotIt .OR. IsParallel ) THEN
                IF(Var % DOFs == 1) THEN
                  WRITE(Name,'("value: Re Eigen ",I0," ",A," in element ",I0)') j,TRIM(Var % Name),l
                ELSE
                  WRITE(Name,'("value: Re Eigen ",I0," ",A," ",I0," in element ",I0)') j,TRIM(Var % Name),i,l
                END IF
                CALL AddToSaveList(TRIM(Name), Val,.FALSE.,ParOper)

                IF( ComplexEigenVectors ) THEN
                  IF(Var % DOFs == 1) THEN
                    WRITE(Name,'("value: Im Eigen ",I0," ",A," in element ",I0)') j,TRIM(Var % Name),l
                  ELSE
                    WRITE(Name,'("value: Im Eigen ",I0," ",A," ",I0," in element ",I0)') j,TRIM(Var % Name),i,l
                  END IF
                  CALL AddToSaveList(TRIM(Name), Val2,.FALSE.,ParOper)                  
                END IF
              END IF
            END DO
          END DO
        END IF

      ELSE IF(Var % Dofs == 1) THEN          

        Val = -HUGE( Val )
        GotIt = .FALSE.
        
        IF( ElementalVar ) THEN
          NodeIndexes => CurrentElement % DgIndexes
        ELSE
          NodeIndexes => CurrentElement % NodeIndexes 
        END IF

        IF( l == 0 ) THEN

        ELSE IF( ASSOCIATED( Var % Perm ) ) THEN
          IF( ALL(Var % Perm(NodeIndexes(1:n)) > 0)) THEN            
            ElementValues(1:n) = Var % Values(Var % Perm(NodeIndexes(1:n)))
            Val = SUM( CoordinateBasis(1:n) * ElementValues(1:n) ) 
            GotIt = .TRUE.
          END IF
        ELSE
          ElementValues(1:n) = Var % Values(NodeIndexes(1:n))
          Val = SUM( CoordinateBasis(1:n) * ElementValues(1:n) ) 
          GotIt = .TRUE.
        END IF

        IF(GotIt .OR. IsParallel) THEN
          WRITE(Name,'("value: ",A," in element ",I0)') TRIM(Var % Name),l
          CALL AddToSaveList(TRIM(Name), Val,.FALSE.,ParOper)                   
        END IF
      END IF

      Var => Var % Next      

    END DO
  END DO
  IF( NoElements > 0 ) THEN
    WRITE (Message,'(A,I0,A)') 'Tabulated points within ',NoElements,' elements'
    CALL Info('SaveScalars',Message)
  END IF

  !------------------------------------------------------------------------------
  ! Update time elapsed from start if requested
  !------------------------------------------------------------------------------
  IF( ListGetLogical( Model % Simulation,'Simulation Timing',GotIt ) ) THEN
    CALL Info('SaveScalars','Adding information on simulation timing',Level=10)
    CT = CPUTime() - ListGetConstReal( Model % Simulation,'cputime0',GotIt)
    RT = RealTime() - ListGetConstReal( Model % Simulation,'realtime0',GotIt)
    CALL AddToSaveList('simulation: cpu time (s)',CT )
    CALL AddToSaveList('simulation: real time (s)',RT )
  END IF

  !------------------------------------------------------------------------------
  ! Read the scalars defined in other modules
  !------------------------------------------------------------------------------
  Lst => ListHead(Model % Simulation)
  l = 0
  DO WHILE( ASSOCIATED( Lst ) )    
    IF ( Lst % Name(1:4) == TRIM(ResultPrefix) ) THEN
      CALL AddToSaveList(Lst % Name, Lst % Fvalues(1,1,1))
      l = l + 1
    END IF
    Lst => Lst % Next
  END DO
  IF(l > 0 ) THEN
    WRITE (Message,'(A,I0,A)') 'Found ',l,' result scalars in simulation section'
    CALL Info('SaveScalars',Message,Level=7)
  END IF

  !------------------------------------------------------------------------------
  !   Add free form expressions with GetCReal
  !------------------------------------------------------------------------------
  NoVar = 0
  GotVar = .TRUE.
  DO WHILE( GotVar )
    NoVar = NoVar + 1    
    WRITE (Name,'(A,I0)') 'Expression ',NoVar
    Val = ListGetCReal( Params, TRIM(Name), GotVar )
    IF( GotVar ) CALL AddToSaveList(TRIM(Name),Val)
  END DO

  !------------------------------------------------------------------------------
  ! Add results in Components
  !------------------------------------------------------------------------------
  IF( ListGetLogical( Params,'Save Component Results',GotIt ) ) THEN
    CALL Info('SaveScalars','Saving results from component',Level=10)
    l = 0
    DO i = 1, Model % NumberOfComponents
      Lst => ListHead( Model % Components(i) % Values )
      DO WHILE( ASSOCIATED( Lst ) )    
        IF ( Lst % Name(1:4) == TRIM(ResultPrefix) ) THEN
          CALL AddToSaveList('component '//TRIM(I2S(i))//': '//TRIM(Lst % Name), Lst % Fvalues(1,1,1))
          l = l + 1
        END IF
        Lst => Lst % Next
      END DO
    END DO
    IF(l > 0 ) THEN
      WRITE (Message,'(A,I0,A)') 'Found ',l,' result scalars in components section'
      CALL Info('SaveScalars',Message,Level=7)
    END IF
  END IF
  

  !------------------------------------------------------------------------------
  ! If there are no values 
  !------------------------------------------------------------------------------
  IF( NoValues == 0 ) THEN
    CALL Warn('SaveScalars','Found no values to save')
    RETURN
  ELSE
    CALL Info('SaveScalars','Found '//TRIM(I2S(NoValues))//' values to save in total',Level=6)
  END IF


  
  !------------------------------------------------------------------------------
  ! Finally save all the scalars into a file 
  !------------------------------------------------------------------------------
  IF( SaveToFile ) THEN

    LineInd = ListGetInteger( Params,'Line Marker',GotIt)
    PrevNoValues = ListGetInteger( Params,'Save Scalars Dofs',GotIt) 

    IF(WriteCore .AND. NoValues /= PrevNoValues) THEN 
      CALL ListAddInteger( Params,'Save Scalars Dofs',NoValues )

      WRITE( Message, '(A)' ) 'Saving names of values to file: '//TRIM(ScalarNamesFile)
      CALL Info( 'SaveScalars', Message, Level=4 )
      
      IF( Solver % TimesVisited > 0 ) THEN
        WRITE ( Message,'(A,I0,A,I0)') 'Number of scalar values differ from previous time:',&
            NoValues,' vs. ',PrevNoValues
        CALL Warn('SaveScalars',Message)
      END IF
      
      IF(ParallelWrite) CALL Info('SaveScalars','Parallel data is written into separate files',Level=6)
      IF(ParallelReduce) CALL Info('SaveScalars','Parallel data is reduced into one file',Level=6)
      IF(FileAppend) CALL Info('SaveScalars','Data is appended to existing file',Level=6)
      
      OPEN (10, FILE=ScalarNamesFile)
      Message = ListGetString(Model % Simulation,'Comment',GotIt)
      IF( GotIt ) THEN
        WRITE(10,'(A)') TRIM(Message)
      END IF

      Message = ListGetString(Params,'Comment',GotIt)
      IF( GotIt ) THEN
        WRITE(10,'(A)') TRIM(Message)
      END IF

      DateStr = GetVersion()
      WRITE( 10,'(A)') 'Elmer version: '//TRIM(DateStr)     

      DateStr = GetRevision( GotIt )
      IF( GotIt ) THEN
        WRITE( 10,'(A)') 'Elmer revision: '//TRIM(DateStr)
      END IF        

      DateStr = GetCompilationDate( GotIt )
      IF( GotIt ) THEN
        WRITE( 10,'(A)') 'Elmer compilation date: '//TRIM(DateStr)
      END IF

      DateStr = GetSifName( GotIt ) 
      IF( GotIt ) THEN
        WRITE( 10,'(A)') 'Solver input file: '//TRIM(DateStr)
      END IF
            
      DateStr = FormatDate()      
      WRITE( 10,'(A)') 'File started at: '//TRIM(DateStr)

      WRITE(10,'(A)') ' '
      WRITE(10,'(A)') 'Variables in columns of matrix: '//TRIM(ScalarsFile)
      IF( LineInd /= 0 ) THEN
        i = 1
        WRITE(10,'(I4,": ",A)') 1,'Line Marker'
      ELSE
        i = 0
      END IF
      DO No=1,NoValues 
        WRITE(10,'(I4,": ",A)') No+i,TRIM(ValueNames(No))
      END DO
      CLOSE(10)
    END IF
    
    !------------------------------------------------------------------------------
    ! In parallel case save the data in different files
    !------------------------------------------------------------------------------
    
    WRITE( Message,'(A)' ) 'Saving values to file: '// TRIM(ScalarsFile)
    CALL Info( 'SaveScalars', Message, Level=4 )
    
    IF ( ParallelWrite ) THEN
      WRITE( ScalarParFile, '(A,i0)' ) TRIM(ScalarsFile)//'.', ParEnv % MyPE
      
      IF( Solver % TimesVisited > 0 .OR. FileAppend) THEN 
        OPEN (10, FILE=ScalarParFile,POSITION='APPEND')
      ELSE 
        OPEN (10,FILE=ScalarParFile)
      END IF
    ELSE IF( WriteCore ) THEN 
      IF( Solver % TimesVisited > 0 .OR. FileAppend) THEN 
        OPEN (10, FILE=ScalarsFile,POSITION='APPEND')
      ELSE 
        OPEN (10,FILE=ScalarsFile)
      END IF
    END IF



    IF( WriteCore ) THEN
      ! If there are multiple lines it may be a good idea to mark each by an index
      IF( LineInd /= 0) THEN
        WRITE (10,'(I6)',advance='no') LineInd
      END IF
      DO No=1,NoValues-1
        IF( ValuesInteger(No) ) THEN
          IntVal = NINT( Values(No) )
          WRITE (10,'(I10)',advance='no') IntVal
        ELSE
          WRITE (10,'(ES22.12E3)',advance='no') Values(No)
        END IF
      END DO
      IF( ValuesInteger(NoValues) ) THEN
        IntVal = NINT( Values(No) )
        WRITE (10,'(I10)') IntVal
      ELSE
        WRITE (10,'(ES22.12E3)') Values(NoValues)
      END IF
      CLOSE(10)
    
      !------------------------------------------------------------------------------
      ! Save comments by line in a metadata file
      !------------------------------------------------------------------------------
      IF( Solver % TimesVisited == 0 .AND. FileAppend .AND. LineInd /= 0 ) THEN      
        Message = ListGetString(Params,'Comment',GotIt)
        Name = TRIM(ScalarsFile) // '.' // TRIM("marker")
        IF( GotIt ) THEN
          OPEN (10, FILE=Name,POSITION='APPEND')
          WRITE(10,'(I6,A,A)') LineInd,': ',TRIM(Message)
        END IF
      END IF
    
      !------------------------------------------------------------------------------
      ! Save data in livegraph format
      !------------------------------------------------------------------------------
    
      IF(LiveGraph) THEN
        ! Save data as comma-separated-values (cvs-file)
        IF( Solver % TimesVisited > 0 .OR. FileAppend) THEN 
          OPEN (10, FILE=TRIM(ScalarsFile)//'.csv',POSITION='APPEND')      
        ELSE 
          OPEN (10, FILE=TRIM(ScalarsFile)//'.csv')
          DO No=1,NoValues-1
            WRITE (10,'(A)',advance='no') TRIM(ValueNames(No))//","
          END DO
          WRITE (10,'(A)') TRIM(ValueNames(NoValues))
        END IF
      
        DO No=1,NoValues-1
          WRITE (10,'(ES22.12E3,A)',advance='no') Values(No),","
        END DO
        WRITE (10,'(ES22.12E3)') Values(NoValues)
        CLOSE(10)
      END IF
    END IF   
  END IF ! of SaveFile
  
  !------------------------------------------------------------------------------
  ! Echo values if requested, this is the default if no output to file
  !------------------------------------------------------------------------------
  IF( EchoValues ) THEN
    CALL Info( 'SaveScalars','Showing computed results:')
    DO No=1,NoValues
      WRITE(Message,'(I4,": ",A,ES22.12E3)') No,TRIM(ValueNames(No)),Values(No)
      CALL Info('SaveScalars',Message)
    END DO
  END IF

  !------------------------------------------------------------------------------
  ! Add data to the Simulation block for possible future use
  !------------------------------------------------------------------------------
  DO No=1,NoValues
    CALL ListAddConstReal(Model % Simulation,TRIM(ValueNames(No)),Values(No))
  END DO
 
  !------------------------------------------------------------------------------
  ! For consistancy checks one may print out a value imitating ComputeChange
  !------------------------------------------------------------------------------
  NormInd = ListGetInteger( Params,'Show Norm Index',GotIt)
  IF( NormInd > 0 .AND. NormInd <= NoValues ) THEN
    Norm = Values( NormInd )
    Solver % Variable % Values = Values( NormInd )
    Solver % Variable % Norm = ABS( Values ( NormInd ) ) 

    Name = ListGetString( Params, 'Equation',GotIt)
    IF(.NOT. GotIt) Name = 'SaveScalars'

    ! Here the name is ComputeChange in order to get the change also to ElmerGUI
    ! albeit in a very dirt style. One could also edit ElmerGUI....
    WRITE( Message, '(a,g15.8,g15.8,a)') &
        'SS (ITER=1) (NRM,RELC): (',Norm, Change,&
        ' ) :: '//TRIM( Name )
    CALL Info( 'ComputeChange', Message, Level=3 )
  END IF
  
  !------------------------------------------------------------------------------
  ! One may also export desired data as a cost function for FindOptimum solver
  !------------------------------------------------------------------------------
  CostInd = ListGetInteger( Params,'Cost Function Index',GotIt)
  IF( CostInd > 0 .AND. CostInd <= NoValues ) THEN
    CALL ListAddConstReal( Model % Simulation,'Cost Function',Values(CostInd) )
  END IF


  IF( ParEnv % MyPe == 0 ) THEN
    IF( ListGetLogical( Params,'Mark Failed Strategy',GotIt) ) THEN
      LineInd = ListGetInteger( Params,'Line Marker',GotIt)
      IF(.NOT. GotIt) THEN
        CALL Fatal('SaveScalars','Failed strategy marked requires > Line Marker <')
      END IF

      CALL Info('SaveScalars','Saving True marker at end')
      Name = 'FINISHED_MARKER_'//TRIM(I2S(LineInd))
      i = 1
      OPEN(10,FILE=Name,STATUS='Unknown')
      WRITE(10,'(I0)') i
    END IF
  END IF

  
  IF( NoElements > 0 ) THEN
    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z)
  END IF

  n = 1
  n = NINT( ParallelReduction(1.0_dp*n) )
  
  CALL Info('SaveScalars', '-----------------------------------------', Level=7 )


!------------------------------------------------------------------------------

CONTAINS


  FUNCTION OperToParOperMap( LocalOper ) RESULT ( ParOper )

    CHARACTER(LEN=MAX_NAME_LEN) :: LocalOper, ParOper


    ParOper = 'NONE'
    IF(.NOT. ParallelReduce ) RETURN


    SELECT CASE(LocalOper)
      
    CASE('nodes','elements','dofs','sum','sum abs','int','int abs','volume','potential energy','convective energy',&
        'diffusive energy','boundary sum','boundary dofs','boundary int','area','diffusive flux',&
        'convective flux','nans','partition checksum','partition neighbours checksum')
      ParOper = 'sum'
            
    CASE('max','max abs','boundary max','boundary max abs')
      ParOper = 'max'
      
    CASE('min','min abs','boundary min','boundary min abs')
      ParOper = 'min'

      ! These operators should already be more of less parallel
    CASE('partitions','cpu time','wall time','cpu memory','norm','nonlin change','steady state change',&
	'nonlin iter','nonlin converged','steady converged','threads')
      ParOper = 'none'

    CASE DEFAULT
      ParOper = ListGetString( Params,'Default Parallel Operator',GotIt)
      IF( .NOT. GotIt ) THEN
        ParOper = 'none'
        CALL Warn('SaveScalars','Reduction not implemented in parallel:'//TRIM(LocalOper))
      END IF      

    END SELECT

    CALL Info('OperToParOperMap',TRIM(LocalOper)//' -> '//TRIM(ParOper),Level=12)

  END FUNCTION OperToParOperMap




  SUBROUTINE AddToSaveList(Name, Val, ValueIsInteger, ParallelOperator )

    INTEGER :: n
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp) :: Val, ParVal
    LOGICAL, OPTIONAL :: ValueIsInteger
    CHARACTER(LEN=MAX_NAME_LEN), OPTIONAL :: ParallelOperator
    !------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: Str, ParOper
    REAL, ALLOCATABLE :: TmpValues(:)
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: TmpValueNames(:)
    LOGICAL, POINTER :: TmpValuesInteger(:)     
    TYPE(Variable_t), POINTER :: TargetVar
    INTEGER :: MPIOper
    INTEGER :: istat
    LOGICAL :: GotParOper

    SAVE TmpValues, TmpValueNames


    ! For the first time allocate some space
    IF(.NOT. ALLOCATED(Values)) THEN
      n = 20
      ALLOCATE( Values(n), ValueNames(n), ValuesInteger(n), STAT=istat )
      IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 6') 	
    END IF

    n = NoValues

    ! If vectors are too small make some more room in a rather dummy way
    IF(n >= SIZE(Values) ) THEN
      ALLOCATE(TmpValues(n), TmpValueNames(n), TmpValuesInteger(n),STAT=istat)
      IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 8') 		
	
      TmpValues(1:n) = Values(1:n)
      TmpValuesInteger(1:n) = ValuesInteger(1:n)
      TmpValueNames(1:n) = ValueNames(1:n)
      DEALLOCATE(Values,ValueNames,ValuesInteger)

      ALLOCATE(Values(n+10), ValueNames(n+10), ValuesInteger(n+10), STAT=istat )
      IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 9') 		

      Values(1:n) = TmpValues(1:n)
      ValuesInteger(1:n) = TmpValuesInteger(1:n)
      ValueNames(1:n) = TmpValueNames(1:n)
      DEALLOCATE(TmpValues,TmpValueNames,TmpValuesInteger)
    END IF

    n = n + 1
    Values(n) = Val
    ValueNames(n) = TRIM(Name)
    IF( PRESENT( ValueIsInteger ) ) THEN
      ValuesInteger(n) = ValueIsInteger
    ELSE
      ValuesInteger(n) = .FALSE.
    END IF
    NoValues = n

    IF( ParallelReduce ) THEN
      GotParOper = .FALSE.
      IF( PRESENT(ParallelOperator)) THEN
	ParOper = ParallelOperator
        GotParOper = .TRUE.
      ELSE

!------------------------------------------------------------------------------------------------
! This is to ensure that parallel operators may be applied to also other scalars than those 
! computed within this solver with the serial operators. The indexing may not be known in advance
! but may be checked by a trial run. Note that the indexing does not often coinsice with the 
! normal operators but then this loop is never reached.
!------------------------------------------------------------------------------------------------

        WRITE (Str,'(A,I0)') 'Parallel Operator ',n 
        ParOper = ListGetString(Params,TRIM(Str),GotParOper)
      END IF

      IF( GotParOper ) THEN
        SELECT CASE( ParOper )
          
        CASE('max' )
          MPIOper = MPI_MAX
        CASE('min' )
          MPIOper = MPI_MIN
        CASE('sum' )
          MPIOper = MPI_SUM

        CASE DEFAULT
          MPIOper = 0
          
        END SELECT
        
        IF( MPIOper > 0 ) THEN
          CALL MPI_ALLREDUCE(Val,ParVal,1,&
              MPI_DOUBLE_PRECISION,MPIOper,ELMER_COMM_WORLD,ierr)
          Values(n) = ParVal
          IF( MPIOper == MPI_MIN ) THEN
            WRITE( ValueNames(n),'(A)') TRIM( ValueNames(n) )//' : mpi_min'
          ELSE IF( MPIOper == MPI_MAX ) THEN
            WRITE( ValueNames(n),'(A)') TRIM( ValueNames(n) )//' : mpi_max'
          ELSE IF( MPIOper == MPI_SUM ) THEN
            WRITE( ValueNames(n),'(A)') TRIM( ValueNames(n) )//' : mpi_sum'
          END IF
        END IF
      END IF
    END IF


    !------------------------------------------------------------------------------
    ! If requested, create variable of the result
    ! This is performed already here so the variable can be used 
    ! in subsequent definitions within SaveScalars.
    !------------------------------------------------------------------------------
    WRITE (Str,'(A,I0)') 'Target Variable ',n
    VariableName = ListGetString( Params, TRIM(Str), GotIt )
    IF( GotIt ) THEN
      TargetVar => VariableGet( Model % Variables, TRIM(VariableName) )
      IF(.NOT. ASSOCIATED(TargetVar)) THEN
        NULLIFY(WrkPntr)
        ALLOCATE(WrkPntr(1),STAT=istat)
	IF( istat /= 0 ) CALL Fatal('SaveScalars','Memory allocation error 5') 	
 
        CALL VariableAdd( Model % Variables, Mesh, Solver, &
            TRIM(VariableName), 1, WrkPntr )
        TargetVar => VariableGet( Model % Variables, TRIM(VariableName) )       
      END IF
      TargetVar % Values(1) = Values(n)
      CALL Info('SaveScalars','Defining: '//TRIM(VariableName)//' = '//TRIM(ValueNames(n)),Level=8)
    END IF

  END SUBROUTINE AddToSaveList


  ! Create table for masking nodes in statistical operators.
  !-------------------------------------------------------------------------
  SUBROUTINE CreateNodeMask( ) 

    INTEGER :: t
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: MaskList

    IF(.NOT. MaskOper ) RETURN

    CALL Info('SaveScalars','Creating mask for: '//TRIM(MaskName),Level=10)
    
    IF( ALLOCATED( NodeMask ) ) THEN
      IF( SIZE( NodeMask ) /= SIZE( Var % Perm ) ) DEALLOCATE( NodeMask ) 
    END IF
   
    IF(.NOT. ALLOCATED( NodeMask ) ) THEN
      ALLOCATE( NodeMask( SIZE( Var % Perm ) ) ) 
    END IF
    NodeMask = .FALSE.
    
    DO t = 1, Mesh % NumberOfBulkElements 
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes

      IF( ElementalVar ) THEN
        NodeIndexes => Element % DgIndexes
      ELSE        
        NodeIndexes => Element % NodeIndexes
      END IF

      ! If we are masking operators with correct (body, body force, or material) then do it here
      IF( BodyOper ) THEN        
        MaskList => GetBodyParams( Element ) 
        IF(.NOT. ASSOCIATED(MaskList)) CYCLE
      ELSE IF( BodyForceOper ) THEN
        MaskList => GetBodyForce( Element, GotIt )
        IF( .NOT. GotIt ) CYCLE
      ELSE IF( MaterialOper ) THEN
        MaskList => GetMaterial( Element, GotIt ) 
        IF(.NOT. GotIt ) CYCLE
      ELSE
        CALL Fatal('MaskNodes','Unknown mask strategy')
      END IF
      IF( ListGetLogical( MaskList, MaskName, GotIt ) ) THEN
        NodeMask(NodeIndexes(1:n)) = .TRUE.
      END IF
    END DO
    
    t = COUNT( NodeMask )    
    CALL Info('SaveScalars','Created mask of size: '&
        //TRIM(I2S(t)),Level=12)

  END SUBROUTINE CreateNodeMask




  FUNCTION VectorStatistics(Var,OperName,NodalOper) RESULT (operx)

    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName
    LOGICAL, OPTIONAL :: NodalOper
    REAL(KIND=dp) :: operx

    REAL(KIND=dp) :: Minimum, Maximum, AbsMinimum, AbsMaximum, &
        Mean, Variance, sumx, sumxx, sumabsx, x, Variance2
    INTEGER :: Nonodes, i, j, k, l, NoDofs, sumi
    TYPE(NeighbourList_t), POINTER :: nlist(:)
    INTEGER, POINTER :: PPerm(:)
    
    CALL Info('SaveScalars','Computing operator: '//TRIM(OperName),Level=12)

    sumi = 0
    sumx = 0.0
    sumxx = 0.0
    sumabsx = 0.0

    Maximum = -HUGE(x)
    Minimum = HUGE(x)
    AbsMaximum = -HUGE(x)
    AbsMinimum = HUGE(x)

    PPerm => Var % Perm
    IF( Var % TYPE == Variable_on_gauss_points .OR. &
        Var % TYPE == Variable_on_elements ) THEN
      NULLIFY( PPerm )
    END IF

    
    NoDofs = Var % Dofs
    IF(ASSOCIATED (PPerm)) THEN
      Nonodes = SIZE(PPerm) 
    ELSE
      Nonodes = SIZE(Var % Values) / NoDofs
    END IF

    
    IF( MaskOper ) THEN
      IF( NoNodes > SIZE(NodeMask) ) THEN
        CALL Info('SaveScalars','Decreasing operator range to size of mask: '&
            //TRIM(I2S(NoNodes))//' vs. '//TRIM(I2S(SIZE(NodeMask))), Level=8)
        NoNodes = SIZE(NodeMask)
      END IF
    END IF

    nlist => NULL()
    IF(ParEnv % PEs>1) THEN
      IF(ASSOCIATED(Var % Solver)) THEN
        IF(ASSOCIATED(Var % Solver % Matrix)) THEN
          IF(ASSOCIATED(Var % Solver % Matrix % ParallelInfo)) THEN
            nlist => Var % Solver % Matrix % ParallelInfo % NeighbourList
          END IF
        END IF
      END IF
    END IF

    IF( PRESENT( NodalOper ) ) THEN
      IF( NodalOper ) THEN
        FirstInd = 1
        LastInd = Mesh % NumberOfNodes
      ELSE
        FirstInd = Mesh % NumberOfNodes + 1
        LastInd = SIZE( PPerm ) 
      END IF
    ELSE
      FirstInd = 1
      LastInd = NoNodes
    END IF

    
    DO i=FirstInd,LastInd 
      IF( MaskOper ) THEN
        IF( .NOT. NodeMask(i) ) CYCLE
      END IF

      j = i
      IF(ASSOCIATED(PPerm)) j = Var % Perm(i)

      IF(j > 0) THEN
        IF( ParEnv % PEs > 1 .AND. ASSOCIATED(nlist) ) THEN
          IF( nlist(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
        END IF

        IF(NoDofs <= 1) THEN
          x = Var % Values(j)
        ELSE
          x = 0.0d0
          DO l=1,NoDofs
            x = x + Var % Values(NoDofs*(j-1)+l) ** 2
          END DO
          x = SQRT(x)
        END IF

        sumi = sumi + 1
        sumx = sumx + x
        sumxx = sumxx + x*x
        sumabsx = sumabsx + ABS( x )

        Maximum = MAX(x,Maximum)
        Minimum = MIN(x,Minimum)
        IF(ABS(x) > ABS(AbsMaximum) ) AbsMaximum = x
        IF(ABS(x) < ABS(AbsMinimum) ) AbsMinimum = x
      END IF
    END DO

    ! If there are no dofs avoid division by zero
    IF(sumi == 0) THEN
      operx = 0.0d0
      RETURN
    END IF

    Mean = sumx / sumi

    Variance2 = sumxx/sumi-Mean*Mean
    IF(Variance2 > 0.0d0) THEN
      Variance = SQRT(Variance2) 
    ELSE
      Variance = 0.0d0
    END IF

    SELECT CASE(OperName)
      
    CASE ('sum')
      operx = sumx

    CASE ('sum square')
      operx = sumxx 

    CASE ('sum abs')
      operx = sumabsx
      
    CASE ('max')
      operx = Maximum
      
    CASE ('min')
      operx = Minimum
      
    CASE ('range')
      operx = Maximum - Minimum

    CASE ('max abs')
      operx = AbsMaximum
      
    CASE ('min abs')
      operx = AbsMinimum
      
    CASE ('mean')
      operx = Mean

    CASE ('mean square')
      operx = SQRT( sumxx / sumi )
      
    CASE ('mean abs')
      operx = sumabsx / sumi

    CASE ('variance')
      operx = Variance
      
    CASE DEFAULT 
      CALL Warn('SaveScalars','Unknown statistical operator!')

    END SELECT
      
    
    CALL Info('SaveScalars','Finished computing operator',Level=12)


  END FUNCTION VectorStatistics



  
!------------------------------------------------------------------------------

  FUNCTION VectorMeanDeviation(Var,OperName) RESULT (Deviation)

    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName
    REAL(KIND=dp) :: Mean, Deviation
    REAL(KIND=dp) :: sumx, sumdx, x, dx
    INTEGER :: Nonodes, i, j, k, NoDofs, sumi
    INTEGER, POINTER :: PPerm(:)
    
    NoDofs = Var % Dofs

    PPerm => Var % Perm
    IF( Var % TYPE == Variable_on_gauss_points .OR. &
        Var % TYPE == Variable_on_elements ) THEN
      NULLIFY( PPerm )
    END IF

    IF(ASSOCIATED (PPerm)) THEN
      Nonodes = SIZE(PPerm) 
    ELSE
      Nonodes = SIZE(Var % Values) / NoDofs
    END IF

    IF( MaskOper ) THEN
      IF( NoNodes > SIZE(NodeMask) ) THEN
        CALL Info('SaveScalars','Decreasing operator range to size of mask: '&
            //TRIM(I2S(NoNodes))//' vs. '//TRIM(I2S(SIZE(NodeMask))) )
        NoNodes = SIZE(NodeMask)
      END IF
    END IF

    sumi = 0
    sumx = 0.0
    Deviation = 0.0

    DO i=1,Nonodes
      IF( MaskOper ) THEN
        IF( .NOT. NodeMask(i) ) CYCLE
      END IF

      j = i

      IF( ParEnv % PEs > 1 ) THEN
        IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
      END IF

      IF(ASSOCIATED(PPerm)) j = PPerm(i)
      IF(j > 0) THEN
        IF(NoDofs <= 1) THEN
          x = Var % Values(j)
        ELSE
          x = 0.0d0
          DO k=1,NoDofs
            x = x + Var % Values(NoDofs*(j-1)+k) ** 2
          END DO
          x = SQRT(x)
        END IF
        sumi = sumi + 1
        sumx = sumx + x
      END IF
    END DO

    IF(sumi == 0) RETURN

    Mean = sumx / sumi
    sumi = 0
    sumdx = 0.0
    DO i=1,Nonodes
      j = i
      IF(ASSOCIATED(PPerm)) j = PPerm(i)
      IF(j > 0) THEN
        IF(NoDofs <= 1) THEN
          x = Var % Values(j)
        ELSE
          x = 0.0d0
          DO k=1,NoDofs
            x = x + Var % Values(NoDofs*(j-1)+k) ** 2
          END DO
          x = SQRT(x)
        END IF
        dx = ABS(x-Mean)
        sumi = sumi + 1
        sumdx = sumdx + dx
      END IF
    END DO
    Deviation = sumdx / sumi

  END FUNCTION VectorMeanDeviation

!------------------------------------------------------------------------------

  FUNCTION BulkIntegrals(Var, OperName, GotCoeff, CoeffName) RESULT (operx)
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName, CoeffName
    LOGICAL :: GotCoeff
    REAL(KIND=dp) :: operx, vol
    
    INTEGER :: t, hits
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:), PermIndexes(:)
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes), dBasisdx(Model % MaxElementNodes,3)
    REAL(KIND=dp) :: EnergyTensor(3,3,Model % MaxElementNodes),&
        EnergyCoeff(Model % MaxElementNodes) 
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L,C(3,3),x,y,z
    REAL(KIND=dp) :: func, coeff, integral1, integral2, Grad(3), CoeffGrad(3)
    REAL(KIND=DP), POINTER :: Pwrk(:,:,:)
    LOGICAL :: Stat
    TYPE(ValueList_t), POINTER :: MaskList

    INTEGER :: i,j,k,p,q,DIM,NoDofs
    
    TYPE(GaussIntegrationPoints_t) :: IntegStuff

    hits = 0
    integral1 = 0._dp
    integral2 = 0._dp
    vol = 0._dp
    EnergyCoeff = 1.0d0

    NoDofs = Var % Dofs

    DIM = CoordinateSystemDimension()


    DO t = 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IF(t == Mesh % NumberOfBulkElements + 1 .AND. hits > 0) GOTO 10

      Element => Mesh % Elements(t)
      Model % CurrentElement => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes

      IF ( Element % TYPE % ElementCode == 101 ) CYCLE
      
      IF( ElementalVar ) THEN
        PermIndexes => Element % DgIndexes
      ELSE
        PermIndexes => Element % NodeIndexes
      END IF

      IF ( ANY(Var % Perm(PermIndexes) == 0 ) ) CYCLE      
      hits = hits + 1
      

      NodeIndexes => Element % NodeIndexes 
      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
      

      ! If we are masking operators with correct (body, body force, or material) then do it here
      IF( BodyOper ) THEN        
        MaskList => GetBodyParams( Element ) 
      ELSE IF( BodyForceOper ) THEN
        MaskList => GetBodyForce( Element, GotIt )
        IF( .NOT. GotIt ) CYCLE
      ELSE IF( MaterialOper ) THEN
        MaskList => GetMaterial( Element, GotIt ) 
        IF(.NOT. GotIt ) CYCLE
      END IF
      IF( MaskOper ) THEN
        IF( .NOT. ListGetLogical( MaskList, MaskName, GotIt ) ) CYCLE
      END IF


      k = ListGetInteger( Model % Bodies( Element % BodyId ) % Values, &
          'Material', GotIt, minv=1, maxv=Model % NumberOfMaterials )

      IF( OperName == 'diffusive energy' ) THEN 
        EnergyTensor = 0.0d0
        IF(GotCoeff) THEN
          CALL ListGetRealArray( Model % Materials(k) % Values, &
              TRIM(CoeffName), Pwrk, n, NodeIndexes )

          IF ( SIZE(Pwrk,1) == 1 ) THEN
            DO i=1,3
              EnergyTensor( i,i,1:n ) = Pwrk( 1,1,1:n )
            END DO
          ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
            DO i=1,MIN(3,SIZE(Pwrk,1))
              EnergyTensor(i,i,1:n) = Pwrk(i,1,1:n)
            END DO
          ELSE
            DO i=1,MIN(3,SIZE(Pwrk,1))
              DO j=1,MIN(3,SIZE(Pwrk,2))
                EnergyTensor( i,j,1:n ) = Pwrk(i,j,1:n)
              END DO
            END DO
          END IF
        ELSE 
          DO i=1,3          
            EnergyTensor(i,i,1:n) = 1.0d0
          END DO
        END IF
      ELSE
        k = ListGetInteger( Model % Bodies( Element % BodyId ) % Values, &
            'Material', GotIt, minv=1, maxv=Model % NumberOfMaterials )
        IF(GotCoeff) THEN
          EnergyCoeff = ListGetReal( Model % Materials(k) % Values, &
              TRIM(CoeffName), n, NodeIndexes(1:n) )
        ELSE
          EnergyCoeff(1:n) = 1.0d0
        END IF
      END IF

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Element )
            
      DO i=1,IntegStuff % n
        U = IntegStuff % u(i)
        V = IntegStuff % v(i)
        W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
            Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
        s = 1.0
        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
          y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
          z = SUM( ElementNodes % z(1:n)*Basis(1:n) )
          s = 2*PI
        END IF
        
        CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
        
        s = s * SqrtMetric * SqrtElementMetric * IntegStuff % s(i)
        coeff = SUM( EnergyCoeff(1:n) * Basis(1:n))
        vol =  coeff * vol + S

        SELECT CASE(OperName)
          
          CASE ('volume')
          integral1 = integral1 + coeff * S

          CASE ('int','int mean')
          func = SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) )
          integral1 = integral1 + S * coeff * func 

          CASE ('int abs','int abs mean')
          func = ABS( SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) ) )
          integral1 = integral1 + S * coeff * func 

          CASE ('int variance')
          func = SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) )
          integral1 = integral1 + S * coeff * func 
          integral2 = integral2 + S * coeff * func**2 

          CASE ('diffusive energy')
          CoeffGrad = 0.0d0
          DO j = 1, DIM
            Grad(j) = SUM( dBasisdx(1:n,j) *  Var % Values(Var % Perm(PermIndexes)) )
            DO k = 1, DIM
              CoeffGrad(j) = CoeffGrad(j) + SUM( EnergyTensor(j,k,1:n) * Basis(1:n) ) * &
                  SUM( dBasisdx(1:n,k) * Var % Values(Var % Perm(PermIndexes)) )
            END DO
          END DO
          
          integral1 = integral1 + s * SUM( Grad(1:DIM) * CoeffGrad(1:DIM) )

          CASE ('convective energy')
          func = SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) )

          IF(NoDofs == 1) THEN
            func = SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) )
            integral1 = integral1 + s * coeff * func**2
          ELSE
            func = 0.0d0
            DO j=1,MIN(DIM,NoDofs)
              func = SUM( Var % Values(NoDofs*(Var % Perm(PermIndexes)-1)+j) * Basis(1:n) )
              integral1 = integral1 + s * coeff * func**2
            END DO
          END IF

          CASE ('potential energy')

          func = SUM( Var % Values(Var % Perm(PermIndexes)) * Basis(1:n) )
          integral1 = integral1 + s * coeff * func

        CASE DEFAULT 
          CALL Warn('SaveScalars','Unknown statistical OPERATOR')

        END SELECT

      END DO

    END DO


10  CONTINUE 
    
    operx = 0.0d0
    IF(hits == 0) RETURN

    SELECT CASE(OperName)
      
      CASE ('volume')
      operx = integral1

      CASE ('int')
      operx = integral1

      CASE ('int abs')
      operx = integral1
      
      CASE ('int mean')
      operx = integral1 / vol        

      CASE ('int abs mean')
      operx = integral1 / vol        

      CASE ('int variance')
      operx = SQRT(integral2/vol-(integral1/vol)**2)

      CASE ('diffusive energy')
      operx = 0.5d0 * integral1

      CASE ('convective energy')
      operx = 0.5d0 * integral1

      CASE ('potential energy')
      operx = integral1
      
    END SELECT


  END FUNCTION BulkIntegrals
!------------------------------------------------------------------------------


  SUBROUTINE BoundaryIntegrals(Var, OperName, GotCoeff, &
      CoeffName, fluxes, areas, fluxescomputed)

    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName, CoeffName
    LOGICAL :: GotCoeff
    REAL(KIND=dp) :: fluxes(:), areas(:)
    INTEGER :: fluxescomputed(:)
    
    INTEGER :: t, FluxBody, LBody, RBody, NActive
    TYPE(Element_t), POINTER :: Element, Parent    
    TYPE(ValueList_t), POINTER :: Material, BCVal
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3)
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes),dBasisdx(Model % MaxElementNodes,3),&
        ParentBasis(Model % MaxElementNodes),&
        ParentdBasisdx(Model % MaxElementNodes,3),EnergyTensor(3,3,Model % MaxElementNodes),&
        EnergyCoeff(Model % MaxElementNodes)
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,up,vp,wp,S,A,L,C(3,3),x,y,z
    REAL(KIND=dp) :: func, coeff, Normal(3), Flow(3), flux
    REAL(KIND=DP), POINTER :: Pwrk(:,:,:) => Null()
    INTEGER, POINTER :: ParentIndexes(:), PermIndexes(:)
    REAL(KIND=dp) :: LocalVectorSolution(3,35)

    LOGICAL :: Stat, Permutated    
    INTEGER :: i,j,k,p,q,DIM,bc,NoDofs,pn,hits,istat
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    TYPE(Nodes_t) :: ParentNodes

    n = Model % MaxElementNodes

    DIM = CoordinateSystemDimension()

    ALLOCATE(ParentNodes % x(n), ParentNodes % y(n), ParentNodes % z(n), PermIndexes(n), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('BoundaryIntegrals','Memory allocation error') 


    Minimum = HUGE(minimum)
    Maximum = -HUGE(maximum)

    IF( OperName == 'area' ) THEN
      NoDofs = 1
      Permutated = .FALSE.
    ELSE
      NoDofs = Var % Dofs
      Permutated = ASSOCIATED(Var % Perm)
    END IF

    DO i=1,3          
      EnergyTensor(i,i,:) = 1.0d0
    END DO

    SELECT CASE(OperName)
      
      CASE('diffusive flux') 
      IF(NoDofs /= 1) THEN
        CALL Warn('SaveScalars','diffusive flux & NoDofs /= 1?')
        RETURN
      END IF

      CASE ('convective flux')
      IF( NoDofs /= 1 .AND. NoDofs < DIM) THEN
        CALL Warn('SaveScalars','convective flux & NoDofs < DIM?')
        RETURN
      END IF
      
      CASE ('area','boundary int','boundary int mean')

    CASE DEFAULT 
      CALL Warn('SaveScalars','Unknown statistical OPERATOR')

    END SELECT

    fluxes = 0.0_dp
    areas = 0.0_dp
    coeff = 1.0_dp


    NActive = GetNOFBoundaryElements()

    DO t = 1, NActive 

      Element => GetBoundaryElement(t, Var % Solver)
      BCVal => GetBC()
      IF (.NOT. ASSOCIATED(BCVal) ) CYCLE

      IF (NoDOFs > 1) THEN
        CALL GetVectorLocalSolution(LocalVectorSolution, &
            UElement=Element, USolver=Var % Solver, UVariable=Var)
      ELSE
        CALL GetScalarLocalSolution(LocalVectorSolution(1,:), &
            UElement=Element, USolver=Var % Solver, UVariable=Var)
        IF( ComponentVar ) THEN
          IF( ASSOCIATED( Var2 ) ) THEN
            CALL GetScalarLocalSolution(LocalVectorSolution(2,:), &
                UElement=Element, USolver=Var % Solver, UVariable=Var2)
            NoDofs = 2
          ELSE
            LocalVectorSolution(2,:) = 0.0_dp
          END IF
          IF( ASSOCIATED( Var3 ) ) THEN
            CALL GetScalarLocalSolution(LocalVectorSolution(3,:), &
                UElement=Element, USolver=Var % Solver, UVariable=Var3)
            NoDofs = 3
          ELSE
            LocalVectorSolution(3,:) = 0.0_dp
          END IF
        END IF
      END IF

      Model % CurrentElement => Element

      IF ( Element % TYPE % ElementCode == 101 ) CYCLE

      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      IF(Permutated) THEN
        PermIndexes(1:n) = Var % Perm(NodeIndexes(1:n))
        IF (ANY( PermIndexes(1:n) == 0)) CYCLE
      ELSE
        PermIndexes(1:n) = NodeIndexes(1:n)        
      END IF


      DO bc=1, Model % NumberOfBCs

        IF ( Model % BCs(bc) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE

        IF(.NOT. ListGetLogical(Model % BCs(bc) % Values,'Flux Integrate',gotIt ) .AND. &
            .NOT. ListGetLogical(Model % BCs(bc) % Values, MaskName, gotIt ) ) CYCLE
                 
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))



        SELECT CASE(OperName)
          
        CASE('diffusive flux') 
          
          FluxBody = ListGetInteger( Model % BCs(bc) % Values, &
              'Flux Integrate Body', gotIt ) 
          IF ( GotIt ) THEN
            IF ( ASSOCIATED( Element  % BoundaryInfo % Left ) ) &
                Lbody = Element % BoundaryInfo % Left % BodyId
            IF ( ASSOCIATED( Element  % BoundaryInfo % Right ) ) &
                Rbody = Element % BoundaryInfo % Right % BodyId
            
            IF ( Lbody == FluxBody ) THEN
              Parent => Element % BoundaryInfo % Left
            ELSEIF ( Rbody == FluxBody ) THEN
              Parent => Element % BoundaryInfo % Right
            ELSE
              WRITE( Message, * ) 'No such flux integrate body on bc ', &
                  Element % BoundaryInfo % Constraint
              CALL Fatal( 'SaveScalars', Message )
            END IF
          ELSE        
            Parent => ELement % BoundaryInfo % Left
            stat = ASSOCIATED( Parent )
            
            IF(Permutated) THEN
              IF(stat) stat = ALL(Var % Perm(Parent % NodeIndexes) > 0)
              
              IF ( .NOT. stat ) THEN
                Parent => ELement % BoundaryInfo % Right              
                stat = ASSOCIATED( Parent )
                IF(stat) stat = ALL(Var % Perm(Parent % NodeIndexes) > 0)
              END IF
            END IF
            IF ( .NOT. stat )  CALL Fatal( 'SaveScalars',&
                'No solution available for specified boundary' )
          END IF                   
          
          pn = Parent % TYPE % NumberOfNodes
          ParentIndexes => Parent % NodeIndexes
          i = ListGetInteger( Model % Bodies(Parent % BodyId) % Values, 'Material', &
              minv=1, maxv=Model % NumberOFMaterials )
          Material => Model % Materials(i) % Values

          ParentNodes % x(1:pn) = Mesh % Nodes % x(ParentIndexes(1:pn))
          ParentNodes % y(1:pn) = Mesh % Nodes % y(ParentIndexes(1:pn))
          ParentNodes % z(1:pn) = Mesh % Nodes % z(ParentIndexes(1:pn))

          EnergyTensor = 0._dp

          IF(GotCoeff) THEN
            CALL ListGetRealArray( Material, CoeffName, Pwrk, pn, ParentIndexes )

            IF ( SIZE(Pwrk,1) == 1 ) THEN
              DO i=1,3
                EnergyTensor( i,i,1:pn ) = Pwrk( 1,1,1:pn )
              END DO
            ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
              DO i=1,MIN(3,SIZE(Pwrk,1))
                EnergyTensor(i,i,1:pn) = Pwrk(i,1,1:pn)
              END DO
            ELSE
              DO i=1,MIN(3,SIZE(Pwrk,1))
                DO j=1,MIN(3,SIZE(Pwrk,2))
                  EnergyTensor( i,j,1:pn ) = Pwrk(i,j,1:pn)
                END DO
              END DO
            END IF
          ELSE
            DO i=1,3          
              EnergyTensor(i,i,1:pn) = 1.0d0
            END DO
          END IF
          EnergyCoeff(1:n) = 1.0D00
          fluxescomputed(bc) = fluxescomputed(bc) + 1


        CASE ('convective flux')

          FluxBody = ListGetInteger( Model % BCs(bc) % Values, &
              'Flux Integrate Body', gotIt ) 
          IF ( GotIt ) THEN
            IF ( ASSOCIATED( Element  % BoundaryInfo % Left ) ) &
                Lbody = Element % BoundaryInfo % Left % BodyId
            IF ( ASSOCIATED( Element  % BoundaryInfo % Right ) ) &
                Rbody = Element % BoundaryInfo % Right % BodyId
            
            IF ( Lbody == FluxBody ) THEN
              Parent => Element % BoundaryInfo % Left
            ELSEIF ( Rbody == FluxBody ) THEN
              Parent => Element % BoundaryInfo % Right
            ELSE
              WRITE( Message, * ) 'No such flux integrate body on bc ', &
                  Element % BoundaryInfo % Constraint
              CALL Fatal( 'SaveScalars', Message )
            END IF
          ELSE        
            Parent => ELement % BoundaryInfo % Left
            stat = ASSOCIATED( Parent )
            
            IF(Permutated) THEN
              IF(stat) stat = ALL(Var % Perm(Parent % NodeIndexes) > 0)
              
              IF ( .NOT. stat ) THEN
                Parent => ELement % BoundaryInfo % Right              
                stat = ASSOCIATED( Parent )
                IF(stat) stat = ALL(Var % Perm(Parent % NodeIndexes) > 0)
              END IF
            END IF
            IF ( .NOT. stat )  CALL Fatal( 'SaveScalars',&
                'No solution available for specified boundary' )
          END IF                   
          
          pn = Parent % TYPE % NumberOfNodes
          ParentIndexes => Parent % NodeIndexes
          i = ListGetInteger( Model % Bodies(Parent % BodyId) % Values, 'Material', &
              minv=1, maxv=Model % NumberOFMaterials )
          Material => Model % Materials(i) % Values

          ParentNodes % x(1:pn) = Mesh % Nodes % x(ParentIndexes(1:pn))
          ParentNodes % y(1:pn) = Mesh % Nodes % y(ParentIndexes(1:pn))
          ParentNodes % z(1:pn) = Mesh % Nodes % z(ParentIndexes(1:pn))

          IF(GotCoeff) THEN
            EnergyCoeff(1:n) = ListGetReal( Material, CoeffName, n, NodeIndexes )
          ELSE
            EnergyCoeff(1:n) = 1.0d0
          END IF          
          fluxescomputed(bc) = fluxescomputed(bc) + 1

        CASE ('area','boundary int','boundary int mean')
           IF(GotCoeff) THEN
             i = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                 minv=1, maxv=Model % NumberOFMaterials )
             Material => Model % Materials(i) % Values
             EnergyCoeff(1:n) = ListGetReal( Material, CoeffName, n, NodeIndexes )
          ELSE
            EnergyCoeff(1:n) = 1.0d0
          END IF
          fluxescomputed(bc) = fluxescomputed(bc) + 1


        END SELECT

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------
        IntegStuff = GaussPoints( Element )
        
        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
          x = SUM( ElementNodes % x(1:n)*Basis(1:n) )
          y = SUM( ElementNodes % y(1:n)*Basis(1:n) )
          z = SUM( ElementNodes % z(1:n)*Basis(1:n) )

          s = 1.0d0

          IF(.FALSE.) THEN
            IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
              s = 2.0d0 * PI 
            END IF
            CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
            s = s * SqrtMetric * SqrtElementMetric * IntegStuff % s(i)
          ELSE
            IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
              s = 2.0d0 * PI * x 
            END IF
            s = s * SqrtElementMetric * IntegStuff % s(i)
          END IF
          
          
          SELECT CASE(OperName)
            
          CASE ('diffusive flux')
            
            Normal = NormalVector( Element,ElementNodes,U,V,.TRUE. )
            CALL GetParentUVW( Element, n, Parent, pn, Up, Vp, Wp, Basis)
            
            stat = ElementInfo( Parent,ParentNodes,Up,Vp,Wp,SqrtElementMetric, &
                ParentBasis,ParentdBasisdx )
            
            Flow = 0.0d0
            DO j = 1, DIM
              DO k = 1, DIM
                IF(Permutated) THEN
                  Flow(j) = Flow(j) + SUM( EnergyTensor(j,k,1:pn) * ParentBasis(1:pn) ) * &
                      SUM( ParentdBasisdx(1:pn,k) * Var % Values(Var % Perm(ParentIndexes(1:pn))) )
                ELSE
                  Flow(j) = Flow(j) + SUM( EnergyTensor(j,k,1:pn) * ParentBasis(1:pn) ) * &
                      SUM( ParentdBasisdx(1:pn,k) * Var % Values(ParentIndexes(1:pn)) ) 
                END IF
              END DO
            END DO
            
            flux = SUM(Normal(1:DIM) * Flow(1:DIM))
            Minimum = MIN(flux,Minimum)
            Maximum = MAX(flux,Maximum)
            
            fluxes(bc) = fluxes(bc) + s * flux
            
          CASE ('convective flux')
            
            Normal = NormalVector( Element,ElementNodes,u,v,.TRUE. )
            coeff = SUM( EnergyCoeff(1:n) * Basis(1:n))
            
            IF(NoDofs == 1) THEN
              func = SUM( LocalVectorSolution(1,1:n) * Basis(1:n) )
              fluxes(bc) = fluxes(bc) + s * coeff * func
              Minimum = MIN(Minimum, coeff*func)
              Maximum = MAX(Maximum, coeff*func)
            ELSE 
              DO j=1,DIM
                Flow(j) = coeff * SUM(LocalVectorSolution(j,1:n)*Basis(1:n))
              END DO
              fluxes(bc) = fluxes(bc) + s * SUM(Normal * Flow)
              Minimum = MIN(Minimum,  SUM(Normal * Flow))
              Maximum = MAX(Maximum,  SUM(Normal * Flow))
            END IF
            
          CASE ('boundary int','boundary int mean')
            
            coeff = SUM( EnergyCoeff(1:n) * Basis(1:n))
            
            IF(NoDofs == 1) THEN
              func = SUM( LocalVectorSolution(1,1:n) * Basis(1:n) )
              flux = coeff * func 
              fluxes(bc) = fluxes(bc) + s * flux
            ELSE 
              flux = 0.0_dp
              DO j=1,NoDofs
                flux = flux + SUM( Var % Values(NoDofs*(PermIndexes(1:n)-1)+j) * Basis(1:n) )**2
              END DO
              flux = coeff * SQRT( flux )
              fluxes(bc) = fluxes(bc) + s * flux
            END IF
     
          CASE ('area')
            coeff = SUM( EnergyCoeff(1:n) * Basis(1:n))
            fluxes(bc) = fluxes(bc) + s * coeff 
            
         END SELECT
          
         areas(bc) = areas(bc) + s * coeff

        END DO

      END DO

    END DO

    DEALLOCATE( ParentNodes % x, ParentNodes % y, ParentNodes % z, PermIndexes )

  END SUBROUTINE BoundaryIntegrals


!------------------------------------------------------------------------------
! Take nodal sum over given boundaries
!------------------------------------------------------------------------------

  SUBROUTINE BoundaryStatistics(Var, OperName, GotCoeff, &
      CoeffName, fluxes, fluxescomputed)

    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName, CoeffName
    LOGICAL :: GotCoeff, FindMinMax
    REAL(KIND=dp) :: fluxes(:), val
    INTEGER :: fluxescomputed(:)
    LOGICAL, ALLOCATABLE :: nodescomputed(:)    
    TYPE(Element_t), POINTER :: Element, Parent    
    TYPE(ValueList_t), POINTER :: Material
    LOGICAL :: Stat, Permutated    
    INTEGER :: i,j,k,p,q,t,DIM,bc,n,hits,istat


    ALLOCATE(NodesComputed(SIZE(Var % Perm)),STAT=istat)
    IF( istat /= 0 ) CALL Fatal('BoundaryStatistics','Memory allocation error') 
	
    NodesComputed = .FALSE.
    hits = 0
    
    NoDofs = Var % Dofs
    Permutated = ASSOCIATED(Var % Perm)
    FindMinMax = .FALSE.
    
    SELECT CASE(OperName)
      
    CASE('boundary sum','boundary dofs','boundary mean') 
      fluxes = 0.0_dp

    CASE('boundary min','boundary min abs')
      fluxes = HUGE( val )
      FindMinMax = .TRUE.
      
    CASE('boundary max','boundary max abs') 
      fluxes = -HUGE( val )
      FindMinMax = .TRUE.
      
    CASE DEFAULT 
      CALL Warn('SaveScalars','Unknown statistical operator')
      
    END SELECT


    DO t = Mesh % NumberOfBulkElements+1, Mesh % NumberOfBulkElements &
        + Mesh % NumberOfBoundaryElements

      Element => Mesh % Elements(t)
      Model % CurrentElement => Mesh % Elements(t)

      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      DO bc=1, Model % NumberOfBCs

        IF ( Model % BCs(bc) % Tag /= Element % BoundaryInfo % Constraint ) CYCLE
        IF(.NOT. ListGetLogical(Model % BCs(bc) % Values, MaskName, gotIt ) ) CYCLE

        hits = hits + 1
        
        DO i=1,n
          j = NodeIndexes(i)

          IF( ParEnv % PEs > 1 ) THEN
            IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          END IF

          IF( Permutated ) j = Var % Perm(j)
          
          IF( j == 0) CYCLE            
          IF( nodescomputed(j) ) CYCLE
          
          IF( NoDofs == 1 ) THEN
            val = Var % Values(j)
          ELSE
            val = 0.0_dp
            DO k=1,NoDofs
              val = val + Var % Values(NoDofs*(j-1)+k)
            END DO
            val = SQRT( val )
          END IF
            
          
          IF(FindMinMax) THEN
            SELECT CASE(OperName)              

            CASE('boundary min')
              fluxes(bc) = MIN( val, fluxes(bc) )

            CASE('boundary max')
              fluxes(bc) = MAX( val, fluxes(bc) )

            CASE('boundary min abs')
              IF(ABS(fluxes(bc)) < ABS( val )) fluxes(bc) = val

            CASE('boundary max abs')
              IF(ABS(fluxes(bc)) > ABS( val )) fluxes(bc) = val

            END SELECT
          ELSE
            fluxes(bc) = fluxes(bc) + val
          END IF

          nodescomputed(j) = .TRUE.         
          fluxescomputed(bc) = fluxescomputed(bc) + 1          
        END DO        
      END DO
    END DO


    SELECT CASE(OperName)
      
    CASE('boundary dofs') 
      fluxes = 1.0 * fluxescomputed
      
    END SELECT

  END SUBROUTINE BoundaryStatistics
  

!------------------------------------------------------------------------------


  SUBROUTINE PolylineIntegrals(Var, OperName, GotCoeff, &
      CoeffName, fluxes, areas, fluxescomputed)

    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=MAX_NAME_LEN) :: OperName, CoeffName
    LOGICAL :: GotCoeff
    REAL(KIND=dp) :: fluxes(:),areas(:)
    INTEGER :: fluxescomputed(:)
    
    INTEGER :: t
    TYPE(Element_t), TARGET :: SideElement
    TYPE(Element_t), POINTER :: Element, Parent    
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: SqrtMetric,Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3),LocalCoordinates(3),Point(3)
    REAL(KIND=dp) :: Basis(Model % MaxElementNodes),dBasisdx(Model % MaxElementNodes,3),&
        ParentBasis(Model % MaxElementNodes),&
        ParentdBasisdx(Model % MaxElementNodes,3),EnergyTensor(3,3,Model % MaxElementNodes),&
        EnergyCoeff(Model % MaxElementNodes) 
    REAL(KIND=dp) :: SqrtElementMetric,U,V,W,up,vp,wp,S,A,L,C(3,3),x,y,z,dx,dy,dz,ds,dsmax
    REAL(KIND=dp) :: func, coeff, Normal(3), Flow(3), x0, y0, z0, pos(2), flux
    REAL(KIND=DP), POINTER :: Pwrk(:,:,:)
    INTEGER, POINTER :: ParentIndexes(:), PermIndexes(:), SideIndexes(:), OnLine(:,:)

    LOGICAL :: Stat, Permutated, Inside    
    INTEGER :: i,j,k,p,q,DIM,bc,NoDofs,pn,Line, NoSides, Side, NodeNumber, LineNode(2), istat
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    TYPE(Nodes_t) :: ParentNodes, LineNodes, SideNodes


    n = Model % MaxElementNodes
    DIM = CoordinateSystemDimension()

    ALLOCATE(ParentNodes % x(n), ParentNodes % y(n), ParentNodes % z(n), PermIndexes(n), &
	SideNodes % x(n), SideNodes % y(n), SideNodes % z(n), SideIndexes(n), &
	LineNodes % x(n), LineNodes % y(n), LineNodes % z(n), &
	OnLine(Mesh % NumberOfNodes,2), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('PolylineIntegrals','Memory allocation error') 


    CALL Info('PolylineIntegrals','OPERATOR: '//TRIM(OperName))

    areas = 0.0_dp
    coeff = 1.0_dp   

    NoDofs = Var % Dofs
    IF( OperName == 'area' ) THEN
      Permutated = .FALSE.
    ELSE
      Permutated = ASSOCIATED(Var % Perm)
    END IF

    DO i=1,3          
      EnergyTensor(i,i,:) = 1.0d0
    END DO

    SELECT CASE(OperName)
      
      CASE('diffusive flux') 
      IF(NoDofs /= 1) THEN
        CALL Warn('SaveScalars','diffusive flux & NoDofs /= 1?')
        RETURN
      END IF

      CASE ('convective flux')
      IF(NoDofs /= 1 .AND. NoDofs < DIM) THEN
        CALL Warn('SaveScalars','convective flux & NoDofs < DIM?')
        RETURN
      END IF
      
      CASE ('area','boundary int','boundary int mean')

    CASE DEFAULT 
      CALL Warn('SaveScalars','Unknown physical OPERATOR')

    END SELECT

    DIM = CoordinateSystemDimension()
    IF(DIM < 3) THEN
      LineNodes % z = 0.0d0
      SideNodes % z(1:2) = 0.0d0
    END IF

    fluxes = 0.0d0
    SideElement % TYPE => GetElementType( 202, .FALSE.)
    SideElement % Bdofs = 0
    Element => SideElement

!   /* Go through the line segments */
    DO Line = 1,NoLines 

      LineNodes % x(1:2) = LineCoordinates(2*Line-1:2*Line,1) 
      LineNodes % y(1:2) = LineCoordinates(2*Line-1:2*Line,2) 
      IF(DIM == 3) LineNodes % z(1:2) = LineCoordinates(2*Line-1:2*Line,3) 
      OnLine = 0

      DO t = 1, Mesh % NumberOfBulkElements 
        
        Parent => Mesh % Elements(t)
        Model % CurrentElement => Mesh % Elements(t)

        NoSides = Parent % TYPE % ElementCode / 100  
        IF(NoSides < 3 .OR. NoSides > 4) CYCLE
        
        pn = Parent % TYPE % NumberOfNodes
        ParentIndexes => Parent % NodeIndexes

        ParentNodes % x(1:pn) = Mesh % Nodes % x(ParentIndexes(1:pn))
        ParentNodes % y(1:pn) = Mesh % Nodes % y(ParentIndexes(1:pn))
        ParentNodes % z(1:pn) = Mesh % Nodes % z(ParentIndexes(1:pn))          
        
        IF(Permutated) THEN
          PermIndexes(1:pn) = Var % Perm(ParentIndexes(1:pn))
          IF (ANY( PermIndexes(1:pn) == 0)) CYCLE
        ELSE
          PermIndexes(1:pn) = ParentIndexes(1:pn)        
        END IF

        NodeNumber = 0

        DO Side = 1, NoSides

          SideIndexes(1) = ParentIndexes(Side)
          SideIndexes(2) = ParentIndexes(MOD(Side,NoSides)+1)
          
          SideNodes % x(1:2) = Mesh % Nodes % x(SideIndexes(1:2))
          SideNodes % y(1:2) = Mesh % Nodes % y(SideIndexes(1:2))
          IF(DIM == 3) SideNodes % z(1:2) = Mesh % Nodes % z(SideIndexes(1:2))
  
          CALL LineIntersectionCoords(SideNodes,LineNodes,Inside,x0,y0,z0,u)

          IF(.NOT. Inside) CYCLE

          NodeNumber = NodeNumber + 1        
          ElementNodes % x(NodeNumber) = x0
          ElementNodes % y(NodeNumber) = y0
          ElementNodes % z(NodeNumber) = z0
          pos(NodeNumber) = u

        END DO

        IF(NodeNumber == 0) CYCLE

        IF(NodeNumber > 2) THEN
          CALL Warn('PolylineIntergrals','There should not be more than 2 intersections!')
          CYCLE
        END IF

        !---------------------------------------------------------------------------
        ! If there is only one intersection the other end of the node must lie
        ! inside the element. Assuming that the line is long compared to the 
        ! element the correct end of the line segment may be easily deduced.
        !---------------------------------------------------------------------------
        IF(NodeNumber == 1) THEN
          IF(pos(1) < 0.5d0) THEN
            i = 1
            pos(2) = 0.0
          ELSE
            i = 2
            pos(2) = 1.0
          END IF
          x0 = LineNodes % x(i)
          y0 = LineNodes % y(i)
          z0 = LineNodes % z(i)            

          ElementNodes % x(2) = LineNodes % x(i)
          ElementNodes % y(2) = LineNodes % y(i)
          ElementNodes % z(2) = LineNodes % z(i)            
        END IF

        IF(ABS(pos(1)-pos(2)) < 1.0d-8) CYCLE

        !-----------------------------------------------------------------------------
        ! Change the order of nodes so that the normal always points to the same direction          
        !-----------------------------------------------------------------------------
        IF(pos(1) < pos(2)) THEN
          ElementNodes % x(2) = ElementNodes % x(1)
          ElementNodes % y(2) = ElementNodes % y(1)
          ElementNodes % z(2) = ElementNodes % z(1)
          ElementNodes % x(1) = x0
          ElementNodes % y(1) = y0
          ElementNodes % z(1) = z0           
        END IF
        
        !--------------------------------------------------------------------------------
        ! The following avoids the cases where the line goes exactly at the element 
        ! interface and therefore the flux would be computed twice
        !--------------------------------------------------------------------------------
        dx = ElementNodes % x(1) - ElementNodes % x(2)
        dy = ElementNodes % y(1) - ElementNodes % y(2)
        dsmax = SQRT(dx*dx+dy*dy)
        LineNode = 0


        DO i=1,Parent % TYPE % ElementCode / 100 
          DO j=1,2
            dx = ParentNodes % x(i) - ElementNodes % x(j)
            dy = ParentNodes % y(i) - ElementNodes % y(j)
            ds = SQRT(dx*dx+dy*dy)
            IF(ds < 1.0d-4 * dsmax) LineNode(j) = ParentIndexes(i)
          END DO
        END DO

        IF(ALL(LineNode(1:2) > 0)) THEN
          IF(ANY(OnLine(LineNode(1),:) == LineNode(2))) CYCLE
          IF(ANY(OnLine(LineNode(2),:) == LineNode(1))) CYCLE

          IF(OnLine(LineNode(1),1) == 0) THEN
            OnLine(LineNode(1),1) = LineNode(2)
          ELSE IF(OnLine(LineNode(1),2) == 0) THEN
            OnLine(LineNode(1),2) = LineNode(2)
          ELSE
            CALL Warn('PolylineIntegrate','This should never happen')
          END IF
          
          IF(OnLine(LineNode(2),1) == 0) THEN
            OnLine(LineNode(2),1) = LineNode(1)
          ELSE IF(OnLine(LineNode(2),2) == 0) THEN
            OnLine(LineNode(2),2) = LineNode(1)
          ELSE
            CALL Warn('PolylineIntegrate','This should never happen')
          END IF
        END IF


        
        i = ListGetInteger( Model % Bodies(Parent % BodyId) % Values, 'Material', &
            minv=1, maxv=Model % NumberOFMaterials )
        Material => Model % Materials(i) % Values
        fluxescomputed(Line) = fluxescomputed(Line) + 1
        
        
        SELECT CASE(OperName)
          
          CASE('diffusive flux') 
          
          IF(GotCoeff) THEN
            CALL ListGetRealArray( Material, CoeffName, Pwrk, &
                pn, ParentIndexes, gotIt )
            
            EnergyTensor = 0._dp
            IF(GotIt) THEN
              IF ( SIZE(Pwrk,1) == 1 ) THEN
                DO i=1,3
                  EnergyTensor( i,i,1:pn ) = Pwrk( 1,1,1:pn )
                END DO
              ELSE IF ( SIZE(Pwrk,2) == 1 ) THEN
                DO i=1,MIN(3,SIZE(Pwrk,1))
                  EnergyTensor(i,i,1:pn) = Pwrk(i,1,1:pn)
                END DO
              ELSE
                DO i=1,MIN(3,SIZE(Pwrk,1))
                  DO j=1,MIN(3,SIZE(Pwrk,2))
                    EnergyTensor( i,j,1:pn ) = Pwrk(i,j,1:pn)
                  END DO
                END DO
              END IF
            ELSE 
              DO i=1,3          
                EnergyTensor(i,i,1:pn) = 1.0d0
              END DO
            END IF
          END IF
          
          CASE ('convective flux','area')
          EnergyCoeff(1:n) = ListGetReal( Material, CoeffName, pn, ParentIndexes, gotIt )
          IF(.NOT. GotIt) EnergyCoeff(1:pn) = 1.0d0
          
        END SELECT

!------------------------------------------------------------------------------
!    Numerical integration
!------------------------------------------------------------------------------

        IntegStuff = GaussPoints( Element, 2 )
        
        DO i=1,IntegStuff % n
          U = IntegStuff % u(i)
          V = IntegStuff % v(i)
          W = IntegStuff % w(i)

!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
          stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
              Basis,dBasisdx )
          
          x = SUM( ElementNodes % x(1:n) * Basis(1:n) )
          y = SUM( ElementNodes % y(1:n) * Basis(1:n) )
          z = SUM( ElementNodes % z(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------

          s = 1.0d0

          IF(.FALSE.) THEN
            IF(CurrentCoordinateSystem() /= Cartesian ) s = 2.0 * PI 
            CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,x,y,z )
            s = s * SqrtMetric * SqrtElementMetric * IntegStuff % s(i)
          ELSE
            IF(CurrentCoordinateSystem() /= Cartesian ) s = 2.0 * PI * x
            s = s * SqrtElementMetric * IntegStuff % s(i)            
          END IF
     
          Normal = NormalVector( Element,ElementNodes,U,V,.FALSE. )

!------------------------------------------------------------------------------
!      Because the intersection nodes do not really exist the field variables 
!      must be avaluated using the nodes of the parent element.
!------------------------------------------------------------------------------
            
          Point(1) = x
          Point(2) = y
          Point(3) = z
          
          IF( .NOT. PointInElement( Parent, ParentNodes, Point, LocalCoordinates ) ) THEN
            CALL Warn('PolylineIntegrals','The node should be in the element by construction!')
            CYCLE
          END IF
          
          Up = LocalCoordinates(1)
          Vp = LocalCoordinates(2)
          Wp = LocalCoordinates(3)
          
          
          stat = ElementInfo( Parent,ParentNodes,Up,Vp,Wp,SqrtElementMetric, &
              ParentBasis,ParentdBasisdx )
          
  
          SELECT CASE(OperName)
            
            CASE ('diffusive flux')           
              Flow = 0.0d0
              DO j = 1, DIM
                DO k = 1, DIM
                  Flow(j) = Flow(j) + SUM( EnergyTensor(j,k,1:pn) * ParentBasis(1:pn) ) * &
                      SUM( ParentdBasisdx(1:pn,k) * Var % Values(PermIndexes(1:pn)) ) 
                END DO
              END DO
              fluxes(Line) = fluxes(Line) + s * SUM(Normal(1:DIM) * Flow(1:DIM))
            
            
            CASE ('convective flux')            
              coeff = SUM( EnergyCoeff(1:pn) * ParentBasis(1:pn))              
              IF(NoDofs == 1) THEN
                func = SUM( Var % Values(PermIndexes(1:pn)) * ParentBasis(1:pn) )                
                fluxes(Line) = fluxes(Line) + s * coeff * func
              ELSE 
                DO j=1,DIM
                  Flow(j) = coeff * &
                      SUM( Var % Values(NoDofs*(PermIndexes(1:pn)-1)+j) * ParentBasis(1:pn) )
                END DO
                fluxes(Line) = fluxes(Line) + s * coeff * SUM(Normal * Flow)
              END IF
 
           CASE ('boundary int','boundary int mean')            
              coeff = SUM( EnergyCoeff(1:pn) * ParentBasis(1:pn))              
              IF(NoDofs == 1) THEN
                func = SUM( Var % Values(PermIndexes(1:pn)) * ParentBasis(1:pn) )                
                fluxes(Line) = fluxes(Line) + s * coeff * func
              ELSE 
                flux = 0.0_dp
                DO j=1,NoDofs
                  flux = flux + &
                      SUM( Var % Values(NoDofs*(PermIndexes(1:pn)-1)+j) * ParentBasis(1:pn) )**2
                END DO
                flux = coeff * SQRT(flux) 
                fluxes(Line) = fluxes(Line) + s * flux
              END IF
                         
            CASE ('area')                        
              coeff = SUM( EnergyCoeff(1:pn) * Basis(1:pn))
              fluxes(Line) = fluxes(Line) + s * coeff 
              
            END SELECT

            areas(Line) = areas(Line) + s * coeff

        END DO

      END DO

    END DO

    DEALLOCATE( ParentNodes % x, ParentNodes % y, ParentNodes % z, PermIndexes, &
        SideNodes % x, SideNodes % y, SideNodes % z, SideIndexes, &
        LineNodes % x, LineNodes % y, LineNodes % z, OnLine)

  END SUBROUTINE PolylineIntegrals
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> This subroutine tests whether the line segment goes through the current
!> face of the element. 
!------------------------------------------------------------------------------
  SUBROUTINE LineIntersectionCoords(Plane,Line,Inside,x0,y0,z0,frac)

    TYPE(Nodes_t) :: Plane, Line
    LOGICAL :: Inside
    REAL (KIND=dp) :: x0, y0, z0, frac

    REAL (KIND=dp) :: A(3,3),B(3),C(3),eps=1.0d-6,detA,absA

    Inside = .FALSE.
    
    ! In 2D the intersection is between two lines
    A(1,1) = Line % x(2) - Line % x(1)
    A(2,1) = Line % y(2) - Line % y(1)
    A(1,2) = Plane % x(1) - Plane % x(2)
    A(2,2) = Plane % y(1) - Plane % y(2)

    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))
    
    IF(ABS(detA) <= eps * absA + 1.0d-20) RETURN
    
    B(1) = Plane % x(1) - Line % x(1) 
    B(2) = Plane % y(1) - Line % y(1) 
    
    CALL InvertMatrix( A,2 )
    C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
    
    IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) RETURN
    
    Inside = .TRUE.
    frac = C(1)
    X0 = Line % x(1) + C(1) * (Line % x(2) - Line % x(1))
    Y0 = Line % y(1) + C(1) * (Line % y(2) - Line % y(1))
    Z0 = Line % z(1) + C(1) * (Line % z(2) - Line % z(1))
    
  END SUBROUTINE LineIntersectionCoords
  
!------------------------------------------------------------------------------
END SUBROUTINE SaveScalars
!------------------------------------------------------------------------------

!> \}
