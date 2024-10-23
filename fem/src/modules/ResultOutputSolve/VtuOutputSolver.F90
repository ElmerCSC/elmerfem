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
MODULE VtuXMLFile
  USE DefUtils 
  USE MeshUtils
  USE SolverUtils
  USE SaveUtils
  USE MainUtils
  USE ElementDescription
  USE AscBinOutputUtils

  IMPLICIT NONE 
  
CONTAINS


  ! Check whether there is any discontinuous galerkin field to be saved. 
  ! It does not make sense to use discontinuous saving if there are no discontinuous fields.
  ! It will even result to errors since probably there are no DG indexes either. 
  !-----------------------------------------------------------------------------------------
  FUNCTION CheckAnyDGField(Model,Params) RESULT ( HaveAnyDG ) 
    TYPE(Model_t) :: Model
    TYPE(ValueList_t), POINTER :: Params

    
    LOGICAL :: HaveAnyDG
    INTEGER :: Rank, Vari, VarType
    CHARACTER(LEN=1024) :: Txt, FieldName
    TYPE(Variable_t), POINTER :: Solution
    LOGICAL :: Found
    
    HaveAnyDG = .FALSE.

    DO Rank = 0,1
      DO Vari = 1, 999
        IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        
        FieldName = GetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName),ThisOnly=.TRUE.)
        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly=.TRUE.)
        END IF
        IF( .NOT. ASSOCIATED( Solution ) ) CYCLE

        VarType = Solution % Type
        
        IF ( VarType == Variable_on_nodes_on_elements .OR. &
            VarType == Variable_on_gauss_points ) THEN
          HaveAnyDG = .TRUE.
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION CheckAnyDGField



  ! Average fields within bodies. Offers good compromise between file size
  ! and honoring discontinuities. 
  !-----------------------------------------------------------------------
  SUBROUTINE AverageBodyFields( Mesh ) 
    
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(Variable_t), POINTER :: Var, Var1
    INTEGER :: NoAve, i
    LOGICAL :: BodySum    
    
    ! The variables may have been averaged already in vector form. 
    ! Inherit the DgAveraged flag to the components. 
    Var => Mesh % Variables    
    DO WHILE( ASSOCIATED( Var ) )       
      IF ( Var % DOfs > 1 .AND. Var % TYPE == Variable_on_nodes_on_elements ) THEN        
        DO i=1,Var % Dofs
          Var1 => VariableGet( Mesh % Variables,ComponentName(Var % Name,i), ThisOnly = .TRUE.)
          IF( ASSOCIATED(Var1)) Var1 % DgAveraged = Var % DgAveraged 
        END DO
      END IF
      Var => Var % Next
    END DO

    NoAve = 0
    Var => Mesh % Variables    
    DO WHILE( ASSOCIATED( Var ) ) 
      
      ! Skip if variable is not active for saving       
      IF ( .NOT. Var % Output ) THEN
        CONTINUE
      ! Skip if variable is global one
      ELSE IF ( SIZE( Var % Values ) == Var % DOFs ) THEN  
        CONTINUE
      ! Each field is present componentwise and as a vector. 
      ! Only do the components, and this takes care of the vectors as well.
      ELSE IF ( Var % DOFs > 1 ) THEN  
        CONTINUE
      ! And finally do the everaging for remaining DG fields only
      ELSE IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN
        NoAve = NoAve + 1

        ! This is really quite dirty!
        ! For variables that scale with h the operator is more naturally a sum 
        ! than an average. 
        !BodySum = ( INDEX( Var % Name,'nodal force' ) /= 0 .OR. &
            !INDEX( Var % Name,' loads' ) /= 0 )
        CALL CalculateBodyAverage( Mesh, Var, .FALSE. )
      END IF

      Var => Var % Next
    END DO

    CALL Info('AverageBodyFields','Reduced '//I2S(NoAve)//' elemental fields',Level=7)

  END SUBROUTINE AverageBodyFields

  
  
  ! Write the filename for saving .vtu, .pvd, and .pvtu files.
  ! Partname, partition and timestep may be added to the name.
  !------------------------------------------------------------------------------------
  SUBROUTINE VtuFileNaming( BaseFile, VtuFile, Suffix, GroupId, FileIndex, Part, NoPath, ParallelBase ) 
    CHARACTER(LEN=*), INTENT(IN) :: BaseFile, Suffix
    CHARACTER(LEN=*), INTENT(INOUT) :: VtuFile
    INTEGER :: GroupId, FileIndex
    INTEGER, OPTIONAL :: Part
    LOGICAL, OPTIONAL :: NoPath
    LOGICAL, OPTIONAL :: ParallelBase

    CHARACTER(MAX_NAME_LEN) :: GroupName
    INTEGER :: i,j,NameOrder(3),PEs
    LOGICAL :: LegacyMode, ParallelBaseName, NoFileindex, Found 

    
    NameOrder = [1,2,3]
    LegacyMode = .FALSE.
    ParallelBaseName = .FALSE.
    IF( PRESENT( ParallelBase ) ) ParallelBaseName = ParallelBase

    NoFileindex = ListGetLogical( CurrentModel % Solver % Values,'No Fileindex',Found )  

    
    VtuFile = BaseFile

    DO j = 1, 3
    
      ! Append vtu file name with group name      
      SELECT CASE( NameOrder(j) )

      CASE( 1 ) 
        ! If we have groups then the piece is set to include the name of the body/bc. 
        IF( GroupId > 0 ) THEN
          IF( GroupId <= CurrentModel % NumberOfBodies ) THEN
            GroupName = ListGetString( CurrentModel % Bodies(GroupId) % Values,"Name", Found)
            IF(.NOT. Found) GroupName = 'Body'//I2S(GroupId)
          ELSE
            i = GroupId - CurrentModel % NumberOfBodies
            GroupName = ListGetString( CurrentModel % BCs(i) % Values,"Name",Found)
            IF(.NOT. Found) GroupName = 'BC'//I2S(i)
          END IF
          VtuFile = TRIM(VtuFile)//"_"//TRIM(GroupName)
        END IF

      CASE( 2 )         
        PEs = ParEnv % PEs
        IF( PEs == 1 ) CYCLE

        IF( PRESENT( Part ) ) THEN
          ! In parallel the mesh consists of pieces called partitions.
          ! Give each partition a name that includes the partition. 
          IF( LegacyMode ) THEN
            WRITE( VtuFile,'(A,A,I4.4,A)') TRIM((VtuFile)),"_",Part,"par"            
          ELSE
            IF ( PEs < 10) THEN                    
              WRITE( VtuFile,'(A,A,I1.1,A,I1.1)') TRIM((VtuFile)),"_",PEs,"np",Part
            ELSE IF ( PEs < 100) THEN                    
              WRITE( VtuFile,'(A,A,I2.2,A,I2.2)') TRIM((VtuFile)),"_",PEs,"np",Part
            ELSE IF ( PEs < 1000) THEN                    
              WRITE( VtuFile,'(A,A,I3.3,A,I3.3)') TRIM((VtuFile)),"_",PEs,"np",Part
            ELSE
              WRITE( VtuFile,'(A,A,I4.4,A,I4.4)') TRIM((VtuFile)),"_",PEs,"np",Part
            END IF
          END IF
        ELSE
          ! Also add number to the wrapper files as it is difficult otherwise
          ! quickly see on which partitioning they were computed. 
          IF( ParallelBaseName ) THEN
            IF ( PEs < 10) THEN                    
              WRITE( VtuFile,'(A,A,I1.1,A)') TRIM((VtuFile)),"_",PEs,"np"
            ELSE IF ( PEs < 100) THEN                    
              WRITE( VtuFile,'(A,A,I2.2,A)') TRIM((VtuFile)),"_",PEs,"np"
            ELSE IF ( PEs < 1000) THEN                    
              WRITE( VtuFile,'(A,A,I3.3,A)') TRIM((VtuFile)),"_",PEs,"np"
            ELSE
              WRITE( VtuFile,'(A,A,I4.4,A)') TRIM((VtuFile)),"_",PEs,"np"
            END IF
          END IF
        END IF
          
      CASE( 3 )         
        ! This is for adding time (or nonlinear iteration/scanning) to the filename.
        IF( FileIndex > 0 .AND. .NOT. NoFileindex ) THEN
          IF( FileIndex < 10000 ) THEN        
            WRITE(VtuFile,'(A,A,I4.4)') TRIM((VtuFile)),"_t",FileIndex
          ELSE
            WRITE(VtuFile,'(A,A,I0)' ) TRIM((VtuFile)),"_t",FileIndex
          END IF
        END IF     
        
      END SELECT
    END DO

    IF( PRESENT( NoPath ) ) THEN
      IF( NoPath ) THEN    
        j = INDEX( VtuFile,'/',BACK=.TRUE.) 
        IF( j > 0 ) VtuFile = VtuFile(j+1:)
      END IF
    END IF
          
    VtuFile = TRIM( VtuFile)//TRIM(Suffix) 

  END SUBROUTINE VtuFileNaming
  
END MODULE VtuXMLFile


!------------------------------------------------------------------------------
!> Subroutine for saving the results in XML based VTK format (VTU). Both ascii and binary
!> output is available, in single or double precision. The format is understood by 
!> visualization software Paraview and ViSit, for example.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VtuOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE VtuXMLFile
    
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(dp) :: dt
  LOGICAL :: TransientSimulation
  
  INTEGER, SAVE :: nTime = 0
  LOGICAL :: GotIt, Parallel, FixedMesh, DG, DN, DoAve
  CHARACTER(MAX_NAME_LEN) :: FilePrefix
  CHARACTER(MAX_NAME_LEN) :: BaseFile, VtuFile, PvtuFile, PvdFile, DataSetFile
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: i, j, k, l, n, m, Partitions, Part, ExtCount, FileindexOffSet, MeshDim, PrecBits, &
             PrecSize, IntSize, FileIndex
  CHARACTER(:), ALLOCATABLE :: OutputDirectory
  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp) :: DoubleWrk
  REAL :: SingleWrk
  LOGICAL :: BinaryOutput, AsciiOutput, SinglePrec, NoFileindex
  CHARACTER(MAX_NAME_LEN) :: Str
  INTEGER, POINTER :: InvFieldPerm(:)
  INTEGER, ALLOCATABLE, TARGET :: NodePerm(:), InvNodePerm(:), InvDgPerm(:), DgPerm(:)
  INTEGER :: NumberOfGeomNodes, NumberOfDofNodes, NumberOfElements, ParallelNodes, ParallelElements
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(ValueList_t),POINTER :: Params
  INTEGER :: MaxModes, MaxModes2, BCOffset, ElemFirst, ElemLast, &
      OutputMeshes, ParallelDofsNodes, LagN
  INTEGER, POINTER :: ActiveModes(:), ActiveModes2(:)
  LOGICAL :: GotActiveModes, GotActiveModes2, EigenAnalysis, &
      WriteIds, SaveLinear, SaveMetainfo, &
      NoPermutation, SaveElemental, SaveNodal, NoInterp
  LOGICAL, ALLOCATABLE :: ActiveElem(:)
  INTEGER, ALLOCATABLE :: GeometryBodyMap(:),GeometryBCMap(:)

! Parameters for buffered binary output
  INTEGER :: BufferSize

  LOGICAL :: TimeCollection, GroupCollection, ParallelBase
  INTEGER :: GroupId, EigenVectorMode
  CHARACTER(*), PARAMETER :: Caller = 'VtuOutputSolver'

  Params => GetSolverParams()
  Mesh => Model % Mesh
  MeshDim = Mesh % MeshDim

  NoInterp = .TRUE.
  IF( ListGetLogical( Params,'Enable Interpolation',GotIt ) ) THEN
    NoInterp = .FALSE.
  END IF
  
  DG = GetLogical( Params,'Discontinuous Galerkin',GotIt)
  DN = GetLogical( Params,'Discontinuous Bodies',GotIt)
  IF( DG .OR. DN ) THEN    
    IF(.NOT. CheckAnyDGField(Model,Params) ) THEN
      CALL Info(Caller,'No DG or IP fields, omitting discontinuity creation!',Level=6)
      DG = .FALSE. 
      DN = .FALSE.
    END IF

    IF( DG .OR. DN ) THEN
      ! Sometimes we have a request to save in DG format even though no equation has been solved as dg.
      ! Then we need to create the Element % DgIndexes for the saving only. If already done this does nothing.
      CALL CheckAndCreateDGIndexes( Mesh )

      DoAve = GetLogical( Params,'Average Within Materials', GotIt)
      IF(.NOT. GotIt) DoAve = .TRUE. 
      
      IF( DoAve ) CALL AverageBodyFields( Mesh )  
    END IF
  END IF
  
  LagN = GetInteger( Params,'Lagrange Element Degree',GotIt) 

  ExtCount = GetInteger( Params,'Output Count',GotIt)
  IF( GotIt ) THEN
    nTime = ExtCount
  ELSE
    i = GetInteger( Params,'Fileindex step',GotIt)
    IF( GotIt ) THEN
      nTime = nTime + i
    ELSE
      nTime = nTime + 1
    END IF
  END IF

  FileIndexOffset = GetInteger( Params,'Fileindex offset',GotIt)
  FileIndex = nTime + FileIndexOffset

  i = GetInteger( Params, 'Output File cycle',GotIt)
  IF(GotIt) FileIndex = MODULO(FileIndex-1,i)+1
  
  BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
  IF( GotIt ) THEN
    AsciiOutput = .NOT. BinaryOutput
  ELSE
    AsciiOutput = GetLogical( Params,'Ascii Output',GotIt)
    BinaryOutput = .NOT. AsciiOutput
  END IF

  ! Do not save metainfo if we use the results of saving for consistency test.
  SaveMetainfo = .NOT. ( GetLogical( Params,'Skip Metainfo',GotIt ) .OR. &
      ListCheckPresent( Params,'Reference Values') )
  
  ParallelBase = GetLogical( Params,'Partition Numbering',GotIt )
  
  IF( BinaryOutput ) THEN
    BufferSize = GetInteger( Params,'Binary Output Buffer Size',GotIt)
    IF( .NOT. GotIt ) BufferSize = MAX(1000, Mesh % NumberOfNodes )
  END IF
  
  SaveElemental = GetLogical( Params,'Save Elemental Fields',GotIt)
  IF(.NOT. GotIt) SaveElemental = .TRUE.
    
  SaveNodal = GetLogical( Params,'Save Nodal Fields',GotIt) 
  IF(.NOT. GotIt) SaveNodal = .TRUE.

  SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
  IF( SinglePrec ) THEN
    CALL Info(Caller,'Using single precision arithmetics in output!',Level=7)
  END IF

  IF( SinglePrec ) THEN
    PrecBits = 32
    PrecSize = KIND( SingleWrk ) 
  ELSE
    PrecBits = 64
    PrecSize = KIND( DoubleWrk ) 
  END IF
  IntSize = KIND(i)

  OutputMeshes = ListGetInteger(Params,'Number of Output Meshes',GotIt)

  Partitions = ParEnv % PEs
  Part = ParEnv % MyPE
  Parallel = (Partitions > 1) .OR. GetLogical(Params,'Enforce Parallel format',GotIt)

  NoFileindex = GetLogical( Params,'No Fileindex',GotIt)

  SaveLinear = GetLogical( Params,'Save Linear Elements',GotIt)

  FilePrefix = GetString( Params,'Output File Name',GotIt )
  IF ( .NOT.GotIt ) FilePrefix = "Output"
  IF ( Mesh % DiscontMesh ) THEN
    FilePrefix = 'discont_'//TRIM(FilePrefix)    
  ELSE IF( OutputMeshes > 1 ) THEN
    i = INDEX( Mesh % Name,'/',.TRUE.)
    IF( i > 0 ) THEN      
      FilePrefix = TRIM(Mesh % Name(i+1:))//'_'//TRIM(FilePrefix)
    ELSE
      FilePrefix = TRIM(Mesh % Name)//'_'//TRIM(FilePrefix)      
    END IF
  END IF
    
  
  IF ( nTime == 1 ) THEN
    CALL Info(Caller,'Saving results in VTK XML format with prefix: '//TRIM(FilePrefix))
    CALL Info(Caller, 'Saving number of partitions: '//I2S(Partitions))
  END IF

  BaseFile = FilePrefix

  CALL SolverOutputDirectory( Solver, BaseFile, OutputDirectory, UseMeshDir = .TRUE.  )
  BaseFile = TRIM(OutputDirectory)// '/' //TRIM(BaseFile)
  
  CALL Info(Caller,'Full filename base is: '//TRIM(Basefile), Level=10 )
    
  FixedMesh = ListGetLogical(Params,'Fixed Mesh',GotIt)

  TimeCollection = GetLogical( Params,'Vtu Time Collection', GotIt ) 
  IF( TimeCollection ) THEN
    IF( NoFileIndex ) THEN
      CALL Warn(Caller,'Vtu time collection cannot work without file indexes')
      NoFileIndex = .FALSE.
    END IF
    IF( .NOT. TransientSimulation ) THEN
      CALL Warn(Caller,'Vtu time collection requires a transient simulation!')
      TimeCollection = .FALSE.
    END IF
  END IF
  
  GroupCollection = GetLogical( Params,'Vtu Part Collection', GotIt ) 

  GroupId = 0    
200 CONTINUE
  IF( GroupCollection ) THEN
    GroupId = GroupId + 1
    CALL Info(Caller,'Saving group '//I2S(GroupId),Level=8)
  END IF

  !------------------------------------------------------------------------------
  ! Initialize stuff for masked saving
  !------------------------------------------------------------------------------
  CALL GenerateSaveMask(Mesh,Params,Parallel,GroupId,SaveLinear,&
      NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements,&
      ElemFirst,ElemLast)
  
  !------------------------------------------------------------------------------
  ! If we have a discontinuous mesh then create the permutation vectors to deal
  ! with the discontinuities.
  !------------------------------------------------------------------------------
  CALL GenerateSavePermutation(Mesh,DG,DN,LagN,SaveLinear,ActiveElem,NumberOfGeomNodes,&
      NoPermutation,NumberOfDofNodes,DgPerm,InvDgPerm,NodePerm,InvNodePerm)
  
  ! The partition is active for saving if there are any nodes 
  ! to write. There can be no elements nor dofs without nodes.
  CALL ParallelActive( NumberOfDofNodes > 0 )

  IF( nTime == 1 ) THEN
    ParallelElements = ParallelReduction( NumberOfElements ) 

    IF( ParallelElements == 0 ) THEN
      CALL Info(Caller, 'Nothing to save for this selection: ',Level=8)
    ELSE
      CALL Info(Caller, 'Total number of elements to save: '&
          //I2S(ParallelElements),Level=8)
      
      ParallelNodes = ParallelReduction( NumberOfGeomNodes ) 
      CALL Info(Caller, 'Total number of geometry nodes to save: '&
          //I2S(ParallelNodes),Level=8)
      
      ParallelNodes = ParallelReduction( NumberOfDofNodes ) 
      CALL Info(Caller, 'Total number of dof nodes to save: '&
          //I2S(ParallelNodes),Level=8)
    END IF
  END IF

  ! Sometimes we may want to ignore eigenmodes completely.
  !--------------------------------------------------------------
  IF( ListGetLogical( Params,'Ignore Eigenmodes',GotIt) ) THEN
    EigenAnalysis = .FALSE.
    EigenVectorMode = 0
    MaxModes = 0
    MaxModes2 = 0
    GOTO 1 
  END IF
  
  !------------------------------------------------------------------------------
  ! Check whether we have nodes coming from different reasons
  !------------------------------------------------------------------------------       
  ActiveModes => ListGetIntegerArray( Params,'Active EigenModes',GotActiveModes ) 
  IF( GotActiveModes ) THEN
    MaxModes = SIZE( ActiveModes )
  ELSE
    MaxModes = GetInteger( Params,'Number of EigenModes',GotIt)
    IF(.NOT. GotIt) MaxModes = GetInteger( Params,'Eigen System Values',GotIt)
    IF(.NOT. GotIt) THEN
      DO i=1,Model % NumberOfSolvers
        MaxModes = MAX( MaxModes, &
            GetInteger( Model % Solvers(i) % Values,'Eigen System Values', GotIt ) )
        MaxModes = MAX( MaxModes, &
            GetInteger( Model % Solvers(i) % Values,'Harmonic System Values', GotIt ) )       
        IF( ListGetLogical( Model % Solvers(i) % Values,'Save Scanning Modes',GotIt ) ) THEN
          MaxModes = MAX( MaxModes, &
              GetInteger( Model % Solvers(i) % Values,'Scanning Loops', GotIt ) )
        END IF
      END DO
    END IF
  END IF
    
  IF( MaxModes > 0 ) THEN
    CALL Info(Caller,'Maximum number of eigen/harmonic modes: '//I2S(MaxModes),Level=7)
    Str = ListGetString( Params,'Eigen Vector Component', GotIt )
    IF( GotIt ) THEN
      IF( Str == 're') THEN
        CONTINUE
      ELSE IF( Str == 'im' ) THEN
        EigenVectorMode = 1
      ELSE IF( Str == 'abs' ) THEN
        EigenVectorMode = 2
      ELSE IF( Str == 'complex' ) THEN
        EigenVectorMode = 3
      ELSE
        CALL Fatal(Caller,'Invalid value for >Eigen System Mode< :'//TRIM(str))
      END IF
      CALL Info(Caller,'Using eigen vector mode: '//I2S(EigenVectorMode),Level=7)
    END IF
  END IF

  ActiveModes2 => ListGetIntegerArray( Params,'Active Constraint Modes',GotActiveModes2 ) 
  IF( GotActiveModes2 ) THEN
    MaxModes2 = SIZE( ActiveModes2 )
  ELSE
    MaxModes2 = 0
    DO i=1,Model % NumberOfSolvers
      IF( .NOT. ASSOCIATED( Model % Solvers(i) % Variable ) ) CYCLE
      MaxModes2 = MAX( MaxModes2, &
          Model % Solvers(i) % Variable % NumberOfConstraintModes )
    END DO
  END IF
  IF( MaxModes2 > 0 ) THEN
    CALL Info(Caller,'Maximum number of constraint modes: '//I2S(MaxModes2),Level=7)
  END IF

  ! This activates the solution of the modes one for each file
  EigenAnalysis = ListGetLogical( Params,'Eigen Analysis',GotIt) .OR. &
      ListGetLogical( Params,'Constraint Modes Analysis',GotIt) 
  IF( EigenAnalysis ) THEN
    CALL Info(Caller,'Saving each mode to different file')
    FileIndex = 1
  END IF

  ! Let's jump here if we ignore all modes. 
1 CONTINUE
  
  BcOffset = 0
  WriteIds = GetLogical( Params,'Save Geometry Ids',GotIt)  
  IF( WriteIds ) THEN
    ! Create the mapping for body ids, default is unity mapping
    IF(.NOT. ALLOCATED(GeometryBodyMap)) THEN
      ALLOCATE( GeometryBodyMap( CurrentModel % NumberOfBodies ) )
    END IF
    j = ListGetInteger( Params,'Default Body Id',GotIt )
    IF( GotIt ) THEN
      GeometryBodyMap = j
    ELSE
      DO i=1,CurrentModel % NumberOfBodies
        GeometryBodyMap(i) = i
      END DO
    END IF

    ! User given mapping
    DO i=1,CurrentModel % NumberOfBodies
      j = ListGetInteger( CurrentModel % Bodies(i) % Values,'Geometry Id',GotIt)
      IF( GotIt ) GeometryBodyMap(i) = j
    END DO

    ! Create mapping for bc ids, default is unity mapping with offset
    IF( .NOT. ALLOCATED( GeometryBCMap ) ) THEN
      ALLOCATE( GeometryBCMap( CurrentModel % NumberOfBCs ) )
    END IF
    j = ListGetInteger( Params,'Default BC Id',GotIt )
    IF( GotIt ) THEN
      GeometryBCMap = j
    ELSE
      ! Determine a default offset
      BCOffset = ListGetInteger( Params,'BC Id Offset',GotIt )
      IF( .NOT. GotIt ) THEN
        IF( ElemFirst <= Mesh % NumberOfBulkElements ) THEN
          BCOffset = 100
          DO WHILE( BCOffset <= Model % NumberOfBodies ) 
            BCOffset = 10 * BCOffset
          END DO
          CALL Info(Caller,'Setting offset for boundary entities: '&
              //I2S(BCOffset),Level=6)
        END IF
      END IF
      DO i=1,CurrentModel % NumberOfBCs
        GeometryBCMap(i) = i + BCOffSet
      END DO
    END IF

    ! User given bc mapping
    DO i=1,CurrentModel % NumberOfBCs
      j = ListGetInteger( CurrentModel % BCs(i) % Values,'Geometry Id',GotIt)
      IF( GotIt ) GeometryBCMap(i) = j
    END DO

  END IF
  
  ! We may need to jump here to write a new eigenmode
100 CONTINUE

  ParallelDofsNodes = ParallelReduction( NumberOfDofNodes ) 
  
  IF(Parallel) THEN
    ! Generate the filename for saving
    !--------------------------------------------------------------------
    CALL VtuFileNaming( BaseFile, PvtuFile,'.pvtu', GroupId, FileIndex, ParallelBase = ParallelBase )
    CALL Info(Caller,'Writing the pvtu file: '//TRIM(PvtuFile), Level=10)
    CALL WritePvtuFile( PVtuFile, Model )
    CALL Info(Caller,'Finished writing pvtu file',Level=12)
  END IF


  ! Write the Vtu file with all the data
  !--------------------------------------------------------------------------
  IF( NumberOfDofNodes > 0 ) THEN
    CALL VtuFileNaming( BaseFile, VtuFile,'.vtu', GroupId, FileIndex, Part+1 ) 
    CALL Info(Caller,'Writing the vtu file: '//TRIM(VtuFile),Level=7)
    CALL WriteVtuFile( VtuFile, Model, FixedMesh )
    CALL Info(Caller,'Finished writing vtu file',Level=12)
  END IF

  ! For transient simulation or group collections write a holder for individual files
  !-----------------------------------------------------------------------------------
  IF( TimeCollection .OR. GroupCollection ) THEN
    CALL VtuFileNaming( BaseFile, PvdFile,'.pvd', GroupId, FileIndex, ParallelBase = ParallelBase )    
    WRITE( PvdFile,'(A,".pvd")' ) TRIM(BaseFile)
    IF( Parallel ) THEN
      CALL VtuFileNaming( BaseFile, DataSetFile,'.pvtu', GroupId, FileIndex, &
          NoPath = .TRUE., ParallelBase = ParallelBase ) 
    ELSE      
      CALL VtuFileNaming( BaseFile, DataSetFile,'.vtu', GroupId, FileIndex, NoPath = .TRUE. ) 
    END IF

    IF( ParallelDofsNodes == 0 ) THEN
      CALL Info(Caller,'Nothing to write in pvd file: '//TRIM(DataSetFile),Level=10)
    ELSE
      CALL Info(Caller,'Writing the pvd file: '//TRIM(DataSetFile),Level=10)
      CALL WritePvdFile( PvdFile, DataSetFile, FileIndex, Model )
      CALL Info(Caller,'Finished writing pvd file',Level=12)     
    END IF
  END IF


  IF( EigenAnalysis ) THEN
    FileIndex = FileIndex + 1
    IF( FileIndex <= MaxModes + MaxModes2 ) GOTO 100
  END IF

  IF( GroupCollection ) THEN
    IF( GroupId < CurrentModel % NumberOfBodies + CurrentModel % NumberOfBCs ) THEN
      GOTO 200 
    END IF
  END IF

  
  IF( ALLOCATED( NodePerm ) ) DEALLOCATE( NodePerm ) 
  IF( ALLOCATED( ActiveElem ) ) DEALLOCATE( ActiveElem ) 
  IF( ALLOCATED( InvNodePerm ) ) DEALLOCATE( InvNodePerm )
  IF( ALLOCATED( GeometryBodyMap ) ) DEALLOCATE( GeometryBodyMap )
  IF( ALLOCATED( GeometryBCMap ) ) DEALLOCATE( GeometryBcMap ) 
  
  CALL Info(Caller,'All done for now',Level=10)     


CONTAINS


  FUNCTION PickComplex(zval,zmode) RESULT( val ) 
    COMPLEX(KIND=dp) :: zval
    INTEGER :: zmode
    REAL(KIND=dp) :: val
    
    SELECT CASE( zmode ) 
    CASE(0,3) 
      val = REAL( zval )
    CASE(1,4) 
      val = AIMAG( zval ) 
    CASE(2) 
      val = ABS( zval ) 
    END SELECT
  END FUNCTION PickComplex
    
  
  
  ! Writes a single VTU file that can be read by Paraview, ViSiT etc.
  !---------------------------------------------------------------------------------------
  SUBROUTINE WriteVtuFile( VtuFile, Model, RemoveDisp )
    CHARACTER(LEN=*), INTENT(IN) :: VtuFile
    TYPE(Model_t) :: Model 
    LOGICAL, INTENT(IN) :: RemoveDisp
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: i,ii,j,jj,k,dofs,Rank,n,m,dim,vari,sdofs,dispdofs, dispBdofs, Offset, &
        NoFields, NoFields2, IndField, iField, NoModes, NoModes2, NoFieldsWritten, &
        cumn, iostat
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FieldNameB, OutStr
    CHARACTER :: lf
    CHARACTER(*), PARAMETER :: Caller = 'WriteVtuFile'
    LOGICAL :: ScalarsExist, VectorsExist, Found, &
        ComponentVector, ComponentVectorB, ComplementExists, &
        Use2, IsHarmonic, FlipActive
    INTEGER :: DoIm
    LOGICAL :: WriteData, WriteXML, L, Buffered
    TYPE(Variable_t), POINTER :: Solution, Solution2, Solution3, TmpSolDg, TmpSolDg2, TmpSolDg3
    INTEGER, POINTER :: Perm(:), PermB(:), DispPerm(:), DispBPerm(:)
    REAL(KIND=dp), POINTER :: Values(:), Values2(:), Values3(:), DispValues(:)
    REAL(KIND=dp), POINTER :: ValuesB(:), ValuesB2(:), ValuesB3(:), DispBValues(:)
    REAL(KIND=dp) :: x,y,z, val,ElemVectVal(3)
    INTEGER, ALLOCATABLE, TARGET :: ElemInd(:)
    INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729
    INTEGER :: TmpIndexes(MAX_LAGRANGE_NODES), VarType
    
    COMPLEX(KIND=dp), POINTER :: EigenVectors(:,:), EigenVectors2(:,:), EigenVectors3(:,:)
    COMPLEX(KIND=dp), POINTER :: EigenVectorsB(:,:), EigenVectorsB2(:,:), EigenVectorsB3(:,:)
    COMPLEX(KIND=dp) :: zval
    REAL(KIND=dp), POINTER :: ConstraintModes(:,:)
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: CurrentElement, Parent
    TYPE(ValueList_t), POINTER :: Params
    REAL(KIND=dp), POINTER :: TmpArray(:,:)
    REAL(KIND=dp) :: CoordScale(3), CoordOffset(3)
    TYPE(Variable_t), POINTER, SAVE :: LagVarX=>NULL(), LagVarY=>NULL(), LagVarZ=>NULL()
    
    ! Initialize the auxiliary module for buffered writing
    !--------------------------------------------------------------
    CALL AscBinWriteInit( AsciiOutput, SinglePrec, VtuUnit, BufferSize )

    ! Linefeed character
    !-----------------------------------
    lf = CHAR(10)

    Offset = 0
    WriteXML = .TRUE.
    WriteData = AsciiOutput
    Params => GetSolverParams()
    Buffered = .TRUE.
    FlipActive = .FALSE.
    TmpSolDg => NULL()
    TmpSolDg2 => NULL()
    TmpSolDg3 => NULL()
    
    ! we could have huge amount of gauss points
    ALLOCATE( ElemInd(512)) !Model % Mesh % MaxElementDOFS))

    ! This format works both for ascii and binary output
    !-------------------------------------------------------------------------
    OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream', STATUS='replace', IOStat=iostat)
    IF( iostat /= 0 ) THEN
      CALL Fatal(Caller,'Opening of file failed: '//TRIM(VtuFile))
    END IF
        
    Solver => Model % Solver

    ! VTU seemingly only works with 3D cases, so enforce it
    dim = 3


    WRITE( OutStr,'(A)') '<?xml version="1.0"?>'//lf
    CALL AscBinStrWrite( OutStr ) 

    IF ( LittleEndian() ) THEN
      OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
    ELSE
      OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'//lf
    END IF
    CALL AscBinStrWrite( OutStr )

    IF( SaveMetainfo ) THEN
      IF( FileIndex == 1 .AND. ParEnv % MyPe == 0 ) THEN
        Txt = GetVersion()
        WRITE( OutStr,'(A)') '<!-- Elmer version: '//TRIM(Txt)//' -->'//lf     
        CALL AscBinStrWrite( OutStr )

        Txt = GetRevision( Found )
        IF( Found ) THEN
          WRITE( OutStr,'(A)') '<!-- Elmer revision: '//TRIM(Txt)//' -->'//lf
          CALL AscBinStrWrite( OutStr )
        END IF

        Txt = GetCompilationDate( Found )
        IF( Found ) THEN
          WRITE( OutStr,'(A)') '<!-- Elmer compilation date: '//TRIM(Txt)//' -->'//lf
          CALL AscBinStrWrite( OutStr )
        END IF

        Txt = GetSifName( Found) 
        IF( Found ) THEN
          WRITE( OutStr,'(A)') '<!-- Solver input file: '//TRIM(Txt)//' -->'//lf
          CALL AscBinStrWrite( OutStr )
        END IF

        Txt = FormatDate()      
        WRITE( OutStr,'(A)') '<!-- File started at: '//TRIM(Txt)//' -->'//lf
        CALL AscBinStrWrite( OutStr )
      END IF
    END IF
      
    WRITE( OutStr,'(A)') '  <UnstructuredGrid>'//lf
    CALL AscBinStrWrite( OutStr )
    WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NumberOfDofNodes,&
        '" NumberOfCells="',NumberOfElements,'">'//lf
    CALL AscBinStrWrite( OutStr )

    ! nodewise information
    !-------------------------------------
    ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
    VectorFieldName = GetString( Params,'Vector Field 1',VectorsExist)
    IF( .NOT. ( ScalarsExist .OR. VectorsExist) ) THEN
      CALL Warn(Caller,'Are there really no scalars or vectors?')
    END IF

    WRITE( OutStr,'(A)') '      <PointData>'//lf
    CALL AscBinStrWrite( OutStr )

    DispDofs = 0
    DispBDofs = 0
    IF(RemoveDisp) THEN
      Solution => VariableGet( Model % Mesh % Variables, 'Displacement',ThisOnly=NoInterp)
      IF( ASSOCIATED( Solution ) ) THEN
        Solver => Solution % Solver
        L = GetLogical( GetSolverParams(Solver),'Displace Mesh',Found)
        IF(.NOT.Found) L=.NOT.EigenOrHarmonicAnalysis(Solver)
        IF (L) THEN
          DispPerm => Solution % Perm
          DispValues => Solution % Values
          DispDofs = Solution % Dofs
        END IF
      END IF

      Solution => VariableGet( Model % Mesh % Variables, 'Mesh Update',ThisOnly=NoInterp)
      IF( ASSOCIATED( Solution ) ) THEN
        DispBPerm => Solution % Perm
        DispBValues => Solution % Values
        DispBDofs = Solution % Dofs
        CALL Info(Caller,'Automatically complement > Displacement < by > Mesh Update < field',Level=7)
      END IF
    END IF

    CoordScale = 1.0_dp
    IF( ListGetLogical( Params,'Coordinate Scaling Revert', Found ) ) THEN
      TmpArray => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',Found )    
      IF( Found ) THEN            
        DO i=1,Model % Mesh % MaxDim 
          j = MIN( i, SIZE(TmpArray,1) )
          CoordScale(i) = 1.0_dp / TmpArray(j,1)
        END DO
      END IF
    ELSE
      ! Enable local coordinate scaling for convenience. 
      TmpArray => ListGetConstRealArray( Params,'Coordinate Scaling',Found )    
      IF( Found ) THEN            
        DO i=1,Model % Mesh % MaxDim 
          j = MIN( i, SIZE(TmpArray,1) )
          CoordScale(i) = TmpArray(j,1)
        END DO
      END IF
    END IF
      
    CoordOffset = 0.0_dp
    TmpArray => ListGetConstRealArray( Params,'Mesh Translate',Found )    
    IF( Found ) THEN            
      DO i=1,MIN( 3, SIZE(TmpArray,1) )
        CoordOffset(i) = TmpArray(i,1)
      END DO
    ELSE
      DO i=1,3
        CoordOffset(i) = ListGetCReal( Params,'Mesh Translate '//I2S(i),Found )
      END DO
    END IF
    

    ! When the data is 'appended' two loops will be taken and the data will be written
    ! on the second loop. Offset is the position in the appended data after the '_' mark.
    !------------------------------------------------------------------------------------
100 Offset = 0

    IF( SaveNodal ) THEN
      CALL Info(Caller,'Writing nodal fields',Level=10)
      NoFieldsWritten = 0
      DO Rank = 0,2
        DO Vari = 1, 999
          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

          FieldName = GetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT

          IF(Rank == 2) THEN
            CALL Fatal(Caller,'Do the tensors')
          END IF

          !---------------------------------------------------------------------
          ! Find the variable with the given name in the normal manner 
          !---------------------------------------------------------------------
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName),ThisOnly=NoInterp)

          ComponentVector = .FALSE.
          IF(ASSOCIATED(Solution)) THEN
            dofs = Solution % DOFs
            NULLIFY(Solution2)
            NULLIFY(Solution3)
          ELSE 
            !---------------------------------------------------------------------
            ! Some vectors are defined by a set of components (either 2 or 3 is possible!)
            !---------------------------------------------------------------------            
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution)) THEN 
              dofs = 1
              Solution2 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly=NoInterp)
              IF(ASSOCIATED(Solution2)) dofs = 2
              Solution3 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly=NoInterp)
              IF(ASSOCIATED(Solution3)) dofs = 3
              ComponentVector = .TRUE.              
            ELSE
              CALL Warn(Caller,'Nonexistent variable: '//TRIM(FieldName)) 
              CYCLE
            END IF
          END IF

          CALL Info(Caller,'Saving variable: '//TRIM(FieldName),Level=10)
          
          VarType = Solution % Type

          IF( LagN > 0 ) THEN
            ! Let us create a L-field on-the-fly and have the solution vector point to that.
            IF(.NOT. ASSOCIATED(TmpSolDg)) THEN
              ALLOCATE(TmpSolDg,TmpSolDg2,TmpSolDg3)
            END IF
            CALL p2LagrangeSwapper( Mesh, Solution, TmpSolDg, LagN, DgPerm, NumberOfDofNodes )
            Solution => TmpSolDg
            IF( ASSOCIATED(Solution2) ) THEN              
              CALL p2LagrangeSwapper( Mesh, Solution2, TmpSolDg2, LagN, DgPerm, NumberOfDofNodes )
              Solution2 => TmpSolDg2
            END IF
            IF( ASSOCIATED(Solution3) ) THEN              
              CALL p2LagrangeSwapper( Mesh, Solution3, TmpSolDg3, LagN, DgPerm, NumberOfDofNodes )
              Solution3 => TmpSolDg3
            END IF            
          ELSE IF ( VarType == Variable_on_nodes_on_elements ) THEN
            IF( .NOT. ( DG .OR. DN ) ) CYCLE
            IF( .NOT. SaveElemental ) CYCLE
          ELSE IF( VarType == Variable_on_elements ) THEN
            CYCLE
          ELSE IF( VarType == Variable_on_gauss_points ) THEN
            IF ( .NOT. ( DG .OR. DN ) ) CYCLE                       

            CALL Ip2DgSwapper( Mesh, Solution, TmpSolDg, Variable_on_nodes_on_elements )
            IF(DN) CALL CalculateBodyAverage( Mesh, TmpSolDg, .FALSE. )

            Solution => TmpSolDg 
            IF(ASSOCIATED(Solution2)) THEN
              CALL Ip2DgSwapper( Mesh, Solution2, TmpSolDg2, Variable_on_nodes_on_elements )
              IF(DN) CALL CalculateBodyAverage( Mesh, TmpSolDg2, .FALSE. )
              Solution2 => TmpSolDg2 
            END IF
            IF(ASSOCIATED(Solution3)) THEN
              CALL Ip2DgSwapper( Mesh, Solution3, TmpSolDg3, Variable_on_nodes_on_elements )
              IF(DN) CALL CalculateBodyAverage( Mesh, TmpSolDg3, .FALSE. )
              Solution3 => TmpSolDg3 
            END IF
          END IF

          VarType = Solution % TYPE
          Values => Solution % Values
          IF ( ASSOCIATED(Solution2) ) Values2 => Solution2 % Values
          IF ( ASSOCIATED(Solution3) ) Values3 => Solution3 % Values

          IF( MaxModes == 0 ) THEN
            EigenVectors => NULL()
            EigenVectors2 => NULL()
            EigenVectors3 => NULL()            
          ELSE
            EigenVectors => Solution % EigenVectors
            IF(ASSOCIATED(Solution2) ) EigenVectors2 => Solution2 % EigenVectors
            IF(ASSOCIATED(Solution3) ) EigenVectors3 => Solution3 % EigenVectors
          END IF

          IF( MaxModes2 == 0 ) THEN
            ConstraintModes => NULL()
          ELSE
            ConstraintModes => Solution % ConstraintModes
          END IF
            
          ! Default is to save the field only once
          NoFields = 0
          NoFields2 = 0
          NoModes = 0
          NoModes2 = 0

          IsHarmonic = .FALSE.
          IF( ASSOCIATED( Solution % Solver ) ) THEN
            IsHarmonic = ListCheckPresent( Solution % Solver % Values, &
                'Harmonic System Values' )
          END IF

          IF( EigenAnalysis ) THEN
            IF( MaxModes > 0 .AND. FileIndex <= MaxModes .AND. &
                ASSOCIATED(EigenVectors) ) THEN  
              NoModes = SIZE( Solution % EigenVectors, 1 )

              IF( GotActiveModes ) THEN
                IndField = ActiveModes( FileIndex ) 
              ELSE
                IndField = FileIndex
              END IF
              IF( IndField > NoModes ) THEN
                WRITE( Message,'(A,I0,A,I0,A)') 'Too few eigenmodes (',&
                    IndField,',',NoModes,') in '//TRIM(FieldName)       
                CALL Warn(Caller,Message)
                CYCLE
              END IF
              NoModes = 1
              NoFields = 1
            ELSE IF( FileIndex > MaxModes .AND. &
                ASSOCIATED(ConstraintModes) ) THEN

              NoModes2 = Solution % NumberOfConstraintModes
              IF( GotActiveModes2 ) THEN
                IndField = ActiveModes2( FileIndex - MaxModes ) 
              ELSE
                IndField = FileIndex - MaxModes 
              END IF
              IF( IndField > NoModes2 ) THEN
                WRITE( Message,'(A,I0,A,I0,A)') 'Too few constraint modes (',&
                    IndField,',',NoModes,') in '//TRIM(FieldName)       
                CALL Warn(Caller,Message)
                CYCLE
              END IF
              NoModes2 = 1
              NoFields2 = 1
            END IF
          ELSE
            IF( MaxModes > 0 .AND. ASSOCIATED(Solution % EigenVectors) ) THEN  
              NoModes = SIZE( Solution % EigenVectors, 1 )
              IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
              NoFields = NoModes
            END IF

            IF( MaxModes2 > 0 .AND. ASSOCIATED(ConstraintModes) ) THEN
              NoModes2 = Solution % NumberOfConstraintModes
              IF( MaxModes2 > 0 ) NoModes2 = MIN( MaxModes2, NoModes2 )
              NoFields2 = NoModes2
            END IF

            IF( NoModes + NoModes2 == 0 ) NoFields = 1
          END IF
          
          Perm => Solution % Perm
          FlipActive = Solution % PeriodicFlipActive 
          
          !---------------------------------------------------------------------
          ! There may be special complementary variables such as 
          ! displacement & mesh update. 
          !---------------------------------------------------------------------
          ComplementExists = .FALSE.
          IF( .TRUE. ) THEN ! IF( NoModes + NoModes2 == 0 ) THEN
            IF(Rank==0) WRITE(Txt,'(A,I0,A)') 'Scalar Field ',Vari,' Complement'
            IF(Rank==1) WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
            IF(Rank==2) WRITE(Txt,'(A,I0,A)') 'Tensor Field ',Vari,' Complement'

            FieldNameB = GetString( Params, TRIM(Txt), Found )
            IF( Found ) THEN
              Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldNameB),ThisOnly=NoInterp)
              ComponentVectorB = .FALSE.
              IF(.NOT. ASSOCIATED( Solution ) ) THEN
                Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldNameB)//' 1',ThisOnly=NoInterp)
                ComponentVectorB = ASSOCIATED(Solution)
                IF(NoModes>0) EigenVectorsB => Solution % EigenVectors
              END IF
              
              IF( ASSOCIATED(Solution)) THEN 
                ValuesB => Solution % Values
                PermB => Solution % Perm 
                IF( ComponentVectorB ) THEN                  
                  Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldNameB)//' 2',ThisOnly=NoInterp)
                  IF( ASSOCIATED(Solution)) THEN
                    ValuesB2 => Solution % Values
                    IF(NoModes>0) EigenVectorsB2 => Solution % EigenVectors
                  END IF
                  Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldNameB)//' 3',ThisOnly=NoInterp)
                  IF( ASSOCIATED(Solution)) THEN
                    ValuesB3 => Solution % Values
                    IF(NoModes>0) EigenVectorsB3 => Solution % EigenVectors
                  END IF
                END IF
                ComplementExists = .TRUE.
              ELSE
                CALL Warn(Caller,'Complement does not exist:'//TRIM(FieldNameB))
              END IF
            END IF
          END IF
          
          IF( dofs > 1 ) THEN
            sdofs = MAX(dofs,dim)
          ELSE
            sdofs = 1
          END IF

          !---------------------------------------------------------------------
          ! Finally save the field values 
          !---------------------------------------------------------------------
          DoIm = 0
300       DO iField = 1, NoFields + NoFields2          

            IF( ( DG .OR. DN ) .AND. VarType == Variable_on_nodes_on_elements ) THEN
              CALL Info(Caller,'Setting field type to discontinuous',Level=12)
              InvFieldPerm => InvDgPerm
            ELSE IF( ALLOCATED( InvNodePerm ) ) THEN
              CALL Info(Caller,'Setting field type to nodal',Level=14)
              InvFieldPerm => InvNodePerm
            ELSE
              InvFieldPerm => NULL()
            END IF

            IF(.NOT. EigenAnalysis ) THEN
              IF( iField <= NoFields ) THEN
                IF( Nomodes > 0 ) THEN
                  IF( GotActiveModes ) THEN
                    IndField = ActiveModes( iField ) 
                  ELSE
                    IndField = iField
                  END IF
                END IF
              ELSE
                IF( Nomodes2 > 0 ) THEN
                  IF( GotActiveModes2 ) THEN
                    IndField = ActiveModes2( iField - NoFields ) 
                  ELSE
                    IndField = iField - NoFields
                  END IF
                END IF
              END IF
            END IF
            
            IF( WriteXML ) THEN
              NoFieldsWritten = NoFieldsWritten + 1

              IF( NoModes + NoModes2 == 0 .OR. EigenAnalysis ) THEN
                WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
              ELSE IF( iField <= NoFields ) THEN
                IF( IsHarmonic ) THEN
                  WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                      TRIM(FieldName)//' HarmonicMode',IndField
                ELSE
                  WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                      TRIM(FieldName)//' EigenMode',IndField
                END IF
              ELSE
                WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                    TRIM(FieldName)//' ConstraintMode',IndField
              END IF
              CALL AscBinStrWrite( OutStr )

              WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'          
              CALL AscBinStrWrite( OutStr ) 

              IF( AsciiOutput ) THEN
                WRITE( OutStr,'(A)') ' format="ascii">'//lf
                CALL AscBinStrWrite( OutStr ) 
              ELSE
                WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
                CALL AscBinStrWrite( OutStr ) 
              END IF
            END IF

            IF( BinaryOutput ) THEN
              k = NumberOfDofNodes * PrecSize * sdofs
              Offset = Offset + IntSize + k
            END IF

            !---------------------------------------------------------------------
            ! Data may also be appended and then its saved on the second sweep
            !---------------------------------------------------------------------
            IF( WriteData ) THEN

              IF( .NOT. NoPermutation .AND. NumberOfDofNodes > 0 ) THEN
                IF(.NOT. ASSOCIATED( InvFieldPerm ) ) THEN
                  CALL Fatal(Caller,'InvFieldPerm not associated!')
                END IF
              END IF
              
              IF( BinaryOutput ) WRITE( VtuUnit ) k

              DO ii = 1, NumberOfDofNodes

                IF( NoPermutation ) THEN
                  i = ii 
                ELSE
                  i = InvFieldPerm(ii) 
                END IF

                IF( ASSOCIATED( Perm ) .AND. LagN == 0 ) THEN
                  j = Perm(i)
                ELSE
                  j = i
                END IF

                Use2 = .FALSE.
                IF( ComplementExists ) THEN
                  IF( j == 0 ) THEN
                    Use2 = .TRUE. 
                    j = PermB(i)
                  END IF
                END IF
                
                DO k=1,sdofs              
                  IF(j==0 .OR. k > dofs) THEN
                    val = 0.0_dp
                  ELSE IF( NoModes > 0 .AND. iField <= NoFields ) THEN
                    IF( Use2 ) THEN
                      IF( ComponentVectorB ) THEN
                        IF( k == 1 ) zval = EigenVectorsB(IndField,j)
                        IF( k == 2 ) zval = EigenVectorsB2(IndField,j)
                        IF( k == 3 ) zval = EigenVectorsB3(IndField,j)
                      ELSE
                        zval = EigenVectorsB(IndField,dofs*(j-1)+k) 
                      END IF
                    ELSE
                      IF( ComponentVector ) THEN
                        IF( k == 1 ) zval = EigenVectors(IndField,j)
                        IF( k == 2 ) zval = EigenVectors2(IndField,j)
                        IF( k == 3 ) zval = EigenVectors3(IndField,j)
                      ELSE
                        zval = EigenVectors(IndField,dofs*(j-1)+k) 
                      END IF
                    END IF
                    val = PickComplex(zval,EigenVectorMode+DoIm) 

                  ELSE IF( NoModes2 > 0 ) THEN
                    val = ConstraintModes(IndField,dofs*(j-1)+k)
                  ELSE
                    IF( Use2 ) THEN
                      IF( ComponentVectorB ) THEN
                        IF( k == 1 ) val = ValuesB(j)
                        IF( k == 2 ) val = ValuesB2(j)
                        IF( k == 3 ) val = ValuesB3(j)
                      ELSE
                        val = ValuesB(dofs*(j-1)+k)              
                      END IF
                    ELSE
                      IF( ComponentVector ) THEN
                        IF( k == 1 ) val = Values(j)
                        IF( k == 2 ) val = Values2(j)
                        IF( k == 3 ) val = Values3(j)
                      ELSE
                        IF(dofs*(j-1)+k > SIZE(Values) .OR. dofs*(j-1)+k < 1 ) THEN
                          PRINT *,'vtu:',dofs,j,k,SIZE(values),dofs*(j-1)+k
                          call flush(6)
                        END IF
                        val = Values(dofs*(j-1)+k)              
                      END IF
                    END IF
                  END IF

                  IF( FlipActive ) THEN
                    IF( Model % Mesh % PeriodicFlip(i) ) val = -val
                  END IF
                  
                  CALL AscBinRealWrite( val )
                END DO
              END DO

              CALL AscBinRealWrite( 0.0_dp, .TRUE. )

            END IF

            IF( AsciiOutput ) THEN
              WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
              CALL AscBinStrWrite( OutStr ) 
            END IF
          END DO

          IF( NoModes > 0 ) THEN !.AND. iField <= NoFields ) THEN
            ! We have chosen to save both real and imaginary components.
            ! We have done real, now redo with im. This is later patch, hence the dirty GOTO. 
            IF(EigenVectorMode == 3 .AND. DoIm == 0 ) THEN
              DoIm = 1 
              CALL Info(Caller,'Doing the imaginary component of: '//TRIM(FieldName),Level=25)
              FieldName = TRIM(FieldName)//' Im'
              GOTO 300 
            END IF
          END IF
          
        END DO
      END DO
    END IF ! IF( SaveNodal )

    IF( WriteXML ) THEN
      CALL Info(Caller,'Number of nodal fields written: '//I2S(NoFieldsWritten),Level=10)
      WRITE( OutStr,'(A)') '      </PointData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    ! Elementwise information
    !-------------------------------------
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      <CellData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    IF( SaveElemental ) THEN
      CALL Info(Caller,'Writing elemental fields',Level=10)
      NoFieldsWritten = 0
      DO Rank = 0,1
        DO Vari = 1, 999
          
          IF( Rank == 0 ) THEN
            WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          ELSE
            WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          END IF
          FieldName = GetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT

          !---------------------------------------------------------------------
          ! Find the variable with the given name in the normal manner 
          !---------------------------------------------------------------------
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName),ThisOnly=NoInterp)
          ComponentVector = .FALSE.

          ! If we are looking for a vector just one dofs won't do!
          ! This circumvents a problem somewhere else in the code. 
          IF( ASSOCIATED( Solution ) ) THEN
            IF( Rank > 0 .AND. Solution % Dofs <= 1 ) NULLIFY( Solution ) 
          END IF

          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE 
              CALL Warn(Caller,'Nonexistent variable: '//TRIM(FieldName))
              CYCLE
            END IF
          END IF
          
          VarType = Solution % TYPE

          IF( DG .OR. DN ) THEN
            Found = ( VarType == Variable_on_elements )
          ELSE
            Found = ( VarType == Variable_on_nodes_on_elements .OR. &
                VarType == Variable_on_gauss_points  .OR. &
                VarType == Variable_on_elements )            
          END IF
          IF (.NOT. Found ) CYCLE
          
          Perm => Solution % Perm
          Dofs = Solution % DOFs
          Values => Solution % Values

          IF( Solution % PeriodicFlipActive ) THEN
            CALL Warn(Caller,'Cannot yet deal with PeriodicFlip in elemental variables!')
          END IF
          
          !---------------------------------------------------------------------
          ! Some vectors are defined by a set of components (either 2 or 3)
          !---------------------------------------------------------------------
          IF( ComponentVector ) THEN
            dofs = 1 
            Solution2 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution2)) THEN
              Values2 => Solution2 % Values
              dofs = 2
            END IF
            Solution3 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution3)) THEN
              Values3 => Solution3 % Values
              dofs = 3
            END IF
          END IF

          IF( MaxModes == 0 ) THEN
            EigenVectors => NULL()
            EigenVectors2 => NULL()
            EigenVectors3 => NULL()
          ELSE
            EigenVectors => Solution % EigenVectors
            IF(ASSOCIATED(Solution2) ) EigenVectors2 => Solution2 % EigenVectors
            IF(ASSOCIATED(Solution3) ) EigenVectors3 => Solution3 % EigenVectors
          END IF
            
          IF( MaxModes2 == 0 ) THEN
            ConstraintModes => NULL()
          ELSE
            ConstraintModes => Solution % ConstraintModes
          END IF
            
          ! Default is to save the field only once
          NoFields = 0
          NoFields2 = 0
          NoModes = 0
          NoModes2 = 0
          
          IF( EigenAnalysis ) THEN
            IF( MaxModes > 0 .AND. FileIndex <= MaxModes .AND. &
                ASSOCIATED(EigenVectors) ) THEN  
              NoModes = SIZE( EigenVectors, 1 )
              
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( FileIndex ) 
              ELSE
                IndField = FileIndex
              END IF
              IF( IndField > NoModes ) THEN
                WRITE( Message,'(A,I0,A,I0,A)') 'Too few eigenmodes (',&
                    IndField,',',NoModes,') in '//TRIM(FieldName)       
                CALL Warn(Caller,Message)
                CYCLE
              END IF
              NoModes = 1
              NoFields = 1
            ELSE IF( FileIndex > MaxModes .AND. &
                ASSOCIATED(ConstraintModes) ) THEN
              
              NoModes2 = Solution % NumberOfConstraintModes
              IF( GotActiveModes2 ) THEN
                IndField = ActiveModes2( FileIndex - MaxModes ) 
              ELSE
                IndField = FileIndex - MaxModes 
              END IF
              IF( IndField > NoModes2 ) THEN
                WRITE( Message,'(A,I0,A,I0,A)') 'Too few constraint modes (',&
                    IndField,',',NoModes,') in '//TRIM(FieldName)       
                CALL Warn(Caller,Message)
                CYCLE
              END IF
              NoModes2 = 1
              NoFields2 = 1
            END IF
          ELSE
            IF( MaxModes > 0 .AND. ASSOCIATED(EigenVectors) ) THEN  
              NoModes = SIZE( EigenVectors, 1 )
              IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
              NoFields = NoModes
            END IF
            
            IF( MaxModes2 > 0 .AND. ASSOCIATED(ConstraintModes) ) THEN
              NoModes2 = Solution % NumberOfConstraintModes
              IF( MaxModes2 > 0 ) NoModes2 = MIN( MaxModes2, NoModes2 )
              NoFields2 = NoModes2
            END IF
            
            IF( NoModes + NoModes2 == 0 ) NoFields = 1
          END IF

          IF( dofs > 1 ) THEN
            sdofs = MAX(dofs,dim)
          ELSE
            sdofs = 1
          END IF
          
          !---------------------------------------------------------------------
          ! Finally save the field values 
          !---------------------------------------------------------------------
          DoIm = 0
400       DO iField = 1, NoFields + NoFields2          

            IF(.NOT. EigenAnalysis ) THEN
              IF( iField <= NoFields ) THEN
                IF( Nomodes > 0 ) THEN
                  IF( GotActiveModes ) THEN
                    IndField = ActiveModes( iField ) 
                  ELSE
                    IndField = iField
                  END IF
                END IF
              ELSE
                IF( Nomodes2 > 0 ) THEN
                  IF( GotActiveModes2 ) THEN
                    IndField = ActiveModes2( iField - NoFields ) 
                  ELSE
                    IndField = iField - NoFields
                  END IF
                END IF
              END IF
            END IF

            IF( WriteXML ) THEN
              NoFieldsWritten = NoFieldsWritten + 1

              CALL Info(Caller,'Writing variable: '//TRIM(FieldName),Level=20)
              IF( NoModes + NoModes2 == 0 .OR. EigenAnalysis ) THEN
                WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
              ELSE IF( iField <= NoFields ) THEN
                IF( IsHarmonic ) THEN
                  WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                      TRIM(FieldName)//' HarmonicMode',IndField
                ELSE
                  WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                      TRIM(FieldName)//' EigenMode',IndField
                END IF
              ELSE
                WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                    TRIM(FieldName)//' ConstraintMode',IndField
              END IF
              CALL AscBinStrWrite( OutStr )

              WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'          
              CALL AscBinStrWrite( OutStr ) 

              IF( AsciiOutput ) THEN
                WRITE( OutStr,'(A)') ' format="ascii">'//lf
                CALL AscBinStrWrite( OutStr ) 
              ELSE
                WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
                CALL AscBinStrWrite( OutStr ) 
              END IF
            END IF
            
            IF( BinaryOutput ) THEN
              k = PrecSize * sdofs * NumberOfElements
              Offset = Offset + IntSize + k
            END IF
                        
            IF( WriteData ) THEN
              IF( BinaryOutput ) WRITE( VtuUnit ) k

              DO i = ElemFirst, ElemLast
                IF( .NOT. ActiveElem(i) ) CYCLE
                CurrentElement => Model % Elements(i)

                ElemVectVal = 0._dp
                ElemInd = 0

                IF( VarType == Variable_on_nodes_on_elements ) THEN

                  IF( SaveLinear ) THEN
                    n = GetElementCorners( CurrentElement )
                  ELSE
                    n = GetElementNOFNodes( CurrentElement )
                  END IF

                  IF ( ASSOCIATED(CurrentElement % BoundaryInfo) .AND. .NOT. &
                      ASSOCIATED(CurrentElement % DGIndexes) ) THEN

                    Parent => CurrentElement % BoundaryInfo % Left
                    IF (.NOT.ASSOCIATED(Parent) ) &
                        Parent => CurrentElement % BoundaryInfo % Right

                    IF ( ASSOCIATED(Parent) ) THEN
                      IF (ASSOCIATED(Parent % DGIndexes) ) THEN
                        DO j=1,n
                          DO k=1,Parent % TYPE % NumberOfNodes
                            IF(Currentelement % NodeIndexes(j) == Parent % NodeIndexes(k)) &
                                ElemInd(j) = Perm( Parent % DGIndexes(k) )
                          END DO
                        END DO
                      END IF
                    END IF
                  ELSE
                    ElemInd(1:n) = Perm( CurrentElement % DGIndexes(1:n) )
                  END IF

                ELSE IF( VarType == Variable_on_gauss_points ) THEN

                  m = CurrentElement % ElementIndex
                  IF( m < SIZE( Perm ) ) THEN
                    n = Perm(m+1)-Perm(m)
                  ELSE
                    n = 0
                  END IF

                  DO j=1,n
                    ElemInd(j) = Perm(m)+j
                  END DO

                ELSE IF( VarType == Variable_on_elements ) THEN

                  m = CurrentElement % ElementIndex
                  n = 1
                  ElemInd(1) = 0

                  IF( ASSOCIATED( Perm ) ) THEN                  
                    IF( m>SIZE( Perm ) ) THEN
                      j = 0
                      IF( ASSOCIATED( CurrentElement % BoundaryInfo ) ) THEN                                            
                        IF( ASSOCIATED( CurrentElement % BoundaryInfo % Left ) ) THEN
                          j = CurrentElement % BoundaryInfo % Left % ElementIndex
                        END IF
                        IF( j <= 0 ) THEN
                          IF( ASSOCIATED( CurrentElement % BoundaryInfo % Right ) ) THEN
                            j = CurrentElement % BoundaryInfo % Right % ElementIndex
                          END IF
                        END IF
                      END IF

                      IF( j == 0 ) THEN
                        CALL Fatal(Caller,'Cannot define parent cell index for element: '//I2S(m))
                      END IF
                      m = j
                    END IF
                    ElemInd(1) = Perm( m ) 
                  END IF
                END IF

                IF ( ALL(ElemInd(1:n) > 0)) THEN                    
                  DO k=1,sdofs
                    IF( k > dofs ) THEN
                      val = 0.0_dp

                    ELSE IF( NoModes > 0 .AND. iField <= NoFields ) THEN
                      IF( ComponentVector ) THEN
                        IF( k == 1 ) zval = SUM(EigenVectors(IndField,ElemInd(1:n)))/n
                        IF( k == 2 ) zval = SUM(EigenVectors2(IndField,ElemInd(1:n)))/n
                        IF( k == 3 ) zval = SUM(EigenVectors3(IndField,ElemInd(1:n)))/n
                      ELSE
                        zval = SUM(EigenVectors(IndField,dofs*(ElemInd(1:n)-1)+k))/n
                      END IF
                      val = PickComplex(zval,EigenVectorMode+DoIm) 

                    ELSE IF( NoModes2 > 0 ) THEN
                      val = ConstraintModes(IndField,dofs*(j-1)+k)

                    ELSE
                      IF(ComponentVector) THEN
                        IF (k==1) val = SUM(Values(ElemInd(1:n)))/n
                        IF (k==2) val = SUM(Values2(ElemInd(1:n)))/n
                        IF (k==3) val = SUM(Values3(ElemInd(1:n)))/n
                      ELSE
                        val = SUM(Values(dofs*(ElemInd(1:n)-1)+k))/n
                      END IF
                    END IF
                    ElemVectVal(k) = val
                  END DO
                END IF
                
                DO k=1,sdofs
                  CALL AscBinRealWrite( ElemVectVal(k) )
                END DO
              END DO
              
              CALL AscBinRealWrite( 0.0_dp, .TRUE. )

            END IF
            
            IF( AsciiOutput ) THEN
              WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
              CALL AscBinStrWrite( OutStr ) 
            END IF

          END DO
            
          IF( NoModes > 0 ) THEN
            IF(EigenVectorMode == 3 .AND. DoIm == 0 ) THEN
              DoIm = 1 
              FieldName = TRIM(FieldName)//' Im'
              GOTO 400 
            END IF
          END IF
          
        END DO
      END DO
      IF( WriteXML ) THEN
        CALL Info(Caller,'Number of elemental fields written: '//I2S(NoFieldsWritten),Level=10)
      END IF
    END IF  ! IF( SaveElemental )

    !---------------------------------------------------------------------
    ! If requested write the body and bc indexes
    !---------------------------------------------------------------------
    IF( WriteIds ) THEN
      IF( WriteXML ) THEN
        CALL Info(Caller,'Writing entity IDs for bodies and boundaries',Level=10)

        WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="GeometryIds"'
        CALL AscBinStrWrite( OutStr )

        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' format="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF

      IF( BinaryOutput ) THEN
        k = IntSize * NumberOfElements
        Offset = Offset + IntSize + k
      END IF

      IF( WriteData ) THEN        
        IF( BinaryOutput ) WRITE( VtuUnit ) k

        DO i = ElemFirst, ElemLast
          IF(.NOT. ActiveElem(i)) CYCLE

          CurrentElement => Model % Elements(i)

          IF( i <= Mesh % NumberOfBulkElements ) THEN
            j = CurrentElement % BodyId
            IF (j>=1 .AND. j<= SIZE(GeometryBodyMap)) j = GeometryBodyMap( j )     
          ELSE
            j = GetBCId( CurrentElement )
            IF ( j>=1 .AND. j<= SIZE(GeometryBCMap)) THEN
              j = GeometryBCMap( j )
            ELSE
              j = BCOffset
            END IF
          END IF

          CALL AscBinIntegerWrite( j )
        END DO
        CALL AscBinIntegerWrite( 0, .TRUE. )
      END IF

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      </CellData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    ! Coordinates of each point
    !-------------------------------------
    IF( WriteXML ) THEN
      CALL Info(Caller,'Writing coordinates for each used node',Level=10)
      WRITE( OutStr,'(A)') '      <Points>'//lf
      CALL AscBinStrWrite( OutStr ) 

      WRITE( OutStr,'(A,I0,A,I0,A)') '        <DataArray type="Float',&
          PrecBits,'" NumberOfComponents="',dim,'"'
      CALL AscBinStrWrite( OutStr )       

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF


    IF( BinaryOutput ) THEN
      k = dim * NumberOfDofNodes * PrecSize
      Offset = Offset + IntSize + k
    END IF

    ! For higher order L-element we need to create also the coordinates on-the-fly.
    ! It would be convenient if the coordinate would be available as a 3-component variable.
    ! If the coordinates change this should be modified...
    IF( LagN > 0 ) THEN
      GotIt = .FALSE.
      IF(.NOT. ASSOCIATED(LagVarX) ) THEN
        ALLOCATE( LagVarX, LagVarY, LagVarZ )
        GotIt = .TRUE.
      END IF
      ! We need to redo coordinates, if they are new or if they are displaced
      IF( ListGetLogicalAnySolver(Model,'Displace Mesh') ) GotIt = .TRUE.
      IF(GotIt) THEN      
        Solution => VariableGet( Mesh % Variables,'Coordinate 1')
        CALL p2LagrangeSwapper( Mesh, Solution, LagVarX, LagN, DgPerm, NumberOfDofNodes )
        Solution2 => VariableGet( Mesh % Variables,'Coordinate 2')
        CALL p2LagrangeSwapper( Mesh, Solution2, LagVarY, LagN, DgPerm, NumberOfDofNodes )
        Solution3 => VariableGet( Mesh % Variables,'Coordinate 3')
        CALL p2LagrangeSwapper( Mesh, Solution3, LagVarZ, LagN, DgPerm, NumberOfDofNodes )
      END IF
    END IF
    
    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k 

      DO ii = 1, NumberOfDofNodes
        IF( NoPermutation ) THEN
          i = ii 
        ELSE
          i = InvNodePerm(ii)
        END IF

        IF( LagN > 0 ) THEN
          x = LagVarX % Values(i)
          y = LagVarY % Values(i)
          z = LagVarZ % Values(i)
        ELSE          
          x = Model % Mesh % Nodes % x( i )
          y = Model % Mesh % Nodes % y( i )
          z = Model % Mesh % Nodes % z( i )
        END IF
          
        ! If displacement field is active remove the displacement from the coordinates
        IF( dispdofs > 0 .OR. dispBdofs > 0) THEN
          j = 0
          IF(dispdofs > 0) THEN
            j = DispPerm(i)
            IF( j > 0 ) THEN
              x = x - DispValues(dispdofs*(j-1)+1)
              y = y - DispValues(dispdofs*(j-1)+2)
              IF(dispdofs == 3) z = z - DispValues(dispdofs*(j-1)+3)
            END IF
          END IF
          IF(dispBdofs > 0 .AND. j==0) THEN
            j = DispBPerm(i)
            IF( j > 0 ) THEN
              x = x - DispBValues(dispBdofs*(j-1)+1)
              y = y - DispBValues(dispBdofs*(j-1)+2)
              IF(dispBdofs == 3) z = z - DispBValues(dispBdofs*(j-1)+3)
            END IF
          END IF
        END IF

        x = CoordScale(1) * x + CoordOffset(1)
        y = CoordScale(2) * y + CoordOffset(2)
        z = CoordScale(3) * z + CoordOffset(3)
        
        CALL AscBinRealWrite( x )
        CALL AscBinRealWrite( y )
        CALL AscBinRealWrite( z )
      END DO

      CALL AscBinRealWrite( 0.0_dp, .TRUE.)
    END IF

    IF( AsciiOutput ) THEN   
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      </Points>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    ! Write out the mesh
    !-------------------------------------
    IF( WriteXML ) THEN
      CALL Info(Caller,'Writing the elemental connectivity data',Level=10)
      WRITE( OutStr,'(A)') '      <Cells>'//lf
      CALL AscBinStrWrite( OutStr ) 

      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="connectivity"'
      CALL AscBinStrWrite( OutStr ) 

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      ! The offset needs to be summed over all elements, this is just the size factor
      k = 0
      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        IF( LagN > 0 ) THEN
          ! When we call without the indexes argument, we just get number of elemental dofs
          n = GetLagrangeIndexes( Mesh, LagN, CurrentElement )
        ELSE IF( SaveLinear ) THEN
          n = GetElementCorners( CurrentElement ) 
        ELSE
          n = GetElementNOFNodes( CurrentElement )
        END IF
        k = k + n * IntSize
      END DO
      Offset = Offset + k + IntSize
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k

      DO i = ElemFirst, ElemLast
        IF(.NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        
        IF( LagN > 0 ) THEN
          n = GetLagrangeIndexes( Mesh, LagN, CurrentElement, TmpIndexes ) 
        ELSE          
          CALL Elmer2VtkIndexes( CurrentElement, DG .OR. DN, SaveLinear, TmpIndexes )          
          IF( SaveLinear ) THEN
            n = GetElementCorners( CurrentElement ) 
          ELSE
            n = GetElementNOFNodes( CurrentElement )
          END IF
        END IF
          
        DO j=1,n
          IF( DN .OR. DG .OR. LagN > 0 ) THEN
            jj = DgPerm( TmpIndexes(j) )
          ELSE IF( NoPermutation ) THEN
            jj = TmpIndexes(j)
          ELSE
            jj = NodePerm( TmpIndexes(j) )
          END IF
          CALL AscBinIntegerWrite( jj - 1 )
        END DO
      END DO
      CALL AscBinIntegerWrite( 0, .TRUE. ) 
    END IF

    IF( AsciiOutput ) THEN
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF


    ! Offsets for element indexes 
    !-------------------------------------------------------------------

    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="offsets"'
      CALL AscBinStrWrite( OutStr ) 

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      k = NumberOfElements * IntSize
      Offset = Offset + IntSize + k
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k 

      cumn = 0
      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        IF( LagN > 0 ) THEN
          n = GetLagrangeIndexes( Mesh, LagN, CurrentElement )
        ELSE IF( SaveLinear ) THEN
          n = GetElementCorners( CurrentElement ) 
        ELSE
          n = GetElementNOFNodes( CurrentElement )
        END IF
        cumn = cumn + n

        CALL AscBinIntegerWrite( cumn )
      END DO

      CALL AscBinIntegerWrite( 0, .TRUE.)

    END IF


    IF( AsciiOutput ) THEN   
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="types"'
      CALL AscBinStrWrite( OutStr ) 

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' FORMAT="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      k = NumberOfElements * IntSize
      Offset = Offset + IntSize + k
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k

      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        j = CurrentElement % TYPE % ElementCode 
        ! If we have higher order Lagrange elements then the family remains the same but
        ! number of nodes increases. 
        IF( LagN > 0 ) THEN
          n = GetLagrangeIndexes( Mesh, LagN, CurrentElement )
          IF (n > 100) THEN
            j = 1000 * ( j / 100 ) + n
          ELSE
            j = 100 * ( j / 100 ) + n
          END IF
        END IF
        n = Elmer2VTKElement(j, SaveLinear )

        CALL AscBinIntegerWrite( n )
      END DO

      CALL AscBinIntegerWrite( 0, .TRUE. )     
    END IF

    IF( AsciiOutput ) THEN
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXml ) THEN
      WRITE( OutStr,'(A)') '      </Cells>'//lf
      CALL AscBinStrWrite( OutStr ) 
      WRITE( OutStr,'(A)') '    </Piece>'//lf
      CALL AscBinStrWrite( OutStr ) 
      WRITE( OutStr,'(A)') '  </UnstructuredGrid>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    IF( BinaryOutput ) THEN
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '<AppendedData encoding="raw">'//lf                    
        CALL AscBinStrWrite( OutStr )           
        WRITE( VtuUnit ) '_'

        WriteXML = .FALSE.
        WriteData = .TRUE.
        GOTO 100
      ELSE
        WRITE( OutStr,'(A)') lf//'</AppendedData>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    WRITE( OutStr,'(A)') '</VTKFile>'//lf
    CALL AscBinStrWrite( OutStr ) 

    WRITE( OutStr,'(A)') ' '
    CALL AscBinStrWrite( OutStr ) 

    CLOSE( VtuUnit )

    CALL AscBinWriteFree()

    IF( ALLOCATED( ElemInd ) ) DEALLOCATE(ElemInd)

    IF( ASSOCIATED( TmpSolDg ) ) THEN
      CALL ReleaseVariableList( TmpSolDg )
      IF(ASSOCIATED(TmpSolDg2)) CALL ReleaseVariableList( TmpSolDg2 )
      IF(ASSOCIATED(TmpSolDg3)) CALL ReleaseVariableList( TmpSolDg3 )      
    END IF
    
    CALL Info(Caller,'Finished writing file',Level=15)
    
  END SUBROUTINE WriteVtuFile


  ! Write collection file that may include timesteps and/or parts
  !--------------------------------------------------------------
  SUBROUTINE WritePvdFile( PvdFile, DataSetFile, nTime, Model )
    CHARACTER(LEN=*), INTENT(IN) :: PvdFile, DataSetFile
    INTEGER :: nTime, RecLen = 0
    TYPE(Model_t) :: Model     
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: n, nLine = 0, iostat
    REAL(KIND=dp) :: time
    CHARACTER :: lf
    CHARACTER(LEN=MAX_NAME_LEN) :: Str
    LOGICAL :: Found

    SAVE RecLen, nLine

    lf = CHAR(10)

    IF( ParEnv % PEs > 1 ) THEN
      IF( ParEnv % MyPE > 0 ) RETURN
    END IF
    time = GetTime()
    IF( GetLogical( Params,'Vtu time previous',Found) ) THEN
      time = time - GetTimestepSize()
    END IF
    
    IF( nLine == 0 ) THEN
      ! Find the maximum record length (modulo four)
      WRITE( Str,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'
      n = LEN_TRIM( Str ) 
      
      WRITE( Str,'(A,ES16.7,A,I0,A)') '<DataSet timestep="',time,&
        '" group="" part="',GroupId,'" file="'//TRIM(DataSetFile)//'"/>'
      n = MAX( LEN_TRIM( Str ), n ) 

      ! Just long enough
      RecLen = ((n/4)+5)*4
      
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='REPLACE', &
          ACCESS='DIRECT', ACTION='WRITE', RECL=RecLen, IOSTAT=iostat)
      IF( iostat /= 0 ) THEN
        CALL Fatal('WritePvdFile','Opening of file failed: '//TRIM(PvdFile))
      END IF
          
      IF ( LittleEndian() ) THEN
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'
      ELSE
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="BigEndian"><Collection>'
      END IF     
      nLine = 1
    ELSE
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='OLD', &
          ACCESS='DIRECT', ACTION='READWRITE', RECL=RecLen, IOSTAT=iostat)     
      IF( iostat /= 0 ) THEN
        CALL Fatal('WritePvdFile','Opening of file failed: '//TRIM(PvdFile))
      END IF
    END IF

    nLine = nLine + 1
    WRITE( VtuUnit,'(A,ES12.3,A,I0,A)',REC=nLine) lf//'<DataSet timestep="',time,&
        '" group="" part="',GroupId,'" file="'//TRIM(DataSetFile)//'"/>'
    WRITE( VtuUnit,'(A)',REC=nLine+1) lf//'</Collection></VTKFile>'

    CLOSE( VtuUnit )

    Visited = .TRUE.
    
  END SUBROUTINE WritePvdFile


  ! Write parallel holder for serial vtu files.
  !-----------------------------------------------------------------------------
  SUBROUTINE WritePvtuFile( PvtuFile, Model )
    CHARACTER(LEN=*), INTENT(IN) :: PVtuFile
    TYPE(Model_t) :: Model 
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: i,j,k,dofs,Rank,n,dim,vari,sdofs,iostat
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FullName
    LOGICAL :: ScalarsExist, VectorsExist, Found, ComponentVector, AllActive, ThisActive
    LOGICAL, POINTER :: ActivePartition(:)
    TYPE(Variable_t), POINTER :: Solution, Solution2, Solution3
    INTEGER :: Active, NoActive, ierr, NoFields, NoFields2, NoModes, NoModes2, IndField, iField, VarType
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    

    ThisActive = ( NumberOfGeomNodes > 0 ) 
   
    Active = 0
    IF( ThisActive ) Active = 1    	

    NoActive = ParallelReduction( Active ) 
    IF( NoActive == 0 ) THEN
      CALL Info('WritePvtuFile','No active partitions for saving',Level=12)
    END IF
    AllActive = ( NoActive == Partitions )
    
    WRITE( Message,'(A,I0,A,I0,A)') 'Number of active partitions is ',&
        NoActive,' (out of ',Partitions,')'
    CALL Info('WritePvtuFile',Message,Level=10)

    IF(.NOT. AllActive ) THEN
      IF( Part == 0 ) THEN
        ALLOCATE( ActivePartition(Partitions))
        ActivePartition = .FALSE.
        ActivePartition(1) = ThisActive 
      END IF

      DO i=2,Partitions
        IF( i == Part + 1 ) THEN
          CALL MPI_SEND(ThisActive, 1, MPI_LOGICAL, &
	       0, 1000, ELMER_COMM_WORLD, ierr )
        END IF
        IF( Part == 0 ) THEN
          CALL MPI_RECV( ActivePartition(i), 1, MPI_LOGICAL, &
	       i-1, 1000, ELMER_COMM_WORLD, status, ierr )
        END IF
      END DO
    END IF

    IF( Part > 0 ) RETURN
    CALL Info('WritePvtuFile','List of active partitions was composed',Level=12)

    OPEN( UNIT=VtuUnit, FILE=PvtuFile, form = 'formatted', STATUS='REPLACE', IOSTAT=iostat)
    IF( iostat /= 0 ) THEN
      CALL Fatal('WritePvtuFile','Opening of file failed: '//TRIM(PvtuFile))
    END IF
        
    dim = 3

    IF ( LittleEndian() ) THEN
      WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    ELSE
      WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="BigEndian">'
    END IF
    WRITE( VtuUnit,'(A)') '  <PUnstructuredGrid>'
    
    ! nodewise information
    !-------------------------------------
    ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
    VectorFieldName = GetString( Params,'Vector Field 1',VectorsExist)
    
    IF( SaveNodal ) THEN
      IF( ScalarsExist .AND. VectorsExist) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)&
            //'" Vectors="'//TRIM(VectorFieldName)//'">'
      ELSE IF( ScalarsExist ) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)//'">'
      ELSE IF( VectorsExist ) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Vectors="'//TRIM(VectorFieldName)//'">'
      ELSE
        CALL Warn('WritePvtuFile','Are there really no scalars or vectors?')
      END IF


      DO Rank = 0,2
        DO Vari = 1, 999
          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

          FieldName = GetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT

          IF(Rank == 2) THEN
            CALL Warn('WritePvtuFile','Do the tensors')
            EXIT
          END IF

          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName),ThisOnly=NoInterp)
          ComponentVector = .FALSE.

          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE
              CALL Warn('WritePvtuFile','Nonexistent variable: '//TRIM(FieldName)) 
              CYCLE
            END IF
          END IF

          VarType = Solution % TYPE

          IF ( VarType == Variable_on_nodes_on_elements ) THEN
            IF( .NOT. ( ( DG .OR. DN ) .AND. SaveElemental ) ) CYCLE
          ELSE IF( VarType == Variable_on_elements ) THEN
            CYCLE
          ELSE IF( VarType == Variable_on_gauss_points ) THEN
            IF ( DG ) THEN
              CONTINUE
            ELSE
              CYCLE
            END IF
          END IF

          ! Maybe we have some eigenmodes or constraint modes? 
          NoModes = 0
          NoModes2 = 0          
          IF( .NOT. EigenAnalysis ) THEN
            IF( MaxModes > 0 .AND. ASSOCIATED(Solution % EigenVectors) ) THEN  
              NoModes = SIZE( Solution % EigenVectors, 1 )
              IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
            END IF

            IF( MaxModes2 > 0 .AND. ASSOCIATED(Solution % ConstraintModes) ) THEN
              NoModes2 = Solution % NumberOfConstraintModes
              IF( MaxModes2 > 0 ) NoModes2 = MIN( MaxModes2, NoModes2 )              
            END IF
          END IF
          NoFields = MAX(1, MAX(NoModes, NoModes2 ) )
            
          dofs = Solution % DOFs
          IF( ComponentVector ) THEN
            Solution2 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution2)) dofs = 2
            Solution3 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly=NoInterp)
            IF( ASSOCIATED(Solution3)) dofs = 3
          END IF

          IF( dofs > 1 ) THEN
            sdofs = MAX(dofs,3)
          ELSE
            sdofs = 1
          END IF
          
          DO iField = 1, NoFields 
            IF( NoModes + NoModes2 == 0 .OR. EigenAnalysis ) THEN
              FullName = TRIM( FieldName ) 
            ELSE          
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( iField ) 
              ELSE
                IndField = iField
              END IF

              IF( NoModes > 0 ) THEN
                WRITE( FullName,'(A,I0)') TRIM( FieldName )//' EigenMode',IndField
              ELSE IF( NoModes2 > 0 ) THEN
                WRITE( FullName,'(A,I0)') TRIM( FieldName )//' ConstraintMode',IndField
              ELSE
                CALL Fatal('WritePvtuFile','Unknown case!')
              END IF
            END IF

            IF( AsciiOutput ) THEN
              WRITE( VtuUnit,'(A,A,I1,A)') '      <PDataArray type="Float64" Name="'//TRIM(FullName), &
                  '" NumberOfComponents="',sdofs,'" format="ascii"/>'  
            ELSE 
              WRITE( VtuUnit,'(A,I0,A,A,I0,A)') '      <PDataArray type="Float',&
                  PrecBits,'" Name="'//TRIM(FullName), &
                  '" NumberOfComponents="',sdofs,'" format="appended"/>'  
            END IF
          END DO

        END DO
      END DO

      IF( ScalarsExist .OR. VectorsExist) THEN
        WRITE( VtuUnit,'(A)') '    </PPointData>'
      END IF
    END IF


    ! Elementwise information
    !-------------------------------------
    WRITE( VtuUnit,'(A)') '    <PCellData>'

    IF( SaveElemental ) THEN
      IF( ScalarsExist .OR. VectorsExist ) THEN
        DO Rank = 0,2
          DO Vari = 1, 999
            IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
            FieldName = GetString( Params, TRIM(Txt), Found )
            IF(.NOT. Found ) EXIT

            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName),ThisOnly=NoInterp)
            ComponentVector = .FALSE.

            IF(ASSOCIATED(Solution)) THEN
              dofs = Solution % DOFs
            ELSE
              IF( Rank == 0 ) THEN
                CALL Warn('WritePvtuFile','Nonexistent scalar variable: '//TRIM(FieldName))
                CYCLE
              ELSE
                Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly=NoInterp)
                IF( ASSOCIATED(Solution)) THEN 
                  ComponentVector = .TRUE.
                  dofs = 1
                  Solution2 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly=NoInterp)
                  IF( ASSOCIATED(Solution2)) dofs = 2
                  Solution3 => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly=NoInterp)
                  IF( ASSOCIATED(Solution3)) dofs = 3                  
                ELSE 
                  CALL Warn('WritePvtuFile','Nonexistent vector variable: '//TRIM(FieldName))
                  CYCLE
                END IF
              END IF
            END IF

            VarType = Solution % TYPE

            IF( DG .OR. DN ) THEN
              Found = ( VarType == Variable_on_elements )
            ELSE
              Found = ( VarType == Variable_on_nodes_on_elements .OR. &
                  VarType == Variable_on_gauss_points  .OR. &
                  VarType == Variable_on_elements )            
            END IF
            IF (.NOT. Found ) CYCLE

            ! Maybe we have some eigenmodes or constraint modes? 
            NoModes = 0
            NoModes2 = 0          
            IF( .NOT. EigenAnalysis ) THEN
              IF( MaxModes > 0 .AND. ASSOCIATED(Solution % EigenVectors) ) THEN  
                NoModes = SIZE( Solution % EigenVectors, 1 )
                IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
              END IF

              IF( MaxModes2 > 0 .AND. ASSOCIATED(Solution % ConstraintModes) ) THEN
                NoModes2 = Solution % NumberOfConstraintModes
                IF( MaxModes2 > 0 ) NoModes2 = MIN( MaxModes2, NoModes2 )              
              END IF
            END IF
            NoFields = MAX(1, MAX(NoModes, NoModes2 ) )

            IF( dofs > 1 ) THEN
              sdofs = MAX(dofs,3)
            ELSE
              sdofs = 1
            END IF

            DO iField = 1, NoFields
              IF( NoModes + NoModes2 == 0 .OR. EigenAnalysis ) THEN
                FullName = TRIM( FieldName ) 
              ELSE          
                IF( NoModes > 0 ) THEN
                  WRITE( FullName,'(A,I0)') TRIM( FieldName )//' EigenMode',IndField
                ELSE IF( NoModes2 > 0 ) THEN
                  WRITE( FullName,'(A,I0)') TRIM( FieldName )//' ConstraintMode',IndField
                ELSE
                  CALL Fatal('WritePvtuFile','Unknown case!')
                END IF
              END IF

              IF( AsciiOutput ) THEN
                WRITE( VtuUnit,'(A,A,I1,A)') '      <PDataArray type="Float64" Name="'//TRIM(FullName), &
                    '" NumberOfComponents="',sdofs,'" format="ascii"/>'  
              ELSE 
                WRITE( VtuUnit,'(A,I0,A,A,I0,A)') '      <PDataArray type="Float',PrecBits,'" Name="'//TRIM(FullName), &
                    '" NumberOfComponents="',sdofs,'" format="appended"/>'  
              END IF
            END DO

          END DO
        END DO
      END IF
    END IF

    ! Body and BC indexes
    IF( WriteIds ) THEN
      IF( AsciiOutput ) THEN
        WRITE( VtuUnit,'(A)') '      <PDataArray type="Int32" Name="GeometryIds" format="ascii"/>'
      ELSE
        WRITE( VtuUnit,'(A)') '      <PDataArray type="Int32" Name="GeometryIds" format="appended"/>'        
      END IF
    END IF

    WRITE( VtuUnit,'(A)') '    </PCellData>'

    ! Coordinates of each point
    !-------------------------------------
    WRITE( VtuUnit,'(A)') '    <PPoints>'
    IF( AsciiOutput ) THEN
      WRITE( VtuUnit,'(A,I0,A)') '      <DataArray type="Float64" NumberOfComponents="',dim,'" format="ascii"/>'    
    ELSE 
      WRITE( VtuUnit,'(A,I0,A,I0,A)') '      <DataArray type="Float',PrecBits,&
          '" NumberOfComponents="',dim,'" format="appended"/>'    
    END IF
    WRITE( VtuUnit,'(A)') '    </PPoints>' 

 
    ! Write the pieces to the file 
    !-------------------------------------
    DO i=1,Partitions
      IF(.NOT. AllActive ) THEN
        IF( .NOT. ActivePartition(i)) CYCLE
      END IF      
      CALL VtuFileNaming( BaseFile, VtuFile,'.vtu', GroupId, FileIndex, i, NoPath = .TRUE.) 
      WRITE( VtuUnit,'(A)' ) '    <Piece Source="'//TRIM(VtuFile)//'"/>'
    END DO

    WRITE( VtuUnit,'(A)') '  </PUnstructuredGrid>'
    WRITE( VtuUnit,'(A)') '</VTKFile>'

    CLOSE( VtuUnit )

    IF(.NOT. AllActive ) DEALLOCATE( ActivePartition ) 
  
  END SUBROUTINE WritePvtuFile

!------------------------------------------------------------------------------
END SUBROUTINE VtuOutputSolver
!------------------------------------------------------------------------------
  
