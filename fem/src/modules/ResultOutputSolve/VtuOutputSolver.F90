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

!------------------------------------------------------------------------------
!> Subroutine for saving the results in XML based VTK format (VTU). Both ascii and binary
!> output is available, in single or double precision. The format is understood by 
!> visualization softwares Paraview and ViSit, for example.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VtuOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils 
  USE MeshUtils
  USE ElementDescription
  USE AscBinOutputUtils
  
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(dp) :: dt
  LOGICAL :: TransientSimulation
  
  INTEGER, SAVE :: nTime = 0
  LOGICAL :: GotIt, Hit, Parallel, FixedMesh, DG, DN
  CHARACTER(MAX_NAME_LEN) :: FilePrefix
  CHARACTER(MAX_NAME_LEN) :: BaseFile, VtuFile, PvtuFile, PvdFile, DataSetFile
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: Var
  INTEGER :: i, j, k, l, n, m, Partitions, Part, ExtCount, FileindexOffSet, MeshDim, PrecBits, &
             PrecSize, IntSize, FileIndex
  CHARACTER(MAX_NAME_LEN) :: Dir
  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp) :: DoubleWrk
  REAL :: SingleWrk

  LOGICAL :: MaskExists, BinaryOutput, AsciiOutput, SinglePrec, NoFileindex, &
      SkipHalo, SaveOnlyHalo, IsHalo, IsBoundaryElement
  CHARACTER(MAX_NAME_LEN) :: Str, MaskName
  TYPE(Variable_t), POINTER :: MaskVar
  INTEGER, POINTER :: MaskPerm(:), InvFieldPerm(:), NodeIndexes(:)
  INTEGER, ALLOCATABLE, TARGET :: NodePerm(:), InvNodePerm(:), InvDgPerm(:), DgPerm(:)
  INTEGER :: NumberOfGeomNodes, NumberOfDofNodes, NumberOfElements, ParallelNodes, ParallelElements, Sweep
  TYPE(Element_t), POINTER :: CurrentElement, LeftElem, RightElem
  TYPE(ValueList_t),POINTER :: Params
  INTEGER :: MaxModes, MaxModes2, BCOffset, ElemFirst, ElemLast, LeftIndex, RightIndex, &
      discontMesh, OutputMeshes
  INTEGER, POINTER :: ActiveModes(:), ActiveModes2(:), Indexes(:)
  LOGICAL :: GotActiveModes, GotActiveModes2, EigenAnalysis, ConstraintAnalysis, &
      WriteIds, SaveBoundariesOnly, SaveBulkOnly, SaveLinear, &
      GotMaskName, NoPermutation, SaveElemental, SaveNodal, GotMaskCond
  LOGICAL, ALLOCATABLE :: ActiveElem(:)
  INTEGER, ALLOCATABLE :: BodyVisited(:),GeometryBodyMap(:),GeometryBCMap(:)
  REAL(KIND=dp), ALLOCATABLE :: MaskCond(:)

! Parameters for buffered binary output
  INTEGER :: BufferSize


  Params => GetSolverParams()
  Mesh => Model % Mesh
  MeshDim = Mesh % MeshDim

  DG = GetLogical( Params,'Discontinuous Galerkin',GotIt)
  DN = GetLogical( Params,'Discontinuous Bodies',GotIt)

  ExtCount = GetInteger( Params,'Output Count',GotIt)
  IF( GotIt ) THEN
    nTime = ExtCount
  ELSE
    nTime = nTime + 1
  END IF

  FileIndexOffset = GetInteger( Params,'Fileindex offset',GotIt)
  FileIndex = nTime + FileIndexOffset

  BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
  IF( GotIt ) THEN
    AsciiOutput = .NOT. BinaryOutput
  ELSE
    AsciiOutput = GetLogical( Params,'Ascii Output',GotIt)
    BinaryOutput = .NOT. AsciiOutput
  END IF
  
  SaveElemental = GetLogical( Params,'Save Elemental Fields',GotIt)
  IF(.NOT. GotIt) SaveElemental = .TRUE.
    
  SaveNodal = GetLogical( Params,'Save Nodal Fields',GotIt) 
  IF(.NOT. GotIt) SaveNodal = .TRUE.

  SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
  IF( SinglePrec ) THEN
    CALL Info('VtuOutputSolver','Using single precision arithmetics in output!',Level=7)
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
    CALL Info('VtuOutputSolver','Saving results in VTK XML format with prefix: '//TRIM(FilePrefix))
    WRITE( Message,'(A,I0)') 'Saving number of partitions: ',Partitions
    CALL Info('VtuOutputSolver', Message )
  END IF


  BaseFile = FilePrefix
  IF ( .NOT. FileNameQualified(FilePrefix) ) THEN
    Dir = GetString( Params,'Output Directory',GotIt) 
    IF(.NOT. GotIt) Dir = GetString( Model % Simulation,&
        'Output Directory',GotIt)     
    IF( GotIt ) THEN
      IF( LEN_TRIM(Dir) > 0 ) THEN
        BaseFile = TRIM(Dir)// '/' //TRIM(FilePrefix)
        CALL MakeDirectory( TRIM(Dir) // CHAR(0) )
      END IF
    ELSE 
      BaseFile = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/' //TRIM(FilePrefix)
    END IF
  END IF
  CALL Info('VtuOutputSolver','Full filename base is: '//TRIM(Basefile), Level=10 )

  

  
  FixedMesh = ListGetLogical(Params,'Fixed Mesh',GotIt)

  
  !------------------------------------------------------------------------------
  ! Initialize stuff for masked saving
  !------------------------------------------------------------------------------
  ! Halo exists only in parallel
  IF( Parallel ) THEN
    SkipHalo = GetLogical( Params,'Skip Halo Elements', GotIt )
    SaveOnlyHalo = GetLogical( Params,'Save Halo Elements Only', GotIt )
  ELSE
    SkipHalo = .FALSE.
    SaveOnlyHalo = .FALSE.
  END IF

  GotMaskName = .FALSE.
  Str = GetString( Params,'Mask Variable',MaskExists)
  IF( MaskExists ) THEN
    MaskVar => VariableGet(Model % Variables,TRIM(Str))
    IF( ASSOCIATED(MaskVar)) MaskPerm => MaskVar % Perm
    MaskExists = ASSOCIATED(MaskPerm)
    IF( MaskExists ) THEN
      CALL Info('VtuOutputSolver','Using > '// TRIM(Str) // ' < as mask variable')
    END IF
  ELSE
    ! Check if there is an additional mask name given
    IF( MeshDim == 2 ) THEN
      MaskName = GetString( Params,'2D Mask Name',GotIt)    
    ELSE IF( MeshDim == 3 ) THEN  
      MaskName = GetString( Params,'3D Mask Name',GotIt)    
    END IF
    IF(.NOT. GotIt) MaskName = GetString( Params,'Mask Name',GotIt) 
    GotMaskName = GotIt
  END IF

  GotMaskCond = .FALSE.
  IF( .NOT. GotMaskName ) THEN
    MaskName = GetString( Params,'Mask Condition',GotMaskCond)
    IF( GotMaskCond ) THEN
      n = Mesh % MaxElementNodes
      ALLOCATE( MaskCond(n) )
    END IF
  END IF

  SaveBoundariesOnly = GetLogical( Params,'Save Boundaries Only',GotIt ) 
  SaveBulkOnly = GetLogical( Params,'Save Bulk Only',GotIt ) 
  
  NumberOfGeomNodes = Mesh % NumberOfNodes
  IF( MaskExists ) THEN
    NumberOfGeomNodes = COUNT( MaskPerm(1:NumberOfGeomNodes) > 0 ) 
  END IF
  NumberOfElements = 0

  IF( NumberOfGeomNodes > 0 ) THEN
    ElemFirst = HUGE( ElemFirst )
    ElemLast = 0 

    ALLOCATE(NodePerm(Mesh % NumberOfNodes))
    NodePerm = 0

    ALLOCATE(ActiveElem(Mesh % NumberOfBulkElements + & 
        Mesh % NumberOfBoundaryElements))
    ActiveElem = .FALSE.

    ! Count the true number of elements and mark the 1st and last element
    !-----------------------------------------------------------------------
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      
      IsBoundaryElement = ( i > Mesh % NumberOfBulkElements )

      IF( IsBoundaryElement ) THEN
        IF( SaveBulkOnly ) CYCLE
      ELSE
        IF( SaveBoundariesOnly ) CYCLE
      END IF
      
      CurrentElement => Mesh % Elements(i)
      Model % CurrentElement => CurrentElement

      IF( GetElementFamily( CurrentElement ) == 1 ) CYCLE          
      
      IF( SkipHalo .OR. SaveOnlyHalo ) THEN
        IF( IsBoundaryElement ) THEN
          IF( ASSOCIATED( CurrentElement % BoundaryInfo ) ) THEN
            LeftElem => CurrentElement % BoundaryInfo % Left
            IF( ASSOCIATED( LeftElem ) ) THEN
              LeftIndex = LeftElem % ElementIndex
              IF( LeftIndex > 0 ) THEN
                IF( Mesh % Elements(LeftIndex) % PartIndex /= ParEnv % MyPe ) LeftIndex = 0
              END IF
            ELSE
              LeftIndex = 0
            END IF
            RightElem => CurrentElement % BoundaryInfo % Right
            IF( ASSOCIATED( RightElem ) ) THEN
              RightIndex = RightElem % ElementIndex
              IF( RightIndex > 0 ) THEN
                IF( Mesh % Elements(RightIndex) % PartIndex /= ParEnv % MyPe ) RightIndex = 0
              END IF
            ELSE
              RightIndex = 0
            END IF
            IsHalo = ( LeftIndex == 0 .AND. RightIndex == 0 )
          ELSE
            IsHalo = .FALSE.
          END IF
        ELSE
          IsHalo = ( CurrentElement % PartIndex /= ParEnv % MyPe )
        END IF

        IF( IsHalo ) THEN
          IF( SkipHalo ) CYCLE
        ELSE
          IF( SaveOnlyHalo ) CYCLE
        END IF
      END IF


      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(CurrentElement % NodeIndexes) <= 0) ) CYCLE
      END IF

      IF( GotMaskName ) THEN
        Hit = .FALSE.
        IF( i <= Mesh % NumberOfBulkElements ) THEN
          l = CurrentElement % BodyId
          k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            Hit = ListGetLogical( Model % BodyForces(k) % Values, TRIM(MaskName), GotIt)
          END  IF
          IF( .NOT. Hit ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              Hit = ListGetLogical( Model % Equations(k) % Values, TRIM(MaskName), GotIt)
            END IF
          END IF
        ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
            Hit = ListGetLogical(Model % BCs(l) % Values, MaskName, GotIt ) 
            EXIT
          END DO
        END IF
        IF(.NOT. Hit ) CYCLE
      END IF

      IF( GotMaskCond ) THEN
        n = CurrentElement % TYPE % NumberOfNodes
        Indexes => CurrentElement % NodeIndexes

        IF( i <= Mesh % NumberOfBulkElements ) THEN
          l = CurrentElement % BodyId
          k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            MaskCond(1:n) = ListGetReal( Model % BodyForces(k) % Values, TRIM(MaskName), &
                n, Indexes, GotIt)
          END  IF

          IF( .NOT. Hit ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              MaskCond(1:n) = ListGetReal( Model % Equations(k) % Values, TRIM(MaskName), &
                  n, Indexes, GotIt)
            END IF
          END IF
        ELSE
          GotIt = .FALSE.
          IF( ASSOCIATED( CurrentElement % BoundaryInfo ) ) THEN
            DO l=1, Model % NumberOfBCs
              IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
              MaskCond(1:n) = ListGetReal(Model % BCs(l) % Values, MaskName, &
                  n, Indexes, GotIt ) 
              EXIT
            END DO
          END IF
        END IF
        IF( .NOT. GotIt ) CYCLE
        IF( .NOT. ALL(MaskCond(1:n) > 0.0_dp ) ) CYCLE
      END IF
      
      ActiveElem(i) = .TRUE.
      NumberOfElements = NumberOfElements + 1
      ElemFirst = MIN( ElemFirst, i )
      ElemLast = MAX( ElemLast, i )
      
      IF( SaveLinear ) THEN
        m = GetElementCorners( CurrentElement ) 
        NodePerm( CurrentElement % NodeIndexes(1:m) ) = 1
      ELSE          
        NodePerm( CurrentElement % NodeIndexes ) = 1
      END IF
      
    END DO
    
    CALL Info('VtuOutputSolver','Number of active elements '//TRIM(I2S(NumberOfElements))//&
        ' out of '//TRIM(I2S(Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements)),Level=10)

    NumberOfGeomNodes = COUNT( NodePerm > 0 ) 
    
    CALL Info('VtuOutputSolver','Number of geometry nodes '//TRIM(I2S(NumberOfGeomNodes))//&
        ' out of '//TRIM(I2S(Mesh % NumberOfNodes)),Level=10)
  END IF


  NumberOfDofNodes = 0

  IF( DG .OR. DN ) THEN    
    IF(.NOT. CheckAnyElementalField() ) THEN
      CALL Info('VtuOutputSolver','No elemental fields, omitting discontinuity creation!',Level=6)
      DG = .FALSE. 
      DN = .FALSE.
    END IF
  END IF



  ! If we have a discontinuous mesh then create the permutation vectors to deal with the discontinuities.
  IF( DG .OR. DN ) THEN
    NoPermutation = .FALSE.

    IF( DN ) THEN      
      CALL AverageBodyFields( Mesh )  
      ALLOCATE( BodyVisited( Mesh % NumberOfNodes ) )
    END IF

    k = 0
    DO i=1,Mesh % NumberOfBulkElements         
      CurrentElement => Mesh % Elements(i)
      k = k + CurrentElement % TYPE % NumberOfNodes
    END DO
    CALL Info('VtuOutputSolver','Maximum number of dofs in DG: '//TRIM(I2S(k)),Level=12)
    ALLOCATE( DgPerm(k) )
    DgPerm = 0

    DO Sweep=1,2
      l = 0
      IF( DG ) THEN
        DO i=1,Mesh % NumberOfBulkElements         
          IF( .NOT. ActiveElem(i) ) CYCLE
          CurrentElement => Mesh % Elements(i)
          NodeIndexes => CurrentElement % NodeIndexes

          IF( SaveLinear ) THEN
            m = GetElementCorners( CurrentElement )
          ELSE
            m = GetElementNOFNodes( CurrentElement )
          END IF

          DO k=1,m
            IF( NodePerm( NodeIndexes(k) ) == 0 ) CYCLE
            l = l + 1
            IF( Sweep == 2 ) THEN
              InvNodePerm(l) = NodeIndexes(k)
              DgPerm( CurrentElement % DGIndexes(k) ) = l
              InvDgPerm(l) = CurrentElement % DGIndexes(k)
            END IF
          END DO
        END DO
      ELSE      
        DO i=1,Model % NumberOfBodies
          BodyVisited = 0
          DO j=1,Mesh % NumberOfBulkElements         
            IF(.NOT. ActiveElem(i) ) CYCLE
            CurrentElement => Mesh % Elements(j)
            IF( CurrentElement % BodyId /= i ) CYCLE
            NodeIndexes => CurrentElement % NodeIndexes

            IF( SaveLinear ) THEN
              m = GetElementCorners( CurrentElement )
            ELSE
              m = GetElementNOFNodes( CurrentElement )
            END IF

            DO k=1,m
              IF( NodePerm( NodeIndexes(k) ) == 0 ) CYCLE
              IF( BodyVisited( NodeIndexes(k) ) > 0 ) THEN
                DgPerm( CurrentElement % DGIndexes(k) ) = BodyVisited( NodeIndexes(k) )
                CYCLE
              END IF
              l = l + 1
              BodyVisited(NodeIndexes(k)) = l
              IF( Sweep == 2 ) THEN
                InvNodePerm(l) = NodeIndexes(k)
                DgPerm( CurrentElement % DGIndexes(k) ) = l
                InvDgPerm(l) = CurrentElement % DGIndexes(k)
              END IF
            END DO
          END DO
        END DO
      END IF

      IF( Sweep == 1 ) THEN
        CALL Info('VtuOutputSolver','Independent dofs in discontinuous mesh: '//TRIM(I2S(l)),Level=10)
        NumberOfDofNodes = l
        ALLOCATE( InvNodePerm(l), InvDgPerm(l) ) 
        InvNodePerm = 0
        InvDgPerm = 0
      END IF
    END DO

    IF( DN ) DEALLOCATE( BodyVisited ) 

  ELSE
    NoPermutation = ( NumberOfGeomNodes == Mesh % NumberOfNodes )    
    IF( NoPermutation ) THEN
      DEALLOCATE( NodePerm ) 
    ELSE
      CALL Info('VtuOutputSolver','Not saving all nodes, creating permutation!',Level=12)
      ALLOCATE( InvNodePerm( NumberOfGeomNodes ) ) 
      InvNodePerm = 0
      j = 0
      DO i=1,Mesh % NumberOfNodes
        IF( NodePerm(i) > 0 ) THEN
          j = j + 1       
          NodePerm(i) = j
          InvNodePerm(j) = i
        END IF
      END DO
    END IF
    NumberOfDofNodes = NumberOfGeomNodes 
  END IF

  ! The partition is active for saving if there are any nodes 
  ! to write. There can be no elements nor dofs without nodes.
  CALL ParallelActive( NumberOfDofNodes > 0 )

  IF( nTime == 1 ) THEN
    ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfGeomNodes ) )
    WRITE( Message,'(A,I8)') 'Total number of geometry nodes to save:',ParallelNodes
    CALL Info('VtuOutputSolver',Message,Level=6)

    ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfDofNodes ) )
    WRITE( Message,'(A,I8)') 'Total number of dof nodes to save:',ParallelNodes
    CALL Info('VtuOutputSolver',Message,Level=6)

    ParallelElements = NINT( ParallelReduction( 1.0_dp * NumberOfElements ) )
    WRITE( Message,'(A,I8)') 'Total number of elements to save:',ParallelElements
    CALL Info('VtuOutputSolver',Message,Level=6)
  END IF

  IF( BinaryOutput ) THEN
    BufferSize = GetInteger( Params,'Binary Output Buffer Size',GotIt)
    IF( .NOT. GotIt ) BufferSize = MAX( NumberOfDofNodes, NumberOfElements )
  END IF



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
    CALL Info('VtuOutputSolver','Maximum number of eigen/harmonic modes: '//TRIM(I2S(MaxModes)),Level=7)
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
    CALL Info('VtuOutputSolver','Maximum number of constraint modes: '//TRIM(I2S(MaxModes2)),Level=7)
  END IF

  ! This activates the solution of the modes one for each file
  EigenAnalysis = ListGetLogical( Params,'Eigen Analysis',GotIt) .OR. &
      ListGetLogical( Params,'Constraint Modes Analysis',GotIt) 
  IF( EigenAnalysis ) THEN
    CALL Info('VtuOutputSolver','Saving each mode to different file')
    FileIndex = 1
  END IF

  BcOffset = 0
  WriteIds = GetLogical( Params,'Save Geometry Ids',GotIt)  
  IF( WriteIds ) THEN
    ! Create the mapping for body ids, default is unity mapping
    ALLOCATE( GeometryBodyMap( CurrentModel % NumberOfBodies ) )
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
    !PRINT *,'GeometryBodyMap:',GeometryBodyMap

    ! Create mapping for bc ids, default is unity mapping with offset
    ALLOCATE( GeometryBCMap( CurrentModel % NumberOfBCs ) )
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
          CALL Info('VtuOutputSolver','Setting offset for boundary entities: '&
              //TRIM(I2S(BCOffset)),Level=6)
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
    !PRINT *,'GeometryBcMap:',GeometryBcMap

  END IF
  

 100   CONTINUE

  IF(Parallel) THEN
    IF( NoFileindex ) THEN
      WRITE( PvtuFile,'(A,".pvtu")' ) TRIM(BaseFile)
    ELSE IF( FileIndex < 10000 ) THEN
      WRITE( PvtuFile,'(A,I4.4,".pvtu")' ) TRIM(BaseFile),FileIndex
    ELSE   
      WRITE( PvtuFile,'(A,I0,".pvtu")' ) TRIM(BaseFile),FileIndex
    END IF
    CALL Info('VtuOutputSolver','Writing the pvtu file: '//TRIM(PvtuFile), Level=10)
    CALL WritePvtuFile( PvtuFile, Model )
    CALL Info('VtuOutputSolver','Finished writing pvtu file',Level=12)
  END IF


  ! Write the Vtu file with all the data
  !--------------------------------------------------------------------------
  IF( NumberOfDofNodes > 0 ) THEN
    IF ( Parallel ) THEN
      IF( NoFileindex ) THEN
        WRITE( VtuFile,'(A,I4.4,A,".vtu")' ) TRIM(BaseFile),Part+1,"par"
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,A,I4.4,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            FileIndex
      ELSE
        WRITE( VtuFile,'(A,I4.4,A,I0,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            FileIndex
      END IF
    ELSE
      IF( NoFileindex ) THEN
        WRITE( VtuFile,'(A,".vtu")' ) TRIM(BaseFile)
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,".vtu")' ) TRIM(BaseFile),FileIndex
      ELSE
        WRITE( VtuFile,'(A,I0,".vtu")' ) TRIM(BaseFile),FileIndex
      END IF
    END IF

    CALL Info('VtuOutputSolver','Writing the vtu file: '//TRIM(VtuFile),Level=7)
    CALL WriteVtuFile( VtuFile, Model, FixedMesh )
    CALL Info('VtuOutputSolver','Finished writing vtu file',Level=12)
  END IF

  ! For transient simulation write a holder for the timesteps
  !-----------------------------------------------------------
  IF( GetLogical( Params,'Vtu Time Collection', GotIt ) ) THEN
    IF( TransientSimulation .AND. .NOT. NoFileIndex ) THEN
      WRITE( PvdFile,'(A,".pvd")' ) TRIM(BaseFile)
      IF( Parallel ) THEN
        DataSetFile = PvtuFile
      ELSE
        DataSetFile = VtuFile
      END IF
      CALL Info('VtuOutputSolver','Writing the pvd file: '//TRIM(DataSetFile),Level=10)
      CALL WritePvdFile( PvdFile, DataSetFile, FileIndex, Model )
      CALL Info('VtuOutputSolver','Finished writing pvd file',Level=12)     
    END IF
  END IF


  IF( EigenAnalysis ) THEN
    FileIndex = FileIndex + 1
    IF( FileIndex <= MaxModes + MaxModes2 ) GOTO 100
  END IF

  IF( NumberOfDofNodes > 0 ) THEN
    IF( .NOT. NoPermutation ) THEN
      DEALLOCATE( InvNodePerm, NodePerm ) 
    END IF
    DEALLOCATE( ActiveElem ) 
  END IF

  IF( WriteIds ) THEN  
    DEALLOCATE( GeometryBodyMap, GeometryBcMap )
  END IF
 
  
  CALL Info('VtuOutputSolver','All done for now',Level=10)     


CONTAINS



  FUNCTION Elmer2VtkElement( ElmerCode, SaveLinear ) RESULT ( VTKCode )
    INTEGER :: ElmerCode
    LOGICAL :: SaveLinear
    INTEGER :: VTKCode
    
    SELECT CASE (ElmerCode)
    CASE( 101 )
      VTKCode = 1
    CASE( 202 )
      VTKCode = 3
    CASE( 203 )
      VTKCode = 21
    CASE( 303 )
      VTKCode = 5
    CASE( 306 )
      VTKCode = 22
    CASE( 404 )
      VTKCode = 9
    CASE( 408 )
      VTKCode = 23
    CASE( 409 )
      VTKCode = 28
    CASE( 504 )
      VTKCode = 10
    CASE( 510 )
      VTKCode = 24
    CASE( 605 )
      VTKCode = 14
    CASE( 613 )
      VTKCode = 27
    CASE( 706 )
      VTKCode = 13
    CASE( 715 ) 
      VTKCode = 26
    CASE( 808 )
      VTKCode = 12
    CASE( 820 )
      VTKCode = 25
    CASE( 827 )
      VTKCode = 29
    CASE DEFAULT
      WRITE(Message,'(A,I0)') 'Not implemented for elementtype: ',ElmerCode
      CALL Fatal('Elmer2VtkElement',Message)
      
    END SELECT


    ! If requested return the 1st order element corresponding to the higher order elements
    IF( SaveLinear ) THEN
      SELECT CASE (VTKCode)
      CASE( 21 )
        VTKCode = 3
      CASE( 22 )
        VTKCode = 5
      CASE( 23, 28 )
        VTKCode = 9
      CASE( 24 )
        VTKCode = 10
      CASE( 27 )
        VTKCode = 14
      CASE( 26 )
        VTKCode = 13
      CASE( 25, 29 )
        VTKCode = 12
      END SELECT
    END IF

  END FUNCTION Elmer2VtkElement


!  FUNCTION Elmer2VtkIndexes( Element, DgElem, SaveLinear ) RESULT ( NodeIndexes )
  SUBROUTINE Elmer2VtkIndexes( Element, DgElem, SaveLinear, NodeIndexes )
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: DgElem
    LOGICAL :: SaveLinear
    INTEGER :: NodeIndexes(:)

    TYPE(Element_t), POINTER :: Parent
    INTEGER, POINTER :: UseIndexes(:)
    INTEGER, TARGET :: NewIndexes(27),BCIndexes(27)
    INTEGER :: ElmerCode, i,j,k,n,hits
    INTEGER, POINTER :: Order(:)
    INTEGER, TARGET, DIMENSION(20) :: &
        Order820 = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16/)
    INTEGER, TARGET, DIMENSION(27) :: &
        Order827 = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16,24,22,21,23,25,26,27/)
    LOGICAL :: DoReorder
      

    ElmerCode = Element % Type % ElementCode


    IF( DGElem ) THEN
      UseIndexes => NULL()
      IF( ASSOCIATED( Element % DGIndexes ) ) THEN
        UseIndexes => Element % DGIndexes
      ELSE IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
        Parent => Element % BoundaryInfo % Left
        IF (.NOT.ASSOCIATED(Parent) ) THEN
          Parent => Element % BoundaryInfo % Right        
        END IF
        IF ( ASSOCIATED(Parent) ) THEN
          IF (ASSOCIATED(Parent % DGIndexes) ) THEN
            n = Element % TYPE % NumberOfNodes 
            hits = 0
            DO j=1,n
              DO k=1,Parent % TYPE % NumberOfNodes
                IF(Element % NodeIndexes(j) == Parent % NodeIndexes(k)) THEN
                  BCIndexes(j) = Parent % DGIndexes(k) 
                  hits = hits + 1
                  EXIT
                END IF
              END DO
            END DO
            UseIndexes => BCIndexes
            IF( Hits < n ) THEN
              CALL Fatal('VtuOutputSolver','Could not determine DG boundary indexes')
            END IF
          END IF
        END IF
      ENDIF

      IF(.NOT. ASSOCIATED( UseIndexes ) ) THEN
        PRINT *,'Problematic BC elem:',Element % BodyId, Element % ElementIndex, Element % NodeIndexes, &
            ASSOCIATED( Element % DgIndexes ), ASSOCIATED( Element % BoundaryInfo ), DGelem, &
            Element % TYPE % ElementCode
        CALL Fatal('VtuOutputSolver','Could not set indexes for boundary element!')        
      END IF
    ELSE
      UseIndexes => Element % NodeIndexes
    END IF

    n = Element % TYPE % NumberOfNodes 


    ! Linear elements never require reordering 
    IF( .NOT. SaveLinear ) THEN
      SELECT CASE (ElmerCode)
        
      CASE( 820 )
        Order => Order820
        DoReOrder = .TRUE.
        
      CASE( 827 ) 
        Order => Order827
        DoReOrder = .TRUE.

      CASE DEFAULT
        DoReorder = .FALSE.
        
      END SELECT
    ELSE
      DoReOrder = .FALSE.
    END IF

    IF( DoReorder ) THEN
      NodeIndexes(1:n) = UseIndexes( Order(1:n) )
    ELSE
      NodeIndexes(1:n) = UseIndexes(1:n)
    END IF
    
  
  END SUBROUTINE Elmer2VtkIndexes




  ! Check whether there is any elemental field to be saved. 
  ! It does not make sense to use discontinuous saving if there are no discontinuous fields.
  ! It will even result to errors since probably there are no DG indexes either. 
  FUNCTION CheckAnyElementalField() RESULT ( HaveAnyElemental ) 

    LOGICAL :: HaveAnyElemental
    INTEGER :: Rank, Vari, VarType
    CHARACTER(LEN=1024) :: Txt, FieldName
    TYPE(Variable_t), POINTER :: Solution
    LOGICAL :: Found
    
    HaveAnyElemental = .FALSE.

    DO Rank = 0,1
      DO Vari = 1, 999
        IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        
        FieldName = GetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
        END IF
        IF( .NOT. ASSOCIATED( Solution ) ) CYCLE

        VarType = Solution % Type
        
        IF ( VarType == Variable_on_nodes_on_elements .OR. &
            VarType == Variable_on_elements .OR. &
            VarType == Variable_on_gauss_points ) THEN
          HaveAnyElemental = .TRUE.
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION CheckAnyElementalField



  ! Average fields within bodies
  SUBROUTINE AverageBodyFields( Mesh ) 
    
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(Variable_t), POINTER :: Var
    INTEGER :: NoAve
    LOGICAL :: BodySum 

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

    CALL Info('VtuOutputSolver','Reduced '//TRIM(I2S(NoAve))//' elemental fields',Level=7)

  END SUBROUTINE AverageBodyFields



  ! Writes a single VTU file that can be read by Paraview, ViSiT etc.
  !---------------------------------------------------------------------------------------
  SUBROUTINE WriteVtuFile( VtuFile, Model, RemoveDisp )
    CHARACTER(LEN=*), INTENT(IN) :: VtuFile
    TYPE(Model_t) :: Model 
    LOGICAL, INTENT(IN) :: RemoveDisp
    INTEGER, PARAMETER :: VtuUnit = 58
    TYPE(Variable_t), POINTER :: Var,Var1
    CHARACTER(LEN=512) :: str
    INTEGER :: i,ii,j,jj,k,dofs,Rank,cumn,n,m,dim,vari,sdofs,dispdofs, disp2dofs, Offset, &
        NoFields, NoFields2, IndField, iField, NoModes, NoModes2, NoFieldsWritten
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FieldName2, OutStr
    CHARACTER :: lf
    LOGICAL :: ScalarsExist, VectorsExist, Found,&
        ComponentVector, ComplementExists, Use2, IsHarmonic
    LOGICAL :: WriteData, WriteXML, L, Buffered
    TYPE(Variable_t), POINTER :: Solution
    INTEGER, POINTER :: Perm(:), Perm2(:), DispPerm(:), Disp2Perm(:)
    REAL(KIND=dp), POINTER :: Values(:), DispValues(:), Disp2Values(:), Values2(:), Values3(:)
    REAL(KIND=dp) :: x,y,z, val,ElemVectVal(3)
    INTEGER, ALLOCATABLE, TARGET :: ElemInd(:)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: TmpIndexes(27), VarType

    COMPLEX(KIND=dp), POINTER :: EigenVectors(:,:)
    REAL(KIND=dp), POINTER :: ConstraintModes(:,:)
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: CurrentElement, Parent
    TYPE(ValueList_t), POINTER :: Params


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

    ALLOCATE( ElemInd(Model % Mesh % MaxElementDOFS))

    ! This is a hack to ensure that the streamed saving will cover the whole file
    !----------------------------------------------------------------------------
    IF(.TRUE.) THEN
      OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'formatted', STATUS='unknown' )
      WRITE( VtuUnit,'(A)') ' '
      CLOSE( VtuUnit ) 
    END IF


    ! This format works both for ascii and binary output
    !-------------------------------------------------------------------------
    OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream', STATUS='unknown' )

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
      CALL Warn('WriteVtuFile','Are there really no scalars or vectors?')
    END IF

    WRITE( OutStr,'(A)') '      <PointData>'//lf
    CALL AscBinStrWrite( OutStr )

    DispDofs = 0
    Disp2Dofs = 0
    IF(RemoveDisp) THEN
      Solution => VariableGet( Model % Mesh % Variables, 'Displacement')
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

      Solution => VariableGet( Model % Mesh % Variables, 'Mesh Update')
      IF( ASSOCIATED( Solution ) ) THEN
        Disp2Perm => Solution % Perm
        Disp2Values => Solution % Values
        Disp2Dofs = Solution % Dofs
        CALL Info('VtuOutputSolver','Automatically complement > Displacement < by > Mesh Update < field',Level=7)
      END IF
    END IF


    ! When the data is 'appended' two loops will be taken and the data will be written
    ! on the second loop. Offset is the position in the appended data after the '_' mark.
    !------------------------------------------------------------------------------------
100 Offset = 0

    IF( SaveNodal ) THEN
      CALL Info('VtuOutputSolver','Writing nodal fields',Level=10)
      NoFieldsWritten = 0
      DO Rank = 0,2
        DO Vari = 1, 999
          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

          FieldName = GetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT

          IF(Rank == 2) THEN
            CALL Fatal('VtuOutputSolver','Do the tensors')
          END IF

          !---------------------------------------------------------------------
          ! Find the variable with the given name in the normal manner 
          !---------------------------------------------------------------------
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
          ComponentVector = .FALSE.
          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE
              WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
              CALL Warn('WriteVtuXMLFile', Txt)
              CYCLE
            END IF
          END IF

          CALL Info('VtuOutputSolver','Saving variable: '//TRIM(FieldName),Level=10)
          
          VarType = Solution % Type

          IF ( VarType == Variable_on_nodes_on_elements ) THEN
            IF( .NOT. ( ( DG .OR. DN ) .AND. SaveElemental ) ) CYCLE
          ELSE IF( VarType == Variable_on_elements ) THEN
            CYCLE
          ELSE IF( VarType == Variable_on_gauss_points ) THEN
            CYCLE
          END IF

          ! Default is to save the field only once
          NoFields = 0
          NoFields2 = 0
          NoModes = 0
          NoModes2 = 0

          EigenVectors => Solution % EigenVectors
          ConstraintModes => Solution % ConstraintModes

          IsHarmonic = .FALSE.
          IF( ASSOCIATED( Solution % Solver ) ) THEN
            IsHarmonic = ListCheckPresent( Solution % Solver % Values, &
                'Harmonic System Values' )
          END IF

          IF( EigenAnalysis ) THEN
            IF( MaxModes > 0 .AND. FileIndex <= MaxModes .AND. &
                ASSOCIATED(EigenVectors) ) THEN  
              NoModes = SIZE( Solution % EigenValues )

              IF( GotActiveModes ) THEN
                IndField = ActiveModes( FileIndex ) 
              ELSE
                IndField = FileIndex
              END IF
              IF( IndField > NoModes ) THEN
                WRITE( Message,'(A,I0,A,I0,A)') 'Too few eigenmodes (',&
                    IndField,',',NoModes,') in '//TRIM(FieldName)       
                CALL Warn('WriteVtuXMLFile',Message)
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
                CALL Warn('WriteVtuXMLFile',Message)
                CYCLE
              END IF
              NoModes2 = 1
              NoFields2 = 1
            END IF
          ELSE
            IF( MaxModes > 0 .AND. ASSOCIATED(Solution % EigenVectors) ) THEN  
              NoModes = SIZE( Solution % EigenValues )
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
          dofs = Solution % DOFs
          Values => Solution % Values
          VarType = Solution % Type
          
          !---------------------------------------------------------------------
          ! Some vectors are defined by a set of components (either 2 or 3)
          !---------------------------------------------------------------------
          IF( ComponentVector ) THEN

            IF( NoModes + NoModes2 > 0 ) THEN
              CALL Warn('WriteVtuXMLFile','Modes cannot currently be given componentwise!')
              CYCLE
            END IF
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
            IF( ASSOCIATED(Solution)) THEN
              Values2 => Solution % Values
              dofs = 2
            END IF
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
            IF( ASSOCIATED(Solution)) THEN
              Values3 => Solution % Values
              dofs = 3
            END IF
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
          END IF

          !---------------------------------------------------------------------
          ! There may be special complementary variables such as 
          ! displacement & mesh update. These are not implemented for modal output. 
          !---------------------------------------------------------------------
          ComplementExists = .FALSE.
          IF( NoModes + NoModes2 == 0 ) THEN
            IF(Rank==0) WRITE(Txt,'(A,I0,A)') 'Scalar Field ',Vari,' Complement'
            IF(Rank==1) WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
            IF(Rank==2) WRITE(Txt,'(A,I0,A)') 'Tensor Field ',Vari,' Complement'

            FieldName2 = GetString( Params, TRIM(Txt), Found )
            IF( Found ) THEN
              Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName2))
              IF( ASSOCIATED(Solution)) THEN 
                Values2 => Solution % Values
                Perm2 => Solution % Perm 
                ComplementExists = .TRUE.
              ELSE
                CALL Warn('WriteVTUFile','Complement does not exist:'//TRIM(FieldName2))
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
          DO iField = 1, NoFields + NoFields2          

            IF( ( DG .OR. DN ) .AND. VarType == Variable_on_nodes_on_elements ) THEN
              CALL Info('WriteVTUFile','Setting field type to discontinuous',Level=12)
              InvFieldPerm => InvDgPerm
            ELSE
              CALL Info('WriteVTUFile','Setting field type to nodal',Level=14)
              InvFieldPerm => InvNodePerm
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

              IF( BinaryOutput ) WRITE( VtuUnit ) k

              DO ii = 1, NumberOfDofNodes

                IF( NoPermutation ) THEN
                  i = ii 
                ELSE
                  i = InvFieldPerm(ii) 
                END IF

                IF( ASSOCIATED( Perm ) ) THEN
                  j = Perm(i)
                ELSE
                  j = i
                END IF

                Use2 = .FALSE.
                IF( ComplementExists ) THEN
                  IF( j == 0 ) THEN
                    Use2 = .TRUE. 
                    j = Perm2(i)
                  END IF
                END IF

                DO k=1,sdofs              
                  IF(j==0 .OR. k > dofs) THEN
                    val = 0.0_dp
                  ELSE IF( ComponentVector ) THEN
                    IF( k == 1 ) val = Values(j)
                    IF( k == 2 ) val = Values2(j)
                    IF( k == 3 ) val = Values3(j)
                  ELSE IF( Use2 ) THEN
                    val = Values2(dofs*(j-1)+k)              
                  ELSE IF( NoModes > 0 .AND. iField <= NoFields ) THEN
                    val = EigenVectors(IndField,dofs*(j-1)+k)                              
                  ELSE IF( NoModes2 > 0 ) THEN
                    val = ConstraintModes(IndField,dofs*(j-1)+k)
                  ELSE
                    val = Values(dofs*(j-1)+k)              
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
        END DO
      END DO
    END IF ! IF( SaveNodal )

    IF( WriteXML ) THEN
      CALL Info('VtuOutputSolver','Number of nodal fields written: '//TRIM(I2S(NoFieldsWritten)),Level=10)
      WRITE( OutStr,'(A)') '      </PointData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    ! Elementwise information
    !-------------------------------------
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      <CellData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    IF( SaveElemental .AND. .NOT. ( DG .OR. DN ) ) THEN
      CALL Info('VtuOutputSolver','Writing elemental fields',Level=10)
      NoFieldsWritten = 0
      DO Rank = 0,1
        DO Vari = 1, 999

          IF( Rank == 0 ) THEN
            WRITE(Txt,'(A,I0)') 'Scalar Field Elemental ',Vari
          ELSE
            WRITE(Txt,'(A,I0)') 'Vector Field Elemental ',Vari
          END IF
          FieldName = GetString( Params, TRIM(Txt), Found )
          L = Found

          IF(.NOT. Found) THEN
            IF( Rank == 0 ) THEN
              WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            ELSE
              WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
            END IF
            FieldName = GetString( Params, TRIM(Txt), Found )
          END IF

          IF(.NOT. Found) EXIT

          !---------------------------------------------------------------------
          ! Find the variable with the given name in the normal manner 
          !---------------------------------------------------------------------
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
          ComponentVector = .FALSE.

          ! If we are looking for a vector just one dofs wont do!
          ! This circumvents a problem somewhere else in the code. 
          IF( ASSOCIATED( Solution ) ) THEN
            IF( Rank > 0 .AND. Solution % Dofs <= 1 ) NULLIFY( Solution ) 
          END IF

          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE 
              IF( L ) THEN
                WRITE(Txt, '(A,A)') 'Nonexistent elemental variable: ',TRIM(FieldName)
                CALL Warn('WriteVtuXMLFile', Txt)
              END IF
              CYCLE
            END IF
          END IF

          VarType = Solution % Type
          Found = ( VarType == Variable_on_nodes_on_elements .OR. &
              VarType == Variable_on_gauss_points .OR. &
              VarType == Variable_on_elements )
          IF (.NOT. Found ) CYCLE

          Perm => Solution % Perm
          Dofs = Solution % DOFs
          Values => Solution % Values
          
          !---------------------------------------------------------------------
          ! Some vectors are defined by a set of components (either 2 or 3)
          !---------------------------------------------------------------------
          IF( ComponentVector ) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
            IF( ASSOCIATED(Solution)) THEN
              Values2 => Solution % Values
              dofs = 2
            END IF
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
            IF( ASSOCIATED(Solution)) THEN
              Values3 => Solution % Values
              dofs = 3
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
          IF( WriteXML ) THEN
            CALL Info('WriteVtuFile','Writing variable: '//TRIM(FieldName) )
            WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
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
            NoFieldsWritten = NoFieldsWritten + 1
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

                IF ( ALL(ElemInd(1:n) > 0)) THEN
                  IF( sdofs == 1 ) THEN
                    ElemVectVal(1) = SUM(Values(ElemInd(1:n))) / n
                  ELSE
                    DO k=1,sdofs
                      IF( k > dofs ) THEN
                        ElemVectVal(k) = 0.0_dp
                      ELSE IF(ComponentVector) THEN
                        IF (k==1) ElemVectVal(k) = SUM(Values(ElemInd(1:n)))/n
                        IF (k==2) ElemVectVal(k) = SUM(Values2(ElemInd(1:n)))/n
                        IF (k==3) ElemVectVal(k) = SUM(Values3(ElemInd(1:n)))/n
                      ELSE
                        ElemVectVal(k) = SUM(Values(dofs*(ElemInd(1:n)-1)+k))/n
                      END IF
                    END DO
                  END IF
                END IF 
                
              ELSE IF( VarType == Variable_on_gauss_points ) THEN

                m = CurrentElement % ElementIndex
                n = Perm(m+1)-Perm(m)
                IF( n > 0 ) THEN
                  DO j=1,n
                    ElemInd(j) = Perm(m)+j
                  END DO
                  
                  IF( sdofs == 1 ) THEN
                    ElemVectVal(1) = SUM(Values(ElemInd(1:n))) / n
                  ELSE
                    DO k=1,sdofs
                      IF( k > dofs ) THEN
                        ElemVectVal(k) = 0.0_dp
                      ELSE IF(ComponentVector) THEN
                        IF (k==1) ElemVectVal(k) = SUM(Values(ElemInd(1:n)))/n
                        IF (k==2) ElemVectVal(k) = SUM(Values2(ElemInd(1:n)))/n
                        IF (k==3) ElemVectVal(k) = SUM(Values3(ElemInd(1:n)))/n
                      ELSE
                        ElemVectVal(k) = SUM(Values(dofs*(ElemInd(1:n)-1)+k))/n
                      END IF
                    END DO
                  END IF
                END IF
                

              ELSE IF( VarType == Variable_on_elements ) THEN
                
                m = CurrentElement % ElementIndex
                
                IF( ASSOCIATED( Perm ) ) m = Perm( m ) 
                
                IF( sdofs == 1 ) THEN
                  ElemVectVal(1) = Values(m) 
                ELSE
                  DO k=1,sdofs
                    IF( k > dofs ) THEN
                      ElemVectVal(k) = 0.0_dp
                    ELSE IF(ComponentVector) THEN
                      IF (k==1) ElemVectVal(k) = Values(m)
                      IF (k==2) ElemVectVal(k) = Values2(m)
                      IF (k==3) ElemVectVal(k) = Values3(m)
                    ELSE
                      ElemVectVal(k) = Values(dofs*(m-1)+k)
                    END IF
                  END DO
                END IF
                
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
      END DO
      IF( WriteXML ) THEN
        CALL Info('VtuOutputSolver','Number of elemental fields written: '//TRIM(I2S(NoFieldsWritten)),Level=10)
      END IF
    END IF  ! IF( SaveElemental )

    !---------------------------------------------------------------------
    ! If requested write the body and bc indexes
    !---------------------------------------------------------------------
    IF( WriteIds ) THEN
      IF( WriteXML ) THEN
        CALL Info('VtuOutputSolver','Writing entity IDs for bodies and boundaries',Level=10)

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
            j = GeometryBodyMap( j )
          ELSE
            j = GetBCId( CurrentElement ) 
            IF ( j>=1 .AND. j<= SIZE(GeometryBCMap)) j = GeometryBCMap( j )
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
      CALL Info('VtuOutputSolver','Writing coordinates for each used node',Level=10)
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

    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k 

      DO ii = 1, NumberOfDofNodes
        IF( NoPermutation ) THEN
          i = ii 
        ELSE
          i = InvNodePerm(ii)
        END IF

        x = Model % Mesh % Nodes % x( i )
        y = Model % Mesh % Nodes % y( i )
        z = Model % Mesh % Nodes % z( i )


        ! If displacement field is active remove the displacement from the coordinates
        IF( dispdofs > 0 .OR. disp2dofs > 0) THEN
          j = 0
          IF(dispdofs > 0) THEN
            j = DispPerm(i)
            IF( j > 0 ) THEN
              x = x - DispValues(dispdofs*(j-1)+1)
              y = y - DispValues(dispdofs*(j-1)+2)
              IF(dispdofs == 3) z = z - DispValues(dispdofs*(j-1)+3)
            END IF
          END IF
          IF(disp2dofs > 0 .AND. j==0) THEN
            j = Disp2Perm(i)
            IF( j > 0 ) THEN
              x = x - Disp2Values(disp2dofs*(j-1)+1)
              y = y - Disp2Values(disp2dofs*(j-1)+2)
              IF(disp2dofs == 3) z = z - Disp2Values(disp2dofs*(j-1)+3)
            END IF
          END IF
        END IF

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
      CALL Info('VtuOutputSolver','Writing the elemental connectivity data',Level=10)
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
        IF( SaveLinear ) THEN
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

        !          NodeIndexes => Elmer2VtkIndexes( CurrentElement, DG .OR. DN, SaveLinear )
        CALL Elmer2VtkIndexes( CurrentElement, DG .OR. DN, SaveLinear, TmpIndexes )

        IF( SaveLinear ) THEN
          n = GetElementCorners( CurrentElement ) 
        ELSE
          n = GetElementNOFNodes( CurrentElement )
        END IF

        DO j=1,n
          IF( DN .OR. DG ) THEN
            jj = DgPerm( TmpIndexes(j) )
!          ELSE IF( DG ) THEN
!            jj = TmpIndexes(j)
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
        IF( SaveLinear ) THEN
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
        n = Elmer2VTKElement(CurrentElement % TYPE % ElementCode, SaveLinear )

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
    DEALLOCATE(ElemInd)

  END SUBROUTINE WriteVtuFile



  SUBROUTINE WritePvdFile( PvdFile, DataSetFile, nTime, Model )
    CHARACTER(LEN=*), INTENT(IN) :: PvdFile, DataSetFile
    INTEGER :: nTime, RecLen = 0
    TYPE(Model_t) :: Model     
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: n
    REAL(KIND=dp) :: time
    CHARACTER :: lf
    CHARACTER(LEN=MAX_NAME_LEN) :: Str
    LOGICAL :: Found

    SAVE RecLen

    lf = CHAR(10)

    IF( ParEnv % PEs > 1 ) THEN
      IF( ParEnv % MyPE > 0 ) RETURN
    END IF
    time = GetTime()
    IF( GetLogical( Params,'Vtu time previous',Found) ) THEN
      time = time - GetTimestepSize()
    END IF


    IF( nTime == 1 .OR. Reclen == 0 ) THEN
      ! Find the maximum record length (modulo four)
      WRITE( Str,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'
      n = LEN_TRIM( Str ) 

      WRITE( Str,'(A,ES16.7,A)') '<DataSet timestep="',time,&
        '" group="" part="0" file="'//TRIM(DataSetFile)//'"/>'
      n = MAX( LEN_TRIM( Str ), n ) 
      
      RecLen = ((n/4)+1)*4
    END IF

    IF( nTime == 1 ) THEN
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='REPLACE', &
          ACCESS='DIRECT', ACTION='WRITE', RECL=RecLen)

      IF ( LittleEndian() ) THEN
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'//lf
      ELSE
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="BigEndian"><Collection>'//lf
      END IF     
    ELSE
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='OLD', &
          ACCESS='DIRECT', ACTION='READWRITE', RECL=RecLen)     
    END IF

    WRITE( VtuUnit,'(A,ES12.3,A)',REC=nTime+1) '<DataSet timestep="',time,&
        '" group="" part="0" file="'//TRIM(DataSetFile)//'"/>'//lf
    WRITE( VtuUnit,'(A)',REC=nTime+2) '</Collection></VTKFile>'//lf

    CLOSE( VtuUnit )

  END SUBROUTINE WritePvdFile


  SUBROUTINE WritePvtuFile( VtuFile, Model )
    CHARACTER(LEN=*), INTENT(IN) :: VtuFile
    TYPE(Model_t) :: Model 
    INTEGER, PARAMETER :: VtuUnit = 58
    TYPE(Variable_t), POINTER :: Var,Var1
    CHARACTER(LEN=512) :: str
    INTEGER :: i,j,k,dofs,Rank,cumn,n,dim,vari,sdofs
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FullName, ShortName
    LOGICAL :: ScalarsExist, VectorsExist, Found, VeloFlag, ComponentVector, &
               AllActive, ThisActive, L
    LOGICAL, POINTER :: ActivePartition(:)
    TYPE(Variable_t), POINTER :: Solution
    INTEGER, POINTER :: Perm(:)
    INTEGER :: Active, NoActive, ierr, NoFields, NoModes, IndField, iField, VarType
    REAL(KIND=dp), POINTER :: Values(:)
    COMPLEX(KIND=dp), POINTER :: EigenVectors(:,:)
    TYPE(Element_t), POINTER :: CurrentElement

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    

    ThisActive = ( NumberOfGeomNodes > 0 ) 
   
    Active = 0
    IF( ThisActive ) Active = 1    	

    NoActive = NINT( ParallelReduction( 1.0_dp * Active ) )
    IF( NoActive == 0 ) THEN
      CALL Warn('WritePvtuFile','No active partitions for saving')
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


    OPEN( UNIT=VtuUnit, FILE=VtuFile, form = 'formatted', STATUS='UNKNOWN' )
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
        
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
        ComponentVector = .FALSE.

        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE
            WRITE(Txt, '(A,A)') 'Nonexistent variable 2: ',TRIM(FieldName)
            CALL Warn('WriteVtuXMLFile', Txt)
            CYCLE
          END IF
        END IF
        
        VarType = Solution % Type
        IF( VarType == Variable_on_nodes_on_elements ) THEN
          IF( .NOT. ( ( DG .OR. DN ) .AND. SaveElemental ) ) CYCLE
        END IF

        IF( ASSOCIATED(Solution % EigenVectors)) THEN
           NoModes = SIZE( Solution % EigenValues )
           IF( ComponentVector ) THEN
             CALL Warn('WritePvtuXMLFile','Eigenmodes cannot be given componentwise!')
             CYCLE
           ELSE IF( EigenAnalysis ) THEN
             IF( GotActiveModes ) THEN
               IndField = ActiveModes( FileIndex ) 
             ELSE
               IndField = FileIndex
             END IF
             IF( IndField > NoModes ) THEN
               CALL Warn('WriteVtuXMLFile','Too few eigenmodes!')
               CYCLE
             END IF
             NoModes = 1
             NoFields = 1
           ELSE	  
             IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
             NoFields = NoModes
           END IF
        ELSE
          NoModes = 0 
          NoFields = 1
        END IF
  
        dofs = Solution % DOFs
        IF( ComponentVector ) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
          IF( ASSOCIATED(Solution)) dofs = 2
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
          IF( ASSOCIATED(Solution)) dofs = 3
        END IF
        
        IF( dofs > 1 ) THEN
          sdofs = MAX(dofs,3)
        ELSE
          sdofs = 1
        END IF

        DO iField = 1, NoFields

          IF( NoModes == 0 .OR. EigenAnalysis ) THEN
            FullName = TRIM( FieldName ) 
          ELSE          
            IF( GotActiveModes ) THEN
              IndField = ActiveModes( iField ) 
            ELSE
              IndField = iField
            END IF
            WRITE( FullName,'(A,I0)') TRIM( FieldName )//' mode',IndField
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
    WRITE( VtuUnit,'(A)') '    </PPointData>'
  END IF


    ! Elementwise information
    !-------------------------------------
    WRITE( VtuUnit,'(A)') '    <PCellData>'

  IF( SaveElemental  .AND. .NOT. ( DG .OR. DN ) ) THEN
    IF( ScalarsExist .OR. VectorsExist ) THEN
      DO Rank = 0,2
        DO Vari = 1, 999

          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field Elemental ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field Elemental ',Vari
          FieldName = GetString( Params, TRIM(Txt), Found )
          L = Found 
          
          IF(.NOT. Found) THEN          
            IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
            FieldName = GetString( Params, TRIM(Txt), Found )
          END IF
          IF(.NOT. Found ) EXIT

          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
          ComponentVector = .FALSE.

          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE 
              IF( L ) THEN
                WRITE(Txt, '(A,A)') 'Nonexistent elemental variable 2: ',TRIM(FieldName)
                CALL Warn('WriteVtuXMLFile', Txt)
              END IF
              CYCLE
            END IF
          END IF
          
          VarType = Solution % Type
          IF( VarType /= Variable_on_nodes_on_elements ) CYCLE

          IF( ASSOCIATED(Solution % EigenVectors)) THEN
            NoModes = SIZE( Solution % EigenValues )
            IF( ComponentVector ) THEN
              CALL Warn('WritePvtuXMLFile','Eigenmodes cannot be given componentwise!')
              CYCLE
            ELSE IF( EigenAnalysis ) THEN
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( FileIndex ) 
              ELSE
                IndField = FileIndex
              END IF
              IF( IndField > NoModes ) THEN
                CALL Warn('WriteVtuXMLFile','Too few eigenmodes!')
                CYCLE
              END IF
              NoModes = 1
              NoFields = 1
            ELSE	  
              IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
              NoFields = NoModes
            END IF
          ELSE
            NoModes = 0 
            NoFields = 1
          END IF
          
          dofs = Solution % DOFs
          IF( ComponentVector ) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
            IF( ASSOCIATED(Solution)) dofs = 2
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
            IF( ASSOCIATED(Solution)) dofs = 3
          END IF
          
          IF( dofs > 1 ) THEN
            sdofs = MAX(dofs,3)
          ELSE
            sdofs = 1
          END IF


          DO iField = 1, NoFields
            
            IF( NoModes == 0 .OR. EigenAnalysis ) THEN
              FullName = TRIM( FieldName ) 
            ELSE          
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( iField ) 
              ELSE
                IndField = iField
              END IF
              WRITE( FullName,'(A,I0)') TRIM( FieldName )//' mode',IndField
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
    j = INDEX( FilePrefix,'/') 
    IF( j == 0 ) THEN
      ShortName = FilePrefix
    ELSE
      ShortName = FilePrefix(j+1:)
    END IF

    DO i=1,Partitions
      IF(.NOT. AllActive ) THEN
        IF( .NOT. ActivePartition(i)) CYCLE
      END IF

      IF( NoFileindex ) THEN
        WRITE( VtuUnit,'(A,I4.4,A,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",'.vtu"/>'
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuUnit,'(A,I4.4,A,I4.4,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",FileIndex,'.vtu"/>'        
      ELSE
        WRITE( VtuUnit,'(A,I4.4,A,I0,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",FileIndex,'.vtu"/>'        
      END IF
    END DO

    WRITE( VtuUnit,'(A)') '  </PUnstructuredGrid>'
    WRITE( VtuUnit,'(A)') '</VTKFile>'

    CLOSE( VtuUnit )

    IF(.NOT. AllActive ) DEALLOCATE( ActivePartition ) 
  
  END SUBROUTINE WritePvtuFile

!------------------------------------------------------------------------------
END SUBROUTINE VtuOutputSolver
!------------------------------------------------------------------------------
  
