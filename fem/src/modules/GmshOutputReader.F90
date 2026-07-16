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
!> Reads mesh and results in ascii format understood by the pre-/postprocessing software Gmsh.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE GmshOutputReader( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------  
  USE DefUtils
  USE SaveUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp) :: Time
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(ValueList_t), POINTER :: SolverParams
  
  LOGICAL :: Found, UseBBox, UseQuadTree, AllocationsDone
  INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerType, GmshType,body_id,&
      Vari, Rank, NoNodes, NoElems, NoBulkElems, ElemDim, MaxElemDim, &
      MaxElemNodes, InputPartitions, ReadPart, AlignCoord, PassiveCoord, &
      CumNodes, CumElems, iostat, NoVars, MeshDim
  INTEGER, POINTER :: NodeIndexes(:), ElmerIndexes(:), MaskPerm(:)
  REAL(KIND=dp) :: x,y,z,dx,BBTol
  REAL(KIND=dp), ALLOCATABLE :: BBox(:,:), MyBBox(:)
  LOGICAL, ALLOCATABLE :: ActivePart(:)
  INTEGER, PARAMETER :: LENGTH = 1024
  CHARACTER(LEN=LENGTH) :: Txt, FieldName, CompName, str
  CHARACTER(MAX_NAME_LEN) :: BaseFile, InputFile, VarName
  CHARACTER(:), ALLOCATABLE :: InputDirectory
  TYPE(Mesh_t), POINTER :: FromMesh, ToMesh
  TYPE(Variable_t), POINTER :: Var
  INTEGER :: FileUnit=1
  TYPE VariableArray_t
    TYPE(Variable_t), POINTER :: Var => NULL()
  END TYPE VariableArray_t
  TYPE(VariableArray_t) :: VarArray(20)

  CHARACTER(*), PARAMETER :: Caller = 'GmshOutputReader'  
  
  INTERFACE
    SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, NewVariables, &
        UseQuadrantTree, Projector, MaskName, FoundNodes, NewMaskPerm, KeepUnfoundNodes )
      USE Types
      TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
      TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
      LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
      CHARACTER(LEN=*),OPTIONAL :: MaskName
      TYPE(Projector_t), POINTER, OPTIONAL :: Projector
      INTEGER, OPTIONAL, POINTER :: NewMaskPerm(:)  
      LOGICAL, OPTIONAL :: KeepUnfoundNodes  
    END SUBROUTINE InterpolateMeshToMeshQ
  END INTERFACE

       
  SAVE AllocationsDone
  
!------------------------------------------------------------------------------

  CALL Info(Caller,'Reading Gmsh results and interpolating to current mesh!')


  SolverParams => GetSolverParams()
  ToMesh => Solver % Mesh

  dim = CoordinateSystemDimension()
  AlignCoord = ListGetInteger( SolverParams,'Align Coordinate',Found )

  UseQuadTree = ListGetLogical( SolverParams,'Use Quadrant Tree',Found ) 
  
  PassiveCoord = ListGetInteger( SolverParams,'Interpolation Passive Coordinate',Found)
  IF(.NOT. Found) PassiveCoord = ABS(AlignCoord)
  
  BaseFile = ListGetString( SolverParams, 'Filename', UnfoundFatal = .TRUE. ) 
  IF(INDEX(BaseFile,'.') == 0) WRITE( BaseFile,'(A,A)') TRIM(BaseFile),".msh"    
  CALL SolverOutputDirectory( Solver, BaseFile, InputDirectory, UseMeshDir = .TRUE. )
  BaseFile = TRIM(InputDirectory)// '/' //TRIM(BaseFile)

  InputPartitions = ListGetInteger( SolverParams,'Filename Partitions', Found )
  IF(.NOT. Found) InputPartitions = 1

  n = -1
  MaskPerm => NULL()
  Str = ListGetString( SolverParams,'Mask Name',Found) 
  IF( Found ) THEN
    ALLOCATE( MaskPerm( ToMesh % NumberOfNodes ) )
    MaskPerm = 0
    CALL MakePermUsingMask( Model, Solver, ToMesh, Str, .FALSE., &
        MaskPerm, n, RequireLogical = .TRUE. )
    CALL Info(Caller,'Using given mask "'//TRIM(Str)//'" for interpolation!')            
    IF(n==0) DEALLOCATE(MaskPerm)
  ELSE IF( ASSOCIATED( Solver % Variable ) ) THEN
    MaskPerm => Solver % Variable % Perm
    NULLIFY(MaskPerm)
    IF( ASSOCIATED(MaskPerm) ) THEN
      CALL Info(Caller,'Using Solver % Variable % Perm as the mask for interpolation!')        
      n = COUNT(MaskPerm(1:ToMesh % NumberOfNodes) > 0)
    END IF
  END IF

  IF( n == 0 ) THEN
    CALL Info(Caller,'Zero masked nodes, returning without interpolation!')
    RETURN
  ELSE IF( n > 0 ) THEN
    CALL Info(Caller,'Number of masked nodes: '//I2S(n),Level=7)
  END IF
    
  UseBBox = .FALSE.
  IF( InputPartitions > 1 ) THEN
    UseBBox = ListGetLogical( SolverParams,'Use Bounding Box',Found ) 
  END IF
   
  ! We use simple bounding box to avoid reading unnecessary pieces in parallel.
  ! By default all pieces are read and a union mesh is created on-the-fly. 
  IF( UseBBox ) THEN
    ALLOCATE(Bbox(InputPartitions,6),MyBBox(6))
    Bbox(:,1:3) = HUGE(x)   ! initialize min values
    Bbox(:,4:6) = -HUGE(x)  ! initialize max values
    ALLOCATE(ActivePart(InputPartitions))
    ActivePart = .TRUE.

    n = ToMesh % NumberOfNodes 
    IF( ASSOCIATED( MaskPerm ) ) THEN
      MyBbox(1) = MINVAL(ToMesh % Nodes % x(1:n),MaskPerm(1:n)>0)
      MyBbox(2) = MINVAL(ToMesh % Nodes % y(1:n),MaskPerm(1:n)>0)
      MyBbox(3) = MINVAL(ToMesh % Nodes % z(1:n),MaskPerm(1:n)>0)
      MyBbox(4) = MAXVAL(ToMesh % Nodes % x(1:n),MaskPerm(1:n)>0)
      MyBbox(5) = MAXVAL(ToMesh % Nodes % y(1:n),MaskPerm(1:n)>0)
      MyBbox(6) = MAXVAL(ToMesh % Nodes % z(1:n),MaskPerm(1:n)>0)
    ELSE
      MyBbox(1) = MINVAL(ToMesh % Nodes % x(1:n))
      MyBbox(2) = MINVAL(ToMesh % Nodes % y(1:n))
      MyBbox(3) = MINVAL(ToMesh % Nodes % z(1:n))
      MyBbox(4) = MAXVAL(ToMesh % Nodes % x(1:n))
      MyBbox(5) = MAXVAL(ToMesh % Nodes % y(1:n))
      MyBbox(6) = MAXVAL(ToMesh % Nodes % z(1:n))
    END IF
      
    BBtol = ListGetConstReal( SolverParams,'Bounding box tolerance',Found )
    IF(.NOT. Found ) BBtol = 1.0d-6
  END IF
  
  AllocationsDone = .FALSE.
  MaxElemDim = 0
  MaxElemNodes = 0
  
10 CONTINUE

  CumNodes = 0
  CumElems = 0
  
  DO ReadPart = 1, InputPartitions 
    
    NoNodes = 0
    NoElems = 0

    IF( InputPartitions == 1 ) THEN
      InputFile = TRIM(BaseFile)
    ELSE
      IF( AllocationsDone .AND. UseBBox ) THEN
        IF( .NOT. ActivePart(ReadPart) ) CYCLE
      END IF

      ! With this heuristics ensure that not every partition tries to read the same file
      ! In ideal case everybody starts on their own partition.
      i = MODULO(ReadPart-1 + ParEnv % MyPe, InputPartitions)
      InputFile = TRIM(BaseFile)//'_'//I2S(InputPartitions)//'np'//I2S(i+1)
    END IF

    IF( InputPartitions > 1 .OR. .NOT. AllocationsDone ) THEN
      CALL Info(Caller,'Reading mesh from file: '//TRIM(InputFile),Level=7)
      OPEN(UNIT=FileUnit, FILE=InputFile, STATUS='old',iostat=iostat)
      IF( iostat /= 0 ) THEN
        IF( InputPartitions > 1 ) THEN
          IF( UseBBox ) ActivePart(ReadPart) = .FALSE.
          ! We need not have all partition present
          CALL Info(Caller,'Could not open file: '//TRIM(InputFile))
          CYCLE
        ELSE
          CALL Fatal(Caller,'Could not open file: '//TRIM(InputFile))
        END IF
      END IF
    END IF
    
    CALL ReadSingleGmshFile()

    CumNodes = CumNodes + NoNodes
    CumElems = CumElems + NoElems 
    
    IF( InputPartitions > 1 .OR. AllocationsDone ) THEN
      CLOSE( FileUnit )
    ELSE
      REWIND( FileUnit )
    END IF    
  END DO

  IF(.NOT. AllocationsDone ) THEN
    IF( InputPartitions > 1 ) THEN
      CALL Info(Caller,'Number of cumulative nodes to read: '//I2S(CumNodes),Level=5)
      CALL Info(Caller,'Number of cumulative elements to read: '//I2S(CumElems),Level=5)
    END IF
    IF( UseBBox ) THEN
      i = COUNT( ActivePart )
      CALL Info(Caller,'Number of active partitions '&
          //I2S(i)//' out of '//I2S(InputPartitions),Level=10)
    END IF

    CALL Info(Caller,'Maximum element dimension: '//I2S(MaxElemDim),Level=7)
    CALL Info(Caller,'Maximum element nodes: '//I2S(MaxElemNodes),Level=7)    
    FromMesh => AllocateMesh(CumElems,0,CumNodes)

    ! This is temporal only
    FromMesh % MeshDim = 3
    FromMesh % MaxElementNodes = MaxElemNodes
    AllocationsDone = .TRUE.    

    CALL Info(Caller,'Creating variable structure',Level=20)
    FromMesh % Variables => NULL()
    ALLOCATE( Perm( CumNodes ) )
    DO i =1, CumNodes
      Perm(i) = i
    END DO

    GOTO 10 
  END IF

  CALL Info(Caller,'Last bulk element index: '//I2S(NoBulkElems),Level=7)

  FromMesh % MeshDim = MAX( MaxElemDim, MeshDim ) 
  FromMesh % NumberOfBulkElements = NoBulkElems
  FromMesh % NumberOfBoundaryElements = NoElems - NoBulkElems
    
  CALL Info(Caller,'Gmsh reader complete!')

  CALL InterpolateFromGmshFile()
  
  CALL Info(Caller,'Interpolation from Gmsh format complete')
  

CONTAINS

  ! Read a Gmsh file with results. When reading multiple files use the offsets for
  ! CumNodes and CumElems. 
  SUBROUTINE ReadSingleGmshFile()

    REAL(KIND=dp) :: GmshVer
    REAL(KIND=dp) :: coord(3)
    INTEGER :: GmshToElmerType(21), GmshIndexes(27)
    INTEGER :: i,j,k
    
    SAVE GmshVer
    
    GmshToElmerType = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
        510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)

    NoVars = 0
    MeshDim = 0
    
    DO WHILE( .TRUE. )
      READ( FileUnit,'(A)',END=20,ERR=20 ) str    
      IF ( SEQL( str, '$MeshFormat') ) THEN
        READ( FileUnit,'(A)',END=20,ERR=20 ) str    
        IF(.NOT. AllocationsDone ) THEN
          READ( str,*) GmshVer
          WRITE(Message,'(A,ES12.3)') 'Gmsh file version: ',GmshVer
          CALL Info(Caller,Message)
        END IF
      END IF

      IF ( SEQL( str, '$Nodes') ) THEN
        READ( FileUnit,'(A)',END=20,ERR=20 ) str    
        READ( str,*) NoNodes
        IF(.NOT. AllocationsDone ) THEN
          CALL Info(Caller,'Number of nodes in mesh: '//I2S(NoNodes),Level=7)
        END IF
        DO i=1,NoNodes 
          READ( FileUnit,'(A)',END=20,ERR=20 ) str     
          READ( str,*) j, coord
          IF( i /= j ) CALL Fatal(Caller,'Do node permutations!')
          IF( AllocationsDone ) THEN
            ! In parallel we have an offset, in serial it is zero
            k = i + CumNodes
            FromMesh % Nodes % x(k) = coord(1)
            FromMesh % Nodes % y(k) = coord(2)
            FromMesh % Nodes % z(k) = coord(3)
          ELSE IF( UseBBox ) THEN
            Bbox(ReadPart,1:3) = MIN(Bbox(ReadPart,1:3),Coord)
            Bbox(ReadPart,4:6) = MAX(Bbox(ReadPart,4:6),Coord)            
          END IF
        END DO
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! EndNodes  
        
        IF(.NOT. AllocationsDone ) THEN
          IF( UseBBox ) THEN
            ! Check the bounding boxes so that we do not need to go through this
            ! piece again if the bounding boxes do not overlap.
            DO i=1,3
              IF(i==PassiveCoord) CYCLE
              IF(i>dim) CYCLE
              IF( Bbox(ReadPart,i) > MyBbox(i+3) + BBtol ) ActivePart(ReadPart) = .FALSE.
              IF( Bbox(ReadPart,i+3) < MyBbox(i) - BBtol ) ActivePart(ReadPart) = .FALSE.
            END DO
            IF( .NOT. ActivePart(ReadPart) ) THEN
              NoNodes = 0
              NoElems = 0
              EXIT
            END IF
          END IF
      
          IF( NoElems > 0 .AND. NoNodes > 0 ) EXIT
        END IF
      END IF ! Nodes

                       
      IF ( SEQL( str, '$Elements') ) THEN
        READ( FileUnit,'(A)',END=20,ERR=20 ) str    
        READ( str,*) NoElems
        IF(.NOT. AllocationsDone ) THEN
          CALL Info(Caller,'Number of elements in mesh: '//I2S(NoElems),Level=7)
        END IF

        DO i=1,NoElems
          READ( FileUnit,'(A)',END=20,ERR=20 ) str     
          READ(str,*) j,GmshType
          IF( i /= j ) CALL Fatal(Caller,'Do element permutations!')        
          ElmerType = GmshToElmerType(GmshType)

          ElemDim = 0
          IF( ElmerType > 500 ) THEN
            ElemDim = 3
          ELSE IF( ElmerType > 300 ) THEN
            ElemDim = 2
          ELSE IF( ElmerType > 200 ) THEN
            ElemDim = 1
          END IF
                    
          IF( AllocationsDone ) THEN
            k = CumElems + i
            Element => FromMesh % Elements(k)
            Element % TYPE => GetElementType(ElmerType)
            n = Element % TYPE % NumberOfNodes
            CALL AllocateVector( Element % NodeIndexes, n )          
            READ(str,*) j,GmshType,j,j,j, Element % Nodeindexes(1:n)
            IF( CumNodes > 0 ) THEN
              Element % NodeIndexes(1:n) = Element % NodeIndexes(1:n) + CumNodes
            END IF
            IF( ElemDim == MaxElemDim ) NoBulkElems = k
          ELSE
            MaxElemDim = MAX( MaxElemDim, ElemDim ) 
            MaxElemNodes = MAX( MaxElemNodes, MODULO( ElmerType, 100 ) )        
          END IF
        END DO
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! EndElements  

        ! Leave early since we don't have structures to read the results
        IF( .NOT. AllocationsDone ) THEN
          IF( NoElems > 0 .AND. NoNodes > 0 ) EXIT
        END IF
      END IF ! Elements


      IF ( SEQL( str, '$NodeData') ) THEN
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! 1  
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! "name"            
        n = LEN_TRIM(str)

        VarName = str(2:n-1)
        CALL Info(Caller,'Reading gmsh variable: '//TRIM(VarName))

        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! 1  
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! time
        READ( str, * ) time
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! dim
        READ( str, * ) MeshDim 
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! visited
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! dofs
        READ( str, * ) dofs
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! n
        READ( str, * ) n

        IF( n /= NoNodes ) THEN
          CALL Info(Caller,'Mismatch is vector size!')
        END IF

        Var => NULL()
        IF( ASSOCIATED( FromMesh % Variables ) ) THEN
          CALL Info(Caller,'Reusing variable structure',Level=20)
          Var => VariableGet( FromMesh % Variables, VarName, ThisOnly = .TRUE. )
        END IF

        IF(.NOT. ASSOCIATED(Var) ) THEN
          CALL VariableAddVector( FromMesh % Variables, FromMesh, Solver, &
              VarName, dofs = dofs, Perm = Perm )
          Var => VariableGet( FromMesh % Variables, VarName, ThisOnly = .TRUE. )
        END IF
          
        DO i=1,NoNodes
          READ( FileUnit,'(A)',END=20,ERR=20 ) str  
          k = i + CumNodes
          !IF( k > SIZE( Var % Values ) ) CALL Fatal('','k is too large!')
          IF( dofs == 1 ) THEN
            READ( str, * ) j, Var % Values(k)
          ELSE
            READ( str, * ) j, Var % Values(dofs*(k-1)+1:dofs*k) 
          END IF
        END DO
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! EndNodeData      

        NoVars = NoVars + 1
        VarArray(NoVars) % Var => Var
                
        IF( InfoActive(28) ) THEN          
          IF(NoVars==1) CALL Info(Caller,'Initial field ranges:')
          CALL VectorValuesRange(Var % Values,SIZE(Var % Values),'From: '//TRIM(Var % Name))       
        END IF
      END IF ! NodeData
    END DO

    
20  CONTINUE

  END SUBROUTINE ReadSingleGmshFile


  
  ! Interpolate from the loaded Gmsh file to exiting mesh.
  !------------------------------------------------------
  SUBROUTINE InterpolateFromGmshFile()

    REAL(KIND=dp), POINTER :: x1(:), x2(:)
    REAL(KIND=dp) :: minx, maxx
    INTEGER :: n1,n2
    
    IF( AlignCoord /= 0 ) THEN
      k = ABS( AlignCoord ) 

      IF( k == 1 ) THEN
        x1 => FromMesh % Nodes % x
        x2 => ToMesh % Nodes % x
      ELSE IF( k == 2 ) THEN
        x1 => FromMesh % Nodes % y
        x2 => ToMesh % Nodes % y
      ELSE IF( k == 3 ) THEN      
        x1 => FromMesh % Nodes % z
        x2 => ToMesh % Nodes % z
      ELSE
        CALL Fatal(Caller,'Invalid value for "Align Coordinate": '//I2S(AlignCoord))
      END IF

      n1 = FromMesh % NumberOfNodes
      n2 = ToMesh % NumberOfNodes

      IF( AlignCoord > 0 ) THEN
        minx = MINVAL( x2(1:n2) )
        maxx = MAXVAL( x1(1:n1) ) 
        dx = minx - maxx
      ELSE
        minx = MINVAL( x1(1:n1) )
        maxx = MAXVAL( x2(1:n2) ) 
        dx = minx - maxx
      END IF

      WRITE(Message,'(A,ES12.3)') 'Aligning coordinate '//I2S(k)//' with: ',dx
      CALL Info(Caller,Message)

      x1 = x1 + dx
    END IF


    !CALL InspectMesh(FromMesh)
    !CALL InspectMesh(ToMesh)
    
    IF( ASSOCIATED( MaskPerm ) ) THEN
      CALL InterpolateMeshToMeshQ( FromMesh, ToMesh, FromMesh % Variables, ToMesh % Variables, &
          UseQuadrantTree=UseQuadTree,NewMaskPerm = MaskPerm ) 
    ELSE
      CALL InterpolateMeshToMeshQ( FromMesh, ToMesh, FromMesh % Variables, ToMesh % Variables, &
          UseQuadrantTree=UseQuadTree)
    END IF
    
    IF( InfoActive(28) ) THEN
      CALL Info(Caller,'Projected field ranges:')
      DO i=1,NoVars        
        Var => VariableGet( ToMesh % Variables, VarArray(i) % Var % Name, ThisOnly = .TRUE.)
        IF(ASSOCIATED(Var)) THEN
          CALL VectorValuesRange(Var % Values,SIZE(Var % Values),'To: '//TRIM(Var % Name))          
        END IF
      END DO
    END IF
        
    
  END SUBROUTINE InterpolateFromGmshFile
    
  
!------------------------------------------------------------------------------
END SUBROUTINE GmshOutputReader
!------------------------------------------------------------------------------


SUBROUTINE GmshOutputReader_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()
  CALL ListAddNewLogical(Params,'No Matrix',.TRUE.)  

END SUBROUTINE GmshOutputReader_init
