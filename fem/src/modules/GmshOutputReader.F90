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
  INTEGER, POINTER :: NewMaskPerm(:)
  
  LOGICAL :: Found, AllocationsDone = .FALSE.
  
  INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerType, GmshType,body_id,&
      Vari, Rank, NoNodes, NoElems, NoBulkElems, ElemDim, MaxElemDim, &
      MaxElemNodes, AlignCoord
  INTEGER :: GmshToElmerType(21), GmshIndexes(27) 
  INTEGER, POINTER :: NodeIndexes(:), ElmerIndexes(:), MaskPerm(:)
  REAL(KIND=dp) :: x,y,z,GmshVer,dx
  REAL(KIND=dp), POINTER :: x1(:), x2(:)
  INTEGER, PARAMETER :: LENGTH = 1024
  CHARACTER(LEN=LENGTH) :: Txt, FieldName, CompName, str
  CHARACTER(MAX_NAME_LEN) :: InputFile, InputDirectory, VarName
  INTEGER :: FileUnit=1
  TYPE(Mesh_t), POINTER :: FromMesh, ToMesh
  TYPE(Variable_t), POINTER :: Var
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


  GmshToElmerType = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
      510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)


  SolverParams => GetSolverParams()

  InputFile = ListGetString( SolverParams, 'Filename', UnfoundFatal = .TRUE. ) 
  IF(INDEX(InputFile,'.') == 0) WRITE( InputFile,'(A,A)') TRIM(InputFile),".msh"
    
  CALL SolverOutputDirectory( Solver, InputFile, InputDirectory, UseMeshDir = .TRUE. )
  InputFile = TRIM(InputDirectory)// '/' //TRIM(InputFile)
   
  dim = CoordinateSystemDimension()

  CALL Info(Caller,'Reading mesh from file: '//TRIM(InputFile))
  OPEN(UNIT=FileUnit, FILE=InputFile, STATUS='old')


  NoNodes = 0
  NoElems = 0
  MaxElemDim = 0
  MaxElemNodes = 0
  
10 CONTINUE
  
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
      IF(.NOT. AllocationsDone ) THEN
        READ( str,*) NoNodes
        CALL Info(Caller,'Number of nodes in mesh: '//TRIM(I2S(NoNodes)),Level=7)
      END IF
      DO i=1,NoNodes 
        READ( FileUnit,'(A)',END=20,ERR=20 ) str     
        READ( str,*) j,x,y,z
        IF( i /= j ) CALL Fatal(Caller,'Do node permutations!')
        IF( AllocationsDone ) THEN
          FromMesh % Nodes % x(i) = x
          FromMesh % Nodes % y(i) = y
          FromMesh % Nodes % z(i) = z
        END IF
      END DO
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! EndNodes  
      IF(.NOT. AllocationsDone ) THEN
        IF( NoElems > 0 .AND. NoNodes > 0 ) EXIT
      END IF
    END IF ! Nodes


    IF ( SEQL( str, '$Elements') ) THEN
      READ( FileUnit,'(A)',END=20,ERR=20 ) str    
      IF(.NOT. AllocationsDone ) THEN
        READ( str,*) NoElems
        CALL Info(Caller,'Number of elements in mesh: '//TRIM(I2S(NoElems)),Level=7)
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
          Element => FromMesh % Elements(i)
          Element % TYPE => GetElementType(ElmerType)
          n = Element % TYPE % NumberOfNodes
          CALL AllocateVector( Element % NodeIndexes, n )          
          READ(str,*) j,GmshType,k,k,k, Element % Nodeindexes(1:n)
          IF( ElemDim == MaxElemDim ) NoBulkElems = i
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
      CALL Info(Caller,'Reading variable: '//TRIM(VarName))
            
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! 1  
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! time
      READ( str, * ) time
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! dim  
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! visited
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! dofs
      READ( str, * ) dofs
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! n
      READ( str, * ) n

      IF( n /= NoNodes ) THEN
        CALL Info(Caller,'Mismatch is vector size!')
      END IF
            
      IF( ASSOCIATED( FromMesh % Variables ) ) THEN
        CALL Info(Caller,'Creating variable structure',Level=20)
      ELSE
        CALL Info(Caller,'Reusing variable structure',Level=20)
        FromMesh % Variables => Null()
        ALLOCATE( Perm( NoNodes ) )
        DO i =1, NoNodes
          Perm(i) = i
        END DO
      END IF
      
      CALL VariableAddVector( FromMesh % Variables, FromMesh, Solver, &
          VarName, dofs = dofs, Perm = Perm )
      Var => VariableGet( FromMesh % Variables, VarName )
      
      DO i=1,NoNodes
        READ( FileUnit,'(A)',END=20,ERR=20 ) str  
        IF( dofs == 1 ) THEN
          READ( str, * ) j, Var % Values(i)
        ELSE
          READ( str, * ) j, Var % Values(dofs*(i-1)+1:dofs*i) 
        END IF
      END DO
      READ( FileUnit,'(A)',END=20,ERR=20 ) str  ! EndNodeData      
      !PRINT *,'Var Range:',MINVAL( Var % Values ), MAXVAL( Var % Values ) 
    END IF ! NodeData
  END DO

20 CONTINUE
  
  IF(AllocationsDone ) THEN
    CALL Info(Caller,'Last bulk element index: '//TRIM(I2S(NoBulkElems)),Level=7)
    FromMesh % NumberOfBulkElements = NoBulkElems
    FromMesh % NumberOfBoundaryElements = NoElems - NoBulkElems
  ELSE
    CALL Info(Caller,'Maximum element dimension: '//TRIM(I2S(MaxElemDim)),Level=7)
    CALL Info(Caller,'Maximum element nodes: '//TRIM(I2S(MaxElemNodes)),Level=7)
    
    FromMesh => AllocateMesh(NoElems,0,NoNodes)    
    FromMesh % MeshDim = MaxElemDim
    FromMesh % MaxElementNodes = MaxElemNodes
       
    REWIND( FileUnit ) 
    AllocationsDone = .TRUE.    
    GOTO 10 
  END IF
    
  CLOSE(FileUnit)
  
  CALL Info(Caller,'Gmsh reader complete!')


  ToMesh => Solver % Mesh

  AlignCoord = ListGetInteger( SolverParams,'Align Coordinate',Found )
  IF( Found ) THEN
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
      CALL Fatal(Caller,'Invalid value for "Align Coordinate": '//TRIM(I2S(AlignCoord)))
    END IF
      
    IF( AlignCoord > 0 ) THEN
      dx = MINVAL( x2 ) - MAXVAL( x1 ) 
    ELSE
      dx = MINVAL( x1 ) - MAXVAL( x2 ) 
    END IF
    WRITE(Message,'(A,ES12.3)') 'Aligning coordinate '//TRIM(I2S(k))//' with: ',dx
    CALL Info(Caller,Message)
    
    x1 = x1 + dx
  END IF

  
  Str = ListGetString( SolverParams,'Mask Name',Found) 
  IF( Found ) THEN
    ALLOCATE( NewMaskPerm( ToMesh % NumberOfNodes ) )
    NewMaskPerm = 0
    CALL MakePermUsingMask( Model, Solver, ToMesh, Str, .FALSE., &
        NewMaskPerm, n, RequireLogical = .TRUE. )
    IF( n == 0 ) THEN
      CALL Fatal(Caller,'Zero masked nodes')
    ELSE
      CALL Info(Caller,'Number of masked nodes: '//TRIM(I2S(n)))
    END IF

    CALL InterpolateMeshToMeshQ( FromMesh, ToMesh, FromMesh % Variables, ToMesh % Variables, &
        NewMaskPerm = NewMaskPerm ) 
  ELSE
    CALL InterpolateMeshToMeshQ( FromMesh, ToMesh, FromMesh % Variables, ToMesh % Variables )
  END IF

  CALL Info(Caller,'Interpolation complete')
  
!------------------------------------------------------------------------------
END SUBROUTINE GmshOutputReader
!------------------------------------------------------------------------------
  
