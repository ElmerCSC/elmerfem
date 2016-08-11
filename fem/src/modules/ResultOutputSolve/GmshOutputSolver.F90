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
!> Saves results in ascii format understood by the pre-/postprocessing software Gmsh.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE GmshOutputSolver( Model,Solver,dt,TransientSimulation )
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
  REAL(KIND=dp), POINTER :: Values(:),Values2(:),Values3(:)
  REAL(KIND=dp) :: Vector(3), Time
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(Variable_t), POINTER :: Solution, TimeVariable
  TYPE(ValueList_t), POINTER :: SolverParams
  
  LOGICAL :: Found, GotField, FileAppend, AlterTopology, MaskExists
  LOGICAL :: EigenAnalysis = .FALSE., EigenActive, ComponentVector
  
  INTEGER :: VisitedTimes = 0, ExtCount
  INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerCode, GmshCode,body_id, Vari, Rank, truedim
  INTEGER :: Tag, NumberOfAllElements, BCOffSet
  INTEGER, PARAMETER :: MaxElemCode = 827
  INTEGER :: ElmerToGmshType(MaxElemCode), GmshToElmerType(21), GmshIndexes(27) 
  INTEGER, POINTER :: NodeIndexes(:), ElmerIndexes(:), MaskPerm(:)

  INTEGER, PARAMETER :: LENGTH = 1024
  CHARACTER(LEN=LENGTH) :: OutputFile, Txt, FieldName, CompName
    
  SAVE VisitedTimes
  
!------------------------------------------------------------------------------

  CALL Info('GmshOutputSolver','Saving results in Gmsh format')

  ExtCount = ListGetInteger( Solver % Values,'Output Count',Found)
  IF( Found ) THEN
    VisitedTimes = ExtCount
  ELSE
    VisitedTimes = VisitedTimes + 1
  END IF

  GmshToElmerType = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
      510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)
  ElmerToGmshType = 0

  DO i=1,SIZE(GmshToElmerType)
    j = GmshToElmerType(i)
    IF( j > 0 ) ElmerToGmshType(j) = i
  END DO

  SolverParams => GetSolverParams()
  EigenAnalysis = GetLogical( SolverParams, 'Eigen Analysis', Found )
  FileAppend = GetLogical( SolverParams,'File Append',Found)
  IF(.NOT. Found) FileAppend = .TRUE.
  AlterTopology = GetLogical( SolverParams,'Alter Topology',Found)
  
  Txt = ListGetString( SolverParams,'Mask Variable',MaskExists)
  IF( MaskExists ) THEN
    Solution => VariableGet(Model % Variables,TRIM(Txt))
    IF( ASSOCIATED(Solution)) MaskPerm => Solution % Perm
    MaskExists = ASSOCIATED(MaskPerm)
  END IF

  OutputFile = GetString( Solver % Values, 'Output File Name', Found )
  IF( .NOT.Found ) OutputFile = 'Output.msh'

  IF(INDEX(OutputFile,'.') == 0) WRITE( OutputFile,'(A,A)') TRIM(OutputFile),".msh"

  dim = CoordinateSystemDimension()
  IF( VisitedTimes > 1 ) THEN
    IF( AlterTopology ) THEN
      OutputFile = NextFreeFilename( OutputFile )
      CALL Info('GmshOutputSolver','Writing mesh and data to a new file: '//TRIM(OutputFile))
    ELSE IF( FileAppend ) THEN      
      CALL Info('GmshOutputSolver','Appending data to the same file: '//TRIM(OutputFile))
      OPEN(UNIT=10, FILE=OutputFile, POSITION='APPEND' )      
      GOTO 10
    ELSE
      OutputFile = NextFreeFilename( OutputFile )          
      CALL Info('GmshOutputSolver','Writing data to a new file: '//TRIM(OutputFile))
      OPEN(UNIT=10, FILE=OutputFile )
      WRITE(10,'(A)') '$MeshFormat'
      WRITE(10,'(A)') '2.0 0 8'
      WRITE(10,'(A)') '$EndMeshFormat'          
      GOTO 10    
    END IF
  END IF


  ! Save the header
  !-------------------------------------------------
  CALL Info('GsmhOutputSolver','Saving results to file: '//TRIM(OutputFile))
  OPEN(UNIT=10, FILE=OutputFile )
  
  WRITE(10,'(A)') '$MeshFormat'
  WRITE(10,'(A)') '2.0 0 8'
  WRITE(10,'(A)') '$EndMeshFormat'    
  

  ! Save the mesh nodes
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the mesh nodes')
  IF( MaskExists ) THEN
    nsize = MAXVAL( MaskPerm ) 
  ELSE
    nsize = Model % NumberOfNodes
  END IF

  WRITE(10,'(A)') '$Nodes'
  WRITE(10,'(I8)') nsize
  IF( dim == 3 ) THEN
    DO i = 1, Model % NumberOfNodes
      IF( MaskExists ) THEN
        IF( MaskPerm(i) == 0 ) CYCLE
      END IF      
      WRITE(10,'(I8,3ES16.7E3)') i,Model % Nodes % x(i),Model % Nodes % y(i), Model % Nodes % z(i)
    END DO
  ELSE 
    DO i = 1, Model % NumberOfNodes
      IF( MaskExists ) THEN
        IF( MaskPerm(i) == 0 ) CYCLE
      END IF            
      WRITE(10,'(I8,2ES16.7E3,A)') i,Model % Nodes % x(i),Model % Nodes % y(i),' 0.0' 
    END DO
  END IF
  WRITE(10,'(A)') '$EndNodes'

  ! Save the mesh elements
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the mesh elements')
  NumberOfAllElements = Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
  
  IF( MaskExists ) THEN
    nsize = 0
    DO i=1,NumberOfAllElements
      Element => Model % Mesh % Elements(i)
      ElmerIndexes => Element % NodeIndexes
      IF( ANY(MaskPerm(ElmerIndexes) == 0) ) CYCLE
      nsize = nsize  + 1
    END DO
  ELSE
    nsize = NumberOfAllElements
  END IF

  BCOffSet = 100
  DO WHILE( BCOffset <= Model % NumberOfBodies ) 
    BCOffset = 10 * BCOffset
  END DO

  WRITE(10,'(A)') '$Elements'
  WRITE(10,'(I8)') nsize
  DO i = 1, NumberOfAllElements
    Element => Model % Mesh % Elements(i)
    ElmerCode = Element % TYPE % ElementCode
    ElmerIndexes => Element % NodeIndexes
    
    IF( MaskExists ) THEN
      IF( ANY(MaskPerm(ElmerIndexes) == 0) ) CYCLE     
    END IF

    GmshCode = ElmerToGmshType(ElmerCode)
    IF( GmshCode == 0 ) THEN
      CALL Warn('GmshOutputSolver','Gmsh element index not found!')
      CYCLE
    END IF

    IF( i <= Model % NumberOfBulkElements ) THEN
      Tag = Element % BodyId
    ELSE
      Tag = GetBCId( Element ) + BCOffset
    END IF

    WRITE(10,'(I8,I3,I3,I5,I5)',ADVANCE='NO') i,GmshCode,2,Tag,Tag
    k = MOD(ElmerCode,100)

    CALL ElmerToGmshIndex(ElmerCode,ElmerIndexes,GmshIndexes)

    DO j=1,k-1
      WRITE(10,'(I8)',ADVANCE='NO') GmshIndexes(j)
    END DO
    WRITE(10,'(I8)') GmshIndexes(k)
  END DO
  WRITE(10,'(A)') '$EndElements'

  ! With a mask the list of physical entities should be checked
  !-------------------------------------------------------------
  IF(.NOT. MaskExists ) THEN
    nsize = Model % NumberOfBodies + Model % NumberOfBCs
    WRITE(10,'(A)') '$PhysicalNames'
    WRITE(10,'(I8)') nsize
    DO i=1,Model % NumberOfBodies 
      Txt = ListGetString( Model % Bodies(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(10,'(I8,A)') i,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(10,'(I8,A,I0,A)') i,'"Body ',i,'"'       
      END IF
    END DO
    DO i=1,Model % NumberOfBCs
      Txt = ListGetString( Model % BCs(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(10,'(I8,A)') i+BCOffset,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(10,'(I8,A,I0,A)') i+BCOffset,'"Boundary Condition ',i,'"'               
      END IF
    END DO
    WRITE(10,'(A)') '$EndPhysicalNames'
  END IF


10 CONTINUE


  ! Time is needed
  !-------------------------------------------------
  TimeVariable => VariableGet( Model % Variables, 'Time' )        
  Time = TimeVariable % Values(1)
  
  ! Loop over different type of variables
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the nodal data')
  DO Rank = 0,2
    DO Vari = 1, 999
      IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
      IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
      IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

      FieldName = GetString( Solver % Values, TRIM(Txt), Found )
      IF(.NOT. Found) EXIT 
      IF( Rank == 2) THEN
        CALL Warn('GmshOutputSolver','Not implemented yet for tensors!')
        CYCLE
      END IF

      ComponentVector = .FALSE.
      Solution => VariableGet( Model % Mesh % Variables, FieldName )
      IF(ASSOCIATED(Solution)) THEN
        Values => Solution % Values
        Perm => Solution % Perm
        dofs = Solution % DOFs
      ELSE
        IF( Rank == 1 ) THEN
          Solution => VariableGet( Model % Mesh % Variables, FieldName//' 1' )
          IF( ASSOCIATED( Solution ) ) THEN
            ComponentVector = .TRUE.
            Values => Solution % Values
            Perm => Solution % Perm
            dofs = 1
            Solution => VariableGet( Model % Mesh % Variables, FieldName//' 2' )
            IF( ASSOCIATED(Solution)) THEN
              Values2 => Solution % Values
              dofs = 2
            END IF            
            Solution => VariableGet( Model % Mesh % Variables, FieldName//' 3' )
            IF( ASSOCIATED(Solution)) THEN
              Values3 => Solution % Values
              dofs = 3
            END IF
          END IF
        END IF
        IF( .NOT. ASSOCIATED(Solution)) THEN
          CALL Warn('GsmhOutputSolver','Variable not present: '//TRIM(FieldName))
          CYCLE
        END IF
      END IF
      IF( ASSOCIATED(Solution % EigenVectors) ) THEN
        CALL Warn('GmshOutputSolver','Eigenvectors related to field: '//TRIM(FieldName))
        CALL Warn('GmshOutputSolver','Eigenvectors saving yet not supported')
      END IF

      truedim = MIN(dofs, dim)
      IF( MaskExists ) THEN
        nsize = 0
        DO i=1,SIZE(Perm)
          IF( Perm(i)==0 .OR. MaskPerm(i) == 0 ) CYCLE
          nsize = nsize + 1
        END DO
        IF( nsize == 0 ) THEN
          CALL Warn('GmshOutputSolver','No dofs with the current mask for saving: '//TRIM(FieldName))         
        END IF
      ELSE
        nsize = SIZE(Values) / Dofs
      END IF

      
      WRITE(10,'(A)') '$NodeData'
      WRITE(10,'(A)') '1'
      WRITE(10,'(A)') '"'//TRIM(FieldName)//'"'
      WRITE(10,'(A)') '1'

      ! Gmsh starts steady state indexes from zero, hence deductions by one
      IF( TransientSimulation ) THEN
        WRITE(10,'(ES16.7E3)') Time
      ELSE
        WRITE(10,'(ES16.7E3)') Time - 1.0_dp
      END IF
      WRITE(10,'(A)') '3'
      WRITE(10,'(I8)') VisitedTimes-1
      IF(Rank == 0) THEN
        WRITE(10,'(A)') '1'
      ELSE IF(Rank == 1) THEN
        WRITE(10,'(A)') '3'
      ELSE 
        WRITE(10,'(A)') '9'
      END IF     
      WRITE(10,'(I8)') nsize
     
      DO i=1,SIZE(Perm) 
        j = Perm(i)
        IF( j == 0) CYCLE
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        
        IF( Rank == 0 ) THEN
          WRITE(10,'(I8,ES16.7E3)') i,Values(j)
        ELSE IF(Rank == 1) THEN
          IF( ComponentVector ) THEN
            IF( truedim == 2 ) THEN
              WRITE(10,'(I8,2ES16.7E3,A)') i,&
                  Values(j),Values2(j),' 0.0'
            ELSE
              WRITE(10,'(I8,3ES16.7E3)') i,&
                  Values(j),Values2(j),Values3(j)
            END IF
          ELSE
            IF( truedim == 2 ) THEN
              WRITE(10,'(I8,2ES16.7E3,A)') i,&
                  Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),' 0.0'
            ELSE
              WRITE(10,'(I8,3ES16.7E3)') i,&
                  Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),Values(dofs*(j-1)+3)
            END IF           
          END IF
        END IF
      END DO
      WRITE(10,'(A)') '$EndNodeData'

    END DO
  END DO
  
      
  IF(.FALSE.) THEN
    WRITE(10,'(A)') '$ElementData'
    WRITE(10,'(A)') '$EndElementData'
  END IF
  
  IF(.FALSE.) THEN
    WRITE(10,'(A)') '$ElementNodeData'
    WRITE(10,'(A)') '$EndElementNodeData'
  END IF
  
  CLOSE(10)
  
  CALL Info('GmshOutputSolver','Gmsh output complete')

CONTAINS



  SUBROUTINE ElmerToGmshIndex(Code,ElmerIndexes,GmshIndexes)

    INTEGER :: Code
    INTEGER :: ElmerIndexes(:),GmshIndexes(:)
    INTEGER :: i,n
    LOGICAL :: reorder, Visited = .FALSE.

    INTEGER, TARGET :: order510(10),order613(13),order715(15),order820(20)
    INTEGER, POINTER :: order(:)

    SAVE Visited

    IF(.NOT. Visited ) THEN
      order510(:) = (/ 0,1,2,3,4,5,6,7,9,8 /)
      order613(:) = (/ 0,1,2,3,4,5,8,10,6,7,9,11,12 /)
      order715(:) = (/ 0,1,2,3,4,5,6,9,7,8,10,11,12,14,13 /)
      order820(:) = (/ 0,1,2,3,4,5,6,7,8,11,12,9,10,12,14,15,16,18,19,17 /)
      Visited = .TRUE.
    END IF

    reorder = .FALSE.

    SELECT CASE( Code )
      
    CASE (510)
      reorder = .TRUE.
      order => order510
      
    CASE (613)
      reorder = .TRUE.
      order => order613
      
    CASE (715)
      reorder = .TRUE.
      order => order715
      
    CASE (820)
      reorder = .TRUE.
      order => order820
     
    CASE DEFAULT
      
    END SELECT

    n = MOD(Code,100) 
    IF( reorder ) THEN
      DO i=1,n 
        GmshIndexes(order(i)+1) = ElmerIndexes(i)
      END DO
    ELSE
      GmshIndexes(1:n) = ElmerIndexes(1:n)      
    END IF


  END SUBROUTINE ElmerToGmshIndex

!------------------------------------------------------------------------------
END SUBROUTINE GmshOutputSolver
!------------------------------------------------------------------------------
  
