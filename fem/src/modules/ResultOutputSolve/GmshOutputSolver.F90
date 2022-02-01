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
  REAL(KIND=dp), POINTER :: Values(:),Values2(:),Values3(:)
  REAL(KIND=dp) :: Vector(3), Time
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(Variable_t), POINTER :: Solution, TimeVariable
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh
  
  LOGICAL :: Found, GotField, FileAppend, AlterTopology, MaskExists
  LOGICAL :: EigenAnalysis = .FALSE., EigenActive, ComponentVector, Parallel
  
  INTEGER :: VisitedTimes = 0, ExtCount
  INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerCode, GmshCode,body_id, Vari, Rank, truedim
  INTEGER :: Tag, NumberOfAllElements, BCOffSet
  INTEGER, PARAMETER :: MaxElemCode = 827
  INTEGER :: ElmerToGmshType(MaxElemCode), GmshToElmerType(21), &
      ElmerIndexes(27), GmshIndexes(27) 
  INTEGER, POINTER :: NodeIndexes(:)

  INTEGER, ALLOCATABLE :: NodePerm(:),DgPerm(:)
  INTEGER, ALLOCATABLE, TARGET :: InvDgPerm(:), InvNodePerm(:)
  LOGICAL, ALLOCATABLE :: ActiveElem(:)
  LOGICAL :: NoPermutation
  INTEGER :: NumberOfGeomNodes, NumberOfDofNodes,NumberOfElements, ElemFirst, ElemLast
  INTEGER, POINTER :: InvFieldPerm(:)
  
  INTEGER, PARAMETER :: LENGTH = 1024
  CHARACTER(LEN=LENGTH) :: Txt, FieldName, CompName
  CHARACTER(MAX_NAME_LEN) :: OutputFile, OutputDirectory
  INTEGER :: GmshUnit
  CHARACTER(*), PARAMETER :: Caller = 'GmshOutputSolver'
    
  SAVE VisitedTimes
  
!------------------------------------------------------------------------------

  CALL Info(Caller,'Saving results in Gmsh format')

  Mesh => Model % Mesh
  Params => Solver % Values
  Parallel = ( ParEnv % PEs > 1 )
  
  ExtCount = ListGetInteger( Params,'Output Count',Found)
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

  EigenAnalysis = GetLogical( Params, 'Eigen Analysis', Found )
  FileAppend = GetLogical( Params,'File Append',Found)
  IF(.NOT. Found) FileAppend = .TRUE.
  AlterTopology = GetLogical( Params,'Alter Topology',Found)
  
  OutputFile = GetString( Solver % Values, 'Output File Name', Found )
  IF( Found ) THEN
    IF(INDEX(OutputFile,'.') == 0) WRITE( OutputFile,'(A,A)') TRIM(OutputFile),".msh"
  ELSE
    OutputFile = 'Output.msh'
  END IF
    
  CALL SolverOutputDirectory( Solver, OutputFile, OutputDirectory, UseMeshDir = .TRUE. )
  OutputFile = TRIM(OutputDirectory)// '/' //TRIM(OutputFile)
  
  !------------------------------------------------------------------------------
  ! Initialize stuff for masked saving
  !------------------------------------------------------------------------------
  CALL GenerateSaveMask(Mesh,Params,Parallel,0,.FALSE.,&
      NodePerm,ActiveElem,NumberOfGeomNodes,NumberOfElements,&
      ElemFirst,ElemLast)
  
  CALL GenerateSavePermutation(Mesh,.FALSE.,.FALSE.,0,.FALSE.,&
      ActiveElem,NumberOfGeomNodes,NoPermutation,NumberOfDofNodes,&
      DgPerm,InvDgPerm,NodePerm,InvNodePerm)
  
  InvFieldPerm => InvNodePerm
  
  dim = CoordinateSystemDimension()
  IF( VisitedTimes > 1 ) THEN
    IF( AlterTopology ) THEN
      OutputFile = NextFreeFilename( OutputFile )
      CALL Info(Caller,'Writing mesh and data to a new file: '//TRIM(OutputFile))
    ELSE IF( FileAppend ) THEN      
      CALL Info(Caller,'Appending data to the same file: '//TRIM(OutputFile))
      OPEN(NEWUNIT=GmshUnit, FILE=OutputFile, POSITION='APPEND' )      
      GOTO 10
    ELSE
      OutputFile = NextFreeFilename( OutputFile )          
      CALL Info(Caller,'Writing data to a new file: '//TRIM(OutputFile))
      OPEN(NEWUNIT=GmshUnit, FILE=OutputFile )
      WRITE(GmshUnit,'(A)') '$MeshFormat'
      WRITE(GmshUnit,'(A)') '2.0 0 8'
      WRITE(GmshUnit,'(A)') '$EndMeshFormat'          
      GOTO 10    
    END IF
  END IF


  ! Save the header
  !-------------------------------------------------
  CALL Info('GsmhOutputSolver','Saving results to file: '//TRIM(OutputFile))
  OPEN(NEWUNIT=GmshUnit, FILE=OutputFile )
  
  WRITE(GmshUnit,'(A)') '$MeshFormat'
  WRITE(GmshUnit,'(A)') '2.0 0 8'
  WRITE(GmshUnit,'(A)') '$EndMeshFormat'    
  

  ! Save the mesh nodes
  !-------------------------------------------------
  CALL Info(Caller,'Writing the mesh nodes')
  CALL WriteGmshNodes()

  ! Save the mesh elements
  !-------------------------------------------------
  CALL Info(Caller,'Writing the mesh elements')
  CALL WriteGmshElements() 

  ! With a mask the list of physical entities should be checked
  !-------------------------------------------------------------
  IF(.NOT. MaskExists ) THEN
!    CALL WritePhysicalNames() 
  END IF

10 CONTINUE

  CALL Info(Caller,'Writing the nodal data')
  CALL WriteGmshData()
      
  IF(.FALSE.) THEN
    WRITE(GmshUnit,'(A)') '$ElementData'
    WRITE(GmshUnit,'(A)') '$EndElementData'
  END IF
  
  IF(.FALSE.) THEN
    WRITE(GmshUnit,'(A)') '$ElementNodeData'
    WRITE(GmshUnit,'(A)') '$EndElementNodeData'
  END IF
  
  CLOSE(GmshUnit)



  IF(ALLOCATED(DgPerm)) DEALLOCATE(DgPerm)
  IF(ALLOCATED(InvDgPerm)) DEALLOCATE(InvDgPerm)
  IF(ALLOCATED(NodePerm)) DEALLOCATE(NodePerm)
  IF(ALLOCATED(InvNodePerm)) DEALLOCATE(InvNodePerm)
  IF(ALLOCATED(ActiveElem)) DEALLOCATE(ActiveElem)

  
  CALL Info(Caller,'Gmsh output complete')

CONTAINS
  
  SUBROUTINE WriteGmshNodes()

    nsize = NumberOfGeomNodes
    
    WRITE(GmshUnit,'(A)') '$Nodes'
    WRITE(GmshUnit,'(I8)') nsize
    DO i = 1, nsize
      IF( NoPermutation ) THEN
        j = i 
      ELSE
        j = InvNodePerm(i)
      END IF
      
      IF( dim == 3 ) THEN
        WRITE(GmshUnit,'(I8,3ES16.7E3)') i,Mesh % Nodes % x(j),Mesh % Nodes % y(j), Mesh % Nodes % z(j)
      ELSE
        WRITE(GmshUnit,'(I8,2ES16.7E3,A)') i,Mesh % Nodes % x(j),Mesh % Nodes % y(j),' 0.0' 
      END IF
    END DO
    WRITE(GmshUnit,'(A)') '$EndNodes'
  END SUBROUTINE WriteGmshNodes
    
  
  SUBROUTINE WriteGmshElements()

    nsize = NumberOfElements 

    BCOffSet = 100
    DO WHILE( BCOffset <= Model % NumberOfBodies ) 
      BCOffset = 10 * BCOffset
    END DO

    WRITE(GmshUnit,'(A)') '$Elements'
    WRITE(GmshUnit,'(I8)') nsize

    l = 0
    DO i = ElemFirst, ElemLast
      IF(.NOT. ActiveElem(i) ) CYCLE

      l = l + 1
      Element => Mesh % Elements(i)
      ElmerCode = Element % TYPE % ElementCode

      n = Element % Type % NumberOfNodes
      IF( NoPermutation ) THEN
        ElmerIndexes(1:n) = Element % NodeIndexes(1:n)
      ELSE
        ElmerIndexes(1:n) = NodePerm(Element % NodeIndexes(1:n))
      END IF
        
      GmshCode = ElmerToGmshType(ElmerCode)
      IF( GmshCode == 0 ) THEN
        CALL Warn(Caller,'Gmsh element index not found!')
        CYCLE
      END IF

      IF( i <= Model % NumberOfBulkElements ) THEN
        Tag = Element % BodyId
      ELSE
        Tag = GetBCId( Element ) + BCOffset
      END IF

      WRITE(GmshUnit,'(I8,I3,I3,I5,I5)',ADVANCE='NO') l,GmshCode,2,Tag,Tag
      k = MOD(ElmerCode,100)

      CALL ElmerToGmshIndex(ElmerCode,ElmerIndexes,GmshIndexes)

      DO j=1,k-1
        WRITE(GmshUnit,'(I8)',ADVANCE='NO') GmshIndexes(j)
      END DO
      WRITE(GmshUnit,'(I8)') GmshIndexes(k)
    END DO
    WRITE(GmshUnit,'(A)') '$EndElements'
  END SUBROUTINE WriteGmshElements


  SUBROUTINE WritePhysicalNames()
    CALL Info(Caller,'Writing the physical entity names')
    nsize = Model % NumberOfBodies + Model % NumberOfBCs
    WRITE(GmshUnit,'(A)') '$PhysicalNames'
    WRITE(GmshUnit,'(I8)') nsize
    DO i=1,Model % NumberOfBodies 
      Txt = ListGetString( Model % Bodies(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(GmshUnit,'(I8,A)') i,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(GmshUnit,'(I8,A,I0,A)') i,'"Body ',i,'"'       
      END IF
    END DO
    DO i=1,Model % NumberOfBCs
      Txt = ListGetString( Model % BCs(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(GmshUnit,'(I8,A)') i+BCOffset,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(GmshUnit,'(I8,A,I0,A)') i+BCOffset,'"Boundary Condition ',i,'"'               
      END IF
    END DO
    WRITE(GmshUnit,'(A)') '$EndPhysicalNames'
  END SUBROUTINE WritePhysicalNames
    
  
  
  SUBROUTINE WriteGmshData()
    INTEGER :: ii

    
    ! Time is needed
    !-------------------------------------------------
    TimeVariable => VariableGet( Model % Variables, 'Time' )        
    Time = TimeVariable % Values(1)

    ! Loop over different type of variables
    !-------------------------------------------------
    CALL Info(Caller,'Writing the nodal data')
    DO Rank = 0,2
      DO Vari = 1, 999
        IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

        FieldName = GetString( Solver % Values, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT 
        IF( Rank == 2) THEN
          CALL Warn(Caller,'Not implemented yet for tensors!')
          CYCLE
        END IF

        ComponentVector = .FALSE.
        Solution => VariableGet( Mesh % Variables, FieldName )
        IF(ASSOCIATED(Solution)) THEN
          Values => Solution % Values
          Perm => Solution % Perm
          dofs = Solution % DOFs
        ELSE
          IF( Rank == 1 ) THEN
            Solution => VariableGet( Mesh % Variables, FieldName//' 1' )
            IF( ASSOCIATED( Solution ) ) THEN
              ComponentVector = .TRUE.
              Values => Solution % Values
              Perm => Solution % Perm
              dofs = 1
              Solution => VariableGet( Mesh % Variables, FieldName//' 2' )
              IF( ASSOCIATED(Solution)) THEN
                Values2 => Solution % Values
                dofs = 2
              END IF
              Solution => VariableGet( Mesh % Variables, FieldName//' 3' )
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
          CALL Warn(Caller,'Eigenvectors related to field: '//TRIM(FieldName))
          CALL Warn(Caller,'Eigenvectors saving yet not supported')
        END IF

        truedim = MIN(dofs, dim)
        nsize = NumberOfGeomNodes

        WRITE(GmshUnit,'(A)') '$NodeData'
        WRITE(GmshUnit,'(A)') '1'
        WRITE(GmshUnit,'(A)') '"'//TRIM(FieldName)//'"'
        WRITE(GmshUnit,'(A)') '1'
        
        ! Gmsh starts steady state indexes from zero, hence deductions by one
        IF( TransientSimulation ) THEN
          WRITE(GmshUnit,'(ES16.7E3)') Time
        ELSE
          WRITE(GmshUnit,'(ES16.7E3)') Time - 1.0_dp
        END IF
        WRITE(GmshUnit,'(A)') '3'
        WRITE(GmshUnit,'(I8)') VisitedTimes-1
        IF(Rank == 0) THEN
          WRITE(GmshUnit,'(A)') '1'
        ELSE IF(Rank == 1) THEN
          WRITE(GmshUnit,'(A)') '3'
        ELSE 
          WRITE(GmshUnit,'(A)') '9'
        END IF
        WRITE(GmshUnit,'(I8)') nsize

        DO ii = 1, NumberOfGeomNodes
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
          
          IF( Rank == 0 ) THEN
            WRITE(GmshUnit,'(I8,ES16.7E3)') ii,Values(j)
          ELSE IF(Rank == 1) THEN
            IF( j == 0 ) THEN
              WRITE(GmshUnit,'(I8,A)') ii,' 0.0 0.0 0.0'                
            ELSE IF( ComponentVector ) THEN
              IF( truedim == 2 ) THEN
                WRITE(GmshUnit,'(I8,2ES16.7E3,A)') ii,&
                    Values(j),Values2(j),' 0.0'
              ELSE
                WRITE(GmshUnit,'(I8,3ES16.7E3)') ii,&
                    Values(j),Values2(j),Values3(j)
              END IF
            ELSE
              IF( truedim == 2 ) THEN
                WRITE(GmshUnit,'(I8,2ES16.7E3,A)') ii,&
                    Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),' 0.0'
              ELSE
                WRITE(GmshUnit,'(I8,3ES16.7E3)') ii,&
                    Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),Values(dofs*(j-1)+3)
              END IF
            END IF
          END IF
        END DO
        WRITE(GmshUnit,'(A)') '$EndNodeData'

      END DO
    END DO
  END SUBROUTINE WriteGmshData
    

   
!------------------------------------------------------------------------------
END SUBROUTINE GmshOutputSolver
!------------------------------------------------------------------------------
  
