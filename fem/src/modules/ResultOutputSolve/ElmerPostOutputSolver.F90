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
!> A fork of the ElmerPost output utility which also includes use of masks.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE ElmerPostOutputSolver( Model, Solver,dt,TransientSimulation,ONOEfound )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL, OPTIONAL :: ONOEfound
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    INTEGER :: TimeCount
    CHARACTER(LEN=512) :: PostFile, FilePrefix
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(Variable_t), POINTER :: Var,Var1,Displacement=>NULL(),MeshUpdate=>NULL(),MaskVar
    
    CHARACTER(LEN=512) :: Row
    CHARACTER(MAX_NAME_LEN) :: Str, DateStr, VarName, Txt
    
    LOGICAL :: gotIt, FreeSurfaceFlag, MoveBoundary, MeshMoved, MaskExists, Found
    INTEGER :: ScalarFields, VectorFields, ONOECount
    REAL(KIND=dp) :: Time,MeshScale, NodeCoords(3)
    REAL(KIND=dp), POINTER :: Values(:), Values2(:), Values3(:)
    INTEGER :: jj,n1,n2,ii,i,j,k,l,n,m,q,Node,DOFs,TimeStep, On_nodes_on_elements, &
         NumberOfNodes, NumberOfElements, ind, Vari, MeshDim, ExtCount
    INTEGER, POINTER :: MaskPerm(:), MaskOrder(:), TimeSteps(:)
    LOGICAL :: MaskAllocated 
    TYPE(Mesh_t), POINTER :: Mesh

    INTEGER, SAVE :: ParallelNodes, SavedCount = 0
    LOGICAL, SAVE :: Visited = .FALSE.
!------------------------------------------------------------------------------
    
    MeshDim = Model % Mesh % MeshDim
    Mesh => Model % Mesh

    IF ( .NOT.PRESENT(ONOEfound) ) THEN
      ExtCount = GetInteger( Solver % Values,'Output Count',GotIt)
    
      IF( GotIt ) THEN
        SavedCount = ExtCount 
      ELSE
        SavedCount = SavedCount + 1
      END IF
      Visited = ( SavedCount > 1 )
    END IF

    Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
        'Timestep Intervals', GotIt )
    IF( GotIt ) THEN
      TimeCount = SUM(TimeSteps)
    ELSE 
      TimeCount = GetInteger( Model % Simulation,'Steady State Max Iterations')
    END IF
    
    FilePrefix = GetString( Solver % Values,'Output File Name',GotIt )
    IF ( .NOT.GotIt ) FilePrefix = "Output"
    
    IF(INDEX(FilePrefix,'.') == 0) THEN
      WRITE( Postfile,'(A,A)') TRIM(FilePrefix),".ep"
    ELSE
      PostFile = FilePrefix
    END IF

    IF ( PRESENT(ONOEfound) ) PostFile = TRIM(Postfile) // ".el"
    
    IF ( INDEX( PostFile, ':') == 0 .AND. PostFile(1:1) /= '/' .AND. &
        PostFile(1:1) /= Backslash ) THEN
      
      IF ( LEN_TRIM(OutputPath) > 0 ) THEN
        IF ( Visited )  THEN
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
              TRIM(PostFile), POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
              TRIM(PostFile),STATUS='UNKNOWN' )
        END IF
      ELSE
        IF ( Visited ) THEN
          OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
        ENDIF
      END IF
    ELSE
      IF ( Visited  ) THEN
        OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
      ELSE
        OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
      END IF
    END IF


    FreeSurfaceFlag = .FALSE.
    MoveBoundary    = .FALSE.
    DO i=1,CurrentModel % NumberOfBCs
      FreeSurfaceFlag = FreeSurfaceFlag .OR. GetLogical( &
          CurrentModel % BCs(i) % Values,'Free Surface', GotIt )
      IF ( FreeSurfaceFlag ) THEN
        MoveBoundary =  GetLogical( &
            CurrentModel % BCs(i) % Values,'Internal Move Boundary', GotIt )
        
        IF ( .NOT. GotIt ) MoveBoundary = .TRUE.
        
        FreeSurfaceFlag = FreeSurfaceFlag .AND. MoveBoundary
      END IF
      
      IF ( FreeSurfaceFlag ) EXIT
    END DO
    
!------------------------------------------------------------------------------
! Initialize stuff for masked saving
!------------------------------------------------------------------------------
    Str = GetString( Solver % Values,'Mask Variable',MaskExists)
    MaskAllocated = .FALSE.
    IF( MaskExists ) THEN
      MaskVar => VariableGet(Solver % Mesh % Variables,TRIM(Str))
      IF( ASSOCIATED(MaskVar)) THEN
        MaskPerm => MaskVar % Perm
        MaskExists = ASSOCIATED(MaskPerm)
      END IF
    ELSE
      IF( MeshDim == 2 ) THEN
        Str = GetString( Solver % Values,'2D Mask Name',GotIt)    
      ELSE IF( MeshDim == 3 ) THEN  
        Str = GetString( Solver % Values,'3D Mask Name',GotIt)    
      END IF
      IF(.NOT. GotIt) Str = GetString( Solver % Values,'Mask Name',GotIt) 
      IF( GotIt ) THEN
        ALLOCATE( MaskPerm( Model % NumberOfNodes ) ) 
        CALL MakePermUsingMask( Model,Solver,Mesh,Str, &
            .FALSE., MaskPerm, NumberOfNodes )
        ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfNodes ) )
        IF( ParallelNodes > 0 ) THEN
          MaskExists = .TRUE.
          MaskAllocated = .TRUE.
        ELSE
          CALL Warn('ElmerPostOutputSolver','Given mask not active: '//TRIM(Str) )
          DEALLOCATE( MaskPerm )
        END IF
      END IF
    END IF
    

    IF(MaskExists) THEN
      CALL Info('ElmerPostOutputSolver','Using > '// TRIM(Str) // ' < as mask variable')
      NumberOfNodes = MAXVAL(MaskPerm)
      ALLOCATE(MaskOrder(NumberOfNodes))
      DO i=1,SIZE(MaskPerm)
        j = MaskPerm(i)
        IF(j > 0) MaskOrder(j) = i
      END DO
      NumberOfElements = 0
      IF ( PRESENT(ONOEfound) ) NumberOfNodes=0
      DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
        Element => Model % Elements(i)
        IF( ALL(MaskPerm(Element % NodeIndexes)>0)) THEN
          NumberOfElements = NumberOfElements + 1
        END IF
        IF ( PRESENT(ONOEfound) ) &
          NumberOfNodes=NumberOfNodes+Element % TYPE % NumberOfNodes
      END DO
    ELSE
      NumberOfNodes = Model % NumberOfNodes
      NumberOfElements =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
      IF ( PRESENT(ONOEfound) ) THEN
        NumberOfNodes=0
        DO i=1,Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(i)
          NumberOfNodes=NumberOfNodes+Element % TYPE % NumberOfNodes
        END DO
      END IF
    END IF
 
!------------------------------------------------------------------------------
!   Count degrees of freedom to be saved
!------------------------------------------------------------------------------
    
    Dofs = 0
    ScalarFields = 0
    On_nodes_on_elements = 0
    DO Vari = 1, 999
      WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
      VarName=ListGetString( Solver % Values, Txt, Found )
      IF ( .NOT. Found ) EXIT
      Var => VariableGet( Solver % Mesh % Variables, Varname )
      IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
        On_Nodes_on_elements=On_Nodes_on_elements+1
        IF ( PRESENT(ONOEfound) ) Dofs=Dofs+1
      ELSE
        IF ( .NOT.PRESENT(ONOEfound) ) Dofs=Dofs+1
      END IF
      ScalarFields = ScalarFields + 1
    END DO
    
    VectorFields = 0
    DO Vari = 1, 999
      WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
      VarName=ListGetString( Solver % Values, Txt, Found )
      IF ( .NOT. Found ) EXIT
      Var => VariableGet( Solver % Mesh % Variables, Varname )
      IF(.NOT.ASSOCIATED(Var)) & 
        Var => VariableGet( Solver % Mesh % Variables, TRIM(Varname)//' 1' )
      IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
        On_Nodes_on_elements=On_Nodes_on_elements+1
        IF ( PRESENT(ONOEfound) ) Dofs=Dofs+3
      ELSE
        IF ( .NOT.PRESENT(ONOEfound) ) Dofs=Dofs+3
      END IF
      VectorFields = VectorFields+1
    END DO
    
    IF( ScalarFields + VectorFields == 0 ) GOTO 10
    
!------------------------------------------------------------------------------
! Write header to output
!------------------------------------------------------------------------------
    IF ( .NOT. Visited ) THEN
      WRITE(PostFileUnit,'(i10,i10,i7,i7)',ADVANCE='NO' ) NumberOfNodes, &
          NumberOfElements, DOFs, TimeCount
      
      DO Vari = 1, ScalarFields
        WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        VarName = ListGetString( Solver % Values, TRIM(Txt) )
        Var => VariableGet( Solver % Mesh % Variables,VarName )
        IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
          IF ( .NOT.PRESENT(ONOEfound) ) CYCLE
        ELSE
          IF ( PRESENT(ONOEfound) ) CYCLE
        END IF

	VarName(1:1) = CHAR(ICHAR(VarName(1:1))-ICHAR('a')+ICHAR('A'))
        k = LEN_TRIM(VarName)
        DO j=1,k
          IF ( VarName(j:j) == ' ' ) VarName(j:j) = '.'
        END DO
       
        WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' scalar: '//TRIM(VarName)
      END DO
      
      DO Vari = 1, VectorFields
        WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
        Var => VariableGet( Solver % Mesh % Variables,VarName )
        IF(.NOT.ASSOCIATED(Var)) & 
          Var => VariableGet( Solver % Mesh % Variables, TRIM(Varname)//' 1' )
        IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
          IF ( .NOT.PRESENT(ONOEfound) ) CYCLE
        ELSE
          IF ( PRESENT(ONOEfound) ) CYCLE
        END IF


	VarName(1:1) = CHAR(ICHAR(VarName(1:1))-ICHAR('a')+ICHAR('A'))
        k = LEN_TRIM(VarName)
        DO j=1,k
          IF ( VarName(j:j) == ' ' ) VarName(j:j) = '.'
        END DO

        WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: '//TRIM(VarName)
      END DO
      
      WRITE(PostFileUnit,'()')
      DateStr = FormatDate()
      WRITE( PostFileUnit, '("#File started at: ",A)' ) TRIM(DateStr)
!------------------------------------------------------------------------------
!   Coordinates
!------------------------------------------------------------------------------
!
      MeshScale = 1.0_dp
      DO i=1,Model % NumberOfSolvers
        IF ( Model % Solvers(i) % Variable % NameLen <= 0 ) CYCLE
        
        IF (Model % Solvers(i) % Variable % Name &
            (1:Model % Solvers(i) % Variable % NameLen) == 'displacement') THEN
          MeshMoved = ListGetLogical( Model % Solvers(i) % Values, &
              'Displace Mesh', Gotit )

          IF ( .NOT. GotIt ) &
            MeshMoved=.NOT.EigenOrHarmonicAnalysis(Model % Solvers(i))

          IF ( .NOT.MeshMoved ) MeshScale = 0.0_dp
        END IF
      END DO
      
      IF ( PRESENT(ONOEfound) ) THEN
        DO i=1,Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(i)
          IF(MaskExists) THEN
            IF( ANY(MaskPerm(Element % NodeIndexes)<=0)) CYCLE
          END IF

          DO j=1,Element % TYPE % NumberOfNodes
            NodeCoords(1) = Model % Nodes % x(Element % NodeIndexes(j))
            NodeCoords(2) = Model % Nodes % y(Element % NodeIndexes(j))
            NodeCoords(3) = Model % Nodes % z(Element % NodeIndexes(j))

            IF( MeshDim <= 2 ) THEN           
              WRITE(PostFileUnit,'(2ES16.7E3,A)') NodeCoords(1:2), ' 0.0'
            ELSE
              WRITE(PostFileUnit,'(3ES16.7E3)') NodeCoords(1:3)
            END IF
          END DO
        END DO
      ELSE
        DO ii=1,NumberOfNodes
          i = ii
          IF(MaskExists) i = MaskOrder(i)

          NodeCoords(1) = Model % Nodes % x(i)
          NodeCoords(2) = Model % Nodes % y(i)
          NodeCoords(3) = Model % Nodes % z(i)
        
          IF ( ASSOCIATED(Displacement) ) THEN
            k = Displacement % Perm(i)
          
            IF ( k > 0 ) THEN
              k = Displacement % DOFs * (k-1)
              NodeCoords(1) = NodeCoords(1) - MeshScale*Displacement % Values(k+1)
              IF( Displacement % DOFs >= 2) THEN
                NodeCoords(2) = NodeCoords(2) - MeshScale*Displacement % Values(k+2)
                END IF
              IF( Displacement % DOFs == 3) THEN
                NodeCoords(3) = NodeCoords(3) - MeshScale*Displacement % Values(k+3)
              END IF
            ELSE
              IF ( ASSOCIATED( MeshUpdate ) ) k = MeshUpdate % Perm(i)
            
              k = MeshUpdate % DOFs * (k-1)
              NodeCoords(1) = NodeCoords(1) - MeshUpdate % Values(k+1)
              IF( MeshUpdate % DOFs >= 2) THEN
                NodeCoords(2) = NodeCoords(2) - MeshUpdate % Values(k+2)
              END IF
              IF( MeshUpdate % DOFs == 3) THEN
                NodeCoords(3) = NodeCoords(3) - MeshUpdate % Values(k+3)
              END IF
            END IF
          END IF

          IF( MeshDim <= 2 ) THEN           
            WRITE(PostFileUnit,'(2ES16.7E3,A)') NodeCoords(1:2), ' 0.0'
          ELSE
            WRITE(PostFileUnit,'(3ES16.7E3)') NodeCoords(1:3)
          END IF
        END DO
      END IF

!------------------------------------------------------------------------------
! Elements
!------------------------------------------------------------------------------
      WRITE(PostFileUnit,'(a)') '#group all'
      ONOECount=0
      DO i=1,Model % NumberOfBulkElements
        Element => Model % Elements(i)
        
        IF(MaskExists) THEN
          IF( ANY(MaskPerm(Element % NodeIndexes)<=0)) CYCLE
        END IF
        
        k = Element % BodyId
        gotIt = .FALSE.
        IF ( k >= 1 .AND. k <= Model % NumberOfBodies ) THEN
          Str = ListGetString( Model % Bodies(k) % Values,'Name',gotIt )
        END IF
        
        IF ( gotIt ) THEN
          k = LEN_TRIM(Str)
          DO j=1,k
            IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
          END DO
          
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
        ELSE
          WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) 'body'//TRIM(I2S(k))//' '
        END IF
        
        WRITE(PostFileUnit,'(i5)', ADVANCE='NO') Element % TYPE % ElementCode
        n = 0
        DO j=1,Element % TYPE % NumberOfNodes,4
          DO k=1,MIN(4,Element % TYPE % NumberOfNodes-n)
            n = n + 1
            ind = Element % NodeIndexes(n)
            IF(PRESENT(ONOEfound)) THEN
              ONOECount = ONOECount+1
              ind = ONOECount
            ELSE IF(MaskExists) THEN
              ind = MaskPerm(ind)
            END IF
            WRITE(PostFileUnit, '(i8)', ADVANCE='NO')  ind - 1
          END DO
          WRITE( PostFileUnit,'(a)' ) ''
        END DO
      END DO

      DO i=Model % NumberOfBulkElements + 1,Model % NumberOfBulkElements + &
          Model % NumberOfBoundaryElements
        
        Element => Model % Elements(i)
        
        IF(MaskExists) THEN
          IF( ANY(MaskPerm(Element % NodeIndexes) <= 0)) CYCLE
        END IF
        
        k = Element % BoundaryInfo % Constraint
        
        gotIt = .FALSE.
        IF ( k >= 1 .AND. k <= Model % NumberOfBCs ) THEN
          Str = ListGetString( Model % BCs(k) % Values,'Name',gotIt )
        END IF
        
        IF ( gotIt ) THEN
          k = LEN_TRIM(Str)
          DO j=1,k
            IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
          END DO
          
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
        ELSE
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' ) 'Constraint'//TRIM(I2S(k))//' '
        END IF
        
        WRITE(PostFileUnit,'(i5)', ADVANCE='NO') Element % TYPE % ElementCode
        DO k=1,Element % TYPE % NumberOfNodes
          ind = Element % NodeIndexes(k)
          IF(PRESENT(ONOEfound)) THEN
            ONOECount=ONOECount+1
            ind = ONOECount
          ELSE IF(MaskExists) THEN
            ind = MaskPerm(ind)
          END IF
          WRITE( PostFileUnit, '(i8)', ADVANCE='NO' )  ind-1
        END DO
        WRITE( PostFileUnit,'(a)' ) ''
      END DO
      WRITE(PostFileUnit,'(a)') '#endgroup all'
    END IF
   
!------------------------------------------------------------------------------
!  Save resulst on a timestep (or steady state iteration step)
!------------------------------------------------------------------------------
 
    TimeStep   = SavedCount
    Var => VariableGet( Model % Variables, 'Time' )        
    Time = 1.0d0
    IF ( ASSOCIATED(Var) ) Time = Var % Values(1)

    
    WRITE( PostFileUnit,'(a,i7,i7,ES16.7E3)' ) '#time ',SavedCount,Timestep,Time

    IF ( PRESENT(ONOEfound) ) THEN
       n1=Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
    ELSE
       n1=1
       n2=NumberOfNodes
    END IF

    DO jj=1,n1

      IF( PRESENT(ONOEfound) ) THEN
        Element => Model % Elements(jj)
        IF(MaskExists) THEN
          IF ( ANY(MaskPerm(Element % NodeIndexes)<=0) ) CYCLE
        END IF
        n2=Element % TYPE % NumberOfNodes

        Parent => NULL()
        IF (ASSOCIATED(Element % BoundaryInfo) ) THEN
          Parent => Element % BoundaryInfo % Left
          n = 0
          IF(ASSOCIATED(Parent)) THEN
            DO l=1,GetElementNOFNodes(Parent)
              DO m=1,n2
                IF (Element % NodeIndexes(m)==Parent % NodeIndexes(l)) n=n+1
              END DO
            END DO
          END IF
          IF (n/=n2) Parent => Element % BoundaryInfo % Right
        END IF
      END IF

      DO ii=1,n2
      
        IF ( PRESENT(ONOEfound) ) THEN
          IF ( ASSOCIATED(Parent) ) THEN
            DO l=1,GetElementNOFNOdes(Parent)
              IF (Element % NodeIndexes(ii)==Parent % NodeIndexes(l)) EXIT
            END DO
            i = Parent % DGIndexes(l)
          ELSE
            i = Element % DGIndexes(ii)
          END IF
        ELSE
          i = ii
          IF(MaskExists) i = MaskOrder(i)
        END IF
      
        DO Vari = 1, ScalarFields
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
          Var => VariableGet( Model % Variables, TRIM(VarName) ) 
          IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
            IF (.NOT.PRESENT(ONOEfound)) CYCLE
          ELSE IF ( PRESENT(ONOEfound) ) THEN
            CYCLE
          END IF
          Values => Var % Values

          k = i
          IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
          IF ( k > 0 ) THEN
            WRITE(PostFileUnit,'(ES16.7E3)',ADVANCE='NO') Values(k)
          ELSE
            WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
          END IF
        END DO

        DO Vari = 1, VectorFields

          WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
          Var => VariableGet( Model % Variables, VarName ) 

          IF( ASSOCIATED(Var) ) THEN
            IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
              IF (.NOT.PRESENT(ONOEfound) ) CYCLE
            ELSE IF ( PRESENT(ONOEfound) ) THEN
              CYCLE
            END IF
            k = i
            Dofs = Var% Dofs
            IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
          
            IF( k == 0 ) THEN
              ! Check for the presence of secondary variable (in practice 'mesh update')
              WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
              VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
              IF( GotIt ) THEN
                Var => VariableGet( Model % Variables, VarName ) 
                k = Var % Perm(i)
              END IF
            END IF
          
            IF( k > 0 ) THEN
              DO j=1,dofs
                WRITE(PostFileUnit,'(ES16.7E3)',ADVANCE='NO') Var % Values(dofs*(k-1)+j)
              END DO
              IF( dofs < 3 ) THEN
                WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
              END IF
            ELSE 
              WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
            END IF
          
          ELSE
            ! Check for the presence of component vectors given by its components i=1,2,3
            Var => VariableGet( Model % Variables, TRIM(VarName)//' 1' ) 
            IF( ASSOCIATED(Var)) THEN
              IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
                IF (.NOT.PRESENT(ONOEfound) ) CYCLE
              ELSE IF ( PRESENT(ONOEfound) ) THEN
                CYCLE
              END IF

              k = i
              IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)

              IF( k == 0 ) THEN
                WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
                VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
                IF( GotIt ) THEN
                  Var => VariableGet( Model % Variables, TRIM(VarName)//' 1' ) 
                  k = Var % Perm(i)
                END IF
              END IF
              IF( k > 0 ) THEN
                Values => Var % Values
                Var => VariableGet( Model % Variables, TRIM(VarName)//' 2' ) 
                Values2 => Var % Values
                IF( MeshDim == 2 ) THEN
                  WRITE(PostFileUnit,'(2ES16.7E3,A)',ADVANCE='NO') Values(k),&
                      Values2(k),' 0.0'               
                ELSE
                  Var => VariableGet( Model % Variables, TRIM(VarName)//' 3' ) 
                  Values3 => Var % Values
                  WRITE(PostFileUnit,'(3ES16.7E3)',ADVANCE='NO') Values(k),&
                      Values2(k), Values3(k)
                END IF
              ELSE
                WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'              
              END IF
            END IF
          END IF
        END DO
        WRITE(PostFileUnit,'()')
      END DO
    END DO

!------------------------------------------------------------------------------
!   We are done here close the files and deallocate
!------------------------------------------------------------------------------
10  CONTINUE

    CLOSE(PostFileUnit)
    
    IF(MaskExists) DEALLOCATE(MaskOrder)
    IF(MaskAllocated) DEALLOCATE( MaskPerm )

    IF ( .NOT. PRESENT(ONOEfound) .AND. On_nodes_on_elements>0 ) &
      CALL ElmerPostOutputSolver(Model,Solver,dt,TransientSimulation,.TRUE.)

  END SUBROUTINE ElmerPostOutputSolver
!------------------------------------------------------------------------------

