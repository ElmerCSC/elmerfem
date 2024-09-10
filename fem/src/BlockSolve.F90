!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/


!> \ingroup ElmerLib
!> \{
!-----------------------------------------------------------------------------
!> Module for solution of matrix equation utilizing block strategies.
!-----------------------------------------------------------------------------

MODULE BlockSolve

 USE DefUtils

 IMPLICIT NONE

  TYPE(BlockMatrix_t), POINTER, SAVE :: TotMatrix

  LOGICAL, PRIVATE :: isParallel=.FALSE.

  TYPE(Solver_t), POINTER, SAVE :: SolverRef
  TYPE(Variable_t), POINTER :: SolverVar => Null()
  TYPE(Matrix_t), POINTER :: SolverMatrix => Null(), SaveMatrix

CONTAINS

#define USEPERM 0
  
  !-----------------------------------------------------------------------------------
  !> If a block variable does not exist it will be created. 
  !> Here only normal nodal elements are supported for the moment. 
  !> Then also the creation of permutation vector is straight-forward.
  !> Note that no reordering is currently performed.
  !
  !> There is limitation regarding non-nodal elements which stems partly from the fact
  !> that an elementtype is solver specific while this one solver could have a number of
  !> different elementtypes for different equations.
  !-----------------------------------------------------------------------------------
  FUNCTION CreateBlockVariable( Solver, VariableNo, VarName, ExtDofs, ExtPerm ) RESULT ( Var )
    
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: VariableNo
    CHARACTER(LEN=*) :: VarName
    INTEGER, OPTIONAL :: ExtDofs
    INTEGER, POINTER, OPTIONAL :: ExtPerm(:)
    TYPE(Variable_t), POINTER :: Var
    
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,t,n,ndeg,body_id,body_id_prev,eq_id,solver_id,Dofs,nsize
    INTEGER, POINTER :: VarPerm(:), ActiveVariables(:)
    INTEGER, ALLOCATABLE :: Indexes(:)
    REAL(KIND=dp), POINTER :: Values(:)
    LOGICAL :: Hit, GotIt
    CHARACTER(:), ALLOCATABLE :: str
    
    LOGICAL :: GlobalBubbles, Found
    INTEGER :: MaxNDOFs, MaxDGDOFs, MaxEDOFs, MaxFDOFs, MaxBDOFs


    CALL Info('CreateBlockVariable','Creating block variables',Level=8)
    
    Mesh => Solver % Mesh
    Params => Solver % Values

    IF( PRESENT( ExtDofs ) ) THEN
      Dofs = ExtDofs
    ELSE 
      str = 'Variable '//I2S(VariableNo)//' Dofs'
      Dofs = ListGetInteger( Params,str, GotIt )
      IF(.NOT. GotIt) Dofs = 1
    END IF
    
    IF( PRESENT( ExtPerm ) ) THEN
      nsize = MAXVAL( ExtPerm ) 
      varPerm => ExtPerm
      GOTO 100
    END IF
    
    Ndeg = 0
    MaxNDOFs  = 0
    MaxBDOFs = 0
    MaxDGDOFs = 0
    DO i=1, Mesh % NumberOFBulkElements
      Element => Mesh % Elements(i)
      MaxNDOFs  = MAX( MaxNDOFs,  Element % NDOFs )
      MaxBDOFs  = MAX( MaxBDOFs,  Element % BDOFs )
      MaxDGDOFs = MAX( MaxDGDOFs, Element % DGDOFs )
    END DO
    
    MaxEDOFs = 0
    DO i=1, Mesh % NumberOFEdges
      Element => Solver % Mesh % Edges(i)
      MaxEDOFs  = MAX( MaxEDOFs,  Element % BDOFs )
    END DO
    
    MaxFDOFs = 0
    DO i=1, Mesh % NumberOFFaces
      Element => Solver % Mesh % Faces(i)
      MaxFDOFs  = MAX( MaxFDOFs,  Element % BDOFs )
    END DO
    
    GlobalBubbles = ListGetLogical( Params, 'Bubbles in Global System', Found )
    IF (.NOT.Found) GlobalBubbles = .TRUE.
    
    Ndeg = Ndeg + Mesh % NumberOfNodes
    IF ( MaxEDOFs > 0 ) Ndeg = Ndeg + MaxEDOFs * Mesh % NumberOFEdges
    IF ( MaxFDOFs > 0 ) Ndeg = Ndeg + MaxFDOFs * Mesh % NumberOFFaces
    IF ( GlobalBubbles ) &
        Ndeg = Ndeg + MaxBDOFs * Mesh % NumberOfBulkElements
    IF ( ListGetLogical( Params, 'Discontinuous Galerkin', Found ) ) &
        Ndeg = MAX( NDeg, MaxDGDOFs * ( Mesh % NumberOfBulkElements + &
        Mesh % NumberOfBoundaryElements) )
    
    ALLOCATE( VarPerm(ndeg) )
    VarPerm = 0
    
    solver_id = 0
    DO i = 1, CurrentModel % NumberOfSolvers
      PSolver => CurrentModel % Solvers(i)
      IF( ASSOCIATED( Solver, PSolver ) ) THEN
        solver_id = i
        EXIT
      END IF
    END DO
    
    ALLOCATE(Indexes(Mesh % MaxElementDOFs))
    body_id_prev = -1
    DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
      Element => Mesh % Elements(t)
      CurrentModel % CurrentElement => Element
      
      body_id = Element % BodyId 
      IF( body_id /= body_id_prev ) THEN
        Hit = .FALSE.
        body_id_prev = body_id
        
        IF( body_id < 1 ) CYCLE
        eq_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values,'Equation')
        IF( eq_id < 1 ) CYCLE
        
        str='Active Variables['//i2s(solver_id)//']'            
        ActiveVariables => ListGetIntegerArray(CurrentModel % Equations(eq_id) % Values, str)
        IF(.NOT. ASSOCIATED(ActiveVariables)) THEN
          ActiveVariables => ListGetIntegerArray( CurrentModel % Equations(eq_id) % Values, &
              'Active Variables')
        END IF
        IF( .NOT. ASSOCIATED(ActiveVariables)) CYCLE
        IF( .NOT. ANY(ActiveVariables == VariableNo) )  CYCLE
        Hit = .TRUE.
      END IF
      
      IF( Hit ) THEN
         n=GetElementDOFs(Indexes)
         VarPerm(Indexes(1:n)) = 1
      END IF
    END DO
    
    j = 0
    DO i = 1, SIZE(VarPerm)
      IF( VarPerm(i) > 0 ) THEN
        j = j + 1
        VarPerm(i) = j
      END IF
    END DO
    nsize = j
    
100 IF( nsize == 0 ) THEN
      CALL Fatal('CreateBlockVariable','Variable '//TRIM(VarName)//' cannot be created')      
    END IF

    CALL Info('CreateBlockVariable','Creating variable: '//TRIM(VarName), Level=6 )
    
    CALL VariableAddVector( Mesh % Variables, Mesh, Solver, &
        TRIM(VarName), Dofs, Perm = VarPerm )          
    
    WRITE( Message,'(A,I0,A)') 'Creating variable '//TRIM(VarName)//' with ',Dofs,' dofs'
    CALL Info('CreateBlockVariable', Message)
    
    Var => VariableGet( Mesh % Variables, VarName )         

  END FUNCTION CreateBlockVariable


  !-------------------------------------------------------------------
  !> This subroutine initializes the block matrix structure so that the 
  !> matrices and vectors have a natural location to save.
  !------------------------------------------------------------------
  SUBROUTINE BlockInitMatrix( Solver, BlockMatrix, BlockDofs, FieldDofs, SkipVar )
    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    INTEGER :: BlockDofs
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix
    INTEGER, OPTIONAL :: FieldDofs
    LOGICAL, OPTIONAL :: SkipVar
    
    TYPE(Solver_t), POINTER :: PSolver
    INTEGER, POINTER :: BlockStruct(:), SlaveSolvers(:)
    LOGICAL :: GotBlockStruct, GotSlaveSolvers, Found
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER :: i,j,k,n,Novar
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(:), ALLOCATABLE :: VarName, str
    LOGICAL :: UseSolverMatrix, IsComplex
    CHARACTER(*), PARAMETER :: Caller = 'BlockInitMatrix'

            
    Params => Solver % Values

    IsComplex = ListGetLogical( Params,'Linear System Complex',Found)
    IF( IsComplex ) THEN
      CALL Info(Caller,'Assuming block matrix to be complex!',Level=8)
    ELSE
      CALL Info(Caller,'Assuming block matrix to be real!',Level=20)
    END IF
    
    BlockMatrix => Solver % BlockMatrix
    IF (ASSOCIATED(BlockMatrix)) THEN
      CALL Info(Caller,'Using existing block matrix',Level=10)
      RETURN
    END IF

    CALL Info(Caller,'Initializing block matrix',Level=10)
    
    ALLOCATE(Solver % BlockMatrix)
    BlockMatrix => Solver % BlockMatrix
 
    BlockStruct => ListGetIntegerArray( Params,'Block Structure',GotBlockStruct)
    BlockMatrix % GotBlockStruct = GotBlockStruct
    
    IF( GotBlockStruct ) THEN
      IF( SIZE( BlockStruct ) /= BlockDofs ) THEN
        CALL Fatal(Caller,'Incompatible size of > Block Structure < given!')
      END IF
      IF( MINVAL( BlockStruct ) < 1 .OR. MAXVAL( BlockStruct ) > BlockDofs ) THEN
        CALL Fatal(Caller,'Incompatible values in > Block Structure < given!')          
      END IF
      NoVar = MAXVAL( BlockStruct )
      CALL Info(Caller,'Using given block structure of size: '//I2S(SIZE( BlockStruct)),Level=8)
      BlockMatrix % BlockStruct => BlockStruct

      ALLOCATE( BlockMatrix % InvBlockStruct(NoVar) )
      BlockMatrix % InvBlockStruct = 0
      DO i=1,BlockDofs
        j = BlockStruct(i)
        IF( BlockMatrix % InvBlockStruct(j) == 0 ) THEN
          BlockMatrix % InvBlockStruct(j) = i
        ELSE
          ! Block structure is not bijection for this component
          BlockMatrix % InvBlockStruct(j) = -1
        END IF
      END DO        
    ELSE
      CALL Info(Caller,'Inheriting blocks from variable dofs',Level=8)
      NoVar = BlockDofs
    END IF    


    IF( BlockMatrix % NoVar == NoVar ) THEN
      CALL Info(Caller,'Reusing existing blockmatrix',Level=6)
      RETURN
    ELSE IF( BlockMatrix % Novar /= 0 ) THEN
      CALL Fatal(Caller,'Previous blockmatrix was of different size?')
    ELSE
      CALL Info(Caller,'Block matrix will be of size '//I2S(NoVar),Level=6)
    END IF
    
    BlockMatrix % Solver => Solver
    BlockMatrix % NoVar = NoVar

    
    ALLOCATE( BlockMatrix % SubMatrix(NoVar,NoVar) )
    DO i=1,NoVar
      DO j=1,NoVar
        DO k=1,2
          Amat => NULL()
          Amat => AllocateMatrix()
          Amat % ListMatrix => NULL()
          Amat % FORMAT = MATRIX_LIST      
          Amat % NumberOfRows = 0
          AMat % COMPLEX = IsComplex
          IF( k==1) THEN
            ! First cycle creates the holder for standard matrix
            BlockMatrix % Submatrix(i,j) % Mat => Amat
          ELSE
            ! Second cycle creates the optional holder for preconditioning matrix
            BlockMatrix % Submatrix(i,j) % PrecMat => Amat
          END IF
        END DO
      END DO
    END DO
    
    ALLOCATE( BlockMatrix % SubMatrixActive(NoVar,NoVar) )
    BlockMatrix % SubMatrixActive = .FALSE.

    ALLOCATE( BlockMatrix % SubMatrixTranspose(NoVar,NoVar) )
    BlockMatrix % SubMatrixTranspose = .FALSE.

    ALLOCATE( BlockMatrix % SubVector(NoVar))
    DO i=1,NoVar
      BlockMatrix % Subvector(i) % Var => NULL()
    END DO

    ALLOCATE( BlockMatrix % Offset(NoVar+1))
    BlockMatrix % Offset = 0
    BlockMatrix % maxsize = 0

    
    IF( PRESENT( SkipVar ) ) THEN
      IF( SkipVar ) THEN
        BlockMatrix % ParentMatrix => Solver % Matrix
        RETURN
      END IF
    END IF


    ! We may have different size of block matrix than the number of actual components.
    ! For example, when we have a projector of a scalar field our block size is (2,2)
    ! but we can only create the (1,1) from the initial matrix system. 
    IF( PRESENT( FieldDofs ) ) THEN
      CALL Info(Caller,'Number of field components: '//I2S(FieldDofs))
      IF( Novar /= FieldDofs ) CALL Info(Caller,'Number of fields and blocks ('&
          //I2S(NoVar)//') differ!')
      IF(.NOT. GotBlockStruct ) NoVar = FieldDofs
    END IF

    
    ! If we have just one variable and also one matrix then no need to look further
    ! This would probably just happen for testing purposes. 
    UseSolverMatrix = (NoVar == 1 )
    IF( UseSolverMatrix ) THEN
      UseSolverMatrix = ASSOCIATED( Solver % Variable )
    END IF
    IF( UseSolverMatrix ) THEN
      UseSolverMatrix = ASSOCIATED( Solver % Matrix )
    END IF        
    IF( UseSolverMatrix ) THEN
      BlockMatrix % SubVector(1) % Var => Solver % Variable
      n = SIZE( Solver % Variable % Values )
      BlockMatrix % Offset(1) = 0
      BlockMatrix % Offset(2) = n      
      BlockMatrix % MaxSize = n
      BlockMatrix % TotSize = n
      RETURN
    END IF

    SlaveSolvers =>  ListGetIntegerArray( Params, &
        'Block Solvers', GotSlaveSolvers )

   
    DO i = 1,NoVar
      IF( GotSlaveSolvers ) THEN
        j = SlaveSolvers(i)

        CALL Info(Caller,'Associating block '//I2S(i)//' with solver: '//I2S(j),Level=10)

        PSolver => CurrentModel % Solvers(j)
        Var => PSolver % Variable 
        VarName = TRIM(Var % Name)

        BlockMatrix % SubVector(i) % Solver => PSolver
        BlockMatrix % SubMatrix(i,i) % Mat => PSolver % Matrix        

      ELSE
        str = 'Variable '//I2S(i)

        VarName = ListGetString( Params, str, Found )
        IF(.NOT. Found ) THEN       
          IF( BlockMatrix % GotBlockStruct ) THEN
           VarName = 'BlockVar '//I2S(i)
          ELSE
            VarName = ComponentName(Solver % Variable % Name,i)            
          END IF
        END IF
        Var => VariableGet( Solver % Mesh % Variables, VarName )
      END IF
      
      !-----------------------------------------------------------------------------------
      ! If variable does not exist it will be created. 
      ! Here it is assumed that all components have the same number of dofs
      ! described by the same permutation vector. If the components are
      ! accounted in normal manner [1,2,3,...] then it suffices just to have 
      ! pointers to the components of the full vector.
      !-----------------------------------------------------------------------------------
      IF(ASSOCIATED( Var ) ) THEN
        CALL Info(Caller,'Using existing variable > '//VarName//' <')		
      ELSE		
        CALL Info(Caller,'Variable > '//VarName//' < does not exist, creating')
        PSolver => Solver
        IF( BlockMatrix % GotBlockStruct ) THEN
          j = COUNT( BlockMatrix % BlockStruct == i ) 
          IF( j == 0 ) THEN
            CALL Fatal(Caller,'Invalid > Block Structure < given!')
          END IF
          Var => CreateBlockVariable(PSolver, i, VarName, j, Solver % Variable % Perm )
        ELSE
          Var => CreateBlockVariable(PSolver, i, VarName )
        END IF
      END IF
      
      BlockMatrix % SubVector(i) % Var => Var
      n = SIZE(Var % Values)
      
      BlockMatrix % Offset(i+1) = BlockMatrix % Offset(i) + n
      BlockMatrix % MaxSize = MAX( BlockMatrix % MaxSize, n )
    END DO
    
    BlockMatrix % TotSize = BlockMatrix % Offset( NoVar + 1 )

    CALL Info(Caller,'All done',Level=12)
      
  END SUBROUTINE BlockInitMatrix
    


  !-------------------------------------------------------------------
  !> This subroutine creates the missing component variables.
  !------------------------------------------------------------------
  SUBROUTINE BlockInitVar( Solver, BlockMatrix )
    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix
    
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER :: i,j,k,n,Novar
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(:), ALLOCATABLE :: VarName, str
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: Vals(:)
    
    Params => Solver % Values
    Mesh => Solver % Mesh
    NoVar = BlockMatrix % NoVar
    
    DO i=1,NoVar
      Amat => BlockMatrix % Submatrix(i,i) % Mat 
      n = Amat % NumberOfRows
      
      BlockMatrix % Offset(i+1) = BlockMatrix % Offset(i) + n
      BlockMatrix % MaxSize = MAX( BlockMatrix % MaxSize, n )
      
      VarName = ComponentName("Block variable",i)            
      Var => VariableGet( Mesh % Variables, VarName )
      IF(.NOT. ASSOCIATED( Var ) ) THEN
        CALL Info('BlockInitVar','Variable > '//VarName//' < does not exist, creating')
        PSolver => Solver
        NULLIFY( Vals )
        ALLOCATE( Vals(n) )
        Vals = 0.0_dp
        
        CALL VariableAdd( Mesh % Variables,Mesh,PSolver,VarName,1,Vals,&
            Output = .FALSE. )
        !Perm,Output,Secondary, TYPE )
        Var => VariableGet( Mesh % Variables, VarName )
      END IF
      BlockMatrix % SubVector(i) % Var => Var

      ! Take the monolithic solution as initial guess
      !DO j=1,n
      !  k = Amat % InvPerm(j)
      !  Var % Values(j) = Solver % Variable % Values(k)
      !END DO

    END DO
        
    BlockMatrix % TotSize = BlockMatrix % Offset( NoVar + 1 )

    CALL Info('BlockInitVar','All done',Level=12)
      
  END SUBROUTINE BlockInitVar


  !-------------------------------------------------------------------
  !> This subroutine copies back the full vector from its components.
  !------------------------------------------------------------------
  SUBROUTINE BlockBackCopyVar( Solver, BlockMatrix )
    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix
    
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER :: i,j,k,n,m,Novar
    TYPE(Variable_t), POINTER :: Var
    
    CALL Info('BlockBackCopyVar','Copying values back to monolithic solution vector',Level=10)

    NoVar = BlockMatrix % NoVar

    m = SIZE( Solver % Variable % Values ) 
   
    DO i=1,NoVar
      Amat => BlockMatrix % Submatrix(i,i) % Mat 
      n = Amat % NumberOfRows
      Var => BlockMatrix % SubVector(i) % Var 
      
      ! Copy the block part to the monolithic solution
      DO j=1,n
        k = Amat % InvPerm(j)
        IF( k < 1 .OR. k > m ) THEN
          PRINT *,'ijk:',i,j,k
          CYCLE
        END IF
        Solver % Variable % Values(k) = Var % Values(j)
      END DO

    END DO
        
    BlockMatrix % TotSize = BlockMatrix % Offset( NoVar + 1 )

    CALL Info('BlockBackCopyVar','All done',Level=15)
      
  END SUBROUTINE BlockBackCopyVar

  
  
  !-------------------------------------------------------------------------------------
  !> Picks the components of a full matrix to the submatrices of a block matrix.
  !> On choice, the user may have exactly the same block matrix structure than 
  !> a leading component.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrix( Solver, NoVar )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar

    INTEGER :: RowVar, ColVar
    TYPE(Matrix_t), POINTER :: SolverMatrix
    TYPE(Matrix_t), POINTER :: Amat
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ReuseMatrix, Found, EliminateZero
    INTEGER::i,j,k,l,n
    REAL(KIND=DP) :: SumAbsMat
    
    CALL Info('BlockPickMatrix','Picking block matrix of size '//I2S(NoVar)//' from monolithic one',Level=10)

    SolverMatrix => Solver % Matrix 
    Params => Solver % Values
        
    ReuseMatrix = ListGetLogical( Params,'Block Matrix Reuse',Found)
    EliminateZero = ListGetLogical( Params, &
                         'Block Eliminate Zero Submatrices', Found )

    DO RowVar=1,NoVar
      DO ColVar=1,NoVar            
        Amat => TotMatrix % Submatrix(RowVar,ColVar) % Mat          
        IF( TotMatrix % GotBlockStruct) THEN
          ! A generic picking method for submatrices
          !----------------------------------------------------------------------
          CALL Info('BlockPickMatrix','Picking generic block matrix ('&
              //I2S(RowVar)//','//I2S(ColVar)//')',Level=20)
          CALL CRS_BlockMatrixPick2(SolverMatrix,Amat,TotMatrix % BlockStruct,RowVar,ColVar)
        ELSE
          ! Picking of standard submatrices of size one.
          !----------------------------------------------------------------------
          IF( ReuseMatrix ) THEN
            IF( RowVar + ColVar > 2 .AND. Amat % NumberOfRows == 0 ) THEN
              CALL Info('BlockPickMatrix','Copying block matrix topology ('&
                  //I2S(RowVar)//','//I2S(ColVar)//')',Level=20)
              CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(1,1) % Mat, Amat )
            END IF
          END IF
          CALL Info('BlockPickMatrix','Picking simple block matrix ('&
              //I2S(RowVar)//','//I2S(ColVar)//')',Level=20)          
          CALL CRS_BlockMatrixPick(SolverMatrix,Amat,NoVar,RowVar,ColVar,RowVar == ColVar )          

          IF( EliminateZero ) THEN
            IF( Amat % NumberOfRows > 0 ) THEN
              SumAbsMat = SUM( ABS( Amat % Values ) )
              IF( SumAbsMat < SQRT( TINY( SumAbsMat ) ) ) THEN
                CALL Info('BlockPickMatrix','Matrix is actually all zero, eliminating it!',Level=12)
                DEALLOCATE( Amat % Values ) 
                IF( .NOT. ReuseMatrix ) THEN
                  DEALLOCATE( Amat % Rows, Amat % Cols )
                  IF( RowVar == ColVar ) DEALLOCATE( Amat % Diag, Amat % rhs ) 
                END IF
                Amat % NumberOfRows = 0
              END IF
            END IF
          END IF

        END IF

!        CALL CRS_SortMatrix( Amat, .TRUE. )        
      END DO
    END DO

    BLOCK
      INTEGER, POINTER :: BlockStruct(:),BlockPerm(:)
      INTEGER ::  nl,nk,nv,n0,nj
      
      CALL Info('BlockPickMatrix','Creating permutation to map between block and mono vectors',Level=12)
      
      n = SolverMatrix % NumberOfRows
      IF(.NOT. ASSOCIATED( TotMatrix % BlockPerm ) ) THEN
        ALLOCATE( TotMatrix % BlockPerm(n) )
      END IF
      BlockPerm => TotMatrix % BlockPerm
      BlockPerm = 0
      
      IF( TotMatrix % GotBlockStruct ) THEN
        BlockStruct => TotMatrix % BlockStruct

        ! This is permutation for nontrivial block structure, for example (1 1 1 2)
        nl = SIZE(BlockStruct)   ! number of original components
        nk = MAXVAL(BlockStruct) ! number of blocks
        nv = n / nl              ! size of each field

        DO k=1,nk
          n0 = COUNT(BlockStruct < k )   ! number of components prior to this block
          nj = COUNT(BlockStruct == k )  ! number of components in this block
          j = 0
          DO l=1,nl
            IF(BlockStruct(l) /= k) CYCLE
            j = j+1            
            DO i=1,nv
              BlockPerm(nv*n0+nj*(i-1)+j) = nl*(i-1)+l
            END DO
          END DO
        END DO        
      ELSE
        ! This is trivial numbering for default block structure (1 2 3 4 ...) 
        nv = n/NoVar
        DO j=1,NoVar
          DO i=1,nv
            BlockPerm((j-1)*nv+i) = Novar*(i-1) + j 
          END DO
        END DO
      END IF
    END BLOCK

  END SUBROUTINE BlockPickMatrix


  !-------------------------------------------------------------------------------------
  !> Picks the components of a full matrix to given domains or bodies.
  !> The rest stays in 1st domain.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickDofsPhysical( Solver, BlockIndex, NoVar )
    
    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER :: BlockIndex(:)
    INTEGER :: Novar
    
    INTEGER::i,j,k,t,n,MinBlock,MaxBlock,body_id,bf_id,bc_id,n_bf
    TYPE(ValueList_t), POINTER :: List
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Found
    INTEGER :: ElemPerm(27)
    INTEGER, POINTER :: Perm(:)
    
    
    CALL Info('BlockPickDofsPhysical','Picking block matrix of size '//I2S(NoVar)//' from monolithic one',Level=10)

    n_bf = CurrentModel % NumberOfBodyForces

    MinBlock = HUGE(MinBlock)
    MaxBlock = 0
    DO i=1,n_bf + CurrentModel % NumberOfBCs
      IF( i <= n_bf ) THEN
        List => CurrentModel % BodyForces(i) % Values
      ELSE
        List => CurrentModel % BCs(i-n_bf) % Values
      END IF        
      j = ListGetInteger( List,'Block Index',Found )
      IF( Found ) THEN
        MinBlock = MIN(j,MinBlock)
        MaxBlock = MAX(j,MaxBlock)      
      END IF
    END DO
    
    IF( MaxBlock == 0 ) THEN
      CALL Fatal('BlockPickDofsPhysical','Cannot create a physical block structure as no >Block Index< given!')
    END IF

    Mesh => Solver % Mesh 
    Perm => Solver % Variable % Perm 
    n = MAXVAL( Perm ) 
    BlockIndex = 0
        
    DO t=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      IF( t <= Mesh % NumberOfBulkElements ) THEN
        body_id = Element % BodyId
        bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values,'Body Force', Found )
        IF( bf_id == 0 ) CYCLE
        List => CurrentModel % BodyForces(bf_id) % Values
      ELSE
        IF(.NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE             
        DO bc_id=1,CurrentModel % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
        END DO               
        IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE        
        List => CurrentModel % BCs(bc_id) % Values
      END IF
      
      j = ListGetInteger( List,'Block Index',Found )
      IF( .NOT. Found ) CYCLE
      
      n = Element % Type % NumberOfNodes
      ElemPerm(1:n) = Perm( Element % NodeIndexes(1:n) )
      IF( ANY(ElemPerm(1:n) == 0 ) ) CYCLE
      
      BlockIndex( ElemPerm(1:n) ) = j
    END DO
    
    n = COUNT( BlockIndex == 0 )
    IF( n > 0 ) THEN
      CALL Info('BlockPickDofsPhysical','Number of indexes without block matrix index: '//I2S(n),Level=7)
      IF( MinBlock > 1 ) THEN
        k = 1
      ELSE
        MaxBlock = MaxBlock + 1
        k = MaxBlock
      END IF
      WHERE( BlockIndex == 0 ) BlockIndex = k
    ELSE
      CALL Info('BlockPickDofsPhysical','All physical domains given block index',Level=10)
    END IF
    
    MaxBlock = ParallelReduction(MaxBlock, 2 ) 
    NoVar = MaxBlock

  END SUBROUTINE BlockPickDofsPhysical
    


  !-------------------------------------------------------------------------------------
  !> Arranges the DOFs of a H(div) approximation into groups
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickHdiv( Solver, BlockIndex, NoVar )
    
    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER :: BlockIndex(:)
    INTEGER :: Novar
    
    INTEGER :: i,j,n,nn,ne,nf,nb,nnis,neis,nfis,nbis
    INTEGER :: nncount,necount,nfcount,nbcount
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found
    INTEGER, POINTER :: Perm(:)
    
    
    CALL Info('BlockPickHdiv','Picking block matrix for mixed hdiv solver',Level=10)

    Mesh => Solver % Mesh
    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges
    nf = Mesh % NumberOfFaces
    nb = Mesh % NumberOfBulkElements

    ! true/false flags whether dof type exists
    nnis = 0
    neis = 0
    nfis = 0
    nbis = 0

    ! counter of types of dofs
    nncount = 0
    necount = 0
    nfcount = 0
    nbcount = 0
    
    Perm => Solver % Variable % Perm 
    n = SIZE( Perm ) 

    DO i=1,n
      j = Perm(i)
      IF( j == 0 ) CYCLE

      IF( i <= nn ) THEN
        nnis = 1
        nncount = nncount + 1
        BlockIndex(j) = 1
      ELSE IF( i <= nn + ne ) THEN
        neis = 1
        necount = necount + 1
        BlockIndex(j) = nnis + 1
      ELSE IF( i <= nn + ne + nf ) THEN
        nfis = 1
        nfcount = nfcount + 1
        BlockIndex(j) = nnis + neis + 1
      ELSE
        nbis = 1
        nbcount = nbcount + 1
        BlockIndex(j) = nnis + neis + nfis + 1
      END IF
    END DO

    IF( nncount > 0 ) CALL Info('BlockPickHdiv','Number of nodal dofs: '//I2S(nncount),Level=8)
    IF( necount > 0 ) CALL Info('BlockPickHdiv','Number of edge dofs: '//I2S(necount),Level=8)
    IF( nfcount > 0 ) CALL Info('BlockPickHdiv','Number of face dofs: '//I2S(nfcount),Level=8)
    IF( nbcount > 0 ) CALL Info('BlockPickHdiv','Number of elemental dofs: '//I2S(nbcount),Level=8)
       
    NoVar = nnis + neis + nfis + nbis

    CALL Info('BlockPickHdiv','Found dofs related to '//I2S(NoVar)//' groups',Level=6)
    
  END SUBROUTINE BlockPickHdiv
  


  
  !-------------------------------------------------------------------------------------
  !> Picks the components of a full matrix when blockindex table is given.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrixPerm( Solver, BlockIndex, NoVar )

    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: BlockIndex(:)
    INTEGER :: Novar

    INTEGER :: bcol,brow,bi,bk,i,k,j,n
    TYPE(Matrix_t), POINTER :: A, B
    INTEGER, ALLOCATABLE :: BlockNumbering(:), rowcount(:), offset(:)
    
    CALL Info('BlockPickMatrixPerm','Picking domainwise block matrix from monolithic one',Level=10)

    A => Solver % Matrix 
    
    n = A % NumberOfRows
    ALLOCATE( BlockNumbering( n ), rowcount(NoVar), offset(NoVar+1) )
    BlockNumbering = 0
    RowCount = 0
    offset = 0

    IF(.NOT. ASSOCIATED(TotMatrix % BlockPerm) ) THEN
      ALLOCATE(TotMatrix % BlockPerm(n))
    END IF
    TotMatrix % BlockPerm = 0 
    
    
    DO i=1,n
      brow = BlockIndex(i)
      rowcount(brow) = rowcount(brow) + 1
      BlockNumbering(i) = rowcount(brow)
    END DO

    
    DO i = 1, NoVar
      B => TotMatrix % SubMatrix(i,i) % Mat
      n = rowcount(i)
      
      ALLOCATE(B % Rhs(n))
      B % rhs = 0.0_dp
      
      ALLOCATE(B % InvPerm(n))
      B % InvPerm = 0 
      ! Add the (n,n) entry since this helps to create most efficiently the full ListMatrix
      ! CALL AddToMatrixElement(B,n,n,0.0_dp)      

      offset(i+1) = offset(i) + n
    END DO
    
    DO i=1,A % NumberOfRows 
      
      brow = BlockIndex(i)
      bi = BlockNumbering(i)
     
      B => TotMatrix % SubMatrix(brow,brow) % Mat
      B % Rhs(bi) = A % Rhs(i)

      B % InvPerm(bi) = i
      
      TotMatrix % BlockPerm(offset(brow)+bi) = i
      
      DO j=A % Rows(i+1)-1,A % Rows(i),-1

        k = A % Cols(j)
        
        bcol = BlockIndex(k)
        bk = BlockNumbering(k)
        
        B => TotMatrix % SubMatrix(brow,bcol) % Mat       
        CALL AddToMatrixElement(B,bi,bk,A % Values(j))
      END DO
    END DO
    
    DO i = 1, NoVar
      DO j = 1, NoVar
        B => TotMatrix % SubMatrix(i,j) % Mat
        CALL List_toCRSMatrix(B)
      END DO
    END DO
              
  END SUBROUTINE BlockPickMatrixPerm

    
  
  !-------------------------------------------------------------------------------------
  !> Picks the components of a full matrix to the submatrices of a block matrix assuming AV solver.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrixAV( Solver, NoVar )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar

    INTEGER :: RowVar, ColVar
    TYPE(Matrix_t), POINTER :: SolverMatrix
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER::i,j,k,i_aa,i_vv,i_av,i_va,n;
    TYPE(Matrix_t), POINTER :: B_aa,B_av,B_va,B_vv,C_aa,C_vv,A,CM
    REAL(KIND=DP) :: SumAbsMat
    
    CALL Info('BlockPickMatrixAV','Picking block matrix from monolithic one',Level=10)

    SolverMatrix => Solver % Matrix 
    
    A => SolverMatrix
    i_aa=0; i_vv=0; i_av=0; i_va=0;
    n = Solver % Mesh % NumberOfNodes

    B_vv => TotMatrix % SubMatrix(1,1) % Mat
    B_va => TotMatrix % SubMatrix(1,2) % Mat
    B_av => TotMatrix % SubMatrix(2,1) % Mat
    B_aa => TotMatrix % SubMatrix(2,2) % Mat

    IF(ASSOCIATED(B_aa % Values)) B_aa % Values=0._dp
    IF(ASSOCIATED(B_av % Values)) B_av % Values=0._dp
    IF(ASSOCIATED(B_va % Values)) B_va % Values=0._dp
    IF(ASSOCIATED(B_vv % Values)) B_vv % Values=0._dp

    DO i=1,SIZE(Solver % Variable % Perm)
      j = Solver % Variable % Perm(i)
      IF(j<=0) CYCLE
      IF (i<=n) THEN
        i_vv=i_vv+1
      ELSE
        i_aa=i_aa+1
      END IF
      DO k=A % Rows(j+1)-1,A % Rows(j),-1
        !         IF(A % Cols(k) /= j.AND.A % Values(k)==0) CYCLE

        IF(i<=n.AND.A % Cols(k)<=n) THEN
          CALL AddToMatrixElement(B_vv,i_vv,A % Cols(k),A % Values(k))
        ELSE IF (i<=n.AND.A % Cols(k)>n) THEN
          CALL AddToMatrixElement(B_va,i_vv,A % Cols(k)-n,A % Values(k))
        ELSE IF (i>n.AND.A % Cols(k)<=n) THEN
          CALL AddToMatrixElement(B_av,i_aa,A % Cols(k),A % Values(k))
        ELSE
          CALL AddToMatrixElement(B_aa,i_aa,A % Cols(k)-n,A % Values(k))
        END IF
      END DO
    END DO

    IF (B_aa % Format == MATRIX_LIST) THEN
      CALL List_toCRSMatrix(B_aa)
      CALL List_toCRSMatrix(B_av)
      CALL List_toCRSMatrix(B_va)
      CALL List_toCRSMatrix(B_vv)
    END IF

    IF(ASSOCIATED(B_aa % Rhs)) THEN
      DEALLOCATE(B_aa % Rhs, B_vv % Rhs)
    END IF

    ALLOCATE(B_aa % Rhs(B_aa % NumberOfRows))
    ALLOCATE(B_vv % Rhs(B_vv % NumberOfRows))

    i_vv=0; i_aa=0
    DO i=1,SIZE(Solver % Variable % Perm)
      j = Solver % Variable % Perm(i)
      IF(j<=0) CYCLE
      IF (i<=n) THEN
        i_vv=i_vv+1
        B_vv % Rhs(i_vv) = A % Rhs(j)
      ELSE 
        i_aa=i_aa+1
        B_aa % Rhs(i_aa) = A % Rhs(j)
      END IF
    END DO


    ! Also inherit the constraints, if any
    ! If the constraints are treated as block matrix also the
    ! pointer should not be associated. 
    CM => A % ConstraintMatrix
    DO WHILE(ASSOCIATED(CM))
      C_aa=>AllocateMatrix(); C_aa % Format=MATRIX_LIST
      C_vv=>AllocateMatrix(); C_vv % Format=MATRIX_LIST
      i_aa=0; i_vv=0;
      DO i=1,CM % NumberOFRows
        IF(CM % Cols(CM % Rows(i))<=n) THEN
          i_vv=i_vv+1
        ELSE
          i_aa=i_aa+1
        END IF
        DO j=CM % Rows(i),CM % Rows(i+1)-1
          IF (CM % Values(j)==0._dp) CYCLE

          IF (CM % Cols(j)<=n) THEN
            CALL AddToMatrixElement(C_vv,i_vv,CM % Cols(j),CM % Values(j))
          ELSE
            CALL AddToMatrixElement(C_aa,i_aa,CM % Cols(j)-n,CM % Values(j))
          END IF
        END DO
      END DO
      C_aa % ConstraintMatrix => Null()!B_aa % ConstraintMatrix XXXXXXX
      B_aa % ConstraintMatrix => C_aa

      C_vv % ConstraintMatrix => Null()!B_vv % ConstraintMatrix YYYYYYY
      B_vv % ConstraintMatrix => C_vv

      CALL List_toCRSMatrix(C_vv)
      CALL List_toCRSMatrix(C_aa)
      ALLOCATE(C_aa % Rhs(C_aa % NumberOfRows)); C_aa % Rhs=0._dp
      ALLOCATE(C_vv % Rhs(C_vv % NumberOfRows)); C_vv % Rhs=0._dp
      CM => CM % ConstraintMatrix
      IF(c_vv%numberofrows<=0) b_vv%constraintmatrix=>null()
    END DO
    
  END SUBROUTINE BlockPickMatrixAV


  !-------------------------------------------------------------------------------------
  !> Picks vertical and horizontal components of a full matrix.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrixHorVer( Solver, NoVar, Cart )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar
    LOGICAL :: Cart

    INTEGER::i,j,k,n,t,ne,dofs,nd,ni,nn,ndir(3),ic,kc
    TYPE(Matrix_t), POINTER :: A,B
    TYPE(Nodes_t), SAVE :: Nodes, EdgeNodes
    INTEGER :: ActiveCoordinate
    INTEGER, ALLOCATABLE :: DTag(:), DPerm(:)
    INTEGER, POINTER :: Indexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Wlen, Wproj, Wtol, u, v, w, DetJ, Normal(3), MaxCoord, MinCoord
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: PiolaVersion, Found, Stat
    TYPE(Element_t), POINTER :: Element, Edge
    REAL(KIND=dp), POINTER :: Coord(:)
    
    REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:), RotWBasis(:,:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)

    n = 28 ! currently just large enough
    ALLOCATE( WBasis(n,3), RotWBasis(n,3), Basis(n), dBasisDx(n,3), Indexes(n) )
    
    
    CALL Info('BlockPickMatrixHorVer','Dividing matrix into vertical and horizontal dofs',Level=10)


    n = MAXVAL(Solver % Variable % Perm)
    Mesh => Solver % Mesh 
    
    A => Solver % Matrix
    dofs = Solver % Variable % Dofs
    
    n = A % NumberOfRows / dofs
    
    ALLOCATE( DTag(n), DPerm(n*dofs)  ) 
    DTag = 0
    DPerm = 0
        
    PiolaVersion = ListGetLogical( Solver % Values,'Use Piola Transform', Found )
    ActiveCoordinate = ListGetInteger( Solver % Values,'Active Coordinate',Found )
    IF(.NOT. Found ) ActiveCoordinate = 3
    Normal = 0.0_dp
    Normal(ActiveCoordinate) = 1.0_dp
    
    Wtol = 1.0e-3
      
    
    DO t=1,Solver % NumberOfActiveElements
      Element => Mesh % Elements( Solver % ActiveElements(t) )
      nn = Element % TYPE % NumberOfNodes

      nd = GetElementDOFs( Indexes, Element, Solver)  
      CALL GetElementNodes( Nodes, Element )


      ! Both strategies give exactly the same set of vertical and horizontal dofs!
      ! Both strategies give exactly the same set of vertical and horizontal dofs!
      IF( Cart ) THEN
        DO ActiveCoordinate = 1, 3
          IF( ActiveCoordinate == 1 ) THEN
            Coord => Nodes % x
          ELSE IF( ActiveCoordinate == 2 ) THEN
            Coord => Nodes % y
          ELSE
            Coord => Nodes % z
          END IF
          
          MinCoord = MINVAL( Coord(1:nn) )
          MaxCoord = MAXVAL( Coord(1:nn) )
          Wlen = MaxCoord - MinCoord 
          
          DO i=1,nd
            j = Solver % Variable % Perm(Indexes(i))
            
            IF( i <= Element % TYPE % NumberOfEdges ) THEN
              Edge => Mesh % Edges( Element % EdgeIndexes(i) )
              CALL GetElementNodes( EdgeNodes, Edge )
              ne = Edge % TYPE % NumberOfNodes
              
              IF( ActiveCoordinate == 1 ) THEN
                Coord => EdgeNodes % x
              ELSE IF( ActiveCoordinate == 2 ) THEN
                Coord => EdgeNodes % y
              ELSE
                Coord => EdgeNodes % z
              END IF
              
              MinCoord = MINVAL( Coord(1:ne) )
              MaxCoord = MAXVAL( Coord(1:ne) )
            ELSE            
              CALL Fatal('BlockPickMatrixHorVer','Cannot do faces yet!')
            END IF

            Wproj = ( MaxCoord - MinCoord ) / Wlen
            
            IF( WProj > 1.0_dp - Wtol ) DTag(j) = ActiveCoordinate
          END DO
        END DO

      ELSE IF(.TRUE.) THEN
        IF( ActiveCoordinate == 1 ) THEN
          Coord => Nodes % x
        ELSE IF( ActiveCoordinate == 2 ) THEN
          Coord => Nodes % y
        ELSE
          Coord => Nodes % z
        END IF

        MinCoord = MINVAL( Coord(1:nn) )
        MaxCoord = MAXVAL( Coord(1:nn) )
        Wlen = MaxCoord - MinCoord 

        DO i=1,nd
          j = Solver % Variable % Perm(Indexes(i))

          IF( i <= Element % Type % NumberOfEdges ) THEN
            Edge => Mesh % Edges( Element % EdgeIndexes(i) )
            CALL GetElementNodes( EdgeNodes, Edge )
            ne = Edge % Type % NumberOfNodes

            IF( Indexes(i) /= Mesh % NumberOfNodes + Element % EdgeIndexes(i) ) THEN
              PRINT *,'ind com:',Indexes(i), Mesh % NumberOfNodes + Element % EdgeIndexes(i), &
                  Mesh % NumberOfNodes 
            END IF

            IF( ActiveCoordinate == 1 ) THEN
              Coord => EdgeNodes % x
            ELSE IF( ActiveCoordinate == 2 ) THEN
              Coord => EdgeNodes % y
            ELSE
              Coord => EdgeNodes % z
            END IF

            MinCoord = MINVAL( Coord(1:ne) )
            MaxCoord = MAXVAL( Coord(1:ne) )
          ELSE            
            ! jj = 2 * ( Element % ElementIndex - 1) + ( i - noedges ) 
            CALL Fatal('BlockPickMatrixHorVer','Cannot do faces yet!')
          END IF

          Wproj = ( MaxCoord - MinCoord ) / Wlen

          IF( WProj > 1.0_dp - Wtol ) THEN  
            DTag(j) = 1  
          ELSE IF( Wproj < Wtol ) THEN
            DTag(j) = 2  
          END IF
        END DO

      ELSE      
        IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion)

        u = SUM( IP % u ) / IP % n
        v = SUM( IP % v ) / IP % n
        w = IP % w(k)

        IF (PiolaVersion) THEN
          stat = EdgeElementInfo( Element, Nodes, u, v, w, &
              DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
              RotBasis = RotWBasis, dBasisdx = dBasisdx, &
              ApplyPiolaTransform = .TRUE.)
        ELSE
          stat = ElementInfo( Element, Nodes, u, v, w, &
              detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element, WBasis, RotWBasis, Basis, dBasisdx)
        END IF

        DO i=1,nd
          j = Solver % Variable % Perm(Indexes(i))
          Wlen = SQRT( SUM( WBasis(i,:)**2 ) )
          IF( Wlen < EPSILON( Wlen ) ) CYCLE

          Wproj = ABS( SUM( WBasis(i,:) * Normal ) ) / Wlen 

          IF( WProj > 1.0_dp - Wtol ) THEN  
            IF( DTag(j) == 2 ) PRINT *,'Vertical edge '//I2S(j)//' is also horizontal?'
            DTag(j) = 1  ! set to be vertical
          ELSE IF( Wproj < Wtol ) THEN
            IF( DTag(j) == 1 ) PRINT *,'Horizontal edge '//I2S(j)//' is also vertical?'
            DTag(j) = 2  ! set to be horizontal
          ELSE
            PRINT *,'Edge '//I2S(j)//' direction undefined: ',Wproj
          END IF
        END DO

      END IF
    END DO
    

    ! Number vertical and horizontal (or all cartesian) dofs separately.
    ndir = 0
    DO i=1,n
      DO j=1,dofs
        k = dofs*(i-1)+j
        ndir(DTag(i)) = ndir(DTag(i)) + 1
        DPerm(k) = ndir(DTag(i))
      END DO
    END DO

    !PRINT *,'Cartesian dofs:',ndir(1:NoVar)

    i = n - SUM( ndir ) 
    IF( i > 0 ) THEN      
      CALL Fatal('BlockPickMatrixHorVer','Could not determine all nodes: '&
          //I2S(i))
    END IF


    ! Allocate vectors if not present
    DO i=1,NoVar
      DO j=1,NoVar
        B => TotMatrix % SubMatrix(i,j) % Mat
        IF( ASSOCIATED( B % Values ) ) B % Values = 0.0_dp
      END DO
      B => TotMatrix % SubMatrix(i,i) % Mat      
      IF(.NOT. ASSOCIATED( B % InvPerm ) ) ALLOCATE( B % InvPerm(ndir(i)) )
      IF(.NOT. ASSOCIATED( B % Rhs) ) ALLOCATE(B % Rhs(ndir(i)) )
      !PRINT *,'a complex', a % complex
      !B % COMPLEX = A % COMPLEX
    END DO
    

    DO i=1,A % NumberOfRows
      ic = (i-1)/dofs+1
      
      DO j=A % Rows(i+1)-1,A % Rows(i),-1
        k = A % Cols(j)
        kc = (k-1)/dofs+1
        
        IF( DTag(ic) < 1 .OR. DTag(ic) > NoVar ) THEN
          PRINT *,'i:',i,ic,Dtag(ic)
        END IF
        
        IF( DTag(kc) < 1 .OR. DTag(kc) > NoVar ) THEN
          PRINT *,'k:',k,kc,Dtag(kc)
        END IF
        
        B => TotMatrix % SubMatrix(DTag(ic),DTag(kc)) % Mat
        
        IF( Dperm(i) < 1 .OR. DPerm(k) < 1 ) THEN
          PRINT *,'ik',Dperm(i),Dperm(k)
          STOP EXIT_ERROR
        END IF
        CALL AddToMatrixElement(B,Dperm(i),DPerm(k),A % Values(j))
      END DO

      B => TotMatrix % SubMatrix(DTag(ic),DTag(ic)) % Mat      
      B % Rhs(Dperm(i)) = A % Rhs(i)          
      B % InvPerm(Dperm(i)) = i          
    END DO
    
    DO i=1,NoVar
      DO j=1,NoVar
        B => TotMatrix % SubMatrix(i,j) % Mat        
        IF (B % FORMAT == MATRIX_LIST) THEN
          CALL List_toCRSMatrix(B)
        END IF
      END DO
    END DO

    
    IF( ASSOCIATED( A % ConstraintMatrix ) ) THEN
      CALL Warn('BlockPickMatrixHorVer','Cannot deal with constraints')
    END IF
    
  END SUBROUTINE BlockPickMatrixHorVer


  !-------------------------------------------------------------------------------------
  !> Makes a quadratic H(curl) approximation to have a block structure
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrixHcurl( Solver, NoVar, DoCmplx )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar
    LOGICAL :: DoCmplx

    INTEGER :: i,j,k,n,m,n0,dofs,ic,kc
    TYPE(Matrix_t), POINTER :: A,B
    INTEGER, ALLOCATABLE :: BlockTag(:), BlockPerm(:), nblock(:)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element, Edge

    Mesh => Solver % Mesh

    IF(.NOT. ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Fatal('BlockPickMatrixHcurl','This subroutine needs Edges!')
    END IF
    IF(.NOT. ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Fatal('BlockPickMatrixHcurl','This subroutine needs Faces!')
    END IF
        
    CALL Info('BlockPickMatrixHcurl','Arranging a quadratic H(curl) approximation into blocks',Level=10)

    A => Solver % Matrix
    dofs = Solver % Variable % Dofs

    m = A % NumberOfRows    
    n = m
    IF(.NOT. DoCmplx) n = n / dofs
      
    ALLOCATE( BlockTag(n), BlockPerm(m), nblock(NoVar) )

    ! Set the default blocks
    IF( DoCmplx ) THEN
      ! Synopsis: the tags of (4x4) block structure
      !   - Re, lowest-order = 1
      !   - Im, lowest-order = 2
      !   - Re, higher-order = 3
      !   - Im, higher-order = 4
      ! while (3x3) structure corresponds to
      !   - Re and Im, lowest-order = 1
      !   - Re, higher-order = 2
      !   - Im, higher-order = 3
      BlockTag(1::2) = 3
      BlockTag(2::2) = 4
    ELSE
      ! Synopsis:
      !   - lowest-order = 1
      !   - higher-order = 2
      !   - higher-order DOFs which are not associated with edges = 3 (optional)
      BlockTag = 2
    END IF
        
    n0 = Mesh % NumberOfNodes    
    DO i=1, Mesh % NumberOfEdges
      ! This corresponds to the lowest-order DOF over an edge
      j = n0 + 2*i-1
      k = Solver % Variable % Perm(j)
      IF(k==0) CYCLE
      IF( DoCmplx ) THEN
        BlockTag(2*k-1) = 1
        BlockTag(2*k) = 2 
      ELSE
        ! If BlockTag array were created for all DOFs with DOFs>1,
        ! it would have a repeated entries occuring in clusters of
        ! size DOFs. Therefore a smaller array can be used to tag DOFs.
        BlockTag(k) = 1
      END IF
    END DO

    IF(NoVar == 3) THEN
      IF( DoCmplx ) THEN
        WHERE( BlockTag > 1 )
          BlockTag = BlockTag - 1
        END WHERE
      ELSE
        DO j = n0 + 2*Mesh % NumberOfEdges + 1, SIZE(Solver % Variable % Perm)
          k = Solver % Variable % Perm(j)
          IF(k>0) BlockTag(k) = 3
        END DO
      END IF
    END IF
    
    ! Number each group of DOFs separately
    BlockPerm = 0
    nblock = 0
    IF( DoCmplx ) THEN
      DO i=1,n
        nblock(BlockTag(i)) = nblock(BlockTag(i)) + 1
        BlockPerm(i) = nblock(BlockTag(i))
      END DO
    ELSE
      DO i=1,n
        n0 = dofs*(i-1)
        DO j=1,dofs
          nblock(BlockTag(i)) = nblock(BlockTag(i)) + 1
          k = n0 + j
          BlockPerm(k) = nblock(BlockTag(i))
        END DO
      END DO
    END IF

    DO i=1,NoVar
      CALL Info('BlockPickMatrixHcurl','Block vector '//I2S(i)//' size: '//I2S(nblock(i)),Level=6)
    END DO

    ! Allocate vectors if not present
    DO i=1,NoVar
      DO j=1,NoVar
        B => TotMatrix % SubMatrix(i,j) % Mat
        IF( ASSOCIATED( B % Values ) ) B % Values = 0.0_dp
      END DO
      B => TotMatrix % SubMatrix(i,i) % Mat      
      IF(.NOT. ASSOCIATED( B % InvPerm ) ) ALLOCATE( B % InvPerm(nblock(i)) )
      IF(.NOT. ASSOCIATED( B % Rhs) ) ALLOCATE(B % Rhs(nblock(i)) )
      !PRINT *,'a complex', a % complex
      !B % COMPLEX = A % COMPLEX
    END DO
    

    DO i=1,A % NumberOfRows
      IF( DoCmplx ) THEN
        ic = i
      ELSE       
        ic = (i-1)/dofs+1
      END IF
        
      DO j=A % Rows(i+1)-1,A % Rows(i),-1
        k = A % Cols(j)

        IF( DoCmplx ) THEN
          kc = k
        ELSE
          kc = (k-1)/dofs+1
        END IF
          
        IF( BlockTag(ic) < 1 .OR. BlockTag(ic) > NoVar ) THEN
          PRINT *,'i:',i,ic,BlockTag(ic)
        END IF
        
        IF( BlockTag(kc) < 1 .OR. BlockTag(kc) > NoVar ) THEN
          PRINT *,'k:',k,kc,BlockTag(kc)
        END IF
        
        B => TotMatrix % SubMatrix(BlockTag(ic),BlockTag(kc)) % Mat
        
        IF( BlockPerm(i) < 1 .OR. BlockPerm(k) < 1 ) THEN
          PRINT *,'ik',BlockPerm(i),BlockPerm(k)
          STOP EXIT_ERROR
        END IF
        CALL AddToMatrixElement(B,BlockPerm(i),BlockPerm(k),A % Values(j))
      END DO

      B => TotMatrix % SubMatrix(BlockTag(ic),BlockTag(ic)) % Mat      
      B % Rhs(BlockPerm(i)) = A % Rhs(i)          
      B % InvPerm(BlockPerm(i)) = i          
    END DO
    
    DO i=1,NoVar
      DO j=1,NoVar
        B => TotMatrix % SubMatrix(i,j) % Mat        
        IF (B % FORMAT == MATRIX_LIST) THEN
          CALL List_toCRSMatrix(B)
        END IF
        CALL Info('BlockPickMatrixHcurl','Matrix '//I2S(10*i+j)//' nonzeros: '//I2S(SIZE(B % Values)),Level=6)
      END DO
    END DO
    
    IF( ASSOCIATED( A % ConstraintMatrix ) ) THEN
      CALL Warn('BlockPickMatrixHcurl','Cannot deal with constraints')
    END IF
    
  END SUBROUTINE BlockPickMatrixHcurl
  

  !-------------------------------------------------------------------------------------
  !> Picks the components of a full matrix to the submatrices of a block matrix assuming AV solver.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickMatrixNodal( Solver, NoVar )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar

    INTEGER :: RowVar, ColVar, Dofs
    TYPE(Matrix_t), POINTER :: SolverMatrix
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER::i,j,k,ii,l,ll,i_aa,i_vv,i_av,i_va,n_a,n_v,ntot,rdof,cdof
    TYPE(Matrix_t), POINTER :: B_aa,B_av,B_va,B_vv,C_aa,C_vv,A,CM
    REAL(KIND=DP) :: SumAbsMat, val
    
    CALL Info('BlockPickMatrixNodal','Picking nodal and non-nodal block matrices from monolithic one',Level=10)

    SolverMatrix => Solver % Matrix 
    
    A => SolverMatrix
    i_aa=0; i_vv=0; i_av=0; i_va=0;
    i = Solver % Mesh % NumberOfNodes
    n_v = MAXVAL( Solver % Variable % Perm(1:i) )
    ntot = MAXVAL( Solver % Variable % Perm ) 
    n_a = ntot - n_v    
    
    dofs = Solver % Variable % Dofs

    !PRINT *,'Dofs, n_v, n_a: ',dofs,n_v,n_a
    
    B_vv => TotMatrix % SubMatrix(1,1) % Mat
    B_va => TotMatrix % SubMatrix(1,2) % Mat
    B_av => TotMatrix % SubMatrix(2,1) % Mat
    B_aa => TotMatrix % SubMatrix(2,2) % Mat

    IF(ASSOCIATED(B_aa % Values)) B_aa % Values=0._dp
    IF(ASSOCIATED(B_av % Values)) B_av % Values=0._dp
    IF(ASSOCIATED(B_va % Values)) B_va % Values=0._dp
    IF(ASSOCIATED(B_vv % Values)) B_vv % Values=0._dp

    IF( .NOT. ASSOCIATED( B_aa % Rhs ) ) ALLOCATE(B_aa % Rhs(dofs*n_a))
    IF( .NOT. ASSOCIATED( B_vv % Rhs ) ) ALLOCATE(B_vv % Rhs(dofs*n_v))

    
    DO i=1,A % NumberOfRows
      
      ii = (i-1)/dofs+1
      rdof = (i-1)/ntot+1

      IF( ii < n_v ) THEN
        B_vv % Rhs(i-(rdof-1)*n_a) = A % Rhs(i)
      ELSE
        B_aa % Rhs(i-rdof*n_v) = A % Rhs(i)
      END IF
      
      DO k=A % Rows(i+1)-1,A % Rows(i),-1
        l = A % Cols(k)
        
        ll = (l-1)/dofs+1
        cdof = (l-1)/ntot+1
        val = A % Values(k)
                    
        IF( ii <= n_v ) THEN
          IF( ll <= n_v ) THEN
            CALL AddToMatrixElement(B_vv,i-(rdof-1)*n_a,l-(cdof-1)*n_a,val)
          ELSE
            IF( ABS( val ) > TINY( val ) ) THEN
              CALL AddToMatrixElement(B_va,i-(rdof-1)*n_a,l-cdof*n_v,val)
            END IF
          END IF
        ELSE
          IF( ll <= n_v ) THEN
            IF( ABS( val ) > TINY( val ) ) THEN
              CALL AddToMatrixElement(B_av,i-rdof*n_v,l-(cdof-1)*n_a,val)
            END IF
          ELSE
            CALL AddToMatrixElement(B_aa,i-rdof*n_v,l-cdof*n_v,val)
          END IF
        END IF
      END DO
    END DO

    IF (B_aa % Format == MATRIX_LIST) THEN
      CALL List_toCRSMatrix(B_aa)
      CALL List_toCRSMatrix(B_av)
      CALL List_toCRSMatrix(B_va)
      CALL List_toCRSMatrix(B_vv)
    END IF
   
    CM => A % ConstraintMatrix
    IF( ASSOCIATED( CM ) ) THEN
      CALL Fatal('BlockPickMatrixNodal','There would be some constraints to pick too!')
    END IF
    
  END SUBROUTINE BlockPickMatrixNodal



  !-------------------------------------------------------------------------------------
  !> Picks the components of the constraint matrix.
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPickConstraint( Solver, NoVar )
    
    TYPE(Solver_t), TARGET :: Solver
    INTEGER :: Novar

    INTEGER :: RowVar, ColVar
    TYPE(Matrix_t), POINTER :: SolverMatrix
    TYPE(Matrix_t), POINTER :: Amat
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found

    LOGICAL :: BlockAV
    INTEGER::i,j,k,l,n,i1,i2,i3,rowi,colj,NoCon,rb,cb
    TYPE(Matrix_t), POINTER :: A,CM,C1,C2,C3,C1prec,C2prec,C3prec
    REAL(KIND=dp) :: PrecCoeff,val
    INTEGER, POINTER :: ConsPerm(:)
    INTEGER :: DoPrec
    CHARACTER(:), ALLOCATABLE :: VarName
    TYPE(Variable_t), POINTER :: Var
    TYPE(Solver_t), POINTER :: PSolver
    LOGICAL :: InheritCM, PrecTrue 
    
    CALL Info('BlockPickConstraint','Picking constraints to block matrix',Level=10)

    
    SolverMatrix => Solver % Matrix 
    Params => Solver % Values
       
    BlockAV = ListGetLogical( Params,'Block A-V System', Found)

    A => SolverMatrix
    n = Solver % Mesh % NumberOfNodes
    
    PrecCoeff = ListGetConstReal( Params,'Block Diag Coeff',Found)
    IF(.NOT. Found ) PrecCoeff = 1.0_dp

    PrecTrue = ListGetLogical( Params,'Block Diag True',Found ) 
    
    
    ! temporarily be generic
    NoCon = 0
    CM => A % ConstraintMatrix
    DO WHILE(ASSOCIATED(CM))
      NoCon = NoCon + 1
      CM => CM % ConstraintMatrix
    END DO

    CALL Info('BlockPickConstraint','Number of constraint matrices: '//I2S(NoCon),Level=10)
    
    InheritCM = (NoVar == 1 ) .AND. (NoCon == 1 )
    
    
    IF( NoVar == 1 ) THEN
      IF( InheritCM ) THEN
        TotMatrix % Submatrix(NoVar+1,1) % Mat => A % ConstraintMatrix
        CALL Info('BlockPickConstraint','Using constraint matrix as is!',Level=10)
        IF( PrecTrue ) THEN
          C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % Mat                  
        ELSE
          C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % PrecMat        
        END IF
      ELSE      
        C1 => TotMatrix % Submatrix(NoVar+1,1) % Mat
        IF( PrecTrue ) THEN
          C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % Mat        
        ELSE
          C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % PrecMat        
        END IF
      END IF
      TotMatrix % SubMatrixTranspose(NoVar+1,1) = .TRUE.
      
    ELSE IF(BlockAV) THEN          
      C1 => TotMatrix % SubMatrix(NoVar+1,1) % Mat
      C2 => TotMatrix % Submatrix(NoVar+2,2) % Mat
      IF( PrecTrue ) THEN
        C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % Mat
        C2prec => TotMatrix % Submatrix(NoVar+2,NoVar+2) % Mat
      ELSE
        C1prec => TotMatrix % Submatrix(NoVar+1,NoVar+1) % PrecMat
        C2prec => TotMatrix % Submatrix(NoVar+2,NoVar+2) % PrecMat
      END IF
      TotMatrix % SubMatrixTranspose(NoVar+1,1) = .TRUE.
      TotMatrix % SubMatrixTranspose(NoVar+2,2) = .TRUE.
    ELSE
      CALL Fatal('BlockPickConstraint','Not done for vectors!')
    END IF
    
    n = Solver % Mesh % NumberOfNodes
    CM => A % ConstraintMatrix
    DO WHILE(ASSOCIATED(CM)) 
      n = MAX( n, MAXVAL( CM % InvPerm ) )
      CM => CM % ConstraintMatrix
    END DO

    
    ALLOCATE( ConsPerm( n ) ) 
    ConsPerm = 0
    

    DO DoPrec = 0, 1
      
      i1 = 0
      i2 = 0
      
      CM => A % ConstraintMatrix
      DO WHILE(ASSOCIATED(CM))         

        DO i=1,CM % NumberOFRows

          rowi = CM % Cols(CM % Rows(i))

          rb = 1          
          IF( BlockAV ) THEN
            IF( rowi > n ) rb = 2
          END  IF

          IF( rb == 1 ) THEN
            i1 = i1 + 1
          ELSE
            i2 = i2 + 1
          END IF

          ! First round initialize the ConsPerm
          IF( DoPrec == 0 ) THEN
            j = CM % InvPerm(i)
            ConsPerm( j ) = i
          END IF
                      
          DO j=CM % Rows(i),CM % Rows(i+1)-1
            IF (CM % Values(j)==0._dp) CYCLE

            colj = CM % Cols(j) 
            cb = 1
            IF( BlockAV ) THEN
              IF( colj > n ) THEN
                cb = 2
                colj = colj - n 
              END IF
            END IF

            val = CM % Values(j)

            IF( DoPrec == 1 ) THEN
              IF( ConsPerm( colj ) > 0 ) THEN
                ! The sign -1 is needed for consistency
                val = -PrecCoeff * val
                IF ( cb == 1 ) THEN
                  CALL AddToMatrixElement(C1prec,i1,ConsPerm(colj),val)
                ELSE
                  CALL AddToMatrixElement(C2prec,i2,ConsPerm(colj),val)
                END IF
              END IF                
            ELSE IF( .NOT. InheritCM ) THEN
              IF ( cb == 1 ) THEN
                CALL AddToMatrixElement(C1,i1,colj,val)
             ELSE
                CALL AddToMatrixElement(C2,i2,colj,val)
              END IF
            END IF
          END DO
        END DO

        CM => CM % ConstraintMatrix 
      END DO

      ! It is more efficient to set the last entry of the list matrix first      
      IF( DoPrec == 0 ) THEN
        CALL AddToMatrixElement(C1prec,i1,i1,0.0_dp)
      END IF
               
    END DO
      
    CALL Info('BlockPickConstraint','Setting format of constraint blocks to CRS',Level=20)
    IF(.NOT. InheritCM ) THEN
      CALL List_toCRSMatrix(C1)
    END IF
    CALL List_toCRSMatrix(C1prec)

    IF( BlockAV ) THEN
      CALL List_toCRSMatrix(C2)    
      CALL List_toCRSMatrix(C2prec)
    END IF

    
    IF( ListGetLogical( Solver % Values,'Save Prec Matrix', Found ) ) THEN   
      CALL SaveProjector(C1prec,.TRUE.,"CM")
    END IF
    
    
    VarName = "lambda"
    Var => VariableGet( Solver % Mesh % Variables, VarName )
    IF(ASSOCIATED( Var ) ) THEN
      CALL Info('BlockPickConstraint','Using existing variable > '//VarName//' <')		
    ELSE		
      CALL Info('BlockPickConstraint','Variable > '//VarName//' < does not exist, creating')
      PSolver => Solver
      
      n = i1
      Var => CreateBlockVariable(PSolver, NoVar+1, VarName, 1, ConsPerm )      

      TotMatrix % SubVector(NoVar+1) % Var => Var      
      TotMatrix % Offset(NoVar+2) = TotMatrix % Offset(NoVar+1) + n
      TotMatrix % MaxSize = MAX( TotMatrix % MaxSize, n )
      TotMatrix % TotSize = TotMatrix % TotSize + n
    END IF
    

  END SUBROUTINE BlockPickConstraint


  
  !-------------------------------------------------------------------------------------
  !> The block preconditioning matrix need not be directly derived from the full 
  !> matrix. Some or all the components may also be derived from a basic operator
  !> such as the Laplacian. 
  !-------------------------------------------------------------------------------------
  SUBROUTINE BlockPrecMatrix( Solver, NoVar )

    TYPE(Solver_t) :: Solver
    INTEGER :: Novar

    INTEGER :: i, RowVar, ColVar, CopyVar
    CHARACTER(:), ALLOCATABLE :: str
    REAL(KIND=dp) :: Coeff
    LOGICAL :: GotIt, GotIt2, DoIt
    INTEGER, POINTER :: VarPerm(:)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Amat, PMat
    TYPE(Variable_t), POINTER :: AVar
    
    CALL Info('BlockPrecMatrix','Checking for tailored preconditioning matrices',Level=6)
    
    Params => Solver % Values
    
    ! The user may give a user defined preconditioner matrix
    !-----------------------------------------------------------
    DO RowVar=1,NoVar
      i = TotMatrix % Submatrix(RowVar,RowVar) % PrecMat % NumberOfRows 

      IF( i > 0 ) CYCLE
      
      str = 'Prec Matrix Diffusion '//I2S(RowVar)
      Coeff = ListGetCReal( Params, str, GotIt)
      
      str = 'Prec Matrix Density '//I2S(RowVar)
      Coeff = ListGetCReal( Params, str, GotIt2)
      
      IF( GotIt .OR. GotIt2 ) THEN        
        CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(RowVar,RowVar) % Mat, &
            TotMatrix % Submatrix(RowVar,RowVar) % PrecMat )   
        
        Amat => TotMatrix % Submatrix(RowVar,RowVar) % PrecMat
        VarPerm => TotMatrix % Subvector(RowVar) % Var % Perm
        IF( GotIt ) THEN
          CALL Info('BlockPrecMatrix','Creating simple preconditioning Laplace matrix',Level=8)
          CALL LaplaceMatrixAssembly( Solver, VarPerm, Amat )
          Amat % Values = Coeff * Amat % Values
        ELSE 
          CALL Info('BlockPrecMatrix','Creating simple preconditioning mass matrix',Level=8)
          CALL MassMatrixAssembly( Solver, VarPerm, Amat )
          Amat % Values = Coeff * Amat % Values
        END IF
        Amat % ParallelInfo => TotMatrix % Submatrix(RowVar,RowVar) % Mat % ParallelInfo
      END IF
      
      str = 'Prec Matrix Complex Coeff '//I2S(RowVar)      
      Coeff = ListGetCReal( Params, str, GotIt )

      IF(.NOT. GotIt) THEN
        str = 'Prec Matrix Complex Coeff'
        Coeff = ListGetCReal( Params, str, GotIt )        
      END IF
      
      IF(.NOT. GotIt) THEN
        GotIt = ListGetLogical( Params,'Block Split Complex', GotIt )  
        Coeff = 1.0_dp
      END IF     
      
      IF( GotIt ) THEN
        IF( NoVar /= 2 .AND. NoVar /= 4 ) THEN
          CALL Fatal('BlockPrecMatrix','Assuming 2 or 4 blocks for the complex preconditioner!')
        END IF

        CALL Info('BlockPrecMatrix','Creating preconditioning matrix from block sums',Level=8)       
        CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(RowVar,RowVar) % Mat, &
            TotMatrix % Submatrix(RowVar,RowVar) % PrecMat )   
        Amat => TotMatrix % Submatrix(RowVar,RowVar) % PrecMat        
        IF( ASSOCIATED( TotMatrix % Submatrix(RowVar,RowVar) % Mat % PrecValues ) ) THEN
          AMat % Values = TotMatrix % Submatrix(RowVar,RowVar) % Mat % PrecValues                
          DEALLOCATE( TotMatrix % Submatrix(RowVar,RowVar) % Mat % PrecValues )
        ELSE
          AMat % Values = TotMatrix % Submatrix(RowVar,RowVar) % Mat % Values                
        END IF
          
        IF( RowVar == 1 .OR. RowVar == 3 ) THEN
          ColVar = RowVar + 1
        ELSE
          ColVar = RowVar - 1
          Coeff = -Coeff
        END IF
        IF( SIZE( Amat % Values ) /= SIZE( TotMatrix % Submatrix(RowVar,ColVar) % Mat % Values ) ) THEN
          CALL Fatal('BlockPrecMatrix','Mismatch in matrix size!')
        END IF
        
        AMat % Values = Amat % Values + &
            Coeff * TotMatrix % Submatrix(RowVar,ColVar) % Mat % Values                
      END IF
      
      str = 'Prec Matrix Diagonal '//I2S(RowVar)
      Coeff = ListGetCReal( Params, str, GotIt)
      IF( GotIt ) THEN
        CopyVar = NoVar+1 - RowVar
        PMat => TotMatrix % Submatrix(RowVar,CopyVar) % Mat
        Amat => TotMatrix % Submatrix(RowVar,RowVar) % PrecMat 
        CALL Info('BlockPrecMatrix','Creating preconditioner from matrix ('&
            //I2S(RowVar)//','//I2S(CopyVar)//')',Level=6)
        PRINT *,'proj matrix sum:',SUM( ABS( PMat % Values ) ) 
        PRINT *,'orig diag matrix sum:',SUM( ABS( TotMatrix % Submatrix(RowVar,RowVar) % Mat % Values ) )
        
        CALL DiagonalMatrixSumming( Solver, PMat, Amat )
        Amat % Values = Coeff * Amat % Values
      END IF
    END DO
    
    str = ListGetString( Params,'Block Matrix Schur Variable', GotIt)      
    IF( GotIt ) THEN
      AVAr => VariableGet( Solver % Mesh % Variables, str )
      IF( .NOT. ASSOCIATED( AVar ) ) THEN
        CALL Fatal('BlockPrecMatrix','Schur variable does not exist: '//str)
      END IF            
      IF( .NOT. ASSOCIATED( AVar % Solver ) ) THEN
        CALL Fatal('BlockPrecMatrix','Schur solver does not exist for: '//str)
      END IF
      IF( .NOT. ASSOCIATED( AVar % Solver % Matrix ) ) THEN
        CALL Fatal('BlockPrecMatrix','Schur matrix does not exist for: '//str)
      END IF
      CALL Info('BlockPrecMatrix','Using Schur matrix to precondition block '//I2S(NoVar))
      TotMatrix % Submatrix(NoVar,NoVar) % PrecMat => AVar % Solver % Matrix
    END IF  

    ! When we have an inner-outer iteration, we could well have a different matrix
    ! assembled for the purpose of preconditioning. Use it here, if available.
    IF(ListGetLogical( Params,'Block Nested System',GotIt ) ) THEN
      Amat => TotMatrix % Submatrix(1,1) % Mat
      IF( ASSOCIATED( Amat % PrecValues ) ) THEN
        PMat => TotMatrix % Submatrix(1,1) % PrecMat
        IF (.NOT. ASSOCIATED(PMat % Values)) THEN
          CALL Info('BlockPrecMatrix','Moving PrecValues to PrecMat!')
          CALL CRS_CopyMatrixTopology( AMat, PMat )
        ELSE
          ! Make a partial check that PrecMat has been derived from the right template:
          IF (.NOT. ASSOCIATED(AMat % Rows, PMat % Rows)) &
              CALL Fatal('BlockPrecMatrix', 'Inconsistent matrix structures')
        END IF
        PMat % Values => Amat % PrecValues
        NULLIFY(Amat % PrecValues)
      END IF
    END IF
    
  END SUBROUTINE BlockPrecMatrix


  ! Check if matrix C that is a coupling block in the block matrix system
  ! couples dofs at parallel interfaces. Assume that C := C_ab where a is
  ! associated to A and b to B.
  !------------------------------------------------------------------------
  FUNCTION CheckParallelCoupling( A, B, C ) RESULT ( Coupled ) 
    TYPE(Matrix_t), POINTER :: A, B, C
    LOGICAL :: Coupled
    
    LOGICAL :: Acoupled, Bcoupled
    INTEGER :: i,j,k
    REAL(KIND=dp) :: Eps
    
    Coupled = .FALSE.
    IF(.NOT. ASSOCIATED( A % ParallelInfo ) ) THEN
      CALL Fatal('CheckParallelCoupling','Matrix A does not have ParallelInfo!')
    END IF
    IF(.NOT. ASSOCIATED( B % ParallelInfo ) ) THEN
      CALL Fatal('CheckParallelCoupling','Matrix B does not have ParallelInfo!')
    END IF
    
    DO i=1,C % NumberOfRows
      DO j=C % Rows(i), C % Rows(i+1)-1
        k = C % Cols(j)
        IF( ABS( C % Values(j) ) < EPSILON( Eps ) ) CYCLE  
        IF ( ASSOCIATED(A % ParallelInfo % NeighbourList(i) % Neighbours) ) THEN
          IF ( SIZE(A % ParallelInfo % NeighbourList(i) % Neighbours) > 1 ) Coupled = .TRUE.
        END IF
        IF ( ASSOCIATED(B % ParallelInfo % NeighbourList(k) % Neighbours) ) THEN       
          IF ( SIZE(B % ParallelInfo % NeighbourList(k) % Neighbours) > 1 ) Coupled = .TRUE.
        END IF
        IF( Coupled ) EXIT
      END DO
    END DO

    IF( Coupled ) THEN
      CALL Info('CheckParallelCoupling','Coupling matrix has parallel connections!',Level=10)
    ELSE
      CALL Info('CheckParallelCoupling','Coupling matrix does not have parallel connections!',Level=10)
    END IF
      
  END FUNCTION CheckParallelCoupling
  

  !> Create the coupling blocks for a linear FSI coupling among various types of
  !> elasticity and fluid solvers.
  !--------------------------------------------------------------------------------
  SUBROUTINE FsiCouplingBlocks( Solver )

    TYPE(Solver_t) :: Solver
    INTEGER :: i, j, k, Novar

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A_fs, A_sf, A_s, A_f
    TYPE(Variable_t), POINTER :: FVar, SVar
    INTEGER, POINTER :: ConstituentSolvers(:)
    LOGICAL :: Found
    LOGICAL :: IsPlate, IsShell, IsNs, IsPres
    CHARACTER(*), PARAMETER :: Caller = 'FsiCouplingBlocks'
    
    Params => Solver % Values
    ConstituentSolvers => ListGetIntegerArray(Params, 'Block Solvers', Found)

    IsPlate = .FALSE.
    IsShell = .FALSE.
    IsNS = .FALSE.
    IsPres = .FALSE.
    
    i = ListGetInteger( Params,'Structure Solver Index',Found)

    IF ( Found ) THEN
      IsPlate = ListGetLogical( CurrentModel % Solvers(i) % Values,&
          'Plate Solver', Found )
      IsShell = ListGetLogical( CurrentModel % Solvers(i) % Values,&
          'Shell Solver', Found )      
    ELSE
      i = ListGetInteger( Params,'Plate Solver Index',IsPlate)
      IF(.NOT. IsPlate ) THEN
        i = ListGetInteger( Params,'Shell Solver Index',IsShell)
      END IF
    END IF
    
    ! The first and second entries in the "Block Solvers" list define 
    ! the solver sections to assemble the (1,1)-block and (2,2)-block, 
    ! respectively. The following check is needed as the solver section
    ! numbers may not index TotMatrix % Submatrix(:,:) directly.
    !
    IF (i > 0) THEN
      Found = .FALSE.
      DO k=1,SIZE(ConstituentSolvers)
        IF (i == ConstituentSolvers(k)) THEN
          Found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. Found) THEN
        CALL Fatal(Caller, &
            'Structure/Plate/Shell Solver Index is not an entry in the Block Solvers array')
      ELSE
        i = k
      END IF
    END IF
    IF (i > 2) CALL Fatal(Caller, &
        'Use the first two entries of Block Solvers to define FSI coupling')
      
    j = ListGetInteger( Params,'Fluid Solver Index',Found)
    IF(.NOT. Found ) THEN
      j = ListGetInteger( Params,'NS Solver Index', IsNs )
      IF( .NOT. IsNs ) THEN
        j = ListGetInteger( Params,'Pressure Solver Index',IsPres)
      END IF
    END IF

    IF (j > 0) THEN
      Found = .FALSE.
      DO k=1,SIZE(ConstituentSolvers)
        IF (j == ConstituentSolvers(k)) THEN
          Found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. Found) THEN
        CALL Fatal(Caller, &
            'Fluid/NS/Pressure Solver Index is not an entry in the Block Solvers array')
      ELSE
        j = k
      END IF
    END IF
    IF (j > 2) CALL Fatal(Caller, &
        'Use the first two entries of Block Solvers to define FSI coupling')

    IF( j == 0 ) THEN
      IF( i > 1 .AND. TotMatrix % NoVar == 2 ) j = 3 - i 
    END IF
    IF( i == 0 ) THEN
      IF( j > 1 .AND. TotMatrix % NoVar == 2 ) i = 3 - j
    END IF      
    
    IF(i<=0 .OR. j<=0) THEN
      IF( i > 0 ) CALL Warn(Caller,'Structure solver given but not fluid!')
      IF( j > 0 ) CALL Warn(Caller,'Fluid solver given but not structure!')
      RETURN
    END IF
    
!    IF (i > TotMatrix % NoVar .OR. j > TotMatrix % NoVar) &
!        CALL Fatal(Caller,'Use solver sections 1 and 2 to define FSI coupling') 
  
    A_fs => TotMatrix % Submatrix(j,i) % Mat
    A_sf => TotMatrix % Submatrix(i,j) % Mat
    
    IF(.NOT. ASSOCIATED( A_fs ) ) THEN
      CALL Fatal(Caller,'Fluid-structure coupling matrix not allocated!')
    END IF
    IF(.NOT. ASSOCIATED( A_sf ) ) THEN
      CALL Fatal(Caller,'Structure-fluid coupling matrix not allocated!')
    END IF
       
    SVar => TotMatrix % Subvector(i) % Var
    FVar => TotMatrix % Subvector(j) % Var

    A_s => TotMatrix % Submatrix(i,i) % Mat
    A_f => TotMatrix % Submatrix(j,j) % Mat
    
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL Fatal(Caller,'Fluid variable not present!')
    END IF
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL Fatal(Caller,'Structure variable not present!')
    END IF

    IF(.NOT. (IsNs .OR. IsPres ) ) THEN
      IsPres = ( FVar % Dofs <= 2 )
      IsNs = .NOT. IsPres
    END IF
    
    CALL FsiCouplingAssembly( Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
        IsPlate, IsShell, IsNS )

    IF( ParEnv % PEs > 1 ) THEN    
      TotMatrix % Submatrix(i,j) % ParallelSquareMatrix = .FALSE.
      TotMatrix % Submatrix(j,i) % ParallelSquareMatrix = .FALSE.
      
      TotMatrix % Submatrix(i,j) % ParallelIsolatedMatrix = &
          .NOT. CheckParallelCoupling(A_s, A_f, A_sf )  
      TotMatrix % Submatrix(j,i) % ParallelIsolatedMatrix = &
          .NOT. CheckParallelCoupling(A_f, A_s, A_fs )  
    END IF
      
       
  END SUBROUTINE FsiCouplingBlocks
    

  !> Create the coupling between elasticity solvers of various types.
  !--------------------------------------------------------------------------------
  SUBROUTINE StructureCouplingBlocks( Solver )

    TYPE(Solver_t) :: Solver
    
    INTEGER :: i,j,k,ind1,ind2,Novar,Nsol
    INTEGER, POINTER :: ConstituentSolvers(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params, ShellParams
    TYPE(Matrix_t), POINTER :: A_fs, A_sf, A_s, A_f
    TYPE(Variable_t), POINTER :: FVar, SVar
    LOGICAL :: IsPlate, IsShell, IsBeam, IsSolid, GotBlockSolvers
    LOGICAL :: DrillingDOFs
    TYPE(Solver_t), POINTER :: PSol
    CHARACTER(*), PARAMETER :: Caller = 'StructureCouplingBlocks'

    
    Params => Solver % Values
    ConstituentSolvers => ListGetIntegerArray(Params, 'Block Solvers', GotBlockSolvers)
    IF(.NOT. GotBlockSolvers ) THEN
      CALL Fatal(Caller,'We need "Block Solvers" defined!')
    END IF

    ! Currently we simply assume the master solver to be listed as the first entry in
    ! the 'Block Solvers' array.
    i = 1
    SVar => TotMatrix % Subvector(i) % Var
    IF(.NOT. ASSOCIATED( SVar ) ) THEN
      CALL Fatal(Caller,'Master structure variable not present!')
    END IF
    A_s => TotMatrix % Submatrix(i,i) % Mat
    
    Nsol = SIZE( ConstituentSolvers )

    
    DO j = 1, Nsol
      ! No need to couple to one self!
      IF(j==1) CYCLE
      IF (j > size(ConstituentSolvers)) CALL Fatal(Caller, &
          'Solid/Plate/Shell/Beam Solver Index larger than Block Solvers array')

      k = ConstituentSolvers(j)
      PSol => CurrentModel % Solvers(k)
      
      IsSolid = ListGetLogical( Psol % Values,'Solid Solver',IsSolid)
      IsPlate = ListGetLogical( Psol % Values,'Plate Solver',IsPlate)
      IsShell = ListGetLogical( Psol % Values,'Shell Solver',IsShell)
      IsBeam = ListGetLogical( Psol % Values,'Beam Solver',IsBeam)

      ind1 = ConstituentSolvers(i)
      ind2 = ConstituentSolvers(j)
      CALL Info(Caller,'Generating coupling between solvers '&
          //I2S(ind1)//' and '//I2S(ind2))

      
      A_fs => TotMatrix % Submatrix(j,i) % Mat
      A_sf => TotMatrix % Submatrix(i,j) % Mat
      
      !SVar => TotMatrix % Subvector(i) % Var
      FVar => TotMatrix % Subvector(j) % Var
      IF(.NOT. ASSOCIATED( FVar ) ) THEN
        CALL Fatal(Caller,'Slave structure variable not present!')
      END IF
      
      !A_s => TotMatrix % Submatrix(i,i) % Mat
      A_f => TotMatrix % Submatrix(j,j) % Mat

      IF (IsShell) THEN
        ShellParams => Fvar % Solver % Values
        DrillingDOFs = GetLogical(ShellParams, 'Drilling DOFs', Found)
      ELSE
        DrillingDOFs = .FALSE.
      END IF
      
      CALL StructureCouplingAssembly( Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
          IsSolid, IsPlate, IsShell, IsBeam, DrillingDOFs)
            
      IF( ParEnv % PEs > 1 ) THEN    
        TotMatrix % Submatrix(i,j) % ParallelSquareMatrix = .FALSE.
        TotMatrix % Submatrix(j,i) % ParallelSquareMatrix = .FALSE.
        
        TotMatrix % Submatrix(i,j) % ParallelIsolatedMatrix = &
            .NOT. CheckParallelCoupling(A_s, A_f, A_sf )  
        TotMatrix % Submatrix(j,i) % ParallelIsolatedMatrix = &
            .NOT. CheckParallelCoupling(A_f, A_s, A_fs )  
      END IF

    END DO
    
  END SUBROUTINE StructureCouplingBlocks
  

  ! This is tailored L2 norm for the many use types of the block solver.
  !---------------------------------------------------------------------
  FUNCTION CompNorm( x, n, npar, A) RESULT ( nrm ) 
    REAL(KIND=dp) :: x(:)
    INTEGER :: n
    INTEGER, OPTIONAL :: npar
    TYPE(Matrix_t), POINTER, OPTIONAL :: A
    REAL(KIND=dp) :: nrm

    INTEGER :: i,m
    REAL(KIND=dp) :: s,stot,ntot

    IF( ParEnv % PEs > 1 .AND. PRESENT( A ) ) THEN
      m = 0
      s = 0.0_dp
      DO i=1,A % NumberOfRows
        IF (Parenv % MyPE /= A % ParallelInfo % NeighbourList(i) % Neighbours(1)) CYCLE
        m = m+1
        s = s + x(i)**2
      END DO
    ELSE
      IF( ParEnv % PEs > 1 .AND. PRESENT( npar ) ) THEN
        m = npar
      ELSE
        m = n
      END IF
      s = SUM(x(1:m)**2)
    END IF
          
    stot = ParallelReduction(s)
    ntot = ParallelReduction(m)
    
    nrm = SQRT( stot / ntot )
    
  END FUNCTION CompNorm
  
  
  !------------------------------------------------------------------------------          
  !> Compute the rhs for the block matrix system which is solved
  !> accounting only the diagonal entries i.e. subtract the non-diagonal 
  !> matrix-vector results from the original r.h.s. vectors.
  !> After this the block diagonal problem Ax=b may be solved.
  !----------------------------------------------------------------------------------
  SUBROUTINE BlockUpdateRhs( BlockMatrix, ThisRow )
    
    TYPE(BlockMatrix_t), TARGET :: BlockMatrix
    INTEGER, OPTIONAL :: ThisRow
    
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: n, NoRow,NoCol, NoVar
    REAL(KIND=dp), POINTER :: x(:),rtmp(:),rhs(:)
    REAL(KIND=dp) :: bnorm
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: GotRhs, Trans
    
    CALL Info('BlockUpdateRhs','Computing block r.h.s',Level=8)

    NoVar = BlockMatrix % NoVar
    
    ! The residual is used only as a temporary vector
    ALLOCATE( rtmp(BlockMatrix % MaxSize) )
    
    
    DO NoRow = 1,NoVar 
      
      ! Optionally only one diagonal block may be updated for
      IF( PRESENT( ThisRow ) ) THEN
        IF( NoRow /= ThisRow ) CYCLE 
      END IF
      
      Var => BlockMatrix % SubVector(NoRow) % Var
      x => Var % Values
      n = SIZE( x )
      
      ! The r.h.s. of the initial system is stored in the Matrix
      !-----------------------------------------------------------
      IF(.NOT. ALLOCATED( BlockMatrix % SubVector(NoRow) % rhs )) THEN
        ALLOCATE( BlockMatrix % SubVector(NoRow) % rhs(n) )
        BlockMatrix % SubVector(NoRow) % rhs = 0.0_dp
        CALL Info('BlockUpdateRhs','Creating rhs for component: '//I2S(NoRow),Level=12)
      END IF
      rhs => BlockMatrix % SubVector(NoRow) % rhs
      
      A => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat
      GotRhs = .FALSE.
      IF( ASSOCIATED( A ) ) THEN
        IF( ASSOCIATED( A % rhs ) ) THEN
          GotRhs = .TRUE.
          rhs = A % rhs
        END IF
      END IF

      IF( .NOT. GotRhs ) rhs = 0.0_dp
      
      DO NoCol = 1,NoVar           
        ! This ensures that the diagonal itself is not subtracted
        ! before computing the bnorm used to estimate the convergence.
        IF( NoCol == NoRow ) CYCLE
        
        Var => BlockMatrix % SubVector(NoCol) % Var
        x => Var % Values

        Trans = .FALSE.        
        A => BlockMatrix % SubMatrix( NoRow, NoCol ) % Mat
        IF( A % NumberOfRows == 0 ) THEN
          IF( BlockMatrix % SubMatrixTranspose(NoCol,NoRow) ) THEN         
            A => TotMatrix % SubMatrix(NoCol,NoRow) % Mat
            Trans = .TRUE.
          END IF
        END IF
                       
        IF( A % NumberOfRows == 0 ) CYCLE

        IF(.NOT. Trans) THEN
          CALL CRS_MatrixVectorMultiply( A, x, rtmp)              
        ELSE 
          CALL CRS_TransposeMatrixVectorMultiply( A, x, rtmp )
        END IF
        rhs(1:n) = rhs(1:n) - rtmp(1:n) 
      END DO

      
      bnorm = CompNorm(rhs,n)
      BlockMatrix % SubVector(NoRow) % bnorm = bnorm
      
      ! Finally deduct the diagonal entry so that we can solve for the residual
      NoCol = NoRow
      Var => BlockMatrix % SubVector(NoCol) % Var
      x => Var % Values
      A => BlockMatrix % SubMatrix( NoRow, NoCol ) % Mat
      IF( A % NumberOfRows > 0 ) THEN
        CALL CRS_MatrixVectorMultiply( A, x, rtmp)              
        rhs(1:n) = rhs(1:n) - rtmp(1:n)
      END IF
      
    END DO
    
    DEALLOCATE( rtmp )
    
  END SUBROUTINE BlockUpdateRhs
  

  !------------------------------------------------------------------------------
  !> Perform matrix-vector product v=Au for block matrices.
  !> Has to be callable outside the module by Krylov methods.
  !------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixVectorProd( u,v,ipar )
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(in) :: u(*)
    REAL(KIND=dp), INTENT(out) :: v(*)
    INTEGER, INTENT(in) :: ipar(*)
    
    INTEGER :: n,i,j,k,NoVar,i1,i2,j1,j2,ll,kk
    REAL(KIND=dp), ALLOCATABLE :: s(:)
    INTEGER :: maxsize,ndofs
    INTEGER, POINTER :: Offset(:)
    REAL(KIND=dp) :: nrm
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:)

    LOGICAL :: Trans, Isolated
    LOGICAL :: DoSum , DoAMGXMV, Found

    DoAMGXMV = ListGetLogical( SolverRef % Values, 'Block AMGX M-V', Found)
    
    CALL Info('BlockMatrixVectorProd','Starting block matrix multiplication',Level=20)

    NoVar = TotMatrix % NoVar
    MaxSize = TotMatrix % MaxSize
    ALLOCATE( s(MaxSize) )

    IF( isParallel ) THEN
      Offset => TotMatrix % ParOffset
    ELSE
      Offset => TotMatrix % Offset
    END IF
    
    v(1:offset(NoVar+1)) = 0
    

    DO i=1,NoVar
      DO j=1,NoVar
        s = 0._dp

        j1 = offset(j)+1
        j2 = offset(j+1)

        Trans = .FALSE.
        Isolated = .FALSE.
        A => TotMatrix % SubMatrix(i,j) % Mat
        
        IF( A % NumberOfRows == 0 ) THEN
          IF( TotMatrix % SubMatrixTranspose(j,i) ) THEN
            A => TotMatrix % SubMatrix(j,i) % Mat
            IF( A % NumberOfRows >  0 ) THEN
              Trans = .TRUE.
              CALL Info('BlockMatrixVectorProd','Multiplying with transpose of submatrix ('&
                  //I2S(j)//','//I2S(i)//')',Level=10)
              IF( isParallel ) THEN
                Isolated = TotMatrix % SubMatrix(j,i) % ParallelIsolatedMatrix
                IF( .NOT. Isolated ) THEN
                  CALL Fatal('BlockMatrixVectorProd','Only isolated matrices may be transposed in parallel!')
                END IF
              END IF
            END IF
          END IF
        ELSE
          Isolated = TotMatrix % SubMatrix(i,j) % ParallelIsolatedMatrix 
        END IF
        IF( A % NumberOfRows == 0) CYCLE
        
        CALL Info('BlockMatrixVectorProd','Multiplying with submatrix ('&
            //I2S(i)//','//I2S(j)//')',Level=15)          
        
        IF (isParallel) THEN
          IF( ASSOCIATED( A % ParMatrix ) ) THEN
            CALL ParallelMatrixVector( A, u(j1:j2), s  )
          ELSE IF( Isolated ) THEN
            IF(.NOT. Trans ) THEN
              CALL CRS_MatrixVectorMultiply( A, u(j1:j2), s )
            ELSE
              CALL CRS_TransposeMatrixVectorMultiply( A, u(j1:j2), s )              
            END IF
          ELSE
            CALL Fatal('BlockMatrixVectorProd','Cannot make the matric-vector product in parallel!')
          END IF
        ELSE
          IF( .NOT. Trans ) THEN
            IF ( DoAMGXMV ) THEN
              CALL AMGXMatrixVectorMultiply(A, u(j1:j2), s, SolverRef )
            ELSE
              CALL CRS_MatrixVectorMultiply( A, u(j1:j2), s )
            END IF
          ELSE
            CALL CRS_TransposeMatrixVectorMultiply( A, u(j1:j2), s )
          END IF
        END IF
          
        IF( InfoActive( 25 ) ) THEN
          PRINT *,'MatVecProdNorm u:',i,j,&
              SQRT(SUM(u(j1:j2)**2)),SUM( u(j1:j2) ), MINVAL( u(j1:j2) ), MAXVAL( u(j1:j2) ) 
          PRINT *,'MatVecProdNorm s:',i,j,&
              SQRT(SUM(s**2)), SUM( s ), MINVAL( s ), MAXVAL( s ) 
        END IF

        v(offset(i)+1:offset(i+1)) = v(offset(i)+1:offset(i+1)) + s(1:offset(i+1)-offset(i))
      END DO
      
      IF( InfoActive( 25 ) ) THEN
        i1 = offset(i)+1
        i2 = offset(i+1)
        b => TotMatrix % Submatrix(i,i) % Mat % Rhs
        IF( ASSOCIATED( b ) ) THEN
          PRINT *,'MatVecProdNorm b:',i,&
              SQRT(SUM(b**2)),SUM( b ), MINVAL( b ), MAXVAL( b )
        END IF
        PRINT *,'MatVecProdNorm v:',i,&
            SQRT(SUM(v(i1:i2)**2)), SUM( v(i1:i2) ), MINVAL( v(i1:i2) ), MAXVAL( v(i1:i2) ) 
      END IF
    END DO
    
    IF( InfoActive( 25 ) ) THEN
      n = offset(NoVar+1)
      nrm = CompNorm(v(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'Mv result norm: ',nrm
      CALL Info('BlockMatrixVectorProd',Message )
    END IF
      
    CALL Info('BlockMatrixVectorProd','Finished block matrix multiplication',Level=20)
!------------------------------------------------------------------------------
  END SUBROUTINE BlockMatrixVectorProd
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> Given a permutation between monolithic and block matrix solutions
!> create a map between the owned dofs for the same in parallel. 
!------------------------------------------------------------------------------
  SUBROUTINE ParallelShrinkPerm()
    TYPE(Matrix_t), POINTER :: A, Adiag    
    INTEGER, POINTER :: BlockPerm(:), ParBlockPerm(:), ParPerm(:)    
    INTEGER :: i,j,k,l,n
    INTEGER, ALLOCATABLE :: ShrinkPerm(:),RenumPerm(:)

    n = TotMatrix % TotSize
           
    ! Example of blockperm: block solution u to monolithic solution v
    ! v(1:n) = u(BlockPerm(1:n)) 

    ! Dense numbering the monolithic dofs
    A => TotMatrix % ParentMatrix 
    IF(.NOt. ASSOCIATED(A) ) A => SolverMatrix

    IF( ASSOCIATED( A ) ) THEN
      CALL Info('ParallelShrinkPerm','Using the monolithic matrix to deduce permutation',Level=6)
    ELSE
      CALL Info('ParallelShrinkPerm','Using the block matrix to deduce permutation',Level=6)
    END IF

    ALLOCATE(ShrinkPerm(n))

    DO j=1,TotMatrix % Novar
      Adiag => TotMatrix % Submatrix(j,j) % Mat
      IF(.NOT. ASSOCIATED( Adiag % ParallelInfo ) ) CYCLE
      k = 0
      l = 0
      ShrinkPerm(1:Adiag % NumberOfRows) = 0 
      DO i=1,Adiag % NumberOfRows
        k = k + 1
        IF (Parenv % MyPE /= Adiag % ParallelInfo % NeighbourList(i) % Neighbours(1)) CYCLE
        l = l+1
        ShrinkPerm(k) = l
      END DO

      IF(.NOT. ASSOCIATED(TotMatrix % Submatrix(j,j) % ParPerm ) ) THEN
        ALLOCATE( TotMatrix % Submatrix(j,j) % ParPerm(l) )
      END IF
      ParPerm => TotMatrix % Submatrix(j,j) % ParPerm
      ParPerm = 0
      DO i=1,Adiag % NumberOfRows
        k = ShrinkPerm(i)
        IF(k==0) CYCLE
        ParPerm(k) = i
      END DO      
    END DO

    
    ShrinkPerm = 0    
    IF( ASSOCIATED( A ) ) THEN
      l = 0
      DO i=1,A % NumberOfRows
        IF (Parenv % MyPE /= A % ParallelInfo % NeighbourList(i) % Neighbours(1)) CYCLE
        l = l+1
        ShrinkPerm(i) = l
      END DO
    END IF
    
    IF(.NOT. ASSOCIATED(TotMatrix % ParPerm) ) THEN
      ALLOCATE( TotMatrix % ParPerm(l) )
    END IF
    ParPerm => TotMatrix % ParPerm
    ParPerm = 0
    DO i=1,n
      j = ShrinkPerm(i)
      IF(j==0) CYCLE
      ParPerm(j) = i
    END DO

    
    IF(.NOT. ASSOCIATED(TotMatrix % ParOffset) ) THEN
      ALLOCATE( TotMatrix % ParOffset(TotMatrix % NoVar+1))
    END IF
    TotMatrix % ParOffset = 0
    k = 0
    l = 0      
    DO j=1,TotMatrix % Novar
      A => TotMatrix % Submatrix(j,j) % Mat
      DO i=1,A % NumberOfRows
        k = k + 1
        IF (Parenv % MyPE /= A % ParallelInfo % NeighbourList(i) % Neighbours(1)) CYCLE
        l = l+1
      END DO
      TotMatrix % ParOffset(j+1) = l 
    END DO
      
    
    CALL Info('ParallelShrinkPerm','Number of parallel dofs in this partition: '// &
        I2S(l)//' / '//I2S(n), Level=6)    
    
    ! We can only make the ParBlockPerm if also BlockPerm exists!
    BlockPerm => TotMatrix % BlockPerm 
    IF(.NOT. ASSOCIATED(BlockPerm) ) RETURN

    ! Dense numbering for the block system dofs
    ALLOCATE(RenumPerm(n))
    RenumPerm = 1    
    DO i=1,n
      j = BlockPerm(i)
      IF( ShrinkPerm(j) == 0 ) RenumPerm(i) = 0
    END DO
    l = 0
    DO i=1,n
      IF(RenumPerm(i) > 0) THEN
        l=l+1
        RenumPerm(i) = l
      END IF
    END DO

    ! Allocate the parallel permutation, if not already done
    IF(.NOT. ASSOCIATED(TotMatrix % ParBlockPerm) ) THEN
      ALLOCATE(TotMatrix % ParBlockPerm(l))
    END IF
    ParBlockPerm => TotMatrix % ParBlockPerm
    ParBlockPerm = 0

    ! And finally create the renumbering scheme.
    ! This was tough to settle. Don't doubt it!
    DO i=1,n
      IF(RenumPerm(i)==0) CYCLE
      ParBlockPerm(RenumPerm(i)) = ShrinkPerm(BlockPerm(i))
    END DO
    
  END SUBROUTINE ParallelShrinkPerm

  
!------------------------------------------------------------------------------
! If the block matrix is just a remake of the monolithic one we may actually
! use the monolithic version when we make the necessary permutations.
! The advantage may be that we have an alternative (more simple) routine for
! debugging and it may also be faster...
! Note that here u and v only include the owned dofs for each partition.
! Hence if we want to make a standard matrix-vector product we need to expand
! these to include all dofs.
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixVectorProdMono( u,v,ipar )
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(in) :: u(*)
    REAL(KIND=dp), INTENT(out) :: v(*)
    INTEGER, INTENT(in) :: ipar(*)
    
    INTEGER :: n,m,i,j,k,NoVar,i1,i2,j1,j2
    REAL(KIND=dp), ALLOCATABLE :: s(:)
    INTEGER :: maxsize,ndofs
    INTEGER, POINTER :: Offset(:)
    REAL(KIND=dp) :: nrm
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:)
    LOGICAL :: DoSum, GotBlockStruct 
    REAL(KIND=dp), POINTER :: utmp(:),vtmp(:)
    INTEGER, POINTER :: BlockPerm(:)
    
    CALL Info('BlockMatrixVectorProdMono','Starting monolithic matrix multiplication',Level=20)

    IF(.NOT.ASSOCIATED(SolverMatrix)) THEN
      CALL Fatal('BlockMatrixVectorProdMono','No matrix to apply.')
    END IF
    
    NoVar = TotMatrix % NoVar
    MaxSize = TotMatrix % MaxSize
    GotBlockStruct = TotMatrix % GotBlockStruct

    n = ipar(3)
    
    IF(isParallel) THEN
      BlockPerm => TotMatrix % ParBlockPerm
      IF(.NOT. ASSOCIATED( BlockPerm ) ) THEN
        BlockPerm => TotMatrix % ParPerm 
      END IF
      IF(.NOT. ASSOCIATED(BlockPerm) ) THEN
        CALL Fatal('BlockMatrixVectorProdMono','How come there is no permutation in parallel?')
      END IF
      Offset => TotMatrix % ParOffset
    ELSE
      BlockPerm => TotMatrix % BlockPerm
      Offset => TotMatrix % Offset
    END IF

    IF( ASSOCIATED( BlockPerm ) ) THEN
      IF(InfoActive(20)) THEN
        i = MAXVAL( BlockPerm )
        j = MINVAL( BlockPerm )        
        IF( j <= 0 .OR. i > n ) THEN
          CALL Fatal('BlockMatrixVectorProdMono','Invalid sizes!')
        END IF
      END IF
    END IF


    IF(isParallel) THEN
      ! Reorder & copy block variable to monolithic variable
      m = SolverMatrix % NumberOfRows
      ALLOCATE(utmp(m),vtmp(m))
      utmp(1:m) = 0.0_dp
      utmp(BlockPerm(1:n)) = u(1:n)
      CALL ParallelMatrixVector(SolverMatrix, utmp(1:m), vtmp(1:m) )
      v(1:n) = vtmp(BlockPerm(1:n))
      DEALLOCATE(utmp,vtmp)
    ELSE
      IF( ASSOCIATED( BlockPerm ) ) THEN
        ! Only reorder between block and monolithic ordering
        ALLOCATE(utmp(n))
        utmp(BlockPerm(1:n)) = u(1:n)
        CALL CRS_MatrixVectorMultiply( SolverMatrix, utmp, v )
        v(1:n) = v(BlockPerm(1:n))              
      ELSE
        CALL CRS_MatrixVectorMultiply( SolverMatrix, u(1:n), v(1:n) )
      END IF
    END IF
               
    IF( InfoActive( 25 ) ) THEN
      nrm = CompNorm(v(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'Mv result norm: ',nrm
      CALL Info('BlockMatrixVectorProdMono',Message )
    END IF
      
    CALL Info('BlockMatrixVectorProdMono','Finished block matrix multiplication',Level=20)
!------------------------------------------------------------------------------
  END SUBROUTINE BlockMatrixVectorProdMono
!------------------------------------------------------------------------------

  

!> Create the vectors needed for block matrix scaling. Currently only
!> real and complex valued row equilibration is supported. Does not perform
!> the actual scaling.
!------------------------------------------------------------------------------
  SUBROUTINE CreateBlockMatrixScaling( )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: i,j,k,l,n,m,NoVar,istat
    REAL(KIND=dp) :: nrm, tmp, blocknrm
    TYPE(Matrix_t), POINTER :: A, Atrans
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    LOGICAL :: ComplexMatrix, GotIt, DiagOnly, PrecScale
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: SolverVar
    CHARACTER(*), PARAMETER :: Caller = 'CreateBlockMatrixScaling'
    
    CALL Info(Caller,'Starting block matrix row equilibration',Level=20)

    NoVar = TotMatrix % NoVar
    
    Params => CurrentModel % Solver % Values 
    DiagOnly = ListGetLogical( Params,'Block Scaling Diagonal',Found ) 
    IF( DiagOnly ) THEN
      CALL Info(Caller,'Considering only diagonal matrices in scaling',Level=20)      
    END IF

    PrecScale = ListGetLogical( Params,'Block Scaling PrecMatrix',Found ) 
    
    m = 0
    DO k=1,NoVar
      A => TotMatrix % SubMatrix(k,k) % Mat
      n = A % NumberOfRows
      IF(.NOT. ASSOCIATED(TotMatrix % Subvector )) THEN
        CALL Fatal(Caller,'Subvector not associated!')
      END IF
      IF( .NOT. ALLOCATED( Totmatrix % SubVector(k) % DiagScaling ) ) THEN
        m = m + 1
        ALLOCATE( TotMatrix % SubVector(k) % DiagScaling(n), STAT=istat )
        IF( istat /= 0 ) THEN
          CALL Fatal(Caller,'Cannot allocate scaling vectors '//I2S(k)//' of size: '//I2S(n))
        END IF
      END IF
    END DO
    IF( m > 0 ) THEN
      CALL Info(Caller,'Allocated '//I2S(m)//' scaling vectors for rhs!',Level=8)
    END IF
    

    blocknrm = 0.0_dp
    DO k=1,NoVar
      GotIt = .FALSE.
      A => TotMatrix % SubMatrix(k,k) % Mat
      IF( ASSOCIATED( A ) ) THEN
        IF( A % NumberOfRows > 0 ) THEN
          GotIt = .TRUE.
          ComplexMatrix = A % COMPLEX
        END IF
      END IF
      IF(.NOT. GotIt) CALL Warn(Caller,'Improve complex matrix detection!')
        
      IF( ComplexMatrix ) THEN
        m = 2
        CALL Info(Caller,'Assuming complex matrix block: '//I2S(k),Level=20)
      ELSE
        m = 1
        CALL Info(Caller,'Assuming real valued matrix block: '//I2S(k),Level=20)
      END IF     
      
      n = A % NumberOfRows 
      Diag => TotMatrix % SubVector(k) % DiagScaling
      Diag = 0.0_dp

      
      DO l=1,NoVar
        IF( DiagOnly ) THEN
          IF( k /= l ) CYCLE
        END IF
        
        Found = .FALSE.
        IF(k==l .AND. PrecScale ) THEN
          A => TotMatrix % Submatrix(k,k) % PrecMat
          Found = ( A % NumberOfRows > 0 )
        END IF

        IF(.NOT.Found) THEN
          A => TotMatrix % Submatrix(k,l) % Mat          
          IF( A % NumberOfRows == 0 ) THEN
            IF( TotMatrix % SubMatrixTranspose(l,k) ) THEN                    
              CALL Info(Caller,'Creating the transpose for real!',Level=12)
              Atrans => CRS_Transpose( TotMatrix % SubMatrix(l,k) % Mat ) 
              TotMatrix % Submatrix(k,l) % Mat => Atrans
              TotMatrix % SubMatrixTranspose(l,k) = .FALSE.
              A => Atrans
            END IF
          END IF
        END IF
        
        IF(.NOT. ASSOCIATED( A  ) ) CYCLE
        IF( A % NumberOfRows ==  0 ) CYCLE

        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values
              
        !---------------------------------------------
        ! Compute 1-norm of each row
        !---------------------------------------------
        DO i=1,n,m
          tmp = 0.0_dp

          IF( ComplexMatrix ) THEN
            DO j=Rows(i),Rows(i+1)-1,2
              tmp = tmp + SQRT( Values(j)**2 + Values(j+1)**2 )
            END DO
          ELSE
            DO j=Rows(i),Rows(i+1)-1        
              tmp = tmp + ABS(Values(j))          
            END DO
          END IF

          ! Compute the sum to the real component, scaling for imaginary will be the same
          Diag(i) = Diag(i) + tmp
        END DO
      END DO

      IF (ParEnv % PEs > 1) THEN
        A => TotMatrix % SubMatrix(k,k) % Mat      
        CALL ParallelSumVector(A, Diag)
      END IF
      
      nrm = MAXVAL(Diag(1:n))
      IF( ParEnv % PEs > 1 ) THEN
        nrm = ParallelReduction(nrm,2)
      END IF
      blocknrm = MAX(blocknrm,nrm)
      
      ! Define the actual scaling vector (for real component)
      DO i=1,n,m
        IF (Diag(i) > TINY( nrm ) ) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
      END DO

      ! Scaling of complex component
      IF( ComplexMatrix ) Diag(2::2) = Diag(1::2)
      
      WRITE( Message,'(A,ES12.5)') 'Unscaled matrix norm for block '//I2S(k)//': ', nrm    
      CALL Info(Caller, Message, Level=10 )      
    END DO

    WRITE( Message,'(A,ES12.5)') 'Unscaled matrix norm: ', blocknrm    
    CALL Info(Caller, Message, Level=7 )
    
  END SUBROUTINE CreateBlockMatrixScaling
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixInfo()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: i,j,k,l,n,m,NoVar
    REAL(KIND=dp) :: nrm, tmp, blocknrm
    TYPE(Matrix_t), POINTER :: A, Atrans
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    LOGICAL :: ComplexMatrix, GotIt, DiagOnly
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: Found
    
    
    CALL Info('BlockMatrixInfo','Showing some ranges of block matrix stuff',Level=10)
    
    NoVar = TotMatrix % NoVar
    m = 0
    
    PRINT *,'BlockInfo:',NoVar   
    
    DO k=1,NoVar
      DO l=1,NoVar
        
        A => TotMatrix % SubMatrix(k,l) % Mat
        IF( .NOT. ASSOCIATED( A ) ) CYCLE
        IF( A % NumberOfRows == 0 ) CYCLE
        
        n = TotMatrix % offset(k+1) - TotMatrix % offset(k)
        IF(isParallel ) THEN
          m = TotMatrix % ParOffset(k+1) -TotMatrix % ParOffset(k)
        END IF
          
        PRINT *,'BlockInfo:',k,l,A % NumberOfRows, n, m, A % COMPLEX

        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values

        PRINT *,'BlockInfo: A range',SUM( Values ), MINVAL( Values ), MAXVAL( Values ) 
      END DO
    END DO
    
  END SUBROUTINE BlockMatrixInfo
!------------------------------------------------------------------------------



  
!> Performs the actual forward or reverse scaling. Optionally the scaling may be
!> applied to only one matrix with an optional r.h.s. The idea is that for
!> block preconditioning we may revert to the original symmetric matrix but
!> still use the optimal row equilibration scaling for the block system. 
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixScaling( reverse, blockrow, blockcol, bext, SkipMatrixScale  )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: reverse
    INTEGER, OPTIONAL :: blockrow, blockcol
    REAL(KIND=dp), POINTER, OPTIONAL :: bext(:)
    LOGICAL, OPTIONAL :: SkipMatrixScale

    INTEGER :: i,j,k,l,n,m,NoVar
    REAL(KIND=dp) :: nrm, tmp
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: backscale
    CHARACTER(*), PARAMETER :: Caller = 'BlockMatrixScaling'
                  
    IF( PRESENT( Reverse ) ) THEN
      BackScale = Reverse
    ELSE
      BackScale = .FALSE.
    END IF
    IF (BackScale) THEN
      CALL Info(Caller,'Performing block matrix reverse row equilibration',Level=10)
    ELSE
      CALL Info(Caller,'Performing block matrix row equilibration',Level=10)
    END IF

    NoVar = TotMatrix % NoVar   
    DO k=1,NoVar
      
      IF( PRESENT( blockrow ) ) THEN
        IF( blockrow /= k ) CYCLE
      END IF
      
      Diag => TotMatrix % SubVector(k) % DiagScaling
      IF( .NOT. ASSOCIATED( Diag ) ) THEN
        CALL Fatal(Caller,'Diag for scaling not associated!')
      END IF
      n = SIZE(Diag) 
      
      IF( BackScale ) Diag = 1.0_dp / Diag 
            
      DO l=1,NoVar        

        ! If we use unscaled special preconditioning matrix we don't need to scale it
        IF( PRESENT( SkipMatrixScale ) ) THEN
          IF( SkipMatrixScale ) CYCLE
        END IF
        
        IF( PRESENT( blockcol ) ) THEN
          IF( blockcol /= l ) CYCLE
        END IF
        
        A => TotMatrix % SubMatrix(k,l) % Mat
        IF( A % NumberOfRows == 0 ) CYCLE
                          
        Rows   => A % Rows
        Cols   => A % Cols
        Values => A % Values
        
        DO i=1,n    
          DO j=Rows(i),Rows(i+1)-1
            Values(j) = Values(j) * Diag(i)
          END DO
        END DO

#if 0
        ! This does not seem to be necessary but actually harmfull.
        A => TotMatrix % SubMatrix(k,l) % PrecMat
        IF( A % NumberOfRows == 0 ) CYCLE
        DO i=1,n    
          DO j=A % Rows(i),A % Rows(i+1)-1
            A % Values(j) = A % Values(j) * Diag(i)
          END DO
        END DO
#endif
      END DO
        
      IF( PRESENT( bext ) ) THEN
        b => bext
      ELSE        
        b => TotMatrix % Submatrix(k,k) % Mat % Rhs
      END IF
      
      IF( ASSOCIATED( b ) ) THEN
        b(1:n) = Diag(1:n) * b(1:n)
      END IF
      
      IF( BackScale ) Diag = 1.0_dp / Diag       
    END DO

    IF( BackScale ) THEN
      CALL Info(Caller,'Finished block matrix reverse row equilibration',Level=10)           
    ELSE
      CALL Info(Caller,'Finished block matrix row equilibration',Level=10)           
    END IF
      
  END SUBROUTINE BlockMatrixScaling
!------------------------------------------------------------------------------


!> Deallocates the block matrix scaling vectors.   
!------------------------------------------------------------------------------
  SUBROUTINE DestroyBlockMatrixScaling()
!------------------------------------------------------------------------------
    INTEGER :: k,NoVar
    
    CALL Info('DestroyBlockMatrixScaling','Deallocating the vectors for block system scaling',Level=10)
              
    NoVar = TotMatrix % NoVar   
    DO k=1,NoVar            
      IF( ALLOCATED( TotMatrix % SubVector(k) % DiagScaling ) ) THEN
        DEALLOCATE( TotMatrix % SubVector(k) % DiagScaling )
      END IF
    END DO

  END SUBROUTINE DestroyBlockMatrixScaling
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Perform block preconditioning for Au=v by solving all the individual diagonal problems.
!> Has to be called outside the module by Krylov methods.
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixPrec( u,v,ipar )    
    IMPLICIT NONE
    REAL(KIND=dp), TARGET, INTENT(out) :: u(*)
    REAL(KIND=dp), TARGET, INTENT(in) :: v(*)
    INTEGER :: ipar(*)
!---------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: rtmp(:),vtmp(:),xtmp(:),btmp(:),diagtmp(:),b(:),x(:), a_rhs_save(:)
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:)
    INTEGER :: i,j,k,l,n,m,NoVar,nc,kk,ll
    TYPE(Solver_t), POINTER :: Solver, Solver_save, ASolver
    INTEGER, POINTER :: Offset(:), BlockPerm(:),ParPerm(:)
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER :: BlockOrder(:)
    TYPE(Matrix_t), POINTER :: A, Aij, mat_save
    TYPE(Variable_t), POINTER :: Var, Var_save
    REAL(KIND=dp) :: nrm
    LOGICAL :: GotOrder, BlockGS, Found, NS, ScaleSystem, DoSum, &
        IsComplex, BlockScaling, DoDiagScaling, DoPrecScaling, UsePrecMat, Trans, &
        Isolated, NoNestedScaling, DoAMGXmv, CalcLoads
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc

    CALL Info('BlockMatrixPrec','Starting block matrix preconditioning',Level=8)

    DoAMGXMV = ListGetLogical( SolverRef % Values, 'Block AMGX M-V', Found)
    
    n = ipar(3)
    
    IF( InfoActive(25) ) THEN
      nrm = CompNorm(v(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'V start norm: ',nrm
      CALL Info('BlockMatrixPrec',Message,Level=10)
    END IF
      
    Solver => CurrentModel % Solver
    Params => Solver % Values
    
    ! Enable user defined order for the solution of blocks
    !---------------------------------------------------------------
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotOrder)
    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    
    NoVar = TotMatrix % NoVar
    Solver => TotMatrix % Solver

    TotMatrix % NoIters = TotMatrix % NoIters + 1
   
    IF( isParallel ) THEN
      offset => TotMatrix % ParOffset
    ELSE
      offset => TotMatrix % Offset
    END IF
      
    IF( n /= Offset(NoVar+1) ) THEN
      CALL Fatal('BlockMatrixPrec','There is a mismatch between sizes!')
    END IF
    
    ! Save the initial solver stuff
    solver_save => Solver
    var_save => Solver % Variable
    mat_save => Solver % Matrix
    rhs_save => Solver % Matrix % RHS

    BlockScaling = ListGetLogical( Params,'Block Scaling',Found )
    
    DoDiagScaling = .FALSE.
    IF( ASSOCIATED( Solver % Matrix ) ) THEN
      DoDiagScaling = ASSOCIATED( Solver % Matrix % diagscaling ) 
    END IF
    IF( DoDiagScaling ) THEN
      CALL Info('BlockMatrixPrec','External diagonal scaling is active!',Level=20)
    ELSE
      CALL Info('BlockMatrixPrec','External diagonal scaling is not active!',Level=30)
    END IF
    IF( BlockScaling ) THEN
      CALL Info('BlockMatrixPrec','Block matrix scaling is active!',Level=20)
    ELSE
      CALL Info('BlockMatrixPrec','Block matrix scaling is not active!',Level=30)
    END IF
      
    IF(DoDiagscaling) THEN
      NoNestedScaling = ListGetLogical( Params,'Eliminate Nested Scaling',Found )
      IF(.NOT. Found) NoNestedScaling = .TRUE.      
    ELSE
      NoNestedScaling = .FALSE.
    END IF

    ! Always treat the inner iterations as truly complex if they are
    CALL ListAddLogical( Params,'Linear System Skip Complex',.FALSE.) 
    
    IF (isParallel) THEN
      ALLOCATE( x(TotMatrix % MaxSize), b(TotMatrix % MaxSize) )
    END IF

    ! Initial guess:
    !-----------------------------------------
    u(1:n) = 0.0_dp
    
    IF( BlockGS ) THEN
      ALLOCATE( vtmp(n), rtmp(n), xtmp(n))
      vtmp(1:n) = v(1:n)
    END IF
    
    CALL ListPushNameSpace('block:')
    
    DO j=1,NoVar
      IF( GotOrder ) THEN
        i = BlockOrder(j)
      ELSE
        i = j
      END IF
      
      WRITE(Message,'(A,I0)') 'Solving block: ',i
      CALL Info('BlockMatrixPrec',Message,Level=8)

      CALL ListPushNameSpace('block '//i2s(i)//i2s(i)//':')
      
      ! Set pointers to the new linear system
      !-------------------------------------------------------------------
      Var => TotMatrix % SubVector(i) % Var

      UsePrecMat = .FALSE.
      IF( ASSOCIATED( TotMatrix % Subvector(i) % Solver ) ) THEN
        ASolver => TotMatrix % SubVector(i) % Solver
        A => ASolver % Matrix
      ELSE
        A => TotMatrix % Submatrix(i,i) % PrecMat
        IF( A % NumberOfRows == 0 ) THEN
          A => TotMatrix % Submatrix(i,i) % Mat
        ELSE
          UsePrecMat = .TRUE.
          CALL Info('BlockMatrixPrec','Using specialized (Schur) preconditioning block',Level=9)
        END IF      
        ASolver => Solver
      END IF      
        
      IF (isParallel) THEN
        ! copy part of full solution to block solution
        x = 0.0_dp
        b = 0.0_dp
        ParPerm => TotMatrix % Submatrix(i,i) % ParPerm
        x(ParPerm) = u(offset(i)+1:offset(i+1))
        IF( BlockGS ) THEN
          b(ParPerm) = vtmp(offset(i)+1:offset(i+1))
        ELSE
          b(ParPerm) = v(offset(i)+1:offset(i+1))
        END IF
      ELSE
        x => u(offset(i)+1:offset(i+1))
        IF( BlockGS ) THEN
          b => vtmp(offset(i)+1:offset(i+1))
        ELSE      
          b => v(offset(i)+1:offset(i+1))
        END IF
      END IF

      IF( InfoActive(25) ) THEN
        ! l is uninitialized!
        !nrm = CompNorm(b,offset(i+1)-offset(i),npar=l)
        nrm = CompNorm(b,offset(i+1)-offset(i))
        WRITE( Message,'(A,ES12.5)') 'Rhs '//I2S(i)//' norm: ',nrm
        CALL Info('BlockMatrixPrec',Message,Level=10)
      END IF
        
      ! Reuse block preconditioner from the first block to other components
      !--------------------------------------------------------------------
      IF( ListGetLogical( Params,'Block Prec Reuse',Found) ) THEN
        DO k = 1, NoVar
          IF( k == i ) CYCLE
          IF( CRS_CopyMatrixPrec( TotMatrix % Submatrix(k,k) % Mat, A ) ) EXIT
        END DO
      END IF
      
      ! We do probably not want to compute the change within each iteration
      CALL ListAddLogical( Asolver % Values,'Skip Advance Nonlinear iter',.TRUE.)         
      CALL ListAddLogical( Asolver % Values,'Skip Compute Nonlinear Change',.TRUE.)         
       
      IF( BlockScaling ) CALL BlockMatrixScaling(.TRUE.,i,i,b,UsePrecMat)

      ! The special preconditioning matrices have not been scaled with the monolithic system.
      ! So we need to transfer the (x,b) of this block to the unscaled system before going
      ! going to solve it. It is probably desirable to use separate scaling for this system. 
      DoPrecScaling = DoDiagScaling .AND. UsePrecMat
      IF( DoPrecScaling ) THEN
        n = A % NumberOfRows
        ALLOCATE( btmp(n), diagtmp(n) )

        IF( TotMatrix % GotBlockStruct ) THEN
          k = TotMatrix % InvBlockStruct(i)
          IF( k <= 0 ) THEN
            CALL Fatal('BlockMatrixPrec','Cannot define the originating block '&
                //I2S(i)//' for scaling!')
          ELSE
            CALL Info('BlockMatrixPrec','Using initial block '//I2S(k)//&
                ' in scaling of prec matrix '//I2S(i),Level=12)
          END IF
        ELSE
          k = i
        END IF

        l = Solver % Variable % DOFs
        diagtmp(1:n) = Solver % Matrix % DiagScaling(k::l)

        ! Scale x & b to the unscaled system of the tailored preconditioning matrix for given block.
        x(1:n) = x(1:n) * diagtmp(1:n)
        btmp(1:n) = b(1:n) / diagtmp(1:n) * Solver % Matrix % RhsScaling**2
      ELSE
        btmp => b
        IF( NoNestedScaling ) THEN
          CALL Info('BlockMatrixPrec','Eliminating scaling for block as outer scaling already done!',Level=25)
          CALL ListAddLogical( Params,'Linear System Skip Scaling',.TRUE.)
        END IF
      END IF
              

      IF( InfoActive( 25 ) ) THEN
        CALL BlockMatrixInfo()
      END IF

      IF( A % COMPLEX ) THEN
        nc = 2
      ELSE
        nc = 1
      END IF
      

      IF(DoAMGXMv) THEN
        ScaleSystem = ListGetLogical( Params,'Linear System Scaling', Found )
        IF(.NOT. Found) ScaleSystem = .TRUE.
        IF ( ScaleSystem ) CALL ScaleLinearSystem(ASolver, A,btmp,x )
        CALL AMGXSolver( A, x, btmp, ASolver )
        IF( ScaleSystem ) CALL BackScaleLinearSystem(ASolver,A,btmp,x)
      ELSE
        CalcLoads = ListGetLogical( ASolver % Values, 'Calculate Loads', Found )
        CALL ListAddLogical( ASolver % Values, 'Calculate Loads', .FALSE.)
        CALL SolveLinearSystem( A, btmp, x, Var % Norm, Var % DOFs, ASolver )
        IF (CalcLoads) CALL ListAddLogical( ASolver % Values, 'Calculate Loads', .TRUE.)        
      END IF

      ! If this was a special preconditioning matrix then update the solution in the scaled system. 
      IF( DoPrecScaling ) THEN
        x(1:n) = x(1:n) / diagtmp(1:n)
        DEALLOCATE( btmp, diagtmp )
      ELSE IF( NoNestedScaling ) THEN
        CALL ListAddLogical( Params,'Linear System Skip Scaling',.FALSE.)
      END IF
        
      IF( InfoActive(20) ) THEN
        nrm = CompNorm(x,offset(i+1)-offset(i),A=A)
        WRITE( Message,'(A,ES12.5)') 'Linear system '//I2S(i)//' norm: ',nrm
        CALL Info('BlockMatrixPrec',Message)
      END IF
        
      IF( BlockScaling ) CALL BlockMatrixScaling(.FALSE.,i,i,b,UsePrecMat)

      IF (isParallel) THEN
        x(1:offset(i+1)-offset(i)) = x(ParPerm) 
        u(offset(i)+1:offset(i+1)) = x(1:offset(i+1)-offset(i))
      END IF
 
      !---------------------------------------------------------------------
      IF( BlockGS ) THEN        
        CALL Info('BlockMatrixPrec','Updating block r.h.s',Level=9)
      
        DO l=j+1,NoVar
          IF( GotOrder ) THEN
            k = BlockOrder(l)
          ELSE
            k = l
          END IF

          str = 'Block Gauss-Seidel Passive '//I2S(k)//I2S(i)
          IF( ListGetLogical( Params, str, Found ) ) CYCLE

          CALL Info('BlockMatrixPrec','Updating r.h.s for component '//I2S(k),Level=15)
          
          ! The residual is used only as a temporary vector
          !-------------------------------------------------------------
          Trans = .FALSE.
          Isolated = .FALSE.          
          Aij => TotMatrix % SubMatrix(k,i) % Mat 

          ! We may use transpose that is actually never created as only its action is used!
          IF( Aij % NumberOfRows == 0 ) THEN
            IF( TotMatrix % SubMatrixTranspose(i,k) ) THEN
              Aij => TotMatrix % SubMatrix(i,k) % Mat
              IF( Aij % NumberOfRows >  0 ) THEN
                Trans = .TRUE.
                CALL Info('BlockMatrixPrec','Multiplying with transpose of submatrix ('&
                    //I2S(i)//','//I2S(k)//')',Level=8)
                IF(isParallel ) THEN
                  IF( TotMatrix % SubMatrix(i,k) % ParallelIsolatedMatrix ) THEN
                    Isolated = .TRUE.
                  ELSE
                    CALL Fatal('BlockMatrixPrec','Only isolated matrices may be transposed in parallel!')
                  END IF
                END IF
              END IF
            END IF
          ELSE
            Isolated = TotMatrix % SubMatrix(k,i) % ParallelIsolatedMatrix 
          END IF
          IF( Aij % NumberOfRows == 0) CYCLE
          

          IF (isParallel) THEN
            IF(ASSOCIATED(Aij % ParMatrix)) THEN
              ! x is packed, r is full 
              CALL ParallelMatrixVector(Aij,x,rtmp )
              
            ELSE IF( Isolated ) THEN
              ! If our matrix is not active on shared nodes we may apply serial Mv
              ! and pack the results to include only the dofs owned by the partition.
              IF( .NOT. Trans ) THEN
                CALL CRS_MatrixVectorMultiply(Aij,x,rtmp)
              ELSE
                CALL CRS_TransposeMatrixVectorMultiply( Aij, x, rtmp )
              END IF

              ParPerm => TotMatrix % Submatrix(i,i) % ParPerm

#if 0
              rtmp(1:offset(i+1)-offset(i)) = rtmp(ParPerm)
#else
              ll = 0
              DO kk=1,A % NumberofRows
                IF (Parenv % MyPE /= A % ParallelInfo % NeighbourList(kk) % Neighbours(1)) CYCLE
                ll = ll+1
                rtmp(ll) = rtmp(kk)
                IF(parperm(ll) /= kk) PRINT *,'Problem:',ll,kk,parperm(ll)
              END DO
#endif
              
            ELSE IF (ASSOCIATED(SolverMatrix)) THEN
              ! Here we don't have the luxury that the block matrix would either have parallel
              ! communication initiated, or not have interface dofs. The last resort is to use
              ! the initial monolithic matrix to perform Mv also for a given block by setting
              ! other dofs to zero. 
              xtmp = 0.0_dp              
              xtmp(offset(i)+1:offset(i+1)) = x(offset(i)+1:offset(i+1))

              BlockPerm => TotMatrix % ParBlockPerm

              xtmp(BlockPerm(1:n)) = xtmp(1:n)              
              CALL ParallelMatrixVector(SolverMatrix,xtmp,rtmp)
              rtmp(1:n) = rtmp(BlockPerm(1:n))

              rtmp(1:offset(k+1)-offset(k)) = rtmp(offset(k)+1:offset(k+1))
              
            ELSE
              CALL Fatal('BlockMatrixPrec','Do not know how to apply parallel matrix!')
            END IF
          ELSE
            Aij => TotMatrix % Submatrix(k,i) % Mat            
            IF( .NOT. Trans ) THEN
                IF(DoAMGXMV) THEN
                  CALL AMGXMatrixVectorMultiply(Aij, x, rtmp, SolverRef )
                ELSE
                  CALL CRS_MatrixVectorMultiply(Aij,x,rtmp )
                END IF
            ELSE
              CALL CRS_TransposeMatrixVectorMultiply( Aij, x, rtmp )
            END IF
          END IF

          ! Up-date the off-diagonal entries to the r.h.s. 
          vtmp(offset(k)+1:offset(k+1)) = vtmp(offset(k)+1:offset(k+1)) &
              - rtmp(1:offset(k+1)-offset(k))

        END DO ! l=j+1,NoVar
        
      END IF

      CALL ListPopNameSpace() ! block ij:
      
    END DO ! j=1,NoVar

    IF (isParallel) DEALLOCATE(x,b)

    CALL ListPopNameSpace('block:') ! block:

    CALL ListAddLogical( Params,'Linear System Refactorize',.FALSE. )
    CALL ListAddLogical( Asolver % Values,'Skip Advance Nonlinear iter',.FALSE.)
    CALL ListAddLogical( Asolver % Values,'Skip Compute Nonlinear Change',.FALSE.)
    
    Solver => Solver_save
    Solver % Matrix => mat_save
    Solver % Matrix % RHS => rhs_save
    Solver % Variable => Var_save

    IF( BlockGS ) THEN
      DEALLOCATE( vtmp, rtmp, xtmp ) 
    END IF

    IF( InfoActive(20) ) THEN
      nrm = CompNorm(v(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'V fin norm: ',nrm
      CALL Info('BlockMatrixPrec',Message,Level=10)
      
      nrm = CompNorm(u(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'U fin norm: ',nrm
      CALL Info('BlockMatrixPrec',Message,Level=10)
    END IF
      
    CALL Info('BlockMatrixPrec','Finished block matrix preconditioning',Level=8)
    
  END SUBROUTINE BlockMatrixPrec



  !> This call takes care of Jacobi & Gauss Seidel block methods. 
  !-----------------------------------------------------------------
  SUBROUTINE BlockStandardIter( Solver, MaxChange )

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: MaxChange

    INTEGER :: i,j,NoVar,RowVar,iter,LinIter,MinIter
    INTEGER, POINTER :: BlockOrder(:)
    LOGICAL :: GotIt, GotBlockOrder, BlockGS
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:), b(:)
    REAL(KIND=dp), POINTER :: dx(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, SolverVar
    REAL(KIND=dp) :: LinTol, TotNorm, dxnorm, xnorm, Relax
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ScaleSystem, BlockScaling, Found

    NoVar = TotMatrix % NoVar
    Params => Solver % Values
    SolverVar => Solver % Variable

    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotBlockOrder)
    LinIter = ListGetInteger( Params,'Linear System Max Iterations',GotIt)
    MinIter = ListGetInteger( Params,'Linear System Min Iterations',GotIt)
    LinTol = ListGetConstReal( Params,'Linear System Convergence Tolerance',GotIt)
    BlockScaling = ListGetLogical( Params,'Block Scaling',GotIt)
    
    CALL ListPushNamespace('block:')

    ! We don't want compute change externally
    CALL ListAddNewLogical( Params,'Skip compute nonlinear change',.TRUE.)
    
    Relax = 1.0_dp
    
    DO iter = 1, LinIter

      ! Store the iteration count
      TotMatrix % NoIters = iter
      
      ! In block Jacobi the r.h.s. is not updated during the iteration cycle
      !----------------------------------------------------------------------
      IF( BlockGS ) THEN
        WRITE( Message,'(A,I0)') 'Block Gauss-Seidel iteration: ',iter
      ELSE
        WRITE( Message,'(A,I0)') 'Block Jacobi iteration: ',iter
        CALL BlockUpdateRhs(TotMatrix)
      END IF
      CALL Info('BlockStandardIter',Message,Level=5)
      MaxChange = 0.0_dp
      TotNorm = 0.0_dp
      
      IF( iter == 2 ) THEN
        CALL ListAddLogical( Params,'No Precondition Recompute',.TRUE.)
      END IF
      
      DO i=1,NoVar
        IF( GotBlockOrder ) THEN
          RowVar = BlockOrder(i)
        ELSE
          RowVar = i
        END IF
        
        ! In gauss-seidel the partial update is immediately taken into account
        !---------------------------------------------------------------------
        IF( BlockGS ) THEN
          CALL BlockUpdateRhs(TotMatrix,RowVar)
        END IF
        
        IF( ListGetLogical( Params,'Block Prec Reuse',GotIt) ) THEN
          DO j = 1, NoVar
            IF( j == RowVar ) CYCLE
            IF( CRS_CopyMatrixPrec( TotMatrix % Submatrix(j,j) % Mat, A ) ) EXIT
          END DO
        END IF
        
        b => TotMatrix % SubVector(RowVar) % rhs
        
        IF( InfoActive( 15 ) ) THEN
          PRINT *,'rhs'//I2S(i)//':',SQRT( SUM(b**2) ), MINVAL( b ), MAXVAL( b ), SUM( b )
        END IF

        Var => TotMatrix % SubVector(RowVar) % Var
        Solver % Variable => Var
        
        A => TotMatrix % Submatrix(i,i) % PrecMat
        IF( A % NumberOfRows == 0 ) THEN
          A => TotMatrix % Submatrix(i,i) % Mat
        ELSE
          CALL Info('BlockStandardIter','Using preconditioning block: '//I2S(i),Level=8)
        END IF
        
        !Solver % Matrix => A

        ! Use the newly computed residual rather than original r.h.s. to solve the equation!!
        rhs_save => A % rhs ! Solver % Matrix % RHS
        A % RHS => b
        
        ! Solving the subsystem
        !-----------------------------------
        ALLOCATE( dx( SIZE( Var % Values ) ) )
        dx = 0.0_dp

        CALL ListPushNamespace('block '//i2s(RowVar)//i2s(RowVar)//':')          

        IF( BlockScaling ) CALL BlockMatrixScaling(.TRUE.,i,i,b)
              
        !IF( ListGetLogical( Solver % Values,'Linear System Complex', Found ) ) A % Complex = .TRUE.

        !ScaleSystem = ListGetLogical( Solver % Values,'block: Linear System Scaling', Found )
        !IF(.NOT. Found) ScaleSystem = .TRUE.

        
        CALL SolveLinearSystem( A, b, dx, Var % Norm, Var % DOFs, Solver )

        IF( BlockScaling ) CALL BlockMatrixScaling(.FALSE.,i,i,b)

        CALL ListPopNamespace()

        ! Revert back to original r.h.s.
        A % RHS => rhs_save
        !Solver % Matrix => mat_save

        IF( iter > 1 ) THEN
          Var % Values = Var % Values + Relax * dx
        ELSE
          Var % Values = Var % Values + dx
        END IF

        dxnorm = CompNorm(dx,A % NumberOfRows, A=A)
        xnorm = CompNorm(Var % Values, A % NumberOfRows, A=A)

        Var % Norm = xnorm
        Var % NonlinChange = dxnorm / xnorm

        WRITE(Message,'(A,2ES12.3)') 'Block '//I2S(RowVar)//' norms: ',xnorm, dxnorm / xnorm
        CALL Info('BlockStandardIter',Message,Level=5)
        
        IF( InfoActive( 20 ) ) THEN
          PRINT *,'dx'//I2S(i)//':',SQRT( SUM(dx**2) ), MINVAL( dx ), MAXVAL( dx ), SUM( dx ), SUM( ABS( dx ) )
        END IF
      
        DEALLOCATE( dx )
          
        TotNorm = TotNorm + Var % Norm
        MaxChange = MAX( MaxChange, Var % NonlinChange )        
      END DO

      WRITE(Message,'(A,2ES12.3)') 'Sum of norms: ',TotNorm, MaxChange
      CALL Info('BlockStandardIter',Message,Level=4)

      IF( MaxChange < LinTol .AND. iter >= MinIter ) THEN
        CALL Info('BlockStandardIter','Converged after iterations: '//I2S(iter),Level=5)
        EXIT
      END IF
      
    END DO
    CALL ListPopNamespace('block:')

    CALL ListAddLogical( Params,'No Precondition Recompute',.FALSE.)
        
    Solver % Variable => SolverVar

  END SUBROUTINE BlockStandardIter


  !---------------------------------------------------------------------------
  !> This call takes care of the iterative Krylov methods for block systems
  !> which can still be preconditioned by block Jacobi or Gauss-Seidel methods
  !---------------------------------------------------------------------------
  SUBROUTINE BlockKrylovIter( Solver, MaxChange )

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: MaxChange

    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc
    INTEGER(KIND=AddrInt) :: iterProc,precProc, mvProc,dotProc,nmrProc, zero=0
    REAL(KIND=dp) :: dpar(20), xnorm,prevxnorm
    REAL(KIND=dp), ALLOCATABLE :: x(:),b(:),r(:)
    
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: SVar
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoVar, ndim, maxsize
    LOGICAL :: Converged, Diverged
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    INTEGER, POINTER :: Offset(:),poffset(:),BlockStruct(:),ParPerm(:)
    INTEGER :: i,j,k,l,ia,ib,istat
    LOGICAL :: LS, BlockAV,Found, UseMono
    CHARACTER(*), PARAMETER :: Caller = 'BlockKrylovIter'
 
    
    CALL Info(Caller,'Starting block system iteration',Level=8)
    
    !CALL ListPushNameSpace('outer:')
    Params => Solver % Values
    
    BlockAV = ListGetLogical(Params,'Block A-V System', Found)

    ndim = TotMatrix % TotSize 
    NoVar = TotMatrix % NoVar

    TotMatrix % NoIters = 0
    
    offset => TotMatrix % Offset
    IF(isParallel) THEN
      poffset => TotMatrix % ParOffset
    ELSE
      poffset => offset
    END IF

    ! Just do some error checks
    DO i=1,NoVar
      IF( .NOT. ASSOCIATED( TotMatrix % Subvector(i) % Var ) ) THEN
        CALL Fatal(Caller,'Subvector '//I2S(i)//' not associated!')
      END IF
      A => TotMatrix % SubMatrix(i,i) % Mat
      IF( .NOT. ASSOCIATED( A ) ) THEN
        CALL Fatal(Caller,'Submatrix '//I2S(11*i)//' not associated!')
      END IF
      IF( .NOT. ASSOCIATED( A % Rhs ) ) THEN
        CALL Warn(Caller,'Submatrix rhs '//I2S(11*i)//' not associated!')
      END IF
    END DO
          
    CALL Info(Caller,'Allocating temporal vectors for block system of size: '&
        //I2S(ndim),Level=15)

    ALLOCATE(x(ndim), b(ndim),r(ndim),STAT=istat)
    IF( istat /= 0 ) THEN
      CALL Fatal(Caller,'Cannot allocate temporal vectors of size: '//I2S(ndim))
    END IF
    
    x = 0.0_dp
    b = 0.0_dp
    r = 0.0_dp
    
    IF (isParallel) THEN
      CALL Info(Caller,'Performing parallel initializations!',Level=18)
      DO i=1,NoVar
        DO j=1,NoVar
          A => TotMatrix % SubMatrix(i,j) % Mat          
          ! ParallelInitSolve expects full vectors
          IF ( i /= j ) THEN
            IF(ASSOCIATED(A % ParMatrix)) CALL ParallelInitSolve(A,r,r,r)
          ELSE
            IF (ASSOCIATED(A % ParMatrix % SplittedMatrix % InsideMatrix % PrecValues)) THEN
              IF (.NOT. ASSOCIATED(A % PrecValues)) & 
                  NULLIFY(A % ParMatrix % SplittedMatrix % InsideMatrix % PrecValues)
            END IF
            CALL ParallelInitSolve(A, TotMatrix % Subvector(i) % Var % Values, A % rhs, r )
            IF( ASSOCIATED(SolverMatrix)) THEN
              x(offset(i)+1:offset(i+1)) = TotMatrix % SubVector(i) % Var % Values        
              IF(ASSOCIATED(A % rhs)) b(offset(i)+1:offset(i+1)) = A % rhs
            END IF
          END IF
        END DO
      END DO
      IF(ASSOCIATED(SolverMatrix)) THEN
        CALL ParallelInitSolve( SolverMatrix, x, b, r )      
        x = 0.0_dp
        b = 0.0_dp
      END IF
    END IF
      
    CALL Info(Caller,'Initializing monolithic system vectors',Level=18)
    
    DO i=1,NoVar
      A => TotMatrix % SubMatrix(i,i) % Mat

      IF (.NOT.isParallel) THEN
        x(offset(i)+1:offset(i+1)) = TotMatrix % SubVector(i) % Var % Values        
        IF(ASSOCIATED(A % rhs)) b(offset(i)+1:offset(i+1)) = A % rhs
      ELSE 
        ParPerm => TotMatrix % SubMatrix(i,i) % ParPerm
        x(poffset(i)+1:poffset(i+1)) = TotMatrix % SubVector(i) % Var % Values(ParPerm)        

        ! This is a little dirty as it uses internal stuff from the parallel structure directly.
        ! However, only this r.h.s. vector seems to be up-to-date.
        IF(ASSOCIATED(A % rhs)) b(poffset(i)+1:poffset(i+1)) = A % rhs(ParPerm) !A % ParMatrix % SplittedMatrix % InsideMatrix % Rhs
      END IF
    END DO
    
    ! Parallel block system only solves for its own variables.
    IF (isParallel) THEN
      ndim = poffset(NoVar+1)
    END IF
    
    !----------------------------------------------------------------------
    ! Solve matrix equation solver with the redefined block matrix operations
    !----------------------------------------------------------------------
    CALL ListAddLogical(Params,'Linear System Free Factorization',.FALSE.)

    precProc = AddrFunc(BlockMatrixPrec)

    UseMono = ListGetLogical(Params,'Block MV Monolithic',Found )
    IF(isParallel .AND. .NOT. Found ) THEN
      ! This is a little dangerous logic since it separates serial and parallel operation!      
      UseMono = .NOT. ( ASSOCIATED(TotMatrix % SubMatrix(1,NoVar) % Mat % ParMatrix) .OR. &
          TotMatrix % SubMatrix(1,NoVar) % ParallelIsolatedMatrix )      
    END IF
          
    IF( UseMono ) THEN
      mvProc = AddrFunc(BlockMatrixVectorProdMono)       
    ELSE
      mvProc = AddrFunc(BlockMatrixVectorProd)       
    END IF
      
    prevXnorm = CompNorm(b,ndim)
    WRITE( Message,'(A,ES12.5)') 'Rhs norm at start: ',PrevXnorm
    CALL Info(Caller,Message,Level=10)

    prevXnorm = CompNorm(x,ndim)
    WRITE( Message,'(A,ES12.5)') 'Solution norm at start: ',PrevXnorm
    CALL Info(Caller,Message,Level=10)

    CALL Info(Caller,'Start of blocks system iteration',Level=18)

    ! Always treat the block system as a real valued system and complex
    ! arithmetics only at the inner level.
    CALL ListAddLogical( Params,'Linear System Skip Complex',.TRUE.) 

    IF(ASSOCIATED(SolverMatrix)) THEN
      A => SolverMatrix
    ELSE
      A => TotMatrix % SubMatrix(1,1) % Mat
    END IF

    !IF( ListGetLogical( Solver % Values,'Linear System Complex', Found ) ) A % COMPLEX = .TRUE.

    IF (isParallel) THEN
      ! Note that at this stage we work with the packed vectors x and b
      A => A % ParMatrix % SplittedMatrix % InsideMatrix
      CALL IterSolver( A,x,b,&
          Solver,ndim=ndim,MatvecF=mvProc,PrecF=precProc,&
          DotF=AddrFunc(SParDotProd), NormF=AddrFunc(SParNorm))

      IF (BlockAV) THEN        
        IF(ASSOCIATED(SolverMatrix)) THEN
          ! Communicate the packed solution among all partitions on the shared nodes
          SolverMatrix % ParMatrix % SplittedMatrix % TmpXvec = x(1:poffset(NoVar+1))
          CALL ParallelUpdateResult(SolverMatrix,x,r)
        ELSE
          CALL Fatal(Caller,'No matrix to apply.')
        END IF
      ELSE
        DO i=1,NoVar
          TotMatrix % SubMatrix(i,i) % Mat % ParMatrix % SplittedMatrix % &
              TmpXvec = x(poffset(i)+1:poffset(i+1))
        END DO

        ! Communicate blockwise information.
        ! After this the packed x is again a full vector with redundant shared dofs
        DO i=1,NoVar
          CALL ParallelUpdateResult(TotMatrix % SubMatrix(i,i) % Mat, &
              x(offset(i)+1:offset(i+1)), r )
        END DO
      END IF
    ELSE
      IF( ListGetLogical( Params,'Linear System test',Found ) ) THEN
        CALL IterSolver( A,x,b,&
            Solver,ndim=ndim,MatvecF=mvProc,PrecF=precProc,&
            DotF=AddrFunc(PseudoZDotProd) )
      ELSE
        CALL IterSolver( A,x,b,&
            Solver,ndim=ndim,MatvecF=mvProc,PrecF=precProc) 
      END IF
    END IF
    CALL info(Caller,'Finished block system iteration',Level=18)
    
    CALL ListAddLogical(Params,'Linear System Refactorize',.TRUE.)
    CALL ListAddLogical(Params,'Linear System Free Factorization',.TRUE.)

    !CALL ListPopNamespace()
    
    Xnorm = CompNorm(x,ndim)
    WRITE( Message,'(A,ES12.5)') 'Solution norm: ',Xnorm
    CALL Info(Caller,Message,Level=8)

    MaxChange = 2*ABS(Xnorm-PrevXnorm)/(Xnorm+PrevXnorm)
    PrevXNorm = Xnorm

    WRITE( Message,'(A,ES12.5)') 'Relative change: ',MaxChange
    CALL Info(Caller,Message,Level=8)
    
    DO i=1,NoVar
      TotMatrix % SubVector(i) % Var % Values(1:offset(i+1)-offset(i)) = & 
          x(offset(i)+1:offset(i+1))
    END DO
      
    ! Copy values back since for nontrivial block-matrix structure the
    ! components do not build the whole solution.
    !-----------------------------------------------------------------
    IF( TotMatrix % GotBlockStruct ) THEN
      SVar => CurrentModel % Solver % Variable
      Svar % Values(TotMatrix % BlockPerm) = x
    ELSE IF (BlockAV) THEN
      SolverVar % Values = x
    END IF

    CALL Info(Caller,'Finished block krylov iteration',Level=20)
   
  END SUBROUTINE blockKrylovIter

  
  
  !> This makes the system monolithic. If it was initially monolithic
  !> and then made block, it does not make any sense. However, for
  !> multiphysics coupled cases it may be a good strategy. 
  !-----------------------------------------------------------------
  SUBROUTINE BlockMonolithicSolve( Solver, MaxChange )

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: MaxChange

    INTEGER :: i,j,k,m,n,nc,mc,c,NoVar,NoCol,NoRow,NoEigen,vdofs,comp
    INTEGER, POINTER :: BlockOrder(:)
    LOGICAL :: GotIt, DampedEigen
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:), b(:)
    REAL(KIND=dp), POINTER :: CollX(:), rhs(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, CompVar, SolverVar
    TYPE(Variable_t), TARGET :: MonolithicVar
    REAL(KIND=dp) :: TotNorm
    CHARACTER(:), ALLOCATABLE :: CompName
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: CollMat
    LOGICAL :: Found, HaveMass, HaveDamp, SaveImag, Visited = .FALSE.
    CHARACTER(*), PARAMETER :: Caller = 'BlockMonolithicSolve'
    
    SAVE Visited, CollMat, CollX, HaveMass, HaveDamp, SaveImag
    
    CALL Info(Caller,'Solving block matrix as monolithic!',Level=6)
    
    NoVar = TotMatrix % NoVar
    Params => Solver % Values
    SolverVar => Solver % Variable
    Solver % Variable => MonolithicVar

        
    IF(.NOT. Visited ) THEN
      n = 0
      m = 0

      CollMat => AllocateMatrix()
      HaveMass = .FALSE.
      HaveDamp = .FALSE.
      
      DO NoRow = 1,NoVar 
        rhs => TotMatrix % SubVector(NoRow) % rhs      
        A => TotMatrix % SubMatrix( NoRow, NoRow ) % Mat
        n = n + A % NumberOfRows
      
        DO NoCol = 1,NoVar
          A => TotMatrix % SubMatrix( NoRow, NoCol ) % Mat
          IF(.NOT. ASSOCIATED(A) ) CYCLE
          IF(.NOT. ASSOCIATED(A % Values) ) CYCLE

          IF(InfoActive(20)) THEN
            CALL VectorValuesRange(A % Values,SIZE(A % Values),&
                'A'//I2S(10*NoRow+NoCol),.TRUE.)       
            IF( ASSOCIATED( A % MassValues ) ) THEN
              CALL VectorValuesRange(A % MassValues,SIZE(A % MassValues),&
                  'M'//I2S(10*NoRow+NoCol),.TRUE.)       
            END IF
          END IF
          
          m = m + SIZE( A % Values )
          IF( ASSOCIATED( A % MassValues ) ) HaveMass = .TRUE.
          IF( ASSOCIATED( A % DampValues ) ) HaveDamp = .TRUE.
        END DO
      END DO

      IF( HaveMass ) THEN
        DO NoRow = 1,NoVar 
          A => TotMatrix % SubMatrix( NoRow, NoRow ) % Mat
          IF(.NOT. ASSOCIATED( A % MassValues ) ) THEN
            CALL Warn(Caller,'MassValues are missing for block: '//I2S(11*NoRow))
          END IF
        END DO
        CALL Info(Caller,'Treating MassValues of block matrix too!',Level=20)
      END IF

      IF( HaveDamp ) THEN
        CALL Info(Caller,'Treating DampValues of block matrix too!',Level=20)
      END IF
        
      NoEigen = Solver %  NOFEigenValues

      DampedEigen = ListGetLogical(Solver % Values,'Eigen System Complex',Found )  
      IF( DampedEigen ) THEN
        CALL Info(Caller,'Creating complex system for eigen values!',Level=6)
      ELSE
        CALL Info(Caller,'Creating real valued system for eigen values!',Level=8)        
      END IF
      
      SaveImag = ListGetLogical(Solver % Values,'Pick Im Component',Found )  
    
      ! The matrix sizes depend on whether we create a complex or real valued system. 
      IF(DampedEigen) THEN
        nc = 2*n
        mc = 4*m
      ELSE
        nc = n
        mc = m
      END IF

      
      CollMat % NumberOfRows = nc
      CALL Info(Caller,'Size of monolithic matrix: '//I2S(nc),Level=7)
      CALL Info(Caller,'Estimated number of nonzeros in monolithic matrix: '//I2S(mc),Level=7)
      
      ALLOCATE( CollMat % rhs(nc), CollMat % Diag(nc), CollMat % Rows(nc+1), CollX(nc) )
      CollMat % rhs = 0.0_dp
      CollMat % Diag = 0
      CollMat % Rows = 0
      CollX = 0.0_dp

      CollMat % Complex = DampedEigen
      
      ALLOCATE( CollMat % Values(mc), CollMat % Cols(mc+1) )
      CollMat % Values = 0.0_dp
      CollMat % Cols = 0

      IF( HaveMass ) THEN
        ALLOCATE( CollMat % MassValues(mc) )
        CollMat % MassValues = 0.0_dp
      END IF
      IF( HaveDamp .AND. .NOT. DampedEigen ) THEN
        ALLOCATE( CollMat % DampValues(mc) )
        CollMat % DampValues = 0.0_dp
      END IF
    END IF

    k = 0
    CollMat % Rows(1) = 1
  

    IF(DampedEigen) THEN
      DO NoRow = 1,NoVar
        n = TotMatrix % Offset(NoRow)
        
        A => TotMatrix % SubMatrix( NoRow, NoRow ) % Mat
        m = A % NumberOfRows

        DO i=1,m

          ! Loop over real and imaginary rows
          DO c=1,2
          
            DO NoCol = 1,NoVar
              A => TotMatrix % SubMatrix( NoRow, NoCol ) % Mat
              IF( .NOT. ASSOCIATED( A ) ) CYCLE
              IF( .NOT. ASSOCIATED( A % Rows ) ) CYCLE
              
              DO j=A % Rows(i),A % Rows(i+1)-1
                ! If we have the imaginary row add the multiplier of imaginary value first
                IF( c == 2 ) THEN
                  k = k + 1
                  CollMat % Cols(k) = 2*TotMatrix % Offset(NoCol) + 2*A % Cols(j)-1
                  IF( ASSOCIATED( A % DampValues) ) THEN
                    CollMat % Values(k) = A % DampValues(j)
                  END IF
                END IF

                k = k + 1
                CollMat % Values(k) = A % Values(j)
                CollMat % Cols(k) = 2*TotMatrix % Offset(NoCol) + 2*(A % Cols(j)-1)+c
                IF( ASSOCIATED( A % MassValues ) ) THEN
                  CollMat % MassValues(k) = A % MassValues(j)
                END IF

                ! If we have the real row add the multiplier of imaginary value last
                IF( c == 1 ) THEN
                  k = k + 1
                  CollMat % Cols(k) = 2*TotMatrix % Offset(NoCol) + 2*A % Cols(j)
                  IF( ASSOCIATED( A % DampValues) ) THEN
                    CollMat % Values(k) = -A % DampValues(j)
                  END IF
                END IF
              END DO
            END DO
            IF(.NOT. Visited ) THEN
              CollMat % Rows(2*((n+i)-1)+c+1) = k+1
            END IF
          END DO
        END DO
      END DO

    ELSE
      DO NoRow = 1,NoVar
        n = TotMatrix % Offset(NoRow)

        A => TotMatrix % SubMatrix( NoRow, NoRow ) % Mat
        m = A % NumberOfRows

        DO i=1,m
          DO NoCol = 1,NoVar
            A => TotMatrix % SubMatrix( NoRow, NoCol ) % Mat
            IF( .NOT. ASSOCIATED( A ) ) CYCLE
            IF( .NOT. ASSOCIATED( A % Rows ) ) CYCLE
            IF( SIZE(A % Rows) < i+1 ) CYCLE
            
            DO j=A % Rows(i),A % Rows(i+1)-1
              k = k + 1
              CollMat % Values(k) = A % Values(j)
              IF(.NOT. Visited ) THEN
                CollMat % Cols(k) = TotMatrix % Offset(NoCol) + A % Cols(j)
              END IF
              IF( HaveMass ) THEN
                IF( ASSOCIATED( A % MassValues ) ) THEN
                  CollMat % MassValues(k) = A % MassValues(j)
                END IF
              END IF
              IF( HaveDamp ) THEN
                IF( ASSOCIATED( A % DampValues ) ) THEN
                  CollMat % DampValues(k) = A % DampValues(j)
                END IF
              END IF
            END DO
            IF( ASSOCIATED( A % rhs ) ) THEN
              CollMat % rhs(n+i) = A % rhs(i)
            END IF
          END DO
          IF(.NOT. Visited ) THEN
            CollMat % Rows(n+i+1) = k+1
          END IF
        END DO
      END DO
    END IF


    IF( ParEnv % PEs > 1 ) THEN
      IF( NoVar /= 2 ) THEN
        CALL Fatal(Caller,'Parallel operation currently assumes just 2x2 blocks!')
      END IF

      CALL Info(Caller,'Merging parallel info of matrices')
      CALL ParallelMergeMatrix( Solver, CollMat, &
          TotMatrix % SubMatrix(1,1) % Mat, &
          TotMatrix % SubMatrix(2,2) % Mat )
    END IF

    
    IF(InfoActive(20)) THEN
      !CALL CRS_CheckSymmetricTopo(CollMat)
      !CALL CRS_CheckComplexTopo(CollMat)
      CALL VectorValuesRange(CollMat % Values,SIZE(CollMat % Values),'Atot',.TRUE.)
      IF( ASSOCIATED( CollMat % MassValues ) ) THEN
        CALL VectorValuesRange(CollMat % MassValues,SIZE(CollMat % MassValues),'Mtot',.TRUE.)
      END IF
    END IF
          
    CALL Info(Caller,'True number of nonzeros in monolithic matrix: '//I2S(k),Level=7)
 
    IF(.NOT. Visited ) THEN
      MonolithicVar % Name = '' ! Some name needed to avoid an uninitialised value error
      MonolithicVar % Values => CollX
      MonolithicVar % Dofs = 1
      MonolithicVar % Perm => NULL()
      
      NoEigen = Solver %  NOFEigenValues
      IF( NoEigen > 0 ) THEN
        n = CollMat % NumberOfRows
        ALLOCATE( MonolithicVar % EigenValues(NoEigen), MonolithicVar % EigenVectors(NoEigen,n) )
        MonolithicVar % EigenValues = 0.0_dp
        MonolithicVar % EigenVectors = 0.0_dp
      END IF
      Visited = .TRUE.      
    END IF

    IF(.NOT. DampedEigen) THEN        
      CALL Info(Caller,'Copying block solution to monolithic vector',Level=12)
      DO i=1,NoVar
        Var => TotMatrix % SubVector(i) % Var 
        n = SIZE( Var % Values )       
        m = TotMatrix % Offset(i)
        CollX(m+1:m+n) = Var % Values(1:n) 
      END DO
    END IF
      
    ! Solve monolithic matrix equation. 
    CALL SolveLinearSystem( CollMat, CollMat % rhs, CollX, TotNorm, 1, Solver )

    CALL Info(Caller,'Copying monolithic vector to block solutions',Level=12)

    ! Copy the 1st eigenmode because this will be used for norms etc. 
    NoEigen = Solver % NOFEigenValues
    IF( NoEigen > 0 ) THEN
      MonolithicVar % Values = REAL( MonolithicVar % EigenVectors(1,:) ) 
    END IF
    
    DO i=1,NoVar
      Var => TotMatrix % SubVector(i) % Var 
      n = SIZE( Var % Values ) 
      m = TotMatrix % Offset(i)
      Var % Values(1:n) = CollX(m+1:m+n) 
      
      IF( NoEigen > 0 ) THEN
        IF(.NOT. ASSOCIATED( Var % EigenValues ) ) THEN
          IF( ASSOCIATED( Var % Solver ) ) THEN
            Var % Solver % NOFEigenValues = NoEigen
          END IF
          ALLOCATE( Var % EigenValues(NoEigen), Var % EigenVectors(NoEigen,n) )
          Var % EigenValues = 0.0_dp
          Var % EigenVectors = 0.0_dp

          vdofs = Var % Dofs
          IF( vdofs > 1 ) THEN
            DO comp=1,vdofs
              CompName = ComponentName(Var % Name,comp)
              CompVar => VariableGet( Solver % Mesh % Variables, Compname )
              CompVar % EigenValues => Var % EigenValues
              CompVar % EigenVectors => Var % EigenVectors(:,comp::vdofs)
            END DO
          END IF
        END IF

        DO k=1,NoEigen                    
          Var % EigenValues(k) = MonolithicVar % EigenValues(k)
          Var % EigenVectors(k,1:n) = MonolithicVar % EigenVectors(k,m+1:m+n)
        END DO
        
      END IF
    END DO

    Solver % Variable => SolverVar
    
  END SUBROUTINE BlockMonolithicSolve


!------------------------------------------------------------------------------
!> An alternative handle for the block solvers to be used by the legacy matrix
!> type. 
!------------------------------------------------------------------------------
  SUBROUTINE BlockSolveInt(A,x,b,Solver)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp), TARGET :: x(:)
    REAL(KIND=dp), TARGET CONTIG :: b(:)
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i,j,k,l,n,nd,NonLinIter,tests,NoTests,iter
    LOGICAL :: GotIt, GotIt2, BlockPrec, BlockAV, &
        BlockHdiv, BlockHcurl, BlockHorVer, BlockCart, BlockNodal, BlockDomain, &
        BlockDummy, BlockComplex 
    INTEGER :: ColVar, RowVar, NoVar, BlockDofs, VarDofs
    
    REAL(KIND=dp) :: NonlinearTol, Norm, PrevNorm, Residual, PrevResidual, &
        TotNorm, MaxChange, alpha, beta, omega, rho, oldrho, s, r, PrevTotNorm, &
        Coeff
    REAL(KIND=dp), POINTER :: SaveValues(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
    LOGICAL :: Robust, LinearSearch, ErrorReduced, IsProcedure, ScaleSystem,&
        LS, BlockScaling, BlockMonolithic, Found
    INTEGER :: HaveConstraint, HaveAdd
    INTEGER, POINTER :: VarPerm(:)
    INTEGER, POINTER :: BlockPerm(:)
    INTEGER, POINTER :: SlaveSolvers(:)
    LOGICAL :: GotSlaveSolvers, SkipVar
    
    
    TYPE(Matrix_t), POINTER :: Amat, SaveCM
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params

    INTEGER, POINTER :: BlockIndex(:)
    

    CALL Info('BlockSolveInt','---------------------------------------',Level=5)

    Params => Solver % Values
    Mesh => Solver % Mesh
    PSolver => Solver

    SolverRef => Solver

    isParallel = ParEnv % PEs > 1
    
    
    ! Determine some parameters related to the block strategy
    !------------------------------------------------------------------------------
    BlockPrec = ListGetLogical( Params,'Block Preconditioner',GotIt)
    IF(.NOT. GotIt) THEN
      CALL Info('BlockSolveInt','Using block preconditioning mode by default')
      BlockPrec = .TRUE.
    END IF

    BlockMonolithic = ListGetLogical( Params,'Block Monolithic',GotIt)
    
    BlockScaling = ListGetLogical( Params,'Block Scaling',GotIt)

    ! Different strategies on how to split the initial monolithic matrix into blocks
    BlockAV = ListGetLogical( Params,'Block A-V System', GotIt)
    BlockHcurl = ListGetLogical( Params,'Block Quadratic Hcurl System', GotIt)   
    BlockHdiv = ListGetLogical( Params,'Block Hdiv system',GotIt)
    BlockNodal = ListGetLogical( Params,'Block Nodal System', GotIt)
    BlockHorVer = ListGetLogical( Params,'Block Hor-Ver System', GotIt)
    BlockCart = ListGetLogical( Params,'Block Cartesian System', GotIt)
    BlockDomain = ListGetLogical( Params,'Block Domain System',GotIt) 
    BlockDummy = ListGetLogical( Params,'Block Nested System',GotIt)
    BlockComplex = ListGetLogical( Params,'Block Complex System',GotIt) 
    
    SlaveSolvers =>  ListGetIntegerArray( Params, &
         'Block Solvers', GotSlaveSolvers )

    SkipVar = .FALSE.
    IF( BlockDomain ) THEN
      n = MAXVAL( Solver % Variable % Perm )
      ALLOCATE( BlockIndex(n) )
      BlockDofs = 0
      CALL BlockPickDofsPhysical( PSolver, BlockIndex, BlockDofs )  
      SkipVar = .TRUE.
    ELSE IF( BlockHdiv ) THEN
      n = MAXVAL( Solver % Variable % Perm )
      ALLOCATE( BlockIndex(n) )
      BlockDofs = 0
      CALL BlockPickHdiv( PSolver, BlockIndex, BlockDofs )  
      SkipVar = .TRUE.
    ELSE IF( BlockAV .OR. BlockNodal .OR. BlockHorVer .OR. BlockHcurl ) THEN
      BlockDofs = 2
      IF( ListGetLogical( Params,'Block Quadratic Hcurl Faces',Found ) ) BlockDofs = 3
      IF(BlockComplex) THEN
        BlockDofs = 2 * BlockDofs 
        IF( ListGetLogical( Params,'Block Quadratic Hcurl semicomplex',Found ) ) BlockDofs = 3
      END IF
      SkipVar = .TRUE.
    ELSE IF( BlockCart ) THEN
      BlockDofs = 3
      SkipVar = .TRUE.
    ELSE IF( GotSlaveSolvers ) THEN
      BlockDofs = SIZE( SlaveSolvers )
    ELSE IF( BlockDummy ) THEN
      BlockDofs = 1
    ELSE
      BlockDofs = Solver % Variable % Dofs      
    END IF
    VarDofs = BlockDofs

    
    HaveConstraint = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  HaveConstraint = 1
    HaveConstraint = ParallelReduction(HaveConstraint)
     
    HaveAdd = 0
    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows > 0 ) HaveAdd = 1
    END IF
    HaveAdd = ParallelReduction(HaveAdd)

    IF( HaveConstraint > 0 ) BlockDofs = BlockDofs + 1
    IF( HaveAdd > 0 ) BlockDofs = BlockDofs + 1    
   
    CALL BlockInitMatrix( Solver, TotMatrix, BlockDofs, VarDofs, SkipVar )
    
    NoVar = TotMatrix % NoVar
    TotMatrix % Solver => Solver
    
    SaveMatrix => Solver % Matrix
    SolverMatrix => A
    Solver % Matrix => A

    SolverVar => Solver % Variable
    SaveValues => SolverVar % Values
    SolverVar % Values => x

    SaveRHS => SolverMatrix % RHS
    SolverMatrix % RHS => b
    
    IF( .NOT. GotSlaveSolvers ) THEN
      CALL Info('BlockSolveInt','Splitting monolithic matrix into pieces',Level=10)
      IF( BlockDomain .OR. BlockHdiv ) THEN
        CALL BlockPickMatrixPerm( Solver, BlockIndex, VarDofs )
        DEALLOCATE( BlockIndex ) 
      ELSE IF( BlockAV ) THEN
        CALL BlockPickMatrixAV( Solver, VarDofs )
      ELSE IF( BlockHcurl ) THEN
        CALL BlockPickMatrixHcurl( Solver, VarDofs, BlockComplex )
        IF(VarDofs == 3 ) BlockComplex = .FALSE.
      ELSE IF( BlockHorVer .OR. BlockCart ) THEN
        CALL BlockPickMatrixHorVer( Solver, VarDofs, BlockCart )       
      ELSE IF( BlockNodal ) THEN
        CALL BlockPickMatrixNodal( Solver, VarDofs )        
      ELSE IF( BlockDummy .OR. VarDofs == 1 ) THEN
        CALL Info('BlockSolveInt','Using the original matrix as the (1,1) block!',Level=10)
        TotMatrix % SubMatrix(1,1) % Mat => SolverMatrix        
        TotMatrix % SubMatrix(1,1) % Mat % Complex = ListGetLogical(Params,'Linear System Complex',Found)        
      ELSE
        CALL BlockPickMatrix( Solver, NoVar ) !VarDofs )
        VarDofs = NoVar
      END IF

      IF( SkipVar ) THEN
        CALL BlockInitVar( Solver, TotMatrix )
      END IF

      CALL BlockPrecMatrix( Solver, VarDofs ) 
      CALL Info('BlockSolveInt','Block matrix system created',Level=12)
    END IF
    
    ! Currently we cannot have both structure-structure and fluid-structure couplings!
    IF( ListGetLogical( Solver % Values,'Structure-Structure Coupling',Found ) ) THEN
      CALL StructureCouplingBlocks( Solver )
    ELSE
      CALL FsiCouplingBlocks( Solver )
    END IF
        
    IF( HaveConstraint > 0 ) THEN
      CALL BlockPickConstraint( Solver, VarDofs )
      ! Storing pointer to CM so that the SolverMatrix won't be treated with the constraints
      SaveCM => Solver % Matrix % ConstraintMatrix
      Solver % Matrix % ConstraintMatrix => NULL()
    END IF

    IF (isParallel) THEN
      CALL Info('BlockSolveInt','Initializing parallel block matrices',Level=12)
      DO RowVar=1,NoVar
        DO ColVar=1,NoVar

          Amat => TotMatrix % SubMatrix(RowVar,ColVar) % Mat

          ! It is reasonable that there is no parallel communication for parallel
          ! matrix if some partition is not participating in the process. 
          IF(Amat % NumberOfRows == 0) CYCLE

          Amat % Comm = Solver % Matrix % Comm
          Parenv % ActiveComm = Amat % Comm
          Solver % Variable => TotMatrix % SubVector(ColVar) % Var
          
          ! This is a coupling matrix that should by construction not lie at the interface
          IF(TotMatrix % Submatrix(RowVar,ColVar) % ParallelIsolatedMatrix ) CYCLE

          ! The parallel solution and hence also initialization works only for
          ! standard square matrices. 
          GotIt = (Amat % NumberOfRows == MAXVAL(Amat % Cols)) 
          TotMatrix % Submatrix(RowVar,ColVar) % ParallelSquareMatrix = GotIt            

          IF(GotIt) THEN
            IF (.NOT.ASSOCIATED(Amat % ParMatrix)) THEN
              CALL ParallelInitMatrix(Solver,Amat)
            END IF
            IF(ASSOCIATED(Amat % ParMatrix )) THEN
              Amat % ParMatrix % ParEnv % ActiveComm = Amat % Comm
              ParEnv => Amat % ParMatrix % ParEnv
            END IF
            !CALL SParIterActiveBarrier()
          END IF
        END DO
      END DO

      CALL ParallelShrinkPerm()
        
      Solver % Variable  => SolverVar
      CALL Info('BlockSolveInt','Initialization of block matrix finished',Level=20)
    END IF
    
    !------------------------------------------------------------------------------
    ! Finally solve the system using 'outer: ' as the optional namespace
    ! for the linear system setting.
    !------------------------------------------------------------------------------          
      
    TotNorm = 0.0_dp
    MaxChange = 0.0_dp

    IF( BlockScaling ) THEN   
      CALL CreateBlockMatrixScaling()
      CALL BlockMatrixScaling(.FALSE.)
    END IF

    CALL ListPushNamespace('outer:')

    IF (BlockScaling) THEN
      ! This simplifies writing a consistent sif file:
      CALL ListAddLogical(Solver % Values, 'Linear System Row Equilibration', .TRUE.)      
    END IF
    
    ! The case with one block is mainly for testing and developing features
    ! related to nonlinearity and assembly.
    !----------------------------------------------------------------------
    IF( NoVar == 1 .AND. .NOT. BlockDummy ) THEN
      CALL Info('BlockSolveInt','Solving in standard manner',Level=6)
      
      Solver % Variable => TotMatrix % SubVector(1) % Var
      Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat
      
      TotNorm = DefaultSolve()
      MaxChange = Solver % Variable % NonlinChange 
    ELSE IF( BlockMonolithic ) THEN
      CALL Info('BlockSolveInt','Using monolithic strategy for the block',Level=6)        
      CALL BlockMonolithicSolve( Solver, MaxChange )      
    ELSE IF( BlockPrec ) THEN
      CALL Info('BlockSolveInt','Using block preconditioning strategy',Level=6)        
      CALL BlockKrylovIter( Solver, MaxChange )
    ELSE
      Solver % Variable => TotMatrix % SubVector(1) % Var
      Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat
      
      CALL Info('BlockSolveInt','Using block solution strategy',Level=6)
      CALL BlockStandardIter( Solver, MaxChange )
    END IF

    CALL ListPopNamespace('outer:')

    IF( BlockScaling ) THEN
      CALL BlockMatrixScaling(.TRUE.)
      CALL DestroyBlockMatrixScaling()
    END IF

    ! For legacy matrices do the backmapping 
    !------------------------------------------
    SolverMatrix % RHS => SaveRHS
    Solver % Matrix => SaveMatrix
    Solver % Variable => SolverVar
    Solver % Variable % Values => SaveValues
       
    IF( HaveConstraint > 0 ) THEN
      ! Restore the pointer to the SolverMatrix
      Solver % Matrix % ConstraintMatrix => SaveCM 
    END IF

    IF( BlockHorVer .OR. BlockCart .OR. BlockDomain .OR. BlockHcurl ) THEN
      CALL BlockBackCopyVar( Solver, TotMatrix )
    END IF

    IF( ListGetLogical( Solver % Values,'Block Save Iterations',Found ) ) THEN
      CALL ListAddInteger(CurrentModel % Simulation,'res: block iterations',TotMatrix % NoIters)
    END IF   
    
    CALL Info('BlockSolveInt','All done')
    CALL Info('BlockSolveInt','-------------------------------------------------',Level=5)

  END SUBROUTINE BlockSolveInt
END MODULE BlockSolve


!------------------------------------------------------------------------------
!> Just a handle for SolveLinearSystem():
!------------------------------------------------------------------------------
SUBROUTINE BlockSolveExt(A,x,b,Solver)
!------------------------------------------------------------------------------
    USE BlockSolve

    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: x(:),b(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, bm
!------------------------------------------------------------------------------

    ! Eliminate recursion for block solvers. 
    bm = ListGetLogical(  Solver % Values, 'Linear System Block Mode', Found)
    IF(Found) &
      CALL ListAddLogical(Solver % Values,'Linear System Block Mode',.FALSE.)

    CALL BlockSolveInt(A,x,b,Solver)

    IF(Found) &
      CALL ListAddLogical(Solver % Values,'Linear System Block Mode',bm )
!------------------------------------------------------------------------------
END SUBROUTINE BlockSolveExt
!------------------------------------------------------------------------------

!> \}
