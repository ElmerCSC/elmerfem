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

  TYPE(Variable_t), POINTER :: SolverVar => Null()
  TYPE(Matrix_t), POINTER :: SolverMatrix => Null()

CONTAINS


  
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
    CHARACTER(LEN=max_name_len) :: VarName
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
    CHARACTER(LEN=max_name_len) :: str, eq
    
    LOGICAL :: GlobalBubbles, Found
    INTEGER :: MaxNDOFs, MaxDGDOFs, MaxEDOFs, MaxFDOFs, MaxBDOFs


    CALL Info('BlockSolver','Creating block variables',Level=8)
    
    Mesh => Solver % Mesh
    Params => Solver % Values

    IF( PRESENT( ExtDofs ) ) THEN
      Dofs = ExtDofs
    ELSE 
      WRITE (str,'(A,I0,A)') 'Variable ',VariableNo,' Dofs'
      Dofs = ListGetInteger( Params, TRIM(str), GotIt )
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
        
        str='Active Variables['//TRIM(i2s(solver_id))//']'            
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
    CALL Info('BlockSolver',Message)
    
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
    TYPE(Matrix_t), POINTER :: Amat, Bmat
    INTEGER :: i,j,n,Novar
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=max_name_len) :: VarName, str
    LOGICAL :: UseSolverMatrix, IsComplex
    
            
    Params => Solver % Values

    IsComplex = ListGetLogical( Params,'Linear System Complex',Found)

    BlockMatrix => Solver % BlockMatrix
    IF (ASSOCIATED(BlockMatrix)) THEN
      CALL Info('BlockInitMatrix','Using existing block matrix',Level=10)
      RETURN
    END IF

    CALL Info('BlockInitMatrix','Initializing block matrix',Level=10)
    
    ALLOCATE(Solver % BlockMatrix)
    BlockMatrix => Solver % BlockMatrix
 
    BlockStruct => ListGetIntegerArray( Params,'Block Structure',GotBlockStruct)
    IF( GotBlockStruct ) THEN
      IF( SIZE( BlockStruct ) /= BlockDofs ) THEN
        CALL Fatal('BlockInitMatrix','Incompatible size of > Block Structure < given!')
      END IF
      IF( MINVAL( BlockStruct ) < 1 .OR. MAXVAL( BlockStruct ) > BlockDofs ) THEN
        CALL Fatal('BlockInitMatrix','Incompatible values in > Block Structure < given!')          
      END IF
      NoVar = MAXVAL( BlockStruct )
      BlockMatrix % BlockStruct => BlockStruct
    ELSE
      NoVar = BlockDofs
    END IF    
    BlockMatrix % GotBlockStruct = GotBlockStruct


    IF( BlockMatrix % NoVar == NoVar ) THEN
      CALL Info('BlockInitMatrix','Reusing existing blockmatrix',Level=6)
      RETURN
    ELSE IF( BlockMatrix % Novar /= 0 ) THEN
      CALL Fatal('BlockInitMatrix','Previous blockmatrix was of different size?')
    ELSE
      CALL Info('BlockInitMatrix','Starting',Level=6)
    END IF
    
    BlockMatrix % Solver => Solver
    BlockMatrix % NoVar = NoVar

    ALLOCATE( BlockMatrix % SubMatrix(NoVar,NoVar) )
    DO i=1,NoVar
      DO j=1,NoVar
        Amat => AllocateMatrix()
        Amat % ListMatrix => NULL()
        Amat % FORMAT = MATRIX_LIST      
        Amat % NumberOfRows = 0
        AMat % Complex = IsComplex
        BlockMatrix % Submatrix(i,j) % Mat => Amat

        Bmat => AllocateMatrix()
        Bmat % ListMatrix => NULL()
        Bmat % FORMAT = MATRIX_LIST      
        Bmat % NumberOfRows = 0
        BMat % Complex = IsComplex
        BlockMatrix % Submatrix(i,j) % PrecMat => Bmat
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
      IF( SkipVar ) RETURN
    END IF
    
    IF( PRESENT( FieldDofs ) ) THEN
      NoVar = FieldDofs
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

        CALL Info('BlockInitMatrix','Associating block '//TRIM(I2S(i))//' with solver: '//TRIM(I2S(j)),Level=10)

        PSolver => CurrentModel % Solvers(j)
        Var => PSolver % Variable 
        VarName = Var % Name

        BlockMatrix % SubVector(i) % Solver => PSolver
        BlockMatrix % SubMatrix(i,i) % Mat => PSolver % Matrix        

      ELSE
        WRITE (str,'(A,I0)') 'Variable ',i

        VarName = ListGetString( Params, TRIM(str), Found )
        IF(.NOT. Found ) THEN       
          IF( BlockMatrix % GotBlockStruct ) THEN
            WRITE (VarName,'(A,I0)') 'BlockVar ',i
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
        CALL Info('BlockInitMatrix','Using existing variable > '//TRIM(VarName)//' <')		
      ELSE		
        CALL Info('BlockInitMatrix','Variable > '//TRIM(VarName)//' < does not exist, creating')
        PSolver => Solver
        IF( BlockMatrix % GotBlockStruct ) THEN
          j = COUNT( BlockMatrix % BlockStruct == i ) 
          IF( j == 0 ) THEN
            CALL Fatal('BlockInitMatrix','Invalid > Block Structure < given!')
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


    CALL Info('BlockInitMatrix','All done',Level=12)
      
  END SUBROUTINE BlockInitMatrix
    


  !-------------------------------------------------------------------
  !> This subroutine creates the minssing component variables.
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
    CHARACTER(LEN=max_name_len) :: VarName, str
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
        CALL Info('BlockInitMatrix','Variable > '//TRIM(VarName)//' < does not exist, creating')
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
    LOGICAL :: ReuseMatrix, Found
    INTEGER::i,j,k,i_aa,i_vv,i_av,i_va,n;
    REAL(KIND=DP) :: SumAbsMat
    
    CALL Info('BlockSolver','Picking block matrix from monolithic one',Level=10)

    SolverMatrix => Solver % Matrix 
    Params => Solver % Values
        
    ReuseMatrix = ListGetLogical( Params,'Block Matrix Reuse',Found)

    DO RowVar=1,NoVar
      DO ColVar=1,NoVar            
        Amat => TotMatrix % Submatrix(RowVar,ColVar) % Mat          
        IF( TotMatrix % GotBlockStruct) THEN
          ! A generic picking method for submatrices
          !----------------------------------------------------------------------
          CALL Info('BlockSolver','Picking generic block matrix ('&
              //TRIM(I2S(RowVar))//','//TRIM(I2S(ColVar))//')',Level=20)
          CALL CRS_BlockMatrixPick2(SolverMatrix,Amat,TotMatrix % BlockStruct,RowVar,ColVar)
        ELSE
          ! Picking of standard submatrices of size one.
          !----------------------------------------------------------------------
          IF( ReuseMatrix ) THEN
            IF( RowVar + ColVar > 2 .AND. Amat % NumberOfRows == 0 ) THEN
              CALL Info('BlockSolver','Copying block matrix topology ('&
                  //TRIM(I2S(RowVar))//','//TRIM(I2S(ColVar))//')',Level=20)
              CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(1,1) % Mat, Amat )
            END IF
          END IF
          CALL Info('BlockSolver','Picking simple block matrix ('&
              //TRIM(I2S(RowVar))//','//TRIM(I2S(ColVar))//')',Level=20)          
          CALL CRS_BlockMatrixPick(SolverMatrix,Amat,NoVar,RowVar,ColVar)          
            
          IF( Amat % NumberOfRows > 0 ) THEN
            SumAbsMat = SUM( ABS( Amat % Values ) )
            IF( SumAbsMat < SQRT( TINY( SumAbsMat ) ) ) THEN
              CALL Info('BlockSolver','Matrix is actually all zero, eliminating it!',Level=20)
              DEALLOCATE( Amat % Values ) 
              IF( .NOT. ReuseMatrix ) THEN
                DEALLOCATE( Amat % Rows, Amat % Cols )
                IF( RowVar == ColVar ) DEALLOCATE( Amat % Diag )
              END IF
              Amat % NumberOfRows = 0
            END IF
          END IF

        END IF
      END DO
    END DO

  END SUBROUTINE BlockPickMatrix


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
    
    CALL Info('BlockSolverAV','Picking block matrix from monolithic one',Level=10)

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
    ! pointer should not be assicoated. 
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
    
    
    CALL Info('BlockPickMatrixHorVer','Dividing matrix in vertical and horizontal dofs',Level=10)


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
            IF( DTag(j) == 2 ) PRINT *,'Vertical edge '//TRIM(I2S(j))//' is also horizontal?'
            DTag(j) = 1  ! set to be vertical
          ELSE IF( Wproj < Wtol ) THEN
            IF( DTag(j) == 1 ) PRINT *,'Horizontal edge '//TRIM(I2S(j))//' is also vertical?'
            DTag(j) = 2  ! set to be horizontal
          ELSE
            PRINT *,'Edge '//TRIM(I2S(j))//' direction undefined: ',Wproj
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

    PRINT *,'Cartesian dofs:',ndir(1:NoVar)

    i = n - SUM( ndir ) 
    IF( i > 0 ) THEN      
      CALL Fatal('BlockPickMatrixHorVer','Could not determine all nodes: '&
          //TRIM(I2S(i)))
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
          STOP
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
    LOGICAL :: ReuseMatrix, Found

    LOGICAL :: BlockAV
    INTEGER::i,j,k,l,n,i1,i2,i3,rowi,colj,NoCon,rb,cb
    TYPE(Matrix_t), POINTER :: A,CM,C1,C2,C3,C1prec,C2prec,C3prec
    REAL(KIND=dp) :: PrecCoeff,val
    INTEGER, POINTER :: ConsPerm(:)
    INTEGER :: DoPrec
    CHARACTER(LEN=max_name_len) :: VarName
    TYPE(Variable_t), POINTER :: Var
    TYPE(Solver_t), POINTER :: PSolver
    LOGICAL :: InheritCM, PrecTrue 
    
    CALL Info('BlockSolver','Picking constraints to block matrix',Level=10)

    
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

    CALL Info('BlockPickConstraint','Number of constraint matrices: '//TRIM(I2S(NoCon)),Level=10)
    
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
      
    CALL Info('BlockSolver','Setting format of constraint blocks to CRS',Level=20)
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
      CALL Info('BlockSolver','Using existing variable > '//TRIM(VarName)//' <')		
    ELSE		
      CALL Info('BlockSolver','Variable > '//TRIM(VarName)//' < does not exist, creating')
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

    INTEGER :: i, RowVar, ColVar
    CHARACTER(LEN=max_name_len) :: str
    REAL(KIND=dp) :: Coeff
    LOGICAL :: GotIt, GotIt2
    INTEGER, POINTER :: VarPerm(:)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Amat

    Params => Solver % Values
 
    ! The user may give a user defined preconditioner matrix
    !-----------------------------------------------------------
    DO RowVar=1,NoVar
      i = TotMatrix % Submatrix(RowVar,RowVar) % PrecMat % NumberOfRows 
      IF( i > 0 ) CYCLE
      
      WRITE (str,'(A,I0)') 'Prec Matrix Diffusion ',RowVar
      Coeff = ListGetCReal( Params, TRIM(str), GotIt)
      
      WRITE (str,'(A,I0)') 'Prec Matrix Density ',RowVar
      Coeff = ListGetCReal( Params, TRIM(str), GotIt2)
      
      IF( GotIt .OR. GotIt2 ) THEN
        CALL Info('BlockPrecMatrix','Creating simple preconditioning matrix')
        
        CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(RowVar,RowVar) % Mat, &
            TotMatrix % Submatrix(RowVar,RowVar) % PrecMat )   
        
        Amat => TotMatrix % Submatrix(RowVar,RowVar) % PrecMat
        VarPerm => TotMatrix % Subvector(RowVar) % Var % Perm
        IF( GotIt ) THEN
          CALL LaplaceMatrixAssembly( Solver, VarPerm, Amat )
          Amat % Values = Coeff * Amat % Values
        ELSE 
          CALL MassMatrixAssembly( Solver, VarPerm, Amat )
          Amat % Values = Coeff * Amat % Values
        END IF
      END IF
    END DO
  END SUBROUTINE BlockPrecMatrix


  !> Create the coupling blocks for a linear FSI coupling among various types of
  !> elasticity and fluid solvers.
  !--------------------------------------------------------------------------------
  SUBROUTINE FsiCouplingBlocks( Solver )

    TYPE(Solver_t) :: Solver
    INTEGER :: i,j,Novar

    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A_fs, A_sf, A_s, A_f
    TYPE(Variable_t), POINTER :: FVar, SVar
    LOGICAL :: IsPlate, IsShell, IsNs, IsPres
    
    Params => Solver % Values

    IsPlate = .FALSE.
    IsShell = .FALSE.
    IsNS = .FALSE.
    IsPres = .FALSE.
    
    i = ListGetInteger( Params,'Structure Solver Index',Found)
    IF( Found ) THEN
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
      
    j = ListGetInteger( Params,'Fluid Solver Index',Found)
    IF(.NOT. Found ) THEN
      j = ListGetInteger( Params,'NS Solver Index', IsNs )
      IF( .NOT. IsNs ) THEN
        j = ListGetInteger( Params,'Pressure Solver Index',IsPres)
      END IF
    END IF
    IF( j == 0 ) THEN
      IF( i > 1 .AND. TotMatrix % NoVar == 2 ) j = 3 - i 
    END IF
    IF( i == 0 ) THEN
      IF( j > 1 .AND. TotMatrix % NoVar == 2 ) i = 3 - j
    END IF      
    
    IF(i<=0 .OR. j<=0) THEN
      IF( i > 0 ) CALL Warn('FsiCouplingBlocks','Structure solver given but not fluid!')
      IF( j > 0 ) CALL Warn('FsiCouplingBlocks','Fluid solver given but not structure!')
      RETURN
    END IF
      
    A_fs => TotMatrix % Submatrix(j,i) % Mat
    A_sf => TotMatrix % Submatrix(i,j) % Mat
    
    IF(.NOT. ASSOCIATED( A_fs ) ) THEN
      CALL Fatal('FsiCouplingBlocks','Fluid-structure coupling matrix not allocated!')
    END IF
    IF(.NOT. ASSOCIATED( A_sf ) ) THEN
      CALL Fatal('FsiCouplingBlocks','Structure-fluid coupling matrix not allocated!')
    END IF
       
    SVar => TotMatrix % Subvector(i) % Var
    FVar => TotMatrix % Subvector(j) % Var

    A_s => TotMatrix % Submatrix(i,i) % Mat
    A_f => TotMatrix % Submatrix(j,j) % Mat
    
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL Fatal('FsiCouplingBlocks','Fluid variable not present!')
    END IF
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL Fatal('FsiCouplingBlocks','Structure variable not present!')
    END IF

    IF(.NOT. (IsNs .OR. IsPres ) ) THEN
      IsPres = ( FVar % Dofs <= 2 )
      IsNs = .NOT. IsPres
    END IF
    
    CALL FsiCouplingAssembly( Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
        IsPlate, IsShell, IsNS )
      
  END SUBROUTINE FsiCouplingBlocks
    

  
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
    LOGICAL :: GotRhs, DoSum
    
    CALL Info('BlockUpdateRhs','Computing block r.h.s',Level=5)

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
        CALL Info('BlockUpdateRhs','Creating rhs for component: '//TRIM(I2S(NoRow)),Level=12)
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
        ! befor computing the bnorm used to estimate the convergence.
        IF( NoCol == NoRow ) CYCLE
        
        Var => BlockMatrix % SubVector(NoCol) % Var
        x => Var % Values

        DoSum = .FALSE.
        A => BlockMatrix % SubMatrix( NoRow, NoCol ) % Mat
        IF( A % NumberOfRows > 0 ) THEN        
          CALL CRS_MatrixVectorMultiply( A, x, rtmp)              
          DoSum = .TRUE.
        ELSE IF( BlockMatrix % SubMatrixTranspose(NoCol,NoRow) ) THEN
          A => TotMatrix % SubMatrix(NoCol,NoRow) % Mat
          IF( A % NumberOfRows >  0 ) THEN
            CALL CRS_TransposeMatrixVectorMultiply( A, x, rtmp )
            DoSum = .TRUE.
          END IF
        END IF
        IF(DoSum) rhs(1:n) = rhs(1:n) - rtmp(1:n) 
      END DO

      
      bnorm = SQRT( SUM( rhs**2 ) / n )
      BlockMatrix % SubVector(NoRow) % bnorm = bnorm
      
      ! Finally duduct the diagonal entry so that we can solve for the residual
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
    
    INTEGER :: i,j,k,NoVar,i1,i2,j1,j2
    REAL(KIND=dp), ALLOCATABLE :: s(:)
    INTEGER :: maxsize,ndofs
    INTEGER, POINTER :: Offset(:)
    REAL(KIND=dp) :: nrm
    TYPE(Matrix_t), POINTER :: A

    REAL(KIND=dp), POINTER :: b(:)
    LOGICAL :: DoSum 
    
    CALL Info('BlockMatrixVectorProd','Starting block matrix multiplication',Level=20)
    
    NoVar = TotMatrix % NoVar
    Offset => TotMatrix % Offset
    MaxSize = TotMatrix % MaxSize

    IF (isParallel) THEN
      IF(.NOT.ASSOCIATED(TotMatrix % SubMatrix(1,NoVar) % Mat % ParMatrix)) THEN
        IF(.NOT.ASSOCIATED(SolverMatrix)) THEN
          CALL Fatal('BlockMatrixVectorProd','No matrix to apply.')
        ELSE
          CALL ParallelMatrixVector(SolverMatrix, u(1:ipar(3)), v(1:ipar(3)))
        END IF
        RETURN
      END IF
    END IF

    ALLOCATE( s(MaxSize) )

    v(1:offset(NoVar+1)) = 0
    

    DO i=1,NoVar
      DO j=1,NoVar
        s = 0._dp

        j1 = offset(j)+1
        j2 = offset(j+1)

        DoSum = .FALSE.
        
        A => TotMatrix % SubMatrix(i,j) % Mat
        IF( A % NumberOfRows > 0 ) THEN
          
          IF( MAXVAL( A % Cols ) > offset(j+1)-offset(j) ) THEN
            CALL Fatal('BlockMatrixVectorProd','Wrong max column index: '&
                //TRIM(I2S(MAXVAL( A % Cols )))//' vs. '//TRIM(I2S(offset(j+1)-offset(j))))
          END IF
          IF( A % NumberofRows > offset(i+1)-offset(i)) THEN         
            CALL Fatal('BlockMatrixVectorProd','Wrong max column index: '&
                //TRIM(I2S( A % NumberOfRows ))//' vs. '//TRIM(I2S(offset(i+1)-offset(i))))
          END IF
          
          CALL Info('BlockMatrixVectorProd','Multiplying with submatrix ('&
              //TRIM(I2S(i))//','//TRIM(I2S(j))//')',Level=8)          
          
          IF (isParallel) THEN
            CALL ParallelMatrixVector( A, u(j1:j2), s  )
          ELSE
            CALL CRS_MatrixVectorMultiply( A, u(j1:j2), s )
          END IF
          DoSum = .TRUE.
          
        ELSE IF( TotMatrix % SubMatrixTranspose(j,i) ) THEN          
          A => TotMatrix % SubMatrix(j,i) % Mat
          IF( A % NumberOfRows > 0 ) THEN            
            CALL Info('BlockMatrixVectorProd','Multiplying with transpose of submatrix ('&
                //TRIM(I2S(j))//','//TRIM(I2S(i))//')',Level=12)
            
            IF (isParallel) THEN
              CALL Fatal('BlockMatrixVectorProd','Do transpose in parallel!')
            ELSE          
              CALL CRS_TransposeMatrixVectorMultiply( A, u(j1:j2), s )
            END IF
            DoSum = .TRUE.
          END IF
        END IF
          
        
        IF( InfoActive( 15 ) ) THEN
          PRINT *,'MatVecProdNorm u:',i,j,&
              SQRT(SUM(u(j1:j2)**2)),SUM( u(j1:j2) ), MINVAL( u(j1:j2) ), MAXVAL( u(j1:j2) ) 
          PRINT *,'MatVecProdNorm s:',i,j,&
              SQRT(SUM(s**2)), SUM( s ), MINVAL( s ), MAXVAL( s ) 
        END IF

        IF( DoSum ) THEN
          DO k=1,offset(i+1)-offset(i)
            v(k+offset(i)) = v(k+offset(i)) + s(k)
          END DO
        END IF
      END DO
      
      IF( InfoActive( 15 ) ) THEN
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
      
    IF( InfoActive( 15 ) ) THEN
      i = offset(NoVar+1)
      nrm = SQRT( SUM( v(1:i)**2 ) )
      PRINT *,'MatVecProdNorm full:',nrm, MINVAL( v(1:i) ), MAXVAL( v(1:i)), SUM( v(1:i) )
    END IF
      
    CALL Info('BlockMatrixVectorProd','Finished block matrix multiplication',Level=20)
!------------------------------------------------------------------------------
  END SUBROUTINE BlockMatrixVectorProd
!------------------------------------------------------------------------------
  

!> Create the vectors needed for block matrix scaling. Currently only
!> real and complex valued row equilibriation is supported. Does not perform
!> the actual scaling.
!------------------------------------------------------------------------------
  SUBROUTINE CreateBlockMatrixScaling( )
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n,m,NoVar
    REAL(KIND=dp) :: nrm, tmp, blocknrm
    TYPE(Matrix_t), POINTER :: A, Atrans
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    LOGICAL :: ComplexMatrix, GotIt, DiagOnly
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: Found !, IsComplex
    TYPE(ValueList_t), POINTER :: Params
    
    
    CALL Info('CreateBlockMatrixScaling','Starting block matrix row equilibriation',Level=10)
    
    NoVar = TotMatrix % NoVar
    
    Params => CurrentModel % Solver % Values 
    DiagOnly = ListGetLogical( Params,'Block Scaling Diagonal',Found ) 
    IF( DiagOnly ) THEN
       CALL Info('CreateBlockMatrixScaling',&
            'Considering only diagonal matrices in scaling',Level=20)      
    END IF

    !IsComplex = ListGetLogical( Params,'Linear System Complex',Found ) 
    
    
    DO k=1,NoVar
      GotIt = .FALSE.
      A => TotMatrix % SubMatrix(k,k) % Mat
      IF( ASSOCIATED( A ) ) THEN
        IF( A % NumberOfRows > 0 ) THEN
          GotIt = .TRUE.
          ComplexMatrix = A % COMPLEX
        END IF
      END IF
      IF(.NOT. GotIt) CALL Warn('CreateBlockMatrixScaling','Improve complex matrix detection!')
        
      IF( ComplexMatrix ) THEN
        m = 2
        CALL Info('CreateBlockMatrixScaling',&
            'Assuming complex matrix block: '//TRIM(I2S(k)),Level=20)
      ELSE
        m = 1
        CALL Info('CreateBlockMatrixScaling',&
            'Assuming real valued matrix block: '//TRIM(I2S(k)),Level=20)
      END IF     
      
      n = TotMatrix % offset(k+1) - TotMatrix % offset(k)

      IF( .NOT. ALLOCATED( Totmatrix % SubVector(k) % DiagScaling ) ) THEN
        ALLOCATE( TotMatrix % SubVector(k) % DiagScaling(n) )
      END IF
      
      Diag => TotMatrix % SubVector(k) % DiagScaling
      Diag = 0.0_dp

      
      DO l=1,NoVar

        IF( DiagOnly ) THEN
          IF( k /= l ) CYCLE
        END IF
        
        A => TotMatrix % SubMatrix(k,l) % Mat

        IF( A % NumberOfRows ==  0 ) THEN
          IF( TotMatrix % SubMatrixTranspose(l,k) ) THEN                    
            CALL Info('CreateBlockMatrixScaling','Creating the transpose for real!')
            Atrans => CRS_Transpose( TotMatrix % SubMatrix(l,k) % Mat ) 
            TotMatrix % Submatrix(k,l) % Mat => Atrans
            TotMatrix % SubMatrixTranspose(l,k) = .FALSE.
            A => Atrans
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
        blocknrm = 0.0_dp
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

          blocknrm = MAX( blocknrm, tmp ) 
          
          ! Compute the sum to the real component, scaling for imaginary will be the same
          Diag(i) = Diag(i) + tmp
        END DO

        ! PRINT *,'BlockNorm:',k,l,blocknrm
        
      END DO
      
      IF (ParEnv % PEs > 1) THEN
        CALL ParallelSumVector(A, Diag)
      END IF
      
      nrm = MAXVAL( Diag(1:n) ) 

      IF( ParEnv % PEs > 1 ) THEN
        nrm = ParallelReduction(nrm,2)
      END IF
      
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

    END DO
    
    WRITE( Message, * ) 'Unscaled matrix norm: ', nrm    
    CALL Info('CreateBlockMatrixScaling', Message, Level=7 )
    
    CALL Info('CreateBlockMatrixScaling','Finished block matrix row equilibriation',Level=20)           
    
  END SUBROUTINE CreateBlockMatrixScaling
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixInfo()
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n,m,NoVar
    REAL(KIND=dp) :: nrm, tmp, blocknrm
    TYPE(Matrix_t), POINTER :: A, Atrans
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    LOGICAL :: ComplexMatrix, GotIt, DiagOnly
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: Found
    
    
    CALL Info('BlockMatrixInfo','Showing some ranges of block matrix stuff',Level=10)
    
    NoVar = TotMatrix % NoVar

    PRINT *,'BlockInfo:',NoVar
    
    
    DO k=1,NoVar
      DO l=1,NoVar
        
        A => TotMatrix % SubMatrix(k,l) % Mat
        IF( .NOT. ASSOCIATED( A ) ) CYCLE
        IF( A % NumberOfRows == 0 ) CYCLE
        
        n = TotMatrix % offset(k+1) - TotMatrix % offset(k)

        PRINT *,'BlockInfo:',k,l,A % NumberOfRows, n, A % COMPLEX

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
!> still use the optimal row equilibriation scaling for the block system. 
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixScaling( reverse, blockrow, blockcol, bext )
!------------------------------------------------------------------------------
    LOGICAL, OPTIONAL :: reverse
    INTEGER, OPTIONAL :: blockrow, blockcol
    REAL(KIND=dp), POINTER, OPTIONAL :: bext(:)

    INTEGER :: i,j,k,l,n,m,NoVar
    REAL(KIND=dp) :: nrm, tmp
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:), Diag(:), Values(:)
    INTEGER, POINTER :: Rows(:), Cols(:)
    LOGICAL :: backscale
    
    
    CALL Info('BlockMatrixScaling','Starting block matrix row equilibriation',Level=10)
           
    IF( PRESENT( Reverse ) ) THEN
      CALL Info('BlockMatrixScaling','Performing reverse scaling',Level=20)
      BackScale = Reverse
    ELSE
      CALL Info('BlockMatrixScaling','Performing forward scaling',Level=20)
      BackScale = .FALSE.
    END IF

    
    NoVar = TotMatrix % NoVar   
    DO k=1,NoVar
      
      IF( PRESENT( blockrow ) ) THEN
        IF( blockrow /= k ) CYCLE
      END IF
      
      n = TotMatrix % offset(k+1) - TotMatrix % offset(k)
      Diag => TotMatrix % SubVector(k) % DiagScaling
      IF( .NOT. ASSOCIATED( Diag ) ) THEN
        CALL Fatal('BlockMatrixScaling','Diag for scaling not associated!')
      END IF
      
      IF( BackScale ) Diag = 1.0_dp / Diag 
            
      DO l=1,NoVar        
        
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
        
    CALL Info('ForwardBlockMatrixScaling','Finished block matrix row equilibriation',Level=20)           

  END SUBROUTINE BlockMatrixScaling
!------------------------------------------------------------------------------


!> Deallocates the block matrix scaling vectors.   
!------------------------------------------------------------------------------
  SUBROUTINE DestroyBlockMatrixScaling()
!------------------------------------------------------------------------------
    INTEGER :: k,NoVar
    
    CALL Info('DestroyBlockMatrixScaling','Starting block matrix row equilibriation',Level=10)
              
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
    REAL(KIND=dp), TARGET, INTENT(out) :: u(*)
    REAL(KIND=dp), TARGET, INTENT(in) :: v(*)
    INTEGER :: ipar(*)
!---------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: rtmp(:),vtmp(:),xtmp(:),b(:), x(:), a_rhs_save(:)
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:)
    INTEGER :: i,j,k,l,NoVar
    TYPE(Solver_t), POINTER :: Solver, Solver_save, ASolver
    INTEGER, POINTER :: Offset(:)
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER :: BlockOrder(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, Var_save

    LOGICAL :: GotOrder, BlockGS, Found, NS, ScaleSystem, DoSum, IsComplex, BlockScaling
    CHARACTER(LEN=MAX_NAME_LEN) :: str
#ifndef USE_ISO_C_BINDINGS
    INTEGER(KIND=AddrInt) :: AddrFunc
#else
    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc
#endif

    CALL Info('BlockMatrixPrec','Starting block matrix preconditioning',Level=6)
    
    Solver => CurrentModel % Solver
    Params => Solver % Values
    
    ! Enable user defined order for the solution of blocks
    !---------------------------------------------------------------
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotOrder)
    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    
    NoVar = TotMatrix % NoVar
    Solver => TotMatrix % Solver
    offset => TotMatrix % Offset

    ! Save the initial solver stuff
    solver_save => Solver
    var_save => Solver % Variable
    mat_save => Solver % Matrix
    rhs_save => Solver % Matrix % RHS

    ! Always treat the inner iterations as truly complex if they are
    CALL ListAddLogical( Params,'Linear System Skip Complex',.FALSE.) 
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.FALSE.) 

    BlockScaling = ListGetLogical( Params,'Block Scaling',Found )
    
#define SOLSYS
#ifdef SOLSYS
    IF (isParallel) &
      ALLOCATE( x(TotMatrix % MaxSize), b(TotMatrix % MaxSize) )
#endif

    ! Initial guess 
    !-----------------------------------------
    u(1:offset(NoVar+1)) = v(1:offset(NoVar+1))

    IF( BlockGS ) THEN
      ALLOCATE( vtmp(offset(NoVar+1)), rtmp(offset(NoVar+1)), xtmp(offset(NoVar+1)))
      vtmp(1:offset(NoVar+1)) = v(1:offset(NoVar+1))
    END IF
    
    CALL ListPushNameSpace('block:')
    
    DO j=1,NoVar
      IF( GotOrder ) THEN
        i = BlockOrder(j)
      ELSE
        i = j
      END IF
      
      WRITE(Message,'(A,I0)') 'Solving block: ',i
      CALL Info('BlockMatrixPrec',Message,Level=6)

      ! Set pointers to the new linear system
      !-------------------------------------------------------------------
      Var => TotMatrix % SubVector(i) % Var

      IF( ASSOCIATED( TotMatrix % Subvector(i) % Solver ) ) THEN
        ASolver => TotMatrix % SubVector(i) % Solver
        A => ASolver % Matrix
      ELSE
        A => TotMatrix % Submatrix(i,i) % PrecMat
        IF( A % NumberOfRows == 0 ) THEN
          A => TotMatrix % Submatrix(i,i) % Mat
        ELSE
          PRINT *,'Using specialized preconditioning block'
        END IF      
        ASolver => Solver
      END IF
      
        
#ifndef SOLSYS
      x => u(offset(i)+1:offset(i+1))
      IF( BlockGS ) THEN
        b => vtmp(offset(i)+1:offset(i+1))
      ELSE      
        b => v(offset(i)+1:offset(i+1))
      END IF
#else
      IF (isParallel) THEN
        x = 0; b=0
        l = 0
        DO k=1,A % NumberofRows
          IF (Parenv % MyPE==A % ParallelInfo % NeighbourList(k) % Neighbours(1)) THEN
            l = l+1
            x(k) = u(offset(i)+l)
            IF( BlockGS ) THEN
              b(k) = vtmp(offset(i)+l)
            ELSE      
              b(k) = v(offset(i)+l)
            END IF
           END IF
        END DO
      ELSE
        x => u(offset(i)+1:offset(i+1))
        IF( BlockGS ) THEN
          b => vtmp(offset(i)+1:offset(i+1))
        ELSE      
          b => v(offset(i)+1:offset(i+1))
        END IF
      END IF
#endif
      !A_rhs_save => A_rhs	
      !A % RHS => b

      ! Reuse block preconditioner from the first block to other components
      !--------------------------------------------------------------------
      IF( ListGetLogical( Params,'Block Prec Reuse',Found) ) THEN
        DO k = 1, NoVar
          IF( k == i ) CYCLE
          IF( CRS_CopyMatrixPrec( TotMatrix % Submatrix(k,k) % Mat, A ) ) EXIT
        END DO
      END IF

      CALL ListPushNameSpace('block '//TRIM(i2s(i))//TRIM(i2s(i))//':')

      ! We do probably not want to compute the change within each iteration
      CALL ListAddNewLogical( Asolver % Values,'Skip Compute Nonlinear Change',.TRUE.)         
        
      ! Revert back if the matrix was set not complex
      !IF( ListGetLogical( ASolver % Values,'Linear System Complex', Found ) ) A % COMPLEX = .TRUE.
        
      IF( InfoActive( 15 ) ) THEN
        PRINT *,'Range pre:',i,TRIM(Var % Name), MINVAL(x),MAXVAL(x)
      END IF

      IF( BlockScaling ) CALL BlockMatrixScaling(.TRUE.,k,k,b)


      IF( InfoActive( 15 ) ) THEN
        CALL BlockMatrixInfo()
      END IF

      
      IF (isParallel) THEN
#ifndef SOLSYS
        GlobalData => A % ParMatrix
        GlobalMatrix => GlobalData % SplittedMatrix % InsideMatrix
        GlobalMatrix % MatVecSubr = A % MatVecSubr
        GlobalMatrix % Ematrix => A
        GlobalMatrix % COMPLEX = A % COMPLEX
 
        CALL IterSolver( GlobalMatrix, x,b, &
            ASolver,MatvecF=AddrFunc(SParMatrixVector), &
            DotF=AddrFunc(SParDotProd), NormF=AddrFunc(SParNorm))
#else
        !CALL SolveSystem( A, ParMatrix, b, x, Var % Norm, Var % DOFs, ASolver )
        CALL SolveLinearSystem( A, b, x, Var % Norm, Var % DOFs, ASolver )
#endif
      ELSE
        
        !ScaleSystem = ListGetLogical( Params,'block: Linear System Scaling', Found )
        !IF(.NOT. Found) ScaleSystem = .TRUE.
        !IF ( ScaleSystem ) CALL ScaleLinearSystem(ASolver, A,b,x )
        
        CALL SolveLinearSystem( A, b, x, Var % Norm, Var % DOFs, ASolver )

        !IF( ScaleSystem ) CALL BackScaleLinearSystem(ASolver,A,b,x)       
      END IF

      IF( BlockScaling ) CALL BlockMatrixScaling(.FALSE.,k,k,b)

      
      IF( InfoActive( 15 ) ) THEN
        PRINT *,'Range post:',i,TRIM(Var % Name), MINVAL(x),MAXVAL(x)
      END IF
        

#ifdef SOLSYS
      IF (isParallel) THEN
        l = 0
        DO k=1,A % NumberofRows
          IF (Parenv % MyPE==A % ParallelInfo % NeighbourList(k) % Neighbours(1)) THEN
            l = l+1
            x(l) = x(k)
            u(offset(i)+l) = x(l)
          END IF
        END DO
      END IF
#endif

 
      !---------------------------------------------------------------------
      IF( BlockGS ) THEN        
        CALL Info('BlockMatrixPrec','Updating block r.h.s',Level=5)
      
        DO l=j+1,NoVar
          IF( GotOrder ) THEN
            k = BlockOrder(l)
          ELSE
            k = l
          END IF

          WRITE( str,'(A,I0,I0)') 'Block Gauss-Seidel Passive ',k,i
          IF( ListGetLogical( Params, str, Found ) ) CYCLE
        
          CALL Info('BlockMatrixPrec','Updating r.h.s for component '//TRIM(I2S(k)),Level=15)

          !IF( ASSOCIATED( TotMatrix % Subvector(i) % Solver ) ) THEN
          !  ASolver => TotMatrix % SubVector(i) % Solver
          !  A => ASolver % Matrix
          !ELSE
          !  A => TotMatrix % Submatrix(i,i) % Mat
          !END IF
          
            ! The residual is used only as a temporary vector
          !-------------------------------------------------------------
          DoSum = .FALSE.

          IF (isParallel) THEN
            IF(ASSOCIATED(TotMatrix % SubMatrix(k,i) % Mat % ParMatrix)) THEN
              CALL ParallelMatrixVector(TotMatrix % SubMatrix(k,i) % Mat,x,rtmp )
              DoSum = .TRUE.
            ELSE IF (ASSOCIATED(SolverMatrix)) THEN
              xtmp=0
              xtmp(offset(i)+1:offset(i+1))=x(offset(i)+1:offset(i+1))
              CALL ParallelMatrixVector(SolverMatrix,xtmp,rtmp)
              rtmp(1:offset(k+1)-offset(k))=rtmp(offset(k)+1:offset(k+1))
              DoSum = .TRUE.
            ELSE
              CALL Fatal('BlockMatrixPrec','No matrix to apply.')
            END IF

            IF( TotMatrix % SubMatrixTranspose(i,k) ) THEN
              CALL Fatal('','Do transpose in parallel!')
            END IF

          ELSE
            !IF( ASSOCIATED( TotMatrix % Subvector(i) % Solver ) ) THEN
            !  ASolver => TotMatrix % SubVector(i) % Solver
            !  A => ASolver % Matrix
            !ELSE
            !  A => TotMatrix % Submatrix(i,i) % Mat
            !END IF

            IF(.NOT.ASSOCIATED(TotMatrix % SubMatrix(k,i) % Mat) ) THEN
            !  CALL Fatal('','what is this')
            !  xtmp=0
            !  xtmp(offset(i)+1:offset(i+1))=x(offset(i)+1:offset(i+1))
            !  CALL CRS_MatrixVectorMultiply(SolverMatrix,x,rtmp )
            !  rtmp(1:offset(k+1)-offset(k))=rtmp(offset(k)+1:offset(k+1))
            END IF
            A => TotMatrix % Submatrix(k,i) % Mat            
            IF( A % NumberOfRows > 0 ) THEN
              CALL CRS_MatrixVectorMultiply(A,x,rtmp )
              DoSum = .TRUE.
            ELSE IF( TotMatrix % SubMatrixTranspose(i,k) ) THEN
              A => TotMatrix % SubMatrix(i,k) % Mat
              IF( A % NumberOfRows >  0 ) THEN
                CALL Info('BlockMatrixPrec','Multiplying with transpose of submatrix ('&
                    //TRIM(I2S(i))//','//TRIM(I2S(k))//')',Level=8)
                CALL CRS_TransposeMatrixVectorMultiply( A, x, rtmp )
                DoSum = .TRUE.
              END IF
            END IF
          END IF

          IF( DoSum ) THEN
            ! Up-date the off-diagonal entries to the r.h.s. 
            vtmp(offset(k)+1:offset(k+1)) = vtmp(offset(k)+1:offset(k+1)) &
                - rtmp(1:offset(k+1)-offset(k))
          END IF                 

        END DO ! l=j+1,NoVar
        
      END IF

      CALL ListPopNameSpace() ! block ij:
      
    END DO ! j=1,NoVar

#ifdef SOLSYS
    IF (isParallel) DEALLOCATE(x,b)
#undef SOLSYS
#endif

    CALL ListPopNameSpace('block:') ! block:

    CALL ListAddLogical( Params,'Linear System Refactorize',.FALSE. )

      
    Solver => Solver_save
    Solver % Matrix => mat_save
    Solver % Matrix % RHS => rhs_save
    Solver % Variable => Var_save

    IF( BlockGS ) THEN
      DEALLOCATE( vtmp, rtmp ) 
    END IF

    CALL Info('BlockMatrixPrec','Finished block matrix preconditioning',Level=6)
    
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
    
    !mat_Save => Solver % Matrix
    !Solver % Matrix => A
    !rhs_save => Solver % Matrix % RHS
    !Solver % Matrix % RHS => b

    CALL ListPushNamespace('block:')

    ! We don't want compute change externally
    CALL ListAddNewLogical( Params,'Skip compute nonlinear change',.TRUE.)
    
    Relax = 1.0_dp
    
    DO iter = 1, LinIter
      
      ! In block Jacobi the r.h.s. is not updated during the iteration cycle
      !----------------------------------------------------------------------
      IF( BlockGS ) THEN
        WRITE( Message,'(A,I0)') 'Block Gauss-Seidel iteration: ',iter
      ELSE
        WRITE( Message,'(A,I0)') 'Block Jacobi iteration: ',iter
        CALL BlockUpdateRhs(TotMatrix)
      END IF
      CALL Info('BlockSolver',Message,Level=6)
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
          PRINT *,'rhs'//TRIM(I2S(i))//':',SQRT( SUM(b**2) ), MINVAL( b ), MAXVAL( b ), SUM( b )
        END IF

        Var => TotMatrix % SubVector(RowVar) % Var
        Solver % Variable => Var
        
        A => TotMatrix % Submatrix(i,i) % PrecMat
        IF( A % NumberOfRows == 0 ) THEN
          A => TotMatrix % Submatrix(i,i) % Mat
        ELSE
          CALL Info('BlockSolver','Using preconditioning block: '//TRIM(I2S(i)))
        END IF
        
        !Solver % Matrix => A

        ! Use the newly computed residual rather than original r.h.s. to solve the equation!!
        rhs_save => A % rhs ! Solver % Matrix % RHS
        A % RHS => b
        
        ! Solving the subsystem
        !-----------------------------------
        ALLOCATE( dx( SIZE( Var % Values ) ) )
        dx = 0.0_dp

        CALL ListPushNamespace('block '//TRIM(i2s(RowVar))//TRIM(i2s(RowVar))//':')          

        IF( BlockScaling ) CALL BlockMatrixScaling(.TRUE.,i,i,b)
              
        !IF( ListGetLogical( Solver % Values,'Linear System Complex', Found ) ) A % Complex = .TRUE.

        !ScaleSystem = ListGetLogical( Solver % Values,'block: Linear System Scaling', Found )
        !IF(.NOT. Found) ScaleSystem = .TRUE.

        !IF ( ScaleSystem ) CALL ScaleLinearSystem(Solver, A, b, dx )
        
        CALL SolveLinearSystem( A, b, dx, Var % Norm, Var % DOFs, Solver )

        IF( BlockScaling ) CALL BlockMatrixScaling(.FALSE.,i,i,b)

        
        !IF( ScaleSystem ) CALL BackScaleLinearSystem(Solver, A, b, dx)       

        CALL ListPopNamespace()

        ! Revert back to original r.h.s.
        A % RHS => rhs_save
        !Solver % Matrix => mat_save

        IF( iter > 1 ) THEN
          Var % Values = Var % Values + Relax * dx
        ELSE
          Var % Values = Var % Values + dx
        END IF
          
        dxnorm = SQRT( SUM(dx**2) )
        xnorm = SQRT( SUM( Var % Values**2 ) )

        Var % Norm = xnorm
        Var % NonlinChange = dxnorm / xnorm
        
        IF( InfoActive( 15 ) ) THEN
          PRINT *,'dx'//TRIM(I2S(i))//':',SQRT( SUM(dx**2) ), MINVAL( dx ), MAXVAL( dx ), SUM( dx ), SUM( ABS( dx ) )
        END IF
      
        DEALLOCATE( dx )
          
        TotNorm = TotNorm + Var % Norm
        MaxChange = MAX( MaxChange, Var % NonlinChange )
      END DO

      IF( InfoActive( 15 ) ) THEN
        PRINT *,'GS Norm:',iter, MaxChange, TotNorm
      END IF

      IF( MaxChange < LinTol .AND. iter >= MinIter ) EXIT

    END DO
    CALL ListPopNamespace('block:')

    !Solver % Matrix % RHS => rhs_save
    !Solver % Matrix => mat_save

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
#ifdef USE_ISO_C_BINDINGS
    EXTERNAL :: AddrFunc
#endif
    INTEGER(KIND=AddrInt) :: iterProc,precProc, mvProc,dotProc,nmrProc, zero=0
    REAL(KIND=dp) :: dpar(20), xnorm,prevxnorm
    REAL(KIND=dp), ALLOCATABLE :: x(:),b(:),r(:)
    
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: SVar
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoVar, ndim, maxsize
    LOGICAL :: Converged, Diverged
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    INTEGER, POINTER :: Offset(:),poffset(:),BlockStruct(:)
    INTEGER :: i,j,k,l,ia,ib,istat
    LOGICAL :: LS, BlockAV,Found

    CALL Info('BlockKrylovIter','Starting block system iteration',Level=8)
    
    !CALL ListPushNameSpace('outer:')
    Params => Solver % Values
    
    BlockAV = ListGetLogical(Params,'Block A-V System', Found)

    Offset => TotMatrix % Offset
    ndim = TotMatrix % TotSize 
    NoVar = TotMatrix % NoVar

    CALL Info('BlockKrylovIter','Allocating temporal vectors for block system of size: '&
        //TRIM(I2S(ndim)),Level=8)

    ALLOCATE(x(ndim), b(ndim),r(ndim),STAT=istat)
    IF( istat /= 0 ) THEN
      CALL Fatal('BlockKrylovIter','Cannot allocate temporal vectors of size: '//TRIM(I2S(ndim)))
    END IF
    
    x=0;b=0;r=0
    
    IF (isParallel) THEN
      CALL Info('BlockKrylovIter','Preforming parallel initializations!',Level=18)
      DO i=1,NoVar
        DO j=1,NoVar
          IF ( i /= j ) THEN
            IF(ASSOCIATED(TotMatrix % SubMatrix(i,j) % Mat % ParMatrix)) &
              CALL ParallelInitSolve( TotMatrix % SubMatrix(i,j) % Mat,r,r,r)
          ELSE 
            x(offset(i)+1:offset(i+1)) = TotMatrix % SubVector(i) % Var % Values
            b(offset(i)+1:offset(i+1)) = TotMatrix % SubMatrix(i,i) % Mat % RHS
            
            CALL ParallelInitSolve( TotMatrix % SubMatrix(i,j) % Mat, &
                x(offset(j)+1:offset(j+1)), b(offset(i)+1:offset(i+1)),r )
          END IF
        END DO
      END DO
      IF(ASSOCIATED(SolverMatrix)) CALL ParallelInitSolve( SolverMatrix, x, b, r )
    END IF
    
    k = 0
    x = 0
    b = 0

    IF (isParallel) ALLOCATE(poffset(NoVar+1))

    CALL Info('BlockKrylovIter','Initializing monolithic system vectors',Level=18)
    
    DO i=1,NoVar

      IF( .NOT. ASSOCIATED( TotMatrix % Subvector(i) % Var ) ) THEN
        CALL Fatal('BlockKrylovIter','Subvector '//TRIM(I2S(i))//' not associated!')
      END IF
      IF( .NOT. ASSOCIATED( TotMatrix % Submatrix(i,i) % Mat ) ) THEN
        CALL Fatal('BlockKrylovIter','Submatrix '//TRIM(I2S(i))//' not associated!')
      END IF
      IF( .NOT. ASSOCIATED( TotMatrix % Submatrix(i,i) % Mat % Rhs ) ) THEN
        CALL Warn('BlockKrylovIter','Submatrix rhs '//TRIM(I2S(i))//' not associated!')
      END IF
      
      IF (.NOT.isParallel) THEN
        x(offset(i)+1:offset(i+1)) = TotMatrix % SubVector(i) % Var % Values        

        IF( ASSOCIATED( TotMatrix % Submatrix(i,i) % Mat % Rhs ) ) THEN          
          b(offset(i)+1:offset(i+1)) = TotMatrix % SubMatrix(i,i) % Mat % rhs
        END IF
        
      ELSE
        A => TotMatrix % SubMatrix(i,i) % Mat
        poffset(i) = k
        DO j=1,offset(i+1)-offset(i)
          IF ( A % ParallelInfo % NeighbourList(j) % Neighbours(1) == ParEnv % Mype ) THEN
            k=k+1
            x(k) = TotMatrix % SubVector(i) % Var % Values(j)
          END IF
        END DO
        poffset(i+1) = k
      END IF
    END DO
    
    
    IF (isParallel) THEN
      DO i=1,NoVar
        A => TotMatrix % SubMatrix(i,i) % Mat
        b(poffset(i)+1:poffset(i+1)) = A % ParMatrix % SplittedMatrix % InsideMatrix % Rhs
      END DO
      
      ndim = poffset(NoVar+1)
      TotMatrix % Offset => poffset
    END IF
    
    !----------------------------------------------------------------------
    ! Solve matrix equation solver with the redefined block matrix operations
    !----------------------------------------------------------------------
    CALL ListAddLogical(Params,'Linear System Free Factorization',.FALSE.)

    precProc = AddrFunc(BlockMatrixPrec)
    mvProc = AddrFunc(BlockMatrixVectorProd)       

    prevXnorm = SQRT( SUM( b**2 ) )
    WRITE( Message,'(A,ES12.3)') 'Rhs norm at start: ',PrevXnorm
    CALL Info('BlockKrylovIter',Message,Level=10)

    prevXnorm = SQRT( SUM( x**2 ) )
    WRITE( Message,'(A,ES12.3)') 'Solution norm at start: ',PrevXnorm
    CALL Info('BlockKrylovIter',Message,Level=10)

    CALL Info('BlockKrylovIter','Start of blocks system iteration',Level=18)

    ! Always treat the block system as a real valued system and complex
    ! arithmetics only at the inner level.
    CALL ListAddLogical( Params,'Linear System Skip Complex',.TRUE.) 

    ! Skip the scaling for block level system as the default routines
    ! would not perform it properly.
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.TRUE.) 
     
    IF(ASSOCIATED(SolverMatrix)) THEN
      A => SolverMatrix
    ELSE
      A => TotMatrix % SubMatrix(1,1) % Mat
    END IF

    !IF( ListGetLogical( Solver % Values,'Linear System Complex', Found ) ) A % COMPLEX = .TRUE.

    IF (isParallel) THEN
        A => A % ParMatrix % SplittedMatrix % InsideMatrix
      CALL IterSolver( A,x,b,&
          Solver,ndim=ndim,MatvecF=mvProc,PrecF=precProc,&
          DotF=AddrFunc(SParDotProd), NormF=AddrFunc(SParNorm))
      
      IF (BlockAV) THEN
        IF(ASSOCIATED(SolverMatrix)) THEN
          SolverMatrix % ParMatrix % SplittedMatrix % TmpXvec = x(1:poffset(NoVar+1))
          CALL ParallelUpdateResult(SolverMatrix,x,r)
        ELSE
          CALL Fatal('BlockKrylovIter','No matrix to apply.')
        END IF
      ELSE
        DO i=1,NoVar
            TotMatrix % SubMatrix(i,i) % Mat % ParMatrix % SplittedMatrix % &
              TmpXvec = x(poffset(i)+1:poffset(i+1))
        END DO
      
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
    CALL info('BlockKrylovIter','Finished block system iteration',Level=18)
    
    CALL ListAddLogical(Params,'Linear System Refactorize',.TRUE.)
    CALL ListAddLogical(Params,'Linear System Free Factorization',.TRUE.)

    !CALL ListPopNamespace()
    Xnorm = SQRT( SUM( x**2) )
    
    WRITE( Message,'(A,ES12.3)') 'Solution norm: ',Xnorm
    CALL Info('BlockKrylovIter',Message,Level=8)


    MaxChange = 2*ABS(Xnorm-PrevXnorm)/(Xnorm+PrevXnorm)
    PrevXNorm = Xnorm

    WRITE( Message,'(A,ES12.3)') 'Relative change: ',MaxChange
    CALL Info('BlockKrylovIter',Message,Level=8)
    
    DO i=1,NoVar
      TotMatrix % SubVector(i) % Var % Values(1:offset(i+1)-offset(i)) = & 
          x(offset(i)+1:offset(i+1))
    END DO
    TotMatrix % Offset => Offset
      
    ! Copy values back since for nontrivial block-matrix structure the
    ! components do not build the whole solution.
    !-----------------------------------------------------------------
    IF( TotMatrix % GotBlockStruct ) THEN
      SVar => CurrentModel % Solver % Variable

      BlockStruct => TotMatrix % BlockStruct
      l = SIZE(BlockStruct) 
      DO j=1,l
        k = BlockStruct(j)
        ia = COUNT( k == BlockStruct(1:l) ) 
        ib = COUNT( k == BlockStruct(1:j) )
        
        DO i=1,SIZE(x)/l
          SVar % Values(l*(i-1)+j) = &
              TotMatrix % SubVector(k) % Var % Values( ia*(i-1) + ib )
        END DO
      END DO
    ELSE IF (BlockAV) THEN
      SolverVar % Values = x
    END IF

    CALL Info('BlockKrylovIter','Finished block krylov iteration',Level=20)
   
  END SUBROUTINE blockKrylovIter



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
    LOGICAL :: GotIt, GotIt2, BlockPrec, BlockGS, BlockJacobi, BlockAV, &
        BlockHorVer, BlockCart, BlockNodal
    INTEGER :: ColVar, RowVar, NoVar, BlockDofs, VarDofs
    
    REAL(KIND=dp) :: NonlinearTol, Norm, PrevNorm, Residual, PrevResidual, &
        TotNorm, MaxChange, alpha, beta, omega, rho, oldrho, s, r, PrevTotNorm, &
        Coeff
    REAL(KIND=dp), POINTER :: SaveValues(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
    CHARACTER(LEN=max_name_len) :: str, VarName, ColName, RowName
    LOGICAL :: Robust, LinearSearch, ErrorReduced, IsProcedure, ScaleSystem,&
        ReuseMatrix, LS, BlockScaling
    INTEGER :: HaveConstraint, HaveAdd
    INTEGER, POINTER :: VarPerm(:)
    INTEGER, POINTER :: SlaveSolvers(:)
    LOGICAL :: GotSlaveSolvers, SkipVar
    
    
    TYPE(Matrix_t), POINTER :: Amat, SaveMatrix, SaveCM
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params


    CALL Info('BlockSolverInt','---------------------------------------',Level=5)

    Params => Solver % Values
    Mesh => Solver % Mesh
    PSolver => Solver

    isParallel = ParEnv % PEs > 1
    
    
    ! Determine some parameters related to the block strategy
    !------------------------------------------------------------------------------
    BlockPrec = ListGetLogical( Params,'Block Preconditioner',GotIt)
    IF(.NOT. GotIt) THEN
      CALL Info('BlockSolver','Using block preconditioning mode by default')
      BlockPrec = .TRUE.
    END IF

    BlockScaling = ListGetLogical( Params,'Block Scaling',GotIt)

    ! Block iteration style: jacobi vs. gauss-seidel
    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',GotIt)    
    BlockJacobi = ListGetLogical( Params,'Block Jacobi',GotIt)

    ! Different strategies on how to split the initial monolithic matrix into blocks
    BlockAV = ListGetLogical( Params,'Block A-V System', GotIt)
    BlockNodal = ListGetLogical( Params,'Block Nodal System', GotIt)
    BlockHorVer = ListGetLogical( Params,'Block Hor-Ver System', GotIt)
    BlockCart = ListGetLogical( Params,'Block Cartesian System', GotIt)
    
    SlaveSolvers =>  ListGetIntegerArray( Params, &
         'Block Solvers', GotSlaveSolvers )

    SkipVar = .FALSE.
    IF( BlockAV .OR. BlockNodal .OR. BlockHorVer ) THEN
      BlockDofs = 2
      SkipVar = .TRUE.
    ELSE IF( BlockCart ) THEN
      BlockDofs = 3
      SkipVar = .TRUE.
    ELSE IF( GotSlaveSolvers ) THEN
      BlockDofs = SIZE( SlaveSolvers )
    ELSE
      BlockDofs = Solver % Variable % Dofs      
    END IF
    VarDofs = BlockDofs
    
    HaveConstraint = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  HaveConstraint = 1
    HaveConstraint = ParallelReduction(HaveConstraint*1._dp)
     
    HaveAdd = 0
    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows > 0 ) HaveAdd = 1
    END IF
    HaveAdd = ParallelReduction(HaveAdd*1._dp)

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
      IF( BlockAV ) THEN
        CALL BlockPickMatrixAV( Solver, VarDofs )
      ELSE IF( BlockHorVer .OR. BlockCart ) THEN
        CALL BlockPickMatrixHorVer( Solver, VarDofs, BlockCart )       
      ELSE IF( BlockNodal ) THEN
        CALL BlockPickMatrixNodal( Solver, VarDofs )        
      ELSE IF( VarDofs > 1 ) THEN
        CALL BlockPickMatrix( Solver, VarDofs )
      ELSE
        CALL Info('BlockSolver','Using the original matrix as the (1,1) block!',Level=10)
        TotMatrix % SubMatrix(1,1) % Mat => SolverMatrix        
      END IF

      IF( SkipVar ) THEN
        CALL BlockInitVar( Solver, TotMatrix )
      END IF

      CALL BlockPrecMatrix( Solver, VarDofs ) 
    END IF

    CALL FsiCouplingBlocks( Solver )
    
    IF( HaveConstraint > 0 ) THEN
      CALL BlockPickConstraint( Solver, VarDofs )
      ! Storing pointer to CM so that the SolverMatrix won't be treated with the constraints
      SaveCM => Solver % Matrix % ConstraintMatrix
      Solver % Matrix % ConstraintMatrix => NULL()
    END IF

    IF (isParallel) THEN
      DO RowVar=1,NoVar
        DO ColVar=1,NoVar
          CALL ParallelActive( .TRUE.)
          Amat => TotMatrix % SubMatrix(RowVar,ColVar) % Mat
          Amat % Comm = Solver % Matrix % Comm
          Parenv % ActiveComm = Amat % Comm
          Solver % Variable => TotMatrix % SubVector(ColVar) % Var

          IF(Amat % NumberOfRows>0) THEN
            IF(Amat % NumberOfRows==MAXVAL(Amat % Cols)) THEN
              IF (.NOT.ASSOCIATED(Amat % ParMatrix)) &          
                CALL ParallelInitMatrix(Solver,Amat)
            END IF

            IF(ASSOCIATED(Amat % ParMatrix )) THEN
              Amat % ParMatrix % ParEnv % ActiveComm = &
                Amat % Comm
              ParEnv = Amat % ParMatrix % ParEnv
            END IF
          END IF
        END DO
      END DO
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
    
    ! The case with one block is mainly for testing and developing features
    ! related to nonlinearity and assembly.
    !----------------------------------------------------------------------
    IF( NoVar == 1 ) THEN
      CALL Info('BlockSolverInt','Solving in standard manner',Level=6)
      
      Solver % Variable => TotMatrix % SubVector(1) % Var
      Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat
      
      TotNorm = DefaultSolve()
      MaxChange = Solver % Variable % NonlinChange 
      
    ELSE IF( BlockPrec ) THEN
      CALL Info('BlockSolverInt','Using block preconditioning strategy',Level=6)        
      CALL BlockKrylovIter( Solver, MaxChange )
    ELSE
      Solver % Variable => TotMatrix % SubVector(1) % Var
      Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat
      
      CALL Info('BlockSolverInt','Using block solution strategy',Level=6)
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

    IF( BlockHorVer .OR. BlockCart ) THEN
      CALL BlockBackCopyVar( Solver, TotMatrix )
    END IF
      
    CALL Info('BlockSolverInt','All done')
    CALL Info('BlockSolverInt','-------------------------------------------------',Level=5)

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

    ! Eliminate recursion for block solvers. 
    CALL ListAddLogical(Solver % Values,'Linear System Block Mode',.FALSE.)
    CALL BlockSolveInt(A,x,b,Solver)
    CALL ListAddLogical(Solver % Values,'Linear System Block Mode',.TRUE.)

!------------------------------------------------------------------------------
END SUBROUTINE BlockSolveExt
!------------------------------------------------------------------------------

!> \}
