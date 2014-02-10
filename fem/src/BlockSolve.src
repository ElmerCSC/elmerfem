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
  SUBROUTINE BlockInitMatrix( Solver, BlockMatrix, BlockDofs )
    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    INTEGER :: BlockDofs
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix

    TYPE(Solver_t), POINTER :: PSolver
    INTEGER, POINTER :: BlockStruct(:)
    LOGICAL :: GotBlockStruct, Found
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER :: i,j,n,Novar
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(LEN=max_name_len) :: VarName, str

    Params => Solver % Values

    BlockMatrix => Solver % BlockMatrix
    IF (ASSOCIATED(BlockMatrix)) RETURN

    ALLOCATE(Solver % BlockMatrix)
    BlockMatrix => Solver % BlockMatrix
 
    BlockStruct => ListGetIntegerArray( Params,'Block Structure',GotBlockStruct)
    IF( GotBlockStruct ) THEN
      IF( SIZE( BlockStruct ) /= BlockDofs ) THEN
        CALL Fatal('BlockSolver','Incompatible size of > Block Structure < given!')
      END IF
      IF( MINVAL( BlockStruct ) < 1 .OR. MAXVAL( BlockStruct ) > BlockDofs ) THEN
        CALL Fatal('BlockSolver','Incompatible values in > Block Structure < given!')          
      END IF
      NoVar = MAXVAL( BlockStruct )
      BlockMatrix % BlockStruct => BlockStruct
    ELSE
      NoVar = BlockDofs
    END IF    
    BlockMatrix % GotBlockStruct = GotBlockStruct


    IF( BlockMatrix % NoVar == NoVar ) THEN
      CALL Info('InitializeBlockMatrix','Reusing existing blockmatrix',Level=6)
      RETURN
    ELSE IF( BlockMatrix % Novar /= 0 ) THEN
      CALL Fatal('InitializeBlockMatrix','Previous blockmatrix was of different size?')
    ELSE
      CALL Info('InitializeBlockMatrix','Starting',Level=6)
    END IF
    
    BlockMatrix % Solver => Solver
    BlockMatrix % NoVar = NoVar

    ALLOCATE( BlockMatrix % SubMatrix(NoVar,NoVar) )
    DO i=1,NoVar
      DO j=1,NoVar
        Amat => AllocateMatrix()
        CALL ClearMatrix(Amat)
        Amat % ListMatrix => NULL()
        Amat % FORMAT = MATRIX_LIST      
        Amat % NumberOfRows = 0
        BlockMatrix % Submatrix(i,j) % Mat => Amat
        
        ALLOCATE( BlockMatrix % Submatrix(i,j) % PrecMat )
        Amat => BlockMatrix % Submatrix(i,j) % PrecMat
        CALL ClearMatrix( Amat )
      END DO
    END DO
    
    ALLOCATE( BlockMatrix % SubMatrixActive(NoVar,NoVar) )
    BlockMatrix % SubMatrixActive = .FALSE.
    
    ALLOCATE( BlockMatrix % SubVector(NoVar))
    ALLOCATE( BlockMatrix % Offset(NoVar+1))
    BlockMatrix % Offset = 0
    BlockMatrix % maxsize = 0
    
    DO i = 1,NoVar
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
      
      !-----------------------------------------------------------------------------------
      ! If variable does not exist it will be created. 
      ! Here it is assumed that all components have the same number of dofs
      ! described by the same permutation vector. If the components are
      ! accounted in normal manner [1,2,3,...] then it suffices just to have 
      ! pointers to the components of the full vector.
      !-----------------------------------------------------------------------------------
      IF(ASSOCIATED( Var ) ) THEN
        CALL Info('BlockSolver','Using existing variable > '//TRIM(VarName)//' <')		
      ELSE		
        CALL Info('BlockSolver','Variable > '//TRIM(VarName)//' < does not exist, creating')
        PSolver => Solver
        IF( BlockMatrix % GotBlockStruct ) THEN
          j = COUNT( BlockMatrix % BlockStruct == i ) 
          IF( j == 0 ) THEN
            CALL Fatal('InitializeBlockMatrix','Invalid > Block Structure < given!')
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
    
    
  END SUBROUTINE BlockInitMatrix
    


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

    LOGICAL :: BlockAV
    INTEGER::i,j,k,i_aa,i_vv,i_av,i_va,n;
    TYPE(Matrix_t), POINTER :: B_aa,B_av,B_va,B_vv,C_aa,C_vv,A,CM

    SolverMatrix => Solver % Matrix 
    Params => Solver % Values

    BlockAV = ListGetLogical( Params,'Block A-V System', Found)
    IF(BlockAV) THEN
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

#if 1
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
if(c_vv%numberofrows<=0) b_vv%constraintmatrix=>null()
      END DO
#endif

      RETURN
    END IF

    ReuseMatrix = ListGetLogical( Params,'Block Matrix Reuse',Found)

    DO RowVar=1,NoVar
      DO ColVar=1,NoVar            
        Amat => TotMatrix % Submatrix(RowVar,ColVar) % Mat          
        IF( TotMatrix % GotBlockStruct) THEN
          ! A generic picking method for submatrices
          !----------------------------------------------------------------------
          CALL CRS_BlockMatrixPick2(SolverMatrix,Amat,TotMatrix % BlockStruct,RowVar,ColVar)
        ELSE
          ! Picking of standard submatrices of size one.
          !----------------------------------------------------------------------
          IF( ReuseMatrix ) THEN
            IF( RowVar + ColVar > 2 .AND. Amat % NumberOfRows == 0 ) THEN
              CALL CRS_CopyMatrixTopology( TotMatrix % Submatrix(1,1) % Mat, Amat )
            END IF
          END IF
          CALL CRS_BlockMatrixPick(SolverMatrix,Amat,NoVar,RowVar,ColVar)          
        END IF
      END DO
    END DO

  END SUBROUTINE BlockPickMatrix

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
    TYPE(Variable_t), POINTER :: Var
    
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
      END IF
      rhs => BlockMatrix % SubVector(NoRow) % rhs
      
      A => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat
      rhs = A % rhs
      
      DO NoCol = 1,NoVar           
        ! This ensures that the diagonal itself is not subtracted
        ! Otherwise the linear system should be solved for dx rather than x
        IF( NoCol == NoRow ) CYCLE
        
        Var => BlockMatrix % SubVector(NoCol) % Var
        x => Var % Values
        A => BlockMatrix % SubMatrix( NoRow, NoCol ) % Mat
        IF( A % NumberOfRows == 0 ) CYCLE
        
        rtmp = 0._dp
        CALL CRS_MatrixVectorMultiply( A, x, rtmp)      
        
        rhs(1:n) = rhs(1:n) - rtmp(1:n)
      END DO
    END DO
    
    DEALLOCATE( rtmp )
    
  END SUBROUTINE BlockUpdateRhs
  

  !------------------------------------------------------------------------------
  !> Perform matrix-vector product for block matrices.
  !> Has to be callable outside the module by Krylov methods.
  !------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixVectorProd( u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,NoVar,ipar(*)
    REAL(KIND=dp) :: u(*),v(*)
    REAL(KIND=dp), ALLOCATABLE :: s(:)
    INTEGER :: maxsize,ndofs
    INTEGER, POINTER :: Offset(:)
    
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
        IF (isParallel) THEN
          CALL ParallelMatrixVector( TotMatrix % SubMatrix(i,j) % Mat, &
                   u(offset(j)+1:offset(j+1)), s  )
        ELSE
          CALL CRS_MatrixVectorMultiply( TotMatrix % SubMatrix(i,j) % Mat, &
                   u(offset(j)+1:offset(j+1)), s )
        END IF

        DO k=1,offset(i+1)-offset(i)
          v(k+offset(i)) = v(k+offset(i)) + s(k)
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE BlockMatrixVectorProd
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> Perform block preconditioning by solving all the individual diagonal problems.
!> Has to be called outside the module by Krylov methods.
!------------------------------------------------------------------------------
  SUBROUTINE BlockMatrixPrec( u,v,ipar )    
    REAL(KIND=dp), TARGET :: u(*), v(*)
    INTEGER :: ipar(*)
!---------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: rtmp(:),vtmp(:),xtmp(:),b(:), x(:), rhs_save(:)
    INTEGER :: i,j,k,l,NoVar
    TYPE(Solver_t), POINTER :: Solver
    INTEGER, POINTER :: Offset(:)
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER :: BlockOrder(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, SolverVar

    LOGICAL :: GotOrder, BlockGS, Found, NS, SkipCompChange
    TYPE(Varying_string) :: namesp
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER(KIND=AddrInt) :: AddrFunc

    CALL Info('BlockMatrixPrec','Starting preconditioning',Level=6)
    
    Solver => CurrentModel % Solver
    Params => Solver % Values
    
    ! Enable user defined order for the solution of blocks
    !---------------------------------------------------------------
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotOrder)
    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    
    NoVar = TotMatrix % NoVar
    Solver => TotMatrix % Solver
    offset => TotMatrix % Offset
    SolverVar => Solver % Variable

    NS = ListGetNameSpace(namesp)

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
      A => TotMatrix % Submatrix(i,i) % PrecMat
      IF( A % NumberOfRows == 0 ) THEN
        A => TotMatrix % Submatrix(i,i) % Mat
      ELSE
        PRINT *,'Using specialized preconditioning block'
      END IF

      Var => TotMatrix % SubVector(i) % Var
      Solver % Variable => Var

      mat_save => Solver % Matrix
      Solver % Matrix => A

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
      rhs_save => Solver % Matrix % RHS
      Solver % Matrix % RHS => b
    
      ! Reuse block preconditioner from the first block to other components
      !--------------------------------------------------------------------
      IF( ListGetLogical( Params,'Block Prec Reuse',Found) ) THEN
        DO k = 1, NoVar
          IF( k == i ) CYCLE
          IF( CRS_CopyMatrixPrec( TotMatrix % Submatrix(k,k) % Mat, A ) ) EXIT
        END DO
      END IF

      CALL ListSetNameSpace('block '//TRIM(i2s(i))//TRIM(i2s(i))//':')
      SkipCompChange = ListGetLogical( Params,'Skip Compute Nonlinear Change',Found)
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)

      IF (isParallel) THEN
#ifndef SOLSYS
        GlobalData => A % ParMatrix
        GlobalMatrix => GlobalData % SplittedMatrix % InsideMatrix
        GlobalMatrix % MatVecSubr = A % MatVecSubr
        GlobalMatrix % Ematrix => A
        GlobalMatrix % COMPLEX = A % COMPLEX
 
        CALL IterSolver( GlobalMatrix, x,b, &
            Solver,MatvecF=AddrFunc(SParMatrixVector), &
                DotF=AddrFunc(SParDotProd), NormF=AddrFunc(SParNorm))
#else
        CALL SolveSystem( A, ParMatrix, b, x, Var % Norm, Var % DOFs, Solver )
#endif
      ELSE
        CALL SolveSystem( A, ParMatrix, b, x, Var % Norm, Var % DOFs, Solver )
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

      Solver % Matrix % RHS => rhs_save
      Solver % Matrix => mat_save

      !---------------------------------------------------------------------
      IF( BlockGS ) THEN        
        CALL Info('BlockSolver','Computing block r.h.s',Level=5)
      
        DO l=j+1,NoVar
          IF( GotOrder ) THEN
            k = BlockOrder(l)
          ELSE
            k = l
          END IF
        
          WRITE( str,'(A,I0,I0)') 'Block Gauss-Seidel Passive ',k,i
          IF( ListGetLogical( Params, str, Found ) ) CYCLE
        
          ! The residual is used only as a temporary vector
          !-------------------------------------------------------------
          IF (isParallel) THEN
            IF(ASSOCIATED(TotMatrix % SubMatrix(k,i) % Mat % ParMatrix)) THEN
              CALL ParallelMatrixVector(TotMatrix % SubMatrix(k,i) % Mat,x,rtmp )
            ELSE IF (ASSOCIATED(SolverMatrix)) THEN
              xtmp=0
              xtmp(offset(i)+1:offset(i+1))=x(offset(i)+1:offset(i+1))
              CALL ParallelMatrixVector(SolverMatrix,xtmp,rtmp)
              rtmp(1:offset(k+1)-offset(k))=rtmp(offset(k)+1:offset(k+1))
            ELSE
              CALL Fatal('BlockKrylovIter','No matrix to apply.')
            END IF
          ELSE
            IF(.NOT.ASSOCIATED(TotMatrix % SubMatrix(k,i) % Mat) ) THEN
              xtmp=0
              xtmp(offset(i)+1:offset(i+1))=x(offset(i)+1:offset(i+1))
              CALL CRS_MatrixVectorMultiply(SolverMatrix,x,rtmp )
              rtmp(1:offset(k+1)-offset(k))=rtmp(offset(k)+1:offset(k+1))
            ELSE
              CALL CRS_MatrixVectorMultiply(TotMatrix % SubMatrix(k,i) % Mat,x,rtmp )
            END IF
          END IF
          vtmp(offset(k)+1:offset(k+1)) = vtmp(offset(k)+1:offset(k+1)) &
                   - rtmp(1:offset(k+1)-offset(k))
        END DO
      END IF
    END DO

#ifdef SOLSYS
    IF (isParallel) DEALLOCATE(x,b)
#undef SOLSYS
#endif


    CALL ListSetNameSpace(CHAR(namesp))

    CALL ListAddLogical( Params,'Linear System Refactorize',.FALSE. )
    CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',SkipCompChange)
    Solver % Variable => SolverVar

    IF( BlockGS ) THEN
      DEALLOCATE( vtmp, rtmp ) 
    END IF

    CALL Info('BlockMatrixPrec','Finished preconditioning',Level=6)
    
  END SUBROUTINE BlockMatrixPrec


  !-----------------------------------------------------------------
  !> This call takes care of Jacobi & Gauss Seidel block methods. 
  !-----------------------------------------------------------------
  SUBROUTINE BlockStandardIter( Solver, MaxChange )

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: MaxChange

    INTEGER :: i,j,NoVar,RowVar,iter,LinIter
    INTEGER, POINTER :: BlockOrder(:)
    LOGICAL :: GotIt, GotBlockOrder, SkipCompChange, BlockGS
    REAL(KIND=dp), POINTER :: b(:), rhs_save(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, SolverVar
    REAL(KIND=dp) :: LinTol, TotNorm
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found

    NoVar = TotMatrix % NoVar
    Params => Solver % Values
    SolverVar => Solver % Variable

    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotBlockOrder)
    LinIter = ListGetInteger( Params,'Linear System Max Iterations',GotIt)
    LinTol = ListGetConstReal( Params,'Linear System Convergence Tolerance',GotIt)

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
        
        A => TotMatrix % Submatrix(RowVar,RowVar) % Mat
        b => TotMatrix % SubVector(RowVar) % rhs
        Var => TotMatrix % SubVector(RowVar) % Var
        Solver % Variable => Var
        
        mat_Save => Solver % Matrix
        Solver % Matrix => A
        rhs_save => Solver % Matrix % RHS
        Solver % Matrix % RHS => b
        
        ! Solving the subsystem
        !-----------------------------------
        CALL ListSetNameSpace('block '//TRIM(i2s(RowVar))//TRIM(i2s(RowVar))//':')          
        CALL SolveSystem( A, ParMatrix, b, &
            Var % Values, Var % Norm, Var % DOFs, Solver )
        
        Solver % Matrix % RHS => rhs_save
        Solver % Matrix => mat_save
        
        TotNorm = TotNorm + Var % Norm
        MaxChange = MAX( MaxChange, Var % NonlinChange )
      END DO      
    END DO

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

    INTEGER(KIND=AddrInt) :: AddrFunc, iterProc,precProc, mvProc,dotProc,nmrProc, zero=0
    REAL(KIND=dp) :: dpar(20), xnorm,prevxnorm
    REAL(KIND=dp), ALLOCATABLE :: x(:),b(:),r(:)
    
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: SVar
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoVar, ndim, maxsize
    LOGICAL :: Converged, Diverged
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    INTEGER, POINTER :: Offset(:),poffset(:),BlockStruct(:)
    INTEGER :: i,j,k,l,ia,ib
    LOGICAL :: LS, BlockAV,Found
    TYPE(Varying_string) :: namesp

    Params => Solver % Values

    BlockAV = ListGetLogical(Params,'Block A-V System', Found)

    Offset => TotMatrix % Offset
    ndim = TotMatrix % TotSize 
    NoVar = TotMatrix % NoVar
    ALLOCATE(x(ndim), b(ndim),r(ndim))
    x=0;b=0;r=0
    
    IF (isParallel) THEN
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
    
    k=0
    x=0;b=0
    IF (isParallel) ALLOCATE(poffset(NoVar+1))

    DO i=1,NoVar
      IF (.NOT.isParallel) THEN
        x(offset(i)+1:offset(i+1)) = TotMatrix % SubVector(i) % Var % Values
        b(offset(i)+1:offset(i+1)) = TotMatrix % SubMatrix(i,i) % Mat % rhs
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

    LS = ListGetNameSpace(namesp)
    CALL ListSetNameSpace('outer:')
    
    prevXnorm = SQRT( SUM( x**2 ) )

    IF(ASSOCIATED(SolverMatrix)) THEN
      A => SolverMatrix
    ELSE
      A => TotMatrix % SubMatrix(1,1) % Mat
    END IF
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
      CALL IterSolver( A,x,b,&
          Solver,ndim=ndim,MatvecF=mvProc,PrecF=precProc )
    END IF
    
    CALL ListAddLogical(Params,'Linear System Refactorize',.TRUE.)
    CALL ListAddLogical(Params,'Linear System Free Factorization',.TRUE.)

    CALL ListSetNameSpace(CHAR(namesp))
    Xnorm = SQRT( SUM( x**2) )
    
    MaxChange = 2*ABS(Xnorm-PrevXnorm)/(Xnorm+PrevXnorm)
    PrevXNorm = Xnorm
    
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
    
  END SUBROUTINE blockKrylovIter



!------------------------------------------------------------------------------
!> An alternative handle for the block solvers to be used by the legacy matrix
!> type. 
!------------------------------------------------------------------------------
  SUBROUTINE BlockSolveInt(A,x,b,Solver)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp), TARGET :: x(:),b(:)
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i,j,k,l,n,nd,NonLinIter,tests,NoTests,iter
    LOGICAL :: GotIt, GotIt2, BlockPrec, BlockGS, BlockAV
    INTEGER :: ColVar, RowVar, NoVar, BlockDofs

    REAL(KIND=dp) :: NonlinearTol, Norm, PrevNorm, Residual, PrevResidual, &
        TotNorm, MaxChange, alpha, beta, omega, rho, oldrho, s, r, PrevTotNorm, &
        Coeff
    REAL(KIND=dp), POINTER :: SaveValues(:), SaveRHS(:)
    CHARACTER(LEN=max_name_len) :: str, VarName, ColName, RowName
    LOGICAL :: Robust, LinearSearch, ErrorReduced, IsProcedure, ScaleSystem,&
        ReuseMatrix, LS
    INTEGER, POINTER :: VarPerm(:)
    
    TYPE(Matrix_t), POINTER :: Amat, SaveMatrix
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Varying_string) :: namesp

    CALL Info('BlockSolver','---------------------------------------',Level=5)

    Params => Solver % Values
    Mesh => Solver % Mesh
    PSolver => Solver

    isParallel = ParEnv % PEs > 1
    
         
    ! Determine some parameters related to the block strategy
    !------------------------------------------------------------------------------
    BlockPrec = ListGetLogical( Params,'Block Preconditioner',GotIt)
    IF(.NOT. GotIt ) BlockPrec = .TRUE.
    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',GotIt)    

    BlockAV = ListGetLogical( Params,'Block A-V System', GotIt)

    IF(BlockAV) THEN
      BlockDofs = 2
    ELSE
      BlockDofs = Solver % Variable % Dofs      
    END IF
    CALL BlockInitMatrix( Solver, TotMatrix, BlockDofs )

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
      
    CALL BlockPickMatrix( Solver, NoVar )
    CALL BlockPrecMatrix( Solver, Novar ) 

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
    
    LS = ListGetNameSpace(namesp)
    CALL ListSetNameSpace('outer:')
    
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
      CALL Info('BlockSolverInt','Using block precontioning strategy',Level=6)        
      CALL BlockKrylovIter( Solver, MaxChange )
    ELSE
      Solver % Variable => TotMatrix % SubVector(1) % Var
      Solver % Matrix => TotMatrix % Submatrix(1,1) % Mat
      
      CALL Info('BlockSolverInt','Using block solution strategy',Level=6)
      CALL BlockStandardIter( Solver, MaxChange )
    END IF
    CALL ListSetNameSpace('')

    ! For legacy matrices do the backmapping 
    !------------------------------------------
    SolverMatrix % RHS => SaveRHS
    Solver % Matrix => SaveMatrix
    Solver % Variable => SolverVar
    Solver % Variable % Values => SaveValues
       
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
    CALL BlockSolveInt(A,x,b,Solver)
!------------------------------------------------------------------------------
END SUBROUTINE BlockSolveExt
!------------------------------------------------------------------------------

!> \}
