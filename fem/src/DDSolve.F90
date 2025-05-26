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

MODULE DDSolve

 USE ParallelUtils 
 USE Integration
 USE ListMatrix
 USE ElementUtils, ONLY : FreeMatrix
 USE MatrixAssembly, ONLY : AddToMatrixElement, ZeroRow, SetMatrixElement
 USE IterativeMethods, ONLY : PseudoZDotProd
 USE IterSolve, ONLY : IterSolver
 USE ElementDescription, ONLY : ElementInfo, EdgeElementInfo
 USE SolverUtils, ONLY : SolveLinearSystem, ScaleLinearSystem, &
     BackScaleLinearSystem, AMGXMatrixVectorMultiply, AMGXSolver, DiagonalMatrixSumming, &
     StructureCouplingAssembly, FSICouplingAssembly, SaveLinearSystem, &
     MassMatrixAssembly, VectorValuesRange, LaplaceMatrixAssembly
 USE MeshUtils, ONLY : SaveProjector   
 USE DefUtils, ONLY : DefaultSolve, GetElementDOFs, GetElementNodes, GetLogical
 
 IMPLICIT NONE

  TYPE(BlockMatrix_t), POINTER, SAVE :: TotMatrix

  LOGICAL, PRIVATE :: isParallel=.FALSE.

  TYPE(Solver_t), POINTER, SAVE :: SolverRef
  TYPE(Variable_t), POINTER :: SolverVar => Null()
  TYPE(Matrix_t), POINTER :: SolverMatrix => Null(), SaveMatrix


  INTEGER, ALLOCATABLE, SAVE :: OwnerPart(:)
  INTEGER, ALLOCATABLE, SAVE :: NoOwners(:)
  
  INTEGER, SAVE :: DDMode = 0 
  REAL(KIND=dp), SAVE :: RobinC
  
  
CONTAINS
  

  !-------------------------------------------------------------------
  !> This subroutine creates the missing component variables.
  !------------------------------------------------------------------
  SUBROUTINE DDInitVar( Solver, BlockMatrix, A, x )    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: x(:)
    
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Matrix_t), POINTER :: Amat
    INTEGER :: i,j,k,n,m,Novar,ni,ndir,NoRow
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    CHARACTER(:), ALLOCATABLE :: VarName, str
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: Vals(:)
    INTEGER, POINTER :: VarPerm(:)
    TYPE(Matrix_t), POINTER :: Ai
    REAL(KIND=dp), POINTER :: xi(:),Wsum(:)
    REAL(KIND=dp) :: s, dval
    
    LOGICAL :: AddVector
    
    Params => Solver % Values
    Mesh => Solver % Mesh
    NoVar = BlockMatrix % NoVar

    IF(.NOT. ASSOCIATED(BlockMatrix % SubVector)) THEN
      ALLOCATE(BlockMatrix % SubVector(NoVar))
    END IF
    
    !A => Solver % Matrix
    !x => Solver % Variable % Values
    n = A % NumberOfRows


    ! Create "OnwerPart" such that maximum partition always owns the dof.
    IF(.NOT. ALLOCATED(OwnerPart) ) THEN
      ALLOCATE(OwnerPart(n),NoOwners(n))
      OwnerPart = 0
      NoOwners = 0

      ALLOCATE(Wsum(n))
      Wsum = 0.0_dp

      DO NoRow = 1,NoVar
        Ai => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat        
        ni = Ai % NumberOfRows 

        IF(.NOT. ASSOCIATED(Ai % InvPerm)) THEN
          ALLOCATE(Ai % InvPerm(ni))
          Ai % InvPerm = 0
          DO j=1,n
            k = Ai % Perm(j)
            IF(k>0) Ai % InvPerm(k) = j
          END DO
        END IF
                
        OwnerPart(Ai % InvPerm) = NoRow
        NoOwners(Ai % InvPerm) = NoOwners(Ai % InvPerm) + 1
        
        IF(.NOT. ALLOCATED( BlockMatrix % SubVector(NoRow) % rhs )) THEN
          ALLOCATE( BlockMatrix % SubVector(NoRow) % rhs(ni) )
          BlockMatrix % SubVector(NoRow) % rhs(ni) = 0.0_dp
        END IF

        ! Sum up all the weights.
        Wsum(Ai % InvPerm) = Wsum(Ai % InvPerm) + Ai % Rhs
      END DO

      DO NoRow = 1,NoVar
        Ai => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat               
        ! Normalize rhs such that it presents the fraction how to distribute it. 
        Ai % rhs = Ai % Rhs / Wsum(Ai % InvPerm)
        !PRINT *,'Owned dofs:',NoRow,COUNT(OwnerPart==NoRow),MINVAL(Ai % Rhs),MAXVAL(Ai % Rhs)
      END DO
      !PRINT *,'Number of Owners range:',MINVAL(NoOwners),MAXVAL(NoOwners)
      DEALLOCATE(Wsum)
    END IF
       
    
    DO NoRow=1,NoVar
      Ai => BlockMatrix % Submatrix(NoRow,NoRow) % Mat 
      ni = Ai % NumberOfRows
      
      Var => BlockMatrix % SubVector(NoRow) % Var 
      IF(.NOT. ASSOCIATED(Var) ) THEN
        ALLOCATE(Var)
        ALLOCATE(Var % Values(ni))
        Var % Dofs = 1
        Var % Values = 0.0_dp
        BlockMatrix % SubVector(NoRow) % Var => Var
      END IF

      ! Copy the initial guess. 
      xi => Var % Values
      xi = x(Ai % InvPerm)      

      ! If we set Dirichlet conditions we need original BulkValues to compute loads.
      IF(.NOT. ASSOCIATED(Ai % BulkValues)) THEN
        ALLOCATE( Ai % BulkValues(SIZE(Ai % Values)))
        Ai % BulkValues = Ai % Values
      END IF

      ! Set the Dirichlet conditions - the matrix side only!
      ! We assume that the loads are distributed evenly among participants. 
      ndir = 0
      DO j=1,ni
        k = Ai % InvPerm(j)        
        IF(A % ConstrainedDOF(k) ) THEN
          ! We scale the diagonal so that we can scale any r.h.s. load vector. 
          s = Ai % Rhs(j) * A % Values(A % Diag(k)) 
          CALL ZeroRow(Ai,j)
          CALL SetMatrixElement(Ai,j,j,s)
          ndir = ndir + 1
        END IF
      END DO
      PRINT *,'Initial Dirichlet cond:',NoRow,ndir,Ai % NumberOfRows
      
      BlockMatrix % MaxSize = MAX( BlockMatrix % MaxSize, ni )
    END DO
    
    CALL Info('DDInitVar','DD variables initialized!',Level=12)
      
  END SUBROUTINE DDInitVar


  !-------------------------------------------------------------------
  !> This subroutine copies back the full vector from its components.
  !------------------------------------------------------------------
  SUBROUTINE DDCopyDomains( Solver, BlockMatrix, x, FromPieces, ToPieces, Relax )
    
    IMPLICIT NONE
    
    TYPE(Solver_t), TARGET :: Solver
    TYPE(BlockMatrix_t), POINTER :: BlockMatrix
    REAL(KIND=dp) :: x(:)
    LOGICAL :: FromPieces, ToPieces
    REAL(KIND=dp), OPTIONAL :: Relax
    INTEGER :: NoRow
    
    
    TYPE(Matrix_t), POINTER :: Ai
    INTEGER :: i,j,k,n,ni,m,Novar
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp), POINTER :: xi(:)
    REAL(KIND=dp) :: c
    

    NoVar = BlockMatrix % NoVar

    IF( FromPieces ) THEN
      CALL Info('DDCopyDomains','Copying from pieces to global!',Level=12)
      IF(PRESENT(Relax)) THEN
        c = Relax
        x = (1-c) * x
      ELSE
        c = 1.0_dp
        x = 0.0_dp
      END IF
      
      DO NoRow=1,NoVar
        Var => TotMatrix % SubVector(NoRow) % Var
        Ai => TotMatrix % Submatrix(NoRow,NoRow) % Mat

        x(Ai % InvPerm) = x(Ai % InvPerm) + c * Var % Values / NoOwners(Ai % InvPerm) 
      END DO
    END IF
     
    IF( ToPieces ) THEN
      CALL Info('DDCopyDomains','Copying from global to pieces!',Level=12)
      DO NoRow=1,NoVar
        Ai => BlockMatrix % Submatrix(NoRow,NoRow) % Mat         
        Var => BlockMatrix % SubVector(NoRow) % Var 
        xi => Var % Values
        
        ! Copy the initial guess. 
        xi = x(Ai % InvPerm)      
      END DO
    END IF
      
  END SUBROUTINE DDCopyDomains

   

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

    s = SUM(x(1:n)**2)          
    nrm = SQRT( s / n )
    
  END FUNCTION CompNorm
  
  
  !------------------------------------------------------------------------------          
  !> Compute the rhs for the block matrix system which is solved
  !> accounting only the diagonal entries i.e. subtract the non-diagonal 
  !> matrix-vector results from the original r.h.s. vectors.
  !> After this the block diagonal problem Ax=b may be solved.
  !----------------------------------------------------------------------------------
  SUBROUTINE DDUpdateRhs( BlockMatrix, A, x, b, ResMode, ThisRow )
    
    TYPE(BlockMatrix_t), TARGET :: BlockMatrix
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: x(:), b(:)
    LOGICAL :: ResMode 
    INTEGER, OPTIONAL :: ThisRow    

    TYPE(Matrix_t), POINTER :: Ar,Ac
    INTEGER :: n, NoRow,NoCol, NoVar, nr, nc, i, j, k, ic, ir, ndir
    REAL(KIND=dp), POINTER :: rtmp(:),ones(:),rhs(:),xr(:),xc(:),TmpValues(:)
    REAL(KIND=dp) :: bnorm, s
    LOGICAL :: Found, FirstTime = .TRUE.
    TYPE(Variable_t), POINTER :: Var

    SAVE FirstTime
    
    CALL Info('DDUpdateRhs','Computing block r.h.s',Level=8)

    IF( FirstTime ) THEN
      IF( ListGetLogical(CurrentModel % Solver % Values,'DD Robin', Found ) ) THEN
        DDMode = 1
        RobinC = ListGetCReal(CurrentModel % Solver % Values,'DD Robin Coefficient')
        PRINT *,'RobinC:',RobinC
      END IF
    END IF
      
    NoVar = BlockMatrix % NoVar
    
    !A => CurrentModel % Solver % Matrix
    !x => CurrentModel % Solver % Variable % Values
    n = A % NumberOfRows
    
    ! This must be just big enough
    ALLOCATE( rtmp(n) ) 
    rtmp = 0.0_dp
    
    ! Initialilze r.h.s. from original system where "Ar % rhs" has been
    ! defined to be a scaling factor. 
    DO NoRow = 1,NoVar
      Ar => BlockMatrix % SubMatrix(NoRow,NoRow) % Mat
      rhs => BlockMatrix % SubVector(NoRow) % rhs
      rhs = b(Ar % InvPerm) * ( Ar % rhs )
    END DO

    ! DDMode 
    ! 1) Dirichlet-Neumann
    ! - some partitions sets the Dirichlet to others
    ! - the others return the Neumann
    !
    ! 2) Robin-Robin
    ! - all partitions set the Neumann for others
    ! - all partitions set the Robin with penalty
    
    
    ! For Neumann condition set the r.h.s. force from calculated load.
    DO NoRow = 1,NoVar       
      ! Optionally only one diagonal block may be updated for
      IF( PRESENT( ThisRow ) ) THEN
        IF( NoRow /= ThisRow ) CYCLE 
      END IF

      Ar => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat            
      Var => BlockMatrix % SubVector(NoRow) % Var
      xr => Var % Values
      nr = SIZE( xr )      

      !xr = x(Ar % InvPerm)      
      TmpValues => Ar % Values
      Ar % Values => Ar % BulkValues 
      CALL CRS_MatrixVectorMultiply( Ar, xr, rtmp)
      
      Ar % Values => TmpValues

      ! Update the nodal force i.e. Neumann.
      IF( DDMode == 0 ) THEN
        DO ir=1,nr
          i = Ar % InvPerm(ir)

          ! The owner piece gets all the loads from non-owner pieces.
          NoCol = OwnerPart(i)
          IF(NoCol /= NoRow) THEN
            Ac => BlockMatrix % SubMatrix( NoCol, NoCol ) % Mat
            rhs => BlockMatrix % SubVector(NoCol) % rhs            
            ic = Ac % Perm(i)
            rhs(ic) = rhs(ic) - rtmp(ir)
          END IF          
        END DO
      ELSE
        DO ir=1,nr
          i = Ar % InvPerm(ir)

          ! All pieces get the loads from other pieces. 
          DO NoCol=1,NoVar 
            IF(NoCol /= NoRow) THEN
              Ac => BlockMatrix % SubMatrix( NoCol, NoCol ) % Mat
              rhs => BlockMatrix % SubVector(NoCol) % rhs            
              ic = Ac % Perm(i)
              IF(ic>0) rhs(ic) = rhs(ic) - rtmp(ir)
            END IF
          END DO
        END DO
      END IF
    END DO

    
    IF( DDMode == 0 ) THEN
      DO NoRow = 1,NoVar 

        ! Optionally only one diagonal block may be updated for
        IF( PRESENT( ThisRow ) ) THEN
          IF( NoRow /= ThisRow ) CYCLE 
        END IF

        Ar => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat      
        nr = Ar % NumberOfRows

        Var => BlockMatrix % SubVector(NoRow) % Var
        xr => Var % Values

        rhs => BlockMatrix % SubVector(NoRow) % rhs

        ndir = 0
        DO ir=1,nr
          i = Ar % InvPerm(ir)
          k = OwnerPart(i)

          ! Dirichlet condition set from owner part to all non-owner parts.
          IF(k /= NoRow ) THEN
            s = 1.0_dp
            CALL ZeroRow(Ar,ir)
            CALL SetMatrixElement(Ar,ir,ir,s)
            rhs(ir) = s * x(i)
            ndir = ndir + 1
          END IF
        END DO
        !PRINT *,'Number of Dirichlet dofs:',NoRow,ndir
      END DO
    ELSE IF( DDMode == 1 ) THEN
      DO NoRow = 1,NoVar 
        Ar => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat      
        nr = Ar % NumberOfRows

        Var => BlockMatrix % SubVector(NoRow) % Var
        xr => Var % Values

        rhs => BlockMatrix % SubVector(NoRow) % rhs

        DO NoCol = 1, NoVar
          IF(NoCol == NoRow) CYCLE

          ! Set Robin BCs between all domains. 
          Ac => BlockMatrix % SubMatrix( NoRow, NoCol) % Mat        
          IF(.NOT. ASSOCIATED(Ac)) CYCLE
          IF(Ac % NumberOfRows == 0) CYCLE
          
          xc => BlockMatrix % SubVector(NoCol) % Var % Values
          rtmp = 0.0_dp
          ! Project the unknowns from other domains and penalize them. 
          CALL CRS_MatrixVectorMultiply( Ac, xc, rtmp)              
          rhs(1:nr) = rhs(1:nr) + robinC * rtmp(1:nr)

          IF(FirstTime) THEN
            nc = SIZE(xc)
            ALLOCATE(ones(nc))
            ones = 1.0_dp
            ! In the first time add diagonal entries to correspond for the implicit Robin part.
            CALL CRS_MatrixVectorMultiply( Ac, ones, rtmp)                         
            DEALLOCATE( ones )

            DO j=1,nr
              k = Ar % InvPerm(j)        
              ! Do not spoil the real Dirichlet BC's setting. 
              IF(A % ConstrainedDOF(k) ) CYCLE
              Ar % Values(Ar % Diag(j)) = Ar % Values(Ar % Diag(j)) + robinC * rtmp(j)
            END DO
          END IF
        END DO
      END DO
    END IF

    ! Set the Dirichlet conditions.
    DO NoRow = 1,NoVar
      Ar => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat
      rhs => BlockMatrix % SubVector(NoRow) % rhs    
      nr = Ar % NumberOfRows
      
      DO j=1,nr        
        k = Ar % InvPerm(j)        
        ! We know that the digonal has been scaled by the line below, so also scale the residual.
        ! s = Ai % Rhs(j) * A % Values(A % Diag(k)) 
        IF(A % ConstrainedDOF(k) ) THEN
          rhs(j) = Ar % Rhs(j) * b(k)  
        END IF
      END DO
    END DO

    
    IF( ResMode ) THEN
      DO NoRow=1, NoVar
        Ar => BlockMatrix % SubMatrix( NoRow, NoRow ) % Mat
        nr = Ar % NumberOfRows 

        rhs => BlockMatrix % SubVector(NoRow) % rhs

        bnorm = CompNorm(rhs,nr)
        BlockMatrix % SubVector(NoRow) % bnorm = bnorm

        ! Finally deduct the diagonal entry so that we can solve for the residual
        Var => BlockMatrix % SubVector(NoRow) % Var
        xr => Var % Values      

        CALL CRS_MatrixVectorMultiply( Ar, xr, rtmp)              
        rhs(1:nr) = rhs(1:nr) - rtmp(1:nr)
      END DO
    END IF
    DEALLOCATE( rtmp ) 
    
    FirstTime = .FALSE.
    
  END SUBROUTINE DDUpdateRhs

  
  
!------------------------------------------------------------------------------
!> Perform block preconditioning for Au=v by solving all the individual diagonal problems.
!> Has to be called outside the module by Krylov methods.
!------------------------------------------------------------------------------
  SUBROUTINE DDMatrixPrec( u,v,ipar )    
    IMPLICIT NONE
    REAL(KIND=dp), TARGET, INTENT(out) :: u(*)
    REAL(KIND=dp), TARGET, INTENT(in) :: v(*)
    INTEGER :: ipar(*)
!---------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: rtmp(:),vtmp(:),xtmp(:),btmp(:),diagtmp(:),b(:),xi(:)
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:)
    INTEGER :: i,j,k,l,n,m,NoVar,nc,kk,ll,istat
    TYPE(Solver_t), POINTER :: Solver, Solver_save
    INTEGER, POINTER :: Offset(:), ParPerm(:)
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER :: BlockOrder(:)
    TYPE(Matrix_t), POINTER :: A, Aij, mat_save
    TYPE(Variable_t), POINTER :: Var, Var_save
    REAL(KIND=dp) :: nrm
    LOGICAL :: GotOrder, BlockGS, Found, NS, ScaleSystem, DoSum, &
        IsComplex, UsePrecMat, CalcLoads, ResMode 
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER(KIND=AddrInt) :: AddrFunc
    EXTERNAL :: AddrFunc

    CALL Info('DDMatrixPrec','Starting block matrix preconditioning',Level=8)

    n = ipar(3)
    
    IF( InfoActive(25) ) THEN
      nrm = CompNorm(v(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'V start norm: ',nrm
      CALL Info('DDMatrixPrec',Message,Level=10)
    END IF
      
    Solver => CurrentModel % Solver
    Params => Solver % Values
    NoVar = TotMatrix % NoVar
    Solver => TotMatrix % Solver
    TotMatrix % NoIters = TotMatrix % NoIters + 1

    ! Enable user defined order for the solution of blocks
    !---------------------------------------------------------------
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotOrder)
    IF(GotOrder) THEN
      IF(SIZE(BlockOrder) /= NoVar) THEN
        CALL Fatal('DDMatrixPrec','Block Order size should be '//I2S(NoVar))
      END IF
    END IF

    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)    
         
    ! Save the initial solver stuff
    solver_save => Solver
    var_save => Solver % Variable
    mat_save => Solver % Matrix
    rhs_save => Solver % Matrix % RHS


    ! Always treat the inner iterations as truly complex if they are
    CALL ListAddLogical( Params,'Linear System Skip Complex',.FALSE.) 
    CALL ListAddLogical( Params,'Linear System Skip Loads',.TRUE.)

    ! Initial guess:
    !-----------------------------------------
    u(1:n) = 0.0_dp

    IF( BlockGS ) THEN
      ALLOCATE( vtmp(n), rtmp(n), xtmp(n),STAT=istat)
      IF ( istat /= 0 ) CALL Fatal('DDMatrixPrec','Memory allocation error for wrk space!')
      vtmp(1:n) = v(1:n)
    END IF
    
    CALL ListPushNameSpace('inner:')
    ResMode = .FALSE.
    CALL DDUpdateRhs(TotMatrix,SolverMatrix,u(1:n),v(1:n),ResMode)

    
    DO j=1,NoVar
      IF( GotOrder ) THEN
        i = BlockOrder(j)
      ELSE
        i = j
      END IF
      
      WRITE(Message,'(A,I0)') 'Solving block: ',i
      CALL Info('DDMatrixPrec',Message,Level=8)

      ! We do probably not want to compute the change within each iteration
      CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.TRUE.)         
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)         
      
      ! Set pointers to the new linear system
      !-------------------------------------------------------------------
      b => TotMatrix % SubVector(i) % rhs        
      Var => TotMatrix % SubVector(i) % Var
      Solver % Variable => Var      
      A => TotMatrix % Submatrix(i,i) % Mat
      xi => Var % Values
      
        
      IF( InfoActive(25) ) THEN
        nrm = CompNorm(b,offset(i+1)-offset(i))
        WRITE( Message,'(A,ES12.5)') 'Rhs '//I2S(i)//' norm: ',nrm
        CALL Info('DDMatrixPrec',Message,Level=10)
      END IF
                     
      btmp => b

      IF( A % COMPLEX ) THEN
        nc = 2
      ELSE
        nc = 1
      END IF

      IF( InfoActive( 15 ) ) THEN
        PRINT *,'rhs'//I2S(i)//':',SQRT( SUM(b**2) ), MINVAL( b ), MAXVAL( b ), SUM( b )
        PRINT *,'x'//I2S(i)//':',SQRT( SUM(xi**2) ), MINVAL( xi ), MAXVAL( xi ), SUM( xi )
        PRINT *,'Asize'//I2S(i)//':',A % NumberOfRows, SIZE(A % Cols), SIZE(A % Rows)
        PRINT *,'Avals'//I2S(i)//':',SQRT( SUM(A % Values**2) ), MINVAL( A % Values ), &
            MAXVAL( A % Values ), SUM( A % Values )
        PRINT *,'Sizes:',A % NumberOfRows, SIZE(b), SIZE(xi)
      END IF
      
      CALL SolveLinearSystem( A, b, xi, Var % Norm, Var % DOFs, Solver )
      
      IF( InfoActive(20) ) THEN
        nrm = CompNorm(xi,A % NumberOfRows,A=A) 
        WRITE( Message,'(A,ES12.5)') 'Linear system '//I2S(i)//' norm: ',nrm
        CALL Info('DDMatrixPrec',Message)
      END IF
      
    END DO ! j=1,NoVar

    ! From domains take average and create an average global solution. 
    CALL DDCopyDomains( Solver, TotMatrix, u(1:n), FromPieces=.TRUE., ToPieces=.TRUE.)
        
    CALL ListPopNameSpace('inner:') ! inner:

    CALL ListAddLogical( Params,'Linear System Refactorize',.FALSE. )
    CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.FALSE.)
    CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.FALSE.)
    CALL ListAddLogical( Params,'Linear System Skip Loads',.FALSE.)

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
      CALL Info('DDMatrixPrec',Message,Level=10)
      
      nrm = CompNorm(u(1:n),n)
      WRITE( Message,'(A,ES12.5)') 'U fin norm: ',nrm
      CALL Info('DDMatrixPrec',Message,Level=10)
    END IF
      
    CALL Info('DDMatrixPrec','Finished block matrix preconditioning',Level=8)
    
  END SUBROUTINE DDMatrixPrec



  !> This call takes care of Jacobi & Gauss Seidel block methods. 
  !-----------------------------------------------------------------
  SUBROUTINE BlockStandardIter( Solver, MaxChange )

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: MaxChange

    INTEGER :: i,j,NoVar,RowVar,iter,LinIter,MinIter
    INTEGER, POINTER :: BlockOrder(:)
    LOGICAL :: GotIt, GotBlockOrder, BlockGS
    REAL(KIND=dp), POINTER CONTIG :: rhs_save(:), b(:)
    REAL(KIND=dp), POINTER :: dx(:), x(:), xi(:)
    TYPE(Matrix_t), POINTER :: A, mat_save
    TYPE(Variable_t), POINTER :: Var, SolverVar
    REAL(KIND=dp) :: LinTol, TotNorm, dxnorm, xnorm, Relax
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ScaleSystem, Found, ResMode 

    NoVar = TotMatrix % NoVar
    Params => Solver % Values
    SolverVar => Solver % Variable
    x => SolverVar % Values

    BlockGS = ListGetLogical( Params,'Block Gauss-Seidel',Found)
    BlockOrder => ListGetIntegerArray( Params,'Block Order',GotBlockOrder)
    LinIter = ListGetInteger( Params,'Linear System Max Iterations',GotIt)
    MinIter = ListGetInteger( Params,'Linear System Min Iterations',GotIt)
    LinTol = ListGetConstReal( Params,'Linear System Convergence Tolerance',GotIt)

    ResMode = ListGetLogical( Params,'DD Residual Mode',Found )
    IF(.NOT. Found) ResMode = .TRUE.
    
    CALL ListPushNamespace('inner:')

    ! We don't want compute change externally
    CALL ListAddNewLogical( Params,'Skip compute nonlinear change',.TRUE.)

    Relax = ListGetConstReal( Params,'DD Relaxation Factor',GotIt)    
    IF(.NOT. GotIt) Relax = 1.0_dp

    DO iter = 1, LinIter
      ! Store the iteration count
      TotMatrix % NoIters = iter
      
      ! In block Jacobi the r.h.s. is not updated during the iteration cycle
      !----------------------------------------------------------------------
      IF( BlockGS ) THEN
        WRITE( Message,'(A,I0)') 'Block Gauss-Seidel iteration: ',iter
      ELSE
        WRITE( Message,'(A,I0)') 'Block Jacobi iteration: ',iter
      END IF

      A => Solver % Matrix
      CALL DDUpdateRhs(TotMatrix, A, x, A % rhs, ResMode)

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
        !IF( BlockGS ) THEN
        !  CALL DDUpdateRhs(TotMatrix,RowVar)
        !END IF
        
        !IF( ListGetLogical( Params,'Block Prec Reuse',GotIt) ) THEN
        !  DO j = 1, NoVar
        !    IF( j == RowVar ) CYCLE
        !    IF( CRS_CopyMatrixPrec( TotMatrix % Submatrix(j,j) % Mat, A ) ) EXIT
        !  END DO
        !END IF
        
        b => TotMatrix % SubVector(RowVar) % rhs
        

        Var => TotMatrix % SubVector(RowVar) % Var
        Solver % Variable => Var
        
        A => TotMatrix % Submatrix(RowVar,RowVar) % Mat

        ! Use the newly computed residual rather than original r.h.s. to solve the equation!!
        rhs_save => A % rhs ! Solver % Matrix % RHS
        A % RHS => b
        xi => Var % Values
        
        ! Solving the subsystem
        !-----------------------------------
        IF( InfoActive( 15 ) ) THEN
          PRINT *,'rhs'//I2S(i)//':',SQRT( SUM(b**2) ), MINVAL( b ), MAXVAL( b ), SUM( b )
          PRINT *,'x'//I2S(i)//':',SQRT( SUM(xi**2) ), MINVAL( xi ), MAXVAL( xi ), SUM( xi )          
          PRINT *,'Asize'//I2S(i)//':',A % NumberOfRows, SIZE(A % Cols), SIZE(A % Rows)
          PRINT *,'Avals'//I2S(i)//':',SQRT( SUM(A % Values**2) ), MINVAL( A % Values ), &
              MAXVAL( A % Values ), SUM( A % Values )
          PRINT *,'Sizes:',A % NumberOfRows, SIZE(b), SIZE(dx)
        END IF

        
        ALLOCATE( dx( SIZE( Var % Values ) ) )
        IF(ResMode) THEN
          dx = 0.0_dp
          CALL SolveLinearSystem( A, b, dx, Var % Norm, Var % DOFs, Solver )
          xi = xi + Relax * dx
        ELSE
          dx = xi
          CALL SolveLinearSystem( A, b, dx, Var % Norm, Var % DOFs, Solver )
          xi = (1.0_dp-Relax) * xi + Relax * dx          
        END IF
               
        ! Revert back to original r.h.s.
        A % RHS => rhs_save
        !Solver % Matrix => mat_save
        
        
#if 1
        dxnorm = SQRT(SUM(dx*dx))
        xnorm = SQRT(SUM(Var % Values**2))
#else   
        dxnorm = CompNorm(dx,A % NumberOfRows, A=A)
        xnorm = CompNorm(Var % Values, A % NumberOfRows, A=A)
#endif
        
        Var % Norm = xnorm
        Var % NonlinChange = dxnorm / xnorm

        WRITE(Message,'(A,2ES12.3)') 'Block '//I2S(RowVar)//' norms: ',xnorm, dxnorm / xnorm
        CALL Info('BlockStandardIter',Message,Level=12)
        
        IF( InfoActive( 20 ) ) THEN
          PRINT *,'dx'//I2S(i)//':',SQRT( SUM(dx**2) ), MINVAL( dx ), MAXVAL( dx ), SUM( dx ), SUM( ABS( dx ) )
        END IF
      
        DEALLOCATE( dx )
          
        TotNorm = TotNorm + Var % Norm
        MaxChange = MAX( MaxChange, Var % NonlinChange )        
      END DO

      ! From domains take average and create an average global solution. 
      CALL DDCopyDomains( Solver, TotMatrix, x, FromPieces=.TRUE., ToPieces=.TRUE.)
      
      WRITE(Message,'(A,2ES12.3)') 'Sum of norms: ',TotNorm, MaxChange
      CALL Info('BlockStandardIter',Message,Level=4)

      IF( MaxChange < LinTol .AND. iter >= MinIter ) THEN
        CALL Info('BlockStandardIter','Converged after iterations: '//I2S(iter),Level=5)
        EXIT
      END IF
      
    END DO
    CALL ListPopNamespace('inner:')

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
    INTEGER(KIND=AddrInt) :: iterProc,precProc, dotProc,nmrProc, zero=0
    REAL(KIND=dp) :: dpar(20), xnorm,prevxnorm
    REAL(KIND=dp), POINTER :: x(:),b(:),r(:)
    
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: SVar
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoVar, ndim, maxsize
    LOGICAL :: Converged, Diverged
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    INTEGER, POINTER :: BlockStruct(:),ParPerm(:)
    INTEGER :: i,j,k,l,ia,ib,istat
    LOGICAL :: LS, BlockAV,Found, UseMono
    CHARACTER(:), ALLOCATABLE :: VarName
    CHARACTER(*), PARAMETER :: Caller = 'BlockKrylovIter'
 
    
    CALL Info(Caller,'Starting block system iteration',Level=8)
    
    !CALL ListPushNameSpace('outer:')
    Params => Solver % Values
    
    ndim = TotMatrix % TotSize 
    NoVar = TotMatrix % NoVar

    TotMatrix % NoIters = 0

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

    ALLOCATE(r(ndim),STAT=istat)
    IF( istat /= 0 ) THEN
      CALL Fatal(Caller,'Cannot allocate temporal vectors of size: '//I2S(ndim))
    END IF
    
    x => SolverVar % Values 
    A => SolverMatrix
    b => A % RHS 

    r = 0.0_dp
          
    CALL Info(Caller,'Initializing monolithic system vectors',Level=18)

   
    
    !----------------------------------------------------------------------
    ! Solve matrix equation solver with the redefined block matrix operations
    !----------------------------------------------------------------------
    CALL ListAddLogical(Params,'Linear System Free Factorization',.FALSE.)

    precProc = AddrFunc(DDMatrixPrec)

    prevXnorm = CompNorm(b,ndim)
    WRITE( Message,'(A,ES12.5)') 'Rhs norm at start: ',PrevXnorm
    CALL Info(Caller,Message,Level=10)    
    IF( PrevXNorm < EPSILON( PrevXNorm ) ) THEN
      CALL Warn(Caller,"With zero'ish r.h.s. it does not make sense to call the solver!")
      RETURN
    END IF

    
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

    CALL IterSolver( A,x,b,Solver,ndim=ndim,PrecF=precProc) 
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
    
    CALL Info(Caller,'Finished block krylov iteration',Level=20)
   
  END SUBROUTINE blockKrylovIter

  
 
!------------------------------------------------------------------------------
!> An alternative handle for the block solvers to be used by the legacy matrix
!> type. 
!------------------------------------------------------------------------------
  SUBROUTINE DDSolveInt(A,x,b,Solver)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp), TARGET :: x(:)
    REAL(KIND=dp), TARGET CONTIG :: b(:)
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Variable_t), POINTER :: Var, DDVar
    INTEGER :: i,j,k,l,n,nr
    LOGICAL :: GotIt, BlockPrec, &
        BlockReIm, BlockDummy, BlockComplex, SkipPrec
    INTEGER :: ColVar, RowVar, NoVar, BlockDofs, VarDofs
    
    REAL(KIND=dp) :: TotNorm, MaxChange
    REAL(KIND=dp), POINTER :: SaveValues(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveRHS(:)
    LOGICAL :: Found, DoScaling 
    INTEGER, POINTER :: VarPerm(:)
    LOGICAL :: OldMatrix
    
    TYPE(Matrix_t), POINTER :: Amat, SaveCM
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params

    INTEGER, POINTER :: BlockIndex(:)
    CHARACTER(*), PARAMETER :: Caller = 'DDSolveInt'

    
    CALL Info(Caller,'---------------------------------------',Level=5)

    Params => Solver % Values
    Mesh => Solver % Mesh
    PSolver => Solver

    SolverRef => Solver
    
    OldMatrix = ASSOCIATED( Solver % BlockMatrix )
    IF( OldMatrix ) THEN
      CALL Info('BlockSolveInit','Using existing DD matrix structures!',Level=12)
      TotMatrix => Solver % BlockMatrix
    END IF
        
    ! Determine some parameters related to the block strategy
    !------------------------------------------------------------------------------
    BlockPrec = ListGetLogical( Params,'DD Preconditioner',GotIt)
    IF(BlockPrec) THEN
      CALL Info(Caller,'Using DD in preconditioning mode',Level=12)
    ELSE
      CALL Info(Caller,'Using DD in standard mode',Level=12)
    END IF

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
       
    CALL DDInitVar( Solver, TotMatrix, A, x )
    CALL DDCopyDomains( Solver, TotMatrix, x, FromPieces = .FALSE., ToPieces = .TRUE. )
    CALL Info(Caller,'DD matrix system created',Level=12)
    
    !------------------------------------------------------------------------------
    ! Finally solve the system using 'outer: ' as the optional namespace
    ! for the linear system setting.
    !------------------------------------------------------------------------------                
    TotNorm = 0.0_dp
    MaxChange = 0.0_dp

    DoScaling = ListGetLogical( Params,'Linear System Scaling',Found )
    IF( DoScaling ) THEN
      CALL Info(Caller,'Performing scaling on domains',Level=7)
      DO RowVar=1,NoVar      
        Amat => TotMatrix % SubMatrix(RowVar,RowVar) % Mat
        nr = Amat % NumberOfRows
        IF(.NOT. ASSOCIATED(Amat % DiagScaling) ) THEN
          ALLOCATE(Amat % DiagScaling(nr))
          Amat % DiagScaling = 0.0_dp
        END IF
        Amat % DiagScaling(1:nr) = A % DiagScaling(Amat % InvPerm)
      END DO
      CALL ScaleLinearSystem(Solver,Amat,b,x,Amat % DiagScaling )
      CALL ListAddLogical( Params,'Linear System Skip Scaling',.TRUE. )
    END IF
      
    ! Solve the DD system either as a solver strategy, or as preconditioner.
    !----------------------------------------------------------------------
    CALL ListPushNamespace('outer:')
    IF( BlockPrec ) THEN
      CALL Info(Caller,'Using block preconditioning strategy',Level=6)        
      CALL BlockKrylovIter( Solver, MaxChange )
    ELSE
      CALL Info(Caller,'Using block solution strategy',Level=6)
      CALL BlockStandardIter( Solver, MaxChange )
    END IF
    CALL ListPopNamespace('outer:')

    IF( DoScaling ) THEN
      CALL ListAddLogical( Params,'Linear System Skip Scaling',.FALSE. )
    END IF
          
    ! Create auxialiry variables for debugging mainly.
    DDVar => VariableGet( Solver % Mesh % Variables,'DD',ThisOnly=.TRUE.)
    IF(ASSOCIATED(DDVar)) THEN
      IF(DDVar % Dofs < NoVar ) THEN
        CALL Info(Caller,'Correct number of component in "DD"',Level=20)
      ELSE IF(DDVar % Dofs < NoVar ) THEN
        CALL Info(Caller,'Will only save '//I2S(DDVar % dofs)//' dofs for "DD"',Level=10)
        DDVar % Values = 0.0_dp
      ELSE IF(DDVar % Dofs > NoVar ) THEN
        CALL Fatal(Caller,'Too large size '//I2S(DDVar % dofs)//' for "DD"!')
      END IF
      DO RowVar=1,DDVar % Dofs
        Amat => TotMatrix % SubMatrix(RowVar,RowVar) % Mat
        Var => TotMatrix % SubVector(RowVar) % Var
        IF(.NOT. ASSOCIATED(Var)) CYCLE
        IF(.NOT. ASSOCIATED(Amat % InvPerm)) CYCLE
        DDVar % Values(DDVar % Dofs*(Amat % InvPerm-1)+RowVar) = Var % Values
      END DO
    END IF
    
    Var => VariableGet( Solver % Mesh % Variables,'NoOwners',ThisOnly=.TRUE.)
    IF(ASSOCIATED(Var)) THEN
      Var % Values = 1.0_dp * NoOwners
    END IF
    Var => VariableGet( Solver % Mesh % Variables,'OwnerPart',ThisOnly=.TRUE.)
    IF(ASSOCIATED(Var)) THEN
      Var % Values = 1.0_dp * OwnerPart
    END IF     
    
    ! For legacy matrices do the backmapping 
    !------------------------------------------
    SolverMatrix % RHS => SaveRHS
    Solver % Matrix => SaveMatrix
    Solver % Variable => SolverVar
    Solver % Variable % Values => SaveValues
       
    IF( ListGetLogical( Solver % Values,'Block Save Iterations',Found ) ) THEN
      CALL ListAddInteger(CurrentModel % Simulation,'res: block iterations',TotMatrix % NoIters)
    END IF   
    
    CALL Info(Caller,'All done')
    CALL Info(Caller,'-------------------------------------------------',Level=5)

  END SUBROUTINE DDSolveInt
  
END MODULE DDSolve




!------------------------------------------------------------------------------
!> Just a handle for SolveLinearSystem():
!------------------------------------------------------------------------------
SUBROUTINE DDSolveExt(A,x,b,Solver)
!------------------------------------------------------------------------------
    USE Types
    USE ParallelUtils
    USE DDSolve, ONLY: DDSolveInt 
    USE Lists, ONLY : ListGetLogical, ListAddLogical
    IMPLICIT NONE
    
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: x(:),b(:)
!------------------------------------------------------------------------------
    LOGICAL :: Found, bm
!------------------------------------------------------------------------------

    ! Eliminate recursion for block solvers. 
    bm = ListGetLogical(  Solver % Values, 'Linear System DD Mode', Found)
    IF(Found) &
      CALL ListAddLogical(Solver % Values,'Linear System DD Mode',.FALSE.)

    IF( ParEnv % PEs == 1 ) THEN
      CALL DDSolveInt(A,x,b,Solver)
    ELSE
      CALL Fatal('DDSolveExt','Implemented currently only in serial!')
    END IF
      
    IF(Found) &
      CALL ListAddLogical(Solver % Values,'Linear System DD Mode',bm )
!------------------------------------------------------------------------------
  END SUBROUTINE DDSolveExt
!------------------------------------------------------------------------------

!> \}
