!/******************************************************************************
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
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Linear system assembly, solve, and control utilities.
!>  Extracted from SolverUtils.
!------------------------------------------------------------------------------

#include "../config.h"

MODULE SolveCore

    USE SolverBasics
    USE BoundaryConditionUtils
    USE IterSolve, ONLY : NumericalError
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE NonlinearAcceleration(A,x,b,Solver,PreSolve,NoSolve)    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) CONTIG :: b(:),x(:)
    TYPE(Solver_t) :: Solver
    LOGICAL :: PreSolve
    LOGICAL, OPTIONAL :: NoSolve
    !------------------------------------------------------------------------------
    ! We have a special structure for the iterates and residuals so that we can
    ! cycle over the pointers instead of the values. 
    TYPE AndersonVect_t
      LOGICAL :: Additive
      REAL(KIND=dp), POINTER :: Iterate(:), Residual(:), Ax(:)
      INTEGER :: tag
    END TYPE AndersonVect_t
    TYPE(AndersonVect_t), ALLOCATABLE :: AndersonBasis(:), AndersonTmp
    INTEGER :: AndersonInterval, ItersCnt, AndersonVecs, VecsCnt, iter, n,i,j,k
    TYPE(Variable_t), POINTER :: iterV, Svar
    REAL(KIND=dp), ALLOCATABLE :: Alphas(:),AxTable(:,:),TmpVec(:) 
    REAL(KIND=dp) :: Nrm, AndersonRelax
    LOGICAL :: Found, DoRelax, KeepBasis, Visited = .FALSE., Parallel    
    INTEGER :: LinInterval
    INTEGER :: PrevSolverId = -1
    
    SAVE AndersonBasis, TmpVec, Alphas, ItersCnt, AndersonInterval, VecsCnt, AndersonVecs, &
        PrevSolverId, AxTable, AndersonRelax, DoRelax, Visited, KeepBasis, LinInterval
        
    IF( PreSolve ) THEN
      CALL Info('NonlinearAcceleration','Performing pre-solution steps',Level=8)
    ELSE
      CALL Info('NonlinearAcceleration','Performing post-solution steps',Level=8)
    END IF

    Parallel = ( ParEnv % PEs > 1 ) 
        
    iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter',UnfoundFatal=.TRUE.)
    iter = NINT(iterV % Values(1))

    IF(PRESENT(NoSolve)) NoSolve = .FALSE.
    
    n = A % NumberOfRows
          
    IF(.NOT. Visited ) THEN
      PrevSolverId = Solver % SolverId
      CALL Info('NonlinearAcceleration','Allocating structures for solution history',Level=8)

      AndersonInterval = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration Interval',Found)      
      LinInterval = ListGetInteger( Solver % Values,&
          'Linear System Acceleration Interval',Found)      

      AndersonVecs = MAX( AndersonInterval, LinInterval )
      IF( AndersonVecs == 0 ) THEN
        CALL Fatal('NonlinearAcceleration','Both acceleration intervals are zero!')
      END IF
            
      AndersonRelax = ListGetCReal( Solver % Values,&
          'Nonlinear System Acceleration Relaxation',DoRelax)
      KeepBasis = ListGetLogical( Solver % Values,&
          'Nonlinear System Acceleration Keep Vectors',Found)            

      ItersCnt = 0    ! relates to "AndersonInterval"
      VecsCnt = 0     ! relates to "AndersonVecs"
      
      IF(.NOT. ALLOCATED( AndersonBasis ) ) THEN
        ALLOCATE( AndersonBasis(AndersonVecs) )
        DO i=1,AndersonVecs
          ALLOCATE( AndersonBasis(i) % Residual(n), &
              AndersonBasis(i) % Iterate(n) )
          AndersonBasis(i) % Residual = 0.0_dp
          AndersonBasis(i) % Iterate = 0.0_dp          
        END DO
        ALLOCATE( TmpVec(n), Alphas(AndersonVecs) )
      END IF
      Visited = .TRUE.
    END IF
    
    IF( PrevSolverId /= Solver % SolverId ) THEN
      CALL Fatal('NonlinearAcceleration','Current implementation only supports one solver!')
    END IF
      
    
    IF( PreSolve ) THEN           
      IF( iter == 1 ) THEN
        ItersCnt = 0
        IF( .NOT. KeepBasis ) VecsCnt = 0
      END IF

      ItersCnt = ItersCnt + 1
      VecsCnt = VecsCnt + 1

      ! Calculate the residual of the matrix equation
      ! Here 'x' comes before being modified hence A(x) is consistent. 
      CALL MatrixVectorMultiply( A, x, TmpVec )
      TmpVec = TmpVec - b

      ! Add the iterate and residual to the basis vectors.
      ! This is fast as we operate with pointers mainly.
      AndersonTmp = AndersonBasis(AndersonVecs)        
      DO i=AndersonVecs,2,-1
        AndersonBasis(i) = AndersonBasis(i-1)
      END DO
      AndersonBasis(1) = AndersonTmp
      AndersonBasis(1) % Residual = TmpVec
      AndersonBasis(1) % Iterate = x 

      ! Pure Anderson sweep is done every AndersonInterval iterations if we have full basis.
      IF(.NOT. DoRelax .AND. AndersonInterval > 0 ) THEN
        IF( VecsCnt >= AndersonVecs .AND. ItersCnt >= AndersonInterval ) THEN
          CALL AndersonMinimize( )
          ItersCnt = 0
          IF(PRESENT(NoSolve)) NoSolve = .TRUE.
          RETURN
        END IF
      END IF
      
      IF( LinInterval > 0 ) THEN
        CALL AndersonGuess()
      END IF
    ELSE
      ! Relaxation strategy is done after each linear solve.
      IF( DoRelax ) THEN
        CALL Info('NonlinearAcceleration','Minimizing residual using history data',Level=6)
        CALL AndersonMinimize( )
      END IF
    END IF

  CONTAINS 


    !------------------------------------------------------------------------------
    FUNCTION Mydot( n, x, y ) RESULT(s)
      !------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp)  :: s
      REAL(KIND=dp) CONTIG :: x(:)
      REAL(KIND=dp) CONTIG, OPTIONAL :: y(:)
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        IF( PRESENT( y ) ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = DOT_PRODUCT( x(1:n), x(1:n) )
        END IF
      ELSE
        IF( PRESENT( y ) ) THEN
          s = ParallelDot( n, x, y )
        ELSE
          s = ParallelDot( n, x, x )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END FUNCTION Mydot
    !------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    SUBROUTINE Mymv( A, x, b, Update )
      !------------------------------------------------------------------------------
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: A
      LOGICAL, OPTIONAL :: Update
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        CALL CRS_MatrixVectorMultiply( A, x, b )
      ELSE
        IF ( PRESENT( Update ) ) THEN
          CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
        ELSE
          CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END SUBROUTINE Mymv
    !------------------------------------------------------------------------------

    
    ! Given set of basis vectors and residuals find a new suggestion for solution.
    ! Either use as such or combine it to solution when relaxation is used.
    ! This is applied to boost nonlinear iteration. 
    !------------------------------------------------------------------------------
    SUBROUTINE AndersonMinimize()
      INTEGER ::m, n, AndersonMinn
      REAL(KIND=dp) :: rr, rb
      
      m = MIN( ItersCnt, AndersonInterval )      
      
      AndersonMinN = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration First Iteration',Found )
      IF(.NOT. (Found .OR. DoRelax)) AndersonMinN = AndersonInterval
            
      ! Nothing to do 
      IF( m < AndersonMinN ) RETURN
      
      ! If size of our basis is just one, there is not much to do...
      ! We can only perform classical relaxation. 
      IF( m == 1 ) THEN
        x = AndersonRelax * x + (1-AndersonRelax) * AndersonBasis(1) % Iterate
        RETURN
      END IF
      
      ! If we are converged then the solution should already be the 1st component.
      ! Hence use that as the basis. 
      Alphas(1) = 1.0_dp     
      TmpVec = AndersonBasis(1) % Residual
      
      ! Minimize the residual
      n = SIZE( AndersonBasis(1) % Residual ) 
      DO k=2,m
        rr = MyDot( n, AndersonBasis(k) % Residual ) 
        rb = MyDot( n, AndersonBasis(k) % Residual, TmpVec )         
        Alphas(k) = -rb / rr 
        TmpVec = TmpVec + Alphas(k) * AndersonBasis(k) % Residual
      END DO

      ! Normalize the coefficients such that the sum equals unity
      ! This way for example, Dirichlet BCs will be honored.
      Alphas = Alphas / SUM( Alphas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Alpha(',i,') = ',Alphas(i)
          CALL Info('NonlinearAcceleration',Message)
        END DO
      END IF
              
      ! Create the new suggestion for the solution vector
      ! We take part of the suggested new solution vector 'x' and
      ! part of minimized residual that was used in anderson acceleration.
      IF( DoRelax ) THEN
        Alphas = Alphas * (1-AndersonRelax)
        x = AndersonRelax * x
        DO k=1,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      ELSE
        x = Alphas(1) * AndersonBasis(1) % Iterate
        DO k=2,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      END IF
        
    END SUBROUTINE AndersonMinimize

    
    ! Given set of basis vectors and a linear system
    ! find a combincation of the vectors that minimizes the norm of the linear
    ! system. This may be used to provide a better initial guess for a linear system.
    !--------------------------------------------------------------------------------
    SUBROUTINE AndersonGuess()
      INTEGER :: AndersonMinN

      REAL(KIND=dp), POINTER, SAVE ::Betas(:), Ymat(:,:)
      LOGICAL, SAVE :: AllocationsDone = .FALSE.
      INTEGER :: i,j,m
      
      IF(.NOT. AllocationsDone ) THEN
        m = LinInterval
        DO i=1,LinInterval
          ALLOCATE( AndersonBasis(i) % Ax(n) )
          AndersonBasis(i) % Ax = 0.0_dp
        END DO
        ALLOCATE(Betas(m),Ymat(m,m))
        AllocationsDone = .TRUE.
      END IF
      
      m = MIN( VecsCnt, LinInterval )      

      ! Calculate the residual of the matrix equation
      DO i=1,m
        CALL Mymv( A, AndersonBasis(i) % Iterate, AndersonBasis(i) % Ax )
      END DO

      DO i=1,m
        DO j=i,m
          Ymat(i,j) = SUM( AxTable(:,i) * AxTable(:,j) )
          Ymat(j,i) = Ymat(i,j)
        END DO
        Betas(i) = SUM( AxTable(:,i) * b )
      END DO
      
      CALL LUSolve(m, YMat(1:m,1:m), Betas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Beta(',i,') = ',Betas(i)
          CALL Info('NonLinearAcceleration',Message)
        END DO
      END IF
                                
      x = Betas(m) * AndersonBasis(m) % Iterate
      DO i=1,m-1
        x = x + Betas(i) * AndersonBasis(i) % Iterate
      END DO

    END SUBROUTINE AndersonGuess
    
  END SUBROUTINE NonlinearAcceleration
!------------------------------------------------------------------------------

  

  

  SUBROUTINE CalculateLoads( Solver, Aaid, x, DOFs, UseBulkValues, NodalLoads, NodalValues ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER  :: Aaid
    REAL(KIND=dp) CONTIG :: x(:)
    INTEGER :: DOFs
    LOGICAL :: UseBulkValues
    TYPE(Variable_t), POINTER, OPTIONAL :: NodalLoads
    REAL(KIND=dp), POINTER, OPTIONAL :: NodalValues(:)
    
    INTEGER :: i,j,k,l,m,ii,This,DOF
    REAL(KIND=dp), POINTER :: TempRHS(:), TempVector(:), Rhs(:), TempX(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
    REAL(KIND=dp) :: Energy, Energy_im
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Found, Rotated
    REAL(KIND=dp), ALLOCATABLE :: BoundarySum(:), BufReal(:)
    INTEGER, ALLOCATABLE :: BoundaryShared(:),BoundaryActive(:),DofSummed(:),BufInteg(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: bc, ind, NoBoundaryActive, NoBCs, ierr
    LOGICAL :: OnlyGivenBCs, IgnorePeriodic
    LOGICAL :: UseVar, Parallel
    LOGICAL, SAVE :: SharedWarned = .FALSE.


    Parallel = Solver % Parallel
      
    UseVar = .FALSE.
    IF(PRESENT( NodalLoads ) ) THEN
      UseVar = ASSOCIATED( NodalLoads )
      IF(.NOT. UseVar ) THEN
        CALL Warn('CalculateLoads','Load variable not associated!')
        RETURN
      END IF
    ELSE IF( PRESENT( NodalValues ) ) THEN
      IF(.NOT. ASSOCIATED( NodalValues ) ) THEN
        CALL Warn('CalculateLoads','Load values not associated!')
        RETURN
      END IF
    ELSE
      CALL Fatal('CalculateLoads','Give either loads variable or values as parameter!')
    END IF
    
    ALLOCATE( TempVector(Aaid % NumberOfRows) )

    IF( UseBulkValues ) THEN
      IF(.NOT. ASSOCIATED(Aaid % BulkValues)) THEN
        CALL Fatal('CalculateLoads','"BulkValues" are not associated!')
      END IF
      SaveValues => Aaid % Values      
      Aaid % Values => Aaid % BulkValues
      Rhs => Aaid % BulkRHS
    ELSE
      Rhs => Aaid % Rhs
    END IF
    
    IF ( Parallel ) THEN
      IF( ASSOCIATED( Rhs ) ) THEN
        ALLOCATE(TempRHS(SIZE(Rhs)))
        TempRHS = Rhs 
        CALL ParallelInitSolve( Aaid, x, TempRHS, Tempvector )
      END IF
      CALL ParallelMatrixVector( Aaid, x, TempVector, .TRUE. )
    ELSE
      CALL MatrixVectorMultiply( Aaid, x, TempVector )
    END IF

    IF(ListGetLogical(Solver % Values, 'Calculate Energy Inner Product', Found) .OR. &
        ListGetLogical(Solver % Values, 'Calculate Energy Norm', Found) ) THEN
      Energy = 0._dp
      IF( ListGetLogical(Solver % Values, 'Linear System Complex', Found) ) THEN
        Energy_im = 0._dp
        DO i = 1, (Aaid % NumberOfRows / 2)
          IF ( Parallel ) THEN
            IF ( Aaid% ParMatrix % ParallelInfo % &
                NeighbourList(2*(i-1)+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          END IF
          Energy    = Energy    + x(2*(i-1)+1) * TempVector(2*(i-1)+1) + x(2*(i-1)+2) * TempVector(2*(i-1)+2)
          Energy_im = Energy_im + x(2*(i-1)+1) * TempVector(2*(i-1)+2) - x(2*(i-1)+2) * TempVector(2*(i-1)+1) 
        END DO
        Energy    = ParallelReduction(Energy)
        Energy_im = ParallelReduction(Energy_im)

        CALL ListAddConstReal( Solver % Values, 'Energy inner product', Energy)
        CALL ListAddConstReal( Solver % Values, 'Energy inner product im', Energy_im)

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product im'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy_im )

        WRITE( Message, * ) 'Energy inner product: ', Energy, Energy_im
        CALL Info( 'CalculateLoads', Message, Level=5)
      ELSE 
        DO i=1,Aaid % NumberOfRows
          IF ( Parallel ) THEN
            IF ( Aaid % ParMatrix % ParallelInfo % &
                NeighbourList(i) % Neighbours(1) /= Parenv % MyPE ) CYCLE
          END IF
          Energy = Energy + x(i)*TempVector(i)
        END DO
        Energy = ParallelReduction(Energy)
        CALL ListAddConstReal( Solver % Values, 'Energy inner product', Energy )

        WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy inner product'
        CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

        WRITE( Message, * ) 'Energy inner product: ', Energy
        CALL Info( 'CalculateLoads', Message, Level=5)
      END IF
    END IF

    IF( ASSOCIATED( Rhs ) ) THEN
      IF ( Parallel ) THEN
        DO i=1,Aaid % NumberOfRows
          IF ( AAid % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % Mype ) THEN
            TempVector(i) = TempVector(i) - TempRHS(i)
          ELSE
            TempVector(i) = 0
          END IF
        END DO
        CALL ParallelSumVector( AAid, Tempvector )
        DEALLOCATE( TempRhs ) 
      ELSE
        TempVector = TempVector - RHS
      END IF
    END IF
          
    IgnorePeriodic = ListGetLogical( Solver % Values,'Calculate Loads Ignore Periodic',Found )
    
    NoBCs = CurrentModel % NumberOfBCs

    IF( IgnorePeriodic ) THEN
      DO This=1,NoBCs
        Projector => CurrentModel  % BCs(This) % PMatrix
        IF (ASSOCIATED(Projector)) THEN
          DO DOF=1,DOFs
            DO i=1,Projector % NumberOfRows
              ii = Projector % InvPerm(i)
              IF( ii == 0 ) CYCLE
              k = Solver % Variable % Perm(ii)
              IF(k<=0) CYCLE
              k = DOFs * (k-1) + DOF
              TempVector(k)=0

              DO l = Projector % Rows(i), Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = Solver % Variable % Perm( Projector % Cols(l) )
                IF ( m > 0 ) THEN
                  m = DOFs * (m-1) + DOF
                  TempVector(k) = TempVector(k) + Projector % Values(l)*TempVector(m)
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END IF
      
    IF( UseVar ) THEN
      DO i=1,SIZE( NodalLoads % Perm )
        j = NodalLoads % Perm(i)
        k = Solver % Variable % Perm(i)
        IF ( j>0 .AND. k>0 ) THEN
          DO dof=1,DOFs
            NodalLoads % Values(DOFs*(j-1)+dof) = TempVector(DOFs*(k-1)+dof)
          END DO
        END IF
      END DO
    ELSE
      NodalValues = TempVector
    END IF
      
    IF( ListGetLogical( Solver % Values,'Calculate Boundary Fluxes',Found ) ) THEN
      CALL Info('CalculateLoads','Computing boundary fluxes from nodal loads',Level=6)

      IF( Solver % Mesh % MaxEdgeDofs > 1 .OR. Solver % Mesh % MaxFaceDOFs > 1 ) THEN
        CALL Warn('CalculateLoads','Boundary flux computation implemented only for nodes for now!')
      END IF

      IF(.NOT. UseVar ) THEN
        CALL Fatal('CalculateLoads','Boundary flux computation needs the variable parameter!')        
      END IF
      
      ALLOCATE( BoundarySum( NoBCs * DOFs ), &
          BoundaryActive( NoBCs ), &
          BoundaryShared( NoBCs ), &
          DofSummed( MAXVAL( NodalLoads % Perm ) ) )
      BoundarySum = 0.0_dp
      BoundaryActive = 0
      BoundaryShared = 0
      DofSummed = 0

      OnlyGivenBCs = ListCheckPresentAnyBC( CurrentModel,'Calculate Boundary Flux')
      
      k = Solver % Mesh % NumberOfBulkElements
      DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
        Element => Solver % Mesh % Elements(i)
        bc = Element % BoundaryInfo % Constraint
           
        IF( bc == 0 ) CYCLE

        IF( OnlyGivenBCs ) THEN
          IF (.NOT. ListGetLogical( CurrentModel % BCs(bc) % Values,&
              'Calculate Boundary Flux',Found) ) CYCLE
        END IF

        DO j=1,Element % TYPE % NumberOfNodes
          ind = NodalLoads % Perm( Element % NodeIndexes(j) )
          IF( ind == 0 ) CYCLE

          ! In this partition sum up only the true owners
          IF ( Parallel ) THEN
            IF ( AAid % ParallelInfo % NeighbourList(ind) % Neighbours(1) &
                /= ParEnv % Mype ) CYCLE
          END IF

          ! Only sum each entry once. If there is a conflict we cannot 
          ! really resolve it with the chosen method so just warn. 
          IF( DofSummed(ind) == 0 ) THEN
            BoundarySum( DOFs*(bc-1)+1 :DOFs*bc ) = BoundarySum( DOFs*(bc-1)+ 1:DOFs*bc ) + &
                NodalLoads % Values( DOFs*(ind-1) + 1: DOFs * ind )
            DofSummed( ind ) = bc
            BoundaryActive( bc ) = 1
          ELSE IF( bc /= DofSummed(ind) ) THEN
            BoundaryShared(bc) = 1
            BoundaryShared(DofSummed(ind)) = 1
          END IF
        END DO
      END DO
      
      NoBoundaryActive = 0
      IF( Parallel ) THEN
        ALLOCATE( BufInteg( NoBCs ), BufReal( NoBCs * DOFs ) )

        BufInteg = BoundaryActive
        CALL MPI_ALLREDUCE( BufInteg, BoundaryActive, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        BufInteg = BoundaryShared
        CALL MPI_ALLREDUCE( BufInteg, BoundaryShared, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        BufReal = BoundarySum 
        CALL MPI_ALLREDUCE( BufReal, BoundarySum, DOFs * NoBCs, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        DEALLOCATE( BufInteg, BufReal ) 
      END IF

      DO i=1,CurrentModel % NumberOfBCs 
        IF( BoundaryActive(i) == 0 ) CYCLE
        IF( BoundaryShared(i) > 0) THEN
          ! This is a property of the mesh and the BCs, so it does not change
          ! between the calls -- and there is one call per nonlinear iteration
          ! and timestep. Say it once, then keep it for the curious only.
          IF( SharedWarned ) THEN
            CALL Info('CalculateLoads','Boundary '//I2S(i)//&
                ' includes inseparable dofs!',Level=8)
          ELSE
            CALL Warn('CalculateLoads','Boundary '//I2S(i)//' includes inseparable dofs!')
          END IF
        END IF
        NoBoundaryActive = NoBoundaryActive + 1

        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over BC '//I2S(i)
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//I2S(j)//' Flux over BC '//I2S(i)
          END IF
          CALL ListAddConstReal( CurrentModel % Simulation, 'res: '//TRIM(Message), &
              BoundarySum(DOFs*(i-1)+j) )
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',BoundarySum(DOFs*(i-1)+j)
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END DO
      SharedWarned = .TRUE.
      
      IF( NoBoundaryActive > 1 ) THEN
        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over all BCs'
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//I2S(j)//' Flux over all BCs'
          END IF
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',SUM(BoundarySum(j::DOFs))
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END IF
      
      DEALLOCATE( DofSummed, BoundaryShared, BoundaryActive, BoundarySum )      


      IF( ListGetLogical( Solver % Values,'Calculate Boundary Weights', Found ) ) THEN
        BLOCK
          CHARACTER(MAX_NAME_LEN) :: Name   
          TYPE(Variable_t), POINTER :: WVar, FVar
          REAL(KIND=dp) :: eps

          CALL Info('CalculateLoads','Computing fluxes on boundaries!',Level=10)
          
          Name = GetVarName(Solver % Variable) // ' Boundary Weights'
          WVar => VariableGet( Solver % Mesh % Variables, Name )

          IF(.NOT. ASSOCIATED(WVar) ) THEN
            CALL Fatal('CalculateLoads','Weight variable is not available!')
          END IF

          Name = GetVarName(Solver % Variable) // ' Boundary Flux'        
          FVar => VariableGet( Solver % Mesh % Variables, Name ) 
          IF(.NOT. ASSOCIATED(FVar) ) THEN
            CALL VariableAddVector( Solver % Mesh % Variables,&
                Solver % Mesh, Solver, Name, Solver % Variable % DOFs, Secondary = .TRUE., &
                Perm = WVar % Perm )
            FVar => VariableGet( Solver % Mesh % Variables, Name ) 
          END IF
          IF(.NOT. ASSOCIATED(FVar) ) THEN
            CALL Fatal('CalculateLoads','Flux variable is not available!')
          END IF

          eps = EPSILON(eps)
          FVar % Values = 0.0_dp
          
          DO i=1,SIZE( WVar % Perm )
            j = WVar % Perm(i)
            IF(j==0) CYCLE
            k = Solver % Variable % Perm(i)
            IF(k==0) CYCLE

            DO dof=1,DOFs
              IF( WVar % Values(j) < eps ) CYCLE 
              FVar % Values(DOFs*(j-1)+dof) = TempVector(DOFs*(k-1)+dof) / WVar % Values(j)
            END DO
          END DO

        END BLOCK
      END IF
    END IF
    
    DEALLOCATE( TempVector )
    
    IF( UseBulkValues ) THEN
      Aaid % Values => SaveValues
    END IF

  END SUBROUTINE CalculateLoads





  ! Create a boundary matrix and at calculate step compute the boundary loads
  ! for one given body. This is not called by default but the user needs to
  ! include it in the code, both at assembly and after solution.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsAssembly( Solver, Element, LocalMatrix, LocalRhs )

    TYPE(Solver_t) :: Solver
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: LocalMatrix(:,:)
    REAL(KIND=dp) :: LocalRhs(:)

    LOGICAL :: FirstStep
    INTEGER :: i,j,k,l,n,Row,Col,Dofs,ElemNo,TargetBody=-1
    TYPE(Matrix_t), POINTER :: BCMat
    REAL(KIND=dp) :: Val
    LOGICAL :: Found
    INTEGER, POINTER :: Perm(:), BCPerm(:)
    CHARACTER(:), ALLOCATABLE :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    SAVE :: BCMat, TargetBody, BCPerm, Perm, Dofs


    FirstStep = ( Solver % ActiveElements(1) == Element % ElementIndex )

    IF( FirstStep ) THEN
      CALL Info('BCLoadsAssembly','Visiting first element',Level=6)
 
      BCMat => Solver % Matrix % EMatrix
      IF(.NOT. ASSOCIATED( BCMat ) ) THEN
        TargetBody = ListGetInteger( Solver % Values,'Boundary Loads Target Body',Found )
        IF( Found ) THEN
          CALL Info('BCLoadsAssembly','Target body set to: '//I2S(TargetBody),Level=6)       
        ELSE
          TargetBody = -1
          RETURN
        END IF

        CALL Info('BCLoadsAssembly','Allocating structures for load computation',Level=8)
        IF ( ParEnv % PEs > 1 ) THEN
          CALL Warn('BCLoadsAssembly','Not implemented in parallel')
        END IF

        ! Mark the active nodes
        ALLOCATE( BCPerm( Solver % Mesh % NumberOfNodes ) )
        BCPerm = 0

        ElemNo = 0
        k = Solver % Mesh % NumberOfBulkElements
        DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
          Element => Solver % Mesh % Elements(i)
          Found = .FALSE.             
          IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
            Found = ( Element % BoundaryInfo % Left % BodyId == TargetBody )
          END IF
          IF(.NOT. Found ) THEN
            IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              Found = ( Element % BoundaryInfo % Right % BodyId == TargetBody )
            END IF
          END IF
          IF( Found ) THEN
            ElemNo = ElemNo + 1
            BCPerm( Element % NodeIndexes ) = 1
          END IF
        END DO

        CALL Info('BCLoadsAssembly','Number of related boundary elements: '//I2S(ElemNo),Level=8)

        n = 0
        DO i=1,Solver % Mesh % NumberOfNodes
          IF( BCPerm(i) > 0 ) THEN
            n = n + 1
            BCPerm(i) = n
          END IF
        END DO
        CALL Info('BCLoadsAssembly','Number of active nodes: '//I2S(n),Level=8)

        ! Create the list matrix 
        BCMat => AllocateMatrix()
        BCMat % Format = MATRIX_LIST           
        CALL AddToMatrixElement( BCMat, n, n, 0.0_dp )
        Solver % Matrix % EMatrix => BCMat

        ALLOCATE( BCMat % Rhs(n) )
        BCMat % Rhs = 0.0_dp
      END IF

      ! When visiting the routine after the 1st iteration the matrix for is already CRS 
      IF( BCMat % Format == MATRIX_CRS ) THEN
        BCMat % Values = 0.0_dp
        BCMat % Rhs = 0.0_dp
      END IF

      Dofs = Solver % Variable % Dofs
      Perm => Solver % Variable % Perm

      Name = TRIM(Solver % Variable % Name)//' BCLoads'
      BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
      IF(.NOT. ASSOCIATED( BCVar ) ) THEN
        CALL Info('CalculateBCLoads','Creating variable: '//TRIM(Name),Level=7)
        CALL VariableAddVector( Solver % Mesh % Variables,&
            Solver % Mesh, Solver, Name, DOFs, Perm = BCPerm )
      END IF
      
    END IF

    IF( Element % BodyId == TargetBody ) THEN       
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        IF ( BCPerm( Element % NodeIndexes(i) ) == 0 ) CYCLE
        DO k=0,Dofs-1
          Row = Dofs * BCPerm( Element % NodeIndexes(i) ) - k
          BCMat % rhs(Row) = BCMat % rhs(Row) + LocalRhs(Dofs*i-k)
          DO j=1,n
            DO l=0,Dofs-1
              Col = Dofs * Perm( Element % NodeIndexes(j) ) - l
              Val = LocalMatrix(Dofs*i-k,Dofs*j-l)
              CALL AddToMatrixElement(BCMat,Row,Col,Val)
            END DO
          END DO
        END DO
      END DO
    END IF


  END SUBROUTINE BCLoadsAssembly


  ! Calculate the boundary loads resulting from the action of boundary matrix.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsComputation( Solver )

    TYPE(Solver_t) :: Solver

    TYPE(Matrix_t), POINTER :: BCMat
    CHARACTER(:), ALLOCATABLE :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    BCMat => Solver % Matrix % EMatrix
    IF(.NOT. ASSOCIATED( BCMat ) ) THEN
      CALL Fatal('BCLoadsComputation','We should have the boundary matrix!')
    END IF
        
    CALL Info('BCLoadsComputation','Computing boundary loads',Level=6)
    IF( BCMat % FORMAT == MATRIX_LIST ) THEN
      CALL List_ToCRSMatrix( BCMat )
      CALL Info('BCLoadsComputation','Matrix format changed to CRS',Level=8)
    END IF

    Name = TRIM(Solver % Variable % Name)//' BCLoads'
    BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
    IF(.NOT. ASSOCIATED( BCVar ) ) THEN
      CALL Fatal('BCLoadsComputation','Variable not present: '//TRIM(Name))
    END IF
    
    CALL MatrixVectorMultiply( BCMat, Solver % Variable % Values, BCVar % Values )
    BCVar % Values = BCVar % Values - BCMat % rhs

    CALL Info('BCLoadsComputation','All done',Level=12)

  END SUBROUTINE BCLoadsComputation


    
!------------------------------------------------------------------------------
!> Prints the values of the CRS matrix to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintMatrix( A, Parallel, CNumbering,SaveMass, SaveDamp, SaveStiff, SkipZeros )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A            !< Structure holding matrix
    LOGICAL :: Parallel    !< are we in parallel mode?
    LOGICAL :: CNumbering  !< Continuous numbering ?
    LOGICAL, OPTIONAL :: SaveMass  !< Should we save the mass matrix
    LOGICAL, OPTIONAL :: SaveDamp  !< Should we save the damping matrix
    LOGICAL, OPTIONAL :: SaveStiff !< Should we save the stiffness matrix
    LOGICAL, OPTIONAL :: SkipZeros !< Should we write zeros or not
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,IndMass,IndDamp,IndStiff,IndMax,row,col
    LOGICAL :: DoMass, DoDamp, DoStiff, Found, Skip0
    REAL(KIND=dp) :: Vals(3), val
    INTEGER, ALLOCATABLE :: Owner(:)

    DoMass = .FALSE.
    IF( PRESENT( SaveMass ) ) DoMass = SaveMass
    IF( DoMass .AND. .NOT. ASSOCIATED( A % MassValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting mass matrix')
      DoMass = .FALSE. 
    END IF

    DoDamp = .FALSE.
    IF( PRESENT( SaveDamp ) ) DoDamp = SaveDamp
    IF( DoDamp .AND. .NOT. ASSOCIATED( A % DampValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting damp matrix')
      DoDamp = .FALSE. 
    END IF

    DoStiff = .TRUE.
    IF( PRESENT( SaveStiff ) ) DoStiff = SaveStiff
    IF( DoStiff .AND. .NOT. ASSOCIATED( A % Values ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting stiff matrix')
      DoStiff = .FALSE. 
    END IF

    Skip0 = .FALSE.
    IF(PRESENT(SkipZeros)) Skip0 = SkipZeros
        
    IF(.NOT. (DoStiff .OR. DoDamp .OR. DoMass ) ) THEN
      CALL Warn('CRS_PrintMatrix','Saving just the topology!')
    END IF
    
    IndStiff = 0
    IndDamp = 0
    IndMass = 0

    IF( DoStiff ) IndStiff = 1
    IF( DoDamp ) IndDamp = IndStiff + 1
    IF( DoMass ) IndMass = MAX( IndStiff, IndDamp ) + 1
    IndMax = MAX( IndStiff, IndDamp, IndMass )

    IF (Parallel.AND.Cnumbering) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)
  
      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, &
          A % Perm, A % Gorder, Owner )
    END IF

    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF
      DO j = A % Rows(i),A % Rows(i+1)-1

        col = A % Cols(j)
        IF(Parallel) THEN
          IF(Cnumbering) THEN
            col = A % Gorder(col)
          ELSE 
            col = A % ParallelInfo % GlobalDOFs(col)
          END IF
        END IF

        IF( DoStiff ) THEN
          Vals(IndStiff) = A % Values(j)
        END IF
        IF( DoDamp ) THEN
          Vals(IndDamp) = A % DampValues(j)
        END IF
        IF( DoMass ) THEN
          Vals(IndMass) = A % MassValues(j)
        END IF

        IF( Skip0 ) THEN
          IF(SUM(ABS(Vals(1:IndMax))) < EPSILON(val)) CYCLE
        END IF
          
        WRITE(1,'(I0,A,I0,A)',ADVANCE='NO') row,' ',col,' '
        IF( IndMax > 0 ) THEN
          WRITE(1,*) Vals(1:IndMax)          
        ELSE
          WRITE(1,'(A)') ' '
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE  PrintMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Prints the values of the right-hand-side vector to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintRHS( A, Parallel, CNumbering, SaveSum )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding matrix
    LOGICAL :: Parallel, CNumbering, SaveSum
!------------------------------------------------------------------------------
    INTEGER :: i, j, row
    REAL(KIND=dp) :: val, asum, rsum
    LOGICAL :: SaveRhs

    SaveRhs = ASSOCIATED(A % rhs)
    
    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF

      IF(SaveRhs) val = A % rhs(i)

      IF(SaveSum) THEN
        rsum = 0.0_dp
        asum = 0.0_dp
        DO j = A % Rows(i),A % Rows(i+1)-1
          rsum = rsum + A % Values(j)
          asum = asum + ABS(A % Values(j))
        END DO
      END IF
      
      WRITE(1,'(I0,A)',ADVANCE='NO') row,' '
      IF( SaveRhs .AND. SaveSum ) THEN
        WRITE(1,*) Val, rsum, asum
      ELSE IF( .NOT. SaveRhs .AND. SaveSum ) THEN
        WRITE(1,*) rsum, asum
      ELSE IF( ABS( Val ) <= TINY( val ) ) THEN
        WRITE(1,'(A)') '0.0'
      ELSE
        WRITE(1,*) Val
      END IF

    END DO

  END SUBROUTINE PrintRHS
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Solves a linear system and also calls the necessary preconditioning routines.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveLinearSystem( A, b, &
       x, Norm, DOFs, Solver, BulkMatrix )
!------------------------------------------------------------------------------
    USE EigenSolve

    REAL(KIND=dp) CONTIG :: b(:), x(:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: DOFs
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Matrix_t), OPTIONAL, POINTER :: BulkMatrix
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Relax,GotIt,Stat,ScaleSystem, EigenAnalysis, HarmonicAnalysis,&
               BackRotation, ApplyRowEquilibration, ApplyLimiter, Parallel, &
               SkipZeroRhs, SkipLoads, ComplexSystem, ComputeChangeScaled, &
               RecursiveAnalysis, CalcLoads, NotParallel
    INTEGER :: n,i,j,k,l,ii,m,DOF,istat,this,mn, AllocStat
    CHARACTER(:), ALLOCATABLE :: Method, Prec, SaveSlot
    TYPE(C_FUNPTR) :: Proc
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Px(:), &
                TempRHS(:), NonlinVals(:)
    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: s,Relaxation,Beta,Gamma,bnorm,Energy,xn,bn
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Aaid, Projector, MP
    REAL(KIND=dp), POINTER :: mx(:), mb(:), mr(:)
    TYPE(Variable_t), POINTER :: IterV
    LOGICAL :: NormalizeToUnity, AndersonAcc, AndersonScaled, NoSolve, Found
    REAL(KIND=dp), POINTER :: pv(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    CHARACTER(*), PARAMETER :: Caller = 'SolveLinearSystem'

    
    TARGET b, x 
    
    INTERFACE 
       SUBROUTINE VankaCreate(A,Solver)
          USE Types
          TYPE(Matrix_t) :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE VankaCreate

       SUBROUTINE CircuitPrecCreate(A,Solver)
          USE Types
          TYPE(Matrix_t), TARGET :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE CircuitPrecCreate

       SUBROUTINE FetiSolver(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE FetiSolver
 
       SUBROUTINE BlockSolveExt(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE BlockSolveExt
    END INTERFACE
!------------------------------------------------------------------------------

    Params => Solver % Values
     
    IF( ListGetLogical( Params,'Linear System Skip Complex',GotIt ) ) THEN
      CALL Info(Caller,'This time skipping complex treatment',Level=20)
      A % COMPLEX = .FALSE.
      ComplexSystem = .FALSE.
    ELSE
      ComplexSystem = ListGetLogical( Params, 'Linear System Complex', GotIt )
      IF ( GotIt ) A % COMPLEX = ComplexSystem
    END IF

    IF( ListGetLogical( Params,'Linear System Skip Scaling',GotIt ) ) THEN     
      CALL Info(Caller,'This time skipping scaling',Level=20)
      ScaleSystem = .FALSE.
    ELSE
      ScaleSystem = ListGetLogical( Params, 'Linear System Scaling', GotIt )
      IF ( .NOT. GotIt  ) ScaleSystem = .TRUE.
    END IF

    SkipLoads = ListGetLogical( Params,'Linear System Skip Loads',GotIt)
    
    
    IF( A % COMPLEX ) THEN
      CALL Info(Caller,'Assuming complex valued linear system',Level=6)
    ELSE
      CALL Info(Caller,'Assuming real valued linear system',Level=8)
    END IF

    Parallel = Solver % Parallel 
      
!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
    IF ( Parallel  ) THEN
      IF( .NOT. ASSOCIATED(A % ParMatrix) ) THEN
        CALL Info(Caller,'Creating parallel matrix structures',Level=8)
        A % Solver => Solver
        CALL ParallelInitMatrix( Solver, A )
        IF(A % ParallelInfo % NothingShared ) THEN
          CALL Info(Caller,'No dofs shared in parallel matrix!',Level=6)
        END IF
      ELSE
        CALL Info(Caller,'Using previously created parallel matrix structures!',Level=15)
      END IF      
      Parallel = ASSOCIATED(A % ParMatrix)       
    END IF

    IF( Parallel ) THEN
      CALL Info(Caller,'Assuming parallel linear system',Level=8)
    ELSE
      CALL Info(Caller,'Assuming serial linear system',Level=8)
    END IF  
        
    IF ( ListGetLogical( Solver % Values, 'Linear System Save',GotIt )) THEN
      saveslot = ListGetString( Solver % Values,'Linear System Save Slot', GotIt )
      IF(SaveSlot == 'linear solve') CALL SaveLinearSystem( Solver, A )
    END IF

!------------------------------------------------------------------------------

    n = A % NumberOfRows

    BackRotation = ListGetLogical(Params,'Back Rotate N-T Solution',GotIt)
    IF (.NOT.GotIt) BackRotation=.TRUE.
    BackRotation = BackRotation .AND. ASSOCIATED(Solver % Variable % Perm)

    IF ( Solver % Matrix % Lumped .AND. Solver % TimeOrder == 1 ) THEN
       Method = ListGetString( Params, 'Timestepping Method', GotIt)
       IF (  Method == 'runge-kutta' .OR. Method == 'explicit euler' ) THEN
         ALLOCATE(Diag(n), TempRHS(n))

         TempRHS= b(1:n)
         Diag = A % Values(A % Diag)

         IF( Parallel ) THEN
           CALL ParallelSumVector(A,Diag)
           CALL ParallelSumVector(A,TempRHS)
         END IF

         DO i=1,n
            IF ( ABS(Diag(i)) /= 0._dp ) x(i) = TempRHS(i) / Diag(i)
         END DO

         DEALLOCATE(Diag, TempRHS)

         IF (BackRotation) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
         Norm = ComputeNorm(Solver, n, x) 
         RETURN
       END IF
    END IF
    
!------------------------------------------------------------------------------
!  These definitions are needed if changing the iterative solver on-the-fly

    Solver % MultiGridSolver = ( ListGetString( Params, &
        'Linear System Solver', GotIt ) == 'multigrid' )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'MG Levels', GotIt, minv=1 ) )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'Multigrid Levels', GotIt, minv=1 ) )
    Solver % MultiGridLevel = Solver % MultigridTotal
!------------------------------------------------------------------------------

    EigenAnalysis = .FALSE.
    IF(.NOT. ListGetLogical( Params,'Skip Eigen Analysis',GotIt)) THEN
      EigenAnalysis = Solver % NOFEigenValues > 0 .AND. &
          ListGetLogical( Params, 'Eigen Analysis',GotIt )
    END IF
          
    HarmonicAnalysis = ( Solver % NOFEigenValues > 0 ) .AND. &
        ListGetLogical( Params, 'Harmonic Analysis',GotIt )
    
    ! These analyses types may require recursive strategies and may also have zero rhs
    RecursiveAnalysis = HarmonicAnalysis .OR. EigenAnalysis     

    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt ) 
    SkipZeroRhs = ListGetLogical( Params,'Skip Zero Rhs Test',GotIt )
#ifdef HAVE_FETI4I
    IF ( C_ASSOCIATED(A % PermonMatrix) ) THEN
      ScaleSystem = .FALSE.
      SkipZeroRhs = .TRUE.
    END IF
#endif

    IF ( .NOT. ( RecursiveAnalysis .OR. ApplyLimiter .OR. SkipZeroRhs ) ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))      
      IF ( bnorm <= TINY( bnorm) ) THEN
        CALL Info(Caller,'Solution trivially zero!',Level=5)
        x = 0.0d0

        ! Increase the nonlinear counter since otherwise some stuff may stagnate
        ! Normally this is done within ComputeChange
        iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        IF(ASSOCIATED(iterV)) THEN
          Solver % Variable % NonlinIter = iterV % Values(1)
          iterV % Values(1) = iterV % Values(1) + 1 
        END IF
        Solver % Variable % Norm = 0.0_dp
        Solver % Variable % NonlinConverged = 1
     
        RETURN
      END IF
    END IF

    IF ( Solver % MultiGridLevel == -1  ) RETURN
    
    ! Set the flags to false to allow recursive strategies for these analysis types, little dirty...
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.FALSE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.FALSE.)
    END IF


!------------------------------------------------------------------------------
!   If solving harmonic analysis go there:
!   --------------------------------------
    IF ( HarmonicAnalysis ) THEN
      CALL SolveHarmonicSystem( A, Solver )
    END IF


!   If solving eigensystem go there:
!   --------------------------------
    IF ( EigenAnalysis ) THEN
      IF ( ScaleSystem ) CALL ScaleLinearSystem(Solver, A )

      CALL SolveEigenSystem( &
          A, Solver %  NOFEigenValues, &
          Solver % Variable % EigenValues,       &
          Solver % Variable % EigenVectors, Solver )
      
      IF ( ScaleSystem ) CALL BackScaleLinearSystem( Solver, A, EigenScaling = .TRUE. ) 
      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )

      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      NormalizeToUnity = ListGetLogical( Solver % Values, &
          'Eigen System Normalize To Unity',Stat)         

      IF(NormalizeToUnity .OR. ListGetLogical( Solver % Values,  &
          'Eigen System Mass Normalize', Stat ) ) THEN

        CALL ScaleEigenVectors( A, Solver % Variable % EigenVectors, &
            SIZE(Solver % Variable % EigenValues), NormalizeToUnity ) 
      END IF

      CALL Info(Caller, 'Repointing '//I2S(Solver % Variable % DOFs)//&
          ' eigenvalue components for: '//TRIM(Solver % Variable % Name))

      IF (A % Complex) THEN
        DO k=1,Solver % Variable % DOFs
          str = ComponentName( Solver % Variable % Name, k )
          Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )

          IF( ASSOCIATED( Var ) ) THEN
            IF (Solver % Variable % DOFs == 2) THEN
              CALL Info('SolveEigenSystem', 'Eigenvalue component ' &
                  //I2S(k)//': '//TRIM(str))
              Var % EigenValues => Solver % Variable % EigenValues
              IF (.NOT. ASSOCIATED(Var % EigenVectors)) THEN
                ALLOCATE(Var % EigenVectors(SIZE(Var % EigenValues), &
                    SIZE(Solver % Variable % EigenVectors,2)), STAT=AllocStat)
              END IF

              DO i=1,SIZE(Var % EigenVectors,2)
                SELECT CASE(k)
                CASE(1)
                  Var % EigenVectors(:,i) = &
                      CMPLX(REAL(Solver % Variable % EigenVectors(:,i)), 0.0_dp, kind=dp)
                CASE(2)
                  ! This is the imaginary component as a real-valued array: 
                  Var % EigenVectors(:,i) = &
                      CMPLX(AIMAG(Solver % Variable % EigenVectors(:,i)), 0.0_dp, kind=dp)
                END SELECT
              END DO
            ELSE
              CALL Warn(Caller, 'A complex-valued system should have 2 DOFs')
            END IF
          END IF
        END DO
      END IF
      
      ! This is temporal (?) fix for a glitch where the complex eigen vector
      ! is expanded to one where real and complex parts follow each other. 
      IF( ListGetLogical( Solver % Values,'Expand Eigen Vectors', Stat ) ) THEN
        CALL ExpandEigenVectors( A, Solver % Variable % EigenVectors, &
            Solver % NOFEigenValues, Solver % Variable % dofs )
      END IF
        
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF
    
    ! We have solved {harmonic,eigen,constraint} system and no need to continue further
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.TRUE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.TRUE.)
      RETURN
    END IF


! Check whether b=0 since then equation Ax=b has only the trivial solution, x=0. 
! In case of a limiter one still may need to check the limiter for contact.
!-----------------------------------------------------------------------------
    IF( Parallel ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
    ELSE
      bnorm = SQRT(SUM(b(1:n)**2))
    END IF

    IF ( bnorm <= TINY( bnorm) .AND..NOT.SkipZeroRhs) THEN
      CALL Info(Caller,'Solution trivially zero!',Level=5)
      x = 0.0d0

      ! Increase the nonlinear counter since otherwise some stuff may stagnate
      ! Normally this is done within ComputeChange
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      IF(ASSOCIATED(iterV)) THEN
        Solver % Variable % NonlinIter = iterV % Values(1)
        iterV % Values(1) = iterV % Values(1) + 1
      END IF
      Solver % Variable % Norm = 0.0_dp
      Solver % Variable % NonlinConverged = 1

      RETURN
    END IF

    AndersonAcc = ListGetLogical( Params,'Nonlinear System Acceleration',GotIt ) 
    AndersonScaled = ListgetLogical( Params,'Nonlinear System Acceleration Scaled',GotIt ) 
    
    IF( AndersonAcc .AND. .NOT. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF(NoSolve) GOTO 120
    END IF
    
!   Convert rhs & initial value to the scaled system:
!   -------------------------------------------------
    IF ( ScaleSystem ) THEN
      CALL ScaleLinearSystem(Solver, A, b, x, &
          RhsScaling = (bnorm/=0._dp), ConstraintScaling=.TRUE. )
    END IF

    ComputeChangeScaled = ListGetLogical(Params,&
        'Nonlinear System Compute Change in Scaled System',GotIt)
    IF(.NOT.GotIt) ComputeChangeScaled = .FALSE.

    IF(ComputeChangeScaled) THEN
      ALLOCATE(NonlinVals(SIZE(x)))
      NonlinVals = x
      IF (ASSOCIATED(Solver % Variable % Perm)) & 
          CALL RotateNTSystemAll(NonlinVals, Solver % Variable % Perm, DOFs)
    END IF

    IF( AndersonAcc .AND. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF( NoSolve ) GOTO 110
    END IF
    
    IF( ListGetLogical( Params,'Linear System Normalize Guess',GotIt ) ) THEN
      CALL NormalizeInitialGuess() 
    ELSE IF( ListGetLogical( Params,'Linear System Nullify Guess',GotIt ) ) THEN
      CALL Info(Caller,'Nullifying initial guess!',Level=30)
      x(1:n) = 0.0_dp
    ELSE IF( ListGetLogical( Params,'Linear System Nullify First Guess',GotIt ) ) THEN
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter', UnfoundFatal=.TRUE.)
      i = NINT(iterV % Values(1))
      IF(i<=1) THEN
        CALL Info(Caller,'Nullifying first initial guess!',Level=30)
        x(1:n) = 0.0_dp
      END IF        
    END IF
    
    Method = ListGetString(Params,'Linear System Solver',GotIt)
    IF(.NOT. GotIt) THEN
      CALL Fatal(Caller,'Give "Linear System Solver", e.g. "iterative" or "direct"')
    END IF
    
    IF (Method=='multigrid' .OR. Method=='iterative' ) THEN
      Prec = ListGetString(Params,'Linear System Preconditioning',GotIt)
      IF( GotIt ) THEN
        CALL Info(Caller,'Linear System Preconditioning: '//TRIM(Prec),Level=8)
        CALL ResetTimer("Prec0-"//TRIM(Prec))
        IF( SEQL(Prec,'vanka') ) THEN
          IF(LEN(Prec)>=6) THEN
            i = ICHAR(Prec(6:6)) - ICHAR('0')
            CALL ListAddNewInteger( Params,'Vanka Mode',i) 
          END IF
          CALL VankaCreate(A,Solver)
        ELSE IF ( Prec=='circuit' ) THEN
          CALL CircuitPrecCreate(A,Solver)
#if 0 
          IF( ListGetLogical(Params,'Linear System Save', GotIt) ) THEN
            IF( ASSOCIATED( A % CircuitMatrix ) ) THEN
              CALL SaveLinearSystem(Solver, Ain = A % CircuitMatrix, LinSysName = "circuit")
            END IF
            IF( ASSOCIATED( Solver % Matrix % AddMatrix ) ) THEN
              CALL SaveLinearSystem(Solver, Ain = Solver % Matrix % AddMatrix, LinSysName = "addmatrix")
            END IF
          END IF
#endif
        END IF
        CALL CheckTimer("Prec0-"//TRIM(Prec),Level=8,Delete=.TRUE.)                  
      END IF
    END IF

    IF (Method=='iterative' ) THEN
      k = ListGetInteger(Params,'Linear System Dense Dofs',GotIt)
      IF(GotIt) A % ndeg = k

      IF( A % ndeg > 1 ) THEN
        IF( CRS_CheckStructuredDofs( A, A % ndeg) ) THEN
          CALL Fatal(Caller,'CRS matrix failed the dense test of size '//I2S(A % ndeg))
          A % ndeg = 0
        ELSE
          CALL Info(Caller,'CRS matrix passed the dense test of size '//I2S(A % ndeg),Level=12)
        END IF
      END IF
      IF( Parallel .AND. ASSOCIATED(A % ParMatrix) ) THEN
        IF( ASSOCIATED(A % ParMatrix % SplittedMatrix) ) &
          A % ParMatrix % SplittedMatrix % InsideMatrix % Ndeg = A % Ndeg
      END IF
    END IF
     
    IF( InfoActive(30) ) THEN
      CALL VectorValuesRange(A % values,SIZE(A % values),'A')
      pv => b
      CALL VectorValuesRange(pv,SIZE(pv),'b')
    END IF


    IF(ListGetLogical(Params, 'Linear System Use Rocalution', Found)) &
      Method = 'rocalution'

    NotParallel = .NOT. Parallel
    IF ( Parallel ) THEN
       IF ( A % ParallelInfo % NothingShared ) NotParallel = .TRUE.
    END IF

    IF ( NotParallel ) THEN
      IF(ListGetLogical(Params, 'Linear System Use Hypre', Found)) Method = 'hypre'

      CALL Info(Caller,'Serial linear System Solver: '//TRIM(Method),Level=8)
      
      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL IterSolver( A, x, b, Solver )
      CASE('feti')
        CALL Fatal(Caller,'Feti solver available only in parallel.')
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
      CASE('amgx')
        CALL AMGXSolver( A, x, b, Solver )
      CASE('rocalution')
        CALL ROCSolver( A, x, b, Solver )
      CASE('hypre')
        CALL SolveHypre( A, x, b, Solver )
      CASE('direct')
        CALL DirectSolver( A, x, b, Solver )        
      CASE DEFAULT        
        CALL Fatal(Caller,'Unknown "Linear System Solver": '//TRIM(Method))
      END SELECT
    ELSE
      CALL Info(Caller,'Parallel linear System Solver: '//TRIM(Method),Level=8)

      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL ParallelIter( A, A % ParallelInfo, DOFs, &
            x, b, Solver, A % ParMatrix )
      CASE('feti')
        CALL FetiSolver( A, x, b, Solver )
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
      CASE('amgx')
        CALL AMGXSolver( A, x, b, Solver )
      CASE('rocalution')
        CALL ROCSolver( A, x, b, Solver )
      CASE('direct')
        CALL DirectSolver( A, x, b, Solver )
      CASE DEFAULT        
        CALL Fatal(Caller,'Unknown "Linear System Solver": '//TRIM(Method))
      END SELECT
    END IF

    IF( InfoActive(30) ) THEN
      pv => x
      CALL VectorValuesRange(pv,SIZE(pv),'x')
    END IF
    
110 IF( AndersonAcc .AND. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF

    IF(ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, NonlinVals, Matrix=A, RHS=b )
      DEALLOCATE(NonlinVals)
    END IF

    IF ( ScaleSystem ) THEN
      CALL BackScaleLinearSystem( Solver, A, b, x, ConstraintScaling=.TRUE. )
    END IF

120 IF( AndersonAcc .AND. .NOT. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF
    
    Aaid => A
    IF (PRESENT(BulkMatrix)) THEN
      IF (ASSOCIATED(BulkMatrix) ) Aaid=>BulkMatrix
    END IF

    NodalLoads => NULL()
    IF(.NOT. SkipLoads ) THEN
      NodalLoads => VariableGet( Solver % Mesh % Variables, &
          GetVarName(Solver % Variable) // ' Loads' )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        ! Nodal loads may be allocated but the user may have toggled
        ! the 'calculate loads' flag such that no load computation should be performed.
        CalcLoads = ListGetLogical( Solver % Values,'Calculate Loads',GotIt )
        IF( .NOT. GotIt ) CalcLoads = .TRUE.
        IF( CalcLoads ) THEN
          CALL Info(Caller,'Calculating nodal loads for: '//&
              GetVarName(Solver % Variable),Level=6)
          CALL CalculateLoads( Solver, Aaid, x, Dofs, .TRUE., NodalLoads ) 
        END IF
      END IF
    END IF
      
    IF (BackRotation) THEN
      CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        CALL BackRotateNTSystem(NodalLoads % Values,NodalLoads % Perm,DOFs)
      END IF
    END IF

!------------------------------------------------------------------------------
    
!------------------------------------------------------------------------------
! Compute the change of the solution with different methods 
!------------------------------------------------------------------------------
    IF(.NOT.ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, Matrix=A, RHS=b )
    END IF
    Norm = Solver % Variable % Norm

    IF(.NOT. SkipLoads ) THEN
      NodalLoads => VariableGet( Solver % Mesh % Variables, &
          GetVarName(Solver % Variable) // ' Residual' )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        CalcLoads = ListGetLogical( Solver % Values,'Calculate Residual',GotIt )
        IF( .NOT. GotIt ) CalcLoads = .TRUE.
        IF( CalcLoads ) THEN
          CALL Info(Caller,'Calculating nodal residual',Level=6)
          CALL CalculateLoads( Solver, Aaid, x, Dofs, .FALSE., NodalLoads ) 
        END IF
      END IF
    END IF
    
!------------------------------------------------------------------------------
 
   Solver % Variable % PrimaryMesh => Solver % Mesh
   CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
         GetVarName(Solver % Variable) )
   
   IF ( ASSOCIATED( NodalLoads ) ) THEN
     NodalLoads % PrimaryMesh => Solver % Mesh
     CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
                  GetVarName(NodalLoads) )
   END IF

!------------------------------------------------------------------------------
! In order to be able to change the preconditioners or solvers the old matrix structures
! must be deallocated on request.

   IF( ListGetLogical( Params, 'Linear System Preconditioning Deallocate', GotIt) ) THEN
     ! ILU preconditioning
     IF( ASSOCIATED(A % ILUValues) ) THEN
       IF(  SIZE( A % ILUValues) /= SIZE(A % Values) ) &
           DEALLOCATE(A % ILUCols, A % ILURows, A % ILUDiag)
       DEALLOCATE(A % ILUValues)
     END IF

     ! Multigrid solver / preconditioner
     IF( Solver % MultigridLevel > 0 ) THEN
       Aaid => A 
       IF(ASSOCIATED( Aaid % Parent) ) THEN
         DO WHILE( ASSOCIATED( Aaid % Parent ) )
           Aaid => Aaid % Parent
         END DO
         DO WHILE( ASSOCIATED( Aaid % Child) )
           Aaid => Aaid % Child
           IF(ASSOCIATED(Aaid % Parent)) DEALLOCATE(Aaid % Parent )
           IF(ASSOCIATED(Aaid % Ematrix)) DEALLOCATE(Aaid % Ematrix )
         END DO
       END IF
     END IF
   END IF

 CONTAINS


   ! Sometimes the r.h.s. may abruptly diminish in value resulting to significant 
   ! convergence issues or it may be that the system scales linearly with the source. 
   ! This flag tries to improve on the initial guess of the linear solvers, and may 
   ! sometimes even result to the exact solution.
   !--------------------------------------------------------------------------------
   SUBROUTINE NormalizeInitialGuess()
     REAL(KIND=dp) :: xn, bn
     REAL(KIND=dp), ALLOCATABLE, TARGET :: TempVector(:)

     
     CALL Info(Caller,'Normalizing initial guess!',Level=30)

     ALLOCATE( TempVector(A % NumberOfRows) )

     IF ( Parallel ) THEN
       IF( .NOT. ALLOCATED( TempRHS ) ) THEN
         ALLOCATE( TempRHS(A % NumberOfRows) ); TempRHS=0._dp
       END IF

       Tempvector = 0._dp
       TempRHS(1:n) = b(1:n)
       CALL ParallelInitSolve( A, x, TempRHS, Tempvector )

       MP => ParallelMatrix(A,mx,mb,mr)
       mn = MP % NumberOfRows

       TempVector = 0._dp
       CALL ParallelMatrixVector( A, mx, TempVector )

       bn = ParallelDot( mn, TempVector, mb )
       xn = ParallelDot( mn, TempVector, TempVector )
       DEALLOCATE( TempRHS )
     ELSE
       CALL MatrixVectorMultiply( A, x, TempVector )
       xn = SUM( TempVector(1:n)**2 )
       bn = SUM( TempVector(1:n) * b(1:n) )
     END IF

     IF( xn > TINY( xn ) ) THEN
       x(1:n) = x(1:n) * ( bn / xn )
       WRITE( Message,'(A,ES12.3)') 'Linear System Normalizing Factor: ',bn/xn
       CALL Info(Caller,Message,Level=6) 
     END IF
     DEALLOCATE( TempVector )

   END SUBROUTINE NormalizeInitialGuess


   
!------------------------------------------------------------------------------
  END SUBROUTINE SolveLinearSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AMGXMatrixVectorMultiply( A, u, v, Solver )
!------------------------------------------------------------------------------
    USE iso_c_binding, only: C_INTPTR_T, C_CHAR, C_NULL_CHAR

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: u(*), v(*)

#ifdef HAVE_AMGX
    INTERFACE
      SUBROUTINE AMGXmv(AMGX, n, rows, cols, vals, b, x, nonlin_update ) BIND(C, Name="AMGXmv")

         USE Types
         USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

         IMPLICIT NONE

         INTEGER(KIND=C_INTPTR_T) :: AMGX
         REAL(KIND=dp) :: vals(*), b(*), x(*)
         INTEGER :: rows(*), cols(*), nonlin_update, n
      END SUBROUTINE AMGXmv
    END INTERFACE

    LOGICAL :: Found
    INTEGER :: nonlin_update, i

    nonlin_update = 0
    CALL AMGXmv( A % AMGXMV, A % NumberOfRows, A % Rows-1, A % Cols-1, &
                  A % Values, u, v, nonlin_update )
#else
    CALL Fatal('AMGXSolver', "AMGX doesn't seem to be included.")
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE AMGXMatrixVectorMultiply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Deterministic probe value for a global dof index, used by the
!> 'AMGX Verify Gather' check in AMGXSolver().
!>
!> Every partition must obtain exactly the same value for the same global index,
!> so this is an integer hash (the MINSTD multiplier) rather than anything that
!> would depend on the libm in use. The product stays well inside INTEGER(8):
!> global indices are below 2**31 and 48271*2**31 is about 1.0e14.
!------------------------------------------------------------------------------
  PURE FUNCTION AMGXProbeValue( g ) RESULT ( v )
!------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: g
    REAL(KIND=dp) :: v
!------------------------------------------------------------------------------
    INTEGER(KIND=8), PARAMETER :: Mult = 48271_8, Modulus = 2147483647_8
!------------------------------------------------------------------------------
    v = REAL( MOD( INT(g,8)*Mult, Modulus ), dp ) / REAL( Modulus, dp ) - 0.5_dp
!------------------------------------------------------------------------------
  END FUNCTION AMGXProbeValue
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AMGXSolver( A, x, b, Solver )
!------------------------------------------------------------------------------
    USE iso_c_binding, only: C_INTPTR_T, C_CHAR, C_NULL_CHAR

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: x(:), b(:)

#ifdef HAVE_AMGX
    INTERFACE
      SUBROUTINE AMGXSolve(AMGX, n, rows, cols, vals, b, x, &
              nonlin_update, config, comm, ng, part_vec, bnrm, solve_status ) BIND(C, Name="AMGXSolve")

         USE Types
         USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

         IMPLICIT NONE

         INTEGER(KIND=C_INTPTR_T) :: AMGX
         REAL(KIND=dp) :: vals(*), b(*), x(*), bnrm
         CHARACTER(KIND=C_CHAR) :: config(*)
         INTEGER :: rows(*), cols(*), nonlin_update, n, comm, ng, part_vec(*), solve_status
      END SUBROUTINE AMGXSolve
    END INTERFACE
#endif

    ! AMGX_SOLVE_STATUS values, mirrored from https://github.com/NVIDIA/AMGX/tree/main/include
    INTEGER, PARAMETER :: AMGX_SOLVE_SUCCESS = 0, AMGX_SOLVE_FAILED = 1, &
        AMGX_SOLVE_DIVERGED = 2, AMGX_SOLVE_NOT_CONVERGED = 3

    CHARACTER(KIND=C_CHAR) :: cfg(MAX_NAME_LEN)
    CHARACTER(LEN=MAX_NAME_LEN) :: config
    LOGICAL :: found, isparallel , nlin, DoFatal, VerifyGather, VerifyOnly, Rebuild
    INTEGER :: nonlin_update, i, j, n, lrow, me, you, solve_status
    REAL(KIND=dp)  :: bnrm

    TYPE(Matrix_t), POINTER :: Bm

    ! The collection state belongs to the matrix, not to this routine: see
    ! AMGXCollection_t in Types. > AC < is only a shorthand for it, so that
    ! several solvers using AMGX on different matrices keep their own.
    TYPE(AMGXCollection_t), POINTER :: AC
    INTEGER :: ng

    ! The collection of the partition-local matrices into whole owned rows below
    ! is plain Elmer + MPI: it needs no AMGX. Its correctness can therefore be
    ! checked on any machine, GPU or not, which is what 'AMGX Verify Gather'
    ! does. Without the library linked in that check is the only thing this
    ! routine can do, and it stops before the solve.
    VerifyGather = ListGetLogical( Solver % Values, 'AMGX Verify Gather', Found )
    VerifyOnly = .FALSE.
#ifndef HAVE_AMGX
    IF( .NOT. VerifyGather ) THEN
      CALL Fatal('AMGXSolver', "AMGX doesn't seem to be included.")
    END IF
    VerifyOnly = .TRUE.
    CALL Info('AMGXSolver','AMGX is not linked in: verifying the matrix gather only.',Level=4)
#endif

    solve_status = AMGX_SOLVE_SUCCESS

    nlin = ListGetLogical( Solver % Values, 'Linear System Refactorize', Found )
    IF ( nlin .OR. .NOT. Found ) THEN
      nonlin_update = 1
    ELSE
      nonlin_update = 0
    END IF

    config = ListGetString( Solver % Values, 'AMGX Config', Found )
    IF(.NOT. Found .AND. .NOT. VerifyOnly ) THEN
      CALL Fatal('AMGXSolver','Keyword > AMGX Config < is required!')
    END IF
    DO i=1,LEN_TRIM(config)
      cfg(i) = config(i:i)
    END DO
    cfg(i) = C_NULL_CHAR

    me =  Parenv % MyPe
    isParallel = Parenv % PEs>1

    IF(isParallel) THEN

      IF( .NOT. ASSOCIATED( A % AMGXColl ) ) ALLOCATE( A % AMGXColl )
      AC => A % AMGXColl
      ng = AC % ng

      BLOCK
        INTEGER, ALLOCATABLE :: part_vec_tmp(:)
        INTEGER, ALLOCATABLE :: Owner(:)

        INTEGER, ALLOCATABLE :: rRows(:),rSize(:),cBuf(:)
        REAL(KIND=dp), ALLOCATABLE :: vBuf(:)
        INTEGER :: status(MPI_STATUS_SIZE),ierr,i,j,k,l,lrow,you,rcnt,proc,totcnt,totnnz

        TYPE SendStuff_t
          INTEGER, ALLOCATABLE :: Size(:), Rows(:)
        END TYPE SendStuff_t
        TYPE(SendStuff_t), ALLOCATABLE :: SendStuff(:)

        ! Received (Bm row, global column) per entry, in the order the values
        ! arrive. Turned into > AC % RecvIdx < once Bm has its CRS structure.
        TYPE(IdxList_t), ALLOCATABLE :: RecvRow(:), RecvCol(:)
        INTEGER :: icnt, slot, nsend, nrecv
        REAL(KIND=dp), ALLOCATABLE :: sbuf(:), rbuf(:)
 
        IF(.NOT. ASSOCIATED(A % CollectionMatrix)) THEN
          ! Re-entering for this same matrix, its collection matrix having been
          ! dropped: throw the stale pattern away.
          IF(ALLOCATED(AC % APerm))THEN
            DEALLOCATE(AC % APerm, AC % iLPerm, AC % part_vec, AC % SendTo, &
                AC % GlobalToLocal)
            IF(ALLOCATED(AC % LocalMap)) DEALLOCATE(AC % LocalMap)
            IF(ALLOCATED(AC % SendIdx)) DEALLOCATE(AC % SendIdx)
            IF(ALLOCATED(AC % RecvIdx)) DEALLOCATE(AC % RecvIdx)
          END IF

          ! This is a different matrix from the one the cached pattern, if any,
          ! was built for.
          AC % PatternReady = .FALSE.
          AC % nnzA = -1
          AC % nrowsA = -1

          n = SIZE(A % ParallelInfo % GlobalDOFs)
          ALLOCATE( Owner(n), AC % APerm(n) )
          CALL ContinuousNumbering(A % ParallelInfo,A % Perm,AC % APerm,Owner)

          ng = 0
          DO i=1, A % NumberOfRows
            IF(A % ParallelInfo % NeighbourList(i) % Neighbours(1) == me) ng=ng+1
          END DO
          ng = NINT(ParallelReduction(1._dp*ng))
          AC % ng = ng

          ALLOCATE(AC % part_vec(ng), part_vec_tmp(ng));
          part_vec_tmp=-1; AC % part_vec=-1
          DO i=1,A % NumberOfRows
            part_vec_tmp(AC % APerm(i)) = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          END DO
          CALL MPI_ALLREDUCE( part_vec_tmp, AC % part_vec, ng, MPI_INTEGER, MPI_MAX, &
                            ELMER_COMM_WORLD, ierr )
          DEALLOCATE(part_vec_tmp)

          Bm => AllocateMatrix(); Bm % Format = MATRIX_LIST
          ALLOCATE(AC % SendTo(ParEnv % Pes),AC % GlobalToLocal(ng),AC % iLPerm(A % NumberOfRows))
        ELSE
          Bm => A % CollectionMatrix
        END IF

        ! The pattern has to be rebuilt on the first pass, if the structure of A
        ! changed under us, or if CRS_AddToMatrixElement() hit a column the
        ! gathered matrix does not have and flipped it back to list format.
        Rebuild = .NOT. AC % PatternReady .OR. Bm % Format == MATRIX_LIST .OR. &
            AC % nnzA /= SIZE(A % Values) .OR. AC % nrowsA /= A % NumberOfRows

        ! The decision has to be unanimous. The two branches below use different
        ! message patterns, so partitions disagreeing would deadlock -- and every
        ! term above is local: CRS_AddToMatrixElement() flips > Bm % Format <
        ! back to list format on the one partition that met an unknown column.
        i = 0
        IF( Rebuild ) i = 1
        Rebuild = ( NINT( ParallelReduction( 1._dp*i, 2 ) ) > 0 )

        IF( Rebuild ) THEN
          AC % iLPerm = 0
          LRow = 0
          AC % SendTo = 0
          AC % GlobalToLocal = 0
 
          IF (Bm % Format == MATRIX_CRS  ) Bm % Values = 0._dp
          DO i=1,A % NumberofRows
            you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
            IF ( you == me ) THEN
              LRow = LRow+1
              AC % iLPerm(LRow) = i
              DO j=A % Rows(i+1)-1, A % Rows(i),-1
                CALL AddToMatrixElement(Bm,LRow,AC % APerm(A  % Cols(j)), A % Values(j))
              END DO
              AC % GlobalToLocal(AC % APerm(i)) = LRow
            ELSE
              AC % SendTo(you+1) = AC % SendTo(you+1)+1
            END IF
          END DO

          ALLOCATE(SendStuff(ParEnv % Pes))
          DO i=1,ParEnv % PEs
            IF( i-1==me ) CYCLE
            IF(.NOT.ParEnv % IsNeighbour(i))  CYCLE
  
            ALLOCATE( SendStuff(i) % Rows(AC % SendTo(i)) )
            ALLOCATE( SendStuff(i) % Size(AC % SendTo(i)) )
          END DO
 
          AC % SendTo = 0
          DO i=1,a % NumberOfRows
            you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
            IF ( you /= me ) THEN
              AC % SendTo(you+1) = AC % SendTo(you+1)+1
              SendStuff(you+1) % Size(AC % SendTo(you+1))  = A % Rows(i+1)-A % Rows(i)
              SendStuff(you+1) % Rows(AC % SendTo(you+1))  = i
            END IF
          END DO

          ! Cache the flat list of value indices per neighbour, in exactly the
          ! order the rows are shipped below, so that later solves can send the
          ! values alone.
          IF(ALLOCATED(AC % SendIdx)) DEALLOCATE(AC % SendIdx)
          IF(ALLOCATED(AC % RecvIdx)) DEALLOCATE(AC % RecvIdx)
          ALLOCATE(AC % SendIdx(ParEnv % PEs), AC % RecvIdx(ParEnv % PEs))
          ALLOCATE(RecvRow(ParEnv % PEs), RecvCol(ParEnv % PEs))
          DO i=1,ParEnv % PEs
            IF( i-1==me .OR. .NOT. ParEnv % IsNeighbour(i) .OR. AC % SendTo(i)==0 ) THEN
              ALLOCATE(AC % SendIdx(i) % Ind(0))
              CYCLE
            END IF
            ALLOCATE(AC % SendIdx(i) % Ind(SUM(SendStuff(i) % Size)))
            l = 0
            DO j=1,AC % SendTo(i)
              k = SendStuff(i) % Rows(j)
              DO you=A % Rows(k), A % Rows(k+1)-1
                l = l+1
                AC % SendIdx(i) % Ind(l) = you
              END DO
            END DO
          END DO

          ! Ensure the MPI buffered-send pool is large enough for this loop's
          ! traffic: call CheckBuffer, similarily to RocSolver
          totcnt = SUM(AC % SendTo)
          totnnz = 0
          DO i=1,ParEnV % PEs
            IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE
            IF(AC % SendTo(i)>0) totnnz = totnnz + SUM(SendStuff(i) % Size)
          END DO
          CALL CheckBuffer( ParEnv % PEs*(1+MPI_BSEND_OVERHEAD) + 2*totcnt + 3*totnnz + &
                     (3*COUNT(AC % SendTo/=0) + 2*totcnt)*MPI_BSEND_OVERHEAD )

          DO i=1,ParEnV % PEs
            IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE

            CALL MPI_BSEND(AC % SendTo(i),1,MPI_INTEGER,i-1,1200,ELMER_COMM_WORLD, ierr)
            IF(AC % SendTo(i)==0) CYCLE

            CALL MPI_BSEND(AC % APerm(SendStuff(i) % Rows),AC % SendTo(i),MPI_INTEGER,i-1, &
                          1201,ELMER_COMM_WORLD,ierr )
  
            CALL MPI_BSEND( SendStuff(i) % Size,AC % SendTo(i),MPI_INTEGER,i-1, &
                          1202,ELMER_COMM_WORLD,ierr )
            DO j=1,AC % SendTo(i)
              k = SendStuff(i) % Rows(j)
              CALL MPI_BSEND(AC % APerm(A % Cols(A % Rows(k):A % Rows(k+1)-1)),SendStuff(i) % Size(j), &
                         MPI_INTEGER,i-1, 1203,ELMER_COMM_WORLD, ierr )

              CALL MPI_BSEND(A % Values(A % Rows(k):A % Rows(k+1)-1),SendStuff(i) % Size(j), &
                      MPI_DOUBLE_PRECISION,i-1,1204,ELMER_COMM_WORLD, ierr )
            END DO
          END DO

          DO i=1,ParEnV % PEs
            IF(i-1==me .OR. .NOT. ParEnv % IsNeighbour(i)) CYCLE

            CALL MPI_RECV(rcnt,1,MPI_INTEGER,MPI_ANY_SOURCE,1200,ELMER_COMM_WORLD,status,ierr)
            IF(rcnt==0) CYCLE

            proc = status(MPI_SOURCE)
            ALLOCATE( rRows(rcnt), rSize(rcnt) )

            CALL MPI_RECV(rRows,rcnt,MPI_INTEGER,proc,1201,ELMER_COMM_WORLD,status,ierr)
            CALL MPI_RECV(rSize,rcnt,MPI_INTEGER,proc,1202,ELMER_COMM_WORLD,status,ierr)

            ALLOCATE( RecvRow(proc+1) % Ind(SUM(rSize)), &
                      RecvCol(proc+1) % Ind(SUM(rSize)) )
            icnt = 0

            DO j=1,rcnt
              k = AC % GlobalToLocal(rRows(j))
              IF ( k==0 ) THEN
                PRINT*,Parenv % MyPE,proc, 'not mine then ?', rRows(j)
                ! Account for the values anyway: the neighbour ships them and
                ! the cached flat ordering has to line up with what it sends.
                DO l=1,rSize(j)
                  icnt = icnt+1
                  RecvRow(proc+1) % Ind(icnt) = 0
                  RecvCol(proc+1) % Ind(icnt) = 0
                END DO
                CYCLE
              END IF
              ALLOCATE(cBuf(rSize(j)), vBuf(rSize(j)))
              CALL MPI_RECV(cBuf,rSize(j),MPI_INTEGER,proc,1203, &
                        ELMER_COMM_WORLD,status,ierr )

              CALL MPI_RECV(vBuf,rSize(j),MPI_DOUBLE_PRECISION,proc, &
                      1204,ELMER_COMM_WORLD, status,ierr)

              DO l=1,rSize(j)
                CAll AddToMatrixElement(Bm,k,cBuf(l),vBuf(l))
                icnt = icnt+1
                RecvRow(proc+1) % Ind(icnt) = k
                RecvCol(proc+1) % Ind(icnt) = cBuf(l)
              END DO
              DEALLOCATE(cbuf, vbuf)
            END DO

            DEALLOCATE( rRows, rSize )
          END DO
          CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
 
          IF( Bm % Format == MATRIX_LIST ) THEN
            CALL List_toCRSMatrix(Bm)
            A % CollectionMatrix => Bm
          END IF

          ! Bm now has its final CRS structure, so the (row,column) pairs above
          ! can be resolved to value slots once and for all.
          IF( Bm % Format == MATRIX_CRS ) THEN
            IF(ALLOCATED(AC % LocalMap)) DEALLOCATE(AC % LocalMap)
            ALLOCATE(AC % LocalMap(SIZE(A % Values)))
            AC % LocalMap = 0

            DO lrow=1,Bm % NumberOfRows
              i = AC % iLPerm(lrow)
              IF( i<=0 ) CYCLE
              DO j=A % Rows(i), A % Rows(i+1)-1
                slot = CRS_Search( Bm % Rows(lrow+1)-Bm % Rows(lrow), &
                    Bm % Cols(Bm % Rows(lrow):Bm % Rows(lrow+1)-1), AC % APerm(A % Cols(j)) )
                IF( slot>0 ) AC % LocalMap(j) = slot + Bm % Rows(lrow) - 1
              END DO
            END DO

            DO i=1,ParEnv % PEs
              IF( .NOT. ALLOCATED(RecvRow(i) % Ind) ) THEN
                ALLOCATE(AC % RecvIdx(i) % Ind(0))
                CYCLE
              END IF
              nrecv = SIZE(RecvRow(i) % Ind)
              ALLOCATE(AC % RecvIdx(i) % Ind(nrecv))
              AC % RecvIdx(i) % Ind = 0
              DO l=1,nrecv
                k = RecvRow(i) % Ind(l)
                IF( k<=0 ) CYCLE
                slot = CRS_Search( Bm % Rows(k+1)-Bm % Rows(k), &
                    Bm % Cols(Bm % Rows(k):Bm % Rows(k+1)-1), RecvCol(i) % Ind(l) )
                IF( slot>0 ) AC % RecvIdx(i) % Ind(l) = slot + Bm % Rows(k) - 1
              END DO
            END DO

            AC % nnzA = SIZE(A % Values)
            AC % nrowsA = A % NumberOfRows
            AC % PatternReady = .TRUE.
            CALL Info('AMGXSolver','Collection pattern cached; later solves exchange values only.',Level=6)
          END IF

        ELSE IF( nonlin_update == 1 ) THEN
          ! Values-only refresh through the cached pattern: one contiguous
          ! message of doubles per neighbour, no indices, no per-row messages.
          Bm % Values = 0._dp

          DO j=1,SIZE(A % Values)
            IF( AC % LocalMap(j)>0 ) &
                Bm % Values(AC % LocalMap(j)) = Bm % Values(AC % LocalMap(j)) + A % Values(j)
          END DO

          totnnz = 0
          k = 0
          DO i=1,ParEnv % PEs
            nsend = SIZE(AC % SendIdx(i) % Ind)
            IF(nsend>0) THEN
              totnnz = totnnz + nsend
              k = k+1
            END IF
          END DO
          CALL CheckBuffer( 2*totnnz + (1+k)*MPI_BSEND_OVERHEAD )

          DO i=1,ParEnv % PEs
            nsend = SIZE(AC % SendIdx(i) % Ind)
            IF( nsend==0 ) CYCLE
            ALLOCATE(sbuf(nsend))
            sbuf = A % Values(AC % SendIdx(i) % Ind)
            CALL MPI_BSEND(sbuf,nsend,MPI_DOUBLE_PRECISION,i-1,1210, &
                ELMER_COMM_WORLD,ierr)
            DEALLOCATE(sbuf)
          END DO

          DO i=1,ParEnv % PEs
            nrecv = SIZE(AC % RecvIdx(i) % Ind)
            IF( nrecv==0 ) CYCLE
            ALLOCATE(rbuf(nrecv))
            CALL MPI_RECV(rbuf,nrecv,MPI_DOUBLE_PRECISION,i-1,1210, &
                ELMER_COMM_WORLD,status,ierr)
            DO l=1,nrecv
              IF( AC % RecvIdx(i) % Ind(l)>0 ) &
                  Bm % Values(AC % RecvIdx(i) % Ind(l)) = &
                      Bm % Values(AC % RecvIdx(i) % Ind(l)) + rbuf(l)
            END DO
            DEALLOCATE(rbuf)
          END DO
        END IF
      END BLOCK

      ! Verify the gather above against the parallel matrix itself.
      !
      ! The parallel matrix is a sum of partition-local contributions,
      ! A = SUM_p A_p, because every element is assembled on exactly one
      ! partition. Hence for any vector u on which all partitions agree,
      ! SUM_p (A_p u) is the exact global product A*u, and ParallelSumVector
      ! performs that sum. This gives a reference for the gathered matrix Bm
      ! that needs neither AMGX nor the SplittedMatrix glue.
      !
      ! Note that Elmer splits the parallel matrix by *column* ownership
      ! (SplitMatrix in SParIterSolver.F90), which lets the parallel
      ! matrix-vector product be formed from local u alone, whereas AMGX wants
      ! whole *rows* on the owning rank -- hence the gather. An orientation
      ! error between the two cancels for a symmetric matrix and does not for a
      ! general one, so run this on a strongly unsymmetric case, with a
      ! symmetric one as the control.
      IF( VerifyGather ) THEN
        BLOCK
          REAL(KIND=dp), ALLOCATABLE :: vref(:), vtst(:)
          REAL(KIND=dp) :: s, emax, vmax, relerr
          INTEGER :: ii, kk

          ALLOCATE( vref(A % NumberOfRows), vtst(Bm % NumberOfRows) )

          ! Reference: the local product of every partition, then summed.
          DO ii=1,A % NumberOfRows
            s = 0._dp
            DO kk=A % Rows(ii), A % Rows(ii+1)-1
              s = s + A % Values(kk) * AMGXProbeValue( AC % APerm(A % Cols(kk)) )
            END DO
            vref(ii) = s
          END DO
          CALL ParallelSumVector( A, vref )

          ! Candidate: the same product from the gathered whole rows. Bm % Cols
          ! already carry the global (AC % APerm) numbering.
          DO ii=1,Bm % NumberOfRows
            s = 0._dp
            DO kk=Bm % Rows(ii), Bm % Rows(ii+1)-1
              s = s + Bm % Values(kk) * AMGXProbeValue( Bm % Cols(kk) )
            END DO
            vtst(ii) = s
          END DO

          emax = 0._dp; vmax = 0._dp
          DO ii=1,Bm % NumberOfRows
            emax = MAX( emax, ABS( vtst(ii) - vref(AC % iLPerm(ii)) ) )
            vmax = MAX( vmax, ABS( vref(AC % iLPerm(ii)) ) )
          END DO
          emax = ParallelReduction( emax, 2 )
          vmax = ParallelReduction( vmax, 2 )

          relerr = emax
          IF( vmax > 0._dp ) relerr = emax / vmax

          WRITE(Message,'(A,ES12.5)') 'AMGX gather vs. parallel matrix-vector, relative error: ',relerr
          IF( relerr > 1.0e-9_dp ) THEN
            CALL Warn('AMGXSolver',Message)
            CALL Warn('AMGXSolver','The gathered matrix does not reproduce the parallel product!')
          ELSE
            CALL Info('AMGXSolver',Message,Level=4)
          END IF

          DEALLOCATE( vref, vtst )
        END BLOCK

        IF( VerifyOnly ) RETURN
      END IF

      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: bb(:),xb(:), r(:)

        n = Bm % NumberOfRows
        j = A  % NumberOfRows
        ALLOCATE(bb(n), xb(n), r(j) )

        r = b(1:j)
        CALL ParallelSumVector(A, r)

        DO i=1,n
          bb(i) = r(AC % iLPerm(i))
          xb(i) = x(AC % iLPerm(i))
        END DO
        bnrm = SQRT(ParallelReduction(SUM(bb**2)))

#ifdef HAVE_AMGX
        CALL AMGXSolve( A % AMGX, n, Bm % Rows-1, Bm % Cols-1, Bm % Values,  &
          bb, xb, nonlin_update, cfg, ELMER_COMM_WORLD, ng, AC % part_vec, bnrm, solve_status )
#endif

        x = 0
        DO i=1,n
          x(AC % iLPerm(i)) = xb(i)
        END DO
        CALL ParallelSumVector(A, x)
      END BLOCK
    ELSE
      IF( VerifyGather ) THEN
        CALL Info('AMGXSolver','Serial run: there is no gather to verify.',Level=4)
        IF( VerifyOnly ) RETURN
      END IF

      n = A % NumberOfRows
      bnrm = SQRT(SUM(b**2))

#ifdef HAVE_AMGX
      CALL AMGXSolve( A % AMGX, n, A % Rows-1, A % Cols-1, &
        A % Values, b, x, nonlin_update, cfg, ELMER_COMM_WORLD, n, &
           A % Diag, bnrm, solve_status ) ! <--- a % diag  for dummy
#endif
    END IF

    IF( solve_status == AMGX_SOLVE_SUCCESS ) THEN
      IF( ASSOCIATED( Solver % Variable ) ) Solver % Variable % LinConverged = 1
    ELSE
      CALL Info('AMGXSolver','Returned status: '//I2S(solve_status),Level=15)
      SELECT CASE( solve_status )
      CASE( AMGX_SOLVE_FAILED )
        CALL NumericalError('AMGXSolver','AMGX solver failed.')
      CASE( AMGX_SOLVE_DIVERGED )
        CALL NumericalError('AMGXSolver','AMGX solve diverged.')
      CASE( AMGX_SOLVE_NOT_CONVERGED )
        DoFatal = ListGetLogical( Solver % Values, 'Linear System Abort Not Converged', Found )
        IF(.NOT. Found ) DoFatal = .TRUE.
        IF( DoFatal ) THEN
          CALL NumericalError('AMGXSolver','Too many iterations were needed.')
        ELSE
          CALL Info('AMGXSolver','AMGX solve did not converge to tolerance',Level=6)
        END IF
      CASE DEFAULT
        CALL Warn('AMGXSolver','AMGX solver returned unrecognized status '//I2S(solve_status))
      END SELECT
      IF( ASSOCIATED( Solver % Variable ) ) Solver % Variable % LinConverged = 0
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AMGXSolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE ROCSolver( A, x, b, Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: x(:), b(:)

#ifdef HAVE_ROCALUTION

    ! Interfaces to c++ code calling ROCalution library
    ! -------------------------------------------------
    INTERFACE
      SUBROUTINE ROCSerialSolve(n, rows, cols, vals, b, x, nonlin_update, &
            imethod, prec, maxiter, tol, schur_n, schur_rows, schur_cols, &
            schur_vals, dofs ) BIND(C, Name="ROCSerialSolve")
        USE Types
        USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

        IMPLICIT NONE
        REAL(KIND=dp) :: vals(*), b(*), x(*), tol, schur_vals(*)
        INTEGER :: rows(*), cols(*), nonlin_update, n, imethod, prec, maxiter
        INTEGER :: schur_n, schur_rows(*), schur_cols(*), dofs
      END SUBROUTINE ROCSerialSolve


      SUBROUTINE ROCParallelSolve(gn, n, rows, cols, vals, b, x, bnrm, goffset, &
         fcomm, imethod, prec, maxiter, tol) BIND(C, Name="ROCParallelSolve")
        USE Types
        USE ISO_C_BINDING, ONLY: C_CHAR, C_INTPTR_T

        IMPLICIT NONE
        REAL(KIND=dp) :: vals(*), b(*), x(*), bnrm, tol
        INTEGER :: gn, n, rows(*), cols(*), goffset(*), fcomm, imethod, prec, maxiter
      END SUBROUTINE ROCParallelSolve
    END INTERFACE

    ! local variables:
    ! ----------------
    LOGICAL :: Found, isParallel, Refactorize
    INTEGER :: nonlin_update, i, j, k,l, m, n, p,q,gn, me

    TYPE(Matrix_t), POINTER ::Rmatrix

    INTEGER, POINTER ::  aPerm(:), iLperm(:), gOffset(:), iRows(:), iCols(:)
    REAL(KIND=dp), POINTER :: iVals(:)
    REAL(KIND=dp), ALLOCATABLE :: dBuf(:)
    INTEGER, ALLOCATABLE :: Owner(:), SendTo(:), iBuf(:), tOffset(:), rRows(:), rSize(:)

    INTEGER :: status(MPI_STATUS_SIZE),ierr,lrow,you,rcnt,proc

    INTEGER :: buf_size, procs

    TYPE(Matrix_arr_t), POINTER :: iMatrix(:)

    ! Define these to get somewhat shorter MPI subroutine calls:
    ! -----------------------------------------------------------
    INTEGER :: xmpi_comm
    INTEGER :: xmpi_src  = MPI_ANY_SOURCE
    INTEGER :: xmpi_max  = MPI_MAX
    INTEGER :: xmpi_int  = MPI_INTEGER
    INTEGER :: xmpi_dbl  = MPI_DOUBLE_PRECISION


    REAL(KIND=dp) :: TOL
    INTEGER :: Imethod, Prec, MaxIter, ILULevel, own_n

    TYPE(ValueList_t), POINTER :: Params

    TYPE SendStuff_t
      INTEGER, ALLOCATABLE :: Size(:), Rows(:)
    END TYPE SendStuff_t
    TYPE(SendStuff_t), ALLOCATABLE :: SendStuff(:)

    REAL(KIND=dp)  :: rt

    Params => Solver % Values

    ! Extract some controls to ROCalution from the simulation control info:
    ! ---------------------------------------------------------------------
    nonlin_update = 1
    Refactorize = ListGetLogical( Params, 'Linear System Refactorize', Found )
    IF ( Found .AND. .NOT. Refactorize ) nonlin_update = 0

    SELECT CASE(ListGetString(Params,'Linear System Iterative Method',Found))
      CASE('cg')
         Imethod = 0;
      CASE('bicgstab')
         Imethod = 1;
      CASE('bicgstabl')
         Imethod = 2;
      CASE('gmres')
         Imethod = 3;
      CASE('fgmres')
         Imethod = 4;
      CASE DEFAULT
         Imethod = 0;
    END SELECT

    SELECT CASE(ListGetString(Params,'Linear System Preconditioning',Found))
      CASE('jacobi')
        Prec = 0;
      CASE('sgs')
        Prec = 1;
      CASE('ilu')
        Prec = 2; ILULevel = 0
      CASE('ilu0')
       Prec = 2; ILULevel = 0
      CASE('ilu1')
       Prec = 2; ILULevel = 1
      CASE('ilu2')
       Prec = 2; ILULevel = 2
      CASE('schur')
       Prec = 3; ILULevel = 0
      CASE DEFAULT
       Prec = 0;
    END SELECT

    MaxIter = ListGetInteger(Params,'Linear System Max Iterations')
    TOL = ListGetCReal(Params,'Linear System Convergence Tolerance')

    ! ---------------------------------------------------------------------

    n = A % NumberOfRows

    procs =  Parenv % PEs
    isParallel = procs>1

    IF(isParallel) THEN
            rt = RealTime()
      me    =  Parenv % MyPe
      xmpi_comm = ELMER_COMM_WORLD

      IF (.NOT.ASSOCIATED(A % ParMatrix)) THEN
        A % Solver => Solver
        CALL ParallelInitMatrix(Solver,A)
      END IF

      A % Solver => Solver
      CALL SetMatrixParEnv( A )

      ! Enforce continuous ascending numbering as required by ROCalution
      ! ----------------------------------------------------------------
      IF (.NOT. ASSOCIATED(A % RocParams % Rmatrix)) THEN
        n = SIZE(A % ParallelInfo % GlobalDOFs)
        ALLOCATE( Owner(n), Aperm(n), tOffset(0:procs), gOffset(0:procs) )

        tOffset = -1
        CALL ContinuousNumbering(A % ParallelInfo,A % Perm,aPerm,Owner,gStart=tOffset(me))

        gOffset = -1
        CALL MPI_ALLREDUCE(tOffset,gOffset,procs,xmpi_int,xmpi_max,xmpi_comm,ierr)

        own_n = SUM(Owner)
        gOffset(procs) = ParallelReduction(own_n)

        ALLOCATE(iLPerm(A % NumberOfRows))

        Rmatrix => AllocateMatrix();
        Rmatrix % Format = MATRIX_LIST
        Rmatrix % ListMatrix => List_AllocateMatrix(own_n)

        ALLOCATE(iMatrix(0:ParEnv % PEs-1))
        DO proc=0,ParEnv % PEs-1
          iMatrix(proc) % M => AllocateMatrix()
          iMatrix(proc) % M % Format = MATRIX_LIST
          iMatrix(proc) % M % ListMatrix => List_AllocateMatrix(own_n)
        END DO

        A % RocParams % Rmatrix => Rmatrix
        A % RocParams % CntPerm => aPerm
        A % RocParams % LocPerm => iLperm
        A % RocParams % gOffset => gOffset
        A % RocParams % iMatrix => iMatrix
      ELSE
        Rmatrix => A % RocParams % Rmatrix
        aPerm   => A % RocParams % CntPerm
        iLPerm  => A % RocParams % LocPerm
        gOffset => A % RocParams % gOffset
        iMatrix => A % RocParams % iMatrix
      END IF

      ! Complete the matrix rows such that each partition has full rows of the 'owned' dofs
      ! -----------------------------------------------------------------------------------
      IF( Rmatrix % Format == MATRIX_LIST .OR. nonlin_update==1 ) THEN
 
        IF (Rmatrix % Format == MATRIX_CRS  ) Rmatrix % Values = 0._dp

        ! Create inside matrix + count rows with values to send for each neighbour
        ! -------------------------------------------------------------------------
        ALLOCATE(SendTo(0:procs-1))
        iLPerm = 0
        LRow = 0
        SendTo = 0
        DO i=1,A % NumberofRows
          you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          IF ( you == me ) THEN
            lRow = lRow + 1
            iLPerm(lRow) = i

            l = A % Rows(i+1) - A % Rows(i)
            dbuf = A % Values(A % Rows(i):A % Rows(i+1)-1)
            ibuf = aPerm(A % Cols(A % Rows(i):A % Rows(i+1)-1))
            CALL SortF(l, ibuf, dbuf)

            IF ( Rmatrix % Format == MATRIX_LIST ) THEN
              CALL List_AddMatrixRow(Rmatrix % ListMatrix, lRow, l, &
                       ibuf, dbuf, SortedInput=.TRUE. )
            ELSE
!             DO j=A % Rows(i+1)-1, A % Rows(i),-1
!               CALL AddToMatrixElement(Rmatrix, lRow, aPerm(A  % Cols(j)), A % Values(j))
!             END DO
              l = 0
              k = RMatrix % Rows(lRow)
              DO j = A % Rows(i), A % Rows(i+1)-1
                l = l + 1
                DO WHILE( ibuf(l) /= RMatrix % Cols(k) )
                  k = k + 1
                END DO
                RMatrix % Values(k) = dbuf(l)
                k = k + 1
              END DO
            ENDIF
          ELSE
            SendTo(you) = SendTo(you)+1
          END IF
        END DO

        ALLOCATE(SendStuff(0:procs-1))
        DO proc=0,procs-1
          IF( proc==me .OR. .NOT. ParEnv % IsNeighbour(proc+1) ) CYCLE
          ALLOCATE( SendStuff(proc) % Rows(SendTo(proc)) )
          ALLOCATE( SendStuff(proc) % Size(SendTo(proc)) )
        END DO
 
        ! Count number of columns of each neighbour's rows to be sent
        ! -----------------------------------------------------------
        SendTo   = 0
        buf_size = 0
        DO i=1,a % NumberOfRows
          you = A % ParallelInfo % NeighbourList(i) % Neighbours(1)
          IF ( you /= me ) THEN
            SendTo(you) = SendTo(you)+1
            SendStuff(you) % Size(Sendto(you))  = A % Rows(i+1)-A % Rows(i)
            SendStuff(you) % Rows(Sendto(you))  = i
            buf_size = buf_size + A % Rows(i+1) - A % Rows(i)
          END IF
        END DO
        CALL CheckBuffer(ParEnv % PEs*(4+4*MPI_BSEND_OVERHEAD) + 4*buf_size)

        ! Send data to neighbours
        ! -----------------------
        DO proc=0,procs-1
          IF(proc==me .OR. .NOT. ParEnv % IsNeighbour(proc+1)) CYCLE

          CALL MPI_BSEND(SendTo(proc),1,xmpi_int,proc,1200,xmpi_comm,ierr)
          IF(Sendto(proc)==0) CYCLE

          ibuf = aPerm(SendStuff(proc) % Rows)
          CALL MPI_BSEND(ibuf,SendTo(proc),xmpi_int,proc,1201,xmpi_comm,ierr)

          ibuf = SendStuff(proc) % Size
          CALL MPI_BSEND(ibuf,SendTo(proc),xmpi_int,proc,1202,xmpi_comm,ierr)

          DO j=1,SendTo(proc)
            k = SendStuff(proc) % Rows(j)
            l = SendStuff(proc) % Size(j)
            dBuf = A % Values(A % Rows(k):A % Rows(k+1)-1)
            iBuf = aPerm(A % Cols(A % Rows(k):A % Rows(k+1)-1))
            CALL SortF(l, ibuf, dbuf)
            CALL MPI_BSEND(dBuf,l,xmpi_dbl,proc,1203,xmpi_comm,ierr)
            IF (Rmatrix % Format == MATRIX_LIST ) THEN
              CALL MPI_BSEND(iBuf,l,xmpi_int,proc,1204,xmpi_comm,ierr)
            END IF
          END DO
        END DO

        ! receive data from neighbours
        ! ----------------------------
        ALLOCATE( rRows(100), rSize(100) )
        DO i=0,procs-1
          IF(i==me .OR. .NOT. ParEnv % IsNeighbour(i+1)) CYCLE

          CALL MPI_RECV(rcnt,1,xmpi_int,xmpi_src,1200,xmpi_comm,status,ierr)
          IF(rcnt==0) CYCLE

          IF(rcnt > SIZE(rRows)) THEN
            DEALLOCATE(rRows,Rsize)
            ALLOCATE( rRows(rcnt), rSize(rcnt) )
          END IF

          proc = status(MPI_SOURCE)
          CALL MPI_RECV(rRows,rcnt,xmpi_int,proc,1201,xmpi_comm,status,ierr)
          CALL MPI_RECV(rSize,rcnt,xmpi_int,proc,1202,xmpi_comm,status,ierr)
          DO j=1,rcnt
            k = rRows(j) - gOffset(me)

            IF ( k<=0 .OR. k>gOffset(me+1)-gOffset(me) ) THEN
              PRINT*,Parenv % MyPE,proc, 'not mine then ?', rRows(j), gOffset(me), gOffset(me+1)
              CYCLE
            END IF

            IF(rSize(j) > SIZE(iBuf)) THEN
              DEALLOCATE(iBuf,dBuf)
              ALLOCATE( iBuf(rSize(j)), dBuf(rSize(j)) )
            END IF
            CALL MPI_RECV(dBuf,rSize(j),xmpi_dbl,proc,1203,xmpi_comm,status,ierr)

            IF ( RMatrix % Format == MATRIX_LIST ) THEN
              CALL MPI_RECV(iBuf,rSize(j),xmpi_int,proc,1204,xmpi_comm,status,ierr)
              CALL List_AddMatrixRow(iMatrix(proc) % M % ListMatrix,k, &
                       rSize(j),iBuf,dBuf,SortedInput=.TRUE.)
            ELSE
              q = 0
              DO l=iMatrix(proc) % M % Rows(k), iMatrix(proc) % M % Rows(k+1)-1
                q = q + 1
                m = iMatrix(proc) % M % Cols(l)
                Rmatrix % Values(m) = RMatrix % Values(m) + dBuf(q)
              END DO
            END IF
          END DO
        END DO

        ! ----------

        IF(Rmatrix % Format == MATRIX_LIST) THEN

          DO proc=0,procs-1
            IF ( proc==me .OR. .NOT. ParEnv % IsNeighbour(proc+1)) CYCLE

            CALL List_toCRSMatrix(iMatrix(proc) % M)
            DO i=1,iMatrix(proc) % M % NumberOfRows
              iRows => iMatrix(proc) % M % Rows
              iCols => iMatrix(proc) % M % Cols
              iVals => iMatrix(proc) % M % Values

              l = iRows(i+1) - iRows(i)
              IF (l>0) THEN
                CALL List_AddMatrixRow( Rmatrix % ListMatrix,i,l, &
                   iCols(iRows(i):iRows(i+1)-1), iVals(iRows(i):iRows(i+1)-1), SortedInput=.TRUE.)
              END IF
            END DO
          END DO

          CALL List_toCRSMatrix(Rmatrix)

          DO proc=0,procs-1
            IF ( iMatrix(proc) % M % NumberOfRows <= 0 ) CYCLE

            iRows => iMatrix(proc) % M % Rows
            iCols => iMatrix(proc) % M % Cols

            DO i=1,iMatrix(proc) % M % NumberOfRows
              l = RMatrix % Rows(i)
              DO j=iRows(i), iRows(i+1)-1
                DO WHILE(Rmatrix % Cols(l) /= iCols(j))
                   l=l+1 
                END DO
                iCols(j) = l
              END DO
            END DO
          END DO
        END IF

!       print*,'ct time: ', realtime()-rt
      END IF

      n = Rmatrix % NumberOfRows
      gn = ParallelReduction(n);


      !  the linear equation solver
      ! ---------------------------
      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: pb(:),px(:), r(:)
        REAL(KIND=dp) :: bnrm

        n = Rmatrix % NumberOfRows
        j = A  % NumberOfRows
        ALLOCATE(pb(n), px(n), r(j) )

        ! complete the RHS of the partition dofs
        !---------------------------------------
        r = b(1:j)
        CALL ParallelSumVector(A, r)

        ! Extract 'owned' dofs
        !---------------------
        DO i=1,n
          pb(i) = r(iLPerm(i))
          px(i) = x(iLPerm(i))
        END DO

        bnrm = SQRT(ParallelReduction(SUM(pb**2)))
        IF(bnrm <AEPS) bnrm=1;

        !  call ROCalution
        ! ----------------
        CALL ROCParallelSolve( gn, n, Rmatrix % Rows-1, Rmatrix % Cols-1, &
             Rmatrix % Values, pb, px, bnrm, gOffset, A % comm, imethod, prec, maxiter, tol )

        ! distribute the result such that all dofs present (even shared and not 'owned')
        ! in a partition have consistent values
        ! ------------------------------------------------------------------------------
        x = 0
        DO i=1,n
          x(iLPerm(i)) = px(i)
        END DO
        CALL ParallelSumVector(A, x)
      END BLOCK

      ! Cleanup, remains to be reconsidered for optimizations
      ! -----------------------------------------------------
      !CALL FreeMatrix(Rmatrix);
      !DEALLOCATE(APerm,ILperm,gOffset)

      !A % RocParams % Rmatrix => Null()
      !A % RocParams % CntPerm => Null()
      !A % RocParams % LocPerm => Null()
      !A % RocParams % gOffset => Null()
    ELSE
      ! Serial case: call the linear solver
      ! -----------------------------------
      BLOCK
        TYPE(Variable_t), POINTER :: SchurV
        TYPE(Matrix_t), POINTER :: Schur 
        REAL(KIND=dp) :: ddum(1)
        INTEGER :: i, j, k, l, dofs, idum(1)

        Schur => NULL()
        dofs = Solver % Variable % DOFs
        IF (prec==3) THEN
          IF ( ListGetLogical( Solver % Values, 'Create Schur Approximation Matrix', Found) ) THEN
            Schur => XCreateSchurApproximation(A)
          ELSE
            SchurV => VariableGet( Solver % Mesh % Variables, 'Schur' )          
            IF ( ASSOCIATED(SchurV) ) Schur => SchurV % Solver % Matrix
          END IF
        END IF

        IF ( ASSOCIATED(Schur) ) THEN
          CALL ROCSerialSolve( n, A % Rows-1, A % Cols-1, A % Values, b, x, &
              nonlin_update, imethod, prec, maxiter, tol, Schur % numberOfRows, &
              Schur % Rows-1, Schur % cols-1, Schur % Values, dofs )
          CALL FreeMatrix( Schur)
        ELSE
          CALL ROCSerialSolve( n, A % Rows-1, A % Cols-1, A % Values, b, x, &
              nonlin_update, imethod, prec, maxiter, tol, 0, idum, idum, ddum, dofs)
        END IF
      END BLOCK
    END IF
#else
    CALL Fatal('ROCSolver', "Rocalution doesn't seem to be included.")
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE ROCSolver
!------------------------------------------------------------------------------


  ! Create matrix S=P((diag(A))^-1)Q
  !------------------------------------------------------------------------  
  FUNCTION XCreateSchurApproximation(A) RESULT ( S ) 

    TYPE(Matrix_t) :: A
    TYPE(Matrix_t), POINTER :: P, Q
    TYPE(Matrix_t), POINTER :: S

    INTEGER :: n, nc, i, j, k, l, j2, k2
    REAL(KIND=dp) :: val
    LOGICAL :: Found
    
    CALL Info('CreateSchurApproximation','Creating Shcur complement for preconditioning!',Level=20)

    NULLIFY(S)
!   IF(.NOT. ASSOCIATED(P) .OR. .NOT. ASSOCIATED(Q)) THEN
!     CALL Info('CreateSchurApproximation','Constraint matrix not associated!')
!     RETURN
!   END IF
    S => AllocateMatrix()
    
    nc = CoordinateSystemDimension() + 1
    n = A % NumberOfRows / nc
    IF(n == 0) THEN
      CALL Info('CreateSchurApproximation','No rows in Constraint matrix!')
      RETURN
    END IF

    S % FORMAT = MATRIX_LIST
      
    ! Add the corner entry to give the max size for list.  
    CALL List_AddToMatrixElement(S % ListMatrix, n, n, 0.0_dp ) 

    l = 0
    DO i=nc,n*nc,nc
      l = l + 1
      DO j=A % Rows(i),A % Rows(i+1)-1
        k = A % Cols(j)
        IF (MOD(k,nc)==0) CYCLE

        val = A % Values(j) / A % Values(A % Diag(k))
        DO j2=A % Rows(k)+nc-1,A % Rows(k+1)-1,nc
          k2 = A % Cols(j2)
          CALL List_AddToMatrixElement(S % ListMatrix, l, (k2-1)/nc+1, -val * A % Values(j2) )
        END DO
      END DO
    END DO

    CALL List_toCRSMatrix(S)
    
    val = 1.0_dp ! SIZE(S % Values) / SIZE(P % Values)
    WRITE(Message,*) 'Schur matrix increase factor: ',val, S % NumberOfrows, SUM(S % Values)
    CALL Info('CreateSchurApproximation',Message)
  END FUNCTION XCreateSchurApproximation





!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  SUBROUTINE LinearSystemResidual( A, b, x, r )
!------------------------------------------------------------------------------

    REAL(KIND=dp) CONTIG :: b(:)   
    REAL(KIND=dp) CONTIG :: x(:)   
    TYPE(Matrix_t), POINTER :: A   
    REAL(KIND=dp), POINTER :: r(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec

    INTEGER :: i,n,nn 

    n = A % NumberOfRows 

    IF (Parenv % Pes>1) THEN
      CALL ParallelInitSolve(A,x,b,r)
      CALL ParallelMatrixVector(A,x,r,.TRUE.)
    ELSE
      CALL MatrixVectorMultiply( A, x, r)
    END IF

    DO i=1,n
      r(i) = b(i) - r(i)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LinearSystemResidual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  FUNCTION LinearSystemMaskedResidualNorm( A, b, x, ActiveRow, ActiveCol ) RESULT ( Nrm )

    REAL(KIND=dp) CONTIG :: b(:)
    REAL(KIND=dp) CONTIG :: x(:)
    TYPE(Matrix_t) :: A
    LOGICAL, DIMENSION(:) :: ActiveRow(:), ActiveCol(:)
    REAL(KIND=dp) :: Nrm
    
    REAL(KIND=dp), ALLOCATABLE :: r(:)
    INTEGER :: i,n,totn
    REAL(KIND=dp) :: r2sum

    n = A % NumberOfRows 

    ALLOCATE(r(n))
   
    IF (Parenv % Pes>1) THEN
      CALL Fatal('LinearSystemMaskedResidualNorm','Not implemented in parallel yet!')
!      CALL ParallelMatrixVector(A, x, r, .TRUE.)
    ELSE
      CALL MaskedMatrixVectorMultiply( A, x, r, ActiveRow, ActiveCol )
    END IF

    DO i=1,n
      IF( ActiveRow(i) ) THEN
        r(i) = b(i) - r(i)
      END IF
    END DO

    totn = ParallelReduction(n)

    r2sum = SUM( r**2 )
    Nrm = SQRT( ParallelReduction(r2sum) / totn )

    DEALLOCATE( r ) 
    
  END FUNCTION LinearSystemMaskedResidualNorm



  FUNCTION HaveRestrictionMatrix( A ) RESULT( HaveConstraint ) 

    TYPE(Matrix_t), POINTER :: A
    LOGICAL :: HaveConstraint

    INTEGER :: n

    IF( .NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('HaveRestrictionMatrix','Matrix A not associated!')
    END IF

    n = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  THEN
      IF ( A % ConstraintMatrix % NumberOFRows > 0 ) n = n + 1 
    END IF
         
    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows > 0 ) n = n + 1
    END IF
    
    n = ParallelReduction(n)
    HaveConstraint = ( n > 0 ) 
    
  END FUNCTION HaveRestrictionMatrix

  
  
!------------------------------------------------------------------------------
!> Solve a system. Various additional utilities are included and 
!> naturally a call to the linear system solver.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveSystem( A,ParA,b,x,Norm,DOFs,Solver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) CONTIG, TARGET :: b(:)   !< The RHS vector
    REAL(KIND=dp) CONTIG :: x(:)   !< Previous solution on entry, new solution on exit (hopefully)
    REAL(KIND=dp) :: Norm          !< L2 Norm of solution
    TYPE(Matrix_t), POINTER :: A   !< The coefficient matrix
    INTEGER :: DOFs                !< Number of degrees of freedom per node for this equation
    TYPE(Solver_t), TARGET :: Solver                 !< Holds various solver options.
    TYPE(SParIterSolverGlobalD_t), POINTER :: ParA   !< holds info for parallel solver, 
                                                     !< if not executing in parallel this is just a dummy.
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh, SaveMEsh
    LOGICAL :: Relax, Found, NeedPrevSol, Timing, ResidualMode, &
        RestrictionMode, BlockMode, GloNum, FirstLoop
    INTEGER :: n,i,j,k,l,m,istat,nrows,ncols,colsj,rowoffset
    CHARACTER(:), ALLOCATABLE :: Method, VariableName
    TYPE(C_FUNPTR) :: Proc
    REAL(KIND=dp) :: Relaxation,Alpha,Beta,Gamma
    REAL(KIND=dp), ALLOCATABLE :: Diag(:), TempVector(:)
    REAL(KIND=dp), POINTER :: bb(:),Res(:)
    REAL(KIND=dp) :: t0,rt0,rst,st,ct
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NMode, LinModes
    CHARACTER(*), PARAMETER :: Caller = 'SolveSystem'
    
    INTERFACE
      SUBROUTINE BlockSolveExt(A,x,b,Solver)
        USE Types
        TYPE(Matrix_t), POINTER :: A
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp) :: x(:), b(:)
      END SUBROUTINE BlockSolveExt
    END INTERFACE   


!------------------------------------------------------------------------------
    Params => Solver % Values

    CALL Info(Caller,'Solving linear system',Level=10)

    Timing = ListCheckPrefix(Params,'Linear System Timing')
    IF( Timing ) THEN
      t0 = CPUTime(); rt0 = RealTime()
    END IF

    n = A % NumberOfRows

    IF(A % FORMAT < 1 .OR. A % FORMAT > 3 ) THEN
      CALL Fatal( Caller,'Not implemented for matrix format: '//I2S(A % format))
    END IF
      
    
    RestrictionMode = HaveRestrictionMatrix( A ) 

    ResidualMode = ListGetLogical( Params,'Linear System Residual Mode',Found )      
      
    BlockMode = ListGetLogical( Params,'Linear System Block Mode',Found ) 
    
!------------------------------------------------------------------------------
! The allocation of previous values has to be here in order to 
! work properly with the Dirichlet elimination.
!------------------------------------------------------------------------------
    NeedPrevSol = ResidualMode
    
    IF(.NOT. NeedPrevSol ) THEN
      Relaxation = ListGetCReal( Params, &
          'Nonlinear System Relaxation Factor', Found )
      IF( Found ) NeedPrevSol = (Relaxation /= 1.0_dp)
    END IF

    IF(.NOT. NeedPrevSol ) THEN
      Method = ListGetString( Params, &
        'Nonlinear System Convergence Measure', Found ) 
      NeedPrevSol = ( Method == 'residual' .OR. Method == 'solution' )
    END IF

    IF( NeedPrevSol ) THEN
      CALL Info(Caller,'Previous solution must be stored before system is solved',Level=10)
      Found = ASSOCIATED(Solver % Variable % NonlinValues)
      IF( Found ) THEN
        IF ( SIZE(Solver % Variable % NonlinValues) /= n) THEN
          DEALLOCATE(Solver % Variable % NonlinValues)
          Found = .FALSE.
        END IF
      END IF
      IF(.NOT. Found) THEN
        ALLOCATE( Solver % Variable % NonlinValues(n), STAT=istat ) 
        IF ( istat /= 0 ) CALL Fatal( Caller, 'Memory allocation error.' )
      END IF
      Solver % Variable % NonlinValues = x(1:n)
    END IF

    IF ( C_ASSOCIATED(Solver % LinBeforeProc) ) THEN
      CALL Info(Caller,'Calling procedure before solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinBeforeProc,CurrentModel,Solver, &
                       A, b, x, n, DOFs, Norm )
       IF ( istat /= 0 ) GOTO 10
    END IF

    ! If residual mode is requested make change of variables:
    ! Ax=b -> Adx = b-Ax0 = r
    IF( ResidualMode ) THEN
      CALL Info(Caller,'Changing the equation to residual based mode',Level=10)
      ALLOCATE( Res(n) ) 

      ! If needed move the current solution to N-T coordinate system
      ! before computing the residual.
      IF (ASSOCIATED(Solver % Variable % Perm)) &
          CALL RotateNTSystemAll(x, Solver % Variable % Perm, DOFs)

      CALL LinearSystemResidual( A, b, x, res )
      bb => res
      ! Set the initial guess for the residual system to zero
      x = 0.0_dp
    ELSE
      bb => b
    END IF

    FirstLoop = .TRUE.
    Nmode = 0

    BLOCK
      LOGICAL :: LFact, FreeFact, ConstraintMatrixConstant
      ! Hoisted out of the constraint-elimination branch below: nested BLOCKs
      ! are miscompiled by some compilers (Intel 20.0), and the extra scope
      ! bought nothing here.
      TYPE(Matrix_t), POINTER :: Acoll

      ConstraintMatrixConstant = ListGetLogical( Solver % Values,  &
          'Constraint Modes Constant Matrix', Found)

      IF(ConstraintMatrixConstant) THEN
        LFact = ListGetLogical( Solver % Values, 'Linear System Refactorize', Found )
        IF(.NOT. Found ) LFact = .TRUE.
      END IF

20    CONTINUE
 
      CALL ConstraintModesDriver( A, x, b, Solver, .TRUE., Nmode, LinModes, FirstLoop = FirstLoop )  

      IF ( LinModes > 0 .AND. ConstraintMatrixConstant ) THEN
        FreeFact = ListGetLogical( Solver % Values, 'Linear System Free Factorization', Found )
        IF (.NOT. Found ) FreeFact = .TRUE.
        CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', .FALSE. )
      END IF
    
      IF( BlockMode ) THEN
        CALL Info(Caller,'Solving linear system with block strategy',Level=10)
        ! Here activate constraint solve only if constraints are not treated as blocks
        IF( RestrictionMode .AND. &
            ListGetLogical( Params, 'Eliminate Linear Constraints', Found) ) THEN
          Acoll  => AllocateMatrix()
          Acoll % FORMAT = MATRIX_LIST
          CALL Info(Caller,'Eliminating constraints before going into block matrix!')
          CALL EliminateLinearRestriction( A, bb, A % ConstraintMatrix, Acoll, Solver, .TRUE. )
          CALL List_ToCRSMatrix(Acoll)

          Acoll % Comm = A % Comm
          Acoll % AddMatrix => A % AddMatrix
          CALL ParallelInitMatrix(Solver, Acoll)

          CALL BlockSolveExt( Acoll, x, Acoll % rhs, Solver )

          CALL Info(Caller,'Freeing collection matrix after solution',Level=10)
          NULLIFY( Acoll % AddMatrix )

          CALL FreeMatrix(Acoll)
          ! The environment of the freed collection matrix is gone with it,
          ! return to the one of the system matrix.
          CALL SetMatrixParEnv( A )

          Acoll => NULL()
        ELSE
          CALL BlockSolveExt( A, x, bb, Solver )
        END IF
      ELSE IF ( RestrictionMode ) THEN
        CALL Info(Caller,'Solving linear system with linear restrictions!',Level=10)
        IF( ListGetLogical( Params,'Save Constraint Matrix',Found ) ) THEN
          GloNum = ListGetLogical( Params,'Save Constraint Matrix Global Numbering',Found )
          CALL SaveProjector(A % ConstraintMatrix,.TRUE.,'cm',Parallel=GloNum)
        END IF
        CALL SolveWithLinearRestriction( A,bb,x,Norm,DOFs,Solver )
      ELSE ! standard mode
        CALL Info(Caller,'Solving linear system in standard way',Level=12)
        CALL SolveLinearSystem( A,bb,x,Norm,DOFs,Solver )
      END IF
      CALL Info(Caller,'System solved',Level=12)
    
      IF( LinModes > 0 .OR. Nmode > 0 ) THEN
        CALL ConstraintModesDriver( A, x, b, Solver, .FALSE., FirstLoop = FirstLoop ) 

        IF (ConstraintMatrixConstant) THEN
          CALL ListAddLogical( Solver % Values, 'Linear System Constant Matrix', .TRUE.)
          CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
        END IF

        FirstLoop = .FALSE.
        IF( Nmode < LinModes ) GOTO 20

        IF ( ConstraintMatrixConstant ) THEN
          CALL ListAddLogical( Solver % Values, 'Linear System Constant Matrix', .FALSE.)
          CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', LFact )
          CALL ListAddLogical( Solver % Values, 'Linear System Free Factorization', FreeFact )
        END IF
      END IF         
    END BLOCK
    
    ! Even in the residual mode the system is reverted back to complete vectors 
    ! and we may forget about the residual.
    IF( ResidualMode ) DEALLOCATE( Res ) 

!------------------------------------------------------------------------------

10  CONTINUE

    IF ( C_ASSOCIATED(Solver % LinAfterProc) ) THEN
      CALL Info(Caller,'Calling procedure after solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinAfterProc, CurrentModel, Solver, &
              A, b, x, n, DOFs, Norm )
    END IF

    IF ( Solver % TimeOrder == 2 ) THEN
      CALL Info(Caller,'Setting up PrevValues for 2nd order transient equations',Level=12)

      IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        CALL Update2ndOrder(n,Solver % dt,x, &
            Solver % Variable % PrevValues,Solver % Alpha,Solver % Beta)
      END IF
    END IF

    IF( Timing ) THEN
      st  = CPUTime() - t0;
      rst = RealTime() - rt0

      WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time (CPU,REAL) for '&
          //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
      CALL Info(Caller,Message,Level=4)    
      
      IF( ListGetLogical(Params,'Linear System Timing',Found)) THEN
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF
      
      IF( ListGetLogical(Params,'Linear System Timing Cumulative',Found)) THEN
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),Found)
        st = st + ct
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),Found)
        rst = rst + ct

        WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time cumulative (CPU,REAL) for '&
            //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
        CALL Info(Caller,Message,Level=7)    
        
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF

    END IF

    CALL Info(Caller,'Finished solving the system',Level=12)

!------------------------------------------------------------------------------
END SUBROUTINE SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve a linear eigen system.
!------------------------------------------------------------------------------
SUBROUTINE SolveEigenSystem( StiffMatrix, NOFEigen, &
    EigenValues, EigenVectors,Solver )
!------------------------------------------------------------------------------
    USE EigenSolve
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), TARGET :: EigenValues(:),EigenVectors(:,:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: StiffMatrix
    INTEGER :: NOFEigen
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    COMPLEX(KIND=dp), POINTER :: dvecs(:,:), evecs(:,:)
    REAL(KIND=dp), POINTER :: p(:)
    INTEGER :: i,j,k,n, AllocStat
    TYPE(Matrix_t), POINTER :: A
    TYPE(Variable_t), POINTER :: Var 
    LOGICAL :: Damped, Direct, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    !------------------------------------------------------------------------------
    n = StiffMatrix % NumberOfRows
    EVecs => EigenVectors

    A => StiffMatrix

    Damped = ListGetLogical( Solver % Values, 'Eigen System Damped', Found)
    Direct = ListGetString( Solver % Values, 'Linear System Solver', Found) == 'direct'

    IF( Damped .AND. Direct ) THEN
      A => GenerateStateEquationSystem(Solver,A,n)
      ALLOCATE(EVecs(NOFEigen,n)); EVecs=0
      CALL ListAddLogical( Solver % Values,'Eigen System Damped', .FALSE. )
    END IF

    IF ( .NOT. A % Complex ) THEN
      CALL Info('SolveEigenSystem','Solving real valued eigen system of size: '//I2S(n),Level=8)

      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolve( Solver, A, n, NOFEigen, EigenValues, evecs )
      ELSE
        CALL ParallelArpackEigenSolve( Solver, A, n, NOFEigen, EigenValues, evecs )
      END IF
    ELSE
      CALL Info('SolveEigenSystem','Solving complex valued eigen system of size: '//I2S(n/2),Level=8)
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolveComplex( Solver, A, n/2, NOFEigen, EigenValues, evecs )
      ELSE
        CALL ParallelArpackEigenSolveComplex( Solver, A, n/2, NOFEigen, EigenValues, EVecs )
      END IF
    END IF

    IF (Damped.AND.Direct) THEN
      n = n/2
      DO i=1,NOFEigen
        DO j=1,n
          EigenVectors(i,j) = EVecs(i,j+n)
        END DO
      END DO
      DEALLOCATE(EVecs)
      CALL ListAddLogical( Solver % Values,'Eigen System Damped', .TRUE. )
    END IF

 CONTAINS

     FUNCTION GenerateStateEquationSystem(Solver,A,n) RESULT(B)

       TYPE(Matrix_t), POINTER :: A,B
       TYPE(Solver_t), TARGET:: Solver

       INTEGER :: i,j,k,l,n, max_gdofs

! Make the system matrices for the damped eigensystem

!      (M 0)           lambda *  ( C K)
!      (0 M) [y x]^T =           (-M 0) [y x]^T

       B => AllocateMatrix()
       B % Format = MATRIX_CRS
       B % Complex = A % Complex

       B % NumberOfRows = 2*n
       l = 4*SIZE(A % Cols)
       ALLOCATE(B % Cols(l), B % Rows(2*n+1), B % Values(l), B % MassValues(l))
       B % Values = 0; B % MassValues = 0

       B % Rows(1) = 1
       l = 0
       DO i=1,n
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)
           B % Values(l) = A % DampValues(k)
           B % MassValues(l) = A % MassValues(k)
         END DO
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)+n
           B % Values(l) = A % Values(k)
         END DO
         B % Rows(i+1) = l+1
       END DO

       DO i=1,n
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)
           B % Values(l) = -A % MassValues(k)
         END DO
 
         DO k=A % Rows(i), A % Rows(i+1)-1
           l = l + 1
           B % Cols(l) = A % Cols(k)+n
           B % MassValues(l) = A % MassValues(k)
         END DO
         B % Rows(n+i+1) = l+1
       END DO

       IF(ParEnv % PEs>1) THEN
         IF(.NOT.ASSOCIATED(A % ParMatrix)) CALL ParallelInitMatrix(Solver,A)

         ALLOCATE(B % ParallelInfo)
         ALLOCATE(B % Perm(2*n))
         ALLOCATE(B % ParallelInfo % GInterface(2*n))
         ALLOCATE(B % ParallelInfo % GlobalDOFs(2*n))
         ALLOCATE(B % ParallelInfo % NeighbourList(2*n))
 
         max_gdofs = ParallelReduction( MAXVAL(A % ParallelInfo % GlobalDOFs), 2 )
         DO i=1,n
           B % ParallelInfo % NeighbourList(i) % Neighbours => A % ParallelInfo % NeighbourList(i) % Neighbours
           B % ParallelInfo % NeighbourList(i+n) % Neighbours => A % ParallelInfo % NeighbourList(i) % Neighbours

           B % ParallelInfo % Ginterface(i) = A % ParallelInfo % GInterface(i)
           B % ParallelInfo % Ginterface(i+n) = A % ParallelInfo % GInterface(i)

           B % ParallelInfo % GlobalDOFs(i) = A % ParallelInfo % GlobalDOFs(i)
           B % ParallelInfo % GlobalDOFs(i+n) = A % ParallelInfo % GlobalDOFs(i)+max_gdofs

           B % Perm(i) = A % Perm(i)
           B % Perm(i+n) = A % Perm(i)+n
         END DO
         B % Solver => Solver
         B % Parmatrix => ParInitMatrix(B, B % ParallelInfo)
       END IF
       n = 2*n
     END FUNCTION GenerateStateEquationSystem
    
!------------------------------------------------------------------------------
END SUBROUTINE SolveEigenSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Compute lumped fluxes, for example for capacitance or impedance matrices. 
!------------------------------------------------------------------------------
SUBROUTINE StoreLumpedFluxes( Solver, NoModes, iMode, FluxesRow, FluxesRowIm, &
    FluxesRhs, FluxesRhsIm, ImpRe, ImpIm ) 
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  INTEGER :: NoModes, iMode
  REAL(KIND=dp) :: FluxesRow(:)
  REAL(KIND=dp), OPTIONAL :: FluxesRowIm(:)
  REAL(KIND=dp), OPTIONAL :: FluxesRhs, FluxesRhsIm
  REAL(KIND=dp), OPTIONAL :: ImpRe, ImpIm
!------------------------------------------------------------------------------
  REAL(KIND=dp), POINTER :: FluxesMatrix(:,:), FluxesMatrixIm(:,:)
  INTEGER :: i,j,k,n,Nmode,Mmode
  REAL(KIND=dp), POINTER :: PValues(:)
  TYPE(LumpedModel_t), POINTER :: Lumped
  CHARACTER(*), PARAMETER :: Caller = 'StoreLumpedFluxes'
      
  CALL Info(Caller,'Storing lumped fluxes',Level=10)

  Lumped => Solver % Lumped
  IF(.NOT. ASSOCIATED( Lumped ) ) THEN
    CALL Info(Caller,'Allocating lumped model of size: '//I2S(NoModes), Level=5)
    ALLOCATE( Solver % Lumped )
    Lumped => Solver % Lumped

    Lumped % NoModes = NoModes
    ALLOCATE( Lumped % CMatrix( NoModes, NoModes ) )
    Lumped % CMatrix = 0.0_dp

    IF( PRESENT(FluxesRowIm) ) THEN
      CALL Info(Caller,'Storing lumped fluxes imaginary',Level=20)
      Lumped % IsComplex = .TRUE.
      ALLOCATE( Lumped % CMatrixIm( NoModes, NoModes ) )
      Lumped % CMatrixIm = 0.0_dp
    END IF
    IF( PRESENT(FluxesRhs) ) THEN
      ALLOCATE( Lumped % CRhs( NoModes ) )
      Lumped % Crhs = 0.0_dp
    END IF
    IF( PRESENT(FluxesRhsIm) ) THEN
      Lumped % IsComplex = .TRUE.
      ALLOCATE( Lumped % CRhsIm( NoModes ) )
      Lumped % CrhsIm = 0.0_dp
    END IF
    IF( PRESENT(ImpRe) ) THEN
      ALLOCATE( Lumped % ImpRe( NoModes ) )
      Lumped % ImpRe = 0.0_dp
    END IF
    IF( PRESENT(ImpIm) ) THEN
      ALLOCATE( Lumped % ImpIm( NoModes ) )
      Lumped % ImpIm = 0.0_dp
    END IF
  END IF

  IF( Lumped % IsComplex ) FluxesMatrixIm => Solver % Lumped % CMatrixIm

  Lumped % CMatrix(iMode,1:NoModes) = FluxesRow(1:NoModes)
  IF(PRESENT(FluxesRowIm) ) THEN
    Lumped % CMatrixIm(iMode,1:NoModes) = FluxesRowIm(1:NoModes)
  END IF
  IF(PRESENT(FluxesRhs) ) THEN
    Lumped % CRhs(iMode) = FluxesRhs
  END IF
  IF(PRESENT(FluxesRhsIm) ) THEN
    Lumped % CRhsIm(iMode) = FluxesRhsIm
  END IF
  IF( PRESENT(ImpRe) ) THEN
    Lumped % ImpRe(iMode) = ImpRe
  END IF
  IF( PRESENT(ImpIm) ) THEN
    Lumped % ImpIm(iMode) = ImpIm
  END IF
        
  Lumped % CntModes = iMode

  BLOCK
    INTEGER :: ADepth
    IF( ListGetLogicalAnySolver( CurrentModel,'Adaptive Mesh Refinement') ) THEN      
      Adepth = CurrentModel % Mesh % AdaptiveDepth
      WRITE(Message,*) 'Row'//I2S(iMode)//' Depth'//I2S(Adepth)//':',FluxesRow(1:NoModes)
      CALL Info('StoreLumpedFluxes',Message,Level=4)
    END IF
  END BLOCK
  
END SUBROUTINE StoreLumpedFluxes


SUBROUTINE FinalizeLumpedMatrix( Solver )

  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
  REAL(KIND=dp), POINTER :: FluxesMatrix(:,:), FluxesMatrixIm(:,:)
  INTEGER :: i,j,k,NoModes
  LOGICAL :: Symmetric, IsComplex, EmWaveMode, CoilMode, Found, DoIt
  REAL(KIND=dp) :: nrm
  CHARACTER(:), ALLOCATABLE :: MatrixFile
  CHARACTER(*), PARAMETER :: Caller = 'FinalizeLumpedMatrix'
  COMPLEX(KIND=dp) :: cx, cb
  REAL(KIND=dp), ALLOCATABLE :: Checksum(:)
  TYPE(LumpedModel_t), POINTER :: Lumped

  CALL Info(Caller,'Finalizing lumped matrix',Level=8)
  
  IF(.NOT. ASSOCIATED( Solver % Lumped ) ) THEN
    CALL Fatal(Caller,'We should not be here without the lumped system!')
  END IF

  EmWaveMode = ListGetLogical( Solver % Values,'Constraint Modes EM Wave', Found )
  CoilMode = ListGetLogical( Solver % Values,'Constraint Modes Coils', Found )
  
  Lumped => Solver % Lumped  
  NoModes = Lumped % NoModes
  IF( Lumped % CntModes /= NoModes ) THEN
    CALL Fatal(Caller,'Trying to deduce '//I2S(NoMOdes)//&
        ' rows with '//I2S(Lumped % CntModes)//' lines of data!')
  END IF
  
  IsComplex = Lumped % IsComplex

  FluxesMatrix => Lumped % CMatrix
  IF( IsComplex ) FluxesMatrixIm => Lumped % CMatrixIm

  ! Energies and Impedances are related by factor two when currents are "1".   
  IF(CoilMode) THEN
    BLOCK
      REAL(KIND=dp) :: DesiredCurr(NoModes)  
      ! Normalize the inductance in case it has been computed with non-unity currents!
      DO i=1,CurrentModel % NumberOfComponents
        j = ListGetInteger( CurrentModel % Components(i) % Values,'Constraint Mode',Found )
        IF(j > 0) THEN
          DesiredCurr(j) = ListGetCReal( CurrentModel % Components(i) % Values,'Desired Coil Current',Found )
          IF(.NOT. Found) DesiredCurr(j) = 1.0_dp
        END IF
      END DO
      DO i=1,NoModes
        DO j=1,NoModes
          FluxesMatrix(i,j) = FluxesMatrix(i,j) / (DesiredCurr(i)*DesiredCurr(j)) 
        END DO
      END DO
    END BLOCK
  END IF

  ! Normalize the S-parameter matrix
  IF(EmWaveMode) THEN    
    ! Normalize by the source    
    BLOCK
      
      LOGICAL :: FixIt, NoNormalize
      
      FixIt =  ListGetLogical( Solver % Values,'Enforce Unity rowsum',Found )
      NoNormalize = ListGetLogical( Solver % values, 'Skip Normalize fluxes', Found )

      IF( InfoActive(20) ) THEN        
        CALL Info( Caller,'Showing matrix before normalization!')
        DO i=1,NoModes
          DO j=1,NoModes
            WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,FluxesMatrix(i,j),FluxesMatrixIm(i,j)
            CALL Info( Caller, Message )
          END DO
        END DO
        DO i=1,NoModes
          WRITE( Message,*) 'Normalization vector '//I2S(i)//':',Lumped % Crhs(i)
          CALL Info( Caller, Message )
        END DO
      END IF      
   
      IF (.NOT. NoNormalize) THEN
        DO i=1,NoModes
          DO j=1,NoModes         
            nrm = SQRT(Lumped % Crhs(j) * Lumped % Crhs(i))                               
            FluxesMatrix(i,j) = FluxesMatrix(i,j) / nrm
            FluxesMatrixIm(i,j) = FluxesMatrixIm(i,j) / nrm
          END DO
        END DO
      END IF
        
      IF( FixIt ) THEN
        DO i=1,NoModes
          nrm = 2*FluxesMatrix(i,i) / SUM(FluxesMatrix(i,:)**2 + FluxesMatrixIm(i,:)**2) 
          FluxesMatrix(i,:) = FluxesMatrix(i,:) * nrm
          FluxesMatrixIm(i,:) = FluxesMatrixIm(i,:) * nrm
          
          WRITE(Message,*) 'Fixing Multiplier '//I2S(i)//':',nrm
          CALL Info( Caller, Message, Level=6)          
        END DO
      END IF

      DO i=1,NoModes
        FluxesMatrix(i,i) = FluxesMatrix(i,i) - 1.0_dp
      END DO
    END BLOCK    
  END IF
  
  Symmetric = ListGetLogical( Solver % Values,&
      'Constraint Modes Fluxes Symmetric', Found ) 
  IF(.NOT. Found) Symmetric = ListGetLogical( Solver % Values,&
      'Constraint Modes Matrix Symmetric', Found ) 
  IF( Symmetric ) THEN        
    CALL Info(Caller,'Enforcing symmetry of reduced system!',Level=8)
    
    IF( InfoActive(10) ) THEN
      CALL Info( Caller,'Showing asymmetry of matrix before enforced symmetry!')
      DO i=1,NoModes
        DO j=i+1,NoModes
          IF( IsComplex ) THEN
            WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,&
                FluxesMatrix(i,j)-FluxesMatrix(j,i),FluxesMatrixIm(i,j)-FluxesMatrixIm(j,i)
          ELSE                
            WRITE( Message, '(I3,I3,ES17.9)' ) i,j,FluxesMatrix(i,j)-FluxesMatrix(j,i)
          END IF
          CALL Info( Caller, Message )
        END DO
      END DO
    END IF
    
    FluxesMatrix = 0.5_dp * ( FluxesMatrix + TRANSPOSE( FluxesMatrix ) )
    IF( IsComplex ) THEN
      FluxesMatrixIm = 0.5_dp * ( FluxesMatrixIm + TRANSPOSE( FluxesMatrixIm ) )
    END IF
  END IF
  
  CALL Info( Caller,'Lumped Matrix', Level=5 )
  DO i=1,NoModes
    DO j=1,NoModes
      IF( Symmetric .AND. j < i ) CYCLE
      IF( IsComplex ) THEN
        WRITE( Message, '(I3,I3,2ES17.9)' ) i,j,FluxesMatrix(i,j),FluxesMatrixIm(i,j)
      ELSE
        WRITE( Message, '(I3,I3,ES17.9)' ) i,j,FluxesMatrix(i,j)
      END IF
      CALL Info( Caller, Message, Level=5 )
    END DO
  END DO

  MatrixFile = ListGetString(Solver % Values,'Constraint Modes Fluxes Filename',Found )
  IF(.NOT. Found) MatrixFile = ListGetString(Solver % Values,'Constraint Modes Matrix Filename',Found )
  IF( Found ) THEN
    ! Find the lowest active partition working here.
    k = ParallelReduction(ParEnv % MyPe,1)     
    IF( k == ParEnv % MyPe ) THEN
      OPEN(10, FILE=MatrixFile)
      DO i=1,NoModes
        WRITE (10,*) FluxesMatrix(i,:)
      END DO
      CLOSE(10)

      IF( IsComplex ) THEN
        OPEN( 11, FILE=TRIM(MatrixFile)//'_im')
        DO i=1,NoModes
          WRITE (11,*) FluxesMatrixIm(i,:) 
        END DO
        CLOSE(11)
        
        ALLOCATE(CheckSum(NoModes))
        OPEN( 11, FILE=TRIM(MatrixFile)//'_abs')
        DO i=1,NoModes
          WRITE (11,*) SQRT(FluxesMatrix(i,:)**2+FluxesMatrixIm(i,:)**2)           
          CheckSum(i) = SUM(FluxesMatrix(i,:)**2+FluxesMatrixIm(i,:)**2) 
        END DO
        CLOSE(11)

        OPEN( 11, FILE=TRIM(MatrixFile)//'_angle')
        DO i=1,NoModes
          WRITE (11,*) ( 180.0_dp / PI ) * ATAN2(FluxesMatrixIm(i,:),FluxesMatrix(i,:))           
        END DO
        CLOSE(11)

        IF( ASSOCIATED( Lumped % ImpRe ) ) THEN
          OPEN( 11, FILE=TRIM(MatrixFile)//'_Z')
          DO i=1,NoModes
            WRITE (11,*) Lumped % ImpRe(i), Lumped % ImpIm(i) 
          END DO
          CLOSE(11)
        END IF

        WRITE(Message,*) 'Normalization checksum: ',CheckSum
        CALL Info(Caller,Message) 
      END IF
      CALL Info( Caller,'Constraint modes fluxes was saved to file '//TRIM(MatrixFile),Level=5)
    END IF
  END IF
  
  nrm = 0.0_dp
  DO i=1,NoModes
    DO j=1,NoModes
      nrm = nrm + ABS(FluxesMatrix(i,j))
      IF(IsComplex) nrm = nrm + ABS(FluxesMatrixIm(i,j))
    END DO
  END DO

  WRITE(Message,'(A,ES12.5)') TRIM(Solver % Variable % Name)//' lumped matrix norm: ',Nrm
  CALL Info(Caller, Message, Level=6)

  DoIt = ListGetLogical( Solver % Values,'Constraint Modes Fluxes Norm', Found )
  IF(.NOT. Found) DoIt = ListGetLogical( Solver % Values,'Constraint Modes Matrix Norm', Found )
  IF(DoIt) CALL ListAddConstReal( CurrentModel % Simulation,&
      'res: '//TRIM(Solver % Variable % Name)//' lumped matrix norm',nrm)
  
  IF( ListGetLogical( Solver % Values,'Constraint Modes Matrix Results', Found ) ) THEN
    CALL Info(Caller,'Adding Constraint Modes Fluxes with "res:" to list',Level=5)
    DO i=1,NoModes
      DO j=1,NoModes
        CALL ListAddConstReal( CurrentModel % Simulation,'res: CMF '//I2S(i)//' '//I2S(j),FluxesMatrix(i,j))
      END DO
    END DO
    IF( IsComplex ) THEN
      DO i=1,NoModes
        DO j=1,NoModes
          CALL ListAddConstReal( CurrentModel % Simulation,'res: CMF Im '//I2S(i)//' '//I2S(j),FluxesMatrixIm(i,j))
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE FinalizeLumpedMatrix


! Given an edge element and associated field variable computes the circulation over the
! element. The consistency of the circulation is checked 1st ensuring the each edge is
! accounted for correctly and the element is consistently accounted for with the reference
! normal.
!------------------------------------------------------------------------------------------
SUBROUTINE EdgeElementCirculation( Circ, Mesh, Element, Nodes, Solver, RefNormal, Var )
  REAL(KIND=dp) :: Circ(:)
  TYPE(Mesh_t) :: Mesh
  TYPE(Element_t) :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(Solver_t) :: Solver
  TYPE(Variable_t), POINTER, OPTIONAL :: Var
  REAL(KIND=dp), OPTIONAL :: RefNormal(3)

  REAL(KIND=dp) :: EdgeVector(3), EdgeVector1(3), Normal(3), DotProd, s
  INTEGER :: dofs, dofIndexes(20), i, j, k, n, i1, i2, j1, j2, nb
  INTEGER :: ownParent, totParent
  TYPE(Variable_t), POINTER :: pVar

  Circ = 0.0_dp
  nb = mGetElementDOFs( dofIndexes, Element, Solver )
  n = Element % Type % NumberOfEdges
  
  IF(PRESENT(Var)) THEN
    pVar => Var
  ELSE
    pVar => Solver % Variable 
  END IF

  IF( ParEnv % PEs > 1 ) THEN
    TotParent = 0
    OwnParent = 0   
    IF(ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
      TotParent = TotParent + 1
      IF(Element % Boundaryinfo % Left % PartIndex == ParEnv % MyPe ) ownParent = ownParent + 1
    END IF
    IF(ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
      TotParent = TotParent + 1
      IF(Element % Boundaryinfo % Right % PartIndex == ParEnv % MyPe ) ownParent = ownParent + 1
    END IF
    IF(ownParent == 0) RETURN
  END IF
  
  DO i=1,nb
    i1 = i
    i2 = MODULO(i,n)+1

    ! Vector in the direction of the edge
    EdgeVector(1) = Nodes % x(i2) - Nodes % x(i1)
    EdgeVector(2) = Nodes % y(i2) - Nodes % y(i1)
    EdgeVector(3) = Nodes % z(i2) - Nodes % z(i1)

    ! Among different elements define a consistent direction for circulation
    ! Store 1st tangent vector and compute the normal from the cross product. 
    IF( PRESENT(RefNormal) ) THEN
      IF(i==1) THEN
        Edgevector1 = -EdgeVector
      ELSE IF(i==2) THEN
        Normal = CrossProduct( EdgeVector1, EdgeVector )
        DotProd = SUM(Normal * RefNormal) 
      END IF
    END IF

    ! Among edges for one element define consistent direction for circulation
    ! Edges are (1,2), (2,3), (3,1) for triangle and likewise for quad.
    j1 = Element % NodeIndexes(i1)
    j2 = Element % NodeIndexes(i2) 
    IF( ParEnv % PEs > 1 ) THEN                            
      j1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
      j2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
    END IF

    ! Integration length and direction
    s = SQRT(SUM(EdgeVector**2))
    IF( j1 < j2) s = -s

    j = pVar % Perm(DofIndexes(i))
    DO k=1,pVar % Dofs
      Circ(k) = Circ(k) + s * pVar % Values(pVar % Dofs*(j-1)+k)
    END DO
  END DO

  ! If the element normal is pointing to different direction reverse the flux.
  IF( PRESENT(RefNormal) ) THEN
    IF( DotProd < 0.0_dp ) Circ = -Circ
  END IF

  ! If we are at interface of partitions the face element may appear twice. Then take just half. 
  IF( ParEnv % PEs > 1 ) THEN
    Circ = Circ * ownParent / totParent    
  END IF
  
END SUBROUTINE EdgeElementCirculation


SUBROUTINE BoundaryCirculation( BCCirc, BcInd, Solver, Var )
  REAL(KIND=dp) :: BCCirc(:)
  INTEGER :: BcInd
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Variable_t), POINTER :: Var
  
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp) :: Normal0(3), ElemCirc(3), parCoeff
  INTEGER :: i,t,k,n,nb,elem,dofs,ownParent,TotParent
  LOGICAL :: NormalSet
  TYPE(Solver_t), POINTER :: pSolver
  
  Mesh => Solver % Mesh
  IF( Mesh % MeshDim < 3 ) THEN
    CALL Fatal('BoundaryCirculation','Currently only available in 3D!')
  END IF
  pSolver => Solver
  dofs = Var % dofs
  BCCirc = 0.0_dp
  NormalSet = .FALSE.
  
  DO elem=Mesh % NumberOfBulkElements + 1, &
      Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        
    
    Element => Mesh % Elements(elem)    

    IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
    IF(Element % BoundaryInfo % Constraint /= CurrentModel % BCs(bcInd) % Tag) CYCLE
    
    n = Element % Type % NumberOfNodes
    CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)

    ! Compute a reference normal to which all other face normals are compared against.
    IF( .NOT. NormalSet ) THEN
      Normal0 = NormalVector(Element,Nodes)
      NormalSet = .TRUE.
    END IF
      
    CALL EdgeElementCirculation( ElemCirc, Mesh, Element, Nodes, pSolver, Normal0, Var )     
    BcCirc(1:dofs) = BcCirc(1:dofs) + ElemCirc(1:dofs)
  END DO

  IF( ParEnv % PEs > 1 ) THEN
    DO i=1,dofs
      BcCirc(i) = ParallelReduction(BCCirc(i))
    END DO
  END IF
    
END SUBROUTINE BoundaryCirculation



!------------------------------------------------------------------------------
!> Solve a linear system with permutated constraints.
!------------------------------------------------------------------------------
SUBROUTINE ConstraintModesDriver( A, x, b, Solver, PreSolve, ThisMode, LinSysModes, FirstLoop )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t), TARGET :: Solver
    REAL(KIND=dp) CONTIG :: x(:),b(:)
    LOGICAL :: PreSolve
    INTEGER, OPTIONAL :: ThisMode, LinSysModes
    LOGICAL, OPTIONAL :: FirstLoop
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i,j,k,n,NoModes,Mmode,ierr
    LOGICAL :: PrecRecompute, Stat, Found, ComputeFluxes, ComputeLinkage, Symmetric, &
        IsComplex, Parallel, ConsiderP, RhsMode, EmWaveMode, CoilMode, GaussLaw
    REAL(KIND=dp), ALLOCATABLE :: Fluxes(:), b0(:), A0(:), TempRhs(:)
    REAL(KIND=dp), ALLOCATABLE :: FluxesRow(:), FluxesRowIm(:)
    REAL(KIND=dp) :: FluxesRhs, FluxesRhsIm, ImpRe, ImpIm
    LOGICAL, ALLOCATABLE :: ConstrainedDOF0(:)
    REAL(KIND=dp) :: flux
    CHARACTER(:), ALLOCATABLE :: MatrixFile, BCName
    INTEGER :: NMode = 0, dof
    TYPE(Variable_t), POINTER :: pVar
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: LinsysMode, EigenMode, GotBC, LumpedMode 
    CHARACTER(*), PARAMETER :: Caller = 'ConstraintModesDriver'

    SAVE FluxesRow, FluxesRowIm, Fluxes, TempRhs, A0, b0, ConstrainedDOF0, LinsysMode, NMode

    
    !------------------------------------------------------------------------------
    IF(PRESENT(LinSysModes)) LinSysModes = 0
    NoModes = Solver % NumberOfConstraintModes 

    IF(NoModes <= 0) RETURN    
    Params => Solver % Values

    ! We can also have a combination of standard analysis + constraint modes
    ! analysis of the frozen state. Then the default solution slot should really do
    ! the standard analysis.
    IF( ListGetLogical( Params,'Constraint Modes Analysis Frozen',Found ) ) THEN
      Solver % Variable % FrozenMode = .TRUE.
      RETURN
    END IF

    EigenMode = ListgetLogical( Solver % Values, 'Eigen Analysis', Found ) 
    
    Var => Solver % Variable
    n = A % NumberOfRows        
    Parallel = Solver % Parallel
    
    IsComplex = ListGetLogical( Params,'Linear System Complex',Found)


    ! If the mode is nodal it is not lumped
    ! If it relates to whole boundary it is. 
    LumpedMode = ListGetLogical( Params,'Constraint Modes Lumped',Found )
    
    ! This is to my understanding not needed. To estimate the fluxes we
    ! basically integrate over basis functions that estimate unity.
    ! For p-elements this means using the linear nodal basis only, not any
    ! of the fake fluxes associated to p-degrees of freedom. Anyways,
    ! we leave this option here for testing etc. 
    ConsiderP = ListGetLogical( Params,'Consider P Fluxes',Found ) 

    ComputeFluxes = ListGetLogical( Params,'Constraint Modes Fluxes',Found) 
    ComputeLinkage = ListGetLogical( Params,'Constraint Modes Linkage',Found )
    RhsMode = ListGetLogical( Params,'Constraint Modes rhs',Found )
    EmWaveMode = ListGetLogical( Params,'Constraint Modes EM Wave', Found ) 
    CoilMode = ListGetLogical( Params,'Constraint Modes Coils',Found ) 

    IF(EmWaveMode) IsComplex = .TRUE.
    GaussLaw = .FALSE.
    IF(EMWaveMode) GaussLaw = ListGetLogical( Params,'Use Gauss Law',Found ) 
    
    IF( EmWaveMode .OR. CoilMode ) THEN
      ComputeFluxes = .TRUE.
      RhsMode = .TRUE.
    END IF

    
    IF( PreSolve ) THEN
      CALL Info(Caller,'Number of constraint modes is: '//I2S(NoModes),Level=8)
          
      ! We loop over the mode if it is not given in some external loop.
      !---------------------------------------------------------------------
      pVar => NULL()
      IF( ListGetLogical( Params,'Nonlinear System Constraint Modes', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'nonlin iter')    
      IF( ListGetLogical( Params,'Steady State Constraint Modes', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'coupled iter')
      IF( ListGetLogical( Params,'Run Control Constraint Modes', Found ) .OR. &
          ListGetLogical( CurrentModel % Control,'Constraint Modes Analysis', Found ) ) &
          pVar => VariableGet( Solver % Mesh % Variables,'run')    
      LinSysMode = .NOT. ASSOCIATED(pVar)
      
      IF(LinSysMode) THEN
        LinSysModes = NoModes
        ! If we combined eigen analysis base + constraint modes base do an empty
        ! cycle that does the eigenmodes first without updating the counter. 
        IF( PRESENT( FirstLoop ) ) THEN
          IF( FirstLoop .AND. EigenMode ) RETURN
        END IF        
        Nmode = ThisMode + 1
      ELSE
        Nmode = NINT( pVar % Values(1) ) 
        LinSysModes = 0
      END IF
      ThisMode = Nmode
      
      IF(EigenMode) CALL ListAddLogical( Solver % Values,'Skip Eigen Analysis',.TRUE.)

      IF( SIZE(x) /= n ) THEN
        CALL Fatal(Caller,'Conflicting sizes for matrix and variable! ('//I2S(SIZE(x))//','//I2S(n)//')')
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        ALLOCATE( FluxesRow(NoModes), Fluxes(n) )       
        IF( IsComplex ) ALLOCATE( FluxesRowIm(NoModes) )         
        
        IF( Parallel ) THEN
          ALLOCATE(TempRHS(SIZE(A % BulkRhs)))
          TempRhs = 0.0_dp
        END IF
        
        IF( IsComplex .OR. CoilMode) THEN
          ALLOCATE( A0(SIZE(A % Values)), b0(n), ConstrainedDOF0(n) )
          A0 = A % Values
          b0 = A % Rhs
          ConstrainedDOF0 = A % ConstrainedDOF
        END IF
      END IF
      
      IF(LinSysMode .AND. NMode == 2 ) THEN
        CALL ListAddLogical( Params,'No Precondition Recompute',.TRUE.)
      END IF

      CALL Info(Caller,'Setting up constrained mode: '//I2S(NMode),Level=6)
      i = Nmode


      ! By default constraint modes are set to 0/1.
      ! However, we can also set the BC's in some other way using prefix "mode 1:" etc.  
      GotBC = .FALSE.
      IF( LumpedMode ) THEN
        DO dof=1,Var % dofs
          BcName = 'mode '//I2S(Nmode)//': '//ComponentName(Var % name,dof)
          IF(ListCheckPresentAnyBC(CurrentModel, BcName ) ) GotBC = .TRUE.
        END DO
      END IF
        
      ! The matrix has been manipulated already before. This ensures
      ! that the system has values 1 at the constraint mode i.
      IF( CoilMode ) THEN                
        b = b0
        WHERE( Var % ConstraintModesIndeces > 0 .AND. Var % ConstraintModesIndeces /= Nmode ) 
          b = 0.0_dp
        END WHERE
      ELSE IF( RhsMode ) THEN
        IF( IsComplex ) THEN       
          A % Values = A0
          WHERE( Var % ConstraintModesIndeces > 0 )
            b = 0.0_dp
          END WHERE
          WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 .OR. &
              Var % ConstraintModesIndeces == 2*Nmode )
            ! Revert rhs to previous
            b = b0
          END WHERE
        ELSE       
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            WHERE( Var % ConstraintModesIndeces == Nmode-1 ) 
              b = 0.0_dp
            END WHERE
          END IF
          WHERE( Var % ConstraintModesIndeces == Nmode ) 
            b = Var % ConstraintModesWeights
          END WHERE
        END IF
      ELSE        
        IF( IsComplex ) THEN
          ! Quick and a little dirty fix for complex capacitance matrix.
          IF( ListGetLogical( Solver % Values,'Calculate Capacitance Matrix',Found ) ) THEN
            WHERE( Var % ConstraintModesIndeces /= 0 ) 
              A % ConstrainedDOF = .TRUE.
              A % DValues = 0.0_dp
            END WHERE
            WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 )
              A % DValues = 1.0_dp
            END WHERE
          ELSE
            A % Values = A0
            b = b0
            A % ConstrainedDOF = ConstrainedDOF0

            WHERE( Var % ConstraintModesIndeces == 2*Nmode-1 )
              A % ConstrainedDOF = .TRUE.
              A % DValues = 1.0_dp
            END WHERE
            WHERE( Var % ConstraintModesIndeces == 2*Nmode ) 
              A % ConstrainedDOF = .TRUE.
              A % DValues = 0.0_dp
            END WHERE
          END IF
          CALL EnforceDirichletConditions( Solver, A, b )

        ELSE IF( GotBC ) THEN
          
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            DO dof=1,Var % dofs
              WHERE( Var % ConstraintModesIndeces == Var % Dofs*(Nmode-2)+dof ) 
                A % DValues = 0.0_dp
              END WHERE
            END DO
          END IF
          
          DO dof=1,Var % dofs
            BcName = 'mode '//I2S(Nmode)//': '//ComponentName(Var % name,dof)
            IF(ListCheckPresentAnyBC(CurrentModel, BcName ) ) THEN            
              CALL Info(Caller,"Setting constraint for: "//TRIM(BCName),Level=7)
              CALL SetDirichletBoundaries( CurrentModel, A, b, &
                  BcName, dof, Var % DOFs, Var % Perm )
            END IF
          END DO

          CALL EnforceDirichletConditions( Solver, A, b )
          
        ELSE
          
          IF( Nmode > 1 .AND. LinSysMode ) THEN
            WHERE( Var % ConstraintModesIndeces == Nmode-1 ) 
              A % DValues = 0.0_dp
            END WHERE
          END IF

          WHERE( Var % ConstraintModesIndeces == Nmode ) 
            A % DValues = 1.0_dp
          END WHERE

          CALL EnforceDirichletConditions( Solver, A, b )                    
        END IF
      END IF
      CALL ListAddLogical( Params,'Skip Zero Rhs Test',.TRUE. )
    END IF
      
    IF( .NOT. PreSolve ) THEN 
      IF( PRESENT( FirstLoop ) ) THEN
        IF( FirstLoop .AND. EigenMode ) RETURN
      END IF

      CALL Info(Caller,'Mode '//I2S(NMode)//' computed, doing some postprocessing',Level=10)

      IF( .NOT. ( IsComplex .OR. CoilMode ) ) THEN
        WHERE( Var % ConstraintModesIndeces == Nmode ) b = 0.0_dp
      END IF
      
      IF( NMode <= Var % NumberOfConstraintModes ) THEN
        Var % ConstraintModes(NMode,:) = x
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        CALL Info(Caller,'Computing lumped fluxes',Level=8)
        
        IF( CoilMode ) THEN
          CALL MagneticEnergies()
        ELSE IF(EmWaveMode ) THEN
          CALL EMWaveFluxes()
        ELSE IF( ComputeFluxes ) THEN
          CALL ConstraintModesFluxes(EmWaveMode)
        ELSE IF( ComputeLinkage ) THEN
          CALL ConstraintModesLinkage()
        END IF
        
        ! Do parallel communication here at one sweep, not before!
        IF(.NOT. EMWaveMode ) THEN
          CALL CommunicateConstraintModesFluxes()
        END IF
        
        IF( IsComplex ) THEN
          IF(EmWaveMode ) THEN
            CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow, FluxesRowIm, FluxesRhs, FluxesRhsIm, ImpRe, ImpIm)
          ELSE
            CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow, FluxesRowIm, FluxesRhs, FluxesRhsIm)
          END IF
        ELSE
          CALL StoreLumpedFluxes(Solver, NoModes, NMode, FluxesRow ) 
        END IF
      END IF
            
      IF(LinSysMode .AND. NMode == NoModes ) THEN
        IF( ComputeFluxes .OR. ComputeLinkage .OR. CoilMode ) THEN
          CALL FinalizeLumpedMatrix( Solver )            
        END IF
        CALL ListAddLogical( Params,'No Precondition Recompute',.FALSE.)
      END IF

      IF( ComputeFluxes .OR. ComputeLinkage ) THEN
        DEALLOCATE( Fluxes, FluxesRow )
        IF( IsComplex ) DEALLOCATE( FluxesRowIm )
        IF( Parallel ) DEALLOCATE(TempRHS)
        IF( IsComplex .OR. CoilMode) THEN
          A % Values = A0
          A % Rhs = b0
          A % ConstrainedDOF = ConstrainedDOF0           
          DEALLOCATE( A0, b0, ConstrainedDOF0 )
        END IF
      END IF

      CALL ListAddLogical( Params,'Skip Zero Rhs Test',.FALSE. )

      IF(EigenMode) THEN
        CALL ListAddLogical( Solver % Values,'Skip Eigen Analysis',.FALSE.)
      END IF


    END IF


    
  CONTAINS

    SUBROUTINE MagneticEnergies()
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: CVar, AVar
      TYPE(ValueList_t), POINTER :: Vlist
      INTEGER, POINTER :: MasterEntities(:)
      REAL(KIND=dp) :: Nrm, FL
      
      Mesh => Solver % Mesh
      CVar => VariableGet( Mesh % Variables,'CoilCurrent e',ThisOnly=.TRUE.)
      IF(.NOT. ASSOCIATED(CVar)) THEN
        CVar => VariableGet( Mesh % Variables,'CoilCurrent',ThisOnly=.TRUE.)
      END IF
      AVar => Solver % Variable

      FluxesRow = 0.0_dp      
      DO j=1,CurrentModel % NumberOfComponents
        Vlist => CurrentModel % Components(j) % Values       
        k = ListGetInteger( Vlist,'Constraint Mode', Found )
        IF(.NOT. Found) CYCLE

        IF( ListGetLogical( Vlist,'Flux Linkage', Found ) ) THEN
          FL = ComponentStokesTheorem(CurrentModel, Mesh, Vlist, AVar, Surf = .TRUE.)
          WRITE(Message,'(A,ES12.3)') 'Flux Linkage '//I2S(NMode)//' '//I2S(k)//':', FL
          CALL Info(Caller,Message)
        END IF

        MasterEntities => ListGetIntegerArray( Vlist,'Master Bodies',Found )
        FluxesRow(k) = ComponentCoilEnergy(CurrentModel, Mesh, MasterEntities, avar, cvar )

        WRITE(Message,'(A,ES12.3)') 'Coil Inductance '//I2S(Nmode)//' '//I2S(k)//':', FluxesRow(k)
        CALL Info(Caller,Message)
      END DO
            
      Nrm = SUM(FluxesRow)
      WRITE(Message,'(A,ES12.3)') 'Energy norm of Impedance matrix row '//I2S(NMode)//': ',Nrm
      CALL Info(Caller,Message,Level=8)

    END SUBROUTINE MagneticEnergies


    SUBROUTINE EMWaveFluxes()
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: AVar
      TYPE(ValueList_t), POINTER :: Vlist
      INTEGER, POINTER :: MasterEntities(:)
      COMPLEX(KIND=dp) :: OutFlux,InFlux,PortFlux,InImp,PortImp
      INTEGER :: i,j,k,n,port,alloc     
      
      CALL Info(Caller,'Using <Ej,Ej> for lumping',Level=10)
      
      Mesh => Solver % Mesh
      AVar => Solver % Variable      
      FluxesRow = 0.0_dp
      FluxesRowIm = 0.0_dp

      FluxesRhs = 1.0_dp 
      FluxesRhsIm = 0.0_dp

      DO port = 1, NoModes
        DO alloc = 0, 1
          n = 0
          DO j=1,CurrentModel % NumberOfBCs        
            Vlist => CurrentModel % BCs(j) % Values       
            k = ListGetInteger( Vlist,'Constraint Mode', Found )
            IF(k == port) THEN
              n = n+1
              IF(alloc == 1 ) MasterEntities(n) = j
            END IF
          END DO
          IF(alloc == 1 ) EXIT
          ALLOCATE(MasterEntities(n))
        END DO

        OutFlux = BoundaryWaveFlux(CurrentModel, Mesh, MasterEntities, Avar, PortFlux, PortImp, port==NMode )
          
        ! Memorize the coefficient for normalization: <Ec,Ej>/<Ei,Ei>                
        ! Real and imag part of: <Ec,Ej>
        FluxesRow(port) = REAL(OutFlux)
        FluxesRowIm(port) = AIMAG(OutFlux) 

        ! Memorize diagonal entry <Ej,Ej*> for future normalization
        IF(port==NMode) THEN
          InFlux = PortFlux
          InImp = PortImp
        END IF
               
        DEALLOCATE(MasterEntities)                
      END DO

      FluxesRhs = REAL(InFlux)      
      FluxesRhsIm = AIMAG(InFlux)
      ImpRe = REAL(InImp)
      ImpIm = AIMAG(InImp)
      
    END SUBROUTINE EMWaveFluxes

      
    
    SUBROUTINE ConstraintModesFluxes(EmWaveMode)
      LOGICAL :: EmWaveMode

      REAL(KIND=dp), POINTER CONTIG :: PValues(:)
      INTEGER :: poffset
      REAL(KIND=dp) :: flux, w
      COMPLEX(KIND=dp) :: cflux, cx, cmult, cb, crhs
      
      ! Use the initial bulk values that do not include Dirichlet conditions
      PValues => A % Values
      IF( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
        CALL Fatal(Caller,'BulkValues not associated!')
      END IF
      A % Values => A % BulkValues
            
      Fluxes = 0.0_dp
      IF ( Parallel ) THEN
        TempRHS = A % BulkRhs 
        CALL ParallelInitSolve( A, x, TempRHS, Fluxes )
        CALL ParallelMatrixVector( A, x, Fluxes, .TRUE. )
      ELSE
        CALL MatrixVectorMultiply( A, x, Fluxes )
      END IF

      ! Revert pointer back 
      A % Values => PValues

      poffset = 2*(NoModes + 1)
      FluxesRow = 0.0_dp
      IF(IsComplex) THEN
        FluxesRowIm = 0.0_dp
        FluxesRhs = 0.0_dp
        FluxesRhsIm = 0.0_dp
      END IF
      
      IF( EmWaveMode ) THEN
        w = ListGetAngularFrequency( Found = Found )
        IF(.NOT. Found) CALL Fatal(Caller,'Energy mode requires "Angular Frequency"!')
        cmult = 1.0/(2*w*CMPLX(0.0_dp,1.0_dp,KIND=dp)) 
      END IF
      
      DO j=1,n
        k = Var % ConstraintModesIndeces(j)
        
        IF( ConsiderP ) THEN
          ! P dofs are associated with negative index as they are not included in ConstraintModesAnalysis.
          IF( k < -1 ) k = k + poffset
        END IF

        IF( k > 0 ) THEN          
          IF( IsComplex ) THEN
            Mmode = (k+1)/2
            IF( MOD(k,2) == 1 ) THEN                
              IF( EmWaveMode ) THEN
                cx = CMPLX(x(j),x(j+1),KIND=dp)
                cflux = cmult * cx * CONJG(CMPLX(Fluxes(j),Fluxes(j+1),KIND=dp))
                crhs = cmult * cx * CONJG(CMPLX(b(j),b(j+1), KIND=dp))
              ELSE
                cflux = CMPLX(Fluxes(j),Fluxes(j+1),KIND=dp)
                crhs = CMPLX(b(j),b(j+1),KIND=dp)
              END IF
              IF( Nmode /= Mmode ) THEN
                FluxesRow(Mmode) = FluxesRow(Mmode) - REAL(cflux)
                FluxesRowIm(Mmode) = FluxesRowIm(Mmode) - AIMAG(cflux)
              ELSE
                FluxesRow(Nmode) = FluxesRow(Nmode) + REAL(cflux)
                FluxesRowIm(Nmode) = FluxesRowIm(Nmode) + AIMAG(cflux)
                FluxesRhs = FluxesRhs + REAL(crhs)
                FluxesRhsIm = FluxesRhsIm + AIMAG(crhs)
              END IF
            END IF
          ELSE
            flux = Fluxes(j)
            Mmode = k 
            IF( Nmode /= Mmode ) THEN                
              FluxesRow(Mmode) = FluxesRow(Mmode) - flux
            END IF
            FluxesRow(Nmode) = FluxesRow(Nmode) + flux
          END IF
        END IF
      END DO

      IF( EmWaveMode ) THEN
        IF(InfoActive(20)) THEN
          PRINT *,'ConstrainModesFluxes:',FluxesRow
          PRINT *,'ConstrainModesFluxesIm:',FluxesRowIm
          PRINT *,'ConstrainModesFluxes:',FluxesRhs,FluxesRhsIm        
        END IF
      END IF
        
    END SUBROUTINE ConstraintModesFluxes

   
    SUBROUTINE ConstraintModesCirculation()
      TYPE(Element_t), POINTER :: Element
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Solver_t), POINTER :: pSolver
      INTEGER, POINTER :: dofIndexes(:), nodeIndexes(:), Perm(:)
      REAL(KIND=dp), POINTER :: Values(:)
      REAL(KIND=dp) :: Normal(3), Normal0(3), Circ(3)
      TYPE(Nodes_t), SAVE :: Nodes
      INTEGER :: i,t,k,n,nb,elem
      LOGICAL :: NormalSet
            
      Mesh => Solver % Mesh
      IF( Mesh % MeshDim < 3 ) THEN
        CALL Fatal('ConstraintModesCirculation','Currently only available in 2D!')
      END IF
      pSolver => Solver

      FluxesRow = 0.0_dp
      IF( IsComplex ) FluxesRowIm = 0.0_dp
      
      DO elem=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        

        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes
        nodeIndexes => Element % NodeIndexes

        CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Element % NodeIndexes)
                      
        nb = mGetElementDOFs( dofIndexes, Element, Solver )

        k = Solver % Variable % ConstraintModesIndeces(Perm(dofIndexes(1)))
        IF(k==0) CYCLE
              
        IF(ANY(Var % ConstraintModesIndeces(Perm(dofIndexes(1:nb))) /= k)) CYCLE        
        
        ! Compute a reference normal to which all other face normals are compared against.
        Normal0 = NormalVector(Element,Nodes)

        CALL EdgeElementCirculation( Circ, Mesh, Element, Nodes, pSolver, Normal0, Var ) 
         
        ! Store elemental fluxes to the total sums. 
        IF( IsComplex ) THEN            
          Mmode = (k+1)/2
          IF( MOD(k,2) == 1 ) THEN                
            FluxesRow(Mmode) = FluxesRow(Mmode) + Circ(1)
          ELSE
            FluxesRowIm(Mmode) = FluxesRowIm(Mmode) + Circ(2)
          END IF
        ELSE
          Mmode = k 
          FluxesRow(Mmode) = FluxesRow(Mmode) + Circ(1)
        END IF
      END DO

    END SUBROUTINE ConstraintModesCirculation

    
    SUBROUTINE ConstraintModesLinkage()
      TYPE(Element_t), POINTER :: Element
      TYPE(Mesh_t), POINTER :: Mesh
      INTEGER, POINTER :: Indexes(:), Perm(:)
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t), SAVE :: Nodes
      REAL(KIND=dp) :: detJ, flux, s, PotAtIp, Basis(12), pot(12)
      INTEGER :: t,k,n,elem

      
      Mesh => Solver % Mesh
      IF( Mesh % MeshDim == 3 ) THEN
        CALL Fatal('ConstraintModesLinkage','Currently only available in 2D!')
      END IF

      FluxesRow = 0.0_dp
      IF( IsComplex ) FluxesRowIm = 0.0_dp
      Perm => Solver % Variable % Perm
      
      DO elem=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements        

        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes
        Indexes => Element % NodeIndexes        
        k = Var % ConstraintModesIndeces(Perm(Indexes(1)))

        IF(k==0) CYCLE
        
        IF(ANY(Var % ConstraintModesIndeces(Perm(Indexes(1:n))) /= k)) CYCLE        
        pot(1:n) = Solver % Variable % Values(Solver % Variable % Perm(Indexes(1:n)))

        CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, n, Indexes)
        IP = GaussPoints( Element )
        
        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )          
          s = IP % s(t) * DetJ
          PotAtIp = SUM(Basis(1:n)*pot(1:n))

          flux = s * PotAtIp 
          
          IF( IsComplex ) THEN
            Mmode = (k+1)/2
            IF( MOD(k,2) == 1 ) THEN                
              FluxesRow(Mmode) = FluxesRow(Mmode) + flux
            ELSE
              FluxesRowIm(Mmode) = FluxesRowIm(Mmode) + flux
            END IF
          ELSE
            Mmode = k 
            FluxesRow(Mmode) = FluxesRow(Mmode) + flux
          END IF
        END DO        
      END DO

    END SUBROUTINE ConstraintModesLinkage
      

    SUBROUTINE CommunicateConstraintModesFluxes()

      REAL(KIND=dp), ALLOCATABLE :: tmpFluxesRow(:)

      IF( Parallel ) THEN
        ALLOCATE(tmpFluxesRow(NoModes))
        tmpFluxesRow = FluxesRow
        CALL MPI_ALLREDUCE(tmpFluxesRow, FluxesRow, NoModes, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ELMER_COMM_WORLD, ierr)
        IF( IsComplex ) THEN
          tmpFluxesRow = FluxesRowIm
          CALL MPI_ALLREDUCE(tmpFluxesRow, FluxesRowIm, NoModes, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
          CALL MPI_ALLREDUCE(FluxesRhs, FluxesRhs, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
          CALL MPI_ALLREDUCE(FluxesRhsIm, FluxesRhsIm, 1, MPI_DOUBLE_PRECISION, &
              MPI_SUM, ELMER_COMM_WORLD, ierr)        
        END IF
        DEALLOCATE(tmpFluxesRow)
      END IF

    END SUBROUTINE CommunicateConstraintModesFluxes
    
    
!------------------------------------------------------------------------------
  END SUBROUTINE ConstraintModesDriver

SUBROUTINE SolveHarmonicSystem( G, Solver )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: G
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: BMatrix, A => NULL()
    INTEGER :: i,j,k,n, kr, ki, DOFs, ne, niter
    LOGICAL :: stat, Found, OptimizeBW, Real_given,Imag_given
    REAL(KIND=dp) :: Omega, norm, s
    REAL(KIND=dp), POINTER :: Freqv(:,:)
    REAL(KIND=dp), ALLOCATABLE :: x(:)
    REAL(KIND=dp), POINTER :: b(:)
    REAL(KIND=dp) :: frequency
    INTEGER :: Nfrequency
    TYPE(ValueList_t), POINTER :: BC
    CHARACTER(:), ALLOCATABLE :: Name

    CALL Info( 'SolveHarmonicSystem', 'Solving initially transient style system as harmonic one', Level=5)
    
    n = Solver % Matrix % NumberofRows
    DOFs = Solver % Variable % DOFs * 2

    A => G
    DO WHILE( ASSOCIATED(A) )
      BMatrix => A
      A => A % EMatrix
      IF ( ASSOCIATED(A) ) THEN
        IF ( A % COMPLEX ) THEN
          CALL Info('SolveHarmonicSystem','Reusing existing harmonic system',Level=10)
          EXIT
        END IF
      END IF
    END DO

    IF ( .NOT. ASSOCIATED(A) ) THEN      
      CALL Info('SolveHarmonicSystem','Creating new matrix for harmonic system',Level=10)      

      OptimizeBW = ListGetLogical(Solver % Values, 'Optimize Bandwidth', Found)
      IF ( .NOT. Found ) OptimizeBW = .TRUE.
      
      A => CreateMatrix( CurrentModel, Solver, Solver % Mesh,   &
              Solver % Variable % Perm, DOFs, MATRIX_CRS, OptimizeBW, &
              ListGetString( Solver % Values, 'Equation') )
      A % COMPLEX = .TRUE.
      BMatrix % EMatrix => A
      ALLOCATE( A % rhs(2*n) )
      
      DO j=1,Solver % Variable % DOFs
        Name = ComponentName( Solver % Variable % Name, j ) 
        DO i=1,CurrentModel % NumberOFBCs
          BC => CurrentModel % BCs(i) % Values
          real_given = ListCheckPresent( BC, Name )
          imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )
          
          IF ( real_given .AND. .NOT. imag_given ) THEN
            CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
          ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
            CALL ListAddConstReal( BC, Name, 0._dp )
          END IF
        END DO
      END DO
    END IF

    b => A % rhs
    ALLOCATE( x(2*n) )
    x = 0
    
    b(1:2*n:2) = G % RHS(1:n)
    b(2:2*n:2) = G % RHS_im(1:n)

    
    Nfrequency = ListGetInteger( Solver % Values,'Harmonic System Values',Found )
    IF( Nfrequency > 1 ) THEN
      freqv => ListGetConstRealArray( Solver % Values, 'Frequency' )
    ELSE
      Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
      IF( .NOT. Found ) THEN
        CALL Fatal( 'SolveHarmonicSystem', '> Frequency < must be given for harmonic analysis.' )
      END IF
      
      Nfrequency = 1
      ! Add the number of frequencies even for case of one for some postprocessing stuff to work 
      CALL ListAddInteger( Solver % Values,'Harmonic System Values',Nfrequency )
    END IF
    
    niter = MIN(Nfrequency,Solver % NOFEigenValues)
    ne=Solver % NofEigenValues
    Solver % NofEigenValues=0

    DO i=1,niter
      IF( Nfrequency > 1 ) THEN
        Frequency = freqv(i,1)
        WRITE( Message, '(a,i5,e12.3)' ) 'Frequency sweep: ', i, frequency
      ELSE
        WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
      END IF
      CALL Info( 'SolveHarmonicSystem', Message, Level=4 )

      omega = 2 * PI * Frequency
      DO k=1,n
        kr = A % Rows(2*(k-1)+1)
        ki = A % Rows(2*(k-1)+2)
        DO j=G % Rows(k),G % Rows(k+1)-1
          A % Values(kr)   =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(kr) = &
              A % Values(kr) - omega**2*G % MassValues(j)
          IF (ASSOCIATED(G % DampValues)) THEN
            A % Values(kr+1) = -G % Dampvalues(j) * omega
            A % Values(ki)   =  G % Dampvalues(j) * omega
          ELSE
            ! The damping is what couples the real and imaginary blocks. Without it
            ! the coupling is zero -- but it must still be written, as A is freshly
            ! allocated here and these are the only entries the loop would leave
            ! untouched.
            A % Values(kr+1) = 0.0_dp
            A % Values(ki)   = 0.0_dp
          END IF
          A % Values(ki+1) =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(ki+1) = &
            A % Values(ki+1) - omega**2*G % MassValues(j)
          kr = kr + 2
          ki = ki + 2
        END DO
      END DO

      
      DO j=1,Solver % Variable % DOFs
        Name = ComponentName( Solver % Variable % Name, j ) 

        CALL SetDirichletBoundaries( CurrentModel, A, b, Name, &
                2*j-1, DOFs, Solver % Variable % Perm )

        CALL SetDirichletBoundaries( CurrentModel, A, b, TRIM(Name) // ' im', &
                2*j, DOFs, Solver % Variable % Perm )
      END DO

      CALL EnforceDirichletConditions( Solver, A, b )
 
      
      CALL SolveLinearSystem( A, b, x, Norm, DOFs, Solver )
      
      DO j=1,n
        Solver % Variable % EigenVectors(i,j) = &
                 CMPLX( x(2*(j-1)+1),x(2*(j-1)+2),KIND=dp )
      END DO
    END DO

    Solver % NOFEigenValues = ne

    DEALLOCATE( x )
!------------------------------------------------------------------------------
 END SUBROUTINE SolveHarmonicSystem
!------------------------------------------------------------------------------



 

!------------------------------------------------------------------------------
!> Just toggles the initial system to harmonic one and back
!------------------------------------------------------------------------------
SUBROUTINE ChangeToHarmonicSystem( Solver, BackToReal )
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  LOGICAL, OPTIONAL :: BackToReal
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: Are => NULL(), Aharm => NULL(), SaveMatrix 
  INTEGER :: i,j,k,n, kr, ki, DOFs, TimeOrder
  LOGICAL :: stat, Found, OptimizeBW, Real_given, Imag_given
  CHARACTER(:), ALLOCATABLE :: Name
  REAL(KIND=dp) :: Omega, s, val
  REAL(KIND=dp), POINTER :: b(:), TmpVals(:)
  REAL(KIND=dp) :: frequency
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: TmpVar, ReVar, HarmVar, SaveVar
  LOGICAL :: ToReal, ParseName, AnyDirichlet, Diagonal, HarmonicReal, EigenMode
  CHARACTER(*), PARAMETER :: Caller = 'ChangeToHarmonicSystem'

  
  IF( .NOT. ASSOCIATED( Solver % Variable ) ) THEN
    CALL Warn(Caller,'Not applicable without a variable')
    RETURN    
  END IF

  IF( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Warn(Caller,'Not applicable without a matrix')
    RETURN    
  END IF

  EigenMode = ListgetLogical( Solver % Values, 'Eigen Analysis', Found )

  ToReal = .FALSE.
  IF( PRESENT( BackToReal ) ) ToReal = BackToReal

  IF( ToReal ) THEN
    IF( ASSOCIATED( Solver % Variable % Evar ) ) THEN
      IF( Solver % Variable % Evar % Dofs < Solver % Variable % Dofs ) THEN
        CALL Info(Caller,'Changing the harmonic results back to real system!',Level=6)

        SaveVar => Solver % Variable
        SaveMatrix => Solver % Matrix 

        Solver % Variable => Solver % Variable % Evar
        Solver % Variable % Evar => SaveVar

        Solver % Matrix => Solver % Matrix % EMatrix
        Solver % Matrix % Ematrix => SaveMatrix

        ! Eliminate cyclic dependence that is a bummer when deallocating stuff
        NULLIFY( Solver % Matrix % EMatrix % Ematrix )
      END IF
    END IF
    RETURN
  END IF


  CALL Info(Caller,'Changing the real transient system to harmonic one!',Level=6)

  SaveMatrix => Solver % Matrix
  SaveVar => Solver % Variable     

  n = Solver % Matrix % NumberofRows
  DOFs = SaveVar % Dofs
  Are => Solver % Matrix

  CALL Info(Caller,'Number of real system rows: '//I2S(n),Level=16)

  ! Obtain the frequency, it may depend on iteration step etc. 
  Omega = 0._dp
  IF (.NOT. EigenMode) THEN
    Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
    IF( .NOT. Found ) THEN
      CALL Fatal( Caller, '> Frequency < must be given for harmonic analysis.' )
    END IF
    WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
    CALL Info( Caller, Message, Level=5 )

     omega = 2 * PI * Frequency
     CALL ListAddConstReal( CurrentModel % Simulation, 'res: frequency', Frequency )
  END IF

  
  HarmonicReal = ListGetLogical( Solver % Values,'Harmonic Mode Real',Found ) 
  IF( HarmonicReal ) THEN
    CALL Info(Caller,'Enforcing harmonic system to be real valued',Level=8)
    IF (ASSOCIATED(Are % MassValues)) THEN
      Are % Values = Are % Values - omega**2* Are % MassValues
    ELSE
      CALL Fatal(Caller,'Harmonic system requires mass!')
    END IF
    ! This is set outside so that it can be called more flexibilly
    CALL EnforceDirichletConditions( Solver, Are, Are % rhs  )
    RETURN
  END IF

  TimeOrder = MAX( Solver % TimeOrder, &
      ListGetInteger(Solver % Values,'Time Derivative Order', Found ) ) 

  IF( TimeOrder == 2 ) THEN
    Diagonal = ListGetLogical( Solver % Values,'Harmonic Mode Block Diagonal',Found )  
    IF(.NOT. Found ) Diagonal = .NOT. ASSOCIATED(Are % DampValues)
    IF( Diagonal ) THEN
      CALL Info(Caller,'2nd order undamped system is assumed to be block diagonal',Level=8)
    END IF
  ELSE
    Diagonal = .FALSE.
    CALL Info(Caller,'1st order system is always assumed to be truly complex',Level=8)
  END IF
    
  
  ! Find whether the matrix already exists
  Aharm => Are % EMatrix
  IF( ASSOCIATED( Aharm ) ) THEN
    CALL Info(Caller,'Found existing harmonic system',Level=10)
    IF( ALLOCATED( Aharm % ConstrainedDOF ) ) Aharm % ConstrainedDOF = .FALSE.
  ELSE    
    ! Create the matrix if it does not
    
    Aharm => CreateChildMatrix( Are, Dofs, 2*Dofs, CreateRhs = .TRUE., Diagonal = Diagonal )

    Aharm % COMPLEX = ListGetLogical( Solver % Values,'Linear System Complex', Found ) 
    IF( .NOT. Found ) Aharm % COMPLEX = .NOT. Diagonal !TRUE. 
  END IF


  ! Set the harmonic system r.h.s
  b => Aharm % rhs
  
  IF( ASSOCIATED( Are % Rhs ) ) THEN
    b(1:2*n:2) = Are % RHS(1:n)
  ELSE
    b(1:2*n:2) = 0.0_dp
  END IF
  
  IF( ASSOCIATED( Are % Rhs_im ) ) THEN
    b(2:2*n:2) = Are % RHS_im(1:n)            
  ELSE
    b(2:2*n:2) = 0.0_dp
  END IF

  ! Mass matrix is always needed, both for 1st and 2nd order systems!
  ! It is always the leading time derivative. 
  IF( .NOT. ASSOCIATED(Are % MassValues) ) THEN
    CALL Fatal(Caller,'We do not have mass matrix values!')
  END IF

  IF( TimeOrder == 2 ) THEN
    IF( ASSOCIATED(Are % DampValues) ) THEN
      CALL Info(Caller,'We have damp matrix values',Level=12)
      IF( Diagonal ) THEN
        CALL Fatal(Caller,'Damping matrix cannot be block diagonal!')
      END IF
    ELSE
      CALL Info(Caller,'We do not have damp matrix values',Level=12)
    END IF
  END IF
    
  ! Set the harmonic system matrix
  IF( EigenMode ) THEN
    ALLOCATE(Aharm % MassValues(SIZE(Aharm % Values)))
    Aharm % MassValues = 0.0_dp

    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        Aharm % Values(kr) = Are % Values(j)
        Aharm % Values(ki+1) = Are % Values(j)
        
        IF( TimeOrder == 2 ) THEN
          IF (ASSOCIATED(Are % DampValues)) THEN
            Aharm % Values(kr+1) = -Are % Dampvalues(j)
            Aharm % Values(ki)   =  Are % Dampvalues(j)
          END IF
          
          Aharm % MassValues(kr) = Are % MassValues(j)
          Aharm % MassValues(ki+1) = Are % MassValues(j)
        ELSE
          Aharm % Values(kr+1) = -Are % Massvalues(j)
          Aharm % Values(ki)   =  Are % Massvalues(j)
        END IF
                    
        kr = kr + 2
        ki = ki + 2
      END DO
    END DO
  ELSE IF( Diagonal ) THEN
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)
        val = val - omega**2* Are % MassValues(j)
        
        Aharm % Values(kr) = val 
        Aharm % Values(ki) = val 
        kr = kr + 1
        ki = ki + 1
      END DO
    END DO
  ELSE
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)

        IF(TimeOrder == 2 ) THEN
          val = val - omega**2* Are % MassValues(j)        
          IF (ASSOCIATED(Are % DampValues)) THEN
            Aharm % Values(kr+1) = -Are % Dampvalues(j) * omega
            Aharm % Values(ki)   =  Are % Dampvalues(j) * omega
          END IF
        ELSE
          Aharm % Values(kr+1) = -Are % Massvalues(j) * omega
          Aharm % Values(ki)   =  Are % Massvalues(j) * omega
        END IF
                    
        Aharm % Values(kr) = val
        Aharm % Values(ki+1) = val     
        
        kr = kr + 2
        ki = ki + 2
      END DO
    END DO
  END IF
    
  AnyDirichlet = .FALSE.
  
  ! Finally set the Dirichlet conditions for the solver    
  DO j=1,DOFs
    Name = ComponentName( Solver % Variable % Name, j ) 
    DO i=1,CurrentModel % NumberOFBCs
      BC => CurrentModel % BCs(i) % Values
      real_given = ListCheckPresent( BC, Name )
      imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )

      IF( real_given .OR. imag_given ) AnyDirichlet = .TRUE.

      IF ( real_given .AND. .NOT. imag_given ) THEN
        CALL Info(Caller,'Setting zero >'//TRIM(Name)//' im< on BC '//I2S(i),Level=12)
        CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
      ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
        CALL Info(Caller,'Setting zero >'//TRIM(Name)//'< on BC '//I2S(i),Level=12)
        CALL ListAddConstReal( BC, Name, 0._dp )
      END IF
    END DO
  END DO


  IF( AnyDirichlet ) THEN
    DO j=1,DOFs
      Name = ComponentName( SaveVar % Name, j ) 
      
      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, Name, &
          2*j-1, 2*DOFs, SaveVar % Perm )

      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, TRIM(Name) // ' im', &
          2*j, 2*DOFs, SaveVar % Perm )
    END DO
  END IF

  
  ! Create the new fields, the total one and the imaginary one
  !-------------------------------------------------------------
  k = INDEX( SaveVar % name, '[' )
  ParseName = ( k > 0 ) 

  ! Name of the full complex variable not used for postprocessing
  IF( ParseName ) THEN
    Name = TRIM(SaveVar % Name(1:k-1))//' complex'
  ELSE
    Name = TRIM( SaveVar % Name )//' complex'
  END IF

  CALL Info(Caller,'Harmonic system full name: '//TRIM(Name),Level=12)


  HarmVar => VariableGet( Solver % Mesh % Variables, Name )
  IF( ASSOCIATED( HarmVar ) ) THEN
    CALL Info(Caller,'Reusing full system harmonic dofs',Level=12)
  ELSE
    CALL Info(Caller,'Creating full system harmonic dofs',Level=12)
    CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name,2*DOFs,Perm=SaveVar % Perm,Output=.FALSE.)
    HarmVar => VariableGet( Solver % Mesh % Variables, Name )
    IF(.NOT. ASSOCIATED( HarmVar ) ) CALL Fatal(Caller,'New created variable should exist!')

    ! Repoint the values of the original solution vector
    HarmVar % Values(1:2*n:2) = SaveVar % Values(1:n)

    ! It beats me why this cannot be deallocated without some NaNs later
    !DEALLOCATE( SaveVar % Values )
    SaveVar % Values => HarmVar % Values(1:2*n:2)
    SaveVar % Secondary = .TRUE.

    ! Repoint the components of the original solution
    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVar => VariableGet( Solver % Mesh % Variables, ComponentName( SaveVar % Name, i ) )
        IF( ASSOCIATED( TmpVar ) ) THEN
          TmpVar % Values => HarmVar % Values(2*i-1::HarmVar % Dofs)
        ELSE
          CALL Fatal(Caller,'Could not find re component '//I2S(i))
        END IF
      END DO
    END IF

    IF( ParseName ) THEN
      Name = ListGetString( Solver % Values,'Imaginary Variable',Found )
      IF(.NOT. Found ) THEN
        CALL Fatal(Caller,'We need > Imaginary Variable < to create harmonic system!')
      END IF
    ELSE
      Name = TRIM( SaveVar % Name )//' im'
      CALL Info(Caller,'Using derived name for imaginary component: '//TRIM(Name),Level=12)
    END IF

    TmpVals => HarmVar % Values(2:2*n:2)
    CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name, DOFs,TmpVals, Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        

    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVals => HarmVar % Values(2*i:2*n:2*Dofs)
        CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
            ComponentName(Name,i),1,TmpVals,Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        
      END DO
    END IF
    
  END IF

  IF ( EigenMode ) THEN 
     IF ( ASSOCIATED( Solver % Variable % EigenValues ) ) THEN
       HarmVar % EigenValues  => Solver % Variable % Eigenvalues
       HarmVar % EigenVectors => Solver % Variable % EigenVectors
     END IF
  END IF


  ! Now change the pointers such that when we visit the linear solver
  ! the system will automatically be solved as complex
  Solver % Variable => HarmVar
  Solver % Matrix => Aharm

  IF(AnyDirichlet) THEN
    IF(ParEnv % PEs>1) CALL ParallelInitMatrix(Solver, Aharm)
    CALL EnforceDirichletConditions( Solver, Aharm, b )
  END IF

  ! Save the original matrix and variable in Ematrix and Evar
  Solver % Matrix % Ematrix => SaveMatrix
  Solver % Variable % Evar => SaveVar    

  ! Eliminate cyclic dependence that is a bummer when deallocating stuff
  ! We are toggling {Are,Aharm} in {Solver % Matrix, Solver % Matrix % Ematrix}
  NULLIFY( Solver % Matrix % EMatrix % Ematrix )

!------------------------------------------------------------------------------
END SUBROUTINE ChangeToHarmonicSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Eliminate linear restriction only using the ListMatrix structure.
!> This is placed in a separate routine such that it can be called
!> when solving with and without block matrix being active.
!------------------------------------------------------------------------------
SUBROUTINE EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
    CollectionMatrix, Solver, CopyStiffMatrix, ExportUsePerm, ExportUseIPerm, ExportUseDiag )
  IMPLICIT NONE
  TYPE(Matrix_t) :: StiffMatrix
  REAL(KIND=dp) :: ForceVector(:) 
  TYPE(Matrix_t), POINTER :: RestMatrix
  TYPE(Matrix_t) :: CollectionMatrix  
  TYPE(Solver_t) :: Solver
  LOGICAL, OPTIONAL :: CopyStiffMatrix
  INTEGER, POINTER, OPTIONAL :: ExportUsePerm(:), ExportUseIPerm(:)
  REAL(KIND=dp), POINTER, OPTIONAL :: ExportUseDiag(:)
  
  INTEGER :: m,n,i,j,k,l,ix,p,q,Loop
  INTEGER, ALLOCATABLE, TARGET :: SlavePerm(:),MasterPerm(:),SlaveIPerm(:),MasterIPerm(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: SlaveDiag(:), MasterDiag(:), DiagDiag(:)
  INTEGER, POINTER :: UsePerm(:), UseIPerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:)
  REAL(KIND=dp) :: scl, val
  REAL(KIND=dp), POINTER :: TVals(:), Vals(:)
  REAL(KIND=dp), POINTER :: CollectionVector(:), RestVector(:)
  TYPE(ListMatrix_t), POINTER :: Lmat(:)
  TYPE(Matrix_t), POINTER :: Xmat, Tmat
  TYPE(ListMatrixEntry_t), POINTER :: cTmp
  LOGICAL :: Found, EliminateSlave, EliminateFromMaster, UseTranspose
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(*), PARAMETER :: Caller = 'EliminateLinearRestriction'
  

  CALL Info(Caller,'Eliminating Constraints from CollectionMatrix',Level=12)

  Params => Solver % Values


  EliminateSlave = ListGetLogical( Params, 'Eliminate Slave',Found )
  EliminateFromMaster = ListGetLogical( Params, 'Eliminate From Master',Found )

  UseTranspose = ListGetLogical(Params, 'Use Transpose values', Found)
  IF( UseTranspose ) THEN
    CALL Info(Caller,'Using transpose values in elimination',Level=15)            
  END IF

  
  n = StiffMatrix % NumberOfRows
  m = RestMatrix % NumberOfRows

  RestVector => NULL()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  IF(.NOT. ASSOCIATED(CollectionMatrix % Rhs) ) THEN
    ALLOCATE(CollectionMatrix % Rhs(n) )
    CollectionMatrix % Rhs = 0.0_dp
  END IF    
  CollectionVector => CollectionMatrix % RHS  

  ! We may optionally ask that the stiffness matrix is copied to the base.
  IF( PRESENT(CopyStiffMatrix)) THEN
    IF(CopyStiffMatrix) THEN
      DO i=StiffMatrix % NumberOfRows,1,-1
        DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
          CALL AddToMatrixElement( CollectionMatrix, &
              i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
        END DO
        CollectionVector(i) = CollectionVector(i) + ForceVector(i)
      END DO
    END IF
  END IF
    
  
  ALLOCATE(SlaveDiag(m),MasterDiag(m),SlavePerm(n),MasterPerm(n),&
      SlaveIPerm(m),MasterIPerm(m),DiagDiag(m))
  SlavePerm  = 0; SlaveIPerm  = 0; 
  MasterPerm = 0; MasterIPerm = 0
  SlaveDiag = 0.0_dp; MasterDiag = 0.0_dp
  DiagDiag = 0.0_dp
  
  Tvals => RestMatrix % TValues
  IF (.NOT.ASSOCIATED(Tvals)) Tvals => RestMatrix % Values 

  ! Extract diagonal entries for constraints:
  !------------------------------------------
  CALL Info(Caller,'Extracting diagonal entries for constraints',Level=15)


  DO i=1, RestMatrix % NumberOfRows
    m = RestMatrix % InvPerm(i)

    IF( m == 0 ) THEN
      CALL Warn(Caller,'InvPerm is zero for row: '//I2S(i))      
      CYCLE
    END IF

    m = MOD(m-1,n) + 1
    SlavePerm(m)  = i
    SlaveIperm(i) = m

    DO j=RestMatrix % Rows(i), RestMatrix % Rows(i+1)-1
      k = RestMatrix % Cols(j)
      val = Tvals(j)

      IF(k>n) THEN
        IF(ABS(val) > ABS(DiagDiag(i))) DiagDiag(i) = val
        CYCLE
      END IF

      ! Don't really really remember/understand the logic here but it seems better to
      ! choose the biggest value in case there are many of them. 
      IF(k == RestMatrix % InvPerm(i)) THEN
        IF(ABS(val) > ABS(SlaveDiag(i))) THEN
          SlaveDiag(i) = val
        END IF
      ELSE
        IF(ABS(val) > ABS(MasterDiag(i))) THEN
          MasterDiag(i) = val
          MasterPerm(k)  = i
          MasterIperm(i) = k
        END IF
      END IF
    END DO

    ! This is less conservative complaint than the original. 
    IF(ABS(SlaveDiag(i)) < TINY(val) .OR. ABS(MasterDiag(i)) < TINY(val)) THEN
      PRINT *,'Diagvals too small',ParEnv % MyPe,i,SlaveDiag(i),MasterDiag(i)
    END IF        
  END DO

  IF(InfoActive(25)) THEN
    PRINT *,'SlaveSum:',SUM(SlaveDiag)
    PRINT *,'MasterSum:',SUM(MasterDiag) 
    PRINT *,'SlaveSum abs:',SUM(ABS(SlaveDiag))
    PRINT *,'MasterSum abs:',SUM(ABS(MasterDiag))
  END IF

  IF(EliminateFromMaster) THEN
    CALL Info(Caller,'Eliminating from master',Level=15)      
    UsePerm  => MasterPerm 
    UseDiag  => MasterDiag
    UseIPerm => MasterIPerm 
  ELSE
    CALL Info(Caller,'Eliminating from slave',Level=15)            
    UsePerm  => SlavePerm
    UseDiag  => SlaveDiag
    UseIPerm => SlaveIPerm
  END IF
      
  IF(UseTranspose) THEN
    Vals => Tvals
  ELSE
    Vals => RestMatrix % Values
  END IF

  ! The rest is done in List Matrix format so move to that in case not yet!
  IF( CollectionMatrix % FORMAT /= MATRIX_LIST ) THEN
    CALL List_ToListMatrix(CollectionMatrix)
  END IF

  ! Replace elimination equations by the constraints (could done be as a postprocessing
  ! step, if eq's totally eliminated from linsys.)
  ! ----------------------------------------------------------------------------------
  CALL Info(Caller,'Deleting rows from equation to be eliminated',Level=15)

  Lmat => CollectionMatrix % ListMatrix
  DO m=1,RestMatrix % NumberOfRows
    i = UseIPerm(m)
    CALL List_DeleteRow(Lmat, i, Keep=.TRUE.)
  END DO

  CALL Info(Caller,'Copying rows from constraint matrix to eliminate dofs',Level=15)
  DO m=1,RestMatrix % NumberOfRows
    i = UseIPerm(m)
    DO l=RestMatrix % Rows(m+1)-1, RestMatrix % Rows(m), -1
      j = RestMatrix % Cols(l)

      ! skip l-coefficient entries, handled separately afterwards:
      ! --------------------------------------------------------
      IF(j > n) CYCLE
      CALL List_AddToMatrixElement( Lmat, i, j, Vals(l) )
    END DO
    CollectionVector(i) = RestVector(m)
  END DO

  ! Eliminate slave dof cycles:
  ! ---------------------------
  Xmat => RestMatrix
  Found = .TRUE.
  Loop = 0
  DO WHILE(Found)
    DO i=Xmat % NumberofRows,1,-1
      q = 0
      DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
        k = Xmat % Cols(j)
        IF(k>n) CYCLE
        IF(UsePerm(k)>0 .AND. ABS(TVals(j))>AEPS) q=q+1
      END DO
      IF(q>1) EXIT
    END DO
    Found = (q>1)

    Tmat => Xmat
    IF(Found) THEN
      Loop = Loop + 1
      CALL Info(Caller,'Recursive elimination round: '//I2S(Loop),Level=15)

      Tmat => AllocateMatrix()
      Tmat % Format = MATRIX_LIST

      DO i=Xmat % NumberofRows,1,-1
        DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
          k = Xmat % Cols(j)
          IF ( ABS(Tvals(j))>AEPS ) &
              CALL List_AddToMatrixElement(Tmat % ListMatrix, i, k, TVals(j))
        END DO
      END DO

      DO m=1,Xmat % NumberOfRows
        i = UseIPerm(m)
        DO j=Xmat % Rows(m), Xmat % Rows(m+1)-1
          k = Xmat % Cols(j)

          ! The size of SlavePerm is often exceeded but I don't really undersrtand the operation...
          ! so this is just a dirty fix.
          IF( k > SIZE( SlavePerm ) ) CYCLE

          l = SlavePerm(k)

          IF(l>0 .AND. k/=i) THEN
            IF(ABS(Tvals(j))<AEPS) CYCLE
            scl = -TVals(j) / SlaveDiag(l)

            CALL List_DeleteMatrixElement( Tmat % ListMatrix, m, k )

            DO q=Xmat % Rows(l+1)-1, Xmat % Rows(l),-1
              IF(ABS(Tvals(q))<AEPS) CYCLE
              ix = Xmat % Cols(q)
              IF ( ix/=k ) &
                  CALL List_AddToMatrixElement( Tmat % ListMatrix, m, ix, scl * TVals(q) )
            END DO
          END IF
        END DO
      END DO

      CALL List_ToCRSMatrix(Tmat)
      Tvals => Tmat % Values
      IF(.NOT.ASSOCIATED(Xmat,RestMatrix)) CALL FreeMatrix(Xmat)
    END IF
    Xmat => TMat
  END DO

  ! Eliminate Lagrange Coefficients:
  ! --------------------------------
  CALL Info(Caller,'Eliminating Lagrange Coefficients',Level=15)

  DO m=1,Tmat % NumberOfRows
    i = UseIPerm(m)
    IF( ABS( UseDiag(m) ) < TINY( 1.0_dp ) ) THEN
      PRINT *,'UseDiag too small:',m,ParEnv % MyPe,UseDiag(m)
      CYCLE
    END IF

    DO j=TMat % Rows(m), TMat % Rows(m+1)-1
      k = TMat % Cols(j)
      IF(k<=n) THEN
        IF(UsePerm(k)/=0) CYCLE
        scl = -Tvals(j) / UseDiag(m)
      ELSE
        k = UseIPerm(k-n)
        scl = -Tvals(j) / UseDiag(m)
      END IF

      DO l=StiffMatrix % Rows(i+1)-1, StiffMatrix % Rows(i),-1
        CALL List_AddToMatrixElement( Lmat, k, &
            StiffMatrix % Cols(l), scl * StiffMatrix % Values(l) )
      END DO
      CollectionVector(k) = CollectionVector(k) + scl * ForceVector(i)
    END DO
  END DO

  IF ( .NOT.ASSOCIATED(Tmat, RestMatrix ) ) CALL FreeMatrix(Tmat)

  ! Eliminate slave dofs, using the constraint equations:
  ! -----------------------------------------------------
  IF ( EliminateSlave ) THEN
    CALL Info(Caller,'Eliminate slave dofs using constraint equations',Level=15)

    CALL List_ToCRSMatrix(CollectionMatrix)
    Tmat => AllocateMatrix()
    Tmat % Format = MATRIX_LIST

    DO i=1,StiffMatrix % NumberOfRows
      IF(UsePerm(i)/=0) CYCLE

      DO m = CollectionMatrix % Rows(i), CollectionMatrix % Rows(i+1)-1
        j = SlavePerm(CollectionMatrix % Cols(m))

        IF(j==0) THEN
          CYCLE
        END IF
        IF( ABS( SlaveDiag(j) ) < TINY( 1.0_dp ) ) THEN
          PRINT *,'SlaveDiag too small:',j,ParEnv % MyPe,SlaveDiag(j)
          CYCLE
        END IF

        scl = -CollectionMatrix % Values(m) / SlaveDiag(j)
        CollectionMatrix % Values(m) = 0._dp

        ! ... and add replacement values:
        ! -------------------------------
        k = UseIPerm(j)
        DO p=CollectionMatrix % Rows(k+1)-1, CollectionMatrix % Rows(k), -1
          l = CollectionMatrix % Cols(p)
          IF ( l /= SlaveIPerm(j) ) &
              CALL List_AddToMatrixElement( Tmat % listmatrix, i, l, scl*CollectionMatrix % Values(p) )
        END DO
        CollectionVector(i) = CollectionVector(i) + scl * CollectionVector(k)
      END DO
    END DO

    CALL List_ToListMatrix(CollectionMatrix)
    Lmat => CollectionMatrix % ListMatrix

    CALL List_ToCRSMatrix(Tmat)
    DO i=TMat % NumberOfRows,1,-1
      DO j=TMat % Rows(i+1)-1,TMat % Rows(i),-1
        CALL List_AddToMatrixElement( Lmat, i, TMat % cols(j), TMat % Values(j) )
      END DO
    END DO
    CALL FreeMatrix(Tmat)
  END IF

  ! Optimize bandwidth, if needed:
  ! ------------------------------
  IF(EliminateFromMaster) THEN
    CALL Info(Caller,'Optimizing bandwidth after elimination',Level=15)
    DO i=1,RestMatrix % NumberOfRows
      j = SlaveIPerm(i)
      k = MasterIPerm(i)

      Ctmp => Lmat(j) % Head
      Lmat(j) % Head => Lmat(k) % Head
      Lmat(k) % Head => Ctmp

      l = Lmat(j) % Degree
      Lmat(j) % Degree = Lmat(k) % Degree
      Lmat(k) % Degree = l

      scl = CollectionVector(j)
      CollectionVector(j) = CollectionVector(k)
      CollectionVector(k) = scl
    END DO
  END IF

  IF( PRESENT(ExportUsePerm) ) THEN
    CALL Info(Caller,'Export UsePerm outside elimination',Level=20) 
    ALLOCATE(ExportUsePerm(SIZE(UsePerm)))
    ExportUsePerm = UsePerm
  END IF
  IF( PRESENT(ExportUseIPerm) ) THEN
    CALL Info(Caller,'Export UseIPerm outside elimination',Level=20) 
    ALLOCATE(ExportUseIPerm(SIZE(UseIPerm)))
    ExportUseIPerm = UseIPerm
  END IF
  IF( PRESENT(ExportUseDiag) ) THEN
    CALL Info(Caller,'Export UseDiag outside elimination',Level=20) 
    ALLOCATE(ExportUseDiag(SIZE(UseDiag)))
    ExportUseDiag = UseDiag
  END IF

  IF(PRESENT(CopyStiffMatrix)) THEN
    IF(CopyStiffMatrix) CALL List_ToCRSMatrix(CollectionMatrix)
  END IF
      
  CALL Info(Caller,'Finished Eliminating Restrictions',Level=12)

END SUBROUTINE EliminateLinearRestriction

  

!------------------------------------------------------------------------------
!>  This subroutine will solve the system with some linear restriction.
!>  The restriction matrix is assumed to be in the ConstraintMatrix-field of 
!>  the StiffMatrix. The restriction vector is the RHS-field of the
!>  ConstraintMatrix.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE SolveWithLinearRestriction( StiffMatrix, ForceVector, &
    Solution, Norm, DOFs, Solver )
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
  REAL(KIND=dp),TARGET :: ForceVector(:) !< The right hand side of the linear equation
  REAL(KIND=dp),TARGET :: Solution(:)    !< Previous solution as input, new solution as output.
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedom of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPointer
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, RestMatrixTranspose 
  REAL(KIND=dp), POINTER CONTIG :: CollectionVector(:), RestVector(:), AddVector(:) 
  REAL(KIND=dp), POINTER  :: MultiplierValues(:), pSol(:),DiagScaling(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat, NoEmptyRows 
  INTEGER :: i, j, k, l, m, n, p,q, ix, Loop, colj, nIter
  TYPE(Variable_t), POINTER :: MultVar, iterV
  REAL(KIND=dp) :: scl, rowsum, Relax, val
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet, &
      NonEmptyRow, ComplexSystem, ConstraintScaling, UseTranspose, EliminateConstraints, &
      SkipConstraints, ResidualMode
  INTEGER, POINTER :: UseIPerm(:), UsePerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:) 
  LOGICAL  :: Parallel, UseTreeGauge, NeedMassDampValues, DoOwnScaling
  LOGICAL, ALLOCATABLE :: TrueDof(:)
  INTEGER, ALLOCATABLE :: Iperm(:)
  REAL(KIND=dp) :: t0,rt0,st,rst
  CHARACTER(:), ALLOCATABLE :: str,MultiplierName
  TYPE(ValueList_t), POINTER :: Params
  CHARACTER(*), PARAMETER :: Caller = 'SolveWithLinearRestriction'

  TYPE(ParEnv_t), POINTER :: ParEnvSave

  SAVE MultiplierValues, SolverPointer
  
!------------------------------------------------------------------------------
  CALL Info( Caller, ' ', Level=12 )
  ParEnvSave => ParEnv

  SolverPointer => Solver  
  Params => Solver % Values

  t0 = CPUTime()
  rt0 = RealTime()
    
  Parallel = Solver % Parallel

  ResidualMode = ListGetLogical( Params,'Restriction System Residual Mode',Found )    
  iterV => VariableGet(Solver % Mesh % Variables,'nonlin iter',UnfoundFatal=.TRUE.)
  nIter = NINT(iterV % Values(1))
    
  NotExplicit = ListGetLogical(Params,'No Explicit Constrained Matrix',Found)
  IF(.NOT. Found) NotExplicit=.FALSE.

  NeedMassDampValues = ListGetLogical(Params, 'Eigen Analysis', Found )

  RestMatrix => NULL()
  IF(.NOT.NotExplicit) RestMatrix => StiffMatrix % ConstraintMatrix
  RestVector => Null()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  AddMatrix => StiffMatrix % AddMatrix
  AddVector => NULL()
  IF(ASSOCIATED(AddMatrix)) AddVector => AddMatrix % RHS

  EliminateConstraints = ListGetLogical( Params, 'Eliminate Linear Constraints', Found)
  
  NumberOfRows = StiffMatrix % NumberOfRows
  
  CollectionMatrix => StiffMatrix % CollectionMatrix
  Refactorize = ListGetLogical(Params,'Linear System Refactorize',Found)
  IF(.NOT.Found) THEN
    Refactorize = .NOT. ( ResidualMode .AND. nIter > 1) 
  END IF

  IF(ASSOCIATED(CollectionMatrix)) THEN
    IF(Refactorize .AND. .NOT.NotExplicit) THEN
      CALL Info( Caller,'Freeing previous collection matrix structures',Level=10)
      CALL FreeMatrix(CollectionMatrix)
      CollectionMatrix => NULL()
      CALL Info( Caller,'Refactoring requested, creating fully new matrix',Level=10)
    ELSE
      CALL Info( Caller,'Trying to keep previous collection matrix structures',Level=10)
    END IF
  END IF
  ! No need to release anything from > Solver % ParEnv < here: it is a mirror
  ! only, the collection matrix created below gets an environment of its own.

  IF(.NOT.ASSOCIATED(CollectionMatrix)) THEN
    CollectionMatrix => AllocateMatrix()
    CollectionMatrix % FORMAT = MATRIX_LIST
  ELSE
    DEALLOCATE(CollectionMatrix % RHS)
    CollectionMatrix % Values = 0.0_dp

    IF(NeedMassDampValues) THEN
      IF(ASSOCIATED(CollectionMatrix % MassValues)) CollectionMatrix % MassValues = 0.0_dp
      IF(ASSOCIATED(CollectionMatrix % DampValues)) CollectionMatrix % DampValues = 0.0_dp
    END IF
  END IF
  IF(NotExplicit) CollectionMatrix % ConstraintMatrix => StiffMatrix % ConstraintMatrix  
  
  NumberOfRows = StiffMatrix % NumberOfRows
  IF(ASSOCIATED(AddMatrix)) NumberOfRows = MAX(NumberOfRows,AddMatrix % NumberOfRows)
  IF(ASSOCIATED(RestMatrix)) THEN
    IF(.NOT.EliminateConstraints) NumberOfRows = NumberOFRows + RestMatrix % NumberOfRows
  END IF

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
       CollectionSolution( NumberOfRows ), STAT = istat )
  IF ( istat /= 0 ) CALL Fatal( Caller, 'Memory allocation error.' )

  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

  ComplexSystem = StiffMatrix % COMPLEX .OR. &
      ListGetLogical(Params,'Linear System Complex', Found )
  
!------------------------------------------------------------------------------
! If multiplier should be exported,  allocate memory and export the variable.
!------------------------------------------------------------------------------

  ExportMultiplier = ListGetLogical(Params, 'Export Lagrange Multiplier', Found )
  IF ( ExportMultiplier ) THEN
    MultiplierName = LagrangeMultiplierName( Solver )
    MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
    j = 0
    IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberofRows
    IF(ASSOCIATED(AddMatrix))  j = j+MAX(0,AddMatrix % NumberofRows-StiffMatrix % NumberOfRows)

    IF ( .NOT. ASSOCIATED(MultVar) ) THEN
      CALL Info(Caller,'Creating variable for Lagrange multiplier',Level=8)
      ALLOCATE( MultiplierValues(j), STAT=istat )
      IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error.')

      MultiplierValues = 0.0_dp
      IF( ComplexSystem ) THEN
        CALL VariableAddVector(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
            MultiplierName, 2, MultiplierValues)               
      ELSE
        CALL VariableAdd(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
            MultiplierName, 1, MultiplierValues)      
      END IF
      MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
    END IF

    IF( InfoActive( 20 ) ) THEN
      CALL VectorValuesRange(MultVar % Values,SIZE(MultVar % Values),TRIM(MultVar % Name))
    END IF

    MultiplierValues => MultVar % Values

    IF (j > SIZE(MultiplierValues)) THEN
      CALL Info(Caller,'Increasing Lagrange multiplier size to: '//I2S(j),Level=8)
      ALLOCATE(MultiplierValues(j)); MultiplierValues=0._dp       
      MultiplierValues(1:SIZE(MultVar % Values)) = MultVar % Values     

      ! If the Lagrange variable includes history also change its size.
      IF( ASSOCIATED( MultVar % PrevValues ) ) THEN
        MultVar % Values = MultVar % PrevValues(:,1)
        DEALLOCATE( MultVar % PrevValues )
        ALLOCATE( MultVar % PrevValues(j,1) )
        MultVar % PrevValues = 0.0_dp
        MultVar % PrevValues(:,1) = MultVar % Values
      END IF

      DEALLOCATE(MultVar % Values)
      MultVar % Values => MultiplierValues
    END IF

    IF( InfoActive(25) ) THEN
      CALL VectorValuesRange(MultVar % values,SIZE(MultVar % values),'MultVar')
    END IF
  ELSE
    MultiplierValues => NULL()
  END IF

  UseTreeGauge = ListGetlogical(Params, 'Use Tree Gauge', Found )

!------------------------------------------------------------------------------
! Put the RestMatrix to lower part of CollectionMatrix
!------------------------------------------------------------------------------

  EnforceDirichlet = ListGetLogical(Params, 'Enforce Exact Dirichlet BCs',Found)
  IF(.NOT.Found) EnforceDirichlet = .TRUE.
  EnforceDirichlet = EnforceDirichlet .AND. ALLOCATED(StiffMatrix % ConstrainedDOF)

  UseTranspose = ListGetLogical(Params, 'Use Transpose values', Found)
  IF( UseTranspose ) THEN
    CALL Info(Caller,'Using transpose values in elimination',Level=15)            
  END IF
    
  CALL Info(Caller,'Number of Rows / Nonzeros in original matrix: '&
      //I2S(StiffMatrix % NumberOfRows)//' / '//I2S(SIZE(StiffMatrix % Values)),Level=22) 
  
  IF(ASSOCIATED(RestMatrix) .AND. .NOT. EliminateConstraints) THEN

    CALL Info(Caller,'Adding ConstraintMatrix into CollectionMatrix',Level=10)

    CALL Info(Caller,'Number of Rows / Nonzeros in constraint matrix: '&
      //I2S(RestMatrix % NumberOfRows)//' / '//I2S(SIZE(RestMatrix % Values)),Level=12) 

    NoEmptyRows = 0
    ConstraintScaling = ListGetLogical(Params, 'Constraint Scaling',Found)
    IF(ConstraintScaling) THEN
      rowsum = ListGetConstReal(Params, 'Constraint Scale', Found)
      IF(Found) RestMatrix % Values = RestMatrix % Values * rowsum
    END IF

    ALLOCATE( iperm(SIZE(Solver % Variable % Perm)) )
    iperm = 0
    DO i=1,SIZE(Solver % Variable % Perm)
      IF ( Solver % Variable % Perm(i)>0) Iperm(Solver % Variable % Perm(i))=i
    END DO

    DO i=RestMatrix % NumberOfRows,1,-1

      k=StiffMatrix % NumberOfRows
      IF(ASSOCIATED(AddMatrix)) k=MAX(k,AddMatrix % NumberOfRows)
      k=k+i

      CALL AddToMatrixElement( CollectionMatrix,k,k,0._dp )
      IF(ComplexSystem) THEN
        IF(MOD(k,2)==0) THEN
          CALL AddToMatrixElement( CollectionMatrix,k,k-1,0._dp )
        ELSE
          CALL AddToMatrixElement( CollectionMatrix,k,k+1,0._dp )
        END IF
      END IF
      NonEmptyRow = .FALSE.

      rowsum = 0._dp
      l = -1
      DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
        IF(RestMatrix % Cols(j)==k) l=j
        rowsum = rowsum + ABS(RestMatrix % Values(j))
      END DO

      IF(rowsum>EPSILON(1._dp)) THEN
        IF(ConstraintScaling) THEN
          IF(l<=0.OR.l>0.AND.RestMatrix % Values(l)==0) THEN
            DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
              RestMatrix % Values(j) = RestMatrix % values(j)/rowsum
            END DO
            RestMatrix % RHS(i) = RestMatrix % RHS(i) / rowsum
          END IF
        END IF

        DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
          ! Skip non-positive column indexes, why should there be any? 
          colj = RestMatrix % Cols(j)
          IF( colj <= 0 ) CYCLE

          ! Complex system requires that we have exact 2x2 block structure. Don't spoil that.
          IF ( .NOT. ComplexSystem ) THEN
            IF( ABS(RestMatrix % Values(j)) < EPSILON(1._dp)*rowsum ) CYCLE
          END IF

          ! If we have Dirichlet condition set for the matrix use that directly and do not add
          ! stuff to the row that would spoil the condition. 
          Found = .TRUE.
          IF (EnforceDirichlet .AND. colj <= StiffMatrix % NumberOfRows) THEN
            Found = .NOT. StiffMatrix % ConstrainedDOF(colj)
          END IF
            
          IF(Found) THEN
            IF (ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE
              val = RestMatrix % Values(j)
            END IF              
            CALL AddToMatrixElement( CollectionMatrix, colj, k, val ) 

            ! Add the Transpose part
            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE
              val = RestMatrix % Values(j)
            END IF

            ! Only add the transpose when it is associated to the unknowns of the initial matrix.
            ! Otherwise the entries related to largrange multipliers would be multiplied by factor 2!
            IF( colj <= StiffMatrix % NumberOfRows ) THEN            
              CALL AddToMatrixElement( CollectionMatrix, k, colj, val ) 
              NonEmptyRow = NonEmptyRow .OR. val /= 0
            END IF
          ELSE
            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              val = RestMatrix % TValues(j)
            ELSE 
              val = RestMatrix % Values(j)
            END IF
            ! Use the value of the Dirichlet condition from "Dvalues"
            CollectionVector(k) = CollectionVector(k) - val * StiffMatrix % DValues(colj)  
          END IF
        END DO
      END IF
 
      Found = .TRUE.
      IF (EnforceDirichlet) THEN
        IF(ASSOCIATED(RestMatrix % InvPerm)) THEN
          l = RestMatrix % InvPerm(i)
          IF(l>0) THEN
            l = MOD(l-1,StiffMatrix % NumberOfRows)+1
            IF(StiffMatrix % ConstrainedDOF(l)) THEN
              l = iperm((l-1)/Solver % Variable % DOFs+1) 
              IF (l<=Solver % Mesh % NumberOfNodes) THEN
                Found = .FALSE.
                CALL ZeroRow(CollectionMatrix,k)
                CollectionVector(k) = 0
                CALL SetMatrixElement(CollectionMatrix,k,k,1._dp)
              END IF
            END IF
          END IF
        END IF
      END IF

      ! If there is no matrix entry, there can be no non-zero r.h.s.
      IF ( Found ) THEN
        IF( .NOT.NonEmptyRow ) THEN
          NoEmptyRows = NoEmptyRows + 1
          CollectionVector(k) = 0._dp
!          might not be the right thing to do in parallel!!
          IF(UseTreeGauge) THEN
            CALL SetMatrixElement( CollectionMatrix,k,k,1._dp )
          END IF 
        ELSE
          IF( ASSOCIATED( RestVector ) ) CollectionVector(k) = CollectionVector(k) + RestVector(i)
        END IF
      END IF
    END DO

    IF( NoEmptyRows > 0 ) THEN
      CALL Info(Caller,'Constraint Matrix has '//I2S(NoEmptyRows)// &
          ' empty rows out of '//I2S(RestMatrix % NumberOfRows),Level=6 )
    END IF

    CALL Info(Caller,'Finished Adding ConstraintMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the AddMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  IF(ASSOCIATED(AddMatrix)) THEN
    CALL Info(Caller,'Adding AddMatrix into CollectionMatrix',Level=12)

    CALL Info(Caller,'Number of Rows / Nonzeros in additional matrix: '&
      //I2S(AddMatrix % NumberOfRows)//' / '//I2S(SIZE(AddMatrix % Values)),Level=12) 
    
    DO i=AddMatrix % NumberOfRows,1,-1

      Found = .TRUE.
      IF (EnforceDirichlet .AND. i<=StiffMatrix % NumberOFRows) THEN
        Found = .NOT. StiffMatrix % ConstrainedDOF(i)
      END IF
        
      IF(Found) THEN
        Found = .FALSE.
        DO j=AddMatrix % Rows(i+1)-1,AddMatrix % Rows(i),-1
            CALL AddToMatrixElement( CollectionMatrix, &
               i, AddMatrix % Cols(j), AddMatrix % Values(j))
            IF (i == AddMatrix % Cols(j)) Found = .TRUE.
        END DO
        
        IF( ASSOCIATED(AddVector)) THEN
          CollectionVector(i) = CollectionVector(i) + AddVector(i)
        END IF
        IF (.NOT.Found) THEN
          CALL AddToMatrixElement( CollectionMatrix, i, i, 0._dp )
          IF(ComplexSystem) THEN
            IF(MOD(i,2)==0) THEN
              CALL AddToMatrixElement( CollectionMatrix,i,i-1,0._dp )
            ELSE
              CALL AddToMatrixElement( CollectionMatrix,i,i+1,0._dp )
            END IF
          END IF
        END IF
      END IF
    END DO
    CALL Info(Caller,'Finished Adding AddMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the StiffMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  CALL Info(Caller,'Adding Stiffness Matrix into CollectionMatrix',Level=12)

  DO i=StiffMatrix % NumberOfRows,1,-1
    DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
      CALL AddToMatrixElement( CollectionMatrix, &
        i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
    END DO
    CollectionVector(i) = CollectionVector(i) + ForceVector(i)
  END DO

!------------------------------------------------------------------------------
! Eliminate constraints instead of adding the Lagrange coefficient equations.
! Assumes biorthogonal basis for Lagrange coefficient interpolation, but not
! necessarily biorthogonal constraint equation test functions.
!------------------------------------------------------------------------------
  IF (ASSOCIATED(RestMatrix) .AND. EliminateConstraints) THEN
    IF ( ExportMultiplier ) THEN
      ! With the multiplier active we need to use it also for elimination in case the
      ! constraint shares some dofs with the multiplier. 
      CALL EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
          CollectionMatrix, Solver, ExportUseIPerm = UseIPerm, ExportUseDiag = UseDiag )
    ELSE
      CALL EliminateLinearRestriction( StiffMatrix, ForceVector, RestMatrix, &
        CollectionMatrix, Solver )
    END IF
  END IF
  
  IF(CollectionMatrix % FORMAT==MATRIX_LIST) THEN
    CALL Info(Caller,'Reverting CollectionMatrix back to CRS matrix',Level=10)
    CALL List_toCRSMatrix(CollectionMatrix)
  END IF
    
  ! CRS-format matrix needed here
  IF ( NeedMassDampValues ) THEN  ! Doesn't work with constraints, "AddMatrix" only !!
    CALL CopyMassDampValues(CollectionMatrix, StiffMatrix, AddMatrix)
  END IF
  
  CALL Info( Caller, 'CollectionMatrix done', Level=12 )

!------------------------------------------------------------------------------
! Assign values to CollectionVector
!------------------------------------------------------------------------------

  j = StiffMatrix % NumberOfRows  
  CollectionSolution(1:j) = Solution(1:j)
  
  i = StiffMatrix % NumberOfRows+1
  j = SIZE(CollectionSolution)
  IF( j >= i) THEN
    CollectionSolution(i:j) = 0._dp
    IF(ExportMultiplier) CollectionSolution(i:j) = MultiplierValues(1:j-i+1)
  END IF

  IF( InfoActive(30) ) THEN
    pSol => CollectionSolution
    CALL VectorValuesRange(pSol,j,'CollectionSolution')           
  END IF
  
  CollectionMatrix % ExtraDOFs = CollectionMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows

  CollectionMatrix % ParallelDOFs = 0
  IF(ASSOCIATED(AddMatrix)) &
    CollectionMatrix % ParallelDOFs = MAX(AddMatrix % NumberOfRows - &
            StiffMatrix % NumberOfRows,0)
    
  CALL Info( Caller, 'CollectionVector done', Level=12 )

!------------------------------------------------------------------------------
! Solve the Collection-system 
!------------------------------------------------------------------------------

! Collectionmatrix % Complex = StiffMatrix % Complex


  CollectionMatrix % Comm = StiffMatrix % Comm

  
  st  = CPUTime() - t0;
  rst = RealTime() - rt0  
  WRITE(Message,'(a,f8.2,f8.2,a)') 'Collection matrix creation time (CPU,REAL): ',st,rst,' (s)'
  CALL Info(Caller,Message,Level=6)    
  
  i = CollectionMatrix % NumberOfRows - StiffMatrix % NumberOfRows
  j = SIZE(CollectionMatrix % Values) - SIZE(StiffMatrix % Values )
  CALL Info(Caller,'Collection matrix increased with '//I2S(i)//&
      ' rows and '//I2S(j)//' non-zeros',Level=8)
    
  IF( InfoActive( 30 ) ) THEN
    CALL VectorValuesRange(CollectionMatrix % Values,SIZE(CollectionMatrix % Values),'A')       
    CALL VectorValuesRange(CollectionMatrix % rhs,SIZE(CollectionMatrix % rhs),'b')       
  END IF
      
  IF( ResidualMode ) THEN
    BLOCK
      REAL(KIND=dp), POINTER :: Res(:)      
      ! If residual mode is requested make change of variables:
      ! Ax=b -> Adx = b-Ax0 = r
      IF( niter > 1 ) THEN
        CALL Info(Caller,'Changing the equation to residual based mode',Level=10)
        ALLOCATE( Res(SIZE(CollectionSolution)) ) 
        CALL LinearSystemResidual( CollectionMatrix, CollectionVector, CollectionSolution, res )
        CollectionVector = Res
        CollectionSolution = 0.0_dp
        DEALLOCATE(Res)
      END IF
    END BLOCK
  END IF
  
  ! We may want to skip ComputeChange including the constraints if we use certain other options
  SkipConstraints = ResidualMode .OR. &
      ListGetLogical( Params, 'Nonlinear System Convergence Without Constraints',Found ) .OR. &
      ListGetLogical( Params, 'NonLinear System Consistent Norm',Found )   
  str = ListGetString( Params, 'NonLinear System Convergence Measure',Found )
  IF( str == 'solution' ) THEN
    SkipConstraints = .TRUE.
    CALL Info(Caller,&
        'Nonlinear system convergence measure == "solution" must skip constraints',Level=10)
  END IF
  IF( SkipConstraints ) THEN
    CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.TRUE.)
  END IF

  DoOwnScaling = ListGetLogical( Params,'Linear System Scaling',Found)
  IF(.NOT. Found) DoOwnScaling = .TRUE.
  IF(.NOT. ResidualMode) DoOwnScaling = .FALSE.
  IF(DoOwnScaling) THEN   
    CALL Info(Caller,'Performing special scaling with constraints',Level=10)
    DiagScaling => CollectionMatrix % DiagScaling
    IF(Niter == 1 ) THEN
      IF(.NOT. ASSOCIATED(DiagScaling) ) THEN
        ALLOCATE( DiagScaling(SIZE(CollectionVector)))
        CollectionMatrix % DiagScaling => DiagScaling
      END IF

      ! Should we scale only part or the full matrix? 
      IF(.FALSE.) THEN
        DiagScaling = 1.0_dp
        StiffMatrix % DiagScaling => DiagScaling
        ! Just build the scaling matrix using only the original stiffness matrix.
        CALL ScaleLinearSystem(Solver,StiffMatrix,ApplyScaling=.FALSE.)     
        CollectionMatrix % ScalingMethod = StiffMatrix % ScalingMethod    
        StiffMatrix % DiagScaling => NULL()      
      ELSE
        CALL ScaleLinearSystem(Solver,CollectionMatrix,ApplyScaling=.FALSE.)     
      END IF        
    END IF
    
    CALL ScaleLinearSystem(Solver,CollectionMatrix,CollectionVector,&
        CollectionSolution,DiagScaling=CollectionMatrix % DiagScaling)
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.TRUE. ) 
  END IF
  
  
  !IF( ListGetLogical( Params,'Linear System Save',Found ) ) THEN        
  !  CALL SaveLinearSystem( Solver, CollectionMatrix,'RestrictedMat')
  !END IF
  

  CALL Info(Caller,'Now solving the linear system with constraints!',Level=10)
  Collectionmatrix % DGMatrix = StiffMatrix %  DGMatrix

  CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
     CollectionSolution, Norm, DOFs, Solver, StiffMatrix )


  IF(DoOwnScaling) THEN
    CALL BackScaleLinearSystem( Solver,CollectionMatrix,CollectionVector,&
        CollectionSolution,CollectionMatrix % DiagScaling)
    CALL ListAddLogical( Params,'Linear System Skip Scaling',.FALSE. ) 
  END IF

  !-------------------------------------------------------------------------------
  ! For restricted systems study the norm without some block components.
  ! For example, excluding gauge constraints may give valuable information
  ! of the real accuracy of the unconstrained system. Currently just for info.
  !-------------------------------------------------------------------------------
  IF( ListGetLogical( Params,'Restricted System Norm',Found ) ) THEN
    ALLOCATE( TrueDof( CollectionMatrix % NumberOfRows ) )
    TrueDof = .TRUE.
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the original system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    IF( ListGetLogical( Params,'Restricted System Norm Skip Nodes',Found ) ) THEN
      i = 1
      j = MAXVAL( Solver % Variable % Perm(1:Solver % Mesh % NumberOfNodes) )
      CALL Info(Caller,'Skipping nodal dof range: '//I2S(i)//'-'//I2S(j),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF

    IF( ListGetLogical( Params,'Restricted System Norm Skip Constraints',Found ) ) THEN
      i = StiffMatrix % NumberOfRows + 1
      j = CollectionMatrix % NumberOfRows      
      CALL Info(Caller,'Skipping constraints dof range: '&
          //I2S(i)//'-'//I2S(j),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the masked system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    DEALLOCATE( TrueDof )
  END IF
    
  
!------------------------------------------------------------------------------
! Separate the solution from CollectionSolution
!------------------------------------------------------------------------------
    CALL Info(Caller,'Picking solution from collection solution',Level=10)

    j = StiffMatrix % NumberOfRows

    IF( ResidualMode .AND. nIter > 1) THEN
      Solution(1:j) = Solution(1:j) + CollectionSolution(1:j)
    ELSE
      Solution(1:j) = CollectionSolution(1:j)
    END IF
    
    IF ( ExportMultiplier ) THEN
      CALL Info(Caller,'Separating Lagrange multiplier from collection solution',Level=10)
      
      IF(ASSOCIATED(RestMatrix) .AND. EliminateConstraints) THEN        
        ! Compute eliminated l-coefficient values:
        ! ---------------------------------------
        MultiplierValues = 0.0_dp
        DO i=1,RestMatrix % NumberOfRows
          scl = 1._dp / UseDiag(i)
          m = UseIPerm(i)
          MultiplierValues(i) = scl * ForceVector(m)
          DO j=StiffMatrix % Rows(m), StiffMatrix % Rows(m+1)-1
            MultiplierValues(i) = MultiplierValues(i) - &
                scl * StiffMatrix % Values(j) * Solution(StiffMatrix % Cols(j))
          END DO
        END DO

        DEALLOCATE( UseIPerm, UseDiag )
      ELSE
        i = StiffMatrix % NumberOfRows
        j=0
        IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberOfRows
        IF(ASSOCIATED(AddMatrix)) &
            j=j+MAX(0,AddMatrix % NumberOfRows - StiffMatrix % NumberOFRows)

        Relax = ListGetCReal( Params,'Lagrange Multiplier Relaxation Factor', Found )
        IF( ResidualMode .AND. nIter > 1 ) THEN
          IF( Found ) THEN          
            MultiplierValues(1:j) = MultiplierValues(1:j) + &
                Relax * CollectionSolution(i+1:i+j)
          ELSE
            MultiplierValues(1:j) = MultiplierValues(1:j) + CollectionSolution(i+1:i+j)
          END IF
        ELSE
          IF( Found ) THEN          
            MultiplierValues(1:j) = (1-Relax) * MultiplierValues(1:j) + &
                Relax * CollectionSolution(i+1:i+j)
          ELSE       
            MultiplierValues(1:j) = CollectionSolution(i+1:i+j)
          END IF
        END IF
                
      END IF      
    END IF

      
!------------------------------------------------------------------------------

    IF( SkipConstraints ) THEN
      CALL ListAddLogical( Params,'Skip Advance Nonlinear iter',.FALSE.)
      CALL ListAddLogical( Params,'Skip Compute Nonlinear Change',.FALSE.)
      CALL ComputeChange(Solver,.FALSE.,StiffMatrix % NumberOfRows,Matrix=StiffMatrix,Rhs=ForceVector)
    END IF
        
    DEALLOCATE(CollectionSolution)
    CollectionMatrix % ConstraintMatrix => NULL()
    StiffMatrix % CollectionMatrix => CollectionMatrix

    ParEnv => ParEnvSave
    
    CALL Info( Caller, 'All done', Level=10 )
CONTAINS


!------------------------------------------------------------------------------
   SUBROUTINE CopyMassDampValues(A,B,C)
!------------------------------------------------------------------------------
     TYPE(Matrix_t)  :: A,B,C
!------------------------------------------------------------------------------
     INTEGER :: i,j,n
!------------------------------------------------------------------------------
     REAL(KIND=dp), POINTER :: svals(:)
!------------------------------------------------------------------------------
     n = SIZE(A % Values)

     IF(ASSOCIATED(B % MassValues) .OR. ASSOCIATED(C % MassValues)) THEN
       IF(.NOT.ASSOCIATED(A % MassValues)) THEN
         ALLOCATE(A % MassValues(n)); A % MassValues = 0._dp
       END IF
     END IF

     IF(ASSOCIATED(B % DampValues) .OR. ASSOCIATED(C % DampValues)) THEN
       IF(.NOT.ASSOCIATED(A % DampValues)) THEN
         ALLOCATE(A % DampValues(n)); A % DampValues = 0._dp
       END IF
     END IF

     svals => A % Values

     IF(ASSOCIATED(B % MassValues)) THEN
       A % Values => A % MassValues
       DO i=B % NumberOfRows,1,-1
         DO j=B % Rows(i+1)-1,B % Rows(i),-1
           CALL AddToMatrixElement( A, i, B % Cols(j), B % Values(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(C % MassValues)) THEN
       A % Values => A % MassValues
       DO i=C % NumberOfRows,1,-1
         DO j=C % Rows(i+1)-1,C % Rows(i),-1
           CALL AddToMatrixElement( A, i, C % Cols(j), C % Values(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(B % DampValues)) THEN
       A % Values => A % DampValues
       DO i=B % NumberOfRows,1,-1
         DO j=B % Rows(i+1)-1,B % Rows(i),-1
           CALL AddToMatrixElement( A, i, B % Cols(j), B % DampValues(j) )
         END DO
       END DO
     END IF
 
     IF(ASSOCIATED(C % DampValues)) THEN
       A % Values => A % DampValues
       DO i=C % NumberOfRows,1,-1
         DO j=C % Rows(i+1)-1,C % Rows(i),-1
           CALL AddToMatrixElement( A, i, C % Cols(j), C % DampValues(j) )
         END DO
       END DO
     END IF

     A % Values => svals
!------------------------------------------------------------------------------
   END SUBROUTINE CopyMassDampValues
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE SolveWithLinearRestriction
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Get the node from on which the controlled value should be set. 
!------------------------------------------------------------------------------
  FUNCTION GetControlNode(Mesh,Perm,Params,iControl) RESULT ( ControlNode )
    TYPE(Mesh_t) :: Mesh
    INTEGER, POINTER :: Perm(:)
    TYPE(ValueList_t), POINTER :: Params    
    INTEGER :: iControl
    INTEGER :: ControlNode

    INTEGER :: i,m
    REAL(KIND=dp) :: Coord(3), MinDist
    REAL(KIND=dp), POINTER :: RealWork(:,:)
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'GetControlNode'

    str = 'Control Node Index '//I2S(iControl)                
    ControlNode = ListGetInteger( Params,str,Found )
    IF(.NOT. Found .AND. iControl == 1 ) THEN
      str = 'Control Node Index'
      ControlNode = ListGetInteger( Params,str,Found )    
    END IF
   
    IF(.NOT. Found ) THEN        
      ControlNode = -1

      Coord = 0.0_dp
      str = 'Control Node Coordinates'
      RealWork => ListGetConstRealArray( Params,str,Found )           
      IF(Found) THEN
        i = iControl
      ELSE
        str = TRIM(str)//' '//I2S(iControl)        
        RealWork => ListGetConstRealArray( Params,str,Found )                         
        i = 1
      END IF
      
      IF( Found ) THEN
        IF(SIZE(RealWork,2)==1) THEN
          m = SIZE(RealWork,1)
          Coord(1:m) = RealWork(1:m,i)
        ELSE
          m = SIZE(RealWork,2)
          Coord(1:m) = RealWork(i,1:m)
        END IF
        
        CALL FindClosestNode(Mesh,Coord,MinDist,ControlNode,ParEnv % PEs>1,Perm=Perm)
        CALL Info(Caller,'Control Node located to index: '//I2S(ControlNode),Level=6)

        ! Add the index for future rounds since it takes time to make the search every time!
        str = 'Control Node Index '//I2S(iControl)  
        CALL ListAddInteger( Params, str, ControlNode )
      END IF
    END IF

  END FUNCTION GetControlNode


  FUNCTION GetControlValue(Mesh,Params,iControl,Var,dof) RESULT ( val )
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params  
    INTEGER :: iControl
    TYPE(Variable_t), POINTER, OPTIONAL :: Var
    INTEGER, OPTIONAL :: dof
    REAL(KIND=dp) :: val
    
    TYPE(Variable_t), POINTER :: pVar
    INTEGER :: i,j
    INTEGER :: dof0
    REAL(KIND=dp) :: val0
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: str, varname
    CHARACTER(*), PARAMETER :: Caller = 'GetControlValue'

    IF(.NOT. ASSOCIATED(Mesh)) CALL Fatal(Caller,'Mesh not associated!')
    IF(.NOT. ASSOCIATED(Params)) CALL Fatal(Caller,'Params not associated!')
    
    str = 'Control Variable'
    varname = ListGetString( Params, str, Found )
    IF(.NOT. Found ) THEN
      str = 'Control Variable '//I2S(iControl)
      varname = ListGetString( Params, str, Found )
    END IF
    IF(.NOT. Found ) THEN
      IF( PRESENT(Var) ) THEN
        pVar => Var
      ELSE
        CALL Fatal(Caller,'Could not find keyword for: '//TRIM(str))
      END IF
    ELSE
      pVar => VariableGet( Mesh % Variables, varname )
      IF(.NOT. ASSOCIATED(pVar) ) THEN
        CALL Fatal(Caller,'Could not find control variable: '//TRIM(varname))
      END IF
    END IF

    i = GetControlNode(Mesh,pVar % Perm,Params,iControl)
    IF(i==-1) THEN
      CALL Fatal(Caller,'Could not find control node!')
    END IF
        
    dof0 = 1   
    IF(PRESENT(dof)) THEN
      dof0 = dof
    ELSE IF( pVar % Dofs > 1) THEN
      dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
    END IF

    IF( i == 0 ) THEN
      val = -HUGE(val)
    ELSE
      j = pVar % dofs*(pVar % Perm(i)-1)+dof0
      val = pVar % Values(j) 
    END IF

    val = ParallelReduction(val,2)
    
    str = 'Control Target Value'        
    val0 = ListGetCReal( Params, str, Found )
    IF(.NOT. Found ) THEN
      str = 'Control Target Value '//I2S(iControl)
      val0 = ListGetCReal( Params, str, Found )
    END IF    
    val = val - val0 
  
    !PRINT *,'Control value:',val,val0,i,j

  END FUNCTION GetControlValue
    

  SUBROUTINE ApplyExplicitControl(Solver)
    TYPE(Solver_t) :: Solver

    INTEGER :: i, n
    REAL(KIND=dp), POINTER :: Fvec(:)
    LOGICAL :: Found
    TYPE(Variable_t), POINTER :: FVar
    
    n = ListGetInteger(Solver % Values,'Number Of Controls',Found )
    IF( n == 0 ) THEN
      CALL Warn('ApplyExplicitControl','Explicit control points requested but no controls available!')
      RETURN
    END IF
    
    FVar => VariableGet( Solver % Mesh % Variables,'cpar' )
    IF(.NOT. ASSOCIATED( FVar ) ) THEN
      CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver,'cpar',n,Global=.TRUE.)
      FVar => VariableGet( Solver % Mesh % Variables,'cpar' )
    END IF

    Fvec => FVar % Values    
    DO i=1,n
      Fvec(i) = GetControlValue(Solver % Mesh,Solver % Values,i,Solver % Variable) 
    END DO

    !PRINT *,'Control values:',Fvec
       
  END SUBROUTINE ApplyExplicitControl

  
!------------------------------------------------------------------------------
!> Given the operation point and an additional r.h.s. source vector find the
!> amplitude for the latter one such that the control problem is resolved.
!> We can request a field value at given point, for example. This tries to
!> mimic some ideas of the "Smart Heater Control" of "HeatSolver" available
!> long ago. This would hopefully be applicable to wider set of modules. 
!------------------------------------------------------------------------------
  SUBROUTINE ControlLinearSystem(Solver,PreSolve)
    TYPE(Solver_t) :: Solver
    LOGICAL :: PreSolve

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: x0(:),b(:),BulkRhsSave(:),dr(:),r0(:),dy(:),y0(:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: dx(:),f(:,:)
    INTEGER, POINTER :: Perm(:)
    INTEGER :: dofs, i, j, nsize, ControlNode, dof0, nControl, iControl,jControl
    REAL(KIND=dp) :: Nrm, val, cand, mincand, Relax
    LOGICAL :: GotF, Found, UseLoads, ExtremumMode, DiagControl    
    REAL(KIND=dp), ALLOCATABLE :: cAmp(:), cTarget(:), cVal(:), dc(:), cSens(:,:)
    INTEGER, ALLOCATABLE :: cDof(:)
    
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'ControlLinearSystem'

    
    SAVE f, cAmp, cTarget, cSens, cVal, dc, cDof

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal(Caller,'Controlling of source terms implemented only in serial!')
    END IF
    
    Params => Solver % Values
    Mesh => Solver % Mesh 
    A => Solver % Matrix
    Var => Solver % Variable    
    b => A % RHS
    x0 => Var % Values
    dofs = Var % Dofs
    Perm => Var % Perm
    nsize = SIZE(x0)


    nControl = ListGetInteger(Params,'Number of Controls',Found ) 
    IF(.NOT. Found ) nControl = 1


    IF( PreSolve ) THEN
      CALL Info(Caller,'Applying controlled sources',Level=7)     
      ALLOCATE(f(nsize,nControl),cAmp(nControl),cTarget(nControl),cVal(nControl),&
          dc(nControl),cSens(nControl,nControl),cDof(nControl))      
      cAmp = 0.0_dp; cTarget = 0.0_dp; cVal = 0.0_dp; dc = 0.0_dp; cSens = 0.0_dp; cDof = 0
      
      f = 0.0_dp
      DO iControl = 1, nControl
       ! This is inherited from previous control iterations.
       str = 'Control Amplitude'
       IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)
       cAmp(iControl) = ListGetCReal( Params, str, Found )

       IF(.NOT. Found ) THEN
         str = 'Initial Control Amplitude'
         IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)
         cAmp(iControl) = ListGetCReal( Params, str, Found )          
        END IF
      END DO
      
      DO iControl = 1, Ncontrol            
        ! Default name for controlled source term
        str = TRIM(Var % Name)//' Control'
        IF(Ncontrol>1) str = TRIM(str)//' '//I2S(iControl)
              
        ! We need to add the control source here in order to be able to use
        ! standard means for convergence monitoring. 
        CALL Info(Caller,'Computing source term for: '//TRIM(str),Level=7)
        CALL SetNodalSources( CurrentModel,Mesh,str,dofs, Perm, GotF, f(:,iControl) )

        ! The additional source needs to be nullified for Dirichlet conditions
        IF( ALLOCATED( A % ConstrainedDOF ) ) THEN
          WHERE( A % ConstrainedDOF ) f(:,iControl) = 0.0_dp
        END IF

        IF(InfoActive(10)) THEN
          DO i=1,dofs
            PRINT *,'ranges b:',i,MINVAL(b(i::dofs)),MAXVAL(b(i::dofs)),SUM(b(i::dofs))
            PRINT *,'ranges f:',i,MINVAL(f(i::dofs,iControl)),&
                MAXVAL(f(i::dofs,iControl)),SUM(f(i::dofs,iControl))
          END DO
        END  IF
       
        IF( ABS(cAmp(iControl)) > 1.0e-20 ) THEN
          b(1:nsize) = b(1:nsize) + cAmp(iControl) * f(1:nsize,iControl)
        END IF
      END DO
    END IF

      
    IF(.NOT. PreSolve ) THEN
      CALL Info(Caller,'Dertermining source term amplitude',Level=7)     
      
      CALL ListPushNameSpace('control:')
      CALL ListAddLogical( Params,'control: Skip Compute Nonlinear Change',.TRUE.)
      CALL ListAddLogical( Params,'control: Skip Advance Nonlinear iter',.TRUE.)

      ALLOCATE(dx(nsize))      
      UseLoads = ListGetLogical( Params,'Control Use Loads', Found )
      IF(UseLoads) THEN        
        ALLOCATE(r0(nsize),dr(nsize))      
      END IF

      DiagControl = ListGetLogical( Params,'Control Diagonal', Found )
      
      dof0 = 1
      IF( dofs > 1) THEN
        dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
      END IF
            
      ! Get the target values for control
      DO iControl = 1, Ncontrol            
        str = 'Control Target Value'
        IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
        val = ListGetCReal( Params,str,UnfoundFatal=.TRUE.)
        cTarget(iControl) = val

        i = GetControlNode(Mesh,Perm,Params,iControl) 

        IF(i>0) THEN
          i = dofs*(Perm(i)-1)+dof0
          cDof(iControl) = i
        END IF
      END DO
      
      ! The possibility to use control for extremum temperature is here included.
      ExtremumMode = .FALSE.
      IF( ANY(cDof==0) ) THEN
        IF( Ncontrol == 1 ) THEN
          ExtremumMode = .TRUE.
        ELSE
          CALL Fatal(Caller,'Extremum control cannot be used with '//I2S(Ncontrol)//' controls!')
        END IF
      END IF
      
      
      DO iControl = 1, Ncontrol            
        ! We already know the sources, now compute their affect
        dx = 0.0_dp
        CALL SolveSystem(A,ParMatrix,f(:,iControl),dx,Nrm,dofs,Solver)
        
        ! We use either solution or reaction force for (y0,dy) so that we can
        ! generalize the control procedures for both. 
        IF( UseLoads ) THEN
          ! Nodal loads with the base case
          IF(iControl==1) THEN
            CALL CalculateLoads( Solver, A, x0, dofs, .TRUE., NodalValues = r0 ) 
          END IF

          ! We we use loads then compute the effect of the controlled source to the
          ! reaction force. Hence some hassle with the temporal pointers.          
          BulkRhsSave => A % BulkRhs
          A % BulkRhs => f(:,iControl)
          CALL CalculateLoads( Solver, A, dx, dofs, .TRUE., NodalValues = dr ) 
          A % BulkRhs => BulkRhsSave
          y0 => r0
          dy => dr
        ELSE
          y0 => x0
          dy => dx
        END IF
       
        val = cTarget(iControl) 
        ControlNode = cDof(iControl)

        IF( ExtremumMode ) THEN
          ! We basically do tuning here already but for the sake of uniformity lets just
          ! register the sensitivity and current value. 
          mincand = HUGE(mincand)
          DO i=1,nsize
            j = dofs*(i-1)+dof0          
            IF(ABS(dy(j)) < TINY(dy(j))) CYCLE
            cand = (val-y0(j))/dy(j)
            IF( ABS(cand) < ABS(mincand) ) THEN
              mincand = cand
              cSens(1,1) = dy(j)
              cVal(iControl) = y0(j)
            END IF
          END DO
          CALL Info(Caller,'Extremum value is easiest found in dof: '//I2S(ControlNode),Level=7)    
        ELSE                       
          DO jControl=1,nControl           
            IF(DiagControl .AND. jControl /= iControl) CYCLE
            cSens(jControl,iControl) = dy(cDof(jControl))
          END DO
          cVal(iControl) = y0(cDof(iControl))
        END IF
      END DO

      IF( InfoActive(20) ) THEN
        PRINT *,'cVal:',cVal
        PRINT *,'cTarget:',cTarget
        
        DO i=1,NControl
          PRINT *,'Sens',i,':',cSens(i,:)
        END DO
      END IF
                  
      ! Here we solve the control equation without any assumption of diagonal dominance etc. 
      dc = cTarget - cVal      
      CALL LuSolve(nControl,cSens,dc)

      Relax = ListGetCReal( Params,'Control Relaxation Factor', Found ) 
      IF( Found ) dc = Relax * dc
      

      DO iControl = 1, Ncontrol                    
        str = 'Control Amplitude'
        IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
        cAmp(iControl) = ListGetCReal( Params,str,Found)
        
        cAmp(iControl) = cAmp(iControl) + dc(iControl)
        CALL ListAddConstReal( Params, str, cAmp(iControl) )
        
        ! Apply control, this always to the solution - not to load
        x0(1:nsize) = x0(1:nsize) + dc(iControl) * dx(1:nsize)

        WRITE(Message,'(A,ES15.6)') 'Applied '//TRIM(str)//': ',cAmp(iControl)      
        CALL Info(Caller,Message,Level=5)
      END DO
        
      CALL ListPopNamespace()
      
      DEALLOCATE(f,dx,cAmp,cTarget,cVal,dc,cSens,cDof)
      IF(UseLoads) DEALLOCATE(dr,r0)      
    END IF

    CALL Info(Caller,'All done for now',Level=15)
    
  END SUBROUTINE ControlLinearSystem



!------------------------------------------------------------------------------
!> Given the operation point and an additional r.h.s. source vector find the
!> amplitude for the latter one such that the control problem is resolved.
!> We can request a field value at given point, for example. This tries to
!> mimic some ideas of the "Smart Heater Control" of "HeatSolver" available
!> long ago. This would hopefully be applicable to wider set of modules. 
!------------------------------------------------------------------------------
  SUBROUTINE ControlNonlinearSystem(Solver,PreSolve)
    TYPE(Solver_t) :: Solver
    LOGICAL :: PreSolve

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: A    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: x0(:),b(:),dr(:),r0(:),dy(:),y0(:),prevvalues(:),x(:),dx(:,:)
    INTEGER, POINTER :: Perm(:)
    INTEGER :: dofs, i, j, nsize, ControlNode, dof0, nControl, iControl=0,jControl
    REAL(KIND=dp) :: Nrm, val, cand, mincand, Relax, Eps
    LOGICAL :: GotF, Found, UseLoads, ExtremumMode, DiagControl, Multiply    
    REAL(KIND=dp), ALLOCATABLE :: cAmp(:), cTarget(:), cVal(:), dc(:), cSens(:,:)
    INTEGER, ALLOCATABLE :: cDof(:)
    TYPE(Model_t), POINTER :: Model    
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'ControlNonlinearSystem'
    
    SAVE cAmp, cTarget, cSens, cVal, dc, cDof, iControl, &
        UseLoads, DiagControl, ExtremumMode, Eps, &
        dy, dx, x0, r0, dr, prevvalues 

    IF( ParEnv % PEs > 1 ) THEN
      CALL Fatal(Caller,'Controlling of source terms implemented only in serial!')
    END IF

    Model => CurrentModel
    Params => Solver % Values
    Mesh => Solver % Mesh 
    A => Solver % Matrix
    Var => Solver % Variable    
    b => A % RHS
    x => Var % Values
    dofs = Var % Dofs
    Perm => Var % Perm
    nsize = SIZE(x)


    nControl = ListGetInteger(Params,'Number of Controls',Found ) 
    IF(.NOT. Found ) nControl = 1

    Multiply = .TRUE.
    
    IF( PreSolve ) THEN
      IF( iControl == 0 ) THEN
        CALL Info(Caller,'Applying controlled sources',Level=7)     
        nsize = SIZE(x)
        ALLOCATE(x0(nsize),dx(nsize,nControl),prevvalues(nsize),&
            cAmp(nControl),cTarget(nControl),cVal(nControl),&
            dc(nControl),cSens(nControl,nControl),cDof(nControl))      
        cAmp = 1.0_dp; cTarget = 0.0_dp; cVal = 0.0_dp; dc = 0.0_dp; cSens = 0.0_dp; cDof = 0

        ! Save previous values
        prevvalues = x

        UseLoads = ListGetLogical( Params,'Control Use Loads', Found )
        IF( UseLoads ) THEN
          ALLOCATE(r0(nsize),dr(nsize))
        END IF
                  
        dof0 = 1
        IF( dofs > 1) THEN
          dof0 = ListGetInteger( Params,'Control Target Component',UnfoundFatal=.TRUE.)
        END IF

        ! Get the target values for control
        DO jControl = 1, Ncontrol            
          str = 'Control Target Value'
          IF(nControl > 1) str = TRIM(str)//' '//I2S(jControl)        
          val = ListGetCReal( Params,str,UnfoundFatal=.TRUE.)
          cTarget(jControl) = val
          !i = GetControlNode(jControl)

          i = GetControlNode(Mesh,Perm,Params,jControl) 

          IF(i>0) i = dofs*(Perm(i)-1)+dof0
          cDof(jControl) = i 
        END DO

        ! The possibility to use control for extremum temperature is here included.
        ExtremumMode = .FALSE.
        IF( ANY(cDof==0) ) THEN
          IF( Ncontrol == 1 ) THEN
            ExtremumMode = .TRUE.
          ELSE
            CALL Fatal(Caller,'Extremum control cannot be used with multiple controls!')
          END IF
        END IF

        DiagControl = ListGetLogical( Params,'Control Diagonal', Found ) 

        Eps = ListGetCReal( Params,'Control Epsilon',Found )
        IF(.NOT. Found ) Eps = 0.01_dp
      ELSE IF( iControl == 1 ) THEN        
        CALL ListPushNameSpace('control:')
        CALL ListAddLogical( Params,'control: Skip Compute Nonlinear Change',.TRUE.)
        CALL ListAddLogical( Params,'control: Skip Advance Nonlinear iter',.TRUE.)        
      END IF
    END IF
    

    IF(.NOT. PreSolve) THEN
      IF(iControl == 0 ) THEN
        x0 = Var % Values

        IF(UseLoads) THEN        
          ! Reaction force for the base case
          CALL CalculateLoads( Solver, A, x, dofs, .TRUE., NodalValues = r0 )
        END IF
      ELSE
        ! Remove variation of the parameters
        val = 1.0/(1.0_dp + eps)
        CALL ListSetParameters( Model, iControl, val, multiply, Found )            

        dx(:,iControl) = x - x0
        
        IF(UseLoads) THEN
          ! Reaction force for the variation
          y0 => r0
          CALL CalculateLoads( Solver, A, x, dofs, .TRUE., NodalValues = dr )
          dr = dr - r0
          dy => dr
        ELSE
          y0 => x0
          dy => dx(:,iControl)
        END IF

        val = cTarget(iControl) 
        ControlNode = cDof(iControl)

        IF( ExtremumMode ) THEN
          ! We basically do tuning here already but for the sake of uniformity lets just
          ! register the sensitivity and current value. 
          mincand = HUGE(mincand)
          DO i=1,nsize
            j = dofs*(i-1)+dof0          
            IF(ABS(dy(j)) < TINY(dy(j))) CYCLE
            cand = (val-y0(j))/dy(j)
            IF( ABS(cand) < ABS(mincand) ) THEN
              mincand = cand
              cSens(1,1) = dy(j) / eps
              cVal(iControl) = y0(j)
            END IF
          END DO
          CALL Info(Caller,'Extremum value is easiest found in dof: '//I2S(ControlNode),Level=7)    
        ELSE                       
          DO jControl=1,nControl           
            IF(DiagControl .AND. jControl /= iControl) CYCLE
            cSens(jControl,iControl) = dy(cDof(jControl)) / eps
          END DO
          cVal(iControl) = y0(cDof(iControl))
        END IF
      END IF
                
      IF(iControl == nControl ) THEN
        CALL Info(Caller,'Dertermining source term amplitude',Level=7)     
        
        ! Here we solve the control equation without any assumption of diagonal dominance etc. 
        dc = cTarget - cVal      
        CALL LuSolve(nControl,cSens,dc)

        IF( InfoActive(20) ) THEN
          PRINT *,'cVal:',cVal
          PRINT *,'cTarget:',cTarget          
          DO i=1,NControl
            PRINT *,'Sens',i,':',cSens(i,:)
          END DO
        END IF
        
        Relax = ListGetCReal( Params,'Control Relaxation Factor', Found ) 
        IF( Found ) dc = Relax * dc        
        
        x = x0
        
        DO jControl = 1, Ncontrol                    
          str = 'Control Amplitude'
          IF(nControl > 1) str = TRIM(str)//' '//I2S(iControl)        
          cAmp(jControl) = ListGetCReal( Params,str,Found)
          IF(.NOT. Found) cAmp(jControl) = 1.0_dp

          val = 1.0_dp + dc(jControl)
          cAmp(jControl) = val * cAmp(jControl) 

          IF( .NOT. multiply ) val = cAmp(jControl)             
          CALL ListSetParameters( Model, jControl, val, multiply, Found )            
            
          CALL ListAddConstReal( Params, str, cAmp(jControl) )
        
          ! Apply control, this always to the solution - not to load
          x = x + dc(jControl) * dx(:,jControl)

          WRITE(Message,'(A,ES15.6)') 'Applied '//TRIM(str)//': ',cAmp(iControl)      
          CALL Info(Caller,Message,Level=5)
        END DO
        
        DEALLOCATE(prevvalues,x0,dx)
        IF(UseLoads) DEALLOCATE(r0,dr)
        DEALLOCATE(cAmp,cTarget,cVal,dc,cSens,cDof)
        
        CALL ListPopNamespace()
        iControl = 0
      ELSE
        iControl = iControl + 1
        ! Add variation from the parameters
        val = (1.0_dp + eps)
        CALL ListSetParameters( Model, iControl, val, multiply, Found )                    

        ! Start from the same base case with the matrix assembly
        x = prevvalues
      END IF

    END IF
      
    CALL Info(Caller,'All done for now',Level=15)
    
  END SUBROUTINE ControlNonlinearSystem

  SUBROUTINE SaveLinearSystem( Solver, Ain, LinSysName, OffsetInd )
!------------------------------------------------------------------------------
    TYPE( Solver_t ) :: Solver
    TYPE(Matrix_t), POINTER, OPTIONAL :: Ain
    CHARACTER(LEN=*), OPTIONAL :: LinSysName
    INTEGER, OPTIONAL :: OffsetInd
!------------------------------------------------------------------------------    
    TYPE(Matrix_t), POINTER :: A
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(:), ALLOCATABLE :: dumpfile, dumpprefix
    INTEGER, POINTER :: Perm(:)
    REAL(KIND=dp), POINTER :: Sol(:)
    INTEGER :: i
    LOGICAL :: SaveMass, SaveDamp, SavePerm, SaveSol, Found , Parallel, CNumbering, SkipZeros, SaveSum, &
        SaveAdios2, SaveStiff
    CHARACTER(*), PARAMETER :: Caller = 'SaveLinearSystem'
!------------------------------------------------------------------------------

    CALL Info(Caller,'Saving linear system',Level=6)

    Parallel = ParEnv % PEs > 1

    Params => Solver % Values
    IF(.NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'Parameter list not associated!')
    END IF

    CNumbering = ListGetLogical(Params, 'Linear System Save Continuous Numbering',Found)

    IF( PRESENT(Ain)) THEN
      A => Ain
    ELSE
      A => Solver % Matrix
    END IF

    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal(Caller,'Matrix not associated!')
    END IF

    SaveSum = ListGetLogical( Params,'Linear System Save Sum',Found)

    SaveMass = ListGetLogical( Params,'Linear System Save Mass',Found)

    SaveDamp = ListGetLogical( Params,'Linear System Save Damp',Found)   

    SkipZeros = ListGetLogical( Params,'Linear System Save Skip Zeros', Found )

    SaveAdios2 = ListGetLogical( Params,'Linear System Save Adios2', Found )

    SaveStiff = ListGetLogical( Params,'Linear System Save Stiff', Found )
    IF(.NOT. Found ) SaveStiff = .TRUE.

    SaveSol = ListGetLogical( Params,'Linear System Save Solution',Found)
    IF( SaveSol ) THEN
      Sol => Solver % Variable % Values
      IF( .NOT. ASSOCIATED( Sol ) ) THEN
        CALL Warn(Caller,'Solution not associated!')
        SaveSol = .FALSE.
      END IF
    END IF

    IF( PRESENT( LinSysName ) ) THEN
      dumpprefix = TRIM(LinSysName) 
    ELSE
      dumpprefix = ListGetString( Params, 'Linear System Save Prefix', Found)
      IF(.NOT. Found ) dumpprefix = 'linsys'
    END IF
          
    IF(.NOT. SaveAdios2) THEN

    dumpfile = TRIM(dumpprefix)//'_a.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
    CALL Info(Caller,'Saving matrix to: '//TRIM(dumpfile),Level=5)
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintMatrix(A,Parallel,Cnumbering,SaveMass=SaveMass,SaveDamp=SaveDamp,SkipZeros=SkipZeros)
    CLOSE(1)

    IF( ASSOCIATED(A % rhs) .OR. SaveSum ) THEN
      dumpfile = TRIM(dumpprefix)//'_b.dat'
      IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
      CALL Info(Caller,'Saving matrix rhs to: '//TRIM(dumpfile),Level=5)
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      CALL PrintRHS(A, Parallel, CNumbering, SaveSum )
      CLOSE(1)
    END IF
      
    SavePerm = ListGetLogical( Params,'Linear System Save Perm',Found)
    IF( SavePerm ) THEN
      Perm => Solver % Variable % Perm
      IF( .NOT. ASSOCIATED( Perm ) ) THEN
        CALL Warn(Caller,'Permuation not associated!')
        SavePerm = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_perm.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
        CALL Info(Caller,'Saving permutation to: '//TRIM(dumpfile),Level=5)
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Perm)
          WRITE(1,'(I0,A,I0)') i,' ',Perm(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF


    IF( SaveSol ) THEN
      dumpfile = TRIM(dumpprefix)//'_sol.dat'
      IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
      CALL Info(Caller,'Saving solution to: '//TRIM(dumpfile),Level=5)
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      DO i=1,SIZE(Sol)
        WRITE(1,'(I0,ES15.6)') i,Sol(i)
      END DO
      CLOSE( 1 )
    END IF
    
    
    dumpfile = TRIM(dumpprefix)//'_sizes.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)
    CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile),Level=5)
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    WRITE(1,*) A % NumberOfRows
    WRITE(1,*) SIZE(A % Values)
    i = 0
    IF( SavePerm ) i = SIZE( Perm ) 
    WRITE(1,*) i        
    WRITE(1,*) MINVAL(A % Cols)
    WRITE(1,*) MAXVAL(A % Cols)
    IF(PRESENT(OffsetInd)) WRITE(1,*) OffsetInd
    CLOSE(1)

    IF(Parallel) THEN
      dumpfile = TRIM(dumpprefix)//'_sizes.dat'
      CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile),Level=6)
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      WRITE(1,*) ParallelReduction(A % ParMatrix % &
                           SplittedMatrix % InsideMatrix % NumberOfRows)
      WRITE(1,*) ParallelReduction(SIZE(A % Values))
      IF( SavePerm ) WRITE(1,*) ParallelReduction(SIZE( Perm ))
      CLOSE(1)
    END IF

    ELSE
      dumpfile = TRIM(dumpprefix)//'.bp'
      CALL SaveLinearSystemAdios2( A, dumpfile, Parallel, CNumbering, SaveMass, SaveDamp, SaveStiff, SaveSol, Sol )
    END IF

    IF( ListGetLogical( Params,'Linear System Save and Stop',Found ) ) THEN
      CALL Info(Caller,'Just saved matrix and stopped!',Level=4)
      STOP EXIT_OK
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SaveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Saves the linear system (CRS matrix, optionally mass/damp matrices, and rhs)
!> using ADIOS2Utils. Global dofs from A % parallelinfo % globaldofs is saved to the output file
!> or if CNumbering is enabled, then the output of ContinuousNumbering A % Gorder is saved.
!> Note! The ADIOS2 file is finalized at the end of the call so multiple timesteps will cause problems.
!------------------------------------------------------------------------------
  SUBROUTINE SaveLinearSystemAdios2( A, dumpfile, Parallel, CNumbering, SaveMass, SaveDamp, SaveStiff, SaveSol, Sol )
!------------------------------------------------------------------------------
#ifdef HAVE_ADIOS2
    USE ADIOS2Utils
#endif
    TYPE(Matrix_t) :: A
    CHARACTER(LEN=*) :: dumpfile
    LOGICAL :: Parallel, CNumbering, SaveMass, SaveDamp, SaveStiff, SaveSol
    REAL(KIND=dp), POINTER :: Sol(:)
!------------------------------------------------------------------------------
#ifdef HAVE_ADIOS2
    TYPE(AdiosWriter_t) :: Writer
    INTEGER, ALLOCATABLE :: Owner(:)
    INTEGER, ALLOCATABLE :: GlobalRows(:), GlobalCols(:), RowNnz(:)
    INTEGER :: i, j, n, nnz, ierr
    CHARACTER(*), PARAMETER :: Caller = 'SaveLinearSystemAdios2'

    CALL Info(Caller,'Saving linear system using Adios2 to: '//TRIM(dumpfile),Level=5)

    n = A % NumberOfRows
    nnz = A % Rows(n+1) - 1

    ALLOCATE( GlobalRows(n),  )


    IF( Parallel ) THEN
      IF( CNumbering ) THEN
        IF( .NOT. ASSOCIATED( A % Gorder ) ) THEN
          ALLOCATE( A % Gorder(n), Owner(n) )
          CALL ContinuousNumbering( A % ParallelInfo, A % Perm, A % Gorder, Owner )
        END IF
        GlobalRows(1:n) = A % Gorder(1:n)
      ELSE
        GlobalRows(1:n) = A % ParallelInfo % GlobalDOFs(1:n)
      END IF
    ELSE
      DO i=1,n
        GlobalRows(i) = i
      END DO
    END IF

    ierr = Writer % init( dumpfile, array_kind = ADIOS2_ARRAY_GLOBAL )

    CALL Writer % write_data( 'global_dofs', GlobalRows )
    CALL Writer % write_data( 'rows', A % rows )
    CALL Writer % write_data( 'cols', A % cols )

    IF( SaveStiff ) THEN
      CALL Writer % write_data( 'values', A % Values(1:nnz) )
    END IF

    IF( SaveMass ) THEN
      CALL Writer % write_data( 'mass_values', A % MassValues(1:nnz) )
    END IF

    IF( SaveDamp ) THEN
      CALL Writer % write_data( 'damp_values', A % DampValues(1:nnz) )
    END IF

    IF( ASSOCIATED( A % rhs ) ) THEN
      CALL Writer % write_data( 'rhs', A % rhs(1:n) )
    END IF

    IF( SaveSol ) THEN
      CALL Writer % write_data( 'solution', Sol )
    END IF

    ierr = Writer % finalize()

    DEALLOCATE( GlobalRows )
#else
    CALL Fatal('SaveLinearSystemAdios2','Elmer was not compiled with ADIOS2 support!')
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE SaveLinearSystemAdios2
!------------------------------------------------------------------------------

  SUBROUTINE LinearSystemMultiply( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar,CoeffVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff,Coeff2
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(:), ALLOCATABLE :: str, VarName

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values
        
    UpdateRhs = ListGetLogical( Params,'Linear System Multiply Consistent',Found)
    Symmetric = ListGetLogical( Params,'Linear System Multiply Symmetric',Found)

    ! First, Multiply the k:th piece of the r.h.s. vector if requested
    !-----------------------------------------------------------
    DO k=1,Dofs
      Mode = 0
      
      str = 'Linear System Rhs Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the rhs with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the rhs with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      IF( Mode == 0 ) THEN
        str = 'Linear System Rhs Variable'
        VarName = ListGetString( Params,str,Found ) 
        NULLIFY( CoeffVar ) 
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          VarName = ListGetString( Params,str//' '//I2S(k),Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the rhs with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )

          !PRINT *,'Range:',Mode,MINVAL(CoeffVar % Values),MAXVAL(CoeffVar % Values)
        END IF
      END IF
      IF( Mode == 0 ) CYCLE
 
      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          Rhs = Coeff * Rhs
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
          Rhs( jk ) = Coeff * Rhs( jk )
        END DO
      END IF
    END DO
    ! End of r.h.s. multiplication

    ! Secondly, Multiply the kl block of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      DO l=1,Dofs
        Mode = 0
        str = 'Linear System Matrix Factor'
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 1
          WRITE( Message,'(A,ES12.3)') 'Multiplying the matrix with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        ELSE
          Coeff = ListGetCReal( Params, str//' '//I2S(k)//I2S(l), Found )
          IF( Found ) THEN
            Mode = 2
            WRITE( Message,'(A,I0,I0,A,ES12.3)') 'Multiplying block (',k,l,') of the matrix with ',Coeff
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF
        IF( Mode == 0 ) THEN
          str = 'Linear System Matrix Variable'
          VarName = ListGetString( Params,str,Found )
          NULLIFY( CoeffVar ) 
          IF( Found ) THEN
            CoeffVar => VariableGet( Mesh % Variables, str )                                    
          ELSE
            VarName = ListGetString( Params,str//' '//I2S(k)//I2S(l),Found )
            IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
          END IF
          IF( ASSOCIATED( CoeffVar ) ) THEN
            IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
              CALL Fatal('LinearSystemMultiply','Permutations should be the same')
            END IF
            Mode = 3
            WRITE( Message,'(A,I0,I0,A)') 'Multiplying block (',k,l,') of the matrix with > '//TRIM(VarName)//' <'
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF

        IF( Mode == 0 ) CYCLE

        IF( Mode == 1 ) THEN
          IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
            Values = Coeff * Values
          END IF
        ELSE 
          DO j=1,SIZE( Perm ) 
            jk = Dofs*(j-1)+k
            IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
            DO i=Rows(jk),Rows(jk+1)-1 
              IF( MODULO( Cols(i), Dofs ) == MODULO( l, Dofs ) ) THEN
                IF( Mode == 3 .AND. Symmetric ) THEN          
                  j2 = (Cols(i)-1)/Dofs + 1 
                  Coeff2 = CoeffVar % Values(j2)
                  Values( i ) = SQRT( Coeff * Coeff2 ) * Values( i )
                ELSE
                  Values( i ) = Coeff * Values( i )
                END IF
              END IF
            END DO
          END DO
        END IF
      END DO
      IF( Mode == 1 ) EXIT
    END DO
    ! end of matrix multiplication


    ! Finally, Multiply the diagonal of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the diagonal with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the matrix diagonal with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )          
        END IF
      END IF

      IF( Mode == 0 ) THEN
        str = 'Linear System Diagonal Variable'
        VarName = ListGetString( Params,str,Found )
        NULLIFY( CoeffVar )
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          VarName = ListGetString( Params,str//' '//I2S(k),Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the diagonal with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE

      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          IF( UpdateRhs ) Rhs = Rhs + ( Coeff - 1 ) * Values( A % Diag ) * ThisVar % Values
          Values( A % Diag ) = Coeff * Values( A % Diag )
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)

          IF( UpdateRhs ) Rhs( jk ) = Rhs( jk ) + (Coeff - 1) * Values(A % Diag(jk)) * ThisVar % Values(jk)
          Values( A % Diag(jk) ) = Coeff * Values( A % Diag(jk) )
        END DO
      END IF
    END DO
    ! end of diagonal multiplication


  END SUBROUTINE LinearSystemMultiply






!---------------------------------------------------------------------------------
!> Set the diagonal entry to given minimum.
!----------------------------------------------------------------------------------
  SUBROUTINE LinearSystemMinDiagonal( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(:), ALLOCATABLE :: str
    INTEGER :: NoSet
    REAL(KIND=dp) :: DiagSum, val, DiagMax

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values

    ! Set the minimum value for each component, only nodel dofs considered
    !---------------------------------------------------------------------
    NoSet = 0
    DiagMax = 0.0_dp
    DiagSum = 0.0_dp
    n = MAXVAL( Perm ( 1:Mesh % NumberOfNodes ) )

    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Min'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Setting minimum of the diagonal to ',Coeff
        CALL Info('LinearSystemMinDiagonal',Message, Level=6 )
      ELSE
        Coeff = ListGetCReal( Params, str//' '//I2S(k), Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Setting minimum of diagonal component ',k,' to ',Coeff
          CALL Info('LinearSystemMinDiagonal',Message, Level=6 )          
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE
      
      DO j=1,n
        jk = Dofs*(j-1)+k
        l = A % Diag(jk) 
        IF( l == 0 ) CYCLE

        ! Only add the diagonal to the owned dof
        IF( ParEnv % PEs > 1 ) THEN
          IF( A % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
        END IF

        val = ABS( Values( l ) )
        DiagSum = DiagSum + val
        DiagMax = MAX( DiagMax, val )
        IF( val < Coeff ) THEN
          Values( A % Diag(jk) ) = Coeff
          NoSet = NoSet + 1
        END IF
      END DO
    END DO

    CALL Info('LinearSystemMinDiagonal',&
        'Number of diagonal values set to minimum: '//I2S(NoSet),Level=5)
    WRITE( Message,'(A,ES12.3)') 'Average abs(diagonal) entry: ',DiagSum / n
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )

    WRITE( Message,'(A,ES12.3)') 'Maximum abs(diagonal) entry: ',DiagMax
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )


  END SUBROUTINE LinearSystemMinDiagonal





  !----------------------------------------------------------------------
  !> Make the high-order flux corrected transport (FCT) correction after 
  !> the low order approximation has been solved. 
  !
  !> For more information see, for example, 
  !> Dmitri Kuzmin (2008): "Explicit and implicit FEM-FCT algorithms with flux linearization"
  !----------------------------------------------------------------------
  SUBROUTINE FCT_Correction( Solver )

    TYPE(Solver_t) :: Solver

    TYPE(Valuelist_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n,i,j,k,k2
    INTEGER, POINTER :: Rows(:),Cols(:),Perm(:)
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: BulkValues(:),u(:),M_L(:),udot(:), &
        pp(:),pm(:),qp(:),qm(:),corr(:),ku(:),ulow(:)
    REAL(KIND=dp), ALLOCATABLE :: mc_udot(:), fct_u(:)
    REAL(KIND=dp), POINTER CONTIG :: M_C(:), SaveValues(:)
    REAL(KIND=dp) :: rsum, Norm,m_ij,k_ij,du,d_ij,f_ij,c_ij,Ceps,CorrCoeff,&
        rmi,rmj,rpi,rpj,dt
    TYPE(Variable_t), POINTER :: Var, Variables
    LOGICAL :: Found, Symmetric, SaveFields, SkipCorrection
    CHARACTER(:), ALLOCATABLE :: VarName, TmpVarName

    REAL(KIND=dp), POINTER :: mmc(:), mmc_h(:), fct_d(:)
    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: ActiveNodes(:)

    Params => Solver % Values

    SkipCorrection = ListGetLogical( Params,'FCT Correction Skip',Found )
    Symmetric = ListGetLogical( Params,'FCT Correction Symmetric',Found )
    SaveFields = ListGetLogical( Params,'FCT Correction Save',Found )

    IF( SkipCorrection .AND. .NOT. SaveFields ) THEN
      CALL Info('FCT_Correction','Skipping the computation of FCT correction',Level=5)
    END IF

    CALL Info('FCT_Correction','Computing corrector for the low order solution',Level=5)

    ! PRINT *,'FCT Norm Before Correction:',SQRT( SUM( Solver % Variable % Values**2) )
 
    Mesh => Solver % Mesh
    Variables => Solver % Mesh % Variables
 
    ! Set pointers 
    A => Solver % Matrix
    n = A % NumberOfRows
    Rows => A % Rows
    Cols => A % Cols

    BulkValues => A % BulkValues
    M_C => A % MassValues
    Perm => Solver % Variable % Perm

    M_L => A % MassValuesLumped 
    IF (ParEnv % PEs>1) CALL ParallelSumVector(A,M_L)
    
    Var => VariableGet( Variables,'timestep size')
    dt = Var % Values(1) 

    ! low order solution at the start, high order in the end
    u => Solver % Variable % Values
    VarName = GetVarName(Solver % Variable)

    ! Here a bunch of vectors are stored for visualization and debugging purposes
    !----------------------------------------------------------------------------

    ! low order solution is the solution without corrections
    ! This is created and saved only if requested
    !---------------------------------------------------------------------------
    IF( SaveFields ) THEN
      TmpVarName = TRIM( VarName )//' fctlow'    
      Var => VariableGet( Variables, TmpVarName )
      IF( .NOT. ASSOCIATED(Var) ) THEN
        CALL VariableAddVector( Variables, Mesh, Solver,&
            TmpVarName, Perm = Perm, Output = SaveFields )
        Var => VariableGet( Variables, TmpVarName )
      END IF
      ulow => Var % Values
      ulow = u
    END IF

    ! Create auxiliary vectors for the fct algorithm
    !---------------------------------------------------------------------------

    ! r.h.s. term
    TmpVarName = TRIM( VarName )//' fctku'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    ku => Var % Values

    ! time derivative from lower order analysis
    TmpVarName = TRIM( VarName )//' fctudot'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    udot => Var % Values

    ! Fields related to flux limiters
    TmpVarName = TRIM( VarName )//' fctpp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pp => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctpm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pm => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctqp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qp => Var % Values

    TmpVarName = TRIM( VarName )//' fctqm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qm => Var % Values

    TmpVarName = TRIM( VarName )//' fctmm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    Var % Values = M_L

    ! higher order correction 
    TmpVarName = TRIM( VarName )//' fctcorr'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    corr => Var % Values


    ! 1) Compute the nodal time derivatives
    ! M_C*udot=K*ulow  (M_C is the consistent mass matrix)
    !----------------------------------------------------------------------
    CALL Info('FCT_Correction','Compute nodal time derivatives',Level=10)
    ! Compute: ku = K*ulow
#if 0
    DO i=1,n
      rsum = 0.0_dp
      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        K_ij = BulkValues(k)
        rsum = rsum + K_ij * u(j) 
      END DO
      ku(i) = rsum
    END DO
    ! Solve the linear system for udot
    ! The stiffness matrix is momentarily replaced by the consistent mass matrix M_C
    ! Also the namespace is replaced to 'fct:' so that different strategies may 
    ! be applied to the mass matrix solution.
    CALL ListPushNameSpace('fct:')
    CALL ListAddLogical( Params,'fct: Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'fct: Skip Advance Nonlinear iter',.TRUE.)
    SaveValues => A % Values
    A % Values => M_C
    CALL SolveLinearSystem( A, ku, udot, Norm, 1, Solver )
    A % Values => SaveValues
    CALL ListPopNamespace()
#else

  BLOCK
    REAL(KIND=dp), ALLOCATABLE :: TmpRhsVec(:), TmpXVec(:)

    SaveValues => A % Values

    IF (Parenv % PEs>1) THEN
      ALLOCATE(TmpRHSVec(n), TmpXVec(n))
      TmpxVec = 0._dp; tmpRHSVec = 0._dp

      A % Values => BulkValues
      CALL ParallelInitSolve(A,TmpXVec,TmpRhsVec,u)
      CALL ParallelVector(A,TmpRhsvec,u)

      CALL ParallelMatrixVector(A,TmpRhsvec,TmpXVec)

      CALL PartitionVector(A,Ku,TmpXVec)
      DEALLOCATE(TmpRhsVec, TmpXVec)
    ELSE
      DO i=1,n
        rsum = 0._dp
        DO k=Rows(i),Rows(i+1)-1
          j = Cols(k)
          K_ij = BulkValues(k)
          rsum = rsum + K_ij * u(j) 
        END DO
        ku(i) = rsum
      END DO
    END IF

    CALL ListPushNameSpace('fct:')
    CALL ListAddLogical( Params,'fct: Skip Compute Nonlinear Change',.TRUE.)
    CALL ListAddLogical( Params,'fct: Skip Advance Nonlinear iter',.TRUE.)
  
    A % Values => M_C
    udot = 0._dp
    CALL SolveLinearSystem(A,Ku,Udot,Norm,1,Solver)
    CALL ListPopNamespace()

    A % Values => SaveValues
  END BLOCK
#endif

    ! Computation of correction factors (Zalesak's limiter)
    ! Code derived initially from Kuzmin's subroutine   
    !---------------------------------------------------------
    CALL Info('FCT_Correction','Compute correction factors',Level=10)
    pp = 0 
    pm = 0
    qp = 0 
    qm = 0

    IF(ParEnv % PEs>1) THEN
      fct_d => A % FCT_D
      mmc    => A % MassValues
      mmc_h  => A % HaloMassValues

      ALLOCATE(ActiveNodes(n)); activeNodes=.FALSE.
      DO i=1,Solver % NumberOfActiveElements
        Element => Solver % Mesh % Elements(Solver % ActiveElements(i))
        IF ( Element % PartIndex /= ParEnv % MyPE ) CYCLE
        Activenodes(Solver % Variable % Perm(Element % NodeIndexes)) = .TRUE.
      END DO
    ELSE
      fct_d => A % FCT_D
      mmc => A % MassValues
    END IF
    DO i=1,n
      IF (ParEnv % PEs > 1 ) THEN
        IF ( .NOT. ActiveNodes(i) ) CYCLE
      end if

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)

        IF (ParEnv % PEs>1) THEN
          IF ( .NOT.ActiveNodes(j) ) CYCLE
        END IF

        ! Compute the raw antidiffusive fluxes
        ! f_ij=m_ij*[udot(i)-udot(j)]+d_ij*[ulow(i)-ulow(j)]
        !-----------------------------------------------------
        ! d_ij and m_ij are both symmetric
        ! Hence F_ji = -F_ij
           
        f_ij = mmc(k)*(udot(i)-udot(j)) + fct_d(k)*(u(i)-u(j))
        IF ( ParEnv % PEs>1 ) f_ij=f_ij+mmc_h(k)*(udot(i)-udot(j))
        ! Compared to Kuzmin's paper F_ij=-F_ij since d_ij and 
        ! udot have different signs. 
        f_ij = -f_ij

        ! Antidiffusive fluxes to be limited
        du = u(j)-u(i)

        ! Prelimiting of antidiffusive fluxes i.e. du and the flux have different signs
        IF (f_ij*du >= TINY(du)) THEN
          f_ij = 0._dp
        ELSE        
          ! Positive/negative edge contributions
          pp(i) = pp(i) + MAX(0._dp,f_ij)
          pm(i) = pm(i) + MIN(0._dp,f_ij)
        END IF

        ! Maximum/minimum solution increments
        qp(i) = MAX(qp(i),du)
        qm(i) = MIN(qm(i),du)
      END DO
    END DO

    ! Computation of nodal correction factors
    ! DO i=1,n
    !  IF( pp(i) > Ceps ) THEN
    !    rp(i) = MIN( 1.0_dp, M_L(i)*qp(i)/pp(i) )
    !  END IF
    !  IF( pm(i) < -Ceps ) THEN
    !    rm(i) = MIN( 1.0_dp, M_L(i)*qm(i)/pm(i) )
    !  END IF
    ! END DO

    ! Correct the low-order solution
    ! (M_L*ufct)_i=(M_L*ulow)_i+dt*sum(alpha_ij*f_ij)
    !-------------------------------------------------
    ! Symmetric flux limiting
    ! Correction of the right-hand side


    !   IF (ParEnv % PEs>1) THEN
    !     CALL ParallelSumVector(A,pm)
    !     CALL ParallelSumVector(A,pp)
    !     CALL ParallelSumVector(A,qm)
    !     CALL ParallelSumVector(A,qp)
    !   END IF
    
    CorrCoeff = ListGetCReal( Params,'FCT Correction Coefficient',Found )
    IF( .NOT. Found ) CorrCoeff = 1._dp

    Ceps = ListGetCReal( Params,'FCT Correction Epsilon',Found )
    IF(.NOT. Found ) Ceps = EPSILON( Ceps )
    corr = 0._dp
    DO i=1,n
      IF (ParEnv % PEs>1) THEN
        IF( .NOT. ActiveNodes(i)) CYCLE
      END IF

      IF( pp(i) > Ceps ) THEN
        rpi = MIN( 1._dp, M_L(i)*qp(i)/pp(i) )
      ELSE
        rpi = 0._dp
      END IF

      IF( pm(i) < -Ceps ) THEN
        rmi = MIN( 1._dp, M_L(i)*qm(i)/pm(i) )
      ELSE
        rmi = 0._dp
      END IF

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        IF(ParEnv % PEs>1) THEN
          IF(.NOT.ActiveNodes(j)) CYCLE
        END IF

        f_ij = mmc(k)*(udot(i)-udot(j))  + fct_d(k)*(u(i)-u(j))
        IF (ParEnv % PEs>1) f_ij = f_ij + mmc_h(k)*(udot(i)-udot(j))
        f_ij = -f_ij

        IF (f_ij > 0) THEN 
          IF( pm(j) < -Ceps ) THEN
            rmj = MIN( 1.0_dp, M_L(j)*qm(j)/pm(j) )
          ELSE
            rmj = 0._dp
          END IF
          c_ij = MIN(rpi,rmj)
        ELSE 
          IF( pp(j) > Ceps ) THEN
            rpj = MIN( 1._dp, M_L(j)*qp(j)/pp(j) )
          ELSE
            rpj = 0._dp
          END IF
          c_ij = MIN(rmi,rpj)
        END IF
        corr(i) = corr(i) + c_ij * f_ij
      END DO
      corr(i) = CorrCoeff * corr(i) / M_L(i)
    END DO

    IF (ParEnv % PEs>1) THEN
!     CALL ParallelSumVector(A,corr)
      DEALLOCATE(A % HaloValues, A % HaloMassValues)
      A % HaloValues => Null(); A % HaloMassValues => Null()
    END IF

    ! Optionally skip applying the correction, just for debugging purposes
    IF( SkipCorrection ) THEN
      CALL Info('FCT_Correction','Skipping Applying corrector',Level=6)
    ELSE
      CALL Info('FCT_Correction','Applying corrector for the low order solution',Level=10)

      u = u + corr

      ! PRINT *,'FCT Norm After Correction:',SQRT( SUM( Solver % Variable % Values**2) )
    END IF

  END SUBROUTINE FCT_Correction


   SUBROUTINE DetermineSoftLimiter( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterV, LimitVar
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,t,ind,dofs, dof, bf, bc, Upper, Removed, Added, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), &
         ElemLimit(:), ElemActive(:), ElemWrk(:)
     REAL(KIND=dp) :: LimitSign, EqSign, ValEps, LoadEps, ValEps0, LoadEps0, &
         MaxValue, MaxLoad, val
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF, GotActive
     LOGICAL, ALLOCATABLE :: LimitDone(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params, Entity
     LOGICAL, ALLOCATABLE :: InterfaceDof(:)
     INTEGER :: ConservativeAfterIters, NonlinIter, CoupledIter, timeIter, DownStreamDirection
     LOGICAL :: Conservative, ConservativeAdd, ConservativeRemove, RelativeEps, &
         DoAdd, DoRemove, DirectionActive, FirstTime, DownStreamRemove, LimitFreeze, &
         AllActive, AllPassive
     TYPE(Mesh_t), POINTER :: Mesh

     CHARACTER(:), ALLOCATABLE :: Name,LimitName, InitName, ActiveName
     CHARACTER(*), PARAMETER :: Caller = 'DetermineSoftLimiter'

     Model => CurrentModel
     Var => Solver % Variable
     Mesh => Solver % Mesh
     
     ! The variable to be constrained by the soft limiters
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     Params => Solver % Values

     ! Check the iterations counts and determine whether this is the first 
     ! time with this solver. 
     !------------------------------------------------------------------------
     iterV => VariableGet( Mesh % Variables,'nonlin iter')
     IF( ASSOCIATED( iterV ) ) THEN
       NonlinIter =  NINT( iterV % Values(1) ) 
     ELSE
       NonlinIter = 1
     END IF

     iterV => VariableGet( Mesh % Variables,'coupled iter')
     IF( ASSOCIATED( iterV ) ) THEN
       CoupledIter = NINT( iterV % Values(1) )
     ELSE
       CoupledIter = 1
     END IF

     iterV => VariableGet( Mesh % Variables,'timestep')
     IF( ASSOCIATED( iterV ) ) THEN
       timeIter = NINT( iterV % Values(1) )
     ELSE
       timeIter = 1
     END IF

     FirstTime = (nonliniter <= 1) .AND. (coupledIter <= 1) .AND. (timeIter == 1)

     ! Always freeze the contact set when going to new timestep since the residual values
     ! on the new timestep are not reliable before a new timestep has been solved. 
     LimitFreeze = (timeIter > 1) .AND. (nonliniter <= 1) .AND. (coupledIter <= 1)

     ! Optionally freeze the contact set for the whole 1st timestep
     IF(timeIter == 1 .AND. .NOT. FirstTime ) THEN
       LimitFreeze = ListGetLogical( Params,'Limiter Passive First Timestep', Found )
     END IF

     ! We can turn optionally the contact set fully active/passive using a global condition.      
     AllActive = .FALSE.
     val = ListGetCReal( Params,'Limiter Global Active Condition',Found )
     IF(Found) AllActive = (val > 0.0_dp) 

     AllPassive = .FALSE.
     val = ListGetCReal( Params,'Limiter Global Passive Condition',Found )
     IF(Found) AllPassive = (val > 0.0_dp) 

     IF(AllActive .AND. AllPassive) THEN
       CALL Fatal(Caller,'Limiter cannot be both AllActive and AllPassive!')
     END IF

     
     IF( FirstTime ) THEN
       CALL Info(Caller,'Initializing soft limiter for solver',Level=7)
     END IF
     IF( LimitFreeze ) THEN
       CALL Info(Caller,'Keeping soft limiter fixed during this cycle!',Level=7)
     END IF
     
     ! Determine variable for computing the contact load used to determine the 
     ! soft limit set.
     !------------------------------------------------------------------------
     CALL Info(Caller,'Determining soft limiter problems',Level=8)
     LoadVar => VariableGet( Model % Variables, &
         GetVarName(Var) // ' Contact Load',ThisOnly = .TRUE. )
     CALL CalculateLoads( Solver, Solver % Matrix, Var % Values, Var % DOFs, .FALSE., LoadVar ) 

     IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
       CALL Fatal(Caller, &
           'No Loads associated with variable '//GetVarName(Var) )
       RETURN
     END IF
     LoadValues => LoadVar % Values


     
     ConservativeAdd = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',Conservative ) 
     IF( Conservative ) THEN
       ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeRemove = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',Found )      
     IF( Found ) THEN
       Conservative = .TRUE.  
       ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF

     DownStreamRemove = ListGetLogical( Params,'Apply Limiter Remove Downstream',Found)
     IF( DownStreamRemove ) THEN
       CALL Info(Caller,'Removing contact dofs only in downstream',Level=8)      
       ConservativeRemove = .TRUE.
       Conservative = .TRUE.
       DownStreamDirection = ListGetInteger( Params,'Apply Limiter Downstream Direction',Found)
       IF(.NOT. Found ) DownStreamDirection = 1
     END IF
       
     LoadEps0 = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps0 = 1.0e-8_dp
     LoadEps = LoadEps0
     
     ValEps0 = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps0 = 1.0e-8_dp
     ValEps = ValEps0

     RelativeEps = ListGetLogical(Params,'Limiter Relative Tolerance', Found ) 
     
     ! The user may want to toggle the sign for various kinds of equations
     ! The default sign that come from standard formulation of Laplace equation.
     !---------------------------------------------------------------------------       
     IF( ListGetLogical( Params,'Limiter Load Sign Negative',Found) ) THEN
       EqSign = -1.0_dp
     ELSE
       EqSign = 1.0_dp
     END IF

     ! Loop through upper and lower limits     
     !------------------------------------------------------------------------
     DO Upper=0,1

       DirectionActive = .FALSE.

       ! If we have both upper and lower limiter then these logical vectors need to be 
       ! reinitialized for the 2nd sweep.
       IF( ALLOCATED( LimitDone) ) LimitDone = .FALSE.
       IF( ALLOCATED( InterfaceDof ) ) InterfaceDof = .FALSE. 

       ! Upper and lower limits have different sign for testing
       !----------------------------------------------------------------------       
       IF( Upper == 0 ) THEN
         LimitSign = -EqSign
       ELSE
         LimitSign = EqSign
       END IF       
       
       ! Go through the components of the field, if many
       !-------------------------------------------------
       DO DOF = 1,dofs

         Name = TRIM(Var % name)
         IF ( Var % DOFs > 1 ) name = ComponentName(name,DOF)

         ! The keywords for the correct lower or upper limit of the variable
         !------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           LimitName = TRIM(name)//' Lower Limit'           
           InitName = TRIM(name)//' Lower Initial'
           ActiveName = TRIM(name)//' Lower Active'
         ELSE
           LimitName = TRIM(name)//' Upper Limit' 
           InitName = TRIM(name)//' Upper Initial' 
           ActiveName = TRIM(name)//' Upper Active' 
         END IF

         AnyLimitBC = ListCheckPresentAnyBC( Model, LimitName )
         AnyLimitBF = ListCheckPresentAnyBodyForce( Model, LimitName )

         ! If there is no active keyword then there really is nothing to do
         !----------------------------------------------------------------
         IF( .NOT. ( AnyLimitBC .OR. AnyLimitBF ) ) CYCLE
         DirectionActive = .TRUE.
         
         CALL Info(Caller,'Applying limit: '//TRIM(LimitName),Level=8)

         ! OK: Do contact for a particular dof and only upper or lower limit
         !------------------------------------------------------------------------

         ! Define the range of elements for which the limiters are active
         !---------------------------------------------------------------
         ElemFirst = Model % NumberOfBulkElements + 1           
         ElemLast = Model % NumberOfBulkElements 
        
         IF( AnyLimitBF ) ElemFirst = 1
         IF( AnyLimitBC ) ElemLast = Model % NumberOfBulkElements + &
             Model % NumberOfBoundaryElements 

         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemActive(n), ElemWrk(n) )
           LimitDone = .FALSE.
         END IF

         ! Check that active set vectors for limiters exist, otherwise allocate
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           IF(ASSOCIATED(Var % LowerLimitActive)) THEN
             IF(SIZE(Var % LowerLimitActive) /= totsize) THEN
               DEALLOCATE(Var % LowerLimitActive)
               DEALLOCATE(Var % LowerLimit)
             END IF
           END IF
           IF( .NOT. ASSOCIATED(Var % LowerLimitActive ) ) THEN
             ALLOCATE( Var % LowerLimitActive( totsize ) )
             Var % LowerLimitActive = .FALSE.
             ALLOCATE( Var % LowerLimit( totsize ) )
             Var % LowerLimit = -HUGE(val)
           END IF
           LimitActive => Var % LowerLimitActive
         ELSE
           IF(ASSOCIATED(Var % UpperLimitActive)) THEN
             IF(SIZE(Var % UpperLimitActive) /= totsize) THEN
               DEALLOCATE(Var % UpperLimitActive)
               DEALLOCATE(Var % UpperLimit)
             END IF
           END IF
           IF( .NOT. ASSOCIATED( Var % UpperLimitActive ) ) THEN
             ALLOCATE( Var % UpperLimitActive( totsize ) )
             Var % UpperLimitActive = .FALSE.
             ALLOCATE( Var % UpperLimit( totsize ) )
             Var % UpperLimit = HUGE(val)
           END IF
           LimitActive => Var % UpperLimitActive
         END IF
 
         Removed = 0
         Added = 0        
         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemActive(n), ElemWrk(n) )
           LimitDone = .FALSE.
         END IF

         IF( RelativeEps ) THEN         
           IF( dofs == 1 ) THEN
             MaxLoad = MAXVAL( ABS( LoadValues ) )
             MaxValue = MAXVAL( ABS( Var % Values ) ) 
           ELSE
             MaxLoad = MAXVAL( ABS( LoadValues(dof::dofs ) ) )
             MaxValue = MAXVAL( ABS( Var % Values(dof::dofs) ) ) 
           END IF
           MaxLoad = ParallelReduction(MaxLoad)
           MaxValue = ParallelReduction(MaxValue)
           IF( InfoActive(20) ) THEN
             WRITE(Message,'(A,ES12.3)') 'Using load maximum for scaling: ',MaxLoad
             CALL Info(Caller,Message)
             WRITE(Message,'(A,ES12.3)') 'Using value maximum for scaling: ',MaxValue
             CALL Info(Caller,Message)
           END IF
           LoadEps = LoadEps0 * MaxLoad
           ValEps = ValEps0 * MaxValue 
         END IF

         IF( FirstTime ) THEN
           ! In the first time set the initial set 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element

             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE
               IF(Element % BodyId == 0) CYCLE
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE               
               Entity => Model % BodyForces(bf) % Values               
             END IF

             ElemLimit(1:n) = ListGetReal( Entity, LimitName, n, NodeIndexes, Found)             
             IF(.NOT. Found) CYCLE

             ElemActive(1:n) = ListGetReal( Entity, ActiveName, n, NodeIndexes, GotActive)

             ! For the 1st time we can use different limit and different active sets, if given.
             IF(FirstTime) THEN
               ElemWrk(1:n) = ListGetReal( Entity, InitName, n, NodeIndexes, Found)
               IF(Found) ElemLimit(1:n) = ElemWrk(1:n)

               ElemWrk(1:n) = ListGetReal( Entity, TRIM(ActiveName)//' Initial', n, NodeIndexes, Found)
               IF(Found) THEN
                 ElemActive(1:n) = ElemWrk(1:n)
                 GotActive = .TRUE.
               END IF

               ! We can also give the initial state in a normal manner in initial condition section.
               IF(Element % BodyId > 0 ) THEN
                 bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                     'Initial Condition', Found)
                 IF(bf>0) THEN
                   Entity => Model % ICs(bf) % Values               
                   ElemWrk(1:n) = ListGetReal( Entity, LimitName, n, NodeIndexes, Found)
                   IF(Found) ElemLimit(1:n) = ElemWrk(1:n)                   
                   ElemWrk(1:n) = ListGetReal( Entity, ActiveName, n, NodeIndexes, Found)
                   IF(Found) THEN
                     ElemActive(1:n) = ElemWrk(1:n)
                     GotActive = .TRUE.
                   END IF
                 END IF
               END IF
             END IF
               
             
             DO i=1,n
               j = FieldPerm( NodeIndexes(i) )
               IF( j == 0 ) CYCLE
               ind = Dofs * ( j - 1) + Dof

               IF( LimitDone(ind) ) CYCLE
             
               ! Go through the active set and free nodes with wrong sign in contact force
               !--------------------------------------------------------------------------       
               IF( AllActive ) THEN
                 IF(.NOT. LimitActive(ind)) THEN
                   added = added + 1
                   LimitActive(ind) = .TRUE.
                 END IF
               ELSE IF( AllPassive ) THEN
                 IF(LimitActive(ind)) THEN
                   removed = removed + 1
                   LimitActive(ind) = .FALSE.
                 END IF                 
               ELSE IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
                 IF(.NOT. LimitActive(ind)) THEN
                   added = added + 1
                   LimitActive(ind) = .TRUE.
                 END IF
               ELSE 
                 val = Var % Values(ind) 
                 IF( Upper == 0 ) THEN
                   DoAdd = ( val < ElemLimit(i) - ValEps )
                 ELSE
                   DoAdd = ( val > ElemLimit(i) + ValEps )
                 END IF
                 IF(DoAdd) THEN
                   added = added + 1
                   LimitActive(ind) = .TRUE.
                 ELSE
                   LimitActive(ind) = .FALSE.
                 END IF
               END IF
                 
               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               !--------------------------------------------------------------
               IF( LimitActive(ind) ) THEN
                 ! Set the Dirichlet conditions already here!
                 Solver % Matrix % DValues(ind) = ElemLimit(i)
                 Solver % Matrix % ConstrainedDOF(ind) = .TRUE.
                 
                 IF( Upper == 0 ) THEN
                   Var % Values(ind) = MAX( Var % Values(ind), ElemLimit(i) )
                 ELSE
                   Var % Values(ind) = MIN( Var % Values(ind), ElemLimit(i) )
                 END IF
               END IF

               LimitDone(ind) = .TRUE.             
             END DO
           END DO

           CYCLE
         END IF ! FirstTime


         IF( Conservative ) THEN
           IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
             ALLOCATE( InterfaceDof( totsize ) )
             InterfaceDof = .FALSE. 
           END IF
           
           ! Mark limited and unlimited neighbours and thereby make a 
           ! list of interface dofs. 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element
             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE             
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE
               Entity => Model % BodyForces(bf) % Values
             END IF          

             ElemLimit(1:n) = ListGetReal( Entity, &
                 LimitName, n, NodeIndexes, Found)             
             IF(.NOT. Found) CYCLE
             
             IF( DownStreamRemove ) THEN
               ! This includes only interface dofs donwstream from
               ! non-contact zone.
               BLOCK
                 REAL(kind=DP) :: r1(3),r2(3),dr(3),reps=1.0d-6
                 
                 DO i=1,n
                   j = FieldPerm( NodeIndexes(i) )
                   IF( j == 0 ) CYCLE
                   ind = Dofs * ( j - 1) + Dof
                   
                   ! Downstream of non-contact zone
                   IF(LimitActive(ind)) CYCLE
                                      
                   DO i2 = i,n
                     IF( i2 == i ) CYCLE                   
                     j2 = FieldPerm( NodeIndexes(i2) )
                     IF( j2 == 0 ) CYCLE
                     ind2 = Dofs * ( j2 - 1) + Dof
                     
                     IF( LimitActive(ind2) ) THEN
                       r2(1) =  Mesh % Nodes % x(NodeIndexes(i2))
                       r2(2) =  Mesh % Nodes % y(NodeIndexes(i2))
                       r2(3) =  Mesh % Nodes % z(NodeIndexes(i2))
                       
                       r1(1) = Mesh % Nodes % x(NodeIndexes(i))
                       r1(2) = Mesh % Nodes % y(NodeIndexes(i))
                       r1(3) = Mesh % Nodes % z(NodeIndexes(i))

                       k = DownStreamDirection 
                       IF( k > 0 ) THEN
                         dr = r2 - r1
                       ELSE
                         dr = r1 - r2
                         k = -k
                       END IF
                       
                       IF( dr(k) < reps ) CYCLE
                       
                       IF( dr(k) > 0.5*SQRT(SUM(dr*dr)) ) THEN
                         InterfaceDof(ind2) = .TRUE.
                         !PRINT *,'downstream coord:',dr
                       END IF
                     END IF
                   END DO
                 END DO
               END BLOCK
             ELSE
               ! This includes all interface dofs
               DO i=1,n
                 j = FieldPerm( NodeIndexes(i) )
                 IF( j == 0 ) CYCLE
                 ind = Dofs * ( j - 1) + Dof
                 
                 DO i2 = i+1,n
                   j2 = FieldPerm( NodeIndexes(i2) )
                   IF( j2 == 0 ) CYCLE
                   ind2 = Dofs * ( j2 - 1) + Dof
                   
                   IF( LimitActive(ind) .NEQV. LimitActive(ind2) ) THEN
                     InterfaceDof(ind) = .TRUE.
                     InterfaceDof(ind2) = .TRUE.
                   END IF
                 END DO
               END DO
             END IF
           END DO

           CALL Info(Caller,&
               'Number of interface dofs: '//I2S(COUNT(InterfaceDof)),Level=8)
         END IF ! Conservative

         IF( DownStreamRemove ) THEN
           t = COUNT(InterfaceDof)
           CALL Info(Caller,'Downstream contact set dofs:'//I2S(t),Level=8)
         END IF
         
       
         ! Add and release dofs from the contact set:
         ! If it is removed it cannot be added. 
         !----------------------------------------------------------------------
         DO t = ElemFirst, ElemLast

           Element => Model % Elements(t)
           Model % CurrentElement => Element
           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           
           Found = .FALSE.
           IF( t > Model % NumberOfBulkElements ) THEN
             DO bc = 1,Model % NumberOfBCs
               IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                 Found = .TRUE.
                 Entity => Model % BCs(bc) % Values
                 EXIT
               END IF
             END DO
             IF(.NOT. Found ) CYCLE
           ELSE             
             IF(Element % BodyId == 0) CYCLE
             bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                 'Body Force', Found)
             IF(.NOT. Found ) CYCLE
             Entity => Model % BodyForces(bf) % Values
           END IF
           
           ElemLimit(1:n) = ListGetReal( Entity, &
               LimitName, n, NodeIndexes, Found)             
           IF(.NOT. Found) CYCLE
           
           ElemActive(1:n) = ListGetReal( Entity, &
               ActiveName, n, NodeIndexes, GotActive)

           DO i=1,n
             j = FieldPerm( NodeIndexes(i) )
             IF( j == 0 ) CYCLE
             ind = Dofs * ( j - 1) + Dof

             IF( LimitDone(ind) ) CYCLE
             
             ! Go through the active set and free nodes with wrong sign in contact force
             !--------------------------------------------------------------------------       
             IF( LimitFreeze ) THEN
               CONTINUE
             ELSE IF( AllActive ) THEN
               IF(.NOT. LimitActive( ind ) ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE. 
               END IF
             ELSE IF( AllPassive ) THEN
               IF(LimitActive( ind ) ) THEN
                 removed = removed + 1
                 LimitActive(ind) = .FALSE.
               END IF               
             ELSE IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
               IF(.NOT. LimitActive( ind ) ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE. 
               END IF
             ELSE IF( LimitActive( ind ) ) THEN
               DoRemove = ( LimitSign * LoadValues(ind) > LimitSign * LoadEps ) 
               IF( DoRemove ) THEN
                 ! In the conservative mode only release nodes from contact set 
                 ! when they are adjacent to dofs that previously was not in the set.
                 ! This means that set is released only at the boundaries. 
                 IF( ConservativeRemove ) DoRemove = InterfaceDof( ind ) 
                 IF( DoRemove ) THEN
                   IF(LimitActive(ind)) THEN                     
                     removed = removed + 1
                     LimitActive(ind) = .FALSE.
                   END IF
                 END IF
               END IF
             ELSE
               ! Go through the dofs that are beyond the contact surface.
               !-----------------------------------------------------------
               val = Var % Values(ind) 
               IF( Upper == 0 ) THEN
                 DoAdd = ( val < ElemLimit(i) - ValEps )
               ELSE
                 DoAdd = ( val > ElemLimit(i) + ValEps )
               END IF
               
               IF( DoAdd ) THEN
                 IF( ConservativeAdd ) DoAdd = InterfaceDof( ind ) 
                 IF( DoAdd ) THEN
                   IF( .NOT. LimitActive(ind) ) THEN
                     added = added + 1
                     LimitActive(ind) = .TRUE.
                   END IF
                 END IF
               END IF
             END IF

             IF( LimitActive(ind) ) THEN
               ! Set the Dirichlet conditions already here!
               Solver % Matrix % DValues(ind) = ElemLimit(i)
               Solver % Matrix % ConstrainedDOF(ind) = .TRUE.

               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               IF( Upper == 0 ) THEN
                 Var % Values(ind) = MAX( Var % Values(ind), ElemLimit(i) )
               ELSE
                 Var % Values(ind) = MIN( Var % Values(ind), ElemLimit(i) )
               END IF
             END IF
             
             LimitDone(ind) = .TRUE.             
           END DO
         END DO
       END DO

       IF( DirectionActive ) THEN      
         ! Output some information before exiting
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           WRITE(Message,'(A,I0)') 'Number of lower limited dofs for '&
               //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         ELSE
           WRITE(Message,'(A,I0)') 'Number of upper limited dofs for '&
               //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         END IF
         CALL Info(Caller,Message,Level=5)
         
         IF(added + removed >= 0) THEN
           CALL Info(Caller,'Added '//I2S(added)//' and removed '&
               //I2S(removed)//' dofs in contact set',Level=6)
         END IF
       END IF
         
     END DO
                
     ! Optionally save the limiters as a field variable so that 
     ! lower limit is given value -1.0 and upper limit value +1.0.
     IF( ListGetLogical( Params,'Save Limiter',Found ) ) THEN
       
       LimitVar => VariableGet( Model % Variables, &
           GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       IF(.NOT. ASSOCIATED( LimitVar ) ) THEN
         CALL Info(Caller,'Creating field for contact: '//TRIM(GetVarName(Var)),Level=7)
         CALL VariableAddVector( Model % Variables, Solver % Mesh, Solver,&
             GetVarName(Var) //' Contact Active', Var % Dofs, Perm = FieldPerm )
         LimitVar => VariableGet( Model % Variables, &
             GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       END IF
              
       LimitVar % Values = 0.0_dp
       IF( ASSOCIATED( Var % LowerLimitActive ) ) THEN
         IF(SIZE(LimitVar % Values) /= SIZE(Var % LowerLimitActive)) THEN
           CALL Fatal(Caller,'Mismatch in size for LimitVar values!')
         END IF
         WHERE( Var % LowerLimitActive ) 
           LimitVar % Values = -1.0_dp
         END WHERE           
       END IF
       IF( ASSOCIATED( Var % UpperLimitActive ) ) THEN
         IF(SIZE(LimitVar % Values) /= SIZE(Var % UpperLimitActive)) THEN
           CALL Fatal(Caller,'Mismatch in size for LimitVar values!')
         END IF
         WHERE( Var % UpperLimitActive ) 
           LimitVar % Values = 1.0_dp
         END WHERE           
       END IF
     END IF

     IF( ALLOCATED( LimitDone ) ) THEN
       DEALLOCATE( LimitDone, ElemLimit, ElemActive, ElemWrk ) 
     END IF
     
     IF( ALLOCATED( InterfaceDof ) ) THEN
       DEALLOCATE( InterfaceDof )
     END IF

     CALL Info(Caller,'All done',Level=12)

   END SUBROUTINE DetermineSoftLimiter

END MODULE SolveCore
