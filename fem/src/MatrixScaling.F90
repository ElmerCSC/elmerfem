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
! *  Authors: Juha Ruokolainen, Peter Raback
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
!>  Linear system scaling and nodal/entity weight utilities.
!>  Extracted from SolverUtils.
!------------------------------------------------------------------------------

MODULE MatrixScaling

    USE Types,         ONLY : dp, Solver_t, Matrix_t
    USE Messages,      ONLY : Info, Warn, Fatal, Message
    USE Lists,         ONLY : ListGetLogical, ListGetString, ListGetCReal
    USE GeneralUtils,  ONLY : I2S
    USE ParallelUtils, ONLY : ParallelReduction, ParallelSumVector, &
                              ParallelInitSolve, ParallelMatrixVector, ParallelVector
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Compute the 1-norms of the rows of a CRS matrix. For a complex-valued matrix
!> (stored in the real-valued representation) the norm of the complex-valued row
!> is returned for both rows of the real-valued representation. In the parallel
!> case the absolute values of the partition contributions are summed, so that 
!> the values are the same in all partitions sharing a row, but each shared 
!> value gives only an upper bound of the norm of the truly assembled row.
!> That is, the computation of the true norms would require summing the partition
!  contributions A_{ij}^(p) before summing the absolute values row-wise.  
!------------------------------------------------------------------------------
  SUBROUTINE RowwiseOneNorms(A, RowNorm, ComplexMatrix, Parallel)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A          !< The matrix
    REAL(KIND=dp), INTENT(OUT) :: RowNorm(:) !< The 1-norms of the rows
    LOGICAL, INTENT(IN) :: ComplexMatrix     !< Is the matrix complex-valued
    LOGICAL, INTENT(IN) :: Parallel          !< Sum the partition contributions
!------------------------------------------------------------------------------
    INTEGER :: i, j, n
    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
    n = A % NumberOfRows

    IF ( ComplexMatrix ) THEN
      ! The odd row of the real-valued representation contains the pairs
      ! (Re(a_ij), -Im(a_ij)) for the complex-valued entries a_ij of the row.
      !$OMP PARALLEL DO &
      !$OMP SHARED(RowNorm, A, N) &
      !$OMP PRIVATE(i, j, s) &
      !$OMP DEFAULT(NONE)
      DO i=1,n,2
        s = 0.0_dp
        DO j=A % Rows(i), A % Rows(i+1)-1, 2
          s = s + SQRT(A % Values(j)**2 + A % Values(j+1)**2)
        END DO
        RowNorm(i) = s
        RowNorm(i+1) = s
      END DO
      !$OMP END PARALLEL DO
    ELSE
      !$OMP PARALLEL DO &
      !$OMP SHARED(RowNorm, A, N) &
      !$OMP PRIVATE(i) &
      !$OMP DEFAULT(NONE)
      DO i=1,n
        RowNorm(i) = SUM(ABS(A % Values(A % Rows(i):A % Rows(i+1)-1)))
      END DO
      !$OMP END PARALLEL DO
    END IF

    IF ( Parallel ) CALL ParallelSumVector(A, RowNorm)
!------------------------------------------------------------------------------
  END SUBROUTINE RowwiseOneNorms
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Scale linear system with different strategies.
!------------------------------------------------------------------------------
  SUBROUTINE ScaleLinearSystem(Solver,A,b,x,DiagScaling, &
          ApplyScaling,RhsScaling,ConstraintScaling,ScalingStr,Shift)
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:), x(:)
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)
    LOGICAL, OPTIONAL :: ApplyScaling, RhsScaling,ConstraintScaling
    CHARACTER(*), OPTIONAL :: ScalingStr
    COMPLEX(KIND=dp), OPTIONAL :: Shift  !< The shift sigma when A % Values = K - sigma*M with M = A % MassValues
    INTEGER :: n

    CHARACTER(:), ALLOCATABLE :: str
    LOGICAL :: Parallel, Found, ComplexMatrix

    n = A % NumberOfRows
    Parallel = Solver % Parallel

    IF( ListGetLogical( Solver % Values,'Linear System Pseudo Complex',Found ) ) THEN
      ComplexMatrix = .TRUE.
    ELSE
      ComplexMatrix = A % COMPLEX
    END IF

    IF( PRESENT( ScalingStr ) ) THEN
      str = ScalingStr
    ELSE IF( ListGetLogical(Solver % Values,'Linear System Row Equilibration',Found ) ) THEN
      str = 'row equilibration'
    ELSE
      str = ListGetString( Solver % Values,'Linear System Scaling Method',Found )
      IF(.NOT. Found) str = 'diagonal'
    END IF

    SELECT CASE( str )
    CASE('diagonal')
      CALL ScaleLinearSystemDiagonal()

    CASE('row equilibration','rowsum')
      CALL RowEquilibration(Solver, A, b, Parallel, ApplyScaling )

    CASE('constant')
      CALL ScaleLinearSystemConstant()

    CASE('none')
      CALL Info('ScaleLinearSystem','No scaling will be applied!',Level=12)
      RETURN

    CASE DEFAULT
      CALL Fatal('ScaleLinearSystem','Unknown scaling method: '//TRIM(str))
    END SELECT


  CONTAINS

    !-------------------------------------------------------------
    !> Scale system Ax = b as (DAD)y = Db, where
    !> D = 1/SQRT(Diag(A)) and y = D^-1 x. By default, if the scaled
    !> system is rewritten as A'y = b', this subroutine also performs
    !> an additional scaling so that the final form is given by
    !> A'(y/|b'|) = b'/|b'|. Whether the last step is taken depends
    !> on the optional argument RhsScaling.
    !-------------------------------------------------------------
    SUBROUTINE ScaleLinearSystemDiagonal()

      INTEGER :: i,j,nRepl
      REAL(KIND=dp) :: bnorm,s,r,m,DiagTol
      COMPLEX(KIND=dp) :: DiagC, MassC, ShiftC
      LOGICAL :: DoRHS, DoCM, Found, Shifted
      REAL(KIND=dp), POINTER  :: Diag(:)
      REAL(KIND=dp), ALLOCATABLE :: RowSum(:), MDiag(:)
      TYPE(Matrix_t), POINTER :: CM


      A % ScalingMethod = 1

      IF( PRESENT( DiagScaling ) ) THEN
        CALL Info('ScaleLinearSystem','Reusing existing > DiagScaling < vector',Level=12)
        Diag => DiagScaling
      ELSE
        IF(.NOT. ASSOCIATED(A % DiagScaling)) THEN
        CALL Info('ScaleLinearSystem','Creating > DiagScaling < vector of size '//I2S(n),Level=10)
          ALLOCATE( A % DiagScaling(n) )
        ELSE
          CALL Info('ScaleLinearSystem','Recomputing > DiagScaling < vector of size '//I2S(n),Level=12)
        END IF
        Diag => A % DiagScaling
        Diag(1:n) = 0._dp

        ! If the matrix has been obtained by shifting as K - sigma*M, a diagonal entry
        ! k_ii - sigma*m_ii may become negligible owing to cancellation, which would give
        ! a huge scaling factor. Such entries are detected by comparing with the upper
        ! bound |k_ii| + |sigma*m_ii| of the diagonal entry, which is then used in the scaling
        ! instead of the diagonal entry. This test does not depend on how the units of 
        ! the field variables have been chosen.
        !
        Shifted = PRESENT(Shift) .AND. ASSOCIATED(A % MassValues)
        IF ( Shifted ) Shifted = SIZE(A % MassValues) == SIZE(A % Values) .AND. ABS(Shift) > 0.0_dp
        ShiftC = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        DiagTol = 0.0_dp
        IF ( Shifted ) THEN
          ShiftC = Shift
          DiagTol = ListGetCReal(Solver % Values, 'Linear System Scaling Diagonal Tolerance', &
              Found, minv = 0.0_dp, DefValue = 1.0e-10_dp)
          ALLOCATE(MDiag(n))
          MDiag(1:n) = 0._dp
        END IF

        IF ( ComplexMatrix ) THEN
          CALL Info('ScaleLinearSystem','Assuming complex matrix while scaling',Level=20)

          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, MDiag, Shifted, A, N) &
          !$OMP PRIVATE(i, j) &
          !$OMP DEFAULT(NONE)
          DO i=1,n,2
            j = A % Diag(i)
            IF(j>0) THEN
              Diag(i)   = A % Values(j)
              Diag(i+1) = A % Values(j+1)
              IF ( Shifted ) THEN
                MDiag(i)   = A % MassValues(j)
                MDiag(i+1) = A % MassValues(j+1)
              END IF
            ELSE
              Diag(i) = 0._dp
              Diag(i+1) = 0._dp
            END IF
          END DO
          !$OMP END PARALLEL DO
        ELSE
          CALL Info('ScaleLinearSystem','Assuming real valued matrix while scaling',Level=25)

          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, MDiag, Shifted, A, N) &
          !$OMP PRIVATE(i, j) &
          !$OMP DEFAULT(NONE)
          DO i=1,n
            j = A % Diag(i)
            IF (j>0) THEN
              Diag(i) = A % Values(j)
              IF ( Shifted ) MDiag(i) = A % MassValues(j)
            END IF
          END DO
          !$OMP END PARALLEL DO
        END IF

        IF ( Parallel ) THEN
          CALL Info('ScaleLinearSystem','Performing parallel summation of > DiagScaling < vector',Level=20)
          CALL ParallelSumVector(A, Diag)
          IF ( Shifted ) CALL ParallelSumVector(A, MDiag)
        END IF

        ! Replace the diagonal entries by their magnitudes, with the cancelled entries
        ! replaced as described above:
        nRepl = 0
        IF ( ComplexMatrix ) THEN
          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, MDiag, Shifted, ShiftC, DiagTol, N) &
          !$OMP PRIVATE(i, s, r, DiagC, MassC) &
          !$OMP REDUCTION(+:nRepl) &
          !$OMP DEFAULT(NONE)
          DO i=1,n,2
            DiagC = CMPLX(Diag(i), -Diag(i+1), KIND=dp)
            s = ABS(DiagC)
            IF ( Shifted .AND. s > TINY(s) ) THEN
              MassC = ShiftC * CMPLX(MDiag(i), -MDiag(i+1), KIND=dp)
              r = ABS(DiagC + MassC) + ABS(MassC)
              IF ( s <= DiagTol * r ) THEN
                s = r
                nRepl = nRepl + 1
              END IF

            END IF
            Diag(i) = s
            Diag(i+1) = s
          END DO
          !$OMP END PARALLEL DO
        ELSE
          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, MDiag, Shifted, ShiftC, DiagTol, N) &
          !$OMP PRIVATE(i, s, r, m) &
          !$OMP REDUCTION(+:nRepl) &
          !$OMP DEFAULT(NONE)
          DO i=1,n
            s = ABS(Diag(i))
            IF ( Shifted .AND. s > TINY(s) ) THEN
              m = REAL(ShiftC) * MDiag(i)
              r = ABS(Diag(i) + m) + ABS(m)
              IF ( s <= DiagTol * r ) THEN
                s = r
                nRepl = nRepl + 1
              END IF
            END IF
            Diag(i) = s
          END DO
          !$OMP END PARALLEL DO
        END IF

        IF ( Shifted ) THEN
          DEALLOCATE(MDiag)
          IF ( nRepl > 0 ) CALL Info('ScaleLinearSystem', &
              'Diagonal entries replaced owing to cancellation: '//I2S(nRepl),Level=10)
        END IF

        ! A zero diagonal entry (as in saddle point problems) is replaced by the 1-norm
        ! of the row:
        s = 0.0_dp
        IF ( ANY(Diag(1:n) <= TINY(s)) ) s = 1.0_dp
        IF ( Parallel ) s = ParallelReduction(s,2)
        IF ( s > 0.0_dp ) THEN
          ALLOCATE(RowSum(n))
          CALL RowwiseOneNorms(A, RowSum, ComplexMatrix, Parallel)
          WHERE( Diag(1:n) <= TINY(s) ) Diag(1:n) = RowSum(1:n)
          DEALLOCATE(RowSum)
        END IF

        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, N) &
        !$OMP PRIVATE(i) &
        !$OMP DEFAULT(NONE)
        DO i=1,n
          IF ( Diag(i) > TINY(Diag(i)) ) THEN
            Diag(i) = 1.0_dp / SQRT(Diag(i))
          ELSE
            Diag(i) = 1.0_dp
          END IF
        END DO
        !$OMP END PARALLEL DO
      END IF


      ! Optionally we may just create the diag and leave the scaling undone
      !--------------------------------------------------------------------
      IF( PRESENT( ApplyScaling ) ) THEN
        IF(.NOT. ApplyScaling ) THEN
          CALL Info('ScaleLinearSystem','Application of scaling skipped!',Level=20)
          RETURN
        END IF
      END IF

      CALL Info('ScaleLinearSystem','Scaling matrix values',Level=20)

      !$OMP PARALLEL &
      !$OMP SHARED(Diag, A, N) &
      !$OMP PRIVATE(i,j) &
      !$OMP DEFAULT(NONE)

      !$OMP DO
      DO i=1,n
        DO j = A % Rows(i), A % Rows(i+1)-1
          A % Values(j) = A % Values(j) * &
              ( Diag(i) * Diag(A % Cols(j)) )
        END DO
      END DO
      !$OMP END DO NOWAIT

      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
          CALL Info('ScaleLinearSystem','Scaling PrecValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % PrecValues(j) = A % PrecValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          CALL Info('ScaleLinearSystem','Scaling MassValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % MassValues(j) = A % MassValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          CALL Info('ScaleLinearSystem','Scaling DampValues',Level=20)
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % DampValues(j) = A % DampValues(j) * &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      !$OMP END PARALLEL

      DoCM = .FALSE.
      IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling

      IF(doCM) THEN
        CM => A % ConstraintMatrix
        IF (ASSOCIATED(CM)) THEN
          CALL Info('ScaleLinearSystem','Scaling Constraints',Level=20)
          !$OMP PARALLEL DO &
          !$OMP SHARED(Diag, CM) &
          !$OMP PRIVATE(i,j) &
          !$OMP DEFAULT(NONE)
          DO i=1,CM % NumberOFRows
            DO j=CM % Rows(i), CM % Rows(i+1)-1
              CM % Values(j) = CM % Values(j) * Diag(CM % Cols(j))
            END DO
          END DO
          !$OMP END PARALLEL DO
        END IF
      END IF

      ! Scale r.h.s. and initial guess
      !--------------------------------
      A % RhsScaling=1._dp
      ! TODO: Add threading
      IF( PRESENT( b ) ) THEN
        CALL Info('ScaleLinearSystem','Scaling Rhs vector',Level=20)

        b(1:n) = b(1:n) * Diag(1:n)
        DoRHS = .TRUE.
        IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
        IF (DoRHS) THEN
          IF( Parallel ) THEN
            BLOCK
              REAL(KIND=dp), ALLOCATABLE :: s(:), r(:)
              ALLOCATE(s(n), r(n))
              s = 0
              r = b(1:n)
              CALL ParallelSUMVector(A, r)
              CALL ParallelVector(A, s, r)
              bnorm = SQRT(ParallelReduction(SUM(s**2)))
            END BLOCK
          ELSE
            bnorm = SQRT(SUM(b(1:n)**2))
          END IF

          IF( bnorm < SQRT( TINY( bnorm ) ) ) THEN
            CALL Info('ScaleLinearSystem','Rhs vector is almost zero, skipping rhs scaling!',Level=20)
          ELSE
            A % RhsScaling = bnorm
            b(1:n) = b(1:n) / bnorm
          END IF
        END IF
      END IF

      IF( PRESENT(x) ) THEN
        x(1:n) = x(1:n) / (Diag(1:n) * A % RhsScaling)
      END IF

    END SUBROUTINE ScaleLinearSystemDiagonal


    !-------------------------------------------------------------
    !>  Scale system Ax = b as (A/Ascl)y = b, with y = Ascl x.
    !>  By default, if the scaled system is rewritten as A'y = b',
    !>  this subroutine also perform an additional scaling, so that
    !>  the final form is given by A'(y/bscl) = b'/bscl, i.e.
    !>  (A/Ascl)(Ascl*x/bscl) = (b/bscl). Whether the last step is
    !>  taken depends on the optional argument RhsScaling.
    !-------------------------------------------------------------
    SUBROUTINE ScaleLinearSystemConstant()

      INTEGER :: i,j,nSum
      REAL(KIND=dp) :: bnorm,s
      COMPLEX(KIND=dp) :: DiagC
      LOGICAL :: DoRHS, DoCM
      REAL(KIND=dp) :: Ascl, bscl, Xscl, bsum, DiagSum
      TYPE(Matrix_t), POINTER :: CM

      IF( Parallel ) THEN
        CALL Info('ScaleLinearSystem','Scaling matrix entries by constant in parallel',Level=10)
      ELSE
        CALL Info('ScaleLinearSystem','Scaling matrix entries by constant in serial',Level=10)
      END IF

      CALL Info('ScaleLinearSystem','Computing > AveScaling < constant',Level=12)

      A % ScalingMethod = 3

      DiagSum = 0.0_dp
      nSum = n

      DO i=1,n
        j = A % Diag(i)
        IF (j>0) THEN
          DiagSum = DiagSum + ABS(A % Values(j))
        END IF
      END DO

      IF ( Parallel ) THEN
        DiagSum = ParallelReduction( DiagSum )
        nSum = ParallelReduction( nSum )
      END IF

      Ascl = DiagSum / nSum
      A % AveScaling = Ascl

      WRITE( Message,'(A,ES12.3)') 'Average diagonal entry: ', Ascl
      CALL Info( 'ScaleLinearSystemConstant', Message, Level=8 )

      ! Optionally we may just create the diag and leave the scaling undone
      !--------------------------------------------------------------------
      IF( PRESENT( ApplyScaling ) ) THEN
        IF(.NOT. ApplyScaling ) RETURN
      END IF

      A % Values = A % Values / Ascl
      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
          A % PrecValues = A % PrecValues / Ascl
        END IF
      END IF
      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          A % MassValues = A % MassValues / Ascl
        END IF
      END IF
      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          A % DampValues = A % DampValues / Ascl
        END IF
      END IF

      DoCM = .FALSE.
      IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling
      IF(doCM) THEN
        CM => A % ConstraintMatrix
        IF (ASSOCIATED(CM)) THEN
          CM % Values = CM % Values / Ascl
        END IF
      END IF

      ! Scale r.h.s. and initial guess
      !--------------------------------
      A % RhsScaling = 1._dp

      IF( PRESENT( b ) ) THEN
        DoRHS = .TRUE.
        IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
        IF (DoRHS) THEN
          bsum = SUM( ABS( b(1:n) ) )
          nSum = n

          IF ( Parallel ) THEN
            bSum = ParallelReduction( bSum )
            nSum = ParallelReduction( nSum )
          END IF

          bscl = bsum / nSum
          b = b / bscl

          A % RhsScaling = bscl

          WRITE( Message,'(A,ES12.3)') 'Average rhs entry: ', bscl
          CALL Info( 'ScaleLinearSystemConstant', Message, Level=7 )
        END IF
      END IF

      IF( PRESENT(x) ) THEN
        Xscl = A % RhsScaling / Ascl
        x(1:n) = x(1:n) / Xscl
      END IF

    END SUBROUTINE ScaleLinearSystemConstant

!-----------------------------------------------------------------------------
  END SUBROUTINE ScaleLinearSystem
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!>   Equilibrate the rows of the coefficient matrix A to
!>   minimize the condition number. The associated rhs vector f is also scaled.
!------------------------------------------------------------------------------
  SUBROUTINE RowEquilibration( Solver, A, f, Parallel, ApplyScaling )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: A
    REAL(KIND=dp), OPTIONAL :: f(:)
    LOGICAL :: Parallel
    LOGICAL, OPTIONAL :: ApplyScaling
!-----------------------------------------------------------------------------
    LOGICAL :: ComplexMatrix, Found
    INTEGER :: i, j, n
    REAL(kind=dp) :: norm
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-------------------------------------------------------------------------

    CALL Info('RowEquilibration','Scaling system such that abs rowsum is unity',Level=15)

    n = A % NumberOfRows
    ComplexMatrix = A % COMPLEX

    IF (ComplexMatrix) CALL Info('RowEquilibration','Using complex matrix norm',Level=15)

    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    IF( .NOT. ASSOCIATED(A % DiagScaling) ) THEN
      ALLOCATE( A % DiagScaling(n) )
    END IF
    Diag => A % DiagScaling

    Diag = 0.0d0
    norm = 0.0d0

    A % ScalingMethod = 2

    !---------------------------------------------
    ! Compute 1-norm of each row
    !---------------------------------------------
    IF (Parallel .AND. ListGetLogical(Solver % Values, 'Row Equilibration Use M-V', Found)) THEN
      IF (ComplexMatrix) THEN
        CALL Warn('RowEquilibration', '"Row Equilibration Use M-V" is not available for complex matrices')
      END IF
    END IF
    IF (Parallel .AND. .NOT. ComplexMatrix .AND. &
        ListGetLogical(Solver % Values, 'Row Equilibration Use M-V', Found, DefValue = .TRUE.)) THEN
      ! The matrix-vector product |A|*1 is evaluated using the splitted parallel matrix
      ! where the partition contributions have been glued, so that each entry of the
      ! global matrix is stored once and has its final value. The result is thus the
      ! row 1-norm of the global matrix (the same as obtained in serial execution, up
      ! to the order of summation). It is however available only for the rows owned
      ! by the partition; the other rows are zeroed and then obtained from the owners.
      BLOCK
        REAL(KIND=dp), ALLOCATABLE :: x(:),y(:),z(:)
        TYPE(Matrix_t), POINTER :: ap
        ALLOCATE(x(n),y(n),z(n))
        ap => a
        x = 1; y=0; z=0
        CALL ParallelInitSolve( ap, x, y, z )
        CALL ParallelMatrixVector( ap, x, Diag, Update=.TRUE., ZeroNotOwned=.TRUE., UseABS=.TRUE. )
      END BLOCK
      CALL ParallelSUMVector(A,Diag)
    ELSE
      CALL RowwiseOneNorms(A, Diag, ComplexMatrix, Parallel)
    END IF

    norm = MAXVAL(Diag(1:n))
    IF( Parallel ) THEN
      norm = ParallelReduction(norm,2)
    END IF

    WRITE( Message, * ) 'Unscaled matrix norm: ', norm
    CALL Info( 'RowEquilibration', Message, Level=5 )

    !--------------------------------------------------
    ! Now, define the scaling matrix by inversion and
    ! perform the actual scaling of the linear system
    !--------------------------------------------------
    IF (ComplexMatrix) THEN
      DO i=1,n,2
        IF (Diag(i) > TINY(norm) ) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
        Diag(i+1) = Diag(i)
      END DO
    ELSE
      DO i=1,n
        IF (Diag(i) > TINY(norm)) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
      END DO
    END IF

    IF( PRESENT( ApplyScaling ) ) THEN
      IF(.NOT. ApplyScaling) THEN
        CALL Info('RowEquilibration', 'Application of scaling skipped!', Level=20)
        RETURN
      END IF
    END IF

    DO i=1,n
      DO j=Rows(i),Rows(i+1)-1
        Values(j) = Values(j) * Diag(i)
      END DO
    END DO

    IF (PRESENT(f)) THEN
      DO i=1,n
        f(i) = Diag(i) * f(i)
      END DO
    END IF

    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) * Diag(i)
          END DO
        END DO
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE RowEquilibration
!------------------------------------------------------------------------------



!--------------------------------------------------------------
!>  Scale the system back to original.
!--------------------------------------------------------------
  SUBROUTINE BackScaleLinearSystem( Solver,A,b,x,DiagScaling,&
      ConstraintScaling, EigenScaling )

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    LOGICAL, OPTIONAL :: ConstraintScaling, EigenScaling
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)

    CHARACTER(:), ALLOCATABLE :: str
    LOGICAL :: Found
    INTEGER :: n

    CALL Info('BackScaleLinearSystem','Scaling back to original scale',Level=14)

    n = A % NumberOfRows

    SELECT CASE( A % ScalingMethod )

    CASE( 1 )
      CALL BackScaleLinearSystemDiagonal()

    CASE( 2 )
      CALL ReverseRowEquilibration( A, b )

    CASE( 3 )
      CALL BackScaleLinearSystemConstant()

    CASE DEFAULT
      CALL Info('BackScaleLinearSystem','No scaling was applied!',Level=20)

    END SELECT

  CONTAINS

    SUBROUTINE BackScaleLinearSystemDiagonal()

      REAL(KIND=dp), POINTER :: Diag(:)
      INTEGER :: i,j
      LOGICAL :: doCM
      TYPE(Matrix_t), POINTER :: CM

      IF( PRESENT( DiagScaling ) ) THEN
        Diag => DiagScaling
      ELSE
        Diag => A % DiagScaling
      END IF

      IF(.NOT. ASSOCIATED( Diag ) ) THEN
        CALL Warn('BackScaleLinearSystem','Diag not associated!')
        RETURN
      END IF

      IF( SIZE( Diag ) /= n ) THEN
        CALL Fatal('BackScaleLinearSystem','Diag of wrong size!')
      END IF

      ! TODO: Add threading
      !
      !      Solve x:  INV(D)x = y
      !      -------------------------------------------
      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Diag(1:n) * A % RhsScaling
      END IF

      IF( PRESENT( b ) ) THEN
        b(1:n) = b(1:n) / Diag(1:n) * A % RhsScaling
      END IF

      IF( PRESENT( EigenScaling ) ) THEN
        IF( EigenScaling ) THEN
          ! TODO: Add threading
          DO i=1,Solver % NOFEigenValues
            !
            !           Solve x:  INV(D)x = y
            !           --------------------------
            IF ( Solver % Matrix % COMPLEX ) THEN
              Solver % Variable % EigenVectors(i,1:n/2) = &
                  Solver % Variable % EigenVectors(i,1:n/2) * Diag(1:n:2)
            ELSE
              Solver % Variable % EigenVectors(i,1:n) = &
                  Solver % Variable % EigenVectors(i,1:n) * Diag(1:n)
            END IF
          END DO
        END IF
      END IF

      !$OMP PARALLEL &
      !$OMP SHARED(Diag, A, N) &
      !$OMP PRIVATE(i, j) &
      !$OMP DEFAULT(NONE)

      !$OMP DO
      DO i=1,n
        DO j=A % Rows(i), A % Rows(i+1)-1
          A % Values(j) = A % Values(j) / (Diag(i) * Diag(A % Cols(j)))
        END DO
      END DO
      !$OMP END DO NOWAIT

      IF ( ASSOCIATED( A % PrecValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % PrecValues(j) = A % PrecValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % MassValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % MassValues(j) = A % MassValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      IF ( ASSOCIATED( A % DampValues ) ) THEN
        IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
          !$OMP DO
          DO i=1,n
            DO j=A % Rows(i), A % Rows(i+1)-1
              A % DampValues(j) = A % DampValues(j) / &
                  ( Diag(i) * Diag(A % Cols(j)) )
            END DO
          END DO
          !$OMP END DO NOWAIT
        END IF
      END IF

      !$OMP END PARALLEL

      ! TODO: Add threading
      doCM=.FALSE.
      IF(PRESENT(ConstraintScaling)) doCM=ConstraintScaling
      IF(doCM) THEN
        CM => A % ConstraintMatrix
        IF (ASSOCIATED(CM)) THEN
          DO i=1,CM % NumberOFRows
            DO j=CM % Rows(i), CM % Rows(i+1)-1
              CM % Values(j) = CM % Values(j) / ( Diag(CM % Cols(j)) )
            END DO
          END DO
        END IF
      END IF

      A % RhsScaling=1._dp

      IF(.NOT. PRESENT(DiagScaling) ) THEN
        DEALLOCATE(A % DiagScaling)
        A % DiagScaling=>NULL()
      END IF

    END SUBROUTINE BackScaleLinearSystemDiagonal


    SUBROUTINE BackScaleLinearSystemConstant()

      REAL(KIND=dp) :: Ascl, Bscl, Xscl
      LOGICAL :: doCM
      TYPE(Matrix_t), POINTER :: CM

      Ascl = A % AveScaling
      Bscl = A % RhsScaling
      Xscl = Bscl / Ascl

      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Xscl
      END IF

      IF( PRESENT( b ) ) THEN
        b(1:n) = Bscl * b(1:n)
      END IF

      IF( PRESENT( EigenScaling ) ) THEN
        IF( EigenScaling ) THEN
          Solver % Variable % EigenVectors = &
              Solver % Variable % EigenVectors * Xscl
        END IF
      END IF

      IF( Ascl > 0.0 .AND. ABS(Ascl-1.0_dp) > EPSILON(Ascl) ) THEN
        A % Values = Ascl * A % Values
        IF ( ASSOCIATED( A % PrecValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
            A % PrecValues = AScl * A % PrecValues
          END IF
        END IF
        IF ( ASSOCIATED( A % MassValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
            A % MassValues = Ascl * A % MassValues
          END IF
        END IF
        IF ( ASSOCIATED( A % DampValues ) ) THEN
          IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
            A % DampValues = Ascl * A % DampValues
          END IF
        END IF

        doCM = .FALSE.
        IF(PRESENT(ConstraintScaling)) doCM=ConstraintScaling
        IF(doCM) THEN
          CM => A % ConstraintMatrix
          IF( ASSOCIATED(CM ) ) THEN
            CM % Values = Ascl * CM % Values
          END IF
        END IF
      END IF

      A % RhsScaling = 1._dp

    END SUBROUTINE BackScaleLinearSystemConstant


  END SUBROUTINE BackScaleLinearSystem


!------------------------------------------------------------------------------
!> Scale the linear system back to original when the linear
!> system scaling has been done by row equilibration.
!------------------------------------------------------------------------------
  SUBROUTINE ReverseRowEquilibration( A, f )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: f(:)
!-----------------------------------------------------------------------------
    INTEGER :: i, j, n
    INTEGER, POINTER :: Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-----------------------------------------------------------------------------
    CALL Info('ReverseRowEquilibration','Scaling back to original scale',Level=14)

    n = A % NumberOfRows
    Diag => A % DiagScaling
    Values => A % Values
    Rows => A % Rows

    IF(.NOT. ASSOCIATED( Diag ) ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag not associated!')
    END IF
    IF( SIZE( Diag ) /= n ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag of wrong size!')
    END IF

    IF (PRESENT(f)) f(1:n) = f(1:n) / Diag(1:n)
    DO i=1,n
      DO j = Rows(i), Rows(i+1)-1
        Values(j) = Values(j) / Diag(i)
      END DO
    END DO

    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) / Diag(i)
          END DO
        END DO
      END IF
    END IF

    DEALLOCATE(A % DiagScaling)
    A % DiagScaling => NULL()

!------------------------------------------------------------------------------
  END SUBROUTINE ReverseRowEquilibration
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> A simplified subroutine for scaling vectors so that they are transformed
!> in the same way as the RHS and the solution vector of the linear system
!> associated with the given matrix A.
!------------------------------------------------------------------------------
  SUBROUTINE ScaleLinearSystemVectors(A, b, n, x, BackScaling)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A
    REAL(KIND=dp), INTENT(INOUT) :: b(:)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: x(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: BackScaling
!------------------------------------------------------------------------------
    LOGICAL :: Backwards
!------------------------------------------------------------------------------
    IF (PRESENT(BackScaling)) THEN
      Backwards = BackScaling
    ELSE
      Backwards = .FALSE.
    END IF

    IF (Backwards) THEN
      !
      ! Perform back-scaling
      !
      SELECT CASE(A % ScalingMethod)
      CASE(1)
        b(1:n) = b(1:n) / A % DiagScaling(1:n) * A % RhsScaling
        IF (PRESENT(x)) x(1:n) = x(1:n) * A % DiagScaling(1:n) * A % RhsScaling
      CASE(2)
        b(1:n) = b(1:n) / A % DiagScaling(1:n)
      CASE(3)
        b(1:n) = A % RhsScaling * b(1:n)
        IF (PRESENT(x)) x(1:n) = (A % RhsScaling / A % AveScaling) * x(1:n)
      CASE DEFAULT
        CALL Fatal('ScaleLinearSystemVectors', 'Unknown method for back-scaling')
      END SELECT
    ELSE
      !
      ! Transform the vectors so that they correspond to the scaled version of linear system
      !
      SELECT CASE(A % ScalingMethod)
      CASE(1)
        b(1:n) = A % DiagScaling(1:n) * b(1:n) / A % RhsScaling
        IF (PRESENT(x)) x(1:n) = x(1:n) / (A % DiagScaling(1:n) * A % RhsScaling)
      CASE(2)
        b(1:n) = A % DiagScaling(1:n) * b(1:n)
      CASE(3)
        b(1:n) = b(1:n) / A % RhsScaling
        IF (PRESENT(x)) x(1:n) = A % AveScaling * x(1:n) / A % RhsScaling
      CASE DEFAULT
        CALL Fatal('ScaleLinearSystemVectors', 'Unknown method for scaling')
      END SELECT
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ScaleLinearSystemVectors
!------------------------------------------------------------------------------

END MODULE MatrixScaling
