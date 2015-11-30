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
! *  Original Date: 01 Oct 1996
! *
! ****************************************************************************/

#include "huti_fdefs.h"

!> \ingroup ElmerLib 
!> \{

!-----------------------------------------------------------------------------
!>  Module defining utility routines & matrix storage for sparse
!>  matrix in Compressed Row Storage (CRS) format.
!-----------------------------------------------------------------------------
MODULE CRSMatrix

  USE Lists
#ifdef USE_ISO_C_BINDINGS
  USE LoadMod
#endif

  IMPLICIT NONE

CONTAINS


!-----------------------------------------------------------------------------
!> Search CRS matrix for what?
!-----------------------------------------------------------------------------
  FUNCTION CRS_Search( N,Array,VALUE ) RESULT ( Index )
!-----------------------------------------------------------------------------
    INTEGER :: N,VALUE,Array(:)
!-----------------------------------------------------------------------------
    INTEGER :: Lower, Upper,Lou,Index
!-----------------------------------------------------------------------------
    Index = 0 
    Upper = N
    Lower = 1

    ! Handle the special case

    IF ( Upper == 0 ) RETURN

    DO WHILE( .TRUE. )
      IF ( Array(Lower) == VALUE ) THEN
        Index = Lower
        EXIT
      ELSE IF ( Array(Upper) == VALUE ) THEN
        Index = Upper
        EXIT
      END IF

      IF ( (Upper-Lower)>1 ) THEN
        Lou = ISHFT((Upper+Lower), -1)
        IF ( Array(Lou) < VALUE ) THEN
          Lower = Lou
        ELSE
          Upper = Lou
        END IF
      ELSE
        EXIT
      END IF
    END DO
    
    RETURN

  END FUNCTION CRS_Search
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Zero a CRS format matrix
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ZeroMatrix(A)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding the matrix
!------------------------------------------------------------------------------
    A % Values = 0.0d0
  END SUBROUTINE CRS_ZeroMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Zero given row from a CRS format matrix
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ZeroRow( A,n )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A !< Structure holding the matrix
    INTEGER :: n        !< Row number to be zeroed
!------------------------------------------------------------------------------
 
    INTEGER :: i

    LOGICAL :: isMass, isDamp, EigenAnalysis, HarmonicAnalysis, Found

    EigenAnalysis=.FALSE.; HarmonicAnalysis=.FALSE.
    IF(ASSOCIATED(A % Solver)) THEN
       EigenAnalysis = A % Solver % NOFEigenValues > 0 .AND. &
           ListGetLogical( A % Solver % Values, 'Eigen Analysis',Found)

       HarmonicAnalysis = A % Solver % NOFEigenValues>0 .AND. &
          ListGetLogical( A % Solver % Values, 'Harmonic Analysis',Found)
    END IF

    isMass = (EigenAnalysis.OR.HarmonicAnalysis).AND.ASSOCIATED(A % MassValues)
    IF ( isMass ) &
      isMass = isMass .AND. SIZE(A % MassValues) == SIZE(A % Values)

    isDamp = (EigenAnalysis.OR.HarmonicAnalysis).AND.ASSOCIATED(A % DampValues)
    IF ( isDamp ) &
      isDamp = isDamp .AND. SIZE(A % DampValues) == SIZE(A % Values)

    DO i=A % Rows(n),A % Rows(n+1)-1
      A % Values(i) = 0.0_dp
    END DO

    IF ( isMass ) THEN
      DO i=A % Rows(n),A % Rows(n+1)-1
         A % MassValues(i) = 0.0_dp
      END DO
    END IF

    IF ( isDamp )  THEN
      DO i=A % Rows(n),A % Rows(n+1)-1
        A % DampValues(i) = 0.0_dp
      END DO
    END IF

  END SUBROUTINE CRS_ZeroRow
!------------------------------------------------------------------------------
  


!------------------------------------------------------------------------------
!>    Sort columns to ascending order for rows of a CRS format matrix-
!>    Optionally also sort the corresponging values of the matrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_SortMatrix( A, ValuesToo )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A   !< Structure holding the matrix
    LOGICAL, OPTIONAL, INTENT(IN) :: ValuesToo  !< Should the values be sorted as well?
!------------------------------------------------------------------------------

    INTEGER :: i,j,n
    LOGICAL :: SortValues
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)

    SortValues = .FALSE.
    IF ( PRESENT(ValuesToo) ) SortValues=ValuesToo

    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    IF ( SortValues ) Values => A % Values
    n = A % NumberOfRows

    IF ( .NOT. A % Ordered ) THEN
      IF ( SortValues ) THEN
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
        !$OMP SHARED(Rows, Cols, Values, N) &
        !$OMP PRIVATE(i)
        DO i=1,N
          CALL SortF( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),Values(Rows(i):Rows(i+1)-1) )
        END DO
        !$OMP END PARALLEL DO
      ELSE
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
        !$OMP SHARED(Rows, Cols, N) &
        !$OMP PRIVATE(i)
        DO i=1,N
          CALL Sort( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1) )
        END DO
        !$OMP END PARALLEL DO
      END IF

      IF ( ASSOCIATED(Diag) ) THEN
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
        !$OMP SHARED(Diag, Rows, Cols, N) &
        !$OMP PRIVATE(i,j)
        DO i=1,N
          DO j=Rows(i),Rows(i+1)-1
            IF ( Cols(j) == i ) THEN
              Diag(i) = j
              EXIT
            END IF
          END DO
        END DO
        !$OMP END PARALLEL DO
      END IF

      A % Ordered = .TRUE.
    END IF

  END SUBROUTINE CRS_SortMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Sort columns to ascending order for rows of a CRS format matrix-
!>    Optionally also sort the corresponging values of the matrix.
!>    This operates just on a basic matrix type. 
!------------------------------------------------------------------------------
  SUBROUTINE CRS_SortBasicMatrix( A, ValuesToo )
!------------------------------------------------------------------------------
    TYPE(BasicMatrix_t), TARGET :: A  !< Structure holding the matrix
    LOGICAL, OPTIONAL, INTENT(IN) :: ValuesToo    !< Should the values be sorted as well?
!------------------------------------------------------------------------------

    INTEGER :: i,j,n

    LOGICAL :: SortValues
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)

    SortValues = .FALSE.
    IF ( PRESENT(ValuesToo) ) SortValues=ValuesToo

    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    IF ( SortValues ) Values => A % Values
    n = A % NumberOfRows

    IF ( SortValues ) THEN
      DO i=1,N
        CALL SortF( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),Values(Rows(i):Rows(i+1)-1) )
      END DO
    ELSE
      DO i=1,N
          CALL Sort( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1) )
      END DO
    END IF

    IF ( ASSOCIATED(Diag) ) THEN
      DO i=1,N
        DO j=Rows(i),Rows(i+1)-1
          IF ( Cols(j) == i ) THEN
            Diag(i) = j
            EXIT
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE CRS_SortBasicMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Fill in the column number to a CRS format matrix (values are not 
!>    affected in any way).
!------------------------------------------------------------------------------
  SUBROUTINE CRS_MakeMatrixIndex( A,i,j )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding matrix
    INTEGER, INTENT(IN) :: i         !< row number of the matrix element
    INTEGER, INTENT(IN) :: j         !< column number of the matrix element
!------------------------------------------------------------------------------

    INTEGER :: k,n
    INTEGER, POINTER :: Cols(:),Rows(:)

    Rows   => A % Rows
    Cols   => A % Cols

    n = Rows(i)
    DO k=Rows(i),Rows(i+1)-1
      IF ( Cols(k) == j ) THEN
        RETURN
      ELSE IF ( Cols(k) < 1 ) THEN
        n = k
        EXIT
      END IF
    END DO

    IF ( Cols(n) >= 1 ) THEN
      WRITE( Message, * ) 'Trying to access non-existent column:',n,Cols(n)
      CALL Error( 'MakeMatrixIndex', Message )
      RETURN
    END IF

    Cols(n) = j
  END SUBROUTINE CRS_MakeMatrixIndex
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Add a given value to an element of a  CRS format matrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_AddToMatrixElement( A,i,j,VALUE )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A     !< Structure holding the matrix
    INTEGER, INTENT(IN) :: i         !< row number of the matrix element
    INTEGER, INTENT(IN) :: j         !< column number of the matrix element
    REAL(KIND=dp), INTENT(IN) :: VALUE   !< value to be added to the matrix element
 !------------------------------------------------------------------------------
    INTEGER :: k
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------
    IF(i>A % NumberOfRows) THEN
      A % FORMAT=MATRIX_LIST; RETURN
    END IF

    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag
    Values => A % Values

    IF ( .NOT.ASSOCIATED(Diag) .OR. i /= j .OR. .NOT. A % Ordered ) THEN
      k = CRS_Search( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),j )
      IF ( k==0 .AND. VALUE/=0 ) A % FORMAT = MATRIX_LIST
      IF ( k==0 ) RETURN
      k = k + Rows(i) - 1
    ELSE
      k = Diag(i)
    END IF
!$omp atomic
    Values(k) = Values(k) + VALUE
  END SUBROUTINE CRS_AddToMatrixElement
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Set a given value to an element of a  CRS format matrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_SetMatrixElement( A,i,j,VALUE )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A     !< Structure holding the matrix
    INTEGER, INTENT(IN) :: i         !< row number of the matrix element
    INTEGER, INTENT(IN) :: j         !< column number of the matrix element
    REAL(KIND=dp), INTENT(IN) :: VALUE   !< new value of the matrix element
!------------------------------------------------------------------------------ 
    INTEGER :: k
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------

    IF(i>A % NumberOfRows) THEN
      A % FORMAT=MATRIX_LIST; RETURN
    END IF

    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag
    Values => A % Values

    IF ( .NOT.ASSOCIATED(Diag) .OR. i /= j .OR. .NOT. A % Ordered ) THEN
       k = CRS_Search( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),j )
      IF ( k==0 ) THEN
         A % FORMAT = MATRIX_LIST
         RETURN
       END IF
       k = k + Rows(i) - 1
    ELSE
       k = Diag(i)
    END IF
    Values(k) = VALUE
  END SUBROUTINE CRS_SetMatrixElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Get a given matrix entry from CRS format matrix.
!------------------------------------------------------------------------------
  FUNCTION CRS_GetMatrixElement( A,i,j ) RESULT ( VALUE )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN):: A     !< Structure holding the matrix
    INTEGER, INTENT(IN) :: i         !< row number of the matrix element
    INTEGER, INTENT(IN) :: j         !< column number of the matrix element
    REAL(KIND=dp) :: VALUE   !< obtained value of the matrix element
!------------------------------------------------------------------------------ 
    INTEGER :: k
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------
    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag
    Values => A % Values

    Value = REAL(0,dp)
    IF ( .NOT.ASSOCIATED(Diag).OR.i /= j .OR. .NOT. A % Ordered ) THEN
       k = CRS_Search( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),j )
      IF ( k==0 ) THEN
         PRINT*,'Trying to get value to nonexistent matrix element: ', i,j
         RETURN
       END IF
       k = k + Rows(i) - 1
    ELSE
       k = Diag(i)
    END IF
    VALUE = Values(k)

  END FUNCTION CRS_GetMatrixElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Get a given matrix entry from CRS format matrix and replace it with a new value
!------------------------------------------------------------------------------
  FUNCTION CRS_ChangeMatrixElement( A,i,j, NewValue ) RESULT ( OldValue )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN):: A     !< Structure holding the matrix
    INTEGER, INTENT(IN) :: i         !< row number of the matrix element
    INTEGER, INTENT(IN) :: j         !< column number of the matrix element
    REAL(KIND=dp), INTENT(IN) :: NewValue  !< Value to be set   
    REAL(KIND=dp) :: OldValue !< Value to be gotten  
!------------------------------------------------------------------------------
 
    INTEGER :: k
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------

    Rows   => A % Rows
    Cols   => A % Cols
    Diag   => A % Diag
    Values => A % Values

    OldValue = REAL(0, dp)
    IF ( .NOT.ASSOCIATED(Diag).OR.i /= j .OR. .NOT. A % Ordered ) THEN
       k = CRS_Search( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1),j )
      IF ( k==0 ) THEN
         PRINT*,'Trying to change value of a nonexistent matrix element: ', i,j,NewValue
         RETURN
       END IF
       k = k + Rows(i) - 1
    ELSE
       k = Diag(i)
    END IF
!$omp critical
    OldValue = Values(k)
    Values(k) = NewValue
!$omp end critical

  END FUNCTION CRS_ChangeMatrixElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Add a row together with another row of a CRS matrix, and thereafter zero it.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_MoveRow( A,n1,n2,coeff )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A    !< Structure holding the matrix
    INTEGER, INTENT(IN) :: n1         !< Row number to be copied and zerod
    INTEGER, INTENT(IN) :: n2         !< Row number to be added
    REAL(KIND=dp),OPTIONAL :: coeff   !< Optional coefficient to multiply the row to be copied with
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: VALUE, c
    INTEGER :: i,j

    IF( PRESENT(Coeff)) THEN
      c = coeff
    ELSE
      c = 1.0_dp
    END IF

    DO i=A % Rows(n1),A % Rows(n1+1)-1
      j = A % Cols(i)
      VALUE = c * A % Values(i) 
      A % Values(i) = 0.0_dp
      CALL CRS_AddToMatrixElement( A,n2,j,VALUE )      
    END DO

  END SUBROUTINE CRS_MoveRow
!------------------------------------------------------------------------------
  


!------------------------------------------------------------------------------
!>    Add a set of values (.i.e. element stiffness matrix) to a CRS format
!>    matrix. 
!------------------------------------------------------------------------------
  SUBROUTINE CRS_GlueLocalMatrix( A,N,Dofs,Indeces,LocalMatrix )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding matrix
     REAL(KIND=dp), INTENT(IN) :: LocalMatrix(:,:)  !< A (N x Dofs) x ( N x Dofs) matrix holding the values to be added to the CRS format matrix
     INTEGER, INTENT(IN) :: N             !< Number of nodes in element
     INTEGER, INTENT(IN) :: Dofs          !< Number of degrees of freedom for one node
     INTEGER, INTENT(IN) :: Indeces(:)    !< Maps element node numbers to global (or partition) node numbers 
     INTEGER::jc,ic
	                                      !! (to matrix rows and columns, if Dofs = 1)
!------------------------------------------------------------------------------ 
     INTEGER :: i,j,k,l,c,Row,Col
     INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
     REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------

     Diag   => A % Diag
     Rows   => A % Rows
     Cols   => A % Cols
     Values => A % Values

     IF ( Dofs == 1 ) THEN
       DO i=1,N
         Row = Indeces(i)
         IF ( Row <=0 ) CYCLE
         DO j=1,N
           Col = Indeces(j)
           IF ( Col <= 0 ) CYCLE
           IF ( Col >= Row ) THEN
             DO c=Diag(Row),Rows(Row+1)-1
               IF ( Cols(c) == Col ) THEN
!$omp atomic
                 Values(c) = Values(c) + LocalMatrix(i,j)
                 EXIT
               END IF
             END DO
           ELSE
             DO c=Rows(Row),Diag(Row)-1
               IF ( Cols(c) == Col ) THEN
!$omp atomic
                 Values(c) = Values(c) + LocalMatrix(i,j)
                 EXIT
               END IF
             END DO
           END IF
         END DO
       END DO
     ELSE
       DO i=1,N
          DO k=0,Dofs-1
             IF ( Indeces(i) <= 0 ) CYCLE
             Row  = Dofs * Indeces(i) - k
             DO j=1,N
                DO l=0,Dofs-1
                   IF ( Indeces(j) <= 0 ) CYCLE
                   Col  = Dofs * Indeces(j) - l
                   IF ( Col >= Row ) THEN
                     DO c=Diag(Row),Rows(Row+1)-1
                        IF ( Cols(c) == Col ) THEN
!$omp atomic
                           Values(c) = Values(c) + LocalMatrix(Dofs*i-k,Dofs*j-l)
                           EXIT
                        END IF
                     END DO
                   ELSE
                     DO c=Rows(Row),Diag(Row)-1
                        IF ( Cols(c) == Col ) THEN
!$omp atomic
                           Values(c) = Values(c) + LocalMatrix(Dofs*i-k,Dofs*j-l)
                           EXIT
                        END IF
                     END DO
                   END IF
                END DO
             END DO
          END DO
       END DO
     END IF

  END SUBROUTINE CRS_GlueLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Add a set of values (.i.e. element stiffness matrix) to a CRS format
!>    matrix. For this matrix the entries are ordered so that first for one
!>    dof you got all nodes, and then for second etc. There may be on offset
!>    to the entries making the subroutine suitable for coupled monolithic
!>    matrix assembly.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_GlueLocalSubMatrix( A,row0,col0,Nrow,Ncol,RowInds,ColInds,&
                  RowDofs,ColDofs,LocalMatrix )
!------------------------------------------------------------------------------ 
     REAL(KIND=dp), INTENT(IN) :: LocalMatrix(:,:)  !< A (Nrow x RowDofs) x ( Ncol x ColDofs) matrix holding the values to be
                                                    !!            added to the CRS format matrix
     TYPE(Matrix_t) :: A           !< Structure holding matrix
     INTEGER, INTENT(IN) :: Nrow   !< Number of nodes for the row variable
     INTEGER, INTENT(IN) :: Ncol   !< Number of nodes for the column variable
     INTEGER, INTENT(IN) :: RowDofs  !< Number of dofs for row variable
     INTEGER, INTENT(IN) :: ColDofs  !< Number of dofs for column variable
     INTEGER, INTENT(IN) :: Col0     !< Offset for column variable
	 INTEGER, INTENT(IN) :: Row0     !< Offset for row variable
     INTEGER, INTENT(IN) :: RowInds(:)  !< Permutation of the row dofs
	 INTEGER, INTENT(IN) :: ColInds(:)  !< Permutation of the column dofs
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,c,Row,Col
     INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
     REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------
     Diag   => A % Diag
     Rows   => A % Rows
     Cols   => A % Cols
     Values => A % Values
     
     DO i=1,Nrow
        DO k=0,RowDofs-1
           IF ( RowInds(i) <= 0 ) CYCLE
           Row  = Row0 + RowDofs * RowInds(i) - k

           DO j=1,Ncol
              DO l=0,ColDofs-1
                 IF ( ColInds(j) <= 0 ) CYCLE
                 Col  = Col0 + ColDofs * ColInds(j) - l

! If Diag does not exist then one cannot separate the gluing into two parts
! In fact this cannot be guarateed for the off-diagonal block matrices and hence 
! this condition is set active.
                IF( .TRUE. ) THEN
                    DO c=Rows(Row),Rows(Row+1)-1
                       IF ( Cols(c) == Col ) THEN
!$omp atomic
                           Values(c) = Values(c) + LocalMatrix(RowDofs*i-k,ColDofs*j-l)
                           EXIT
                        END IF
                     END DO
                     IF( Cols(c) /= Col ) PRINT *,'NO HIT 1',Row,Col
                 ELSE IF ( Col >= Row ) THEN
                    DO c=Diag(Row),Rows(Row+1)-1
                       IF ( Cols(c) == Col ) THEN
!$omp atomic
                          Values(c) = Values(c) + LocalMatrix(RowDofs*i-k,ColDofs*j-l)
                          EXIT
                       END IF
                    END DO
                    IF( Cols(c) /= Col ) PRINT *,'NO HIT 2',Row,Col
                 ELSE
                    DO c=Rows(Row),Diag(Row)-1
                       IF ( Cols(c) == Col ) THEN
!$omp atomic
                           Values(c) = Values(c) + LocalMatrix(RowDofs*i-k,ColDofs*j-l)
                           EXIT
                        END IF
                     END DO
                     IF( Cols(c) /= Col ) PRINT *,'NO HIT 3',Row,Col
                  END IF
               END DO
            END DO
         END DO
      END DO
    END SUBROUTINE CRS_GlueLocalSubMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> When Dirichlet conditions are set by zeroing the row except for setting 
!> the diagonal entry to one, the matrix symmetry is broken. This routine 
!> maintains the symmetric structure of the matrix equation.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_SetSymmDirichlet( A,b,n,val )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A       !< Structure holding matrix
    INTEGER, INTENT(IN) :: n              !< Index of the dofs to be fixed   
    REAL(KIND=dp) :: b(:)     !< right-hand-side of the matrix equation
    REAL(KIND=dp), INTENT(IN) :: val      !< Dirichlet value to be set
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,k1,k2
    REAL(KIND=dp) :: t
    LOGICAL :: isMass, isDamp

    isMass = ASSOCIATED(A % MassValues)
    IF ( isMass ) &
      isMass = isMass .AND. SIZE(A % MassValues) == SIZE(A % Values)

    isDamp = ASSOCIATED(A % DampValues)
    IF ( isDamp ) &
      isDamp = isDamp .AND. SIZE(A % DampValues) == SIZE(A % Values)

    IF(.NOT.A % NoDirichlet) THEN
      DO l=A % Rows(n),A % Rows(n+1)-1
         i = A % Cols(l)
         IF ( i == n ) CYCLE

         IF ( n > i ) THEN
           k1 = A % Diag(i)+1
           k2 = A % Rows(i+1)-1
         ELSE 
           k1 = A % Rows(i)
           k2 = A % Diag(i)-1
         END IF

         k = k2 - k1 + 1
         IF ( k <= 30 ) THEN
           DO j = k1, k2
             IF ( A % Cols(j) == n ) THEN
               b(i) = b(i) - A % Values(j) * val
               A % Values(j) = 0.0_dp
               IF ( isMass ) A % MassValues(j) = 0._dp
               IF ( isDamp ) A % DampValues(j) = 0._dp
               EXIT
             ELSE IF ( A % Cols(j) > n ) THEN
               EXIT
             END IF
           END DO
         ELSE
           j = CRS_Search( k,A % Cols(k1:k2),n )
           IF ( j > 0 ) THEN
             j = j + k1 - 1
             b(i) = b(i) - A % Values(j) * val
             A % Values(j) = 0.0_dp
             IF ( isMass ) A % MassValues(j) = 0._dp
             IF ( isDamp ) A % DampValues(j) = 0._dp
           END IF
         END IF
      END DO
      CALL CRS_ZeroRow(A,n)
      A % Values(A % Diag(n)) = 1._dp
      b(n) = val
    END IF

    IF(ALLOCATED(A % Dvalues)) THEN
      A % DValues(n) = val
    ELSE
      b(n) = val
    END IF
    IF(ALLOCATED(A % ConstrainedDOF)) A % ConstrainedDOF(n) = .TRUE.
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_SetSymmDirichlet
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes the rowsoum of a given row in a CRS matrix.
!------------------------------------------------------------------------------
FUNCTION CRS_RowSum( A,k ) RESULT(rsum)
!------------------------------------------------------------------------------
   TYPE(Matrix_t), INTENT(IN) :: A       !< Structure holding matrix
   INTEGER, INTENT(IN) :: k              !< Row index 
   REAL(KIND=dp) :: rsum                 !< Sum of the entries on the row
!------------------------------------------------------------------------------
   INTEGER :: i
   
   rsum = 0.0D0
   DO i=A % Rows(k), A % Rows(k+1)-1
     rsum  = rsum + A % Values( i )
   END DO
!------------------------------------------------------------------------------
END FUNCTION CRS_RowSum
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes information on the matrix rowsums.
!------------------------------------------------------------------------------
SUBROUTINE CRS_RowSumInfo( A, Values ) 
!------------------------------------------------------------------------------
   TYPE(Matrix_t), INTENT(IN) :: A       
   REAL(KIND=dp), POINTER, OPTIONAL :: Values(:)
!------------------------------------------------------------------------------
   REAL(KIND=dp), POINTER :: PValues(:)
   INTEGER :: i,j,k              
   REAL(KIND=dp) :: val,rsum,absrsum     
   REAL(KIND=dp) :: minrsum,maxrsum,minabsrsum,maxabsrsum
!------------------------------------------------------------------------------

   minrsum = HUGE(minrsum)
   maxrsum = -HUGE(maxrsum)
   minabsrsum = HUGE(minabsrsum)
   maxabsrsum = 0.0_dp
   
   IF( PRESENT( Values ) ) THEN
     PValues => Values
   ELSE
     PValues => A % Values
   END IF

   DO i=1,A % NumberOfRows
     rsum = 0.0_dp
     absrsum = 0.0_dp

     DO j=A % Rows(i), A % Rows(i+1)-1
       val = PValues( j )
       rsum = rsum + val
       absrsum = absrsum + ABS( val ) 
     END DO

     minrsum = MIN( minrsum, rsum ) 
     maxrsum = MAX( maxrsum, rsum ) 
     minabsrsum = MIN( minabsrsum, absrsum ) 
     maxabsrsum = MAX( maxabsrsum, absrsum ) 
   END DO

   WRITE( Message,'(A,ES12.4)') 'Total sum:',SUM( PValues )
   CALL Info( 'CRS_RowSumInfo', Message ) 
   WRITE( Message,'(A,2ES12.4)') 'Rowsum range:',minrsum,maxrsum
   CALL Info( 'CRS_RowSumInfo', Message ) 
   WRITE( Message,'(A,2ES12.4)') 'Absolute rowsum range:',minabsrsum,maxabsrsum
   CALL Info( 'CRS_RowSumInfo', Message ) 
   
!------------------------------------------------------------------------------
 END SUBROUTINE CRS_RowSumInfo
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Create the structures required for a CRS format matrix.
!------------------------------------------------------------------------------
  FUNCTION CRS_CreateMatrix( N,Total,RowNonzeros,Ndeg,Reorder,AllocValues ) RESULT(A)
!------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: N    !< Number of rows in the matrix
    INTEGER, INTENT(IN) :: Total  !< Total number of nonzero entries in the matrix
    INTEGER, INTENT(IN) :: Ndeg   !< Negrees of freedom
    INTEGER, INTENT(IN) :: RowNonzeros(:)  !< Number of nonzero entries in rows of the matrix
    INTEGER, INTENT(IN) :: Reorder(:)      !< Permutation index for bandwidth reduction
    LOGICAL, INTENT(IN) :: AllocValues     !< Should the values arrays be allocated ?
    TYPE(Matrix_t), POINTER :: A  !>  Pointer to the created Matrix_t structure.
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,istat
    INTEGER, POINTER :: InvPerm(:)
!------------------------------------------------------------------------------

    A => AllocateMatrix()

    k = Ndeg*Ndeg*Total

    ALLOCATE( A % Rows(n+1),A % Diag(n),A % Cols(k),STAT=istat )

    IF ( istat == 0 .AND. AllocValues ) THEN
      ALLOCATE( A % Values(k), STAT=istat )
    END IF

    IF ( istat /= 0 ) THEN
      CALL Fatal( 'CreateMatrix', 'Memory allocation error.' )
    END IF

    NULLIFY( A % ILUValues )
    NULLIFY( A % CILUValues )

    InvPerm => A % Diag ! just available memory space...
    j = 0
    DO i=1,SIZE(Reorder)
       IF ( Reorder(i) > 0 ) THEN
          j = j + 1
          InvPerm(Reorder(i)) = j
       END IF
    END DO

    A % NumberOfRows = N
    A % Rows(1) = 1
    DO i=2,n
       j = InvPerm((i-2)/Ndeg+1)
       A % Rows(i) = A % Rows(i-1) + Ndeg*RowNonzeros(j)
    END DO
    j = InvPerm((n-1)/ndeg+1)
    A % Rows(n+1) = A % Rows(n)  +  Ndeg*RowNonzeros(j)

    A % Cols = 0
    A % Diag = 0

    A % Ordered = .FALSE.
  END FUNCTION CRS_CreateMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Matrix vector product (v = Au) for a matrix given in CRS format.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_MatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: u   !< Vector to be multiplied
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: v  !< Result vector
    TYPE(Matrix_t), INTENT(IN) :: A                !< Structure holding matrix
!------------------------------------------------------------------------------
     INTEGER, POINTER  CONTIG :: Cols(:),Rows(:)
     REAL(KIND=dp), POINTER  CONTIG :: Values(:)

     INTEGER :: i,j,n
     REAL(KIND=dp) :: rsum
#ifdef HAVE_MKL
	INTERFACE
		SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
	 		USE Types
	 		CHARACTER :: transa
	 		INTEGER :: m
	 		REAL(KIND=dp) :: a(*)
	 		INTEGER :: ia(*), ja(*)
	 		REAL(KIND=dp) :: x(*), y(*)
	 	END SUBROUTINE mkl_dcsrgemv
	END INTERFACE
#endif

!------------------------------------------------------------------------------

     n = A % NumberOfRows
     Rows   => A % Rows
     Cols   => A % Cols
     Values => A % Values

    IF  ( A % MatvecSubr /= 0 ) THEN
#ifdef USE_ISO_C_BINDINGS
      CALL MatVecSubrExt(A % MatVecSubr,A % SpMV, n,Rows,Cols,Values,u,v,0)
#else
      CALL MatVecSubr(A % MatVecSubr,A % SpMV, n,Rows,Cols,Values,u,v,0)
#endif
      RETURN
   END IF

	! Use MKL to perform mvp if it is available
#ifdef HAVE_MKL
	CALL mkl_dcsrgemv('N', n, Values, Rows, Cols, u, v)
#else
!$omp parallel do private(j,rsum)
     DO i=1,n
        rsum = 0.0d0
        DO j=Rows(i),Rows(i+1)-1
           rsum = rsum + u(Cols(j)) * Values(j)
        END DO
        v(i) = rsum
     END DO
!$omp end parallel do
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_MatrixVectorMultiply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Matrix-vector product v = |A|u with A a matrix in the CRS format and
!>  |.| the matrix function giving the absolute values of the argument 
!>  components. This special mv subroutine may be needed in connection with
!>  certain stopping criteria for iterative linear solvers. 
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ABSMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: u   !< The vector u
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: v  !< The result vector v
    TYPE(Matrix_t), INTENT(IN) :: A                !< The structure holding the matrix A
!------------------------------------------------------------------------------
    INTEGER, POINTER  CONTIG :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER  CONTIG :: Values(:), Abs_Values(:)

    
    INTEGER :: i,j,n
    REAL(KIND=dp) :: rsum
#ifdef HAVE_MKL
    INTERFACE
      SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
        USE Types
        CHARACTER :: transa
        INTEGER :: m
        REAL(KIND=dp) :: a(*)
        INTEGER :: ia(*), ja(*)
        REAL(KIND=dp) :: x(*), y(*)
      END SUBROUTINE mkl_dcsrgemv
    END INTERFACE
#endif
!------------------------------------------------------------------------------

    n = A % NumberOfRows
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    IF  ( A % MatvecSubr /= 0 ) THEN
#ifdef USE_ISO_C_BINDINGS
      ALLOCATE(Abs_Values(SIZE(A % Values)))
      Abs_Values = ABS(Values)
      CALL MatVecSubrExt(A % MatVecSubr,A % SpMV, n,Rows,Cols,Abs_Values,u,v,0) ! TODO: (bug) must be ABS(Values)
      DEALLOCATE(Abs_Values)
#else
      CALL MatVecSubr(A % MatVecSubr,A % SpMV, n,Rows,Cols,ABS(Values),u,v,0)
#endif
      RETURN
    END IF

    ! Use MKL to perform mvp if it is available
#ifdef HAVE_MKL
    CALL mkl_dcsrgemv('N', n, ABS(Values), Rows, Cols, u, v)
#else
!$omp parallel do private(j,rsum)
    DO i=1,n
      rsum = 0.0d0
      DO j=Rows(i),Rows(i+1)-1
        rsum = rsum + u(Cols(j)) * ABS(Values(j))
      END DO
      v(i) = rsum
    END DO
!$omp end parallel do
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_ABSMatrixVectorMultiply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Calculate transpose of A in CRS format: B = A^T
!------------------------------------------------------------------------------
     FUNCTION CRS_Transpose( A ) RESULT(B)
!------------------------------------------------------------------------------
       IMPLICIT NONE
       
       TYPE(Matrix_t), POINTER :: A, B
       
       INTEGER, ALLOCATABLE :: Row(:)
       INTEGER :: NVals
       INTEGER :: i,j,k,istat,n

       B => AllocateMatrix()
       
       NVals = SIZE( A % Values )
       B % NumberOfRows = MAXVAL( A % Cols )

       ALLOCATE( B % Rows( B % NumberOfRows +1 ), B % Cols( NVals ), &
           B % Values( Nvals ), B % Diag( B % NumberOfRows ), Row( B % NumberOfRows ), &
           STAT=istat )
       IF ( istat /= 0 )  CALL Fatal( 'CRS_Transpose', &
           'Memory allocation error.' )

       B % Diag = 0
       Row = 0       
       DO i = 1, NVals
         Row( A % Cols(i) ) = Row( A % Cols(i) ) + 1
       END DO
       
       B % Rows(1) = 1
       DO i = 1, B % NumberOfRows
         B % Rows(i+1) = B % Rows(i) + Row(i)
       END DO
       B % Cols = 0
       
       DO i = 1, B % NumberOfRows
         Row(i) = B % Rows(i)
       END DO
      
       DO i = 1, A % NumberOfRows

         DO j = A % Rows(i), A % Rows(i+1) - 1
           k = A % Cols(j)

           IF ( Row(k) < B % Rows(k+1) ) THEN 
             B % Cols( Row(k) ) = i
             B % Values( Row(k) ) = A % Values(j)
             Row(k) = Row(k) + 1
           ELSE
             WRITE( Message, * ) 'Trying to access non-existent column', i,k,j
             CALL Error( 'CRS_Transpose', Message )
             RETURN
           END IF
         END DO
       END DO

       DEALLOCATE( Row )

!------------------------------------------------------------------------------
     END FUNCTION CRS_Transpose
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Matrix vector product (v = A^T u) for a matrix given in CRS format.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_TransposeMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: u   !< Vector to be multiplied
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: v  !< Result vector
    TYPE(Matrix_t), INTENT(IN) :: A                !< Structure holding matrix
!------------------------------------------------------------------------------
     INTEGER, POINTER  CONTIG :: Cols(:),Rows(:)
     REAL(KIND=dp), POINTER  CONTIG :: Values(:)

     INTEGER :: i,j,k,n
     REAL(KIND=dp) :: rsum
#ifdef HAVE_MKL
	INTERFACE
		SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
	 		USE Types
	 		CHARACTER :: transa
	 		INTEGER :: m
	 		REAL(KIND=dp) :: a(*)
	 		INTEGER :: ia(*), ja(*)
	 		REAL(KIND=dp) :: x(*), y(*)
	 	END SUBROUTINE mkl_dcsrgemv
	END INTERFACE
#endif
!------------------------------------------------------------------------------

     n = A % NumberOfRows
     Rows   => A % Rows
     Cols   => A % Cols
     Values => A % Values

	! Use MKL to perform mvp if it is available
#ifdef HAVE_MKL
	CALL mkl_dcsrgemv('T', n, Values, Rows, Cols, u, v)
#else
     v(1:n) = 0.0_dp
     DO i=1,n
       DO j=Rows(i),Rows(i+1)-1
         k = Cols(j)
         v(k) = v(k) + u(i) * Values(j)
       END DO
     END DO
#endif
!------------------------------------------------------------------------------
   END SUBROUTINE CRS_TransposeMatrixVectorMultiply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Add another matrix B to matrix A, and eliminate B
!------------------------------------------------------------------------------
  SUBROUTINE CRS_MergeMatrix( A,B,PermA,PermB )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A             !< Structure holding the master matrix
    TYPE(Matrix_t), POINTER :: B             !< Structure holding the slave matrix
    INTEGER, POINTER, OPTIONAL :: PermA(:)   !< Permutation of the master dofs
    INTEGER, POINTER, OPTIONAL :: PermB(:)   !< Permutation of the slave dofs
!------------------------------------------------------------------------------
     INTEGER, POINTER  CONTIG :: ColsA(:),RowsA(:),ColsB(:),RowsB(:),Rows(:),Cols(:)
     REAL(KIND=dp), POINTER  CONTIG :: ValuesA(:),ValuesB(:),Values(:)
     INTEGER :: i,j,k,n,nA,nB,kb,iA,iB
     LOGICAL :: Set,UsePerm
!------------------------------------------------------------------------------
     REAL(kind=dp) :: sumA, sumB

     CALL Info('CRS_MergeMatrix','Merging two matrices',Level=9)

     IF(.NOT. ASSOCIATED(A)) THEN
       CALL Fatal('CRS_MergeMatrix','A not associated')
     ELSE IF(.NOT. ASSOCIATED(B)) THEN
       CALL Fatal('CRS_MergeMatrix','B not associated')
     END IF

     UsePerm = PRESENT( PermA ) 

     IF( UsePerm ) THEN
       IF(.NOT. PRESENT( PermB ) ) THEN
         CALL Fatal('CRS_MergeMatrix','Either both PermA and PermB or neither')
       END IF
        
       n = SIZE(PermA)
       IF( SIZE(PermB) /= n ) THEN
         CALL Fatal('CRS_MergeMatrix','Mismatch in perm size')
       END IF
     ELSE
       n = A % NumberOfRows
       IF( n /= B % NumberOfRows ) THEN
         CALL Fatal('CRS_MergeMatrix','Mismatch in matrix size')
       END IF
     END IF

     RowsA   => A % Rows
     ColsA   => A % Cols
     ValuesA => A % Values

     RowsB   => B % Rows
     ColsB   => B % Cols
     ValuesB => B % Values

     Set = .FALSE.
     sumA = 0.0_dp
     sumB = 0.0_dp

100  kb = 0
     IF( UsePerm ) THEN
       DO i=1,n
         iA = PermA(i)
         iB = PermB(i)
         IF( iA > 0 .AND. iB > 0 ) THEN
           WRITE (Message,'(A,I0,I0)') 'Code the possibility to merge rows: ',iA,iB
           CALL Fatal('CRS_MergeMatrix',Message)
         END IF

         IF( iA > 0 ) THEN
           nA = RowsA(iA+1)-RowsA(iA) 
         ELSE
           nA = 0
         END IF
         IF( iB > 0 ) THEN
           nB = RowsB(iB+1)-RowsB(iB)
         ELSE
           nB = 0
         END IF

         IF( nA > 0 ) THEN           
           DO j=RowsA(iA),RowsA(iA+1)-1
             kb = kb + 1
             IF( Set ) THEN
               Cols(kb) = ColsA(j)
               Values(kb) = ValuesA(j)
               sumA = sumA + ValuesA(j)
             END IF
           END DO
         ELSE IF( nB > 0 ) THEN
           DO j=RowsB(iB),RowsB(iB+1)-1
             kb = kb + 1
             IF( Set ) THEN
               Cols(kb) = ColsB(j)
               Values(kb) = ValuesB(j)
               sumB = sumB + ValuesB(j)
             END IF
           END DO
         END IF
         IF( Set ) THEN
           Rows(i+1) = kb+1
         END IF
       END DO
     ELSE
       DO i=1,n
         nA = RowsA(i+1)-RowsA(i) 
         nB = RowsB(i+1)-RowsB(i)
         IF( nA > 0 .AND. nB > 0 ) THEN
           WRITE (Message,'(A,I0,I0)') 'Code the possibility to merge rows: ',iA,iB
           CALL Fatal('CRS_MergeMatrix',Message)
         ELSE IF( nA > 0 ) THEN
           DO j=RowsA(i),RowsA(i+1)-1
             kb = kb + 1
             IF( Set ) THEN
               Cols(kb) = ColsA(j)
               Values(kb) = ValuesA(j)
             END IF
           END DO
         ELSE IF( nB > 0 ) THEN
           DO j=RowsB(i),RowsB(i+1)-1
             kb = kb + 1
             IF( Set ) THEN
               Cols(kb) = ColsB(j)
               Values(kb) = ValuesB(j)
             END IF
           END DO
         END IF
         IF( Set ) THEN
           Rows(i+1) = kb+1
         END IF
       END DO
     END IF

     IF( kb == 0 ) THEN
       CALL Fatal('CRS_MergeMatrix','Union size is zero?')
     END IF
     
     IF(.NOT. Set) THEN
       ALLOCATE( Rows(n+1), Cols(kb), Values(kb) )
       Rows(1) = 1
       Set = .TRUE.
       CALL Info('CRS_MergeMatrix','Done Allocating and going now really',Level=9)
       GOTO 100
     END IF
     
     DEALLOCATE( RowsA, RowsB, ColsA, ColsB, ValuesA, ValuesB )     
     B % NumberOfRows = 0

     A % Rows => Rows
     A % Cols => Cols
     A % Values => Values
     A % NumberOfRows = n

     CALL Info('CRS_MergeMatrix','Merging of matrices finisged',Level=9)

!------------------------------------------------------------------------------
   END SUBROUTINE CRS_MergeMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Matrix vector product (v = Au) for a matrix given in CRS format
!>    assuming complex valued matrix equation.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN) :: u   !< Vector to be multiplied
    COMPLEX(KIND=dp), DIMENSION(*), INTENT(OUT) :: v  !< Result vector
    TYPE(Matrix_t), INTENT(IN) :: A                !< Structure holding matrix
!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: i,j,n
    COMPLEX(KIND=dp) :: s,rsum
!------------------------------------------------------------------------------
    n = A % NumberOfRows / 2
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

!$omp parallel do private(rsum,j,s)
    DO i=1,n
       rsum = CMPLX( 0.0d0, 0.0d0,KIND=dp )
       DO j=Rows(2*i-1),Rows(2*i)-1,2
          s = CMPLX( Values(j), -Values(j+1), KIND=dp )
          rsum = rsum + s * u((Cols(j)+1)/2)
       END DO
       v(i) = rsum
    END DO
!$omp end parallel do
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_ComplexMatrixVectorMultiply
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  SUBROUTINE CRS_ApplyProjector( PMatrix, u, uperm, v, vperm, Trans )
!-------------------------------------------------------------------------------
    TYPE(Matrix_t) :: PMatrix
    REAL(KIND=dp) :: u(:),v(:)
    LOGICAL, OPTIONAL :: Trans
    INTEGER, POINTER :: uperm(:), vperm(:)
!-------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n
    REAL(KIND=dp), POINTER :: Values(:)
    LOGICAL :: LTrans
    INTEGER, POINTER :: Rows(:), Cols(:)
!-------------------------------------------------------------------------------
    LTrans = .FALSE.
    IF ( PRESENT( Trans ) ) LTrans = Trans

    n = PMatrix % NumberOfRows
    Rows   => PMatrix % Rows
    Cols   => PMatrix % Cols
    Values => PMatrix % Values

    IF ( ASSOCIATED( uperm ) .AND. ASSOCIATED( vperm ) ) THEN
       IF ( LTrans ) THEN
          DO i=1,n
             k = uperm(i)
             IF ( k > 0 ) THEN
                DO j=Rows(i),Rows(i+1)-1
                   l = vperm(Cols(j))
                   IF ( l > 0 ) v(l) = v(l) + u(k) * Values(j)
                END DO
             END IF
          END DO
       ELSE
          DO i=1,n
             l = vperm(i)
             IF ( l > 0 ) THEN
               IF ( ANY(Values(Rows(i):Rows(i+1)-1)/=0) ) v(l)=0
             END IF
          END DO

          DO i=1,n
             l = vperm(i)
             IF ( l > 0 ) THEN
                DO j = Rows(i), Rows(i+1)-1
                   k = uperm(Cols(j))
                   IF ( k>0 ) v(l) = v(l) + u(k) * Values(j)
                END DO
             END IF
          END DO
       END IF
    ELSE
       IF ( LTrans ) THEN
          DO i=1,n
             DO j=Rows(i),Rows(i+1)-1
                v(Cols(j)) = v(Cols(j)) + u(i) * Values(j)
             END DO
          END DO
       ELSE
          DO i=1,n
             DO j = Rows(i), Rows(i+1)-1
                v(i) = v(i) + u(Cols(j)) * Values(j)
             END DO
          END DO
       END IF
    END IF
!-------------------------------------------------------------------------------
  END SUBROUTINE CRS_ApplyProjector
!-------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!>    Diagonal preconditioning of a CRS format matrix. Matrix is accessed
!>    from a global variable GlobalMatrix. Note that if the matrix has been 
!> scaled so that the diagonal entries are already ones, this subroutine is obsolite.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_DiagPrecondition( u,v,ipar )
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: u  !< Resulting approximate solution after preconditioning
    REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: v  !< Given right-hand-side
    INTEGER, DIMENSION(*) :: ipar   !< structure holding info from (HUTIter-iterative solver package)
!------------------------------------------------------------------------------
    INTEGER :: i,j,n
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    REAL(KIND=dp), POINTER :: Values(:)

    Diag   => GlobalMatrix % Diag
    Rows   => GlobalMatrix % Rows
    Cols   => GlobalMatrix % Cols
    Values => GlobalMatrix % Values

    n = GlobalMatrix % NumberOfRows

    IF ( .NOT. GlobalMatrix % Ordered ) THEN
       DO i=1,N
          CALL SortF( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1), &
                   Values(Rows(i):Rows(i+1)-1) )
       END DO
       DO i=1,N
          DO j=Rows(i),Rows(i+1)-1
             IF ( Cols(j) == i ) THEN
                Diag(i) = j
                EXIT
             END IF
          END DO
       END DO
       GlobalMatrix % Ordered = .TRUE.
    END IF

    DO i=1,n
       IF  ( ABS( Values(Diag(i))) > AEPS ) THEN
           u(i) = v(i) / Values(Diag(i))
       ELSE
           u(i) = v(i)
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_DiagPrecondition
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Diagonal preconditioning of a CRS format matrix for complex valued matrix equations. 
!>    Matrix is accessed from a global variable GlobalMatrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexDiagPrecondition( u,v,ipar )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), DIMENSION(*), INTENT(OUT) :: u  !< Resulting approximate solution after preconditioning
    COMPLEX(KIND=dp), DIMENSION(*), INTENT(IN) :: v  !< Given right-hand-side
    INTEGER, DIMENSION(*) :: ipar   !< structure holding info from (HUTIter-iterative solver package)
!------------------------------------------------------------------------------
    INTEGER :: i,j,n
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    COMPLEX(KIND=dp) :: A
    REAL(KIND=dp), POINTER :: Values(:)

    Diag   => GlobalMatrix % Diag
    Rows   => GlobalMatrix % Rows
    Cols   => GlobalMatrix % Cols
    Values => GlobalMatrix % Values

    n = GlobalMatrix % NumberOfRows

    IF ( .NOT. GlobalMatrix % Ordered ) THEN
       DO i=1,N
          CALL SortF( Rows(i+1)-Rows(i),Cols(Rows(i):Rows(i+1)-1), &
                   Values(Rows(i):Rows(i+1)-1) )
       END DO

       DO i=1,N
          DO j=Rows(i),Rows(i+1)-1
             IF ( Cols(j) == i ) THEN
                Diag(i) = j
                EXIT
             END IF
          END DO
       END DO
       GlobalMatrix % Ordered = .TRUE.
    END IF

    DO i=1,n/2
       A = CMPLX( Values(Diag(2*i-1)), -Values(Diag(2*i-1)+1), KIND=dp )
       u(i) = v(i) / A
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE CRS_ComplexDiagPrecondition
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Pics the block diagonal entries from matrix A to build matrix B.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_BlockDiagonal(A,B,Blocks) 
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A      !< The initial matrix
    TYPE(Matrix_t) :: B  !< The block diagonal matrix
    INTEGER, INTENT(IN) :: Blocks        !< Number of blocks used in the decomposition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,kb,n

    IF(Blocks <= 1) RETURN

    N = A % NumberOfRows
    B % NumberOfRows = N
    
    kb = 0
    DO i=1,N
      DO k= A % Rows(i), A % Rows(i+1)-1
        l = A % Cols(k)
        IF( MOD(i,Blocks) == MOD(l,Blocks)) kb = kb + 1
      END DO
    END DO
    ALLOCATE(B % Rows(N+1),B % Cols(kb), B % Values(kb), B % Diag(n))
      
    kb = 1
    DO i=1,N
      B % Rows(i) = kb
      DO k = A % Rows(i), A % Rows(i+1)-1
        l = A % Cols(k)
        IF( MOD(i,Blocks) == MOD(l,Blocks) ) THEN
          B % Values(kb) = A % Values(k)
          B % Cols(kb) = A % Cols(k)
          IF( B % Cols(kb) == i) B % Diag(i) = kb
          kb = kb + 1  
        END IF
      END DO
    END DO
    B % Rows(N+1) = kb

!------------------------------------------------------------------------------
  END SUBROUTINE CRS_BlockDiagonal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Removes zeros from the matrix structure.
!> This might be done in order to save memory, or to speed up the matrix 
!> operations. One must be carefull since the fact the an entry is zero
!> does not always imply that it would be zero throughout the simulation.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_RemoveZeros( A  )
!-------------------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A     !< The matrix which will be returned with the non-zeros removed
!-------------------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,kb,kb0,n,rowkb
    INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
    REAL(KIND=DP) :: val
    REAL(KIND=DP), POINTER :: Values(:)

    N = A % NumberOfRows
    
    ! Count the number of nonzeros
    ! The diagonal entry is assumed always to exist.
    kb = 0
    DO i=1,N
      DO k= A % Rows(i), A % Rows(i+1)-1
        l = A % Cols(k)
        val = A % Values(k)
        IF( i == l .OR. ABS( val ) > TINY( val) ) THEN
          kb = kb + 1
        END IF
      END DO
    END DO
    
    kb0 = SIZE( A % Values )
    IF( kb == kb0 ) THEN
      CALL Info('CRS_RemoveZeros','There are no zeros to remove',Level=6)
      RETURN
    END IF

    WRITE( Message,'(A,F8.3,A)') 'Fraction of zeros to remove: ',100.0*(1.0-1.0*kb/kb0),' %'
    CALL Info('CRS_RemoveZeros', Message)

    ! These are new
    ALLOCATE(Cols(kb), Values(kb))
   
    ! These are overwritten
    Diag => A % Diag
    Rows => A % Rows


    kb = 1
    DO i=1,N
      ! Memorize this as this should not be set before the next loop
      rowkb = kb
      DO k = A % Rows(i), A % Rows(i+1)-1
        l = A % Cols(k)
        val = A % Values(k) 

        IF( i == l ) THEN
          Diag(i) = kb
        ELSE IF( .NOT. ( ABS(val) > TINY(val) ) ) THEN
          CYCLE
        END IF

        ! Set the new entry to the matrix
        Values(kb) = A % Values(k)
        Cols(kb) = l
        kb = kb + 1

      END DO
      Rows(i) = rowkb
    END DO
    Rows(N+1) = kb

    DEALLOCATE( A % Values, A % Cols ) 
    A % Values => Values
    A % Cols => Cols

!------------------------------------------------------------------------------
  END SUBROUTINE CRS_RemoveZeros
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
!> Makes a algebraic lower order scheme assuming steady state advection-diffusion equation.
!> This can be applied together with flux corrected transport (FCT) scheme.
!> Also creates a lumped mass to MassValuesLumped and saves the original stiffness
!> matrix values to BulkValues.
!
!> For more information see, for example, 
!> Dmitri Kuzmin (2008): "Explicit and implicit FEM-FCT algorithms with flux linearization"
!------------------------------------------------------------------------------
  SUBROUTINE CRS_FCTLowOrder( A )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A           !< Initial higher order matrix
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,k2,l,kb,n
    REAL(KIND=dp) :: Aij,Aji,Aii,Dij,dFij,msum
    REAL(KIND=dp), POINTER :: ML(:)
    INTEGER, POINTER :: Rows(:), Cols(:),Diag(:)
    LOGICAL :: Found, Positive
    INTEGER :: pcount

    CALL Info('CRS_FCTLowOrder','Making low order FCT correction to matrix',Level=5)

    N = A % NumberOfRows
    Rows => A % Rows
    Cols => A % Cols
    Diag => A % Diag

    IF(.NOT. ASSOCIATED(A % FCT_D) ) THEN
      ALLOCATE( A % FCT_D(SIZE( A % Values) ) ) 
    END IF
    A % FCT_D = 0.0_dp    

    IF(.NOT. ASSOCIATED(A % BulkValues) ) THEN
      ALLOCATE( A % BulkValues(SIZE( A % Values) ) ) 
    END IF
    A % BulkValues = A % Values

    pcount = 0
    Positive = .TRUE.

    DO i=1,N

      Aii = A % Values(Diag(i))
      IF( Aii > 0.0_dp ) pcount = pcount + 1

      DO k = Rows(i), Rows(i+1)-1
        j = Cols(k)

        ! Go through the lower symmetric part (ij) and find the corresponding entry (ji)
        IF( i >= j ) CYCLE

        ! First find entry (j,i)
        Found = .FALSE.
        DO k2 = Rows(j), Rows(j+1)-1
          IF( Cols(k2) == i ) THEN
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. Found ) THEN
          CALL Fatal('CRS_FCTLowOrder','Entry not found, matrix topology not symmetric!?')
        END IF

        ! Formula (30) in Kuzmin's paper
        ! Modified so that it also works similarly to A and -A
        Aij = A % Values(k)
        Aji = A % Values(k2)

        ! Positive = Aii > 0.0_dp
        ! In Kuzmin's paper matrix K is -K compared to Elmer convention.
        ! Hence also the condition here is opposite. 
        IF( Positive ) THEN
          Dij = MIN( -Aij, -Aji, 0.0_dp )
        ELSE
          Dij = MAX( -Aij, -Aji, 0.0_dp )
        END IF

        IF(.FALSE.) THEN
          PRINT *,'ij',i,j,Cols(k2),Cols(k)
          PRINT *,'Diag',Cols(Diag(i)),Cols(Diag(j))
          PRINT *,'A',Aij,Aji,Aii,Dij
        END IF

        ! Formula (32) in Kuzmin's paper
        IF( ABS( Dij ) > 0.0_dp ) THEN
          A % FCT_D(k) = A % FCT_D(k) + Dij
          A % FCT_D(k2) = A % FCT_D(k2) + Dij
          A % FCT_D(Diag(i)) = A % FCT_D(Diag(i)) - Dij
          A % FCT_D(Diag(j)) = A % FCT_D(Diag(j)) - Dij          
        END IF
        
      END DO
    END DO

    WRITE( Message,'(A,I0,A,I0,A)') 'Positive diagonals ',pcount,' (out of ',n,')'
    CALL Info('CRS_FCTLowOrder',Message,Level=8)

    A % Values = A % Values + A % FCT_D
    
    ! Just some optional stuff for debugging purposes
    IF(.FALSE.) THEN
      CALL CRS_RowSumInfo( A, A % BulkValues ) 
      CALL CRS_RowSumInfo( A, A % Values ) 
      CALL CRS_RowSumInfo( A, A % FCT_D ) 
    END IF

    ! Create a lumped mass matrix by computing the rowsums of the 
    ! initial mass matrix.
    CALL Info('CRS_FCTLowOder','Creating lumped mass matrix',Level=10)
    IF( .NOT. ASSOCIATED( A % MassValuesLumped ) ) THEN         
      ALLOCATE( A % MassValuesLumped( n ) )
    END IF
    ML => A % MassValuesLumped     
    DO i=1,n
      msum = 0.0_dp
      DO j=Rows(i),Rows(i+1)-1
        msum = msum + A % MassValues(j)
      END DO
      ML(i) = msum 
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE CRS_FCTLowOrder
!------------------------------------------------------------------------------

 

!------------------------------------------------------------------------------
!> Copies the matrix topology from matrix A to build matrix B. Note that the 
!> topology is really reused, so that if matrix A is destroyed also matrix B 
!> becomes unusable.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_CopyMatrixTopology(A,B)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A  !< Initial matrix
    TYPE(Matrix_t) :: B              !< New matrix with same matrix topology
!------------------------------------------------------------------------------
    INTEGER :: kb,n,istat
    N = A % NumberOfRows

    IF ( n == 0 ) THEN
      CALL Fatal('CRS_CopyMatrixTopology','The first matrix is assumed to exist')
    END IF
    
    IF ( A % FORMAT /= MATRIX_CRS ) THEN
      CALL Fatal('CRS_CopyMatrixTopology','The matrix structure should be CRS!')
    END IF
    IF ( B % NumberOfRows /= 0 ) THEN
      CALL Fatal('CRS_CopyMatrixTopology','The other matrix is assumed not to exist')
    END IF

    CALL Info('CRS_CopyMatrixTopology','Reusing matrix topology',Level=9)


    B % NumberOfRows = n
    B % ListMatrix => NULL()
    B % FORMAT = A % FORMAT

    B % Rows => A % Rows
    B % Cols => A % Cols
    IF( ASSOCIATED( A % Diag ) ) THEN
      B % Diag => A % Diag
    END IF

!    IF( ASSOCIATED( A % rhs ) ) THEN
!      ALLOCATE( B % rhs(n), STAT = istat )
!      IF( istat /= 0 ) CALL Fatal('CRS_CopyMatrixTopology','memory allocation error 1')
!      B % rhs = 0.0_dp
!    END IF

    kb = SIZE( A % Values )
    ALLOCATE( B % Values(kb), STAT=istat )
    IF( istat /= 0 ) CALL Fatal('CRS_CopyMatrixTopology','memory allocation error 2')
    B % Values = 0.0_dp

!------------------------------------------------------------------------------
  END SUBROUTINE CRS_CopyMatrixTopology
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Pics a block from matrix A to build matrix B. It is assumed that the 
!> matrix is split into given number of equally sized blocks.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_BlockMatrixPick(A,B,Blocks,Nrow,Ncol)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A   !< Initial matrix
    TYPE(Matrix_t) :: B   !< Submatrix picked from the larger matrix
    INTEGER, INTENT(IN) :: Blocks     !< Number of blocks in the initial matrix
    INTEGER, INTENT(IN) :: Nrow       !< Row to be picked
    INTEGER, INTENT(IN) :: Ncol       !< Column to be picked
!------------------------------------------------------------------------------    
	INTEGER :: i,j,k,l,kb,n,Nrow0,Ncol0,nsub
    INTEGER :: lsub,isub,istat,modNcol
    LOGICAL :: NewMatrix, Diagonal

    IF(Blocks <= 1) THEN
      CALL Fatal('CRS_BlockMatrixPick','No applicable to just one block!')
      RETURN
    END IF

    N = A % NumberOfRows
    Nsub = N / Blocks
    modNcol = MOD( Ncol,Blocks)

    NewMatrix = ( B % NumberOfRows == 0 ) 
    Diagonal = ( Nrow == Ncol ) 

    IF( NewMatrix ) THEN
      B % ListMatrix => NULL()
      B % FORMAT = MATRIX_CRS

      B % NumberOfRows = Nsub    
      kb = 0
      
      DO isub=1,Nsub
        i = Blocks * ( isub - 1 ) + Nrow 
        DO k= A % Rows(i), A % Rows(i+1)-1
          l = A % Cols(k)
          IF( MOD(l,Blocks) == modNcol ) THEN
            kb = kb + 1
          END IF
        END DO
      END DO
      
      IF( kb == 0 ) THEN
        PRINT *,'Nrow:',Nrow,'Ncol:',Ncol
        CALL Warn('CRS_BlockMatrixPick','No matrix entries in submatrix')
        RETURN
      END IF

      ALLOCATE(B % Rows(nsub+1),B % Cols(kb), B % Values(kb),STAT=istat )
      IF( istat /= 0 ) CALL Fatal('CRS_BlockMatrixPick','memory allocation error 1')
    END IF

    IF( Diagonal ) THEN
      IF( .NOT. ASSOCIATED( B % Diag ) ) THEN
        ALLOCATE( B % Diag(nsub), STAT=istat)
        IF( istat /= 0 ) CALL Fatal('CRS_BlockMatrixPick','memory allocation error 2')      
      END IF
      IF( .NOT. ASSOCIATED( B % Rhs ) ) THEN
        ALLOCATE( B % rhs(nsub), STAT=istat)
        IF( istat /= 0 ) CALL Fatal('CRS_BlockMatrixPick','memory allocation error 3')      
      END IF
    END IF


    kb = 1
    DO isub=1,Nsub

      IF( NewMatrix ) B % Rows(isub) = kb
      i = Blocks * ( isub - 1 ) + Nrow 

      DO k = A % Rows(i), A % Rows(i+1)-1

        l = A % Cols(k)
        IF( MOD( l, Blocks ) == modNcol ) THEN
          lsub = ( l - 1) / Blocks + 1
          B % Values(kb) = A % Values(k)
          IF( NewMatrix ) THEN
            B % Cols(kb) = lsub
            IF( Diagonal .AND. isub == lsub ) B % Diag(isub) = kb
          END IF
          kb = kb + 1  
        END IF
      END DO
      
      IF( Diagonal ) B % rhs(isub) = A % rhs(i)

    END DO
    IF( NewMatrix ) B % Rows(Nsub+1) = kb

!------------------------------------------------------------------------------
  END SUBROUTINE CRS_BlockMatrixPick
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Pics a block from matrix A to build matrix B. 
!> This subroutine enables the use of 
!> nontrivial block decompositions. 
!------------------------------------------------------------------------------
  SUBROUTINE CRS_BlockMatrixPick2(A,B,BlockStruct,Nrow,Ncol)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A   !< Initial matrix
    TYPE(Matrix_t) :: B   !< Submatrix picked from the larger matrix
    INTEGER, POINTER :: BlockStruct(:)     !< Block decomposition structure of the initial matrix
    INTEGER, INTENT(IN) :: Nrow       !< Row to be picked
    INTEGER, INTENT(IN) :: Ncol       !< Column to be picked
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,kb,n,Nrow0,Ncol0,nsub,Mrow,Mcol,mr,mc,imsub,lmsub
    INTEGER :: lsub,isub,istat,modNcol,Blocks
    LOGICAL :: NewMatrix, Allocated, Diagonal, Hit
    INTEGER, ALLOCATABLE :: Irow(:), Icol(:)
    
    Blocks = SIZE( BlockStruct )

    IF(Blocks <= 1) THEN
      CALL Fatal('CRS_BlockMatrixPick','No applicable to just one block!')
      RETURN
    END IF

    N = A % NumberOfRows

    Mrow = 0
    Mcol = 0
    ALLOCATE( Irow(Blocks), Icol(Blocks) ) 
    Irow = 0
    Icol = 0

    DO i=1,Blocks
      IF( BlockStruct(i) == Nrow ) THEN
        Mrow = Mrow + 1
        Irow(Mrow) = i
      END IF
      IF( BlockStruct(i) == Ncol ) THEN
        Mcol = Mcol + 1
        Icol(Mcol) = i
      END IF
    END DO

    IF( Mrow == 0 .OR. Mcol == 0 ) THEN
      CALL Fatal('CRS_BlockMatrixPick','Nothing to pick!')
    END IF

    Nsub = N / Blocks
    modNcol = MOD( Ncol,Blocks)

    NewMatrix = ( B % NumberOfRows == 0 ) 
    Allocated = .NOT. NewMatrix
    Diagonal = ( Nrow == Ncol ) 

    IF( .NOT. Allocated ) THEN
      PRINT *,'block rows no:',Mrow,' inds:',Irow(1:Mrow)
      PRINT *,'block cols no:',Mcol,' inds:',Icol(1:Mcol)
      B % ListMatrix => NULL()
      B % FORMAT = MATRIX_CRS
      B % NumberOfRows = Mrow *  Nsub    
    END IF

100 kb = 1      
    DO isub=1,Nsub

      DO mr=1,Mrow
        imsub = Mrow*(isub-1)+mr
        IF( Allocated .AND. NewMatrix ) B % Rows( imsub ) = kb
        i = Blocks * ( isub - 1 ) + Irow(mr) 
        
        DO k= A % Rows(i), A % Rows(i+1)-1
          l = A % Cols(k)
          Hit = .FALSE.

          DO mc = 1, Mcol
            IF( MOD(l,Blocks) == MOD( Icol(mc), Blocks) ) THEN
              Hit = .TRUE.
              EXIT
            END IF
          END DO

          IF( Hit ) THEN
            IF( Allocated ) THEN
              lmsub = Mcol * ( ( l - 1) / Blocks ) + mc
              B % Values(kb) = A % Values(k)
              IF( NewMatrix ) THEN
                B % Cols(kb) = lmsub
                IF( Diagonal ) THEN
                  IF( imsub == lmsub ) B % Diag(imsub) = kb
                END IF
              END IF            
              IF( Diagonal ) B % rhs(imsub) = A % rhs(i)              
            END IF
            kb = kb + 1
          END IF

        END DO
      END DO
    END DO
    
    IF( .NOT. Allocated ) THEN
      IF( kb == 1 ) THEN
        CALL Warn('CRS_BlockMatrixPick','No matrix entries in submatrix')
        RETURN
      END IF

      ALLOCATE(B % Rows(Mrow*nsub+1),B % Cols(kb-1), B % Values(kb-1),STAT=istat )
      IF( istat /= 0 ) CALL Fatal('CRS_BlockMatrixPick','memory allocation error 1')
      
      B % Rows(Mrow*Nsub+1) = kb
      
      IF( Diagonal ) THEN
        ALLOCATE( B % Diag(Mrow*nsub), B % rhs(Mrow*nsub), STAT=istat)
        IF( istat /= 0 ) CALL Fatal('CRS_BlockMatrixPick','memory allocation error 2')      
      END IF

      IF( A % COMPLEX ) THEN
        IF( MOD( Mrow, 2) == 0 .AND. MOD( Mcol, 2) == 0 ) THEN
          B % COMPLEX = .TRUE.
        END IF
      END IF

      Allocated = .TRUE.
      GOTO 100
    END IF


!------------------------------------------------------------------------------
  END SUBROUTINE CRS_BlockMatrixPick2
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Copies some preconditioning structures from matrix A to B. 
!> The intent is to allow saving of memory and CPU time for cases
!> where similar preconditioning could be used for many components.
!------------------------------------------------------------------------------
  FUNCTION CRS_CopyMatrixPrec(A,B) RESULT(Status)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A  !< Initial matrix
    TYPE(Matrix_t) :: B  !< New matrix with the same preconditioner
    LOGICAL :: Status  !< Returns true if copying was successful or unnecessary
!------------------------------------------------------------------------------
    INTEGER :: n

    Status = .FALSE.

    IF( ASSOCIATED( B % IluValues ) ) THEN
!      CALL Info('CRS_CopyMatrixPrec','ILU preconditioner already exists')
      Status = .TRUE.
    END IF

    IF( ASSOCIATED( B % Ematrix ) ) THEN
!      CALL Info('CRS_CopyMatrixPrec','AMG preconditioner already exists')
      Status = .TRUE.      
    END IF
    
    IF( Status ) RETURN

    IF( SIZE( A % Values ) /= SIZE( B % Values ) ) THEN
      PRINT *,'sizes',SIZE( A % Values ), SIZE( B % Values )
      CALL Info('CRS_CopyMatrixPrec','Mismatch in size, returning')            
      RETURN
    END IF
    
    IF( ASSOCIATED( A % IluValues ) ) THEN
      CALL Info('CRS_CopyMatrixPrec','Reusing ILU preconditioner topology',Level=9)      
      B % IluRows => A % IluRows
      B % IluCols => A % IluCols
      B % IluDiag => A % IluDiag

      n = SIZE( A % ILUValues ) 
      ALLOCATE( B % IluValues(n) )
      B % IluValues = 0.0_dp
      Status = .TRUE.
      RETURN
    END IF

    RETURN

    ! This should be still worked on....
    !------------------------------------------------
    IF( ASSOCIATED( A % Ematrix ) ) THEN
      CALL Info('CRS_CopyMatrixPrec','Reusing AMG preconditioner topology',Level=9)      
      B % Ematrix => A % Ematrix
    END IF

  END FUNCTION CRS_CopyMatrixPrec



!------------------------------------------------------------------------------
!>    Builds an incomplete (ILU(n)) factorization for a iterative solver
!>    preconditioner. Real matrix version.
!------------------------------------------------------------------------------
  FUNCTION CRS_IncompleteLU(A,ILUn) RESULT(Status)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A          !< Structure holding input matrix, will also hold the factorization on exit.
    INTEGER, INTENT(IN) :: ILUn  !< Order of fills allowed 0-9
    LOGICAL :: Status            !< Whether or not the factorization succeeded.
!------------------------------------------------------------------------------
    LOGICAL :: Warned
    INTEGER :: i,j,k,l,m,n,istat

    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    REAL(KIND=dp), POINTER :: ILUValues(:), Values(:)

#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: st, tx
#else
    REAL(KIND=dp) :: CPUTime, st, tx
#endif

    LOGICAL, ALLOCATABLE :: C(:), D(:)
    REAL(KIND=dp), ALLOCATABLE ::  S(:), T(:)

    INTEGER, POINTER :: ILUCols(:),ILURows(:),ILUDiag(:)

    TYPE(Matrix_t), POINTER :: A1
!------------------------------------------------------------------------------
    WRITE(Message,'(a,i1,a)')  &
         'ILU(',ILUn,') (Real), Starting Factorization:'
    CALL Info( 'CRS_IncompleteLU', Message, Level = 5 )
    st = CPUTime()

    N = A % NumberOfRows

    IF(N == 0) THEN
       A % ILURows => A % Rows
       A % ILUCols => A % Cols
       A % ILUDiag => A % Diag

       Status = .TRUE.
       RETURN
    END IF

    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    IF(ASSOCIATED(A % PrecValues)) THEN
      Values => A % PrecValues
    ELSE
      Values => A % Values
    END IF

    IF ( .NOT. ASSOCIATED(A % ILUValues) ) THEN

       IF ( ILUn == 0 ) THEN
          A % ILURows => A % Rows
          A % ILUCols => A % Cols
          A % ILUDiag => A % Diag
       ELSE
          CALL InitializeILU1( A,n )

          IF ( ILUn > 1 ) THEN
             ALLOCATE( A1 )

             DO i=1,ILUn-1
                CALL Info('CRS_IncompleLU','Recursive round: '//I2S(i))

                A1 % Cols => A % ILUCols
                A1 % Rows => A % ILURows
                A1 % Diag => A % ILUDiag

                CALL InitializeILU1( A1,n )

                A % ILUCols => A1 % ILUCols
                A % ILURows => A1 % ILURows
                A % ILUDiag => A1 % ILUDiag

                DEALLOCATE( A1 % Cols, A1 % Rows, A1 % Diag )
             END DO

             DEALLOCATE(A1)
          END IF
       END IF

       m = A % ILURows(N+1)-1
       ALLOCATE( A % ILUValues(A % ILURows(N+1)-1), STAT=istat )

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_IncompleteLU', 'Memory allocation error.' )
       ELSE
         CALL Info('CRS_IncompleteLU','Allocated LU matrix of size: '//i2s(m),Level=10 ) 
       END IF
    END IF

    ILURows   => A % ILURows
    ILUCols   => A % ILUCols
    ILUDiag   => A % ILUDiag
    ILUValues => A % ILUValues
!
!   Allocate space for storing one full row:
!   ----------------------------------------
    ALLOCATE( C(n), S(n) )
    C = .FALSE.
    S =  0.0d0

    IF ( A % Cholesky ) THEN
      CALL Info('CRS_IncompleteLU','Performing incomplete Cholesky',Level=12)

      ALLOCATE( T(n) )
      T =  0._dp
     !
     ! The factorization row by row:
     ! -----------------------------
     Warned = .FALSE.
     DO i=1,N

       ! Convert current row to full form for speed,
       ! only flagging the nonzero entries:
       ! -------------------------------------------
       DO k=Rows(i), Diag(i)
         j = Cols(k)
         T(j) = Values(k)
       END DO

       DO k=ILURows(i), ILUDiag(i)
         j = ILUCols(k)
         C(j) = .TRUE.
         S(j) = ILUValues(k)
       END DO

       ! This is the factorization part for the current row:
       ! ---------------------------------------------------
       S(i) = T(i)
       DO m=ILURows(i),ILUDiag(i)-1
         j = ILUCols(m)
         S(j) = T(j)
         DO l = ILURows(j),ILUDiag(j)-1
           k = ILUCols(l)
           S(j) = S(j) - S(k) * ILUValues(l)
         END DO
         S(j) = S(j) * ILUValues(ILUDiag(j))
         S(i) = S(i) - S(j)**2
       END DO

       IF ( S(i) <= AEPS ) THEN
         S(i) = 1._dp
         IF ( .NOT. Warned )  THEN
           CALL Warn( 'Cholesky factorization:', &
               'Negative diagonal: not pos.def. or badly conditioned matrix' )
           Warned = .TRUE.
         END IF
       ELSE
         S(i) = 1._dp / SQRT(S(i))
       END IF

       ! Convert the row back to  CRS format:
       ! ------------------------------------
       DO k=Rows(i), Diag(i)
         j = Cols(k) 
         T(j) = 0._dp
       END DO

       DO k=ILURows(i), ILUDiag(i)
         j = ILUCols(k)
         ILUValues(k) = S(j)
         S(j) =  0._dp
         C(j) = .FALSE.
       END DO
     END DO

    ELSE
      CALL Info('CRS_IncompleteLU','Performing incomplete LU',Level=12)

     ! The factorization row by row:
     ! -----------------------------
     DO i=1,N

       ! Convert current row to full form for speed,
       ! only flagging the nonzero entries:
       ! -------------------------------------------
       DO k=Rows(i), Rows(i+1)-1
          S(Cols(k)) = Values(k)
       END DO

       DO k = ILURows(i), ILURows(i+1)-1
          C(ILUCols(k)) = .TRUE.
       END DO
!
!      This is the factorization part for the current row:
!      ---------------------------------------------------
       DO m=ILURows(i),ILUDiag(i)-1
         k = ILUCols(m)
         IF ( S(k) == 0._dp ) CYCLE

         IF ( ABS(ILUValues(ILUDiag(k))) > AEPS ) &
           S(k) = S(k) / ILUValues(ILUDiag(k)) 

         DO l = ILUDiag(k)+1, ILURows(k+1)-1
           j = ILUCols(l)
           IF ( C(j) ) THEN
             S(j) = S(j) - S(k) * ILUValues(l)
           END IF
         END DO
       END DO

!
!      Convert the row back to  CRS format:
!      ------------------------------------
       DO k=ILURows(i), ILURows(i+1)-1
         IF ( C(ILUCols(k)) ) THEN
           ILUValues(k)  = S(ILUCols(k))
           S(ILUCols(k)) =  0.0d0
           C(ILUCols(k)) = .FALSE.
         END IF
       END DO
     END DO

     !
     ! Prescale the diagonal for the LU solve:
     ! ---------------------------------------
     DO i=1,N
       IF ( ABS(ILUValues(ILUDiag(i))) < AEPS ) THEN
         ILUValues(ILUDiag(i)) = 1.0d0
       ELSE
         ILUValues(ILUDiag(i)) = 1.0d0 / ILUValues(ILUDiag(i))
       END IF
     END DO
    END IF


!------------------------------------------------------------------------------
    WRITE(Message,'(a,i1,a,i9)') 'ILU(', ILUn, &
        ') (Real), NOF nonzeros: ',ILURows(n+1)
    CALL Info( 'CRS_IncompleteLU', Message, Level=5 )

    WRITE(Message,'(a,i1,a,i9)') 'ILU(', ILUn, &
        ') (Real), filling (%) : ',   &
         FLOOR(ILURows(n+1)*(100.0d0/Rows(n+1)))
    CALL Info( 'CRS_IncompleteLU', Message, Level=5 )

    WRITE(Message,'(A,I1,A,F8.2)') 'ILU(',ILUn, &
        ') (Real), Factorization ready at (s): ', CPUTime()-st
    CALL Info( 'CRS_IncompleteLU', Message, Level=5 )

    Status = .TRUE.
!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE InitializeILU1( A, n )
!------------------------------------------------------------------------------
      TYPE(Matrix_t) :: A
      INTEGER :: n
      INTEGER :: i,j,k,l,istat,RowMin,RowMax,Nonzeros
      INTEGER :: C(n)
      INTEGER, POINTER :: Cols(:),Rows(:),Diag(:), &
           ILUCols(:),ILURows(:),ILUDiag(:)
!------------------------------------------------------------------------------

      Diag => A % Diag
      Rows => A % Rows
      Cols => A % Cols

      ALLOCATE( A % ILURows(N+1),A % ILUDiag(N),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_IncompleteLU', 'Memory allocation error.' )
      END IF

      ILURows => A % ILURows
      ILUDiag => A % ILUDiag
!
!     Count fills, row by row:
!     ------------------------
      NonZeros = Rows(N+1) - 1

      C = 0
      DO i=1,n
         DO k=Rows(i), Rows(i+1)-1
            C(Cols(k)) = 1
         END DO

         DO k = Cols(Rows(i)), i-1
            IF ( C(k) /= 0 ) THEN
               DO l=Diag(k)+1, Rows(k+1)-1
                  j = Cols(l)
                  IF ( C(j) == 0 ) THEN
                    Nonzeros = Nonzeros + 1

                    IF( Nonzeros == HUGE( NonZeros ) ) THEN
                      CALL Error('CRS_IncompleteLU','Number of nonzeros larger than HUGE(Integer)')
                      CALL Fatal('CRS_IncompleteLU','Try some cheaper preconditioner!')                     
                    END IF

                  END IF
               END DO
            END IF
         END DO

         DO k = Rows(i), Rows(i+1)-1
            C(Cols(k)) = 0
         END DO
      END DO

      CALL Info('CRS_IncompleteLU','Number of nonzeros: '//TRIM(I2S(NonZeros)),Level=12)

!------------------------------------------------------------------------------

      ALLOCATE( A % ILUCols(Nonzeros),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_IncompleteLU', 'Memory allocation error.' )
      END IF
      ILUCols => A % ILUCols

!------------------------------------------------------------------------------

!
!     Update row nonzero structures: 
!     ------------------------------
      C = 0
      ILURows(1) = 1
      DO i=1,n
         DO k=Rows(i), Rows(i+1)-1
            C(Cols(k)) = 1
         END DO

         RowMin = Cols(Rows(i))
         RowMax = Cols(Rows(i+1)-1)

         DO k=RowMin, i-1
            IF ( C(k) == 1 ) THEN
               DO l=Diag(k)+1,Rows(k+1)-1
                  j = Cols(l)
                  IF ( C(j) == 0 ) THEN
                     C(j) = 2
                     RowMax = MAX( RowMax, j )
                  END IF
               END DO
            END IF
         END DO

         j = ILURows(i) - 1
         DO k = RowMin, RowMax 
            IF ( C(k) > 0 ) THEN
               j = j + 1
               C(k) = 0
               ILUCols(j) = k
               IF ( k == i ) ILUDiag(i) = j
            END IF
         END DO
         ILURows(i+1) = j + 1
      END DO

      CALL Info('CRS_IncompleteLU','Updated nonzero elements',Level=20)

!------------------------------------------------------------------------------
    END SUBROUTINE InitializeILU1
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  END FUNCTION CRS_IncompleteLU
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Buids an incomplete (ILU(n)) factorization for an iterative solver
!>    preconditioner. Complex matrix version.
!------------------------------------------------------------------------------
  FUNCTION CRS_ComplexIncompleteLU(A,ILUn) RESULT(Status)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Strcture holding input matrix, will also hold the factorization on exit.
    INTEGER, INTENT(IN) :: ILUn   !< Order of fills allowed 0-9
    LOGICAL :: Status  !< Whether or not the factorization succeeded.
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,istat

    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    REAL(KIND=dp), POINTER ::  Values(:)
    COMPLEX(KIND=dp), POINTER :: ILUValues(:)

    INTEGER, POINTER :: ILUCols(:),ILURows(:),ILUDiag(:)

    TYPE(Matrix_t), POINTER :: A1

#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: st
#else
    REAL(KIND=dp) :: st, CPUTime
#endif

    LOGICAL, ALLOCATABLE :: C(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: S(:), T(:)
!------------------------------------------------------------------------------

    WRITE(Message,'(a,i1,a)') 'ILU(',ILUn,') (Complex), Starting Factorization:'
    CALL Info( 'CRS_ComplexIncompleteLU', Message, Level=5 )
    st = CPUTime()

    N = A % NumberOfRows
    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    IF(ASSOCIATED(A % PrecValues)) THEN
      Values => A % PrecValues
    ELSE
      Values => A % Values
    END IF

    IF ( .NOT.ASSOCIATED(A % CILUValues) ) THEN

       ALLOCATE( A1 )
       A1 % NumberOfRows = N / 2

       ALLOCATE( A1 % Rows(n/2+1) )
       ALLOCATE( A1 % Diag(n/2) )
       ALLOCATE( A1 % Cols(SIZE(A % Cols) / 4) )

       A1 % Rows(1) = 1
       k = 0
       DO i=1,n,2
          DO j=A % Rows(i),A % Rows(i+1)-1,2
             k = k + 1
             A1 % Cols(k) = (A % Cols(j)+1) / 2
             IF ( A % Cols(j) == i ) A1 % Diag((i+1)/2) = k
          END DO
          A1 % Rows((i+1)/2+1) = k+1
       END DO

       IF ( ILUn == 0 ) THEN
          A % ILUCols => A1 % Cols
          A % ILURows => A1 % Rows
          A % ILUDiag => A1 % Diag
       ELSE
          CALL InitializeComplexILU1( A1, n/2 )

          A % ILUCols => A1 % ILUCols
          A % ILURows => A1 % ILURows
          A % ILUDiag => A1 % ILUDiag

          DEALLOCATE( A1 % Cols,A1 % Rows,A1 % Diag )
       END IF

       DEALLOCATE( A1 )

       IF ( ILUn > 1 ) THEN
          ALLOCATE( A1 )
          A1 % NumberOfRows = N / 2

          DO i=1,ILUn-1
             A1 % Cols => A % ILUCols
             A1 % Rows => A % ILURows
             A1 % Diag => A % ILUDiag

             CALL InitializeComplexILU1( A1, n/2 )

             A % ILUCols => A1 % ILUCols
             A % ILURows => A1 % ILURows
             A % ILUDiag => A1 % ILUDiag

             DEALLOCATE( A1 % Cols,A1 % Rows,A1 % Diag )
          END DO

          DEALLOCATE(A1)
       END IF

       ALLOCATE( A % CILUValues(A % ILURows(N/2+1)),STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'CRS_ComplexIncompleteLU', 'Memory allocation error.' )
       END IF
    END IF

    ILURows   => A % ILURows
    ILUCols   => A % ILUCols
    ILUDiag   => A % ILUDiag
    ILUValues => A % CILUValues

!
!   Allocate space for storing one full row:
!   ----------------------------------------
    ALLOCATE( C(n/2), S(n/2) )
    C = .FALSE.
    S =  0.0d0

    IF ( A % Cholesky ) THEN
      ALLOCATE( T(n/2) )
      T =  0.0d0
     !
     ! The factorization row by row:
     ! -----------------------------
     DO i=1,N/2

       ! Convert current row to full form for speed,
       ! only flagging the nonzero entries:
       ! -------------------------------------------
       DO k = Rows(2*i-1), Rows(2*i)-1,2
          T((Cols(k)+1)/2) = CMPLX( Values(k), -Values(k+1), KIND=dp )
       END DO

       DO j=ILURows(i), ILUDiag(i)
          C(ILUCols(j)) = .TRUE.
          S(ILUCols(j)) = ILUValues(j)
       END DO

       ! This is the factorization part for the current row:
       ! ---------------------------------------------------
       S(i) = T(i)
       DO m=ILURows(i),ILUDiag(i)-1
         j = ILUCols(m)
         S(j) = T(j)
         DO l = ILURows(j),ILUDiag(j)-1
           k = ILUCols(l)
           S(j) = S(j) - S(k) * CONJG(ILUValues(l))
         END DO
         S(j) = S(j) * ILUValues(ILUDiag(j))
         S(i) = S(i) - S(j)*CONJG(S(j))
       END DO

       S(i) = 1._dp / SQRT(S(i))

       ! Convert the row back to  CRS format:
       ! ------------------------------------
       DO k = Rows(2*i-1), Rows(2*i)-1,2
         T((Cols(k)+1)/2) =  0._dp
       END DO

       DO k=ILURows(i), ILUDiag(i)
         ILUValues(k)  = S(ILUCols(k))
         S(ILUCols(k)) =  0._dp
         C(ILUCols(k)) = .FALSE.
       END DO
     END DO

    ELSE

     ! The factorization row by row:
     ! -----------------------------
     DO i=1,N/2

       ! Convert the current row to full form for speed,
       ! only flagging the nonzero entries:
       ! -----------------------------------------------
       DO k = ILURows(i), ILURows(i+1)-1
         C(ILUCols(k)) = .TRUE.
       END DO

       DO k = Rows(2*i-1), Rows(2*i)-1,2
         S((Cols(k)+1)/2) = CMPLX( Values(k), -Values(k+1), KIND=dp )
       END DO

       ! This is the factorization part for the current row:
       ! ---------------------------------------------------
       DO m=ILURows(i),ILUDiag(i)-1
         k = ILUCols(m)
         IF ( S(k) == 0._dp ) CYCLE

         IF ( ABS(ILUValues(ILUDiag(k))) > AEPS ) &
           S(k) = S(k) / ILUValues(ILUDiag(k)) 

         DO l = ILUDiag(k)+1, ILURows(k+1)-1
           j = ILUCols(l)
           IF ( C(j) ) THEN
             S(j) = S(j) - S(k) * ILUValues(l)
           END IF
         END DO
       END DO

       ! Convert the row back to  CRS format:
       ! ------------------------------------
       DO k=ILURows(i), ILURows(i+1)-1
         ILUValues(k)  = S(ILUCols(k))
         S(ILUCols(k)) =  0._dp
         C(ILUCols(k)) = .FALSE.
       END DO
     END DO

     ! Prescale the diagonal for the LU solve:
     ! ---------------------------------------
     DO i=1,n/2
       IF ( ABS(ILUValues(ILUDiag(i))) < AEPS ) THEN
         ILUValues(ILUDiag(i)) = 1._dp
       ELSE
         ILUValues(ILUDiag(i)) = 1._dp / ILUValues(ILUDiag(i))
       END IF
     END DO
    END IF

!------------------------------------------------------------------------------

    WRITE(Message,'(a,i1,a,i9)') 'ILU(', ILUn, &
        ') (Complex), NOF nonzeros: ',ILURows(n/2+1)
    CALL Info( 'CRS_ComplexIncompleteLU', Message, Level=5 )

    WRITE(Message,'(a,i1,a,i9)') 'ILU(', ILUn, &
        ') (Complex), filling (%) : ',   &
         FLOOR(ILURows(n/2+1)*(400.0d0/Rows(n+1)))
    CALL Info( 'CRS_ComplexIncompleteLU', Message, Level=5 )

    WRITE(Message,'(A,I1,A,F8.2)') 'ILU(',ILUn, &
        ') (Complex), Factorization ready at (s): ', CPUTime()-st
    CALL Info( 'CRS_ComplexIncompleteLU', Message, Level=5 )

    Status = .TRUE.
!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE InitializeComplexILU1( A, n )
!------------------------------------------------------------------------------
      TYPE(Matrix_t) :: A
      INTEGER :: n

      INTEGER :: i,j,k,l,istat,RowMin,RowMax,Nonzeros

      INTEGER :: C(n)
      INTEGER, POINTER :: Cols(:),Rows(:),Diag(:), &
           ILUCols(:),ILURows(:),ILUDiag(:)
!------------------------------------------------------------------------------

      Diag => A % Diag
      Rows => A % Rows
      Cols => A % Cols

      ALLOCATE( A % ILURows(N+1),A % ILUDiag(N),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ComplexIncompleteLU', 'Memory allocation error.' )
      END IF

      ILURows => A % ILURows
      ILUDiag => A % ILUDiag

!
!     Count fills, row by row:
!     ------------------------
      NonZeros = Rows(N+1) - 1
      C = 0
      DO i=1,n
         DO k=Rows(i), Rows(i+1)-1
            C(Cols(k)) = 1
         END DO

         DO k = Cols(Rows(i)), i-1
            IF ( C(k) /= 0 ) THEN
               DO l=Diag(k)+1, Rows(k+1)-1
                  j = Cols(l)
                  IF ( C(j) == 0 ) Nonzeros = Nonzeros + 1
               END DO
            END IF
         END DO

         DO k = Rows(i), Rows(i+1)-1
            C(Cols(k)) = 0
         END DO
      END DO

!------------------------------------------------------------------------------

      ALLOCATE( A % ILUCols(Nonzeros),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ComplexIncompleteLU', 'Memory allocation error.' )
      END IF
      ILUCols => A % ILUCols

!------------------------------------------------------------------------------

!
!     Update row nonzero structures: 
!     ------------------------------
      C = 0
      ILURows(1) = 1
      DO i=1,n
         DO k=Rows(i), Rows(i+1)-1
            C(Cols(k)) = 1
         END DO

         RowMin = Cols(Rows(i))
         RowMax = Cols(Rows(i+1)-1)

         DO k=RowMin, i-1
            IF ( C(k) == 1 ) THEN
               DO l=Diag(k)+1,Rows(k+1)-1
                  j = Cols(l)
                  IF ( C(j) == 0 ) THEN
                     C(j) = 2
                     RowMax = MAX( RowMax, j )
                  END IF
               END DO
            END IF
         END DO

         j = ILURows(i) - 1
         DO k = RowMin, RowMax 
            IF ( C(k) > 0 ) THEN
               j = j + 1
               C(k) = 0
               ILUCols(j) = k
               IF ( k == i ) ILUDiag(i) = j
            END IF
         END DO
         ILURows(i+1) = j + 1
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE InitializeComplexILU1
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  END FUNCTION CRS_ComplexIncompleteLU
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Buids an incomplete (ILUT) factorization for an iterative solver
!>    preconditioner. Real matrix version.
!------------------------------------------------------------------------------
  FUNCTION CRS_ILUT(A,TOL) RESULT(Status)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding input matrix, will also hold the factorization on exit.
    REAL(KIND=dp), INTENT(IN) :: TOL  !< Drop toleranece: if ILUT(i,j) <= NORM(A(i,:))*TOL the value is dropped.
    LOGICAL :: Status    !< Whether or not the factorization succeeded.
!------------------------------------------------------------------------------
    INTEGER :: n
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: t
#else
    REAL(KIND=dp) :: CPUTime, t
#endif
!------------------------------------------------------------------------------

    CALL Info( 'CRS_ILUT', 'Starting factorization:', Level=5 )
    t = CPUTime()

    n = A % NumberOfRows

    IF ( ASSOCIATED( A % ILUValues ) ) THEN
       DEALLOCATE( A % ILURows, A % ILUDiag, A % ILUCols, A % ILUValues )
    END IF
!
!   ... and then to the point:
!   --------------------------
    CALL ComputeILUT( A, n, TOL )
! 
    WRITE( Message, * ) 'ILU(T) (Real), NOF nonzeros: ',A % ILURows(N+1)
    CALL Info( 'CRS_ILUT', Message, Level=5 )
    WRITE( Message, * ) 'ILU(T) (Real), filling (%): ', &
         FLOOR(A % ILURows(N+1)*(100.0d0/A % Rows(N+1)))
    CALL Info( 'CRS_ILUT', Message, Level=5 )
    WRITE(Message,'(A,F8.2)') 'ILU(T) (Real), Factorization ready at (s): ', CPUTime()-t
    CALL Info( 'CRS_ILUT', Message, Level=5 )

    Status = .TRUE.
!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE ComputeILUT( A,n,TOL )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: TOL
      INTEGER :: n
      TYPE(Matrix_t) :: A
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: WORKN = 128

      INTEGER :: i,j,k,l,istat, RowMin, RowMax
      REAL(KIND=dp) :: NORMA

      REAL(KIND=dp), POINTER :: Values(:), ILUValues(:), CWork(:)

      INTEGER, POINTER :: Cols(:), Rows(:), Diag(:), &
           ILUCols(:), ILURows(:), ILUDiag(:), IWork(:)

      LOGICAL :: C(n)
#ifdef USE_ISO_C_BINDINGS
      REAL(KIND=dp) :: S(n), cptime, ttime, t
#else
      REAL(KIND=dp) :: S(n), CPUTime, cptime, ttime, t
#endif
!------------------------------------------------------------------------------

      ttime  = CPUTime()
      cptime = 0.0d0

      Diag => A % Diag
      Rows => A % Rows
      Cols => A % Cols
      IF(ASSOCIATED(A % PrecValues)) THEN
        Values => A % PrecValues
      ELSE
        Values => A % Values
      END IF

      ALLOCATE( A % ILURows(N+1),A % ILUDiag(N),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ILUT', 'Memory allocation error.' )
      END IF

      ILURows => A % ILURows
      ILUDiag => A % ILUDiag

      ALLOCATE( ILUCols( WORKN*N ),  ILUValues( WORKN*N ), STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ILUT', 'Memory allocation error.' )
      END IF
!
!     The factorization row by row:
!     -----------------------------
      ILURows(1) = 1
      S =  0.0d0
      C = .FALSE.

      DO i=1,n
!
!        Convert the current row to full form for speed,
!        only flagging the nonzero entries:
!        -----------------------------------------------
         DO k=Rows(i), Rows(i+1) - 1
            C(Cols(k)) = .TRUE.
            S(Cols(k)) = Values(k)
         END DO
!
!        Check bandwidth for speed, bandwidth optimization
!        helps here ALOT, use it!
!        -------------------------------------------------
         RowMin = Cols(Rows(i))
         RowMax = Cols(Rows(i+1)-1)
!
!        Here is the factorization part for the current row:
!        ---------------------------------------------------
         DO k=RowMin,i-1
            IF ( C(k) ) THEN
               IF ( ABS(ILUValues(ILUDiag(k))) > AEPS ) &
                 S(k) = S(k) / ILUValues(ILUDiag(k)) 
              
               DO l=ILUDiag(k)+1, ILURows(k+1)-1
                  j = ILUCols(l)
                  IF ( .NOT. C(j) ) THEN
                     C(j) = .TRUE.
                     RowMax = MAX( RowMax,j )
                  END IF
                  S(j) = S(j) - S(k) * ILUValues(l)
               END DO
            END IF
         END DO
!
!        This is the ILUT part, drop element ILU(i,j), if
!        ABS(ILU(i,j)) <= NORM(A(i,:))*TOL:
!        -------------------------------------------------
         NORMA = SQRT( SUM( ABS(Values(Rows(i):Rows(i+1)-1))**2 ) )

         j = ILURows(i)-1
         DO k=RowMin, RowMax
            IF ( C(k) ) THEN
               IF ( ABS(S(k)) >= TOL*NORMA .OR. k==i ) THEN
                  j = j + 1
                  ILUCols(j)   = k
                  ILUValues(j) = S(k)
                  IF ( k == i ) ILUDiag(i) = j
               END IF
               S(k) =  0.0d0
               C(k) = .FALSE.
            END IF
         END DO
         ILURows(i+1) = j + 1
!
!        Preparations for the next row:
!        ------------------------------
         IF ( i < N ) THEN
!
!           Check if still enough workspace:
!           --------------------------------
            IF ( SIZE(ILUCols) < ILURows(i+1) + N ) THEN

               t = CPUTime()
!              k = ILURows(i+1) + MIN( WORKN, n-i ) * n
               k = ILURows(i+1) + MIN( 0.75d0*ILURows(i+1), (n-i)*(1.0d0*n) )
               ALLOCATE( IWork(k), STAT=istat )
               IF ( istat /= 0 ) THEN
                  CALL Fatal( 'CRS_ILUT', 'Memory allocation error.' )
               END IF
               DO j=1,ILURows(i+1)-1
                 IWork(j) = ILUCols(j)
               END DO
               DEALLOCATE( ILUCols )

               ALLOCATE( CWork(k), STAT=istat )
               IF ( istat /= 0 ) THEN
                  CALL Fatal( 'CRS_ILUT', 'Memory allocation error.' )
               END IF
               DO j=1,ILURows(i+1)-1
                 CWork(j) = ILUValues(j)
               END DO
               DEALLOCATE( ILUValues )

               ILUCols   => IWork
               ILUValues => CWork
               NULLIFY( IWork, CWork )
               cptime = cptime + ( CPUTime() - t )
!              PRINT*,'tot: ', CPUTime()-ttime, 'copy: ', cptime
            END IF
         END IF
      END DO
!
!     Prescale the diagonal for the LU solve:
!     ---------------------------------------
      DO i=1,n
         IF ( ABS(ILUValues(ILUDiag(i))) < AEPS ) THEN
            ILUValues(ILUDiag(i)) = 1.0d0
         ELSE
            ILUValues(ILUDiag(i)) = 1.0d0 / ILUValues(ILUDiag(i))
         END IF
      END DO

      A % ILUCols   => ILUCols
      A % ILUValues => ILUValues
!------------------------------------------------------------------------------
    END SUBROUTINE ComputeILUT
!------------------------------------------------------------------------------
  END FUNCTION CRS_ILUT
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Buids an incomplete (ILUT) factorization for an iterative solver
!>    preconditioner. Complex matrix version.
!------------------------------------------------------------------------------
  FUNCTION CRS_ComplexILUT(A,TOL) RESULT(Status)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding input matrix, will also hold the factorization on exit.
    REAL(KIND=dp), INTENT(IN) :: TOL  !< Drop toleranece: if ILUT(i,j) <= NORM(A(i,:))*TOL the value is dropped.
    LOGICAL :: Status    !< Whether or not the factorization succeeded.
!------------------------------------------------------------------------------
    INTEGER :: n
#ifdef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: t
#else
    REAL(KIND=dp) :: CPUTime, t
#endif
!------------------------------------------------------------------------------

    CALL Info( 'CRS_ComplexILUT', 'ILU(T) (Complex), Starting factorization: ', Level=5 )
    t = CPUTime()

    n = A % NumberOfRows / 2

    IF ( ASSOCIATED( A % CILUValues ) ) THEN
       DEALLOCATE( A % ILURows, A % ILUCols, A % ILUDiag, A % CILUValues )
    END IF
!
!   ... and then to the point:
!   --------------------------
    CALL ComplexComputeILUT( A, n, TOL )
 
!------------------------------------------------------------------------------
    
    WRITE( Message, * ) 'ILU(T) (Complex), NOF nonzeros: ',A % ILURows(n+1)
    CALL Info( 'CRS_ComplexILUT', Message, Level=5 )
    WRITE( Message, * ) 'ILU(T) (Complex), filling (%): ', &
         FLOOR(A % ILURows(n+1)*(400.0d0/A % Rows(2*n+1)))
    CALL Info( 'CRS_ComplexILUT', Message, Level=5 )
    WRITE(Message,'(A,F8.2)') 'ILU(T) (Complex), Factorization ready at (s): ', CPUTime()-t
    CALL Info( 'CRS_ComplexILUT', Message, Level=5 )

    Status = .TRUE.
!------------------------------------------------------------------------------

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE ComplexComputeILUT( A,n,TOL )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: TOL
      INTEGER :: n
      TYPE(Matrix_t) :: A
!------------------------------------------------------------------------------
      INTEGER, PARAMETER :: WORKN = 128

      INTEGER :: i,j,k,l,istat,RowMin,RowMax
      REAL(KIND=dp) :: NORMA

      REAL(KIND=dp), POINTER :: Values(:)
      COMPLEX(KIND=dp), POINTER :: ILUValues(:), CWork(:)

      INTEGER, POINTER :: Cols(:), Rows(:), Diag(:), &
           ILUCols(:), ILURows(:), ILUDiag(:), IWork(:)

      LOGICAL :: C(n)
      COMPLEX(KIND=dp) :: S(n)
!------------------------------------------------------------------------------

      Diag => A % Diag
      Rows => A % Rows
      Cols => A % Cols
      IF (ASSOCIATED(A % PrecValues)) THEN
        Values => A % PrecValues
      ELSE
        Values => A % Values
      END IF

      ALLOCATE( A % ILURows(n+1),A % ILUDiag(n),STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ComplexILUT', 'Memory allocation error.' )
      END IF

      ILURows => A % ILURows
      ILUDiag => A % ILUDiag

      ALLOCATE( ILUCols( WORKN*N ),  ILUValues( WORKN*N ), STAT=istat )
      IF ( istat /= 0 ) THEN
         CALL Fatal( 'CRS_ComplexILUT', 'Memory allocation error.' )
      END IF
!
!     The factorization row by row:
!     -----------------------------
      ILURows(1) = 1
      C = .FALSE.
      S = CMPLX( 0.0d0, 0.0d0, KIND=dp )

      DO i=1,n
!
!        Convert the current row to full form for speed,
!        only flagging the nonzero entries:
!        -----------------------------------------------
         DO k=Rows(2*i-1), Rows(2*i)-1,2
            C((Cols(k)+1) / 2) = .TRUE.
            S((Cols(k)+1) / 2) = CMPLX( Values(k), -Values(k+1), KIND=dp )
         END DO
!
!        Check bandwidth for speed, bandwidth optimization
!        helps here ALOT, use it!
!        -------------------------------------------------
         RowMin = (Cols(Rows(2*i-1)) + 1) / 2
         RowMax = (Cols(Rows(2*i)-1) + 1) / 2
!
!        Here is the factorization part for the current row:
!        ---------------------------------------------------
         DO k=RowMin,i-1
            IF ( C(k) ) THEN
               IF ( ABS(ILUValues(ILUDiag(k))) > AEPS ) &
                 S(k) = S(k) / ILUValues(ILUDiag(k)) 
              
               DO l=ILUDiag(k)+1, ILURows(k+1)-1
                  j = ILUCols(l)
                  IF ( .NOT. C(j) ) THEN
                     C(j) = .TRUE.
                     RowMax = MAX( RowMax,j )
                  END IF
                  S(j) = S(j) - S(k) * ILUValues(l)
               END DO
            END IF
         END DO

!
!        This is the ILUT part, drop element ILUT(i,j), if
!        ABS(ILUT(i,j)) <= NORM(A(i,:))*TOL:
!        -------------------------------------------------
         NORMA = 0.0d0
         DO k = Rows(2*i-1), Rows(2*i)-1, 2
            NORMA = NORMA + Values(k)**2 + Values(k+1)**2
         END DO
         NORMA = SQRT(NORMA)

         j = ILURows(i)-1
         DO k=RowMin, RowMax
            IF ( C(k) ) THEN
               IF ( ABS(S(k)) > TOL*NORMA .OR. k==i ) THEN
                  j = j + 1
                  ILUCols(j)   = k
                  ILUValues(j) = S(k)
                  IF ( k == i ) ILUDiag(i) = j
               END IF
               C(k) = .FALSE.
               S(k) = CMPLX( 0.0d0, 0.0d0, KIND=dp )
            END IF
         END DO
         ILURows(i+1) = j + 1
!
!        Preparations for the next row:
!        ------------------------------
         IF ( i < N ) THEN
!
!           Check if still enough workspace:
!           --------------------------------
            IF ( SIZE(ILUCols) < ILURows(i+1) + n ) THEN
!              k = ILURows(i+1) + MIN( WORKN, n-i ) * n
               k = ILURows(i+1) + MIN( 0.75d0*ILURows(i+1), (n-i)*(1.0d0*n) )

               ALLOCATE( IWork(k), STAT=istat )
               IF ( istat /= 0 ) THEN
                  CALL Fatal( 'CRS_ComplexILUT', 'Memory allocation error.' )
               END IF
               DO j=1,ILURows(i+1)-1
                 IWork(j) = ILUCols(j)
               END DO
               DEALLOCATE( ILUCols )

               ALLOCATE( CWork(k), STAT=istat )
               IF ( istat /= 0 ) THEN
                  CALL Fatal( 'CRS_ComplexILUT', 'Memory allocation error.' )
               END IF
               DO j=1,ILURows(i+1)-1
                 CWork(j) = ILUValues(j)
               END DO
               DEALLOCATE( ILUValues )

               ILUCols   => IWork
               ILUValues => CWork
            END IF
         END IF
      END DO
!
!     Prescale the diagonal for the LU solve:
!     ---------------------------------------
      DO i=1,n
         IF ( ABS(ILUValues(ILUDiag(i))) < AEPS ) THEN
            ILUValues(ILUDiag(i)) = 1.0d0
         ELSE
            ILUValues(ILUDiag(i)) = 1.0d0 / ILUValues(ILUDiag(i))
         END IF
      END DO

      A % ILUCols    => ILUCols
      A % CILUValues => ILUValues
      NULLIFY( ILUCols, ILUValues )
!------------------------------------------------------------------------------
    END SUBROUTINE ComplexComputeILUT
!------------------------------------------------------------------------------
  END FUNCTION CRS_ComplexILUT
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Incomplete factorization preconditioner solver for a CRS format matrix.
!>    Matrix is accessed from a global variable GlobalMatrix. The real version.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_LUPrecondition( u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(*), INTENT(IN) :: ipar  !< structure holding info from (HUTIter-iterative solver package)
    REAL(KIND=dp), DIMENSION(HUTI_NDIM), INTENT(IN) :: v   !< Right-hand-side vector
    REAL(KIND=dp), DIMENSION(HUTI_NDIM), INTENT(OUT) :: u   !< Solution vector

    u = v
    CALL CRS_LUSolve( HUTI_NDIM,GlobalMatrix,u )
  END SUBROUTINE CRS_LUPrecondition
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Incomplete factorization preconditioner solver for a CRS format matrix.
!>    Matrix is accessed from a global variable GlobalMatrix. The Complex version.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexLUPrecondition( u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(*), INTENT(IN) :: ipar  !< structure holding info from (HUTIter-iterative solver package)
    COMPLEX(KIND=dp), DIMENSION(HUTI_NDIM), INTENT(IN) :: v !< Right-hand-side vector
    COMPLEX(KIND=dp), DIMENSION(HUTI_NDIM), INTENT(OUT) :: u !< Solution vector

    u = v
    CALL CRS_ComplexLUSolve( HUTI_NDIM,GlobalMatrix,u )
  END SUBROUTINE CRS_ComplexLUPrecondition
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Solve a system (Ax=b) after factorization A=LUD has been done. This
!>    routine is meant as a part of  a preconditioner for an iterative solver.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_LUSolve( N,A,b )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(IN) :: A  !< Structure holding input matrix
    INTEGER, INTENT(IN) :: N   !< Size of the system
    DOUBLE PRECISION :: b(n)   !< on entry the RHS vector, on exit the solution vector.
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,row,col,nn
    DOUBLE PRECISION :: s
    DOUBLE PRECISION, POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------

    Diag => A % ILUDiag
    Rows => A % ILURows
    Cols => A % ILUCols
    Values => A % ILUValues

!
!   if no ilu provided do diagonal solve:
!   -------------------------------------
    IF ( .NOT. ASSOCIATED( Values ) ) THEN
       DO i=1,A % NumberOfRows
         b(i) = b(i) / A % Values( A % Diag(i) )
       END DO
       RETURN
    END IF

!--------------------------------------------------------------------
! The following #ifdefs  seem really necessery, if speed is an issue:
! SGI compiler optimizer  wants to know the sizes of the arrays very
! explicitely, while DEC compiler seems to make a copy of some of the
! arrays on the subroutine call (destroying performance).
!--------------------------------------------------------------------
#ifndef SGI
    IF ( A % Cholesky ) THEN
      !
      ! Forward substitute (solve z from Lz = b)
      DO i=1,n
        s = b(i)
        DO j=Rows(i),Diag(i)-1
           s = s - Values(j) * b(Cols(j))
        END DO
        b(i) = s * Values(Diag(i))
      END DO

      !
      ! Backward substitute (solve x from L^Tx = z)
      DO i=n,1,-1
        b(i) = b(i) * Values(Diag(i))
        DO j=Rows(i),Diag(i)-1
           b(Cols(j)) = b(Cols(j)) - Values(j) * b(i)
        END DO
      END DO
    ELSE
      !
      ! Forward substitute (solve z from Lz = b)
      DO i=1,n
         s = b(i)
         DO j=Rows(i),Diag(i)-1
            s = s - Values(j) * b(Cols(j))
         END DO
         b(i) = s 
      END DO

      !
      ! Backward substitute (solve x from UDx = z)
      DO i=n,1,-1
         s = b(i)
         DO j=Diag(i)+1,Rows(i+1)-1
            s = s - Values(j) * b(Cols(j))
         END DO
         b(i) = Values(Diag(i)) * s
      END DO
    END IF
#else
    CALL LUSolve( n,SIZE(Cols),Rows,Cols,Diag,Values,b )

  CONTAINS

    SUBROUTINE LUSolve( n,m,Rows,Cols,Diag,Values,b )
      INTEGER :: n,m,Rows(n+1),Cols(m),Diag(n)
      REAL(KIND=dp) :: Values(m),b(n)

      INTEGER :: i,j

      !
      ! Forward substitute (solve z from Lz = b)
      DO i=1,n
         DO j=Rows(i),Diag(i)-1
            b(i) = b(i) - Values(j) * b(Cols(j))
         END DO
      END DO

      !
      ! Backward substitute (solve x from UDx = z)
      DO i=n,1,-1
         DO j=Diag(i)+1,Rows(i+1)-1
            b(i) = b(i) - Values(j) * b(Cols(j))
         END DO
         b(i) = Values(Diag(i)) * b(i)
      END DO
    END SUBROUTINE LUSolve
#endif

  END SUBROUTINE CRS_LUSolve
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Solve a complex system Ax=b after factorization A=LUD has been
!>    done. This routine is meant as a part of  a preconditioner for an
!>    iterative solver.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexLUSolve( N,A,b )
!------------------------------------------------------------------------------ 
    TYPE(Matrix_t), INTENT(IN) :: A   !< Structure holding input matrix
    INTEGER, INTENT(IN) :: N          !< Size of the system
    COMPLEX(KIND=dp) :: b(N)          !< on entry the RHS vector, on exit the solution vector.
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), POINTER :: Values(:)
    INTEGER :: i,j
    COMPLEX(KIND=dp) :: x, s
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
!------------------------------------------------------------------------------

    Diag => A % ILUDiag
    Rows => A % ILURows
    Cols => A % ILUCols
    Values => A % CILUValues

!
!   if no ilu provided do diagonal solve:
!   -------------------------------------
    IF ( .NOT. ASSOCIATED( Values ) ) RETURN

!---------------------------------------------------------------------
! The following #ifdefs  seem really necessery, if speed is an issue:
! SGI compiler optimizer  wants to know the sizes of the arrays very
! explicitely, while DEC compiler seems to make a copy of some of the
! arrays on the subroutine call (destroying performance).
!--------------------------------------------------------------------
#ifndef SGI
    IF ( A % Cholesky ) THEN
      !
      ! Forward substitute
      DO i=1,n
         s = b(i)
         DO j=Rows(i),Diag(i)-1
            s = s - Values(j) * b(Cols(j))
         END DO
         b(i) = s * Values(Diag(i))
      END DO

      !
      ! Backward substitute
      DO i=n,1,-1
         b(i) = b(i) * Values(Diag(i))
         DO j=Rows(i),Diag(i)-1
           b(Cols(j)) = b(Cols(j)) - Values(j) * b(i)
         END DO
      END DO
    ELSE
      !
      ! Forward substitute
      DO i=1,n
         s = b(i)
         DO j=Rows(i),Diag(i)-1
            s = s - Values(j) * b(Cols(j))
         END DO
         b(i) = s
      END DO

      !
      ! Backward substitute
      DO i=n,1,-1
         s = b(i)
         DO j=Diag(i)+1,Rows(i+1)-1
            s = s - Values(j) * b(Cols(j))
         END DO
         b(i) = Values(Diag(i)) * s
      END DO
    END IF
#else
    CALL ComplexLUSolve( n,SIZE(Cols),Rows,Cols,Diag,Values,b )

  CONTAINS

    SUBROUTINE ComplexLUSolve( n,m,Rows,Cols,Diag,Values,b )
      INTEGER :: n,m,Rows(n+1),Cols(m),Diag(n)
      COMPLEX(KIND=dp) :: Values(m),b(n)

      INTEGER :: i,j

      !
      ! Forward substitute
      DO i=1,n
         DO j=Rows(i),Diag(i)-1
            b(i) = b(i) - Values(j) * b(Cols(j))
         END DO
      END DO

      !
      ! Backward substitute
      DO i=n,1,-1
         DO j=Diag(i)+1,Rows(i+1)-1
            b(i) = b(i) - Values(j) * b(Cols(j))
         END DO
         b(i) = Values(Diag(i)) * b(i)
      END DO
    END SUBROUTINE ComplexLUSolve
#endif

  END SUBROUTINE CRS_ComplexLUSolve
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Matrix vector product (v = Au) for a matrix given in CRS format. The
!>    matrix is accessed from a global variable GlobalMatrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_MatrixVectorProd( u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(*), INTENT(IN) :: ipar      !< Structure holding info HUTIter-iterative solver package
    REAL(KIND=dp), INTENT(IN) :: u(HUTI_NDIM)  !< vector to multiply u
    REAL(KIND=dp) :: v(HUTI_NDIM)  !< result vector

!------------------------------------------------------------------------------
    INTEGER, POINTER  CONTIG :: Cols(:),Rows(:)
!   INTEGER, POINTER :: Cols(:),Rows(:)
    REAL(KIND=dp), POINTER  CONTIG :: Values(:)
!   REAL(KIND=dp), POINTER :: Values(:)

#ifdef HAVE_MKL
	INTERFACE
		SUBROUTINE mkl_dcsrgemv(transa, m, a, ia, ja, x, y)
	 		USE Types
	 		CHARACTER :: transa
	 		INTEGER :: m
	 		REAL(KIND=dp) :: a(*)
	 		INTEGER :: ia(*), ja(*)
	 		REAL(KIND=dp) :: x(*), y(*)
	 	END SUBROUTINE mkl_dcsrgemv
	END INTERFACE
#endif

    INTEGER :: i,j,n
    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------

    n = GlobalMatrix % NumberOfRows
    Rows   => GlobalMatrix % Rows
    Cols   => GlobalMatrix % Cols
    Values => GlobalMatrix % Values

    IF  ( GlobalMatrix % MatVecSubr /= 0 ) THEN
#ifdef USE_ISO_C_BINDINGS
      CALL MatVecSubrExt(GlobalMatrix % MatVecSubr, &
                  GlobalMatrix % SpMV, n,Rows,Cols,Values,u,v,0)
#else
      CALL MatVecSubr(GlobalMatrix % MatVecSubr, &
                  GlobalMatrix % SpMV, n,Rows,Cols,Values,u,v,0)
#endif
      RETURN
   END IF
!--------------------------------------------------------------------
! The following #ifdefs  seem really necessery, if speed is an issue:
! SGI compiler optimizer  wants to know the sizes of the arrays very
! explicitely, while DEC compiler seems to make a copy of some of the
! arrays on the subroutine call (destroying performance).
!--------------------------------------------------------------------
#ifndef SGI
    IF ( HUTI_EXTOP_MATTYPE == HUTI_MAT_NOTTRPSED ) THEN
#ifdef HAVE_MKL
    CALL mkl_dcsrgemv('N', n, Values, Rows, Cols, u, v)
#else
!$omp parallel do private(j,s)
       DO i=1,n
          s = 0.0d0
          DO j=Rows(i),Rows(i+1)-1
             s = s + Values(j) * u(Cols(j))
          END DO
          v(i) = s
       END DO
!$omp end parallel do
#endif
    ELSE
       v(1:n) = 0.0d0
       DO i=1,n
          s = u(i)
          DO j=Rows(i),Rows(i+1)-1
             v(Cols(j)) = v(Cols(j)) + s * Values(j)
          END DO
       END DO
    END IF
!    IF ( ASSOCIATED( GlobalMatrix % EMatrix ) ) THEN
!       n = GlobalMatrix % EMatrix % NumberOFRows
!       Rows   => GlobalMatrix % EMatrix % Rows
!       Cols   => GlobalMatrix % EMatrix % Cols
!       Values => GlobalMatrix % EMatrix % Values
!
!       allocate( w(n) )
!       DO i=1,n
!          s = 0.0d0
!          DO j=Rows(i),Rows(i+1)-1
!             s = s + Values(j) * u(Cols(j))
!          END DO
!          w(i) = s
!       END DO
!
!       DO i=1,n
!          s = w(i)
!          DO j=Rows(i),Rows(i+1)-1
!             v(Cols(j)) = v(Cols(j)) + s * Values(j)
!          END DO
!       END DO
!       deallocate( w )
!    END IF
#else
    CALL MatVec( n,SIZE(Cols),Rows,Cols,Values,u,v )

  CONTAINS

    SUBROUTINE MatVec( n,m,Rows,Cols,Values,u,v )
      INTEGER :: n,m
      INTEGER :: Rows(n+1),Cols(m)
      REAL(KIND=dp) :: Values(m),u(n),v(n)

      INTEGER :: i,j

      IF ( HUTI_EXTOP_MATTYPE == HUTI_MAT_NOTTRPSED ) THEN
         v(1:n) = 0.0d0
         DO i=1,n
            DO j=Rows(i),Rows(i+1)-1
               v(i) = v(i) + Values(j) * u(Cols(j))
            END DO
         END DO
      ELSE
         v(1:n) = 0.0d0
         DO i=1,n
            s = u(i)
            DO j=Rows(i),Rows(i+1)-1
               v(Cols(j)) = v(Cols(j)) + s * Values(j)
            END DO
         END DO
      END IF
    END SUBROUTINE MatVec
#endif

  END SUBROUTINE CRS_MatrixVectorProd
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Complex matrix vector product (v = Au) for a matrix given in
!>    CRS format. The matrix is accessed from a global variable
!>    GlobalMatrix.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_ComplexMatrixVectorProd( u,v,ipar )
!------------------------------------------------------------------------------

    INTEGER, DIMENSION(*), INTENT(IN) :: ipar      !< Structure holding info HUTIter-iterative solver package
    COMPLEX(KIND=dp), INTENT(IN) :: u(HUTI_NDIM)      !< vector to multiply u
    COMPLEX(KIND=dp) :: v(HUTI_NDIM)  !< result vector

!------------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)
    INTEGER :: i,j,n
    COMPLEX(KIND=dp) :: s,rsum
    REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------

    n = HUTI_NDIM
    Rows   => GlobalMatrix % Rows
    Cols   => GlobalMatrix % Cols
    Values => GlobalMatrix % Values

!----------------------------------------------------------------------
! The following #ifdefs  seem really necessery, if speed is an issue:
! SGI compiler optimizer  wants to know the sizes of the arrays very
! explicitely, while DEC compiler seems to make a copy of some of the
! arrays on the subroutine call (destroying performance).
!----------------------------------------------------------------------
#ifndef SGI
    IF ( HUTI_EXTOP_MATTYPE == HUTI_MAT_NOTTRPSED ) THEN
!$omp parallel do private(rsum,j,s)
       DO i=1,n
          rsum = CMPLX( 0.0d0, 0.0d0, KIND=dp )
          DO j=Rows(2*i-1),Rows(2*i)-1,2
             s = CMPLX( Values(j), -Values(j+1), KIND=dp )
             rsum = rsum + s * u((Cols(j)+1)/2)
          END DO
          v(i) = rsum
       END DO
!$omp end parallel do
    ELSE
       v = CMPLX( 0.0d0, 0.0d0, KIND=dp )
       DO i=1,n
          rsum = u(i)
          DO j=Rows(2*i-1),Rows(2*i)-1,2
             s = CMPLX( Values(j), -Values(j+1), KIND=dp )
             v((Cols(j)+1)/2) = v((Cols(j)+1)/2) + s * rsum
          END DO
       END DO
    END IF
#else
    CALL ComplexMatVec( n,SIZE(Cols),Rows,Cols,Values,u,v )

  CONTAINS

    SUBROUTINE ComplexMatVec( n,m,Rows,Cols,Values,u,v )
      INTEGER :: n,m
      INTEGER :: Rows(2*n+1),Cols(m)
      REAL(KIND=dp) :: Values(m)
      COMPLEX(KIND=dp) :: u(n),v(n)

      INTEGER :: i,j
      COMPLEX(KIND=dp) :: s, rsum

      IF ( HUTI_EXTOP_MATTYPE == HUTI_MAT_NOTTRPSED ) THEN
         DO i=1,n
            rsum = CMPLX( 0.0d0, 0.0d0, KIND=dp )
            DO j=Rows(2*i-1),Rows(2*i)-1,2
               s = CMPLX( Values(j), -Values(j+1), KIND=dp )
               rsum = rsum + s * u((Cols(j)+1)/2)
            END DO
            v(i) = rsum
         END DO
      ELSE
         v = CMPLX( 0.0d0, 0.0d0, KIND=dp )
         DO i=1,n
            rsum = u(i)
            DO j=Rows(2*i-1),Rows(2*i)-1,2
               s = CMPLX( Values(j), -Values(j+1), KIND=dp )
               v((Cols(j)+1)/2) = v((Cols(j)+1)/2) + s * rsum
            END DO
         END DO
      END IF
    END SUBROUTINE ComplexMatVec
#endif

  END SUBROUTINE CRS_ComplexMatrixVectorProd
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Check the matrix for correctness.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_InspectMatrix( A )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A   !< Structure holding the matrix
!------------------------------------------------------------------------------

    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    INTEGER :: i,j,k,k2,n,ColMin,ColMax,RowMin,RowMax,ColN,RowN,ValN,UnSym
    REAL(KIND=dp) :: ValMin,ValMax,TotalSum,DiagSum
    LOGICAL :: Hit

    Diag   => A % Diag
    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    n = A % NumberOfRows
    RowMin = MINVAL( Rows )
    RowMax = MAXVAL( Rows )
    RowN = SIZE( Rows )

    PRINT *,'Rows:'
    PRINT *,'size:',RowN
    PRINT *,'min:',RowMin
    PRINT *,'max:',RowMax

    IF( RowMin < 1 ) THEN
      PRINT *,'Outliers:'
      j =  0
      DO i=1,RowN
        IF( Rows(i) < 1 ) THEN
          j = j + 1
          IF( j < 10 ) PRINT *,'i',i,Rows(i)
        END IF
      END DO
      IF( j > 0 ) THEN
        PRINT *,'Number of row outliers:',j
      END IF
    END IF

    ColMin = MINVAL( Cols )
    ColMax = MAXVAL( Cols )
    ColN = SIZE( Cols )

    PRINT *,'Cols:'
    PRINT *,'size:',ColN
    PRINT *,'min:',ColMin
    PRINT *,'max:',ColMax

    IF( ColMin < 1 ) THEN
      PRINT *,'Outliers:'
      j = 0
      DO i=1,ColN
        IF( Cols(i) < 1 ) THEN
          j = j + 1
          IF( j < 10 )PRINT *,'i',i,Cols(i)
        END IF
      END DO
      IF( j > 0 ) THEN
        PRINT *,'Number of column outliers:',j
      END IF
    END IF

    ValMin = MINVAL( Values )
    ValMax = MAXVAL( Values )
    ValN = SIZE( Values )
    TotalSum = SUM( Values )

    PRINT *,'Values:'
    PRINT *,'size:',ValN
    PRINT *,'min:',ValMin
    PRINT *,'max:',ValMax
    PRINT *,'sum:',TotalSum
    
    IF( ColN /= RowMax - 1  ) THEN
      PRINT *,'Conflicting max row index :',n,RowMax-1
    END IF
    IF( N /= ColMax  ) THEN
      PRINT *,'Conflicting max row index:',N,ColMax
    END IF
    IF( ColN /= ValN  ) THEN
      PRINT *,'Conflicting size of Values vs. Cols:',ValN,ColN
    END IF

    PRINT *,'Checking Diag vector:'
    j = 0
    DO i=1,n
      IF( i /= Cols(Diag(i)) ) THEN
        j = j + 1
        IF( j < 10 ) PRINT *,'diag:',i,Cols(Diag(i)) 
      END IF
    END DO
    IF( j > 0 ) THEN
      PRINT *,'Number of errornous diag entries:',j
    ELSE
      DiagSum = SUM( Values(Diag) )
      PRINT *,'diagonal sum:',DiagSum

! ratios are used to study cluster MG method, for example
      PRINT *,'ratios:',TotalSum/DiagSum,(TotalSum-DiagSum)/DiagSum
      
      TotalSum = SUM( ABS( Values ) )
      DiagSum = SUM( ABS( Values(Diag) ) )
      PRINT *,'abs ratios:',TotalSum/DiagSum,(TotalSum-DiagSum)/DiagSum

    END IF
    
    UnSym = 0
    DO i=1,n
      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        Hit = .FALSE.
        DO k2=Rows(j),Rows(j+1)-1
          IF( Cols(k2) == i ) THEN
            Hit = .TRUE.
            EXIT
          END IF
        END DO
        IF( .NOT. Hit ) THEN
          UnSym = UnSym + 1
          IF( UnSym < 10 ) PRINT *,'unsym:',i,j
        END IF
      END DO
    END DO

    PRINT *,'Number of unsymmetric entries in matrix topology:',UnSym


  END SUBROUTINE CRS_InspectMatrix


!------------------------------------------------------------------------------
!> Eliminate redundant entries in matrix A by summing same (i,j) combinations
!> together. It is assumed that the same entries are on the same CRS row.
!------------------------------------------------------------------------------
  SUBROUTINE CRS_PackMatrix( A )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A     !< Structure holding the matrix
!------------------------------------------------------------------------------ 
    INTEGER :: i,j,k,k2,n,rowi,nofs0
    INTEGER, ALLOCATABLE :: LocalCols(:), ColIndex(:)
    INTEGER, POINTER :: Cols(:),Rows(:),Diag(:)
    REAL(KIND=dp), POINTER :: Values(:), TValues(:)
    REAL(KIND=dp), ALLOCATABLE :: LocalValues(:), LocalTValues(:)
!------------------------------------------------------------------------------

    IF(A % NumberOfRows == 0 ) RETURN

    Rows    => A % Rows
    Cols    => A % Cols
    Diag    => A % Diag
    Values  => A % Values
    TValues => A % TValues

    nofs0 = Rows(A % NumberOfRows+1)-1

    n = 0
    DO i=1,A % NumberOfRows
      n = MAX( n, Rows(i+1)+1-Rows(i) )
    END DO

    ALLOCATE( LocalCols(n), LocalValues(n), LocalTValues(n) )
    LocalCols = 0
    LocalValues = 0.0_dp

    ALLOCATE(ColIndex(MAXVAL(Cols)))
    ColIndex = 0

    k2 = 0
    DO i=1,A % NumberOfRows
      Rowi = k2+1

      ! Memorize the matrix row
      DO k=1,Rows(i+1)-Rows(i)
        LocalCols(k)   = Cols(Rows(i)+k-1)
        LocalValues(k) = Values(Rows(i)+k-1)
        IF(ASSOCIATED(TValues)) &
          LocalTValues(k) = TValues(Rows(i)+k-1)
      END DO

      ! Pack the matrix row 
      DO k=1,Rows(i+1)-Rows(i)
        j = LocalCols(k) 
        IF( ColIndex(j) == 0 ) THEN
          k2 = k2 + 1
          ColIndex(j) = k2
          Cols(k2) = Cols(Rows(i)+k-1)
          Values(k2) = LocalValues(k)
          IF(ASSOCIATED(TValues)) &
            TValues(k2) = LocalTValues(k)
        ELSE
          Values(ColIndex(j)) = Values(ColIndex(j)) + LocalValues(k)
          IF(ASSOCIATED(TValues)) &
            TValues(ColIndex(j)) = TValues(ColIndex(j)) + LocalTValues(k)
        END IF
      END DO
      
      ! Nullify the index table
      DO k=1,Rows(i+1)-Rows(i)
        j = LocalCols(k) 
        ColIndex(j) = 0 
      END DO
      Rows(i) = Rowi
    END DO
    Rows(i) = k2+1

    CALL Info('CRS_PackMatrix','Number of summed-up matrix entries: '&
        //TRIM(I2S(nofs0-k2)),Level=8)

  END SUBROUTINE CRS_PackMatrix
!------------------------------------------------------------------------------



END MODULE CRSMatrix
!------------------------------------------------------------------------------

!> \} ElmerLib
