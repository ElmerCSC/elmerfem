!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{


MODULE ListMatrix

    USE CRSMatrix
    USE GeneralUtils

    IMPLICIT NONE

    INTEGER, PARAMETER :: LISTMATRIX_GROWTH = 1000

CONTAINS

!-------------------------------------------------------------------------------
!> Returns a handle to an allocated list matrix.
!-------------------------------------------------------------------------------
  FUNCTION List_AllocateMatrix(N) RESULT(Matrix)
!-------------------------------------------------------------------------------
    INTEGER :: i,n,istat
    TYPE(ListMatrix_t), POINTER :: Matrix(:)

    ALLOCATE( Matrix(n), STAT=istat )
    IF( istat /= 0 ) THEN
      CALL Fatal('List_AllocateMatrix','Allocation error for ListMatrix of size: '//I2S(n))
    END IF

    !$OMP PARALLEL
    !$OMP DO
    DO i=1,n
      Matrix(i) % Head => NULL()
      Matrix(i) % Level = 0
      Matrix(i) % Degree = 0
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
!-------------------------------------------------------------------------------
  END FUNCTION List_AllocateMatrix
!-------------------------------------------------------------------------------
 

!-------------------------------------------------------------------------------
!> Frees a list matrix.
!-------------------------------------------------------------------------------
   SUBROUTINE List_FreeMatrix( N, List )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: N
!-------------------------------------------------------------------------------

     TYPE(ListMatrixEntry_t), POINTER :: p,p1
     INTEGER :: i
!-------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED(List) ) RETURN

     !$OMP PARALLEL DO &
     !$OMP SHARED(List,N) &
     !$OMP PRIVATE(p, p1) DEFAULT(NONE)
     DO i=1,N
        p => List(i) % Head
        DO WHILE( ASSOCIATED(p) )
           p1 => p % Next
           DEALLOCATE( p )
           p => p1 
        END DO
     END DO
     !$OMP END PARALLEL DO
     DEALLOCATE( List )
!-------------------------------------------------------------------------------
   END SUBROUTINE List_FreeMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Enlarge a list matrix so that in can take in new entries.
!-------------------------------------------------------------------------------
  FUNCTION List_EnlargeMatrix(Matrix,N) RESULT(NewMatrix)
!-------------------------------------------------------------------------------
    INTEGER :: i,n
    TYPE(ListMatrix_t), POINTER :: Matrix(:), NewMatrix(:)

    NewMatrix => List_AllocateMatrix(n)
    IF ( ASSOCIATED(Matrix) ) THEN
       DO i=1,SIZE(Matrix)
        NewMatrix(i)=Matrix(i)
      END DO
      DEALLOCATE(Matrix)
    END IF
!-------------------------------------------------------------------------------
  END FUNCTION List_EnlargeMatrix
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Transfer the flexible list matrix to the more efficient CRS matrix that is 
!> used in most places of the code. Here the target is the rows and columns of the matrix.
!-------------------------------------------------------------------------------
  SUBROUTINE List_ToCRS(L,Rows,Cols,Diag)
!-------------------------------------------------------------------------------
    TYPE(ListMatrix_t) :: L(:)
    INTEGER :: i,j,n
    TYPE(Matrix_t), POINTER :: A
    TYPE(ListMatrixEntry_t), POINTER :: P
    INTEGER, POINTER CONTIG :: Rows(:),Cols(:),Diag(:)

    DO n=SIZE(L),1,-1
      IF ( L(n) % Degree>0 ) EXIT
    END DO

    ALLOCATE( Rows(n+1), Diag(n) )
    Rows(1) = 1
    DO i=1,n
      Rows(i+1) = Rows(i) + L(i) % Degree
    END DO
    ALLOCATE( Cols(Rows(i+1)-1) )
    j = 0
    DO i=1,n
      P => L(i) % Head
      DO WHILE(ASSOCIATED(P))
        j = j + 1
        Cols(j) = P % Index
        P => P % Next
      END DO
    END DO

    CALL Info('List_ToCRS',&
        'Number of entries in CRS matrix: '//I2S(Rows(n+1)-1),Level=8)

    A => AllocateMatrix()
    A % NumberOfRows = n
    A % Rows => Rows
    A % Diag => Diag
    A % Cols => Cols
    CALL CRS_SortMatrix(A)
    DEALLOCATE(A)
!-------------------------------------------------------------------------------
  END SUBROUTINE List_ToCRS
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Transfer the flexible list matrix to the more efficient CRS matrix that is 
!> used in most places of the code. The matrix structure can accommodate both forms.
!-------------------------------------------------------------------------------
  SUBROUTINE List_ToCRSMatrix(A)
!-------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    
    TYPE(ListMatrix_t), POINTER :: L(:)   
    INTEGER :: i,j,n
    TYPE(ListMatrixEntry_t), POINTER :: P
    INTEGER, POINTER CONTIG :: Rows(:),Cols(:),Diag(:)
    REAL(KIND=dp), POINTER CONTIG :: Values(:)

    IF( A % FORMAT /= MATRIX_LIST ) THEN
      CALL Warn('List_ToCRSMatrix','The initial matrix type is not List')
      RETURN
    END IF
    
    L => A % ListMatrix

    IF( .NOT. ASSOCIATED( L ) ) THEN
!     CALL Warn('ListToCRSMatrix','List not associated')
      A % FORMAT = MATRIX_CRS      
      A % NumberOfRows = 0
      RETURN
    END IF 
    
    DO n=SIZE(L),1,-1
      IF ( L(n) % Degree > 0 ) EXIT
    END DO
    
    ALLOCATE( Rows(n+1), Diag(n) )
    Diag = 0
    Rows(1) = 1
    DO i=1,n
      Rows(i+1) = Rows(i) + L(i) % Degree
    END DO

    CALL Info('List_ToCRSMatrix',&
        'Changing matrix type with number of non-zeros: '//I2S(Rows(n+1)-1),Level=8)

    ALLOCATE( Cols(Rows(n+1)-1)) 
    ALLOCATE( Values(Rows(n+1)-1) )

    j = 0
    DO i=1,n
      P => L(i) % Head
      DO WHILE(ASSOCIATED(P))
        j = j + 1
        Cols(j)   = P % Index
        Values(j) = P % Val
        P => P % Next
      END DO
    END DO
    
    A % NumberOfRows = n
    A % Rows => Rows
    A % Diag => Diag
    A % Cols => Cols
    A % Values => Values  
  
    A % Ordered=.FALSE.
    CALL CRS_SortMatrix( A )

    CALL List_FreeMatrix( SIZE(L), L )
    A % ListMatrix => NULL()

    A % FORMAT = MATRIX_CRS
    CALL Info('List_ToCRSMatrix','Matrix format changed from List to CRS', Level=12)

!-------------------------------------------------------------------------------
  END SUBROUTINE List_ToCRSMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Convert CRS matrix to list matrix
!-------------------------------------------------------------------------------
  SUBROUTINE List_ToListMatrix(A,Truncate)
!-------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    LOGICAL, OPTIONAL :: Truncate
    
    INTEGER :: i,j,n
    LOGICAL :: Trunc
    TYPE(ListMatrixEntry_t), POINTER :: CList

    Trunc=.FALSE.
    IF(PRESENT(Truncate)) Trunc=Truncate

    A % ListMatrix => List_AllocateMatrix(A % NumberOfRows)

    DO i=1,A % NumberOfRows
      A % ListMatrix(i) % Level  = 0
      A % ListMatrix(i) % Degree = 0

      IF(A % Rows(i) == A % Rows(i+1)) THEN
        A % ListMatrix(i) % Head => NULL()
        CYCLE
      END IF

      ALLOCATE(A % ListMatrix(i) % Head)
      Clist => A % ListMatrix(i) % Head
      Clist % Next => NULL()

      DO j=A % Rows(i), A % Rows(i+1)-1
        IF(Trunc) THEN
          IF (A % Cols(j) > A % NumberOfRows) EXIT
        END IF

        IF (j>A % Rows(i)) THEN
          IF ( Clist % Index >= A % Cols(j) ) THEN
            CALL Warn( 'List_ToListMatrix()', 'Input matrix not ordered ? ')
            GOTO 100
          END IF
          ALLOCATE(Clist % Next)
          Clist => Clist % Next
          CList % Next => NULL()
        END IF

        CList % Val = A % Values(j)
        CList % Index = A % Cols(j)
        A % ListMatrix(i) % Degree = A % ListMatrix(i) % Degree + 1
      END DO
    END DO

    GOTO 200

100 CONTINUE

    ! If not ordered input ...

    CALL List_FreeMatrix(i,A % ListMatrix)
    A % ListMatrix => Null()

    DO i=1,A % NumberOfRows
      DO j=A % Rows(i+1)-1,A % Rows(i),-1
        IF(Trunc) THEN
          IF (A % Cols(j) > A % NumberOfRows) CYCLE
        END IF
        CALL List_SetMatrixElement(A % ListMatrix,i,A % Cols(j),A % Values(j))
      END DO
    END DO

200 CONTINUE

    A % FORMAT = MATRIX_LIST

    IF( ASSOCIATED( A % Rows ) ) DEALLOCATE( A % Rows )
    IF( ASSOCIATED( A % Cols ) ) DEALLOCATE( A % Cols )
    IF( ASSOCIATED( A % Diag ) ) DEALLOCATE( A % Diag )
    IF( ASSOCIATED( A % Values ) ) DEALLOCATE( A % Values )

    A % Rows => Null()  
    A % Cols => Null()  
    A % Diag => Null()  
    A % Values => NULL()

    ! If the CRS matrix had a specific structure it is probably spoiled when going into
    ! free form matrix structure.
    A % ndeg = -1 
    
    CALL Info('List_ToListMatrix','Matrix format changed from CRS to List', Level=7)
!-------------------------------------------------------------------------------
  END SUBROUTINE List_ToListMatrix
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   FUNCTION List_GetMatrixIndex(List,k1,k2 ) RESULT(Entry)
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
     TYPE(ListMatrixEntry_t), POINTER :: CList,Prev, Entry, Dummy
!-------------------------------------------------------------------------------

     INTEGER :: i, istat

     
     IF ( .NOT. ASSOCIATED(List) ) List=>List_AllocateMatrix(k1)

     IF ( k1>SIZE(List) ) THEN
       List => List_EnlargeMatrix(List,MAX(k1, &
             SIZE(List)+LISTMATRIX_GROWTH) )
     END IF

     Clist => List(k1) % Head

     IF ( .NOT. ASSOCIATED(Clist) ) THEN
        Dummy => NULL()
        Entry => List_GetMatrixEntry(k2, Dummy )

        List(k1) % Degree = 1
        List(k1) % Head => Entry
        RETURN
     END IF

     NULLIFY( Prev )
     DO WHILE( ASSOCIATED(CList) )
        IF ( Clist % INDEX >= k2 ) EXIT
        Prev  => Clist
        CList => CList % Next
     END DO

     IF ( ASSOCIATED( CList ) ) THEN
        IF ( CList % INDEX == k2 ) THEN
          Entry => Clist
          RETURN
        END IF
     END IF

     Entry => List_GetMatrixEntry(k2, CList)
     IF ( ASSOCIATED( Prev ) ) THEN
         Prev % Next => Entry
     ELSE
        List(k1) % Head => Entry
     END IF
 
     List(k1) % Degree = List(k1) % Degree + 1
!-------------------------------------------------------------------------------
   END FUNCTION List_GetMatrixIndex
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE List_AddMatrixIndexes(List,k1,nk2,Ind)
   ! Add an array of sorted indeces to a row in ListMatrix_t. "ind" may
   ! contain duplicate entries.
!-------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(ListMatrix_t) :: List(:)
     INTEGER, INTENT(IN) :: k1, nk2
     INTEGER, INTENT(IN) :: Ind(nk2)

     TYPE(ListMatrixEntry_t), POINTER :: RowPtr, PrevPtr, Entry, Dummy
!-------------------------------------------------------------------------------
     INTEGER :: i,k2,k2i,j, k,prevind

     IF (k1>SIZE(List)) THEN
       CALL Fatal('List_AddMatrixIndexes','Row index out of bounds: '//TRIM(I2S(k1)))
     END IF
     
     ! Add each element in Ind to the row list
     RowPtr => List(k1) % Head
    
     ! First element needs special treatment as it may modify 
     ! the list starting point
     IF (.NOT. ASSOCIATED(RowPtr)) THEN
       Dummy => NULL() 
       Entry => List_GetMatrixEntry(Ind(1),Dummy)
       List(k1) % Degree = 1
       List(k1) % Head => Entry
       k2i = 2
       prevind = ind(1)
     ELSE IF (RowPtr % Index > Ind(1)) THEN
       Entry => List_GetMatrixEntry(Ind(1),RowPtr)
       List(k1) % Degree = List(k1) % Degree + 1
       List(k1) % Head => Entry
       k2i = 2
       prevind = ind(1)
     ELSE IF (RowPtr % Index == Ind(1)) THEN
        k2i = 2
        prevind = ind(1)
     ELSE
       k2i = 1
       prevind = -1
     END IF

     PrevPtr => List(k1) % Head
     RowPtr  => List(k1) % Head % Next

     DO i=k2i,nk2
       k2=Ind(i)
       if (k2 == prevind) cycle

       ! Find a correct place place to add index to
       DO WHILE( ASSOCIATED(RowPtr) )
         IF (RowPtr % Index >= k2) EXIT
         PrevPtr => RowPtr
         RowPtr  => RowPtr % Next
       END DO
       
       IF (ASSOCIATED(RowPtr)) THEN
         ! Do not add duplicates
         IF (RowPtr % Index /= k2) THEN
           ! Create new element between PrevPtr and RowPtr
           Entry => List_GetMatrixEntry(k2,RowPtr)
           PrevPtr % Next => Entry
           List(k1) % Degree = List(k1) % Degree + 1

           ! Advance to next element in list
           PrevPtr => Entry
!          RowPtr  => 
         ELSE
           ! Advance to next element in list
           PrevPtr => RowPtr
           RowPtr  => RowPtr % Next
         END IF
       ELSE
         EXIT
       END IF

       prevind = k2
     END DO

     DO j=i,nk2
       k2 = Ind(j)
       if (k2 == prevind) cycle
       prevind = k2

       Dummy => NULL()
       Entry => List_GetMatrixEntry(k2,Dummy)
       PrevPtr % Next => Entry
       PrevPtr => PrevPtr % Next
       List(k1) % Degree = List(k1) % Degree + 1
     END DO
!-------------------------------------------------------------------------------
   END SUBROUTINE List_AddMatrixIndexes
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_AddMatrixRow(List,k1,nk2,Ind,Vals)
   ! Add an array of sorted indeces to a row in ListMatrix_t. "ind" may
   ! contain duplicate entries.
!-------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(ListMatrix_t) :: List(:)
     INTEGER, INTENT(IN) :: k1, nk2
     INTEGER, INTENT(INOUT) :: Ind(nk2)
     REAL(KIND=dp), INTENT(INOUT) :: Vals(nk2)

     TYPE(ListMatrixEntry_t), POINTER :: RowPtr, PrevPtr, Entry, Dummy
!-------------------------------------------------------------------------------
     INTEGER :: i,k2,k2i,j, k,prevind

     IF (k1>SIZE(List)) THEN
       CALL Fatal('List_AddMatrixIndexes','Row index out of bounds: '//TRIM(I2S(k1)))
     END IF

     CALL SortF(nk2, Ind, Vals)
     
     ! Add each element in Ind to the row list
     RowPtr => List(k1) % Head
    
     ! First element needs special treatment as it may modify 
     ! the list starting point
     IF (.NOT. ASSOCIATED(RowPtr)) THEN
       Dummy => NULL() 
       Entry => List_GetMatrixEntry(Ind(1),Dummy)
       Entry % Val = Vals(1)
       List(k1) % Degree = 1
       List(k1) % Head => Entry
       k2i = 2
       prevind = ind(1)
     ELSE IF (RowPtr % Index > Ind(1)) THEN
       Entry => List_GetMatrixEntry(Ind(1),RowPtr)
       Entry % Val = Vals(1)
       List(k1) % Degree = List(k1) % Degree + 1
       List(k1) % Head => Entry
       k2i = 2
       prevind = ind(1)
     ELSE IF (RowPtr % Index == Ind(1)) THEN
        k2i = 2
        prevind = ind(1)
     ELSE
       k2i = 1
       prevind = -1
     END IF

     PrevPtr => List(k1) % Head
     RowPtr  => List(k1) % Head % Next

     DO i=k2i,nk2
       k2 = Ind(i)
       if (k2 == prevind) cycle

       ! Find a correct place place to add index to
       DO WHILE( ASSOCIATED(RowPtr) )
         IF (RowPtr % Index >= k2) EXIT
         PrevPtr => RowPtr
         RowPtr  => RowPtr % Next
       END DO
       
       IF (ASSOCIATED(RowPtr)) THEN
         ! Do not add duplicates
         IF (RowPtr % Index /= k2) THEN
           ! Create new element between PrevPtr and RowPtr
           Entry => List_GetMatrixEntry(k2,RowPtr)
           Entry % Val = Vals(i)
           PrevPtr % Next => Entry
           List(k1) % Degree = List(k1) % Degree + 1

           ! Advance to next element in list
           PrevPtr => Entry
         ELSE
           ! Advance to next element in list
           RowPtr % Val = RowPtr % Val + Vals(i)
           PrevPtr => RowPtr
           RowPtr  => RowPtr % Next
         END IF
       ELSE
         EXIT
       END IF

       prevind = k2
     END DO

     DO j=i,nk2
       k2 = Ind(j)
       if (k2 == prevind) cycle
       prevind = k2

       Dummy => NULL()
       Entry => List_GetMatrixEntry(k2,Dummy)
       Entry % Val = Vals(j)
       PrevPtr % Next => Entry
       PrevPtr => PrevPtr % Next
       List(k1) % Degree = List(k1) % Degree + 1
     END DO
!-------------------------------------------------------------------------------
   END SUBROUTINE List_AddMatrixRow
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
   FUNCTION List_GetMatrixEntry(ind, next) RESULT(ListEntry)
!-------------------------------------------------------------------------------
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: ind
     TYPE(ListMatrixEntry_t), POINTER, INTENT(IN) :: next
     TYPE(ListMatrixEntry_t), POINTER :: ListEntry

     INTEGER :: istat
     
     ALLOCATE(ListEntry, STAT=istat)
     IF( istat /= 0 ) THEN
        CALL Fatal('List_GetMatrixEntry','Could not allocate entry!')
     END IF

     ListEntry % Val = REAL(0,dp)
     ListEntry % INDEX = ind
     ListEntry % Next => next
!-------------------------------------------------------------------------------
   END FUNCTION List_GetMatrixEntry     
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE List_DeleteMatrixElement(List,k1,k2)
!-------------------------------------------------------------------------------
     INTEGER :: k1,k2
     TYPE(ListMatrix_t) :: List(:)
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: Clist,Prev

     Prev => NULL()
     Clist => List(k1) % Head
     DO WHILE(ASSOCIATED(Clist))
       IF (Clist % Index >= k2) EXIT
       Prev  => Clist
       Clist => Clist % Next
     END DO
     IF (.NOT.ASSOCIATED(Clist)) RETURN

     IF (Clist % Index /= k2) RETURN

     IF (ASSOCIATED(Prev)) THEN
       Prev % Next => Clist % Next
     ELSE
       List(k1) % Head => Clist % Next
     END IF
     DEALLOCATE(Clist)
     List(k1) % Degree = MAX(List(k1) % Degree-1,0)
!-------------------------------------------------------------------------------
   END SUBROUTINE List_DeleteMatrixElement
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_DeleteRow(List,k1,Keep)
!-------------------------------------------------------------------------------
     INTEGER :: k1,k2
     LOGICAL, OPTIONAL :: Keep
     TYPE(ListMatrix_t) :: List(:)
!-------------------------------------------------------------------------------
     LOGICAL :: lKeep
     INTEGER::n
     TYPE(ListMatrixEntry_t), POINTER :: Clist,Next

     n = SIZE(List)
     IF(k1<=0.OR.k1>n) RETURN

     Clist=>List(k1) % Head
     DO WHILE(ASSOCIATED(Clist))
       Next=>Clist % Next
       DEALLOCATE(Clist)
       Clist=>Next
     END DO

     lKeep = .FALSE.
     IF(PRESENT(Keep)) lKeep = Keep
     
     IF(lKeep) THEN
       List(k1) % Degree=0
       List(k1) % Head=>NULL()
     ELSE
       List(k1:n-1)=List(k1+1:n)
       List(n) % Degree=0
       List(n) % Head=>NULL()
     END IF
!-------------------------------------------------------------------------------
   END SUBROUTINE List_DeleteRow
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_DeleteCol(List,k1)
!-------------------------------------------------------------------------------
     INTEGER :: k1
     TYPE(ListMatrix_t) :: List(:)
!-------------------------------------------------------------------------------
     INTEGER::i,n
     TYPE(ListMatrixEntry_t), POINTER :: Clist,Prev

     n=SIZE(List)

     DO i=1,n
       Prev => NULL()
       Clist => List(i) % Head
       DO WHILE(ASSOCIATED(Clist))
         IF(Clist % Index>=k1) EXIT
         Prev  => Clist
         Clist => Clist % Next
       END DO

       IF(.NOT.ASSOCIATED(Clist)) CYCLE

       IF (Clist % Index==k1) THEN
         IF(ASSOCIATED(Prev)) THEN
           Prev % Next => Clist % Next
         ELSE
           List(i) % Head => Clist % Next
         END IF
         List(i) % Degree = MAX(List(i) % Degree-1,0)
         DEALLOCATE(Clist)
       END IF
     END DO
!-------------------------------------------------------------------------------
   END SUBROUTINE List_DeleteCol
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_AddToMatrixElement( List,k1,k2,Val,SetVal )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
     REAL(KIND=dp) :: Val
     LOGICAL, OPTIONAL :: SetVal 
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: CList,Prev, Entry
     LOGICAL :: Set     

     Set = .FALSE.
     IF( PRESENT(SetVal)) Set = SetVal

     Entry => List_GetMatrixIndex(List,k1,k2)
     IF ( Set ) THEN
       Entry % Val = Val
     ELSE
       Entry % Val = Entry % Val + Val
     END IF
!-------------------------------------------------------------------------------
   END SUBROUTINE List_AddToMatrixElement
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE List_AddMatrixIndex( List,k1,k2  )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: CList,Prev, Entry
     LOGICAL :: Set     

     Entry => List_GetMatrixIndex(List,k1,k2)
!-------------------------------------------------------------------------------
   END SUBROUTINE List_AddMatrixIndex
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_SetMatrixElement( List,k1,k2,Val,SetVal )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
     TYPE(ListMatrixEntry_t), POINTER :: CList,Prev, Entry
     REAL(KIND=dp) :: Val
     LOGICAL, OPTIONAL :: SetVal 

     CALL List_AddToMatrixElement( List,k1,k2,Val,.TRUE.)
!-------------------------------------------------------------------------------
   END SUBROUTINE List_SetMatrixElement
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   FUNCTION List_GetMatrixElement( List,k1,k2 ) RESULT ( Val )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
     TYPE(ListMatrixEntry_t), POINTER :: CList,Prev, Entry
     REAL(KIND=dp) :: Val
!-------------------------------------------------------------------------------

     Val = 0.0_dp

     IF ( .NOT. ASSOCIATED(List) ) RETURN
     IF ( k1>SIZE(List) ) RETURN
     Clist => List(k1) % Head
     IF ( .NOT. ASSOCIATED(Clist) ) RETURN

     NULLIFY( Prev )
     DO WHILE( ASSOCIATED(CList) )
        IF ( Clist % INDEX == k2 ) Val = CList % Val
        IF ( Clist % INDEX >= k2 ) RETURN
        Prev  => Clist
        CList => CList % Next
     END DO
!-------------------------------------------------------------------------------
   END FUNCTION List_GetMatrixElement
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_ZeroRow( List,k1 )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: CList
     
     IF ( .NOT. ASSOCIATED(List) ) THEN
       CALL Warn('List_ZeroRow','No List matrix present!')
       RETURN
     END IF
     
     IF ( k1 > SIZE(List) ) THEN
       CALL Warn('List_ZeroRow','No such row!')
       RETURN
     END IF
     
     Clist => List(k1) % Head
     IF ( .NOT. ASSOCIATED(Clist) ) THEN
       CALL Warn('List_ZeroRow','Row not associated!')
       RETURN
     END IF
     
     DO WHILE( ASSOCIATED(CList) )
       Clist % Val = 0.0_dp
       CList => CList % Next
     END DO
!-------------------------------------------------------------------------------
   END SUBROUTINE List_ZeroRow
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_MoveRow( List,n1,n2,coeff,staycoeff )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: n1, n2
     REAL(KIND=dp), OPTIONAL :: coeff, staycoeff
!-------------------------------------------------------------------------------
     INTEGER :: k2
     REAL(KIND=dp) :: val, c, d
     TYPE(ListMatrixEntry_t), POINTER :: CList

     IF( PRESENT(coeff)) THEN
       c = coeff
     ELSE
       c = 1.0_dp
     END IF

     IF( PRESENT(staycoeff)) THEN
       d = staycoeff
     ELSE
       d = 0.0_dp
     END IF
              
     IF ( .NOT. ASSOCIATED(List) ) THEN
       CALL Warn('List_MoveRow','No List matrix present!')
       RETURN
     END IF
     
     IF ( n1 > SIZE(List) ) THEN
       CALL Warn('List_MoveRow','No row to move!')
       RETURN
     END IF
     
     Clist => List(n1) % Head
     IF ( .NOT. ASSOCIATED(Clist) ) THEN
       CALL Warn('List_MoveRow','Row not associated!')
       RETURN
     END IF
     
     DO WHILE( ASSOCIATED(CList) )
       k2 = Clist % Index
       Val = Clist % Val
       Clist % VAL = d * Val 

! This could be made more optimal as all the entries are for the same row!
       CALL List_AddToMatrixElement(List,n2,k2,c*Val)

       CList => CList % Next
     END DO

!-------------------------------------------------------------------------------
   END SUBROUTINE List_MoveRow
!-------------------------------------------------------------------------------


! Exchange row structure between two matrix rows.
! Currently this is not optimal since we copy the structure back-and-forth.
!-------------------------------------------------------------------------------
   SUBROUTINE List_ExchangeRowStructure( List,n1,n2 )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: n1, n2
!-------------------------------------------------------------------------------
     INTEGER :: k1, k2
     TYPE(ListMatrixEntry_t), POINTER :: CList1, CList2, Lptr
              
     IF ( .NOT. ASSOCIATED(List) ) THEN
       CALL Warn('List_ExchangeRowStructure','No List matrix present!')
       RETURN
     END IF
         
     Clist1 => List(n1) % Head
     IF ( .NOT. ASSOCIATED(Clist1) ) THEN
       CALL Warn('List__ExchangeRowStructure','Row1 not associated!')
       RETURN
     END IF

     Clist2 => List(n2) % Head
     IF ( .NOT. ASSOCIATED(Clist2) ) THEN
       CALL Warn('List__ExchangeRowStructure','Row2 not associated!')
       RETURN
     END IF
     
     DO WHILE( ASSOCIATED(CList1) )
       k1 = Clist1 % Index
       Lptr => List_GetMatrixIndex( List,n2,k1 )
       CList1 => CList1 % Next
     END DO
     
     DO WHILE( ASSOCIATED(CList2) )
       k2 = Clist2 % Index
       Lptr => List_GetMatrixIndex( List,n1,k2 )
       CList2 => CList2 % Next
     END DO
     
!-------------------------------------------------------------------------------
   END SUBROUTINE List_ExchangeRowStructure
!-------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Add the entries of a local matrix to a list-format matrix.    
!------------------------------------------------------------------------------
  SUBROUTINE List_GlueLocalMatrix( A,N,Dofs,Indexes,LocalMatrix )
!------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: A(:)
     INTEGER :: N,DOFs, Indexes(:)
     REAL(KIND=dp) :: LocalMatrix(:,:)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Val
     INTEGER :: i,j,k,l,c,Row,Col
     
     DO i=1,n
       IF (Indexes(i)<=0) CYCLE
       DO k=0,Dofs-1
         Row = Dofs*Indexes(i)-k
         DO j=1,n
           IF (Indexes(j)<=0) CYCLE
           DO l=0,Dofs-1
             Col = Dofs * Indexes(j) - l
             Val = LocalMatrix(Dofs*i-k,Dofs*j-l)
             CALL List_AddToMatrixElement(A,Row,Col,Val)
           END DO
         END DO

       END DO
     END DO
   END SUBROUTINE List_GlueLocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Add the entries of a local matrix to a list-format matrix by allowing
!>    offsets
!------------------------------------------------------------------------------
   SUBROUTINE List_GlueLocalSubMatrix( List,row0,col0,Nrow,Ncol, &
          RowInds,ColInds,RowDofs,ColDofs,LocalMatrix )
!------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:) 
     INTEGER :: Nrow,Ncol,RowDofs,ColDofs,Col0,Row0,RowInds(:),ColInds(:)
     REAL(KIND=dp) :: LocalMatrix(:,:)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Val
     INTEGER :: i,j,k,l,c,Row,Col
     
     DO i=1,Nrow
       DO k=0,RowDofs-1
         IF ( RowInds(i) <= 0 ) CYCLE
         Row = Row0 + RowDofs * RowInds(i) - k
         
         DO j=1,Ncol
           DO l=0,ColDofs-1
             IF ( ColInds(j) <= 0 ) CYCLE
             Col  = Col0 + ColDofs * ColInds(j) - l
             Val = LocalMatrix(RowDofs*i-k,ColDofs*j-l)
             CALL List_AddToMatrixElement(List,Row,Col,Val)
           END DO
         END DO

       END DO
     END DO
   END SUBROUTINE List_GlueLocalSubMatrix
!------------------------------------------------------------------------------

END MODULE ListMatrix

!> \} ElmerLib


