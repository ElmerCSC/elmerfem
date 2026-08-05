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

    ! Matrix entries are handed out from a pool of large chunks instead of being
    ! allocated one at a time. This saves one allocation per matrix nonzero and,
    ! more importantly, keeps entries created in sequence close to each other in
    ! memory, which is what the row walks mostly pay for. Deleted entries go on
    ! a free list and are recycled, so the pool stays at the high water mark of
    ! simultaneously live entries rather than growing with the total ever built.
    ! Once every entry has been returned all chunks but the newest are released.
    !
    ! NOTE: like the rest of this module the pool assumes that a list matrix is
    ! built by one thread at a time. Threaded assembly has its own per-thread
    ! pools in ListMatrixArray.
    INTEGER, PARAMETER, PRIVATE :: LISTMATRIX_CHUNKMIN = 4096
    INTEGER, PARAMETER, PRIVATE :: LISTMATRIX_CHUNKMAX = 262144

    TYPE(ListMatrixEntryPool_t), POINTER, PRIVATE, SAVE :: EntryChunks => NULL()
    TYPE(ListMatrixEntry_t), POINTER, PRIVATE, SAVE :: FreeEntries => NULL()
    INTEGER(KIND=8), PRIVATE, SAVE :: EntriesInUse = 0

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

     TYPE(ListMatrixEntry_t), POINTER :: p,Head,Tail,RowTail
     INTEGER :: i
     INTEGER(KIND=8) :: cnt
!-------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED(List) ) RETURN

     ! Chain the rows of each thread together, then hand the whole chain back to
     ! the entry pool in a single splice. Nothing is freed per entry any more.
     !$OMP PARALLEL &
     !$OMP SHARED(List,N) &
     !$OMP PRIVATE(i,p,Head,Tail,RowTail,cnt) DEFAULT(NONE)
     Head => NULL()
     Tail => NULL()
     cnt = 0
     !$OMP DO
     DO i=1,N
        p => List(i) % Head
        IF ( .NOT. ASSOCIATED(p) ) CYCLE

        RowTail => p
        cnt = cnt + 1
        DO WHILE( ASSOCIATED( RowTail % Next ) )
           RowTail => RowTail % Next
           cnt = cnt + 1
        END DO

        ! The first row processed provides the tail of the whole chain
        IF ( .NOT. ASSOCIATED(Tail) ) Tail => RowTail

        RowTail % Next => Head
        Head => p
     END DO
     !$OMP END DO

     ! Unnamed !$OMP CRITICAL on purpose, as in GaussPointsInit: a named critical
     ! section allocates its per-name mutex lazily on first encounter, and that
     ! lazy init crashes inside MPI-spawned processes with GOMP on Windows+MSMPI.
     ! The unnamed form uses a pre-allocated global lock. Do not reintroduce a
     ! name here. The wider lock costs nothing, this being entered once per
     ! thread per matrix rather than once per entry.
     !$OMP CRITICAL
     CALL List_ReturnEntryChain( Head, Tail, cnt )
     !$OMP END CRITICAL
     !$OMP END PARALLEL

     DEALLOCATE( List )
     CALL List_ShrinkEntryPool()
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
    INTEGER :: i,j,n,istat
    TYPE(Matrix_t), POINTER :: A
    TYPE(ListMatrixEntry_t), POINTER :: P
    INTEGER, POINTER CONTIG :: Rows(:),Cols(:),Diag(:)

    DO n=SIZE(L),1,-1
      IF ( L(n) % Degree>0 ) EXIT
    END DO

    ALLOCATE( Rows(n+1), Diag(n), STAT=istat)
    IF(istat /= 0 ) THEN
      CALL Fatal('List_ToCRS','Could not allocate memory for CRS Rows of size '//I2S(n))
    END IF

    
    Rows(1) = 1
    DO i=1,n
      Rows(i+1) = Rows(i) + L(i) % Degree
    END DO
    ALLOCATE( Cols(Rows(i+1)-1), Stat=istat)
    IF(istat /= 0 ) THEN
      CALL Fatal('List_ToCRS','Could not allocate memory for CRS Cols of size '//I2S(Rows(i+1)-1))
    END IF


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
    INTEGER :: i,j,k,n,m,istat
    TYPE(ListMatrixEntry_t), POINTER :: P
    INTEGER, POINTER CONTIG :: Rows(:),Cols(:),Diag(:)
    REAL(KIND=dp), POINTER CONTIG :: Values(:)

    IF( A % FORMAT /= MATRIX_LIST ) THEN
      CALL Warn('List_ToCRSMatrix','The initial matrix type is not List')
      RETURN
    END IF
    
    L => A % ListMatrix

    IF( .NOT. ASSOCIATED( L ) ) THEN
      A % FORMAT = MATRIX_CRS      
      A % NumberOfRows = 0
      RETURN
    END IF 
    
    DO n=SIZE(L),1,-1
      IF ( L(n) % Degree > 0 ) EXIT
    END DO
    CALL Info('List_ToCRSMatrix','List size '//I2S(SIZE(L))//' vs. active rows '//I2S(n),Level=25)
    
    ALLOCATE( Rows(n+1), Diag(n), STAT=istat)
    IF(istat /= 0 ) THEN
      CALL Fatal('List_ToCRSMatrix','Could not allocate memory for CRS Rows of size '//I2S(n))
    END IF
    
    Diag = 0
    Rows(1) = 1
    DO i=1,n
      Rows(i+1) = Rows(i) + L(i) % Degree
    END DO

    m = Rows(n+1)-1
    CALL Info('List_ToCRSMatrix',&
        'Changing matrix type with number of non-zeros: '//I2S(m),Level=8)

    ALLOCATE( Cols(m),Values(m),STAT=istat) 
    IF(istat /= 0 ) THEN
      CALL Fatal('List_ToCRS','Could not allocate memory for CRS Cols & Values of size '//I2S(m))
    END IF

    j = 0
    DO i=1,n
      P => L(i) % Head
      DO WHILE(ASSOCIATED(P))
        j = j + 1
        Cols(j) = P % Index
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
    TYPE(ListMatrixEntry_t), POINTER :: CList, Dummy

    Trunc=.FALSE.
    IF(PRESENT(Truncate)) Trunc=Truncate
    Dummy => NULL()

    A % ListMatrix => List_AllocateMatrix(A % NumberOfRows)

    DO i=1,A % NumberOfRows
      A % ListMatrix(i) % Level  = 0
      A % ListMatrix(i) % Degree = 0

      IF(A % Rows(i) == A % Rows(i+1)) THEN
        A % ListMatrix(i) % Head => NULL()
        CYCLE
      END IF

      Clist => NULL()

      DO j=A % Rows(i), A % Rows(i+1)-1
        IF(Trunc) THEN
          IF (A % Cols(j) > A % NumberOfRows) EXIT
        END IF

        IF ( ASSOCIATED(Clist) ) THEN
          IF ( Clist % Index >= A % Cols(j) ) THEN
            CALL Warn( 'List_ToListMatrix()', 'Input matrix not ordered ? ')
            GOTO 100
          END IF
          Clist % Next => List_GetMatrixEntry( A % Cols(j), Dummy )
          Clist => Clist % Next
        ELSE
          A % ListMatrix(i) % Head => List_GetMatrixEntry( A % Cols(j), Dummy )
          Clist => A % ListMatrix(i) % Head
        END IF

        CList % Val = A % Values(j)
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

     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER, INTENT(IN) :: k1, nk2
     INTEGER, INTENT(IN) :: Ind(nk2)

     TYPE(ListMatrixEntry_t), POINTER :: RowPtr, PrevPtr, Entry, Dummy
!-------------------------------------------------------------------------------
     INTEGER :: i,k2,k2i,j,prevind

     IF (k1>SIZE(List)) THEN
       List => List_EnlargeMatrix(List,MAX(k1, &
             SIZE(List)+LISTMATRIX_GROWTH) )
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
   SUBROUTINE List_AddMatrixRow(List,k1,nk2,Ind,Vals,SortedInput,KeepOrder)
   ! Add an array of sorted indeces to a row in ListMatrix_t. "ind" may
   ! contain duplicate entries.
!-------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(ListMatrix_t), INTENT(INOUT), POINTER :: List(:)
     INTEGER, INTENT(IN) :: k1, nk2
     INTEGER, INTENT(INOUT) :: Ind(nk2)
     REAL(KIND=dp), INTENT(INOUT) :: Vals(nk2)
     LOGICAL, INTENT(IN), OPTIONAL :: SortedInput
     LOGICAL, INTENT(IN), OPTIONAL :: KeepOrder
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: RowPtr, PrevPtr, ENTRY, Dummy
     INTEGER, ALLOCATABLE :: OrigInd(:)
     REAL(KIND=dp), ALLOCATABLE :: OrigVals(:)
     LOGICAL :: DoSort, DoOrder
     INTEGER :: i,k2,k2i,j,prevind

     DoSort = .TRUE.
     DoOrder = .FALSE.
     IF( PRESENT(SortedInput) ) DoSort = ( .NOT. SortedInput  )

     ! Note that sorting spoils the order of Ind and Vals!
     IF(DoSort) THEN
       IF(PRESENT(KeepOrder)) DoOrder = KeepOrder
       IF( DoOrder ) THEN
         OrigInd = Ind(1:nk2)
         OrigVals = Vals(1:nk2)
       END IF
       CALL SortF(nk2, Ind, Vals)
     END IF

       
     IF(ASSOCIATED(List)) THEN
       i = SIZE(List)
     ELSE
       i = 0
     END IF
     
     IF (k1>i) THEN
       List => List_EnlargeMatrix(List,MAX(k1, &
             i+LISTMATRIX_GROWTH) )
     END IF
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
       prevind = IND(1)
     ELSE IF (RowPtr % Index == Ind(1)) THEN
        k2i = 2
        prevind = ind(1)
        RowPtr % Val = RowPtr % Val + Vals(1)
     ELSE
       k2i = 1
       prevind = -1
     END IF

     PrevPtr => List(k1) % Head
     RowPtr  => List(k1) % Head % Next

     DO i=k2i,nk2
       k2 = Ind(i)
       IF(k2 == prevind) THEN
          CYCLE
       END IF

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
       IF (k2 == prevind) THEN
         cycle
       END IF
       prevind = k2

       Dummy => NULL()
       Entry => List_GetMatrixEntry(k2,Dummy)
       Entry % Val = Vals(j)
       PrevPtr % Next => Entry
       PrevPtr => PrevPtr % Next
       List(k1) % Degree = List(k1) % Degree + 1
     END DO

     IF( DoOrder ) THEN
       Ind(1:nk2) = OrigInd 
       Vals(1:nk2) = OrigVals 
     END IF
     
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

     IF( ASSOCIATED( FreeEntries ) ) THEN
       ! Recycle a previously deleted entry
       ListEntry => FreeEntries
       FreeEntries => ListEntry % Next
     ELSE
       IF( .NOT. ASSOCIATED( EntryChunks ) ) THEN
         CALL List_NewEntryChunk()
       ELSE IF( EntryChunks % NextIndex > SIZE( EntryChunks % Entries ) ) THEN
         CALL List_NewEntryChunk()
       END IF
       ListEntry => EntryChunks % Entries( EntryChunks % NextIndex )
       EntryChunks % NextIndex = EntryChunks % NextIndex + 1
     END IF

     EntriesInUse = EntriesInUse + 1

     ListEntry % Val = REAL(0,dp)
     ListEntry % INDEX = ind
     ListEntry % Next => next
!-------------------------------------------------------------------------------
   END FUNCTION List_GetMatrixEntry
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Add a new chunk of entries to the pool. The chunks grow geometrically so that
!> even a very large matrix needs only a handful of them.
!-------------------------------------------------------------------------------
   SUBROUTINE List_NewEntryChunk()
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntryPool_t), POINTER :: Chunk
     INTEGER :: n, istat

     n = LISTMATRIX_CHUNKMIN
     IF( ASSOCIATED( EntryChunks ) ) THEN
       n = MIN( 2*SIZE( EntryChunks % Entries ), LISTMATRIX_CHUNKMAX )
     END IF

     ALLOCATE( Chunk, STAT=istat )
     IF( istat == 0 ) ALLOCATE( Chunk % Entries(n), STAT=istat )
     IF( istat /= 0 ) THEN
       CALL Fatal('List_NewEntryChunk', &
           'Could not allocate entry chunk of size: '//I2S(n))
     END IF

     Chunk % NextIndex = 1
     Chunk % Next => EntryChunks
     EntryChunks => Chunk
!-------------------------------------------------------------------------------
   END SUBROUTINE List_NewEntryChunk
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Return a single entry to the pool. The caller must have unlinked it from its
!> row already, and must not read the entry afterwards.
!-------------------------------------------------------------------------------
   SUBROUTINE List_ReturnEntry( Entry )
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: Entry

     Entry % Next => FreeEntries
     FreeEntries => Entry
     EntriesInUse = EntriesInUse - 1
!-------------------------------------------------------------------------------
   END SUBROUTINE List_ReturnEntry
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Return a whole chain of "cnt" entries, from Head to Tail, in one splice.
!-------------------------------------------------------------------------------
   SUBROUTINE List_ReturnEntryChain( Head, Tail, cnt )
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntry_t), POINTER :: Head, Tail
     INTEGER(KIND=8) :: cnt

     IF( .NOT. ASSOCIATED( Head ) ) RETURN

     Tail % Next => FreeEntries
     FreeEntries => Head
     EntriesInUse = EntriesInUse - cnt
!-------------------------------------------------------------------------------
   END SUBROUTINE List_ReturnEntryChain
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> If no list matrix holds an entry any more then every chunk is free, so give
!> them back to the system. The newest and largest chunk is kept so that the
!> next matrix does not have to start from scratch.
!-------------------------------------------------------------------------------
   SUBROUTINE List_ShrinkEntryPool()
!-------------------------------------------------------------------------------
     TYPE(ListMatrixEntryPool_t), POINTER :: Chunk, NextChunk

     IF( EntriesInUse /= 0 ) RETURN
     IF( .NOT. ASSOCIATED( EntryChunks ) ) RETURN

     FreeEntries => NULL()

     Chunk => EntryChunks % Next
     DO WHILE( ASSOCIATED( Chunk ) )
       NextChunk => Chunk % Next
       DEALLOCATE( Chunk % Entries )
       DEALLOCATE( Chunk )
       Chunk => NextChunk
     END DO

     EntryChunks % Next => NULL()
     EntryChunks % NextIndex = 1
!-------------------------------------------------------------------------------
   END SUBROUTINE List_ShrinkEntryPool
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
     CALL List_ReturnEntry(Clist)
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
       CALL List_ReturnEntry(Clist)
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
         CALL List_ReturnEntry(Clist)
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
     TYPE(ListMatrixEntry_t), POINTER :: Entry
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
     TYPE(ListMatrixEntry_t), POINTER :: Entry

     Entry => List_GetMatrixIndex(List,k1,k2)
!-------------------------------------------------------------------------------
   END SUBROUTINE List_AddMatrixIndex
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
   SUBROUTINE List_SetMatrixElement( List,k1,k2,Val,SetVal )
!-------------------------------------------------------------------------------
     TYPE(ListMatrix_t), POINTER :: List(:)
     INTEGER :: k1,k2
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
     TYPE(ListMatrixEntry_t), POINTER :: CList
     REAL(KIND=dp) :: Val
!-------------------------------------------------------------------------------

     Val = 0.0_dp

     IF ( .NOT. ASSOCIATED(List) ) RETURN
     IF ( k1>SIZE(List) ) RETURN
     Clist => List(k1) % Head

     DO WHILE( ASSOCIATED(CList) )
        IF ( Clist % INDEX >= k2 ) THEN
          IF ( Clist % INDEX == k2 ) Val = CList % Val
          RETURN
        END IF
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
     TYPE(ListMatrixEntry_t), POINTER :: p1, p2, prev2, Entry

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

     p1 => List(n1) % Head
     IF ( .NOT. ASSOCIATED(p1) ) THEN
       CALL Warn('List_MoveRow','Row not associated!')
       RETURN
     END IF

     ! Ensure row n2 is within the list bounds
     IF ( n2 > SIZE(List) ) THEN
       List => List_EnlargeMatrix(List, MAX(n2, SIZE(List)+LISTMATRIX_GROWTH))
     END IF

     ! Merge-walk: both rows are sorted so p2 only advances forward — O(d1+d2)
     NULLIFY(prev2)
     p2 => List(n2) % Head

     DO WHILE( ASSOCIATED(p1) )
       k2  = p1 % Index
       val = p1 % Val
       p1 % Val = d * val

       ! Advance p2 to the first entry with Index >= k2
       DO WHILE( ASSOCIATED(p2) )
         IF ( p2 % Index >= k2 ) EXIT
         prev2 => p2
         p2    => p2 % Next
       END DO

       IF ( ASSOCIATED(p2) .AND. p2 % Index == k2 ) THEN
         ! Entry already exists in row n2 — accumulate and advance past it
         p2 % Val = p2 % Val + c * val
         prev2 => p2
         p2    => p2 % Next
       ELSE
         ! Insert before p2 (or append if p2 is NULL)
         Entry => List_GetMatrixEntry(k2, p2)
         Entry % Val = c * val
         IF ( ASSOCIATED(prev2) ) THEN
           prev2 % Next => Entry
         ELSE
           List(n2) % Head => Entry
         END IF
         List(n2) % Degree = List(n2) % Degree + 1
         prev2 => Entry
       END IF

       p1 => p1 % Next
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
     INTEGER :: i, d1, d2
     INTEGER, ALLOCATABLE :: Ind1(:), Ind2(:)
     TYPE(ListMatrixEntry_t), POINTER :: CList1, CList2

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

     ! Take a copy of both index sets before either row is touched, then merge
     ! each into the other in one pass. Both rows are sorted so no sorting is
     ! needed. The outcome is the union of the two structures for both rows.
     d1 = 0
     DO WHILE( ASSOCIATED(CList1) )
       d1 = d1 + 1
       CList1 => CList1 % Next
     END DO

     d2 = 0
     DO WHILE( ASSOCIATED(CList2) )
       d2 = d2 + 1
       CList2 => CList2 % Next
     END DO

     ALLOCATE( Ind1(d1), Ind2(d2) )

     i = 0
     Clist1 => List(n1) % Head
     DO WHILE( ASSOCIATED(CList1) )
       i = i + 1
       Ind1(i) = Clist1 % Index
       CList1 => CList1 % Next
     END DO

     i = 0
     Clist2 => List(n2) % Head
     DO WHILE( ASSOCIATED(CList2) )
       i = i + 1
       Ind2(i) = Clist2 % Index
       CList2 => CList2 % Next
     END DO

     CALL List_AddMatrixIndexes( List,n2,d1,Ind1 )
     CALL List_AddMatrixIndexes( List,n1,d2,Ind2 )

     DEALLOCATE( Ind1, Ind2 )

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
     INTEGER :: i,k,nc,nu,Row
     INTEGER :: Cols(N*Dofs),UCols(N*Dofs),Slot(N*Dofs),Lcol(N*Dofs)
     REAL(KIND=dp) :: Vals(N*Dofs)

     ! Every row of the element shares the same column structure, so it is
     ! sorted and compressed just once. Each row is then glued with a single
     ! merge pass over the row list rather than one list walk per entry.
     CALL List_LocalColStructure( N,Dofs,Indexes,0,nc,nu,Cols,UCols,Slot,Lcol )
     IF( nu == 0 ) RETURN

     DO i=1,n
       IF (Indexes(i)<=0) CYCLE
       DO k=0,Dofs-1
         Row = Dofs*Indexes(i)-k

         CALL List_GatherLocalRow( LocalMatrix,Dofs*i-k,nc,nu,Slot,Lcol,Vals )
         CALL List_AddMatrixRow( A,Row,nu,UCols,Vals,SortedInput=.TRUE.)
       END DO
     END DO
   END SUBROUTINE List_GlueLocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Build the sorted and compressed column structure shared by all rows of a
!>    local (element) matrix. "Cols" holds the global column of each local
!>    column in order, "UCols(1:nu)" the sorted distinct ones, "Slot" maps each
!>    of the "nc" local columns to its place in UCols and "Lcol" back to the
!>    column index of the local matrix. Duplicate global columns, which arise
!>    if the index list repeats an entry, land in a single slot.
!------------------------------------------------------------------------------
   SUBROUTINE List_LocalColStructure( Ncol,ColDofs,ColInds,Col0,nc,nu, &
          Cols,UCols,Slot,Lcol )
!------------------------------------------------------------------------------
     INTEGER, INTENT(IN) :: Ncol,ColDofs,Col0
     INTEGER, INTENT(IN) :: ColInds(:)
     INTEGER, INTENT(OUT) :: nc,nu
     INTEGER, INTENT(OUT) :: Cols(:),UCols(:),Slot(:),Lcol(:)
!------------------------------------------------------------------------------
     INTEGER :: j,l,p

     nc = 0
     nu = 0
     DO j=1,Ncol
       IF( ColInds(j) <= 0 ) CYCLE
       DO l=0,ColDofs-1
         nc = nc + 1
         Cols(nc) = Col0 + ColDofs * ColInds(j) - l
         Lcol(nc) = ColDofs*j-l
       END DO
     END DO
     IF( nc == 0 ) RETURN

     UCols(1:nc) = Cols(1:nc)
     CALL Sort( nc, UCols )

     nu = 1
     DO p=2,nc
       IF( UCols(p) /= UCols(nu) ) THEN
         nu = nu + 1
         UCols(nu) = UCols(p)
       END IF
     END DO

     DO p=1,nc
       Slot(p) = SearchI( nu, UCols, Cols(p) )
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE List_LocalColStructure
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Gather one row of a local matrix into the compressed column order given
!>    by List_LocalColStructure().
!------------------------------------------------------------------------------
   SUBROUTINE List_GatherLocalRow( LocalMatrix,Lrow,nc,nu,Slot,Lcol,Vals )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(IN) :: LocalMatrix(:,:)
     INTEGER, INTENT(IN) :: Lrow,nc,nu
     INTEGER, INTENT(IN) :: Slot(:),Lcol(:)
     REAL(KIND=dp), INTENT(OUT) :: Vals(:)
!------------------------------------------------------------------------------
     INTEGER :: p

     Vals(1:nu) = 0.0_dp
     DO p=1,nc
       Vals(Slot(p)) = Vals(Slot(p)) + LocalMatrix(Lrow,Lcol(p))
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE List_GatherLocalRow
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
     INTEGER :: i,k,nc,nu,Row
     INTEGER :: Cols(Ncol*ColDofs),UCols(Ncol*ColDofs)
     INTEGER :: Slot(Ncol*ColDofs),Lcol(Ncol*ColDofs)
     REAL(KIND=dp) :: Vals(Ncol*ColDofs)

     ! As in List_GlueLocalMatrix: sort the shared column structure once, then
     ! add each row with a single merge pass.
     CALL List_LocalColStructure( Ncol,ColDofs,ColInds,Col0,nc,nu, &
             Cols,UCols,Slot,Lcol )
     IF( nu == 0 ) RETURN

     DO i=1,Nrow
       IF ( RowInds(i) <= 0 ) CYCLE
       DO k=0,RowDofs-1
         Row = Row0 + RowDofs * RowInds(i) - k

         CALL List_GatherLocalRow( LocalMatrix,RowDofs*i-k,nc,nu,Slot,Lcol,Vals )
         CALL List_AddMatrixRow( List,Row,nu,UCols,Vals,SortedInput=.TRUE.)
       END DO
     END DO
   END SUBROUTINE List_GlueLocalSubMatrix
!------------------------------------------------------------------------------

END MODULE ListMatrix

!> \} ElmerLib


