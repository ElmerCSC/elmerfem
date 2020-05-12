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
!
!/******************************************************************************
! *
! *  Authors: Mikko Byckling
! *
! *  Original Date: 25 September 2017
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{


MODULE ListMatrixArray
  USE Messages
  USE Types
  USE GeneralUtils, ONLY : I2S
  
  IMPLICIT NONE
CONTAINS
  
  !-------------------------------------------------------------------------------
  !> Allocates an empty array list matrix.
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_Allocate(ListMatrixArray, N, PoolSize, Atomic)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    INTEGER,INTENT(IN) :: N
    INTEGER, OPTIONAL :: PoolSize
    LOGICAL, OPTIONAL :: Atomic
    
    INTEGER :: i,istat, nthr, TID, psize
    LOGICAL :: InitLocks
    
    psize = 1024
    IF (PRESENT(PoolSize)) psize = PoolSize

    InitLocks = .FALSE.
    IF (PRESENT(Atomic)) InitLocks = Atomic
    
    ! Allocate ListMatrix and associated pools
    nthr = 1
    !$ nthr = omp_get_max_threads()
    ALLOCATE( ListMatrixArray % Rows(n), &
              ListMatrixArray % Pool(nthr), STAT=istat )
    IF( istat /= 0 ) THEN
      CALL Fatal('ListMatrixArray_AllocateMatrix',&
                 'Allocation error for ListMatrix of size: '//TRIM(I2S(n)))
    END IF
    IF (InitLocks) CALL ListMatrixArray_InitializeAtomic(ListMatrixArray)
    
    !$OMP PARALLEL &
    !$OMP SHARED(ListMatrixArray, N, psize) &
    !$OMP PRIVATE(i, TID) DEFAULT(NONE)
    
    TID = 1
    !$ TID = omp_get_thread_num()+1

    CALL ListMatrixPool_Initialize(ListMatrixArray % Pool(TID), psize)
    
    !$OMP DO
    DO i=1,N
      ListMatrixArray % Rows(i) % Head => NULL()
      ListMatrixArray % Rows(i) % Level = 0
      ListMatrixArray % Rows(i) % Degree = 0
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
  END SUBROUTINE ListMatrixArray_Allocate
 
  !-------------------------------------------------------------------------------
  !> Free an array list matrix.
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_Free( ListMatrixArray )
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray

    TYPE(ListMatrixEntryPool_t), POINTER :: p, p1
    INTEGER :: N,TID
    
    N = SIZE(ListMatrixArray % Pool)
    !$OMP PARALLEL &
    !$OMP SHARED(ListMatrixArray, N) &
    !$OMP PRIVATE(p, p1, TID) DEFAULT(NONE)

    !$OMP DO
    DO TID=1,N
       CALL ListMatrixPool_Free(ListMatrixArray % Pool(TID))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    CALL ListMatrixArray_FreeAtomic(ListMatrixArray)
    
    DEALLOCATE(ListMatrixArray % Rows, ListMatrixArray % Pool)
  END SUBROUTINE ListMatrixArray_Free

  SUBROUTINE ListMatrixArray_InitializeAtomic(ListMatrixArray)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray

    INTEGER :: i, N, istat
    
#ifdef _OPENMP
    N = SIZE(ListMatrixArray % Rows)
    
    ALLOCATE( ListMatrixArray % RowLocks(n), STAT=istat )
    IF( istat /= 0 ) THEN
      CALL Fatal('ListMatrixArray_InitializeAtomic',&
            'Allocation error for ListMatrix row locks of size: '//TRIM(I2S(n)))
    END IF
      
    !$OMP PARALLEL DO &
    !$OMP SHARED(ListMatrixArray,N) &
    !$OMP PRIVATE(i) DEFAULT(NONE)
    DO i=1,N
      CALL omp_init_lock(ListMatrixArray % RowLocks(i))
    END DO
    !$OMP END PARALLEL DO
#endif
  END SUBROUTINE ListMatrixArray_InitializeAtomic

  SUBROUTINE ListMatrixArray_FreeAtomic(ListMatrixArray)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray

    INTEGER :: i, N

#ifdef _OPENMP
    IF (ALLOCATED(ListMatrixArray % RowLocks)) THEN
      N = SIZE(ListMatrixArray % RowLocks)
      
      !$OMP PARALLEL DO &
      !$OMP SHARED(ListMatrixArray,N) &
      !$OMP PRIVATE(i) DEFAULT(NONE)
      DO i=1,N
        CALL omp_destroy_lock(ListMatrixArray % RowLocks(i))
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE(ListMatrixArray % RowLocks)
    END IF
#endif
  END SUBROUTINE ListMatrixArray_FreeAtomic

  SUBROUTINE ListMatrixArray_LockRow(ListMatrixArray, row, Atomic)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    INTEGER, INTENT(IN) :: row
    LOGICAL, OPTIONAL :: Atomic

#ifdef _OPENMP
    IF (PRESENT(ATOMIC)) THEN
      IF (Atomic) CALL omp_set_lock(ListMatrixArray % RowLocks(row))
    END IF
#endif
  END SUBROUTINE ListMatrixArray_LockRow

  SUBROUTINE ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    INTEGER, INTENT(IN) :: row
    LOGICAL, OPTIONAL :: Atomic

#ifdef _OPENMP
    IF (PRESENT(ATOMIC)) THEN
      IF (Atomic) CALL omp_unset_lock(ListMatrixArray % RowLocks(row))
    END IF
#endif
  END SUBROUTINE ListMatrixArray_UnlockRow
  
  !-------------------------------------------------------------------------------
  !> Transfer sparsity pattern of the array list matrix format to a graph format,
  !> used in most places of the code. 
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_ToGraph( ListMatrixArray, Graph)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    TYPE(Graph_t) :: Graph
    
    ! TODO
    CALL Fatal('ListMatrixArray_ToGraph','Not implemented yet!')
  END SUBROUTINE ListMatrixArray_ToGraph

  !-------------------------------------------------------------------------------
  !> Transfer the flexible list matrix to the more efficient CRS matrix that is 
  !> used in most places of the code. The matrix structure can accommodate both forms.
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_ToCRSMatrix( ListMatrixArray, CRSMatrix )
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    TYPE(Matrix_t) :: CRSMatrix
    
    ! TODO
    CALL Fatal('ListMatrixArray_ToCRSMatrix','Not implemented yet!')
  END SUBROUTINE ListMatrixArray_ToCRSMatrix

  SUBROUTINE ListMatrixArray_FromCRSMatrix( ListMatrixArray, CRSMatrix )
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    TYPE(Matrix_t) :: CRSMatrix
    
    ! TODO
    CALL Fatal('ListMatrixArray_FromCRSMatrix','Not implemented yet!')
  END SUBROUTINE ListMatrixArray_FromCRSMatrix

  !-------------------------------------------------------------------------------
  !> Add index (row,col) to the matrix sparsity structure 
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_AddEntry(ListMatrixArray, row, col, val, Atomic)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    INTEGER, INTENT(IN) :: row, col
    REAL(KIND=dp), OPTIONAL :: val
    LOGICAL, OPTIONAL :: Atomic

    TYPE(ListMatrixEntry_t), POINTER :: CEntryPtr, PEntryPtr, NEntryPtr
    INTEGER :: TID

    TID = 1
    !$ TID = omp_get_thread_num() + 1

    CALL ListMatrixArray_LockRow(ListMatrixArray, row, Atomic)
    
    CEntryPtr => ListMatrixArray % Rows(row) % Head
    IF (.NOT. ASSOCIATED(CEntryPtr)) THEN
       ! Empty matrix row, add entry and return
       ListMatrixArray % Rows(row) % Head => &
            ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, NULL())
       ListMatrixArray % Rows(row) % Degree = 1
       CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
       RETURN
    ELSE IF (CEntryPtr % Index == col) THEN
       ! Do not add duplicates, nothing to do!
       CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
       RETURN
    ELSE IF (CEntryPtr % Index > col) THEN
       ! Add a new entry to the Head of list
       ListMatrixArray % Rows(row) % Head => &
            ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, CEntryPtr)
       ListMatrixArray % Rows(row) % Degree = & 
            ListMatrixArray % Rows(row) % Degree + 1
       CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
       RETURN
    END IF
    
    ! Search a correct place for the element
    PEntryPtr => CEntryPtr
    CEntryPtr => CEntryPtr % Next
    
    DO WHILE( ASSOCIATED(CEntryPtr) )
       ! Do not add duplicates
       IF (CEntryPtr % Index == col) THEN
         CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
         RETURN
       END IF
       ! Place found, exit search loop
       IF (CEntryPtr % Index > col) EXIT
       
       PEntryPtr => CEntryPtr
       CEntryPtr => CEntryPtr % Next
    END DO

    ! Add entry to the correct place in the list
    PEntryPtr % Next => ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, CEntryPtr)
    ListMatrixArray % Rows(row) % Degree = & 
          ListMatrixArray % Rows(row) % Degree + 1
    CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
  END SUBROUTINE ListMatrixArray_AddEntry

  !-------------------------------------------------------------------------------
  !> Add indexes on a single row to the matrix sparsity structure.
  !-------------------------------------------------------------------------------
  SUBROUTINE ListMatrixArray_AddEntries(ListMatrixArray, row, nentry, Indexes, Perm, Atomic)
    IMPLICIT NONE
    TYPE(ListMatrixArray_t) :: ListMatrixArray
    INTEGER, INTENT(IN) :: row, nentry
    INTEGER, INTENT(IN) :: Indexes(nentry), Perm(nentry)
    LOGICAL, OPTIONAL :: Atomic
    
    TYPE(ListMatrixEntry_t), POINTER :: CEntryPtr, PEntryPtr, NEntryPtr
    INTEGER :: TID, centry, sentry, rentry, col
        
    TID = 1
    !$ TID = omp_get_thread_num() + 1

    CALL ListMatrixArray_LockRow(ListMatrixArray, row, Atomic)
    
    CEntryPtr => ListMatrixArray % Rows(row) % Head
    sentry = 1
    col = Indexes(Perm(1))
    IF (.NOT. ASSOCIATED(CEntryPtr)) THEN
       ! Empty matrix row, add entry and continue
       CEntryPtr => ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, NULL())
       ListMatrixArray % Rows(row) % Head => CEntryPtr
       ListMatrixArray % Rows(row) % Degree = 1
       sentry = 2
    ELSE IF (CEntryPtr % Index == col) THEN
       ! Do not add duplicates, continue with next element!
       sentry = 2
    ELSE IF (CEntryPtr % Index > col) THEN
       ! Add a new entry to the Head of list
       NEntryPtr => ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, CEntryPtr)
       CEntryPtr => NEntryPtr
       ListMatrixArray % Rows(row) % Head => CEntryPtr
       ListMatrixArray % Rows(row) % Degree = & 
            ListMatrixArray % Rows(row) % Degree + 1
       sentry = 2
    END IF
    
    ! Search a correct place for the element
    PEntryPtr => CEntryPtr
    CEntryPtr => CEntryPtr % Next

    DO centry=sentry,nentry
       col=Indexes(Perm(centry))

       ! Find a correct place to add index to
       DO WHILE( ASSOCIATED(CEntryPtr) )
         IF (CEntryPtr % Index >= col) EXIT
         PEntryPtr => CEntryPtr
         CEntryPtr => PEntryPtr % Next
       END DO
       
       IF (ASSOCIATED(CEntryPtr)) THEN
         ! Do not add duplicates
         IF (CEntryPtr % Index /= col) THEN
           ! Create new element between PEntryPtr and CEntryPtr
           NEntryPtr => ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, CEntryPtr)
           PEntryPtr % Next => NEntryPtr
           ListMatrixArray % Rows(row) % Degree = & 
                ListMatrixArray % Rows(row) % Degree + 1

           ! Advance to next element in list
           PEntryPtr => NEntryPtr
           CEntryPtr => NEntryPtr % Next
         ELSE
           ! Advance to next element in list
           PEntryPtr => CEntryPtr
           CEntryPtr => CEntryPtr % Next
         END IF
       ELSE
         ! List matrix row contains no more entries
         EXIT
       END IF
     END DO

     ! Add rest of the entries in Indexes to list matrix row (if any)
     DO rentry=centry,nentry
       col=Indexes(Perm(rentry))
       NEntryPtr => ListMatrixPool_GetListEntry(ListMatrixArray % Pool(TID), col, NULL())
       PEntryPtr % Next => NEntryPtr
       PEntryPtr => NEntryPtr
       ListMatrixArray % Rows(row) % Degree = & 
            ListMatrixArray % Rows(row) % Degree + 1
     END DO

     CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
   END SUBROUTINE ListMatrixArray_AddEntries

   !-------------------------------------------------------------------------------
   !> Delete entry (row,col) from the matrix sparsity structure 
   !-------------------------------------------------------------------------------
   SUBROUTINE ListMatrixArray_DeleteEntry(ListMatrixArray, row, col, Atomic)
     IMPLICIT NONE
     TYPE(ListMatrixArray_t) :: ListMatrixArray
     INTEGER, INTENT(IN) :: row, col
     LOGICAL, OPTIONAL :: Atomic
     
     TYPE(ListMatrixEntry_t), POINTER :: CEntryPtr, PEntryPtr
     INTEGER :: TID
     LOGICAL :: NotFound
     
     TID = 1
     !$ TID = omp_get_thread_num() + 1

     CALL ListMatrixArray_LockRow(ListMatrixArray, row, Atomic)
     
     ! Search for element from the list
     PEntryPtr => NULL()
     CEntryPtr => ListMatrixArray % Rows(row) % Head     
     DO WHILE( ASSOCIATED(CEntryPtr) )
       IF (CEntryPtr % Index >= col) EXIT
       
       PEntryPtr => CEntryPtr
       CEntryPtr => CEntryPtr % Next
     END DO

     IF (ASSOCIATED(CEntryPtr)) THEN
       IF (CEntryPtr % Index == col) THEN
         ! Element found, delete it from list and add it to
         ! the pooled list of deleted entries
         IF (ASSOCIATED(PEntryPtr)) THEN
           PEntryPtr % Next => CEntryPtr % Next
         ELSE
           ListMatrixArray % Rows(row) % Head => CEntryPtr % Next
         END IF
         CALL ListMatrixPool_AddDeletedEntry(ListMatrixArray % Pool(TID), CEntryPtr)
         
         ListMatrixArray % Rows(row) % Degree = &
           MAX(ListMatrixArray % Rows(row) % Degree - 1, 0)
       END IF
     END IF
     
     CALL ListMatrixArray_UnlockRow(ListMatrixArray, row, Atomic)
   END SUBROUTINE ListMatrixArray_DeleteEntry
   
   !-------------------------------------------------------------------------------
   !> ListMatrixPool support routines
   !-------------------------------------------------------------------------------
   SUBROUTINE ListMatrixPool_Initialize(Pool, PoolSize)
     IMPLICIT NONE
     
     TYPE(ListMatrixPool_t) :: Pool
     INTEGER, INTENT(IN) :: PoolSize
     
     Pool % EntryPool => NULL()
     Pool % Deleted => NULL()
     Pool % PoolSize = PoolSize
     CALL ListMatrixPool_EnLarge(Pool)
   END SUBROUTINE ListMatrixPool_Initialize

   SUBROUTINE ListMatrixPool_Enlarge(Pool)
     IMPLICIT NONE

     TYPE(ListMatrixPool_t) :: Pool

     TYPE(ListMatrixEntryPool_t), POINTER :: EntryPool

     INTEGER :: astat
     
     ALLOCATE(EntryPool, STAT=astat)
     IF (astat == 0) ALLOCATE(EntryPool % Entries(Pool % PoolSize), STAT=astat)
     IF (astat /= 0) THEN
        CALL Fatal('ListMatrixPool_Enlarge','Pool allocation failed')
     END IF

     EntryPool % NextIndex = 1

     EntryPool % Next => Pool % EntryPool
     Pool % EntryPool => EntryPool
   END SUBROUTINE ListMatrixPool_Enlarge

   SUBROUTINE ListMatrixPool_Free(Pool)
     IMPLICIT NONE

     TYPE(ListMatrixPool_t) :: Pool
     TYPE(ListMatrixEntryPool_t), POINTER :: EntryPool, EntryPoolNext

     EntryPool => Pool % EntryPool
     DO WHILE (ASSOCIATED(EntryPool))
        EntryPoolNext => EntryPool % Next
        DEALLOCATE(EntryPool % Entries)
        DEALLOCATE(EntryPool)
        EntryPool => EntryPoolNext
     END DO
   END SUBROUTINE ListMatrixPool_Free

   FUNCTION ListMatrixPool_GetListEntry(Pool, ind, Next) RESULT(ListEntry)
     IMPLICIT NONE

     TYPE(ListMatrixPool_t) :: Pool
     INTEGER, INTENT(IN) :: ind
     TYPE(ListMatrixEntry_t), POINTER :: Next

     TYPE(ListMatrixEntry_t), POINTER :: ListEntry

     ! Check if deleted entries are available
     IF (ASSOCIATED(Pool % Deleted)) THEN
        ListEntry => Pool % Deleted
        Pool % Deleted => ListEntry % Next
     ELSE
        ! No deleted entries available, allocate a new pool if necessary
        IF (Pool % PoolSize < Pool % EntryPool % NextIndex) THEN
           CALL ListMatrixPool_Enlarge(Pool)
        END IF
        
        ! Get next element from pool
        ListEntry => Pool % EntryPool % Entries(Pool % EntryPool % NextIndex)
        Pool % EntryPool % NextIndex = Pool % EntryPool % NextIndex + 1
     END IF
     
     ListEntry % Index = ind
     ListEntry % Next => Next
   END FUNCTION ListMatrixPool_GetListEntry

   SUBROUTINE ListMatrixPool_AddDeletedEntry(Pool, DEntry)
     IMPLICIT NONE

     TYPE(ListMatrixPool_t) :: Pool
     TYPE(ListMatrixEntry_t), POINTER :: DEntry

     ! Add new deleted entry to the head of the deleted entries list
     DEntry % Next => Pool % Deleted
     Pool % Deleted => DEntry
   END SUBROUTINE ListMatrixPool_AddDeletedEntry

END MODULE ListMatrixArray

!> \} ElmerLib
