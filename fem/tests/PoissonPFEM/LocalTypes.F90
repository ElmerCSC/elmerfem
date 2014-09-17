MODULE LocalTypes
	USE Types
    IMPLICIT NONE

	TYPE ElementCache_t
		TYPE(ElementType_t), POINTER :: Type => NULL()
		INTEGER :: nc = 0
		REAL(Kind=dp), POINTER :: U(:)=>NULL(),V(:)=>NULL(),W(:)=>NULL()
	END TYPE ElementCache_t

    TYPE VertexMap_t
        INTEGER, ALLOCATABLE :: vlist(:)
    END TYPE VertexMap_t

    TYPE IntegerList_t
        INTEGER :: nelem = 0
        INTEGER, ALLOCATABLE :: entries(:)
    END TYPE IntegerList_t

    TYPE IntegerHashBucket_t
        TYPE(IntegerList_t), ALLOCATABLE :: list
        TYPE(IntegerList_t), ALLOCATABLE :: entryloc
    END TYPE IntegerHashBucket_t

    TYPE IntegerHashSet_t
        TYPE(IntegerHashBucket_t), ALLOCATABLE :: set(:)
        TYPE(IntegerList_t), ALLOCATABLE :: entries
    END TYPE IntegerHashSet_t

    INTEGER, PARAMETER :: INTEGERLIST_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_CHAIN_DEFAULT_SIZE = 8
    REAL(KIND=dp), PARAMETER :: INTEGERHASHSET_FILLRATIO = 0.75

CONTAINS

    SUBROUTINE IntegerHashSetInit(iset, isize)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN), OPTIONAL :: isize

        INTEGER :: i, n, cn, allocstat

        n = INTEGERHASHSET_DEFAULT_SIZE
        IF (PRESENT(isize)) n = isize

        IF (ALLOCATED(iset % set)) DEALLOCATE(iset % set)

        ALLOCATE(iset % set(FLOOR(n/INTEGERHASHSET_FILLRATIO)), &
                 iset % entries, STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('IntegerHashSetInit', 'Memory allocation error!')

        CALL IntegerListInit(iset % entries, n)
    END SUBROUTINE IntegerHashSetInit

    SUBROUTINE IntegerHashSetAdd(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: hkey, lind, allocstat

        ! Check for overfill and rehash if necessary
        ! TODO!

        ! Add entry to Hash table
        hkey = IntegerHashSetHashFunction(iset, key)
        IF (.NOT. ALLOCATED(iset % set(hkey) % list)) THEN
            ! Add new chain and entry
            ALLOCATE(iset % set(hkey) % list, iset % set(hkey) % entryloc, STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerHashSetAdd', 'Memory allocation error!')

            CALL IntegerListInit(iset % set(hkey) % list, INTEGERHASHSET_CHAIN_DEFAULT_SIZE)
            CALL IntegerListInit(iset % set(hkey) % entryloc, INTEGERHASHSET_CHAIN_DEFAULT_SIZE)

            CALL IntegerListAdd(iset % set(hkey) % list, key)
            CALL IntegerListAdd(iset % entries, key)
            CALL IntegerListAdd(iset % set(hkey) % entryloc, iset % entries % nelem)
        ELSE
            ! Old chain, just add entry
            ASSOCIATE(chain => iset % set(hkey))
                lind = IntegerListFind(chain % list, key)

                ! A set can only contain one reference to each item
                IF (lind < 0) THEN
                    CALL IntegerListAdd(chain % list, key)
                    CALL IntegerListAdd(iset % entries, key)
                    CALL IntegerListAdd(chain % entryloc, iset % entries % nelem)
                END IF
            END ASSOCIATE
        END IF
    END SUBROUTINE IntegerHashSetAdd

    SUBROUTINE IntegerHashSetDelete(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: hkey, lind, ekey

        hkey = IntegerHashSetHashFunction(iset, key)
        ASSOCIATE(chain => iset % set(hkey))

            IF (ALLOCATED(chain % list)) THEN
                lind = IntegerListFind(chain % list, key)

                IF (lind > 0) THEN
                    ! Delete values from chain and store index of entry
                    CALL IntegerListDeleteAt(chain % list, lind)
                    ekey = IntegerListAt(chain % entryloc, lind)
                    CALL IntegerListDeleteAt(chain % entryloc, lind)
                    ! Delete entry from entry list
                    CALL IntegerListDeleteAt(iset % entries, ekey)
                END IF
            END IF
        END ASSOCIATE
    END SUBROUTINE IntegerHashSetDelete

    SUBROUTINE IntegerHashSetDeleteAll(iset)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER :: ent, hkey

        DO ent=1, iset % entries % nelem
            hkey = IntegerHashSetHashFunction(iset, iset % entries % entries(ent))

            ! Elements of chain may not be allocated
            ! if they have been removed before
            IF (ALLOCATED(iset % set(hkey) % list)) THEN
                CALL IntegerListDeleteAll(iset % set(hkey) % list)
                CALL IntegerListDeleteAll(iset % set(hkey) % entryloc)
                DEALLOCATE(iset % set(hkey) % list, iset % set(hkey) % entryloc)
            END IF
        END DO

        CALL IntegerListDeleteAll(iset % entries)
    END SUBROUTINE IntegerHashSetDeleteAll

    FUNCTION IntegerHashSetFind(iset, key) RESULT(found)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t), INTENT(IN) :: iset
        INTEGER, INTENT(IN) :: key
        LOGICAL :: found

        INTEGER :: hkey, ind

        found = .FALSE.
        hkey = IntegerHashSetHashFunction(iset, key)
        IF (ALLOCATED(iset % set(hkey) % list)) THEN
            ind = IntegerListFind(iset % set(hkey) % list, key)
            IF (ind > 0) found = .TRUE.
        END IF
    END FUNCTION IntegerHashSetFind

    FUNCTION IntegerHashSetHashFunction(iset, key) RESULT(hkey)
        IMPLICIT NONE
        TYPE(IntegerHashSet_t), INTENT(IN) :: iset
        INTEGER, INTENT(IN) :: key

        INTEGER :: hkey

        ! NIY
        hkey = 0
    END FUNCTION IntegerHashSetHashFunction

    SUBROUTINE IntegerListInit(ilist, isize)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, OPTIONAL, INTENT(IN) :: isize

        INTEGER :: n, allocstat

        n = INTEGERLIST_DEFAULT_SIZE
        IF (PRESENT(isize)) n = isize

        IF (ALLOCATED(ilist % entries)) DEALLOCATE(ilist % entries)
        ALLOCATE(ilist % entries(n), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('IntegerListInit', 'Memory allocation error!')

        ilist % nelem = 0
    END SUBROUTINE IntegerListInit

    FUNCTION IntegerListAt(ilist, ind) RESULT(key)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: ind
        INTEGER :: key

        key = 0
        IF (ind < 0 .AND. ind > SIZE(ilist % entries)) RETURN
        key = ilist % entries(ind)
    END FUNCTION IntegerListAt

    SUBROUTINE IntegerListAdd(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER, ALLOCATABLE :: elementsnew(:)
        INTEGER :: n, allocstat

        n = size(ilist % entries)
        IF (ilist % nelem + 1 > n) THEN
            ! Reallocate list structure with double the size
            ALLOCATE(elementsnew(2*n), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerListAdd', 'Memory allocation error!')

            ! Copy elements, deallocate old element vector and move allocation
            elementsnew(1:n)=ilist % entries(1:n)
            DEALLOCATE(ilist % entries)
            CALL MOVE_ALLOC(elementsnew, ilist % entries)
        END IF

        ! Add new element
        ilist % nelem = ilist % nelem + 1
        ilist % entries(ilist % nelem)=item
    END SUBROUTINE IntegerListAdd

    SUBROUTINE IntegerListDelete(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER :: i, ind

        ind = IntegerListFind(ilist, item)
        CALL IntegerListDeleteAt(ilist, ind)
    END SUBROUTINE IntegerListDelete

    SUBROUTINE IntegerListDeleteAt(ilist, ind)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: ind

        INTEGER :: i

        IF (ind < 0 .AND. ind > SIZE(ilist % entries)) RETURN

        DO i=ind, ilist % nelem-1
            ilist % entries(i) = ilist % entries(i+1)
        END DO
        ilist % entries(ilist % nelem) = 0
        ilist % nelem = ilist % nelem - 1
    END SUBROUTINE IntegerListDeleteAt

    SUBROUTINE IntegerListDeleteAll(ilist)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist

        INTEGER :: allocstat

        ilist % entries = 0
        ilist % nelem = 0
    END SUBROUTINE IntegerListDeleteAll

    FUNCTION IntegerListFind(ilist, item) RESULT(ind)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item
        INTEGER :: ind

        INTEGER :: i

        ind = -1
        DO i=1,ilist % nelem
            IF (ilist % entries(i) == item) THEN
                ind = i
                EXIT
            END IF
        END DO
    END FUNCTION IntegerListFind

END MODULE LocalTypes
