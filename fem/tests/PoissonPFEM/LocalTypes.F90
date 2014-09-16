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

    TYPE IntegerHashSet_t
        INTEGER :: nelem = 0
        TYPE(IntegerList_t), ALLOCATABLE :: set(:)
        INTEGER, ALLOCATABLE :: entries(:)
    END TYPE IntegerHashSet_t

    INTEGER, PARAMETER :: INTEGERLIST_DEFAULT_SIZE = 64
    INTEGER, PARAMETER :: INTEGERHASHSET_DEFAULT_SIZE = 64
    REAL(KIND=dp), PARAMETER :: INTEGERHASHSET_FILLRATIO = 0.75

CONTAINS

    SUBROUTINE IntegerHashSetInit(iset, isize)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN), OPTIONAL :: isize

        INTEGER :: i, n, allocstat

        n = INTEGERHASHSET_DEFAULT_SIZE
        IF (PRESENT(isize)) n = isize

        IF (ALLOCATED(iset % set)) THEN
            DEALLOCATE(iset % set)
            DEALLOCATE(iset % entries)
        END IF
        ALLOCATE(iset % set(FLOOR(n/INTEGERHASHSET_FILLRATIO)), &
                 iset % entries(n), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('IntegerHashSetInit', 'Memory allocation error!')

        DO i=1,n
            CALL IntegerListInit(iset % set(i), 8)
        END DO

        iset % entries = 0
        iset % nelem = 0
    END SUBROUTINE IntegerHashSetInit

    SUBROUTINE IntegerHashSetAdd(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        ! NIY
    END SUBROUTINE IntegerHashSetAdd

    SUBROUTINE IntegerHashSetDelete(iset, key)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset
        INTEGER, INTENT(IN) :: key

        ! NIY
    END SUBROUTINE IntegerHashSetDelete

    SUBROUTINE IntegerHashSetDeleteAll(iset)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t) :: iset

        ! NIY
    END SUBROUTINE IntegerHashSetDeleteAll

    SUBROUTINE IntegerHashSetFind(iset, key, found)
        IMPLICIT NONE

        TYPE(IntegerHashSet_t), INTENT(IN) :: iset
        INTEGER, INTENT(IN) :: key
        LOGICAL, INTENT(OUT) :: found

        ! NIY
    END SUBROUTINE IntegerHashSetFind

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

    SUBROUTINE IntegerListAdd(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER, ALLOCATABLE :: elementsnew(:)
        INTEGER :: n, allocstat

        n = size(ilist % elements)
        IF (ilist % nelem + 1 > n) THEN
            ! Reallocate list structure with double the size
            ALLOCATE(elementsnew(2*n), STAT=allocstat)
            IF (allocstat /= 0) CALL Fatal('IntegerListAdd', 'Memory allocation error!')

            ! Copy elements, deallocate old element vector and move allocation
            elementsnew(1:n)=ilist % elements(1:n)
            DEALLOCATE(ilist % elements)
            CALL MOVE_ALLOC(elementsnew, ilist % elements)
        END IF

        ! Add new element
        ilist % nelem = ilist % nelem + 1
        ilist % elements(ilist % nelem)=item
    END SUBROUTINE IntegerListAdd

    SUBROUTINE IntegerListDelete(ilist, item)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item

        INTEGER :: i, ind

        CALL IntegerListFind(ilist, item, ind)
        IF (ind > 0) THEN
            DO i=ind, ilist % nelem-1
                ilist % elements(i) = ilist % elements(i+1)
            END DO
            ilist % elements(ilist % nelem) = 0
            ilist % nelem = ilist % nelem - 1
        END IF
    END SUBROUTINE IntegerListDelete

    SUBROUTINE IntegerListDeleteAll(ilist)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist

        INTEGER :: allocstat

        ilist % entries = 0
        ilist % nelem = 0
    END SUBROUTINE IntegerListDeleteAll

    SUBROUTINE IntegerListFind(ilist, item, ind)
        IMPLICIT NONE

        TYPE(IntegerList_t) :: ilist
        INTEGER, INTENT(IN) :: item
        INTEGER, INTENT(OUT) :: ind

        ind = -1
        DO i=1,ilist % nelem
            IF (ilist % elements(i) == item) THEN
                ind = i
                EXIT
            END IF
        END DO
    END SUBROUTINE IntegerListFind

END MODULE LocalTypes
