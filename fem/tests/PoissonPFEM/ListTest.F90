PROGRAM ListTest
    USE LocalTypes
    IMPLICIT NONE

    INTEGER, PARAMETER :: n=20

    TYPE(IntegerList_t) :: ilist
    INTEGER :: i

    CALL IntegerListInit(ilist, 8)

    ! Add numbers to list
    DO i=1, n
        CALL IntegerListAdd(ilist, i)
    END DO

    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)

    ! Remove every other element
    DO i=1,n,2
        CALL IntegerListDelete(ilist, i)
    END DO

    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)

    ! Try to remove first element
    i = IntegerListFind(ilist, 2)
    IF (i>0) CALL IntegerListDeleteAt(ilist, i)

    ! Try to remove last element
    i = IntegerListFind(ilist, n)
    IF (i>0) CALL IntegerListDeleteAt(ilist, i)

    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)

END PROGRAM ListTest
