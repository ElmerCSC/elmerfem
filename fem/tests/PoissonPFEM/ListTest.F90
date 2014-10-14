PROGRAM ListTest
    USE LocalTypes
    IMPLICIT NONE

    INTEGER, PARAMETER :: n=20

    TYPE(IntegerList_t) :: ilist
    INTEGER :: i

    CALL IntegerListInit(ilist, 8)

    WRITE (*,'(A,I0,A,I0)') 'Add numbers from ', 1, ' to ', n
    ! Add numbers to list
    DO i=1, n
        CALL IntegerListAdd(ilist, i)
    END DO

    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)
    
    WRITE (*,'(A,I0,A,I0,A,I0)') 'Remove numbers from ', 1, ' to ', n, ' by ', 2
    ! Remove every other element
    DO i=1,n,2
        CALL IntegerListDelete(ilist, i)
    END DO

    WRITE (*,'(A,I0)') 'Find elements'
    DO i=1,n
        WRITE (*,'(A,I0,A,I0)') 'Find ', i, ' got ', IntegerListFind(ilist, i)
    END DO
 
    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)

    WRITE (*,'(A,I0)') 'Find nonexisting elements'
    IF (IntegerListFind(ilist, -10) >= 0) THEN
        WRITE (*,*) 'Error: expect element -10 not in list' 
    END IF
    IF (IntegerListAt(ilist, -10) /= 0) THEN
        WRITE (*,*) 'Error: expect element -10 not in list' 
    END IF
   

    WRITE (*,'(A)') 'Remove first element'
    ! Try to remove first element
    i = IntegerListFind(ilist, 2)
    IF (i>0) CALL IntegerListDeleteAt(ilist, i)

    WRITE (*,'(A)') 'Remove middle element'
    ! Try to remove first element
    i = IntegerListFind(ilist, 10)
    IF (i>0) CALL IntegerListDeleteAt(ilist, i)

    WRITE (*,'(A)') 'Remove last element'
    ! Try to remove last element
    i = IntegerListFind(ilist, n)
    IF (i>0) CALL IntegerListDeleteAt(ilist, i)

    WRITE (*,'(A)') 'Remove nonexisting elements'
    CALL IntegerListDelete(ilist, -10)
    CALL IntegerListDeleteAt(ilist, -10)

    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)
    
    WRITE (*,*) 'Remove all elements'
    CALL IntegerListDeleteAll(ilist)
    
    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)

    WRITE (*,*) 'Re-insert all elements'
    ! Add numbers to list
    DO i=1, n
        CALL IntegerListAdd(ilist, i)
    END DO
    
    WRITE (*,*) 'nelem=', ilist % nelem
    WRITE (*,*) ilist % entries(1:ilist % nelem)
    
    CALL IntegerListDeleteAll(ilist)
END PROGRAM ListTest
