PROGRAM SetTest
    USE LocalTypes
    IMPLICIT NONE

    INTEGER, PARAMETER :: n=16, isetsize=8

    TYPE(IntegerHashSet_t) :: iset
    INTEGER :: i, hval
    LOGICAL :: hashInRange

    CALL IntegerHashSetInit(iset, isetsize)

    ! Hash function values for keys
    WRITE (*,'(A,I0)') 'Hash table size is ', SIZE(iset % set)
    WRITE (*,'(A,I0)') 'Hash keys from 1 to ', 4*isetsize
    hashInRange = .TRUE.
    DO i=1,4*isetsize
        hval = IntegerHashSetHashFunction(iset, i)
        IF (hval < 1 .OR. hval > SIZE(iset % set)) THEN
            WRITE (*,'(A,I0,A,I0)') 'Out of range: key=',i,' , value=', hval
            hashInRange = .FALSE.
        ELSE
            WRITE (*,'(A,I0,A,I0)') 'key=',i,' , value=', hval
        END IF
    END DO
    IF (hashInRange) THEN
        WRITE (*,*) 'Hash function passes correctness check.'
    ELSE
        WRITE (*,*) 'ERROR: Incorrect hash function. Aborting test.'
        STOP
    END IF

    WRITE (*,'(A,I0,A,I0,A,I0)') 'Add numbers from ', 1, ' to ', n, ' by ', 4
    ! Add elements to set
    DO i=1, n, 4
        CALL IntegerHashSetAdd(iset, i)
    END DO
    WRITE (*,'(A,I0,A,I0,A,I0)') 'Find numbers from ', 1, ' to ', n
    DO i=1,n
        IF (IntegerHashSetFind(iset, i)) WRITE (*,'(A,I0)') 'Found ', i
    END DO
    CALL CheckIntegerHashSetConsistency(iset)
    ! CALL WriteIntegerHashSetContents(iset)
    
    WRITE (*,'(A,I0,A,I0,A,I0)') 'Insert numbers from ', 1, ' to ', n
    DO i=1,n
        CALL IntegerHashSetAdd(iset, i)
    END DO
    WRITE (*,'(A,I0,A,I0)') 'Find numbers from ', 1, ' to ', n
    DO i=1,n
        IF (IntegerHashSetFind(iset, i)) WRITE (*,'(A,I0)') 'Found ', i
    END DO
    ! CALL WriteIntegerHashSetContents(iset)
    CALL CheckIntegerHashSetConsistency(iset)

    WRITE (*,'(A,I0,A,I0,A,I0)') 'Delete numbers from ', 1, ' to ', n, ' by ', 2
    DO i=1,n,2
        CALL IntegerHashSetDelete(iset, i)
    END DO
    WRITE (*,'(A,I0,A,I0,A,I0)') 'Find numbers from ', 1, ' to ', n
    DO i=1,n
        IF (IntegerHashSetFind(iset, i)) WRITE (*,'(A,I0)') 'Found ', i
    END DO
    CALL CheckIntegerHashSetConsistency(iset)
    ! CALL WriteIntegerHashSetContents(iset)

    WRITE (*,'(A,I0,A,I0,A,I0)') 'Delete numbers ', 1,' and ', 5
    CALL IntegerHashSetDelete(iset, 1)
    CALL IntegerHashSetDelete(iset, 5)
    IF (IntegerHashSetFind(iset, 1)) WRITE (*,'(A,I0)') 'Error, found ', 1
    IF (IntegerHashSetFind(iset, 5)) WRITE (*,'(A,I0)') 'Error ,found ', 5
    CALL CheckIntegerHashSetConsistency(iset)

    WRITE (*,'(A)') 'Remove all elements'
    CALL IntegerHashSetDeleteAll(iset)
    WRITE (*,'(A,I0,A,I0,A,I0)') 'Find numbers from ', 1, ' to ', n
    DO i=1,n
        IF (IntegerHashSetFind(iset, i)) WRITE (*,'(A,I0)') 'Found ', i
    END DO
    CALL CheckIntegerHashSetConsistency(iset)

    CONTAINS

        SUBROUTINE CheckIntegerHashSetConsistency(iset)
            IMPLICIT NONE
            
            TYPE(IntegerHashSet_t) :: iset

            INTEGER :: i, j, ent, lind
            LOGICAL :: consistent
            
            consistent = .TRUE.
            ! Check that mappings from buckets to elements is consistent
            DO i=1,SIZE(iset % set)
                IF (ALLOCATED(iset % set(i) % list)) THEN
                    DO j=1,iset % set(i) % list % nelem
                        ent = iset % set(i) % list % entries(j)
                        lind = IntegerListFind(iset % entries, ent)
                        
                        ! Each element must be present in entry list
                        IF (lind <= 0) THEN
                            WRITE (*,*) 'set=', i, ', entry=', j, ' val=', ent, &
                                  ' not found from entry list although in set'
                            consistent = .FALSE.
                        END IF
                    END DO
                END IF
            END DO

            ! Check that mapping from elements to buckets is consistent
            DO i=1,iset % entries % nelem
                ent = iset % entries% entries(i)
                IF (.NOT. IntegerHashSetFind(iset, ent)) THEN
                    WRITE (*,*) 'entries(i)=', ent, ' not found from set although in list'
                    consistent = .FALSE.
                END IF
            END DO

            IF (consistent) THEN 
                WRITE (*,*) 'Map consistent'
            ELSE
                WRITE (*,*) 'ERROR: Map not consistent'
            END IF

        END SUBROUTINE CheckIntegerHashSetConsistency

END PROGRAM SetTest
