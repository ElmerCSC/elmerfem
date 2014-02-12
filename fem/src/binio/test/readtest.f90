PROGRAM ReadTest
    
    USE BinIO

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = selected_real_kind(15)
    INTEGER, PARAMETER :: fu = 10

    IF( TestFile( "writetest.checkL" ) .AND. TestFile( "writetest.checkB" ) )THEN
        PRINT *, 'OK'
    ELSE
        PRINT *, 'DAMN!'
    END IF

CONTAINS

    LOGICAL FUNCTION TestFile( FName )
        CHARACTER(*) :: FName
        INTEGER :: a, stat, i
        REAL(dp) :: v(3), v2(3)
        CHARACTER(50) :: s, s2
        CHARACTER(256) :: Error
        CHARACTER :: E

        CALL BinOpen( fu, FName, "read", stat )
        IF ( stat /= 0 ) STOP 1

        CALL BinReadString( fu,E,stat )
        IF ( stat /= 0 ) STOP "1.5"
        CALL BinSetInputEndianess( fu,E )

        CALL BinReadInt4( fu, a, stat )
        IF ( stat /= 0 ) STOP 2

        CALL BinReadString( fu, s, stat )
        IF ( stat /= 0 ) STOP 3

        DO i = 1, SIZE( v )
            CALL BinReadDouble( fu, v(i), stat )
            IF ( stat /= 0 ) STOP 4
        END DO

        CALL BinReadString( fu, s2(1:8), stat )
        IF ( stat /= 0 ) STOP 5

        DO i = 1, SIZE( v2 )
            CALL BinReadDouble( fu, v2(i), stat )
            IF ( stat /= 0 ) STOP 6
        END DO

        ! At this point, we should have reached the end.
        CALL BinReadInt4( fu, i, stat )
        IF ( stat >= 0 ) STOP 7

        CALL BinClose( fu, stat )
        IF ( stat /= 0 ) STOP 8

        IF ( a == 5 .AND. TRIM(s) == "Hello World!"  &
                    .AND. ALL( v == (/ -7.13_dp, 0.01_dp, HUGE(v) /) ) &
                    .AND. ALL( v2 == (/ -7.13_dp, 0.01_dp, HUGE(v2) /) ) &
                    .AND. s2(2:7) == "Humhum" .AND. IACHAR(s2(1:1)) == 10 &
                                              .AND. IACHAR(s2(8:8)) == 10 ) then
            TestFile = .TRUE.
        ELSE
            !PRINT *, FName
            !PRINT *, a
            !PRINT *, "'", TRIM(s), "'"
            !PRINT *, v
            !PRINT *, "'", TRIM(s2(2:7)), "'"
            !PRINT *, IACHAR(s2(1:1)), IACHAR(s2(8:8))
            !PRINT *, v2
            TestFile = .FALSE.
        END IF
    END FUNCTION TestFile

END PROGRAM ReadTest
