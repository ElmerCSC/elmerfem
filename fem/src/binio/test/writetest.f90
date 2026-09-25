PROGRAM WriteTest

    USE BinIO

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: fu = 10
    INTEGER :: a, stat, i, iounit
    REAL(dp) :: v(3)
    CHARACTER(1) :: E

    a = 5
    v = (/ -7.13_DP, 0.01_DP, HUGE(v) /)

    CALL BinOpen( fu, "writetest.out", "write", stat )
    IF ( stat /= 0 ) STOP 1

    CALL BinEndianess( E )
    CALL BinWriteChar( fu,E,stat )
    IF ( stat /= 0 ) STOP 2

    CALL BinWriteInt4( fu, a, stat )
    IF ( stat /= 0) STOP 3

    CALL BinWriteString( fu, "Hello World!", stat )
    IF ( stat /= 0 ) STOP 4

    DO i = 1, SIZE( v )
        CALL BinWriteDouble( fu, v(i), stat )
        IF ( stat /= 0 ) STOP 5
    END DO

    CALL BinClose( fu, stat )
    IF ( stat /= 0 ) STOP 6

    OPEN( NEWUNIT=iounit, file="writetest.out", position="append" )
    WRITE( iounit, '(A)' ) ''
    WRITE( iounit, '(A)' ) "Humhum"
    CLOSE( iounit )

    CALL BinOpen( fu, "writetest.out", "append", stat )
    IF ( stat /= 0 ) STOP 7

    DO i = 1, SIZE( v )
        CALL BinWriteDouble( fu, v(i), stat )
        IF ( stat /= 0 ) STOP 8
    END DO

    CALL BinClose( fu, stat )
    IF ( stat /= 0 ) STOP 9

END PROGRAM WriteTest
