!  Copyright (c) December 5 2006 - , CSC - IT Center for Science Ltd., Finland
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program (in file fem/GPL-2); if not, write to the
!  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
!  Boston, MA 02110-1301, USA.
!
!
! Wrappers for a bunch of C functions defined in binio.c to read/write
! from/to a binary stream in a portable (independent of endianness) manner.
!
! The OPTIONAL 'Status' argument, if PRESENT, will be set to > 0 for error (use
! StrErrorF to create meaningful error messages), 0 for success and (for the
! BinRead* procedures) -1 for end-of-file.  If .NOT.PRESENT(Status), then, on
! error, an error message will be written and the program terminated.
!
MODULE BinIO

    USE Kinds

    IMPLICIT NONE

    ! TODO: Use some kind of autoconf magic to set these to the C parameters
    ! SEEK_SET etc. That way we could pass the directly onwards to fseeko in the
    ! C code.
    INTEGER, PARAMETER :: BIN_SEEK_SET = 0, BIN_SEEK_CUR = 1, BIN_SEEK_END = 2

    INTERFACE
        FUNCTION BinFTell( Unit ) BIND(C, NAME="binftell_c")
            USE Kinds
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: Unit
            INTEGER(IntOff_k) :: BinFTell
        END FUNCTION BinFTell

        SUBROUTINE BinFSeek( Unit, Offset, Whence ) BIND(C, NAME="binfseek_c")
            USE Kinds
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT), INTENT(IN) :: Unit
            INTEGER(IntOff_k), INTENT(IN) :: Offset
            INTEGER, INTENT(IN) :: Whence
        END SUBROUTINE BinFSeek

        SUBROUTINE BinEndianess(e) BIND(C, NAME="binendianess_c")
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(KIND=C_CHAR), INTENT(OUT) :: e(*)
        END SUBROUTINE BinEndianess

        SUBROUTINE BinSetInputEndianess( Unit, e ) BIND(C, NAME="binsetinputendianess_c")
          USE, INTRINSIC :: ISO_C_BINDING
          INTEGER(C_INT) :: Unit
          CHARACTER(KIND=C_CHAR) :: e(*)
        END SUBROUTINE BinSetInputEndianess

        SUBROUTINE BinOpen_C(Unit,File,FileLen,Action,Status) BIND(C, NAME='binopen_c')
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(KIND=C_CHAR) :: File(*), Action(*)
          INTEGER(C_INT) :: Unit,FileLen,Status
        END SUBROUTINE BinOpen_C

        SUBROUTINE BinClose_C(unit, stat) BIND(C, NAME="binclose_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: unit, stat
        END SUBROUTINE BinClose_C

        SUBROUTINE BinWriteChar_C( Unit, c, Status ) BIND(C, NAME='binwritechar_c')
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(C_CHAR) :: c(*)
          INTEGER(C_INT) :: Unit, Status
        END SUBROUTINE BinWriteChar_C

        SUBROUTINE BinReadInt4_C(Unit, a, Status_) BIND(C, NAME="binreadint4_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: Unit, a, Status_
        END SUBROUTINE BinReadInt4_C

        SUBROUTINE BinWriteInt4_C( Unit, a, Status_ ) BIND(C, NAME="binwriteint4_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: Unit, Status_, a ! Note that this assumes LP64!
        END SUBROUTINE BinWriteInt4_c

        SUBROUTINE BinReadInt8_C(unit, a, stat) BIND(C, NAME="binreadint8_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: unit, stat
            INTEGER(C_INT_LEAST64_T) :: a
        END SUBROUTINE BinReadInt8_C

        SUBROUTINE BinWriteInt8_C(unit, a, stat) BIND(C, NAME="binwriteint8_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: unit, stat
            INTEGER(C_INT_LEAST64_T) :: a
        END SUBROUTINE BinWriteInt8_C

        SUBROUTINE BinReadDouble_C(unit, a, stat) BIND(C, NAME="binreaddouble_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: unit, stat
            REAL(C_DOUBLE) :: a
        END SUBROUTINE BinReadDouble_C

        SUBROUTINE BinWriteDouble_C(unit, a, stat) BIND(C, NAME="binwritedouble_c")
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT) :: unit, stat
            REAL(C_DOUBLE) :: a
        END SUBROUTINE BinWriteDouble_C

        SUBROUTINE BinWriteString_C( Unit, s, len, Status ) BIND(C, NAME='binwritestring_c')
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(C_CHAR) :: s(*)
          INTEGER(C_INT) :: Unit, len, Status
        END SUBROUTINE BinWriteString_C

        SUBROUTINE BinReadString_C( Unit, s, len, Status ) BIND(C, NAME='binreadstring_c')
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(C_CHAR) :: s(*)
          INTEGER(C_INT) :: Unit, len, Status
        END SUBROUTINE BinReadString_C


        SUBROUTINE StrErrorF_C( e,s,len) BIND(C, NAME="strerrorf_c")
          USE, INTRINSIC :: ISO_C_BINDING
          CHARACTER(C_CHAR) :: s(*)
          INTEGER(C_INT) :: e,len
        END SUBROUTINE StrErrorF_C
    END INTERFACE

    PRIVATE :: HandleStatus

CONTAINS

    SUBROUTINE HandleStatus( Status, Status_, MsgPrefix )
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER, INTENT(IN) :: Status_
        CHARACTER(LEN=*), INTENT(IN) :: MsgPrefix
        CHARACTER(LEN=100) :: Msg
        INTEGER, PARAMETER :: STDERR = 0

        IF ( PRESENT( Status ) ) THEN
            Status = Status_
        ELSE
            IF ( Status_ > 0 ) THEN
                CALL StrErrorF( Status_, Msg )
                WRITE( STDERR, * ) TRIM(MsgPrefix) // ": " // TRIM(Msg)
                STOP
            END IF
        END IF
    END SUBROUTINE HandleStatus


    SUBROUTINE BinOpen( Unit, File, Action, Status )
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER, INTENT(IN) :: Unit
        CHARACTER(*), INTENT(IN) :: File
        CHARACTER(*), INTENT(IN) :: Action ! "write", "append" or "read"
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinOpen_C( Unit, TRIM(File)//C_NULL_CHAR, LEN(TRIM(File))+1, Action, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Can't open file " &
                                            // TRIM(File) )
    END SUBROUTINE BinOpen


    SUBROUTINE BinClose( Unit, Status )
        INTEGER, INTENT(IN) :: Unit
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinClose_C( Unit, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Can't close file" )
    END SUBROUTINE BinClose


    SUBROUTINE BinWriteInt4( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        INTEGER(Int4_k), INTENT(IN) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinWriteInt4_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error writing Int4" )
    END SUBROUTINE BinWriteInt4


    SUBROUTINE BinReadInt4( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        INTEGER(Int4_k), INTENT(OUT) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinReadInt4_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error reading Int4" )
    END SUBROUTINE BinReadInt4


    SUBROUTINE BinWriteInt8( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        INTEGER(Int8_k), INTENT(IN) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinWriteInt8_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error writing Int8" )
    END SUBROUTINE BinWriteInt8


    SUBROUTINE BinReadInt8( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        INTEGER(Int8_k), INTENT(OUT) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinReadInt8_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error reading Int8" )
    END SUBROUTINE BinReadInt8


    SUBROUTINE BinWriteDouble( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        DOUBLE PRECISION, INTENT(IN) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinWriteDouble_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error writing Double" )
    END SUBROUTINE BinWriteDouble


    SUBROUTINE BinReadDouble( Unit, a, Status )
        INTEGER, INTENT(IN) :: Unit
        DOUBLE PRECISION, INTENT(OUT) :: a
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinReadDouble_C( Unit, a, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error reading Double" )
    END SUBROUTINE BinReadDouble


    ! Write a CHARACTER(1).  (Note: to read a CHARACTER(1), with no '\0' at the
    ! end, just use BinReadString with a CHARACTER(1) as argument.
    SUBROUTINE BinWriteChar( UNIT, c, Status )
        INTEGER, INTENT(IN) :: Unit
        CHARACTER, INTENT(IN) :: c
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinWriteChar_C( Unit, c, Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error writing char" )
    END SUBROUTINE BinWriteChar


    ! Write 's' to file pointed to by 'unit', and append a '\0'.
    SUBROUTINE BinWriteString( UNIT, s, Status )
        INTEGER, INTENT(IN) :: Unit
        CHARACTER(*), INTENT(IN) :: s
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinWriteString_C( Unit, s, LEN(s), Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error writing string" )
    END SUBROUTINE BinWriteString

    ! Read bytes up to next NULL byte (but no more than len(s)) into 's'
    ! from file pointed to by 'unit', and pad with spaces. The NULL byte
    ! will be removed from the stream, but not added to 's'.
    SUBROUTINE BinReadString( Unit, s, Status )
        INTEGER, INTENT(IN) :: Unit
        CHARACTER(*), INTENT(OUT) :: s
        INTEGER, OPTIONAL, INTENT(OUT) :: Status
        INTEGER :: Status_

        CALL BinReadString_C( Unit, s, LEN(s), Status_ )
        CALL HandleStatus( Status, Status_, "BINIO: Error reading string" )
    END SUBROUTINE BinReadString

    ! Return a string representation of the error code 'e' (as returned in a
    ! status variable) in 's'. (Actually just a wrapper to the C function
    ! strerror().)
    SUBROUTINE StrErrorF( e, s )
        INTEGER, INTENT(IN) :: e
        CHARACTER(*), INTENT(OUT) :: s

        CALL StrErrorF_C( e, s, LEN(s) )
    END SUBROUTINE StrErrorF

END MODULE BinIO
