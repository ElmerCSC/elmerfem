!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! *****************************************************************************/

! #ifndef USE_ISO_C_BINDINGS
#include "../config.h"
! #endif

!> \ingroup ElmerLib
!> \}

!-----------------------------------------------------------------------------
!>  Miscallenous utilities.
!-----------------------------------------------------------------------------
MODULE GeneralUtils

USE Types
#ifdef USE_ISO_C_BINDINGS
USE LoadMod
#endif

#ifdef HAVE_LUA
USE, INTRINSIC :: ISO_C_BINDING
#endif

IMPLICIT NONE

INTERFACE AllocateVector
  MODULE PROCEDURE AllocateRealVector, AllocateIntegerVector, &
                   AllocateComplexVector, AllocateLogicalVector, &
                   AllocateElementVector
END INTERFACE

INTERFACE AllocateArray
  MODULE PROCEDURE AllocateRealArray, AllocateIntegerArray, &
                   AllocateComplexArray, AllocateLogicalArray
END INTERFACE

INTERFACE ComponentName
   MODULE PROCEDURE ComponentNameStr, ComponentNameVar
END INTERFACE

    REAL(KIND=dp), PRIVATE :: AdvanceTime1, AdvanceTime2

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE StartAdvanceOutput( SolverName, OutputType )
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: SolverName, OutputType
!------------------------------------------------------------------------------
#ifndef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: RealTime
#endif
     AdvanceTime1 = RealTime()
     AdvanceTime2 = RealTime()
     CALL Info( SolverName, OutputType, Level=5 )
!------------------------------------------------------------------------------
  END SUBROUTINE StartAdvanceOutput
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AdvanceOutput(t,n,dot_t,percent_t)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: t,n
     REAL(KIND=dp), OPTIONAL :: dot_t,percent_t
!------------------------------------------------------------------------------
     INTEGER :: i
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: d_t, p_t
#else
     REAL(KIND=dp) :: RealTime, d_t, p_t
#endif
!------------------------------------------------------------------------------
     d_t = 1._dp
     p_t = 20._dp
     IF ( PRESENT(dot_t) ) d_t = dot_t
     IF ( PRESENT(percent_t) ) p_t = percent_t

     IF ( RealTime() - AdvanceTime1 > d_t ) THEN
       CALL Info( '', '.', Level=5, noAdvance=.TRUE. )

       IF ( RealTime() - AdvanceTime2 > p_t ) THEN
         i = NINT(t*100.0/n)
         WRITE(Message, '(i3,a)' ) i, '%'
         CALL Info( '', Message, Level=5 )
         AdvanceTime2 = RealTime()
       END IF
       AdvanceTime1 = RealTime()
     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AdvanceOutput
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   PURE FUNCTION lentrim(str) RESULT(n)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*), INTENT(IN) :: str
     INTEGER :: n
     DO n=LEN(str),1,-1
       IF ( str(n:n) /= ' ' ) EXIT
     END DO
!------------------------------------------------------------------------------
   END FUNCTION lentrim
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compare equality of start of s1 to (in most uses string literal) s2.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PURE FUNCTION SEQL(s1,s2) RESULT(L)
!------------------------------------------------------------------------------
    LOGICAL :: L
    CHARACTER(LEN=*), INTENT(IN) :: s1,s2
!------------------------------------------------------------------------------
    INTEGER :: n
!------------------------------------------------------------------------------
    L = .FALSE.
    n = LEN(s2)
    IF(LEN(s1) < n) RETURN
    IF (s1(1:n)==s2) L=.TRUE.
!------------------------------------------------------------------------------
  END FUNCTION SEQL
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Converts integer to string. Handy when writing output with integer data.
!------------------------------------------------------------------------------
  PURE FUNCTION i2s(ival) RESULT(str)
!------------------------------------------------------------------------------
    INTEGER, INTENT(in) :: ival
    CHARACTER(LEN=12) :: str
!------------------------------------------------------------------------------
    INTEGER :: i,j,n,m,t,v
    CHARACTER, PARAMETER :: DIGITS(0:9)=['0','1','2','3','4','5','6','7','8','9']
!------------------------------------------------------------------------------
     str = ' '

     IF ( ival >= 0 ) THEN
       j=0
       v=ival
     ELSE
       str(1:1)='-'
       j=1
       v=-ival
     END IF

     IF (v<10) THEN
       str(j+1:j+1)=DIGITS(v)
     ELSE
       n=2
       m=10
       DO WHILE(10*m<=v)
         n=n+1
         m=m*10
       END DO

       DO i=j+1,j+n
         t = v / m
         str(i:i) = DIGITS(t)
         v = v - t*m
         m = m / 10
       END DO
     END IF
!------------------------------------------------------------------------------
  END FUNCTION i2s
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Converts string of length n to an integer number. 
!------------------------------------------------------------------------------
  PURE FUNCTION s2i(str,n) RESULT(ival)
!------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    INTEGER :: ival
    CHARACTER(LEN=n), INTENT(IN) :: str
!------------------------------------------------------------------------------
    LOGICAL :: neg
    INTEGER :: j,k    
    INTEGER, PARAMETER :: ic0 = ICHAR('0')
    
    neg = str(1:1)=='-'
    k=1
    IF ( neg ) k=2
    
    ival = 0
    DO j=k,n
      ival = 10*ival + ICHAR(str(j:j)) - ic0
    END DO
    IF(neg) ival=-ival
  END FUNCTION s2i
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> Converts a string into a number of integer numbers
!> It is assumed that the integers may also be separated by 
!> the given separator. 
!------------------------------------------------------------------------------
  FUNCTION str2ints(str,ints,sep) RESULT(n)
!------------------------------------------------------------------------------
    INTEGER, INTENT(out) :: ints(:)
    CHARACTER(LEN=*), INTENT(in) :: str
    CHARACTER, OPTIONAL, INTENT(in) :: sep

    INTEGER :: i,k,l,m,n,ic, icsep
    INTEGER, PARAMETER :: ic0 = ICHAR('0'), ic9 = ICHAR('9'), icm = ICHAR('-'), &
        ics = ICHAR(' ')

    IF( PRESENT( sep ) ) THEN
      icsep = ICHAR(sep)
    ELSE
      icsep = ics
    END IF


    k = LEN_TRIM(str)
    l = 1
    n = 0
    DO WHILE(l<=k.AND.n<SIZE(ints))
      DO WHILE(l<=k)
        ic = ICHAR(str(l:l))
        IF( ic == ics .OR. ic == icsep ) THEN
          CONTINUE
        ELSE
          EXIT
        END IF
        l=l+1
      END DO
      IF(l>k) EXIT
      IF(.NOT.(ic==icm .OR. ic>=ic0 .AND. ic<=ic9)) EXIT

      m = l+1
      DO WHILE(m<=k)
        ic = ICHAR(str(m:m))
        IF(ic<ic0 .OR. ic>ic9) EXIT
        m=m+1
      END DO

      n = n + 1
      ints(n) = s2i(str(l:m-1),m-l)
      l = m
    END DO
  END FUNCTION str2ints
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE SystemCommand( cmd ) 
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: cmd
    CALL SystemC( TRIM(cmd) // CHAR(0) )
!------------------------------------------------------------------------------
  END SUBROUTINE SystemCommand
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  PURE FUNCTION FileNameQualified(file) RESULT(L)
!------------------------------------------------------------------------------
    LOGICAL :: L
    CHARACTER(*), INTENT(IN) :: file
!------------------------------------------------------------------------------
    L = INDEX(file,':')>0 .OR. file(1:1)=='/' .OR. file(1:1)==Backslash
!------------------------------------------------------------------------------
  END FUNCTION FileNameQualified
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   PURE FUNCTION LittleEndian() RESULT(L)
!------------------------------------------------------------------------------
     LOGICAL :: L
!------------------------------------------------------------------------------
     INTEGER(1) :: s(2)
     INTEGER(2), PARAMETER :: t = 256*7+8
!------------------------------------------------------------------------------
     s = TRANSFER(t,s)
     L=s(1)==8
!------------------------------------------------------------------------------
   END FUNCTION LittleEndian
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION FormatDate() RESULT( date )
!------------------------------------------------------------------------------
    CHARACTER( LEN=20 ) :: date
    INTEGER :: dates(8)

    CALL DATE_AND_TIME( VALUES=dates )
    WRITE( date, &
     '(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)' ) &
                dates(1),dates(2),dates(3),dates(5),dates(6),dates(7)
!------------------------------------------------------------------------------
  END FUNCTION FormatDate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sort an array of integer values. 
!------------------------------------------------------------------------------
   PURE SUBROUTINE Sort( n,a )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in)  :: n
     INTEGER, INTENT(inout) :: a(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,ra
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )
        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
        ELSE
         ra = a(ir)
         a(ir) = a(1)
         ir = ir - 1
         IF ( ir == 1 ) THEN
           a(1) = ra
           RETURN
         END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF

          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            i = j
            j =  j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Sort
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Sort an integer array a, together with an another integer array.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortI( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: a(:),b(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,ra,rb
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )
        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
         ra = a(ir)
         rb = b(ir)
         a(ir) = a(1)
         b(ir) = b(1)
         ir = ir - 1
         IF ( ir == 1 ) THEN
           a(1) = ra
           b(1) = rb
           RETURN
         END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
             IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j =  j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortI
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sort an index array, and change the order of an real array accordingly.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortF( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: a(:)
     REAL(KIND=dp), INTENT(inout) :: b(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,ra
     REAL(KIND=dp) :: rb
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortF
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sort an real array, and change the order of an index array accordingly.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortD( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: b(:)
     REAL(KIND=dp), INTENT(inout) :: a(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,rb
     REAL(KIND=dp) :: ra
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortD
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sort an complex array, and organize an index table accordingly.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortC( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: b(:)
     COMPLEX(KIND=dp), INTENT(inout) :: a(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,rb
     COMPLEX(KIND=dp) :: ra
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )
        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir ) THEN
             IF ( ABS(a(j))<ABS(a(j+1)) ) j = j+1
          END IF
          IF ( ABS(ra)<ABS(a(j)) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Order real components in b in a decreasing order and return the new order
!> of indexes in a.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortR( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: a(:)
     REAL(KIND=dp), INTENT(inout) :: b(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,ra
     REAL(KIND=dp) :: rb
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
             IF ( b(j) > b(j+1) ) j = j+1
          END IF
          IF ( rb > b(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Search an integer value in an ordered array. 
!------------------------------------------------------------------------------
   PURE FUNCTION SearchI( N,Array,Val ) RESULT ( Idx )
!------------------------------------------------------------------------------
    INTEGER, INTENT(in) :: N,Val,Array(:)
!------------------------------------------------------------------------------
    INTEGER :: Lower, Upper,Lou,Idx
!------------------------------------------------------------------------------

    Idx = 0 
    Upper = N
    Lower = 1

    ! Handle the special case

    IF ( Upper == 0 ) RETURN

    DO WHILE( .TRUE. )
      IF ( Array(Lower) == Val) THEN
         Idx = Lower
         EXIT
      ELSE IF ( Array(Upper) == Val ) THEN
         Idx = Upper
         EXIT
      END IF

      IF ( (Upper-Lower)>1 ) THEN
        Lou = ISHFT((Upper + Lower), -1)
        IF ( Array(Lou) < Val ) THEN
          Lower = Lou
        ELSE
          Upper = Lou
        END IF
      ELSE
        EXIT
      END IF
    END DO
    
    RETURN

!------------------------------------------------------------------------------
  END FUNCTION SearchI
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Search a real value in an ordered array. 
!------------------------------------------------------------------------------
  PURE FUNCTION SearchR( N,Array,Val ) RESULT ( Idx )
!------------------------------------------------------------------------------

    INTEGER, INTENT(in) :: N
    INTEGER :: Idx
    REAL(KIND=dP), INTENT(in) :: Val,Array(:)
!------------------------------------------------------------------------------
    INTEGER :: Lower, Upper,Lou
!------------------------------------------------------------------------------

    Idx = 0
    Upper = N
    Lower = 1

    ! Handle the special case
    IF ( Upper == 0 ) RETURN

    DO WHILE( .TRUE. )
      IF ( ABS( Array(Lower) - Val) < TINY(Val)  ) THEN
        Idx = Lower
        EXIT
      ELSE IF ( ABS( Array(Upper) - Val ) < TINY(Val) ) THEN
        Idx = Upper
        EXIT
      END IF

      IF ( (Upper-Lower) > 1 ) THEN
        Lou = ISHFT((Upper + Lower), -1)
        IF ( Array(Lou) < Val ) THEN
          Lower = Lou
        ELSE
          Upper = Lou
        END IF
      ELSE
        EXIT
      END IF
    END DO
    
    RETURN

!------------------------------------------------------------------------------
  END FUNCTION SearchR
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE OpenIncludeFile( Unit, FileName, IncludePath )
!------------------------------------------------------------------------------
    INTEGER :: Unit
    CHARACTER(LEN=*) :: FileName, IncludePath
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,k0,k1,l
    CHARACTER(LEN=1024) :: name, TmpName
!------------------------------------------------------------------------------

    i = 1
    name = FileName
    DO WHILE( name(i:i) == ' ' .OR. name(i:i)=='"')
      i = i + 1
    END DO
    j = LEN_TRIM(name)
    IF ( name(j:j) == '"' ) j=j-1
    name = TRIM(name(i:j))

    IF ( INDEX(name,':') == 0 .AND. name(1:1) /= '/' .AND. &
              name(1:1) /= Backslash ) THEN
       k0 = 1
       DO WHILE( IncludePath(k0:k0) == '"' )
         k0 = k0+1
       END DO
       k1 = INDEX( IncludePath, ';' )

       DO WHILE( k1 >= k0 )
         DO k = k1-1,k0,-1
           IF ( IncludePath(k:k) /= ' ' .AND. IncludePath(k:k)/='"' ) EXIT
         END DO 
         IF ( IncludePath(k:k) == '"' ) k=k-1
         IF ( k >= k0 ) THEN
           WRITE( tmpName,'(a,a,a)' ) IncludePath(k0:k), '/', TRIM(name)
           OPEN( Unit, FILE=TRIM(tmpName), STATUS='OLD',ERR=10 )
           RETURN
         END IF
10       CONTINUE
         k0 = k1+1
         k1 = INDEX( IncludePath(k0:), ';' ) + k0 - 1
       END DO

       IF ( LEN_TRIM(IncludePath(k0:))>0 ) THEN
         k1 = INDEX( IncludePath(k0:), '"' ) + k0 - 2
         IF ( k1 < k0 ) k1=LEN_TRIM(IncludePath)
         tmpName = TRIM(IncludePath(k0:k1)) //  '/' // TRIM(name)
         OPEN( Unit, FILE=TRIM(TmpName), STATUS='OLD',ERR=20 )
         RETURN
       END IF

20     CONTINUE
       OPEN( Unit, FILE=TRIM(name), STATUS='OLD' )
    ELSE
       OPEN( Unit, FILE=TRIM(name), STATUS='OLD' )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE OpenIncludeFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Read a (logical) line from FORTRAN device Unit and remove leading, trailing,
!>  and multiple blanks between words. Also convert uppercase characters to
!>  lowercase.The logical line can continue the several physical lines by adding
!>  the backslash (\) mark at the end of a physical line. 
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ReadAndTrim( Unit,str,echo,literal ) RESULT(l)
!------------------------------------------------------------------------------
!******************************************************************************
!
!
!  ARGUMENTS:
!
!     INTEGER :: Unit
!       INPUT: Fortran unit number to read from
!
!     CHARACTER :: str
!       OUTPUT: The string read from the file
!
!  FUNCTION RESULT:
!      LOGICAL :: l
!        Success of the read operation
!
!******************************************************************************
     INTEGER, PARAMETER :: MAXLEN = 16384

     INTEGER :: Unit
     CHARACTER(LEN=:), ALLOCATABLE :: str

     LOGICAL, OPTIONAL :: Echo, literal

     LOGICAL :: l

     CHARACTER(LEN=:), ALLOCATABLE :: temp
     CHARACTER(LEN=12) :: tmpstr
     CHARACTER(LEN=MAXLEN) :: readstr = ' ', copystr = ' ', matcstr=' ' , IncludePath=' '

     LOGICAL :: InsideQuotes, OpenSection=.FALSE.
     INTEGER :: i,j,k,m,ValueStarts=0,inlen,ninlen,outlen,IncludeUnit=28,IncludeUnitBase=28

     CHARACTER(LEN=MAX_NAME_LEN) :: Prefix = '  '

     INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A'),Tab=9
     CHARACTER(LEN=MAXLEN) :: tmatcstr, tcmdstr
     INTEGER :: tninlen

     SAVE ReadStr, ValueStarts, Prefix, OpenSection

     IF ( PRESENT(literal) ) literal=.FALSE.
     l = .TRUE.

     IF(.NOT.ALLOCATED(str)) ALLOCATE(CHARACTER(512)::str)
     outlen = LEN(str)

     IF ( ValueStarts==0 .AND. OpenSection ) THEN
       str = 'end'
       ValueStarts = 0
       OpenSection = .FALSE.
       RETURN
     END IF

     IF ( ValueStarts == 0 ) THEN
        tmpstr = ' '
        DO WHILE( .TRUE. )
          IF ( IncludeUnit < IncludeUnitBase ) THEN
            READ( IncludeUnit,'(A)',END=1,ERR=1 ) readstr
            GO TO 2
1           CLOSE(IncludeUnit)
            IncludeUnit = IncludeUnit+1
            READ( Unit,'(A)',END=10,ERR=10 ) readstr
2           CONTINUE
          ELSE
            READ( Unit,'(A)',END=10,ERR=10 ) readstr
          END IF

          readstr = ADJUSTL(readstr)

          DO k=1,12
            j = ICHAR(readstr(k:k))
            IF ( j >= A .AND. j<= Z ) THEN
              Tmpstr(k:k) = CHAR(j+U2L)
            ELSE
              tmpstr(k:k) = readstr(k:k)
            END IF
          END DO

          IF ( SEQL(Tmpstr, 'include path') ) THEN
            k = LEN_TRIM(readstr)
            IncludePath(1:k-13) = readstr(14:k)
            Tmpstr = ''
          ELSE
            EXIT
          END IF
        END DO

        IF ( SEQL(tmpstr, 'include ') ) THEN
          IncludeUnit = IncludeUnit-1
          CALL OpenIncludeFile( IncludeUnit, TRIM(readstr(9:)), IncludePath )
          READ( IncludeUnit,'(A)',END=3,ERR=3 ) readstr
          GO TO 4
3         CLOSE(IncludeUnit)
          IncludeUnit = IncludeUnit+1
          READ( Unit,'(A)',END=10,ERR=10 ) readstr
4         CONTINUE
        END IF
        ninlen = LEN_TRIM(readstr)
     ELSE
        inlen = LEN_TRIM(readstr)
        ninlen = inlen-ValueStarts+1
        IF ( Prefix == ' ' ) THEN
           readstr = readstr(ValueStarts:inlen)
        ELSE IF ( Prefix == '::' ) THEN
           readstr = readstr(ValueStarts:inlen)
           OpenSection = .TRUE.
           Prefix = ' '
        ELSE
           DO i=ValueStarts,inlen
              IF ( readstr(i:i) ==  ')' ) THEN
                 readstr(i:i) = ' '
                 EXIT
              ELSE IF ( readstr(i:i) == ',' ) THEN
                 readstr(i:i) = ' '
              END IF
           END DO
           ninlen = ninlen + LEN_TRIM(Prefix) + 1
           readstr = TRIM(Prefix) // ' ' // readstr(ValueStarts:inlen)
        END IF
     END IF

     ValueStarts = 0
     InsideQuotes  = .FALSE.

     i = INDEX( readstr(1:ninlen), '!' )
     IF ( i>0 ) ninlen=i-1

     i = 1
     inlen = ninlen
     DO WHILE( i <= inlen )
       IF ( readstr(i:i) == '"' ) InsideQuotes = .NOT.InsideQuotes
       IF ( .NOT. InsideQuotes .AND. readstr(i:i) == CHAR(92) .AND. i==inlen ) THEN
          readstr(i:i) = ' '
          IF ( IncludeUnit < IncludeUnitBase ) THEN
            READ( IncludeUnit,'(A)',END=10,ERR=10 ) readstr(i+1:MAXLEN)
          ELSE
            READ( Unit,'(A)',END=10,ERR=10 ) readstr(i+1:MAXLEN)
          END IF
          DO j=LEN(readstr),i+1,-1
             IF ( readstr(j:j) /= ' ' ) EXIT
          END DO
          inlen = inlen + j-i
       END IF
       i = i + 1
     END DO

#ifdef HAVE_LUA

     block 
       integer :: lstat
       character(kind=c_char, len=:), pointer :: lua_result
       integer :: result_len
       logical :: closed_region, first_bang
       closed_region = .false.
       first_bang = .true.
       i = INDEX( readstr(1:inlen), '#' )

       IF ( i>0 .AND. i<inlen ) THEN
         m = i
         copystr(i:inlen) = readstr(i:inlen)
         DO WHILE(i<=inlen)
           IF ( copystr(i:i) == '#' ) THEN
             DO j=i+1,inlen-1
               IF ( copystr(j:j) == '#' ) EXIT
             END DO
             ninlen = j - i

             ! Initialize variables for each copy of Lua interpreter separately

             !$OMP PARALLEL DEFAULT(NONE) &
             !$OMP SHARED(copystr, i, matcstr, ninlen, inlen, closed_region, first_bang, j) &
             !$OMP PRIVATE(tcmdstr, tninlen, lstat, result_len, lua_result) 

             tninlen = ninlen
             tcmdstr = copystr(i+1:inlen)

             IF(tcmdstr(tninlen:tninlen) == '#') then 
               closed_region = .TRUE.
             ELSE
               closed_region = .FALSE.
             END IF

             IF(closed_region) THEN
               lstat = lua_dostring( LuaState, &
                   'return tostring('// tcmdstr(1:tninlen-1) // ')'//c_null_char, 1)
             ELSE
               IF (i == 1 .and. first_bang .and. j == inlen) THEN  ! ' # <luacode>' case, dont do 'return tostring(..)'.
                                                                   ! Instead, just execute the line in the lua interpreter
                 lstat = lua_dostring( LuaState, tcmdstr(1:tninlen) // c_null_char, 1)
               ELSE ! 'abc = # <luacode>' case, oneliners only
                 lstat = lua_dostring( LuaState, &
                     'return tostring('// tcmdstr(1:tninlen) // ')'//c_null_char, 1)
               END IF
             END IF
             lua_result => lua_popstring(LuaState, result_len)

             !$OMP SINGLE 
             matcstr(1:result_len) = lua_result(1:result_len)
             ninlen = result_len
             !$OMP END SINGLE

             !$OMP END PARALLEL

             DO k=1,ninlen
               readstr(m:m) = matcstr(k:k)
               m = m + 1
             END DO
             i = j+1
           ELSE
             readstr(m:m) = copystr(i:i)
             i = i + 1
             m = m + 1
           END IF
           first_bang = .false.
         END DO
         IF ( m <= inlen ) readstr(m:inlen) = ' '
         inlen = m-1
       END IF
     end block
#endif

     i = INDEX( readstr(1:inlen), '$' )
     IF ( i>0 .AND. i<inlen ) THEN
       m = i
       copystr(i:inlen) = readstr(i:inlen)
       DO WHILE(i<=inlen)
         IF ( copystr(i:i) == '$' ) THEN
            DO j=i+1,inlen-1
              IF ( copystr(j:j) == '$' ) EXIT
            END DO
            ninlen = j - i

            ! Initialize variables for each copy of MATC separately

            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP SHARED(copystr, i, matcstr, ninlen, inlen) &
            !$OMP PRIVATE(tcmdstr, tmatcstr, tninlen)

            tninlen = ninlen
            tcmdstr = copystr(i+1:inlen)
            CALL MATC( tcmdstr, tmatcstr, tninlen )
            !$OMP BARRIER

            !$OMP SINGLE
            matcstr(1:tninlen) = tmatcstr(1:tninlen)
            ninlen = tninlen
            !$OMP END SINGLE

            !$OMP END PARALLEL

            DO k=1,ninlen
              readstr(m:m) = matcstr(k:k)
              m = m + 1
            END DO
            i = j+1
         ELSE
            readstr(m:m) = copystr(i:i)
            i = i + 1
            m = m + 1
         END IF
       END DO
       IF ( m <= inlen ) readstr(m:inlen) = ' '
       inlen = m-1
     END IF

     IF ( PRESENT( Echo ) ) THEN
        IF ( Echo ) WRITE( 6, '(a)' ) readstr(1:inlen)
     END IF

     i = 1
     DO WHILE(i <= inlen )
        IF (readstr(i:i) /= ' ' .AND. ICHAR(readstr(i:i))/=Tab ) EXIT
        i = i + 1
     END DO

     InsideQuotes = .FALSE.

     IF ( PRESENT(literal) ) THEN
       IF ( readstr(i:i) == '"' ) literal=.TRUE.
     END IF

     k = 1
     DO WHILE( i<=inlen )
        IF ( readstr(i:i) == '"' ) THEN
          InsideQuotes = .NOT.InsideQuotes
          i=i+1
          IF ( i>inlen ) EXIT
        END IF

        IF ( .NOT.InsideQuotes ) THEN
           IF ( readstr(i:i) == '!' .OR. readstr(i:i) == '#' .OR. &
                readstr(i:i) == '=' .OR. readstr(i:i) == '(' .OR. &
                readstr(i:i) == ';' .OR. readstr(i:i+1) == '::' ) EXIT 
           IF (ICHAR( readstr(i:i))<32.AND.ICHAR(readstr(i:i))/=Tab) EXIT
        END IF

        DO WHILE( i <= inlen )
          IF ( readstr(i:i) == '"'  ) THEN
            InsideQuotes = .NOT.InsideQuotes
            i=i+1
            IF ( i>inlen ) EXIT
          END IF

          IF ( .NOT.InsideQuotes ) THEN
             IF ( readstr(i:i) == ' ' .OR. readstr(i:i) == '=' .OR. &
                  readstr(i:i) == ';' .OR. readstr(i:i) == '(' .OR. &
                  readstr(i:i+1) == '::' ) EXIT 
             IF ( ICHAR( readstr(i:i))<32 ) EXIT
          END IF

          IF ( k>outlen ) THEN
             temp = str
             DEALLOCATE(str)
             outlen=LEN(temp)+512
             ALLOCATE(CHARACTER(outlen)::str)
             str(1:LEN(temp))=temp; str(LEN(temp)+1:)=''
             DEALLOCATE(temp)
          END IF

          j = ICHAR( readstr(i:i) )
          IF ( .NOT.InsideQuotes .AND. j>=A .AND. j<=Z ) THEN
            str(k:k) = CHAR(j+U2L)
          ELSE IF ( .NOT.InsideQuotes .AND. j==Tab ) THEN
            str(k:k) = ' '
          ELSE
            str(k:k) = readstr(i:i)
          ENDIF

          i = i + 1
          k = k + 1
        END DO

        IF ( k <= outlen ) str(k:k) = ' '
        k = k + 1

        DO WHILE( i<=inlen )
          IF ( readstr(i:i) /= ' ' .AND. ICHAR(readstr(i:i))/=Tab ) EXIT
          i = i + 1
        END DO
     END DO
     str(k:)=' '

     IF ( i <= inlen ) THEN
       Prefix = ' '
       IF ( ReadStr(i:i) == '=' ) THEN
         ValueStarts = i + 1
       ELSE IF ( ReadStr(i:i) == ';' ) THEN
         ValueStarts = i + 1
       ELSE IF ( ReadStr(i:i) == '(' ) THEN
         ValueStarts = i + 1
         Prefix = 'Size'
       ELSE IF ( ReadStr(i:i+1) == '::' ) THEN
         ValueStarts = i + 2
         Prefix = '::'
       ELSE IF ( ICHAR(readstr(i:i)) < 32 ) THEN
         DO WHILE( i <= inlen )
           IF ( ICHAR(readstr(i:i)) >= 32 ) EXIT
           i = i + 1
         END DO
         IF ( i <= inlen ) ValueStarts = i
       END IF
     END IF
     RETURN
10   CONTINUE
     l = .FALSE.
!------------------------------------------------------------------------------
   END FUNCTION ReadAndTrim
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  PURE FUNCTION GetVarName(Var) RESULT(str)
!------------------------------------------------------------------------------
    TYPE(Variable_t), INTENT(in) :: Var
    CHARACTER(LEN=Var % NameLen) :: str
!------------------------------------------------------------------------------
    str = Var % Name(1:Var % NameLen)
!------------------------------------------------------------------------------
  END FUNCTION GetVarName
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION ComponentNameVar( Var, Component ) RESULT(str)
!------------------------------------------------------------------------------
    TYPE(Variable_t),INTENT(in) :: Var
    INTEGER, OPTIONAL,INTENT(in) :: Component
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    IF ( Var % Name(1:Var % NameLen) == 'flow solution' ) THEN
      str='flow solution'
      IF ( .NOT. PRESENT(Component) ) RETURN
      IF ( Component == Var % DOFs ) THEN
        str = 'pressure'
        RETURN
      ELSE
        str = 'velocity ' // TRIM(i2s(Component))
      END IF
    ELSE
      str = ComponentName(Var % Name, Component)
    END IF
!------------------------------------------------------------------------------
END FUNCTION ComponentNameVar
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION ComponentNameStr( BaseName, Component_arg ) RESULT(str)
!------------------------------------------------------------------------------
    INTEGER, OPTIONAL, INTENT(in) :: Component_arg
    CHARACTER(LEN=*), INTENT(in) :: BaseName
!------------------------------------------------------------------------------
    INTEGER :: ind, ind1, DOFsTot, DOFs, Component
    CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
    ind = INDEX( BaseName,'[' )

    Component = 0
    IF ( PRESENT(Component_arg) ) Component=Component_arg

    IF ( ind<=0 ) THEN
      str = BaseName
      IF ( Component > 0 ) THEN
        str = TRIM(str) // ' ' // TRIM(i2s(Component) )
      END IF
    ELSE IF( Component == 0 ) THEN
      str = BaseName(1:ind-1)
    ELSE
      DOFsTot = 0
      DO WHILE( .TRUE. )
        ind1 = INDEX( BaseName(ind+1:),':' )+ind
        IF ( ind1 <= ind ) THEN
           CALL Fatal( 'ComponentName', 'Syntax error in variable definition.' )
        END IF
        READ(BaseName(ind1+1:),'(i1)') DOFs
        DOFsTot = DOFsTot+DOFs
        IF ( DOFsTot>=Component ) EXIT
        ind = ind1+2
      END DO
      str = BaseName(ind+1:ind1-1)
      IF ( DOFs>1 ) THEN
        DOFs = Component - DOFsTot + DOFs
        str = TRIM(str) // ' ' // TRIM(i2s(DOFs) )
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ComponentNameStr
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a tridiagonal linear system. 
!------------------------------------------------------------------------------
    PURE SUBROUTINE SolveTriDiag( n, y, h, r )
!------------------------------------------------------------------------------
       INTEGER, INTENT(in) :: n
       REAL(KIND=dp), INTENT(out) :: r(:)
       REAL(KIND=dp), INTENT(in)  :: y(:), h(:)

       REAL(KIND=dp) :: s,b(n)
       INTEGER :: i

       DO i=2,n-1
         b(i) = 2 * ( h(i-1) + h(i) )
         r(i) = 3 * ( h(i)   * ( y(i)-y(i-1) ) / h(i-1) + &
                      h(i-1) * ( y(i+1)-y(i) ) / h(i) )
       END DO

       r(2) = r(2) - h(2) * r(1)
       DO i=2,n-2
         s = -h(i+1) / b(i)
         r(i+1) = r(i+1) + s * r(i)
         b(i+1) = b(i+1) + s * h(i-1)
       END DO

       DO i=n-1,2,-1
          r(i) = (r(i) - h(i-1) * r(i+1)) / b(i)
       END DO
!------------------------------------------------------------------------------
    END SUBROUTINE SolveTriDiag
!------------------------------------------------------------------------------


    FUNCTION CheckMonotone(n,x) RESULT ( Monotone )
      REAL(KIND=dp), INTENT(in) :: x(:)
      INTEGER, INTENT(in) :: n
      LOGICAL :: Monotone
      
      INTEGER :: i
      
      Monotone = .TRUE.
      DO i=1,n-1
        IF( x(i+1) <= x(i) ) THEN
          Monotone = .FALSE.
          WRITE (Message,'(E14.7,A,E14.7)')  x(i),'>=',x(i+1)
          CALL WARN('CheckMonotone', Message)
          EXIT
        END IF
      END DO           
      
    END FUNCTION CheckMonotone

!------------------------------------------------------------------------------
!> Solver for the coefficients of a cubic spline.
!------------------------------------------------------------------------------
    PURE SUBROUTINE CubicSpline( n,x,y,r, monotone )
!------------------------------------------------------------------------------
      REAL(KIND=dp), INTENT(in)  :: x(:),y(:)
      REAL(KIND=dp), INTENT(out) :: r(:)
      INTEGER, INTENT(in) :: n
      LOGICAL, OPTIONAL, INTENT(in) :: monotone

      REAL(KIND=dp) ::  t,h(n),tau, alpha, beta
      INTEGER :: i
      LOGICAL :: mono

      DO i=1,n-1
        h(i) = x(i+1) - x(i)
      END DO

      r(1) = (y(2) - y(1) )  / h(1)
      r(n) = (y(n) - y(n-1) ) / h(n-1)

      mono = .FALSE.
      IF(PRESENT(monotone)) mono = Monotone

      IF (mono) THEN
        DO i=1,n-1
          h(i) = (y(i+1) - y(i) ) / h(i)
        END DO

        DO i=2,n-1
          r(i) = (h(i-1) + h(i))/2
        END DO

        DO i=1,n-1
          IF(ABS(h(i))<10*AEPS) THEN
            r(i) = 0._dp; r(i+1) = 0._dp
            CYCLE
          END IF

          alpha = r(i) / h(i); beta = r(i+1) / h(i)
          IF ( alpha < 0._dp .OR. beta < 0._dp ) THEN
            r(i) = 0._dp;
            CYCLE
          END IF
 
          tau = SQRT(alpha**2 + beta**2)
          IF(tau > 3) THEN
            tau = 3._dp / tau
            r(i) = alpha*tau*h(i); r(i+1)=beta*tau*h(i)
          END IF
        END DO
      ELSE
        CALL SolveTriDiag( n,y,h,r )
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE CubicSpline
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Evalulate a cubic spline.
!------------------------------------------------------------------------------
   PURE FUNCTION CubicSplineVal(x,y,r,t) RESULT(s)
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s
      REAL(KIND=dp), INTENT(in) :: x(:),y(:),r(:),t
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: a,b,c,d,h,lt

      h = x(2)-x(1)
      a = -2 * ( y(2) - y(1) ) + (   r(1) + r(2) ) * h
      b =  3 * ( y(2) - y(1) ) - ( 2*r(1) + r(2) ) * h
      c = r(1) * h
      d = y(1)

      lt = (t - x(1)) / h
      s = ((a*lt + b) * lt + c) * lt + d
!------------------------------------------------------------------------------
   END FUNCTION CubicSplineVal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Evalulate derivative of cubic spline.
!------------------------------------------------------------------------------
   PURE FUNCTION CubicSplinedVal(x,y,r,t) RESULT(s)
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: s
      REAL(KIND=dp), INTENT(in) :: x(:),y(:),r(:),t
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: a,b,c,h,lt

      h = x(2)-x(1)
      a = -2 * ( y(2) - y(1) ) + (   r(1) + r(2) ) * h
      b =  3 * ( y(2) - y(1) ) - ( 2*r(1) + r(2) ) * h
      c = r(1) * h

      lt = (t - x(1)) / h
      s = ((3*a*lt + 2*b) * lt + c)/h
!------------------------------------------------------------------------------
   END FUNCTION CubicSplinedVal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Search array index such that tval(i) <= t < tval(i+1)
!------------------------------------------------------------------------------
   PURE FUNCTION SearchInterval( tval, t ) RESULT(i)
!------------------------------------------------------------------------------
      INTEGER :: i
      REAL(KIND=dp), INTENT(in) :: tval(:), t
!------------------------------------------------------------------------------
      INTEGER :: n,n0,n1
!------------------------------------------------------------------------------

      n = SIZE(tval)

      IF (t < tval(2)) THEN
        i = 1
      ELSE IF (t>=tval(n-1)) THEN
        i = n-1
      ELSE
        n0 = 1
        n1 = n
        i = (n0+n1)/2
        DO WHILE(.TRUE.)
          IF  ( tval(i) <= t .AND. tval(i+1)>t ) EXIT

          IF ( tval(i) >  t ) THEN
            n1 = i-1 
          ELSE
            n0 = i+1
          END IF
          i = (n0+n1)/2
        END DO
      END IF
      IF(i>n-1) i=n-1
      
!------------------------------------------------------------------------------
   END FUNCTION SearchInterval
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> As SearchInterval, but doesn't assume we'll find the value in the interval
   !------------------------------------------------------------------------------
   PURE FUNCTION SearchIntPosition( tval, t ) RESULT(i)
     !------------------------------------------------------------------------------
     INTEGER :: i
     INTEGER, INTENT(in) :: tval(:), t
     !------------------------------------------------------------------------------
     INTEGER :: n,n0,n1
     !------------------------------------------------------------------------------

     n = SIZE(tval)

     IF (t < tval(1)) THEN
       i = 0
     ELSE IF (t>=tval(n)) THEN
       i = n
     ELSE
       n0 = 1
       n1 = n
       i = (n0+n1)/2
       DO WHILE(.TRUE.)
         IF  ( tval(i) <= t .AND. tval(i+1)>t ) EXIT

         IF ( tval(i) >  t ) THEN
           n1 = i-1 
         ELSE
           n0 = i+1
         END IF
         i = (n0+n1)/2
       END DO
     END IF
     IF(i>n) i=n

   !------------------------------------------------------------------------------
   END FUNCTION SearchIntPosition
   !------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Interpolate values in a curve given by linear table or splines.
!------------------------------------------------------------------------------
   PURE FUNCTION InterpolateCurve( TValues,FValues,T, CubicCoeff) RESULT( F )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: F
     REAL(KIND=dp), INTENT(iN) :: TValues(:),FValues(:),T
     REAL(KIND=dp), OPTIONAL, POINTER, INTENT(in) :: CubicCoeff(:)
!------------------------------------------------------------------------------
     INTEGER :: i,n
     LOGICAL :: Cubic
!------------------------------------------------------------------------------

     n = SIZE(TValues)

     ! This is a misuse of the interpolation in case of standard dependency
     ! of type y=a*x.  
     IF( n == 1 ) THEN
       F = FValues(1) * T
       RETURN
     END IF

     i = SearchInterval( Tvalues, t )

     Cubic = PRESENT(CubicCoeff)
     Cubic = Cubic .AND. T>=Tvalues(1) .AND. T<=Tvalues(n)
     IF ( Cubic ) Cubic = Cubic.AND.ASSOCIATED(CubicCoeff)

     IF ( Cubic ) THEN
       F = CubicSplineVal(Tvalues(i:i+1),FValues(i:i+1),CubicCoeff(i:i+1),T)
     ELSE
       F = (T-TValues(i)) / (TValues(i+1)-TValues(i))
       F = (1-F)*FValues(i) + F*FValues(i+1)
     END IF
   END FUNCTION InterpolateCurve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Derivate a curve given by linear table or splines.
!------------------------------------------------------------------------------
   PURE FUNCTION DerivateCurve( TValues,FValues,T,CubicCoeff ) RESULT( F )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: F
     REAL(KIND=dp), INTENT(in) :: TValues(:),FValues(:),T
     REAL(KIND=dp), OPTIONAL, POINTER, INTENT(in) :: CubicCoeff(:)
!------------------------------------------------------------------------------
     INTEGER :: i,n
     LOGICAL :: Cubic
!------------------------------------------------------------------------------
     n = SIZE(TValues)

     i = SearchInterval( Tvalues, t )

     Cubic = PRESENT(CubicCoeff)
     Cubic = Cubic .AND. T>=Tvalues(1) .AND. T<=Tvalues(n)
     IF ( Cubic ) Cubic = Cubic.AND.ASSOCIATED(CubicCoeff)

     IF (Cubic) THEN
       F = CubicSplinedVal(Tvalues(i:i+1),FValues(i:i+1),CubicCoeff(i:i+1),T)
     ELSE
       F = (FValues(i+1)-FValues(i)) / (TValues(i+1)-TValues(i))
     END IF
!------------------------------------------------------------------------------
   END FUNCTION DerivateCurve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Integrate a curve given by linear table or splines.
!------------------------------------------------------------------------------
   PURE SUBROUTINE CumulativeIntegral(TValues,FValues,CubicCoeff,Cumulative)
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(in)  :: TValues(:),FValues(:)
     REAL(KIND=dp), INTENT(out) :: Cumulative(:)
     REAL(KIND=dp), OPTIONAL, POINTER, INTENT(in) :: CubicCoeff(:)
!------------------------------------------------------------------------------
     INTEGER :: i,n
     LOGICAL :: Cubic
     REAL(KIND=dp) :: t(2), y(2), r(2), h, a, b, c, d
!------------------------------------------------------------------------------
     n = SIZE(TValues)

     Cubic = PRESENT(CubicCoeff)
     IF ( Cubic ) Cubic = Cubic.AND.ASSOCIATED(CubicCoeff)

     ! here only complete intervals:
     ! -----------------------------
     Cumulative(1) = 0._dp
     IF ( Cubic ) THEN
       DO i=1,n-1
         t(1) = Tvalues(i)
         t(2) = Tvalues(i+1)

         y(1) = FValues(i)
         y(2) = FValues(i+1)

         r(1) = CubicCoeff(i)
         r(2) = CubicCoeff(i+1)

         h  = t(2) - t(1)

         a = (-2 * (y(2) - y(1)) + (  r(1) + r(2)) * h)/4
         b = ( 3 * (y(2) - y(1)) - (2*r(1) + r(2)) * h)/3
         c = (r(1) * h)/2
         d = y(1)
         Cumulative(i+1) = Cumulative(i) + h*(a+b+c+d)
       END DO
     ELSE
       DO i=1,n-1
         t(1) = Tvalues(i)
         t(2) = Tvalues(i+1)

         y(1) = FValues(i)
         y(2) = FValues(i+1)

         h  = t(2) - t(1)
         c = (y(2)-y(1))/2
         d = y(1)
         Cumulative(i+1) = Cumulative(i) + h*(c+d)
       END DO
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CumulativeIntegral
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Integrate a curve given by linear table or splines.
!------------------------------------------------------------------------------
   PURE FUNCTION IntegrateCurve(TValues,FValues,CubicCoeff,T0,T1,Cumulative) RESULT(sumf)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: sumf

     REAL(KIND=dp), INTENT(in) :: TValues(:),FValues(:)
     REAL(KIND=dp), OPTIONAL, INTENT(in) :: T0, T1
     REAL(KIND=dp), OPTIONAL, INTENT(in) :: Cumulative(:)
     REAL(KIND=dp), OPTIONAL, POINTER, INTENT(in) :: CubicCoeff(:)
!------------------------------------------------------------------------------
     INTEGER :: i,n,i0,i1
     LOGICAL :: Cubic
     REAL(KIND=dp) :: t(2), y(2), r(2), h, a, b, c, d, s0, s1, tt0, tt1
!------------------------------------------------------------------------------
     n = SIZE(TValues)

     tt0 = TValues(1)
     IF(PRESENT(t0)) tt0=t0

     tt1 = TValues(n)
     IF(PRESENT(t1)) tt1=t1

     sumf = 0._dp
     IF(tt0>=tt1) RETURN

     ! t0 < first, t1 <= first
     IF(tt1<=Tvalues(1)) THEN
       t(1) = Tvalues(1)
       t(2) = Tvalues(2)

       y(1) = FValues(1)
       y(2) = FValues(2)

       h  = t(2) - t(1)
       s0 = (tt0 - t(1)) / h
       s1 = (tt1 - t(1)) / h
       c = (y(2) - y(1)) / 2
       d = y(1)
       sumf = sumf + h * ((c*s1 + d)*s1 - (c*s0 + d)*s0)
       RETURN
     END IF

     ! t0 >= last, t1 > last
     IF(tt0>=Tvalues(n)) THEN
       t(1) = Tvalues(n-1)
       t(2) = Tvalues(n)

       y(1) = FValues(n-1)
       y(2) = FValues(n)

       h  = t(2) - t(1)
       s0 = (tt0 - t(1)) / h
       s1 = (tt1 - t(1)) / h
       c = (y(2) - y(1)) / 2
       d = y(1)
       sumf = sumf + h * ((c*s1 + d)*s1 - (c*s0 + d)*s0)
       RETURN
     END IF

     ! first interval outside 
     IF(tt0<Tvalues(1)) THEN
       t(1) = Tvalues(1)
       t(2) = Tvalues(2)

       y(1) = FValues(1)
       y(2) = FValues(2)

       h  = t(2) - t(1)
       s0 = (tt0 - t(1)) / h
       c = (y(2) - y(1)) / 2
       d = y(1)
       sumf = sumf - h * (c*s0 + d)*s0
       tt0 = Tvalues(1)
     END IF

     ! last interval outside 
     IF(tt1>Tvalues(n)) THEN
       t(1) = Tvalues(n-1)
       t(2) = Tvalues(n)

       y(1) = FValues(n-1)
       y(2) = FValues(n)

       h  = t(2) - t(1)
       s1 = (tt1 - t(1)) / h
       c = (y(2) - y(1)) / 2
       d = y(1)
       sumf = sumf + h * ( (c*s1 + d)*s1 - (c+d) )
       tt1 = Tvalues(n)
     END IF

     IF(tt0 >= tt1) RETURN

     Cubic = PRESENT(CubicCoeff)
     IF ( Cubic ) Cubic = Cubic.AND.ASSOCIATED(CubicCoeff)

     i0 = SearchInterval( Tvalues, tt0 )

     ! first (possibly partial) interval:
     ! -------------------------------------
     t(1) = Tvalues(i0)
     t(2) = Tvalues(i0+1)

     h  = t(2) - t(1)
     s0 = (tt0-t(1))/h
     s1 = MIN((tt1-t(1))/h,1._dp)

     IF(s0>0 .OR. s1<1) THEN
       y(1) = FValues(i0)
       y(2) = FValues(i0+1)

       IF(Cubic) THEN
         r(1) = CubicCoeff(i0)
         r(2) = CubicCoeff(i0+1)

         a = (-2 * (y(2) - y(1)) + (  r(1) + r(2)) * h)/4
         b = ( 3 * (y(2) - y(1)) - (2*r(1) + r(2)) * h)/3
         c = (r(1) * h)/2
         d = y(1)
         sumf = sumf + h * ( (((a*s1 + b)*s1 + c)*s1 + d)*s1 - &
                   (((a*s0 + b)*s0 + c)*s0 + d)*s0 )

       ELSE
         c = (y(2)-y(1))/2
         d = y(1)
         sumf = sumf + h * ( (c*s1 + d)*s1 - (c*s0 + d)*s0 )
       END IF
       i0 = i0 + 1 
       tt0 = Tvalues(i0)
       IF(tt0 >= tt1) RETURN
     END IF

     i1 = SearchInterval( Tvalues, tt1 )

     ! last (possibly partial) interval:
     ! ------------------------------------
     t(1) = Tvalues(i1)
     t(2) = Tvalues(i1+1)

     h  = t(2) - t(1)

     s0 = MAX((tt0-t(1))/h, 0.0_dp)
     s1 = (tt1-t(1))/h

     IF(s0>0 .OR. s1<1) THEN
       y(1) = FValues(i1)
       y(2) = FValues(i1+1)

       IF(Cubic) THEN
         r(1) = CubicCoeff(i1)
         r(2) = CubicCoeff(i1+1)

         a = (-2 * (y(2) - y(1)) + (  r(1) + r(2)) * h)/4
         b = ( 3 * (y(2) - y(1)) - (2*r(1) + r(2)) * h)/3
         c = (r(1) * h)/2
         d = y(1)
         sumf = sumf + h * ( (((a*s1 + b)*s1 + c)*s1 + d)*s1 - &
                 (((a*s0 + b)*s0 + c)*s0 + d)*s0 )
       ELSE
         c = (y(2)-y(1))/2
         d = y(1)
         sumf = sumf + h * ( (c*s1 + d)*s1 - (c*s0 + d)*s0 )
       END IF
       i1 = i1 - 1 
       tt1 = Tvalues(i1+1)
       IF(tt0 >= tt1) RETURN
     END IF

     ! here only complete intervals:
     ! -----------------------------

     IF(PRESENT(Cumulative)) THEN
       sumf = sumf + Cumulative(i1+1) - Cumulative(i0)
       RETURN
     END IF

     IF ( Cubic ) THEN
       DO i=i0,i1
         t(1) = Tvalues(i)
         t(2) = Tvalues(i+1)

         y(1) = FValues(i)
         y(2) = FValues(i+1)

         r(1) = CubicCoeff(i)
         r(2) = CubicCoeff(i+1)

         h  = t(2) - t(1)

         a = (-2 * (y(2) - y(1)) + (  r(1) + r(2)) * h)/4
         b = ( 3 * (y(2) - y(1)) - (2*r(1) + r(2)) * h)/3
         c = (r(1) * h)/2
         d = y(1)
         sumf = sumf + h * (a+b+c+d)
       END DO
     ELSE
       DO i=i0,i1
         t(1) = Tvalues(i)
         t(2) = Tvalues(i+1)

         y(1) = FValues(i)
         y(2) = FValues(i+1)

         h  = t(2) - t(1)
         c = (y(2)-y(1))/2
         d = y(1)
         sumf = sumf + h * (c+d)
       END DO
     END IF
!------------------------------------------------------------------------------
   END FUNCTION IntegrateCurve
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a 2 x 2 linear system.
!------------------------------------------------------------------------------
   SUBROUTINE SolveLinSys2x2( A, x, b )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(out) :: x(:)
     REAL(KIND=dp), INTENT(in)  :: A(:,:),b(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: detA
!------------------------------------------------------------------------------
     detA = A(1,1) * A(2,2) - A(1,2) * A(2,1)

     IF ( detA == 0.0d0 ) THEN
       WRITE( Message, * ) 'Singular matrix, sorry!'
       CALL Error( 'SolveLinSys2x2', Message )
       RETURN
     END IF

     detA = 1.0d0 / detA
     x(1) = detA * (A(2,2) * b(1) - A(1,2) * b(2))
     x(2) = detA * (A(1,1) * b(2) - A(2,1) * b(1))
!------------------------------------------------------------------------------
   END SUBROUTINE SolveLinSys2x2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a 3 x 3 linear system.
!------------------------------------------------------------------------------
   SUBROUTINE SolveLinSys3x3( A, x, b )
!------------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(out) :: x(:)
     REAL(KIND=dp), INTENT(in)  :: A(:,:),b(:)
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: C(2,2),y(2),g(2),s,t,q
!------------------------------------------------------------------------------

     IF ( ABS(A(1,1))>ABS(A(1,2)) .AND. ABS(A(1,1))>ABS(A(1,3)) ) THEN
       q = 1.0d0 / A(1,1)
       s = q * A(2,1)
       t = q * A(3,1)
       C(1,1) = A(2,2) - s * A(1,2)
       C(1,2) = A(2,3) - s * A(1,3)
       C(2,1) = A(3,2) - t * A(1,2)
       C(2,2) = A(3,3) - t * A(1,3)

       g(1) = b(2) - s * b(1)
       g(2) = b(3) - t * b(1)
       CALL SolveLinSys2x2( C,y,g )
       
       x(2) = y(1)
       x(3) = y(2)
       x(1) = q * ( b(1) - A(1,2) * x(2) - A(1,3) * x(3) )
     ELSE IF ( ABS(A(1,2)) > ABS(A(1,3)) ) THEN
       q = 1.0d0 / A(1,2)
       s = q * A(2,2)
       t = q * A(3,2)
       C(1,1) = A(2,1) - s * A(1,1)
       C(1,2) = A(2,3) - s * A(1,3)
       C(2,1) = A(3,1) - t * A(1,1)
       C(2,2) = A(3,3) - t * A(1,3)
       
       g(1) = b(2) - s * b(1)
       g(2) = b(3) - t * b(1)
       CALL SolveLinSys2x2( C,y,g )

       x(1) = y(1)
       x(3) = y(2)
       x(2) = q * ( b(1) - A(1,1) * x(1) - A(1,3) * x(3) )
     ELSE
       q = 1.0d0 / A(1,3)
       s = q * A(2,3)
       t = q * A(3,3)
       C(1,1) = A(2,1) - s * A(1,1)
       C(1,2) = A(2,2) - s * A(1,2)
       C(2,1) = A(3,1) - t * A(1,1)
       C(2,2) = A(3,2) - t * A(1,2)

       g(1) = b(2) - s * b(1)
       g(2) = b(3) - t * b(1)
       CALL SolveLinSys2x2( C,y,g )

       x(1) = y(1)
       x(2) = y(2)
       x(3) = q * ( b(1) - A(1,1) * x(1) - A(1,2) * x(2) )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE SolveLinSys3x3
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE ClearMatrix( Matrix ) 
     TYPE(Matrix_t), POINTER, INTENT(in) :: Matrix
INCLUDE "mpif.h"
  
     Matrix % FORMAT = MATRIX_CRS

      NULLIFY( Matrix % Child )
      NULLIFY( Matrix % Parent )
      NULLIFY( Matrix % EMatrix )
      NULLIFY( Matrix % ConstraintMatrix )

      NULLIFY( Matrix % Perm )
      NULLIFY( Matrix % InvPerm )

      NULLIFY( Matrix % Cols )
      NULLIFY( Matrix % Rows )
      NULLIFY( Matrix % Diag )
 
      NULLIFY( Matrix % RHS )
      NULLIFY( Matrix % Force )
      NULLIFY( Matrix % RHS_im )

      NULLIFY( Matrix % Values )
      NULLIFY( Matrix % ILUValues )
      NULLIFY( Matrix % MassValues )
      NULLIFY( Matrix % DampValues )

      NULLIFY( Matrix % BulkRHS )
      NULLIFY( Matrix % BulkValues )

      NULLIFY( Matrix % BulkResidual )

      NULLIFY( Matrix % ILUCols )
      NULLIFY( Matrix % ILURows )
      NULLIFY( Matrix % ILUDiag )

      NULLIFY( Matrix % CRHS )
      NULLIFY( Matrix % CForce )

      NULLIFY( Matrix % ParMatrix )

      NULLIFY( Matrix % CValues )
      NULLIFY( Matrix % CILUValues )
      NULLIFY( Matrix % CMassValues )
      NULLIFY( Matrix % CDampValues )

!     NULLIFY( Matrix % GRows )
!     NULLIFY( Matrix % RowOwner )
      NULLIFY( Matrix % GOrder )
      NULLIFY( Matrix % EPerm )

      NULLIFY( Matrix % ParMatrix )

      NULLIFY( Matrix % ParallelInfo )
#ifdef HAVE_UMFPACK
      Matrix % UMFPack_Numeric = 0
#endif

      Matrix % Cholesky  = .FALSE.
      Matrix % Lumped    = .FALSE.
      Matrix % Ordered   = .FALSE. 
      Matrix % COMPLEX   = .FALSE.
      Matrix % Symmetric = .FALSE.
      Matrix % SolveCount   = 0
      Matrix % NumberOfRows = 0

      Matrix % ProjectorBC = 0
      Matrix % ProjectorType = PROJECTOR_TYPE_DEFAULT
      
      Matrix % Solver => NULL()

      Matrix % DGMatrix = .FALSE.
      Matrix % Comm = ELMER_COMM_WORLD

   END SUBROUTINE ClearMatrix


!------------------------------------------------------------------------------
   FUNCTION AllocateMatrix() RESULT(Matrix)
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: Matrix
!------------------------------------------------------------------------------
      ALLOCATE( Matrix )
      CALL ClearMatrix( Matrix )
!------------------------------------------------------------------------------
   END FUNCTION AllocateMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE FreeQuadrantTree( Root )
!------------------------------------------------------------------------------
    TYPE(Quadrant_t), POINTER :: Root

    INTEGER :: i

    IF ( .NOT. ASSOCIATED( Root ) ) RETURN

    IF ( ASSOCIATED(Root % Elements) ) DEALLOCATE( Root % Elements )

    IF ( ASSOCIATED( Root % ChildQuadrants ) ) THEN
       DO i=1,SIZE(Root % ChildQuadrants)
          CALL FreeQuadrantTree( Root % ChildQuadrants(i) % Quadrant )
       END DO
       DEALLOCATE( Root % ChildQuadrants )
       NULLIFY( Root % ChildQuadrants )
    END IF

    DEALLOCATE( Root )
!------------------------------------------------------------------------------
  END SUBROUTINE FreeQuadrantTree
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateRealVector( F, n, From, FailureMessage )
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: F(:)
    INTEGER :: n
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n > 0 ) THEN
       ALLOCATE( F(n), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
       IF ( PRESENT( FailureMessage  ) ) THEN
          WRITE( Message, * )'Unable to allocate ', n, ' element real array.'
          CALL Error( 'AllocateRealVector', Message )
          IF ( PRESENT( From ) ) THEN
             WRITE( Message, * )'Requested From: ', TRIM(From)
             CALL Error( 'AllocateRealVector', Message )
          END IF
          IF ( PRESENT( FailureMessage ) ) THEN
             CALL Fatal( 'AllocateRealVector', FailureMessage )
          END IF
       END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateRealVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateComplexVector( f, n, From, FailureMessage )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), POINTER :: f(:)
    INTEGER :: n
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n > 0 ) THEN
       ALLOCATE( f(n), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
       IF ( PRESENT( FailureMessage  ) ) THEN
          WRITE( Message, * )'Unable to allocate ', n, ' element real array.'
          CALL Error( 'AllocateComplexVector', Message )
          IF ( PRESENT( From ) ) THEN
             WRITE( Message, * )'Requested From: ', TRIM(From)
             CALL Error( 'AllocateComplexVector', Message )
          END IF
          IF ( PRESENT( FailureMessage ) ) THEN
             CALL Fatal( 'AllocateComplexVector', FailureMessage )
          END IF
       END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateComplexVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateIntegerVector( f, n, From, FailureMessage )
!------------------------------------------------------------------------------
    INTEGER, POINTER :: f(:)
    INTEGER :: n
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n > 0 ) THEN
       ALLOCATE( f(n), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
       IF ( PRESENT( FailureMessage  ) ) THEN
          WRITE( Message, * )'Unable to allocate ', n, ' element integer array.'
          CALL Error( 'AllocateIntegerVector', Message )
          IF ( PRESENT( From ) ) THEN
             WRITE( Message, * )'Requested From: ', TRIM(From)
             CALL Error( 'AllocateIntegerVector', Message )
          END IF
          IF ( PRESENT( FailureMessage ) ) THEN
             CALL Fatal( 'AllocateIntegerVector', FailureMessage )
          END IF
       END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateIntegerVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateLogicalVector( f, n, From, FailureMessage )
!------------------------------------------------------------------------------
    LOGICAL, POINTER :: f(:)
    INTEGER :: n
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n > 0 ) THEN
       ALLOCATE( f(n), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
       IF ( PRESENT( FailureMessage  ) ) THEN
          WRITE( Message, * )'Unable to allocate ', n, ' element integer array.'
          CALL Error( 'AllocateLogicalVector', Message )
          IF ( PRESENT( From ) ) THEN
             WRITE( Message, * )'Requested From: ', TRIM(From)
             CALL Error( 'AllocateLogicalVector', Message )
          END IF
          IF ( PRESENT( FailureMessage ) ) THEN
             CALL Fatal( 'AllocateLogicalVector', FailureMessage )
          END IF
       END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateLogicalVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateElementVector( f, n, From, FailureMessage )
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: f(:)
    INTEGER :: n
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n > 0 ) THEN
       ALLOCATE( f(n), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
      WRITE( Message, * )'Unable to allocate ', n, ' element integer array.'
      CALL Error( 'AllocateElementVector', Message )
      IF ( PRESENT( From ) ) THEN
        WRITE( Message, * )'Requested From: ', TRIM(From)
        CALL Error( 'AllocateElementVector', Message )
      END IF
      IF ( PRESENT( FailureMessage ) ) THEN
        CALL Fatal( 'AllocateElementVector', FailureMessage )
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateElementVector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateRealArray( f, n1, n2, From, FailureMessage )
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: f(:,:)
    INTEGER :: n1,n2
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n1 > 0 .AND. n2 > 0 ) THEN
       ALLOCATE( f(n1,n2), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
      WRITE( Message, * )'Unable to allocate ', n1, ' by ', n2, ' element real matrix.'
      CALL Error( 'AllocateRealArray', Message )
      IF ( PRESENT( From ) ) THEN
        WRITE( Message, * )'Requested From: ', TRIM(From)
        CALL Error( 'AllocateRealArray', Message )
      END IF
      IF ( PRESENT( FailureMessage ) ) THEN
        CALL Fatal( 'AllocateRealArray', FailureMessage )
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE  AllocateRealArray
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AllocateComplexArray( f, n1, n2, From, FailureMessage )
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp), POINTER :: f(:,:)
    INTEGER :: n1,n2
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n1 > 0 .AND. n2 > 0 ) THEN
       ALLOCATE( f(n1,n2), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
      WRITE( Message, * )'Unable to allocate ', n1, ' by ', n2, ' element real matrix.'
      CALL Error( 'AllocateComplexArray', Message )
      IF ( PRESENT( From ) ) THEN
        WRITE( Message, * )'Requested From: ', TRIM(From)
        CALL Error( 'AllocateComplexArray', Message )
      END IF
      IF ( PRESENT( FailureMessage ) ) THEN
        CALL Fatal( 'AllocateComplexArray', FailureMessage )
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE  AllocateComplexArray
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateIntegerArray( f, n1, n2, From, FailureMessage )
!------------------------------------------------------------------------------
    INTEGER, POINTER :: f(:,:)
    INTEGER :: n1,n2
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n1 > 0 .AND. n2 > 0 ) THEN
       ALLOCATE( f(n1,n2), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
      WRITE( Message, * )'Unable to allocate ', n1, ' by ', n2, ' element integer matrix.'
      CALL Error( 'AllocateIntegerArray', Message )
      IF ( PRESENT( From ) ) THEN
        WRITE( Message, * )'Requested From: ', TRIM(From)
        CALL Error( 'AllocateIntegerArray', Message )
      END IF
      IF ( PRESENT( FailureMessage ) ) THEN
        CALL Fatal( 'AllocateIntegerArray', FailureMessage )
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE  AllocateIntegerArray
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE AllocateLogicalArray( f, n1, n2, From, FailureMessage )
!------------------------------------------------------------------------------
    LOGICAL, POINTER :: f(:,:)
    INTEGER :: n1,n2
    CHARACTER(LEN=*), OPTIONAL :: From, FailureMessage
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

    istat = -1
    IF ( n1 > 0 .AND. n2 > 0 ) THEN
       ALLOCATE( f(n1,n2), STAT=istat )
    END IF
    IF ( istat /=  0 ) THEN
      WRITE( Message, * )'Unable to allocate ', n1, ' by ', n2, ' element integer matrix.'
      CALL Error( 'AllocateLogicalArray', Message )
      IF ( PRESENT( From ) ) THEN
        WRITE( Message, * )'Requested From: ', TRIM(From)
        CALL Error( 'AllocateLogicalArray', Message )
      END IF
      IF ( PRESENT( FailureMessage ) ) THEN
        CALL Fatal( 'AllocateLogicalArray', FailureMessage )
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE  AllocateLogicalArray
!------------------------------------------------------------------------------

  ! Pad given integer value to be the next largest multiple of nbyte
  FUNCTION IntegerNBytePad(val, nbyte) RESULT(padval)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: val, nbyte
    INTEGER :: padval
    ! Parameters and variables
    INTEGER, PARAMETER :: bytesinint = KIND(val)
    INTEGER :: nbytesinint

    ! Compute number of nbytes in int
    nbytesinint = nbyte/bytesinint
    ! Compute value padded to multiples of n-byte
    padval=((val-1)/nbytesinint)*nbytesinint+nbytesinint
  END FUNCTION IntegerNBytePad

  ! Pad given value to be the next largest multiple of nbyte
  FUNCTION NBytePad(val, bytesinelem, nbyte) RESULT(padval)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: val, bytesinelem, nbyte
    INTEGER :: padval
    ! Variables
    INTEGER :: nbytesinelem

    ! Compute number of nbytes in a single element
    nbytesinelem = nbyte/bytesinelem
    ! Compute value padded to multiples of n-byte
    padval=((val-1)/nbytesinelem)*nbytesinelem+nbytesinelem
  END FUNCTION NBytePad

!------------------------------------------------------------------------------
!> Given the filename0 (and suffix0) find the 1st free filename
!> that does not exist in the current working directory
!------------------------------------------------------------------------------

  FUNCTION NextFreeFilename(Filename0,Suffix0,LastExisting) RESULT (Filename)

    CHARACTER(LEN=MAX_NAME_LEN) :: Filename0
    CHARACTER(LEN=MAX_NAME_LEN), OPTIONAL :: Suffix0 
    LOGICAL, OPTIONAL :: LastExisting
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    CHARACTER(LEN=MAX_NAME_LEN) :: Prefix, Suffix, PrevFilename
    LOGICAL :: FileIs
    INTEGER :: No, ind, len
    
    ind = INDEX( FileName0,'.',.TRUE. )
    len = LEN_TRIM(Filename0)
    IF(ind > 0) THEN
      Prefix = Filename0(1:ind-1)
      Suffix = Filename0(ind:len)
    ELSE
      Prefix = Filename0(1:len)
      IF(PRESENT(Suffix0)) THEN
        Suffix = '.'//TRIM(Suffix0)
      ELSE
        Suffix = '.dat'
      END IF
    END IF

    DO No = 1,9999
      IF( No > 0 ) PrevFilename = Filename
      IF( No < 10) THEN
        WRITE( FileName,'(A,I1,A)') TRIM(Prefix),No,TRIM(Suffix)
      ELSE IF( No < 100) THEN
        WRITE( FileName,'(A,I2,A)') TRIM(Prefix),No,TRIM(Suffix)
      ELSE IF( No < 1000) THEN
        WRITE( FileName,'(A,I3,A)') TRIM(Prefix),No,TRIM(Suffix)
      ELSE IF( No < 10000) THEN
        WRITE( FileName,'(A,I4,A)') TRIM(Prefix),No,TRIM(Suffix)
      END IF
      INQUIRE( FILE=Filename, EXIST=FileIs )
      IF(.NOT. FileIs) EXIT
    END DO

    IF( PRESENT(LastExisting)) THEN
      IF( LastExisting ) Filename = PrevFilename
    END IF

!------------------------------------------------------------------------------
  END FUNCTION NextFreeFilename
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Given the filename0 add a string related to the partitioning.
!------------------------------------------------------------------------------

  FUNCTION AddFilenameParSuffix(Filename0,Suffix0,Parallel,MyPe) RESULT (Filename)

    CHARACTER(LEN=MAX_NAME_LEN) :: Filename0
    CHARACTER(LEN=*), OPTIONAL :: Suffix0 
    LOGICAL :: Parallel
    INTEGER :: MyPe
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    CHARACTER(LEN=MAX_NAME_LEN) :: Prefix, Suffix
    INTEGER :: No, ind, len


    ind = INDEX( FileName0,'.',.TRUE. )
    len = LEN_TRIM(Filename0)
    IF(ind > 0) THEN
      Prefix = Filename0(1:ind-1)
      Suffix = Filename0(ind:len)
    ELSE
      Prefix = Filename0(1:len)
      IF(PRESENT(Suffix0)) THEN
        Suffix = '.'//TRIM(Suffix0)
      ELSE
        Suffix = '.dat'
      END IF
    END IF

    IF( Parallel ) THEN
      No = MyPe + 1
      IF( No < 10000 ) THEN
        WRITE( FileName,'(A,I4.4,A)') TRIM(Prefix),No,TRIM(Suffix)
      ELSE
        WRITE( FileName,'(A,I0,A)') TRIM(Prefix),No,TRIM(Suffix)
      END IF
    ELSE
      FileName = TRIM(Prefix)//TRIM(Suffix)
    END IF

!------------------------------------------------------------------------------
  END FUNCTION AddFilenameParSuffix
!------------------------------------------------------------------------------


  !---------------------------------------------------------<
  !> Returns values from a normal distribution to be used in 
  !> thermal velocity distribution, for example.
  !---------------------------------------------------------
  FUNCTION NormalRandom() RESULT ( normalrand ) 
    
    REAL(KIND=dp) :: normalrand,mean
    INTEGER :: flag = 0
    REAL(KIND=dp) :: fac,gsave,rsq,r1,r2 
    
    SAVE flag,gsave 
    
    IF (flag == 0) THEN 
      rsq=2.0_dp 
      
      DO WHILE(rsq >= 1.0_dp .OR. rsq == 0.0_dp ) 
        CALL RANDOM_NUMBER(r1)
        CALL RANDOM_NUMBER(r2)
        r1 = 2.0_dp * r1 - 1.0_dp
        r2 = 2.0_dp * r2 - 1.0_dp
        rsq = r1*r1 + r2*r2 
      ENDDO
      
      fac = SQRT(-2.0_dp * LOG(rsq) / rsq) 
      gsave = r1 * fac 
      normalrand = r2 * fac 
      flag = 1 
    ELSE 
      normalrand = gsave 
      flag = 0 
    ENDIF
    
  END FUNCTION NormalRandom


  !---------------------------------------------------------
  !> Returns values from a even distribution [0,1]
  !---------------------------------------------------------
  FUNCTION EvenRandom() RESULT ( rand )     
    REAL(KIND=dp) :: rand
    CALL RANDOM_NUMBER(rand)
  END FUNCTION EvenRandom
   

  SUBROUTINE ForceLoad
    CALL MPI_SEND()
  END SUBROUTINE ForceLoad


END MODULE GeneralUtils


!---------------------------------------------------------
!> Module mainly for writing xml based vtk files. 
!> The idea is that same routines save both the ascii 
!> and binary format. 
!---------------------------------------------------------
MODULE AscBinOutputUtils
  
  
  USE Types
  IMPLICIT NONE
  
  LOGICAL, PRIVATE :: AsciiOutput, SinglePrec
  INTEGER, PRIVATE :: VtuUnit = 0, BufferSize = 0
  REAL, POINTER, PRIVATE :: FVals(:)
  REAL(KIND=dp), POINTER, PRIVATE :: DVals(:)
  INTEGER, POINTER, PRIVATE :: IVals(:)
  INTEGER, PRIVATE :: INoVals, NoVals

  SAVE :: AsciiOutput, SinglePrec, VtuUnit,  BufferSize, &
      FVals, DVals, IVals
  


CONTAINS


  ! Initialize the buffer for writing, choose mode etc.
  !-----------------------------------------------------------------
  SUBROUTINE AscBinWriteInit( IsAscii, IsSingle, UnitNo, BufSize )
    
    LOGICAL :: IsAscii, IsSingle
    INTEGER :: UnitNo, BufSize
    
    AsciiOutput =  IsAscii
    SinglePrec = IsSingle
    VtuUnit = UnitNo
    BufferSize = BufSize

    CALL Info('AscBinWriteInit','Initializing buffered ascii/binary writing',Level=8)
    IF( AsciiOutput ) THEN
      CALL Info('AscBinWriteInit','Writing in ascii',Level=10)
    ELSE
      CALL Info('AscBinWriteInit','Writing in binary',Level=10)
    END IF

    IF( SinglePrec ) THEN
      CALL Info('AscBinWriteInit','Writing in single precision',Level=10)
    ELSE
      CALL Info('AscBinWriteInit','Writing in double precision',Level=10)
    END IF

    WRITE(Message,'(A,I0)')  'Writing to unit number: ',VtuUnit
    CALL Info('AscBinWriteInit',Message,Level=10)

    IF(.NOT. AsciiOutput ) THEN
      WRITE(Message,'(A,I0)')  'Size of buffer is: ',BufferSize
      CALL Info('AscBinWriteInit',Message,Level=10)
      
      ALLOCATE( Ivals( BufferSize ) ) 
      IF( SinglePrec ) THEN
        ALLOCATE( FVals( BufferSize ) ) 
      ELSE
        ALLOCATE( Dvals( BufferSize ) ) 
      END IF
      
      INoVals = 0 
      NoVals = 0
    END IF

  END SUBROUTINE AscBinWriteInit


  ! Free the buffer, next buffer can be different in size and type
  !---------------------------------------------------------------
  SUBROUTINE AscBinWriteFree()

    CALL Info('AscBinWriteFree','Terminating buffered ascii/binary writing',Level=10)

    IF( AsciiOutput ) RETURN

    IF( SinglePrec ) THEN
      DEALLOCATE( FVals )
    ELSE
      DEALLOCATE( DVals ) 
    END IF
    DEALLOCATE( IVals ) 

    BufferSize = 0
    VtuUnit = 0
    
  END SUBROUTINE AscBinWriteFree



  ! The writing of xml strings is done here to allow easier modification
  ! of output strategies.
  !-------------------------------------------------------------------------
  SUBROUTINE AscBinStrWrite( Str )
    
    CHARACTER(LEN=1024) :: Str 
    INTEGER, PARAMETER :: VtuUnit = 58
    
    WRITE( VtuUnit ) TRIM(Str)        
    
  END SUBROUTINE AscBinStrWrite
  

  ! Write a binary value, either in single or double precision
  !------------------------------------------------------------------------
  SUBROUTINE AscBinRealWrite( val, EmptyBuffer )
    
    INTEGER, PARAMETER :: VtuUnit = 58
    REAL(KIND=dp) :: val
    LOGICAL, OPTIONAL :: EmptyBuffer
    LOGICAL :: Empty
    CHARACTER(LEN=1024) :: Str 

    IF( VtuUnit == 0 ) THEN
      CALL Fatal('AscBinRealWrite','Buffer not initialized for writing')
    END IF

    IF( PRESENT( EmptyBuffer ) ) THEN
      Empty = EmptyBuffer
    ELSE
      Empty = .FALSE.
    END IF

    ! Ascii output is not buffered
    IF( AsciiOutput ) THEN
      IF( Empty ) RETURN
      IF( ABS( val ) <= TINY ( val ) ) THEN
        WRITE(Str,'(A)') " 0.0"
      ELSE IF( SinglePrec ) THEN
        WRITE( Str,'(ES12.3E3)') val       
      ELSE
        WRITE( Str,'(ES16.7E3)') val
      END IF
      WRITE( VtuUnit ) TRIM(Str)        
      RETURN
    END IF

    ! Buffered binary output
    IF( Empty .OR. NoVals == BufferSize ) THEN
      IF( NoVals == 0 ) THEN
        RETURN
      ELSE IF( SinglePrec ) THEN
        WRITE( VtuUnit ) Fvals(1:NoVals)
      ELSE
        WRITE( VtuUnit ) DVals(1:NoVals) 
      END IF
      NoVals = 0
      IF( Empty ) RETURN 
    END IF
    
    ! Save values in the buffer (either single or double prec.)
    NoVals = NoVals + 1
    IF( SinglePrec ) THEN
      Fvals(NoVals) = val
    ELSE
      DVals(NoVals) = val
    END IF


  END SUBROUTINE AscBinRealWrite


  ! Write an integer value 
  !-------------------------------------------------
  SUBROUTINE AscBinIntegerWrite( ival, EmptyBuffer )
    
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: ival
    LOGICAL, OPTIONAL :: EmptyBuffer
    LOGICAL :: Empty
    CHARACTER(LEN=1024) :: Str 

    IF( VtuUnit == 0 ) THEN
      CALL Fatal('AscBinIntegerWrite','Buffer not initialized for writing')
    END IF

    IF( PRESENT( EmptyBuffer ) ) THEN
      Empty = EmptyBuffer
    ELSE
      Empty = .FALSE.
    END IF

    IF( AsciiOutput ) THEN
      IF( Empty ) RETURN
      WRITE( Str, '(" ",I0)') ival
      WRITE( VtuUnit ) TRIM(Str)        
      RETURN
    END IF

    IF( Empty .OR. INoVals == BufferSize ) THEN
      IF( INoVals == 0 ) THEN
        RETURN
      ELSE 
        WRITE( VtuUnit ) Ivals(1:INoVals)
      END IF
      INoVals = 0
      IF( Empty ) RETURN 
    END IF

    INoVals = INoVals + 1
    Ivals(INoVals) = ival
        
  END SUBROUTINE AscBinIntegerWrite
  
  
END MODULE AscBinOutputUtils


!> \}

