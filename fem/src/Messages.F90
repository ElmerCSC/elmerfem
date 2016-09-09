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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 2002
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{


!-------------------------------------------------------------------------------
!>  Message output routines for sending information and dealing with exceptions.
!-------------------------------------------------------------------------------
MODULE Messages

   IMPLICIT NONE
   CHARACTER(LEN=512) :: Message = ' '
   INTEGER, PRIVATE :: i
   LOGICAL :: OutputPrefix=.FALSE., OutputCaller=.TRUE.
   LOGICAL :: OutputLevelMask(0:31) = (/ (.TRUE.,i=1,32) /)
   INTEGER :: MaxOutputLevel=32, MinOutputLevel=0, OutputPE = 0
   INTEGER :: MaxOutputPE = 0, MinOutputPE = 0

CONTAINS

!-----------------------------------------------------------------------
!> Prints information on the standard output if the requested or 
!> default output level does not surpass the maximum output level.
!-----------------------------------------------------------------------
   SUBROUTINE Info( Caller, String, noAdvance, Level )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     INTEGER, OPTIONAL :: Level
     LOGICAL, OPTIONAL :: noAdvance
!-----------------------------------------------------------------------

     LOGICAL :: nadv, nadv1 = .FALSE.
     INTEGER :: n
     INTEGER, PARAMETER :: DefLevel = 4
     SAVE nadv1

!-----------------------------------------------------------------------

     IF ( OutputPE < 0 ) RETURN

     IF ( PRESENT( Level ) ) THEN
       IF ( .NOT. OutputLevelMask(Level) ) RETURN
     ELSE
       ! The default level of info
       !-------------------------------------------
       IF( .NOT. OutputLevelMask(DefLevel) ) RETURN
     END IF

     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance

     IF(.NOT. nadv1 ) THEN
       IF ( OutputPrefix ) THEN
         WRITE( *,'(A)', ADVANCE = 'NO' ) 'INFO:: '
       END IF

       IF ( OutputCaller ) THEN
         WRITE( *,'(A)', ADVANCE = 'NO' ) TRIM(Caller) // ': '
       END IF

       ! If there are several partitions to be saved than plot the partition too
       IF ( MaxOutputPE > 0 ) THEN
         WRITE( *,'(A,I0,A)', ADVANCE = 'NO' ) 'Part',OutputPE,': '
       END IF
     END IF


     IF ( nadv ) THEN
        WRITE( *,'(A)', ADVANCE = 'NO' )  TRIM(String)
     ELSE
        WRITE( *,'(A)', ADVANCE = 'YES' ) TRIM(String)
     END IF
     nadv1 = nadv

     CALL FLUSH(6)
!-----------------------------------------------------------------------
   END SUBROUTINE Info
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> May be used to skip computation that only relates to printing info.
!-----------------------------------------------------------------------
   FUNCTION InfoActive( Level ) RESULT( Show ) 
!-----------------------------------------------------------------------
     INTEGER, OPTIONAL :: Level
     LOGICAL :: Show
!-----------------------------------------------------------------------
     INTEGER, PARAMETER :: DefLevel = 4
!-----------------------------------------------------------------------

     IF ( PRESENT( Level ) ) THEN
       Show = OutputLevelMask(Level)
     ELSE
       Show = OutputLevelMask(DefLevel) 
     END IF

!-----------------------------------------------------------------------
   END FUNCTION InfoActive
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!> When a suspicious incident takes place this subroutine may be used
!> to inform the user.
!-----------------------------------------------------------------------
   SUBROUTINE Warn( Caller, String, noAdvance )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     LOGICAL, OPTIONAL :: noAdvance
!-----------------------------------------------------------------------

     LOGICAL :: nadv, nadv1 = .FALSE.
     SAVE nadv1

!-----------------------------------------------------------------------
     IF ( .NOT. OutputLevelMask(2) ) RETURN

     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance

     IF ( nadv ) THEN
        WRITE( *, '(A,A,A,A)', ADVANCE='NO' ) &
          'WARNING:: ', TRIM(Caller), ': ', TRIM(String)
     ELSE
        IF ( .NOT. nadv1 ) THEN
           WRITE( *, '(A,A,A,A)', ADVANCE='YES' ) &
             'WARNING:: ', TRIM(Caller), ': ', TRIM(String)
        ELSE
           WRITE( *, '(A)', ADVANCE='YES' ) TRIM(String)
        END IF
     END IF
     nadv1 = nadv
     CALL FLUSH(6)
!-----------------------------------------------------------------------
   END SUBROUTINE Warn
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!> This routine may be used to inform the user of an error.
!-----------------------------------------------------------------------
   SUBROUTINE Error( Caller, String, noAdvance )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     LOGICAL, OPTIONAL :: noAdvance
!-----------------------------------------------------------------------

     LOGICAL :: nadv, nadv1 = .FALSE.
     SAVE nadv1

!-----------------------------------------------------------------------
     IF ( .NOT. OutputLevelMask(1) ) RETURN

     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance

     IF ( nadv ) THEN
        WRITE( *, '(A,A,A,A)', ADVANCE='NO' ) &
          'ERROR:: ', TRIM(Caller), ': ', TRIM(String )
     ELSE
        IF ( .NOT. nadv1 ) THEN
           WRITE( *, '(A,A,A,A)', ADVANCE='YES' ) &
             'ERROR:: ', TRIM(Caller), ': ', TRIM(String)
        ELSE
           WRITE( *, '(A)', ADVANCE='YES' ) TRIM(String)
        END IF
     END IF
     nadv1 = nadv
     CALL FLUSH(6)
!-----------------------------------------------------------------------
   END SUBROUTINE Error
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> This routine may be used to terminate the program in the case of an error.
!-----------------------------------------------------------------------
   SUBROUTINE Fatal( Caller, String, noAdvance )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     LOGICAL, OPTIONAL :: noAdvance
!-----------------------------------------------------------------------

     LOGICAL :: nadv, nadv1 = .FALSE.
     SAVE nadv1

!-----------------------------------------------------------------------
     IF ( .NOT. OutputLevelMask(0) ) STOP

     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance

     IF ( nadv ) THEN
        WRITE( *, '(A,A,A,A)', ADVANCE='NO' ) &
          'ERROR:: ', TRIM(Caller), ': ', TRIM(String )
     ELSE
        IF ( .NOT. nadv1 ) THEN
           WRITE( *, '(A,A,A,A)', ADVANCE='YES' ) &
             'ERROR:: ', TRIM(Caller), ': ', TRIM(String)
        ELSE
           WRITE( *, '(A)', ADVANCE='YES' ) TRIM(String)
        END IF
        STOP
     END IF
     nadv1 = nadv
     CALL FLUSH(6)
!-----------------------------------------------------------------------
   END SUBROUTINE Fatal
!-----------------------------------------------------------------------

END MODULE Messages

!> \}
