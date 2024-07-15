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

  ! In Fortran standard 2003 this should work, else activate the 2nd
#if 1
   USE, INTRINSIC :: iso_fortran_env, ONLY : stdout=>output_unit
#else
#define stdout 6
#endif


#ifdef HAVE_XIOS
  USE XIOS
#endif

   IMPLICIT NONE
   
   CHARACTER(LEN=512) :: Message = ' '
   INTEGER, PRIVATE :: i
   LOGICAL :: OutputPrefix=.FALSE., OutputCaller=.TRUE.
   LOGICAL :: OutputLevelMask(0:31) = .TRUE.
   INTEGER :: MaxOutputLevel=31, MinOutputLevel=0, OutputPE = 0
   INTEGER :: MaxOutputThread=0, MaxThreads=1
   INTEGER :: MaxOutputPE = 0, MinOutputPE = 0
   INTEGER :: InfoOutUnit = stdout
   INTEGER, PARAMETER :: InfoToFileUnit = 33   
   LOGICAL :: InfoToFile = .FALSE.
   
   INTEGER, PARAMETER :: EXIT_OK=0, EXIT_ERROR=1

#ifdef HAVE_XIOS
   LOGICAL :: USE_XIOS = .FALSE. 
#endif


CONTAINS

!-----------------------------------------------------------------------
!> Prints information on the standard output if the requested or 
!> default output level does not surpass the maximum output level.
!-----------------------------------------------------------------------
   SUBROUTINE Info( Caller, String, noAdvance, Level )
!-----------------------------------------------------------------------
#define MEMDEBUG 0     
#if MEMDEBUG
     INTERFACE
       FUNCTION cpumemory() RESULT(dbl) BIND(C,name='cpumemory')
         USE, INTRINSIC :: ISO_C_BINDING
         REAL(C_DOUBLE) :: dbl
       END FUNCTION cpumemory
     END INTERFACE
     INTEGER(KIND=8) :: CurrUse
#endif
     
     CHARACTER(LEN=*) :: Caller, String
     INTEGER, OPTIONAL :: Level
     LOGICAL, OPTIONAL :: noAdvance
!-----------------------------------------------------------------------
     LOGICAL :: nadv, nadv1 = .FALSE.
     INTEGER :: n
     INTEGER, PARAMETER :: DefLevel = 4
     LOGICAL :: StdoutSet = .FALSE.
     INTEGER :: nthread, omp_get_thread_num

     
     SAVE nadv1

!-----------------------------------------------------------------------          
     IF ( OutputPE < 0 ) RETURN

     IF ( PRESENT( Level ) ) THEN
       if (Level > MaxOutputLevel) RETURN
       IF ( .NOT. OutputLevelMask(Level) ) RETURN
     ELSE
       ! The default level of info
       !-------------------------------------------
       IF( .NOT. OutputLevelMask(DefLevel) ) RETURN
     END IF

     IF(.NOT. StdoutSet ) THEN
       StdoutSet = .TRUE.
     END IF

     nthread = -1
     IF( MaxOutputThread > 0 ) THEN
       !$ nthread = omp_get_thread_num()+1
       IF(nthread > MaxOutputThread ) RETURN
     END IF
     
     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance
     
     IF(.NOT. nadv1 ) THEN
       IF ( OutputPrefix ) THEN
         WRITE( InfoOutUnit,'(A)', ADVANCE = 'NO' ) 'INFO:: '
       END IF

       IF ( OutputCaller ) THEN
         WRITE( InfoOutUnit,'(A)', ADVANCE = 'NO' ) TRIM(Caller) // ': '
       END IF
     END IF

     IF ( nadv ) THEN
       IF( MaxOutputPE > 0 .AND. .NOT. InfoToFile ) THEN
         IF( MaxOutputThread > 1 ) THEN
           WRITE( InfoOutUnit,'(A,I0,A,I0,A)', ADVANCE = 'NO' ) 'Part',OutputPE,' Thread',nthread,': '//TRIM(String)
         ELSE
           WRITE( InfoOutUnit,'(A,I0,A)', ADVANCE = 'NO' ) 'Part',OutputPE,': '//TRIM(String)
         END IF
       ELSE
         IF( MaxOutputThread > 1 ) THEN 
           WRITE( InfoOutUnit,'(A,I0,A)', ADVANCE = 'NO' ) 'Thread',nthread,': '//TRIM(String)
         ELSE
           WRITE( InfoOutUnit,'(A)', ADVANCE = 'NO' ) TRIM(String)
         END IF
       END IF
     ELSE
#if MEMDEBUG
       CurrUse = NINT( CPUMemory() ) 
       IF( MaxOutputPE > 0 .AND. .NOT. InfoToFile ) THEN
         WRITE( InfoOutUnit,'(A,I0,A,A,T50,A,I0)', ADVANCE = 'YES' ) 'Part',OutputPE,': ',TRIM(String), &
             'MEM: ',CurrUse
       ELSE
         WRITE( InfoOutUnit,'(A,T50,A,I0)', ADVANCE = 'YES' ) TRIM(String),'MEM: ',CurrUse
       END IF
#else
       IF( MaxOutputPE > 0 .AND. .NOT. InfoToFile ) THEN
         IF( MaxOutputThread > 1 ) THEN
           WRITE( InfoOutUnit,'(A,I0,A,I0,A)', ADVANCE = 'YES' ) 'Part',OutputPE,' Thread',nthread,': '//TRIM(String)
         ELSE
           WRITE( InfoOutUnit,'(A,I0,A)', ADVANCE = 'YES' ) 'Part',OutputPE,': '//TRIM(String)
         END IF
       ELSE
         IF( MaxOutputThread > 1 ) THEN 
           WRITE( InfoOutUnit,'(A,I0,A)', ADVANCE = 'YES' ) 'Thread',nthread,': '//TRIM(String)
         ELSE
           WRITE( InfoOutUnit,'(A)', ADVANCE = 'YES' ) TRIM(String)
         END IF
       END IF
#endif
     END IF
     nadv1 = nadv

     CALL FLUSH(InfoOutUnit)

          
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
       IF ( MaxOutputPE > 0 ) THEN
         WRITE( InfoOutUnit, '(A,A,A,I0,A,A)', ADVANCE='NO' ) &
             'WARNING:: ', TRIM(Caller), ': Part',OutputPE,':', TRIM(String)
       ELSE
         WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='NO' ) &
             'WARNING:: ', TRIM(Caller), ': ', TRIM(String)
       END IF
     ELSE
       IF ( .NOT. nadv1 ) THEN
         IF( MaxOutputPE > 0 ) THEN
           WRITE( InfoOutUnit, '(A,A,A,I0,A,A)', ADVANCE='YES' ) &
               'WARNING:: ', TRIM(Caller), ': Part',OutputPE,':', TRIM(String)
         ELSE
           WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='YES' ) &
               'WARNING:: ', TRIM(Caller), ': ', TRIM(String)
         END IF
       ELSE
         WRITE( InfoOutUnit, '(A)', ADVANCE='YES' ) TRIM(String)
       END IF
     END IF
     nadv1 = nadv
     CALL FLUSH(InfoOutUnit)
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
        WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='NO' ) &
          'ERROR:: ', TRIM(Caller), ': ', TRIM(String )
     ELSE
        IF ( .NOT. nadv1 ) THEN
           WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='YES' ) &
             'ERROR:: ', TRIM(Caller), ': ', TRIM(String)
        ELSE
           WRITE( InfoOutUnit, '(A)', ADVANCE='YES' ) TRIM(String)
        END IF
     END IF
     nadv1 = nadv
     CALL FLUSH(InfoOutUnit)
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

     IF ( .NOT. OutputLevelMask(0) ) STOP EXIT_ERROR

     nadv = .FALSE.
     IF ( PRESENT( noAdvance ) ) nadv = noAdvance

     IF ( nadv ) THEN
        WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='NO' ) &
          'ERROR:: ', TRIM(Caller), ': ', TRIM(String )
     ELSE
        IF ( .NOT. nadv1 ) THEN
           WRITE( InfoOutUnit, '(A,A,A,A)', ADVANCE='YES' ) &
             'ERROR:: ', TRIM(Caller), ': ', TRIM(String)
        ELSE
           WRITE( InfoOutUnit, '(A)', ADVANCE='YES' ) TRIM(String)
        END IF
        STOP EXIT_ERROR
     END IF
     nadv1 = nadv
     CALL FLUSH(InfoOutUnit)

#ifdef HAVE_XIOS
     IF (USE_XIOS) THEN
       CALL xios_context_finalize()
       CALL xios_finalize()
     ENDIF
#endif 
!-----------------------------------------------------------------------
   END SUBROUTINE Fatal
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> This routine may be used to terminate the program in the case of an error.
!-----------------------------------------------------------------------
   SUBROUTINE Assert(Condition, Caller, ErrorMessage)
!-----------------------------------------------------------------------
     CHARACTER(LEN=*), OPTIONAL :: Caller, ErrorMessage
     LOGICAL :: Condition
!-----------------------------------------------------------------------
     IF ( .NOT. OutputLevelMask(0) ) STOP EXIT_ERROR

     IF(Condition) RETURN !Assertion passed

     WRITE( Message, '(A)') 'ASSERTION ERROR'

     IF(PRESENT(Caller)) THEN
       WRITE( Message, '(A,A,A)') TRIM(Message),': ',TRIM(Caller)
     END IF

     IF(PRESENT(ErrorMessage)) THEN
       WRITE( Message, '(A,A,A)') TRIM(Message),': ',TRIM(ErrorMessage)
     END IF

     WRITE( *, '(A)', ADVANCE='YES' ) Message

     !Provide a stack trace if no caller info provided
#ifdef __GFORTRAN__
     IF(.NOT.PRESENT(Caller)) CALL BACKTRACE
#endif

     STOP EXIT_ERROR
!-----------------------------------------------------------------------
   END SUBROUTINE Assert
!-----------------------------------------------------------------------

   
END MODULE Messages

!> \}
