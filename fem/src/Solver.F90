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

!> \ingroup ElmerLib

!------------------------------------------------------------------------------
!> A caller for the Elmer main program.
!------------------------------------------------------------------------------


#include "../config.h"

PROGRAM Solver
   USE Types
   USE GeneralUtils
   USE ParallelUtils

   IMPLICIT NONE

   REAL(KIND=dp) :: CT, RT
   INTEGER, PARAMETER :: Initialize=0
   INTEGER :: tlen
   LOGICAL :: Silent
   CHARACTER(:), ALLOCATABLE :: DateStr
   CHARACTER(LEN=MAX_NAME_LEN) :: toutput

   INTEGER :: iargc, nargs, arglen
   CHARACTER(:), ALLOCATABLE :: buf
   TYPE(ArgStr_t), ALLOCATABLE :: args(:)

   INTERFACE
     SUBROUTINE ElmerSolver(initialize, args, NoArgs)
       USE Types
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: initialize
       INTEGER, INTENT(IN) :: NoArgs
       TYPE(ArgStr_t), INTENT(IN) :: args(:)
     END SUBROUTINE ElmerSolver
   END INTERFACE

   CALL envir( 'ELMERSOLVER_OUTPUT_TOTAL', toutput, tlen )
   Silent = toutput(1:1)=='0' .OR. toutput(1:5)=='false'

   CT = CPUtime()
   RT = RealTime()

   IF ( .NOT. Silent ) THEN
     DateStr = FormatDate()
     WRITE( *,'(A,A)' ) "ELMER SOLVER (v " // ELMER_FEM_VERSION // ") STARTED AT: ", TRIM(DateStr)
     CALL FLUSH(6)
   END IF

   ! Get number of command line arguments
   nargs = COMMAND_ARGUMENT_COUNT()

   ! Collect command line arguments
   IF( nargs > 0 ) THEN 
     ALLOCATE(args(nargs))
     ALLOCATE(CHARACTER(MAX_PATH_LEN)::buf)

     iargc = 0
     DO WHILE( iargc < nargs )
       iargc = iargc + 1 
       CALL GET_COMMAND_ARGUMENT(iargc, buf, length=arglen)
       args(iargc) % astr = buf(1:arglen)
     END DO

     DEALLOCATE(buf)
   END IF

   CALL ElmerSolver(Initialize, args, nargs)

   IF ( .NOT. Silent ) THEN
     IF ( ParEnv % myPE == 0 ) THEN
       WRITE( *,'(a,F12.2,F12.2)' ) 'SOLVER TOTAL TIME(CPU,REAL): ', &
                   CPUTime()-CT, RealTime()-RT
       DateStr = FormatDate()
       WRITE( *,'(A,A)' ) 'ELMER SOLVER FINISHED AT: ', TRIM(DateStr)
       CALL FLUSH(6)
     END IF
   END IF
   
END PROGRAM Solver

! ******************************************************************************

!> \}
