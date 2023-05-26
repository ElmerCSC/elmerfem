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

   REAL(KIND=dp) :: CT, RT
   INTEGER, PARAMETER :: Initialize=0
   INTEGER :: tlen
   LOGICAL :: Silent
   CHARACTER(:), ALLOCATABLE :: DateStr
   CHARACTER(LEN=MAX_NAME_LEN) :: toutput

   CALL envir( 'ELMERSOLVER_OUTPUT_TOTAL', toutput, tlen )
   Silent = toutput(1:1)=='0' .OR. toutput(1:5)=='false'

   CT = CPUtime()
   RT = RealTime()

   IF ( .NOT. Silent ) THEN
     DateStr = FormatDate()
     WRITE( *,'(A,A)' ) "ELMER SOLVER (v " // VERSION // ") STARTED AT: ", TRIM(DateStr)
     CALL FLUSH(6)
   END IF

   CALL ElmerSolver(Initialize)

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
