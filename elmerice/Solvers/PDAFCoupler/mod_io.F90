!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! *****************************************************************************
! * COUPLING WITH PDAF (Parallel Data Assimilation Framework)
! *  https://pdaf.awi.de/trac/wiki/WikiStart
! *****************************************************************************
! * Generic routines for printing info and IO
MODULE mod_IO
    USE, INTRINSIC :: iso_fortran_env, ONLY : stdout=>output_unit
    USE netcdf
    USE mod_assimilation, ONLY : screen
    USE PDAF, ONLY : PDAF_abort

    IMPLICIT NONE
    INTEGER :: OutUnit = stdout

CONTAINS

!**********************************************************************
!* Print message from Caller to stdout (as all PDAF messages)       ***
!**********************************************************************        
 SUBROUTINE LOCAL_INFO(Caller, String, Level)
   !-----------------------------------------------------------------------
   CHARACTER(LEN=*) :: Caller, String
   INTEGER :: Level
   !-----------------------------------------------------------------------

   IF (Level > screen) RETURN

   WRITE( OutUnit , '(A,A,A,A)') &
          'PDAF_ELMER:: ', TRIM(Caller), ': ', TRIM(String )

   CALL FLUSH(OutUnit)

 END SUBROUTINE LOCAL_INFO

!**********************************************************************
!* Print message from caller and call PDAF_abort                    ***
!**********************************************************************  
 SUBROUTINE PDAF_FATAL(Caller, String)
   !-----------------------------------------------------------------------
   CHARACTER(LEN=*) :: Caller, String
   !-----------------------------------------------------------------------

   WRITE( OutUnit , '(A,A,A,A)') &
          'PDAF_ELMER:: ERROR:: ', TRIM(Caller), ': ', TRIM(String )
   CALL FLUSH(OutUnit)
   CALL PDAF_abort(1)

 END SUBROUTINE PDAF_FATAL

!**********************************************************************
!* Generic status check for netcf functions                         ***
!*   Write error and abort if the function return an error          ***
!**********************************************************************
 SUBROUTINE NFCHECK(status)
    ! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then
       WRITE( OutUnit , *) trim(nf90_strerror(status))
       CALL PDAF_abort(1)
    end if

  END SUBROUTINE NFCHECK

END MODULE mod_IO  
