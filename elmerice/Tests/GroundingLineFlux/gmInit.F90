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
! A simple user function to initialise the grouned mask
!  mask=0 for  x=0
!  mask=-1 for x>0
!  mask=+1 for x<0
FUNCTION gmInit( Model, nodenumber, x) RESULT(gm)
  USE DefUtils
  IMPLICIT NONE

  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: x,gm
  INTEGER       :: nodenumber
  REAL(KIND=dp),PARAMETER :: EEPS=1.0e-6


  gm=0._dp
  IF (x.gt.EEPS) THEN
    gm=-1.0_dp
  ELSE IF (x.lt.-EEPS) THEN
    gm=+1.0_dp
  END IF

  End
