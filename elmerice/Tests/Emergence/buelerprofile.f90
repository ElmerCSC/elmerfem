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
! ******************************************************************************
FUNCTION Bueler ( Model, nodenumber, coord) RESULT(elevation)
   USE Types
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber
   REAL(KIND=dp) :: coord(3), elevation
  
   REAL(KIND=dp), PARAMETER :: L = 361.25d03, H0 = 2.06d03, n = 3.0_dp, minh=100.0_dp
   REAL(KIND=dp) :: T1, T2, H, R

   R = SQRT(coord(1)*coord(1) + coord(2)*coord(2))
   T1 = ( H0/((n - 1.0_dp)**(n/(2.0_dp*n + 2.0_dp))) )
   T2 = (n + 1.0_dp)*(R/L) - n*((R/L)**((n + 1.0_dp)/n)) + &
        n*((1.0_dp - (R/L))**((n + 1.0_dp)/n) ) - 1.0
   H =  T1 * ( T2**(n/(2.0*n + 2.0)) )
   elevation = MAX(H,minh)
 END FUNCTION Bueler
