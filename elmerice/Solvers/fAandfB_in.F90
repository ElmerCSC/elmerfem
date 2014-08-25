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
! *
! *  Authors: Olivier Gagliardini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> Function a(D) and b(D) from Gagliardini and Meyssonier, 1997.
!> modified to fulfill the condition 3x a(D) >= 2x b(D) for D > 0.1
!> for n=3 only
FUNCTION ParameterA ( D ) RESULT(a)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   REAL(KIND=dp) :: a, D, DD     

    IF (D >= 1.0_dp) THEN
      a = 1.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.4 ) DD = 0.4_dp
      a = EXP( 13.22240_dp - 15.78652_dp * DD )

    ELSE
      a =  1.0_dp  + 2.0/3.0 * (1.0_dp - D) 
      a = a / ( D**1.5 )

    END IF
END FUNCTION ParameterA 


FUNCTION ParameterB ( D ) RESULT(b)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   REAL(KIND=dp) :: b, D, DD 

    IF (D >= 1.0_dp) THEN
      b = 0.0_dp           

    ELSE IF (D <= 0.81) THEN
      DD = D
      IF (D < 0.4 ) DD = 0.4_dp
      b = EXP( 15.09371_dp - 20.46489_dp * DD )

    ELSE
      b = (1.0_dp - D)**(1.0/3.0) 
      b = 3.0/4.0 * ( b / (3.0 * (1 - b)) )**1.5

    END IF
END FUNCTION ParameterB 
