! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st Apsril 1995 - , CSC - IT Center for Science Ltd., Finland
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
! *  along with this program (in fle fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 5 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *  
! *****************************************************************************/
!
!
!------------------------------------------------------------------------------
SUBROUTINE LinearSliding( Beta1, Beta2, n1,n2,n3, txz,tyz,p, uB1,uB2,uB3)
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  ! Compute basal velocity from stresses and sliding coefficients
  !
  !  REAL(KIND=dp) :: Beta1, Beta2
  !     INPUT: Slip coefficents
  !
  !  REAL(KIND=dp) :: n1, n2, n3
  !     INPUT: normal vectors
  !
  !  REAL(KIND=dp) :: txz, tyz, p
  !     INPUT: shear stresses, pressure
  !
  !  REAL(KIND=dp) :: uB1, uB2, uB3
  !     OUTPUT: basal velocity compoenents
  !
  !******************************************************************************

 
  !------------------------------------------------------------------------------
  ! This routine is called for each bedrock node
  !------------------------------------------------------------------------------

  USE ElementUtils

  IMPLICIT NONE

  REAL(KIND=dp) ::  Beta1,Beta2 
  REAL(KIND=dp) :: n1,n2,n3
  REAL(KIND=dp) :: txz,tyz,p
  REAL(KIND=dp) :: Tangent(3),Tangent2(3)
  REAL(KIND=dp) :: uB_t1,uB_t2,TauB_t1,TauB_t2
  REAL(KIND=dp), INTENT(OUT) :: uB1,uB2,uB3

 !------------------------------------------------------------------------------
 !    set the bottom velocity for current node, 
 !       input: beta, normal vec, nodenumber. return: uB
 !------------------------------------------------------------------------------

  CALL TangentDirections((/ n1, n2, n3/), Tangent, Tangent2 ) 


  TauB_t1 = Tangent(1)*(-n1*p+n3*txz)+Tangent(2)*(-n2*p+n3*tyz)+ &
       Tangent(3)*(n1*txz+n2*tyz-n3*p)

  TauB_t2 = Tangent2(1)*(-n1*p+n3*txz)+Tangent2(2)*(-n2*p+n3*tyz)+ &
       Tangent2(3)*(n1*txz+n2*tyz-n3*p)

  uB_t1=-TauB_t1/Beta1
  uB_t2=-TauB_t2/Beta2

  uB1=uB_t1*Tangent(1)+uB_t2*Tangent2(1)
  uB2=uB_t1*Tangent(2)+uB_t2*Tangent2(2)
  uB3=uB_t1*Tangent(3)+uB_t2*Tangent2(3)


  !------------------------------------------------------------------------------
END SUBROUTINE LinearSliding
!------------------------------------------------------------------------------

