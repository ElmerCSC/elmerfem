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

! Auxiliary routines for HUTI
!
! $Id: huti_aux.src,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $

#include "huti_intdefs.h"
MAKE_INCLUDE(#,include, "huti_fdefs.h" )

!*************************************************************************
!*************************************************************************
!
! This is a dummy preconditioning routine copying only the vector

subroutine MAKE_SUBRN( huti_,dummy_pcondfun ) ( u, v, ipar )

  implicit none

  F_PRECISION_TYPE, dimension(*) :: u, v
  integer, dimension(*) :: ipar

  integer :: i

  !*************************************************************************
  
  do i = 1, HUTI_NDIM
     u(i) = v(i)
  end do

  return

end subroutine MAKE_SUBRN( huti_,dummy_pcondfun )

!*************************************************************************

!*************************************************************************
!*************************************************************************
!
! This routine fills a vector with pseudo random numbers

subroutine MAKE_SUBRN( huti_,randvec ) ( u, ipar )

  implicit none

  F_PRECISION_TYPE, dimension(*) :: u
  integer, dimension(*) :: ipar

  integer :: i
  real :: harvest

  !*************************************************************************

  do i = 1, HUTI_NDIM
     call random_number( harvest )
     u(i) = harvest
#if defined(C_PRE) || defined(Z_PRE)
     call random_number( harvest )
     u(i) = cmplx( 0, harvest )
#endif
  end do

  return

end subroutine MAKE_SUBRN( huti_,randvec )

!*************************************************************************
!*************************************************************************
!
! Construct LU decompostion of the given matrix and solve LUu = v
!

subroutine MAKE_SUBRN( huti_,lusolve ) ( n, lumat, u, v )

  implicit none

  integer :: n
  F_PRECISION_TYPE, dimension(n,n) :: lumat
  F_PRECISION_TYPE, dimension(n) :: u, v

  integer :: i, j, k

  !*************************************************************************
  
  !
  ! This is from Saad''s book, Algorithm 10.4
  !

  do i = 2,n
     do k = 1, i - 1

	! Check for small pivot

	if ( abs(lumat(k,k)) .lt. 1.0e-16 ) then
	   print *, '(libhuti.a) GMRES: small pivot', lumat(k,k)
	end if

	! Compute a_ik = a_ik / a_kk

	lumat(i,k) = lumat(i,k) / lumat(k,k)

	do j = k + 1, n

	   ! Compute a_ij = a_ij - a_ik * a_kj

	   lumat(i,j) = lumat(i,j) - lumat(i,k) * lumat(k,j)
	end do

     end do
  end do

  ! Forward solve, Lu = v

  do i = 1, n

     ! Compute u(i) = v(i) - sum L(i,j) u(j)

     u(i) = v(i)
     do k = 1, i - 1
        u(i) = u(i) - lumat(i,k) * u(k)
     end do
  end do

  ! Backward solve, u = inv(U) u

  do i = n, 1, -1

     ! Compute u(i) = u(i) - sum U(i,j) u(j)

     do k = i + 1, n
	u(i) = u(i) - lumat(i,k) * u(k)
     end do

     ! Compute u(i) = u(i) / U(i,i)

     u(i) = u(i) / lumat(i,i) 
  end do

  return

end subroutine MAKE_SUBRN( huti_,lusolve )

