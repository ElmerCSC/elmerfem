
! Auxiliary routines for HUTI
!
! $Id: huti_aux_Z.F90,v 1.5 2005/05/04 09:57:27 vierinen Exp $



































#include  "huti_fdefs.h" 

!*************************************************************************
!*************************************************************************
!
! This is a dummy preconditioning routine copying only the vector

subroutine  huti_zdummy_pcondfun  ( u, v, ipar )

  implicit none

  double complex, dimension(*) :: u, v
  integer, dimension(*) :: ipar

  integer :: i

  !*************************************************************************
  
  do i = 1, HUTI_NDIM
     u(i) = v(i)
  end do

  return

end subroutine  huti_zdummy_pcondfun 

!*************************************************************************

!*************************************************************************
!*************************************************************************
!
! This routine fills a vector with pseudo random numbers

subroutine  huti_zrandvec  ( u, ipar )

  implicit none

  double complex, dimension(*) :: u
  integer, dimension(*) :: ipar

  integer :: i
  real :: harvest

  !*************************************************************************

  do i = 1, HUTI_NDIM
     call random_number( harvest )
     u(i) = harvest

     call random_number( harvest )
     u(i) = cmplx( 0, harvest )

  end do

  return

end subroutine  huti_zrandvec 

!*************************************************************************
!*************************************************************************
!
! Construct LU decompostion of the given matrix and solve LUu = v
!

subroutine  huti_zlusolve  ( n, lumat, u, v )

  implicit none

  integer :: n
  double complex, dimension(n,n) :: lumat
  double complex, dimension(n) :: u, v

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

end subroutine  huti_zlusolve 

