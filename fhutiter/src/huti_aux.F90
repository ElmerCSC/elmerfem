module huti_aux

  !
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
  ! **************************************************************************/

  ! Auxiliary routines for HUTI
  !

#include  "huti_fdefs.h" 

contains

  !*************************************************************************
  !*************************************************************************
  !
  ! This is a dummy preconditioning routine copying only the vector

  subroutine  huti_sdummy_pcondfun  ( u, v, ipar )

    implicit none

    real, dimension(*) :: u, v
    integer, dimension(*) :: ipar

    integer :: i

    !*************************************************************************

    do i = 1, HUTI_NDIM
       u(i) = v(i)
    end do

    return

  end subroutine  huti_sdummy_pcondfun

  !*************************************************************************

  !*************************************************************************
  !*************************************************************************
  !
  ! This routine fills a vector with pseudo random numbers

  subroutine  huti_srandvec  ( u, ipar )

    implicit none

    real, dimension(*) :: u
    integer, dimension(*) :: ipar

    integer :: i
    real :: harvest

    !*************************************************************************

    do i = 1, HUTI_NDIM
       call random_number( harvest )
       u(i) = harvest
    end do

    return

  end subroutine  huti_srandvec

  !*************************************************************************
  !*************************************************************************
  !
  ! Construct LU decompostion of the given matrix and solve LUu = v
  !

  subroutine  huti_slusolve  ( n, lumat, u, v )

    implicit none

    integer :: n
    real, dimension(n,n) :: lumat
    real, dimension(n) :: u, v

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

  end subroutine  huti_slusolve


  !*************************************************************************
  !*************************************************************************
  !
  ! This is a dummy preconditioning routine copying only the vector

  subroutine  huti_ddummy_pcondfun  ( u, v, ipar )

    implicit none

    double precision, dimension(*) :: u, v
    integer, dimension(*) :: ipar

    integer :: i

    !*************************************************************************
    !$OMP PARALLEL DO
    do i = 1, HUTI_NDIM
       u(i) = v(i)
    end do
    !$OMP END PARALLEL DO 

    return

  end subroutine  huti_ddummy_pcondfun

  !*************************************************************************

  !*************************************************************************
  !*************************************************************************
  !
  ! This routine fills a vector with pseudo random numbers

  subroutine  huti_drandvec  ( u, ipar )

    implicit none

    double precision, dimension(*) :: u
    integer, dimension(*) :: ipar

    integer :: i
    real :: harvest

    !*************************************************************************

    do i = 1, HUTI_NDIM
       call random_number( harvest )
       u(i) = harvest
    end do

    return

  end subroutine  huti_drandvec

  !*************************************************************************
  !*************************************************************************
  !
  ! Construct LU decompostion of the given matrix and solve LUu = v
  !

  subroutine  huti_dlusolve  ( n, lumat, u, v )

    implicit none

    integer :: n
    double precision, dimension(n,n) :: lumat
    double precision, dimension(n) :: u, v

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

  end subroutine  huti_dlusolve


  !*************************************************************************
  !*************************************************************************
  !
  ! This is a dummy preconditioning routine copying only the vector

  subroutine  huti_cdummy_pcondfun  ( u, v, ipar )

    implicit none

    complex, dimension(*) :: u, v
    integer, dimension(*) :: ipar

    integer :: i

    !*************************************************************************

    do i = 1, HUTI_NDIM
       u(i) = v(i)
    end do

    return

  end subroutine  huti_cdummy_pcondfun

  !*************************************************************************

  !*************************************************************************
  !*************************************************************************
  !
  ! This routine fills a vector with pseudo random numbers

  subroutine  huti_crandvec  ( u, ipar )

    implicit none

    complex, dimension(*) :: u
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

  end subroutine  huti_crandvec

  !*************************************************************************
  !*************************************************************************
  !
  ! Construct LU decompostion of the given matrix and solve LUu = v
  !

  subroutine  huti_clusolve  ( n, lumat, u, v )

    implicit none

    integer :: n
    complex, dimension(n,n) :: lumat
    complex, dimension(n) :: u, v

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

  end subroutine  huti_clusolve


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

end module huti_aux
