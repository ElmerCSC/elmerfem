#include "huti_fdefs.h"

!******************************************************************
!******************************************************************
!
! Purpose of this example is to show how to use HUTIter library
!

program hutiexample

  use globals_module
  implicit none

  external huti_d_tfqmr, huti_d_cg
  external own_matvec

  ! local variables

  character (len=64) :: filename
  integer :: ndim, i, j
  double precision :: value

  ! Variables needed for HUTI

  double precision, dimension(:,:), allocatable, target :: A
  double precision, dimension(:), allocatable :: b
  double precision, dimension(:), allocatable :: x
  double precision, dimension(:,:), allocatable :: work

  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar

  integer :: workdim
  parameter ( workdim = HUTI_TFQMR_WORKSIZE )

  !******************************************************************

  write (6, *) 'Enter name of the matrix file :'
  read (5, *) filename
  write (6, *) 'Enter leading dimension of the matrix :'
  read (5, *) ndim

  allocate( A(ndim, ndim) )
  A_ptr => A

  allocate( b(ndim) )
  allocate( x(ndim) )
  allocate( work(ndim, workdim) )

  A = 0; b = 1

  open ( 11, file=filename, status='old', err=20 )

10 continue
  read ( 11, *, end=20 ) i, j, value
  A(i, j) = value
  go to 10

20 close(11)

  HUTI_NDIM = ndim
  HUTI_WRKDIM = workdim
  HUTI_DBUGLVL = HUTI_ITEROUTPUT
  HUTI_MAXIT = HUTI_DFLTMAXIT
  HUTI_TOLERANCE = HUTI_DFLTTOLERANCE

  call HUTI_D_TFQMR( x, b, ipar, dpar, work, own_matvec, 0, 0, 0, 0, 0 )

  print *, 'HUTI returned ', HUTI_INFO

  open( 11, file=trim( filename ) // '.out' )
  write( 11, '(I6,A,E25.17)') ( i,' ', x(i), i=1,ndim )
  close( 11 )

end program hutiexample

!******************************************************************
!******************************************************************
!
! This is an example of the user supplied matvec routine for HUTI
! solvers.
!

subroutine own_matvec ( u, v, ipar )

  use globals_module
  implicit none

  ! Parameters
  double precision, dimension(*) :: u, v
  integer, dimension(*) :: ipar

  ! Local variables
  integer :: i,j
  double precision :: dsum

  do i = 1,HUTI_NDIM
     dsum = 0
     do j = 1,HUTI_NDIM
	dsum = dsum + A_ptr(i,j) * u(j)
     end do
     v(i) = dsum
  end do

end subroutine own_matvec

!******************************************************************

