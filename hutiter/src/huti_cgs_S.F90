
!
! Subroutines to implement Conjugate Gradient Squared iteration
!
! $Id: huti_cgs_S.F90,v 1.8 2005/06/02 15:35:27 vierinen Exp $



































#include  "huti_fdefs.h" 

!*************************************************************************
!*************************************************************************
!
! This subroutine is based on a book by Barret et al.:
! "Templates for the Solution of Linear Systems: Building Blocks for
!  Iterative Methods", 1993.
!
! All matrix-vector operations are done externally, so we do not need
! to know about the matrix structure (sparse or dense). Memory allocation
! for the working arrays has also been done externally.

!*************************************************************************
! Work array is used in the following order:
! work(:,1) = r tilde (zero)
! work(:,2) = p
! work(:,3) = q
! work(:,4) = u
! work(:,5) = t1v (temporary)
! work(:,6) = t2v (temporary)
! work(:,7) = r
!
!*************************************************************************
! Definitions to make the code more understandable and to make it look
! like the pseudo code
!

#define  X  xvec 
#define  B  rhsvec 

#define  RTLD  work(:,1) 
#define  RTLD_ind  1 
#define  P  work(:,2) 
#define  P_ind  2 
#define  Q  work(:,3) 
#define  Q_ind  3 
#define  U  work(:,4) 
#define  U_ind  4 
#define  T1V  work(:,5) 
#define  T1V_ind  5 
#define  T2V  work(:,6) 
#define  T2V_ind  6 
#define  R  work(:,7) 
#define  R_ind  7 
  
!*************************************************************************
!*************************************************************************
! Single precision version
!*************************************************************************
!*************************************************************************

subroutine  huti_scgssolv  ( ndim, wrkdim, xvec, rhsvec, ipar,&
                            dpar, work, matvecsubr, pcondlsubr, pcondrsubr, &
                            dotprodfun, normfun, stopcfun )


  implicit none

  external matvecsubr, pcondlsubr, pcondrsubr
  external dotprodfun, normfun, stopcfun
  real :: dotprodfun
  real :: normfun
  real :: stopcfun

  ! Parameters

  integer :: ndim, wrkdim
  real, dimension(ndim) :: xvec, rhsvec
  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
  real, dimension(ndim,wrkdim) :: work

  ! Local variables

  real :: rho, oldrho, alpha, beta
  integer :: iter_count

  real :: residual, rhsnorm, precrhsnorm

  !
  ! End of variable declarations
  !*********************************************************************

  !*********************************************************************
  ! The actual CGS begins here (look the pseudo code in the
  ! "Templates..."-book, page 26)
  !
  ! First the initialization part
  !

  oldrho = 0
  iter_count = 1

  ! The following applies for all matrix operations in this solver

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED

  ! Norms of right-hand side vector are used in convergence tests

  if ( HUTI_STOPC .eq. HUTI_TRESID_SCALED_BYB .or. & 
       HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYB ) then
     rhsnorm = normfun( HUTI_NDIM, B, 1 )
  end if
  if ( HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYPRECB ) then
     call pcondlsubr( P, B, ipar )
     precrhsnorm = normfun( HUTI_NDIM, P, 1 )
  end if

  ! Generate vector X if needed

  if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
     call  huti_srandvec   ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  call matvecsubr( X, R, ipar )
  R = B - R
  RTLD = R

  !
  ! This is where the loop starts (that is we continue from here after
  ! the first iteration)
  !

300 continue

  rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
  if ( rho .eq. 0 ) then
     HUTI_INFO = HUTI_CGS_RHO
     go to 1000
  end if

  if ( iter_count .eq. 1 ) then
     U = R
     P = U
  else
     beta = rho / oldrho
     U = R + beta * Q
     P = U + beta * Q + beta * beta * P
  end if

  call pcondlsubr( T2V, P, ipar )
  call pcondrsubr( T1V, T2V, ipar )

  call matvecsubr( T1V, T2V, ipar )

  alpha = rho / dotprodfun( HUTI_NDIM, RTLD, 1, T2V, 1 )
  Q = U - alpha * T2V

  T2V = U + Q

  call pcondlsubr( U, T2V, ipar )
  call pcondrsubr( T1V, U, ipar )
  X = X + alpha * T1V

  call matvecsubr( T1V, T2V, ipar )
  R = R - alpha * T2V

  !
  ! Check the convergence against selected stopping criterion
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call matvecsubr( X, T1V, ipar )
     T1V = T1V - B
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call matvecsubr( X, T1V, ipar )
     T1V = T1V - B
     residual = normfun( HUTI_NDIM, T1V, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     T1V = alpha * T1V
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call matvecsubr( X, T1V, ipar )
     T1V = T1V - B
     residual = normfun( HUTI_NDIM, T1V, 1 )
  end select

  !
  ! Print debugging output if desired
  !

  if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
     if ( mod(iter_count, HUTI_DBUGLVL) .eq. 0 ) then
        write (*, '(I8, E11.4)') iter_count, residual
     end if
  end if

  if ( residual .lt. HUTI_TOLERANCE ) then
     HUTI_INFO = HUTI_CONVERGENCE
     go to 1000
  end if

  oldrho = rho

  !
  ! Return next time back to the iteration loop (without initialization)
  !

  iter_count = iter_count + 1
  if ( iter_count .gt. HUTI_MAXIT ) then
     HUTI_INFO = HUTI_MAXITER
     go to 1000
  end if

  go to 300

  !
  ! This is where we exit last time (after enough iterations or breakdown)
  !

1000 continue
  if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
     write (*, '(I8, E11.4)') iter_count, residual
  end if

  HUTI_ITERS = iter_count
  return

  ! End of execution
  !*********************************************************************

end subroutine  huti_scgssolv 

!*************************************************************************
