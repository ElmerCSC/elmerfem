
!
! Subroutine to implement QMR iterative method (double complex)
!
! $Id: huti_qmr_C.F90,v 1.8 2005/06/02 15:35:27 vierinen Exp $



































#include  "huti_fdefs.h" 

!*************************************************************************
!
! This subroutine is based on a paper by Freund and Nachtigal:
! "An Implementation of the QMR Method Based on Coupled Two-Term
!  Recurrences", 1994 (SIAM J. Sci. Comput, March 1994)
! and a book by Barret et al.:
! "Templates for the Solution of Linear Systems: Building Blocks for
!  Iterative Methods", 1993.
! 
! All matrix-vector operations are done externally, so we do not need
! to know about the matrix structure (sparse or dense). Memory allocation
! for the working arrays has also been done externally.
!
!*************************************************************************
! Work array is used in the following order:
! work(:,1) = v
! work(:,2) = v tilde
! work(:,3) = y
! work(:,4) = y tilde
! work(:,5) = w
! work(:,6) = w tilde
! work(:,7) = z
! work(:,8) = z tilde
! work(:,9) = p
! work(:,10) = p tilde
! work(:,11) = q
! work(:,12) = d
! work(:,13) = s
! work(:,14) = r
!
!*************************************************************************
! Definitions to make the code more understandable and to make it look
! like the pseudo code
!
#define  X  xvec 
#define  B  rhsvec 

#define  V  work(:,1) 
#define  V_ind  1 
#define  VTLD  work(:,2) 
#define  VTLD_ind  2 
#define  Y  work(:,3) 
#define  Y_ind  3 
#define  YTLD  work(:,4) 
#define  YTLD_ind  4 
#define  W  work(:,5) 
#define  W_ind  5 
#define  WTLD  work(:,6) 
#define  WTLD_ind  6 
#define  Z  work(:,7) 
#define  Z_ind  7 
#define  ZTLD  work(:,8) 
#define  ZTLD_ind  8 
#define  P  work(:,9) 
#define  P_ind  9 
#define  PTLD  work(:,10) 
#define  PTLD_ind  10 
#define  Q  work(:,11) 
#define  Q_ind  11 
#define  D  work(:,12) 
#define  D_ind  12 
#define  S  work(:,13) 
#define  S_ind  13 
#define  R  work(:,14) 
#define  R_ind  14 
!*************************************************************************
  
!*************************************************************************
!*************************************************************************
! Complex version
!*************************************************************************
!*************************************************************************

subroutine  huti_cqmrsolv  ( ndim, wrkdim, xvec, rhsvec, &
                          ipar, dpar, work, matvecsubr, pcondlsubr, &
                          pcondrsubr, dotprodfun, normfun, stopcfun )

  implicit none

  ! Parameters

  external matvecsubr, pcondlsubr, pcondrsubr
  external dotprodfun, normfun, stopcfun
  complex :: dotprodfun
  real :: normfun
  real :: stopcfun

  integer :: ndim, wrkdim
  complex, dimension(ndim) :: xvec, rhsvec
  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
  complex, dimension(ndim,wrkdim) :: work

  ! Local variables

  complex :: beta, gamma, oldgamma, delta, rho, rhonext
  complex :: psi, theta, oldtheta, eta, epsilon
  integer iter_count
  real :: residual, rhsnorm, precrhsnorm

  !
  ! End of variable declarations
  !*********************************************************************

  !*********************************************************************
  ! The actual QMR begins here (look the pseudo code in the
  ! "Templates.."-book on page 24)
  !
  ! First the initialization part
  !

  oldgamma = 0
  oldtheta = 0
  iter_count = 1

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
     call  huti_crandvec   ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call matvecsubr( X, R, ipar )
  R = B - R
  VTLD = R
  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call pcondlsubr( Y, VTLD, ipar )

  WTLD = R
  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call pcondrsubr( Z, WTLD, ipar )

  rho = normfun( HUTI_NDIM, Y, 1 )
  psi = normfun( HUTI_NDIM, Z, 1 )
  oldgamma = 1
  eta = -1

  !
  ! This is where the loop starts (that is we continue from here after
  ! the first iteration)
  !

300 continue

  if (( rho .eq. 0 ) .or. ( psi .eq. 0 )) then
     HUTI_INFO = HUTI_QMR_RHO_PSI
     go to 1000
  end if

  V = VTLD / rho
  Y = Y / rho
  W = WTLD / psi
  Z = Z / psi
  delta = dotprodfun( HUTI_NDIM, Z, 1, Y, 1 )

  if ( delta .eq. 0 ) then
     HUTI_INFO = HUTI_QMR_DELTA
     go to 1000
  end if

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call pcondrsubr( YTLD, Y, ipar )

  HUTI_EXTOP_MATTYPE = HUTI_MAT_TRPSED
  call pcondlsubr( ZTLD, Z, ipar )

  if ( iter_count .eq. 1 ) then
     P = YTLD
     Q = ZTLD
  else
     P = YTLD - ((psi * delta)/epsilon) * P
     Q = ZTLD - ((rho * delta)/epsilon) * Q
  end if

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call matvecsubr( P, PTLD, ipar )
  epsilon = dotprodfun( HUTI_NDIM, Q, 1, PTLD, 1 )
  if ( epsilon .eq. 0 ) then
     HUTI_INFO = HUTI_QMR_EPSILON
     go to 1000
  end if

  beta = epsilon / delta
  if ( beta .eq. 0 ) then
     HUTI_INFO = HUTI_QMR_BETA
     go to 1000
  end if

  VTLD = PTLD - beta * V

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED
  call pcondlsubr( Y, VTLD, ipar )

  rhonext = normfun( ndim, Y, 1 )

  HUTI_EXTOP_MATTYPE = HUTI_MAT_TRPSED
  call matvecsubr( Q, WTLD, ipar )

  WTLD = WTLD - beta * W

  HUTI_EXTOP_MATTYPE = HUTI_MAT_TRPSED
  call pcondrsubr( Z, WTLD, ipar )

  psi = normfun( HUTI_NDIM, Z, 1 )
  theta = rhonext / (oldgamma * abs( beta ))
  gamma = 1 / sqrt( 1 + theta * theta )
  if ( gamma .eq. 0 ) then
     HUTI_INFO = HUTI_QMR_GAMMA
     go to 1000
  end if

  eta = -1 * ( eta * rho * gamma * gamma ) / ( beta * oldgamma * oldgamma )

  if ( iter_count .eq. 1 ) then
     D = eta * P
     S = eta * PTLD
  else
     D = eta * P + ( oldtheta * gamma ) * ( oldtheta * gamma ) * D
     S = eta * PTLD + ( oldtheta * gamma ) * ( oldtheta * gamma ) * S
  end if
  X = X + D
  R = R - S

  !
  ! Check the convergence against selected stopping criterion
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call matvecsubr( X, YTLD, ipar )
     YTLD = YTLD - B
     residual = normfun( HUTI_NDIM, YTLD, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call matvecsubr( X, YTLD, ipar )
     YTLD = YTLD - B
     residual = normfun( HUTI_NDIM, YTLD, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     residual = normfun( HUTI_NDIM, D, 1 )
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call matvecsubr( X, YTLD, ipar )
     YTLD = YTLD - B
     residual = normfun( HUTI_NDIM, YTLD, 1 )
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

  rho = rhonext
  oldgamma = gamma
  oldtheta = theta

  !
  ! Return back to the iteration loop (without initialization)
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

end subroutine  huti_cqmrsolv 

!*********************************************************************
