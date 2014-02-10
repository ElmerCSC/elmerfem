
!
! Subroutine to implement BiConjugate Gradient Stabilised iteration
!
! $Id: huti_bicgstab_Z.F90,v 1.5 2005/05/04 09:57:33 vierinen Exp $



































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
! work(:,3) = p tilde
! work(:,4) = v
! work(:,5) = s
! work(:,6) = s tilde
! work(:,7) = t
! work(:,8) = r
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
#define  T1V  work(:,3) 
#define  T1V_ind  3 
#define  V  work(:,4) 
#define  V_ind  4 
#define  S  work(:,5) 
#define  S_ind  5 
#define  T2V  work(:,6) 
#define  T2V_ind  6 
#define  T  work(:,7) 
#define  T_ind  7 
#define  R  work(:,8) 
#define  R_ind  8 
  
!*************************************************************************
!*************************************************************************
! Double complex version
!*************************************************************************
!*************************************************************************

subroutine  huti_zbicgstabsolv  ( ndim, wrkdim, xvec, rhsvec, &
                            ipar, dpar, work, matvecsubr, pcondlsubr, &
                            pcondrsubr, dotprodfun, normfun, stopcfun )


  implicit none

  external matvecsubr, pcondlsubr, pcondrsubr
  external dotprodfun, normfun, stopcfun
  double complex :: dotprodfun
  double precision :: normfun
  double precision :: stopcfun

  ! Parameters

  integer :: ndim, wrkdim
  double complex, dimension(ndim) :: xvec, rhsvec
  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
  double complex, dimension(ndim,wrkdim) :: work

  ! Local variables

  double complex :: rho, oldrho, alpha, beta, omega
  integer :: iter_count

  double precision :: residual, rhsnorm, precrhsnorm

  !
  ! End of variable declarations
  !*********************************************************************

  !*********************************************************************
  ! The actual BiCGSTAB begins here (look the pseudo code in the
  ! "Templates..."-book, page 27)
  !
  ! First the initialization part
  !

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
     call  huti_zrandvec   ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  call matvecsubr( X, R, ipar )
  R = B - R
  RTLD = R
  P = 0; V = 0
  oldrho = 1; omega = 1; alpha = 0

  !
  ! This is where the loop starts (that is we continue from here after
  ! the first iteration)
  !

300 continue

  rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
  if ( rho .eq. 0 ) then
     HUTI_INFO = HUTI_BICGSTAB_RHO
     go to 1000
  end if

  beta = ( rho * alpha ) / ( oldrho * omega )
  P = R + beta * ( P - omega * V )

  call pcondlsubr( V, P, ipar )
  call pcondrsubr( T1V, V, ipar )
  call matvecsubr( T1V, V, ipar )

  alpha = rho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
  S = R - alpha * V

  residual = normfun( HUTI_NDIM, S, 1 )
  if ( residual .lt. HUTI_EPSILON ) then
     X = X + alpha * T1V
     HUTI_INFO = HUTI_BICGSTAB_SNORM
     go to 1000
  end if

  call pcondlsubr( T, S, ipar )
  call pcondrsubr( T2V, T, ipar )
  call matvecsubr( T2V, T, ipar )

  omega = ( dotprodfun( HUTI_NDIM, T, 1, S, 1 ) ) / &
          ( dotprodfun( HUTI_NDIM, T, 1, T, 1 ) )
  X = X + alpha * T1V + omega * T2V
  R = S - omega * T

  !
  ! Check the convergence against selected stopping criterion
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call matvecsubr( X, T2V, ipar )
     T1V = T2V - B
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call matvecsubr( X, T2V, ipar )
     T1V = T2V - B
     residual = normfun( HUTI_NDIM, T1V, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     T1V = alpha * T1V + omega * T2V
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call matvecsubr( X, T2V, ipar )
     T1V = T2V - B
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

  if ( omega .eq. 0 ) then
     HUTI_INFO = HUTI_BICGSTAB_OMEGA
     go to 1000
  end if

  oldrho = rho

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

end subroutine  huti_zbicgstabsolv 

!*************************************************************************
