module huti_bicgstab_2
  use huti_aux
  implicit none
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
! *****************************************************************************/

!
! Subroutine to implement BiConjugate Gradient Stabilised (2) iteration

#include  "huti_fdefs.h" 

!*************************************************************************
!*************************************************************************
!
! This subroutine is based on a paper by Henk A. Van der Vorst:
! "Parallel Iterative Solution Methods for Linear Systems arising from
!  Discretized PDE's". This is the Bi-CGSTAB(2) version.
!
! All matrix-vector operations are done externally, so we do not need
! to know about the matrix structure (sparse or dense). Memory allocation
! for the working arrays has also been done externally.

!*************************************************************************
! Work array is used in the following order:
! work(:,1) = r tilde (zero)
! work(:,2) = u
! work(:,3) = t1v
! work(:,4) = v
! work(:,5) = s
! work(:,6) = w
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
#define  U  work(:,2) 
#define  U_ind  2 
#define  T1V  work(:,3) 
#define  T1V_ind  3 
#define  V  work(:,4) 
#define  V_ind  4 
#define  S  work(:,5) 
#define  S_ind  5 
#define  W  work(:,6) 
#define  W_ind  6 
#define  T  work(:,7) 
#define  T_ind  7 
#define  R  work(:,8) 
#define  R_ind  8 

!*************************************************************************

contains

  !*************************************************************************
  !*************************************************************************
  ! Single precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_sbicgstab_2solv  ( ndim, wrkdim, xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )

    use huti_interfaces
    implicit none

    procedure( mv_iface_s ), pointer :: matvecsubr
    procedure( pc_iface_s ), pointer :: pcondlsubr 
    procedure( pc_iface_s ), pointer :: pcondrsubr 
    procedure( dotp_iface_s ), pointer :: dotprodfun 
    procedure( norm_iface_s ), pointer :: normfun 
    procedure( stopc_iface_s ), pointer :: stopcfun 

    ! Parameters

    integer :: ndim, wrkdim
    real, dimension(ndim) :: xvec, rhsvec
    integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
    double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
    real, dimension(ndim,wrkdim) :: work

    ! Local variables

    real :: rho, oldrho, alpha, beta, omega1, omega2
    real :: tau, delta, myy
    integer :: iter_count

    real :: residual, rhsnorm, precrhsnorm

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
       call pcondlsubr( T1V, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, T1V, 1 )
    end if

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_srandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
       X = 1
    end if

    call pcondrsubr( U, X, ipar )
    call matvecsubr( U, R, ipar )
    U = B - R
    call pcondlsubr( R, U, ipar )
    RTLD = R
    U = 0
    oldrho = 1; omega2 = 1; alpha = 0

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    oldrho = -omega2 * oldrho

    !
    ! This is the even BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    U = R - beta * U

    call pcondrsubr( V, U, ipar )
    call matvecsubr( V, T1V, ipar )
    call pcondlsubr( V, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
    R = R - alpha * V

    call pcondrsubr( S, R, ipar )
    call matvecsubr( S, T1V, ipar )
    call pcondlsubr( S, T1V, ipar )

    X = X + alpha * U

    !
    ! This is the odd BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, S, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    V = S - beta * V

    call pcondrsubr( W, V, ipar )
    call matvecsubr( W, T1V, ipar )
    call pcondlsubr( W, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, W, 1 )
    U = R - beta * U
    R = R - alpha * V
    S = S - alpha * W

    call pcondrsubr( T, S, ipar )
    call matvecsubr( T, T1V, ipar )
    call pcondlsubr( T, T1V, ipar )

    !
    ! This is the GCR(2) part
    !

    omega1 = dotprodfun( HUTI_NDIM, R, 1, S, 1 )
    myy = dotprodfun( HUTI_NDIM, S, 1, S, 1 )
    delta = dotprodfun( HUTI_NDIM, S, 1, T, 1 )
    tau = dotprodfun( HUTI_NDIM, T, 1, T, 1 )
    omega2 = dotprodfun( HUTI_NDIM, R, 1, T, 1 )

    tau = tau - ( delta * delta ) / myy
    omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau
    omega1 = ( omega1 - delta * omega2 ) / myy

    X = X + omega1 * R + omega2 * S + alpha * U
    R = R - omega1 * S - omega2 * T

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
       T1V = omega1 * R + omega2 * S + alpha * U
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
    U = U - omega1 * V - omega2 * W

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

  end subroutine  huti_sbicgstab_2solv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Double precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_dbicgstab_2solv  ( ndim, wrkdim, xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )

    use huti_interfaces
    implicit none

    procedure( mv_iface_d ), pointer :: matvecsubr 
    procedure( pc_iface_d ), pointer :: pcondlsubr 
    procedure( pc_iface_d ), pointer :: pcondrsubr 
    procedure( dotp_iface_d ), pointer :: dotprodfun 
    procedure( norm_iface_d ), pointer :: normfun 
    procedure( stopc_iface_d ), pointer :: stopcfun 

    ! Parameters

    integer :: ndim, wrkdim
    double precision, dimension(ndim) :: xvec, rhsvec
    integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
    double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
    double precision, dimension(ndim,wrkdim) :: work

    ! Local variables

    double precision :: rho, oldrho, alpha, beta, omega1, omega2
    double precision :: tau, delta, myy
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
       call pcondlsubr( T1V, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, T1V, 1 )
    end if

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_drandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
       X = 1
    end if

    call pcondrsubr( U, X, ipar )
    call matvecsubr( U, R, ipar )
    U = B - R
    call pcondlsubr( R, U, ipar )
    RTLD = R
    U = 0
    oldrho = 1; omega2 = 1; alpha = 0

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    oldrho = -omega2 * oldrho

    !
    ! This is the even BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    U = R - beta * U

    call pcondrsubr( V, U, ipar )
    call matvecsubr( V, T1V, ipar )
    call pcondlsubr( V, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
    R = R - alpha * V

    call pcondrsubr( S, R, ipar )
    call matvecsubr( S, T1V, ipar )
    call pcondlsubr( S, T1V, ipar )

    X = X + alpha * U

    !
    ! This is the odd BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, S, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    V = S - beta * V

    call pcondrsubr( W, V, ipar )
    call matvecsubr( W, T1V, ipar )
    call pcondlsubr( W, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, W, 1 )
    U = R - beta * U
    R = R - alpha * V
    S = S - alpha * W

    call pcondrsubr( T, S, ipar )
    call matvecsubr( T, T1V, ipar )
    call pcondlsubr( T, T1V, ipar )

    !
    ! This is the GCR(2) part
    !

    omega1 = dotprodfun( HUTI_NDIM, R, 1, S, 1 )
    myy = dotprodfun( HUTI_NDIM, S, 1, S, 1 )
    delta = dotprodfun( HUTI_NDIM, S, 1, T, 1 )
    tau = dotprodfun( HUTI_NDIM, T, 1, T, 1 )
    omega2 = dotprodfun( HUTI_NDIM, R, 1, T, 1 )

    tau = tau - ( delta * delta ) / myy
    omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau
    omega1 = ( omega1 - delta * omega2 ) / myy

    X = X + omega1 * R + omega2 * S + alpha * U
    R = R - omega1 * S - omega2 * T

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
       T1V = omega1 * R + omega2 * S + alpha * U
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
    U = U - omega1 * V - omega2 * W

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

  end subroutine  huti_dbicgstab_2solv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_cbicgstab_2solv  ( ndim, wrkdim, xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )

    use huti_interfaces
    implicit none

    procedure( mv_iface_c ), pointer :: matvecsubr 
    procedure( pc_iface_c ), pointer :: pcondlsubr 
    procedure( pc_iface_c ), pointer :: pcondrsubr 
    procedure( dotp_iface_c ), pointer :: dotprodfun 
    procedure( norm_iface_c ), pointer :: normfun 
    procedure( stopc_iface_c ), pointer :: stopcfun 

    ! Parameters

    integer :: ndim, wrkdim
    complex, dimension(ndim) :: xvec, rhsvec
    integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
    double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
    complex, dimension(ndim,wrkdim) :: work

    ! Local variables

    complex :: rho, oldrho, alpha, beta, omega1, omega2
    complex :: tau, delta, myy
    integer :: iter_count

    real :: residual, rhsnorm, precrhsnorm

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
       call pcondlsubr( T1V, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, T1V, 1 )
    end if

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_crandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
       X = 1
    end if

    call pcondrsubr( U, X, ipar )
    call matvecsubr( U, R, ipar )
    U = B - R
    call pcondlsubr( R, U, ipar )
    RTLD = R
    U = 0
    oldrho = 1; omega2 = 1; alpha = 0

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    oldrho = -omega2 * oldrho

    !
    ! This is the even BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    U = R - beta * U

    call pcondrsubr( V, U, ipar )
    call matvecsubr( V, T1V, ipar )
    call pcondlsubr( V, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
    R = R - alpha * V

    call pcondrsubr( S, R, ipar )
    call matvecsubr( S, T1V, ipar )
    call pcondlsubr( S, T1V, ipar )

    X = X + alpha * U

    !
    ! This is the odd BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, S, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    V = S - beta * V

    call pcondrsubr( W, V, ipar )
    call matvecsubr( W, T1V, ipar )
    call pcondlsubr( W, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, W, 1 )
    U = R - beta * U
    R = R - alpha * V
    S = S - alpha * W

    call pcondrsubr( T, S, ipar )
    call matvecsubr( T, T1V, ipar )
    call pcondlsubr( T, T1V, ipar )

    !
    ! This is the GCR(2) part
    !

    omega1 = dotprodfun( HUTI_NDIM, R, 1, S, 1 )
    myy = dotprodfun( HUTI_NDIM, S, 1, S, 1 )
    delta = dotprodfun( HUTI_NDIM, S, 1, T, 1 )
    tau = dotprodfun( HUTI_NDIM, T, 1, T, 1 )
    omega2 = dotprodfun( HUTI_NDIM, R, 1, T, 1 )

    tau = tau - ( delta * delta ) / myy
    omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau
    omega1 = ( omega1 - delta * omega2 ) / myy

    X = X + omega1 * R + omega2 * S + alpha * U
    R = R - omega1 * S - omega2 * T

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
       T1V = omega1 * R + omega2 * S + alpha * U
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF


    
    U = U - omega1 * V - omega2 * W

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

  end subroutine  huti_cbicgstab_2solv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Double complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_zbicgstab_2solv  ( ndim, wrkdim, xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )

    use huti_interfaces
    implicit none

    procedure( mv_iface_z ), pointer :: matvecsubr 
    procedure( pc_iface_z ), pointer :: pcondlsubr
    procedure( pc_iface_z ), pointer :: pcondrsubr 
    procedure( dotp_iface_z ), pointer :: dotprodfun 
    procedure( norm_iface_z ), pointer :: normfun 
    procedure( stopc_iface_z ), pointer :: stopcfun 

    ! Parameters

    integer :: ndim, wrkdim
    double complex, dimension(ndim) :: xvec, rhsvec
    integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
    double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
    double complex, dimension(ndim,wrkdim) :: work

    ! Local variables

    double complex :: rho, oldrho, alpha, beta, omega1, omega2
    double complex :: tau, delta, myy
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
       call pcondlsubr( T1V, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, T1V, 1 )
    end if

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_zrandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
       X = 1
    end if

    call pcondrsubr( U, X, ipar )
    call matvecsubr( U, R, ipar )
    U = B - R
    call pcondlsubr( R, U, ipar )
    RTLD = R
    U = 0
    oldrho = 1; omega2 = 1; alpha = 0

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    oldrho = -omega2 * oldrho

    !
    ! This is the even BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, R, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    U = R - beta * U

    call pcondrsubr( V, U, ipar )
    call matvecsubr( V, T1V, ipar )
    call pcondlsubr( V, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
    R = R - alpha * V

    call pcondrsubr( S, R, ipar )
    call matvecsubr( S, T1V, ipar )
    call pcondlsubr( S, T1V, ipar )

    X = X + alpha * U

    !
    ! This is the odd BiCG step
    !

    rho = dotprodfun( HUTI_NDIM, RTLD, 1, S, 1 )
    if ( rho .eq. 0 ) then
       HUTI_INFO = HUTI_BICGSTAB_2_RHO
       go to 1000
    end if

    beta = ( rho * alpha ) / oldrho
    oldrho = rho
    V = S - beta * V

    call pcondrsubr( W, V, ipar )
    call matvecsubr( W, T1V, ipar )
    call pcondlsubr( W, T1V, ipar )

    alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, W, 1 )
    U = R - beta * U
    R = R - alpha * V
    S = S - alpha * W

    call pcondrsubr( T, S, ipar )
    call matvecsubr( T, T1V, ipar )
    call pcondlsubr( T, T1V, ipar )

    !
    ! This is the GCR(2) part
    !

    omega1 = dotprodfun( HUTI_NDIM, R, 1, S, 1 )
    myy = dotprodfun( HUTI_NDIM, S, 1, S, 1 )
    delta = dotprodfun( HUTI_NDIM, S, 1, T, 1 )
    tau = dotprodfun( HUTI_NDIM, T, 1, T, 1 )
    omega2 = dotprodfun( HUTI_NDIM, R, 1, T, 1 )

    tau = tau - ( delta * delta ) / myy
    omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau
    omega1 = ( omega1 - delta * omega2 ) / myy

    X = X + omega1 * R + omega2 * S + alpha * U
    R = R - omega1 * S - omega2 * T

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
       T1V = omega1 * R + omega2 * S + alpha * U
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call pcondrsubr( S, X, ipar )
       call matvecsubr( S, T1V, ipar )
       T1V = T1V - B
       call pcondlsubr( S, T1V, ipar )
       residual = normfun( HUTI_NDIM, S, 1 )
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
    U = U - omega1 * V - omega2 * W

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

  end subroutine  huti_zbicgstab_2solv

  !*************************************************************************

end module huti_bicgstab_2
