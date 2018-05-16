module huti_bicgstab
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
  ! Subroutine to implement BiConjugate Gradient Stabilised iteration
  !
  ! $Id: huti_bicgstab.src,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $


#include  "huti_fdefs.h" 

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


contains


  !*************************************************************************
  !*************************************************************************
  ! Single precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_sbicgstabsolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    real :: rho, oldrho, alpha, beta, omega
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
       !HUTI_INFO = HUTI_BICGSTAB_SNORM
       HUTI_INFO = HUTI_CONVERGENCE
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
    
    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
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

  end subroutine  huti_sbicgstabsolv

  !*************************************************************************



  !*************************************************************************
  !*************************************************************************
  ! Double precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_dbicgstabsolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    double precision :: rho, oldrho, alpha, beta, omega
    integer :: iter_count, i

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
       call  huti_drandvec   ( X, ipar )
    else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
#ifdef _OPENMP
       !$OMP PARALLEL DO
       do i=1,ndim
          xvec(i) = 1
       end do
       !$OMP END PARALLEL DO
#else
       X = 1
#endif
    end if

    call matvecsubr( X, R, ipar )
#ifdef _OPENMP
    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i=1,ndim
       work(i,R_ind) = rhsvec(i) - work(i,R_ind)
       work(i,RTLD_ind) = work(i,R_ind)
    end do
    !$OMP END DO NOWAIT
    !$OMP DO
    do i=1,ndim
       work(i,P_ind) = 0
       work(i,V_ind) = 0
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
#else
    R = B - R
    RTLD = R
    P = 0; V = 0
#endif
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

#ifdef _OPENMP
    !$OMP PARALLEL DO
    do i=1,ndim
       work(i,P_ind) = work(i,R_ind) + beta * ( work(i,P_ind) - omega * work(i,V_ind) )
    end do
    !$OMP END PARALLEL DO
#else
    P = R + beta * ( P - omega * V )
#endif

    call pcondlsubr( V, P, ipar )
    call pcondrsubr( T1V, V, ipar )
    call matvecsubr( T1V, V, ipar )

    alpha = rho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
#ifdef _OPENMP
    !$OMP PARALLEL DO
    do i=1,ndim
       work(i,S_ind) = work(i,R_ind) - alpha * work(i,V_ind)
    end do
    !$OMP END PARALLEL DO
#else
    S = R - alpha * V
#endif

    residual = normfun( HUTI_NDIM, S, 1 )
    if ( residual .lt. HUTI_EPSILON ) then
#ifdef _OPENMP
    !$OMP PARALLEL DO
    do i=1,ndim
       xvec(i) = xvec(i) + alpha * work(i,T1V_ind)
    end do
    !$OMP END PARALLEL DO
#else
       X = X + alpha * T1V
#endif
       HUTI_INFO = HUTI_CONVERGENCE
       !HUTI_INFO = HUTI_BICGSTAB_SNORM
       go to 1000
    end if

    call pcondlsubr( T, S, ipar )
    call pcondrsubr( T2V, T, ipar )
    call matvecsubr( T2V, T, ipar )

    omega = ( dotprodfun( HUTI_NDIM, T, 1, S, 1 ) ) / &
         ( dotprodfun( HUTI_NDIM, T, 1, T, 1 ) )

#ifdef _OPENMP
    !$OMP PARALLEL DO
    do i=1,ndim
       xvec(i) = xvec(i) + alpha * work(i,T1V_ind) + omega * work(i,T2V_ind)
       work(i,R_ind) = work(i,S_ind) - omega * work(i,T_ind)
    end do
    !$OMP END PARALLEL DO
#else
    X = X + alpha * T1V + omega * T2V
    R = S - omega * T
#endif

    !
    ! Check the convergence against selected stopping criterion
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call matvecsubr( X, T2V, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       do i=1,ndim
          work(i,T1V_ind) = work(i,T2V_ind) - rhsvec(i)
       end do
       !$OMP END PARALLEL DO
#else
       T1V = T2V - B
#endif
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call matvecsubr( X, T2V, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       do i=1,ndim
          work(i,T1V_ind) = work(i,T2V_ind) - rhsvec(i)
       end do
       !$OMP END PARALLEL DO
#else
       T1V = T2V - B
#endif
       residual = normfun( HUTI_NDIM, T1V, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       residual = normfun( HUTI_NDIM, R, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
#ifdef _OPENMP
       !$OMP PARALLEL DO
       do i=1,ndim
          work(i,T1V_ind) = alpha * work(i,T1V_ind) + omega * work(i,T2V_ind)
       end do
       !$OMP END PARALLEL DO
#else
       T1V = alpha * T1V + omega * T2V
#endif
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call matvecsubr( X, T2V, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       do i=1,ndim
          work(i,T1V_ind) = work(i,T2V_ind) - rhsvec(i)
       end do
       !$OMP END PARALLEL DO
#else
       T1V = T2V - B
#endif
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
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

  end subroutine  huti_dbicgstabsolv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_cbicgstabsolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    complex :: rho, oldrho, alpha, beta, omega
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
       call pcondlsubr( P, B, ipar )
       precrhsnorm = normfun( HUTI_NDIM, P, 1 )
    end if

    ! Generate vector X if needed

    if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
       call  huti_crandvec   ( X, ipar )
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
       HUTI_INFO = HUTI_CONVERGENCE
       !HUTI_INFO = HUTI_BICGSTAB_SNORM
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
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

  end subroutine  huti_cbicgstabsolv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Double complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_zbicgstabsolv  ( ndim, wrkdim, xvec, rhsvec, &
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
       !HUTI_INFO = HUTI_BICGSTAB_SNORM
       HUTI_INFO = HUTI_CONVERGENCE
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

    IF( residual /= residual .OR. residual > HUTI_MAXTOLERANCE ) THEN
      HUTI_INFO = HUTI_DIVERGENCE
      GOTO 1000
    END IF
    
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


end module huti_bicgstab
