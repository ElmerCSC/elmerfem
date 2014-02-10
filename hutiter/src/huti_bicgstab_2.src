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

!
! Subroutine to implement BiConjugate Gradient Stabilised (2) iteration
!
! $Id: huti_bicgstab_2.src,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $

#include "huti_intdefs.h"
MAKE_INCLUDE(#,include, "huti_fdefs.h" )

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

MAKE_DEFINE(#,define, X, xvec )
MAKE_DEFINE(#,define, B, rhsvec )

MAKE_DEFINE(#,define, RTLD, work(:,1) )
MAKE_DEFINE(#,define, RTLD_ind, 1 )
MAKE_DEFINE(#,define, U, work(:,2) )
MAKE_DEFINE(#,define, U_ind, 2 )
MAKE_DEFINE(#,define, T1V, work(:,3) )
MAKE_DEFINE(#,define, T1V_ind, 3 )
MAKE_DEFINE(#,define, V, work(:,4) )
MAKE_DEFINE(#,define, V_ind, 4 )
MAKE_DEFINE(#,define, S, work(:,5) )
MAKE_DEFINE(#,define, S_ind, 5 )
MAKE_DEFINE(#,define, W, work(:,6) )
MAKE_DEFINE(#,define, W_ind, 6 )
MAKE_DEFINE(#,define, T, work(:,7) )
MAKE_DEFINE(#,define, T_ind, 7 )
MAKE_DEFINE(#,define, R, work(:,8) )
MAKE_DEFINE(#,define, R_ind, 8 )
  
!*************************************************************************
!*************************************************************************
! PRECISION_COMMENT version
!*************************************************************************
!*************************************************************************

subroutine MAKE_SUBRN( huti_,bicgstab_2solv ) ( ndim, wrkdim, xvec, rhsvec, &
                            ipar, dpar, work, matvecsubr, pcondlsubr, &
                            pcondrsubr, dotprodfun, normfun, stopcfun )


  implicit none

  external matvecsubr, pcondlsubr, pcondrsubr
  external dotprodfun, normfun, stopcfun
  F_PRECISION_TYPE :: dotprodfun
  NORMFUN_PREC_TYPE :: normfun
  NORMFUN_PREC_TYPE :: stopcfun

  ! Parameters

  integer :: ndim, wrkdim
  F_PRECISION_TYPE, dimension(ndim) :: xvec, rhsvec
  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
  F_PRECISION_TYPE, dimension(ndim,wrkdim) :: work

  ! Local variables

  F_PRECISION_TYPE :: rho, oldrho, alpha, beta, omega1, omega2
  F_PRECISION_TYPE :: tau, delta, myy
  integer :: iter_count

  NORMFUN_PREC_TYPE :: residual, rhsnorm, precrhsnorm

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
     call MAKE_SUBRN( huti_,randvec )  ( X, ipar )
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

end subroutine MAKE_SUBRN( huti_,bicgstab_2solv )

!*************************************************************************
