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
! Subroutines to implement Transpose Free QMR iteration
!
! $Id: huti_tfqmr.src,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $

#include "huti_intdefs.h"
MAKE_INCLUDE(#,include, "huti_fdefs.h" )

!*************************************************************************
!*************************************************************************
!
! These subroutines are based on a paper by Roland W. Freund:
! "A Transpose-Free Quasi-Minimal Residual Algorithm for Non-Hermitian
!  Linear Systems", 1993 (SIAM J. Sci. Comput, March 1993)
!
! All matrix-vector operations are done externally, so we do not need
! to know about the matrix structure (sparse or dense). Memory allocation
! for the working arrays has also been done externally.

!*************************************************************************
! Work array is used in the following order:
! work(:,1) = v
! work(:,2) = y
! work(:,3) = y new
! work(:,4) = r tilde (zero)
! work(:,5) = t1v (temporary for matrix-vector operations)
! work(:,6) = t2v (temporary for matrix-vector operations)
! work(:,7) = w
! work(:,8) = d
! work(:,9) = r
! work(:,10) = trv (temporary vector for residual computations)
!
!*************************************************************************
! Definitions to make the code more understandable and to make it look
! like the pseudo code
!
MAKE_DEFINE(#,define, X, xvec )
MAKE_DEFINE(#,define, B, rhsvec )

MAKE_DEFINE(#,define, V, work(:,1) )
MAKE_DEFINE(#,define, V_ind, 1 )
MAKE_DEFINE(#,define, Y, work(:,2) )
MAKE_DEFINE(#,define, Y_ind, 2 )
MAKE_DEFINE(#,define, YNEW, work(:,3) )
MAKE_DEFINE(#,define, YNEW_ind, 3 )
MAKE_DEFINE(#,define, RTLD, work(:,4) )
MAKE_DEFINE(#,define, RTLD_ind, 4 )
MAKE_DEFINE(#,define, T1V, work(:,5) )
MAKE_DEFINE(#,define, T1V_ind, 5 )
MAKE_DEFINE(#,define, T2V, work(:,6) )
MAKE_DEFINE(#,define, T2V_ind, 6 )
MAKE_DEFINE(#,define, W, work(:,7) )
MAKE_DEFINE(#,define, W_ind, 7 )
MAKE_DEFINE(#,define, D, work(:,8) )
MAKE_DEFINE(#,define, D_ind, 8 )
MAKE_DEFINE(#,define, R, work(:,9) )
MAKE_DEFINE(#,define, R_ind, 9 )
MAKE_DEFINE(#,define, TRV, work(:,10) )
MAKE_DEFINE(#,define, TRV_ind, 10 )

! This is the magic ratio for upperb and tolerance used in upper bound
! convergence test
MAKE_DEFINE(#,define, UPPERB_TOL_RATIO, 10.0 )

!*************************************************************************
!*************************************************************************
! PRECISION_COMMENT version
!*************************************************************************
!*************************************************************************

subroutine MAKE_SUBRN( huti_,tfqmrsolv ) ( ndim, wrkdim, xvec, rhsvec, ipar,&
                            dpar, work, matvecsubr, pcondlsubr, pcondrsubr, &
                            dotprodfun, normfun, stopcfun )


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

  F_PRECISION_TYPE :: rho, oldrho, eta, tau, gamma, oldgamma, alpha
  F_PRECISION_TYPE :: beta, c
  integer :: iter_count

  NORMFUN_PREC_TYPE :: residual, upperb, rhsnorm, precrhsnorm

  !
  ! End of variable declarations
  !*********************************************************************

  !*********************************************************************
  ! The actual TFQMR begins here (look the pseudo code in the
  ! "A Transpose-Free..."-paper, algorithm 5.1)
  !
  ! First the initialization part
  !

  iter_count = 1

  ! The following applies for all matrix operations in this solver

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED

  ! Norms of right-hand side vector are used in convergence tests

  if ( HUTI_STOPC .eq. HUTI_TRESID_SCALED_BYB .or. & 
       HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYB .or. &
       HUTI_STOPC .eq. HUTI_UPPERB_STOPC ) then
     rhsnorm = normfun( HUTI_NDIM, B, 1 )
  end if
  if ( HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYPRECB ) then
     call pcondlsubr( D, B, ipar )
     precrhsnorm = normfun( HUTI_NDIM, D, 1 )
  end if

  !
  ! Part 1A - 1C
  !

  ! Generate vector X if needed

  if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
     call MAKE_SUBRN( huti_,randvec )  ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  call pcondrsubr( D, X, ipar )
  call matvecsubr( D, R, ipar )
  D = B - R
  call pcondlsubr( R, D, ipar )

  Y = R; W = R
  call pcondrsubr( V, Y, ipar )
  call matvecsubr( V, D, ipar )
  call pcondlsubr( V, D, ipar )
  T2V = V

  D = 0
  tau = normfun( HUTI_NDIM, R, 1 )
  oldgamma = 0; gamma = 0; eta = 0

  RTLD = R
  oldrho = dotprodfun ( HUTI_NDIM, RTLD, 1, R, 1 )
  if ( oldrho .eq. 0 ) then
     HUTI_INFO = HUTI_TFQMR_RHO
     go to 1000
  end if

  !
  ! This is where the loop starts (that is we continue from here after
  ! the first iteration)
  !
  !
  ! Part 2A
  !

300 continue

  alpha = oldrho / dotprodfun( HUTI_NDIM, RTLD, 1, V, 1 )
  YNEW = Y - alpha * V

  !
  ! Part 2B
  !
  !
  ! This is the inner loop from 2n-1 to 2n
  !

  ! First the 2n-1 case

  ! Note: We have already MATRIX * Y in T2V

  W = W - alpha * T2V
  gamma = ( normfun( HUTI_NDIM, W, 1 )) / tau
  c = 1 / sqrt( 1 + gamma * gamma )
  tau = tau * gamma * c

  D = Y + ((oldgamma * oldgamma * eta) / alpha) * D
  eta = c * c * alpha
  X = X + eta * D

  oldgamma = gamma

  !
  ! Check the convergence against selected stopping criterion
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, ipar )
     residual = normfun( HUTI_NDIM, TRV, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, ipar )
     residual = normfun( HUTI_NDIM, TRV, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, ipar )
     residual = normfun( HUTI_NDIM, TRV, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     R = eta * D
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_UPPERB_STOPC)
     upperb = real( sqrt( 2.0 * iter_count ) * tau / rhsnorm)
     if ( ( upperb / HUTI_TOLERANCE ) .lt. UPPERB_TOL_RATIO ) then
	call pcondrsubr( TRV, X, ipar )
        call matvecsubr( TRV, R, ipar )
        TRV = R - B
	call pcondlsubr( R, TRV, ipar )
        residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
     else
        residual = upperb
     end if
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  end select

  if ( residual .lt. HUTI_TOLERANCE ) then
     HUTI_INFO = HUTI_CONVERGENCE
     go to 1000
  end if

  !
  ! And then the 2n case
  !

  call pcondrsubr( T1V, YNEW, ipar )
  call matvecsubr( T1V, R, ipar )
  call pcondlsubr( T1V, R, ipar )

  W = W - alpha * T1V
  gamma = ( normfun( HUTI_NDIM, W, 1 )) / tau
  c = 1 / sqrt( 1 + gamma * gamma )
  tau = tau * gamma * c

  D = YNEW + ((oldgamma * oldgamma * eta) / alpha) * D
  eta = c * c * alpha
  X = X + eta * D

  oldgamma = gamma

  !
  ! Check the convergence against selected stopping criterion
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, ipar )
     residual = normfun( HUTI_NDIM, TRV, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, ipar )
     residual = normfun( HUTI_NDIM, TRV, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( TRV, R, 1 )
     residual = normfun( HUTI_NDIM, TRV, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     R = eta * D
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_UPPERB_STOPC)
     upperb = real( sqrt( 2.0 * iter_count ) * tau / rhsnorm)
     if ( ( upperb / HUTI_TOLERANCE ) .lt. UPPERB_TOL_RATIO ) then
	call pcondrsubr( TRV, X, ipar )
        call matvecsubr( TRV, R, ipar )
        TRV = R - B
	call pcondlsubr( R, TRV, ipar )
        residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
     else
        residual = upperb
     end if
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call pcondrsubr( TRV, X, ipar )
     call matvecsubr( TRV, R, ipar )
     TRV = R - B
     call pcondlsubr( R, TRV, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  end select

  if ( residual .lt. HUTI_TOLERANCE ) then
     HUTI_INFO = HUTI_CONVERGENCE
     go to 1000
  end if

  !
  ! Produce debugging output if desired
  !

  if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
     if ( mod(iter_count, HUTI_DBUGLVL) .eq. 0 ) then
        if ( HUTI_STOPC .eq. HUTI_UPPERB_STOPC ) then
           write (*, '(I8, 2E17.7)') iter_count, residual, upperb
        else
           write (*, '(I8, E17.7)') iter_count, residual
        end if
     end if
  end if

  !
  ! Part 2C
  !

  rho = dotprodfun( HUTI_NDIM, RTLD, 1, W, 1 )
  beta = rho / oldrho
  YNEW = W + beta * YNEW
  call pcondrsubr( T2V, YNEW, ipar )
  call matvecsubr( T2V, R, ipar )
  call pcondlsubr( T2V, R, ipar )

  ! Note: we still have MATRIX * YNEW in T1V

  V = T2V + beta * T1V + beta * beta * V
  Y = YNEW

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

  ! Compute the unpreconditioned X

  call pcondrsubr( TRV, X, ipar )
  X = TRV

  if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
     if ( HUTI_STOPC .eq. HUTI_UPPERB_STOPC ) then
	write (*, '(I8, 2E17.7)') iter_count, residual, upperb
     else
	write (*, '(I8, E17.7)') iter_count, residual
     end if
  end if
  HUTI_ITERS = iter_count

  return

  ! End of execution
  !*********************************************************************

end subroutine MAKE_SUBRN( huti_,tfqmrsolv )

!*************************************************************************
