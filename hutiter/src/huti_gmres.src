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
! Subroutines to implement Generalized Minimum Residual iterative method
!
! $Id: huti_gmres.src,v 1.1.1.1 2005/04/15 10:31:18 vierinen Exp $

#include "huti_intdefs.h"
MAKE_INCLUDE(#,include, "huti_fdefs.h" )

!*************************************************************************
!*************************************************************************
!
! These subroutines are based on a book by Barret et al.:
! "Templates for the Solution of Linear Systems: Building Blocks for
!  Iterative Methods", 1993.
!
! All matrix-vector operations are done externally, so we do not need
! to know about the matrix structure (sparse or dense). So has the
! memory allocation for the working arrays done also externally.

!*************************************************************************
! Work array is used in the following order:
! work(:,1) = z
! work(:,2) = p
! work(:,3) = q
! work(:,4) = r
!
! Definitions to make the code more understandable and to make it look
! like the pseudo code (these are commond to all precisions)
!

MAKE_DEFINE(#,define, X, xvec )
MAKE_DEFINE(#,define, B, rhsvec )

MAKE_DEFINE(#,define, W, work(:,1) )
MAKE_DEFINE(#,define, W_ind, 1 )
MAKE_DEFINE(#,define, R, work(:,2) )
MAKE_DEFINE(#,define, R_ind, 2 )
MAKE_DEFINE(#,define, S, work(:,3) )
MAKE_DEFINE(#,define, S_ind, 3 )
MAKE_DEFINE(#,define, VTMP, work(:,4) )
MAKE_DEFINE(#,define, VTMP_ind, 4 )
MAKE_DEFINE(#,define, T1V, work(:,5) )
MAKE_DEFINE(#,define, T1V_ind, 5 )
MAKE_DEFINE(#,define, V, work(:,6) )
MAKE_DEFINE(#,define, V_ind, 6 )

#define Varr2(ind1,ind2) work(:,V_ind+ind1-1:V_ind+ind2-1)
#define Varr(ind) work(:,V_ind+ind-1)


!*************************************************************************
  
!*************************************************************************
!*************************************************************************
! PRECISION_COMMENT version
!*************************************************************************
!*************************************************************************

subroutine MAKE_SUBRN( huti_,gmressolv ) ( ndim, wrkdim, xvec, rhsvec, &
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

  integer :: iter_count
  NORMFUN_PREC_TYPE :: residual, rhsnorm, precrhsnorm
  NORMFUN_PREC_TYPE :: bnrm, alpha, beta

  INTEGER :: reit, info, i, j, k, l, m, n

  F_PRECISION_TYPE :: temp, temp2, error, gamma

  ! Local arrays

  F_PRECISION_TYPE, DIMENSION(HUTI_GMRES_RESTART+1,HUTI_GMRES_RESTART+1) :: H
  F_PRECISION_TYPE, &
       DIMENSION((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1)) :: HLU
  F_PRECISION_TYPE, DIMENSION(HUTI_GMRES_RESTART+1) :: CS, SN, Y

  !
  ! End of variable declarations
  !*********************************************************************

  !*********************************************************************
  ! The actual GMRES begins here (look the pseudo code in the
  ! "Templates..."-book, page 20)
  !
  ! First the initialization part
  !

  iter_count = 1

  bnrm = normfun( HUTI_NDIM, B, 1 )

  ! Norms of right-hand side vector are used in convergence tests

  if ( HUTI_STOPC .eq. HUTI_TRESID_SCALED_BYB .or. & 
       HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYB ) then
     rhsnorm = bnrm
  end if
  if ( HUTI_STOPC .eq. HUTI_PRESID_SCALED_BYPRECB ) then
     call pcondlsubr( T1V, B, ipar )
     precrhsnorm = normfun( HUTI_NDIM, T1V, 1 )
  end if

  ! The following applies for all matrix operations in this solver

  HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED

  ! Generate vector X if needed

  if ( HUTI_INITIALX .eq. HUTI_RANDOMX ) then
     call MAKE_SUBRN( huti_,randvec )  ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  call pcondrsubr( T1V, X, ipar )
  call matvecsubr( T1V, R, ipar )
  T1V = B - R
  call pcondlsubr( R, T1V, ipar )

  m = HUTI_GMRES_RESTART
  Varr2(1,m+1) = 0
  H = 0
  CS = 0
  SN = 0
  VTMP = 0
  work(1,VTMP_ind) = 1.0

  !
  ! This is where the loop starts (that is we continue from here after
  ! the first iteration)
  !

300 continue

  call pcondrsubr( T1V, X, ipar )
  call matvecsubr( T1V, R, ipar )
  T1V = B - R
  call pcondlsubr( R, T1V, ipar )

  alpha = normfun( HUTI_NDIM, R, 1 )
  if ( alpha .eq. 0 ) then
     HUTI_INFO = HUTI_GMRES_ALPHA
     go to 1000
  end if

  Varr(1) = R / alpha
  S = alpha * VTMP

  !
  ! Construct orthonormal
  !

  DO i = 1, m
     call pcondrsubr( W, Varr(i), ipar )
     call matvecsubr( W, T1V, ipar )
     call pcondlsubr( W, T1V, ipar )

     DO k = 1, i
	H(k,i) = dotprodfun( HUTI_NDIM, W, 1, Varr(k), 1 )
	W = W - H(k,i) * Varr(k)
     END DO

     beta = normfun( HUTI_NDIM, W, 1 )
     if ( beta .eq. 0 ) then
	HUTI_INFO = HUTI_GMRES_BETA
	go to 1000
     end if

     H(i+1,i) = beta
     Varr(i+1) = W / H(i+1, i)

     !
     ! Compute the Givens rotation
     !

     DO k = 1, i-1
	temp = CS(k) * H(k,i) + SN(k) * H(k+1,i)
	H(k+1,i) = -1 * SN(k) * H(k,i) + CS(k) * H(k+1,i)
	H(k,i) = temp
     END DO

     IF ( H(i+1,i) .eq. 0 ) THEN
	CS(i) = 1; SN(i) = 0
     ELSE
	IF ( abs( H(i+1,i) ) .gt. abs( H(i,i) ) ) THEN
	   temp2 = H(i,i) / H(i+1,i)
	   SN(i) = 1 / sqrt( 1 + ( temp2 * temp2 ))
	   CS(i) = temp2 * SN(i)
	ELSE
	   temp2 = H(i+1,i) / H(i,i)
	   CS(i) = 1 / sqrt( 1 + ( temp2 * temp2 ))
	   SN(i) = temp2 * CS(i)
	END IF
     END IF

     temp = CS(i) * work(i,S_ind)
     work(i+1,S_ind) = -1 * SN(i) * work(i,S_ind)
     work(i,S_ind) = temp

     H(i,i) = ( CS(i) * H(i,i) ) + ( SN(i) * H(i+1,i) )
     H(i+1,i) = 0

     error = abs( work(i+1,S_ind) ) / bnrm 
     IF ( REAL( error ) .lt.  HUTI_TOLERANCE ) THEN

	HLU = 0; j = 1
	do k = 1, i
	   do l = 1, i
	      HLU(j) = H(l,k)
	      j = j + 1
	   end do
	end do

	call MAKE_SUBRN( huti_,lusolve ) ( i, HLU, Y, work(:,S_ind) )

	X = X + MATMUL( Varr2(1,i), Y(1:i) )

	EXIT
     END IF

  END DO

  IF ( REAL( error ) .lt. HUTI_TOLERANCE ) THEN
     GOTO 500
  END IF

  HLU = 0; j = 1
  do k = 1, m
     do l = 1, m
	HLU(j) = H(l,k)
	j = j + 1
     end do
  end do
  
  call MAKE_SUBRN( huti_,lusolve ) ( m, HLU, Y, work(:,S_ind) )

  X = X + MATMUL( Varr2(1,m), Y(1:m) )

500 CONTINUE

  !
  ! Check the convergence
  !

  select case (HUTI_STOPC)
  case (HUTI_TRUERESIDUAL)
     call pcondrsubr( T1V, X, ipar )
     call matvecsubr( T1V, R, ipar )
     T1V = B - R
     call pcondlsubr( R, T1V, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  case (HUTI_TRESID_SCALED_BYB)
     call pcondrsubr( T1V, X, ipar )
     call matvecsubr( T1V, R, ipar )
     T1V = B - R
     call pcondlsubr( R, T1V, ipar )
     residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
  case (HUTI_PSEUDORESIDUAL)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( T1V, R, ipar )
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_PRESID_SCALED_BYB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( T1V, R, ipar )
     residual = normfun( HUTI_NDIM, T1V, 1 ) / rhsnorm
  case (HUTI_PRESID_SCALED_BYPRECB)
     call matvecsubr( X, R, ipar )
     R = R - B
     call pcondlsubr( T1V, R, ipar )
     residual = normfun( HUTI_NDIM, T1V, 1 ) / precrhsnorm
  case (HUTI_XDIFF_NORM)
     T1V = MATMUL( Varr2(1,m), Y(1:m) )
     residual = normfun( HUTI_NDIM, T1V, 1 )
  case (HUTI_USUPPLIED_STOPC)
     residual = stopcfun( X, B, R, ipar, dpar )
  case default
     call pcondrsubr( T1V, X, ipar )
     call matvecsubr( T1V, R, ipar )
     T1V = R - B
     call pcondlsubr( R, T1V, ipar )
     residual = normfun( HUTI_NDIM, R, 1 )
  end select
     
  work(m+1,S_ind) = normfun( HUTI_NDIM, R, 1 )

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

end subroutine MAKE_SUBRN( huti_,gmressolv )

!*************************************************************************
