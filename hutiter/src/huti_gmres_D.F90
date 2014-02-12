
!
! Subroutines to implement Generalized Minimum Residual iterative method
!
! $Id: huti_gmres_D.F90,v 1.5 2005/05/04 09:57:39 vierinen Exp $



































#include  "huti_fdefs.h" 

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

#define  X  xvec 
#define  B  rhsvec 

#define  W  work(:,1) 
#define  W_ind  1 
#define  R  work(:,2) 
#define  R_ind  2 
#define  S  work(:,3) 
#define  S_ind  3 
#define  VTMP  work(:,4) 
#define  VTMP_ind  4 
#define  T1V  work(:,5) 
#define  T1V_ind  5 
#define  V  work(:,6) 
#define  V_ind  6 





!*************************************************************************
  
!*************************************************************************
!*************************************************************************
! Double precision version
!*************************************************************************
!*************************************************************************

subroutine  huti_dgmressolv  ( ndim, wrkdim, xvec, rhsvec, &
                          ipar, dpar, work, matvecsubr, pcondlsubr, &
                          pcondrsubr, dotprodfun, normfun, stopcfun )


  implicit none

  external matvecsubr, pcondlsubr, pcondrsubr
  external dotprodfun, normfun, stopcfun
  double precision :: dotprodfun
  double precision :: normfun
  double precision :: stopcfun

  ! Parameters

  integer :: ndim, wrkdim
  double precision, dimension(ndim) :: xvec, rhsvec
  integer, dimension(HUTI_IPAR_DFLTSIZE) :: ipar
  double precision, dimension(HUTI_DPAR_DFLTSIZE) :: dpar
  double precision, dimension(ndim,wrkdim) :: work

  ! Local variables

  integer :: iter_count
  double precision :: residual, rhsnorm, precrhsnorm
  double precision :: bnrm, alpha, beta

  INTEGER :: reit, info, i, j, k, l, m, n

  double precision :: temp, temp2, error, gamma

  ! Local arrays

  double precision, DIMENSION(HUTI_GMRES_RESTART+1,HUTI_GMRES_RESTART+1) :: H
  double precision, &
       DIMENSION((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1)) :: HLU
  double precision, DIMENSION(HUTI_GMRES_RESTART+1) :: CS, SN, Y

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
     call  huti_drandvec   ( X, ipar )
  else if ( HUTI_INITIALX .ne. HUTI_USERSUPPLIEDX ) then
     X = 1
  end if

  call pcondrsubr( T1V, X, ipar )
  call matvecsubr( T1V, R, ipar )
  T1V = B - R
  call pcondlsubr( R, T1V, ipar )

  m = HUTI_GMRES_RESTART
  work(:,V_ind+1-1:V_ind+m+1-1) = 0
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

  work(:,V_ind+1-1) = R / alpha
  S = alpha * VTMP

  !
  ! Construct orthonormal
  !

  DO i = 1, m
     call pcondrsubr( W, work(:,V_ind+i-1), ipar )
     call matvecsubr( W, T1V, ipar )
     call pcondlsubr( W, T1V, ipar )

     DO k = 1, i
	H(k,i) = dotprodfun( HUTI_NDIM, W, 1, work(:,V_ind+k-1), 1 )
	W = W - H(k,i) * work(:,V_ind+k-1)
     END DO

     beta = normfun( HUTI_NDIM, W, 1 )
     if ( beta .eq. 0 ) then
	HUTI_INFO = HUTI_GMRES_BETA
	go to 1000
     end if

     H(i+1,i) = beta
     work(:,V_ind+i+1-1) = W / H(i+1, i)

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

	call  huti_dlusolve  ( i, HLU, Y, work(:,S_ind) )

	X = X + MATMUL( work(:,V_ind+1-1:V_ind+i-1), Y(1:i) )

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
  
  call  huti_dlusolve  ( m, HLU, Y, work(:,S_ind) )

  X = X + MATMUL( work(:,V_ind+1-1:V_ind+m-1), Y(1:m) )

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
     T1V = MATMUL( work(:,V_ind+1-1:V_ind+m-1), Y(1:m) )
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

end subroutine  huti_dgmressolv 

!*************************************************************************
