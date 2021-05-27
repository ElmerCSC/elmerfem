module huti_gmres
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
  ! Subroutines to implement Generalized Minimum Residual iterative method
  !

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

contains

  !*************************************************************************

  !*************************************************************************
  !*************************************************************************
  ! Single precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_sgmressolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    integer :: iter_count
    real :: residual, rhsnorm, precrhsnorm
    real :: bnrm, alpha, beta

    INTEGER :: reit, info, i, j, k, l, m, n

    real :: temp, temp2, error, gamma

    ! Local arrays

    real, DIMENSION(HUTI_GMRES_RESTART+1,HUTI_GMRES_RESTART+1) :: H
    real, &
         DIMENSION((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1)) :: HLU
    real, DIMENSION(HUTI_GMRES_RESTART+1) :: CS, SN, Y

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
       call  huti_srandvec   ( X, ipar )
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

          call  huti_slusolve  ( i, HLU, Y, work(:,S_ind) )

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

    call  huti_slusolve  ( m, HLU, Y, work(:,S_ind) )

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
          WRITE (*, '(A, I8, E11.4)') '   gmres:',iter_count, residual
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
       WRITE (*, '(A, I8, E11.4)') '   gmres:',iter_count, residual
    end if

    HUTI_ITERS = iter_count
    return

    ! End of execution
    !*********************************************************************

  end subroutine  huti_sgmressolv

  !*************************************************************************

  !*************************************************************************
  !*************************************************************************
  ! Double precision version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_dgmressolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    integer :: iter_count
    double precision :: residual, rhsnorm, precrhsnorm
    double precision :: bnrm, alpha, beta

    INTEGER :: reit, info, i, j, k, l, m, n, ii, jj

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
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         xvec(ii) = 1
       END DO
       !$OMP END PARALLEL DO
#else   
        X = 1
#endif
    end if

    call pcondrsubr( T1V, X, ipar )
    call matvecsubr( T1V, R, ipar )
#ifdef _OPENMP
    !$OMP PARALLEL DO
    DO ii=1,ndim
      work(ii,T1V_ind) = rhsvec(ii) - work(ii,R_ind)
    END DO
    !$OMP END PARALLEL DO
#else
    T1V = B - R
#endif
    call pcondlsubr( R, T1V, ipar )

    m = HUTI_GMRES_RESTART
#ifdef _OPENMP
    !$OMP PARALLEL PRIVATE(j)
    DO jj=V_ind+1-1,V_ind+m+1-1
      !$OMP DO
      DO ii=1,ndim
        work(ii,jj) = 0
      END DO
      !$OMP END DO
    END DO
    !$OMP DO
    DO ii=1,ndim
      work(ii,VTMP_ind) = 0
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
#else
    work(:,V_ind+1-1:V_ind+m+1-1) = 0
    VTMP = 0
#endif
    H = 0
    CS = 0
    SN = 0
    work(1,VTMP_ind) = 1.0

    !
    ! This is where the loop starts (that is we continue from here after
    ! the first iteration)
    !

300 continue

    call pcondrsubr( T1V, X, ipar )
    call matvecsubr( T1V, R, ipar )
#ifdef _OPENMP
    !$OMP PARALLEL DO
    DO ii=1,ndim
      work(ii,T1V_ind) = rhsvec(ii) - work(ii,R_ind)
    END DO
    !$OMP END PARALLEL DO
#else
    T1V = B - R
#endif
    call pcondlsubr( R, T1V, ipar )

    alpha = normfun( HUTI_NDIM, R, 1 )
    if ( alpha .eq. 0 ) then
       HUTI_INFO = HUTI_GMRES_ALPHA
       go to 1000
    end if

#ifdef _OPENMP
    !$OMP PARALLEL
    !$OMP DO
    DO ii=1,ndim
      work(ii,V_ind+1-1) = work(ii,R_ind) / alpha
    END DO
    !$OMP END DO
    !$OMP DO
    DO ii=1,ndim
      work(ii,S_ind) = alpha * work(ii,VTMP_ind)
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
#else
    work(:,V_ind+1-1) = R / alpha
    S = alpha * VTMP
#endif

    !
    ! Construct orthonormal
    !

    DO i = 1, m
       call pcondrsubr( W, work(:,V_ind+i-1), ipar )
       call matvecsubr( W, T1V, ipar )
       call pcondlsubr( W, T1V, ipar )

       DO k = 1, i
          H(k,i) = dotprodfun( HUTI_NDIM, W, 1, work(:,V_ind+k-1), 1 )
#ifdef _OPENMP
          !$OMP PARALLEL DO
          DO ii=1,ndim
            work(ii,W_ind) = work(ii,W_ind) - H(k,i) * work(ii,V_ind+k-1)
          END DO
          !$OMP END PARALLEL DO
#else
          W = W - H(k,i) * work(:,V_ind+k-1)
#endif
       END DO

       beta = normfun( HUTI_NDIM, W, 1 )
       if ( beta .eq. 0 ) then
          HUTI_INFO = HUTI_GMRES_BETA
          go to 1000
       end if

       H(i+1,i) = beta
#if _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,V_ind+i+1-1) = work(ii,W_ind) / H(i+1, i)
       END DO
       !$OMP END PARALLEL DO
#else
       work(:,V_ind+i+1-1) = W / H(i+1, i)
#endif
       
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

#ifdef _OPENMP
          CALL DGEMV('N',ndim,i,1D0,work(1,V_ind),ndim,Y,1,1D0,xvec,1)
#else
          X = X + MATMUL( work(:,V_ind+1-1:V_ind+i-1), Y(1:i) )
#endif
          
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

#ifdef _OPENMP
    CALL DGEMV('N',ndim,m,1D0,work(1,V_ind),ndim,Y,1,1D0,xvec,1)
#else
    X = X + MATMUL( work(:,V_ind+1-1:V_ind+m-1), Y(1:m) )
#endif


500 CONTINUE

    !
    ! Check the convergence
    !

    select case (HUTI_STOPC)
    case (HUTI_TRUERESIDUAL)
       call pcondrsubr( T1V, X, ipar )
       call matvecsubr( T1V, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,T1V_ind) = rhsvec(ii) - work(ii,R_ind)
       END DO
       !$OMP END PARALLEL DO
#else
       T1V = B - R
#endif
       call pcondlsubr( R, T1V, ipar )
       residual = normfun( HUTI_NDIM, R, 1 )
    case (HUTI_TRESID_SCALED_BYB)
       call pcondrsubr( T1V, X, ipar )
       call matvecsubr( T1V, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,T1V_ind) = rhsvec(ii) - work(ii,R_ind)
       END DO
       !$OMP END PARALLEL DO
#else
       T1V = B - R
#endif
       call pcondlsubr( R, T1V, ipar )
       residual = normfun( HUTI_NDIM, R, 1 ) / rhsnorm
    case (HUTI_PSEUDORESIDUAL)
       call matvecsubr( X, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,R_ind) = work(ii,R_ind) - rhsvec(ii)
       END DO
       !$OMP END PARALLEL DO
#else
       R = R - B
#endif
       call pcondlsubr( T1V, R, ipar )
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_PRESID_SCALED_BYB)
       call matvecsubr( X, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,R_ind) = work(ii,R_ind) - rhsvec(ii)
       END DO
       !$OMP END PARALLEL DO
#else
       R = R - B
#endif
       call pcondlsubr( T1V, R, ipar )
       residual = normfun( HUTI_NDIM, T1V, 1 ) / rhsnorm
    case (HUTI_PRESID_SCALED_BYPRECB)
       call matvecsubr( X, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,R_ind) = work(ii,R_ind) - rhsvec(ii)
       END DO
       !$OMP END PARALLEL DO
#else
       R = R - B
#endif
       call pcondlsubr( T1V, R, ipar )
       residual = normfun( HUTI_NDIM, T1V, 1 ) / precrhsnorm
    case (HUTI_XDIFF_NORM)
#ifdef _OPENMP
       CALL DGEMV('N',ndim,m,1D0,work(1,V_ind),ndim,Y,1,0D0,work(1,T1V_ind),1)
#else
       T1V = MATMUL( work(:,V_ind+1-1:V_ind+m-1), Y(1:m) )
#endif
       residual = normfun( HUTI_NDIM, T1V, 1 )
    case (HUTI_USUPPLIED_STOPC)
       residual = stopcfun( X, B, R, ipar, dpar )
    case default
       call pcondrsubr( T1V, X, ipar )
       call matvecsubr( T1V, R, ipar )
#ifdef _OPENMP
       !$OMP PARALLEL DO
       DO ii=1,ndim
         work(ii,T1V_ind) = work(ii,R_ind) - rhsvec(ii)
       END DO
       !$OMP END PARALLEL DO
#else
       T1V = R - B
#endif
       call pcondlsubr( R, T1V, ipar )
       residual = normfun( HUTI_NDIM, R, 1 )
    end select

    work(m+1,S_ind) = normfun( HUTI_NDIM, R, 1 )

    !
    ! Print debugging output if desired
    !

    if ( HUTI_DBUGLVL .ne. HUTI_NO_DEBUG ) then
       if ( mod(iter_count, HUTI_DBUGLVL) .eq. 0 ) then
          WRITE (*, '(A, I8, E11.4)') '   gmres:',iter_count, residual
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
       WRITE (*, '(A, I8, E11.4)') '   gmres:',iter_count, residual
    end if

    HUTI_ITERS = iter_count
    return

    ! End of execution
    !*********************************************************************

  end subroutine  huti_dgmressolv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_cgmressolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    integer :: iter_count
    real :: residual, rhsnorm, precrhsnorm
    real :: bnrm, alpha, beta

    INTEGER :: reit, info, i, j, k, l, m, n

    complex :: temp, temp2, error, gamma

    ! Local arrays

    complex, DIMENSION(HUTI_GMRES_RESTART+1,HUTI_GMRES_RESTART+1) :: H
    complex, &
         DIMENSION((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1)) :: HLU
    complex, DIMENSION(HUTI_GMRES_RESTART+1) :: CS, SN, Y

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
       call  huti_crandvec   ( X, ipar )
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

          call  huti_clusolve  ( i, HLU, Y, work(:,S_ind) )

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

    call  huti_clusolve  ( m, HLU, Y, work(:,S_ind) )

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
          WRITE (*, '(A, I8, E11.4)') '   gmresz:',iter_count, residual
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
       WRITE (*, '(A, I8, E11.4)') '   gmresz:', iter_count, residual
    end if

    HUTI_ITERS = iter_count
    return

    ! End of execution
    !*********************************************************************

  end subroutine  huti_cgmressolv

  !*************************************************************************


  !*************************************************************************
  !*************************************************************************
  ! Double complex version
  !*************************************************************************
  !*************************************************************************

  subroutine  huti_zgmressolv  ( ndim, wrkdim, xvec, rhsvec, &
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

    integer :: iter_count
    double precision :: residual, rhsnorm, precrhsnorm
    double precision :: bnrm, alpha, beta

    INTEGER :: reit, info, i, j, k, l, m, n

    double complex :: temp, temp2, error, gamma

    ! Local arrays

    double complex, DIMENSION(HUTI_GMRES_RESTART+1,HUTI_GMRES_RESTART+1) :: H
    double complex, &
         DIMENSION((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1)) :: HLU
    double complex, DIMENSION(HUTI_GMRES_RESTART+1) :: CS, SN, Y

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
       call  huti_zrandvec   ( X, ipar )
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

          call  huti_zlusolve  ( i, HLU, Y, work(:,S_ind) )

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

    call  huti_zlusolve  ( m, HLU, Y, work(:,S_ind) )

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
          WRITE (*, '(A, I8, E11.4)') '   gmresz:',iter_count, residual
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
       WRITE (*, '(A,I8, E11.4)') '   gmresz:',iter_count, residual
    end if

    HUTI_ITERS = iter_count
    return

    ! End of execution
    !*********************************************************************

  end subroutine  huti_zgmressolv

  !*************************************************************************


end module huti_gmres
