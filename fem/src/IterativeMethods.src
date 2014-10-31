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
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20.9.2007
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Module containing iterative methods. Uses the calling procedure of the 
!>  HUTIter package for similar interfacing. The idea is that the future 
!>  development of iterative methods could be placed in this module. 
!------------------------------------------------------------------------------


#include "huti_fdefs.h"

! if using old huti_fdefs.h, later obsolite
#ifndef HUTI_MAXTOLERANCE
#define HUTI_MAXTOLERANCE dpar(2)
#endif
#ifndef HUTI_SGSPARAM
#define HUTI_SGSPARAM dpar(3)
#endif
#ifndef HUTI_BICGSTABL_L
#define HUTI_BICGSTABL_L ipar(16)
#endif
#ifndef HUTI_DIVERGENCE
#define HUTI_DIVERGENCE 3
#endif
#ifndef HUTI_GCR_RESTART
#define HUTI_GCR_RESTART ipar(17)
#endif

MODULE IterativeMethods
  
  USE Types
  USE CRSMatrix  

  IMPLICIT NONE
  
  INTEGER :: nc
  LOGICAL :: Constrained

  TYPE(Matrix_t), POINTER, PRIVATE :: CM
  
CONTAINS
  

!------------------------------------------------------------------------------
!> Symmetric Gauss-Seidel iterative method for linear systems. This is not really of practical
!> use but may be used for testing, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_sgs( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
        
    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    REAL(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp) :: &
      xvec(HUTI_NDIM),rhsvec(HUTI_NDIM),work(HUTI_WRKDIM,HUTI_NDIM)

    INTEGER :: ndim
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual, Omega
    LOGICAL :: Converged, Diverged

    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL 
    Omega = HUTI_SGSPARAM

    CALL sgs(ndim, GlobalMatrix, xvec, rhsvec, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval, Omega)

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER
    
  CONTAINS 
 
!------------------------------------------------------------------------------
    SUBROUTINE SGS( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval, Omega )
 !------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      INTEGER :: i,j,k,n
      REAL(KIND=dp) :: x(n),b(n)
      REAL(KIND=dp) :: Omega
      INTEGER, POINTER :: Cols(:),Rows(:)
      REAL(KIND=dp), POINTER :: Values(:)
      REAL(KIND=dp), ALLOCATABLE :: R(:)
      LOGICAL :: Converged, Diverged
      INTEGER :: OutputInterval
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual, bnorm,rnorm,w,s

      Rows   => A % Rows
      Cols   => A % Cols
      Values => A % Values
      
      ALLOCATE( R(n) )
      
      CALL matvecsubr( x, r, ipar )
     
      r(1:n) = b(1:n) - r(1:n)
      bnorm = normfun(n, b, 1)
      rnorm = normfun(n, r, 1) 

      Residual = rnorm / bnorm
      Converged = (Residual < MinTolerance) 
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
      IF( Converged .OR. Diverged) RETURN

      DO k=1,Rounds
        DO i=1,n
          s = 0.0d0
          DO j=Rows(i),Rows(i+1)-1
            s = s + x(Cols(j)) * Values(j)
          END DO
          x(i) = x(i) + Omega * (b(i)-s) / Values(A % Diag(i))
        END DO
        
        DO i=n,1,-1
          s = 0.0d0
          DO j=Rows(i),Rows(i+1)-1
            s = s + x(Cols(j)) * Values(j)
          END DO
          x(i) = x(i) + Omega * (b(i)-s) / Values(A % Diag(i))
        END DO
        
        CALL matvecsubr( x, r, ipar )
        r(1:n) = b(1:n) - r(1:n)
        rnorm = normfun(n, r, 1)
        
        Residual = rnorm / bnorm
        IF( MOD(k,OutputInterval) == 0) THEN
          WRITE (*, '(I8, 2E11.4)') k, rnorm, residual
        END IF
        
        Converged = (Residual < MinTolerance) 
        Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
        IF( Converged .OR. Diverged) RETURN
        
      END DO
    END SUBROUTINE SGS
!------------------------------------------------------------------------------
  END SUBROUTINE itermethod_sgs
 !------------------------------------------------------------------------------
 


!------------------------------------------------------------------------------
!> Jacobi iterative method for linear systems. This is not really of practical
!> use but may be used for testing, for example. 
!> Note that if the scaling is performed so that the diagonal entry is one
!> the division by it is unnecessary. Hence for this method scaling is not 
!> needed. 
!------------------------------------------------------------------------------
 SUBROUTINE itermethod_jacobi( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
       
    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    REAL(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp) :: &
       xvec(HUTI_NDIM),rhsvec(HUTI_NDIM),work(HUTI_WRKDIM,HUTI_NDIM)

    INTEGER :: ndim
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged
    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL 
       
    CALL jacobi(ndim, GlobalMatrix, xvec, rhsvec, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval )

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER   

  CONTAINS 
    
    
    SUBROUTINE Jacobi( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      REAL(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: OutputInterval
      REAL(KIND=dp) :: bnorm,rnorm
      REAL(KIND=dp), ALLOCATABLE :: R(:)
!------------------------------------------------------------------------------
      INTEGER :: i,j,n
!------------------------------------------------------------------------------
      
      Converged = .FALSE.
      Diverged = .FALSE.
      
      ALLOCATE( R(n) )
      
      CALL matvecsubr( x, r, ipar )
      r(1:n) = b(1:n) - r(1:n)
      
      bnorm = normfun(n, b, 1)
      rnorm = normfun(n, r, 1)
      
      Residual = rnorm / bnorm
      Converged = (Residual < MinTolerance) 
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)    
      IF( Converged .OR. Diverged) RETURN
      
      DO i=1,Rounds
        DO j=1,n
          x(j) = x(j) + r(j) / A % Values(A % diag(j))
        END DO
        CALL matvecsubr( x, r, ipar )
        
        r(1:n) = b(1:n) - r(1:n)
        rnorm = normfun(n, r, 1)
        
        Residual = rnorm / bnorm
        
        IF( MOD(i,OutputInterval) == 0) THEN
          WRITE (*, '(I8, 2E11.4)') i, rnorm, residual
        END IF
        
        Converged = (Residual < MinTolerance) 
        Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
        IF( Converged .OR. Diverged) EXIT
      END DO
      
      DEALLOCATE( R )
      
    END SUBROUTINE Jacobi
    
!------------------------------------------------------------------------------
  END SUBROUTINE itermethod_jacobi
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> Richardson iterative method for linear systems. This may of actual use for 
!> mass matrices. Actually this is not the simple Richardson iteration method
!> as it is preconditioned with the lumped mass matrix. 
!> Note that if scaling is performed by the "row equilibrium" method then
!> lumped mass is by construction unity (assuming all-positive entries).
!> So for this method scaling is not needed. 
!------------------------------------------------------------------------------
 SUBROUTINE itermethod_richardson( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
       
    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    REAL(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp) :: &
       xvec(HUTI_NDIM),rhsvec(HUTI_NDIM),work(HUTI_WRKDIM,HUTI_NDIM)

    INTEGER :: ndim
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged
    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL 
       
    CALL richardson(ndim, GlobalMatrix, xvec, rhsvec, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval )

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER   

  CONTAINS 
    
    
    SUBROUTINE Richardson( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      REAL(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: OutputInterval
      REAL(KIND=dp) :: bnorm,rnorm,s,q
      INTEGER, POINTER :: Cols(:),Rows(:)
      REAL(KIND=dp), POINTER :: Values(:)
      REAL(KIND=dp), ALLOCATABLE :: R(:), M(:)
!------------------------------------------------------------------------------
      INTEGER :: i,j,k,n
!------------------------------------------------------------------------------

      Rows   => A % Rows
      Cols   => A % Cols
      Values => A % Values
      
      Converged = .FALSE.
      Diverged = .FALSE.
      
      ALLOCATE( R(n), M(n) )
      
      CALL matvecsubr( x, r, ipar )
      r(1:n) = b(1:n) - r(1:n)
      
      bnorm = normfun(n, b, 1)
      rnorm = normfun(n, r, 1)

      Residual = rnorm / bnorm
      Converged = (Residual < MinTolerance) 
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
      IF( Converged .OR. Diverged) RETURN

      ! Perform preconditioning by mass lumping
      DO i=1,n
        s = 0.0_dp
        DO j=Rows(i),Rows(i+1)-1
          s = s + Values( j )
        END DO
        M(i) = s 
      END DO

      DO k=1,Rounds
        DO i=1,n
          IF( k == 1 ) THEN
            x(i) = b(i) / M(i) 
          ELSE
            x(i) = x(i) + r(i) / M(i)
          END IF
        END DO
        
        CALL matvecsubr( x, r, ipar )

        r(1:n) = b(1:n) - r(1:n)
        rnorm = normfun(n, r, 1)
        
        Residual = rnorm / bnorm
        
        IF( MOD(k,OutputInterval) == 0) THEN
          WRITE (*, '(I8, 2E11.4)') k, rnorm, residual
        END IF
        
        Converged = (Residual < MinTolerance) 
        Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
        IF( Converged .OR. Diverged) EXIT
      END DO
      
      DEALLOCATE( R, M )
      
    END SUBROUTINE Richardson
    
!------------------------------------------------------------------------------
  END SUBROUTINE itermethod_richardson
!------------------------------------------------------------------------------
  

!-----------------------------------------------------------------------------------
    SUBROUTINE C_matvec(u,v,ipar,matvecsubr)
!-----------------------------------------------------------------------------------
      INTEGER :: ipar(*)
      REAL(KIND=dp) :: u(*),v(*)

      EXTERNAL matvecsubr
!-----------------------------------------------------------------------------------
      INTEGER :: i,j,k,l,ndim
!-----------------------------------------------------------------------------------
      ndim = HUTI_NDIM
      CALL matvecsubr( u, v, ipar )
      IF (Constrained) THEN
         DO i=1,CM % NumberOfRows
           k = ndim+i
           v(k) = 0._dp
           DO j = CM % Rows(i),CM % Rows(i+1)-1
             l = CM % Cols(j)
             v(l) = v(l) + CM % Values(j)*u(k)
             v(k) = v(k) + CM % Values(j)*u(l)
           END DO
         END DO
       END IF
!-----------------------------------------------------------------------------------
    END SUBROUTINE C_matvec
!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------
    SUBROUTINE C_lpcond(u,v,ipar,pcondlsubr)
!-----------------------------------------------------------------------------------
      INTEGER :: ipar(*)
      REAL(KIND=dp) :: u(*),v(*)

      EXTERNAL pcondlsubr
!-----------------------------------------------------------------------------------
      INTEGER :: ndim
!-----------------------------------------------------------------------------------
      ndim = HUTI_NDIM
      IF(Constrained) HUTI_NDIM = ndim+nc
      CALL pcondlsubr(u,v,ipar)
      IF(Constrained) HUTI_NDIM = ndim
!-----------------------------------------------------------------------------------
    END SUBROUTINE C_lpcond
!-----------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>   This routine solves real linear systems Ax = b by using the BiCGStab(l) algorithm 
!>   with l >= 2 and the right-oriented ILU(n) preconditioning. 
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_bicgstabl( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
    

    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    REAL(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), TARGET :: &
       xvec(HUTI_NDIM),rhsvec(HUTI_NDIM),work(HUTI_WRKDIM,HUTI_NDIM)

    INTEGER :: ndim,i,j,k
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, UseStopCFun

    TYPE(Matrix_t),POINTER :: A

    REAL(KIND=dp), POINTER :: x(:),b(:)

    A => GlobalMatrix
    CM => A % ConstraintMatrix
    Constrained = ASSOCIATED(CM)
    
    ndim = HUTI_NDIM

    x => xvec
    b => rhsvec
    nc = 0
    IF (Constrained) THEN
      nc = CM % NumberOfRows
      Constrained = nc>0
      IF(Constrained) THEN
        ALLOCATE(x(ndim+nc),b(ndim+nc))
        IF(.NOT.ALLOCATED(CM % ExtraVals))THEN
          ALLOCATE(CM % ExtraVals(nc)); CM % extraVals=0._dp
        END IF
        b(1:ndim) = rhsvec; b(ndim+1:) = CM % RHS
        x(1:ndim) = xvec; x(ndim+1:) = CM % extraVals
      END IF
    END IF

    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL 
    PolynomialDegree = HUTI_BICGSTABL_L 
    UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC

    CALL RealBiCGStabl(ndim+nc, A,x,b, Rounds, MinTol, MaxTol, &
        Converged, Diverged, OutputInterval, PolynomialDegree )

    IF(Constrained) THEN
      xvec=x(1:ndim)
      rhsvec=b(1:ndim)
      CM % extraVals = x(ndim+1:ndim+nc)
      DEALLOCATE(x,b)
    END IF

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER

  CONTAINS



!-----------------------------------------------------------------------------------
!>   The subroutine has been written using as a starting point the work of D.R. Fokkema 
!>   (subroutine zbistbl v1.1 1998). Dr. Fokkema has given the right to distribute
!>   the derived work under GPL and hence the original more conservative 
!> copyright notice of the subroutine has been removed accordingly.  
!-----------------------------------------------------------------------------------
    SUBROUTINE RealBiCGStabl( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
        Diverged, OutputInterval, l, StoppingCriterionType )
!----------------------------------------------------------------------------------- 
      INTEGER :: l   ! polynomial degree
      INTEGER :: n, MaxRounds, OutputInterval   
      LOGICAL :: Converged, Diverged
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp) :: x(n), b(n)
      REAL(KIND=dp) :: Tol, MaxTol
      INTEGER, OPTIONAL :: StoppingCriterionType 
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: zero, one, t(n), kappa0, kappal 
      REAL(KIND=dp) :: dnrm2, rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
          delta = 1.0d-2, bnrm, bw_errorind, tottime
      INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
          IluOrder
      REAL(KIND=dp) :: alpha, beta, omega, rho0, rho1, sigma, ddot, varrho, hatgamma
      LOGICAL rcmp, xpdt, GotIt, BackwardError, EarlyExit
      CHARACTER(LEN=MAX_NAME_LEN) :: str
      REAL(KIND=dp), ALLOCATABLE :: work(:,:), rwork(:,:)
!------------------------------------------------------------------------------
    
      IF ( l < 2) CALL Fatal( 'RealBiCGStabl', 'Polynomial degree < 2' )
      
      IF ( ALL(x == 0.0d0) ) x = b

      zero = 0.0d0
      one =  1.0d0

      ALLOCATE( work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)) )
      work = 0.0d0
      rwork = 0.0d0
    
      rr = 1
      r = rr+1
      u = r+(l+1)
      xp = u+(l+1)
      bp = xp+1
    
      z = 1
      zz = z+(l+1)
      y0 = zz+(l+1)
      yl = y0+1
      y = yl+1
    
      CALL C_matvec(x,work(:,r),ipar,matvecsubr)

      work(1:n,r) = b(1:n) - work(1:n,r)
      bnrm  = normfun(n, b(1:n), 1 )
      rnrm0 = normfun(n, work(1:n,r), 1 )

      !-------------------------------------------------------------------
      ! Check whether the initial guess is already converged, diverged or NaN
      !--------------------------------------------------------------------
      Diverged = .FALSE.
      IF (bnrm /= bnrm) THEN
        CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: bnrm = NaN.' )
      ENDIF
      IF(rnrm0 /= rnrm0 ) THEN
        CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: nrm0 = NaN.' )
      END IF

      errorind = rnrm0 / bnrm
      IF(errorind /= errorind ) THEN
        CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: errorind = NaN.' )
      END IF
     
      Converged = (errorind < Tol)
      Diverged = (errorind > MaxTol) 

      IF( Converged .OR. Diverged) RETURN

      EarlyExit = .FALSE.

      work(1:n,rr) = work(1:n,r) 
      work(1:n,bp) = work(1:n,r)
      work(1:n,xp) = x(1:n)

      rnrm = rnrm0
      mxnrmx = rnrm0
      mxnrmr = rnrm0      
      x(1:n) = zero      
      alpha = zero
      omega = one
      sigma = one
      rho0 = one

      DO Round=1,MaxRounds
        !-------------------------
        ! --- The BiCG part ---
        !-------------------------
        rho0 = -omega*rho0
      
        DO k=1,l
          rho1 = dotprodfun(n, work(1:n,rr), 1, work(1:n,r+k-1), 1 )
          IF (rho0 == zero) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: rho0 == zero.' )
          ENDIF
          IF (rho1 /= rho1) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: rho1 == NaN.' )
          ENDIF
         
          beta = alpha*(rho1/rho0)
          rho0 = rho1
          DO j=0,k-1
            work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
          ENDDO

          CALL C_lpcond( t, work(:,u+k-1), ipar, pcondlsubr )
          CALL C_matvec( t, work(:,u+k), ipar, matvecsubr )
          sigma = dotprodfun(n, work(1:n,rr), 1, work(1:n,u+k), 1 )
          
          IF (sigma == zero) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: sigma == zero.' )
          ENDIF
          IF (sigma /= sigma) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: sigma == NaN.' )
          ENDIF

          alpha = rho1/sigma
          x(1:n) = x(1:n) + alpha * work(1:n,u)
          DO j=0,k-1
            work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
          ENDDO

          CALL C_lpcond( t, work(:,r+k-1), ipar,pcondlsubr )
          CALL C_matvec( t, work(:,r+k), ipar, matvecsubr )

          rnrm = normfun(n, work(1:n,r), 1 )
          IF (rnrm /= rnrm) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: rnrm == NaN.' )
          ENDIF

          mxnrmx = MAX (mxnrmx, rnrm)
          mxnrmr = MAX (mxnrmr, rnrm)

          !----------------------------------------------------------------------
          ! In some simple cases, a few BiCG updates may already be enough to
          ! obtain the solution. The following is for handling this special case. 
          !----------------------------------------------------------------------
          errorind = rnrm / bnrm

!         IF( OutputInterval /= 0) THEN
!           WRITE (*, '(I8, 2E11.4)') 0, rnrm, errorind
!         END IF

          Converged = (errorind < Tol) 
          Diverged = (errorind /= errorind)

          IF (Converged .OR. Diverged) THEN
             EarlyExit = .TRUE.
             EXIT
          END IF
        END DO

        IF (EarlyExit) EXIT

        !--------------------------------------
        ! --- The convex polynomial part ---
        !--------------------------------------
        
        DO i=1,l+1
          DO j=1,i
            rwork(i,j) = dotprodfun(n, work(1:n,r+i-1), 1, work(1:n,r+j-1), 1 ) 
          END DO
        END DO
        DO j=2,l+1
          rwork(1:j-1,j) = rwork(j,1:j-1)
        END DO
          
        rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
        CALL dgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
            iwork, stat)
      
        ! --- tilde r0 and tilde rl (small vectors)
      
        rwork(1,y0) = -one
        rwork(2:l,y0) = rwork(2:l,z) 
        CALL dgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
            rwork(2:l,y0), l-1, stat)
        rwork(l+1,y0) = zero
        
        rwork(1,yl) = zero
        rwork(2:l,yl) = rwork(2:l,z+l) 
        CALL dgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
            rwork(2:l,yl), l-1, stat)
        rwork(l+1,yl) = -one
      
        ! --- Convex combination
      
        CALL dsymv ('u', l+1, one, rwork(1:l+1,z:z+l), l+1, &
            rwork(1:l+1,y0), 1, zero, rwork(1:l+1,y), 1)
        kappa0 = SQRT( ddot(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1) )
        CALL dsymv ('u', l+1, one, rwork(1:l+1,z:z+l), l+1, &
            rwork(1:l+1,yl), 1, zero, rwork(1:l+1,y), 1)
        kappal = SQRT( ddot(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) )
        CALL dsymv ('u', l+1, one, rwork(1:l+1,z:z+l), l+1, &
          rwork(1:l+1,y0), 1, zero, rwork(1:l+1,y), 1)
        varrho = ddot(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &
            (kappa0*kappal)
        hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
            kappa0/kappal
        rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)
        
        !  --- Update
        
        omega = rwork(l+1,y0)
        DO j=1,l
          work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
          x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
          work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
        ENDDO
    
        CALL dsymv ('u', l+1, one, rwork(1:l+1,z:z+l), l+1, &
            rwork(1:l+1,y0), 1, zero, rwork(1:l+1,y), 1)
        rnrm = SQRT( ddot(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1) )
        
        !---------------------------------------
        !  --- The reliable update part ---
        !---------------------------------------
        
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          CALL C_lpcond( t, x, ipar,pcondlsubr )
          CALL C_matvec( t, work(:,r), ipar, matvecsubr )

          work(1:n,r) = work(1:n,bp) - work(1:n,r)
          mxnrmr = rnrm
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            work(1:n,xp) = work(1:n,xp) + t(1:n)
            x(1:n) = zero
            work(1:n,bp) = work(1:n,r)
            mxnrmx = rnrm
          ENDIF
        ENDIF
        
        IF (rcmp) THEN
          IF (xpdt) THEN       
            t(1:n) = work(1:n,xp)
          ELSE
            t(1:n) = t(1:n) + work(1:n,xp)  
          END IF
        ELSE
          CALL C_lpcond( t, x, ipar,pcondlsubr )
          t(1:n) = t(1:n)+work(1:n,xp)
        END IF
      
        errorind = rnrm / bnrm

        IF( MOD(Round,OutputInterval) == 0) THEN
          WRITE (*, '(I8, 2E11.4)') Round, rnrm, errorind
        END IF
    
        Converged = (errorind < Tol) 
        Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
        IF( Converged .OR. Diverged) EXIT    
      END DO

      IF(OutputInterval /= HUGE(OutputInterval)) THEN
        WRITE (*, '(I8, 2E11.4)') Round, rnrm, errorind
      END IF
    
      !------------------------------------------------------------
      ! We have solved z = P*x, with P the preconditioner, so finally 
      ! solve the true unknown x
      !------------------------------------------------------------
      t(1:n) = x(1:n)
      CALL C_lpcond( x, t, ipar,pcondlsubr )
      x(1:n) = x(1:n) + work(1:n,xp)        
!------------------------------------------------------------------------------
    END SUBROUTINE RealBiCGStabl
!------------------------------------------------------------------------------

  END SUBROUTINE itermethod_bicgstabl
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   This routine solves real linear systems Ax = b by using the GCR algorithm 
!> (Generalized Conjugate Residual).
!------------------------------------------------------------------------------
 SUBROUTINE itermethod_gcr( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
       
    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    REAL(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), TARGET :: &
       xvec(HUTI_NDIM),rhsvec(HUTI_NDIM),work(HUTI_WRKDIM,HUTI_NDIM)

    INTEGER :: ndim, RestartN
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, UseStopCFun

    TYPE(Matrix_t),POINTER::A

    REAL(KIND=dp), POINTER :: x(:),b(:)
    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    RestartN = HUTI_GCR_RESTART 
    UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC

    x => xvec
    b => rhsvec
    nc = 0

    A => GlobalMatrix
    CM => A % ConstraintMatrix
    Constrained = ASSOCIATED(CM)

    IF (Constrained) THEN
      nc = CM % NumberOfRows
      Constrained = nc>0
      IF(Constrained) THEN
        ALLOCATE(x(ndim+nc),b(ndim+nc))
        IF(.NOT.ALLOCATED(CM % ExtraVals))THEN
          ALLOCATE(CM % ExtraVals(nc)); CM % extraVals=0._dp
        END IF
        b(1:ndim) = rhsvec; b(ndim+1:) = CM % RHS
        x(1:ndim) = xvec; x(ndim+1:) = CM % extraVals
      END IF
    END IF
       
    CALL GCR(ndim+nc, GlobalMatrix, x, b, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval, RestartN )

    IF(Constrained) THEN
      xvec=x(1:ndim)
      rhsvec=b(1:ndim)
      CM % extraVals = x(ndim+1:ndim+nc)
      DEALLOCATE(x,b)
    END IF

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER   

  CONTAINS 
    
    
    SUBROUTINE GCR( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval, m) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      REAL(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: n, OutputInterval, m
      REAL(KIND=dp) :: bnorm,rnorm
      REAL(KIND=dp), ALLOCATABLE :: R(:)

      REAL(KIND=dp), ALLOCATABLE :: S(:,:), V(:,:), T1(:), T2(:),TT(:)

!------------------------------------------------------------------------------
      INTEGER :: i,j,k
      REAL(KIND=dp) :: alpha, beta, trueres(n), trueresnorm, normerr
!------------------------------------------------------------------------------
      
      Converged = .FALSE.
      Diverged = .FALSE.
      
      ALLOCATE( R(n), T1(n), T2(n),TT(n) )
      IF ( m > 1 ) THEN
         ALLOCATE( S(n,m-1), V(n,m-1) )
         V(1:n,1:m-1) = 0.0d0	
         S(1:n,1:m-1) = 0.0d0
      END IF	
      
      CALL C_matvec( x, r, ipar, matvecsubr )
      r(1:n) = b(1:n) - r(1:n)
      
      bnorm = normfun(n, b, 1)
      rnorm = normfun(n, r, 1)
      
      IF (UseStopCFun) THEN
        Residual = stopcfun(x,b,r,ipar,dpar)
      ELSE
        Residual = rnorm / bnorm
      END IF
      Converged = (Residual < MinTolerance) 
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)
      IF( Converged .OR. Diverged) RETURN
       
      DO k=1,Rounds
	 !----------------------------------------------
	 ! Check for restarting
         !----------------------------------------------
         IF ( MOD(k,m)==0 ) THEN
            j = m
         ELSE
            j = MOD(k,m)
            !--------------------------------------------
            ! Compute the true residual when restarting:
            !--------------------------------------------
            IF ( (j==1) .AND. (k>1) ) THEN
               CALL C_matvec( x, r, ipar, matvecsubr )
               r(1:n) = b(1:n) - r(1:n)
            END IF
         END IF

         !----------------------------------------------------------
         ! Perform the preconditioning...
         !---------------------------------------------------------------
         CALL C_lpcond( T1, r, ipar,pcondlsubr )
         CALL C_matvec( T1, T2, ipar, matvecsubr )

         !--------------------------------------------------------------
         ! Perform the orthogonalization of the search directions....
         !--------------------------------------------------------------
         DO i=1,j-1
            beta = dotprodfun(n, V(1:n,i), 1, T2(1:n), 1 )
            T1(1:n) = T1(1:n) - beta * S(1:n,i)
            T2(1:n) = T2(1:n) - beta * V(1:n,i)        
         END DO

         alpha = normfun(n, T2(1:n), 1 )
         T1(1:n) = 1.0d0/alpha * T1(1:n)
         T2(1:n) = 1.0d0/alpha * T2(1:n)

         !-------------------------------------------------------------
         ! The update of the solution and save the search data...
         !------------------------------------------------------------- 
         beta = dotprodfun(n, T2(1:n), 1, r(1:n), 1 )
         x(1:n) = x(1:n) + beta * T1(1:n)      
         r(1:n) = r(1:n) - beta * T2(1:n)
	 IF ( j /= m ) THEN
            S(1:n,j) = T1(1:n)
            V(1:n,j) = T2(1:n)
	 END IF       

         !--------------------------------------------------------------
         ! Check whether the convergence criterion is met 
         !--------------------------------------------------------------
         rnorm = normfun(n, r, 1)
         !CALL C_matvec( x, trueres, ipar, matvecsubr )
         !trueres(1:n) = b(1:n) - trueres(1:n)
         !trueresnorm = normfun(n, trueres, 1)
         IF (UseStopCFun) THEN
           Residual = stopcfun(x,b,r,ipar,dpar)
           IF( MOD(k,OutputInterval) == 0) THEN
             WRITE (*, '(I8, 3E11.4)') k, rnorm / bnorm, residual
           END IF           
         ELSE
           Residual = rnorm / bnorm
           IF( MOD(k,OutputInterval) == 0) THEN
             WRITE (*, '(I8, E11.4)') k, residual
           END IF
         END IF
           
         Converged = (Residual < MinTolerance) 
         !-----------------------------------------------------------------
         ! Make an additional check that the true residual agrees with 
         ! the iterated residual:
         !-----------------------------------------------------------------
         IF (Converged ) THEN
            CALL C_matvec( x, trueres, ipar, matvecsubr )
            trueres(1:n) = b(1:n) - trueres(1:n)
            TrueResNorm = normfun(n, trueres, 1)
            NormErr = ABS(TrueResNorm - rnorm)/TrueResNorm
            IF ( NormErr > 1.0d-1 ) THEN
               CALL Info('WARNING', 'Iterated GCR solution may not be accurate', Level=2)
               WRITE( Message, * ) 'Iterated GCR residual norm = ', rnorm
               CALL Info('WARNING', Message, Level=2)
               WRITE( Message, * ) 'True residual norm = ', TrueResNorm
               CALL Info('WARNING', Message, Level=2)   
            END IF
         END IF
         Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)    
         IF( Converged .OR. Diverged) EXIT
        
      END DO
      
      DEALLOCATE( R, T1, T2 )
      IF ( m > 1 ) DEALLOCATE( S, V)
      
    END SUBROUTINE GCR
    
!------------------------------------------------------------------------------
  END SUBROUTINE itermethod_gcr
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
!> This routine provides the complex version to the GCR linear solver.
!------------------------------------------------------------------------------
 SUBROUTINE itermethod_z_gcr( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
!------------------------------------------------------------------------------
       
    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    COMPLEX(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun
    
    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp) :: &
       xvec(2*HUTI_NDIM),rhsvec(2*HUTI_NDIM),work(HUTI_WRKDIM,2*HUTI_NDIM)

    COMPLEX(KIND=dp) :: y(HUTI_NDIM),f(HUTI_NDIM)
    INTEGER :: ndim, RestartN, i
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged
    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    RestartN = HUTI_GCR_RESTART 

    !----------------------------------------------------------------------------
    ! Transform the solution vector and the right-hand side vector to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------    
    DO i=1,ndim
      y(i) = CMPLX( xvec(2*i-1), xvec(2*i), kind=dp )
      f(i) = CMPLX( rhsvec(2*i-1), rhsvec(2*i), kind=dp )
    END DO

       
    CALL GCR_Z(ndim, GlobalMatrix, y, f, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval, RestartN )

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER
   
    !----------------------------------------------
    ! Return the solution as a real vector...
    !----------------------------------------------
    DO i=1,ndim
      xvec( 2*i-1 ) = REAL( y(i) )
      xvec( 2*i ) = AIMAG( y(i) )
    END DO



  CONTAINS 
    
    
    SUBROUTINE GCR_Z( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval, m) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      COMPLEX(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: n, OutputInterval, m
      REAL(KIND=dp) :: bnorm,rnorm
      COMPLEX(KIND=dp), ALLOCATABLE :: R(:)

      COMPLEX(KIND=dp), ALLOCATABLE :: S(:,:), V(:,:), T1(:), T2(:)

!------------------------------------------------------------------------------
      INTEGER :: i,j,k
      COMPLEX(KIND=dp) :: beta
      REAL(KIND=dp) :: alpha, trueresnorm, normerr
      COMPLEX(KIND=dp) :: trueres(n)
!------------------------------------------------------------------------------
      
      Converged = .FALSE.
      Diverged = .FALSE.
      
      ALLOCATE( R(n), T1(n), T2(n) )
      IF ( m > 1 ) THEN
         ALLOCATE( S(n,m-1), V(n,m-1) )
         V(1:n,1:m-1) = CMPLX( 0.0d0, 0.0d0, kind=dp)
         S(1:n,1:m-1) = CMPLX( 0.0d0, 0.0d0, kind=dp)
      END IF	
      
      CALL matvecsubr( x, r, ipar )
      r(1:n) = b(1:n) - r(1:n)
      
      bnorm = normfun(n, b, 1)
      rnorm = normfun(n, r, 1)
      
      Residual = rnorm / bnorm
      Converged = (Residual < MinTolerance) 
      Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)    
      IF( Converged .OR. Diverged) RETURN
      
      DO k=1,Rounds
	 !----------------------------------------------
	 ! Check for restarting
         !--------------------------------------------- 
         IF ( MOD(k,m)==0 ) THEN
            j = m
         ELSE
            j = MOD(k,m)
            !--------------------------------------------
            ! Compute the true residual when restarting:
            !--------------------------------------------
            IF ( (j==1) .AND. (k>1) ) THEN
               CALL matvecsubr( x, r, ipar ) 
               r(1:n) = b(1:n) - r(1:n)
            END IF
         END IF
         !----------------------------------------------------------
         ! Perform the preconditioning...
         !---------------------------------------------------------------
         CALL pcondlsubr( T1, r, ipar )         
         CALL matvecsubr( T1, T2, ipar )
         !--------------------------------------------------------------
         ! Perform the orthogonalization of the search directions....
         !--------------------------------------------------------------
         DO i=1,j-1
            beta = dotprodfun(n, V(1:n,i), 1, T2(1:n), 1 )
            T1(1:n) = T1(1:n) - beta * S(1:n,i)
            T2(1:n) = T2(1:n) - beta * V(1:n,i)        
         END DO

         alpha = normfun(n, T2(1:n), 1 )
         T1(1:n) = 1.0d0/alpha * T1(1:n)
         T2(1:n) = 1.0d0/alpha * T2(1:n)

         !-------------------------------------------------------------
         ! The update of the solution and save the search data...
         !------------------------------------------------------------- 
         beta = dotprodfun(n, T2(1:n), 1, r(1:n), 1 )
         x(1:n) = x(1:n) + beta * T1(1:n)      
         r(1:n) = r(1:n) - beta * T2(1:n)
	 IF ( j /= m ) THEN
            S(1:n,j) = T1(1:n)
            V(1:n,j) = T2(1:n)
	 END IF       

         !--------------------------------------------------------------
         ! Check whether the convergence criterion is met 
         !--------------------------------------------------------------
         rnorm = normfun(n, r, 1)
         Residual = rnorm / bnorm
        
         IF( MOD(k,OutputInterval) == 0) THEN
            WRITE (*, '(I8, E11.4)') k, residual
         END IF
        
         Converged = (Residual < MinTolerance)
         !-----------------------------------------------------------------
         ! Make an additional check that the true residual agrees with 
         ! the iterated residual:
         !-----------------------------------------------------------------
         IF (Converged ) THEN
            CALL matvecsubr( x, trueres, ipar )
            trueres(1:n) = b(1:n) - trueres(1:n)
            TrueResNorm = normfun(n, trueres, 1)
            NormErr = ABS(TrueResNorm - rnorm)/TrueResNorm
            IF ( NormErr > 1.0d-1 ) THEN
               CALL Info('WARNING', 'Iterated GCR solution may not be accurate', Level=2)
               WRITE( Message, * ) 'Iterated GCR residual norm = ', rnorm
               CALL Info('WARNING', Message, Level=2)
               WRITE( Message, * ) 'True residual norm = ', TrueResNorm
               CALL Info('WARNING', Message, Level=2)   
            END IF
         END IF 
         Diverged = (Residual > MaxTolerance) .OR. (Residual /= Residual)    
         IF( Converged .OR. Diverged) EXIT
        
      END DO
      
      DEALLOCATE( R, T1, T2 )
      IF ( m > 1 ) DEALLOCATE( S, V)
      
    END SUBROUTINE GCR_Z
    
!------------------------------------------------------------------------------
  END SUBROUTINE itermethod_z_gcr
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> This routine provides a complex version of the BiCGstab(l) solver for linear systems.
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_z_bicgstabl( xvec, rhsvec, &
       ipar, dpar, work, matvecsubr, pcondlsubr, &
       pcondrsubr, dotprodfun, normfun, stopcfun )
 !------------------------------------------------------------------------------

    IMPLICIT NONE

    EXTERNAL matvecsubr, pcondlsubr, pcondrsubr
    EXTERNAL dotprodfun, normfun, stopcfun
    COMPLEX(KIND=dp) :: dotprodfun
    REAL(KIND=dp) :: normfun
    REAL(KIND=dp) :: stopcfun

    ! Parameters
    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp) :: &
         xvec(2*HUTI_NDIM),rhsvec(2*HUTI_NDIM),work(HUTI_WRKDIM,2*HUTI_NDIM)

    COMPLEX(KIND=dp) :: y(HUTI_NDIM),f(HUTI_NDIM)
    INTEGER :: ndim, i, PolynomialDegree
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol
    LOGICAL :: Converged, Diverged
    !--------------------------------------------------------------------------------    

    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    PolynomialDegree = HUTI_BICGSTABL_L 
    !----------------------------------------------------------------------------
    ! Transform the solution vector and the right-hand side vector to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------    
    DO i=1,ndim
       y(i) = CMPLX( xvec(2*i-1), xvec(2*i), kind=dp )
       f(i) = CMPLX( rhsvec(2*i-1), rhsvec(2*i), kind=dp )
    END DO

    CALL ComplexBiCGStabl(ndim, GlobalMatrix, y, f, Rounds, MinTol, MaxTol, &
         Converged, Diverged, OutputInterval, PolynomialDegree)

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER

    !----------------------------------------------
    ! Return the solution as a real vector...
    !----------------------------------------------
    DO i=1,ndim
       xvec( 2*i-1 ) = REAL( y(i) )
       xvec( 2*i ) = AIMAG( y(i) )
    END DO


  CONTAINS 

    !-----------------------------------------------------------------------------------
    SUBROUTINE ComplexBiCGStabl( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
         Diverged, OutputInterval, l, StoppingCriterionType )
    !-----------------------------------------------------------------------------------
    !   This subroutine solves complex linear systems Ax = b by using the BiCGStab(l) algorithm 
    !   with l >= 2 and the right-oriented preconditioning. 
    !
    !   The subroutine has been written using as a starting point the work of D.R. Fokkema 
    !   (subroutine zbistbl v1.1 1998). Dr. Fokkema has given the right to distribute
    !   the derived work under GPL and hence the original copyright notice of the subroutine
    !   has been removed accordingly.  
    !
    !----------------------------------------------------------------------------------- 
      INTEGER :: l   ! polynomial degree
      INTEGER :: n, MaxRounds, OutputInterval   
      LOGICAL :: Converged, Diverged
      TYPE(Matrix_t), POINTER :: A
      COMPLEX(KIND=dp) :: x(n), b(n)
      REAL(KIND=dp) :: Tol, MaxTol
      INTEGER, OPTIONAL :: StoppingCriterionType 
      !------------------------------------------------------------------------------
      COMPLEX(KIND=dp) :: zzero, zone, t(n), kappa0, kappal 
      REAL(KIND=dp) :: rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
           delta = 1.0d-2, bnrm
      INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round
      COMPLEX(KIND=dp) :: alpha, beta, omega, rho0, rho1, sigma, zdotc, varrho, hatgamma
      COMPLEX(KIND=dp), ALLOCATABLE :: work(:,:), rwork(:,:)
      LOGICAL rcmp, xpdt, EarlyExit
      !------------------------------------------------------------------------------

      IF ( l < 2) CALL Fatal( 'RealBiCGStabl', 'Polynomial degree < 2' )

      IF ( ALL(x == CMPLX(0.0d0,0.0d0,kind=dp)) ) x = b

      zzero = CMPLX( 0.0d0,0.0d0, KIND=dp)
      zone =  CMPLX( 1.0d0,0.0d0, KIND=dp)

      ALLOCATE( work(n,3+2*(l+1)), rwork(l+1,3+2*(l+1)) )
      work = CMPLX( 0.0d0, 0.0d0, KIND=dp )
      rwork = CMPLX( 0.0d0, 0.0d0, KIND=dp )

      rr = 1
      r = rr+1
      u = r+(l+1)
      xp = u+(l+1)
      bp = xp+1

      z = 1
      zz = z+(l+1)
      y0 = zz+(l+1)
      yl = y0+1
      y = yl+1

      CALL matvecsubr( x, work(1:n,r), ipar )
      work(1:n,r) = b(1:n) - work(1:n,r)
      bnrm = normfun(n, b(1:n), 1)
      rnrm0 = normfun(n, work(1:n,r), 1)

      !-------------------------------------------------------------------
      ! Check whether the initial guess satisfies the stopping criterion
      !--------------------------------------------------------------------
      errorind = rnrm0 / bnrm
      Converged = (errorind < Tol)
      Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)

      IF( Converged .OR. Diverged) RETURN
      EarlyExit = .FALSE.

      work(1:n,rr) = work(1:n,r) 
      work(1:n,bp) = work(1:n,r)
      work(1:n,xp) = x(1:n)

      rnrm = rnrm0
      mxnrmx = rnrm0
      mxnrmr = rnrm0
      x(1:n) = zzero    
      alpha = zzero
      omega = zone
      sigma = zone
      rho0 = zone

      DO Round=1,MaxRounds 
         !-------------------------
         ! --- The BiCG part ---
         !-------------------------
         rho0 = -omega*rho0

         DO k=1,l
            rho1 = dotprodfun(n, work(1:n,rr), 1, work(1:n,r+k-1), 1)
            IF (rho0 == zzero) THEN
               CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
            ENDIF
            beta = alpha*(rho1/rho0)
            rho0 = rho1
            DO j=0,k-1
               work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
            ENDDO
            CALL pcondlsubr( t, work(1:n,u+k-1), ipar )
            CALL matvecsubr( t, work(1:n,u+k),   ipar )

            sigma = dotprodfun(n, work(1:n,rr), 1, work(1:n,u+k), 1)
            IF (sigma == zzero) THEN
               CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error.' )
            ENDIF
            alpha = rho1/sigma
            x(1:n) = x(1:n) + alpha * work(1:n,u)
            DO j=0,k-1
               work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
            ENDDO
            CALL pcondlsubr( t, work(1:n,r+k-1), ipar )
            CALL matvecsubr( t, work(1:n,r+k),   ipar )
            rnrm = normfun(n, work(1:n,r), 1)
            mxnrmx = MAX (mxnrmx, rnrm)
            mxnrmr = MAX (mxnrmr, rnrm)

            !----------------------------------------------------------------------
            ! In some simple cases, a few BiCG updates may already be enough to
            ! obtain the solution. The following is for handling this special case. 
            !----------------------------------------------------------------------
            errorind = rnrm / bnrm
            Converged = (errorind < Tol) 
            IF (Converged) THEN
               EarlyExit = .TRUE.
               EXIT
            END IF

         ENDDO
         
         IF (EarlyExit) EXIT        

         !--------------------------------------
         ! --- The convex polynomial part ---
         !--------------------------------------

         DO i=1,l+1
            DO j=1,i
               rwork(i,j) = dotprodfun(n, work(1:n,r+i-1), 1, work(1:n,r+j-1),1 ) 
            END DO
         END DO
         DO j=2,l+1
            rwork(1:j-1,j) = CONJG( rwork(j,1:j-1) )
         END DO

         rwork(1:l+1,zz:zz+l) = rwork(1:l+1,z:z+l)
         CALL zgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
              iwork, stat)

         ! --- tilde r0 and tilde rl (small vectors)

         rwork(1,y0) = -zone
         rwork(2:l,y0) = rwork(2:l,z) 
         CALL zgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
              rwork(2:l,y0), l-1, stat)
         rwork(l+1,y0) = zzero

         rwork(1,yl) = zzero
         rwork(2:l,yl) = rwork(2:l,z+l) 
         CALL zgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
              rwork(2:l,yl), l-1, stat)
         rwork(l+1,yl) = -zone

         ! --- Convex combination

         CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
              rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
         kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) ) ! replace zdotc
         CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
              rwork(1:l+1,yl), 1, zzero, rwork(1:l+1,y), 1)
         kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )  ! replace zdotc
         CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
              rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
         varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &            ! replace zdotc
              (kappa0*kappal)
         hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
              kappa0/kappal
         rwork(1:l+1,y0) = rwork(1:l+1,y0) - hatgamma * rwork(1:l+1,yl)

         !  --- Update

         omega = rwork(l+1,y0)
         DO j=1,l
            work(1:n,u) = work(1:n,u) - rwork(j+1,y0) * work(1:n,u+j)
            x(1:n) = x(1:n) + rwork(j+1,y0) * work(1:n,r+j-1)
            work(1:n,r) = work(1:n,r) - rwork(j+1,y0) * work(1:n,r+j)
         ENDDO

         CALL zhemv ('u', l+1, zone, rwork(1:l+1,z:z+l), l+1, &
              rwork(1:l+1,y0), 1, zzero, rwork(1:l+1,y), 1)
         rnrm = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )

         !---------------------------------------
         !  --- The reliable update part ---
         !---------------------------------------

         mxnrmx = MAX (mxnrmx, rnrm)
         mxnrmr = MAX (mxnrmr, rnrm)
         xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
         rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
         IF (rcmp) THEN
            ! PRINT *, 'Performing residual update...'
            CALL pcondlsubr( t, x, ipar )
            CALL matvecsubr( t, work(1:n,r), ipar )
            work(1:n,r) = work(1:n,bp) - work(1:n,r)
            mxnrmr = rnrm
            IF (xpdt) THEN
               ! PRINT *, 'Performing solution update...'
               work(1:n,xp) = work(1:n,xp) + t(1:n)
               x(1:n) = zzero
               work(1:n,bp) = work(1:n,r)
               mxnrmx = rnrm
            ENDIF
         ENDIF

         IF (rcmp) THEN
            IF (xpdt) THEN       
               t(1:n) = work(1:n,xp)
            ELSE
               t(1:n) = t(1:n) + work(1:n,xp)  
            END IF
         ELSE
            CALL pcondlsubr( t, x, ipar )
            t(1:n) = t(1:n)+work(1:n,xp)
         END IF

         errorind = rnrm/bnrm
         IF( MOD(Round,OutputInterval) == 0) THEN
            WRITE (*, '(I8, E11.4)') Round, errorind
         END IF

         Converged = (errorind < Tol) 
         Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)    
         IF( Converged .OR. Diverged) EXIT    
      END DO

      IF( EarlyExit .AND. (OutputInterval/=HUGE(OutputInterval)) ) THEN
         WRITE (*, '(I8, E11.4)') Round, errorind         
      END IF

      !------------------------------------------------------------
      ! We have solved z = P*x, so finally solve the true unknown x
      !------------------------------------------------------------
      t(1:n) = x(1:n)
      CALL pcondlsubr( x, t, ipar )
      x(1:n) = x(1:n) + work(1:n,xp)      

    !----------------------------------------------------------
    END SUBROUTINE ComplexBiCGStabl
    !----------------------------------------------------------

!--------------------------------------------------------------
  END SUBROUTINE itermethod_z_bicgstabl
!--------------------------------------------------------------

END MODULE IterativeMethods

!> \} ElmerLib

