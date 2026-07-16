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
! *  Authors: Juha Ruokolainen, Peter Råback, Mika Malinen, Martin van Gijzen
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

! if using old huti_fdefs.h, later obsolete
#ifndef HUTI_MAXTOLERANCE
#define HUTI_MAXTOLERANCE dpar(2)
#endif
#ifndef HUTI_SGSPARAM
#define HUTI_SGSPARAM dpar(3)
#endif
#ifndef HUTI_PSEUDOCOMPLEX
#define HUTI_PSEUDOCOMPLEX ipar(7)
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
#ifndef HUTI_IDRS_S
#define HUTI_IDRS_S ipar(18)
#endif


MODULE IterativeMethods
  
  USE Types
  USE CRSMatrix  
  USE SParIterComm
  
  IMPLICIT NONE
  
  INTEGER :: nc
  LOGICAL :: Constrained

  TYPE(Matrix_t), POINTER, PRIVATE :: CM
  
CONTAINS
  

  ! When treating a complex system with iterative solver, norm and matrix-vector product are
  ! similar for real-valued and complex-valued systems. However, the inner product is different.
  ! For pseudo-complex systems this routine generates also the complex part of the product.
  ! This may have a favourable effect on convergence.
  !
  ! This routine has same API as the fully real-valued system but every second call returns
  ! the missing complex part.
  !
  ! This routine assumes that in x and y the values follow each other. 
  !-----------------------------------------------------------------------------------
  FUNCTION PseudoZDotProd( ndim, x, xind, y, yind ) RESULT( d )
  !-----------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: ndim, xind, yind
    REAL(KIND=dp) :: x(*)
    REAL(KIND=dp) :: y(*)
    REAL(KIND=dp) :: d
        
    INTEGER :: i, callcount = 0
    REAL(KIND=dp) :: a,b
    
    SAVE callcount, a, b

    IF( callcount == 0 ) THEN    
      ! z = x^H*y = (x_re-i*x_im)(y_re+i*y_im)       
      ! =>  z_re = x_re*y_re + x_im*y_im
      !     z_im = x_re*y_im - x_im*y_re
      
      a = SUM( x(1:ndim) * y(1:ndim) )
      b = SUM( x(1:ndim:2) * y(2:ndim:2) - x(2:ndim:2) * y(1:ndim:2) )

      IF (ParEnv % PEs > 1) THEN
        CALL SParActiveSUM(a,0)
        CALL SParActiveSUM(b,0)
      END IF

      d = a 
      callcount = callcount + 1
    ELSE
      d = b
      callcount = 0
    END IF
      
    !-----------------------------------------------------------------------------------
  END FUNCTION PseudoZDotProd
  !-----------------------------------------------------------------------------------

  
  ! As the previous but assumes that the real and complex values are ordered blockwise. 
  !-----------------------------------------------------------------------------------
  FUNCTION PseudoZDotProd2( ndim, x, xind, y, yind ) RESULT( d )
  !-----------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: ndim, xind, yind
    REAL(KIND=dp) :: x(*)
    REAL(KIND=dp) :: y(*)
    REAL(KIND=dp) :: d
        
    INTEGER :: i, callcount = 0
    REAL(KIND=dp) :: a,b
    
    SAVE callcount, a, b

    IF( callcount == 0 ) THEN    
      a = SUM( x(1:ndim) * y(1:ndim) )
      b = SUM( x(1:ndim/2) * y(ndim/2+1:ndim) - x(ndim/2+1:ndim) * y(1:ndim/2) )

      IF (ParEnv % PEs > 1) THEN
        CALL SParActiveSUM(a,0)
        CALL SParActiveSUM(b,0)
      END IF
      
      d = a 
      callcount = callcount + 1
    ELSE
      d = b
      callcount = 0
    END IF
      
!-----------------------------------------------------------------------------------
  END FUNCTION PseudoZDotProd2
!-----------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> Symmetric Gauss-Seidel iterative method for linear systems. This is not really of practical
!> use but may be used for testing, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_sgs( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )

    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

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
          CALL FLUSH(6)
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

    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

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
          CALL FLUSH(6)
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

    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

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
          CALL FLUSH(6)
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
      USE huti_interfaces
      IMPLICIT NONE
      PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
      INTEGER :: ipar(*)
      REAL(KIND=dp) :: u(*),v(*)
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
    RECURSIVE SUBROUTINE C_rpcond(u,v,ipar,pcondrsubr)
!-----------------------------------------------------------------------------------
      USE huti_interfaces
      IMPLICIT NONE
      PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
      INTEGER :: ipar(*)
      REAL(KIND=dp) :: u(*),v(*)

!-----------------------------------------------------------------------------------
      INTEGER :: ndim
!-----------------------------------------------------------------------------------
      ndim = HUTI_NDIM
      IF(Constrained) HUTI_NDIM = ndim+nc
      CALL pcondrsubr(u,v,ipar)
      IF(Constrained) HUTI_NDIM = ndim
!-----------------------------------------------------------------------------------
    END SUBROUTINE C_rpcond
!-----------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>   This routine solves real linear systems Ax = b by using the BiCGStab(l) algorithm 
!>   with l >= 2 and the right-oriented ILU(n) preconditioning. 
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_bicgstabl( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )

    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    ! DOUBLE PRECISION, DIMENSION(HUTI_NDIM), TARGET :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work
    ! DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    INTEGER :: ndim,i,j,k
    INTEGER :: Rounds, OutputInterval, PolynomialDegree
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, Halted, UseStopCFun, PseudoComplex

    TYPE(Matrix_t),POINTER :: A

    REAL(KIND=dp), POINTER CONTIG :: x(:),b(:)

    ! Variables related to robust mode
    LOGICAL :: Robust 
    INTEGER :: BestIter,BadIterCount,MaxBadIter, RobustStart
    REAL(KIND=dp) :: BestNorm,RobustStep,RobustTol,RobustMaxTol
    REAL(KIND=dp), ALLOCATABLE :: Bestx(:)

    
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

    PseudoComplex = ( HUTI_PSEUDOCOMPLEX > 0 )  
    
    Converged = .FALSE.
    Diverged = .FALSE.
    Halted = .FALSE.
    
    Robust = ( HUTI_ROBUST == 1 )
    IF( Robust ) THEN
      RobustTol = HUTI_ROBUST_TOLERANCE
      RobustStep = HUTI_ROBUST_STEPSIZE
      RobustMaxTol = HUTI_ROBUST_MAXTOLERANCE
      MaxBadIter = HUTI_ROBUST_MAXBADIT
      RobustStart = HUTI_ROBUST_START
      BestNorm = SQRT(HUGE(BestNorm))
      BadIterCount = 0
      BestIter = 0      
      ALLOCATE( BestX(ndim))
    END IF
    
    CALL RealBiCGStabl(ndim+nc, A,x,b, Rounds, MinTol, MaxTol, &
         Converged, Diverged, Halted, OutputInterval, PolynomialDegree )

    IF(Constrained) THEN
      xvec=x(1:ndim)
      rhsvec=b(1:ndim)
      CM % extraVals = x(ndim+1:ndim+nc)
      DEALLOCATE(x,b)
    END IF
    
    IF( Robust ) THEN
      DEALLOCATE( BestX )
    END IF
      
    IF(Converged) THEN
      HUTI_INFO = HUTI_CONVERGENCE
    ELSE IF(Diverged) THEN
      HUTI_INFO = HUTI_DIVERGENCE
    ELSE IF(Halted) THEN
      HUTI_INFO = HUTI_HALTED
    ELSE
      HUTI_INFO = HUTI_MAXITER
    END IF
      
  CONTAINS

!-----------------------------------------------------------------------------------
!>   The subroutine has been written using as a starting point the work of D.R. Fokkema 
!>   (subroutine zbistbl v1.1 1998). Dr. Fokkema has given the right to distribute
!>   the derived work under GPL and hence the original more conservative 
!> copyright notice of the subroutine has been removed accordingly.  
!-----------------------------------------------------------------------------------
    SUBROUTINE RealBiCGStabl( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
        Diverged, Halted, OutputInterval, l)
!----------------------------------------------------------------------------------- 
      INTEGER :: l   ! polynomial degree
      INTEGER :: n, MaxRounds, OutputInterval   
      LOGICAL :: Converged, Diverged, Halted
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp) :: x(n), b(n)
      REAL(KIND=dp) :: Tol, MaxTol
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: zero, one, t(n), kappa0, kappal 
      REAL(KIND=dp) :: dnrm2, rnrm0, rnrm, mxnrmx, mxnrmr, errorind, &
          delta = 1.0d-2, bnrm, bw_errorind, tottime
      INTEGER :: i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, iwork(l-1), stat, Round, &
          IluOrder
      REAL(KIND=dp) :: alpha, beta, omega, rho0, rho1, sigma, ddot, varrho, hatgamma
      LOGICAL rcmp, xpdt, GotIt, EarlyExit
      REAL(KIND=dp), ALLOCATABLE :: work(:,:)
      REAL(KIND=dp) :: rwork(l+1,3+2*(l+1))
      REAL(KIND=dp) :: tmpmtr(l-1,l-1), tmpvec(l-1)
      REAL(KIND=dp) :: beta_im
!------------------------------------------------------------------------------
    
      IF ( l < 2) CALL Fatal( 'RealBiCGStabl', 'Polynomial degree < 2' )
      
      IF ( ALL(x == 0.0d0) ) x = b

      zero = 0.0d0
      one =  1.0d0

      ALLOCATE( work(n,3+2*(l+1)) )
      !$OMP PARALLEL PRIVATE(j)
      DO j=1,3+2*(l+1)
         !$OMP DO SCHEDULE(STATIC)
         DO i=1,n
            work(i,j) = 0.0d0
         END DO
         !$OMP END DO
      END DO
      !$OMP END PARALLEL
      DO j=1,3+2*(l+1)
         DO i=1,l+1
            rwork(i,j) = 0.0d0
         END DO
      END DO

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
    
      ! CALL C_matvec(x,work(:,r),ipar,matvecsubr)
      CALL C_matvec(x,work(1,r),ipar,matvecsubr)
      
      !$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i=1,n
         work(i,r) = b(i) - work(i,r)
      END DO
      !$OMP END PARALLEL DO
      ! bnrm  = normfun(n, b(1:n), 1 )
      ! rnrm0 = normfun(n, work(1:n,r), 1 )
      bnrm  = normfun(n, b(1), 1 )
      rnrm0 = normfun(n, work(1,r), 1 )

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

      IF( Converged .OR. Diverged ) RETURN

      EarlyExit = .FALSE.

      !$OMP PARALLEL 
      !$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         work(i,rr) = work(i,r) 
         work(i,bp) = work(i,r)
      END DO
      !$OMP END DO NOWAIT
      !$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         work(i,xp) = x(i)
      END DO
      !$OMP END DO NOWAIT
      !$OMP DO SCHEDULE(STATIC)
      DO i=1,n
         x(i) = zero
      END DO
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL

      rnrm = rnrm0
      mxnrmx = rnrm0
      mxnrmr = rnrm0  
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
          ! rho1 = dotprodfun(n, work(1:n,rr), 1, work(1:n,r+k-1), 1 )
          rho1 = dotprodfun(n, work(1,rr), 1, work(1,r+k-1), 1 )
          IF (rho0 == zero) THEN
            CALL Warn( 'RealBiCGStab(l)', 'Iteration halted: rho0 == zero.' )
            Halted = .TRUE.
            GOTO 100
          ENDIF
          IF (rho1 /= rho1) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: rho1 == NaN.' )
          ENDIF
         
          beta = alpha*(rho1/rho0)
          rho0 = rho1
          !$OMP PARALLEL PRIVATE(j)
          DO j=0,k-1
             !$OMP DO SCHEDULE(STATIC)
             DO i=1,n
                work(i,u+j) = work(i,r+j) - beta*work(i,u+j)
             END DO
             !$OMP END DO
          ENDDO
          !$OMP END PARALLEL

          ! CALL C_rpcond( t, work(:,u+k-1), ipar, pcondrsubr )
          CALL C_rpcond( t, work(1,u+k-1), ipar, pcondrsubr )
          ! CALL C_matvec( t, work(:,u+k), ipar, matvecsubr )
          CALL C_matvec( t, work(1,u+k), ipar, matvecsubr )
          ! sigma = dotprodfun(n, work(1:n,rr), 1, work(1:n,u+k), 1 )
          sigma = dotprodfun(n, work(1,rr), 1, work(1,u+k), 1 )
          
          IF (sigma == zero) THEN
            CALL Warn( 'RealBiCGStab(l)', 'Iteration halted: sigma == zero.' )
            Halted = .TRUE.
            GOTO 100
          ENDIF
          IF (sigma /= sigma) THEN
            CALL Fatal( 'RealBiCGStab(l)', 'Breakdown error: sigma == NaN.' )
          ENDIF
          
          alpha = rho1/sigma

          !$OMP PARALLEL PRIVATE(j)
          !$OMP DO SCHEDULE(STATIC)
          DO i=1,n
             x(i) = x(i) + alpha * work(i,u)
          END DO
          !$OMP END DO NOWAIT
          DO j=0,k-1
             !$OMP DO SCHEDULE(STATIC)
             DO i=1,n
                work(i,r+j) = work(i,r+j) - alpha * work(i,u+j+1)
             END DO
             !$OMP END DO
          ENDDO
          !$OMP END PARALLEL

          ! CALL C_rpcond( t, work(:,r+k-1), ipar,pcondrsubr )
          ! CALL C_matvec( t, work(:,r+k), ipar, matvecsubr )
          CALL C_rpcond( t, work(1,r+k-1), ipar,pcondrsubr )
          CALL C_matvec( t, work(1,r+k), ipar, matvecsubr )

          ! rnrm = normfun(n, work(1:n,r), 1 )
          rnrm = normfun(n, work(1,r), 1 )

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
             ! rwork(i,j) = dotprodfun(n, work(1:n,r+i-1), 1, work(1:n,r+j-1), 1 ) 
             rwork(i,j) = dotprodfun(n, work(1,r+i-1), 1, work(1,r+j-1), 1 ) 
          END DO
        END DO
        DO j=2,l+1
           DO i=1,j-1
              rwork(i,j) = rwork(j,i)
           END DO
        END DO
        DO j=0,l-1
           DO i=1,l+1
              rwork(i,zz+j) = rwork(i,z+j)
           END DO
        END DO
        DO j=1,l-1
           DO i=1,l-1
              tmpmtr(i,j) = rwork(i+1,zz+j)
           END DO
        END DO
        ! CALL dgetrf (l-1, l-1, rwork(2:l,zz+1:zz+l-1), l-1, &
        !     iwork, stat)
        CALL dgetrf (l-1, l-1, tmpmtr, l-1, &
             iwork, stat)
      
        ! --- tilde r0 and tilde rl (small vectors)
        
        rwork(1,y0) = -one
        DO i=2,l
           rwork(i,y0) = rwork(i,z)
        END DO
        DO i=1,l-1
           tmpvec(i) = rwork(i+1,y0)
        END DO
        ! CALL dgetrs('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
        !     rwork(2:l,y0), l-1, stat)
        CALL dgetrs('n', l-1, 1, tmpmtr, l-1, iwork, &
             tmpvec, l-1, stat)
        DO i=1,l-1
           rwork(i+1,y0) = tmpvec(i)
        END DO
        rwork(l+1,y0) = zero
        
        rwork(1,yl) = zero
        DO i=1,l-1
           rwork(i+1,yl) = rwork(i+1,z+l)
           tmpvec(i) = rwork(i+1,yl)
        END DO
        ! CALL dgetrs ('n', l-1, 1, rwork(2:l,zz+1:zz+l-1), l-1, iwork, &
        !     rwork(2:l,yl), l-1, stat)
        CALL dgetrs ('n', l-1, 1, tmpmtr, l-1, iwork, &
             tmpvec, l-1, stat)
        DO i=1,l-1
           rwork(i+1,yl) = tmpvec(i)
        END DO
        rwork(l+1,yl) = -one
      
        ! --- Convex combination          
        
        CALL dsymv ('u', l+1, one, rwork(1,z), l+1, &
            rwork(1,y0), 1, zero, rwork(1,y), 1)
        kappa0 = ddot(l+1, rwork(1,y0), 1, rwork(1,y), 1)

        ! If untreated this would result to NaN's
        IF( kappa0 <= 0.0 ) THEN
          CALL Warn('RealBiCGStab(l)','kappa0^2 is non-positive, iteration halted')
          Halted = .TRUE.
          GOTO 100
        END IF
        kappa0 = SQRT( kappa0 ) 

        CALL dsymv ('u', l+1, one, rwork(1,z), l+1, &
            rwork(1,yl), 1, zero, rwork(1,y), 1)
        kappal = ddot(l+1, rwork(1,yl), 1, rwork(1,y), 1 )
        
        ! If untreated this would result to NaN's
        IF( kappal <= 0.0 ) THEN
          CALL Warn('RealBiCGStab(l)','kappal^2 is non-positive, iteration halted')
          Halted = .TRUE.
          GOTO 100 
        END IF
        kappal = SQRT( kappal )

        CALL dsymv ('u', l+1, one, rwork(1,z), l+1, &
          rwork(1,y0), 1, zero, rwork(1,y), 1)

        varrho = ddot(l+1, rwork(1,yl), 1, rwork(1,y), 1) / &
            (kappa0*kappal)
        
        hatgamma = varrho/ABS(varrho) * MAX(ABS(varrho),7d-1) * &
            kappa0/kappal
        DO i=1,l+1
           rwork(i,y0) = rwork(i,y0) - hatgamma * rwork(i,yl)
        END DO
        !  --- Update
        
        omega = rwork(l+1,y0)
        !$OMP PARALLEL PRIVATE(j,i) FIRSTPRIVATE(rwork)
        DO j=1,l
           !$OMP DO SCHEDULE(STATIC)
           DO i=1,n
              work(i,u) = work(i,u) - rwork(j+1,y0) * work(i,u+j)
           END DO
           !$OMP END DO
           !$OMP DO SCHEDULE(STATIC)
           DO i=1,n
              x(i) = x(i) + rwork(j+1,y0) * work(i,r+j-1)
           END DO
           !$OMP END DO
           !$OMP DO SCHEDULE(STATIC)
           DO i=1,n
              work(i,r) = work(i,r) - rwork(j+1,y0) * work(i,r+j)
           END DO
           !$OMP END DO
        ENDDO
        !$OMP END PARALLEL
    
        CALL dsymv ('u', l+1, one, rwork(1,z), l+1, &
            rwork(1,y0), 1, zero, rwork(1,y), 1)
        rnrm = ddot(l+1, rwork(1,y0), 1, rwork(1,y), 1)

        IF( rnrm < 0.0 ) THEN
          CALL Warn('RealBiCGStab(l)','rnrm^2 is negative, iteration halted')
          Halted = .TRUE.
          GOTO 100 
        END IF        
        rnrm = SQRT( rnrm ) 
        
        !---------------------------------------
        !  --- The reliable update part ---
        !---------------------------------------
        
        mxnrmx = MAX (mxnrmx, rnrm)
        mxnrmr = MAX (mxnrmr, rnrm)
        xpdt = (rnrm < delta*rnrm0 .AND. rnrm0 < mxnrmx)
        rcmp = ((rnrm < delta*mxnrmr .AND. rnrm0 < mxnrmr) .OR. xpdt)
        IF (rcmp) THEN
          ! PRINT *, 'Performing residual update...'
          CALL C_rpcond( t, x, ipar,pcondrsubr )
          ! CALL C_matvec( t, work(:,r), ipar, matvecsubr )
          CALL C_matvec( t, work(1,r), ipar, matvecsubr )

          mxnrmr = rnrm
          !$OMP PARALLEL DO SCHEDULE(STATIC)
          DO i=1,n
             work(i,r) = work(i,bp) - work(i,r)
          END DO
          !$OMP END PARALLEL DO
          IF (xpdt) THEN
            ! PRINT *, 'Performing solution update...'
            !$OMP PARALLEL DO SCHEDULE(STATIC)
             DO i=1,n
                work(i,xp) = work(i,xp) + t(i)
                x(i) = zero
                work(i,bp) = work(i,r)
             END DO
             !$OMP END PARALLEL DO

             mxnrmx = rnrm
          ENDIF
        ENDIF
        
        IF (rcmp) THEN
          IF (xpdt) THEN       
             !$OMP PARALLEL DO SCHEDULE(STATIC)
             DO i=1,n
                t(i) = work(i,xp)
             END DO
             !$OMP END PARALLEL DO
          ELSE
             !$OMP PARALLEL DO SCHEDULE(STATIC)
             DO i=1,n
                t(i) = t(i) + work(i,xp)  
             END DO
             !$OMP END PARALLEL DO
          END IF
        ELSE
          CALL C_rpcond( t, x, ipar,pcondrsubr )
          !$OMP PARALLEL DO
          DO i=1,n
             t(i) = t(i)+work(i,xp)
          END DO
          !$OMP END PARALLEL DO
        END IF
      
        errorind = rnrm / bnrm

        IF( MOD(Round,OutputInterval) == 0) THEN
          WRITE (*, '(I8, 2E11.4)') Round, rnrm, errorind
          CALL FLUSH(6)
        END IF
        
        IF( Robust ) THEN
          IF (Round>=RobustStart ) THEN
            IF( errorInd < RobustStep * BestNorm ) THEN
              BestIter = Round
              BestNorm = errorInd
              Bestx = x
              BadIterCount = 0
            ELSE
              BadIterCount = BadIterCount + 1
            END IF

            IF( BestNorm <  RobustTol .AND. &
                  ( errorInd > RobustMaxTol .OR. BadIterCount > MaxBadIter ) ) THEN
              EXIT
            END IF
          END IF
        END IF
               
        Converged = (errorind < Tol) 
        Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
        IF( Converged .OR. Diverged) EXIT    
      END DO

100   IF( Robust ) THEN
        IF( BestNorm < RobustTol ) THEN
          Converged = .TRUE.
        END IF
        IF( BestNorm < errorInd ) THEN
          x = Bestx
        END IF
        IF(OutputInterval /= HUGE(OutputInterval)) THEN
          WRITE(*,'(A,I8,E11.4,I8,2E11.4)') 'BiCGStabl robust: ',&
              MIN(MaxRounds,Round), BestNorm, BestIter, rnrm, errorind
          CALL FLUSH(6)
        END IF
      ELSE
        IF(OutputInterval /= HUGE(OutputInterval)) THEN
          WRITE (*, '(A, I8, 2E11.4)') 'BiCGStabl: ', MIN(MaxRounds,Round), rnrm, errorind
          CALL FLUSH(6)
        END IF
      END IF
            
      !------------------------------------------------------------
      ! We have solved z = P*x, with P the preconditioner, so finally 
      ! solve the true unknown x
      !------------------------------------------------------------
      !$OMP PARALLEL DO
      DO i=1,n
         t(i) = x(i)
      END DO
      !$OMP END PARALLEL DO
      CALL C_rpcond( x, t, ipar,pcondrsubr )
      !$OMP PARALLEL DO
      DO i=1,n
         x(i) = x(i) + work(i,xp)
      END DO
      !$OMP END PARALLEL DO
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

    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    INTEGER :: ndim, RestartN
    INTEGER :: Rounds, MinIter, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, UseStopCFun
    LOGICAL :: PseudoComplex

    TYPE(Matrix_t),POINTER::A

    REAL(KIND=dp), POINTER :: x(:),b(:)

    CALL Info('Itermetod_gcr','Starting GCR iteration',Level=25)

    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinIter = HUTI_MINIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    RestartN = HUTI_GCR_RESTART 
    UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC

    Converged = .FALSE.
    Diverged = .FALSE.
    PseudoComplex = ( HUTI_PSEUDOCOMPLEX > 0 )  
      
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
        Converged, Diverged, OutputInterval, RestartN, MinIter )

    
    IF(Constrained) THEN
      xvec = x(1:ndim)
      rhsvec = b(1:ndim)
      CM % extraVals = x(ndim+1:ndim+nc)
      DEALLOCATE(x,b)
    END IF

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER   

  CONTAINS 
    
    
    SUBROUTINE GCR( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval, m, MinIter) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds, MinIter
      REAL(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: n, OutputInterval, m
      REAL(KIND=dp) :: bnorm,rnorm
      REAL(KIND=dp), ALLOCATABLE :: R(:)

      REAL(KIND=dp), ALLOCATABLE :: S(:,:), V(:,:), T1(:), T2(:)

!------------------------------------------------------------------------------
      INTEGER :: i,j,k,ksum=0
      REAL(KIND=dp) :: alpha, beta, trueres(n), trueresnorm, normerr
      REAL(KIND=dp) :: beta_im
!------------------------------------------------------------------------------
      INTEGER :: allocstat
        
      ALLOCATE( R(n), T1(n), T2(n), STAT=allocstat )
      IF( allocstat /= 0 ) THEN
        CALL Fatal('GCR','Failed to allocate memory of size: '//I2S(n))
      END IF

      IF ( m > 1 ) THEN
        ALLOCATE( S(n,m-1), V(n,m-1), STAT=allocstat )
        IF( allocstat /= 0 ) THEN
          CALL Fatal('GCR','Failed to allocate memory of size: '&
              //I2S(n)//' x '//I2S(m-1))
        END IF
        
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
      Converged = (Residual < MinTolerance) .AND. ( MinIter <= 0 )
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
         CALL C_rpcond( T1, r, ipar, pcondrsubr )
         CALL C_matvec( T1, T2, ipar, matvecsubr )

         !--------------------------------------------------------------
         ! Perform the orthogonalization of the search directions....
         !--------------------------------------------------------------
         DO i=1,j-1
           beta = dotprodfun(n, V(1:n,i), 1, T2(1:n), 1 )

           IF( PseudoComplex ) THEN
             ! The even call is for the complex part of beta
             ! This has to be before the T1 and T2 vectors are tampered
             ! For convenience we subtract 
             beta_im = dotprodfun(n, V(1:n,i), 1, T2(1:n), 1 )

             IF( HUTI_PSEUDOCOMPLEX == 2 ) THEN
               T1(1:n/2) = T1(1:n/2) + beta_im * S(n/2+1:n,i) 
               T1(n/2+1:n) = T1(n/2+1:n) - beta_im * S(1:n/2,i)                    
               
               T2(1:n/2) = T2(1:n/2) + beta_im * V(1+n/2:n,i)
               T2(1+n/2:n) = T2(1+n/2:n) - beta_im * V(1:n/2,i)                                
             ELSE
               T1(1:n:2) = T1(1:n:2) + beta_im * S(2:n:2,i) 
               T1(2:n:2) = T1(2:n:2) - beta_im * S(1:n:2,i)                    
               
               T2(1:n:2) = T2(1:n:2) + beta_im * V(2:n:2,i)
               T2(2:n:2) = T2(2:n:2) - beta_im * V(1:n:2,i)
             END IF
           END IF
           
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

         IF( PseudoComplex ) THEN
           beta_im = dotprodfun(n, T2(1:n), 1, r(1:n), 1 )

           IF( HUTI_PSEUDOCOMPLEX == 2 ) THEN
             x(1:n/2) = x(1:n/2) - beta_im * T1(1+n/2:n)
             x(1+n/2:n) = x(1+n/2:n) + beta_im * T1(1:n/2)                    
             r(1:n/2) = r(1:n/2) + beta_im * T2(1+n/2:n)
             r(1+n/2:n) = r(1+n/2:n) - beta_im * T2(1:n/2)                    
           ELSE
             x(1:n:2) = x(1:n:2) - beta_im * T1(2:n:2)
             x(2:n:2) = x(2:n:2) + beta_im * T1(1:n:2)                    
             r(1:n:2) = r(1:n:2) + beta_im * T2(2:n:2)
             r(2:n:2) = r(2:n:2) - beta_im * T2(1:n:2)                    
           END IF
         END IF
                      
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

         IF (UseStopCFun) THEN
           Residual = stopcfun(x,b,r,ipar,dpar)
           IF( MOD(k,OutputInterval) == 0) THEN
             WRITE (*, '(A, I6, 2E12.4)') '   gcr:',k, rnorm / bnorm, residual
             CALL FLUSH(6)
           END IF
         ELSE
           Residual = rnorm / bnorm
           IF( MOD(k,OutputInterval) == 0) THEN
             IF( PseudoComplex ) THEN
               WRITE (*, '(A, I6, 3E12.4, A)') '   gcr:',k, residual, beta, beta_im,'i'
             ELSE
               WRITE (*, '(A, I6, 2E12.4)') '   gcr:',k, residual, beta
             END IF
             CALL FLUSH(6)
           END IF
         END IF
           
         Converged = (Residual < MinTolerance) .AND. ( k >= MinIter )
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
             CALL Warn('IterMethod_GCR','Iterated GCR solution may not be accurate')
             i = 4
           ELSE
             i = 8
           END IF
           WRITE( Message,'(A,I0,A,ES12.3)') 'Iterated residual norm after ',k,' iters:', rnorm
           CALL Info('IterMethod_GCR', Message, Level=i)
           WRITE( Message,'(A,ES12.3)') 'True residual norm::', TrueResNOrm
           CALL Info('IterMethod_GCR', Message, Level=i)            

           IF( InfoActive(20) ) THEN
             ksum = ksum + k
             CALL Info('IterMethod_GCR','Total number of GCR iterations: '//I2S(ksum))           
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


   
!-----------------------------------------------------------------------------------
!>  This subroutine solves real linear systems Ax = b by using the IDR(s) algorithm
!>  with s >= 1 and the right-oriented preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_idrs( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
!------------------------------------------------------------------------------
    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_d ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_d ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_d ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_d ), POINTER :: normfun
    PROCEDURE( stopc_iface_d ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    REAL(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    ! DOUBLE PRECISION, DIMENSION(HUTI_NDIM), TARGET :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    REAL(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work
    ! DOUBLE PRECISION, DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work
    INTEGER :: ndim,i,j,k
    INTEGER :: Rounds, OutputInterval, s
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, UseStopCFun

    TYPE(Matrix_t), POINTER :: A

    REAL(KIND=dp), POINTER :: x(:),b(:)

    ! Variables related to robust mode
    LOGICAL :: Robust 
    INTEGER :: BestIter,BadIterCount,MaxBadIter
    REAL(KIND=dp) :: BestNorm,RobustStep,RobustTol,RobustMaxTol
    REAL(KIND=dp), ALLOCATABLE :: Bestx(:)

    LOGICAL :: Smoothing 

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
    s = HUTI_IDRS_S
    UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC
    
    Robust = ( HUTI_ROBUST == 1 )
    IF( Robust ) THEN
      RobustTol = HUTI_ROBUST_TOLERANCE
      RobustStep = HUTI_ROBUST_STEPSIZE
      RobustMaxTol = HUTI_ROBUST_MAXTOLERANCE
      MaxBadIter = HUTI_ROBUST_MAXBADIT
      BestNorm = SQRT(HUGE(BestNorm))
      BadIterCount = 0
      BestIter = 0      
      ALLOCATE( BestX(ndim))
    END IF

    Smoothing = ( HUTI_SMOOTHING == 1) 

    Converged = .FALSE.
    Diverged = .FALSE.
    
    CALL RealIDRS(ndim+nc, A,x,b, Rounds, MinTol, MaxTol, &
         Converged, Diverged, OutputInterval, s )

    IF(Constrained) THEN
      xvec=x(1:ndim)
      rhsvec=b(1:ndim)
      CM % extraVals = x(ndim+1:ndim+nc)
      DEALLOCATE(x,b)
    END IF

    IF( Robust ) THEN
      DEALLOCATE( BestX )
    END IF
    

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER

  CONTAINS

!-----------------------------------------------------------------------------------
!   The subroutine RealIDRS solves real linear systems Ax = b by using the IDR(s) 
!   algorithm with s >= 1 and the right-oriented preconditioning.
!
!   The subroutine RealIDRS has been written by M.B. van Gijzen
!----------------------------------------------------------------------------------- 
    SUBROUTINE RealIDRS( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
        Diverged, OutputInterval, s)
!----------------------------------------------------------------------------------- 
      INTEGER :: s   ! IDR parameter
      INTEGER :: n, MaxRounds, OutputInterval   
      LOGICAL :: Converged, Diverged
      TYPE(Matrix_t), POINTER :: A
      REAL(KIND=dp) :: x(n), b(n)
      REAL(KIND=dp) :: Tol, MaxTol

      ! Local arrays:
      REAL(kind=dp) :: P(n,s)
      REAL(kind=dp) :: G(n,s)
      REAL(kind=dp) :: U(n,s)
      REAL(kind=dp) :: r(n)
      REAL(kind=dp) :: v(n)
      REAL(kind=dp) :: t(n)
      REAL(kind=dp) :: M(s,s), f(s), mu(s)
      REAL(kind=dp) :: alpha(s), beta(s), gamma(s)

      REAL(kind=dp) :: om, tr, tr_s, tt
      REAL(kind=dp) :: nr, nt, rho, kappa

      REAL(kind=dp), ALLOCATABLE :: r_s(:), x_s(:)
      REAL(kind=dp) :: theta
      
      INTEGER :: iter                         ! number of iterations
      INTEGER :: ii                           ! inner iterations index
      INTEGER :: jj                           ! G-space index
      REAL(kind=dp) :: normb, normr, errorind ! for tolerance check
      INTEGER :: i,j,k,l                      ! loop counters
!----------------------------------------------------------------------------------- 
      
      ! Compute initial residual
      normb = normfun(n,b,1)
      CALL C_matvec( x, t, ipar, matvecsubr )
      r = b - t
      
      !-------------------------------------------------------------------
      ! Check whether the initial guess satisfies the stopping criterion
      !--------------------------------------------------------------------
      IF (UseStopCFun) THEN
        errorind = stopcfun(x,b,r,ipar,dpar)
      ELSE
        normr = normfun(n,r,1)
        errorind = normr / normb
      END IF
      Converged = (errorind < Tol)
      Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)

      IF ( Converged .OR. Diverged ) RETURN

      U = 0.0d0
      IF ( Smoothing ) THEN
        ALLOCATE( r_s(n), x_s(n) )
        x_s = x
        r_s = r 
      END IF

      ! Define P(n,s) and kappa
#if 1
      CALL RANDOM_NUMBER(P)
#else
      ! this is alternative generation of initial basis vectors
      ! it is deterministic but not as good...
      l = 0
      k = 2        
      DO j=1,s
        DO i=1,n
          P(i,j) = MODULO(i+l,k) / (1.0*(k-1)) 
        END DO
        l = k
        k = 2*k + 1
      END DO
#endif
              
      DO j = 1,s
        DO k = 1,j-1
          alpha(k) = dotprodfun(n, P(:,k), 1, P(:,j), 1 )
          P(:,j) = P(:,j) - alpha(k)*P(:,k)
        END DO
        P(:,j) = P(:,j)/normfun(n,P(:,j),1)
      END DO
      kappa = 0.7d0

      ! Initialize local variables:
      M = 0.0d0
      om = 1.0d0
      iter = 0
      jj = 0
      ii = 0

      
      ! This concludes the initialisation phase    
      
      ! Main iteration loop, build G-spaces:
      DO WHILE ( (.NOT. Converged) .AND. (.NOT. Diverged) ) 

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Generate s vectors in G_j
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! New right-hand side for small system:
        DO k = 1,s
          f(k) = dotprodfun(n, P(:,k), 1, r, 1 )
        END DO

        DO k = 1,s

          ! Update inner iteration counter
          ii = ii + 1

          ! Compute new v
          v = r
          IF ( jj > 0 ) THEN

            ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            DO i = k,s
              gamma(i) = f(i)
              DO j = k,i-1
                gamma(i) = gamma(i) - M(i,j)*gamma(j)
              END DO
              gamma(i) = gamma(i)/M(i,i)
              v = v - gamma(i)*G(:,i)
            END DO

            ! Compute new U(:,k)
            CALL C_rpcond( t, v, ipar, pcondrsubr ) 
            t = om*t
            DO i = k,s
              t = t + gamma(i)*U(:,i)
            END DO
            U(:,k) = t

          ELSE

            ! Updates for the first s iterations (in G_0):
            CALL C_rpcond( U(:,k), v, ipar, pcondrsubr )

          END IF

          ! Compute new G(:,k), G(:,k) is in space G_j
          CALL C_matvec( U(:,k), G(:,k), ipar, matvecsubr )

          ! Bi-Orthogonalise the new basis vectors:
          DO i = 1,s
            mu(i) = dotprodfun(n, P(:,i), 1, G(:,k), 1 )
          END DO
          DO i = 1,k-1
            alpha(i) = mu(i)
            DO j = 1, i-1
              alpha(i) = alpha(i) - M(i,j)*alpha(j)
            END DO
            alpha(i) = alpha(i)/M(i,i)
            G(:,k) = G(:,k) - G(:,i)*alpha(i)
            U(:,k) = U(:,k) - U(:,i)*alpha(i)
            mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
          END DO
          M(k:s,k) = mu(k:s)

          ! Break down?
          IF ( ABS(M(k,k)) <= TINY(tol) ) THEN
            Diverged = .TRUE.
            EXIT
          END IF

          ! Make r orthogonal to p_i, i = 1..k, update solution and residual
          beta(k) = f(k)/M(k,k)
          r = r - beta(k)*G(:,k)
          x = x + beta(k)*U(:,k)
          
          ! New f = P'*r (first k  components are zero)
          IF ( k < s ) THEN
            f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
          END IF

          IF (Smoothing) THEN
            t = r_s - r
            tr_s = dotprodfun(n, t, 1, r_s, 1 )
            tt = dotprodfun(n, t, 1, t, 1 )
            theta = tr_s / tt
            
            r_s = r_s - theta * t
            x_s = x_s - theta * (x_s - x)
          END IF

          ! Check for convergence
          iter = iter + 1
          IF (UseStopCFun) THEN
            IF (Smoothing) THEN
              errorind = stopcfun(x_s,b,r_s,ipar,dpar)
            ELSE
              errorind = stopcfun(x,b,r,ipar,dpar)
            END IF
          ELSE
            IF (Smoothing) THEN
              normr = normfun(n,r_s,1)  
            ELSE
              normr = normfun(n,r,1)
            END IF
            errorind = normr/normb
          END IF

          IF( MOD(iter,OutputInterval) == 0) THEN
            WRITE (*, '(I8, E11.4)') iter, errorind
          END IF

          Converged = (errorind < Tol)
          Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
          IF ( Converged .OR. Diverged ) EXIT
          IF (iter == MaxRounds) EXIT
          
         
        END DO ! Now we have computed s+1 vectors in G_j
        IF ( Converged .OR. Diverged ) EXIT
        IF (iter == MaxRounds) EXIT

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Compute first residual in G_j+1
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! Update G-space counter
        jj = jj + 1

        ! Compute first residual in G_j+1
        ! Note: r is already perpendicular to P so v = r

        ! Preconditioning:
        CALL C_rpcond( v, r, ipar, pcondrsubr )
        ! Matrix-vector multiplication:
        CALL C_matvec( v, t, ipar, matvecsubr )

        ! Computation of a new omega
        ! 'Maintaining the convergence':
        nr = normfun(n,r,1)
        nt = normfun(n,t,1)
        tr = dotprodfun(n, t, 1, r, 1 )
        rho = ABS(tr/(nt*nr))
        om=tr/(nt*nt)
        IF ( rho < kappa ) THEN
          om = om*kappa/rho
        END IF

        IF ( ABS(om) <= EPSILON(tol) ) THEN
          Diverged = .TRUE.
          EXIT
        END IF

        ! Update solution and residual
        r = r - om*t
        x = x + om*v
        
        IF (Smoothing) THEN
          t = r_s - r
          tr_s = dotprodfun(n, t, 1, r_s, 1 )
          tt = dotprodfun(n, t, 1, t, 1 )
          theta = tr_s / tt
          r_s = r_s - theta * t
          x_s = x_s - theta * (x_s - x)
        END IF

        ! Check for convergence
        iter = iter + 1
        IF (UseStopCFun) THEN
          IF (Smoothing) THEN
            errorind = stopcfun(x_s,b,r_s,ipar,dpar)
          ELSE
            errorind = stopcfun(x,b,r,ipar,dpar)
          END IF
        ELSE
          IF (Smoothing) THEN
            normr = normfun(n,r_s,1)  
          ELSE
            normr = normfun(n,r,1)
          END IF
          errorind = normr/normb
        END IF
        
        IF( MOD(iter,OutputInterval) == 0) THEN
          WRITE (*, '(I8, E11.4)') iter, errorind
          CALL FLUSH(6)
        END IF

        IF( Robust ) THEN
          ! Always store the best solution so far (with some small margin)
          IF( errorInd < RobustStep * BestNorm ) THEN
            BestIter = iter
            BestNorm = errorInd
            IF (Smoothing) THEN
              Bestx = x_s
            ELSE
              Bestx = x
            END IF
            BadIterCount = 0
          ELSE
            BadIterCount = BadIterCount + 1
          END IF

          ! If we have diverged too much and have found already a good candidate, then take it
          IF( BestNorm <  RobustTol .AND. &
              ( errorInd > RobustMaxTol .OR. BadIterCount > MaxBadIter ) ) THEN
            EXIT
          END IF
          
        END IF
                        
        Converged = (errorind < Tol)
        Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
        IF (iter == MaxRounds) EXIT
      END DO ! end of while loop

      IF( Smoothing ) x = x_s
      
      IF( Robust ) THEN
        IF( BestNorm < RobustTol ) THEN
          Converged = .TRUE.
        END IF
        IF( BestNorm < errorInd ) THEN
          x = Bestx
        END IF        
        IF(OutputInterval /= HUGE(OutputInterval)) THEN
          WRITE(*,'(A,I8,E11.4,I8,E11.4)') 'Idrs robust: ',&
              iter, BestNorm, BestIter, errorind
          CALL FLUSH(6)
        END IF
      ELSE
        IF(OutputInterval /= HUGE(OutputInterval)) THEN
          WRITE(*,'(A,I8,E11.4)') 'Idrs: ',&
              iter, errorind
          CALL FLUSH(6)
        END IF
      END IF
      
    !----------------------------------------------------------
    END SUBROUTINE RealIDRS
    !----------------------------------------------------------

!--------------------------------------------------------------
  END SUBROUTINE itermethod_idrs
!--------------------------------------------------------------


!------------------------------------------------------------------------------
!> This routine provides the complex version to the GCR linear solver.
!------------------------------------------------------------------------------
 SUBROUTINE itermethod_z_gcr( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
!------------------------------------------------------------------------------
    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    COMPLEX(KIND=dp) :: y(HUTI_NDIM),f(HUTI_NDIM)
    INTEGER :: ndim, RestartN, i
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol, Residual
    LOGICAL :: Converged, Diverged, UseStopCFun

    CALL Info('Itermetod_z_gcr','Starting GCR iteration',Level=25)
    
    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    RestartN = HUTI_GCR_RESTART 
    UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC
    
    Converged = .FALSE.
    Diverged = .FALSE.
    
    !----------------------------------------------------------------------------
    ! Transform the solution vector and the right-hand side vector to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------
    DO i=1,ndim
      y(i)=xvec(i)
      f(i)=rhsvec(i)
    END DO
       
    CALL GCR_Z(ndim, GlobalMatrix, y, f, Rounds, MinTol, MaxTol, Residual, &
        Converged, Diverged, OutputInterval, RestartN )

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER
   
    DO i=1,ndim
      xvec(i) = y(i)
    END DO

  CONTAINS 
    
    
!------------------------------------------------------------------------------  
    SUBROUTINE GCR_Z( n, A, x, b, Rounds, MinTolerance, MaxTolerance, Residual, &
        Converged, Diverged, OutputInterval, m) 
!------------------------------------------------------------------------------
      TYPE(Matrix_t), POINTER :: A
      INTEGER :: Rounds
      COMPLEX(KIND=dp) :: x(n),b(n)
      LOGICAL :: Converged, Diverged
      REAL(KIND=dp) :: MinTolerance, MaxTolerance, Residual
      INTEGER :: n, OutputInterval, m
!------------------------------------------------------------------------------      
      REAL(KIND=dp) :: bnorm,rnorm
      COMPLEX(KIND=dp), ALLOCATABLE :: R(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: S(:,:), V(:,:), T1(:), T2(:)
      INTEGER :: i,j,k,allocstat
      COMPLEX(KIND=dp) :: beta, czero
      REAL(KIND=dp) :: alpha, trueresnorm, normerr
      COMPLEX(KIND=dp), ALLOCATABLE :: trueres(:)
!------------------------------------------------------------------------------
            
      ALLOCATE( R(n), T1(n), T2(n), trueres(n), STAT=allocstat )
      IF( allocstat /= 0) &
          CALL Fatal('GCR_Z','Failed to allocate memory of size: '//I2S(n))
      IF ( m > 1 ) THEN        
         ALLOCATE( S(n,m-1), V(n,m-1), STAT=allocstat )
         IF ( allocstat /= 0 ) THEN
           CALL Fatal('GCR_Z','Failed to allocate memory of size: '&
               //I2S(n)//' x '//I2S(m-1))
         END IF
         czero = CMPLX( 0.0_dp, 0.0_dp )
         V(1:n,1:m-1) = czero
         S(1:n,1:m-1) = czero
      END IF	
      
      CALL matvecsubr( x, r, ipar )
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
         CALL pcondrsubr( T1, r, ipar )         
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

         IF (UseStopCFun) THEN
           Residual = stopcfun(x,b,r,ipar,dpar)
           IF( MOD(k,OutputInterval) == 0) THEN
             WRITE (*, '(A, I6, 2E12.4)') '   gcr:',k, rnorm / bnorm, residual
             CALL FLUSH(6)
           END IF
         ELSE
           Residual = rnorm / bnorm
           IF( MOD(k,OutputInterval) == 0) THEN
             WRITE (*, '(A, I8, 3ES12.4,A)') '   gcrz:',k, residual, beta,'i'
             CALL FLUSH(6)
           END IF
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
               CALL FLUSH(6)
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
    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

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
      y(i)=xvec(i)
      f(i)=rhsvec(i)
    END DO

    CALL ComplexBiCGStabl(ndim, GlobalMatrix, y, f, Rounds, MinTol, MaxTol, &
         Converged, Diverged, OutputInterval, PolynomialDegree)

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER

    DO i=1,ndim
      xvec(i) = y(i)
    END DO

  CONTAINS 

    !-----------------------------------------------------------------------------------
    SUBROUTINE ComplexBiCGStabl( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
         Diverged, OutputInterval, l)
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
               CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error 1.' )
            ENDIF
            beta = alpha*(rho1/rho0)
            rho0 = rho1
            DO j=0,k-1
               work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
            ENDDO
            CALL pcondrsubr( t, work(1:n,u+k-1), ipar )
            CALL matvecsubr( t, work(1:n,u+k),   ipar )

            sigma = dotprodfun(n, work(1:n,rr), 1, work(1:n,u+k), 1)
            IF (sigma == zzero) THEN
               CALL Fatal( 'ComplexBiCGStab(l)', 'Breakdown error 2.' )
            ENDIF
            alpha = rho1/sigma
            x(1:n) = x(1:n) + alpha * work(1:n,u)
            DO j=0,k-1
               work(1:n,r+j) = work(1:n,r+j) - alpha * work(1:n,u+j+1)
            ENDDO
            CALL pcondrsubr( t, work(1:n,r+k-1), ipar )
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
         call zmv( rwork(1:l+1,z:z+l), rwork(1:l+1,y0), rwork(1:l+1,y), l+1 )
         kappa0 = SQRT( ABS(zdotc(l+1, rwork(1:l+1,y0), 1, rwork(1:l+1,y), 1)) )  ! replace zdotc

         call zmv( rwork(1:l+1,z:z+l), rwork(1:l+1,yl), rwork(1:l+1,y), l+1 )
         kappal = SQRT( ABS(zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1)) )  ! replace zdotc

         call zmv( rwork(1:l+1,z:z+l), rwork(1:l+1,y0), rwork(1:l+1,y), l+1 )
         varrho = zdotc(l+1, rwork(1:l+1,yl), 1, rwork(1:l+1,y), 1) / &           ! replace zdotc
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

         call zmv( rwork(1:l+1,z:z+l), rwork(1:l+1,y0), rwork(1:l+1,y), l+1 )
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
            CALL pcondrsubr( t, x, ipar )
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
            CALL pcondrsubr( t, x, ipar )
            t(1:n) = t(1:n)+work(1:n,xp)
         END IF

         errorind = rnrm/bnrm
         IF( MOD(Round,OutputInterval) == 0) THEN
           WRITE (*, '(I8, E11.4)') Round, errorind
           CALL FLUSH(6)
         END IF

         Converged = (errorind < Tol) 
         Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)    
         IF( Converged .OR. Diverged) EXIT    
      END DO

      IF( EarlyExit .AND. (OutputInterval/=HUGE(OutputInterval)) ) THEN
        WRITE (*, '(I8, E11.4)') Round, errorind
        CALL FLUSH(6)
      END IF

      !------------------------------------------------------------
      ! We have solved z = P*x, so finally solve the true unknown x
      !------------------------------------------------------------
      t(1:n) = x(1:n)
      CALL pcondrsubr( x, t, ipar )
      x(1:n) = x(1:n) + work(1:n,xp)      

    !----------------------------------------------------------
    END SUBROUTINE ComplexBiCGStabl
    !----------------------------------------------------------

!--------------------------------------------------------------

    SUBROUTINE zmv(a,x,y,n)
      integer, INTENT(in) :: n
      complex(kind=dp), INTENT(in)  ::a(n,n), x(n)
      complex(kind=dp), INTENT(out) ::y(n)
 
      complex(kind=dp), parameter :: zone  = cmplx(1._dp, 0._dp,kind=dp)
      complex(kind=dp), parameter :: zzero = cmplx(0._dp, 0._dp,kind=dp)
    
      IF(n>8) THEN
        CALL zhemv ('u', n, zone, a, n, x, 1, zzero, y, 1)
        RETURN
      END IF

      y(1) = a(1,1)*x(1) + a(1,2)*x(2) + a(1,3)*x(3)
      y(2) = a(2,1)*x(1) + a(2,2)*x(2) + a(2,3)*x(3)
      y(3) = a(3,1)*x(1) + a(3,2)*x(2) + a(3,3)*x(3)
      IF(n==3) RETURN

      y(1) = y(1) + a(1,4)*x(4)
      y(2) = y(2) + a(2,4)*x(4)
      y(3) = y(3) + a(3,4)*x(4)
      y(4) = a(4,1)*x(1)+a(4,2)*x(2)+a(4,3)*x(3)+a(4,4)*x(4)
      IF(n==4) RETURN

      y(1) = y(1) + a(1,5)*x(5)
      y(2) = y(2) + a(2,5)*x(5)
      y(3) = y(3) + a(3,5)*x(5)
      y(4) = y(4) + a(4,5)*x(5)
      y(5) = a(5,1)*x(1)+a(5,2)*x(2)+a(5,3)*x(3)+a(5,4)*x(4)+ &
             a(5,5)*x(5)
      IF(n==5) RETURN

      y(1) = y(1) + a(1,6)*x(6)
      y(2) = y(2) + a(2,6)*x(6)
      y(3) = y(3) + a(3,6)*x(6)
      y(4) = y(4) + a(4,6)*x(6)
      y(5) = y(5) + a(5,6)*x(6)
      y(6) = a(6,1)*x(1)+a(6,2)*x(2)+a(6,3)*x(3)+a(6,4)*x(4)+ &
             a(6,5)*x(5)+a(6,6)*x(6)
      IF(n==6) RETURN

      y(1) = y(1) + a(1,7)*x(7)
      y(2) = y(2) + a(2,7)*x(7)
      y(3) = y(3) + a(3,7)*x(7)
      y(4) = y(4) + a(4,7)*x(7)
      y(5) = y(5) + a(5,7)*x(7)
      y(6) = y(6) + a(6,7)*x(7)
      y(7) = a(7,1)*x(1)+a(7,2)*x(2)+a(7,3)*x(3)+a(7,4)*x(4)+ &
             a(7,5)*x(5)+a(7,6)*x(6)+a(7,7)*x(7)
      IF(n==7) RETURN

      y(1) = y(1) + a(1,8)*x(8)
      y(2) = y(2) + a(2,8)*x(8)
      y(3) = y(3) + a(3,8)*x(8)
      y(4) = y(4) + a(4,8)*x(8)
      y(5) = y(5) + a(5,8)*x(8)
      y(6) = y(6) + a(6,8)*x(8)
      y(7) = y(7) + a(7,8)*x(8)
      y(8) = a(8,1)*x(1)+a(8,2)*x(2)+a(8,3)*x(3)+a(8,4)*x(4)+ &
             a(8,5)*x(5)+a(8,6)*x(6)+a(8,7)*x(7)+a(8,8)*x(8)

    END SUBROUTINE zmv

!--------------------------------------------------------------
  END SUBROUTINE itermethod_z_bicgstabl
!--------------------------------------------------------------


!------------------------------------------------------------------------------
!> This routine provides a complex version of the IDR(s) solver for linear systems.
!------------------------------------------------------------------------------
  SUBROUTINE itermethod_z_idrs( xvec, rhsvec, &
      ipar, dpar, work, matvecsubr, pcondlsubr, &
      pcondrsubr, dotprodfun, normfun, stopcfun )
!------------------------------------------------------------------------------
    USE huti_interfaces
    IMPLICIT NONE
    PROCEDURE( mv_iface_z ), POINTER :: matvecsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondlsubr
    PROCEDURE( pc_iface_z ), POINTER :: pcondrsubr
    PROCEDURE( dotp_iface_z ), POINTER :: dotprodfun
    PROCEDURE( norm_iface_z ), POINTER :: normfun
    PROCEDURE( stopc_iface_z ), POINTER :: stopcfun

    INTEGER, DIMENSION(HUTI_IPAR_DFLTSIZE) :: ipar
    COMPLEX(KIND=dp), TARGET, DIMENSION(HUTI_NDIM) :: xvec, rhsvec
    DOUBLE PRECISION, DIMENSION(HUTI_DPAR_DFLTSIZE) :: dpar
    COMPLEX(KIND=dp), DIMENSION(HUTI_WRKDIM,HUTI_NDIM) :: work

    COMPLEX(KIND=dp) :: y(HUTI_NDIM),f(HUTI_NDIM)
    INTEGER :: ndim, i, s
    INTEGER :: Rounds, OutputInterval
    REAL(KIND=dp) :: MinTol, MaxTol
    LOGICAL :: Converged, Diverged
    !--------------------------------------------------------------------------------    

    ndim = HUTI_NDIM
    Rounds = HUTI_MAXIT
    MinTol = HUTI_TOLERANCE
    MaxTol = HUTI_MAXTOLERANCE
    OutputInterval = HUTI_DBUGLVL
    s = HUTI_IDRS_S 
    !----------------------------------------------------------------------------
    ! Transform the solution vector and the right-hand side vector to 
    ! complex-valued vectors y and f
    !---------------------------------------------------------------------------
    DO i=1,ndim
      y(i)=xvec(i)
      f(i)=rhsvec(i)
    END DO
    CALL ComplexIDRS(ndim, GlobalMatrix, y, f, Rounds, MinTol, MaxTol, &
         Converged, Diverged, OutputInterval, s )

    IF(Converged) HUTI_INFO = HUTI_CONVERGENCE
    IF(Diverged) HUTI_INFO = HUTI_DIVERGENCE
    IF ( (.NOT. Converged) .AND. (.NOT. Diverged) ) HUTI_INFO = HUTI_MAXITER

    DO i=1,ndim
      xvec(i) = y(i)
    END DO

  CONTAINS 

!----------------------------------------------------------------------------------- 
!   This subroutine solves complex linear systems Ax = b by using the IDR(s) algorithm 
!   with s >= 1 and the right-oriented preconditioning. 
!
!   The subroutine ComplexIDRS has been written by M.B. van Gijzen
!-----------------------------------------------------------------------------------
    SUBROUTINE ComplexIDRS( n, A, x, b, MaxRounds, Tol, MaxTol, Converged, &
        Diverged, OutputInterval, s )
!-----------------------------------------------------------------------------------
      INTEGER :: s  
      INTEGER :: n, MaxRounds, OutputInterval   
      LOGICAL :: Converged, Diverged, UseStopCFun
      TYPE(Matrix_t), POINTER :: A
      COMPLEX(KIND=dp) :: x(n), b(n)
      REAL(KIND=dp) :: Tol, MaxTol
!------------------------------------------------------------------------------

      ! Local arrays:
      REAL(kind=dp) :: Pr(n,s), Pi(n,s) 
      COMPLEX(kind=dp) :: P(n,s)
      COMPLEX(kind=dp) :: G(n,s)
      COMPLEX(kind=dp) :: U(n,s)
      COMPLEX(kind=dp) :: r(n) 
      COMPLEX(kind=dp) :: v(n)   
      COMPLEX(kind=dp) :: t(n)  
      COMPLEX(kind=dp) :: M(s,s), f(s), mu(s)
      COMPLEX(kind=dp) :: alpha(s), beta(s), gamma(s)

      COMPLEX(kind=dp) :: om, tr    
      REAL(kind=dp) :: nr, nt, rho, kappa

      INTEGER :: iter                         ! number of iterations
      INTEGER :: ii                           ! inner iterations index
      INTEGER :: jj                           ! G-space index
      REAL(kind=dp) :: normb, normr, errorind ! for tolerance check
      INTEGER :: i,j,k,l                      ! loop counters

      UseStopCFun = HUTI_STOPC == HUTI_USUPPLIED_STOPC
      
      U = 0.0d0

      ! Compute initial residual, set absolute tolerance
      normb = normfun(n,b,1)
      CALL matvecsubr( x, t, ipar )
      r = b - t
      IF (UseStopCFun) THEN
        errorind = stopcfun(x,b,r,ipar,dpar)
      ELSE
        normr = normfun(n,r,1)
        errorind = normr / normb
      END IF
      
      !-------------------------------------------------------------------
      ! Check whether the initial guess satisfies the stopping criterion
      !--------------------------------------------------------------------
      Converged = (errorind < Tol)
      Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)

      IF ( Converged .OR. Diverged ) RETURN

      ! Define P and kappa 
      CALL RANDOM_NUMBER(Pr)
      CALL RANDOM_NUMBER(Pi)
      P = Pr + (0.,1.)*Pi

      DO j = 1,s
        DO k = 1,j-1
          alpha(k) = dotprodfun(n, P(:,k), 1, P(:,j), 1 )
          P(:,j) = P(:,j) - alpha(k)*P(:,k)
        END DO
        P(:,j) = P(:,j)/normfun(n,P(:,j),1)
      END DO
      kappa = 0.7d0

      ! Initialize local variables:
      M = 0.0d0
      om = 1.0d0
      iter = 0
      jj = 0
      ii = 0
      ! This concludes the initialisation phase

      ! Main iteration loop, build G-spaces:

      DO WHILE ( .NOT. Converged .AND. .NOT. Diverged )  ! start of iteration loop

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Generate s vectors in G_j
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! New right-hand side for small system:
        DO k = 1,s
          f(k) = dotprodfun(n, P(:,k), 1, r, 1 )
        END DO

        DO k = 1,s

          ! Update inner iteration counter
          ii = ii + 1

          ! Compute new v
          v = r 
          IF ( jj > 0 ) THEN

            ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
            DO i = k,s
              gamma(i) = f(i)
              DO j = k,i-1
                gamma(i) = gamma(i) - M(i,j)*gamma(j)
              END DO
              gamma(i) = gamma(i)/M(i,i)
              v = v - gamma(i)*G(:,i)
            END DO

            ! Compute new U(:,k)
            CALL pcondrsubr( t, v, ipar )
            t = om*t
            DO i = k,s
              t = t + gamma(i)*U(:,i)
            END DO
            U(:,k) = t

          ELSE 

            ! Updates for the first s iterations (in G_0):
            CALL pcondrsubr( U(:,k), v, ipar )

          END IF

          ! Compute new G(:,k), G(:,k) is in space G_j
          CALL matvecsubr( U(:,k), G(:,k), ipar )

          ! Bi-Orthogonalise the new basis vectors: 
          DO i = 1,s
            mu(i) = dotprodfun(n, P(:,i), 1, G(:,k), 1 )
          END DO
          DO i = 1,k-1
            alpha(i) = mu(i)
            DO j = 1, i-1
              alpha(i) = alpha(i) - M(i,j)*alpha(j)
            END DO
            alpha(i) = alpha(i)/M(i,i)
            G(:,k) = G(:,k) - G(:,i)*alpha(i)
            U(:,k) = U(:,k) - U(:,i)*alpha(i)
            mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
          END DO
          M(k:s,k) = mu(k:s)

          ! Break down?
          IF ( ABS(M(k,k)) <= TINY(tol) ) THEN
            Diverged = .TRUE.
            EXIT
          END IF

          ! Make r orthogonal to p_i, i = 1..k, update solution and residual 
          beta(k) = f(k)/M(k,k)
          r = r - beta(k)*G(:,k)
          x = x + beta(k)*U(:,k)

          ! New f = P'*r (first k  components are zero)
          IF ( k < s ) THEN
            f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
          END IF

          ! Check for convergence
          IF (UseStopCFun) THEN
            errorind = stopcfun(x,b,r,ipar,dpar)
          ELSE
            normr = normfun(n,r,1)
            errorind = normr/normb
          END IF
          iter = iter + 1
          
          IF( MOD(iter,OutputInterval) == 0) THEN
            WRITE (*, '(I8, E11.4)') iter, errorind
            CALL FLUSH(6)
          END IF

          Converged = (errorind < Tol)
          Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
          IF ( Converged .OR. Diverged ) EXIT
          IF (iter == MaxRounds) EXIT

        END DO ! Now we have computed s+1 vectors in G_j
        IF ( Converged .OR. Diverged ) EXIT
        IF (iter == MaxRounds) EXIT

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Compute first residual in G_j+1
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! Update G-space counter
        jj = jj + 1

        ! Compute first residual in G_j+1
        ! Note: r is already perpendicular to P so v = r

        ! Preconditioning:
        CALL pcondrsubr( v, r, ipar )
        ! Matrix-vector multiplication:
        CALL matvecsubr( v, t, ipar )

        ! Computation of a new omega
        ! 'Maintaining the convergence':
        nr = normfun(n,r,1)
        nt = normfun(n,t,1)
        tr = dotprodfun(n, t, 1, r, 1 )
        rho = ABS(tr/(nt*nr))
        om=tr/(nt*nt)
        IF ( rho < kappa ) THEN
          om = om*kappa/rho
        END IF

        IF ( ABS(om) <= EPSILON(tol) ) THEN 
          Diverged = .TRUE.
          EXIT
        END IF

        ! Update solution and residual
        r = r - om*t 
        x = x + om*v 

        ! Check for convergence
        IF (UseStopCFun) THEN
          errorind = stopcfun(x,b,r,ipar,dpar)
        ELSE
          normr = normfun(n,r,1)
          errorind = normr/normb
        END IF
        iter = iter + 1

        IF( MOD(iter,OutputInterval) == 0) THEN
          WRITE (*, '(I8, E11.4)') iter, errorind
          CALL FLUSH(6)
        END IF

        Converged = (errorind < Tol)
        Diverged = (errorind > MaxTol) .OR. (errorind /= errorind)
        IF (iter == MaxRounds) EXIT

      END DO ! end of while loop

    !----------------------------------------------------------
    END SUBROUTINE ComplexIDRS
    !----------------------------------------------------------

!--------------------------------------------------------------
  END SUBROUTINE itermethod_z_idrs
!--------------------------------------------------------------

END MODULE IterativeMethods

!> \} ElmerLib
