!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Module containing the direct solvers for linear systems given in CRS format.
!> Included are Lapack band matrix solver, multifrontal Umfpack, MUMPS, SuperLU, 
!> and Pardiso. Note that many of these are linked in with ElmerSolver only 
!> if they are made available at the compilation time. 
!------------------------------------------------------------------------------

MODULE DirectSolve

   USE CRSMatrix
   USE Lists
   USE BandMatrix
   USE SParIterSolve
   USE SparIterGlobals

   IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Solver the complex linear system using direct band matrix solver from of Lapack.
!------------------------------------------------------------------------------
   SUBROUTINE ComplexBandSolver( A,x,b, Free_fact )
!------------------------------------------------------------------------------

     LOGICAL, OPTIONAL :: Free_Fact
     TYPE(Matrix_t) :: A
     REAL(KIND=dp) :: x(*),b(*)
!------------------------------------------------------------------------------

   
     INTEGER :: i,j,k,istat,Subband,N
     COMPLEX(KIND=dp), ALLOCATABLE :: BA(:,:)

     REAL(KIND=dp), POINTER CONTIG :: Values(:)
     INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

     SAVE BA
!------------------------------------------------------------------------------

     IF ( PRESENT(Free_Fact) ) THEN
       IF ( Free_Fact ) THEN
         IF ( ALLOCATED(BA) ) DEALLOCATE(BA)
         RETURN
       END IF
     END IF

     Rows => A % Rows
     Cols => A % Cols
     Diag => A % Diag
     Values => A % Values

     n = A % NumberOfRows
     x(1:n) = b(1:n)
     n = n / 2

     IF ( A % Format == MATRIX_CRS .AND. .NOT. A % Symmetric ) THEN
       Subband = 0
       DO i=1,N
         DO j=Rows(2*i-1),Rows(2*i)-1,2
           Subband = MAX(Subband,ABS((Cols(j)+1)/2-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(3*SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= 3*Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(3*SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(2*i-1),Rows(2*i)-1,2
           k = i - (Cols(j)+1)/2 + 2*Subband + 1
           BA(k,(Cols(j)+1)/2) = CMPLX(Values(j), -Values(j+1), KIND=dp )
         END DO
       END DO

       CALL SolveComplexBandLapack( N,1,BA,x,Subband,3*Subband+1 )

     ELSE IF ( A % Format == MATRIX_CRS ) THEN

       Subband = 0
       DO i=1,N
         DO j=Rows(2*i-1),Diag(2*i-1)
           Subband = MAX(Subband,ABS((Cols(j)+1)/2-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Direct solver memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(2*i-1),Diag(2*i-1)
           k = i - (Cols(j)+1)/2 + 1
           BA(k,(Cols(j)+1)/2) = CMPLX(Values(j), -Values(j+1), KIND=dp )
         END DO
       END DO

       CALL SolveComplexSBandLapack( N,1,BA,x,Subband,Subband+1 )

     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ComplexBandSolver 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solver the real linear system using direct band matrix solver from of Lapack.
!------------------------------------------------------------------------------
   SUBROUTINE BandSolver( A,x,b,Free_Fact )
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Free_Fact
     TYPE(Matrix_t) :: A
     REAL(KIND=dp) :: x(*),b(*)
!------------------------------------------------------------------------------

     INTEGER :: i,j,k,istat,Subband,N
     REAL(KIND=dp), ALLOCATABLE :: BA(:,:)

     REAL(KIND=dp), POINTER CONTIG :: Values(:)
     INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

     SAVE BA
!------------------------------------------------------------------------------
     IF ( PRESENT(Free_Fact) ) THEN
       IF ( Free_Fact ) THEN
         IF ( ALLOCATED(BA) ) DEALLOCATE(BA)
         RETURN
       END IF
     END IF

     N = A % NumberOfRows

     x(1:n) = b(1:n)

     Rows => A % Rows
     Cols => A % Cols
     Diag => A % Diag
     Values => A % Values

     IF ( A % Format == MATRIX_CRS ) THEN ! .AND. .NOT. A % Symmetric ) THEN
        Subband = 0
        DO i=1,N
          DO j=Rows(i),Rows(i+1)-1
            Subband = MAX(Subband,ABS(Cols(j)-i))
          END DO
        END DO

        IF ( .NOT.ALLOCATED( BA ) ) THEN

          ALLOCATE( BA(3*SubBand+1,N),stat=istat )

          IF ( istat /= 0 ) THEN
            CALL Fatal( 'BandSolver', 'Memory allocation error.' )
          END IF

        ELSE IF ( SIZE(BA,1) /= 3*Subband+1 .OR. SIZE(BA,2) /= N ) THEN

          DEALLOCATE( BA )
          ALLOCATE( BA(3*SubBand+1,N),stat=istat )

          IF ( istat /= 0 ) THEN
            CALL Fatal( 'BandSolver', 'Memory allocation error.' )
          END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(i),Rows(i+1)-1
           k = i - Cols(j) + 2*Subband + 1
           BA(k,Cols(j)) = Values(j)
         END DO
       END DO

       CALL SolveBandLapack( N,1,BA,x,Subband,3*Subband+1 )

     ELSE IF ( A % Format == MATRIX_CRS ) THEN

       Subband = 0
       DO i=1,N
         DO j=Rows(i),Diag(i)
           Subband = MAX(Subband,ABS(Cols(j)-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'BandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'BandSolver', 'Memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(i),Diag(i)
           k = i - Cols(j) + 1
           BA(k,Cols(j)) = Values(j)
         END DO
       END DO

       CALL SolveSBandLapack( N,1,BA,x,Subband,Subband+1 )

     ELSE IF ( A % Format == MATRIX_BAND ) THEN
       CALL SolveBandLapack( N,1,Values,x,Subband,3*Subband+1 )
     ELSE IF ( A % Format == MATRIX_SBAND ) THEN
       CALL SolveSBandLapack( N,1,Values,x,Subband,Subband+1 )
     END IF

!------------------------------------------------------------------------------
  END SUBROUTINE BandSolver 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using Umfpack multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE UMFPack_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
    LOGICAL, OPTIONAL :: Free_Fact
    TYPE(Matrix_t) :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), TARGET :: x(*), b(*)

    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

#include "../config.h"
#ifdef HAVE_UMFPACK
#ifdef USE_ISO_C_BINDINGS

#ifdef ARCH_32_BITS
#define CAddrInt c_int32_t
#else
#define CAddrInt c_int64_t
#endif
  ! Standard int version
  INTERFACE
    SUBROUTINE umf4def( control ) &
       BIND(C,name='umf4def')
       USE, INTRINSIC :: ISO_C_BINDING
       REAL(C_DOUBLE) :: control(*)
    END SUBROUTINE umf4def

    SUBROUTINE umf4sym( m,n,rows,cols,values,symbolic,control,iinfo ) &
       BIND(C,name='umf4sym')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: m,n,rows(*),cols(*)
       INTEGER(CAddrInt) ::  symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4sym

    SUBROUTINE umf4num( rows,cols,values,symbolic,numeric, control,iinfo ) &
       BIND(C,name='umf4num')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: rows(*),cols(*)
       INTEGER(CAddrInt) ::  numeric, symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4num

    SUBROUTINE umf4sol( sys, x, b, numeric, control, iinfo ) &
       BIND(C,name='umf4sol')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: sys
       INTEGER(CAddrInt) :: numeric
       REAL(C_DOUBLE) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4sol

    SUBROUTINE umf4fsym(symbolic) &
        BIND(C,name='umf4fsym')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: symbolic
    END SUBROUTINE umf4fsym

    SUBROUTINE umf4fnum(numeric) &
        BIND(C,name='umf4fnum')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: numeric
    END SUBROUTINE umf4fnum

  END INTERFACE

  ! Long int version
  INTERFACE
    SUBROUTINE umf4_l_def( control ) &
       BIND(C,name='umf4_l_def')
       USE, INTRINSIC :: ISO_C_BINDING
       REAL(C_DOUBLE) :: control(*)
    END SUBROUTINE umf4_l_def

    SUBROUTINE umf4_l_sym( m,n,rows,cols,values,symbolic,control,iinfo ) &
       BIND(C,name='umf4_l_sym')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) ::m,n,rows(*),cols(*) 
       INTEGER(C_LONG) :: m,n,rows(*),cols(*) !TODO: m,n of are called with AddrInt kind
       INTEGER(CAddrInt) ::  symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_sym

    SUBROUTINE umf4_l_num( rows,cols,values,symbolic,numeric, control,iinfo ) &
       BIND(C,name='umf4_l_num')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) :: rows(*),cols(*)
       INTEGER(C_LONG) :: rows(*),cols(*)
       INTEGER(CAddrInt) ::  numeric, symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_num

    SUBROUTINE umf4_l_sol( sys, x, b, numeric, control, iinfo ) &
       BIND(C,name='umf4_l_sol')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) :: sys
       INTEGER(C_LONG) :: sys
       INTEGER(CAddrInt) :: numeric
       REAL(C_DOUBLE) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4_l_sol

    SUBROUTINE umf4_l_fnum(numeric) &
        BIND(C,name='umf4_l_fnum')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: numeric
    END SUBROUTINE umf4_l_fnum

    SUBROUTINE umf4_l_fsym(symbolic) &
        BIND(C,name='umf4_l_fsym')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: symbolic
    END SUBROUTINE umf4_l_fsym
  END INTERFACE
#else
    INTERFACE
    SUBROUTINE umf4def( control )
       USE Types
       REAL(KIND=dp) :: control(*) 
    END SUBROUTINE umf4def

    SUBROUTINE umf4sym( m,n,rows,cols,values,symbolic,control,iinfo )
       USE Types
       INTEGER :: m,n,rows(*),cols(*)
       INTEGER(KIND=AddrInt) ::  symbolic
       REAL(KIND=dp) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4sym

    SUBROUTINE umf4num( rows,cols,values,symbolic,numeric, control,iinfo )
       USE Types
       INTEGER :: rows(*),cols(*)
       INTEGER(KIND=AddrInt) ::  numeric, symbolic
       REAL(KIND=dp) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4num

    SUBROUTINE umf4sol( sys, x, b, numeric, control, iinfo )
       USE Types
       INTEGER :: sys
       INTEGER(KIND=AddrInt) :: numeric
       REAL(KIND=dp) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4sol

    SUBROUTINE umf4_l_def( control )
       USE Types
       REAL(KIND=dp) :: control(*) 
    END SUBROUTINE umf4_l_def

    SUBROUTINE umf4_l_sym( m,n,rows,cols,values,symbolic,control,iinfo )
       USE Types
       INTEGER(KIND=AddrInt) :: m,n,rows(*),cols(*)
       INTEGER(KIND=AddrInt) ::  symbolic
       REAL(KIND=dp) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_sym

    SUBROUTINE umf4_l_num( rows,cols,values,symbolic,numeric, control,iinfo )
       USE Types
       INTEGER(KIND=AddrInt) :: rows(*),cols(*)
       INTEGER(KIND=AddrInt) ::  numeric, symbolic
       REAL(KIND=dp) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_num

    SUBROUTINE umf4_l_sol( sys, x, b, numeric, control, iinfo )
       USE Types
       INTEGER(KIND=AddrInt) :: sys
       INTEGER(KIND=AddrInt) :: numeric
       REAL(KIND=dp) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4_l_sol
  END INTERFACE
#endif

  INTEGER :: i, n, status, sys
  REAL(KIND=dp) :: iInfo(90), Control(20)
  INTEGER(KIND=AddrInt) :: symbolic, zero=0
  INTEGER(KIND=C_LONG) :: ln, lsys
  INTEGER(KIND=C_LONG), ALLOCATABLE :: LRows(:), LCols(:)

  SAVE iInfo, Control
 
  LOGICAL :: Factorize, FreeFactorize, stat, BigMode

  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % UMFPack_Numeric/=0 ) THEN
        CALL umf4fnum(A % UMFPack_Numeric)
        A % UMFPack_Numeric = 0
      END IF
      RETURN
    END IF
  END IF

  BigMode = ListGetString( Solver % Values, &
      'Linear System Direct Method' ) == 'big umfpack'

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  n = A % NumberofRows
  Rows => A % Rows
  Cols => A % Cols
  Diag => A % Diag
  Values => A % Values

  IF ( Factorize .OR. A% UmfPack_Numeric==0 ) THEN
    IF ( A % UMFPack_Numeric /= 0 ) THEN
      IF( BigMode ) THEN
        CALL umf4_l_fnum( A % UMFPack_Numeric )
      ELSE
        CALL umf4fnum( A % UMFPack_Numeric )
      END IF
      A % UMFPack_Numeric = 0
    END IF

    IF ( BigMode ) THEN
      ALLOCATE( LRows(SIZE(Rows)), LCols(SIZE(Cols)) )
      DO i=1,n
        LRows(i) = Rows(i)-1
        LCols(i) = Cols(i)-1
      END DO
      ln = n ! TODO: Kludge: ln is AddrInt and n is regular INTEGER
      CALL umf4_l_def( Control )
      CALL umf4_l_sym( ln,ln, LRows, LCols, Values, Symbolic, Control, iInfo )
    ELSE
      Rows = Rows-1
      Cols = Cols-1
      CALL umf4def( Control )
      CALL umf4sym( n,n, Rows, Cols, Values, Symbolic, Control, iInfo )
    END IF

    IF (iinfo(1)<0) THEN
      PRINT *, 'Error occurred in umf4sym: ', iinfo(1)
      STOP
    END IF

    IF ( BigMode ) THEN
      CALL umf4_l_num(LRows, LCols, Values, Symbolic, A % UMFPack_Numeric, Control, iInfo )
    ELSE
      CALL umf4num( Rows, Cols, Values, Symbolic, A % UMFPack_Numeric, Control, iInfo )
    END IF

    IF (iinfo(1)<0) THEN
      PRINT*, 'Error occurred in umf4num: ', iinfo(1)
      STOP
    ENDIF

    IF ( BigMode ) THEN
      DEALLOCATE( LRows, LCols )
      CALL umf4_l_fsym( Symbolic )
    ELSE
      A % Rows = A % Rows+1
      A % Cols = A % Cols+1
      CALL umf4fsym( Symbolic )
    END IF
  END IF

  IF ( BigMode ) THEN
    lsys = 2
    CALL umf4_l_sol( lsys, x, b, A % UMFPack_Numeric, Control, iInfo )
  ELSE
    sys = 2
    CALL umf4sol( sys, x, b, A % UMFPack_Numeric, Control, iInfo )
  END IF

  IF (iinfo(1)<0) THEN
    PRINT*, 'Error occurred in umf4sol: ', iinfo(1)
    STOP
  END IF
 
  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    IF ( BigMode ) THEN
      CALL umf4_l_fnum(A % UMFPack_Numeric)
    ELSE
      CALL umf4fnum(A % UMFPack_Numeric)
    END IF
    A % UMFPack_Numeric = 0
  END IF
#else
   CALL Fatal( 'UMFPack_SolveSystem', 'UMFPACK Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE UMFPack_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using Cholmod multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE Cholmod_SolveSystem( Solver,A,x,b,Free_fact)
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: x(*), b(*)

  INTEGER(KIND=AddrInt) :: cholmod_ffactorize

  LOGICAL :: Factorize, FreeFactorize, Found

  REAL(KIND=dp), POINTER CONTIG :: Vals(:)
  INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

#ifdef HAVE_CHOLMOD
  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % Cholmod/=0 ) THEN
        CALL cholmod_ffree(A % cholmod)
        A % cholmod = 0
      END IF
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', Found )
  IF ( .NOT. Found ) Factorize = .TRUE.

  IF ( Factorize .OR. A% cholmod==0 ) THEN
    IF ( A % cholmod/=0 ) THEN
      CALL cholmod_ffree(A % cholmod)
      A % cholmod = 0
    END IF

    Rows => A % Rows
    Cols => A % Cols
    Vals => A % Values

    Rows=Rows-1; Cols=Cols-1 ! c numbering
    A % Cholmod=cholmod_ffactorize(A % NumberOfRows, Rows, Cols, Vals)
    Rows=Rows+1; Cols=Cols+1 ! fortran numbering
  END IF

  CALL cholmod_fsolve(A % cholmod, A % NumberOfRows, x, b);

  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', Found )
  IF ( .NOT. Found ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    CALL cholmod_ffree(A % cholmod)
    A % cholmod = 0
  END IF
#else
   CALL Fatal( 'Cholmod_SolveSystem', 'Cholmod Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Cholmod_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using SuiteSparseQR multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE SPQR_SolveSystem( Solver,A,x,b,Free_fact)
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: x(*), b(*)

  INTEGER(KIND=AddrInt) :: spqr_ffactorize

  INTEGER :: i,spqr_ffree
  LOGICAL :: Factorize, FreeFactorize, Found

  REAL(KIND=dp), POINTER CONTIG :: Vals(:)
  INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

#ifdef HAVE_CHOLMOD
  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % Cholmod/=0 ) THEN
        IF(spqr_ffree(A % cholmod)==0) A % Cholmod=0
      END IF
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', Found )
  IF ( .NOT. Found ) Factorize = .TRUE.

  IF ( Factorize .OR. A% cholmod==0 ) THEN
    IF ( A % cholmod/=0 ) THEN
      i=spqr_ffree(A % cholmod)
      A % cholmod = 0
    END IF

    Rows => A % Rows
    Cols => A % Cols
    Vals => A % Values

    Rows=Rows-1; Cols=Cols-1 ! c numbering
    A % Cholmod=spqr_ffactorize(A % NumberOfRows, Rows, Cols, Vals)
    Rows=Rows+1; Cols=Cols+1 ! fortran numbering
  END IF

  CALL spqr_fsolve(A % cholmod, A % NumberOfRows, x, b);

  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', Found )
  IF ( .NOT. Found ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    i=spqr_ffree(A % cholmod)
    A % cholmod = 0
  END IF
#else
   CALL Fatal( 'SPQR_SolveSystem', 'SPQR Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE SPQR_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using MUMPS direct solver. This is a legacy solver
!> with complicated dependencies. This is only available in parallel. 
!------------------------------------------------------------------------------
  SUBROUTINE Mumps_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
 
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_MUMPS
  INCLUDE 'mpif.h'

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,ip,ierr,icntlft,nzloc
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER, ALLOCATABLE :: memb(:)
  INTEGER :: Comm_active, Group_active, Group_world

  REAL(KIND=dp), ALLOCATABLE :: dbuf(:)


  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( ASSOCIATED(A % MumpsID) ) THEN
        DEALLOCATE( A % MumpsID % irn_loc, &
           A % MumpsID % jcn_loc, A % MumpsID % rhs,  &
              A % MumpsID % isol_loc, A % MumpsID % sol_loc, A % Gorder)
        A % Gorder=>Null()

        A % MumpsID % job = -2
        CALL DMumps(A % MumpsID)
        DEALLOCATE(A % MumpsID)
      END IF
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.ASSOCIATED(A % MumpsID) ) THEN
    IF ( ASSOCIATED(A % MumpsID) ) THEN
      DEALLOCATE( A % MumpsID % irn_loc, &
         A % MumpsID % jcn_loc, A % MumpsID % Rhs,  &
            A % MumpsID % isol_loc, A % MumpsID % sol_loc, A % Gorder)
      A % MumpsID % job = -2
      CALL DMumps(A % MumpsID)
      DEALLOCATE(A % MumpsID)
    END IF

    ALLOCATE(A % MumpsID)

    A % MumpsID % Comm = A % Comm
    A % MumpsID % par  =  1
    A % MumpsID % job  = -1
    A % MumpsID % Keep =  0

    matsym = ListGetLogical( Solver % Values, 'Linear System Symmetric', stat)
    matspd = ListGetLogical( Solver % Values, 'Linear System Positive Definite', stat)

    ! force unsymmetric mode when "row equilibration" is used
    scaled = ListGetLogical( Solver % Values, 'Linear System Scaling', stat)
    IF(.NOT.stat) scaled = .TRUE.
    IF(scaled) THEN
      IF(ListGetLogical( Solver % Values, 'Linear System Row Equilibration',stat)) matsym=.FALSE.
    END IF

    IF(matsym) THEN
      IF ( matspd) THEN
        A % MumpsID % sym = 1
      ELSE
        A % MumpsID % sym = 0 ! 2=symmetric, but unsymmetric solver seems faster, at least in a few
                              ! simple cases...  more testing needed...
      END IF
    ELSE
      A % MumpsID % sym = 0
    END IF

    CALL DMumps(A % MumpsID)

    IF(ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)

    IF(ASSOCIATED(A % ParallelInfo)) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)

      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, &
          A % Perm, A % Gorder, Owner )

      CALL MPI_ALLREDUCE( SUM(Owner), A % MumpsID % n, &
         1, MPI_INTEGER, MPI_SUM, A % MumpsID % Comm, ierr )
      DEALLOCATE(Owner)
    ELSE
      CALL MPI_ALLREDUCE( A % NumberOfRows, A % MumpsId % n, &
            1, MPI_INTEGER, MPI_MAX, A % Comm, ierr )

      ALLOCATE(A % Gorder(A % NumberOFrows))
      DO i=1,A % NumberOFRows
        A % Gorder(i) = i
      END DO
    END IF

   ! Set matrix for Mumps (unsymmetric case)
    IF (A % mumpsID % sym == 0) THEN
      A % MumpsID % nz_loc = A % Rows(A % NumberOfRows+1)-1

      ALLOCATE( A % MumpsID % irn_loc(A % MumpsID % nz_loc) )
      DO i=1,A % NumberOfRows
        ip = A % Gorder(i)
        DO j=A % Rows(i),A % Rows(i+1)-1
          A % MumpsID % irn_loc(j) = ip
        END DO
      END DO

      ALLOCATE( A % MumpsID % jcn_loc(A % MumpsId % nz_loc) )
      DO i=1,A % MumpsID % nz_loc
        A % MumpsID % jcn_loc(i) = A % Gorder(A % Cols(i))
      END DO
      A % MumpsID % a_loc   => A % values
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
        END DO
      END DO

      A % MumpsID % nz_loc = nzloc

      ALLOCATE( A % MumpsID % irn_loc(A % MumpsID % nz_loc) )
      ALLOCATE( A % MumpsID % jcn_loc(A % MumpsId % nz_loc) )
      ALLOCATE( A % MumpsID % A_loc(A % MumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
          A % mumpsID % IRN_loc(nzloc) = A % Gorder(i)
          A % mumpsID % A_loc(nzloc) = A % Values(j)
          A % mumpsID % JCN_loc(nzloc) = A % Gorder(A % Cols(j))
        END DO
      END DO
    END IF


    ALLOCATE(A % MumpsID % rhs(A % MumpsId % n))

    A % MumpsID % icntl(2)  = 0 ! suppress printing of diagnostics and warnings
    A % MumpsID % icntl(3)  = 0 ! suppress statistics
    A % MumpsID % icntl(4)  = 1 ! the same as the two above, but doesn't seem to work.
    A % MumpsID % icntl(5)  = 0 ! matrix format 'assembled'

    icntlft = ListGetInteger(Solver % Values, &
          'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % MumpsID % icntl(14) = icntlft
    END IF
    A % MumpsID % icntl(18) = 3 ! 'distributed' matrix 
    A % MumpsID % icntl(21) = 1 ! 'distributed' solution phase

    A % MumpsID % job = 4
    CALL DMumps(A % MumpsID)
    CALL Flush(6)

    A % MumpsID % lsol_loc = A % mumpsid % info(23)
    ALLOCATE(A % MumpsID % sol_loc(A % MumpsId % lsol_loc))
    ALLOCATE(A % MumpsID % isol_loc(A % MumpsId % lsol_loc))
  END IF

 ! sum the rhs from all procs. Could be done
 ! for neighbours only (i guess):
 ! ------------------------------------------
  A % MumpsID % RHS = 0
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    A % MumpsId % RHS(ip) = b(i)
  END DO
  ALLOCATE( dbuf(A % MumpsID % n) )
  dbuf = A % MumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % MumpsID % RHS, &
    A % MumpsID % n, MPI_DOUBLE_PRECISION, MPI_SUM, A % MumpsID % Comm, ierr )

 ! Solution:
 ! ---------
  A % MumpsID % job = 3
  CALL DMumps(A % MumpsID)

 ! Distribute the solution to all:
 ! -------------------------------
  A % MumpsId % Rhs = 0
  DO i=1,A % MumpsID % lsol_loc
    A % MumpsID % RHS(A % MumpsID % isol_loc(i)) = &
            A % MumpsID % sol_loc(i)
  END DO
  dbuf = A % MumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % MumpsID % RHS, &
    A % MumpsID % N, MPI_DOUBLE_PRECISION, MPI_SUM,A %  MumpsID % Comm, ierr )

  DEALLOCATE(dbuf)

 ! Select the values which belong to us:
 ! -------------------------------------
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    x(i) = A % MumpsId % RHS(ip)
  END DO

  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    DEALLOCATE( A % MumpsID % irn_loc, &
       A % MumpsID % jcn_loc, A % MumpsID % Rhs,  &
         A % MumpsID % isol_loc, A % MumpsID % sol_loc, A % Gorder)

    IF ( A % MumpsID % Sym/=0 ) DEALLOCATE( A % MumpsID % A_loc )

    A % Gorder=>Null()

    A % MumpsID % job = -2
    CALL DMumps(A % MumpsID)
    DEALLOCATE(A % MumpsID)
    A % MumpsId => NULL()
  END IF

#else
   CALL Fatal( 'Mumps_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Mumps_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves local linear system using MUMPS direct solver. If the solved system
!> is singular, optionally one possible solution is returned.
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_SolveSystem( Solver, A, x, b, Free_Fact )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp), TARGET :: x(*), b(*)
     LOGICAL, OPTIONAL :: Free_Fact

     INTEGER :: i
     LOGICAL :: Factorize, FreeFactorize, stat

#ifdef HAVE_MUMPS
     ! Free local Mumps instance if requested
     IF (PRESENT(Free_Fact)) THEN
       IF (Free_Fact) THEN
         CALL MumpsLocal_Free(A)
         RETURN
       END IF
     END IF

     ! Refactorize local matrix if needed
     Factorize = ListGetLogical( Solver % Values, &
                                'Linear System Refactorize', stat )
     IF (.NOT. stat) Factorize = .TRUE.

     IF (Factorize .OR. .NOT. ASSOCIATED(A % mumpsIDL)) THEN
       CALL MumpsLocal_Factorize(Solver, A)
     END IF

     ! Set RHS
     A % mumpsIDL % NRHS = 1
     A % mumpsIDL % LRHS = A % mumpsIDL % n
     DO i=1,A % NumberOfRows
       A % mumpsIDL % RHS(i) = b(i)
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, b, 1, A % mumpsIDL % RHS, 1)

     ! SOLUTION PHASE
     A % mumpsIDL % job = 3
     CALL DMumps(A % mumpsIDL)

     ! TODO: If solution is not local, redistribute the solution vector here

     ! Set local solution
     DO i=1,A % NumberOfRows
       x(i)=A % mumpsIDL % RHS(i)
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, A % mumpsIDL % RHS, 1, x, 1)

     FreeFactorize = ListGetLogical( Solver % Values, &
                                  'Linear System Free Factorization', stat )
     IF (.NOT. stat) FreeFactorize = .TRUE.

     IF (Factorize .AND. FreeFactorize) CALL MumpsLocal_Free(A)
#else
     CALL Fatal( 'MumpsLocal_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Factorize local matrix with Mumps
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_Factorize(Solver, A)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Solver_t) :: Solver
       TYPE(Matrix_t) :: A

#ifdef HAVE_MUMPS
       INCLUDE 'mpif.h'

     INTEGER :: i, j, n, nz, allocstat, icntlft, ptype, nzloc
     LOGICAL :: matpd, matsym, nullpiv, stat

     ! INTEGER :: myrank, ierr
     ! CHARACTER(len=32) :: buf

    IF ( ASSOCIATED(A % mumpsIDL) ) THEN
         CALL MumpsLocal_Free(A)
    END IF

    ALLOCATE(A % mumpsIDL)

    ! INITIALIZATION PHASE

    ! Initialize local instance of Mumps
    ! TODO (Hybridization): change this if local system needs to be solved
    ! with several cores
    A % mumpsIDL % COMM = MPI_COMM_SELF
    A % mumpsIDL % PAR = 1 ! Host (=self) takes part in factorization

    ! Check if matrix is symmetric or spd
    matsym = ListGetLogical(Solver % Values, &
                            'Linear System Symmetric', stat)
    matpd = ListGetLogical(Solver % Values, &
                           'Linear System Positive Definite', stat)

    A % mumpsIDL % SYM = 0
    IF (matsym) THEN
      IF (matpd) THEN
        ! Matrix is symmetric positive definite
        A % mumpsIDL % SYM = 1
      ELSE
        ! Matrix is symmetric
        A % mumpsIDL % SYM = 2
     END IF
    ELSE
      ! Matrix is unsymmetric
      A % mumpsIDL % SYM = 0
    END IF
    A % mumpsIDL % JOB  = -1 ! Initialize
    CALL DMumps(A % mumpsIDL)

    ! FACTORIZE PHASE

    ! Set stdio parameters
    A % mumpsIDL % ICNTL(1)  = 6  ! Error messages to stdout
    A % mumpsIDL % ICNTL(2)  = -1 ! No diagnostic and warning messages
    A % mumpsIDL % ICNTL(3)  = -1 ! No statistic messages
    A % mumpsIDL % ICNTL(4)  = 1  ! Print only error messages

    ! Set matrix format
    A % mumpsIDL % ICNTL(5)  = 0  ! Assembled matrix format
    A % mumpsIDL % ICNTL(18) = 0 ! Centralized matrix
    A % mumpsIDL % ICNTL(21) = 0 ! Centralized dense solution phase

    ! Check if solution of singular systems is ok
    A % mumpsIDL % ICNTL(24) = 0
    nullpiv = ListGetLogical(Solver % Values, &
                            'Mumps Solve Singular', stat)
    IF (nullpiv) THEN
      A % mumpsIDL % ICNTL(24) = 1
      A % mumpsIDL % CNTL(1) = 1D-2     ! Pivoting threshold
      A % mumpsIDL % CNTL(3) = 1D-9     ! Null pivot detection threshold
      A % mumpsIDL % CNTL(5) = 1D6      ! Fixation value for null pivots
      A % mumpsIDL % CNTL(13) = 1       ! Do not use ScaLAPACK on the root node
      ! TODO: if needed, here set CNTL(3) and CNTL(5) as parameters for
      ! more accurate null pivot detection
    END IF

    ! Set permutation strategy for Mumps
    ptype = ListGetInteger(Solver % Values, &
                                'Mumps Permutation Type', stat)
    IF (stat) THEN
      A % mumpsIDL % ICNTL(6) = ptype
    END IF

    ! TODO: Change this if system is larger than local.
    ! For larger than local systems define global->local numbering
    n = A % NumberofRows
    nz = A % Rows(A % NumberOfRows+1)-1
    A % mumpsIDL % N  = n
    A % mumpsIDL % NZ = nz
    ! A % mumpsIDL % nz_loc = nz

    ! Allocate rows and columns for MUMPS
    ALLOCATE( A % mumpsIDL % IRN(nz), &
          A % mumpsIDL % JCN(nz), &
          A % mumpsIDL % A(nz), STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('MumpsLocal_Factorize', &
            'Memory allocation for MUMPS row and column indices failed.')
    END IF

    ! Set matrix for Mumps (unsymmetric case)
    IF (A % mumpsIDL % sym == 0) THEN
      DO i=1,A % NumberOfRows
         DO j=A % Rows(i),A % Rows(i+1)-1
            A % mumpsIDL % IRN(j) = i
         END DO
       END DO

       ! Set columns and values
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % JCN(i) = A % Cols(i)
       END DO
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % A(i) = A % Values(i)
       END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows
        DO j=A % Rows(i),A % Rows(i+1)-1
          ! Only output lower triangular part to Mumps
          IF (i<=A % Cols(j)) THEN
            nzloc = nzloc + 1
            A % mumpsIDL % IRN(nzloc) = i
            A % mumpsIDL % JCN(nzloc) = A % Cols(j)
            A % mumpsIDL % A(nzloc) = A % Values(j)
          END IF
        END DO
      END DO
    END IF

    icntlft = ListGetInteger(Solver % Values, 'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % mumpsIDL % ICNTL(14) = icntlft
    END IF

    A % mumpsIDL % JOB = 1 ! Perform analysis
    CALL DMumps(A % mumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('MumpsLocal_Factorize','Mumps analysis phase failed')
    END IF

    A % mumpsIDL % JOB = 2 ! Perform factorization
    CALL DMumps(A % mumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('MumpsLocal_Factorize','Mumps factorize phase failed')
    END IF

    ! Allocate RHS
    ALLOCATE(A % mumpsIDL % RHS(A % mumpsIDL % N), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('MumpsLocal_Factorize', &
                 'Memory allocation for MUMPS solution vector and RHS failed.' )
    END IF
#else
   CALL Fatal( 'MumpsLocal_Factorize', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_Factorize
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve local nullspace using MUMPS direct solver. On exit, z will be
!> allocated and will hold the jth local nullspace vectors as z(j,:).
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_SolveNullSpace(Solver, A, z, nz)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Solver_t) :: Solver
      TYPE(Matrix_t) :: A
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:), TARGET :: z
      INTEGER :: nz, nrhs

#ifdef HAVE_MUMPS
      INCLUDE 'mpif.h'

      INTEGER :: j,k,n, allocstat
      LOGICAL :: Factorize, FreeFactorize, stat

      REAL(KIND=dp), DIMENSION(:), POINTER :: rhs_tmp, nspace

      Factorize = ListGetLogical( Solver % Values,&
                                  'Linear System Refactorize', stat )
      IF (.NOT. stat) Factorize = .TRUE.

      ! Refactorize local matrix if needed
      IF ( Factorize .OR. (.NOT. ASSOCIATED(A % mumpsIDL)) ) THEN
            CALL MumpsLocal_Factorize(Solver, A)
      END IF

      ! Check symmetry condition
      IF (.NOT. (A % mumpsIDL % sym == 1 .OR. A % mumpsIDL % sym == 2)) THEN
         CALL Fatal( 'MumpsLocal_SolveNullSpace', &
           'Mumps null space computation not supported for unsymmetric matrices.')
      END IF

      ! Get the size of the local nullspace and allocate memory
      ! TODO: this is incorrect, should be number of columns
      n = A % NumberOfRows
      nz = A % mumpsIDL % INFOG(28)

      IF (nz == 0) THEN
            ALLOCATE(z(nz,n))
          RETURN
      END IF

      ALLOCATE(nspace(n*nz), STAT=allocstat)
      IF (allocstat /= 0) THEN
         CALL Fatal( 'MumpsLocal_SolveNullSpace', &
               'Storage allocation for local nullspace failed.')
      END IF

      rhs_tmp => A % mumpsIDL % RHS
      nrhs = A % mumpsIDL % NRHS
      A % mumpsIDL % RHS => nspace
      A % mumpsIDL % NRHS = nz

      ! Solve for the complete local nullspace
      A % mumpsIDL % JOB = 3
      A % mumpsIDL % ICNTL(25) = -1
      CALL DMumps(A % mumpsIDL)
      ! Check return status
      IF (A % mumpsIDL % INFO(1)<0) THEN
          CALL Fatal('MumpsLocal_SolveNullSpace','Mumps nullspace solution failed')
      END IF

      ! Restore pointer to local solution vector
      A % mumpsIDL % ICNTL(25) = 0
      A % mumpsIDL % RHS => rhs_tmp
      A % mumpsIDL % NRHS = nrhs

      FreeFactorize = ListGetLogical( Solver % Values, &
            'Linear System Free Factorization', stat )
      IF (.NOT. stat) FreeFactorize = .TRUE.

      IF (Factorize .AND. FreeFactorize) THEN
          CALL MumpsLocal_Free(A)
       END IF

      ALLOCATE(z(nz,n), STAT=allocstat)
      IF (allocstat /= 0) THEN
          CALL Fatal( 'MumpsLocal_SolveNullSpace', &
                           'Storage allocation for local nullspace failed.')
      END IF

      ! Copy computed nullspace to z and deallocate nspace
      DO j=1,nz
            ! z is row major, cannot use DCOPY here
        DO k=1,n
              z(j,k)=nspace(n*(j-1)+k)
        END DO
      END DO
      DEALLOCATE(nspace)
#else
   CALL Fatal( 'MumpsLocal_SolveNullSpace', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_SolveNullSpace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Free local Mumps variables and solver internal allocations
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_Free(A)
!------------------------------------------------------------------------------
        IMPLICIT NONE

        TYPE(Matrix_t) :: A

#ifdef HAVE_MUMPS
     IF (ASSOCIATED(A % mumpsIDL)) THEN
          ! Deallocate Mumps structures
          DEALLOCATE( A % mumpsIDL % irn, A % mumpsIDL % jcn, &
             A % mumpsIDL % a, A % mumpsIDL % rhs)

          ! Free Mumps internal allocations
          A % mumpsIDL % job = -2
          CALL DMumps(A % mumpsIDL)
          DEALLOCATE(A % mumpsIDL)

          A % mumpsIDL => Null()
        END IF
#else
     CALL Fatal( 'MumpsLocal_Free', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_Free
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solves a linear system using SuperLU direct solver.
!------------------------------------------------------------------------------
  SUBROUTINE SuperLU_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_SUPERLU
      LOGICAL :: stat, Factorize, FreeFactorize
      integer :: n, nnz, nrhs, iinfo, iopt, nprocs

      interface
        subroutine solve_superlu( iopt, nprocs, n, nnz, nrhs, values, cols, &
                     rows, b, ldb, factors, iinfo )
            use types
            integer :: iopt, nprocs, n, nnz, nrhs, cols(*), rows(*), ldb, iinfo
            real(kind=dp) :: values(*), b(*)
            integer(kind=addrint) :: factors
        end subroutine solve_superlu
      end interface

      IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
          IF ( A % SuperLU_Factors/= 0 ) THEN
            iopt = 3
            CALL Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, &
               A % Cols, A % Rows, x, n, A % SuperLU_Factors, iinfo )
            A % SuperLU_Factors = 0
          END IF
          RETURN
        END IF
      END IF

      n = A % NumberOfRows
      nrhs = 1
      x(1:n) = b(1:n)
      nnz = A % Rows(n+1)-1

      nprocs = ListGetInteger( Solver % Values, & 
              'Linear System Number of Threads', stat )
      IF ( .NOT. stat ) nprocs = 1
!
      Factorize = ListGetLogical( Solver % Values, &
         'Linear System Refactorize', stat )
      IF ( .NOT. stat ) Factorize = .TRUE.

      IF ( Factorize .OR. A % SuperLU_Factors==0 ) THEN

        IF ( A % SuperLU_Factors/= 0 ) THEN
          iopt = 3
          call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
                     A % Rows, x, n, A % SuperLU_Factors, iinfo )
          A % SuperLU_Factors=0
        END IF

        ! First, factorize the matrix. The factors are stored in *factors* handle.
        iopt = 1
        call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
              A % Rows, x, n, A % SuperLU_Factors, iinfo )
 
        if (iinfo .eq. 0) then
           write (*,*) 'Factorization succeeded'
        else
           write(*,*) 'INFO from factorization = ', iinfo
        endif
      END IF

      !
      ! Second, solve the system using the existing factors.
      iopt = 2
      call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
           A % Rows,  x, n, A % SuperLU_Factors, iinfo )
!
      if (iinfo .eq. 0) then
         write (*,*) 'Solve succeeded'
      else
         write(*,*) 'INFO from triangular solve = ', iinfo
      endif

! Last, free the storage allocated inside SuperLU
      FreeFactorize = ListGetLogical( Solver % Values, &
          'Linear System Free Factorization', stat )
      IF ( .NOT. stat ) FreeFactorize = .TRUE.

      IF ( Factorize .AND. FreeFactorize ) THEN
        iopt = 3
        call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
                    A % Rows, x, n, A % SuperLU_Factors, iinfo )
        A % SuperLU_Factors = 0
      END IF
#endif      
!------------------------------------------------------------------------------
  END SUBROUTINE SuperLU_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Permon solver
!------------------------------------------------------------------------------
  SUBROUTINE Permon_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
#ifdef HAVE_FETI4I
   use feti4i
#endif
 
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_FETI4I
  INCLUDE 'mpif.h'

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,nd,ip,ierr,icntlft,nzloc
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER :: n_dof_partition

  INTEGER, ALLOCATABLE :: memb(:), DirichletInds(:), Neighbours(:), IntOptions(:)
  REAL(KIND=dp), ALLOCATABLE :: DirichletVals(:), RealOptions(:)
  INTEGER :: Comm_active, Group_active, Group_world

  REAL(KIND=dp), ALLOCATABLE :: dbuf(:)


  n_dof_partition = A % NumberOfRows

!  INTERFACE
!     FUNCTION Permon_InitSolve(n, gnum, nd, dinds, dvals, n_n, n_ranks) RESULT(handle) BIND(c,name='permon_initsolve') 
!        USE, INTRINSIC :: ISO_C_BINDING
!        TYPE(C_PTR) :: handle
!        INTEGER(C_INT), VALUE :: n, nd, n_n
!        REAL(C_DOUBLE) :: dvals(*)
!        INTEGER(C_INT) :: gnum(*), dinds(*), n_ranks(*)
!     END FUNCTION Permon_Initsolve
!
!     SUBROUTINE Permon_Solve( handle, x, b ) BIND(c,name='permon_solve')
!        USE, INTRINSIC :: ISO_C_BINDING
!        REAL(C_DOUBLE) :: x(*), b(*)
!        TYPE(C_PTR), VALUE :: handle
!     END SUBROUTINE Permon_solve
!  END INTERFACE

  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.C_ASSOCIATED(A % PermonSolverInstance) ) THEN
    IF ( C_ASSOCIATED(A % PermonSolverInstance) ) THEN
       CALL Fatal( 'Permon', 're-entry not implemented' )
    END IF

    nd = COUNT(A % ConstrainedDOF)
    ALLOCATE(DirichletInds(nd), DirichletVals(nd))
    j = 0
    DO i=1,A % NumberOfRows
      IF(A % ConstrainedDOF(i)) THEN
        j = j + 1
        DirichletInds(j) = i; DirichletVals(j) = A % Dvalues(i)
      END IF
    END DO

    !TODO sequential case not working
    n = 0
    ALLOCATE(neighbours(Parenv % PEs))
    DO i=1,ParEnv % PEs
      IF( ParEnv % IsNeighbour(i) .AND. i-1/=ParEnv % myPE) THEN
        n = n + 1
        neighbours(n) = i-1
      END IF
    END DO

    !A % PermonSolverInstance = Permon_InitSolve( SIZE(A % ParallelInfo % GlobalDOFs), &
    !     A % ParallelInfo % GlobalDOFs, nd,  DirichletInds, DirichletVals, n, neighbours )

    IF( n_dof_partition /=  SIZE(A % ParallelInfo % GlobalDOFs) ) THEN
      CALL Fatal( 'Permon', &
        'inconsistency: A % NumberOfRows /=  SIZE(A % ParallelInfo % GlobalDOFs' )
    END IF

    ALLOCATE(IntOptions(10))
    ALLOCATE(RealOptions(1))

    CALL FETI4ISetDefaultIntegerOptions(IntOptions)
    CALL FETI4ISetDefaultRealOptions(RealOptions)

    CALL FETI4ICreateInstance(A % PermonSolverInstance, A % PermonMatrix, &
      A % NumberOfRows, b, A % ParallelInfo % GlobalDOFs, &
      n, neighbours, &
      nd, DirichletInds, DirichletVals, &
      IntOptions, RealOptions)
  END IF

  !CALL Permon_Solve( A % PermonSolverInstance, x, b )
  CALL FETI4ISolve(A % PermonSolverInstance, n_dof_partition, x)
#else
   CALL Fatal( 'Permon_SolveSystem', 'Permon Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Permon_SolveSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solves a linear system using Pardiso direct solver (from either MKL or
!> official Pardiso distribution. If possible, MKL-version is used).
!------------------------------------------------------------------------------
  SUBROUTINE Pardiso_SolveSystem( Solver,A,x,b,Free_fact )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), TARGET :: x(*), b(*)
    LOGICAL, OPTIONAL :: Free_fact

! MKL version of Pardiso (interface is different)
#if defined(HAVE_MKL)
    INTERFACE
      SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
                           values, rows, cols, perm, nrhs, iparm, msglvl, b, x, ierror)
        USE Types
        IMPLICIT NONE
        REAL(KIND=dp) :: values(*), b(*), x(*)
        INTEGER(KIND=AddrInt) :: pt(*)
        INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
        INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
      END SUBROUTINE pardiso

      SUBROUTINE pardisoinit(pt, mtype, iparm)
        USE Types
        IMPLICIT NONE
        INTEGER(KIND=AddrInt) :: pt(*)
        INTEGER :: mtype
        INTEGER :: iparm(*)
      END SUBROUTINE pardisoinit
    END INTERFACE

    INTEGER maxfct, mnum, mtype, phase, n, nrhs, ierror, msglvl
    INTEGER, POINTER :: Iparm(:)
    INTEGER i, j, k, nz, idum(1), nzutd
    LOGICAL :: Found, matsym, matpd
    REAL*8  :: ddum(1)

    LOGICAL :: Factorize, FreeFactorize
    INTEGER :: tlen, allocstat
    CHARACTER(LEN=MAX_NAME_LEN) :: threads, mat_type

    REAL(KIND=dp), POINTER :: values(:)
    INTEGER, POINTER  :: rows(:), cols(:)

    ! Check if system needs to be refactorized
    Factorize = ListGetLogical( Solver % Values, &
                  'Linear System Refactorize', Found )
    IF ( .NOT. Found ) Factorize = .TRUE.

    ! Set matrix type for Pardiso
    mat_type = ListGetString( Solver % Values, 'Linear System Matrix Type', Found )
    
    IF (Found) THEN
      SELECT CASE(mat_type)
      CASE('positive definite')
        mtype = 2
      CASE('symmetric indefinite')
        mtype = -2
      CASE('structurally symmetric')
        mtype = 1
      CASE('nonsymmetric', 'general')
        mtype = 11
      CASE DEFAULT
        mtype = 11
      END SELECT
    ELSE
      ! Check if matrix is symmetric or spd
      matsym = ListGetLogical(Solver % Values, &
                            'Linear System Symmetric', Found)

      matpd = ListGetLogical(Solver % Values, &
                    'Linear System Positive Definite', Found)

      IF (matsym) THEN
        IF (matpd) THEN
          ! Matrix is symmetric positive definite
           mtype = 2
         ELSE
          ! Matrix is structurally symmetric (can't handle indefinite systems!!!!)
          mtype = 1
        END IF
      ELSE
        ! Matrix is unsymmetric
        mtype = 11
      END IF
    END IF

    ! Free factorization if requested
    IF ( PRESENT(Free_Fact) ) THEN
      IF ( Free_Fact ) THEN
        IF(ASSOCIATED(A % PardisoId)) THEN
          phase     = -1           ! release internal memory
          n = A % NumberOfRows
          maxfct = 1
          mnum   = 1
          nrhs   = 1
          msglvl = 0
          CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
              ddum, idum, idum, idum, nrhs, A % PardisoParam, msglvl, ddum, ddum, ierror)
          DEALLOCATE(A % PardisoId, A % PardisoParam)
          A % PardisoId => NULL()
          A % PardisoParam => NULL()
        END IF
        RETURN
      END IF
    END IF

    ! Get number of rows and number of nonzero elements
    n = A % Numberofrows
    nz = A % Rows(n+1)-1

    ! Copy upper triangular part of symmetric positive definite system
    IF ( ABS(mtype) == 2 ) THEN
      nzutd = 0
      DO i=1,n
        nzutd = nzutd + A % Rows(i+1)-A % Diag(i)
      END DO
      
      ALLOCATE( values(nzutd), cols(nzutd), rows(n+1), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('Pardiso_SolveSystem', &
           'Memory allocation for row and column indices failed')
      END IF

      ! Copy upper triangular part of A to Pardiso structure
      Rows(1)=1
      DO i=1,n
        ! Set up row pointers and copy values
        nzutd = A % Rows(i+1)-A % Diag(i)
        Rows(i+1) = Rows(i)+nzutd
        DO j=0,nzutd-1
          Cols(Rows(i)+j)=A % Cols(A%Diag(i)+j)
          Values(Rows(i)+j)=A % Values(A%Diag(i)+j)
        END DO
      END DO
      
    ELSE
      Cols => A % Cols
      Rows => A % Rows
      Values => A % Values
    END IF

    ! Set up parameters
    msglvl    = 0 ! Do not write out any info
    maxfct    = 1 ! Set up space for 1 matrix at most
    mnum      = 1 ! Matrix to use in the solution phase (1st and only one)
    nrhs      = 1 ! Use only one RHS
    iparm => A % PardisoParam

    ! Compute factorization if necessary
    IF ( Factorize .OR. .NOT.ASSOCIATED(A % PardisoID) ) THEN
      ! Free factorization
      IF (ASSOCIATED(A % PardisoId)) THEN
        phase = -1 ! release internal memory
        CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
                    ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
        DEALLOCATE(A % PardisoId, A % PardisoParam)
        A % PardisoId => NULL()
        A % PardisoParam => NULL()
      END IF

      ! Allocate pardiso internal structures
      ALLOCATE(A % PardisoId(64), A % PardisoParam(64), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('Pardiso_SolveSystem', &
                     'Memory allocation for Pardiso failed')
      END IF
      iparm => A % PardisoParam
      iparm = 0
      A % PardisoId = 0
 
      ! Set up scaling values for solver based on matrix type
      CALL pardisoinit(A % PardisoId, mtype, iparm)

      ! Set up rest of parameters explicitly
      iparm(1)=1      ! Do not use solver default parameters
      iparm(2)=2      ! Minimize fill-in with nested dissection from Metis
      iparm(3)=0      ! Reserved
      iparm(4)=0      ! Always compute factorization
      iparm(5)=0      ! No user input permutation
      iparm(6)=0      ! Write solution vector to x
      iparm(8)=5      ! Number of iterative refinement steps
      iparm(18)=-1    ! Report nnz(L) and nnz(U)
      iparm(19)=0     ! Do not report Mflops
      iparm(27)=0     ! Do not check sparse matrix representation
      iparm(35)=0     ! Use Fortran style indexing
      iparm(60)=0     ! Use in-core version of Pardiso

      ! Perform analysis
      phase = 11            ! Analysis
      CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
            values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      IF (ierror .NE. 0) THEN
        WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
        CALL Fatal('Pardiso_SolveSystem','Error during analysis phase')
      END IF

      ! Perform factorization
      phase = 22            ! Factorization
      CALL pardiso (A % pardisoId, maxfct, mnum, mtype, phase, n, &
           values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      IF (ierror .NE. 0) THEN
        WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
        CALL Fatal('Pardiso_SolveSystem','Error during factorization phase')
      END IF
    END IF ! Compute factorization

    ! Perform solve
    phase = 33            ! Solve, iterative refinement
    CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
                values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror)

    IF (ierror .NE. 0) THEN
      WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
      CALL Fatal('Pardiso_SolveSystem','Error during solve phase')
    END IF

    ! Release memory if needed
    FreeFactorize = ListGetLogical( Solver % Values, &
                       'Linear System Free Factorization', Found )
    IF ( .NOT. Found ) FreeFactorize = .TRUE.

    IF ( Factorize .AND. FreeFactorize ) THEN
      phase     = -1           ! release internal memory
      CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
            ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      DEALLOCATE(A % PardisoId, A % PardisoParam)
      A % PardisoId => NULL()
      A % PardisoParam => NULL()
    END IF
    IF (ABS(mtype) == 2) DEALLOCATE(Values, Rows, Cols)

! Distribution version of Pardiso
#elif defined(HAVE_PARDISO)
!..   All other variables

      INTERFACE
        SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
          values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror, dparm)
          USE Types
          REAL(KIND=dp) :: values(*), b(*), x(*), dparm(*)
          INTEGER(KIND=AddrInt) :: pt(*)
          INTEGER :: idum(*), nrhs, iparm(*), msglvl, ierror
          INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
        END SUBROUTINE pardiso

        SUBROUTINE pardisoinit(pt,mtype,solver,iparm,dparm,ierror)
          USE Types
          INTEGER :: mtype, iparm(*),ierror,solver
          REAL(KIND=dp) :: dparm(*)
          INTEGER(KIND=AddrInt) :: pt(*)
        END SUBROUTINE pardisoinit
      END INTERFACE

      INTEGER maxfct, mnum, mtype, phase, n, nrhs, ierror, msglvl
      INTEGER, POINTER :: Iparm(:)
      INTEGER i, j, k, nz, idum(1)
      LOGICAL :: Found, Symm, Posdef
      REAL*8  waltime1, waltime2, ddum(1), dparm(64)

      LOGICAL :: Factorize, FreeFactorize
      INTEGER :: tlen
      CHARACTER(LEN=MAX_STRING_LEN) :: threads

      REAL(KIND=dp), POINTER :: values(:)
      INTEGER, POINTER  :: rows(:), cols(:)

      IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
          IF(ASSOCIATED(A % PardisoId)) THEN
            phase     = -1           ! release internal memory
            n = A % NumberOfRows
            mtype  = 11
            maxfct = 1
            mnum   = 1
            nrhs   = 1
            msglvl = 0
            CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
             ddum, idum, idum, idum, nrhs, A % PardisoParam, msglvl, ddum, ddum, ierror,dparm)
            DEALLOCATE(A % PardisoId, A % PardisoParam)
            A % PardisoId => NULL()
            A % PardisoParam => NULL()
          END IF
          RETURN
        END IF
      END IF

!  .. Setup Pardiso control parameters und initialize the solvers
!     internal address pointers. This is only necessary for the FIRST
!     call of the PARDISO solver.
!
      Factorize = ListGetLogical( Solver % Values, &
         'Linear System Refactorize', Found )
      IF ( .NOT. Found ) Factorize = .TRUE.

      symm = ListGetLogical( Solver % Values, &
         'Linear System Symmetric', Found )

      posdef = ListGetLogical( Solver % Values, &
         'Linear System Positive Definite', Found )

      mtype = 11

      cols => a % cols
      rows => a % rows
      values => a % values
      n = A % Numberofrows

      IF ( posdef ) THEN
        nz = A % rows(n+1)-1
        allocate( values(nz), cols(nz), rows(n+1) )
        k = 1
        do i=1,n
         rows(i)=k
         do j=a % rows(i), a % rows(i+1)-1
           if ( a % cols(j)>=i .and. a % values(j) /= 0._dp ) then
             cols(k) = a % cols(j)
             values(k) = a % values(j)
             k = k + 1
           end if
         end do
        end do
        rows(n+1)=k
        mtype     = 2
      ELSE IF ( symm ) THEN
        mtype = 1
      END IF

      msglvl    = 0       ! with statistical information
      maxfct    = 1
      mnum      = 1
      nrhs      = 1
      iparm => A % PardisoParam


      IF ( Factorize .OR. .NOT.ASSOCIATED(A % PardisoID) ) THEN
        IF(ASSOCIATED(A % PardisoId)) THEN
          phase     = -1           ! release internal memory
          CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
           ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror,dparm)
          DEALLOCATE(A % PardisoId, A % PardisoParam)
          A % PardisoId => NULL()
          A % PardisoParam => NULL()
        END IF

        ALLOCATE(A % PardisoId(64), A % PardisoParam(64))
        iparm => A % PardisoParam
        CALL pardisoinit(A % PardisoId, mtype, 0, iparm, dparm, ierror )

!..     Reordering and Symbolic Factorization, This step also allocates
!       all memory that is necessary for the factorization
!
        phase     = 11      ! only reordering and symbolic factorization
        msglvl    = 0       ! with statistical information
        maxfct    = 1
        mnum      = 1
        nrhs      = 1

!  ..   Numbers of Processors ( value of OMP_NUM_THREADS )
#ifdef USE_ISO_C_BINDINGS
        CALL envir( 'OMP_NUM_THREADS', threads, tlen )
#else
        CALL envir( 'OMP_NUM_THREADS'//char(0), threads, tlen )
#endif
        iparm(3) = 1
        IF ( tlen>0 ) &
          READ(threads(1:tlen),*) iparm(3)
        IF (iparm(3)<=0) Iparm(3) = 1

        CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
          values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        IF (ierror .NE. 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', ierror
          STOP
        END IF

!..     Factorization.
        phase     = 22  ! only factorization
        CALL pardiso (A % pardisoId, maxfct, mnum, mtype, phase, n, &
         values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        IF (ierror .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', ierror
          STOP
        ENDIF
      END IF

!..   Back substitution and iterative refinement
      phase     = 33  ! only factorization
      iparm(8)  = 0   ! max number of iterative refinement steps
                      ! (if perturbations occur in the factorization
                      ! phase, two refinement steps are taken)

      CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
       values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror, dparm)

!..   Termination and release of memory
      FreeFactorize = ListGetLogical( Solver % Values, &
          'Linear System Free Factorization', Found )
      IF ( .NOT. Found ) FreeFactorize = .TRUE.

      IF ( Factorize .AND. FreeFactorize ) THEN
        phase     = -1           ! release internal memory
        CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
          ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        DEALLOCATE(A % PardisoId, A % PardisoParam)
        A % PardisoId => NULL()
        A % PardisoParam => NULL()
      END IF

      IF(posdef) DEALLOCATE( values, rows, cols )
#else
      CALL Fatal( 'Parsido_SolveSystem', 'Pardiso solver has not been installed.' )
#endif

!------------------------------------------------------------------------------
  END SUBROUTINE Pardiso_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a linear system using Cluster Pardiso direct solver from MKL
!------------------------------------------------------------------------------
  SUBROUTINE CPardiso_SolveSystem( Solver,A,x,b,Free_fact )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), TARGET :: x(*), b(*)
    LOGICAL, OPTIONAL :: Free_fact

! Cluster Pardiso
#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
              values, rows, cols, perm, nrhs, iparm, msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    INTEGER :: phase, n, ierror
    INTEGER, POINTER :: Iparm(:)
    INTEGER i, j, k, nz, nzutd, idum(1), nl, nt
    LOGICAL :: Found, matsym, matpd
    REAL(kind=dp) :: ddum(1)
    REAL(kind=dp), POINTER, DIMENSION(:) :: dbuf

    LOGICAL :: Factorize, FreeFactorize
    INTEGER :: tlen, allocstat

    REAL(KIND=dp), POINTER CONTIG :: values(:)
    INTEGER, POINTER CONTIG :: rows(:), cols(:)

    ! Free factorization if requested
    IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
            CALL CPardiso_Free(A)
        END IF

        RETURN
    END IF

    ! Check if system needs to be refactorized
    Factorize = ListGetLogical( Solver % Values, &
                                'Linear System Refactorize', Found )
    IF ( .NOT. Found ) Factorize = .TRUE.

    ! Compute factorization if necessary
    IF ( Factorize .OR. .NOT.ASSOCIATED(A % CPardisoID) ) THEN
        CALL CPardiso_Factorize(Solver, A)
    END IF

    ! Get global start and end of domain
    nl = A % CPardisoId % iparm(41)
    nt = A % CPardisoId % iparm(42)

    ! Gather RHS
    A % CPardisoId % rhs = 0D0
    DO i=1,A % NumberOfRows
        A % CPardisoId % rhs(A % Gorder(i)-nl+1) = b(i)
    END DO

    ! Perform solve
    phase = 33      ! Solve, iterative refinement
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype,  phase, A % CPardisoId % n, &
          A % CPardisoID % aa, A % CPardisoID % ia, A % CPardisoID % ja, idum, &
          A % CPardisoId % nrhs, A % CPardisoID % iparm, &
          A % CPardisoId % msglvl, A % CPardisoId % rhs, &
          A % CPardisoId % x,  A % Comm, ierror)

    IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during solve phase')
    END IF

    ! Distribute solution
    DO i=1,A % NumberOfRows
        x(i)=A % CPardisoId % x(A % Gorder(i)-nl+1)
    END DO

    ! Release memory if needed
    FreeFactorize = ListGetLogical( Solver % Values, &
                                    'Linear System Free Factorization', Found )
    IF ( .NOT. Found ) FreeFactorize = .TRUE.

    IF ( Factorize .AND. FreeFactorize ) THEN
        CALL CPardiso_Free(A)
    END IF

#else
    CALL Fatal( 'CParsido_SolveSystem', 'Cluster Pardiso solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE CPardiso_SolveSystem
!------------------------------------------------------------------------------

#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
  SUBROUTINE CPardiso_Factorize(Solver, A)
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A

    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
              values, rows, cols, perm, nrhs, iparm, msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    LOGICAL :: matsym, matpd, Found
    INTEGER :: i, j, k, rind, lrow, rptr, rsize, lind, tind
    INTEGER :: allocstat, n, nz, nzutd, nl, nt, nOwned, nhalo, ierror
    INTEGER :: phase, idum(1)
    REAL(kind=dp) :: ddum(1)
    REAL(kind=dp), DIMENSION(:), POINTER :: aa
    INTEGER, DIMENSION(:), POINTER CONTIG :: iparm, ia, ja, owner, dsize, iperm, Order

    INTEGER :: fid
    CHARACTER(LEN=MAX_NAME_LEN) :: mat_type

    ! Free old factorization if necessary
    IF (ASSOCIATED(A % CPardisoId)) THEN
        CALL CPardiso_Free(A)
    END IF

    ! Allocate Pardiso structure
    ALLOCATE(A % CPardisoId, STAT=allocstat)
    IF (allocstat /= 0) THEN
    CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso failed')
    END IF
    ! Allocate control structures
    ALLOCATE(A % CPardisoId% ID(64), A % CPardisoId % IParm(64), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso failed')
    END IF

    iparm => A % CPardisoId % IParm
    ! Initialize control parameters and solver Id's
    DO i=1,64
        iparm(i)=0
    END DO
    DO i=1,64
        A % CPardisoId % ID(i)=0
    END DO

    ! Set matrix type for CPardiso
    mat_type = ListGetString( Solver % Values, 'Linear System Matrix Type', Found )
    IF (Found) THEN
      SELECT CASE(mat_type)
      CASE('positive definite')
        A % CPardisoID % mtype = 2
      CASE('symmetric indefinite')
        A % CPardisoID % mtype = -2
      CASE('structurally symmetric')
        A % CPardisoID % mtype = 1
      CASE('nonsymmetric', 'general')
        A % CPardisoID % mtype = 11
      CASE DEFAULT
        A % CPardisoID % mtype = 11
      END SELECT
    ELSE
      ! Check if matrix is symmetric or spd
      matsym = ListGetLogical(Solver % Values, &
                            'Linear System Symmetric', Found)

      matpd = ListGetLogical(Solver % Values, &
                    'Linear System Positive Definite', Found)

      IF (matsym) THEN
        IF (matpd) THEN
          ! Matrix is symmetric positive definite
          A % CPardisoID % mtype = 2
        ELSE
          ! Matrix is structurally symmetric
          A % CPardisoID % mtype = 1
        END IF
      ELSE
        ! Matrix is nonsymmetric
        A % CPardisoID % mtype = 11
      END IF
    END IF

    ! Set up continuous numbering for the whole computation domain
    n = SIZE(A % ParallelInfo % GlobalDOFs)
    ALLOCATE(A % Gorder(n), Owner(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('CPardiso_Factorize', &
                    'Memory allocation for CPardiso global numbering failed')
    END IF
    CALL ContinuousNumbering(A % ParallelInfo, A % Perm, A % Gorder, Owner, nOwn=nOwned)

    ! Compute the number of global dofs
    CALL MPI_ALLREDUCE(nOwned, A % CPardisoId % n, &
                       1, MPI_INTEGER, MPI_SUM, A % Comm, ierror)
    DEALLOCATE(Owner)

    ! Find bounds of domain
    nl = A % Gorder(1)
    nt = A % Gorder(1)
    DO i=2,n
        ! NOTE: Matrix is structurally symmetric
        rind = A % Gorder(i)
        nl = MIN(rind, nl)
        nt = MAX(rind, nt)
    END DO

    ! Allocate temp storage for global numbering
    ALLOCATE(Order(n), iperm(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                    'Memory allocation for CPardiso global numbering failed')
    END IF

    ! Sort global numbering to build matrix
    Order(1:n) = A % Gorder(1:n)
    DO i=1,n
        iperm(i)=i
    END DO
    CALL SortI(n, Order, iperm)

    ! Allocate storage for CPardiso matrix
    nhalo = (nt-nl+1)-n
    nz = A % Rows(A % NumberOfRows+1)-1
    IF (ABS(A % CPardisoID % mtype) == 2) THEN
        nzutd = ((nz-n)/2)+1 + n
        ALLOCATE(A % CPardisoID % ia(nt-nl+2), &
                 A % CPardisoId % ja(nzutd+nhalo), &
                 A % CPardisoId % aa(nzutd+nhalo), &
                 STAT=allocstat)
    ELSE
        ALLOCATE(A % CPardisoID % ia(nt-nl+2), &
                 A % CPardisoId % ja(nz+nhalo), &
                 A % CPardisoId % aa(nz+nhalo), &
                STAT=allocstat)
    END IF
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso matrix failed')
    END IF
    ia => A % CPardisoId % ia
    ja => A % CPardisoId % ja
    aa => A % CPardisoId % aa

    ! Build distributed CRS matrix
    ia(1) = 1
    lrow = 1      ! Next row to add
    rptr = 1      ! Pointer to next row to add, equals ia(lrow)
    lind = Order(1)-1 ! Row pointer for the first round
    
    ! Add rows of matrix 
    DO i=1,n
      ! Add empty rows until the beginning of the row to add
      ! (first round adds nothing due to choice of lind)
      tind = Order(i)
      rsize = (tind-lind)-1
      
      ! Put zeroes to the diagonal
      DO j=1,rsize
        ia(lrow+j)=rptr+j
        ja(rptr+(j-1))=lind+j
        aa(rptr+(j-1))=0D0
      END DO
      ! Set up row pointers
      rptr = rptr + rsize
      lrow = lrow + rsize
      
      ! Add next row
      rind = iperm(i)
      lind = A % rows(rind)
      tind = A % rows(rind+1)
      IF (ABS(A % CPardisoId % mtype) == 2) THEN
        ! Add only upper triangular elements for symmetric matrices
        rsize = 0
        DO j=lind, tind-1
          IF (A % Gorder(A % Cols(j)) >= Order(i)) THEN
            ja(rptr+rsize)=A % Gorder(A % Cols(j))
            aa(rptr+rsize)=A % values(j)
            rsize = rsize + 1
          END IF
        END DO

      ELSE
        rsize = tind-lind
        DO j=lind, tind-1
          ja(rptr+(j-lind))=A % Gorder(A % Cols(j))
          aa(rptr+(j-lind))=A % values(j)
        END DO
      END IF
        
      ! Sort column indices
      CALL SortF(rsize, ja(rptr:rptr+rsize), aa(rptr:rptr+rsize))
        
      ! Set up row pointers
      rptr = rptr + rsize
      lrow = lrow + 1
      ia(lrow) = rptr

      lind = Order(i) ! Store row index for next round
    END DO

    ! Deallocate temp storage
    DEALLOCATE(Order, iperm)

    ! Set up parameters
    A % CPardisoId % msglvl    = 0 ! Do not write out = 0 / write out = 1 info
    A % CPardisoId % maxfct    = 1 ! Set up space for 1 matrix at most
    A % CPardisoId % mnum      = 1 ! Matrix to use in the solution phase (1st and only one)
    A % CPardisoId % nrhs      = 1 ! Use only one RHS
    ALLOCATE(A % CPardisoId % rhs(nt-nl+1), &
             A % CPardisoId % x(nt-nl+1), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso rhs and solution vector x failed')
    END IF

    ! Set up parameters explicitly
    iparm(1)=1          ! Do not use = 1 / use = 0 solver default parameters
    iparm(2)=2          ! Minimize fill-in with OpenMP nested dissection
    iparm(5)=0          ! No user input permutation
    iparm(6)=0          ! Write solution vector to x
    iparm(8)=0          ! Number of iterative refinement steps
    IF (A % CPardisoID % mtype ==11 .OR. &
        A % CPardisoID % mtype == 13) THEN
      iparm(10)=13      ! Perturbation value 10^-iparm(10) in case of small pivots
      iparm(11)=1       ! Use scalings from symmetric weighted matching
      iparm(13)=1       ! Use permutations from nonsymmetric weighted matching
    ELSE
      iparm(10)=8       ! Perturbation value 10^-iparm(10) in case of small pivots
      iparm(11)=0       ! Do not use scalings from symmetric weighted matching
      iparm(13)=0       ! Do not use permutations from symmetric weighted matching
    END IF
    
    iparm(21)=1         ! Do not use Bunch Kaufman pivoting
    iparm(27)=0         ! Do not check sparse matrix representation
    iparm(28)=0         ! Use double precision
    iparm(35)=0         ! Use Fortran indexing

    ! CPardiso matrix input format
    iparm(40) = 2       ! Distributed solution phase, distributed solution vector
    iparm(41) = nl      ! Beginning of solution domain
    iparm(42) = nt      ! End of solution domain

    ! Perform analysis
    phase = 11      ! Analysis
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype, phase, A % CPardisoId % n, &
          aa, ia, ja, idum, A % CPardisoId % nrhs, iparm, &
          A % CPardisoId % msglvl, &
          ddum, ddum, A % Comm, ierror)
    IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during analysis phase')
    END IF

    ! Perform factorization
    phase = 22      ! Factorization
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype, phase, A % CPardisoId % n, &
          aa, ia, ja, idum, A % CPardisoId % nrhs, iparm, &
          A % CPardisoId % msglvl, &
          ddum, ddum,  A % Comm, ierror)
    IF (ierror .NE. 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during factorization phase')
    END IF
  END SUBROUTINE CPardiso_Factorize


  SUBROUTINE CPardiso_Free(A)
    IMPLICIT NONE

    TYPE(Matrix_t) :: A
    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
                           values, rows, cols, perm, nrhs, iparm, &
                           msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    INTEGER :: ierror, phase, idum(1)
    REAL(kind=dp) :: ddum(1)

    IF(ASSOCIATED(A % CPardisoId)) THEN
        phase     = -1           ! release internal memory
        CALL cluster_sparse_solver(A % CPardisoId % ID, &
              A % CPardisoID % maxfct, A % CPardisoID % mnum, &
              A % CPardisoID % mtype, phase, A % CPardisoID % n, &
              A % CPardisoId % aa, A % CPardisoId % ia, A % CPardisoId % ja, &
              idum, A % CPardisoID % nrhs, A % CPardisoID % IParm, &
              A % CPardisoID % msglvl, ddum, ddum, A % Comm, ierror)

        ! Deallocate global ordering
        DEALLOCATE(A % Gorder)

        ! Deallocate control structures
        DEALLOCATE(A % CPardisoId % ID)
        DEALLOCATE(A % CPardisoID % IParm)
        ! Deallocate matrix and permutation
        DEALLOCATE(A % CPardisoId % ia, A % CPardisoId % ja, &
                   A % CPardisoId % aa, A % CPardisoID % rhs, &
                   A % CPardisoID % x)
        ! Deallocate CPardiso structure
        DEALLOCATE(A % CPardisoID)
    END IF
  END SUBROUTINE CPardiso_Free
#endif


!------------------------------------------------------------------------------
  SUBROUTINE DirectSolver( A,x,b,Solver,Free_Fact )
!------------------------------------------------------------------------------

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: x(*),b(*)
    TYPE(Matrix_t) :: A
    LOGICAL, OPTIONAL :: Free_Fact
!------------------------------------------------------------------------------

    LOGICAL :: GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: Method
!------------------------------------------------------------------------------

    IF ( PRESENT(Free_Fact) ) THEN
      IF ( Free_Fact ) THEN
        CALL BandSolver( A, x, b, Free_Fact )
        CALL ComplexBandSolver( A, x, b, Free_Fact )
#ifdef HAVE_MUMPS
        CALL Mumps_SolveSystem( Solver, A, x, b, Free_Fact )
        CALL MumpsLocal_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#if defined(HAVE_MKL) || defined(HAVE_PARDISO)
        CALL Pardiso_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
        CALL CPardiso_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_SUPERLU
        CALL SuperLU_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_UMFPACK
        CALL Umfpack_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_CHOLMOD
        CALL SPQR_SolveSystem( Solver, A, x, b, Free_Fact )
        CALL Cholmod_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_FETI4I
        CALL Permon_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
        RETURN
        RETURN
      END IF
    END IF

    Method = ListGetString( Solver % Values, 'Linear System Direct Method',GotIt )
    IF ( .NOT. GotIt ) Method = 'banded'
    
    
    CALL Info('DirectSolver','Using direct method: '//TRIM(Method),Level=9)

    SELECT CASE( Method )
      CASE( 'banded', 'symmetric banded' )
        IF ( .NOT. A % Complex ) THEN
           CALL BandSolver( A, x, b )
        ELSE
           CALL ComplexBandSolver( A, x, b )
        END IF

      CASE( 'umfpack', 'big umfpack' )
        CALL Umfpack_SolveSystem( Solver, A, x, b )

      CASE( 'cholmod' )
        CALL Cholmod_SolveSystem( Solver, A, x, b )

      CASE( 'spqr' )
        CALL SPQR_SolveSystem( Solver, A, x, b )

      CASE( 'mumps' )
        CALL Mumps_SolveSystem( Solver, A, x, b )

      CASE( 'mumpslocal' )
        CALL MumpsLocal_SolveSystem( Solver, A, x, b )

      CASE( 'superlu' )
        CALL SuperLU_SolveSystem( Solver, A, x, b )

      CASE( 'permon' )
        CALL Permon_SolveSystem( Solver, A, x, b )

      CASE( 'pardiso' )
        CALL Pardiso_SolveSystem( Solver, A, x, b )

      CASE( 'cpardiso' )
        CALL CPardiso_SolveSystem( Solver, A, x, b )

      CASE DEFAULT
        CALL Fatal( 'DirectSolver', 'Unknown direct solver method.' )
    END SELECT

!------------------------------------------------------------------------------
  END SUBROUTINE DirectSolver
!------------------------------------------------------------------------------

END MODULE DirectSolve

!> \} ElmerLib
