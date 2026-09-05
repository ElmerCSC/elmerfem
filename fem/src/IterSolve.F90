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
! ******************************************************************************
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
! ****************************************************************************/

#include "huti_fdefs.h"

!> \ingroup ElmerLib 
!> \{




!------------------------------------------------------------------------------
!> Module containing control of the iterative solvers for linear systems
!> that come with the Elmer suite.
!------------------------------------------------------------------------------
MODULE IterSolve

!$ USE omp_lib ! conditionally, for the thread ids in the dot products below

   USE Lists
   USE BandMatrix
   USE IterativeMethods
   USE huti_sfe

   IMPLICIT NONE

   !/*
   ! * Iterative method selection
   ! */
   INTEGER, PARAMETER, PRIVATE :: ITER_BiCGStab     =           320
   INTEGER, PARAMETER, PRIVATE :: ITER_TFQMR        =           330
   INTEGER, PARAMETER, PRIVATE :: ITER_CG           =           340
   INTEGER, PARAMETER, PRIVATE :: ITER_CGS          =           350
   INTEGER, PARAMETER, PRIVATE :: ITER_GMRES        =           360
   INTEGER, PARAMETER, PRIVATE :: ITER_BiCGStab2    =           370
   INTEGER, PARAMETER, PRIVATE :: ITER_SGS          =           380
   INTEGER, PARAMETER, PRIVATE :: ITER_JACOBI       =           390
   INTEGER, PARAMETER, PRIVATE :: ITER_RICHARDSON   =           391
   INTEGER, PARAMETER, PRIVATE :: ITER_BICGSTABL    =           400
   INTEGER, PARAMETER, PRIVATE :: ITER_GCR          =           410
   INTEGER, PARAMETER, PRIVATE :: ITER_IDRS         =           420
   INTEGER, PARAMETER, PRIVATE :: ITER_MPRGP        =           430

   !/*
   ! * Preconditioning type code
   ! */
   INTEGER, PARAMETER, PRIVATE :: PRECOND_NONE      =           500
   INTEGER, PARAMETER, PRIVATE :: PRECOND_DIAGONAL  =           510
   INTEGER, PARAMETER, PRIVATE :: PRECOND_ILUn      =           520
   INTEGER, PARAMETER, PRIVATE :: PRECOND_ILUT      =           530
   INTEGER, PARAMETER, PRIVATE :: PRECOND_MG        =           540
   INTEGER, PARAMETER, PRIVATE :: PRECOND_BILUn     =           550
   INTEGER, PARAMETER, PRIVATE :: PRECOND_Vanka     =           560
   INTEGER, PARAMETER, PRIVATE :: PRECOND_Circuit   =           570
   INTEGER, PARAMETER, PRIVATE :: PRECOND_Slave     =           580

   INTEGER, PARAMETER :: stack_max=64
   INTEGER :: stack_pos=0
   LOGICAL :: FirstCall(stack_max)

   REAL(KIND=dp), POINTER :: fm_Diag(:), fm_G(:,:)

CONTAINS


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
#ifndef HUTI_MPRGP_GAMMA
#define HUTI_MPRGP_GAMMA dpar(6)
#endif
#ifndef HUTI_MPRGP_TOLFACTOR
#define HUTI_MPRGP_TOLFACTOR dpar(7)
#endif
!#ifndef HUTI_MPRGP_BOUND
!#define HUTI_MPRGP_BOUND ipar(19)
!#endif
#ifndef HUTI_MPRGP_ADAPT
#define HUTI_MPRGP_ADAPT ipar(19)
#endif

!------------------------------------------------------------------------------
!> Dummy preconditioner, if linear system scaling is active this corresponds
!> to diagonal preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE pcond_dummy(u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM)
    INTEGER :: i
!------------------------------------------------------------------------------
    !$OMP PARALLEL DO
    DO i=1,HUTI_NDIM
       u(i) = v(i)
    END DO
    !$OMP END PARALLEL DO
!------------------------------------------------------------------------------
  END SUBROUTINE pcond_dummy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Complex dummy preconditioner, if linear system scaling is active this corresponds
!> to diagonal preconditioning.
!------------------------------------------------------------------------------
  SUBROUTINE pcond_dummy_cmplx(u,v,ipar )
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(HUTI_NDIM), v(HUTI_NDIM)
!------------------------------------------------------------------------------
    u = v 
!------------------------------------------------------------------------------
  END SUBROUTINE pcond_dummy_cmplx
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE fm_DiagPrec( u,v,ipar )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: u(*),v(*)
    INTEGER :: ipar(*)

    INTEGER :: n

    n = HUTI_NDIM
    u(1:n) = v(1:n)*fm_diag(1:n)
!------------------------------------------------------------------------------
  END SUBROUTINE fm_DiagPrec
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE fm_MatVec( u,v,ipar )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(*),v(*), ct, rsum, cumt=0, s

    INTEGER :: i,j,n

    n = HUTI_NDIM
#if 1
!   CALL DSYMV('U',n,1.0_dp,Jacobian,n,u,1,0.0_dp,v,1)
    CALL DGEMV('N',n,n,1.0_dp,fm_G,n,u,1,0.0_dp,v,1)
#else
    v(1:n) = 0
!$omp parallel do private(i,j,s) shared(n,u,v,fM_G)
    DO i=1,n
       s = 0._dp
       DO j=1,n
         s = s + fm_G(i,j) * u(j)
       END DO
       v(i) = s
    END DO
!$omp end parallel do
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE fm_MatVec
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !> Create mask for skipping edges on a given boundary. 
  !------------------------------------------------------------------------------
  SUBROUTINE CreateEdgeSkipMask(SkipMask)

    LOGICAL, POINTER :: SkipMask(:)
    INTEGER :: t,n0,e0,t0,bc_id,i,j
    LOGICAL :: Found, Piola
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: pVar    

    Mesh => CurrentModel % Mesh

    NULLIFY(pVar)
    DO t=1,CurrentModel % NumberOfSolvers
      IF(ListGetLogical(CurrentModel % Solvers(t) % Values,'Edge Basis',Found ) ) THEN
        pVar => CurrentModel % Solvers(t) % Variable
        EXIT
      END IF
    END DO
    IF(.NOT. ASSOCIATED(pVar)) THEN
      CALL Fatal('CreateEdgeSkipMask','Could not find "Edge Basis" defined in any Solver!')
    END IF

    n0 = Mesh % NumberOfNodes
    t0 = Mesh % NumberOfBulkElements
    e0 = Mesh % NumberOfEdges

    Piola = ListGetLogicalAnySolver( CurrentModel,'Use Piola Transform' ) 

    SkipMask = .FALSE.

    DO t=t0+1,t0+Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)

      IF(.NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      DO bc_id=1,CurrentModel % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE
      BC => CurrentModel % BCs(bc_id) % Values

      IF(ListGetLogical(BC,'Edge Skip Mask',Found ) ) THEN
        SkipMask(pVar % Perm(n0+Element % EdgeIndexes)) = .TRUE.
      END IF
    END DO

    i = COUNT(SkipMask)
    CALL Info('CreateEdgeSkipMask','Mask includes edges on BC: '//I2S(i)//' (out of '//I2S(e0)//')',Level=7)   

    
    ! It is not self-evident that we should include the additional Piola nodes
    ! in the set of nodes to be skipped in smoothing / krylov iteration.
    ! Numerical evidence seems to suggest that this is a good idea. 
    IF(Piola) THEN
      IF(SIZE(pVar % Perm) < n0+e0+2*Mesh % NumberOfFaces) THEN
        CALL Fatal('CreateEdgeSkipMask','Size of Perm too small for Piola!')
      END IF
      
      DO t=1, Mesh % NumberOfFaces
        Element => Mesh % Faces(t)

        ! Only for quads do we have the extra dofs related to Piola transformed edge elements.
        IF(Element % TYPE % ElementCode / 100 == 4 ) THEN
          IF(ALL(SkipMask(pVar % Perm(n0+Element % EdgeIndexes)))) THEN
            DO i=0,1
              j = pVar % Perm(n0+e0+2*t-i)
              IF(j>0) SkipMask(j) = .TRUE.
            END DO
          END IF
        END IF
      END DO
      
      i = COUNT(SkipMask)
      CALL Info('CreateEdgeSkipMask','Mask includes total dofs on BC: '//I2S(i), Level=7)
    END IF
    
  END SUBROUTINE CreateEdgeSkipMask

  
  !------------------------------------------------------------------------------
  !> Create mask for skipping nodes on a given boundary. 
  !------------------------------------------------------------------------------
  SUBROUTINE CreateNodeSkipMask(SkipMask, pVar )

    LOGICAL, POINTER :: SkipMask(:)
    TYPE(Variable_t), POINTER :: pVar    

    INTEGER :: t,n0,e0,t0,bc_id
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh

    IF(.NOT. ListGetLogicalAnyBC(CurrentModel,'Edge Skip Mask' ) ) RETURN
    
    Mesh => CurrentModel % Mesh      
    t0 = Mesh % NumberOfBulkElements
    SkipMask = .FALSE.
    
    DO t=t0+1,t0+Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)

      IF(.NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      DO bc_id=1,CurrentModel % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE
      BC => CurrentModel % BCs(bc_id) % Values

      IF(ListGetLogical(BC,'Edge Skip Mask',Found ) ) THEN
        WHERE(pVar % Perm(Element % NodeIndexes) > 0) 
          SkipMask(pVar % Perm(Element % NodeIndexes)) = .TRUE.
        END WHERE
      END IF
    END DO

    n0 = COUNT(SkipMask)
    CALL Info('CreateNodeSkipMask','Created mask for skipping nodes: '//I2S(n0),Level=7)
    
  END SUBROUTINE CreateNodeSkipMask

  

!> Computed masked dot product.
!----------------------------------------------------------------------
!> The inner products below are the reason a threaded solve is not
!> reproducible run to run unless they are written this way. With
!> REDUCTION(+:s) the per-thread partial sums are combined in the order the
!> threads happen to arrive, so the last bits of every dot product are a
!> function of the scheduling, not of the data. A Krylov method turns that
!> into a different iterate: BiCGStab on ContactBlunt2Djump stopped between
!> 65 and 67 iterations from the same bitwise identical matrix and rhs, and
!> the solution moved by ~1e-6 relative -- far above the 1e-8 linear
!> tolerance, since the contact system is ill conditioned. The contact
!> iteration then took anywhere from 23 to 32 nonlinear steps against a
!> limit of 40, i.e. the test was passing on luck.
!>
!> Accumulating into per-thread slots under SCHEDULE(STATIC) and summing the
!> slots in thread order instead makes the result bitwise reproducible for a
!> given thread count. It does NOT make it agree with the serial sum, and is
!> not meant to: reproducibility across runs is what a regression test needs.
!>
!> END DO NOWAIT matters here. A thread writes only its own slot, and the join
!> at END PARALLEL already orders those writes against the serial sum that
!> follows, so the barrier at END DO buys nothing -- and it is not free: on a
!> 5141 long vector, the length of the ContactBlunt2Djump system, keeping it
!> cost 21% at four threads and 15% at two against the REDUCTION this replaced,
!> while with NOWAIT the same measurement is within 1-2% of it. Longer vectors
!> are memory bound and show no difference either way; at one thread this path
!> skips the parallel region altogether and is ~11% faster than the original.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
FUNCTION MaskedDotProd( ndim, x, xind, y, yind ) RESULT(dres)
!----------------------------------------------------------------------
  IMPLICIT NONE

  ! Parameters
  INTEGER :: ndim, xind, yind
  REAL(KIND=dp) :: x(*)
  REAL(KIND=dp) :: y(*)
  REAL(KIND=dp) :: dres

  ! Local variables
  REAL(KIND=dp) :: s
  INTEGER :: i

  LOGICAL, POINTER, SAVE :: SkipMask(:) => NULL()
  
  IF(.NOT. ASSOCIATED(SkipMask)) THEN
    ALLOCATE(SkipMask(ndim))
    CALL CreateEdgeSkipMask(SkipMask)
  END IF
    
  BLOCK
    REAL(KIND=dp), ALLOCATABLE :: part(:)
    REAL(KIND=dp) :: s
    INTEGER :: nthr, thr
    nthr = 1
!$  nthr = omp_get_max_threads()

    IF( nthr <= 1 ) THEN
      dres = 0
      DO i=1,ndim
      IF(SkipMask(i)) CYCLE
        dres = dres + y(i) * x(i)
      END DO
    ELSE
      ALLOCATE( part(nthr) )
      part = 0
!$OMP PARALLEL PRIVATE(i,thr,s) SHARED(part) NUM_THREADS(nthr)
      thr = 1
!$    thr = omp_get_thread_num() + 1
      s = 0
!$OMP DO SCHEDULE(STATIC)
      DO i=1,ndim
      IF(SkipMask(i)) CYCLE
        s = s + y(i) * x(i)
      END DO
!$OMP END DO NOWAIT
      part(thr) = s
!$OMP END PARALLEL
      dres = 0
      DO i=1,nthr
        dres = dres + part(i)
      END DO
      DEALLOCATE( part )
    END IF
  END BLOCK
!!!CALL SParActiveSUM(dres,0)

!----------------------------------------------------------------------
 END FUNCTION MaskedDotProd
!----------------------------------------------------------------------


!> Compute global 2-norm of vector x
!----------------------------------------------------------------------
FUNCTION MaskedNorm( ndim, x, xind ) RESULT(dres)
!----------------------------------------------------------------------
  IMPLICIT NONE

  ! Parameters

  INTEGER :: ndim, xind
  REAL(KIND=dp) :: x(*)
  REAL(KIND=dp) :: dres

  ! Local variables
  INTEGER :: i
  LOGICAL, POINTER, SAVE :: SkipMask(:) => NULL()
  
  IF(.NOT. ASSOCIATED(SkipMask)) THEN
    ALLOCATE(SkipMask(ndim))
    CALL CreateEdgeSkipMask(SkipMask)
  END IF

  BLOCK
    REAL(KIND=dp), ALLOCATABLE :: part(:)
    REAL(KIND=dp) :: s
    INTEGER :: nthr, thr
    nthr = 1
!$  nthr = omp_get_max_threads()

    IF( nthr <= 1 ) THEN
      dres = 0
      DO i=1,ndim
      IF(SkipMask(i)) CYCLE
        dres = dres + x(i)*x(i)
      END DO
    ELSE
      ALLOCATE( part(nthr) )
      part = 0
!$OMP PARALLEL PRIVATE(i,thr,s) SHARED(part) NUM_THREADS(nthr)
      thr = 1
!$    thr = omp_get_thread_num() + 1
      s = 0
!$OMP DO SCHEDULE(STATIC)
      DO i=1,ndim
      IF(SkipMask(i)) CYCLE
        s = s + x(i)*x(i)
      END DO
!$OMP END DO NOWAIT
      part(thr) = s
!$OMP END PARALLEL
      dres = 0
      DO i=1,nthr
        dres = dres + part(i)
      END DO
      DEALLOCATE( part )
    END IF
  END BLOCK
!!!CALL SParActiveSUM(dres,0)
  dres = SQRT(dres)

!----------------------------------------------------------------------
END FUNCTION MaskedNorm
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  FUNCTION Otmp_ddot( ndim, x, xind, y, yind ) RESULT(dres)
!----------------------------------------------------------------------
    IMPLICIT NONE

    ! Parameters
    INTEGER :: ndim, xind, yind
    REAL(KIND=dp) :: x(*)
    REAL(KIND=dp) :: y(*)
    REAL(KIND=dp) :: dres

    INTEGER :: i

    IF(  xind/=1 .OR. yind /=1 ) THEN
       dres = ddot(ndim,x,xind,y,yind)
       RETURN
    END IF

    BLOCK
    REAL(KIND=dp), ALLOCATABLE :: part(:)
    REAL(KIND=dp) :: s
    INTEGER :: nthr, thr
    nthr = 1
!$  nthr = omp_get_max_threads()

    IF( nthr <= 1 ) THEN
      dres = 0
      DO i=1,ndim
        dres = dres + x(i) * y(i)
      END DO
    ELSE
      ALLOCATE( part(nthr) )
      part = 0
!$OMP PARALLEL PRIVATE(i,thr,s) SHARED(part) NUM_THREADS(nthr)
      thr = 1
!$    thr = omp_get_thread_num() + 1
      s = 0
!$OMP DO SCHEDULE(STATIC)
      DO i=1,ndim
        s = s + x(i) * y(i)
      END DO
!$OMP END DO NOWAIT
      part(thr) = s
!$OMP END PARALLEL
      dres = 0
      DO i=1,nthr
        dres = dres + part(i)
      END DO
      DEALLOCATE( part )
    END IF
    END BLOCK

!----------------------------------------------------------------------
  END FUNCTION Otmp_ddot
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  FUNCTION Otmp_zdotc( ndim, x, xind, y, yind ) RESULT(zres)
!----------------------------------------------------------------------
    IMPLICIT NONE

    ! Parameters
    INTEGER :: ndim, xind, yind
    COMPLEX(KIND=dp) :: x(*)
    COMPLEX(KIND=dp) :: y(*)
    COMPLEX(KIND=dp) :: zres

    INTEGER :: i

    IF(  xind/=1 .OR. yind /=1 ) THEN
       zres = zdotc(ndim,x,xind,y,yind)
       RETURN
    END IF

    BLOCK
    COMPLEX(KIND=dp), ALLOCATABLE :: part(:)
    COMPLEX(KIND=dp) :: s
    INTEGER :: nthr, thr
    nthr = 1
!$  nthr = omp_get_max_threads()

    IF( nthr <= 1 ) THEN
      zres = 0
      DO i=1,ndim
        zres = zres + DCONJG(x(i)) * y(i)
      END DO
    ELSE
      ALLOCATE( part(nthr) )
      part = 0
!$OMP PARALLEL PRIVATE(i,thr,s) SHARED(part) NUM_THREADS(nthr)
      thr = 1
!$    thr = omp_get_thread_num() + 1
      s = 0
!$OMP DO SCHEDULE(STATIC)
      DO i=1,ndim
        s = s + DCONJG(x(i)) * y(i)
      END DO
!$OMP END DO NOWAIT
      part(thr) = s
!$OMP END PARALLEL
      zres = 0
      DO i=1,nthr
        zres = zres + part(i)
      END DO
      DEALLOCATE( part )
    END IF
    END BLOCK

!----------------------------------------------------------------------
  END FUNCTION Otmp_zdotc
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> As Otmp_zdotc but WITHOUT conjugation, i.e. the bilinear form x^T y.
!> CG is the one Krylov method whose correctness requires the operator to be
!> self-adjoint in the inner product being used, and the complex systems
!> assembled by Elmer are complex symmetric (A = A^T) rather than Hermitian.
!> Methods that only need an inner product for orthogonality or residual
!> minimization -- BiCGStab, GMRES, GCR, TFQMR, CGS -- keep using Otmp_zdotc.
!> NOTE: switching CG to this form is necessary but on its own has not been
!> shown to be sufficient; see the handover notes.
!----------------------------------------------------------------------
  FUNCTION Otmp_zdotu( ndim, x, xind, y, yind ) RESULT(zres)
!----------------------------------------------------------------------
    IMPLICIT NONE

    ! Parameters
    INTEGER :: ndim, xind, yind
    COMPLEX(KIND=dp) :: x(*)
    COMPLEX(KIND=dp) :: y(*)
    COMPLEX(KIND=dp) :: zres

    INTEGER :: i

    IF(  xind/=1 .OR. yind /=1 ) THEN
       zres = zdotu(ndim,x,xind,y,yind)
       RETURN
    END IF

    BLOCK
    COMPLEX(KIND=dp), ALLOCATABLE :: part(:)
    COMPLEX(KIND=dp) :: s
    INTEGER :: nthr, thr
    nthr = 1
!$  nthr = omp_get_max_threads()

    IF( nthr <= 1 ) THEN
      zres = 0
      DO i=1,ndim
        zres = zres + x(i) * y(i)
      END DO
    ELSE
      ALLOCATE( part(nthr) )
      part = 0
!$OMP PARALLEL PRIVATE(i,thr,s) SHARED(part) NUM_THREADS(nthr)
      thr = 1
!$    thr = omp_get_thread_num() + 1
      s = 0
!$OMP DO SCHEDULE(STATIC)
      DO i=1,ndim
        s = s + x(i) * y(i)
      END DO
!$OMP END DO NOWAIT
      part(thr) = s
!$OMP END PARALLEL
      zres = 0
      DO i=1,nthr
        zres = zres + part(i)
      END DO
      DEALLOCATE( part )
    END IF
    END BLOCK

!----------------------------------------------------------------------
  END FUNCTION Otmp_zdotu
!----------------------------------------------------------------------

    
!------------------------------------------------------------------------------
!> The routine that decides which linear system solver to call, and calls it.
!> There are two main sources of iterations within Elmer.
!> 1) The old HUTiter library that includes the most classic iterative Krylov
!>    methods.
!> 2) The internal MODULE IterativeMethods that includes some classic iterative
!>    methods and also some more recent Krylov methods. 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!> May the scalar 2N values of a complex matrix be released for the duration of a
!> solve? Factored out of IterSolver so the test reads as one thing and each
!> clause can carry the reason it is there.
!>
!> IT IS A WHITELIST, and deliberately so: the consumers reach the matrix through
!> GlobalMatrix and cannot be enumerated by grep, so this admits only what has
!> been checked by experiment rather than excluding what happened to be noticed.
!> "Linear System Poison Scalar Values" is that experiment -- it fills the array
!> with 1e300 over the same window and restores the exact bits after -- and every
!> clause below traces to something it caught or confirmed.
!------------------------------------------------------------------------------
  FUNCTION FreeEligible( A, Params, PCondType, StopcProc, BlockCRS, &
      PoisonVals, MatvecF, MatvecReadsNoValues ) RESULT( OK )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: PCondType
    INTEGER(KIND=AddrInt) :: StopcProc
    LOGICAL :: BlockCRS, PoisonVals, OK
    INTEGER(KIND=AddrInt), OPTIONAL :: MatvecF
    LOGICAL, OPTIONAL :: MatvecReadsNoValues
    LOGICAL :: Found
!------------------------------------------------------------------------------
    OK = BlockCRS .AND. .NOT. PoisonVals
    IF( OK .AND. .NOT. ASSOCIATED( A % Values ) ) OK = .FALSE.

    ! A caller-supplied product may read anything, so it has to say otherwise.
    IF( OK .AND. PRESENT( MatvecF ) ) THEN
      OK = .FALSE.
      IF( PRESENT( MatvecReadsNoValues ) ) OK = MatvecReadsNoValues
    END IF

    ! The backward-error criteria read the coefficients directly, e.g.
    ! SUM(GlobalMatrix % Values**2) in BackwardError.F90.
    IF( OK ) OK = ( StopcProc == 0 )

    ! Every preconditioner here reads the block view or the compact factors.
    ! Anything else -- multigrid, Vanka, block, circuit, the auxiliary space
    ! solver -- has not been probed and stays out.
    IF( OK ) OK = ( PCondType == PRECOND_NONE .OR. &
        PCondType == PRECOND_ILUn .OR. PCondType == PRECOND_ILUT .OR. &
        ( PCondType == PRECOND_DIAGONAL .AND. ASSOCIATED( A % BDiag ) ) )

    ! ParallelUtils used to point InsideMatrix % Values at MassValues, and
    ! DEALLOCATE on a pointer that was pointer-assigned is undefined. That idiom
    ! is gone, but ASSOCIATED(a,b) is the one thing testable here, so test it.
    IF( OK ) OK = &
        .NOT. ASSOCIATED( A % Values, A % MassValues ) .AND. &
        .NOT. ASSOCIATED( A % Values, A % DampValues ) .AND. &
        .NOT. ASSOCIATED( A % Values, A % PrecValues ) .AND. &
        .NOT. ASSOCIATED( A % Values, A % BulkValues ) .AND. &
        .NOT. ASSOCIATED( A % Values, A % ILUValues )

    IF( OK ) OK = ListGetLogical( Params, &
        'Linear System Free Scalar Values', Found, DefValue = .TRUE. )
!------------------------------------------------------------------------------
  END FUNCTION FreeEligible
!------------------------------------------------------------------------------


  RECURSIVE SUBROUTINE IterSolver( A,x,b,Solver,ndim,DotF, &
              NormF,MatvecF,PrecF,StopcF,MatvecReadsNoValues,DotFU )
!------------------------------------------------------------------------------
    USE huti_sfe
    USE ListMatrix
    USE SParIterGlobals
    USE GeneralUtils, ONLY : ComplexValues
    IMPLICIT NONE

!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), DIMENSION(:), TARGET CONTIG :: x,b
    TYPE(Matrix_t), TARGET :: A
    !> Set by a caller supplying MatvecF to assert that its product does not
    !> read A % Values. Only the caller can know this, and it decides whether
    !> the scalar values may be released across the iteration.
    LOGICAL, OPTIONAL :: MatvecReadsNoValues
    INTEGER, OPTIONAL :: ndim
    INTEGER(KIND=AddrInt), OPTIONAL :: DotF, NormF, MatVecF, PrecF, StopcF
    !> The unconjugated (bilinear) counterpart of DotF, for a caller that
    !> supplies its own complex inner product. CG needs this form and every
    !> other complex method here needs DotF; a caller cannot pick between them
    !> without duplicating the keyword parsing that decides IterType below, so
    !> it hands over both and the choice is made here. Ignored for real
    !> systems and for any method other than CG.
    INTEGER(KIND=AddrInt), OPTIONAL :: DotFU
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: Adiag,CM,PrecMat,SaveGlobalM

    REAL(KIND=dp) :: dpar(HUTI_DPAR_DFLTSIZE),stopfun
!   external stopfun
    REAL(KIND=dp), ALLOCATABLE :: work(:,:)
    INTEGER :: i,j,k,N,ipar(HUTI_IPAR_DFLTSIZE),wsize,istat,IterType,PCondType,ILUn,Blocks
    LOGICAL :: PoisonVals, FreeVals, ValsFreed
    INTEGER :: nScalarVals
    REAL(KIND=dp), ALLOCATABLE :: SaveVals(:)
    LOGICAL :: Internal, NullEdges
    LOGICAL :: ComponentwiseStopC, NormwiseStopC, RowEquilibration
    LOGICAL :: Condition,GotIt, Refactorize,Found,GotDiagFactor,Robust
    LOGICAL :: ComplexSystem, PseudoComplexSystem, DoFatal, LeftOriented
    LOGICAL :: BlockCRS
    
    REAL(KIND=dp) :: ILUT_TOL, DiagFactor

    TYPE(ValueList_t), POINTER :: Params

    CHARACTER(:), ALLOCATABLE :: str

    EXTERNAL MultigridPrec
    EXTERNAL NormwiseBackwardError, ComponentwiseBackwardError
    EXTERNAL NormwiseBackwardErrorGeneralized
    EXTERNAL NormwiseBackwardError_Z
    
    INTEGER(KIND=Addrint) :: dotProc, normProc, pcondProc, &
        pconddProc, mvProc, iterProc, StopcProc
    INTEGER(KIND=Addrint) :: AddrFunc
    COMPLEX(KIND=dp), POINTER :: xC(:), bC(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: workC(:,:)
    EXTERNAL :: AddrFunc    

    INTERFACE
      SUBROUTINE VankaCreate(A,Solver)
        USE Types
        TYPE(Matrix_t) :: A
        TYPE(Solver_t) :: Solver
      END SUBROUTINE VankaCreate
      
      SUBROUTINE VankaPrec(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        REAL(KIND=dp) :: u(*),v(*)
      END SUBROUTINE VankaPrec

      SUBROUTINE CircuitPrec(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        REAL(KIND=dp) :: u(*),v(*)
      END SUBROUTINE CircuitPrec

      SUBROUTINE CircuitPrecComplex(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        COMPLEX(KIND=dp) :: u(*),v(*)
      END SUBROUTINE CircuitPrecComplex

      SUBROUTINE SlavePrec(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        REAL(KIND=dp) :: u(*),v(*)
      END SUBROUTINE SlavePrec

      SUBROUTINE SlavePrecComplex(u,v,ipar)
        USE Types
        INTEGER :: ipar(*)
        COMPLEX(KIND=dp) :: u(*),v(*)
      END SUBROUTINE SlavePrecComplex      
    END INTERFACE
!------------------------------------------------------------------------------
    N = A % NumberOfRows
    IF ( PRESENT(ndim) ) n=ndim
    
    ipar = 0
    dpar = 0.0D0
    pconddProc = 0
!------------------------------------------------------------------------------
    Params => Solver % Values
    str = ListGetString( Params,'Linear System Iterative Method',Found )
    IF( .NOT. Found ) THEN
      CALL Warn('IterSolver','> Linear System Iterative Method < not found, using BiCGstab')
      str = 'bicgstab'      
    ELSE
      CALL Info('IterSolver','Using iterative method: '//TRIM(str),Level=9)
    END IF
    
    IF( ListGetLogical( Params,'Linear System Skip Complex',GotIt ) ) THEN
      CALL Info('IterSolver','This time skipping complex treatment',Level=20)
      A % COMPLEX = .FALSE.
      ComplexSystem = .FALSE.
    ELSE
      ComplexSystem = ListGetLogical( Params,'Linear System Complex',Found ) 
      IF( .NOT. Found ) ComplexSystem = A % COMPLEX 
    END IF
    
    PseudoComplexSystem = ListGetLogical( Params,'Linear System Pseudo Complex',Found ) 

    IF( ComplexSystem ) THEN
      CALL Info('IterSolver','Matrix is complex valued',Level=10)
    ELSE IF( PseudoComplexSystem ) THEN
      CALL Info('IterSolver','Matrix is pseudo complex valued',Level=10)
    ELSE    
      CALL Info('IterSolver','Matrix is real valued',Level=12)
    END IF

    
    SELECT CASE(str)
    CASE('bicgstab2')
      ! NOTE:
      ! BiCGStab2 should be nearly the same as BiCGStabl with the parameter l=2, but
      ! the implementation of BiCGStabl uses the right-oriented preconditioning, while
      ! BiCGStab2 works as expected only with the left-oriented preconditioning. Due to
      ! the difference in the preconditioning the convergence histories may be quite different.
      ! The complex BiCGStab2 does not convince, use BiCGStabl instead:
      IF (ComplexSystem ) THEN
        IterType = ITER_BICGstabl
      ELSE
        IterType = ITER_BiCGStab2
      END IF
    CASE('bicgstabl')
      IterType = ITER_BICGstabl
    CASE('bicgstab')
      IterType = ITER_BiCGStab
    CASE('tfqmr')
      IterType = ITER_TFQMR
    CASE('cgs')
      IterType = ITER_CGS
    CASE('cg')
      IterType = ITER_CG
    CASE('gmres')
      IterType = ITER_GMRES
    CASE('sgs')
      IterType = ITER_SGS
    CASE('jacobi')
      IterType = ITER_jacobi
    CASE('richardson')
      IterType = ITER_richardson
    CASE('gcr')
      IterType = ITER_GCR
    CASE('idrs')
      IterType = ITER_IDRS
    CASE DEFAULT
      IterType = ITER_BiCGStab
    CASE('mprgp')
      IterType = ITER_MPRGP
    END SELECT
    
!------------------------------------------------------------------------------

    HUTI_WRKDIM = 0
    HUTI_PSEUDOCOMPLEX = 0
    IF( PseudoComplexSystem ) THEN
      HUTI_PSEUDOCOMPLEX = 1     
      IF ( ListGetLogical( Params,'Block Split Complex',Found ) ) HUTI_PSEUDOCOMPLEX = 2
    END IF
    Internal = .FALSE.
    
    SELECT CASE ( IterType )
      
      ! Solvers from HUTiter
      !-------------------------------------------------------       
    CASE (ITER_BiCGStab)
      HUTI_WRKDIM = HUTI_BICGSTAB_WORKSIZE
      
    CASE (ITER_BiCGStab2)
      HUTI_WRKDIM = HUTI_BICGSTAB_2_WORKSIZE
      
    CASE (ITER_TFQMR)
      HUTI_WRKDIM = HUTI_TFQMR_WORKSIZE
      
    CASE (ITER_CG)
      HUTI_WRKDIM = HUTI_CG_WORKSIZE
      
    CASE (ITER_CGS)
      HUTI_WRKDIM = HUTI_CGS_WORKSIZE
      
    CASE (ITER_GMRES)
      HUTI_GMRES_RESTART = ListGetInteger( Params,&
          'Linear System GMRES Restart',  GotIt ) 
      IF ( .NOT. GotIT ) HUTI_GMRES_RESTART = 10
      HUTI_WRKDIM = HUTI_GMRES_WORKSIZE + HUTI_GMRES_RESTART
      
      ! Solvers from IterativeMethods.src
      !-------------------------------------------------------       
    CASE (ITER_SGS)
      HUTI_WRKDIM = 1
      HUTI_SGSPARAM = ListGetConstReal( Params,'SGS Overrelaxation Factor',&
          GotIt,minv=0.0_dp,maxv=2.0_dp)
      IF(.NOT. GotIt) HUTI_SGSPARAM = 1.8_dp
      Internal = .TRUE.
      
    CASE (ITER_Jacobi, ITER_Richardson)
      HUTI_WRKDIM = 1
      Internal = .TRUE.
      
    CASE (ITER_GCR)
      HUTI_WRKDIM = 1
      HUTI_GCR_RESTART = ListGetInteger( Params, &
          'Linear System GCR Restart',  GotIt ) 
      IF ( .NOT. GotIt ) THEN
        i = ListGetInteger( Params,'Linear System Max Iterations', minv=1 )
        IF( i > 200 ) THEN
          i = 200
          CALL Info('IterSolver','"Linear System GCR Restart" not given, setting it to '//I2S(i),Level=4)
        END IF
        HUTI_GCR_RESTART = i
      END IF
      Internal = .TRUE.
      
    CASE (ITER_BICGSTABL)
      HUTI_WRKDIM = 1
      HUTI_BICGSTABL_L = ListGetInteger( Params,'BiCGstabl polynomial degree',&
          GotIt,minv=2)
      IF(.NOT. GotIt) HUTI_BICGSTABL_L = 2
      Internal = .TRUE.

    CASE (ITER_IDRS)
      HUTI_WRKDIM = 1
      HUTI_IDRS_S = ListGetInteger( Params,'IDRS parameter',GotIt,minv=1)
      IF(.NOT. GotIt) HUTI_IDRS_S = 4
      Internal = .TRUE.
    
    CASE (ITER_MPRGP)
      HUTI_WRKDIM = 1
      Internal = .TRUE.
      HUTI_MPRGP_GAMMA = ListGetConstReal( Params, 'Linear System MPRGP Gamma', GotIt )
      IF(.NOT. GotIt) HUTI_MPRGP_GAMMA = 1.0_dp
      HUTI_MPRGP_TOLFACTOR = ListGetConstReal( Params, 'Linear System MPRGP TolFactor', GotIt )
      IF(.NOT. GotIt) HUTI_MPRGP_TOLFACTOR = 5.0_dp
      !HUTI_MPRGP_BOUND = ListGetString( Params, 'Linear System MPRGP Bound Type', GotIt )
      !IF(.NOT. GotIt) HUTI_MPRGP_BOUND = 'lower' ! TODO: should write error if no bounds      
      HUTI_MPRGP_ADAPT = 1
      IF( ListGetLogical( Params, 'Linear System MPRGP Adaptive', GotIt ) ) THEN
        HUTI_MPRGP_ADAPT = 1
      ELSE
        IF(GotIt) HUTI_MPRGP_ADAPT = 0
      END IF
        
    END SELECT
!------------------------------------------------------------------------------
    
    wsize = HUTI_WRKDIM
    
    StopcProc = 0
    IF (PRESENT(StopcF)) THEN
       StopcProc = StopcF
       HUTI_STOPC = HUTI_USUPPLIED_STOPC
    ELSE
       ComponentwiseStopC = ListGetLogical(Params,'Linear System Componentwise Backward Error',GotIt)
       IF (ComponentwiseStopC) THEN
          IF (ComplexSystem) THEN
            CALL Info('IterSolver', 'Linear System Componentwise Backward Error is active')
            CALL Fatal('IterSolver', 'Error computation does not support Linear System Complex = True')
          END IF
          StopcProc = AddrFunc(ComponentwiseBackwardError)
          HUTI_STOPC = HUTI_USUPPLIED_STOPC
       ELSE
          NormwiseStopC = ListGetLogical(Params,'Linear System Normwise Backward Error',GotIt)
          IF (NormwiseStopC) THEN
             RowEquilibration = ListGetLogical(Params,'Linear System Row Equilibration',GotIt)
             IF (RowEquilibration) THEN
               IF (ComplexSystem) THEN
                 StopcProc = AddrFunc(NormwiseBackwardError_Z)
               ELSE
                 StopcProc = AddrFunc(NormwiseBackwardError)
               END IF
             ELSE
               IF (ComplexSystem) THEN
                 CALL Info('IterSolver', 'Linear System Normwise Backward Error is active')
                 CALL Fatal('IterSolver', 'Error computation needs Linear System Row Equilibration = True')
               ELSE
                 StopcProc = AddrFunc(NormwiseBackwardErrorGeneralized)
               END IF
             END IF
             HUTI_STOPC = HUTI_USUPPLIED_STOPC
          ELSE
             HUTI_STOPC = HUTI_TRESID_SCALED_BYB
          END IF
       END IF
    END IF
    HUTI_NDIM  = N
    
    HUTI_DBUGLVL  = ListGetInteger( Params, &
        'Linear System Residual Output', GotIt )
    IF ( .NOT.Gotit ) HUTI_DBUGLVL = 1
    
    IF ( Parenv % myPE /= 0 ) HUTI_DBUGLVL=0
    
    HUTI_MAXIT = ListGetInteger( Params, &
        'Linear System Max Iterations', minv=1 )
    
    HUTI_MINIT = ListGetInteger( Params, &
        'Linear System Min Iterations', GotIt )
    
    IF( ComplexSystem ) THEN
      ALLOCATE(workC(N/2,wsize), stat=istat)
      IF ( istat /= 0 ) THEN
        CALL Fatal( 'IterSolve', 'Memory allocation failure.' )
      END IF
      workC = cmplx(0,0,dp)
    ELSE
      ALLOCATE(work(N,wsize), stat=istat)
      IF ( istat /= 0 ) THEN
        CALL Fatal( 'IterSolve', 'Memory allocation failure.' )
      END IF
      !$OMP PARALLEL PRIVATE(j)
      DO j=1,wsize
        !$OMP DO
        DO i=1,N
          work(i,j) = real(0,dp)
        END DO
        !$OMP END DO
      END DO
      !$OMP END PARALLEL
    END IF

    IF ( (IterType == ITER_BiCGStab2 .OR. IterType == ITER_BiCGStabL .OR. &
         IterType == ITER_BiCGStab ) .AND. ALL(x == 0.0) ) x = 1.0d-8

    HUTI_INITIALX = HUTI_USERSUPPLIEDX
    
    HUTI_TOLERANCE = ListGetCReal( Params, &
        'Linear System Convergence Tolerance' )
    
    HUTI_MAXTOLERANCE = ListGetCReal( Params, &
        'Linear System Divergence Limit', GotIt)
    IF(.NOT. GotIt) HUTI_MAXTOLERANCE = 1.0d20
    
    IF( ListGetLogical( Params,'Linear System Robust',GotIt) ) THEN
      HUTI_ROBUST = 1
      HUTI_ROBUST_TOLERANCE = ListGetCReal( Params,'Linear System Robust Tolerance',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_TOLERANCE = HUTI_TOLERANCE**(2.0/3.0)
      HUTI_ROBUST_MAXTOLERANCE = ListGetCReal( Params,'Linear System Robust Limit',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_MAXTOLERANCE = SQRT( HUTI_TOLERANCE )      
      HUTI_ROBUST_STEPSIZE = ListGetCReal( Params,'Linear System Robust Margin',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_STEPSIZE = 1.1_dp
      HUTI_ROBUST_MAXBADIT = ListGetInteger( Params,'Linear System Robust Max Iterations',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_MAXBADIT = HUTI_MAXIT / 2
      HUTI_ROBUST_START = ListGetInteger( Params,'Linear System Robust Start Iteration',GotIt)
      IF(.NOT. GotIt ) HUTI_ROBUST_START = 1
    ELSE
      HUTI_ROBUST = 0
    END IF


    IF( ListGetLogical( Params,'IDRS Smoothing',GotIt) ) THEN
      HUTI_SMOOTHING = 1
    ELSE
      HUTI_SMOOTHING = 0
    END IF
      
    
!------------------------------------------------------------------------------

    ! By default the right-oriented preconditioning is applied, but BiCGStab2,
    ! GMRES and TFQMR are called with the left-oriented preconditining since 
    ! the right-oriented preconditioning does not work as expected. The methods
    ! from the module IterativeMethods use always the right-oriented preconditioning:
    !
    SELECT CASE ( IterType )
    CASE (ITER_BiCGStab2, ITER_GMRES, ITER_TFQMR)
      LeftOriented = .TRUE.
    CASE DEFAULT
      LeftOriented = ListGetLogical(Params, 'Linear System Left Preconditioning', GotIt)
      IF (Internal) LeftOriented = .FALSE.
    END SELECT

    
    ! Build the block view here rather than at the matvec selection further
    ! down: the preconditioner is set up in between, and the complex ILU reads
    ! the view when it is present. Refreshing it afterwards would have the
    ! factorization see the values of the previous solve.
    ! In parallel A is the SplittedMatrix InsideMatrix and MatvecF is present,
    ! so the product below stays SParCMatrixVector -- but that routine hands the
    ! local part straight to CRS_ComplexMatrixVectorMultiply, and the local ILU
    ! factorizes this same matrix, so both want the view built here just as in
    ! serial.
    ! On by default. The view costs roughly 40% of the matrix on top of the
    ! scalar form it is derived from, so "Linear System Block CRS = False"
    ! turns it off where memory is tighter than time.
    BlockCRS = ComplexSystem
    IF( BlockCRS ) BlockCRS = ListGetLogical( Params,'Linear System Block CRS', &
        Found, DefValue = .TRUE. )
    IF( BlockCRS ) CALL CRS_BuildBlockCRS( A )

    ! PROBE. Poison the scalar values from here -- the moment the view exists --
    ! through to the end of the iteration, and restore the exact bits after. The
    ! window deliberately covers the PRECONDITIONER SETUP as well as the Krylov
    ! loop, because that is where the next increment of peak memory is: the
    ! factors are allocated while the matrix is still resident. Anything in
    ! either phase that still reads A % Values turns the solve into Inf/NaN and
    ! says so, which is how the release's eligibility list was arrived at and is
    ! what any widening of it has to be justified with.
    ValsFreed = .FALSE.
    PoisonVals = .FALSE.
    IF( BlockCRS ) PoisonVals = ListGetLogical( Params, &
        'Linear System Poison Scalar Values', Found )
    IF( PoisonVals ) THEN
      ALLOCATE( SaveVals(SIZE(A % Values)) )
      SaveVals = A % Values
      A % Values = 1.0e300_dp
      CALL Info('IterSolver', &
          'PROBE: scalar values poisoned from the view onwards',Level=5)
    END IF

    IF ( .NOT. PRESENT(PrecF) ) THEN
      str = ListGetString( Params, 'Linear System Preconditioning',gotit )
      IF ( .NOT.gotit ) str = 'none'
      
      A % Cholesky = ListGetLogical( Params,'Linear System Symmetric ILU', Gotit )
      
      ILUn = -1
      IF ( str == 'none' ) THEN
        PCondType = PRECOND_NONE

      ELSE IF ( str == 'diagonal' ) THEN
        PCondType = PRECOND_DIAGONAL

      ELSE IF ( str == 'ilut' ) THEN
        ILUT_TOL = ListGetCReal( Params, &
            'Linear System ILUT Tolerance',GotIt )
        PCondType = PRECOND_ILUT

      ELSE IF ( SEQL(str, 'ilu') ) THEN
        ILUn = ListGetInteger( Params, &
            'Linear System ILU Order', gotit )
        IF ( .NOT.gotit ) THEN
          IF(LEN(str)>=4) ILUn = ICHAR(str(4:4)) - ICHAR('0')
        END IF
        IF ( ILUn  < 0 .OR. ILUn > 9 ) ILUn = 0
        PCondType = PRECOND_ILUn

      ELSE IF ( SEQL(str, 'bilu') ) THEN
        ILUn = 0
        IF(LEN(str)>=5) ILUn = ICHAR(str(5:5)) - ICHAR('0')
        IF ( ILUn  < 0 .OR. ILUn > 9 ) ILUn = 0
        IF( Solver % Variable % Dofs == 1) THEN
          CALL Warn('IterSolver','BILU for one dofs is equal to ILU!')
          PCondType = PRECOND_ILUn
        ELSE
          PCondType = PRECOND_BILUn
        END IF

      ELSE IF ( str == 'multigrid' ) THEN
        PCondType = PRECOND_MG

      ELSE IF ( SEQL(str,'vanka') ) THEN
        PCondType = PRECOND_VANKA
        
      ELSE IF ( str == 'auxiliary space solver' .OR. str == 'slave' ) THEN
        PCondType = PRECOND_SLAVE
        
      ELSE IF ( str == 'circuit' ) THEN
        ILUn = ListGetInteger( Params, 'Linear System ILU Order', gotit )
        IF(.NOT.Gotit ) ILUn=-1
        PCondType = PRECOND_Circuit

      ELSE
        PCondType = PRECOND_NONE
        CALL Warn( 'IterSolve', 'Unknown preconditioner type: '//TRIM(str)//', feature disabled.' )
      END IF

      IF ( .NOT. ListGetLogical( Params, 'No Precondition Recompute',GotIt ) ) THEN
        CALL ResetTimer("Prec-"//TRIM(str))

        n = ListGetInteger( Params, 'Linear System Precondition Recompute', GotIt )
        IF ( n <= 0 ) n = 1
        
        Refactorize = ListGetLogical( Params, 'Linear System Refactorize', Gotit )
        IF ( .NOT. Gotit ) Refactorize = .TRUE.
        
        IF (.NOT.(ASSOCIATED(A % ILUValues).OR.ASSOCIATED(A % CILUValues)).OR. &
                  (Refactorize.AND.MOD(A % SolveCount, n)==0) ) THEN


          IF ( A % FORMAT == MATRIX_CRS ) THEN

            ! Optionally one may emphasize the diagonal entries in the linear system
            ! to make the preconditioning more stable.
            !-------------------------------------------------------------------------          
            DiagFactor = ListGetCReal( Params,'Linear System ILU Factor',GotIt ) 
            GotDiagFactor = ( DiagFactor > EPSILON( DiagFactor ) ) 
            IF( GotDiagFactor ) THEN
              CALL Info('IterSolver','Applying diagonal relaxation for ILU', Level=8)
              DiagFactor = DiagFactor + 1.0_dp
              A % Values( A % Diag ) = DiagFactor * A % Values( A % Diag )      
            END IF

            IF ( ComplexSystem ) THEN
              IF ( PCondType == PRECOND_ILUn .OR. (PCondType==PRECOND_Circuit.AND.ILUn>=0) ) THEN
                NullEdges = ListGetLogical(Params, 'Edge Basis', GotIt)
                CM => A % ConstraintMatrix
                IF(NullEdges.OR.ASSOCIATED(CM)) THEN
                  CALL Info('IterSolver','Omitting edge dofs from being target of ILUn',Level=20)

                  IF(ASSOCIATED(A % ILURows)) DEALLOCATE(A % ILURows)
                  IF(ASSOCIATED(A % ILUCols)) DEALLOCATE(A % ILUCols)
                  IF(ASSOCIATED(A % ILUDiag)) DEALLOCATE(A % ILUDiag)
                  IF(ASSOCIATED(A % CILUValues)) DEALLOCATE(A % CILUValues)

                  PrecMat => AllocateMatrix()
                  PrecMat % FORMAT = MATRIX_LIST
                  PrecMat % CIluValues => NULL()

                  IF(ASSOCIATED(CM)) THEN
                    DO i=CM % NumberOfRows,1,-1
                      k = i + A % NumberOfRows
                      CALL List_AddMatrixIndex( PrecMat % ListMatrix,k,k)
                      IF(MOD(k,2)==0) THEN
                        CALL List_AddMatrixIndex(PrecMat % ListMatrix, k, k-1)
                      ELSE
                        CALL List_AddMatrixIndex(PrecMat % ListMatrix, k, k+1)
                      END IF

                      DO j=CM % Rows(i+1)-1,CM % Rows(i),-1
                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, CM % Cols(j), CM % Values(j))

                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             CM % Cols(j), i + A % NumberOfRows, CM % Values(j))
                      END DO
                    END DO
                  END IF

                  k = A % NumberOfRows - A % ExtraDOFs
                  DO i=A % NumberOfRows,1,-1
                    IF(i>k) THEN
                       CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i)
                       IF(MOD(i,2)==0) THEN
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i-1)
                       ELSE
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i+1)
                       END IF
                    ELSE IF (NullEdges) THEN
                       CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, i, 1._dp)
                       IF(MOD(i,2)==0) THEN
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i-1)
                       ELSE
                         CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i+1)
                       END IF
                    END IF

                    DO j=A % Rows(i+1)-1,A % Rows(i),-1
                      IF (i>k .OR. A % Cols(j)>k .OR. .NOT.NullEdges) THEN
                        CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, A % Cols(j), A % Values(j))
                      END IF
                    END DO
                  END DO

                  CALL List_ToCRSMatrix(PrecMat)
                  Condition = CRS_ComplexIncompleteLU(PrecMat,ILUn)

                  A % ILURows => PrecMat % IluRows
                  A % ILUCols => PrecMat % IluCols
                  A % ILUDiag => PrecMat % IluDiag
                  A % CILUvalues => PrecMat % CIluValues

                  DEALLOCATE(PrecMat % Values)
                  IF(.NOT.ASSOCIATED(A % ILURows,PrecMat % Rows)) DEALLOCATE(PrecMat % Rows)
                  IF(.NOT.ASSOCIATED(A % ILUCols,PrecMat % Cols)) DEALLOCATE(PrecMat % Cols)
                  IF(.NOT.ASSOCIATED(A % ILUDiag,PrecMat % Diag)) DEALLOCATE(PrecMat % Diag)
                  DEALLOCATE(PrecMat)
                ELSE
                  Condition = CRS_ComplexIncompleteLU(A,ILUn)
                END IF
              ELSE IF ( PCondType == PRECOND_ILUT ) THEN
                Condition = CRS_ComplexILUT( A,ILUT_TOL )
              END IF
            ELSE IF (ILUn>=0 .OR. PCondType == PRECOND_ILUT) THEN  ! Not ComplexSystem
              SELECT CASE(PCondType)
              CASE(PRECOND_ILUn, PRECOND_Circuit)
                NullEdges = ListGetLogical(Params, 'Edge Basis', GotIt)
                CM => A % ConstraintMatrix
                IF(NullEdges.OR.ASSOCIATED(CM)) THEN

                  IF(ASSOCIATED(A % ILURows)) DEALLOCATE(A % ILURows)
                  IF(ASSOCIATED(A % ILUCols)) DEALLOCATE(A % ILUCols)
                  IF(ASSOCIATED(A % ILUDiag)) DEALLOCATE(A % ILUDiag)
                  IF(ASSOCIATED(A % ILUValues)) DEALLOCATE(A % ILUValues)

                  PrecMat => AllocateMatrix()
                  PrecMat % FORMAT = MATRIX_LIST

                  IF(ASSOCIATED(CM)) THEN
                    DO i=CM % NumberOfRows,1,-1
                      CALL List_AddMatrixIndex( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, i + A % NumberOFrows)

                      DO j=CM % Rows(i+1)-1,CM % Rows(i),-1
                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             i + A % NumberOfRows, CM % Cols(j), CM % Values(j))

                        CALL List_AddToMatrixElement( PrecMat % ListMatrix, &
                             CM % Cols(j), i + A % NumberOfRows, CM % Values(j))
                      END DO
                    END DO
                  END IF

                  k = A % NumberOfRows - A % ExtraDOFs
                  DO i=A % NumberOfRows,1,-1
                    IF(i>k) THEN
                       CALL List_AddMatrixIndex(PrecMat % ListMatrix, i, i)
                    ELSE IF (NullEdges) THEN
                       CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, i, 1._dp)
                    END IF

                    DO j=A % Rows(i+1)-1,A % Rows(i),-1
                      IF (i>k .OR. A % Cols(j)>k .OR. .NOT.NullEdges) THEN
                        CALL List_AddToMatrixElement(PrecMat % ListMatrix, i, A % Cols(j), A % Values(j))
                      END IF
                    END DO
                  END DO

                  CALL List_ToCRSMatrix(PrecMat)
                  Condition = CRS_IncompleteLU(PrecMat,ILUn,Params)

                  A % ILURows => PrecMat % IluRows
                  A % ILUCols => PrecMat % IluCols
                  A % ILUDiag => PrecMat % IluDiag
                  A % ILUvalues => PrecMat % IluValues

                  DEALLOCATE(PrecMat % Values)
                  IF(.NOT.ASSOCIATED(A % ILURows,PrecMat % Rows)) DEALLOCATE(PrecMat % Rows)
                  IF(.NOT.ASSOCIATED(A % ILUCols,PrecMat % Cols)) DEALLOCATE(PrecMat % Cols)
                  IF(.NOT.ASSOCIATED(A % ILUDiag,PrecMat % Diag)) DEALLOCATE(PrecMat % Diag)
                  DEALLOCATE(PrecMat)
                ELSE
                  Condition = CRS_IncompleteLU(A,ILUn,Params)
                END IF
              CASE(PRECOND_ILUT)
                Condition = CRS_ILUT( A,ILUT_TOL )
              CASE(PRECOND_BILUn)
                Blocks = Solver % Variable % Dofs
                IF ( Blocks <= 1 ) THEN
                  Condition = CRS_IncompleteLU(A,ILUn,Params)
                ELSE
                  IF( .NOT. ASSOCIATED( A % ILUValues ) ) THEN
                    Adiag => AllocateMatrix()
                    CALL CRS_BlockDiagonal(A,Adiag,Blocks)
                    Condition = CRS_IncompleteLU(Adiag,ILUn,Params)
                    A % ILURows   => Adiag % ILURows
                    A % ILUCols   => Adiag % ILUCols
                    A % ILUValues => Adiag % ILUValues
                    A % ILUDiag   => Adiag % ILUDiag                 
                    IF (ILUn > 0) THEN
                      DEALLOCATE(Adiag % Rows,Adiag % Cols, Adiag % Diag, Adiag % Values)
                    END IF
                    DEALLOCATE( Adiag )
                  ELSE
                    Condition = CRS_IncompleteLU(A,ILUn,Params)
                  END IF
                END IF
              CASE(PRECOND_VANKA)
                !                  CALL VankaCreate( A, SolverParam )
              END SELECT
            END IF

            IF( GotDiagFactor ) THEN
              CALL Info('IterSolver','Reverting diagonal relaxation for ILU', Level=10)
               A % Values( A % Diag ) = A % Values( A % Diag ) / DiagFactor
            END IF

          ELSE
            IF ( PCondType == PRECOND_ILUn ) THEN
              CALL Warn( 'IterSolve', 'No ILU Preconditioner for Band Matrix format,' )
              CALL Warn( 'IterSolve', 'using Diagonal preconditioner instead...' )
              PCondType = PRECOND_DIAGONAL
            END IF
          END IF
        END IF
        CALL CheckTimer("Prec-"//TRIM(str),Level=8,Delete=.TRUE.)                  
      END IF
    END IF
    
    A % SolveCount = A % SolveCount + 1
!------------------------------------------------------------------------------

    IF ( PRESENT(MatvecF) ) THEN
      mvProc = MatvecF
    ELSE
      IF ( .NOT. ComplexSystem ) THEN
        mvProc = AddrFunc( CRS_MatrixVectorProd )
      ELSE
        ! A complex matrix is stored fourfold redundantly as 2N real rows of
        ! 2x2 blocks. Taking the product against a compact block view of it
        ! instead is worth roughly 1.8x on the product itself, at the cost of
        ! carrying the view alongside the scalar form. The view itself was
        ! built above, ahead of the preconditioner.
        IF( BlockCRS ) THEN
          mvProc = AddrFunc( CRS_BlockComplexMatrixVectorProd )
        ELSE
          mvProc = AddrFunc( CRS_ComplexMatrixVectorProd )
        END IF
      END IF
    END IF
    
    IF ( PRESENT(dotF) ) THEN
      dotProc = dotF
    ELSE
      dotProc = 0
    END IF
    
    IF ( PRESENT(normF) ) THEN
      normProc = normF
    ELSE
      normProc = 0
    END IF
    
    
    IF ( PRESENT(PrecF) ) THEN
      pcondProc = PrecF
    ELSE
      SELECT CASE( PCondType )
      CASE (PRECOND_NONE)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( pcond_dummy )
        ELSE
          pcondProc = AddrFunc( pcond_dummy_cmplx  )
        END IF
        
      CASE (PRECOND_DIAGONAL)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CRS_DiagPrecondition )
        ELSE
          pcondProc = AddrFunc( CRS_ComplexDiagPrecondition )
        END IF
        
      CASE (PRECOND_ILUn, PRECOND_ILUT, PRECOND_BILUn )
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CRS_LUPrecondition )
        ELSE
          pcondProc = AddrFunc( CRS_ComplexLUPrecondition )
        END IF

      CASE (PRECOND_MG)
        pcondProc = AddrFunc( MultiGridPrec )
        
      CASE (PRECOND_VANKA)
        pcondProc = AddrFunc( VankaPrec )

      CASE (PRECOND_Slave)
        IF(ListGetLogical( Solver % Values,'Linear System Refactorize First',Found ) ) THEN
          CALL LIstAddLogical( Solver % Values,'Linear System Refactorize',.TRUE.)
        END IF        
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( SlavePrec )
        ELSE
          pcondProc = AddrFunc( SlavePrecComplex )
        END IF
        
      CASE (PRECOND_Circuit)
        IF ( .NOT. ComplexSystem ) THEN
          pcondProc = AddrFunc( CircuitPrec )
        ELSE
          pcondProc = AddrFunc( CircuitPrecComplex )
        END IF
        
      CASE DEFAULT
        pcondProc = 0
      END SELECT
    END IF
    

    IF ( .NOT. ComplexSystem ) THEN
      SELECT CASE ( IterType )

       ! Solvers from HUTiter library 
       !-------------------------------------------------------       
      CASE (ITER_BiCGStab)
        iterProc = AddrFunc( HUTI_D_BICGSTAB )
      CASE (ITER_BiCGStab2)
        iterProc = AddrFunc( HUTI_D_BICGSTAB_2 )
      CASE (ITER_TFQMR)
        iterProc = AddrFunc( HUTI_D_TFQMR )
      CASE (ITER_CG)
        iterProc = AddrFunc( HUTI_D_CG )
      CASE (ITER_CGS)
        iterProc = AddrFunc( HUTI_D_CGS )
      CASE (ITER_GMRES)
        iterProc = AddrFunc( HUTI_D_GMRES )
        
        ! Solvers from IterativeMethods.src 
        !-------------------------------------------------------
      CASE (ITER_SGS)
        iterProc = AddrFunc( itermethod_sgs )
      CASE (ITER_JACOBI)
        iterProc = AddrFunc( itermethod_jacobi )
      CASE (ITER_RICHARDSON)
        iterProc = AddrFunc( itermethod_richardson )
      CASE (ITER_GCR)        
        iterProc = AddrFunc( itermethod_gcr )
      CASE (ITER_BICGSTABL)
        iterProc = AddrFunc( itermethod_bicgstabl )
      CASE (ITER_IDRS)
        iterProc = AddrFunc( itermethod_idrs )
      CASE (ITER_MPRGP)
        iterProc = AddrFunc( itermethod_mprgp )
        
      END SELECT
      
      IF( Internal ) THEN
        IF( ListGetLogical( Params,'Linear System Skip Mask',Found ) ) THEN
          CALL Info('IterSolver','Using edge skip mask for linear system solver!')
          dotProc = AddrFunc(MaskedDotProd)
          normproc = AddrFunc(MaskedNorm)        
        ELSE IF( PseudoComplexSystem ) THEN
          IF( HUTI_PSEUDOCOMPLEX == 1 ) THEN
            CALL Info('IterSolver','Setting dot product function to: PseudoZDotProd',Level=15)
            dotProc = AddrFunc( PseudoZDotProd )
          ELSE
            CALL Info('IterSolver','Setting dot product function to: PseudoZDotProd2',Level=15)
            dotProc = AddrFunc( PseudoZDotProd2 )             
          END IF
        ELSE        
!         IF ( dotProc  == 0 ) dotProc = AddrFunc(ddot)
        END IF
        IF ( normProc == 0 ) normproc = AddrFunc(dnrm2)
        IF( HUTI_DBUGLVL == 0) HUTI_DBUGLVL = HUGE( HUTI_DBUGLVL )        
      END IF

      IF ( dotProc  == 0 ) dotProc = AddrFunc(Otmp_ddot)
      
    ELSE
      HUTI_NDIM = HUTI_NDIM / 2
      SELECT CASE ( IterType )

        ! Solvers from HUTiter library 
        !-------------------------------------------------------       
      CASE (ITER_BiCGStab)
        iterProc = AddrFunc( HUTI_Z_BICGSTAB )
      CASE (ITER_BiCGStab2)
        iterProc = AddrFunc( HUTI_Z_BICGSTAB_2 )
      CASE (ITER_TFQMR)
        iterProc = AddrFunc( HUTI_Z_TFQMR )
      CASE (ITER_CG)
        iterProc = AddrFunc( HUTI_Z_CG )
      CASE (ITER_CGS)
        iterProc = AddrFunc( HUTI_Z_CGS )
      CASE (ITER_GMRES)
        iterProc = AddrFunc( HUTI_Z_GMRES )
        
        ! Solvers from IterativeMethods.src 
        !-------------------------------------------------------
      CASE (ITER_GCR)
        iterProc = AddrFunc( itermethod_z_gcr )
      CASE (ITER_BICGSTABL)
        iterProc = AddrFunc( itermethod_z_bicgstabl )
      CASE (ITER_IDRS)
        iterProc = AddrFunc( itermethod_z_idrs )
      CASE DEFAULT
        CALL Fatal('IterSolver', 'Complex arithmetic version of the given linear solver is not available')
      END SELECT
      
      IF( Internal ) THEN
!       IF ( dotProc  == 0 ) dotProc = AddrFunc(zdotc)
        IF ( normProc == 0 ) normproc = AddrFunc(dznrm2)
        IF( HUTI_DBUGLVL == 0) HUTI_DBUGLVL = HUGE( HUTI_DBUGLVL )
      END IF

      IF ( IterType == ITER_CG ) THEN
        ! CG needs the unconjugated bilinear form on these complex symmetric
        ! systems; the other complex methods want the Hermitian product. This
        ! has to be honoured for a caller-supplied product too, not just the
        ! default: the parallel path hands in SParCDotProd, and with the
        ! conjugating form CG has no self-adjoint operator to work with and
        ! diverges outright -- HelmholtzFEM's CG pass ran the residual up to
        ! 5.5e+01 over 301 iterations at np=2 and np=4 alike before this.
        IF ( PRESENT( DotFU ) ) THEN
          dotProc = DotFU
        ELSE IF ( dotProc == 0 ) THEN
          dotProc = AddrFunc(Otmp_zdotu)
        END IF
      ELSE IF ( dotProc == 0 ) THEN
        dotProc = AddrFunc(Otmp_zdotc)
      END IF

    END IF
    
!------------------------------------------------------------------------------

    stack_pos = stack_pos+1
    IF(stack_pos>stack_max) THEN
      CALL Fatal('IterSolver', 'Recursion too deep ('//I2S(stack_pos)//' vs '//I2S(stack_max)//')')
    ELSE IF(stack_pos<=0) THEN
      CALL Fatal('IterSolver', 'eh')
    END IF
    FirstCall(stack_pos) = .TRUE.

    SaveGlobalM => GlobalMatrix
    GlobalMatrix => A
    
    IF ( ComplexSystem ) THEN
      ! x and b already hold the complex vectors as consecutive (Re,Im) pairs,
      ! so alias them instead of copying into complex temporaries. The solution
      ! is then written straight back into x and needs no copy-out either.
      xC => ComplexValues( x, HUTI_NDIM )
      bC => ComplexValues( b, HUTI_NDIM )

      CALL Info('IterSolver','Calling complex iterative solver',Level=32)

      ! Release the scalar values across the iteration. Once the block view
      ! carries the product and the ILU factors are built, nothing in the Krylov
      ! loop reads A % Values -- and that array is two thirds of the matrix,
      ! dropped exactly where the footprint peaks, with factors, matrix and
      ! Krylov work vectors all resident at once. It is rebuilt from the view on
      ! the way out, exactly, so nothing downstream can tell.
      !
      ! ELIGIBILITY IS A WHITELIST, and deliberately so: the readers cannot be
      ! enumerated by grep (they reach the matrix through GlobalMatrix), so the
      ! rule admits only what has been checked rather than excluding what has
      ! been noticed. Known in-window readers that must keep the array:
      !   CRS_ComplexDiagPrecondition   used to read Values every iteration and
      !                                 excluded diagonal preconditioning
      !                                 outright; it now takes the diagonal
      !                                 block from the view, so it is admitted
      !                                 whenever BDiag exists. Re-checked with
      !                                 the poison probe below, which failed on
      !                                 this case before the conversion and
      !                                 passes after it.
      !   the backward-error stop criteria, e.g. SUM(GlobalMatrix % Values**2)
      !                                 in BackwardError.F90 -- hence StopcProc==0
      !   a caller-supplied product, which may read anything, so it has to
      !                                 assert otherwise through
      !                                 MatvecReadsNoValues. SParCMatrixVector
      !                                 does, on the strength of the poison
      !                                 probe run at np4.
      ! Established by poisoning the array and watching what breaks, not by
      ! reading: "Linear System Poison Scalar Values" below is that probe, kept
      ! because it is what any widening of this list has to be justified with.
      !
      ! The aliasing guard is not decoration. ParallelUtils points
      ! InsideMatrix % Values at InsideMatrix % MassValues for part of the eigen
      ! path, and DEALLOCATE on a pointer that was pointer-assigned rather than
      ! allocated is undefined. ASSOCIATED(a,b) is the one thing that can be
      ! tested here, so test it.
      ! Released for the iteration and not earlier, deliberately. Releasing before
      ! the factorization instead was built and measured and gained NOTHING: the
      ! peak of a run like this is set by mesh handling and assembly, well before
      ! the solve, so all a release can do is clip the solve phase back under
      ! that ceiling. Once it is under, freeing sooner or freeing more is
      ! invisible. See section 7p of complex-storage-estimate.txt.
      IF( .NOT. ValsFreed ) THEN
        FreeVals = FreeEligible( A, Params, PCondType, StopcProc, &
            BlockCRS, PoisonVals, MatvecF, MatvecReadsNoValues )
        IF( FreeVals ) THEN
          nScalarVals = SIZE( A % Values )
          DEALLOCATE( A % Values )
          A % Values => NULL()
          ValsFreed = .TRUE.
          CALL Info('IterSolver','Released the scalar values across the iteration: '// &
              I2S(nScalarVals)//' reals',Level=8)
        END IF
      END IF

      IF (LeftOriented) THEN
        CALL IterCall( iterProc, xC, bC, ipar, dpar, workC, &
            mvProc, pcondProc, pconddProc, dotProc, normProc, stopcProc )
      ELSE
        CALL IterCall( iterProc, xC, bC, ipar, dpar, workC, &
            mvProc, pconddProc, pcondProc, dotProc, normProc, stopcProc )
      END IF

      IF( ValsFreed ) THEN
        ALLOCATE( A % Values(nScalarVals) )
        CALL CRS_ExpandBlockCRS( A )
        ValsFreed = .FALSE.
      END IF

      IF( PoisonVals ) THEN
        A % Values = SaveVals
        DEALLOCATE( SaveVals )
      END IF

      ! No copy-back needed: xC aliases x.
      xC => NULL()
      bC => NULL()
    ELSE
      CALL Info('IterSolver','Calling real-valued iterative solver',Level=32)

      IF (LeftOriented) THEN
        CALL IterCall( iterProc, x, b, ipar, dpar, work, &
            mvProc, pcondProc, pconddProc, dotProc, normProc, stopcProc )          
      ELSE
        CALL IterCall( iterProc, x, b, ipar, dpar, work, &
            mvProc, pconddProc, pcondProc, dotProc, normProc, stopcProc )
      END IF
    ENDIF

    GlobalMatrix => SaveGlobalM
    
    stack_pos=stack_pos-1
    
    IF ( ComplexSystem ) HUTI_NDIM = HUTI_NDIM * 2

    !------------------------------------------------------------------------------
    IF ( HUTI_INFO == HUTI_CONVERGENCE ) THEN
      IF( ASSOCIATED( Solver % Variable ) ) THEN
        Solver % Variable % LinConverged = 1
      END IF
    ELSE
      CALL Info('IterSolve','Returned return code: '//I2S(HUTI_INFO),Level=15)
      IF( HUTI_INFO == HUTI_DIVERGENCE ) THEN
        CALL NumericalError( 'IterSolve', 'System diverged over maximum tolerance.')
      ELSE IF( HUTI_INFO == HUTI_MAXITER ) THEN                
        DoFatal = ListGetLogical( Params,'Linear System Abort Not Converged',Found )
        IF(.NOT. Found ) DoFatal = .TRUE.
        IF( DoFatal ) THEN
          CALL NumericalError('IterSolve','Too many iterations were needed.')
        ELSE
          CALL Info('IterSolve','Linear iteration did not converge to tolerance',Level=6)
        END IF
      ELSE IF( HUTI_INFO == HUTI_HALTED ) THEN
        CALL Warn('IterSolve','Iteration halted due to problem in algorithm, trying to continue')
      END IF
      IF( ASSOCIATED( Solver % Variable ) ) THEN
        Solver % Variable % LinConverged = 0
      END IF
    END IF
!------------------------------------------------------------------------------
    IF ( ComplexSystem ) THEN
      DEALLOCATE( workC )
    ELSE
      DEALLOCATE( work )
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE IterSolver
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> This routine may be used to either inform user or terminate following
!> convergence/numerical issues, based on a flag in the SIF. Default
!> behaviour terminates execution.
!-----------------------------------------------------------------------
   SUBROUTINE NumericalError( Caller, String, IsFatal )
!-----------------------------------------------------------------------
     CHARACTER(LEN=*) :: Caller, String
     LOGICAL, OPTIONAL :: IsFatal
!-----------------------------------------------------------------------
     LOGICAL :: DoFatal, Found
!-----------------------------------------------------------------------

     !Fatality logic:
     ! 1) Respect calling routine's wishes if present
     ! 2) Respect solver specific option if present
     ! 3) Respect global abort flag if present
     ! 4) Otherwise fatal (backwards compatibility)

     IF(PRESENT(IsFatal)) THEN
       DoFatal = IsFatal
     ELSE
       DoFatal = ListGetLogical(CurrentModel % Simulation,&
           'Global Abort Not Converged',Found)
       IF(.NOT. Found ) DoFatal = .TRUE.
     END IF

     IF(DoFatal) THEN
       CALL Fatal(Caller,'Numerical Error: '//TRIM(String))
     ELSE
       CALL Warn(Caller,'Numerical Error: '//TRIM(String))
     END IF

!-----------------------------------------------------------------------
   END SUBROUTINE NumericalError
!-----------------------------------------------------------------------


END MODULE IterSolve

!> \} ElmerLib
! a trivial whitespace
