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
! *  Original Date: 01 Oct 1996
! *
! ******************************************************************************/

!--------------------------------------------------------------------------------
!>  Module defining element type and operations. The most basic FEM routines
!>  are here, handling the basis functions, global derivatives, etc...
!--------------------------------------------------------------------------------
!> \ingroup ElmerLib
!> \{

#include "../config.h"

!--------------------------------------------------------------------------------
!>  Element information assembly: ElementInfo, EdgeElementInfo, FaceElementInfo,
!>  metric computation, global derivative computation. StabParam also lives here.
!--------------------------------------------------------------------------------
MODULE ElemInfo
   USE ElementBasis
   USE ElementGeometry
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE H1Basis
   USE Lists
!$ USE omp_lib ! Include module conditionally (for omp_in_parallel below)

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: GetElementType, StabParam, &
             EdgeElementStyle, ElementInfo, ElementInfoVec, &
             ElementInfoVec_ComputePElementBasis, ElementInfoVec_ElementBasisToGlobal, &
             ElementSize, &
             FaceElementInfo, PiolaTransformationData, &
             FaceElementOrientation, FaceElementBasisOrdering, PickActiveFace, &
             EdgeElementInfo, &
             TriangleFaceDofsOrdering, TriangleFaceDofsOrdering2nd, &
             TriangleFaceDofsOrdering2, SquareFaceDofsOrdering, &
             ReorderingAndSignReversionsData, &
             EdgeWhitneyComponents2D, FaceWhitneyComponents2D, WeightedWhitneyForms, &
             GetEdgeBasis, mGetElementDOFs, &
             ElementMetric, ElementMetricVec, &
             GlobalFirstDerivativesInternal, GlobalFirstDerivatives, &
             GlobalSecondDerivatives
#ifdef HAVE_QP
   PUBLIC :: CheckMetric, ElementMetricQP
#endif

CONTAINS

   FUNCTION GetElementType( code,CompStabFlag ) RESULT(element)
!------------------------------------------------------------------------------
      INTEGER :: code
      LOGICAL, OPTIONAL :: CompStabFlag
      TYPE(ElementType_t), POINTER :: element
!------------------------------------------------------------------------------
!     Local variables
!------------------------------------------------------------------------------
      TYPE(Nodes_t) :: Nodes
      INTEGER :: sdim
      TYPE(Element_t), POINTER :: Elm
!------------------------------------------------------------------------------
      element => ElementTypeList

      DO WHILE( ASSOCIATED(element) )
        IF ( code == element % ElementCode ) EXIT
        element => element % NextElementType
      END DO

      IF ( .NOT. ASSOCIATED( element ) ) THEN
        WRITE( message, * ) &
             'Element type code ',code,' not found. Ignoring element.'
        CALL Warn( 'GetElementType', message )
        RETURN
      END IF

      IF ( PRESENT( CompStabFlag ) ) THEN
        IF ( .NOT. CompStabFlag ) RETURN
      END IF

      IF ( Element % StabilizationMK == 0.0d0 ) THEN
        ALLOCATE( Elm )
        Elm % TYPE => element
        Elm % BDOFs  = 0
        Elm % DGDOFs = 0
        NULLIFY( Elm % PDefs )
        NULLIFY( Elm % DGIndexes )
        NULLIFY( Elm % EdgeIndexes )
        NULLIFY( Elm % FaceIndexes )
        NULLIFY( Elm % BubbleIndexes )
        Nodes % x => Element % NodeU
        Nodes % y => Element % NodeV
        Nodes % z => Element % NodeW

        sdim = CurrentModel % Dimension
        CurrentModel % Dimension = Element % Dimension
        CALL StabParam( Elm, Nodes, Element % NumberOfNodes, &
                 Element % StabilizationMK )
        CurrentModel % Dimension = sdim

        DEALLOCATE(Elm)
      END IF

   END FUNCTION GetElementType
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute convection diffusion equation stab. parameter  for each and every
!> element of the model by solving the largest eigenvalue of
!
!> Lu = \lambda Gu,
!
!> L = (\nablda^2 u,\nabla^ w), G = (\nabla u,\nabla w)
!------------------------------------------------------------------------------
   SUBROUTINE StabParam(Element,Nodes,n,mK,hK,UseLongEdge)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t) :: Element
      INTEGER :: n
      TYPE(Nodes_t) :: Nodes
      REAL(KIND=dp) :: mK
      REAL(KIND=dp), OPTIONAL :: hK
      LOGICAL, OPTIONAL :: UseLongEdge
!------------------------------------------------------------------------------
      INTEGER :: info,p,q,i,j,t,dim
      REAL(KIND=dp) :: EIGR(n),EIGI(n),Beta(n),s,ddp(3),ddq(3),dNodalBasisdx(n,n,3)
      REAL(KIND=dp) :: u,v,w,L(n-1,n-1),G(n-1,n-1),Work(16*n)
      REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),ddBasisddx(n,3,3),detJ

      LOGICAL :: stat
      TYPE(GaussIntegrationPoints_t) :: IntegStuff

      IF ( Element % TYPE % BasisFunctionDegree <= 1 ) THEN
         SELECT CASE( Element % TYPE % ElementCode ) 
           CASE( 202, 303, 404, 504, 605, 706  )
              mK = 1.0d0 / 3.0d0
           CASE( 808 )
              mK = 1.0d0 / 6.0d0
         END SELECT
         IF ( PRESENT( hK ) ) hK = ElementDiameter( Element, Nodes, UseLongEdge)
         RETURN
      END IF

      dNodalBasisdx = 0._dp
      DO p=1,n
        u = Element % TYPE % NodeU(p)
        v = Element % TYPE % NodeV(p)
        w = Element % TYPE % NodeW(p)
        stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
        dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
      END DO

      dim = CoordinateSystemDimension()
      IntegStuff = GaussPoints( Element )
      L = 0.0d0
      G = 0.0d0
      DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis, &
                dBasisdx )

        s = detJ * IntegStuff % s(t)

        DO p=2,n
          DO q=2,n
            ddp = 0.0d0
            ddq = 0.0d0
            DO i=1,dim
              G(p-1,q-1) = G(p-1,q-1) + s * dBasisdx(p,i) * dBasisdx(q,i)
              ddp(i) = ddp(i) + SUM( dNodalBasisdx(p,1:n,i) * dBasisdx(1:n,i) )
              ddq(i) = ddq(i) + SUM( dNodalBasisdx(q,1:n,i) * dBasisdx(1:n,i) )
            END DO
            L(p-1,q-1) = L(p-1,q-1) + s * SUM(ddp) * SUM(ddq)
          END DO
        END DO
      END DO

      IF ( ALL(ABS(L) < AEPS) ) THEN
        mK = 1.0d0 / 3.0d0
        IF ( PRESENT(hK) ) THEN
          hK = ElementDiameter( Element,Nodes,UseLongEdge)
        END IF
        RETURN
      END IF


      CALL DSYGV( 1,'N','U',n-1,L,n-1,G,n-1,EIGR,Work,12*n,info )
      mK = EIGR(n-1)

      IF ( mK < 10*AEPS ) THEN
        mK = 1.0d0 / 3.0d0
        IF ( PRESENT(hK) ) THEN
          hK = ElementDiameter( Element,Nodes,UseLongEdge )
        END IF
        RETURN
      END IF

      IF ( PRESENT( hK ) ) THEN
        hK = SQRT( 2.0d0 / (mK * Element % TYPE % StabilizationMK) )
        mK = MIN( 1.0d0 / 3.0d0, Element % TYPE % StabilizationMK )
      ELSE
        SELECT CASE(Element % TYPE % ElementCode / 100)
        CASE(2,4,8) 
          mK = 4 * mK
        END SELECT
        mK = MIN( 1.0d0/3.0d0, 2/mK )
      END IF

!------------------------------------------------------------------------------
   END SUBROUTINE StabParam
   SUBROUTINE EdgeElementStyle(VList, PiolaVersion, SecondFamily, QuadraticApproximation, &
       BasisDegree, GradientVersion, Check ) 

     TYPE(ValueList_t), POINTER :: VList
     LOGICAL :: PiolaVersion
     LOGICAL, OPTIONAL :: SecondFamily
     LOGICAL, OPTIONAL :: QuadraticApproximation
     INTEGER, OPTIONAL :: BasisDegree
     LOGICAL, OPTIONAL :: GradientVersion
     LOGICAL, OPTIONAL :: Check
     
     LOGICAL :: Found, Quadratic, Cubic, Second 
     
     Quadratic = ListGetLogical(VList,'Quadratic Approximation', Found )
     Cubic = ListGetLogical(VList,'Cubic Approximation', Found )
     
     Second = ListGetLogical(Vlist,'Second Kind Basis', Found )
     IF( Quadratic .OR. Cubic) THEN
       PiolaVersion = .TRUE.
     ELSE       
       IF(Second) THEN
         PiolaVersion = .TRUE.
       ELSE    
         PiolaVersion = ListGetLogical(Vlist,'Use Piola Transform', Found )
       END IF
     END IF

     IF(PRESENT(SecondFamily)) THEN
       SecondFamily = Second
     END IF
     
     IF(PRESENT(BasisDegree)) THEN
       BasisDegree = 1
       IF(Quadratic) THEN
         BasisDegree = 2
       ELSE IF (Cubic) THEN
         BasisDegree = 3
       END IF
     END IF

     IF(PRESENT(QuadraticApproximation)) THEN
       QuadraticApproximation = Quadratic
     END IF

     IF (PRESENT(GradientVersion)) THEN
       GradientVersion = ListGetLogical(VList, 'Gradient Basis Functions', Found) .OR. &
           ListGetLogical(VList, 'Simplicial Mesh', Found)
     END IF
     
     ! When initializing the consistency of the keywords may be checked.
     ! Also always add the Piola flag since it determines the type of IPs.
     IF( PRESENT(Check)) THEN
       IF(Check) THEN
         IF(PiolaVersion) THEN
           IF(.NOT. ListCheckPresent(Vlist,'Use Piola Transform')) THEN
             IF(Quadratic .OR. Cubic) THEN
               CALL Info('EdgeElementStyle','"Quadratic/Cubic Approximation" requested without Piola. ' &
                   //'Setting "Use Piola Transform = True"')
             ELSE IF( Second ) THEN
               CALL Info('EdgeElementStyle','"Second Kind Basis" requested without Piola. ' &
                   //'Setting "Use Piola Transform = True"')
             END IF
             CALL ListAddLogical(Vlist,'Use Piola Transform',.TRUE.)
           END IF
         END IF
       END IF
     END IF    
     
   END SUBROUTINE EdgeElementStyle

   
!------------------------------------------------------------------------------
!>  Return the referential description b(f(p)) of the basis function b(x),
!>  with f mapping points p on a reference element to points x on a physical
!>  element. The referential description of the spatial gradient field grad b 
!>  and, if requested, the second spatial derivatives may also be returned.
!>  Also return the square root of the determinant of the metric tensor
!>  (=sqrt(det(J^TJ))) related to the mapping f.
!------------------------------------------------------------------------------
   RECURSIVE FUNCTION ElementInfo( Element, Nodes, u, v, w, detJ, &
       Basis, dBasisdx, ddBasisddx, SecondDerivatives, Bubbles, BasisDegree, &
       EdgeBasis, RotBasis, USolver, ip_index ) RESULT(stat)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element             !< Element structure
     TYPE(Nodes_t)   :: Nodes                       !< Element nodal coordinates.
     REAL(KIND=dp) :: u                             !< 1st local coordinate at which to calculate the basis function.
     REAL(KIND=dp) :: v                             !< 2nd local coordinate.
     REAL(KIND=dp) :: w                             !< 3rd local coordinate.
     REAL(KIND=dp) :: detJ                          !< Square root of determinant of element coordinate system metric
     REAL(KIND=dp) :: Basis(:)                      !< Basis function values at p=(u,v,w)
     REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)       !< Global first derivatives of basis functions at (u,v,w)
     REAL(KIND=dp), OPTIONAL :: ddBasisddx(:,:,:)   !< Global second derivatives of basis functions at (u,v,w) if requested
     LOGICAL, OPTIONAL :: SecondDerivatives         !< Are the second derivatives needed? (still present for historical reasons)
     LOGICAL, OPTIONAL :: Bubbles                   !< Are the bubbles to be evaluated.
     INTEGER, OPTIONAL :: BasisDegree(:)            !< Degree of each basis function in Basis(:) vector. 
                                                    !< May be used with P element basis functions
     REAL(KIND=dp), OPTIONAL :: EdgeBasis(:,:)      !< If present, the values of H(curl)-conforming basis functions B(f(p))
     REAL(KIND=dp), OPTIONAL :: RotBasis(:,:)       !< The referential description of the spatial curl of B
     TYPE(Solver_t), POINTER, OPTIONAL :: USolver   !< The solver used to call the basis functions.
     INTEGER, OPTIONAL, INTENT(IN) :: ip_index      !< Integration point index for O(1) cache lookup.
     LOGICAL :: Stat                                !< If .FALSE. element is degenerate.
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PSolver => NULL(), PrevSolver => NULL()

     LOGICAL :: Compute2ndDerivatives

     INTEGER :: i, j, k, q, n, dim, cdim, ip_slot
     REAL(KIND=dp) :: ElmMetric(3,3), LtoGMap(3,3)

     REAL(KIND=dp) :: NodalBasis(Element % TYPE % NumberOfNodes), &
             dLBasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3)

     REAL(KIND=dp), ALLOCATABLE :: ddlBasisddx(:,:,:)

     INTEGER :: EdgeBasisDegree
     LOGICAL :: SecondFamily, GradVersion
     LOGICAL :: PerformPiolaTransform, Found, SerendipityPBasis
     
     SAVE PrevSolver, EdgeBasisDegree, PerformPiolaTransform, SecondFamily, GradVersion
     !$OMP THREADPRIVATE(PrevSolver, EdgeBasisDegree, PerformPiolaTransform, SecondFamily, &
     !$OMP               GradVersion)
!------------------------------------------------------------------------------

     IF( PRESENT( USolver ) ) THEN
       pSolver => USolver
     ELSE
       pSolver => CurrentModel % Solver
     END IF
     
     IF(PRESENT(EdgeBasis)) THEN       
       IF( .NOT. ASSOCIATED( PrevSolver, PSolver ) ) THEN
         PrevSolver => pSolver                  
         CALL EdgeElementStyle(pSolver % Values, PerformPiolaTransform, SecondFamily, &
             BasisDegree = EdgeBasisDegree, GradientVersion = GradVersion)
       END IF
       IF( PerformPiolaTransform ) THEN
         stat = EdgeElementInfo(Element,Nodes,u,v,w,detF=Detj,Basis=Basis, &
             EdgeBasis=EdgeBasis,RotBasis=RotBasis,dBasisdx=dBasisdx,&
             SecondFamily = SecondFamily, BasisDegree = EdgeBasisDegree, &
             ApplyPiolaTransform = PerformPiolaTransform, &
             GradientVersion = GradVersion)
       ELSE
         IF(Element % Type % ElementCode == 504 .AND. ANY([u,v,w] < 0.0) ) THEN
           PRINT *,'Negative local coordinates for tet:',u,v,w
         END IF
         stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )
         CALL GetEdgeBasis(Element,EdgeBasis,RotBasis,Basis,dBasisdx)         
       END IF
       RETURN
     END IF

     stat = .TRUE.
     n    = Element % Type % NumberOfNodes
     dim  = Element % Type % Dimension
     cdim = CoordinateSystemDimension()

     IF ( Element % Type % ElementCode == 101 ) THEN
        detJ = 1.0d0
        Basis(1) = 1.0d0
        IF ( PRESENT(dBasisdx) ) dBasisdx(1,:) = 0.0d0
        RETURN
     END IF

     Compute2ndDerivatives = PRESENT(SecondDerivatives) .AND. PRESENT(ddBasisddx)
     IF(Compute2ndDerivatives) Compute2ndDerivatives = SecondDerivatives

     IF (Compute2ndDerivatives) &
       CALL EvalSecondDerivativesRef(Element, pSolver, u, v, w, n, dim, Basis, &
           MAX(SIZE(Nodes % x),SIZE(ddBasisddx)), ddLBasisddx)

     ip_slot = 0
     IF (PRESENT(ip_index)) ip_slot = ip_index
     ! Evaluate reference basis functions and local gradients at (u,v,w).
     ! Reference basis cache lookup/fill. The cache lives on Element % Type,
     ! which is SHARED by every element of this type across all OpenMP threads
     ! (a mesh has only a handful of distinct Type objects). A single shared
     ! cache filled from inside a parallel assembly loop was an unsynchronized
     ! data race: concurrent threads that cache-miss grabbed the same next slot
     ! (ip = BasisCacheCount+1) and tore each other's (U,V,W)-key vs payload
     ! writes, so a later lookup could match a key while reading a different
     ! point's basis values -> wrong local matrix -> intermittently, locally
     ! wrong solution (the ~1% ElastElstatBeamNodal wrong-norm flake on Windows
     ! CI; invisible to Valgrind memcheck since every access is in bounds -- a
     ! race, not a memory error). Each thread has its own cache column (last
     ! index tid), so both lookup and fill are lock-free and race-free and the
     ! cache stays active during threaded assembly. Same per-thread pattern as
     ! DefUtils' index stores.
     RefBasisBlock: BLOCK
       INTEGER :: ip, tid
       tid = 1
       !$ tid = omp_get_thread_num() + 1
       IF ( ip_slot > 0 ) THEN
         ! O(1) direct slot lookup
         IF ( ip_slot <= Element % TYPE % BasisCacheCount(tid) ) THEN
           Basis(1:n)       = Element % Type % BasisCache(ip_slot, 1:n, tid)
           dLBasisdx(1:n,:) = Element % Type % dBasisCache(ip_slot, 1:n, :, tid)
           EXIT RefBasisBlock
         END IF
       ELSE
         ! Linear scan by (u,v,w) coordinates
         DO ip = 1, Element % TYPE % BasisCacheCount(tid)
           IF ( Element % Type % BasisCacheU(ip, tid) == u .AND. &
                Element % Type % BasisCacheV(ip, tid) == v .AND. &
                Element % Type % BasisCacheW(ip, tid) == w ) THEN
             Basis(1:n)       = Element % TYPE % BasisCache(ip, 1:n, tid)
             dLBasisdx(1:n,:) = Element % TYPE % dBasisCache(ip, 1:n, :, tid)
             EXIT RefBasisBlock
           END IF
         END DO
       END IF

       CALL NodalBasisFunctions(n, Basis, element, u, v, w, pSolver)
       CALL NodalFirstDerivatives(n, dLBasisdx, element, u, v, w, pSolver)

       ! Store in this thread's cache column if space available.
       IF ( ip_slot > 0 ) THEN
         ip = ip_slot
       ELSE
         ip = Element % Type % BasisCacheCount(tid) + 1
       END IF

       IF ( ip <= ELEM_BASIS_CACHE_SIZE ) THEN
         Element % Type % BasisCacheU(ip, tid) = u
         Element % Type % BasisCacheV(ip, tid) = v
         Element % Type % BasisCacheW(ip, tid) = w
         Element % Type % BasisCache(ip, 1:n, tid)     = Basis(1:n)
         Element % Type % dBasisCache(ip, 1:n, :, tid) = dLBasisdx(1:n, :)
         Element % Type % BasisCacheCount(tid) = MAX(Element % Type % BasisCacheCount(tid), ip)
       END IF
     END BLOCK RefBasisBlock

     q = n
     CALL EvalPElementBasis(Element, pSolver, u, v, w, n, q, Basis, dLBasisdx, &
             Compute2ndDerivatives, ddLBasisddx, BasisDegree)

!------------------------------------------------------------------------------

     ! Element (contravariant) metric and square root of determinant
     !--------------------------------------------------------------
#ifdef HAVE_QP
     IF(Element % Status==0) THEN
       stat = CheckMetric(q, Element, Nodes, dLBasisdx)
       Element % Status = MERGE(1, 2, stat)
     END IF
#endif
     stat = ElementMetric(q,Element,Nodes,ElmMetric,detJ,dLBasisdx,LtoGMap)
     IF ( .NOT. stat ) RETURN
     IF ( PRESENT(dBasisdx) ) THEN
       dBasisdx = 0.0d0
       DO k=1,dim
         DO j=1,cdim
           DO i=1,q
             dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k) * LtoGMap(j,k)
           END DO
         END DO
       END DO
     END IF

     ! Get matrix of second derivatives, if needed:
     !---------------------------------------------
     IF (Compute2ndDerivatives) &
       CALL GlobalSecondDerivatives(Element,Nodes,ddBasisddx,u,v,w,ElmMetric,dLBasisdx,ddLBasisddx,q)

!------------------------------------------------------------------------------
     IF ( PRESENT(Bubbles) .AND. .NOT. isActivePElement(Element, pSolver) ) THEN
       CALL EvalBubbleBasis(Element, Nodes, u, v, w, detJ, n, cdim, &
           Basis, dBasisdx, Bubbles, stat)
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ElementInfo


   SUBROUTINE EvalSecondDerivativesRef(Element, pSolver, u, v, w, &
       n, dim, Basis, nalloc, ddLBasisddx)
     IMPLICIT NONE
     TYPE(Element_t), TARGET, INTENT(IN) :: Element
     TYPE(Solver_t), POINTER, INTENT(IN) :: pSolver
     REAL(KIND=dp), INTENT(IN) :: u, v, w
     INTEGER, INTENT(IN) :: n, dim, nalloc
     REAL(KIND=dp), INTENT(INOUT) :: Basis(:)
     REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ddLBasisddx(:,:,:)
     INTEGER :: i
     ALLOCATE(ddLBasisddx(nalloc,3,3))
     Basis = 0
     ddLBasisddx = 0._dp
     DO i=1,n
       Basis(i) = 1
       SELECT CASE(dim)
       CASE(1)
         ddLBasisddx(i,1,1) = SecondDerivatives1D(element,basis,u)
       CASE(2)
         ddLBasisddx(i,1:2,1:2) = SecondDerivatives2D(element,basis,u,v)
       CASE(3)
         SELECT CASE(Element % Type % ElementCode)
         CASE(605)
           IF(isActivePElement(Element,pSolver)) THEN
             ddLBasisddx(i,:,:) = ddPyramidNodalPBasis(i,u,v,w)
           ELSE
             ddLBasisddx(i,:,:) = SecondDerivatives3D(element,basis,u,v,w)
           END IF
         CASE(706)
           IF(isActivePElement(element,pSolver)) THEN
             ddLBasisddx(i,:,:) = ddWedgeNodalPBasis(i,u,v,w)
           ELSE
             ddLBasisddx(i,:,:) = SecondDerivatives3D(element,basis,u,v,w)
           END IF
         CASE DEFAULT
           ddLBasisddx(i,:,:) = SecondDerivatives3D(element,basis,u,v,w)
         END SELECT
       END SELECT
       Basis(i) = 0
     END DO
   END SUBROUTINE EvalSecondDerivativesRef


   SUBROUTINE EvalBubbleBasis(Element, Nodes, u, v, w, detJ, n, cdim, &
       Basis, dBasisdx, Bubbles, stat)
     IMPLICIT NONE
     TYPE(Element_t), TARGET, INTENT(IN) :: Element
     TYPE(Nodes_t), INTENT(IN) :: Nodes
     REAL(KIND=dp), INTENT(IN) :: u, v, w
     REAL(KIND=dp), INTENT(INOUT) :: detJ
     INTEGER, INTENT(IN) :: n, cdim
     REAL(KIND=dp), INTENT(INOUT) :: Basis(:), dBasisdx(:,:)
     LOGICAL, INTENT(IN) :: Bubbles
     LOGICAL, INTENT(INOUT) :: stat
     REAL(KIND=dp) :: BubbleValue, LinBasis(8), dLinBasisdx(8,3)
     TYPE(Element_t) :: Bubble
     INTEGER :: i, j

       Bubble % BDOFs = 0
       NULLIFY( Bubble % PDefs )
       NULLIFY( Bubble % EdgeIndexes )
       NULLIFY( Bubble % FaceIndexes )
       NULLIFY( Bubble % BubbleIndexes )

       IF ( Bubbles .AND. SIZE(Basis) >= 2*n ) THEN

         SELECT CASE(Element % TYPE % ElementCode / 100)
           CASE(2)

              IF ( Element % TYPE % ElementCode == 202 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(202)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                          LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1)
                END DO
              END DO

           CASE(3)

              IF ( Element % TYPE % ElementCode == 303 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(303)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF
  
              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(3,j) * LinBasis(1) * LinBasis(2)
                END DO
              END DO

           CASE(4)

              IF ( Element % TYPE % ElementCode == 404 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(404)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                             LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(1,j) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(3,j) * LinBasis(1)
                END DO
              END DO

           CASE(5)

              IF ( Element % TYPE % ElementCode == 504 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(504)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3) * LinBasis(4)
              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(1,j) * &
                                    LinBasis(2) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(2,j) * &
                                    LinBasis(1) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(3,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(4,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(3)
                END DO
              END DO

           CASE(8)

              IF ( Element % TYPE % ElementCode == 808 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(808)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                  LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(7)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(1,j) * LinBasis(7)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(7,j) * LinBasis(1)
                END DO
              END DO

         CASE DEFAULT
 
              WRITE( Message, '(a,i4,a)' ) 'Bubbles for element: ', &
               Element % TYPE % ElementCode, ' are not implemented.'
              CALL Error( 'ElementInfo', Message )
              CALL Fatal( 'ElementInfo', 'Please use p-element basis instead.' )

         END SELECT
       END IF
   END SUBROUTINE EvalBubbleBasis


   SUBROUTINE EvalPElementBasis(Element, pSolver, u, v, w, n, q, Basis, dLBasisdx, &
       Compute2ndDerivatives, ddLBasisddx, BasisDegree)
     IMPLICIT NONE
     TYPE(Element_t), TARGET, INTENT(IN) :: Element
     TYPE(Solver_t), POINTER, INTENT(IN) :: pSolver
     REAL(KIND=dp), INTENT(IN) :: u, v, w
     INTEGER, INTENT(IN) :: n
     INTEGER, INTENT(INOUT) :: q
     REAL(KIND=dp), INTENT(INOUT) :: Basis(:), dLBasisdx(:,:)
     LOGICAL, INTENT(IN) :: Compute2ndDerivatives
     REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: ddLBasisddx(:,:,:)
     INTEGER, OPTIONAL, INTENT(INOUT) :: BasisDegree(:)
     LOGICAL :: degrees, invert, SerendipityPBasis
     INTEGER :: i, j, k, f, p, nb, EDOFs, BDOFs, BodyId, tetraType, locali, localj
     INTEGER :: tmp(4), direction(4)
     INTEGER :: GIndexes(Element % Type % NumberOfNodes)
     TYPE(Element_t), POINTER :: Parent, Edge, Face
      ! P ELEMENT CODE:
      ! ---------------
       !
       ! Check whether the polynomial degree of each basis functions is asked
       ! and, if so, initialize by the degree of linear basis:
       ! ---------------------------------------------------
     IF ( isActivePElement(element,pSolver) ) THEN
       degrees = .FALSE.
       IF ( PRESENT(BasisDegree)) THEN 
         degrees = .TRUE.
         BasisDegree = 0
         BasisDegree(1:n) = 1
       END IF

       BodyId = Element % BodyId
       IF (BodyId==0 .AND. ASSOCIATED(Element % BoundaryInfo)) THEN
         Parent => Element % PDefs % LocalParent         
         IF(ASSOCIATED(Parent)) BodyId = Parent % BodyId
         IF( BodyId == 0 ) THEN
           Parent => Element % BoundaryInfo % Left
           IF( ASSOCIATED(Parent)) BodyId = Parent % BodyId
         END IF
         IF(BodyId == 0) THEN
           Parent => Element % BoundaryInfo % Right
           IF( ASSOCIATED(Parent)) BodyId = Parent % BodyId
         END IF
       END IF

       IF (BodyId==0) THEN
         CALL Warn('ElementInfo', 'Element '//I2S(Element % ElementIndex)//' of type '//&
             I2S(Element % TYPE % ElementCode)//' has 0 BodyId, assuming index 1')
         BodyId = 1
       END IF

       ! If running in parallel use global indexing in orienting degrees of freedom
       GIndexes = Element % NodeIndexes
       IF (ASSOCIATED(pSolver % Mesh % ParallelInfo % GlobalDOFs)) &
         GIndexes = pSolver % Mesh % ParallelInfo % GlobalDOFs(GIndexes)

       SerendipityPBasis = Element % PDefs % Serendipity

!------------------------------------------------------------------------------
      SELECT CASE( Element % TYPE % ElementCode ) 
!------------------------------------------------------------------------------

      ! P element code for line element:
      ! --------------------------------
      CASE(202)
         ! Get element p
         p = pSolver % Def_Dofs(2,BodyId,6)
         BDOFs = MAX(GetBubbleDOFs(Element, p), pSolver % Def_Dofs(2,BodyId,5))

         ! Bubbles of line element
         IF (BDOFs > 0) THEN
            ! For boundary element integration check direction
            invert = .FALSE.
            IF ( Element % PDefs % isEdge .AND. &
                    GIndexes(1)>GIndexes(2) ) invert = .TRUE.

            ! For each bubble get the value of basis function
            DO i=1, BDOFs
               IF (q >= SIZE(Basis)) EXIT
               q = q + 1
               
               Basis(q) = LineBubblePBasis(i+1,u,invert)
               dLBasisdx(q,1) = dLineBubblePBasis(i+1,u,invert)
               IF(Compute2ndDerivatives) THEN
                 ddLBasisddx(q,1,1) = ddLineBubblePBasis(i+1,u,invert)
               END IF
               
               ! Polynomial degree of basis function to vector
               IF (degrees) BasisDegree(q) = 1+i
            END DO
         END IF

!------------------------------------------------------------------------------
      ! P element code for triangles:
      CASE(303)
         EDOFs = GetEdgeDOFs(Element, pSolver % Def_Dofs(3,BodyId,6))
         ! Edges of triangle
         IF ( ASSOCIATED( Element % EdgeIndexes ) .AND. EDOFs > 0) THEN
            
            ! For each edge calculate the value of edge basis function
            edges_triangle: DO i=1,3
               Edge => pSolver % Mesh % Edges( Element % EdgeIndexes(i) )

               ! Get local number of edge start and endpoint nodes
               tmp(1:2) = getTriangleEdgeMap(i)
               locali = tmp(1)
               localj = tmp(2)

               ! Invert edge for parity if needed
               invert = .FALSE.
               IF ( GIndexes(locali)>GIndexes(localj) ) invert=.TRUE.

               ! For each edge DOF get the value of p-basis function
               ! NOTE: Edges may not have correct information about the count of DOFs
               !        per edge, so the following would not work:
               !       EDOFs = GetEdgeDOFs(Edge, pSolver % Def_Dofs(2,BodyId,6))
               !
               DO k=1,EDOFs
                  IF (q >= SIZE(Basis)) EXIT edges_triangle
                  q = q + 1
                  
                  ! Value of basis functions for edge=i and i=k+1 by parity
                  Basis(q) = TriangleEdgePBasis(i, k+1, u, v, invert)
                  dLBasisdx(q,1:2) = dTriangleEdgePBasis(i, k+1, u, v, invert)
                  IF(Compute2ndDerivatives) THEN
                    ddLBasisddx(q,1:2,1:2) = ddTriangleEdgePBasis(i,k+1,u,v,invert)
                  END IF
                  
                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO
            END DO edges_triangle
         END IF

         ! Bubbles of p triangle      

         ! Get element p
         p = pSolver % Def_Dofs(3,BodyId,6)
         nb = pSolver % Def_Dofs(3,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)

         IF (BDOFs > 0) THEN
            p = getEffectiveBubbleP(element,p,bdofs)
            
            ! For boundary element direction needs to be calculated
            IF (Element % PDefs % isEdge) THEN
               direction = 0
               ! Get direction of this face (mask for face = boundary element nodes)
               direction(1:3) = getTriangleFaceDirection(Element, [ 1,2,3 ], GIndexes)
            END IF

            bubbles_triangle: DO i = 0,p-3
               DO j = 0,p-i-3
                  IF ( q >= SIZE(Basis) ) EXIT bubbles_triangle
                  q = q + 1

                  ! Get bubble basis functions and their derivatives
                  ! 3d Boundary element has a direction
                  IF (Element % PDefs % isEdge) THEN
                     Basis(q) = TriangleEBubblePBasis(i,j,u,v,direction) 
                     dLBasisdx(q,1:2) = dTriangleEBubblePBasis(i,j,u,v,direction)

                     IF(Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = ddTriangleEBubblePBasis(i,j,u,v,direction)
                     END IF
                  ELSE
                  ! 2d element bubbles have no direction
                     Basis(q) = TriangleBubblePBasis(i,j,u,v) 
                     dLBasisdx(q,1:2) = dTriangleBubblePBasis(i,j,u,v)

                     IF(Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = ddTriangleBubblePBasis(i,j,u,v)
                     END IF
                  END IF
                  
                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 3+i+j
               END DO
            END DO bubbles_triangle
         END IF
!------------------------------------------------------------------------------
      ! P element code for quads:
      CASE(404)
         ! Edges of p quadrilateral
         EDOFs = GetEdgeDOFs(Element, pSolver % Def_Dofs(4,BodyId,6))
         IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
            ! For each edge calculate the values of edge basis functions 
            edges_quad: DO i=1,4
               Edge => pSolver % Mesh % Edges( Element % EdgeIndexes(i) )

               ! Choose correct parity by global edge dofs
               tmp(1:2) = getQuadEdgeMap(i)
               locali = tmp(1)
               localj = tmp(2)
               
               ! Invert parity if needed
               invert = .FALSE.
               IF (GIndexes(locali) > GIndexes(localj)) invert = .TRUE. 

               ! For each DOF in edge calculate the value of p-basis function
               DO k=1,EDOFs
                  IF ( q >= SIZE(Basis) ) EXIT edges_quad
                  q = q + 1

                  ! Get values of basis functions for edge=i and i=k+1 by parity
                  IF (SerendipityPBasis) THEN
                    Basis(q) = SD_QuadEdgePBasis(i,k+1,u,v,invert)
                    ! Get value of derivatives of basis functions
                    dLBasisdx(q,1:2) = SD_dQuadEdgePBasis(i,k+1,u,v,invert)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,1:2,1:2) = SD_ddQuadEdgePBasis(i,k+1,u,v,invert)
                    END IF
                  ELSE
                    Basis(q) = QuadEdgePBasis(i,k+1,u,v,invert)
                    ! Get value of derivatives of basis functions
                    dLBasisdx(q,1:2) = dQuadEdgePBasis(i,k+1,u,v,invert)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,1:2,1:2) = ddQuadEdgePBasis(i,k+1,u,v,invert)
                    END IF
                  END IF
                  
                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO              
            END DO edges_quad
         END IF

         ! Bubbles of p quadrilateral, the number of which may have been defined explicitly or
         ! be determined by the specified degree of approximation. However, we never omit bubbles
         ! which are part of the FE space of the specified degree
  
         ! Get the specified element P:
         p = pSolver % Def_Dofs(4,BodyId,6)
         nb = pSolver % Def_Dofs(4,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb) 

         IF (BDOFs > 0) THEN
            p = getEffectiveBubbleP(element,p,bdofs)

            ! For boundary element direction needs to be calculated
            IF (Element % PDefs % isEdge) THEN
               direction = getSquareFaceDirection(Element, [ 1,2,3,4 ], GIndexes )
            END IF
           
            ! For each bubble calculate the value of p basis function
            ! and its derivatives for index pairs i,j>=2, i+j=4,...,p
            IF(SerendipityPBasis) THEN
              SD_bubbles_quad: DO i=2,p-2
                 DO j=2,p-i
                   IF ( q >= SIZE(Basis) ) EXIT SD_bubbles_quad
                   q = q + 1
                  
                   ! Get values of bubble functions
                   ! 3D boundary elements have a direction
                   IF (Element % PDefs % isEdge) THEN
                     Basis(q) = SD_QuadBubblePBasis(i,j,u,v,direction)
                     dLBasisdx(q,1:2) = SD_dQuadBubblePBasis(i,j,u,v,direction)
                     IF (Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = SD_ddQuadBubblePBasis(i,j,u,v)
                     END IF
                   ELSE
                     ! 2d element bubbles have no direction
                     Basis(q) = SD_QuadBubblePBasis(i,j,u,v)
                     dLBasisdx(q,1:2) = SD_dQuadBubblePBasis(i,j,u,v)
                     IF (Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = SD_ddQuadBubblePBasis(i,j,u,v)
                     END IF
                   END IF
                   ! Polynomial degree of basis function to vector
                   IF (degrees) BasisDegree(q) = i+j
                 END DO
              END DO SD_bubbles_quad
            ELSE
              bubbles_quad: DO i=0,p-2
                 DO j=0,p-2
                   IF ( q >= SIZE(Basis) ) EXIT bubbles_quad
                   q = q + 1
                  
                   ! Get values of bubble functions
                   ! 3D boundary elements have a direction
                   IF (Element % PDefs % isEdge) THEN
                     Basis(q) = QuadBubblePBasis(i,j,u,v,direction)
                     dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v,direction)
                     IF (Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = ddQuadBubblePBasis(i,j,u,v)
                     END IF
                   ELSE
                     ! 2d element bubbles have no direction
                     Basis(q) = QuadBubblePBasis(i,j,u,v)
                     dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v)
                     IF (Compute2ndDerivatives) THEN
                       ddLBasisddx(q,1:2,1:2) = ddQuadBubblePBasis(i,j,u,v)
                     END IF
                   END IF
                   ! Polynomial degree of basis function to vector
                   IF (degrees) BasisDegree(q) = 2+i+j
                 END DO
              END DO bubbles_quad
            END IF
         END IF
!------------------------------------------------------------------------------
      ! P element code for tetrahedra:
      CASE(504)
         p = pSolver % Def_Dofs(5,BodyId,6)
         EDOFs = GetEdgeDOFs(Element, p)
         tetraType = Element % PDefs % TetraType

         ! Edges of p tetrahedron
         IF ( ASSOCIATED( Element % EdgeIndexes ) .AND. EDOFs > 0) THEN   
            ! For each edge i calculate the values of edge functions
            edges_tetrahedron: DO i=1,6
               Edge => pSolver % Mesh % Edges (Element % EdgeIndexes(i))
               
               ! For each edge DOF k calculate the value of edge function
               ! and its derivatives
               DO k=1, EDOFs
                  IF (q >= SIZE(Basis)) EXIT edges_tetrahedron
                  q = q + 1

                  Basis(q) = TetraEdgePBasis(i,k+1,u,v,w,tetraType)
                  dLBasisdx(q,:) = dTetraEdgePBasis(i,k+1,u,v,w,tetraType)
                  IF(Compute2ndDerivatives) THEN
                     ddLBasisddx(q,:,:) = ddTetraEdgePBasis(i,k+1,u,v,w,tetraType)
                  END IF

                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO
            END DO edges_tetrahedron
         END IF

         ! Faces of p tetrahedron
         IF ( ASSOCIATED( Element % FaceIndexes )) THEN
            ! For each face calculate values of face functions
            faces_tetrahedron: DO F=1,4
               Face => pSolver % Mesh % Faces (Element % FaceIndexes(F))

               ! Get face p 
               !p = MAX(pSolver % Def_Dofs(5,BodyId,6), Face % PDefs % P)

               ! Do not solve face DOFs if there is not any
               !IF (GetFaceDOFs(Element, p, F) <= 0) CYCLE

               tmp(1:3) = getTetraFaceMap(F,tetraType)
               direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3), GIndexes )

               ! For each DOF in face calculate values of face function and 
               ! its derivatives for index pairs 
               ! i,j=0,..,p-3, i+j=0,..,p-3
               DO i=0,p-3
                  DO j=0,p-i-3
                     IF (q >= SIZE(Basis)) EXIT faces_tetrahedron
                     q = q + 1 
                   
                     Basis(q) = TetraFacePBasis(F,i,j,u,v,w, tetraType )
                     dLBasisdx(q,:) = dTetraFacePBasis(F,i,j,u,v,w, tetraType )
                     IF(Compute2ndDerivatives) THEN
                       ddLBasisddx(q,:,:) = ddTetraFacePBasis(F,i,j,u,v,w,tetraType )
                     END IF

                     ! Polynomial degree of basis function to vector
                     IF (degrees) BasisDegree(q) = 3+i+j
                  END DO
               END DO
            END DO faces_tetrahedron
         END IF

         ! Bubbles of p tetrahedron
         nb = pSolver % Def_Dofs(5,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb) 
         IF ( BDOFs > 0 ) THEN
            p = getEffectiveBubbleP(element,p,bdofs)

            ! For each bubble DOF calculate the value of bubble function
            ! and its derivatives for index pairs
            ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
            bubbles_tetrahedron: DO i=0,p-4
               DO j=0,p-i-4
                  DO k=0,p-i-j-4
                     IF (q >= SIZE(Basis)) EXIT bubbles_tetrahedron
                     q = q + 1

                     Basis(q) = TetraBubblePBasis(i,j,k,u,v,w)
                     dLBasisdx(q,:) = dTetraBubblePBasis(i,j,k,u,v,w)
                     IF(Compute2ndDerivatives) THEN
                       ddLBasisddx(q,:,:) = ddTetraBubblePBasis(i,j,k,u,v,w)
                     END IF
                     ! Polynomial degree of basis function to vector
                     IF (degrees) BasisDegree(q) = 4+i+j+k
                  END DO
               END DO
            END DO bubbles_tetrahedron
            
         END IF
!------------------------------------------------------------------------------
      ! P element code for pyramids:
      CASE(605)

         IF(SerendipityPBasis) THEN
           CALL Fatal('ElementInfo', 'p-Pyramid not implemented for serendipity scheme, ' // &
                       'please use the full scheme instead.')
         END IF

         ! Edges of P Pyramid
         p = pSolver % Def_Dofs(6,BodyId,6)
         EDOFs = GetEdgeDOFs(Element, p)
         IF (ASSOCIATED( Element % EdgeIndexes ) .AND. EDOFs > 0) THEN
            ! For each edge calculate values of edge functions
            edges_pyramid: DO i=1,8
               Edge => pSolver % Mesh % Edges( Element % EdgeIndexes(i) )

               ! Get local indexes of current edge
               tmp(1:2) = getPyramidEdgeMap(i)
               locali = tmp(1)
               localj = tmp(2)

               ! Determine edge direction
               invert = .FALSE.
               
               ! Invert edge if local first node has greater global index than second one
               IF ( GIndexes(locali) > GIndexes(localj) ) invert = .TRUE.

               ! For each edge DOF k calculate the value of edge function
               ! and its derivatives
               DO k=1,EDOFs
                  IF ( q >= SIZE(Basis) ) EXIT edges_pyramid
                  q = q + 1

                  ! Get values of edge basis functions and their derivatives
                  Basis(q) = PyramidEdgePBasis(i,k+1,u,v,w,invert)
                  dLBasisdx(q,:) = dPyramidEdgePBasis(i,k+1,u,v,w,invert)
                  IF (Compute2ndDerivatives) THEN
                    ddLBasisddx(q,:,:) = ddPyramidEdgePBasis(i,k+1,u,v,w,invert)
                  END IF
                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO
            END DO edges_pyramid
         END IF

         
         ! Faces of P Pyramid
         IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
            ! For each face in pyramid, calculate the values of face functions
            faces_pyramid: DO F=1,5
               Face => pSolver % Mesh % Faces( Element % FaceIndexes(F) )
               
               ! Get face p
               !p = MAX(pSolver % Def_Dofs(6,BodyId,6), Face % PDefs % P) 

               ! Do not solve face dofs, if there is not any
               !IF (GetFaceDOFs(Element, p, F) <= 0) CYCLE
               
               ! Handle triangle and square faces separately
               SELECT CASE(F)
               CASE (1)
                  direction = 0; invert=.FALSE.
                  ! Get global direction vector for enforcing parity
                  tmp(1:4) = getPyramidFaceMap(F)
                  direction(1:4) = getSquareFaceDirection( Element, tmp(1:4), GIndexes )

                  ! For each face calculate the values of functions for index
                  ! pairs i,j=2,..,p-2 i+j=4,..,p

!                DO i=0,p-2
!                   DO j=0,p-i-2
                  DO i=0,p-2
                     DO j=0,p-2
                        IF ( q >= SIZE(Basis) ) EXIT faces_pyramid
                        q = q + 1
                        
                        Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                        dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)
                        IF (Compute2ndDerivatives) THEN
                          ddLBasisddx(q,:,:) = ddPyramidFacePBasis(F,i,j,u,v,w,direction)
                        END IF
                        
                        ! Polynomial degree of basis function to vector
                        IF (degrees) BasisDegree(q) = 2+i+j
                     END DO
                  END DO

               CASE (2,3,4,5)
                  direction = 0
                  ! Get global direction vector for enforcing parity
                  tmp(1:4) = getPyramidFaceMap(F) 
                  direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3), GIndexes )

                  ! For each face calculate the values of functions for index
                  ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                  DO i=0,p-3
                     DO j=0,p-i-3
                        IF ( q >= SIZE(Basis) ) EXIT faces_pyramid
                        q = q + 1

                        Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                        dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)
                        IF (Compute2ndDerivatives) THEN
                          ddLBasisddx(q,:,:) = ddPyramidFacePBasis(F,i,j,u,v,w,direction)
                        END IF

                        ! Polynomial degree of basis function to vector
                        IF (degrees) BasisDegree(q) = 3+i+j
                     END DO
                  END DO
               END SELECT    
            END DO faces_pyramid
         END IF

         ! Bubbles of P Pyramid
         nb = pSolver % Def_Dofs(6,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb) 
         IF ( BDOFs > 0 ) THEN
            p = getEffectiveBubbleP(element,p,bdofs)
  
            ! Calculate the values of bubble functions for indexes
            ! i,j,k=0,..,p-3 i+j+k=0,..,p-3
            bubbles_pyramid: DO i=0,p-3
               DO j=0,p-i-3
                  DO k=0,p-i-j-3
                     IF ( q >= SIZE(Basis)) EXIT bubbles_pyramid
                     q = q + 1

                     Basis(q) = PyramidBubblePBasis(i,j,k,u,v,w)
                     dLBasisdx(q,:) = dPyramidBubblePBasis(i,j,k,u,v,w)
                     IF (Compute2ndDerivatives) THEN
                       ddLBasisddx(q,:,:) = ddPyramidBubblePBasis(i,j,k,u,v,w)
                     END IF
                     
                     ! Polynomial degree of basis function to vector
                     IF (degrees) BasisDegree(q) = 3+i+j+k
                  END DO
               END DO
            END DO bubbles_pyramid
         END IF
         
!------------------------------------------------------------------------------
      ! P element code wedges:
      CASE(706)
         p = pSolver % Def_Dofs(7,BodyId,6)
         EDOFs = GetEdgeDOFs(Element, p)
         ! Edges of P Wedge
         IF (ASSOCIATED( Element % EdgeIndexes ) .AND. EDOFs > 0) THEN
            ! For each edge i calculate the values of edge functions
            edges_prism: DO i=1,9
               Edge => pSolver % Mesh % Edges( Element % EdgeIndexes(i) )
               
               ! Get local indexes of current edge
               tmp(1:2) = getWedgeEdgeMap(i)
               locali = tmp(1)
               localj = tmp(2)

               ! Determine edge direction
               invert = .FALSE.
               ! Invert edge if local first node has greater global index than second one
               IF ( GIndexes(locali) > GIndexes(localj) ) invert = .TRUE.
        
               ! For each edge DOF k calculate the value of edge function
               ! and its derivatives
               DO k=1,EDOFs
                  IF ( q >= SIZE(Basis) ) EXIT edges_prism
                  q = q + 1

                  ! Get values of edge basis functions and their derivatives
                  IF(SerendipityPBasis) THEN
                    Basis(q) = SD_WedgeEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,:) = SD_dWedgeEdgePBasis(i,k+1,u,v,w,invert)
                    IF(Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = SD_ddWedgeEdgePBasis(i,k+1,u,v,w,invert)
                    END IF
                  ELSE
                    Basis(q) = WedgeEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,:) = dWedgeEdgePBasis(i,k+1,u,v,w,invert)
                    IF(Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = ddWedgeEdgePBasis(i,k+1,u,v,w,invert)
                    END IF
                  END IF

                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO
            END DO edges_prism
         END IF

         ! The faces of p-wedge 
         IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
            ! For each face in wedge, calculate the values of face functions
            faces_prism: DO F=1,5
               Face => pSolver % Mesh % Faces( Element % FaceIndexes(F) )

               !p = MAX(pSolver % Def_Dofs(7,BodyId,6), Face % PDefs % P) 

               ! Do not solve face dofs, if there is not any
               !IF (GetFaceDOFs(Element, p, F) <= 0) CYCLE
               
               ! Handle triangle and square faces separately
               SELECT CASE(F)
               CASE (1,2)
                  direction = 0
                  ! Get global direction vector for enforcing parity
                  tmp(1:4) = getWedgeFaceMap(F) 
                  direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3), GIndexes )
                  
                  ! For each face calculate the values of functions for index
                  ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                  DO i=0,p-3
                     DO j=0,p-i-3
                        IF ( q >= SIZE(Basis) ) EXIT faces_prism
                        q = q + 1

                        IF(SerendipityPBasis) THEN
                          Basis(q) = SD_WedgeFacePBasis(F,i,j,u,v,w,direction)
                          dLBasisdx(q,:) = SD_dWedgeFacePBasis(F,i,j,u,v,w,direction)
                          IF(Compute2ndDerivatives) THEN
                             ddLBasisddx(q,:,:) = SD_ddWedgeFacePBasis(F,i,j,u,v,w,direction)
                          END IF
                        ELSE
                          Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)
                          IF(Compute2ndDerivatives) THEN
                             ddLBasisddx(q,:,:) = ddWedgeFacePBasis(F,i,j,u,v,w,direction)
                          END IF
                        END IF

                        ! Polynomial degree of basis function to vector
                        IF (degrees) BasisDegree(q) = 3+i+j
                     END DO
                  END DO
               CASE (3,4,5)
                  direction = 0
                  ! Get global direction vector for enforcing parity
                  invert = .FALSE.
                  tmp(1:4) = getWedgeFaceMap(F)
                  direction(1:4) = getSquareFaceDirection( Element, tmp(1:4), GIndexes )
                  
                  ! First and second node must form a face in upper or lower triangle
                  IF (.NOT. wedgeOrdering(direction)) THEN
                     invert = .TRUE.
                     tmp(1) = direction(2)
                     direction(2) = direction(4)
                     direction(4) = tmp(1)
                  END IF

                  ! For each face calculate values of functions from index
                  ! pairs i,j=2,..,p-2 i+j=4,..,p
                  IF(SerendipityPBasis) THEN
                    DO i=2,p-2
                       DO j=2,p-i
                          IF ( q >= SIZE(Basis) ) EXIT faces_prism
                          q = q + 1

                          IF (.NOT. invert) THEN
                             Basis(q) = SD_WedgeFacePBasis(F,i,j,u,v,w,direction)
                             dLBasisdx(q,:) = SD_dWedgeFacePBasis(F,i,j,u,v,w,direction)
                             IF(Compute2ndDerivatives) THEN
                                ddLBasisddx(q,:,:) = SD_ddWedgeFacePBasis(F,i,j,u,v,w,direction)
                             END IF
                          ELSE
                             Basis(q) = SD_WedgeFacePBasis(F,j,i,u,v,w,direction)
                             dLBasisdx(q,:) = SD_dWedgeFacePBasis(F,j,i,u,v,w,direction)
                             IF(Compute2ndDerivatives) THEN
                                ddLBasisddx(q,:,:) = SD_ddWedgeFacePBasis(F,j,i,u,v,w,direction)
                             END IF
                          END IF
                          ! Polynomial degree of basis function to vector
                          IF (degrees) BasisDegree(q) = i+j
                       END DO
                    END DO
                  ELSE
                    DO i=0,p-2
                       DO j=0,p-2
                          IF ( q >= SIZE(Basis) ) EXIT faces_prism
                          q = q + 1

                          Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)
                          IF(Compute2ndDerivatives) THEN
                             ddLBasisddx(q,:,:) = ddWedgeFacePBasis(F,i,j,u,v,w,direction)
                          END IF
   
                          ! Polynomial degree of basis function to vector
                          IF (degrees) BasisDegree(q) = 2+i+j
                       END DO
                    END DO
                  END IF
               END SELECT
            END DO faces_prism
         END IF

         ! Bubbles of P Wedge
         nb = pSolver % Def_Dofs(7,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb) 
         IF ( BDOFs > 0 ) THEN

            p = getEffectiveBubbleP(element,p,bdofs)
            
            IF(SerendipityPBasis) THEN
              ! For each bubble calculate the value of basis function and its derivative
              ! for index pairs i,j=0,..,p-5 k=2,..,p-3 i+j+k=2,..,p-3
              SD_bubbles_prism: DO i=0,p-5
                 DO j=0,p-5-i
                    DO k=2,p-3-i-j
                       IF ( q >= SIZE(Basis) ) EXIT SD_bubbles_prism
                       q = q + 1

                       Basis(q) = SD_WedgeBubblePBasis(i,j,k,u,v,w)
                       dLBasisdx(q,:) = SD_dWedgeBubblePBasis(i,j,k,u,v,w)
                       IF(Compute2ndDerivatives) THEN
                         ddLBasisddx(q,:,:) = SD_ddWedgeBubblePBasis(i,j,k,u,v,w)
                       END IF

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 3+i+j+k
                    END DO
                 END DO
              END DO SD_bubbles_prism
            ELSE
              bubbles_prism: DO i=0,p-3
                 DO j=0,p-i-3
                    DO k=0,p-2
                       IF ( q >= SIZE(Basis) ) EXIT bubbles_prism
                       q = q + 1

                       Basis(q) = WedgeBubblePBasis(i,j,k,u,v,w)
                       dLBasisdx(q,:) = dWedgeBubblePBasis(i,j,k,u,v,w)
                       IF(Compute2ndDerivatives) THEN
                         ddLBasisddx(q,:,:) = ddWedgeBubblePBasis(i,j,k,u,v,w)
                       END IF

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 2+i+j+k
                    END DO
                 END DO
              END DO bubbles_prism
            END IF
         END IF

!------------------------------------------------------------------------------
      ! P element code for bricks:
      CASE(808) 
         p = pSolver % Def_Dofs(8,BodyId,6)
         EDOFs = GetEdgeDOFs(Element, p)
         ! Edges of P brick
         IF ( ASSOCIATED( Element % EdgeIndexes ) .AND. EDOFs > 0) THEN
            ! For each edge i calculate the values of edge functions 
            edges_brick: DO i=1,12
               Edge => pSolver % Mesh % Edges( Element % EdgeIndexes(i) )
               
               ! Get local indexes of current edge
               tmp(1:2) = getBrickEdgeMap(i)
               locali = tmp(1)
               localj = tmp(2)
               
               ! Determine edge direction
               invert = .FALSE.
               
               ! Invert edge if local first node has greater global index than second one
               IF (GIndexes(locali)>GIndexes(localj)) invert = .TRUE.
               
               ! For each edge DOF k calculate the values of edge function
               ! and its derivatives
               DO k=1,EDOFs
                  IF ( q >= SIZE(Basis) ) EXIT edges_brick
                  q = q + 1

                  ! Get values of edge basis functions and their derivatives
                  IF(SerendipityPBasis) THEN
                    Basis(q) = SD_BrickEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,:) = SD_dBrickEdgePBasis(i,k+1,u,v,w,invert)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = SD_ddBrickEdgePBasis(i,k+1,u,v,w,invert)
                    END IF
                  ELSE
                    Basis(q) = BrickEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,:) = dBrickEdgePBasis(i,k+1,u,v,w,invert)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = ddBrickEdgePBasis(i,k+1,u,v,w,invert)
                    END IF
                  END IF

                  ! Polynomial degree of basis function to vector
                  IF (degrees) BasisDegree(q) = 1+k
               END DO
            END DO edges_brick
         END IF

         ! Faces of P brick
         IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in brick, calculate values of face functions
           faces_brick: DO F=1,6
              Face => pSolver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Get p for face
              !p = MAX(pSolver % Def_Dofs(8,BodyId,6), Face % PDefs % P)
                           
              ! Do not calculate face values if no dofs
              !IF (GetFaceDOFs(Element, p, F)<= 0) CYCLE
               
              ! Generate direction vector for this face
              tmp(1:4) = getBrickFaceMap(F)
              direction(1:4) = getSquareFaceDirection(Element, tmp, GIndexes)

              ! For each face calculate the values of functions for index
              ! pairs i,j=2,..,p-2 i+j=4,..,p
              IF(SerendipityPBasis) THEN
                DO i=2,p-2
                  DO j=2,p-i
                    IF ( q >= SIZE(Basis) ) EXIT faces_brick

                    q = q + 1
                    Basis(q) = SD_BrickFacePBasis(F,i,j,u,v,w,direction)
                    dLBasisdx(q,:) = SD_dBrickFacePBasis(F,i,j,u,v,w,direction)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = SD_ddBrickFacePBasis(F,i,j,u,v,w,direction)
                    END IF
                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j
                   END DO
                END DO
              ELSE
                DO i=0,p-2
                  DO j=0,p-2
                    IF ( q >= SIZE(Basis) ) EXIT faces_brick
  
                    q = q + 1
                    Basis(q) = BrickFacePBasis(F,i,j,u,v,w,direction)
                    dLBasisdx(q,:) = dBrickFacePBasis(F,i,j,u,v,w,direction)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = ddBrickFacePBasis(F,i,j,u,v,w,direction)
                    END IF
                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 2+i+j
                  END DO
                END DO
              END IF
           END DO faces_brick
         END IF

         ! Bubbles of p brick
         nb = pSolver % Def_Dofs(8,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb) 
         IF ( BDOFs > 0 ) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           IF(SerendipityPBasis) THEN
             SD_bubbles_brick: DO i=2,p-4
               DO j=2,p-i-2
                 DO k=2,p-i-j
                    IF ( q >= SIZE(Basis)) EXIT SD_bubbles_brick
                    q = q + 1

                    Basis(q) = SD_BrickBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = SD_dBrickBubblePBasis(i,j,k,u,v,w)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = SD_ddBrickBubblePBasis(i,j,k,u,v,w)
                    END IF
                     
                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j+k
                 END DO
               END DO
             END DO SD_bubbles_brick
           ELSE 
             bubbles_brick: DO i=0,p-2
               DO j=0,p-2
                 DO k=0,p-2
                   IF ( q >= SIZE(Basis)) EXIT bubbles_brick
                   q = q + 1

                    Basis(q) = BrickBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dBrickBubblePBasis(i,j,k,u,v,w)
                    IF (Compute2ndDerivatives) THEN
                      ddLBasisddx(q,:,:) = ddBrickBubblePBasis(i,j,k,u,v,w)
                    END IF
                     
                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 2+i+j+k
                 END DO
               END DO
             END DO bubbles_brick
           END IF
         END IF

      END SELECT
      END IF
     END SUBROUTINE EvalPElementBasis

!------------------------------------------------------------------------------
   
   ! SUBROUTINE ElementInfoVec_InitWork(m, n)
   !   IMPLICIT NONE

   !   INTEGER, INTENT(IN) :: m, n
   !   INTEGER :: allocstat

   !   allocstat = 0
   !   IF (.NOT. ALLOCATED(BasisWrk)) THEN
   !     ALLOCATE(BasisWrk(m,n), &
   !             dBasisdxWrk(m,n,3), &
   !             LtoGMapsWrk(m,3,3), &
   !             DetJWrk(m), &
   !             uWrk(m), vWrk(m), wWrk(m), STAT=allocstat)
   !   ELSE IF (SIZE(BasisWrk,1) /= m .OR. SIZE(BasisWrk,2) /= n) THEN
   !     DEALLOCATE(BasisWrk, dBasisdxWrk, LtoGMapsWrk, DetJWrk, uWrk, vWrk, wWrk)
   !     ALLOCATE(BasisWrk(m,n), &
   !             dBasisdxWrk(m,n,3), &
   !             LtoGMapsWrk(m,3,3), &
   !             DetJWrk(m), &
   !             uWrk(m), vWrk(m), wWrk(m), STAT=allocstat)
   !   END IF

   !   ! Check memory allocation status
   !   IF (allocstat /= 0) THEN
   !     CALL Error('ElementInfo_InitWork','Storage allocation for local element basis failed')
   !   END IF
   ! END SUBROUTINE ElementInfoVec_InitWork

   ! SUBROUTINE ElementInfoVec_FreeWork()
   !   IMPLICIT NONE

   !   IF (ALLOCATED(BasisWrk)) THEN
   !     DEALLOCATE(BasisWrk, dBasisdxWrk, LtoGMapsWrk, DetJWrk, uWrk, vWrk, wWrk)
   !   END IF
   ! END SUBROUTINE ElementInfoVec_FreeWork

!
!------------------------------------------------------------------------------
   FUNCTION ElementInfoVec( Element, Nodes, nc, u, v, w, detJ, nbmax, &
               Basis, dBasisdx, USolver ) RESULT(retval)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element    !< Element structure
     TYPE(Nodes_t)   :: Nodes              !< Element nodal coordinates.
     INTEGER, INTENT(IN) :: nc             !< Number of local coordinates to compute values of the basis function
     REAL(KIND=dp), POINTER CONTIG :: u(:)  !< 1st local coordinates at which to calculate the basis function.
     REAL(KIND=dp), POINTER CONTIG :: v(:)  !< 2nd local coordinates.
     REAL(KIND=dp), POINTER CONTIG :: w(:)  !< 3rd local coordinates.
     REAL(KIND=dp) CONTIG, INTENT(OUT) :: detJ(:) !< Square roots of determinants of element coordinate system metric at coordinates
     INTEGER, INTENT(IN) :: nbmax          !< Maximum number of basis functions to compute
     REAL(KIND=dp) CONTIG :: Basis(:,:)    !< Basis function values at (u,v,w)
     REAL(KIND=dp) CONTIG, OPTIONAL :: dBasisdx(:,:,:)    !< Global first derivatives of basis functions at (u,v,w)
     TYPE(Solver_t), TARGET, OPTIONAL :: USolver
     LOGICAL :: retval                             !< If .FALSE. element is degenerate. or if local storage allocation fails

     ! Internal work arrays (always needed)
     REAL(KIND=dp) :: uWrk(VECTOR_BLOCK_LENGTH), vWrk(VECTOR_BLOCK_LENGTH), wWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: BasisWrk(VECTOR_BLOCK_LENGTH,nbmax)
     REAL(KIND=dp) :: dBasisdxWrk(VECTOR_BLOCK_LENGTH,nbmax,3)
     REAL(KIND=dp) :: DetJWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: LtoGMapsWrk(VECTOR_BLOCK_LENGTH,3,3)

     TYPE(Solver_t), POINTER :: pSolver
     
     INTEGER :: i, l, n, dim, cdim, ll, ncl, lln, nbp
     LOGICAL :: elem
!DIR$ ATTRIBUTES ALIGN:64::uWrk, vWrk, wWrk, BasisWrk, dBasisdxWrk, DetJWrk, LtoGMapsWrk
     
     !------------------------------------------------------------------------------
     ! Special case, Element: POINT
     IF (Element % TYPE % ElementCODE == 101) THEN
       DetJ(1:nc) = REAL(1, dp)
       Basis(1:nc,1) = REAL(1, dp)
       IF (PRESENT(dBasisdx)) THEN
         DO i=1,nc
           dBasisdx(i,1,1) = REAL(0, dp)
         END DO
       END IF
       retval = .TRUE.
       RETURN
     END IF
     
     ! Set up workspace arrays 
     ! CALL ElementInfoVec_InitWork(VECTOR_BLOCK_LENGTH, nbmax)
     IF ( nbmax < Element % TYPE % NumberOfNodes ) THEN
       CALL Fatal('ElementInfoVec','Not enough storage to compute local element basis')
     END IF

     IF(PRESENT(dBasisdx))  &
       dBasisdx = 0._dp ! avoid uninitialized stuff depending on coordinate dimension...

     IF( isActivePelement(Element) ) THEN
       retval =  ElementInfoVec_ComputePElementBasis(Element,Nodes,nc,u,v,w,detJ,nbmax,Basis,&
             uWrk,vWrk,wWrk,BasisWrk,dBasisdxWrk,DetJWrk,LtoGmapsWrk,dBasisdx,USolver)
     ELSE
       retval = .TRUE.
       n    = Element % TYPE % NumberOfNodes
       dim  = Element % TYPE % DIMENSION
       cdim = CoordinateSystemDimension()

       IF (nc <= VECTOR_SMALL_THRESH) THEN
         DO l=1,nc
           CALL NodalBasisFunctions(n, BasisWrk(l,:), Element, u(l), v(l), w(l))
           CALL NodalFirstDerivatives(n, dBasisdxWrk(l,:,:), Element, u(l), v(l), w(l))
         END DO
         elem = ElementMetricVec(Element, Nodes, nc, n, DetJ(1:nc), &
               nbmax, dBasisdxWrk, LtoGMapsWrk)
         IF (.NOT. elem) THEN; retval = .FALSE.; RETURN; END IF
         Basis(1:nc, 1:n) = BasisWrk(1:nc, 1:n)
         IF (PRESENT(dBasisdx)) &
           CALL ElementInfoVec_ElementBasisToGlobal(nc, n, nbmax, dBasisdxWrk, &
               dim, cdim, LtoGMapsWrk, 1, dBasisdx)
         RETURN
       END IF

       DO ll=1,nc,VECTOR_BLOCK_LENGTH
         lln = MIN(ll+VECTOR_BLOCK_LENGTH-1,nc)
         ncl = lln-ll+1

         ! Block copy input
         uWrk(1:ncl) = u(ll:lln)
         IF (cdim > 1) THEN
           vWrk(1:ncl) = v(ll:lln)
         END IF
         IF (cdim > 2) THEN
           wWrk(1:ncl) = w(ll:lln)
         END IF

         ! H1Basis vectorized dispatch: only safe for minimum-node (linear) elements
         ! (type >= nodes means no higher-order nodes to fill).
         ! Higher-order elements fall through to the scalar loop.
         IF (Element % TYPE % ElementCode/100 >= MODULO(Element % TYPE % ElementCode, 100)) THEN
           SELECT CASE(Element % TYPE % ElementCode / 100)
           CASE(2)
             nbp = 0; CALL H1Basis_LineNodal(ncl, uWrk, nbmax, BasisWrk, nbp)
             nbp = 0; CALL H1Basis_dLineNodal(ncl, uWrk, nbmax, dBasisdxWrk, nbp)
           CASE(3)
             CALL H1Basis_TriangleNodal(ncl, uWrk, vWrk, nbmax, BasisWrk)
             CALL H1Basis_dTriangleNodal(ncl, uWrk, vWrk, nbmax, dBasisdxWrk)
           CASE(4)
             nbp = 0; CALL H1Basis_QuadNodal(ncl, uWrk, vWrk, nbmax, BasisWrk, nbp)
             nbp = 0; CALL H1Basis_dQuadNodal(ncl, uWrk, vWrk, nbmax, dBasisdxWrk, nbp)
           CASE(5)
             CALL H1Basis_TetraNodal(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk)
             CALL H1Basis_dTetraNodal(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk)
           CASE(7)
             CALL H1Basis_WedgeNodal(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk)
             CALL H1Basis_dWedgeNodal(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk)
           CASE(8)
             nbp = 0; CALL H1Basis_BrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
             nbp = 0; CALL H1Basis_dBrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbp)
           CASE DEFAULT
             DO l=1,ncl
               CALL NodalBasisFunctions(n, BasisWrk(l,:), element, uWrk(l), vWrk(l), wWrk(l))
               CALL NodalFirstDerivatives(n, dBasisdxWrk(l,:,:), element, uWrk(l), vWrk(l), wWrk(l))
             END DO
           END SELECT
         ELSE
           DO l=1,ncl
             CALL NodalBasisFunctions(n, BasisWrk(l,:), element, uWrk(l), vWrk(l), wWrk(l))
             CALL NodalFirstDerivatives(n, dBasisdxWrk(l,:,:), element, uWrk(l), vWrk(l), wWrk(l))
           END DO
         END IF
         Basis(ll:lln, 1:n) = BasisWrk(1:ncl, 1:n)

         ! Element (contravariant) metric and square root of determinant
         !--------------------------------------------------------------
         elem = ElementMetricVec( Element, Nodes, ncl, n, DetJWrk, &
                nbmax, dBasisdxWrk, LtoGMapsWrk )

         IF (.NOT. elem) THEN
           retval = .FALSE.
           RETURN
           END IF

         !_ELMER_OMP_SIMD
         DO i=1,ncl
           DetJ(i+ll-1)=DetJWrk(i)
         END DO

         ! Get global basis functions
         !--------------------------------------------------------------
         ! First derivatives
         IF (PRESENT(dBasisdx)) THEN
!DIR$ FORCEINLINE
           CALL ElementInfoVec_ElementBasisToGlobal(ncl, n, nbmax, dBasisdxWrk, dim, cdim, LtoGMapsWrk, ll, dBasisdx)
         END IF
       END DO
     END IF
   END FUNCTION ElementInfoVec
     
   FUNCTION ElementInfoVec_ComputePElementBasis(Element, Nodes, nc, u, v, w, DetJ, nbmax, Basis, &
      uWrk, vWrk, wWrk, BasisWrk, dBasisdxWrk, DetJWrk, LtoGmapsWrk, dBasisdx, USolver) RESULT(retval)
     IMPLICIT NONE
     TYPE(Element_t), TARGET :: Element    !< Element structure
     TYPE(Nodes_t)   :: Nodes              !< Element nodal coordinates.
     INTEGER, INTENT(IN) :: nc             !< Number of local coordinates to compute values of the basis function
     REAL(KIND=dp), POINTER CONTIG :: u(:)  !< 1st local coordinates at which to calculate the basis function.
     REAL(KIND=dp), POINTER CONTIG :: v(:)  !< 2nd local coordinates.
     REAL(KIND=dp), POINTER CONTIG :: w(:)  !< 3rd local coordinates.
     REAL(KIND=dp) CONTIG, INTENT(OUT) :: detJ(:) !< Square roots of determinants of element coordinate system metric at coordinates
     INTEGER, INTENT(IN) :: nbmax          !< Maximum number of basis functions to compute
     REAL(KIND=dp) CONTIG :: Basis(:,:)    !< Basis function values at (u,v,w)
     ! Internal work arrays
     REAL(KIND=dp) :: uWrk(VECTOR_BLOCK_LENGTH), vWrk(VECTOR_BLOCK_LENGTH), wWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: BasisWrk(VECTOR_BLOCK_LENGTH,nbmax)
     REAL(KIND=dp) :: dBasisdxWrk(VECTOR_BLOCK_LENGTH,nbmax,3)
     REAL(KIND=dp) :: DetJWrk(VECTOR_BLOCK_LENGTH)
     REAL(KIND=dp) :: LtoGMapsWrk(VECTOR_BLOCK_LENGTH,3,3)
     REAL(KIND=dp) CONTIG, OPTIONAL :: dBasisdx(:,:,:)    !< Global first derivatives of basis functions at (u,v,w)
     TYPE(Solver_t), TARGET, OPTIONAL :: USolver
     LOGICAL :: retval                    !< If .FALSE. element is degenerate. or if local storage allocation fails


     !------------------------------------------------------------------------------
     !    Local variables
     !------------------------------------------------------------------------------
     INTEGER :: EdgeDegree(H1Basis_MaxPElementEdges), &
           FaceDegree(H1Basis_MaxPElementFaces), &
           EdgeDirection(H1Basis_MaxPElementEdgeNodes,H1Basis_MaxPElementEdges), &
           FaceDirection(H1Basis_MaxPElementFaceNodes,H1Basis_MaxPElementFaces)

     INTEGER :: cdim, dim, i, j, k, l, ll, lln, ncl, ip, n, p, nb, bdofs, &
           nbp, nbq, nbdxp, allocstat, ncpad, EdgeMaxDegree, FaceMaxDegree, BodyId

     TYPE(Solver_t), POINTER :: pSolver
     TYPE(Element_t), POINTER :: Parent

     LOGICAL :: invertBubble, elem, SerendipityPBasis
 
!DIR$ ATTRIBUTES ALIGN:64::EdgeDegree, FaceDegree
!DIR$ ATTRIBUTES ALIGN:64::EdgeDirection, FaceDirection
!DIR$ ASSUME_ALIGNED uWrk:64, vWrk:64, wWrk:64, BasisWrk:64, dBasisdxWrk:64, DetJWrk:64, LtoGMapsWrk:64

     IF( PRESENT( USolver ) ) THEN
       pSolver => USolver
     ELSE
       pSolver => CurrentModel % Solver
     END IF

     BodyId = Element % BodyId
     IF( isActivePElement(Element)) THEN
       IF (BodyId==0 .AND. ASSOCIATED(Element % BoundaryInfo)) THEN
         Parent => Element % PDefs % LocalParent
         IF(ASSOCIATED(Parent)) BodyId = Parent % BodyId
       END IF
       SerendipityPBasis = Element % PDefs % Serendipity
     ELSE
       IF (BodyId==0 .AND. ASSOCIATED(Element % BoundaryInfo)) THEN
         Parent => Element % BoundaryInfo % Left
         IF(ASSOCIATED(Parent)) BodyId = Parent % BodyId
       END IF
     END IF

     IF (BodyId==0) THEN
       CALL Warn('ElementInfoVec', 'Element '//I2S(Element % ElementIndex)//' of type '//&
           I2S(Element % TYPE % ElementCode)//' has 0 BodyId, assuming index 1')
       BodyId = 1
     END IF

     retval = .TRUE.
     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     dBasisdxWrk = 0._dp ! avoid uninitialized stuff depending on coordinate dimension...

     ! Block the computation for large values of input points
     DO ll=1,nc,VECTOR_BLOCK_LENGTH
       lln = MIN(ll+VECTOR_BLOCK_LENGTH-1,nc)
       ncl = lln-ll+1

       ! Set number of computed basis functions
       nbp = 0
       nbdxp = 0

       ! Block copy input
       uWrk(1:ncl) = u(ll:lln)
       IF (cdim > 1) THEN
         vWrk(1:ncl) = v(ll:lln)
       END IF
       IF (cdim > 2) THEN
         wWrk(1:ncl) = w(ll:lln)
       END IF

       ! Compute local p element basis
       SELECT CASE (Element % Type % ElementCode)
         ! Element: LINE
       CASE (202)
         ! Compute nodal basis
         CALL H1Basis_LineNodal(ncl, uWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dLineNodal(ncl, uWrk, nbmax, dBasisdxWrk, nbdxp)

         ! Element bubble functions
         p = pSolver % Def_Dofs(2,BodyId,6)
         nb = pSolver % Def_Dofs(2,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)

         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           ! For first round of blocked loop, compute edge direction
           IF (ll==1) THEN
             IF (Element % PDefs % isEdge .AND. &
                   Element % NodeIndexes(1)> Element % NodeIndexes(2)) THEN
               invertBubble = .TRUE.
             ELSE
               invertBubble = .FALSE.
             END IF
           END IF

           CALL H1Basis_LineBubbleP(ncl, uWrk, P, nbmax, BasisWrk, nbp, invertBubble)
           CALL H1Basis_dLineBubbleP(ncl, uWrk, P, nbmax, dBasisdxWrk, nbdxp, invertBubble)
         END IF

         ! Element: TRIANGLE
       CASE (303)
         ! Compute nodal basis
         CALL H1Basis_TriangleNodalP(ncl, uWrk, vWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dTriangleNodalP(ncl, uWrk, vWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes)) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree>1 ) THEN
             nbq = nbp + SUM(EdgeDegree(1:3)-1)
             IF(nbmax >= nbq ) THEN
               CALL H1Basis_TriangleEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, BasisWrk, &
                     nbp, EdgeDirection)
               CALL H1Basis_dTriangleEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, dBasisdxWrk, &
                     nbdxp, EdgeDirection)
             END IF
           END IF
         END IF

         ! Element bubble functions
         p = pSolver % Def_Dofs(3,BodyId,6)
         nb = pSolver % Def_Dofs(3,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)

         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             IF (Element % PDefs % isEdge) THEN
               ! Get 2D face direction
               CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                     1, Element % NodeIndexes, FaceDirection)
             END IF
           END IF
           IF (Element % PDefs % isEdge) THEN
             CALL H1Basis_TriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp, &
                   FaceDirection(1:3,1))
             CALL H1Basis_dTriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp, &
                   FaceDirection(1:3,1))
           ELSE
             CALL H1Basis_TriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_dTriangleBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp)
           END IF
         END IF

         ! QUADRILATERAL
       CASE (404)
         ! Compute nodal basis
         CALL H1Basis_QuadNodal(ncl, uWrk, vWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dQuadNodal(ncl, uWrk, vWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             nbq = nbp + SUM(EdgeDegree(1:4)-1)
             IF(nbmax >= nbq) THEN
               IF(SerendipityPBasis) THEN
                 CALL H1Basis_SD_QuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_SD_dQuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
               ELSE
                 CALL H1Basis_QuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_dQuadEdgeP(ncl, uWrk, vWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
               END IF
             END IF
           END IF
         END IF

         ! Element bubble functions
         p = pSolver % Def_Dofs(4,BodyId,6)
         nb = pSolver % Def_Dofs(4,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)

         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           IF(nbmax-nbp<getBubbleDOFs(Element,p)) THEN
             IF(SerendipityPBasis) THEN
               CALL Fatal("ElementInfoVec", &
                 "Not enough space for storing bubble basis, check your #bubbles: i*(i-1)/2 (0,1,3,6,10,15,...)")
             ELSE
               CALL Fatal("ElementInfoVec", &
                 "Not enough space for storing bubble basis, check your #bubbles: i^2 (0,1,4,9,16,25,...)")
             END IF
           END IF

           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             IF (Element % PDefs % isEdge) THEN
               ! Get 2D face direction
               CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                     1,  Element % NodeIndexes,  FaceDirection)
             END IF
           END IF

           IF (Element % PDefs % isEdge) THEN
             IF(SerendipityPBasis) THEN
               CALL H1Basis_SD_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp, &
                     FaceDirection(1:4,1))
               CALL H1Basis_SD_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp, &
                     FaceDirection(1:4,1))
             ELSE
               CALL H1Basis_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp, &
                     FaceDirection(1:4,1))
               CALL H1Basis_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp, &
                     FaceDirection(1:4,1))
             END IF
           ELSE
             IF(SerendipityPBasis) THEN
               CALL H1Basis_SD_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp)
               CALL H1Basis_SD_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp)
             ELSE
               CALL H1Basis_QuadBubbleP(ncl, uWrk, vWrk, P, nbmax, BasisWrk, nbp)
               CALL H1Basis_dQuadBubbleP(ncl, uWrk, vWrk, P, nbmax, dBasisdxWrk, nbdxp)
             END IF
           END IF
         END IF

         ! TETRAHEDRON
       CASE (504)
         ! Compute nodal basis
         CALL H1Basis_TetraNodalP(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)

         ! Compute local first derivatives
         CALL H1Basis_dTetraNodalP(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             ! Get polynomial degree of each edge
             EdgeMaxDegree = 0
             IF( CurrentModel % Solver % Mesh % MaxEdgeDofs == 0 ) THEN
               CONTINUE             
             ELSE
               DO i=1,6
                 EdgeDegree(i) = CurrentModel % Solver % &
                       Mesh % Edges( Element % EdgeIndexes(i) ) % BDOFs + 1
                 EdgeMaxDegree = MAX(EdgeDegree(i),EdgeMaxDegree)
               END DO
             END IF

             ! Tetrahedral directions are enforced by tetra element types
             IF (EdgeMaxDegree > 1) THEN
               CALL H1Basis_GetTetraEdgeDirection(Element % PDefs % TetraType, EdgeDirection)
             END IF
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             nbq = nbp + SUM(EdgeDegree(1:6)-1)
             IF(nbmax >= nbq) THEN
               CALL H1Basis_TetraEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                     EdgeDirection)
               CALL H1Basis_dTetraEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                     EdgeDirection)
             END IF
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             ! Get polynomial degree of each face
             FaceMaxDegree = 0

             IF( CurrentModel % Solver % Mesh % MaxFaceDofs == 0 ) THEN
               CONTINUE             
             ELSE IF (CurrentModel % Solver % Mesh % MinFaceDOFs == &
                   CurrentModel % Solver % Mesh % MaxFaceDOFs) THEN
               FaceMaxDegree = CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(1) ) % PDefs % P
               FaceDegree(1:Element % Type % NumberOfFaces) = FaceMaxDegree
             ELSE
               DO i=1,4
                 IF (CurrentModel % Solver % Mesh % &
                       Faces( Element % FaceIndexes(i) ) % BDOFs /= 0) THEN
                   FaceDegree(i) = CurrentModel % Solver % Mesh % &
                         Faces( Element % FaceIndexes(i) ) % PDefs % P
                   FaceMaxDegree = MAX(FaceDegree(i), FaceMaxDegree)
                 ELSE
                   FaceDegree(i) = 0
                 END IF
               END DO
             END IF

             IF (FaceMaxDegree > 1) THEN
               CALL H1Basis_GetTetraFaceDirection(Element % PDefs % TetraType, FaceDirection)
             END IF
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree>1 ) THEN
             nbq = nbp
             DO i=1,4
               DO j=0,FaceDegree(i)
                  nbq = nbq + MAX(FaceDegree(i)-j-2,0)
               END DO
             END DO
  
             IF (nbmax >= nbq ) THEN
               CALL H1Basis_TetraFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                     FaceDirection)
               CALL H1Basis_dTetraFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                     FaceDirection)
             END IF
           END IF
         END IF

         ! Element bubble functions
         p = pSolver % Def_Dofs(5,BodyId,6)
         nb = pSolver % Def_Dofs(5,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)
           CALL H1Basis_TetraBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
           CALL H1Basis_dTetraBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
         END IF

       ! Pyramid
       CASE (605)
         IF(SerendipityPBasis) THEN
           CALL Fatal('ElementInfoVec', 'p-Pyramid not available for serendipity scheme, ' // &
                  'please use full polynomial scheme instead.' )
         END IF

         ! Compute nodal basis
         CALL H1Basis_PyramidNodalP(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dPYramidNodalP(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1)THEN
             nbq = nbp+SUM(EdgeDegree(1:8)-1)
             IF(nbmax >= nbq) THEN
               CALL H1Basis_PyramidEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                     EdgeDirection)

               CALL H1Basis_dPyramidEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                     EdgeDirection)
             END IF
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             CALL GetElementMeshFaceInfo(CurrentModel % Solver % Mesh, &
                   Element, FaceDegree, FaceDirection, FaceMaxDegree)
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1 ) THEN
             nbq = nbp
             ! Square faces
             DO i=1,1
               DO j=0,FaceDegree(i)-2
                 nbq = nbq + MAX(FaceDegree(i)-1,0)
               END DO
             END DO

             ! Triangle faces
             DO i=2,5
               DO j=0,FaceDegree(i)-3
                 nbq = nbq + MAX(FaceDegree(i)-j-2,0)
               END DO
             END DO
             
             IF(nbmax >= nbq) THEN
               CALL H1Basis_PyramidFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                     FaceDirection)
               CALL H1Basis_dPyramidFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                     FaceDirection)
             END IF
           END IF
         END IF

         ! Element bubble functions
         p = pSolver % Def_Dofs(6,BodyId,6)
         nb = pSolver % Def_Dofs(6,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           CALL H1Basis_PyramidBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
           CALL H1Basis_dPyramidBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
         END IF


         ! WEDGE
       CASE (706)
         ! Compute nodal basis
         CALL H1Basis_WedgeNodalP(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dWedgeNodalP(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1)THEN
             nbq = nbp+SUM(EdgeDegree(1:9)-1)
             IF(nbmax >= nbq) THEN
               IF(SerendipityPBasis) THEN
                 CALL H1Basis_SD_WedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_SD_dWedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
               ELSE
                 CALL H1Basis_WedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_dWedgeEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
               END IF
             END IF
           END IF
         END IF

         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             CALL GetElementMeshFaceInfo(CurrentModel % Solver % Mesh, &
                   Element, FaceDegree, FaceDirection, FaceMaxDegree)
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1 ) THEN
             nbq = nbp
             ! Triangle faces
             DO i=1,2
               DO j=0,FaceDegree(i)-3
                 nbq = nbq + MAX(FaceDegree(i)-j-2,0)
               END DO
             END DO
             ! Square faces
             DO i=3,5
               IF(SerendipityPBasis) THEN
                 DO j=2,FaceDegree(i)-2
                   nbq = nbq + MAX(FaceDegree(i)-j-1,0)
                 END DO
               ELSE
                 DO j=0,FaceDegree(i)-2
                   nbq = nbq + MAX(FaceDegree(i)-1,0)
                 END DO
               END IF
             END DO
             
             IF(nbmax >= nbq) THEN
               IF(SerendipityPBasis) THEN
                 CALL H1Basis_SD_WedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                       FaceDirection)
                 CALL H1Basis_SD_dWedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                       FaceDirection)
               ELSE
                 CALL H1Basis_WedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                       FaceDirection)
                 CALL H1Basis_dWedgeFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                       FaceDirection)
               END IF
             END IF
           END IF
         END IF

         ! Element bubble functions
         p = pSolver % Def_Dofs(7,BodyId,6)
         nb = pSolver % Def_Dofs(7,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)
           IF(SerendipityPBasis) THEN
             CALL H1Basis_SD_WedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_SD_dWedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
           ELSE
             CALL H1Basis_WedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_dWedgeBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
           END IF
         END IF

         ! HEXAHEDRON
       CASE (808)
         ! Compute local basis
         CALL H1Basis_BrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, BasisWrk, nbp)
         ! Compute local first derivatives
         CALL H1Basis_dBrickNodal(ncl, uWrk, vWrk, wWrk, nbmax, dBasisdxWrk, nbdxp)

         IF (ASSOCIATED( Element % EdgeIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! edge directions
           IF (ll==1) THEN
             CALL GetElementMeshEdgeInfo(CurrentModel % Solver % Mesh, &
                   Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
           END IF

           ! Compute basis function values
           IF (EdgeMaxDegree > 1) THEN
             nbq = nbp + SUM(EdgeDegree(1:12)-1)
             IF(nbmax >= nbq) THEN
               IF(SerendipityPBasis) THEN
                 CALL H1Basis_SD_BrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_SD_dBrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
               ELSE
                 CALL H1Basis_BrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, BasisWrk, nbp, &
                       EdgeDirection)
                 CALL H1Basis_dBrickEdgeP(ncl, uWrk, vWrk, wWrk, EdgeDegree, nbmax, dBasisdxWrk, nbdxp, &
                       EdgeDirection)
              END IF
             END IF
           END IF
         END IF


         IF (ASSOCIATED( Element % FaceIndexes )) THEN
           ! For first round of blocked loop, compute polynomial degrees and 
           ! face directions
           IF (ll==1) THEN
             CALL GetElementMeshFaceInfo(CurrentModel % Solver % Mesh, &
                   Element, FaceDegree, FaceDirection, FaceMaxDegree)
           END IF

           ! Compute basis function values
           IF (FaceMaxDegree > 1) THEN
             nbq = nbp
             DO i=1,6
               DO j=2,FaceDegree(i)
                 nbq = nbq + MAX(FaceDegree(i)-j-1,0)
               END DO
             END DO

             IF(nbmax >= nbq) THEN
               IF(SerendipityPBasis) THEN
                 CALL H1Basis_SD_BrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                       FaceDirection)
                 CALL H1Basis_SD_dBrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                       FaceDirection)
               ELSE
                 CALL H1Basis_BrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, BasisWrk, nbp, &
                       FaceDirection)
                 CALL H1Basis_dBrickFaceP(ncl, uWrk, vWrk, wWrk, FaceDegree, nbmax, dBasisdxWrk, nbdxp, &
                       FaceDirection)
               END IF
             END IF
           END IF
         END IF

         
         ! Element bubble functions
         p = pSolver % Def_Dofs(8,BodyId,6)
         nb = pSolver % Def_Dofs(8,BodyId,5)
         BDOFs = MAX(GetBubbleDOFs(Element, p), nb)
         IF (BDOFs > 0) THEN
           p = getEffectiveBubbleP(element,p,bdofs)

           IF(nbmax-nbp<getBubbleDOFs(Element,p)) THEN
             IF(SerendipityPBasis) THEN
               CALL Fatal("ElementInfoVec", &
                 "Not enough space for storing bubble basis, check your #bubbles: i*(i-1)*(i-1)/2 (0,1,4,10,16,...)")
             ELSE
               CALL Fatal("ElementInfoVec", &
                 "Not enough space for storing bubble basis, check your #bubbles: i^3: (0,1,8,27,64,...)")
             END IF
           END IF

           IF(SerendipityPBasis) THEN
             CALL H1Basis_SD_BrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_SD_dBrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
           ELSE
             CALL H1Basis_BrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, BasisWrk, nbp)
             CALL H1Basis_dBrickBubbleP(ncl, uWrk, vWrk, wWrk, P, nbmax, dBasisdxWrk, nbdxp)
           END IF
         END IF

         
       CASE DEFAULT
         WRITE( Message, '(a,i4,a)' ) 'Vectorized basis for element: ', &
               Element % TYPE % ElementCode, ' not implemented.'
         CALL Error( 'ElementInfoVec', Message )
         CALL Fatal( 'ElementInfoVec', 'ElementInfoVec is still does not include pyramids.' )
       END SELECT

       ! Copy basis function values to global array
       Basis(ll:lln, 1:nbp) = BasisWrk(1:ncl, 1:nbp)

       !--------------------------------------------------------------
       ! Element (contravariant) metric and square root of determinant
       !--------------------------------------------------------------
       elem = ElementMetricVec( Element, Nodes, ncl, nbp, DetJWrk, &
             nbmax, dBasisdxWrk, LtoGMapsWrk )
       IF (.NOT. elem) THEN
         retval = .FALSE.
         RETURN
       END IF

       !_ELMER_OMP_SIMD
       DO i=1,ncl
         DetJ(i+ll-1)=DetJWrk(i)
       END DO

       ! Get global basis functions
       !--------------------------------------------------------------
       ! First derivatives
       IF (PRESENT(dBasisdx)) THEN
!DIR$ FORCEINLINE
         CALL ElementInfoVec_ElementBasisToGlobal(ncl, nbp, nbmax, dBasisdxWrk, dim, cdim, LtoGMapsWrk, ll, dBasisdx)
       END IF
     END DO ! Block over Gauss points

  CONTAINS
   
     SUBROUTINE GetElementMeshEdgeInfo(Mesh, Element, EdgeDegree, EdgeDirection, EdgeMaxDegree)
       IMPLICIT NONE
       
       TYPE(Mesh_t), INTENT(IN) :: Mesh
       TYPE(Element_t), INTENT(IN) :: Element
       INTEGER, INTENT(OUT) :: EdgeDegree(H1Basis_MaxPElementEdges), &
               EdgeDirection(H1Basis_MaxPElementEdgeNodes,H1Basis_MaxPElementEdges)
       INTEGER, INTENT(OUT) :: EdgeMaxDegree
       INTEGER :: i

       EdgeMaxDegree = 0

       IF( Mesh % MaxEdgeDofs == 0 ) THEN
         CONTINUE             

       ELSE IF (Mesh % MinEdgeDOFs == Mesh % MaxEdgeDOFs) THEN
          EdgeDegree(1:Element % Type % NumberOfEdges) = Mesh % MaxEdgeDOFs + 1
          EdgeMaxDegree = Mesh % MaxEdgeDOFs + 1
       ELSE
       ! Get polynomial degree of each edge separately
!DIR$ LOOP COUNT MAX=12
          DO i=1,Element % Type % NumberOfEdges
             EdgeDegree(i) = Mesh % Edges( Element % EdgeIndexes(i) ) % BDOFs + 1
             EdgeMaxDegree = MAX(EdgeDegree(i), EdgeMaxDegree)
          END DO
       END IF

       ! Get edge directions if needed
       IF (EdgeMaxDegree > 1) THEN
         CALL H1Basis_GetEdgeDirection(Element % Type % ElementCode, &
                                       Element % Type % NumberOfEdges, &
                                       Element % NodeIndexes, &
                                       EdgeDirection)
       END IF
     END SUBROUTINE GetElementMeshEdgeInfo
     
     SUBROUTINE GetElementMeshFaceInfo(Mesh, Element, FaceDegree, FaceDirection, FaceMaxDegree)
       IMPLICIT NONE
       
       TYPE(Mesh_t), INTENT(IN) :: Mesh
       TYPE(Element_t), INTENT(IN) :: Element
       INTEGER, INTENT(OUT) :: FaceDegree(H1Basis_MaxPElementFaces), &
               FaceDirection(H1Basis_MaxPElementFaceNodes,H1Basis_MaxPElementFaces)
       INTEGER, INTENT(OUT) :: FaceMaxDegree
       INTEGER :: i

       ! Get polynomial degree of each face
       FaceMaxDegree = 0
       
       IF( Mesh % MaxFaceDofs == 0 ) THEN
         CONTINUE              

       ELSE IF (Mesh % MinFaceDOFs == Mesh % MaxFaceDOFs) THEN
          FaceMaxDegree = Mesh % Faces( Element % FaceIndexes(1) ) % PDefs % P
          FaceDegree(1:Element % Type % NumberOfFaces) = FaceMaxDegree
       ELSE
!DIR$ LOOP COUNT MAX=6
          DO i=1,Element % Type % NumberOfFaces
             IF (Mesh % Faces( Element % FaceIndexes(i) ) % BDOFs /= 0) THEN
                FaceDegree(i) = Mesh % Faces( Element % FaceIndexes(i) ) % PDefs % P
                FaceMaxDegree = MAX(FaceDegree(i), FaceMaxDegree)
             ELSE
                FaceDegree(i) = 0
             END IF
          END DO
       END IF

       ! Get face directions
       IF (FaceMaxDegree > 1) THEN
         CALL H1Basis_GetFaceDirection(Element % Type % ElementCode, &
                                       Element % Type % NumberOfFaces, &
                                       Element % NodeIndexes, &
                                       FaceDirection)
       END IF
     END SUBROUTINE GetElementMeshFaceInfo     
!------------------------------------------------------------------------------
   END FUNCTION ElementInfoVec_ComputePElementBasis
!------------------------------------------------------------------------------
   
   SUBROUTINE ElementInfoVec_ElementBasisToGlobal(npts, nbasis, nbmax, dLBasisdx, dim, cdim, LtoGMap, offset, dBasisdx)
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: npts
     INTEGER, INTENT(IN) :: nbasis
     INTEGER, INTENT(IN) :: nbmax
     REAL(KIND=dp), INTENT(IN) :: dLBasisdx(VECTOR_BLOCK_LENGTH,nbmax,3)
     INTEGER, INTENT(IN) :: dim
     INTEGER, INTENT(IN) :: cdim
     REAL(KIND=dp), INTENT(IN) :: LtoGMap(VECTOR_BLOCK_LENGTH,3,3)
     INTEGER, INTENT(IN) :: offset
     REAL(KIND=dp) CONTIG :: dBasisdx(:,:,:)

     INTEGER :: i, j, l
!DIR$ ASSUME_ALIGNED dLBasisdx:64, LtoGMap:64

     ! Map local basis function to global
     SELECT CASE (dim)
     CASE(1)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)
           END DO
         END DO
       END DO
     CASE(2)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             ! Map local basis function to global
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)+ &
                   dLBasisdx(l,i,2)*LtoGMap(l,j,2)
           END DO
         END DO
       END DO
     CASE(3)
       !DIR$ LOOP COUNT MAX=3
       DO j=1,cdim
         DO i=1,nbasis
           !_ELMER_OMP_SIMD
           DO l=1,npts
             ! Map local basis function to global
             dBasisdx(l+offset-1,i,j) = dLBasisdx(l,i,1)*LtoGMap(l,j,1)+ &
                   dLBasisdx(l,i,2)*LtoGMap(l,j,2)+ &
                   dLBasisdx(l,i,3)*LtoGMap(l,j,3)
           END DO
         END DO
       END DO
     END SELECT

   END SUBROUTINE ElementInfoVec_ElementBasisToGlobal

   
!------------------------------------------------------------------------------
!>  Returns just the size of the element at its center.
!>  providing a more economical way than calling ElementInfo. 
!------------------------------------------------------------------------------
   FUNCTION ElementSize( Element, Nodes ) RESULT ( detJ )

     TYPE(Element_t) :: Element
     TYPE(Nodes_t) :: Nodes
     REAL(KIND=dp) :: detJ

     REAL(KIND=dp) :: u,v,w
     REAL(KIND=dp), ALLOCATABLE :: Basis(:)
     INTEGER :: n,family
     LOGICAL :: Stat


     family = Element % TYPE % ElementCode / 100
     n = Element % TYPE % NumberOfNodes
     ALLOCATE( Basis(n) )

     SELECT CASE ( family )
       
       CASE ( 1 ) ! node
         DetJ = 1.0_dp
         RETURN

       CASE ( 2 ) ! line
         u = 0.0_dp
         v = 0.0_dp
         w = 0.0_dp

       CASE ( 3 ) ! tri
         u = 0.5_dp
         v = 0.5_dp
         w = 0.0_dp
         
       CASE ( 4 ) ! quad
         u = 0.0_dp
         v = 0.0_dp
         w = 0.0_dp

       CASE ( 5 ) ! tet
         u = 0.5_dp
         v = 0.5_dp
         w = 0.5_dp

       CASE ( 6 ) ! pyramid
         u = 0.0_dp
         v = 0.0_dp
         w = 0.0_dp

       CASE ( 7 ) ! wedge
         u = 0.5_dp
         v = 0.5_dp
         w = 0.0_dp

       CASE ( 8 ) ! hex
         u = 0.0_dp
         v = 0.0_dp
         w = 0.0_dp
         
       CASE DEFAULT
         CALL Fatal('ElementSize','Not implemented for elementtype')

       END SELECT

       Stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

     END FUNCTION ElementSize
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
!>  Return H(div)-conforming face element basis function values and their divergence 
!>  with respect to the reference element coordinates at a given point on the
!>  reference element. Here the basis for a real element K is constructed by  
!>  transforming the basis functions defined on the reference element k via the 
!>  Piola transformation. The data for performing the Piola transformation is also returned.
!>  Note that the reference element is chosen as in the p-approximation so that
!>  the reference element edges/faces have the same length/area. This choice simplifies 
!>  the associated assembly procedure.
!>     With giving the optional argument ApplyPiolaTransform = .TRUE., this function
!>  also performs the Piola transform, so that the basis functions and their spatial
!>  div as defined on the physical element are returned.
!>    The implementation is not yet complete as all element shapes are not supported. 
!---------------------------------------------------------------------------------
     RECURSIVE FUNCTION FaceElementInfo( Element, Nodes, u, v, w, F, detF, &
         Basis, FBasis, DivFBasis, dBasisdx, BDM, Dual, BasisDegree, &
         ApplyPiolaTransform, LeftHanded) RESULT(stat)
!------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element        !< Element structure
       TYPE(Nodes_t) :: Nodes                    !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: u                        !< 1st reference element coordinate at which the basis functions are evaluated
       REAL(KIND=dp) :: v                        !< 2nd reference element coordinate
       REAL(KIND=dp) :: w                        !< 3rd reference element coordinate
       REAL(KIND=dp), OPTIONAL :: F(3,3)         !< The gradient F=Grad f, with f the element map f:k->K
       REAL(KIND=dp) :: detF                     !< The absolute value of the determinant of the gradient matrix F
       REAL(KIND=dp) :: Basis(:)                 !< Standard nodal basis functions evaluated at (u,v,w)
       REAL(KIND=dp) :: FBasis(:,:)              !< Face element basis functions b spanning the reference element space   
       REAL(KIND=dp), OPTIONAL :: DivFBasis(:)   !< The divergence of basis functions with respect to the local coordinates
       REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)  !< The first derivatives of the H1-conforming basis functions at (u,v,w)
       LOGICAL, OPTIONAL :: BDM                  !< If .TRUE., a basis for BDM space is constructed
       LOGICAL, OPTIONAL :: Dual                 !< If .TRUE., create an alternate dual basis
       INTEGER, OPTIONAL :: BasisDegree          !< This has limited functionality at the moment
       LOGICAL, OPTIONAL :: ApplyPiolaTransform  !< If  .TRUE., perform the Piola transform so that, instead of b
                                                 !< and Div b, return  B(f(p)) and (div B)(f(p)) with B(x) the basis 
                                                 !< functions on the physical element and div the spatial divergence operator.
       LOGICAL, OPTIONAL :: LeftHanded           !< Indicates whether detF is negative
       LOGICAL :: Stat                           !< Should be .FALSE. for a degenerate element but this is not yet checked
!-----------------------------------------------------------------------------------------------------------------
!      Local variables
!------------------------------------------------------------------------------------------------------------
       INTEGER, PARAMETER :: MaxDOFs = 48 ! The largest DOF count handled, revise when new elements are added

       TYPE(Mesh_t), POINTER :: Mesh
       INTEGER, POINTER :: EdgeMap(:,:), FaceMap(:,:)
       INTEGER :: SquareFaceMap(4)
       INTEGER :: DOFs
       INTEGER :: n, dim, cdim, q, i, j, k, I1, I2
       INTEGER :: FDofMap(6,4), DofsPerFace, FaceIndices(4)
       INTEGER :: Family, RTDegree, GIndexes(27)
       REAL(KIND=dp) :: LF(3,3), LG(3,3)
       REAL(KIND=dp) :: DivBasis(MaxDOFs)
       REAL(KIND=dp) :: dLbasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3), S, D1, D2, fun, dfun, wfun(2)
       REAL(KIND=dp) :: WorkBasis(24,3), WorkDivBasis(24)

       LOGICAL :: ReverseSign(6), CreateBDMBasis, Parallel
       LOGICAL :: CreateDualBasis
       LOGICAL :: PerformPiolaTransform
!-----------------------------------------------------------------------------------------------------
       Mesh => CurrentModel % Solver % Mesh

       Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)

       !--------------------------------------------------------------------
       ! Check whether BDM or dual basis functions should be created and 
       ! whether the Piola transform is already applied within this function.
       !---------------------------------------------------------------------
       CreateBDMBasis = .FALSE.
       IF ( PRESENT(BDM) ) CreateBDMBasis = BDM
       RTDegree = 0
       IF (PRESENT(BasisDegree)) THEN
         RTDegree = BasisDegree - 1
         IF (BasisDegree > 2) CALL Fatal('ElementDescription::FaceElementInfo', 'Unsupported element degree')
       END IF
       CreateDualBasis = .FALSE.
       IF ( PRESENT(Dual) ) CreateDualBasis = Dual
       PerformPiolaTransform = .FALSE.
       IF ( PRESENT(ApplyPiolaTransform) ) PerformPiolaTransform = ApplyPiolaTransform       
       !-----------------------------------------------------------------------------------------------------
       stat = .TRUE.
       Basis = 0.0d0
       FBasis = 0.0d0
       IF (PRESENT(DivFBasis)) DivFBasis = 0.0d0
       DivBasis = 0.0d0
       LF = 0.0d0

       dLbasisdx = 0.0d0      
       n = Element % TYPE % NumberOfNodes
       dim = Element % TYPE % DIMENSION
       cdim = CoordinateSystemDimension()
       
       IF ( Element % TYPE % ElementCode == 101 ) THEN
          detF = 1.0d0
          Basis(1) = 1.0d0
          IF (PRESENT(dBasisdx)) dBasisdx(1,:) = 0.0d0
          RETURN
       END IF

       !-----------------------------------------------------------------------
       ! The standard nodal basis functions on the reference element and
       ! their derivatives with respect to the local coordinates. These define 
       ! the mapping of the reference element to an actual element on the 
       ! background mesh but are not the basis functions for face element approximation.
       ! Remark: Using reference elements having the faces of the same area
       ! simplifies the implementation of element assembly procedures.
       !-----------------------------------------------------------------------
       Family = Element % TYPE % ElementCode / 100
       SELECT CASE(Family)
       CASE(2)
         DO q=1,2
           Basis(q) = LineNodalPBasis(q, u)
           dLBasisdx(q,1) = dLineNodalPBasis(q, u)
         END DO
         IF (RTDegree == 1) THEN
           DOFs = 3
           ! Basis(3) = LineBubblePBasis(2, u)
           ! dLBasisdx(q,1) = dLineBubblePBasis(2, u)
         ELSE
           DOFs = 2
         END IF
       CASE(3)
          DO q=1,n
             Basis(q) = TriangleNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
          END DO
       CASE(4)
          DO q=1,n
             Basis(q) = QuadNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v) 
          END DO
       CASE(5)
          DO q=1,n
             Basis(q) = TetraNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
          END DO
       CASE(8)
         DO q=1,n
           Basis(q) = BrickNodalPBasis(q, u, v, w)
           dLBasisdx(q,1:3) = dBrickNodalPBasis(q, u, v, w)
         END DO
       CASE DEFAULT
          CALL Fatal('ElementDescription::FaceElementInfo','Unsupported element type')
       END SELECT          

       
       GIndexes(1:n) = Element % NodeIndexes(1:n)
       IF( Parallel ) GIndexes(1:n) = Mesh % ParallelInfo % GlobalDOFs(GIndexes(1:n))             
       
       !-----------------------------------------------------------------------
       ! Get data for performing the Piola transformation...
       !-----------------------------------------------------------------------
       stat = PiolaTransformationData(n, Element, Nodes, LF, detF, dLBasisdx) 
       !------------------------------------------------------------------------
       ! ... in order to define the basis for the element space X(K) via 
       ! applying the Piola transformation as
       !    X(K) = { B | B = 1/(det F) F b(f^{-1}(x)) }
       ! with b giving the face element basis function on the reference element k,
       ! f mapping k to the actual element K, i.e. K = f(k) and F = Grad f. This 
       ! function returns the local basis functions b and their divergence (with respect
       ! to local coordinates) evaluated at the integration point. The effect of 
       ! the Piola transformation need to be considered when integrating, so we 
       ! shall return also the values of F and det F.
       !
       ! The construction of face element bases could be done in an alternate way for 
       ! triangles and tetrahedra, while the chosen approach has the benefit that
       ! it generalizes to other cases. For example general quadrilaterals may now 
       ! be handled in the same way.
       !---------------------------------------------------------------------------
       IF (PRESENT(dBasisdx) .AND. cdim == dim) THEN
         LG = 0.0d0
         IF (cdim == dim) THEN
           SELECT CASE(Element % TYPE % ElementCode / 100)
           CASE(3,4)
             LG(1,1) = 1.0d0/detF * LF(2,2)
             LG(1,2) = -1.0d0/detF * LF(1,2)
             LG(2,1) = -1.0d0/detF * LF(2,1)
             LG(2,2) = 1.0d0/detF * LF(1,1)
           CASE(5,6,7,8)
             CALL InvertMatrix3x3(LF,LG,detF)       
           CASE DEFAULT
             CALL Fatal('ElementDescription::FaceElementInfo','Unsupported element type')
           END SELECT
           LG(1:dim,1:dim) = TRANSPOSE( LG(1:dim,1:dim) )
         END IF
       END IF
       
       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(2)
         ! TO DO: Implement possible sign reversions
         FBasis(1,1) = -Basis(1)
         DivBasis(1) = -dLBasisdx(q,1)
         FBasis(2,1) = Basis(2)
         DivBasis(2) = -dLBasisdx(q,2)
         IF (RTDegree > 0) THEN
           FBasis(3,1) = 4.0d0 * Basis(1) * Basis(2)
           DivBasis(2) = 4.0d0 * dLBasisdx(1,1) * Basis(2) + 4.0d0 * Basis(1) * dLBasisdx(2,1)
         END IF
        
       CASE(3)
          !----------------------------------------------------------------
          ! Note that the global orientation of face normal is taken to be
          ! n = t x e_z where the tangent vector t is aligned with
          ! the element edge and points towards the node that has
          ! a larger global index.
          !---------------------------------------------------------------
          EdgeMap => GetEdgeMap(3)
          !EdgeMap => GetEdgeMap(GetElementFamily(Element))

          !-----------------------------------------------------------------------------------
          ! Check first whether a sign reversion will be needed as face dofs have orientation.
          !-----------------------------------------------------------------------------------
          CALL FaceElementOrientation(Element, ReverseSign)

          IF (CreateBDMBasis) THEN
             !----------------------------------------------------------------------------
             ! This is for the BDM space of degree k=1.
             !----------------------------------------------------------------------------
             DOFs = 6
             DofsPerFace = 2
             !----------------------------------------------------------------------------
             ! First tabulate the basis functions in the default order.
             !----------------------------------------------------------------------------
             ! Two basis functions defined on face 12:
             !-------------------------------------------------
             FBasis(1,1) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + u + v)             
             FBasis(1,2) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) + 3.0d0 * u + v)
             DivBasis(1) = sqrt(3.0d0)/3.0d0
             
             FBasis(2,1) = sqrt(3.0d0)/6.0d0 * (sqrt(3.0d0) + u - v)             
             FBasis(2,2) = sqrt(3.0d0)/6.0d0 * (-sqrt(3.0d0) - 3.0d0 * u + v)
             DivBasis(2) = sqrt(3.0d0)/3.0d0

             ! Two basis functions defined on face 23:
             
             FBasis(3,1) = 1.0d0/(3.0d0+sqrt(3.0d0)) * (2.0d0+sqrt(3.0d0)+(2.0d0+sqrt(3.0d0))*u-(1.0d0+sqrt(3.0d0))*v)
             FBasis(3,2) = 1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
             DivBasis(3) = sqrt(3.0d0)/3.0d0

             FBasis(4,1) = 1.0d0/6.0d0 * (-3.0d0+sqrt(3.0d0)+(-3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
             FBasis(4,2) = 1.0d0/6.0d0 * ( 3.0d0+sqrt(3.0d0) ) * v
             DivBasis(4) = sqrt(3.0d0)/3.0d0


             ! Two basis functions defined on face 31:

             FBasis(5,1) = 1.0d0/( 3.0d0+sqrt(3.0d0) ) * ( 1.0d0 - u - v - sqrt(3.0d0)*v ) 
             FBasis(5,2) = ( 3.0d0+2.0d0*sqrt(3.0d0) ) * v /(3.0d0*(1.0d0+sqrt(3.0d0)))
             DivBasis(5) = sqrt(3.0d0)/3.0d0

             FBasis(6,1) = 1.0d0/6.0d0 * (-3.0d0-sqrt(3.0d0)+(3.0d0+sqrt(3.0d0))*u + 2.0d0*sqrt(3.0d0)*v)
             FBasis(6,2) = 1.0d0/6.0d0 * ( -3.0d0+sqrt(3.0d0) ) * v
             DivBasis(6) = sqrt(3.0d0)/3.0d0

             !-----------------------------------------------------
             ! Now do the reordering and sign reversion:
             !-----------------------------------------------------
             DO q=1,3
               IF (ReverseSign(q)) THEN
                 DO j=1,DofsPerFace
                   i = (q-1)*DofsPerFace + j
                   WorkBasis(j,1:2) = FBasis(i,1:2)
                   WorkDivBasis(j) = DivBasis(i)
                 END DO
                 i = 2*q - 1
                 FBasis(i,1:2) = -WorkBasis(2,1:2)
                 DivBasis(i) = -WorkDivBasis(2)
                 i = 2*q
                 FBasis(i,1:2) = -WorkBasis(1,1:2)
                 DivBasis(i) = -WorkDivBasis(1)
               END IF
             END DO

          ELSE
             SELECT CASE (RTDegree)
             CASE(0) 
               DOFs = 3

               FBasis(1,1) = SQRT(3.0d0)/6.0d0 * u
               FBasis(1,2) = -0.5d0 + SQRT(3.0d0)/6.0d0 * v
               DivBasis(1) =  SQRT(3.0d0)/3.0d0
               IF (ReverseSign(1)) THEN
                 FBasis(1,:) = -FBasis(1,:)
                 DivBasis(1) = -DivBasis(1)
               END IF

               FBasis(2,1) = SQRT(3.0d0)/6.0d0 * (1.0d0 + u)
               FBasis(2,2) = SQRT(3.0d0)/6.0d0 * v
               DivBasis(2) =  SQRT(3.0d0)/3.0d0        
               IF (ReverseSign(2)) THEN
                 FBasis(2,:) = -FBasis(2,:)
                 DivBasis(2) = -DivBasis(2)
               END IF

               FBasis(3,1) = SQRT(3.0d0)/6.0d0 * (-1.0d0 + u)
               FBasis(3,2) = SQRT(3.0d0)/6.0d0 * v
               DivBasis(3) =  SQRT(3.0d0)/3.0d0          
               IF (ReverseSign(3)) THEN
                 FBasis(3,:) = -FBasis(3,:)
                 DivBasis(3) = -DivBasis(3)
               END IF

             CASE(1)
               !
               ! We use a non-hierarchic basis which is motivated by flux reconstruction.
               ! The degrees of freedom associated with the element faces (edges) can be
               ! integrated for a given flux q as (q.n,w) where the weights are the Lagrange
               ! basis functions.
               DOFs = 8
               !-------------------------------------------------
               ! Two basis functions defined on the face 12.
               !-------------------------------------------------
               WorkBasis(3,1) = SQRT(3.0d0)/6.0d0 * u
               WorkBasis(3,2) = -0.5d0 + SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(3) = SQRT(3.0d0)/3.0d0
               IF (ReverseSign(1)) THEN
                 WorkBasis(3,:) = -WorkBasis(3,:)
                 WorkDivBasis(3) = -WorkDivBasis(3)
               END IF

               wfun(1) = 4.0d0 * Basis(1) - 2.0d0 * Basis(2)
               wfun(2) = 4.0d0 * Basis(2) - 2.0d0 * Basis(1)
               WorkBasis(1,1:2) = wfun(1) * WorkBasis(3,1:2)
               WorkBasis(2,1:2) = wfun(2) * WorkBasis(3,1:2)
               WorkDivBasis(1) = wfun(1) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(1,1:2) - 2.0d0 * dLBasisdx(2,1:2)))
               WorkDivBasis(2) = wfun(2) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(2,1:2) - 2.0d0 * dLBasisdx(1,1:2)))
               
               i = EdgeMap(1,1)
               j = EdgeMap(1,2)
               IF (GIndexes(j)<GIndexes(i)) THEN
                 FBasis(1,1:2) = WorkBasis(2,1:2)
                 DivBasis(1) = WorkDivBasis(2) 
                 FBasis(2,1:2) = WorkBasis(1,1:2)
                 DivBasis(2) = WorkDivBasis(1)  
               ELSE
                 FBasis(1,1:2) = WorkBasis(1,1:2)
                 DivBasis(1) = WorkDivBasis(1)
                 FBasis(2,1:2) = WorkBasis(2,1:2)
                 DivBasis(2) = WorkDivBasis(2)
               END IF

               !-------------------------------------------------
               ! Two basis functions defined on the face 23.
               !-------------------------------------------------
               WorkBasis(3,1) = SQRT(3.0d0)/6.0d0 * (1.0d0 + u)
               WorkBasis(3,2) = SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(3) =  SQRT(3.0d0)/3.0d0        
               IF (ReverseSign(2)) THEN
                 WorkBasis(3,:) = -WorkBasis(3,:)
                 WorkDivBasis(3) = -WorkDivBasis(3)
               END IF

               wfun(1) = 4.0d0 * Basis(2) - 2.0d0 * Basis(3)
               wfun(2) = 4.0d0 * Basis(3) - 2.0d0 * Basis(2)
               WorkBasis(1,1:2) = wfun(1) * WorkBasis(3,1:2)
               WorkBasis(2,1:2) = wfun(2) * WorkBasis(3,1:2)
               WorkDivBasis(1) = wfun(1) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(2,1:2) - 2.0d0 * dLBasisdx(3,1:2)))
               WorkDivBasis(2) = wfun(2) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(3,1:2) - 2.0d0 * dLBasisdx(2,1:2)))

               i = EdgeMap(2,1)
               j = EdgeMap(2,2)
               IF (GIndexes(j)<GIndexes(i)) THEN
                 FBasis(3,1:2) = WorkBasis(2,1:2)
                 DivBasis(3) = WorkDivBasis(2) 
                 FBasis(4,1:2) = WorkBasis(1,1:2)
                 DivBasis(4) = WorkDivBasis(1)  
               ELSE
                 FBasis(3,1:2) = WorkBasis(1,1:2)
                 DivBasis(3) = WorkDivBasis(1)
                 FBasis(4,1:2) = WorkBasis(2,1:2)
                 DivBasis(4) = WorkDivBasis(2)
               END IF
               
               !-------------------------------------------------
               ! Two basis functions defined on the face 31.
               !-------------------------------------------------
               WorkBasis(3,1) = SQRT(3.0d0)/6.0d0 * (-1.0d0 + u)
               WorkBasis(3,2) = SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(3) =  SQRT(3.0d0)/3.0d0          
               IF (ReverseSign(3)) THEN
                 WorkBasis(3,:) = -WorkBasis(3,:)
                 WorkDivBasis(3) = -WorkDivBasis(3)
               END IF
               
               wfun(1) = 4.0d0 * Basis(3) - 2.0d0 * Basis(1)
               wfun(2) = 4.0d0 * Basis(1) - 2.0d0 * Basis(3)
               WorkBasis(1,1:2) = wfun(1) * WorkBasis(3,1:2)
               WorkBasis(2,1:2) = wfun(2) * WorkBasis(3,1:2)
               WorkDivBasis(1) = wfun(1) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(3,1:2) - 2.0d0 * dLBasisdx(1,1:2)))
               WorkDivBasis(2) = wfun(2) * WorkDivBasis(3) + &
                   SUM(WorkBasis(3,1:2) * (4.0d0 * dLBasisdx(1,1:2) - 2.0d0 * dLBasisdx(3,1:2)))
               
               i = EdgeMap(3,1)
               j = EdgeMap(3,2)
               IF (GIndexes(j)<GIndexes(i)) THEN
                 FBasis(5,1:2) = WorkBasis(2,1:2)
                 DivBasis(5) = WorkDivBasis(2) 
                 FBasis(6,1:2) = WorkBasis(1,1:2)
                 DivBasis(6) = WorkDivBasis(1)  
               ELSE
                 FBasis(5,1:2) = WorkBasis(1,1:2)
                 DivBasis(5) = WorkDivBasis(1)
                 FBasis(6,1:2) = WorkBasis(2,1:2)
                 DivBasis(6) = WorkDivBasis(2)
               END IF

               !-------------------------------------------------
               ! Two basis functions defined on the interior 123.
               ! Note: The ordering of these functions is not specified,
               !       although the choice is made unique.
               !-------------------------------------------------               
               WorkBasis(1,1) = SQRT(3.0d0)/6.0d0 * u
               WorkBasis(1,2) = -0.5d0 + SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(1) = Basis(3) * SQRT(3.0d0)/3.0d0 + SUM(WorkBasis(1,1:2) * dLBasisdx(3,1:2))
               WorkBasis(1,1:2) = Basis(3) * WorkBasis(1,1:2)
               
               WorkBasis(2,1) = SQRT(3.0d0)/6.0d0 * (1.0d0 + u)
               WorkBasis(2,2) = SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(2) = Basis(1) * SQRT(3.0d0)/3.0d0 + SUM(WorkBasis(2,1:2) * dLBasisdx(1,1:2))
               WorkBasis(2,1:2) = Basis(1) * WorkBasis(2,1:2)
               
               WorkBasis(3,1) = SQRT(3.0d0)/6.0d0 * (-1.0d0 + u)
               WorkBasis(3,2) = SQRT(3.0d0)/6.0d0 * v
               WorkDivBasis(3) = Basis(2) * SQRT(3.0d0)/3.0d0 + SUM(WorkBasis(3,1:2) * dLBasisdx(2,1:2))
               WorkBasis(3,1:2) = Basis(2) * WorkBasis(3,1:2)
               
               FaceIndices(1:3) = GIndexes(1:3)
               IF ( FaceIndices(1) < FaceIndices(2) ) THEN
                 k = 1
               ELSE
                 k = 2
               END IF
               IF ( FaceIndices(k) > FaceIndices(3) ) THEN
                 k = 3
               END IF

               SELECT CASE(k)
               CASE(1)
                 FBasis(7,1:2) = WorkBasis(1,1:2)
                 DivBasis(7) = WorkDivBasis(1)
                 FBasis(8,1:2) = WorkBasis(3,1:2)
                 DivBasis(8) = WorkDivBasis(3)
               CASE(2)
                 FBasis(7,1:2) = WorkBasis(1,1:2)
                 DivBasis(7) = WorkDivBasis(1)
                 FBasis(8,1:2) = WorkBasis(2,1:2)
                 DivBasis(8) = WorkDivBasis(2)
               CASE(3)
                 FBasis(7,1:2) = WorkBasis(2,1:2)
                 DivBasis(7) = WorkDivBasis(2)
                 FBasis(8,1:2) = WorkBasis(3,1:2)
                 DivBasis(8) = WorkDivBasis(3)
               END SELECT
               
             END SELECT
          END IF
          
       CASE(4)
          DOFs = 6
          !--------------------------------------------------------------------
          ! Quadrilateral Arnold-Boffi-Falk (ABF) element basis of degree k=0
          !--------------------------------------------------------------------
          EdgeMap => GetEdgeMap(4)
          SquareFaceMap(:) = (/ 1,2,3,4 /)          

          IF (.NOT. CreateDualBasis) THEN
             !-------------------------------------------------
             ! Four basis functions defined on the edges
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             FBasis(1,1) = 0.0d0
             FBasis(1,2) = -((-1.0d0 + v)*v)/4.0d0
             DivBasis(1) = (1.0d0 - 2*v)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(1,:) = -FBasis(1,:)
               DivBasis(1) = -DivBasis(1)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             FBasis(2,1) = (u*(1.0d0 + u))/4.0d0
             FBasis(2,2) = 0.0d0
             DivBasis(2) = (1 + 2.0d0*u)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(2,:) = -FBasis(2,:)
               DivBasis(2) = -DivBasis(2)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             FBasis(3,1) = 0.0d0
             FBasis(3,2) = (v*(1.0d0 + v))/4.0d0
             DivBasis(3) = (1.0d0 + 2.0d0*v)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(3,:) = -FBasis(3,:)
               DivBasis(3) = -DivBasis(3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             FBasis(4,1) = -((-1.0d0 + u)*u)/4.0d0
             FBasis(4,2) = 0.0d0
             DivBasis(4) = (1.0d0 - 2.0d0*u)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(4,:) = -FBasis(4,:)
               DivBasis(4) = -DivBasis(4)
             END IF

             !--------------------------------------------------------------------
             ! Additional two basis functions associated with the element interior
             !-------------------------------------------------------------------
             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkDivBasis(:) = 0.0d0

             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = (-1.0d0 + v**2)/2.0d0
             WorkDivBasis(1) = v

             WorkBasis(2,1) = (1.0d0 - u**2)/2.0d0
             WorkBasis(2,2) = 0.0d0
             WorkDivBasis(2) = -u

             FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             FBasis(5,:) = D1 * WorkBasis(I1,:)
             DivBasis(5) = D1 * WorkDivBasis(I1)
             FBasis(6,:) = D2 * WorkBasis(I2,:)
             DivBasis(6) = D2 * WorkDivBasis(I2)   
          ELSE
             !---------------------------------------------------------------------------
             ! Create alternate basis functions for the ABF space so that these basis
             ! functions are dual to the standard basis functions when the mesh is regular.
             ! First four basis functions which are dual to the standard edge basis 
             ! functions:
             !----------------------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             FBasis(1,1) = 0.0d0
             FBasis(1,2) = (-3.0d0*(-1.0d0 - 2.0d0*v + 5.0d0*v**2))/4.0d0
             DivBasis(1) = (-3.0d0*(-1.0d0 + 5.0d0*v))/2.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(1,:) = -FBasis(1,:)
               DivBasis(1) = -DivBasis(1)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             FBasis(2,1) = (3.0d0*(-1.0d0 + 2.0d0*u + 5.0d0*u**2))/4.0d0
             FBasis(2,2) = 0.0d0
             DivBasis(2) = (3.0d0*(1.0d0 + 5.0d0*u))/2.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivBasis(2) = -DivBasis(2)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             FBasis(3,1) = 0.0d0
             FBasis(3,2) = (3.0d0*(-1.0d0 + 2.0d0*v + 5.0d0*v**2))/4.0d0
             DivBasis(3) = (3.0d0*(1.0d0 + 5.0d0*v))/2.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(3,:) = -FBasis(3,:)
               DivBasis(3) = -DivBasis(3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             FBasis(4,1) = (-3.0d0*(-1.0d0 - 2.0d0*u + 5.0d0*u**2))/4.0d0
             FBasis(4,2) = 0.0d0
             DivBasis(4) = (-3.0d0*(-1.0d0 + 5.0d0*u))/2.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               FBasis(4,:) = -FBasis(4,:)
               DivBasis(4) = -DivBasis(4)
             END IF

             !-------------------------------------------------------------------------
             ! Additional two dual basis functions associated with the element interior
             !-------------------------------------------------------------------------
             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkDivBasis(:) = 0.0d0

             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = (3.0d0*(-3.0d0 + 5.0d0*v**2))/8.0d0
             WorkDivBasis(1) = 15.0d0*v/4.0d0

             WorkBasis(2,1) = (3.0d0*(3.0d0 - 5.0d0*u**2))/8.0d0
             WorkBasis(2,2) = 0.0d0
             WorkDivBasis(2) = -15.0d0*u/4.0d0

             FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))              
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             FBasis(5,:) = D1 * WorkBasis(I1,:)
             DivBasis(5) = D1 * WorkDivBasis(I1)
             FBasis(6,:) = D2 * WorkBasis(I2,:)
             DivBasis(6) = D2 * WorkDivBasis(I2)
          END IF

       CASE(5)
          !-----------------------------------------
          ! This branch is for handling tetrahedra
          !-----------------------------------------------------------------------------------
          ! Check first whether a sign reversion will be needed as face dofs have orientation.
          ! If the sign is not reversed, the positive value of the degree of freedom produces
          ! positive outward flux from the element through the face handled.
          !-----------------------------------------------------------------------------------
          CALL FaceElementOrientation(Element, ReverseSign)

          IF (CreateBDMBasis) THEN
             DOFs = 12
             DofsPerFace = 3 ! This choice is used for the BDM space of degree k=1
             !----------------------------------------------------------------------------
             ! Create a table of BDM basis functions in the default order
             !----------------------------------------------------------------------------
             ! Face {213}:
             WorkBasis(1,1) = (3*Sqrt(6.0d0) + 2*Sqrt(6.0d0)*u - 3*Sqrt(2.0d0)*v - 3*w)/12.0_dp
             WorkBasis(1,2) = (-2*Sqrt(2.0d0) - 3*Sqrt(2.0d0)*u + Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(1,3) = (-8 - 12*u + 4*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0_dp

             WorkBasis(2,1) = (2*Sqrt(6.0d0)*u + 3*(-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w))/12.0_dp
             WorkBasis(2,2) = (-2*Sqrt(2.0d0) + 3*Sqrt(2.0d0)*u + Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(2,3) = u + (-8 + 4*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0_dp

             WorkBasis(3,1) = -u/(2.0*Sqrt(6.0d0))
             WorkBasis(3,2) = (Sqrt(2.0d0) + 3*Sqrt(6.0d0)*v - 2*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(3,3) = (4 - 8*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w)/12.0_dp

             ! Face {124}:
             WorkBasis(4,1) = (2*Sqrt(6.0d0)*u + 3*(-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w))/12.0_dp
             WorkBasis(4,2) = (-6*Sqrt(2.0d0) + 9*Sqrt(2.0d0)*u + 2*Sqrt(6.0d0)*v + 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(4,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(5,1) = (3*Sqrt(6.0d0) + 2*Sqrt(6.0d0)*u - 3*Sqrt(2.0d0)*v - 3*w)/12.0_dp
             WorkBasis(5,2) = (-6*Sqrt(2.0d0) - 9*Sqrt(2.0d0)*u + 2*Sqrt(6.0d0)*v + 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(5,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(6,1) = -u/(2.0*Sqrt(6.0d0))
             WorkBasis(6,2) = (3*Sqrt(2.0d0) - Sqrt(6.0d0)*v - 6*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(6,3) = (5*w)/(2.0*Sqrt(6.0d0))

             ! Face {234}:
             WorkBasis(7,1) = (5*Sqrt(6.0d0) + 5*Sqrt(6.0d0)*u - 6*Sqrt(2.0d0)*v - 6*w)/12.0_dp
             WorkBasis(7,2) = -v/(2.0*Sqrt(6.0d0))
             WorkBasis(7,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(8,1) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u + 6*Sqrt(2.0d0)*v - 3*w)/12.0_dp
             WorkBasis(8,2) = (5*Sqrt(6.0)*v - 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(8,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(9,1) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u + 9*w)/12.0_dp
             WorkBasis(9,2) = (-(Sqrt(6.0d0)*v) + 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(9,3) = (5*w)/(2.0*Sqrt(6.0d0))

             ! Face {314}:             
             WorkBasis(10,1) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 6*Sqrt(2.0d0)*v + 3*w)/12.0_dp
             WorkBasis(10,2) = (5*Sqrt(6.0d0)*v - 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(10,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(11,1) = (-5*Sqrt(6.0d0) + 5*Sqrt(6.0d0)*u + 6*Sqrt(2.0d0)*v + 6*w)/12.0_dp
             WorkBasis(11,2) = -v/(2.0*Sqrt(6.0d0))
             WorkBasis(11,3) = -w/(2.0*Sqrt(6.0d0))
             WorkBasis(12,1) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 9*w)/12.0_dp
             WorkBasis(12,2) = (-(Sqrt(6.0d0)*v) + 3*Sqrt(3.0d0)*w)/12.0_dp
             WorkBasis(12,3) = (5*w)/(2.0*Sqrt(6.0d0))

             !----------------------------------------------------------------------
             ! Find out how face basis functions must be ordered so that the global
             ! indexing convention is respected. 
             !-----------------------------------------------------------------------
             CALL FaceElementBasisOrdering(Element, FDofMap(1:4,1:3))

             !-----------------------------------------------------
             ! Now do the actual reordering and sign reversion
             !-----------------------------------------------------
             DO q=1,4
                IF (ReverseSign(q)) THEN
                   S = -1.0d0
                ELSE
                   S = 1.0d0
                END IF

                DO j=1,DofsPerFace
                   k = FDofMap(q,j)
                   i = (q-1)*DofsPerFace + j
                   FBasis(i,:) = S * WorkBasis((q-1)*DofsPerFace+k,:)
                   DivBasis(i) = S * sqrt(3.0d0)/(2.0d0*sqrt(2.0d0))
                END DO
             END DO

          ELSE
             DOFs = 4
             !-------------------------------------------------------------------------
             ! The basis functions that define RT space on reference element
             !-----------------------------------------------------------------------
             FBasis(1,1) = SQRT(2.0d0)/4.0d0 * u
             FBasis(1,2) = -SQRT(6.0d0)/12.0d0 + SQRT(2.0d0)/4.0d0 * v
             FBasis(1,3) = -1.0d0/SQRT(3.0d0) + SQRT(2.0d0)/4.0d0 * w
             DivBasis(1) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( ReverseSign(1) ) THEN
                FBasis(1,:) = -FBasis(1,:)
                DivBasis(1) = -DivBasis(1)
             END IF

             FBasis(2,1) = SQRT(2.0d0)/4.0d0 * u
             FBasis(2,2) = -SQRT(6.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * v
             FBasis(2,3) = SQRT(2.0d0)/4.0d0 * w
             DivBasis(2) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( ReverseSign(2) ) THEN
                FBasis(2,:) = -FBasis(2,:)
                DivBasis(2) = -DivBasis(2)
             END IF

             FBasis(3,1) = SQRT(2.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * u
             FBasis(3,2) = SQRT(2.0d0)/4.0d0 * v
             FBasis(3,3) = SQRT(2.0d0)/4.0d0 * w
             DivBasis(3) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( ReverseSign(3) ) THEN
                FBasis(3,:) = -FBasis(3,:)
                DivBasis(3) = -DivBasis(3)
             END IF

             FBasis(4,1) = -SQRT(2.0d0)/4.0d0 + SQRT(2.0d0)/4.0d0 * u
             FBasis(4,2) = SQRT(2.0d0)/4.0d0 * v
             FBasis(4,3) = SQRT(2.0d0)/4.0d0 * w
             DivBasis(4) = 3.0d0*SQRT(2.0d0)/4.0d0
             IF ( ReverseSign(4) ) THEN
                FBasis(4,:) = -FBasis(4,:)
                DivBasis(4) = -DivBasis(4)
             END IF
          END IF
       CASE(8)
         !--------------------------------------------------------------
         ! This branch is for handling brick elements
         !--------------------------------------------------------------  
         ! Check first whether a sign reverse will be needed.
         ! If the sign is not reversed, the positive value of the degree of freedom produces
         ! positive outward flux from the element through the face handled.
         !-----------------------------------------------------------------------------------
         CALL FaceElementOrientation(Element, ReverseSign)

         DOFs = 48   ! 4 DOFs per face and 24 elementwise DOFs
         DofsPerFace = 4
         WorkBasis = 0.0d0

         !
         ! Face 2143:
         !
         SquareFaceMap(:) = (/ 2,1,4,3 /)
         DO q=1,4
           WorkBasis(q,3) = -1.0d0 * QuadNodalPBasis(SquareFaceMap(q), u, v) * LineNodalPBasis(1, w)
           WorkDivBasis(q) = -1.0d0 * QuadNodalPBasis(SquareFaceMap(q), u, v) * dLineNodalPBasis(1, w)
         END DO

         !
         ! Face 5678:
         !
         DO q=1,4
           WorkBasis(4+q,3) = QuadNodalPBasis(q, u, v) * LineNodalPBasis(2, w)
           WorkDivBasis(4+q) = QuadNodalPBasis(q, u, v) * dLineNodalPBasis(2, w)
         END DO
         
         !
         ! Face 1265:
         !         
         DO q=1,4
           WorkBasis(8+q,2) = -1.0d0 * QuadNodalPBasis(q, u, w) * LineNodalPBasis(1, v)
           WorkDivBasis(8+q) = -1.0d0 * QuadNodalPBasis(q, u, w) * dLineNodalPBasis(1, v)
         END DO

         !
         ! Face 2376:
         !         
         DO q=1,4
           WorkBasis(12+q,1) = QuadNodalPBasis(q, v, w) * LineNodalPBasis(2, u)
           WorkDivBasis(12+q) = QuadNodalPBasis(q, v, w) * dLineNodalPBasis(2, u)
         END DO

         !
         ! Face 3487:
         !
         SquareFaceMap(:) = (/ 2,1,4,3 /)
         DO q=1,4
           WorkBasis(16+q,2) = QuadNodalPBasis(SquareFaceMap(q), u, w) * LineNodalPBasis(2, v)
           WorkDivBasis(16+q) = QuadNodalPBasis(SquareFaceMap(q), u, w) * dLineNodalPBasis(2, v)
         END DO

         !
         ! Face 4152:
         !
         SquareFaceMap(:) = (/ 2,1,4,3 /)
         DO q=1,4
           WorkBasis(20+q,1) = -1.0d0 * QuadNodalPBasis(SquareFaceMap(q), v, w) * LineNodalPBasis(1, u)
           WorkDivBasis(20+q) = -1.0d0 * QuadNodalPBasis(SquareFaceMap(q), v, w) * dLineNodalPBasis(1, u)
         END DO

         !----------------------------------------------------------------------
         ! Find out how face basis functions must be ordered so that the global
         ! indexing convention is respected. 
         !-----------------------------------------------------------------------
         CALL FaceElementBasisOrdering(Element, FDofMap(1:6,1:4))

         !-----------------------------------------------------
         ! Now do the actual reordering and sign reverses
         !-----------------------------------------------------
         DO q=1,6
           IF (ReverseSign(q)) THEN
             S = -1.0d0
           ELSE
             S = 1.0d0
           END IF

           DO j=1,DofsPerFace
             k = FDofMap(q,j)
             i = (q-1)*DofsPerFace + j
             FBasis(i,:) = S * WorkBasis((q-1)*DofsPerFace+k,:)
             DivBasis(i) = S * WorkDivBasis((q-1)*DofsPerFace+k)
           END DO
         END DO

         !
         ! 24 interior basis functions (8 per coordinate direction)
         !
         k = 24
         DO j=1,2
           SELECT CASE(j)
           CASE(1)
             fun = 1.0d0
             dfun = 0.0d0
           CASE(2)
             fun = 2.0d0 * u
             dfun = 2.0d0
           END SELECT
           DO q=1,4
             k = k + 1
             FBasis(k,1) = QuadNodalPBasis(q, v, w) * LineNodalPBasis(1, u) * LineNodalPBasis(2, u) * fun
             DivBasis(k) = QuadNodalPBasis(q, v, w) * ( dLineNodalPBasis(1, u) * LineNodalPBasis(2, u) * fun + &
                 LineNodalPBasis(1, u) * dLineNodalPBasis(2, u) * fun + &
                 LineNodalPBasis(1, u) * LineNodalPBasis(2, u) * dfun )
           END DO
         END DO

         DO j=1,2
           SELECT CASE(j)
           CASE(1)
             fun = 1.0d0
             dfun = 0.0d0
           CASE(2)
             fun = 2.0d0 * v
             dfun = 2.0d0
           END SELECT
           DO q=1,4
             k = k + 1
             FBasis(k,2) = QuadNodalPBasis(q, u, w) * LineNodalPBasis(1, v) * LineNodalPBasis(2, v) * fun
             DivBasis(k) = QuadNodalPBasis(q, u, w) * ( dLineNodalPBasis(1, v) * LineNodalPBasis(2, v) * fun + &
                 LineNodalPBasis(1, v) * dLineNodalPBasis(2, v) * fun + &
                 LineNodalPBasis(1, v) * LineNodalPBasis(2, v) * dfun )
           END DO
         END DO

         DO j=1,2
           SELECT CASE(j)
           CASE(1)
             fun = 1.0d0
             dfun = 0.0d0
           CASE(2)
             fun = 2.0d0 * w
             dfun = 2.0d0
           END SELECT
           DO q=1,4
             k = k + 1
             FBasis(k,3) = QuadNodalPBasis(q, u, v) * LineNodalPBasis(1, w) * LineNodalPBasis(2, w) * fun
             DivBasis(k) = QuadNodalPBasis(q, u, v) * ( dLineNodalPBasis(1, w) * LineNodalPBasis(2, w) * fun + &
                 LineNodalPBasis(1, w) * dLineNodalPBasis(2, w) * fun + &
                 LineNodalPBasis(1, w) * LineNodalPBasis(2, w) * dfun )
           END DO
         END DO

       CASE DEFAULT
          CALL Fatal('ElementDescription::FaceElementInfo','Unsupported element type')
       END SELECT

       IF (PerformPiolaTransform) THEN
         DO j=1,DOFs
           DO k=1,dim
             WorkBasis(1,k) = SUM( LF(k,1:dim) * FBasis(j,1:dim) )
           END DO
           FBasis(j,1:dim) = 1.0d0/DetF * WorkBasis(1,1:dim)
           
           DivBasis(j) = 1.0d0/DetF * DivBasis(j)
         END DO
         ! Make the returned value DetF to act as a metric term for integration
         ! over the volume of the element:
         IF (PRESENT(LeftHanded)) LeftHanded = detF < 0.0d0
         DetF = ABS(DetF)
       END IF

       ! ----------------------------------------------------------------------
       ! Get global first derivatives of the nodal basis functions if wanted:
       ! ----------------------------------------------------------------------
       IF ( PRESENT(dBasisdx) ) THEN
         dBasisdx = 0.0d0
         IF (cdim == dim) THEN       
           DO i=1,n
             DO j=1,dim
               DO k=1,dim
                 dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LG(j,k)
               END DO
             END DO
           END DO
         ELSE
           CALL Warn('ElementDescription::FaceElementInfo', &
               'Cannot return gradient for elements embedded in a higher-dimensional space')
         END IF
       END IF
       
       IF (PRESENT(F)) F = LF
       IF (PRESENT(DivFBasis)) DivFBasis(1:DOFs) = DivBasis(1:DOFs)
!-----------------------------------------------------------------------------
     END FUNCTION FaceElementInfo
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------
!> This function returns data for performing the Piola transformation 
!------------------------------------------------------------------------------------------------
     FUNCTION PiolaTransformationData(nn,Element,Nodes,F,DetF,dLBasisdx) RESULT(Success)
!-------------------------------------------------------------------------------------------------
       INTEGER :: nn                   !< The number of classic nodes used in the element mapping
       TYPE(Element_t) :: Element      !< Element structure
       TYPE(Nodes_t) :: Nodes          !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: F(:,:)         !< The gradient of the element mapping
       REAL(KIND=dp) :: DetF           !< The determinant of the gradient matrix (or the Jacobian matrix)
       REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of nodal basis functions with respect to local coordinates
       LOGICAL :: Success              !< Could and should return .FALSE. if the element is degenerate
!-----------------------------------------------------------------------------------------------------
!      Local variables
!-------------------------------------------------------------------------------------------------
       REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
       INTEGER :: cdim,dim,n,i
!-------------------------------------------------------------------------------------------------
       x => Nodes % x
       y => Nodes % y
       z => Nodes % z     

       ! cdim = CoordinateSystemDimension()
       n = MIN( SIZE(x), nn )
       dim  = Element % TYPE % DIMENSION

       !------------------------------------------------------------------------------
       ! The gradient of the element mapping K = f(k), with k the reference element
       !------------------------------------------------------------------------------
       F = 0.0d0
       DO i=1,dim
          F(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
          F(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
          !IF (dim == 3) &
          ! In addition to the case dim = 3, the following entries may be useful  
          ! with dim=2 when natural BCs in 3-D are handled. 
          F(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
       END DO

       SELECT CASE( dim )
       CASE(1)
          DetF = sqrt(SUM(F(1:3,1)**2))
       CASE (2)
          DetF = F(1,1)*F(2,2) - F(1,2)*F(2,1)
       CASE(3)
          DetF = F(1,1) * ( F(2,2)*F(3,3) - F(2,3)*F(3,2) ) + &
               F(1,2) * ( F(2,3)*F(3,1) - F(2,1)*F(3,3) ) + &
               F(1,3) * ( F(2,1)*F(3,2) - F(2,2)*F(3,1) )
       END SELECT

       success = .TRUE.
!------------------------------------------------
     END FUNCTION PiolaTransformationData
!------------------------------------------------

!-----------------------------------------------------------------------------------
!> Get information about whether a sign reversion will be needed to obtain right
!> DOFs for face (vector) elements. If the sign is not reversed, the positive value of 
!> the degree of freedom produces positive outward flux from the element through 
!> the face handled.
!-----------------------------------------------------------------------------------
SUBROUTINE FaceElementOrientation(Element, ReverseSign, FaceIndex, Nodes)
!-----------------------------------------------------------------------------------
  IMPLICIT NONE

  TYPE(Element_t), INTENT(IN) :: Element       !< A 3-D/2-D element having 2-D/1-D faces 
  LOGICAL, INTENT(OUT) :: ReverseSign(:)       !< Face-wise information about the sign reversions
  INTEGER, OPTIONAL, INTENT(IN) :: FaceIndex   !< Check just one face that is specified here
  TYPE(Nodes_t), OPTIONAL :: Nodes             !< An inactive variable related to code verification
!-----------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: Parallel
  
  INTEGER, POINTER :: FaceMap(:,:)
  INTEGER, TARGET :: TetraFaceMap(4,3), BrickFaceMap(6,4)
  INTEGER :: FaceIndices(4), GIndexes(27)
  INTEGER :: j, q, first_face, last_face

  ! Some inactive variables that were used in the code verification
  LOGICAL :: ReverseSign2(4), CheckSignReversions
  INTEGER :: n, i, k, A, B, C, D, I1, I2
  REAL(KIND=dp) :: t1(3), t2(3), m(3), e(3), D1, D2
!-----------------------------------------------------------------------------------
  ReverseSign(:) = .FALSE.

  IF (PRESENT(FaceIndex)) THEN
    first_face = FaceIndex
    last_face = FaceIndex
  ELSE
    first_face = 1
  END IF

  Mesh => CurrentModel % Solver % Mesh
  Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)

  n = Element % Type % NumberOfNodes
  GIndexes(1:n) = Element % NodeIndexes(1:n)
  IF( Parallel ) GIndexes(1:n) = Mesh % ParallelInfo % GlobalDOFs(GIndexes(1:n))
  
  SELECT CASE(Element % TYPE % ElementCode / 100)
  CASE(3)
    FaceMap => GetEdgeMap(3) 

    IF (.NOT. PRESENT(FaceIndex)) last_face = 3
    IF (SIZE(ReverseSign) < last_face) CALL Fatal('FaceElementOrientation', &
        'Too small array for listing element faces')
    
    DO q=first_face,last_face
      FaceIndices(1:2) = GIndexes((FaceMap(q,1:2)))
      IF (FaceIndices(2) < FaceIndices(1)) ReverseSign(q) = .TRUE.
    END DO

  CASE(4)
    FaceMap => GetEdgeMap(4)

    IF (.NOT. PRESENT(FaceIndex)) last_face = 4
    IF (SIZE(ReverseSign) < last_face) CALL Fatal('FaceElementOrientation', &
        'Too small array for listing element faces')
    
    DO q=first_face,last_face
      FaceIndices(1:2) = GIndexes((FaceMap(q,1:2)))
      IF (FaceIndices(2) < FaceIndices(1)) ReverseSign(q) = .TRUE.
    END DO

  CASE(5)
    TetraFaceMap(1,:) = (/ 2, 1, 3 /)
    TetraFaceMap(2,:) = (/ 1, 2, 4 /)
    TetraFaceMap(3,:) = (/ 2, 3, 4 /) 
    TetraFaceMap(4,:) = (/ 3, 1, 4 /)

    FaceMap => TetraFaceMap

    IF (.NOT. PRESENT(FaceIndex)) last_face = 4
    IF (SIZE(ReverseSign) < last_face) CALL Fatal('FaceElementOrientation', &
        'Too small array for listing element faces')

    DO q=first_face,last_face
      FaceIndices(1:3) = GIndexes(FaceMap(q,1:3))
      IF ( (FaceIndices(1) < FaceIndices(2)) .AND. (FaceIndices(1) < FaceIndices(3)) ) THEN
        IF (FaceIndices(3) < FaceIndices(2)) THEN
          ReverseSign(q) = .TRUE.
        END IF
      ELSE IF ( ( FaceIndices(2) < FaceIndices(1) ) .AND. ( FaceIndices(2) < FaceIndices(3) ) ) THEN
        IF ( FaceIndices(1) < FaceIndices(3) ) THEN
          ReverseSign(q) = .TRUE.
        END IF
      ELSE  
        IF ( FaceIndices(2) < FaceIndices(1) ) THEN
          ReverseSign(q) = .TRUE.
        END IF
      END IF
    END DO

    !----------------------------------------------------------------------
    ! Another way for finding sign reversions in the case of tetrahedron. 
    ! This code is retained here, although it was used for verification purposes...
    !----------------------------------------------------------------------
    CheckSignReversions = .FALSE.
    IF (CheckSignReversions) THEN
      DO q=1,4
        ReverseSign2(q) = .FALSE.
        i = FaceMap(q,1)
        j = FaceMap(q,2)
        k = FaceMap(q,3)

        IF ( ( GIndexes(i) < GIndexes(j) ) .AND. ( GIndexes(i) < GIndexes(k) ) ) THEN
          A = i
          IF (GIndexes(j) < GIndexes(k)) THEN
            B = j
            C = k
          ELSE
            B = k
            C = j
          END IF
        ELSE IF ( ( GIndexes(j) < GIndexes(i) ) .AND. ( GIndexes(j) < GIndexes(k) ) ) THEN
          A = j
          IF (GIndexes(i) < GIndexes(k)) THEN
            B = i
            C = k
          ELSE
            B = k
            C = i
          END IF
        ELSE
          A = k
          IF (GIndexes(i) < GIndexes(j)) THEN
            B = i
            C = j
          ELSE
            B = j
            C = i
          END IF
        END IF

        t1(1) = Nodes % x(B) - Nodes % x(A)
        t1(2) = Nodes % y(B) - Nodes % y(A)              
        t1(3) = Nodes % z(B) - Nodes % z(A)

        t2(1) = Nodes % x(C) - Nodes % x(A)
        t2(2) = Nodes % y(C) - Nodes % y(A)              
        t2(3) = Nodes % z(C) - Nodes % z(A)

        m(1:3) = CrossProduct(t1,t2)

        SELECT CASE(q)
        CASE(1)
          D = 4
        CASE(2)
          D = 3 
        CASE(3)
          D = 1
        CASE(4)
          D = 2                   
        END SELECT

        e(1) = Nodes % x(D) - Nodes % x(A)
        e(2) = Nodes % y(D) - Nodes % y(A)                
        e(3) = Nodes % z(D) - Nodes % z(A)  

        IF ( SUM(m(1:3) * e(1:3)) > 0.0d0 ) ReverseSign2(q) = .TRUE.

      END DO

      IF ( ANY(ReverseSign(1:4) .NEQV. ReverseSign2(1:4)) ) THEN
        PRINT *, 'CONFLICTING SIGN REVERSIONS SUGGESTED'
        PRINT *, ReverseSign(1:4)
        PRINT *, ReverseSign2(1:4)
        STOP EXIT_ERROR
      END IF
    END IF

  CASE(8)
    !
    ! Write the face map such that by default the normal points outwards
    ! from the brick:
    !
    BrickFaceMap(1,:) = (/ 2, 1, 4, 3 /)
    BrickFaceMap(2,:) = (/ 5, 6, 7, 8 /)
    BrickFaceMap(3,:) = (/ 1, 2, 6, 5 /)
    BrickFaceMap(4,:) = (/ 2, 3, 7, 6 /)
    BrickFaceMap(5,:) = (/ 3, 4, 8, 7 /)
    BrickFaceMap(6,:) = (/ 4, 1, 5, 8 /)

    FaceMap => BrickFaceMap

    IF (.NOT. PRESENT(FaceIndex)) last_face = 6
    IF (SIZE(ReverseSign) < last_face) CALL Fatal('FaceElementOrientation', &
        'Too small array for listing element faces')

    DO q=first_face,last_face
      FaceIndices(1:4) = GIndexes(FaceMap(q,1:4))    
      CALL SquareFaceDofsOrdering(I1, I2, D1, D2, FaceIndices(1:4), ReverseSign(q))
    END DO

  CASE DEFAULT
    CALL Fatal('FaceElementOrientation', 'Unsupported element family')
  END SELECT
!-----------------------------------------------------------------------------------
END SUBROUTINE FaceElementOrientation
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
!> This subroutine produces information about how the basis functions of face (vector)
!> elements have to be reordered to conform with the global ordering convention.
!-----------------------------------------------------------------------------------
SUBROUTINE FaceElementBasisOrdering(Element, FDofMap, FaceIndex, ReverseSign)
!-----------------------------------------------------------------------------------
  IMPLICIT NONE

  TYPE(Element_t), INTENT(IN) :: Element       !< A 3-D element having 2-D faces
  INTEGER, INTENT(OUT) :: FDofMap(:,:)         !< Face-wise information for the basis permutation  
  INTEGER, OPTIONAL, INTENT(IN) :: FaceIndex   !< Check just one face that is specified here
  LOGICAL, OPTIONAL, INTENT(OUT) :: ReverseSign(:) !< For bricks face-wise information about the sign reversions
!-----------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh 
  LOGICAL :: Parallel
  LOGICAL :: ReverseNormal(6)
  INTEGER, POINTER :: FaceMap(:,:)
  INTEGER, TARGET :: TetraFaceMap(4,3), BrickFaceMap(6,4), FaceIndices(4), GIndexes(27)
  INTEGER :: n, i, j, k, l, q, first_face, last_face
!-----------------------------------------------------------------------------------
  FDofMap = 0
  ReverseNormal(:) = .FALSE.

  IF (PRESENT(FaceIndex)) THEN
    first_face = FaceIndex
    last_face = FaceIndex
  ELSE
    first_face = 1
  END IF

  Mesh => CurrentModel % Solver % Mesh
  Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)
  
  n = Element % TYPE % NumberOfNodes
  GIndexes(1:n) = Element % NodeIndexes(1:n)
  IF( Parallel ) GIndexes(1:n) = Mesh % ParallelInfo % GlobalDOFs(GIndexes(1:n))
  

  SELECT CASE(Element % TYPE % ElementCode / 100)
  CASE(5)
    !
    ! This handles the tetrahedron of Nedelec's second family
    !
    TetraFaceMap(1,:) = (/ 2, 1, 3 /)
    TetraFaceMap(2,:) = (/ 1, 2, 4 /)
    TetraFaceMap(3,:) = (/ 2, 3, 4 /) 
    TetraFaceMap(4,:) = (/ 3, 1, 4 /)

    FaceMap => TetraFaceMap

    IF (.NOT. PRESENT(FaceIndex)) last_face = 4

    DO q=first_face,last_face
      FaceIndices(1:3) = GIndexes(FaceMap(q,1:3))
      IF ( ( FaceIndices(1) < FaceIndices(2) ) .AND. ( FaceIndices(1) < FaceIndices(3) ) ) THEN
        FDofMap(q,1) = 1
        IF (FaceIndices(2) < FaceIndices(3)) THEN
          FDofMap(q,2) = 2
          FDofMap(q,3) = 3                      
        ELSE
          FDofMap(q,2) = 3
          FDofMap(q,3) = 2
        END IF
      ELSE IF ( ( FaceIndices(2) < FaceIndices(1) ) .AND. ( FaceIndices(2) < FaceIndices(3) ) ) THEN
        FDofMap(q,1) = 2
        IF (FaceIndices(1) < FaceIndices(3)) THEN
          FDofMap(q,2) = 1
          FDofMap(q,3) = 3
        ELSE
          FDofMap(q,2) = 3
          FDofMap(q,3) = 1
        END IF
      ELSE
        FDofMap(q,1) = 3
        IF (FaceIndices(1) < FaceIndices(2)) THEN
          FDofMap(q,2) = 1
          FDofMap(q,3) = 2 
        ELSE
          FDofMap(q,2) = 2
          FDofMap(q,3) = 1 
        END IF
      END IF
    END DO

  CASE(8)
    !
    ! Write the face map such that by default the normal points outwards
    ! from the brick:
    !
    BrickFaceMap(1,:) = (/ 2, 1, 4, 3 /)
    BrickFaceMap(2,:) = (/ 5, 6, 7, 8 /)
    BrickFaceMap(3,:) = (/ 1, 2, 6, 5 /)
    BrickFaceMap(4,:) = (/ 2, 3, 7, 6 /)
    BrickFaceMap(5,:) = (/ 3, 4, 8, 7 /)
    BrickFaceMap(6,:) = (/ 4, 1, 5, 8 /)

    FaceMap => BrickFaceMap

    IF (.NOT. PRESENT(FaceIndex)) last_face = 6

    DO q=first_face,last_face
      FaceIndices(1:4) = GIndexes(FaceMap(q,1:4))
    
!      CALL SquareFaceDofsOrdering(I1, I2, D1, D2, FaceIndices(1:4), ReverseSign(q))

      i = 1
      j = 2
      IF ( FaceIndices(i) < FaceIndices(j) ) THEN
        k = i
      ELSE
        k = j
      END IF
      i = 4
      j = 3 
      IF ( FaceIndices(i) < FaceIndices(j) ) THEN
        l = i
      ELSE
        l = j
      END IF
      IF ( FaceIndices(k) > FaceIndices(l) ) THEN
        k = l
      END IF
!      A = k

      SELECT CASE(k)
      CASE(1)
        FDofMap(q,1) = 1
        FDofMap(q,3) = 3
        IF ( FaceIndices(2) < FaceIndices(4) ) THEN
          FDofMap(q,2) = 2
          FDofMap(q,4) = 4
        ELSE
          FDofMap(q,2) = 4
          FDofMap(q,4) = 2
          ReverseNormal(q) = .TRUE.
        END IF
      CASE(2)
        FDofMap(q,2) = 1
        FDofMap(q,4) = 3
        IF ( FaceIndices(3) < FaceIndices(1) ) THEN
          FDofMap(q,1) = 4
          FDofMap(q,3) = 2
        ELSE
          FDofMap(q,1) = 2
          FDofMap(q,3) = 4
          ReverseNormal(q) = .TRUE.
        END IF
      CASE(3)
        FDofMap(q,3) = 1
        FDofMap(q,1) = 3
        IF ( FaceIndices(4) < FaceIndices(2) ) THEN
          FDofMap(q,2) = 4
          FDofMap(q,4) = 2
        ELSE
          FDofMap(q,2) = 2
          FDofMap(q,4) = 4
          ReverseNormal(q) = .TRUE.
        END IF
      CASE(4)
        FDofMap(q,4) = 1
        FDofMap(q,2) = 3
        IF ( FaceIndices(1) < FaceIndices(3) ) THEN
          FDofMap(q,1) = 2
          FDofMap(q,3) = 4
        ELSE
          FDofMap(q,1) = 4
          FDofMap(q,3) = 2
          ReverseNormal(q) = .TRUE.
        END IF
      CASE DEFAULT
        CALL Fatal('ElementDescription::FaceElementBasisOrdering','Erratic square face Indices')
      END SELECT

    END DO

    IF (PRESENT(ReverseSign)) ReverseSign(1:6) = ReverseNormal(1:6)

  CASE DEFAULT
    CALL Fatal('FaceElementBasisOrdering', 'Unsupported element family')
  END SELECT
!-----------------------------------------------------------------------------------
END SUBROUTINE FaceElementBasisOrdering
!-----------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Here the given element can be supposed to be some face of its parent element.
!> The index of the face in reference to the parent element and pointer
!> to the face are returned. The given element and the face returned are thus
!> representations of the same entity but they may still be indexed differently.
!------------------------------------------------------------------------------
SUBROUTINE PickActiveFace(Mesh, Parent, Element, Face, ActiveFaceId)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(Mesh_t), INTENT(IN) :: Mesh  
  TYPE(Element_t), INTENT(IN) :: Parent, Element
  TYPE(Element_t), POINTER, INTENT(OUT) :: Face
  INTEGER, INTENT(OUT) :: ActiveFaceId
!------------------------------------------------------------------------------
  INTEGER :: matches, k, l
!------------------------------------------------------------------------------
  SELECT CASE(Element % TYPE % ElementCode / 100)
  CASE(2)
    IF ( ASSOCIATED(Parent % EdgeIndexes) ) THEN
      DO ActiveFaceId=1,Parent % TYPE % NumberOfEdges
        Face => Mesh % Edges(Parent % EdgeIndexes(ActiveFaceId))
        matches = 0
        DO k=1,Element % TYPE % NumberOfNodes
          DO l=1,Face % TYPE % NumberOfNodes
            IF (Element % NodeIndexes(k) == Face % NodeIndexes(l)) &
                matches=matches+1
          END DO
        END DO
        IF (matches==Element % TYPE % NumberOfNodes) EXIT
      END DO
    ELSE
      matches = 0
    END IF
  CASE(3,4)
    IF ( ASSOCIATED(Parent % FaceIndexes) ) THEN
      DO ActiveFaceId=1,Parent % TYPE % NumberOfFaces
        Face => Mesh % Faces(Parent % FaceIndexes(ActiveFaceId))
        IF ((Element % TYPE % ElementCode / 100) /= (Face % TYPE % ElementCode / 100)) CYCLE
        matches = 0
        DO k=1,Element % TYPE % NumberOfNodes
          DO l=1,Face % TYPE % NumberOfNodes
            IF (Element % NodeIndexes(k) == Face % NodeIndexes(l)) &
                matches=matches+1
          END DO
        END DO
        IF (matches == Element % TYPE % NumberOfNodes ) EXIT
      END DO
    ELSE
      matches = 0
    END IF
  CASE DEFAULT
    CALL Fatal('PickActiveFace', 'Element variable is of a wrong dimension')
  END SELECT

  IF (matches /= Element % TYPE % NumberOfNodes) THEN
    Face => NULL()
    ActiveFaceId = 0
    CALL Warn('PickActiveFace', 'The element is not a face of given parent')
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE PickActiveFace
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Perform the cross product of two vectors
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
!>  Return H(curl)-conforming edge element basis function values and their Curl  
!>  with respect to the reference element coordinates at a given point on the
!>  reference element. Here the basis for a real element K is constructed by  
!>  transforming the basis functions defined on the reference element k via a version
!>  of the Piola transformation designed for functions in H(curl). This construction
!>  differs from the approach taken in the alternate subroutine GetEdgeBasis, which
!>  does not make reference to the Piola transformation and hence may have limitations
!>  in its extendability. The data for performing the Piola transformation is also returned.
!>  Note that the reference element is chosen as in the p-approximation so that
!>  the reference element edges/faces have the same length/area. This choice simplifies
!>  the associated assembly procedure.
!>     With giving the optional argument ApplyPiolaTransform = .TRUE., this function
!>  also performs the Piola transform, so that the basis functions and their spatial
!>  curl as defined on the physical element are returned.
!>     In the lowest-order case this function returns the basis functions belonging
!>  to the optimal family which is not subject to degradation of convergence on
!>  meshes consisting of non-affine physical elements. The second-order elements
!>  are members of the Nedelec's first family and are constructed in a hierarchic
!>  fashion (the lowest-order basis functions give a partial construction of
!>  the second-order basis).
!---------------------------------------------------------------------------------
     FUNCTION EdgeElementInfo( Element, Nodes, u, v, w, F, G, detF, &
          Basis, EdgeBasis, RotBasis, dBasisdx, SecondFamily, BasisDegree, &
          ApplyPiolaTransform, ReadyEdgeBasis, ReadyRotBasis, &
          TangentialTrMapping, GradientVersion) RESULT(stat)
!------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element        !< Element structure
       TYPE(Nodes_t) :: Nodes                    !< Data corresponding to the classic element nodes
       REAL(KIND=dp) :: u                        !< 1st reference element coordinate at which the basis functions are evaluated
       REAL(KIND=dp) :: v                        !< 2nd local coordinate
       REAL(KIND=dp) :: w                        !< 3rd local coordinate
       REAL(KIND=dp), OPTIONAL :: F(3,3)         !< The gradient F=Grad f, with f the element map f:k->K
       REAL(KIND=dp), OPTIONAL :: G(3,3)         !< The transpose of the inverse of the gradient F
       REAL(KIND=dp) :: detF                     !< The determinant of the gradient matrix F
       REAL(KIND=dp) :: Basis(:)                 !< H1-conforming basis functions evaluated at (u,v,w)
       REAL(KIND=dp) :: EdgeBasis(:,:)           !< The basis functions b spanning the reference element space
       REAL(KIND=dp), OPTIONAL :: RotBasis(:,:)  !< The Curl of the edge basis functions with respect to the local coordinates
       REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)  !< The first derivatives of the H1-conforming basis functions at (u,v,w)
       LOGICAL, OPTIONAL :: SecondFamily         !< If .TRUE., a Nedelec basis of the second kind is returned (only simplicial elements)
       INTEGER, OPTIONAL :: BasisDegree          !< The approximation degree 2 (or even 3 in some cases) is also supported
       LOGICAL, OPTIONAL :: ApplyPiolaTransform  !< If  .TRUE., perform the Piola transform so that, instead of b
                                                 !< and Curl b, return  B(f(p)) and (curl B)(f(p)) with B(x) the basis 
                                                 !< functions on the physical element and curl the spatial curl operator.
                                                 !< In this case the absolute value of detF is returned.
       REAL(KIND=dp), OPTIONAL :: ReadyEdgeBasis(:,:) !< A pretabulated edge basis function can be given
       REAL(KIND=dp), OPTIONAL :: ReadyRotBasis(:,:)  !< The preretabulated Curl of the edge basis function
       LOGICAL, OPTIONAL :: TangentialTrMapping  !< To return b x n, with n=(0,0,1) the normal to the 2D reference element.
                                                 !< The Piola transform is then the usual div-conforming version.    
       LOGICAL, OPTIONAL :: GradientVersion      !< Use an alternate basis of the first kind, lacking support for pyramids and bricks
       LOGICAL :: Stat                           !< .FALSE. for a degenerate element
!-----------------------------------------------------------------------------------------------------------------
!      Local variables
!------------------------------------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh
       TYPE(Element_t), POINTER :: Parent, Face, pElement
       INTEGER :: n, dim, cdim, q, i, j, k, l, A, I1, I2, I3, I4, FaceIndices(4), A0, B0, C0, J1, J2, J3
       REAL(KIND=dp) :: dLbasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3), WorkBasis(4,3), WorkCurlBasis(4,3)
       REAL(KIND=dp) :: D1, D2, B(3), curlB(3), GT(3,3), LG(3,3), LF(3,3)
       REAL(KIND=dp) :: ElmMetric(3,3), detJ, CurlBasis(54,3)
       REAL(KIND=dp) :: t(3), s(3), v1, v2, v3, h1, h2, h3, dh1, dh2, dh3, grad(2), grad_i(2), grad_j(2), grad_k(2)
       REAL(KIND=dp) :: LBasis(Element % TYPE % NumberOfNodes), Beta(4), EdgeSign(16)
       REAL(KIND=dp) :: fs1, fs2
       REAL(KIND=dp) :: sfun, tfun, hfun, gfun, bfun
       REAL(KIND=dp) :: grad_sfun(3), grad_tfun(3), grad_hfun(3), grad_gfun(3)
       REAL(KIND=dp) :: svec(3), tvec(3), hvec(3), gvec(3)
       REAL(KIND=dp) :: grad_svec(3,3), grad_tvec(3,3), grad_hvec(3,3), grad_gvec(3,3)
       REAL(KIND=dp) :: WorkWeight(2), grad_weight(2,1:3)
       REAL(KIND=dp) :: Wrk(4,3), WrkCurl(4,3)
       LOGICAL :: Create2ndKindBasis, PerformPiolaTransform, UsePretabulatedBasis, Parallel
       LOGICAL :: SecondOrder, ThirdOrder, ApplyTraceMapping, Found
       LOGICAL :: ReverseSign(4)
       LOGICAL :: ScaleFaceBasis, RedefineFaceBasis, UseWForms = .false.
       LOGICAL :: GradVersion
       LOGICAL :: InformAboutWForms = .TRUE.
       INTEGER, POINTER :: EdgeMap(:,:)
       INTEGER :: TriangleFaceMap(3), SquareFaceMap(4), BrickFaceMap(6,4), PrismSquareFaceMap(3,4), DOFs, GIndexes(27)
       INTEGER :: ActiveFaceId, EDOFs, FDOFs, BDOFs
!----------------------------------------------------------------------------------------------------------
       IF (InformAboutWForms) THEN
         IF (UseWForms) CALL Info('EdgeElementInfo', 'Expressing basis by using Whitney forms')
         InformAboutWForms = .FALSE.
       END IF
       RedefineFaceBasis = .TRUE. ! Left as an emergency switch to revert to the original (ill-conditioned) basis
       ScaleFaceBasis = .TRUE.
       fs1 = 28.0d0
       fs2 = 84.0d0
       
       Mesh => CurrentModel % Solver % Mesh
       Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)

       stat = .TRUE.
       Basis = 0.0d0
       EdgeBasis = 0.0d0
       WorkBasis = 0.0d0
       CurlBasis = 0.0d0
       LG = 0.0d0
       !--------------------------------------------------------------------------------------------
       ! Check whether ready edge basis function values are available to reduce computation.
       ! If they are available, this function is used primarily to obtain the Piola transformation.
       !--------------------------------------------------------------------------------------------
       UsePretabulatedBasis = .FALSE.
       IF ( PRESENT(ReadyEdgeBasis) .AND. PRESENT(ReadyRotBasis) ) UsePretabulatedBasis = .TRUE.
       !------------------------------------------------------------------------------------------
       ! Check whether the Nedelec basis functions of the second kind or higher order basis
       ! functions should be created and whether the Piola transform is already applied within 
       ! this function.
       !------------------------------------------------------------------------------------------
       Create2ndKindBasis = .FALSE.
       IF ( PRESENT(SecondFamily) ) Create2ndKindBasis = SecondFamily
       IF (Create2ndKindBasis .AND. .NOT.(Element % TYPE % ElementCode / 100 == 2 .OR. &
           Element % TYPE % ElementCode / 100 == 3 .OR. &
           Element % TYPE % ElementCode / 100 == 5)) THEN
         CALL Fatal('EdgeElementInfo', 'Second Kind Basis = True is not supported for the given element shape')
       END IF
       
       SecondOrder = .FALSE.
       ThirdOrder = .FALSE.
       IF ( PRESENT(BasisDegree) ) THEN
         SecondOrder = BasisDegree == 2
         IF (.NOT. SecondOrder) ThirdOrder = BasisDegree == 3 
       END IF
       PerformPiolaTransform = .FALSE.
       IF ( PRESENT(ApplyPiolaTransform) ) PerformPiolaTransform = ApplyPiolaTransform
       
       ApplyTraceMapping = .FALSE.
       IF ( PRESENT(TangentialTrMapping) ) ApplyTraceMapping = TangentialTrMapping

       GradVersion = .FALSE.
       IF ( PRESENT(GradientVersion) ) GradVersion = GradientVersion
       IF (GradVersion .AND. .NOT.(Element % TYPE % ElementCode / 100 == 2 .OR. &
           Element % TYPE % ElementCode / 100 == 3 .OR. &
           Element % TYPE % ElementCode / 100 == 4 .OR. &
           Element % TYPE % ElementCode / 100 == 5 .OR. &
           Element % TYPE % ElementCode / 100 == 7)) THEN
         CALL Fatal('EdgeElementInfo', 'Gradient Basis Functions = True is not supported for the given element shape')
       END IF
           
       !-------------------------------------------------------------------------------------------
       dLbasisdx = 0.0d0      
       n = Element % TYPE % NumberOfNodes
       dim = Element % TYPE % DIMENSION
       cdim = CoordinateSystemDimension()

       IF ( Element % TYPE % ElementCode == 101 ) THEN
         detF = 1.0d0
         Basis(1) = 1.0d0
         IF ( PRESENT(dBasisdx) ) dBasisdx(1,:) = 0.0d0
         RETURN
       END IF

       GIndexes(1:n) = Element % NodeIndexes(1:n)
       IF( Parallel ) GIndexes(1:n) = Mesh % ParallelInfo % GlobalDOFs(GIndexes(1:n))
            
       !IF (cdim == 3 .AND. dim==1) THEN
       !  CALL Warn('EdgeElementInfo', 'Traces of 2-D edge elements have not been implemented yet')
       !  RETURN
       !END IF

       !-----------------------------------------------------------------------
       ! The standard nodal basis functions on the reference element and
       ! their derivatives with respect to the local coordinates. These define 
       ! the mapping of the reference element to an actual element on the background 
       ! mesh but are not the basis functions for the edge element approximation.
       ! Remark: Using reference elements having the edges of the same length
       ! simplifies the implementation of element assembly procedures.
       !-----------------------------------------------------------------------
       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(2)
         IF (SecondOrder .AND. n==3) CALL Fatal('EdgeElementInfo', &
             'The lowest-order background mesh needed for trace evaluation over an edge')
         IF (Create2ndKindBasis) CALL Fatal('EdgeElementInfo', &
             'Traces of 2-D edge elements (the 2nd family) have not been implemented yet')
         IF (SecondOrder) THEN
           DOFs = 2
         ELSE
           DOFs = 1
         END IF
         DO q=1,2
           Basis(q) = LineNodalPBasis(q, u)
           dLBasisdx(q,1) = dLineNodalPBasis(q, u)
         END DO
       CASE(3)
         IF (SecondOrder .OR. ThirdOrder) THEN
           ! DOFs is the number of H(curl)-conforming basis functions:
           IF (SecondOrder) THEN
             IF (Create2ndKindBasis) THEN
               DOFs = 12
             ELSE
               DOFs = 8
             END IF
             IF (.NOT.(n==3 .OR. n==6)) CALL Fatal('EdgeElementInfo', 'A 3-node or 6-node background element expected')
           ELSE
             IF (Create2ndKindBasis) THEN
               DOFs = 20
             ELSE
               DOFs = 15
             END IF
             IF (.NOT. n==3) CALL Fatal('EdgeElementInfo', 'A 3-node background element expected')
           END IF
             
           IF (n == 6) THEN
             ! Here the element of the background mesh is of type 306.
             ! The Lagrange interpolation basis on the p-approximation reference element:
             Basis(1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v))/6.0d0
             dLBasisdx(1,1) = -0.5d0 + u + v/Sqrt(3.0d0)
             dLBasisdx(1,2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
             Basis(2) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2.0d0*Sqrt(3.0d0)*v))/6.0d0
             dLBasisdx(2,1) = 0.5d0 + u - v/Sqrt(3.d0)
             dLBasisdx(2,2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
             Basis(3) = (v*(-Sqrt(3.0d0) + 2.0d0*v))/3.0d0
             dLBasisdx(3,1) = 0.0d0
             dLBasisdx(3,2) =  -(1.0d0/Sqrt(3.0d0)) + (4.0d0*v)/3.0d0
             Basis(4) = (3.0d0 - 3.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + v**2)/3.0d0
             dLBasisdx(4,1) = -2.0d0*u
             dLBasisdx(4,2) = (-2.0d0*(Sqrt(3.0d0) - v))/3.0d0
             Basis(5) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
             dLBasisdx(5,1) =  (2.0d0*v)/Sqrt(3.0d0)
             dLBasisdx(5,2) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2.0d0*v))/3.0d0
             Basis(6) = (-2.0d0*v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0           
             dLBasisdx(6,1) = (-2.0d0*v)/Sqrt(3.0d0)
             dLBasisdx(6,2) = (-2.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 2.0d0*v))/3.0d0
           ELSE
             ! Here the element of the background mesh is of type 303:
             DO q=1,3
               Basis(q) = TriangleNodalPBasis(q, u, v)
               dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
             END DO
           END IF
         ELSE
           DO q=1,n
             Basis(q) = TriangleNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dTriangleNodalPBasis(q, u, v) 
           END DO
           IF (Create2ndKindBasis) THEN
             DOFs = 6
           ELSE
             DOFs = 3
           END IF
         END IF
       CASE(4)
         IF (SecondOrder) THEN
           ! The second-order quad from the Nedelec's first family: affine physical elements may be needed
           DOFs = 12
         ELSE
           ! The lowest-order quad from the optimal family (ABF_0)
           DOFs = 6
         END IF
         IF (n>4) THEN
           ! Here the background mesh is supposed to be of type 408/409
           CALL NodalBasisFunctions2D(Basis, Element, u, v)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, w)
         ELSE
           ! Here the background mesh is of type 404           
           DO q=1,4
             Basis(q) = QuadNodalPBasis(q, u, v)
             dLBasisdx(q,1:2) = dQuadNodalPBasis(q, u, v) 
           END DO
         END IF
       CASE(5)
         IF (SecondOrder .OR. ThirdOrder) THEN
           IF (SecondOrder) THEN
             IF (Create2ndKindBasis) THEN
               DOFs = 30
             ELSE
               DOFs = 20
             END IF
           ELSE
             IF (Create2ndKindBasis) THEN
               CALL Fatal('EdgeElementInfo', 'A cubic element of the 2nd kind is not yet available')
               DOFs = 60
             ELSE
               DOFs = 45
             END IF
             IF (.NOT. n==4) CALL Fatal('EdgeElementInfo', 'A 4-node background element expected')             
           END IF
           
           IF (n == 10) THEN
             ! Here the element of the background mesh is of type 510.
             ! The Lagrange interpolation basis on the p-approximation reference element:
             Basis(1) = (6.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - Sqrt(6.0d0)*w + 2.0d0*Sqrt(2.0d0)*v*w + &
                 w**2 + 2.0d0*u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/12.0d0
             dLBasisdx(1,1) = -0.5d0 + u + v/Sqrt(3.0d0) + w/Sqrt(6.0d0)
             dLBasisdx(1,2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v + Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(1,3) = (-Sqrt(6.0d0) + 2.0d0*Sqrt(6.0d0)*u + 2.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(2) = (6.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - Sqrt(6.0d0)*w + 2.0d0*Sqrt(2.0d0)*v*w + &
                 w**2 - 2.0d0*u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/12.0d0
             dLBasisdx(2,1) = 0.5d0 + u - v/Sqrt(3.0d0) - w/Sqrt(6.0d0)
             dLBasisdx(2,2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v + Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(2,3) = (-Sqrt(6.0d0) - 2.0d0*Sqrt(6.0d0)*u + 2.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(3) =  (8.0d0*v**2 + w*(Sqrt(6.0d0) + w) - 4.0d0*v*(Sqrt(3.0d0) + Sqrt(2.0d0)*w))/12.0d0
             dLBasisdx(3,1) = 0.0d0
             dLBasisdx(3,2) = (-Sqrt(3.0d0) + 4.0d0*v - Sqrt(2.0d0)*w)/3.0d0
             dLBasisdx(3,3) = (Sqrt(6.0d0) - 4.0d0*Sqrt(2.0d0)*v + 2.0d0*w)/12.0d0
             Basis(4) = (w*(-Sqrt(6.0d0) + 3.0d0*w))/4.0d0
             dLBasisdx(4,1) = 0.0d0
             dLBasisdx(4,2) = 0.0d0
             dLBasisdx(4,3) = (-Sqrt(6.0d0) + 6.0d0*w)/4.0d0
             Basis(5) =  (6.0d0 - 6.0d0*u**2 - 4.0d0*Sqrt(3.0d0)*v + 2.0d0*v**2 - 2.0d0*Sqrt(6.0d0)*w + &
                 2.0d0*Sqrt(2.0d0)*v*w + w**2)/6.0d0
             dLBasisdx(5,1) = -2.0d0*u
             dLBasisdx(5,2) = (-2.0d0*Sqrt(3.0d0) + 2.0d0*v + Sqrt(2.0d0)*w)/3.0d0
             dLBasisdx(5,3) = (-Sqrt(6.0d0) + Sqrt(2.0d0)*v + w)/3.0d0
             Basis(6) =  (-4.0d0*v**2 + w*(-Sqrt(6.0d0) - Sqrt(6.0d0)*u + w) + v*(4.0d0*Sqrt(3.0d0) + &
                 4.0d0*Sqrt(3.0d0)*u - Sqrt(2.0d0)*w))/6.0d0
             dLBasisdx(6,1) = (2.0d0*v)/Sqrt(3.0d0) - w/Sqrt(6.0d0)
             dLBasisdx(6,2) = (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 8.0d0*v - Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(6,3) = (-Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v + 2.0d0*w)/6.0d0
             Basis(7) =  (-4.0d0*v**2 + w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + w) - &
                 v*(-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + Sqrt(2.0d0)*w))/6.0d0
             dLBasisdx(7,1) = (-2.0d0*v)/Sqrt(3.0d0) + w/Sqrt(6.0d0)
             dLBasisdx(7,2) = (4.0d0*Sqrt(3.0d0) - 4.0d0*Sqrt(3.0d0)*u - 8.0d0*v - Sqrt(2.0d0)*w)/6.0d0
             dLBasisdx(7,3) = (-Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v + 2.0d0*w)/6.0d0
             Basis(8) = -(w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + w))/2.0d0
             dLBasisdx(8,1) = -(Sqrt(1.5d0)*w)
             dLBasisdx(8,2) = -(w/Sqrt(2.0d0))
             dLBasisdx(8,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 2.0d0*w)/2.0d0
             Basis(9) = ((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - w)*w)/2.0d0
             dLBasisdx(9,1) = Sqrt(1.5d0)*w
             dLBasisdx(9,2) = -(w/Sqrt(2.0d0))
             dLBasisdx(9,3) = (Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 2.0d0*w)/2.0d0
             Basis(10) = Sqrt(2.0d0)*v*w - w**2/2.0d0
             dLBasisdx(10,1) = 0.0d0
             dLBasisdx(10,2) = Sqrt(2.0d0)*w
             dLBasisdx(10,3) = Sqrt(2.0d0)*v - w
           ELSE
             ! Here the element of the background mesh is of type 504: 
             DO q=1,4
               Basis(q) = TetraNodalPBasis(q, u, v, w)
               dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
             END DO
           END IF
         ELSE
           DO q=1,n
             Basis(q) = TetraNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dTetraNodalPBasis(q, u, v, w)
           END DO
           IF (Create2ndKindBasis) THEN
             DOFs = 12
           ELSE
             DOFs = 6
           END IF
         END IF
       CASE(6)
         IF (SecondOrder) THEN
           ! The second-order pyramid from the Nedelec's first family
           DOFs = 31
         ELSE
           ! The lowest-order pyramid from the optimal family
           DOFs = 10
         END IF

         IF (n==13) THEN
           ! Here the background mesh is supposed to be of type 613. The difference between the standard
           ! reference element and the p-reference element can be taken into account by a simple scaling:
           CALL NodalBasisFunctions3D(Basis, Element, u, v, sqrt(2.0d0)*w)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, sqrt(2.0d0)*w)
           dLBasisdx(1:n,3) = sqrt(2.0d0) * dLBasisdx(1:n,3)
         ELSE
           ! Background mesh elements of the type 605:
           DO q=1,n
             Basis(q) = PyramidNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dPyramidNodalPBasis(q, u, v, w)
           END DO
         END IF

       CASE(7)
         IF (SecondOrder) THEN
           ! The second-order prism from the Nedelec's first family: affine physical elements may be needed
           DOFs = 36
         ELSE
           ! The lowest-order prism from the optimal family
           DOFs = 15
         END IF

         IF (n==15) THEN
           ! Here the background mesh is of type 715.
           ! The Lagrange interpolation basis on the p-approximation reference element:

           h1 = -0.5d0*w + 0.5d0*w**2
           h2 = 0.5d0*w + 0.5d0*w**2
           h3 = 1.0d0 - w**2
           dh1 = -0.5d0 + w
           dh2 = 0.5d0 + w
           dh3 = -2.0d0 * w
           
           WorkBasis(1,1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2.0d0*Sqrt(3.0d0)*v))/6
           grad(1) = -0.5d0 + u + v/Sqrt(3.0d0)
           grad(2) = (-Sqrt(3.0d0) + 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
           Basis(1) = WorkBasis(1,1) * h1
           dLBasisdx(1,1:2) = grad(1:2) * h1
           dLBasisdx(1,3) = WorkBasis(1,1) * dh1
           Basis(4) = WorkBasis(1,1) * h2
           dLBasisdx(4,1:2) = grad(1:2) * h2
           dLBasisdx(4,3) = WorkBasis(1,1) * dh2
           Basis(13) = WorkBasis(1,1) * h3
           dLBasisdx(13,1:2) = grad(1:2) * h3
           dLBasisdx(13,3) = WorkBasis(1,1) * dh3

           WorkBasis(1,1) = (3.0d0*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2.0d0*Sqrt(3.0d0)*v))/6.0d0
           grad(1) = 0.5d0 + u - v/Sqrt(3.d0)
           grad(2) = (-Sqrt(3.0d0) - 2.0d0*Sqrt(3.0d0)*u + 2.0d0*v)/6.0d0
           Basis(2) = WorkBasis(1,1) * h1
           dLBasisdx(2,1:2) = grad(1:2) * h1
           dLBasisdx(2,3) = WorkBasis(1,1) * dh1
           Basis(5) = WorkBasis(1,1) * h2
           dLBasisdx(5,1:2) = grad(1:2) * h2
           dLBasisdx(5,3) = WorkBasis(1,1) * dh2
           Basis(14) = WorkBasis(1,1) * h3
           dLBasisdx(14,1:2) = grad(1:2) * h3
           dLBasisdx(14,3) = WorkBasis(1,1) * dh3

           WorkBasis(1,1) = (v*(-Sqrt(3.0d0) + 2.0d0*v))/3.0d0
           grad(1) = 0.0d0
           grad(2) = -(1.0d0/Sqrt(3.0d0)) + (4.0d0*v)/3.0d0
           Basis(3) = WorkBasis(1,1) * h1
           dLBasisdx(3,1:2) = grad(1:2) * h1
           dLBasisdx(3,3) = WorkBasis(1,1) * dh1
           Basis(6) = WorkBasis(1,1) * h2
           dLBasisdx(6,1:2) = grad(1:2) * h2
           dLBasisdx(6,3) = WorkBasis(1,1) * dh2
           Basis(15) = WorkBasis(1,1) * h3
           dLBasisdx(15,1:2) = grad(1:2) * h3
           dLBasisdx(15,3) = WorkBasis(1,1) * dh3

           h1 = 0.5d0 * (1.0d0 - w)
           dh1 = -0.5d0
           h2 = 0.5d0 * (1.0d0 + w)
           dh2 = 0.5d0

           WorkBasis(1,1) = (3.0d0 - 3.0d0*u**2 - 2.0d0*Sqrt(3.0d0)*v + v**2)/3.0d0
           grad(1) = -2.0d0*u
           grad(2) = (-2.0d0*(Sqrt(3.0d0) - v))/3.0d0
           Basis(7) = WorkBasis(1,1) * h1
           dLBasisdx(7,1:2) = grad(1:2) * h1
           dLBasisdx(7,3) = WorkBasis(1,1) * dh1
           Basis(10) = WorkBasis(1,1) * h2
           dLBasisdx(10,1:2) = grad(1:2) * h2
           dLBasisdx(10,3) = WorkBasis(1,1) * dh2

           WorkBasis(1,1) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
           grad(1) = (2.0d0*v)/Sqrt(3.0d0)
           grad(2) = (2.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2.0d0*v))/3.0d0
           Basis(8) = WorkBasis(1,1) * h1
           dLBasisdx(8,1:2) = grad(1:2) * h1
           dLBasisdx(8,3) = WorkBasis(1,1) * dh1
           Basis(11) = WorkBasis(1,1) * h2
           dLBasisdx(11,1:2) = grad(1:2) * h2
           dLBasisdx(11,3) = WorkBasis(1,1) * dh2

           WorkBasis(1,1) = (-2.0d0*v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0
           grad(1) = (-2.0d0*v)/Sqrt(3.0d0)
           grad(2) = (-2.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 2.0d0*v))/3.0d0
           Basis(9) = WorkBasis(1,1) * h1
           dLBasisdx(9,1:2) = grad(1:2) * h1
           dLBasisdx(9,3) = WorkBasis(1,1) * dh1
           Basis(12) = WorkBasis(1,1) * h2
           dLBasisdx(12,1:2) = grad(1:2) * h2
           dLBasisdx(12,3) = WorkBasis(1,1) * dh2
         ELSE
           ! Here the background mesh is of type 706
           DO q=1,n
             Basis(q) = WedgeNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dWedgeNodalPBasis(q, u, v, w)
           END DO
         END IF
       CASE(8)
         IF (SecondOrder) THEN
           ! The second-order brick from the Nedelec's first family: affine physical elements may be needed
           DOFs = 54
         ELSE
           ! The lowest-order brick from the optimal family
           DOFs = 27
         END IF
         IF (n>8) THEN
           ! Here the background mesh is supposed to be of type 820/827
           CALL NodalBasisFunctions3D(Basis, Element, u, v, w)
           CALL NodalFirstDerivatives(n, dLBasisdx, Element, u, v, w) 
         ELSE
           ! Here the background mesh is of type 808
           DO q=1,n
             Basis(q) = BrickNodalPBasis(q, u, v, w)
             dLBasisdx(q,1:3) = dBrickNodalPBasis(q, u, v, w)
           END DO
         END IF
       CASE DEFAULT
         CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
       END SELECT

       !-----------------------------------------------------------------------
       ! Get data for performing the Piola transformation...
       !-----------------------------------------------------------------------
       stat = PiolaTransformationData(n, Element, Nodes, LF, detF, dLBasisdx) 
       !------------------------------------------------------------------------
       ! ... in order to define the basis for the element space X(K) via 
       ! applying a version of the Piola transformation as
       !    X(K) = { B | B = F^{-T}(f^{-1}(x)) b(f^{-1}(x)) }
       ! with b giving the edge basis function on the reference element k,
       ! f mapping k to the actual element K, i.e. K = f(k) and F = Grad f. This 
       ! function returns the local basis functions b and their Curl (with respect
       ! to local coordinates) evaluated at the integration point. The effect of 
       ! the Piola transformation need to be considered when integrating, so we 
       ! shall return also the values of F, G=F^{-T} and det F.
       !
       ! It should be noted that the case of 2-D surface elements embedded in
       ! the three-dimensional space is handled as a special case. Then F^{-T}
       ! is replaced by the transpose of the pseudoinverse of F. The Piola 
       ! transformation then maps a 2-component field to a 3-component vector
       ! field which is tangential to the 2-D surface.
       !
       ! The construction of edge element bases could be done in an alternate way for 
       ! triangles and tetrahedra, while the chosen approach has the benefit that
       ! it generalizes to other cases. For example general quadrilaterals may now 
       ! be handled in the same way.
       !---------------------------------------------------------------------------
       IF (cdim == dim) THEN
          SELECT CASE(Element % TYPE % ElementCode / 100)
          CASE(3,4)
             LG(1,1) = 1.0d0/detF * LF(2,2)
             LG(1,2) = -1.0d0/detF * LF(1,2)
             LG(2,1) = -1.0d0/detF * LF(2,1)
             LG(2,2) = 1.0d0/detF * LF(1,1)
          CASE(5,6,7,8)
             CALL InvertMatrix3x3(LF,LG,detF)       
          CASE DEFAULT
             CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
          END SELECT
          LG(1:dim,1:dim) = TRANSPOSE( LG(1:dim,1:dim) )
       END IF

       IF (UsePretabulatedBasis) THEN
         DO i=1,DOFs
           EdgeBasis(i,1:3) = ReadyEdgeBasis(i,1:3)
           CurlBasis(i,1:3) = ReadyRotBasis(i,1:3)
         END DO
       ELSE
         SELECT CASE(Element % TYPE % ElementCode / 100)
         CASE(2)
           !--------------------------------------------------------------
           ! This is a special case to return the tangential components 
           ! trace of 2D elements
           !--------------------------------------------------------------
           !
           ! The sign reversion of basis must be checked via the parent element:
           !
           Parent => Element % BoundaryInfo % Left
           IF (.NOT. ASSOCIATED(Parent)) THEN
             Parent => Element % BoundaryInfo % Right
           END IF

           IF (.NOT. ASSOCIATED(Parent)) THEN
             CALL Warn('EdgeElementInfo', 'cannot create curl-conforming basis functions, zeros returned')
             RETURN
           END IF
           !
           ! Identify the edge representing the element among the edges of 
           ! the parent element:
           !
           pElement => Element 
           CALL PickActiveFace(Mesh, Parent, pElement, Face, ActiveFaceId)
           IF (ActiveFaceId == 0) RETURN
           !
           ! Use the parent element to check whether sign reversions are needed:
           !
           CALL FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)
           
           IF (ReverseSign(ActiveFaceId)) THEN
             EdgeBasis(1,1) = -0.5d0
           ELSE
             EdgeBasis(1,1) = 0.5d0
           END IF
           IF (SecondOrder) THEN
             EdgeBasis(2,1) = 1.5d0 * u
           END IF
           CurlBasis(1:DOFs,:) = 0.0d0

         CASE(3)
           !--------------------------------------------------------------
           ! This branch is for handling triangles. Note that
           ! the global orientation of the edge tangent t is defined such that
           ! t points towards the node that has a larger global index.
           !--------------------------------------------------------------
           EdgeMap => GetEdgeMap(3)
           !EdgeMap => GetEdgeMap(GetElementFamily(Element))

           IF (Create2ndKindBasis .OR. GradVersion .AND. &
               (SecondOrder .OR. ThirdOrder)) THEN

             IF (Create2ndKindBasis) THEN
               ! This construction follows Sun, Lee, Cendes. SIAM J. Sci. Comput. 23(4):1053-1076.
               ! The first basis function associated with an edge is the Whitney form, while 
               ! the second basis function corresponds to a gradient field.

               IF (SecondOrder) THEN
                 EDOFs = 3
                 FDOFs = 3
               ELSE
                 EDOFs = 2
                 FDOFs = 0
               END IF
             ELSE
               !
               ! An alternate first-family basis of degree 2 or 3 for faster solution with iterative methods. 
               !
               IF (SecondOrder) THEN
                 EDOFs = 2
                 FDOFs = 2
               ELSE
                 ! The case of third-order basis 
                 EDOFs = 3
                 FDOFs = 6
               END IF
             END IF
               
             DO k=1,3
               
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)

               CALL EdgeWhitneyComponents2D(WorkBasis(1:2,:), WorkCurlBasis(1:2,:), i, j, u, v)
               
               IF (Create2ndKindBasis .AND. SecondOrder .OR. &
                   GradVersion .AND. ThirdOrder) THEN
                 WorkWeight(1) = 2.0d0*Basis(i) - Basis(j)
                 WorkWeight(2) = 2.0d0*Basis(j) - Basis(i)

                 grad_weight(1,1:2) = 2.0d0*dLBasisdx(i,1:2) - dLBasisdx(j,1:2)
                 grad_weight(2,1:2) = 2.0d0*dLBasisdx(j,1:2) - dLBasisdx(i,1:2)
               END IF

               IF (GIndexes(j) < GIndexes(i)) THEN
                 I1 = 2
                 I2 = 1
               ELSE
                 I1 = 1
                 I2 = 2
               END IF

               DO l=1,EDOFs
                 SELECT CASE(l)
                 CASE(1)
                   sfun = -1.0d0
                   tfun = 1.0d0
                 CASE(2)
                   sfun = 1.0d0
                   tfun = 1.0d0
                 CASE(3)
                   sfun = -WorkWeight(I1)
                   tfun = WorkWeight(I2)
                   grad_sfun(1:2) = -grad_weight(I1,1:2)
                   grad_tfun(1:2) = grad_weight(I2,1:2)
                 CASE DEFAULT
                   CALL Fatal('ElementDescription::EdgeElementInfo','sfun/tfun not defined')
                 END SELECT

                 EdgeBasis(EDOFs*(k-1)+l,1:2) = sfun * WorkBasis(I1,1:2) + tfun * WorkBasis(I2,1:2)
                 CurlBasis(EDOFs*(k-1)+l,3) = sfun * WorkCurlBasis(I1,3) + tfun * WorkCurlBasis(I2,3)

                 IF (l > 2) THEN
                   CurlBasis(EDOFs*(k-1)+l,3) = CurlBasis(EDOFs*(k-1)+l,3) + &
                       grad_sfun(1)*WorkBasis(I1,2) + grad_tfun(1)*WorkBasis(I2,2) - &
                       grad_sfun(2)*WorkBasis(I1,1) - grad_tfun(2)*WorkBasis(I2,1)
                 END IF
               END DO
             END DO

             ! The basis functions associated with the faces for higher-order cases
             IF (FDOFs > 0) THEN
               TriangleFaceMap(:) = (/ 1,2,3 /)

               CALL FaceWhitneyComponents2D(WorkBasis(1:3,:), WorkCurlBasis(1:3,:), u, v)
               
               ! Create permutation:
               FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
               CALL TriangleFaceDofsOrdering2nd(I1,I2,I3,FaceIndices(1:3))

               ! Create the basis:
               DO l=1,FDOFs

                 SELECT CASE(l)
                 CASE(1)
                   sfun = 1.0d0
                   tfun = 1.0d0
                   hfun = -2.0d0
                 CASE(2)
                   sfun = 1.0d0
                   tfun = -1.0d0
                   hfun = 0.0d0
                 CASE(3)
                   sfun = 1.0d0
                   tfun = 1.0d0
                   hfun = 1.0d0                                  
                 CASE(4)
                   sfun = Basis(I2) - Basis(I3)
                   tfun = Basis(I3) - Basis(I1)
                   hfun = Basis(I1) - Basis(I2)

                   grad_sfun(1:2) = dLBasisdx(I2,1:2) - dLBasisdx(I3,1:2)
                   grad_tfun(1:2) = dLBasisdx(I3,1:2) - dLBasisdx(I1,1:2)
                   grad_hfun(1:2) = dLBasisdx(I1,1:2) - dLBasisdx(I2,1:2)
                 CASE(5)
                   sfun = 393.0d0 * Basis(I1) + 80.0d0 * Basis(I2) - 212.0d0 * Basis(I3)
                   tfun = -393.0d0 * Basis(I2) - 80.0d0 * Basis(I1) + 212.0d0 * Basis(I3)
                   hfun = -313.0d0 * Basis(I1) + 313.0d0 * Basis(I2)

                   grad_sfun(1:2) = 393.0d0 * dLBasisdx(I1,1:2) + 80.0d0 * dLBasisdx(I2,1:2) - 212.0d0 * dLBasisdx(I3,1:2)
                   grad_tfun(1:2) = -393.0d0 * dLBasisdx(I2,1:2) - 80.0d0 * dLBasisdx(I1,1:2) + 212.0d0 * dLBasisdx(I3,1:2)
                   grad_hfun(1:2) = -313.0d0 * dLBasisdx(I1,1:2) + 313.0d0 * dLBasisdx(I2,1:2)
                 CASE(6)
                   sfun = -131.0d0 * Basis(I1) + 168.0d0 * Basis(I2) - 124.0d0 * Basis(I3)
                   tfun = -131.0d0 * Basis(I2) + 168.0d0 * Basis(I1) - 124.0d0 * Basis(I3)
                   hfun = -37.0d0 * Basis(I1) - 37.0d0 * Basis(I2) + 248.0d0 * Basis(I3)

                   grad_sfun(1:2) = -131.0d0 * dLBasisdx(I1,1:2) + 168.0d0 * dLBasisdx(I2,1:2) - 124.0d0 * dLBasisdx(I3,1:2)
                   grad_tfun(1:2) = -131.0d0 * dLBasisdx(I2,1:2) + 168.0d0 * dLBasisdx(I1,1:2) - 124.0d0 * dLBasisdx(I3,1:2)
                   grad_hfun(1:2) = -37.0d0 * dLBasisdx(I1,1:2) - 37.0d0 * dLBasisdx(I2,1:2) + 248.0d0 * dLBasisdx(I3,1:2)                 
                 END SELECT

                 EdgeBasis(3*EDOFs + l,1:2) = sfun * WorkBasis(I1,1:2) + tfun * WorkBasis(I2,1:2) + &
                     hfun * WorkBasis(I3,1:2)
                 CurlBasis(3*EDOFs + l,3) = sfun * WorkCurlBasis(I1,3) + tfun * WorkCurlBasis(I2,3) + &
                     hfun * WorkCurlBasis(I3,3)

                 IF (l > 3) THEN
                   CurlBasis(3*EDOFs+l,3) = CurlBasis(3*EDOFs+l,3) + &
                       grad_sfun(1)*WorkBasis(I1,2) + grad_tfun(1)*WorkBasis(I2,2) + grad_hfun(1)*WorkBasis(I3,2) - &
                       grad_sfun(2)*WorkBasis(I1,1) - grad_tfun(2)*WorkBasis(I2,1) - grad_hfun(2)*WorkBasis(I3,1)
                 END IF

               END DO
             END IF

           ELSE
             
             !------------------------------------------------------------
             ! The optimal/Nedelec basis functions of the first kind. We employ
             ! a hierarchic basis, so the lowest-order basis functions are
             ! also utilized in the construction of the second-order basis. 
             ! First the edge 12 ...
             !------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0
             EdgeBasis(1,2) = u/(2.0d0*Sqrt(3.0d0))
             CurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,3) = -CurlBasis(1,3)
             END IF
             IF (SecondOrder) THEN
               EdgeBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0
               EdgeBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
               CurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0                     
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the edge 23:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 3
               EdgeBasis(4,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0
               EdgeBasis(4,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0
               CurlBasis(4,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0
             ELSE
               k = 2
             END IF
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(k,1) = -v/(2.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (1 + u)/(2.0d0*Sqrt(3.0d0))
             CurlBasis(k,3) =  1.0d0/Sqrt(3.0d0)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,3) = -CurlBasis(k,3)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the edge 31:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 5
               EdgeBasis(6,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
               EdgeBasis(6,2) = -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0
               CurlBasis(6,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0                     
             ELSE
               k = 3
             END IF
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(k,1) = -v/(2.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0))
             CurlBasis(k,3) = 1.0d0/Sqrt(3.0d0)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,3) = -CurlBasis(k,3)
             END IF

             IF (SecondOrder) THEN
               !-------------------------------------------------
               ! Two basis functions defined on the face 123:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 1,2,3 /)          
               FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               WorkBasis(1,1) = ((Sqrt(3.0d0) - v)*v)/6.0d0
               WorkBasis(1,2) = (u*v)/6.0d0
               WorkCurlBasis(1,3) = (-Sqrt(3.0d0) + 3.0d0*v)/6.0d0
               WorkBasis(2,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0))
               WorkBasis(2,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(2,3) =(-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0
               WorkBasis(3,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkBasis(3,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(3,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

               IF (RedefineFaceBasis) THEN
                 EdgeBasis(7,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
                 CurlBasis(7,3) = 0.5d0 * D1 * WorkCurlBasis(I1,3) + 0.5d0 * D2 * WorkCurlBasis(I2,3)
                 EdgeBasis(8,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
                 CurlBasis(8,3) = 0.5d0 * D2 * WorkCurlBasis(I2,3) - 0.5d0 * D1 * WorkCurlBasis(I1,3)
               ELSE
                 EdgeBasis(7,:) = D1 * WorkBasis(I1,:)
                 CurlBasis(7,3) = D1 * WorkCurlBasis(I1,3)
                 EdgeBasis(8,:) = D2 * WorkBasis(I2,:)
                 CurlBasis(8,3) = D2 * WorkCurlBasis(I2,3)  
               END IF
               
               ! Finally, scale to reduce ill-conditioning:
               IF (ScaleFaceBasis) THEN
                 EdgeBasis(7,:) = sqrt(fs1) * EdgeBasis(7,:)
                 EdgeBasis(8,:) = sqrt(fs2) * EdgeBasis(8,:)
                 CurlBasis(7,3) = sqrt(fs1) * CurlBasis(7,3)
                 CurlBasis(8,3) = sqrt(fs2) * CurlBasis(8,3)
               END IF
             END IF
           END IF

         CASE(4)
           !--------------------------------------------------------------
           ! This branch is for handling quadrilaterals
           !--------------------------------------------------------------
           EdgeMap => GetEdgeMap(4)
           IF (SecondOrder) THEN
             IF (GradVersion) THEN
               !
               ! An alternate basis which is compatible with the basis originally constructed for
               ! simplicial elements when GradientVersion = .TRUE.. Here the basis functions are  
               ! defined in terms of the Lobatto shape functions Phi(k,.) and the Legendre
               ! polynomials LegendreP(1,.)
               !
               EDOFs = 2
               FDOFs = 4

               DO k=1,4
                 i = EdgeMap(k,1)
                 j = EdgeMap(k,2)

                 SELECT CASE(k)
                 CASE(1)
                   ! (u,v) -> 1/2 * P0 * L1(v) * e1
                   I1 = 1
                   WorkBasis(1,1) = 0.1D1 / 0.4D1 - v / 0.4D1
                   WorkCurlBasis(1,3) = 0.1D1 / 0.4D1

                   ! (u,v) -> grad( -1/sqrt(6) * phi_2(u) * L1(v) ) so that the tangential
                   ! components trace on the edge is given by -1/2 P1(u) * e1 = -1/2 u * e1

                   WorkBasis(2,1) = -1.0d0/sqrt(6.0d0) * dPhi(2,u) * LineNodalPBasis(1,v)
                   WorkBasis(2,2) = -1.0d0/sqrt(6.0d0) * Phi(2,u) * dLineNodalPBasis(1,v)
                   WorkCurlBasis(2,3) = 0.0d0

                 CASE(2)
                   ! (u,v) -> 1/2 * P0 * L2(u) * e2
                   I1 = 2
                   WorkBasis(1,2) = 0.1D1 / 0.4D1 + u / 0.4D1
                   WorkCurlBasis(1,3) = 0.1D1 / 0.4D1

                   ! (u,v) -> grad( -1/sqrt(6) * phi_2(v) * L2(u) ) so that the tangential
                   ! components trace on the edge is given by -1/2 P1(v) * e2 = -1/2 v * e2

                   WorkBasis(2,1) = -1.0d0/sqrt(6.0d0) * Phi(2,v) * dLineNodalPBasis(2,u)
                   WorkBasis(2,2) = -1.0d0/sqrt(6.0d0) * dPhi(2,v) * LineNodalPBasis(2,u)
                   WorkCurlBasis(2,3) = 0.0d0

                 CASE(3)
                   ! (u,v) -> -1/2 * P0 * L2(v) * e1
                   I1 = 1
                   WorkBasis(1,1) = -0.1D1 / 0.4D1 - v / 0.4D1
                   WorkCurlBasis(1,3) = 0.1D1 / 0.4D1

                   ! (u,v) -> grad( -1/sqrt(6) * phi_2(u) * L2(v) ) ; cf. CASE(1)
                   WorkBasis(2,1) = -1.0d0/sqrt(6.0d0) * dPhi(2,u) * LineNodalPBasis(2,v)
                   WorkBasis(2,2) = -1.0d0/sqrt(6.0d0) * Phi(2,u) * dLineNodalPBasis(2,v)
                   WorkCurlBasis(2,3) = 0.0d0

                 CASE(4)
                   ! (u,v) -> -1/2 * P0 * L1(u) * e2
                   I1 = 2
                   WorkBasis(1,2) = -0.1D1 / 0.4D1 + u / 0.4D1
                   WorkCurlBasis(1,3) = 0.1D1 / 0.4D1

                   ! (u,v) -> grad( -1/sqrt(6) * phi_2(v) * L1(u) ) ; cf. CASE(2)

                   WorkBasis(2,1) = -1.0d0/sqrt(6.0d0) * Phi(2,v) * dLineNodalPBasis(1,u)
                   WorkBasis(2,2) = -1.0d0/sqrt(6.0d0) * dPhi(2,v) * LineNodalPBasis(1,u)
                   WorkCurlBasis(2,3) = 0.0d0

                 END SELECT

                 IF (GIndexes(j) < GIndexes(i)) THEN
                   WorkBasis(1,I1) = -WorkBasis(1,I1)
                   WorkCurlBasis(1,3) = -WorkCurlBasis(1,3)
                 END IF

                 EdgeBasis(EDOFs*(k-1)+1,I1) = WorkBasis(1,I1)
                 CurlBasis(EDOFs*(k-1)+1,3) = WorkCurlBasis(1,3)

                 !DO l=2,EDOFs 
                 EdgeBasis(EDOFs*(k-1)+2,1:2) = WorkBasis(2,1:2)
                 CurlBasis(EDOFs*(k-1)+2,3) = WorkCurlBasis(2,3)
                 !END DO
               END DO

               SquareFaceMap(:) = (/ 1,2,3,4 /)
               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               WorkBasis(1:4,1:2) = 0.0d0

               ! (u,v) ->  P0/2 * (-2) * sqrt(2.0d0/3.0d0) * phi_2(v) e1
               !        =  1/2 P0 * 4 L_1(v) L_2(v) e1
               WorkBasis(1,1) = -sqrt(2.0d0/3.0d0) * Phi(2,v)
               WorkCurlBasis(1,3) = sqrt(2.0d0/3.0d0) * dPhi(2,v)

               ! (u,v) ->  (-2) * sqrt(2.0d0/3.0d0) * phi_2(u) * P0/2 e2
               !        =  1/2 P0 * 4 L_1(u) L_2(u) e2               
               WorkBasis(2,2) = -sqrt(2.0d0/3.0d0) * Phi(2,u)
               WorkCurlBasis(2,3) = -sqrt(2.0d0/3.0d0) * dPhi(2,u)

               ! (u,v) ->  D_u [-1/sqrt(6) * phi_2(u) * (-2) * sqrt(2/3) * phi_2(v) ] e_1
               !        = -1/2 P1(u) * [(-2) * sqrt(2/3) * phi_2(v)] e_1
               !        = -1/2 P1(u) * [4 L_1(v) L_2(v)] e_1
               WorkBasis(3,1) = -1.0d0/2.0d0 * LegendreP(1,u) * (-2.0d0) * sqrt(2.0d0/3.0d0) * Phi(2,v)
               WorkCurlBasis(3,3) = 1.0d0/2.0d0 * LegendreP(1,u) * (-2.0d0) * sqrt(2.0d0/3.0d0) * dPhi(2,v)

               ! (u,v) -> D_v [-1/sqrt(6) * phi_2(v) * (-2) * sqrt(2/3) * phi_2(u) ] e_2
               !        = -1/2 P1(v) * [(-2) * sqrt(2/3) * phi_2(u)] e_2
               !        = -1/2 P1(v) * [4 L_1(u) L_2(u)] e_2
               WorkBasis(4,2) = (-2.0d0) * sqrt(2.0d0/3.0d0) * Phi(2,u) * (-1.0d0/2.0d0) * LegendreP(1,v)
               WorkCurlBasis(4,3) = (-2.0d0) * sqrt(2.0d0/3.0d0) * dPhi(2,u) * (-1.0d0/2.0d0) * LegendreP(1,v)

               DO l=1,FDOFs

                 SELECT CASE(l)
                 CASE(1)
                   ! (u,v) -> -sqrt(2/3) * P0 * phi_2(v) e1
                   !        = (1/2 P0) * [-2 * sqrt(2/3) * phi_2(v)] e1
                   !        = (1/2 P0) * 4 L_1(v) L_2(v) e1
                   !
                   sfun = 1.0d0
                   ! tfun = 0.0d0
                   EdgeBasis(4*EDOFs + l,1:2) = sfun * D1 * WorkBasis(I1,1:2)
                   CurlBasis(4*EDOFs + l,3) = sfun * D1 * WorkCurlBasis(I1,3)
                 CASE(2)
                   ! (u,v) -> -sqrt(2/3) * phi_2(u) * P0 e2
                   !        = (1/2 P0) * [-2 * sqrt(2/3) * phi_2(u)] e2
                   !        = (1/2 P0) * 4 L_1(u) L_2(u) e2
                   !
                   !sfun = 0.0d0
                   tfun = 1.0d0
                   EdgeBasis(4*EDOFs + l,1:2) = tfun * D2 * WorkBasis(I2,1:2)
                   CurlBasis(4*EDOFs + l,3) = tfun * D2 * WorkCurlBasis(I2,3)
                 CASE(3)
                   ! (u,v) ->  -1/2 P1(u) * [4 L_1(v) L_2(v)] e_1,  or -1/2 P1(v) * [4 L_1(u) L_2(u)] e_2
                   sfun = 1.0d0
                   tfun = 0.0d0
                   q = 2
                   ! Note that sign changes never happen
                   EdgeBasis(4*EDOFs + l,1:2) = sfun * WorkBasis(q+I1,1:2)
                   CurlBasis(4*EDOFs + l,3) = sfun * WorkCurlBasis(q+I1,3)
                 CASE(4)
                   ! (u,v) -> grad( -1/sqrt(6) * (-2 * sqrt(2/3)) * phi_2(u) * phi_2(v) )
                   !        = -1/2 P1(u) * [4 L_1(v) L_2(v)] e_1 - 1/2 P1(v) * [4 L_1(u) L_2(u)] e_2
                   !
                   sfun = 1.0d0
                   tfun = 1.0d0
                   q = 2
                   ! Note that sign changes never happen
                   EdgeBasis(4*EDOFs + l,1:2) = sfun * WorkBasis(q+I1,1:2) + tfun * WorkBasis(q+I2,1:2)
                   CurlBasis(4*EDOFs + l,3) = 0.0d0
                 END SELECT
               END DO
             ELSE
               !---------------------------------------------------------------
               ! The second-order element from the Nedelec's first family with
               ! a hierarchic basis. This element may not be optimally accurate
               ! if the physical element is not affine.
               ! First, the eight basis functions associated with the edges:
               !--------------------------------------------------------------
               i = EdgeMap(1,1)
               j = EdgeMap(1,2)
               EdgeBasis(1,1) = 0.1D1 / 0.4D1 - v / 0.4D1
               CurlBasis(1,3) = 0.1D1 / 0.4D1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(1,:) = -EdgeBasis(1,:)
                 CurlBasis(1,3) = -CurlBasis(1,3)
               END IF
               EdgeBasis(2,1) = 0.3D1 * u * (0.1D1 / 0.4D1 - v / 0.4D1)
               CurlBasis(2,3) = 0.3D1 / 0.4D1 * u

               i = EdgeMap(2,1)
               j = EdgeMap(2,2)
               EdgeBasis(3,2) = 0.1D1 / 0.4D1 + u / 0.4D1 
               CurlBasis(3,3) = 0.1D1 / 0.4D1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(3,:) = -EdgeBasis(3,:)
                 CurlBasis(3,3) = -CurlBasis(3,3)
               END IF
               EdgeBasis(4,2) = 0.3D1 * v * (0.1D1 / 0.4D1 + u / 0.4D1)
               CurlBasis(4,3) = 0.3D1 / 0.4D1 * v

               i = EdgeMap(3,1)
               j = EdgeMap(3,2)
               EdgeBasis(5,1) = -0.1D1 / 0.4D1 - v / 0.4D1
               CurlBasis(5,3) = 0.1D1 / 0.4D1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(5,:) = -EdgeBasis(5,:)
                 CurlBasis(5,3) = -CurlBasis(5,3)
               END IF
               EdgeBasis(6,1) = -0.3D1 * u * (-0.1D1 / 0.4D1 - v / 0.4D1)
               CurlBasis(6,3) = -0.3D1 / 0.4D1 * u

               i = EdgeMap(4,1)
               j = EdgeMap(4,2)
               EdgeBasis(7,2) = -0.1D1 / 0.4D1 + u / 0.4D1
               CurlBasis(7,3) = 0.1D1 / 0.4D1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(7,:) = -EdgeBasis(7,:)
                 CurlBasis(7,3) = -CurlBasis(7,3)
               END IF
               EdgeBasis(8,2) = -0.3D1 * v * (-0.1D1 / 0.4D1 + u / 0.4D1)
               CurlBasis(8,3) = -0.3D1 / 0.4D1 * v

               !--------------------------------------------------------------------
               ! Additional four basis functions associated with the element interior
               !-------------------------------------------------------------------
               SquareFaceMap(:) = (/ 1,2,3,4 /)          
               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = 0.2D1 * (0.1D1 / 0.2D1 - v / 0.2D1) * (0.1D1 / 0.2D1 + v / 0.2D1)
               WorkCurlBasis(1,3) = v
               WorkBasis(2,1) = 0.12D2 * u * (0.1D1 / 0.2D1 - v / 0.2D1) * (0.1D1 / 0.2D1 + v / 0.2D1)
               WorkCurlBasis(2,3) = 0.6D1 * u * (0.1D1 / 0.2D1 + v / 0.2D1) - &
                   0.6D1 * u * (0.1D1 / 0.2D1 - v / 0.2D1)

               WorkBasis(3,2) = 0.2D1 * (0.1D1 / 0.2D1 - u / 0.2D1) * (0.1D1 / 0.2D1 + u / 0.2D1)
               WorkCurlBasis(3,3) = -u
               WorkBasis(4,2) = 0.12D2 * v * (0.1D1 / 0.2D1 - u / 0.2D1) * (0.1D1 / 0.2D1 + u / 0.2D1)
               WorkCurlBasis(4,3) = -0.6D1 * v * (0.1D1 / 0.2D1 + u / 0.2D1) + &
                   0.6D1 * v * (0.1D1 / 0.2D1 - u / 0.2D1)

               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(9,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(9,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(10,:) = 0.5d0 * WorkBasis(2*(I1-1)+2,:)
               CurlBasis(10,:) = 0.5d0 * WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(11,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(11,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(12,:) = 0.5d0 * WorkBasis(2*(I2-1)+2,:)
               CurlBasis(12,:) = 0.5d0 * WorkCurlBasis(2*(I2-1)+2,:)
             END IF
           ELSE
             !------------------------------------------------------
             ! The Arnold-Boffi-Falk element of degree k=0 which is
             ! a member of the optimal edge element family. 
             ! First, four basis functions defined on the edges
             !-------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = ((-1.0d0 + v)*v)/4.0d0
             EdgeBasis(1,2) = 0.0d0
             CurlBasis(1,3) = (1.0d0 - 2*v)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,3) = -CurlBasis(1,3)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1.0d0 + u))/4.0d0
             CurlBasis(2,3) = (1.0d0 + 2*u)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,3) = -CurlBasis(2,3)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(3,1) = -(v*(1.0d0 + v))/4.0d0
             EdgeBasis(3,2) = 0.0d0
             CurlBasis(3,3) = (1.0d0 + 2*v)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,3) = -CurlBasis(3,3)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = -((-1 + u)*u)/4.0d0
             CurlBasis(4,3) = (1.0d0 - 2*u)/4.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,3) = -CurlBasis(4,3)
             END IF

             !--------------------------------------------------------------------
             ! Additional two basis functions associated with the element interior
             !-------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)          

             WorkBasis(1,:) = 0.0d0
             WorkBasis(2,:) = 0.0d0
             WorkCurlBasis(1,:) = 0.0d0
             WorkCurlBasis(2,:) = 0.0d0         

             WorkBasis(1,1) = (1.0d0 - v**2)/2.0d0
             WorkBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = v

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = (1.0d0 - u**2)/2.0d0
             WorkCurlBasis(2,3) = -u

             FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(5,:) = D1 * WorkBasis(I1,:)
             CurlBasis(5,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(6,:) = D2 * WorkBasis(I2,:)
             CurlBasis(6,:) = D2 * WorkCurlBasis(I2,:)
           END IF

         CASE(5)
           !--------------------------------------------------------------
           ! This branch is for handling tetrahedra
           !--------------------------------------------------------------
           EdgeMap => GetEdgeMap(5)

           IF (Create2ndKindBasis .OR. GradVersion .AND. &
               (SecondOrder .OR. ThirdOrder)) THEN

             BDOFs = 0
             IF (Create2ndKindBasis) THEN
               !
               ! This construction follows Sun, Lee, Cendes. SIAM J. Sci. Comput. 23(4):1053-1076.
               ! The first basis function associated with an edge is always the Whitney form, while 
               ! the second basis function corresponds to a gradient field.
               !
               IF (SecondOrder) THEN
                 EDOFs = 3
                 FDOFs = 3
               ELSE
                 EDOFs = 2
                 FDOFs = 0
               END IF
             ELSE
               !
               ! An alternate first-family basis of degree 2 or 3 for faster solution with iterative methods. 
               !
               IF (SecondOrder) THEN
                 EDOFs = 2
                 FDOFs = 2
               ELSE
                 ! The case of third-order basis
                 EDOFs = 3
                 FDOFs = 6
                 BDOFs = 3
               END IF
             END IF
             
             DO k=1,6
               
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)

               tvec(1:3) = Basis(i) * dLBasisdx(j,1:3)
               svec(1:3) = Basis(j) * dLBasisdx(i,1:3)

               grad_svec(1,2) = dLBasisdx(j,2) * dLBasisdx(i,1)
               grad_svec(1,3) = dLBasisdx(j,3) * dLBasisdx(i,1)
               grad_svec(2,1) = dLBasisdx(j,1) * dLBasisdx(i,2)
               grad_svec(2,3) = dLBasisdx(j,3) * dLBasisdx(i,2)
               grad_svec(3,1) = dLBasisdx(j,1) * dLBasisdx(i,3)
               grad_svec(3,2) = dLBasisdx(j,2) * dLBasisdx(i,3)

               grad_tvec(1,2) = dLBasisdx(i,2) * dLBasisdx(j,1)
               grad_tvec(1,3) = dLBasisdx(i,3) * dLBasisdx(j,1)
               grad_tvec(2,1) = dLBasisdx(i,1) * dLBasisdx(j,2)
               grad_tvec(2,3) = dLBasisdx(i,3) * dLBasisdx(j,2)
               grad_tvec(3,1) = dLBasisdx(i,1) * dLBasisdx(j,3)
               grad_tvec(3,2) = dLBasisdx(i,2) * dLBasisdx(j,3)

               WorkBasis(1,1:3) = svec(1:3)
               WorkBasis(2,1:3) = tvec(1:3)

               WorkCurlBasis(1,1) = grad_svec(3,2) - grad_svec(2,3)
               WorkCurlBasis(1,2) = grad_svec(1,3) - grad_svec(3,1)
               WorkCurlBasis(1,3) = grad_svec(2,1) - grad_svec(1,2)

               WorkCurlBasis(2,1) = grad_tvec(3,2) - grad_tvec(2,3)
               WorkCurlBasis(2,2) = grad_tvec(1,3) - grad_tvec(3,1)
               WorkCurlBasis(2,3) = grad_tvec(2,1) - grad_tvec(1,2)

               IF (Create2ndKindBasis .AND. SecondOrder .OR. &
                   GradVersion .AND. ThirdOrder) THEN
                 WorkWeight(1) = 2.0d0*Basis(i) - Basis(j)
                 WorkWeight(2) = 2.0d0*Basis(j) - Basis(i)

                 grad_weight(1,1:3) = 2.0d0*dLBasisdx(i,1:3) - dLBasisdx(j,1:3)
                 grad_weight(2,1:3) = 2.0d0*dLBasisdx(j,1:3) - dLBasisdx(i,1:3)
               END IF

               IF (GIndexes(j) < GIndexes(i)) THEN
                 I1 = 2
                 I2 = 1
               ELSE
                 I1 = 1
                 I2 = 2
               END IF

               DO l=1,EDOFs
                 SELECT CASE(l)
                 CASE(1)
                   sfun = -1.0d0
                   tfun = 1.0d0
                 CASE(2)
                   sfun = 1.0d0
                   tfun = 1.0d0
                 CASE(3)
                   sfun = -WorkWeight(I1)
                   tfun = WorkWeight(I2)
                   grad_sfun(1:3) = -grad_weight(I1,1:3)
                   grad_tfun(1:3) = grad_weight(I2,1:3)
                 CASE DEFAULT
                   CALL Fatal('ElementDescription::EdgeElementInfo','sfun/tfun not defined')
                 END SELECT

                 EdgeBasis(EDOFs*(k-1)+l,1:3) = sfun * WorkBasis(I1,1:3) + tfun * WorkBasis(I2,1:3)
                 CurlBasis(EDOFs*(k-1)+l,1:3) = sfun * WorkCurlBasis(I1,1:3) + tfun * WorkCurlBasis(I2,1:3)

                 IF (l > 2) THEN
                   CurlBasis(EDOFs*(k-1)+l,1) = CurlBasis(EDOFs*(k-1)+l,1) + &
                       grad_sfun(2)*WorkBasis(I1,3) + grad_tfun(2)*WorkBasis(I2,3) - &
                       grad_sfun(3)*WorkBasis(I1,2) - grad_tfun(3)*WorkBasis(I2,2)

                   CurlBasis(EDOFs*(k-1)+l,2) = CurlBasis(EDOFs*(k-1)+l,2) + &
                       grad_sfun(3)*WorkBasis(I1,1) + grad_tfun(3)*WorkBasis(I2,1) - &
                       grad_sfun(1)*WorkBasis(I1,3) - grad_tfun(1)*WorkBasis(I2,3)

                   CurlBasis(EDOFs*(k-1)+l,3) = CurlBasis(EDOFs*(k-1)+l,3) + &
                       grad_sfun(1)*WorkBasis(I1,2) + grad_tfun(1)*WorkBasis(I2,2) - &
                       grad_sfun(2)*WorkBasis(I1,1) - grad_tfun(2)*WorkBasis(I2,1)
                 END IF
               END DO
             END DO

             ! The basis functions associated with the faces for higher-order cases
             IF (FDOFs > 0) THEN
               DO k=1,4
                 SELECT CASE(k)
                 CASE(1)
                   TriangleFaceMap(:) = (/ 2,1,3 /)
                 CASE(2)
                   TriangleFaceMap(:) = (/ 1,2,4 /)
                 CASE(3)
                   TriangleFaceMap(:) = (/ 2,3,4 /)
                 CASE(4)
                   TriangleFaceMap(:) = (/ 3,1,4 /)
                 END SELECT

                 I1 = TriangleFaceMap(1)
                 I2 = TriangleFaceMap(2)
                 I3 = TriangleFaceMap(3)

                 WorkBasis(1,1:3) = Basis(I2) * Basis(I3) * dLBasisdx(I1,1:3)
                 WorkBasis(2,1:3) = Basis(I1) * Basis(I3) * dLBasisdx(I2,1:3)
                 WorkBasis(3,1:3) = Basis(I1) * Basis(I2) * dLBasisdx(I3,1:3)

                 ! The gradient of each row of WorkBasis:

                 grad_svec(1,2) = (dLBasisdx(I2,2) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,2)) * dLBasisdx(I1,1)
                 grad_svec(1,3) = (dLBasisdx(I2,3) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,3)) * dLBasisdx(I1,1)
                 grad_svec(2,1) = (dLBasisdx(I2,1) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,1)) * dLBasisdx(I1,2)
                 grad_svec(2,3) = (dLBasisdx(I2,3) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,3)) * dLBasisdx(I1,2)
                 grad_svec(3,1) = (dLBasisdx(I2,1) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,1)) * dLBasisdx(I1,3)
                 grad_svec(3,2) = (dLBasisdx(I2,2) * Basis(I3) + &
                     Basis(I2) * dLBasisdx(I3,2)) * dLBasisdx(I1,3)

                 grad_tvec(1,2) = (dLBasisdx(I1,2) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,2)) * dLBasisdx(I2,1)
                 grad_tvec(1,3) = (dLBasisdx(I1,3) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,3)) * dLBasisdx(I2,1)
                 grad_tvec(2,1) = (dLBasisdx(I1,1) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,1)) * dLBasisdx(I2,2)
                 grad_tvec(2,3) = (dLBasisdx(I1,3) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,3)) * dLBasisdx(I2,2)
                 grad_tvec(3,1) = (dLBasisdx(I1,1) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,1)) * dLBasisdx(I2,3)
                 grad_tvec(3,2) = (dLBasisdx(I1,2) * Basis(I3)  + &
                     Basis(I1) * dLBasisdx(I3,2)) * dLBasisdx(I2,3)

                 grad_hvec(1,2) = (dLBasisdx(I1,2) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,2)) * dLBasisdx(I3,1)
                 grad_hvec(1,3) = (dLBasisdx(I1,3) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,3)) * dLBasisdx(I3,1)
                 grad_hvec(2,1) = (dLBasisdx(I1,1) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,1)) * dLBasisdx(I3,2)
                 grad_hvec(2,3) = (dLBasisdx(I1,3) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,3)) * dLBasisdx(I3,2)
                 grad_hvec(3,1) = (dLBasisdx(I1,1) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,1)) * dLBasisdx(I3,3)
                 grad_hvec(3,2) = (dLBasisdx(I1,2) * Basis(I2)  + &
                     Basis(I1) * dLBasisdx(I2,2)) * dLBasisdx(I3,3)

                 WorkCurlBasis(1,1) = grad_svec(3,2) - grad_svec(2,3)
                 WorkCurlBasis(1,2) = grad_svec(1,3) - grad_svec(3,1)
                 WorkCurlBasis(1,3) = grad_svec(2,1) - grad_svec(1,2)

                 WorkCurlBasis(2,1) = grad_tvec(3,2) - grad_tvec(2,3)
                 WorkCurlBasis(2,2) = grad_tvec(1,3) - grad_tvec(3,1)
                 WorkCurlBasis(2,3) = grad_tvec(2,1) - grad_tvec(1,2)

                 WorkCurlBasis(3,1) = grad_hvec(3,2) - grad_hvec(2,3)
                 WorkCurlBasis(3,2) = grad_hvec(1,3) - grad_hvec(3,1)
                 WorkCurlBasis(3,3) = grad_hvec(2,1) - grad_hvec(1,2)

                 ! Create permutation:
                 FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
                 CALL TriangleFaceDofsOrdering2nd(I1,I2,I3,FaceIndices(1:3))

                 ! Create the basis:
                 DO l=1,FDOFs

                   IF (UseWForms) THEN
                     !
                     ! This allows for an alternate implementation where the basis functions
                     ! are expressed in terms of Whitney forms. For the third basis function
                     ! of gradient type the representation in terms of Whitney forms is not
                     ! yet available.
                     !
                     IF (l == 3) THEN
                       ! Revert to the values before overwriting within this loop:
                       J1 = TriangleFaceMap(1)
                       J2 = TriangleFaceMap(2)
                       J3 = TriangleFaceMap(3)

                       WorkBasis(1,1:3) = Basis(J2) * Basis(J3) * dLBasisdx(J1,1:3)
                       WorkBasis(2,1:3) = Basis(J1) * Basis(J3) * dLBasisdx(J2,1:3)
                       WorkBasis(3,1:3) = Basis(J1) * Basis(J2) * dLBasisdx(J3,1:3)
                       
                       WorkCurlBasis(1,1) = grad_svec(3,2) - grad_svec(2,3)
                       WorkCurlBasis(1,2) = grad_svec(1,3) - grad_svec(3,1)
                       WorkCurlBasis(1,3) = grad_svec(2,1) - grad_svec(1,2)

                       WorkCurlBasis(2,1) = grad_tvec(3,2) - grad_tvec(2,3)
                       WorkCurlBasis(2,2) = grad_tvec(1,3) - grad_tvec(3,1)
                       WorkCurlBasis(2,3) = grad_tvec(2,1) - grad_tvec(1,2)
                       
                       WorkCurlBasis(3,1) = grad_hvec(3,2) - grad_hvec(2,3)
                       WorkCurlBasis(3,2) = grad_hvec(1,3) - grad_hvec(3,1)
                       WorkCurlBasis(3,3) = grad_hvec(2,1) - grad_hvec(1,2)

                       CALL TriangleFaceDofsOrdering2nd(I1,I2,I3,FaceIndices(1:3))
                     ELSE
                       CALL WeightedWhitneyForms(WorkBasis(1:3,:), WorkCurlBasis(1:3,:), k, u, v, w)
                       CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices,A0,B0,C0)
                     END IF
                   END IF
                   
                   SELECT CASE(l)
                   CASE(1)
                     !
                     ! This creates the basis function L_C W_{AB} - 2 L_B W_{AC} 
                     !
                     IF (UseWForms) THEN
                       sfun = 1.0d0
                       tfun = -2.0d0
                       hfun = 0.0d0
                     ELSE
                       sfun = 1.0d0
                       tfun = 1.0d0
                       hfun = -2.0d0
                     END IF
                   CASE(2)
                     !
                     ! This creates the basis function -L_C W_{AB}
                     !
                     IF (UseWForms) THEN
                       sfun = -1.0d0
                       tfun = 0.0d0
                       hfun = 0.0d0
                     ELSE
                       sfun = 1.0d0
                       tfun = -1.0d0
                       hfun = 0.0d0
                     END IF
                   CASE(3)
                     ! This corresponds to the second-order gradient:
                     sfun = 1.0d0
                     tfun = 1.0d0
                     hfun = 1.0d0
                   CASE(4)
                     !
                     ! This creates the basis function (L_A - L_B) L_B W_{AC} + (L_C - L_A) L_C W_{AB} 
                     ! 
                     IF (UseWForms) THEN
                       sfun = Basis(C0) - Basis(A0)
                       tfun = Basis(A0) - Basis(B0)
                       hfun = 0.0d0
                       grad_sfun(1:3) = dLBasisdx(C0,1:3) - dLBasisdx(A0,1:3)
                       grad_tfun(1:3) = dLBasisdx(A0,1:3) - dLBasisdx(B0,1:3)
                       grad_hfun(1:3) = 0.0d0
                     ELSE
                       grad_sfun(1:3) = dLBasisdx(I2,1:3) - dLBasisdx(I3,1:3)
                       grad_tfun(1:3) = dLBasisdx(I3,1:3) - dLBasisdx(I1,1:3)
                       grad_hfun(1:3) = dLBasisdx(I1,1:3) - dLBasisdx(I2,1:3)
                     END IF
                   CASE(5)
                     !
                     ! This creates the basis function (-80 L_A - 393 L_B + 212 L_C) L_C W_{AB} + (-313 L_A + 313 L_B) L_B W_{AC} 
                     !
                     IF (UseWForms) THEN
                       sfun = -393.0d0 * Basis(B0) - 80.0d0 * Basis(A0) + 212.0d0 * Basis(C0)
                       tfun = -313.0d0 * Basis(A0) + 313.0d0 * Basis(B0)
                       hfun = 0.0d0
                       
                       grad_sfun(1:3) = -393.0d0 * dLBasisdx(B0,1:3) - 80.0d0 * dLBasisdx(A0,1:3) + 212.0d0 * dLBasisdx(C0,1:3)
                       grad_tfun(1:3) = -313.0d0 * dLBasisdx(A0,1:3) + 313.0d0 * dLBasisdx(B0,1:3)
                       grad_hfun(1:3) = 0.0d0
                     ELSE
                       ! Note that hfun contains a correction to the basis proposed by Sun, Lee, Cendes
                       ! so that the weight functions form a partition of the unity.
                       sfun = 393.0d0 * Basis(I1) + 80.0d0 * Basis(I2) - 212.0d0 * Basis(I3)
                       tfun = -393.0d0 * Basis(I2) - 80.0d0 * Basis(I1) + 212.0d0 * Basis(I3)
                       hfun = -313.0d0 * Basis(I1) + 313.0d0 * Basis(I2)

                       grad_sfun(1:3) = 393.0d0 * dLBasisdx(I1,1:3) + 80.0d0 * dLBasisdx(I2,1:3) - 212.0d0 * dLBasisdx(I3,1:3)
                       grad_tfun(1:3) = -393.0d0 * dLBasisdx(I2,1:3) - 80.0d0 * dLBasisdx(I1,1:3) + 212.0d0 * dLBasisdx(I3,1:3)
                       grad_hfun(1:3) = -313.0d0 * dLBasisdx(I1,1:3) + 313.0d0 * dLBasisdx(I2,1:3)
                     END IF
                   CASE(6)
                     !
                     ! This creates the basis function (168 L_A - 131 L_B - 124 L_C) L_C W_{AB} + (-37 L_A - 37 L_B + 248 L_C) L_B W_{AC} 
                     !
                     IF (UseWForms) THEN
                       sfun = -131.0d0 * Basis(B0) + 168.0d0 * Basis(A0) - 124.0d0 * Basis(C0)
                       tfun = -37.0d0 * Basis(A0) - 37.0d0 * Basis(B0) + 248.0d0 * Basis(C0)
                       hfun = 0.0d0

                       grad_sfun(1:3) = -131.0d0 * dLBasisdx(B0,1:3) + 168.0d0 * dLBasisdx(A0,1:3) - 124.0d0 * dLBasisdx(C0,1:3)
                       grad_tfun(1:3) = -37.0d0 * dLBasisdx(A0,1:3) - 37.0d0 * dLBasisdx(B0,1:3) + 248.0d0 * dLBasisdx(C0,1:3)
                       grad_hfun(1:3) = 0.0d0
                     ELSE
                       ! Note that hfun contains a correction to the basis proposed by Sun, Lee, Cendes
                       ! so that the weight functions form a partition of the unity
                       sfun = -131.0d0 * Basis(I1) + 168.0d0 * Basis(I2) - 124.0d0 * Basis(I3)
                       tfun = -131.0d0 * Basis(I2) + 168.0d0 * Basis(I1) - 124.0d0 * Basis(I3)
                       hfun = -37.0d0 * Basis(I1) - 37.0d0 * Basis(I2) + 248.0d0 * Basis(I3)

                       grad_sfun(1:3) = -131.0d0 * dLBasisdx(I1,1:3) + 168.0d0 * dLBasisdx(I2,1:3) - 124.0d0 * dLBasisdx(I3,1:3)
                       grad_tfun(1:3) = -131.0d0 * dLBasisdx(I2,1:3) + 168.0d0 * dLBasisdx(I1,1:3) - 124.0d0 * dLBasisdx(I3,1:3)
                       grad_hfun(1:3) = -37.0d0 * dLBasisdx(I1,1:3) - 37.0d0 * dLBasisdx(I2,1:3) + 248.0d0 * dLBasisdx(I3,1:3)
                     END IF
                     
                   END SELECT

                   IF (UseWForms .AND. l /= 3) THEN
                     EdgeBasis(6*EDOFs + FDOFs*(k-1)+l,1:3) = sfun * D1 * WorkBasis(I1,1:3) + tfun * D2 * WorkBasis(I2,1:3)
                     CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1:3) = sfun * D1 * WorkCurlBasis(I1,1:3) + tfun * D2 * WorkCurlBasis(I2,1:3)
                   ELSE
                     EdgeBasis(6*EDOFs + FDOFs*(k-1)+l,1:3) = sfun * WorkBasis(I1,1:3) + tfun * WorkBasis(I2,1:3) + &
                         hfun * WorkBasis(I3,1:3)
                     CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1:3) = sfun * WorkCurlBasis(I1,1:3) + tfun * WorkCurlBasis(I2,1:3) + &
                         hfun * WorkCurlBasis(I3,1:3)
                   END IF
                   
                   IF (l > 3) THEN
                     IF (UseWForms) THEN
                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1) + &
                           grad_sfun(2)*D1*WorkBasis(I1,3) + grad_tfun(2)*D2*WorkBasis(I2,3) - &
                           grad_sfun(3)*D1*WorkBasis(I1,2) - grad_tfun(3)*D2*WorkBasis(I2,2)

                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,2) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,2) + &
                           grad_sfun(3)*D1*WorkBasis(I1,1) + grad_tfun(3)*D2*WorkBasis(I2,1) - &
                           grad_sfun(1)*D1*WorkBasis(I1,3) - grad_tfun(1)*D2*WorkBasis(I2,3)

                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,3) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,3) + &
                           grad_sfun(1)*D1*WorkBasis(I1,2) + grad_tfun(1)*D2*WorkBasis(I2,2) - &
                           grad_sfun(2)*D1*WorkBasis(I1,1) - grad_tfun(2)*D2*WorkBasis(I2,1)
                       
                     ELSE
                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,1) + &
                           grad_sfun(2)*WorkBasis(I1,3) + grad_tfun(2)*WorkBasis(I2,3) + grad_hfun(2)*WorkBasis(I3,3) - &
                           grad_sfun(3)*WorkBasis(I1,2) - grad_tfun(3)*WorkBasis(I2,2) - grad_hfun(3)*WorkBasis(I3,2)

                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,2) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,2) + &
                           grad_sfun(3)*WorkBasis(I1,1) + grad_tfun(3)*WorkBasis(I2,1) + grad_hfun(3)*WorkBasis(I3,1) - &
                           grad_sfun(1)*WorkBasis(I1,3) - grad_tfun(1)*WorkBasis(I2,3) - grad_hfun(1)*WorkBasis(I3,3)

                       CurlBasis(6*EDOFs + FDOFs*(k-1)+l,3) = CurlBasis(6*EDOFs + FDOFs*(k-1)+l,3) + &
                           grad_sfun(1)*WorkBasis(I1,2) + grad_tfun(1)*WorkBasis(I2,2) + grad_hfun(1)*WorkBasis(I3,2) - &
                           grad_sfun(2)*WorkBasis(I1,1) - grad_tfun(2)*WorkBasis(I2,1) - grad_hfun(2)*WorkBasis(I3,1)
                     END IF
                   END IF


                 END DO
               END DO

             END IF

             IF (BDOFs > 0) THEN
               I1 = 1
               I2 = 2
               I3 = 3
               I4 = 4
               
               WorkBasis(1,1:3) = Basis(I2) * Basis(I3) * Basis(I4) * dLBasisdx(I1,1:3)
               WorkBasis(2,1:3) = Basis(I1) * Basis(I3) * Basis(I4) * dLBasisdx(I2,1:3)
               WorkBasis(3,1:3) = Basis(I1) * Basis(I2) * Basis(I4) * dLBasisdx(I3,1:3)
               WorkBasis(4,1:3) = Basis(I1) * Basis(I2) * Basis(I3) * dLBasisdx(I4,1:3)

               ! The gradient of each row of WorkBasis (TO DO: Restructure as a loop)

               grad_svec(1,2) = (dLBasisdx(I2,2) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,2) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,2)) * dLBasisdx(I1,1)
               grad_svec(1,3) = (dLBasisdx(I2,3) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,3) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,3)) * dLBasisdx(I1,1)
               grad_svec(2,1) = (dLBasisdx(I2,1) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,1) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,1)) * dLBasisdx(I1,2)
               grad_svec(2,3) = (dLBasisdx(I2,3) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,3) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,3)) * dLBasisdx(I1,2)
               grad_svec(3,1) = (dLBasisdx(I2,1) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,1) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,1)) * dLBasisdx(I1,3)
               grad_svec(3,2) = (dLBasisdx(I2,2) * Basis(I3) * Basis(I4) + &
                   Basis(I2) * dLBasisdx(I3,2) * Basis(I4) + &
                   Basis(I2) * Basis(I3) * dLBasisdx(I4,2)) * dLBasisdx(I1,3)
               
               grad_tvec(1,2) = (dLBasisdx(I1,2) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,2) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,2)) * dLBasisdx(I2,1)
               grad_tvec(1,3) = (dLBasisdx(I1,3) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,3) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,3)) * dLBasisdx(I2,1)
               grad_tvec(2,1) = (dLBasisdx(I1,1) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,1) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,1)) * dLBasisdx(I2,2)
               grad_tvec(2,3) = (dLBasisdx(I1,3) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,3) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,3)) * dLBasisdx(I2,2)
               grad_tvec(3,1) = (dLBasisdx(I1,1) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,1) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,1)) * dLBasisdx(I2,3)
               grad_tvec(3,2) = (dLBasisdx(I1,2) * Basis(I3) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I3,2) * Basis(I4) + &
                   Basis(I1) * Basis(I3) * dLBasisdx(I4,2)) * dLBasisdx(I2,3)
               
               grad_hvec(1,2) = (dLBasisdx(I1,2) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,2) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,2)) * dLBasisdx(I3,1)
               grad_hvec(1,3) = (dLBasisdx(I1,3) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,3) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,3)) * dLBasisdx(I3,1)
               grad_hvec(2,1) = (dLBasisdx(I1,1) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,1) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,1)) * dLBasisdx(I3,2)
               grad_hvec(2,3) = (dLBasisdx(I1,3) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,3) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,3)) * dLBasisdx(I3,2)
               grad_hvec(3,1) = (dLBasisdx(I1,1) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,1) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,1)) * dLBasisdx(I3,3)
               grad_hvec(3,2) = (dLBasisdx(I1,2) * Basis(I2) * Basis(I4) + &
                   Basis(I1) * dLBasisdx(I2,2) * Basis(I4) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I4,2)) * dLBasisdx(I3,3)
               
               grad_gvec(1,2) = (dLBasisdx(I1,2) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,2) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,2)) * dLBasisdx(I4,1)
               grad_gvec(1,3) = (dLBasisdx(I1,3) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,3) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,3)) * dLBasisdx(I4,1)
               grad_gvec(2,1) = (dLBasisdx(I1,1) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,1) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,1)) * dLBasisdx(I4,2)
               grad_gvec(2,3) = (dLBasisdx(I1,3) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,3) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,3)) * dLBasisdx(I4,2)
               grad_gvec(3,1) = (dLBasisdx(I1,1) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,1) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,1)) * dLBasisdx(I4,3)
               grad_gvec(3,2) = (dLBasisdx(I1,2) * Basis(I2) * Basis(I3) + &
                   Basis(I1) * dLBasisdx(I2,2) * Basis(I3) + &
                   Basis(I1) * Basis(I2) * dLBasisdx(I3,2)) * dLBasisdx(I4,3)

               WorkCurlBasis(I1,1) = grad_svec(3,2) - grad_svec(2,3)
               WorkCurlBasis(I1,2) = grad_svec(1,3) - grad_svec(3,1)
               WorkCurlBasis(I1,3) = grad_svec(2,1) - grad_svec(1,2)

               WorkCurlBasis(I2,1) = grad_tvec(3,2) - grad_tvec(2,3)
               WorkCurlBasis(I2,2) = grad_tvec(1,3) - grad_tvec(3,1)
               WorkCurlBasis(I2,3) = grad_tvec(2,1) - grad_tvec(1,2)

               WorkCurlBasis(I3,1) = grad_hvec(3,2) - grad_hvec(2,3)
               WorkCurlBasis(I3,2) = grad_hvec(1,3) - grad_hvec(3,1)
               WorkCurlBasis(I3,3) = grad_hvec(2,1) - grad_hvec(1,2)

               WorkCurlBasis(I4,1) = grad_gvec(3,2) - grad_gvec(2,3)
               WorkCurlBasis(I4,2) = grad_gvec(1,3) - grad_gvec(3,1)
               WorkCurlBasis(I4,3) = grad_gvec(2,1) - grad_gvec(1,2)
               

               DO l=1,BDOFs

                 SELECT CASE(l)
                 CASE(1)
                   sfun = 1.0d0
                   tfun = 1.0d0
                   hfun = -1.0d0
                   gfun = -1.0d0
                 CASE(2)
                   sfun = 0.0d0
                   tfun = 0.0d0
                   hfun = 1.0d0
                   gfun = -1.0d0
                 CASE(3)
                   sfun = 1.0d0
                   tfun = -1.0d0
                   hfun = 0.0d0
                   gfun = 0.0d0
                 CASE DEFAULT
                   CALL Fatal('ElementDescription::EdgeElementInfo','Bubble count exceeds the current ability')                   
                 END SELECT

                 EdgeBasis(6*EDOFs + 4*FDOFs + l,1:3) = sfun * WorkBasis(I1,1:3) + tfun * WorkBasis(I2,1:3) + &
                       hfun * WorkBasis(I3,1:3) + gfun * WorkBasis(I4,1:3)
                 CurlBasis(6*EDOFs + 4*FDOFs + l,1:3) = sfun * WorkCurlBasis(I1,1:3) + tfun * WorkCurlBasis(I2,1:3) + &
                       hfun * WorkCurlBasis(I3,1:3) + gfun * WorkCurlBasis(I4,1:3)
               END DO

             END IF

           ELSE
             !-------------------------------------------------------------
             ! The optimal/Nedelec basis functions of the first kind. We employ
             ! a hierarchic basis, so the lowest-order basis functions are
             ! also utilized in the construction of the second-order basis. 
             ! The first the edge ...
             !-------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = (6.0d0 - 2.0d0*Sqrt(3.0d0)*v - Sqrt(6.0d0)*w)/24.0d0
             EdgeBasis(1,2) = u/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(1,3) = u/(4.0d0*Sqrt(6.0d0))            
             CurlBasis(1,1) = 0.0d0
             CurlBasis(1,2) = -1.0d0/(2.0d0*Sqrt(6.0d0))
             CurlBasis(1,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF
             IF (SecondOrder) THEN
               EdgeBasis(2,1) = -(u*(-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/4.0d0
               EdgeBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
               EdgeBasis(2,3) = (Sqrt(1.5d0)*u**2)/2.0d0
               CurlBasis(2,1) = 0.0d0
               CurlBasis(2,2) = (-3.0d0*Sqrt(1.5d0)*u)/2.0d0
               CurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0                   
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the second edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 3
               EdgeBasis(4,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*(4.0d0*v - Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(4,2) = -((1.0d0 + u - Sqrt(3.0d0)*v)*&
                   (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(4,3) = -((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*&
                   (-1.0d0 - u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,1) = (-9.0d0*(1.0d0 + u - Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,2) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(4,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0
             ELSE
               k = 2
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(k,1) = (-4.0d0*v + Sqrt(2.0d0)*w)/(16.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w)/48.0d0
             EdgeBasis(k,3) = -(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)/(24.0d0*Sqrt(2.0d0))
             CurlBasis(k,1) = 1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 1.0d0/(4.0d0*Sqrt(6.0d0))
             CurlBasis(k,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the third edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 5
               EdgeBasis(6,1) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)*&
                   (4.0d0*v - Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(6,2) = -((-1.0d0 + u + Sqrt(3.0d0)*v)*&
                   (-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w))/16.0d0
               EdgeBasis(6,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)*&
                   (-1.0d0 + u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,1) = (9.0d0*(-1.0d0 + u + Sqrt(3.0d0)*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,2) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(6,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
             ELSE
               k = 3
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(k,1) = (-4.0d0*v + Sqrt(2.0d0)*w)/(16.0d0*Sqrt(3.0d0))
             EdgeBasis(k,2) = (-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w)/48.0d0
             EdgeBasis(k,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - 3.0d0*Sqrt(2.0d0)*v)/48.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 1.0d0/(4.0d0*Sqrt(6.0d0))
             CurlBasis(k,3) = 1.0d0/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the fourth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 7
               EdgeBasis(8,1) = (3.0d0*w*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 4.0d0*w))/16.0d0
               EdgeBasis(8,2) = (w*(-3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + &
                   4.0d0*Sqrt(3.0d0)*w))/16.0d0
               EdgeBasis(8,3) = -((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)*&
                   (-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v + 2.0d0*Sqrt(6.0d0)*w))/(8.0d0*Sqrt(2.0d0))
               CurlBasis(8,1) = (-3.0d0*(-3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u + &
                   Sqrt(6.0d0)*v + 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               CurlBasis(8,2) = (9.0d0*(-Sqrt(6.0d0) + Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 4.0d0*w))/16.0d0
               CurlBasis(8,3) = 0.0d0
             ELSE
               k = 4
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(k,1) = (Sqrt(1.5d0)*w)/8.0d0
             EdgeBasis(k,2) = w/(8.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = (Sqrt(6.0d0) - Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)/16.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = Sqrt(1.5d0)/4.0d0
             CurlBasis(k,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the fifth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 9
               EdgeBasis(10,1) = (3.0d0*(Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 4.0d0*w)*w)/16.0d0
               EdgeBasis(10,2) = (w*(-3.0d0*Sqrt(2.0d0) - 3.0d0*Sqrt(2.0d0)*u + &
                   Sqrt(6.0d0)*v + 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               EdgeBasis(10,3) = ((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
                   (-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v + 2.0d0*Sqrt(6.0d0)*w))/16.0d0
               CurlBasis(10,1) = (3.0d0*(3.0d0*Sqrt(2.0d0) + 3.0d0*Sqrt(2.0d0)*u - &
                   Sqrt(6.0d0)*v - 4.0d0*Sqrt(3.0d0)*w))/16.0d0
               CurlBasis(10,2) = (9.0d0*(Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 4.0d0*w))/16.0d0
               CurlBasis(10,3) = 0.0d0
             ELSE
               k = 5
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             EdgeBasis(k,1) = -(Sqrt(1.5d0)*w)/8.0d0
             EdgeBasis(k,2) = w/(8.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = (Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)/16.0d0
             CurlBasis(k,1) = -1.0d0/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = -Sqrt(1.5d0)/4.0d0
             CurlBasis(k,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             !-------------------------------------------------
             ! Basis functions associated with the sixth edge:
             !-------------------------------------------------
             IF (SecondOrder) THEN
               k = 11
               EdgeBasis(12,1) = 0.0d0
               EdgeBasis(12,2) = (Sqrt(3.0d0)*(Sqrt(2.0d0)*v - 2.0d0*w)*w)/4.0d0
               EdgeBasis(12,3) = (Sqrt(1.5d0)*v*(-v + Sqrt(2.0d0)*w))/2.0d0
               CurlBasis(12,1) = (-3.0d0*(Sqrt(6.0d0)*v - 2.0d0*Sqrt(3.0d0)*w))/4.0d0
               CurlBasis(12,2) = 0.0d0
               CurlBasis(12,3) = 0.0d0
             ELSE
               k = 6
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             EdgeBasis(k,1) = 0.0d0
             EdgeBasis(k,2) = -w/(4.0d0*Sqrt(2.0d0))
             EdgeBasis(k,3) = v/(4.0d0*Sqrt(2.0d0))
             CurlBasis(k,1) = 1.0d0/(2.0d0*Sqrt(2.0d0))
             CurlBasis(k,2) = 0.0d0
             CurlBasis(k,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
             END IF

             ! -------------------------------------------------------------
             ! Finally scale the lowest-order basis functions so that 
             ! (b,t) = 1 when the integration is done over the element edge.
             ! -------------------------------------------------------------
             IF (SecondOrder) THEN
               DO k=1,6
                 EdgeBasis(2*(k-1)+1,:) = 2.0d0 * EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = 2.0d0 * CurlBasis(2*(k-1)+1,:)
               END DO
             ELSE
               DO k=1,6
                 EdgeBasis(k,:) = 2.0d0 * EdgeBasis(k,:)
                 CurlBasis(k,:) = 2.0d0 * CurlBasis(k,:)
               END DO
             END IF

             IF (SecondOrder) THEN

               DO k=1,4
                 !
                 ! Two additional basis functions on each face
                 !
                 SELECT CASE(k)
                 CASE(1)
                   TriangleFaceMap(:) = (/ 2,1,3 /)
                 CASE(2)
                   TriangleFaceMap(:) = (/ 1,2,4 /)
                 CASE(3)
                   TriangleFaceMap(:) = (/ 2,3,4 /)
                 CASE(4)
                   TriangleFaceMap(:) = (/ 3,1,4 /)
                 END SELECT

                 FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
                 CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

                 CALL WeightedWhitneyForms(WorkBasis(1:3,1:3), WorkCurlBasis(1:3,1:3), k, u, v, w)
                 
                 IF (RedefineFaceBasis) THEN
                   EdgeBasis(12+2*(k-1)+1,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
                   CurlBasis(12+2*(k-1)+1,:) = 0.5d0 * D1 * WorkCurlBasis(I1,:) + 0.5d0 * D2 * WorkCurlBasis(I2,:)
                   EdgeBasis(12+2*k,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
                   CurlBasis(12+2*k,:) = 0.5d0 * D2 * WorkCurlBasis(I2,:) - 0.5d0 * D1 * WorkCurlBasis(I1,:)
                 ELSE
                   EdgeBasis(12+2*(k-1)+1,:) = D1 * WorkBasis(I1,:)
                   CurlBasis(12+2*(k-1)+1,:) = D1 * WorkCurlBasis(I1,:)
                   EdgeBasis(12+2*k,:) = D2 * WorkBasis(I2,:)
                   CurlBasis(12+2*k,:) = D2 * WorkCurlBasis(I2,:)  
                 END IF
               END DO
                 
               ! Finally, scale to reduce ill-conditioning:
               IF (ScaleFaceBasis) THEN
                 EdgeBasis(13:20:2,:) = sqrt(fs1) * EdgeBasis(13:20:2,:)
                 CurlBasis(13:20:2,:) = sqrt(fs1) * CurlBasis(13:20:2,:)
                 EdgeBasis(14:20:2,:) = sqrt(fs2) * EdgeBasis(14:20:2,:)
                 CurlBasis(14:20:2,:) = sqrt(fs2) * CurlBasis(14:20:2,:)                 
               END IF
             END IF
           END IF

         CASE(6)
           !--------------------------------------------------------------
           ! This branch is for handling pyramidic elements
           !--------------------------------------------------------------         
           EdgeMap => GetEdgeMap(6)

           IF (SecondOrder) THEN
             EdgeSign = 1.0d0

             LBasis(1) = 0.1D1 / 0.4D1 - u / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 + &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(2) = 0.1D1 / 0.4D1 + u / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 - &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(3) = 0.1D1 / 0.4D1 + u / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 + &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(4) = 0.1D1 / 0.4D1 - u / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1 - &
                 u * v / ( (0.1D1 - w * sqrt(0.2D1) / 0.2D1) * 0.4D1 )
             LBasis(5) = w * sqrt(0.2D1) / 0.2D1

             Beta(1) = 0.1D1 / 0.2D1 - u / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(2) = 0.1D1 / 0.2D1 - v / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(3) = 0.1D1 / 0.2D1 + u / 0.2D1 - w * sqrt(0.2D1) / 0.4D1
             Beta(4) = 0.1D1 / 0.2D1 + v / 0.2D1 - w * sqrt(0.2D1) / 0.4D1

             ! Edge 12:
             !--------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = 0.1D1 / 0.4D1 - v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = sqrt(0.2D1) * u * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ((w * sqrt(0.2D1) - 0.2D1) * 0.8D1)
             CurlBasis(1,1) = sqrt(0.2D1) * u / ((w * sqrt(0.2D1) - 0.2D1) * 0.4D1)
             CurlBasis(1,2) = -sqrt(0.2D1) / 0.8D1 - sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(1,3) = 0.1D1 / 0.4D1
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
               EdgeSign(1) = -1.0d0
             END IF

             EdgeBasis(2,1:3) = 3.0d0 * u * EdgeBasis(1,1:3)
             CurlBasis(2,1) = 0.3D1 / 0.4D1 * u ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(2,2) = -0.3D1 / 0.8D1 * u * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) + &
                 4.0D0 * v - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(2,3) = 0.3D1 / 0.4D1 * u

             ! Edge 23:
             !--------------------------------------------------------------
             k = 3 ! k=2 for first-order
             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(k,1) = 0.0d0
             EdgeBasis(k,2) = 0.1D1 / 0.4D1 + u / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,3) = sqrt(0.2D1) * v * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(k,1) = sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 ) + sqrt(0.2D1) /  0.8D1
             CurlBasis(k,2) = sqrt(0.2D1) * v / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,3) = 0.1D1 / 0.4D1
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * v * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * v * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) - & 
                 4.0D0 * u - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = 0.3D1 / 0.4D1 * v ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = 0.3D1 / 0.4D1 * v

             ! Edge 43:
             !--------------------------------------------------------------
             k = 5 ! k=3 for first-order
             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(k,1) = 0.1D1 / 0.4D1 + v / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,2) = 0.0d0
             EdgeBasis(k,3) = sqrt(0.2D1) * u * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )

             CurlBasis(k,1) = -sqrt(0.2D1) * u / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,2) = -sqrt(0.2D1) / 0.8D1 - sqrt(0.2D1) * (w * sqrt(0.2D1) - &
                 2.0D0 * v - 0.2D1) / ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             CurlBasis(k,3) = -0.1D1 / 0.4D1
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * u * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = -0.3D1 / 0.4D1 * u ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * u * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) - &
                 4.0D0 * v - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = -0.3D1 / 0.4D1 * u


             ! Edge 14:
             !--------------------------------------------------------------
             k = 7 ! k=4 for first-order
             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(k,1) = 0.0d0 
             EdgeBasis(k,2) = 0.1D1 / 0.4D1 - u / 0.4D1 - w * sqrt(0.2D1) / 0.8D1
             EdgeBasis(k,3) = sqrt(0.2D1) * v * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / & 
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )

             CurlBasis(k,1) = sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / ( (w * &
                 sqrt(0.2D1) - 0.2D1) * 0.8D1 ) + sqrt(0.2D1) / 0.8D1
             CurlBasis(k,2) = -sqrt(0.2D1) * v / ( (w * sqrt(0.2D1) - 0.2D1) * 0.4D1 )
             CurlBasis(k,3) = -0.1D1 / 0.4D1
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * v * EdgeBasis(k,1:3)
             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * v * sqrt(0.2D1) * (0.3D1 * w * sqrt(0.2D1) + &
                 4.0D0 * u - 0.6D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,2) = -0.3D1 / 0.4D1 * v ** 2 * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1)
             CurlBasis(k+1,3) = -0.3D1 / 0.4D1 * v


             ! Edge 15:
             !--------------------------------------------------------------
             k = 9 ! k=5 for first-order             
             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             EdgeBasis(k,1) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,2) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,3) = -sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w - &
                 0.2D1 * sqrt(0.2D1) * u * w - &
                 0.2D1 * sqrt(0.2D1) * v * w + u * w ** 2 + v * w ** 2 + 0.2D1 * w * sqrt(0.2D1) - &
                 0.2D1 * u * v - w ** 2 + 0.2D1 * u + 0.2D1 * v - 0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2 

             CurlBasis(k,1) = (-sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * &
                 u * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = -(-sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * &
                 v * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(1)+LBasis(3) )

             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * (-0.9D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w - &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + &
                 0.24D2 * u * w + 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * (-0.3D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) + &
                 0.4D1 * sqrt(0.2D1) * v ** 2 + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u* v * w - &
                 0.4D1 * v ** 2 * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + &
                 0.12D2 * u * w + 0.24D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,3) = 0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u - v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 25:
             !--------------------------------------------------------------
             k = 11 ! k=6 for first-order  
             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             EdgeBasis(k,1) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,2) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             EdgeBasis(k,3) = sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w - 0.2D1 * &
                 sqrt(0.2D1) * u * w + 0.2D1 * sqrt(0.2D1) * v * w + u * w ** 2 - v * w ** 2 - &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v + w ** 2 + 0.2D1 * u - 0.2D1 * v + 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = (-sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * & 
                 v * w - 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(2)+LBasis(4) )

             CurlBasis(k+1,1) = 0.3D1 / 0.8D1 * (0.9D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u** 2 * w + &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - &
                 0.6D1 * v * sqrt(0.2D1) - 0.24D2 * u * w + 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = -0.3D1 / 0.8D1 * (-0.3D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - &
                 0.4D1 * sqrt(0.2D1) * v ** 2 - 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.4D1 * v ** 2 * w + 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) + 0.12D2 * u * w - 0.24D2 * v * w - 0.2D1 * sqrt(0.2D1) + &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)** 2
             CurlBasis(k+1,3) = 0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u + v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 35:
             !--------------------------------------------------------------
             k = 13 ! k=7 for first-order  
             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             EdgeBasis(k,1) = -w * sqrt(0.2D1)/ 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) 
             EdgeBasis(k,2) = -w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / & 
                 (w * sqrt(0.2D1) - 0.2D1)
             EdgeBasis(k,3) = -sqrt(0.2D1)/ 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w + 0.2D1 * &
                 sqrt(0.2D1) * u * w + 0.2D1 * sqrt(0.2D1) * v * w - u * w ** 2 - v * w ** 2 + &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v - w ** 2 - 0.2D1 * u - 0.2D1 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,2) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * &
                 v * w + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) ** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(3)+LBasis(1) )

             CurlBasis(k+1,1) = -0.3D1 / 0.8D1 * (0.9D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w - &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) - 0.24D2 * u * w - 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,2) = 0.3D1 / 0.8D1 * (0.3D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) + &
                 0.4D1 * sqrt(0.2D1) * v ** 2 + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u *v * w - &
                 0.4D1 * v ** 2 * w - 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) - &
                 0.12D2 * u * w - 0.24D2 * v * w + 0.2D1 * sqrt(0.2D1) - 0.16D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             CurlBasis(k+1,3) = -0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u - v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Edge 45:
             !--------------------------------------------------------------
             k = 15 ! k=8 for first-order  
             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             EdgeBasis(k,1) = w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) 
             EdgeBasis(k,2) = -w * sqrt(0.2D1) / 0.8D1 * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1)
             EdgeBasis(k,3) = sqrt(0.2D1) / 0.4D1 * (0.2D1 * sqrt(0.2D1) * u * v * w + &
                 0.2D1 * sqrt(0.2D1) * u * w - 0.2D1 * sqrt(0.2D1) * v * w - u * w ** 2 + v * w ** 2 - &
                 0.2D1 * w * sqrt(0.2D1) - 0.2D1 * u * v + w ** 2 - 0.2D1 * u + 0.2D1 * v + 0.2D1) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2 
             CurlBasis(k,1) = -(-sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.2D1 * u * w - &
                 0.2D1 * sqrt(0.2D1) + 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1)** 2 * 0.2D1 )
             CurlBasis(k,2) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.2D1 * v * w + &
                 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1)** 2 * 0.2D1 )
             CurlBasis(k,3) = 0.0d0 
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(k,:) = -EdgeBasis(k,:)
               CurlBasis(k,:) = -CurlBasis(k,:)
               EdgeSign(k) = -1.0d0
             END IF

             EdgeBasis(k+1,1:3) = 3.0d0 * EdgeSign(k) * EdgeBasis(k,1:3) * ( LBasis(5)-LBasis(4)+LBasis(2) )

             CurlBasis(k+1,1) = -0.3D1 / 0.8D1 * (-0.9D1 * sqrt(0.2D1) * u * w ** 2 + &
                 0.3D1 * sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u** 2 * w + &
                 0.8D1 * u * v * w - 0.6D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + &
                 0.6D1 * v * sqrt(0.2D1) + 0.24D2 * u * w - 0.12D2 * v * w + 0.2D1 * sqrt(0.2D1) - &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             CurlBasis(k+1,2) = 0.3D1 / 0.8D1 * (0.3D1 * sqrt(0.2D1) * u * w ** 2 - &
                 0.9D1 * sqrt(0.2D1) * v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - &
                 0.4D1 * sqrt(0.2D1) * v ** 2 - 0.13D2 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u *v * w + &
                 0.4D1 * v ** 2 * w + 0.6D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - &
                 0.6D1 * v * sqrt(0.2D1) - 0.12D2 * u * w + 0.24D2 * v * w - 0.2D1 * sqrt(0.2D1) + &
                 0.16D2 * w) / (w * sqrt(0.2D1) - 0.2D1)**2
             CurlBasis(k+1,3) = -0.3D1 / 0.8D1 * w * sqrt(0.2D1) * (u + v) / (w * sqrt(0.2D1) - 0.2D1)


             ! Square face:
             ! ------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)

             WorkBasis(1,1:3) = 2.0d0 * ( EdgeSign(1) * EdgeBasis(1,1:3) * Beta(4) + &
                 EdgeSign(5) * EdgeBasis(5,1:3) * Beta(2) ) / (1.0d0 - LBasis(5))
             WorkCurlBasis(1,1) = -0.2D1 * u * v * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(1,2) = -(sqrt(0.2D1) * w ** 2 + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / & 
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(1,3) = -0.2D1 * v / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(2,1:3) = 3.0d0 * WorkBasis(1,1:3) * u
             WorkCurlBasis(2,1) = -0.6D1 * u ** 2 * sqrt(0.2D1) * v / (w * sqrt(0.2D1) - 0.2D1)** 2
             WorkCurlBasis(2,2) = 0.3D1 / 0.2D1 * u * (0.2D1 * sqrt(0.2D1) * v ** 2 - &
                 0.3D1 * sqrt(0.2D1) * w ** 2 - 0.6D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(2,3) = -0.6D1 * u * v / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(3,1:3) = 2.0d0 * ( EdgeSign(3) * EdgeBasis(3,1:3) * Beta(1) + &
                 EdgeSign(7) * EdgeBasis(7,1:3) * Beta(3) ) / (1.0d0 - LBasis(5))
             WorkCurlBasis(3,1) = (sqrt(0.2D1) * w ** 2 + 0.2D1 * sqrt(0.2D1) - 0.4D1 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(3,2) = 0.2D1 * u * v * sqrt(0.2D1) / (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(3,3) = 0.2D1 * u / (w * sqrt(0.2D1) - 0.2D1)

             WorkBasis(4,1:3) = 3.0d0 * WorkBasis(3,1:3) * v
             WorkCurlBasis(4,1) = -0.3D1 / 0.2D1 * v * (0.2D1 * sqrt(0.2D1) * u ** 2 - &
                 0.3D1 * sqrt(0.2D1) * w ** 2 - 0.6D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (w * sqrt(0.2D1) - 0.2D1) ** 2
             WorkCurlBasis(4,2) = 0.6D1 * sqrt(0.2D1) * v ** 2 * u / (w * sqrt(0.2D1) - 0.2D1)**2
             WorkCurlBasis(4,3) = 0.6D1 * u * v / (w * sqrt(0.2D1) - 0.2D1)

             ! -------------------------------------------------------------------
             ! Finally apply an order change and sign reversions if needed. 
             ! -------------------------------------------------------------------
             FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(17,:) = D1 * WorkBasis(2*(I1-1)+1,:)
             CurlBasis(17,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
             EdgeBasis(18,:) = WorkBasis(2*(I1-1)+2,:)
             CurlBasis(18,:) = WorkCurlBasis(2*(I1-1)+2,:)
             EdgeBasis(19,:) = D2 * WorkBasis(2*(I2-1)+1,:)
             CurlBasis(19,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
             EdgeBasis(20,:) = WorkBasis(2*(I2-1)+2,:)
             CurlBasis(20,:) = WorkCurlBasis(2*(I2-1)+2,:) 

             
             !-------------------------------------------------
             ! Two basis functions defined on the face 125:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 1,2,5 /)           
             FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             WorkBasis(1,1:3) = LBasis(5) * EdgeSign(1) * EdgeBasis(1,1:3)
             WorkCurlBasis(1,1) = w * u / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,2) = (-0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - & 
                 0.4D1 * v * w - 0.2D1 * sqrt(0.2D1) + 0.8D1 * w) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(3) * EdgeSign(9) * EdgeBasis(9,1:3)
             WorkCurlBasis(2,1) = (sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 - &
                 0.2D1 * u * w - 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,2) = -(-0.3D1 * sqrt(0.2D1) * u * w ** 2 + 0.2D1 * sqrt(0.2D1) * &
                 v * w ** 2 + 0.6D1 * u * v * sqrt(0.2D1) - 0.7D1 * sqrt(0.2D1) * w ** 2 - &
                 0.8D1 * u * v * w + 0.3D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + 0.2D1 * v * sqrt(0.2D1) + &
                 0.12D2 * u * w - 0.6D1 * v * w - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)**2 )
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.16D2 )

             WorkBasis(3,1:3) = Beta(1) * EdgeSign(11) * EdgeBasis(11,1:3)
             WorkCurlBasis(3,1) = (-sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 + &
                 0.2D1 * u * w - 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)** 2 ) 
             WorkCurlBasis(3,2) = -(-0.3D1 * sqrt(0.2D1) * u * w ** 2 - 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * u * v * sqrt(0.2D1) + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) - 0.2D1 * v * sqrt(0.2D1) + 0.12D2 * u * w + &
                 0.6D1 * v * w + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1)**2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             IF (RedefineFaceBasis) THEN
               EdgeBasis(21,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
               CurlBasis(21,:) = 0.5d0 * D1 * WorkCurlBasis(I1,:) + 0.5d0 * D2 * WorkCurlBasis(I2,:)
               EdgeBasis(22,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
               CurlBasis(22,:) = 0.5d0 * D2 * WorkCurlBasis(I2,:) - 0.5d0 * D1 * WorkCurlBasis(I1,:)
             ELSE             
               EdgeBasis(21,:) = D1 * WorkBasis(I1,:)
               CurlBasis(21,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(22,:) = D2 * WorkBasis(I2,:)
               CurlBasis(22,:) = D2 * WorkCurlBasis(I2,:)              
             END IF
               
             !-------------------------------------------------
             ! Two basis functions defined on the face 235:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 2,3,5 /)          
             FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             WorkBasis(1,1:3) = LBasis(5) * EdgeSign(3) * EdgeBasis(3,1:3)
             WorkCurlBasis(1,1) = (0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - 0.4D1 * u * w + &
                 0.2D1 * sqrt(0.2D1) - 0.8D1 * w) / ( (w * sqrt(0.2D1) - 0.2D1) * 0.8D1 )
             WorkCurlBasis(1,2) = w * v / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(4) * EdgeSign(11) * EdgeBasis(11,1:3)
             WorkCurlBasis(2,1) = -(0.2D1 * sqrt(0.2D1) * u * w ** 2 + 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 + 0.2D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) - 0.6D1 * u * w - &
                 0.12D2 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2) 
             WorkCurlBasis(2,2) = (sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 )
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(2) * EdgeSign(13) * EdgeBasis(13,1:3)
             WorkCurlBasis(3,1) = -(-0.2D1 * sqrt(0.2D1) * u * w ** 2 + 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 - 0.2D1 * u * sqrt(0.2D1) + 0.6D1 * v * sqrt(0.2D1) + 0.6D1 * u * w - &
                 0.12D2 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = (-sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 ( (w * sqrt(0.2D1) - 0.2D1) * 0.16D2 )

             IF (RedefineFaceBasis) THEN
               EdgeBasis(23,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
               CurlBasis(23,:) = 0.5d0 * D1 * WorkCurlBasis(I1,:) + 0.5d0 * D2 * WorkCurlBasis(I2,:)
               EdgeBasis(24,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
               CurlBasis(24,:) = 0.5d0 * D2 * WorkCurlBasis(I2,:) - 0.5d0 * D1 * WorkCurlBasis(I1,:)
             ELSE
               EdgeBasis(23,:) = D1 * WorkBasis(I1,:)
               CurlBasis(23,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(24,:) = D2 * WorkBasis(I2,:)
               CurlBasis(24,:) = D2 * WorkCurlBasis(I2,:)
             END IF

             !-------------------------------------------------
             ! Two basis functions defined on the face 345:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 3,4,5 /)           
             FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             WorkBasis(1,1:3) = -LBasis(5) * EdgeSign(5) * EdgeBasis(5,1:3)
             WorkCurlBasis(1,1) = w * u / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,2) = (0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * v * sqrt(0.2D1) - 0.4D1 * w * v + &
                 0.2D1 * sqrt(0.2D1) - 0.8D1 * w) / (0.8D1 * (w * sqrt(0.2D1)- 0.2D1) )
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(1) * EdgeSign(13) * EdgeBasis(13,1:3)
             WorkCurlBasis(2,1) = -(-sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * u * w - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,2) = (0.3D1 * sqrt(0.2D1) * u * w ** 2 - 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) - 0.2D1 * v * sqrt(0.2D1) - 0.12D2 * u * w + &
                 0.6D1 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(3) * EdgeSign(15) * EdgeBasis(15,1:3)
             WorkCurlBasis(3,1) = -(sqrt(0.2D1) * u * w ** 2 + 0.4D1 * sqrt(0.2D1) * u ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * u * w - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = (0.3D1 * sqrt(0.2D1) * u * w ** 2 + 0.2D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 + 0.6D1 * u * sqrt(0.2D1) + 0.2D1 * v * sqrt(0.2D1) - 0.12D2 * u * w - &
                 0.6D1 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * u - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             IF (RedefineFaceBasis) THEN
               EdgeBasis(25,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
               CurlBasis(25,:) = 0.5d0 * D1 * WorkCurlBasis(I1,:) + 0.5d0 * D2 * WorkCurlBasis(I2,:)
               EdgeBasis(26,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
               CurlBasis(26,:) = 0.5d0 * D2 * WorkCurlBasis(I2,:) - 0.5d0 * D1 * WorkCurlBasis(I1,:)
             ELSE
               EdgeBasis(25,:) = D1 * WorkBasis(I1,:)
               CurlBasis(25,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(26,:) = D2 * WorkBasis(I2,:)
               CurlBasis(26,:) = D2 * WorkCurlBasis(I2,:)              
             END IF
               
             !-------------------------------------------------
             ! Two basis functions defined on the face 415:
             !-------------------------------------------------
             TriangleFaceMap(:) = (/ 4,1,5 /)          
             FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
             CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             WorkBasis(1,1:3) = -LBasis(5) * EdgeSign(7) * EdgeBasis(7,1:3)
             WorkCurlBasis(1,1) = (-0.3D1 * sqrt(0.2D1) * w ** 2 + 0.2D1 * u * sqrt(0.2D1) - &
                 0.4D1 * u * w - 0.2D1 * sqrt(0.2D1) + 0.8D1 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) )
             WorkCurlBasis(1,2) = w * v / (w * sqrt(0.2D1) - 0.2D1) / 0.4D1
             WorkCurlBasis(1,3) = w * sqrt(0.2D1) / 0.8D1

             WorkBasis(2,1:3) = Beta(2) * EdgeSign(15) * EdgeBasis(15,1:3)
             WorkCurlBasis(2,1) = (-0.2D1 * sqrt(0.2D1) * u * w ** 2 - 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v + 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w - &
                 0.3D1 * w ** 3 - 0.2D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) + 0.6D1 * u * w + &
                 0.12D2 * w * v + 0.2D1 * sqrt(0.2D1) - 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 )
             WorkCurlBasis(2,2) = -(-sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 + 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(2,3) = w * sqrt(0.2D1) * (w * sqrt(0.2D1) - 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             WorkBasis(3,1:3) = Beta(4) * EdgeSign(9) * EdgeBasis(9,1:3)
             WorkCurlBasis(3,1) = (0.2D1 * sqrt(0.2D1) * u * w ** 2 - 0.3D1 * sqrt(0.2D1) * v * w ** 2 + &
                 0.6D1 * sqrt(0.2D1) * u * v - 0.7D1 * sqrt(0.2D1) * w ** 2 - 0.8D1 * u * v * w + &
                 0.3D1 * w ** 3 + 0.2D1 * u * sqrt(0.2D1) - 0.6D1 * v * sqrt(0.2D1) - 0.6D1 * u * w + &
                 0.12D2 * w * v - 0.2D1 * sqrt(0.2D1) + 0.10D2 * w) / &
                 (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,2) = -(sqrt(0.2D1) * v * w ** 2 + 0.4D1 * sqrt(0.2D1) * v ** 2 - &
                 0.8D1 * sqrt(0.2D1) * w ** 2 - 0.4D1 * v ** 2 * w + 0.3D1 * w ** 3 - 0.2D1 * w * v - &
                 0.4D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.8D1 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             WorkCurlBasis(3,3) = -w * sqrt(0.2D1) * (w * sqrt(0.2D1) + 2.0D0 * v - 0.2D1) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 

             IF (RedefineFaceBasis) THEN
               EdgeBasis(27,:) = 0.5d0 * D1 * WorkBasis(I1,:) + 0.5d0 * D2 * WorkBasis(I2,:)
               CurlBasis(27,:) = 0.5d0 * D1 * WorkCurlBasis(I1,:) + 0.5d0 * D2 * WorkCurlBasis(I2,:)
               EdgeBasis(28,:) = 0.5d0 * D2 * WorkBasis(I2,:) - 0.5d0 * D1 * WorkBasis(I1,:)
               CurlBasis(28,:) = 0.5d0 * D2 * WorkCurlBasis(I2,:) - 0.5d0 * D1 * WorkCurlBasis(I1,:)
             ELSE
               EdgeBasis(27,:) = D1 * WorkBasis(I1,:)
               CurlBasis(27,:) = D1 * WorkCurlBasis(I1,:)
               EdgeBasis(28,:) = D2 * WorkBasis(I2,:)
               CurlBasis(28,:) = D2 * WorkCurlBasis(I2,:)              
             END IF

             ! Finally three interior basis functions:
             ! -----------------------------------------------------------------------------------
             EdgeBasis(29,1:3) = LBasis(5) * Beta(4) * EdgeSign(1) * EdgeBasis(1,1:3)
             CurlBasis(29,1) = u * v * w / (0.4D1 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(29,2) = (0.2D1 * sqrt(0.2D1) * v ** 2 - 0.9D1 * sqrt(0.2D1) * w ** 2 - &
                 0.4D1 * v ** 2 * w + 0.4D1 * w ** 3 - 0.2D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(29,3) = sqrt(0.2D1) * v * w / 0.8D1

             EdgeBasis(30,1:3) = LBasis(5) * Beta(3) * EdgeSign(7) * EdgeBasis(7,1:3)
             CurlBasis(30,1) = -(0.2D1 * sqrt(0.2D1) * u ** 2 - 0.9D1 * sqrt(0.2D1) * w **2 - &
                 0.4D1 * u ** 2 * w + 0.4D1 * w ** 3 - 0.2D1 * sqrt(0.2D1) + 0.12D2 * w) / &
                 (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(30,2) = -u * v * w / (0.4D1* (w * sqrt(0.2D1) - 0.2D1) ) 
             CurlBasis(30,3) = -sqrt(0.2D1) * u * w / 0.8D1

             EdgeBasis(31,1:3) = Beta(3) * Beta(4) * EdgeSign(9) * EdgeBasis(9,1:3)
             CurlBasis(31,1) = (0.2D1 * sqrt(0.2D1) * u ** 2 * w ** 2 + 0.2D1 * sqrt(0.2D1) * u * v * w ** 2 -&
                 0.2D1 * sqrt(0.2D1) * w ** 4 + 0.6D1 * sqrt(0.2D1) * u ** 2 * v - &
                 0.11D2 * sqrt(0.2D1) * v * w ** 2 - 0.8D1 * u ** 2 * v * w + 0.4D1 * v * w ** 3 + &
                 0.2D1 * sqrt(0.2D1) * u ** 2 - 0.15D2 * sqrt(0.2D1) * w ** 2 - 0.6D1 * u ** 2 * w - &
                 0.4D1 * u * v * w + 0.13D2 * w ** 3 - 0.6D1 * v * sqrt(0.2D1) + 0.20D2 * w * v - &
                 0.2D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             CurlBasis(31,2) = -(0.2D1 * sqrt(0.2D1) * u * v * w ** 2 + 0.2D1 * sqrt(0.2D1) * v ** 2 * w**2 - &
                 0.2D1 * sqrt(0.2D1) * w ** 4 + 0.6D1 * sqrt(0.2D1) * u * v ** 2 - &
                 0.11D2 * sqrt(0.2D1) * u * w ** 2 - 0.8D1 * u * v ** 2 * w + 0.4D1 * u * w ** 3 + &
                 0.2D1 * sqrt(0.2D1) * v ** 2 - 0.15D2 * sqrt(0.2D1) * w ** 2 - 0.4D1 * u * v * w - &
                 0.6D1 * v ** 2 * w + 0.13D2 * w ** 3 - 0.6D1 * u * sqrt(0.2D1) + 0.20D2 * u *w - &
                 0.2D1 * sqrt(0.2D1) + 0.14D2 * w) / (0.16D2 * (w * sqrt(0.2D1) - 0.2D1) ** 2 ) 
             CurlBasis(31,3) = -(u - v) * w * sqrt(0.2D1) / 0.16D2

             ! Finally, scale to reduce ill-conditioning:
             IF (ScaleFaceBasis) THEN
               EdgeBasis(21:27:2,:) = sqrt(fs1) * EdgeBasis(21:27:2,:)
               CurlBasis(21:27:2,:) = sqrt(fs1) * CurlBasis(21:27:2,:)
               EdgeBasis(22:28:2,:) = sqrt(fs2) * EdgeBasis(22:28:2,:)
               CurlBasis(22:28:2,:) = sqrt(fs2) * CurlBasis(22:28:2,:)                 

               EdgeBasis(29:30,:) = sqrt(506.9d0) * EdgeBasis(29:30,:)
               CurlBasis(29:30,:) = sqrt(506.9d0) * CurlBasis(29:30,:)
               EdgeBasis(31,:) = sqrt(167.8d0) * EdgeBasis(31,:)
               CurlBasis(31,:) = sqrt(167.8d0) * CurlBasis(31,:)
             END IF
             
           ELSE
             !-----------------------------------------------------------------------------------------
             ! The lowest-order pyramid from the optimal family. Now these basis functions are 
             ! also contained in the set of hierarchic basis functions, so this branch could be
             ! removed by making some code modifications (to do?).
             !-----------------------------------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = (v*(-1 + (2*v)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = (u*v*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,1) = (u*(-Sqrt(2.0d0) + 2*Sqrt(2.0d0)*v + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,2) = (v*(Sqrt(2.0d0) - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(1,3) = (-2 + 4*v + Sqrt(2.0d0)*w)/(-8 + 4*Sqrt(2.0d0)*w)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1 + (2*u)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(2,3) = (u*v*(Sqrt(2.0d0) + Sqrt(2.0d0)*u - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,1) = (u*(Sqrt(2.0d0) - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,2) = -(v*(Sqrt(2.0d0) + 2*Sqrt(2.0d0)*u - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(2,3) = (2 + 4*u - Sqrt(2.0d0)*w)/(8 - 4*Sqrt(2.0d0)*w)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(3,1) = (v*(1 + (2*v)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(3,2) = 0.0d0
             EdgeBasis(3,3) = (u*v*(Sqrt(2.0d0) + Sqrt(2.0d0)*v - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,1) = (u*(Sqrt(2.0d0) + 2*Sqrt(2.0d0)*v - w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,2) = (v*(-Sqrt(2.0d0) + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(3,3) = (2 + 4*v - Sqrt(2.0d0)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = (u*(-1 + (2*u)/(2 - Sqrt(2.0d0)*w)))/4.0d0
             EdgeBasis(4,3) = (u*v*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,1) = (u*(-Sqrt(2.0d0) + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,2) = -(v*(-Sqrt(2.0d0) + 2*Sqrt(2.0d0)*u + w))/(2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(4,3) = (2 - 4*u - Sqrt(2.0d0)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             EdgeBasis(5,1) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(5,2) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(5,3) = (u*(-2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) + 4*w - Sqrt(2.0d0)*w**2) - &
                 (-1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2))/(4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,1) = (-2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) + 4*w - Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,2) = (2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*v - 4*w + 2*v*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(5,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             EdgeBasis(6,1) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*v + w))/(8.0d0 - 4*Sqrt(2.0d0)*w)
             EdgeBasis(6,2) = (w*(-Sqrt(2.0d0) - Sqrt(2.0d0)*u + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(6,3) = (-((-1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2)) + & 
                 u*(2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*v - 4*w + 4*v*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(6,1) = -(2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(6,2) = (-2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) + 4*w - Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2) 
             CurlBasis(6,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             EdgeBasis(7,1) = ((Sqrt(2.0d0) + Sqrt(2.0d0)*v - w)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(7,2) = ((Sqrt(2.0d0) + Sqrt(2.0d0)*u - w)*w)/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(7,3) = ((1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2) + &
                 u*(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) - 4*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,1) = (2*Sqrt(2.0d0) + 2*u*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,2) = -(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2 + Sqrt(2.0d0)*w)**2)
             CurlBasis(7,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             EdgeBasis(8,1) = (w*(-Sqrt(2.0d0) - Sqrt(2.0d0)*v + w))/(-8.0d0 + 4*Sqrt(2.0d0)*w)
             EdgeBasis(8,2) = (w*(-Sqrt(2.0d0) + Sqrt(2.0d0)*u + w))/(8.0d0 - 4*Sqrt(2.0d0)*w)
             EdgeBasis(8,3) = ((1 + v)*(2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2) - &
                 u*(2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - 2*w) - 4*w + Sqrt(2.0d0)*w**2))/ &
                 (4.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,1) = (2*Sqrt(2.0d0) - 2*Sqrt(2.0d0)*u - 4*w + 2*u*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,2) = (2*Sqrt(2.0d0) + 2*v*(Sqrt(2.0d0) - w) - 4*w + Sqrt(2.0d0)*w**2)/ &
                 (2.0d0*(-2.0d0 + Sqrt(2.0d0)*w)**2)
             CurlBasis(8,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             ! ------------------------------------------------------------------
             ! The last two basis functions are associated with the square face.
             ! We first create the basis function in the default order without
             ! sign reversions.
             ! ------------------------------------------------------------------
             SquareFaceMap(:) = (/ 1,2,3,4 /)

             WorkBasis(1,1) = (2.0d0 - 2*v**2 - 2*Sqrt(2.0d0)*w + w**2)/(4.0d0 - 2*Sqrt(2.0d0)*w)
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = (u*(1.0d0 - (4*v**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2))/(2.0d0*Sqrt(2.0d0))
             WorkCurlBasis(1,1) = (-2*Sqrt(2.0d0)*u*v)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(1,2) = (-2*Sqrt(2.0d0) + 4*w - Sqrt(2.0d0)*w**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(1,3) = (2.0d0*v)/(2.0d0 - Sqrt(2.0d0)*w)

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = (2.0d0 - 2*u**2 - 2*Sqrt(2.0d0)*w + w**2)/(4.0d0 - 2*Sqrt(2.0d0)*w)
             WorkBasis(2,3) = (v*(1.0d0 - (4*u**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2))/(2.0d0*Sqrt(2.0d0))
             WorkCurlBasis(2,1) = (2*Sqrt(2.0d0) - 4*w + Sqrt(2.0d0)*w**2)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(2,2) = (2*Sqrt(2.0d0)*u*v)/(-2.0d0 + Sqrt(2.0d0)*w)**2
             WorkCurlBasis(2,3) = (2*u)/(-2.0d0 + Sqrt(2.0d0)*w)

             ! -------------------------------------------------------------------
             ! Finally apply an order change and sign reversions if needed. 
             ! -------------------------------------------------------------------
             FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(9,:) = D1 * WorkBasis(I1,:)
             CurlBasis(9,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(10,:) = D2 * WorkBasis(I2,:)
             CurlBasis(10,:) = D2 * WorkCurlBasis(I2,:)          
           END IF

         CASE(7)
           !--------------------------------------------------------------
           ! This branch is for handling prismatic (or wedge) elements
           !--------------------------------------------------------------           
           EdgeMap => GetEdgeMap(7)

           IF (SecondOrder) THEN
             GRADIENT_VERSION_PRISM: IF (GradVersion) THEN
               EDOFs = 2               
               !
               ! First handle the edges which bound a triangular face
               !
               EDGES_BOUNDING_TRIANGLE: DO k=1,3
                 !
                 ! We utilize the Nedelec basis for the triangle, so we
                 ! compute component functions to create the Whitney forms
                 ! associated with the edges of triangle.
                 
                 EdgeMap => GetEdgeMap(3)
                 i = EdgeMap(k,1)
                 j = EdgeMap(k,2)

                 CALL EdgeWhitneyComponents2D(Wrk(1:2,:), WrkCurl(1:2,:), i, j, u, v)
                 bfun = TriangleNodalPBasis(i,u,v) * TriangleNodalPBasis(j,u,v)

                 EdgeMap => GetEdgeMap(7)

                 TRIANGULAR_ENDS: DO q=1,2
                   !
                   ! If k=1, handle the first and fourth edge when q=1 and q=2, respectively.
                   ! If k=2, handle the second and fifth edge.
                   ! If k=3, handle the third and sixth edge.
                   !
                   SELECT CASE(q)
                   CASE(1)
                     ! Pick a blending function:
                     h1 = 0.5d0 * (1-w)
                     dh1 = -0.5d0
                   CASE(2)
                     i = EdgeMap(3+k,1)
                     j = EdgeMap(3+k,2)
                     ! Pick a blending function:
                     h1 = 0.5d0 * (1+w)
                     dh1 = 0.5d0                     
                   END SELECT

                   IF (GIndexes(j) < GIndexes(i)) THEN
                     I1 = 2
                     I2 = 1
                   ELSE
                     I1 = 1
                     I2 = 2
                   END IF

                   DO l=1,EDOFs
                     SELECT CASE(l)
                     CASE(1)
                       sfun = -1.0d0
                       tfun = 1.0d0
                     CASE(2)
                       sfun = 1.0d0
                       tfun = 1.0d0
                     CASE DEFAULT
                       CALL Fatal('ElementDescription::EdgeElementInfo','sfun/tfun not defined')
                     END SELECT

                     B(1:2) = sfun * Wrk(I1,1:2) + tfun * Wrk(I2,1:2)
                     CurlB(3) = sfun * WrkCurl(I1,3) + tfun * WrkCurl(I2,3)
                     EdgeBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 1:2) = B(1:2) * h1
                     
                     SELECT CASE(l)
                     CASE(1)
                       CurlBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 1) = -B(2) * dh1
                       CurlBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 2) = B(1) * dh1
                       CurlBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 3) = CurlB(3) * h1 
                     CASE(2)
                       ! The basis function obtained as grad(bfun*h1)
                       EdgeBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 3) = bfun * dh1
                       CurlBasis((q-1)*3*EDOFs + EDOFs*(k-1) + l, 1:3) = 0.0d0
                     END SELECT
                   END DO
                 END DO TRIANGULAR_ENDS
               END DO EDGES_BOUNDING_TRIANGLE

               AXIAL_EDGES: DO k=1,3

                 i = EdgeMap(6+k,1)
                 j = EdgeMap(6+k,2)
                 grad(1:2) = dTriangleNodalPBasis(k, u, v)

                 WorkBasis(1,3) = 0.5d0 * TriangleNodalPBasis(k, u, v)
                 WorkCurlBasis(1,1) = 0.5d0* grad(2)
                 WorkCurlBasis(1,2) = -0.5d0* grad(1)

                 ! (u,v,w) -> grad( -1/sqrt(6) * L_k(u,v) * phi_2(w) )
                 !          = grad( L_k(u,v) * L_1(w) * L_2(w) )
                 WorkBasis(2,1) = -1.0d0/sqrt(6.0d0) * grad(1) * Phi(2,w)
                 WorkBasis(2,2) = -1.0d0/sqrt(6.0d0) * grad(2) * Phi(2,w)
                 WorkBasis(2,3) = -1.0d0/sqrt(6.0d0) * TriangleNodalPBasis(k, u, v) * dPhi(2,w)

                 IF (GIndexes(j) < GIndexes(i)) THEN
                   WorkBasis(1,3) = -WorkBasis(1,3)
                   WorkCurlBasis(1,1:2) = -WorkCurlBasis(1,1:2)
                 END IF
                 
                 EdgeBasis(6*EDOFs+(k-1)*EDOFs+1,3) = WorkBasis(1,3)
                 CurlBasis(6*EDOFs+(k-1)*EDOFs+1,1:2) = WorkCurlBasis(1,1:2)
                 
                 EdgeBasis(6*EDOFs+(k-1)*EDOFs+2,1:3) = WorkBasis(2,1:3)
                 CurlBasis(6*EDOFs+(k-1)*EDOFs+2,1:3) = 0.0d0 
               END DO AXIAL_EDGES

               ! The triangular faces and two internal (bubble) functions
               !
               FDOFs = 2
               TRIANGULAR_FACES: DO k=1,3
                 !
                 ! We utilize the basis functions of the triangle
                 !
                 CALL FaceWhitneyComponents2D(Wrk(1:3,:), WrkCurl(1:3,:), u, v)
                 
                 SELECT CASE(k)
                 CASE(1)
                   TriangleFaceMap(:) = (/ 1,2,3 /)
                   h1 = 0.5d0 * (1-w)
                   dh1 = -0.5d0
                 CASE(2)
                   TriangleFaceMap(:) = (/ 4,5,6 /)
                   h1 = 0.5d0 * (1+w)
                   dh1 = 0.5d0
                 CASE(3)
                   h1 = 0.25d0 * (1-w**2)
                   dh1 = -0.5d0 * w
                 END SELECT

                 IF (k == 3) THEN
                   I1 = 1
                   I2 = 2
                   I3 = 3
                   ! Check the following offset if the order is higher than 2
                   c0 = 9*EDOFs + 2*FDOFs + 3*4
                 ELSE
                   FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
                   CALL TriangleFaceDofsOrdering2nd(I1,I2,I3,FaceIndices(1:3))
                   c0 = 9*EDOFs + (k-1)*FDOFs
                 END IF
                   
                 DO l=1,FDOFs
                   
                   SELECT CASE(l)
                   CASE(1)
                     sfun = 1.0d0
                     tfun = 1.0d0
                     hfun = -2.0d0
                   CASE(2)
                     sfun = 1.0d0
                     tfun = -1.0d0
                     hfun = 0.0d0
                   END SELECT

                   B(1:2) = sfun * Wrk(I1,1:2) + tfun * Wrk(I2,1:2) + hfun * Wrk(I3,1:2)
                   CurlB(3) = sfun * WrkCurl(I1,3) + tfun * WrkCurl(I2,3) + &
                       hfun * WrkCurl(I3,3)

                   EdgeBasis(c0+l,1:2) = B(1:2) * h1
                   CurlBasis(c0+l,1) = -B(2) * dh1
                   CurlBasis(c0+l,2) = B(1) * dh1
                   CurlBasis(c0+l,3) = CurlB(3) * h1
                 END DO
               END DO TRIANGULAR_FACES
                 
               ! The quadrilateral faces
               !
               FDOFs = 4
               c0 = 9*EDOFs + 2*2
               QUAD_FACES_PRISM: DO k=1,3

                 SELECT CASE(k)
                 CASE(1)
                   SquareFaceMap(:) = (/ 1,2,5,4 /)
                 CASE(2)
                   SquareFaceMap(:) = (/ 2,3,6,5 /)
                 CASE(3)
                   SquareFaceMap(:) = (/ 3,1,4,6 /)
                 END SELECT

                 h1 = 1.0d0 - w**2
                 dh1 = -2.0d0 * w
                 
                 i = SquareFaceMap(1)
                 j = SquareFaceMap(2)

                 CALL EdgeWhitneyComponents2D(Wrk(1:2,:), WrkCurl(1:2,:), i, j, u, v)

                 WorkBasis(:,:) = 0.0d0
                 WorkCurlBasis(:,:) = 0.0d0

                 ! The case where blending is done in the direction of the axis of prism:
                 DO l=1,FDOFs/2
                   
                   SELECT CASE(l)
                   CASE(1)
                     sfun = -1.0d0
                     tfun = 1.0d0
                   CASE(2)
                     sfun = 1.0d0
                     tfun = 1.0d0
                   CASE DEFAULT
                     CALL Fatal('ElementDescription::EdgeElementInfo','sfun/tfun not defined')
                   END SELECT

                   B(1:2) = sfun * Wrk(1,1:2) + tfun * Wrk(2,1:2)
                   CurlB(3) = sfun * WrkCurl(1,3) + tfun * WrkCurl(2,3)

                   WorkBasis(2*(l-1)+1,1:2) = B(1:2) * h1
                   
                   WorkCurlBasis(2*(l-1)+1,1) = -B(2) * dh1
                   WorkCurlBasis(2*(l-1)+1,2) = B(1) * dh1
                   WorkCurlBasis(2*(l-1)+1,3) = CurlB(3) * h1 
                 END DO

                 grad_i = dTriangleNodalPBasis(i,u,v)
                 grad_j = dTriangleNodalPBasis(j,u,v)
                 bfun = TriangleNodalPBasis(i,u,v) * TriangleNodalPBasis(j,u,v)
                 
                 WorkBasis(2,3) = 2.0d0 * bfun
                 WorkCurlBasis(2,1) = 2.0d0 * (grad_i(2)*TriangleNodalPBasis(j,u,v) + &
                     TriangleNodalPBasis(i,u,v)*grad_j(2))
                 WorkCurlBasis(2,2) = -2.0d0 * (grad_i(1)*TriangleNodalPBasis(j,u,v) + &
                     TriangleNodalPBasis(i,u,v)*grad_j(1))
                 
                 WorkBasis(4,3) = dh1 * bfun
                 WorkCurlBasis(4,1) = -0.5d0 * w * 4.0d0 * (grad_i(2)*TriangleNodalPBasis(j,u,v) + &
                     TriangleNodalPBasis(i,u,v)*grad_j(2))
                 WorkCurlBasis(4,2) = 0.5d0 * w * 4.0d0 * (grad_i(1)*TriangleNodalPBasis(j,u,v) + &
                     TriangleNodalPBasis(i,u,v)*grad_j(1))

                 FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
                 CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

                 EdgeBasis(c0 + (k-1)*FDOFs + 1,1:3) = D1 * WorkBasis(I1,1:3)
                 CurlBasis(c0 + (k-1)*FDOFs + 1,1:3) = D1 * WorkCurlBasis(I1,1:3)
                 EdgeBasis(c0 + (k-1)*FDOFs + 2,1:3) = D2 * WorkBasis(I2,1:3)
                 CurlBasis(c0 + (k-1)*FDOFs + 2,1:3) = D2 * WorkCurlBasis(I2,1:3)
                 ! Note that sign changes never happen:
                 EdgeBasis(c0 + (k-1)*FDOFs + 3,1:3) = WorkBasis(2+I1,1:3)
                 CurlBasis(c0 + (k-1)*FDOFs + 3,1:3) = WorkCurlBasis(2+I1,1:3)
                 EdgeBasis(c0 + (k-1)*FDOFs + 4,1:3) = WorkBasis(2+I1,1:3) + WorkBasis(2+I2,1:3)
                 CurlBasis(c0 + (k-1)*FDOFs + 4,1:3) = 0.0d0
                 !CurlBasis(c0 + (k-1)*FDOFs + 4,1:3) = WorkCurlBasis(2+I1,1:3) + WorkCurlBasis(2+I2,1:3)

               END DO QUAD_FACES_PRISM
               
             ELSE
               !---------------------------------------------------------------
               ! The second-order element from the Nedelec's first family 
               ! (note that the lowest-order prism element is from a different 
               ! family). This element may not be optimally accurate if 
               ! the physical element is not affine.
               !--------------------------------------------------------------             
               h1 = 0.5d0 * (1-w)
               dh1 = -0.5d0
               h2 = 0.5d0 * (1+w)
               dh2 = 0.5d0
               h3 = h1 * h2
               dh3 = -0.5d0 * w

               ! ---------------------------------------------------------
               ! The first and fourth edges ...
               !--------------------------------------------------------
               ! The corresponding basis functions for the triangle:
               !--------------------------------------------------------
               WorkBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0
               WorkBasis(1,2) = u/(2.0d0*Sqrt(3.0d0))
               WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
               WorkBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0
               WorkBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0
               WorkCurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0

               i = EdgeMap(1,1)
               j = EdgeMap(1,2)
               EdgeBasis(1,1:2) = WorkBasis(1,1:2) * h1
               CurlBasis(1,1) = -WorkBasis(1,2) * dh1
               CurlBasis(1,2) = WorkBasis(1,1) * dh1
               CurlBasis(1,3) = WorkCurlBasis(1,3) * h1
               EdgeBasis(2,1:2) = WorkBasis(2,1:2) * h1
               CurlBasis(2,1) = -WorkBasis(2,2) * dh1
               CurlBasis(2,2) = WorkBasis(2,1) * dh1
               CurlBasis(2,3) = WorkCurlBasis(2,3) * h1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(1,1:2) = -EdgeBasis(1,1:2)
                 CurlBasis(1,1:3) = -CurlBasis(1,1:3)
               END IF

               i = EdgeMap(4,1)
               j = EdgeMap(4,2)
               EdgeBasis(7,1:2) = WorkBasis(1,1:2) * h2
               CurlBasis(7,1) = -WorkBasis(1,2) * dh2
               CurlBasis(7,2) = WorkBasis(1,1) * dh2
               CurlBasis(7,3) = WorkCurlBasis(1,3) * h2
               EdgeBasis(8,1:2) = WorkBasis(2,1:2) * h2
               CurlBasis(8,1) = -WorkBasis(2,2) * dh2
               CurlBasis(8,2) = WorkBasis(2,1) * dh2
               CurlBasis(8,3) = WorkCurlBasis(2,3) * h2
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(7,1:2) = -EdgeBasis(7,1:2)
                 CurlBasis(7,1:3) = -CurlBasis(7,1:3)
               END IF

               ! ---------------------------------------------------------
               ! The second and fifth edges ...
               !--------------------------------------------------------
               ! The corresponding basis functions for the triangle:
               !--------------------------------------------------------
               WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0))
               WorkBasis(1,2) = (1 + u)/(2.0d0*Sqrt(3.0d0))
               WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0)
               WorkBasis(2,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0
               WorkBasis(2,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0
               WorkCurlBasis(2,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0

               i = EdgeMap(2,1)
               j = EdgeMap(2,2)
               EdgeBasis(3,1:2) = WorkBasis(1,1:2) * h1
               CurlBasis(3,1) = -WorkBasis(1,2) * dh1
               CurlBasis(3,2) = WorkBasis(1,1) * dh1
               CurlBasis(3,3) = WorkCurlBasis(1,3) * h1
               EdgeBasis(4,1:2) = WorkBasis(2,1:2) * h1
               CurlBasis(4,1) = -WorkBasis(2,2) * dh1
               CurlBasis(4,2) = WorkBasis(2,1) * dh1
               CurlBasis(4,3) = WorkCurlBasis(2,3) * h1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(3,1:2) = -EdgeBasis(3,1:2)
                 CurlBasis(3,1:3) = -CurlBasis(3,1:3)
               END IF

               i = EdgeMap(5,1)
               j = EdgeMap(5,2)
               EdgeBasis(9,1:2) = WorkBasis(1,1:2) * h2
               CurlBasis(9,1) = -WorkBasis(1,2) * dh2
               CurlBasis(9,2) = WorkBasis(1,1) * dh2
               CurlBasis(9,3) = WorkCurlBasis(1,3) * h2
               EdgeBasis(10,1:2) = WorkBasis(2,1:2) * h2
               CurlBasis(10,1) = -WorkBasis(2,2) * dh2
               CurlBasis(10,2) = WorkBasis(2,1) * dh2
               CurlBasis(10,3) = WorkCurlBasis(2,3) * h2
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(9,1:2) = -EdgeBasis(9,1:2)
                 CurlBasis(9,1:3) = -CurlBasis(9,1:3)
               END IF

               ! ---------------------------------------------------------
               ! The third and sixth edges ...
               !--------------------------------------------------------
               ! The corresponding basis functions for the triangle:
               !--------------------------------------------------------
               WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0))
               WorkBasis(1,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0))
               WorkCurlBasis(1,3) =  1.0d0/Sqrt(3.0d0)
               WorkBasis(2,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0
               WorkBasis(2,2) = -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0
               WorkCurlBasis(2,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0

               i = EdgeMap(3,1)
               j = EdgeMap(3,2)
               EdgeBasis(5,1:2) = WorkBasis(1,1:2) * h1
               CurlBasis(5,1) = -WorkBasis(1,2) * dh1
               CurlBasis(5,2) = WorkBasis(1,1) * dh1
               CurlBasis(5,3) = WorkCurlBasis(1,3) * h1
               EdgeBasis(6,1:2) = WorkBasis(2,1:2) * h1
               CurlBasis(6,1) = -WorkBasis(2,2) * dh1
               CurlBasis(6,2) = WorkBasis(2,1) * dh1
               CurlBasis(6,3) = WorkCurlBasis(2,3) * h1
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(5,1:2) = -EdgeBasis(5,1:2)
                 CurlBasis(5,1:3) = -CurlBasis(5,1:3)
               END IF

               i = EdgeMap(6,1)
               j = EdgeMap(6,2)
               EdgeBasis(11,1:2) = WorkBasis(1,1:2) * h2
               CurlBasis(11,1) = -WorkBasis(1,2) * dh2
               CurlBasis(11,2) = WorkBasis(1,1) * dh2
               CurlBasis(11,3) = WorkCurlBasis(1,3) * h2
               EdgeBasis(12,1:2) = WorkBasis(2,1:2) * h2
               CurlBasis(12,1) = -WorkBasis(2,2) * dh2
               CurlBasis(12,2) = WorkBasis(2,1) * dh2
               CurlBasis(12,3) = WorkCurlBasis(2,3) * h2
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(11,1:2) = -EdgeBasis(11,1:2)
                 CurlBasis(11,1:3) = -CurlBasis(11,1:3)
               END IF

               ! -------------------------------------------------------
               ! The edges 14, 25 and 36
               !--------------------------------------------------------
               DO q = 1,3
                 i = EdgeMap(6+q,1)
                 j = EdgeMap(6+q,2)
                 grad(1:2) = dTriangleNodalPBasis(q, u, v)
                 EdgeBasis(12+(q-1)*2+1,3) = 0.5d0 * TriangleNodalPBasis(q, u, v)
                 CurlBasis(12+(q-1)*2+1,1) = 0.5d0* grad(2)
                 CurlBasis(12+(q-1)*2+1,2) = -0.5d0* grad(1)
                 EdgeBasis(12+(q-1)*2+2,3) = 3.0d0 * EdgeBasis(12+(q-1)*2+1,3) * w
                 CurlBasis(12+(q-1)*2+2,1) = 1.5d0 * grad(2) * w
                 CurlBasis(12+(q-1)*2+2,2) = -1.5d0 * grad(1) * w

                 IF(GIndexes(j)<GIndexes(i)) THEN
                   EdgeBasis(12+(q-1)*2+1,3) = -EdgeBasis(12+(q-1)*2+1,3)
                   CurlBasis(12+(q-1)*2+1,1:2) = -CurlBasis(12+(q-1)*2+1,1:2)
                 END IF
               END DO

               !-------------------------------------------------
               ! Two basis functions defined on the face 123:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 1,2,3 /)
               FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               WorkBasis(1,1) = ((Sqrt(3.0d0) - v)*v)/6.0d0
               WorkBasis(1,2) = (u*v)/6.0d0
               WorkCurlBasis(1,3) = (-Sqrt(3.0d0) + 3.0d0*v)/6.0d0
               WorkBasis(2,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0))
               WorkBasis(2,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(2,3) =(-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0
               WorkBasis(3,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkBasis(3,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0))
               WorkCurlBasis(3,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

               IF (RedefineFaceBasis) THEN
                 EdgeBasis(19,1:2) = (D1 * WorkBasis(I1,1:2) + D2 * WorkBasis(I2,1:2)) * 0.5d0 * h1
                 EdgeBasis(20,1:2) = (-D1 * WorkBasis(I1,1:2) + D2 * WorkBasis(I2,1:2)) * 0.5d0 * h1

                 CurlBasis(19,1) = -(D1 * WorkBasis(I1,2) + D2 * WorkBasis(I2,2)) * 0.5d0 * dh1
                 CurlBasis(19,2) = (D1 * WorkBasis(I1,1) + D2 * WorkBasis(I2,1)) * 0.5d0 * dh1
                 CurlBasis(19,3) = (D1 * WorkCurlBasis(I1,3) + D2 * WorkCurlBasis(I2,3)) * 0.5d0 * h1

                 CurlBasis(20,1) = -(-D1 * WorkBasis(I1,2) + D2 * WorkBasis(I2,2)) * 0.5d0 * dh1
                 CurlBasis(20,2) = (-D1 * WorkBasis(I1,1) + D2 * WorkBasis(I2,1)) * 0.5d0 * dh1
                 CurlBasis(20,3) = (-D1 * WorkCurlBasis(I1,3) + D2 * WorkCurlBasis(I2,3)) * 0.5d0 * h1

               ELSE
                 EdgeBasis(19,1:2) = D1 * WorkBasis(I1,1:2) * h1
                 CurlBasis(19,1) = -D1 * WorkBasis(I1,2) * dh1
                 CurlBasis(19,2) = D1 * WorkBasis(I1,1) * dh1
                 CurlBasis(19,3) = D1 * WorkCurlBasis(I1,3) * h1

                 EdgeBasis(20,1:2) = D2 * WorkBasis(I2,1:2) * h1
                 CurlBasis(20,1) = -D2 * WorkBasis(I2,2) * dh1
                 CurlBasis(20,2) = D2 * WorkBasis(I2,1) * dh1
                 CurlBasis(20,3) = D2 * WorkCurlBasis(I2,3) * h1
               END IF

               !-------------------------------------------------
               ! Two basis functions defined on the face 456:
               !-------------------------------------------------
               TriangleFaceMap(:) = (/ 4,5,6 /)
               FaceIndices(1:3) = GIndexes(TriangleFaceMap(1:3))
               CALL TriangleFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               IF (RedefineFaceBasis) THEN
                 EdgeBasis(21,1:2) = (D1 * WorkBasis(I1,1:2) + D2 * WorkBasis(I2,1:2)) * 0.5d0 * h2
                 EdgeBasis(22,1:2) = (-D1 * WorkBasis(I1,1:2) + D2 * WorkBasis(I2,1:2)) * 0.5d0 * h2

                 CurlBasis(21,1) = -(D1 * WorkBasis(I1,2) + D2 * WorkBasis(I2,2)) * 0.5d0 * dh2
                 CurlBasis(21,2) = (D1 * WorkBasis(I1,1) + D2 * WorkBasis(I2,1)) * 0.5d0 * dh2
                 CurlBasis(21,3) = (D1 * WorkCurlBasis(I1,3) + D2 * WorkCurlBasis(I2,3)) * 0.5d0 * h2

                 CurlBasis(22,1) = -(-D1 * WorkBasis(I1,2) + D2 * WorkBasis(I2,2)) * 0.5d0 * dh2
                 CurlBasis(22,2) = (-D1 * WorkBasis(I1,1) + D2 * WorkBasis(I2,1)) * 0.5d0 * dh2
                 CurlBasis(22,3) = (-D1 * WorkCurlBasis(I1,3) + D2 * WorkCurlBasis(I2,3)) * 0.5d0 * h2

               ELSE
                 EdgeBasis(21,1:2) = D1 * WorkBasis(I1,1:2) * h2
                 CurlBasis(21,1) = -D1 * WorkBasis(I1,2) * dh2
                 CurlBasis(21,2) = D1 * WorkBasis(I1,1) * dh2
                 CurlBasis(21,3) = D1 * WorkCurlBasis(I1,3) * h2

                 EdgeBasis(22,1:2) = D2 * WorkBasis(I2,1:2) * h2
                 CurlBasis(22,1) = -D2 * WorkBasis(I2,2) * dh2
                 CurlBasis(22,2) = D2 * WorkBasis(I2,1) * dh2
                 CurlBasis(22,3) = D2 * WorkCurlBasis(I2,3) * h2
               END IF

               ! scale to reduce ill-conditioning:
               IF (ScaleFaceBasis) THEN
                 EdgeBasis(19:21:2,:) = sqrt(fs1) * EdgeBasis(19:21:2,:)
                 CurlBasis(19:21:2,:) = sqrt(fs1) * CurlBasis(19:21:2,:)
                 EdgeBasis(20:22:2,:) = sqrt(fs2) * EdgeBasis(20:22:2,:)
                 CurlBasis(20:22:2,:) = sqrt(fs2) * CurlBasis(20:22:2,:)
               END IF

               !-------------------------------------------------
               ! Four basis functions defined on the face 1254:
               !-------------------------------------------------              
               SquareFaceMap(:) = (/ 1,2,5,4 /)          
               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = (3.0d0 - Sqrt(3.0d0)*v)/6.0d0 * 4.0d0 * h3
               WorkBasis(1,2) = u/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
               WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
               WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
               WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
               WorkBasis(2,1) = -(u*(-3.0d0 + Sqrt(3.0d0)*v))/2.0d0 * 4.0d0 * h3
               WorkBasis(2,2) = (Sqrt(3.0d0)*u**2)/2.0d0 * 4.0d0 * h3
               WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
               WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
               WorkCurlBasis(2,3) = (3.0d0*Sqrt(3.0d0)*u)/2.0d0 * 4.0d0 * h3

               WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(1, u, v) * TriangleNodalPBasis(2, u, v)
               grad(1:2) = dTriangleNodalPBasis(1, u, v) * TriangleNodalPBasis(2, u, v) + &
                   TriangleNodalPBasis(1, u, v) * dTriangleNodalPBasis(2, u, v)
               WorkCurlBasis(3,1) = 2.0d0 * grad(2)
               WorkCurlBasis(3,2) = -2.0d0 * grad(1)
               WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
               WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
               WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(23,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(23,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(24,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(24,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(25,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(25,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(26,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(26,:) = WorkCurlBasis(2*(I2-1)+2,:)            

               !-------------------------------------------------
               ! Four basis functions defined on the face 2365:
               !-------------------------------------------------              
               SquareFaceMap(:) = (/ 2,3,6,5 /)          
               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
               WorkBasis(1,2) = (1 + u)/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
               WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
               WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
               WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
               WorkBasis(2,1) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*v)/4.0d0 * 4.0d0 * h3
               WorkBasis(2,2) = (Sqrt(3.0d0)*(1.0d0 + u)*(-1.0d0 - u + Sqrt(3.0d0)*v))/4.0d0 * 4.0d0 * h3
               WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
               WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
               WorkCurlBasis(2,3) = (-3.0d0*(Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v))/4.0d0 * 4.0d0 * h3

               WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(2, u, v) * TriangleNodalPBasis(3, u, v)
               grad(1:2) = dTriangleNodalPBasis(2, u, v) * TriangleNodalPBasis(3, u, v) + &
                   TriangleNodalPBasis(2, u, v) * dTriangleNodalPBasis(3, u, v)
               WorkCurlBasis(3,1) = 2.0d0 * grad(2)
               WorkCurlBasis(3,2) = -2.0d0 * grad(1)
               WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
               WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
               WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(27,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(27,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(28,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(28,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(29,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(29,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(30,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(30,:) = WorkCurlBasis(2*(I2-1)+2,:)  

               !-------------------------------------------------
               ! Four basis functions defined on the face 3146:
               !-------------------------------------------------              
               SquareFaceMap(:) = (/ 3,1,4,6 /)          
               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = -v/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
               WorkBasis(1,2) = (-1 + u)/(2.0d0*Sqrt(3.0d0)) * 4.0d0 * h3
               WorkCurlBasis(1,1) = -WorkBasis(1,2)/h3 * dh3 
               WorkCurlBasis(1,2) = WorkBasis(1,1)/h3 * dh3 
               WorkCurlBasis(1,3) = 1.0d0/Sqrt(3.0d0) * 4.0d0 * h3
               WorkBasis(2,1) = (v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0 * 4.0d0 * h3
               WorkBasis(2,2) =  -(Sqrt(3.0d0)*(-1.0d0 + u)*(-1.0d0 + u + Sqrt(3.0d0)*v))/4.0d0 * 4.0d0 * h3
               WorkCurlBasis(2,1) = -WorkBasis(2,2)/h3 * dh3 
               WorkCurlBasis(2,2) = WorkBasis(2,1)/h3 * dh3
               WorkCurlBasis(2,3) = (-3.0d0*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v))/4.0d0 * 4.0d0 * h3

               WorkBasis(3,3) = 2.0d0 * TriangleNodalPBasis(3, u, v) * TriangleNodalPBasis(1, u, v)
               grad(1:2) = dTriangleNodalPBasis(3, u, v) * TriangleNodalPBasis(1, u, v) + &
                   TriangleNodalPBasis(3, u, v) * dTriangleNodalPBasis(1, u, v)
               WorkCurlBasis(3,1) = 2.0d0 * grad(2)
               WorkCurlBasis(3,2) = -2.0d0 * grad(1)
               WorkBasis(4,3) = 3.0d0 * WorkBasis(3,3) * w
               WorkCurlBasis(4,1) = 6.0d0 * grad(2) * w
               WorkCurlBasis(4,2) = -6.0d0 * grad(1) * w

               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(31,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(31,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(32,:) = WorkBasis(2*(I1-1)+2,:)
               CurlBasis(32,:) = WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(33,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(33,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(34,:) = WorkBasis(2*(I2-1)+2,:)
               CurlBasis(34,:) = WorkCurlBasis(2*(I2-1)+2,:)  

               !-------------------------------------------------
               ! Two basis functions associated with the interior
               !-------------------------------------------------    
               EdgeBasis(35,1) = (v*(1.0d0 + u - v/Sqrt(3.0d0)))/(4.0d0*Sqrt(3.0d0)) * h3
               EdgeBasis(35,2) = ((-1.0d0 + u)*(-3.0d0 - 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
               CurlBasis(35,1) = -EdgeBasis(35,2)/h3 * dh3
               CurlBasis(35,2) = EdgeBasis(35,1)/h3 * dh3
               CurlBasis(35,3) = (-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0 * h3

               EdgeBasis(36,1) = (v*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
               EdgeBasis(36,2) = -((1.0d0 + u)*(-3.0d0 + 3.0d0*u + Sqrt(3.0d0)*v))/(12.0d0*Sqrt(3.0d0)) * h3
               CurlBasis(36,1) = -EdgeBasis(36,2)/h3 * dh3
               CurlBasis(36,2) = EdgeBasis(36,1)/h3 * dh3
               CurlBasis(36,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0 * h3

               IF (ScaleFaceBasis) THEN
                 EdgeBasis(35:36,1:2) = sqrt(150.0d0) * EdgeBasis(35:36,1:2)
                 CurlBasis(35:36,1:3) = sqrt(150.0d0) * CurlBasis(35:36,1:3)
               END IF
             END IF GRADIENT_VERSION_PRISM
           ELSE
             !--------------------------------------------------------------
             ! The lowest-order element from the optimal family. The optimal
             ! accuracy is obtained also for non-affine meshes.
             ! -------------------------------------------------------------
             ! First nine basis functions associated with the edges
             ! -------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = -((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + w)*w)/12.0d0
             EdgeBasis(1,2) = (u*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(1,3) = 0.0d0
             CurlBasis(1,1) = (u*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(1,2) = -((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + 2*w))/12.0d0
             CurlBasis(1,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(2,1) = -(v*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(2,2) = ((1.0d0 + u)*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0)) 
             EdgeBasis(2,3) = 0.0d0
             CurlBasis(2,1) = ((1.0d0 + u)*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(2,2) = (v*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(2,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(3,1) = -(v*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(3,2) = ((-1.0d0 + u)*(-1.0d0 + w)*w)/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(3,3) = 0.0d0
             CurlBasis(3,1) = ((-1.0d0 + u)*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(3,2) = (v*(1.0d0 - 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(3,3) = ((-1.0d0 + w)*w)/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(4,1) = -((-3.0d0 + Sqrt(3.0d0)*v)*w*(1.0d0 + w))/12.0d0
             EdgeBasis(4,2) = (u*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(4,3) = 0.0d0
             CurlBasis(4,1) = -(u*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(4,2) = -((-3.0d0 + Sqrt(3.0d0)*v)*(1.0d0 + 2.0d0*w))/12.0d0
             CurlBasis(4,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             EdgeBasis(5,1) = -(v*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(5,2) = ((1.0d0 + u)*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(5,3) = 0.0d0
             CurlBasis(5,1) = -((1.0d0 + u)*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(5,2) = -(v*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(5,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             EdgeBasis(6,1) = -(v*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(6,2) = ((-1.0d0 + u)*w*(1.0d0 + w))/(4.0d0*Sqrt(3.0d0))
             EdgeBasis(6,3) = 0.0d0
             CurlBasis(6,1) = -((-1.0d0 + u)*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(6,2) = -(v*(1.0d0 + 2.0d0*w))/(4.0d0*Sqrt(3.0d0))
             CurlBasis(6,3) = (w*(1.0d0 + w))/(2.0d0*Sqrt(3.0d0))
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             EdgeBasis(7,1) = 0.0d0
             EdgeBasis(7,2) = 0.0d0
             EdgeBasis(7,3) = (3*u**2 + v*(-Sqrt(3.0d0) + v) + u*(-3.0d0 + 2*Sqrt(3.0d0)*v))/12.0d0
             CurlBasis(7,1) = (-Sqrt(3.0d0) + 2*Sqrt(3.0d0)*u + 2*v)/12.0d0
             CurlBasis(7,2) = (3.0d0 - 6*u - 2*Sqrt(3.0d0)*v)/12.0d0
             CurlBasis(7,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             EdgeBasis(8,1) = 0.0d0
             EdgeBasis(8,2) = 0.0d0
             EdgeBasis(8,3) = (3*u**2 + v*(-Sqrt(3.0d0) + v) + u*(3.0d0 - 2*Sqrt(3.0d0)*v))/12.0d0
             CurlBasis(8,1) = (-Sqrt(3.0d0) - 2*Sqrt(3.0d0)*u + 2*v)/12.0d0
             CurlBasis(8,2) = (-3.0d0 - 6*u + 2*Sqrt(3.0d0)*v)/12.0d0 
             CurlBasis(8,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             i = EdgeMap(9,1)
             j = EdgeMap(9,2)
             EdgeBasis(9,1) = 0.0d0
             EdgeBasis(9,2) = 0.0d0
             EdgeBasis(9,3) = (v*(-Sqrt(3.0d0) + 2*v))/6.0d0
             CurlBasis(9,1) = (-Sqrt(3.0d0) + 4*v)/6.0d0
             CurlBasis(9,2) = 0.0d0
             CurlBasis(9,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(9,:) = -EdgeBasis(9,:)
               CurlBasis(9,:) = -CurlBasis(9,:)
             END IF

             ! ---------------------------------------------------------------------
             ! Additional six basis functions on the square faces (two per face).
             ! ---------------------------------------------------------------------         
             PrismSquareFaceMap(1,:) = (/ 1,2,5,4 /)
             PrismSquareFaceMap(2,:) = (/ 2,3,6,5 /)
             PrismSquareFaceMap(3,:) = (/ 3,1,4,6 /)

             ! The first square face:
             WorkBasis(1,1) = ((-3.0d0 + Sqrt(3.0d0)*v)*(-1.0d0 + w**2))/6.0d0
             WorkBasis(1,2) = -(u*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = (u*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,2) = (-1.0d0 + v/Sqrt(3.0d0))*w
             WorkCurlBasis(1,3) = -((-1.0d0 + w**2)/Sqrt(3.0d0)) 

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = (3.0d0 - 3*u**2 - 2*Sqrt(3.0d0)*v + v**2)/6.0d0
             WorkCurlBasis(2,1) = (-Sqrt(3.0d0) + v)/3.0d0
             WorkCurlBasis(2,2) = u
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(PrismSquareFaceMap(1,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(10,:) = D1 * WorkBasis(I1,:)
             CurlBasis(10,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(11,:) = D2 * WorkBasis(I2,:)
             CurlBasis(11,:) = D2 * WorkCurlBasis(I2,:) 

             ! The second square face:
             WorkBasis(1,1) = (v*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,2) = -((1.0d0 + u)*(-1.0d0 + w**2))/(2.0d0*Sqrt(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((1.0d0 + u)*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,2) = (v*w)/Sqrt(3.0d0)
             WorkCurlBasis(1,3) = -((-1.0d0 + w**2)/Sqrt(3.0d0))

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - v)*v)/3.0d0
             WorkCurlBasis(2,1) = (Sqrt(3.0d0) + Sqrt(3.0d0)*u - 2*v)/3.0d0
             WorkCurlBasis(2,2) = -(v/Sqrt(3.0d0))
             WorkCurlBasis(2,3) = 0.0d0 

             FaceIndices(1:4) = GIndexes(PrismSquareFaceMap(2,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(12,:) = D1 * WorkBasis(I1,:)
             CurlBasis(12,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(13,:) = D2 * WorkBasis(I2,:)
             CurlBasis(13,:) = D2 * WorkCurlBasis(I2,:) 

             ! The third square face:
             WorkBasis(1,1) = (v*(-1.0d0 + w**2))/(2.0d0*SQRT(3.0d0))
             WorkBasis(1,2) = -((-1.0d0 + u)*(-1.0d0 + w**2))/(2.0d0*SQRT(3.0d0))
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((-1.0d0 + u)*w)/SQRT(3.0d0)
             WorkCurlBasis(1,2) = (v*w)/SQRT(3.0d0)
             WorkCurlBasis(1,3) = -(-1.0d0 + w**2)/SQRT(3.0d0)

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -(v*(-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v))/3.0d0
             WorkCurlBasis(2,1) = (Sqrt(3.0d0) - Sqrt(3.0d0)*u - 2*v)/3.0d0
             WorkCurlBasis(2,2) = v/Sqrt(3.0d0)
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(PrismSquareFaceMap(3,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(14,:) = D1 * WorkBasis(I1,:)
             CurlBasis(14,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(15,:) = D2 * WorkBasis(I2,:)
             CurlBasis(15,:) = D2 * WorkCurlBasis(I2,:) 
           END IF

         CASE(8)
           !--------------------------------------------------------------
           ! This branch is for handling brick elements
           !--------------------------------------------------------------           
           EdgeMap => GetEdgeMap(8)
           
           IF (SecondOrder) THEN
             !---------------------------------------------------------------
             ! The second-order element from the Nedelec's first family 
             ! (note that the lowest-order brick element is from a different 
             ! family). This element may not be optimally accurate if 
             ! the physical element is not affine.
             !--------------------------------------------------------------             
  
             ! Edges 12 and 43 ...
             DO q=1,2
               k = 2*q-1 ! Edge number k: 1 ~ 12 and 3 ~ 43 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(1,w) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = 0.5d0 * (-0.5d0) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,3) = -0.5d0 * LineNodalPBasis(1,w) * dLineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(1,w) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = 1.5d0 * (-0.5d0) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,3) = -1.5d0 * LineNodalPBasis(1,w) * u * dLineNodalPBasis(q,v)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO

             ! Edges 56 and 87 ...
             DO q=1,2
               k = 4 + 2*q-1 ! Edge number k: 5 ~ 56 and 7 ~ 87 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(2,w) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = 0.5d0 * 0.5d0 * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,3) = -0.5d0 * LineNodalPBasis(2,w) * dLineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(2,w) * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = 1.5d0 * 0.5d0 * u * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,3) = -1.5d0 * LineNodalPBasis(2,w) * u * dLineNodalPBasis(q,v)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO

             ! Edges 23 and 14 ...
             DO q=1,2
               k = 2*q ! Edge number k: 2 ~ 23 and 4 ~ 14 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,2) = 0.5d0 * LineNodalPBasis(1,w) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,1) = -0.5d0 * (-0.5d0) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(1,w) * dLineNodalPBasis(3-q,u)
               EdgeBasis(2*(k-1)+2,2) = 1.5d0 * LineNodalPBasis(1,w) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,1) = -1.5d0 * (-0.5d0) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(1,w) * v * dLineNodalPBasis(3-q,u)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO            

             ! Edges 67 and 58 ...
             DO q=1,2
               k = 4+2*q ! Edge number k: 6 ~ 67 and 8 ~ 58 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,2) = 0.5d0 * LineNodalPBasis(2,w) * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,1) = -0.5d0 * 0.5d0 * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(2,w) * dLineNodalPBasis(3-q,u)
               EdgeBasis(2*(k-1)+2,2) = 1.5d0 * LineNodalPBasis(2,w) * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,1) = -1.5d0 * 0.5d0 * v * LineNodalPBasis(3-q,u) 
               CurlBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(2,w) * v * dLineNodalPBasis(3-q,u)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO          

             ! Edges 15 and 48 ...
             DO q=1,2
               k = 8+3*(q-1)+1 ! Edge number k: 9 ~ 15 and 12 ~ 48 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(1,u) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(1,u) * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = -0.5d0 * dLineNodalPBasis(1,u) * LineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(1,u) * w * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(1,u) * w * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = -1.5d0 * dLineNodalPBasis(1,u) * w * LineNodalPBasis(q,v)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO         

             ! Edges 26 and 37 ...
             DO q=1,2
               k = 9+q ! Edge number k: 10 ~ 26 and 11 ~ 37 
               i = EdgeMap(k,1)
               j = EdgeMap(k,2)
               EdgeBasis(2*(k-1)+1,3) = 0.5d0 * LineNodalPBasis(2,u) * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,1) = 0.5d0 * LineNodalPBasis(2,u) * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+1,2) = -0.5d0 * dLineNodalPBasis(2,u) * LineNodalPBasis(q,v)
               EdgeBasis(2*(k-1)+2,3) = 1.5d0 * LineNodalPBasis(2,u) * w * LineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,1) = 1.5d0 * LineNodalPBasis(2,u) * w * dLineNodalPBasis(q,v) 
               CurlBasis(2*(k-1)+2,2) = -1.5d0 * dLineNodalPBasis(2,u) * w * LineNodalPBasis(q,v)
               IF(GIndexes(j)<GIndexes(i)) THEN
                 EdgeBasis(2*(k-1)+1,:) = -EdgeBasis(2*(k-1)+1,:)
                 CurlBasis(2*(k-1)+1,:) = -CurlBasis(2*(k-1)+1,:)
               END IF
             END DO     

             ! ---------------------------------------------------------------------
             ! Additional basis functions on the square faces (four per face).
             ! ---------------------------------------------------------------------         

             ! Faces 1234 and 5678:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,2,3,4 /)
               CASE(2)
                 SquareFaceMap(:) = (/ 5,6,7,8 /)
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * LineNodalPBasis(q,w)
               WorkCurlBasis(1,2) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * dLineNodalPBasis(q,w)
               WorkCurlBasis(1,3) = v * LineNodalPBasis(q,w)

               WorkBasis(2,1) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * u * LineNodalPBasis(q,w)
               WorkCurlBasis(2,2) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * u * dLineNodalPBasis(q,w)
               WorkCurlBasis(2,3) = -12.0d0 * (-0.5d0 * v) * u * dLineNodalPBasis(q,w) 

               WorkBasis(3,2) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * LineNodalPBasis(q,w)
               WorkCurlBasis(3,1) = -2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * dLineNodalPBasis(q,w)
               WorkCurlBasis(3,3) = -u * LineNodalPBasis(q,w)
               
               WorkBasis(4,2) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * v * LineNodalPBasis(q,w)
               WorkCurlBasis(4,1) = -12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * v * dLineNodalPBasis(q,w)
               WorkCurlBasis(4,3) = 12.0d0 * (-0.5d0 * u) * v * LineNodalPBasis(q,w)
               
               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               k = 24
               EdgeBasis(k+4*(q-1)+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+4*(q-1)+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+4*(q-1)+2,:) = 0.5d0 * WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+4*(q-1)+2,:) = 0.5d0 * WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+4*(q-1)+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+4*(q-1)+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4*(q-1)+4,:) = 0.5d0 * WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4*(q-1)+4,:) = 0.5d0 * WorkCurlBasis(2*(I2-1)+2,:)
             END DO

             ! Faces 1265 and 4378:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,2,6,5 /)
                 k = 32
               CASE(2)
                 SquareFaceMap(:) = (/ 4,3,7,8 /)
                 k = 40
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,1) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * LineNodalPBasis(q,v)
               WorkCurlBasis(1,2) = 2.0d0 * (-0.5d0 * w) * LineNodalPBasis(q,v)
               WorkCurlBasis(1,3) = -2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * dLineNodalPBasis(q,v)

               WorkBasis(2,1) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * LineNodalPBasis(q,v)
               WorkCurlBasis(2,2) = 12.0d0 * (-0.5d0 * w) * u * LineNodalPBasis(q,v)
               WorkCurlBasis(2,3) = -12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * dLineNodalPBasis(q,v)

               WorkBasis(3,3) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * LineNodalPBasis(q,v)
               WorkCurlBasis(3,1) = 2.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * dLineNodalPBasis(q,v)
               WorkCurlBasis(3,2) = u * LineNodalPBasis(q,v)
               
               WorkBasis(4,3) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * w * LineNodalPBasis(q,v)
               WorkCurlBasis(4,1) = 12.0d0 * LineNodalPBasis(1,u) * LineNodalPBasis(2,u) * w * dLineNodalPBasis(q,v)
               WorkCurlBasis(4,2) = -12.0d0 * (-0.5d0 * u) * w * LineNodalPBasis(q,v)
               
               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(k+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+2,:) = 0.5d0 * WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+2,:) = 0.5d0 * WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4,:) = 0.5d0 * WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4,:) = 0.5d0 * WorkCurlBasis(2*(I2-1)+2,:)
             END DO
             
             ! Faces 2376 and 1485:
             DO q=1,2
               SELECT CASE(q)
               CASE(1)
                 SquareFaceMap(:) = (/ 1,4,8,5 /)
                 k = 44
               CASE(2)
                 SquareFaceMap(:) = (/ 2,3,7,6 /)
                 k = 36
               END SELECT

               WorkBasis = 0.0d0
               WorkCurlBasis = 0.0d0

               WorkBasis(1,2) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * LineNodalPBasis(q,u)
               WorkCurlBasis(1,1) = -2.0d0 * (-0.5d0 * w) * LineNodalPBasis(q,u)
               WorkCurlBasis(1,3) = 2.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * dLineNodalPBasis(q,u)

               WorkBasis(2,2) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * LineNodalPBasis(q,u)
               WorkCurlBasis(2,1) = -12.0d0 * (-0.5d0 * w) * v * LineNodalPBasis(q,u)
               WorkCurlBasis(2,3) = 12.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * dLineNodalPBasis(q,u)

               WorkBasis(3,3) = 2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * LineNodalPBasis(q,u)
               WorkCurlBasis(3,1) = 2.0d0 * (-0.5d0 * v) * LineNodalPBasis(q,u)
               WorkCurlBasis(3,2) = -2.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * dLineNodalPBasis(q,u)
               
               WorkBasis(4,3) = 12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * LineNodalPBasis(q,u)
               WorkCurlBasis(4,1) = 12.0d0 * (-0.5d0 * v) * w * LineNodalPBasis(q,u)
               WorkCurlBasis(4,2) = -12.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * dLineNodalPBasis(q,u)
               
               FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
               CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

               EdgeBasis(k+1,:) = D1 * WorkBasis(2*(I1-1)+1,:)
               CurlBasis(k+1,:) = D1 * WorkCurlBasis(2*(I1-1)+1,:)
               EdgeBasis(k+2,:) = 0.5d0 * WorkBasis(2*(I1-1)+2,:)
               CurlBasis(k+2,:) = 0.5d0 * WorkCurlBasis(2*(I1-1)+2,:)
               EdgeBasis(k+3,:) = D2 * WorkBasis(2*(I2-1)+1,:)
               CurlBasis(k+3,:) = D2 * WorkCurlBasis(2*(I2-1)+1,:)
               EdgeBasis(k+4,:) = 0.5d0 * WorkBasis(2*(I2-1)+2,:)
               CurlBasis(k+4,:) = 0.5d0 * WorkCurlBasis(2*(I2-1)+2,:)
             END DO

             ! Interior basis functions, two per coordinate direction:

             EdgeBasis(49,1) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * &
                 LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(49,2) = 8.0d0 * (-0.5d0 * w) * LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(49,3) = -8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * (-0.5d0 * v)

             EdgeBasis(50,1) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u * &
                 LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(50,2) = 24.0d0 * (-0.5d0 * w) * u * LineNodalPBasis(1,v) * LineNodalPBasis(2,v)
             CurlBasis(50,3) = -24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * u *  (-0.5d0 * v)

 
             EdgeBasis(51,2) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(51,1) = -8.0d0 * (-0.5d0 * w) * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(51,3) = 8.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * (-0.5d0 * u)

             EdgeBasis(52,2) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(52,1) = -24.0d0 * (-0.5d0 * w) * v * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(52,3) = 24.0d0 * LineNodalPBasis(1,w) * LineNodalPBasis(2,w) * v * (-0.5d0 * u)
            
             EdgeBasis(53,3) = 8.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(53,1) = 8.0d0 * (-0.5d0 * v) * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(53,2) = -8.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * (-0.5d0 * u)

             EdgeBasis(54,3) = 24.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * &
                 LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(54,1) = 24.0d0 * (-0.5d0 * v) * w * LineNodalPBasis(1,u) * LineNodalPBasis(2,u)
             CurlBasis(54,2) = -24.0d0 * LineNodalPBasis(1,v) * LineNodalPBasis(2,v) * w * (-0.5d0 * u)

           ELSE
             !--------------------------------------------------------------
             ! The lowest-order element from the optimal family. The optimal
             ! accuracy is obtained also for non-affine meshes.
             ! -------------------------------------------------------------
             ! First twelve basis functions associated with the edges
             ! -------------------------------------------------------------
             i = EdgeMap(1,1)
             j = EdgeMap(1,2)
             EdgeBasis(1,1) = ((-1.0d0 + v)*v*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(1,2) = 0.0d0
             EdgeBasis(1,3) = 0.0d0
             CurlBasis(1,1) = 0.0d0
             CurlBasis(1,2) = ((-1.0d0 + v)*v*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(1,3) = -((-1.0d0 + 2*v)*(-1.0d0 + w)*w)/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(1,:) = -EdgeBasis(1,:)
               CurlBasis(1,:) = -CurlBasis(1,:)
             END IF

             i = EdgeMap(2,1)
             j = EdgeMap(2,2)
             EdgeBasis(2,1) = 0.0d0
             EdgeBasis(2,2) = (u*(1.0d0 + u)*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(2,3) = 0.0d0
             CurlBasis(2,1) = -(u*(1.0d0 + u)*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(2,2) = 0.0d0
             CurlBasis(2,3) = ((1.0d0 + 2*u)*(-1.0d0 + w)*w)/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(2,:) = -EdgeBasis(2,:)
               CurlBasis(2,:) = -CurlBasis(2,:)
             END IF

             i = EdgeMap(3,1)
             j = EdgeMap(3,2)
             EdgeBasis(3,1) = (v*(1.0d0 + v)*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(3,2) = 0.0d0
             EdgeBasis(3,3) = 0.0d0
             CurlBasis(3,1) = 0.0d0
             CurlBasis(3,2) = (v*(1.0d0 + v)*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(3,3) = -((1.0d0 + 2*v)*(-1.0d0 + w)*w)/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(3,:) = -EdgeBasis(3,:)
               CurlBasis(3,:) = -CurlBasis(3,:)
             END IF

             i = EdgeMap(4,1)
             j = EdgeMap(4,2)
             EdgeBasis(4,1) = 0.0d0
             EdgeBasis(4,2) = ((-1.0d0 + u)*u*(-1.0d0 + w)*w)/8.0d0
             EdgeBasis(4,3) = 0.0d0
             CurlBasis(4,1) = -((-1.0d0 + u)*u*(-1.0d0 + 2*w))/8.0d0
             CurlBasis(4,2) = 0.0d0
             CurlBasis(4,3) = ((-1.0d0 + 2*u)*(-1.0d0 + w)*w)/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(4,:) = -EdgeBasis(4,:)
               CurlBasis(4,:) = -CurlBasis(4,:)
             END IF

             i = EdgeMap(5,1)
             j = EdgeMap(5,2)
             EdgeBasis(5,1) = ((-1.0d0 + v)*v*w*(1.0d0 + w))/8.0d0
             EdgeBasis(5,2) = 0.0d0
             EdgeBasis(5,3) = 0.0d0
             CurlBasis(5,1) = 0.0d0
             CurlBasis(5,2) = ((-1.0d0 + v)*v*(1.0d0 + 2*w))/8.0d0 
             CurlBasis(5,3) = -((-1.0d0 + 2*v)*w*(1.0d0 + w))/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(5,:) = -EdgeBasis(5,:)
               CurlBasis(5,:) = -CurlBasis(5,:)
             END IF

             i = EdgeMap(6,1)
             j = EdgeMap(6,2)
             EdgeBasis(6,1) = 0.0d0
             EdgeBasis(6,2) = (u*(1.0d0 + u)*w*(1.0d0 + w))/8.0d0
             EdgeBasis(6,3) = 0.0d0
             CurlBasis(6,1) = -(u*(1.0d0 + u)*(1.0d0 + 2*w))/8.0d0
             CurlBasis(6,2) = 0.0d0
             CurlBasis(6,3) = ((1.0d0 + 2*u)*w*(1.0d0 + w))/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(6,:) = -EdgeBasis(6,:)
               CurlBasis(6,:) = -CurlBasis(6,:)
             END IF

             i = EdgeMap(7,1)
             j = EdgeMap(7,2)
             EdgeBasis(7,1) = (v*(1.0d0 + v)*w*(1.0d0 + w))/8.0d0
             EdgeBasis(7,2) = 0.0d0
             EdgeBasis(7,3) = 0.0d0
             CurlBasis(7,1) = 0.0d0
             CurlBasis(7,2) = (v*(1.0d0 + v)*(1.0d0 + 2*w))/8.0d0
             CurlBasis(7,3) = -((1.0d0 + 2*v)*w*(1.0d0 + w))/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(7,:) = -EdgeBasis(7,:)
               CurlBasis(7,:) = -CurlBasis(7,:)
             END IF

             i = EdgeMap(8,1)
             j = EdgeMap(8,2)
             EdgeBasis(8,1) = 0.0d0
             EdgeBasis(8,2) = ((-1.0d0 + u)*u*w*(1.0d0 + w))/8.0d0
             EdgeBasis(8,3) = 0.0d0
             CurlBasis(8,1) = -((-1.0d0 + u)*u*(1.0d0 + 2*w))/8.0d0
             CurlBasis(8,2) = 0.0d0
             CurlBasis(8,3) = ((-1.0d0 + 2*u)*w*(1.0d0 + w))/8.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(8,:) = -EdgeBasis(8,:)
               CurlBasis(8,:) = -CurlBasis(8,:)
             END IF

             i = EdgeMap(9,1)
             j = EdgeMap(9,2)
             EdgeBasis(9,1) = 0.0d0
             EdgeBasis(9,2) = 0.0d0
             EdgeBasis(9,3) = ((-1.0d0 + u)*u*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(9,1) = ((-1.0d0 + u)*u*(-1.0d0 + 2*v))/8.0d0
             CurlBasis(9,2) = -((-1.0d0 + 2*u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(9,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(9,:) = -EdgeBasis(9,:)
               CurlBasis(9,:) = -CurlBasis(9,:)
             END IF

             i = EdgeMap(10,1)
             j = EdgeMap(10,2)
             EdgeBasis(10,1) = 0.0d0
             EdgeBasis(10,2) = 0.0d0
             EdgeBasis(10,3) = (u*(1.0d0 + u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(10,1) = (u*(1.0d0 + u)*(-1.0d0 + 2*v))/8.0d0
             CurlBasis(10,2) = -((1.0d0 + 2*u)*(-1.0d0 + v)*v)/8.0d0
             CurlBasis(10,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(10,:) = -EdgeBasis(10,:)
               CurlBasis(10,:) = -CurlBasis(10,:)
             END IF

             i = EdgeMap(11,1)
             j = EdgeMap(11,2)
             EdgeBasis(11,1) = 0.0d0
             EdgeBasis(11,2) = 0.0d0
             EdgeBasis(11,3) = (u*(1.0d0 + u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(11,1) = (u*(1.0d0 + u)*(1.0d0 + 2*v))/8.0d0
             CurlBasis(11,2) = -((1.0d0 + 2*u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(11,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(11,:) = -EdgeBasis(11,:)
               CurlBasis(11,:) = -CurlBasis(11,:)
             END IF

             i = EdgeMap(12,1)
             j = EdgeMap(12,2)
             EdgeBasis(12,1) = 0.0d0
             EdgeBasis(12,2) = 0.0d0
             EdgeBasis(12,3) = ((-1.0d0 + u)*u*v*(1.0d0 + v))/8.0d0
             CurlBasis(12,1) = ((-1.0d0 + u)*u*(1.0d0 + 2*v))/8.0d0
             CurlBasis(12,2) = -((-1.0d0 + 2*u)*v*(1.0d0 + v))/8.0d0
             CurlBasis(12,3) = 0.0d0
             IF(GIndexes(j)<GIndexes(i)) THEN
               EdgeBasis(12,:) = -EdgeBasis(12,:)
               CurlBasis(12,:) = -CurlBasis(12,:)
             END IF

             ! ---------------------------------------------------------------------
             ! Additional twelve basis functions on the square faces (two per face).
             ! ---------------------------------------------------------------------         
             BrickFaceMap(1,:) = (/ 1,2,3,4 /)          
             BrickFaceMap(2,:) = (/ 5,6,7,8 /)
             BrickFaceMap(3,:) = (/ 1,2,6,5 /)
             BrickFaceMap(4,:) = (/ 2,3,7,6 /)
             BrickFaceMap(5,:) = (/ 4,3,7,8 /)
             BrickFaceMap(6,:) = (/ 1,4,8,5 /)

             ! The first face:
             WorkBasis(1,1) = -((-1.0d0 + v**2)*(-1.0d0 + w)*w)/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v**2)*(-1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(1,3) = (v*(-1.0d0 + w)*w)/2.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = -((-1.0d0 + u**2)*(-1.0d0 + w)*w)/4.0d0
             WorkBasis(2,3) = 0.0d0
             WorkCurlBasis(2,1) = ((-1.0d0 + u**2)*(-1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(2,2) = 0.0d0
             WorkCurlBasis(2,3) = -(u*(-1.0d0 + w)*w)/2.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(1,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(13,:) = D1 * WorkBasis(I1,:)
             CurlBasis(13,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(14,:) = D2 * WorkBasis(I2,:)
             CurlBasis(14,:) = D2 * WorkCurlBasis(I2,:) 

             ! The second face:
             WorkBasis(1,1) = -((-1.0d0 + v**2)*w*(1.0d0 + w))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v**2)*(1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(1,3) = (v*w*(1.0d0 + w))/2.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = -((-1.0d0 + u**2)*w*(1.0d0 + w))/4.0d0
             WorkBasis(2,3) = 0.0d0
             WorkCurlBasis(2,1) = ((-1.0d0 + u**2)*(1.0d0 + 2*w))/4.0d0
             WorkCurlBasis(2,2) = 0.0d0
             WorkCurlBasis(2,3) = -(u*w*(1.0d0 + w))/2.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(2,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(15,:) = D1 * WorkBasis(I1,:)
             CurlBasis(15,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(16,:) = D2 * WorkBasis(I2,:)
             CurlBasis(16,:) = D2 * WorkCurlBasis(I2,:) 

             ! The third face:
             WorkBasis(1,1) = -((-1.0d0 + v)*v*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -((-1.0d0 + v)*v*w)/2.0d0
             WorkCurlBasis(1,3) = ((-1.0d0 + 2*v)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u**2)*(-1.0d0 + v)*v)/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u**2)*(-1.0d0 + 2*v))/4.0d0
             WorkCurlBasis(2,2) = (u*(-1.0d0 + v)*v)/2.0d0
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(3,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(17,:) = D1 * WorkBasis(I1,:)
             CurlBasis(17,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(18,:) = D2 * WorkBasis(I2,:)
             CurlBasis(18,:) = D2 * WorkCurlBasis(I2,:) 

             ! The fourth face:
             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = -(u*(1.0d0 + u)*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = (u*(1.0d0 + u)*w)/2.0d0
             WorkCurlBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = -((1.0d0 + 2*u)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -(u*(1.0d0 + u)*(-1 + v**2))/4.0d0
             WorkCurlBasis(2,1) = -(u*(1.0d0 + u)*v)/2.0d0
             WorkCurlBasis(2,2) = ((1.0d0 + 2*u)*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(4,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(19,:) = D1 * WorkBasis(I1,:)
             CurlBasis(19,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(20,:) = D2 * WorkBasis(I2,:)
             CurlBasis(20,:) = D2 * WorkCurlBasis(I2,:) 

             ! The fifth face:
             WorkBasis(1,1) = -(v*(1.0d0 + v)*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,2) = 0.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = 0.0d0
             WorkCurlBasis(1,2) = -(v*(1.0d0 + v)*w)/2.0d0
             WorkCurlBasis(1,3) = ((1.0d0 + 2*v)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u**2)*v*(1.0d0 + v))/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u**2)*(1.0d0 + 2*v))/4.0d0
             WorkCurlBasis(2,2) = (u*v*(1.0d0 + v))/2.0d0
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(5,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(21,:) = D1 * WorkBasis(I1,:)
             CurlBasis(21,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(22,:) = D2 * WorkBasis(I2,:)
             CurlBasis(22,:) = D2 * WorkCurlBasis(I2,:) 

             ! The sixth face:
             WorkBasis(1,1) = 0.0d0
             WorkBasis(1,2) = -((-1.0d0 + u)*u*(-1.0d0 + w**2))/4.0d0
             WorkBasis(1,3) = 0.0d0
             WorkCurlBasis(1,1) = ((-1.0d0 + u)*u*w)/2.0d0
             WorkCurlBasis(1,2) = 0.0d0
             WorkCurlBasis(1,3) = -((-1.0d0 + 2*u)*(-1.0d0 + w**2))/4.0d0

             WorkBasis(2,1) = 0.0d0
             WorkBasis(2,2) = 0.0d0
             WorkBasis(2,3) = -((-1.0d0 + u)*u*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,1) = -((-1.0d0 + u)*u*v)/2.0d0
             WorkCurlBasis(2,2) = ((-1.0d0 + 2*u)*(-1.0d0 + v**2))/4.0d0
             WorkCurlBasis(2,3) = 0.0d0

             FaceIndices(1:4) = GIndexes(BrickFaceMap(6,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)

             EdgeBasis(23,:) = D1 * WorkBasis(I1,:)
             CurlBasis(23,:) = D1 * WorkCurlBasis(I1,:)
             EdgeBasis(24,:) = D2 * WorkBasis(I2,:)
             CurlBasis(24,:) = D2 * WorkCurlBasis(I2,:) 

             ! ------------------------------------------------------------------------
             ! Additional basis functions on the element interior (three per element)
             ! -----------------------------------------------------------------------
             EdgeBasis(25,1) = ((-1.0d0 + v**2)*(-1.0d0 + w**2))/2.0d0
             EdgeBasis(25,2) = 0.0d0
             EdgeBasis(25,3) = 0.0d0
             CurlBasis(25,1) = 0.0d0
             CurlBasis(25,2) = (-1.0d0 + v**2)*w
             CurlBasis(25,3) = v - v*w**2

             EdgeBasis(26,1) = 0.0d0
             EdgeBasis(26,2) = ((-1.0d0 + u**2)*(-1.0d0 + w**2))/2.0d0
             EdgeBasis(26,3) = 0.0d0
             CurlBasis(26,1) = w - u**2*w
             CurlBasis(26,2) = 0.0d0
             CurlBasis(26,3) = u*(-1 + w**2)

             EdgeBasis(27,1) = 0.0d0
             EdgeBasis(27,2) = 0.0d0
             EdgeBasis(27,3) = ((-1.0d0 + u**2)*(-1.0d0 + v**2))/2.0d0
             CurlBasis(27,1) = (-1.0d0 + u**2)*v
             CurlBasis(27,2) = u - u*v**2
             CurlBasis(27,3) = 0.0d0
           END IF

         CASE DEFAULT
           CALL Fatal('ElementDescription::EdgeElementInfo','Unsupported element type')
         END SELECT
       END IF

       IF (cdim == dim) THEN
          !--------------------------------------------------------------------------------
          ! To optimize computation, this branch avoids calling the ElementMetric function
          ! since all necessary data has already been found via PiolaTransformationData.
          !-------------------------------------------------------------------------------
          IF (PerformPiolaTransform) THEN
             DO j=1,DOFs
                DO k=1,dim
                   B(k) = SUM( LG(k,1:dim) * EdgeBasis(j,1:dim) )
                END DO
                EdgeBasis(j,1:dim) = B(1:dim)

                IF (dim == 2) THEN
                   CurlBasis(j,3) = 1.0d0/DetF * CurlBasis(j,3)
                ELSE
                   DO k=1,dim
                      B(k) = 1.0d0/DetF * SUM( LF(k,1:dim) * CurlBasis(j,1:dim) )
                   END DO
                   CurlBasis(j,1:dim) = B(1:dim)
                END IF
             END DO
             ! Make the returned value DetF to act as a metric term for integration
             ! over the volume of the element: 
             DetF = ABS(DetF)
          END IF

          ! ----------------------------------------------------------------------
          ! Get global first derivatives of the nodal basis functions if wanted:
          ! ----------------------------------------------------------------------
          IF ( PRESENT(dBasisdx) ) THEN
             dBasisdx = 0.0d0
             DO i=1,n
                DO j=1,dim
                   DO k=1,dim
                      dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LG(j,k)
                   END DO
                END DO
             END DO
          END IF
       ELSE
          ! ----------------------------------------------------------------------
          ! We should enter this branch in the case of 2-D elements (dim=2) 
          ! embedded in the three-dimensional space (cdim=3). The following function
          ! defines LG to be the transpose of the pseudoinverse of F = LF.
          ! ----------------------------------------------------------------------       
          IF (PerformPiolaTransform .OR. PRESENT(dBasisdx) .OR. ApplyTraceMapping) THEN
             IF ( .NOT. ElementMetric( n, Element, Nodes, &
                  ElmMetric, detJ, dLBasisdx, LG ) ) THEN
                stat = .FALSE.
                RETURN
             END IF
          END IF

          IF (ApplyTraceMapping) THEN
            ! Perform operation b -> b x n. The resulting field transforms under the usual 
            ! Piola transform (like div-conforming field). For a general surface element
            ! embedded in 3D we return B(f(p))=1/sqrt(a) F(b x n) where a is the determinant of
            ! the metric tensor, F=[a1 a2] with a1 and a2 surface basis vectors and (b x n) is
            ! considered to be 2-vector (the trivial component ignored). Note that asking simultaneously 
            ! for the curl of the basis is not an expected combination.
            DO j=1,DOFs
              WorkBasis(1,1:2) = EdgeBasis(j,1:2)
              EdgeBasis(j,1) = WorkBasis(1,2)
              EdgeBasis(j,2) = -WorkBasis(1,1)
            END DO
            IF (PerformPiolaTransform) THEN
              DO j=1,DOFs 
                DO k=1,cdim
                  B(k) = SUM( LF(k,1:dim) * EdgeBasis(j,1:dim) ) / DetJ
                END DO
                EdgeBasis(j,1:cdim) = B(1:cdim)                
              END DO
            END IF
          ELSE
            IF (PerformPiolaTransform) THEN
              DO j=1,DOFs
                DO k=1,cdim
                  B(k) = SUM( LG(k,1:dim) * EdgeBasis(j,1:dim) )
                END DO
                EdgeBasis(j,1:cdim) = B(1:cdim)
                ! The returned spatial curl in the case cdim=3 and dim=2 handled here
                ! has limited usability. This handles only either a transformation of
                ! the type x_3 = p_3 or the normal component of curl for an arbitrarily
                ! oriented surface element. Note that the normal component is returned
                ! as the third entry, so this value has to be multiplied with the normal
                ! vector to get the vector representation of the normal component with
                ! respect to the coordinate axes. 
                CurlBasis(j,3) = 1.0d0/DetJ * CurlBasis(j,3)
              END DO
            END IF
          END IF

          ! Make the returned value DetF to act as a metric term for integration
          ! over the volume of the element: 
          DetF = DetJ

          ! ----------------------------------------------------------------------
          ! Get global first derivatives of the nodal basis functions if wanted:
          ! ----------------------------------------------------------------------
          IF ( PRESENT(dBasisdx) ) THEN
             dBasisdx = 0.0d0
             DO i=1,n
                DO j=1,cdim
                   DO k=1,dim
                      dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LG(j,k)
                   END DO
                END DO
             END DO
          END IF

       END IF

       IF(PRESENT(F)) F = LF
       IF(PRESENT(G)) G = LG
       IF(PRESENT(RotBasis)) RotBasis(1:DOFs,:) = CurlBasis(1:DOFs,:)
!-----------------------------------------------------------------------------
     END FUNCTION EdgeElementInfo
!------------------------------------------------------------------------------

     
     
!----------------------------------------------------------------------------
     SUBROUTINE TriangleFaceDofsOrdering(I1,I2,D1,D2,Ind,A,B,C)       
!-----------------------------------------------------------------------------
! This is used for selecting what additional basis functions are associated
! with a triangular face in the case of second-order approximation in H(curl).
! Given a triangular face [ijk] this subroutine can be used to pick two basis
! functions from an array of three candidate functions (for Nedelec's first family)
!
!    b_1 = L_k W_{ij}
!    b_2 = L_j W_{ik}
!    b_3 = L_i W_{jk}       
!
! such that the two basis functions are L_C W_{AB} and L_B W_{AC}. Here W_{ij}
! denotes the Whitney form and {A,B,C} are the global node indices such that
! A < B < C. D1 and D2 indicate whether sign reversions must be applied to
! the pre-tabulated basis functions. The indices corresponding to A, B and C
! may also be returned.      
! ----------------------------------------------------------------------------
       INTEGER, INTENT(OUT) :: I1, I2
       REAL(KIND=dp), INTENT(OUT) :: D1, D2
       INTEGER, INTENT(IN) :: Ind(4)
       INTEGER, OPTIONAL, INTENT(OUT) :: A, B, C
!---------------------------------------------------------------------------
       INTEGER ::  i, j, k
! --------------------------------------------------------------------------
       D1 = 1.0d0
       D2 = 1.0d0
       IF ( Ind(1) < Ind(2) ) THEN
          i = 1
       ELSE
          i = 2
       END IF
       IF ( Ind(i) > Ind(3) ) THEN
          i = 3
       END IF

       SELECT CASE(i)
       CASE(1)
          IF (Ind(3) > Ind(2)) THEN
             j = 2
             k = 3
             I1 = 1
             I2 = 2
          ELSE
             j = 3
             k = 2
             I1 = 2
             I2 = 1             
          END IF
       CASE(2)
         IF (Ind(3) > Ind(1)) THEN
             j = 1
             k = 3
             I1 = 1
             I2 = 3
             D1 = -1.0d0
          ELSE
             j = 3
             k = 1
             I1 = 3
             I2 = 1
             D2 = -1.0d0             
          END IF
       CASE(3)
          IF (Ind(2) > Ind(1)) THEN
             j = 1
             k = 2
             I1 = 2
             I2 = 3
          ELSE
             j = 2
             k = 1
             I1 = 3
             I2 = 2
          END IF
          D1 = -1.0d0
          D2 = -1.0d0          
       CASE DEFAULT
          CALL Fatal('ElementDescription::TriangleFaceDofsOrdering','Erratic triangular face Indices')
       END SELECT
       IF (PRESENT(A)) A = i
       IF (PRESENT(B)) B = j
       IF (PRESENT(C)) C = k
!---------------------------------------------------------
     END SUBROUTINE TriangleFaceDofsOrdering
!-----------------------------------------------------------

!----------------------------------------------------------------------------
     SUBROUTINE TriangleFaceDofsOrdering2nd(I1,I2,I3,Ind)       
!-----------------------------------------------------------------------------
! This is used for selecting the order of additional basis functions associated
! with a triangular face in the case of a higher-order approximation in H(curl) when
! the Nedelec second family is used. Given a triangular face [ijk] this subroutine
! can be used to permute an array of three candidate basis functions
!
!    b_1 = L_j L_k grad L_i
!    b_2 = L_i L_k grad L_j
!    b_3 = L_i L_j grad L_k
!
! such that the basis functions are  {L_B L_C grad L_A, L_A L_C grad L_B,
! L_A L_B grad L_C}. Here {A,B,C} are the global node indices such that
! A < B < C. These indices are returned as I1 = A, I2 = B and I3 = C.
! ----------------------------------------------------------------------------
       INTEGER, INTENT(OUT) :: I1, I2, I3
       INTEGER, INTENT(IN) :: Ind(3)
!---------------------------------------------------------------------------

       IF ( Ind(1) < Ind(2) ) THEN
          I1 = 1
       ELSE
          I1 = 2
       END IF
       IF ( Ind(I1) > Ind(3) ) THEN
          I1 = 3
       END IF

       SELECT CASE(I1)
       CASE(1)
          IF (Ind(3) > Ind(2)) THEN
             I2 = 2
             I3 = 3
          ELSE
             I2 = 3
             I3 = 2
          END IF
       CASE(2)
         IF (Ind(3) > Ind(1)) THEN
             I2 = 1
             I3 = 3
          ELSE
             I2 = 3
             I3 = 1
          END IF
       CASE(3)
          IF (Ind(2) > Ind(1)) THEN
             I2 = 1
             I3 = 2
          ELSE
             I2 = 2
             I3 = 1
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::TriangleFaceDofsOrdering2nd','Erratic triangular face Indices')
       END SELECT
!---------------------------------------------------------
     END SUBROUTINE TriangleFaceDofsOrdering2nd
!-----------------------------------------------------------

     

!-------------------------------------------------------------
     SUBROUTINE TriangleFaceDofsOrdering2(t,s,Ind)       
!-------------------------------------------------------------------------------
! Returns two unit vectors t and s for spanning constant vector fields
! defined on a triangular face. As a rule for orientation, the vector t is defined 
! as t = Grad L_B - Grad L_A where L_A and L_B are the Lagrange basis functions
! associated with the nodes that has the smallest global indices A and B (A<B).
! Then s = Sqrt(3)* grad L_C, with C corresponding to the largest global index.
!-------------------------------------------------------------------------------
       INTEGER ::  Ind(4)
       REAL(KIND=dp) :: t(3), s(3)
!----------------------------------------------------------
       INTEGER ::  k, A
! -------------------------------------------------------------------
       t = 0.0d0
       s = 0.0d0

       IF ( Ind(1) < Ind(2) ) THEN
          k = 1
       ELSE
          k = 2
       END IF
       IF ( Ind(k) > Ind(3) ) THEN
          k = 3
       END IF
       A = k

       SELECT CASE(A)
       CASE(1)
          IF ( Ind(2) < Ind(3) ) THEN ! B=2, tangent = AB = 12
             t(1) = 1.0d0
             t(2) = 0.0
             s(1) = 0.0d0
             s(2) = 1.0d0
          ELSE ! B=3, tangent = AB = 13
             t(1) = 0.5d0
             t(2) = Sqrt(3.0d0)/2.0d0
             s(1) = Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0
          END IF
       CASE(2)     
          IF ( Ind(1) < Ind(3) ) THEN ! B=1, tangent = AB = 21
             t(1) = -1.0d0
             t(2) = 0.0
             s(1) = 0.0d0
             s(2) = 1.0d0
          ELSE ! B=3, tangent = AB = 23
             t(1) = -0.5d0
             t(2) = Sqrt(3.0d0)/2.0d0
             s(1) = -Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0
          END IF
       CASE(3)
          IF ( Ind(1) < Ind(2) ) THEN ! B=1, tangent = AB = 31
             t(1) = -0.5d0
             t(2) = -Sqrt(3.0d0)/2.0d0
             s(1) = Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0          
          ELSE ! B=2, tangent = AB = 32
             t(1) = 0.5d0
             t(2) = -Sqrt(3.0d0)/2.0d0            
             s(1) = -Sqrt(3.0d0)/2.0d0
             s(2) = -0.5d0       
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::TriangleFaceDofsOrdering','Erratic square face Indices')
       END SELECT
!---------------------------------------------------------
     END SUBROUTINE TriangleFaceDofsOrdering2
!-----------------------------------------------------------

!---------------------------------------------------------
!> This subroutine can be used to create a unique parametrization of
!> quadrilateral faces so that different elements sharing the same
!> face can list the basis functions associated with the face in
!> a unique order. If the face of the reference element is represented
!> by default using two basis vectors e(1,:) and e(2,:), the unique
!> parametrization uses the basis E1 = D1 * e(I1,:) and 
!> E2 = D2 * e(I2,:). 
!----------------------------------------------------------------------
     SUBROUTINE SquareFaceDofsOrdering(I1, I2, D1, D2, Ind, ReverseSign)       
!----------------------------------------------------------------------
       INTEGER, INTENT(OUT) ::  I1, I2      !< Permutation info about coordinate directions
       REAL(KIND=dp), INTENT(OUT) :: D1, D2 !< Sign reversion info related to coordinate directions  
       INTEGER, INTENT(IN) :: Ind(4)        !< The global indices of quadrilateral face
       LOGICAL, OPTIONAL, INTENT(OUT) :: ReverseSign   ! Is e(1,:) x e(2,:) /=  E1 x E2
!----------------------------------------------------------
       INTEGER ::  i, j, k, l, A
       LOGICAL :: ReverseNormal
! -------------------------------------------------------------------
!  Find input for applying an order change and sign reversions to two
!  basis functions associated with a square face. To this end, 
!  find nodes A, B, C such that A has the minimal global index,
!  AB and AC are edges, with C having the largest global index. 
!  Then AB gives the positive direction for the first face DOF and
!  AC gives the positive direction for the second face DOF.
!  REMARK: This convention must be followed when creating basis
!  functions for other element types which are intended to be compatible
!  with the element type to which this rule is applied.
! -------------------------------------------------------------------
       i = 1
       j = 2
       IF ( Ind(i) < Ind(j) ) THEN
          k = i
       ELSE
          k = j
       END IF
       i = 4
       j = 3 
       IF ( Ind(i) < Ind(j) ) THEN
          l = i
       ELSE
          l = j
       END IF
       IF ( Ind(k) > Ind(l) ) THEN
          k = l
       END IF
       A = k

       ReverseNormal = .FALSE.
       
       SELECT CASE(A)
       CASE(1)
          IF ( Ind(2) < Ind(4) ) THEN
             I1 = 1
             I2 = 2
             D1 = 1.0d0
             D2 = 1.0d0
          ELSE
             I1 = 2
             I2 = 1
             D1 = 1.0d0
             D2 = 1.0d0
             ReverseNormal = .TRUE.
          END IF
       CASE(2)
          IF ( Ind(3) < Ind(1) ) THEN
             I1 = 2
             I2 = 1
             D1 = 1.0d0
             D2 = -1.0d0
          ELSE
             I1 = 1
             I2 = 2
             D1 = -1.0d0
             D2 = 1.0d0
             ReverseNormal = .TRUE.
          END IF
       CASE(3)
          IF ( Ind(4) < Ind(2) ) THEN
             I1 = 1
             I2 = 2
             D1 = -1.0d0
             D2 = -1.0d0
          ELSE
             I1 = 2
             I2 = 1
             D1 = -1.0d0
             D2 = -1.0d0
             ReverseNormal = .TRUE.
          END IF
       CASE(4)
          IF ( Ind(1) < Ind(3) ) THEN
             I1 = 2
             I2 = 1
             D1 = -1.0d0
             D2 = 1.0d0
          ELSE
             I1 = 1
             I2 = 2
             D1 = 1.0d0
             D2 = -1.0d0
             ReverseNormal = .TRUE.
          END IF
       CASE DEFAULT
          CALL Fatal('ElementDescription::SquareFaceDofsOrdering','Erratic square face Indices')
       END SELECT

       IF (PRESENT(ReverseSign)) ReverseSign = ReverseNormal
!----------------------------------------------------------
     END SUBROUTINE SquareFaceDofsOrdering
!----------------------------------------------------------

!----------------------------------------------------------------------------------
!>  Returns data for rearranging H(curl)-conforming basis functions so that 
!>  compatibility with the convention for defining global DOFs is attained.
!>  If n basis function value have already been tabulated in the default order
!>  as BasisArray(1:n,:), then SignVec(1:n) * BasisArray(PermVec(1:n),:) gives
!>  the basis vector values corresponding to the global DOFs.
!>  TO DO: support for second-order basis functions, triangles and quads missing
!------------------------------------------------------------------------------------
     SUBROUTINE ReorderingAndSignReversionsData(Element,Nodes,PermVec,SignVec)
!-------------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Element_t), TARGET :: Element        !< Element structure
       TYPE(Nodes_t) :: Nodes                    !< Data corresponding to the classic element nodes
       INTEGER :: PermVec(:)                     !< At exit the permutation vector for performing reordering
       REAL(KIND=dp) :: SignVec(:)               !< At exit the vector for performing sign changes
!---------------------------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh       
       INTEGER, POINTER :: EdgeMap(:,:)
       INTEGER :: SquareFaceMap(4), BrickFaceMap(6,4), PrismSquareFaceMap(3,4), GIndexes(27), DOFs, i, j, k
       INTEGER :: FaceIndices(4), I1, I2, n
       REAL(KIND=dp) :: D1, D2
       LOGICAL :: Parallel
!---------------------------------------------------------------------------------------------------
       Mesh => CurrentModel % Solver % Mesh

       Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)
       
       SignVec = 1.0d0
       
       n = Element % TYPE % NumberOfNodes
       GIndexes(1:n) = Element % NodeIndexes(1:n)
       IF(Parallel) GIndexes(1:n) = Mesh % ParallelInfo % GlobalDofs(GIndexes(1:n))

       SELECT CASE( Element % TYPE % ElementCode / 100 )
       !CASE(3) needs to be done

       !CASE(4) needs to be done

       CASE(5)
         ! NOTE: The Nedelec second family is not yet supported
         EdgeMap => GetEdgeMap(5)
         DO k=1,6
           i = EdgeMap(k,1)
           j = EdgeMap(k,2)
           IF (GIndexes(j)<GIndexes(i)) SignVec(k) = -1.0d0
           PermVec(k) = k
         END DO

       CASE(6)
          EdgeMap => GetEdgeMap(6)
          DO k=1,8
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             IF (GIndexes(j)<GIndexes(i)) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! -----------------------------------------------------
          ! Additional two basis functions on the square face
          ! -----------------------------------------------------
          SquareFaceMap(:) = (/ 1,2,3,4 /)
          FaceIndices(1:4) = GIndexes(SquareFaceMap(1:4))
          CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)
          i = 8
          PermVec(i+1) = i+I1 
          PermVec(i+2) = i+I2
          SignVec(i+1) = D1
          SignVec(i+2) = D2
 
       CASE(7)
          EdgeMap => GetEdgeMap(7)
          DO k=1,9
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             IF (GIndexes(j)<GIndexes(i)) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! ---------------------------------------------------------------------
          ! Additional six basis functions on the square faces (two per face).
          ! ---------------------------------------------------------------------         
          PrismSquareFaceMap(1,:) = (/ 1,2,5,4 /)
          PrismSquareFaceMap(2,:) = (/ 2,3,6,5 /)
          PrismSquareFaceMap(3,:) = (/ 3,1,4,6 /)
          DO k=1,3
             FaceIndices(1:4) = GIndexes(PrismSquareFaceMap(k,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)
             i = 9+(k-1)*2
             PermVec(i+1) = i+I1 
             PermVec(i+2) = i+I2
             SignVec(i+1) = D1
             SignVec(i+2) = D2 
          END DO

       CASE(8)
          EdgeMap => GetEdgeMap(8)
          DO k=1,12
             i = EdgeMap(k,1)
             j = EdgeMap(k,2)
             IF (GIndexes(j)<GIndexes(i)) SignVec(k) = -1.0d0
             PermVec(k) = k
          END DO
          ! ---------------------------------------------------------------------
          ! Additional twelve basis functions on the square faces (two per face).
          ! ---------------------------------------------------------------------         
          BrickFaceMap(1,:) = (/ 1,2,3,4 /)          
          BrickFaceMap(2,:) = (/ 5,6,7,8 /)
          BrickFaceMap(3,:) = (/ 1,2,6,5 /)
          BrickFaceMap(4,:) = (/ 2,3,7,6 /)
          BrickFaceMap(5,:) = (/ 4,3,7,8 /)
          BrickFaceMap(6,:) = (/ 1,4,8,5 /)
          DO k=1,6
             FaceIndices(1:4) = GIndexes(BrickFaceMap(k,1:4))
             CALL SquareFaceDofsOrdering(I1,I2,D1,D2,FaceIndices)
             i = 12+(k-1)*2
             PermVec(i+1) = i+I1 
             PermVec(i+2) = i+I2
             SignVec(i+1) = D1
             SignVec(i+2) = D2 
          END DO
          PermVec(25) = 25
          PermVec(26) = 26         
          PermVec(27) = 27
           
       CASE DEFAULT
          CALL Fatal('ElementDescription::ReorderingAndSignReversionsData','Unsupported element type')
       END SELECT
!----------------------------------------------------------
     END SUBROUTINE ReorderingAndSignReversionsData
!----------------------------------------------------------

!------------------------------------------------------------------------     
!    Given an edge [ij] of a triangle this subroutine returns
!
!    Workbasis(1,1:2) = L_j grad L_i
!    Workbasis(2,1:2) = L_i grad L_j
!
!    and the values of their curl at a given point (u,v). Suitable linear
!    combinations of these functions then give the basic Whitney forms. 
!------------------------------------------------------------------------
     SUBROUTINE EdgeWhitneyComponents2D(WorkBasis, WorkCurlBasis, i, j, u, v)
!------------------------------------------------------------------------
       
       REAL(KIND=dp), INTENT(OUT) :: WorkBasis(2,3), WorkCurlBasis(2,3)
       INTEGER, INTENT(IN) :: i, j
       REAL(KIND=dp), INTENT(IN) :: u, v
!------------------------------------------------------------------------
       REAL(KIND=dp) :: grad_i(2), grad_j(2)
       REAL(KIND=dp) :: grad_svec(2,2), grad_tvec(2,2)
!------------------------------------------------------------------------
       grad_i = dTriangleNodalPBasis(i,u,v)
       grad_j = dTriangleNodalPBasis(j,u,v)

       grad_svec(1,2) = grad_j(2) * grad_i(1)
       grad_svec(2,1) = grad_j(1) * grad_i(2)

       grad_tvec(1,2) = grad_i(2) * grad_j(1)
       grad_tvec(2,1) = grad_i(1) * grad_j(2)

       WorkBasis(1:2,3) = 0.0d0
       WorkBasis(1,1:2) = TriangleNodalPBasis(j,u,v) * grad_i
       WorkBasis(2,1:2) = TriangleNodalPBasis(i,u,v) * grad_j

       WorkCurlBasis(1:2,1:2) = 0.0d0
       WorkCurlBasis(1,3) = grad_svec(2,1) - grad_svec(1,2)
       WorkCurlBasis(2,3) = grad_tvec(2,1) - grad_tvec(1,2)
!------------------------------------------------------------------------
     END SUBROUTINE EdgeWhitneyComponents2D
!------------------------------------------------------------------------


!------------------------------------------------------------------------     
!    Given a face [ijk] of a triangle this subroutine returns
!
!    Workbasis(1,1:2) = L_j L_k grad L_i
!    Workbasis(2,1:2) = L_i L_k grad L_j
!    Workbasis(3,1:2) = L_i L_j grad L_k     
!
!    and the values of their curl at a given point (u,v).
!------------------------------------------------------------------------
     SUBROUTINE FaceWhitneyComponents2D(WorkBasis, WorkCurlBasis, u, v)
!------------------------------------------------------------------------
       
       REAL(KIND=dp), INTENT(OUT) :: WorkBasis(3,3), WorkCurlBasis(3,3)
!       INTEGER, INTENT(IN) :: i, j
       REAL(KIND=dp), INTENT(IN) :: u, v
!------------------------------------------------------------------------
       REAL(KIND=dp) :: grad_i(2), grad_j(2), grad_k(2)
       REAL(KIND=dp) :: grad_svec(2,2), grad_tvec(2,2), grad_hvec(2,2)
!------------------------------------------------------------------------
       WorkBasis(:,3) = 0.0d0
       WorkBasis(1,1:2) = TriangleNodalPBasis(2,u,v) * TriangleNodalPBasis(3,u,v) * dTriangleNodalPBasis(1,u,v) 
       WorkBasis(2,1:2) = TriangleNodalPBasis(1,u,v) * TriangleNodalPBasis(3,u,v) * dTriangleNodalPBasis(2,u,v)
       WorkBasis(3,1:2) = TriangleNodalPBasis(1,u,v) * TriangleNodalPBasis(2,u,v) * dTriangleNodalPBasis(3,u,v)
       
       grad_i = dTriangleNodalPBasis(1,u,v)
       grad_j = dTriangleNodalPBasis(2,u,v)
       grad_k = dTriangleNodalPBasis(3,u,v)
                 
       grad_svec(1,2) = (grad_j(2) * TriangleNodalPBasis(3,u,v) + &
           TriangleNodalPBasis(2,u,v) * grad_k(2)) * grad_i(1)
       grad_svec(2,1) = (grad_j(1) * TriangleNodalPBasis(3,u,v) + &
           TriangleNodalPBasis(2,u,v) * grad_k(1)) * grad_i(2)

       grad_tvec(1,2) = (grad_i(2) * TriangleNodalPBasis(3,u,v)  + &
           TriangleNodalPBasis(1,u,v) * grad_k(2)) * grad_j(1)
       grad_tvec(2,1) = (grad_i(1) * TriangleNodalPBasis(3,u,v)  + &
           TriangleNodalPBasis(1,u,v) * grad_k(1)) * grad_j(2)

       grad_hvec(1,2) = (grad_i(2) * TriangleNodalPBasis(2,u,v)  + &
           TriangleNodalPBasis(1,u,v) * grad_j(2)) * grad_k(1)
       grad_hvec(2,1) = (grad_i(1) * TriangleNodalPBasis(2,u,v)  + &
           TriangleNodalPBasis(1,u,v) * grad_j(1)) * grad_k(2)

       WorkCurlBasis(1:3,1:2) = 0.0d0
       WorkCurlBasis(1,3) = grad_svec(2,1) - grad_svec(1,2)
       WorkCurlBasis(2,3) = grad_tvec(2,1) - grad_tvec(1,2)
       WorkCurlBasis(3,3) = grad_hvec(2,1) - grad_hvec(1,2)
!------------------------------------------------------------------------
     END SUBROUTINE FaceWhitneyComponents2D
!------------------------------------------------------------------------
     
!------------------------------------------------------------------------
     SUBROUTINE WeightedWhitneyForms(WorkBasis, WorkCurlBasis, k, u, v, w)
!------------------------------------------------------------------------
!    Given the kth face of a tetrahedron this subroutine returns
!
!    b_1 = L_k W_{ij}
!    b_2 = L_j W_{ik}
!    b_3 = L_i W_{jk}       
!
!    and the values of their curl at a given point, with W_{ij} denoting
!    the Whitney forms. Here the triangular faces [ijk] are indexed as
!     
!    k=1  [213]
!    k=2  [124]
!    k=3  [234]
!    k=4  [314]     
!------------------------------------------------------------------------
     REAL(KIND=dp), INTENT(OUT) :: WorkBasis(3,3), WorkCurlBasis(3,3)
     INTEGER, INTENT(IN) :: k
     REAL(KIND=dp), INTENT(IN) :: u, v, w
!------------------------------------------------------------------------
     SELECT CASE(k)
     CASE(1) !213
       WorkBasis(1,1) = ((4.0d0*v - Sqrt(2.0d0)*w)*&
           (-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(3.0d0))
       WorkBasis(1,2) = -(u*(4.0d0*v - Sqrt(2.0d0)*w))/24.0d0
       WorkBasis(1,3) = (u*(-2.0d0*Sqrt(2.0d0)*v + w))/24.0d0
       WorkCurlBasis(1,1) = -u/(4.0d0*Sqrt(2.0d0))
       WorkCurlBasis(1,2) = (Sqrt(6.0d0) + 3.0d0*Sqrt(2.0d0)*v - 3.0d0*w)/24.0d0
       WorkCurlBasis(1,3) = (Sqrt(3.0d0) - 3.0d0*v)/6.0d0

       WorkBasis(2,1) = ((4.0d0*v - Sqrt(2.0d0)*w)*(-6.0d0 + 6.0d0*u + &
           2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(96.0d0*Sqrt(3.0d0))
       WorkBasis(2,2) = -((4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - 3.0d0*Sqrt(2.0d0)*w)*&
           (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/288.0d0
       WorkBasis(2,3) = ((Sqrt(3.0d0) + Sqrt(3.0d0)*u - 3.0d0*v)*&
           (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(144.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,1) = -(-6.0d0 + 2.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
           Sqrt(6.0d0)*w)/(16.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,2) = (2.0d0*Sqrt(3.0d0) - 6.0d0*Sqrt(3.0d0)*u + &
           6.0d0*v - 3.0d0*Sqrt(2.0d0)*w)/(48.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,3) = (Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u - 3.0d0*v)/12.0d0

       WorkBasis(3,1) = -((4.0d0*v - Sqrt(2.0d0)*w)*(-6.0d0 - 6.0d0*u + &
           2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(96.0d0*Sqrt(3.0d0))
       WorkBasis(3,2) = ((-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + 3.0d0*Sqrt(2.0d0)*w)* &
           (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/288.0d0
       WorkBasis(3,3) = -((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + 3.0d0*v)* &
           (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(144.0d0*Sqrt(2.0d0))
       WorkCurlBasis(3,1) = -(-6.0d0 - 2.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
           Sqrt(6.0d0)*w)/(16.0d0*Sqrt(2.0d0))
       WorkCurlBasis(3,2) = (-2.0d0*Sqrt(3.0d0) - 6.0d0*Sqrt(3.0d0)*u - 6.0d0*v + &
           3.0d0*Sqrt(2.0d0)*w)/(48.0d0*Sqrt(2.0d0))
       WorkCurlBasis(3,3) = (-Sqrt(3.0d0) - 3.0d0*Sqrt(3.0d0)*u + 3.0d0*v)/12.0d0
       
     CASE(2) !124
       WorkBasis(1,1) = -(w*(-6.0d0 + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(8.0d0*Sqrt(6.0d0))
       WorkBasis(1,2) = (u*w)/(4.0d0*Sqrt(2.0d0))
       WorkBasis(1,3) = (u*w)/8.0d0
       WorkCurlBasis(1,1) = -u/(4.0d0*Sqrt(2.0d0))
       WorkCurlBasis(1,2) = (Sqrt(6.0d0) - Sqrt(2.0d0)*v - 3.0d0*w)/8.0d0
       WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

       WorkBasis(2,1) = -(w*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + &
           Sqrt(6.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(2,2) = (w*(1.0d0 + u - v/Sqrt(3.0d0) - w/Sqrt(6.0d0)))/(8.0d0*Sqrt(2.0d0))
       WorkBasis(2,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)* &
           (-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,1) = (-3.0d0*Sqrt(2.0d0) - Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(2,2) = (Sqrt(6.0d0) + 3.0d0*Sqrt(6.0d0)*u - Sqrt(2.0d0)*v - 3.0d0*w)/16.0d0
       WorkCurlBasis(2,3) =  w/(4.0d0*Sqrt(2.0d0))

       WorkBasis(3,1) = (w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(3,2) = -(w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(48.0d0*Sqrt(2.0d0))
       WorkBasis(3,3) = -((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
           (-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/96.0d0
       WorkCurlBasis(3,1) = (-3.0d0*Sqrt(2.0d0) + Sqrt(2.0d0)*u + Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(3,2) = (-Sqrt(6.0d0) + 3.0d0*Sqrt(6.0d0)*u + Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
       WorkCurlBasis(3,3) = -w/(4.0d0*Sqrt(2.0d0))
       
     CASE(3) ! 234
       WorkBasis(1,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
       WorkBasis(1,2) = (w*(4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u - &
           3.0d0*Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(1,3) = -((1.0d0 + u - Sqrt(3.0d0)*v)*w)/16.0d0
       WorkCurlBasis(1,1) = (-2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u + 3.0d0*Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(1,2) = (-2.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
       WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

       WorkBasis(2,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
       WorkBasis(2,2) = -(w*(-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(2,3) = -((Sqrt(6.0d0) + Sqrt(6.0d0)*u - Sqrt(2.0d0)*v)*&
           (-4.0d0*v + Sqrt(2.0d0)*w))/(32.0d0*Sqrt(3.0d0))
       WorkCurlBasis(2,1) = (2.0d0*Sqrt(2.0d0) + 2.0d0*Sqrt(2.0d0)*u - &
           2.0d0*Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(2,2) = (-4.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
       WorkCurlBasis(2,3) = w/(4.0d0*Sqrt(2.0d0))

       WorkBasis(3,1) = 0.0d0
       WorkBasis(3,2) = (w*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
       WorkBasis(3,3) = -(v*(-6.0d0 - 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
       WorkCurlBasis(3,1) = (2.0d0*Sqrt(2.0d0) + 2.0d0*Sqrt(2.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/8.0d0
       WorkCurlBasis(3,2) = -v/(4.0d0*Sqrt(2.0d0))
       WorkCurlBasis(3,3) = -w/(4.0d0*Sqrt(2.0d0))

     CASE(4) ! 314
       WorkBasis(1,1) = (w*(-2.0d0*Sqrt(2.0d0)*v + w))/16.0d0
       WorkBasis(1,2) = (w*(-4.0d0*Sqrt(3.0d0) + 4.0d0*Sqrt(3.0d0)*u + &
           3.0d0*Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(1,3) = -((-1.0d0 + u + Sqrt(3.0d0)*v)*w)/16.0d0
       WorkCurlBasis(1,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - 3.0d0*Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(1,2) = (-2.0d0*Sqrt(2.0d0)*v + 3.0d0*w)/16.0d0
       WorkCurlBasis(1,3) = w/(2.0d0*Sqrt(2.0d0))

       WorkBasis(2,1) = 0.0d0
       WorkBasis(2,2) = (w*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
       WorkBasis(2,3) = -(v*(-6.0d0 + 6.0d0*u + 2.0d0*Sqrt(3.0d0)*v + Sqrt(6.0d0)*w))/(24.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - Sqrt(6.0d0)*v - Sqrt(3.0d0)*w)/8.0d0
       WorkCurlBasis(2,2) = v/(4.0d0*Sqrt(2.0d0))
       WorkCurlBasis(2,3) =  w/(4.0d0*Sqrt(2.0d0))

       WorkBasis(3,1) = ((2.0d0*Sqrt(2.0d0)*v - w)*w)/16.0d0
       WorkBasis(3,2) = -(w*(-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkBasis(3,3) = ((-Sqrt(3.0d0) + Sqrt(3.0d0)*u + v)*&
           (-4.0d0*v + Sqrt(2.0d0)*w))/(16.0d0*Sqrt(6.0d0))
       WorkCurlBasis(3,1) = (2.0d0*Sqrt(2.0d0) - 2.0d0*Sqrt(2.0d0)*u - &
           2.0d0*Sqrt(6.0d0)*v + Sqrt(3.0d0)*w)/16.0d0
       WorkCurlBasis(3,2) = (4.0d0*Sqrt(2.0d0)*v - 3.0d0*w)/16.0d0
       WorkCurlBasis(3,3) =  -w/(4.0d0*Sqrt(2.0d0))

     CASE DEFAULT
       CALL Fatal('WeightedWhitneyForms', 'A wrong face index')       
     END SELECT
!------------------------------------------------------------------------
   END SUBROUTINE WeightedWhitneyForms
!------------------------------------------------------------------------
     

! --------------------------------------------------------------------------------------
!> This subroutine contains an older design for providing edge element basis functions
!> of the lowest-degree. Obtaining optimal accuracy with these elements may require that 
!> the element map is affine, while the edge basis functions given by the newer design 
!> (the function EdgeElementInfo) should also work on general meshes. 
!------------------------------------------------------------------------
   SUBROUTINE GetEdgeBasis( Element, WBasis, RotWBasis, Basis, dBasisdx )
!------------------------------------------------------------------------
     TYPE(Element_t),TARGET :: Element
     REAL(KIND=dp) :: WBasis(:,:), RotWBasis(:,:), Basis(:), dBasisdx(:,:)
!------------------------------------------------------------------------
     TYPE(Element_t),POINTER :: Edge
     TYPE(Mesh_t), POINTER :: Mesh
     REAL(KIND=dp) :: u,v,w,dudx(3,3),du(3),Base,dBase(3),tBase(3), &
                rBase(3),triBase(3),dtriBase(3,3), G(3,3), F(3,3), detF, detG, &
                EdgeBasis(8,3), CurlBasis(8,3)
     LOGICAL :: Parallel,stat
     INTEGER :: i,j,k,n,nj,nk,i1,i2
     INTEGER, POINTER :: EdgeMap(:,:)
!------------------------------------------------------------------------
     Mesh => CurrentModel % Solver % Mesh

     Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)
     
     IF (Element % TYPE % BasisFunctionDegree>1) THEN
       CALL Fatal('GetEdgeBasis',"Can't handle but linear elements, sorry.") 
     END IF

     SELECT CASE(Element % TYPE % ElementCode / 100)
     CASE(4,7,8)
       n = Element % TYPE % NumberOfNodes
       u = SUM(Basis(1:n)*Element % TYPE % NodeU(1:n))
       v = SUM(Basis(1:n)*Element % TYPE % NodeV(1:n))
       w = SUM(Basis(1:n)*Element % TYPE % NodeW(1:n))

       dudx(1,:) = MATMUL(Element % TYPE % NodeU(1:n),dBasisdx(1:n,:))
       dudx(2,:) = MATMUL(Element % TYPE % NodeV(1:n),dBasisdx(1:n,:))
       dudx(3,:) = MATMUL(Element % TYPE % NodeW(1:n),dBasisdx(1:n,:))

       triBase(1) = 1-u-v
       triBase(2) = u
       triBase(3) = v

       dtriBase(1,:) = -dudx(1,:)-dudx(2,:) 
       dtriBase(2,:) =  dudx(1,:)
       dtriBase(3,:) =  dudx(2,:)
     CASE(6)
       n = Element % TYPE % NumberOfNodes
       u = SUM(Basis(1:n)*Element % TYPE % NodeU(1:n))
       v = SUM(Basis(1:n)*Element % TYPE % NodeV(1:n))
       w = SUM(Basis(1:n)*Element % TYPE % NodeW(1:n))

       G(1,:) = MATMUL(Element % TYPE % NodeU(1:n),dBasisdx(1:n,:))
       G(2,:) = MATMUL(Element % TYPE % NodeV(1:n),dBasisdx(1:n,:))
       G(3,:) = MATMUL(Element % TYPE % NodeW(1:n),dBasisdx(1:n,:))            

       detG =  G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
                  G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
                  G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )
       detF = 1.0d0/detG
       CALL InvertMatrix3x3(G,F,detG)
       
       !------------------------------------------------------------
       ! The basis functions spanning the reference element space and
       ! their Curl with respect to the local coordinates
       ! ------------------------------------------------------------
       EdgeBasis(1,1) = (1.0d0 - v - w)/4.0d0
       EdgeBasis(1,2) = 0.0d0
       EdgeBasis(1,3) = (u*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,1) = u/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,2) = -(-2.0d0 + v + 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(1,3) = 0.25d0

       EdgeBasis(2,1) = 0.0d0
       EdgeBasis(2,2) = (1.0d0 + u - w)/4.0d0
       EdgeBasis(2,3) = (v*(1.0d0 + u - w))/(4.0d0 - 4.0d0*w)
       CurlBasis(2,1) = (2.0d0 + u - 2.0d0*w)/(4.0d0 - 4.0d0*w)
       CurlBasis(2,2) = v/(4.0d0*(-1.0d0 + w))
       CurlBasis(2,3) = 0.25d0       

       EdgeBasis(3,1) = (1.0d0 + v - w)/4.0d0
       EdgeBasis(3,2) = 0.0d0
       EdgeBasis(3,3) = (u*(1.0d0 + v - w))/(4.0d0 - 4.0d0*w)
       CurlBasis(3,1) = u/(4.0d0 - 4.0d0*w)
       CurlBasis(3,2) = (2.0d0 + v - 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(3,3) = -0.25d0

       EdgeBasis(4,1) = 0.0d0
       EdgeBasis(4,2) = (1.0d0 - u - w)/4.0d0
       EdgeBasis(4,3) = (v*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       CurlBasis(4,1) = (-2.0d0 + u + 2.0d0*w)/(4.0d0*(-1.0d0 + w))
       CurlBasis(4,2) = v/(4.0d0 - 4.0d0*w)
       CurlBasis(4,3) = -0.25d0

       EdgeBasis(5,1) = (w*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(5,2) = (w*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(5,3) = (-((-1.0d0 + v)*(-1.0d0 + w)**2) + u*(v - (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(5,1) = -(-1.0d0 + u + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(5,2) = (-1.0d0 + v + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(5,3) = 0.0d0

       EdgeBasis(6,1) = -(w*(-1.0d0 + v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(6,2) = (w*(-1.0d0 - u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(6,3) = (-((-1.0d0 + v)*(-1.0d0 + w)**2) + u*((-1.0d0 + w)**2 + v*(-1.0d0 + 2.0d0*w)))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(6,1) = (1.0d0 + u - w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(6,2) = -(-1.0d0 + v + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(6,3) = 0.0d0    

       EdgeBasis(7,1) = ((1.0d0 + v - w)*w)/(4.0d0*(-1.0d0 + w))
       EdgeBasis(7,2) = ((1.0d0 + u - w)*w)/(4.0d0*(-1.0d0 + w))
       EdgeBasis(7,3) = ((1.0d0 + v)*(-1.0d0 + w)**2 + u*(v + (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(7,1) = (1.0d0 + u - w)/(2.0d0 - 2.0d0*w)
       CurlBasis(7,2) = (1.0d0 + v - w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(7,3) = 0.0d0

       EdgeBasis(8,1) = (w*(-1.0d0 - v + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(8,2) = -(w*(-1.0d0 + u + w))/(4.0d0*(-1.0d0 + w))
       EdgeBasis(8,3) = ((1.0d0 + v)*(-1.0d0 + w)**2 - u*(v + (-1.0d0 + w)**2 - 2.0d0*v*w))/&
            (4.0d0*(-1.0d0 + w)**2)
       CurlBasis(8,1) = (-1.0d0 + u + w)/(2.0d0*(-1.0d0 + w))
       CurlBasis(8,2) = (1.0d0 + v - w)/(2.0d0 - 2.0d0*w)
       CurlBasis(8,3) = 0.0d0

     END SELECT

     EdgeMap => GetEdgeMap(Element % TYPE % ElementCode / 100)
     DO i=1,SIZE(Edgemap,1)
       j = EdgeMap(i,1); k = EdgeMap(i,2)

       nj = Element % Nodeindexes(j)
       nk = Element % Nodeindexes(k)
       IF (Parallel) THEN
         nj=Mesh % ParallelInfo % GlobalDOFs(nj)
         nk=Mesh % ParallelInfo % GlobalDOFs(nk)
       END IF
         
       SELECT CASE(Element % TYPE % ElementCode / 100)
       CASE(3,5)
         WBasis(i,:) = Basis(j)*dBasisdx(k,:) - Basis(k)*dBasisdx(j,:)

         RotWBasis(i,1) = 2.0_dp * ( dBasisdx(j,2) * dBasisdx(k,3) - &
                       dBasisdx(j,3) * dBasisdx(k,2) )
         RotWBasis(i,2) = 2.0_dp * ( dBasisdx(j,3) * dBasisdx(k,1) - &
                       dBasisdx(j,1) * dBasisdx(k,3) )
         RotWBasis(i,3) = 2.0_dp * ( dBasisdx(j,1) * dBasisdx(k,2) - &
                       dBasisdx(j,2) * dBasisdx(k,1) )

       CASE(6)
          !-----------------------------------------------------------------------
          ! Create the referential description of basis functions and their 
          ! spatial curl on the physical element via applying the Piola transform:
          !-----------------------------------------------------------------------
          DO k=1,3
             WBasis(i,k) = SUM( G(1:3,k) * EdgeBasis(i,1:3) )
          END DO
          DO k=1,3
             RotWBasis(i,k) = 1.0d0/DetF * SUM( F(k,1:3) * CurlBasis(i,1:3) )
          END DO

       CASE(7)
         SELECT CASE(i)
          CASE(1)
            j=1;k=2; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(2)
            j=2;k=3; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(3)
            j=3;k=1; Base=(1-w)/2; dBase=-dudx(3,:)/2
          CASE(4)
            j=1;k=2; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(5)
            j=2;k=3; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(6)
            j=3;k=1; Base=(1+w)/2; dBase= dudx(3,:)/2
          CASE(7)
            Base=triBase(1); dBase=dtriBase(1,:); du=dudx(3,:)/2
          CASE(8)
            Base=triBase(2); dBase=dtriBase(2,:); du=dudx(3,:)/2
          CASE(9)
            Base=triBase(3); dBase=dtriBase(3,:); du=dudx(3,:)/2
         END SELECT

         IF(i<=6) THEN
            tBase = (triBase(j)*dtriBase(k,:)-triBase(k)*dtriBase(j,:))
            rBase(1) = 2*Base*(dtriBase(j,2)*dtriBase(k,3)-dtriBase(k,2)*dtriBase(j,3)) + &
                              dBase(2)*tBase(3) - dBase(3)*tBase(2)

            rBase(2) = 2*Base*(dtriBase(j,3)*dtriBase(k,1)-dtriBase(k,3)*dtriBase(j,1)) + &
                              dBase(3)*tBase(1) - dBase(1)*tBase(3)

            rBase(3) = 2*Base*(dtriBase(j,1)*dtriBase(k,2)-dtriBase(k,1)*dtriBase(j,2)) + &
                              dBase(1)*tBase(2) - dBase(2)*tBase(1)

            RotWBasis(i,:)=rBase
            WBasis(i,:)=tBase*Base
         ELSE
            WBasis(i,:)=Base*du
            RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))
            RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))
            RotWBasis(i,3)=(dBase(1)*du(2) - dBase(2)*du(1))
         END IF
       CASE(4)
         SELECT CASE(i)
          CASE(1)
             du=dudx(1,:); Base=(1-v)*(1-w)
             dBase(:)=-dudx(2,:)*(1-w)-(1-v)*dudx(3,:)
          CASE(2)
             du=dudx(2,:); Base=(1+u)*(1-w)
             dBase(:)= dudx(1,:)*(1-w)-(1+u)*dudx(3,:)
          CASE(3)
             du=-dudx(1,:); Base=(1+v)*(1-w)
             dBase(:)= dudx(2,:)*(1-w)-(1+v)*dudx(3,:)
          CASE(4)
             du=-dudx(2,:); Base=(1-u)*(1-w)
             dBase(:)=-dudx(1,:)*(1-w)-(1-u)*dudx(3,:)
         END SELECT

         wBasis(i,:) = Base*du/n
         RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))/n
         RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))/n
         RotWBasis(i,3) = (dBase(1)*du(2) - dBase(2)*du(1))/n
       CASE(8)
         SELECT CASE(i)
          CASE(1)
             du=dudx(1,:); Base=(1-v)*(1-w)
             dBase(:)=-dudx(2,:)*(1-w)-(1-v)*dudx(3,:)
          CASE(2)
             du=dudx(2,:); Base=(1+u)*(1-w)
             dBase(:)= dudx(1,:)*(1-w)-(1+u)*dudx(3,:)
          CASE(3)
             du=dudx(1,:); Base=(1+v)*(1-w)
             dBase(:)= dudx(2,:)*(1-w)-(1+v)*dudx(3,:)
          CASE(4)
             du=dudx(2,:); Base=(1-u)*(1-w)
             dBase(:)=-dudx(1,:)*(1-w)-(1-u)*dudx(3,:)
          CASE(5)
             du=dudx(1,:); Base=(1-v)*(1+w)
             dBase(:)=-dudx(2,:)*(1+w)+(1-v)*dudx(3,:)
          CASE(6)
             du=dudx(2,:); Base=(1+u)*(1+w)
             dBase(:)= dudx(1,:)*(1+w)+(1+u)*dudx(3,:)
          CASE(7)
             du=dudx(1,:); Base=(1+v)*(1+w)
             dBase(:)= dudx(2,:)*(1+w)+(1+v)*dudx(3,:)
          CASE(8)
             du=dudx(2,:); Base=(1-u)*(1+w)
             dBase(:)=-dudx(1,:)*(1+w)+(1-u)*dudx(3,:)
          CASE(9)
             du=dudx(3,:); Base=(1-u)*(1-v)
             dBase(:)=-dudx(1,:)*(1-v)-(1-u)*dudx(2,:)
          CASE(10)
             du=dudx(3,:); Base=(1+u)*(1-v)
             dBase(:)= dudx(1,:)*(1-v)-(1+u)*dudx(2,:)
          CASE(11)
             du=dudx(3,:); Base=(1+u)*(1+v)
             dBase(:)= dudx(1,:)*(1+v)+(1+u)*dudx(2,:)
          CASE(12)
             du=dudx(3,:); Base=(1-u)*(1+v)
             dBase(:)=-dudx(1,:)*(1+v)+(1-u)*dudx(2,:)
         END SELECT

         wBasis(i,:)=Base*du/n
         RotWBasis(i,1)=(dBase(2)*du(3) - dBase(3)*du(2))/n
         RotWBasis(i,2)=(dBase(3)*du(1) - dBase(1)*du(3))/n
         RotWBasis(i,3)=(dBase(1)*du(2) - dBase(2)*du(1))/n
       CASE DEFAULT
         CALL Fatal( 'Edge Basis', 'Not implemented for this element type.')
       END SELECT

       IF( nk < nj ) THEN
         WBasis(i,:) = -WBasis(i,:); RotWBasis(i,:) = -RotWBasis(i,:)
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE GetEdgeBasis
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Return the elementwise number of degrees of freedom and their indexes for
!> a particular solver
!------------------------------------------------------------------------------
   FUNCTION mGetElementDOFs( Indexes, UElement, USolver, NotDG, UMesh ) RESULT(nd)
!------------------------------------------------------------------------------
     INTEGER :: Indexes(:)
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     LOGICAL, OPTIONAL :: NotDG
     TYPE(Mesh_t), OPTIONAL, TARGET :: UMesh
     INTEGER :: nd     
!------------------------------------------------------------------------------
     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent, Face
     TYPE(Mesh_t), POINTER :: Mesh

     LOGICAL :: Found, GB, DGDisable, NeedEdges, Bubbles
     INTEGER :: i,j,k,id, nb, p, NDOFs, MaxNDOFs, EDOFs, MaxEDOFs, FDOFs, MaxFDOFs, BDOFs
     INTEGER :: Ind, ElemFamily, ParentFamily, face_type, face_id
     INTEGER :: NodalIndexOffset, EdgeIndexOffset, FaceIndexOffset
!------------------------------------------------------------------------------
     IF ( PRESENT( USolver ) ) THEN
       Solver => USolver
     ELSE
       Solver => CurrentModel % Solver
     END IF
     
     nd = 0

     IF (.NOT. ASSOCIATED(Solver)) THEN
       CALL Warn('mGetElementDOFS', 'Cannot return DOFs data without knowing solver')
       RETURN
     END IF
     
     IF( PRESENT( UMesh ) ) THEN
       Mesh => UMesh
     ELSE
       Mesh => Solver % Mesh
     END IF
            
     IF ( PRESENT( UElement ) ) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF
     ElemFamily = Element % TYPE % ElementCode / 100

     DGDisable=.FALSE.
     IF (PRESENT(NotDG)) DGDisable=NotDG

     IF ( .NOT. DGDisable .AND. Solver % DG ) THEN
       DO i=1,Element % DGDOFs
         nd = nd + 1
         Indexes(nd) = Element % DGIndexes(i)
       END DO

       IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
         IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
           DO i=1,Element % BoundaryInfo % Left % DGDOFs
             nd = nd + 1
             Indexes(nd) = Element % BoundaryInfo % Left % DGIndexes(i)
           END DO
         END IF
         IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
           DO i=1,Element % BoundaryInfo % Right % DGDOFs
             nd = nd + 1
             Indexes(nd) = Element % BoundaryInfo % Right % DGIndexes(i)
           END DO
         END IF
       END IF

       IF ( nd > 0 ) RETURN
     END IF

     id = Element % BodyId
     IF ( Id==0 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) &
           id = Element % BoundaryInfo % Left % BodyId

       IF (id == 0 .OR. id > CurrentModel % NumberOfBodies) THEN
         IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) &
             id = Element % BoundaryInfo % Right % BodyId
       END IF
     END IF
     !
     ! In some cases it may happen that this function
     ! is called although the BodyId of the element structure hasn't
     ! been set. The following "guess" would be risky if the element
     ! definition depended on body index. It's desirable that
     ! the caller takes care of the creation of the body index so that
     ! the following row need not be considered.
     IF (id==0) id=1


     IF (SIZE(Solver % Def_Dofs,2) < id) CALL Fatal('mGetElementDOFS', &
         'Indexing outside array bounds: '//I2S(SIZE(Solver % Def_Dofs,2))//' vs. '//I2S(id))
     
     IF (.NOT.ASSOCIATED(Mesh)) THEN
       IF ( Solver % Def_Dofs(ElemFamily,id,1)>0 ) THEN  
         CALL Warn('mGetElementDOFS', &
             'Solver mesh unknown, the node indices are returned')
         MaxNDOFs = 1
       ELSE
         CALL Warn('mGetElementDOFS', &
             'Solver mesh unknown, no indices returned')
         RETURN
       END IF
     ELSE
       MaxNDOFs = Mesh % MaxNDOFs
     END IF
     NodalIndexOffset = MaxNDOFs * Mesh % NumberOfNodes     

     NDOFs = Solver % Def_Dofs(ElemFamily,id,1)
     IF (NDOFs > 0) THEN
       DO i=1,Element % TYPE % NumberOfNodes
         DO j=1,NDOFs
           nd = nd + 1
           Indexes(nd) = MaxNDOFs * (Element % NodeIndexes(i)-1) + j
         END DO
       END DO
     END IF

     ! The DOFs of advanced elements cannot be returned without knowing mesh
     ! ---------------------------------------------------------------------
     IF (.NOT.ASSOCIATED(Mesh)) RETURN

     NeedEdges = .FALSE.
     DO i=2,SIZE(Solver % Def_Dofs,3)
       IF (Solver % Def_Dofs(ElemFamily, id, i)>=0) THEN
         NeedEdges = .TRUE.
         EXIT
       END IF
     END DO

     IF (.NOT. NeedEdges) THEN
       !
       ! Check whether face DOFs have been generated by "-quad_face b: ..." or
       ! "-tri_face b: ..."
       !
       IF (ElemFamily == 3 .OR. ElemFamily == 4) THEN
         IF (Solver % Def_Dofs(6+ElemFamily, id, 5)>=0) NeedEdges = .TRUE.
       ELSE
         !
         ! Check finally if 3-D faces are associated with face bubbles
         !
         IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           DO j=1,Element % TYPE % NumberOfFaces
             Face => Mesh % Faces(Element % FaceIndexes(j))
             face_type = Face % TYPE % ElementCode/100
             k = 0
             IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
               face_id  = Face % BoundaryInfo % Left % BodyId
               k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF
             IF (k == 0) THEN
               IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
                 face_id = Face % BoundaryInfo % Right % BodyId
                 k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
               END IF
             END IF
             IF (k > 0) THEN
               NeedEdges = .TRUE.
               EXIT
             END IF
           END DO
         END IF
       END IF
     END IF

     IF ( .NOT. NeedEdges ) RETURN

     MaxFDOFs = Mesh % MaxFaceDOFs
     MaxEDOFs = Mesh % MaxEdgeDOFs
     EdgeIndexOffset = MaxEDOFs * Mesh % NumberOfEdges
     FaceIndexOffset = MaxFDOFs * Mesh % NumberOfFaces

BLOCK
  LOGICAL  :: EdgesDone, FacesDone
  TYPE(Element_t), POINTER :: Edge

       EdgesDone = .FALSE.
       FacesDone = .FALSE.

       IF ( ASSOCIATED(Element % EdgeIndexes) ) THEN
         EdgesDone = .TRUE.
         DO j=1,Element % TYPE % NumberOfEdges
           Edge => Mesh % Edges( Element % EdgeIndexes(j) )
           IF( Edge % Type % ElementCode == Element % Type % ElementCode) THEN
             IF ( .NOT. (Solver % GlobalBubbles .AND. &
                   Element % BodyId>0.AND.ASSOCIATED(Element % BoundaryInfo)) ) THEN
               EdgesDone = .FALSE.
               CYCLE
             END IF
           END IF

           EDOFs = 0
           IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
             EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
           ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect when p varies over mesh; cf. what is done in InitialPermutation
             EDOFs = getEdgeDOFs(Element, Solver % Def_Dofs(ElemFamily,id,6))
           END IF

           DO i=1,EDOFs
             nd = nd + 1
             Indexes(nd) = MaxEDOFs*(Element % EdgeIndexes(j)-1) + &
                 i + NodalIndexOffset
           END DO
         END DO
       END IF

       IF ( ASSOCIATED(Element % FaceIndexes) ) THEN
         FacesDone = .TRUE.
         DO j=1,Element % TYPE % NumberOfFaces
           Face => Mesh % Faces( Element % FaceIndexes(j) )

           IF (Face % Type % ElementCode == Element % Type % ElementCode) THEN
             IF ( .NOT. (Solver % GlobalBubbles .AND. &
                 Element % BodyId>0.AND.ASSOCIATED(Element % BoundaryInfo)) ) THEN
               FacesDone = .FALSE.
               CYCLE
             END IF
           END IF

           k = MAX(0,Solver % Def_Dofs(ElemFamily,id,3))
           IF (k == 0) THEN
             !
             ! NOTE: This depends on what face dofs have been introduced
             ! by using the construct "-quad_face b: ..." and
             ! "-tri_face b: ..."
             !
             face_type = Face % TYPE % ElementCode/100
             IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
               face_id  = Face % BoundaryInfo % Left % BodyId
               k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
             END IF
             IF (k == 0) THEN
               IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
                 face_id = Face % BoundaryInfo % Right % BodyId
                 k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
               END IF
             END IF
           END IF

           FDOFs = 0
           IF (k > 0) THEN
             FDOFs = k
           ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect when p varies over mesh; cf. what is done in InitialPermutation
             FDOFs = getFaceDOFs(Element,Solver % Def_Dofs(ElemFamily,id,6),j,Face)
           END IF

           DO i=1,FDOFs
             nd = nd + 1
             Indexes(nd) = MaxFDOFs*(Element % FaceIndexes(j)-1) + i + &
                 NodalIndexOffset + EdgeIndexOffset
           END DO
         END DO
       END IF

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN

       IF (isActivePelement(Element, Solver)) THEN
         Parent => Element % pDefs % LocalParent
       ELSE
         Parent => Element % BoundaryInfo % Left
         IF (.NOT.ASSOCIATED(Parent) ) &
             Parent => Element % BoundaryInfo % Right
       END IF
       IF (.NOT.ASSOCIATED(Parent) ) RETURN
       ParentFamily = Parent % TYPE % ElementCode / 100

       SELECT CASE(ElemFamily)
       CASE(2)
         IF ( .NOT. EdgesDone .AND. ASSOCIATED(Parent % EdgeIndexes) ) THEN
           IF ( isActivePElement(Element, Solver) ) THEN
             Ind=Element % PDefs % LocalNumber
           ELSE
             DO Ind=1,Parent % TYPE % NumberOfEdges
               Edge => Mesh % Edges(Parent % EdgeIndexes(ind))
               k = 0
               DO i=1,Edge % TYPE % NumberOfNodes
                 DO j=1,Element % TYPE % NumberOfNodes
                   IF ( Edge % NodeIndexes(i)==Element % NodeIndexes(j) ) k=k+1
                 END DO
               END DO
               IF ( k==Element % TYPE % NumberOfNodes) EXIT
             END DO
           END IF

           EDOFs = 0
           IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
             EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
           ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
             EDOFs = getEdgeDOFs(Parent, Solver % Def_Dofs(ParentFamily,id,6))
           END IF

           DO i=1,EDOFs
             nd = nd + 1
             Indexes(nd) = MaxEDOFs*(Parent % EdgeIndexes(Ind)-1) + &
                 i + NodalIndexOffset
           END DO
         END IF

       CASE(3,4)
         IF ( .NOT. FacesDone .AND. ASSOCIATED( Parent % FaceIndexes ) ) THEN

           IF ( isActivePElement(Element, Solver) ) THEN
             Ind=Element % PDefs % LocalNumber
           ELSE
             DO Ind=1,Parent % TYPE % NumberOfFaces
               Face => Mesh % Faces(Parent % FaceIndexes(ind))
               k = 0
               DO i=1,Face % TYPE % NumberOfNodes
                 DO j=1,Element % TYPE % NumberOfNodes
                   IF ( Face % NodeIndexes(i)==Element % NodeIndexes(j)) k=k+1
                 END DO
               END DO
               IF ( k==Face % TYPE % NumberOfNodes) EXIT
             END DO
           END IF

           IF (Ind >= 1 .AND. Ind <= Parent % Type % NumberOfFaces) THEN

             IF (ASSOCIATED(Element % FaceIndexes).AND. isActivePelement(Element, Solver) ) THEN
               Face => Mesh % Faces(Element % PDefs % localParent % Faceindexes(Ind))
             ELSE
               Face => Element
             END IF

             IF (.NOT.EdgesDone .AND. ASSOCIATED(Face % EdgeIndexes)) THEN
               DO j=1,Face % TYPE % NumberOFEdges
                 Edge => Mesh % Edges(Face % EdgeIndexes(j))

                 EDOFs = 0
                 IF (Solver % Def_Dofs(ElemFamily,id,2) >= 0) THEN
                   EDOFs = Solver % Def_Dofs(ElemFamily,id,2)
                 ELSE IF (Solver % Def_Dofs(ElemFamily,id,6) > 1) THEN
! TO DO: This is not yet perfect when p varies over mesh; cf. what is done in InitialPermutation
                   EDOFs = getEdgeDOFs(Element, Solver % Def_Dofs(ElemFamily,id,6))
                 END IF

                 DO i=1,EDOFs
                   nd = nd + 1
                   Indexes(nd) = MaxEDOFs*(Face % EdgeIndexes(j)-1) + &
                       i + NodalIndexOffset                   
                 END DO
               END DO
             END IF
             
             FDOFs = 0
             IF (Solver % Def_Dofs(ParentFamily,id,6) > 1) THEN
               FDOFs = getFaceDOFs(Parent,Solver % Def_Dofs(ParentFamily,id,6),Ind,Face)
             ELSE
               k = MAX(0,Solver % Def_Dofs(ElemFamily,id,3))
               IF (k == 0) THEN
                 !
                 ! NOTE: This depends on what dofs have been introduced
                 ! by using the construct "-quad_face b: ..." and
                 ! "-tri_face b: ..."
                 !
                 face_type = Face % TYPE % ElementCode/100
                 IF (ASSOCIATED(Face % BoundaryInfo % Left)) THEN
                   face_id  = Face % BoundaryInfo % Left % BodyId
                   k = MAX(0,Solver % Def_Dofs(face_type+6,face_id,5))
                 END IF
                 IF (k == 0) THEN
                   IF (ASSOCIATED(Face % BoundaryInfo % Right)) THEN
                     face_id = Face % BoundaryInfo % Right % BodyId
                     k = MAX(k,Solver % Def_Dofs(face_type+6,face_id,5))
                   END IF
                 END IF
               END IF

               IF (k > 0) THEN
                 FDOFs = k
               END IF
             END IF

             DO i=1,FDOFs
               nd = nd + 1
               Indexes(nd) = MaxFDOFs*(Parent % FaceIndexes(Ind)-1) + i + &
                   NodalIndexOffset + EdgeIndexOffset
             END DO
           END IF
         END IF
       END SELECT
     ELSE
       IF (ASSOCIATED(Element % BubbleIndexes) .AND. Solver % GlobalBubbles) THEN
         BDOFs = 0
         nb = Solver % Def_Dofs(ElemFamily,id,5)
         p = Solver % Def_Dofs(ElemFamily,id,6)
         IF (nb >= 0 .OR. p >= 1) THEN
           IF (p > 1) BDOFs = GetBubbleDOFs(Element, p)
           BDOFs = MAX(nb, BDOFs)
         ELSE
           IF (ASSOCIATED(Solver % Values)) THEN
             Bubbles = ListGetLogical(Solver % Values, 'Bubbles', Found )
             ! The following is not a right way to obtain the bubble count
             ! in order to support solverwise definitions
             IF (Bubbles) BDOFs = SIZE(Element % BubbleIndexes)
           END IF
         END IF
         DO i=1,BDOFs
           nd = nd + 1
           Indexes(nd) = NodalIndexOffset + EdgeIndexOffset + FaceIndexOffset + &
               Element % BubbleIndexes(i)
         END DO
       END IF
     END IF
   END BLOCK

!------------------------------------------------------------------------------
  END FUNCTION mGetElementDOFs
!------------------------------------------------------------------------------

#ifdef HAVE_QP
!------------------------------------------------------------------------------
!>    Check element by comparing determinants of the metric tensor computed
!>    in double and quad precision.
!------------------------------------------------------------------------------
   FUNCTION CheckMetric(nDOFs,Elm,Nodes,dLBasisdx) RESULT(Success)
!------------------------------------------------------------------------------
     INTEGER :: nDOFs                !< Number of active nodes in element
     TYPE(Element_t)  :: Elm         !< Element structure
     TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
     REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
     LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: GeomId     
     INTEGER :: cdim,dim,i,j,k,n,imin,jmin
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z

     INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(24)     

     REAL(KIND=dp) :: dp_dx(3,3),dp_G(3,3),dp_GI(3,3),dp_s, dp_DetG
     REAL(KIND=qp) :: qp_dx(3,3),qp_G(3,3),qp_GI(3,3),qp_s, qp_DetG, eps
!------------------------------------------------------------------------------
     success = .TRUE.

     x => Nodes % x
     y => Nodes % y
     z => Nodes % z

     cdim = CoordinateSystemDimension()
     n = MIN( SIZE(x), nDOFs )
     dim  = elm % TYPE % DIMENSION

     eps = 1.0d-6
!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       dp_dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
       dp_dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
       dp_dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )

       qp_dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
       qp_dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
       qp_dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
     END DO
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
     DO i=1,dim
       DO j=1,dim
         dp_s = 0.0_dp
         qp_s = 0.0_dp
         DO k=1,cdim
           dp_s = dp_s + dp_dx(k,i)*dp_dx(k,j)
           qp_s = qp_s + qp_dx(k,i)*qp_dx(k,j)
         END DO
         dp_G(i,j) = dp_s
         qp_G(i,j) = qp_s
       END DO
     END DO

!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
     CASE (1)
       dp_DetG  = dp_G(1,1)
       qp_DetG  = qp_G(1,1)

!------------------------------------------------------------------------------
!      Surface elements
!------------------------------------------------------------------------------
     CASE (2)
       dp_DetG = ( dp_G(1,1)*dp_G(2,2) - dp_G(1,2)*dp_G(2,1) )
       qp_DetG = ( qp_G(1,1)*qp_G(2,2) - qp_G(1,2)*qp_G(2,1) )

!------------------------------------------------------------------------------
!      Volume elements
!------------------------------------------------------------------------------
     CASE (3)
       dp_DetG = dp_G(1,1) * ( dp_G(2,2)*dp_G(3,3) - dp_G(2,3)*dp_G(3,2) ) + &
                 dp_G(1,2) * ( dp_G(2,3)*dp_G(3,1) - dp_G(2,1)*dp_G(3,3) ) + &
                 dp_G(1,3) * ( dp_G(2,1)*dp_G(3,2) - dp_G(2,2)*dp_G(3,1) )

       qp_DetG = qp_G(1,1) * ( qp_G(2,2)*qp_G(3,3) - qp_G(2,3)*qp_G(3,2) ) + &
                 qp_G(1,2) * ( qp_G(2,3)*qp_G(3,1) - qp_G(2,1)*qp_G(3,3) ) + &
                 qp_G(1,3) * ( qp_G(2,1)*qp_G(3,2) - qp_G(2,2)*qp_G(3,1) )
     END SELECT
     
     Success = ABS(dp_detG-qp_detG) <= eps*ABS(qp_DetG)
!------------------------------------------------------------------------------
   END FUNCTION CheckMetric
!------------------------------------------------------------------------------
#endif
   
!------------------------------------------------------------------------------
!>    Compute contravariant metric tensor (=J^TJ)^-1 of element coordinate
!>    system, and square root of determinant of covariant metric tensor
!>    (=sqrt(det(J^TJ)))
!------------------------------------------------------------------------------
   FUNCTION ElementMetric(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap) RESULT(Success)
!------------------------------------------------------------------------------
     INTEGER :: nDOFs                !< Number of active nodes in element
     TYPE(Element_t)  :: Elm         !< Element structure
     TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
     REAL(KIND=dp) :: Metric(:,:)    !< Contravariant metric tensor
     REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
     REAL(KIND=dp) :: DetG           !< SQRT of determinant of metric tensor
     REAL(KIND=dp) :: LtoGMap(3,3)   !< Transformation to obtain the referential description of the spatial gradient
     LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3),s,smin,eps=0
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     INTEGER :: GeomId     
     INTEGER :: cdim,dim,i,j,k,n,imin,jmin
!------------------------------------------------------------------------------
     success = .TRUE.

     x => Nodes % x
     y => Nodes % y
     z => Nodes % z

     cdim = CoordinateSystemDimension()
     n = MIN( SIZE(x), nDOFs )
     dim  = elm % TYPE % DIMENSION

#ifdef HAVE_QP
     IF(Elm % Status == 2) THEN
       IF (ElementMetricQP(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap)) RETURN
       GOTO 100
     END IF
#endif

     eps = (EPSILON(eps))**dim
!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
       dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
       dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
     END DO
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
     DO i=1,dim
       DO j=1,dim
         s = 0.0_dp
         DO k=1,cdim
           s = s + dx(k,i)*dx(k,j)
         END DO
         G(i,j) = s
       END DO
     END DO
!    G(1:dim,1:dim) = MATMUL( TRANSPOSE(dx(1:cdim,1:dim)),dx(1:cdim,1:dim) )
!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
     CASE (1)
       DetG  = G(1,1)

       IF ( DetG <= eps ) GOTO 100

       Metric(1,1) = 1.0d0 / DetG
       DetG  = SQRT( DetG )

!------------------------------------------------------------------------------
!      Surface elements
!------------------------------------------------------------------------------
     CASE (2)
       DetG = ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )

       IF ( DetG <= eps ) GOTO 100

       Metric(1,1) =  G(2,2) / DetG
       Metric(1,2) = -G(1,2) / DetG
       Metric(2,1) = -G(2,1) / DetG
       Metric(2,2) =  G(1,1) / DetG
       DetG = SQRT(DetG)

!------------------------------------------------------------------------------
!      Volume elements
!------------------------------------------------------------------------------
     CASE (3)
       DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
              G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
              G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )

       IF ( DetG <= eps ) GOTO 100

       CALL InvertMatrix3x3(G,GI,detG)
       Metric = GI
       DetG = SQRT(DetG)
     END SELECT
     
!--------------------------------------------------------------------------------------
!    Construct a transformation X = LtoGMap such that (grad B)(f(p)) = X(p) Grad b(p),
!    with Grad the gradient with respect to the reference element coordinates p and 
!    the referential description of the spatial field B(x) satisfying B(f(p)) = b(p).
!    If cdim > dim (e.g. a surface embedded in the 3-dimensional space), X is
!    the transpose of the pseudo-inverse of Grad f.
!-------------------------------------------------------------------------------
     DO i=1,cdim
       DO j=1,dim
         s = 0.0d0
         DO k=1,dim
           s = s + dx(i,k) * Metric(k,j)
         END DO
         LtoGMap(i,j) = s
       END DO
     END DO
!    LtoGMap(1:cdim,1:dim) = MATMUL(dx(1:cdim,1:dim), Metric(1:dim,1:dim) )

     ! Return here also implies success = .TRUE.
     RETURN

100  CONTINUE

#ifdef HAVE_QP
     ! Try recursively with quadratic precision.
     ! With just double precision for very flat elements the DetJ may be poorly evaluated. 
     IF( Elm % Status /= 2) THEN
       Success = ElementMetricQP(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap) 
       IF( Success ) RETURN
     END IF
#endif
     
     WRITE( Message,'(A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex
     CALL Error( 'ElementMetric', Message )
     
     IF( ASSOCIATED( Elm % BoundaryInfo ) ) THEN
       WRITE( Message,'(A,I0,A,ES14.6)') 'Boundary Id: ',Elm % BoundaryInfo % Constraint,' DetG:',DetG
     ELSE
       WRITE( Message,'(A,I0,A,ES14.6)') 'Body Id: ',Elm % BodyId,' DetG:',DetG
     END IF
     CALL Info( 'ElementMetric', Message, Level=3 )

     DO i=1,n
       WRITE( Message,'(A,I0,A,3ES14.6)') 'Node: ',i,' Coord:',x(i),y(i),z(i)       
       CALL Info( 'ElementMetric', Message, Level=3 )
     END DO

     ! Find the two nodes closest to each other:
     smin = HUGE(smin)
     DO i=1,n
       DO j=i+1,n
         s = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
         IF( s < smin ) THEN
           imin = i
           jmin = j
           smin = s           
         END IF
       END DO
     END DO
     smin = SQRT(smin)

     WRITE( Message,'(A,I0,A,I0,A,I0,A,I0,A,ES14.6)') 'Closest distance: ',imin,'-',jmin,&
         ' (',Elm % NodeIndexes(imin),'-',Elm % NodeIndexes(jmin),') |dCoord|:',smin
     CALL Info( 'ElementMetric', Message, Level=3 )

     IF ( cdim < dim ) THEN
       WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
       CALL Info( 'ElementMetric', Message, Level=3 )
     END IF

!------------------------------------------------------------------------------
   END FUNCTION ElementMetric
!------------------------------------------------------------------------------

#ifdef HAVE_QP
!------------------------------------------------------------------------------
! Quadratic precision version of the previous that is called when the DetJ appear
! to be close to zero or negative. 
!------------------------------------------------------------------------------
   FUNCTION ElementMetricQP(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap) RESULT(Success)
!------------------------------------------------------------------------------
     INTEGER :: nDOFs                !< Number of active nodes in element
     TYPE(Element_t)  :: Elm         !< Element structure
     TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
     REAL(KIND=dp) :: Metric(:,:)    !< Contravariant metric tensor
     REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
     REAL(KIND=dp) :: DetG           !< SQRT of determinant of metric tensor
     REAL(KIND=dp) :: LtoGMap(3,3)   !< Transformation to obtain the referential description of the spatial gradient
     LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     INTEGER :: GeomId     
     INTEGER :: cdim,dim,i,j,k,n

! Local Quadratic precision variables     
     INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(24)     
     REAL(KIND=qp) :: dx(3,3),G(3,3),GI(3,3),s,DetGqp
!------------------------------------------------------------------------------
     success = .FALSE.

     x => Nodes % x
     y => Nodes % y
     z => Nodes % z

     cdim = CoordinateSystemDimension()
     n = MIN( SIZE(x), nDOFs )
     dim  = elm % TYPE % DIMENSION
     DetG = 0.0_dp

!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
       dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
       dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
     END DO
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
     DO i=1,dim
       DO j=1,dim
         s = 0.0d0
         DO k=1,cdim
           s = s + dx(k,i)*dx(k,j)
         END DO
         G(i,j) = s
       END DO
     END DO
!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
     CASE (1)
       DetGqp  = G(1,1)

       IF ( DetGqp <= TINY( DetG ) ) RETURN

       Metric(1,1) = 1.0d0 / DetGqp

!------------------------------------------------------------------------------
!      Surface elements
!------------------------------------------------------------------------------
     CASE (2)
       DetGqp = ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )

       IF ( DetGqp <= TINY( DetG ) ) RETURN

       Metric(1,1) =  G(2,2) / DetGqp
       Metric(1,2) = -G(1,2) / DetGqp
       Metric(2,1) = -G(2,1) / DetGqp
       Metric(2,2) =  G(1,1) / DetGqp

!------------------------------------------------------------------------------
!      Volume elements
!------------------------------------------------------------------------------
     CASE (3)
       DetGqp = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
           G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
           G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )

       IF ( DetGqp <= TINY( DetG ) ) RETURN

       CALL InvertMatrix3x3QP( G,GI,detGqp )
       Metric = GI
     END SELECT

     DetG = SQRT(DetGqp)     
     Success = .TRUE.
     
!--------------------------------------------------------------------------------------
     DO i=1,cdim
       DO j=1,dim
         s = 0.0d0
         DO k=1,dim
           s = s + dx(i,k) * Metric(k,j)
         END DO
         LtoGMap(i,j) = s
       END DO
     END DO
     
!------------------------------------------------------------------------------
   END FUNCTION ElementMetricQP
!------------------------------------------------------------------------------
#endif
   
!------------------------------------------------------------------------------
   FUNCTION ElementMetricVec( Elm, Nodes, nc, ndof, DetJ, nbmax, dLBasisdx, LtoGMap) RESULT(AllSuccess)
!------------------------------------------------------------------------------
     TYPE(Element_t)  :: Elm                                 !< Element structure
     TYPE(Nodes_t)    :: Nodes                               !< element nodal coordinates
     INTEGER, INTENT(IN) :: nc                               !< Number of points to map
     INTEGER :: ndof                                         !< Number of active nodes in element
     REAL(KIND=dp) :: DetJ(VECTOR_BLOCK_LENGTH)              !< SQRT of determinant of element coordinate metric at each point
     INTEGER, INTENT(IN) :: nbmax                            !< Maximum total number of basis functions in local basis
     REAL(KIND=dp) :: dLBasisdx(VECTOR_BLOCK_LENGTH,nbmax,3) !< Derivatives of element basis function with 
                                                             !<  respect to local coordinates at each point
     REAL(KIND=dp) :: LtoGMap(VECTOR_BLOCK_LENGTH,3,3)       !< Mapping between local and global coordinates
     LOGICAL :: AllSuccess                  !< Returns .FALSE. if some point in element is degenerate
!------------------------------------------------------------------------------
!       Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: dx(VECTOR_BLOCK_LENGTH,3,3)
     REAL(KIND=dp) :: Metric(VECTOR_BLOCK_LENGTH,6), &
             G(VECTOR_BLOCK_LENGTH,6)       ! Symmetric Metric(nc,3,3) and G(nc,3,3)

     REAL(KIND=dp) :: s
     INTEGER :: cdim,dim,i,j,k,l,n,ip, jj, kk
     INTEGER :: ldbasis, ldxyz, utind
!DIR$ ATTRIBUTES ALIGN:64::Metric
!DIR$ ATTRIBUTES ALIGN:64::dx
!DIR$ ATTRIBUTES ALIGN:64::G
!DIR$ ASSUME_ALIGNED dLBasisdx:64, LtoGMap:64, DetJ:64
     !------------------------------------------------------------------------------
     AllSuccess = .TRUE.

     ! Coordinates (single array)
     n = MIN( SIZE(Nodes % x, 1), ndof )

     ! Dimensions (coordinate system and element)
     cdim = CoordinateSystemDimension()
     dim  = elm % TYPE % DIMENSION

     ! Leading dimensions for local basis and coordinate arrays
     ldbasis = SIZE(dLBasisdx, 1)
     ldxyz = SIZE(Nodes % xyz, 1)

     ! For linear, extruded and otherwise regular elements mapping has to be computed
     ! only once, the problem is to identify these cases...
     !------------------------------------------------------------------------------
     !       Partial derivatives of global coordinates with respect to local coordinates
     !------------------------------------------------------------------------------
     ! Avoid DGEMM calls for nc small
     IF (nc < VECTOR_SMALL_THRESH) THEN
       DO l=1,dim
         DO j=1,3
           dx(1:nc,j,l)=REAL(0,dp)
           DO k=1,n
!DIR$ UNROLL
             DO i=1,nc
               dx(i,j,l)=dx(i,j,l)+dLBasisdx(i,k,l)*Nodes % xyz(k,j)
             END DO
           END DO
         END DO
       END DO
     ELSE
       DO i=1,dim
         CALL DGEMM('N','N',nc, 3, n, &
                 REAL(1,dp), dLbasisdx(1,1,i), ldbasis, &
                 Nodes % xyz, ldxyz, REAL(0, dp), dx(1,1,i), VECTOR_BLOCK_LENGTH)
       END DO
     END IF
     !------------------------------------------------------------------------------
     !       Compute the covariant metric tensor of the element coordinate system (symmetric)
     !------------------------------------------------------------------------------
     ! Linearized upper triangular indices for accesses to G
     ! | (1,1) (1,2) (1,3) | = | 1 2 4 |
     ! |       (2,2) (2,3) |   |   3 5 |
     ! |             (3,3) |   |     6 |
     ! G is symmetric, compute only the upper triangular part of G=dx^Tdx
!DIR$ LOOP COUNT MAX=3
     DO j=1,dim
!DIR$ LOOP COUNT MAX=3
       DO i=1,j
!DIR$ INLINE
         utind = GetSymmetricIndex(i,j)
         SELECT CASE (cdim)
         CASE(1)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)
           END DO
         CASE(2)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)
           END DO
         CASE(3)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)+dx(l,3,i)*dx(l,3,j)
           END DO
         END SELECT
       END DO
     END DO

     !------------------------------------------------------------------------------
     !       Convert the metric to contravariant base, and compute the SQRT(DetG)
     !------------------------------------------------------------------------------
     SELECT CASE( dim )
       !------------------------------------------------------------------------------
       !       Line elements
       !------------------------------------------------------------------------------
     CASE (1)
       ! Determinants
       ! DetJ(1:nc)  = G(1:nc,1,1)
       DetJ(1:nc)  = G(1:nc,1)

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         !_ELMER_OMP_SIMD
         DO i=1,nc
           ! Metric(i,1,1) = REAL(1,dp)/DetJ(i)
           Metric(i,1) = REAL(1,dp)/DetJ(i)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT( DetJ(i))
         END DO
       END IF


       !------------------------------------------------------------------------------
       !       Surface elements
       !------------------------------------------------------------------------------
     CASE (2)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = ( G(i,1,1)*G(i,2,2) - G(i,1,2)*G(i,2,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*G(i,3)-G(i,2)*G(i,2)
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp)/DetJ(i)
           ! G is symmetric
           ! All in one go, with redundancies eliminated
           Metric(i,1) =  s*G(i,3)
           Metric(i,2) = -s*G(i,2)
           Metric(i,3) =  s*G(i,1)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
       !------------------------------------------------------------------------------
       !       Volume elements
       !------------------------------------------------------------------------------
     CASE (3)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = G(i,1,1) * ( G(i,2,2)*G(i,3,3) - G(i,2,3)*G(i,3,2) ) + &
         !           G(i,1,2) * ( G(i,2,3)*G(i,3,1) - G(i,2,1)*G(i,3,3) ) + &
         !           G(i,1,3) * ( G(i,2,1)*G(i,3,2) - G(i,2,2)*G(i,3,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*(G(i,3)*G(i,6)-G(i,5)*G(i,5)) + &
                 G(i,2)*(G(i,5)*G(i,4)-G(i,2)*G(i,6)) + &
                 G(i,4)*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp) / DetJ(i)
           ! Metric(i,1,1) =  s * (G(i,2,2)*G(i,3,3) - G(i,3,2)*G(i,2,3))
           ! Metric(i,2,1) = -s * (G(i,2,1)*G(i,3,3) - G(i,3,1)*G(i,2,3))
           ! Metric(i,3,1) =  s * (G(i,2,1)*G(i,3,2) - G(i,3,1)*G(i,2,2))
           ! G is symmetric

           ! All in one go, with redundancies eliminated
           Metric(i,1)= s*(G(i,3)*G(i,6)-G(i,5)*G(i,5))
           Metric(i,2)=-s*(G(i,2)*G(i,6)-G(i,4)*G(i,5))
           Metric(i,3)= s*(G(i,1)*G(i,6)-G(i,4)*G(i,4))
           Metric(i,4)= s*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
           Metric(i,5)=-s*(G(i,1)*G(i,5)-G(i,2)*G(i,4))
           Metric(i,6)= s*(G(i,1)*G(i,3)-G(i,2)*G(i,2))
         END DO

         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
     END SELECT

     IF (AllSuccess) THEN
       SELECT CASE(dim)
       CASE(1)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1)
           END DO
         END DO
       CASE(2)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3)
           END DO
         END DO
       CASE(3)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2) + dx(l,i,3)*Metric(l,4)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3) + dx(l,i,3)*Metric(l,5)
             LtoGMap(l,i,3) = dx(l,i,1)*Metric(l,4) + dx(l,i,2)*Metric(l,5) + dx(l,i,3)*Metric(l,6)
           END DO
         END DO
       END SELECT
     ELSE

       ! Degenerate element!
       WRITE( Message,'(A,I0,A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex, ', pt=', i
       CALL Error( 'ElementMetricVec', Message )
       WRITE( Message,'(A,G10.3)') 'DetG:',DetJ(i)
       CALL Info( 'ElementMetricVec', Message, Level=3 )
       DO i=1,cdim
         WRITE( Message,'(A,I0,A,3G10.3)') 'Dir: ',i,' Coord:',Nodes % xyz(i,1),&
                 Nodes % xyz(i,2), Nodes % xyz(i,3)
         CALL Info( 'ElementMetricVec', Message, Level=3 )
       END DO
       IF (cdim < dim) THEN
         WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
         CALL Info( 'ElementMetricVec', Message, Level=3 )
       END IF
     END IF

   CONTAINS

     FUNCTION GetSymmetricIndex(i,j) RESULT(utind)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: i, j
       INTEGER :: utind

       IF (i>j) THEN
         utind = i*(i-1)/2+j
       ELSE
         utind = j*(j-1)/2+i
       END IF
     END FUNCTION GetSymmetricIndex
!------------------------------------------------------------------------------
   END FUNCTION ElementMetricVec
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>    Given element structure return value of the first partial derivatives with
!>    respect to global coordinates of a quantity x given at element nodes at
!>    local coordinate point u,v,w inside the element. Element basis functions
!>    are used to compute the value. This is internal version, and shouldn't
!>    usually be called directly by the user, but through the wrapper routine
!>    GlobalFirstDerivatives.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalFirstDerivativesInternal( elm,nodes,df,gx,gy,gz, &
                       Metric,dLBasisdx )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: element
!      INPUT: element structure
!
!    Type(Nodes_t) :: nodes
!      INPUT: element nodal coordinate arrays
!     
!     REAL(KIND=dp) :: f(:)
!      INPUT: Nodal values of the quantity whose partial derivative we want to know
!
!     REAL(KIND=dp) :: gx = @f(u,v)/@x, gy = @f(u,v)/@y, gz = @f(u,v)/@z
!      OUTPUT: Values of the partial derivatives
!
!     REAL(KIND=dp) :: Metric(:,:)
!      INPUT: Contravariant metric tensor of the element coordinate system
!
!     REAL(KIND=dp), OPTIONAL :: dLBasisdx(:,:)
!      INPUT: Values of partial derivatives with respect to local coordinates
!
!   FUNCTION VALUE:
!      .TRUE. if element is ok, .FALSE. if degenerated
!
!------------------------------------------------------------------------------
   !
   ! Return value of first derivatives of a quantity f in global
   ! coordinates at point (u,v) in gx,gy and gz.
   !
     TYPE(Element_t) :: elm
     TYPE(Nodes_t) :: nodes
 
     REAL(KIND=dp) :: df(:),Metric(:,:)
     REAL(KIND=dp) :: gx,gy,gz
     REAL(KIND=dp) :: dLBasisdx(:,:)

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
     REAL(KIND=dp) :: dx(3,3),dfc(3),s

     INTEGER :: cdim,dim,i,j,n,NB
!------------------------------------------------------------------------------

     n    = elm % TYPE % NumberOfNodes
     dim  = elm % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     x => nodes % x
     y => nodes % y
     z => nodes % z
!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local, and
!    partial derivatives of the quantity given, also with respect to local
!    coordinates
!------------------------------------------------------------------------------
     SELECT CASE(cdim)
       CASE(1)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
         END DO

       CASE(2)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
            dx(2,i) = SUM( y(1:n)*dLBasisdx(1:n,i) )
         END DO

       CASE(3)
         DO i=1,dim
            dx(1,i) = SUM( x(1:n)*dLBasisdx(1:n,i) )
            dx(2,i) = SUM( y(1:n)*dLBasisdx(1:n,i) )
            dx(3,i) = SUM( z(1:n)*dLBasisdx(1:n,i) )
         END DO
     END SELECT
!------------------------------------------------------------------------------
!    Contravariant components of partials in element coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       s = 0.0d0
       DO j=1,dim
         s = s + Metric(i,j) * df(j)
       END DO
       dfc(i) = s
     END DO
!------------------------------------------------------------------------------
!    Transform partials to space coordinates
!------------------------------------------------------------------------------
     gx = 0.0d0
     gy = 0.0d0
     gz = 0.0d0
     SELECT CASE(cdim)
       CASE(1)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )

       CASE(2)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )
         gy = SUM( dx(2,1:dim) * dfc(1:dim) )

       CASE(3)
         gx = SUM( dx(1,1:dim) * dfc(1:dim) )
         gy = SUM( dx(2,1:dim) * dfc(1:dim) )
         gz = SUM( dx(3,1:dim) * dfc(1:dim) )
     END SELECT

   END SUBROUTINE GlobalFirstDerivativesInternal
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of the first partial derivative with
!>   respect to global coordinates of a quantity f given at element nodes at
!>   local coordinate point u,v,w inside the element. Element basis functions
!>   are used to compute the value.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalFirstDerivatives( Elm, Nodes, df, gx, gy, gz, &
                    Metric, dLBasisdx )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!   Type(Element_t) :: element
!     INPUT: element structure
!
!   Type(Nodes_t) :: nodes
!     INPUT: element nodal coordinate arrays
!     
!   REAL(KIND=dp) :: f(:)
!     INPUT: Nodal values of the quantity whose partial derivatives we want
!            to know
!
!   REAL(KIND=dp) :: gx=@f(u,v,w)/@x, gy=@f(u,v,w)/@y, gz=@f(u,v,w)/@z
!     OUTPUT: Values of the partial derivatives
!
!   REAL(KIND=dp) :: u,v,w
!     INPUT: Point at which to evaluate the partial derivative
!
!   REAL(KIND=dp)L :: dLBasisdx(:,:)
!     INPUT: Values of partial derivatives of basis functions with respect to
!            local coordinates
!
!   REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
!     INPUT: Values of partial derivatives of basis functions with respect to
!            global coordinates can be given here, if known, otherwise they
!            will be computed from the element basis functions.
!
!------------------------------------------------------------------------------

     TYPE(Element_t) :: elm
     TYPE(Nodes_t) :: nodes

     REAL(KIND=dp) :: gx,gy,gz
     REAL(KIND=dp) :: dLBasisdx(:,:),Metric(:,:),df(:)

!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: n
!------------------------------------------------------------------------------

    CALL GlobalFirstDerivativesInternal( Elm, Nodes, df, &
              gx, gy, gz, Metric, dLBasisdx )

   END SUBROUTINE GlobalFirstDerivatives
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>   Given element structure return value of a quantity x given at element nodes
!>   at local coordinate point u inside the element. Element basis functions are
!>   used to compute the value. This is just a wrapper routine and will call the
!>   real function according to element dimension.   
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>          Compute elementwise matrix of second partial derivatives
!>          at given point u,v,w in global coordinates.
!------------------------------------------------------------------------------
   SUBROUTINE GlobalSecondDerivatives(elm,nodes,values,u,v,w,Metric,&
                     dBasisdx,ddLBasisddx,nd)
!------------------------------------------------------------------------------
!  
!       Parameters:
!  
!           Input:   (Element_t) structure describing the element
!                    (Nodes_t)   element nodal coordinates
!                    (double precision) F nodal values of the quantity
!                    (double precision) u,v point at which to evaluate
!  
!           Output:   3x3 matrix (values) of partial derivatives
!  
!------------------------------------------------------------------------------

     TYPE(Nodes_t)   :: nodes
     TYPE(Element_t) :: elm

     INTEGER :: nd
 
     REAL(KIND=dp) :: u,v,w
     REAL(KIND=dp) ::  Metric(:,:)
     REAL(KIND=dp) ::  values(:,:,:)
     REAL(KIND=dp) :: dBasisdx(:,:), ddLBasisddx(:,:,:)
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,l,n,q,dim,cdim

     REAL(KIND=dp), DIMENSION(3,3,3) :: C1,C2,ddx
     REAL(KIND=dp) :: df(3), cddf(3,3),ddf(3,3),dx(3,3)

     REAL(KIND=dp) :: s
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
!------------------------------------------------------------------------------
#if 0
#if 1
!
! This is actually not quite correct...
!
     IF ( elm % TYPE % BasisFunctionDegree <= 1 ) RETURN
#else
!
! this is ...
!
     IF ( elm % TYPE % ElementCode <= 202 .OR. &
          elm % TYPE % ElementCode == 303 .OR. &
          elm % TYPE % ElementCode == 504 ) RETURN
#endif
#endif

     n  = elm % TYPE % NumberOfNodes
     x => nodes % x
     y => nodes % y
     z => nodes % z

     dim  = elm % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()


!------------------------------------------------------------------------------
!    Partial derivatives of the basis functions are given, just
!    sum for the first partial derivatives...
!------------------------------------------------------------------------------
     dx = 0.0d0
     SELECT CASE( cdim )
       CASE(1)
         DO i=1,dim
           dx(1,i) = SUM( x(1:nd)*dBasisdx(1:nd,i) )
         END DO

       CASE(2)
         DO i=1,dim
           dx(1,i) = SUM( x(1:nd)*dBasisdx(1:nd,i) )
           dx(2,i) = SUM( y(1:nd)*dBasisdx(1:nd,i) )
         END DO

       CASE(3)
         DO i=1,dim
           dx(1,i) = SUM( x(1:nd)*dBasisdx(1:nd,i) )
           dx(2,i) = SUM( y(1:nd)*dBasisdx(1:nd,i) )
           dx(3,i) = SUM( z(1:nd)*dBasisdx(1:nd,i) )
         END DO
     END SELECT
!------------------------------------------------------------------------------
!     Get second partial derivatives with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
       DO j=1,dim
         ddx(1,i,j) = SUM(ddLBasisddx(1:nd,i,j)*x(1:nd) )
         ddx(2,i,j) = SUM(ddLBasisddx(1:nd,i,j)*y(1:nd) )
         ddx(3,i,j) = SUM(ddLBasisddx(1:nd,i,j)*z(1:nd) )
       END DO
     END DO
!
!------------------------------------------------------------------------------
!    Christoffel symbols of the second kind of the element coordinate system
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          DO k=1,dim
            s = 0.0d0
            DO l=1,cdim
              s = s + ddx(l,i,j)*dx(l,k)
            END DO
            C2(i,j,k) = s
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
!    Christoffel symbols of the first kind
!------------------------------------------------------------------------------
      DO i=1,dim
        DO j=1,dim
          DO k=1,dim
            s = 0.0d0
            DO l=1,dim
              s = s + Metric(k,l)*C2(i,j,l)
            END DO
            C1(i,j,k) = s
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
!     First add ordinary partials (change of the quantity with coordinates)...
!------------------------------------------------------------------------------
     Values = 0.0d0
     DO q=1,nd
       df  = dBasisdx(q,:)
       ddf = ddLBasisddx(q,:,:)

!------------------------------------------------------------------------------
!     ... then add change of coordinates
!------------------------------------------------------------------------------
        DO i=1,dim
          DO j=1,dim
            s = 0.0d0
            DO k=1,dim
              s = s - C1(i,j,k)*df(k)
            END DO
            ddf(i,j) = ddf(i,j) + s
          END DO
        END DO
!------------------------------------------------------------------------------
!       Convert to contravariant base
!------------------------------------------------------------------------------
        DO i=1,dim
          DO j=1,dim
            s = 0.0d0
            DO k=1,dim
              DO l=1,dim
                s = s + Metric(i,k)*Metric(j,l)*ddf(k,l)
              END DO
            END DO
            cddf(i,j) = s
          END DO
        END DO
!------------------------------------------------------------------------------
!      And finally transform to global coordinates 
!------------------------------------------------------------------------------
        DO i=1,cdim
          DO j=1,cdim
            s = 0.0d0
            DO k=1,dim
              DO l=1,dim
                s = s + dx(i,k)*dx(j,l)*cddf(k,l)    
              END DO
            END DO
            Values(q,i,j) = s
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE GlobalSecondDerivatives

END MODULE ElemInfo
