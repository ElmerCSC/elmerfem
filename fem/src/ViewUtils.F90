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
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16 Apr 2024
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \}

#include "../config.h"
!------------------------------------------------------------------------------
!>  Utility routines for the Elmer main program.
!------------------------------------------------------------------------------
MODULE ViewUtils

     USE MeshBasics

CONTAINS

!------------------------------------------------------------------------------
  ! Find planar ares and reduce those, if found, to fewer elements
!------------------------------------------------------------------------------
  FUNCTION PlanarReduce( n, Normals, Coord, Mesh ) RESULT(MeshOut)
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    TYPE(Mesh_t), POINTER :: MeshOut
    INTEGER :: n
    REAL(KIND=dp) :: Normals(:), Coord(:)
!------------------------------------------------------------------------------
    INTEGER :: i, j, k, Usedn, Setn,  pn,nn
    LOGICAL, ALLOCATABLE :: Used(:)
    INTEGER, ALLOCATABLE :: Set(:), Ref(:)

    LOGICAL  :: problems

    REAL(KIND=dp) :: t0,t1,t2
    TYPE(Mesh_t), POINTER :: Mesh2
!------------------------------------------------------------------------------
    !t0 = cputime()

    CALL FindMeshEdges2D(Mesh)

    ! Elements may carry their original global ElementIndex (e.g. BulkElements+i).
    ! Traverse uses ElementIndex as an index into Used(n)/Set(n)/Normals(3*n), so
    ! it must equal the local position 1..n here.
    DO i=1,n
      Mesh % Elements(i) % ElementIndex = i
    END DO

    ALLOCATE( Set(n), Used(n), Ref(SIZE(Coord)/3) )
    ALLOCATE(Mesh2, MeshOut)
    Mesh2   = Mesh
    MeshOut = Mesh
    ALLOCATE(MeshOut % Elements(n))

    ! Find connected planar areas, if any ...
    i=1
    Used  = .FALSE.
    Usedn = 0;
    pn = 0
    nn = 0
    DO WHILE( Usedn < N )
      IF ( Used(i) ) THEN
        i=i+1; CYCLE
      END IF

      Used(i) = .TRUE.
      Usedn = Usedn+1

      Setn = 1
      Set(1) = i
      CALL Traverse( N, i, Normals, Set, Setn, Used, Usedn, Mesh )
      pn = pn + 1

      IF ( Setn<=2 ) THEN
         DO j=1,Setn
           nn = nn + 1

           MeshOut % Elements(nn) = Mesh % Elements(Set(j))
           k = Mesh % Elements(Set(j)) % Type % NumberOfNodes

           ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(k))
           MeshOut % Elements(nn) % NodeIndexes = &
                      Mesh % Elements(Set(j)) % NodeIndexes
         END DO
         CYCLE
      END IF


      ! Given one planar area, reduce it to one "element"
      CALL ReducePlanarSet()

      i = 1
    END DO

    IF(nn == Mesh % NumberOfBulkElements ) THEN
      CALL Info('PlanaReduce','No reduction in the element count! Using the orignal mesh!')
      DO i=1,nn
        DEALLOCATE(MeshOut % Elements(i) % NodeIndexes)
      END DO
      DEALLOCATE(MeshOut % Elements)
      DEALLOCATE(MeshOut)
      MeshOut => NULL()
      RETURN     
    END IF
    
    MeshOut % NumberOfBulkElements = nn
    MeshOut % NumberOfBoundaryElements = 0

    CALL Info('PlanarReduce','Reduced mesh number of elements: '//I2S(nn))
   
    !t0 = cputime() - t0
    !PRINT*,' Shadow elments,and time spent ', nn, t0

CONTAINS

  ! ...  find connected elements with equal "normals"
  RECURSIVE SUBROUTINE Traverse(n, i, Normals, Set, Setn, Used, Usedn, Mesh )
     INTEGER :: N, i, UsedN, Setn, Set(:)
     LOGICAL :: Used(:)
     REAL(KIND=dp) :: Normals(:)
     TYPE(Mesh_t) :: Mesh

     REAL(KIND=dp), PARAMETER :: eps=1d-4
     LOGICAL :: Eqn
     INTEGER :: j, l, m, e
     TYPE(Element_t), POINTER:: el, pel, ed

     IF ( Usedn >= n) RETURN

     el => Mesh % Elements(i)

     DO l = 1,el % Type % NumberOfEdges
       ed => Mesh % Edges(el % EdgeIndexes(l)) 

       pel => ed % BoundaryInfo % Left
       IF ( ASSOCIATED(pel, el) ) THEN
         pel => ed % BoundaryInfo % Right
       END IF
       IF (.NOT. ASSOCIATED(pel)) CYCLE

       j = pel % ElementIndex
       IF ( Used(j) ) CYCLE

       eqn = .TRUE.
       DO m=1,3
         eqn = eqn .AND. ABS(Normals(3*(i-1)+m) - Normals(3*(j-1)+m)) < eps
       END DO

       IF ( eqn ) THEN
         Setn  = Setn + 1
         Set(Setn) = j
         Used(j) = .TRUE.
         Usedn = Usedn + 1
         CALL Traverse(N, j, Normals, Set, Setn, Used, Usedn, Mesh )
         IF (Usedn>=n) EXIT
       END IF
     END DO
   END SUBROUTINE Traverse


   ! ...  find "edges" of the 1d line elements (used to traverse connected elements)
   SUBROUTINE FindEdges0(Mesh)
     TYPE(Mesh_t) :: Mesh
     INTEGER :: i, j, maxi, n2
     TYPE(Element_t), POINTER :: el

     maxi = 0
     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % Type => GetElementType(202)
       Mesh % Elements(i) % ElementIndex = i
       maxi = MAX(maxi,Mesh % Elements(i) % NodeIndexes(1))
       maxi = MAX(maxi,Mesh % Elements(i) % NodeIndexes(2))
     END DO

     ALLOCATE(Mesh % Edges(maxi)) 

     DO i=1,maxi
       ALLOCATE( Mesh % Edges(i) % BoundaryInfo )
       Mesh % Edges(i) % Type => GetElementType(101)
       Mesh % Edges(i) % BoundaryInfO % Left  => Null()
       Mesh % Edges(i) % BoundaryInfO % Right => Null()
     END DO

     DO i=1,Mesh % NumberOfBulkElements
       ALLOCATE( Mesh % Elements(i) % EdgeIndexes(2) )
       Mesh % Elements(i) % Type % NumberOfEdges = 2
       DO j=1,Mesh % Elements(i) % Type % NumberOfNodes
         Mesh % Elements(i) % EdgeIndexes(j) = Mesh % Elements(i) % NodeIndexes(j)

         el => Mesh % Edges(Mesh % Elements(i) % NodeIndexes(j))

         IF ( .NOT. ASSOCIATED( el % BoundaryInfo % Left ) ) THEN
           el % BoundaryInfo % Left => Mesh % Elements(i)
           CYCLE
         END IF

         IF ( .NOT. ASSOCIATED( el % BoundaryInfo % Right ) ) THEN
           el % BoundaryInfo % Right => Mesh % Elements(i)
           CYCLE
         END IF

         STOP 'k'
       END DO
     END DO
     Mesh % NumberOfEdges = maxi
   END SUBROUTINE FindEdges0

!------------------------------------------------------------------------------
!> Reduce the connected planar set Set(1:Setn) into as few output elements as
!> possible, and append them to MeshOut.
!>
!> This used to be an inline BLOCK construct with three further BLOCKs nested
!> inside it (and a GOTO jumping out of them). Intel Fortran 20.0 miscompiles
!> nested BLOCKs that declare allocatable/automatic arrays, so the whole thing
!> is a contained procedure instead; the mesh and bookkeeping variables come in
!> through host association exactly as they did through the BLOCK scoping.
!------------------------------------------------------------------------------
   SUBROUTINE ReducePlanarSet()
!------------------------------------------------------------------------------
     LOGICAL :: handled
     INTEGER :: i2, j, k, l, m, mm, Usedn2, Setn2, np2, pn2
     INTEGER :: ind(128,2), ind2(128,2)
     REAL(KIND=dp) :: c(3), d(3), e(3)

     TYPE(Element_t), POINTER :: el, ed

     LOGICAL, ALLOCATABLE :: Used2(:)
     INTEGER, ALLOCATABLE :: Set2(:)
     REAL(KIND=dp), ALLOCATABLE :: DirVec2(:)

     ! ... hole bridging (was a nested BLOCK)
     INTEGER :: ntri_h, ej_h
     INTEGER, ALLOCATABLE :: tris_h(:,:)

     ! ... ear clipping (was a nested BLOCK)
     INTEGER :: entri, ej
     INTEGER, ALLOCATABLE :: etris(:,:)
     REAL(KIND=dp) :: eplnorm(3)

     ! ... circle / fan triangulation (was a nested BLOCK)
     LOGICAL :: Circle
     REAL(KIND=dp) :: cx, cy, cz, r0, cc(3), dd(3)
#ifdef use_circle_detection
     REAL(KIND=dp) :: r
#endif
!------------------------------------------------------------------------------
     ALLOCATE(Set2(4*Setn), Used2(4*Setn), DirVec2(4*3*Setn))

     ALLOCATE( Mesh2 % Elements(n) )
     DO j=1,Setn
       Mesh2 % Elements(j) = Mesh % Elements(Set(j))
       Mesh2 % Elements(j) % ElementIndex = j
       k = Mesh % Elements(Set(j)) % Type % NumberOfNodes
       ALLOCATE(Mesh2 % Elements(j) % NodeIndexes(k))
       Mesh2 % Elements(j) % NodeIndexes = Mesh % Elements(Set(j)) % NodeIndexes
     END DO
     Mesh2 % NumberOfBulkElements = Setn
     Mesh2 % NumberOfBoundaryElements = 0


     ! Find outer edges, then reduce edge lines to #n 1d-sets
     CALL FindMeshEdges2D(Mesh2)

     np2 = 0
     DO j=1,Mesh2 % NumberOfEdges
       IF ( ASSOCIATED(Mesh2 % Edges(j) % BoundaryInfo % Right) ) CYCLE

       ed => Mesh2 % Edges(j)
       el => ed % BoundaryInfo % Left
       np2 = np2 + 1

       ! normal vector of the edge in the plane of the parent element
       !
       ! outer edge direction
       k = 3*(ed % NodeIndexes(1)-1)+1
       l = 3*(ed % NodeIndexes(2)-1)+1
       c = Coord(k:k+2) - Coord(l:l+2)

       ! 1st parent  edge
       k = 3*(el % NodeIndexes(1)-1)+1
       l = 3*(el % NodeIndexes(2)-1)+1
       d = Coord(k:k+2) - Coord(l:l+2)

       ! 2nd parent edge
       l = 3*(el % NodeIndexes(3)-1)+1
       e = Coord(k:k+2) - Coord(l:l+2)

       ! d normal to parent
       d = CrossProduct(d,e)

       ! normal to edge in plane of the parent
       d = CrossProduct(c,d)
       d = d/SQRT(SUM(d**2))

       ! parent "center point"
       c = 0
       DO l=1,3
         mm = 3*(el % NodeIndexes(l)-1)+1
         c = c + Coord(mm:mm+2)
       END DO

       ! edge centerpoint
       e = 0
       DO l=1,2
         mm = 3*(ed % NodeIndexes(l)-1)+1
         e = e + Coord(mm:mm+2)
       END DO
       ! direction vector from edge center to parent center 
       c = c/3 - e/2

       ! "outer" normal
       IF ( SUM(d*c)>0 ) d=-d

       k = 3*(np2-1)+1
       DirVec2(k:k+2) = d
     END DO

     np2 = 0
     DO j=1,Mesh2 % NumberOfEdges
       IF ( ASSOCIATED(Mesh2 % Edges(j) % BoundaryInfo % Right) ) CYCLE

       np2 = np2 + 1
       ed => Mesh2 % Edges(j)
       el => Mesh2 % Elements(np2)
       IF (np2<=Setn) DEALLOCATE(el % NodeIndexes)

       el % ElementIndex = np2
       ALLOCATE(el % NodeIndexes(2))
       el % NodeIndexes = ed % NodeIndexes
     END DO
     Mesh2 % NumberOfBulkElements = np2

     ! -----
     DO j=1,Mesh2 % NumberOfEdges
       DEALLOCATE(Mesh2 % Edges(j) % NodeIndexes)
     END DO
     DEALLOCATE(Mesh2 % Edges)
     Mesh2 % Edges => NULL()
     ! ----

     CALL FindEdges0(Mesh2)

     ! ... edge lines ...
     pn2 = 0
     i2 = 1
     Usedn2 = 0;
     Used2  = .FALSE.
     np2 = Mesh2 % numberofbulkelements

     DO WHILE( Usedn2 < np2 )
       IF ( Used2(i2) ) THEN
         i2=i2+1; CYCLE
       END IF

       Used2(i2) = .TRUE.
       Usedn2 = Usedn2 + 1
       Setn2 = 1
       Set2(1) = i2

       CALL Traverse( np2, i2, DirVec2, Set2, Setn2, Used2, Usedn2, Mesh2 )
       pn2 = pn2 + 1

       ! ... pick exterme nodes of the 1d sets ...
       Ref = 0
       DO j=1,Setn2
         el => Mesh2 % Elements(Set2(j))
         Ref(el % NodeIndexes(1:2)) = Ref(el % NodeIndexes(1:2)) + 1
       END DO

       l = 0
       DO j=1,Setn2
         el => Mesh2 % Elements(Set2(j))
         DO k=1,2
           IF ( Ref(el % NodeIndexes(k))==1 ) THEN
              l = l + 1
              ind(pn2,l) = el % NodeIndexes(k)
           END IF
         END DO
         IF ( l>= 2 ) EXIT
       END DO

       i2 = 1
     END DO

     ! ... order the points to a closed (hopefully) polygonal shape
     ind2 = -1
     ind2(1,:) = ind(1,:)
     m = 1
     ind(1,:)  = 0
     DO k=2,pn2
       j = ind2(k-1,2)
       DO l=2,pn2
         IF ( ind(l,1) == j ) THEN
           ind2(k,:) = ind(l,:)
           ind(l,:)  = 0
           m = m + 1
           EXIT
         ELSE IF ( ind(l,2) == j ) THEN
           ind2(k,1) = ind(l,2)
           ind2(k,2) = ind(l,1)
           ind(l,:)  = 0
           m = m + 1
           EXIT
         END IF
       END DO
     END DO

     IF ( ind2(m,2) /= ind2(1,1) ) THEN
       ! Open chain = truly broken topology, cannot repair
       CALL Info('PlanarReduce','Could not construct superelement? Using original elements.',Level=10)
       DO j=1,Setn
         nn = nn + 1

         MeshOut % Elements(nn) = Mesh % Elements(Set(j)) 
         k = Mesh % Elements(Set(j)) % TYPE % NumberOfNodes
         MeshOut % Elements(nn) % TYPE => Mesh % Elements(Set(j)) % Type

         ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(k))
         MeshOut % Elements(nn) % NodeIndexes = &
             Mesh % Elements(Set(j)) % NodeIndexes
       END DO

     ELSE
       ! Multiple closed loops: outer boundary + hole(s)
       handled = .FALSE.
       IF ( m < pn2 ) THEN
         ALLOCATE(tris_h(3, pn2+32))
         eplnorm = Normals(3*(Set(1)-1)+1:3*(Set(1)-1)+3)
         CALL BridgeHolesAndTriangulate(Coord, eplnorm,&
                ind2(1:m,1), m, ind(1:pn2,:), pn2, &
                tris_h, ntri_h)
         DO ej_h = 1, ntri_h
           nn = nn + 1
           MeshOut % Elements(nn) % Type => GetElementType(303)
           ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(3))
           MeshOut % Elements(nn) % NodeIndexes(1:3) = tris_h(1:3,ej_h)
         END DO
         DEALLOCATE(tris_h)
         handled = .TRUE.
       END IF

       IF ( .NOT. handled ) THEN
         IF ( .NOT. IsConvex(Coord, ind2(1:pn2,1), pn2) ) THEN
           ! Non-convex simple polygon: triangulate by ear clipping
           ALLOCATE(etris(3,pn2))
           eplnorm = Normals(3*(Set(1)-1)+1:3*(Set(1)-1)+3)
           CALL EarClipTriangulate(Coord, ind2(1:pn2,1), pn2, eplnorm, etris, entri)
           DO ej=1,entri
             nn = nn + 1
             MeshOut % Elements(nn) % Type => GetElementType(303)
             ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(3))
             MeshOut % Elements(nn) % NodeIndexes(1:3) = etris(1:3,ej)
           END DO
           DEALLOCATE(etris)

         ELSE IF ( pn2 == 4 ) THEN
           !  construct one quad
           nn = nn + 1

           MeshOut % Elements(nn) % Type => GetElementType(404)
           ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(4))
           MeshOut % Elements(nn) % NodeIndexes = ind2(1:4,1)

         ELSE
           ! construct a single circle or pn2-2 triangles (assumes convex area). ...
           cx = SUM(Coord(3*(ind2(1:pn2,1)-1)+1))/pn2
           cy = SUM(Coord(3*(ind2(1:pn2,1)-1)+2))/pn2
           cz = SUM(Coord(3*(ind2(1:pn2,1)-1)+3))/pn2

           j = 3*(ind2(1,1)-1)
           r0 = 0
           r0 = r0 + (Coord(j+1) - cx)**2
           r0 = r0 + (Coord(j+2) - cy)**2
           r0 = r0 + (Coord(j+3) - cz)**2

#ifdef use_circle_detection
           Circle = .TRUE.
           DO j=2,pn2
             k = 3*(ind2(j,1)-1)
             r = 0
             r = r + (Coord(k+1) - cx)**2
             r = r + (Coord(k+2) - cy)**2
             r = r + (Coord(k+3) - cz)**2
             IF ( ABS(r-r0) > 1d-8 ) Circle = .FALSE.
           END DO
#else
           Circle = .FALSE.
#endif

           IF ( Circle ) THEN
             k=3*(ind2(1,1)-1)+1
             l=3*(ind2(2,1)-1)+1
             cc = Coord(l:l+2) - Coord(k:k+2)

             l=3*(ind2(3,1)-1)+1
             dd = Coord(l:l+2) - Coord(k:k+2)

             cc = CrossProduct(cc,dd)
             cc = cc / SUM(SQRT(cc**2))

             nn = nn + 1

             MeshOut % Elements(nn) % Type => GetElementType(101)

             ALLOCATE(MeshOut % Elements(nn) % PropertyData)
             MeshOut % Elements(nn) % PropertyData % Name = "circle"
             ALLOCATE(MeshOut % Elements(nn) % PropertyData % Values(8))

             MeshOut % Elements(nn) % PropertyData % Values(1) = 0        ! circle segment inner radius (not used)
             MeshOut % Elements(nn) % PropertyData % Values(2) = SQRT(r0) ! outer radius...
             MeshOut % Elements(nn) % PropertyData % Values(3) = cx       ! center point
             MeshOut % Elements(nn) % PropertyData % Values(4) = cy
             MeshOut % Elements(nn) % PropertyData % Values(5) = cz
             MeshOut % Elements(nn) % PropertyData % Values(6) = cc(1)    ! normal vector
             MeshOut % Elements(nn) % PropertyData % Values(7) = cc(2)
             MeshOut % Elements(nn) % PropertyData % Values(8) = cc(3)
           ELSE
             DO j=2,pn2-1
               nn = nn + 1

               ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(3))
               MeshOut % Elements(nn) % Type => GetElementType(303)
               MeshOut % Elements(nn) % NodeIndexes(1) = ind2(1,1)
               MeshOut % Elements(nn) % NodeIndexes(2) = ind2(j,1)
               MeshOut % Elements(nn) % NodeIndexes(3) = ind2(j+1,1)
             END DO
           END IF
         END IF
       END IF  ! handled
     END IF

     ! ----
     DO j=1,Mesh2 % NumberOfBulkElements
       DEALLOCATE(Mesh2 % Elements(j) % NodeIndexes)
       DEALLOCATE(Mesh2 % Elements(j) % EdgeIndexes)
     END DO
     DEALLOCATE(Mesh2 % Edges)
     DEALLOCATE(Mesh2 % Elements)
     ! ----
!------------------------------------------------------------------------------
   END SUBROUTINE ReducePlanarSet
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END FUNCTION PlanarReduce
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 FUNCTION IsConvex(Coord,iv,n) RESULT(cnvx)
    INTEGER :: n, iv(n)
    LOGICAL :: cnvx
    REAL(KIND=dp) :: Coord(:)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nrm(3),a(3),b(3),c(3),ln,tol=1.d-12
    INTEGER :: i,j,k,l,im,ip

    cnvx = .FALSE.

   ! Compute reference normal
   i = 3*(iv(1)-1)+1
   a = Coord(i:i+2)

   i = 3*(iv(2)-1)+1
   b = Coord(i:i+2)

   i = 3*(iv(3)-1)+1
   c = Coord(i:i+2)

   nrm = CrossProduct(b-a, c-b)
   ln = SUM(SQRT(nrm**2))
   IF (ln<tol) RETURN
   nrm = nrm / ln

   DO i = 1, n
    ip = MOD(i,n)+1
    im = MOD(i-2+n,n) + 1

    j = 3*(iv(i)-1)+1
    k = 3*(iv(im)-1)+1
    l = 3*(iv(ip)-1)+1

    a = Coord(j:j+2)  - Coord(k:k+2)
    b = Coord(l:l+2)  - Coord(j:j+2)
    c = CrossProduct(a, b)
    IF (SUM(c*nrm) < -tol) RETURN
   END DO

   cnvx = .TRUE.
!------------------------------------------------------------------------------
 END FUNCTION IsConvex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE BridgeHolesAndTriangulate(Coord, planeNorm, outerV, nOuter, &
       holeSegs, nSegs, tris, ntri)
!------------------------------------------------------------------------------
!   Collect hole loops from the non-zero entries of holeSegs, bridge each hole
!   to the outer polygon via a rightward-ray visibility cut, then ear-clip.
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(IN)  :: Coord(:), planeNorm(3)
    INTEGER,       INTENT(IN)  :: nOuter, outerV(nOuter)
    INTEGER,       INTENT(IN)  :: nSegs, holeSegs(nSegs,2)
    INTEGER,       INTENT(OUT) :: ntri, tris(3,*)
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MAXV = 64, MAXH = 16

    INTEGER       :: holePoly(MAXV, MAXH), holeLen(MAXH), nHoles
    INTEGER       :: merged(MAXV*2), nMerged
    INTEGER       :: curOuter(MAXV*2), nCurOuter
    REAL(KIND=dp) :: opx(MAXV*2), opy(MAXV*2)
    REAL(KIND=dp) :: hpxA(MAXV), hpyA(MAXV)
    REAL(KIND=dp) :: uu(3), vv(3), ref(3), tmp3(3)
    INTEGER       :: ii, s, t, hh, hMax, bIdx, hlen, cur
    LOGICAL       :: segUsed(nSegs), found, reverseHole
    REAL(KIND=dp) :: hpxMax, hpyRay, xi, minXi, tpar, denom
    REAL(KIND=dp) :: x1, y1, x2, y2, oArea, hArea
!------------------------------------------------------------------------------
    ! Build 2D coordinate frame in the plane
    ref = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    IF (ABS(planeNorm(1)) > 0.9_dp) ref = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    uu = ref - SUM(ref*planeNorm)*planeNorm
    uu = uu / SQRT(SUM(uu**2))
    vv = CrossProduct(planeNorm, uu)

    ! Collect hole loops from non-zero (unused) segments
    segUsed(1:nSegs) = (holeSegs(1:nSegs,1) == 0)
    nHoles = 0
    DO s = 1, nSegs
      IF (segUsed(s)) CYCLE
      IF (nHoles >= MAXH) EXIT
      nHoles = nHoles + 1
      segUsed(s) = .TRUE.
      holePoly(1, nHoles) = holeSegs(s,1)
      holePoly(2, nHoles) = holeSegs(s,2)
      hlen = 2
      cur  = holeSegs(s,2)
      DO WHILE (cur /= holePoly(1,nHoles) .AND. hlen < MAXV)
        found = .FALSE.
        DO t = 1, nSegs
          IF (segUsed(t)) CYCLE
          IF (holeSegs(t,1) == cur) THEN
            segUsed(t) = .TRUE.; cur = holeSegs(t,2)
            hlen = hlen + 1; holePoly(hlen, nHoles) = cur
            found = .TRUE.; EXIT
          ELSE IF (holeSegs(t,2) == cur) THEN
            segUsed(t) = .TRUE.; cur = holeSegs(t,1)
            hlen = hlen + 1; holePoly(hlen, nHoles) = cur
            found = .TRUE.; EXIT
          END IF
        END DO
        IF (.NOT. found) EXIT
      END DO
      holeLen(nHoles) = hlen - 1   ! last vertex == first; don't repeat
    END DO

    ! Bridge each hole into the growing outer polygon
    curOuter(1:nOuter) = outerV(1:nOuter)
    nCurOuter = nOuter

    DO hh = 1, nHoles
      ! Project current outer polygon to 2D
      DO ii = 1, nCurOuter
        s = 3*(curOuter(ii)-1)+1
        tmp3 = Coord(s:s+2)
        opx(ii) = SUM(tmp3*uu); opy(ii) = SUM(tmp3*vv)
      END DO

      ! Project hole vertices to 2D
      DO ii = 1, holeLen(hh)
        s = 3*(holePoly(ii,hh)-1)+1
        tmp3 = Coord(s:s+2)
        hpxA(ii) = SUM(tmp3*uu); hpyA(ii) = SUM(tmp3*vv)
      END DO

      ! Signed areas to detect winding direction
      oArea = 0.0_dp
      DO ii = 1, nCurOuter
        t = MOD(ii, nCurOuter) + 1
        oArea = oArea + opx(ii)*opy(t) - opx(t)*opy(ii)
      END DO
      hArea = 0.0_dp
      DO ii = 1, holeLen(hh)
        t = MOD(ii, holeLen(hh)) + 1
        hArea = hArea + hpxA(ii)*hpyA(t) - hpxA(t)*hpyA(ii)
      END DO
      ! Holes should wind opposite to the outer polygon; reverse if same
      reverseHole = (oArea * hArea > 0.0_dp)

      ! Rightmost hole vertex — anchor for the rightward-ray bridge
      hMax = 1
      DO ii = 2, holeLen(hh)
        IF (hpxA(ii) > hpxA(hMax)) hMax = ii
      END DO
      hpxMax = hpxA(hMax); hpyRay = hpyA(hMax)

      ! Find the nearest outer edge intersected by the rightward ray from hMax
      bIdx  = 1
      minXi = HUGE(1.0_dp)
      DO ii = 1, nCurOuter
        t  = MOD(ii, nCurOuter) + 1
        x1=opx(ii); y1=opy(ii); x2=opx(t); y2=opy(t)
        denom = y2 - y1
        IF (ABS(denom) < 1.0d-14) CYCLE
        tpar = (hpyRay - y1) / denom
        IF (tpar < 0.0_dp .OR. tpar > 1.0_dp) CYCLE
        xi = x1 + tpar*(x2 - x1)
        IF (xi < hpxMax - 1.0d-12) CYCLE
        IF (xi < minXi) THEN
          minXi = xi
          ! Bridge to the rightmost endpoint of this edge
          IF (opx(ii) >= opx(t)) THEN; bIdx = ii
          ELSE;                         bIdx = t
          END IF
        END IF
      END DO

      ! Merged polygon: outer(bIdx..bIdx-1), outer(bIdx) [bridge],
      !                 hole(hMax..hMax-+) [reversed if needed], hole(hMax) [bridge back]
      nMerged = 0
      DO ii = 0, nCurOuter - 1
        nMerged = nMerged + 1
        merged(nMerged) = curOuter(MOD(bIdx-1+ii, nCurOuter)+1)
      END DO
      nMerged = nMerged + 1
      merged(nMerged) = curOuter(bIdx)                   ! bridge: repeat bIdx

      DO ii = 0, holeLen(hh) - 1
        nMerged = nMerged + 1
        IF (reverseHole) THEN
          merged(nMerged) = holePoly( &
               MOD(hMax-1-ii+2*holeLen(hh), holeLen(hh))+1, hh)
        ELSE
          merged(nMerged) = holePoly( &
               MOD(hMax-1+ii, holeLen(hh))+1, hh)
        END IF
      END DO
      nMerged = nMerged + 1
      merged(nMerged) = holePoly(hMax, hh)               ! bridge back

      curOuter(1:nMerged) = merged(1:nMerged)
      nCurOuter = nMerged
    END DO

    ! Ear-clip the fully merged polygon
    CALL EarClipTriangulate(Coord, merged(1:nMerged), nMerged, planeNorm, tris, ntri)
!------------------------------------------------------------------------------
  END SUBROUTINE BridgeHolesAndTriangulate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE EarClipTriangulate(Coord, iv, nv, planeNorm, tris, ntri)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(IN)  :: Coord(:), planeNorm(3)
    INTEGER,       INTENT(IN)  :: nv, iv(nv)
    INTEGER,       INTENT(OUT) :: ntri, tris(3,*)
!------------------------------------------------------------------------------
    INTEGER       :: ii, a, b, c, m, cnt, steps
    INTEGER       :: nxt(nv), prv(nv)
    REAL(KIND=dp) :: px(nv), py(nv), uu(3), vv(3), ref(3)
    REAL(KIND=dp) :: ax, ay, bx, by, cx, cy, mx, my
    REAL(KIND=dp) :: cross, area2, d1, d2, d3
    LOGICAL       :: isEar, hasNeg, hasPos
!------------------------------------------------------------------------------
    ! Build 2D coordinate frame in the plane
    ref = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    IF (ABS(planeNorm(1)) > 0.9_dp) ref = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    uu = ref - SUM(ref*planeNorm)*planeNorm
    uu = uu / SQRT(SUM(uu**2))
    vv = CrossProduct(planeNorm, uu)

    ! Project polygon vertices to 2D
    DO ii = 1, nv
      m = 3*(iv(ii)-1)+1
      ref = Coord(m:m+2)
      px(ii) = SUM(ref*uu)
      py(ii) = SUM(ref*vv)
    END DO

    ! Signed area (Shoelace) to determine polygon winding
    area2 = 0.0_dp
    DO ii = 1, nv
      a = MOD(ii,nv)+1
      area2 = area2 + px(ii)*py(a) - px(a)*py(ii)
    END DO

    ! Circular doubly-linked list
    DO ii = 1, nv
      nxt(ii) = MOD(ii,nv)+1
      prv(ii) = MOD(ii-2+nv,nv)+1
    END DO

    ntri  = 0
    cnt   = nv
    ii    = 1
    steps = 0

    DO WHILE (cnt > 3)
      IF (steps > cnt) EXIT   ! degenerate polygon, give up

      a = prv(ii); b = ii; c = nxt(ii)
      ax=px(a); ay=py(a); bx=px(b); by=py(b); cx=px(c); cy=py(c)

      ! Triangle signed area (must match polygon winding to be an ear)
      cross = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax)

      IF ((area2 >= 0.0_dp .AND. cross >= 0.0_dp) .OR. &
          (area2 <  0.0_dp .AND. cross <= 0.0_dp)) THEN
        isEar = .TRUE.
        m = nxt(c)
        DO WHILE (m /= a)
          mx=px(m); my=py(m)
          ! Sign of each edge vs. query point
          d1 = (bx-ax)*(my-ay) - (by-ay)*(mx-ax)
          d2 = (cx-bx)*(my-by) - (cy-by)*(mx-bx)
          d3 = (ax-cx)*(my-cy) - (ay-cy)*(mx-cx)
          hasNeg = (d1 < 0.0_dp) .OR. (d2 < 0.0_dp) .OR. (d3 < 0.0_dp)
          hasPos = (d1 > 0.0_dp) .OR. (d2 > 0.0_dp) .OR. (d3 > 0.0_dp)
          ! Same sign (or zero) = inside or on boundary -> not an ear
          IF (.NOT.(hasNeg .AND. hasPos)) THEN
            isEar = .FALSE.; EXIT
          END IF
          m = nxt(m)
        END DO
      ELSE
        isEar = .FALSE.
      END IF

      IF (isEar) THEN
        ntri = ntri + 1
        tris(1,ntri) = iv(a); tris(2,ntri) = iv(b); tris(3,ntri) = iv(c)
        nxt(a) = c; prv(c) = a
        cnt   = cnt - 1
        steps = 0
        ii    = c
      ELSE
        steps = steps + 1
        ii    = nxt(ii)
      END IF
    END DO

    ! Last triangle
    IF (cnt == 3) THEN
      ntri = ntri + 1
      tris(1,ntri) = iv(prv(ii)); tris(2,ntri) = iv(ii); tris(3,ntri) = iv(nxt(ii))
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE EarClipTriangulate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 FUNCTION LoadShadowMesh(Name) RESULT(ShadowMesh)
!------------------------------------------------------------------------------
   CHARACTER(LEN=*) :: Name
   TYPE(Mesh_t), POINTER :: ShadowMesh

   INTEGER :: i, j, ncnt, id, body, code, node(8), iostat
   REAL(KIND=dp)  :: x,y,z
   CHARACTER(LEN=256) :: Line

   ShadowMesh => AllocateMesh()

   OPEN( 33, File=TRIM(name)//'/mesh.header', STATUS='OLD', IOSTAT=iostat )
   IF (iostat /= 0 ) CALL Error( "ViewUtils:", "Error reading shadow mesh.")

   print*,TRIM(Name)//'/mesh.header'
   READ(33, *) NofNodes, NofElements
   CLOSE(33)

   ALLOCATE( ShadowMesh % Nodes % x(NofNodes), &
             ShadowMesh % Nodes % y(NofNodes), &
             ShadowMesh % Nodes % z(NofNodes), &
             ShadowMesh % Elements(NofElements) )

   OPEN( 33, File=TRIM(name)//'/mesh.nodes', STATUS='OLD', IOSTAT=iostat  )
   IF (iostat /= 0 ) CALL Error( "ViewUtils:", "Error reading shadow mesh.")
   DO i=1, NofNodes
     READ(33, *) id, body, x, y, z
     ShadowMesh % Nodes % x(i) = x
     ShadowMesh % Nodes % y(i) = y
     ShadowMesh % Nodes % z(i) = z
   END DO
   CLOSE(33)
   ShadowMesh % NumberOfNodes = NofNodes

   OPEN( 33, File=TRIM(name)//'/mesh.elements', STATUS='OLD', IOSTAT=iostat )
   IF (iostat /= 0 ) CALL Error( "ViewUtils:", "Error reading shadow mesh.")
   DO i=1, NofElements
     READ(33, '(a)' ) line
     READ(line, *) id, body, code

     ncnt = code - 100 * (code / 100)
     READ(line, *) id, body, code, (node(j),j=1,ncnt)

     ALLOCATE( ShadowMesh % Elements(i) % NodeIndexes(ncnt) )
     ShadowMesh % Elements(i) % NodeIndexes = node(1:ncnt)
     ShadowMesh % Elements(i) % Bodyid = body
     ShadowMesh % Elements(i) % ElementIndex = i
     ShadowMesh % Elements(i) % Type => GetElementType(code)
   END DO
   ShadowMesh % NumberOfBulkElements = Nofelements
   CLOSE(33)

!------------------------------------------------------------------------------
 END FUNCTION LoadShadowMesh
!------------------------------------------------------------------------------


!> \} ElmerLib

END MODULE ViewUtils
