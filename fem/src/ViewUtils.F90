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

     USE MeshUtils

CONTAINS

!------------------------------------------------------------------------------
  ! Find planar ares and reduce those, if found, to fewer elements
!------------------------------------------------------------------------------
  FUNCTION PlanarReduce( n, Normals, Coord, Mesh ) RESULT(MeshOut)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh, MeshOut
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
      BLOCK
        REAL(KIND=dp) ::  c(3),r,eps=1d-8
        INTEGER :: i2,j,k,l,Usedn2, Setn2, np2, pn2, ind(32,2), ind2(32,2)

        TYPE(Element_t), POINTER :: el, ed

        LOGICAL, ALLOCATABLE :: Used2(:)
        INTEGER, ALLOCATABLE :: Set2(:)
        REAL(KIND=dp), ALLOCATABLE :: DirVec2(:)

        ALLOCATE(Set2(2*Setn), Used2(2*Setn), DirVec2(2*3*Setn))
       
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

          ! normal vector of  the edge in the plane of the parent element
          BLOCK
          REAL(KIND=dp) :: rr(3,3), d(3), e(3)

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
            m = 3*(el % NodeIndexes(l)-1)+1
            c = c + Coord(m:m+2)
          END DO

          ! edge centerpoint
          e = 0
          DO l=1,2
            m = 3*(ed % NodeIndexes(l)-1)+1
            e = e + Coord(m:m+2)
          END DO
          ! direction vector from edge center to parent center 
          c = c/3 - e/2

          ! "outer" normal
          IF ( SUM(d*c)>0 ) d=-d

          k = 3*(np2-1)+1
          DirVec2(k:k+2) = d
          END BLOCK
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
        ! ----

        CALL FindEdges0(Mesh2)

!       DO i2=1,Mesh2 % NumberOfBulkElements
!         k = 3*(Mesh2 % Elements(i2) % NodeIndexes(1)-1)+1
!         l = 3*(Mesh2 % Elements(i2) % NodeIndexes(2)-1)+1
!         c = ABS( Coord(k:k+2) - Coord(l:l+2) )
!         r = SQRT(SUM(c**2))
!         DirVec2(3*(i2-1)+1:3*(i2-1)+3)  = c / r
!       END DO

        ! ... edge lines ...
        pn2 = 0
        i2 = 1
        Usedn2 = 0;
        Used2  = .FALSE.

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

        problems = .TRUE.
        IF ( m==pn2 .AND. ind2(m,2)==ind2(1,1) ) &
                 problems = .NOT. IsConvex(Coord,Ind2(:,1), m)

        IF ( problems ) THEN
          PRINT*,'no luck, maybe hole(s) in geometry, polygon not convex, or something? Using the orignal mesh.'
          DO i=1,nn
            DEALLOCATE(MeshOut % Elements(i) % NodeIndexes)
          END DO
          DEALLOCATE(MeshOut % Elements)
          DEALLOCATE(MeshOut)
          MeshOut => Null()
          RETURN
        END IF

        ! ----
        DO j=1,Mesh2 % NumberOfBulkElements
          DEALLOCATE(Mesh2 % Elements(j) % NodeIndexes)
          DEALLOCATE(Mesh2 % Elements(j) % EdgeIndexes)
        END DO
        DEALLOCATE(Mesh2 % Edges)
        DEALLOCATE(Mesh2 % Elements)
        ! ----

        IF ( pn2 == 4 ) THEN
          !  construct one quad
          nn = nn + 1

          MeshOut % Elements(nn) % Type => GetElementType(404)
          ALLOCATE(MeshOut % Elements(nn) % NodeIndexes(4))
          MeshOut % Elements(nn) % NodeIndexes = ind2(1:4,1)
        ELSE

          ! construct a single circle or pn2-2 triangles (assumes convex area). ...
          BLOCK
            REAL(KIND=dp) :: cx, cy, cz, r, c(3),d(3)
            LOGICAL :: Circle

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
            DO i=2,pn2
              j = 3*(ind2(i,1)-1)
              r = 0
              r = r + (Coord(j+1) - cx)**2
              r = r + (Coord(j+2) - cy)**2
              r = r + (Coord(j+3) - cz)**2
              IF ( ABS(r-r0) > 1d-8 ) Circle = .FALSE.
            END DO
#else
            Circle = .FALSE.
#endif

            IF ( Circle ) THEN
              k=3*(ind2(1,1)-1)+1
              l=3*(ind2(2,1)-1)+1
              c = Coord(l:l+2) - Coord(k:k+2)

              l=3*(ind2(3,1)-1)+1
              d = Coord(l:l+2) - Coord(k:k+2)
 
              c = CrossProduct(c,d)
              c = c / SUM(SQRT(c**2))

              nn = nn + 1

              MeshOut % Elements(nn) % Type => GetElementType(101)

              ALLOCATE(MeshOut % Elements(nn) % PropertyData)
              MeshOut % Elements(nn) % PropertyData % Name = "circle"
              ALLOCATE(MeshOut % Elements(nn) % PropertyData % Values(8))

              MeshOut % Elements(nn) % PropertyData % Values(1) = 0        ! circle segment inner radius (not used)
              MeshOut % Elements(nn) % PropertyData % Values(2) = SQRT(r)  ! outer radius...
              MeshOut % Elements(nn) % PropertyData % Values(3) = cx       ! center point
              MeshOut % Elements(nn) % PropertyData % Values(4) = cy
              MeshOut % Elements(nn) % PropertyData % Values(5) = cz
              MeshOut % Elements(nn) % PropertyData % Values(6) = c(1)     ! normal vector
              MeshOut % Elements(nn) % PropertyData % Values(7) = c(2)
              MeshOut % Elements(nn) % PropertyData % Values(8) = c(3)
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
          END BLOCK
        END IF
      END BLOCK
      i = 1
    END DO

    MeshOut % NumberOfBulkElements = nn
    MeshOut % NumberOfBoundaryElements = 0

    !t0 = cputime() - t0
    !PRINT*,' Shadow elments,and time spent ', nn, t0

CONTAINS

  ! ...  find connected elements with equal "normals"
  RECURSIVE SUBROUTINE Traverse(n, i, Normals, Set, Setn, Used, Usedn, Mesh )
     INTEGER :: N, i, UsedN, Setn, Set(:)
     LOGICAL :: Used(:)
     REAL(KIND=dp) :: Normals(:)
     TYPE(Mesh_t), POINTER :: Mesh

     REAL(KIND=dp), PARAMETER :: eps=1d-8
     LOGICAL :: Eqn
     INTEGER :: j, l, m, e
     TYPE(Element_t), POINTER:: el, pel, ed

     IF ( Usedn >= n) RETURN

     el => Mesh % Elements(i)

     DO l = 1,el % Type % NumberOfEdges
       ed => Mesh % Edges(el % EdgeIndexes(l)) 

       pel => ed % BoundaryInfo %  Left
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
         Setn  = Setn+1
         Set(Setn) = j
         Used(j) = .TRUE.
         Usedn = Usedn+1
         CALL Traverse(N, j, Normals, Set, Setn, Used, Usedn, Mesh )
       END IF
     END DO
   END SUBROUTINE Traverse


   ! ...  find "edges" of the 1d line elements (used to traverse connected elements)
   SUBROUTINE FindEdges0(Mesh)
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i, j, maxi, n2
     TYPE(Element_t), POINTER :: el

     maxi = 0
     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % Type => GetElementType(202)
       Mesh % Elements(i) % ElementIndex = i
       maxi = MAX(maxi,Mesh % Elements(i) % NodeIndexes(1))
       maxi = MAX(maxi,Mesh % Elements(i) % NodeIndexes(2))
     END DO

     ALLOCATE( Mesh % Edges(maxi) ) 

     DO i=1,maxi
       ALLOCATE( mesh % edges(i) % BoundaryInfo )
       Mesh % Edges(i) % Type => GetElementType(101)
       Mesh % Edges(i) % BoundaryInfO % left => null()
       Mesh % Edges(i) % BoundaryInfO % right => null()
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

         stop 'k'

       END DO
     END DO
     Mesh % NumberOfEdges = maxi
   END SUBROUTINE FindEdges0

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
