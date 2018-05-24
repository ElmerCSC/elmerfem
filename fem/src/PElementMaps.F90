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
! *  Authors: Mikko Byckling, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 23 Aug 2004
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!-------------------------------------------------------------------------------
!>  Module defining mappings for p elements. These include nodal points 
!>  contained by faces and edges, element boundary maps (edges for 2d elements,
!>  faces for 3d) and mappings from faces to edge numbers. Mappings defined in 
!>  this module are compatible with basis functions defined in module 
!>  PElementBase.
!-------------------------------------------------------------------------------

MODULE PElementMaps
  USE Types

  IMPLICIT NONE

  ! Private mappings. For access use get[Element][Type]Map(i)
  PRIVATE QuadEdgeMap, TriangleEdgeMap, &
       TetraEdgeMap1, TetraFaceMap1, TetraFaceEdgeMap1, &
       TetraEdgeMap2, TetraFaceMap2, TetraFaceEdgeMap2, &
       BrickEdgeMap, BrickFaceMap, BrickFaceEdgeMap, &
       WedgeEdgeMap, WedgeFaceMap, WedgeFaceEdgeMap, &
       PyramidEdgeMap, PyramidFaceMap, PyramidFaceEdgeMap, &
       MInit
  ! Mappings
  INTEGER, TARGET, SAVE :: QuadEdgeMap(4,2), TriangleEdgeMap(3,2), &
       TetraEdgeMap1(6,2), TetraFaceMap1(4,3), TetraFaceEdgeMap1(4,3), &
       TetraEdgeMap2(6,2), TetraFaceMap2(4,3), TetraFaceEdgeMap2(4,3),&
       BrickEdgeMap(12,2), BrickFaceMap(6,4), BrickFaceEdgeMap(6,4), &
       WedgeEdgeMap(9,2), WedgeFaceMap(5,4), WedgeFaceEdgeMap(5,4), &
       PyramidEdgeMap(8,2), PyramidFaceMap(5,4), PyramidFaceEdgeMap(5,4)

  LOGICAL, SAVE :: MInit = .FALSE.
  !$OMP THREADPRIVATE(MInit, QuadEdgeMap, TriangleEdgeMap, &
  !$OMP&              TetraEdgeMap1, TetraFaceMap1, TetraFaceEdgeMap1, &
  !$OMP&              TetraEdgeMap2, TetraFaceMap2, TetraFaceEdgeMap2,&
  !$OMP&              BrickEdgeMap, BrickFaceMap, BrickFaceEdgeMap, &
  !$OMP&              WedgeEdgeMap, WedgeFaceMap, WedgeFaceEdgeMap, &
  !$OMP&              PyramidEdgeMap, PyramidFaceMap, PyramidFaceEdgeMap)
CONTAINS

  ! MAPPINGS

  ! First some direct mappings to elements. These should not be used directly 
  ! unless element type is implicitly known from context. Better way is to use
  ! getElement[Boundary,Edge,Face]Map -routines.

    ! Call: localEdge = getQuadEdge(i)
    !
    ! Function returns mapping from edge number to edge endpoints 

    FUNCTION getQuadEdgeMap(i) RESULT(localEdge)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(2) :: localEdge

      IF (.NOT. MInit) CALL InitializeMappings()
      
      localEdge(:) = QuadEdgeMap(i,:)
    END FUNCTION getQuadEdgeMap

    ! Call: localEdge = getTriangleEdge(i)
    ! 
    ! Function returns mapping from edge number to edge endpoints

    FUNCTION getTriangleEdgeMap(i) RESULT(localEdge)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(2) :: localEdge

      IF (.NOT. MInit) CALL InitializeMappings()

      localEdge(:)=TriangleEdgeMap(i,:)
    END FUNCTION getTriangleEdgeMap
    
    ! Call: localEdge = getBrickEdgeMap(i)
    ! 
    ! Function returns mapping from edge number to edge endpoints

    FUNCTION getBrickEdgeMap(i) RESULT(localEdge)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(2) :: localEdge

      IF (.NOT. MInit) CALL InitializeMappings()

      localEdge(:) = BrickEdgeMap(i,:)
    END FUNCTION getBrickEdgeMap
    
    ! Call: localFace = getBrickFaceMap(i)
    ! 
    ! Function returns mapping from face number to face nodes

    FUNCTION getBrickFaceMap(i) RESULT(localFace)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(4) :: localFace

      IF (.NOT. MInit) CALL InitializeMappings()

      localFace(:) = BrickFaceMap(i,:)
    END FUNCTION getBrickFaceMap

    ! Call: localEdge = getFaceEdgeMap(face, localNode)
    !
    ! getFaceEdgeMap returns number of local edge when given face and 
    ! its local node number. Node number is treated as edges beginning point

    FUNCTION getBrickFaceEdgeMap(face, localNode) RESULT(localEdge)
      IMPLICIT NONE

      ! Parameters 
      INTEGER, INTENT(IN) :: face, localNode
      ! Variables
      INTEGER :: localEdge

      IF (.NOT. MInit) CALL InitializeMappings()

      localEdge = BrickFaceEdgeMap(face,localNode)

      IF (localEdge == 0) THEN
         WRITE (*,'(A,I2,I3)') 'Unknown combination node for (face,node)', face,localNode 
         STOP
      END IF
    END FUNCTION getBrickFaceEdgeMap


    FUNCTION getTetraEdgeMap(i,TYPE) RESULT(edge)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN), OPTIONAL :: TYPE
      INTEGER :: t
      INTEGER, DIMENSION(2) :: edge

      IF (.NOT. MInit) CALL InitializeMappings()

      ! If type not present use default (1)
      t = 1
      IF (PRESENT(TYPE)) t = TYPE

      ! Select edge map by tetra type
      SELECT CASE (t)
      CASE (1)
         edge(:) = TetraEdgeMap1(i,:)
      CASE (2)
         edge(:) = TetraEdgeMap2(i,:)
      CASE DEFAULT
         CALL Fatal('PElementMaps::getTetraEdgeMap','Unknown tetra type')
      END SELECT
    END FUNCTION getTetraEdgeMap


    FUNCTION getTetraFaceMap(i,TYPE) RESULT(face)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN), OPTIONAL :: TYPE
      INTEGER :: t
      INTEGER, DIMENSION(3) :: face

      IF (.NOT. MInit) CALL InitializeMappings()
      
      ! If type not present use default (1)
      t = 1
      IF (PRESENT(TYPE)) t = TYPE

      ! Select face map by tetra type 
      SELECT CASE(t)
      CASE (1)
         face(:) = TetraFaceMap1(i,:)
      CASE (2)
         face(:) = TetraFaceMap2(i,:)
      CASE DEFAULT 
         CALL Fatal('PElementMaps::getTetraFaceMap','Unknown tetra type')
      END SELECT
    END FUNCTION getTetraFaceMap

    FUNCTION getWedgeEdgeMap(i) RESULT(edge)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(2) :: edge

      IF (.NOT. MInit) CALL InitializeMappings()

      edge(:) = WedgeEdgeMap(i,:)
    END FUNCTION getWedgeEdgeMap


    FUNCTION getWedgeFaceMap(i) RESULT(face)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(4) :: face

      IF (.NOT. MInit) CALL InitializeMappings()

      face(:) = WedgeFaceMap(i,:)
    END FUNCTION getWedgeFaceMap


    FUNCTION getPyramidEdgeMap(i) RESULT(edge)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(2) :: edge

      IF (.NOT. MInit) CALL InitializeMappings()

      edge(:) = PyramidEdgeMap(i,:)
    END FUNCTION getPyramidEdgeMap


    FUNCTION getPyramidFaceMap(i) RESULT(face)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: i
      INTEGER, DIMENSION(4) :: face

      IF (.NOT. MInit) CALL InitializeMappings()

      face(:) = PyramidFaceMap(i,:)
    END FUNCTION getPyramidFaceMap


!------------------------------------------------------------------------------
!>     Mapping from element local edge or face number to nodes contained in 
!>     that edge or face. 
!------------------------------------------------------------------------------
    FUNCTION getElementBoundaryMap(Element, i) RESULT(map)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get map for
!
!    INTEGER, INTENT(IN) :: i
!      INPUT: Local number of elements edge or face
!
!  FUNCTION VALUE:
!    INTEGER :: map(4)
!       Map containing local node numbers of given local edge or face
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t) :: Element
      INTEGER, INTENT(IN) :: i
      
      INTEGER :: map(4)

      IF (.NOT. MInit) CALL InitializeMappings()
      map = 0

      ! Function is not defined for non p elements
      !IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
      !   CALL Warn('PElementMaps::getElementBoundaryMap','Element not p element') 
      !   RETURN
      !END IF

      SELECT CASE(Element % TYPE % ElementCode / 100)
      CASE (3)
         map(1:2) = getTriangleEdgeMap(i)
      CASE (4)
         map(1:2) = getQuadEdgeMap(i)
      CASE (5)
         map(1:3) = getTetraFaceMap(i,Element % PDefs % TetraType)
      CASE (6)
         map(1:4) = getPyramidFaceMap(i)
      CASE (7)
         map(1:4) = getWedgeFaceMap(i)
      CASE (8)
         map(1:4) = getBrickFaceMap(i)
      CASE DEFAULT
         CALL Fatal('PElementMaps::getElementBoundaryMap','Unsupported element type')
      END SELECT
    END FUNCTION getElementBoundaryMap


!------------------------------------------------------------------------------
!>     Mapping from element local face to local edges contained in face. Given
!>     element and local face number this routine returns numbers of local edges
!>     on face. 
!------------------------------------------------------------------------------
    FUNCTION getFaceEdgeMap( Element, i) RESULT(map)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get map for
!
!    INTEGER, INTENT(IN) :: i
!      INPUT: Local number of element face
!
!  FUNCTION VALUE:
!    INTEGER :: map(4)
!       Map containing local numbers of edges on face
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t) :: Element
      INTEGER, INTENT(IN) :: i

      INTEGER :: elementCode, map(4)

      elementCode = Element % TYPE % ElementCode

      IF (.NOT. MInit) CALL InitializeMappings()

      ! Function is not defined for non p elements
      !IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
      !   CALL Warn('PElementMaps::getFaceEdgeMap','Element not p element') 
      !   map = 0
      !   RETURN
      !END IF

      SELECT CASE(elementCode / 100)
      CASE (5)
         map = 0
         SELECT CASE (Element % PDefs % TetraType)
         CASE (1)
            map(1:3) = TetraFaceEdgeMap1(i,:)
         CASE (2)
            map(1:3) = TetraFaceEdgeMap2(i,:)
         CASE DEFAULT
            CALL Fatal('PElementMaps::getFaceEdgeMap','Unknown tetra type')
         END SELECT
      CASE (6)
         map(1:4) = PyramidFaceEdgeMap(i,:)
      CASE (7)
         map(1:4) = WedgeFaceEdgeMap(i,:)
      CASE (8)
         map(1:4) = BrickFaceEdgeMap(i,:)
      CASE DEFAULT
         CALL Fatal('PElementMaps::getFaceEdgeMap','Unsupported element type')
      END SELECT
    END FUNCTION getFaceEdgeMap


!------------------------------------------------------------------------------
!>     Get mappings for given element to element edges and their nodes. Given 
!>     element, this routine returns a map containing nodes (endpoints) of
!>     elements edges. 
!------------------------------------------------------------------------------
    SUBROUTINE GetElementEdgeMap( Element, map )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get map for
!
!    INTEGER :: map(:,:)
!       OUTPUT: Map containing local node numbers of local edges
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t) :: Element
      INTEGER,  POINTER :: map(:,:)

      IF (.NOT. MInit) CALL InitializeMappings()

      ! Function is not defined for non p elements
      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementMaps::GetElementEdgeMap','Element not p element') 
         map = 0
         RETURN
      END IF

      SELECT CASE (Element % TYPE % ElementCode / 100)
      CASE (3)
         map => TriangleEdgeMap
      CASE (4)
         map => QuadEdgeMap
      CASE (5)
         SELECT CASE( Element % PDefs % TetraType )
         CASE (1)
            map => TetraEdgeMap1
         CASE (2)
            map => TetraEdgeMap2
         CASE DEFAULT
            CALL Fatal('PElementMaps::GetElementEdgeMap','Unknown tetra type for p element')
         END SELECT
      CASE (6)
         map => PyramidEdgeMap
      CASE (7)
         map => WedgeEdgeMap
      CASE (8)
         map => BrickEdgeMap
      CASE DEFAULT
         CALL Fatal('PElementMaps::GetElementEdgeMap','Unsupported element type')
      END SELECT
    END SUBROUTINE GetElementEdgeMap
   

!------------------------------------------------------------------------------
!>     Get mappings for given element to element faces and their nodes. Given 
!>     element, this routine returns a map containing nodes (endpoints) of
!>     elements face. 
!------------------------------------------------------------------------------
    SUBROUTINE GetElementFaceMap( Element, faceMap )
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get map for
!
!    INTEGER :: map(:,:)
!       OUTPUT: Map containing local node numbers of local faces
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      TYPE(Element_t) :: Element
      INTEGER, POINTER :: faceMap(:,:)

      IF (.NOT. MInit) CALL InitializeMappings()

      ! Function is not defined for non p elements
      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementMaps::GetElementFaceMap','Element not p element') 
         NULLIFY(faceMap)
         RETURN
      END IF

      SELECT CASE (Element % TYPE % ElementCode / 100)
!      CASE (3)
!         facemap => TriangleEdgeMap
!      CASE (4)
!         facemap => QuadEdgeMap
      CASE (5)
         SELECT CASE( Element % PDefs % TetraType )
         CASE (1)
            faceMap => TetraFaceMap1
         CASE (2)
            faceMap => TetraFaceMap2
         CASE DEFAULT
            CALL Fatal('PElementMaps::GetElementFaceMap','Unknown tetra type for p element')
         END SELECT
      CASE (6)
         faceMap => PyramidFaceMap
      CASE (7)
         faceMap => WedgeFaceMap
      CASE (8)
         faceMap => BrickFaceMap
      CASE DEFAULT
         CALL Fatal('PElementMaps::GetElementFaceMap','Unsupported element type')
      END SELECT
    END SUBROUTINE GetElementFaceMap


!------------------------------------------------------------------------------    
!>     Get mappings for given element to elements faces and their edge. Given 
!>     element, this routine returns a map containing local edge numbers of
!>     elements faces. 
!------------------------------------------------------------------------------    
    SUBROUTINE GetElementFaceEdgeMap( Element, faceEdgeMap )
!------------------------------------------------------------------------------    
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to get map for
!
!    INTEGER :: map(:,:)
!       OUTPUT: Map containing local edge numbers of local faces
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t) :: Element
      INTEGER, POINTER :: faceEdgeMap(:,:)
      
      IF (.NOT. MInit) CALL InitializeMappings()

      ! Function is not defined for non p elements
      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementMaps::GetElementFaceEdgeMap','Element not p element') 
         NULLIFY(faceEdgeMap)
         RETURN
      END IF

      SELECT CASE (Element % TYPE % ElementCode / 100)
      CASE (5)
         SELECT CASE( Element % PDefs % TetraType )
         CASE (1)
            faceEdgeMap => TetraFaceEdgeMap1
         CASE (2)
            faceEdgeMap => TetraFaceEdgeMap2
         CASE DEFAULT
            CALL Fatal('PElementMaps::GetElementFaceEdgeMap','Unknown tetra type for p element')
         END SELECT
      CASE (6)
         faceEdgeMap => PyramidFaceEdgeMap
      CASE (7)
         faceEdgeMap => WedgeFaceEdgeMap
      CASE (8)
         faceEdgeMap => BrickFaceEdgeMap
      CASE DEFAULT
         CALL Fatal('PElementMaps::GetElementFaceEdgeMap','Unsupported element type')
      END SELECT
    END SUBROUTINE GetElementFaceEdgeMap


!------------------------------------------------------------------------------
!>   This subroutine initializes element mappings.
!------------------------------------------------------------------------------
    SUBROUTINE InitializeMappings() 
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      CALL Info('PElementMaps::InitializeMappings','Initializing mappings for elements',Level=10)

      ! Quad edge mappings
      QuadEdgeMap(1,:) = (/ 1,2 /)
      QuadEdgeMap(2,:) = (/ 2,3 /)
      QuadEdgeMap(3,:) = (/ 4,3 /)
      QuadEdgeMap(4,:) = (/ 1,4 /)

      ! Triangle edge mappings
      TriangleEdgeMap(1,:) = (/ 1,2 /)
      TriangleEdgeMap(2,:) = (/ 2,3 /)
      TriangleEdgeMap(3,:) = (/ 3,1 /)

      ! Brick edge mappings
      BrickEdgeMap(1,:) = (/ 1,2 /)
      BrickEdgeMap(2,:) = (/ 2,3 /)
      BrickEdgeMap(3,:) = (/ 4,3 /)
      BrickEdgeMap(4,:) = (/ 1,4 /)
      BrickEdgeMap(5,:) = (/ 5,6 /)
      BrickEdgeMap(6,:) = (/ 6,7 /)
      BrickEdgeMap(7,:) = (/ 8,7 /)
      BrickEdgeMap(8,:) = (/ 5,8 /)
      BrickEdgeMap(9,:) = (/ 1,5 /)
      BrickEdgeMap(10,:) = (/ 2,6 /)
      BrickEdgeMap(11,:) = (/ 3,7 /)
      BrickEdgeMap(12,:) = (/ 4,8 /)

      ! Brick face mappings
      BrickFaceMap(1,:) = (/ 1,2,3,4 /) ! xi,eta
      BrickFaceMap(2,:) = (/ 5,6,7,8 /) ! xi,eta
      BrickFaceMap(3,:) = (/ 1,2,6,5 /) ! xi,zeta
      BrickFaceMap(4,:) = (/ 2,3,7,6 /) ! eta,zeta
      ! BrickFaceMap(5,:) = (/ 3,4,8,7 /)
      BrickFaceMap(5,:) = (/ 4,3,7,8 /)
      ! BrickFaceMap(6,:) = (/ 4,1,5,8 /) 
      BrickFaceMap(6,:) = (/ 1,4,8,5 /)

      BrickFaceEdgeMap(1,:) = (/ 1,2,3,4 /)
      BrickFaceEdgeMap(2,:) = (/ 5,6,7,8 /)    
      BrickFaceEdgeMap(3,:) = (/ 1,10,5,9 /)
      BrickFaceEdgeMap(4,:) = (/ 2,11,6,10 /)
      ! BrickFaceEdgeMap(5,:) = (/ 3,12,7,11 /)
      BrickFaceEdgeMap(5,:) = (/ 3,11,7,12 /)
      ! BrickFaceEdgeMap(6,:) = (/ 4,9,8,12 /)
      BrickFaceEdgeMap(6,:) = (/ 4,12,8,9 /)

      ! Tetra edge mappings (not needed for enforcing parity!)
      ! Type 1
      TetraEdgeMap1(1,:) = (/ 1,2 /)
      TetraEdgeMap1(2,:) = (/ 2,3 /)
      TetraEdgeMap1(3,:) = (/ 1,3 /)
      TetraEdgeMap1(4,:) = (/ 1,4 /)
      TetraEdgeMap1(5,:) = (/ 2,4 /)
      TetraEdgeMap1(6,:) = (/ 3,4 /)
      ! Type 2
      TetraEdgeMap2(1,:) = (/ 1,2 /)
      TetraEdgeMap2(2,:) = (/ 3,2 /)
      TetraEdgeMap2(3,:) = (/ 1,3 /)
      TetraEdgeMap2(4,:) = (/ 1,4 /)
      TetraEdgeMap2(5,:) = (/ 2,4 /)
      TetraEdgeMap2(6,:) = (/ 3,4 /)

      ! Tetra face mappings (not needed for enforcing parity!)
      ! Type 1
      TetraFaceMap1(1,:) = (/ 1,2,3 /)
      TetraFaceMap1(2,:) = (/ 1,2,4 /)
      TetraFaceMap1(3,:) = (/ 2,3,4 /)
      TetraFaceMap1(4,:) = (/ 1,3,4 /)
      ! Type 2 
      TetraFaceMap2(1,:) = (/ 1,3,2 /)
      TetraFaceMap2(2,:) = (/ 1,2,4 /)
      TetraFaceMap2(3,:) = (/ 3,2,4 /)
      TetraFaceMap2(4,:) = (/ 1,3,4 /)

      ! Type 1 
      TetraFaceEdgeMap1(1,:) = (/ 1,2,3 /)
      TetraFaceEdgeMap1(2,:) = (/ 1,5,4 /)
      TetraFaceEdgeMap1(3,:) = (/ 2,6,5 /)
      TetraFaceEdgeMap1(4,:) = (/ 3,6,4 /)
      ! Type 2 
      TetraFaceEdgeMap2(1,:) = (/ 3,2,1 /)
      TetraFaceEdgeMap2(2,:) = (/ 1,5,4 /)
      TetraFaceEdgeMap2(3,:) = (/ 2,5,6 /)
      TetraFaceEdgeMap2(4,:) = (/ 3,6,4 /)

      ! Wedge edge mappings
      WedgeEdgeMap(1,:) = (/ 1,2 /)
      WedgeEdgeMap(2,:) = (/ 2,3 /)
      WedgeEdgeMap(3,:) = (/ 3,1 /)
      WedgeEdgeMap(4,:) = (/ 4,5 /)
      WedgeEdgeMap(5,:) = (/ 5,6 /)
      WedgeEdgeMap(6,:) = (/ 6,4 /)
      WedgeEdgeMap(7,:) = (/ 1,4 /)
      WedgeEdgeMap(8,:) = (/ 2,5 /)
      WedgeEdgeMap(9,:) = (/ 3,6 /)

      ! Wedge face mappings
      WedgeFaceMap(1,:) = (/ 1,2,3,0 /)
      WedgeFaceMap(2,:) = (/ 4,5,6,0 /)
      WedgeFaceMap(3,:) = (/ 1,2,5,4 /)
      WedgeFaceMap(4,:) = (/ 2,3,6,5 /)
      WedgeFaceMap(5,:) = (/ 3,1,4,6 /)

      WedgeFaceEdgeMap(1,:) = (/ 1,2,3,0 /)
      WedgeFaceEdgeMap(2,:) = (/ 4,5,6,0 /)
      WedgeFaceEdgeMap(3,:) = (/ 1,8,4,7 /)
      WedgeFaceEdgeMap(4,:) = (/ 2,9,5,8 /)
      WedgeFaceEdgeMap(5,:) = (/ 3,7,6,9 /)
      
      ! Pyramid edge mappings 
      PyramidEdgeMap(1,:) = (/ 1,2 /)
      PyramidEdgeMap(2,:) = (/ 2,3 /)
      PyramidEdgeMap(3,:) = (/ 4,3 /)
      PyramidEdgeMap(4,:) = (/ 1,4 /)
      PyramidEdgeMap(5,:) = (/ 1,5 /)
      PyramidEdgeMap(6,:) = (/ 2,5 /)
      PyramidEdgeMap(7,:) = (/ 3,5 /)
      PyramidEdgeMap(8,:) = (/ 4,5 /)

      ! Pyramid face mappings
      PyramidFaceMap(1,:) = (/ 1,2,3,4 /)
      PyramidFaceMap(2,:) = (/ 1,2,5,0 /)
      PyramidFaceMap(3,:) = (/ 2,3,5,0 /)
      PyramidFaceMap(4,:) = (/ 3,4,5,0 /)
      PyramidFaceMap(5,:) = (/ 4,1,5,0 /)

      PyramidFaceEdgeMap(1,:) = (/ 1,2,3,4 /)
      PyramidFaceEdgeMap(2,:) = (/ 1,6,5,0 /)
      PyramidFaceEdgeMap(3,:) = (/ 2,7,6,0 /)
      PyramidFaceEdgeMap(4,:) = (/ 3,8,7,0 /)
      PyramidFaceEdgeMap(5,:) = (/ 4,5,8,0 /)

      MInit = .TRUE.
    END SUBROUTINE InitializeMappings

!------------------------------------------------------------------------------
  FUNCTION getEdgeDOFs( Element, p ) RESULT(EdgeDOFs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: EdgeDOFs
    INTEGER, INTENT(IN) :: p

    IF (.NOT. ASSOCIATED(Element % PDefs) ) THEN
       EdgeDOFs = 0
       RETURN
    END IF

    EdgeDOFs = MAX(0, p-1)
!------------------------------------------------------------------------------
  END FUNCTION getEdgeDOFs
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>     Based on element face polynomial degree p, return degrees of freedom for
!>     given face. 
!------------------------------------------------------------------------------
  FUNCTION getFaceDOFs( Element, p, faceNumber ) RESULT(faceDOFs)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t), POINTER :: Element
!      INPUT: Element to get face dofs to 
!
!    INTEGER :: p
!      INPUT: Face polynomial degree p
!
!    INTEGER :: faceNumber
!      INPUT: Local number of face for element (important for wedges and 
!        pyramids).
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: faceDOFs
!       number of face dofs for Element
!    
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Element_t) :: Element
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN), OPTIONAL :: faceNumber
    INTEGER :: faceDOFs

    ! This function is not defined for non p elements
    IF (.NOT. ASSOCIATED(Element % PDefs) ) THEN
       faceDOFs = 0
       RETURN
    END IF

    faceDOFs = 0
    SELECT CASE(Element % TYPE % ElementCode / 100)
    ! Tetrahedron
    CASE (5)
       IF (p >= 3) faceDOFs = (p-1)*(p-2)/2
    ! Pyramid
    CASE (6)
       SELECT CASE(faceNumber)
          CASE (1)
             IF (p >= 4) faceDOFs = (p-2)*(p-3)/2
          CASE (2:5)
             IF (p >= 3) faceDOFs = (p-1)*(p-2)/2
       END SELECT
    ! Wedge
    CASE (7)
       SELECT CASE(faceNumber)
       CASE (1,2)
          IF (p >= 3) faceDOFs = (p-1)*(p-2)/2
       CASE (3:5)
          IF (p >= 4) faceDOFs = (p-2)*(p-3)/2
       END SELECT
    ! Brick   
    CASE (8)
       IF (p >= 4) faceDOFs = (p-2)*(p-3)/2
    CASE DEFAULT
       CALL Warn('MeshUtils::getFaceDOFs','Unsupported p element type')
       faceDOFs = p
    END SELECT

    faceDOFs = MAX(0, faceDOFs)
  END FUNCTION getFaceDOFs


!------------------------------------------------------------------------------
!>     Based on element bubble polynomial degree p, return degrees of freedom for
!>     given elements bubbles. 
!------------------------------------------------------------------------------
  FUNCTION getBubbleDOFs( Element, p) RESULT(bubbleDOFs)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t), POINTER :: Element
!      INPUT: Element to get bubble dofs to 
!
!    INTEGER :: p
!      INPUT: Element polynomial degree p
!
!  FUNCTION VALUE:
!    REAL(KIND=dp) :: bubbleDOFs
!       number of bubble dofs for Element
!    
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Element_t) :: Element
    INTEGER, INTENT(IN) :: p
    INTEGER :: bubbleDOFs
    
    ! This function is not defined for non p elements
    IF (.NOT. ASSOCIATED(Element % PDefs) ) THEN
       bubbleDOFs = 0
       RETURN
    END IF

    ! Select by element type
    bubbleDOFs = 0
    SELECT CASE (Element % TYPE % ElementCode / 100)
    ! Line 
    CASE (2)
      IF (p >= 2) bubbleDOFs = p - 1
    ! Triangle
    CASE (3)
      IF (p >= 3) bubbleDOFs = (p-1)*(p-2)/2
    ! Quad
    CASE (4)
       IF (p >= 4) bubbleDOFs = (p-2)*(p-3)/2
    ! Tetrahedron
    CASE (5)
       IF (p >= 4) bubbleDOFs = (p-1)*(p-2)*(p-3)/6
    ! Pyramid
    CASE (6)
       IF (p >= 4) bubbleDOFs = (p-1)*(p-2)*(p-3)/6
    ! Wedge
    CASE (7)
       IF (p >= 5) bubbleDOFs = (p-2)*(p-3)*(p-4)/6
    ! Brick
    CASE (8)
       IF (p >= 6) bubbleDOFs = (p-3)*(p-4)*(p-5)/6
    CASE DEFAULT
       CALL Warn('MeshUtils::getBubbleDOFs','Unsupported p element type')
       bubbleDOFs = p
    END SELECT

    bubbleDOFs = MAX(0, bubbleDOFs)
  END FUNCTION getBubbleDOFs



!------------------------------------------------------------------------------
  FUNCTION isActivePElement(Element) RESULT(retVal)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t), INTENT(IN) :: Element
    INTEGER :: c
    LOGICAL :: retVal

    retVal = isPelement(Element)

    IF(ASSOCIATED(CurrentModel % Solver)) THEN
      IF(ALLOCATED(CurrentModel % Solver % Def_Dofs)) THEN
        c = Element % Type % ElementCode / 100
        retVal = retVal.AND.ANY(CurrentModel % Solver % Def_Dofs(c,:,6)>0)
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION isActivePElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Checks if given element is a p element.   
!------------------------------------------------------------------------------
    FUNCTION isPElement( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p triangle, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: Element
      LOGICAL :: retVal
      retVal = ASSOCIATED(Element % PDefs)
!------------------------------------------------------------------------------
    END FUNCTION isPElement
!------------------------------------------------------------------------------
    

!------------------------------------------------------------------------------
!>    Function checks if given element is p element triangle
!------------------------------------------------------------------------------
    FUNCTION isPTriangle( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p triangle, .FALSE. otherwise
!
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: Element
      LOGICAL :: retVal

      ! Check elementcode and p element flag
      retVal = Element % TYPE % ElementCode/100==3 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPTriangle
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Function checks if given element is p element quad
!------------------------------------------------------------------------------
    FUNCTION isPQuad( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p quad, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: Element
      LOGICAL :: retVal

      ! Check elementcode and p element flag
      retVal = Element % TYPE % ElementCode/100==4 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPQuad
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Function checks if given element is p element tetra
!------------------------------------------------------------------------------
    FUNCTION isPTetra( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p tetra, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: Element
      LOGICAL :: retVal

      retVal = Element % TYPE % ElementCode/100==5 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPTetra
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Function checks if given element is p element wedge
!------------------------------------------------------------------------------
    FUNCTION isPWedge( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p wedge, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: Element
      LOGICAL :: retVal

      retVal = Element % TYPE % ElementCode/100==7 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPWedge
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Function checks if given element is p element pyramid
!------------------------------------------------------------------------------
    FUNCTION isPPyramid( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p pyramid, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: ELement
      LOGICAL :: retVal

      retVal = Element % TYPE % ElementCode/100==6 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPPyramid
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>    Function checks if given element is p element brick
!------------------------------------------------------------------------------
    FUNCTION isPBrick( Element ) RESULT(retVal)
!------------------------------------------------------------------------------
!
!  ARGUMENTS:
!    Type(Element_t) :: Element
!      INPUT: Element to check
!
!  FUNCTION VALUE:
!    LOGICAL :: retVal
!       .TRUE. if given element is a p brick, .FALSE. otherwise
!    
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), INTENT(IN) :: ELement
      LOGICAL :: retVal

      retVal = Element % TYPE % ElementCode/100==8 .AND. isPElement(Element)
!------------------------------------------------------------------------------
    END FUNCTION isPBrick
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Get the number of gauss points for P-elements.
!------------------------------------------------------------------------------
    FUNCTION getNumberOfGaussPoints( Element, Mesh ) RESULT(ngp)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Mesh_t) :: Mesh 
      TYPE(Element_t) :: Element
      INTEGER :: ngp
!------------------------------------------------------------------------------
      INTEGER :: edgeP, faceP, bubbleP, nb, maxp

      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementBase::getNumberOfGaussPoints','Element not p element')
         ngp = 0
         RETURN
      END IF

      ! Max p of edges
      edgeP = 0
      IF ( Element % TYPE % DIMENSION == 2 .OR. &
           Element % TYPE % DIMENSION == 3) THEN
         edgeP = getEdgeP( Element, Mesh )
      END IF
      
      ! Max p of faces
      faceP = 0
      IF ( Element % TYPE % DIMENSION == 3 ) THEN
         faceP = getFaceP(Element, Mesh )
      END IF
      
      ! Element bubble p
      bubbleP = 0
      IF (Element % BDOFs > 0) THEN
         bubbleP = Element % PDefs % P
         
         SELECT CASE( Element % TYPE % ElementCode / 100 )
         CASE(3)
             nb = MAX( GetBubbleDOFs( Element, bubbleP ), Element % BDOFs )
             bubbleP = CEILING( ( 3.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )

         CASE(4)
             nb = MAX( GetBubbleDOFs( Element, bubbleP ), Element % BDOFs )
             bubbleP = CEILING( ( 5.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )

         CASE(5)
             nb = MAX( GetBubbleDOFs(Element, bubbleP ), Element % BDOFs )
             bubbleP=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0 / &
                    (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

         CASE(6)
             nb = MAX( GetBubbleDOFs(Element, bubbleP ), Element % BDOFs )
             bubbleP=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0 / &
                    (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

         CASE(7)
             nb = MAX( GetBubbleDOFs( Element, bubbleP ), Element % BDOFs )
             bubbleP=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0 / &
                    (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+3)

         CASE(8)
             nb = MAX( GetBubbleDOFs(Element, bubbleP ), Element % BDOFs )
             bubbleP=CEILING(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0 / &
                    (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+4)
         END SELECT
      END IF

      ! Get number of Gauss points. Number needed is 2*max(p)=2*(max(p)+1)/2
      maxp = MAX(1, edgeP, faceP, bubbleP) + 1
      ngp = maxp ** Element % TYPE % DIMENSION
!------------------------------------------------------------------------------
    END FUNCTION getNumberOfGaussPoints
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
    FUNCTION getEdgeP( Element, Mesh ) RESULT(edgeP)
!------------------------------------------------------------------------------
      IMPLICIT NONE
      
      TYPE(Mesh_t) :: Mesh
      TYPE(Element_t) :: Element 
      TYPE(Element_t), POINTER :: Edge
      
      INTEGER :: edgeP, i

      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementBase::getEdgeP','Element not p element')
         edgeP = 0
         RETURN
      END IF

      ! Get max p of edges of element if any
      edgeP = 0
      IF (ASSOCIATED(Element % EdgeIndexes)) THEN
         DO i=1, Element % TYPE % NumberOfEdges
            Edge => Mesh % Edges(Element % EdgeIndexes(i))
            ! Here if edge has no dofs it effectively has degree of 0
            IF (Edge % BDOFs <= 0) CYCLE
            edgeP = MAX(edgeP,Edge % PDefs % P)
         END DO
      END IF
!------------------------------------------------------------------------------
    END FUNCTION getEdgeP
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    FUNCTION getFaceP( Element, Mesh ) RESULT(faceP)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t) :: Element
      TYPE(Element_t), POINTER :: Face
      INTEGER :: faceP, i
      TYPE(Mesh_t) :: Mesh
      
      IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
         CALL Warn('PElementBase::getFaceP','Element not p element')
         faceP = 0
         RETURN
      END IF

      faceP = 0
      IF (ASSOCIATED(Element % FaceIndexes)) THEN
         DO i=1,Element % TYPE % NumberOfFaces
            Face => Mesh % Faces( Element % FaceIndexes(i) )
            ! Here if face has no dofs it effectively has degree of 0
            IF (Face % BDOFs <= 0) CYCLE
            faceP = MAX(faceP, Face % PDefs % P)
         END DO
      END IF
!------------------------------------------------------------------------------
    END FUNCTION getFaceP
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    FUNCTION getNumberOfGaussPointsFace( Face, Mesh ) RESULT(ngp)
!------------------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(Element_t), POINTER :: Face, Edge
      INTEGER :: ngp
      TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
      INTEGER :: edgeP, i, maxp

      ! Get max p of edges contained in face
      edgeP = 0
      DO i=1, Face % TYPE % NumberOfEdges
         Edge => Mesh % Edges ( Face % EdgeIndexes(i) )
         edgeP = MAX(edgeP, Edge % PDefs % P)
      END DO

      ! If no face dofs use max edge dofs as number of points
      IF (Face % BDOFs <= 0) THEN
         ngp = (edgeP + 1) ** 2 
         RETURN
      END IF

      ! Get number of Gauss points
      maxp = MAX(edgeP, Face % PDefs % P) + 1
      ngp = maxp ** 2
!------------------------------------------------------------------------------
    END FUNCTION
!------------------------------------------------------------------------------
 
  
!------------------------------------------------------------------------------
!>     Subroutine for getting reference p element nodes (because these are NOT
!>     yet defined in element description files)
!------------------------------------------------------------------------------
  SUBROUTINE GetRefPElementNodes(Element, U, V, W, PerformCheck)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    REAL(KIND=dp), POINTER CONTIG :: U(:), V(:), W(:)
    LOGICAL, OPTIONAL :: PerformCheck
!--------------------------------------------------------------------------------
    INTEGER :: n
!--------------------------------------------------------------------------------
    LOGICAL :: PElementRequired
!--------------------------------------------------------------------------------    
    PElementRequired = .TRUE.
    IF ( PRESENT(PerformCheck) ) PElementRequired = PerformCheck
    IF (PElementRequired) THEN
       ! If element is not p element return
       IF ( .NOT. isPElement(Element) ) THEN 
          CALL Warn('PElementMaps::GetRefPElementNodes','Element given not a p element')
          RETURN
       END IF
    END IF

    ! Reserve space for element nodes
    n = Element % TYPE % NumberOfNodes

    ! Select by element type given
    SELECT CASE(Element % TYPE % ElementCode / 100)
    ! Line
    CASE(2)
       U(1:n) = (/ -1d0,1d0 /)
    ! Triangle
    CASE(3)
       U(1:n) = (/ -1d0,1d0,0d0 /)
       V(1:n) = (/ 0d0,0d0,SQRT(3.0d0) /)
    ! Quad
    CASE(4)
       U(1:n) = (/ -1d0,1d0,1d0,-1d0 /)
       V(1:n) = (/ -1d0,-1d0,1d0,1d0 /)
    ! Tetrahedron
    CASE(5)
       U(1:n) = (/ -1d0,1d0,0d0,0d0 /)
       V(1:n) = (/ 0d0,0d0,SQRT(3.0d0),1.0d0/SQRT(3.0d0) /)
       W(1:n) = (/ 0d0,0d0,0d0,2*SQRT(2.0d0/3.0d0) /)
    ! Pyramid
    CASE(6)
       U(1:n) = (/ -1d0,1d0,1d0,-1d0,0d0 /)
       V(1:n) = (/ -1d0,-1d0,1d0,1d0,0d0 /)
       W(1:n) = (/ 0d0,0d0,0d0,0d0,SQRT(2.0d0) /)
    ! Wedge
    CASE(7)
       U(1:n) = (/ -1d0,1d0,0d0,-1d0,1d0,0d0 /)
       V(1:n) = (/ 0d0,0d0,SQRT(3.0d0),0d0,0d0,SQRT(3.0d0) /)
       W(1:n) = (/ -1d0,-1d0,-1d0,1d0,1d0,1d0 /)
    ! Brick
    CASE(8)
       U(1:n) = (/ -1d0,1d0,1d0,-1d0,-1d0,1d0,1d0,-1d0 /)
       V(1:n) = (/ -1d0,-1d0,1d0,1d0,-1d0,-1d0,1d0,1d0 /)
       W(1:n) = (/ -1d0,-1d0,-1d0,-1d0,1d0,1d0,1d0,1d0 /)
    CASE DEFAULT
       CALL Warn('PElementMaps::GetRefPElementNodes','Unknown element type')
    END SELECT
!------------------------------------------------------------------------------
    END SUBROUTINE GetRefPElementNodes
!------------------------------------------------------------------------------


END MODULE

!> \}
