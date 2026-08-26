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
!>  Element geometry: normal/surface vectors, coordinate mapping, inside tests,
!>  cutting/splitting.
!--------------------------------------------------------------------------------
MODULE ElementGeometry
   USE ElementBasis
   USE Messages
   USE Integration
   USE LinearAlgebra
   USE CoordinateSystems
   USE PElementMaps
   USE PElementBase
   USE Lists

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TriangleInside, QuadInside, TetraInside, BrickInside, &
             CheckPassiveElement, CheckNormalDirection, CheckNormalDirectionParent, &
             NormalVector, SurfaceVector, &
             LineFaceIntersection, LineFaceIntersection2, PointFaceDistance, &
             GlobalToLocal, &
             getTriangleFaceDirection, getSquareFaceDirection, wedgeOrdering, &
             ComputeRotationMatrix, CutSingleElement, SplitSingleElement

CONTAINS

!------------------------------------------------------------------------------


  

!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a triangle, whose node
!>     coordinates are given in nx,ny,nz. Method: Invert the basis
!>     functions....
!------------------------------------------------------------------------------
  FUNCTION TriangleInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:) !< Node coordinate arrays
    REAL(KIND=dp) :: x,y,z             !< point which to consider
    LOGICAL :: inside                  !< result of the in/out test
!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: a00,a01,a10,a11,b00,b01,b10,b11,detA,px,py,u,v
!------------------------------------------------------------------------------

    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y ) RETURN

    A00 = nx(2) - nx(1)
    A01 = nx(3) - nx(1)
    A10 = ny(2) - ny(1)
    A11 = ny(3) - ny(1)

    detA = A00*A11 - A01*A10
    IF ( ABS(detA) < AEPS ) RETURN

    detA = 1 / detA

    B00 =  A11*detA
    B01 = -A01*detA
    B10 = -A10*detA
    B11 =  A00*detA

    px = x - nx(1)
    py = y - ny(1)
    u = 0.0d0
    v = 0.0d0

    u = B00*px + B01*py
    IF ( u < 0.0d0 .OR. u > 1.0d0 ) RETURN

    v = B10*px + B11*py
    IF ( v < 0.0d0 .OR. v > 1.0d0 ) RETURN

    inside = (u + v <=  1.0d0)
!------------------------------------------------------------------------------
   END FUNCTION TriangleInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a quadrilateral, whose
!>     node coordinates are given in nx,ny,nz. Method: Invert the
!>     basis functions....
!------------------------------------------------------------------------------
   FUNCTION QuadInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:) !< Node coordinate arrays
    REAL(KIND=dp) :: x,y,z             !< point which to consider
    LOGICAL :: inside                  !< result of the in/out test
!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: r,a,b,c,d,ax,bx,cx,dx,ay,by,cy,dy,px,py,u,v
!------------------------------------------------------------------------------
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y ) RETURN

    ax = 0.25*(  nx(1) + nx(2) + nx(3) + nx(4) )
    bx = 0.25*( -nx(1) + nx(2) + nx(3) - nx(4) )
    cx = 0.25*( -nx(1) - nx(2) + nx(3) + nx(4) )
    dx = 0.25*(  nx(1) - nx(2) + nx(3) - nx(4) )

    ay = 0.25*(  ny(1) + ny(2) + ny(3) + ny(4) )
    by = 0.25*( -ny(1) + ny(2) + ny(3) - ny(4) )
    cy = 0.25*( -ny(1) - ny(2) + ny(3) + ny(4) )
    dy = 0.25*(  ny(1) - ny(2) + ny(3) - ny(4) )

    px = x - ax
    py = y - ay

    a = cy*dx - cx*dy
    b = bx*cy - by*cx + dy*px - dx*py
    c = by*px - bx*py

    u = 0.0d0
    v = 0.0d0

    IF ( ABS(a) < AEPS ) THEN
      r = -c / b
      IF ( r < -1.0d0 .OR. r > 1.0d0 ) RETURN

      v = r
      u = (px - cx*r)/(bx + dx*r)
      inside = (u >= -1.0d0 .AND. u <= 1.0d0)
      RETURN
    END IF

    d = b*b - 4*a*c
    IF ( d < 0.0d0 ) RETURN

    d = SQRT(d)
    IF ( b>0 ) THEN
      r = -2*c/(b+d)
    ELSE
      r = (-b+d)/(2*a)
    END IF
    IF ( r >= -1.0d0 .AND. r <= 1.0d0 ) THEN
      v = r
      u = (px - cx*r)/(bx + dx*r)
        
      IF ( u >= -1.0d0 .AND. u <= 1.0d0 ) THEN
        inside = .TRUE.
        RETURN
      END IF
    END IF

    IF ( b>0 ) THEN
      r = -(b+d)/(2*a)
    ELSE
      r = 2*c/(-b+d)
    END IF
    IF ( r >= -1.0d0 .AND. r <= 1.0d0 ) THEN
      v = r
      u = (px - cx*r)/(bx + dx*r)
      inside = u >= -1.0d0 .AND. u <= 1.0d0
      RETURN
    END IF
!------------------------------------------------------------------------------
  END FUNCTION QuadInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a tetrahedron, whose
!>     node coordinates are given in nx,ny,nz. Method: Invert the
!>     basis functions....
!------------------------------------------------------------------------------
  FUNCTION TetraInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:) !< Node coordinate arrays
    REAL(KIND=dp) :: x,y,z             !< point which to consider
    LOGICAL :: inside                  !< result of the in/out test
!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: A00,A01,A02,A10,A11,A12,A20,A21,A22,detA
    REAL(KIND=dp) :: B00,B01,B02,B10,B11,B12,B20,B21,B22
    REAL(KIND=dp) :: px,py,pz,u,v,w
!------------------------------------------------------------------------------
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y .OR. MAXVAL(nz) < z ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y .OR. MINVAL(nz) > z ) RETURN

    A00 = nx(2) - nx(1)
    A01 = nx(3) - nx(1)
    A02 = nx(4) - nx(1)

    A10 = ny(2) - ny(1)
    A11 = ny(3) - ny(1)
    A12 = ny(4) - ny(1)

    A20 = nz(2) - nz(1)
    A21 = nz(3) - nz(1)
    A22 = nz(4) - nz(1)

    detA =        A00*(A11*A22 - A12*A21)
    detA = detA + A01*(A12*A20 - A10*A22)
    detA = detA + A02*(A10*A21 - A11*A20)
    IF ( ABS(detA) < AEPS ) RETURN

    detA = 1 / detA

    px = x - nx(1)
    py = y - ny(1)
    pz = z - nz(1)

    B00 = (A11*A22 - A12*A21)*detA
    B01 = (A21*A02 - A01*A22)*detA
    B02 = (A01*A12 - A11*A02)*detA

    u = B00*px + B01*py + B02*pz
    IF ( u < 0.0d0 .OR. u > 1.0d0 ) RETURN


    B10 = (A12*A20 - A10*A22)*detA
    B11 = (A00*A22 - A20*A02)*detA
    B12 = (A10*A02 - A00*A12)*detA

    v = B10*px + B11*py + B12*pz
    IF ( v < 0.0d0 .OR. v > 1.0d0 ) RETURN


    B20 = (A10*A21 - A11*A20)*detA
    B21 = (A01*A20 - A00*A21)*detA
    B22 = (A00*A11 - A10*A01)*detA

    w = B20*px + B21*py + B22*pz
    IF ( w < 0.0d0 .OR. w > 1.0d0 ) RETURN

    inside = (u + v + w) <= 1.0d0
!------------------------------------------------------------------------------
  END FUNCTION TetraInside
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!>     Figure out if given point x,y,z is inside a brick, whose node coordinates
!>     are given in nx,ny,nz. Method: Divide to tetrahedrons.
!------------------------------------------------------------------------------
  FUNCTION BrickInside( nx,ny,nz,x,y,z ) RESULT(inside)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: nx(:),ny(:),nz(:) !< Node coordinate arrays
    REAL(KIND=dp) :: x,y,z             !< point which to consider
    LOGICAL :: inside                  !< result of the in/out test
!------------------------------------------------------------------------------
!   Local variables
!------------------------------------------------------------------------------
    INTEGER :: i,j
    REAL(KIND=dp) :: px(4),py(4),pz(4),r,s,t,maxx,minx,maxy,miny,maxz,minz
    INTEGER :: map(3,12)
!------------------------------------------------------------------------------
    map = RESHAPE( [ 0,1,2,   0,2,3,   4,5,6,   4,6,7,   3,2,6,   3,6,7,  &
     1,5,6,   1,6,2,   0,4,7,   0,7,3,   0,1,5,   0,5,4 ], [ 3,12 ] ) + 1
    
    inside = .FALSE.

    IF ( MAXVAL(nx) < x .OR. MAXVAL(ny) < y .OR. MAXVAL(nz) < z ) RETURN
    IF ( MINVAL(nx) > x .OR. MINVAL(ny) > y .OR. MINVAL(nz) > z ) RETURN

    px(1) = 0.125d0 * SUM(nx)
    py(1) = 0.125d0 * SUM(ny)
    pz(1) = 0.125d0 * SUM(nz)

    DO i=1,12
      px(2:4) = nx(map(1:3,i))
      py(2:4) = ny(map(1:3,i))
      pz(2:4) = nz(map(1:3,i))

      IF ( TetraInside( px,py,pz,x,y,z ) ) THEN
        inside = .TRUE.
        RETURN
      END IF
    END DO
!------------------------------------------------------------------------------
  END FUNCTION BrickInside
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Check if the current element has been defined passive.
!> This is done by inspecting a looking an the values of "varname Passive"
!> in the Body Force section. It is determined to be passive if it has 
!> more positive than negative hits in an element.
!------------------------------------------------------------------------------
  FUNCTION CheckPassiveElement( UElement )  RESULT( IsPassive )
    !------------------------------------------------------------------------------
    TYPE(Element_t), OPTIONAL, TARGET :: UElement
    LOGICAL :: IsPassive
    !------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element,tmp
    REAL(KIND=dp), ALLOCATABLE :: Passive(:)
    INTEGER :: body_id, bf_id, nlen, NbrNodes, PassNodes
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: PassName
    LOGICAL :: NoPassiveElements = .FALSE.    
    TYPE(Solver_t), POINTER :: pSolver, PrevSolver => NULL()
    TYPE(ValueList_t), POINTER :: BodyForce => NULL()
    INTEGER :: ActiveMin = -1, PassiveMin = -1, prev_body_id = -1
    LOGICAL :: DoCheck = .FALSE.
    
    SAVE Passive, NoPassiveElements, PrevSolver, PassName, prev_body_id, &
        BodyForce, ActiveMin, PassiveMin, DoCheck
    !$OMP THREADPRIVATE(Passive, NoPassiveElements, PrevSolver, PassName, prev_body_id, &
    !$OMP               BodyForce, ActiveMin, PassiveMin, DoCheck ) 
    !------------------------------------------------------------------------------
    IsPassive = .FALSE.
    pSolver => CurrentModel % Solver
    
    IF( .NOT. ASSOCIATED( pSolver, PrevSolver ) ) THEN
      PrevSolver => pSolver          
      nlen = CurrentModel % Solver % Variable % NameLen
      PassName = GetVarName(CurrentModel % Solver % Variable) // ' Passive'     
      NoPassiveElements = .NOT. ListCheckPresentAnyBodyForce(CurrentModel, PassName)

      ! Nullify the BodyForce memories also if we have new solver.
      prev_body_id = -1
    END IF
    
    IF( NoPassiveElements ) RETURN       

    IF (PRESENT(UElement)) THEN
      tmp => CurrentModel % CurrentElement
      Element => UElement
      CurrentModel % CurrentElement => Element
    ELSE
#ifdef _OPENMP
      IF (omp_in_parallel()) THEN
        CALL Fatal('CheckPassiveElement', &
             'Need an element to update inside a threaded region')
      END IF
#endif
      Element => CurrentModel % CurrentElement
    END IF

    body_id = Element % BodyId 
    IF ( body_id <= 0 )  RETURN   ! body_id == 0 for boundary elements

    ! Do some mundane list operations if we have different body than previously. 
    IF(body_id /= prev_body_id ) THEN
      prev_body_id = body_id
      
      bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values, &
          'Body Force', DoCheck , minv=1,maxv=CurrentModel % NumberOfBodyForces )
      IF(DoCheck) THEN
        BodyForce => CurrentModel % BodyForces(bf_id) % Values
        DoCheck = ListCheckPresent( BodyForce, PassName)
      END IF
      IF(DoCheck) THEN
        PassiveMin = ListGetInteger( pSolver % Values,'Passive Element Min Nodes',Found )      
        IF(.NOT. Found) PassiveMin = ListGetInteger( BodyForce,'Passive Element Min Nodes',Found )              
        ActiveMin = ListGetInteger( pSolver % Values,'Active Element Min Nodes',Found )               
        IF(.NOT. Found) ActiveMin = ListGetInteger( BodyForce,'Active Element Min Nodes',Found )
      END IF
    END IF
    
    IF(DoCheck) THEN 
      NbrNodes = Element % TYPE % NumberOfNodes
      IF ( ALLOCATED(Passive) ) THEN
        IF ( SIZE(Passive) < NbrNodes ) THEN
          DEALLOCATE(Passive)
          ALLOCATE( Passive(NbrNodes) )
        END IF
      ELSE
        ALLOCATE( Passive(NbrNodes) )
      END IF
      Passive(1:NbrNodes) = ListGetReal( BodyForce, PassName, NbrNodes, Element % NodeIndexes )
      PassNodes = COUNT(Passive(1:NbrNodes)>0)

      ! Go through the extremum cases first, and if the element is not either fully 
      ! active or passive, then check for some possible given criteria for determining 
      ! the element active / passive. 
      !------------------------------------------------------------------------------
      IF( PassNodes == 0 ) THEN
        CONTINUE
      ELSE IF( PassNodes == NbrNodes ) THEN
        IsPassive = .TRUE.
      ELSE
        IF( PassiveMin > 0 ) THEN
          IsPassive = ( PassNodes >= PassiveMin )
        ELSE IF( ActiveMin > 0 ) THEN
          IsPassive = ( PassNodes > NbrNodes - ActiveMin )
        ELSE
          IsPassive = ( 2*PassNodes > NbrNodes )
        END IF
      END IF
    END IF

    IF (PRESENT(UElement)) THEN
       CurrentModel % CurrentElement => tmp
    END IF
!------------------------------------------------------------------------------
  END FUNCTION CheckPassiveElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Normal will point into body with lower body ID.
!> or outwards, if no elements on the other side.
!------------------------------------------------------------------------------
  SUBROUTINE CheckNormalDirection( Boundary,Normal,x,y,z,turn )
!------------------------------------------------------------------------------

    TYPE(Element_t) :: Boundary
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Normal(3),x,y,z
    LOGICAL, OPTIONAL :: turn
!------------------------------------------------------------------------------

    TYPE (Element_t), POINTER :: Element,LeftElement,RightElement

    INTEGER :: LMat,RMat,n,k

    REAL(KIND=dp) :: u,v,w,dCoord(3)
    REAL(KIND=dp), ALLOCATABLE :: nx(:),ny(:),nz(:)
    REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES)
    LOGICAL :: LPassive
!------------------------------------------------------------------------------
    
    IF(.NOT. ASSOCIATED( Boundary % BoundaryInfo ) )  RETURN
    
    k = Boundary % BoundaryInfo % OutBody

    LeftElement => Boundary % BoundaryInfo % Left

    Element => Null()
    IF ( ASSOCIATED(LeftELement) ) THEN
       RightElement => Boundary % BoundaryInfo % Right
       IF ( ASSOCIATED( RightElement ) ) THEN ! we have a body-body boundary        
         IF ( k > 0 ) THEN ! declared outbody 
           IF ( LeftElement % BodyId == k ) THEN
             Element => RightElement
           ELSE
             Element => LeftElement
           END IF
         ELSE IF (LeftElement % BodyId > RightElement % BodyId) THEN ! normal pointing into body with lower body ID
             Element => LeftElement
         ELSE IF (LeftElement % BodyId < RightElement % BodyId) THEN! normal pointing into body with lower body ID
           Element => RightElement
         ELSE ! active/passive boundary
           LPassive = CheckPassiveElement( LeftElement )
           IF (LPassive .NEQV. CheckPassiveElement( RightElement )) THEN 
             IF(LPassive) THEN
               Element => RightElement
             ELSE
               Element => LeftElement
             END IF
           END IF
         END IF
       ELSE ! body-vacuum boundary from left->right
         Element => LeftElement
       END IF
    ELSE! body-vacuum boundary from right->left
       Element => Boundary % BoundaryInfo % Right
    END IF

    IF ( .NOT. ASSOCIATED(Element) ) RETURN

    n = Element % TYPE % NumberOfNodes

    ALLOCATE( nx(n), ny(n), nz(n) )

    nx(1:n) = CurrentModel % Nodes % x(Element % NodeIndexes)
    ny(1:n) = CurrentModel % Nodes % y(Element % NodeIndexes)
    nz(1:n) = CurrentModel % Nodes % z(Element % NodeIndexes)

    SELECT CASE( Element % TYPE % ElementCode / 100 )

    CASE(2,4,8)
      u = 0.0_dp
      v = 0.0_dp
      w = 0.0_dp
    CASE(3)
      u = 1.0d0/3
      v = 1.0d0/3
      w = 0.0d0
    CASE(5)
      u = 1.0d0/4
      v = 1.0d0/4
      w = 1.0d0/4
    CASE(6)
      u = 0.0
      v = 0.0
      w = 1.0d0/3
    CASE(7)
      u = 1.0d0/3
      v = 1.0d0/3
      w = 0.0d0
    CASE DEFAULT
      CALL Fatal('CheckNormalDirection','Invalid elementcode for parent element!')   
      
    END SELECT

    CALL NodalBasisFunctions( n, Basis, Element, u, v, w )
    dCoord(1) = DOT_PRODUCT( Basis(1:n), nx(1:n) ) - x
    dCoord(2) = DOT_PRODUCT( Basis(1:n), ny(1:n) ) - y
    dCoord(3) = DOT_PRODUCT( Basis(1:n), nz(1:n) ) - z
  
    IF ( PRESENT(turn) ) turn = .FALSE.
    IF ( SUM( dCoord * Normal ) > 0 ) THEN
       IF ( Element % BodyId /= k ) THEN
          Normal = -Normal
          IF ( PRESENT(turn) ) turn = .TRUE.
       END IF
    ELSE IF (  Element % BodyId == k ) THEN
       Normal = -Normal
       IF ( PRESENT(turn) ) turn = .TRUE.
    END IF
    DEALLOCATE( nx,ny,nz )
!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalDirection
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Normal will point out from the parent.
!------------------------------------------------------------------------------
  SUBROUTINE CheckNormalDirectionParent( Boundary,Normal,x,y,z,Element,turn )
!------------------------------------------------------------------------------

    TYPE(Element_t) :: Boundary
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Normal(3),x,y,z
    TYPE(Element_t), POINTER :: Element
    LOGICAL, OPTIONAL :: turn
!------------------------------------------------------------------------------
    INTEGER :: n,k
    REAL(KIND=dp) :: u,v,w,x1,y1,z1
    REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES)
    REAL(KIND=dp), ALLOCATABLE :: nx(:),ny(:),nz(:)
    LOGICAL :: LPassive
!------------------------------------------------------------------------------

    IF( PRESENT( turn ) ) turn = .FALSE.

    IF ( .NOT. ASSOCIATED(Element) ) RETURN

    n = Element % TYPE % NumberOfNodes

    ALLOCATE( nx(n), ny(n), nz(n) )

    nx(1:n) = CurrentModel % Nodes % x(Element % NodeIndexes)
    ny(1:n) = CurrentModel % Nodes % y(Element % NodeIndexes)
    nz(1:n) = CurrentModel % Nodes % z(Element % NodeIndexes)

    SELECT CASE( Element % TYPE % ElementCode / 100 )
    CASE(2,4,8); u = 0.0_dp;    v = 0.0_dp;    w = 0.0_dp
    CASE(3);     u = 1.0d0/3;   v = 1.0d0/3;   w = 0.0_dp
    CASE(5);     u = 1.0d0/4;   v = 1.0d0/4;   w = 1.0d0/4
    CASE(6);     u = 0.0_dp;    v = 0.0_dp;    w = 1.0d0/3
    CASE(7);     u = 1.0d0/3;   v = 1.0d0/3;   w = 0.0_dp
    CASE DEFAULT
      CALL Fatal('CheckNormalDirection','Invalid elementcode for parent element!')
    END SELECT

    CALL NodalBasisFunctions( n, Basis, Element, u, v, w )
    x1 = DOT_PRODUCT( Basis(1:n), nx(1:n) ) - x
    y1 = DOT_PRODUCT( Basis(1:n), ny(1:n) ) - y
    z1 = DOT_PRODUCT( Basis(1:n), nz(1:n) ) - z
    
    ! Swap the sign if the tentative normal points to the center, it should point outward
    IF ( x1*Normal(1) + y1*Normal(2) + z1*Normal(3) > 0 ) THEN
      Normal = -Normal
      IF ( PRESENT(turn) ) turn = .TRUE.
    END IF

    DEALLOCATE( nx,ny,nz )
!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalDirectionParent
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!> Gives the normal vector of a boundary element.
!> For noncurved elements the normal vector does not depend on the local coordinate
!> while otherwise it does. There are different uses of the function where some
!> do not have the luxury of knowing the local coordinates and hence the center
!> point is used as default.
!------------------------------------------------------------------------------
  RECURSIVE FUNCTION NormalVector( Boundary,BoundaryNodes,u0,v0,Check,Parent,Turn) RESULT(Normal)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Boundary
    TYPE(Nodes_t)   :: BoundaryNodes
    REAL(KIND=dp), OPTIONAL :: u0,v0
    LOGICAL, OPTIONAL :: Check
    TYPE(Element_t), POINTER, OPTIONAL :: Parent
    LOGICAL, OPTIONAL :: Turn
    REAL(KIND=dp) :: Normal(3)
!------------------------------------------------------------------------------
    LOGICAL :: CheckBody, CheckParent
    TYPE(ElementType_t),POINTER :: elt
    REAL(KIND=dp) :: u,v,Auu,Auv,Avu,Avv,detA,x,y,z
    REAL(KIND=dp) :: dxdu,dxdv,dydu,dydv,dzdu,dzdv
    REAL(KIND=dp), DIMENSION(:), POINTER :: nx,ny,nz
    REAL(KIND=dp) :: Tangent1(3), Tangent2(3)
    REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES), dLBasisdx(MAX_ELEMENT_NODES, 3)
    TYPE(Nodes_t) :: ParentNodes
    TYPE(Element_t), POINTER :: pParent
    INTEGER :: n, meshDim, elemDim
    
!------------------------------------------------------------------------------

    nx => BoundaryNodes % x
    ny => BoundaryNodes % y
    nz => BoundaryNodes % z   

    elemDim = Boundary % TYPE % DIMENSION

    IF(ASSOCIATED( CurrentModel % Mesh ) ) THEN
      meshDim = CurrentModel % Mesh % MeshDim
    ELSE
      meshDim = CurrentModel % dimension
    END IF
      
    SELECT CASE ( elemDim )

    CASE ( 0 ) 
      Normal(1) = 1.0_dp
      Normal(2:3) = 0.0_dp

    CASE ( 1 )
      IF( meshDim == 3 ) THEN
        ! We have 1D element but 3D mesh
        ! Define the normal in the plane defined by the 2D parent element.
        IF( PRESENT( u0 ) ) THEN
          u = u0
        ELSE
          u = 0.0_dp
        END IF

        ! 1st tangent vector is defined by the edge direction
        n = Boundary % TYPE % NumberOfNodes
        CALL NodalFirstDerivatives( n, dLBasisdx, Boundary, u, 0.0_dp, 0.0_dp )
        dxdu = DOT_PRODUCT( dLBasisdx(1:n,1), nx(1:n) )
        dydu = DOT_PRODUCT( dLBasisdx(1:n,1), ny(1:n) )
        dzdu = DOT_PRODUCT( dLBasisdx(1:n,1), nz(1:n) )
        
        detA = dxdu*dxdu + dydu*dydu + dzdu*dzdu
        IF ( detA <= 0._dp ) THEN
          Normal = 0._dp
          RETURN
        END IF
        detA = 1.0_dp / SQRT(detA)
        Tangent1(1) = dxdu * detA
        Tangent1(2) = dydu * detA
        Tangent1(3) = dzdu * detA

        ! The 2nd tangent element is the normal vector of the parent element
        IF( PRESENT( Parent ) ) THEN
          pParent => Parent
        ELSE
          pParent => Boundary % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED(pParent) ) THEN
            pParent => Boundary % BoundaryInfo % Right
          END IF          
        END IF

        n = pParent % TYPE % NumberOfNodes
        ALLOCATE( ParentNodes % x(n), ParentNodes % y(n), ParentNodes % z(n) )        
        ParentNodes % x(1:n) = CurrentModel % Nodes % x(pParent % NodeIndexes)
        ParentNodes % y(1:n) = CurrentModel % Nodes % y(pParent % NodeIndexes)
        ParentNodes % z(1:n) = CurrentModel % Nodes % z(pParent % NodeIndexes)
        Tangent2 = NormalVector( pParent, ParentNodes) 
        DEALLOCATE( ParentNodes % x, ParentNodes % y, ParentNodes % z)
        
        Normal = CrossProduct( Tangent1, Tangent2 )         
      ELSE        
        IF( PRESENT( u0 ) ) THEN
          u = u0
        ELSE
          u = 0.0_dp
        END IF

        n = Boundary % TYPE % NumberOfNodes
        CALL NodalFirstDerivatives( n, dLBasisdx, Boundary, u, 0.0_dp, 0.0_dp )
        dxdu = DOT_PRODUCT( dLBasisdx(1:n,1), nx(1:n) )
        dydu = DOT_PRODUCT( dLBasisdx(1:n,1), ny(1:n) )

        detA = dxdu*dxdu + dydu*dydu
        IF ( detA <= 0._dp ) THEN
          Normal = 0._dp
          RETURN
        END IF
        detA = 1.0_dp / SQRT(detA)
        Normal(1) = -dydu * detA
        Normal(2) =  dxdu * detA
        Normal(3) =  0.0d0
      END IF
        
    CASE ( 2 ) 
      IF( PRESENT( u0 ) ) THEN
        u = u0
        v = v0
      ELSE
        IF( Boundary % TYPE % ElementCode / 100 == 3 ) THEN
          u = 1.0_dp/3
          v = 1.0_dp/3
        ELSE
          u = 0.0_dp
          v = 0.0_dp
        END IF
      END IF

      n = Boundary % TYPE % NumberOfNodes
      CALL NodalFirstDerivatives( n, dLBasisdx, Boundary, u, v, 0.0_dp )
      dxdu = DOT_PRODUCT( dLBasisdx(1:n,1), nx(1:n) )
      dydu = DOT_PRODUCT( dLBasisdx(1:n,1), ny(1:n) )
      dzdu = DOT_PRODUCT( dLBasisdx(1:n,1), nz(1:n) )
      dxdv = DOT_PRODUCT( dLBasisdx(1:n,2), nx(1:n) )
      dydv = DOT_PRODUCT( dLBasisdx(1:n,2), ny(1:n) )
      dzdv = DOT_PRODUCT( dLBasisdx(1:n,2), nz(1:n) )

      Auu = dxdu*dxdu + dydu*dydu + dzdu*dzdu
      Auv = dxdu*dxdv + dydu*dydv + dzdu*dzdv
      Avv = dxdv*dxdv + dydv*dydv + dzdv*dzdv

      detA = 1.0d0 / SQRT(Auu*Avv - Auv*Auv)

      Normal(1) = (dydu * dzdv - dydv * dzdu) * detA
      Normal(2) = (dxdv * dzdu - dxdu * dzdv) * detA
      Normal(3) = (dxdu * dydv - dxdv * dydu) * detA
    
    CASE DEFAULT
      CALL Fatal('NormalVector','No normal for '&
          //I2S(Boundary % TYPE % ElementCode)//' in '//I2S(meshDim)//'dim mesh!')
      
    END SELECT


    CheckParent = .FALSE.
    IF( PRESENT( Parent ) ) CheckParent = ASSOCIATED( Parent ) 
    
    CheckBody = .FALSE.
    IF ( PRESENT(Check) ) CheckBody = Check

    IF ( .NOT. ( CheckBody .OR. CheckParent ) ) RETURN
   
    SELECT CASE( Boundary % TYPE % ElementCode / 100 ) 

    CASE(1)
      x = nx(1)
      y = nx(1)
      z = nz(1)

    CASE(2,4)
      n = Boundary % TYPE % NumberOfNodes
      CALL NodalBasisFunctions( n, Basis, Boundary, 0.0d0, 0.0d0, 0.0d0 )
      x = DOT_PRODUCT( Basis(1:n), nx(1:n) )
      y = DOT_PRODUCT( Basis(1:n), ny(1:n) )
      z = DOT_PRODUCT( Basis(1:n), nz(1:n) )

    CASE(3)
      n = Boundary % TYPE % NumberOfNodes
      CALL NodalBasisFunctions( n, Basis, Boundary, 1.0d0/3, 1.0d0/3, 0.0d0 )
      x = DOT_PRODUCT( Basis(1:n), nx(1:n) )
      y = DOT_PRODUCT( Basis(1:n), ny(1:n) )
      z = DOT_PRODUCT( Basis(1:n), nz(1:n) )
    END SELECT

    IF( CheckParent ) THEN
      CALL CheckNormalDirectionParent( Boundary, Normal, x, y, z, Parent,Turn )   
    ELSE
      CALL CheckNormalDirection( Boundary,Normal,x,y,z,Turn )
    END IF

!------------------------------------------------------------------------------
  END FUNCTION NormalVector
!------------------------------------------------------------------------------

#if 0
!------------------------------------------------------------------------------
!> More economical normal vector computation assuming linear geometry description.
!------------------------------------------------------------------------------
  RECURSIVE FUNCTION NormalVectorLinear( Boundary,BoundaryNodes,Parent) RESULT(Normal)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Boundary
    TYPE(Nodes_t) :: BoundaryNodes
    TYPE(Element_t), POINTER, OPTIONAL :: Parent
    REAL(KIND=dp) :: Normal(3)
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    REAL(KIND=dp) :: vec0(3), vec1(3), vec2(3), vec3(3) 
    TYPE(Element_t), POINTER :: pParent
    INTEGER :: i,i1,i2,i3,i4,n,m,ElemDim,MeshDim
    
!------------------------------------------------------------------------------

    x => CurrentModel % Nodes % x
    y => CurrentModel % Nodes % y
    z => CurrentModel % Nodes % z

    IF( PRESENT( Parent ) ) THEN
      pParent => Parent
    ELSE IF( ASSOCIATED( Boundary % BoundaryInfo ) ) THEN
      pParent => Boundary % BoundaryInfo % Left
      IF(.NOT. ASSOCIATED(pParent) ) THEN
        pParent => Boundary % BoundaryInfo % Right
      END IF
    END IF

    ElemDim = Boundary % Type % Dimension 
    MeshDim = CurrentModel % Mesh % MeshDim 
    
    IF(ElemDim <= MeshDim-1 .OR. .NOT. (ASSOCIATED(pParent)) ) THEN
      SELECT CASE ( ElemDim ) 
        
      CASE ( 0 ) 
        Normal(1) = 1.0_dp
        Normal(2:3) = 0.0_dp

      CASE ( 1 )
        i1 = Boundary % NodeIndexes(1)
        i2 = Boundary % NodeIndexes(2)

        vec1(1) = x(i2) - x(i1)
        vec1(2) = y(i2) - y(i1)
        vec1(3) = 0.0_dp

        Normal(1) = -vec1(2)
        Normal(2) = vec1(1)
        Normal(3) = 0.0_dp

        Normal = Normal / SQRT(SUM(Normal**2))

      CASE( 2 ) 
        n = Boundary % TYPE % ElementCode / 100 

        i1 = Boundary % NodeIndexes(1)
        IF(n==4) THEN
          i2 = Boundary % NodeIndexes(2)
          i3 = Boundary % NodeIndexes(3)
          i4 = Boundary % NodeIndexes(4)
        ELSE
          i2 = Boundary % NodeIndexes(2)
          i3 = Boundary % NodeIndexes(3)
          i4 = i1
        END IF
        
        vec1(1) = x(i3) - x(i1)
        vec1(2) = y(i3) - y(i1)
        vec1(3) = z(i3) - z(i1)
        
        vec2(1) = x(i4) - x(i2)
        vec2(2) = y(i4) - y(i2)
        vec2(3) = z(i4) - z(i2)
          
        Normal = CrossProduct( vec1, vec2 )
        Normal = Normal / SQRT(SUM(Normal**2))

      CASE DEFAULT
        CALL Fatal('NormalVector','Invalid dimension for determining normal!')

      END SELECT
      
    ELSE 

      SELECT CASE ( ElemDim ) 
        
      CASE ( 0 )                
        i1 = pParent % NodeIndexes(1)
        i2 = pParent % NodeIndexes(2)
        
        Normal(1) = x(i2) - x(i1)
        Normal(2) = y(i2) - y(i1)
        Normal(3) = 0.0_dp

        Normal = Normal / SQRT(SUM(Normal**2))
        IF( i1 == Boundary % NodeIndexes(1) ) THEN
          Normal = -Normal
        END IF
                       
      CASE ( 1 )
        i1 = Boundary % NodeIndexes(1)
        i2 = Boundary % NodeIndexes(2)

        vec1(1) = x(i1)
        vec1(2) = y(i1)
        vec1(3) = z(i1)

        vec2(1) = x(i2)
        vec2(2) = y(i2)
        vec2(3) = z(i2)
               
        vec0 = vec1-vec2
        vec0 = vec0 / SQRT(SUM(vec0**2))
        
        n = pParent % TYPE % ElementCode / 100 

        vec2 = 0.0_dp
        DO i=1,n
          i3 = pParent % NodeIndexes(i)
          IF(i3 == i1 .OR. i3 == i2 ) CYCLE

          ! Vector stretching from edge center to the other nodes
          ! of the parent element. 
          vec2(1) = vec3(1) + x(i3) 
          vec3(1) = vec3(1) + x(i3) 
          vec3(1) = vec3(1) + x(i3) 
        END DO
        ! Subtract the average 
        vec3 = vec3 - (n-2)*(vec1+vec2)/2 
        
        ! Remove projection in the direction of the line
        Normal = vec3 - SUM(vec0*vec3)*vec0
        Normal = -Normal / SQRT(SUM(Normal**2))

      CASE( 2 ) 
        n = Boundary % TYPE % ElementCode / 100 
        
        i1 = Boundary % NodeIndexes(1)
        IF(n==4) THEN
          i2 = Boundary % NodeIndexes(2)
          i3 = Boundary % NodeIndexes(3)
          i4 = Boundary % NodeIndexes(4)
        ELSE
          i2 = Boundary % NodeIndexes(2)
          i3 = Boundary % NodeIndexes(3)
          i4 = i1
        END IF
          
        vec1(1) = x(i3) - x(i1)
        vec1(2) = y(i3) - y(i1)
        vec1(3) = z(i3) - z(i1)
        
        vec2(1) = x(i4) - x(i2)
        vec2(2) = y(i4) - y(i2)
        vec2(3) = z(i4) - z(i2)
          
        Normal = CrossProduct( vec1, vec2 )
        Normal = Normal / SQRT(SUM(Normal**2))

        m = pParent % TYPE % ElementCode / 100 
        vec1 = 0.0_dp
        vec2 = 0.0_dp
        DO i=1,m
          i1 = pParent % NodeIndexes(i)
          IF( ANY( Boundary % NodeIndexes == i1 ) ) THEN
            vec1(1) = vec1(1) + x(i1)
            vec1(2) = vec1(2) + y(i1)
            vec1(3) = vec1(3) + z(i1)            
          ELSE
            vec2(1) = vec2(1) + x(i1)
            vec2(2) = vec2(2) + y(i1)
            vec2(3) = vec2(3) + z(i1)
          END IF
        END DO

        vec1 = vec1 / n
        vec2 = vec2 / (m-n)

        IF( SUM( (vec1-vec2)*Normal ) < 0.0_dp ) THEN
          Normal = -Normal
        END IF
        
      CASE DEFAULT
        CALL Fatal('NormalVector','Invalid dimension for determining normal!')
        
      END SELECT
    END IF
      
!------------------------------------------------------------------------------
  END FUNCTION NormalVectorLinear
!------------------------------------------------------------------------------
#endif


  
!------------------------------------------------------------------------------
!> Returns a point that is most importantly supposed to be on the surface
!> For noncurved elements this may simply be the mean while otherwise
!> there may be a need to find the surface node using the local coordinates.
!> Hence the optional parameters. Typically the NormalVector and SurfaceVector
!> should be defined at the same position.
!------------------------------------------------------------------------------
  FUNCTION SurfaceVector( Boundary,BoundaryNodes,u,v ) RESULT(Surface)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Boundary
    TYPE(Nodes_t)   :: BoundaryNodes
    REAL(KIND=dp),OPTIONAL :: u,v
    REAL(KIND=dp) :: Surface(3)
!------------------------------------------------------------------------------
    REAL(KIND=dp), DIMENSION(:), POINTER :: nx,ny,nz
    REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES)
    INTEGER :: i,n
!------------------------------------------------------------------------------

    nx => BoundaryNodes % x
    ny => BoundaryNodes % y
    nz => BoundaryNodes % z
    n = Boundary % TYPE % NumberOfNodes

    IF( .NOT. PRESENT( u ) ) THEN
      Surface(1) = SUM( nx ) / n
      Surface(2) = SUM( ny ) / n
      Surface(3) = SUM( nz ) / n
    ELSE
      IF( Boundary % TYPE % DIMENSION == 1 ) THEN
        CALL NodalBasisFunctions( n, Basis, Boundary, u, 0.0_dp, 0.0_dp )
      ELSE
        CALL NodalBasisFunctions( n, Basis, Boundary, u, v, 0.0_dp )
      END IF
      Surface(1) = DOT_PRODUCT( Basis(1:n), nx(1:n) )
      Surface(2) = DOT_PRODUCT( Basis(1:n), ny(1:n) )
      Surface(3) = DOT_PRODUCT( Basis(1:n), nz(1:n) )
    END IF

!------------------------------------------------------------------------------
  END FUNCTION SurfaceVector
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> This subroutine tests where the intersection between the line defined by two 
!> points and a plane (or line) defined by a boundary element meet. There is
!> an intersection if ( 0 < Lambda < 1 ). Of all intersections the first one is 
!> that with the smallest positive lambda. 
!---------------------------------------------------------------------------
  FUNCTION LineFaceIntersection(FaceElement,FaceNodes,&
      Rinit,Rfin,u,v) RESULT ( Lambda )
!---------------------------------------------------------------------------
    TYPE(Nodes_t) :: FaceNodes
    TYPE(Element_t) :: FaceElement
    REAL(KIND=dp) :: Rinit(3),Rfin(3)
    REAL(KIND=dp),OPTIONAL :: u,v
    REAL(KIND=dp) :: Lambda

    REAL (KIND=dp) :: Surface(3),t1(3),t2(3),Normal(3),Rproj
    REAL (KIND=dp) :: Lambda0
    INTEGER :: third

    third = 3

100 CONTINUE

    ! For higher order elements this may be a necessity
    IF( PRESENT( u ) .AND. PRESENT(v) ) THEN
      Surface = SurfaceVector( FaceElement, FaceNodes, u, v )
      Normal = NormalVector( FaceElement, FaceNodes, u, v )

    ELSE IF( FaceElement % TYPE % DIMENSION == 2 ) THEN
      ! Any point known to be at the surface, even corner node
      Surface(1) = FaceNodes % x(1)
      Surface(2) = FaceNodes % y(1)
      Surface(3) = FaceNodes % z(1)

      ! Tangent vector, nor normalized to unity!
      t1(1) = FaceNodes % x(2) - Surface(1)
      t1(2) = FaceNodes % y(2) - Surface(2)
      t1(3) = FaceNodes % z(2) - Surface(3)

      t2(1) = FaceNodes % x(third) - Surface(1)
      t2(2) = FaceNodes % y(third) - Surface(2)
      t2(3) = FaceNodes % z(third) - Surface(3)

      ! Normal vector obtained from the cross product of tangent vectors
      ! This is not normalized to unity as value of lambda does not depend on its magnitude
      Normal(1) = t1(2)*t2(3) - t1(3)*t2(2)
      Normal(2) = t1(3)*t2(1) - t1(1)*t2(3)
      Normal(3) = t1(1)*t2(2) - t1(2)*t2(1)
    ELSE
      Surface(1) = FaceNodes % x(1)
      Surface(2) = FaceNodes % y(1)
      Surface(3) = 0.0_dp

      Normal(1) = Surface(2) - FaceNodes % y(2)
      Normal(2) = FaceNodes % x(2) - Surface(1)
      Normal(3) = 0.0_dp      
    END IF

    ! Project of the line to the face normal
    Rproj = SUM( (Rfin - Rinit) * Normal )
    
    IF( ABS( Rproj ) < TINY( Rproj ) ) THEN
      ! if the intersection cannot be defined make it an impossible one
      Lambda = -HUGE( Lambda ) 
    ELSE
      Lambda = SUM( ( Surface - Rinit ) * Normal ) / Rproj
    END IF

    IF( FaceElement % Type % NumberOfNodes == 4 ) THEN
      IF( third == 3 ) THEN
        third = 4
        Lambda0 = Lambda
        GOTO 100
      END IF
      IF( ABS( Lambda0 ) < ABS( Lambda) ) THEN
        Lambda = Lambda0 
      END IF
   END IF


  END FUNCTION LineFaceIntersection
  

!---------------------------------------------------------------------------
!> This subroutine performs a similar test as above using slightly different 
!> strategy.
!---------------------------------------------------------------------------
  FUNCTION LineFaceIntersection2(FaceElement,FaceNodes,Rinit,Rfin,Intersect) RESULT ( Lambda ) 

    TYPE(Nodes_t) :: FaceNodes
    TYPE(Element_t) :: FaceElement
    REAL(KIND=dp) :: Rinit(3), Rfin(3),Lambda
    LOGICAL :: Intersect
!----------------------------------------------------------------------------
    REAL (KIND=dp) :: A(3,3),B(3),C(3),Eps,Eps2,Eps3,detA,absA,ds
    INTEGER :: split, i, n, notriangles, triangle, ElemDim

    Eps = EPSILON( Eps )
    Eps2 = SQRT(TINY(Eps2))    
    Eps3 = 1.0d-12
    Lambda = -HUGE( Lambda )
    Intersect = .FALSE.
    ElemDim = FaceElement % TYPE % DIMENSION 

    ! Then solve the exact points of intersection from a 3x3 or 2x2 linear system
    !--------------------------------------------------------------------------
    IF( ElemDim == 2 ) THEN
      n = FaceElement % Type % NumberOfNodes
      ! In 3D rectangular faces are treated as two triangles
      IF( n == 4 .OR. n == 8 .OR. n == 9 ) THEN
        notriangles = 2
      ELSE
        notriangles = 1
      END IF

      DO triangle=1,notriangles
          
        A(1:3,1) = Rfin(1:3) - Rinit(1:3)
        
        IF(triangle == 1) THEN
          A(1,2) = FaceNodes % x(1) - FaceNodes % x(2)
          A(2,2) = FaceNodes % y(1) - FaceNodes % y(2)
          A(3,2) = FaceNodes % z(1) - FaceNodes % z(2)
        ELSE 
          A(1,2) = FaceNodes % x(1) - FaceNodes % x(4)
          A(2,2) = FaceNodes % y(1) - FaceNodes % y(4)
          A(3,2) = FaceNodes % z(1) - FaceNodes % z(4)
        END IF

        A(1,3) = FaceNodes % x(1) - FaceNodes % x(3)
        A(2,3) = FaceNodes % y(1) - FaceNodes % y(3)
        A(3,3) = FaceNodes % z(1) - FaceNodes % z(3)
        
        ! Check for linearly dependent vectors
        detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
             - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
             + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        absA = SUM(ABS(A(1,1:3))) * SUM(ABS(A(2,1:3))) * SUM(ABS(A(3,1:3))) 

        IF(ABS(detA) <= eps * absA + Eps2) CYCLE
!        print *,'detA',detA

        B(1) = FaceNodes % x(1) - Rinit(1)
        B(2) = FaceNodes % y(1) - Rinit(2)
        B(3) = FaceNodes % z(1) - Rinit(3)
        
        CALL InvertMatrix( A,3 )
        C(1:3) = MATMUL( A(1:3,1:3),B(1:3) )
        
        IF( ANY(C(2:3) < -Eps3) .OR. ANY(C(2:3) > 1.0_dp + Eps3 ) ) CYCLE
        IF( C(2)+C(3) > 1.0_dp + Eps3 ) CYCLE

        ! Relate the point of intersection to local coordinates
        !IF(corners < 4) THEN
        !  u = C(2)
        !  v = C(3)
        !ELSE IF(corners == 4 .AND. split == 0) THEN
        !  u = 2*(C(2)+C(3))-1
        !  v = 2*C(3)-1
        !ELSE 
        !  ! For the 2nd split of the rectangle the local coordinates switched
        !  v = 2*(C(2)+C(3))-1
        !  u = 2*C(3)-1        
        !END IF
        
        Intersect = .TRUE.
        Lambda = C(1)
        EXIT
 
      END DO
    ELSE
      ! In 2D the intersection is between two lines
      
      A(1:2,1) = Rfin(1:2) - Rinit(1:2)
      A(1,2) = FaceNodes % x(1) - FaceNodes % x(2)
      A(2,2) = FaceNodes % y(1) - FaceNodes % y(2)

      detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))

      ! Lines are almost parallel => no intersection possible
      IF(ABS(detA) <= eps * absA + Eps2) RETURN

      B(1) = FaceNodes % x(1) - Rinit(1)
      B(2) = FaceNodes % y(1) - Rinit(2)

      CALL InvertMatrix( A,2 )
      C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
     
      IF(C(2) < -Eps3 .OR. C(2) > 1.0_dp + Eps3 ) RETURN

      Intersect = .TRUE.
      Lambda = C(1)

!      u = -1.0d0 + 2.0d0 * C(2)

    END IF

!    IF(.NOT. Inside) RETURN

!    stat = ElementInfo( Element, FaceNodes, U, V, W, SqrtElementMetric, &
!        Basis, dBasisdx )
    
!    Weights(1:n) = Basis(1:n)
!    MaxInd = 1
!    DO i=2,n
!      IF(Weights(MaxInd) < Weights(i)) MaxInd = i
!    END DO

  END FUNCTION LineFaceIntersection2
  
 

!---------------------------------------------------------------------------
!> This subroutine computes the signed distance of a point from a surface.
!---------------------------------------------------------------------------
  FUNCTION PointFaceDistance(BoundaryElement,BoundaryNodes,&
      Coord,Normal,u0,v0) RESULT ( Dist )
!---------------------------------------------------------------------------
    TYPE(Nodes_t) :: BoundaryNodes
    TYPE(Element_t) :: BoundaryElement
    REAL(KIND=dp) :: Coord(3),Normal(3)
    REAL(KIND=dp),OPTIONAL :: u0,v0
    REAL(KIND=dp) :: Dist

    REAL (KIND=dp) :: Surface(3),t1(3),t2(3),u,v

    ! For higher order elements this may be a necessity
    IF( PRESENT( u0 ) .AND. PRESENT(v0) ) THEN
      u = u0
      v = v0
      Surface = SurfaceVector( BoundaryElement, BoundaryNodes, u, v )
    ELSE
      u = 0.0_dp
      v = 0.0_dp

      ! Any point known to be at the surface, even corner node
      Surface(1) = BoundaryNodes % x(1)
      Surface(2) = BoundaryNodes % y(1)
      Surface(3) = BoundaryNodes % z(1)
    END IF

    Normal = NormalVector( BoundaryElement, BoundaryNodes, u, v, .TRUE. )

    ! Project of the line to the face normal
    Dist = SUM( (Surface - Coord ) * Normal ) 
END FUNCTION PointFaceDistance



!------------------------------------------------------------------------------
!> Convert global coordinates x,y,z inside element to local coordinates
!> u,v,w of the element.
!> @todo Change to support p elements
!------------------------------------------------------------------------------
  SUBROUTINE GlobalToLocal( u,v,w,x,y,z,Element,ElementNodes )
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: x,y,z,u,v,w
    TYPE(Element_t) :: Element
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MaxIter = 50
    INTEGER :: i,n
    REAL(KIND=dp) :: r,s,t,delta(3),prevdelta(3),J(3,3),J1(3,2),det,swap,acc,err,scl,eps
    REAL(KIND=dp) :: Basis(MAX_ELEMENT_NODES), dLBasisdx(MAX_ELEMENT_NODES,3)
    LOGICAL :: Converged
!------------------------------------------------------------------------------

    u = 0._dp
    v = 0._dp
    w = 0._dp
    IF (Element % TYPE % DIMENSION==0) RETURN

    n = Element % TYPE % NumberOfNodes
    scl = MAXVAL(ElementNodes % x(1:n)) - MINVAL(ElementNodes % x(1:n)) + &
        MAXVAL(ElementNodes % y(1:n)) - MINVAL(ElementNodes % y(1:n)) + &
        MAXVAL(ElementNodes % z(1:n)) - MINVAL(ElementNodes % z(1:n))
        
    
    ! @todo Not supported yet
!   IF (ASSOCIATED(Element % PDefs)) THEN
!      CALL Fatal('GlobalToLocal','P elements not supported yet!')
!   END IF

    eps = EPSILON(eps)
    acc = eps * scl**2
    Converged = .FALSE.

     delta = 0._dp

!------------------------------------------------------------------------------
    DO i=1,Maxiter
!------------------------------------------------------------------------------
      CALL NodalBasisFunctions( n, Basis, Element, u, v, w )
      r = DOT_PRODUCT( Basis(1:n), ElementNodes % x(1:n) ) - x
      s = DOT_PRODUCT( Basis(1:n), ElementNodes % y(1:n) ) - y
      t = DOT_PRODUCT( Basis(1:n), ElementNodes % z(1:n) ) - z

      err = r**2 + s**2 + t**2

      IF ( err < acc ) THEN
        Converged = .TRUE.
        EXIT
      END IF

      prevdelta = delta
      delta = 0.d0

      CALL NodalFirstDerivatives( n, dLBasisdx, Element, u, v, w )

      SELECT CASE( Element % TYPE % DIMENSION )
      CASE(1)

        J(1,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % x(1:n) )
        J(2,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % y(1:n) )
        J(3,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % z(1:n) )

        det = SUM( J(1:3,1)**2 )
        delta(1) = (r*J(1,1)+s*J(2,1)+t*J(3,1))/det

      CASE(2)

         J(1,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % x(1:n) )
         J(1,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % x(1:n) )
         J(2,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % y(1:n) )
         J(2,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % y(1:n) )

        SELECT CASE( CoordinateSystemDimension() )
           CASE(3)
              J(3,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % z(1:n) )
              J(3,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % z(1:n) )

              delta(1) = r
              delta(2) = s
              delta(3) = t
              delta(1:2) = MATMUL( TRANSPOSE(J(1:3,1:2)), delta )
              r = delta(1)
              s = delta(2)

              J(1:2,1:2) = MATMUL( TRANSPOSE(J(1:3,1:2)), J(1:3,1:2) )
              delta(3)   = 0.0d0
         END SELECT

         CALL SolveLinSys2x2( J(1:2,1:2), delta(1:2), [ r, s] )

      CASE(3)
        J(1,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % x(1:n) )
        J(1,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % x(1:n) )
        J(1,3) = DOT_PRODUCT( dLBasisdx(1:n,3), ElementNodes % x(1:n) )

        J(2,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % y(1:n) )
        J(2,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % y(1:n) )
        J(2,3) = DOT_PRODUCT( dLBasisdx(1:n,3), ElementNodes % y(1:n) )

        J(3,1) = DOT_PRODUCT( dLBasisdx(1:n,1), ElementNodes % z(1:n) )
        J(3,2) = DOT_PRODUCT( dLBasisdx(1:n,2), ElementNodes % z(1:n) )
        J(3,3) = DOT_PRODUCT( dLBasisdx(1:n,3), ElementNodes % z(1:n) )

        CALL SolveLinSys3x3( J, delta, [ r, s, t ] )

      END SELECT

      IF( i > 10 ) THEN
        ! If the same values is suggested over and over again, then exit
        ! This may be a sign that the node is off-plane and cannot be 
        ! described within the element.
        IF( SUM( ABS( delta - prevdelta ) ) < eps ) EXIT

        ! Use sloppier criteria when iteration still unsuccessful
        IF( i > 20 ) THEN
          IF( SUM( ABS( delta - prevdelta ) ) < 1.0e-8 ) EXIT
        END IF

        ! If the iteration does not proceed try with some relaxation
        delta = 0.5_dp * delta 
      END IF

      u = u - delta(1)
      v = v - delta(2)
      w = w - delta(3)


!------------------------------------------------------------------------------
    END DO
!------------------------------------------------------------------------------

    IF ( .NOT. Converged ) THEN        
      IF( err > 1.0e-8 ) THEN
        IF( i > MaxIter ) THEN
          CALL Warn( 'GlobalToLocal', 'did not converge.')
          PRINT *,'rst',i,r,s,t
          PRINT *,'err',err,acc,eps
          PRINT *,'delta',delta,prevdelta
          PRINT *,'uvw',u,v,w
          PRINT *,'code',Element % TYPE % ElementCode
          PRINT *,'x:',x,ElementNodes % x(1:n)
          PRINT *,'y:',y,ElementNodes % y(1:n)
          PRINT *,'z:',z,ElementNodes % z(1:n)
        ELSE
!          CALL Warn( 'GlobalToLocal', 'Node may be out of element')
!          PRINT *,'rst',i,r,s,t,acc
        END IF
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GlobalToLocal
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
!>     Given element and its face map (for some triangular face of element ), 
!>     this routine returns global direction of triangle face so that 
!>     functions are continuous over element boundaries
!------------------------------------------------------------------------------
  FUNCTION getTriangleFaceDirection( Element, FaceMap, Indexes ) RESULT(globalDir)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Element_t) :: Element   !< Element to get direction to
    INTEGER :: FaceMap(3)        !< Element triangular face map
    INTEGER :: Indexes(:)
    INTEGER :: globalDir(3)      !< Global direction of triangular face as local node numbers.
!------------------------------------------------------------------------------
    INTEGER :: i, nodes(3)  
    
    ! Put global nodes of face into sorted order
    nodes(1:3) = Indexes( FaceMap )
    CALL sort(3, nodes)
    
    globalDir = 0
    ! Find local numbers of sorted nodes. These local nodes 
    ! span continuous functions over element boundaries
    DO i=1,Element % TYPE % NumberOfNodes
       IF (nodes(1) == Indexes(i)) THEN
          globalDir(1) = i
       ELSE IF (nodes(2) == Indexes(i)) THEN
          globalDir(2) = i
       ELSE IF (nodes(3) == Indexes(i)) THEN
          globalDir(3) = i
       END IF
    END DO
  END FUNCTION getTriangleFaceDirection


!------------------------------------------------------------------------------
!>     Given element and its face map (for some square face of element ), 
!>     this routine returns global direction of square face so that 
!>     functions are continuous over element boundaries
!------------------------------------------------------------------------------
  FUNCTION getSquareFaceDirection( Element, FaceMap, Indexes ) RESULT(globalDir)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element   !< Element to get direction to
    INTEGER :: FaceMap(:)        !< Element square face map
    INTEGER :: Indexes(:)
    INTEGER :: globalDir(4)      !< Global direction of square face as local node numbers.
!------------------------------------------------------------------------------
    INTEGER :: i, A,B,C,D, nodes(4), minGlobal

    ! Get global nodes 
    nodes(1:4) = Indexes( FaceMap )

    ! Find min global node
    minGlobal = nodes(1)
    A = 1
    DO i=2,4
       IF (nodes(i) < minGlobal) THEN
          A = i
          minGlobal = nodes(i)
       END IF
    END DO

    ! Now choose node B as the smallest node NEXT to min node
    B = MOD(A,4)+1
    C = MOD(A+3,4)
    IF (C == 0) C = 4
    D = MOD(A+2,4)
    IF (D == 0) D = 4
    IF (nodes(B) > nodes(C)) THEN
       i = B
       B = C
       C = i
    END IF

    ! Finally find local numbers of nodes A,B and C. They uniquely
    ! define a global face so that basis functions are continuous 
    ! over element boundaries
    globalDir = 0
    DO i=1,Element % TYPE % NumberOfNodes
       IF (nodes(A) == Indexes(i)) THEN
          globalDir(1) = i
       ELSE IF (nodes(B) == Indexes(i)) THEN
          globalDir(2) = i
       ELSE IF (nodes(C) == Indexes(i)) THEN
          globalDir(4) = i
       ELSE IF (nodes(D) == Indexes(i)) THEN
          globalDir(3) = i
       END IF
    END DO
  END FUNCTION getSquareFaceDirection


!------------------------------------------------------------------------------
!>     Function checks if given local numbering of a square face
!>     is legal for wedge element
!------------------------------------------------------------------------------
  FUNCTION wedgeOrdering( ordering ) RESULT(retVal)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, DIMENSION(4), INTENT(IN) :: ordering  !< Local ordering of a wedge square face
    LOGICAL :: retVal                              !< .TRUE. iff given ordering is legal for wedge square face.

    retVal = .FALSE.
    IF ((ordering(1) >= 1 .AND. ordering(1) <= 3 .AND.&
         ordering(2) >= 1 .AND. ordering(2) <= 3) .OR. &
       (ordering(1) >= 4 .AND. ordering(1) <= 6 .AND.&
       ordering(2) >= 4 .AND. ordering(2) <= 6)) THEN
       retVal = .TRUE.
    END IF
  END FUNCTION wedgeOrdering

  !---------------------------------------------------------
  !> Computes the 3D rotation matrix for a given 
  !> surface normal vector
  !---------------------------------------------------------
  FUNCTION ComputeRotationMatrix(PlaneVector) RESULT ( RotMat )

    REAL(KIND=dp) :: PlaneVector(3), RotMat(3,3), ex(3), ey(3), ez(3)
    INTEGER :: i, MinIndex, MidIndex, MaxIndex

    !Ensure PlaneVector is the unit normal
    PlaneVector = PlaneVector / SQRT( SUM(PlaneVector ** 2) )
    
    !The new z-axis is normal to the defined surface
    ez = PlaneVector

    MaxIndex = MAXLOC(ABS(ez),1)
    MinIndex = MINLOC(ABS(ez),1)

    !Special case when calving front perfectly aligned to either
    ! x or y axis. In this case, make minindex = 3 (ex points upwards)
    IF(ABS(ez(3)) == ABS(ez(2)) .OR. ABS(ez(3)) == ABS(ez(1))) &
         MinIndex = 3

    DO i=1,3
       IF(i == MaxIndex .OR. i == MinIndex) CYCLE
       MidIndex = i
    END DO

    ex(MinIndex) = 1.0
    ex(MidIndex) = 0.0
    
    ex(MaxIndex) = -ez(MinIndex)/ez(MaxIndex)
    ex = ex / SQRT( SUM(ex ** 2) )

    !The new y-axis is orthogonal to new x and z axes
    ey = CrossProduct(ez, ex)
    ey = ey / SQRT( SUM(ey ** 2) ) !just in case...

    RotMat(1,:) = ex
    RotMat(2,:) = ey
    RotMat(3,:) = ez

  END FUNCTION ComputeRotationMatrix



  ! Observe the cuts on one single element.
  ! This could be used to improve on the integration rules if we know where the
  ! element should be split.
  !-----------------------------------------------------------------------------
  SUBROUTINE CutSingleElement(Element, ElemNodes, ElemPhi, ElemCut )

    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: ElemNodes
    REAL(KIND=dp) :: ElemPhi(:)
    LOGICAL :: ElemCut(:)

    INTEGER :: i,i2,n
    REAL(KIND=dp) :: h1,h2,hprod,r
    REAL(KIND=dp), PARAMETER :: Eps=1.0e-3
    
    n = Element % TYPE % ElementCode / 100
    ElemCut(1:2*n) = .FALSE.
    
    h1 = MINVAL(ElemPhi(1:n))
    h2 = MAXVAL(ElemPhi(1:n))
    IF(h1*h2 >= 0.0_dp) RETURN

    IF( (SIZE(ElemNodes % x) < 2*n) ) THEN
      CALL Fatal('CutSingleElement','ElemNodes too small!')
    END IF
    
    DO i=1, n
      i2 = MODULO(i,n)+1
      h1 = ElemPhi(i)
      h2 = ElemPhi(i2)
      hprod = h1*h2            

      ! First mark the cut nodes.        
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))        
        IF( r <= Eps ) THEN
          ElemCut(i2) = .TRUE.
        ELSE IF((1.0-r < Eps) ) THEN
          ElemCut(i) = .TRUE.
        ELSE
          ElemCut(n+i) = .TRUE.

          ! We update nodes so that the element on-the-fly can point to then using NodeIndexes. 
          ElemNodes % x(n+i) = (1-r) * ElemNodes % x(i2) + r * ElemNodes % x(i)
          ElemNodes % y(n+i) = (1-r) * ElemNodes % y(i2) + r * ElemNodes % y(i)
          ElemNodes % z(n+i) = (1-r) * ElemNodes % z(i2) + r * ElemNodes % z(i)
        END IF
      ELSE IF( ABS(hprod) < 1.0d-20 ) THEN
        IF(ABS(h1) < 1.0e-20) ElemCut(i) = .TRUE. 
        IF(ABS(h2) < 1.0e-20) ElemCut(i2) = .TRUE.
      END IF
    END DO

  END SUBROUTINE CutSingleElement


  ! Given a single element and a list of node and edge cuts create a list of
  ! pieces coming from the split.
  !---------------------------------------------------------------------------
  SUBROUTINE SplitSingleElement(Element, ElemCut, ElemNodes, m, &
      IsCut, IsMore, LocalInds, SgnNode )

    TYPE(Element_t) :: Element
    LOGICAL :: ElemCut(:)
    TYPE(Nodes_t) :: ElemNodes
    INTEGER :: m
    LOGICAL :: IsCut, IsMore
    INTEGER :: LocalInds(:)
    INTEGER :: SgnNode
    
    
    INTEGER :: n,n_split,n_cut,ElemType,SplitCase,iCase,subcase
    INTEGER :: j,j2,j3,j4,mmax
    REAL(KIND=dp) :: s1,s2
    
    SAVE :: subcase, j, j2, j3, j4, mmax, s1, s2
    !$OMP THREADPRIVATE(subcase, j, j2, j3, j4, mmax, s1, s2)

    ElemType = Element % TYPE % ElementCode
    n = ElemType / 100
        
    n_split = COUNT( ElemCut(n+1:2*n) )
    n_cut = COUNT( ElemCut(1:n) )
    
    IsMore = .FALSE.   
    IsCut = (n_split > 0)

    ! Nothing to do, use original element.
    IF(.NOT. IsCut) RETURN

    ! This allows use case to deal with element types, edge splits and node splits at the same time. 
    ! It is a matter of taste if this is ok or not...
    SplitCase = 100 * ElemType + 10 * n_split + n_cut
    iCase = 0
    LocalInds = 0
    
    SELECT CASE( SplitCase ) 

      
    CASE( 30320, 30321 ) 
      ! Triangle being cut on two edges.
      IF( m == 1 ) THEN
        ! Find the only edge that is not cut
        DO j=1,3
          IF( .NOT. ElemCut( n + j ) ) EXIT
        END DO
        j2 = MODULO(j,3)+1
        j3 = MODULO(j+1,3)+1
        mmax = 3
        
        ! There are two ways to split the triangle.
        ! Choose the one with shorter diameter.
        s1 = (ElemNodes % x(j) - ElemNodes % x(n + j2))**2 + &
            (ElemNodes % y(j) - ElemNodes % y(n + j2))**2 + &
            (ElemNodes % z(j) - ElemNodes % z(n + j2))**2 
        s2 = (ElemNodes % x(j2) - ElemNodes % x(n + j3))**2 + &
            (ElemNodes % y(j2) - ElemNodes % y(n + j3))**2 + &
            (ElemNodes % z(j2) - ElemNodes % z(n + j3))**2 

        LocalInds(1) = j
        LocalInds(2) = j2                 
        IF( s1 < s2 ) THEN
          LocalInds(3) = n + j2
        ELSE
          LocalInds(3) = n + j3
        END IF
        SgnNode = 1
        iCase = 1
      ELSE IF(m==2) THEN
        IF( s1 < s2 ) THEN
          LocalInds(1) = j
        ELSE
          LocalInds(1) = j2       
        END IF
        LocalInds(2) = n + j2
        LocalInds(3) = n + j3

        SgnNode = 1
        iCase = 2
      ELSE IF(m==3) THEN
        LocalInds(1) = n + j3
        LocalInds(2) = n + j2
        LocalInds(3) = j3

        SgnNode = 3
        iCase = 3
      END IF

    CASE( 30311 ) 
      ! Triangle being cut on one edge and one node. 
      IF( m == 1 ) THEN
        ! Find the only edge that is cut
        DO j=1,3
          IF( ElemCut( n + j ) ) EXIT
        END DO
        j2 = MODULO(j,3)+1
        j3 = MODULO(j+1,3)+1
      END IF
      
      ! One cut result to splitted elements only if the opposing node is cut through
      IF( ElemCut(j3) ) THEN
        IF(m==1) THEN
          LocalInds(1) = n + j
          LocalInds(2) = j2
          LocalInds(3) = j3
          
          SgnNode = 2
          iCase = 4
          mmax = 2
        ELSE IF(m==2) THEN
          LocalInds(1) = n + j
          LocalInds(2) = j3
          LocalInds(3) = j
          
          sgnNode = 3
          iCase = 5
        END IF
      ELSE IF(ElemCut(j) .OR. ElemCut(j2)) THEN
        LocalInds(1:3) = [1,2,3]          
        
        iCase = 6
        SgnNode = j3          
        mmax = 1
      END IF

    CASE( 40420, 40421 ) 
      ! Quadrilateral being cut on two edges. 
      
      IF( m == 1 ) THEN
        subcase = 0
        IF( ElemCut(n+1) .AND. ElemCut(n+3) ) THEN
          subcase = 1
          j = 1
          mmax = 2
        ELSE IF( ElemCut(n+2) .AND. ElemCut(n+4) ) THEN
          subcase = 1
          j = 2
          mmax = 2
        ELSE
          DO j=1,4
            j2 = MODULO(j,4)+1
            IF( ElemCut(n+j) .AND. ElemCut(n+j2) ) THEN
              subcase = 2 
              mmax = 3
              EXIT
            END IF
          END DO
        END IF
        IF( subcase == 0 ) THEN
          CALL Fatal('SplitSingleElement','This case not treated yet for 404!')
        END IF
      END IF

      
      IF( subcase == 1 ) THEN        
        mmax = 2
        
        IF( m == 1 ) THEN
          j2 = MODULO(j,4)+1
          j3 = MODULO(j+1,4)+1
          j4 = MODULO(j+2,4)+1
          
          LocalInds(1) = j
          LocalInds(2) = n + j
          LocalInds(3) = n + j3
          LocalInds(4) = j4          
          
          SgnNode = 1
          iCase = 7
        ELSE IF(m==2) THEN
          LocalInds(1) = j2
          LocalInds(2) = j3
          LocalInds(3) = n + j3
          LocalInds(4) = n + j
          
          SgnNode = 1
          iCase = 8
        END IF

      ELSE IF( subcase == 2 ) THEN
        mmax = 4

        IF( m == 1 ) THEN
          j2 = MODULO(j,4)+1
          j3 = MODULO(j+1,4)+1
          j4 = MODULO(j+2,4)+1

          LocalInds(1) = n + j
          LocalInds(2) = j2
          LocalInds(3) = n + j2

          SgnNode = 2
          iCase = 9
        ELSE IF(m==2) THEN
          LocalInds(1) = j
          LocalInds(2) = n + j
          LocalInds(3) = j4

          SgnNode = 3
          iCase = 10
        ELSE IF(m==3) THEN
          LocalInds(1) = n + j
          LocalInds(2) = n + j2
          LocalInds(3) = j4

          SgnNode = 3
          iCase = 11
        ELSE IF(m==4) THEN
          LocalInds(1) = n + j2
          LocalInds(2) = j3
          LocalInds(3) = j4

          SgnNode = 3
          iCase = 12
        END IF

      END IF

    CASE( 40411 ) 
      ! Quadrilateral being cut on one edge and one node.  

      ! Find the only edge that is cut
      DO j=1,4
        IF( ElemCut( n + j ) ) EXIT
      END DO
      j2 = MODULO(j,4)+1
      j3 = MODULO(j+1,4)+1
      j4 = MODULO(j+2,4)+1

      ! IF we cut node associated to the same edge, we don't really have a split element,
      IF(ElemCut(j) .OR. ElemCut(j2)) THEN
        LocalInds(1:4) = [1,2,3,4]
        iCase = 13
        SgnNode = j3          
        mmax = 1
      ELSE
        mmax = 2
        IF( ElemCut(j3) ) THEN
          IF(m==1) THEN
            LocalInds(1) = n + j
            LocalInds(2) = j2
            LocalInds(3) = j3
            LocalInds(4) = j4

            iCase = 14
            SgnNode = 3
          ELSE IF(m==2) THEN
            LocalInds(1) = j
            LocalInds(2) = n + j
            LocalInds(3) = j4

            iCase = 15
            SgnNode = 1
          END IF

        ELSE IF( ElemCut(j4)) THEN
          IF(m==1) THEN
            LocalInds(1) = j
            LocalInds(2) = n + j
            LocalInds(3) = j3
            LocalInds(4) = j4

            iCase = 16
            SgnNode = 4
          ELSE IF(m==2) THEN
            LocalInds(1) = n + j
            LocalInds(2) = j2
            LocalInds(3) = j3

            iCase = 17
            SgnNode = 2
          END IF
        END IF
      END IF
      
    CASE DEFAULT
      PRINT *,'ElemCut:',ElemCut(1:n*n)
      CALL Fatal('SplitSingleElement','Unknown split case in element divisions: '//I2S(SplitCase))
    END SELECT

    IsMore = (m < mmax )
    !IF(iCase>0) nCase(iCase) = nCase(iCase) + 1
    
  END SUBROUTINE SplitSingleElement

END MODULE ElementGeometry
