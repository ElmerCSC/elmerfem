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
! *  Authors: Juha Ruokolainen, Ville Savolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 01 Oct 1996
! *
! ****************************************************************************/

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Module containing interpolation and quadrant tree routines.
!------------------------------------------------------------------------------

MODULE Interpolation

   USE Types
   USE SParIterGlobals
   USE CoordinateSystems
   USE ElementDescription, ONLY : GlobalToLocal, ElementInfo
   USE PElementMaps, ONLY : isActivePElement
   USE Integration, ONLY : GaussIntegrationPoints_t, GaussPoints

   IMPLICIT NONE
   
 CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE FindLeafElements( Point, dim, RootQuadrant, LeafQuadrant )
!------------------------------------------------------------------------------
     REAL(KIND=dp), DIMENSION(3) :: Point
     REAL(KIND=dp), DIMENSION(6) :: GeometryBoundingBox
     INTEGER :: dim
     TYPE(Quadrant_t), POINTER :: RootQuadrant, LeafQuadrant
!------------------------------------------------------------------------------

     LeafQuadrant => RootQuadrant
     GeometryBoundingBox = RootQuadrant % BoundingBox
!
!    Find recursively the last generation
!    quadrant the point belongs to:
!    -------------------------------------
     CALL FindPointsQuadrant(Point, dim, LeafQuadrant)

   CONTAINS

!------------------------------------------------------------------------------
     RECURSIVE SUBROUTINE FindPointsQuadrant(Point, dim, MotherQuadrant)
!------------------------------------------------------------------------------

     REAL(KIND=dp), DIMENSION(3) :: Point
     INTEGER :: dim
     TYPE(Quadrant_t), POINTER :: MotherQuadrant
!------------------------------------------------------------------------------
     TYPE(Quadrant_t), POINTER :: ChildQuadrant
     INTEGER :: i
     REAL(KIND=dp) :: BBox(6), eps3
     REAL(KIND=dp), PARAMETER :: eps2=0.0_dp !!!!!!! *** !!!!!!
!------------------------------------------------------------------------------
     
!    Loop over ChildQuadrants:
!    -------------------------
     DO i=1, 2**dim
        ChildQuadrant => MotherQuadrant % ChildQuadrants(i) % Quadrant
        BBox = ChildQuadrant % BoundingBox
!
!       ******** NOTE: eps2 set to zero at the moment **********
!
        eps3 = eps2 * MAXVAL(BBox(4:6) - BBox(1:3))
        BBox(1:3) = BBox(1:3) - eps3
        BBox(4:6) = BBox(4:6) + eps3
        !
        ! Is the point in ChildQuadrant(i)?
        ! ----------------------------------
        IF ( ( Point(1) >= BBox(1)) .AND. (Point(1) <= BBox(4) ) .AND. &
             ( Point(2) >= BBox(2)) .AND. (Point(2) <= BBox(5) ) .AND. &
             ( Point(3) >= BBox(3)) .AND. (Point(3) <= BBox(6) ) ) EXIT
     END DO

     IF ( i > 2**dim ) THEN
!       PRINT*,'Warning: point not found in any of the quadrants ?'
        NULLIFY( MotherQuadrant )
        RETURN
     END IF

     MotherQuadrant => ChildQuadrant
!
!    Are we already in the LeafQuadrant ?
!    If not, search for the next generation
!    ChildQuadrants for the point:
!    ---------------------------------------
     IF ( ASSOCIATED ( MotherQuadrant % ChildQuadrants ) )THEN
        CALL FindPointsQuadrant( Point, dim, MotherQuadrant )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE FindPointsQuadrant
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE FindLeafElements
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!>    Checks whether a given point belongs to a given bulk element.
!>    If it does, returns the local coordinates in the bulk element
!------------------------------------------------------------------------------
     FUNCTION PointInElement( Element, ElementNodes, Point, &
	  LocalCoordinates, GlobalEps, LocalEps, NumericEps, &
          GlobalDistance, LocalDistance, EdgeBasis, &
          USolver ) RESULT(IsInElement)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element  !< Bulk element we are checking
    TYPE(Nodes_t) :: ElementNodes        !< The nodal points of the bulk element
    LOGICAL :: IsInElement               !< Whether the node lies within the element
    REAL(KIND=dp), DIMENSION(:) :: Point      !< Point under study.
    REAL(KIND=dp), DIMENSION(:) :: LocalCoordinates  !< Local coordinates corresponding to the global ones.
    REAL(KIND=dp), OPTIONAL :: GlobalEps !< Required accuracy of global coordinates
    REAL(KIND=dp), OPTIONAL :: LocalEps  !< Required accuracy of local coordinates
    REAL(KIND=dp), OPTIONAL :: NumericEps !< Accuracy of numberical operations
    REAL(KIND=dp), OPTIONAL :: GlobalDistance !< Returns the distance from the element in global coordinates.
    REAL(KIND=dp), OPTIONAL :: LocalDistance  !< Returns the distance from the element in local coordinates.
    LOGICAL, OPTIONAL :: EdgeBasis
    TYPE(Solver_t), POINTER, OPTIONAL :: USolver 
!------------------------------------------------------------------------------
    INTEGER :: n
    INTEGER :: i
    LOGICAL :: ComputeDistance, trans
    REAL(KIND=dp) :: ug,vg,wg,xdist,ydist,zdist,sumdist,eps0,eps1,eps2,escale,&
        minx,maxx,miny,maxy,minz,maxz
!------------------------------------------------------------------------------


!   Initialize:
!   -----------

    IsInElement = .FALSE.
    n = Element % TYPE % NumberOfNodes

    ! The numeric precision 
    IF ( PRESENT(NumericEps) ) THEN
      eps0 = NumericEps
    ELSE
      eps0 = EPSILON( eps0 )
    END IF
    
    ! The rough check, used for global coordinates
    IF ( PRESENT(GlobalEps) ) THEN
      Eps1 = GlobalEps
    ELSE
      Eps1 = 1.0e-4
    END IF 
    
    ! The more detailed condition, used for local coordinates
    IF ( PRESENT(LocalEps) ) THEN
      Eps2 = LocalEps
    ELSE
      Eps2 = 1.0e-10
    END IF

    IF( PRESENT( LocalDistance ) ) THEN
      LocalDistance = HUGE( LocalDistance ) 
    END IF
    
   IF( Eps1 < 0.0_dp ) THEN
      CONTINUE
   ELSE IF( PRESENT( GlobalDistance ) ) THEN
      ! When distance has to be computed all coordinate directions need to be checked
      
      minx = MINVAL( ElementNodes % x(1:n) )
      maxx = MAXVAL( ElementNodes % x(1:n) )

      miny = MINVAL( ElementNodes % y(1:n) )
      maxy = MAXVAL( ElementNodes % y(1:n) )

      minz = MINVAL( ElementNodes % z(1:n) )
      maxz = MAXVAL( ElementNodes % z(1:n) )
      
      xdist = MAX( MAX( Point(1) - maxx, 0.0_dp ), minx - Point(1) )
      ydist = MAX( MAX( Point(2) - maxy, 0.0_dp ), miny - Point(2) )
      zdist = MAX( MAX( Point(3) - maxz, 0.0_dp ), minz - Point(3) )
      
      GlobalDistance = SQRT( xdist**2 + ydist**2 + zdist**2)
      
      IF( xdist > eps0 + eps1 * (maxx - minx) ) RETURN 
      IF( ydist > eps0 + eps1 * (maxy - miny) ) RETURN 
      IF( zdist > eps0 + eps1 * (maxz - minz) ) RETURN 
    ELSE
      ! Otherwise make decision independently after each coordinate direction
      
      minx = MINVAL( ElementNodes % x(1:n) )
      maxx = MAXVAL( ElementNodes % x(1:n) )
      xdist = MAX( MAX( Point(1) - maxx, 0.0_dp ), minx - Point(1) )
      IF( xdist > eps0 + eps1 * (maxx - minx) ) RETURN 
      
      miny = MINVAL( ElementNodes % y(1:n) )
      maxy = MAXVAL( ElementNodes % y(1:n) )
      ydist = MAX( MAX( Point(2) - maxy, 0.0_dp ), miny - Point(2) )
      IF( ydist > eps0 + eps1 * (maxy - miny) ) RETURN 
      
      minz = MINVAL( ElementNodes % z(1:n) )
      maxz = MAXVAL( ElementNodes % z(1:n) )
      zdist = MAX( MAX( Point(3) - maxz, 0.0_dp ), minz - Point(3) )
      IF( zdist > eps0 + eps1 * (maxz - minz) ) RETURN 
    END IF

!   Get element local coordinates from global
!   coordinates of the point:
!   -----------------------------------------
    CALL GlobalToLocal( ug, vg, wg, Point(1), Point(2), Point(3), &
        Element, ElementNodes )

! Currently the eps of global coordinates is mixed with the eps of local
! coordinates which is a bit disturbin. There could be sloppier global
! coordinate search and a more rigorous local coordinate search.

    SELECT CASE ( Element % TYPE % ElementCode / 100 )
      CASE(1)
        sumdist = ug

      CASE(2)
        sumdist = MAX( ug - 1.0, MAX( -ug - 1.0, 0.0 ) )

      CASE(3)
        sumdist = MAX( -ug, 0.0 ) + MAX( -vg, 0.0 ) 
        sumdist = sumdist + MAX( ug + vg - 1.0_dp, 0.0 ) 

      CASE(4)
        sumdist = MAX( ug - 1.0, MAX( -ug -1.0, 0.0 ) )
        sumdist = sumdist + MAX( vg - 1.0, MAX( -vg - 1.0, 0.0 ) )

      CASE(5)
        sumdist = MAX( -ug, 0.0 ) + MAX( -vg, 0.0 ) + MAX( -wg, 0.0 ) 
        sumdist = sumdist + MAX( ug + vg + wg - 1.0, 0.0 ) 
        
      CASE(7)
        sumdist = MAX( -ug, 0.0 ) + MAX( -vg, 0.0 ) 
        sumdist = sumdist + MAX( ug + vg - 1.0_dp, 0.0 ) 
        sumdist = sumdist + MAX( wg - 1.0, MAX( -wg - 1.0, 0.0 ) )

      CASE(8)
        sumdist = MAX( ug - 1.0, MAX( -ug -1.0, 0.0 ) )
        sumdist = sumdist + MAX( vg - 1.0, MAX( -vg - 1.0, 0.0 ) )
        sumdist = sumdist + MAX( wg - 1.0, MAX( -wg - 1.0, 0.0 ) )
        
      CASE DEFAULT
        WRITE( Message,'(A,I4)') 'Not implemented for element code',&
            Element % TYPE % ElementCode 
        CALL Warn('PointInElement',Message)
      END SELECT


      IF( sumdist < eps0 + eps2 ) THEN
        IsInElement = .TRUE.
      END IF

      IF( PRESENT( LocalDistance ) ) THEN
        LocalDistance = sumdist
      END IF
        

    trans = PRESENT(EdgeBasis)
    IF(trans) trans=EdgeBasis
    trans = trans .OR. isActivePElement(Element,USolver)

    IF (trans) THEN
      SELECT CASE(Element % Type % ElementCode/100)
      CASE(3)
         ug = 2*ug + vg - 1
         vg = SQRT(3._dp)*vg
      CASE(5)
         ug = 2*ug + vg + wg - 1
         vg = SQRT(3._dp)*vg + 1/SQRT(3._dp)*wg
         wg = 2*SQRT(2/3._dp)*wg
      CASE(6)
         wg = SQRT(2._dp)*wg
      CASE(7)
         ug = 2*ug + vg - 1
         vg = SQRT(3._dp)*vg
      END SELECT
    END IF

    LocalCoordinates(1) = ug
    LocalCoordinates(2) = vg
    LocalCoordinates(3) = wg
!------------------------------------------------------------------------------
    END FUNCTION PointInElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>    Builds a tree hierarchy recursively bisectioning the geometry
!>    bounding box, and partitioning the bulk elements in the
!>    last level of the tree hierarchy
!------------------------------------------------------------------------------
    SUBROUTINE BuildQuadrantTree(Mesh, BoundingBox, RootQuadrant)
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh                        !< Finite element mesh
    REAL(KIND=dp), DIMENSION(6) :: BoundingBox  !< XMin, YMin, ZMin, XMax, YMax, ZMax
    TYPE(Quadrant_t), POINTER :: RootQuadrant   !< Quadrant tree structure root
!------------------------------------------------------------------------------
    INTEGER :: dim, Generation, i
    REAL(KIND=dp) :: XMin, XMax, YMin, YMax, ZMin, ZMax
    TYPE(Quadrant_t), POINTER :: MotherQuadrant
    INTEGER :: MaxLeafElems

    dim = MAX( Mesh % MeshDim, CoordinateSystemDimension() )
        
    IF ( dim == 3 ) THEN
      MaxLeafElems = 16
    ELSE
      MaxLeafElems = 8
    END IF

    Generation = 0

    XMin = BoundingBox(1)
    XMax = BoundingBox(4)
    IF ( dim >= 2 ) THEN
      YMin = BoundingBox(2)
      YMax = BoundingBox(5)
    ELSE
      YMin = 0.d0
      YMax = 0.d0
    END IF
    IF ( dim == 3) THEN
      ZMin = BoundingBox(3)
      ZMax = BoundingBox(6)
    ELSE
      ZMin = 0.d0
      ZMax = 0.d0
    END IF

! Create Mother of All Quadrants
    ALLOCATE ( RootQuadrant )

    RootQuadrant % BoundingBox = [ XMin, YMin, ZMin, XMax, YMax, ZMax ]
    RootQuadrant % NElemsInQuadrant = Mesh % NumberOfBulkElements

    ALLOCATE ( RootQuadrant % Elements( Mesh % NumberOfBulkElements ) )
    RootQuadrant % Elements = [ (i, i=1,Mesh % NumberOfBulkElements) ]

! Start building the quadrant tree
    CALL Info( 'BuildQuandrantTree', 'Start', Level=4 )
    MotherQuadrant => RootQuadrant
    CALL CreateChildQuadrants( MotherQuadrant, dim )
    CALL Info( 'BuildQuandrantTree', 'Ready', Level=4 )

  CONTAINS

!-------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE CreateChildQuadrants( MotherQuadrant, dim )
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Model is automatically available (internal subroutine)
!-------------------------------------------------------------------------------
    TYPE(Quadrant_t), POINTER :: MotherQuadrant
    INTEGER :: i, dim, n
    TYPE(QuadrantPointer_t) :: ChildQuadrant(8)
    REAL(KIND=dp) :: XMin, XMax, YMin, YMax, ZMin, ZMax

!-------------------------------------------------------------------------------
! Create 2**dim child quadrants
!-------------------------------------------------------------------------------
    n = 2**dim
    ALLOCATE ( MotherQuadrant % ChildQuadrants(n) )
    DO i=1, n
       ALLOCATE( MotherQuadrant % ChildQuadrants(i) % Quadrant )
       ChildQuadrant(i) % Quadrant => &
            MotherQuadrant % ChildQuadrants(i) % Quadrant
       ChildQuadrant(i) % Quadrant % NElemsInQuadrant = 0
       NULLIFY ( ChildQuadrant(i) % Quadrant % Elements )
       NULLIFY ( ChildQuadrant(i) % Quadrant % ChildQuadrants )
    END DO
!-------------------------------------------------------------------------------
    XMin = MotherQuadrant % BoundingBox(1)
    YMin = MotherQuadrant % BoundingBox(2)
    ZMin = MotherQuadrant % BoundingBox(3)

    XMax = MotherQuadrant % BoundingBox(4)
    YMax = MotherQuadrant % BoundingBox(5)
    ZMax = MotherQuadrant % BoundingBox(6)
    MotherQuadrant % Size = MAX ( MAX( XMax-XMin, YMax-YMin), ZMax-ZMin )

    ChildQuadrant(1) % Quadrant % BoundingBox = [ XMin, YMin, ZMin, &
      (XMin + XMax)/2.d0, (YMin + YMax)/2.d0, (ZMin + ZMax)/2.d0 ]

    ChildQuadrant(2) % Quadrant % BoundingBox = [ (XMin+XMax)/2.d0, &
        YMin, ZMin, XMax, (YMin+YMax)/2.d0, (ZMin+ZMax)/2.d0 ]

    IF ( dim >= 2 ) THEN
       ChildQuadrant(3) % Quadrant % BoundingBox = [ XMin, (YMin+YMax)/2.d0, &
               ZMin, (XMin+XMax)/2.d0, YMax, (ZMin+ZMax)/2.d0 ]

       ChildQuadrant(4) % Quadrant % BoundingBox = [ (XMin+XMax)/2.d0, &
           (YMin+YMax)/2.d0, ZMin, XMax, YMax, (ZMin+ZMax)/2.d0 ]
    END IF

    IF ( dim == 3 ) THEN
       ChildQuadrant(5) % Quadrant % BoundingBox = [ XMin, YMin, &
          (ZMin+ZMax)/2.d0, (XMin+XMax)/2.d0, (YMin+YMax)/2.d0, ZMax ]

       ChildQuadrant(6) % Quadrant % BoundingBox = [ (XMin+XMax)/2.d0, YMin, &
               (ZMin+ZMax)/2.d0, XMax, (YMin+YMax)/2.d0, ZMax ]

       ChildQuadrant(7) % Quadrant % BoundingBox = [ XMin, (YMin+YMax)/2.d0, &
               (ZMin+ZMax)/2.d0, (XMin+XMax)/2.d0, YMax, ZMax ]

       ChildQuadrant(8) % Quadrant % BoundingBox = [ (XMin+XMax)/2.d0, &
             (YMin+YMax)/2.d0, (ZMin+ZMax)/2.d0, XMax, YMax, ZMax ]
    END IF
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Loop over all elements in the mother quadrant,
! placing them in one of the 2^dim child quadrants
!-------------------------------------------------------------------------------
    CALL PutElementsInChildQuadrants( ChildQuadrant, MotherQuadrant, dim )

!-------------------------------------------------------------------------------
! Check whether we need to branch for the next level
!-------------------------------------------------------------------------------
    DO i=1,n
       ChildQuadrant(i) % Quadrant % Size = MotherQuadrant % Size / 2
       IF ( ChildQuadrant(i) % Quadrant % NElemsInQuadrant > MaxLeafElems ) THEN
          IF ( ChildQuadrant(i) % Quadrant % Size > &
                    ChildQuadrant(i) % Quadrant % MinElementSize ) THEN
             IF ( Generation <= 8 ) THEN
                Generation = Generation + 1
                CALL CreateChildQuadrants( ChildQuadrant(i) % Quadrant, dim )
                Generation = Generation - 1
             END IF
          END IF
       END IF
    END DO

    DEALLOCATE ( MotherQuadrant % Elements )
    NULLIFY ( MotherQuadrant % Elements )
!-------------------------------------------------------------------------------
    END SUBROUTINE CreateChildQuadrants
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Loop over all elements in the MotherQuadrant, placing them
!> in one of the 2^dim child quadrants.
!-------------------------------------------------------------------------------
    RECURSIVE SUBROUTINE PutElementsInChildQuadrants( ChildQuadrant, &
                   MotherQuadrant, dim )
!-------------------------------------------------------------------------------
      TYPE(QuadrantPointer_t) :: ChildQuadrant(8)
      TYPE(Quadrant_t), POINTER :: MotherQuadrant
      INTEGER :: dim
      REAL(KIND=dp) :: eps3
      REAL(KIND=dp), PARAMETER :: eps2=2.5d-2
!-------------------------------------------------------------------------------

      TYPE(Element_t), POINTER :: CurrentElement
      INTEGER :: i, j, t, n
      INTEGER, POINTER :: NodeIndexes(:)
      INTEGER :: ElementList(2**dim, MotherQuadrant % NElemsInQuadrant)

      LOGICAL :: ElementInQuadrant
      REAL(KIND=dp) :: BBox(6), XMin, XMax, YMin, YMax, ZMin, ZMax, ElementSize

!-------------------------------------------------------------------------------

      DO i=1,2**dim
         ChildQuadrant(i) % Quadrant % NElemsInQuadrant = 0
         ChildQuadrant(i) % Quadrant % MinElementSize   = 1.0d20
      END DO

!-------------------------------------------------------------------------------
      DO t=1, MotherQuadrant % NElemsInQuadrant
!-------------------------------------------------------------------------------
         CurrentElement => Mesh % Elements( MotherQuadrant % Elements(t) )
         n = CurrentElement % TYPE % NumberOfNodes
         NodeIndexes => CurrentElement % NodeIndexes

! Get element coordinates
         XMin = MINVAL( Mesh % Nodes % x(NodeIndexes) )
         XMax = MAXVAL( Mesh % Nodes % x(NodeIndexes) )
         YMin = MINVAL( Mesh % Nodes % y(NodeIndexes) )
         YMax = MAXVAL( Mesh % Nodes % y(NodeIndexes) )
         ZMin = MINVAL( Mesh % Nodes % z(NodeIndexes) )
         ZMax = MAXVAL( Mesh % Nodes % z(NodeIndexes) )
         ElementSize = MAX( MAX( XMax-XMin, YMax-YMin ), ZMax-ZMin )

!-------------------------------------------------------------------------------
! Is the element in one of the child quadrants?:
! Check whether element bounding box crosses any of the child quadrant
! bounding boxes:
!-------------------------------------------------------------------------------
         DO i=1, 2**dim ! loop over child quadrants

            BBox = ChildQuadrant(i) % Quadrant % BoundingBox

            eps3 = 0.0d0
            eps3 = MAX( eps3, BBox(4) - BBox(1) )
            eps3 = MAX( eps3, BBox(5) - BBox(2) )
            eps3 = MAX( eps3, BBox(6) - BBox(3) )
            eps3 = eps2 * eps3

            BBox(1:3) = BBox(1:3) - eps3
            BBox(4:6) = BBox(4:6) + eps3

            ElementInQuadrant = .TRUE.
            IF ( XMax < BBox(1) .OR. XMin > BBox(4) .OR. &
                 YMax < BBox(2) .OR. YMin > BBox(5) .OR. &
                 ZMax < BBox(3) .OR. ZMin > BBox(6) ) ElementInQuadrant = .FALSE.

!-------------------------------------------------------------------------------

            IF ( ElementInQuadrant ) THEN
               ChildQuadrant(i) % Quadrant % NElemsInQuadrant = &
                   ChildQuadrant(i) % Quadrant % NElemsInQuadrant + 1

               ChildQuadrant(i) % Quadrant % MinElementSize = &
                 MIN(ElementSize, ChildQuadrant(i) % Quadrant % MinElementSize)

               ! We allocate and store also the midlevels temporarily
               ! (for the duration of the construction routine):
               ! ----------------------------------------------------
               ElementList(i,ChildQuadrant(i) % Quadrant % NElemsInQuadrant) = &
                               MotherQuadrant % Elements(t) 
            END IF
!-------------------------------------------------------------------------------
         END DO
!-------------------------------------------------------------------------------
      END DO

!-------------------------------------------------------------------------------

      DO i=1,2**dim
         IF ( ChildQuadrant(i) % Quadrant % NElemsInQuadrant /= 0 ) THEN
            ALLOCATE ( ChildQuadrant(i) % Quadrant % Elements ( &
               ChildQuadrant(i) % Quadrant % NElemsInQuadrant ) )

            ChildQuadrant(i) % Quadrant % Elements (1: &
                ChildQuadrant(i) % Quadrant % NElemsInQuadrant) = &
                ElementList(i,1:ChildQuadrant(i) % Quadrant % NElemsInQuadrant)
         END IF
      END DO

!-------------------------------------------------------------------------------
    END SUBROUTINE PutElementsInChildQuadrants
!-------------------------------------------------------------------------------
  END SUBROUTINE BuildQuadrantTree
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Returns the local coordinate values from the given mesh
!> structure for the given set of nodal Indexes.
!-------------------------------------------------------------------------------
  SUBROUTINE CopyElementNodesFromMesh(ElementNodes, Mesh, n, Indexes)
!-------------------------------------------------------------------------------    
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n
    INTEGER, POINTER :: Indexes(:)
!-------------------------------------------------------------------------------
    INTEGER :: m
!-------------------------------------------------------------------------------    
    IF ( .NOT. ASSOCIATED( ElementNodes % x ) ) THEN
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    ELSE
      m = SIZE(ElementNodes % x)
      IF ( m < n ) THEN
        DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z)
        ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
      ELSE IF( m > n ) THEN
        ElementNodes % x(n+1:m) = 0.0_dp
        ElementNodes % y(n+1:m) = 0.0_dp
        ElementNodes % z(n+1:m) = 0.0_dp
      END IF
    END IF

    ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
    ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
    ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))
!-------------------------------------------------------------------------------
  END SUBROUTINE CopyElementNodesFromMesh
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Create a matrix representation of the Nedelec interpolation operator which 
!> operates on a vector field expressed in terms of the nodal basis functions
!> and gives the values of DOFs for obtaining its vector element (Nedelec)
!> interpolant. The current implementation assumes that all DOFs are associated
!> with edges, so that the geometric domain of the finite element given as input
!> is supposed to be one-dimensional.
!------------------------------------------------------------------------------
  SUBROUTINE NodalToNedelecPiMatrix(PiMat, Edge, Mesh, dim, SecondFamily)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(OUT) :: PiMat(2,6)      !< The interpolation operator as a matrix 
    TYPE(Element_t), POINTER, INTENT(IN) :: Edge  !< The element for which the operator is created
    TYPE(Mesh_t), POINTER, INTENT(IN) :: Mesh     !< The Edge should belong to the mesh given
    INTEGER, INTENT(IN) :: dim                    !< The number of components of the vector field  
    LOGICAL, OPTIONAL, INTENT(IN) :: SecondFamily !< To select the Nedelec family    
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: SecondKindBasis, stat
    INTEGER, ALLOCATABLE, SAVE :: Ind(:)
    
    INTEGER :: EDOFs, i, k, p, i1, i2, j1, j2, n
    REAL(KIND=dp) :: Basis(2), detJ, s, e(3), t(3), fun(3), u, v
!------------------------------------------------------------------------------
    IF ((Edge % Type % ElementCode / 100) /= 2) THEN
      CALL Warn('NodalToNedelecPiMatrix', 'A 1-dimensional element expected')
      RETURN
    END IF
        
    IF (.NOT. ASSOCIATED(Mesh % Edges)) THEN
      CALL Fatal('NodalToNedelecPiMatrix', 'Mesh edges are not associated!')
    END IF

    n = Edge % Type % NumberOfNodes
    IF (n /= 2) CALL Fatal('NodalToNedelecPiMatrix', &
        'A 2-node element expected ')

    IF (PRESENT(SecondFamily)) THEN
      SecondKindBasis = SecondFamily
    ELSE
      SecondKindBasis = .FALSE.
    END IF

    IF (SecondKindBasis) THEN
      EDOFs = 2
    ELSE
      EDOFs = 1  
    END IF

    CALL CopyElementNodesFromMesh(Nodes, Mesh, n,  Edge % NodeIndexes)

    t(1) = Nodes % x(2) - Nodes % x(1)
    t(2) = Nodes % y(2) - Nodes % y(1)
    t(3) = Nodes % z(2) - Nodes % z(1)
      
    i1 = Edge % NodeIndexes(1)
    i2 = Edge % NodeIndexes(2)
    IF (ParEnv % PEs > 1) THEN                            
      j1 = Mesh % ParallelInfo % GlobalDOFs(i1)             
      j2 = Mesh % ParallelInfo % GlobalDOFs(i2)             
    ELSE
      j1 = i1
      j2 = i2
    END IF

    IF (j2 < j1) t = -t      
    t = t/SQRT(SUM(t**2))

    PiMat = 0.0_dp
    IP = GaussPoints(Edge)
    DO p=1,IP % n
      stat = ElementInfo(Edge, Nodes, IP % u(p), IP % v(p), IP % w(p), DetJ, Basis)
      s = IP % s(p) * DetJ        

      DO k=1,dim
        e(:) = 0.0_dp
        e(k) = 1.0_dp
        DO i=1,n
          fun(:) = Basis(i) * e(:)
          IF (SecondKindBasis) THEN
            u = IP % u(p)
            v = 0.5d0*(1.0d0-sqrt(3.0d0)*u)
            PiMat(1,3*(i-1)+k) = PiMat(1,3*(i-1)+k) + s * SUM(fun*t)*v
            v = 0.5d0*(1.0d0+sqrt(3.0d0)*u)
            PiMat(2,3*(i-1)+k) = PiMat(2,3*(i-1)+k) + s * SUM(fun*t)*v
          ELSE
            PiMat(1,3*(i-1)+k) = PiMat(1,3*(i-1)+k) + s * SUM(fun*t)  
          END IF
        END DO
      END DO
    END DO    
!------------------------------------------------------------------------------
  END SUBROUTINE NodalToNedelecPiMatrix
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
END MODULE Interpolation
!-------------------------------------------------------------------------------

!> \} ElmerLib
