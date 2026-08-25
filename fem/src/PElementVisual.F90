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
!> \ingroup ElmerLib 
!> \{
!
!------------------------------------------------------------------------------
!> High-order visualization of p-element fields: turn a field given in the
!> hierarchic p-basis into an elementwise Lagrange interpolant.
!>
!> This lives in a module of its own so that its caller (p2LagrangeSwapper in
!> SolverBasics) gets a compiler-checked interface. It used to be an external
!> procedure in InterpolateMeshToMesh.F90, with a hand-written INTERFACE block
!> at the call site, and when the dummy argument PElement lost its POINTER
!> attribute the two definitions drifted apart unnoticed, leaving the whole
!> feature silently dead. Unlike the rest of that file, this routine needs
!> nothing above MeshBasics, so nothing forces it to be external.
!------------------------------------------------------------------------------
MODULE PElementVisual

  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Given a finite element field expressed in terms of the hierarchic p-basis
!> find its elementwise Lagrange interpolant of a specified degree. At the moment
!> this subroutine has limited functionality, since the coordinates of the Lagrange
!> nodes are not yet defined for all element types in the case of high p.
!------------------------------------------------------------------------------
  SUBROUTINE HierarchicPToLagrange(PElement, Degree, PSol, LSol, DOFs, PSolver)
!------------------------------------------------------------------------------
    USE MeshBasics, ONLY: AllocateElement
    USE ElementDescription
    IMPLICIT NONE

    TYPE(Element_t), TARGET, INTENT(IN) :: PElement          !< An element structure capable for p-approximation
    INTEGER, INTENT(IN) :: Degree                            !< The desired order of the Lagrange interpolant
    REAL(KIND=dp), INTENT(IN) :: PSol(:,:)                   !< The DOFs of p-solution
    REAL(KIND=dp), INTENT(INOUT) :: LSol(:,:)                !< The DOFs for the Lagrange interpolation
    INTEGER, OPTIONAL, INTENT(IN) :: DOFs                    !< The first dimension of PSol
    TYPE(Solver_t), POINTER, OPTIONAL, INTENT(IN) :: PSolver !< For seeking p-definitions from this solver
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MAX_LINE_DEGREE = 8
    INTEGER, PARAMETER :: MAX_TRIANGLE_DEGREE = 8
    INTEGER, PARAMETER :: MAX_TRIANGLE_NODES = 45
    INTEGER, PARAMETER :: MAX_QUAD_DEGREE = 8
    INTEGER, PARAMETER :: MAX_TETRA_DEGREE = 4
    INTEGER, PARAMETER :: MAX_TETRA_NODES = 35
    INTEGER, PARAMETER :: MAX_PYRAMID_DEGREE = 2
    INTEGER, PARAMETER :: MAX_PYRAMID_NODES = 15
    INTEGER, PARAMETER :: MAX_PRISM_DEGREE = 4
    INTEGER, PARAMETER :: MAX_PRISM_NODES = 75
    INTEGER, PARAMETER :: MAX_BRICK_DEGREE = 8
    INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729

    INTEGER, PARAMETER :: MAX_DEGREE = 8  ! The maximal polynomial degree over any element shape 

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t) :: PNodes
    TYPE(Element_t), POINTER :: LElement => NULL()
    LOGICAL ::  AllocationsDone = .FALSE.
    LOGICAL :: PVersion, stat
    INTEGER, ALLOCATABLE :: BasisDegree(:)
    INTEGER :: Family, Fields, MaxFields, i, j, k, n, t, nd
    INTEGER :: i_start, i_end
    INTEGER :: Offsets(MAX_DEGREE)
    INTEGER :: Edge(2), Face(4)
    REAL(KIND=dp), ALLOCATABLE :: PBasis(:)
    REAL(KIND=dp), POINTER :: NodesU(:), NodesV(:), NodesW(:)
    REAL(KIND=dp) :: ut, vt, wt, detJ
    REAL(KIND=dp) :: r(3), delta

    REAL(KIND=dp), TARGET, SAVE :: VTKLineU(MAX_LINE_DEGREE+1,MAX_LINE_DEGREE+1)

    REAL(KIND=dp), TARGET, SAVE :: VTKTriangleU(MAX_TRIANGLE_DEGREE,MAX_TRIANGLE_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKTriangleV(MAX_TRIANGLE_DEGREE,MAX_TRIANGLE_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKTriangleW(MAX_TRIANGLE_NODES)
    INTEGER, SAVE :: TriangleNodeCounts(MAX_TRIANGLE_DEGREE)

    REAL(KIND=dp), TARGET, SAVE :: VTKQuadU(MAX_QUAD_DEGREE,(MAX_QUAD_DEGREE+1)**2)
    REAL(KIND=dp), TARGET, SAVE :: VTKQuadV(MAX_QUAD_DEGREE,(MAX_QUAD_DEGREE+1)**2)
    REAL(KIND=dp), TARGET, SAVE :: VTKQuadW((MAX_QUAD_DEGREE+1)**2)
    INTEGER, SAVE :: QuadNodeCounts(MAX_QUAD_DEGREE)

    REAL(KIND=dp), TARGET, SAVE :: VTKTetraU(MAX_TETRA_DEGREE,MAX_TETRA_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKTetraV(MAX_TETRA_DEGREE,MAX_TETRA_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKTetraW(MAX_TETRA_DEGREE,MAX_TETRA_NODES)
    INTEGER, SAVE :: TetraNodeCounts(MAX_TETRA_DEGREE)
    INTEGER :: VTKTetraFaceMap(4,3)

    REAL(KIND=dp), TARGET, SAVE :: VTKPyramidU(MAX_PYRAMID_DEGREE,MAX_PYRAMID_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKPyramidV(MAX_PYRAMID_DEGREE,MAX_PYRAMID_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKPyramidW(MAX_PYRAMID_DEGREE,MAX_PYRAMID_NODES)
    INTEGER, SAVE :: PyramidNodeCounts(MAX_PYRAMID_DEGREE)

    REAL(KIND=dp), TARGET, SAVE :: VTKPrismU(MAX_PRISM_DEGREE,MAX_PRISM_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKPrismV(MAX_PRISM_DEGREE,MAX_PRISM_NODES)
    REAL(KIND=dp), TARGET, SAVE :: VTKPrismW(MAX_PRISM_DEGREE,MAX_PRISM_NODES)
    INTEGER, SAVE :: PrismNodeCounts(MAX_PRISM_DEGREE)

    REAL(KIND=dp), TARGET, SAVE :: VTKBrickU(MAX_BRICK_DEGREE,(MAX_BRICK_DEGREE+1)**3)
    REAL(KIND=dp), TARGET, SAVE :: VTKBrickV(MAX_BRICK_DEGREE,(MAX_BRICK_DEGREE+1)**3)
    REAL(KIND=dp), TARGET, SAVE :: VTKBrickW(MAX_BRICK_DEGREE,(MAX_BRICK_DEGREE+1)**3)
    INTEGER, SAVE :: BrickNodeCounts(MAX_BRICK_DEGREE)
    INTEGER :: VTKBrickFaceMap(6,4)

    CHARACTER(*), PARAMETER :: Caller = 'HierarchicPToLagrange'

    SAVE PNodes, LElement, AllocationsDone, PBasis, BasisDegree
!------------------------------------------------------------------------------
    LSol = 0.0d0

    MaxFields = MIN(SIZE(PSol,1), SIZE(LSol,1))
    IF (PRESENT(DOFs)) THEN
      IF (DOFs > MaxFields) THEN
        CALL Warn(Caller, 'Too small arrays to write all fields')
        Fields = MaxFields
      ELSE
        Fields = DOFs
      END IF
    ELSE
      Fields = MaxFields ! Write as many fields as seems to be possible
    END IF

    Family = PElement % TYPE % ElementCode / 100

    IF (Degree == 1) THEN
      DO k=1,Fields
        LSol(k,1:Family) = PSol(k,1:Family)
      END DO
      RETURN
    END IF

    Mesh => CurrentModel % Solver % Mesh
    IF (PRESENT(PSolver)) THEN
      IF(ASSOCIATED(PSolver)) Mesh => PSolver % Mesh
    END IF

    IF (.NOT. AllocationsDone) THEN
      LElement => AllocateElement()
      n = Mesh % MaxElementDOFs
      ALLOCATE(PBasis(n), BasisDegree(n), PNodes % x(n), &
          PNodes % y(n), PNodes % z(n))

      ! --------------------------------------------------------------------------
      ! Create the nodes for line elements (works for any given degree):
      ! --------------------------------------------------------------------------
      VTKLineU(1:MAX_LINE_DEGREE,1) = -1.0d0
      VTKLineU(1:MAX_LINE_DEGREE,2) = 1.0d0
      DO k=2,MAX_LINE_DEGREE
        delta = 2.0d0/k
        ut = -1.0d0
        DO i=1,k-1
          ut = ut + delta 
          VTKLineU(k,i+2) = ut
        END DO
      END DO
      VTKLineU(MAX_LINE_DEGREE+1,:) = 0.0d0  ! One extra row to get zero values for non-active parameters

      ! --------------------------------------------------------------------------
      ! Create the nodes for triangles (basically works for any given degree):
      ! --------------------------------------------------------------------------
      VTKTriangleU(1:MAX_TRIANGLE_DEGREE,1) = -1.0d0
      VTKTriangleU(1:MAX_TRIANGLE_DEGREE,2) = 1.0d0
      VTKTriangleU(1:MAX_TRIANGLE_DEGREE,3) = 0.0d0
      VTKTriangleV(1:MAX_TRIANGLE_DEGREE,1) = 0.0d0
      VTKTriangleV(1:MAX_TRIANGLE_DEGREE,2) = 0.0d0
      VTKTriangleV(1:MAX_TRIANGLE_DEGREE,3) = SQRT(3.0d0)
      VTKTriangleW = 0.0d0

      ! The edges:
      TriangleNodeCounts = 3
      DO k=1,3
        Edge = GetTriangleEdgeMap(k)
        Offsets(1:MAX_TRIANGLE_DEGREE) = TriangleNodeCounts(1:MAX_TRIANGLE_DEGREE)
        CALL NodesOnEdge(VTKTriangleU, VTKTriangleV, Edge(1), Edge(2), MAX_TRIANGLE_DEGREE, &
            Offsets, TriangleNodeCounts)
      END DO

      ! The bubbles:
      VTKTriangleU(3,10) = 0.0d0
      VTKTriangleV(3,10) = SQRT(3.0d0)/3.0d0
      TriangleNodeCounts(3) = TriangleNodeCounts(3) + 1

      DO j=4,MAX_TRIANGLE_DEGREE
        !
        ! We can reuse the nodes of (j-3)th order triangle:
        !
        i_start = 3*j + 1
        SELECT CASE(j)
        CASE(4)
          n = 3
        CASE(5)
          n = 6
        CASE(6)
          n = 10
        CASE(7)
          n = 15
        CASE(8)
          n = 21
        CASE DEFAULT
          CALL Fatal(Caller, 'Triangle interior nodes supported up to degree 8')
        END SELECT
        i_end = 3*j + n

        VTKTriangleU(j,i_start:i_end) = (j-3.0d0)/j * VTKTriangleU(j-3,1:n)
        vt = SQRT(3.0d0)/j
        VTKTriangleV(j,i_start:i_end) = vt + (j-3.0d0)/j * VTKTriangleV(j-3,1:n)
        TriangleNodeCounts(j) = TriangleNodeCounts(j) + n
      END DO

      ! --------------------------------------------------------------------------
      ! Create the nodes for quads (works for any given degree):
      ! --------------------------------------------------------------------------
      VTKQuadU(1:MAX_QUAD_DEGREE,1) = -1.0d0
      VTKQuadU(1:MAX_QUAD_DEGREE,2:3) = 1.0d0
      VTKQuadU(1:MAX_QUAD_DEGREE,4) = -1.0d0
      VTKQuadV(1:MAX_QUAD_DEGREE,1:2) = -1.0d0
      VTKQuadV(1:MAX_QUAD_DEGREE,3:4) = 1.0d0

      ! The edges:
      QuadNodeCounts = 4
      DO k=1,4
        Edge = GetQuadEdgeMap(k)
        Offsets(1:MAX_QUAD_DEGREE) = QuadNodeCounts(1:MAX_QUAD_DEGREE)
        CALL NodesOnEdge(VTKQuadU, VTKQuadV, Edge(1), Edge(2), MAX_QUAD_DEGREE, &
            Offsets, QuadNodeCounts)
      END DO

      ! The bubbles:
      Offsets(1:MAX_QUAD_DEGREE) = QuadNodeCounts(1:MAX_QUAD_DEGREE)
      DO k=2,MAX_QUAD_DEGREE
        delta = 2.0d0/k
        vt = -1.0d0
        DO i=1,k-1
          vt = vt + delta
          ut = -1.0d0
          DO j=1,k-1
            ut = ut + delta 
            VTKQuadU(k,j+Offsets(k)) = ut
            VTKQuadV(k,j+Offsets(k)) = vt         
            QuadNodeCounts(k) = QuadNodeCounts(k) + 1
          END DO
          Offsets(k) = QuadNodeCounts(k)
        END DO
      END DO
      VTKQuadW = 0.0d0

      ! --------------------------------------------------------------------------
      ! Create the nodes for tetrahedra:
      ! --------------------------------------------------------------------------
      LElement % Type => GetElementType(504, .FALSE.)
      CALL GetRefPElementNodes(LElement % Type, VTKTetraU(1,1:4), &
          VTKTetraV(1,1:4), VTKTetraW(1,1:4)) 

      DO k=2,MAX_TETRA_DEGREE
        VTKTetraU(k,1:4) = VTKTetraU(1,1:4)
        VTKTetraV(k,1:4) = VTKTetraV(1,1:4)
        VTKTetraW(k,1:4) = VTKTetraW(1,1:4)
      END DO

      ! The edges:
      TetraNodeCounts = 4
      DO k=1,6
        Edge = GetTetraEdgeMap(k)
        Offsets(1:MAX_TETRA_DEGREE) = TetraNodeCounts(1:MAX_TETRA_DEGREE)
        IF (k == 3) THEN
          CALL NodesOnEdge(VTKTetraU, VTKTetraV, Edge(2), Edge(1), MAX_TETRA_DEGREE, &
              Offsets, TetraNodeCounts, VTKTetraW)
        ELSE
          CALL NodesOnEdge(VTKTetraU, VTKTetraV, Edge(1), Edge(2), MAX_TETRA_DEGREE, &
              Offsets, TetraNodeCounts, VTKTetraW)
        END IF
      END DO
      
      ! The faces:
      Offsets(1:MAX_TETRA_DEGREE) = TetraNodeCounts(1:MAX_TETRA_DEGREE)
      VTKTetraFaceMap(1,:) = (/ 1,2,4 /)
      VTKTetraFaceMap(2,:) = (/ 3,4,2 /)
      VTKTetraFaceMap(3,:) = (/ 1,4,3 /)
      VTKTetraFaceMap(4,:) = (/ 1,3,2 /)

      DO j=3,MAX_TETRA_DEGREE
        SELECT CASE(j)
        CASE(3)
          i_start = 10
          i_end = 10
        CASE(4)
          i_start = 13
          i_end = 15
        CASE DEFAULT
          CALL Fatal(Caller, 'Tetra face nodes supported up to degree 4')
        END SELECT
        n = i_end - i_start + 1 ! The number of new nodes per face

        DO k=1,4
          DO i=1,n
            ! Pick the positions of face nodes from a lower dimensional element 
            r(1) = VTKTriangleU(j,i+i_start-1)
            r(2) = VTKTriangleV(j,i+i_start-1)
            ut = 0.0d0
            vt = 0.0d0
            wt = 0.0d0
            DO t=1,3
              PBasis(1) = TriangleNodalPBasis(t, r(1), r(2))
              ut = ut + PBasis(1) * VTKTetraU(1,VTKTetraFaceMap(k,t))
              vt = vt + PBasis(1) * VTKTetraV(1,VTKTetraFaceMap(k,t))
              wt = wt + PBasis(1) * VTKTetraW(1,VTKTetraFaceMap(k,t))
            END DO
            VTKTetraU(j,i+Offsets(j)) = ut
            VTKTetraV(j,i+Offsets(j)) = vt
            VTKTetraW(j,i+Offsets(j)) = wt
            TetraNodeCounts(j) = TetraNodeCounts(j) + 1
          END DO
          Offsets(j) = TetraNodeCounts(j)
        END DO
      END DO

      ! Bubbles:
      !DO j=4,MAX_TETRA_DEGREE
      r(1) = SUM(VTKTetraU(1,1:4))/4
      r(2) = SUM(VTKTetraV(1,1:4))/4
      r(3) = SUM(VTKTetraW(1,1:4))/4
      VTKTetraU(4,1+Offsets(4)) = r(1)
      VTKTetraV(4,1+Offsets(4)) = r(2)
      VTKTetraW(4,1+Offsets(4)) = r(3)
      TetraNodeCounts(4) = TetraNodeCounts(4) + 1
      !END DO

      ! --------------------------------------------------------------------------
      ! Create the nodes for pyramids:
      ! --------------------------------------------------------------------------
      LElement % Type => GetElementType(605, .FALSE.)
      CALL GetRefPElementNodes(LElement % Type, VTKPyramidU(1,1:5), &
          VTKPyramidV(1,1:5), VTKPyramidW(1,1:5)) 

      DO k=2,MAX_PYRAMID_DEGREE
        VTKPyramidU(k,1:5) = VTKPyramidU(1,1:5)
        VTKPyramidV(k,1:5) = VTKPyramidV(1,1:5)
        VTKPyramidW(k,1:5) = VTKPyramidW(1,1:5)
      END DO

      ! The edges:
      PyramidNodeCounts = 5
      DO k=1,8
        Edge = GetPyramidEdgeMap(k)
        Offsets(1:MAX_PYRAMID_DEGREE) = PyramidNodeCounts(1:MAX_PYRAMID_DEGREE)
        CALL NodesOnEdge(VTKPyramidU, VTKPyramidV, Edge(1), Edge(2), MAX_PYRAMID_DEGREE, &
            Offsets, PyramidNodeCounts, VTKPyramidW)
      END DO
      Offsets(1:MAX_PYRAMID_DEGREE) = PyramidNodeCounts(1:MAX_PYRAMID_DEGREE)

      ! The quad face:
      Face = GetPyramidFaceMap(1)
      DO j=2,MAX_PYRAMID_DEGREE
        SELECT CASE(j)
        CASE(2)
          i_start = 9
          i_end = 9
        CASE(3)
          i_start = 13
          i_end = 16
        CASE DEFAULT
          CALL Fatal(Caller, 'Pyramid (quad) face nodes supported up to degree 3')
        END SELECT
        n = i_end - i_start + 1 ! The number of new nodes per face

        DO i=1,n
          ! Pick the positions of face nodes from a lower dimensional element 
          r(1) = VTKQuadU(j,i+i_start-1)
          r(2) = VTKQuadV(j,i+i_start-1)
          ut = 0.0d0
          vt = 0.0d0
          wt = 0.0d0
          DO t=1,4
            PBasis(1) = QuadNodalPBasis(t, r(1), r(2))
            ut = ut + PBasis(1) * VTKPyramidU(1,Face(t))
            vt = vt + PBasis(1) * VTKPyramidV(1,Face(t))
            wt = wt + PBasis(1) * VTKPyramidW(1,Face(t))
          END DO
          VTKPyramidU(j,i+Offsets(j)) = ut
          VTKPyramidV(j,i+Offsets(j)) = vt
          VTKPyramidW(j,i+Offsets(j)) = wt
          PyramidNodeCounts(j) = PyramidNodeCounts(j) + 1
        END DO
        Offsets(j) = PyramidNodeCounts(j) 
      END DO

      ! The triangular faces (TO DO: support for degrees p > 3):

      ! The bubbles (TO DO: support for degrees p > 2):
      r(1) = SUM(VTKPyramidU(1,1:5))/5
      r(2) = SUM(VTKPyramidV(1,1:5))/5
      r(3) = SUM(VTKPyramidW(1,1:5))/5
      VTKPyramidU(2,1+Offsets(2)) = r(1)
      VTKPyramidV(2,1+Offsets(2)) = r(2)
      VTKPyramidW(2,1+Offsets(2)) = r(3)
      PyramidNodeCounts(2) = PyramidNodeCounts(2) + 1

      ! --------------------------------------------------------------------------
      ! Create the nodes for prisms:
      ! --------------------------------------------------------------------------
      LElement % Type => GetElementType(706, .FALSE.)
      CALL GetRefPElementNodes(LElement % Type, VTKPrismU(1,1:6), &
          VTKPrismV(1,1:6), VTKPrismW(1,1:6)) 

      DO k=2,MAX_PRISM_DEGREE
        VTKPrismU(k,1:6) = VTKPrismU(1,1:6)
        VTKPrismV(k,1:6) = VTKPrismV(1,1:6)
        VTKPrismW(k,1:6) = VTKPrismW(1,1:6)
      END DO

      ! The edges:
      PrismNodeCounts = 6
      DO k=1,9
        Edge = GetWedgeEdgeMap(k)
        Offsets(1:MAX_PRISM_DEGREE) = PrismNodeCounts(1:MAX_PRISM_DEGREE)
        CALL NodesOnEdge(VTKPrismU, VTKPrismV, Edge(1), Edge(2), MAX_PRISM_DEGREE, &
            Offsets, PrismNodeCounts, VTKPrismW)
      END DO
      Offsets(1:MAX_PRISM_DEGREE) = PrismNodeCounts(1:MAX_PRISM_DEGREE)

      ! The triangular faces:
      DO j=3,MAX_PRISM_DEGREE
        SELECT CASE(j)
        CASE(3)
          i_start = 10
          i_end = 10
        CASE(4)
          i_start = 13
          i_end = 15
        CASE DEFAULT
          CALL Fatal(Caller, 'Prism (triangular) face nodes supported up to degree 4')
        END SELECT
        n = i_end - i_start + 1 ! The number of new nodes per face
        DO k=1,2
          Face = GetWedgeFaceMap(k)
          DO i=1,n
            ! Pick the positions of face nodes from a lower dimensional element 
            r(1) = VTKTriangleU(j,i+i_start-1)
            r(2) = VTKTriangleV(j,i+i_start-1)
            ut = 0.0d0
            vt = 0.0d0
            wt = 0.0d0
            DO t=1,3
              PBasis(1) = TriangleNodalPBasis(t, r(1), r(2))
              ut = ut + PBasis(1) * VTKPrismU(1,Face(t))
              vt = vt + PBasis(1) * VTKPrismV(1,Face(t))
              wt = wt + PBasis(1) * VTKPrismW(1,Face(t))
            END DO
            VTKPrismU(j,i+Offsets(j)) = ut
            VTKPrismV(j,i+Offsets(j)) = vt
            VTKPrismW(j,i+Offsets(j)) = wt
            PrismNodeCounts(j) = PrismNodeCounts(j) + 1
          END DO
          Offsets(j) = PrismNodeCounts(j)
        END DO
      END DO

      ! The quad faces:
      DO j=2,MAX_PRISM_DEGREE
        i_start = 4*j + 1
        i_end = (j+1)**2
        n = (j-1)**2  ! The number of new nodes per face n = i_end - i_start + 1

        DO k=3,5
          Face = GetWedgeFaceMap(k)
          DO i=1,n
            ! Pick the positions of face nodes from a lower dimensional element 
            r(1) = VTKQuadU(j,i+i_start-1)
            r(2) = VTKQuadV(j,i+i_start-1)
            ut = 0.0d0
            vt = 0.0d0
            wt = 0.0d0
            DO t=1,4
              PBasis(1) = QuadNodalPBasis(t, r(1), r(2))
              ut = ut + PBasis(1) * VTKPrismU(1,Face(t))
              vt = vt + PBasis(1) * VTKPrismV(1,Face(t))
              wt = wt + PBasis(1) * VTKPrismW(1,Face(t))
            END DO
            VTKPrismU(j,i+Offsets(j)) = ut
            VTKPrismV(j,i+Offsets(j)) = vt
            VTKPrismW(j,i+Offsets(j)) = wt
            PrismNodeCounts(j) = PrismNodeCounts(j) + 1
          END DO
          Offsets(j) = PrismNodeCounts(j) 
        END DO
      END DO

      DO j=3,MAX_PRISM_DEGREE
        i_start = 3*j + 1
        n = (j-1)*(j-2)/2
        i_end = 3*j + n
        DO i=1,j-1
          VTKPrismU(j,Offsets(j)+1:Offsets(j)+n) = VTKTriangleU(j,i_start:i_end)
          VTKPrismV(j,Offsets(j)+1:Offsets(j)+n) = VTKTriangleV(j,i_start:i_end)
          VTKPrismW(j,Offsets(j)+1:Offsets(j)+n) = VTKLineU(j,2+i)
          PrismNodeCounts(j) = PrismNodeCounts(j) + n
          Offsets(j) = PrismNodeCounts(j)
        END DO
      END DO

      ! --------------------------------------------------------------------------
      ! Create the nodes for hexahedra:
      ! --------------------------------------------------------------------------
      LElement % Type => GetElementType(808, .FALSE.)
      CALL GetRefPElementNodes(LElement % Type, VTKBrickU(1,1:8), &
          VTKBrickV(1,1:8), VTKBrickW(1,1:8)) 

      DO k=2,MAX_BRICK_DEGREE
        VTKBrickU(k,1:8) = VTKBrickU(1,1:8)
        VTKBrickV(k,1:8) = VTKBrickV(1,1:8)
        VTKBrickW(k,1:8) = VTKBrickW(1,1:8)
      END DO
      BrickNodeCounts = 8

      ! The edges:
      ! It seems that VTK cell types 72 and 12/29 are not interchangeable.
      ! We first pick the edge map which works for 72 and then re-run
      ! the same code up to the degree 2 with a different edge map to obtain consistency.
      DO k=1,12
        ! The following is needed for 72:
        SELECT CASE(k)
        CASE(11)
          Edge = GetBrickEdgeMap(12)
        CASE(12)
          Edge = GetBrickEdgeMap(11)
        CASE DEFAULT
          Edge = GetBrickEdgeMap(k)
        END SELECT
        Offsets(1:MAX_BRICK_DEGREE) = BrickNodeCounts(1:MAX_BRICK_DEGREE)
        CALL NodesOnEdge(VTKBrickU, VTKBrickV, Edge(1), Edge(2), MAX_BRICK_DEGREE, &
            Offsets, BrickNodeCounts, VTKBrickW)
      END DO
      ! The re-run needs a reset:
      BrickNodeCounts(2) = 8
      DO k=1,12
        ! The following is needed for 29:
        Edge = GetBrickEdgeMap(k)
        Offsets(2) = BrickNodeCounts(2)
        CALL NodesOnEdge(VTKBrickU, VTKBrickV, Edge(1), Edge(2), 2, &
            Offsets, BrickNodeCounts, VTKBrickW)
      END DO

      ! The faces:
      Offsets(1:MAX_BRICK_DEGREE) = BrickNodeCounts(1:MAX_BRICK_DEGREE)
      VTKBrickFaceMap(1,:) = (/ 1,4,8,5 /)
      VTKBrickFaceMap(2,:) = (/ 2,3,7,6 /)
      VTKBrickFaceMap(3,:) = (/ 1,2,6,5 /)
      VTKBrickFaceMap(4,:) = (/ 4,3,7,8 /)
      VTKBrickFaceMap(5,:) = (/ 1,2,3,4 /)
      VTKBrickFaceMap(6,:) = (/ 5,6,7,8 /)

      DO j=2,MAX_BRICK_DEGREE
        i_start = 4*j + 1
        i_end = (j+1)**2
        n = (j-1)**2  ! The number of new nodes per face n = i_end - i_start + 1

        DO k=1,6
          DO i=1,n
            ! Pick the positions of face nodes from a lower dimensional element 
            r(1) = VTKQuadU(j,i+i_start-1)
            r(2) = VTKQuadV(j,i+i_start-1)
            ut = 0.0d0
            vt = 0.0d0
            wt = 0.0d0
            DO t=1,4
              PBasis(1) = QuadNodalPBasis(t, r(1), r(2))
              ut = ut + PBasis(1) * VTKBrickU(1,VTKBrickFaceMap(k,t))
              vt = vt + PBasis(1) * VTKBrickV(1,VTKBrickFaceMap(k,t))
              wt = wt + PBasis(1) * VTKBrickW(1,VTKBrickFaceMap(k,t))
            END DO
            VTKBrickU(j,i+Offsets(j)) = ut
            VTKBrickV(j,i+Offsets(j)) = vt
            VTKBrickW(j,i+Offsets(j)) = wt
            BrickNodeCounts(j) = BrickNodeCounts(j) + 1
          END DO
          Offsets(j) = BrickNodeCounts(j)
        END DO
      END DO

      ! Bubbles:
      VTKBrickU(2,1+Offsets(2)) = 0.0d0
      VTKBrickV(2,1+Offsets(2)) = 0.0d0
      VTKBrickW(2,1+Offsets(2)) = 0.0d0
      BrickNodeCounts(2) = BrickNodeCounts(2) + 1

      DO j=3,MAX_BRICK_DEGREE
        i_start = 4*j + 1
        i_end = (j+1)**2
        n = (j-1)**2
        DO i=1,j-1
          VTKBrickU(j,Offsets(j)+1:Offsets(j)+n) = VTKQuadU(j,i_start:i_end)
          VTKBrickV(j,Offsets(j)+1:Offsets(j)+n) = VTKQuadV(j,i_start:i_end)
          VTKBrickW(j,Offsets(j)+1:Offsets(j)+n) = VTKLineU(j,2+i)
          BrickNodeCounts(j) = BrickNodeCounts(j) + n
          Offsets(j) = BrickNodeCounts(j)
        END DO
      END DO
      AllocationsDone = .TRUE.
    END IF

    IF (PRESENT(PSolver)) THEN
      PVersion = IsActivePElement(PElement, PSolver)
    ELSE
      PVersion = IsPElement(PElement)
    END IF
    
    IF (.NOT. PVersion) THEN
      CALL Warn(Caller, 'The input element is not a p-element, returning')
      RETURN
    END IF
    
    PNodes % x(:) = 0.0d0
    PNodes % y(:) = 0.0d0
    PNodes % z(:) = 0.0d0
    !
    ! NOTE: We shall not need the derivatives of basis functions or
    ! the determinant of the element mapping. Hence using the following
    ! nodal coordinates is sufficient, although in principle we
    ! might disregard some values 
    !
    n = PElement % TYPE % NumberOfNodes
    PNodes % x(1:n) = Mesh % Nodes % x(PElement % NodeIndexes(1:n))
    PNodes % y(1:n) = Mesh % Nodes % y(PElement % NodeIndexes(1:n))
    PNodes % z(1:n) = Mesh % Nodes % z(PElement % NodeIndexes(1:n))

    SELECT CASE(Family)
    CASE(2)
      IF (Degree > MAX_LINE_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (lines) exceeds the supported maximum')
      END IF
      n = Degree + 1
      NodesU => VTKLineU(Degree,1:n)
      NodesV => VTKLineU(MAX_LINE_DEGREE+1,1:n)
      NodesW => VTKLineU(MAX_LINE_DEGREE+1,1:n)

    CASE(3)
      IF (Degree > MAX_TRIANGLE_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (triangles) exceeds the supported maximum')
      END IF
      n = TriangleNodeCounts(Degree)
      NodesU => VTKTriangleU(Degree,1:n)
      NodesV => VTKTriangleV(Degree,1:n)
      NodesW => VTKTriangleW(1:n)

    CASE(4)
      IF (Degree > MAX_QUAD_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (quads) exceeds the supported maximum')
      END IF
      n = QuadNodeCounts(Degree)
      NodesU => VTKQuadU(Degree,1:n)
      NodesV => VTKQuadV(Degree,1:n)
      NodesW => VTKQuadW(1:n)

    CASE(5)
      IF (Degree > MAX_TETRA_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (tetrahedra) exceeds the supported maximum')
      END IF
      n = TetraNodeCounts(Degree)
      NodesU => VTKTetraU(Degree,1:n)
      NodesV => VTKTetraV(Degree,1:n)
      NodesW => VTKTetraW(Degree,1:n)

    CASE(6)
      IF (Degree > MAX_PYRAMID_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (pyramids) exceeds the supported maximum')
      END IF
      n = PyramidNodeCounts(Degree)
      NodesU => VTKPyramidU(Degree,1:n)
      NodesV => VTKPyramidV(Degree,1:n)
      NodesW => VTKPyramidW(Degree,1:n)

    CASE(7)
      IF (Degree > MAX_PRISM_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (prisms) exceeds the supported maximum')
      END IF
      n = PrismNodeCounts(Degree)
      NodesU => VTKPrismU(Degree,1:n)
      NodesV => VTKPrismV(Degree,1:n)
      NodesW => VTKPrismW(Degree,1:n)

    CASE(8)
      IF (Degree > MAX_BRICK_DEGREE) THEN
        CALL Fatal(Caller, 'Lagrange Element Degree (bricks) exceeds the supported maximum')
      END IF
      n = BrickNodeCounts(Degree)
      NodesU => VTKBrickU(Degree,1:n)
      NodesV => VTKBrickV(Degree,1:n)
      NodesW => VTKBrickW(Degree,1:n)

    CASE DEFAULT
      CALL Fatal(Caller, 'An unknown element family')
    END SELECT

    BasisDegree(:) = 0

    DO t = 1,n
      ut = NodesU(t)
      vt = NodesV(t)
      wt = NodesW(t)

      !PRINT *, 'CALLING AT POINT ', t,ut,vt,wt      

      ! It seems that the p-basis for pyramids cannot be evaluated at the fifth node,
      ! so we make a small "error":
      IF (Family == 6 .AND. t==5) wt = wt - 1.0d-8

      stat = ElementInfo(PElement, PNodes, ut, vt, wt, detJ, PBasis, BasisDegree=BasisDegree, &
          USolver=PSolver)
      IF (t == 1) THEN
        nd = COUNT(BasisDegree > 0) 
        IF (nd == 0) CALL Fatal(Caller, 'p-basis needed but the classical basis returned')
        ! PRINT *, 'p-basis dimension ', nd, SUM( PBasis(1:nd)), SUM(PBasis(1:PElement % Type % NumberOfNodes))
      END IF
        
      DO k = 1,Fields
        LSol(k,t) = SUM(PBasis(1:nd) * PSol(k,1:nd))
      END DO
    END DO

  CONTAINS

    ! --------------------------------------------------------------------------------
    ! Given the vertices of a 2-node edge in the first rows of the coordinate arrays
    ! NodesU, NodesV and NodesW, create additional nodes corresponding to the Lagrange
    ! interpolation up to the polynomial degree Maxdegree. The edge is defined by
    ! specifying its start and end indices and is assumed to be of length 2.
    ! The array Offsets gives the node counts before adding new nodes on the edge,
    ! while the array NodeCounts gives the counts after this action.
    ! --------------------------------------------------------------------------------
    SUBROUTINE NodesOnEdge(NodesU, NodesV, ind_start, ind_end, MaxDegree, Offsets, &
        NodeCounts, NodesW)

      IMPLICIT NONE
      REAL(KIND=dp), INTENT(INOUT):: NodesU(:,:), NodesV(:,:) 
      INTEGER, INTENT(IN) :: ind_start, ind_end
      INTEGER, INTENT(IN) :: MaxDegree
      INTEGER, INTENT(IN) :: OffSets(:)
      INTEGER, INTENT(INOUT) :: NodeCounts(:)
      REAL(KIND=dp), OPTIONAL, INTENT(INOUT):: NodesW(:,:)
      ! --------------------------------------------------------------
      REAL(KIND=dp), PARAMETER :: L = 2.0_dp
      LOGICAL :: CreateW
      INTEGER :: i, k
      REAL(KIND=dp) :: r(3), ut, vt, wt, delta
      ! --------------------------------------------------------------
      CreateW = PRESENT(NodesW)
      r(1) =  NodesU(1,ind_end) - NodesU(1,ind_start)
      r(2) =  NodesV(1,ind_end) - NodesV(1,ind_start)
      IF (CreateW) THEN
        r(3) =  NodesW(1,ind_end) - NodesW(1,ind_start)
        k = 3
      ELSE
        k = 2
      END IF
      r(1:k) = r(1:k)/SQRT(SUM(r(1:k)**2))

      DO k=2,MaxDegree
        delta = L/k
        ut = NodesU(1,ind_start)
        vt = NodesV(1,ind_start)
        IF (CreateW) wt = NodesW(1,ind_start)
        DO i=1,k-1
          ut = ut + delta * r(1)
          vt = vt + delta * r(2)
          NodesU(k,i+Offsets(k)) = ut
          NodesV(k,i+Offsets(k)) = vt
          IF (CreateW) THEN
            wt = wt + delta * r(3)
            NodesW(k,i+Offsets(k)) = wt
          END IF
          NodeCounts(k) = NodeCounts(k) + 1
        END DO
      END DO
    END SUBROUTINE NodesOnEdge

!------------------------------------------------------------------------------
  END SUBROUTINE HierarchicPToLagrange
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END MODULE PElementVisual
!------------------------------------------------------------------------------

!> \} ElmerLib
