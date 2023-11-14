!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Joe Todd
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 18/10/19
! *
! *****************************************************************************
!>  Computes anisotropic target element size based on distance from calving front
SUBROUTINE GlacierMeshMetricAniso(Model, nodenumber, y, TargetLength)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: y, TargetLength(:)
  !----------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Variable_t), POINTER :: TimeVar
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp) :: xx,yy,zz,t,told, this_dist, MinDist2, Dist,&
       lc_mindist, lc_maxdist, lc_min, lc_max,s,dx,dz,NodeDepth,NodeElev
  REAL(KIND=dp), DIMENSION(3) :: Point
  INTEGER, POINTER :: FrontPerm(:),BottomPerm(:),SurfacePerm(:),Nodeindexes(:)
  INTEGER :: i,NNodes, NFront, layers, nbottom, n, BCTag, nsurface
  LOGICAL :: NewTime,FirstTime=.TRUE., Debug, Found, ExtrudeLayers, ThisBC
  LOGICAL, ALLOCATABLE :: BottomElemMask(:), SurfElemMask(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FrontMaskName,BottomMaskName,SurfaceMaskName

  SAVE :: FirstTime, told, FrontPerm, Mesh, NNodes, nfront, lc_maxdist, lc_mindist,&
       lc_max, lc_min, dz, newtime, ExtrudeLayers, layers, BottomElemMask, BottomPerm,&
       SurfElemMask, SurfacePerm

  Debug = .FALSE.
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)
  Solver => Model % Solver
  FrontMaskName="Calving Front Mask"
  BottomMaskName="Bottom Surface Mask"
  SurfaceMaskName='Top Surface Mask'

  IF (FirstTime) THEN
    FirstTime = .FALSE.
    NewTime = .TRUE.
    told = t

    Mesh => Model % Mesh
    Material => GetMaterial(Mesh % Elements(1)) !TODO, this is not generalised

    lc_maxdist = ListGetConstReal(Material, "GlacierMeshMetric Max Distance",  DefValue=2000.0_dp)
    lc_mindist = ListGetConstReal(Material, "GlacierMeshMetric Min Distance",  DefValue=200.0_dp)
    lc_max = ListGetConstReal(Material, "GlacierMeshMetric Max LC",  DefValue=500.0_dp)
    lc_min = ListGetConstReal(Material, "GlacierMeshMetric Min LC",  DefValue=75.0_dp)
    ExtrudeLayers = ListGetLogical(Material,"GlacierMeshMetric Vertical From Layers", DefValue=.FALSE.)
    IF(ExtrudeLayers) THEN
      layers = ListGetInteger(Material, "GlacierMeshMetric Vertical Layers", Found=Found)
      IF(.NOT. Found) CALL FATAL('GlacierMeshMetric', 'Vertical requested by layers but number of layers not given')
    ELSE
      dz = ListGetConstReal(Material, "GlacierMeshMetric Vertical LC",  DefValue=50.0_dp)
    END IF
  END IF

  IF(t > told) NewTime = .TRUE. !TODO - replace this with Samuel's mesh counter logic
  IF(NewTime) THEN
    told = t
    NewTime = .FALSE.

    Mesh => Model % Mesh
    NNodes = Mesh % NumberOfNodes

    !mask is the most computationally expense element of this routine
    IF(ASSOCIATED(FrontPerm)) DEALLOCATE(FrontPerm)
    CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
         .FALSE., FrontPerm, nfront, ParallelComm = .FALSE.)

    IF(ExtrudeLayers) THEN

      ! create top and bottom perms
      IF(ASSOCIATED(BottomPerm)) DEALLOCATE(BottomPerm)
      CALL MakePermUsingMask( Model, Solver, Mesh, BottomMaskName, &
      .FALSE., BottomPerm, nbottom, ParallelComm = .FALSE.)
      IF(ASSOCIATED(SurfacePerm)) DEALLOCATE(SurfacePerm)
      CALL MakePermUsingMask( Model, Solver, Mesh, SurfaceMaskName, &
      .FALSE., SurfacePerm, nsurface, ParallelComm = .FALSE.)

      !create mask of elems as with an unstructred mesh all nodes can be in mask but not elem
      DO i=1,Model % NumberOfBCs
        ThisBC = ListGetLogical(Model % BCs(i) % Values, BottomMaskName, Found)
        IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
        BCtag =  Model % BCs(i) % Tag
        EXIT
      END DO

      IF(ALLOCATED(BottomElemMask)) DEALLOCATE(BottomElemMask)
      ALLOCATE(BottomElemMask(Mesh % NumberOfBulkElements &
           + Mesh % NumberOfBoundaryElements))
      BottomElemMask = .TRUE.
      DO i=Mesh % NumberOfBulkElements+1, &
           Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
         IF(Mesh % Elements(i) % BoundaryInfo % constraint == BCTag) &
              BottomElemMask(i) = .FALSE.
      END DO

      !create mask of elems as with an unstructred mesh all nodes can be in mask but not elem
      DO i=1,Model % NumberOfBCs
        ThisBC = ListGetLogical(Model % BCs(i) % Values, SurfaceMaskName, Found)
        IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
        BCtag =  Model % BCs(i) % Tag
        EXIT
      END DO

      IF(ALLOCATED(SurfElemMask)) DEALLOCATE(SurfElemMask)
      ALLOCATE(SurfElemMask(Mesh % NumberOfBulkElements &
           + Mesh % NumberOfBoundaryElements))
      SurfElemMask = .TRUE.
      DO i=Mesh % NumberOfBulkElements+1, &
           Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
         IF(Mesh % Elements(i) % BoundaryInfo % constraint == BCTag) &
              SurfElemMask(i) = .FALSE.
      END DO
    END IF

  END IF

  xx = Mesh % Nodes % x(nodenumber)
  yy = Mesh % Nodes % y(nodenumber)
  zz = Mesh % Nodes % z(nodenumber)

  MinDist2 = HUGE(MinDist2)
  DO i=1,NNodes
    IF(FrontPerm(i) == 0) CYCLE
    !Note - sqrt is an expensive FLOP so only do it once...
    this_dist = ((xx - Mesh % Nodes % x(i))**2.0) + &
         ((yy - Mesh % Nodes % y(i))**2.0) + &
         ((zz - Mesh % Nodes % z(i))**2.0)
    MinDist2 = MIN(this_dist, MinDist2)
  END DO
  Dist = SQRT(MinDist2)

  !Apply caps to this distance:
  Dist = MAX(Dist,lc_mindist)
  Dist = MIN(Dist,lc_maxdist)

  s = (Dist - lc_mindist) / (lc_maxdist - lc_mindist)
  dx = s*lc_max + (1-s)*lc_min

  TargetLength(1) = dx
  TargetLength(2) = dx

  IF(ExtrudeLayers) THEN

    Point(1) = xx
    Point(2) = yy
    Point(3) = 0.0_dp

    NodeDepth = GetZFromMask(Mesh, Point, BottomPerm, BottomElemMask)
    NodeElev = GetZFromMask(Mesh, Point, SurfacePerm, SurfElemMask)

    dz = (NodeElev-NodeDepth)/(layers-1)

    TargetLength(3) = dz
  ELSE
    TargetLength(3) = dz
  END IF

  IF(Debug) PRINT *,'Node ',nodenumber,'TargetLength: ',TargetLength,' dist: ',Dist,' s: ',s

  CONTAINS

  FUNCTION GetZFromMask(Mesh, Point, NodePerm, ElemMask) RESULT(NodeDepth)

    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: NodePerm(:)
    LOGICAL, ALLOCATABLE :: ElemMask(:)
    REAL(KIND=dp) :: NodeDepth
    !-----------------------------------------
    TYPE(Element_t),POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp), POINTER :: ElementValues(:)
    REAL(KIND=dp), DIMENSION(3) :: NodeCoord, LocalCoordinates
    REAL(KIND=dp), DIMENSION(2) :: Distances, Depths, Weights
    REAL(KIND=dp) :: eps_global_limit, eps_local_limit,&
        eps_global, eps_local, eps_numeric
    REAL(KIND=dp) :: MinDist,MinDist2,dist
    INTEGER :: i,n,NodeIndx,NodeIndx2
    LOGICAL :: Found

    eps_global = 2.0e-10_dp
    eps_local = 1.0e-10_dp
    eps_numeric = 1.0e-10_dp
    eps_global_limit = 1.0E-2_dp
    eps_local_limit = 1.0E-2_dp

    n = Mesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), &
         ElementNodes % z(n), ElementValues(n) )

    DO WHILE(.TRUE.)
      DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

        IF(ElemMask(i)) CYCLE

        Element => Mesh % Elements(i)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z = 0.0_dp
        Found = PointInElement( Element, ElementNodes, &
            Point, LocalCoordinates, eps_global, eps_local, eps_numeric)
        IF( Found ) EXIT
      END DO
      IF( Found ) EXIT

      eps_global = eps_global * 10.0_dp
      eps_local = eps_local * 10.0_dp
      IF(eps_global > eps_global_limit) EXIT
      IF(eps_local > eps_local_limit) EXIT
    END DO

    IF(.NOT. Found) THEN

      MinDist = HUGE(1.0_dp); MinDist2 = HUGE(1.0_dp)
      DO i=1,NNodes

        IF(NodePerm(i) == 0) CYCLE

        NodeCoord(1) = Mesh % Nodes % x(i)
        NodeCoord(2) = Mesh % Nodes % y(i)
        NodeCoord(3) = 0.0_dp

        Dist = SUM( ( Point(:2) - NodeCoord(:2) )**2 )
        IF( Dist < MinDist ) THEN
          IF(MinDist < MinDist2) THEN
            MinDist2 = MinDist
            NodeIndx2 = NodeIndx
          END IF
          MinDist = Dist
          NodeIndx = i
        ELSE IF( Dist < MinDist2) THEN
          MinDist2 = Dist
          NodeIndx2 = i
        END IF
      END DO

      Distances(1) = MinDist ** 0.5
      Distances(2) = MinDist2 ** 0.5

      Weights(1) = Distances(1) ** (-1.0_dp)
      Weights(2) = Distances(2) ** (-1.0_dp)

      Depths(1) = Mesh % Nodes % z(NodeIndx) * Weights(1)
      Depths(2) = Mesh % Nodes % z(NodeIndx2) * Weights(2)

      NodeDepth = SUM(Depths)/SUM(Weights)
    ELSE
      ElementValues(1:n) = Mesh % Nodes % z(NodeIndexes)

      NodeDepth = InterpolateInElement( &
          Element, ElementValues, LocalCoordinates(1), &
          LocalCoordinates(2), LocalCoordinates(3) )
    END IF
  END FUNCTION GetZFromMask

END SUBROUTINE GlacierMeshMetricAniso
