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
! *  FreeSurface utilities
! !  TODO: Not quite finished
! !  TODO1: outdated should be removed... (12.12.2003, Juha)
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Jun 1997
! *
! ****************************************************************************/

!> Internal free surface utilities. 
!> \deprecated
!> \ingroup ElmerLib
!> \{

MODULE FreeSurface

   USE DirectSolve
   USE IterSolve
   USE ElementUtils

   IMPLICIT NONE

CONTAINS
!-------------------------------------------------------------------------------
!
! 
   SUBROUTINE MeanCurvature( Model )
!-------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
!-------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,t,BC,NodeIndexes(16)
    INTEGER, POINTER :: Reorder(:),Visited(:)
    LOGICAL :: L

    REAL(KIND=dp) :: ddxddu,ddyddu,ddzddu,ddxdudv,ddydudv,ddzdudv, &
      ddxddv,ddyddv,ddzddv,Auu,Auv,Avv,Buu,Buv,Bvv,detA,u,v,x,y,z

    REAL(KIND=dp), ALLOCATABLE :: dxdu(:),dydu(:),dzdu(:)
    REAL(KIND=dp), ALLOCATABLE :: dxdv(:),dydv(:),dzdv(:)
 
    REAL(KIND=dp), TARGET :: nx(16),ny(16),nz(16),Nrm(3)
    REAL(KIND=dp), POINTER :: Curvature(:)

    REAL(KIND=dp) :: Basis(16),dBasisdx(16,3),ddBasisddx(16,3,3)
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Boundary

    LOGICAL :: stat

    REAL(KIND=dp) :: Metric(3,3),SqrtMetric,Symbols(3,3,3),dSymbols(3,3,3,3)
!-------------------------------------------------------------------------------
    ALLOCATE( Visited(Model % NumberOfNodes) )

    IF ( .NOT.ASSOCIATED( Model % FreeSurfaceNodes) ) THEN
      ALLOCATE( Model % FreeSurfaceNodes( Model % NumberOfNodes ) )
      Model % FreeSurfaceNodes = 0
    END IF
    Reorder => Model % FreeSurfaceNodes

    ! Count nodes on free boundaries, and make a permutation table for nodes,
    ! should perhaps be saved instead of recomputed every time:
    !------------------------------------------------------------------------
    Visited = 0
    Reorder = 0
    k = 0
    DO BC = 1,Model % NumberOfBCs
      IF (.NOT. ListGetLogical(Model % BCs(BC) % Values, 'Free Surface',L)) CYCLE

      DO t=Model % NumberOfBulkElements+1,Model % NumberOfBulkElements + &
                   Model % NumberOfBoundaryElements

        Boundary => Model % Elements(t)
        IF ( Boundary % BoundaryInfo % Constraint /= BC ) CYCLE
  
        n = Boundary % Type % NumberOfNodes
        DO i=1,n
          j = Boundary % NodeIndexes(i)
          IF ( Visited(j) == 0 ) THEN
            k = k + 1
            Reorder(j) = k
          END IF
          Visited(j) = Visited(j) + 1
        END DO
      END DO
    END DO 
!-------------------------------------------------------------------------------
! If no free surfaces, return
!-------------------------------------------------------------------------------
    IF ( k == 0 ) THEN
      DEALLOCATE( Visited )
      RETURN
    END IF
!-------------------------------------------------------------------------------
! Allocate some memories...
!-------------------------------------------------------------------------------
    IF ( .NOT.ASSOCIATED( Model % BoundaryCurvatures ) ) THEN
      ALLOCATE( Model % BoundaryCurvatures(k) )
    END IF
    Curvature => Model % BoundaryCurvatures
    Curvature = 0.0D0

    ALLOCATE( dxdu(k),dydu(k),dzdu(k),dxdv(k),dydv(k),dzdv(k) )
    dxdu = 0.0D0
    dydu = 0.0D0
    dzdu = 0.0D0
    dxdv = 0.0D0
    dydv = 0.0D0
    dzdv = 0.0D0
!-------------------------------------------------------------------------------
!   Compute sum of first derivatives on nodes of the boundaries
!-------------------------------------------------------------------------------
    DO BC = 1,Model % NumberOfBCs

      IF (.NOT. ListGetLogical(Model % BCs(BC) % Values, 'Free Surface',L)) CYCLE

      DO t=Model % NumberOfBulkElements+1,Model % NumberOfBulkElements + &
                   Model % NumberOfBoundaryElements

        Boundary => Model % Elements(t)
        IF ( Boundary % BoundaryInfo % Constraint /= BC ) CYCLE

        n = Boundary % Type % NumberOfNodes
        nx(1:n) = Model % Nodes % x(Boundary % NodeIndexes)
        ny(1:n) = Model % Nodes % y(Boundary % NodeIndexes)
        nz(1:n) = Model % Nodes % z(Boundary % NodeIndexes)
#if 1
nodes % x => nx(1:n)
nodes % y => ny(1:n)
nodes % z => nz(1:n)
#endif

        DO i=1,n
          j = Reorder(Boundary % NodeIndexes(i))

          IF ( Boundary % Type % Dimension == 1 ) THEN
            u = Boundary % Type % NodeU(i)

#if 0
            dxdu(j) = dxdu(j) + FirstDerivative1D( Boundary,nx(1:n),u )
            dydu(j) = dydu(j) + FirstDerivative1D( Boundary,ny(1:n),u )
#else
stat = ElementInfo( Boundary, Nodes, u, 0.0d0, 0.0d0, detA, &
     Basis, dBasisdx )
dydu(j) = dydu(j) + SUM( dBasisdx(1:n,1) * ny(1:n) )
#endif
          ELSE
            u = Boundary % Type % NodeU(i)
            v = Boundary % Type % NodeV(i)

            dxdu(j) = dxdu(j) + FirstDerivativeInU2D( Boundary,nx(1:n),u,v )
            dydu(j) = dydu(j) + FirstDerivativeInU2D( Boundary,ny(1:n),u,v )
            dzdu(j) = dzdu(j) + FirstDerivativeInU2D( Boundary,nz(1:n),u,v )

            dxdv(j) = dxdv(j) + FirstDerivativeInV2D( Boundary,nx(1:n),u,v )
            dydv(j) = dydv(j) + FirstDerivativeInV2D( Boundary,ny(1:n),u,v )
            dzdv(j) = dzdv(j) + FirstDerivativeInV2D( Boundary,nz(1:n),u,v )
          END IF
        END DO
      END DO
    END DO

!-------------------------------------------------------------------------------
!   Average of the derivatives at nodes...
!-------------------------------------------------------------------------------
    DO BC=1,Model % NumberOfBCs
      IF (.NOT. ListGetLogical(Model % BCs(BC) % Values, 'Free Surface',L)) CYCLE

      DO t=Model % NumberOfBulkElements+1,Model % NumberOfBulkElements + &
                   Model % NumberOfBoundaryElements

        Boundary => Model % Elements(t)
        IF ( Boundary % BoundaryInfo % Constraint /= BC ) CYCLE

        n = Boundary % Type % NumberOfNodes
        DO i=1,n
          j = Boundary % NodeIndexes(i)
          IF ( Visited(j) > 0 ) THEN
            dxdu(Reorder(j)) = dxdu(Reorder(j)) / Visited(j)
            dydu(Reorder(j)) = dydu(Reorder(j)) / Visited(j)
            dzdu(Reorder(j)) = dzdu(Reorder(j)) / Visited(j)
            dxdv(Reorder(j)) = dxdu(Reorder(j)) / Visited(j)
            dydv(Reorder(j)) = dydv(Reorder(j)) / Visited(j)
            dzdv(Reorder(j)) = dzdv(Reorder(j)) / Visited(j)
            Visited(j) = -Visited(j)
          END IF
        END DO
      END DO
    END DO

    Visited = ABS( Visited )

!-------------------------------------------------------------------------------
!   The curvature computation begins
!-------------------------------------------------------------------------------
    DO BC=1,Model % NumberOfBCs

      IF (.NOT. ListGetLogical(Model % BCs(BC) % Values, 'Free Surface',L)) CYCLE

      DO t=Model % NumberOfBulkElements+1,Model % NumberOfBulkElements + &
                   Model % NumberOfBoundaryElements

        Boundary => Model % Elements(t)

        IF ( Boundary % BoundaryInfo % Constraint /= BC ) CYCLE

        n = Boundary % Type % NumberOfNodes
        NodeIndexes(1:n) = Boundary % NodeIndexes
!-------------------------------------------------------------------------------
!       Go through element nodal points
!-------------------------------------------------------------------------------
#if 1
nx(1:n) = Model % Nodes % x(NodeIndexes(1:n))
ny(1:n) = Model % Nodes % y(NodeIndexes(1:n))
nz(1:n) = Model % Nodes % z(NodeIndexes(1:n))
nodes % x => nx(1:n)
nodes % y => ny(1:n)
nodes % z => nz(1:n)
#endif
        DO i=1,n
          j = Reorder( NodeIndexes(i) )

          u = Boundary % Type % NodeU(i)
          v = 0.0D0
          IF ( Boundary % Type % Dimension > 1 ) v = Boundary % Type % NodeV(i)

          SqrtMetric = 1.0D0
          x = Model % Nodes % x(NodeIndexes(i))
          y = Model % Nodes % y(NodeIndexes(i))
          z = Model % Nodes % z(NodeIndexes(i))
  
          IF  (CurrentCoordinateSystem() /= Cartesian ) THEN
            CALL CoordinateSystemInfo( Metric,SqrtMetric,Symbols,dSymbols,X,Y,Z )
          END IF
!-------------------------------------------------------------------------------
!       2D case, compute the curvature
!-------------------------------------------------------------------------------
          IF ( Boundary % Type % Dimension == 1 ) THEN
#if 0
!-------------------------------------------------------------------------------
!           Second partial derivatives of the space coordinates with respect to
!           curve coordinate
!-------------------------------------------------------------------------------
            ddxddu = FirstDerivative1D( Boundary,dxdu(Reorder(NodeIndexes(1:n))),u )
            ddyddu = FirstDerivative1D( Boundary,dydu(Reorder(NodeIndexes(1:n))),u )
!-------------------------------------------------------------------------------
!           curve 'metric'
!-------------------------------------------------------------------------------
            Auu  = dxdu(j)*dxdu(j) + dydu(j)*dydu(j)
            detA = 1.0d0 / SQRT(Auu)

            Nrm(1) = -dydu(j)
            Nrm(2) =  dxdu(j)
            Nrm(3) =  0.0D0
            Nrm = Nrm / SQRT( SUM( Nrm**2 ) )
            CALL CheckNormalDirection( Boundary,Nrm,x,y,z )
!-------------------------------------------------------------------------------
!           and finally the curvature
!-------------------------------------------------------------------------------
            Curvature(j) = Curvature(j) + &
                0.5d0 * Auu * ( ddxddu*Nrm(1) + ddyddu*Nrm(2) )
#else
stat = ElementInfo( Boundary, Nodes, u, 0.0d0, 0.0d0, detA, &
     Basis, dBasisdx, ddBasisddx, .TRUE. )
ddyddu = SUM( dBasisdx(1:n,1) * dydu(Reorder(NodeIndexes(1:n))) )
Curvature(j) = Curvature(j) + ( y * ddyddu - (1+dydu(j)**2) ) / ( y * ( 1+dydu(j)**2 )**(3.0d0/2.0d0) )
#endif


          ELSE
!-------------------------------------------------------------------------------
!        3D case, compute the curvature
!-------------------------------------------------------------------------------
            u = Boundary % Type % NodeU(i)
            v = Boundary % Type % NodeV(i)
!-------------------------------------------------------------------------------
!           Second partial derivatives of the space coordinates with respect to
!           surface coordinates
!-------------------------------------------------------------------------------
            ddxddu  = FirstDerivativeInU2D( Boundary, dxdu(NodeIndexes(1:n)),u,v )
            ddyddu  = FirstDerivativeInU2D( Boundary, dydu(NodeIndexes(1:n)),u,v )
            ddzddu  = FirstDerivativeInU2D( Boundary, dzdu(NodeIndexes(1:n)),u,v )

            ddxdudv = FirstDerivativeInU2D( Boundary, dxdv(NodeIndexes(1:n)),u,v )
            ddydudv = FirstDerivativeInU2D( Boundary, dydv(NodeIndexes(1:n)),u,v )
            ddzdudv = FirstDerivativeInU2D( Boundary, dzdv(NodeIndexes(1:n)),u,v )

            ddxddv  = FirstDerivativeInV2D( Boundary, dxdv(NodeIndexes(1:n)),u,v )
            ddyddv  = FirstDerivativeInV2D( Boundary, dydv(NodeIndexes(1:n)),u,v )
            ddzddv  = FirstDerivativeInV2D( Boundary, dzdv(NodeIndexes(1:n)),u,v )
!-------------------------------------------------------------------------------
!           Surface metric
!-------------------------------------------------------------------------------
            Auu = dxdu(k)*dxdu(k) + dydu(k)*dydu(k) + dzdu(k)*dzdu(k)
            Auv = dxdu(k)*dxdv(k) + dydu(k)*dydv(k) + dzdu(k)*dzdv(k)
            Avv = dxdv(k)*dxdv(k) + dydv(k)*dydv(k) + dzdv(k)*dzdv(k)

            detA = 1.0D0 / SQRT(Auu*Avv - Auv*Auv)
!-------------------------------------------------------------------------------
!           Change metric to contravariant form
!-------------------------------------------------------------------------------
            u = Auu
            Auu =  Avv * detA
            Auv = -Auv * detA
            Avv =    u * detA
!-------------------------------------------------------------------------------
!           Normal vector to surface
!-------------------------------------------------------------------------------
            Nrm(1) = (dydu(k) * dzdv(k) - dydv(k) * dzdu(k)) * detA
            Nrm(2) = (dxdv(k) * dzdu(k) - dxdu(k) * dzdv(k)) * detA
            Nrm(3) = (dxdu(k) * dydv(k) - dxdv(k) * dydu(k)) * detA
            Nrm = Nrm / SQRT( SUM( Nrm**2 ) )

            CALL CheckNormalDirection( Boundary,Nrm,x,y,z )
!-------------------------------------------------------------------------------
!           The second fundamental form of the surface
!-------------------------------------------------------------------------------
            Buu = ddxddu  * Nrm(1) + ddyddu  * Nrm(2) + ddzddu  * Nrm(3)
            Buv = ddxdudv * Nrm(1) + ddydudv * Nrm(2) + ddzdudv * Nrm(3)
            Bvv = ddxddv  * Nrm(1) + ddyddv  * Nrm(2) + ddzddv  * Nrm(3)
!-------------------------------------------------------------------------------
!           And finally, the curvature
!-------------------------------------------------------------------------------
            Curvature(j) = Curvature(j) + (Auu*Buu + 2*Auv*Buv + Avv*Bvv)
          END IF
        END DO
      END DO
    END DO
visited=abs(visited)
!-------------------------------------------------------------------------------
! Average to nodes
!-------------------------------------------------------------------------------
    DO BC=1,Model % NumberOfBCs
      IF (.NOT. ListGetLogical(Model % BCs(BC) % Values, 'Free Surface',L)) CYCLE

      DO t=Model % NumberOfBulkElements+1,Model % NumberOfBulkElements + &
                   Model % NumberOfBoundaryElements

        Boundary => Model % Elements(t)
        IF ( Boundary % BoundaryInfo % Constraint /= BC ) CYCLE

        n = Boundary % Type % NumberOfNodes
        DO i=1,n
          j = Boundary % NodeIndexes(i)
          IF ( Visited(j) > 0 ) THEN
            Curvature(Reorder(j)) = Curvature(Reorder(j)) / Visited(j)
write(*,*) j,curvature(Reorder(j))
            Visited(j) = -Visited(j)
          END IF
        END DO
      END DO
    END DO
!-------------------------------------------------------------------------------
!   We are done here.
!-------------------------------------------------------------------------------
    DEALLOCATE( Visited,dxdu,dydu,dzdu,dxdv,dydv,dzdv )
print*,'----------------------'
!-------------------------------------------------------------------------------
  END SUBROUTINE MeanCurvature
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
  SUBROUTINE MoveBoundary( Model,Relax )
!-------------------------------------------------------------------------------
    USE DefUtils
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: Relax
!-------------------------------------------------------------------------------
    INTEGER :: i,ii,j,k,m,n,t,Which, FIter
    LOGICAL, ALLOCATABLE :: Visited(:),Turned(:)
    LOGICAL :: L

    REAL(KIND=dp) :: Auu,Auv,Avv,Buu,Buv,Bvv,detA,u,v,x,y,z,S,R,Ux,Uy,Uz

    REAL(KIND=dp) :: dxdu,dydu,dzdu, FEPS
    REAL(KIND=dp) :: dxdv,dydv,dzdv,x1,x2,y1,y2,z1,z2
 
    REAL(KIND=dp), TARGET :: N1,N2,N3,Nrm(3)
    REAL(KIND=dp), POINTER :: Curvature(:)

    TYPE(Element_t), POINTER :: Boundary, Element
    TYPE(Nodes_t) :: BoundaryNodes, ElementNodes

    REAL(KIND=dp), ALLOCATABLE :: XCoord(:),YCoord(:),ZCoord(:), AvarageNormal(:,:)
    LOGICAL :: XMoved, YMoved, ZMoved

    TYPE( Variable_t), POINTER :: Velocity1,Velocity2,Velocity3

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: BC

    REAL(KIND=dp) :: SqrtMetric, dLBasisdx(16,2),NodalBasis(16)

    SAVE BoundaryNodes, ElementNodes
!-------------------------------------------------------------------------------

    FEPS = ListGetConstReal( Model % Solver % Values, &
          'Free Surface Convergence Tolerance', L )
    IF ( .NOT. L ) FEPS = 1.0d-12

    FIter = ListGetInteger( Model % Solver % Values, &
          'Free Surface Max Iterations', L )
    IF ( .NOT. L ) FIter = 100

    Velocity1 => VariableGet( Model % Variables, 'Velocity 1' )
    Velocity2 => VariableGet( Model % Variables, 'Velocity 2' )
    Velocity3 => VariableGet( Model % Variables, 'Velocity 3' )
!-------------------------------------------------------------------------------
    ALLOCATE( Visited(Model % NumberOfNodes),Turned(Model % NumberOfNodes) )
    Visited = .FALSE.
    Turned  = .FALSE.
!-------------------------------------------------------------------------------

    ALLOCATE( XCoord(Model % NumberOfNodes) )
    XMoved = .FALSE.
    XCoord = Model % Nodes % x

    ALLOCATE( YCoord(Model % NumberOfNodes) )
    YMoved = .FALSE.
    YCoord = Model % Nodes % y

    ALLOCATE( ZCoord(Model % NumberOfNodes) )
    ZMoved = .FALSE.
    ZCoord = Model % Nodes % z

    ALLOCATE( AvarageNormal(Model % NumberOfBCs,3) )
    AvarageNormal = 0

    Ux = 0; Uy = 0; Uz = 0;
!-------------------------------------------------------------------------------
!   Check normal direction first in a separate loop, cause within the node
!   moving loop, the parent elements might be a little funny...!
!-------------------------------------------------------------------------------
    Mesh => Model % Solver % Mesh
    DO t = 1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement() ) CYCLE

      BC => GetBC()
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE
      IF ( .NOT.GetLogical(BC, 'Free Surface',L) ) CYCLE

      n = GetElementNOFNodes()
      CALL GetElementNodes( BoundaryNodes )
!-------------------------------------------------------------------------------
!     Go through boundary element nodes
!-------------------------------------------------------------------------------
      DO i=1,n
        k = Boundary % NodeIndexes(i)

        IF ( Boundary % Type % Dimension == 1 ) THEN
          IF ( i==1 ) THEN
             j=2
          ELSE
             j=1
          END IF
          j = Boundary % NodeIndexes(j)
          x1  = XCoord(j) - XCoord(k)
          y1  = YCoord(j) - YCoord(k)

          Ux = Velocity1 % Values( Velocity1 % Perm(k) )
          Uy = Velocity2 % Values( Velocity2 % Perm(k) )

          IF ( Ux*x1 + Uy*y1 > 0 ) CYCLE
        END IF

!-------------------------------------------------------------------------------
!       Should not move the same node twice,so check if already done
!-------------------------------------------------------------------------------
        IF ( .NOT.Visited(k) ) THEN
            Visited(k) = .TRUE.
!-------------------------------------------------------------------------------
!           2D case, compute normal
!-------------------------------------------------------------------------------
            IF ( Boundary % Type % Dimension == 1 ) THEN
              u = Boundary % Type % NodeU(i)
!------------------------------------------------------------------------------
!             Basis function derivatives with respect to local coordinates
!------------------------------------------------------------------------------
              NodalBasis(1:n) = 0.0d0
              DO m=1,n
                NodalBasis(m)  = 1.0d0
                dLBasisdx(m,1) = FirstDerivative1D( Boundary,NodalBasis,u )
                NodalBasis(m)  = 0.0d0
              END DO

              Nrm(1) = -SUM( BoundaryNodes % y(1:n)*dLBasisdx(1:n,1) )
              Nrm(2) =  SUM( BoundaryNodes % x(1:n)*dLBasisdx(1:n,1) )
              Nrm(3) =  0.0D0
            ELSE
!-------------------------------------------------------------------------------
!            3D case, compute normal
!-------------------------------------------------------------------------------
              u = Boundary % Type % NodeU(i)
              v = Boundary % Type % NodeV(i)

              NodalBasis(1:n) = 0.0D0
              DO m=1,N
                NodalBasis(m)  = 1.0D0
                dLBasisdx(m,1) = FirstDerivativeInU2D( Boundary,NodalBasis,u,v )
                dLBasisdx(m,2) = FirstDerivativeInV2D( Boundary,NodalBasis,u,v )
                NodalBasis(m)  = 0.0D0
              END DO

              dxdu = SUM( BoundaryNodes % x(1:n) * dLBasisdx(1:n,1) )
              dydu = SUM( BoundaryNodes % y(1:n) * dLBasisdx(1:n,1) )
              dzdu = SUM( BoundaryNodes % z(1:n) * dLBasisdx(1:n,1) )

              dxdv = SUM( BoundaryNodes % x(1:n) * dLBasisdx(1:n,2) )
              dydv = SUM( BoundaryNodes % y(1:n) * dLBasisdx(1:n,2) )
              dzdv = SUM( BoundaryNodes % z(1:n) * dLBasisdx(1:n,2) )

              Nrm(1) = dydu * dzdv - dydv * dzdu
              Nrm(2) = dxdv * dzdu - dxdu * dzdv
              Nrm(3) = dxdu * dydv - dxdv * dydu
            END IF
!-------------------------------------------------------------------------------
!           Turn the normal to point outwards, or towards less dense material
!-------------------------------------------------------------------------------
            x = Model % Nodes % x(k)
            y = Model % Nodes % y(k)
            z = Model % Nodes % z(k)
            CALL CheckNormalDirection( Boundary,Nrm,x,y,z,Turned(k) )

            m = Boundary % BoundaryInfo % Constraint
            R = SQRT(SUM(Nrm**2))
            AvarageNormal(m,1:3) = AvarageNormal(m,1:3) + ABS(Nrm)/R
          END IF
       END DO
    END DO

!-------------------------------------------------------------------------------
!   Iterate until convergence
!-------------------------------------------------------------------------------
    DO ii=1,FIter
       S = 0.0d0
       Visited = .FALSE.

       DO t=1,Mesh % NumberOfBoundaryElements
         Boundary => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE
   
         BC => GetBC()
         IF (.NOT.GetLogical(BC, 'Free Surface',L)) CYCLE

         i = CoordinateSystemDimension()
         Which = ListGetInteger( BC, 'Free Coordinate', L,minv=1, maxv=i )
         IF ( .NOT. L ) THEN
           i =  Boundary % BoundaryInfo % Constraint
           Nrm = AvarageNormal(i,:)
           Which = 1
           DO i=2,3
             IF ( Nrm(i) > Nrm(Which) ) Which=i
           END DO
         END IF

         n = GetElementNOFNOdes()
         CALL GetElementNodes( BoundaryNodes )

         DO i=1,n
           k = Boundary % NodeIndexes(i)

           IF ( Boundary % Type % Dimension == 1 ) THEN
             IF ( i==1 ) THEN
                j=2
             ELSE
               j=1
             END IF
             j = Boundary % NodeIndexes(j)
             x1  = XCoord(j) - XCoord(k)
             y1  = YCoord(j) - YCoord(k)

             Ux = Velocity1 % Values( Velocity1 % Perm(k) )
             Uy = Velocity2 % Values( Velocity2 % Perm(k) )

             IF ( Ux*x1 + Uy*y1 > 0 ) CYCLE
           END IF

!-------------------------------------------------------------------------------
!          Should not move the same node twice,so check if already done
!-------------------------------------------------------------------------------
           IF ( .NOT.Visited(k) ) THEN
             Visited(k) = .TRUE.
!-------------------------------------------------------------------------------
!            2D case, compute normal
!-------------------------------------------------------------------------------
             IF ( Boundary % Type % Dimension == 1 ) THEN
               u = Boundary % Type % NodeU(i)
!------------------------------------------------------------------------------
!              Basis function derivatives with respect to local coordinates
!------------------------------------------------------------------------------
               NodalBasis(1:n) = 0.0D0
               DO m=1,n
                 NodalBasis(m)  = 1.0D0
                 dLBasisdx(m,1) = FirstDerivative1D( Boundary,NodalBasis,u )
                 NodalBasis(m)  = 0.0D0
               END DO

               Nrm(1) = -SUM( BoundaryNodes % y(1:n)*dLBasisdx(1:n,1) )
               Nrm(2) =  SUM( BoundaryNodes % x(1:n)*dLBasisdx(1:n,1) )
               Nrm(3) =  0.0D0
             ELSE
!----------------------------------------------------------------------------
!             3D case, compute normal
!----------------------------------------------------------------------------
               u = Boundary % Type % NodeU(i)
               v = Boundary % Type % NodeV(i)

               NodalBasis(1:n) = 0.0D0
               DO j=1,n
                 NodalBasis(j)  = 1.0D0
                 dLBasisdx(j,1) = FirstDerivativeInU2D( Boundary,NodalBasis,u,v )
                 dLBasisdx(j,2) = FirstDerivativeInV2D( Boundary,NodalBasis,u,v )
                 NodalBasis(j)  = 0.0D0
               END DO

               dxdu = SUM( BoundaryNodes % x(1:n) * dLBasisdx(1:n,1) )
               dydu = SUM( BoundaryNodes % y(1:n) * dLBasisdx(1:n,1) )
               dzdu = SUM( BoundaryNodes % z(1:n) * dLBasisdx(1:n,1) )

               dxdv = SUM( BoundaryNodes % x(1:n) * dLBasisdx(1:n,2) )
               dydv = SUM( BoundaryNodes % y(1:n) * dLBasisdx(1:n,2) )
               dzdv = SUM( BoundaryNodes % z(1:n) * dLBasisdx(1:n,2) )

               Nrm(1) = dydu * dzdv - dydv * dzdu
               Nrm(2) = dxdv * dzdu - dxdu * dzdv
               Nrm(3) = dxdu * dydv - dxdv * dydu
             END IF
!-------------------------------------------------------------------------------
!            Turn the normal to point outwards, or towards less dense material
!-------------------------------------------------------------------------------
             IF ( Turned(k) ) Nrm = -Nrm
!-------------------------------------------------------------------------------
!            Now then, lets move the node so that u.n will be reduced.
!-------------------------------------------------------------------------------
             IF ( Boundary % Type % Dimension == 1 ) THEN
!-------------------------------------------------------------------------------
!              TODO ::  This won't handle the three node line
!              2D case, move the nodes..
!-------------------------------------------------------------------------------
               R = Ux * Nrm(1) + Uy * Nrm(2)
               IF ( Which == 2 ) THEN
                 IF ( ABS(Ux) > AEPS ) THEN
                   IF ( Turned(k) ) THEN
                       Model % Nodes % y(k) = Model % Nodes % y(k) - &
                                 R / ( Ux*dLBasisdx(i,1) )
                   ELSE
                       Model % Nodes % y(k) = Model % Nodes % y(k) + &
                                 R / ( Ux*dLBasisdx(i,1) )
                   END IF
                   YMoved = .TRUE.
                 END IF
               ELSE
                 IF ( ABS(Uy) > AEPS ) THEN
                   IF ( Turned(k) ) THEN
                      Model % Nodes % x(k) = Model % Nodes % x(k) + &
                                R / ( Uy*dLBasisdx(i,1) )
                   ELSE
                      Model % Nodes % x(k) = Model % Nodes % x(k) - &
                                R / ( Uy*dLBasisdx(i,1) )
                   END IF
                   XMoved = .TRUE.
                 END IF
               END IF
             ELSE
!-------------------------------------------------------------------------------
!              3D case, move the nodes..
!              TODO :: This is just guesswork, no testing done...
!----------------------------------------------------------------------------
               Ux = Velocity1 % Values( Velocity1 % Perm(k) )
               Uy = Velocity2 % Values( Velocity2 % Perm(k) )
               Uz = Velocity3 % Values( Velocity3 % Perm(k) )

               R = Ux * Nrm(1) + Uy * Nrm(2) + Uz * Nrm(3)
               IF ( Which == 1 ) THEN
                 IF ( ABS(Uy) > AEPS .OR. ABS(Uz) > AEPS ) THEN
                   Model % Nodes % x(k) = Model % Nodes % x(k) + R / &
                    ( (dzdu*dLBasisdx(i,2) - dzdv*dLBasisdx(i,1))*Uy + &
                      (dydv*dLBasisdx(i,1) - dydu*dLBasisdx(i,2))*Uz )
                   XMoved = .TRUE.
                 END IF
               ELSE IF ( Which == 2 ) THEN
                 IF ( ABS(Ux) > AEPS .OR. ABS(Uz) > AEPS ) THEN
                   Model % Nodes % y(k) = Model % Nodes % y(k) + R / &
                    ( (dzdv*dLBasisdx(i,1) - dzdu*dLBasisdx(i,2))*Ux + &
                      (dxdu*dLBasisdx(i,2) - dxdv*dLBasisdx(i,1))*Uz )
                   YMoved = .TRUE.
                 END IF
               ELSE
                 IF ( ABS(Ux) > AEPS .OR. ABS(Uy) > AEPS ) THEN
                   Model % Nodes % z(k) = Model % Nodes % z(k) + R / &
                    ( (dydu*dLBasisdx(i,2) - dydv*dLBasisdx(i,1))*Ux + &
                      (dxdv*dLBasisdx(i,1) - dxdu*dLBasisdx(i,2))*Uy )
                   ZMoved = .TRUE.
                 END IF
               END IF
             END IF
!-------------------------------------------------------------------------------
             S = S + ( (Ux*Nrm(1)+Uy*Nrm(2)+Uz*Nrm(3))/SQRT(SUM(Nrm**2)) )**2
!-------------------------------------------------------------------------------
           END IF
        END DO
      END DO
!-------------------------------------------------------------------------------
      S = SQRT(S)
      WRITE( Message, '(a,i4,a,e21.12)' ) 'Iter: ', ii, ' Free surface Residual: ',S
      CALL Info( 'FreeSurface', Message, Level=4 )
      IF ( S < FEPS ) EXIT
    END DO
!-------------------------------------------------------------------------------
    IF ( XMoved ) THEN
      CALL PoissonSolve( Model,XCoord,YCoord,ZCoord,Model % Nodes % x,1 )
      IF ( Relax /= 1.0d0 ) &
        Model % Nodes % x = Relax * Model % Nodes % x + (1-Relax)*XCoord
    END IF

    IF ( YMoved ) THEN
      CALL PoissonSolve( Model,XCoord,YCoord,ZCoord,Model % Nodes % y,2 )
      IF ( Relax /= 1.0d0 ) &
        Model % Nodes % y = Relax * Model % Nodes % y + (1-Relax)*YCoord
    END IF

    IF ( ZMoved ) THEN
      CALL PoissonSolve( Model,XCoord,YCoord,ZCoord,Model % Nodes % z,3 )
      IF ( Relax /= 1.0d0 ) &
        Model % Nodes % z = Relax * Model % Nodes % z + (1-Relax)*ZCoord
    END IF

    DEALLOCATE( Visited,Turned,XCoord,YCoord,ZCoord,AvarageNormal )
 
    ALLOCATE( ElementNodes % x(Model % MaxElementNodes) )
    ALLOCATE( ElementNodes % y(Model % MaxElementNodes) )
    ALLOCATE( ElementNodes % z(Model % MaxElementNodes) )

    DO i=1,Model % NumberOfBulkElements
      Element => Model % Elements(i)
      n = Element % Type % NumberOfNodes
      ElementNodes % x(1:n) = Model % Nodes % x(Element % NodeIndexes)
      ElementNodes % y(1:n) = Model % Nodes % y(Element % NodeIndexes)
      ElementNodes % z(1:n) = Model % Nodes % z(Element % NodeIndexes)
      CALL StabParam( Element, ElementNodes, n, &
               Element % StabilizationMK, Element % hK )
    END DO
   
    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z)
!-------------------------------------------------------------------------------
  END SUBROUTINE MoveBoundary
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
  SUBROUTINE PoissonSolve( Model,NX,NY,NZ,Solution,Moved )
!-------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Moved
    REAL(KIND=dp) :: NX(:),NY(:),NZ(:),Solution(:)
!-------------------------------------------------------------------------------

    TYPE(Matrix_t), POINTER :: CMatrix
    INTEGER, POINTER :: CPerm(:)

    LOGICAL :: FirstTime = .TRUE.
    REAL(KIND=dp), ALLOCATABLE :: ForceVector(:)

    REAL(KIND=dp) :: Basis(16),dBasisdx(16,3),ddBasisddx(16,3,3)
    REAL(KIND=dp) :: SqrtElementMetric,u,v,w,s,A
    LOGICAL :: Stat

    INTEGER :: i,j,k,n,q,p,t,dim

    TYPE(Solver_t), POINTER :: Solver
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element

    INTEGER :: N_Integ

    REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:), &
                       LocalMatrix(:,:)

    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

    SAVE CMatrix,CPerm,FirstTime,ForceVector,Solver,Nodes,LocalMatrix
!-------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      FirstTime = .FALSE.

      ALLOCATE( Solver )
      Solver % Values => ListAllocate()
      Solver % Mesh => CurrentModel % Mesh

      CALL ListAddString( Solver % Values, &
                      'Linear System Iterative Method', 'CGS' )
      CALL ListAddInteger( Solver % Values, &
                          'Linear System Max Iterations', 500 )
      CALL ListAddConstReal( Solver % Values, &
                'Linear System Convergence Tolerance', 1.0D-9 )
      CALL ListAddString( Solver % Values, &
                      'Linear System Preconditioning', 'ILU0' )
      CALL ListAddInteger( Solver % Values, &
                           'Linear System Residual Output', 1 )

      ALLOCATE( ForceVector( Model % NumberOfNodes ),  &
            CPerm( Model % NumberOfNodes) )

      CMatrix => CreateMatrix( CurrentModel,Solver,Solver % Mesh,CPerm,1,MATRIX_CRS,.FALSE. )

      ALLOCATE( LocalMatrix( Model % MaxElementNodes,Model % MaxElementNodes ) )
      ALLOCATE( Nodes % x(Model % MaxElementNodes) )
      ALLOCATE( Nodes % y(Model % MaxElementNodes) )
      ALLOCATE( Nodes % z(Model % MaxElementNodes) )

      Solver % TimeOrder = 0
    END IF

!------------------------------------------------------------------------------
 
    CALL CRS_ZeroMatrix( CMatrix )
    ForceVector = 0.0D0

!------------------------------------------------------------------------------
    DO i=1,Model % NumberOfBulkElements
!------------------------------------------------------------------------------
      Element => Model % Elements(i)

      dim = Element % Type % Dimension
      n = Element % Type % NumberOfNodes

      Nodes % x(1:n) = nx( Element % NodeIndexes )
      Nodes % y(1:n) = ny( Element % NodeIndexes )
      Nodes % z(1:n) = nz( Element % NodeIndexes )

      IntegStuff = GaussPoints( Element )
      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!------------------------------------------------------------------------------
      LocalMatrix = 0.0D0
      DO t=1,N_Integ
        u = U_Integ(t)
        v = V_Integ(t)
        w = W_Integ(t)
!------------------------------------------------------------------------------
!       Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                   Basis,dBasisdx )

        s = SqrtElementMetric * S_Integ(t)
        DO p=1,N
          DO q=1,N
            A = dBasisdx(p,Moved) * dBasisdx(q,Moved)
            LocalMatrix(p,q) = LocalMatrix(p,q) + s*A
          END DO
        END DO
      END DO

      CALL CRS_GlueLocalMatrix( CMatrix,n,1,Element % NodeIndexes,LocalMatrix )
    END DO

!-------------------------------------------------------------------------------
    DO t = Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
                     Model % NumberOfBoundaryElements
!-------------------------------------------------------------------------------
      Element => Model % Elements(t)
      n = Element % Type % NumberOfNodes

      DO i=1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
          IF ( .NOT.ListGetLogical(Model % BCs(i) % Values,'Free Moving',stat) ) THEN
            DO j=1,n
              k = Element % NodeIndexes(j)

              ForceVector(k) = Solution(k)
              CALL CRS_ZeroRow( CMatrix,k )
              CALL CRS_SetMatrixElement( CMatrix,k,k,1.0D0 )
            END DO
          END IF
        END IF
      END DO
!-------------------------------------------------------------------------------
    END DO
!-------------------------------------------------------------------------------

    CALL IterSolver( CMatrix,Solution,ForceVector,Solver )

!-------------------------------------------------------------------------------
  END SUBROUTINE PoissonSolve
!-------------------------------------------------------------------------------

END MODULE FreeSurface

!> \}
