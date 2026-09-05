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
! *  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA  02110-1301  USA
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Model lumping: reduce an elastic body to the 6x6 spring matrix seen at one
!> boundary, plus the rigid-body mass matrix.
!>
!> Six load cases are run -- three translations and three rotations -- and each
!> contributes one row or column of the lumped stiffness. Which of the two, and
!> how the case is imposed, depends on "Fix Displacement":
!>
!>   .TRUE.  (the default) prescribes a displacement on the lumping boundary and
!>           reads the reaction off the bulk matrix;
!>   .FALSE. applies a pure force or moment and integrates the response.
!>
!> The two are independent approximations of the same beam and agree to a couple
!> of percent, which is what fem/tests/beam-springs and beam-springs-fixdisp
!> check between them.
!>
!> Lifted out of StressSolve so that it is not one solver's private property --
!> nothing in it is specific to how the displacement was obtained. It needs the
!> solved displacement, the solver's bulk matrix, and a boundary flagged with
!> "Model Lumping Boundary"; anything that supplies those can use it.
!------------------------------------------------------------------------------

MODULE ModelLumping

  USE DefUtils
  USE LinearAlgebra

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ModelLumping_t, ModelLumpingInit, ModelLumpingLoads, &
      ModelLumpingDisplacements, ModelLumpingSprings

!------------------------------------------------------------------------------
!> State of one lumping run. It has to persist across the six load cases: the
!> stiffness accumulates a row or column at a time, and the geometry integrals
!> are computed once at the start and read by every case after.
!------------------------------------------------------------------------------
  TYPE :: ModelLumping_t
    REAL(KIND=dp) :: Area = 0.0_dp
    REAL(KIND=dp) :: Center(3) = 0.0_dp
    REAL(KIND=dp) :: Moments(3,3) = 0.0_dp
    REAL(KIND=dp) :: Kmat(6,6) = 0.0_dp
    REAL(KIND=dp) :: KmatMin(6,6) = 0.0_dp
    REAL(KIND=dp), ALLOCATABLE :: NodalLoads(:)
    LOGICAL, ALLOCATABLE :: NodeVisited(:)
    TYPE(Nodes_t) :: ParentNodes
    LOGICAL :: FixDisplacement = .TRUE.
  END TYPE ModelLumping_t

CONTAINS

!------------------------------------------------------------------------------
!> The file the lumped matrices are written to. One place, since five different
!> outputs are built from it by suffix.
!------------------------------------------------------------------------------
  FUNCTION LumpingFile( Solver ) RESULT( KmatFile )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    CHARACTER(LEN=MAX_NAME_LEN) :: KmatFile
!------------------------------------------------------------------------------
    LOGICAL :: Found
!------------------------------------------------------------------------------
    KmatFile = ListGetString( Solver % Values, 'Model Lumping Filename', Found )
    IF ( .NOT. Found ) KmatFile = "Kmat.dat"
!------------------------------------------------------------------------------
  END FUNCTION LumpingFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Start a lumping run: read the method, and compute the geometry of the lumping
!> boundary and the rigid body mass. Call once, before the load cases.
!------------------------------------------------------------------------------
  SUBROUTINE ModelLumpingInit( Lump, Solver, Model )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
    LOGICAL :: Found
!------------------------------------------------------------------------------
    Lump % FixDisplacement = ListGetLogical( Solver % Values, &
        'Fix Displacement', Found, DefValue = .TRUE. )

    IF ( Lump % FixDisplacement ) THEN
      CALL Info( 'ModelLumping', &
          'Using six fixed displacement to compute the spring matrix', Level=5 )
    ELSE
      CALL Info( 'ModelLumping', &
          'Using six pure forces and moments to compute the spring matrix', Level=5 )
    END IF

    Lump % Kmat = 0.0_dp
    Lump % KmatMin = 0.0_dp

    CALL CoordinateIntegrals( Lump, Model )
    CALL CartesianMass( Lump, Solver, Model )
!------------------------------------------------------------------------------
  END SUBROUTINE ModelLumpingInit
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Area, centre of area and the second moments of the lumping boundary.
!>
!> Two passes: the first gives the area and its centre, the second the squared
!> deviations about that centre.
!------------------------------------------------------------------------------
  SUBROUTINE CoordinateIntegrals( Lump, Model )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    REAL(KIND=dp) :: Coords(3), detJ, u, v, w, s
    INTEGER :: power, t, k, i, j, n, N_Integ, dim, maxnodes
    LOGICAL :: FoundBoundary, Found, stat

    SAVE Nodes
!------------------------------------------------------------------------------
    Mesh => Model % Mesh
    dim = CoordinateSystemDimension()
    maxnodes = Model % MaxElementNodes
    ALLOCATE( Basis(maxnodes), dBasisdx(maxnodes,3) )

    FoundBoundary = .FALSE.
    Lump % Area = 0.0_dp
    Lump % Center = 0.0_dp
    Lump % Moments = 0.0_dp

    ! On the first round compute area and center of area.
    ! On the second round compute the square deviations from the mean.

    DO power = 1,2

      DO t=1,Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        BC => GetBC()
        IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
!------------------------------------------------------------------------------
        IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE

        FoundBoundary = .TRUE.
        n = GetElementNOFNodes()
        CALL GetElementNodes( Nodes )

        IntegStuff = GaussPoints( Element )
        U_Integ => IntegStuff % u
        V_Integ => IntegStuff % v
        W_Integ => IntegStuff % w
        S_Integ => IntegStuff % s
        N_Integ =  IntegStuff % n

        DO k=1,N_Integ
          u = U_Integ(k)
          v = V_Integ(k)
          w = W_Integ(k)

          ! Basis function values & derivatives at the integration point:
          !--------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
              Basis, dBasisdx )

          s = detJ * S_Integ(k)
          IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
              CurrentCoordinateSystem() == CylindricSymmetric ) THEN
            s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
          END IF

          Coords(1) = SUM(Basis(1:n) * Nodes % x(1:n))
          IF (dim > 1) THEN
            Coords(2) =  SUM(Basis(1:n) * Nodes % y(1:n))
          END IF
          IF (dim > 2) THEN
            Coords(3) =  SUM(Basis(1:n) * Nodes % z(1:n))
          END IF

          IF(power == 1) THEN
            Lump % Area = Lump % Area + s
            Lump % Center(1:dim) = Lump % Center(1:dim) + s * Coords(1:dim)
          ELSE
            Coords(1:dim) = Coords(1:dim) - Lump % Center(1:dim)
            DO i = 1,dim
              DO j = 1,dim
                Lump % Moments(i,j) = Lump % Moments(i,j) + s * Coords(i) * Coords(j)
              END DO
            END DO
          END IF

        END DO
      END DO

      IF(.NOT. FoundBoundary) THEN
        CALL Fatal('ModelLumping','Model lumping boundary must be defined')
      END IF

      IF(power == 1) Lump % Center(1:dim) = Lump % Center(1:dim) / Lump % Area
    END DO

    DEALLOCATE( Basis, dBasisdx )
!------------------------------------------------------------------------------
  END SUBROUTINE CoordinateIntegrals
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The load of one case, as pure force or pure moment, over one boundary element.
!>
!> Pure moments may only be computed under certain conditions that should be
!> valid for boundaries with normal in the direction of some axis.
!------------------------------------------------------------------------------
  SUBROUTINE ModelLumpingLoads( Lump, Permutation, Nodes, n, Forces )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    INTEGER :: Permutation, n
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Forces(:,:)
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: y(:), z(:)
    REAL(KIND=dp) :: c, Eps
    LOGICAL :: isy, isz
    INTEGER :: ix,iy,iz
!------------------------------------------------------------------------------
    Forces = 0.0d0
    Eps = 1.0d-6

    IF(Permutation <= 3) THEN
      Forces(Permutation,1:n) = 1.0 / Lump % Area
    ELSE IF(Permutation <= 6) THEN
      ix = MOD(Permutation - 4, 3) + 1
      iy = MOD(Permutation - 3, 3) + 1
      iz = MOD(Permutation - 2, 3) + 1

      IF(Permutation == 4) THEN
        z => Nodes % Z
        y => Nodes % Y
      ELSE IF(Permutation == 5) THEN
        z => Nodes % X
        y => Nodes % Z
      ELSE IF(Permutation == 6) THEN
        z => Nodes % Y
        y => Nodes % X
      END IF

      isy = (ABS(Lump % Moments(iy,ix)) < Eps * Lump % Moments(iy,iy))
      isz = (ABS(Lump % Moments(iz,ix)) < Eps * Lump % Moments(iz,iz))

      IF(isy) THEN
        c = 1.0 / Lump % Moments(iy,iy)
        Forces(iz,1:n) = c * (y(1:n) - Lump % Center(iy))
      ELSE IF(isz) THEN
        c = -1.0 / Lump % Moments(iz,iz)
        Forces(iy,1:n) = c * (z(1:n) - Lump % Center(iz))
      ELSE
        c = 1.0 / (Lump % Moments(iy,iy) + Lump % Moments(iz,iz) )
        Forces(iy,1:n) = -c * (z(1:n) - Lump % Center(iz))
        Forces(iz,1:n) =  c * (y(1:n) - Lump % Center(iy))
        CALL Warn('ModelLumping','Moment matrix not diagonalazible!')
        PRINT *,Lump % Moments(iy,ix),Lump % Moments(iz,ix), &
            Lump % Moments(iy,iy),Lump % Moments(iz,iz)
      END IF
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ModelLumpingLoads
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Impose one load case as a prescribed displacement -- a pure translation or a
!> pure rotation of the lumping boundary.
!>
!> Called after the assembly and before the Dirichlet conditions, so that this is
!> what the boundary ends up constrained to.
!------------------------------------------------------------------------------
  SUBROUTINE ModelLumpingDisplacements( Lump, Solver, Model, Permutation )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    INTEGER :: Permutation
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Matrix_t), POINTER :: StiffMatrix
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp), POINTER :: ForceVector(:)
    INTEGER, POINTER :: Perm(:), NodeIndexes(:)
    INTEGER :: j,k,l,n,t,dim
    LOGICAL :: Found
    REAL(KIND=dp) :: Coords(3), dCoords(3), dFii, dx
!------------------------------------------------------------------------------
    Mesh => Model % Mesh
    dim = CoordinateSystemDimension()

    StiffMatrix => Solver % Matrix
    ForceVector => StiffMatrix % RHS
    Perm => Solver % Variable % Perm

    dX   = 1.0d-2*SQRT(Lump % Area)
    dFii = 1.0d-2

    DO t = 1, Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( .NOT. ActiveBoundaryElement()) CYCLE
      n = GetElementNOFNodes()

      BC => GetBC()
      IF ( .NOT.ASSOCIATED( BC ) ) CYCLE

      IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE

      NodeIndexes => Element % NodeIndexes

      DO j=1,n
        k = Perm(NodeIndexes(j))
        IF(k == 0) CYCLE

        dCoords = 0.0d0
        IF(Permutation <= 3) THEN
          dCoords(Permutation) = dX
        ELSE
          Coords(1) = Mesh % Nodes % x(NodeIndexes(j))
          Coords(2) = Mesh % Nodes % y(NodeIndexes(j))
          Coords(3) = Mesh % Nodes % z(NodeIndexes(j))
          Coords = Coords - Lump % Center
          IF (Permutation == 4) THEN
            dCoords(2) = -dFii * Coords(3)
            dCoords(3) = dFii * Coords(2)
          ELSE IF(Permutation == 5) THEN
            dCoords(1) = dFii * Coords(3)
            dCoords(3) = -dFii * Coords(1)
          ELSE IF(Permutation == 6) THEN
            dCoords(1) = -dFii * Coords(2)
            dCoords(2) = dFii * Coords(1)
          END IF

        END IF

        DO l=1,dim
          CALL SetDirichletPoint( StiffMatrix, ForceVector, l, dim, Perm, &
              NodeIndexes(j), dCoords(l) )
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ModelLumpingDisplacements
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> At the end of each load case assemble one line of the stiffness matrix, and
!> after the last one invert it and write it out. The displacements and the
!> springs are taken to be the average values on the surface.
!------------------------------------------------------------------------------
  SUBROUTINE ModelLumpingSprings( Lump, Solver, Model, Permutation )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    INTEGER :: Permutation
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), LocalDisp(:,:), &
        xp(:), yp(:), zp(:)
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    REAL(KIND=dp), POINTER :: Displacement(:)
    REAL(KIND=dp), POINTER CONTIG :: PValues(:)
    REAL(KIND=dp) :: detJ, u, v, w, s, up, vp, wp, dFii, Dx, &
        ForceAtIp(3), MomentAtIp(3), Coord(3)
    INTEGER, POINTER :: Indexes(:), DisplPerm(:)
    INTEGER :: N_Integ, pn, n, i, j, k, t, dim, maxnodes, STDOFs
    LOGICAL :: stat, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: KmatFile

    SAVE Nodes
!------------------------------------------------------------------------------
    Mesh => Model % Mesh
    dim = CoordinateSystemDimension()
    maxnodes = Model % MaxElementNodes

    Displacement => Solver % Variable % Values
    DisplPerm => Solver % Variable % Perm
    STDOFs = Solver % Variable % DOFs

    ALLOCATE( Basis(maxnodes), dBasisdx(maxnodes,3), LocalDisp(dim,maxnodes), &
        xp(maxnodes), yp(maxnodes), zp(maxnodes) )

    ! Allocated once and kept: this is called six times a run, and the earlier
    ! version leaked a set of these arrays on every one of them.
    IF ( .NOT. ASSOCIATED( Lump % ParentNodes % x ) ) THEN
      ALLOCATE( Lump % ParentNodes % x(maxnodes), Lump % ParentNodes % y(maxnodes), &
          Lump % ParentNodes % z(maxnodes) )
    END IF

    dFii = 1.0d-2
    dX = 1.0d-2*SQRT(Lump % Area)

    IF (Permutation == 1) THEN
      Lump % Kmat = 0.0d0
      Lump % KmatMin = HUGE(Lump % KmatMin)
    END IF

    IF( Lump % FixDisplacement ) THEN
      IF(Permutation == 1) THEN
        n = SIZE( Displacement ) / STDOFs
        IF ( ALLOCATED( Lump % NodalLoads ) ) DEALLOCATE( Lump % NodalLoads )
        IF ( ALLOCATED( Lump % NodeVisited ) ) DEALLOCATE( Lump % NodeVisited )
        ALLOCATE( Lump % NodalLoads( STDOFs * n ), Lump % NodeVisited( n ) )
      END IF

      Lump % NodalLoads = 0.0d0
      PValues => Solver % Matrix % Values
      Solver % Matrix % Values => Solver % Matrix % BulkValues
      CALL MatrixVectorMultiply( Solver % Matrix, Displacement, Lump % NodalLoads)
      Solver % Matrix % Values => PValues

      Lump % NodeVisited = .FALSE.

      DO t = 1, Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        BC => GetBC()
        IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
        IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE

        n = GetElementNOFNodes()
        CALL GetElementNodes( Nodes )
        Indexes => Element % NodeIndexes

        DO i=1,n
          j = DisplPerm( Indexes(i) )
          IF(Lump % NodeVisited(j)) CYCLE
          Lump % NodeVisited(j) = .TRUE.

          Coord(1) = Nodes % x(i)
          Coord(2) = Nodes % y(i)
          Coord(3) = Nodes % z(i)
          Coord = Coord - Lump % Center

          DO k=1,dim
            ForceAtIP(k) = Lump % NodalLoads(3*(j-1)+k)
          END DO

          MomentAtIp(1) = -ForceAtIp(2) * Coord(3) + ForceAtIp(3) * Coord(2)
          MomentAtIp(2) = -ForceAtIp(3) * Coord(1) + ForceAtIp(1) * Coord(3)
          MomentAtIp(3) = -ForceAtIp(1) * Coord(2) + ForceAtIp(2) * Coord(1)

          Lump % Kmat(1:3,Permutation) = Lump % Kmat(1:3,Permutation) + ForceAtIp
          Lump % Kmat(4:6,Permutation) = Lump % Kmat(4:6,Permutation) + MomentAtIp
        END DO
      END DO


    ELSE
      DO t = 1, Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF ( .NOT. ActiveBoundaryElement() ) CYCLE

        BC => GetBC()
        IF ( .NOT.ASSOCIATED( BC ) ) CYCLE
        IF(.NOT. GetLogical( BC, 'Model Lumping Boundary',Found )) CYCLE

        n = GetElementNOFNodes()
        CALL GetElementNodes( Nodes )

        ! Get parent element & nodes:
        ! ---------------------------
        Parent => Element % BoundaryInfo % Left
        stat = ASSOCIATED( Parent )
        ! Guarded on stat, not on .NOT.stat: the test is only meaningful once
        ! Parent is known to exist, and the negated form dereferenced a null
        ! pointer in exactly the case it was meant to protect against, while
        ! never applying the permutation test when Left was there. The Right
        ! branch below has always had this the right way round.
        IF ( stat ) stat = ALL(DisplPerm(Parent % NodeIndexes) > 0)
        IF ( .NOT. stat ) THEN
          Parent => Element % BoundaryInfo % Right
          stat = ASSOCIATED( Parent )
          IF ( stat ) stat = ALL(DisplPerm(Parent % NodeIndexes) > 0)
          IF ( .NOT. stat ) CALL Fatal( 'ModelLumping', &
              'Cannot find proper parent for side element' )
        END IF
        pn = GetElementNOFNodes( Parent )
        CALL GetElementNodes( Lump % ParentNodes, Parent )
        CALL GetVectorLocalSolution( LocalDisp, UElement=Parent )

        ! Get boundary nodal points in parent local coordinates:
        ! ------------------------------------------------------
        DO i = 1,n
          DO j = 1,pn
            IF ( Element % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
              xp(i) = Parent % TYPE % NodeU(j)
              yp(i) = Parent % TYPE % NodeV(j)
              zp(i) = Parent % TYPE % NodeW(j)
              EXIT
            END IF
          END DO
        END DO

        IntegStuff = GaussPoints( Element )

        U_Integ => IntegStuff % u
        V_Integ => IntegStuff % v
        W_Integ => IntegStuff % w
        S_Integ => IntegStuff % s
        N_Integ =  IntegStuff % n

        DO k=1,N_Integ
          u = U_Integ(k)
          v = V_Integ(k)
          w = W_Integ(k)

          ! Basis function values & derivatives at the integration point:
          !--------------------------------------------------------------
          stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
              Basis, dBasisdx )

          s = detJ * S_Integ(k)
          IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
              CurrentCoordinateSystem() == CylindricSymmetric ) THEN
            s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
          END IF

          ! The plane  elements only include the  derivatives in the direction
          ! of the plane. Therefore compute the derivatives of the displacement
          ! field from the parent element:
          ! -------------------------------------------------------------------
          Up = SUM( xp(1:n) * Basis(1:n) )
          Vp = SUM( yp(1:n) * Basis(1:n) )
          Wp = SUM( zp(1:n) * Basis(1:n) )

          stat = ElementInfo( Parent, Lump % ParentNodes, Up, Vp, Wp, detJ, &
              Basis, dBasisdx )

          DO i=1,dim
            ForceAtIP(i) = SUM( Basis(1:pn) * LocalDisp(i,1:pn) )
          END DO

          MomentAtIP(1) = 0.5 * &
              ( SUM( dBasisdx(1:pn,2) * LocalDisp(3,1:pn)) &
              - SUM( dBasisdx(1:pn,3) * LocalDisp(2,1:pn)) )
          MomentAtIp(2) = 0.5 * &
              ( SUM( dBasisdx(1:pn,3) * LocalDisp(1,1:pn)) &
              - SUM( dBasisdx(1:pn,1) * LocalDisp(3,1:pn)) )
          MomentAtIp(3) = 0.5 * &
              ( SUM( dBasisdx(1:pn,1) * LocalDisp(2,1:pn)) &
              - SUM( dBasisdx(1:pn,2) * LocalDisp(1,1:pn)) )

          Lump % Kmat(Permutation,1:3) = Lump % Kmat(Permutation,1:3) + s * ForceAtIp
          Lump % Kmat(Permutation,4:6) = Lump % Kmat(Permutation,4:6) + s * MomentAtIp

          DO i = 1,dim
            IF(ABS(Lump % KmatMin(Permutation,i)) > ABS(ForceAtIp(i))) THEN
              Lump % KmatMin(Permutation,i) = ForceAtIp(i)
            END IF
            IF(ABS(Lump % KmatMin(Permutation,i+3)) > ABS(MomentAtIp(i))) THEN
              Lump % KmatMin(Permutation,i+3) = MomentAtIp(i)
            END IF
          END DO
        END DO
      END DO
    END IF



    IF(Permutation == 6) THEN
      KmatFile = LumpingFile( Solver )

      WRITE( Message, * ) 'Saving lumped elastic spring to file ', TRIM(KmatFile)
      CALL Info( 'ModelLumping', Message, Level=4 )

      IF (Lump % FixDisplacement) THEN
        Lump % Kmat(:,1:3) = Lump % Kmat(:,1:3) / dX
        Lump % Kmat(:,4:6) = Lump % Kmat(:,4:6) / dFii

        IF( ListGetLogical(Solver % Values,'Symmetrisize',stat)) THEN
          Lump % Kmat = (Lump % Kmat + TRANSPOSE(Lump % Kmat)) / 2.0d0
        END IF
      ELSE
        Lump % Kmat = Lump % Kmat / Lump % Area

        ! Save the Kmatrix prior to inversion to external file
        OPEN (10, FILE= TRIM(KmatFile) // ".inv")
        DO i=1,Permutation
          WRITE(10,'(6ES17.8E3)') Lump % Kmat(i,:)
        END DO
        CLOSE(10)

        OPEN (10, FILE= TRIM(KmatFile) // ".min-inv")
        DO i=1,Permutation
          WRITE(10,'(6ES17.8E3)') Lump % KmatMin(i,:)
        END DO
        CLOSE(10)

        IF(ListGetLogical(Solver % Values,'Symmetrisize',stat)) THEN
          Lump % Kmat = (Lump % Kmat + TRANSPOSE(Lump % Kmat)) / 2.0d0
          Lump % KmatMin = (Lump % KmatMin + TRANSPOSE(Lump % KmatMin)) / 2.0d0
        END IF

        CALL InvertMatrix(Lump % Kmat,Permutation)
        CALL InvertMatrix(Lump % KmatMin,Permutation)

        OPEN (10, FILE= TRIM(KmatFile) // ".min" )
        DO i=1,Permutation
          WRITE(10,'(6ES17.8E3)') Lump % KmatMin(i,:)
        END DO
        CLOSE(10)
      END IF

      ! Save the Kmatrix to an external file
      OPEN (10, FILE=KmatFile)
      DO i=1,Permutation
        WRITE(10,'(6ES17.8E3)') Lump % Kmat(i,:)
      END DO
      CLOSE(10)

      ! Save the area center to an external file
      OPEN (10, FILE= TRIM(KmatFile) // ".center")
      WRITE(10,'(3ES17.8E3)') Lump % Center
      CLOSE(10)
    END IF

    IF(Lump % FixDisplacement .AND. Permutation == 6) THEN
      DEALLOCATE( Lump % NodalLoads, Lump % NodeVisited )
    END IF

    DEALLOCATE( Basis, dBasisdx, LocalDisp, xp, yp, zp )
!------------------------------------------------------------------------------
  END SUBROUTINE ModelLumpingSprings
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Generalized cartesian lumped mass matrix: the total mass and the moments of
!> inertia of the body about the lumping centre, written out beside the springs.
!------------------------------------------------------------------------------
  SUBROUTINE CartesianMass( Lump, Solver, Model )
!------------------------------------------------------------------------------
    TYPE(ModelLumping_t) :: Lump
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), Density(:)
    REAL(KIND=dp) :: vol, SqrtElementMetric, Dens
    REAL(KIND=dp) :: x, y, z, U, V, W, S
    REAL(KIND=dp) :: Moment0, Moment1(3), Moment2(3,3), Center(3), MassMatrix(6,6)
    INTEGER :: i, n, t, mat_id, body_id, maxnodes
    LOGICAL :: stat
    CHARACTER(LEN=MAX_NAME_LEN) :: KmatFile

    SAVE Nodes
!------------------------------------------------------------------------------
    Mesh => Model % Mesh
    maxnodes = Mesh % MaxElementNodes
    ALLOCATE( Basis(maxnodes), dBasisdx(maxnodes,3), Density(maxnodes) )

!------------------------------------------------------------------------------
! Do some initialization stuff
!------------------------------------------------------------------------------

    vol = 0.0d0
    Moment0 = 0.0d0
    Moment1 = 0.0d0
    Moment2 = 0.0d0
    Center = Lump % Center

!------------------------------------------------------------------------------
! Integrate the lumped mass over the volume/area
!------------------------------------------------------------------------------

    DO t = 1, Solver % NumberOfActiveElements
      Element => Mesh % Elements( Solver % ActiveElements( t ) )
      Model % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      CALL GetElementNodes( Nodes, Element )

      body_id = Element % BodyId
      mat_id = ListGetInteger( Model % Bodies( body_id ) % Values, &
          'Material', minv=1,maxv=Model % NumberOfMaterials )
      Material => Model % Materials(mat_id) % Values
      Density(1:n) = ListGetReal( Material, 'Density', n, NodeIndexes(1:n) )

      IntegStuff = GaussPoints( Element )

      DO i=1,IntegStuff % n

        U = IntegStuff % u(i)
        V = IntegStuff % v(i)
        W = IntegStuff % w(i)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,U,V,W,&
            SqrtElementMetric,Basis,dBasisdx)
!------------------------------------------------------------------------------
!      Coordinatesystem dependent info
!------------------------------------------------------------------------------
        s = SqrtElementMetric * IntegStuff % s(i)
        x = SUM(Nodes % x(1:n) * Basis(1:n)) - Center(1)
        y = SUM(Nodes % y(1:n) * Basis(1:n)) - Center(2)
        z = SUM(Nodes % z(1:n) * Basis(1:n)) - Center(3)

        IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
          s = 2.0 * PI * x * s
        END IF
        vol =  vol + S
        dens = SUM(Basis(1:n) * Density(1:n) )

        Moment0 = Moment0 + s * dens

        Moment1(1) = Moment1(1) + s * x * dens
        Moment1(2) = Moment1(2) + s * y * dens
        Moment1(3) = Moment1(3) + s * z * dens

        Moment2(1,1) = Moment2(1,1) + s * ( y*y + z*z)  * dens
        Moment2(2,2) = Moment2(2,2) + s * ( x*x + z*z )  * dens
        Moment2(3,3) = Moment2(3,3) + s * ( x*x + y*y ) * dens

        Moment2(1,2) = Moment2(1,2) - s * x * y * dens
        Moment2(1,3) = Moment2(1,3) - s * x * z * dens
        Moment2(2,3) = Moment2(2,3) - s * y * z * dens
      END DO
    END DO

    IF(Vol < AEPS) THEN
      DEALLOCATE( Basis, dBasisdx, Density )
      RETURN
    END IF

    Moment2(2,1) = Moment2(1,2)
    Moment2(3,1) = Moment2(1,3)
    ! Symmetrising the other way round, as the two lines above do. Only the
    ! upper triangle is accumulated, so the old form copied the never-written
    ! (3,2) onto the computed (2,3) -- zeroing it, which made the reported
    ! "res: Moment of inertia YZ" and the mass matrix's yz coupling always 0.
    Moment2(3,2) = Moment2(2,3)

    CALL ListAddConstReal(Model % Simulation,'res: Mass',Moment0)

    CALL ListAddConstReal(Model % Simulation,'res: Lumped Center X',Center(1))
    CALL ListAddConstReal(Model % Simulation,'res: Lumped Center Y',Center(2))
    CALL ListAddConstReal(Model % Simulation,'res: Lumped Center Z',Center(3))

    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XX',Moment2(1,1))
    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia YY',Moment2(2,2))
    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia ZZ',Moment2(3,3))
    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XY',Moment2(1,2))
    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia XZ',Moment2(1,3))
    CALL ListAddConstReal(Model % Simulation,'res: Moment of inertia YZ',Moment2(2,3))

    MassMatrix = 0.0d0
    DO i= 1,3
      MassMatrix(i,i) = Moment0
    END DO
    MassMatrix(4:6,4:6) = Moment2

    ! Save the mass matrix to an external file
    KmatFile = LumpingFile( Solver )
    OPEN (10, FILE= TRIM(KmatFile) // ".mass")
    DO i=1,6
      WRITE(10,'(6ES17.8E3)') MassMatrix(i,:)
    END DO
    CLOSE(10)

    DEALLOCATE( Basis, dBasisdx, Density )
!------------------------------------------------------------------------------
  END SUBROUTINE CartesianMass
!------------------------------------------------------------------------------

END MODULE ModelLumping
!------------------------------------------------------------------------------
