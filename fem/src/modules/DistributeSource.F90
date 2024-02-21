!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
!/******************************************************************************
! *
! *  Module for computing source current over a nonconforming mesh.
! *
! *  Authors: Peter Råback, Juhani Kataja
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: ~2015
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Solver for computing sources given in local mesh to global mesh.
!> Computes both the source and the fixing potential source related to it.
!-------------------------------------------------------------------------------
SUBROUTINE DistributeSource(Model, Solver, dt, Transient)
!-------------------------------------------------------------------------------
  USE DefUtils  
  IMPLICIT NONE
!-------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!-------------------------------------------------------------------------------
  LOGICAL :: Found
  TYPE(Variable_t), POINTER :: FixSol
  TYPE(Mesh_t), POINTER :: MeshLocal, MeshGlobal
  TYPE(Element_t), POINTER :: Element
!-------------------------------------------------------------------------------

!  INTERFACE
!    SUBROUTINE InterpolateMeshToMeshQ( OldMesh, NewMesh, OldVariables, &
!        NewVariables, UseQuadrantTree, Projector, MaskName, FoundNodes )
!      USE Types
!      TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
!      TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
!      LOGICAL, OPTIONAL :: UseQuadrantTree,FoundNodes(:)
!      CHARACTER(LEN=*),OPTIONAL :: MaskName
!      TYPE(Projector_t), POINTER, OPTIONAL :: Projector
!    END SUBROUTINE InterpolateMeshToMeshQ
!  END INTERFACE


  MeshLocal => GetMesh()
  MeshGlobal => Model % Meshes
  
  IF( ASSOCIATED( MeshLocal, MeshGlobal ) ) THEN
    CALL Fatal('DistributeSource','Meshes are the same!')
  END IF

  CALL InterpolateAndIntegrateSource( Model, Solver, MeshLocal, MeshGlobal, 'Local Source','g', .TRUE. )

CONTAINS
  
  !------------------------------------------------------------------------------
  !>    Interpolates and integrates functions of one mesh in another mesh.
  !------------------------------------------------------------------------------
  SUBROUTINE InterpolateAndIntegrateSource( Model, Solver, LocalMesh, GlobalMesh, &
      LocalSourceName, GlobalRhsName, UseQuadrantTree )
    !------------------------------------------------------------------------------
    USE DefUtils
    !-------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t), TARGET  :: LocalMesh   
    TYPE(Mesh_t), TARGET  :: GlobalMesh  
    CHARACTER(LEN=*) :: LocalSourceName
    CHARACTER(LEN=*) :: GlobalRhsName
    LOGICAL, OPTIONAL :: UseQuadrantTree  
    !------------------------------------------------------------------------------
    INTEGER :: dim
    TYPE(Nodes_t) :: ElementNodes, ElementNodesB
    INTEGER :: i, j, k, l, n, m, nb, ne, NoFails, NoFound, IndexesB(100)
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    TYPE(Variable_t), POINTER :: RhsSol
    TYPE(Quadrant_t), POINTER :: LeafQuadrant
    TYPE(Element_t),POINTER :: Element, ElementB
    TYPE(ValueList_t), POINTER :: BodyForce
    REAL(KIND=dp), ALLOCATABLE :: Basis(:),BasisB(:),RotWBasis(:,:), &
        WBasis(:,:), dBasisdx(:,:), dBAsisdxb(:,:)
    REAL(KIND=dp), ALLOCATABLE :: LocalSource(:), LocalCurrDens(:,:)
    REAL(KIND=dp) :: BoundingBox(6), detJ, detJb,u,v,w,ub,vb,wb,CurrDensAtip(3)
    LOGICAL :: UseQTree, TryQTree, Stat, EdgeBasis, PiolaT    
    TYPE(Quadrant_t), POINTER :: RootQuadrant
    INTEGER, POINTER :: NodeIndexes(:), NodeIndexesB(:)
    INTEGER :: EdgeBasisDegree 
    LOGICAL :: Found
    INTEGER :: eps_tries
    REAL(KIND=dp) :: eps1 = 0.1, eps2, eps_global, eps_local, eps_numeric
    INTEGER, ALLOCATABLE :: RhsInd(:), FixInd(:)
    TYPE(Solver_t), POINTER :: GSolver
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(*), PARAMETER :: Caller = 'IntegrateSource'
        
    
    SAVE :: ElementNodes, ElementNodesB 

    !------------------------------------------------------------------------------
    CALL Info(Caller,'Projecting current sources between meshes')

    Params => GetSolverParams(Solver)
    
    ! Check if using the spatial division hierarchy for the search:
    ! -------------------------------------------------------------
    dim = CoordinateSystemDimension()
    
    IF ( .NOT. PRESENT( UseQuadrantTree ) ) THEN
      UseQTree = .TRUE.
    ELSE
      UseQTree = UseQuadrantTree
    ENDIF

    IF ( UseQTree ) THEN
      RootQuadrant => GlobalMesh % RootQuadrant
      IF ( .NOT.ASSOCIATED( RootQuadrant ) ) THEN
        BoundingBox(1) = MINVAL(GlobalMesh % Nodes % x)
        BoundingBox(2) = MINVAL(GlobalMesh % Nodes % y)
        BoundingBox(3) = MINVAL(GlobalMesh % Nodes % z)
        BoundingBox(4) = MAXVAL(GlobalMesh % Nodes % x)
        BoundingBox(5) = MAXVAL(GlobalMesh % Nodes % y)
        BoundingBox(6) = MAXVAL(GlobalMesh % Nodes % z)
        
        eps2 = 0.1_dp * MAXVAL(BoundingBox(4:6)-BoundingBox(1:3))
        BoundingBox(1:3) = BoundingBox(1:3) - eps2
        BoundingBox(4:6) = BoundingBox(4:6) + eps2
        
        CALL BuildQuadrantTree( GlobalMesh,BoundingBox,GlobalMesh % RootQuadrant)
        RootQuadrant => GlobalMesh % RootQuadrant

        CALL Info(Caller,'Quadrant tree build for interpolation',Level=10)
      END IF
    END IF
    
    TryQTree = ASSOCIATED(RootQuadrant) .AND. UseQTree 
       
    !------------------------------------------------------------------------------
    n = MAX( GlobalMesh % MaxElementDOFs, LocalMesh % MaxElementDOFs )
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    ALLOCATE( ElementNodesB % x(n), ElementNodesB % y(n), ElementNodesB % z(n) )
    ALLOCATE( Basis(n), BasisB(n), LocalSource(n), LocalCurrDens(3,n), FixInd(n), &
             RhsInd(n), WBasis(n,3), RotWBasis(n,3), dBasisdx(n,3), dBasisdxb(n,3) )
    
    eps_global = ListGetConstReal( Model % Simulation,  &
        'Interpolation Global Epsilon', Stat)
    IF(.NOT. Stat) eps_global = 2.0d-10

    eps_local = ListGetConstReal( Model % Simulation,  &
        'Interpolation Local Epsilon', Stat )
    IF(.NOT. Stat) eps_local = 1.0d-10
    
    eps_tries = ListGetInteger( Model % Simulation,  &
        'Interpolation Max Iterations', Stat )
    IF(.NOT. Stat) eps_tries = 12

    eps_numeric = ListGetConstReal( Model % Simulation, &
        'Interpolation Numeric Epsilon', Stat)
    IF(.NOT. Stat) eps_numeric = 1.0e-10

    NoFails = 0
    NoFound = 0

        
    EdgeBasis = GetLogical( Params,'Edge Basis', Stat )

    RhsSol => VariableGet( GlobalMesh % Variables,GlobalRhsName, .TRUE. )
    IF(.NOT. ASSOCIATED( RhsSol) ) THEN
      CALL Fatal(Caller,'Could not find variable: '//TRIM(GlobalRhsName))
    END IF
    DO i=1,Model % NumberOfSolvers
      IF( ASSOCIATED( RhsSol % Solver, Model % Solvers(i) ) ) THEN
        CALL Info(Caller,'Source will be projected to solver id: '//I2S(i),Level=8)
        GSolver => Model % Solvers(i)
        EXIT
      END IF
    END DO
    IF(.NOT. ASSOCIATED( GSolver ) ) THEN
      CALL Fatal(Caller,'Could not associate r.h.s. vector to any solver!')
    END IF
    
    EdgeBasis = .FALSE.
    PiolaT = .FALSE.
    n = GlobalMesh % NumberOfNodes
    m = SIZE( RhsSol % Perm )
    IF( m > n ) EdgeBasis = ANY( RhsSol % Perm(n+1:m) > 0 )
    IF( EdgeBasis ) THEN
      CALL Info(Caller,'Projecting source to edge basis functions')
      FixSol => VariableGet( GlobalMesh % Variables, GlobalRhsName // ' fix', .TRUE. )
      IF(.NOT. ASSOCIATED(FixSol) ) THEN
        CALL Warn(Caller,'Could not find variable: '//TRIM(GlobalRhsName)//' fix')
      END IF
      CALL EdgeElementStyle(GSolver % Values, PiolaT, BasisDegree = EdgeBasisDegree ) 
    ELSE
      CALL Fatal(Caller,'Currently we assume edge basis!')
    END IF

    ! Initialize the source term and fixing potential source term
    RhsSol % Values = 0.0_dp
    IF(ASSOCIATED(FixSol)) FixSol % Values = 0.0_dp

    ! Create the nodal coordinates for all Gaussian integration points
    !-----------------------------------------------------------------    
    CALL Info(Caller,'Integrating over bulk elements of local mesh',Level=8)

    DO i=1, LocalMesh % NumberOfBulkElements
      Element => LocalMesh % Elements(i)
      
      Model % CurrentElement => Element
      n = Element % TYPE % NumberOfNodes        
      NodeIndexes => Element % NodeIndexes

      ElementNodes % x(1:n) = LocalMesh % Nodes % x(NodeIndexes(1:n))
      ElementNodes % y(1:n) = LocalMesh % Nodes % y(NodeIndexes(1:n))
      ElementNodes % z(1:n) = LocalMesh % Nodes % z(NodeIndexes(1:n))
      
      BodyForce => GetBodyForce(Element)
      IF(.NOT. ASSOCIATED( BodyForce ) ) CYCLE
      
      CALL GetRealVector( BodyForce, LocalCurrDens(1:3,1:n), 'Local Current Density',Found)
      IF(.NOT.Found) CYCLE

      IF( EdgeBasis ) THEN
        IntegStuff = GaussPoints( Element, EdgeBasis = EdgeBasis, PReferenceElement=PiolaT,&
            EdgeBasisDegree=EdgeBasisDegree )
      ELSE
        IntegStuff = GaussPoints( Element )
      END IF

      DO j=1,IntegStuff % n
        
        u = IntegStuff % u(j)
        v = IntegStuff % v(j)
        w = IntegStuff % w(j)

        Stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, Basis, dBasisdx )

        Point(1) = SUM( Basis(1:n) * ElementNodes % x(1:n) )
        Point(2) = SUM( Basis(1:n) * ElementNodes % y(1:n) )
        Point(3) = SUM( Basis(1:n) * ElementNodes % z(1:n) )

        CurrDensAtIp(1) = SUM( Basis(1:n) * LocalCurrDens(1,1:n) )
        CurrDensAtIp(2) = SUM( Basis(1:n) * LocalCurrDens(2,1:n) )
        CurrDensAtIp(3) = SUM( Basis(1:n) * LocalCurrDens(3,1:n) )

        !PRINT *,'CurrentDensity',CurrDensAtIp,'at coordinate',Point
        
        !------------------------------------------------------------------------------
        ! Find in which old mesh bulk element the point belongs to
        !------------------------------------------------------------------------------
        Found = .FALSE.

        IF( TryQTree ) THEN
          !------------------------------------------------------------------------------
          ! Find the last existing quadrant that the point belongs to
          ! Using octree search.
          !------------------------------------------------------------------------------
          Found = .FALSE.
          CALL FindLeafElements(Point, dim, RootQuadrant, LeafQuadrant)
          
          IF ( ASSOCIATED(LeafQuadrant) ) THEN
            ! Go through the bulk elements in the last ChildQuadrant
            ! only.  Try to find matching element with progressively
            ! sloppier tests. Allow at most 100 % of slack:
            ! -------------------------------------------------------
            Eps1 = eps_global
            Eps2 = eps_local

            DO k=1, LeafQuadrant % NElemsInQuadrant
              ElementB => GlobalMesh % Elements(LeafQuadrant % Elements(k))

              nb = ElementB % TYPE % NumberOfNodes
              NodeIndexesB  => ElementB % NodeIndexes

              ne = Elementb % Type % NumberOfEdges
              IndexesB(1:ne) = Elementb % EdgeIndexes + GlobalMesh % NumberOfNodes

              ElementNodesB % x(1:nB) = GlobalMesh % Nodes % x(NodeIndexesB)
              ElementNodesB % y(1:nB) = GlobalMesh % Nodes % y(NodeIndexesB)
              ElementNodesB % z(1:nB) = GlobalMesh % Nodes % z(NodeIndexesB)

              Found = PointInElement( ElementB, ElementNodesB, Point, LocalCoordinates )
              IF ( Found ) EXIT
            END DO
          ELSE
            ! Dummy N^2 search
            DO k=1,GlobalMesh % NumberOfBulkElements
              ElementB => GlobalMesh % Elements(k)

              nb = ElementB % TYPE % NumberOfNodes
              NodeIndexesB  => ElementB % NodeIndexes

              ne = Elementb % Type % NumberOfEdges
              IndexesB(1:ne) = Elementb % EdgeIndexes + GlobalMesh % NumberOfNodes

              ElementNodesB % x(1:nb) = GlobalMesh % Nodes % x(NodeIndexesB)
              ElementNodesB % y(1:nb) = GlobalMesh % Nodes % y(NodeIndexesB)
              ElementNodesB % z(1:nb) = GlobalMesh % Nodes % z(NodeIndexesB)

              Found = PointInElement( ElementB, ElementNodesB, Point, LocalCoordinates  ) 
              IF( Found ) EXIT
            END DO
          END IF
        END IF

        IF(.NOT. Found ) THEN
          CALL Info( Caller,'Point '//I2s(i)//' was not found in any of the elements!',Level=15)
          NoFails = NoFails + 1
          CYCLE
        END IF

        RhsInd(1:ne) = RhsSol % Perm( IndexesB(1:ne) )
        !IF( ANY( RhsInd(1:ne) == 0 ) .OR. ANY( FixInd(1:nb) == 0) ) CYCLE

        NoFound = NoFound + 1

        vb = LocalCoordinates(1)
        ub = LocalCoordinates(2)
        wb = LocalCoordinates(3)

        ! Evaluate basis functions at integration point
        stat = ElementInfo(ElementB, ElementNodesB, ub, vb, wb, detJb, Basisb, dBasisdxb,&
            EdgeBasis = Wbasis, RotBasis = RotWBasis, USolver = GSolver ) 

        ! Note that we use the metric of the local mesh but populate it on the global mesh
        IF( ALL( RhsInd(1:ne) > 0 ) ) THEN
          RhsSol % Values( RhsInd(1:ne) ) = RhsSol % Values( RhsInd(1:ne) ) + &
              MATMUL( WBasis(1:ne,:), CurrDensAtIP ) * detJ * IntegStuff % s(j)
        END IF
          
        IF( ASSOCIATED( FixSol ) ) THEN
          FixInd(1:nb) = FixSol % Perm( NodeIndexesB )
          IF( ALL( FixInd(1:nb) > 0 ) ) THEN
            FixSol % Values( FixInd(1:nb) ) = FixSol % Values( FixInd(1:nb) ) + &
                MATMUL( dBasisdxb(1:nb,:), CurrDensAtIp ) * detJ * IntegStuff % s(j)
          END IF
        END IF
      END DO
    END DO
    
    CALL Info(Caller,'Integration finished',Level=15)
    
    IF( NoFails == 0 ) THEN
      CALL Info(Caller,'Found all nodes '//I2S(NoFound)//&
          ' in the target mesh',Level=6)
    ELSE
      CALL Warn(Caller,'Points not found: '//I2S(NoFails)//' (found '&
          //I2S(NoFound)//')')
    END IF

        
    PRINT *,'Global source sum:',SUM( RhsSol % Values ), SUM( ABS( RhsSol % Values ) ) 
    IF( ASSOCIATED( FixSol ) ) THEN
      PRINT *,'Global fixing sum:',SUM( FixSol % Values ), SUM( ABS( FixSol % Values ) )
    END IF
      
  END SUBROUTINE InterpolateAndIntegrateSource


END SUBROUTINE DistributeSource
