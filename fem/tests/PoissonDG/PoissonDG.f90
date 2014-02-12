!/******************************************************************************
! *
! *       ELMER, A Computational Fluid Dynamics Program.
! *
! *       Copyright 1st April 1995 - , CSC - IT Center for Science Ltd.,
! *                                    Finland.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! ******************************************************************************
! *
! *                    Author:       Juha Ruokolainen
! *
! *                    Address: CSC - IT Center for Science Ltd.
! *                                Keilaranta 14, P.O. BOX 405
! *                                  02101 Espoo, Finland
! *                                  Tel. +358 0 457 2723
! *                                Telefax: +358 0 457 2302
! *                              EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by: Mikko Lyly
! *
! *       Date of modification: 02 Dec 2003
! *
! *****************************************************************************/
 
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: PoissonSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solves the Poisson equation with discontinous Galerkin method!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  DOUBLE PRECISION :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver
     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: BC, BodyForce

     TYPE(Element_t), POINTER :: Element, Edge, &
       ParentElement, LeftParent, RightParent

     LOGICAL :: AllocationsDone = .FALSE., GotIt, Stat
     INTEGER :: k, n, t, istat, i, j, Indexes(128)
     REAL(KIND=dp) :: Norm, at, st, CPUTime, Gamma
     REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), &
              FORCE(:), EpsilonBoundary(:)

     SAVE STIFF, LOAD, FORCE, AllocationsDone, EpsilonBoundary
!*******************************************************************************

     TYPE( Element_t ), POINTER :: Edges(:), Faces(:)
     TYPE(Nodes_t) :: ElementNodes

     TYPE(Mesh_t), POINTER :: Mesh

     REAL(KIND=dp), POINTER :: LocalSolution(:)
     INTEGER :: Active, DOFs
!*******************************************************************************

     Mesh => GetMesh()

     Edges => Mesh % Edges
     Faces => Mesh % Faces
!------------------------------------------------------------------------------
!    Initialize & allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone ) THEN
        N = 2 * Mesh % MaxElementDOFs
        ALLOCATE( FORCE(N), STIFF(2*N,2*N), LOAD(N), EpsilonBoundary(N), &
             LocalSolution(N), STAT = istat )
        
       IF ( istat /= 0 ) CALL FATAL('PoissonDG','Memory allocation error.' )

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------
!    Assembly of the bulk elements
!------------------------------------------------------------------------------
     at = CPUTime()
     Active = GetNOFActive()
     CALL DefaultInitialize()
     DO t = 1, Active
        Element => GetActiveElement( t ) 
        n = GetElementNOfNodes( Element )
        
        BodyForce => GetBodyForce( Element )
        LOAD(1:n) = GetReal( BodyForce, 'Source', GotIt )
        
        CALL LocalMatrix( STIFF, FORCE, LOAD, Element, n ) 
        CALL DefaultUpdateEquations( STIFF, FORCE )
     END DO
!------------------------------------------------------------------------------
!    Assembly of the edge terms
!------------------------------------------------------------------------------
    Gamma = ListGetConstReal( Solver % Values, 'gamma1', GotIt )
    IF( .NOT. GotIt ) Gamma = 1.0d3

    FORCE = 0.0d0
    DO t = 1, Mesh % NumberOfEdges
       Edge => Edges(t)
       IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
       LeftParent  => Edge % BoundaryInfo % Left
       RightParent => Edge % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED(RightParent) ) CYCLE

       n = GetElementNOFnodes( Edge )
       CALL LocalJumps( STIFF, Edge, n, LeftParent, RightParent, gamma )

       CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
    END DO
!------------------------------------------------------------------------------
!   Loop over the boundary elements (Nitsche boundary stab. terms)
!------------------------------------------------------------------------------
    Gamma = ListGetConstReal( Solver % Values, 'gamma2', GotIt )
    IF( .NOT. GotIt ) Gamma = 1.0d-3
    
    DO t = 1, Mesh % NumberOfBoundaryElements
       Element => GetBoundaryElement(t)
       IF( .NOT. ActiveBoundaryElement() )  CYCLE
       IF( GetElementFamily(Element) == 1 ) CYCLE
       
       BC => GetBC()
       IF ( ASSOCIATED(BC) ) THEN
          ParentElement => Element % BoundaryInfo % Left
          IF ( .NOT. ASSOCIATED( ParentElement ) ) &
             ParentElement => Element % BoundaryInfo % Right
          
          n = GetElementNOFNodes( Element )
          k = GetElementNOFnodes( ParentElement )

          LOAD(1:n) = GetReal( BC, 'g', GotIt )
          EpsilonBoundary(1:n) = GetReal( BC, 'e', GotIt )
          
          CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
            Element, ParentElement, n, k, EpsilonBoundary, Gamma )
          
          CALL DefaultUpdateEquations( STIFF, FORCE )
       END IF
     END DO
     CALL DefaultFinishAssembly()
     PRINT*,'Assembly (s): ',CPUTime()-at

     st = CPUTime()
     Norm = DefaultSolve()
     PRINT*,'Solve (s):    ',CPUTime()-st
!------------------------------------------------------------------------------
!    Write the post file (so far works only with triangles in 2d):
!------------------------------------------------------------------------------
     PRINT *,'INFO: Writing post file "dg.ep"'
     OPEN( 10, FILE='dg.ep', STATUS='unknown' )
     
     WRITE( 10, '(i8,i8,i3,i6,a)' ) &
          3*Active, Active, 1, 1, ' scalar: Somequantity'
     
     DO t=1,Active
        Element => GetActiveElement(t) 
        n = GetElementNOFNodes( Element )
        
        DO i = 1,n
           WRITE( 10, '(3e20.12)') &
                Mesh % Nodes % x(Element % NodeIndexes(i)), &
                Mesh % Nodes % y(Element % NodeIndexes(i)), &
                Mesh % Nodes % z(Element % NodeIndexes(i))
        END DO
     END DO

     DO i=1,Active
        WRITE( 10,'(a,i6)',advance='no' ) 'body ', &
             Mesh % Elements(i) % Type % ElementCode
  
        n = Mesh % Elements(i) % Type % NumberOfNodes
        DO j=1,n
           WRITE( 10,'(i8)', ADVANCE='no' ) n*(i-1) + j - 1
        END DO
        WRITE( 10,'(a)' ) ''
     END DO
     
     WRITE( 10,* ) '#time ',1,1,1

     DO i=1,ACtive
       Element => GetActiveElement(i)
       n = GetElementDOFs( Indexes )
       DO k=1,n
          DO j=1, Solver % Variable % DOFs
             WRITE( 10, '(1e30.15)') Solver % Variable % Values( &
               Solver % Variable % DOFs*( &
                   Solver % Variable % Perm(Indexes(k))-1)+j )
          END DO
       END DO
     END DO

     CLOSE(10)

   CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW( EdgeNodes, nEdge, &
         ParentNodes, nParent, U, V, W, EdgeBasis )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE( Nodes_t ) :: EdgeNodes, ParentNodes
      INTEGER :: nEdge, nParent
      REAL( KIND=dp ) :: U, V, W, EdgeBasis(:)
!------------------------------------------------------------------------------
      INTEGER :: i, j, Check
      REAL(KIND=dp) :: Dist, DistTolerance
      REAL(KIND=dp) :: NodalParentU(nEdge), &
           NodalParentV(nEdge), NodalParentW(nEdge)
!------------------------------------------------------------------------------
      DistTolerance = 1.0d-12

      Check = 0
      DO i = 1,nEdge
         DO j = 1,nParent
            Dist = (EdgeNodes % x(i) - ParentNodes % x(j))**2 & 
                 + (EdgeNodes % y(i) - ParentNodes % y(j))**2 & 
                 + (EdgeNodes % z(i) - ParentNodes % z(j))**2

            IF( Dist < DistTolerance ) THEN
               Check = Check+1
               NodalParentU(i) = RightParent % Type % NodeU(j)
               NodalParentV(i) = RightParent % Type % NodeV(j)
               NodalParentW(i) = RightParent % Type % NodeW(j)
            END IF

         END DO
      END DO
      IF( Check /= nEdge ) STOP 'Error' 

      U = SUM( EdgeBasis(1:nEdge) * NodalParentU(1:nEdge) )
      V = SUM( EdgeBasis(1:nEdge) * NodalParentV(1:nEdge) )
      W = SUM( EdgeBasis(1:nEdge) * NodalParentW(1:nEdge) )
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
     SUBROUTINE LocalJumps( STIFF,Edge,n,LeftParent,RightParent,gamma )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), gamma
       INTEGER :: n
       TYPE(Element_t), POINTER :: Edge, LeftParent, RightParent

!------------------------------------------------------------------------------
       REAL(KIND=dp) :: EdgeBasis(n), EdgedBasisdx(n,3)
       REAL(KIND=dp) :: LeftBasis(n+1), LeftdBasisdx(n+1,3)
       REAL(KIND=dp) :: RightBasis(n+1), RightdBasisdx(n+1,3)
       REAL(KIND=dp) :: LeftdBasisdn(n+1), RightdBasisdn(n+1)
       REAL(KIND=dp) :: Jump(2*(n+1)), AverageFlux(2*(n+1))
       REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S
       LOGICAL :: Stat
       INTEGER :: i, j, p, q, dim, t, nEdge, nParent
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       REAL(KIND=dp) :: hE, Normal(3), LeftOut(3)

       TYPE(Nodes_t) :: EdgeNodes, LeftParentNodes, RightParentNodes
       SAVE EdgeNodes, LeftParentNodes, RightParentNodes 
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       nEdge = n
       nParent = n+1
       STIFF = 0.0d0

       CALL GetElementNodes( EdgeNodes, Edge )
       CALL GetElementNodes( LeftParentNodes, LeftParent )
       CALL GetElementNodes( RightParentNodes, RightParent )
       hE = ElementDiameter( Edge, EdgeNodes )

       LeftOut(1) = SUM(LeftParentNodes % x(1:nParent)) / nParent
       LeftOut(2) = SUM(LeftParentNodes % y(1:nParent)) / nParent
       LeftOut(3) = SUM(LeftParentNodes % z(1:nParent)) / nParent
       LeftOut(1) = SUM(EdgeNodes % x(1:n)) / n - LeftOut(1)
       LeftOut(2) = SUM(EdgeNodes % y(1:n)) / n - LeftOut(2)
       LeftOut(3) = SUM(EdgeNodes % z(1:n)) / n - LeftOut(3)

!------------------------------------------------------------------------------
!      Numerical integration over the edge
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Edge )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Edge, EdgeNodes, U, V, W, SqrtElementMetric, &
              EdgeBasis, EdgedBasisdx )

         S = S * SqrtElementMetric

         Normal = NormalVector( Edge, EdgeNodes, U, V, .FALSE. )
         IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

!        Find basis functions for the parent elements:
!        ----------------------------------------------
         CALL FindParentUVW( EdgeNodes, nEdge, LeftParentNodes, &
              nParent, U, V, W, EdgeBasis )

         stat = ElementInfo( LeftParent, LeftParentNodes, &
              U, V, W, SqrtElementMetric, LeftBasis, LeftdBasisdx )

         CALL FindParentUVW( EdgeNodes, nEdge, RightParentNodes, &
              nParent, U, V, W, EdgeBasis )

         stat = ElementInfo( RightParent, RightParentNodes, &
              U, V, W, SqrtElementMetric, RightBasis, RightdBasisdx )

!        Integrate jump terms:
!        ---------------------
         Jump(1:nParent) = LeftBasis(1:nParent)
         Jump(nParent+1:2*nParent) = -RightBasis(1:nParent)

         DO i = 1,nParent
            LeftdBasisdn(i)  = SUM( LeftdBasisdx(i,:)  * Normal(:) )
            RightdBasisdn(i) = SUM( RightdBasisdx(i,:) * Normal(:) )
         END DO

         AverageFlux(1:nParent) = LeftdBasisdn(1:nParent) / 2.0d0
         AverageFlux(nParent+1:2*nParent) = RightdBasisdn(1:nParent) / 2.0d0

         DO p = 1,2*nParent
            DO q = 1,2*nParent
               STIFF(p,q) = STIFF(p,q) + (gamma/hE)*Jump(p)*Jump(q) * s
               STIFF(p,q) = STIFF(p,q) - AverageFlux(p) * Jump(q)   * s
               STIFF(p,q) = STIFF(p,q) + Jump(p) * AverageFlux(q)   * s
            END DO
         END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------      
     SUBROUTINE LocalMatrix( STIFF, FORCE, LOAD, Element, n  )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
       INTEGER :: n
       TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S,A,L
       LOGICAL :: Stat

       INTEGER :: i,p,q,t,dim

       TYPE(Nodes_t) :: Nodes
       SAVE Nodes
 
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
 
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       FORCE = 0.0d0
       STIFF = 0.0d0
       CALL GetElementNodes( Nodes, Element )
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                    Basis, dBasisdx )

         S = S * SqrtElementMetric

         L = SUM( LOAD(1:n) * Basis(1:n) )

!------------------------------------------------------------------------------
!        The Poisson equation
!------------------------------------------------------------------------------
         STIFF(1:n,1:n) = STIFF(1:n,1:n) +  s * &
                MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
         FORCE(1:n) = FORCE(1:n) + s*L*Basis(1:n)
!------------------------------------------------------------------------------
      END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix

!==============================================================================

!------------------------------------------------------------------------------
     SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
          Element, ParentElement, n, k, EpsilonBoundary, Gamma )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), &
            LOAD(:), EpsilonBoundary(:), Gamma
       INTEGER :: n, k
       TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
       REAL(KIND=dp) :: ParentBasis(k), ParentdBasisdx(k,3), &
             ParentdBasisdn(k)
       INTEGER :: i,j,p,q,t,dim

       REAL(KIND=dp) :: Normal(3), ParentU, ParentV, ParentW
       REAL(kind=dp) :: e, g, h, c1, c2, c3, L

       REAL(KIND=dp) :: ParentNodalU(n), parentNodalV(n), ParentNodalW(n)

       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       REAL(KIND=dp) :: SqrtElementMetric,U,V,W,S
       LOGICAL :: Stat

       TYPE(Nodes_t) :: Nodes, ParentNodes
       SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       FORCE = 0.0d0
       STIFF = 0.0d0

       CALL GetElementNodes( Nodes, Element )
       CALL GetElementNodes( ParentNodes, ParentElement )
!------------------------------------------------------------------------------
!      Numerical integration
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints( Element )
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

         Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                    Basis,dBasisdx )
         S = S * SqrtElementMetric

         ParentNodalU = 0.0d0
         ParentNodalV = 0.0d0
         ParentNodalW = 0.0d0
         DO i = 1,n
            DO j = 1,k
               IF( Element % NodeIndexes(i) == ParentElement % NodeIndexes(j) ) THEN
                  ParentNodalU(i) = ParentElement % Type % NodeU(j)
                  ParentNodalV(i) = ParentElement % Type % NodeV(j)
                  ParentNodalW(i) = ParentElement % Type % NodeW(j)
               END IF
            END DO
         END DO
         ParentU = SUM( Basis(1:n) * ParentNodalU(1:n) )
         ParentV = SUM( Basis(1:n) * ParentNodalV(1:n) )
         ParentW = SUM( Basis(1:n) * ParentNodalW(1:n) )

         stat = ElementInfo( ParentElement, ParentNodes, &
              ParentU, ParentV, ParentW, SqrtElementMetric, &
              ParentBasis, ParentdBasisdx )

         L = SUM( LOAD(1:n) * Basis(1:n) )
         e = SUM( EpsilonBoundary(1:n) * Basis(1:n) )

         DO p = 1,k
            ParentdBasisdn(p) = SUM( Normal(1:3) * ParentdBasisdx(p,1:3) )
         END DO
!------------------------------------------------------------------------------
!        The stabilization terms
!------------------------------------------------------------------------------
         h = ElementDiameter( Element, Nodes )
         c1 = Gamma*h/(e + Gamma*h)
         c2 = 1.0d0/(e + Gamma*h)
         c3 = (e * Gamma*h)/(e + Gamma*h)
         DO p = 1,k
           DO q = 1,k
              STIFF(p,q) = STIFF(p,q) - c1 * ParentdBasisdn(p) * ParentBasis(q) * s
              STIFF(p,q) = STIFF(p,q) - c1 * ParentBasis(p) * ParentdBasisdn(q) * s
              STIFF(p,q) = STIFF(p,q) + c2 * ParentBasis(p) * ParentBasis(q) * s
              STIFF(p,q) = STIFF(p,q) - c3 * ParentdBasisdn(p) * ParentdBasisdn(q) * s
           END DO
           FORCE(p) = FORCE(p) - c1 * ParentdBasisdn(p) * L * s
           FORCE(p) = FORCE(p) + c2 * ParentBasis(p) * L * s
         END DO
!------------------------------------------------------------------------------
       END DO
!------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------
   END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------
