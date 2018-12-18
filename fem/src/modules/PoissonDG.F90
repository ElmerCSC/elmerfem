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

     TYPE(Element_t), POINTER :: Element, &
       ParentElement, LeftParent, RightParent

     LOGICAL :: AllocationsDone = .FALSE., GotIt, Stat
     INTEGER :: k, n, np, nr, nl, t, istat, i, j, Indexes(128)
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, st
#else
     REAL(KIND=dp) :: at, st, CPUTime
#endif
     REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), &
              FORCE(:), EpsilonBoundary(:)
     REAL(KIND=dp) :: Beta, Gamma, Norm
     
     SAVE STIFF, LOAD, FORCE, AllocationsDone, EpsilonBoundary
!*******************************************************************************

     TYPE( Element_t ), POINTER :: Edges(:), Faces(:), p(:), elm
     TYPE(Nodes_t) :: ElementNodes

     TYPE(Mesh_t), POINTER :: Mesh

     LOGICAL :: found

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

    Beta = ListGetConstReal( Solver % Values, 'beta', GotIt )
    IF( .NOT. GotIt ) Beta = Gamma

    IF(CoordinateSystemDimension() == 2 ) THEN
      np = Mesh % NumberOfEdges
      p => Edges
    ELSE
      np = Mesh % NumberOfFaces
      p => Faces
    END IF

    FORCE = 0.0d0
    DO t = 1, np
       elm => p(t)
       IF ( .NOT. ActiveBoundaryElement(elm, DGBoundary=.TRUE.) ) CYCLE
       
       LeftParent  => elm % BoundaryInfo % Left
       RightParent => elm  % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED(RightParent) ) CYCLE

       n = GetElementNOFnodes(elm)
       nl = GetElementNOFDOFs(LeftParent)
       nr = GetElementNOFDOFs(RightParent)

       BC => GetBC(elm)
       IF (  ASSOCIATED(BC) .AND. GetLogical( BC, 'Discont BC', Found) ) THEN
         CALL LocalJumpsDiscontBC( STIFF, Elm, n, LeftParent, nl, RightParent, nr, beta )
       ELSE
         CALL LocalJumps( STIFF, Elm, n, LeftParent, nl, RightParent, nr, gamma )
       END IF

       CALL DefaultUpdateEquations( STIFF, FORCE, elm )
    END DO
!------------------------------------------------------------------------------
!   Loop over the boundary elements (Nitsche boundary stab. terms)
!------------------------------------------------------------------------------
    Gamma = ListGetConstReal( Solver % Values, 'gamma2', GotIt )
    IF( .NOT. GotIt ) Gamma = 1.0d-3
    
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
!    PRINT*,'Assembly (s): ',CPUTime()-at

     st = CPUTime()
     Norm = DefaultSolve()
!    PRINT*,'Solve (s):    ',CPUTime()-st

   CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW( Nodes, n, &
         ParentNodes, nParent, U, V, W, Basis )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE( Nodes_t ) :: Nodes, ParentNodes
      INTEGER :: n, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i, j, Check
      REAL(KIND=dp) :: Dist, DistTolerance
      REAL(KIND=dp) :: NodalParentU(n), &
           NodalParentV(n), NodalParentW(n)
!------------------------------------------------------------------------------
      DistTolerance = 1.0d-12

      Check = 0
      DO i = 1,n
         DO j = 1,nParent
            Dist = (Nodes % x(i) - ParentNodes % x(j))**2 & 
                 + (Nodes % y(i) - ParentNodes % y(j))**2 & 
                 + (Nodes % z(i) - ParentNodes % z(j))**2

            IF( Dist < DistTolerance ) THEN
               Check = Check+1
               NodalParentU(i) = RightParent % Type % NodeU(j)
               NodalParentV(i) = RightParent % Type % NodeV(j)
               NodalParentW(i) = RightParent % Type % NodeW(j)
            END IF

         END DO
      END DO
      IF( Check /= n ) STOP 'Error' 

      U = SUM( Basis(1:n) * NodalParentU(1:n) )
      V = SUM( Basis(1:n) * NodalParentV(1:n) )
      W = SUM( Basis(1:n) * NodalParentW(1:n) )
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
     SUBROUTINE LocalJumps( STIFF,Elm,n,LeftParent,nl,RightParent,nr,gamma )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), gamma
       INTEGER :: n,nl,nr
       TYPE(Element_t), POINTER :: Elm, LeftParent, RightParent

!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
       REAL(KIND=dp) :: LeftBasis(nl), LeftdBasisdx(nl,3)
       REAL(KIND=dp) :: RightBasis(nr), RightdBasisdx(nr,3)
       REAL(KIND=dp) :: LeftdBasisdn(nl), RightdBasisdn(nr)
       REAL(KIND=dp) :: Jump(nl+nr), AverageFlux(nl+nr)
       REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S
       LOGICAL :: Stat
       INTEGER :: i, j, p, q, dim, t
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       REAL(KIND=dp) :: hE, Normal(3), LeftOut(3)

       TYPE(Nodes_t) ::Nodes, LeftParentNodes, RightParentNodes
       SAVE Nodes, LeftParentNodes, RightParentNodes 
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       STIFF = 0.0d0

       CALL GetElementNodes( Nodes, elm )
       CALL GetElementNodes( LeftParentNodes, LeftParent )
       CALL GetElementNodes( RightParentNodes, RightParent )
       hE = ElementDiameter( elm, Nodes )

       LeftOut(1) = SUM(LeftParentNodes % x(1:nl)) / nl
       LeftOut(2) = SUM(LeftParentNodes % y(1:nl)) / nl
       LeftOut(3) = SUM(LeftParentNodes % z(1:nl)) / nl

       LeftOut(1) = SUM(Nodes % x(1:n)) / n - LeftOut(1)
       LeftOut(2) = SUM(Nodes % y(1:n)) / n - LeftOut(2)
       LeftOut(3) = SUM(Nodes % z(1:n)) / n - LeftOut(3)

!------------------------------------------------------------------------------
!      Numerical integration over the edge
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints(elm)
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( elm, Nodes, U, V, W, SqrtElementMetric, &
              Basis, dBasisdx )

         S = S * SqrtElementMetric

         Normal = NormalVector( Elm, Nodes, U, V, .FALSE. )
         IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

!        Find basis functions for the parent elements:
!        ----------------------------------------------
         CALL FindParentUVW( Nodes, n, LeftParentNodes, &
              nl, U, V, W, Basis )

         stat = ElementInfo( LeftParent, LeftParentNodes, &
              U, V, W, SqrtElementMetric, LeftBasis, LeftdBasisdx )

         CALL FindParentUVW( Nodes, n, RightParentNodes, &
              nr, U, V, W, Basis )

         stat = ElementInfo( RightParent, RightParentNodes, &
              U, V, W, SqrtElementMetric, RightBasis, RightdBasisdx )

!        Integrate jump terms:
!        ---------------------
         Jump(1:nl) = -LeftBasis(1:nl)
         Jump(nl+1:nl+nr) = RightBasis(1:nr)

         DO i = 1,nl
            LeftdBasisdn(i)  = SUM( LeftdBasisdx(i,:)  * Normal(:) )
         END DO

         DO i = 1,nr
            RightdBasisdn(i) = SUM( RightdBasisdx(i,:) * Normal(:) )
         END DO

         AverageFlux(1:nl) = LeftdBasisdn(1:nl) / 2.0d0
         AverageFlux(nl+1:nl+nr) = RightdBasisdn(1:nr) / 2.0d0

         DO p = 1,nl+nr
            DO q = 1,nl+nr
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
     SUBROUTINE LocalJumpsDiscontBC( STIFF,Elm,n,LeftParent,nl,RightParent,nr,gamma )
!------------------------------------------------------------------------------
       REAL(KIND=dp) :: STIFF(:,:), gamma
       INTEGER :: n,nl,nr
       TYPE(Element_t), POINTER :: Elm, LeftParent, RightParent

!------------------------------------------------------------------------------
       REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
       REAL(KIND=dp) :: LeftBasis(nl), LeftdBasisdx(nl,3)
       REAL(KIND=dp) :: RightBasis(nr), RightdBasisdx(nr,3)
       REAL(KIND=dp) :: LeftdBasisdn(nl), RightdBasisdn(nr)
       REAL(KIND=dp) :: Jump(nl+nr), AverageFlux(nl+nr)
       REAL(KIND=dp) :: SqrtElementMetric, U, V, W, S
       LOGICAL :: Stat
       INTEGER :: i, j, p, q, dim, t
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       REAL(KIND=dp) :: hE, Normal(3), LeftOut(3)

       TYPE(Nodes_t) ::Nodes, LeftParentNodes, RightParentNodes
       SAVE Nodes, LeftParentNodes, RightParentNodes 
!------------------------------------------------------------------------------
       dim = CoordinateSystemDimension()
       STIFF = 0.0d0

       CALL GetElementNodes( Nodes, elm )
       CALL GetElementNodes( LeftParentNodes, LeftParent )
       CALL GetElementNodes( RightParentNodes, RightParent )
       hE = ElementDiameter( elm, Nodes )

       LeftOut(1) = SUM(LeftParentNodes % x(1:nl)) / nl
       LeftOut(2) = SUM(LeftParentNodes % y(1:nl)) / nl
       LeftOut(3) = SUM(LeftParentNodes % z(1:nl)) / nl

       LeftOut(1) = SUM(Nodes % x(1:n)) / n - LeftOut(1)
       LeftOut(2) = SUM(Nodes % y(1:n)) / n - LeftOut(2)
       LeftOut(3) = SUM(Nodes % z(1:n)) / n - LeftOut(3)

!------------------------------------------------------------------------------
!      Numerical integration over the edge
!------------------------------------------------------------------------------
       IntegStuff = GaussPoints(elm)
 
       DO t=1,IntegStuff % n
         U = IntegStuff % u(t)
         V = IntegStuff % v(t)
         W = IntegStuff % w(t)
         S = IntegStuff % s(t)

!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( elm, Nodes, U, V, W, SqrtElementMetric, &
              Basis, dBasisdx )

         S = S * SqrtElementMetric

         Normal = NormalVector( Elm, Nodes, U, V, .FALSE. )
         IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

!        Find basis functions for the parent elements:
!        ----------------------------------------------
         CALL FindParentUVW( Nodes, n, LeftParentNodes, &
              nl, U, V, W, Basis )

         stat = ElementInfo( LeftParent, LeftParentNodes, &
              U, V, W, SqrtElementMetric, LeftBasis, LeftdBasisdx )

         CALL FindParentUVW( Nodes, n, RightParentNodes, &
              nr, U, V, W, Basis )

         stat = ElementInfo( RightParent, RightParentNodes, &
              U, V, W, SqrtElementMetric, RightBasis, RightdBasisdx )

!        Integrate jump terms:
!        ---------------------
         Jump(1:nl) = -LeftBasis(1:nl)
         Jump(nl+1:nl+nr) = RightBasis(1:nr)

         DO i = 1,nl
            LeftdBasisdn(i)  = SUM( LeftdBasisdx(i,:)  * Normal(:) )
         END DO

         DO i = 1,nr
            RightdBasisdn(i) = SUM( RightdBasisdx(i,:) * Normal(:) )
         END DO

         AverageFlux(1:nl) = LeftdBasisdn(1:nl) / 2.0d0
         AverageFlux(nl+1:nl+nr) = RightdBasisdn(1:nr) / 2.0d0

         DO p = 1,nl+nr
            DO q = 1,nl+nr
               STIFF(p,q) = STIFF(p,q) + (gamma/hE)*Jump(p)*Jump(q) * s
               STIFF(p,q) = STIFF(p,q) - AverageFlux(p) * Jump(q)   * s
               STIFF(p,q) = STIFF(p,q) + Jump(p) * AverageFlux(q)   * s
            END DO
         END DO
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE LocalJumpsDisContBC
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
