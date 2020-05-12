!
! IMPORTANT NOTE: The example implementation given here is obsolete and should 
! not be considered as a model for implementing solvers based on edge finite 
! elements. For the vector Helmholtz solver see instead the solver module 
!
!   .../fem/src/modules/VectorHelmholtz.F90
!
SUBROUTINE HrotSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the vector Helmholtz equation
! 
!   -rot rot u + A u = f
!
!  using edge elements (Nedelec/Whitney basis of lowest degree).
!
!  Dirichlet b.c. is given for u.t. Robin b.c. given in the form
!
!     rot u = g - B u.t
!
!  The variariational form discretized by the fem is precisely
!
!   (rot u, rot v) + (A u, v) + <B u.t, v.t> = (f, v) + <g, v.t> 
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element, Edge, Parent
  TYPE(Nodes_t) :: Nodes, EdgeNodes
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: Var

  REAL(KIND=dp) :: Norm, Normal(3), Normal2(3), &
       u, v, w, NormalSign, Edgeh

  INTEGER :: n, nb, nd, t, istat, i, j, k, l, active

  REAL(KIND=dp),    ALLOCATABLE :: LOAD(:,:), Acoef(:,:), Bcoef(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), FORCE(:)

  TYPE(Mesh_t), POINTER :: Mesh
  LOGICAL :: stat, EigenAnalysis

  SAVE STIFF, LOAD, MASS, FORCE, Acoef, Bcoef, &
       AllocationsDone, Nodes, EdgeNodes
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(2*N), LOAD(2,N), STIFF(2*N,2*N), &
          MASS(2*N,2*N), Acoef(2,N), Bcoef(2,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  Solver % Matrix % Complex = .TRUE.
  EigenAnalysis = GetLogical( GetSolverParams(), 'Eigen Analysis' )
  
  !System assembly:
  !----------------
  active = GetNOFActive()
  CALL DefaultInitialize()
  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes() ! kulmat
     nd = GetElementNOFDOFs()  ! vapausasteet

     LOAD = 0.0d0
     BodyForce => GetBodyForce()
     IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1,1:n) = GetReal( BodyForce, 'f 1', Found )
        Load(2,1:n) = GetReal( BodyForce, 'f 2', Found )
     END IF

     Acoef = 0.0d0
     Material => GetMaterial( Element )
     IF ( ASSOCIATED(material) ) THEN
        Acoef(1,1:n) = GetReal( Material, 'A 1', Found )
        Acoef(2,1:n) = GetReal( Material, 'A 2', Found )
     END IF

     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix(  STIFF, MASS, FORCE, LOAD, Acoef, Element, n, nd )

     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )
     IF ( EigenAnalysis ) CALL DefaultUpdateMass( MASS )
  END DO

!------------------------------------------------------------------------------
!
!  Neumann & Newton BCs:
!  ---------------------
   Var => Solver % Variable

   DO t=1, Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF ( GetElementFamily() == 1 .OR. .NOT.ActiveBoundaryElement() ) CYCLE
      
      BC => GetBC(Element)
      IF ( .NOT. ASSOCIATED(BC) ) CYCLE
         
      n = GetElementNOFNodes(Element)
      
      Load(1,1:n) = GetReal( BC, 'g 1', Found )
      Load(2,1:n) = GetReal( BC, 'g 2', Found )
      
      Bcoef(1,1:n) = GetReal( BC, 'B 1', Found )
      Bcoef(2,1:n) = GetReal( BC, 'B 2', Found )
         
!     The boundary integral terms are computed for the parent element:
!     ----------------------------------------------------------------
      Parent => Element % BoundaryInfo % Left
      IF ( .NOT. ASSOCIATED( Parent ) ) THEN
         Parent => Element % BoundaryInfo % Right
      END IF
      IF ( .NOT. ASSOCIATED( Parent ) ) CYCLE

      STIFF = 0.0d0
      FORCE = 0.0d0
      DO j=1,Parent % Type % NumberOfEdges
         Edge => Mesh % Edges( Parent % EdgeIndexes(j) )

         n = 0
         DO k=1,Element % Type % NumberOfNodes
            DO l=1,Edge % Type % NumberOfNodes
               IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) n = n + 1
            END DO
         END DO

         IF ( n == Edge % Type % NumberOfNodes ) THEN
            Edgeh = (Model % Mesh % Nodes % x(Edge % NodeIndexes(2)) &
                    -Model % Mesh % Nodes % x(Edge % NodeIndexes(1)))**2 &
                   +(Model % Mesh % Nodes % y(Edge % NodeIndexes(2)) &
                    -Model % Mesh % Nodes % y(Edge % NodeIndexes(1)))**2 
            Edgeh = SQRT( Edgeh )

            FORCE(j) = DCMPLX( SUM( Load(1,1:2))/2.0d0, &
                              SUM( Load(2,1:2))/2.0d0 )
            
            STIFF(j,j) =  DCMPLX( SUM(Bcoef(1,1:2))/2.0d0 / Edgeh, &
                                 SUM(Bcoef(2,1:2))/2.0d0 / Edgeh )
         END IF
      END DO
      CALL DefaultUpdateEquations( STIFF, FORCE, Parent )
   END DO

!------------------------------------------------------------------------------
   CALL DefaultFinishAssembly()
   CALL MyDirichletBCs()

   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()

!*****************************************************************************
   CALL WriteResults
!*****************************************************************************

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, MASS, FORCE, LOAD, Acoef, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: LOAD(:,:), Acoef(:,:)
    COMPLEX(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,L1,L2
    REAL(KIND=dp) :: WhitneyBasis(3,2), RotWhitneyBasis(3), A1, A2
    LOGICAL :: Stat
    INTEGER :: t, i, j, k
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )


    STIFF = 0.0d0
    FORCE = 0.0d0
    MASS = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
!      Haetaan ensin normaalit solmukantafunktiot:
!      -------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t),  detJ, Basis, dBasisdx )

       L1 = SUM( Basis(1:n) * LOAD(1,1:n) )
       L2 = SUM( Basis(1:n) * LOAD(2,1:n) )

       A1 = SUM( Basis(1:n) * Acoef(1,1:n) )
       A2 = SUM( Basis(1:n) * Acoef(2,1:n) )

!      Haetaan Raviart-Thomas kantafunktiot:
!      -------------------------------------
       CALL EdgeBasis( WhitneyBasis, RotWhitneyBasis, &
            IP % u(t), IP % v(t), IP % w(t), Element, n )

!      Compute element stiffness matrix and force vector:
!      ---------------------------------------------------
       DO i = 1,3
         FORCE(i) = FORCE(i) + (L1*WhitneyBasis(i,1) + L2*WhitneyBasis(i,2)) * detJ * IP % s(t)
         DO j = 1,3
           STIFF(i,j) = STIFF(i,j) + RotWhitneyBasis(i) * RotWhitneyBasis(j) * detJ * IP % s(t)
           STIFF(i,j) = STIFF(i,j) &
                + DCMPLX( A1, A2 ) * SUM(WhitneyBasis(i,:)*WhitneyBasis(j,:)) * detJ * IP % s(t)
           MASS(i,j) = MASS(i,j) +  SUM(WhitneyBasis(i,:)*WhitneyBasis(j,:)) * detJ * IP % s(t)
         END DO
       END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

  SUBROUTINE WriteResults
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: u, v, w, Ratkaisu(2), Ratkaisu2(2), WhitneyBasis(3,2), &
         RotWhitneyBasis(3), Edgeh, Rot, Rot2

    INTEGER :: Active, Indexes(6)

    OPEN(unit=10, file='ratkaisu.ep')

   Active = GetNOFActive()
   t = Active * 3 ! nurkkasolmujen lukumäärä

   write(10,'(4i10,A)') t, Active, 8, 1, &
        ' vector: ReU vector: ImU scalar: ReRotU scalar: ImRotU'

!  Kirjoitetaan solmukoordinaatit:
!  -------------------------------
   Active = GetNOFActive()
   DO i = 1, Active
      Element => GetActiveElement(i)
      DO j = 1,3
         k = Element % NodeIndexes(j)
         write(10,*) Model % Mesh % Nodes % x(k), &
              Model % Mesh % Nodes % y(k), Model % Mesh % Nodes % z(k)
      END DO
   END DO

!  Kirjoitetaan elementit:
!  -----------------------
   WRITE(10,'(a)') '#group all'
   DO i = 1, Active
      Element => GetActiveElement(i)
      WRITE(10,*) 'body1   303', 3*(i-1)+0,3*(i-1)+1,3*(i-1)+2
   END DO
   WRITE(10,'(a)') '#endgroup all'

!  Kirjoitetaan tulokset:
!  ----------------------
   Var => VariableGet( Model % Variables, 'tangent component' )

   WRITE(10,'(a)') '#time 1   1  1  1.0'
   DO i = 1, Active

      Element => GetActiveElement(i)

!     Luuppi nurkkapisteiden yli:
!     ---------------------------
      w = 0.0d0
      DO j = 1,3

         ! Nurkan koordinaatit:
         ! --------------------
         SELECT CASE(j)
         CASE(1)
            u = 0.0d0
            v = 0.0d0
         CASE(2)
            u = 1.0d0
            v = 0.0d0
         CASE(3)
            u = 0.0d0
            v = 1.0d0
         END SELECT

         ! Kantafunktioiden arvot nurkkapisteessä:
         ! ---------------------------------------
         CALL EdgeBasis( WhitneyBasis, RotWhitneyBasis, u, v, w, Element, n )
         k = GetElementDOFs( Indexes )

         ratkaisu = 0.0d0
         ratkaisu2 = 0.0d0
         Rot = 0.0d0
         Rot2 = 0.0d0
         DO k = 1,3
           ratkaisu(1) = ratkaisu(1) + WhitneyBasis(k,1) * Var % Values(2*Var % Perm( Indexes(k))-1)
           ratkaisu(2) = ratkaisu(2) + WhitneyBasis(k,2) * Var % Values(2*Var % Perm( Indexes(k))-1)

           ratkaisu2(1) = ratkaisu2(1) + WhitneyBasis(k,1) * Var % Values(2*Var % Perm( Indexes(k)))
           ratkaisu2(2) = ratkaisu2(2) + WhitneyBasis(k,2) * Var % Values(2*Var % Perm( Indexes(k)))
            
           Rot = Rot + RotWhitneyBasis(k) * Var % Values(2*Var % Perm(Indexes(k))-1 )
           Rot2 = Rot2 + RotWhitneyBasis(k) * Var % Values(2*Var % Perm(Indexes(k)) )
         END DO
         WRITE(10,'(6E20.10)') Ratkaisu(1), Ratkaisu(2), 0.0, Ratkaisu2(1), Ratkaisu2(2), 0.0, Rot, Rot2
      END DO
   END DO

   CLOSE(10)
!------------------------------------------------------------------------
 END SUBROUTINE WriteResults
!------------------------------------------------------------------------


!------------------------------------------------------------------------
 SUBROUTINE MyDirichletBCs( USolver )
!------------------------------------------------------------------------
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver

     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Solver_t), POINTER :: Solver
     REAL(KIND=dp), POINTER    :: b(:)
     INTEGER :: i,j, k, l, n, nb, DOF
     REAL(KIND=dp) :: xx
     LOGICAL :: Flag,Found
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement

     CHARACTER(LEN=MAX_NAME_LEN) :: name

!    For 2d elements:
!    ----------------
     REAL(KIND=dp) :: BoundaryData(2), Edgeh

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => Model % Solver
     END IF

     A => Solver % Matrix
     x => Solver % Variable
     b => A % RHS

     DO DOF=1,x % DOFs
        name = x % name
        IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)
        CALL SetDirichletBoundaries( Model, A, b, Name, DOF, x % DOFs, x % Perm )

        SaveElement => Model % CurrentElement

!       Dirichlet BCs for face & edge DOFs:
!       -----------------------------------
        DO i=1,Mesh % NumberOfBoundaryElements

           Element => GetBoundaryElement(i)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           BC => GetBC()
           IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
           IF ( .NOT. ListCheckPresent( BC, TRIM(Name) ) ) CYCLE

           BoundaryData(1:2) = GetReal( BC, TRIM(name), Found )

           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) THEN
               Parent => Element % BoundaryInfo % Right
           END IF
           IF ( .NOT. ASSOCIATED( Parent ) ) CYCLE

           IF ( ASSOCIATED( Mesh % Edges ) ) THEN
              DO j=1,Parent % Type % NumberOfEdges
                 Edge => Mesh % Edges( Parent % EdgeIndexes(j) )

                 Edgeh = (Model % Mesh % Nodes % x(Edge % NodeIndexes(2)) &
                         -Model % Mesh % Nodes % x(Edge % NodeIndexes(1)))**2 &
                        +(Model % Mesh % Nodes % y(Edge % NodeIndexes(2)) &
                         -Model % Mesh % Nodes % y(Edge % NodeIndexes(1)))**2 
                 Edgeh = SQRT( Edgeh )

                 n = 0
                 DO k=1,Element % Type % NumberOfNodes
                    DO l=1,Edge % Type % NumberOfNodes
                       IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) n = n + 1
                    END DO
                 END DO

                 IF ( n == Edge % Type % NumberOfNodes ) THEN

                    DO k=1,Edge % BDOFs

                       n = Mesh % NumberofNodes + &
                           (Parent % EdgeIndexes(j)-1) * Mesh % MaxEdgeDOFs + k

                       n = x % Perm( n )
                       IF ( n <= 0 ) CYCLE

                          DO l=1,x % DOFs
                          nb = x % DOFs * (n-1) + l
                          CALL CRS_ZeroRow( A, nb )

!                         Check if this could be done any better:
!                         ---------------------------------------
                          IF( l == 1 .AND. TRIM(name) == 'tangent component 1' ) THEN
                             IF( Edge % NodeIndexes(2) > Edge % NodeIndexes(1) ) THEN
                                A % RHS(nb) = +SUM( BoundaryData(1:2) ) / 2.0d0 * Edgeh
                             ELSE
                                A % RHS(nb) = -SUM( BoundaryData(1:2) ) / 2.0d0 * Edgeh
                             END IF
                          END IF

                          IF( l == 2 .AND. TRIM(name) == 'tangent component 2' ) THEN
                             IF( Edge % NodeIndexes(2) > Edge % NodeIndexes(1) ) THEN
                                A % RHS(nb) = +SUM( BoundaryData(1:2) ) / 2.0d0 * Edgeh
                             ELSE
                                A % RHS(nb) = -SUM( BoundaryData(1:2) ) / 2.0d0 * Edgeh
                             END IF
                          END IF

                          A % Values( A % Diag(nb) ) = 1.0d0
                       END DO
                    END DO
                 END IF
              END DO
           END IF

           IF ( ASSOCIATED( Mesh % Faces ) ) THEN
              DO j=1,Parent % Type % NumberOfFaces
                 Face => Mesh % Faces( Parent % FaceIndexes(j) )
                 n = 0
                 DO k=1,Element % Type % NumberOfNodes
                    DO l=1,Face % Type % NumberOfNodes
                       IF ( Face % NodeIndexes(l) == Element % NodeIndexes(k) ) n = n + 1
                    END DO
                 END DO

                 IF ( n == Face % Type % NumberOfNodes ) THEN
                    DO k=1,Face % BDOFs
                       n = Mesh % Numberofnodes +  &
                         Mesh % MaxEdgeDOFs * Mesh % NumberOfEdges + &
                            Mesh % MaxFaceDOFs * (Parent % FaceIndexes(j)-1) + k

                       n = x % Perm( n )
                       IF ( n <= 0 ) CYCLE
                       DO l=1,x % DOFs
                          nb = x % DOFs * (n-1) + l
                          CALL CRS_ZeroRow( A, nb )
                          A % RHS(nb) = 0.0d0
                          A % Values( A % Diag(nb) ) = 1.0d0
                       END DO
                    END DO
                    EXIT
                 END IF
              END DO
           END IF
        END DO
        Model % CurrentElement => SaveElement
     END DO
!------------------------------------------------------------------------
   END SUBROUTINE MyDirichletBCs
!------------------------------------------------------------------------


!------------------------------------------------------------------------
   SUBROUTINE EdgeBasis( WhitneyBasis, RotWhitneyBasis, u, v, w, Element, n )
!------------------------------------------------------------------------
     REAL(KIND=dp) :: WhitneyBasis(:,:), RotWhitneyBasis(:), u, v, w
     TYPE(Element_t) :: Element

     INTEGER :: n
!------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes
     SAVE Nodes

     REAL(KIND=dp) :: detJ, Basis(3), dBasisdx(3,3)
     LOGICAL :: stat
!------------------------------------------------------------------------
     CALL GetElementNodes( Nodes )
     
     stat = ElementInfo( Element, Nodes, u, v, w, &
          detJ, Basis, dBasisdx )

     WhitneyBasis(1,1) = Basis(1) * dBasisdx(2,1) - Basis(2) * dBasisdx(1,1)
     WhitneyBasis(1,2) = Basis(1) * dBasisdx(2,2) - Basis(2) * dBasisdx(1,2)
     WhitneyBasis(2,1) = Basis(2) * dBasisdx(3,1) - Basis(3) * dBasisdx(2,1)
     WhitneyBasis(2,2) = Basis(2) * dBasisdx(3,2) - Basis(3) * dBasisdx(2,2)
     WhitneyBasis(3,1) = Basis(3) * dBasisdx(1,1) - Basis(1) * dBasisdx(3,1)
     WhitneyBasis(3,2) = Basis(3) * dBasisdx(1,2) - Basis(1) * dBasisdx(3,2)

     RotWhitneyBasis(1) = 2.0d0 * ( dBasisdx(1,2) * dBasisdx(2,1) - dBasisdx(2,2) * dBasisdx(1,1) )
     RotWhitneyBasis(2) = 2.0d0 * ( dBasisdx(2,2) * dBasisdx(3,1) - dBasisdx(3,2) * dBasisdx(2,1) )
     RotWhitneyBasis(3) = 2.0d0 * ( dBasisdx(3,2) * dBasisdx(1,1) - dBasisdx(1,2) * dBasisdx(3,1) )

     IF( Element % NodeIndexes(2) < Element % NodeIndexes(1) ) THEN
        WhitneyBasis(1,:)  = -WhitneyBasis(1,:)
        RotWhitneyBasis(1) = -RotWhitneyBasis(1)
     END IF
     
     IF( Element % NodeIndexes(3) < Element % NodeIndexes(2) ) THEN
        WhitneyBasis(2,:)  = -WhitneyBasis(2,:)
        RotWhitneyBasis(2) = -RotWhitneyBasis(2)
     END IF
     
     IF( Element % NodeIndexes(1) < Element % NodeIndexes(3) ) THEN
        WhitneyBasis(3,:)  = -WhitneyBasis(3,:)
        RotWhitneyBasis(3) = -RotWhitneyBasis(3)
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE EdgeBasis
!------------------------------------------------------------------------------
 END SUBROUTINE HrotSolver
!------------------------------------------------------------------------------
