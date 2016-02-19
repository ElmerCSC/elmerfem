SUBROUTINE Wsolve( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation!
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
  LOGICAL :: AllocationsDone = .FALSE., Found, PosEl, NegEl
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, istat, active, i, j
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, BC, BodyParams, Material
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
                                RotM(:,:,:), Tcoef(:,:,:)
  REAL (KIND=DP), POINTER :: Cwrk(:,:,:)
  CHARACTER(LEN=MAX_NAME_LEN):: CoilType
  LOGICAL :: CoilBody


  SAVE STIFF, LOAD, FORCE, Tcoef, RotM, Cwrk, AllocationsDone
!------------------------------------------------------------------------------

  IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), Tcoef(3,3,N), RotM(3,3,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'Wsolve', 'Memory allocation error.' )
     END IF

     NULLIFY( Cwrk )

     AllocationsDone = .TRUE.
  END IF

   !System assembly:
   !----------------
   Active = GetNOFActive()
   CALL DefaultInitialize()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      LOAD = 0.0d0
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) &
         Load(1:n) = GetReal( BodyForce, 'Source', Found )

      BodyParams => GetBodyParams( Element )
      IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('Wsolve', 'Body Parameters not found')

      CoilBody = .False.
      CoilType = GetString(BodyParams, 'Coil Type', Found)
      IF (.NOT. Found) THEN
        CoilType = ''
      ELSE
        SELECT CASE (CoilType)
        CASE ('stranded')
          CoilBody = .True.
        CASE ('massive')
          CoilBody = .True.
        CASE ('foil winding')
          CoilBody = .True.
          CALL GetElementRotM(Element, RotM, n)
        CASE DEFAULT
          CALL Fatal ('Poisson', 'Non existent Coil Type Chosen!')
        END SELECT
      END IF

      Tcoef = 0.0d0
      Material => GetMaterial( Element )
      IF ( ASSOCIATED(Material) ) THEN
!------------------------------------------------------------------------------
!      Read conductivity values (might be a tensor)
!------------------------------------------------------------------------------

        CALL ListGetRealArray( Material, &
               'Electric Conductivity', Cwrk, n, Element % NodeIndexes, Found )

        IF (Found) THEN
           IF ( SIZE(Cwrk,1) == 1 ) THEN
              DO i=1,3
                 Tcoef( i,i,1:n ) = Cwrk( 1,1,1:n )
              END DO
           ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
              DO i=1,MIN(3,SIZE(Cwrk,1))
                 Tcoef(i,i,1:n) = Cwrk(i,1,1:n)
              END DO
           ELSE
              DO i=1,MIN(3,SIZE(Cwrk,1))
                 DO j=1,MIN(3,SIZE(Cwrk,2))
                    Tcoef( i,j,1:n ) = Cwrk(i,j,1:n)
                 END DO
              END DO
           END IF
        END IF
      END IF
 
      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, CoilType, Tcoef, RotM, n, nd+nb )
      CALL LCondensate( nd, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   DO t=1, Solver % Mesh % NumberOfBoundaryElements
     ! get element and BC info
     ! -----------------------
     Element => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement() ) CYCLE
     n = GetElementNOFNodes()
     ! no evaluation of Neumann BC’s on points
     IF ( GetElementFamily() == 1 ) CYCLE
     BC => GetBC()
     FORCE = 0.0d00
     STIFF = 0.0d00

     PosEl = .False.
     NegEl = .False.
     PosEl = GetLogical(BC,'Positive Electrode', Found)
     NegEl = GetLogical(BC,'Negative Electrode', Found)

     LOAD = 0._dp
     IF (PosEl) LOAD = 1._dp
     IF (NegEl) LOAD = -1._dp
     
     IF (Solver % Variable % name == 'w') CALL BoundaryCondition(LOAD, FORCE, Element, n)

     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   CALL DefaultFinishAssembly()
   CALL DefaultDirichletBCs()

   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, CoilType, Tcoef, RotM, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:), Tcoef(:,:,:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP, C(3,3), &
                     RotMLoc(3,3), RotM(3,3,n)
    CHARACTER(LEN=MAX_NAME_LEN):: CoilType
    LOGICAL :: Stat
    INTEGER :: i,j,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    C(1,1:3) = [1,0,0]
    C(2,1:3) = [0,1,0]
    C(3,1:3) = [0,0,1]

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      ! Compute the conductivity tensor
      ! -------------------------------
      DO i=1,3
        DO j=1,3
          C(i,j) = SUM( Tcoef(i,j,1:n) * Basis(1:n) )
          IF(CoilType == 'foil winding') RotMLoc(i,j) = SUM( RotM(i,j,1:n) * Basis(1:n) )
        END DO
      END DO

      ! Transform the conductivity tensor (in case of a foil winding):
      ! --------------------------------------------------------------
      IF (CoilType == 'foil winding') C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      DO i=1,nd
        DO j=1,nd
          ! Finally, the elemental matrix & vector:
          !----------------------------------------
          STIFF(i,j) = STIFF(i,j) + SUM(MATMUL(C, dBasisdx(i,:)) * dBasisdx(j,:))*detJ*IP % s(t)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = &
         K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!----------------------------------------------------------------
  SUBROUTINE BoundaryCondition(LOAD, FORCE, Element, n)
!----------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!----------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
    LocalHeatCapacity, LocalDensity
    LOGICAL :: stat, getSecondDerivatives
    INTEGER :: t,j
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!----------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    FORCE = 0.0d0
    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    !-----------------------------------------------------------------
    ! Loop over Gauss-points (boundary element Integration)
    !-----------------------------------------------------------------
    DO t=1,IP % n
      !Basis function values & derivatives at the integration point:
      !-------------------------------------------------------------
      getSecondDerivatives = .FALSE.
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
      IP % W(t), detJ, Basis, dBasisdx)

      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      DO j=1,n
        FORCE(j) = FORCE(j) + IP % s(t)*DetJ*LoadAtIP*Basis(j)
      END DO
    END DO
  END SUBROUTINE BoundaryCondition

!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
   INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
   INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]

   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE Wsolve
!------------------------------------------------------------------------------
