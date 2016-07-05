
SUBROUTINE WaveSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the wave equation 
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
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, istat, active
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, Material
  REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), DAMP(:,:), STIFF(:,:), LOAD(:), FORCE(:), &
      Speed(:)
  LOGICAL :: EigenAnalysis

  SAVE MASS, DAMP, STIFF, LOAD, FORCE, Speed, AllocationsDone
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), Speed(N), STIFF(N,N), DAMP(n,n), MASS(n,n), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  EigenAnalysis = GetLogical( GetSolverParams(),'Eigen Analysis',Found )


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

    Material => GetMaterial()
    Speed(1:n) = GetReal( Material,'Wave Speed')

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(  MASS, DAMP, STIFF, FORCE, LOAD, Speed, Element, n, nd+nb )
    IF( TransientSimulation ) THEN
      CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
    END IF
    CALL DefaultUpdateEquations( STIFF, FORCE )
    IF ( EigenOrHarmonicAnalysis() ) CALL DefaultUpdateMass( MASS )
  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()


  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  MASS, DAMP, STIFF, FORCE, LOAD, Speed, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: MASS(:,:), DAMP(:,:), STIFF(:,:), FORCE(:), LOAD(:), Speed(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,SpeedAtIp
    LOGICAL :: Stat
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    MASS = 0.0d0
    DAMP = 0.0d0
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! The wave speed at the integration point:
      !------------------------------------------
      SpeedAtIp = SUM( Basis(1:n) * Speed(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      DO p=1,nd
        DO q=1,nd
          MASS(p,q) = MASS(p,q) + IP % s(t) * DetJ * Basis(p) * Basis(q)
        END DO
      END DO

      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
          SpeedAtIp * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
    !------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE WaveSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE WaveSolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )

!------------------------------------------------------------------------------
END SUBROUTINE WaveSolver_Init
!------------------------------------------------------------------------------
