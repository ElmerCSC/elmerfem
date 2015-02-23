SUBROUTINE PressureSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!
!  A Navier-Stokes solver with pressure projection, the pressure part
!  15/12/2003, Juha
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
  LOGICAL :: AllocationsDone = .FALSE., Newton = .FALSE., Found
  TYPE(Element_t),POINTER :: Element

  INTEGER :: n, nb, nd, t, istat, dim, active
  REAL(KIND=dp) :: Norm = 0, PrevNorm, RelC, PressureRelax

  TYPE(Variable_t), POINTER :: PressureVariable
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:,:), FORCE(:), Velocity(:,:)

  SAVE STIFF, LOAD, FORCE, Velocity, AllocationsDone
!------------------------------------------------------------------------------

  dim = CoordinateSystemDimension()
  Mesh => GetMesh()

  ! Allocate some permanent storage, this is done first time only:
  !---------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), STIFF(N,N), Velocity(dim,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PressureSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

  ! Initialize the system and do the assembly:
  !-------------------------------------------
  active = GetNOFActive()
  CALL DefaultInitialize()
  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes()
     nd = GetElementNOFDOFs()
     nb = GetElementNOFBDOFs()

     ! Get previous elementwise velocity iterate:
     !-------------------------------------------
     CALL GetVectorLocalSolution( Velocity, 'VelocityTot' )

     ! Get element local matrix and rhs vector:
     !-----------------------------------------
     CALL LocalMatrix(  STIFF, FORCE,  Velocity, &
            Element, n, nd, nd+nb, dim )
     IF ( nb>0 ) CALL LCondensate( nd, nb, STIFF, FORCE )

     ! Update global matrix and rhs vector from local matrix & vector:
     !----------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()
  
!------------------------------------------------------------------------------

  ! Solve the system:
  !------------------
  Norm = DefaultSolve()

  PressureRelax = GetConstReal( GetSolverParams(), 'Pressure Relax', Found )
  IF ( .NOT. Found ) PressureRelax = 1.0d0

  PressureVariable => VariableGet( Mesh % Variables, 'PressureTot')
  PressureVariable % Values = PressureVariable % Values + &
       PressureRelax * Solver % Variable % Values

!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, NodalVelo, Element, n, nd, ntot, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp), TARGET :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: NodalVelo(:,:)
    INTEGER :: dim, n, nd, ntot
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(ntot),dBasisdx(ntot,3), &
                  DetJ, divU
    REAL(KIND=dp), POINTER :: A(:,:), F(:)
    LOGICAL :: Stat
    INTEGER :: t, i, j, k, l, p, q
    TYPE(GaussIntegrationPoints_t) :: IP

    REAL(KIND=dp) :: s, c

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element,ntot )
    DO t=1,IP % n
       ! Basis function values & derivatives at the integration point:
       !--------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
        IP % W(t),  detJ, Basis, dBasisdx )

       s = IP % s(t) * detJ

       ! Previous velocity at the integration point:
       !--------------------------------------------
       DivU = 0.0d0
       DO i=1,dim
          DivU = DivU + SUM( NodalVelo(i,1:nd) * dBasisdx(1:nd,i) )
       END DO

       ! Finally, the elemental matrix & vector:
       !----------------------------------------
       STIFF(1:ntot,1:ntot) = STIFF(1:ntot,1:ntot) + s * &
           MATMUL( dBasisdx, TRANSPOSE(dBasisdx) )
       FORCE(1:ntot) = FORCE(1:ntot) - s * DivU * Basis
    END DO
!   DO i=n+1,nd
!     FORCE(i)   = 0
!     STIFF(i,i) = 1
!   END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE LCondensate( N, nb, K, F )
!------------------------------------------------------------------------------
      USE LinearAlgebra
      INTEGER :: N, nb, dim
      REAL(KIND=dp) :: K(:,:),F(:), Kbb(Nb,Nb), &
       Kbl(nb,n),Klb(n,nb),Fb(nb)

      INTEGER :: m, i, j, l, p, Cdofs(n), Bdofs(nb)

      DO i=1,n
        CDOFs(i) = i
      END DO

      DO i=1,nb
        BDOFs(i) = i+n
      END DO

      Kbb = K(Bdofs,Bdofs)
      Kbl = K(Bdofs,Cdofs)
      Klb = K(Cdofs,Bdofs)
      Fb  = F(Bdofs)

      CALL InvertMatrix( Kbb,Nb )

      F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
      K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )
!------------------------------------------------------------------------------
    END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE PressureSolver
!------------------------------------------------------------------------------
