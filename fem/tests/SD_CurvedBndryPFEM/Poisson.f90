SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
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
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm,  Energy
  INTEGER :: i,j,n, nb, nd, t, istat, Active

  TYPE(ValueList_t), POINTER :: BodyForce,BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
               NEWX(:), NEWY(:), POT(:)

  TYPE(Matrix_t), POINTER :: A
  TYPE(Mesh_t),   POINTER :: Mesh

  SAVE STIFF, LOAD, FORCE, POT, AllocationsDone
!------------------------------------------------------------------------------

  ! Allocate some permanent storage, this is done first time only:
  ! --------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), POT(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

   n = SIZE( Mesh % Nodes % x )
   ALLOCATE( NEWX(n),NEWY(n) )
   NEWX = Mesh % Nodes % x
   NEWY = Mesh % Nodes % y

   LOAD=0
   Active = GetNOFActive()
   CALL DefaultInitialize()
   DO t=1,Active
     Element => GetActiveElement(t)
     nd = GetElementNOFDOFs(  Element )
     n  = GetElementNOFNodes( Element )
     CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   CALL DefaultFinishAssembly()

   Element => GetBoundaryElement(1)
   BC => GetBC()

   CALL ListAddConstReal( BC,'Potential',0.0d0,GetProcAddr("./Poisson circx") )
   CALL DefaultDirichletBCs()
   Norm = DefaultSolve()
   NEWX = Solver % Variable % Values(Solver % Variable % Perm)

   CALL DefaultInitialize()
   DO t=1,Active
     Element => GetActiveElement(t)
     nd = GetElementNOFDOFs(  Element )
     n  = GetElementNOFNodes( Element )
     CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   CALL DefaultFinishAssembly()
   CALL ListAddConstReal( BC,'Potential',0.0d0,GetProcAddr("./Poisson circy") )
   CALL DefaultDirichletBCs()
   Norm = DefaultSolve()
   NEWY = Solver % Variable % Values(Solver % Variable % Perm)

   Mesh % Nodes % x = NEWX
   Mesh % Nodes % y = NEWY
   DEALLOCATE( NEWX,NEWY )

   ! System assembly:
   ! ----------------
   CALL DefaultInitialize()
   DO t=1,active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()

      BodyForce => GetBodyForce()
      LOAD = 0.0d0
      IF ( ASSOCIATED(BodyForce) ) &
        Load(1:n) = GetReal( BodyForce, 'Source', Found )

      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   CALL DefaultFinishAssembly()
   CALL ListAddConstReal( BC,'Potential',0.0d0)
   CALL DefaultDirichletBCs()
   Norm = DefaultSolve()

   Energy  = 0.0d0
   LOAD    = 0.0d0
   DO t=1,active
      Element => GetActiveElement(t)
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      CALL LocalMatrix( STIFF, FORCE, LOAD, Element, n, nd )
      CALL GetScalarLocalSolution(POT)
      Energy = Energy + SUM( POT*MATMUL(STIFF,POT) )
   END DO
   PRINT*,Energy
   A => GetMatrix()
   PRINT*,'DOFs, Error in energy: ',A % NumberOfRows, ABS(Energy-PI/2.0d0)

CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),ddBasisddx(nd,3,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
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

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL( dBasisdx, TRANSPOSE(dBasisdx) )
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION CircX(Model,x,y,z) RESULT(s)
!------------------------------------------------------------------------------
  USE Types
  TYPE(Model_t) :: Model
  REAL(KIND=dp)::x,y,z,s,r

  r = SQRT(x**2 + y**2 )
  s = 0.0d0
  IF ( r /= 0.0d0 ) s = SQRT(2.0d0)*x/r
!------------------------------------------------------------------------------
END FUNCTION CircX
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION CircY(Model,x,y,z) RESULT(s)
!------------------------------------------------------------------------------
  USE Types
  TYPE(Model_t) :: Model
  REAL(KIND=dp)::x,y,z,s,r

  r = SQRT(x**2 + y**2 )
  s = 0.0d0
  IF ( r /= 0.0d0 ) s = SQRT(2.0d0)*y/r
!------------------------------------------------------------------------------
END FUNCTION CircY
!------------------------------------------------------------------------------
