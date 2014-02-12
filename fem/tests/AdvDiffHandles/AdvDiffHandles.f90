!------------------------------------------------------------------------------
!> Simplified advection-diffusion solver with constant coefficients.
!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  INTEGER :: n, nb, nd, t, istat, NoActive
  TYPE(Matrix_t), POINTER :: A
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), MASS(:,:), LOAD(:), FORCE(:)
  SAVE STIFF, MASS, LOAD, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  CALL Info( 'AdvDiffSolver','')
  CALL Info( 'AdvDiffSolver','----------------------------------------------------------')
  CALL Info( 'AdvDiffSolver','Solving advection-diffusion equation with material handles')

  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
    N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
    ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), MASS(n,n), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'AdvDiffSolver', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
   END IF
   
   ! System assembly:
   !-----------------
   CALL ResetTimer('AdvDiffSolver::StiffAssembly')
   CALL DefaultInitialize()

   NoActive = GetNOFActive()
   DO t=1,NoActive

      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      ! Get element local matrix and rhs vector:
      ! Usually Elmer assembles stiffness and mass matrices and 
      ! adds them directly into global matrix. Here they are
      ! kept separately for later use.
      !------------------------------------------------------------
      CALL LocalMatrix( STIFF, MASS, FORCE, LOAD, Element, n, nd+nb )

      CALL DefaultUpdateEquations( STIFF, FORCE )
      IF( TransientSimulation ) THEN
        CALL DefaultUpdateMass( MASS )
      END IF
   END DO

   ! Includes 'Linear System FCT'
   CALL DefaultFinishAssembly()
   CALL CheckTimer('AdvDiffSolver::StiffAssembly',Delete=.TRUE.)

   IF (TransientSimulation) THEN     
     CALL ResetTimer('AdvDiffSolver::MassAssembly')
     IF(.TRUE.) THEN
       A => GetMatrix()      
       CALL Add1stOrderTime_CRS( A, A % rhs, dt, Solver )
     ELSE
       CALL Default1stOrderTimeGlobal()
     END IF
     CALL CheckTimer('AdvDiffSolver::MassAssembly',Delete=.TRUE.)
   END  IF

   CALL DefaultDirichletBCs()
   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()
   
   CALL Info( 'AdvDiffSolver','Add done')
   CALL Info( 'AdvDiffSolver','----------------------------------------------------------')


CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, MASS, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), MASS(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: Velo(3,nd), Velo_AtIp(3)
    REAL(KIND=dp) :: D_AtIp, C_AtIp, M_AtIp, S_AtIp
    LOGICAL :: Stat
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t) :: S_handle, D_handle, C_handle, M_handle
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName

    SAVE Nodes, S_handle, D_handle, C_handle, M_handle

!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    VarName = ComponentName( Solver % Variable )

    STIFF = 0._dp
    MASS = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
    Velo = 0._dp
    CALL GetScalarLocalSolution(Velo(1,:),'Velocity 1')
    CALL GetScalarLocalSolution(Velo(2,:),'Velocity 2')
    CALL GetScalarLocalSolution(Velo(3,:),'Velocity 3')

    BodyForce => GetBodyForce()
    Material => GetMaterial()

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      S_AtIp = ListGetRealAtIp( S_Handle,BodyForce,Basis, &
          TRIM(VarName)//' Source')         

      ! Material parameters at integration point
      !--------------------------------------------------
      D_AtIp = ListGetRealAtIp( D_Handle,Material,Basis,&
          TRIM(VarName)//' Diffusivity')
!      C_AtIp = ListGetRealAtIp( C_Handle,Material,Basis,&
!          TRIM(ComponentName(Solver % Variable))//' Capacity')
      C_AtIp = 1.0_dp

!      M_AtIp = ListGetRealAtIp( M_Handle,Material,Basis, &
!          TRIM(ComponentName(Solver % Variable))//' Density')
      M_AtIp = 1.0_dp

      ! Velocity at integration point
      !------------------------------
      DO i=1,dim
        Velo_AtIp(i) = SUM(Velo(i,1:n) * Basis(1:n))
      END DO

      ! Finally, the elemental matrix & vector:
      ! First the diffusive part
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             D_AtIp * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      ! Convective part
      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + IP % s(t) * DetJ * &
             C_AtIp * SUM(Velo_AtIp * dBasisdx(q,:)) * Basis(p)
        END DO
      END DO

      ! Mass matrix
      IF( TransientSimulation ) THEN
        DO p=1,nd
          DO q=1,nd
            MASS(p,q) = MASS(p,q) + IP % s(t) * DetJ * &
                M_AtIp * Basis(q) * Basis(p)
          END DO
        END DO
      END IF

      ! Source term
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * S_AtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

END SUBROUTINE AdvDiffSolver
