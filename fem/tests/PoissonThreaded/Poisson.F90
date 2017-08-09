SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!  Solve the Poisson equation with OpenMP threading for assembly and solution
!  phases implemented. Original case by Mikko B. 
!------------------------------------------------------------------------------
  USE Types
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
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: nnz, nthreads, n, nb, nd, t, istat, active
  LOGICAL :: Found
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce
  REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), LOAD(:), FORCE(:)
  INTEGER, POINTER :: Indexes(:)
  INTEGER :: nind, ii, jj
  INTEGER, SAVE :: maxdofs
  LOGICAL, SAVE :: AllocationsDone = .FALSE.
  !$OMP THREADPRIVATE(STIFF, LOAD, FORCE, AllocationsDone, maxdofs)
  ! Variables related to graph colouring
  INTEGER :: col, cli, cti
  TYPE(Graph_t) :: DualGraph
  TYPE(GraphColour_t) :: GraphColouring
  TYPE(Graph_t), POINTER :: ColourIndexList
  INTEGER, ALLOCATABLE :: dualptr(:), dualind(:), colours(:), &
        cptr(:), cind(:)
  LOGICAL :: VecAsm, MCAsm

  REAL(kind=dp), POINTER :: LoadPtr(:)
  !------------------------------------------------------------------------------

  Mesh => GetMesh()

  !System assembly:
  !----------------
  Active = GetNOFActive()
  CALL DefaultInitialize()

  ! Check that mesh colouring is actually in use
  VecAsm = .TRUE.
  MCAsm = ListGetLogical( Solver % Values, 'MultiColour Solver', Found )
  IF (VecAsm .AND. MCAsm .AND. Found) THEN
    CALL Info('PoissonSolver', 'Vectorized assembly in use')

    ! Solver will construct colouring from dual mesh as follows
    ! and store it to Solver % ColourIndexList -variable
    ! (see MainUtils::AddEquationsBasics for details) 

    ! Construct the dual graph from Elmer mesh
    ! CALL ElmerMeshToDualGraph(Mesh, DualGraph)
    
    ! Colour the dual graph
    ! CALL ElmerGraphColour(DualGraph, GraphColouring)
    
    ! Deallocate dual graph as it is no longer needed
    ! CALL Graph_Deallocate(DualGraph)
    
    ! Construct colour lists
    ! CALL ElmerColouringToGraph(GraphColouring, ColourIndexList)
    ! CALL Colouring_Deallocate(GraphColouring)
    ColourIndexList => Solver % ColourIndexList
  ELSE
    CALL Fatal('PoissonSolver', 'Vectorized assembly not in use')
  END IF
  
  CALL ResetTimer('ThreadedAssembly')
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(Solver, Mesh, Active,nthreads, ColourIndexList, cli, cti, VecAsm) &
  !$OMP PRIVATE(BodyForce, Element, col, n, nd, nb, t, istat, Found, Norm, &
  !$OMP         Indexes, nind, LoadPtr)

  !Allocate some permanent storage per thread:
  !--------------------------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    n = Mesh % MaxElementDOFs  ! just big enough for elemental arrays
    maxdofs = n
    ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  ! Perform FE assembly one colour at a time (thread-safe)
  ! Element indices are listed in packed (CRS) structure  ColourIndexList
  DO col=1,ColourIndexList % n

    !$OMP SINGLE
    CALL Info('PoissonSolve','Assembly of colour: '//TRIM(I2S(col)),Level=10)
    cli = ColourIndexList % ptr(col)
    cti = ColourIndexList % ptr(col+1)-1
    ! Set current colour to Solver (normally done via call to GetNOFActive())
    Solver % CurrentColour = col
    !$OMP END SINGLE

    !$OMP DO SCHEDULE(STATIC)
    DO t=cli, cti
      ! Element => GetActiveElement(ColourIndexList % ind(t))
      Element => GetActiveElement(t-cli+1)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)

      LOAD(1:n) = 0.0d0
      BodyForce => GetBodyForce(Element)

      IF ( ASSOCIATED(BodyForce) ) THEN
        CALL GetRealValues( BodyForce, 'Source', Load, Found, UElement=Element )
      END IF

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix( STIFF, FORCE, LOAD, Element, n, nd+nb )
      CALL LCondensate( nd, nb, STIFF, FORCE, maxdofs )

      CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element, VecAssembly=VecAsm )
    END DO ! Element loop
    !$OMP END DO

  END DO ! Colour loop
  !$OMP END PARALLEL

  CALL CheckTimer('ThreadedAssembly',Delete=.TRUE.)

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  Norm = DefaultSolve()

  CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    USE LinearForms
    REAL(KIND=dp) CONTIG :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    LOGICAL :: Stat
    INTEGER :: i, dim, ldbasis, ngp, allocstat
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:), dBasisdx(:,:,:), DetJ(:), LoadAtIPs(:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJ, LoadAtIPs
    !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, DetJ, LoadAtIPs)
!------------------------------------------------------------------------------
    CALL GetElementNodesVec( Nodes, UElement=Element )
    STIFF = 0.0d0
    FORCE = 0.0d0

    ! Get integration points
    IP = GaussPoints( Element )
    ! IP = GaussPoints( Element, 64 )
    ! IP = GaussPoints( Element, 512 )
    ngp = IP % n

    ! Reserve workspace
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), &
              DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LocalMatrix','Storage allocation for local element basis failed')
      END IF
    ELSE IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) THEN
      DEALLOCATE(Basis, dBasisdx, DetJ, LoadAtIPs)
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), &
              DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LocalMatrix','Storage allocation for local element basis failed')
      END IF
    END IF
    
    ldbasis = SIZE(Basis,1)
    ! Compute values of all basis functions at all integration points
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, DetJ, SIZE(Basis,2), Basis, dBasisdx )
    ! Compute actual integration weights (recycle memory space of DetJ)
    DetJ(1:ngp) = IP % s(1:ngp)*Detj(1:ngp)

    ! STIFF=STIFF+(grad u, grad u)
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJ, STIFF)

    ! Source terms at IPs
    !------------------------------------------
    ! LoadAtIPs(1:ngp) = MATMUL( Basis(1:ngp,1:n), LOAD(1:n) )
    CALL LinearForms_ProjectToU(ngp, n, Basis, LOAD, LoadAtIPs)

    ! FORCE=FORCE+(u,f)
    CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, LoadAtIPs, FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

  
! Condensate bubble functions from the local stiffness matrix and force vector.
! Uses BLAS and LAPACK
!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F, kld )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    IMPLICIT NONE

    ! Parameters
    INTEGER, INTENT(IN) :: N, Nb
    REAL(KIND=dp) :: K(kld, kld)
    REAL(KIND=dp) :: F(kld)
    INTEGER, INTENT(IN) :: kld
    
    ! Variables
    REAL(KIND=dp) :: Kbb(Nb,Nb), &
            Kbl(Nb,N), Klb(N,Nb), Fb(Nb), DNRM2
    
    INTEGER :: m, i, j, l, p, nfo, PIV(Nb)
    ! INTEGER :: Ldofs(N), Bdofs(Nb)
    
    ! BLAS interfaces
    INTERFACE
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        IMPLICIT NONE
        DOUBLE PRECISION ALPHA,BETA
        INTEGER INCX,INCY,LDA,M,N
        CHARACTER TRANS
        DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      END SUBROUTINE DGEMV
    END INTERFACE
    
    INTERFACE
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        IMPLICIT NONE
        DOUBLE PRECISION ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,M,N
        CHARACTER TRANSA,TRANSB
        DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      END SUBROUTINE DGEMM
    END INTERFACE
    
    ! LAPACK interfaces
    INTERFACE
      SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
        IMPLICIT NONE
        INTEGER :: INFO, LDA, M, N
        INTEGER :: IPIV( * )
        DOUBLE PRECISION :: A( LDA, * )
      END SUBROUTINE DGETRF
    END INTERFACE
    
    INTERFACE
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        IMPLICIT NONE
        CHARACTER ::TRANS
        INTEGER :: INFO, LDA, LDB, N, NRHS
        INTEGER :: IPIV(*)
        DOUBLE PRECISION :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DGETRS
    END INTERFACE
    
    IF ( Nb <= 0 ) RETURN

    ! Ldofs = (/ (i, i=1,n) /)
    ! Bdofs = (/ (i, i=n+1,n+nb) /)

    ! Can we use blas/lapack without copying the matrices?
    ! Kbb = K(Bdofs,Bdofs)
    DO j=n+1,n+nb
      DO i=n+1,n+nb
        Kbb(i-n,j-n) = K(i,j)
      END DO
    END DO
    ! Kbl = K(Bdofs,Ldofs)
    DO j=1,n
      DO i=n+1,n+nb
        Kbl(i-n,j) = K(i,j)
      END DO
    END DO
    ! Klb = K(Ldofs,Bdofs)
    DO j=n+1,n+nb
      DO i=1,n
        Klb(i,j-n) = K(i,j)
      END DO
    END DO
    ! Fb  = F(Bdofs)
    DO i=n+1,n+nb
      Fb(i-n) = F(i)
    END DO

    ! CALL InvertMatrix( Kbb,nb )
    CALL DGETRF(nb, nb, Kbb, nb, PIV, NFO)

    ! F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    CALL DGETRS('N', Nb, 1, Kbb, Nb, PIV, Fb, Nb, NFO)
    CALL DGEMV('N', n, nb, -1D0, Klb, n, Fb, 1, 1D0, F, 1)

    ! K(1:n,1:n) = &
    !        K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
    CALL DGETRS('N', Nb, n, Kbb, Nb, PIV, Kbl, Nb, NFO)
    CALL DGEMM('N', 'N', n, n, nb, -1D0, Klb, n, Kbl, nb, 1D0, K, kld)
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

END SUBROUTINE PoissonSolver  
!------------------------------------------------------------------------------
! END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------
