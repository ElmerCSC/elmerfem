SUBROUTINE H1BasisEvaluation( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Unit test for new H1 basis function routines in Elmer
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
    USE ISO_C_BINDING
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver

    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
    REAL(kind=dp), PARAMETER :: tol1d = 1E-12, tol3d=1E-12
    INTEGER :: nerror, netest

    nerror = 0

    ! 1D tests
    netest = TestLineElement(Solver, tol1d)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Line element contained errors')
    END IF
    nerror = nerror + netest
    
    ! 2D tests 
    netest = TestTriangleElement(Solver, tol1d, .FALSE.)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Triangle element contained errors')
    END IF
    nerror = nerror + netest

    ! 2D element in a 3D mesh
    netest = TestTriangleElement(Solver, tol1d, .TRUE.)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Triangle face element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestQuadElement(Solver, tol1d, .FALSE.)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Quad element contained errors')
    END IF
    nerror = nerror + netest

    ! 2D element in a 3D mesh
    netest = TestQuadElement(Solver, tol1d, .TRUE.)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Quad face element contained errors')
    END IF
    nerror = nerror + netest

    ! 3D tests
    netest = TestTetraElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Tetra element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestWedgeElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Wedge element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestBrickElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1BasisEvaluation','Brick element contained errors')
    END IF
    nerror = nerror + netest

    ! Build solution norm for error checking
    Solver % Variable % Norm = REAL(1+nerror,dp)
    Solver % Variable % Values = REAL(1+nerror,dp)
    
CONTAINS

  FUNCTION TestLineElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim, perm, q, tag
    INTEGER, PARAMETER :: P = 6, NREP = 100, BubblePerm = 2
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    LOGICAL :: Invert(BubblePerm)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 202, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF

    GP = GaussPoints(Element)

    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof * BubblePerm
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF

    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    Invert(1:2) = [.FALSE., .TRUE.]

    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions1D( Basis(i,1:nndof), Element, UWrk(i))
        CALL NodalFirstDerivatives1D( dBasisdx(i,1:nndof,1:3), Element, UWrk(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis (with and without inversion)
      q = 2
      t_start_tmp=ftimer()
      DO perm=1,BubblePerm
        DO j=1, Element % BDOFs
          q = q + 1
          DO i=1,ngp
            Basis(i,q) = LineBubblePBasis(j+1, UWrk(i),Invert(perm))
            dBasisdx(i,q,1) = dLineBubblePBasis(j+1, UWrk(i), Invert(perm))
          END DO
        END DO
      END DO
      t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
    END DO
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_LineNodal(ngp, UWrk, nbasisvec, BasisVec)
      CALL H1Basis_dLineNodal(ngp, UWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      DO perm=1,BubblePerm
        CALL H1Basis_LineBubbleP(ngp, UWrk, P, nbasisvec, BasisVec, Invert(perm))
        CALL H1Basis_dLineBubbleP(ngp, UWrk, P, ndbasisdxvec, dBasisdxVec, &
              Invert(perm))
      END DO
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()

    CALL PrintTestData(Element, ngp, nrep, nbdof, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestLineElement

  FUNCTION TestTriangleElement(Solver, tol, isEdge) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    LOGICAL :: isEdge

    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nedof, nfdof, &
            nbdof, allocstat, nbasisvec, ndbasisdxvec, rep, dim, perm, tag
    INTEGER, PARAMETER :: P = 6, NREP = 100, EdgePerm = 2, BubblePerm=3
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, t_tot_e, t_totvec_e, &
            t_tot_f, t_totvec_f, t_tot_b, t_totvec_b
    
    INTEGER :: EdgeDir(H1Basis_MaxPElementEdgeNodes,&
                       H1Basis_MaxPElementEdges,&
                       EdgePerm), EdgeP(H1Basis_MaxPElementEdges), &
               BubbleDir(H1Basis_MaxPElementFaceNodes, &
                       BubblePerm)
    LOGICAL :: InvertEdge(3, EdgePerm)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 303, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF
    GP = GaussPoints(Element,4)
    
    nndof = Element % Type % NumberOfNodes
    nedof = GetEdgeDOFs( Element, P )
    IF (isEdge) THEN
      nbdof = Element % BDofs*BubblePerm
    ELSE
      nbdof = Element % BDofs
    END IF
    nbasis = nndof + 3*nedof*EdgePerm + nbdof
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF
    
    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    DO i=1,3
      EdgeDir(1:2,i,1)=getTriangleEdgeMap(i)
    END DO
    ! Invert direction
    DO i=1,3
      EdgeDir(2:1:-1,i,2)=EdgeDir(1:2,i,1)
    END DO
    InvertEdge(:,1) = .FALSE.
    InvertEdge(:,2) = .TRUE.

    ! Initialize bubble mappings if needed
    IF (isEdge) THEN
      BubbleDir(1:3,1)=[ 1,2,3 ]
      DO i=2,BubblePerm
        BubbleDir(1:3,i)=CSHIFT(BubbleDir(1:3,1),i-1)
      END DO
    END IF

    EdgeP = P

    t_tot_n = REAL(0,dp)
    t_tot_e = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = TriangleNodalPBasis(ndof, UWrk(i), VWrk(i))
          dBasisdx(i,ndof,1:2) = dTriangleNodalPBasis(ndof, UWrk(i), VWrk(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Edge basis
      t_start_tmp=ftimer()
      q = 3
      DO perm=1,EdgePerm
        DO i=1,3
          DO ndof=1,nedof
            q=q+1
            DO k=1,ngp
              Basis(k,q) = TriangleEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), &
                                              InvertEdge(i,perm))
              dBasisdx(k,q,1:2) = dTriangleEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), &
                                              InvertEdge(i,perm))
            END DO
          END DO
        END DO
      END DO
      t_tot_e=t_tot_e+(ftimer()-t_start_tmp)
        
      ! Bubble basis
      IF (.NOT. isEdge) THEN
        t_start_tmp=ftimer()
        DO i = 0,p-3
          DO j = 0,p-i-3
            q = q + 1
            DO k = 1, ngp
              Basis(k, q) = TriangleBubblePBasis(i,j,UWrk(k), GP % v(k))
              dBasisdx(k, q, 1:2) = dTriangleBubblePBasis(i,j,GP % u(k), GP % v(k))  
            END DO
          END DO
        END DO
        t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
      ELSE
        t_start_tmp=ftimer()
        DO perm=1,BubblePerm
          DO i = 0,p-3
            DO j = 0,p-i-3
              q = q + 1
              DO k = 1, ngp
                Basis(k, q) = TriangleEBubblePBasis(i,j,UWrk(k), GP % v(k), &
                        BubbleDir(1:3,perm))
                dBasisdx(k, q, 1:2) = dTriangleEBubblePBasis(i,j,GP % u(k), GP % v(k), &
                        BubbleDir(1:3,perm))  
              END DO
            END DO
          END DO
        END DO
        t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
      END IF
    END DO ! NREP
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_e = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_TriangleNodalP(ngp, UWrk, VWrk, nbasisvec, BasisVec)
      CALL H1Basis_dTriangleNodalP(ngp, UWrk, VWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      DO perm=1,EdgePerm
        CALL H1Basis_TriangleEdgeP(ngp, UWrk, VWrk, EdgeP, &
                nbasisvec, BasisVec, EdgeDir(:,:,perm))
        CALL H1Basis_dTriangleEdgeP(ngp, UWrk, VWrk, EdgeP, &
                ndbasisdxvec, dBasisdxVec, EdgeDir(:,:,perm))
      END DO
      t_totvec_e=t_totvec_e+(ftimer()-t_start_tmp)
      
      IF (.NOT. isEdge) THEN
        t_start_tmp=ftimer()
        CALL H1Basis_TriangleBubbleP(ngp, UWrk, VWrk, P, &
                nbasisvec, BasisVec)
        CALL H1Basis_dTriangleBubbleP(ngp, UWrk, VWrk, P, &
                ndbasisdxvec, dBasisdxVec)
        t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
      ELSE
        t_start_tmp=ftimer()
        DO perm=1,BubblePerm
          CALL H1Basis_TriangleBubbleP(ngp, UWrk, VWrk, P, &
                  nbasisvec, BasisVec, BubbleDir(1:3,perm))
          CALL H1Basis_dTriangleBubbleP(ngp, UWrk, VWrk, P, &
                  ndbasisdxvec, dBasisdxVec, BubbleDir(1:3,perm))
        END DO
        t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
      END IF

    END DO
    t_endvec = ftimer()

    CALL PrintTestData(Element, ngp, nrep, 3*nedof+nbdof, &
            t_tot_n, t_tot_e+t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_e+t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestTriangleElement
  
  FUNCTION TestQuadElement(Solver, tol, isEdge) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    LOGICAL :: isEdge

    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nedof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim, perm, tag
    INTEGER, PARAMETER :: P = 6, NREP = 100, EdgePerm = 2, BubblePerm = 4
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_e, t_totvec_e, t_tot_b, t_totvec_b
    INTEGER :: EdgeDir(H1Basis_MaxPElementEdgeNodes,&
                       H1Basis_MaxPElementEdges,&
                       EdgePerm), EdgeP(H1Basis_MaxPElementEdges),&
               BubbleDir(H1Basis_MaxPElementFaceNodes,BubblePerm)
    LOGICAL :: InvertEdge(4, EdgePerm)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec, EdgeDir
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 404, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF

    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nedof = getEdgeDOFs( Element, P )
    IF (isEdge) THEN
      nbdof = Element % BDofs*BubblePerm
    ELSE
      nbdof = Element % BDofs
    END IF
    nbasis = nndof + 4*nedof*EdgePerm + nbdof
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF
    
    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    DO i=1,4
      EdgeDir(1:2,i,1)=getQuadEdgeMap(i)
    END DO
    ! Invert direction
    DO i=1,4
      EdgeDir(2:1:-1,i,2)=EdgeDir(1:2,i,1)
    END DO
    InvertEdge(:,1) = .FALSE.
    InvertEdge(:,2) = .TRUE.

    ! Initialize face mappings if needed
    IF (isEdge) THEN
      BubbleDir(1:4,1)=[ 1,2,3,4 ]
      DO i=2,BubblePerm
        BubbleDir(1:4,i)=CSHIFT(BubbleDir(1:4,1),i-1)
      END DO
    END IF

    EdgeP = P
    
    t_tot_n = REAL(0,dp)
    t_tot_e = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions2D( Basis(i,1:nndof), Element, UWrk(i), VWrk(i))
        CALL NodalFirstDerivatives2D( dBasisdx(i,1:nndof,1:3), Element, UWrk(i), VWrk(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Edge basis
      t_start_tmp=ftimer()
      q = 4
      DO perm=1,EdgePerm
        DO i=1,4
          DO ndof=1,nedof
            q=q+1
            DO k=1,ngp
              Basis(k,q) = QuadEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), InvertEdge(i,perm))
              dBasisdx(k,q,1:2) = dQuadEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), InvertEdge(i,perm))
            END DO
          END DO
        END DO
      END DO
      t_tot_e=t_tot_e+(ftimer()-t_start_tmp)

      ! Bubble basis 
      IF (.NOT. isEdge) THEN
        t_start_tmp=ftimer()
        DO i=2,(p-2)
          DO j=2,(p-i)
            q = q + 1
            DO k = 1, ngp
              Basis(k, q) = QuadBubblePBasis(i,j,UWrk(k),VWrk(k))
              dBasisdx(k, q, 1:2) = dQuadBubblePBasis(i,j,UWrk(k),VWrk(k))
            END DO
          END DO
        END DO
        t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
      ELSE
        t_start_tmp=ftimer()
        DO perm=1,BubblePerm
          DO i=2,(p-2)
            DO j=2,(p-i)
              q = q + 1
              DO k = 1, ngp
                Basis(k, q) = QuadBubblePBasis(i,j,UWrk(k),VWrk(k),&
                                               BubbleDir(1:4,perm))
                dBasisdx(k, q, 1:2) = dQuadBubblePBasis(i,j,UWrk(k),VWrk(k),&
                                               BubbleDir(1:4,perm))
              END DO
            END DO
          END DO
        END DO
        t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
      END IF
    END DO
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_e = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_QuadNodal(ngp, UWrk, VWrk, nbasisvec, BasisVec)
      CALL H1Basis_dQuadNodal(ngp, UWrk, VWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      
      ! Edge basis
      t_start_tmp=ftimer()
      DO perm=1,EdgePerm
        CALL H1Basis_QuadEdgeP(ngp, UWrk, VWrk, EdgeP, &
                nbasisvec, BasisVec, EdgeDir(:,:,perm))
        CALL H1Basis_dQuadEdgeP(ngp, UWrk, VWrk, EdgeP, &
                ndbasisdxvec, dBasisdxVec, EdgeDir(:,:,perm))
      END DO
      t_totvec_e=t_totvec_e+(ftimer()-t_start_tmp)

      IF (.NOT. isEdge) THEN
        t_start_tmp=ftimer()
        CALL H1Basis_QuadBubbleP(ngp, UWrk, VWrk, P, &
                nbasisvec, BasisVec)
        
        CALL H1Basis_dQuadBubbleP(ngp, UWrk, VWrk, P, &
                ndbasisdxvec, dBasisdxVec)
        t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
      ELSE
        t_start_tmp=ftimer()
        DO perm=1,BubblePerm
          CALL H1Basis_QuadBubbleP(ngp, UWrk, VWrk, P, &
                  nbasisvec, BasisVec, BubbleDir(1:4,perm))
          
          CALL H1Basis_dQuadBubbleP(ngp, UWrk, VWrk, P, &
                  ndbasisdxvec, dBasisdxVec, BubbleDir(1:4,perm))
        END DO
        t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
      END IF
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, 4*nedof+nbdof,&
            t_tot_n, t_tot_e+t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_e+t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestQuadElement

  FUNCTION TestTetraElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, k, l, q, ndof, nedof, nfdof, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim, tag, perm
    INTEGER, PARAMETER :: P = 6, NREP = 100, EdgePerm=2, FacePerm=2
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
          t_start_tmp, t_tot_n, t_totvec_n, &
          t_tot_e, t_totvec_e, t_tot_f, t_totvec_f, t_tot_b, t_totvec_b
    INTEGER :: EdgeDir(H1Basis_MaxPElementEdgeNodes,&
                       H1Basis_MaxPElementEdges,&
                       EdgePerm), EdgeP(H1Basis_MaxPElementEdges), &
               FaceDir(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces, &
                       FacePerm), FaceP(H1Basis_MaxPElementFaces), &
               TetraType(H1Basis_MaxPElementEdges, EdgePerm)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec, EdgeDir, FaceDir
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 504, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF

    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nedof = getEdgeDOFs( Element, P )
    nfdof = getFaceDOFs( Element, P )
    nbdof = Element % BDofs
    nbasis = nndof + 6*nedof*EdgePerm + 4*nfdof*FacePerm + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF
    
    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    ! Tetra type 1 and 2
    DO perm=1,EdgePerm
      DO i=1,6
        EdgeDir(1:2,i,perm)=getTetraEdgeMap(i,perm)
      END DO
    END DO
    FaceDir = 0
    DO perm=1,FacePerm
      DO i=1,4
        FaceDir(1:3,i,perm)=getTetraFaceMap(i,perm)
      END DO
    END DO
    TetraType(:,1) = 1
    TetraType(:,2) = 2
    EdgeP = P
    FaceP = P
    
    t_tot_n = REAL(0,dp)
    t_tot_e = REAL(0,dp)
    t_tot_f = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = TetraNodalPBasis(ndof, UWrk(i), VWrk(i), WWrk(i))
          dBasisdx(i,ndof,1:3) = dTetraNodalPBasis(ndof, UWrk(i), VWrk(i), WWrk(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Edge basis
      t_start_tmp=ftimer()
      q = 4
      DO perm=1,EdgePerm
        DO i=1,6
          DO ndof=1,nedof
            q= q + 1
            DO k=1,ngp
              Basis(k,q) = TetraEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), TetraType(i,perm))
              dBasisdx(k,q,1:3) = dTetraEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), TetraType(i,perm))
            END DO
          END DO
        END DO
      END DO
      t_tot_e=t_tot_e+(ftimer()-t_start_tmp)
      
      ! Face basis
      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        DO i=1,4
          DO j=0,p-3
            DO k=0,p-j-3
              q= q + 1
              DO l=1,ngp
                Basis(l,q) = TetraFacePBasis(i, j, k, UWrk(l), VWrk(l), WWrk(l), TetraType(i,perm))
                dBasisdx(l,q,1:3) = dTetraFacePBasis(i, j, k, UWrk(l), VWrk(l), WWrk(l), TetraType(i,perm))
              END DO
            END DO
          END DO
        END DO
      END DO
      t_tot_f=t_tot_e+(ftimer()-t_start_tmp)
      
      ! Bubble basis 
      t_start_tmp=ftimer()
      DO i=0,p-4
        DO j=0,p-i-4
          DO k=0,p-i-j-4
            q = q + 1
            DO l = 1, ngp
              Basis(l, q) = TetraBubblePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l))
              dBasisdx(l, q, 1:3) = dTetraBubblePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l))  
            END DO
          END DO
        END DO
      END DO
      t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
    END DO
    t_end = ftimer()

    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_e = REAL(0,dp)
    t_totvec_f = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_TetraNodalP(ngp, UWrk, VWrk, WWrk, nbasisvec, BasisVec)
      CALL H1Basis_dTetraNodalP(ngp, UWrk, VWrk, WWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      DO perm=1,EdgePerm
        CALL H1Basis_TetraEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              nbasisvec, BasisVec, EdgeDir(:,:,perm))
        CALL H1Basis_dTetraEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              ndbasisdxvec, dBasisdxVec, EdgeDir(:,:,perm))
      END DO
      t_totvec_e=t_totvec_e+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        CALL H1Basis_TetraFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              nbasisvec, BasisVec, FaceDir(:,:,perm))
        CALL H1Basis_dTetraFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              ndbasisdxvec, dBasisdxVec, FaceDir(:,:,perm))
      END DO
      t_totvec_f=t_totvec_f+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      CALL H1Basis_TetraBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              nbasisvec, BasisVec)
      CALL H1Basis_dTetraBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()

    CALL PrintTestData(Element, ngp, nrep, 6*nedof+4*nfdof+nbdof, &
            t_tot_n, t_tot_e+t_tot_f+t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_e+t_totvec_f+t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestTetraElement

  
  FUNCTION TestWedgeElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, k, l, q, ngp, nerror, nbasis, nndof, nedof, &
          nfdof, nbdof, allocstat, perm, &
          nbasisvec, ndbasisdxvec, rep, dim, tag, ndof
    INTEGER, PARAMETER :: P = 6, NREP = 100, EdgePerm=2, FacePerm=4
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_e, t_totvec_e, t_tot_f, t_totvec_f, t_tot_b, t_totvec_b
    INTEGER :: EdgeDir(H1Basis_MaxPElementEdgeNodes,&
                       H1Basis_MaxPElementEdges,&
                       EdgePerm), EdgeP(H1Basis_MaxPElementEdges), &
               FaceDir(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces, &
                       FacePerm), FaceP(H1Basis_MaxPElementFaces)
    INTEGER :: direction(H1Basis_MaxPElementFaceNodes), tmp
    LOGICAL :: InvertEdge(9, EdgePerm), invert
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec, EdgeDir, FaceDir
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 706, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF

    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nedof = getEdgeDOFs( Element, P )
    nfdof = 2*getFaceDofs( Element, P, 1)+3*getFaceDOFs( Element, P, 3 )
    nbdof = Element % BDofs
    nbasis = nndof + 9*nedof*EdgePerm + nfdof*FacePerm + nbdof
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF
    
    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    DO i=1,9
      EdgeDir(1:2,i,1)=getWedgeEdgeMap(i)
    END DO
    ! Invert direction
    DO i=1,9
      EdgeDir(2:1:-1,i,2)=EdgeDir(1:2,i,1)
    END DO
    InvertEdge(:,1) = .FALSE.
    InvertEdge(:,2) = .TRUE.

    DO i=1,5
      FaceDir(1:4,i,1)=getWedgeFaceMap(i)
    END DO
    DO j=2,FacePerm
      DO i=1,2
        FaceDir(1:3,i,j)=CSHIFT(FaceDir(1:3,i,1), j-1)
      END DO
      DO i=3,5
        FaceDir(1:4,i,j)=CSHIFT(FaceDir(1:4,i,1), j-1)
      END DO
    END DO

    EdgeP = P
    FaceP = P

    t_tot_n = REAL(0,dp)
    t_tot_e = REAL(0,dp)
    t_tot_f = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = WedgeNodalPBasis(ndof, UWrk(i), VWrk(i), WWrk(i))
          dBasisdx(i,ndof,1:3) = dWedgeNodalPBasis(ndof, UWrk(i), VWrk(i), WWrk(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Edge basis
      t_start_tmp=ftimer()
      q = 6
      DO perm=1,EdgePerm
        DO i=1,9
          DO ndof=1,nedof
            q=q+1
            DO k=1,ngp
              Basis(k,q) = WedgeEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), &
                      InvertEdge(i,perm))
              dBasisdx(k,q,1:3) = dWedgeEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), &
                      InvertEdge(i,perm))
            END DO
          END DO
        END DO
      END DO
      t_tot_e=t_tot_e+(ftimer()-t_start_tmp)

      ! Face basis
      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        ! Triangle faces
        DO i=1,2
          DO j=0,p-3
            DO k=0,p-j-3
              q = q + 1
              DO l=1,ngp
                Basis(l,q) = WedgeFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                        FaceDir(1:4,i,perm))
                dBasisdx(l,q,1:3) = dWedgeFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                        FaceDir(1:4,i,perm))
              END DO
            END DO
          END DO
        END DO
        ! Square faces
        DO i=3,5
          ! First and second node must form an edge in the upper or lower triangle
          direction(1:4) = FaceDir(1:4,i,perm)
          invert=.FALSE. 
          IF (.NOT. wedgeOrdering(direction)) THEN
            invert = .TRUE.
            tmp = direction(2)
            direction(2) = direction(4)
            direction(4) = tmp
          END IF

          DO j=2,p-2
            DO k=2,p-j
              q = q + 1
              IF (.NOT. invert) THEN
                DO l=1,ngp
                  Basis(l,q) = WedgeFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                          direction)
                  dBasisdx(l,q,1:3) = dWedgeFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                          direction)
                END DO
              ELSE
                DO l=1,ngp
                  Basis(l,q) = WedgeFacePBasis(i,k,j,UWrk(l), VWrk(l), WWrk(l), &
                          direction)
                  dBasisdx(l,q,1:3) = dWedgeFacePBasis(i,k,j,UWrk(l), VWrk(l), WWrk(l), &
                          direction)
                END DO
              END IF
            END DO
          END DO
        END DO
      END DO
      t_tot_f=t_tot_f+(ftimer()-t_start_tmp)
      
      ! Bubble basis 
      t_start_tmp=ftimer()
      DO i=0,p-5
        DO j=0,p-5-i
          DO k=2,p-3-i-j
            q = q + 1
            DO l = 1, ngp
              Basis(l, q) = WedgeBubblePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l))
              dBasisdx(l, q, 1:3) = dWedgeBubblePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l))  
            END DO
          END DO
        END DO
      END DO
      t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
    END DO
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_e = REAL(0,dp)
    t_totvec_f = REAL(0,dp)
    t_totvec_b = REAL(0,dp)    
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_WedgeNodalP(ngp, UWrk, VWrk, WWrk, nbasisvec, BasisVec)
      CALL H1Basis_dWedgeNodalP(ngp, UWrk, VWrk, WWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      t_start_tmp=ftimer()
      DO perm=1,EdgePerm
        CALL H1Basis_WedgeEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              nbasisvec, BasisVec, EdgeDir(:,:,perm))
        CALL H1Basis_dWedgeEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              ndbasisdxvec, dBasisdxVec, EdgeDir(:,:,perm))
      END DO
      t_totvec_e=t_totvec_e+(ftimer()-t_start_tmp)    
      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        CALL H1Basis_WedgeFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              nbasisvec, BasisVec, FaceDir(:,:,perm))
        CALL H1Basis_dWedgeFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              ndbasisdxvec, dBasisdxVec, FaceDir(:,:,perm))
      END DO
      t_totvec_f=t_totvec_f+(ftimer()-t_start_tmp)    
      t_start_tmp=ftimer()
      CALL H1Basis_WedgeBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              nbasisvec, BasisVec)
      CALL H1Basis_dWedgeBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, 9*nedof+nfdof+nbdof, &
          t_tot_n, t_tot_e+t_tot_f+t_tot_b, t_end-t_start, &
          t_totvec_n, t_totvec_e+t_totvec_f+t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestWedgeElement  
  
  FUNCTION TestBrickElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:), UWrk(:), VWrk(:), WWrk(:)

    INTEGER :: i, j, k, l, q, ngp, nerror, nbasis, nndof, nedof, &
            nfdof, nbdof, allocstat, perm, &
            nbasisvec, ndbasisdxvec, rep, dim, tag, ndof
    INTEGER, PARAMETER :: P = 6, NREP = 100, EdgePerm=2, FacePerm=4
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_e, t_totvec_e, t_tot_f, t_totvec_f, t_tot_b, t_totvec_b
    INTEGER :: EdgeDir(H1Basis_MaxPElementEdgeNodes,&
                       H1Basis_MaxPElementEdges,&
                       EdgePerm), EdgeP(H1Basis_MaxPElementEdges), &
               FaceDir(H1Basis_MaxPElementFaceNodes, &
                       H1Basis_MaxPElementFaces, &
                       FacePerm), FaceP(H1Basis_MaxPElementFaces)
    LOGICAL :: InvertEdge(12, EdgePerm)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec, EdgeDir, FaceDir
!DIR$ ATTRIBUTES ALIGN:64 :: UWrk, VWrk, WWrk

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 808, P)
    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = Element % Type % ElementCode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF

    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nedof = getEdgeDOFs( Element, P )
    nfdof = getFaceDOFs( Element, P )
    nbdof = Element % BDofs
    nbasis = nndof + 12*nedof*EdgePerm + 6*nfdof*FacePerm + nbdof
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            UWrk(ngp), VWrk(ngp), WWrk(ngp), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1BasisEvaluation',&
              'Storage allocation for local element basis failed')
    END IF
    
    ! Copy Gauss points to local arrays
    UWrk(1:ngp) = GP % U(1:ngp)
    VWrk(1:ngp) = GP % V(1:ngp)
    WWrk(1:ngp) = GP % W(1:ngp)

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    DO i=1,12
      EdgeDir(1:2,i,1)=getBrickEdgeMap(i)
    END DO
    ! Invert direction
    DO i=1,12
      EdgeDir(2:1:-1,i,2)=EdgeDir(1:2,i,1)
    END DO
    InvertEdge(:,1) = .FALSE.
    InvertEdge(:,2) = .TRUE.

    DO i=1,6
      FaceDir(1:4,i,1)=getBrickFaceMap(i)
    END DO
    DO j=2,FacePerm
      DO i=1,6
        FaceDir(1:4,i,j)=CSHIFT(FaceDir(1:4,i,1), j-1)
      END DO
    END DO

    EdgeP = P
    FaceP = P

    t_tot_n = REAL(0,dp)
    t_tot_e = REAL(0,dp)
    t_tot_f = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions3D( Basis(i,1:nndof), Element, &
                UWrk(i), VWrk(i), WWrk(i))
        CALL NodalFirstDerivatives3D( dBasisdx(i,1:nndof,1:3), Element, &
                UWrk(i), VWrk(i), WWrk(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Edge basis
      t_start_tmp=ftimer()
      q = 8
      DO perm=1,EdgePerm
        DO i=1,12
          DO ndof=1,nedof
            q=q+1
            DO k=1,ngp
              Basis(k,q) = BrickEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), &
                      InvertEdge(i,perm))
              dBasisdx(k,q,1:3) = dBrickEdgePBasis(i, ndof+1, UWrk(k), VWrk(k), WWrk(k), &
                      InvertEdge(i,perm))
            END DO
          END DO
        END DO
      END DO
      t_tot_e=t_tot_e+(ftimer()-t_start_tmp)

      ! Face basis
      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        DO i=1,6
          DO j=2,FaceP(i)-2
            DO k=2,FaceP(i)-j
              q = q + 1
              DO l=1,ngp
                Basis(l,q) = BrickFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                        FaceDir(1:4,i,perm))
                dBasisdx(l,q,1:3) = dBrickFacePBasis(i,j,k,UWrk(l), VWrk(l), WWrk(l), &
                        FaceDir(1:4,i,perm))
              END DO
            END DO
          END DO
        END DO
      END DO
      t_tot_f=t_tot_f+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      DO i=2,p-4
        DO j=2,p-i-2
          DO k=2,p-i-j
            q = q + 1
            DO l=1,ngp
              Basis(l,q) = BrickBubblePBasis(i,j,k, &
                      UWrk(l), VWrk(l), WWrk(l))
              dBasisdx(l,q,:) = dBrickBubblePBasis(i,j,k, &
                      UWrk(l), VWrk(l), WWrk(l))
            END DO
          END DO
        END DO
      END DO
      t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
    END DO
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_e = REAL(0,dp)
    t_totvec_f = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      nbasisvec = 0
      ndbasisdxvec = 0
      t_start_tmp=ftimer()
      CALL H1Basis_BrickNodal(ngp, UWrk, VWrk, WWrk, nbasisvec, BasisVec)
      CALL H1Basis_dBrickNodal(ngp, UWrk, VWrk, WWrk, ndbasisdxvec, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)

      t_start_tmp=ftimer()
      DO perm=1,EdgePerm
        CALL H1Basis_BrickEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              nbasisvec, BasisVec, EdgeDir(:,:,perm))
        CALL H1Basis_dBrickEdgeP(ngp, UWrk, VWrk, WWrk, EdgeP, &
              ndbasisdxvec, dBasisdxVec, EdgeDir(:,:,perm))
      END DO
      t_totvec_e=t_totvec_e+(ftimer()-t_start_tmp)

      t_start_tmp=ftimer()
      DO perm=1,FacePerm
        CALL H1Basis_BrickFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              nbasisvec, BasisVec, FaceDir(:,:,perm))
        CALL H1Basis_dBrickFaceP(ngp, UWrk, VWrk, WWrk, FaceP, &
              ndbasisdxvec, dBasisdxVec, FaceDir(:,:,perm))
      END DO
      t_totvec_f=t_totvec_f+(ftimer()-t_start_tmp)

      t_start_tmp=ftimer()
      CALL H1Basis_BrickBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              nbasisvec, BasisVec)
      CALL H1Basis_dBrickBubbleP(ngp, UWrk, VWrk, WWrk, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()

    CALL PrintTestData(Element, ngp, nrep, 12*nedof+6*nfdof+nbdof, &
            t_tot_n, t_tot_e+t_tot_f+t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_e+t_totvec_f+t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec, UWrk, VWrk, WWrk)
  END FUNCTION TestBrickElement

  FUNCTION TestBasis(ngp, nbasis, ndim, Basis1, Basis2, dBasisdx1, dBasisdx2, tol) RESULT(nerror)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ngp, nbasis, ndim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: Basis1(:,:), Basis2(:,:), &
            dBasisdx1(:,:,:), dBasisdx2(:,:,:)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror

    INTEGER :: i, j, dim

    WRITE (*,'(A)') 'TestBasis: Testing basis functions versus reference implementation.'
    WRITE (*,'(3(A,I0))') 'TestBasis: ngp=', ngp, ', nbasis=', nbasis, &
            ', ndim=', ndim

    nerror = 0
    ! Test basis
    DO j=1,nbasis
      DO i=1,ngp
        IF (ABS(Basis1(i,j)-Basis2(i,j)) >= tol) THEN
          nerror = nerror + 1
          WRITE (*,*) 'Basis:', i,j,Basis1(i,j), Basis2(i,j)
        END IF
      END DO
    END DO
    ! Test derivatives
    DO dim=1, ndim
      DO j=1,nbasis
        DO i=1,ngp
          IF (ABS(dBasisdx1(i,j,dim)-dBasisdx2(i,j,dim)) >= tol) THEN
            nerror = nerror + 1
            WRITE (*,*) 'dBasisdx:', i,j,dim, dBasisdx1(i,j,dim), dBasisdx2(i,j,dim)
          END IF
        END DO
      END DO
    END DO

    IF (nerror == 0) THEN
      WRITE (*,'(A,ES12.3)') 'TestBasis: Passed without errors, tol=', tol
    ELSE
      WRITE (*,'(A,ES12.3)') 'TestBasis: Failed with errors, tol=', tol
    END IF
  END FUNCTION TestBasis

  SUBROUTINE PrintTestData(Element, ngp, nrep, nndofs, t_n1, t_b1, t_tot1, &
          t_n2, t_b2, t_tot2)
  
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER, INTENT(IN) :: ngp, nrep, nndofs
    REAL(kind=dp), INTENT(IN) :: t_n1, t_b1, t_tot1, t_n2, t_b2, t_tot2

    WRITE (*,'(A,I0)') 'Element type=', Element % TYPE % ElementCode
    IF (ASSOCIATED(Element % PDefs)) THEN
      WRITE (*,'(A,I0)') 'Element polynomial degree=', Element % PDefs % P
    END IF
    WRITE (*,'(A,I0)') 'Element number of nodes=', Element % TYPE % NumberOfNodes
    WRITE (*,'(A,I0)') 'Element number of nonnodal dofs=', nndofs
    WRITE (*,'(A,I0)') 'Element number of bubble dofs=', GetElementNOFBDOFs()
    WRITE (*,'(A,I0)') 'Number of Gauss points=', ngp
    WRITE (*,'(A,I0)') 'Number of repetitions=', nrep
    WRITE (*,'(A,3F12.9)') 'Nodal  basis, nodal/bubble/total t(s):', t_n1, &
            t_b1, t_tot1
    WRITE (*,'(A,3F12.9)') 'Vector basis, nodal/bubble/total t(s):', t_n2, &
            t_b2, t_tot2

  END SUBROUTINE PrintTestData
  
  FUNCTION AllocatePElement(Mesh, ElementCode, P) RESULT(PElement)
    IMPLICIT NONE
    
    TYPE(Mesh_t) :: Mesh
    INTEGER, INTENT(IN) :: ElementCode, P
    TYPE(Element_t), POINTER :: PElement


    ! Construct P element
    PElement => AllocateElement()
    PElement % Type => GetElementType( ElementCode )
    CALL AllocatePDefinitions(PElement)

    PElement % BDofs = GetBubbleDofs( PElement, P )
    PElement % PDefs % P = P
    PElement % PDefs % GaussPoints = GetNumberOfGaussPoints(PElement, &
            Mesh)
  END FUNCTION AllocatePElement

  SUBROUTINE DeallocatePElement(PElement)
    IMPLICIT NONE
    
    TYPE(Element_t), POINTER :: PElement

    DEALLOCATE(PElement % PDefs)
    DEALLOCATE(PElement)
  END SUBROUTINE DeallocatePElement

  ! Portable wall-clock timer
  FUNCTION ftimer() RESULT(timerval)
    IMPLICIT NONE
    
    REAL(KIND=dp) :: timerval
    INTEGER(KIND=8) :: t, rate
    
#ifdef _OPENMP
    timerval = OMP_GET_WTIME()
#else
    CALL SYSTEM_CLOCK(t,count_rate=rate)
    timerval = REAL(t,dp)/rate
#endif
  END FUNCTION ftimer

!------------------------------------------------------------------------------
END SUBROUTINE H1BasisEvaluation
!------------------------------------------------------------------------------
