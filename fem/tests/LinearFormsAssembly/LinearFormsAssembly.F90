SUBROUTINE LinearFormsAssembly( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Unit test for linear form / vectorized basis function computation
!  routines in Elmer
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
    USE LinearForms
    USE ISO_C_BINDING
#ifdef _OPENMP
    USE omp_lib
#endif
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
    REAL(KIND=dp), PARAMETER :: tol1d = 1D-12, tol2d=1D-12, tol3d=1D-12
    REAL(KIND=dp) :: float_P
    INTEGER :: nerror, netest, P
    LOGICAL :: Found
    
    nerror = 0
    float_P = ListGetCReal(GetSolverParams(), 'P', Found)
    IF (.NOT. Found) float_P = 6.0_dp
    P = NINT(float_P)
    
    ! 1D tests
    netest = TestLineElement(Solver, P, tol1d)
    IF (netest /= 0) THEN
       CALL Warn('LinearFormsAssembly','Line element contained errors')
    END IF
    nerror = nerror + netest
    
    ! ! 2D tests 
    netest = TestTriangleElement(Solver, P, tol2d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Triangle element contained errors')
    END IF
    nerror = nerror + netest
    
    netest = TestQuadElement(Solver, P, tol2d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Quad element contained errors')
    END IF
    nerror = nerror + netest

    ! ! 3D tests
    netest = TestTetraElement(Solver, P, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Tetra element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestWedgeElement(Solver, P, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Wedge element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestBrickElement(Solver, P, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Brick element contained errors')
    END IF
    nerror = nerror + netest

    ! Build solution norm for error checking
    Solver % Variable % Norm = REAL(1+nerror,dp)
    Solver % Variable % Values = REAL(1+nerror,dp)
    
CONTAINS

  FUNCTION TestLineElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 202, P, tol)
  END FUNCTION TestLineElement
  
  FUNCTION TestTriangleElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 303, P, tol)
  END FUNCTION TestTriangleElement

  FUNCTION TestQuadElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 404, P, tol)
  END FUNCTION TestQuadElement

  FUNCTION TestTetraElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 504, P, tol)
  END FUNCTION TestTetraElement

  FUNCTION TestWedgeElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 706, P, tol)
  END FUNCTION TestWedgeElement

  FUNCTION TestBrickElement(Solver, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 808, P, tol)
  END FUNCTION TestBrickElement
  
  FUNCTION TestElement(Solver, ecode, P, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: ecode
    INTEGER, INTENT(IN) :: P
    REAL(kind=dp), INTENT(IN) :: tol

    TYPE(Element_t), POINTER :: Element, SingleElement
    TYPE(Mesh_t), POINTER :: NewMesh, OldMesh
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
            STIFFvec(:,:), FORCEvec(:)

    INTEGER :: i, j, k, l, q, nerror, nbasis, nndof, allocstat, tag, nthr, &
            nbasisvec, ndbasisdxvec, rep, dim, lm_eval, lm_eval_vec, NumGP

    INTEGER, PARAMETER :: NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_tot, t_startvec, t_endvec, t_tot_vec
    TYPE(GaussIntegrationPoints_t) :: Quadrature
    
    nerror = 0
    lm_eval = 0
    lm_eval_vec = 0
    t_tot = REAL(0,dp)
    t_tot_vec = REAL(0,dp)

    ! Create a mesh with a single element
    OldMesh => Solver % Mesh
    NewMesh => NULL()
    CALL AllocateMeshAndPElement(NewMesh, ecode, P, SingleElement)
    Solver % Mesh => NewMesh

    ! Insert P element definitions to Solver mapping (sets P elements as "active")
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = ecode / 100
      Solver % Def_Dofs(tag,1,6) = P
    END IF
    IF (ecode==404) THEN
      Quadrature = GaussPointsQuad((P+1)**2)
    ELSE
      Quadrature = GaussPoints(SingleElement, PReferenceElement=.TRUE.)
    END IF
    NumGP = Quadrature % N
    
    !$OMP PARALLEL SHARED(Solver, SingleElement, ecode, tol) &
    !$OMP PRIVATE(STIFF, FORCE, STIFFvec, FORCEvec, &
    !$OMP         LOAD, Element, nndof, nbasis, &
    !$OMP         allocstat, rep, t_start, t_end, &
    !$OMP         t_startvec, t_endvec) &
    !$OMP REDUCTION(+:nerror,t_tot,t_tot_vec,lm_eval,lm_eval_vec) &
    !$OMP DEFAULT(NONE)
    
    ! Construct a temporary element based on the one created previously
    Element => ClonePElement(SingleElement)

    nndof = Element % Type % NumberOfNodes
    nbasis = nndof + GetElementNOFDOFs( Element )
    
    ! Reserve workspace
    ALLOCATE(STIFF(nbasis, nbasis), FORCE(nbasis), &
            STIFFvec(nbasis, nbasis), FORCEvec(nbasis), &
            LOAD(nndof), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('LinearForms',&
              'Storage allocation for local matrices failed')
    END IF
    
    ! Initialize artificial load vector
    LOAD = REAL(1,dp)
    
    ! Warmup
    CALL LocalMatrix( STIFF, FORCE, LOAD, Element, nndof, nbasis)
    !$OMP BARRIER
    t_start = ftimer()
    DO rep=1,NREP
      ! Construct local matrix
!DIR$ NOINLINE
      CALL LocalMatrix( STIFF, FORCE, LOAD, Element, nndof, nbasis)
    END DO
    t_end = ftimer()
    lm_eval = NREP
    
    ! Warmup
    CALL LocalMatrixVec( STIFFvec, FORCEvec, LOAD, Element, nndof, nbasis)
    !$OMP BARRIER
    t_startvec = ftimer()
    DO rep=1,NREP
      ! Construct local matrix
!DIR$ NOINLINE
      CALL LocalMatrixVec( STIFFvec, FORCEvec, LOAD, Element, nndof, nbasis)
    END DO
    t_endvec = ftimer()
    lm_eval_vec = NREP

    nerror = TestLocalMatrix(nbasis, STIFF, STIFFvec, tol)
    nerror = nerror + TestLocalForce(nbasis, FORCE, FORCEvec, tol)

    t_tot = t_end - t_start
    t_tot_vec = t_endvec - t_startvec

    CALL DeallocatePElement(Element)
    DEALLOCATE(STIFF, FORCE, STIFFvec, FORCEvec, LOAD)

    !$OMP END PARALLEL

    ! Normalize the times / thread
    nthr = 1
    !$ nthr = omp_get_max_threads()
    t_tot = t_tot / nthr
    t_tot_vec = t_tot_vec / nthr
    CALL PrintTestData(SingleElement, NumGP, t_tot, lm_eval, &
            t_tot_vec, lm_eval_vec)

    CALL DeallocateTemporaryMesh(NewMesh)
    Solver % Mesh => OldMesh
    CALL DeallocatePElement(SingleElement)

    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = ecode / 100
      Solver % Def_Dofs(tag,1,6) = 0
    END IF
  END FUNCTION TestElement
  
  SUBROUTINE LocalMatrix( STIFF, FORCE, LOAD, Element, n, nd )
    IMPLICIT NONE
    REAL(KIND=dp) CONTIG :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) CONTIG, INTENT(IN) :: LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP,Weight
    LOGICAL :: Stat
    INTEGER :: i,j,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetReferenceElementNodes( Nodes, Element )
    STIFF = 0.0d0
    FORCE = 0.0d0

    dim = Element % Type % Dimension
    !Numerical integration:
    !----------------------
    IF (Element % TYPE % ElementCode / 100 == 4) THEN
      IP = GaussPointsQuad((Element % PDefs % P+1)**2)
    ELSE
      IP = GaussPoints( Element, PReferenceElement=.TRUE. )
    END IF

    DO t=1,IP % n

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                           IP % W(t), detJ, Basis, dBasisdx)
      
      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      Weight = IP % s(t) * DetJ
      ! STIFF=STIFF+(grad u, grad v)
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
            MATMUL( dBasisdx(1:nd,1:dim), TRANSPOSE( dBasisdx(1:nd,1:dim) ) )
      
      DO p=1,nd
        DO q=1,nd
          ! STIFF=STIFF+(grad u,v)
          ! -----------------------------------
          STIFF (p,q) = STIFF(p,q) + Weight * SUM(dBasisdx(q,1:dim)) * Basis(p)

          ! STIFF=STIFF+(u,v)
          STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
        END DO
      END DO
      
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
  END SUBROUTINE LocalMatrix

  SUBROUTINE LocalMatrixVec( STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) CONTIG :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) CONTIG, INTENT(IN) :: LOAD(:)
    INTEGER, INTENT(IN) :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    LOGICAL :: Stat
    INTEGER :: i, j, ngp, allocstat, gp
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:), dBasisdx(:,:,:), &
            DetJ(:), LoadAtIPs(:)
    !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, DetJ, LoadAtIPs)
!------------------------------------------------------------------------------
    CALL GetReferenceElementNodes( Nodes, Element )

    STIFF = REAL(0, dp)
    FORCE = REAL(0, dp)

    ! Get integration points
    IF (Element % TYPE % ElementCode / 100 == 4) THEN
      IP = GaussPointsQuad((Element % PDefs % P+1)**2)
    ELSE
      IP = GaussPoints( Element, PReferenceElement=.TRUE. )
    END IF
    ngp = IP % n

    ! Reserve workspace
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), &
              DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LocalMatrixVec',&
                   'Storage allocation for local element basis failed')
      END IF
    ELSE IF (SIZE(Basis,1) /= ngp .OR. SIZE(Basis,2) /= nd) THEN
      DEALLOCATE(Basis, dBasisdx, DetJ, LoadAtIPs)
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), &
              DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('LocalMatrixVec',&
                'Storage allocation for local element basis failed')
      END IF
    END IF
    
    ! Compute values of all basis functions at all integration points
    stat = ElementInfoVec( Element, Nodes, ngp, &
            IP % U, IP % V, IP % W, DetJ, SIZE(Basis,2), Basis, dBasisdx )

        
    ! Compute actual integration weights (recycle memory space of DetJ)
    DO i=1,ngp
       DetJ(i) = Ip % s(i)*Detj(i)
    END DO

    ! STIFF=STIFF+(grad u, grad u)
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, &
             dBasisdx, DetJ, STIFF)
    ! STIFF=STIFF+(u,u)
    CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, STIFF)
    ! STIFF=STIFF+(grad u,v)
    CALL LinearForms_GradUdotU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, Basis, DetJ, STIFF)
    
    ! Source terms at IPs
    !------------------------------------------
    ! LoadAtIPs(1:ngp) = MATMUL( Basis(1:ngp,1:n), LOAD(1:n) )
    CALL LinearForms_ProjectToU(ngp, n, Basis, LOAD, LoadAtIPs)

    ! FORCE=FORCE+(u,f)
    CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, LoadAtIPs, FORCE)
  END SUBROUTINE LocalMatrixVec

  FUNCTION TestLocalMatrix(nbasis, STIFF1, STIFF2, tol) RESULT(nerror)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbasis
    REAL(KIND=dp) CONTIG, INTENT(IN) :: STIFF1(:,:), STIFF2(:,:)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror

    INTEGER :: i, j, dim

    nerror = 0
    ! Test element of local matrix 
    DO j=1,nbasis
      DO i=1,nbasis
        IF (ABS(STIFF1(i,j)-STIFF2(i,j)) >= tol) THEN
          nerror = nerror + 1
          WRITE (*,*) 'STIFF:', i,j,STIFF1(i,j), STIFF2(i,j)
        END IF
      END DO
    END DO
  END FUNCTION TestLocalMatrix
  
  FUNCTION TestLocalForce(nbasis, FORCE1, FORCE2, tol) RESULT(nerror)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbasis
    REAL(KIND=dp) CONTIG, INTENT(IN) :: FORCE1(:), FORCE2(:)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror

    INTEGER :: i, dim

    nerror = 0
    ! Test element of local force vector
    DO i=1,nbasis
      IF (ABS(FORCE1(i)-FORCE2(i)) >= tol) THEN
        nerror = nerror + 1
        WRITE (*,*) 'FORCE:', i, FORCE1(i), FORCE2(i)
      END IF
    END DO
  END FUNCTION TestLocalForce

  SUBROUTINE GetReferenceElementNodes( ElementNodes, Element )
     TYPE(Nodes_t), TARGET :: ElementNodes
     TYPE(Element_t) :: Element

     INTEGER :: i, n, padn, sz, astat

     n = Element % TYPE % NumberOfNodes
     padn = n

     ! TODO: Implement padding
     IF (.NOT. ALLOCATED( ElementNodes % xyz)) THEN
       ! Deallocate old storage
       IF (ASSOCIATED(ElementNodes % x)) DEALLOCATE(ElementNodes % x) 
       IF (ASSOCIATED(ElementNodes % y)) DEALLOCATE(ElementNodes % y) 
       IF (ASSOCIATED(ElementNodes % z)) DEALLOCATE(ElementNodes % z) 

       ! Allocate new storage
       ALLOCATE(ElementNodes % xyz(padn,3))
       ElementNodes % xyz = REAL(0,dp)
       ElementNodes % x => ElementNodes % xyz(1:n,1)
       ElementNodes % y => ElementNodes % xyz(1:n,2)
       ElementNodes % z => ElementNodes % xyz(1:n,3)
     ELSE IF (SIZE(ElementNodes % xyz, 1)<padn) THEN
       DEALLOCATE(ElementNodes % xyz)
       ALLOCATE(ElementNodes % xyz(padn,3))
       ElementNodes % xyz = REAL(0,dp)
       ElementNodes % x => ElementNodes % xyz(1:n,1)
       ElementNodes % y => ElementNodes % xyz(1:n,2)
       ElementNodes % z => ElementNodes % xyz(1:n,3)
     ELSE
       ElementNodes % x => ElementNodes % xyz(1:n,1)
       ElementNodes % y => ElementNodes % xyz(1:n,2)
       ElementNodes % z => ElementNodes % xyz(1:n,3)
       sz = SIZE(ElementNodes % xyz,1)
       SELECT CASE(Element % TYPE % DIMENSION)
       CASE(1)
         DO i=n+1,sz
           ElementNodes % xyz(i,1) = REAL(0,dp)
         END DO
         DO i=1,sz
           ElementNodes % xyz(i,2) = REAL(0,dp)
           ElementNodes % xyz(i,3) = REAL(0,dp)
         END DO
       CASE(2)
         DO i=n+1,sz
           ElementNodes % xyz(i,1) = REAL(0,dp)
           ElementNodes % xyz(i,2) = REAL(0,dp)
         END DO
         DO i=1,sz
           ElementNodes % xyz(i,3) = REAL(0,dp)
         END DO
       CASE(3)
         DO i=n+1,sz
           ElementNodes % xyz(i,1) = REAL(0,dp)
           ElementNodes % xyz(i,2) = REAL(0,dp)
           ElementNodes % xyz(i,3) = REAL(0,dp)
         END DO
       CASE DEFAULT
         CALL Fatal('GetReferenceElementNodes','Unsupported element dimension')
       END SELECT
     END IF

     IF (isPElement(Element)) THEN
       CALL GetRefPElementNodes(Element, &
               ElementNodes % x, &
               ElementNodes % y, & 
               ElementNodes % z)
     ELSE
       IF (ASSOCIATED(Element % Type % NodeU)) THEN
         !DIR$ IVDEP
         DO i=1,n
           ElementNodes % x(i) = Element % Type % NodeU(i)
         END DO
       END IF
       IF (ASSOCIATED(Element % Type % NodeV)) THEN
         !DIR$ IVDEP
         DO i=1,n
           ElementNodes % y(i) = Element % Type % NodeV(i)
         END DO
       END IF
       IF (ASSOCIATED(Element % Type % NodeW)) THEN
         !DIR$ IVDEP
         DO i=1,n
           ElementNodes % z(i) = Element % Type % NodeW(i)
         END DO
       END IF
     END IF
  END SUBROUTINE GetReferenceElementNodes

  SUBROUTINE PrintTestData(Element, ngp, t_n1, evals1, t_n2, evals2)
  
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER, INTENT(IN) :: ngp, evals1, evals2
    REAL(kind=dp), INTENT(IN) :: t_n1, t_n2

    WRITE (*,'(A,I0)') 'Element type=', Element % TYPE % ElementCode
    IF (ASSOCIATED(Element % PDefs)) THEN
      WRITE (*,'(A,I0)') 'Element polynomial degree=', Element % PDefs % P
    END IF
    WRITE (*,'(A,L1)') 'Active P element=', isActivePElement(Element)
    WRITE (*,'(A,I0)') 'Element number of nodes=', Element % TYPE % NumberOfNodes
    WRITE (*,'(A,I0)') 'Element number of nonnodal dofs=', GetElementNOFDOFs(Element)
    WRITE (*,'(A,I0)') 'Element number of bubble dofs=', GetElementNOFBDOFs()
    WRITE (*,'(A,I0)') 'Number of Gauss points=', ngp
    WRITE (*,'(A,I0)') 'Nodal basis, number of local matrix evaluations=', evals1
    WRITE (*,'(A,F12.9)') 'Nodal  basis, local matrix assembly t(s):', t_n1
    WRITE (*,'(A,F12.2)') 'Nodal  basis, local matrix evaluations/sec:', evals1/t_n1
    WRITE (*,'(A,I0)') 'Vector basis, number of local matrix evaluations=', evals2 
    WRITE (*,'(A,F12.9)') 'Vector basis, local matrix assembly t(s):', t_n2
    WRITE (*,'(A,F12.2)') 'Vector basis, local matrix evaluations/sec:', evals2/t_n2

  END SUBROUTINE PrintTestData
  
  SUBROUTINE AllocateMeshAndPElement(Mesh, ElementCode, P, PElement)
    IMPLICIT NONE
    
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, INTENT(IN) :: ElementCode, P
    TYPE(Element_t), POINTER :: PElement
    INTEGER :: node, edge, face, astat

    ! Allocate a mesh with a single element
    Mesh => AllocateMesh()
    ALLOCATE(Mesh % Elements(1), STAT=astat)
    IF (astat /= 0) THEN
      CALL Fatal('AllocateMeshAndElement','Allocation of mesh element failed')
    END IF
    Mesh % NumberOfBulkElements = 1

    ! Construct P element
    PElement => AllocateElement()
    PElement % ElementIndex = 1
    PElement % BodyId = 1
    PElement % Type => GetElementType( ElementCode )
    CALL AllocatePDefinitions(PElement)

    ! Add element node indexes to element
    ALLOCATE(PElement % NodeIndexes(PElement % Type % NumberOfNodes), STAT=astat)
    IF (astat /= 0) THEN
      CALL Fatal('AllocateMeshAndElement','Allocation of node indices failed')
    END IF
    PElement % NodeIndexes(:) = [(node, node=1,PElement % Type % NumberOfNodes)]

    IF (PElement % Type % Dimension > 1) THEN
      ! Add element edge indexes to element
      ALLOCATE(PElement % EdgeIndexes(PElement % Type % NumberOfEdges), STAT=astat)
      IF (astat /= 0) THEN
        CALL Fatal('AllocateMeshAndElement','Allocation of edge indices failed')
      END IF
      PElement % EdgeIndexes(:) = [(edge, edge=1,PElement % Type % NumberOfEdges)]
      
      ! Add all element edges to mesh
      ALLOCATE(Mesh % Edges(PElement % Type % NumberOfEdges), STAT=astat)
      IF (astat /= 0) THEN
        CALL Fatal('AllocateMeshAndElement','Allocation of mesh edges failed')
      END IF
      Mesh % NumberOfEdges = PElement % Type % NumberOfEdges
      ! Mesh % MinEdgeDofs = HUGE(Mesh % MinEdgeDofs)
      Mesh % MinEdgeDofs = 0
      Mesh % MaxEdgeDofs = 0
      DO edge=1,PElement % Type % NumberOfEdges
        CALL InitializePElement(Mesh % Edges(edge))
        ! Allocate edge element (always 202, even in 3D)
        Mesh % Edges(edge) % Type => GetElementType(202, .FALSE.)
        CALL AllocatePDefinitions(Mesh % Edges(edge))
        Mesh % Edges(edge) % PDefs % P = P
        Mesh % Edges(edge) % PDefs % isEdge = .TRUE.
        Mesh % Edges(edge) % PDefs % GaussPoints = (P+1) ** Mesh % Edges(edge) % Type % DIMENSION
        Mesh % Edges(edge) % PDefs % LocalNumber = edge
        Mesh % Edges(edge) % BDOFs = GetBubbleDofs(Mesh % Edges(edge), P)
        Mesh % MinEdgeDofs = MIN(Mesh % MinEdgeDofs, Mesh % Edges(edge) % BDOFs)
        Mesh % MaxEdgeDofs = MAX(Mesh % MaxEdgeDofs, Mesh % Edges(edge) % BDOFs)
      END DO
    END IF

    IF (PElement % Type % Dimension > 2) THEN
      ! Add element face indexes to element
      ALLOCATE(PElement % FaceIndexes(PElement % Type % NumberOfFaces))
      IF (astat /= 0) THEN
        CALL Fatal('AllocateMeshAndElement','Allocation of face indices failed')
      END IF
      PElement % FaceIndexes(:) = [(face, face=1,PElement % Type % NumberOfFaces)]
      
      ! Add element faces to mesh
      ALLOCATE(Mesh % Faces(PElement % Type % NumberOfFaces), STAT=astat)
      IF (astat /= 0) THEN
        CALL Fatal('AllocateMeshAndElement','Allocation of mesh faces failed')
      END IF
      Mesh % NumberOfFaces = PElement % Type % NumberOfFaces
      Mesh % MinFaceDofs = HUGE(Mesh % MinFaceDofs)
      Mesh % MaxFaceDofs = 0
      DO face=1,PElement % Type % NumberOfFaces
        CALL InitializePElement(Mesh % Faces(face))
        ! Allocate edge element (always 303 or 404, depending on element)
        SELECT CASE (ElementCode/100)
        CASE(5)
          Mesh % Faces(face) % Type => GetElementType(303, .FALSE.)
        CASE(6)
          IF ( face == 1 ) THEN
            Mesh % Faces(Face) % Type => GetElementType( 404, .FALSE. )
          ELSE
            Mesh % Faces(Face) % Type => GetElementType( 303, .FALSE. )
          END IF
        CASE(7)
          IF ( face <= 2 ) THEN
            Mesh % Faces(Face) % Type => GetElementType( 303, .FALSE. )
          ELSE
            Mesh % Faces(Face) % Type => GetElementType( 404, .FALSE. )
          END IF
        CASE(8)
          Mesh % Faces(Face) % Type => GetElementType( 404, .FALSE.)
        CASE DEFAULT
          CALL Fatal('AllocateMeshAndElement','Unknown element type')
        END SELECT
        CALL AllocatePDefinitions(Mesh % Faces(face))
        Mesh % Faces(face) % PDefs % P = P
        Mesh % Faces(face) % PDefs % isEdge = .TRUE.
        Mesh % Faces(face) % PDefs % GaussPoints = (P+1) ** Mesh % Faces(face) % Type % DIMENSION
        Mesh % Faces(face) % PDefs % LocalNumber = face
        Mesh % Faces(face) % BDOFs = GetBubbleDofs(Mesh % Faces(face), P)
        Mesh % MinFaceDofs = MAX(Mesh % MinFaceDofs, Mesh % Faces(face) % BDOFs)
        Mesh % MaxFaceDofs = MAX(Mesh % MaxFaceDofs, Mesh % Faces(face) % BDOFs)
      END DO
    END IF 

    PElement % BDofs = GetBubbleDofs( PElement, P )
    PElement % PDefs % P = P
    IF (ElementCode == 504) THEN
      PElement % PDefs % TetraType = 1
    ELSE
      PElement % PDefs % TetraType = 0
    END IF
    PElement % PDefs % isEdge = .FALSE.
    ! Set max element dofs to mesh
    Mesh % MaxElementDOFs = PElement % Type % NumberOfNodes + &
           PElement % Type % NumberOfEdges * Mesh % MaxEdgeDOFs + &
           PElement % Type % NumberOfFaces * Mesh % MaxFaceDOFs + &
           PElement % BDOFs

    PElement % PDefs % GaussPoints = (P+1) ** PElement % Type % DIMENSION
  END SUBROUTINE AllocateMeshAndPElement

  SUBROUTINE InitializePElement(Element)
    IMPLICIT NONE
    TYPE(Element_t) :: Element

    Element % BDOFs    =  0
    Element % NDOFs    =  0
    Element % BodyId   = 1
    Element % Splitted =  0
    Element % hK = 0
    Element % ElementIndex = 0
    Element % StabilizationMk = 0
    NULLIFY( Element % TYPE )
    NULLIFY( Element % PDefs )
    NULLIFY( Element % BubbleIndexes )
    NULLIFY( Element % DGIndexes )
    NULLIFY( Element % NodeIndexes )
    NULLIFY( Element % EdgeIndexes )
    NULLIFY( Element % FaceIndexes )
    NULLIFY( Element % BoundaryInfo )
  END SUBROUTINE InitializePElement

  FUNCTION ClonePElement(Element) RESULT(ClonedElement)
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element, ClonedElement
    
    ALLOCATE(ClonedElement)
    CALL InitializePElement(ClonedElement)
    ClonedElement % BDOFs = Element % BDOFs
    ClonedElement % NDOFs =  Element % NDOFs
    ClonedElement % BodyId = Element % BodyId
    ClonedElement % Splitted = Element % Splitted 
    ClonedElement % hK = Element % hK
    ClonedElement % ElementIndex = Element % ElementIndex
    ClonedElement % StabilizationMk = Element % StabilizationMk
    ClonedElement % Type => Element % Type
    IF (ASSOCIATED(Element % PDefs)) THEN
      CALL AllocatePDefinitions(ClonedElement)
      ClonedElement % PDefs % P = Element % PDefs % P
      ClonedElement % PDefs % TetraType = Element % PDefs % TetraType
      ClonedElement % PDefs % isEdge = Element % PDefs % isEdge
      ClonedElement % PDefs % GaussPoints = Element % PDefs % GaussPoints
      ClonedElement % PDefs % pyramidQuadEdge = Element % PDefs % pyramidQuadEdge
      ClonedElement % PDefs % localNumber = Element % PDefs % localNumber
    END IF
    IF (ASSOCIATED( Element % NodeIndexes )) THEN
      ALLOCATE(ClonedElement % NodeIndexes( SIZE(Element % NodeIndexes)))
      ClonedElement % NodeIndexes(:) = Element % NodeIndexes(:)
    END IF
    IF (ASSOCIATED( Element % EdgeIndexes )) THEN
      ALLOCATE(ClonedElement % EdgeIndexes( SIZE(Element % EdgeIndexes)))
      ClonedElement % EdgeIndexes(:) = Element % EdgeIndexes(:)
    END IF
    IF (ASSOCIATED( Element % FaceIndexes )) THEN
      ALLOCATE(ClonedElement % FaceIndexes( SIZE(Element % FaceIndexes)))
      ClonedElement % FaceIndexes(:) = Element % FaceIndexes(:)
    END IF
  END FUNCTION ClonePElement
  
  SUBROUTINE DeallocatePElement(PElement)
    IMPLICIT NONE
    
    TYPE(Element_t), POINTER :: PElement

    IF (ASSOCIATED(PElement % PDefs)) DEALLOCATE(PElement % PDefs)
    IF (ASSOCIATED(PElement % NodeIndexes)) DEALLOCATE(PElement % NodeIndexes)
    IF (ASSOCIATED(PElement % EdgeIndexes)) DEALLOCATE(PElement % EdgeIndexes)
    IF (ASSOCIATED(PElement % FaceIndexes)) DEALLOCATE(PElement % FaceIndexes)
    DEALLOCATE(PElement)
  END SUBROUTINE DeallocatePElement

  SUBROUTINE DeallocateTemporaryMesh(Mesh)
    IMPLICIT NONE
 
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: edge, face

     IF (ASSOCIATED(Mesh % Elements)) DEALLOCATE(Mesh % Elements)
     DO edge=1,Mesh % NumberOfEdges
       IF (ASSOCIATED(Mesh % Edges(edge) % PDefs)) DEALLOCATE(Mesh % Edges(edge) % PDefs)
     END DO
     IF (ASSOCIATED(Mesh % Edges)) DEALLOCATE(Mesh % Edges)
     DO face=1,Mesh % NumberOfFaces
       IF (ASSOCIATED(Mesh % Faces(face) % PDefs)) DEALLOCATE(Mesh % Faces(face) % PDefs)
     END DO
     IF (ASSOCIATED(Mesh % Faces)) DEALLOCATE(Mesh % Faces)

     DEALLOCATE( Mesh % Nodes )
     DEALLOCATE( Mesh ) 
  END SUBROUTINE DeallocateTemporaryMesh
  
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
END SUBROUTINE LinearFormsAssembly
!------------------------------------------------------------------------------
