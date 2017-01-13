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
    REAL(kind=dp), PARAMETER :: tol1d = 1D-12, tol2d=1D-12, tol3d=1D-12
    INTEGER :: nerror, netest

    nerror = 0

    ! 1D tests
    netest = TestLineElement(Solver, tol1d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Line element contained errors')
    END IF
    nerror = nerror + netest
    
    ! 2D tests 
    netest = TestTriangleElement(Solver, tol2d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Triangle element contained errors')
    END IF
    nerror = nerror + netest
    
    netest = TestQuadElement(Solver, tol2d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Quad element contained errors')
    END IF
    nerror = nerror + netest

    ! 3D tests
    netest = TestTetraElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Tetra element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestWedgeElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Wedge element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestBrickElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('LinearFormsAssembly','Brick element contained errors')
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
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 202, tol)
  END FUNCTION TestLineElement
  
  FUNCTION TestTriangleElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 303, tol)
  END FUNCTION TestTriangleElement

  FUNCTION TestQuadElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 404, tol)
  END FUNCTION TestQuadElement

  FUNCTION TestTetraElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 504, tol)
  END FUNCTION TestTetraElement

  FUNCTION TestWedgeElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 706, tol)
  END FUNCTION TestWedgeElement

  FUNCTION TestBrickElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror
    
    nerror = TestElement(Solver, 808, tol)
  END FUNCTION TestBrickElement
  
  FUNCTION TestElement(Solver, ecode, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: ecode
    REAL(kind=dp), INTENT(IN) :: tol

    TYPE(Element_t), POINTER :: Element    
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:), &
            STIFFvec(:,:), FORCEvec(:)

    INTEGER :: i, j, k, l, q, nerror, nbasis, nndof, nbdof, allocstat, tag, &
            nbasisvec, ndbasisdxvec, rep, dim, lm_eval, lm_eval_vec, NumGP

    INTEGER, PARAMETER :: P = 1, NREP = 100, NumGP1D = 8
    REAL(kind=dp) :: t_start, t_end, t_tot, t_startvec, t_endvec, t_tot_vec
    
    nerror = 0
    lm_eval = 0
    lm_eval_vec = 0
    t_tot = REAL(0,dp)
    t_tot_vec = REAL(0,dp)

    ! Insert P element definitions to Solver mapping
    IF (ALLOCATED(Solver % Def_Dofs)) THEN
      tag = ecode / 100
      Solver % Def_Dofs(tag,1,6) = 1
    END IF

    !$OMP PARALLEL SHARED(Solver, ecode, tol) &
    !$OMP PRIVATE(STIFF, FORCE, STIFFvec, FORCEvec, &
    !$OMP         LOAD, Element, nndof, nbdof, nbasis, &
    !$OMP         NumGP, allocstat, rep, t_start, t_end, &
    !$OMP         t_startvec, t_endvec) &
    !$OMP REDUCTION(+:nerror,t_tot,t_tot_vec,lm_eval,lm_eval_vec) &
    !$OMP DEFAULT(NONE)

    Element => AllocatePElement(Solver, ecode, P)

    NumGP = NumGP1D ** Element % Type % Dimension
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 

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
    CALL LocalMatrix( STIFF, FORCE, LOAD, Element, nndof, nbasis, NumGP)
    !$OMP BARRIER
    t_start = ftimer()
    DO rep=1,NREP
      ! Construct local matrix
!DIR$ NOINLINE
      CALL LocalMatrix( STIFF, FORCE, LOAD, Element, nndof, nbasis, NumGP)
    END DO
    t_end = ftimer()
    lm_eval = NREP
    
    ! Warmup
    CALL LocalMatrixVec( STIFFvec, FORCEvec, LOAD, Element, nndof, nbasis, NumGP)
    !$OMP BARRIER
    t_startvec = ftimer()
    DO rep=1,NREP
      ! Construct local matrix
!DIR$ NOINLINE
      CALL LocalMatrixVec( STIFFvec, FORCEvec, LOAD, Element, nndof, nbasis, NumGP)
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
    Element => AllocatePElement(Solver, ecode, P)
    NumGP = NumGP1D ** Element % Type % Dimension
    CALL PrintTestData(Element, NumGP, t_tot, lm_eval, &
            t_tot_vec, lm_eval_vec)
    CALL DeallocatePElement(Element)
  END FUNCTION TestElement
  
  SUBROUTINE LocalMatrix( STIFF, FORCE, LOAD, Element, n, nd, NumGP )
    IMPLICIT NONE
    REAL(KIND=dp) CONTIG :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) CONTIG, INTENT(IN) :: LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    INTEGER, OPTIONAL :: NumGP
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    CALL GetReferenceElementNodes( Nodes, Element )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IF (PRESENT(NumGP)) THEN
      IP = GaussPoints( Element, NumGP )
    ELSE
      IP = GaussPoints( Element )
    END IF

    DO t=1,IP % n

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                           IP % W(t), detJ, Basis, dBasisdx)
      
      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
                MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO    
  END SUBROUTINE LocalMatrix

  SUBROUTINE LocalMatrixVec( STIFF, FORCE, LOAD, Element, n, nd, NumGP )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) CONTIG :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) CONTIG, INTENT(IN) :: LOAD(:)
    INTEGER, INTENT(IN) :: n, nd
    TYPE(Element_t), POINTER :: Element
    INTEGER, OPTIONAL :: NumGP
!------------------------------------------------------------------------------
    LOGICAL :: Stat
    INTEGER :: i, ngp, allocstat
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
    IF (PRESENT(NumGP)) THEN
      IP = GaussPoints( Element, NumGP )
    ELSE
      IP = GaussPoints( Element )
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
    ELSE IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) THEN
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
            IP % U, IP % V, IP % W, DetJ, Basis, dBasisdx )
    
    ! Compute actual integration weights (recycle memory space of DetJ)
    DO i=1,ngp
       DetJ(i) = Ip % s(i)*Detj(i)
    END DO

    ! STIFF=STIFF+(grad u, grad u)
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, &
            dBasisdx, DetJ, STIFF)

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
     END IF
     
     IF (isPElement(Element)) THEN
       CALL GetRefPElementNodes(Element, &
               ElementNodes % x, &
               ElementNodes % y, & 
               ElementNodes % z)
     ELSE
       IF (ASSOCIATED(Element % Type % NodeU)) THEN
         ElementNodes % x(1:n) = Element % Type % NodeU(1:n)
       END IF
       IF (ASSOCIATED(Element % Type % NodeV)) THEN
         ElementNodes % y(1:n) = Element % Type % NodeV(1:n)
       END IF
       IF (ASSOCIATED(Element % Type % NodeW)) THEN
         ElementNodes % z(1:n) = Element % Type % NodeW(1:n)
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
  
  FUNCTION AllocatePElement(Solver, ElementCode, P) RESULT(PElement)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: ElementCode, P
    TYPE(Element_t), POINTER :: PElement
    
    ! Construct P element
    PElement => AllocateElement()
    PElement % Type => GetElementType( ElementCode )
    CALL AllocatePDefinitions(PElement)
    PElement % BDofs = GetBubbleDofs( PElement, P )
    PElement % PDefs % P = P
    PElement % PDefs % GaussPoints = GetNumberOfGaussPoints(PElement, &
            Solver % Mesh)
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
    INTEGER :: t, rate
    
#ifdef _OPENMP
    timerval = OMP_GET_WTIME()
#else
    CALL SYSTEM_CLOCK(t,count_rate=rate)
    timerval = REAL(t,KIND(dp))/REAL(rate,KIND(dp))
#endif
  END FUNCTION ftimer

!------------------------------------------------------------------------------
END SUBROUTINE LinearFormsAssembly
!------------------------------------------------------------------------------
