SUBROUTINE H1Basis( Model,Solver,dt,TransientSimulation )
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
      CALL Warn('H1Basis','Line element contained errors')
    END IF
    nerror = nerror + netest
    
    ! 2D tests 
    netest = TestTriangleElement(Solver, tol1d)
    IF (netest /= 0) THEN
      CALL Warn('H1Basis','Triangle element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestQuadElement(Solver, tol1d)
    IF (netest /= 0) THEN
      CALL Warn('H1Basis','Quad element contained errors')
    END IF
    nerror = nerror + netest

    ! 3D tests
    netest = TestTetraElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1Basis','Tetra element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestWedgeElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1Basis','Wedge element contained errors')
    END IF
    nerror = nerror + netest

    netest = TestBrickElement(Solver, tol3d)
    IF (netest /= 0) THEN
      CALL Warn('H1Basis','Brick element contained errors')
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
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 6, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 202, P)
    GP = GaussPoints(Element)

    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions1D( Basis(i,1:nndof), Element, GP % U(i))
        CALL NodalFirstDerivatives1D( dBasisdx(i,1:nndof,1:3), Element, GP % U(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      DO j=1, Element % BDOFs
        DO i=1,ngp
          Basis(i,j+2) = LineBubblePBasis(j+1, GP % U(i),.FALSE.)
          dBasisdx(i,j+2,1) = dLineBubblePBasis(j+1, GP % U(i), .FALSE.)
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
      t_start_tmp=ftimer()
      CALL H1LineNodalBasisVec(ngp, GP % U, BasisVec)
      CALL H1dLineNodalBasisVec(ngp, GP % U, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 2
      t_start_tmp=ftimer()
      CALL H1LineBubblePBasisVec(ngp, GP % U, P, nbasisvec, BasisVec, .FALSE.)
      ndbasisdxvec = 2
      CALL H1dLineBubblePBasisVec(ngp, GP % U, P, ndbasisdxvec, dBasisdxVec, &
              .FALSE.)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()

    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestLineElement

  FUNCTION TestTriangleElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 6, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 303, P)
    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = TriangleNodalPBasis(ndof, GP % U(i), GP % V(i))
          dBasisdx(i,ndof,1:2) = dTriangleNodalPBasis(ndof, GP % U(i), GP % V(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      q = 3
      DO i = 0,p-3
        DO j = 0,p-i-3
          q = q + 1
          DO k = 1, ngp
            Basis(k, q) = TriangleBubblePBasis(i,j,GP % U(k), GP % v(k))
            dBasisdx(k, q, 1:2) = dTriangleBubblePBasis(i,j,GP % u(k), GP % v(k))  
          END DO
        END DO
        t_tot_b=t_tot_b+(ftimer()-t_start_tmp)
      END DO
    END DO
    t_end = ftimer()
    
    ! Initialize arrays
    BasisVec = 0
    dBasisdxVec = 0

    t_totvec_n = REAL(0,dp)
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      t_start_tmp=ftimer()
      CALL H1TriangleNodalPBasisVec(ngp, GP % U, GP % V, BasisVec)
      CALL H1dTriangleNodalPBasisVec(ngp, GP % U, GP % V, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 3
      t_start_tmp=ftimer()
      CALL H1TriangleBubblePBasisVec(ngp, GP % U, GP % V, P, &
              nbasisvec, BasisVec)
      ndbasisdxvec = 3
      CALL H1dTriangleBubblePBasisVec(ngp, GP % U, GP % V, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestTriangleElement
  
  FUNCTION TestQuadElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 6, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 404, P)
    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions2D( Basis(i,1:nndof), Element, GP % U(i), GP % V(i))
        CALL NodalFirstDerivatives2D( dBasisdx(i,1:nndof,1:3), Element, GP % U(i), GP % V(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)
    
      ! Bubble basis 
      t_start_tmp=ftimer()
      q=4
      DO i=2,(p-2)
        DO j=2,(p-i)
          q = q + 1
          DO k = 1, ngp
            Basis(k, q) = QuadBubblePBasis(i,j,GP % U(k),GP % V(k))
            dBasisdx(k, q, 1:2) = dQuadBubblePBasis(i,j,GP % U(k),GP % V(k))
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
      t_start_tmp=ftimer()
      CALL H1QuadNodalBasisVec(ngp, GP % U, GP % V, BasisVec)
      CALL H1dQuadNodalBasisVec(ngp, GP % U, GP % V, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 4
      t_start_tmp=ftimer()
      CALL H1QuadBubblePBasisVec(ngp, GP % U, GP % V, P, &
              nbasisvec, BasisVec)
      ndbasisdxvec = 4
      CALL H1dQuadBubblePBasisVec(ngp, GP % U, GP % V, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestQuadElement

  FUNCTION TestTetraElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 6, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 504, P)
    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = TetraNodalPBasis(ndof, GP % U(i), GP % V(i), GP % W(i))
          dBasisdx(i,ndof,1:3) = dTetraNodalPBasis(ndof, GP % U(i), GP % V(i), GP % W(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      q = 4
      DO i=0,p-4
        DO j=0,p-i-4
          DO k=0,p-i-j-4
            q = q + 1
            DO l = 1, ngp
              Basis(l, q) = TetraBubblePBasis(i,j,k,GP % U(l), GP % V(l), GP % W(l))
              dBasisdx(l, q, 1:3) = dTetraBubblePBasis(i,j,k,GP % U(l), GP % V(l), GP % W(l))  
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
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      t_start_tmp=ftimer()
      CALL H1TetraNodalPBasisVec(ngp, GP % U, GP % V, GP % W, BasisVec)
      CALL H1dTetraNodalPBasisVec(ngp, GP % U, GP % V, GP % W, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 4
      t_start_tmp=ftimer()
      CALL H1TetraBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              nbasisvec, BasisVec)
      ndbasisdxvec = 4
      CALL H1dTetraBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestTetraElement

  
  FUNCTION TestWedgeElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, k, l, q, ndof, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 5, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 706, P)
    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO ndof=1,nndof
        DO i=1,ngp
          Basis(i,ndof) = WedgeNodalPBasis(ndof, GP % U(i), GP % V(i), GP % W(i))
          dBasisdx(i,ndof,1:3) = dWedgeNodalPBasis(ndof, GP % U(i), GP % V(i), GP % W(i))
        END DO
      END Do
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      q = 6
      DO i=0,p-5
        DO j=0,p-5-i
          DO k=2,p-3-i-j
            q = q + 1
            DO l = 1, ngp
              Basis(l, q) = WedgeBubblePBasis(i,j,k,GP % U(l), GP % V(l), GP % W(l))
              dBasisdx(l, q, 1:3) = dWedgeBubblePBasis(i,j,k,GP % U(l), GP % V(l), GP % W(l))  
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
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      t_start_tmp=ftimer()
      CALL H1WedgeNodalPBasisVec(ngp, GP % U, GP % V, GP % W, BasisVec)
      CALL H1dWedgeNodalPBasisVec(ngp, GP % U, GP % V, GP % W, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 6
      t_start_tmp=ftimer()
      CALL H1WedgeBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              nbasisvec, BasisVec)
      ndbasisdxvec = 6
      CALL H1dWedgeBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestWedgeElement

  FUNCTION TestBrickElement(Solver, tol) RESULT(nerror)
    IMPLICIT NONE
    
    TYPE(Solver_t) :: Solver
    REAL(kind=dp), INTENT(IN) :: tol
    TYPE(Element_t), POINTER :: Element
    TYPE( GaussIntegrationPoints_t ) :: GP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:,:), dBasisdx(:,:,:), &
            BasisVec(:,:), dBasisdxVec(:,:,:)
!DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, BasisVec, dBasisdxVec

    INTEGER :: i, j, k, l, q, ngp, nerror, nbasis, nndof, nbdof, allocstat, &
            nbasisvec, ndbasisdxvec, rep, dim
    INTEGER, PARAMETER :: P = 7, NREP = 100
    REAL(kind=dp) :: t_start, t_end, t_startvec, t_endvec, &
            t_start_tmp, t_tot_n, t_totvec_n, &
            t_tot_b, t_totvec_b
    

    nerror = 0
    Element => AllocatePElement(Solver % Mesh, 808, P)
    GP = GaussPoints(Element)
    
    nndof = Element % Type % NumberOfNodes
    nbdof = Element % BDofs
    nbasis = nndof + nbdof 
    ngp = GP % N

    ! Reserve workspace for finite element basis
    ALLOCATE(Basis(ngp,nbasis), dBasisdx(ngp,nbasis,3), &
            BasisVec(ngp,nbasis), dBasisdxVec(ngp,nbasis,3), &
            STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('H1Basis',&
              'Storage allocation for local element basis failed')
    END IF

    ! Initialize arrays
    Basis = 0
    dBasisdx = 0
    
    t_tot_n = REAL(0,dp)
    t_tot_b = REAL(0,dp)
    t_start = ftimer()
    DO rep=1,NREP
      ! Nodal basis 
      t_start_tmp=ftimer()
      DO i=1,ngp
        CALL NodalBasisFunctions3D( Basis(i,1:nndof), Element, &
                GP % U(i), GP % V(i), GP % W(i))
        CALL NodalFirstDerivatives3D( dBasisdx(i,1:nndof,1:3), Element, &
                GP % U(i), GP % V(i), GP % W(i))
      END DO
      t_tot_n=t_tot_n+(ftimer()-t_start_tmp)

      ! Bubble basis 
      t_start_tmp=ftimer()
      q = 8
      DO i=2,p-4
        DO j=2,p-i-2
          DO k=2,p-i-j
            q = q + 1
            DO l=1,ngp
              Basis(l,q) = BrickBubblePBasis(i,j,k, &
                      GP % U(l), GP % V(l), GP % W(l))
              dBasisdx(l,q,:) = dBrickBubblePBasis(i,j,k, &
                      GP % U(l), GP % V(l), GP % W(l))
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
    t_totvec_b = REAL(0,dp)
    t_startvec = ftimer()
    DO rep=1,NREP
      t_start_tmp=ftimer()
      CALL H1BrickNodalBasisVec(ngp, GP % U, GP % V, GP % W, BasisVec)
      CALL H1dBrickNodalBasisVec(ngp, GP % U, GP % V, GP % W, dBasisdxVec)
      t_totvec_n=t_totvec_n+(ftimer()-t_start_tmp)
      nbasisvec = 8
      t_start_tmp=ftimer()
      CALL H1BrickBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              nbasisvec, BasisVec)
      ndbasisdxvec = 8
      CALL H1dBrickBubblePBasisVec(ngp, GP % U, GP % V, GP % W, P, &
              ndbasisdxvec, dBasisdxVec)
      t_totvec_b=t_totvec_b+(ftimer()-t_start_tmp)
    END DO
    t_endvec = ftimer()
    
    CALL PrintTestData(Element, ngp, nrep, &
            t_tot_n, t_tot_b, t_end-t_start, &
            t_totvec_n, t_totvec_b, t_endvec-t_startvec)

    nerror = TestBasis(ngp, nbasis, Element % TYPE % DIMENSION, Basis, BasisVec, &
            dBasisdx, dBasisdxVec, tol)

    CALL DeallocatePElement(Element)
    DEALLOCATE(Basis, dBasisdx, BasisVec, dBasisdxVec)
  END FUNCTION TestBrickElement

  FUNCTION TestBasis(ngp, nbasis, ndim, Basis1, Basis2, dBasisdx1, dBasisdx2, tol) RESULT(nerror)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ngp, nbasis, ndim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: Basis1(:,:), Basis2(:,:), &
            dBasisdx1(:,:,:), dBasisdx2(:,:,:)
    REAL(kind=dp), INTENT(IN) :: tol
    INTEGER :: nerror

    INTEGER :: i, j, dim

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
  END FUNCTION TestBasis

  SUBROUTINE PrintTestData(Element, ngp, nrep, t_n1, t_b1, t_tot1, &
          t_n2, t_b2, t_tot2)
  
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER, INTENT(IN) :: ngp, nrep
    REAL(kind=dp), INTENT(IN) :: t_n1, t_b1, t_tot1, t_n2, t_b2, t_tot2

    WRITE (*,'(A,I0)') 'Element type=', Element % TYPE % ElementCode
    IF (ASSOCIATED(Element % PDefs)) THEN
      WRITE (*,'(A,I0)') 'Element polynomial degree=', Element % PDefs % P
    END IF
    WRITE (*,'(A,I0)') 'Element number of nodes=', Element % TYPE % NumberOfNodes
    WRITE (*,'(A,I0)') 'Element number of nonnodal dofs=', GetElementNOFDOFs(Element)
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
    INTEGER :: t, rate
    
#ifdef _OPENMP
    timerval = OMP_GET_WTIME()
#else
    CALL SYSTEM_CLOCK(t,count_rate=rate)
    timerval = REAL(t,KIND(dp))/REAL(rate,KIND(dp))
#endif
  END FUNCTION ftimer

!------------------------------------------------------------------------------
END SUBROUTINE H1Basis
!------------------------------------------------------------------------------
