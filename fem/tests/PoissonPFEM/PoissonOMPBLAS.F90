#ifdef USE_OMP4
#define VECT $OMP SIMD
#define ENDVECT $OMP END SIMD
#else
#define VECT DIR$ SIMD
#define ENDVECT
#endif

SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation with OpenMP threading for assembly and solution
!  phases implemented
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
  USE LocalTypes
  USE ElementBasisFunctions

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
  INTEGER :: nind
  INTEGER, SAVE :: maxdofs
  LOGICAL, SAVE :: AllocationsDone = .FALSE.
  !$OMP THREADPRIVATE(STIFF, LOAD, FORCE, AllocationsDone, maxdofs)

  ! Variables related to graph colouring
  INTEGER :: col, cli, cti, ngd, ngc
  INTEGER, ALLOCATABLE :: dualptr(:), dualind(:), colours(:), &
        cptr(:), cind(:)

  REAL(kind=dp) :: t_start, t_end
  REAL(kind=dp) :: s_start, s_end, ls_start, ls_end, time_s, time_ls

  !------------------------------------------------------------------------------
  Mesh => GetMesh()

  nthreads=1
  !System assembly:
  !----------------
  Active = GetNOFActive()
  CALL DefaultInitialize()

  ! Construct the dual graph from Elmer mesh
  t_start = ftimer()
  CALL ElmerMeshToDualGraph(Mesh, ngd, dualptr, dualind)
  t_end = ftimer()
  WRITE (*,'(A,ES12.3,A)') 'Dual graph creation total: ', t_end - t_start, ' sec.'

  ! Colour the dual graph
  t_start = ftimer()
  CALL ElmerGraphColour(ngd, dualptr, dualind, ngc, colours)
  t_end = ftimer()
  WRITE (*,'(A,ES12.3,A)') 'Graph colouring total: ', t_end - t_start, ' sec.'
  WRITE (*,'(A,I0)') 'Number of colours created ngc=', ngc

  ! Deallocate dual mesh
  DEALLOCATE(dualptr, dualind)

  ! Construct colour lists
  t_start = ftimer()
  CALL ElmerGatherColourLists(ngc, colours, cptr, cind)
  t_end = ftimer()
  WRITE (*,'(A,ES12.3,A)') 'Colour gather total: ', t_end-t_start, ' sec.'

  ! Start timer for solver
  s_start = ftimer()

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(Solver, Mesh, Active,nthreads, ngc, cptr, cind) &
  !$OMP PRIVATE(BodyForce, Element, col, cli, cti, n, nd, nb, t, istat, Found, Norm, &
  !$OMP         Indexes, nind)

#ifdef _OPENMP
  !$OMP SINGLE
  nthreads = omp_get_num_procs()
  WRITE (*,'(A,I0)') 'Number of processors=', nthreads
  nthreads = omp_get_max_threads()
  WRITE (*,'(A,I0)') 'Maximum number of threads=', nthreads
  nthreads = OMP_GET_NUM_THREADS()
  WRITE (*,'(A,I0)') 'Number of threads=', nthreads
  !$OMP END SINGLE
#endif


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

  ! TODO: Perform FE assembly one colour at a time. 
  ! Element indices are listed in CRS structure cptr, cind with ngc colours in total
  DO col=1,ngc
    cli = cptr(col)
    cti = cptr(col+1)-1

    !$OMP DO SCHEDULE(STATIC)
    DO t=cli, cti
      Element => GetActiveElement(cind(t))
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)

      LOAD(1:n) = 0.0d0
      BodyForce => GetBodyForce(Element)
      ! ifort still has a compiler induced race condition bug for the
      ! cases where the contents of a return value pointer have to be copied 
      ! to an array
      !$OMP CRITICAL
      IF ( ASSOCIATED(BodyForce) ) &
            LOAD(1:n) = GetReal( BodyForce, 'Source', Found, UElement=Element )
      !$OMP END CRITICAL
      ! IF (ANY(LOAD(1:n) /= omp_get_thread_num())) THEN
      !      WRITE (*,*) omp_get_thread_num(), ' has a wrong value'
      ! END IF

      !Get element local matrix and rhs vector:
      !----------------------------------------
      ! CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd+nb )
      ! Compare matrices (should be the same)
      ! Norm = DNRM2( (nd+nb)*(nd+nb), STIFF, 1)
      ! WRITE (*,'(A,ES12.5)') '||L||_F=', Norm
      ! Norm = DNRM2( (nd+nb), FORCE, 1)
      ! WRITE (*,'(A,ES12.5)') '||f||_2=', Norm
      ! Vectorized version
      CALL LocalMatrixVec( STIFF, FORCE, LOAD, Element, n, nd+nb )
      ! IF (t>0) STOP
      ! Compare matrices (should be the same)
      ! Norm = DNRM2( (nd+nb)*(nd+nb), STIFF, 1)
      ! WRITE (*,'(A,ES12.5)') '||L||_F=', Norm
      ! Norm = DNRM2( (nd+nb), FORCE, 1)
      ! WRITE (*,'(A,ES12.5)') '||f||_2=', Norm
      ! STOP

      CALL LCondensate( nd, nb, STIFF, FORCE, maxdofs )
      ! CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element )
      ! TEMP, to test vectorized version of glueing process
      Indexes => GetVecIndexStore()
      nind = GetElementDOFs( Indexes, Element, Solver )
      CALL UpdateGlobalEquationsVec( Solver % Matrix, STIFF, Solver % Matrix % RHS, FORCE, nind, &
            Solver % Variable % DOFs, Solver % Variable % Perm(Indexes(1:nind)), &
            UElement=Element, MCAssembly=.FALSE. )
    END DO ! Element loop
    !$OMP END DO

  END DO ! Colour loop
  !$OMP END PARALLEL

  ! Deallocate colouring 
  DEALLOCATE(colours, cptr, cind)

  ! End timer
  s_end = ftimer()
  time_s = (s_end-s_start)

  ! Compute the frobenius norm of the global matrix
  nnz = Solver % Matrix % Rows(Solver % Matrix % NumberOfRows+1)-1
  Norm = DNRM2(nnz, Solver % Matrix % Values, 1)
  WRITE (*,'(A,ES24.16)') '||A||_F=', Norm
  Norm = DNRM2(Solver % Matrix % NumberOfRows, Solver % Matrix % RHS, 1)
  WRITE (*,'(A,ES24.16)') '||b||_2=', Norm

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  ! And finally, solve:
  !--------------------
  WRITE (*,'(A,I0,A,I0)') 'Assembly done, n=', Solver % Matrix % NumberOfRows, ', nelem=', Active
  ! Time the linear solve as well
  ls_start = ftimer()
  Norm = DefaultSolve()
  ls_end = ftimer()

  time_ls = (ls_end-ls_start)

  WRITE (*,'(A, I0)') 'OMP_NUM_THREADS=', nthreads
  WRITE (*,'(A,F12.5)') 'Assembly (s)=', time_s
  WRITE (*,'(A,F12.5)') 'Solve (s)=', time_ls
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !$OMP THREADPRIVATE(Nodes)
    REAL(KIND=dp) :: DASUM

!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    ! IP = GaussPoints( Element, 64 )
    ! IP = GaussPoints( Element, 512 )
    ! print *,Element % Type % ElementCode, IP % n
    ! STOP
    DO t=1,IP % n

      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = LElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                           IP % W(t), detJ, Basis, dBasisdx)

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )
      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
                MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      ! CALL DGEMM('N', 'T', nd, nd, 3, IP % s(t) * DetJ, dBasisdx, nd, dBasisdx, nd, 1D0, STIFF, nd)

      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
      ! CALL DAXPY(nd, IP % s(t) * DetJ * LoadAtIP, Basis, 1, FORCE, 1)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

! Condensate bubble functions from the local stiffness matrix and force vector.
! Uses BLAS and LAPACK
! TODO: get rid of local copies
! TODO: build double complex version via module interface
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

    INTEGER :: m, i, j, l, p, nfo, Ldofs(N), Bdofs(Nb), PIV(Nb)

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

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

	! Can we use blas/lapack without copying the matrices?
    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

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

!------------------------------------------------------------------------------
! Local element routines to make compilation less heavy
   FUNCTION LElementMetric(nDOFs,Elm,Nodes,Metric,DetG,dLBasisdx,LtoGMap) RESULT(Success)
!------------------------------------------------------------------------------
     INTEGER :: nDOFs                !< Number of active nodes in element
     TYPE(Element_t)  :: Elm         !< Element structure
     TYPE(Nodes_t)    :: Nodes       !< element nodal coordinates
     REAL(KIND=dp) :: Metric(3,3)    !< Contravariant metric tensor
     REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
     REAL(KIND=dp) :: DetG           !< SQRT of determinant of element coordinate metric
     REAL(KIND=dp) :: LtoGMap(3,3)   !< Mapping between local and global coordinates
     LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3),s
     REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z

     INTEGER :: cdim,dim,i,j,k,n
!------------------------------------------------------------------------------
     success = .TRUE.

     x => Nodes % x
     y => Nodes % y
     z => Nodes % z
     cdim = CoordinateSystemDimension()
     n = MIN( SIZE(x), nDOFs )
     dim  = elm % TYPE % DIMENSION
!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
     DO i=1,dim
          dx(1,i) = DOT_PRODUCT(dlBasisdx(1:n,i), x(1:n))
          dx(2,i) = DOT_PRODUCT(dlBasisdx(1:n,i), y(1:n))
          dx(3,i) = DOT_PRODUCT(dlBasisdx(1:n,i), z(1:n))
     !    dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
     !    dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
     !    dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
     END DO
     ! The above is actually equivalent to this, but faster
     ! dx(1,1:dim) = MATMUL(x(1:n), dlBasisdx(1:n,1:dim))
     ! dx(2,1:dim) = MATMUL(y(1:n), dlBasisdx(1:n,1:dim))
     ! dx(3,1:dim) = MATMUL(z(1:n), dlBasisdx(1:n,1:dim))
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
     DO i=1,dim
        DO j=1,dim
           s = 0.0d0
           DO k=1,cdim
             s = s + dx(k,i)*dx(k,j)
           END DO
           G(i,j) = s
        END DO
     END DO
!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
       CASE (1)
         DetG  = G(1,1)

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         Metric(1,1) = 1.0d0 / DetG
         DetG  = SQRT( DetG )

!------------------------------------------------------------------------------
!      Surface elements
!------------------------------------------------------------------------------
       CASE (2)
         DetG = ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         Metric(1,1) =  G(2,2) / DetG
         Metric(1,2) = -G(1,2) / DetG
         Metric(2,1) = -G(2,1) / DetG
         Metric(2,2) =  G(1,1) / DetG
         DetG = SQRT(DetG)

!------------------------------------------------------------------------------
!      Volume elements
!------------------------------------------------------------------------------
       CASE (3)
         DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
                G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
                G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )

         IF ( DetG <= TINY( DetG ) ) GOTO 100

         CALL InvertMatrix3x3( G,GI,detG )
         Metric = GI
         DetG = SQRT(DetG)
     END SELECT

     DO i=1,cdim
       DO j=1,dim
         s = 0.0d0
         DO k=1,dim
           s = s + dx(i,k) * Metric(k,j)
         END DO
         LtoGMap(i,j) = s
       END DO
     END DO

!!$      print *,'Element index=', Elm % ElementIndex
!!$      print *,'Basis functions at Gauss point:'
!!$      DO i=1,dim
!!$          print '(8(ES12.4))', dLBasisdx(1:n,i)
!!$      END DO
!!$      print *,'Differentials at Gauss point:'
!!$      DO i=1,3
!!$          print '(3(ES12.4))', dx(i,1:dim)
!!$      END DO
!!$      print *,'Metric tensor at Gauss point:'
!!$      DO i=1,3
!!$          print '(3(ES12.4))', G(i,1:3)
!!$      END DO
!!$      print *,'Inverse of Metric tensor at Gauss point:'
!!$      DO i=1,3
!!$          print '(3(ES12.4))', Metric(i,1:3)
!!$      END DO
!!$      print *,'Determinant=', DetG
!!$      print *,'Local to Global mapping (Jacobian inverse):'
!!$      DO i=1,cdim
!!$          print '(3(ES12.4))', LtoGMap(i,1:dim)
!!$      END DO
! Return here also implies success = .TRUE.
     RETURN


100  Success = .FALSE.

     WRITE( Message,'(A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex
     CALL Error( 'ElementMetric', Message )
     WRITE( Message,'(A,G10.3)') 'DetG:',DetG
     CALL Info( 'ElementMetric', Message, Level=3 )
     DO i=1,cdim
       WRITE( Message,'(A,I0,A,3G10.3)') 'Dir: ',i,' Coord:',x(i),y(i),z(i)
       CALL Info( 'ElementMetric', Message, Level=3 )
     END DO
     IF ( cdim < dim ) THEN
       WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
       CALL Info( 'ElementMetric', Message, Level=3 )
     END IF

!------------------------------------------------------------------------------
   END FUNCTION LElementMetric
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   RECURSIVE FUNCTION LElementInfo( Element, Nodes, u, v, w, detJ, &
     Basis, dBasisdx, ddBasisddx, SecondDerivatives, Bubbles, BasisDegree ) RESULT(stat)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element             !< Element structure
     TYPE(Nodes_t)   :: Nodes                       !< Element nodal coordinates.
     REAL(KIND=dp) :: u                             !< 1st Local coordinate at which to calculate the basis function.
     REAL(KIND=dp) :: v                             !< 2nd local coordinate.
     REAL(KIND=dp) :: w                             !< 3rd local coordinate.
     REAL(KIND=dp) :: detJ                          !< Square root of determinant of element coordinate system metric
     REAL(KIND=dp) :: Basis(:)                      !< Basis function values at (u,v,w)
     REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)       !< Global first derivatives of basis functions at (u,v,w)
     REAL(KIND=dp), OPTIONAL :: ddBasisddx(:,:,:)   !< Global second derivatives of basis functions at (u,v,w) if requested
     INTEGER, OPTIONAL :: BasisDegree(:)            !< Degree of each basis function in Basis(:) vector.
	                                                !! May be used with P element basis functions
     LOGICAL, OPTIONAL :: SecondDerivatives         !< Are the second derivatives needed? (still present for historical reasons)
     LOGICAL, OPTIONAL :: Bubbles                   !< Are the bubbles to be evaluated.
     ! LOGICAL, OPTIONAL :: Debug
     LOGICAL :: Stat                                !< If .FALSE. element is degenerate.
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BubbleValue, dBubbledx(3), t, s, LtoGMap(3,3)
     LOGICAL :: invert, degrees
     INTEGER :: i, j, k, l, q, p, f, n, nb, dim, cdim, locali, localj,  &
          tmp(4), direction(4)
     REAL(KIND=dp) :: LinBasis(8), dLinBasisdx(8,3), ElmMetric(3,3)

     REAL(KIND=dp) :: NodalBasis(Element % TYPE % NumberOfNodes), &
             dLBasisdx(MAX(SIZE(Nodes % x),SIZE(Basis)),3)

     INTEGER :: ldlbasis, ldbasis
     TYPE(Element_t) :: Bubble
     TYPE(Element_t), POINTER :: Edge, Face
     ! LOGICAL :: dbg
!------------------------------------------------------------------------------
     stat = .TRUE.
     n    = Element % TYPE % NumberOfNodes
     dim  = Element % TYPE % DIMENSION
     cdim = CoordinateSystemDimension()

     ! dbg = .FALSE.
	 ! IF (PRESENT(Debug)) dbg = Debug

     IF ( Element % TYPE % ElementCode == 101 ) THEN
        detJ = 1.0d0
        Basis(1) = 1.0d0
        IF ( PRESENT(dBasisdx) ) dBasisdx(1,:) = 0.0d0
        RETURN
     END IF

     Basis = 0.0d0
     CALL NodalBasisFunctions(n, Basis, element, u, v, w)

     dLbasisdx = 0.0d0
     CALL NodalFirstDerivatives(n, dLBasisdx, element, u, v, w)
     ! IF (dbg) THEN
     ! 	WRITE (*,'(3(E14.5))') dLBasisdx(1,1:3)
     ! END IF
     q = n

     ! P ELEMENT CODE:
     ! ---------------
     IF ( isPElement(element) ) THEN

      ! Check for need of P basis degrees and set degree of
      ! linear basis if vector asked:
      ! ---------------------------------------------------
      degrees = .FALSE.
      IF ( PRESENT(BasisDegree)) THEN
        degrees = .TRUE.
        BasisDegree = 0
        BasisDegree(1:n) = 1
      END IF

!------------------------------------------------------------------------------
     SELECT CASE( Element % TYPE % ElementCode )
!------------------------------------------------------------------------------

     ! P element code for line element:
     ! --------------------------------
     CASE(202)
        ! Bubbles of line element
        IF (Element % BDOFs > 0) THEN
           ! For boundary element integration check direction
           invert = .FALSE.
           IF ( Element % PDefs % isEdge .AND. &
                Element % NodeIndexes(1)>Element % NodeIndexes(2) ) invert = .TRUE.

           ! For each bubble in line element get value of basis function
           DO i=1, Element % BDOFs
              IF (q >= SIZE(Basis)) CYCLE
              q = q + 1

              Basis(q) = LineBubblePBasis(i+1,u,invert)
              dLBasisdx(q,1) = dLineBubblePBasis(i+1,u,invert)

              ! Polynomial degree of basis function to vector
              IF (degrees) BasisDegree(q) = 1+i
           END DO
        END IF

!------------------------------------------------------------------------------
! P element code for edges and bubbles of triangle
     CASE(303)
        ! Edges of triangle
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge calculate the value of edge basis function
           DO i=1,3
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Get local number of edge start and endpoint nodes
              tmp(1:2) = getTriangleEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Invert edge for parity if needed
              invert = .FALSE.
              IF ( Element % NodeIndexes(locali)>Element % NodeIndexes(localj) ) invert=.TRUE.

              ! For each dof in edge get value of p basis function
              DO k=1,Edge % BDOFs
                 IF (q >= SIZE(Basis)) CYCLE
                 q = q + 1

                 ! Value of basis functions for edge=i and i=k+1 by parity
                 Basis(q) = TriangleEdgePBasis(i, k+1, u, v, invert)
                 ! Value of derivative of basis function
                 dLBasisdx(q,1:2) = dTriangleEdgePBasis(i, k+1, u, v, invert)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Bubbles of p triangle
        IF ( Element % BDOFs > 0 ) THEN
           ! Get element p
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = NINT( ( 3.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )

           ! For boundary element direction needs to be calculated
           IF (Element % PDefs % isEdge) THEN
              direction = 0
              ! Get direction of this face (mask for face = boundary element nodes)
              direction(1:3) = getTriangleFaceDirection(Element, [ 1,2,3 ])
           END IF

           DO i = 0,p-3
              DO j = 0,p-i-3
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Get bubble basis functions and their derivatives
                 ! 3d Boundary element has a direction
                 IF (Element % PDefs % isEdge) THEN
                    Basis(q) = TriangleEBubblePBasis(i,j,u,v,direction)
                    dLBasisdx(q,1:2) = dTriangleEBubblePBasis(i,j,u,v,direction)
                 ELSE
                 ! 2d element bubbles have no direction
                    Basis(q) = TriangleBubblePBasis(i,j,u,v)
                    dLBasisdx(q,1:2) = dTriangleBubblePBasis(i,j,u,v)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 3+i+j
              END DO
           END DO
        END IF
!------------------------------------------------------------------------------
! P element code for quadrilateral edges and bubbles
     CASE(404)
        ! Edges of p quadrilateral
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge begin node calculate values of edge functions
           DO i=1,4
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Choose correct parity by global edge dofs
              tmp(1:2) = getQuadEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Invert parity if needed
              invert = .FALSE.
              IF (Element % NodeIndexes(locali) > Element % NodeIndexes(localj)) invert = .TRUE.

              ! For each DOF in edge calculate value of p basis function
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! For pyramid square face edges use different basis
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    Basis(q) = QuadPyraEdgePBasis(i,k+1,u,v,invert)
                    dLBasisdx(q,1:2) = dQuadPyraEdgePBasis(i,k+1,u,v,invert)
                 ! Normal case, use basis of quadrilateral
                 ELSE
                    ! Get values of basis functions for edge=i and i=k+1 by parity
                    Basis(q) = QuadEdgePBasis(i,k+1,u,v,invert)
                    ! Get value of derivatives of basis functions
                    dLBasisdx(q,1:2) = dQuadEdgePBasis(i,k+1,u,v,invert)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Bubbles of p quadrilateral
        IF ( Element % BDOFs > 0 ) THEN
          ! Get element P
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p = NINT( ( 5.0d0+SQRT(1.0d0+8.0d0*nb) ) / 2.0d0 )

           ! For boundary element direction needs to be calculated
           IF (Element % PDefs % isEdge) THEN
              direction = 0
              direction = getSquareFaceDirection(Element, [ 1,2,3,4 ])
           END IF

           ! For each bubble calculate value of p basis function
           ! and their derivatives for index pairs i,j>=2, i+j=4,...,p
           DO i=2,(p-2)
              DO j=2,(p-i)
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Get values of bubble functions
                 ! 3D boundary elements have a direction
                 IF (Element % PDefs % isEdge) THEN
                    Basis(q) = QuadBubblePBasis(i,j,u,v,direction)
                    dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v,direction)
                 ELSE
                 ! 2d element bubbles have no direction
                    Basis(q) = QuadBubblePBasis(i,j,u,v)
                    dLBasisdx(q,1:2) = dQuadBubblePBasis(i,j,u,v)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = i+j
              END DO
           END DO
        END IF
!------------------------------------------------------------------------------
! P element code for tetrahedron edges, faces and bubbles
     CASE(504)
        ! Edges of p tetrahedron
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge calculate value of edge functions
           DO i=1,6
              Edge => CurrentModel % Solver % Mesh % Edges (Element % EdgeIndexes(i))

              ! Do not solve edge DOFS if there is not any
              IF (Edge % BDOFs <= 0) CYCLE

              ! For each DOF in edge calculate value of edge functions
              ! and their derivatives for edge=i, i=k+1
              DO k=1, Edge % BDOFs
                 IF (q >= SIZE(Basis)) CYCLE
                 q = q + 1

                 Basis(q) = TetraEdgePBasis(i,k+1,u,v,w, Element % PDefs % TetraType)
                 dLBasisdx(q,1:3) = dTetraEdgePBasis(i,k+1,u,v,w, Element % PDefs % TetraType)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of p tetrahedron
        IF ( ASSOCIATED( Element % FaceIndexes )) THEN
           ! For each face calculate value of face functions
           DO F=1,4
              Face => CurrentModel % Solver % Mesh % Faces (Element % FaceIndexes(F))

              ! Do not solve face DOFs if there is not any
              IF (Face % BDOFs <= 0) CYCLE

              ! Get face p
              p = Face % PDefs % P

              ! For each DOF in face calculate value of face functions and
              ! their derivatives for face=F and index pairs
              ! i,j=0,..,p-3, i+j=0,..,p-3
              DO i=0,p-3
                 DO j=0,p-i-3
                    IF (q >= SIZE(Basis)) CYCLE
                    q = q + 1

                    Basis(q) = TetraFacePBasis(F,i,j,u,v,w, Element % PDefs % TetraType)
                    dLBasisdx(q,1:3) = dTetraFacePBasis(F,i,j,u,v,w, Element % PDefs % TetraType)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 3+i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p tetrahedron
        IF ( Element % BDOFs > 0 ) THEN
           p = Element % PDefs % P

           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=NINT(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

           ! For each DOF in bubbles calculate value of bubble functions
           ! and their derivatives for index pairs
           ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF (q >= SIZE(Basis)) CYCLE
                    q = q + 1

                    Basis(q) = TetraBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,1:3) = dTetraBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 4+i+j+k
                 END DO
              END DO
           END DO

        END IF
!------------------------------------------------------------------------------
! P element code for pyramid edges, faces and bubbles
     CASE(605)
        ! Edges of P Pyramid
        IF (ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in wedge, calculate values of edge functions
           DO i=1,8
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE

              ! Get local indexes of current edge
              tmp(1:2) = getPyramidEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Determine edge direction
              invert = .FALSE.

              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.

              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Get values of edge basis functions and their derivatives
                 Basis(q) = PyramidEdgePBasis(i,k+1,u,v,w,invert)
                 dLBasisdx(q,1:3) = dPyramidEdgePBasis(i,k+1,u,v,w,invert)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of P Pyramid
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in pyramid, calculate values of face functions
           DO F=1,5
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Do not solve face dofs, if there is not any
              IF ( Face % BDOFs <= 0) CYCLE

              ! Get face p
              p = Face % PDefs % P

              ! Handle triangle and square faces separately
              SELECT CASE(F)
              CASE (1)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getPyramidFaceMap(F)
                 direction(1:4) = getSquareFaceDirection( Element, tmp(1:4) )

                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = i+j
                    END DO
                 END DO

              CASE (2,3,4,5)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getPyramidFaceMap(F)
                 direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3) )

                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       Basis(q) = PyramidFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dPyramidFacePBasis(F,i,j,u,v,w,direction)

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              END SELECT
           END DO
        END IF

        ! Bubbles of P Pyramid
        IF (Element % BDOFs >= 0) THEN
           ! Get element p
           p = Element % PDefs % p
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=NINT(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+2)

           ! Calculate value of bubble functions from indexes
           ! i,j,k=0,..,p-4 i+j+k=0,..,p-4
           DO i=0,p-4
              DO j=0,p-i-4
                 DO k=0,p-i-j-4
                    IF ( q >= SIZE(Basis)) CYCLE
                    q = q + 1

                    Basis(q) = PyramidBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dPyramidBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 4+i+j+k
                 END DO
              END DO
           END DO
        END IF

!------------------------------------------------------------------------------
! P element code for wedge edges, faces and bubbles
     CASE(706)
        ! Edges of P Wedge
        IF (ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in wedge, calculate values of edge functions
           DO i=1,9
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE

              ! Get local indexes of current edge
              tmp(1:2) = getWedgeEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Determine edge direction
              invert = .FALSE.
              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.

              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! Use basis compatible with pyramid if neccessary
                 ! @todo Correct this!
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    CALL Fatal('ElementInfo','Pyramid compatible wedge edge basis NIY!')
                 END IF

                 ! Get values of edge basis functions and their derivatives
                 Basis(q) = WedgeEdgePBasis(i,k+1,u,v,w,invert)
                 dLBasisdx(q,1:3) = dWedgeEdgePBasis(i,k+1,u,v,w,invert)

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of P Wedge
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in wedge, calculate values of face functions
           DO F=1,5
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Do not solve face dofs, if there is not any
              IF ( Face % BDOFs <= 0) CYCLE

              p = Face % PDefs % P

              ! Handle triangle and square faces separately
              SELECT CASE(F)
              CASE (1,2)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 tmp(1:4) = getWedgeFaceMap(F)
                 direction(1:3) = getTriangleFaceDirection( Element, tmp(1:3) )

                 ! For each face calculate values of functions from index
                 ! pairs i,j=0,..,p-3 i+j=0,..,p-3
                 DO i=0,p-3
                    DO j=0,p-i-3
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                       dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = 3+i+j
                    END DO
                 END DO
              CASE (3,4,5)
                 direction = 0
                 ! Get global direction vector for enforcing parity
                 invert = .FALSE.
                 tmp(1:4) = getWedgeFaceMap(F)
                 direction(1:4) = getSquareFaceDirection( Element, tmp(1:4) )

                 ! First and second node must form a face in upper or lower triangle
                 IF (.NOT. wedgeOrdering(direction)) THEN
                    invert = .TRUE.
                    tmp(1) = direction(2)
                    direction(2) = direction(4)
                    direction(4) = tmp(1)
                 END IF

                 ! For each face calculate values of functions from index
                 ! pairs i,j=2,..,p-2 i+j=4,..,p
                 DO i=2,p-2
                    DO j=2,p-i
                       IF ( q >= SIZE(Basis) ) CYCLE
                       q = q + 1

                       IF (.NOT. invert) THEN
                          Basis(q) = WedgeFacePBasis(F,i,j,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,i,j,u,v,w,direction)
                       ELSE
                          Basis(q) = WedgeFacePBasis(F,j,i,u,v,w,direction)
                          dLBasisdx(q,:) = dWedgeFacePBasis(F,j,i,u,v,w,direction)
                       END IF

                       ! Polynomial degree of basis function to vector
                       IF (degrees) BasisDegree(q) = i+j
                    END DO
                 END DO
              END SELECT

           END DO
        END IF

        ! Bubbles of P Wedge
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from element
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs( Element, p ), Element % BDOFs )
           p=NINT(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+3)

           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j=0,..,p-5 k=2,..,p-3 i+j+k=2,..,p-3
           DO i=0,p-5
              DO j=0,p-5-i
                 DO k=2,p-3-i-j
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1

                    Basis(q) = WedgeBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dWedgeBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = 3+i+j+k
                 END DO
              END DO
           END DO
        END IF

!------------------------------------------------------------------------------
! P element code for brick edges, faces and bubbles
     CASE(808)
        ! Edges of P brick
        IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
           ! For each edge in brick, calculate values of edge functions
           DO i=1,12
              Edge => CurrentModel % Solver % Mesh % Edges( Element % EdgeIndexes(i) )

              ! Do not solve edge dofs, if there is not any
              IF (Edge % BDOFs <= 0) CYCLE

              ! Get local indexes of current edge
              tmp(1:2) = getBrickEdgeMap(i)
              locali = tmp(1)
              localj = tmp(2)

              ! Determine edge direction
              invert = .FALSE.

              ! Invert edge if local first node has greater global index than second one
              IF ( Element % NodeIndexes(locali) > Element % NodeIndexes(localj) ) invert = .TRUE.

              ! For each DOF in edge calculate values of edge functions
              ! and their derivatives for edge=i and i=k+1
              DO k=1,Edge % BDOFs
                 IF ( q >= SIZE(Basis) ) CYCLE
                 q = q + 1

                 ! For edges connected to pyramid square face, use different basis
                 IF (Edge % PDefs % pyramidQuadEdge) THEN
                    ! Get values of edge basis functions and their derivatives
                    Basis(q) = BrickPyraEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,1:3) = dBrickPyraEdgePBasis(i,k+1,u,v,w,invert)
                 ! Normal case. Use standard brick edge functions
                 ELSE
                    ! Get values of edge basis functions and their derivatives
                    Basis(q) = BrickEdgePBasis(i,k+1,u,v,w,invert)
                    dLBasisdx(q,1:3) = dBrickEdgePBasis(i,k+1,u,v,w,invert)
                 END IF

                 ! Polynomial degree of basis function to vector
                 IF (degrees) BasisDegree(q) = 1+k
              END DO
           END DO
        END IF

        ! Faces of P brick
        IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
           ! For each face in brick, calculate values of face functions
           DO F=1,6
              Face => CurrentModel % Solver % Mesh % Faces( Element % FaceIndexes(F) )

              ! Do not calculate face values if no dofs
              IF (Face % BDOFs <= 0) CYCLE

              ! Get p for face
              p = Face % PDefs % P

              ! Generate direction vector for this face
              tmp(1:4) = getBrickFaceMap(F)
              direction(1:4) = getSquareFaceDirection(Element, tmp)

              ! For each face calculate values of functions from index
              ! pairs i,j=2,..,p-2 i+j=4,..,p
              DO i=2,p-2
                 DO j=2,p-i
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1
                    Basis(q) = BrickFacePBasis(F,i,j,u,v,w,direction)
                    dLBasisdx(q,:) = dBrickFacePBasis(F,i,j,u,v,w,direction)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j
                 END DO
              END DO
           END DO
        END IF

        ! Bubbles of p brick
        IF ( Element % BDOFs > 0 ) THEN
           ! Get p from bubble DOFs
           p = Element % PDefs % P
           nb = MAX( GetBubbleDOFs(Element, p), Element % BDOFs )
           p=NINT(1/3d0*(81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+1d0/ &
                   (81*nb+3*SQRT(-3d0+729*nb**2))**(1/3d0)+4)

           ! For each bubble calculate value of basis function and its derivative
           ! for index pairs i,j,k=2,..,p-4, i+j+k=6,..,p
           DO i=2,p-4
              DO j=2,p-i-2
                 DO k=2,p-i-j
                    IF ( q >= SIZE(Basis) ) CYCLE
                    q = q + 1
                    Basis(q) = BrickBubblePBasis(i,j,k,u,v,w)
                    dLBasisdx(q,:) = dBrickBubblePBasis(i,j,k,u,v,w)

                    ! Polynomial degree of basis function to vector
                    IF (degrees) BasisDegree(q) = i+j+k
                 END DO
              END DO
           END DO
        END IF

     END SELECT
     END IF ! P element flag check
!------------------------------------------------------------------------------

     ! Element (contravariant) metric and square root of determinant
     !--------------------------------------------------------------
     IF ( .NOT. LElementMetric( q, Element, Nodes, &
           ElmMetric, detJ, dLBasisdx, LtoGMap ) ) THEN
        stat = .FALSE.
        RETURN
     END IF

     ! Get global first derivatives:
     !------------------------------
     IF ( PRESENT(dBasisdx) ) THEN
     	! IF (dbg) THEN
     	! WRITE (*,*) 'Metric:'
        ! DO i=1,cdim
		! 	WRITE (*,'(3(E14.5))') LtoGMap(i,1:3)
		! END DO
		! END IF

        dBasisdx = 0.0d0
        ! Use BLAS for p elements with p>1
     	IF (isPElement(element) .AND. q > Element % Type % NumberOfNodes) THEN
    		ldlbasis = SIZE(dLBasisdx,1)
        	ldbasis = SIZE(dBasisdx,1)
        	CALL DGEMM('N', 'T', q, cdim, dim, 1D0, dLBasisdx, ldlbasis, LtoGMap, 3, 0D0, dBasisdx, ldbasis)
     	ELSE
   			! TODO: optimize this
        	DO i=1,q
           		DO j=1,cdim
              		DO k=1,dim
                		dBasisdx(i,j) = dBasisdx(i,j) + dLBasisdx(i,k)*LtoGMap(j,k)
              		END DO
          		END DO
         	END DO

!            print *,'dLBasisdx at Gauss point'
!            DO j=1,cdim
!                print '(8ES12.4)',dLBasisdx(1:q,j)
!            END DO
!            print *,'dBasisdx at Gauss point'
!            DO j=1,cdim
!                print '(8ES12.4)',dBasisdx(1:q,j)
!            END DO
         END IF ! P element flag check
     END IF

     ! Get matrix of second derivatives, if needed:
     !---------------------------------------------
     IF ( PRESENT(ddBasisddx) .AND. PRESENT(SecondDerivatives) ) THEN
       IF ( SecondDerivatives ) THEN
         NodalBasis = 0.0d0
         ddBasisddx(1:n,:,:) = 0.0d0
         DO q=1,n
           NodalBasis(q) = 1.0d0
           CALL GlobalSecondDerivatives(Element,Nodes,NodalBasis, &
               ddBasisddx(q,:,:),u,v,w,ElmMetric,dLBasisdx )
           NodalBasis(q) = 0.0d0
         END DO
       END IF
     END IF

!------------------------------------------------------------------------------
!    Generate bubble basis functions, if requested. Bubble basis is as follows:
!    B_i (=(N_(i+n)) = B * N_i, where N_i:s are the nodal basis functions of
!    the element, and B the basic bubble, i.e. the product of nodal basis
!    functions of the corresponding linear element for triangles and tetras,
!    and product of two diagonally opposed nodal basisfunctions of the
!    correspoding (bi-,tri-)linear element for 1d-elements, quads and hexas.
!------------------------------------------------------------------------------
     IF ( PRESENT( Bubbles ) ) THEN
       Bubble % BDOFs = 0
       NULLIFY( Bubble % PDefs )
       NULLIFY( Bubble % EdgeIndexes )
       NULLIFY( Bubble % FaceIndexes )
       NULLIFY( Bubble % BubbleIndexes )

       IF ( Bubbles .AND. SIZE(Basis) >= 2*n ) THEN

         SELECT CASE(Element % TYPE % ElementCode / 100)
           CASE(2)

              IF ( Element % TYPE % ElementCode == 202 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(202)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                          LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1)
                END DO
              END DO

           CASE(3)

              IF ( Element % TYPE % ElementCode == 303 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(303)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(1,j) * LinBasis(2) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(2,j) * LinBasis(1) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                       dLinBasisdx(3,j) * LinBasis(1) * LinBasis(2)
                END DO
              END DO

           CASE(4)

              IF ( Element % TYPE % ElementCode == 404 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(404)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                             LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(3)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(1,j) * LinBasis(3)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                         dLinBasisdx(3,j) * LinBasis(1)
                END DO
              END DO

           CASE(5)

              IF ( Element % TYPE % ElementCode == 504 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(504)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                            LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(2) * LinBasis(3) * LinBasis(4)
              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(1,j) * &
                                    LinBasis(2) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(2,j) * &
                                    LinBasis(1) * LinBasis(3) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(3,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(4)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * dLinBasisdx(4,j) * &
                                    LinBasis(1) * LinBasis(2) * LinBasis(3)
                END DO
              END DO

           CASE(8)

              IF ( Element % TYPE % ElementCode == 808 ) THEN
                LinBasis(1:n) = Basis(1:n)
                dLinBasisdx(1:n,1:cdim) = dBasisdx(1:n,1:cdim)
              ELSE
                Bubble % TYPE => GetElementType(808)

                stat = ElementInfo( Bubble, nodes, u, v, w, detJ, &
                  LinBasis, dLinBasisdx )
              END IF

              BubbleValue = LinBasis(1) * LinBasis(7)

              DO i=1,n
                Basis(n+i) = Basis(i) * BubbleValue
                DO j=1,cdim
                  dBasisdx(n+i,j) = dBasisdx(i,j) * BubbleValue

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(1,j) * LinBasis(7)

                  dBasisdx(n+i,j) = dBasisdx(n+i,j) + Basis(i) * &
                        dLinBasisdx(7,j) * LinBasis(1)
                END DO
              END DO

         CASE DEFAULT

              WRITE( Message, '(a,i4,a)' ) 'Bubbles for element: ', &
               Element % TYPE % ElementCode, ' are not implemented.'
              CALL Error( 'ElementInfo', Message )
              CALL Fatal( 'ElementInfo', 'Please use p-element basis instead.' )

         END SELECT
       END IF
     END IF
!------------------------------------------------------------------------------
   END FUNCTION LElementInfo
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: IPw(:)
    LOGICAL :: Stat
    INTEGER :: i, dim, ldbasis, ngp, allocstat
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:), dBasisdx(:,:,:), wrk(:,:), DetJ(:), LoadAtIPs(:)
    !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, wrk, DetJ, LoadAtIPs)
!------------------------------------------------------------------------------
    CALL GetElementNodesVec( Nodes, UElement=Element )
    STIFF = 0.0d0
    FORCE = 0.0d0

    ! Get integration points
    IP = GaussPoints( Element )
    ! IP = GaussPoints( Element, 64 )
    ! IP = GaussPoints( Element, 512 )
    IPw => IP % s
    ngp = IP % n

	! Reserve workspace (this should be common for all solvers and handled by the Elmer kernel)
	IF (.NOT. ALLOCATED(Basis)) THEN
		ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), wrk(ngp,nd), &
		         DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
		IF (allocstat /= 0) THEN
		    CALL Fatal('LocalMatrixVec','Storage allocation for local element basis failed')
		END IF
	ELSE IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) THEN
		DEALLOCATE(Basis, dBasisdx, wrk, DetJ, LoadAtIPs)
		ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), wrk(ngp,nd), &
		         DetJ(ngp), LoadAtIPs(ngp), STAT=allocstat)
		IF (allocstat /= 0) THEN
		    CALL Fatal('LocalMatrixVec','Storage allocation for local element basis failed')
		END IF
	END IF

    ldbasis = SIZE(Basis,1)
	! Compute values of all basis functions at all integration points
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, DetJ, nd, Basis, dBasisdx )

	! DO i=1,2
	!     WRITE (*,'(I0,1(ES14.5))'), i, Basis(i,1)
	! 	WRITE (*,'(I0,3(ES14.5))'), i, dBasisdx(i,1,1:3)
	! END DO

	! Compute actual integration weights
	DetJ(1:ngp) = Ipw(1:ngp)*Detj(1:ngp)

	DO dim=1, Element % Type % Dimension
    	! Compute elemental matrix & vector:
      	!----------------------------------------
      	! STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
      	!          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

		! Set wrk as dBasisdx(:,i) * Ip % s(i) * Detj(i)
    	DO i=1,nd
			wrk(1:ngp,i)=DetJ(1:ngp)*dBasisdx(1:ngp,i,dim)
		END DO

   		! Elemental stiffness matrix \grad u \dot \grad u
   		CALL DGEMM('T', 'N', nd, nd, ngp, 1D0, dBasisdx(1,1,dim), ldbasis, wrk, ldbasis, 1D0, STIFF, nd)
   	END DO

    ! Source terms at IPs
    !------------------------------------------
    ! LoadAtIPs(1:ngp) = MATMUL( Basis(1:ngp,1:n), LOAD(1:n) )
    CALL DGEMV('N', ngp, n, 1D0, Basis, ldbasis, LOAD, 1, 0D0, LoadAtIPs, 1)
    ! Integration weights, loads and determinants
    DO i=1,ngp
    	DetJ(i) = DetJ(i)*LoadAtIPs(i)
    	! FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    	! CALL DAXPY(nd, IP % s(ip) * DetJ(ip) * LoadAtIPs(ip), Basis(1,ip), 1, FORCE, 1)
    END DO

    CALL DGEMV('T', ngp, nd, 1D0, Basis, ldbasis, DetJ, 1, 0D0, FORCE, 1)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Local element routines to make compilation less heavy
    FUNCTION ElementMetricVec(ndof, Elm, Nodes, nc, DetJ, dLBasisdx, LtoGMap) RESULT(AllSuccess)
!------------------------------------------------------------------------------
        INTEGER :: ndof                            !< Number of active nodes in element
        TYPE(Element_t)  :: Elm                    !< Element structure
        TYPE(Nodes_t)    :: Nodes                  !< element nodal coordinates
        INTEGER :: nc                              !< Number of points to map
        REAL(KIND=dp) :: DetJ(nc)                  !< SQRT of determinant of element coordinate metric at each point
        REAL(KIND=dp) :: dLBasisdx(:,:,:)          !< Derivatives of element basis function with respect to local coordinates at each point
        REAL(KIND=dp) :: LtoGMap(:,:,:)            !< Mapping between local and global coordinates
        LOGICAL :: AllSuccess                      !< Returns .FALSE. if some point in element is degenerate
!------------------------------------------------------------------------------
!       Local variables
!------------------------------------------------------------------------------
        REAL(KIND=dp) :: Metric(nc,3,3)            ! Workspace
        REAL(KIND=dp) :: dx(nc,3,3),G(nc,3,3)
!DIR$ ATTRIBUTES ALIGN:64::Metric
!DIR$ ATTRIBUTES ALIGN:64::dx
!DIR$ ATTRIBUTES ALIGN:64::G

        REAL(KIND=dp), POINTER :: xyz(:,:)
        REAL(KIND=dp) :: s
        INTEGER :: cdim,dim,i,j,k,l,n,ip
        INTEGER :: ldbasis, ldxyz
!------------------------------------------------------------------------------
        AllSuccess = .TRUE.

        ! Coordinates (single array)
        xyz => Nodes % xyz
        n = MIN( SIZE(Nodes % x, 1), ndof )

        ! Dimensions (coordinate system and element)
        cdim = CoordinateSystemDimension()
        dim  = elm % TYPE % DIMENSION

        ! Leading dimensions for local basis and coordinate arrays
        ldbasis = SIZE(dLBasisdx, 1)
        ldxyz = SIZE(xyz, 1)

        ! For linear, extruded and otherwise regular elements mapping has to be computed
        ! only once, the problem is to identify these cases...
        ! @todo Simplify computation for linear cases
        ! @todo Reduce memory consumption by taking into account that G is symmetric

!------------------------------------------------------------------------------
!       Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
        DO i=1,dim
            CALL DGEMM('N','N',nc, 3, n, &
     	               REAL(1,dp), dLbasisdx(1,1,i), ldbasis, &
     	               xyz, ldxyz, REAL(0, dp), dx(1,1,i), nc)
     	END DO
!------------------------------------------------------------------------------
!       Compute the covariant metric tensor of the element coordinate system (symmetric)
!------------------------------------------------------------------------------
      
      ! G(1:nc,1:dim,1:cdim)=REAL(0,dp)
!DIR$ LOOP COUNT MAX=3
        DO j=1,dim
            ! G is symmetric, compute only the upper triangular part
!DIR$ LOOP COUNT MAX=3
            DO i=1,j
                G(1:nc,i,j) = REAL(0,dp)
                DO k=1,cdim
                    ! @todo This operation can be blocked over l
!VECT
                    DO l=1,nc
                        G(l,i,j)=G(l,i,j)+dx(l,k,i)*dx(l,k,j)
                    END DO
!ENDVECT
                END DO
            END DO
        END DO

!------------------------------------------------------------------------------
!       Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
     	SELECT CASE( dim )
!------------------------------------------------------------------------------
!       Line elements
!------------------------------------------------------------------------------
       	CASE (1)
       	    ! Determinants
            DetJ(1:nc)  = G(1:nc,1,1)

            DO i=1,nc
                IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
                    AllSuccess = .FALSE.
                    EXIT
                END IF
            END DO

            IF (AllSuccess) THEN
!VECT
                DO i=1,nc
                    Metric(i,1,1) = REAL(1,dp)/DetJ(i)
                END DO
!ENDVECT
!VECT
                DO i=1,nc
                    DetJ(i) = SQRT( DetJ(i))
                END DO
!ENDVECT
            END IF


!------------------------------------------------------------------------------
!       Surface elements
!------------------------------------------------------------------------------
        CASE (2)
       	    ! Determinants
!VECT
       	    DO i=1,nc
                ! DetJ(i) = ( G(i,1,1)*G(i,2,2) - G(i,1,2)*G(i,2,1) )
                ! G is symmetric
                DetJ(i) = G(i,1,1)*G(i,2,2)-G(i,1,2)*G(i,1,2)
            END DO
!ENDVECT

            DO i=1,nc
                IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
                    AllSuccess = .FALSE.
                    EXIT
                END IF
            END DO

            IF (AllSuccess) THEN
                ! @todo Since G=G^T, it holds G^{-1}=(G^T)^{-1}
!VECT
                DO i=1,nc
                    ! @todo Move this to a SIMD function InvertMatrix2x2
                    s = REAL(1,dp)/DetJ(i)
                    ! G is symmetric
                    Metric(i,1,1) =  s*G(i,2,2)
                    Metric(i,2,1) = -s*G(i,1,2)

                    Metric(i,1,2) = -s*G(i,1,2)
                    Metric(i,2,2) =  s*G(i,1,1)
                END DO
!ENDVECT
!VECT
                DO i=1,nc
                    DetJ(i) = SQRT(DetJ(i))
                END DO
!ENDVECT
            END IF
!------------------------------------------------------------------------------
!       Volume elements
!------------------------------------------------------------------------------
       	CASE (3)
       	    ! Determinants
!VECT
       	    DO i=1,nc
                ! DetJ(i) = G(i,1,1) * ( G(i,2,2)*G(i,3,3) - G(i,2,3)*G(i,3,2) ) + &
                !           G(i,1,2) * ( G(i,2,3)*G(i,3,1) - G(i,2,1)*G(i,3,3) ) + &
                !           G(i,1,3) * ( G(i,2,1)*G(i,3,2) - G(i,2,2)*G(i,3,1) )
                ! G is symmetric
                DetJ(i) = G(i,1,1)*(G(i,2,2)*G(i,3,3)-G(i,2,3)*G(i,2,3)) + &
                          G(i,1,2)*(G(i,2,3)*G(i,1,3)-G(i,1,2)*G(i,3,3)) + &
                          G(i,1,3)*(G(i,1,2)*G(i,2,3)-G(i,2,2)*G(i,1,3))
            END DO
!ENDVECT

            DO i=1,nc
                IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
                    AllSuccess = .FALSE.
                    EXIT
                END IF
            END DO

            IF (AllSuccess) THEN
                ! @todo Since G=G^T, it holds G^{-1}=(G^T)^{-1}
!VECT
                DO i=1,nc
                    ! @todo Move this to a SIMD function InvertMatrix3x3

                    s = REAL(1,dp) / DetJ(i)
                    ! Metric(i,1,1) =  s * (G(i,2,2)*G(i,3,3) - G(i,3,2)*G(i,2,3))
                    ! Metric(i,2,1) = -s * (G(i,2,1)*G(i,3,3) - G(i,3,1)*G(i,2,3))
                    ! Metric(i,3,1) =  s * (G(i,2,1)*G(i,3,2) - G(i,3,1)*G(i,2,2))
                    ! G is symmetric
                    Metric(i,1,1)= s*(G(i,2,2)*G(i,3,3)-G(i,2,3)*G(i,2,3))
                    Metric(i,2,1)=-s*(G(i,1,2)*G(i,3,3)-G(i,1,3)*G(i,2,3))
                    Metric(i,3,1)= s*(G(i,1,2)*G(i,2,3)-G(i,1,3)*G(i,2,2))

                    ! Metric(i,1,2) = -s * (G(i,1,2)*G(i,3,3) - G(i,3,2)*G(i,1,3))
                    ! Metric(i,2,2) =  s * (G(i,1,1)*G(i,3,3) - G(i,3,1)*G(i,1,3))
                    ! Metric(i,3,2) = -s * (G(i,1,1)*G(i,3,2) - G(i,3,1)*G(i,1,2))
                    ! G is symmetric
                    Metric(i,1,2)=-s*(G(i,1,2)*G(i,3,3)-G(i,2,3)*G(i,1,3))
                    Metric(i,2,2)= s*(G(i,1,1)*G(i,3,3)-G(i,1,3)*G(i,1,3))
                    Metric(i,3,2)=-s*(G(i,1,1)*G(i,2,3)-G(i,1,3)*G(i,1,2))

                    ! Metric(i,1,3) =  s * (G(i,1,2)*G(i,2,3) - G(i,2,2)*G(i,1,3))
                    ! Metric(i,2,3) = -s * (G(i,1,1)*G(i,2,3) - G(i,2,1)*G(i,1,3))
                    ! Metric(i,3,3) =  s * (G(i,1,1)*G(i,2,2) - G(i,2,1)*G(i,1,2))
                    ! G is symmetric
                    Metric(i,1,3)= s*(G(i,1,2)*G(i,2,3)-G(i,2,2)*G(i,1,3))
                    Metric(i,2,3)=-s*(G(i,1,1)*G(i,2,3)-G(i,1,2)*G(i,1,3))
                    Metric(i,3,3)= s*(G(i,1,1)*G(i,2,2)-G(i,1,2)*G(i,1,2))
                END DO
!ENDVECT

!VECT
                DO i=1,nc
                    DetJ(i) = SQRT(DetJ(i))
                END DO
!ENDVECT
            END IF
     	END SELECT

        IF (AllSuccess) THEN
            ! Construct mapping
!DIR$ LOOP COUNT MAX=3
            DO j=1,dim
!DIR$ LOOP COUNT MAX=3
                DO i=1,cdim
                    LtoGMap(1:nc,i,j) = REAL(0,dp)
                    DO k=1,dim
!VECT
                        DO l=1,nc
                            LtoGMap(l,i,j) = LtoGMap(l,i,j) + dx(l,i,k)*Metric(l,k,j)
                        END DO
!ENDVECT
                    END DO
                END DO
            END DO

!!$             DO j=1,nc
!!$                 print *,'Element index=', Elm % ElementIndex
!!$                 print *,'Basis functions at Gauss point:'
!!$                 DO i=1,dim
!!$                     print '(8(ES12.4))', dLBasisdx(j,1:n,i)
!!$                 END DO
!!$                 print *,'Differentials at Gauss point:'
!!$                 DO i=1,3
!!$                     print '(3(ES12.4))', dx(j,i,1:dim)
!!$                 END DO
!!$                 print *,'Metric tensor at Gauss point:'
!!$                 DO i=1,3
!!$                     print '(3(ES12.4))', G(j,i,1:3)
!!$                 END DO
!!$                 print *,'Inverse of metric tensor at Gauss point:'
!!$                 DO i=1,3
!!$                     print '(3(ES12.4))', Metric(j,i,1:3)
!!$                 END DO
!!$                 print *,'Determinant=', DetJ(j)
!!$ 
!!$                 print *,'Local to Global mapping (Jacobian inverse):'
!!$                 DO i=1,cdim
!!$                    print '(3(ES12.4))', LtoGMap(j,i,1:dim)
!!$                END DO
!!$            END DO
        ELSE
!!$            DO j=1,nc
!!$                print *,'Element index=', Elm % ElementIndex
!!$                print *,'Basis functions at Gauss point:'
!!$                DO i=1,dim
!!$                    print '(8(ES12.4))', dLBasisdx(j,1:n,i)
!!$                END DO
!!$                print *,'Differentials at Gauss point:'
!!$                DO i=1,3
!!$                    print '(3(ES12.4))', dx(j,i,1:dim)
!!$                END DO
!!$                print *,'Metric tensor at Gauss point:'
!!$                DO i=1,3
!!$                    print '(3(ES12.4))', G(j,i,1:3)
!!$                END DO
!!$                print *,'Inverse of metric tensor at Gauss point:'
!!$                DO i=1,3
!!$                    print '(3(ES12.4))', Metric(j,i,1:3)
!!$                END DO
!!$                print *,'Determinant=', DetJ(j)
!!$
!!$                print *,'Local to Global mapping (Jacobian inverse):'
!!$                DO i=1,cdim
!!$                    print '(3(ES12.4))', LtoGMap(j,i,1:dim)
!!$                END DO
!!$            END DO

            WRITE( Message,'(A,I0,A,I0,A,I0)') 'Degenerate ',dim,'D element: ',Elm % ElementIndex, ', pt=', i
            CALL Error( 'ElementMetricVec', Message )
            WRITE( Message,'(A,G10.3)') 'DetG:',DetJ(i)
            CALL Info( 'ElementMetricVec', Message, Level=3 )
            DO i=1,cdim
                WRITE( Message,'(A,I0,A,3G10.3)') 'Dir: ',i,' Coord:',xyz(i,1),xyz(i,2),xyz(i,3)
                CALL Info( 'ElementMetricVec', Message, Level=3 )
            END DO
            IF (cdim < dim) THEN
                WRITE( Message,'(A,I0,A,I0)') 'Element dim larger than meshdim: ',dim,' vs. ',cdim
                CALL Info( 'ElementMetricVec', Message, Level=3 )
            END IF
     	END IF
!------------------------------------------------------------------------------
   END FUNCTION ElementMetricVec
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElementInfoVec( Element, Nodes, nc, u, v, w, detJ, nb, Basis, dBasisdx ) RESULT(retval)
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Element_t), TARGET :: Element                     !< Element structure
     TYPE(Nodes_t)   :: Nodes                               !< Element nodal coordinates.
     INTEGER, INTENT(IN) :: nc                              !< Number of local coordinates to compute values of the basis function
     REAL(KIND=dp), INTENT(IN), TARGET :: u(nc)             !< 1st local coordinates at which to calculate the basis function.
     REAL(KIND=dp), INTENT(IN), TARGET :: v(nc)             !< 2nd local coordinates.
     REAL(KIND=dp), INTENT(IN), TARGET :: w(nc)             !< 3rd local coordinates.
     REAL(KIND=dp), INTENT(OUT) :: detJ(nc)                 !< Square roots of determinants of element coordinate system metric at coordinates
     INTEGER, INTENT(IN) :: nb                              !< Maximum number of basis function values to compute at most
     REAL(KIND=dp) :: Basis(nc,nb)                          !< Basis function values at (u,v,w)
     REAL(KIND=dp), OPTIONAL :: dBasisdx(nc,nb,3)           !< Global first derivatives of basis functions at (u,v,w)
     LOGICAL :: retval                                      !< If .FALSE. element is degenerate. or if local storage allocation fails
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: cdim, dim, i, j, k, l, ip, n, nbc, allocstat
	 LOGICAL :: elem

	 ! ELEMENT CACHE
     TYPE(ElementCache_t), SAVE :: ElementCache
	 ! Local storage for basis function values and mapping
     REAL(KIND=dp), ALLOCATABLE, SAVE :: LBasis(:,:), dLBasisdx(:,:,:), &
                    LtoGMaps(:,:,:), dLBasis(:,:), dBasis(:,:)
     !$OMP THREADPRIVATE(LBasis, dLBasisdx, LtoGMaps, dLBasis, dBasis, &
     !$OMP               ElementCache)
!------------------------------------------------------------------------------
	retval = .TRUE.
    n    = Element % TYPE % NumberOfNodes
    dim  = Element % TYPE % DIMENSION
    cdim = CoordinateSystemDimension()

	! Set up local memory
	IF (.NOT. ALLOCATED(LBasis)) THEN
		ALLOCATE(LBasis(nc,nb), dLBasisdx(nc,nb,3), LtoGMaps(nc,3,3), &
		         dLBasis(nb,3), dBasis(nb, 3), STAT=allocstat)
		! Check memory allocation status
		IF (allocstat /= 0) THEN
		    CALL Error('ElementInfoVec','Storage allocation for local element basis failed')
			retval = .FALSE.
			RETURN
		END IF
	ELSE IF (SIZE(LBasis,1) < nc .OR. SIZE(LBasis,2) < nb) THEN
		DEALLOCATE(LBasis, dLBasisdx,LtoGMaps, dLBasis, dBasis)
		ALLOCATE(LBasis(nc,nb), dLBasisdx(nc,nb,3), LtoGMaps(nc,3,3), &
		         dLBasis(nb,3), dBasis(nb,3), STAT=allocstat)
		! Check memory allocation status
		IF (allocstat /= 0) THEN
		    CALL Error('ElementInfoVec','Storage allocation for local element basis failed')
			retval = .FALSE.
			RETURN
		END IF
	END IF

	! P ELEMENT CODE:
    ! ---------------
    IF ( isPElement(element) ) THEN
     	! TODO: implement support
		CALL Error('ElementInfoVec','Vectorized p-elements not implemented yet')
		CALL Fatal( 'ElementInfoVec', 'ElementInfoVec is still experimental.' )
    END IF ! P element flag check

    IF (nb< Element % TYPE % NumberOfNodes) THEN
    	CALL Fatal('ElementInfoVec','Not enough storage to compute local element basis')
    END IF

	! Check if values have been cached (does not work for anything other than regular elements)
	IF (.NOT. isElementBasisInCache( Element, nc, U, V, W, ElementCache ) ) THEN
		! Element local basis not previously computed

		! Compute local nodal basis
 		SELECT CASE (Element % Type % ElementCode)
 		! LINE
 		CASE (202)
 			! Compute nodal basis
 			CALL LineNodalBasisVec(nc, u, LBasis)
 			! Compute local first derivatives
 			CALL dLineNodalBasisVec(nc, u, dLBasisdx)
 			nbc = 2
 		! TRIANGLE
 		CASE (303)
 			! Compute nodal basis
 			CALL TriangleNodalBasisVec(nc, u, v, LBasis)
 			! Compute local first derivatives
 			CALL dTriangleNodalBasisVec(nc, u, v, dLBasisdx)
 			nbc = 3
 		! QUADRILATERAL
 		CASE (404)
 			! Compute nodal basis
 			CALL QuadNodalBasisVec(nc, u, v, LBasis)
 			! Compute local first derivatives
 			CALL dQuadNodalBasisVec(nc, u, v, dLBasisdx)
 			nbc = 4
 		! TETRAHEDRON
 		CASE (504)
 			! Compute nodal basis
 			CALL TetraNodalBasisVec(nc, u, v, w, LBasis)
 			! Compute local first derivatives
 			CALL dTetraNodalBasisVec(nc, u, v, w, dLBasisdx)
 			nbc = 4
 		! WEDGE
 		CASE (706)
 			! Compute nodal basis
 			CALL WedgeNodalBasisVec(nc, u, v, w, LBasis)
 			! Compute local first derivatives
 			CALL dWedgeNodalBasisVec(nc, u, v, w, dLBasisdx)
 			nbc = 6
 		! HEXAHEDRON
    	CASE (808)
    		! Compute local basis
     		CALL BrickNodalBasisVec(nc, u, v, w, LBasis)
     		! Compute local first derivatives
     		CALL dBrickNodalBasisVec(nc, u, v, w, dLBasisdx)
     		nbc = 8
    	CASE DEFAULT
        	WRITE( Message, '(a,i4,a)' ) 'Vectorized basis for element: ', &
              	   Element % TYPE % ElementCode, ' not implemented.'
     		CALL Error( 'ElementInfoVec', Message )
     		CALL Fatal( 'ElementInfoVec', 'ElementInfoVec is still experimental.' )
    	END SELECT

    	! Add element basis to cache
    	CALL CacheElementBasis( Element, nc, U, V, W, ElementCache )
    ELSE
    	! TEMP! This does not work for other than nodal elements
    	nbc = Element % Type % NumberOfNodes
    END IF

!------------------------------------------------------------------------------

    ! Element (contravariant) metric and square root of determinant
    !--------------------------------------------------------------
    elem = ElementMetricVec( nbc, Element, Nodes, nc, detJ, dLBasisdx, LtoGMaps )
    IF (.NOT. elem) THEN
    	retval = .FALSE.
        RETURN
    END IF

    ! Get global basis functions
    !------------------------------
    ! Global basis (copy, use LBasis to implement caching)
    Basis(1:nc,1:nbc)=LBasis(1:nc,1:nbc)

    ! First derivatives
    IF ( PRESENT(dBasisdx) ) THEN
        dBasisdx(1:nc,1:nbc,1:cdim) = REAL(0,dp)
!DIR$ LOOP COUNT MAX=3
        DO j=1,cdim
!DIR$ LOOP COUNT MAX=3
           DO i=1,nbc
!DIR$ LOOP COUNT MAX=3
              DO k=1,dim
! VECT
                    DO l=1,nc
                        ! Map local basis function to global
                        dBasisdx(l,i,j) = dBasisdx(l,i,j) + dLBasisdx(l,i,k)*LtoGMaps(l,j,k)
                    END DO
! ENDVECT
                END DO
            END DO
        END DO

!        DO i=1,nc
!            print *,'dLBasisdx at Gauss point'
!            DO j=1,cdim
!                print '(8ES12.4)',dLBasisdx(i,1:nbc,j)
!            END DO
!            print *,'dBasisdx at Gauss point'
!            DO j=1,cdim
!                print '(8ES12.4)',dBasisdx(i,1:nbc,j)
!            END DO
!        END DO
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ElementInfoVec
!------------------------------------------------------------------------------

   	 FUNCTION isElementBasisInCache( Elem, nc, U, V, W, ElemCache ) RESULT(inCache)
   	 	IMPLICIT NONE
   	 	TYPE(Element_t) :: Elem
   	 	INTEGER, INTENT(IN) :: nc
   	 	REAL(Kind=dp), TARGET :: U(:), V(:), W(:)
   	 	TYPE(ElementCache_t) :: ElemCache
   	 	LOGICAL :: inCache

   	 	inCache = .TRUE.
		IF (ElemCache % nc /= nc .OR. &
		    .NOT. ASSOCIATED(ElemCache % Type, Elem % Type)) THEN
			inCache = .FALSE.
			RETURN
		END IF
		! TODO: this does not currently measure anything, i.e., U, V and W always
		! point to same memory even though contents may vary
		IF (.NOT. ASSOCIATED(ElemCache % U, U) .OR. &
		    .NOT. ASSOCIATED(ElemCache % V, V) .OR. &
		    .NOT. ASSOCIATED(ElemCache % W, W)) THEN

		    inCache = .FALSE.
			RETURN
		END IF
   	 END FUNCTION isElementBasisInCache

   	 SUBROUTINE CacheElementBasis( Elem, nc, U, V, W, ElemCache )
   	 	IMPLICIT NONE
   	 	TYPE(Element_t) :: Elem
   	 	INTEGER, INTENT(IN) :: nc
   	 	REAL(Kind=dp), TARGET :: U(:), V(:), W(:)
   	 	TYPE(ElementCache_t) :: ElemCache

		! Set pointers
		ElemCache % Type => Elem % Type
		ElemCache % nc = nc
		ElemCache % U => U
		ElemCache % V => V
		ElemCache % W => W
   	 END SUBROUTINE

!> Add element local matrices & vectors to global matrices and vectors.
!> Vectorized version, does not support normal or tangential boundary
!> conditions yet.
    SUBROUTINE UpdateGlobalEquationsVec( Gmtr, Lmtr, Gvec, Lvec, n, &
                                         NDOFs, NodeIndexes, RotateNT, UElement, MCAssembly )
        TYPE(Matrix_t), POINTER :: Gmtr         !< The global matrix
        REAL(KIND=dp) :: Lmtr(:,:)              !< Local matrix to be added to the global matrix.
        REAL(KIND=dp) :: Gvec(:)                !< Element local force vector.
        REAL(KIND=dp) :: Lvec(:)                !< The global RHS vector.
        INTEGER :: n                            !< Number of nodes.
        INTEGER :: NDOFs                        !< Number of degrees of free per node.
        INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping.
        LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
        TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
        LOGICAL, OPTIONAL :: MCAssembly   !< Assembly process is multicoloured and guaranteed race condition free 

        INTEGER :: dim, i,j,k
        INTEGER :: Ind(n*NDOFs)
        REAL(KIND=dp) :: Vals(n*NDOFs)
!DIR$ ATTRIBUTES ALIGN:64::Ind, Vals
        TYPE(Element_t), POINTER :: Element
        LOGICAL :: Rotate
        LOGICAL :: ColouredAssembly

        IF (PRESENT(UElement)) THEN
            Element => UElement
        ELSE
            Element => CurrentModel % CurrentElement
        END IF

        IF ( CheckPassiveElement(Element) )  RETURN

        Rotate = .TRUE.
        IF ( PRESENT(RotateNT) ) Rotate = RotateNT

        ColouredAssembly = .FALSE.
        IF ( PRESENT(MCAssembly) ) ColouredAssembly = MCAssembly

        dim = CoordinateSystemDimension()
        ! TEMP
        ! IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN

            ! DO i=1,Element % TYPE % NumberOfNodes
            !     Ind(i) = BoundaryReorder(Element % NodeIndexes(i))
            ! END DO
            ! TODO: See that RotateMatrix is vectorized
            ! CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
            !                    Ind, BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
        IF (Rotate .AND. ndofs >= dim) THEN
            CALL Fatal('UpdateGlobalEquationsVec', &
                       'Normal or tangential boundary conditions not supported yet!')
        END IF

        IF ( ASSOCIATED( Gmtr ) ) THEN
            SELECT CASE( Gmtr % FORMAT )
            CASE( MATRIX_CRS )
                CALL CRS_GlueLocalMatrixVec(Gmtr, n, NDOFs, NodeIndexes, Lmtr, ColouredAssembly)
            CASE DEFAULT
                CALL Fatal('UpdateGlobalEquationsVec','Not implemented for given matrix type')
            END SELECT
        END IF

        ! Check for multicolored assembly
        IF (ColouredAssembly) THEN 
          IF (ANY(NodeIndexes<=0)) THEN
            ! Vector masking needed, no ATOMIC needed
!VECT
            DO i=1,n
              IF (NodeIndexes(i)>0) THEN
!DIR$ LOOP COUNT MIN=1, AVG=3
!DIR$ IVDEP
                DO j=1,NDOFs
                  k = NDOFs*(NodeIndexes(i)-1) + j
                  Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
                END DO
              END IF
            END DO
!ENDVECT
          ELSE
            ! No vector masking needed, no ATOMIC needed
!VECT
            DO i=1,n
!DIR$ LOOP COUNT MIN=1, AVG=3
!DIR$ IVDEP
              DO j=1,NDOFs
                k = NDOFs*(NodeIndexes(i)-1) + j
                Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
              END DO
            END DO
!ENDVECT
          END IF ! Vector masking
        ELSE
          IF (ANY(NodeIndexes<=0)) THEN
            ! Vector masking needed, ATOMIC needed
            DO i=1,n
              IF (NodeIndexes(i)>0) THEN
!DIR$ LOOP COUNT MIN=1, AVG=3
!DIR$ IVDEP
                DO j=1,NDOFs
                  k = NDOFs*(NodeIndexes(i)-1) + j
                  !$OMP ATOMIC
                  Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
                END DO
              END IF
            END DO
          ELSE
            ! No vector masking needed, ATOMIC needed
            DO i=1,n
!DIR$ LOOP COUNT MIN=1, AVG=3
!DIR$ IVDEP
              DO j=1,NDOFs
                k = NDOFs*(NodeIndexes(i)-1) + j
                !$OMP ATOMIC
                Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
              END DO
            END DO
          END IF ! Vector masking
        END IF ! Coloured assembly
    END SUBROUTINE UpdateGlobalEquationsVec


    SUBROUTINE CRS_GlueLocalMatrixVec(Gmtr, N, NDOFs, Indices, Lmtr, MCAssembly)
        TYPE(Matrix_t) :: Gmtr               !< Global matrix
        INTEGER, INTENT(IN) :: N             !< Number of nodes in element
        INTEGER, INTENT(IN) :: NDOFs         !< Number of degrees of freedom for one node
        INTEGER, INTENT(IN) :: Indices(:)    !< Maps element node numbers to global (or partition) node numbers
        REAL(KIND=dp), INTENT(IN) :: Lmtr(:,:)  !< A (N x Dofs) x ( N x Dofs) matrix holding the values to be added
        LOGICAL :: MCAssembly                !< Is the assembly multicolored or not (free of race conditions)

        ! INTEGER, ALLOCATABLE, SAVE :: Lind(:) ! Lind((N*NDOFs)*(N*NDOFs))
        ! REAL(KIND=dp), ALLOCATABLE, SAVE :: Lvals(:) ! Lvals((N*NDOFs)*(N*NDOFs))
        ! !$OMP THREADPRIVATE(Lind, Lvals)
        INTEGER :: Lind((N*NDOFs)*(N*NDOFs))
        REAL(KIND=dp) :: Lvals((N*NDOFs)*(N*NDOFs))

        INTEGER :: i,j, nzind
        INTEGER :: ci, ri, rli, rti, rdof, cdof

        ! INTEGER :: allocstat
        ! LOGICAL :: allocwrk

        INTEGER, POINTER :: gia(:), gja(:)
        REAL(KIND=dp), POINTER :: gval(:)
!DIR$ ATTRIBUTES ALIGN:64::Lind
!DIR$ ATTRIBUTES ALIGN:64::Lvals

!        allocwrk = .FALSE.
!        IF (.NOT. ALLOCATED(Lind)) THEN
!            allocwrk = .TRUE.
!        ELSE IF (SIZE(Lind) < (N*NDOFs)*(N*NDOFs)) THEN
!            DEALLOCATE(Lind, Lvals)
!            allocwrk = .TRUE.
!        END IF
!        IF (allocwrk) THEN
!            ALLOCATE(Lind((N*NDOFs)*(N*NDOFs)), Lvals((N*NDOFs)*(N*NDOFs)), STAT=allocstat)
!            IF (allocstat /= 0) THEN
!                CALL Fatal('CRS_GlueLocalMatrixVec','Workspace allocation failed')
!            END IF
!        END IF

        gia  => Gmtr % Rows
        gja   => Gmtr % Cols
        gval => Gmtr % Values

        ! Check if vector masking is needed
        IF (ANY(Indices<=0)) THEN

            ! Masking and counting needed in the assignment (slower)
            IF (NDOFs == 1) THEN
                ! Separate case for only 1 DOF per node

                ! Construct index array
                nzind = 0
!DIR$ LOOP COUNT MIN=1, AVG=6
                DO i=1,N
                    IF (Indices(i) > 0) THEN
                        ! Row index
                        ri = Indices(i)

                        ! Get row pointers
                        rli = gia(ri)
                        rti = gia(ri+1)-1
!DIR$ LOOP COUNT MIN=1, AVG=6
                        DO j=1,N
                            ! Get global matrix index for entry (ri,Indices(j)).
                            IF (Indices(j) > 0) THEN
                                nzind = nzind + 1
!DIR$ INLINE
                                Lind(nzind)=BinarySearch(gja, Indices(j), rli, rti)
                                Lvals(nzind)=Lmtr(i,j)
                            END IF
                        END DO
                    END IF
                END DO
            ELSE
                ! More than 1 DOF per node
                ! Construct index array
                nzind = 0
!DIR$ LOOP COUNT MIN=1, AVG=6
                DO i=1,N
                    IF (Indices(i) > 0) THEN
!DIR$ LOOP COUNT MIN=2, AVG=3
                        DO rdof=1,NDOFs
                            ! Row index
                            ri = NDOFs*(Indices(i)-1)+rdof

                            ! Get row pointers
                            rli = gia(ri)
                            rti = gia(ri+1)-1
!DIR$ LOOP COUNT MIN=1, AVG=6
                            DO j=1,N
                                IF (Indices(j) > 0) THEN
!DIR$ LOOP COUNT MIN=2, AVG=3
                                    DO cdof=1,NDOFs
                                        ci = NDOFs*(Indices(j)-1)+cdof
                                        ! Get global matrix index for entry (ri,ci).
!DIR$ INLINE
                                        Lind(nzind+cdof)=BinarySearch(gja, ci, rli, rti)
                                        Lvals(nzind+cdof)=Lmtr(NDOFs*(i-1)+rdof,NDOFs*(j-1)+cdof)
                                    END DO
                                    nzind = nzind + cdof
                                END IF
                            END DO
                        END DO
                    END IF
                END DO
            END IF ! NDOFs==1 check

        ELSE
            ! No masking needed (faster)
            IF (NDOFs == 1) THEN
                ! Separate case for only 1 DOF per node

                ! Construct index array
!DIR$ LOOP COUNT MIN=1, AVG=6
                DO i=1,N
                    ! Row index
                    ri = Indices(i)

                    ! Get row pointers
                    rli = gia(ri)
                    rti = gia(ri+1)-1
!DIR$ LOOP COUNT MIN=1, AVG=6
                    DO j=1,N
                        ! Get global matrix index for entry (ri,Indices(j)).
!DIR$ INLINE
                        Lind(N*(i-1)+j)=BinarySearch(gja, Indices(j), rli, rti)
                        Lvals(N*(i-1)+j)=Lmtr(i,j)
                    END DO
                END DO
                nzind = N*N
            ELSE
                ! More than 1 DOF per node

                ! Construct index array
!DIR$ LOOP COUNT MIN=1, AVG=6
                DO i=1,N
!DIR$ LOOP COUNT MIN=2, AVG=3
                    DO rdof=1,NDOFs
                        ! Row index
                        ri = NDOFs*(Indices(i)-1)+rdof

                        ! Get row pointers
                        rli = gia(ri)
                        rti = gia(ri+1)-1
!DIR$ LOOP COUNT MIN=1, AVG=6
                        DO j=1,N
!DIR$ LOOP COUNT MIN=2, AVG=3
                            DO cdof=1,NDOFs
                                ci = NDOFs*(Indices(j)-1)+cdof
                                ! Get global matrix index for entry (ri,ci).
!DIR$ INLINE
                                Lind((NDOFs*N)*(i-1)+NDOFs*(j-1)+cdof)=BinarySearch(gja, ci, rli, rti)
                                Lvals((NDOFs*N)*(i-1)+NDOFs*(j-1)+cdof)=Lmtr(NDOFs*(i-1)+rdof,NDOFs*(j-1)+cdof)
                            END DO
                        END DO
                    END DO

                    nzind = (NDOFs*N)*(NDOFs*N)
                END DO
            END IF ! NDOFs==1 check
        END IF ! Masking check

        ! Sort the indices and the corresponding values
        ! NOTE: this seems to be slow and thus not used, maybe a new sort might help?
!       CALL SortF(nzind, Lind, Lvals)

        ! The actual contribution loop
        IF (MCAssembly) THEN
!DIR$ PREFETCH gval:1:64
!DIR$ PREFETCH gval:0:8
!VECT
          DO i=1,nzind
            gval(Lind(i)) = gval(Lind(i)) + Lvals(i)
          END DO
!ENDVECT
        ELSE
!DIR$ PREFETCH gval:1:64
!DIR$ PREFETCH gval:0:8
          DO i=1,nzind
            !$OMP ATOMIC
            gval(Lind(i)) = gval(Lind(i)) + Lvals(i)
          END DO
        END IF

    END SUBROUTINE CRS_GlueLocalMatrixVec

    PURE FUNCTION BinarySearch(arr, key, lind, tind) RESULT(keyloc)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: arr(:)
        INTEGER, INTENT(IN) :: key, lind, tind

        INTEGER, PARAMETER :: LINSEARCHTHRESH = 8
        INTEGER :: keyloc
        INTEGER :: li, ti, mi

        li = lind
        ti = tind
        DO WHILE ((li+LINSEARCHTHRESH)<ti)
            ! Compute midpoint
            mi = li + ((ti - li) / 2)

            IF (arr(mi)<key) THEN
                li = mi + 1
            ELSE
                ti = mi
            END IF
        END DO

        IF (li<ti) THEN
            keyloc = 0
            DO mi=li,ti
                IF (arr(mi)==key) keyloc = mi
            END DO
        ELSE IF (li == ti .AND. arr(li)==key) THEN
            keyloc = li
        ELSE
            keyloc = 0
        END IF
    END FUNCTION BinarySearch

    FUNCTION GetVecIndexStore() RESULT(ind)
        IMPLICIT NONE
        INTEGER, POINTER :: Ind(:)
        INTEGER :: istat
        INTEGER, ALLOCATABLE, TARGET, SAVE :: VecIndexStore(:)
        !$OMP THREADPRIVATE(VecIndexStore)

        IF ( .NOT. ALLOCATED(VecIndexStore) ) THEN
            ALLOCATE( VecIndexStore(1024), STAT=istat )
            VecIndexStore = 0
            IF ( istat /= 0 ) CALL Fatal( 'GetVecIndexStore', 'Memory allocation error.' )
        END IF

        ind => VecIndexStore( : )
    END FUNCTION GetVecIndexStore

!------------------------------------------------------------------------------
END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------
