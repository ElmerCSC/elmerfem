SUBROUTINE BDMSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the mixed formulation of the Poisson equation by using the BDM
!  basis of degree k = 1
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
  TYPE(Element_t),POINTER :: Element, Edge, Parent
  TYPE(Nodes_t) :: Nodes, EdgeNodes
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: Var

  REAL(KIND=dp) :: Norm, Normal(3), Normal2(3), &
       u, v, w, NormalSign, Edgeh, PresErr, FluxErr, EK(2), SolNorm, FluxNorm

  INTEGER :: n, nb, np, nd, t, istat, i, j, k, l, active, dim

  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), Bcoef(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A

  LOGICAL :: stat

  INTEGER, ALLOCATABLE :: Indeces(:)

  SAVE STIFF, LOAD, FORCE, Acoef, Bcoef, &
       AllocationsDone, Nodes, EdgeNodes, Indeces
!------------------------------------------------------------------------------
  dim = CoordinateSystemDimension()

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs  ! just big enough
     ALLOCATE( FORCE(N), LOAD(6,N), STIFF(N,N), &
          Acoef(N), Bcoef(N), Indeces(N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'BDMSolver', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  Solver % Matrix % Complex = .FALSE.
  A => GetMatrix()

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize()

  DO t=1,active
     Element => GetActiveElement(t)
     n  = GetElementNOFNodes() ! The number of nodes corresponding to the background mesh
     nd = GetElementNOFDOFs()  ! The total number of degrees of freedom
     nb = size(Element % BubbleIndexes(:)) ! The number of elementwise degrees of freedom

     LOAD = 0.0d0
     BodyForce => GetBodyForce()
     IF ( ASSOCIATED(BodyForce) ) THEN
        Load(1,1:n) = GetReal( BodyForce, 'Body Source 1', Found )
        Load(2,1:n) = GetReal( BodyForce, 'Body Source 2', Found )
        Load(3,1:n) = GetReal( BodyForce, 'Body Source 3', Found )
     END IF

     Acoef(1:n) = 1.0d0
     IF (.FALSE.) THEN
        Material => GetMaterial( Element )
        IF ( ASSOCIATED(Material) ) THEN
           Acoef(1:n) = GetReal( Material, 'Conductivity', Found )
           IF (.NOT. Found) CALL Fatal( 'BDMSolver', 'Conductivity must be specified' )
        END IF
     END IF

     !Get element local matrix and rhs vector:
     !----------------------------------------
     CALL LocalMatrix( STIFF, FORCE, LOAD, Acoef, Element, n, nd, nb, dim)

     !Update global matrix and rhs vector from local matrix & vector:
     !---------------------------------------------------------------
     CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO

  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()  

  
  !-------------------------------------------------------------------
  ! Compute the error for the model problem considered:
  !--------------------------------------------------------------------
   if (.true.) then  

      PresErr = 0.0d0
      FluxErr = 0.0d0
      SolNorm = 0.0d0
      FluxNorm = 0.0d0
      DO t=1,Solver % NumberOfActiveElements

         Element => GetActiveElement(t)
         n  = GetElementNOFNodes()
         nd = GetElementDOFs( Indeces )
         nb = size(Element % BubbleIndexes(:))
         np = 0  ! Change to np=4 if using the element type "n:1 f:3 b:1"

         Load(1,1) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(nd)) )

         Load(2,1) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+1)) )
         Load(2,2) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+2)) )
         Load(2,3) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+3)) )
         Load(2,4) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+4)) )
         Load(2,5) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+5)) )
         Load(2,6) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+6)) )  
         Load(2,7) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+7)) )
         Load(2,8) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+8)) )
         Load(2,9) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+9)) )
         Load(2,10) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+10)) )
         Load(2,11) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+11)) )
         Load(2,12) = Solver % Variable % Values( Solver % Variable % &
              Perm(Indeces(np+12)) )  

         !print *, 'Elementwise pressure sol = ', Load(1,1)

         CALL MyComputeError3d(Load, Element, n, nd, nb, dim, EK, SolNorm, FluxNorm)

         PresErr = PresErr + EK(1)
         FluxErr = FluxErr + EK(2)

      end DO      

      Print *, 'L2 Pressure Error = ', sqrt(PresErr)/sqrt(SolNorm)
      Print *, 'L2 Flux Error = ', sqrt(FluxErr)/sqrt(FluxNorm)    

   end if



CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Acoef, Element, n, nd, nb, dim )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:,:), Acoef(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, nb, dim
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: nu
    REAL(KIND=dp) :: RTBasis(nd-nb,3), DivRTBasis(nd-nb), WBasis(dim)
    REAL(KIND=dp) :: Basis(n), DetJ, F(3,3), C(3,3), L(3), M(3), xq, yq, zq, &
         uq, vq, wq, sq
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0

    !-------------------------------------
    ! Numerical integration over element:
    !-------------------------------------
    IF (.true.) THEN
       !-------------------------------------------------------------------------------------------
       ! The reference element is chosen to be that used for p-approximation,
       ! so we need to switch to using a quadrature which would not be used otherwise
       !--------------------------------------------------------------------------------------------
       IF (dim == 2) THEN
          IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
       ELSE
          !----------------------------------------------------------------
          ! The lowest-order BDM approximation on tetrahedra
          !----------------------------------------------------------------
          IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
       END IF
    ELSE
       !-------------------------------------------------------------
       ! This would work if the classic reference element were used.
       !-------------------------------------------------------------
       IP = GaussPoints(Element)
    END IF

    np = 0  ! Set np = n, if nodal dofs are employed; otherwise set np = 0

    DO t=1,IP % n
       stat = FaceElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), F, detJ, Basis, RTBasis, DivRTBasis, BDM = .TRUE.)

       C = MATMUL( TRANSPOSE(F(1:3,1:3)), F(1:3,1:3) )

       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
       zq = SUM( Nodes % z(1:n) * Basis(1:n) )      

       !nu = SUM( Basis(1:n) * Acoef(1:n) )
 
       !----------------------------------------------------------------
       ! The following branch could be used to produce the 
       ! Galerkin projection of the pressure for visualization.
       !------------------------------------------------------------------
       if (np > 0) then
          DO p = 1,n
             DO q = 1,n       
                STIFF(p,q) = STIFF(p,q) + Basis(p) * Basis(q) * detJ * IP % s(t)    
             END DO

             DO q = nd-nb+1,nd
                STIFF(p,q) = STIFF(p,q) - Basis(p) * 1.0d0 * detJ * IP % s(t)            
             END DO

          end DO
       end if

       !--------------------------------------------------------------
       ! The contribution from the variation with the flux variable q
       !---------------------------------------------------------------
       WBasis = 0.0d0
       DO p = 1,nd-np-nb
          DO i=1,dim
             WBasis(i) = SUM( C(i,1:dim) * RTBasis(p,1:dim) )
          END DO
          i = np + p
          DO q = 1,nd-np-nb
             j = np + q
             STIFF(i,j) = STIFF(i,j) + 1.0d0 * &
                  SUM( RTBasis(q,1:dim) * WBasis(1:dim) ) * 1.0d0/detJ * IP % s(t)
          END DO

          DO q = nd-nb+1,nd
             STIFF(i,q) = STIFF(i,q) + 1.0d0 * DivRTBasis(p) * IP % s(t)
          END DO
       END DO

       !--------------------------------------------------
       ! The contribution from the constraint div q = -f
       !--------------------------------------------------
       DO p = nd-nb+1,nd
          DO q = 1,nd-np-nb
             j = np + q
             STIFF(p,j) = STIFF(p,j) + 1.0d0 * DivRTBasis(q) * IP % s(t)
          END DO
          !-------------------------------------------------------------------
          ! The definition of the in-built body source term is contained here:
          !-------------------------------------------------------------------
          FORCE(p) = FORCE(p) - 128.0d0 * ( (yq**2 - 1.0d0/4.0d0)*(zq**2 - 1.0d0/4.0d0) + &
               (xq**2 - 1.0d0/4.0d0)*(zq**2 - 1.0d0/4.0d0) + &  
               (xq**2 - 1.0d0/4.0d0)*(yq**2 - 1.0d0/4.0d0) ) * &
               1.0d0 * detJ * IP % s(t)
       END DO

    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE MyComputeError3d(  LOAD, Element, n, nd, nb, dim, EK, SolNorm, FluxNorm)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Load(:,:), EK(2), SolNorm, FluxNorm
    TYPE(Element_t), POINTER :: Element    
    INTEGER :: n, nd, nb, dim
!------------------------------------------------------------------
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: i, t
    LOGICAL :: Stat
    REAL(KIND=dp) :: RTBasis(nd-nb,3), DivRTBasis(nd-nb), WBasis(dim), Basis(n), &
         uq, vq, wq, sq
    REAL(KIND=dp) :: xq, yq, zq, DetJ, F(3,3), sol, gradsol(dim,1), fluxsol(dim,1), tmp(1,dim)

    TYPE(Nodes_t), SAVE :: Nodes
!---------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    EK = 0.0d0

    IP = GaussPointsTetra(4)
       
    DO i=1,IP % n  
       uq = IP % u(i) 
       vq = IP % v(i)
       wq = IP % w(i)            
       sq = IP % s(i)
       IP % u(i) = -1.0d0 + 2.0d0*uq + vq + wq
       IP % v(i) = SQRT(3.0d0)*vq + 1.0d0/SQRT(3.0d0)*wq
       IP % w(i) = SQRT(8.0d0)/SQRT(3.0d0)*wq
       IP % s(i) = SQRT(8.0d0)*2.0d0*sq
    END DO


    DO t=1,IP % n
       stat = FaceElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), F, detJ, Basis, RTBasis, DivRTBasis, BDM = .TRUE.)
       
       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
       zq = SUM( Nodes % z(1:n) * Basis(1:n) )      
    
       sol = -64.0d0*(xq**2-1.0d0/4.0d0) * (yq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(1,1) = -64.0d0 * 2.0d0 * xq * (yq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(2,1) = -64.0d0 * 2.0d0 * yq * (xq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(3,1) = -64.0d0 * 2.0d0 * zq * (xq**2-1.0d0/4.0d0) * (yq**2-1.0d0/4.0d0)       

       ! For computing pressure error...
       SolNorm = SolNorm + sol**2 * detJ * IP % s(t)
       EK(1) = EK(1) + (sol - Load(1,1))**2 * detJ * IP % s(t) 

       ! For computing the gradient error...
       FluxNorm = FluxNorm + SUM( gradsol(1:3,1)*gradsol(1:3,1) ) * detJ * IP % s(t)

       tmp(1,1:dim) = RTBasis(1,1:dim)*Load(2,1) + RTBasis(2,1:dim)*Load(2,2) + &
            RTBasis(3,1:dim)*Load(2,3) + RTBasis(4,1:dim)*Load(2,4) + &
            RTBasis(5,1:dim)*Load(2,5) + &
            RTBasis(6,1:dim)*Load(2,6) + RTBasis(7,1:dim)*Load(2,7) + &
            RTBasis(8,1:dim)*Load(2,8) + &
            RTBasis(9,1:dim)*Load(2,9) + RTBasis(10,1:dim)*Load(2,10) + &
            RTBasis(11,1:dim)*Load(2,11) + RTBasis(12,1:dim)*Load(2,12)

       fluxsol(1,1) = 1.0d0/detJ * (F(1,1) * tmp(1,1) + F(1,2) * tmp(1,2) + F(1,3) * tmp(1,3))
       fluxsol(2,1) = 1.0d0/detJ * (F(2,1) * tmp(1,1) + F(2,2) * tmp(1,2) + F(2,3) * tmp(1,3)) 
       fluxsol(3,1) = 1.0d0/detJ * (F(3,1) * tmp(1,1) + F(3,2) * tmp(1,2) + F(3,3) * tmp(1,3))           

       EK(2) = EK(2) + SUM( (gradsol(1:dim,1) - fluxsol(1:dim,1))**2 ) * &
            detJ * IP % s(t) 
          
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE MyComputeError3d
!-----------------------------------------------------------------------------


!------------------------------------------------------------------------------
 END SUBROUTINE BDMSolver
!------------------------------------------------------------------------------
