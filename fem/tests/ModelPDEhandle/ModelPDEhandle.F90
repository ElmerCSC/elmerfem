!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation,
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
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
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter, dim
  LOGICAL :: Found
!------------------------------------------------------------------------------
  TYPE(ValueHandle_t) :: Load_h, FieldSource_h, DiffCoeff_h, ReactCoeff_h, ConvCoeff_h, &
      TimeCoeff_h, ConvVelo1_h, ConvVelo2_h, ConvVelo3_h, &
      BCFlux_h, BCCoeff_h, BCExt_h
      
  CALL ListInitElementKeyword( Load_h,'Body Force','Field Source')
  CALL ListInitElementKeyword( DiffCoeff_h,'Material','Diffusion Coefficient')
  CALL ListInitElementKeyword( ReactCoeff_h,'Material','Reaction Coefficient')
  CALL ListInitElementKeyword( ConvCoeff_h,'Material','Convection Coefficient')
  CALL ListInitElementKeyword( TimeCoeff_h,'Material','Time Derivative Coefficient')
  CALL ListInitElementKeyword( ConvVelo1_h,'Material','Convection Velocity 1')
  CALL ListInitElementKeyword( ConvVelo2_h,'Material','Convection Velocity 2')
  CALL ListInitElementKeyword( ConvVelo3_h,'Material','Convection Velocity 3')

  CALL ListInitElementKeyword( BCFlux_h,'Boundary Condition','Field Flux')
  CALL ListInitElementKeyword( BCCoeff_h,'Boundary Condition','Robin Coefficient')
  CALL ListInitElementKeyword( BCExt_h,'Boundary Condition','External Field')
  
  CALL DefaultStart()
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  dim = CoordinateSystemDimension()


  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()
      CALL LocalMatrix(  Element, n, nd+nb )
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        nb = GetElementNOFBDOFs()
        CALL LocalMatrixBC(  Element, n, nd+nb )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()
  
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: D,C,R, rho,Velo(3,n),a(3), Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------


    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    a = 0.0_dp
    
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
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found ) 
      rho = ListGetElementReal( TimeCoeff_h, Basis, Element, Found ) 

      a(1) = ListGetElementReal( ConvVelo1_h, Basis, Element, Found ) 
      a(2) = ListGetElementReal( ConvVelo2_h, Basis, Element, Found ) 
      IF( dim == 3 ) THEN
        a(3) = ListGetElementReal( ConvVelo3_h, Basis, Element, Found ) 
      END IF
        
      D = ListGetElementReal( DiffCoeff_h, Basis, Element, Found ) 
      C = ListGetElementReal( ConvCoeff_h, Basis, Element, Found ) 
      R = ListGetElementReal( ReactCoeff_h, Basis, Element, Found ) 

      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             D * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd
          ! advection term (C*grad(u),v)
          ! -----------------------------------
          STIFF (p,q) = STIFF(p,q) + Weight * &
             C * SUM(a(1:dim)*dBasisdx(q,1:dim)) * Basis(p)

          ! reaction term (R*u,v)
          ! -----------------------------------
          STIFF(p,q) = STIFF(p,q) + Weight * R*Basis(q) * Basis(p)

          ! time derivative (rho*du/dt,v):
          ! ------------------------------
          MASS(p,q) = MASS(p,q) + Weight * rho * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = ListGetElementReal( BCFlux_h, Basis, Element, Found ) 

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = ListGetElementReal( BCCoeff_h, Basis, Element, Found ) 
      Ext = ListGetElementReal( BCExt_h, Basis, Element, Found ) 

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
        END DO
      END DO

      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)
    END DO
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

! Perform static condensation in case bubble dofs are present
!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------
