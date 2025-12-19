!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver_Init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  LOGICAL :: Found
  
  IF( ListGetLogical( Solver % Values,'CutFEM',Found )  ) THEN
    CALL ListAddNewLogical( Solver % Values,'Use Global Mass Matrix',.TRUE.)
    CALL ListAddNewLogical( Solver % Values,'CutFEM Solver',.TRUE.)
  END IF
     
END SUBROUTINE AdvDiffSolver_Init
  

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
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found
!------------------------------------------------------------------------------
  LOGICAL :: DoCut, IsMore, IsCut, isActive
  TYPE(Matrix_t), POINTER :: SaveMatrix
  TYPE(Element_t), POINTER :: pElement
  REAL(KIND=dp) :: MatVol(0:10)
  
  Mesh => Solver % Mesh

  CALL DefaultStart()

  ! In this version we do on-the-fly splitting of elements. 
  DoCut = ListGetLogical( Solver % Values,'CutFEM',Found ) 
  IsMore = .FALSE.

  ! This is just for checking. "0" refers to the length of the line and others are for different
  ! materials. We can check how well line length of sphere or area of spheare is maintained.   
  MatVol = 0.0_dp
  
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1
  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    
    ! Basic System assembly, skip the cutted and empty elements.
    !----------------------------------------------------------
    !CALL DefaultInitialize()
    Active = GetNOFActive(Solver)
    
    DO t=1,Active
      Element => GetActiveElement(t)
      
      ! We need to add this special section to solver where elements are cut on-the-fly.
      IF(DoCut) THEN
        CALL CutinterfaceCheck(Element,isCut,isActive)
        IF(.NOT. isActive) CYCLE

        ! If element is cut then split it on-the-fly until there are no more pieces. 
10      pElement => CutInterfaceBulk(Element,isCut,isMore)        
        IF(isCut) THEN          
          n  = GetElementNOFNodes(pElement)
          nd = GetElementNOFDOFs(pElement)
          nb = GetElementNOFBDOFs(pElement)
          IF( ALL(Solver % Variable % Perm(pElement % NodeIndexes) > 0 ) ) THEN
            CALL LocalMatrix( pElement, n, nd, CutElem = .TRUE. ) 
          END IF
          IF(IsMore) GOTO 10
          CYCLE
        END IF
      END IF
            
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      CALL LocalMatrix( Element, n, nd+nb )          
    END DO

    CALL DefaultFinishBulkAssembly()

    ! This is a special loop over bulk elements where BC's are created on-the-fly if the element is cut.
    ! We need to create both the bulk integrals and boundary integrals at the same time since we may
    ! be eliminating boundary integrals so they should be collected to the same matrix. 
    IF( DoCut ) THEN
      DO t=1,Active
        Element => GetActiveElement(t)

        ! We may have boundary element even though there is no cut since
        ! if two nodes are splitted then this creates a BC element. 
20      pElement => CutInterfaceBC(Element,isCut,isMore)        
        IF(ASSOCIATED(pElement)) THEN
          n = pElement % TYPE % NumberOfNodes
          nd = n
          IF( ALL(Solver % Variable % Perm(pElement % NodeIndexes) > 0 ) ) THEN
            Model % CurrentElement => pElement
            CALL LocalMatrixBC(pElement, n, nd, CutElem = .TRUE. ) 
          END IF
          IF(IsMore) GOTO 20          
        END IF        
      END DO
    END IF
    
    ! Now the real boundary assembly. Here we assume that this is so far from
    ! the cut interface that we need not worry about a conflict.
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)          
      
      IF(DoCut) THEN
        IF( ANY(Solver % Variable % Perm(Element % NodeIndexes) == 0 ) ) CYCLE
      END IF
      
      IF(ActiveBoundaryElement(Element)) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()                    
        CALL LocalMatrixBC(  Element, n, nd )
      END IF
    END DO
    
    IF( DoCut ) THEN
      IF( InfoActive(10) ) THEN
        PRINT *,'CutFEM Interface area:',MatVol(0)
        DO n=1,10
          IF(MatVol(n) > EPSILON(MatVol(n))) PRINT *,'CutFEM body '//I2S(n)//' volume: ',MatVol(n)
        END DO
      END IF
    END IF
      
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( DefaultConverged() ) EXIT    
  END DO
  
  CALL DefaultFinish()
  Solver % Variable % Norm = Norm
  
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, CutElem )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL, OPTIONAL :: CutElem
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: diff_coeff(n), conv_coeff(n), react_coeff(n), &
                     time_coeff(n), D,C,R, rho,Weight, a(3), Velo(3,n)
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,DoCutElem
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes

    SAVE Nodes
!------------------------------------------------------------------------------
    
    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, Element )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
    
    BodyForce => GetBodyForce(Element)
    DoCutElem = .FALSE.
    IF(PRESENT(CutElem)) DoCutElem = .TRUE.
    
    IF ( ASSOCIATED(BodyForce) ) &
       Load(1:n) = GetReal( BodyForce,'field source', Found, Element )
        
    Material => GetMaterial(Element)
    diff_coeff(1:n)=GetReal(Material,'diffusion coefficient',Found,Element)
    react_coeff(1:n)=GetReal(Material,'reaction coefficient',Found,Element)
    conv_coeff(1:n)=GetReal(Material,'convection coefficient',Found,Element)
    time_coeff(1:n)=GetReal(Material,'time derivative coefficient',Found,Element)
        
    Velo = 0._dp
    DO i=1,dim
      Velo(i,1:n)=GetReal(Material,&
          'convection velocity '//I2S(i),Found,Element)
    END DO
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL Info('AdvDiffSolver','Integration points in 1st element: '//I2S(IP % n),Level=8)
    END IF

    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )
      
      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      rho = SUM(Basis(1:n)*time_coeff(1:n))
      a = MATMUL(Velo(:,1:n),Basis(1:n))
      D = SUM(Basis(1:n)*diff_coeff(1:n))
      C = SUM(Basis(1:n)*conv_coeff(1:n))
      R = SUM(Basis(1:n)*react_coeff(1:n))
      
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
      IF(DoCut) MatVol(Element % BodyId) = MatVol(Element % BodyId) + Weight 
    END DO

    IF(TransientSimulation) THEN
      IF( DoCutElem ) THEN
        BLOCK
          TYPE(Matrix_t), POINTER   :: A
          TYPE(Variable_t), POINTER :: x        
          A => Solver % Matrix
          x => Solver % Variable
          CALL UpdateMassMatrix( A, MASS, nd, x % DOFs, &
              x % Perm(Element % NodeIndexes(1:n)), A % MassValues )
        END BLOCK
      ELSE
        CALL Default1stOrderTime(MASS,STIFF,FORCE,Element)
      END IF
    END IF

    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
    IF( DoCutElem ) THEN
      BLOCK
        TYPE(Matrix_t), POINTER   :: A
        TYPE(Variable_t), POINTER :: x
        A => Solver % Matrix
        x => Solver % Variable                    
        CALL UpdateGlobalEquations( A,STIFF,A % rhs,FORCE,nd,x % DOFs, &
             x % Perm(Element % NodeIndexes), UElement=Element )       
      END BLOCK
    ELSE
      CALL DefaultUpdateEquations(STIFF,FORCE,Element )
    END IF
        
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, CutElem ) 
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL, OPTIONAL :: CutElem
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Flux(n), Coeff(n), Ext_t(n), F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    BC => GetBC(Element)

    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, Element )

    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
    
    Flux(1:n)  = GetReal( BC,'field flux', Found, Element )
    Coeff(1:n) = GetReal( BC,'robin coefficient', Found, Element )
    Ext_t(1:n) = GetReal( BC,'external field', Found, Element )
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = SUM(Basis(1:n)*flux(1:n))

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = SUM(Basis(1:n)*coeff(1:n))
      Ext = SUM(Basis(1:n)*ext_t(1:n))

      DO p=1,nd
        DO q=1,nd
          STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
        END DO
      END DO
      
      FORCE(1:nd) = FORCE(1:nd) + Weight * (F + C*Ext) * Basis(1:nd)      
    END DO
    
    IF( PRESENT(CutElem) ) THEN
      BLOCK
        TYPE(Matrix_t), POINTER   :: A
        TYPE(Variable_t), POINTER :: x
        A => Solver % Matrix
        x => Solver % Variable                    
        CALL UpdateGlobalEquations( A,STIFF,A % rhs,FORCE,nd,x % DOFs, &
             x % Perm(Element % NodeIndexes), UElement=Element )       
      END BLOCK
      MatVol(0) = MatVol(0) + Weight 
    ELSE
      CALL DefaultUpdateEquations(STIFF,FORCE,Element )
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------
