!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module for solving heating equation.
! *  Partly vectorized version with handles. 
! *
! *  Authors: Peter Råback & Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 20.01.2020
! * 
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. HeatSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE HeatSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'HeatSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  INTEGER :: dim, n, m
  LOGICAL :: DB, DG
  
  Params => GetSolverParams()
  dim = CoordinateSystemDimension()

  ! Default variable name
  CALL ListAddNewString( Params,'Variable','Temperature')

  ! Tell the matrix structure creation about the need of view factor coupling
  CALL ListAddNewLogical( Params,'Radiation Solver',.TRUE.)
  
  DG = GetLogical( Params,'Discontinuous Galerkin',Found ) 
  DB = GetLogical( Params,'DG Reduced Basis',Found ) 
  
  IF( DG .OR. DB ) THEN
    ! Enforcing indirect nodal connections in parallel for DG just to be sure
    ! The special BCs may require this. 
    CALL ListAddLogical( Params,'DG Indirect Connections',.TRUE.)
  END IF
  
  CALL ListWarnUnsupportedKeyword('body force','Smart Heater Control',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('body force','Integral Heat Source',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('body force','Friction Heat',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('material','Compressibility Model',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('material','Phase Change Model',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('material','Heat Transfer Multiplier',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('equation','Phase Change Model',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('solver','Adaptive Mesh Refinement',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('boundary condition','Phase Change',FatalFound=.TRUE.)

  CALL ListWarnUnsupportedKeyword('boundary condition','Heat Gap',Found )
  IF( Found .AND. .NOT. DG ) THEN
    CALL Fatal(Caller,'Keyword supported only with DG active: "Heat Gap"')
  END IF
   
END SUBROUTINE HeatSolver_Init


!-----------------------------------------------------------------------------
!> A modern version for the heat eqution supporting multi-threading and
!> SIMD friendly ElmerSolver kernels. This tries to be backward compatible
!> with the lagacy HeatSolver but some rarely used features are missing. 
!------------------------------------------------------------------------------
SUBROUTINE HeatSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Radiation
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  INTEGER :: n, nb, nd, t, active, dim
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr
  LOGICAL :: Found, VecAsm, InitHandles, InitDiscontHandles, AxiSymmetric, &
      DG, DB, Newton, HaveFactors, DiffuseGray
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), POINTER :: Temperature(:)
  INTEGER, POINTER :: TempPerm(:)
  REAL(KIND=dp), ALLOCATABLE :: Temps4(:), Emiss(:)
  REAL(KIND=dp) :: Norm, StefBoltz
  CHARACTER(LEN=MAX_NAME_LEN) :: EqName
  CHARACTER(*), PARAMETER :: Caller = 'HeatSolver'

  
  IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving energy equation for temperature')

  ! The View and Gebhardt factors may change if the shape and/or emissivities
  ! have changed. The routine may also affect matrix topology.
  !---------------------------------------------------------------------------
  Mesh => GetMesh()
  CALL RadiationFactors( Solver, .FALSE.) 
  HaveFactors = ListCheckPresentAnyBC( Model,'Radiation')

  IF( HaveFactors ) THEN
    StefBoltz = ListGetConstReal( Model % Constants,&
        'Stefan Boltzmann',UnfoundFatal=HaveFactors)
  END IF

  Temperature => Solver % Variable % Values
  TempPerm => Solver % Variable % Perm
   
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian ) 
  dim = CoordinateSystemDimension() 
  Params => GetSolverParams()  
  EqName = ListGetString( Params,'Equation', Found ) 

  DG = GetLogical( Params,'Discontinuous Galerkin',Found ) 
  DB = GetLogical( Params,'DG Reduced Basis',Found ) 

  maxiter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1
  
  nthr = 1
  !$ nthr = omp_get_max_threads()

  nColours = GetNOFColours(Solver)

  VecAsm = ListGetLogical( Params,'Vector Assembly',Found )
  IF(.NOT. Found ) THEN
    VecAsm = (nColours > 1) .OR. (nthr > 1)
  END IF
  
  IF( VecAsm .AND. AxiSymmetric ) THEN
    CALL Info(Caller,'Vectorized assembly not yet available in axisymmetric case',Level=7)    
    VecAsm = .FALSE.
  END IF
  
  IF( VecAsm ) THEN
    CALL Info(Caller,'Performing vectorized bulk element assembly',Level=7)
  ELSE
    CALL Info(Caller,'Performing non-vectorized bulk element assembly',Level=7)      
  END IF
  
  CALL DefaultStart()

  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    CALL Info(Caller,'Heat solver iteration: '//TRIM(I2S(iter)))

    Newton = GetNewtonActive()
    
    ! Initialize the matrix equation to zero.
    !---------------------------------------
    CALL DefaultInitialize()

    ! For speed compute averaged esissivity and temperature over boundary elements
    ! for diffuse gray radiation.
    !-----------------------------------------------------------------------------
    IF( HaveFactors ) THEN
      IF( iter == 1 ) THEN
        CALL TabulateBoundaryAverages(Mesh, Temps4, Emiss) 
      ELSE
        CALL TabulateBoundaryAverages(Mesh, Temps4 )
      END IF
    END IF
    
    totelem = 0
    
    !$OMP PARALLEL &
    !$OMP SHARED(Solver, Active, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb,col, InitHandles) &
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
   
    DO col=1,nColours
      
      !$OMP SINGLE
      CALL Info( Caller,'Assembly of colour: '//TRIM(I2S(col)),Level=15)
      Active = GetNOFActive(Solver)
      !$OMP END SINGLE
      
      InitHandles = .TRUE.
      !$OMP DO
      DO t=1,Active
        Element => GetActiveElement(t)
        totelem = totelem + 1
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        nb = GetElementNOFBDOFs(Element)
        IF( VecAsm ) THEN
          CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm, InitHandles )
        ELSE
          CALL LocalMatrix(  Element, n, nd+nb, nb, InitHandles )
        END IF
      END DO
      !$OMP END DO
    END DO
    !$OMP END PARALLEL 
    
    totelem = 0
    
    CALL DefaultFinishBulkAssembly()
    
    nColours = GetNOFBoundaryColours(Solver)

    CALL Info(Caller,'Performing boundary element assembly',Level=12)
    
    !$OMP PARALLEL &
    !$OMP SHARED(Active, Solver, nColours, VecAsm, DiffuseGray ) &
    !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    DO col=1,nColours
      !$OMP SINGLE
      CALL Info(Caller,'Assembly of boundary colour: '//TRIM(I2S(col)),Level=10)
      Active = GetNOFBoundaryActive(Solver)
      !$OMP END SINGLE
      
      InitHandles = .TRUE. 
      !$OMP DO
      DO t=1,Active
        Element => GetBoundaryElement(t)
        !WRITE (*,*) Element % ElementIndex
        totelem = totelem + 1
        IF(ActiveBoundaryElement(Element)) THEN
          n  = GetElementNOFNodes(Element)
          nd = GetElementNOFDOFs(Element)
          nb = GetElementNOFBDOFs(Element)

          CALL LocalMatrixBC(  Element, n, nd+nb, nb, VecAsm, DiffuseGray, InitHandles )
          IF( DiffuseGray ) THEN
            CALL LocalMatrixDiffuseGray(  Element, n, nd+nb, nb )
          END IF
        END IF
      END DO
      !$OMP END DO
    END DO
    !$OMP END PARALLEL
    
    IF( DG ) THEN
      BLOCK
        INTEGER :: ElemCount, n1, n2
        TYPE(Element_t), POINTER :: ElemList(:), Parent1, Parent2
        LOGICAL :: BcDone
                
        IF( dim == 2 ) THEN
          ElemCount = Mesh % NumberOfEdges
          ElemList => Mesh % Edges
        ELSE
          ElemCount = Mesh % NumberOfFaces
          ElemList => Mesh % Faces
        END IF

        InitDiscontHandles = .TRUE.      
        DO t = 1, ElemCount
          Element => ElemList(t)
          IF ( .NOT. ActiveBoundaryElement(Element, DGBoundary=.TRUE.) ) CYCLE

          n = GetElementNOFnodes(Element)

          Parent1  => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent1 ) ) CYCLE

          Parent2 => Element  % BoundaryInfo % Right
          IF(.NOT. ASSOCIATED( Parent2 ) ) CYCLE

          n1 = GetElementNOFDOFs(Parent1)
          n2 = GetElementNOFDOFs(Parent2)

          BCDone = .FALSE.
          CALL LocalJumpsDiscontBC( Element, n, Parent1, &
              n1, Parent2, n2, InitDiscontHandles, BCDone )

          IF( .NOT. ( BCDone .OR. DB ) ) THEN          
            CALL LocalJumps( Element, n, Parent1, n1, Parent2, n2 )
          END IF
        END DO
      END BLOCK
    END IF
        
    CALL DefaultFinishBoundaryAssembly()
        
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()
    
    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( DefaultConverged(Solver) ) EXIT

  END DO
  
  CALL DefaultFinish()

CONTAINS 
  

!------------------------------------------------------------------------------
! Assembly of the matrix entries arising from the bulk elements. SIMD version.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm, InitHandles )
!------------------------------------------------------------------------------
    USE LinearForms
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(IN) :: VecAsm
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJVec(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, POINTER  :: CondAtIpVec(:), CpAtIpVec(:), TmpVec(:), &
        SourceAtIpVec(:), RhoAtIpVec(:)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,ngp,allocstat
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
        ConvVelo_h, PerfRate_h, PerfDens_h, PerfCp_h, &
        PerfRefTemp_h
    TYPE(VariableHandle_t), SAVE :: ConvField_h
    
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, DetJVec, &
    !$OMP               MASS, STIFF, FORCE, Nodes, &
    !$OMP               Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
    !$OMP               ConvVelo_h, PerfRate_h, PerfDens_h, PerfCp_h, &
    !$OMP               PerfRefTemp_h, ConvField_h)
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJVec
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Source_h,'Body Force','Heat Source')
      CALL ListInitElementKeyword( Cond_h,'Material','Heat Conductivity')
      CALL ListInitElementKeyword( Cp_h,'Material','Heat Capacity')
      CALL ListInitElementKeyword( Rho_h,'Material','Density')

      CALL ListInitElementKeyword( ConvFlag_h,'Equation','Convection')      
      CALL ListInitElementKeyword( ConvVelo_h,'Material','Convection Velocity',InitVec3D=.TRUE.)

      str = GetString( Params, 'Temperature Convection Field', Found )
      IF(.NOT. Found ) str = 'Flow Solution'
      CALL ListInitElementVariable( ConvField_h, str )
      
      CALL ListInitElementKeyword( PerfRate_h,'Body Force','Perfusion Rate')
      CALL ListInitElementKeyword( PerfDens_h,'Body Force','Perfusion Density')
      CALL ListInitElementKeyword( PerfRefTemp_h,'Body Force','Perfusion Reference Temperature')
      CALL ListInitElementKeyword( PerfCp_h,'Body Force','Perfusion Heat Capacity')
           
      InitHandles = .FALSE.
    END IF
    
    IP = GaussPoints( Element, PReferenceElement = .TRUE.)      
    ngp = IP % n
    
    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
      IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
          DEALLOCATE(Basis, dBasisdx, DetJVec, MASS, STIFF, FORCE, TmpVec )
    END IF
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJVec(ngp), &
          MASS(nd,nd), STIFF(nd,nd), FORCE(nd), &
          TmpVec(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodesVec( Nodes, UElement=Element )
    
    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    
    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJvec, &
        SIZE(Basis,2), Basis, dBasisdx )
    
    ! Compute actual integration weights (recycle the memory space of DetJVec)
    DetJVec(1:ngp) = IP % s(1:ngp) * DetJVec(1:ngp)
    
    RhoAtIpVec => ListGetElementRealVec( Rho_h, ngp, Basis, Element, Found ) 

    ! thermal conductivity term: STIFF=STIFF+(kappa*grad(u),grad(v))
    CondAtIpVec => ListGetElementRealVec( Cond_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, dim, dBasisdx, DetJVec, STIFF, CondAtIpVec )
    END IF

    ! time derivative term: MASS=MASS+(rho*cp*dT/dt,v)
    IF( Transient ) THEN
      CpAtIpVec => ListGetElementRealVec( Cp_h, ngp, Basis, Element, Found ) 
      TmpVec(1:ngp) = CpAtIpVec(1:ngp) * RhoAtIpVec(1:ngp)        
      CALL LinearForms_UdotU(ngp, nd, dim, Basis, DetJVec, MASS, TmpVec )
    END IF
      
    ! source term: FORCE=FORCE+(u,f)
    SourceAtIpVec => ListGetElementRealVec( Source_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      TmpVec(1:ngp) = SourceAtIpVec(1:ngp) * RhoAtIpVec(1:ngp)        
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJVec, TmpVec, FORCE)
    END IF
      
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Assembly of the matrix entries arising from the bulk elements. Not vectorized.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: weight, SourceAtIp, CpAtIp, RhoAtIp, CondAtIp, DetJ, A, VeloAtIp(3)
    REAL(KIND=dp) :: PerfRateAtIp, PerfDensAtIp, PerfCpAtIp, PerfRefTempAtIp, PerfCoeff
    REAL(KIND=dp), POINTER :: CondTensor(:,:)
    LOGICAL :: Stat,Found,ConvComp,ConvConst
    INTEGER :: i,j,t,p,q,m,allocstat,CondRank
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
        ConvVelo_h, PerfRate_h, PerfDens_h, PerfCp_h, PerfRefTemp_h
    TYPE(VariableHandle_t), SAVE :: ConvField_h
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Source_h,'Body Force','Heat Source')
      CALL ListInitElementKeyword( Cond_h,'Material','Heat Conductivity')
      CALL ListInitElementKeyword( Cp_h,'Material','Heat Capacity')
      CALL ListInitElementKeyword( Rho_h,'Material','Density')

      CALL ListInitElementKeyword( ConvFlag_h,'Equation','Convection')
      
      CALL ListInitElementKeyword( ConvVelo_h,'Material','Convection Velocity',InitVec3D=.TRUE.)

      str = GetString( Params, 'Temperature Convection Field', Found )
      IF(.NOT. Found ) str = 'Flow Solution'
      CALL ListInitElementVariable( ConvField_h, str )
      
      CALL ListInitElementKeyword( PerfRate_h,'Body Force','Perfusion Rate')
      CALL ListInitElementKeyword( PerfDens_h,'Body Force','Perfusion Density')
      CALL ListInitElementKeyword( PerfRefTemp_h,'Body Force','Perfusion Reference Temperature')
      CALL ListInitElementKeyword( PerfCp_h,'Body Force','Perfusion Heat Capacity')
     
      InitHandles = .FALSE.
    END IF
    
    IP = GaussPointsAdapt( Element )
      
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(Basis(m), dBasisdx(m,3),&
          MASS(m,m), STIFF(m,m), FORCE(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed in LocalMatrix')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp

    ConvConst = ListCompareElementString( ConvFlag_h,'constant',Element, Found )    
    ConvComp = ListCompareElementString( ConvFlag_h,'computed',Element, Found )
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      IF ( AxiSymmetric ) THEN
        Weight = Weight * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF

      RhoAtIp = ListGetElementReal( Rho_h, Basis, Element, Found )
      
      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      CondAtIp = ListGetElementReal( Cond_h, Basis, Element, Found, &
         GaussPoint = t, Rdim = CondRank, Rtensor = CondTensor ) 
      IF(.NOT. Found ) THEN
        CALL Fatal(Caller,'Required keyword: '//TRIM(Cond_h % Name))
      END IF
      IF( CondRank == 0 ) THEN
        STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
            CondAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )
      ELSE 
        DO p=1,nd
          DO q=1,nd
            A = 0.0_dp
            IF( CondRank == 1 ) THEN
              DO i=1,dim
                A = A + CondTensor(i,1) * dBasisdx(p,i) * dBasisdx(q,i)
              END DO
            ELSE
              DO i=1,dim
                DO j=1,dim
                  A = A + CondTensor(i,j) * dBasisdx(p,i) * dBasisdx(q,j)
                END DO
              END DO
            END IF
            STIFF(p,q) = STIFF(p,q) + Weight * A
          END DO
        END DO
      END IF

      IF( ConvConst .OR. ConvComp .OR. Transient ) THEN
        CpAtIp = ListGetElementReal( Cp_h, Basis, Element, Found )
      END IF
              
      IF( ConvConst .OR. ConvComp ) THEN
        IF( ConvConst ) THEN
          VeloAtIp = ListGetElementReal3D( ConvVelo_h, Basis, Element )
        ELSE
          VeloAtIp = ListGetElementVectorSolution( ConvField_h, Basis, Element, dofs = dim )
        END IF
        
        ! advection term (C*grad(u),v)
        ! -----------------------------------
        DO p=1,nd
          DO q=1,nd
            STIFF (p,q) = STIFF(p,q) + Weight * &
                CpAtIp * RhoAtIp * SUM(VeloAtIp(1:dim)*dBasisdx(q,1:dim)) * Basis(p)
          END DO
        END DO
      END IF
      
      ! reaction term (R*u,v) - perfusion      
      ! -----------------------------------
      PerfRateAtIp = ListGetElementReal( PerfRate_h, Basis, Element, Found )
      IF( Found ) THEN
        PerfDensAtIp = ListGetElementReal( PerfRate_h, Basis, Element, Found )
        PerfCpAtIp = ListGetElementReal( PerfCp_h, Basis, Element, Found )
        PerfRefTempAtIp = ListGetElementReal( PerfRefTemp_h, Basis, Element, Found )
        PerfCoeff = PerfrateAtIp * PerfDensAtIp * PerfCpAtIp 
        DO p=1,nd
          DO q=1,nd        
            STIFF(p,q) = STIFF(p,q) + Weight * PerfCoeff
          END DO
        END DO        
        FORCE(1:nd) = FORCE(1:nd) + Weight * PerfCoeff * PerfRefTempAtIp * Basis(1:nd)
      END IF
                      
      ! Time derivative term
      ! -----------------------------------
      IF( Transient ) THEN
        CpAtIp = ListGetElementReal( Cp_h, Basis, Element, Found )
        RhoAtIp = ListGetElementReal( Rho_h, Basis, Element, Found )
        DO p=1,nd
          MASS(p,1:nd) = MASS(p,1:nd) + Weight * &
                CpAtIp * RhoAtIp * Basis(p) * Basis(1:nd)
        END DO
      END IF

      SourceAtIP = ListGetElementReal( Source_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * RhoAtIp * Basis(1:nd)
      END IF
    END DO
    
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Compute the fraction to be assembled. In serial case it is always one,
! in parallel case only true parents result to assmebly, mixed parents gives
! assembly fraction of 1/2. 
!------------------------------------------------------------------------------
  FUNCTION BCAssemblyFraction( Element ) RESULT ( AssFrac ) 
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: AssFrac

    INTEGER :: NoParents, NoOwners
    
    IF( ParEnv % PEs > 1 ) THEN    
      NoParents = 0; NoOwners = 0
      IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
        IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN   
          NoParents = NoParents + 1
          IF ( Element % BoundaryInfo % Left % PartIndex == ParEnv % myPE ) NoOwners = NoOwners + 1
        END IF
        IF(  ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
          NoParents = NoParents + 1
          IF ( Element % BoundaryInfo % Right % PartIndex == ParEnv % myPE ) NoOwners = NoOwners + 1
        END IF
      END IF
      AssFrac = 1.0_dp * NoOwners / NoParents
    ELSE
      AssFrac = 1.0_dp
    END IF
    
  END FUNCTION BCAssemblyFraction
 !------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------
! Assembly of the matrix entries arising from the Neumann and Robin conditions.
! Also farfield condition and idealized radiation are treated here. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, VecAsm, DiffuseGray, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: VecAsm, DiffuseGray
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,Weight, T0, &
        RadC, RadF, RadText, Text, Emis, AssFrac
    REAL(KIND=dp) :: Basis(nd),DetJ,Coord(3),Normal(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    LOGICAL :: Stat,Found,RobinBC,RadIdeal,RadDiffuse
    INTEGER :: i,j,t,p,q,Indexes(n)
    INTEGER :: NoOwners, NoParents
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: HeatFlux_h, HeatTrans_h, ExtTemp_h, Farfield_h, &
        RadFlag_h, RadExtTemp_h, EmisBC_h, EmisMat_h 
    TYPE(Element_t), POINTER :: Parent

    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes,HeatFlux_h,HeatTrans_h,ExtTemp_h,Farfield_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN
    
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( HeatFlux_h,'Boundary Condition','Heat Flux')
      CALL ListInitElementKeyword( HeatTrans_h,'Boundary Condition','Heat Transfer Coefficient')
      CALL ListInitElementKeyword( ExtTemp_h,'Boundary Condition','External Temperature')
      CALL ListInitElementKeyword( Farfield_h,'Boundary Condition','Farfield Temperature')
      CALL ListInitElementKeyword( RadFlag_h,'Boundary Condition','Radiation')
      CALL ListInitElementKeyword( RadExtTemp_h,'Boundary Condition','Radiation External Temperature')
      CALL ListInitElementKeyword( EmisBC_h,'Boundary Condition','Emissivity')
      CALL ListInitElementKeyword( EmisMat_h,'Material','Emissivity')
            
      InitHandles = .FALSE.
    END IF

    ! In parallel if we have halo the same BC element may occur several times.
    ! Fetch here the fraction of the assembly to be accounted in this occurance. 
    AssFrac = BCAssemblyFraction(Element)
    IF( AssFrac < TINY( AssFrac ) ) RETURN

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
               
    RadIdeal = ListCompareElementString( RadFlag_h,'idealized',Element, Found )    
    RadDiffuse = ListCompareElementString( RadFlag_h,'diffuse gray',Element, Found )

    ! This routine does not do diffuse gray radiation.
    ! Pass on the information to the routine that does. 
    DiffuseGray = RadDiffuse
    Parent => NULL()
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      
      IF ( AxiSymmetric ) THEN
        Weight = Weight * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = ListGetElementReal( HeatFlux_h, Basis, Element, Found )
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
      END IF

      ! Robin condition (c*(T-T_0)):
      ! ---------------------------
      Text = ListGetElementReal( Farfield_h, Basis, Element, RobinBC )      
      IF( RobinBC ) THEN
        Coord(1) = SUM( Nodes % x(1:n)*Basis(1:n) )
        Coord(2) = SUM( Nodes % y(1:n)*Basis(1:n) )
        Coord(3) = SUM( Nodes % z(1:n)*Basis(1:n) )
        Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t), .TRUE. )
        C = SUM( Coord * Normal ) / SUM( Coord * Coord )         
      ELSE
        C = ListGetElementReal( HeatTrans_h, Basis, Element, RobinBC )
        IF(RobinBC) Text = ListGetElementReal( ExtTemp_h, Basis, Element, Found )
      END IF

      IF( RadIdeal ) THEN
        RadText = ListGetElementReal( RadExtTemp_h, Basis, Element, Found )
        IF(.NOT. Found ) THEN
          RadText = ListGetElementReal( ExtTemp_h, Basis, Element, Found )
        END IF

        ! Basis not treated right yet        
        ! Emis = ListGetElementRealParent( EmisMat_h, Basis, Element, Found ) RESULT( RValue ) 
        Emis = ListGetElementRealParent( EmisMat_h, Element = Element, Found = Found )
        IF( .NOT. Found ) THEN
          Emis = ListGetElementReal( EmisBC_h, Basis, Element = Element, Found = Found ) 
        END IF

       IF( DG ) THEN
          T0 = SUM( Basis(1:n) * Temperature(TempPerm(Element % DGIndexes(1:n))))
        ELSE
          T0 = SUM( Basis(1:n) * Temperature(TempPerm(Element % NodeIndexes)))
        END IF
          
        IF( Newton ) THEN
          RadC = StefBoltz * Emis * 4*T0**3
          RadF = StefBoltz * Emis * (3*T0**4+RadText**4) 
        ELSE
          RadC = Emis * StefBoltz * (T0**3 + &
              T0**2*RadText+T0*RadText**2 + RadText**3)
          RadF = RadC * RadText
        END IF
      END IF
        
      IF( RobinBC .OR. RadIdeal) THEN
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * ( C + RadC ) * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * ( C * Text + RadF ) * Basis(1:nd) 
      END IF

    END DO

    IF( ABS(AssFrac-1.0_dp) > TINY( AssFrac ) ) THEN
      FORCE(1:nd) = AssFrac * FORCE(1:nd)
      STIFF(1:nd,1:nd) = AssFrac * STIFF(1:nd,1:nd)
    END IF

    IF( DG ) THEN
      Indexes(1:n) = TempPerm( Element % DGIndexes(1:n) )
      CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
          Solver % Matrix % Rhs, FORCE, n, 1, Indexes(1:n), UElement=Element)      
    ELSE    
      CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! To save some time tabulate the data needed for the diffuse gray radiation. 
! Temps4 is the ^4 averaged temperature over elements.
!------------------------------------------------------------------------------
  SUBROUTINE TabulateBoundaryAverages( Mesh, Temps4, Emiss )
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     REAL(KIND=dp), ALLOCATABLE :: Temps4(:)
     REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: Emiss(:)
 !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
     INTEGER :: bindex, nb, n, j, noactive
     REAL :: NodalEmissivity(12), NodalTemp(12)
     
     nb = Mesh % NumberOfBoundaryElements
     NoActive = 0
     
     DO j=1,nb
       bindex = j + Mesh % NumberOfBulkElements
       Element => Mesh % Elements(bindex)

       BC => GetBC(Element)
       IF(.NOT. ASSOCIATED( BC ) ) CYCLE
       
       IF( ListGetString( BC,'Radiation',Found) /= 'diffuse gray' ) CYCLE
       NoActive = NoActive + 1

       IF(.NOT. ALLOCATED( Temps4 ) ) THEN
         ALLOCATE( Temps4(nb), Emiss(nb) )
         Temps4 = 0.0_dp
         Emiss = 0.0_dp
       END IF
         
       n = GetElementNOFNodes(Element)

       IF( DG ) THEN
         NodalTemp(1:n) = Temperature( TempPerm(Element % DGIndexes) )
       ELSE
         NodalTemp(1:n) = Temperature( TempPerm(Element % NodeIndexes) )
       END IF
       Temps4(j) = ( SUM( NodalTemp(1:n)**4 )/ n )**(1._dp/4._dp)       
       
       IF( PRESENT( Emiss ) ) THEN
         NodalEmissivity(1:n) = GetReal(BC,'Emissivity',Found)
         IF (.NOT. Found) &
             NodalEmissivity(1:n) = GetParentMatProp('Emissivity',Element)
         Emiss(j) = SUM(NodalEmissivity(1:n)) / n
       END IF
     END DO
     
   END SUBROUTINE TabulateBoundaryAverages
!------------------------------------------------------------------------------

   
!------------------------------------------------------------------------------
! Assembly of the matrix entries arising diffuse gray radiation. All the terms
! are treated here. This is a special routine since the view factors create
! additional connections to the matrix. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixDiffuseGray( Element, n, nd, nb )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,T0, Emis, Emis2, RadC, RadF, RadText, Text, Fj, &
        RadLoadAtIp, A1, A2, AngleFraction, Topen, Emis1, AssFrac
    REAL(KIND=dp) :: Basis(nd),DetJ,Coord(3),Normal(3),Atext(12),Base(12),S,RadCoeffAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp), POINTER :: Fact(:)
    TYPE(Element_t), POINTER :: RadElement
    LOGICAL :: Stat,Found,BCOpen
    INTEGER :: i,j,l,t,p,q,Indexes(n),bindex,k,k1,k2,m,nf,nf_imp
    INTEGER :: NoOwners, NoParents
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    INTEGER, POINTER :: ElementList(:),pIndexes(:)
    REAL(KIND=dp), POINTER :: ForceVector(:)
    
    REAL(KIND=dp) :: x,NodalTemp(12),&
        RadCoeff(12),RadLoad(12)
   
    SAVE Nodes
!------------------------------------------------------------------------------
    
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN
    
    IF( ListGetString( BC,'Radiation',Found) /= 'diffuse gray' ) RETURN

    AssFrac = BCAssemblyFraction(Element)
    IF( AssFrac < TINY( AssFrac ) ) RETURN

    n = GetElementNOFNodes(Element)       
    CALL GetElementNodes( Nodes, UElement=Element) 
    n = Element % TYPE % NumberOfNodes
 
    Fact => Element % BoundaryInfo % GebhardtFactors % Factors
    ElementList => Element % BoundaryInfo % GebhardtFactors % Elements

    bindex = Element % ElementIndex - Solver % Mesh % NumberOfBulkElements
    nf = Element % BoundaryInfo % GebhardtFactors % NumberOfFactors
    nf_imp = Element % BoundaryInfo % GebhardtFactors % NumberOfImplicitFactors      
    IF( nf_imp == 0 ) nf_imp = nf

    ForceVector => Solver % Matrix % rhs
    Temperature => Solver % Variable % Values
    TempPerm => Solver % Variable % Perm

    Emis1 = Emiss(bindex)
    
    IP = GaussPoints( Element )       

    BCOpen = GetLogical( BC, 'Radiation Boundary Open', Found)
    IF( BCOpen ) THEN
      AngleFraction = SUM( Fact(1:nf) ) / Emis1
    ELSE
      AngleFraction = 1.0_dp
    END IF

    STIFF(1:n,1:n) = 0.0_dp
    FORCE(1:n) = 0.0_dp      

    IF( DG ) THEN
      NodalTemp(1:n) = Temperature( TempPerm( Element % DGIndexes ) )
    ELSE
      NodalTemp(1:n) = Temperature( TempPerm( Element % NodeIndexes ) )
    END IF
      
    ! Go through surfaces (j) this surface (i) is getting radiated from.
    !------------------------------------------------------------------------------        
    IF ( Newton ) THEN                
      ! Linearization of T^4_i term
      !----------------------------------------------------------------------------
      RadCoeff(1:n) = 4 * Emis1 * NodalTemp(1:n)**3 * StefBoltz
      RadLoad(1:n) = 3 * Emis1 * NodalTemp(1:n)**4 * StefBoltz 
      Base = 0.0_dp

      DO t=1,IP % n
        stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
        s = detJ * IP % s(t)        
        IF ( AxiSymmetric ) THEN
          s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
        END IF

        RadCoeffAtIp = SUM( Basis(1:n) * RadCoeff(1:n) )
        RadLoadAtIp = SUM( Basis(1:n) * RadLoad(1:n) )
                
        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + s * Basis(p) * Basis(q) * RadCoeffAtIp
          END DO
          FORCE(p) = FORCE(p) + s * Basis(p) * RadLoadAtIp
        END DO
                
        Base(1:n) = Base(1:n) + s * Basis(1:n) 
      END DO
      
      ! Linearization of the G_jiT^4_j term
      !------------------------------------------------------------------------------
      DO j=1,nf
        RadElement => Mesh % Elements(ElementList(j))
        k = RadElement % TYPE % NumberOfNodes
        Fj = Fact(j)

        ! Gebhardt factors are given elementwise at the center
        ! of the element, so take average of nodal temperatures
        !-------------------------------------------------------------
        Text = Temps4(j)
        
        IF( j <= nf_imp ) THEN        
          ! Linearization of the G_jiT^4_j term
          !------------------------------------------------------------------------------
          RadLoadAtIp = -3 * Fj * Text**4 * StefBoltz
          RadCoeffAtIp = -4 * Fj * Text**3 * StefBoltz
          
          ! Integrate the contribution of surface j over surface j and add to global matrix
          !------------------------------------------------------------------------------                    
          IF( Dg ) THEN
            DO p=1,n
              k1 = TempPerm( Element % DGIndexes(p) )            
              DO q=1,k
                k2 = TempPerm( RadElement % DGIndexes(q) )
                CALL AddToMatrixElement( Solver % Matrix,k1,k2,RadCoeffAtIp*Base(p)/k )
              END DO
              ForceVector(k1) = ForceVector(k1) + RadLoadAtIp * Base(p)
            END DO
          ELSE
            DO p=1,n
              k1 = TempPerm( Element % NodeIndexes(p) )            
              DO q=1,k
                k2 = TempPerm( RadElement % NodeIndexes(q) )
                CALL AddToMatrixElement( Solver % Matrix,k1,k2,RadCoeffAtIp*Base(p)/k )
              END DO
              ForceVector(k1) = ForceVector(k1) + RadLoadAtIp * Base(p)
            END DO
          END IF
        ELSE
          ! Explicit part, no linearization
          !------------------------------------------------------------------------
          RadLoadAtIp = Fj * Text**4 * StefBoltz
          DO p=1,n
            FORCE(p) = FORCE(p) + RadLoadAtIp * Base(p)
          END DO
        END IF
      END DO
    ELSE
      ! Compute the weighted sum of T^4
      Text = 0.0_dp
      
      DO j=1,nf
        RadElement => Mesh % Elements(ElementList(j))
        k = RadElement % TYPE % NumberOfNodes
        Fj = Fact(j)
        bindex = Element % ElementIndex - Solver % Mesh % NumberOfBulkElements
        Emis2 = Emiss(bindex)

        Text=Text+Emis2*Fj*Temps4(bindex)**4
      END DO

      Text = Text / Emis1**2
    END IF
   

    ! Add the missing part of the incoming radiation in case the boundary is open
    !----------------------------------------------------------------------------
    IF( BCOpen ) THEN
      AText(1:n) = GetReal( BC, 'Radiation External Temperature',Found )
      IF(.NOT. Found) AText(1:n) = GetReal( BC, 'External Temperature' )

      IF( AngleFraction < 1.0_dp ) THEN
        Topen = (SUM( Atext(1:n)**2 ) )**0.25_dp
        IF( Newton ) THEN        
          RadLoadAtIp = (1.0_dp-AngleFraction) * Emis1 * Topen**4 * StefBoltz
          DO p=1,n
            FORCE(p) = FORCE(p) + Base(p) * RadLoadAtIp 
          END DO
        ELSE
          Text = Text + (1.0_dp-AngleFraction) * Topen**4
        END IF
      END IF
    END IF
        
    ! Because we split the product in T^4-T_ext^4 we cannot linearize it before
    ! having computed the complete T_ext^4. So this is done in the end.
    !----------------------------------------------------------------------------
    IF( .NOT. Newton ) THEN      
      Text = Text**0.25
      DO t=1,IP % n
        stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
        s = detJ * IP % s(t)        
        IF ( AxiSymmetric ) THEN
          s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
        END IF

        T0 = SUM( Basis(1:n) * NodalTemp(1:n) )
        RadCoeffAtIp = Emis1 * StefBoltz * &
            (T0**3 + T0**2*Text + T0*Text**2 + Text**3)
                
        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + s * Basis(p) * Basis(q) * RadCoeffAtIp
          END DO
          FORCE(p) = FORCE(p) + s * Basis(p) * RadCoeffAtIp * Text
        END DO
      END DO
    END IF

    
    ! Glue standard local matrix equation to the global matrix
    ! The view factor part has already been glued.
    !-----------------------------------------------------------------
    
    IF( DG ) THEN
      pIndexes => Element % DGIndexes
    ELSE
      pIndexes => Element % NodeIndexes
    END IF

    DO p=1,n
      k1 = TempPerm( pIndexes(p) )
      DO q=1,n
        k2 = TempPerm( pIndexes(q) )
        CALL AddToMatrixElement( Solver % Matrix,k1,k2,STIFF(p,q))
      END DO
      ForceVector(k1) = ForceVector(k1) + FORCE(p)
    END DO
      
  END SUBROUTINE LocalMatrixDiffuseGray
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
! This assembles the local jumps related to standard DG formulation.
! For fully reduced basis this is possibly never needed. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalJumps( Element,n,LeftParent,nl,RightParent,nr)
!------------------------------------------------------------------------------
    INTEGER :: n,nl,nr
    TYPE(Element_t), POINTER :: Element, LeftParent, RightParent
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:),FORCE(:)   
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: LeftBasis(nl), LeftdBasisdx(nl,3)
    REAL(KIND=dp) :: RightBasis(nr), RightdBasisdx(nr,3)
    REAL(KIND=dp) :: LeftdBasisdn(nl), RightdBasisdn(nr)
    REAL(KIND=dp) :: Jump(nl+nr), AverageFlux(nl+nr)
    REAL(KIND=dp) :: detJ, U, V, W, S
    LOGICAL :: Stat
    INTEGER :: i, j, k, p, q, t, m, allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp) :: hE, Normal(3), LeftOut(3), Gamma
    TYPE(Nodes_t) ::Nodes, LeftParentNodes, RightParentNodes
    LOGICAL :: AllocationsDone = .FALSE.

    SAVE Nodes, LeftParentNodes, RightParentNodes, STIFF, FORCE, &
        Basis, dBasisdx, Gamma, AllocationsDone

    !------------------------------------------------------------------------------
    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(Basis(m), dBasisdx(m,3),STIFF(2*m,2*m), FORCE(2*m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed in LocalJumps')
      END IF

      gamma = ListGetCReal( Params,'Dg Continuity Penalty',Found )
      IF(.NOT. Found ) gamma = 0.001

      AllocationsDone = .TRUE.
    END IF

    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    CALL GetElementNodes( Nodes, Element )
    CALL GetElementNodes( LeftParentNodes, LeftParent )
    CALL GetElementNodes( RightParentNodes, RightParent )

    hE = ElementDiameter( Element, Nodes )

    LeftOut(1) = SUM(LeftParentNodes % x(1:nl)) / nl
    LeftOut(2) = SUM(LeftParentNodes % y(1:nl)) / nl
    LeftOut(3) = SUM(LeftParentNodes % z(1:nl)) / nl

    LeftOut(1) = SUM(Nodes % x(1:n)) / n - LeftOut(1)
    LeftOut(2) = SUM(Nodes % y(1:n)) / n - LeftOut(2)
    LeftOut(3) = SUM(Nodes % z(1:n)) / n - LeftOut(3)

    !------------------------------------------------------------------------------
    !      Numerical integration over the edge
    !------------------------------------------------------------------------------
    IP = GaussPoints(Element)
    
    DO k=1,IP % n
      U = IP % u(k)
      V = IP % v(k)
      W = IP % w(k)
      S = IP % s(k)

      !------------------------------------------------------------------------------
      !        Basis function values & derivatives at the integration point
      !------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, U, V, W, detJ, Basis, dBasisdx )

      S = S * detJ
      
      Normal = NormalVector( Element, Nodes, U, V, .FALSE. )
      IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

      ! Find basis functions for the parent elements:
      !----------------------------------------------
      CALL GetParentUVW( Element, n, LeftParent, nl, U, V, W, Basis )

      stat = ElementInfo( LeftParent, LeftParentNodes, &
          U, V, W, detJ, LeftBasis, LeftdBasisdx )
      
      CALL GetParentUVW( Element, n, RightParent, nr, U, V, W, Basis )

      stat = ElementInfo( RightParent, RightParentNodes, &
          U, V, W, detJ, RightBasis, RightdBasisdx )
      
      ! Integrate jump terms:
      !-------------------------
      Jump(1:nl) = LeftBasis(1:nl)
      Jump(nl+1:nl+nr) = -RightBasis(1:nr)
      
      DO i = 1,nl
        LeftdBasisdn(i)  = SUM( LeftdBasisdx(i,:)  * Normal(:) )
      END DO
      
      DO i = 1,nr
        RightdBasisdn(i) = SUM( RightdBasisdx(i,:) * Normal(:) )
      END DO
      
      AverageFlux(1:nl) = LeftdBasisdn(1:nl) / 2.0d0
      AverageFlux(nl+1:nl+nr) = RightdBasisdn(1:nr) / 2.0d0
      
      DO p = 1,nl+nr
        DO q = 1,nl+nr
          STIFF(p,q) = STIFF(p,q) + (gamma/hE)*Jump(p)*Jump(q) * s
          STIFF(p,q) = STIFF(p,q) + AverageFlux(p) * Jump(q)   * s
          STIFF(p,q) = STIFF(p,q) - Jump(p) * AverageFlux(q)   * s
        END DO
      END DO
    END DO
    
    CALL DefaultUpdateEquations( STIFF, FORCE, Element )
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------
! Swap parent elements of a boundary element such that always the Parent1
! belongs to the desired body.
!-------------------------------------------------------------------------
  FUNCTION SwapParentsOnFlag(Parent1, Parent2,FoundJump) RESULT ( Swapped ) 
!-------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent1, Parent2
    TYPE(Element_t), POINTER :: pElem
    LOGICAL :: FoundJump, Swapped    
    TYPE(ValueList_t), POINTER :: Mat
    LOGICAL :: LeftActive, RightActive, HaveJump
    INTEGER :: Body1 = -1, Body2 = -1
        
    SAVE Body1, Body2, LeftActive, HaveJump

    HaveJump = .FALSE.
    Swapped = .FALSE.
    
    ! If we visit subroutine again with same body combination then use also the previous analysis.   
    Mat => GetMaterial( Parent1 )
    LeftActive = ListGetLogical( Mat,'Heat Gap Parent',Found )        
    Mat =>  GetMaterial( Parent2 ) 
    RightActive = ListGetLogical( Mat,'Heat Gap Parent',Found )
      
    IF( LeftActive .AND. RightActive ) THEN
      HaveJump = .FALSE.
    ELSE IF( LeftActive ) THEN
      HaveJump = .TRUE.
    ELSE IF( RightActive ) THEN
      HaveJump = .TRUE.
    ELSE
      HaveJump = .FALSE.
      LeftActive = .TRUE.
    END IF
    
    FoundJump = HaveJump
    IF( .NOT. FoundJump ) RETURN
    
    ! Swicth the reference body always to Parent1
    IF(.NOT. LeftActive ) THEN
      pElem => Parent1
      Parent1 => Parent2
      Parent2 => pElem
      Swapped = .TRUE.
    END IF
    
  END FUNCTION SwapParentsOnFlag
!------------------------------------------------------------------------------

  

!------------------------------------------------------------------------------
! Add jump boundary conditions. These may only occur in conjunction with
! discontinuous Galerkin method. 
!------------------------------------------------------------------------------
  SUBROUTINE LocalJumpsDiscontBC( Element,n,&
      Parent1,n1,Parent2,n2,InitHandles,BCDone)
!------------------------------------------------------------------------------
    INTEGER :: n, n1, n2
    TYPE(Element_t), POINTER :: Element, Parent1, Parent2
    LOGICAL :: InitHandles, BCDone 
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: Basis(n), Jump(n), detJ, S, alpha, beta, AssFrac
    LOGICAL :: Stat
    INTEGER :: i, j, k, p, q, t, i1, i2, JumpOrder, ntmp
    INTEGER :: DgIndexes(2*n)
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t) :: Beta_h, Cond_h, BetaParent_h, CondParent_h
    TYPE(ValueList_t), POINTER :: Mat
    TYPE(Element_t), POINTER :: pElem
    LOGICAL :: DiagJump, Swapped, FoundBodyJump, FoundBCJump
    REAL(KIND=dp) :: Alpha0, Beta0
    LOGICAL :: AllocationsDone = .FALSE.
    INTEGER :: allocstat, m
    
    SAVE Beta_h, Cond_h, BetaParent_h, CondParent_h, Nodes, JumpOrder, DiagJump, &
        Alpha0, Beta0, AllocationsDone, STIFF, FORCE
    
    !------------------------------------------------------------------------------

    BCDone = .FALSE.
    
    ! Both sides need to be active parent elements for a jump condition
    IF( .NOT. CheckElementEquation( Model, Parent1, EqName ) ) RETURN
    IF( .NOT. CheckElementEquation( Model, Parent2, EqName ) ) RETURN

    IF(.NOT. AllocationsDone ) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(STIFF(2*m,2*m), FORCE(2*m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed in LocalJumps')
      END IF
      AllocationsDone = .TRUE.
    END IF
    
    
    IF( InitHandles ) THEN
      JumpOrder = ListGetInteger( Params,'Jump Integration Order',Found )     
      
      DiagJump = ListGetLogical( Params,'Diagonal Jump Glue',Found )
      IF( DiagJump ) THEN
        CALL Info( Caller,'Setting gluing projector to be diagonal',Level=7)
      ELSE
        CALL Info( Caller,'Setting gluing projector to standard Galerkin',Level=7)
      END IF
      
      CALL ListInitElementKeyword( Cond_h,'Boundary Condition','Heat Gap Coefficient')
      CALL ListInitElementKeyword( Beta_h,'Boundary Condition','Heat Gap Flux')      
      CALL ListInitElementKeyword( CondParent_h,'Material','Heat Gap Coefficient')
      CALL ListInitElementKeyword( BetaParent_h,'Material','Heat Gap Flux')      
      
      Alpha0 = ListGetCReal( Params,'Heat Gap Coefficient',Found)
      Beta0 = ListGetCReal( Params,'Heat Gap Flux',Found)
      
      InitHandles = .FALSE.      
    END IF
    
    ! Get the material parameters from the BC 
    FoundBCJump = .FALSE.
    Mat => GetBC( Element ) 
    IF( ASSOCIATED(Mat) ) THEN      
      FoundBCJump =  GetLogical( Mat,'Heat Gap', Found )
    END IF

    ! Or find body jump between two parents
    Swapped = SwapParentsOnFlag( Parent1, Parent2, FoundBodyJump )
    
    IF( .NOT. ( FoundBCJump .OR. FoundBodyJump ) ) RETURN

    AssFrac = BCAssemblyFraction(Element)
    IF( AssFrac < TINY( AssFrac ) ) RETURN    
    
    IF( FoundBCJump ) THEN
      pElem => Element
    ELSE
      pElem => Parent1
    END IF
    
    IF( Swapped ) THEN
      ntmp = n1
      n1 = n2
      n2 = ntmp
    END IF
    
    ! Find the DG indexes for the local assembly 
    !---------------------------------------------
    DgIndexes(1:2*n) = 0
    DO i=1,n
      j = Element % NodeIndexes(i)
      DO i1 = 1, n1
        IF( Parent1 % NodeIndexes( i1 ) == j ) THEN
          DgIndexes(i) = Parent1 % DGIndexes( i1 )
          EXIT
        END IF
      END DO
      DO i2 = 1, n2
        IF( Parent2 % NodeIndexes( i2 ) == j ) THEN
          DgIndexes(n+i) = Parent2 % DGIndexes( i2 )
          EXIT
        END IF
      END DO
    END DO

    IF( ANY( DgIndexes(1:2*n) == 0 ) ) THEN
      CALL Fatal(Caller,'There should not be zero DG indexes!')
    END IF
    
    DgIndexes(1:2*n) = TempPerm( DgIndexes(1:2*n) )
        
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    CALL GetElementNodes( Nodes, Element )
   
    !------------------------------------------------------------------------------
    !      Numerical integration over the edge
    !------------------------------------------------------------------------------
    IP = GaussPoints(Element,RelOrder=JumpOrder)

    DO t=1,IP % n
      !------------------------------------------------------------------------------
      !        Basis function values & derivatives at the integration point
      !------------------------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis )

      S = IP % s(t) * detJ
      IF ( AxiSymmetric ) THEN
        S = S * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF

      IF( FoundBCJump ) THEN
        alpha = ListGetElementReal( Cond_h, Basis, pElem, Found )
      ELSE
        alpha = ListGetElementReal( CondParent_h, Basis, pElem, Found )
      END IF
      IF(.NOT. Found ) alpha = alpha0
      
      IF( FoundBCJump ) THEN
        beta = ListGetElementReal( Beta_h, Basis, pElem, Found ) 
      ELSE
        beta = ListGetElementReal( BetaParent_h, Basis, pElem, Found ) 
      END IF
      IF(.NOT. Found ) beta = beta0
      
      DO p = 1,n                  
        IF( DiagJump ) THEN          
          ! 1st side
          STIFF(p,p) = STIFF(p,p) + alpha * Basis(p) * s 
          STIFF(p,n+p) = STIFF(p,n+p) - alpha * Basis(p) * s 
          ! 2nd side
          STIFF(p+n,p) = STIFF(p+n,p) - alpha * Basis(p) * s 
          STIFF(p+n,n+p) = STIFF(p+n,n+p) + alpha * Basis(p) * s 
        ELSE
          DO q = 1,n
            STIFF(p,q) = STIFF(p,q) + alpha * Basis(p) * Basis(q) * s 
            STIFF(p,n+q) = STIFF(p,n+q) - alpha * Basis(p) * Basis(q) * s 
            STIFF(p+n,q) = STIFF(p+n,q) - alpha * Basis(p) * Basis(q) * s 
            STIFF(p+n,n+q) = STIFF(p+n,n+q) + alpha * Basis(p) * Basis(q) * s 
          END DO
        END IF          
        FORCE(p) = FORCE(p) + beta * Basis(p) * s 
        FORCE(p+n) = FORCE(p+n) - beta * Basis(p) * s 
      END DO
    END DO

    ! In parallel case the contribution will come from both sides.
    ! Hence scale it by half and neglect pure halo contributions. 
    IF( ABS(AssFrac-1.0_dp) > TINY( AssFrac ) ) THEN     
      FORCE(1:2*n) = AssFrac * FORCE(1:2*n) 
      STIFF(1:2*n,1:2*n) = AssFrac * STIFF(1:2*n,1:2*n) 
    END IF
        
    ! We need our own caller since we may have switched order of parents
    ! This results to the need to have our own scaling in parallel.
    CALL UpdateGlobalEquations( Solver % Matrix, STIFF, Solver % Matrix % rhs, FORCE, &
        2*n, 1, DgIndexes(1:2*n), UElement=Element )

    BCDone = .TRUE.

!------------------------------------------------------------------------------
  END SUBROUTINE LocalJumpsDisContBC
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
END SUBROUTINE HeatSolver
!------------------------------------------------------------------------------


