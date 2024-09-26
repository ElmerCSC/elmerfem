!------------------------------------------------------------------------------
!> Topology optimization workflow with the SIMP method
!> Based heavily on the ideas presented in Python code in:
!> "A 165 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013"
!------------------------------------------------------------------------------
SUBROUTINE TopoOpt_init0( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  INTEGER :: i
  
  Params => Solver % Values

  CALL ListAddNewInteger( Params,'Primary Solver Index', 1 )  
  i = ListGetInteger( Params,'Primary Solver Index' )
  IF(ListGetLogical( Params,'Solve Adjoint Problem', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Adjoint')
  END IF
    
!------------------------------------------------------------------------------
END SUBROUTINE TopoOpt_init0
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE TopoOpt_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  Params => Solver % Values

  ! These automatically allocate elemental variables that are then created by library
  ! even before visiting the subroutine below. 
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo rho' )
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo mult' )
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo ce' )
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo dc' )

  IF( ListGetLogical( Params,'Create BW Topology', Found ) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'-elem topo bw' )
  END IF
  
  ! Add a global variable to store the norm to.
  CALL ListAddNewString( Params,'Variable','-nooutput -global topoopt_nrm')
    
!------------------------------------------------------------------------------
END SUBROUTINE TopoOpt_init
!------------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!> Solver for topology optimization. It is assumed that it is paired with a
!> generic elasticity equation.
!------------------------------------------------------------------------------
SUBROUTINE TopoOpt( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found 
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: pVar, uVar, aVar
  TYPE(Matrix_t), POINTER :: Fmat 
  REAL(KIND=dp), ALLOCATABLE :: local_sol_array(:,:), local_sol(:)
  REAL(KIND=dp), POINTER :: ce(:), dc(:), dv(:), dv0(:), bw(:), xTopo(:), xPhys(:), xMult(:)
  INTEGER :: TimesVisited = 0, dim, dofs, Niter, i, j, n, m, Nelems, Nnodes, nsize, cMode
  REAL(KIND=dp) :: volFrac, penal, emin, efrac, gt, obj, val, wmin, Diff(3)
  TYPE(Solver_t), POINTER :: PhysSolver
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  LOGICAL :: PdeFilter, SimpleFilter, Csymmetry, SolveAdj
  INTEGER, POINTER :: ElemPerm(:)
  CHARACTER(:), ALLOCATABLE :: filterMethod, filterType
  CHARACTER(*), PARAMETER :: Caller = 'TopoOpt'

  
  SAVE :: TimesVisited, Fmat, xTopo, xPhys, xMult, Niter, PhysSolver, dim, Mesh, &
      local_sol_array, local_sol, ce, dc, dv, dv0, bw, wmin, FilterMethod, FilterType, &
      gt, Nnodes, Nelems, uVar, aVar, dofs, Nodes, PdeFilter, SimpleFilter, Diff, nsize, &
      ElemPerm, Csymmetry, SolveAdj, obj
  
  
  CALL Info(Caller,'-------------------------------------------')
  CALL Info(Caller,'Updating density for topology optimization')
  
  ! Note: Physical problem should be solved when we come here.
  
  Params => Solver % Values

  IF( TimesVisited == 0) THEN
    Mesh => Solver % Mesh
    dim = Mesh % MeshDim
    Nnodes = Mesh % NumberOfNodes
    Nelems = Mesh % NumberOfBulkElements

    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric )
    
    CALL Info(Caller,'Number of nodes: '//I2S(Nnodes))
    CALL Info(Caller,'Number of bulk elements: '//I2S(Nelems))

    filterMethod = ListGetString( Params,'Filter Method',UnfoundFatal=.TRUE.)
    CALL Info(Caller,'Filter Method: '//TRIM(filterMethod))

    filterType = ListGetString( Params,'Filter Type',UnfoundFatal=.TRUE.)
    CALL Info(Caller,'Filter Type: '//TRIM(filterType))

    PdeFilter = .FALSE.
    SimpleFilter = .FALSE.
    SELECT CASE ( filterType )
    CASE('simple')
      SimpleFilter = .TRUE.
    CASE('pde')
      PdeFilter = .TRUE.
    CASE('distance')
    CASE DEFAULT
      CALL Fatal(Caller,'Unknown filter type: '//TRIM(filterType))
    END SELECT

    IF(PdeFilter .AND. .NOT. ASSOCIATED(Solver % Matrix)) THEN
      CALL Fatal(Caller,'Pde Filter requires field variable & matrix to be associated!')
    END IF
    
    ! This is not generic. We assume the stress solver to be the 1st solver for now. 
    i = ListGetInteger( Params,'Primary Solver Index' )
    PhysSolver => Model % Solvers(i)
    uVar => PhysSolver % Variable
    dofs = uVar % dofs

    ! Get adjoint variable, or if not requested use primary variable for adjoint as well. 
    SolveAdj = ListGetLogical( Params,'Solve Adjoint Problem', Found )
    IF( SolveAdj ) THEN
      aVar => VariableGet( Mesh % Variables,'Adjoint' )
      IF(.NOT. ASSOCIATED(aVar)) CALL Fatal(Caller,'Could not find adjoint variable!')
    ELSE
      aVar => uVar
    END IF
        
    ! These fields are created also for visualization in mind!
    pVar => VariableGet( Mesh % Variables,"topo rho", UnfoundFatal = .TRUE.)
    xPhys => pVar % Values
    ElemPerm => pVar % Perm
    nsize = SIZE(xPhys)
    CALL Info(Caller,'Size of elemental fields: '//I2S(nsize),Level=10)

    pVar => VariableGet( Mesh % Variables,"topo mult", UnfoundFatal = .TRUE.)
    xMult => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo ce", UnfoundFatal = .TRUE.)
    ce => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo dc", UnfoundFatal = .TRUE.)
    dc => pVar % Values
    pVar => VariableGet( Mesh % Variables,"topo bw")
    IF(ASSOCIATED(pVar)) THEN
      bw => pVar % Values
    ELSE
      bw => NULL()
    END IF
    
    ! Allocate full vectors
    ALLOCATE( xTopo(nsize) )
    xTopo = xPhys
    gt = 0.0_dp
    ce = 0.0_dp
    dc = 0.0_dp
            
    ! Calculate elemental volume   
    ALLOCATE( dv0(nsize) ) 
    dv0 = 0.0_dp
    DO i=1,Nelems
      j = ElemPerm(i)
      IF(j==0) CYCLE
      Element => Mesh % Elements(i)
      CALL CopyElementNodesFromMesh( Nodes, Solver % Mesh, &
          Element % TYPE % NumberOfNodes, Element % NodeIndexes )
      dv0(j) = ElementSize( Element, Nodes )      
      IF(Csymmetry) dv0(j) = dv0(j) * SUM(Nodes % x) / Element % TYPE % NumberOfNodes
    END DO
    IF(InfoActive(20)) THEN
      CALL VectorValuesRange(dv0,nsize,'dv0')       
    END IF
    
    ! Allocate elemental stuff
    n = Mesh % MaxElementDofs        
    ALLOCATE(local_sol_array(dofs,n), local_sol(dofs*n))
    
    wmin = ListGetConstReal( Params,'Sensitivity Filter Threshold', Found )
    IF(.NOT. Found) wmin = 1.0e-3

    IF(PdeFilter ) THEN      
      BLOCK 
        REAL(KIND=dp), POINTER :: HWrk(:,:)
        INTEGER :: d1, d2
        Hwrk => ListGetConstRealArray( Params,'PDE Filter Diffusion Constant',UnfoundFatal = .TRUE.)

        d1 = SIZE(Hwrk,1)
        d2 = SIZE(Hwrk,2)
        
        IF (d1 == 1 .AND. d2 == 1 ) THEN
          Diff = Hwrk( 1,1 )
        ELSE IF(d1 == 1 .AND. d2 >= dim ) THEN
          Diff(1:dim) = Hwrk(1,1:dim)
        ELSE IF(d1 >= dim .AND. d2 == 1 ) THEN
          Diff(1:dim) = Hwrk(1:dim,1)
        ELSE
          CALL Fatal(Caller,'Invalid size for "PDE Filter Diffusion Constant": '//I2S(d1)//' x '//I2S(d2))
        END IF        
      END BLOCK
    ELSE
      IF(ParEnv % PEs > 1 ) THEN
        CALL Fatal(Caller,'Only PDE Filter is implemented in parallel!')        
      END IF
      IF( SimpleFilter ) THEN
        FMat => CreateSimpleFilter()
        Niter = MAX(1,ListGetInteger( Params,'Simple Filter Iterations', Found ))
      ELSE
        Fmat => CreateDistanceFilter()
        Niter = 1
      END IF      
      CALL NormalizeFilter(Fmat,.TRUE.)
      val = 1.0_dp * SIZE(Fmat % Cols) / Fmat % NumberOfRows 
      WRITE(Message,'(A,ES12.3)') 'Average number of hits in filter:',val
      CALL Info(Caller,Message)
    END IF
      
    SELECT CASE( FilterMethod )
    CASE('sensitivity')
      dv => dv0 
    CASE('density')
      ALLOCATE(dv(nsize))
      IF( PdeFilter ) THEN
        CALL ApplyPdeFilter( dv0, dv, Diff )
      ELSE
        CALL ApplyFilter( Fmat, dv0, dv, niter, Trans=.TRUE. )
      END IF      
      CALL VectorValuesRange(dv,SIZE(dv),'dv')       
    CASE('none')
      dv => dv0
    CASE DEFAULT
      CALL Fatal(Caller,'Uknown filtering method: '//TRIM(FilterMethod))
    END SELECT
  END IF  ! TimesVisited==0
  
  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(PhysSolver % Matrix % Values,SIZE(PhysSolver % Matrix % Values),'Kmat')       
    CALL VectorValuesRange(uVar % Values,SIZE(uVar % Values),TRIM(uVar % Name))
  END IF

  ! These parameters can depend on time etc. 
  penal = ListGetCReal( Params,'Penalty Exponent',UnfoundFatal=.TRUE.)
  volFrac = ListGetCReal( Params,'Volume Fraction',UnfoundFatal=.TRUE.)
  emin = ListGetCReal( Params,'Minimum Relative Density',UnfoundFatal=.TRUE.)  
  efrac = 1.0_dp - emin

  ! Go to internal density interval [0.0,1.0] (used by the reference code)
  ! xPhys = (xPhys-emin)/efrac
  IF( TimesVisited == 0 ) THEN
    xTopo = xPhys
  END IF
  
  ! Gradients/Sensitivities with respect to the SIMP objective function and
  ! the volume constraint.       
  cMode = GetCycleMode()
  ! 0 - normal
  ! 1 - init cycle
  ! 2 - mid cycle
  ! 3 - end cycle
  
  IF( cMode == 0 .OR. cMode == 1 ) THEN
    obj = 0.0_dp
    dc = 0.0_dp
  END IF
    
  CALL ObjectiveGradients(xPhys,ce,dc,dv0,obj)

  IF( cMode == 1 .OR. cMode == 2 ) THEN
    CALL Info(Caller,'Mid of cycle, finishing early!')
    GOTO 1
  END IF
      
  obj = ParallelReduction( obj ) 

  
  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(xPhys,SIZE(xPhys),'xPhys')       
    CALL VectorValuesRange(ce,SIZE(ce),'ce')       
    CALL VectorValuesRange(dc,SIZE(dc),'dc')       
  END IF
  WRITE(Message,*) 'Objective function: ',obj
  CALL Info(Caller,Message)
  
  ! Pre-filter  
  SELECT CASE( FilterMethod )
  CASE('sensitivity')
    IF( PdeFilter ) THEN
      CALL ApplyPdeFilter( dc, dc, Diff, xTopo, wmin )
    ELSE
      CALL ApplyFilter( Fmat, dc, dc, niter, xTopo, wmin )
    END IF
    IF(InfoActive(20)) THEN
      CALL VectorValuesRange(dc,SIZE(dc),'dc pre')       
    END IF
  CASE('density')
    IF( PdeFilter ) THEN
      CALL ApplyPDEFilter( dc, dc, Diff )
    ELSE
      CALL ApplyFilter( Fmat, dc, dc, niter, Trans=.TRUE. )
    END IF
    IF(InfoActive(20)) THEN
      CALL VectorValuesRange(dc,SIZE(dc),'dc pre')       
    END IF
  CASE('none')
    CALL Info(Caller,'Applying no filtering')
  CASE DEFAULT
    CALL Fatal(Caller,'Uknown filtering method: '//TRIM(FilterMethod))
  END SELECT

  CALL UpdateDensities(xTopo,dc,dv,gt)

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(xTopo,SIZE(xTopo),'xTopo')       
  END IF

  ! Post-filter
  SELECT CASE( FilterMethod )
  CASE('density')
    IF( PdeFilter ) THEN
      CALL ApplyPDEFilter( xTopo, xPhys, Diff )
    ELSE
      CALL ApplyFilter( Fmat, xTopo, xPhys, niter )
    END IF
  CASE DEFAULT
    xPhys = xTopo
  END SELECT

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(xPhys,SIZE(xPhys),'xPhys')       
    CALL VectorValuesRange(ABS(xPhys-xTopo),SIZE(xPhys),'dx')       
  END IF
  
  ! We may pass the objective function as the norm to control convergence
  IF(SIZE(Solver % Variable % Values) == 1 ) THEN
    Solver % Variable % Values = obj
  END IF

  IF(ASSOCIATED(bw)) THEN
    CALL DefineTopologyBW(dv0,xPhys,bw)    
  END IF

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(xPhys,SIZE(xPhys),'xPhys2')       
  END IF
  
  ! Multiplier for local stiffness matrix of the external solver.
1 xMult = emin + efrac * xPhys**penal
    
  TimesVisited = TimesVisited + 1
  
CONTAINS

  FUNCTION GetCycleMode() RESULT ( Mode )
    INTEGER :: Mode
    INTEGER :: nCycle, nT

    Mode = 0
    nCycle = ListGetInteger( Model % Simulation,'Periodic Timesteps',Found )
    nT = TimesVisited 

    ! 1st iteration make a something to create asymmetry
    IF(nCycle == 0 .OR. nT == 0 ) RETURN

    SELECT CASE( MODULO(nT,nCycle) )
    CASE( 0 )
      Mode = 3
    CASE( 1 )
      Mode = 1
    CASE DEFAULT
      Mode = 2
    END SELECT
    
  END FUNCTION GetCycleMode


  
  !---------------------------------------------------------------------
  !> xphys: in [0,1], densities used for scaling the material properties
  !> ce: elementsize energies before scaling by density
  !> dc: gradient with respect to objective function.
  !> dv: gradient with respect to volume constraint.
  !---------------------------------------------------------------------
  SUBROUTINE ObjectiveGradients(x,ce,dc,dv,obj) 
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp), POINTER :: ce(:),dc(:),dv(:)
    REAL(KIND=dp) :: obj

    INTEGER :: i,j,k
    REAL(KIND=dp), ALLOCATABLE:: Stiff(:,:), Force(:)

    
    n = Solver % Mesh % MaxElementNodes * dofs
    ALLOCATE(Stiff(n,n), Force(n) )
    
    DO i=1,Solver % NumberOfActiveElements
      Element => Mesh % Elements(Solver % ActiveElements(i))
      n = Element % TYPE % NumberOfNodes    
      m = dofs*n

      ! Get the local stiffness matrix as saved by the primary solver
      CALL GetLocalMatrixStorage( PhysSolver, m, Stiff, Force, Found, ActiveInd = i )       
      IF(.NOT. Found) CALL Fatal(Caller,'Could not find local stiffness matrix!')

      ! Get the solution from stress solver 
      IF(dofs == 1) THEN
        CALL GetLocalSolution( local_sol,UElement=Element,USolver=PhysSolver) 
      ELSE
        CALL GetLocalSolution( local_sol_array,UElement=Element,USolver=PhysSolver) 
        local_sol(1:m) = RESHAPE( local_sol_array(1:dofs,1:n), [m] )
      END IF

      ! Elemental energy assuming unity multiplier.
      ce(i) = SUM( local_sol(1:m) * MATMUL( Stiff, local_sol(1:m) ) )            
    END DO

    !PRINT *,'Objective:',obj,emin,efrac,penal,SUM(x),SUM(ce)
    
    ! Derivative of elemental energy
    dc = dc - penal*x**(penal-1) * efrac * ce

    ! Objective function
    obj = obj + SUM( (emin + efrac*( x**penal ) ) * ce )
    
  END SUBROUTINE ObjectiveGradients
    

  !--------------------------------------------------------------------------
  !> Optimality criteria method (section 2.2 in paper) for maximum/minimum 
  !> stiffness/compliance. Heuristic updating scheme for the element densities 
  !> to find the Lagrangian multiplier.
  !--------------------------------------------------------------------------
  SUBROUTINE UpdateDensities(x,dc,dv,g)

    REAL(KIND=dp), POINTER :: x(:), dc(:), dv(:)
    REAL(KIND=dp) :: g

    REAL(KIND=dp), ALLOCATABLE :: xnew(:)
    REAL(KIND=dp) :: Vi, l1, l2, lmid, move, err, tol, V0
    INTEGER :: k
    LOGICAL :: Visited = .FALSE.

    SAVE move, tol
    
    l1 = 0.0_dp
    l2 = 1.0e9_dp
    ! maximum update of density

    move = ListGetCReal(Params,'Bisection search max change',Found )
    IF(.NOT. Found) move = 0.2    
    
    tol = ListGetCReal(Params,'Bisection search tolerance',Found )
    IF(.NOT. Found) tol = 1.0e-6
      
    ! Desired total volume
    V0 = volFrac * SUM(dv)
    V0 = ParallelReduction(V0)
    
    ALLOCATE(xnew(SIZE(x)))
    xnew = 0.0_dp
        
    DO k=1,1000            
      lmid = 0.5*(l2+l1)
      
      ! Note: xnew in [0,1]
      ! Suggested new density
      xnew = x*SQRT(-dc/(dv*lmid))

      ! Regulators and limiters
      xnew = MAX(0.0_dp,MAX(x-move,MIN(1.0_dp,MIN(x+move,xnew))))

      ! Volume balance should become zero!
      Vi = SUM(dv*xnew)
      Vi = ParallelReduction(Vi)
     
      IF (Vi > V0) THEN
        l1 = lmid
      ELSE
        l2 = lmid
      END IF
      
      err = (l2-l1)/(l1+l2)
      IF( err < tol ) EXIT      

      IF( InfoActive(20)) THEN
        PRINT *,'Bisection:',k,Vi,l1,l2,err
      END IF
    END DO

    x = xnew 
    g = Vi - V0
    CALL Info(Caller,'Number of bisection iterations: '//I2S(k))
    
  END SUBROUTINE UpdateDensities
                

  !------------------------------------------------------------------
  !> Applies a filter given by CRS matrix with rowsum scaled to unity.
  !> Optionally use a weight vector before and after scaling.
  !------------------------------------------------------------------
  SUBROUTINE ApplyFilter( Fmat, x, xf, niter, w, wmin, Trans )
    TYPE(Matrix_t), POINTER :: Fmat
    REAL(KIND=dp), POINTER :: x(:), xf(:)
    INTEGER, OPTIONAL :: niter
    REAL(KIND=dp), POINTER, OPTIONAL :: w(:)
    REAL(KIND=dp), OPTIONAL :: wmin
    LOGICAL, OPTIONAL :: Trans
    
    REAL(KIND=dp), ALLOCATABLE :: xtmp(:)
    REAL(KIND=dp), POINTER :: SValues(:)
    INTEGER :: n, m, i, j
    LOGICAL :: DoTrans 

    m = 1
    IF(PRESENT(niter)) m = niter

    DoTrans = .FALSE.
    IF(PRESENT(Trans)) DoTrans = Trans

    ALLOCATE(xtmp(SIZE(x)))
    
    IF( PRESENT(w)) THEN
      xtmp = x*w
      IF(.NOT. PRESENT(wmin)) THEN
        CALL Fatal(Caller,'If we have weight we need "wmin" as well!')
      END IF
      CALL Info(Caller,'Scaling filter with weight!',Level=20)
    ELSE
      xtmp = x
    END IF

    DO i=1,m
      IF( DoTrans ) THEN
        CALL TransposeMatrixVectorMultiply( Fmat, xtmp, xf)      
      ELSE
        CALL MatrixVectorMultiply( Fmat, xtmp, xf)
      END IF
      IF(i<m) THEN
        xtmp = xf
      END IF      
    END DO
    
    IF( PRESENT(w)) THEN
      xf = xf/MAX(w,wmin)
    END IF
    DEALLOCATE(xtmp)      

  END SUBROUTINE ApplyFilter


  !------------------------------------------------------------------
  !> Applies a ave/min/max filter that is determined by the topology of
  !> the filter matrix but no values are used.
  !------------------------------------------------------------------
  SUBROUTINE ApplyTopologyFilter( Fmat, mode, x, xf, niter )
    TYPE(Matrix_t), POINTER :: Fmat
    INTEGER :: mode
    REAL(KIND=dp), POINTER :: x(:), xf(:)
    INTEGER, OPTIONAL :: niter

    REAL(KIND=dp), POINTER :: xtmp(:)
    INTEGER, POINTER :: Cols(:), Rows(:)
    INTEGER :: n, m, i, j

    m = 1
    IF(PRESENT(niter)) m = niter

    n = Fmat % NumberOfRows
    Rows => Fmat % Rows
    Cols => Fmat % Cols

    IF(m>1) THEN
      ALLOCATE(xtmp(n))
      xtmp = x
    ELSE
      xtmp => x
    END IF
        
    DO j=1,m                
      DO i=1,n
        SELECT CASE( mode )
        CASE( 0 )
          xf(i) = SUM(xtmp(Cols(Rows(i):Rows(i+1)-1))) / (Rows(i+1)-Rows(i))
        CASE( 1 )
          xf(i) = MINVAL(xtmp(Cols(Rows(i):Rows(i+1)-1)))
        CASE( 2 )
          xf(i) = MAXVAL(xtmp(Cols(Rows(i):Rows(i+1)-1)))
        END SELECT
      END DO
      IF(j<m) xtmp = xf
    END DO
      
    IF(m>1) DEALLOCATE(xtmp)
    
  END SUBROUTINE ApplyTopologyFilter


  
  
  !-------------------------------------------------------------------
  !> Normalize the entries such that the rowsum (or columnsum) is unity
  !-------------------------------------------------------------------
  SUBROUTINE NormalizeFilter(A,TransNorm)
    TYPE(Matrix_t), POINTER :: A
    LOGICAL :: TransNorm

    INTEGER :: i,j,k,n
    REAL(KIND=dp), ALLOCATABLE :: colsum(:)
    REAL(KIND=dp) :: rsum

    n = A % NumberOfRows 

#if 0
    IF( TransNorm ) THEN
      CALL Info('NormalizeFilter','Normalizing filter by columnsum!')
      ! First calculate the colsum
      ALLOCATE(colsum(n))
      colsum = 0.0_dp
      DO i=1, n
        DO j = A % Rows(i), A % Rows(i+1)-1
          k = A % Cols(j)
          colsum(k) = colsum(k) + A % Values(j)
        END DO
      END DO
      
      ! Now create the transposed normalized projector
      IF(.NOT. ASSOCIATED(A % TValues)) THEN
        ALLOCATE(A % TValues(SIZE(A % Values)))
        A % TValues = 0.0_dp
      END IF
      DO i=1, A % NumberOfRows
        DO j = A % Rows(i), A % Rows(i+1)-1      
          k = A % Cols(j)
          A % TValues(j) = A % Values(j) / colsum(k)
        END DO
      END DO
    END IF
#endif
    
    ! Then create the standard projector normalized by rowsum
    CALL Info('NormalizeFilter','Normalizing filter by rowsum!')
    DO i=1, A % NumberOfRows
      rsum = 0.0_dp
      DO j = A % Rows(i), A % Rows(i+1)-1
        rsum = rsum + A % Values(j)
      END DO
      DO j = A % Rows(i), A % Rows(i+1)-1      
        A % Values(j) = A % Values(j) / rsum
      END DO
    END DO
    
  END SUBROUTINE NormalizeFilter


  !----------------------------------------------------------------------------
  !> Create filter that inclues just closest neighbours associated attached by
  !> faces (3D) or edges (2D). This has rather small support and needs to be
  !> typically applied several times. 
  !----------------------------------------------------------------------------
  FUNCTION CreateSimpleFilter() RESULT ( Emat ) 
    TYPE(Matrix_t), POINTER :: Emat    

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Face, ElemA, ElemB
    INTEGER :: i,j,k,k1,k2,NoElems,kcum(27)
    
    CALL Info(Caller,'Creating filter based on element-to-element connectivity')
    Mesh => Solver % Mesh
    IF( Mesh % MeshDim == 3 ) THEN
      CALL FindMeshFaces3D(Mesh)
      NoElems = Mesh % NumberOfFaces
    ELSE
      CALL FindMeshEdges2D(Mesh)
      NoElems = Mesh % NumberOfEdges
    END IF

    ! Create sparse matrix for element-to-element connectivity
    Emat => AllocateMatrix()
    Emat % FORMAT = MATRIX_LIST

    ! Add the max index first because list matrix likes this
    i = Mesh % NumberOfBulkElements
    CALL List_AddToMatrixElement( EMat % ListMatrix,i,i,0.0_dp )

    DO i=1, NoElems 
      IF( Mesh % MeshDim == 3 ) THEN
        Face => Mesh % Faces(i)
      ELSE
        Face => Mesh % Edges(i)
      END IF
      ElemA => Face % BoundaryInfo % Left
      ElemB => Face % BoundaryInfo % Right
      IF(.NOT. ASSOCIATED(ElemA) .OR. .NOT. ASSOCIATED(ElemB) ) CYCLE

      k1 = ElemA % ElementIndex
      k2 = ElemB % ElementIndex

      ! Each neighbour gets the same weight
      CALL List_AddToMatrixElement( Emat % ListMatrix,k1,k2,1.0_dp )
      CALL List_AddToMatrixElement( Emat % ListMatrix,k2,k1,1.0_dp )

      ! Set diagonals too. This way the filter has 0.5 weight for itself. 
      CALL List_AddToMatrixElement( Emat % ListMatrix,k1,k1,1.0_dp )
      CALL List_AddToMatrixElement( Emat % ListMatrix,k2,k2,1.0_dp )
    END DO

    ! Go from list matrix to more efficient CRS matrix
    CALL List_ToCRSMatrix(Emat)

    IF(InfoActive(10)) THEN
      k1 = HUGE(k1)
      k2 = 0
      kcum = 0
      DO i=1, Emat % NumberOfRows
        k = Emat % Rows(i+1) - Emat % Rows(i)
        k1 = MIN(k1, k)
        k2 = MAX(k2, k)
        kcum(k) = kcum(k) + 1
      END DO
      DO i=1,k2
        IF(kcum(i)>0) PRINT *,'Cumulative hits:',i,kcum(i)
      END DO
    END IF
      
    CALL Info(Caller,'Number of hits range for filter ['//I2S(k1)//','//I2S(k2)//']')    
    CALL Info(Caller,'Number of rows in filter: '//TRIM(I2S(Emat % NumberOfRows)))
    CALL Info(Caller,'Number of non-zeros in filter: '//TRIM(I2S(SIZE(Emat % Values))))
    
  END FUNCTION CreateSimpleFilter


  !------------------------------------------------------------------------------------------------
  !> Create filter that includes all elements witing distance smaller than "rmin" between elements.
  !> We use the connectivity of simple filter to find the candidate elements. 
  !------------------------------------------------------------------------------------------------
  FUNCTION CreateDistanceFilter() RESULT ( Rmat ) 
    TYPE(Matrix_t), POINTER :: Rmat

    TYPE(Matrix_t), POINTER :: Emat
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: ElemCenters(:,:)
    REAL(KIND=dp) :: rfilter, rfilter2, rik2
    INTEGER :: NoElems,i,i2,j,k,k1,k2,k3,n,kmax,kmin
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: Inds(:),kcum(:)    
    REAL(KIND=dp), ALLOCATABLE :: Dist(:)
    
    CALL Info(Caller,'Creating filter based on element-to-element distance')
    Mesh => Solver % Mesh

    CALL ResetTimer('DistanceFilter')

    dim = Mesh % MeshDim
    rfilter = ListGetCReal(Solver % Values,'Distance Filter Radius', UnfoundFatal = .TRUE.)     
    rfilter2 = rfilter**2
    NoElems = Mesh % NumberOfBulkElements

    n = 1000
    ALLOCATE(Inds(n), Dist(n) ) 
    Inds = 0
    Dist = 0.0_dp
    
    Emat => CreateSimpleFilter()
    
    ! Compute center of elements for speedier distance computation.
    ALLOCATE(ElemCenters(dim,NoElems))
    DO i=1,NoElems
      Element => Mesh % Elements(i) 
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      ElemCenters(1,i) = SUM(Mesh % Nodes % x(NodeIndexes)) / n
      ElemCenters(2,i) = SUM(Mesh % Nodes % y(NodeIndexes)) / n
      IF(dim==3) ElemCenters(3,i) = SUM(Mesh % Nodes % z(NodeIndexes)) / n
    END DO
    
    ! Create sparse matrix for element-to-element connectivity
    Rmat => AllocateMatrix()
    Rmat % FORMAT = MATRIX_LIST

    ! Add the max index first because list matrix likes this
    CALL List_AddToMatrixElement( RMat % ListMatrix,NoElems,NoElems,0.0_dp )

    kmax = 0
    kmin = HUGE(kmin)
    
    DO i=1, NoElems
      k1 = 1
      k3 = 1
      Inds(1) = i
      Dist(1) = rfilter
      DO WHILE(.TRUE.)
        k2=k3
        DO k=k1,k2
          DO j=Emat % Rows(Inds(k)), Emat % Rows(Inds(k)+1)-1
            i2 = Emat % Cols(j)            
            IF(ANY(Inds(1:k3) == i2)) CYCLE                          

            ! square of distance between element centers
            rik2 = SUM((ElemCenters(:,i)-ElemCenters(:,i2))**2)

            ! Check whether the centerpoints are within filter radius
            IF( rik2 < rfilter2 ) THEN
              k3 = k3 + 1
              IF(k3 > SIZE(Inds)) CALL Fatal(Caller,'Too many neighbours!')
              Inds(k3) = i2
              Dist(k3) = rfilter-SQRT(rik2)
            END IF
          END DO
        END DO
        ! We found no new elements within radius
        IF(k3 == k2) EXIT
        
        ! We have tested neighbours for 'k2' elements already
        k1 = k2+1
      END DO

      ! Assemble the found ones in a new row.
      ! For List_Add the reverse ordering seems to be faster.
      CALL SortF(k3,Inds,Dist)

      DO k1=1,k3
        DO k=k1+1,k2          
          IF(Inds(k1) == Inds(k) ) CALL Fatal(Caller,'Duplicate indeces when creating distance filter!')
        END DO
      END DO

      DO k1=k3,1,-1
        CALL List_AddToMatrixElement( RMat % ListMatrix,i,Inds(k1),Dist(k1))
      END DO     
              
      kmax = MAX(kmax,k3)
      kmin = MIN(kmin,k3)
    END DO
    
    ! Go from list matrix to more efficient CRS matrix
    CALL List_ToCRSMatrix(Rmat)       
    
    ! We do not need the element-to-element connectivity any more. 
    CALL FreeMatrix(Emat)

    IF(InfoActive(10)) THEN
      ALLOCATE(kcum(kmax))
      kcum = 0
      DO i=1, Rmat % NumberOfRows
        k = Rmat % Rows(i+1) - Rmat % Rows(i)
        kcum(k) = kcum(k) + 1
      END DO
      DO i=1,kmax
        IF(kcum(i)>0) PRINT *,'Cumulative hits:',i,kcum(i)
      END DO
    END IF
      
    CALL Info(Caller,'Number of hits range for filter ['//I2S(kmin)//','//I2S(kmax)//']')
    CALL Info(Caller,'Number of rows in filter: '//TRIM(I2S(Rmat % NumberOfRows)))
    CALL Info(Caller,'Number of non-zeros in filter: '//TRIM(I2S(SIZE(Rmat % Values))))
    CALL CheckTimer(Caller,Delete=.TRUE.)
    
  END FUNCTION CreateDistanceFilter


!------------------------------------------------------------------------------
!> Given a topology xPhys create a 0/1 presentation that conserves volume. 
!------------------------------------------------------------------------------
  SUBROUTINE DefineTopologyBW(dv0,xPhys,bw)    
    REAL(KIND=dp), POINTER :: dv0(:), xPhys(:), bw(:)

    REAL(KIND=dp) :: xlow, xup, xmid, h
    REAL(KIND=dp), ALLOCATABLE :: histv(:), cumv(:)

    REAL(KIND=dp) :: Vtot, Vtarget
    INTEGER :: i,j,k,m,iter
    LOGICAL :: Hit
    
    m = 100
    ALLOCATE(histv(0:m),cumv(0:m))
    
    xlow = 0.0_dp
    xup = 1.0_dp

    Vtot = SUM(dv0(1:nsize))
    Vtarget = volFrac * Vtot
    
    DO iter=1,10
      h = (xup-xlow) / m 
      histv = 0.0_dp
      cumv = 0.0_dp
      DO i=1,nsize
        j = MAX(0,MIN(CEILING((xPhys(i)-xlow)/h),m))
        histv(j) = histv(j) + dv0(i)
      END DO
      cumv(0) = histv(0)

      ! For parallel runs communicate the histogram here.
      DO i=1,m        
        cumv(i) = cumv(i-1) + histv(i)
      END DO

      Hit = .FALSE.
      DO i=1,m
        IF(cumv(i-1) < Vtarget .AND. cumv(i) > Vtarget) THEN
          xlow = xlow + (i-1)*h
          xup = xlow + h
          Hit = .TRUE.
          EXIT
        ELSE IF(ABS(cumv(i-1)-Vtarget) < EPSILON(h)) THEN
          xlow = xlow + (i-1)*h
          xup = xlow
          EXIT
        ELSE IF(ABS(cumv(i)-Vtarget) < EPSILON(h)) THEN
          xlow = xlow + i*h
          xup = xlow          
          EXIT
        END IF
      END DO
      IF(.NOT. Hit) EXIT
    END DO
    
    xmid = (xlow+xup)/2.0_dp
    WHERE(xPhys > xmid )
      bw = 1.0_dp
    ELSE WHERE
      bw = 0.0_dp
    END WHERE

    IF(InfoActive(7)) THEN
      PRINT *,'Mass Conserving Limit:',xmid,xlow,xup,iter
    END IF
    
  END SUBROUTINE DefineTopologyBW


!------------------------------------------------------------------------------
!> Assembly of the matrix equation used for PDE filtering.
!> We may assembly both matrix and r.h.s., or just the r.h.s.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, DoMatrix, Diff, x ) 
!------------------------------------------------------------------------------
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: DoMatrix
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp) :: Diff(3)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight, D, xi
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: STIFF(n,n), FORCE(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    D = 1.0_dp
    
    CALL GetElementNodes( Nodes )
    
    ! Separate matrix and force vector integration because we may use lower order
    ! integration scheme for the force vector.
    !----------------------------------------------------------------------------
    IF(DoMatrix ) THEN
      STIFF = 0._dp
      IP = GaussPoints( Element )
      DO t=1,IP % n
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )

        Weight = IP % s(t) * DetJ
        IF(Csymmetry) weight = Weight * SUM(Basis(1:n) * Nodes % x(1:n)) 
        
        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + Weight * ( Basis(p) * Basis(q) + & 
                SUM( Diff(1:dim) * dBasisdx(p,1:dim) * dBasisdx(q,1:dim) ) )
          END DO
        END DO
      END DO
    END IF

    FORCE = 0._dp

    IF( Csymmetry ) THEN
      ! We don't integrate accurately area with one gauss point for cylindrical coordinates. 
      IP = GaussPoints( Element )
    ELSE
      IP = GaussPoints( Element, np=1 )
    END IF
      
    xi = x(ElemPerm(Element % ElementIndex))

    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      IF(Csymmetry) weight = Weight * SUM(Basis(1:n) * Nodes % x(1:n))      

      FORCE(1:n) = FORCE(1:n) + Weight * Basis(1:n) * xi
    END DO

    IF( DoMatrix ) THEN
      CALL DefaultUpdateEquations(STIFF,FORCE)
    ELSE
      CALL DefaultUpdateForce(FORCE)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a diffusion-reaction equation to smooth down given values "x" to "xf".
!> PDE based filtering is ideal since it can use the parallel machinery of Elmer. 
!------------------------------------------------------------------------------
  SUBROUTINE ApplyPDEFilter(x, xf, Diff, w, wmin )
!------------------------------------------------------------------------------

    REAL(KIND=dp), POINTER :: x(:), xf(:)
    REAL(KIND=dp) :: Diff(3)
    REAL(KIND=dp), POINTER, OPTIONAL :: w(:)
    REAL(KIND=dp), OPTIONAL :: wmin
    
    REAL(KIND=dp), POINTER :: xtmp(:)
    INTEGER :: n, t, active
    LOGICAL :: DoMatrix = .TRUE.
    REAL(KIND=dp) :: Norm    

    ! Create weighted elemental field if requested. 
    IF( PRESENT(w)) THEN
      ALLOCATE(xtmp(SIZE(x)))
      xtmp = x*w      
      IF(.NOT. PRESENT(wmin)) THEN
        CALL Fatal(Caller,'If we have weight we need "wmin" as well!')
      END IF
      CALL Info(Caller,'Scaling filter with weight!',Level=20)
    ELSE
      xtmp => x
    END IF

    ! Assembly the matrix equation at the 1st time.
    ! Later just define the r.h.s. vector. 
    IF( DoMatrix ) THEN
      CALL DefaultInitialize()
    ELSE
      Solver % Matrix % Rhs = 0.0_dp
    END IF

    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      CALL LocalMatrix(  Element, n, DoMatrix, Diff, xtmp )
    END DO

    Norm = DefaultSolve()
    pVar => Solver % Variable 
    
    ! After solving the nodal values we need to transfer them back to elemental values. 
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      xf(t) = SUM(pVar % Values(pVar % Perm(Element % NodeIndexes)))/n
    END DO

    ! If weigting was used revert back. 
    IF( PRESENT(w)) THEN
      DEALLOCATE(xtmp)      
      xf = xf/MAX(w,wmin)
    END IF

    ! We have done the matrix. Freeze it and never touch it again.
    DoMatrix = .FALSE.
    
  END SUBROUTINE ApplyPDEFilter



!------------------------------------------------------------------------------
!> Solves a adjoint problem of the primary problem with different rhs.
!------------------------------------------------------------------------------
  SUBROUTINE SolveAdjointProblem(x)
!------------------------------------------------------------------------------

    REAL(KIND=dp), POINTER :: x(:)
    
    TYPE(Solver_t), POINTER :: pSolver
    REAL(KIND=dp), POINTER :: pRhs(:)
    

    INTEGER :: n, t, active
    LOGICAL :: DoMatrix = .TRUE.
    REAL(KIND=dp) :: Norm
    REAL(KIND=dp), POINTER :: aRhs(:)


    pSolver => Model % Solver
    pRhs => PhysSolver % Matrix % rhs 
       
    Model % Solver => PhysSolver
    PhysSolver % Variable => aVar

    IF( DoMatrix ) THEN
      ALLOCATE( aRhs(SIZE(aVar % Values)) )
      PhysSolver % Matrix % Rhs => aRhs
      aRhs = 0.0_dp
      Active = GetNOFActive()
      DO t=1,Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(t)
        n  = GetElementNOFNodes()
        !      CALL LocalMatrix(  Element, n, DoMatrix, Diff, xtmp )
      END DO
      ! We have done the rhs. No need to redo. 
      DoMatrix = .FALSE.
    ELSE
      PhysSolver % Matrix % Rhs => aRhs
    END IF

    ! Solve the adjoint problem with the same matrix equation is the primary matrix.  
    CALL ListAddLogical( PhysSolver % Values,'Skip Compute Change',.TRUE.)    
    Norm = DefaultSolve()
    CALL ListAddLogical( PhysSolver % Values,'Skip Compute Change',.FALSE.)
    
    ! Revert the saved pointers back
    Model % Solver => pSolver
    PhysSolver % Matrix % Rhs => pRhs
    PhysSolver % Variable => uVar
    
  END SUBROUTINE SolveAdjointProblem

  
!------------------------------------------------------------------------------
END SUBROUTINE TopoOpt
!------------------------------------------------------------------------------
