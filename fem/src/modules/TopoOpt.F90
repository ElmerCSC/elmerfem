!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
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
  LOGICAL :: Found, HaveField
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

  HaveField = ( ListGetString( Params,'Filter Type', Found ) == 'pde' ) 

  IF( HaveField ) THEN
    CALL ListAddNewString( Params,'Variable','xNodal') 
  ELSE
    ! Add a global variable to store the norm to if no variable present.
    CALL ListAddNewString( Params,'Variable','-nooutput -global topoopt_nrm')
  END IF

  IF( ListGetLogical( Params,'Create Zero Levelset', Found ) ) THEN
    IF( HaveField ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'topo levelset' )
    ELSE
      CALL Warn('TopoOpt_init','Only PDE filter can create levelset field!')
    END IF
  END IF
  
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
  USE MeshUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  LOGICAL :: Found 
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: pVar, uVar, aVar, pVarExt
  TYPE(Matrix_t), POINTER :: Fmat 
  REAL(KIND=dp), ALLOCATABLE :: local_sol_array(:,:), local_sol(:), local_act(:)
  REAL(KIND=dp), POINTER :: ce(:), dc(:), dv(:), dv0(:), bw(:), zeroset(:), xTopo(:), xPhys(:), xMult(:)
  INTEGER :: TimesVisited = 0, dim, dofs, Niter, i, j, n, m, Nelems, Nnodes, nsize, cMode
  REAL(KIND=dp) :: volFrac, penal, emin, efrac, gt, obj, val, wmin, Diff(3)
  TYPE(Solver_t), POINTER :: PhysSolver
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  LOGICAL :: PdeFilter, SimpleFilter, Csymmetry, SolveAdj, PhysSym, ElemField
  INTEGER, POINTER :: ElemPerm(:)
  INTEGER, ALLOCATABLE, SAVE :: SumPerm(:)
  LOGICAL, ALLOCATABLE, SAVE :: InterfaceNode(:)
  LOGICAL :: SkipInterface
  INTEGER :: nPer
  CHARACTER(:), ALLOCATABLE :: filterMethod, filterType
  CHARACTER(*), PARAMETER :: Caller = 'TopoOpt'

  
  SAVE :: TimesVisited, Fmat, xTopo, xPhys, xMult, Niter, PhysSolver, dim, Mesh, &
      local_sol_array, local_sol, local_act, ce, dc, dv, dv0, bw, zeroset, wmin, FilterMethod, FilterType, &
      gt, Nnodes, Nelems, uVar, aVar, dofs, Nodes, PdeFilter, SimpleFilter, Diff, nsize, &
      ElemPerm, Csymmetry, SolveAdj, obj, nPer, PhysSym, SkipInterface
  
  
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
    IF(.NOT. ListGetLogical(PhysSolver % Values,'Local Matrix Storage',Found ) ) THEN
      CALL Fatal(Caller,'Primary solver should have active "Local Matrix Storage"')
    END IF
      
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

    pVar => VariableGet( Mesh % Variables,"topo levelset")
    IF(ASSOCIATED(pVar)) THEN
      zeroset => pVar % Values
    ELSE
      zeroset => NULL()
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
    
    nPer = ListGetInteger( Params,'Periodic PhysSolver',Found ) 
    IF(nPer > 1) THEN
      ALLOCATE(SumPerm(Solver % Mesh % NumberOfBulkElements))
      PhysSym = ListGetLogical( Params,'Periodic PhysSolver Symmetric',Found )      
      ElemField = .TRUE.
      CALL RotationalPeriodicSumPerm(Solver, Solver % Mesh, 360.0_dp/nPer, &
          Solver % Variable % Perm, SumPerm, ElemField, PhysSym )
    END IF
    
    ! Allocate elemental stuff
    n = Mesh % MaxElementDofs        
    ALLOCATE(local_sol_array(dofs,n), local_sol(dofs*n), local_act(dofs*n))
    
    wmin = ListGetConstReal( Params,'Sensitivity Filter Threshold', Found )
    IF(.NOT. Found) wmin = 1.0e-3
    SkipInterface = .FALSE.
    
    IF(PdeFilter ) THEN      
      BLOCK 
        REAL(KIND=dp), POINTER :: HWrk(:,:)
        INTEGER, ALLOCATABLE :: NodeCount(:)
        INTEGER :: d1, d2, t
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

        SkipInterface = ListGetLogical( Params,'PDE Filter skip Interface',Found )
        IF( SkipInterface ) THEN
          ALLOCATE(NodeCount(Mesh % NumberOfNodes))
          NodeCount = 0
          DO t=1,Mesh % NumberOfBulkElements
            Element => Mesh % Elements(t)
            NodeCount(Element % NodeIndexes) = NodeCount(Element % NodeIndexes) + 1 
          END DO
          ALLOCATE(InterfaceNode(Mesh % NumberOfNodes))
          InterfaceNode = (NodeCount < 4 )
          DEALLOCATE(NodeCount)

          t = COUNT(InterfaceNode)
          CALL Info(Caller,'Number of Interface nodes: '//I2S(t),Level=7)
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
    CALL Info(Caller,'Extracting B&W coloring',Level=8) 
    CALL DefineTopologyBW(dv0,xPhys,bw)    
  END IF

  IF(ASSOCIATED(zeroset)) THEN
    CALL Info(Caller,'Extracting zero levelset function',Level=8) 
    CALL DefineTopologyZeroLevel(Solver % Variable,zeroset)    
  END IF

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(xPhys,SIZE(xPhys),'xPhys2')       
  END IF
  
  ! Multiplier for local stiffness matrix of the external solver.
1 CONTINUE

  xMult = emin + efrac * xPhys**penal

  IF( nPer > 1 ) THEN
    pVar => VariableGet( Mesh % Variables,"topo mult", UnfoundFatal = .TRUE.)
    pVarExt => VariableGet( Mesh % Variables,"topo mult extended", UnfoundFatal = .TRUE.)
    pVarExt % Values = 1.0_dp
    DO i=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(i)
      j = ABS(SumPerm(i))
      IF(j==0) CYCLE
      pVarExt % Values(pVarExt % Perm(i)) = pVar % Values(pVar % Perm(j))
    END DO
  END IF
    
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

    INTEGER :: i,j,k,l,NoModes, NoActive, sgn, sgn1, sgn2
    LOGICAL :: UseAdjoint, Found
    TYPE(Variable_t), POINTER :: AdjSol
    REAL(KIND=dp), ALLOCATABLE:: Stiff(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: spos, sneg
    
    n = Solver % Mesh % MaxElementNodes * dofs
    ALLOCATE(Stiff(n,n), Force(n) )

    UseAdjoint = ListGetLogical(Params,'Use Adjoint Solution',Found)
    IF(UseAdjoint) THEN
      AdjSol => VariableGet(Solver % Mesh % Variables,TRIM(PhysSolver % Variable % Name)//' adjoint')
      IF(.NOT. ASSOCIATED(AdjSol)) CALL Fatal(Caller,'Did not find Adjoint solution!')
      CALL Info(Caller,'Using adjoint solution: '//TRIM(AdjSol % Name))
      IF(InfoActive(20)) THEN
        PRINT *,'Adjoing interval:',MINVAL(AdjSol % Values), MAXVAL(AdjSol % Values), SIZE(AdjSol % Values)
      END IF
    END IF

    sgn1 = LIstGetInteger(Params,'Sign A',Found )
    IF(.NOT. Found) sgn1 = 1
    sgn2 = LIstGetInteger(Params,'Sign B',Found )
    IF(.NOT. Found) sgn2 = 1
    
    NoModes = ListGetInteger(Params,'No Modes',Found )
    IF(.NOT. Found ) THEN
      NoModes = PhysSolver % Variable % NumberOfConstraintModes 
    END IF

    IF(nPer > 1) THEN
      NoActive = PhysSolver % NumberOfActiveElements
    ELSE
      NoActive = Solver % NumberOfActiveElements
    END IF
    
    !ce = 0.0_dp
    !PRINT *,'NULLIFY ce'

    spos = 0.0_dp
    sneg = 0.0_dp
    
    DO i=1,NoActive
      IF(nPer > 1 ) THEN
        j = PhysSolver % ActiveElements(i)
        Element => Mesh % Elements(j)
        j = SumPerm(j)        
        IF(j==0) CYCLE
        IF(j<0) THEN
          j=-j
          sgn=-1
        ELSE
          sgn=1
        END IF
      ELSE
        j = Solver % ActiveElements(i)
        Element => Mesh % Elements(j)
      END IF
      l = ElemPerm(j)

      !IF(j>SIZE(ce)) PRINT *,'Too big j:',j,SIZE(ce)
      
      
      n = Element % TYPE % NumberOfNodes    
      m = dofs*n

      ! Get the local stiffness matrix as saved by the primary solver
      CALL GetLocalMatrixStorage( PhysSolver, m, Stiff, Force, Found, &
          ElemInd = Element % ElementIndex ) 
      IF(.NOT. Found) CALL Fatal(Caller,'Could not find local stiffness matrix!')

      ! Get the solution from stress solver 
      IF(dofs == 1) THEN
        CALL GetLocalSolution( local_sol,UElement=Element,USolver=PhysSolver) 
      ELSE
        CALL GetLocalSolution( local_sol_array,UElement=Element,USolver=PhysSolver) 
        local_sol(1:m) = RESHAPE( local_sol_array(1:dofs,1:n), [m] )
      END IF

      local_act(1:m) = MATMUL( Stiff(1:m,1:m), local_sol(1:m) )
            
      ! Elemental energy assuming unity multiplier.
      IF(UseAdjoint) THEN
        ! Get the solution from stress solver 
        IF(dofs == 1) THEN
          CALL GetLocalSolution( local_sol,UElement=Element,UVariable=AdjSol)
        ELSE
          CALL GetLocalSolution( local_sol_array,UElement=Element,UVariable=AdjSol)
          local_sol(1:m) = RESHAPE( local_sol_array(1:dofs,1:n), [m] )
        END IF
        ! This sign leads to convergence of the bisection iteration. 
        IF(sgn>0) THEN
          spos = spos + SUM( local_sol(1:m) * local_act(1:m) )
          ce(l) = ce(l) + sgn1 * SUM( local_sol(1:m) * local_act(1:m) )
        ELSE
          sneg = sneg + SUM( local_sol(1:m) * local_act(1:m) )
          ce(l) = ce(l) + sgn2 * SUM( local_sol(1:m) * local_act(1:m) )
        END IF          
      ELSE IF(NoModes > 0 ) THEN
        ce(l) = 0.0_dp
        DO k=1,NoModes
          IF(dofs == 1) THEN
            CALL GetLocalConsmode( local_sol,UElement=Element,USolver=PhysSolver,NoMode=k) 
          ELSE
            CALL GetLocalConsmode( local_sol_array,UElement=Element,USolver=PhysSolver,NoMode=k) 
            local_sol(1:m) = RESHAPE( local_sol_array(1:dofs,1:n), [m] )
          END IF
          ce(l) = ce(l) + SUM( local_sol(1:m) * local_act(1:m) )           
        END DO
      ELSE
        ce(l) = SUM( local_sol(1:m) * local_act(1:m) )
      END IF
    END DO
    
    ! Derivative of elemental energy
    dc = dc - penal*x**(penal-1) * efrac * ce

    ! Objective function
    obj = obj + SUM( (emin + efrac*( x**penal ) ) * ce )

    IF(InfoActive(20)) THEN
      PRINT *,'Objective:',obj,spos,sneg
      PRINT *,'ce:',SUM(ce),SUM(ABS(ce)),MINVAL(ce), MAXVAL(ce)
      PRINT *,'dc:',SUM(dc),SUM(ABS(dc)),MINVAL(dc), MAXVAL(dc)
      PRINT *,'x:',SUM(x),SUM(ABS(x)),MINVAL(x), MAXVAL(x)
    END IF
      
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
    REAL(KIND=dp) :: Vi, l1, l2, lmid, move, err, tol, V0, damp
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

    damp = ListGetCReal(Params,'Bisection search damping exponent',Found )
    IF(.NOT. Found) damp = 0.5_dp
      
    ! Desired total volume
    V0 = volFrac * SUM(dv)
    V0 = ParallelReduction(V0)
    
    ALLOCATE(xnew(SIZE(x)))
    xnew = 0.0_dp
        
    DO k=1,1000            
      lmid = 0.5*(l2+l1)
      
      ! Note: xnew in [0,1]
      ! Suggested new density
      xnew = x*(MAX(1.0e-10,-dc/(dv*lmid)))**damp

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
      IF( InfoActive(15)) THEN
        PRINT *,'Bisection:',k,Vi,l1,l2,err
      END IF

      IF( err < tol ) EXIT      
    END DO

    x = xnew 
    g = Vi - V0
    CALL Info(Caller,'Number of bisection iterations: '//I2S(k),Level=7)
    WRITE(Message,'(A,2ES12.3)') 'Volume target and accuracy: ',V0,g
    CALL Info(Caller, Message, Level=7)
    
    
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

    REAL(KIND=dp) :: xlow, xup, xmid, h, q
    REAL(KIND=dp), ALLOCATABLE :: histv(:), cumv(:), tmp_histv(:)

    REAL(KIND=dp) :: Vtot, Vtarget
    INTEGER :: i,j,k,m,iter,ierr
    LOGICAL :: Hit
    
    m = 100
    ALLOCATE(histv(0:m),cumv(0:m))
    IF(ParEnv % MyPe > 1 ) THEN
      ALLOCATE(tmp_histv(0:m))
    END IF
    
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

      IF( ParEnv % PEs > 1 ) THEN
        tmp_histv(0:m) = histv(0:m)
        CALL MPI_ALLREDUCE( tmp_histv, histv, m+1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF              

      cumv(0) = histv(0)
      DO i=1,m        
        cumv(i) = cumv(i-1) + histv(i)
      END DO

      Hit = .FALSE.
      q = 1.0_dp
      
      DO i=1,m
        IF(cumv(i-1) < Vtarget .AND. cumv(i) > Vtarget) THEN
          xlow = xlow + (i-1)*h
          xup = xlow + h

          q = (Vtarget-cumv(i-1))/(cumv(i)-cumv(i-1))

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
    
    xmid = (1-q)*xlow + q*xup

    WHERE(xPhys > xmid )
      bw = 1.0_dp
    ELSE WHERE
      bw = 0.0_dp
    END WHERE

    WRITE(Message,'(A,ES12.3)') 'Mass conserving B&W limit after '//I2S(iter)//' iters: ',xmid
    CALL Info(Caller,Message,Level=7)
    
  END SUBROUTINE DefineTopologyBW


!------------------------------------------------------------------------------
!> Given a nodal topology xPhys find a zero levelset such that the volume
!> constraint is conserved as accurately as possible.
!------------------------------------------------------------------------------
  SUBROUTINE DefineTopologyZeroLevel(xPhysVar,bw)    
    TYPE(Variable_t), POINTER :: xPhysVar
    REAL(KIND=dp), POINTER :: dv0(:), bw(:)

    REAL(KIND=dp) :: xlow, xup, xmid, h, xAtIp, weight, detJ, f, q
    REAL(KIND=dp), ALLOCATABLE :: histv(:), cumv(:), tmp_histv(:), Basis(:)

    REAL(KIND=dp) :: Vtot, Vtarget
    INTEGER :: i,j,k,m,n,t,iter,elem,RelOrder,ierr
    LOGICAL :: Hit, Stat
    TYPE(GaussIntegrationPoints_t) :: IP

    m = 1000
    ALLOCATE(histv(m+1),cumv(m+1))
    IF( ParEnv % PEs > 1 ) THEN
      ALLOCATE(tmp_histv(m+1))
    END IF
      
    xlow = 0.0_dp
    xup = 1.0_dp

    n = Mesh % MaxElementNodes
    ALLOCATE(Basis(n))
    Basis = 0.0_dp

    RelOrder = ListGetInteger( Solver % Values,'Levelset Integration Relative Order',Found)
    IF(.NOT. Found) RelOrder = 1
    
    DO iter=1,1 !0
      h = (xup-xlow) / m 
      histv = 0.0_dp
      cumv = 0.0_dp

      DO elem=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(elem)
        n = Element % Type % NumberOfNodes
        
        IP = GaussPoints(Element, RelOrder=RelOrder)

        DO t=1,IP % n
          stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis )
          xAtIp = SUM(Basis(1:n) * xPhysVar % Values(xPhysVar % Perm(Element % NodeIndexes)))
          weight = detJ * IP % s(t)

          IF(xAtIp <= xlow ) THEN
            histv(1) = histv(1) + weight 
          ELSE IF( xAtIp >= xup ) THEN
            histv(m+1) = histv(m+1) + weight 
          ELSE
            f = (xAtIp-xlow)/h
            j = CEILING(f)
            q = j-f
            histv(j+1) = histv(j+1) + (1-q) * weight 
            histv(j) = histv(j) + q * weight
          END IF
        END DO
      END DO

      IF( ParEnv % PEs > 1 ) THEN
        ! For parallel runs communicate the histogram here.
        tmp_histv = histv
        CALL MPI_ALLREDUCE( tmp_histv, histv, m+1, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF              

      cumv(1) = histv(1)
      DO i=2,m+1        
        cumv(i) = cumv(i-1) + histv(i)
      END DO

      Vtot = cumv(m+1) 
      Vtarget = (1-volFrac) * Vtot
            
      IF( ParEnv % MyPe == 0) THEN
        ! We may optionally save the histogram and its cumulative sum.  
        IF( ListGetLogical( Solver % Values,'Save Density Histogram',Found ) ) THEN
          BLOCK
            INTEGER :: IoUnit
            OPEN(NEWUNIT=IoUnit, FILE="dens_hist.dat")
            DO i=1,m+1
              WRITE(IoUnit,*) i-1, (i-1)*h, histv(i), cumv(i), cumv(i)/Vtot
            END DO
            CLOSE(IoUnit)
          END BLOCK
        END IF
      END IF
              
      Hit = .FALSE.
      q = 1.0_dp
      
      DO i=1,m
        IF(cumv(i) < Vtarget .AND. cumv(i+1) > Vtarget) THEN
          xlow = xlow + (i-1)*h
          xup = xlow + h

          q = (Vtarget-cumv(i))/(cumv(i+1)-cumv(i))
          Hit = .TRUE.
          EXIT
        ELSE IF(ABS(cumv(i)-Vtarget) < EPSILON(h)) THEN
          xlow = xlow + (i-1)*h
          xup = xlow
          EXIT
        ELSE IF(ABS(cumv(i+1)-Vtarget) < EPSILON(h)) THEN
          xlow = xlow + i*h
          xup = xlow          
          EXIT
        END IF
      END DO
      IF(.NOT. Hit) EXIT
    END DO

    ! This is the new approximation of the mid value that gives the desider volume within (x>xmid).
    xmid = (1-q)*xlow + q*xup
             
    IF( ListGetLogical( Solver % Values,'Levelset Symmmetric',Found ) ) THEN
      ! Define levelset as simple offset from nodal density.
      bw = xPhysVar % Values-xmid
    ELSE
      ! Map levelset between [-1,1] such that zero levelset is at desired value. 
      WHERE(xPhysVar % Values > xmid )
        bw = (xPhysVar % Values-xmid)/(1.0_dp-xmid)
      ELSE WHERE
        bw = (xPhysVar % Values-xmid)/xmid      
      END WHERE
    END IF
      
    WRITE(Message,'(A,ES12.3)') 'Mass conserving zero levelset: ',xmid
    CALL Info(Caller,Message,Level=7)
    
  END SUBROUTINE DefineTopologyZeroLevel

  

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
            STIFF(p,q) = STIFF(p,q) + Weight * Basis(p) * Basis(q) 
          END DO
        END DO

        DO p=1,n
          IF(SkipInterface) THEN
            IF(InterfaceNode(Element % NodeIndexes(p))) CYCLE
          END IF
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + Weight * &  
                SUM( Diff(1:dim) * dBasisdx(p,1:dim) * dBasisdx(q,1:dim) )
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
    
    CALL DefaultDirichletBCs()
    
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


#if 0 
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
#endif
  
!------------------------------------------------------------------------------
END SUBROUTINE TopoOpt
!------------------------------------------------------------------------------
