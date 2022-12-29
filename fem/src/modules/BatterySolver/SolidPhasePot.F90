!------------------------------------------------------------------------------
! For copyrights see the BatteryUtils.F90 file in this directory!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initial subroutine for the conservation of charge in solid phase equation
!------------------------------------------------------------------------------
SUBROUTINE SolidPhasePot_Init( Model,Solver,dt,Transient)
  !------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'SolidPhasePot_init'
  TYPE(ValueList_t), POINTER :: Params  
  INTEGER :: dim
  LOGICAL :: Found
  
  Params => GetSolverParams()
  dim = CoordinateSystemDimension()

  CALL ListAddNewString( Params,'Variable','PhiS')

  ! Let this solver allocate stuff for the solver that uses 1D mesh to compute Cs
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Cs' )   
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Cs Ave' )    
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'SOC')

  ! Overpotential
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Eta' )   
    
  ! The flux between solid and electrolyte phase
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Jli' )   

  IF( ListGetLogicalAnySolver( Model,'Save Solid Phase Diff') ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Cs Diff' )    
  END IF  

  IF( ListGetLogicalAnySolver( Model,'Study Jli Balance') ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Jli integral' )    
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Cs init' )    
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'Cs err' )    
  END IF  

  IF( ListGetLogical( Params,'Linearize Flux',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Calculate Phis Sensitivity',.TRUE.)
  END IF

!---------------------------------
END SUBROUTINE SolidPhasePot_Init
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves for the conservation of charge in solid phase equation
!------------------------------------------------------------------------------
SUBROUTINE SolidPhasePot(Model,Solver,dt,Transient)
  !------------------------------------------------------------------------------
  USE DefUtils
  USE BatteryModule

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
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active, dim
  INTEGER :: iter, maxiter
  LOGICAL :: Found, InitHandles, Newton
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(*), PARAMETER :: Caller = 'SolidPhasePot'   
  TYPE(Variable_t), POINTER :: SensVar
  
  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving potential of the solid phase')
  CALL Info(Caller,'------------------------------------------------')

  CALL InitializeBattery()

  CALL DefaultStart()

  Mesh => GetMesh()
  Params => GetSolverParams()

  dim = CoordinateSystemDimension() 

  maxiter = ListGetInteger( Params, &
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  Newton = ListGetLogical( Params,'Linearize Flux',Found )
  IF( Newton ) THEN
    CALL Info(Caller,'Using linearized flux',Level=8)
  END IF

  IF(.NOT. ASSOCIATED( PhisVar, Solver % Variable ) ) THEN
    CALL Fatal(Caller,'This solver should own "Phis"')
  END IF
  
    
  ! Nonlinear iteration loop/solver:
  !--------------------------
  DO iter=1,maxiter
    IF(maxiter>1) CALL Info(Caller,'Nonlinear system iteration: '//I2S(iter),Level=5)
    
    ! Calls new values of flux according to the Butler-Volmer equation    
    CALL ButlerVolmerUpdate(Solver)
    
    IF( Newton ) THEN
      SensVar => VariableGet( Mesh % Variables,'dJli dPhis')
      IF(.NOT. ASSOCIATED( SensVar ) ) THEN
        CALL Fatal(Caller,'Variable "dJli dPhis" not present!')
      END IF
    END IF

    ! System assembly:
    !----------------
    CALL DefaultInitialize() 

    CALL Info(Caller,'Performing bulk element assembly',Level=12)
    Active = GetNOFActive(Solver) 
    InitHandles = .TRUE.
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element) 
      nd = GetElementNOFDOFs(Element)  
      nb = GetElementNOFBDOFs(Element) 
      CALL LocalMatrix(  Element, n, nd+nb, nb, InitHandles )
    END DO
    CALL DefaultFinishBulkAssembly()
    
    CALL Info(Caller,'Performing boundary element assembly',Level=12)
    Active = GetNOFBoundaryActive(Solver) 
    InitHandles = .TRUE. 
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement(Element)) THEN 
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        nb = GetElementNOFBDOFs(Element)
        CALL LocalMatrixBC(  Element, n, nd+nb, nb, InitHandles )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly() 
    CALL DefaultFinishAssembly() 
    CALL DefaultDirichletBCs() 

    ! Solves the system of equations:
    !--------------------
    Norm = DefaultSolve()
    
    CALL VariableRange( PhisVar, 8, AnodeWeight ) 

    IF( Solver % Variable % NonlinConverged == 1 ) EXIT
  END DO

  CALL DefaultFinish() 
  
  CALL Info(Caller,'All done',Level=10)
 
 
CONTAINS

  !------------------------------------------------------------------------------
  ! Assembly of the matrix for solid phase potential in bulk.
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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: ElemSource(:), ElemSens(:), ElemPhis(:)
    REAL(KIND=dp) :: weight, detJ, SourceAtIp, CondAtIp, SensAtIp, PhisAtIp, EpsAtIp
    LOGICAL :: Stat,Found
    INTEGER :: i,j,p,q,t,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: Eps_h, CondCoeff_h
    !------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( CondCoeff_h,'Material',&
          'Solid Phase Electrical Conductivity')
      CALL ListInitElementKeyword( Eps_h,'Material',&
          'Active Particle Volume Fraction')
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()
    IP = GaussPointsAdapt( Element ) 

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(Basis(m), dBasisdx(m,3), MASS(m,m), STIFF(m,m), &
          FORCE(m), ElemSource(m), ElemSens(m), ElemPhis(m), &
          STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    ! Obtaining the source term j_li
    !-----------------------------------------    
    ElemSource(1:n) = JliVar % Values( JliVar % Perm( Element % NodeIndexes ) )
    IF( Newton ) THEN
      ElemSens(1:n) = SensVar % Values( SensVar % Perm( Element % NodeIndexes ) )
      ElemPhis(1:n) = PhisVar % Values( PhisVar % Perm( Element % NodeIndexes ) )
    END IF
    
    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp

    ! Generate weak form for Eq. (3.9) in [1]
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )

      ! Integration weight
      Weight = IP % s(t) * DetJ

      ! Conductivity term at IP + assembly of the K matrix
      CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element ) 
      EpsAtIp = ListGetElementReal( Eps_h, Basis, Element ) 
      
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          EpsAtIp * CondAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )

      ! Source term at IP + assembly of the right hand side matrix
      SourceAtIP = SUM( Basis(1:nd) * ElemSource(1:nd) )       
      
      ! Note that we give negative sign for the source as consistent with Eq. (3.9) in [1]
      FORCE(1:nd) = FORCE(1:nd) - Weight * SourceAtIP * Basis(1:nd)      
      
      IF( Newton ) THEN
        SensAtIp = SUM( Basis(1:n) * ElemSens(1:n) )
        PhisAtIp = SUM( Basis(1:n) * ElemPhis(1:n) ) 
        
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * SensAtIp * Basis(p) * Basis(q)
          END DO
          FORCE(p) = FORCE(p) + Weight * SensAtIp * PhisAtIp * Basis(p)  
        END DO
      END IF

    END DO

    CALL CondensateP( nd-nb, nb, STIFF, FORCE )

    ! Inserts local matrices into global matrix
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element )

  END SUBROUTINE LocalMatrix
  !---------------------------


  ! Assembly of the matrix entries arising from the boundary conditions.
  !------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, InitHandles )
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: Weight, DetJ, CurrDensAtIp, ExtPotAtIp, ExtCondAtIp
    REAL(KIND=dp) :: Basis(nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: CurrDens_h, ExtPot_h, ExtCond_h

    SAVE Nodes
    !------------------------------------------------------------------------------
    BC => GetBC(Element) 
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    ! Receives boundary values from the sif file
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( CurrDens_h,'Boundary Condition','Current Density') 
      CALL ListInitElementKeyword( ExtPot_h,'Boundary Condition','External Potential') 
      CALL ListInitElementKeyword( ExtCond_h,'Boundary Condition','External Conductivity') 
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ! queries the node points for the Element and reset the local matrices and vectors
    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    ! Numerical integration:
    !-----------------------
    IP = GaussPointsAdapt( Element )
    
    !integration loop:
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis )

      !Integration weights
      Weight = IP % s(t) * DetJ
      
      ! Given current density at current collectors:
      CurrDensAtIp = ListGetElementReal( CurrDens_h, Basis, Element, Found )
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * CurrDensAtIp * Basis(1:nd)
      END IF
 
      ExtCondAtIp = ListGetElementReal( ExtCond_h, Basis, Element, Found )
      IF( Found ) THEN
        ExtPotAtIp = ListGetElementReal( ExtPot_h, Basis, Element, Found )
        IF( Found ) THEN
          FORCE(1:nd) = FORCE(1:nd) + Weight * ExtCondAtIp * ExtPotAtIp * Basis(1:nd)
        END IF
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * ExtCondAtip * Basis(q) * Basis(p)
          END DO
        END DO
      END IF
    END DO

    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element)

  END SUBROUTINE LocalMatrixBC
  !------------------------------------------------------------------------------

END SUBROUTINE SolidPhasePot
!------------------------------------------------------------------------------
