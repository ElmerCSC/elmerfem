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
SUBROUTINE HeatSolver_Init0(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------  
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, Serendipity
  
  Params => GetSolverParams()

  IF( ListCheckPresentAnyEquation( Model,'Convection' ) .OR. &
      ListCheckPresentAnyEquation( Model,'Draw Velocity') .OR. &
      ListGetLogical( Params,'Bubbles',Found) ) THEN
    IF( .NOT. ListCheckPresent( Params,'Element') ) THEN
      Serendipity = GetLogical( GetSimulation(), 'Serendipity P Elements', Found)
      IF(.NOT.Found) Serendipity = .TRUE.
      IF(Serendipity) THEN
        CALL ListAddString(Params,'Element', &
            'p:1 -tri b:1 -tetra b:1 -quad b:3 -brick b:4 -prism b:4 -pyramid b:4')
      ELSE
        CALL ListAddString(Params,'Element', &
            'p:1 -tri b:1 -tetra b:1 -quad b:4 -brick b:8 -prism b:4 -pyramid b:4')
      END IF
      CALL ListAddNewLogical(Params,'Bubbles in Global System',.FALSE.)
    END IF
  END IF
  
!------------------------------------------------------------------------------
END SUBROUTINE HeatSolver_Init0
!------------------------------------------------------------------------------


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
  INTEGER :: dim
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
  CALL ListWarnUnsupportedKeyword('solver','Current Control',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('boundary condition','Phase Change',FatalFound=.TRUE.)

  IF(.NOT. ( DG .OR. DB ) ) THEN
    CALL ListWarnUnsupportedKeyword('boundary condition','Heat Gap',Found)
    IF( Found ) THEN
      CALL Fatal(Caller,'Keyword supported only with DG active: "Heat Gap"')
    END IF
  END IF

  ! These use one flag to call library features to compute automatically
  ! a conductivity matrix.
  IF( ListGetLogical(Params,'Calculate Conductance Matrix',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Matrix Symmetric',.TRUE.)
    CALL ListAddNewString( Params,'Constraint Modes Matrix Filename',&
        'ThermalConductanceMatrix.dat',.FALSE.)
    CALL ListRenameAllBC( Model,'Conductivity Body','Constraint Mode Temperature')
  END IF

  ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVIVTY
  CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
  
END SUBROUTINE HeatSolver_Init


!-----------------------------------------------------------------------------
!> A modern version for the heat equation supporting multi-threading and
!> SIMD friendly ElmerSolver kernels. This tries to be backward compatible
!> with the legacy HeatSolver but some rarely used features are missing. 
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
      DG, DB, Newton, HaveFactors, DiffuseGray, Radiosity, Spectral, &
      Converged, PostCalc = .FALSE.
  TYPE(Variable_t), POINTER :: PostWeight, PostFlux, PostAbs, PostEmis, PostTemp
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), POINTER :: Temperature(:)
  INTEGER, POINTER :: TempPerm(:)
  REAL(KIND=dp), ALLOCATABLE :: Temps4(:), Emiss(:), Absorp(:), Reflect(:),RadiatorPowers(:)
  REAL(KIND=dp) :: Norm, StefBoltz
  CHARACTER(LEN=MAX_NAME_LEN) :: EqName
  CHARACTER(*), PARAMETER :: Caller = 'HeatSolver'

  INTERFACE
    FUNCTION HeatSolver_Boundary_Residual( Model,Edge,Mesh,Quant,Perm,Gnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END FUNCTION HeatSolver_Boundary_Residual

    FUNCTION HeatSolver_Edge_Residual( Model,Edge,Mesh,Quant,Perm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END FUNCTION HeatSolver_Edge_Residual

    FUNCTION HeatSolver_Inside_Residual( Model,Element,Mesh,Quant,Perm, Fnorm ) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END FUNCTION HeatSolver_Inside_Residual
  END INTERFACE
  
  IF (.NOT. ASSOCIATED(Solver % Matrix)) RETURN

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving energy equation for temperature')

  ! The View and Gebhart factors may change if the shape and/or emissivities
  ! have changed. The routine may also affect matrix topology.
  !---------------------------------------------------------------------------
  Mesh => GetMesh()
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian ) 
  dim = CoordinateSystemDimension() 
  Params => GetSolverParams()  
  EqName = ListGetString( Params,'Equation', Found ) 

  Radiosity = GetLogical( Params, 'Radiosity Model', Found )
  Spectral = GetLogical( Params,'Spectral Model',Found )
  IF( Spectral ) Radiosity = .TRUE. 
  
  IF(.NOT.Radiosity) CALL RadiationFactors( Solver, .FALSE.,.FALSE.) 

  HaveFactors = ListCheckPresentAnyBC( Model,'Radiation')

  IF( HaveFactors ) THEN
    StefBoltz = ListGetConstReal( Model % Constants,&
        'Stefan Boltzmann',UnfoundFatal=HaveFactors)
  END IF

  Temperature => Solver % Variable % Values
  TempPerm => Solver % Variable % Perm
   
  DB = GetLogical( Params,'DG Reduced Basis',Found ) 
  DG = GetLogical( Params,'Discontinuous Galerkin',Found ) 

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
    CALL Info(Caller,'Heat solver iteration: '//I2S(iter))

    Newton = GetNewtonActive()

100 CONTINUE
    IF(Radiosity) CALL RadiationFactors( Solver, .FALSE., Newton) 
    
    ! Initialize the matrix equation to zero.
    !---------------------------------------
    CALL DefaultInitialize()
    CALL CalculateRadiosityFields(Pre=.TRUE.)
    
    ! For speed compute averaged emissivity and temperature over boundary elements
    ! for diffuse gray radiation.
    !-----------------------------------------------------------------------------
    IF( HaveFactors ) THEN
      CALL TabulateBoundaryAverages(Mesh, Temps4, Emiss, Absorp, Reflect) 
    END IF
    
    totelem = 0
    
    !$OMP PARALLEL &
    !$OMP SHARED(Solver, Active, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb,col, InitHandles) &
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    InitHandles = .TRUE.
    
    DO col=1,nColours
      
      !$OMP SINGLE
      CALL Info( Caller,'Assembly of colour: '//I2S(col),Level=15)
      Active = GetNOFActive(Solver)
      !$OMP END SINGLE
      
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

    BLOCK
      REAL(KIND=dp), POINTER :: RadiatorCoords(:,:)
      TYPE(ValueList_t), POINTER :: RadList

      ! If radiator is in body force section then use it:
      ! This will make it easier to make GUIs etc.
      IF( .NOT. ListCheckPresentAnyBodyForce( Model,'Radiator Coordinates',RadList ) ) &
          RadList => Params

      CALL GetConstRealArray( RadList, RadiatorCoords, 'Radiator Coordinates', Found)

      IF(Found) THEN
        n = SIZE(RadiatorCoords,1)
        ALLOCATE( RadiatorPowers(n))
        DO t=1,n
          RadiatorPowers(t)=GetCReal(RadList, 'Radiator Power '//I2S(t), Found)
        END DO
      END IF
    END BLOCK

    !!OMP PARALLEL &
    !!OMP SHARED(Active, Solver, nColours, VecAsm, DiffuseGray, RadiatorPowers ) &
    !!OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
    !!OMP REDUCTION(+:totelem) DEFAULT(NONE)
    InitHandles = .TRUE. 
    DO col=1,nColours
      !!OMP SINGLE
      CALL Info(Caller,'Assembly of boundary colour: '//I2S(col),Level=10)
      Active = GetNOFBoundaryActive(Solver)
      !!OMP END SINGLE      
      !!OMP DO
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
      !!OMP END DO
    END DO
    !!OMP END PARALLEL
    
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
    
    IF (ALLOCATED(RadiatorPowers)) DEALLOCATE( RadiatorPowers)
        
    CALL DefaultFinishBoundaryAssembly()
        
    CALL DefaultFinishAssembly()

    CALL DefaultDirichletBCs()

    ! Check stepsize for nonlinear iteration
    !------------------------------------------------------------------------------
    IF( DefaultLinesearch( Converged ) ) GOTO 100
    IF( Converged ) EXIT
        
    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( DefaultConverged(Solver) ) EXIT
  END DO
  
  CALL DefaultFinish()
  CALL CalculateRadiosityFields(Pre=.FALSE.)

 IF ( ListGetLogical( Solver % Values, 'Adaptive Mesh Refinement', Found ) ) THEN
   IF( .NOT. ListGetLogical(Params,'Library Adaptivity',Found) ) THEN
     CALL RefineMesh( Model,Solver,Temperature,TempPerm, &
         HeatSolver_Inside_Residual, HeatSolver_Edge_Residual, &
         HeatSolver_Boundary_Residual )
   END IF
 END IF
   
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
        SourceAtIpVec(:), RhoAtIpVec(:),VeloAtIpVec(:,:),ConvVelo(:,:),ConvVelo_i(:)

    LOGICAL :: Stat,Found,ConvComp,ConvConst
    INTEGER :: i,ngp,allocstat
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
        ConvVelo_h(3), PerfRate_h, PerfDens_h, PerfCp_h, &
        PerfRefTemp_h, VolSource_h, OrigMesh_h
    TYPE(VariableHandle_t), SAVE :: ConvField_h
    
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, DetJVec, &
    !$OMP               MASS, STIFF, FORCE, Nodes, ConvVelo, VeloAtIpVec, &
    !$OMP               ConvVelo_i, RhoAtIpVec, SourceAtIpVec, TmpVec, CPAtIpVec, CondAtIpVec, &
    !$OMP               Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
    !$OMP               ConvVelo_h, PerfRate_h, PerfDens_h, PerfCp_h, &
    !$OMP               PerfRefTemp_h, ConvField_h, VolSource_h, OrigMesh_h)
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJVec
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Source_h,'Body Force','Heat Source')
      CALL ListInitElementKeyword( VolSource_h,'Body Force','Volumetric Heat Source')
      CALL ListInitElementKeyword( Cond_h,'Material','Heat Conductivity')
      CALL ListInitElementKeyword( Cp_h,'Material','Heat Capacity')
      CALL ListInitElementKeyword( Rho_h,'Material','Density')

      IF( ListCheckPresentAnyMaterial( Model,'Draw Velocity' ) ) THEN
        CALL Fatal(Caller,'Vectorized assembly not implemented for "Draw Velocity"')
      END IF

      CALL ListInitElementKeyword( ConvFlag_h,'Equation','Convection')      
      DO i=1,3
        CALL ListInitElementKeyword( ConvVelo_h(i),'Material','Convection Velocity '//I2S(i))
      END DO

      str = GetString( Params, 'Temperature Convection Field', Found )
      IF(.NOT. Found ) str = 'Flow Solution'
      CALL ListInitElementVariable( ConvField_h, str )
      
      CALL ListInitElementKeyword( PerfRate_h,'Body Force','Perfusion Rate')
      CALL ListInitElementKeyword( PerfDens_h,'Body Force','Perfusion Density')
      CALL ListInitElementKeyword( PerfRefTemp_h,'Body Force','Perfusion Reference Temperature')
      CALL ListInitElementKeyword( PerfCp_h,'Body Force','Perfusion Heat Capacity')

      CALL ListInitElementKeyword( OrigMesh_h,'Equation','Convection Original Mesh')
      
      InitHandles = .FALSE.
    END IF

    IF( UseLocalMatrixCopy( Solver, Element % ElementIndex ) ) GOTO 10
    
    IP = GaussPointsAdapt(Element)
    ngp = IP % n

    !-----------------------------------------------------------------------------
    ! Output the number of integration points as information.
    ! This in not fully informative if several element types are present.
    !-----------------------------------------------------------------------------    
    IF( Element % ElementIndex == 1 ) THEN
      CALL Info(Caller,'Number of 1st integration points: '//I2S(IP % n), Level=10)
    END IF
        
    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
      IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
          DEALLOCATE(Basis, dBasisdx, DetJVec, MASS, STIFF, FORCE, &
          TmpVec, ConvVelo )
    END IF
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJVec(ngp), &
          MASS(nd,nd), STIFF(nd,nd), FORCE(nd), ConvVelo(ngp,3), &
          TmpVec(ngp), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    IF( ListGetElementLogical( OrigMesh_h ) ) THEN      
      CALL GetElementNodesOrigVec( Nodes, UElement=Element )
    ELSE
      CALL GetElementNodesVec( Nodes, UElement=Element )
    END IF
      
    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp

    ConvConst = ListCompareElementString( ConvFlag_h,'constant',Element, Found )    
    ConvComp = ListCompareElementString( ConvFlag_h,'computed',Element, Found )
    
    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJvec, &
        SIZE(Basis,2), Basis, dBasisdx )
    
    ! Compute actual integration weights (recycle the memory space of DetJVec)
    DetJVec(1:ngp) = IP % s(1:ngp) * DetJVec(1:ngp)

    ! Get pointer to vector including density on all integration points
    RhoAtIpVec => ListGetElementRealVec( Rho_h, ngp, Basis, Element, Found ) 

    ! thermal conductivity term: STIFF=STIFF+(kappa*grad(u),grad(v))
    CondAtIpVec => ListGetElementRealVec( Cond_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, dim, dBasisdx, DetJVec, STIFF, CondAtIpVec )
    END IF

    ! We need heat capacity only if the case is transient or we have convection
    IF( ConvConst .OR. ConvComp .OR. Transient ) THEN
      CpAtIpVec => ListGetElementRealVec( Cp_h, ngp, Basis, Element, Found ) 
      TmpVec(1:ngp) = CpAtIpVec(1:ngp) * RhoAtIpVec(1:ngp)        
    END IF

    ! convection, either constant or computed
    ! STIFF=STIFF+(C*grad(u),v)
    IF( ConvConst .OR. ConvComp ) THEN
      IF( ConvConst ) THEN
        DO i=1,dim
          ConvVelo_i => ListGetElementRealVec( ConvVelo_h(i), ngp, Basis, Element, Found ) 
          IF( Found ) ConvVelo(1:ngp,i) = ConvVelo_i(1:ngp)
        END DO
        VeloAtIpVec => ConvVelo
      ELSE
        VeloAtIpVec => ListGetElementVectorSolutionVec( ConvField_h, ngp, dim, Basis, Element )
      END IF      
      CALL LinearForms_GradUdotU(ngp, nd, dim, dBasisdx, Basis, DetJVec, STIFF, &
          TmpVec, VeloAtIpVec )       
    END IF
            
    ! time derivative term: MASS=MASS+(rho*cp*dT/dt,v)
    IF( Transient ) THEN
      CALL LinearForms_UdotU(ngp, nd, dim, Basis, DetJVec, MASS, TmpVec )
    END IF
      
    ! source term: FORCE=FORCE+(u,f)
    SourceAtIpVec => ListGetElementRealVec( VolSource_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJVec, SourceAtIpVec, FORCE)      
    ELSE
      SourceAtIpVec => ListGetElementRealVec( Source_h, ngp, Basis, Element, Found ) 
      IF( Found ) THEN
        TmpVec(1:ngp) = SourceAtIpVec(1:ngp) * RhoAtIpVec(1:ngp)        
        CALL LinearForms_UdotF(ngp, nd, Basis, DetJVec, TmpVec, FORCE)
      END IF
    END IF
      
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
10  CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------


  ! We have a special type of velocity implemented that follows thin regions
  ! assuming that convection is aligned with the elements, and the elements
  ! are structural ones. Either 404 or 808 type of elements are ok as for now.
  !------------------------------------------------------------------------------
  FUNCTION CalculatePlateTangent(n,Nodes) RESULT ( PlateTan )
    INTEGER :: n    
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: PlateTan(3)

    INTEGER, SAVE :: ActiveCoord = -1
    INTEGER :: i, sgn
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp) :: xmean
    
    IF( ActiveCoord < 1 ) THEN
      ActiveCoord = ListGetInteger( Params,'Active Coordinate',Found )
      IF(.NOT. Found ) THEN
        CALL Fatal('CalculatePlateTangent','Keyword "Draw Velocity" requires "Active Coordinate" to be given!')
      END IF
    END IF
      
    IF(ActiveCoord==1) THEN
      x => Nodes % x
    ELSE IF(ActiveCoord==2) THEN
      x => Nodes % y
    ELSE IF(ActiveCoord==3) THEN
      x => Nodes % z
    ELSE
      CALL Fatal('CalculatePlateTangent','"Active Coordinate" must be either 1, 2 or 3!')
    END IF

    IF( n /= 4 .AND. n /= 8 ) THEN
      CALL Warn('CalculatePlateTangent',&
          'Heuristics is well suited only for structural meshes: '//I2S(n))
    END IF
    
    xmean = SUM(x(1:n)) / n

    PlateTan = 0.0_dp
    DO i=1,n
      IF(x(i) > xmean ) THEN
        sgn = 1
      ELSE
        sgn = -1
      END IF
      PlateTan(1) = PlateTan(1) + sgn * Nodes % x(i)
      PlateTan(2) = PlateTan(2) + sgn * Nodes % y(i)
      PlateTan(3) = PlateTan(3) + sgn * Nodes % z(i)
    END DO

    PlateTan = PlateTan / SQRT( SUM( PlateTan**2 ) )
    
  END FUNCTION CalculatePlateTangent


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
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,nd)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: weight, SourceAtIp, CpAtIp, RhoAtIp, CondAtIp, DetJ, A, VeloAtIp(3)
    REAL(KIND=dp) :: PerfRateAtIp, PerfDensAtIp, PerfCpAtIp, PerfRefTempAtIp, PerfCoeff
    REAL(KIND=dp) :: PlateTangent(3), PlateSpeed
    REAL(KIND=dp), POINTER :: CondTensor(:,:)
    LOGICAL :: Stat,Found,ConvComp,ConvConst
    INTEGER :: i,j,t,p,q,CondRank
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    TYPE(ValueHandle_t), SAVE :: Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, &
        ConvVelo_h, PlateSpeed_h, PerfRate_h, PerfDens_h, PerfCp_h, PerfRefTemp_h, &
        VolSource_h, OrigMesh_h

    TYPE(VariableHandle_t), SAVE :: ConvField_h

!$OMP  THREADPRIVATE(Source_h, Cond_h, Cp_h, Rho_h, ConvFlag_h, ConvVelo_h, PlateSpeed_h, PerfRate_h, &
!$OMP& PerfDens_h, PerfCp_h, PerfRefTemp_h, VolSource_h, OrigMesh_h, ConvField_h, Nodes )

!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Source_h,'Body Force','Heat Source')
      CALL ListInitElementKeyword( VolSource_h,'Body Force','Volumetric Heat Source')
      CALL ListInitElementKeyword( Cond_h,'Material','Heat Conductivity')
      CALL ListInitElementKeyword( Cp_h,'Material','Heat Capacity')
      CALL ListInitElementKeyword( Rho_h,'Material','Density')

      CALL ListInitElementKeyword( ConvFlag_h,'Equation','Convection')
      
      CALL ListInitElementKeyword( ConvVelo_h,'Material','Convection Velocity',InitVec3D=.TRUE.)
      CALL ListInitElementKeyword( PlateSpeed_h,'Material','Draw Velocity')

      str = GetString( Params, 'Temperature Convection Field', Found )
      IF(.NOT. Found ) str = 'Flow Solution'
      CALL ListInitElementVariable( ConvField_h, str )
      
      CALL ListInitElementKeyword( PerfRate_h,'Body Force','Perfusion Rate')
      CALL ListInitElementKeyword( PerfDens_h,'Body Force','Perfusion Density')
      CALL ListInitElementKeyword( PerfRefTemp_h,'Body Force','Perfusion Reference Temperature')
      CALL ListInitElementKeyword( PerfCp_h,'Body Force','Perfusion Heat Capacity')

      CALL ListInitElementKeyword( OrigMesh_h,'Equation','Convection Original Mesh')
      
      InitHandles = .FALSE.
    END IF

    IF( UseLocalMatrixCopy( Solver, Element % ElementIndex ) ) GOTO 20
    
    IP = GaussPointsAdapt( Element )
    IF( Element % ElementIndex == 1 ) THEN
      CALL Info(Caller,'Number of 1st integration points: '//I2S(IP % n), Level=10)
    END IF
      
    IF( ListGetElementLogical( OrigMesh_h ) ) THEN
      CALL GetElementNodesOrig( Nodes, UElement=Element )
    ELSE
      CALL GetElementNodes( Nodes, UElement=Element )
    END IF
      
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
          PlateSpeed = ListGetElementReal( PlateSpeed_h, Basis, Element, Found )
          IF( Found ) THEN
            IF(t==1) PlateTangent = CalculatePlateTangent(n,Nodes)
            VeloAtIp = PlateTangent * PlateSpeed
          ELSE
            VeloAtIp = ListGetElementReal3D( ConvVelo_h, Basis, Element )
          END IF
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

      SourceAtIP = ListGetElementReal( VolSource_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      ELSE
        SourceAtIP = ListGetElementReal( Source_h, Basis, Element, Found ) 
        IF( Found ) THEN
          FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * RhoAtIp * Basis(1:nd)
        END IF
      END IF
    END DO
    
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
20  CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Compute the fraction to be assembled. In serial case it is always one,
! in parallel case only true parents result to assembly, mixed parents gives
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
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), ElemWeight(nd)
    LOGICAL :: Stat,Found,RobinBC,RadIdeal,RadDiffuse,TorBC
    INTEGER :: t,p,q,Indexes(n)
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       

    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: HeatFlux_h, HeatTrans_h, ExtTemp_h, Farfield_h, &
        RadFlag_h, RadExtTemp_h, EmisBC_h, EmisMat_h, TorBC_h 

    !$OMP  THREADPRIVATE(Nodes,HeatFlux_h,HeatTrans_h,ExtTemp_h,Farfield_h,RadFlag_h, &
    !$OMP& RadExtTemp_h, EmisBC_h, EmisMat_h, TorBC_h )
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
      CALL ListInitElementKeyword( TorBC_h,'Boundary Condition','Radiator BC')
      
      InitHandles = .FALSE.
    END IF

    ! In parallel if we have halo the same BC element may occur several times.
    ! Fetch here the fraction of the assembly to be accounted in this occurrence.
    AssFrac = BCAssemblyFraction(Element)
    IF( AssFrac < TINY( AssFrac ) ) RETURN

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    ElemWeight = 0._dp
    
    RadIdeal = ListCompareElementString( RadFlag_h,'idealized',Element, Found )    
    RadDiffuse = ListCompareElementString( RadFlag_h,'diffuse gray',Element, Found )

    IF( DG ) THEN
      CALL DgRadiationIndexes(Element,n,Indexes,.FALSE.)
    END IF
    
    ! This routine does not do diffuse gray radiation.
    ! Pass on the information to the routine that does. 
    DiffuseGray = RadDiffuse
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )

    ! Is this a radiator BC? 
    TorBC = ListGetElementLogical( TorBC_h, Element, Found = Found ) 
    TorBC = TorBC .AND. .NOT. DiffuseGray

        
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      Weight = IP % s(t) * DetJ
      
      IF ( AxiSymmetric ) THEN
        Weight = Weight * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF
      ElemWeight(1:nd) = ElemWeight(1:nd) + Weight * Basis(1:nd) 
      
      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------

      F = ListGetElementReal( HeatFlux_h, Basis, Element, Found )
      IF( TorBC ) THEN
        IF(ALLOCATED(Element % Boundaryinfo % Radiators)) THEN
          Found = .TRUE.
          F = F + SUM(RadiatorPowers*Element % BoundaryInfo % Radiators)
        END IF
      END IF

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

        Emis = ListGetElementRealParent( EmisMat_h, Basis, Element = Element, Found = Found )
        IF( .NOT. Found ) THEN
          Emis = ListGetElementReal( EmisBC_h, Basis, Element = Element, Found = Found ) 
        END IF
        IF(.NOT. Found ) THEN
          CALL Warn(Caller,'Emissivity should be available for radiating BC: '&
              //TRIM(ListGetString(BC,'name')))
          CYCLE
        END IF
        
        IF( DG ) THEN
          T0 = SUM( Basis(1:n) * Temperature(TempPerm(Indexes(1:n))))
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
      ELSE
        RadC = 0; RadF=0;
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

    ! Calculate fluxes on-the-fly
#if 0
    IF( PostCalc ) THEN
      BLOCK
        INTEGER :: ElemPerm(27)
        ElemPerm(1:n) = PostFlux % Perm(Element % NodeIndexes)
        IF(ALL(ElemPerm(1:n) > 0 )) THEN
          PostWeight % Values(ElemPerm(1:n)) = PostWeight % Values(ElemPerm(1:n)) + ElemWeight(1:n)
          PostFlux % Values(ElemPerm(1:n)) = PostFlux % Values(ElemPerm(1:n)) + FORCE(1:n)
        END IF
      END BLOCK
    END IF
#endif
    
    IF( ABS(AssFrac-1.0_dp) > TINY( AssFrac ) ) THEN
      FORCE(1:nd) = AssFrac * FORCE(1:nd)
      STIFF(1:nd,1:nd) = AssFrac * STIFF(1:nd,1:nd)
    END IF

    IF( DG ) THEN
      CALL UpdateGlobalEquations( Solver % Matrix, STIFF, &
          Solver % Matrix % Rhs, FORCE, n, 1, TempPerm(Indexes(1:n)), UElement=Element)      
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
  SUBROUTINE TabulateBoundaryAverages( Mesh, Temps4, Emiss, Absorp, Reflect )
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     REAL(KIND=dp), ALLOCATABLE :: Temps4(:)
     REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: Emiss(:), Absorp(:), Reflect(:)
 !------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
     INTEGER :: bindex, nb, n, j, noactive
     INTEGER :: ElemInds(12)
     REAL(KIND=dp) :: NodalVal(12), NodalTemp(12)
     
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
         ALLOCATE( Temps4(nb), Emiss(nb), Absorp(nb), Reflect(nb) )
         Temps4 = 0.0_dp
         Emiss = 0.0_dp
         Absorp = 0.0_dp
         Reflect = 0.0_dp
       END IF
         
       n = GetElementNOFNodes(Element)

       IF( DG ) THEN
         CALL DgRadiationIndexes(Element,n,ElemInds,.TRUE.)
         NodalTemp(1:n) = Temperature(TempPerm(ElemInds(1:n)))
       ELSE
         NodalTemp(1:n) = Temperature(TempPerm(Element % NodeIndexes))
       END IF
       Temps4(j) = ( SUM( NodalTemp(1:n)**4 )/ n )**(1._dp/4._dp)       

       IF( PRESENT( Emiss ) ) THEN
         NodalVal(1:n) = GetReal(BC,'Emissivity',Found)
         IF (Found) THEN
           Emiss(j) = SUM(NodalVal(1:n)) / n
           NodalVal(1:n) = GetReal(BC,'Absorptivity',Found)
           IF(Found) THEN
             Absorp(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp(j) = Emiss(j)
           END IF
           NodalVal(1:n) = GetReal(BC,'Reflectivity',Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp(j)
           END IF
         ELSE
           NodalVal(1:n) = GetParentMatProp('Emissivity',Element)
           Emiss(j) = SUM(NodalVal(1:n)) / n
           NodalVal(1:n) = GetParentMatProp('Absorptivity',Element, Found)
           IF(Found) THEN
             Absorp(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Absorp(j) = Emiss(j)
           END IF
           NodalVal(1:n) = GetParentMatProp('Reflectivity',Element, Found)
           IF(Found) THEN
             Reflect(j) = SUM(NodalVal(1:n)) / n
           ELSE
             Reflect(j) = 1-Absorp(j)
           END IF
         END IF
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
    REAL(KIND=dp) :: T0,Text, Fj, &
        RadLoadAtIp, AngleFraction, Topen, Emis1, Abso1, Refl1, AssFrac
    REAL(KIND=dp) :: Basis(nd),DetJ,Atext(12),Base(12),S,RadCoeffAtIP
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), TempAtIp
    REAL(KIND=dp), POINTER :: Fact(:) 
    TYPE(Element_t), POINTER :: RadElement
    LOGICAL :: Stat,Found,BCOpen,Radiators
    INTEGER :: j,t,p,q,bindex,k,k1,k2,nf,nf_imp
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    INTEGER, POINTER :: ElementList(:),pIndexes(:)
    REAL(KIND=dp), POINTER :: ForceVector(:)   
    REAL(KIND=dp) :: NodalTemp(12)
    INTEGER, TARGET :: ElemInds(12),ElemInds2(12)
    
    SAVE Nodes

!$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN
    
    IF( ListGetString( BC,'Radiation',Found) /= 'diffuse gray' ) RETURN

    AssFrac = BCAssemblyFraction(Element)
    IF( AssFrac < TINY( AssFrac ) ) RETURN

    CALL GetElementNodes( Nodes, UElement=Element) 
    n = Element % TYPE % NumberOfNodes
    
    IF( .NOT. ASSOCIATED( Element % BoundaryInfo % RadiationFactors ) ) THEN
      CALL Fatal(Caller,'Radiation factors not calculated for boundary!')
    END IF
    
    Fact => Element % BoundaryInfo % RadiationFactors % Factors
    ElementList => Element % BoundaryInfo % RadiationFactors % Elements

    bindex = Element % ElementIndex - Solver % Mesh % NumberOfBulkElements
    nf = Element % BoundaryInfo % RadiationFactors % NumberOfFactors
      
    nf_imp = Element % BoundaryInfo % RadiationFactors % NumberOfImplicitFactors      
    IF( nf_imp == 0 ) nf_imp = nf

    ForceVector => Solver % Matrix % rhs
    Temperature => Solver % Variable % Values
    TempPerm => Solver % Variable % Perm

    Emis1 = Emiss(bindex)
    Refl1 = Reflect(bindex)
    Abso1 = Absorp(bindex)
    
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
      CALL DgRadiationIndexes(Element,n,ElemInds,.TRUE.)
      NodalTemp(1:n) = Temperature( TempPerm( ElemInds(1:n) ) )
    ELSE
      NodalTemp(1:n) = Temperature( TempPerm( Element % NodeIndexes ) )
    END IF
    
    Text  = 0.0_dp

    Radiators = ALLOCATED(Element % BoundaryInfo % Radiators) .AND. &
             ALLOCATED(RadiatorPowers)
        
    IF(Radiosity) THEN 
      IF( BCOpen ) THEN
        CALL Fatal(Caller,'Radiosity model not yet working with open boundaries!')
      END IF

      IF(Refl1 < EPSILON(Refl1) ) THEN
        CALL Fatal(Caller,'Radiosity Model does not work for zero reflectivity (emissivity one)!')
      END IF

      Base = 0.0_dp
      IF(.NOT. Spectral) Emis1 = Emis1 / Refl1     
      TempAtIp = SUM(NodalTemp(1:n))/n

      IF(Newton) THEN
        RadLoadAtIp =  (3 * Emis1 * TempAtIp**3 * StefBoltz - Fact(2)) * TempAtIp &
             + Fact(1) 
        RadCoeffAtIp = 4 * Emis1 * TempAtIp**3 * StefBoltz - Fact(2)
      ELSE
        RadCoeffAtIp = Emis1 * StefBoltz * TempAtIp**3
        RadLoadAtIp = Fact(1)
      END IF
      
      DO t=1,IP % n
        stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
        s = detJ * IP % s(t)        
        IF ( AxiSymmetric ) THEN
          s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
        END IF

        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + s * Basis(p)*Basis(q) * RadCoeffAtIp 
          END DO
          FORCE(p) = FORCE(p) + s * Basis(p) * RadLoadAtIp 
        END DO
        Base(1:n) = Base(1:n) + s * Basis(1:n) 
      END DO
        
    ELSE ! .NOT. Radiosity ) 
      ! Go through surfaces (j) this surface (i) is getting radiated from.
      !------------------------------------------------------------------------------        
      IF ( Newton ) THEN                
        ! Linearization of T^4_i term
        !----------------------------------------------------------------------------
        Base = 0.0_dp

        DO t=1,IP % n
          stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
          s = detJ * IP % s(t)        
          IF ( AxiSymmetric ) THEN
            s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
          END IF

          TempAtIp = SUM(NodalTemp(1:n))/n
          RadCoeffAtIp = 4 * Emis1 * TempAtIp**3 * StefBoltz
          RadLoadAtIp =  3 * Emis1 * TempAtIp**4 * StefBoltz

          DO p=1,n
            DO q=1,n
              STIFF(p,q) = STIFF(p,q) + s * Basis(p)*Basis(q)*RadCoeffAtIp 
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

          ! Gebhart factors are given elementwise at the center
          ! of the element, so take average of nodal temperatures
          !-------------------------------------------------------------
          bindex = ElementList(j) - Solver % Mesh % NumberOfBulkElements
          Text = Temps4(bindex)

          IF( j <= nf_imp ) THEN        
            ! Linearization of the G_jiT^4_j term
            !------------------------------------------------------------------------------
            RadCoeffAtIp = -4 * Fj * Text**3 * StefBoltz
            RadLoadAtIp  = -3 * Fj * Text**4 * StefBoltz

            IF(Radiators) THEN
              IF(ALLOCATED(RadElement % BoundaryInfo % Radiators)) THEN
                RadLoadAtIp = RadLoadAtIp + &
                    Fj * SUM(RadElement % BoundaryInfo % Radiators * RadiatorPowers) * &
                    (1-Emiss(bindex)) / Emiss(bindex)
              END IF
            END IF

            ! Integrate the contribution of surface j over surface j and add to global matrix
            !------------------------------------------------------------------------------                    
            IF( Dg ) THEN
              CALL DgRadiationIndexes(RadElement,k,ElemInds2,.TRUE.)                              

              DO p=1,n
                k1 = TempPerm( ElemInds(p))
                DO q=1,k
                  k2 = TempPerm( ElemInds2(q) )
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

        ! Add radiators in case the radiosity model is not used
        !----------------------------------------------------------------------------
        IF( Radiators ) THEN
          DO p=1,n
            FORCE(p) = FORCE(p) + Base(p) * Emis1 * SUM(Element % BoundaryInfo % &
                Radiators * RadiatorPowers ) 
          END DO
        END IF
        
      ELSE ! .NOT. Newton 
        ! Compute the weighted sum of T^4
        
        Text = 0._dp
        DO j=1,nf
          Fj = Fact(j)

          RadElement => Mesh % Elements(ElementList(j))
          bindex = RadElement % ElementIndex - Solver % Mesh % NumberOfBulkElements
          Text = Text + Fj*Temps4(bindex)**4 / Emis1

          IF(Radiators) THEN
            IF(ALLOCATED(RadElement % BoundaryInfo % Radiators)) THEN
              Text = Text + Fj * &
                  SUM(RadElement % BoundaryInfo % Radiators * RadiatorPowers) * &
                  (1-Emiss(bindex)) / Emiss(bindex) / Emis1 / StefBoltz
            END IF
          END IF
        END DO
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
        Base = 0.0_dp
        Text = Text**0.25_dp
        DO t=1,IP % n
          stat = ElementInfo( Element,Nodes,IP % u(t),IP % v(t),IP % w(t),detJ,Basis )
          s = detJ * IP % s(t)        
          IF ( AxiSymmetric ) THEN
            s = s * SUM( Nodes % x(1:n) * Basis(1:n) )
          END IF

          T0 = SUM( Basis(1:n) * NodalTemp(1:n) )
          RadCoeffAtIp = Emis1 * StefBoltz*(T0**3 + T0**2*Text + T0*Text**2 + Text**3)

          DO p=1,n
            DO q=1,n
              STIFF(p,q) = STIFF(p,q) + s * Basis(p) * Basis(q) * RadCoeffAtIp
            END DO
            FORCE(p) = FORCE(p) + s * Basis(p) * RadCoeffAtIp * Text
          END DO
          Base(1:n) = Base(1:n) + s * Basis(1:n) 
        END DO
      END IF
    END IF ! .NOT. Radiosity
      
    ! Calculate fluxes on-the-fly
    IF( PostCalc ) THEN
      BLOCK
        INTEGER :: ElemPerm(27)
        ElemPerm(1:n) = PostFlux % Perm(Element % NodeIndexes)
        IF(ALL(ElemPerm(1:n) > 0 )) THEN
          PostWeight % Values(ElemPerm(1:n)) = PostWeight % Values(ElemPerm(1:n)) + Base(1:n)
          PostFlux % Values(ElemPerm(1:n)) = PostFlux % Values(ElemPerm(1:n)) + Fact(1) * Base(1:n)
          IF( Spectral ) THEN
            PostEmis % Values(ElemPerm(1:n)) = PostEmis % Values(ElemPerm(1:n)) + Emiss(bindex) * Base(1:n)
            PostAbs % Values(ElemPerm(1:n)) = PostAbs % Values(ElemPerm(1:n)) + Fact(3) * Base(1:n)
            PostTemp % Values(ElemPerm(1:n)) = PostTemp % Values(ElemPerm(1:n)) + Fact(4) * Base(1:n)
          END IF
        END IF
      END BLOCK
    END IF
    
    ! Glue standard local matrix equation to the global matrix
    ! The view factor part has already been glued.
    !-----------------------------------------------------------------
    
    IF( DG ) THEN
      pIndexes => ElemInds
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
    INTEGER :: i, k, p, q, m, allocstat
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
        
    SAVE LeftActive, HaveJump

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
    
    ! Switch the reference body always to Parent1
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
    REAL(KIND=dp) :: Basis(n), detJ, S, alpha, beta, AssFrac
    LOGICAL :: Stat
    INTEGER :: i, j, p, q, t, i1, i2, JumpOrder, ntmp
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
        FORCE(p) = FORCE(p) + beta/2 * Basis(p) * s 
        FORCE(p+n) = FORCE(p+n) + beta/2 * Basis(p) * s 
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

  SUBROUTINE CalculateRadiosityFields(Pre) 
    LOGICAL :: Pre    
    LOGICAL :: Visited = .FALSE., CalcRadiosityFields = .TRUE.
    INTEGER, POINTER :: Perm(:)
    INTEGER :: i,t,nsize = 0
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: c
    
    SAVE Perm, nsize, CalcRadiosityFields, Visited

    IF(.NOT. CalcRadiosityFields ) RETURN
    
    IF(.NOT. Visited ) THEN    
      Visited = .TRUE.
      CalcRadiosityFields = ListGetLogical( Params,'Calculate Radiosity Fields',Found ) 
      IF( CalcRadiosityFields ) THEN
        IF(.NOT. Radiosity ) THEN
          CALL Warn('CalculateRadiosityFields','Radiosity Model is not active, fields omitted!')
          CalcRadiosityFields = .FALSE.
          RETURN
        END IF
        
        ALLOCATE(Perm(Solver % Mesh % NumberOfNodes))
        Perm = 0
        
        CALL Info(Caller,'Creating permutation for radiosity fields',Level=8)
        DO t=Mesh % NumberOfBulkElements+1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(t)
          BC => GetBC(Element)
          IF(.NOT. ASSOCIATED( BC ) ) CYCLE       
          IF( ListCheckPresent( BC,'Radiation') .OR. &
              ListCheckPresent( BC,'Radiator') ) THEN
            Perm(Element % NodeIndexes) = 1
          END IF
        END DO
        
        nsize = 0
        DO i=1,Mesh % NumberOfNodes
          IF(Perm(i) > 0) THEN
            nsize = nsize+1
            Perm(i) = nsize
          END IF
        END DO
        CALL Info(Caller,'Number of active nodes for boundary fields: '//I2S(nsize),Level=10)
        
        CALL DefaultVariableAdd('Radiation Weight',Perm=Perm,Var=PostWeight,Output=.FALSE.)
        CALL DefaultVariableAdd('Radiation Flux',Perm=Perm,Var=PostFlux,Secondary=.TRUE.)
        IF( Radiosity ) THEN
          CALL DefaultVariableAdd('Absorptivity',Perm=Perm,Var=PostAbs,Secondary=.TRUE.)
          CALL DefaultVariableAdd('Emissivity',Perm=Perm,Var=PostEmis,Secondary=.TRUE.)
          CALL DefaultVariableAdd('Radiation Temperature',Perm=Perm,Var=PostTemp,Secondary=.TRUE.)
        END IF
        IF(nsize == 0) CalcRadiosityFields = .FALSE.
      END IF
      PostCalc = CalcRadiosityFields 
    ELSE IF(Pre) THEN
      PostWeight % Values = 0.0_dp
      PostFlux % Values = 0.0_dp
      IF(Spectral) THEN
        PostAbs % Values = 0.0_dp
        PostEmis % Values = 0.0_dp
        PostTemp % Values = 0.0_dp
      END IF
    ELSE
      WHERE( PostWeight % Values > EPSILON(c) )
        PostFlux % Values = PostFlux % Values / PostWeight % Values
      END WHERE
      IF( Spectral ) THEN
        WHERE( PostWeight % Values > EPSILON(c) )
          PostAbs % Values = PostAbs % Values / PostWeight % Values
          PostEmis % Values = PostEmis % Values / PostWeight % Values
          PostTemp % Values = PostTemp % Values / PostWeight % Values
        END WHERE
      END IF      
    END IF
         
  END SUBROUTINE CalculateRadiosityFields

!------------------------------------------------------------------------------
END SUBROUTINE HeatSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION HeatSolver_Boundary_Residual( Model, Edge, Mesh, Quant, Perm,Gnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     USE Radiation

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,l,t,dim,Pn,En,nd
     LOGICAL :: stat, Found
     INTEGER, ALLOCATABLE :: Indexes(:)
     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), ExtTemperature(:), &
         TransferCoeff(:), EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
         dBasisdx(:,:), Temperature(:), Flux(:), NodalEmissivity(:)
     REAL(KIND=dp) :: Conductivity, Emissivity, StefanBoltzmann
     REAL(KIND=dp) :: Normal(3), EdgeLength, gx, gy, gz
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: Residual, ResidualNorm
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     LOGICAL :: First = .TRUE., Dirichlet
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Indicator = 0.0d0
     Gnorm     = 0.0d0

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT
!    
!    ---------------------------------------------

     Element => Edge % BoundaryInfo % Left

     IF ( .NOT. ASSOCIATED( Element ) ) THEN
        Element => Edge % BoundaryInfo % Right
     ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
        Element => Edge % BoundaryInfo % Right
     END IF

     IF ( .NOT. ASSOCIATED( Element ) ) RETURN
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     en = Edge % TYPE % NumberOfNodes
     pn = Element % TYPE % NumberOfNodes

     ALLOCATE(EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     nd = GetElementNOFDOFs(Element)
     ALLOCATE( Temperature(nd), Basis(nd), ExtTemperature(nd), &
        TransferCoeff(en), x(en), y(en), z(en), EdgeBasis(nd), &
        dBasisdx(nd,3), NodalConductivity(nd), Flux(nd), &
        NodalEmissivity(nd), Indexes(nd) ) 

     nd = GetElementDOFs(Indexes,Element)

     ALLOCATE(Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )
     Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))


     DO l = 1,en
       DO k = 1,pn
          IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
             x(l) = Element % TYPE % NodeU(k)
             y(l) = Element % TYPE % NodeV(k)
             z(l) = Element % TYPE % NodeW(k)
             EXIT
          END IF
       END DO
     END DO
!
!    Integrate square of residual over boundary element:
!    ---------------------------------------------------

     Indicator    = 0.0d0
     EdgeLength   = 0.0d0
     ResidualNorm = 0.0d0

     DO j=1,Model % NumberOfBCs
        IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag ) CYCLE

!       IF ( .NOT. ListGetLogical( Model % BCs(j) % Values, &
!                 'Heat Flux BC', Found ) ) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
        s = ListGetConstReal( Model % BCs(j) % Values,'Temperature',Dirichlet )


!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
        Flux(1:en) = ListGetReal( Model % BCs(j) % Values, &
          'Heat Flux', en, Edge % NodeIndexes, Found )

!       ...convective heat transfer:
!       ----------------------------
        TransferCoeff(1:en) =  ListGetReal( Model % BCs(j) % Values, &
          'Heat Transfer Coefficient', en, Edge % NodeIndexes, Found )

        ExtTemperature(1:En) = ListGetReal( Model % BCs(j) % Values, &
          'External Temperature', en, Edge % NodeIndexes, Found )

!       ...black body radiation:
!       ------------------------
        Emissivity      = 0.0d0
        StefanBoltzmann = 0.0d0

        SELECT CASE(ListGetString(Model % BCs(j) % Values,'Radiation',Found))
           !------------------
           CASE( 'idealized' )
           !------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                  ListGetConstReal( Model % Constants,'Stefan Boltzmann',UnfoundFatal=.TRUE. )

           !---------------------
           CASE( 'diffuse gray' )
           !---------------------

              NodalEmissivity(1:en) = ListGetReal( Model % BCs(j) % Values, &
                   'Emissivity', en, Edge % NodeIndexes, Found)
              IF(.NOT. Found) THEN
                 NodalEmissivity(1:en) = GetParentMatProp( 'Emissivity', Edge)
              END IF
              Emissivity = SUM( NodalEmissivity(1:en)) / en

              StefanBoltzMann = &
                    ListGetConstReal( Model % Constants,'Stefan Boltzmann' )

              ExtTemperature(1:en) =  ComputeRadiationLoad( Model, &
                      Mesh, Edge, Quant, Perm, Emissivity )
        END SELECT

!       get material parameters:
!       ------------------------
        k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
                    minv=1, maxv=Model % NumberOFMaterials)

        CALL ListGetRealArray( Model % Materials(k) % Values, &
               'Heat Conductivity', Hwrk, en, Edge % NodeIndexes )
        NodalConductivity(1:en) = Hwrk(1,1,1:en)

!       elementwise nodal solution:
!       ---------------------------
        nd = GetElementDOFs(Indexes,Element)
        Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))

!       do the integration:
!       -------------------
        EdgeLength   = 0.0d0
        ResidualNorm = 0.0d0

        IntegStuff = GaussPoints(Edge)

        DO t=1,IntegStuff % n
           u = IntegStuff % u(t)
           v = IntegStuff % v(t)
           w = IntegStuff % w(t)

           stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
                    EdgeBasis, dBasisdx )
           Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )

           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              s = IntegStuff % s(t) * detJ
           ELSE
              gx = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
              gy = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
              gz = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )
              CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                         Symb, dSymb, gx, gy, gz )
              s = IntegStuff % s(t) * detJ * SqrtMetric
           END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
           u = SUM( EdgeBasis(1:en) * x(1:en) )
           v = SUM( EdgeBasis(1:en) * y(1:en) )
           w = SUM( EdgeBasis(1:en) * z(1:en) )
           stat = ElementInfo(Element,Nodes, u, v, w, detJ,Basis,dBasisdx )

!
!          Heat conductivity at the integration point:
!          --------------------------------------------
           Conductivity = SUM( NodalConductivity(1:en) * EdgeBasis(1:en) )
!
!          given flux at integration point:
!          --------------------------------
           Residual = -SUM( Flux(1:en) * EdgeBasis(1:en) )

!          convective ...:
!          ----------------
           Residual = Residual + SUM(TransferCoeff(1:en) * EdgeBasis(1:en)) * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) )

!          black body radiation...:
!          -------------------------
           Residual = Residual + &
                Emissivity * StefanBoltzmann * &
                     ( SUM( Temperature(1:nd) * Basis(1:nd) ) ** 4 - &
                       SUM( ExtTemperature(1:en) * EdgeBasis(1:en) ) ** 4 )

!          flux given by the computed solution, and 
!          force norm for scaling the residual:
!          -----------------------------------------
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              DO k=1,dim
                 Residual = Residual + Conductivity  * &
                    SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(k)

                 Gnorm = Gnorm + s * (Conductivity * &
                       SUM(dBasisdx(1:nd,k) * Temperature(1:nd)) * Normal(k))**2
              END DO
           ELSE
              DO k=1,dim
                 DO l=1,dim
                    Residual = Residual + Metric(k,l) * Conductivity  * &
                       SUM( dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l)

                    Gnorm = Gnorm + s * (Metric(k,l) * Conductivity * &
                      SUM(dBasisdx(1:nd,k) * Temperature(1:nd) ) * Normal(l))**2
                 END DO
              END DO
           END IF

           EdgeLength   = EdgeLength + s
           IF ( .NOT. Dirichlet ) THEN
              ResidualNorm = ResidualNorm + s * Residual ** 2
           END IF
        END DO
        EXIT
     END DO

     IF ( CoordinateSystemDimension() == 3 ) EdgeLength = SQRT(EdgeLength)

!    Gnorm = EdgeLength * Gnorm
     Indicator = EdgeLength * ResidualNorm
!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Boundary_Residual
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  FUNCTION HeatSolver_Edge_Residual(Model,Edge,Mesh,Quant,Perm) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE

     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2)
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
     TYPE(Nodes_t) :: Nodes, EdgeNodes
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,l,n,t,dim,En,Pn,nd
     INTEGER, ALLOCATABLE :: Indexes(:)
     LOGICAL :: stat
     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), &
            EdgeBasis(:), Basis(:), dBasisdx(:,:), Temperature(:)
     REAL(KIND=dp) :: Grad(3,3), Normal(3), EdgeLength, Jump, Conductivity
     REAL(KIND=dp) :: u, v, w, s, detJ
     REAL(KIND=dp) :: ResidualNorm
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

     !    Initialize:
     !    -----------
     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

     Metric = 0.0d0
     DO i = 1,3
        Metric(i,i) = 1.0d0
     END DO
     Grad = 0.0d0
!
!    ---------------------------------------------

     n = Mesh % MaxElementDOFs
     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

     en = Edge % TYPE % NumberOfNodes
     ALLOCATE( EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en) )

     EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
     EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
     EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

     ALLOCATE( NodalConductivity(en), EdgeBasis(en), Basis(n), &
        dBasisdx(n,3), x(en), y(en), z(en), Temperature(n), Indexes(n) )

!    Integrate square of jump over edge:
!    -----------------------------------
     ResidualNorm = 0.0d0
     EdgeLength   = 0.0d0
     Indicator    = 0.0d0

     IntegStuff = GaussPoints( Edge )

     DO t=1,IntegStuff % n

        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
             EdgeBasis, dBasisdx )

        Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( EdgeBasis(1:en) * EdgeNodes % x(1:en) )
           v = SUM( EdgeBasis(1:en) * EdgeNodes % y(1:en) )
           w = SUM( EdgeBasis(1:en) * EdgeNodes % z(1:en) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        ! 
        ! Compute flux over the edge as seen by elements
        ! on both sides of the edge:
        ! ----------------------------------------------
        DO i = 1,2
           SELECT CASE(i)
              CASE(1)
                 Element => Edge % BoundaryInfo % Left
              CASE(2)
                 Element => Edge % BoundaryInfo % Right
           END SELECT
!
!          Can this really happen (maybe it can...)  ?      
!          -------------------------------------------
           IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
           pn = Element % TYPE % NumberOfNodes

           DO j = 1,en
              DO k = 1,pn
                 IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
                    x(j) = Element % TYPE % NodeU(k)
                    y(j) = Element % TYPE % NodeV(k)
                    z(j) = Element % TYPE % NodeW(k)
                    EXIT
                 END IF
              END DO
           END DO

           u = SUM(EdgeBasis(1:en) * x(1:en))
           v = SUM(EdgeBasis(1:en) * y(1:en))
           w = SUM(EdgeBasis(1:en) * z(1:en))
!
!          Get parent element basis & derivatives at the integration point:
!          -----------------------------------------------------------------
           nd = GetElementDOFs(Indexes,Element)
           Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
           Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
           Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))

           stat = ElementInfo(Element,Nodes,u,v,w,detJ,Basis,dBasisdx)
!
!          Material parameters:
!          --------------------
           k = ListGetInteger( Model % Bodies( &
                    Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOFMaterials )

           CALL ListGetRealArray( Model % Materials(k) % Values, &
                   'Heat Conductivity', Hwrk,en, Edge % NodeIndexes )

           NodalConductivity(1:en) = Hwrk( 1,1,1:en )
           Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
!
!          Temperature at element nodal points:
!          ------------------------------------
           Temperature(1:nd) = Quant( Perm(Indexes(1:nd)) )
!
!          Finally, the flux:
!          ------------------
           DO j=1,dim
              Grad(j,i) = Conductivity * SUM( dBasisdx(1:nd,j) * Temperature(1:nd) )
           END DO
        END DO

!       Compute squre of the flux jump:
!       -------------------------------   
        EdgeLength  = EdgeLength + s
        Jump = 0.0d0
        DO k=1,dim
           IF ( CurrentCoordinateSystem() == Cartesian ) THEN
              Jump = Jump + (Grad(k,1) - Grad(k,2)) * Normal(k)
           ELSE
              DO l=1,dim
                 Jump = Jump + &
                       Metric(k,l) * (Grad(k,1) - Grad(k,2)) * Normal(l)
              END DO
           END IF
        END DO
        ResidualNorm = ResidualNorm + s * Jump ** 2
     END DO

     IF (dim==3) EdgeLength = SQRT(EdgeLength)
     Indicator = EdgeLength * ResidualNorm
!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Edge_Residual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION HeatSolver_Inside_Residual( Model, Element, Mesh, &
        Quant, Perm, Fnorm ) RESULT( Indicator )
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,n,t,dim,nd
     INTEGER, ALLOCATABLE :: Indexes(:)

     LOGICAL :: stat, Found, Compressible, VolSource
     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp), POINTER :: Hwrk(:,:,:)

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)

     REAL(KIND=dp), ALLOCATABLE :: NodalDensity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalCapacity(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
     REAL(KIND=dp), ALLOCATABLE :: Velo(:,:), Pressure(:)
     REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Temperature(:), PrevTemp(:)
     REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:)

     REAL(KIND=dp) :: u, v, w, s, detJ, Density, Capacity

     REAL(KIND=dp) :: SpecificHeatRatio, ReferencePressure, dt
     REAL(KIND=dp) :: Residual, ResidualNorm, Area, Conductivity

     TYPE( ValueList_t ), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: First = .TRUE.
     SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     Fnorm     = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
     IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN

     IF ( First ) THEN
        First = .FALSE.
        NULLIFY( Hwrk )
     END IF

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           dim = 3
        CASE DEFAULT
           dim = CoordinateSystemDimension()
     END SELECT

!    Alllocate local arrays
!    ----------------------
     nd = GetElementNOFDOFs(Element)
     n = GetElementNOFNodes(Element)
     ALLOCATE( NodalDensity(nd), NodalCapacity(nd), NodalConductivity(nd),      &
         Velo(3,nd), Pressure(nd), NodalSource(nd), Temperature(nd), PrevTemp(nd), &
         Basis(nd), dBasisdx(nd,3), ddBasisddx(nd,3,3), Indexes(nd) )
!
!    Element nodal points:
!    ---------------------
     ALLOCATE( Nodes % x(nd), Nodes % y(nd), Nodes % z(nd) )

     nd = GetElementDOFs(Indexes,Element)
     Nodes % x = Mesh % Nodes % x(Indexes(1:nd))
     Nodes % y = Mesh % Nodes % y(Indexes(1:nd))
     Nodes % z = Mesh % Nodes % z(Indexes(1:nd))
!
!    Elementwise nodal solution:
!    ---------------------------
     Temperature(1:nd) = Quant(Perm(Indexes(1:nd)))
!
!    Check for time dep.
!    -------------------
     PrevTemp(1:nd) = Temperature(1:nd)
     dt = Model % Solver % dt
     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Temperature', .TRUE. )
        PrevTemp(1:nd) = Var % PrevValues(Var % Perm(Indexes(1:nd)),1)
     END IF
!
!    Material parameters: conductivity, heat capacity and density
!    -------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     CALL ListGetRealArray(Material,'Heat Conductivity',Hwrk,n,Element % NodeIndexes)
     NodalConductivity(1:n) = Hwrk( 1,1,1:n)
     NodalDensity(1:n) = GetReal(Material, 'Density', Found )
     NodalCapacity(1:n) = GetReal(Material, 'Heat Capacity', Found )
!
!    Check for compressible flow equations:
!    --------------------------------------
     Compressible = .FALSE.

     IF (  ListGetString( Material, 'Compressibility Model', Found ) == &
                 'perfect gas equation 1' ) THEN

        Compressible = .TRUE.

        Pressure = 0.0d0
        Var => VariableGet( Mesh % Variables, 'Pressure', .TRUE. )
        IF ( ASSOCIATED( Var ) ) THEN
           Pressure(1:n) = Var % Values(Var % Perm(Indexes(1:n)))
        END IF

        ReferencePressure = GetConstReal( Material, 'Reference Pressure' )
        SpecificHeatRatio = GetConstReal( Material, 'Specific Heat Ratio' )
        NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
              ( (SpecificHeatRatio - 1) * NodalCapacity(1:n) * Temperature(1:n) )
     END IF
!
!    Get (possible) convection velocity at the nodes of the element:
!    ----------------------------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                minv=1, maxv=Model % NumberOFEquations )

     Velo = 0.0d0
     SELECT CASE( ListGetString( Model % Equations(k) % Values, &
                         'Convection', Found ) )

        !-----------------
        CASE( 'constant' )
        !-----------------

           Velo(1,1:n) = GetReal( Material, 'Convection Velocity 1',Found )
           Velo(2,1:n) = GetReal( Material, 'Convection Velocity 2',Found )
           Velo(3,1:n) = GetReal( Material, 'Convection Velocity 3',Found )

        !-----------------
        CASE( 'computed' )
        !-----------------

           Var => VariableGet( Mesh % Variables, 'Velocity 1', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              IF ( ALL(Var % Perm(Indexes(1:n)) > 0 ) ) THEN
                 Velo(1,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
                 Var => VariableGet( Mesh % Variables, 'Velocity 2', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(2,1:n) = Var % Values(Var % Perm(Indexes(1:n)) )
                 Var => VariableGet( Mesh % Variables, 'Velocity 3', .TRUE. )
                 IF ( ASSOCIATED( Var ) ) &
                    Velo(3,1:n) = Var % Values(Var % Perm(Indexes(1:n)))
              END IF
           END IF
     END SELECT

!
!    Heat source:
!    ------------
!
     k = ListGetInteger( &
         Model % Bodies(Element % BodyId) % Values,'Body Force',Found, &
                 1, Model % NumberOFBodyForces)

     NodalSource = 0.0d0
     IF( k > 0 ) THEN
       NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
           'Volumetric Heat Source',VolSource ) 
       IF( .NOT. VolSource ) THEN
         NodalSource(1:n) = GetReal( Model % BodyForces(k) % Values, &
             'Heat Source',  Found )
       END IF
     END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

     ResidualNorm = 0.0d0
     Area = 0.0d0

     IntegStuff = GaussPoints( Element )
     ddBasisddx = 0

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM(Basis(1:nd) * Nodes % x(1:nd))
           v = SUM(Basis(1:nd) * Nodes % y(1:nd))
           w = SUM(Basis(1:nd) * Nodes % z(1:nd))

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                       Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Capacity     = SUM(NodalCapacity(1:n) * Basis(1:n))
        Density      = SUM(NodalDensity(1:n) * Basis(1:n))
        Conductivity = SUM(NodalConductivity(1:n) * Basis(1:n))
!
!       Residual of the convection-diffusion (heat) equation:
!        R = \rho * c_p * (@T/@t + u.grad(T)) - &
!            div(C grad(T)) + p div(u) - h,
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = \rho * c_p * (@T/@t + u^j T_{,j}) - &
!          g^{jk} (C T_{,j}}_{,k} + p div(u) - h
!       ---------------------------------------------------
!

        IF( VolSource ) THEN
          Residual = -SUM( NodalSource(1:n) * Basis(1:n) )
        ELSE
          Residual = -Density * SUM( NodalSource(1:n) * Basis(1:n) )
        END IF
          
        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           DO j=1,dim
!
!             - grad(C).grad(T):
!             --------------------
!
              Residual = Residual - &
                 SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                 SUM( NodalConductivity(1:n) * dBasisdx(1:n,j) )

!
!             - C div(grad(T)):
!             -------------------
!
              Residual = Residual - Conductivity * &
                 SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,j) )
           END DO
        ELSE
           DO j=1,dim
              DO k=1,dim
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
                 Residual = Residual - Metric(j,k) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) ) * &
                    SUM( NodalConductivity(1:n) * dBasisdx(1:n,k) )

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
                 Residual = Residual - Metric(j,k) * Conductivity * &
                    SUM( Temperature(1:nd) * ddBasisddx(1:nd,j,k) )
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
                 DO l=1,dim
                    Residual = Residual + Metric(j,k) * Conductivity * &
                      Symb(j,k,l) * SUM( Temperature(1:nd) * dBasisdx(1:nd,l) )
                 END DO
              END DO
           END DO
        END IF

!       + \rho * c_p * (@T/@t + u.grad(T)):
!       -----------------------------------
        Residual = Residual + Density * Capacity *  &
           SUM((Temperature(1:nd)-PrevTemp(1:nd))*Basis(1:nd)) / dt

        DO j=1,dim
           Residual = Residual + &
              Density * Capacity * SUM(Velo(j,1:n) * Basis(1:n)) * &
                    SUM( Temperature(1:nd) * dBasisdx(1:nd,j) )
        END DO


        IF ( Compressible ) THEN
!
!          + p div(u) or p u^j_{,j}:
!          -------------------------
!
           DO j=1,dim
              Residual = Residual + &
                 SUM( Pressure(1:n) * Basis(1:n) ) * &
                      SUM( Velo(j,1:n) * dBasisdx(1:n,j) )

              IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
                 DO k=1,dim
                    Residual = Residual + &
                       SUM( Pressure(1:n) * Basis(1:n) ) * &
                           Symb(j,k,j) * SUM( Velo(k,1:n) * Basis(1:n) )
                 END DO
              END IF
           END DO
        END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
!       DO i=1,dim
!          Fnorm = Fnorm + s * (Density *SUM(NodalSource(1:n)*Basis(1:n)))**2
!       END DO
        Area = Area + s
        ResidualNorm = ResidualNorm + s *  Residual ** 2
     END DO

!    Fnorm = Element % hk**2 * Fnorm
     Indicator = Element % hK**2 * ResidualNorm
!------------------------------------------------------------------------------
   END FUNCTION HeatSolver_Inside_Residual
!------------------------------------------------------------------------------
