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
! *  Module for static electric fields. 
! *  Based on the multithreaded and vectorized StatCurrentSolveVec that is based om
! * ModelPDE by Mikko Byckling. Replaces gradually the old StatElecSolver that is
! * not optimized.  
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 20.4.2021
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. StatElecSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, CalculateElemental, CalculateNodal
  INTEGER :: dim
   
  Params => GetSolverParams()
  dim = CoordinateSystemDimension()

  CALL ListAddNewString( Params,'Variable','Potential')

  CalculateElemental = ListGetLogical( Params,'Calculate Elemental Fields',Found )
  CalculateNodal = ListGetLogical( Params,'Calculate Nodal Fields',Found )

  IF(.NOT. (CalculateElemental .OR. CalculateNodal ) ) THEN
    CalculateNodal = .TRUE.
  END IF

  IF (ListGetLogical(Params,'Calculate Electric Energy',Found)) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Electric Energy Density e' )
    IF( CalculateNodal ) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Electric Energy Density' )
  END IF

  IF( ListGetLogical(Params,'Calculate Elecric Flux',Found) ) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Elecric Flux e[Elecric Flux e:'//TRIM(I2S(dim))//']' )
    IF( CalculateNodal ) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Elecric Flux[Elecric Flux:'//TRIM(I2S(dim))//']' )       
  END IF

  IF( ListGetLogical(Params,'Calculate Electric Field',Found) ) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Electric Field e[Electric Field e:'//TRIM(I2S(dim))//']' )
    IF( CalculateNodal ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Electric Field[Electric Field:'//TRIM(I2S(dim))//']' )
  END IF

  ! Nodal fields that may directly be associated as nodal loads
  IF (ListGetLogical(Params,'Calculate Nodal Energy',Found))  THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
        'Nodal Energy Density' )
  END IF
  IF( ListGetLogical(Params,'Calculate Nodal Flux',Found) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Nodal Electric Flux[Nodal Electric Flux:'//TRIM(I2S(dim))//']' )
  END IF

  ! These use one flag to call library features to compute automatically
  ! a capacitance matrix.
  IF( ListGetLogical(Params,'Calculate Capacitance Matrix',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes Symmetric',.TRUE.)
    IF( ListCheckPresent( Params,'Capacitance Matrix Filename') ) THEN
      CALL ListRename( Params,'Capacitance Matrix Filename',&
          'Constraint Modes Fluxes Filename', Found ) 
    ELSE     
      CALL ListAddNewString( Params,'Constraint Modes Fluxes Filename',&
          'CapacitanceMatrix.dat',.FALSE.)
    END IF
    CALL ListRenameAllBC( Model,'Capacitance Body','Constraint Mode Potential')
  END IF

  CALL ListAddInteger( Params,'Time Derivative Order', 0 )
  
  CALL ListWarnUnsupportedKeyword('solver','adaptive mesh redinement',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('body force','piezo material',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('boundary condition','Layer Relative Permittivity',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('boundary condition','infinity bc',FatalFound=.TRUE.)

END SUBROUTINE StatElecSolver_Init


!-----------------------------------------------------------------------------
!> A modern version for static current conduction supporting multithreading and
!> SIMD friendly ElmerSolver kernels. 
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
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
  INTEGER :: n, nb, nd, t, active, dim, RelOrder
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr
  LOGICAL :: Found, VecAsm, InitHandles, AxiSymmetric
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver'
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving static electric field for insulators')

  CALL DefaultStart()

  Mesh => GetMesh()
  Params => GetSolverParams()
  
  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian ) 
  dim = CoordinateSystemDimension() 

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
    CALL Info(Caller,'Vectorized loop not yet available in axisymmetric case',Level=7)    
    VecAsm = .FALSE.
  END IF

  IF( VecAsm ) THEN
    CALL Info(Caller,'Performing vectorized bulk element assembly',Level=7)
  ELSE
    CALL Info(Caller,'Performing non-vectorized bulk element assembly',Level=7)      
  END IF

  RelOrder = GetInteger( Params,'Relative Integration Order',Found ) 

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()

    totelem = 0

    CALL ResetTimer( Caller//'BulkAssembly' )

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

    CALL CheckTimer(Caller//'BulkAssembly',Delete=.TRUE.)
    totelem = 0

    CALL DefaultFinishBulkAssembly()

    nColours = GetNOFBoundaryColours(Solver)

    CALL Info(Caller,'Performing boundary element assembly',Level=12)
    CALL ResetTimer(Caller//'BCAssembly')

    !$OMP PARALLEL &
    !$OMP SHARED(Active, Solver, nColours, VecAsm) &
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
          CALL LocalMatrixBC(  Element, n, nd+nb, nb, VecAsm, InitHandles )
        END IF
      END DO
      !$OMP END DO
    END DO
    !$OMP END PARALLEL

    CALL CheckTimer(Caller//'BCAssembly',Delete=.TRUE.)

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged == 1 ) EXIT

  END DO

  CALL DefaultFinish()

  
CONTAINS

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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, POINTER  :: EpsAtIpVec(:), SourceAtIpVec(:)
    REAL(KIND=dp) :: eps0, weight
    LOGICAL :: Stat,Found, Pref
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, DetJVec, &
    !$OMP               STIFF, FORCE, Nodes, &
    !$OMP               SourceCoeff_h, EpsCoeff_h, &
    !$OMP               SourceAtIpVec, EpsAtIpVec )
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJVec
    !DIR$ ATTRIBUTES ALIGN:64 :: STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Source')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity')
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      END IF
      IF( .NOT. Found ) Eps0 = 8.854187817e-12
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()

    ! Use p-element basis, unless 2nd or 3rd order nodal element
    Pref = isPElement(Element) .OR. Element % Type % BasisFunctionDegree<2
    IF( RelOrder /= 0 ) THEN
      IP = GaussPoints( Element, PReferenceElement = Pref, RelOrder = RelOrder)
    ELSE
      IP = GaussPoints( Element, PReferenceElement = Pref )
    END IF
      
    ngp = IP % n

    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
      IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
            DEALLOCATE(Basis, dBasisdx, DetJVec, STIFF, FORCE )
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJVec(ngp), &
          STIFF(nd,nd), FORCE(nd), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodesVec( Nodes, UElement=Element )

    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp

    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJvec, &
        SIZE(Basis,2), Basis, dBasisdx )

    ! Compute actual integration weights (recycle the memory space of DetJVec)
    DO i=1,ngp
      DetJVec(i) = IP % s(i) * DetJVec(i)
    END DO

    ! diffusivity term: STIFF=STIFF+(eps*grad(u),grad(v))
    EpsAtIpVec => ListGetElementRealVec( EpsCoeff_h, ngp, Basis, Element )
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, &
        DetJVec, STIFF, EpsAtIpVec )
    STIFF(1:nd,1:nd) = Eps0 * STIFF(1:nd,1:nd)

    ! source term: FORCE=FORCE+(u,f)
    SourceAtIpVec => ListGetElementRealVec( SourceCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJVec, SourceAtIpVec, FORCE)
    END IF
      
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )

    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec


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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, DetJ, A
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,p,q,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Source')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity')
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      END IF
      IF( .NOT. Found ) Eps0 = 8.854187817e-12
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()

    IF( RelOrder /= 0 ) THEN
      IP = GaussPoints( Element, RelOrder = RelOrder)
    ELSE
      IP = GaussPoints( Element )
    END IF
      
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDofs
      ALLOCATE(Basis(m), dBasisdx(m,3),&
          STIFF(m,m), FORCE(m), STAT=allocstat)
      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      IF ( AxiSymmetric ) THEN
        Weight = Weight * 2 * PI * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      EpsAtIp = ListGetElementReal( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )      
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
          Eps0 * EpsAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )

      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      END IF
    END DO
    
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, VecAsm, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: VecAsm
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,Ext, Weight,Eps0
    REAL(KIND=dp) :: Basis(nd),DetJ,Coord(3),Normal(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,RobinBC
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Robin_h, Ext_h, Farfield_h, Infty_h

    SAVE Nodes, Eps0
    !$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h,Farfield_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Electric Flux')
      CALL ListInitElementKeyword( Infty_h,'Boundary Condition','Electric Infinity BC')
      CALL ListInitElementKeyword( Farfield_h,'Boundary Condition','Farfield Potential')
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      END IF
      IF( .NOT. Found ) Eps0 = 8.854187817e-12 
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
           

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
        Weight = Weight * 2 * PI * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = ListGetElementReal( Flux_h, Basis, Element, Found )
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
      END IF

      ! Robin type of condition for farfield potential (r*(u-u_0)):
      ! -----------------------------------------------------------
      IF( ListGetElementLogical( Infty_h, Element, Found ) ) THEN
        Coord(1) = SUM( Nodes % x(1:n)*Basis(1:n) )
        Coord(2) = SUM( Nodes % y(1:n)*Basis(1:n) )
        Coord(3) = SUM( Nodes % z(1:n)*Basis(1:n) )
        
        Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t), .TRUE. )
        C = Eps0 * SUM( Coord * Normal ) / SUM( Coord * Coord )         
        
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        Ext = ListGetElementReal( Farfield_h, Basis, Element, Found )
        IF( Found ) THEN
          FORCE(1:nd) = FORCE(1:nd) + Weight * C * Ext * Basis(1:nd)
        END IF
      END IF
    END DO
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE StatElecSolver
!------------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!> A solver performing the postprocessing for the primary solver.
!> This solver is called at the DefaultFinish() slot of the primary solver.
!------------------------------------------------------------------------------
SUBROUTINE StatElecSolver_post( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
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
  INTEGER :: i, dofs, n, nb, nd, t, active
  LOGICAL :: Found, InitHandles
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: WeightVector(:),FORCE(:,:),MASS(:,:),&
      PotInteg(:),PotVol(:)
  INTEGER, POINTER :: WeightPerm(:)
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver_post'
  LOGICAL :: CalcCurrent, CalcField, CalcElectricDisp, NeedScaling, ConstantWeights, &
      Axisymmetric, CalcAvePotential
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) :: EnergyTot, Voltot
  
  TYPE PostVars_t
    TYPE(Variable_t), POINTER :: Var => NULL()
    INTEGER :: FieldType = -1
    LOGICAL :: NodalField = .FALSE.
    LOGICAL :: HaveVar = .FALSE.
  END TYPE PostVars_t
  TYPE(PostVars_t) :: PostVars(8)
  !------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Calculating postprocessing fields')
  
  Mesh => GetMesh()
  Params => GetSolverParams()
    
  ConstantWeights = ListGetLogical( Params,'Constant Weights',Found ) 

  AxiSymmetric = ( CurrentCoordinateSystem() /= Cartesian )     

  CalcAvePotential = ListGetLogical( Params,'Calculate Average Potential',Found )

  IF( CalcAvePotential ) THEN   
    n = Model % NumberOfBodies 
    ALLOCATE( PotInteg(n), PotVol(n) )
    PotInteg = 0.0_dp
    PotVol = 0.0_dp
  END IF
  
  ! Joule losses: type 1, component 1
  PostVars(1) % Var => VariableGet( Mesh % Variables, 'Nodal Energy Density')
  PostVars(1) % NodalField = .TRUE.
  PostVars(2) % Var => VariableGet( Mesh % Variables, 'Electric Energy Density')
  PostVars(3) % Var => VariableGet( Mesh % Variables, 'Electric Energy Density e')
  PostVars(1:3) % FieldType = 1

  ! Electric current: type 2, components 2:4
  PostVars(4) % Var => VariableGet( Mesh % Variables, 'Nodal Electric Flux')
  PostVars(4) % NodalField = .TRUE.
  PostVars(5) % Var => VariableGet( Mesh % Variables, 'Elecric Flux')
  PostVars(6) % Var => VariableGet( Mesh % Variables, 'Elecric Flux e')
  PostVars(4:6) % FieldType = 2 

  ! Electric field: type 3, components 5-7
  PostVars(7) % Var => VariableGet( Mesh % Variables, 'Electric Field')
  PostVars(8) % Var => VariableGet( Mesh % Variables, 'Electric Field e')
  PostVars(7:8) % FieldType = 3 

  ! Do this since the "associated" command cannot handle vectors!
  ! Also initialize the field to zero since some of these are additive
  DO i=1,8
    PostVars(i) % HaveVar = ASSOCIATED( PostVars(i) % Var )
    IF( PostVars(i) % HaveVar ) PostVars(i) % Var % Values = 0.0_dp
  END DO
  
  CalcElectricDisp = ANY( PostVars(1:3) % HaveVar ) 
  CalcCurrent = ANY( PostVars(4:6) % HaveVar ) 
  CalcField = ANY( PostVars(7:8) % HaveVar ) 

  n = COUNT( PostVars(1:8) % HaveVar )
  CALL Info(Caller,'Number of '//TRIM(I2S(n))//' postprocessing fields',Level=8)
    
  ! Only create the nodal weights if we need to scale some nodal field
  NeedScaling = .FALSE.
  DO i=1,8
    IF( .NOT. PostVars(i) % HaveVar ) CYCLE
    IF( PostVars(i) % NodalField ) CYCLE
    IF( PostVars(i) % Var % TYPE == Variable_on_nodes ) THEN
      CALL Info(Caller,'Creating a weighting for scaling purposes from '//TRIM(I2S(i)),Level=10)
      NeedScaling = .TRUE.
      WeightPerm => PostVars(i) % Var % Perm 
      ALLOCATE( WeightVector( MAXVAL( WeightPerm ) ) )
      WeightVector = 0.0_dp
      EXIT
    END IF
  END DO

  n = Mesh % MaxElementDOFs
  ALLOCATE( MASS(n,n), FORCE(8,n) ) ! 1+1+3+3 components for force

  CALL Info(Caller,'Calculating local field values',Level=12) 
  InitHandles = .TRUE.
  EnergyTot = 0.0_dp
  VolTot = 0.0_dp
  DO t = 1, GetNOFActive()
    Element => GetActiveElement(t)
    IF( ParEnv % PEs > 1 ) THEN
      IF( ParEnv % MyPe /= Element % PartIndex ) CYCLE
    END IF
    n  = GetElementNOFNodes(Element)
    CALL LocalPostAssembly( Element, n, InitHandles, MASS, FORCE )
    CALL LocalPostSolve( Element, n, MASS, FORCE )
  END DO
  
  IF( NeedScaling ) THEN
    CALL Info(Caller,'Scaling the field values with weights',Level=12)
    CALL GlobalPostScale()
  END IF

  IF( CalcAvePotential ) THEN
    BLOCK        
      REAL(KIND=dp), ALLOCATABLE:: PotTmp(:)
      INTEGER :: ierr      
      REAL(KIND=dp) :: PotAve
      n = Model % NumberOfBodies
      IF( ParEnv % PEs > 1 ) THEN
        ALLOCATE( PotTmp(n) )        
        CALL MPI_ALLREDUCE(PotVol,PotTmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,ParEnv % ActiveComm,ierr)
        PotVol = PotTmp
        CALL MPI_ALLREDUCE(PotInteg,PotTmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,ParEnv % ActiveComm,ierr)
        PotInteg = PotTmp
        DEALLOCATE( PotTmp ) 
      END IF

      DO i = 1, n
        IF( PotVol(i) < EPSILON( PotVol(i) ) ) CYCLE
        PotAve = PotInteg(i) / PotVol(i)
        WRITE( Message,'(A,ES12.5)') 'Average body'//TRIM(I2S(i))//' potential: ',PotAve
        CALL Info(Caller,Message,Level=7)
        CALL ListAddConstReal( Model % Simulation,&
            'res: Average body'//TRIM(I2S(i))//' potential',PotAve)
      END DO
    END BLOCK
  END IF
    
  CALL Info(Caller,'All done',Level=12)
  

CONTAINS
   
  SUBROUTINE LocalPostAssembly( Element, n, InitHandles, MASS, FORCE )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
    REAL(KIND=dp) :: MASS(:,:), FORCE(:,:)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:),ElementPot(:)
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, DetJ
    REAL(KIND=dp) :: EpsGrad(3), Grad(3), Heat
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,p,q,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
    
!------------------------------------------------------------------------------
    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Source')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity')
      Found = .FALSE.
      IF( ASSOCIATED( Model % Constants ) ) THEN
        Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      END IF
      IF( .NOT. Found ) Eps0 = 8.854187817e-12
      InitHandles = .FALSE.
    END IF

    dim = CoordinateSystemDimension()

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = Mesh % MaxElementDOFs   
      ALLOCATE(Basis(m), dBasisdx(m,3), ElementPot(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )
    CALL GetScalarLocalSolution( ElementPot ) 
    
    ! Initialize
    MASS  = 0._dp
    FORCE = 0._dp

    IP = GaussPoints( Element )
    
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detJ, Basis, dBasisdx )
      Weight = IP % s(t) * DetJ

      IF ( AxiSymmetric ) THEN
        Weight = Weight * 2 * PI * SUM( Nodes % x(1:n)*Basis(1:n) )
      END IF
       
      DO i=1,n
        DO j=1,n
          MASS(i,j) = MASS(i,j) + Weight * Basis(i) * Basis(j)
        END DO
      END DO

      ! Compute the integration weights 
      !----------------------------------------------------------------------------
      FORCE(1,1:n) = FORCE(1,1:n) + Weight * Basis(1:n)

      EpsAtIp = Eps0 * ListGetElementReal( EpsCoeff_h, Basis, Element, Found, &
         GaussPoint = t )
        
      ! Compute the electric field from the potential: E = -grad Phi
      !------------------------------------------------------------------------------
      DO j = 1, DIM
        Grad(j) = SUM( dBasisdx(1:n,j) * ElementPot(1:n) )
      END DO
      IF( CalcField ) THEN
        DO j=1,dim
          Force(5+j,1:n) = Force(5+j,1:n) - Grad(j) * Weight * Basis(1:n)
        END DO
      END IF

      ! Compute the volume current: J = eps (-grad Phi)
      !------------------------------------------------------------------------------
      EpsGrad(1:dim) = EpsAtIp * Grad(1:dim)

      IF( CalcCurrent ) THEN
        DO j=1,dim
          Force(2+j,1:n) = Force(2+j,1:n) - EpsGrad(j) * Weight * Basis(1:n)
        END DO
      END IF

      ! Compute the Joule heating: H,tot = Integral (E . D)dV
      !------------------------------------------------------------------------------
      Heat = SUM( Grad(1:dim) * EpsGrad(1:dim) )      
      IF( CalcElectricDisp ) THEN
        Force(2,1:n) = Force(2,1:n) + Heat * Weight * Basis(1:n)
      END IF

      IF( CalcAvePotential ) THEN
        i = Element % BodyId
        PotVol(i) = PotVol(i) + Weight
        PotInteg(i) = PotInteg(i) + Weight * SUM( Basis(1:n) * ElementPot(1:n) ) 
      END IF
      
      VolTot = VolTot + Weight
      EnergyTot = EnergyTot + Weight * Heat 
    END DO

    IF( NeedScaling ) THEN
      IF( ConstantWeights ) THEN
        WeightVector( WeightPerm( Element % NodeIndexes ) ) = &
            WeightVector( WeightPerm( Element % NodeIndexes ) ) + 1.0_dp
      ELSE
        WeightVector( WeightPerm( Element % NodeIndexes ) ) = &
            WeightVector( WeightPerm( Element % NodeIndexes ) ) + Force(1,1:n)
      END IF
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE LocalPostAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalPostSolve( Element, n, A, b )
!------------------------------------------------------------------------------    
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    REAL(KIND=dp) :: b(:,:), A(:,:)
!------------------------------------------------------------------------------
    INTEGER :: pivot(n),ind(n),i,j,m,dofs,dofcount,FieldType,Vari
    REAL(KIND=dp) :: x(n)
    TYPE(Variable_t), POINTER :: pVar
    LOGICAL :: LocalSolved
!------------------------------------------------------------------------------
    
    CALL LUdecomp(A,n,pivot)

    ! Weight is the 1st column
    dofcount = 1
    DO FieldType = 1, 3

      IF( FieldType == 1 ) THEN
        ! Joule heating has one component
        dofs = 1
      ELSE
        ! Current and electric field has three components
        dofs = 3
      END IF
     
      DO m=1,dofs
        dofcount = dofcount+1
        x = b(dofcount,1:n)
        LocalSolved = .FALSE.
        
        DO Vari = 1, 8
          pVar => PostVars(Vari) % Var
          IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
          IF( PostVars(Vari) % FieldType /= FieldType ) CYCLE
          IF( m > pVar % Dofs ) CYCLE
          
          ! The nodal fields need not be solved for.
          ! Note the nodal field should come before the distributed fields!!
          IF( PostVars(Vari) % NodalField ) THEN
            CONTINUE
          ELSE IF(.NOT. LocalSolved ) THEN
            CALL LUSolve(n,MASS,x,pivot)
            LocalSolved = .TRUE.
          END IF

          ! Note that even though while calling we implicitly assumes elemental
          ! and nodal fields the convention is not assumed here.
          IF( pVar % TYPE == variable_on_nodes_on_elements ) THEN
            ind = pVar % dofs * (pVar % Perm(Element % DGIndexes(1:n))-1)+m
            pVar % Values(ind(1:n)) = x(1:n)          
          ELSE IF( pVar % TYPE == variable_on_nodes ) THEN
            ind = pVar % dofs * (pVar % Perm(Element % NodeIndexes(1:n))-1)+m
            IF( ConstantWeights ) THEN
              pVar % Values(ind(1:n)) = pVar % Values(ind(1:n)) + x(1:n)                    
            ELSE
              pVar % Values(ind(1:n)) = pVar % Values(ind(1:n)) + b(1,1:n) * x(1:n)                               
            END IF
          ELSE IF( pVar % TYPE == variable_on_elements ) THEN
            j = pVar % dofs * ( pVar % Perm( Element % ElementIndex )-1)+m
            pVar % Values(j) = SUM( x(1:n) ) / n
          ELSE
            CALL Warn('LocalPostSolve','Do not know what to do with variable type: '//TRIM(I2S(pVar % TYPE)))
          END IF
        END DO
      END DO
    END DO   
!------------------------------------------------------------------------------
  END SUBROUTINE LocalPostSolve
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GlobalPostScale()
!------------------------------------------------------------------------------
   INTEGER :: dofs,i,j,Vari
   TYPE(Variable_t), POINTER :: pVar
   REAL(KIND=dp), ALLOCATABLE :: tmp(:)
   LOGICAL :: DoneWeight = .FALSE.
   REAL(KIND=dp) :: PotDiff, Capacitance, ControlTarget, ControlScaling, val
   
   VolTot = ParallelReduction(VolTot)
   EnergyTot = ParallelReduction(EnergyTot)
  
   WRITE( Message,'(A,ES12.3)') 'Total Electric Energy:', Energytot
   CALL Info( Caller, Message, Level=6 )
   CALL ListAddConstReal( Model % Simulation,'res: Electric Energy', Energytot )
   
   PotDiff = DirichletDofsRange( Solver )     
   IF( PotDiff > TINY( PotDiff ) ) THEN
     CALL ListAddConstReal( Model % Simulation,'res: Potential Difference', PotDiff )
     Capacitance = 2 * EnergyTot / PotDiff**2 
     WRITE( Message,'(A,ES12.3)') 'Effective Capacitance:', Capacitance
     CALL Info(Caller, Message, Level=6 )
     CALL ListAddConstReal( Model % Simulation,&
         'RES: Effective Capacitance', Capacitance )
   END IF
     
    
   DO Vari = 1, 8
     pVar => PostVars(Vari) % Var
     IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
     IF( PostVars(Vari) % NodalField ) CYCLE
     
     ! This is the only type of variable needing scaling!
     IF( pVar % TYPE /= variable_on_nodes ) CYCLE

     dofs = pVar % Dofs
     
     IF ( ParEnv % PEs > 1) THEN
       ! If we need to scale then also communicate the weight
       IF( .NOT. DoneWeight ) THEN
         CALL ParallelSumVector(Solver % Matrix, WeightVector )
         DoneWeight = .TRUE.
       END IF

       IF( dofs == 1 ) THEN
         CALL ParallelSumVector(Solver % Matrix, pVar % Values )
       ELSE
         IF(.NOT. ALLOCATED( tmp ) ) THEN
           ALLOCATE( tmp( SIZE( WeightVector ) ) )         
         END IF         
         DO i=1,dofs
           tmp = pVar % Values(i::dofs)
           CALL ParallelSumVector(Solver % Matrix, tmp)
           pVar % Values(i::dofs) = tmp
         END DO
       END IF
     END IF
     
     DO i=1,dofs
       WHERE( ABS( WeightVector ) > EPSILON( val ) ) &
           pVar % Values(i::dofs) = pVar % Values(i::dofs) / WeightVector
     END DO
   END DO

   ! Apply physical scaling in the end, if requested
   !------------------------------------------------------------------------
   ControlTarget = GetCReal( Params,'Energy Control',Found)
   IF( Found ) THEN
     CALL Info( Caller,'Scaling energy to desired value',Level=6)
     ControlScaling = SQRT( ControlTarget / EnergyTot )
   END IF

   IF( .NOT. Found ) THEN
     ControlTarget = GetCReal( Params,'Charge Control', Found ) 
     IF( Found ) THEN
       CALL Info( Caller,'Scaling charge to desired value',Level=6)      
       IF( PotDiff < TINY( PotDiff ) ) THEN
         CALL Fatal(Caller,'Charge cannot be controlled without pot. difference')
       END IF
       ControlScaling = ControlTarget / ( EnergyTot / PotDiff )
     END IF
   END IF
     
   IF( Found ) THEN
     WRITE( Message,'(A,ES12.3)') 'Control Scaling:', ControlScaling
     CALL Info(Caller, Message, Level=4 )
     CALL ListAddConstReal( Model % Simulation,'RES: StatElec Scaling',ControlScaling )

     Solver % Variable % Values = ControlScaling * Solver % Variable % Values
          
     DO Vari = 1, 8 
       pVar => PostVars(Vari) % Var
       IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
       IF( PostVars(Vari) % FieldType == 1 ) THEN
         ! Energy density scales quadratically
         pVar % Values = (ControlScaling**2) * pVar % Values
       ELSE
         ! other fields scave linearly         
         pVar % Values = ControlScaling * pVar % Values
       END IF
     END DO
   END IF

#if 0
   DO Vari = 1, 8 
     pVar => PostVars(Vari) % Var
     IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
     ! check this
     CALL InvalidateVariable( Model % Meshes, Solver % Mesh,PostVars(Vari) % Var % Name ) 
   END DO
#endif
   
!------------------------------------------------------------------------------
 END SUBROUTINE GlobalPostScale
!------------------------------------------------------------------------------

!------------------------------------------------------------------------
END SUBROUTINE StatElecSolver_Post
!------------------------------------------------------------------------
