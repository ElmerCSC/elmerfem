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
! *  Based on the multithreaded and vectorized StatCurrentSolveVec that is based on
! *  ModelPDE by Mikko Byckling. Replaces gradually the old StatElecSolver that is
! *  not optimized.  
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
SUBROUTINE StatElecSolver_Init0( Model,Solver,dt,Transient )
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
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver_Init0'
!------------------------------------------------------------------------------
  IF (ListCheckPresentAnyMaterial(Model, 'Relative Permittivity Im')) THEN
    CALL Info(Caller, 'Complex-valued solution triggered by Relative Permittivity Im', Level=7)
    Params => GetSolverParams()
    CALL ListAddNewLogical(Params, 'Linear System Complex', .TRUE.)
    CALL ListAddNewInteger(Params, 'Variable DOFs', 2)
  END IF
!------------------------------------------------------------------------------
END SUBROUTINE StatElecSolver_Init0
!------------------------------------------------------------------------------
  
  
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
  LOGICAL :: Found, CalculateElemental, CalculateNodal, PostActive
  INTEGER :: dim
   
  Params => GetSolverParams()
  dim = CoordinateSystemDimension()

  IF (ListGetInteger(Params, 'Variable DOFs', Found) == 2) THEN
    CALL ListAddNewString( Params,'Variable','Potential[Potential re:1 Potential im:1]')
  ELSE
    CALL ListAddNewString( Params,'Variable','Potential')
  END IF
    
  PostActive = .FALSE.
  
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
    PostActive = .TRUE.
  END IF

  IF( ListGetLogical(Params,'Calculate Elecric Flux',Found) ) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Elecric Flux e[Elecric Flux e:'//I2S(dim)//']' )
    IF( CalculateNodal ) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Elecric Flux[Elecric Flux:'//I2S(dim)//']' )       
    PostActive = .TRUE.
  END IF

  IF( ListGetLogical(Params,'Calculate Electric Field',Found) ) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Electric Field e[Electric Field e:'//I2S(dim)//']' )
    IF( CalculateNodal ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Electric Field[Electric Field:'//I2S(dim)//']' )
    PostActive = .TRUE.
  END IF

  ! Nodal fields that may directly be associated as nodal loads
  IF (ListGetLogical(Params,'Calculate Nodal Energy',Found))  THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
        'Nodal Energy Density' )
    PostActive = .TRUE.
  END IF
  IF( ListGetLogical(Params,'Calculate Nodal Flux',Found) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Nodal Electric Flux[Nodal Electric Flux:'//I2S(dim)//']' )
    PostActive = .TRUE.
  END IF

  ! These use one flag to call library features to compute automatically
  ! a capacitance matrix.
  IF( ListGetLogical(Params,'Calculate Capacitance Matrix',Found ) ) THEN
    CALL Info('StatElecSolver_init','Using Constraint Modes functionality for Capacitance Matrix')
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Matrix Symmetric',.TRUE.)
    IF( ListCheckPresent( Params,'Capacitance Matrix Filename') ) THEN
      CALL ListRename( Params,'Capacitance Matrix Filename',&
          'Constraint Modes Matrix Filename', Found ) 
    ELSE     
      CALL ListAddNewString( Params,'Constraint Modes Matrix Filename',&
          'CapacitanceMatrix.dat',.FALSE.)
    END IF
    CALL ListRenameAllBC( Model,'Capacitance Body','Constraint Mode Potential')
    CALL ListRenameAllBodyForce( Model,'Capacitance Body','Constraint Mode Potential')
    CALL ListAddLogical( Params,'Optimize Bandwidth',.FALSE.)
    CALL Info('StatElecSolver_init','Suppressing bandwidth optimization in Capacitance Matrix computation!')
  END IF

  CALL ListAddInteger( Params,'Time Derivative Order', 0 )
  
  CALL ListWarnUnsupportedKeyword('solver','adaptive mesh redinement',FatalFound=.TRUE.)
  CALL ListWarnUnsupportedKeyword('body force','piezo material',FatalFound=.TRUE.)
  IF( ListCheckPresentAnyBC(Model,'infinity bc') ) THEN
    CALL Fatal('StatElecSolver_init','Use "Elecric Infinity BC" instead of "Infinity BC"')
  END IF
  
  ! If no fields need to be computed do not even call the _post solver!
  CALL ListAddLogical(Params,'PostSolver Active',PostActive)

  ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVIVTY
  CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
  
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
  LOGICAL :: Found, VecAsm, InitHandles, AxiSymmetric, CVersion
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver'

#ifndef LIBRARY_ADAPTIVITY
  INTERFACE
    FUNCTION StatElecSolver_Boundary_Residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END FUNCTION StatElecSolver_Boundary_Residual
    
    FUNCTION StatElecSolver_Edge_Residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END FUNCTION StatElecSolver_Edge_Residual
    
    FUNCTION StatElecSolver_Inside_Residual(Model, Element, Mesh, Quant, Perm, Fnorm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END FUNCTION StatElecSolver_Inside_Residual
  END INTERFACE
#endif
  
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving static electric field for insulators')

  Mesh => GetMesh()
  Params => GetSolverParams()

  CVersion = Solver % Variable % DOFs == 2
  IF (CVersion) CALL Info(Caller, &
      'Performing assembly for a complex-valued system',Level=7)
  
  IF( ListGetLogical( Params,'Follow P Curvature', Found )  ) THEN
    CALL FollowCurvedBoundary( Model, Mesh, .TRUE. ) 
  END IF
      
  CALL DefaultStart()

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

  IF (CVersion) VecAsm = .FALSE.
  
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
    !$OMP SHARED(Solver, Active, nColours, VecAsm, CVersion) &
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
          CALL LocalMatrix(  Element, n, nd+nb, nb, InitHandles, CVersion )
        END IF
      END DO
      !$OMP END DO
    END DO
    !$OMP END PARALLEL 

    CALL CheckTimer(Caller//'BulkAssembly',Delete=.TRUE.)
    totelem = 0

    CALL DefaultFinishBulkAssembly()

    ! If a complex-valued problem, the assembly for boundaries is not yet ready
    IF (CVersion) GOTO 201
    
    nColours = GetNOFBoundaryColours(Solver)

    CALL Info(Caller,'Performing boundary element assembly',Level=12)
    CALL ResetTimer(Caller//'BCAssembly')

    !$OMP PARALLEL &
    !$OMP SHARED(Active, Solver, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    InitHandles = .TRUE. 
    DO col=1,nColours
      !$OMP SINGLE
      CALL Info(Caller,'Assembly of boundary colour: '//I2S(col),Level=10)
      Active = GetNOFBoundaryActive(Solver)
      !$OMP END SINGLE

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
201 CONTINUE
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged == 1 ) EXIT

  END DO

  CALL DefaultFinish()

  IF (ListGetLogical(Solver % Values, 'Adaptive Mesh Refinement', Found)) THEN
    IF(.NOT. ListGetLogical( Solver % Values,'Library Adaptivity',Found )) THEN
      CALL RefineMesh(Model, Solver, Solver % Variable % Values, Solver % Variable % Perm, &
          StatElecSolver_Inside_Residual, StatElecSolver_Edge_Residual, &
          StatElecSolver_Boundary_Residual)
    END IF
  END IF

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
    REAL(KIND=dp) :: eps0
    LOGICAL :: Stat,Found
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

    IF( RelOrder /= 0 ) THEN
      IP = GaussPoints( Element, RelOrder = RelOrder)
    ELSE
      IP = GaussPoints( Element )
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
  SUBROUTINE LocalMatrix( Element, n, nd, nb, InitHandles, CVersion )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
    LOGICAL, INTENT(IN) :: CVersion 
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: STIFF(:,:), FORCE(:)
    COMPLEX(KIND=dp), ALLOCATABLE, SAVE :: CSTIFF(:,:), CFORCE(:)
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, DetJ
    COMPLEX(KIND=dp) :: CEpsAtIp
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Density')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity',InitIm=CVersion)
      
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
          STIFF(m,m), FORCE(m), CSTIFF(m,m), CFORCE(m), STAT=allocstat)
      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    STIFF = 0._dp
    FORCE = 0._dp
    IF (CVersion) THEN
      CSTIFF =  CMPLX(0.0_dp, 0.0_dp, kind=dp)
      CFORCE =  CMPLX(0.0_dp, 0.0_dp, kind=dp)
    END IF
    
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
      IF (CVersion) THEN
        CEpsAtIp = ListGetElementComplex(EpsCoeff_h, Basis, Element, Found, GaussPoint = t)
        CSTIFF(1:nd,1:nd) = CSTIFF(1:nd,1:nd) + Weight * &
            Eps0 * CEpsAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )
      ELSE
        EpsAtIp = ListGetElementReal( EpsCoeff_h, Basis, Element, Found, GaussPoint = t )      
        STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
            Eps0 * EpsAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )
      END IF

      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        IF (CVersion) THEN
          CFORCE(1:nd) = CFORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
        ELSE
          FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
        END IF
      END IF
    END DO

    IF (CVersion) THEN
      CALL CondensateP( nd-nb, nb, CStiff, CForce )
      CALL DefaultUpdateEquations(CStiff, CForce)
    ELSE
      CALL CondensateP( nd-nb, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
    END IF
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
    REAL(KIND=dp) :: Weight,Eps0,Alpha,Beta,Ext
    REAL(KIND=dp) :: Basis(nd),DetJ,Coord(3),Normal(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,GotSome
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Farfield_h, Infty_h, &
        LayerEps_h, LayerH_h, LayerRho_h, LayerV_h
    REAL(KIND=dp) :: LayerEps, LayerV, LayerRho, LayerH    
    SAVE Nodes, Eps0
    !$OMP THREADPRIVATE(Nodes,Flux_h,Farfield_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Electric Flux')
      CALL ListInitElementKeyword( Infty_h,'Boundary Condition','Electric Infinity BC')
      CALL ListInitElementKeyword( Farfield_h,'Boundary Condition','Farfield Potential')
      CALL ListInitElementKeyword( LayerEps_h,'Boundary Condition','Layer Relative Permittivity')
      CALL ListInitElementKeyword( LayerH_h,'Boundary Condition','Layer Thickness')
      CALL ListInitElementKeyword( LayerRho_h,'Boundary Condition','Layer Charge Density')
      CALL ListInitElementKeyword( LayerV_h,'Boundary Condition','Electrode Potential')
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
      ! BC: -epsilon@Phi/@n = -alpha Phi + beta
      !------------------------------------------

      ! Given flux:
      ! -----------
      Alpha = 0.0_dp
      Beta = ListGetElementReal( Flux_h, Basis, Element, GotSome )
      
      IF( ListGetElementLogical( Infty_h, Element, Found ) ) THEN
        ! Robin type of condition for farfield potential (r*(u-u_0)):
        ! -----------------------------------------------------------
        GotSome = .TRUE.
        Coord(1) = SUM( Nodes % x(1:n)*Basis(1:n) )
        Coord(2) = SUM( Nodes % y(1:n)*Basis(1:n) )
        Coord(3) = SUM( Nodes % z(1:n)*Basis(1:n) )
        
        Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t), .TRUE. )
        Ext = ListGetElementReal( Farfield_h, Basis, Element, Found )

        Alpha = Eps0 * SUM( Coord * Normal ) / SUM( Coord * Coord )         
        Beta = Beta + Alpha * Ext
      ELSE
        ! Boundary condition for electrostatic layer on the boundary.
        !------------------------------------------------------------
        LayerEps = ListGetElementReal( LayerEps_h, Basis, Element, Found )
        IF(Found) THEN
          GotSome = .TRUE.
          LayerH = ListGetElementReal( LayerH_h, Basis, Element, Found )
          IF ( .NOT. Found) THEN
            CALL Fatal( Caller,'Charge > Layer thickness < not given!' )
          END IF 
          LayerV = ListGetElementReal( LayerV_h, Basis, Element, Found )
          LayerRho = ListGetElementReal( LayerRho_h, Basis, Element, Found )

          Alpha = LayerEps / LayerH          
          Beta = Beta + Alpha * LayerV + 0.5_dp * LayerRho * LayerH / Eps0
        END IF
      END IF

      IF( GotSome ) THEN      
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * Alpha * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * Beta * Basis(1:nd)
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
  INTEGER :: i, n, nd, t
  LOGICAL :: Found, InitHandles = .TRUE.
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: WeightVector(:),FORCE(:,:),MASS(:,:),&
      PotInteg(:),PotVol(:)
  INTEGER, POINTER :: WeightPerm(:)
  CHARACTER(*), PARAMETER :: Caller = 'StatElecSolver_post'
  LOGICAL :: CalcCurrent, CalcField, CalcElectricDisp, NeedScaling, ConstantWeights, &
      Axisymmetric, CalcAvePotential, DoAve
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
  DoAve = ListGetLogical( Params,'Average Within Materials',Found ) 
  
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
  CALL Info(Caller,'Number of '//I2S(n)//' postprocessing fields',Level=8)
    
  ! Only create the nodal weights if we need to scale some nodal field
  NeedScaling = .FALSE.
  DO i=1,8
    IF( .NOT. PostVars(i) % HaveVar ) CYCLE
    IF( PostVars(i) % NodalField ) CYCLE
    IF( PostVars(i) % Var % TYPE == Variable_on_nodes ) THEN
      CALL Info(Caller,'Creating a weighting for scaling purposes from '//I2S(i),Level=10)
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
    nd = GetElementNOFDOFs(Element)
    CALL LocalPostAssembly( Element, n, nd, InitHandles, MASS, FORCE )
    CALL LocalPostSolve( Element, n, MASS, FORCE )
  END DO

  
  IF( NeedScaling ) THEN
    CALL Info(Caller,'Scaling the field values with weights',Level=12)
    CALL GlobalPostScale()
  END IF

  IF( DoAve ) THEN
    CALl GlobalPostAve()
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
        WRITE( Message,'(A,ES12.5)') 'Average body'//I2S(i)//' potential: ',PotAve
        CALL Info(Caller,Message,Level=7)
        CALL ListAddConstReal( Model % Simulation,&
            'res: Average body'//I2S(i)//' potential',PotAve)
      END DO
    END BLOCK
  END IF
    
  CALL Info(Caller,'All done',Level=12)

  

CONTAINS
   
  SUBROUTINE LocalPostAssembly( Element, n, nd, InitHandles, MASS, FORCE )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd
    TYPE(Element_t), POINTER :: Element
    LOGICAL, INTENT(INOUT) :: InitHandles
    REAL(KIND=dp) :: MASS(:,:), FORCE(:,:)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:),ElementPot(:)
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: EpsAtIp, DetJ
    REAL(KIND=dp) :: EpsGrad(3), Grad(3), Heat
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, EpsCoeff_h
    SAVE Eps0
    
!------------------------------------------------------------------------------
    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Charge Density')
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
      Basis = 0.0_dp; dBasisdx = 0.0_dp; ElementPot = 0.0_dp
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element, USolver=Solver )
    CALL GetScalarLocalSolution( ElementPot, UElement=Element, USolver=Solver) 
    
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
        Grad(j) = SUM( dBasisdx(1:nd,j) * ElementPot(1:nd) )
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
    LOGICAL :: LocalSolved, Erroneous
!------------------------------------------------------------------------------
    
    CALL LUdecomp(A,n,pivot,Erroneous)
    IF (Erroneous) CALL Fatal('LocalPostSolve', 'LU-decomposition fails')

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
            CALL LUSolve(n,A,x,pivot)
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
            CALL Warn('LocalPostSolve','Do not know what to do with variable type: '//I2S(pVar % TYPE))
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
         
   ! If we need to scale then also communicate the weight
   IF (ParEnv % PEs>1 .AND. NeedScaling) THEN
     CALL ParallelSumVector(Solver % Matrix, WeightVector )
   END IF

   DO Vari = 1, 8
     pVar => PostVars(Vari) % Var
     IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
     IF( PostVars(Vari) % NodalField ) CYCLE
     
     ! This is the only type of variable needing scaling!
     IF( pVar % TYPE /= variable_on_nodes ) CYCLE

     dofs = pVar % Dofs
     
     IF ( ParEnv % PEs > 1) THEN
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
       WHERE( ABS( WeightVector ) > TINY( val ) ) &
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

!------------------------------------------------------------------------------
  SUBROUTINE GlobalPostAve()
!------------------------------------------------------------------------------
    INTEGER :: i, Vari
    TYPE(Variable_t), POINTER :: pVar

    DO Vari = 1, 8
      pVar => PostVars(Vari) % Var
      IF( .NOT. ASSOCIATED( pVar ) ) CYCLE

      ! This is the only type of variable needing averaging!
      IF( pVar % TYPE == variable_on_nodes_on_elements ) THEN
        ! We reset the flag since, otherwise this will not be averaged...
        pVar % DgAveraged = .FALSE.
        CALL Info(Caller,'Averaging for field: '//TRIM(pVar % Name),Level=10)
        CALL CalculateBodyAverage(Mesh, pVar, .FALSE.)
      END IF
    END DO

  END SUBROUTINE GlobalPostAve
!------------------------------------------------------------------------------
END SUBROUTINE StatElecSolver_Post
!------------------------------------------------------------------------



!------------------------------------------------------------------------------
! Subroutine for computing residuals for adaptive mesh refinement.
!------------------------------------------------------------------------------
FUNCTION StatElecSolver_boundary_Residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
!------------------------------------------------------------------------------
  USE DefUtils
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
  INTEGER :: i, j, k, n, l, t, dim, Pn, En, nd
  LOGICAL :: stat, Found
  INTEGER, ALLOCATABLE :: Indexes(:)
  REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
  REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
  REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:), &
      EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
      dBasisdx(:, :), Potential(:), Flux(:)  
  REAL(KIND=dp) :: Normal(3), EdgeLength, gx, gy, gz, Permittivity
  REAL(KIND=dp) :: u, v, w, s, detJ
  REAL(KIND=dp) :: Residual, ResidualNorm
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

  LOGICAL :: First = .TRUE., Dirichlet
  SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
  IF (First) THEN
    First = .FALSE.
    NULLIFY (Hwrk)
  END IF

  Indicator = 0.0d0
  Gnorm = 0.0d0

  Metric = 0.0d0
  DO i = 1, 3
    Metric(i, i) = 1.0d0
  END DO

  SELECT CASE (CurrentCoordinateSystem())
  CASE (AxisSymmetric, CylindricSymmetric)
    dim = 3
  CASE DEFAULT
    dim = CoordinateSystemDimension()
  END SELECT
!
!    ---------------------------------------------

  Element => Edge % BoundaryInfo % Left

  IF (.NOT. ASSOCIATED(Element)) THEN
    Element => Edge % BoundaryInfo % Right
  ELSE IF (ANY(Perm(Element % NodeIndexes) <= 0)) THEN
    Element => Edge % BoundaryInfo % Right
  END IF

  IF (.NOT. ASSOCIATED(Element)) RETURN
  IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN

  en = Edge % TYPE % NumberOfNodes
  pn = Element % TYPE % NumberOfNodes

  ALLOCATE (EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en))

  EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
  EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
  EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

  nd = GetElementNOFDOFs(Element)
  ALLOCATE (Potential(nd), Basis(nd), &
            x(en), y(en), z(en), EdgeBasis(nd), &
            dBasisdx(nd, 3), NodalPermittivity(nd), Flux(nd), &
            Indexes(nd))

  nd = GetElementDOFs(Indexes, Element)

  ALLOCATE (Nodes % x(nd), Nodes % y(nd), Nodes % z(nd))
  Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
  Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
  Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))

  DO l = 1, en
    DO k = 1, pn
      IF (Edge % NodeIndexes(l) == Element % NodeIndexes(k)) THEN
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

  Indicator = 0.0d0
  EdgeLength = 0.0d0
  ResidualNorm = 0.0d0

  DO j = 1, Model % NumberOfBCs
    IF (Edge % BoundaryInfo % Constraint /= Model % BCs(j) % Tag) CYCLE

!
!       Check if dirichlet BC given:
!       ----------------------------
    Dirichlet = ListCheckPresent(Model % BCs(j) % Values, &
                                 ComponentName(Model % Solver % Variable))
    IF (.NOT. Dirichlet) THEN
      Dirichlet = ListCheckPrefix(Model % BCs(j) % Values, &
                                  'Constraint Mode')
    END IF
    ! TODO s = ListGetConstReal( Model % BCs(j) % Values,'Potential',Dirichlet )

!       Get various flux bc options:
!       ----------------------------

!       ...given flux:
!       --------------
    Flux(1:en) = ListGetReal(Model % BCs(j) % Values, &
                             'Electric Flux', en, Edge % NodeIndexes, Found)

!       get material parameters:
!       ------------------------
    k = ListGetInteger(Model % Bodies(Element % BodyId) % Values, 'Material', &
                       minv=1, maxv=Model % NumberOFMaterials)

    CALL ListGetRealArray(Model % Materials(k) % Values, &
                          'Relative Permittivity', Hwrk, en, Edge % NodeIndexes, stat)
    IF (.NOT. stat) &
      CALL ListGetRealArray(Model % Materials(k) % Values, &
                            'Permittivity', Hwrk, En, Edge % NodeIndexes)
    NodalPermittivity(1:en) = Hwrk(1, 1, 1:en)

!       elementwise nodal solution:
!       ---------------------------
    nd = GetElementDOFs(Indexes, Element)
    Potential(1:nd) = Quant(Perm(Indexes(1:nd)))

!       do the integration:
!       -------------------
    EdgeLength = 0.0d0
    ResidualNorm = 0.0d0

    IntegStuff = GaussPoints(Edge)

    DO t = 1, IntegStuff % n
      u = IntegStuff % u(t)
      v = IntegStuff % v(t)
      w = IntegStuff % w(t)

      stat = ElementInfo(Edge, EdgeNodes, u, v, w, detJ, &
                         EdgeBasis, dBasisdx)
      Normal = NormalVector(Edge, EdgeNodes, u, v, .TRUE.)

      IF (CurrentCoordinateSystem() == Cartesian) THEN
        s = IntegStuff % s(t) * detJ
      ELSE
        gx = SUM(EdgeBasis(1:en) * EdgeNodes % x(1:en))
        gy = SUM(EdgeBasis(1:en) * EdgeNodes % y(1:en))
        gz = SUM(EdgeBasis(1:en) * EdgeNodes % z(1:en))
        CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                  Symb, dSymb, gx, gy, gz)
        s = IntegStuff % s(t) * detJ * SqrtMetric
      END IF

!
!          Integration point in parent element local
!          coordinates:
!          -----------------------------------------
      u = SUM(EdgeBasis(1:en) * x(1:en))
      v = SUM(EdgeBasis(1:en) * y(1:en))
      w = SUM(EdgeBasis(1:en) * z(1:en))
      stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)

!
!          Electric permittivity at the integration point:
!          --------------------------------------------
      Permittivity = SUM(NodalPermittivity(1:en) * EdgeBasis(1:en))
!
!          given flux at integration point:
!          --------------------------------
      Residual = -SUM(Flux(1:en) * EdgeBasis(1:en))

!          flux given by the computed solution, and
!          force norm for scaling the residual:
!          -----------------------------------------
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        DO k = 1, dim
          Residual = Residual + Permittivity * &
                     SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k)

          Gnorm = Gnorm + s * (Permittivity * &
                               SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k))**2
        END DO
      ELSE
        DO k = 1, dim
          DO l = 1, dim
            Residual = Residual + Metric(k, l) * Permittivity * &
                       SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(l)

            Gnorm = Gnorm + s * (Metric(k, l) * Permittivity * &
                                 SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(l))**2
          END DO
        END DO
      END IF

      EdgeLength = EdgeLength + s
      IF (.NOT. Dirichlet) THEN
        ResidualNorm = ResidualNorm + s * Residual**2
      END IF
    END DO
    EXIT
  END DO

  IF (CoordinateSystemDimension() == 3) EdgeLength = SQRT(EdgeLength)

!    Gnorm = EdgeLength * Gnorm
  Indicator = EdgeLength * ResidualNorm
!------------------------------------------------------------------------------
END FUNCTION StatElecSolver_boundary_residual
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION StatElecSolver_edge_residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
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
  INTEGER :: i, j, k, l, n, t, dim, En, Pn, nd
  INTEGER, ALLOCATABLE :: Indexes(:)
  LOGICAL :: stat
  REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
  REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
  REAL(KIND=dp) :: Permittivity
  REAL(KIND=dp) :: u, v, w, s, detJ
  REAL(KIND=dp) :: Grad(3, 3), Normal(3), EdgeLength, Jump
  REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:), x(:), y(:), z(:), EdgeBasis(:), &
      Basis(:), dBasisdx(:, :), Potential(:)
  REAL(KIND=dp) :: ResidualNorm
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
  TYPE(ValueList_t), POINTER :: Material

  LOGICAL :: First = .TRUE.
  SAVE Hwrk, First
!------------------------------------------------------------------------------

  !    Initialize:
  !    -----------
  IF (First) THEN
    First = .FALSE.
    NULLIFY (Hwrk)
  END IF

  SELECT CASE (CurrentCoordinateSystem())
  CASE (AxisSymmetric, CylindricSymmetric)
    dim = 3
  CASE DEFAULT
    dim = CoordinateSystemDimension()
  END SELECT

  Metric = 0.0d0
  DO i = 1, 3
    Metric(i, i) = 1.0d0
  END DO
  Grad = 0.0d0
!
!    ---------------------------------------------

  n = Mesh % MaxElementDOFs
  ALLOCATE (Nodes % x(n), Nodes % y(n), Nodes % z(n))

  en = Edge % TYPE % NumberOfNodes
  ALLOCATE (EdgeNodes % x(en), EdgeNodes % y(en), EdgeNodes % z(en))

  EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
  EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
  EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)

  ALLOCATE (NodalPermittivity(en), EdgeBasis(en), Basis(n), &
            dBasisdx(n, 3), x(en), y(en), z(en), Potential(n), Indexes(n))

!    Integrate square of jump over edge:
!    -----------------------------------
  ResidualNorm = 0.0d0
  EdgeLength = 0.0d0
  Indicator = 0.0d0

  IntegStuff = GaussPoints(Edge)

  DO t = 1, IntegStuff % n

    u = IntegStuff % u(t)
    v = IntegStuff % v(t)
    w = IntegStuff % w(t)

    stat = ElementInfo(Edge, EdgeNodes, u, v, w, detJ, &
                       EdgeBasis, dBasisdx)

    Normal = NormalVector(Edge, EdgeNodes, u, v, .FALSE.)

    IF (CurrentCoordinateSystem() == Cartesian) THEN
      s = IntegStuff % s(t) * detJ
    ELSE
      u = SUM(EdgeBasis(1:en) * EdgeNodes % x(1:en))
      v = SUM(EdgeBasis(1:en) * EdgeNodes % y(1:en))
      w = SUM(EdgeBasis(1:en) * EdgeNodes % z(1:en))

      CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                Symb, dSymb, u, v, w)
      s = IntegStuff % s(t) * detJ * SqrtMetric
    END IF

    !
    ! Compute flux over the edge as seen by elements
    ! on both sides of the edge:
    ! ----------------------------------------------
    DO i = 1, 2
      SELECT CASE (i)
      CASE (1)
        Element => Edge % BoundaryInfo % Left
      CASE (2)
        Element => Edge % BoundaryInfo % Right
      END SELECT
!
!          Can this really happen (maybe it can...)  ?
!          -------------------------------------------
      IF (ANY(Perm(Element % NodeIndexes) <= 0)) CYCLE
!
!          Next, get the integration point in parent
!          local coordinates:
!          -----------------------------------------
      pn = Element % TYPE % NumberOfNodes

      DO j = 1, en
        DO k = 1, pn
          IF (Edge % NodeIndexes(j) == Element % NodeIndexes(k)) THEN
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
      nd = GetElementDOFs(Indexes, Element)
      Nodes % x(1:nd) = Mesh % Nodes % x(Indexes(1:nd))
      Nodes % y(1:nd) = Mesh % Nodes % y(Indexes(1:nd))
      Nodes % z(1:nd) = Mesh % Nodes % z(Indexes(1:nd))

      stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis, dBasisdx)
!
!          Material parameters:
!          --------------------
      k = ListGetInteger(Model % Bodies( &
                         Element % BodyId) % Values, 'Material', &
                         minv=1, maxv=Model % NumberOFMaterials)

      Material => Model % Materials(k) % Values
      CALL ListGetRealArray(Material, &
                            'Relative Permittivity', Hwrk, en, Edge % NodeIndexes, stat)  ! Should we have stat here?
      IF (.NOT. stat) &
        CALL ListGetRealArray(Material, &
                              'Permittivity', Hwrk, En, Edge % NodeIndexes)

      NodalPermittivity(1:en) = Hwrk(1, 1, 1:en)
      Permittivity = SUM(NodalPermittivity(1:en) * EdgeBasis(1:en))
!
!          Potential at element nodal points:
!          ------------------------------------
      Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
!
!          Finally, the flux:
!          ------------------
      DO j = 1, dim
        Grad(j, i) = Permittivity * SUM(dBasisdx(1:nd, j) * Potential(1:nd))
      END DO
    END DO

!       Compute square of the flux jump:
!       -------------------------------
    EdgeLength = EdgeLength + s
    Jump = 0.0d0
    DO k = 1, dim
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        Jump = Jump + (Grad(k, 1) - Grad(k, 2)) * Normal(k)
      ELSE
        DO l = 1, dim
          Jump = Jump + &
                 Metric(k, l) * (Grad(k, 1) - Grad(k, 2)) * Normal(l)
        END DO
      END IF
    END DO
    ResidualNorm = ResidualNorm + s * Jump**2
  END DO

  IF (dim == 3) EdgeLength = SQRT(EdgeLength)
  Indicator = EdgeLength * ResidualNorm

  DEALLOCATE (Nodes % x, Nodes % y, Nodes % z)

  DEALLOCATE (EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
  DEALLOCATE (x, y, z, NodalPermittivity, EdgeBasis, &
              Basis, dBasisdx, Potential)
!------------------------------------------------------------------------------
END FUNCTION StatElecSolver_Edge_Residual
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION StatElecSolver_Inside_residual(Model, Element, Mesh, &
                                Quant, Perm, Fnorm) RESULT(Indicator)
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
  INTEGER :: i, j, k, l, n, t, dim, nd
  INTEGER, ALLOCATABLE :: Indexes(:)
  LOGICAL :: stat, Found
  TYPE(Variable_t), POINTER :: Var
  REAL(KIND=dp), POINTER :: Hwrk(:, :, :)
  REAL(KIND=dp) :: SqrtMetric, Metric(3, 3), Symb(3, 3, 3), dSymb(3, 3, 3, 3)
  REAL(KIND=dp), ALLOCATABLE :: NodalPermittivity(:)
  REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Potential(:), PrevPot(:)
  REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:, :), ddBasisddx(:, :, :)
  REAL(KIND=dp) :: u, v, w, s, detJ
  REAL(KIND=dp) :: Permittivity, dt 
  REAL(KIND=dp) :: Residual, ResidualNorm, Area
  TYPE(ValueList_t), POINTER :: Material
  TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

  LOGICAL :: First = .TRUE.
  SAVE Hwrk, First
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
  Indicator = 0.0d0
  Fnorm = 0.0d0
!
!    Check if this eq. computed in this element:
!    -------------------------------------------
  IF (ANY(Perm(Element % NodeIndexes) <= 0)) RETURN

  IF (First) THEN
    First = .FALSE.
    NULLIFY (Hwrk)
  END IF

  Metric = 0.0d0
  DO i = 1, 3
    Metric(i, i) = 1.0d0
  END DO

  SELECT CASE (CurrentCoordinateSystem())
  CASE (AxisSymmetric, CylindricSymmetric)
    dim = 3
  CASE DEFAULT
    dim = CoordinateSystemDimension()
  END SELECT

!    Alllocate local arrays
!    ----------------------
  nd = GetElementNOFDOFs(Element)
  n = GetElementNOFNodes(Element)
  ALLOCATE (NodalPermittivity(nd), &
            PrevPot(nd), NodalSource(nd), Potential(nd), &
            Basis(nd), dBasisdx(nd, 3), ddBasisddx(nd, 3, 3), Indexes(nd))
!
!    Element nodal points:
!    ---------------------
  ALLOCATE (Nodes % x(nd), Nodes % y(nd), Nodes % z(nd))

  nd = GetElementDOFs(Indexes, Element)
  Nodes % x = Mesh % Nodes % x(Indexes(1:nd))
  Nodes % y = Mesh % Nodes % y(Indexes(1:nd))
  Nodes % z = Mesh % Nodes % z(Indexes(1:nd))
!
!    Elementwise nodal solution:
!    ---------------------------
  Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
!
!    Check for time dep.
!    -------------------
  PrevPot(1:nd) = Potential(1:nd)
  dt = Model % Solver % dt
  IF (ListGetString(Model % Simulation, 'Simulation Type') == 'transient') THEN
    Var => VariableGet(Model % Variables, 'Potential', .TRUE.)
    PrevPot(1:nd) = Var % PrevValues(Var % Perm(Indexes(1:nd)), 1)
  END IF
!
!    Material parameters: relative permittivity
!    ------------------------------------------
  k = ListGetInteger(Model % Bodies(Element % BodyId) % Values, 'Material', &
                     minv=1, maxv=Model % NumberOfMaterials)

  Material => Model % Materials(k) % Values

  CALL ListGetRealArray(Material, 'Relative Permittivity', Hwrk, n, Element % NodeIndexes, stat)
  IF (.NOT. stat) &
    CALL ListGetRealArray(Material, 'Permittivity', Hwrk, n, Element % NodeIndexes)
  NodalPermittivity(1:n) = Hwrk(1, 1, 1:n)

!
!    Charge density (source):
!    ------------------------
!
  k = ListGetInteger( &
      Model % Bodies(Element % BodyId) % Values, 'Body Force', Found, &
      1, Model % NumberOFBodyForces)

  NodalSource = 0.0d0
  IF (Found .AND. k > 0) THEN
    NodalSource(1:n) = ListGetReal(Model % BodyForces(k) % Values, &
                                   'Charge Density', n, Element % NodeIndexes, stat)
    IF (.NOT. stat) &
      NodalSource(1:n) = ListGetReal(Model % BodyForces(k) % Values, &
                                     'Source', n, Element % NodeIndexes)
  END IF

!
!    Integrate square of residual over element:
!    ------------------------------------------

  ResidualNorm = 0.0d0
  Area = 0.0d0

  IntegStuff = GaussPoints(Element)
  ddBasisddx = 0

  DO t = 1, IntegStuff % n
    u = IntegStuff % u(t)
    v = IntegStuff % v(t)
    w = IntegStuff % w(t)

    stat = ElementInfo(Element, Nodes, u, v, w, detJ, &
                       Basis, dBasisdx, ddBasisddx, .TRUE., .FALSE.)

    IF (CurrentCoordinateSystem() == Cartesian) THEN
      s = IntegStuff % s(t) * detJ
    ELSE
      u = SUM(Basis(1:nd) * Nodes % x(1:nd))
      v = SUM(Basis(1:nd) * Nodes % y(1:nd))
      w = SUM(Basis(1:nd) * Nodes % z(1:nd))

      CALL CoordinateSystemInfo(Metric, SqrtMetric, &
                                Symb, dSymb, u, v, w)
      s = IntegStuff % s(t) * detJ * SqrtMetric
    END IF

    Permittivity = SUM(NodalPermittivity(1:n) * Basis(1:n))
!
!       Residual of the electrostatic equation:
!
!        R = -div(e grad(u)) - s
!       ---------------------------------------------------
!
!       or more generally:
!
!        R = -g^{jk} (C T_{,j}}_{,k} - s
!       ---------------------------------------------------
!
    Residual = -SUM(NodalSource(1:n) * Basis(1:n))

    IF (CurrentCoordinateSystem() == Cartesian) THEN
      DO j = 1, dim
!
!             - grad(e).grad(T):
!             --------------------
!
        Residual = Residual - &
                   SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                   SUM(NodalPermittivity(1:n) * dBasisdx(1:n, j))

!
!             - e div(grad(u)):
!             -------------------
!
        Residual = Residual - Permittivity * &
                   SUM(Potential(1:nd) * ddBasisddx(1:nd, j, j))
      END DO
    ELSE
      DO j = 1, dim
        DO k = 1, dim
!
!                - g^{jk} C_{,k}T_{j}:
!                ---------------------
!
          Residual = Residual - Metric(j, k) * &
                     SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                     SUM(NodalPermittivity(1:n) * dBasisdx(1:n, k))

!
!                - g^{jk} C T_{,jk}:
!                -------------------
!
          Residual = Residual - Metric(j, k) * Permittivity * &
                     SUM(Potential(1:nd) * ddBasisddx(1:nd, j, k))
!
!                + g^{jk} C {_jk^l} T_{,l}:
!                ---------------------------
          DO l = 1, dim
            Residual = Residual + Metric(j, k) * Permittivity * &
                       Symb(j, k, l) * SUM(Potential(1:nd) * dBasisdx(1:nd, l))
          END DO
        END DO
      END DO
    END IF

!
!       Compute also force norm for scaling the residual:
!       -------------------------------------------------
    DO i = 1, dim
      Fnorm = Fnorm + s * (SUM(NodalSource(1:n) * Basis(1:n)))**2
    END DO
    Area = Area + s
    ResidualNorm = ResidualNorm + s * Residual**2
  END DO

!    Fnorm = Element % hk**2 * Fnorm
  Indicator = Element % hK**2 * ResidualNorm
!------------------------------------------------------------------------------
END FUNCTION StatElecSolver_inside_residual
!------------------------------------------------------------------------------
