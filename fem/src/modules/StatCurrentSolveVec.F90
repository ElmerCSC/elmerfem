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
! *  Module for static current conduction.
! *  Based on the multithreaded and vectorized ModelPDE by Mikko Byckling.
! *  Replaces gradually the old StatCurrentSolver that is not optimized.  
! *
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Created: 18.09.2018
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
!> Initialization of the primary solver, i.e. StatCurrentSolver.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'StatCurrentSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, CalculateElemental, CalculateNodal, PostActive 
  INTEGER :: dim  
   
  Params => GetSolverParams()
  dim = CoordinateSystemDimension()
  PostActive = .FALSE.
  
  CALL ListAddNewString( Params,'Variable','Potential')
  
  CalculateElemental = ListGetLogical( Params,'Calculate Elemental Fields',Found )
  CalculateNodal = ListGetLogical( Params,'Calculate Nodal Fields',Found )
  
  IF(.NOT. (CalculateElemental .OR. CalculateNodal ) ) THEN
    CalculateNodal = .TRUE.
  END IF

  IF (ListGetLogical(Params,'Calculate Joule Heating',Found)) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Joule Heating e' )
    IF( CalculateNodal ) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Joule Heating' )
    PostActive = .TRUE.
  END IF
  
  IF( ListGetLogical(Params,'Calculate Volume Current',Found) ) THEN
    IF( CalculateElemental ) & 
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        '-dg Volume Current e[Volume Current e:'//I2S(dim)//']' )
    IF( CalculateNodal ) &
        CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Volume Current[Volume Current:'//I2S(dim)//']' )       
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
  IF (ListGetLogical(Params,'Calculate Nodal Heating',Found))  THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params), &
        'Nodal Joule Heating' )
    PostActive = .TRUE.
  END IF
  IF( ListGetLogical(Params,'Calculate Nodal Current',Found) ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
        'Nodal Current[Nodal Current:'//I2S(dim)//']' )
    PostActive = .TRUE.
  END IF

  ! These use one flag to call library features to compute automatically
  ! a conductivity matrix.
  IF( ListGetLogical(Params,'Calculate Conductivity Matrix',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Constraint Modes Analysis',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Lumped',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Fluxes',.TRUE.)
    CALL ListAddNewLogical( Params,'Constraint Modes Matrix Symmetric',.TRUE.)
    CALL ListAddNewString( Params,'Constraint Modes Matrix Filename',&
        'ConductivityMatrix.dat',.FALSE.)
    CALL ListRenameAllBC( Model,'Conductivity Body','Constraint Mode Potential')
  END IF

  ! If no fields need to be computed do not even call the _post solver!
  CALL ListAddLogical(Params,'PostSolver Active',PostActive)

  ! If library adaptivity is compiled with, use that by default.
#ifdef LIBRARY_ADAPTIVIVTY
  CALL ListAddNewLogical(Params,'Library Adaptivity',.TRUE.)
#endif
  
END SUBROUTINE StatCurrentSolver_Init


!-----------------------------------------------------------------------------
!> A modern version for static current conduction supporting multithreading and
!> SIMD friendly ElmerSolver kernels. 
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Adaptive
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
  CHARACTER(*), PARAMETER :: Caller = 'StatCurrentSolver'
!------------------------------------------------------------------------------

  INTERFACE
    FUNCTION StatCurrentSolver_Boundary_Residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
      INTEGER :: Perm(:)
    END FUNCTION StatCurrentSolver_Boundary_Residual

    FUNCTION StatCurrentSolver_Edge_Residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Edge
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2)
      INTEGER :: Perm(:)
    END FUNCTION StatCurrentSolver_Edge_Residual

    FUNCTION StatCurrentSolver_Inside_Residual(Model, Element, Mesh, Quant, Perm, Fnorm) RESULT(Indicator)
      USE Types
      TYPE(Element_t), POINTER :: Element
      TYPE(Model_t) :: Model
      TYPE(Mesh_t), POINTER :: Mesh
      REAL(KIND=dp) :: Quant(:), Indicator(2), Fnorm
      INTEGER :: Perm(:)
    END FUNCTION StatCurrentSolver_Inside_Residual
  END INTERFACE

!------------------------------------------------------------------------------
  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving static current conduction solver')

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
      CALL Info( Caller,'Assembly of colour: '//I2S(col),Level=15)
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
      CALL Info(Caller,'Assembly of boundary colour: '//I2S(col),Level=10)
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

  IF (ListGetLogical(Params, 'Adaptive Mesh Refinement', Found)) THEN
    IF(.NOT. ListGetLogical(Params,'Library Adaptivity',Found)) THEN
      CALL RefineMesh(Model, Solver, Solver % Variable % Values, Solver % Variable % Perm, &
          StatCurrentSolver_Inside_Residual, StatCurrentSolver_Edge_Residual, &
          StatCurrentSolver_Boundary_Residual)
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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, POINTER  :: CondAtIpVec(:), EpsAtIpVec(:), SourceAtIpVec(:)
    REAL(KIND=dp) :: eps0, weight
    LOGICAL :: Stat,Found, Pref
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, EpsCoeff_h
    SAVE Eps0
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, Eps0, DetJVec, &
    !$OMP               MASS, STIFF, FORCE, Nodes, &
    !$OMP               SourceCoeff_h, CondCoeff_h, EpsCoeff_h, &
    !$OMP               SourceAtIpVec, CondAtIpVec, EpsAtIpVec )
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJVec
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Source')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
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
            DEALLOCATE(Basis, dBasisdx, DetJVec, MASS, STIFF, FORCE )
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJVec(ngp), &
          MASS(nd,nd), STIFF(nd,nd), FORCE(nd), STAT=allocstat)
      
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
    DO i=1,ngp
      DetJVec(i) = IP % s(i) * DetJVec(i)
    END DO

    ! electric conductivity term: STIFF=STIFF+(rho*grad(u),grad(v))
    CondAtIpVec => ListGetElementRealVec( CondCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJVec, STIFF, CondAtIpVec )
    END IF
    
    ! time derivative of potential: MASS=MASS+(eps*grad(u),grad(v))
    IF( Transient ) THEN
      EpsAtIpVec => ListGetElementRealVec( EpsCoeff_h, ngp, Basis, Element, Found ) 
      IF( Found ) THEN
        CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJVec, MASS, EpsAtIpVec )
        MASS(1:nd,1:nd) = Eps0 * MASS(1:nd,1:nd)
      END IF
    END IF
      
    ! source term: FORCE=FORCE+(u,f)
    SourceAtIpVec => ListGetElementRealVec( SourceCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJVec, SourceAtIpVec, FORCE)
    END IF
      
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: eps0, weight
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, CondAtIp, DetJ, A
    REAL(KIND=dp), POINTER :: CondTensor(:,:)
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,p,q,dim,m,allocstat,CondRank
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, EpsCoeff_h

    SAVE Eps0
!------------------------------------------------------------------------------

    !$OMP THREADPRIVATE(Basis, dBasisdx, Eps0, &
    !$OMP               MASS, STIFF, FORCE, Nodes, &
    !$OMP               SourceCoeff_h, CondCoeff_h, EpsCoeff_h )
    
    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Source')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity',UnfoundFatal=.TRUE.)
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
          MASS(m,m), STIFF(m,m), FORCE(m), STAT=allocstat)
      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL GetElementNodes( Nodes, UElement=Element )

    ! Initialize
    MASS  = 0._dp
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
      CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, Found, &
         GaussPoint = t, Rdim = CondRank, Rtensor = CondTensor ) 
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

      IF( Transient ) THEN
        EpsAtIp = Eps0 * ListGetElementReal( EpsCoeff_h, Basis, Element, Found )
        IF( Found ) THEN
          MASS(1:nd,1:nd) = MASS(1:nd,1:nd) + Weight * &
              EpsAtIp * MATMUL( dBasisdx(1:nd,:), TRANSPOSE( dBasisdx(1:nd,:) ) )
        END IF
      END IF

      SourceAtIP = ListGetElementReal( SourceCoeff_h, Basis, Element, Found ) 
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * SourceAtIP * Basis(1:nd)
      END IF
    END DO
    
    IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
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
    REAL(KIND=dp) :: F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),DetJ,Coord(3),Normal(3)
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found,RobinBC
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Robin_h, Ext_h, Farfield_h

    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h,Farfield_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Current Density')
      CALL ListInitElementKeyword( Robin_h,'Boundary Condition','External Conductivity')
      CALL ListInitElementKeyword( Ext_h,'Boundary Condition','External Potential')
      CALL ListInitElementKeyword( Farfield_h,'Boundary Condition','Farfield Potential')
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

      ! Robin condition (r*(u-u_0)):
      ! ---------------------------
      Ext = ListGetElementReal( Farfield_h, Basis, Element, RobinBC )
      IF( RobinBC ) THEN
        Coord(1) = SUM( Nodes % x(1:n)*Basis(1:n) )
        Coord(2) = SUM( Nodes % y(1:n)*Basis(1:n) )
        Coord(3) = SUM( Nodes % z(1:n)*Basis(1:n) )
        Normal = NormalVector( Element, Nodes, IP % u(t), IP % v(t), .TRUE. )
        C = SUM( Coord * Normal ) / SUM( Coord * Coord )         
      ELSE
        C = ListGetElementReal( Robin_h, Basis, Element, RobinBC )
        Ext = ListGetElementReal( Ext_h, Basis, Element, Found )
      END IF
        
      IF( RobinBC ) THEN
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * C * Ext * Basis(1:nd)
      END IF      
    END DO
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE StatCurrentSolver
!------------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!> A solver performing the postprocessing for the primary solver.
!> This solver is called at the DefaultFinish() slot of the primary solver.
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver_post( Model,Solver,dt,Transient )
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
  INTEGER :: i, dofs, n, nb, nd, t, active, CondRank
  LOGICAL :: Found, InitHandles
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp), ALLOCATABLE :: WeightVector(:),MASS(:,:),FORCE(:,:),&
      PotInteg(:),PotVol(:)
  INTEGER, POINTER :: WeightPerm(:)
  CHARACTER(*), PARAMETER :: Caller = 'StatCurrentSolver_post'
  LOGICAL :: CalcCurrent, CalcField, CalcHeating, NeedScaling, ConstantWeights, &
      Axisymmetric, CalcAvePotential
  TYPE(ValueList_t), POINTER :: Params
  REAL(KIND=dp) :: HeatingTot, Voltot
  REAL(KIND=dp), POINTER :: CondTensor(:,:)  
  
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
  PostVars(1) % Var => VariableGet( Mesh % Variables, 'Nodal Joule Heating')
  PostVars(1) % NodalField = .TRUE.
  PostVars(2) % Var => VariableGet( Mesh % Variables, 'Joule Heating')
  PostVars(3) % Var => VariableGet( Mesh % Variables, 'Joule Heating e')
  PostVars(1:3) % FieldType = 1

  ! Electric current: type 2, components 2:4
  PostVars(4) % Var => VariableGet( Mesh % Variables, 'Nodal Current')
  PostVars(4) % NodalField = .TRUE.
  PostVars(5) % Var => VariableGet( Mesh % Variables, 'Volume Current')
  PostVars(6) % Var => VariableGet( Mesh % Variables, 'Volume Current e')
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
  
  CalcHeating = ANY( PostVars(1:3) % HaveVar ) 
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
  HeatingTot = 0.0_dp
  VolTot = 0.0_dp

  !$OMP PARALLEL &
  !$OMP SHARED(Solver, Active) &
  !$OMP PRIVATE(t,Element, n, InitHandles, MASS, FORCE)
  
  !$OMP SINGLE
  Active = GetNOFActive(Solver)
  !$OMP END SINGLE

  !$OMP DO
  DO t = 1, Active
    Element => GetActiveElement(t)
    IF( ParEnv % PEs > 1 ) THEN
      IF( ParEnv % MyPe /= Element % PartIndex ) CYCLE
    END IF
    n  = GetElementNOFNodes(Element)
    CALL LocalPostAssembly( Element, n, InitHandles, MASS, FORCE )
    CALL LocalPostSolve( Element, n, MASS, FORCE )
  END DO
  !$OMP END DO 
  !$OMP END PARALLEL
  
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
        WRITE( Message,'(A,ES12.5)') 'Average body'//I2S(i)//' potential: ',PotAve
        CALL Info(Caller,Message,Level=7)
        CALL ListAddConstReal( Model % Simulation,&
            'res: Average body'//I2S(i)//' potential',PotAve)
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
    REAL(KIND=dp) :: SourceAtIp, EpsAtIp, CondAtIp, DetJ
    REAL(KIND=dp) :: Grad(3), CondGrad(3), Heat
    LOGICAL :: Stat,Found
    INTEGER :: i,j,t,p,q,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, EpsCoeff_h
    SAVE Eps0
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, Eps0, ElementPot, &
    !$OMP               Nodes,SourceCoeff_h, CondCoeff_h, EpsCoeff_h)

    
!------------------------------------------------------------------------------
    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Source')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
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

      CondAtIp = ListGetElementReal( CondCoeff_h, Basis, Element, Found, &
         GaussPoint = t, Rdim = CondRank, Rtensor = CondTensor ) 

      ! EpsAtIp = Eps0 * ListGetElementReal( EpsCoeff_h, Basis, Element, Found )
        
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

      ! Compute the volume current: J = cond (-grad Phi)
      !------------------------------------------------------------------------------
      IF( CondRank == 0 ) THEN
        CondGrad(1:dim) = CondAtIp * Grad(1:dim)
      ELSE IF( CondRank == 1 ) THEN
        CondGrad(1:dim) = CondTensor(1:dim,1) * Grad(1:dim)
      ELSE IF( CondRank == 2 ) THEN
        DO i = 1, DIM
          CondGrad(i) = SUM( CondTensor(i,1:dim) * Grad(1:dim) )
        END DO
      END IF

      IF( CalcCurrent ) THEN
        DO j=1,dim
          Force(2+j,1:n) = Force(2+j,1:n) - CondGrad(j) * Weight * Basis(1:n)
        END DO
      END IF

      ! Compute the Joule heating: H,tot = Integral (E . D)dV
      !------------------------------------------------------------------------------
      Heat = SUM( Grad(1:dim) * CondGrad(1:dim) )      
      IF( CalcHeating ) THEN
        Force(2,1:n) = Force(2,1:n) + Heat * Weight * Basis(1:n)
      END IF

      IF( CalcAvePotential ) THEN
        i = Element % BodyId
        PotVol(i) = PotVol(i) + Weight
        PotInteg(i) = PotInteg(i) + Weight * SUM( Basis(1:n) * ElementPot(1:n) ) 
      END IF
      
      VolTot = VolTot + Weight
      HeatingTot = HeatingTot + Weight * Heat 
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
   LOGICAL :: DoneWeight = .FALSE.
   REAL(KIND=dp) :: PotDiff, Resistance, ControlTarget, ControlScaling, val
   
   VolTot     = ParallelReduction(VolTot)
   HeatingTot = ParallelReduction(HeatingTot)

  
   WRITE( Message, * ) 'Total Heating Power   :', Heatingtot
   CALL Info( Caller, Message, Level=6 )
   CALL ListAddConstReal( Model % Simulation,'RES: Total Joule Heating', Heatingtot )
   
   PotDiff = DirichletDofsRange( Solver )     
   IF( PotDiff > TINY( PotDiff ) ) THEN
     Resistance = PotDiff**2 / HeatingTot
     WRITE( Message, * ) 'Effective Resistance  :', Resistance
     CALL Info(Caller, Message, Level=6 )
     CALL ListAddConstReal( Model % Simulation,'RES: Effective Resistance', Resistance )
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
   ControlTarget = GetCReal( Params,'Power Control',Found)
   IF( Found ) THEN
     CALL Info( Caller,'Scaling power to desired value',Level=6)
     ControlScaling = SQRT( ControlTarget / HeatingTot )
   END IF

   IF( .NOT. Found ) THEN
     ControlTarget = GetCReal( Params,'Current Control', Found ) 
     IF( Found ) THEN
       CALL Info( Caller,'Scaling current to desired value',Level=6)      
       IF( PotDiff < TINY( PotDiff ) ) THEN
         CALL Fatal(Caller,'Current cannot be controlled without pot. difference')
       END IF
       ControlScaling = ControlTarget / ( HeatingTot / PotDiff )
     END IF
   END IF
     
   IF( Found ) THEN
     WRITE( Message, * ) 'Control Scaling       :', ControlScaling
     CALL Info(Caller, Message, Level=4 )
     CALL ListAddConstReal( Model % Simulation, &
         'RES: CurrentSolver Scaling', ControlScaling )
     Solver % Variable % Values = ControlScaling * Solver % Variable % Values
          
     DO Vari = 1, 8 
       pVar => PostVars(Vari) % Var
       IF( .NOT. ASSOCIATED( pVar ) ) CYCLE
       IF( PostVars(Vari) % FieldType == 1 ) THEN
         ! Joule heating scales quadratically
         pVar % Values = (ControlScaling**2) * pVar % Values
       ELSE
         ! other fields save linearly         
         pVar % Values = ControlScaling * pVar % Values
       END IF
     END DO
   END IF
     
   
!------------------------------------------------------------------------------
 END SUBROUTINE GlobalPostScale
!------------------------------------------------------------------------------

!------------------------------------------------------------------------
END SUBROUTINE StatCurrentSolver_Post
!------------------------------------------------------------------------
    
  
  !------------------------------------------------------------------------------
  FUNCTION StatCurrentSolver_boundary_residual(Model, Edge, Mesh, Quant, Perm, Gnorm) RESULT(Indicator)
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
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), &
        EdgeBasis(:), Basis(:), x(:), y(:), z(:), &
        dBasisdx(:, :), Potential(:), Flux(:)
    REAL(KIND=dp) :: Normal(3), EdgeLength, gx, gy, gz, Conductivity
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
              dBasisdx(nd, 3), NodalConductivity(nd), Flux(nd), &
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
                            'Electric Conductivity', Hwrk, en, Edge % NodeIndexes, stat)
      IF (.NOT. stat) THEN
        CALL Fatal('StatCurrentSolver_boundary_residual:','Electric Conductivity not found')
      END IF
      NodalConductivity(1:en) = Hwrk(1, 1, 1:en)
  
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
  !          Conductivity at the integration point:
  !          --------------------------------------
        Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
  !
  !          given flux at integration point:
  !          --------------------------------
        Residual = -SUM(Flux(1:en) * EdgeBasis(1:en))
  
  !          flux given by the computed solution, and
  !          force norm for scaling the residual:
  !          -----------------------------------------
        IF (CurrentCoordinateSystem() == Cartesian) THEN
          DO k = 1, dim
            Residual = Residual + Conductivity * &
                        SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k)
  
            Gnorm = Gnorm + s * (Conductivity * &
                                  SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(k))**2
          END DO
        ELSE
          DO k = 1, dim
            DO l = 1, dim
              Residual = Residual + Metric(k, l) * Conductivity * &
                          SUM(dBasisdx(1:nd, k) * Potential(1:nd)) * Normal(l)
  
              Gnorm = Gnorm + s * (Metric(k, l) * Conductivity * &
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
    END FUNCTION StatCurrentSolver_boundary_residual
    !------------------------------------------------------------------------------
      
  
  FUNCTION StatCurrentSolver_edge_residual(Model, Edge, Mesh, Quant, Perm) RESULT(Indicator)
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
    REAL(KIND=dp) :: Conductivity
    REAL(KIND=dp) :: u, v, w, s, detJ
    REAL(KIND=dp) :: Grad(3, 3), Normal(3), EdgeLength, Jump
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:), x(:), y(:), z(:), EdgeBasis(:), &
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
  
    ALLOCATE (NodalConductivity(en), EdgeBasis(en), Basis(n), &
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
        IF (.NOT. ASSOCIATED(Element)) CYCLE
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
                              'Electric Conductivity', Hwrk, en, Edge % NodeIndexes, stat)
        IF (.NOT. stat) THEN
          CALL Fatal('StatCurrentSolver_edge_residual:', 'Electric Conductivity not found')
        END IF
  
        NodalConductivity(1:en) = Hwrk(1, 1, 1:en)
        Conductivity = SUM(NodalConductivity(1:en) * EdgeBasis(1:en))
  !
  !          Potential at element nodal points:
  !          ------------------------------------
        Potential(1:nd) = Quant(Perm(Indexes(1:nd)))
  !
  !          Finally, the flux:
  !          ------------------
        DO j = 1, dim
          Grad(j, i) = Conductivity * SUM(dBasisdx(1:nd, j) * Potential(1:nd))
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
    DEALLOCATE (x, y, z, NodalConductivity, EdgeBasis, &
                Basis, dBasisdx, Potential)
  !------------------------------------------------------------------------------
    END FUNCTION StatCurrentSolver_edge_residual
    !------------------------------------------------------------------------------
      
  !------------------------------------------------------------------------------
  FUNCTION StatCurrentSolver_inside_residual(Model, Element, Mesh, Quant, Perm, Fnorm) RESULT(Indicator)
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
    REAL(KIND=dp), ALLOCATABLE :: NodalConductivity(:)
    REAL(KIND=dp), ALLOCATABLE :: NodalSource(:), Potential(:), PrevPot(:)
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:, :), ddBasisddx(:, :, :)
    REAL(KIND=dp) :: u, v, w, s, detJ
    REAL(KIND=dp) :: Conductivity, dt
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
  
  !    Allocate local arrays
  !    ----------------------
    nd = GetElementNOFDOFs(Element)
    n = GetElementNOFNodes(Element)
    ALLOCATE (NodalConductivity(nd), &
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
  !    Material parameters: conductivity
  !    ---------------------------------
    k = ListGetInteger(Model % Bodies(Element % BodyId) % Values, 'Material', &
                        minv=1, maxv=Model % NumberOfMaterials)
  
    Material => Model % Materials(k) % Values
  
    CALL ListGetRealArray(Material, 'Electric Conductivity', Hwrk, n, Element % NodeIndexes, stat)
    IF (.NOT. stat) THEN
      CALL Fatal('StatCurrentSolver_inside_residual:', 'Electric Conductivity not found')
    END IF
    NodalConductivity(1:n) = Hwrk(1, 1, 1:n)
  
  !
  !    Current source density (source):
  !    --------------------------------
  !
    k = ListGetInteger( &
        Model % Bodies(Element % BodyId) % Values, 'Body Force', Found, &
        1, Model % NumberOFBodyForces)
  
    NodalSource = 0.0d0
    IF (Found .AND. k > 0) THEN
      NodalSource(1:n) = ListGetReal(Model % BodyForces(k) % Values, &
          'Current Source', n, Element % NodeIndexes, stat)
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
  
      Conductivity = SUM(NodalConductivity(1:n) * Basis(1:n))
  !
  !       Residual of the current conservation equation:
  !
  !        R = -div(σ grad(u)) - s
  !       ---------------------------------------------------
  !
  !       or more generally:
  !
  !        R = -g^{jk} (σ u_{,j}}_{,k}) - s
  !       ---------------------------------------------------
  !
      Residual = -SUM(NodalSource(1:n) * Basis(1:n))
  
      IF (CurrentCoordinateSystem() == Cartesian) THEN
        DO j = 1, dim
  !
  !             - grad(σ).grad(u):
  !             -------------------
          Residual = Residual - &
                      SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                      SUM(NodalConductivity(1:n) * dBasisdx(1:n, j))
  
  !
  !             - σ div(grad(u)):
  !             ------------------
          Residual = Residual - Conductivity * &
                      SUM(Potential(1:nd) * ddBasisddx(1:nd, j, j))
        END DO
      ELSE
        DO j = 1, dim
          DO k = 1, dim
  !
  !                - g^{jk} σ_{,k} u_{,j}:
  !                ------------------------
            Residual = Residual - Metric(j, k) * &
                        SUM(Potential(1:nd) * dBasisdx(1:nd, j)) * &
                        SUM(NodalConductivity(1:n) * dBasisdx(1:n, k))
  
  !
  !                - g^{jk} σ u_{,jk}:
  !                --------------------
            Residual = Residual - Metric(j, k) * Conductivity * &
                        SUM(Potential(1:nd) * ddBasisddx(1:nd, j, k))
  !
  !                + g^{jk} σ Γ_{jk}^l u_{,l}:
  !                ----------------------------
            DO l = 1, dim
              Residual = Residual + Metric(j, k) * Conductivity * &
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
  END FUNCTION StatCurrentSolver_inside_residual
        !------------------------------------------------------------------------------
        
  
  
