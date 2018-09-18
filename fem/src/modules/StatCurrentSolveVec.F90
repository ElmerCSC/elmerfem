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
SUBROUTINE StatCurrentSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'StatCurrentSolver_init'
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found
  INTEGER :: dim
  
  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
    CALL Fatal(Caller,'Implemented only in cartesian coordinates')
  END IF

  IF (ListGetLogical(Params,'Calculate Joule Heating',Found)) &
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
      'Joule Heating' )
  
  IF (ListGetLogical(Params,'Calculate Nodal Heating',Found)) &
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
      'Nodal Joule Heating' )
  
  IF( ListGetLogical(Params,'Calculate Volume Current',Found) ) THEN
    IF( Dim == 2 ) THEN
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          'Volume Current[Volume Current:2]' )
    ELSE
      CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params), &
          'Volume Current[Volume Current:3]' )
    END IF
  END IF

  CALL ListUntreatedFatal( Params,'Power Control',Caller)
  CALL ListUntreatedFatal( Params,'Current Control',Caller)
  CALL ListUntreatedFatal( Params,'Calculate Volume Current',Caller)
  CALL ListUntreatedFatal( Params,'Calculate Volume Current',Caller)

END SUBROUTINE StatCurrentSolver_Init


!-----------------------------------------------------------------------------
!> A modern version for static current conduction supporting multithreading and
!> SIMD friendly ElmerSolver kernels. 
!------------------------------------------------------------------------------
SUBROUTINE StatCurrentSolver( Model,Solver,dt,TransientSimulation )
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
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr
  LOGICAL :: Found, VecAsm, InitHandles
  CHARACTER(*), PARAMETER :: Caller = 'StatCurrentSolver'
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving static current conduction solver')

  CALL DefaultStart()
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  nthr = 1
  !$ nthr = omp_get_max_threads()

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()

    totelem = 0

    nColours = GetNOFColours(Solver)
    VecAsm = (nColours > 1) .OR. (nthr == 1)
    
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
        CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm, InitHandles )
      END DO
      !$OMP END DO
    END DO
    !$OMP END PARALLEL 

    CALL CheckTimer(Caller//'BulkAssembly',Delete=.TRUE.)
    totelem = 0

    CALL DefaultFinishBulkAssembly()

    nColours = GetNOFBoundaryColours(Solver)
    VecAsm = (nColours > 1) .OR. (nthr == 1)

    CALL ResetTimer(Caller//'BCAssembly')
    
    !$OMP PARALLEL &
    !$OMP SHARED(Active, Solver, nColours, VecAsm) &
    !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
    !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
    DO col=1,nColours
       !$OMP SINGLE
       CALL Info('ModelPDEthreaded','Assembly of boundary colour: '//TRIM(I2S(col)),Level=10)
       Active = GetNOFBoundaryActive(Solver)
       !$OMP END SINGLE

       InitHandles = .TRUE. 
       !$OMP DO
       DO t=1,Active
          Element => GetBoundaryElement(t)
          ! WRITE (*,*) Element % ElementIndex
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
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJ(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, POINTER  :: CondCoeff(:), ReactCoeff(:), &
        EpsCoeff(:), SourceCoeff(:)
    REAL(KIND=dp) :: eps0
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes

    TYPE(ValueHandle_t), SAVE :: SourceCoeff_h, CondCoeff_h, ReactCoeff_h, EpsCoeff_h
    
    !$OMP THREADPRIVATE(Basis, dBasisdx, DetJ, &
    !$OMP               MASS, STIFF, FORCE, Nodes, &
    !$OMP               SourceCoeff_h, CondCoeff_h, ReactCoeff_h, EpsCoeff_h, &
    !$OMP               SourceCoeff, CondCoeff, ReactCoeff, EpsCoeff )
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJ
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------

    ! This InitHandles flag might be false on threaded 1st call
    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( SourceCoeff_h,'Body Force','Current Source')
      CALL ListInitElementKeyword( CondCoeff_h,'Material','Electric Conductivity')
      !CALL ListInitElementKeyword( ReactCoeff_h,'Material','Reaction Coefficient')
      CALL ListInitElementKeyword( EpsCoeff_h,'Material','Relative Permittivity')
      Eps0 = ListGetCReal( Model % Constants,'Permittivity Of Vacuum',Found )
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()
    IP = GaussPoints( Element )
    ngp = IP % n

    ! Deallocate storage if needed
    IF (ALLOCATED(Basis)) THEN
      IF (SIZE(Basis,1) < ngp .OR. SIZE(Basis,2) < nd) &
            DEALLOCATE(Basis, dBasisdx, DetJ, MASS, STIFF, FORCE )
    END IF

    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      ALLOCATE(Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJ(ngp), &
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
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJ, &
         SIZE(Basis,2), Basis, dBasisdx )

    ! Compute actual integration weights (recycle the memory space of DetJ)
    DO t=1,ngp
      DetJ(t) = IP % s(t) * Detj(t)
    END DO

    ! electric conductivity term: STIFF=STIFF+(rho*grad(u),grad(v))
    CondCoeff => ListGetElementRealVec( CondCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJ, STIFF, CondCoeff )
    END IF
   
    ! time derivative of potential: MASS=MASS+(eps*grad(u),grad(v))
    EpsCoeff => ListGetElementRealVec( EpsCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % DIMENSION, dBasisdx, DetJ, MASS, EpsCoeff )
      MASS(1:nd,1:nd) = Eps0 * MASS(1:nd,1:nd)
    END IF
      
    ! source term: FORCE=FORCE+(u,f)
    SourceCoeff => ListGetElementRealVec( SourceCoeff_h, ngp, Basis, Element, Found ) 
    IF( Found ) THEN
      CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, SourceCoeff, FORCE)
    END IF

    ! These are waiting for future use
    ! reaction term: STIFF=STIFF+(R*u,v)
    ! CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, STIFF, ReactCoeff )
    ! time derivative term    
    ! CALL LinearForms_UdotU(ngp, nd, Element % TYPE % DIMENSION, Basis, DetJ, MASS, EpsCoeff )
     
    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
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
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BC       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Robin_h, Ext_h

    SAVE Nodes
    !$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Current Density')
      CALL ListInitElementKeyword( Robin_h,'Boundary Condition','Electric Resistivity')
      CALL ListInitElementKeyword( Ext_h,'Boundary Condition','External Potential')
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
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

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
      C = ListGetElementReal( Robin_h, Basis, Element, Found )

      IF( Found ) THEN
        Ext = ListGetElementReal( Ext_h, Basis, Element, Found )
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
