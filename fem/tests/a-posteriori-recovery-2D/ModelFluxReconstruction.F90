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
!/*****************************************************************************/
! *
! * A prototype solver for advection-diffusion-reaction equation,
! * This equation is generic and intended for education purposes
! * but may also serve as a starting point for more complex solvers.
! *
! * Here additional post-processing code is also contained so as to
! * test some flux recovery/equilibration techniques. The exact solution
! * is hard-coded as this is used as a feasibility test. Otherwise
! * the code should be generic but still only in 2D.
! *
! * The author of post-processing part: mika.malinen@csc.fi 
! *
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/

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
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter
  LOGICAL :: Found
!------------------------------------------------------------------------------

  CALL DefaultStart()
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter

    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()
      CALL LocalMatrix(  Element, n, nd+nb )
    END DO

    CALL DefaultFinishBulkAssembly()

    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF(ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes()
        nd = GetElementNOFDOFs()
        CALL LocalMatrixBC(  Element, n, nd )
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()
    IF( DefaultConverged() ) EXIT    

  END DO

  CALL DefaultFinish()

  IF (ListGetLogical(Solver % Values, 'Flux Recovery', Found)) THEN
    !
    ! We might think of replacing the following call by calling RefineMesh with a new
    ! optional argument if CALL ComputeError therein were replaced by
    ! an alternate subroutine call to compute the error indicator
    !
    CALL FluxRecovery(Model, Solver)
  END IF
    
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: diff_coeff(n), conv_coeff(n),react_coeff(n), &
                     time_coeff(n), D,C,R, rho,Velo(3,n),a(3), Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    BodyForce => GetBodyForce()
    IF ( ASSOCIATED(BodyForce) ) &
       Load(1:n) = GetReal( BodyForce,'field source', Found )

    Material => GetMaterial()
    diff_coeff(1:n)=GetReal(Material,'diffusion coefficient',Found)
    react_coeff(1:n)=GetReal(Material,'reaction coefficient',Found)
    conv_coeff(1:n)=GetReal(Material,'convection coefficient',Found)
    time_coeff(1:n)=GetReal(Material,'time derivative coefficient',Found)

    Velo = 0._dp
    DO i=1,dim
      Velo(i,1:n)=GetReal(Material,&
          'convection velocity '//I2S(i),Found)
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
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
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
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp

    Flux(1:n)  = GetReal( BC,'field flux', Found )
    Coeff(1:n) = GetReal( BC,'robin coefficient', Found )
    Ext_t(1:n) = GetReal( BC,'external field', Found )

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
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE FluxRecovery(Model, Solver)
!------------------------------------------------------------------------------
! Flux recovery & a posteriori error estimation
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: RTFlux, NodalLoads
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: RTSolverPars
  TYPE(Element_t), POINTER :: Parent, Face
  INTEGER, ALLOCATABLE :: VisitsCounter(:)
  INTEGER, ALLOCATABLE, SAVE :: Indices(:), RT_Indices(:)
  INTEGER, POINTER :: FaceMap(:,:)
  LOGICAL, SAVE :: AllocationsDone = .FALSE.
  LOGICAL :: Found, UseReactions, Parallel, PostSmoothing, BDM
  LOGICAL :: ReverseSign(3)
  LOGICAL :: OrientationsMatch
  INTEGER :: n, nb, nd, t, active  
  INTEGER :: i, j, k, m, p, ni, nj, nd_rt, istat, ActiveFaceId
  INTEGER :: i_r, j_r
  REAL(KIND=dp), POINTER  :: ErrorIndicator(:)
  REAL(KIND=dp) :: UK(3), s, R1, R2, w1, w2, detw, r_i, r_j, hK
  REAL(KIND=dp) :: LinFun(8)   ! The size corresponds to RT_1(K)
  REAL(KIND=dp) :: Err, SolNorm, SolNormEst, Est, APostEst, Est_K
  REAL(KIND=dp), ALLOCATABLE :: RTFluxPost(:), ReactionWeights(:)
!------------------------------------------------------------------------------  
  
  ! We need a variable that is constructed as an approximation in RT_1/BDM
  !
  Mesh => GetMesh()
  RTFlux => VariableGet(Solver % Mesh % Variables, 'RTFlux')

  IF (ASSOCIATED(RTFlux)) THEN
    IF (.NOT. ASSOCIATED(RTFlux % Solver)) &
        CALL Fatal('FluxEquilibration', 'RTFlux variable is not associated with any solver')

    RTSolverPars => GetSolverParams(RTFlux % Solver)
    BDM = GetLogical(RTSolverPars, 'Second Kind Basis', Found)

    n = SIZE(RTFlux % Values)
    ALLOCATE( VisitsCounter(n), STAT=istat )
    VisitsCounter = 0

    CALL AllocateVector(ErrorIndicator, Solver % Mesh % NumberOfBulkElements)
    
    IF (.NOT. AllocationsDone) THEN
      N = Mesh % MaxElementDOFs
      ALLOCATE(Indices(N), RT_Indices(N), STAT=istat)
      AllocationsDone = .TRUE.
    ELSE
      IF (SIZE(Indices) < Mesh % MaxElementDOFs) THEN
        CALL Fatal('FluxEquilibration', 'mesh changed, the maximun counts of element DOFs are different?')
      END IF            
    END IF
  ELSE
    CALL Fatal('FluxEquilibration', 'RTFlux variable cannot be found')
  END IF

  
  !------------------------------------------------------------------------------
  ! Step I: 
  ! Compute the values of linear functionals (that is, DOFs) to obtain the elementwise
  ! representation of the flux (stress resultant) in RT_1(K) (or BDM(K)). We add the computed values
  ! to entries of the variable which is suitable for defining a conforming approximation
  ! in RT_1. A later averaging may be applied to ensure that a conforming approximation
  ! in RT_1 is obtained.
  !------------------------------------------------------------------------------
  Active = GetNOFActive()
  DO K=1,Active
    Element => GetActiveElement(K)
    IF (.NOT.(GetElementFamily(Element) == 3)) CYCLE

    n  = GetElementNOFNodes(Element)
    nd = GetElementDOFs(Indices, Element, Solver)
    nb = GetElementNOFBDOFs(Element, Solver)
    IF (nb > 0) CALL Fatal('FluxEquilibration', 'Bubbles in Global System = True assumed')
    
    UK(1:nd) = Solver % Variable % Values(Solver % Variable % Perm(Indices(1:nd)))

    nd_rt = GetElementDOFs(RT_Indices, Element, USolver = RTFlux % Solver)

    ! Compute the values of DOFs for the representation in RT_1(K)
    !
    CALL EstimateError(UK, Element, n, nd, LinFun = LinFun)
    
    DO i=1,nd_rt
      j = RTFlux % Solver % Variable % Perm(RT_Indices(i))
      RTFlux % Values(j) = RTFlux % Values(j) + LinFun(i)
      VisitsCounter(j) = VisitsCounter(j) + 1
    END DO
  END DO

  ! Averaging and applying the know constraints related to BCs. First, average:
  !
  DO i=1,SIZE(RTFlux % Values)
    IF (VisitsCounter(i) > 1) THEN
      RTFlux % Values(i) = RTFlux % Values(i) / VisitsCounter(i)
    END IF
  END DO

  NodalLoads => NULL()
  NodalLoads => VariableGet(Solver % Mesh % Variables, &
      GetVarName(Solver % Variable) // ' Loads' )
  
  IF (.NOT. ASSOCIATED(NodalLoads)) THEN
    UseReactions = .FALSE.
  ELSE  
    ! First, check whether there are reactions caused by BCs
    !
    IF (ALLOCATED(Solver % Matrix % ConstrainedDOF)) THEN
      IF (COUNT(Solver % Matrix % ConstrainedDOF) > 0) THEN
        UseReactions = .TRUE.
      ELSE
        UseReactions = .FALSE.
      END IF
    END IF
  END IF

  Apply_reactions: IF (UseReactions) THEN
    !
    ! Create data in order to average reactions 
    !
    ! First we tag DOFs on the model boundary by computing suitable weights
    ! for averaging
    ! TO DO: This should also be done for the DOFs which are associated with
    !        point loads. They need not necessarily be on the boundary.
    !        Smaller arrays could also be allocated by revising the code.
    !
    ALLOCATE(ReactionWeights(SIZE(Solver % Variable % Values)), STAT=istat)
    ReactionWeights = 0
    
    Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)

!    m = 0
    count_shared_nodes: DO K=1, GetNOFBoundaryElements()
      Element => GetBoundaryElement(K)
      IF (.NOT.(GetElementFamily(Element) == 2)) CYCLE
      IF (ActiveBoundaryElement()) THEN
        n = GetElementNOFNodes()
        nd = GetElementDOFs(Indices)
        
        IF (COUNT(Solver % Matrix % ConstrainedDOF(Solver % Variable % Perm(Indices(1:n)))) == n) THEN
!          m = m + 1
          hK = SQRT((Mesh % Nodes % x(Indices(2)) - Mesh % Nodes % x(Indices(1)))**2 + &
              (Mesh % Nodes % y(Indices(2)) - Mesh % Nodes % y(Indices(1)))**2)
          DO i=1,n
            ReactionWeights(Indices(i)) = ReactionWeights(Indices(i)) + hK
          END DO
        END IF
      END IF
    END DO count_shared_nodes

    ! Now, loop again to create the values of DOFs using the reactions. We already know what elements
    ! should be inspected, so this could be made more efficient ...
    !
    w1 = 0.5d0 * (1.0d0 + 1.0d0/sqrt(3.0d0))
    w2 = 0.5d0 * (1.0d0 - 1.0d0/sqrt(3.0d0))
    detw = w1**2 - w2**2
    
    replace_boundary_dofs: DO K=1, GetNOFBoundaryElements()
      Element => GetBoundaryElement(K)
      IF (.NOT.(GetElementFamily(Element) == 2)) CYCLE
      IF (ActiveBoundaryElement()) THEN
        n = GetElementNOFNodes()
        nd = GetElementDOFs(Indices)
        
        IF (COUNT(Solver % Matrix % ConstrainedDOF(Solver % Variable % Perm(Indices(1:n)))) == n) THEN
          !
          ! We need the parent to check the sign reversion
          !
          Parent => Element % BoundaryInfo % Left
          IF (.NOT. ASSOCIATED(Parent)) THEN
            Parent => Element % BoundaryInfo % Right
          END IF
          IF (.NOT. ASSOCIATED(Parent)) Call Fatal('FluxEquilibration', 'A parent element is not defined')  
          !
          ! Identify the face representing the element among the faces of 
          ! the parent element:
          !
          CALL PickActiveFace(Mesh, Parent, Element, Face, ActiveFaceId)
          IF (ActiveFaceId == 0) Call Fatal('FluxEquilibration', 'Cannot determine an element face')

          CALL FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)

          IF (IsLeftHanded(Parent)) THEN
            IF (ReverseSign(ActiveFaceId)) THEN
              s = 1.0d0
            ELSE
              s = -1.0d0
            END IF
          ELSE
            IF (ReverseSign(ActiveFaceId)) THEN
              s = -1.0d0
            ELSE
              s = 1.0d0
            END IF
          END IF
          
          FaceMap => GetEdgeMap(GetElementFamily(Parent))
          i_r = Element % NodeIndexes(1)
          j_r = Element % NodeIndexes(2)
          IF (ReverseSign(ActiveFaceId)) THEN
            i = Parent % NodeIndexes(FaceMap(ActiveFaceId,2))
            ! j = Parent % NodeIndexes(FaceMap(ActiveFaceId,1))
          ELSE
            i = Parent % NodeIndexes(FaceMap(ActiveFaceId,1))
            ! j = Parent % NodeIndexes(FaceMap(ActiveFaceId,2))
          END IF
          OrientationsMatch = i_r == i

          hK = SQRT((Mesh % Nodes % x(Indices(2)) - Mesh % Nodes % x(Indices(1)))**2 + &
              (Mesh % Nodes % y(Indices(2)) - Mesh % Nodes % y(Indices(1)))**2)
          
          IF (OrientationsMatch) THEN
            R1 = s * NodalLoads % Values(Solver % Variable % Perm(i_r)) * hk / ReactionWeights(i_r)
            R2 = s * NodalLoads % Values(Solver % Variable % Perm(j_r)) * hk / ReactionWeights(j_r)
          ELSE
            R1 = s * NodalLoads % Values(Solver % Variable % Perm(j_r)) * hK / ReactionWeights(j_r)
            R2 = s * NodalLoads % Values(Solver % Variable % Perm(i_r)) * hK / ReactionWeights(i_r)
          END IF

          nd_rt = GetElementDOFs(RT_Indices, Parent, USolver = RTFlux % Solver)
!          r_i = NodalLoads % Values(Solver % Variable % Perm(i_r))
!          r_j = NodalLoads % Values(Solver % Variable % Perm(j_r))
!          print *, '===='
!          print *, '=== treating boundary element', K
!          print *, '=== boundary element indexes', i_r, j_r
!          print *, '==== corresponding reactions', r_i, r_j
!          print *, '==== corresponding reactions averaged', r_i * hK / ReactionWeights(i_r), &
!              r_j * hK / ReactionWeights(j_r)
!          print *, '=== parent element indexes', i, j
!          print *, '=== orientations match', OrientationsMatch
!          print *, '=== parent left-handed', IsLeftHanded(Parent), Parent % NodeIndexes(1:3)
!          print *, '=== the first DOF now/before', R1, &
!              RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId - 1)))
!          print *, '=== the second DOF now/before', R2, &
!              RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId)))

          ! Finally use the reactions to replace the values of DOFs on the boundary
          !
          IF (BDM) THEN
            RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId - 1))) = &
                1.0d0/detw * (w1*R1 - w2*R2)
            RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId))) = &
                1.0d0/detw * (-w2*R1 + w1*R2)
          ELSE
            RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId - 1))) = R1
            RTFlux % Values(RTFlux % Solver % Variable % Perm(RT_Indices(2*ActiveFaceId))) = R2
          END IF
        END IF
      END IF
    END DO replace_boundary_dofs
  END IF Apply_reactions

  !------------------------------------------------------------------------------
  ! Step II: 
  ! Equilibrate fluxes (stress resultants). This brings us back to the recovery
  ! which is defined only in the local RT space.
  ! Here the exact solution is also used to study the accuracy.
  !------------------------------------------------------------------------------

  PostSmoothing = .FALSE. ! Eliminated as this doesn't seem to be beneficial
  IF (PostSmoothing) THEN
    ALLOCATE(RTFluxPost(SIZE(RTFlux % Values)))
    RTFluxPost(:) = 0.0d0
    VisitsCounter(:) = 0
  END IF
  
  Err = 0.0d0
  SolNorm = 0.0d0
  Est = 0.0d0
  APostEst = 0.0d0
  SolNormEst = 0.0d0
  
  Elementwise_equilibration: DO K=1,Active
    Element => GetActiveElement(K)
    IF (.NOT. (GetElementFamily(Element) == 3)) CYCLE

    n  = GetElementNOFNodes()
    nd = GetElementDOFs(Indices)
    nb = GetElementNOFBDOFs()
    IF (nb > 0) CALL Fatal('FluxEquilibration', 'Bubbles in Global System = True assumed')
    
    UK(1:nd) = Solver % Variable % Values( Solver % Variable % Perm(Indices(1:nd)) )

    nd_rt = GetElementDOFs(RT_Indices, Element, USolver = RTFlux % Solver)

    ! At the calling time, we still have the RT DOFs respecting the global continuity requirements
    !
    DO i=1,nd_rt
      j = RTFlux % Solver % Variable % Perm(RT_Indices(i))
      LinFun(i) = RTFlux % Values(j)
    END DO
    
    CALL EstimateError(UK, Element, n, nd, Err, SolNorm, Est, APostEst, PostLinFun = LinFun, &
        APostEst_K = Est_K, SolNormEst=SolNormEst)

    ErrorIndicator(K) = SQRT(Est_K)

    IF (PostSmoothing) THEN
      DO i=1,nd_rt
        j = RTFlux % Solver % Variable % Perm(RT_Indices(i))
        RTFluxPost(j) = RTFluxPost(j) + LinFun(i)
        VisitsCounter(j) = VisitsCounter(j) + 1
      END DO
    END IF
    
  END DO Elementwise_equilibration

  IF (SolNormEst > AEPS) ErrorIndicator = (1.0d0 / SQRT(SolNormEst)) * ErrorIndicator
  
  WRITE (*, '(A,E16.8)') 'Solution Norm = ', SQRT(ParallelReduction(SolNorm))
  WRITE (*, '(A,E16.8)') 'Error Norm = ', SQRT(ParallelReduction(Err))/SQRT(ParallelReduction(SolNorm))
  WRITE (*, '(A,E16.8)') 'Recovery Error Norm = ', SQRT(ParallelReduction(Est))/SQRT(ParallelReduction(SolNorm))
  WRITE (*, '(A,E16.8)') 'A posteriori Error = ', SQRT(ParallelReduction(APostEst))/SQRT(ParallelReduction(SolNorm))
  WRITE (*, '(A,E16.8)') 'A posteriori Error Est = ', SQRT(SUM(ErrorIndicator**2))
  WRITE (*, '(A,E16.8)') 'Efficiency Factor = ', SQRT(ParallelReduction(APostEst))/SQRT(ParallelReduction(Err))
!  CALL ShowVectorHistogram(ErrorIndicator,SIZE(ErrorIndicator))
  
  !
  ! The following computation would average again to obtain the recovery field in the global RT_1 space
  !
  post_smoothing: IF (PostSmoothing) THEN

    Err = 0.0d0
    Est = 0.0d0
    APostEst = 0.0d0
    SolNorm = 0.0d0

    ! Post averaging
    !
    DO i=1,SIZE(RTFlux % Values)
      IF (VisitsCounter(i) > 1) THEN
        RTFluxPost(i) = RTFluxPost(i) / VisitsCounter(i)
      END IF
    END DO

    Err = 0.0d0
    Est = 0.0d0
    APostEst = 0.0d0
    SolNorm = 0.0d0

    DO K=1,Active
      Element => GetActiveElement(K)
      IF ( .NOT. (GetElementFamily(Element) == 3) ) CYCLE

      n  = GetElementNOFNodes()
      nd = GetElementDOFs(Indices)
      nb = GetElementNOFBDOFs()
      IF (nb > 0) CALL Fatal('FluxEquilibration', 'Bubbles in Global System = True assumed')

      UK(1:nd) = Solver % Variable % Values( Solver % Variable % Perm(Indices(1:nd)) )

      nd_rt = GetElementDOFs(RT_Indices, Element, USolver = RTFlux % Solver)

      DO i=1,nd_rt
        j = RTFlux % Solver % Variable % Perm(RT_Indices(i))
        LinFun(i) = RTFluxPost(j)
      END DO

      CALL EstimateError(UK, Element, n, nd, Err, SolNorm, Est, APostEst, PostLinFun = LinFun, &
          UseGiven = .TRUE.)

    END DO

    WRITE (*, '(A,E16.8)') 'Solution Norm = ', SQRT(ParallelReduction(SolNorm))
    WRITE (*, '(A,E16.8)') 'Error Norm = ', SQRT(ParallelReduction(Err))/SQRT(ParallelReduction(SolNorm))
    WRITE (*, '(A,E16.8)') 'Estimated Error Norm = ', SQRT(ParallelReduction(Est))/SQRT(ParallelReduction(SolNorm))
    WRITE (*, '(A,E16.8)') 'A posteriori Error Est = ', SQRT(ParallelReduction(APostEst))/SQRT(ParallelReduction(SolNorm))
    WRITE (*, '(A,E16.8)') 'Efficiency Factor = ', SQRT(ParallelReduction(APostEst))/SQRT(ParallelReduction(Err))

  END IF post_smoothing

  DEALLOCATE(ErrorIndicator)
  
CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE EstimateError(UK, Element, n, nd, Err, SolNorm, Est, APostEst, &
      LinFun, PostLinFun, UseGiven, APostEst_K, SolNormEst)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(IN) :: UK(:)
    TYPE(Element_t), POINTER, INTENT(IN) :: Element    
    INTEGER, INTENT(IN) :: n, nd
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: Err, SolNorm, Est, APostEst
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: LinFun(8)
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: PostLinFun(8)
    LOGICAL, OPTIONAL, INTENT(IN) :: UseGiven
    REAL(KIND=dp), OPTIONAL, INTENT(INOUT) :: APostEst_K, SolNormEst
!--------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    
    REAL(KIND=dp) :: FBasis(8,3), DivFBasis(8), Mass(11,11), RHS(11), c(11), d(11), s
    REAL(KIND=dp) :: Basis(nd), dBasis(nd,3), DetJ, xt, uq, vq, weight, EA, f
    REAL(KIND=dp) :: U, gradU(3), Uh, gradUh(3), Nh(3), intf, divN, totflux, dc
    REAL(KIND=dp) :: testfun(2), diff_coeff(n), Load(n)
    REAL(KIND=dp) :: faceflux(3), faceweights(3), facedelta(3)
    REAL(KIND=dp) :: w1, w2
    LOGICAL :: Stat, Found
    LOGICAL :: ReverseSign(3), LeftHanded, Parallel, UseLM
    LOGICAL :: FirstOrderEquilibration 
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: t, i, j, p, q, ni, nj, np, FDOFs, DOFs, ElementOrder
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
!------------------------------------------------------------------------------
    IF (BDM) THEN
      FDOFs = 6
      ElementOrder = 1
      FirstOrderEquilibration = .FALSE.
    ELSE
      FDOFs = 8
      ElementOrder = 2
      FirstOrderEquilibration = .TRUE.
    END IF

    UseLM = .TRUE.
    IF (UseLM) THEN
      IF (BDM) THEN
        DOFs = FDOFs + 1
      ELSE
        DOFs = FDOFs + 3
      END IF
    ELSE
      DOFs = FDOFs
    END IF
       
    CALL GetElementNodes( Nodes )
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE., EdgeBasisDegree=2)

    Material => GetMaterial()
    diff_coeff(1:n) = GetReal(Material, 'diffusion coefficient', Found)

    IF (PRESENT(APostEst_K)) APostEst_K = 0.0d0
    
    c = 0.0d0
    IF (PRESENT(UseGiven) .AND. PRESENT(PostLinFun)) THEN
      IF (UseGiven) THEN
        c(1:FDOFs) = PostLinFun(1:FDOFs)
        GOTO 303
      END IF
    END IF

    Load = 0.0d0
    BodyForce => GetBodyForce()
    IF (ASSOCIATED(BodyForce)) &
       Load(1:n) = GetReal( BodyForce,'field source', Found )

    
    IF (PRESENT(LinFun)) THEN
      !
      ! Compute a local representation of the flux in the elementwise RT_1(K) space.
      ! The elementwise weak formulation is of the form
      !           1/k (N,v) = -(u,div v) + <u,v.n> 
      !
      Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)
      EdgeMap => GetEdgeMap(3)
      CALL FaceElementOrientation(Element, ReverseSign)

      Mass = 0.0_dp
      RHS = 0.0_dp
      IF (BDM) THEN 
        w1 = 0.5d0 * (1.0d0 + 1.0d0/sqrt(3.0d0))
        w2 = 0.5d0 * (1.0d0 - 1.0d0/sqrt(3.0d0))
      !ELSE
      !  w1 = 1.0d0
      !  w2 = 0.0d0
      END IF
      
      DO t=1,IP % n

        stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detF=detJ, Basis=Basis, FBasis=FBasis, &
            DivFBasis=DivFBasis, BDM=BDM, BasisDegree=ElementOrder, ApplyPiolaTransform=.TRUE., &
            LeftHanded=LeftHanded)

        Weight = IP % s(t) * DetJ

        Uh = SUM(UK(1:nd) * Basis(1:nd))
        EA = SUM(diff_coeff(1:n) * Basis(1:n))
        
        DO p=1,FDOFs
          DO q=1,FDOFs
            Mass(p,q) = Mass(p,q) + SUM(FBasis(q,1:2) * FBasis(p,1:2)) * Weight / EA
          END DO
          RHS(p) = RHS(p) - Weight * Uh * DivFBasis(p)
        END DO
        
        IF (UseLM) THEN
          !
          ! Enforce the zeroth-order equilibration by using a Lagrange multiplier
          !
          f = SUM(Load(1:n) * Basis(1:n))
          RHS(FDOFs+1) = RHS(FDOFs+1) + Weight * f

          IF (FirstOrderEquilibration) THEN
            !
            ! Enforce the first-order equilibration by using a Lagrange multiplier
            !
            testfun(1) = Basis(2) - Basis(1)
            testfun(2) = Basis(3) - Basis(1)
            DO p=1,FDOFs
              Mass(10,p)= Mass(10,p) - divFBasis(p) * testfun(1) * Weight
              Mass(11,p)= Mass(11,p) - divFBasis(p) * testfun(2) * Weight
              RHS(10) = RHS(10) + f * testfun(1) * Weight
              RHS(11) = RHS(11) + f * testfun(2) * Weight
              Mass(p,10)= Mass(p,10) - divFBasis(p) * testfun(1) * Weight
              Mass(p,11)= Mass(p,11) - divFBasis(p) * testfun(2) * Weight
            END DO
          END IF
        END IF
      END DO

      ! Assemble the boundary terms:
      ! Loop over all faces (here edges)
      !
      DO p=1,3
        ! Check whether the sign is reversed
        IF (LeftHanded) THEN
          IF (ReverseSign(p)) THEN
            s = 1.0d0
          ELSE
            s = -1.0d0
          END IF
        ELSE
          IF (ReverseSign(p)) THEN
            s = -1.0d0
          ELSE
            s = 1.0d0
          END IF
        END IF
        ! Check the order of basis functions
        i = EdgeMap(p,1)
        j = EdgeMap(p,2)
        ni = Element % NodeIndexes(i)
        IF (Parallel) ni = Mesh % ParallelInfo % GlobalDOFs(ni)             
        nj = Element % NodeIndexes(j)
        IF (Parallel) nj = Mesh % ParallelInfo % GlobalDOFs(nj)

        IF (BDM) THEN
          IF (nj<ni) THEN
            RHS(2*p-1) = RHS(2*p-1) + s * (w2 * UK(i) + w1 * UK(j))
            RHS(2*p) = RHS(2*p) + s * (w1 * UK(i) + w2 * UK(j))
          ELSE
            RHS(2*p-1) = RHS(2*p-1) + s * (w1 * UK(i) + w2 * UK(j))
            RHS(2*p) = RHS(2*p) + s * (w2 * UK(i) + w1 * UK(j))
          END IF
        ELSE
          ! The value of RHS vector is just the nodal DOF up to the sign
          IF (nj<ni) THEN
            ! The first weight function is s * basis(j)
            RHS(2*p-1) = RHS(2*p-1) + s * UK(j)
            RHS(2*p) = RHS(2*p) + s * UK(i)
          ELSE
            RHS(2*p-1) = RHS(2*p-1) + s * UK(i)
            RHS(2*p) = RHS(2*p) + s * UK(j)
          END IF
        END IF
          
        IF (UseLM) THEN
          Mass(FDOFs+1,2*p-1) = -s
          Mass(FDOFs+1,2*p) = -s
          Mass(2*p-1,FDOFs+1) = -s
          Mass(2*p,FDOFs+1) = -s
        END IF
        
      END DO

      ! Finally solve the local representation of the flux
      CALL InvertMatrix(Mass(1:DOFs,1:DOFs), DOFs)
      c(1:DOFs) = MATMUL(MASS(1:DOFs,1:DOFs), RHS(1:DOFs))
      LinFun(1:FDOFs) = c(1:FDOFs)
      RETURN
    END IF

    IF (PRESENT(PostLinFun)) THEN
      c(1:FDOFs) = PostLinFun(1:FDOFs)

      EdgeMap => GetEdgeMap(3)
      CALL FaceElementOrientation(Element, ReverseSign)
      Parallel = ASSOCIATED(Mesh % ParallelInfo % GInterface)
      
      IF (.NOT. BDM) THEN
        !
        ! Perform the zeroth-order equilibration (note that BDM does not seem to benefit
        ! from this):
        !
        intf = 0.0d0
        divN = 0.0d0
        
        DO t=1,IP % n

          stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detF=detJ, Basis=Basis, FBasis=FBasis, &
              DivFBasis=DivFBasis, BDM=BDM, BasisDegree=ElementOrder, ApplyPiolaTransform=.TRUE., &
              LeftHanded=LeftHanded)

          Weight = IP % s(t) * DetJ
          f = SUM(Load(1:n) * Basis(1:n))
          intf = intf + f * Weight
        END DO

        ! Loop over all faces (here edges)
        !
        totflux = 0.0d0
        faceflux = 0.0d0
        DO p=1,3
          IF (LeftHanded) THEN
            IF (ReverseSign(p)) THEN
              s = 1.0d0
            ELSE
              s = -1.0d0
            END IF
          ELSE
            IF (ReverseSign(p)) THEN
              s = -1.0d0
            ELSE
              s = 1.0d0
            END IF
          END IF
          totflux = totflux + s*(c(2*p-1)+c(2*p))
          faceflux(p) = faceflux(p) + s*(c(2*p-1)+c(2*p))
!          IF (c(2*p-1)*c(2*p) < 0.0 .AND. -c(2*p-1)*c(2*p)/(MAXVAL(ABS(c(:))))**2 > 1.0d-16) THEN
!            CALL Warn('FluxEquilibration', 'No consensus on flux sign')
!            print *, 'c1', c(2*p-1)
!            print *, 'c2', c(2*p)
!          END IF
        END DO

        DO p=1,3
          FaceWeights(p) = abs(faceflux(p))/sum(abs(faceflux(:)))
        END DO
        
!        print *, 'total flux out', totflux
!        print *, 'sources tot', -intf
!        print *, 'change flux out by an amount', -intf - totflux
!        dc = (-intf - totflux)/6.0d0
        
        DO p=1,3
          facedelta(p) = FaceWeights(p) * (-intf - totflux)   
        END DO        
        !
        ! The redefinition to ensure the zeroth-order equilibration:
        !
        DO p=1,3
          IF (LeftHanded) THEN
            IF (ReverseSign(p)) THEN
              s = 1.0d0
            ELSE
              s = -1.0d0
            END IF
          ELSE
            IF (ReverseSign(p)) THEN
              s = -1.0d0
            ELSE
              s = 1.0d0
            END IF
          END IF
          
!          d(2*p-1) = c(2*p-1) + s*dc
!          d(2*p) = c(2*p) + s*dc
          d(2*p-1) = c(2*p-1) + s*facedelta(p)*0.5d0
          d(2*p) = c(2*p) + s*facedelta(p)*0.5d0
        END DO

        post_check: IF (.TRUE.) THEN
          totflux = 0.0d0
          DO p=1,3
            IF (LeftHanded) THEN
              IF (ReverseSign(p)) THEN
                s = 1.0d0
              ELSE
                s = -1.0d0
              END IF
            ELSE
              IF (ReverseSign(p)) THEN
                s = -1.0d0
              ELSE
                s = 1.0d0
              END IF
            END IF

            totflux = totflux + s*(d(2*p-1)+d(2*p))
          END DO

          IF (ABS(intf + totflux) > 1.0d1 * AEPS) THEN 
            print *, 'Warning: change flux out by an amount', -intf - totflux
          END IF
        END IF post_check

        c(1:6) = d(1:6)
      END IF

      ! Apply the first-order equilibration
      !
      IF (FirstOrderEquilibration) THEN
        Mass = 0.0_dp
        RHS = 0.0_dp
        
        DO t=1,IP % n
          
          stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detF=detJ, Basis=Basis, FBasis=FBasis, &
              DivFBasis=DivFBasis, BasisDegree=2, ApplyPiolaTransform=.TRUE.)

          Weight = IP % s(t) * DetJ
          f = SUM(Load(1:n) * Basis(1:n))
          
          testfun(1) = Basis(2) - Basis(1)
          testfun(2) = Basis(3) - Basis(1)

          DO p=1,2
            DO q=1,2
              Mass(p,q) = Mass(p,q) - divFBasis(6+q) * testfun(p) * Weight
            END DO
            RHS(p) = RHS(p) + f * testfun(p) * Weight
            DO q=1,6
              RHS(p) = RHS(p) + c(q) * divFBasis(q) * testfun(p) * Weight
            END DO
          END DO
        END DO
        
        CALL InvertMatrix(Mass(1:2,1:2),2)
        c(7:8) = MATMUL(MASS(1:2,1:2), RHS(1:2))
      END IF
      PostLinFun(1:FDOFs) = c(1:FDOFs)
    END IF
    
    ! Compute a posteriori estimate:
303 CONTINUE
    DO t=1,IP % n

      ! Now get the face basis functions so that we can evaluate
      ! the flux in terms of them
      !
      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detF=detJ, Basis=Basis, FBasis=FBasis, &
            DivFBasis=DivFBasis, dBasisdx=dBasis, BDM=BDM, BasisDegree=ElementOrder, &
            ApplyPiolaTransform=.TRUE.)
      
      Weight = IP % s(t) * DetJ

      !  The exact solution:
      IF (PRESENT(Err) .OR. PRESENT(SolNorm)) THEN
        xt = SUM(Basis(1:n) * Nodes % x(1:n))
        U = 1.0d0 - xt**2/9.0d0
        gradU(:) = 0.0
        gradU(1) = -2.0d0 * xt/9.0d0
      END IF

      Uh = SUM(UK(1:nd) * Basis(1:nd))
      DO i=1,2
        gradUh(i) = SUM(UK(1:nd) * dBasis(1:nd,i))
      END DO
      gradUh(3) = 0.0d0
      
      EA = SUM(diff_coeff(1:n) * Basis(1:n))
      Nh = 0.0d0
      DO i=1,2
        Nh(i) = SUM( c(1:FDOFs) * FBasis(1:FDOFs,i) )
      END DO
      
      IF (.TRUE.) THEN
        ! The error of flux
        IF (PRESENT(Err)) Err = Err + &
            SUM((EA * gradU(1:2) - EA * gradUh(1:2))**2) * Weight
        IF (PRESENT(SolNorm)) SolNorm = SolNorm + &
            SUM((EA * gradU(1:2))**2) * Weight
        IF (PRESENT(SolNormEst)) SolNormEst = SolNormEst + &
            SUM((Nh(1:2))**2) * Weight
        
        ! A posteriori estimate based on the recovery:
        !
        IF (PRESENT(Est)) Est = Est + SUM((EA * gradU(1:2) - Nh(1:2))**2) * Weight
        IF (PRESENT(APostEst)) APostEst = APostEst + SUM((EA * gradUh(1:2) - Nh(1:2))**2) * Weight
        IF (PRESENT(APostEst_K)) APostEst_K = APostEst_K + SUM((EA * gradUh(1:2) - Nh(1:2))**2) * Weight
      ELSE
        ! The error of the field
        IF (PRESENT(Err)) Err = Err + &
            (U - Uh)**2 * Weight
        IF (PRESENT(SolNorm)) SolNorm = SolNorm + &
            U**2 * Weight
      END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE EstimateError
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION IsLeftHanded(Element) RESULT(LeftHanded)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    LOGICAL :: LeftHanded
!------------------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: Stat
    REAL(KIND=dp) :: Basis(Element % Type % NumberOfNOdes), DetJ
    REAL(KIND=dp) :: FBasis(8,3)
!------------------------------------------------------------------------------
    CALL GetElementNodes(Nodes, Element)    
    
    stat = FaceElementInfo(Element, Nodes, 0.0d0, SQRT(3.0d0)/3.0d0, &
            0.0d0, detF=detJ, Basis=Basis, FBasis=FBasis, &
            BasisDegree=1, ApplyPiolaTransform=.TRUE., &
            LeftHanded=LeftHanded)
!------------------------------------------------------------------------------    
  END FUNCTION IsLeftHanded
!------------------------------------------------------------------------------  
    
!------------------------------------------------------------------------------
END SUBROUTINE FluxRecovery
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE FluxEquilibration_Dummy(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
END SUBROUTINE FluxEquilibration_Dummy
!------------------------------------------------------------------------------
