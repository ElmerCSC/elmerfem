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
! *
! *  Solve the mixed formulation of the generalized Poisson equation by using 
! *  div-conforming (face) finite elements.
! *
! *  NOTE: It is assumed that the last bubble DOF is used for approximating 
! *        the scalar variable. That is, the scalar variable is approximated 
! *        as an elementwise constant.
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Feb 13, 2019
! *
!******************************************************************************


!------------------------------------------------------------------------------
SUBROUTINE MixedPoisson_Init0(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverPars
  LOGICAL :: Found, SecondFamily
  INTEGER :: dim
  CHARACTER(LEN=MAX_NAME_LEN) :: csys, VarName

!------------------------------------------------------------------------------
  SolverPars => GetSolverParams()
  SecondFamily = GetLogical(SolverPars, 'Second Kind Basis', Found)

  csys = ListGetString(Model % Simulation, 'Coordinate System', Found)
  IF (.NOT. Found) THEN 
    IF (.NOT. ListCheckPresent(SolverPars, 'Element')) &
        CALL Fatal('MixedPoisson_Init0', 'The keyword Element should be specified')
  ELSE
    !
    ! The coordinate system dimension cannot yet be returned by the function
    ! CoordinateSystemDimension due to early execution of this initialization
    ! subroutine;  instead the value of "Coordinate System" is employed.
    !
    SELECT CASE (csys)
    CASE('cartesian 2d')
      IF (SecondFamily) THEN
        CALL ListAddNewString(SolverPars, "Element", "n:0 e:2 b:1")
      ELSE
        CALL ListAddNewString(SolverPars, "Element", "n:0 e:1 -tri b:1 -quad b:3")
      END IF

    CASE('cartesian 3d')
      IF (SecondFamily) THEN
        CALL ListAddNewString(SolverPars, "Element", &
            "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:3")
      ELSE
        CALL ListAddNewString(SolverPars, "Element", &
            "n:0 -tetra b:1 -brick b:25 -quad_face b:4 -tri_face b:1")
      END IF

    CASE DEFAULT
      IF (.NOT. ListCheckPresent(SolverPars, 'Element')) &
          CALL Fatal('MixedPoisson_Init0', 'The keyword Element should be specified')     
    END SELECT
  END IF

  CALL ListAddNewLogical(SolverPars, 'Bubbles in Global System', .TRUE.)
  
  ! Add scalar variable if not present, and get its name
  CALL ListAddNewString(SolverPars,'Potential Variable','mixedpot' )
  VarName = ListGetString(SolverPars,'Potential Variable')
    
  CALL ListAddString( SolverPars,NextFreeKeyword(&
      'Exported Variable',SolverPars),'-elem '//TRIM(VarName))

  CALL ListAddString( SolverPars,NextFreeKeyword(&
      'Exported Variable',SolverPars), TRIM(VarName))
  
  CALL ListAddNewString(SolverPars,'Flux Variable','mixedflux' )
  VarName = ListGetString(SolverPars,'Flux Variable')

  CALL ListAddString( SolverPars,NextFreeKeyword(&
      'Exported Variable',SolverPars),'-dofs 3 -dg '//TRIM(VarName))

  CALL ListAddString( SolverPars,NextFreeKeyword(&
      'Exported Variable',SolverPars), TRIM(VarName))

  CALL ListAddNewString( SolverPars,'Element Integration Points',&
      '-tri 3 -quad 9 -tetra 4 -brick 64')

  ! This solver always needs PostSolver since primary fields are not very intuitive
  CALL ListAddLogical( SolverPars,'PostSolver Active',.TRUE.)
  
!------------------------------------------------------------------------------
END SUBROUTINE MixedPoisson_Init0
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE MixedPoisson(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables:
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: AllocationsDone = .FALSE.
  LOGICAL :: Found, InitHandles, SecondFamily
  INTEGER :: dim, n, nb, nd, t, istat, active
  REAL(KIND=dp), ALLOCATABLE :: Stiff(:,:), Mass(:,:), Force(:)
  REAL(KIND=dp) :: Norm, Cond


  SAVE Stiff, Mass, Force, AllocationsDone
!------------------------------------------------------------------------------
  CALL DefaultStart()

  dim = CoordinateSystemDimension()

  ! Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs  ! just big enough
    ALLOCATE( Force(N), Stiff(N,N), Mass(N,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'MixedPoisson', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF

  !-----------------------
  ! System assembly:
  !----------------------
  SecondFamily = GetLogical(GetSolverParams(), 'Second Kind Basis', Found)

  active = GetNOFActive()
  CALL DefaultInitialize()

  InitHandles = .TRUE.
  DO t=1,active
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() ! Nodes count corresponding to the background mesh
    nd = GetElementNOFDOFs()  ! The total number of degrees of freedom
    nb = SIZE(Element % BubbleIndexes(:)) ! The number of elementwise degrees 
                                          ! of freedom. NOTE: GetElementNOFBDOFs()
                                          ! doesn't return the right value here 

    ! Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(Stiff, Mass, Force, Element, n, nd, nb, dim, SecondFamily, &
        TransientSimulation, InitHandles )

    IF (TransientSimulation) CALL Default1stOrderTime(Mass, Stiff, Force)  

    ! Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations(Stiff, Force)

  END DO

  CALL DefaultFinishBulkAssembly()

  InitHandles = .TRUE.
  active = GetNOFBoundaryElements()
  DO t=1,active
    Element => GetBoundaryElement(t)
    IF (ActiveBoundaryElement()) THEN
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)

      Cond = SUM(GetReal(GetBC(), GetVarName(Solver % Variable)//' Condition', Found))/n
      IF(Cond >= 0) CALL LocalMatrixBC(Element, Mesh, n, nd, SecondFamily, InitHandles)
    END IF
  END DO

  CALL DefaultFinishBoundaryAssembly()

  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()  

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Stiff, Mass, Force, Element, n, nd, nb, dim, &
      SecondFamily, NeedMass, InitHandles )
!------------------------------------------------------------------------------
    USE MGDynMaterialUtils, ONLY : GetTensor

    REAL(KIND=dp) :: Stiff(:,:), Mass(:,:), Force(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: nd  ! The total count of DOFs (nodal, facial and elementwise)
    INTEGER :: nb  ! The number of elementwise DOFs (for the scalar and flux variables)
    INTEGER :: dim
    LOGICAL :: SecondFamily
    LOGICAL :: NeedMass
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MaxFaceBasisDim = 48
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np, mu_rank
    REAL(KIND=dp), POINTER :: mu_tensor(:,:)
    REAL(KIND=dp) :: a, f, v(3), tmp(dim), mu
    REAL(KIND=dp) :: FaceBasis(MaxFaceBasisDim,3), DivFaceBasis(MaxFaceBasisDim)
    REAL(KIND=dp) :: Basis(nd), DetJ, s
    TYPE(ValueHandle_t), SAVE :: SourceField_h, ConvVelo_h, MatPar_h, MatTensor_h
    
!------------------------------------------------------------------------------
    IF(InitHandles ) THEN
      CALL ListInitElementKeyword( SourceField_h, 'Body Force','Source Field')

      ! The original implementation uses inverse for material parameter, but not for material tensor?
      ! Hence we need to separate the two!
      CALL ListInitElementKeyword( MatPar_h,'Material','Material Parameter',DefRValue=1.0_dp)
      CALL ListInitElementKeyword( MatTensor_h,'Material','Material Tensor')
      
      CALL ListInitElementKeyword( ConvVelo_h,'Material','Convection Velocity',InitVec3D=.TRUE.)
      InitHandles = .FALSE.
    END IF


    CALL GetElementNodes( Nodes )

    Stiff = 0.0d0
    Mass = 0.0d0
    Force = 0.0d0

    !------------------------------------------------------------------------
    ! The reference element is chosen to be that used for p-approximation,
    ! so we need to switch to using a quadrature which would not be used 
    ! otherwise
    !------------------------------------------------------------------------
    IF( .FALSE. ) THEN
      ! This rule should be ok, but there seems to be two sets of IPs for
      ! reference element. 
      IP = GaussPointsAdapt( Element, PReferenceElement = .TRUE. )
    ELSE
      SELECT CASE( GetElementFamily(Element) )
      CASE(3)
        IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
      CASE(4)
        IP = GaussPointsQuad(9)
      CASE(5)
        IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
      CASE(8)
        IP = GaussPointsBrick(64)
      CASE DEFAULT
        CALL Fatal('MixedPoisson', 'A non-supported element type')
      END SELECT
    END IF

#if 0
    IF( Element % ElementIndex == 2 ) THEN
      n = IP % n
      PRINT *,'IP:',Element % TYPE % ElementCode, n
      PRINT *,'IP u:',IP % u(1:n)
      PRINT *,'IP v:',IP % v(1:n)
      PRINT *,'IP w:',IP % w(1:n)
      PRINT *,'IP s:',IP % s(1:n)
    END IF
#endif
    
    ! Set np = n, if nodal dofs are employed; otherwise set np = 0:
    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)    

    DO t=1,IP % n
      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)


      s = detJ * IP % s(t)

      !----------------------------------------------------------------
      ! The following branch could be used to produce the 
      ! Galerkin projection of the pressure for visualization.
      !------------------------------------------------------------------
      IF (np > 0) THEN
        DO p = 1,n
          DO q = 1,n       
            Stiff(p,q) = Stiff(p,q) + Basis(p) * Basis(q) * s    
          END DO

          DO q = nd,nd
            Stiff(p,q) = Stiff(p,q) - Basis(p) * 1.0d0 * s            
          END DO
        END DO
      END IF

      !--------------------------------------------------------------
      ! The contribution from the variation with the flux variable q
      !---------------------------------------------------------------
      mu = ListGetElementReal( MatTensor_h, Basis, Element, Found, GaussPoint = t, &
          Rdim = mu_rank, Rtensor = mu_tensor)
      ! Note that if mu is scalar then mu=1/a.
      IF( Found ) THEN
        IF( mu_rank == 0 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              Stiff(i,j) = Stiff(i,j) + mu * &
                  SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * s           
            END DO
          END DO
        ELSE IF (mu_rank == 1 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              tmp = SUM(mu_tensor(1:dim,1)*FaceBasis(q,1:dim))
              Stiff(i,j) = Stiff(i,j) + & 
                  SUM( tmp(1:dim) * FaceBasis(p,1:dim) ) * s
            END DO
          END DO
        ELSE IF (mu_rank == 2 ) THEN
          DO p = 1,nd-np-1
            i = np + p
            DO q = 1,nd-np-1
              j = np + q
              tmp = MATMUL(mu_tensor(1:dim,1:dim),FaceBasis(q,1:dim))
              Stiff(i,j) = Stiff(i,j) + & 
                  SUM( tmp(1:dim) * FaceBasis(p,1:dim) ) * s
            END DO
          END DO
        END IF
      ELSE
        ! If this is not given then it is defaulted to one.
        a = ListGetElementReal( MatPar_h, Basis, Element, Found, GaussPoint = t )      
        DO p = 1,nd-np-1
          i = np + p
          DO q = 1,nd-np-1
            j = np + q
            Stiff(i,j) = Stiff(i,j) + (1.0_dp / a ) * &
                SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * s           
          END DO
        END DO
      END IF

        
      DO p = 1,nd-np-1
        i = np + p
        DO q = nd,nd
          Stiff(i,q) = Stiff(i,q) + 1.0d0 * DivFaceBasis(p) * s
        END DO
      END DO

      !--------------------------------------------------
      ! The contribution from the constraint div q = -f
      !--------------------------------------------------
      DO p = nd,nd
        DO q = 1,nd-np-1
          j = np + q
          Stiff(p,j) = Stiff(p,j) + 1.0d0 * DivFaceBasis(q) * s
        END DO
      END DO

      ! Contribution of convection
      !-----------------------------------------------------
      v = ListGetElementReal3D( ConvVelo_h, Basis, Element, Found, GaussPoint = t )      
      IF ( Found ) THEN
        DO p = nd,nd
          DO q = 1,nd-np-1
            j = np + q
            Stiff(p,j) = Stiff(p,j) - SUM(FaceBasis(q,1:dim) * v(1:dim)) * 1.0d0 * s
          END DO
        END DO
      END IF

      ! Contribution of source term
      !-----------------------------------------------------
      f = ListGetElementReal( SourceField_h, Basis, Element, Found, GaussPoint = t )      
      IF ( Found ) THEN
        DO p = nd,nd
          Force(p) = Force(p) - f * 1.0d0 * s 
        END DO
      END IF

      IF (NeedMass) THEN
        !
        ! Here piecewise constant approximation is assumed. This
        ! could be integrated with one-point quadrature.
        !
        Mass(nd,nd) = Mass(nd,nd) - s
      END IF

    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, Mesh, n, nd, SecondFamily, InitHandles)
!------------------------------------------------------------------------------
    USE ElementDescription, ONLY : PickActiveFace
    IMPLICIT NONE

    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n   ! The number of background element nodes
    INTEGER :: nd  ! The total count of DOFs (nodal and facial)
    LOGICAL :: SecondFamily
    LOGICAL :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Parent, Face
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: ReverseSign(6), stat, AssembleForce, OrientationsMatch
    LOGICAL :: RevertIndices
    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET :: TetraFaceMap(4,3), BrickFaceMap(6,4)
    INTEGER :: FDofMap(6,4)
    INTEGER :: OriginalIndices(4)
    INTEGER :: j, k, l, t, np, p, ActiveFaceId, Family, FDOFs
    INTEGER :: ParentFamily
    REAL(KIND=dp) :: Force(nd), Basis(n), TraceBasis(nd), WorkTrace(nd)
    REAL(KIND=dp) :: detJ, s, u, v, w, g
    TYPE(ValueHandle_t), SAVE :: ScalarField_h
 
    SAVE Nodes
!------------------------------------------------------------------------------
    Family = GetElementFamily(Element)
    IF (Family == 1) RETURN

    BC => GetBC()
    IF (.NOT. ASSOCIATED(BC)) RETURN

    IF (InitHandles) THEN
      CALL ListInitElementKeyword(ScalarField_h, 'Boundary Condition', &
          'Scalar Field')
      InitHandles = .FALSE.
    END IF
    IF (ScalarField_h % NotPresentAnywhere) RETURN

    ! 
    ! The sign reversion of basis will be checked via the parent element:
    ! 
    Parent => Element % BoundaryInfo % Left
    IF (.NOT. ASSOCIATED(Parent)) THEN
      Parent => Element % BoundaryInfo % Right
    END IF
    IF (.NOT. ASSOCIATED(Parent)) RETURN
    ParentFamily = GetElementFamily(Parent)
    !
    ! Identify the face representing the element among the faces of 
    ! the parent element:
    !
    CALL PickActiveFace(Mesh, Parent, Element, Face, ActiveFaceId)
    IF (ActiveFaceId == 0) RETURN
    FDOFs = Face % BDOFs
    IF (FDOFs < 1) RETURN
    !
    ! Use the parent element to check whether sign reversions are needed:
    !
    CALL FaceElementOrientation(Parent, ReverseSign, ActiveFaceId)

    !
    ! In the case of the basis of the second kind the effect of orientation
    ! must be taken into account:
    !
    RevertIndices = .FALSE.
    IF ((ParentFamily == 3 .OR. ParentFamily == 5) .AND. SecondFamily) THEN
      SELECT CASE(Family)
      CASE(2)
        !
        ! Check whether the parametrization of the element conforms with the global positive 
        ! orientation of the edge:
        !
        FaceMap => GetEdgeMap(GetElementFamily(Parent))
        IF (ReverseSign(ActiveFaceId)) THEN
          OrientationsMatch = ( Element % NodeIndexes(1) == Parent % NodeIndexes(FaceMap(ActiveFaceId,2)) )
        ELSE
          OrientationsMatch = ( Element % NodeIndexes(1) == Parent % NodeIndexes(FaceMap(ActiveFaceId,1)) )
        END IF

      CASE(3)
        IF (FDOFs /= 3) CALL Fatal('ModelMixedPoisson', '3-DOF faces expected')
        TetraFaceMap(1,:) = (/ 2, 1, 3 /)
        TetraFaceMap(2,:) = (/ 1, 2, 4 /)
        TetraFaceMap(3,:) = (/ 2, 3, 4 /) 
        TetraFaceMap(4,:) = (/ 3, 1, 4 /)

        !FaceMap => TetraFaceMap 

        CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId)

        IF (ANY(Element % NodeIndexes(1:3) /= Parent % NodeIndexes(TetraFaceMap(ActiveFaceId,1:3)))) THEN
          !
          ! The parent element face is indexed differently, reorder and revert afterwards:
          !
          OriginalIndices(1:3) = Element % NodeIndexes(1:3)
          Element % NodeIndexes(1:3) = Parent % NodeIndexes(TetraFaceMap(ActiveFaceId,1:3))
          RevertIndices = .TRUE.
        END IF

      END SELECT
    ELSE
      SELECT CASE(GetElementFamily(Parent))
      CASE(8)
        IF (FDOFs /= 4) CALL Fatal('ModelMixedPoisson', '4-DOF faces expected')
 
        BrickFaceMap(1,:) = (/ 2, 1, 4, 3 /)
        BrickFaceMap(2,:) = (/ 5, 6, 7, 8 /)
        BrickFaceMap(3,:) = (/ 1, 2, 6, 5 /)
        BrickFaceMap(4,:) = (/ 2, 3, 7, 6 /)
        BrickFaceMap(5,:) = (/ 3, 4, 8, 7 /)
        BrickFaceMap(6,:) = (/ 4, 1, 5, 8 /)

        CALL FaceElementBasisOrdering(Parent, FDofMap, ActiveFaceId)
        
        IF (ANY(Element % NodeIndexes(1:4) /= Parent % NodeIndexes(BrickFaceMap(ActiveFaceId,1:4)))) THEN
          !
          ! The parent element face is indexed differently, reorder and revert afterwards:
          !
          OriginalIndices(1:4) = Element % NodeIndexes(1:4)
          Element % NodeIndexes(1:4) = Parent % NodeIndexes(BrickFaceMap(ActiveFaceId,1:4))
          RevertIndices = .TRUE.        
        END IF
      END SELECT
    END IF

    np = n * Solver % Def_Dofs(GetElementFamily(Parent), Parent % BodyId, 1)

    IF (ReverseSign(ActiveFaceId)) THEN
      s = -1.0d0
    ELSE
      s = 1.0d0
    END IF

    CALL GetElementNodes(Nodes)
    IF (Family == 3) THEN
      !
      ! Integration must be done over p-reference element
      !
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    ELSE
      IP = GaussPoints(Element)
    END IF

    Force = 0.0d0
    DO t=1,IP % n
      !
      ! NOTE: Here the effect of the Piola transformation is taken into account
      !       such that the multiplication with DetJ is not needed
      ! TO CONSIDER: Get the traces of vector-values basis functions 
      !              by calling a subroutine
      !
      SELECT CASE(Family)
      CASE(2)
        !--------------------------------------------------------------
        ! Basis function values at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), DetJ, Basis)
        IF (SecondFamily) THEN
          u = IP % U(t)
          IF (OrientationsMatch) THEN
            TraceBasis(1) = s * 0.5d0 * (1.0d0 - sqrt(3.0d0)*u)
            TraceBasis(2) = s * 0.5d0 * (1.0d0 + sqrt(3.0d0)*u)
          ELSE
            TraceBasis(1) = s * 0.5d0 * (1.0d0 + sqrt(3.0d0)*u)
            TraceBasis(2) = s * 0.5d0 * (1.0d0 - sqrt(3.0d0)*u)
          END IF
        ELSE
          TraceBasis(1) = s * 0.5d0
        END IF
      CASE(3)
        u = IP % U(t)
        v = IP % V(t)
        DO k=1,n
          Basis(k) = TriangleNodalPBasis(k, u, v)
        END DO

        IF (SecondFamily) THEN
          WorkTrace(1) = -1.0d0*(-8 - 12*u + 4*Sqrt(3.0d0)*v)/12.0d0
          WorkTrace(2) = -1.0d0*(u + (-8 + 4*Sqrt(3.0d0)*v)/12.0d0)
          WorkTrace(3) = -1.0d0*(4.0d0 - 8*Sqrt(3.0d0)*v)/12.0d0
          !
          ! Reorder and reverse signs:
          !
          DO j=1,FDOFs
            k = FDofMap(ActiveFaceId,j)
            TraceBasis(j) = s * WorkTrace(k)
          END DO
        ELSE
          TraceBasis(1) = s / SQRT(3.0d0)
        END IF
      CASE(4)
        stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), DetJ, Basis)

        DO j=1,FDOFs
          k = FDofMap(ActiveFaceId,j)
          TraceBasis(j) = s * Basis(k)
        END DO
      END SELECT

      w = IP % s(t) ! NOTE: No need to multiply with DetJ
      g = ListGetElementReal(ScalarField_h, Basis, Element, AssembleForce)

      IF (AssembleForce) THEN
        DO p = 1,nd-np
          j = np + p
          Force(j) = Force(j) + g * TraceBasis(p) * w 
        END DO
      END IF
    END DO

    IF (AssembleForce) CALL DefaultUpdateForce(Force)

    IF (RevertIndices) THEN
      SELECT CASE(Family)
      CASE(3)
        Element % NodeIndexes(1:3) = OriginalIndices(1:3)
      CASE(4)
        Element % NodeIndexes(1:4) = OriginalIndices(1:4)
      END SELECT
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE MixedPoisson
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Postprocessing utility for the main solver. 
!------------------------------------------------------------------------------
SUBROUTINE MixedPoisson_post(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params 
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  LOGICAL :: Found
  TYPE(Variable_t), POINTER :: pVar, Var, fVar
  TYPE(Element_t), POINTER :: Element
  INTEGER, ALLOCATABLE :: Indexes(:)
  INTEGER :: dim, active, t, n, nd, nb, np
  REAL(KIND=dp) :: val

  LOGICAL :: SecondFamily

  REAL(KIND=dp), ALLOCATABLE :: Flux_x(:), Flux_y(:), Flux_z(:)

  Params => GetSolverParams()

  SecondFamily = GetLogical(Params, 'Second Kind Basis', Found)
  
  Mesh => GetMesh()
  Var => Solver % Variable
  dim = CoordinateSystemDimension()

  n = Solver % Mesh % MaxElementDOFs
  ALLOCATE( Indexes(n) )
  ALLOCATE( Flux_x(n), Flux_y(n), Flux_z(n))
  
  ! Get the elemental pressure variable where postprocessing is saved to
  VarName = ListGetString(Params,'Potential Variable')
  pVar => VariableGet( Mesh % Variables, VarName ) 

  VarName = ListGetString(Params,'Flux Variable')
  fVar => VariableGet( Mesh % Variables, VarName ) 
  
  active = GetNOFActive()
  
  DO t=1,active
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() 
    nd = GetElementDOFs( Indexes, Element )  
    nb = SIZE(Element % BubbleIndexes(:))

    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)    

    ! last bubble dofs is the pressure
    val = Var % Values( Var % Perm(Indexes(nd)) )

    ! assign it to be outputted
    pVar % Values(pVar % Perm(Element % ElementIndex)) = val

    CALL GetFlux( n,nd-np-1,Var % Values(Var % Perm(Indexes(np+1:nd-1))) )

    fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+1) = Flux_x(1:n)
    fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+2) = Flux_y(1:n)
    fVar % Values(3*(fVar % Perm(Element % dGIndexes(1:n))-1)+3) = Flux_z(1:n)

  END DO

CONTAINS

  SUBROUTINE GetFlux( n, nval,vals )

    INTEGER :: n, nval
    REAL(KIND=dp) :: vals(:)
    !------------------------------------------------------------------------   
    REAL(KIND=dp) :: MASS(n,n), FORCE_x(n), FORCE_y(n), FORCE_z(n)
    INTEGER, PARAMETER :: MaxFaceBasisDim = 48
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np
    REAL(KIND=dp) :: FaceBasis(MaxFaceBasisDim,3), DivFaceBasis(MaxFaceBasisDim)
    REAL(KIND=dp) :: Basis(n), DetJ, s


    !------------------------------------------------------------------------
    ! The reference element is chosen to be that used for p-approximation,
    ! so we need to switch to using a quadrature which would not be used 
    ! otherwise
    !------------------------------------------------------------------------
    SELECT CASE( GetElementFamily(Element) )
    CASE(3)
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    CASE(4)
      IP = GaussPointsQuad(9)
    CASE(5)
      IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
    CASE(8)
      IP = GaussPointsBrick(64)
    CASE DEFAULT
      CALL Fatal('MixedPoisson', 'A non-supported element type')
    END SELECT

    CALL GetElementNodes(Nodes)

    MASS = 0._dp
    FORCE_x = 0._dp
    FORCE_y = 0._dp
    FORCE_z = 0._dp
    ! Set np = n, if nodal dofs are employed; otherwise set np = 0:
    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1)    

    DO t=1,IP % n
      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)

      s = IP % s(t) * DetJ

      DO i=1,n
        DO j=1,n
          MASS(i,j) = MASS(i,j) + s * Basis(i)*Basis(j)
        END DO
        FORCE_x(i) = FORCE_x(i) + s * SUM(FaceBasis(1:nval,1) * vals(1:nval)) * Basis(i)
        FORCE_y(i) = FORCE_y(i) + s * SUM(FaceBasis(1:nval,2) * vals(1:nval)) * Basis(i)
        FORCE_z(i) = FORCE_z(i) + s * SUM(FaceBasis(1:nval,3) * vals(1:nval)) * Basis(i)
      END DO
    END DO

    CALL InvertMatrix(MASS,n)
    Flux_x(1:n) = MATMUL(MASS, FORCE_x)
    Flux_y(1:n) = MATMUL(MASS, FORCE_y)
    Flux_z(1:n) = MATMUL(MASS, FORCE_z)

  END SUBROUTINE GetFlux


!------------------------------------------------------------------------------        
END SUBROUTINE MixedPoisson_post
!-----------------------------------------------------------------------------
