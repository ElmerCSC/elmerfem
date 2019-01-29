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
! *  Authors: Mika Malinen
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: Jan 21, 2019
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
SUBROUTINE BestApproximationSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!
!  Solve the best approximation of the vector field
! 
!              U = (1,1) in 2D, or          (Test Mode = Integer 1)   
!              U = (1,1,1) in 3D            (Test Mode = Integer 2)
!
!  with respect to the L2 norm (the default) or an energy norm using 
!  H(div)-conforming basis functions. Here the energy norm corresponds to 
!  the operator I + MatPar * div grad, with MatPar a scalar field specified
!  by the user. Additionally, compute the relative error of the solution or 
!  of the div field using the L2 norm. This solver can thus be used for checking
!  that the convergence rate is correct.
!
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
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t), POINTER :: BodyForce, Material, BC
  TYPE(Variable_t), POINTER :: Var

  REAL(KIND=dp) :: Norm, u, v, w, Err, EK, SolNorm, MatPar = 0.0d0

  INTEGER :: n, ne, nf, nb, np, nd, t, istat, i, j, k, l, active, dim, TestMode
  INTEGER :: ElementOrder

  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:)
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Matrix_t), POINTER :: A

  LOGICAL :: stat, SecondFamily, ErrorEstimation
  LOGICAL :: UseDivNorm

  INTEGER, ALLOCATABLE :: Indices(:)

  SAVE STIFF, LOAD, FORCE, Acoef, AllocationsDone, Nodes, Indices
!------------------------------------------------------------------------------

  ElementOrder = 1
  SecondFamily = GetLogical( GetSolverParams(), 'Second Kind Basis', Found)

  ErrorEstimation = GetLogical( GetSolverParams(), 'Error Computation', Found)
  IF (.NOT. Found) ErrorEstimation = .TRUE.

  UseDivNorm = GetLogical( GetSolverParams(), 'Compute Div Error', Found)

  TestMode = ListGetInteger( GetSolverParams(), 'Test Mode', Found)
  IF (.NOT. Found) TestMode = 1

  dim = CoordinateSystemDimension()

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs  ! just big enough
    ALLOCATE( FORCE(N), LOAD(6,N), STIFF(N,N), &
        Indices(N), Acoef(N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'BestApproximationSolver', 'Memory allocation error.' )
    END IF
    AllocationsDone = .TRUE.
  END IF
  
  Solver % Matrix % COMPLEX = .FALSE.
  A => GetMatrix()

  !-----------------------
  ! System assembly:
  !----------------------
  active = GetNOFActive()
  CALL DefaultInitialize()

  DO t=1,active

    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() ! The number of nodes corresponding to the background mesh
    nd = GetElementNOFDOFs()
    nb = GetElementNOFBDOFs()

    Acoef(1:n) = 0.0d0
    Material => GetMaterial( Element )
    IF ( ASSOCIATED(Material) ) THEN
      Acoef(1:n) = GetReal( Material, 'Material Param', Found )
      IF (.NOT. Found) Acoef(1:n) = 0.0d0
    END IF

    !Get element local matrix and rhs vector:
    !----------------------------------------
    CALL LocalMatrix(STIFF, FORCE, Acoef, Element, n, nd+nb, dim)

    !Update global matrix and rhs vector from local matrix & vector:
    !---------------------------------------------------------------
    CALL DefaultUpdateEquations( STIFF, FORCE )

  END DO

  CALL DefaultFinishBulkAssembly()
  CALL DefaultFinishAssembly()

  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()


  !-------------------------------------------------------------------
  ! Compute the norm of the error
  !--------------------------------------------------------------------
  IF (ErrorEstimation) THEN
    Err = 0.0d0
    SolNorm = 0.0d0
    DO t=1,Solver % NumberOfActiveElements

      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementDOFs( Indices )

      Load(1,1:nd) = Solver % Variable % Values( Solver % Variable % &
          Perm(Indices(1:nd)) )

      CALL ComputeApproximationError(Load, Element, n, nd, dim, Err, SolNorm, UseDivNorm)
    END DO

    WRITE (*, '(A,E16.8)') 'Error Norm = ', SQRT(ParallelReduction(Err))/SQRT(ParallelReduction(SolNorm))
  END IF

CONTAINS


!---------------------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(STIFF, FORCE, NodalMatPar, Element, n, nd, dim)
!---------------------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), NodalMatPar(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, dim
    !------------------------------------------------------------------------------
    REAL(KIND=dp) :: MatPar
    REAL(KIND=dp) :: FaceBasis(nd,3), DivFaceBasis(nd)
!    REAL(KIND=dp) :: dBasisdx(n,3)
    REAL(KIND=dp) :: Basis(n), DetJ, xq, yq, zq, uq, vq, wq, sq
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np

    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    STIFF = 0.0d0
    FORCE = 0.0d0

    !------------------------------------------------------------------------------
    ! Numerical integration over element: Use the same quadrature which 
    ! would be used for the curl-conforming elemenents
    !------------------------------------------------------------------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE., &
         EdgeBasisDegree=ElementOrder)    
    
    ! Set np = n, if nodal dofs are employed; otherwise set np = 0:
    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1) 

    DO t=1,IP % n

      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)

      xq = SUM( Nodes % x(1:n) * Basis(1:n) )
      yq = SUM( Nodes % y(1:n) * Basis(1:n) )
      zq = SUM( Nodes % z(1:n) * Basis(1:n) )

      MatPar = SUM( NodalMatPar(1:n) * Basis(1:n) )

      !----------------------------------------------------------------
      ! The following branch could be used to produce the 
      ! Galerkin projection of a solution component for visualization.
      !------------------------------------------------------------------
      IF (np > 0) THEN
        DO p = 1,n
          DO q = 1,n       
            STIFF(p,q) = STIFF(p,q) + Basis(p) * Basis(q) * detJ * IP % s(t)    
          END DO

          DO q = 1,nd-np
            j = np + q
            STIFF(p,j) = STIFF(p,j) - FaceBasis(q,1) * Basis(p) * detJ * IP % s(t)
          END DO
        END DO
      END IF

      !--------------------------------------------------------------
      ! The equations for H(div)-conforming part...
      !---------------------------------------------------------------
      DO p = 1,nd-np
        !----------------------------
        ! The inner product (u,v)_E
        !----------------------------
        i = np + p
        DO q = 1,nd-np
          j = np + q
          STIFF(i,j) = STIFF(i,j) + 1.0d0 * &
              SUM( FaceBasis(q,1:dim) * FaceBasis(p,1:dim) ) * detJ * IP % s(t) + &
              MatPar * DivFaceBasis(q) * DivFaceBasis(p) * detJ * IP % s(t)
        END DO

        !----------------------------------------
        ! RHS corresponding to the exact solution 
        !----------------------------------------
        SELECT CASE(TestMode)
        CASE (1)
          FORCE(i) = FORCE(i) +  (1.0d0) * FaceBasis(p,1) * detJ * IP % s(t) + &
              (1.0d0)* FaceBasis(p,2) * detJ * IP % s(t)
        CASE (2)
          FORCE(i) = FORCE(i) +  (1.0d0) * FaceBasis(p,1) * detJ * IP % s(t) + &
              (1.0d0) * FaceBasis(p,2) * detJ * IP % s(t) + &
              (1.0d0) * FaceBasis(p,3) * detJ * IP % s(t)
        CASE DEFAULT
          CALL Fatal('FaceFEInterpolationTest', 'This Test Mode needs an update')
        END SELECT
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!----------------------------------------------------------------------------------
  SUBROUTINE ComputeApproximationError(Load, Element, n, nd, dim, EK, SolNorm, UseDivNorm)
!----------------------------------------------------------------------------------
    REAL(KIND=dp) :: Load(:,:), EK, SolNorm
    TYPE(Element_t), POINTER :: Element    
    INTEGER :: n, nd, dim
    LOGICAL :: UseDivNorm
!--------------------------------------------------------------------------------
    REAL(KIND=dp) :: FaceBasis(nd,3), DivFaceBasis(nd)
    REAL(KIND=dp) :: Basis(n), DetJ, xq, yq, zq, uq, vq, wq, sq, &
         u(3), divu, sol(3), divsol, e(3), dive
!    REAL(KIND=dp) :: dBasisdx(n,3)
    LOGICAL :: Stat, Found
    INTEGER :: t, i, j, p, q, np
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )

    !-------------------------------------
    ! Numerical integration over element:
    !-------------------------------------
    IP = GaussPoints(Element, EdgeBasis=.TRUE., PReferenceElement=.TRUE., EdgeBasisDegree=ElementOrder) 

    ! Set np = n, if nodal dofs are employed; otherwise set np = 0:
    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyId, 1) 

    DO t=1,IP % n

      stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
          DivFBasis=DivFaceBasis, BDM=SecondFamily, ApplyPiolaTransform=.TRUE.)

      xq = SUM( Nodes % x(1:n) * Basis(1:n) )
      yq = SUM( Nodes % y(1:n) * Basis(1:n) )
      zq = SUM( Nodes % z(1:n) * Basis(1:n) )

      u = 0.0d0
      DO i=1,dim
        u(i) = SUM( Load(1,np+1:nd) * FaceBasis(1:nd-np,i) )
      END DO

      divu = SUM( Load(1,np+1:nd) * DivFaceBasis(1:nd-np) )

      ! Compute the square of the energy norm of the solution and error:
      SELECT CASE(TestMode)
      CASE (1)
        sol(1) = 1.0d0
        sol(2) = 1.0d0
        sol(3) = 0.0d0
        divsol = 0.0d0
      CASE (2)
        sol(1) = 1.0d0
        sol(2) = 1.0d0
        sol(3) = 1.0d0
        divsol = 0.0d0
      CASE DEFAULT
        CALL Fatal('FaceFEInterpolationTest', 'This Test Mode needs an update')
      END SELECT

      e = 0.0d0
      e(1:dim) = sol(1:dim) - u(1:dim)  
      dive = divsol - divu

      IF (UseDivNorm .AND. (TestMode>2)) THEN
        ! Div error in L2:
        !-------------------
        SolNorm = SolNorm + divsol * divsol * detJ
        EK = EK + dive * dive * detJ       
 
        ! Energy norm:
        !--------------
        !SolNorm = SolNorm + (SUM( Sol(1:3) * Sol(1:3) ) + 1.0d0 * divsol * divsol)  * detJ
        !EK = EK + (SUM( e(1:3) * e(1:3) ) + 1.0d0 * dive * dive) * detJ

      ELSE
        ! L2 norm
        SolNorm = SolNorm + SUM( Sol(1:3) * Sol(1:3) ) * detJ
        EK = EK + SUM( e(1:3) * e(1:3) ) * detJ
      END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeApproximationError
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE BestApproximationSolver
!------------------------------------------------------------------------------
