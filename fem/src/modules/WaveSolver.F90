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
! *  Solvers: WaveSolver 
! *  Authors: Juha Ruokolainen, Peter Råback, Mika Malinen
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: ~2013
! *
! *****************************************************************************/


!------------------------------------------------------------------------------
SUBROUTINE WaveSolver_Init(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverParams
!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  CALL ListAddInteger( SolverParams, 'Time derivative order', 2 )

!------------------------------------------------------------------------------
END SUBROUTINE WaveSolver_Init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solves the transient wave equation for a scalar variable using H1-conforming 
!> basis functions and the Galerkin method.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE WaveSolver(Model, Solver, dt, TransientSimulation)
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
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm, Pave
  INTEGER :: dim, maxiter, iter, n, nb, nd, t, active
  LOGICAL :: Found, InitHandles
!------------------------------------------------------------------------------
  
  dim = CoordinateSystemDimension()
  maxiter = ListGetInteger(GetSolverParams(), &
      'Nonlinear System Max Iterations', Found, minv=1)
  IF(.NOT. Found ) maxiter = 1

  CALL DefaultStart()

  DO iter=1,maxiter
    !----------------
    !System assembly:
    !----------------
    CALL DefaultInitialize()
    InitHandles = .TRUE.
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element)
      nb = GetElementNOFBDOFs(Element)
      CALL LocalMatrix(Element, n, nd+nb, dim, InitHandles)
    END DO

    CALL DefaultFinishBulkAssembly()

    InitHandles = .TRUE.
    Active = GetNOFBoundaryElements()
    DO t=1,Active
      Element => GetBoundaryElement(t)
      IF (ActiveBoundaryElement()) THEN
        n  = GetElementNOFNodes(Element)
        nd = GetElementNOFDOFs(Element)
        CALL LocalMatrixBC(Element, n, nd, InitHandles)
      END IF
    END DO

    CALL DefaultFinishBoundaryAssembly()

    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()


    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

    IF( GetLogical( Solver % Values, 'Set Average To Zero', Found ) ) THEN
      Pave = SUM( Solver % Variable % Values) / &
          SIZE( Solver % Variable % Values ) 
      Solver % Variable % Values = Solver % Variable % Values - Pave
    END IF

  END DO

  CALL DefaultFinish()

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(Element, n, nd, dim, InitHandles)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd, dim
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueHandle_t) :: Load_h, SoundSpeed_h, DampingCoeff_h, ReactCoeff_h
    TYPE(Nodes_t) :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    LOGICAL :: Stat, Found, AssembleSource, DampingActive, ReactiveMedium

    INTEGER :: i, t, p, q

    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd), DAMP(nd,nd)
    REAL(KIND=dp) :: c, att, react, LoadAtIP
    REAL(KIND=dp) :: Weight, Basis(nd), dBasisdx(nd,3), DetJ

    SAVE Load_h, SoundSpeed_h, DampingCoeff_h, ReactCoeff_h, Nodes
!------------------------------------------------------------------------------
    IF (InitHandles) THEN
      CALL ListInitElementKeyword(Load_h, 'Body Force', 'Sound source')
      CALL ListInitElementKeyword(SoundSpeed_h, 'Material', 'Sound speed', &
          UnfoundFatal=.TRUE.)
      CALL ListInitElementKeyword(DampingCoeff_h, 'Material', 'Sound damping')
      CALL ListInitElementKeyword(ReactCoeff_h, 'Material', &
          'Sound reaction damping')
      InitHandles = .FALSE.
    END IF

    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    DAMP = 0._dp

    AssembleSource = .FALSE.

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      c = ListGetElementReal(SoundSpeed_h, Basis, Element)
      att = ListGetElementReal(DampingCoeff_h, Basis, Element, DampingActive)
      react = ListGetElementReal(ReactCoeff_h, Basis, Element, ReactiveMedium)

      ! TO DO: Source is now a scalar field div(b). Rather than giving div(b)
      ! it would be better to give the vector b and to apply
      ! integration by parts to get always consistent flux BCs.
      IF (.NOT. Load_h % NotPresentAnywhere) &
          LoadAtIP = ListGetElementReal(Load_h, Basis, Element, AssembleSource) 

      Weight = IP % s(t) * DetJ

      ! The Laplace term:
      ! -----------------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) - Weight * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      ! Damping term to implement "viscous wave equation":
      !---------------------------------------------------
      IF (DampingActive) Damp(1:nd,1:nd) = Damp(1:nd,1:nd) - Weight * &
          (att/c**2) * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

      DO p=1,nd
        DO q=1,nd

          ! Reaction term:
          ! -----------------------------------
          IF (ReactiveMedium) Damp(p,q) = Damp(p,q) - Weight * &
              (react/c**2) * Basis(q) * Basis(p)

          ! The 2nd time derivative:
          ! ------------------------------
          MASS(p,q) = MASS(p,q) - Weight * &
              1.0_dp/ c**2  * Basis(q) * Basis(p)

        END DO
      END DO

      IF (AssembleSource) &
          FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF( TransientSimulation) THEN
       CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
    END IF

    ! Applying static condensation is a risky endeavour since the values of 
    ! the bubble DOFs at previous time steps are not recovered, disable the 
    ! static condensation?
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF, FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC(Element, n, nd, InitHandles)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    INTEGER, INTENT(IN) :: n, nd
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    TYPE(ValueHandle_t) :: Flux_h, SoundSpeed_h
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: AssembleFlux, OutflowBC, Stat, Found
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), DAMP(nd,nd), MASS(nd,nd)
    REAL(KIND=dp) :: c, g
    REAL(KIND=dp) :: Weight, Basis(nd), dBasisdx(nd,3), DetJ
    INTEGER :: t, p, q
 
    SAVE Flux_h, SoundSpeed_h, Nodes
!------------------------------------------------------------------------------
    BC => GetBC()
    IF (.NOT.ASSOCIATED(BC)) RETURN

    IF (InitHandles) THEN
      CALL ListInitElementKeyword(Flux_h, 'Boundary Condition', &
          'Source Acceleration')
      CALL ListInitElementKeyword(SoundSpeed_h, 'Material', 'Sound speed', &
          UnfoundFatal=.TRUE.)
      InitHandles = .FALSE.
    END IF

    OutflowBC = GetLogical(BC, 'Plane Wave BC', Found)
    IF (.NOT. Found) OutflowBC = GetLogical(BC, 'Outflow Boundary', Found)
    IF (.NOT. OutflowBC .AND. Flux_h % NotPresentAnywhere) RETURN

    CALL GetElementNodes( Nodes )

    STIFF = 0._dp
    FORCE = 0._dp
    DAMP = 0._dp
    MASS = 0._dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx)

      Weight = IP % s(t) * DetJ

      IF (OutflowBC) THEN
        c = ListGetElementRealParent(SoundSpeed_h, Basis, Element)
        DO p=1,nd
          DO q=1,nd
            DAMP(p,q) = DAMP(p,q) - Weight * Basis(q) * Basis(p) / c
          END DO
        END DO
      END IF

      g = ListGetElementReal(Flux_h, Basis, Element, AssembleFlux)
      IF (AssembleFlux) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * g * Basis(1:nd)
      END IF

    END DO

    IF( TransientSimulation) THEN
       CALL Default2ndOrderTime( MASS, DAMP, STIFF, FORCE )
    END IF
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE WaveSolver
!------------------------------------------------------------------------------
