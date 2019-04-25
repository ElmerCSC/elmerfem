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
! * Compute the L2 errors of the mixed solution of the Poisson equation
! * in a special case.
! *
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 12, 2019
! *
!******************************************************************************

!------------------------------------------------------------------------------
SUBROUTINE IntegrateError(Model, Solver, dt, TransientSimulation)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: SolverPars, TargetSolverPars
  TYPE(Solver_t), POINTER :: TargetSolver
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Element_t),POINTER :: Element

  LOGICAL :: AllocationsDone = .FALSE., Found, SecondFamily

  CHARACTER(LEN=MAX_NAME_LEN) :: TargetEq, Eq
 
  INTEGER, ALLOCATABLE :: Indices(:)  
  INTEGER :: i, istat, nlen
  INTEGER :: n, nb, np, nd, t, dim

  REAL(KIND=dp), ALLOCATABLE :: SolDOFs(:)
  REAL(KIND=dp) :: PresErr, FluxErr, EK(2), SolNorm, FluxNorm

  SAVE SolDOFs, Indices, AllocationsDone
!------------------------------------------------------------------------------
  SolverPars => GetSolverParams()

  !
  ! Get pointer to the solver the variable of which is used to integrate 
  ! the error:
  !

  TargetEq = ListGetString(SolverPars, 'Target Equation', UnfoundFatal=.TRUE.)
  nlen = LEN_TRIM(TargetEq)
  TargetSolver => NULL()
  DO i=1,Model % NumberOfSolvers
    Eq = ListGetString(Model % Solvers(i) % Values, 'Equation', Found)
    IF (Found) THEN
      IF (nlen == LEN_TRIM(Eq)) THEN
        IF (TargetEq(1:nlen) == Eq(1:nlen)) TargetSolver => Model % Solvers(i)
      END IF
    END IF
  END DO
  IF (.NOT. ASSOCIATED(TargetSolver)) CALL Fatal('ComputeError', &
      'Target Equation does not exist')

  TargetSolverPars => GetSolverParams(TargetSolver)
  SecondFamily = GetLogical(TargetSolverPars, 'Second Kind Basis', Found)
  dim = CoordinateSystemDimension()


  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()
  IF ( .NOT. AllocationsDone ) THEN
     N = Mesh % MaxElementDOFs
     ALLOCATE(SolDOFs(N), Indices(N), STAT=istat)
     IF (istat /= 0) THEN
        CALL Fatal( 'ComputeError', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  
  PresErr = 0.0d0
  FluxErr = 0.0d0
  SolNorm = 0.0d0
  FluxNorm = 0.0d0
  DO t=1,TargetSolver % NumberOfActiveElements

    Element => GetActiveElement(t, TargetSolver)
    n  = GetElementNOFNodes(Element)
    nd = GetElementDOFs(Indices, Element, TargetSolver)
    nb = size(Element % BubbleIndexes(:))

    np = n * TargetSolver % Def_Dofs(GetElementFamily(Element), &
        Element % BodyId, 1)

    DO i=1,nd
      SolDOFs(i) = TargetSolver % Variable % Values(TargetSolver % Variable % &
          Perm(Indices(i)))
    END DO

    CALL ElementwiseError(SolDofs, Element, n, nd, nb, np, dim, EK, SolNorm, &
        FluxNorm, SecondFamily)

    PresErr = PresErr + EK(1)
    FluxErr = FluxErr + EK(2)

  END DO
  WRITE (*, '(A,E16.8)') 'L2 Scalar Error = ', SQRT(ParallelReduction(PresErr))/SQRT(ParallelReduction(SolNorm))
  WRITE (*, '(A,E16.8)') 'L2 Flux Error =', SQRT(ParallelReduction(FluxErr))/SQRT(ParallelReduction(FluxNorm))

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE ElementwiseError(SolDOFs, Element, n, nd, nb, np, dim, EK, &
      SolNorm, FluxNorm, SecondFamily)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: SolDOFs(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n, nd, nb, np, dim
    REAL(KIND=dp) :: EK(2), SolNorm, FluxNorm
    LOGICAL :: SecondFamily
!------------------------------------------------------------------
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: i, p, t
    LOGICAL :: Stat
    REAL(KIND=dp) :: FaceBasis(nd-nb,3), DivFaceBasis(nd-nb), Basis(n)
    REAL(KIND=dp) :: uq, vq, wq, sq
    REAL(KIND=dp) :: xq, yq, zq, DetJ, sol, gradsol(dim), fluxsol(dim)


!---------------------------------------------------------------------------
    IF (nb /= 1) CALL Fatal('ComputeError', 'One bubble DOF assumed')

    CALL GetElementNodes( Nodes )

    EK = 0.0d0

    SELECT CASE( GetElementFamily(Element) )
    CASE(3)
      IP = GaussPointsTriangle(3, PReferenceElement=.TRUE.)
    CASE(5)
      IP = GaussPointsTetra(4, PReferenceElement=.TRUE.)
    CASE DEFAULT
      CALL Fatal('ComputeError', 'A simplicial mesh assumed currently')
    END SELECT 

    DO t=1,IP % n
       stat = FaceElementInfo(Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), detF=detJ, Basis=Basis, FBasis=FaceBasis, &
           DivFBasis=DivFaceBasis, BDM = SecondFamily, ApplyPiolaTransform=.TRUE.)
       
       xq = SUM( Nodes % x(1:n) * Basis(1:n) )
       yq = SUM( Nodes % y(1:n) * Basis(1:n) )
       zq = SUM( Nodes % z(1:n) * Basis(1:n) )      
    
       sol = -64.0d0*(xq**2-1.0d0/4.0d0) * (yq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(1) = -64.0d0 * 2.0d0 * xq * (yq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(2) = -64.0d0 * 2.0d0 * yq * (xq**2-1.0d0/4.0d0) * (zq**2-1.0d0/4.0d0)
       gradsol(3) = -64.0d0 * 2.0d0 * zq * (xq**2-1.0d0/4.0d0) * (yq**2-1.0d0/4.0d0)       

       ! For computing the error of the  scalar field:
       SolNorm = SolNorm + sol**2 * detJ * IP % s(t)
       EK(1) = EK(1) + (sol - SolDOFs(nd))**2 * detJ * IP % s(t) 


       FluxNorm = FluxNorm + SUM( gradsol(1:3)*gradsol(1:3) ) * detJ * IP % s(t)

       FluxSol = 0.0d0
       DO p = 1,nd-np-nb
         i = np + p
         FluxSol(1:dim) = FluxSol(1:dim) + FaceBasis(p,1:dim) * SolDOFs(i)
       END DO

       EK(2) = EK(2) + SUM( (gradsol(1:dim) - fluxsol(1:dim))**2 ) * &
            detJ * IP % s(t) 
          
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ElementwiseError
!-----------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE IntegrateError
!------------------------------------------------------------------------------
