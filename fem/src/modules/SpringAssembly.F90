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
! * Add node-wise specified springs to structural models. This solver can be
! * called as an additional assembly solver by using the keyword "Assembly 
! * Solvers" in the solver section associated with the model which is modified
! * to have the springs. Element-wise spring constraints over higher-
! * dimensional entities are handled in the primary solver code. This solver
! * assumes that the places of the springs are listed using the keyword
! * "Target Nodes".
! *
! *
! *  Authors: Mika Malinen
! *  Email:   mika.malinen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: March 25, 2019
! *
!******************************************************************************

!------------------------------------------------------------------------------
SUBROUTINE SpringAssembler(Model, Solver, dt, TransientSimulation)
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
  TYPE(ValueList_t), POINTER :: SolverParams, ValueList
  TYPE(Solver_t), POINTER :: SolverPtr
  TYPE(Element_t), TARGET :: Element

  LOGICAL, ALLOCATABLE :: Assemble(:) 
  LOGICAL :: Found, NodesFound

  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: BC, DOFs, i, istat, j, n

  REAL(KIND=dp), ALLOCATABLE :: K(:,:), f(:), Work(:,:)

  CHARACTER(LEN=MAX_NAME_LEN) :: UName
!------------------------------------------------------------------------------
  SolverParams => GetSolverParams()

  !
  ! Find the primary solver whose matrix should be modified to include springs:
  !
  UName = GetString(SolverParams, 'Displacement Variable Name', Found)
  IF (.NOT. Found ) UName = 'displacement'
  Found = .FALSE.
  DO i=1,Model % NumberOfSolvers
    SolverPtr => Model % Solvers(i)
    IF (UName == GetVarName(SolverPtr % Variable)) THEN
      Found = .TRUE.
      EXIT
    END IF
  END DO
  IF (.NOT. Found ) THEN
    CALL Fatal('SpringAssembler', 'Solver associated with Displacement Variable Name > '&
        //TRIM(UName)//' < not found!')
  END IF

  DOFs = SolverPtr % Variable % DOFs
  ALLOCATE(K(DOFs,DOFs), F(DOFs), Assemble(DOFs))
  !
  ! Create a point element to serve calling a standard assembly routine:
  !
  Element % Type => GetElementType(101, .FALSE.)
  Element % NDOFs = 1

  !
  ! Loop over all BCs and seek for spring constraints DEFINED NODEWISE
  ! (elementwise conditions are still handled in the solver code):
  !
  DO BC=1,Model % NumberOfBCs
    ValueList => Model % BCs(BC) % Values
    NodesFound = ListCheckPresent(ValueList, 'Target Nodes')
    
    IF (.NOT. NodesFound) CYCLE
 
    NodeIndexes => ListGetIntegerArray(ValueList, 'Target Nodes')
    n = SIZE(NodeIndexes)

    IF (.NOT. ALLOCATED(Work)) THEN 
      ALLOCATE(Work(DOFs,n), STAT=istat)
    ELSE
      IF (SIZE(Work,2) < n) THEN
        DEALLOCATE(Work)
        ALLOCATE(Work(DOFs,n), STAT=istat)
      END IF
    END IF

    Work(1,1:n) = ListGetReal(ValueList, 'Spring', n, NodeIndexes, Found)
    IF (Found) THEN
      CALL Warn('SpringAssembler', 'Define a spring by components (Spring i = ...)')
      CALL Warn('SpringAssembler', 'Skipping a definition (Spring = ...)')
    END IF
    
    Work = 0.0_dp
    DO i=1,DOFs
      Work(i,1:n) = ListGetReal(ValueList, ComponentName('Spring',i), n, &
          NodeIndexes, Assemble(i))
    END DO
    IF (ALL(.NOT. Assemble(:))) CYCLE

    f = 0.0_dp

    DO i=1,n
      Element % NodeIndexes => NodeIndexes(i:i)
      IF (SolverPtr % Variable % Perm(Element % NodeIndexes(1)) <= 0) CYCLE
      K = 0.0_dp
      DO j=1,DOFs
        IF (.NOT. Assemble(j)) CYCLE
        K(j,j) = Work(j,i)
      END DO
      CALL DefaultUpdateEquations(K, f, Element, SolverPtr)
    END DO

  END DO

  IF (ALLOCATED(Work)) DEALLOCATE(Work)
  DEALLOCATE(K, f, Assemble)

!------------------------------------------------------------------------------
END SUBROUTINE SpringAssembler
!------------------------------------------------------------------------------
