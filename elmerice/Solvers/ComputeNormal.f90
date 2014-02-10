!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
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
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Juha Ruokolainen
! *  Email:   Juha.Ruokolainen [at] csc.fi 
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - Scientific Computing Ltd.  
! *               Keilaranta 14                    
! *               02101 Espoo, Finland             
! *                                                 
! *       Original Date: 14 May 2007                
! * 
! *****************************************************************************
!>  Module containing solver for computation of surface normals
SUBROUTINE ComputeNormalSolver( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!Modification have been done to deal with problems where we don't want to compute
!the normal on all the boundaries (problem with the corners)
!
!Keywords are: ComputeAll = True/False computing / or not the normal on every boundary
!              ComputeNormal = True to compute the normal on a given boundary
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: NormalSolution
  TYPE(ValueList_t), POINTER :: BC, SolverParams

  INTEGER :: i, j, k, n, t, DIM
  REAL(KIND=dp) :: u, v, w, s 
  REAL(KIND=dp), POINTER :: Nvector(:)
  INTEGER, POINTER :: Permutation(:)

  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp) :: Bu, Bv, Normal(3), NormalCond(4)

  LOGICAL :: CompAll = .TRUE., CompBC = .TRUE., Found

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'ComputeNormalSolver'

!---------------------------------------
! Setup pointer to the current solution:
!---------------------------------------
  NormalSolution => VariableGet( Solver % Mesh % Variables, 'Normal Vector' )
  IF ( ASSOCIATED( NormalSolution ) ) THEN 
     Nvector => NormalSolution % Values
     Permutation => NormalSolution % Perm
  ELSE
     PRINT *,'FATAL: Unable to set pointer to the current solution'
     STOP
  END IF

  CALL INFO(SolverName, 'Computing Normal Vector for Nodes', level=3)

  SolverParams => GetSolverParams()
  CompAll = GetLogical (SolverParams, 'ComputeAll', Found)

  IF (.NOT.Found) THEN
    WRITE(Message,'(A)') ('ComputeAll not found, Normal is computed for all the boundaries')
    CALL INFO(SolverName, Message, level = 20)
    CompAll = .TRUE.
    CompBC = .TRUE.
  ELSE
    IF (.NOT.CompAll) CompBC = .FALSE.
    IF (CompAll) CompBC = .TRUE.
  END IF

  DIM = CoordinateSystemDimension()
 
  DO t = 1, Solver % Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    IF (ParEnv % myPe .NE. Element % partIndex) CYCLE
    n = GetElementNOFNodes(Element)
    BC => GetBC(Element)
    IF (n == 1) CYCLE

    CALL GetElementNodes( Nodes, Element )

    IF (.NOT.CompAll) THEN
      CompBC = GetLogical ( BC,'ComputeNormal',Found)
      IF (.NOT.Found) THEN
          NormalCond = 0.0
          NormalCond(1:n) = GetReal( BC, 'ComputeNormal Condition', Found )

	! If at least one value in NormalCond > 0, then CompBC=.true.
        IF (COUNT(NormalCond > 0.0) > 0) CompBC = .TRUE.
      END IF
    END IF

    IF (CompBC) THEN

      DO i = 1,n
        IF (NormalCond(i) .LT. 0) CYCLE
        j = Element % NodeIndexes( i )
        k = Permutation(j)
           
        Bu = Element % Type % NodeU(i)
        IF ( Element % Type % Dimension > 1 ) THEN
          Bv = Element % Type % NodeV(i)
        ELSE
          Bv = 0.0D0
        END IF
        Normal = NormalVector(Element, Nodes, Bu, Bv, .TRUE.)
        Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k) +& 
        Normal(1:DIM)
      END DO
    END IF
  END DO
  
  DO i=1,Model % NumberOfNodes
    k = Permutation(i)
    IF ( k > 0 ) THEN
      s = SQRT( SUM( Nvector(DIM*(k-1)+1:DIM*k)**2 ) )

      IF ( s /= 0.0D0 ) THEN
	Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k)/s 
      END IF
    END IF

  END DO

  CALL INFO(SolverName, 'End', level=3)

!------------------------------------------------------------------------------
END SUBROUTINE ComputeNormalSolver
!------------------------------------------------------------------------------
