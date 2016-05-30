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
! *  Authors: Olivier Gagliardini, Gael Durand
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> Solver for creating a mask on whether the lower side of an ice sheet/shelf is
!>  grounded or not. +1=grounded,-1=detached, 0=grounding line (=last grounded node)
SUBROUTINE GroundedSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  For the bottom surface, creates and updates a mask which may be equal to -1, 0 or 1

!  GroundedMask = + 1 if grounded
!               = - 1 if floating
!               = 0   if on the grounding line (also grounded)
!
!  Consequently, a node is grounded if GroundedMask >= 0
!
!  y is the vertical in 2D ; z is the vertical in 3D
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
!******************************************************************************
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
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat,UnFoundFatal=.TRUE.

  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum, bedrockSource
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: z, toler
  REAL(KIND=dp), ALLOCATABLE :: zb(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE AllocationsDone, DIM, SolverName, zb, toler
  !------------------------------------------------------------------------------

!  NULLIFY(bedrockPerm,bedrockVar)

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing grounded mask from geometry', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     mn = Solver % Mesh % MaxElementNodes
     IF (AllocationsDone) DEALLOCATE(zb)     
     ALLOCATE(zb(mn), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error.' )
     END IF
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
     AllocationsDone = .TRUE.
  END IF
  
  SolverParams => GetSolverParams()
  toler = GetConstReal(SolverParams, 'Toler', GotIt)
  IF (.NOT.GotIt) THEN
     CALL FATAL(SolverName, 'No tolerance given for the Grounded Mask.')
  END IF

  bedrockName = GetString(SolverParams, 'Bedrock Variable', GotIt)
  IF (GotIt) THEN
     bedrockSource = VARIABLE
     CALL info(SolverName, 'Bedrock Variable name found', level=8)
  ELSE
     bedrockName = GetString(SolverParams, 'Bedrock Material', GotIt)
     IF (GotIt) THEN
        bedrockSource = MATERIAL_NAMED
        CALL info(SolverName, 'Bedrock Material name found', level=8)
     ELSE
        bedrockSource = MATERIAL_DEFAULT     
        CALL info(SolverName, 'No Bedrock Variable or Material; searching for material \"Min Zs Bottom\".', level=8)
     END IF
  END IF


    
  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     
     SELECT CASE(bedrockSource)
     CASE (VARIABLE)
        bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName,UnFoundFatal=UnFoundFatal)
        bedrockPerm => bedrockVar % Perm
        zb(1:n) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
        NULLIFY(bedrockPerm)
        NULLIFY(bedrockVar)
     CASE (MATERIAL_NAMED)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,bedrockName, n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     CASE (MATERIAL_DEFAULT)
        Material => GetMaterial( Element )
        zb(1:n) = ListGetReal( Material,'Min Zs Bottom',n , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
     END SELECT
     
     CALL GetElementNodes( Nodes )
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        IF (DIM == 2) THEN
           z = Nodes % y( i )
        ELSE IF (DIM == 3) THEN
           z = Nodes % z( i )
        END IF
        
        ! Geometrical condition. Is the node is above the bedrock 
        ! (plus the tolerance)?  Note: zb includes tolerance.
        IF (z > zb(i)) THEN
           VariableValues(Nn) = -1.0_dp
        ELSE
           VariableValues(Nn) = 1.0_dp
        END IF
     END DO
  END DO
  
  !--------------------------------------------------------------
  ! Grounding line loop to label grounded points at grounding Line.
  !--------------------------------------------------------------
  ! Loop over each element:
  !  if the sum of the element masks is lower than the element number 
  !  of nodes minus the number of zeros (i.e. if the element has at 
  !  least one floating node), then each mask equal to 1 is modified 
  !  to 0 (i.e. this node is on the grounding line).  
  DO t = 1, Solver % NumberOfActiveElements     
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     MSum = 0
     ZSum = 0
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        MSum = MSum + VariableValues(Nn)
        IF (VariableValues(Nn) == 0.0_dp) ZSum = ZSum + 1.0_dp
     END DO
     
     IF (MSum + ZSum < n) THEN
        DO i = 1, n
           Nn = Permutation(Element % NodeIndexes(i))
           IF (Nn==0) CYCLE
           IF (VariableValues(Nn) == 1.0_dp) THEN
              VariableValues(Nn) = 0.0_dp
              IF (DIM==2) PRINT *, 'Grounding Line, x', Nodes % x( i )
              IF (DIM==3) PRINT *, 'Grounding Line, (x,y)', Nodes % x( i ), Nodes % y( i )
           END IF
        END DO
     END IF
  END DO
  
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )
 
  CALL INFO( SolverName , 'Done')
 
END SUBROUTINE GroundedSolver 
!------------------------------------------------------------------------------
SUBROUTINE GroundedSolverInit( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
! for Grounded Mask initialisation purpose
! same method than above
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolverInit'

  CALL FATAL( SolverName, 'This solver is deprecated due to code redundancy, &
       please GroundedSolver instead' )

!------------------------------------------------------------------------------
END SUBROUTINE GroundedSolverInit 
!------------------------------------------------------------------------------



