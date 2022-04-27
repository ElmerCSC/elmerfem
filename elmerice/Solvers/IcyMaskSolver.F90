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
! *  Authors: Olivier Gagliardini
! *  Email: olivier.gagliardini@univ-grenoble-alpes.fr   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: Dec 2020 
! * 
! *****************************************************************************
!> Solver for creating a mask on whether there is ice (H> Hin) or not (H<Hmin)
!>  +1= Icy (H>Hmin),-1= Ice Free (H<Hmin), 0=contour of the glacier
SUBROUTINE IcyMaskSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  For the Upper surface, creates and updates a mask which may be equal to -1, 0 or 1
!
!  Icy Mask = + 1 if Icy (glacier)
!           = - 1 if Ice Free 
!           = 0   on the contouur of the glacier (first node with H=Hmin)
!           < -1 for Icy Isolated nodes (useful to move them back to initial elevation) 
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
  TYPE(Variable_t), POINTER :: PointerToVariable, ThickVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat, UnFoundFatal=.TRUE., &
             Isolated = .FALSE., IsolatedEdge = .FALSE., FirstTime = .TRUE., &
             FoundIsolated  

  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum
  INTEGER, POINTER :: Permutation(:), ThickPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), ThickValues(:)
  REAL(KIND=dp) :: toler, Hmin
  REAL(KIND=dp), ALLOCATABLE :: Hice(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'IcyMaskSolver'

  INTEGER, PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3
       
  SAVE AllocationsDone, DIM, SolverName, Hice, toler, Hmin, Isolated, IsolatedEdge, &
       FirstTime               
  !------------------------------------------------------------------------------

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing a mask for the glacier extend', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------
  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     mn = Solver % Mesh % MaxElementNodes
     IF (AllocationsDone) DEALLOCATE(Hice)     
     ALLOCATE(Hice(mn), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error.' )
     END IF
     CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
     AllocationsDone = .TRUE.
  END IF
  
  IF (FirstTime) THEN 
    SolverParams => GetSolverParams()
    toler = GetConstReal(SolverParams, 'Toler', GotIt)
    IF (.NOT.GotIt) THEN
      CALL FATAL(SolverName, 'No tolerance given for the Icy Mask.')
    END IF
    Hmin = GetConstReal(SolverParams, 'Ice Free Thickness', GotIt)
    IF (.NOT.GotIt) THEN
      CALL FATAL(SolverName, 'No Ice Free Thickness given for the Icy Mask.')
    END IF
    Isolated = GetLogical(SolverParams, 'Remove Isolated Points', GotIt)
    IF (Isolated) THEN
      IsolatedEdge = GetLogical(SolverParams, 'Remove Isolated Edges', GotIt)
    END IF
    FirstTime = .FALSE. 
  END IF
  

  ! Read the Thickness variable
   ThickVar => VariableGet(Model % Mesh % Variables, 'Thickness', UnFoundFatal=UnFoundFatal)
   ThickValues => ThickVar % Values
   ThickPerm => ThickVar % Perm
  
  ! Initialise to small value to avoid problem in parallel when using halo elements
  VariableValues = -100.0

  !--------------------------------------------------------------
  ! Icy/Ice Free Mask is based on variable Thickness
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    CALL GetElementNodes( Nodes )

    Hice(1:n) = ThickValues(ThickPerm(Element % NodeIndexes)) - toler 

    DO i = 1, n
      Nn = Permutation(Element % NodeIndexes(i))
      IF (Nn==0) CYCLE

      IF (Hice(i) > Hmin) THEN
        VariableValues(Nn) = 2.0_dp
      ELSE
        VariableValues(Nn) = -1.0_dp
      END IF
    END DO
  END DO
     
  IF (Isolated) THEN 
     ! Look for Isolated Nodes (one node surounded only by Ice Free nodes)
     ! These nodes are tagued with Mask value < -1
    DO t = 1, Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      n = GetElementNOFNodes()
      CALL GetElementNodes( Nodes )
      Hice(1:n) = ThickValues(ThickPerm(Element % NodeIndexes)) - toler 
      
      ! Look only at element with one node with H>Hmin or 2 nodes (IsolatedEdge)
      IF (IsolatedEdge) THEN 
        FoundIsolated = (COUNT(Hice(1:n)>Hmin)==1).OR.(COUNT(Hice(1:n)>Hmin)==2)
      ELSE
        FoundIsolated = COUNT(Hice(1:n)>Hmin)==1
      END IF
      IF (FoundIsolated) THEN
         DO i = 1, n
           Nn = Permutation(Element % NodeIndexes(i))
           IF (Nn==0) CYCLE

           IF (Hice(i) > Hmin) THEN
             IF (ABS(VariableValues(Nn)-1.0)<AEPS) THEN !=1 from previous iter
               VariableValues(Nn) = 1.0_dp ! we let it one
             ELSE IF (ABS(VariableValues(Nn)-2.0)<AEPS) THEN ! first time we see that node  
               VariableValues(Nn) = -1.0_dp
             ELSE
               VariableValues(Nn) = VariableValues(Nn)-1.0_dp
             END IF
           ELSE
             VariableValues(Nn) = -1.0_dp
           END IF
         END DO
      ELSE 
        DO i = 1, n
          Nn = Permutation(Element % NodeIndexes(i))
          IF (Nn==0) CYCLE

          IF (Hice(i) > Hmin) THEN
            VariableValues(Nn) = 1.0_dp
          ELSE
            VariableValues(Nn) = -1.0_dp
          END IF
        END DO
      END IF
    END DO
  END IF 
  
  !--------------------------------------------------------------
  ! Loop to label contour points of the glacier 
  !--------------------------------------------------------------
  ! Loop over each element:
  DO t = 1, Solver % NumberOfActiveElements     
     Element => GetActiveElement(t)
     n = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     MSum = 0
     ZSum = 0
     
     DO i = 1, n
        Nn = Permutation(Element % NodeIndexes(i))
        IF (Nn==0) CYCLE
        IF (VariableValues(Nn)>1.5) VariableValues(Nn)=1.0_dp
        MSum = MSum + VariableValues(Nn)
        IF (ABS(VariableValues(Nn)) < AEPS) ZSum = ZSum - 1.0_dp
     END DO
     
     IF (MSum + ZSum > -n) THEN
        DO i = 1, n
           Nn = Permutation(Element % NodeIndexes(i))
           IF (Nn==0) CYCLE
           IF (ABS(VariableValues(Nn)+1.0_dp)<AEPS) THEN
              VariableValues(Nn) = 0.0_dp
             ! IF (DIM==3) PRINT *, 'Glacier contour line, (x,y)', Nodes % x( i ), Nodes % y( i )
           END IF
        END DO
     END IF

  END DO
  
  ! Take the maximal value for nodes at the boundaty of two partitions 
  ! to account for the case of some margin nodes at the interface between two partitions
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 2 )
 
  CALL INFO( SolverName , 'Done')
 
END SUBROUTINE IcyMaskSolver 
!------------------------------------------------------------------------------
