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
! *  Author: Eef v. Dongen, Thomas Zwinger
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 25-01-2017, based upon GroundedMask
! *****************************************************************************
!> Solver for creating a nodewise mask on whether the lower surface of an ice shelf is
!> far enough away from GL to be modeled with SSA.
!> NOTE: since floatation criterum check is not relevant anymore, the whole coupling could be based on
!> variable gl distance in stead and this solver would not be necessary?
!> +1=FS,-1=SSA, 0=coupling interface (=last FS node, first SSA node)
SUBROUTINE SSAmask( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  For the bottom surface, creates and updates a mask which may be equal to -1, 0 or 1
!
!  SSAmask = + 1 if to be modeled with FS
!               = - 1 if to be modeled with SSA
!               = 0   if on the coupling interface (last node modeled with FS, first SSA)
!
!  Consequently, a node is floating if SSAmask <= 0 (st similar to GroundedMask)
!  An element is allowed to be modeled with SSA if the minimum distance to GL >= GLDistToler
! !note that sea level might change over time! now set at 0.0!
!
! NOTE: needs distance wrt GL as input
!
!  y is the vertical in 2D ; z is the vertical in 3, works for structured meshes only
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
  USE ElementDescription

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
  TYPE(Variable_t), POINTER :: PointerToVariable, GLVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat,UnFoundFatal=.TRUE.
  LOGICAL :: FirstTime = .True., NewTime
  INTEGER :: i, mn, n, t, Nn, istat, DIM, MSum, ZSum, surfSource
  INTEGER, POINTER :: Permutation(:),  GLPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: z, gldisttoler
  REAL(KIND=dp), ALLOCATABLE :: gldist(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'SSAmask', gldistName

  INTEGER,PARAMETER :: MATERIAL_DEFAULT = 1, MATERIAL_NAMED = 2, VARIABLE = 3

  SAVE AllocationsDone, DIM, SolverName
  SAVE gldisttoler, gldist
  !------------------------------------------------------------------------------


  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  CALL INFO(SolverName, 'Computing SSA mask from gl distance', level=3)

  !--------------------------------------------------------------
  ! Allocate some permanent storage:
  !--------------------------------------------------------------


  SolverParams => GetSolverParams()

  gldistName = GetString(SolverParams, 'GL distance Variable', GotIt)
  IF (.NOT. GotIt) CALL FATAL(SolverName, 'No GL distance Variable given')
  gldisttoler = GetConstReal(SolverParams, 'GL distance Toler', GotIt)
  IF (.NOT.GotIt) THEN
    CALL info(SolverName, 'No GL distance Toler specified, set to default 50 km', level=3)
    gldisttoler = 50000_dp
  END IF

  IF ( (.NOT. AllocationsDone) .OR. Solver % Mesh % Changed ) THEN
    DIM = CoordinateSystemDimension()
    mn = Solver % Mesh % MaxElementNodes
    IF (AllocationsDone) DEALLOCATE(GLdist)
    ALLOCATE(gldist(mn),STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL FATAL( SolverName, 'Memory allocation error GLdist.' )
    END IF
    CALL INFO( SolverName, 'Memory allocation done.',Level=1 )
    AllocationsDone = .TRUE.
  END IF

  !--------------------------------------------------------------
  ! find nodes with min distance  GL Toler wrt GL.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    GLVar => VariableGet(Model % Mesh % Variables, gldistName,UnFoundFatal=UnFoundFatal)
    GLPerm => GLVar % Perm
    gldist(1:n) =  GLVar % values(GLPerm(Element % NodeIndexes))
    CALL GetElementNodes( Nodes )

    DO i = 1, n
      Nn = Permutation(Element % NodeIndexes(i))
      IF (Nn==0) CYCLE
      IF(gldist(i)>gldisttoler) THEN
        VariableValues(Nn) = -1.0_dp !far from GL; SSA
      ELSE
        VariableValues(Nn) = 1.0_dp !too close to GL
      END IF
    END DO
  END DO

  !--------------------------------------------------------------
  ! Coupling interface loop to label points at interface.
  !--------------------------------------------------------------
  ! Loop over each element:
  !  if the sum of the element masks is lower than the element number
  !  of nodes minus the number of zeros (i.e. if the element has at
  !  least one SSA node), then each mask equal to 1 is modified
  !  to 0 (i.e. this node is on the coupling interface).
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
          IF (DIM==2) PRINT *, 'Coupling interface, x', Nodes % x( i )
          IF (DIM==3) PRINT *, 'Coupling interface, (x,y)', Nodes % x( i ), Nodes % y( i )
        END IF
      END DO
    END IF
  END DO

  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, 1 )

  CALL INFO( SolverName , 'Done')

END SUBROUTINE SSAmask
!------------------------------------------------------------------------------

!> Solver for computing the areal (line) weights for coupling the SSA residual
!> to the Stokes solution using External Pressure variable for transfer
SUBROUTINE SSAWeights( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Needs a SSA mask variable constructed by solver SSAmask (see above)
!  given by as string by "SSA Mask Name"
!  Assumes linear elements to be used - we just half the line as weight along
!  the interface line at the bedrock
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
  TYPE(Variable_t), POINTER :: PointerToVariable, MaskVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: FirstTime = .TRUE., GotIt, UnFoundFatal=.TRUE.
    INTEGER :: i, i1, i2, k, n, t, DIM, ElemFamily
  INTEGER, POINTER :: Permutation(:),  MaskPerm(:), EdgeMap(:,:)

  REAL(KIND=dp), POINTER :: VariableValues(:), MaskVal(:)
  REAL(KIND=dp) :: ri(2), s2

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'SSAWeight', MaskName


  SAVE FirstTime, DIM, SolverName,MaskName
  !------------------------------------------------------------------------------

  SolverParams => GetSolverParams()

  IF (FirstTime) THEN
    MaskName = GetString(SolverParams, 'SSA Mask Name', GotIt)
    IF (.NOT.GotIt) CALL FATAL(SolverName,'No "SSA Mask Name" found in Solver section')
    DIM = CoordinateSystemDimension()
    FirstTime = .FALSE.
  END IF
  
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  MaskVar => VariableGet(Model % Mesh % Variables,MaskName,UnFoundFatal=UnFoundFatal)
  IF (.NOT.ASSOCIATED(MaskVar)) THEN
    WRITE(Message,'(A,A)') TRIM(MaskName),' not found'
    CALL FATAL(SolverName,Message)
  END IF
  MaskVal  => MaskVar % Values
  MaskPerm => MaskVar % Perm

  VariableValues = 0.0_dp

  ! Loop all elements 
  DO t = 1, Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()

    CALL GetElementNodes( Nodes )

    IF (ANY(MaskVal(MaskPerm(Element % NodeIndexes)) == 0.0_dp)) THEN
      ElemFamily = GetElementFamily( Element )    
      EdgeMap => GetEdgeMap( ElemFamily )
      ri = 0.0_dp
      DO k=1,Element % Type % NumberOfEdges
        i1 = EdgeMap(k,1)
        i2 = EdgeMap(k,2)
        IF ((MaskVal(MaskPerm(Element % NodeIndexes(i1))) == 0) &
             .AND. (MaskVal(MaskPerm(Element % NodeIndexes(i2))) == 0)) THEN
          ri(1) = Nodes % x(i2) - Nodes % x(i1)
          IF (DIM == 3) &
               ri(2) = Nodes % y(i2) - Nodes % y(i1)
          s2 = SQRT( SUM( ri * ri ) )
          VariableValues(Permutation(Element % NodeIndexes(i1))) = &
               VariableValues(Permutation(Element % NodeIndexes(i1))) + 0.25_dp * s2
          VariableValues(Permutation(Element % NodeIndexes(i2))) = &
               VariableValues(Permutation(Element % NodeIndexes(i2))) + 0.25_dp * s2
        END IF
      END DO
    END IF
  END DO
END SUBROUTINE SSAWeights
