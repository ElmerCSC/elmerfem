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
! *  Authors: Olivier Gagliardini, Gael Durand, Rupert Gladstone, Samuel Cook
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!
! Rupert's notes June 2024 (TODO: integrate these into documentation):
!
! Unifying Samuel's code for identifying isolated ungrounded regions with the main
! grounded mask code.
!
! Aim:
!
! 1. All relevant functionality to be accessed through the grounded solver.
! 2. Default behaviour is backward compatible with non-GMvalid grounded solver:
!     one grounded mask that allows isolated ungrounded regions.
! 3. New option to provide a second grounded mask in which isolated ungrounded
!     regions are removed (unlike Samuel's original, this second mask will
!     identify the grounding line itself, i.e. the values of -1, 0, 1 will have
!     the same meaning as the original grounded mask).
! 
! Additional solver option:
!  'Connected mask name = string xxx'
! This needs to correspond to an existing variable, probably an exported variable.
! Samuel's calving front mask also needs to be specified at the appropriate BC.
!
! Example.
!  Add this to the GroundedSolver:
!   Connected mask name = string ConnMask
!   Exported Variable 1 = -dofs 1 "ConnMask"
!  Add this to the front BC:
!   Calving Front Mask = Logical true
!  

!> Solver for creating a mask on whether the lower side of an ice sheet/shelf is
!>  grounded or not. +1=grounded,-1=detached, 0=grounding line (=last grounded node)
SUBROUTINE GroundedSolver( Model,Solver,dt,TransientSimulation )
  !------------------------------------------------------------------------------
  !******************************************************************************
  !
  !  For the bottom surface, creates and updates a mask which may be equal to -1, 0 or 1
  !  
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
  TYPE(Variable_t), POINTER :: PointerToVariable, bedrockVar, FrontVar, LSvar, ConnMaskVar
  TYPE(Nodes_t), SAVE :: Nodes

  LOGICAL :: AllocationsDone = .FALSE., GotIt, stat,UnFoundFatal=.TRUE.,&
       AllGrounded = .FALSE., useLSvar = .FALSE.,                       &
       CheckConn ! check ocean connectivity (creates separate mask without isolated ungrounded regions)

  INTEGER :: ii, mn, en, t, Nn, istat, DIM, MSum, ZSum, bedrockSource
  INTEGER, POINTER :: Permutation(:), bedrockPerm(:), LSvarPerm(:), ConnMaskPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp) :: z, toler
  REAL(KIND=dp), ALLOCATABLE :: zb(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'GroundedSolver', bedrockName,&
       FrontVarName, LSvarName, ConnMaskName

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

  ConnMaskName = ListGetString(SolverParams, 'Connected mask name',GotIt, UnFoundFatal=.FALSE.)
  IF (GotIt) THEN
     CALL INFO( SolverName, '>Connected mask name< found, connectivity will be checked.',Level=5 )
     CheckConn = .TRUE.
     ConnMaskVar => VariableGet(Model % Mesh % Variables, ConnMaskName,UnFoundFatal=.TRUE.)
     ConnMaskPerm => ConnMaskVar % Perm
  ELSE
     CALL INFO( SolverName, '>Connected mask name< not found, not using.',Level=5 )
     CheckConn = .FALSE.
  END IF
     
  !This to enforce all nodes grounded when doing non-calving hydrology to
  !restart a calving simulation from
  AllGrounded = GetLogical(SolverParams, 'All Grounded', GotIt)
  IF(.NOT. GotIt) AllGrounded = .FALSE.
  IF(.NOT. AllGrounded) THEN

    ! If the user gives a lower surface variable name then use this instead of
    ! node coords for the lower ice surface.
    LSvarName = GetString(SolverParams, 'lower surface variable', GotIt)
    IF (GotIt) THEN
      CALL info(SolverName, 'lower surface variable name found', level=8)
      useLSvar = .TRUE.
    ELSE
      useLSvar = .FALSE.
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
  END IF

  !Any variable defined on the calving front
  FrontVarName = GetString(SolverParams, 'Front Variable', GotIt)
  IF(GotIt) THEN
    FrontVar => VariableGet(Model % Mesh % Variables, FrontVarName,UnFoundFatal=UnFoundFatal)
  ELSE
    CALL INFO( SolverName , 'No front variable defined. Some basal frontal nodes may be left with GroundedMask=1')
    FrontVar => NULL()
  END IF
    
  !--------------------------------------------------------------
  ! Grounded/floating loop based on height of base above bedrock.
  !--------------------------------------------------------------
  DO t = 1, Solver % NumberOfActiveElements
     Element => GetActiveElement(t)
     en = GetElementNOFNodes()
     
     IF(.NOT. AllGrounded) THEN

       SELECT CASE(bedrockSource)
       CASE (VARIABLE)
          bedrockVar => VariableGet(Model % Mesh % Variables, bedrockName,UnFoundFatal=UnFoundFatal)
          bedrockPerm => bedrockVar % Perm
          zb(1:en) =  bedrockVar % values(bedrockPerm(Element % NodeIndexes)) + toler
          NULLIFY(bedrockPerm)
          NULLIFY(bedrockVar)
       CASE (MATERIAL_NAMED)
          Material => GetMaterial( Element )
          zb(1:en) = ListGetReal( Material,bedrockName, en , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
       CASE (MATERIAL_DEFAULT)
          Material => GetMaterial( Element )
          zb(1:en) = ListGetReal( Material,'Min Zs Bottom',en , & 
             Element % NodeIndexes, GotIt,UnFoundFatal=UnFoundFatal) + toler
       END SELECT
     END IF
     
     CALL GetElementNodes( Nodes )
     
     DO ii = 1, en
        Nn = Permutation(Element % NodeIndexes(ii))
        IF (Nn==0) CYCLE
        !To enforce grounding
        IF(AllGrounded) THEN
          VariableValues(Nn) = 1.0_dp
          CYCLE
        END IF

        IF (useLSvar) THEN
          LSvar => VariableGet(Model % Mesh % Variables, LSvarName, UnFoundFatal=UnFoundFatal)
          LSvarPerm => LSvar % Perm
          z = LSvar % values( LSvarPerm(Element % NodeIndexes(ii)) )
        ELSE
          IF (DIM == 2) THEN
            z = Nodes % y( ii )
          ELSE IF (DIM == 3) THEN
            z = Nodes % z( ii )
          END IF

        END IF
        
        ! Geometrical condition. Is the node is above the bedrock 
        ! (plus the tolerance)?  Note: zb includes tolerance.
        IF (z > zb(ii)) THEN
          VariableValues(Nn) = -1.0_dp
        ELSE
          VariableValues(Nn) = 1.0_dp
        END IF
     END DO
  END DO

  ! Check connectivity of ungrounded regions to the front (previously GMvalid solver)
  IF (CheckConn) CALL FrontConn( )
  
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
     en = GetElementNOFNodes()
     CALL GetElementNodes( Nodes )
     MSum = 0
     ZSum = 0
     
     DO ii = 1, en
        Nn = Permutation(Element % NodeIndexes(ii))
        IF (Nn==0) CYCLE
        MSum = MSum + VariableValues(Nn)
        IF (ABS(VariableValues(Nn))<AEPS) ZSum = ZSum + 1.0_dp
     END DO

     IF (MSum + ZSum < en) THEN
        DO ii = 1, en
           Nn = Permutation(Element % NodeIndexes(ii))
           IF (Nn==0) CYCLE
           IF (ABS(VariableValues(Nn)-1.0_dp)<AEPS) THEN
              VariableValues(Nn) = 0.0_dp
              IF (DIM==2) PRINT *, 'Grounding Line, x', Nodes % x( ii )
              IF (DIM==3) PRINT *, 'Grounding Line, (x,y)', Nodes % x( ii ), Nodes % y( ii )
           END IF
        END DO
     END IF

     IF (CheckConn) THEN
        MSum = 0
        ZSum = 0
        DO ii = 1, en
           Nn = ConnMaskPerm(Element % NodeIndexes(ii))
           IF (Nn==0) CYCLE
           MSum = MSum + ConnMaskVar % Values(Nn)
           IF (ABS(VariableValues(Nn))<AEPS) ZSum = ZSum + 1.0_dp
        END DO
        IF (MSum + ZSum < en) THEN
           DO ii = 1, en
              Nn = ConnMaskPerm(Element % NodeIndexes(ii))
              IF (Nn==0) CYCLE
              IF (ABS(ConnMaskVar % Values(Nn)-1.0_dp)<AEPS) THEN
                 ConnMaskVar % Values(Nn) = 0.0_dp
              END IF
           END DO
        END IF
     END IF

    !To label all basal frontal nodes not already ungrounded or on GL as on GL -
    !by definition, they're the last grounded node. Also necessary to make
    !plumes work properly.
     IF(ASSOCIATED(FrontVar)) THEN
       DO ii=1, en
         Nn = FrontVar % Perm(Element % NodeIndexes(ii))
         IF(Nn==0) CYCLE
         Nn = Permutation(Element % NodeIndexes(ii))
         IF(Nn==0) CYCLE
         IF (VariableValues(Nn) > 0.0) VariableValues(Nn) = 0.0_dp
       END DO
     END IF
  END DO
  
  IF (CheckConn) THEN
     IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, ConnMaskVar % Values, OPER_MIN )
  END IF
  IF ( ParEnv % PEs>1 ) CALL ParallelSumVector( Solver % Matrix, VariableValues, OPER_MIN )
  
  CALL INFO( SolverName , 'Done')
  
CONTAINS
  
  ! *****************************************************************************/
  ! * An improved version of the routine to calculate basal melt rates on
  ! * ungrounded ice, producing a validity mask instead (1 = ungrounded area
  ! * connected to the ice front; 0 = isolated patch).
  ! ******************************************************************************
  ! *
  ! *  Authors: Samuel Cook
  ! *  Email:   samuel.cook@univ-grenoble-alpes.fr
  ! *  Web:     http://www.csc.fi/elmer
  ! *  Address: CSC - IT Center for Science Ltd.
  ! *           Keilaranta 14
  ! *           02101 Espoo, Finland
  ! *
  ! *  Original Date: 08.2019
  ! *
  ! ****************************************************************************/
  SUBROUTINE FrontConn ()
    USE Types
    USE CoordinateSystems
    USE DefUtils
    USE ElementDescription
    USE CalvingGeometry
    
    IMPLICIT NONE
    
    !-----------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Matrix_t), POINTER :: Matrix
    TYPE(Variable_t), POINTER :: Var, GroundedVar
    TYPE(CrevasseGroup3D_t), POINTER :: FloatGroups, CurrentGroup, DelGroup
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    TYPE(GaussIntegrationPoints_t) :: IntegStuff
    
    REAL(KIND=dp) :: GMCheck, SMeltRate, WMeltRate, SStart, SStop, &
         TotalArea, TotalBMelt, ElemBMelt, s, t, season,&
         SqrtElementMetric,U,V,W,Basis(Model % MaxElementNodes)
    INTEGER :: NoNodes, j, FaceNodeCount, GroupNodeCount, county, &
         Active, ierr, kk, FoundNew, AllFoundNew
    INTEGER, PARAMETER :: FileUnit = 75, MaxFloatGroups = 1000, MaxNeighbours = 20
    INTEGER, POINTER :: Perm(:), InvPerm(:), FrontPerm(:)=>NULL(), Neighbours(:,:), &
         NeighbourHolder(:), NoNeighbours(:), NodeIndexes(:)
    INTEGER, ALLOCATABLE :: AllGroupNodes(:), PartNodeCount(:), AllPartGroupNodes(:), &
         disps(:)
    LOGICAL :: Found, OutputStats, Visited=.FALSE., Debug, stat, Summer
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, GMaskVarName, FrontMaskName, OutfileName, mode
    
    Debug = .FALSE.
    
    SolverName = "GM Front connectivity"
    Mesh => Solver % Mesh
    
    !Identify nodes on the front
    FrontMaskName = "Calving Front Mask"
    CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
         .FALSE., FrontPerm, FaceNodeCount)
    
    !Need the matrix for finding neighbours
    Matrix => Solver % Matrix
    
    IF(.NOT. ASSOCIATED(ConnMaskVar)) CALL Fatal(SolverName, "Front connectivity needs a variable!")
    ConnMaskVar % Values = 1.0_dp
    
    NoNodes = COUNT(ConnMaskPerm > 0)

    !    Model, Solver, dt, TransientSimulation, ConnMaskVar
!    Var => Solver % Variable
!          VariableValues(Nn) = 1.0_dp
!  PointerToVariable => Solver % Variable
!  Permutation  => PointerToVariable % Perm
!  VariableValues => PointerToVariable % Values

!    GMaskVarName = ListGetString(Params, "GroundedMask Variable", Found)
!    IF(.NOT. Found) GMaskVarName = "GroundedMask"
!    GroundedVar => VariableGet(Mesh % Variables, GMaskVarName, .TRUE., UnfoundFatal=.TRUE.)

    GroundedVar => Solver % Variable
    
    GMCheck = -1.0_dp
    
    !Set up inverse perm for FindNodeNeighbours
    InvPerm => CreateInvPerm(Matrix % Perm) !Create inverse perm for neighbour search
    ALLOCATE(Neighbours(Mesh % NumberOfNodes, MaxNeighbours), NoNeighbours(Mesh % NumberOfNodes))
    Neighbours = 0
    
    !Find neighbours for each node on the bed
    DO ii = 1, Mesh % NumberOfNodes
       IF(ConnMaskPerm(ii) <= 0) CYCLE
       
       NeighbourHolder => FindNodeNeighbours(ii, Matrix, &
            Matrix % Perm, 1, InvPerm)
       
       Neighbours(ii,1:SIZE(NeighbourHolder)) = NeighbourHolder
       NoNeighbours(ii) = SIZE(NeighbourHolder)
       DEALLOCATE(NeighbourHolder)
    END DO
    
    !Reuse some old calving code
    !Find groups of connected floating nodes on the base
    FloatGroups => NULL()
    CALL FindCrevasseGroups(Mesh, GroundedVar, Neighbours, &
         -0.5_dp, FloatGroups)
    
    !Check groups are valid (connected to front)
    CurrentGroup => FloatGroups
    DO WHILE(ASSOCIATED(CurrentGroup))
       CurrentGroup % FrontConnected = .FALSE.
       DO ii=1, CurrentGroup % NumberOfNodes
          
          IF(FrontPerm(CurrentGroup % NodeNumbers(ii)) > 0) THEN
             CurrentGroup % FrontConnected = .TRUE.
             EXIT
          END IF
       END DO
       CurrentGroup => CurrentGroup % Next
    END DO
    
    DO kk=1,MaxFloatGroups
       FoundNew = 0
       !Count and gather nodes from all valid groups
       GroupNodeCount = 0
       county = 0
       DO ii=1,2
          IF(ii==2) ALLOCATE(AllGroupNodes(GroupNodeCount))
          CurrentGroup => FloatGroups
          DO WHILE(ASSOCIATED(CurrentGroup))
             IF(CurrentGroup % FrontConnected) THEN
                IF(ii==1) THEN
                   GroupNodeCount = GroupNodeCount + CurrentGroup % NumberOfNodes
                ELSE
                   DO j=1, CurrentGroup % NumberOfNodes
                      county = county + 1
                      AllGroupNodes(county) = Mesh % ParallelInfo % GlobalDOFs(CurrentGroup % NodeNumbers(j))
                   END DO
                END IF
             END IF
             CurrentGroup => CurrentGroup % Next
          END DO
       END DO
       
       !Distribute info to/from all partitions about groups connected to front
       ALLOCATE(PartNodeCount(ParEnv % PEs))
       
       CALL MPI_ALLGATHER(GroupNodeCount, 1, MPI_INTEGER, PartNodeCount, 1, &
            MPI_INTEGER, MPI_COMM_WORLD, ierr)
       
       ALLOCATE(AllPartGroupNodes(SUM(PartNodeCount)), disps(ParEnv % PEs))
       disps(1) = 0
       DO ii=2,ParEnv % PEs
          disps(ii) = disps(ii-1) + PartNodeCount(ii-1)
       END DO
       
       CALL MPI_ALLGATHERV(AllGroupNodes, GroupNodeCount, MPI_INTEGER, &
            AllPartGroupNodes, PartNodeCount, disps, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       
       !Cycle unconnected groups, looking for partition boundary in connected groups
       CurrentGroup => FloatGroups
       DO WHILE(ASSOCIATED(CurrentGroup))
          IF(.NOT. CurrentGroup % FrontConnected) THEN
             DO ii=1,CurrentGroup % NumberOfNodes
                
                IF(ANY(Mesh % ParallelInfo % GlobalDOFs(CurrentGroup % NodeNumbers(ii)) == &
                     AllPartGroupNodes)) THEN
                   CurrentGroup % FrontConnected = .TRUE.
                   FoundNew = 1
                END IF
                
             END DO
          END IF
          CurrentGroup => CurrentGroup % Next
       END DO
       CALL MPI_ALLREDUCE(FoundNew, AllFoundNew, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)
       IF(AllFoundNew == 1) THEN
          DEALLOCATE(AllGroupNodes, PartNodeCount, AllPartGroupNodes, disps)
       ELSE
          EXIT
       END IF
       IF (kk.EQ.MaxFloatGroups) CALL FATAL( SolverName, 'Hard coded loop limit reached; needs recoding!' )
    END DO !k
    
    !Cycle all connected groups, setting melt rate
    CurrentGroup => FloatGroups
    DO WHILE(ASSOCIATED(CurrentGroup))
       IF(CurrentGroup % FrontConnected) THEN
          DO ii=1,CurrentGroup % NumberOfNodes
             ConnMaskVar % Values(ConnMaskVar % Perm(CurrentGroup % NodeNumbers(ii))) = GMCheck
          END DO
       END IF
       CurrentGroup => CurrentGroup % Next
    END DO
    
    !Deallocate floatgroups linked list
    CurrentGroup => FloatGroups
    DO WHILE(ASSOCIATED(CurrentGroup))
       DelGroup => CurrentGroup
       CurrentGroup => CurrentGroup % Next
       
       IF(ASSOCIATED(DelGroup % NodeNumbers)) DEALLOCATE(DelGroup % NodeNumbers)
       IF(ASSOCIATED(DelGroup % FrontNodes)) DEALLOCATE(DelGroup % FrontNodes)
       IF(ASSOCIATED(DelGroup % BoundaryNodes)) DEALLOCATE(DelGroup % BoundaryNodes)
       DEALLOCATE(DelGroup)
    END DO
    
    DEALLOCATE(Neighbours, NoNeighbours, FrontPerm, InvPerm)
  END SUBROUTINE FrontConn
  
END SUBROUTINE GroundedSolver

