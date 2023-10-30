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
! *  Subroutine utilising ParMMG to perform parallel remeshing
! *
! ******************************************************************************
! *
! *  Authors:
! *  Email:   
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 6/4/22
! *
! *****************************************************************************/

! Subroutine to remesh entire mesh in parallel
! Uses MeshRemeshing ParMMG remeshing routines
!------------------------------------------------------------------------------
SUBROUTINE ParallelRemesh( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE MeshPartition
  USE MeshRemeshing
  USE MainUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!-----------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh, OutMesh, FinalMesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: TimeVar
  REAL(kind=dp), POINTER :: WorkReal(:)
  INTEGER, POINTER :: WorkPerm(:)
  INTEGER, ALLOCATABLE :: EdgePairs(:,:)
  LOGICAL :: Parallel, Rebalance, Found, ManAssignEdges, Angle
  LOGICAL, ALLOCATABLE :: Boundaries(:)
  INTEGER :: i,j,k,n, PairCount, ierr, parts, Time
  REAL(kind=dp) :: TimeReal, DetAngle
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = "ParallelRemesh", MeshDir, MeshName, tmp,&
      SaveMeshName
!------------------------------------------------

#ifdef HAVE_PARMMG

  IF(ParEnv % PEs > 1) THEN
      Parallel = .TRUE.
  ELSE
      CALL FATAL(Solvername, 'Needs to be run in parallel')
  END IF

  Mesh => Model % Mesh
  SolverParams => GetSolverParams()

  Rebalance = ListGetLogical(SolverParams, "Rebalance", Found, DefValue = .TRUE.)
  IF(.NOT. Found) CALL INFO(SolverName, "Option 'Rebalance' not found. DefValue is TRUE")
  ManAssignEdges = ListGetLogical(SolverParams, "Manually Assign Edges", Found, DefValue = .FALSE.)
  IF(.NOT. Found) CALL INFO(SolverName, "Option 'Manually Assign Edges' not found. Default is FALSE")
  Angle = ListGetLogical(SolverParams, "Automatic Angle Detection", Found, DefValue = .FALSE.)
  IF(.NOT. Found) CALL INFO(SolverName, "Automatic Angle Detection turned off")

  TimeVar => VariableGet( Model % Mesh % Variables, 'Timestep' )
  TimeReal = TimeVar % Values(1)
  Time = INT(TimeReal)

  SaveMeshName = Mesh % Name

  IF(ManAssignEdges) THEN
    IF(Angle) CALL WARN(SolverName, "Manually assigning edge elems so cannot automatically detect angles")
    ALLOCATE(Boundaries(Model % NumberOFBCs))
    Boundaries = .TRUE.
    CALL GetRemeshEdgeNodes(Mesh, Boundaries, EdgePairs, PairCount)
    CALL DistributedRemeshParMMG(Model, Mesh,OutMesh, EdgePairs, PairCount,Params=SolverParams)
  ELSE

    IF(Angle) THEN
      DetAngle = ListGetConstReal(SolverParams,"Detection Angle", Found)
      IF(.NOT. Found) CALL FATAL(SolverName, "Detection Angle not given!")
      CALL DistributedRemeshParMMG(Model, Mesh,OutMesh, Params=SolverParams, Angle=DetAngle)
    ELSE
      CALL DistributedRemeshParMMG(Model, Mesh,OutMesh, Params=SolverParams)
    END IF
  END IF

  !global element numbers - node numbers already assigned in recovery of PMMG mesh
  CALL RenumberGElems(OutMesh)

  IF(Rebalance) THEN
    CALL Zoltan_Interface( Model, OutMesh, StartImbalanceTol=1.1_dp, TolChange=0.02_dp, MinElems=10 )

    FinalMesh => RedistributeMesh(Model, OutMesh, .TRUE., .FALSE.)
    CALL ReleaseMesh(OutMesh)
  ELSE
    FinalMesh => OutMesh
  END IF

  MeshDir = ListGetString(SolverParams,"Save Mesh Name", Found ) 
  IF(.NOT. Found) MeshDir = SaveMeshName

  WRITE(MeshName, '(A,i0)') TRIM(MeshDir), time  
  OutMesh => NULL()
  FinalMesh % Name = TRIM(MeshName)
  FinalMesh % OutputActive = .TRUE.
  FinalMesh % Changed = .TRUE.
  

  !MakeDirectory seems unable to create multi/level/directories
  !so create the top level first, then the lower
  IF(ParEnv % MyPe==0) CALL MakeDirectory(TRIM(MeshName) // CHAR(0))

  parts = ParEnv % PEs

  tmp = TRIM(MeshName)//"/partitioning."
  MeshDir = TRIM(tmp)

  WRITE (tmp, '(A,I0)') TRIM(MeshDir),parts
  MeshDir = TRIM(tmp)

  IF(ParEnv % MyPe==0) THEN
    PRINT *, 'Save Mesh, creating directory: '//MeshDir
    CALL MakeDirectory(TRIM(MeshDir) // CHAR(0))
  END IF
  CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

  CALL WriteMeshToDisk2(Model, FinalMesh, MeshDir, ParEnv % MyPE)
    
  CALL SwapMesh(Model, Mesh, FinalMesh % Name)
  Model % Mesh % Name = TRIM(SaveMeshName)

#else
  CALL FATAL(SolverName, 'ParMMG needs to be installed to use this solver!')
#endif

  CONTAINS

  ! Cycle through all 303 elements of GatheredMesh, creating lists of those
  ! on the given boundaries
  ! Cycle these lists, identifying elements on different boundaries, which
  ! share nodes (therefore share a theoretical 202 element), construct
  ! list of these 202 elements
  ! These then fed into ParMMG (ParallelRemesh)
  ! copied and altered from GetCalvingRemeshNodes in CalvingGeometry.F90
  SUBROUTINE GetRemeshEdgeNodes(Mesh, Boundaries, Shared, TotalCount)

    USE Types
    USE SParIterComm
    USE MainUtils

    IMPLICIT NONE

    TYPE(Mesh_t),POINTER :: Mesh
    TYPE(Element_t),POINTER :: Element
    LOGICAL :: Boundaries(:)
    !---------------
    INTEGER :: i,j,k, BoundaryNumber, NumNodes, Match, BoundaryID, TotalCount, &
          FirstBdryID, SecondBdryID, CountSoFar, ReqBCs, TotalBCs
    INTEGER, ALLOCATABLE :: ElementNodes(:), Counts(:), BdryNodes(:,:,:), &
          CountPairs(:,:),SharedPairs(:,:,:,:),Shared(:, :)
    LOGICAL :: Debug, Counted, FirstMatch, SecondMatch, ThirdMatch
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    SolverName = 'GetRemeshEdgeNodes'

    ReqBCs = COUNT(Boundaries)
    TotalBCs = SIZE(Boundaries)
    ALLOCATE(Counts(TotalBCs))
    counts = 0
    
    ALLOCATE(BdryNodes(TotalBCs,3,100))

    DO i=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(i)
      
      ElementNodes = Element % NodeIndexes
      BoundaryNumber = Element % BoundaryInfo % constraint
        
      NumNodes = Element % TYPE % NumberOfNodes
      IF (NumNodes /= 3) CALL FATAL(Solvername, "BoundaryElements must be 303s")

      DO BoundaryID=1,TotalBCs
        IF(.NOT. Boundaries(BoundaryID)) CYCLE
        IF (BoundaryNumber == BoundaryID) THEN
            Counts(BoundaryID) = Counts(BoundaryID) + 1
          IF (Counts(BoundaryID) > SIZE(BdryNodes(BoundaryID,1,:))) THEN
            IF(Debug) PRINT *, BoundaryID, 'BdryNodes, doubling array size'
            CALL DoubleBdryArraySize(BdryNodes)
          END IF
          BdryNodes(BoundaryID,:,Counts(BoundaryID)) = ElementNodes
        END IF    
      END DO
    END DO

    !set counts for calving and other boundary shared nodes
    ALLOCATE(CountPairs(TotalBCs-1,TotalBCs-1))
    CountPairs(:,:) = 0
      
    !set allocatables
    ALLOCATE(SharedPairs(TotalBCs-1,TotalBCs-1,2,100))

    ! loop for 1-2, 1-3 ... 1-6, 2,3 ... 5,6
    DO FirstBdryID=1,TotalBCs-1
      IF(.NOT. Boundaries(FirstBdryID)) CYCLE
      IF (Counts(FirstBdryID) /= 0) THEN
        DO i=1, Counts(FirstBdryID)
          DO SecondBdryID=FirstBdryID+1,TotalBCs
            IF(.NOT. Boundaries(SecondBdryID)) CYCLE
            IF (Counts(SecondBdryID) /= 0) THEN
              DO j=1, Counts(SecondBdryID)
                Match = 0
                FirstMatch=.FALSE.
                SecondMatch=.FALSE.
                ThirdMatch=.FALSE.
                DO k=1,3
                  IF (BdryNodes(FirstBdryID,1,i) == BdryNodes(SecondBdryID,k,j)) THEN
                    FirstMatch=.TRUE.
                    Match = Match + 1
                  END IF
                  IF (BdryNodes(FirstBdryID,2,i) == BdryNodes(SecondBdryID,k,j)) THEN
                    SecondMatch=.TRUE.
                    Match = Match + 1
                  END IF
                  IF (BdryNodes(FirstBdryID,3,i) == BdryNodes(SecondBdryID,k,j)) THEN
                    ThirdMatch=.TRUE.
                    Match = Match + 1
                  END IF
                END DO
                IF (Match == 2) THEN
                  CountPairs(FirstBdryID,SecondBdryID-FirstBdryID) = CountPairs(FirstBdryId,SecondBdryID-FirstBdryID) + 1
                  IF (CountPairs(FirstBdryID,SecondBdryID-FirstBdryID) > &
                        SIZE(SharedPairs(FirstBdryID,SecondBdryID-FirstBdryID,1,:))) THEN
                    IF(Debug) PRINT *, 'SharedPairs boundaryIDs-,',FirstBdryID,SecondBdryID,'doubling size of node array.'
                    CALL DoubleEdgeArraySize(SharedPairs)
                  END IF
                  IF (FirstMatch .AND. SecondMatch) THEN
                    SharedPairs(FirstBdryID,SecondBdryID-FirstBdryID,:,CountPairs(FirstBdryID,SecondBdryID-FirstBdryID)) &
                        = BdryNodes(FirstBdryID,1:2,i)
                  ELSE IF (SecondMatch .AND. ThirdMatch) THEN
                    SharedPairs(FirstBdryID,SecondBdryID-FirstBdryID,:,CountPairs(FirstBdryID,SecondBdryID-FirstBdryID)) &
                        = BdryNodes(FirstBdryID,2:3,i)
                  ELSE IF (FirstMatch .AND. ThirdMatch) THEN
                    SharedPairs(FirstBdryID,SecondBdryID-FirstBdryID,1,CountPairs(FirstBdryID,SecondBdryID-FirstBdryID)) &
                        = BdryNodes(FirstBdryID,1,i)
                    SharedPairs(FirstBdryID,SecondBdryID-FirstBdryID,2,CountPairs(FirstBdryID,SecondBdryID-FirstBdryID)) &
                        = BdryNodes(FirstBdryID,3,i)
                  END IF
                ELSE IF (Match == 3) THEN
                  PRINT*, 'BoundaryElement: Duplicated', FirstBdryID,BdryNodes(FirstBdryID,:,i)
                  PRINT*, 'BoundaryElement: Duplicated', SecondBdryID,BdryNodes(SecondBdryID,:,j), j
                  CALL FATAL(Solvername, "BoundaryElement: Duplicated")
                END IF
              END DO
            END IF
          END DO
        END DO
      END IF
    END DO
      
    TotalCount=0
    DO i=1,TotalBCs-1
      IF(.NOT. Boundaries(i)) CYCLE
      DO j=1,TotalBCs-1
        IF(.NOT. Boundaries(j)) CYCLE
        TotalCount=TotalCount+CountPairs(i,j)
      END DO
    END DO

    ALLOCATE(Shared(2, TotalCount))

    CountSoFar=0
    DO i=1,TotalBCs-1
      DO j=1,TotalBCs-1
        Shared(:,1+CountSoFar:CountSoFar+CountPairs(i,j)) = SharedPairs(i,j,:,1:CountPairs(i,j))
        CountSoFar = CountSoFar + CountPairs(i,j)
      END DO
    END DO

  END SUBROUTINE GetRemeshEdgeNodes

  SUBROUTINE DoubleEdgeArraySize(Vec, fill)
    !only doubles size in one dimension
    INTEGER, ALLOCATABLE :: Vec(:,:,:,:)
    INTEGER, OPTIONAL :: fill
    !----------------------------------------
    INTEGER, ALLOCATABLE :: WorkVec(:,:,:,:), D(:)
    
    ALLOCATE(D(3))
    d = SHAPE(Vec)

    ALLOCATE(WorkVec(d(1),d(2),d(3),d(4)))

    WorkVec = Vec
    
    DEALLOCATE(Vec)
    ALLOCATE(Vec(d(1),d(2),d(3),2*d(4)))

    IF(PRESENT(fill)) THEN
      Vec = fill
    ELSE
      Vec = 0
    END IF

    Vec(:,:,:,1:d(4)) = WorkVec

  END SUBROUTINE DoubleEdgeArraySize

  SUBROUTINE DoubleBdryArraySize(Vec, fill)
    !only doubles size in one dimension
    INTEGER, ALLOCATABLE :: Vec(:,:,:)
    INTEGER, OPTIONAL :: fill
    !----------------------------------------
    INTEGER, ALLOCATABLE :: WorkVec(:,:,:), D(:)
    
    ALLOCATE(D(3))
    d = SHAPE(Vec)

    ALLOCATE(WorkVec(d(1), d(2),d(3)))

    WorkVec = Vec
    
    DEALLOCATE(Vec)
    ALLOCATE(Vec(d(1),d(2),2*d(3)))

    IF(PRESENT(fill)) THEN
      Vec = fill
    ELSE
      Vec = 0
    END IF

    Vec(:,:,1:d(3)) = WorkVec

  END SUBROUTINE DoubleBdryArraySize
END SUBROUTINE ParallelRemesh

! Subroutine to create metric for anisotrophic remeshing
! based off a levelset variable
! currently a USF style subroutine as this is how the calving remeshing works
! could change to a subroutine called with remeshing
SUBROUTINE MeshMetricAniso(Model, nodenumber, y, TargetLength)

  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: nodenumber
  REAL(KIND=dp) :: y, TargetLength(:)
  !----------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: TimeVar, LSetVar
  TYPE(Solver_t), POINTER :: Solver
  REAL(KIND=dp) :: t,told, Dist,lc_mindist, lc_maxdist, lc_min, lc_max,s,dx,dz, LevelSet
  INTEGER, POINTER :: LSetPerm(:)
  LOGICAL :: NewTime,FirstTime=.TRUE., Debug, ZInd, Found
  REAL(KIND=dp), POINTER :: LSetValues(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FuncName="MeshMetricAniso", RemeshVarName

  SAVE :: FirstTime, told, Mesh, lc_maxdist, lc_mindist,&
       lc_max, lc_min, dz, newtime, ZInd, RemeshVarName

  Debug = .FALSE.
  Timevar => VariableGet( Model % Variables,'Time')
  t = TimeVar % Values(1)
  Solver => Model % Solver
  SolverParams => Solver % Values

  IF (FirstTime) THEN
    FirstTime = .FALSE.
    NewTime = .TRUE.
    told = t

    Mesh => Model % Mesh

    lc_maxdist = ListGetConstReal(SolverParams, "MeshMetric Max Distance",  Found)
    IF(.NOT. Found) CALL FATAL(FuncName, "Variable: 'MeshMetric Max Distance' not set")
    lc_mindist = ListGetConstReal(SolverParams, "MeshMetric Min Distance",  Found)
    IF(.NOT. Found) CALL FATAL(FuncName, "Variable: 'MeshMetric Min Distance' not set")
    lc_max = ListGetConstReal(SolverParams, "MeshMetric Max LC",  Found)
    IF(.NOT. Found) CALL FATAL(FuncName, "Variable: 'MeshMetric Max LC' not set")
    lc_min = ListGetConstReal(SolverParams, "MeshMetric Min LC",  Found)
    IF(.NOT. Found) CALL FATAL(FuncName, "Variable: 'MeshMetric Min LC' not set")
    ZInd = ListGetLogical(SolverParams, "MeshMetric Z Independent", Found, DefValue = .FALSE.)
    IF(.NOT. Found) CALL INFO(FuncName, "Option 'MeshMetric Z Independent' not found. Default is FALSE")
    IF(ZInd) THEN
      dz = ListGetConstReal(SolverParams, "MeshMetric Vertical LC",  Found)
      IF(.NOT. Found) CALL FATAL(FuncName, "Variable: 'MeshMetric Vertical LC' not set")
    END IF
    RemeshVarName = ListGetString(SolverParams, "Remesh Variable Name", Found)
    IF(.NOT. Found) THEN
      CALL INFO(FuncName, "Remesh Variable Name not found default is Remesh Levelset")
      RemeshVarName = TRIM("Remesh Levelset")
    END IF
  END IF

  IF(t > told) NewTime = .TRUE. !TODO - replace this with Samuel's mesh counter logic
  IF(NewTime) THEN
    told = t
    NewTime = .FALSE.

    Mesh => Model % Mesh

    !mask is the most computationally expense element of this routine
    IF(ASSOCIATED(LSetValues)) NULLIFY(LSetValues)
    IF(ASSOCIATED(LSetPerm)) NULLIFY(LSetPerm)

  END IF

  LSetVar => VariableGet(Mesh % Variables, RemeshVarName, .TRUE., UnfoundFatal=.TRUE.)
  LSetValues => LSetVar % Values
  LSetPerm => LSetVar % Perm

  ! absolute value of the levelset eg distance from zero contour
  LevelSet = ABS(LSetValues(LSetPerm(nodenumber)))

  !Apply caps to this distance:
  Dist = MAX(LevelSet,lc_mindist)
  Dist = MIN(Dist,lc_maxdist)

  s = (Dist - lc_mindist) / (lc_maxdist - lc_mindist)
  dx = s*lc_max + (1-s)*lc_min

  IF(.NOT. ZInd) dz = dx

  TargetLength(1) = dx
  TargetLength(2) = dx
  TargetLength(3) = dz

  IF(Debug) PRINT *,'Node ',nodenumber,'TargetLength: ',TargetLength,' dist: ',Dist,' s: ',s

END SUBROUTINE MeshMetricAniso
