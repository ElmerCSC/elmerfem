
!*****************************************************************************/
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
! * This is a (3D) improvement of FrontalMelt (USF_Frontal) with basic plume shape etc
! * Plume locations are defined by a point in (x,y). When the front moves,
! * plumes are assumed to move normal to the average front parallel.
!
! * Plumes are defined in BC section, e.g.:
! * Plume Count = 1
! * Plume 1 X = 10550.0
! * Plume 1 Y = -96651.0
! *
! * Subsequently modified to make plumes dynamic, taking subglacial discharge
! * and location from GlaDS hydrology solvers as part of HydroCalving system
! ******************************************************************************
! *
! *  Authors: Joe Todd, Samuel Cook
! *  Email:   sc690@cam.ac.uk
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 19.05.2015
! *                 27.06.2018
! *
! ****************************************************************************/

 SUBROUTINE Plume (Model, Solver, dt, TransientSimulation)
   USE Types
   USE CoordinateSystems
   USE DefUtils
   USE ElementDescription

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
   !-----------------------------------
   TYPE(Mesh_t), POINTER :: Mesh, HydroMesh
   TYPE(Solver_t), POINTER :: WorkSolver
   TYPE(Variable_t), POINTER :: ElevVar, ToeCalveVar, WorkVar, WorkVar2,&
                                BMRVar, WorkVar3
   TYPE(ValueList_t), POINTER :: Params, Material
   TYPE(Element_t), POINTER :: Element, Edge
   TYPE(Nodes_t) :: ElementNodes, WorkNodes
   TYPE(GaussIntegrationPoints_t) :: IntegStuff

   REAL(KIND=dp) :: t, s, x, y, m, FrontOrientation(3), PlumeNormal(3), &
        NodeHolder(2,3), PlInt(3), SeaLevel, Melt, prop, hor_dist, plume_z, &
        d0_meltrate, plume_width, BMeltMax, BMeltMME, BMeltDmDz, &
        SqrtElementMetric,U,V,W, &
        Basis(Model % MaxElementNodes), TotalArea, TotalPMelt, TotalBMelt, &
        ElemPMelt, ElemBMelt, ElemToeMelt, Target_PMelt_Average, TotalToeMelt, &
        Target_BMelt_Average, BMelt_Average, PMelt_Average, scale, NodeElev, BMSummerStop, &
        BMSummerStart, Season, aboveMelt, meMelt, Dist, MinDist, ChannelQ,&
        Q0, Plume1MR, Plume2MR, PlProp, Node, NearestNode(3),&
        TargetNode(3), MaxX, MinX, MaxY, MinY, PlDist(2), MeshRes, BMRDist,&
        BMRMinDist, PlDepth, SStart, SStop

   REAL(KIND=dp), ALLOCATABLE :: Xs(:), Ys(:), DwDz(:), W0(:), DmDz(:), MMR(:), MME(:), &
        PlumePoints(:,:,:), PlStart(:,:),PlStop(:,:), PointStore(:),&
        DistArray(:), PlInQ(:), SheetQ(:), PlFinalQ(:), Zi(:), Xi(:), Ta(:),&
        Sa(:), PlAxis(:), PlPos(:,:), NearestFrontNodes(:,:),&
        XArray(:), YArray(:), ZArray(:), PlCoordArray(:,:),&
        Plz(:), PlMR(:), Row(:), PlZArray(:,:), PlMRArray(:,:),&
        TempPlCoordArray(:), TempPlZArray(:), TempPlMRArray(:), MPIArray(:),&
        BMRDistArray(:)
        
   REAL(KIND=dp), POINTER :: PArray(:,:) => NULL(), PArray2(:,:) => NULL(), MeltRate(:),&
        BMeltRate(:) => NULL(), PMeltRate(:) => NULL(),xx(:), yy(:), &
        ToeMeltRate(:), ZPointer(:), MRPointer(:), ZOutput(:), MROutput(:)
   REAL(KIND=dp), ALLOCATABLE, TARGET :: FPOLZ(:), FPOLMR(:)
   INTEGER :: Pl,i,j,k,l,n,dim, PlCount, Active, ierr, ExtrudedLevels, county, &
        Filerows,status(MPI_STATUS_SIZE),TotalNodes,col, counter, aboveIdx, meIdx, &
        NodesPerLevel, TotalPlCount, RowIndex, AxisIndex, myPE, OutputSize,&
        ProcID, MostGL, Output
   INTEGER, ALLOCATABLE :: PlumeOwner(:), idx(:), FNColumns(:), GLNodes(:),&
                           HydroGLNodes(:), SearchIndex(:)
   INTEGER, POINTER :: MeltPerm(:), NodeIndexes(:), ColumnPerm(:), OrderPerm(:)
   INTEGER, PARAMETER :: OutFileUnit = 65, InputFileUnit = 75
   LOGICAL :: Found, Parallel, Boss, Debug, BMeltSwitch, stat,&
              AverageMelt=.FALSE., OutputStats=.FALSE., Visited=.FALSE.,&
              BMFromFile, PlFromFile, RemoveToe, found_intersection, Calving,&
              NoPlume
   LOGICAL, ALLOCATABLE :: PlActive(:), IHavePlume(:), WhoHasPlume(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, PlumeStr, PlumeStrX, &
        PlumeStrY, PlumeStrDwDz, PlumeStrW0, PlumeStrDmDz, PlumeStrMME, PlumeStrMMR, &
        PlumeStrStart, PlumeStrStop, PlumeStrFile, PlumeFile, OutfileName, &
        BMSFile, BMWFile, dumy, ToeCalveVarName, PlMode, BGMode, SaTaZiFile

   TYPE PlumeDefinition_t
      REAL(KIND=dp), POINTER :: z(:), meltrate(:), width(:)
      !NOTE: these are never deallocated, but they don't get reallocated each time, so meh...
   END TYPE PlumeDefinition_t
   TYPE(PlumeDefinition_t), ALLOCATABLE :: ConicalPlumes(:), HydroPlume(:)
   TYPE(PlumeDefinition_t) :: BMeltSummer, BMeltWinter

   SAVE :: Visited, BMeltSummer, BMeltWinter, BMeltSwitch, BMFromFile, BMeltMax, BMeltMME, BMeltDmDz,&
        BMSummerStart, BMSummerStop, ConicalPlumes, DmDz, DwDz, W0, MME, MMR, PlStart, PlStop,&
        PlActive, Xs, Ys, RemoveToe

   INTERFACE
     SUBROUTINE PlumeSolver(Depth, Front, TempA, SalA, Discharge, DepthOutput, MeltOutput, MeshRes, PlDepth)
       USE Types
       REAL(KIND=dp), ALLOCATABLE :: Depth(:), Front(:), TempA(:), SalA(:)
       REAL(KIND=dp) :: Discharge, MeshRes, PlDepth
       REAL(KIND=dp), POINTER :: DepthOutput(:), MeltOutput(:)
     END SUBROUTINE PlumeSolver
   END INTERFACE

   !---------------------------------------------------------
   ! Basic info about solver & simulation
   !---------------------------------------------------------

   Debug = .TRUE.

   SolverName = "Plume"
   Params => GetSolverParams()

   Parallel = (ParEnv % PEs > 1)
   Boss = ParEnv % MyPE == 0

   Mesh => Model % Mesh
   dim = CoordinateSystemDimension()

   Calving = ListGetLogical(Model % Simulation, 'Calving', Found, UnfoundFatal=.FALSE.)
   IF(.NOT. FOUND) Calving = .FALSE.

   t = GetTime()
   season = t - FLOOR(t)

   IF(dim /= 3) CALL Fatal(SolverName, "User Function only works in 3D!")

   PlMode = ListGetString( Params, 'Plume Melt Mode', Found, UnfoundFatal=.FALSE.)
   IF(.NOT. Found) PlMode = 'off'
   PlMode = TRIM(PlMode)

   !Potential to turn melt off
   IF(.NOT. Calving) THEN

     BGMode = ListGetString( Params, 'Background Melt Mode', Found, UnfoundFatal=.TRUE.)
     BGMode = TRIM(BGMode)
   END IF

   IF(.NOT. ASSOCIATED(Solver % Variable)) CALL Fatal(SolverName, "No variable associated!")
   MeltRate => Solver % Variable % Values
   MeltPerm => Solver % Variable % Perm

   ALLOCATE(BMeltRate(SIZE(MeltRate)), PMeltRate(SIZE(MeltRate)))
   BMeltRate = 0.0_dp
   PMeltRate = 0.0_dp
   MeltRate = 0.0_dp

   !Determine mesh levels
   ExtrudedLevels = ListGetInteger(CurrentModel % Simulation,'Extruded Mesh Levels',Found)
   IF(.NOT. Found) ExtrudedLevels = &
        ListGetInteger(CurrentModel % Simulation,'Remesh Extruded Mesh Levels',Found)
   IF(.NOT. Found) CALL Fatal("Remesh",&
        "Unable to find 'Extruded Mesh Levels' or 'Remesh Extruded Mesh Levels'")

   !Modelled melt rates are typically *not* maximum at the base, meaning that
   !progressive melt can result in the formation of a toe. These cause problems
   !for the current calving implementation. Thus, the user can request that these
   !submarine 'toes' be removed by imposing a melt rate which is maximum at the base.
   RemoveToe = ListGetLogical( Params, "Force Toe Calving", Found)
   IF(.NOT. Found) RemoveToe = .FALSE.

   IF(RemoveToe) THEN
     !Get the variable which will hold the 'forced toe calving' distance
     ToeCalveVarName = "Toe Calving"
     ToeCalveVar => VariableGet(Mesh % Variables, ToeCalveVarName, .TRUE.)
     IF(.NOT. ASSOCIATED(ToeCalveVar)) THEN
       WRITE(Message,'(A,A,A)')  "Requested Toe Removal but no variable: '",&
            TRIM(ToeCalveVarName),"' found."
       CALL Fatal(SolverName, Message)
     END IF

     !Remesh gives Exported variables their own perms. Revert this - same perm is useful
     IF(.NOT. ASSOCIATED(ToeCalveVar % Perm, MeltPerm)) THEN
       DEALLOCATE(ToeCalveVar % Perm)
       ToeCalveVar % Perm => MeltPerm
     END IF

     ALLOCATE(ToeMeltRate(SIZE(MeltRate)))
     ToeMeltRate = 0.0_dp

     !Need global mesh structure info
     IF(Parallel) THEN
       !Rather than summing NoNodes from each part, we simply find
       !the maximum global node number
       CALL MPI_AllReduce(MAXVAL(Mesh % ParallelInfo % GlobalDOFs), TotalNodes, &
            1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD,ierr)
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     ELSE
       TotalNodes = Mesh % NumberOfNodes
     END IF

     NodesPerLevel = TotalNodes / ExtrudedLevels

     ALLOCATE(FNColumns(Mesh % NumberOfNodes))
     FNColumns = MOD(Mesh % ParallelInfo % GlobalDOFs, NodesPerLevel)

     !PRINT *,ParEnv % MyPE, ' debug, nodesperlevel, totalnodes, extrudedlevels: ',&
     !     NodesPerLevel, TotalNOdes, ExtrudedLevels
   END IF

   !Get the elevation variable
   ElevVar => VariableGet(Mesh % Variables, "Elevation", .TRUE.)
   IF(.NOT. ASSOCIATED(ElevVar)) CALL Fatal(SolverName,"Couldn't find 'Elevation' variable &
        &needed to compute background melt rate")

   !Get the orientation of the calving front
   !TODO: generalize and link
   IF(.NOT. Calving) THEN
     PArray => ListGetConstRealArray( Params,'Front Orientation', Found, UnfoundFatal=.TRUE.)
     DO i=1,3
       FrontOrientation(i) = PArray(i,1)
     END DO

     !We can define the plume's position by an initial point (Xs,Ys) and a
     !vector (in line with the front) which describes its movement back and forward
     !as the front migrates. In this case, the plume is a vertical line
     !which moves forwards and backwards (with the front) along the FrontOrientation
     !vector. Thus, the normal vector which define's the plume's plane is the reciprocal
     !of the FrontOrientation vector (in 2D):
     PlumeNormal(1) = FrontOrientation(2)
     PlumeNormal(2) = -1.0_dp * FrontOrientation(1)
     PlumeNormal(3) = 0.0_dp
   END IF
  
   Material => GetMaterial()
   SeaLevel = GetCReal(Material, 'Sea Level', Found)
   IF(.NOT. Found) SeaLevel = 0.0_dp
   !PRINT *, 'P1',ParEnv % myPE
   !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

   IF(Calving) THEN
     !Find all GL nodes on hydromesh
     DO i=1,Model % NumberOfSolvers
       WorkSolver => Model % Solvers(i)
       IF(WorkSolver % Variable % Name == 'hydraulic potential') THEN
         HydroMesh => WorkSolver % Mesh
         EXIT
       END IF
     END DO
     
     WorkVar => VariableGet(HydroMesh % Variables, 'groundedmask', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
     j=1
     ALLOCATE(GLNodes(HydroMesh % NumberOfNodes))
     DO i=1,HydroMesh % NumberOfNodes
       IF(WorkVar % Values(WorkVar % Perm(i)) == 0.0) THEN
         GLNodes(j) = i
         j=j+1
       ELSE
         CYCLE
       END IF
     END DO
     IF(j>1) THEN
       ALLOCATE(HydroGLNodes(j-1))
       HydroGLNodes(1:SIZE(HydroGLNodes)) = GLNodes(1:SIZE(HydroGLNodes))
     ELSE
       ALLOCATE(HydroGLNodes(1))
       HydroGLNodes(1) = 0
       PlCount = 0
     END IF
     DEALLOCATE(GLNodes)

     !For each hydro GL node, find nearest calving front node
     MinDist = 100000.0
     BMRMinDist = 100000.0
     ALLOCATE(NearestFrontNodes(SIZE(HydroGLNodes), 5),&
             DistArray(ParEnv % PEs), XArray(ParEnv % PEs),&
             YArray(ParEnv % PEs), ZArray(ParEnv % PEs), &
             BMRDistArray(ParEnv % PEs))
     WorkVar => VariableGet(HydroMesh % Variables, 'zb', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
     BMRVar => VariableGet(Model % Mesh % Variables, 'gmcheck', ThisOnly=.TRUE., UnfoundFatal=.FALSE.)
     XArray = 0.0_dp
     YArray = 0.0_dp
     ZArray = 0.0_dp
     DistArray = 0.0_dp
     BMRDistArray = 0.0_dp
     MostGL = 0
     CALL MPI_ALLREDUCE(SIZE(HydroGLNodes), MostGL, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)

     !Cycle through all nodes in all processes
     DO i=1, ParEnv % PEs
       myPE = i-1
       DO j=1,MostGL !Need to have all processes running same length loop
         MinDist = 100000.0
         BMRMinDist = 100000.0
         IF(j>SIZE(HydroGLNodes)) THEN
           TargetNode(:) = 0.0
         ELSEIF(HydroGLNodes(j) == 0) THEN
           TargetNode(:) = 0.0
         ELSE
           TargetNode(1) = HydroMesh % Nodes % x(HydroGLNodes(j))
           TargetNode(2) = HydroMesh % Nodes % y(HydroGLNodes(j))
           TargetNode(3) = WorkVar % Values(WorkVar % Perm((HydroGLNodes(j))))
         END IF
         CALL MPI_BCAST(TargetNode, 3, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)
         DO k=1,Solver % Mesh % NumberOfNodes
           IF(ALL(TargetNode(1:3)==0.0)) THEN
             NearestNode(1) = 0.0
             NearestNode(2) = 0.0
             NearestNode(3) = 0.0
             MinDist = 99999.0
             EXIT
           END IF
           IF(MeltPerm(k)<=0) CYCLE
           Dist = (Solver % Mesh % Nodes % x(k) - TargetNode(1))**2
           Dist = Dist + (Solver % Mesh % Nodes % y(k) - TargetNode(2))**2
           Dist = Dist + (Solver % Mesh % Nodes % z(k) - TargetNode(3))**2
           Dist = SQRT(Dist)
           IF(Dist<MinDist) THEN
             MinDist = Dist
             NearestNode(1) = Solver % Mesh % Nodes % x(k)
             NearestNode(2) = Solver % Mesh % Nodes % y(k)
             NearestNode(3) = Solver % Mesh % Nodes % z(k)
           END IF
         END DO
         DO k=1,Solver % Mesh % NumberOfNodes
           IF(ALL(TargetNode(1:3)==0.0)) THEN
             BMRDist = 99999.0
             EXIT
           END IF
           IF(.NOT. ASSOCIATED(BMRVar)) THEN
             BMRMinDist = 0
             CYCLE
           END IF
           IF(BMRVar % Perm(k)<=0) CYCLE
           IF(BMRVar % Values(BMRVar % Perm(k)) < 1.0) CYCLE
           Dist = (Solver % Mesh % Nodes % x(k) - TargetNode(1))**2
           Dist = Dist + (Solver % Mesh % Nodes % y(k) - TargetNode(2))**2
           Dist = SQRT(Dist)
           IF(Dist<BMRMinDist) THEN
             BMRMinDist = Dist
           END IF
         END DO

         CALL MPI_GATHER(MinDist, 1, MPI_DOUBLE_PRECISION, DistArray, 1, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)
         CALL MPI_GATHER(BMRMinDist, 1, MPI_DOUBLE_PRECISION, BMRDistArray, 1, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)
         CALL MPI_GATHER(NearestNode(1), 1, MPI_DOUBLE_PRECISION, XArray, 1, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)
         CALL MPI_GATHER(NearestNode(2), 1, MPI_DOUBLE_PRECISION, YArray, 1, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)
         CALL MPI_GATHER(NearestNode(3), 1, MPI_DOUBLE_PRECISION, ZArray, 1, MPI_DOUBLE_PRECISION, myPE, ELMER_COMM_WORLD, ierr)

         IF(ParEnv % myPE == myPE) THEN
           IF(j > SIZE(HydroGLNodes)) CYCLE
           MinDist = 100000
           DO k=1, SIZE(DistArray)
             IF(DistArray(k) < MinDist) THEN
               ProcID = (k-1)
               MinDist = DistArray(k)
             END IF
           END DO
           BMRMinDist = 100000
           DO k=1, SIZE(BMRDistArray)
             IF(BMRDistArray(k) < BMRMinDist) THEN
               BMRMinDist = BMRDistArray(k)
             END IF
           END DO
           NearestFrontNodes(j, 1) = XArray(ProcID+1)
           NearestFrontNodes(j, 2) = YArray(ProcID+1)
           NearestFrontNodes(j, 3) = ZArray(ProcID+1)
           NearestFrontNodes(j, 4) = MinDist
           NearestFrontNodes(j, 5) = BMRMinDist
           !This is for the case where nodes on the two meshes exactly coincide
           !and a 'depth' of 0 is returned
           IF(NearestFrontNodes(j,3) > -2.0) NearestFrontNodes(j,3) = TargetNode(3)
         END IF
       END DO !j
     END DO !i
     DEALLOCATE(DistArray, BMRDistArray, XArray, YArray, ZArray)
     !PRINT *, 'P2',ParEnv % myPE
     !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

     !If distance < mesh res, continue
     !If distance > mesh res, but < (say) 500m, also continue (probably means
     !you've got an ungrounded area that hasn't calved, but fine to add to
     !plume
     !However, ungrounded areas inland effectively remove hydrology between them
     !and actual GL, so need to include channels and sheet feeding them
     !Thus, set critical distance to, say, 4000m to get all water flowing to
     !terminus, but ignore any weird GL artefacts further inland
     MeshRes = GetCReal(Params, 'Mesh Resolution', Found)
     IF(.NOT. Found) THEN
       CALL Info(SolverName, 'Mesh resolution not specified; setting to 250')
       MeshRes = 250.0
     END IF
     DO i=1,SIZE(NearestFrontNodes,1)
       IF(NearestFrontNodes(i,4)>500.0 .AND.&
         NearestFrontNodes(i,5)>(MeshRes+(MeshRes/10))) THEN
         NearestFrontNodes(i,1:5) = 0.0
         HydroGLNodes(i) = 0
       END IF
     END DO
     NULLIFY(BMRVar)

     !Calculate nodal discharge and find GL node with highest discharge
     !(channel + sheet). As things stand, this will double count channels
     !running along the GL itself (i.e., edges with both nodes on the GL),
     !but these a) should be prohibited by the BCs and b) are probably
     !negligible anyway
     ALLOCATE(PlInQ(SIZE(HydroGLNodes)), SheetQ(3))
     WorkVar => VariableGet(HydroMesh % Variables, 'channel flux', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
     WorkVar2 => VariableGet(HydroMesh % Variables, 'sheet discharge', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
     WorkVar3 => VariableGet(HydroMesh % Variables, 'sheet thickness', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)
     j=1
     DO i=1, SIZE(HydroGLNodes)
       SheetQ = 0.0_dp
       ChannelQ = 0.0
       IF(HydroGLNodes(i) == 0) THEN
         PlInQ(i) = 0.0
         CYCLE
       END IF
       DO j=1,HydroMesh % NumberOfEdges
         Edge => HydroMesh % Edges(j)
         IF(ANY(Edge % NodeIndexes(1:2) == HydroGLNodes(i))) THEN
           ChannelQ = ChannelQ + WorkVar % Values(WorkVar % Perm(HydroMesh % NumberOfNodes+j))
         ELSE
           CYCLE
         END IF
       END DO
       DO j=1,2
         SheetQ(j) = SheetQ(j) + (WorkVar2 % Values(2*(WorkVar2 % Perm(HydroGLNodes(i))-1)+j))
       END DO
       SheetQ(3) = SQRT((SheetQ(1)**2)+(SheetQ(2)**2))*WorkVar3 % Values(WorkVar3 % Perm(HydroGLNodes(i)))
       PlInQ(i) = ChannelQ + SheetQ(3)
     END DO

     !Check for multiple entries that have same NearestFrontNode and combine Q
     k=1
     PlCount = 0
     ALLOCATE(PlPos(SIZE(NearestFrontNodes,1),3), PlFinalQ(SIZE(PlInQ)))
     PlPos = 0.0_dp
     PlFinalQ = 0.0_dp
     DO i=1, SIZE(NearestFrontNodes,1)
       IF(HydroGLNodes(i) == 0.0) CYCLE
       DO j=1, SIZE(NearestFrontNodes,1)
         IF(i==j) CYCLE
         IF(HydroGLNodes(j) == 0.0) CYCLE
         IF(ALL(ANINT(NearestFrontNodes(j,1:3)) == ANINT(NearestFrontNodes(i,1:3)))) THEN
           PlInQ(i) = PlInQ(i) + PlInQ(j)
           PlInQ(j) = 0.0
           HydroGLNodes(j) = 0
           NearestFrontNodes(j,1:5) = 0
         END IF
       END DO
       PlPos(k,1:3) = NearestFrontNodes(i,1:3)
       PlFinalQ(k) = PlInQ(i)
       !Discharge needs to be in m2/s, so divide by some guess at the outlet
       !width. Could use channel area from GLaDS, but this assumes semi-circle
       !R channels. Jackson et al. (2017) show line/wedge plumes from low and 
       !broad outlets are likely what Greenlandic tidewater glaciers display.
       !Just need to take care of time units here (m3>m2 done later)
       PlFinalQ(k) = PlFinalQ(k)/(365.25*24*60*60)
       k=k+1
       PlCount = PlCount + 1
     END DO  
     !To stop model bothering if no plume with enough discharge to actually do
     !anything
     !IF(PlCount > 0 .AND. MAXVAL(PlFinalQ)<1E-4) THEN
       !PlCount = 0
     !END IF
     !PRINT *, 'P3',ParEnv % myPE
     !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

     !Link everything up to the output array at relevant points further down
     !Model (at moment) has vertical calving front, so, for plume model, xi is
     !constant.
     !Zi should start at depth of bottom of calving front where plume is sited
     !and then be defined at evenly-spaced increments up to 0 depth (i.e. the 
     !water surface).
     !Ta and Sa read in from file and then use FindPointOnLine to interpolate
     !values for them at that plume's set of zi points
     !Q0 taken straight from hydrology - will be total discharge
     !Then run plume model.
     Output = 0
     CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
     IF(PlCount > 0) THEN
       ALLOCATE(HydroPlume(PlCount))
       SELECT CASE(PlMode)
       CASE("seasonal")
         SStart = ListGetConstReal( Params, "Plume Melt Summer Start", UnfoundFatal=.TRUE.)
         SStop = ListGetConstReal( Params, "Plume Melt Summer Stop", UnfoundFatal=.TRUE.)
         IF(season > SStart .AND. season < SStop) THEN
           SaTaZiFile = ListGetString( Params, "Summer Salinity Temp Depth Input File", UnfoundFatal=.TRUE.)
         ELSE
           SaTaZiFile = ListGetString( Params, "Winter Salinity Temp Depth Input File", UnfoundFatal=.TRUE.)
         END IF
       CASE("constant")
         SaTaZiFile = ListGetString( Params, "Salinity Temp Depth Input File", UnfoundFatal=.TRUE.)
       CASE("off")
         RETURN
       CASE DEFAULT
         PRINT *, 'Please specify what mode you want plumes to operate in.&
                  Options are seasonal, constant or off'
         RETURN
       END SELECT

       OPEN(UNIT=InputFileUnit, FILE=SaTaZiFile, IOSTAT=ierr)

       !Check file length
       Filerows = 0
       DO WHILE(.TRUE.)
         READ(InputFileUnit, *,IOSTAT=ierr) dumy
         IF(ierr /= 0) EXIT
         Filerows = Filerows + 1
       END DO

       REWIND(InputFileUnit)
       ALLOCATE( Zi(Filerows), Sa(Filerows), Ta(Filerows), Xi(Filerows))

       DO i=1,Filerows
         READ(InputFileUnit, *) Zi(i), Sa(i), Ta(i)
         IF(Zi(i) > 0.0) Zi(i) = Zi(i)*(-1)
       END DO
       CLOSE(InputFileUnit)

       DO Pl=1,PlCount
         PlDepth = PlPos(Pl, 3)
 
         !Assume vertical calving front. Internal mesh extrusion kind of imposes
         !it, so unlikely to change this soon.
         Xi(:) = 0.0
         !Divide discharge by nominal frontal mesh resolution to get m2/s
         Q0 = PlFinalQ(Pl)/MeshRes

         ALLOCATE(HydroPlume(Pl) % z(SIZE(Zi)),&
                  HydroPlume(Pl) % meltrate(SIZE(Zi)))
         HydroPlume(Pl) % z = 0.0_dp
         HydroPlume(Pl) % meltrate = 0.0_dp

         ZOutput => HydroPlume(Pl) % z
         MROutput => HydroPlume(Pl) % meltrate
         IF(Q0 .LE. 0.0 .OR. PlDepth > -10.0) THEN
           !This avoids the solver failing on Q=0 and spewing NaNs everywhere
           !and also potential weird geometry edge cases where the PlDepth ends
           !up zero or positive
           MROutput(:) = 0.0
           ZOutput(:) = Zi(:)
         ELSE
           !This solves for a line/wedge plume of assumed width of MeshRes
           !Truncated line/wedge plumes seem better fit for combined channel+sheet
           !discharge and for Jackson et al. (2017)'s observations
           CALL PlumeSolver(Zi, Xi, Ta, Sa, Q0, ZOutput, MROutput, MeshRes, PlDepth)
         END IF
         !DO i=1,SIZE(MROutput)
           !IF(ISNAN(MROutput(i))) PRINT *, 'Found a NaN ',Q0,PlDepth,PlPos(Pl,:)
         !END DO
         IF(SIZE(ZOutput)>Output) THEN
           Output = SIZE(ZOutput)
         END IF
       END DO
       IF(Output == 0) PlCount = 0
     END IF !PlCount > 0
     !PRINT *, 'P4',ParEnv % myPE
     !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
     
     !Work out whether x or y axis is major frontal axis
     !Allocate TotalPlCount and OutputSize to give fixed size for MPI - ensures
     !have enough space in every partition for any number of plumes
     CALL MPI_ALLREDUCE(PlCount, TotalPlCount, 1, MPI_INTEGER, MPI_SUM, ELMER_COMM_WORLD, ierr)
     IF(TotalPlCount > 0) THEN
       CALL MPI_ALLREDUCE(Output, OutputSize, 1, MPI_INTEGER, MPI_MAX, ELMER_COMM_WORLD, ierr)
       ALLOCATE(PlAxis(TotalPlCount), PlZ(OutputSize*TotalPlCount),&
               PlMR(OutputSize*TotalPlCount), PlActive(TotalPlCount))       
       ALLOCATE(TempPlCoordArray(TotalPlCount*ParEnv % PEs))
       ALLOCATE(TempPlZArray(OutputSize*TotalPlCount*ParEnv % PEs))
       ALLOCATE(TempPlMRArray(OutputSize*TotalPlCount*ParEnv % PEs))
     
       PlActive(:) = .TRUE.
       PlZ(:) = 9999.0
       TempPlZArray(:) = 9999.0
       PlMR(:) = -10000.0
       TempPlMRArray(:) = -1.0
       PlAxis(:) = 0.0_dp

       MaxX = -1E16
       MaxY = -1E16
       MinX = 1E16
       MinY = 1E16
       DO i=1, SIZE(PlPos,1)
         IF(PlPos(i,1) == 0.0) CYCLE
         IF(PlPos(i,1)>MaxX) MaxX = PlPos(i,1)
         IF(PlPos(i,1)<MinX) MinX = PlPos(i,1)
         IF(PlPos(i,2)>MaxY) MaxY = PlPos(i,2)
         IF(PlPos(i,2)<MinY) MinY = PlPos(i,2)
       END DO
       x = MaxX - MinX
       y = MaxY - MinY

       IF(PlCount>0) THEN
         IF(x>y) THEN
           PlAxis(1:PlCount) = PlPos(1:PlCount,1)
           AxisIndex = 1
         ELSE
           PlAxis(1:PlCount) = PlPos(1:PlCount,2)
           AxisIndex = 2
         END IF

         DO i=1,PlCount
           ZOutput => HydroPlume(i) % z
           MROutput => HydroPlume(i) % meltrate
           PlZ(1+((i-1)*OutputSize):OutputSize*i) = ZOutput(1:OutputSize)
           PlMR(1+((i-1)*OutputSize):OutputSize*i) = MROutput(1:OutputSize)
         END DO
       ELSE 
         PlAxis(:) = 0.0
         PlZ(:) = 9999.0
         PlMR(:) = -10000.0
       END IF
       !PRINT *, 'P5',ParEnv % myPE
       !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
     
       !MPI stuff to make everything available to all processes - send z and MR 
       !to separate arrays using PlCount as row index, I think
       CALL MPI_GATHER(PlAxis(:), TotalPlCount, MPI_DOUBLE_PRECISION,&
            TempPlCoordArray(:), TotalPlCount, MPI_DOUBLE_PRECISION, 0,&
            ELMER_COMM_WORLD, ierr)
       CALL MPI_GATHER(PlZ(:), OutputSize*TotalPlCount, MPI_DOUBLE_PRECISION,&
            TempPlZArray(:), OutPutSize*TotalPlCount, MPI_DOUBLE_PRECISION, 0,&
            ELMER_COMM_WORLD, ierr)
       CALL MPI_GATHER(PlMR(:), OutputSize*TotalPlCount, MPI_DOUBLE_PRECISION,&
            TempPlMRArray(:), OutPutSize*TotalPlCount, MPI_DOUBLE_PRECISION, 0,&
            ELMER_COMM_WORLD, ierr)

       ALLOCATE(PlCoordArray(TotalPlCount, 2), PlZArray(OutputSize,TotalPlCount),&
               PlMRArray(OutputSize,TotalPlCount))

       IF(ParEnv % myPE == 0) THEN
         j=1
         DO i=1,SIZE(TempPlCoordArray)
           IF(TempPlCoordArray(i)<1E-16 .AND. TempPlCoordArray(i)>-1E-16) CYCLE
           PlCoordArray(j,1) = TempPlCoordArray(i)
           PlCoordArray(j,2) = j
           j=j+1
         END DO
         j=1
         k=1
         l=1
         DO i=1,SIZE(TempPlZArray)
           IF(TempPlZArray(i)==-10000.0 .OR. TempPlZArray(i) == 9999.0) CYCLE
           m = (j-1)/OutputSize
           l = AINT(m)+1
           k = j-((l-1)*OutputSize)
           j = j+1
           PlZArray(k,l) = TempPlZArray(i)
           PlMRArray(k,l) = TempPlMRArray(i)
         END DO
         !CoordArray needs sorting by co-ordinate to order plumes directionally
         !Second column of array is then effectively plume perm, and can be used
         !as index to access correct column in Z and MR arrays, so these can be
         !unsorted
         ALLOCATE(Row(2))
         DO i=1, TotalPlCount
           RowIndex = MINLOC(PlCoordArray(i:TotalPlCount,1), DIM=1)+i-1
           Row(:) = PlCoordArray(i,:)
           PlCoordArray(i,:) = PlCoordArray(RowIndex,:)
           PlCoordArray(RowIndex,:) = Row(:)
         END DO
         DEALLOCATE(Row)
       END IF
       !PRINT *, 'P6',ParEnv % myPE
       !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

       !MPI call to send full final plume arrays to every partition
       CALL MPI_BCAST(PlZArray, OutputSize*TotalPlCount, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
       CALL MPI_BCAST(PlMRArray, OutputSize*TotalPlCount, MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
       CALL MPI_BCAST(PlCoordArray, SIZE(PlCoordArray), MPI_DOUBLE_PRECISION, 0, ELMER_COMM_WORLD, ierr)
       CALL MPI_BCAST(TotalPlCount, 1, MPI_INTEGER, 0, ELMER_COMM_WORLD, ierr)
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
       DEALLOCATE(TempPlCoordArray, TempPlZArray, TempPlMRArray, Plz, PlMR)
     ELSE
       ALLOCATE(PlActive(1))
       PlActive(1) = .FALSE.
     END IF !TotalPlCount > 0
   ELSE 
     PlCount = ListGetInteger( Params, "Plume Count", UnfoundFatal=.TRUE.)
   END IF

   PlFromFile = ListGetLogical( Params, "Plumes From File", Found)
   IF(.NOT. Found) PlFromFile = .FALSE.

   BMFromFile = ListGetLogical( Params, "Background Melt From File", Found)
   IF(.NOT. Found) BMFromFile = .FALSE.
   AverageMelt = GetLogical( Params, "Scale Melt To Average", Found)

   IF(.NOT. Calving) THEN
     IF(.NOT. Visited) THEN
       ALLOCATE(Xs(PlCount), Ys(PlCount), &
          PlActive(PlCount), PlStart(PlCount, 10),&
          PlStop(PlCount, 10))
       PlActive = .FALSE.

       !Set points to NaN - fill them later by checking for potential point
       PlStart = 0.0_dp
       PlStart = PlStart / PlStart
       PlStop = 0.0_dp
       PlStop = PlStop / PlStop

       IF(PlFromFile) THEN
         ALLOCATE(ConicalPlumes(PlCount))
       ELSE
         ALLOCATE(DmDz(PlCount), DwDz(PlCount), MMR(PlCount), MME(PlCount), W0(PlCount))
       END IF
     END IF

     !Leaving this unSAVEd, potential to change number of extruded levels?
     ALLOCATE(PlumePoints(PlCount, ExtrudedLevels, 3))
     !NaN, see above
     PlumePoints = 0.0_dp
     PlumePoints = PlumePoints / PlumePoints

     !Get output file name for stats
     OutfileName = ListGetString(Params,"Melt Stats File", OutputStats, UnfoundFatal=.TRUE.)

     !Check for background melt
     SELECT CASE(BGMode)
     CASE("off")
       BMeltSwitch = .FALSE.
     CASE("seasonal")
       BMeltSwitch = .TRUE.
     CASE Default
       CALL Fatal(SolverName, "Invalid Background Melt Mode, valid options are 'seasonal' and 'off'.")
     END SELECT

     !------------------------------------------------------
     !Determine how the planar (background) plume is defined
     !------------------------------------------------------
     IF(BMeltSwitch .AND. .NOT. Visited) THEN

       IF(BMFromFile) THEN

         BMSFile = ListGetString( Params, "Background Melt Summer Input File", UnfoundFatal=.TRUE.)
         BMWFile = ListGetString( Params, "Background Melt Winter Input File", UnfoundFatal=.TRUE.)

         BMSummerStart = ListGetConstReal( Params, "Background Melt Summer Start", UnfoundFatal=.TRUE.)
         BMSummerStop = ListGetConstReal( Params, "Background Melt Summer Stop", UnfoundFatal=.TRUE.)

         OPEN(UNIT=InputFileUnit, FILE=BMSFile, IOSTAT=ierr)

         !Check file length
         Filerows = 0
         DO WHILE(.TRUE.)
           READ(InputFileUnit, *,IOSTAT=ierr) dumy
           IF(ierr /= 0) EXIT
           Filerows = Filerows + 1
         END DO

         REWIND(InputFileUnit)
         ALLOCATE( BMeltSummer % z(Filerows),&
            BMeltSummer % meltrate(Filerows))

         DO i=1,Filerows
           READ(InputFileUnit, *) BMeltSummer % z(i), BMeltSummer % meltrate(i)
           IF(Boss) PRINT *,TRIM(SolverName),': Summer z, melt: ', BMeltSummer % z(i), &
              BMeltSummer % meltrate(i)
         END DO
         CLOSE(InputFileUnit)

         !Check plume definition monotonically increasing in z
         DO i=2, FileRows
           IF(BMeltSummer % z(i) <= BMeltSummer % z(i-1)) THEN
             CALL Fatal(SolverName, "Background Summer Plume definition should be &
                &monotonically increasing in z.")
           END IF
         END DO

         OPEN(UNIT=InputFileUnit, FILE=BMWFile, IOSTAT=ierr)

         !Check file length
         Filerows = 0
         DO WHILE(.TRUE.)
           READ(InputFileUnit, *,IOSTAT=ierr) dumy
           IF(ierr /= 0) EXIT
           Filerows = Filerows + 1
         END DO

         REWIND(InputFileUnit)
         ALLOCATE( BMeltWinter % z(Filerows),&
            BMeltWinter % meltrate(Filerows))

         DO i=1,Filerows
           READ(InputFileUnit, *) BMeltWinter % z(i), BMeltWinter % meltrate(i)
           IF(Boss) PRINT *,TRIM(SolverName),': Winter z, melt: ', BMeltWinter % z(i),&
              BMeltWinter % meltrate(i)
         END DO
         CLOSE(InputFileUnit)

         !Check plume definition monotonically increasing in z
         DO i=2, FileRows
           IF(BMeltWinter % z(i) <= BMeltWinter % z(i-1)) THEN
             CALL Fatal(SolverName, "Background Winter Plume definition should be &
                &monotonically increasing in z.")
           END IF
         END DO

       ELSE
         BMeltMax = GetConstReal( Params, "Background Melt Max", Found)
         IF(.NOT. Found) CALL Fatal(SolverName, "Background melt requested but no&
            &'Background Melt Max' Found")

         BMeltMME = GetConstReal( Params, "Background Melt Max Melt Elevation", Found)
         IF(.NOT. Found) CALL Fatal(SolverName, "Background melt requested but no&
            &'Background Melt Max Melt Elevation' Found")

         BMeltDmDz = GetConstReal( Params, "Background Melt DmDz", Found)
         IF(.NOT. Found) CALL Fatal(SolverName, "Background melt requested but no&
            &'Background Melt DmDz' Found")
       END IF
     END IF

     AverageMelt = GetLogical( Params, "Scale Melt To Average", Found)
     IF(Found .AND. AverageMelt) THEN
       Target_BMelt_Average = GetConstReal( Params, "Average Background Melt Rate", Found)
       IF(.NOT. Found) CALL Fatal(SolverName, &
          "Requested 'Scale Melt To Average' but no 'Average Background Melt Rate' found.")
       Target_PMelt_Average = GetConstReal( Params, "Average Plume Melt Rate", Found)
       IF(.NOT. Found) CALL Fatal(SolverName, &
          "Requested 'Scale Melt To Average' but no 'Average Plume Melt Rate' found.")
     ELSE
       AverageMelt = .FALSE.
     END IF

     !Find the plume definitions
     IF(.NOT. Visited) THEN
       DO Pl=1, PlCount

         WRITE(PlumeStr,'(a,i0)') "Plume ",Pl

         WRITE(PlumeStrX,'(a,a)') TRIM(PlumeStr), " X"
         WRITE(PlumeStrY,'(a,a)') TRIM(PlumeStr), " Y"
         WRITE(PlumeStrStart,'(a,a)') TRIM(PlumeStr), " Start Times"
         WRITE(PlumeStrStop,'(a,a)') TRIM(PlumeStr), " Stop Times"

         Xs(Pl) = ListGetConstReal( Params, PlumeStrX, Found, UnfoundFatal=.TRUE.)
         Ys(Pl) = ListGetConstReal( Params, PlumeStrY, Found, UnfoundFatal=.TRUE.)

         PArray => ListGetConstRealArray( Params,PlumeStrStart, Found, UnfoundFatal=.TRUE.)
         PArray2 => ListGetConstRealArray( Params,PlumeStrStop, Found, UnfoundFatal=.TRUE.)

         IF(SIZE(PArray) /= SIZE(PArray2)) &
            CALL Fatal(SolverName, "Plume start and stop array size mismatch")

         DO i=1, SIZE(Parray,1)
           PlStart(Pl,i) = PArray(i,1)
           PlStop(Pl,i) = PArray2(i,1)
         END DO

         IF(PlFromFile) THEN
           WRITE(PlumeStrFile,'(a,a)') TRIM(PlumeStr), " Input File"
           PlumeFile = ListGetString(Params, PlumeStrFile, UnfoundFatal=.TRUE.)

           OPEN(UNIT=InputFileUnit, FILE=PlumeFile, IOSTAT=ierr)
           IF(ierr /= 0) CALL Fatal(SolverName, "Failed to open plume input file")

           !Check file length
           Filerows = 0
           DO WHILE(.TRUE.)
             READ(InputFileUnit, *,IOSTAT=ierr) dumy
             IF(ierr /= 0) EXIT
             Filerows = Filerows + 1
           END DO

           REWIND(InputFileUnit)
           ALLOCATE( ConicalPlumes(Pl) % z(Filerows),&
              ConicalPlumes(Pl) % meltrate(Filerows),&
              ConicalPlumes(Pl) % width(Filerows))

           DO i=1,Filerows
             READ(InputFileUnit, *) ConicalPlumes(Pl) % z(i), &
                ConicalPlumes(Pl) % meltrate(i), &
                ConicalPlumes(Pl) % width(i)

             PRINT *,'Debug, plume info: ',i,ConicalPlumes(Pl) % z(i), &
                ConicalPlumes(Pl) % meltrate(i), &
                ConicalPlumes(Pl) % width(i)
           END DO
           CLOSE(InputFileUnit)

           !Check plume definition monotonically increasing in z
           DO i=2, FileRows
             IF(ConicalPlumes(Pl) % z(i) <= ConicalPlumes(Pl) % z(i-1)) THEN
               CALL Fatal(SolverName, "Plume definition should be monotonically increasing in z.")
             END IF
           END DO

         ELSE

           WRITE(PlumeStrDwDz,'(a,a)') TRIM(PlumeStr), " DwDz"
           WRITE(PlumeStrW0,'(a,a)') TRIM(PlumeStr), " Initial Width"
           WRITE(PlumeStrDmDz,'(a,a)') TRIM(PlumeStr), " DmDz"
           WRITE(PlumeStrMMR,'(a,a)') TRIM(PlumeStr), " Max Melt Rate"
           WRITE(PlumeStrMME,'(a,a)') TRIM(PlumeStr), " Max Melt Elevation"

           Dmdz(Pl) = ListGetConstReal( Params, PlumeStrDmdz, Found, UnfoundFatal=.TRUE.)
           DwDz(Pl) = ListGetConstReal( Params, PlumeStrDwDz, Found, UnfoundFatal=.TRUE.)
           W0(Pl) = ListGetConstReal( Params, PlumeStrW0, Found, UnfoundFatal=.TRUE.)
           MME(Pl) = ListGetConstReal( Params, PlumeStrMME, Found, UnfoundFatal=.TRUE.)
           MMR(Pl) = ListGetConstReal( Params, PlumeStrMMR, Found, UnfoundFatal=.TRUE.)
         END IF
       END DO
     END IF

     DO Pl=1,PlCount
       SELECT CASE(PlMode)
       CASE("seasonal")
         DO i=1, SIZE(PlStart,2)
           IF(ISNAN(PlStart(Pl,i)) .OR. ISNAN(PlStop(Pl,i))) EXIT

           IF((season > PlStart(Pl,i)) .AND. (season < PlStop(Pl,i))) THEN
             PlActive(Pl) = .TRUE.
             EXIT
           ELSE
             PlActive(Pl) = .FALSE.
           END IF
         END DO
         PRINT *,'Debug, plume ',Pl,' is active: ',PlActive(Pl)
       CASE("off")
         PlActive(Pl) = .FALSE.
       CASE DEFAULT
         CALL Fatal(SolverName, "Unknown plume melt mode, valid options are 'seasonal' and 'off'.")
       END SELECT
     END DO


     !----------------------------------------------------------
     !  Strategy:
     !
     !  Based on x,y position of plume, find and MPI the projection
     !  of that position onto the front.
     !
     !  In other words, pass a series of coordinates which define the
     !  plumes journey 'up the front'

     IF(Parallel) ALLOCATE(IHavePlume(PlCount), WhoHasPlume(PlCount*ParEnv % PEs), PlumeOwner(PlCount))

     ALLOCATE(NodeIndexes(4)) !Assume working with side of extruded mesh

     IHavePlume = .FALSE.
     WhoHasPlume = .FALSE.
     PlumeOwner = -1

     !------------------------------------------------------------------
     !Find all the element intersections for the each plume's centreline
     !If none, this plume doesn't exist in this partition.
     !------------------------------------------------------------------

     DO Pl=1, PlCount

        Found = .FALSE.
        county = 0

        Active = GetNOFActive()
        DO j=1,Active
           Element => GetActiveElement(j)
           IF(Element % TYPE % ElementCode == 101) CYCLE
           IF(Element % TYPE % ElementCode /= 404) &
              CALL Fatal(SolverName, 'Found a non-404 element, this solver assumes extruded mesh')

           !Determine the two pairs of nodes which each make up a horizontal bar
           !Extruded mesh means lower two node indexes = lower two nodes
           NodeIndexes = Element % NodeIndexes
           CALL Sort( 4, NodeIndexes )

           DO k=1,2

              NodeHolder(1,1) = Mesh % Nodes % x(NodeIndexes(k*2 - 1))
              NodeHolder(1,2) = Mesh % Nodes % y(NodeIndexes(k*2 - 1))
              NodeHolder(1,3) = Mesh % Nodes % z(NodeIndexes(k*2 - 1))

              NodeHolder(2,1) = Mesh % Nodes % x(NodeIndexes(k*2))
              NodeHolder(2,2) = Mesh % Nodes % y(NodeIndexes(k*2))
              NodeHolder(2,3) = Mesh % Nodes % z(NodeIndexes(k*2))

              !Check for plume presence in this element
              CALL PlanePointIntersection((/Xs(Pl), Ys(Pl), 0.0_dp/), PlumeNormal, &
                 NodeHolder(1,:), NodeHolder(2,:),PlInt, found_intersection)

              IF(found_intersection .AND. &
                 ( (PlInt(1) > NodeHolder(1,1)) .NEQV. (PlInt(1) > NodeHolder(2,1)) )) THEN

                 IHavePlume(Pl) = .TRUE.
                 IF(.NOT. ANY(PlumePoints(Pl,:,3) == PlInt(3))) THEN
                    county = county + 1
                    IF(county > ExtrudedLevels) &
                       CALL Fatal(SolverName, "Found too many intersections...")

                    PlumePoints(Pl, county, :) = PlInt(:)
                 END IF
              END IF

           END DO
        END DO
     END DO

     DEALLOCATE(NodeIndexes)

     !Based on the xy position and the frontal geometry, find the plume base (z)

     !MPI comms
     IF(Parallel) THEN

        ! - which part has which plumes
        CALL MPI_ALLGATHER(IHavePlume, PlCount, MPI_LOGICAL, WhoHasPlume, &
           PlCount, MPI_LOGICAL, MPI_COMM_WORLD, ierr)

        DO i=1, ParEnv % PEs
           DO j=1,PlCount
              IF(.NOT. WhoHasPlume((i-1)*PlCount + j)) CYCLE

              IF(PlumeOwner(j) /= -1) THEN
                 WRITE(Message,'(a,i0,a,i0,i0)') "Plume ",j," has multiple owners: ",PlumeOwner(j), i
                 CALL Warn(SolverName, Message)
              END IF
              PlumeOwner(j) = i
           END DO
        END DO

        !Potential for plume to have multiple owners - i.e. plume def
        !is split over multiple partitions. If this is the case, negotiate here,
        !before passing between all parts
        DO Pl=1,PlCount
          !nan check
          IF(((ParEnv % MyPE + 1) == PlumeOwner(Pl)) .AND. &
             ANY(PlumePoints(Pl,:,1) /= PlumePoints(Pl,:,1)) ) THEN

            county = COUNT(PlumePoints(Pl,:,1) /= PlumePoints(Pl,:,1))
            ALLOCATE(idx(county),PointStore(county))
            county = 0
            DO i=1,ExtrudedLevels
              IF(PlumePoints(Pl,i,1) /= PlumePoints(Pl,i,1)) THEN
                county = county + 1
                idx(county) = i
              END IF
            END DO

            DO i=1,3
              CALL MPI_RECV(PointStore, county, MPI_DOUBLE, MPI_ANY_SOURCE, &
                 1000+i, MPI_COMM_WORLD, status, ierr)
              PlumePoints(Pl,idx,i) = PointStore
            END DO

            DEALLOCATE(PointStore,idx)

          ELSE IF(((ParEnv % MyPE + 1) /= PlumeOwner(Pl)) .AND. &
             ANY(PlumePoints(Pl,:,1) == PlumePoints(Pl,:,1)) ) THEN

            county = COUNT(PlumePoints(Pl,:,1) == PlumePoints(Pl,:,1))
            ALLOCATE(idx(county),PointStore(county))
            county = 0
            DO i=1,ExtrudedLevels
              IF(PlumePoints(Pl,i,1) == PlumePoints(Pl,i,1)) THEN
                county = county + 1
                idx(county) = i
              END IF
            END DO

            DO i=1,3
              PointStore = PlumePoints(Pl,idx,i)
              CALL MPI_SEND(PointStore, county, MPI_DOUBLE, (PlumeOwner(Pl)-1),&
                 1000+i, MPI_COMM_WORLD, ierr)
            END DO

            DEALLOCATE(PointStore,idx)
          END IF

        END DO

        !Pass the plume geometry
        DO Pl=1,PlCount
           IF(Debug .AND. Boss) PRINT *, 'Plume ',Pl,' in partition: ',PlumeOwner(Pl)

           !Send node x,y,z
           CALL MPI_BCast(PlumePoints(Pl,:,1), ExtrudedLevels,&
              MPI_DOUBLE, PlumeOwner(Pl)-1, MPI_COMM_WORLD, ierr)
           CALL MPI_BCast(PlumePoints(Pl,:,2), ExtrudedLevels,&
              MPI_DOUBLE, PlumeOwner(Pl)-1, MPI_COMM_WORLD, ierr)
           CALL MPI_BCast(PlumePoints(Pl,:,3), ExtrudedLevels,&
              MPI_DOUBLE, PlumeOwner(Pl)-1, MPI_COMM_WORLD, ierr)
        END DO

     END IF !parallel

     IF(Boss .AND. Debug) THEN
        PRINT *,'Plume x: ', PlumePoints(:,:,1)
        PRINT *,'Plume y: ', PlumePoints(:,:,2)
        PRINT *,'Plume z: ', PlumePoints(:,:,3)
     END IF

     !Cycle nodes on front
     DO i=1, Mesh % NumberOfNodes
        IF(MeltPerm(i) <= 0) CYCLE

        NodeHolder(1,1) = Mesh % Nodes % x(i)
        NodeHolder(1,2) = Mesh % Nodes % y(i)
        NodeHolder(1,3) = Mesh % Nodes % z(i)

        !Not submarine
        IF(NodeHolder(1,3) > SeaLevel) THEN
           CYCLE
        END IF

        !--------------------------------------------
        !Get contribution from background melt rate
        !--------------------------------------------
        IF(BMeltSwitch) THEN
          NodeElev = ElevVar % Values(ElevVar % Perm(i))
          IF(NodeElev < 0.0_dp) NodeElev = 0.0_dp !might be -1.0E-100...

          IF(BMFromFile) THEN
            !1D interp from z, melt file
            IF(season < BMSummerStart .OR. season > BMSummerStop) THEN
              BMeltRate(MeltPerm(i)) = FindPointOnLine(BMeltWinter % z, &
                 BMeltWinter % meltrate, NodeElev)
            ELSE
              BMeltRate(MeltPerm(i)) = FindPointOnLine(BMeltSummer % z, &
                 BMeltSummer % meltrate, NodeElev)
            END IF

            PRINT *,'Node: ',i,' z:', Mesh % Nodes % z(i), NodeElev, &
               ' background melt: ', BMeltRate(MeltPerm(i))
          ELSE
            !Background melt rate increases from zero at elev=0, to BMeltMax @ BMeltMME
            !Then it decreases at a fixed dm/dz, until (if) it reaches zero.
            IF(NodeElev > BMeltMME) THEN
              BMeltRate(MeltPerm(i)) = BMeltMax + BMeltDmDz*(NodeElev - BMeltMME)
              IF(BMeltRate(MeltPerm(i)) < 0.0_dp) BMeltRate(MeltPErm(i)) = 0.0_dp
            ELSE
              BMeltRate(MeltPerm(i)) = BMeltMax * (NodeElev / BMeltMME)
            END IF
          END IF
        END IF

        !----------------------------------
        !Get contributions from each conical plume
        !----------------------------------
        Melt = 0.0_dp

        DO Pl=1,PlCount

           IF(.NOT. PlActive(Pl)) CYCLE !plume isn't active
           IF(NodeHolder(1,3) < PlumePoints(Pl,1,3)) CYCLE !below plume base

           !Find plume points above and below current node
           DO j=2,ExtrudedLevels
              IF(NodeHolder(1,3) < PlumePoints(Pl,j,3)) EXIT
           END DO
           IF(j > ExtrudedLevels) CALL Fatal(SolverName,&
              "Didn't anticipate this, plume definition doesn't reach sea level?")

           !Interpolate the plume definition to current point's vertical coordinate
           prop = (NodeHolder(1,3) - PlumePoints(Pl,j-1,3)) / &
              (PlumePoints(Pl,j,3) - PlumePoints(Pl,j-1,3))

           NodeHolder(2,1) = (PlumePoints(Pl,j-1,1) * (1-prop)) + (PlumePoints(Pl,j,1) * prop)
           NodeHolder(2,2) = (PlumePoints(Pl,j-1,2) * (1-prop)) + (PlumePoints(Pl,j,2) * prop)
           NodeHolder(2,3) = (PlumePoints(Pl,j-1,3) * (1-prop)) + (PlumePoints(Pl,j,3) * prop)

           IF( (NodeHolder(1,3) - NodeHolder(2,3)) > 1.0E-10 ) THEN
              PRINT *,'points: ', NodeHolder(1,3), NodeHolder(2,3)
              CALL Fatal(SolverName, "Error in plume centrepoint calculation")
           END IF
           !How far from the plume centreline is the point?
           hor_dist = ( (NodeHolder(1,1) - NodeHolder(2,1))**2 + &
                (NodeHolder(1,2) - NodeHolder(2,2))**2) ** 0.5

           !How far from the base of the plume?
           plume_z = NodeHolder(2,3) - PlumePoints(Pl,1,3)

           IF(PlFromFile) THEN

             xx => ConicalPlumes(Pl) % z
             yy => ConicalPlumes(Pl) % meltrate
             d0_meltrate = FindPointOnLine(xx, yy, plume_z)

             yy => ConicalPlumes(Pl) % width
             plume_width = FindPointOnLine(xx, yy, plume_z)
           ELSE

             !Width at base, plus increase through z
             plume_width = W0(Pl) + plume_z * DwDz(Pl)

             IF(plume_z < MME(Pl)) THEN
               !Below the plume peak melt elev, varies linearly from zero to MMR
               prop = plume_z / MME(Pl)
               d0_meltrate = prop * MMR(Pl)
             ELSE
               !Above peak melt elev, drops off at rate DmDz
               d0_meltrate = MMR(Pl) +  ((plume_z - MME(Pl)) * DmDz(Pl))
               IF(d0_meltrate < 0.0_dp) d0_meltrate = 0.0_dp
             END IF

           END IF

           Melt = Melt + d0_meltrate * EXP(-(hor_dist/plume_width)**2.0)
         END DO

         PMeltRate(MeltPerm(i)) = Melt
     END DO

     !For every node, check for Plume melt exceeding Background melt, and turn off the latter
     DO i=1,SIZE(PMeltRate)
       IF(PMeltRate(i) > BMeltRate(i)) THEN
         BMeltRate(i) = 0.0_dp
       ELSE
         PMeltRate(i) = 0.0_dp
       END IF
     END DO
   END IF !Not Calving

   IF(Calving .AND. TotalPlCount > 0) THEN
     !PRINT *, 'P7',ParEnv % myPE
     !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)
   !Get model outputs into BMelt variable
     ALLOCATE(SearchIndex(2), FPOLZ(OutputSize), FPOLMR(OutputSize))
     DO i=1, Mesh % NumberOfNodes
       IF(MeltPerm(i) <= 0) CYCLE
       IF(Mesh % Nodes % z(i) .GE. 0.0) THEN
          PMeltRate(MeltPerm(i)) = 0.0
          CYCLE
       END IF
       IF(AxisIndex==1) THEN
         Node = Mesh % Nodes % x(i)
       ELSE
         Node = Mesh % Nodes % y(i)
       END IF   
       SearchIndex = 0.0_dp
       PlDist = 0.0_dp
       MinDist = 100000.0
       NoPlume = .FALSE.
  
       !Find which two plumes between and calculate 1D distance. Bit of a
       !simplification, but don't have non-major-axis co-ordinate here
       !TODO Put both x and y co-ordinates in PlCoordArray?
       DO j=1,TotalPlCount
         IF(Node >= PlCoordArray(j,1)) CYCLE
         !If point less than 1st plume, i.e. off one end of array
         IF(j==1) THEN
           SearchIndex(1) = PlCoordArray(j,2)
           PlDist(1) = 0.0 
           SearchIndex(2) = PlCoordArray(j,2)
           PlDist(2) = ABS(PlCoordArray(j,1) - Node)
           NoPlume = .TRUE.
         ELSE
           SearchIndex(1) = PlCoordArray(j,2)
           PlDist(1) = ABS(PlCoordArray(j,1) - PlCoordArray(j-1,1))
           SearchIndex(2) = PlCoordArray(j-1,2)
           PlDist(2) = ABS(PlCoordArray(j-1,1) - Node)
         END IF
         EXIT
       END DO
       
       IF(ALL(SearchIndex == 0.0)) THEN
         !If cycles through whole array without finding plumes to be between;
         !i.e. point greater than all plumes, so off other end of array
         SearchIndex(1) = PlCoordArray(TotalPlCount,2)
         PlDist(1) = 0.0
         SearchIndex(2) = PlCoordArray(TotalPlCount,2)
         PlDist(2) = ABS(PlCoordArray(TotalPlCount,1) - Node)
         NoPlume = .TRUE.
       END IF

       Node = ABS(Mesh % Nodes % z(i))
       
       FPOLZ(:) = ABS(PlZArray(:,SearchIndex(1)))
       FPOLMR(:) = PlMRArray(:,SearchIndex(1))
       ZPointer => FPOLZ
       MRPointer => FPOLMR
       Plume1MR = FindPointOnLine(ZPointer, MRPointer, Node)
       FPOLZ(:) = ABS(PlZArray(:,SearchIndex(2)))
       FPOLMR(:) = PlMRArray(:,SearchIndex(2))
       ZPointer => FPOLZ
       MRPointer => FPOLMR
       Plume2MR = FindPointOnLine(ZPointer, MRPointer, Node)

       !Then linearly interpolate between the two resulting melt values at that
       !depth and then add Gaussian decay, if relevant (i.e. if off one end of
       !plume list)
       IF(PlDist(1) == 0.0) THEN
         PlProp = 0.0
       ELSE      
         PlProp = (PlDist(2)) / (PlDist(1))
       END IF
       PMeltRate(MeltPerm(i)) = (PlProp * Plume1MR + ((1-PlProp) * Plume2MR))
       IF(NoPlume) THEN
         IF(PlDist(2) .NE. 0.0) THEN
           PMeltRate(MeltPerm(i)) = PMeltrate(MeltPerm(i)) * EXP(-(PlDist(2))**2.0)
         END IF
       END IF
       IF(PMeltRate(MeltPerm(i)) < 0.0) PMeltRate(MeltPerm(i)) = 0.0

       !May also want to do some sort of smoothing here
     END DO
     DEALLOCATE(SearchIndex, FPOLZ, FPOLMR)
     NULLIFY(ZPointer, MRPointer)
   END IF
   !PRINT *, 'P8',ParEnv % myPE
   !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

   IF(RemoveToe) THEN
     !Identify columns of front nodes, and add artificial 'toe' melting to
     !ensure that melt doesn't decrease downwards.
     DO WHILE(.TRUE.) !cycle columns
       !Find a new column
       col = -1
       DO j=1, Mesh % NumberOfNodes
         IF(MeltPerm(j) <= 0) CYCLE
         IF(FNColumns(j) /= -1) THEN !not done
           col = FNColumns(j)
           EXIT
         END IF
       END DO
       IF(col == -1) EXIT !All done

       !Gather front nodes in specified column
       WorkNodes % NumberOfNodes = COUNT((MeltPerm > 0) .AND. (FNColumns == col))
       n = WorkNodes % NumberOfNodes
       ALLOCATE(WorkNodes % z(n), ColumnPerm(n))

       counter = 1
       DO j=1, Mesh % NumberOfNodes
         IF(MeltPerm(j) <= 0) CYCLE
         IF(FNColumns(j) == col) THEN
           WorkNodes % z(counter) = Mesh % Nodes % z(j)
           ColumnPerm(counter) = j
           counter = counter + 1
         END IF
       END DO

       !Order by ascending WorkNodes % z
       ALLOCATE(OrderPerm(n))
       OrderPerm = [(i,i=1,n)]
       CALL SortD( n, WorkNodes % z, OrderPerm )

       !Go down through column, assuring melt rate doesn't decrease
       DO i=n-1,1,-1

         aboveIdx = MeltPerm(ColumnPerm(OrderPerm(i+1)))
         meIdx = MeltPerm(ColumnPerm(OrderPerm(i)))
         IF(WorkNodes % z(i) .GE. 0.0) THEN
           ToeMeltRate(meIdx) = 0.0
           CYCLE
         END IF

         aboveMelt = PMeltRate(aboveIdx) + BMeltRate(aboveIdx) + ToeMeltRate(aboveIdx)
         meMelt = PMeltRate(meIdx) + BMeltRate(meIdx) + ToeMeltRate(meIdx)
         !IF(Debug) PRINT *,ParEnv % MyPE,' Debug, me, above melt: ',meMelt, aboveMelt
         IF(aboveMelt > meMelt) THEN
           ToeMeltRate(meIdx) = aboveMelt - meMelt
         END IF
       END DO

       FNColumns(ColumnPerm) = -1 !mark column nodes as done already

       DEALLOCATE(WorkNodes % z, ColumnPerm, OrderPerm)

     END DO
     DEALLOCATE(FNColumns)
   END IF
   !PRINT *, 'P9',ParEnv % myPE
   !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

   !Possibility to scale plume and background melt contributions to observed average
   IF(AverageMelt .OR. OutputStats) THEN

      IF(AverageMelt) PRINT *,'Plume, scaling melt rates to observed averages'

      IF ( CurrentCoordinateSystem() /= Cartesian ) &
           CALL Fatal(SolverName, "Only cartesian coordinate system supported.")

      n = Mesh % MaxElementNodes
      ALLOCATE(ElementNodes % x(n),&
           ElementNodes % y(n),&
           ElementNodes % z(n))

      !Integrate melt rate and area
      TotalPMelt = 0.0_dp
      TotalBMelt = 0.0_dp
      TotalArea = 0.0_dp

      IF(RemoveToe) THEN
        TotalToeMelt = 0.0_dp
      END IF

      Active = GetNOFActive()
      DO i=1,Active

         Element => GetActiveElement(i)
         ElemBMelt = 0.0_dp
         ElemPMelt = 0.0_dp
         IF(RemoveToe) ElemToeMelt = 0.0_dp

         n = Element % TYPE % NumberOfNodes
         NodeIndexes => Element % NodeIndexes

         ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

         IF(ALL(ElementNodes % z(1:n) > SeaLevel)) CYCLE !Entirely above sea level, disregard

         IntegStuff = GaussPoints( Element )

         DO j=1,IntegStuff % n

            U = IntegStuff % u(j)
            V = IntegStuff % v(j)
            W = IntegStuff % w(j)

            !------------------------------------------------------------------------------
            !        Basis function values at the integration point
            !------------------------------------------------------------------------------
            stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
                 Basis )

            !assume cartesian here
            s = SqrtElementMetric * IntegStuff % s(j)

            ElemBMelt = ElemBMelt + s * SUM(BMeltRate(MeltPerm(NodeIndexes)) * Basis(1:n))
            ElemPMelt = ElemPMelt + s * SUM(PMeltRate(MeltPerm(NodeIndexes)) * Basis(1:n))

            IF(RemoveToe) THEN
              ElemToeMelt = ElemToeMelt + &
                   s * SUM(ToeMeltRate(MeltPerm(NodeIndexes)) * Basis(1:n))
            END IF

            TotalArea = TotalArea + s
         END DO

         TotalPMelt = TotalPMelt + ElemPMelt
         TotalBMelt = TotalBMelt + ElemBMelt

         IF(RemoveToe) TotalToeMelt = TotalToeMelt + ElemToeMelt
      END DO
      DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z)

      IF(Debug) THEN
         PRINT *, 'Plume total submarine area: ',TotalArea
         PRINT *, 'Plume total background melt: ',TotalBMelt
         PRINT *, 'Plume total plume melt: ',TotalPMelt
         IF(RemoveToe) THEN
           PRINT *, 'Plume total toe melt: ', TotalToeMelt
         END IF
      END IF

      CALL MPI_AllReduce(MPI_IN_PLACE, TotalArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_AllReduce(MPI_IN_PLACE, TotalPMelt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_AllReduce(MPI_IN_PLACE, TotalBMelt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      IF(RemoveToe) THEN
        CALL MPI_AllReduce(MPI_IN_PLACE, TotalToeMelt, 1, MPI_DOUBLE, &
             MPI_SUM, MPI_COMM_WORLD, ierr)
      END IF

      BMelt_Average = TotalBMelt / TotalArea
      PMelt_Average = TotalPMelt / TotalArea

      IF(AverageMelt) THEN
        IF(RemoveToe) CALL Fatal(SolverName, &
             "Calving Toe Removal doesn't play with Scale to Average Melt")

        scale = Target_BMelt_Average / BMelt_Average
        BMeltRate = BMeltRate * scale
        !Also scale total value for output
        TotalBMelt = TotalBMelt * scale

        IF(Debug) PRINT *,'Plume, Scaling background melt by factor: ', scale

        scale = Target_PMelt_Average / PMelt_Average
        PMeltRate = PMeltRate * scale
        !Also scale total value for output
        TotalPMelt = TotalPMelt * scale

        IF(Debug) PRINT *,'Plume, Scaling plume melt by factor: ', scale
      END IF
    END IF
    !PRINT *, 'P10',ParEnv % myPE
    !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    IF(OutputStats) THEN
      IF(ParEnv % MyPE == 0) THEN

        IF(.NOT. Visited) THEN
          OPEN( UNIT=OutFileUnit, File=OutfileName, STATUS='UNKNOWN')

          IF(RemoveToe) THEN
            WRITE(OutFileUnit, '(A)') "Timestep, Time, No. Plumes, Submerged Area, &
                 &Total Background Melt (m^3), Total Plume Melt (m^3), &
                 &Total Melt (m^3), Total Melt (actual)(m^3), &
                 &Forced Toe Calving (m^3)"
          ELSE
            WRITE(OutFileUnit, '(A)') "Timestep, Time, No. Plumes, Submerged Area, &
                 &Total Background Melt (m^3), Total Plume Melt (m^3), &
                 &Total Melt (m^3)"
          END IF
        ELSE
          OPEN( UNIT=OutFileUnit, File=OutfileName, STATUS='UNKNOWN', ACCESS='APPEND' )
        END IF

        IF(RemoveToe) THEN
          WRITE(OutFileUnit, &
               '(I0,ES20.11,i5,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11,ES20.11)') &
               GetTimestep(), GetTime(), COUNT(PlActive),TotalArea, &
               TotalBMelt, TotalPMelt, TotalBMelt+TotalPMelt+TotalToeMelt,&
               TotalBMelt+TotalPMelt,TotalToeMelt

        ELSE
          WRITE(OutFileUnit, '(I0,ES20.11,i5,ES20.11,ES20.11,ES20.11,ES20.11)') &
               GetTimestep(), GetTime(), COUNT(PlActive),TotalArea, &
               TotalBMelt, TotalPMelt, TotalBMelt+TotalPMelt
        ENDIF
        CLOSE( OutFileUnit )
      END IF
    END IF
    !PRINT *, 'P11',ParEnv % myPE
    !CALL MPI_BARRIER(ELMER_COMM_WORLD, ierr)

    !Total melt rate is background rate + plume
    MeltRate = PMeltRate + BMeltRate

    IF(RemoveToe) THEN
      MeltRate = MeltRate + ToeMeltRate
      ToeCalveVar % Values = ToeMeltRate
    END IF
    
    NULLIFY(WorkVar, WorkSolver, HydroMesh, Edge, ZOutput, MROutput, WorkVar2,&
           WorkVar3)
    DEALLOCATE(BMeltRate, PMeltRate, HydroGLNodes)
    IF(Calving) DEALLOCATE(PlActive)
    IF(Calving .AND. PlCount > 0) THEN
      DEALLOCATE(Zi, Xi, Sa, Ta)!, ZiInput)
      DO i=1, SIZE(HydroPlume)
        NULLIFY(HydroPlume(i) % z, HydroPlume(i) % meltrate)
      END DO
      DEALLOCATE(HydroPlume)
    END IF
    IF(Calving .AND. TotalPlCount > 0) THEN
      DEALLOCATE(NearestFrontNodes, PlPos, PlFinalQ, PlAxis,&
                PlCoordArray, PlZArray, PlMRArray)
    END IF
    IF(RemoveToe) DEALLOCATE(ToeMeltRate)
    Visited = .TRUE.

    CONTAINS

      FUNCTION FindPointOnLine(xx,yy,x) RESULT(y)
        REAL(KIND=dp) :: x,y,prop
        REAL(KIND=dp), POINTER :: xx(:),yy(:)
        INTEGER :: i, indexx

        IF(SIZE(xx)/=SIZE(yy)) CALL Fatal(SolverName, "Input file has different size columns!")

        indexx = 0

        DO i=1,SIZE(xx)
          IF(x >= xx(i)) CYCLE
          indexx = i
          EXIT
        END DO

        IF(indexx <= 1) THEN !Either below first, or above last, 0 melt
          y = 0.0_dp
        ELSE
          prop = (x - xx(indexx-1)) / (xx(indexx) - xx(indexx-1))
          y = (prop * yy(indexx)) + ((1-prop) * yy(indexx-1))
        END IF

      END FUNCTION FindPointOnLine

      SUBROUTINE PlanePointIntersection ( pp, pnorm, p1, p2, p_intersect, found_intersection )
        !This is defined in the CalvingGeometry module, but, for some reason
        !I have yet to fathom, Elmer kept throwing a paddy about not being able
        !to find it, so I've copied it in here for my own sanity.
        !Get the intersection point between a line and plane in 3D
        ! Plane defined by point "pp" and norm "pnorm", line defined by points "p1" and "p2"
        ! Intersection returned in p_intersect
        !found_intersection = .FALSE. if they happen to be parallel

        REAL(KIND=dp) :: pp(3), pnorm(3), p1(3), p2(3), p_intersect(3)
        LOGICAL :: found_intersection
        !----------------------------
        REAL(KIND=dp) :: pl(3), dist

        pl = p2 - p1

        IF(ABS(DOT_PRODUCT(pl,pnorm)) < EPSILON(1.0_dp)) THEN
          !Line and plane are parallel...
          found_intersection = .FALSE.
          RETURN
        END IF

        dist = DOT_PRODUCT((pp - p1), pnorm) / DOT_PRODUCT(pl,pnorm)

        p_intersect = p1 + dist*pl
        found_intersection = .TRUE.

      END SUBROUTINE PlanePointIntersection
    END SUBROUTINE Plume
!-------------------------------------------------------------------------------
!This routine is an adapatation of the IcePlume package implemented in MITgcm
!by Tom Cowton, which is itself based on work by Donald Slater.
!Solves for a linear plume at all nodes along the calving front
!Works on a m-s standard, rather than the m-a standard used generally in Elmer,
!so time-dependent inputs (discharge) and outputs (melt rate) need to be
!appropriately handled to ensure interfaces correctly with rest of Elmer
!Uses the ODEPack freely-available ODE library, which relies on implicit
!variables, hence lack of IMPLICIT NONE statements here.
!Assumes Zi starts at 0 with depths in minus numbers, last entry being maximum
!depth at which data defined. Xi, Ta and Sa then assumed to be ordered in same
!manner (e.g. last entry in Ta is assumed to be Ta at max depth in Zi). All four
!arrays should therefore be the same size.
  SUBROUTINE PlumeSolver (Depth, Front, TempA, SalA, Discharge, DepthOutput,&
                          MeltOutput, MeshRes, PlDepth)
    USE Types

    EXTERNAL SheetPlume, JEX
       
    REAL(KIND=dp), ALLOCATABLE :: Depth(:), Front(:), TempA(:), SalA(:)
    REAL(KIND=dp) :: Discharge, MeshRes, PlDepth
    REAL(KIND=dp), POINTER :: DepthOutput(:), MeltOutput(:)

    !-----------------------------------------------------------------------

    INTEGER :: i, N, K, Counter, Counter2
    REAL(KIND=dp) :: RHO, Temperature, Salinity, Depth2, Tambient, Sambient,&
                     Rho_Plume, Rho_Ambient
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Z(:), X(:), Temp(:), Sal(:), MR(:),&
                                          MRZ(:)
    REAL(KIND=dp), POINTER :: ZPointer(:), TPointer(:), SPointer(:),&
                              MRPointer(:), MRZPointer(:)
    REAL(KIND=dp), ALLOCATABLE :: ZProf(:), RProf(:), WProf(:), AProf(:),&
                                  MIntProf(:)
    REAL(KIND=dp), ALLOCATABLE :: TProf(:), SProf(:), ZProfAbs(:)
    REAL(KIND=dp) :: TCommon(1000), SCommon(1000), ZCommon(1000), MDCommon(1000,2)
    COMMON TCommon, SCommon, ZCommon, MDCommon, Counter, Counter2

    !For ODEPACK
    INTEGER :: IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,&
               MF, NEQ
    REAL(KIND=dp) :: Y(7), T, TOUT, RTOL, ATOL, RWORK(116)

    !Y is the I/O vector for DLSODE
    !Y(1) = plume thickness/radius
    !Y(2) = plume velocity
    !Y(3) = plume temperature
    !Y(4) = plume salinity
    !Y(5) = plume area
    !Y(6) = area-integrated melt

    NEQ = 6
    LRW = 116
    LIW = 116
    ITOL = 1
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    MF = 10
    IWORK(7) = 2 !Limits number of times repeat error messages are printed
    RTOL = 1.0D-5
    ATOL = 1.0D-5

    !Actual initial velocity and radius values not important, provided they
    !combine to give the desired initial discharge, hence setting velocity to 1
    !and radius to actual discharge
    !No point plumes starting below point where have fjord data (melting below
    !this point, if relevant, dealt with by Toe Calving routine), so initial
    !conditions for temperature and salinity set to those at max data depth
    !Plume initial conditions assumed to be basically fresh, cold water
    ALLOCATE(Z(SIZE(Depth)), Sal(SIZE(SalA)), Temp(SIZE(TempA)), MR(1000),&
            MRZ(1000))
    Z = Depth
    Sal = SalA
    Temp = TempA
    ZPointer => Z
    SPointer => Sal
    TPointer => Temp
    N = SIZE(Z)
    Y(1) = 1.0
    Y(2) = Discharge !0.1
    Y(3) = 1.0D-3
    Y(4) = 1.0D-3
    Y(5) = 0.0
    Y(6) = 0.0
    Y(7) = MeshRes

    !PlDepth = -83.0

    !Set up output profiles
    ALLOCATE(ZProf(N), ZProfAbs(N), RProf(N), WProf(N), TProf(N), SProf(N),&
            AProf(N), MIntProf(N))

    DO i=1, N
      ZProf(i) = Z(i)
      ZProfAbs(i) = ABS(Z(i))
      RProf(i) = 0.0
      WProf(i) = 0.0
      TProf(i) = 0.0
      SProf(i) = 0.0
      AProf(i) = 0.0
      MIntProf(i) = 0.0
    END DO

    ZCommon(1:N) = ZProf(1:N)
    ZCommon(N+1:SIZE(ZCommon)) = -9999.9
    MDCommon(:,:) = 0.0
    MRPointer => MR
    MRZPointer => MRZ
    !MDCommon(N+1:SIZE(MDCommon)) = -9999.9

    !Start at bottom of data range
    T = ZProf(N)
    !Next point to move up to
    TOUT = ZProf(N-1)

    !Set initial conditions
    RProf(N) = Y(1)
    WProf(N) = Y(2)
    TProf(N) = Y(3)
    SProf(N) = Y(4)
    AProf(N) = Y(5)
    MIntProf(N) = Y(6)
    TCommon(1:N) = Temp(1:N)
    SCommon(1:N) = Sal(1:N)
    TCommon(N+1:SIZE(TCommon)) = -9999.9
    SCommon(N+1:SIZE(SCommon)) = -9999.9

    K = N-1
    Counter2 = 1

    !Start at bottom and work up through supplied profiles
    DO IOUT=K,1,-1
      Counter = IOUT+1
      !Cycle through profiles until reach depth plume should start at
      IF(PlDepth > T) THEN
        T=TOUT !Define present depth
        TOUT = ZProf(IOUT-1) !Define next depth 
        RProf(IOUT) = 0.0
        WProf(IOUT) = 0.0
        TProf(IOUT) = 0.0
        SProf(IOUT) = 0.0
        AProf(IOUT) = 0.0
        MIntProf(IOUT) = 0.0
        RProf(IOUT-1) = Y(1)
        WProf(IOUT-1) = Y(2)
        TProf(IOUT-1) = Y(3)
        SProf(IOUT-1) = Y(4)
        AProf(IOUT-1) = Y(5)
        MIntProf(IOUT-1) = Y(6)
        CYCLE
      END IF
        
      !Checks for plume reaching neutral buoyancy in previous layer
      IF(ISTATE > -1) THEN
        Y(7) = MeshRes
        CALL DLSODE(SheetPlume, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,&
             ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)

        !Test for neutral buoyancy now. If solver returns ISTATE = -1, it's not
        !been able to meet required tolerances, which is usually because plume
        !has reached neutral buoyancy, so isn't rising farther. So, no point
        !continuing plume model, so end call.
        !Want to catch plume at point of neutral buoyancy, so need to compare
        !ambient and plume density too. If the plume density is >= the ambient
        !density, we also set ISTATE = -1 to end the call to the plume model.

        Temperature = Y(3)
        Salinity = Y(4)
        Depth2 = T

        !Calculate density
        Rho_Plume = RHO(Temperature, Salinity, Depth2)

        IF(IOUT == 1) THEN
          TAmbient = TProf(1)
          SAmbient = SProf(1)
        ELSE
          TAmbient = 0.5*(Temp(IOUT-1)+Temp(IOUT))
          SAmbient = 0.5*(Sal(IOUT-1)+Sal(IOUT))
        END IF

        !Calculate ambient density
        Rho_Ambient = RHO(TAmbient, SAmbient, Depth2)

        IF(Rho_Plume > Rho_Ambient) ISTATE = -1

        !If ISTATE now < 0, we've hit neutral buoyancy and can zero plume area
        !and velocity.
        IF(ISTATE < 0) THEN
          Y(1) = 0.0_dp
          Y(2) = 0.0_dp
        ELSE
          !If plume not at neutral buoyancy, move up to next depth value and 
          !call plume model again, unless we're at the surface....

          IF(IOUT .NE. 1) THEN

            T=TOUT !Define present depth
            TOUT = ZProf(IOUT-1) !Define next depth 
          END IF
        END IF
      ELSE !For when plume has reached neutral buoyancy
 
        !No plume values
        Y(1) = 0.0
        Y(2) = 0.0
        Y(3) = 0.0
        Y(4) = 0.0
        Y(5) = 0.0
        Y(6) = 0.0
      END IF !ISTATE > -1

      !Save results for that step
      RProf(IOUT) = Y(1)
      WProf(IOUT) = Y(2)
      TProf(IOUT) = Y(3)
      SProf(IOUT) = Y(4)
      AProf(IOUT) = Y(5)

      !This works. For some reason, Y(6) gives stupid melt rates (50 m/d+), but
      !mdot in the SheetPlume subroutine is correct (according to Donald Slater,
      !who compared it against his model), so this is to extract it and get it
      !into the output
      MR = MDCommon(:,2)
      MRZ = MDCommon(:,1)
      MIntProf(IOUT) = FindPointOnLine(MRZPointer, MRPointer, ZProf(IOUT+1))*86400
    
    END DO !IOUT

    !To get from m/d to m/a to fit with rest of Elmer
    DepthOutput = ZProf
    MeltOutput = MIntProf*365.25

    DEALLOCATE(Z, Sal, Temp)
    DEALLOCATE(ZProf, ZProfAbs, RProf, WProf, TProf, SProf, AProf, MIntProf)
    NULLIFY(ZPointer, SPointer, TPointer, MRZPointer, MRPointer)

    CONTAINS
      FUNCTION FindPointOnLine(xx,yy,x) RESULT(y)
        REAL(KIND=dp) :: x,y,prop
        REAL(KIND=dp), POINTER :: xx(:),yy(:)
        INTEGER :: i, indexx

        IF(SIZE(xx)/=SIZE(yy)) CALL Fatal("FindPointOnLine", "Input file has different size columns!")

        indexx = 0

        DO i=1,SIZE(xx)
          IF(x >= xx(i)) CYCLE
          indexx = i
          EXIT
        END DO

        IF(indexx <= 1) THEN !Either below first, or above last, 0 melt
          y = yy(1) !0.0_dp
        ELSE
          prop = (x - xx(indexx-1)) / (xx(indexx) - xx(indexx-1))
          y = (prop * yy(indexx)) + ((1-prop) * yy(indexx-1))
        END IF

      END FUNCTION FindPointOnLine


  END SUBROUTINE PlumeSolver
!-------------------------------------------------------------------------------

  SUBROUTINE SheetPlume(NEQ, T, Y, YDOT)
    USE Types

    INTEGER :: NEQ, Counter, Counter2
    REAL(KIND=dp) :: T, Y(7), YDOT(6), TAmbient, SAmbient, Rho_0, Rho_1,&
                     Tin, mdot, Sb, Tb, a, b, c, RHO
    REAL(KIND=dp) :: E_0, icetemp, Rho_ref, g, c_w, c_i, L, lambda1, lambda2,&
                     lambda3, GamT, GamS, Cd
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Z(:), Temp(:), S(:)
    REAL(KIND=dp), POINTER :: ZPointer(:), TPointer(:), SPointer(:)
    REAL(KIND=dp) :: TCommon(1000), SCommon(1000), ZCommon(1000), MDCommon(1000,2)
    COMMON TCommon, SCommon, ZCommon, MDCommon, Counter, Counter2
 
    !Y(1) = r
    !Y(2) = w
    !Y(3) = T
    !Y(4) = S

    !Interpolate based on Ta and Sa
    j=0
    DO i=1, 1000
      IF(ZCommon(i) .NE. -9999.9) j=j+1
    END DO
    
    ALLOCATE(Z(j), Temp(j), S(j))
    Tin = T
    Z(:) = ZCommon(1:j)
    Temp(:) = TCommon(1:j)
    S(:) = SCommon(1:j)

    ZPointer => Z
    TPointer => Temp
    SPointer => S

    IF(T .LE. ZCommon(j)) THEN
      TAmbient = TCommon(j)
      SAmbient = SCommon(j)
    ELSE IF (T .GE. MAXVAL(ZCommon)) THEN 
      TAmbient = TCommon(MAXLOC(Z,1))
      SAmbient = SCommon(MAXLOC(Z,1))
    ELSE
      TAmbient = FindPointOnLine(ZPointer, TPointer, T)
      SAmbient = FindPointOnLine(ZPointer, SPointer, T)
    END IF

    Rho_1 = RHO(Y(3), Y(4), Z(Counter)) !Plume density
    Rho_0 = RHO(TAmbient, SAmbient, Z(Counter)) !Ambient density

    !Constants
    Rho_ref = 1020.D0
    g       = 9.81D0
    c_w     = 3994.D0
    c_i     = 2009.D0
    L       = 334000D0
    lambda1 = -0.0573D0
    lambda2 = 0.0832D0
    lambda3 = 0.000761D0
    GamT    = 0.011
    GamS    = 0.00031
    Cd      = 0.02  !0.0025 10*higher; see Ezhova et al. (2018)
    icetemp = 0.0
    E_0     = 0.1D0

    !Equations for Sb, Tb and mdot
    a = lambda1*(GamT*c_w-GamS*c_i)
    b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*T+iceTemp-(L/c_i))&
        -GamT*c_w*(Y(3)-lambda2-lambda3*T)
    c = GamS*Y(4)*(c_i*(lambda2+lambda3*T-iceTemp)+L)

    Sb = (1.0/(2.0*a))*(-b-((b**2.0-4.0*a*c)**0.5))
    Tb = lambda1*Sb+lambda2+lambda3*T
    mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb

    !Differential Equations
    !Plume Thickness
    YDOT(1) = 2*E_0+Cd-(g*Y(1)/(Y(2)**2))*(Rho_0-Rho_1)/Rho_ref+2*mdot/Y(2)
 
    !Plume Vertical Velocity
    YDOT(2) = -(Y(2)/Y(1))*(E_0+Cd+mdot/Y(2))+(g/Y(2))*(Rho_0-Rho_1)/Rho_ref

    !Plume Temperature
    YDOT(3) = E_0*TAmbient/Y(1)-(Y(3)/Y(1))*(E_0+mdot/Y(2))+(mdot/(Y(1)*Y(2)))&
              *(Tb-(L/c_w)-(c_i/c_w)*(Tb-iceTemp))

    !Plume Salinity
    YDOT(4) = E_0*SAmbient/Y(1)-(Y(4)/Y(1))*(E_0+mdot/Y(2))

    !Along-integrated Melt Rate and Contact Area
    YDOT(5) = Y(7)
    YDOT(6) = Y(7)*mdot
 
    MDCommon(Counter2,1) = T 
    MDCommon(Counter2,2) = mdot 

    T = Tin
    Counter2 = Counter2+1

    NULLIFY(ZPointer, TPointer, SPointer)
    DEALLOCATE(Z, Temp, S)
    CONTAINS

      FUNCTION FindPointOnLine(xx,yy,x) RESULT(y)
        REAL(KIND=dp) :: x,y,prop
        REAL(KIND=dp), POINTER :: xx(:),yy(:)
        INTEGER :: i, indexx

        IF(SIZE(xx)/=SIZE(yy)) CALL Fatal("FindPointOnLine", "Input file has different size columns!")

        indexx = 0

        DO i=1,SIZE(xx)
          IF(x <= xx(i)) CYCLE !Other way round because depths negative
          indexx = i
          EXIT
        END DO

        IF(indexx <= 1) THEN !Either below first, or above last, 0 melt
          y = yy(1) !0.0_dp
        ELSE
          prop = (x - xx(indexx-1)) / (xx(indexx) - xx(indexx-1))
          y = (prop * yy(indexx)) + ((1-prop) * yy(indexx-1))
        END IF

      END FUNCTION FindPointOnLine

  END SUBROUTINE SheetPlume

  SUBROUTINE JEX()
    RETURN
  END SUBROUTINE JEX

  DOUBLE PRECISION FUNCTION RHO(Temp, S, Z)

    !Equation of state (UNESCO 1983)

    !T = temperature (degrees C)
    !S = salinity (PSU)
    !Z = depth (m)
    USE Types

    REAL(KIND=dp) :: Temp, S, Z

    REAL(KIND=dp) :: Rho_0, g, P, kw, Aw, Bw, k0, bulk_modulus, A, B, C, Rho_w,&
                   Rho_zero

    PARAMETER(rho_0=1027)
    PARAMETER(g=9.81)

    P = Rho_0*g*ABS(z)*1.0E-5

    kw = 19652.21+148.4206*Temp-2.327105*Temp**2+1.360477E-2*(Temp**3)-&
         5.155288E-5*(Temp**4)
    Aw = 3.239908+1.43713E-3*Temp+1.16092E-4*Temp**2-5.77905E-7*Temp**3
    Bw = 8.50935E-5-6.12293E-6*Temp+5.2787E-8*(Temp**2)
    k0 = kw+(54.6746-0.603459*Temp+1.09987E-2*(Temp**2)-6.1670E-5*(Temp**3))*S&
         +(7.944E-2+1.6483E-2*Temp-5.3009E-4*(Temp**2))*(S**1.5)
    A = Aw+(2.2838E-3-1.0981E-5*Temp-1.6078E-6*(Temp**2))*S+1.91075E-4*(S**1.5)
    B = Bw+(-9.9348E-7+2.0816E-8*Temp+9.1697E-10*Temp**2)*S
    bulk_modulus = k0+A*P+B*P**2

    A = 8.24493E-1-4.0899E-3*Temp+7.6438E-5*Temp**2-8.2467E-7*Temp**3&
        +5.3875E-9*Temp**4
    B = -5.72466E-3+1.0227E-4*Temp-1.6546E-6*Temp**2
    C = 4.8314E-4
    Rho_w = 999.842594+6.793952E-2*Temp-9.095290E-3*Temp**2+1.001685E-4*&
            Temp**3-1.120083E-6*Temp**4+6.536336E-9*Temp**5
    Rho_zero = Rho_w+A*S+B*(S**1.5)+C*(S**2)

    RHO = Rho_zero/(1-(P/bulk_modulus))

  END FUNCTION RHO
