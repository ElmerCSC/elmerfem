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
! ******************************************************************************
! *
! *  Authors: Joe Todd
! *  Email:   jat71@cam.ac.uk
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 19.01.2016
! *
! ****************************************************************************/

 SUBROUTINE GMValid (Model, Solver, dt, TransientSimulation)
   USE Types
   USE CoordinateSystems
   USE DefUtils
   USE ElementDescription
   USE CalvingGeometry

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
   !-----------------------------------
   TYPE(Mesh_t), POINTER :: Mesh
   TYPE(Matrix_t), POINTER :: Matrix
   TYPE(Variable_t), POINTER :: Var, GroundedVar
   TYPE(ValueList_t), POINTER :: Params
   TYPE(CrevasseGroup3D_t), POINTER :: FloatGroups, CurrentGroup, DelGroup
   TYPE(Element_t), POINTER :: Element
   TYPE(Nodes_t) :: ElementNodes
   TYPE(GaussIntegrationPoints_t) :: IntegStuff

   REAL(KIND=dp) :: GMCheck, SMeltRate, WMeltRate, SStart, SStop, &
        TotalArea, TotalBMelt, ElemBMelt, s, t, season,&
        SqrtElementMetric,U,V,W,Basis(Model % MaxElementNodes)
   INTEGER :: DIM, NoNodes, i,j,n, FaceNodeCount, GroupNodeCount, county, &
        Active, ierr, k, FoundNew, AllFoundNew
   INTEGER, PARAMETER :: FileUnit = 75
   INTEGER, POINTER :: Perm(:), InvPerm(:), FrontPerm(:)=>NULL(), Neighbours(:,:), &
        NeighbourHolder(:), NoNeighbours(:), NodeIndexes(:)
   INTEGER, ALLOCATABLE :: AllGroupNodes(:), PartNodeCount(:), AllPartGroupNodes(:), &
        disps(:)
   LOGICAL :: Found, OutputStats, Visited=.FALSE., Debug, stat, Summer
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, GMaskVarName, FrontMaskName, OutfileName, mode

   Debug = .FALSE.

   SolverName = "GMValidator"
   Params => Solver % Values
   Mesh => Solver % Mesh

   DIM = CoordinateSystemDimension()
   IF(DIM /= 3) CALL Fatal(SolverName, "This solver only works in 3D!")

   !t = GetTime()
   !season = t - FLOOR(t)

   !mode = ListGetString( Params, 'Basal Melt Mode', Found, UnfoundFatal=.TRUE.)
   !mode = TRIM(mode)

   !Get output file name for stats
   !OutfileName = ListGetString(Params,"Basal Melt Stats File", OutputStats, UnfoundFatal=.TRUE.)

   !Identify nodes on the front
   FrontMaskName = "Calving Front Mask"
   CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
        .FALSE., FrontPerm, FaceNodeCount)

   !Need the matrix for finding neighbours
   Matrix => Solver % Matrix

   Var => Solver % Variable
   IF(.NOT. ASSOCIATED(Var)) CALL Fatal(SolverName, "Solver needs a variable!")
   Perm => Var % Perm
   Var % Values = 0.0_dp

   NoNodes = COUNT(Perm > 0)

   GMaskVarName = ListGetString(Params, "GroundedMask Variable", Found)
   IF(.NOT. Found) GMaskVarName = "GroundedMask"
   GroundedVar => VariableGet(Mesh % Variables, GMaskVarName, .TRUE., UnfoundFatal=.TRUE.)

   !SELECT CASE(mode)
   !CASE("seasonal")
     !SMeltRate = ListGetConstReal(Params, "Basal Melt Summer Rate", UnfoundFatal=.TRUE.)
     !WMeltRate = ListGetConstReal(Params, "Basal Melt Winter Rate", UnfoundFatal=.TRUE.)
     !SStart = ListGetConstReal(Params, "Basal Melt Summer Start", UnfoundFatal=.TRUE.)
     !SStop = ListGetConstReal(Params, "Basal Melt Summer Stop", UnfoundFatal=.TRUE.)

     !Summer = .FALSE.
     !IF(SStop > SStart) THEN
     !  IF(season > SStart .AND. season < SStop) Summer = .TRUE.
     !ELSE
     !  IF(season > SStart .OR. season < SStop) Summer = .TRUE.
     !END IF

     !IF(Summer) THEN
     !  MeltRate = SMeltRate
     !ELSE
     !  MeltRate = WMeltRate
     !END IF
   !CASE("off")
     GMCheck = 1.0_dp
   !CASE DEFAULT
     !CALL Fatal(SolverName, "Unknown basal melt mode, valid options are 'seasonal' and 'off'.")
   !END SELECT

   !Set up inverse perm for FindNodeNeighbours
   InvPerm => CreateInvPerm(Matrix % Perm) !Create inverse perm for neighbour search
   ALLOCATE(Neighbours(Mesh % NumberOfNodes, 10), NoNeighbours(Mesh % NumberOfNodes))
   Neighbours = 0

   !Find neighbours for each node on the bed
   DO i=1, Mesh % NumberOfNodes
     IF(Perm(i) <= 0) CYCLE

     NeighbourHolder => FindNodeNeighbours(i, Matrix, &
          Matrix % Perm, 1, InvPerm)

     Neighbours(i,1:SIZE(NeighbourHolder)) = NeighbourHolder
     NoNeighbours(i) = SIZE(NeighbourHolder)
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
     DO i=1, CurrentGroup % NumberOfNodes

       IF(FrontPerm(CurrentGroup % NodeNumbers(i)) > 0) THEN
         CurrentGroup % FrontConnected = .TRUE.
         EXIT
       END IF
     END DO
     CurrentGroup => CurrentGroup % Next
   END DO

   DO k=1,1000
     FoundNew = 0
     !Count and gather nodes from all valid groups
     GroupNodeCount = 0
     county = 0
     DO i=1,2
       IF(i==2) ALLOCATE(AllGroupNodes(GroupNodeCount))

       CurrentGroup => FloatGroups
       DO WHILE(ASSOCIATED(CurrentGroup))
         IF(CurrentGroup % FrontConnected) THEN

           IF(i==1) THEN
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
     DO i=2,ParEnv % PEs
       disps(i) = disps(i-1) + PartNodeCount(i-1)
     END DO

     CALL MPI_ALLGATHERV(AllGroupNodes, GroupNodeCount, MPI_INTEGER, &
        AllPartGroupNodes, PartNodeCount, disps, MPI_INTEGER, MPI_COMM_WORLD, ierr)

     !Cycle unconnected groups, looking for partition boundary in connected groups
     CurrentGroup => FloatGroups
     DO WHILE(ASSOCIATED(CurrentGroup))
       IF(.NOT. CurrentGroup % FrontConnected) THEN
         DO i=1,CurrentGroup % NumberOfNodes

           IF(ANY(Mesh % ParallelInfo % GlobalDOFs(CurrentGroup % NodeNumbers(i)) == &
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
   END DO !k

   !Cycle all connected groups, setting melt rate
   CurrentGroup => FloatGroups
   DO WHILE(ASSOCIATED(CurrentGroup))
     IF(CurrentGroup % FrontConnected) THEN
       DO i=1,CurrentGroup % NumberOfNodes
         Var % Values(Var % Perm(CurrentGroup % NodeNumbers(i))) = GMCheck
       END DO
     END IF
     CurrentGroup => CurrentGroup % Next
   END DO

   !IF(OutputStats) THEN

   !   IF ( CurrentCoordinateSystem() /= Cartesian ) &
   !        CALL Fatal(SolverName, "Only cartesian coordinate system supported.")

   !   n = Mesh % MaxElementNodes
   !   ALLOCATE(ElementNodes % x(n),&
   !        ElementNodes % y(n),&
   !        ElementNodes % z(n))

      !Integrate melt rate and area
   !   TotalBMelt = 0.0_dp
   !   TotalArea = 0.0_dp

   !   Active = GetNOFActive()
   !   DO i=1,Active

   !      Element => GetActiveElement(i)
   !      ElemBMelt = 0.0_dp

   !      n = Element % TYPE % NumberOfNodes
   !      NodeIndexes => Element % NodeIndexes

         !Not interested in area/melt calcs for grounded elements
   !      IF(ALL(GroundedVar % Values(GroundedVar % Perm(NodeIndexes)) >= 0)) CYCLE

   !      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
   !      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
   !      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))

   !      IntegStuff = GaussPoints( Element )

   !      DO j=1,IntegStuff % n

   !         U = IntegStuff % u(j)
   !         V = IntegStuff % v(j)
   !         W = IntegStuff % w(j)

            !------------------------------------------------------------------------------
            !        Basis function values at the integration point
            !------------------------------------------------------------------------------
   !         stat = ElementInfo( Element,ElementNodes,U,V,W,SqrtElementMetric, &
   !              Basis )

            !assume cartesian here
   !         s = SqrtElementMetric * IntegStuff % s(j)

            !Check here for grounded
   !         ElemBMelt = ElemBMelt + s * SUM(Var % Values(Var % Perm(NodeIndexes)) * Basis(1:n))
   !         TotalArea = TotalArea + s
   !      END DO

   !      TotalBMelt = TotalBMelt + ElemBMelt
   !   END DO
   !   DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z)

   !   IF(Debug) THEN
   !      PRINT *, 'BasalMelt3D: total submarine area: ',TotalArea
   !      PRINT *, 'BasalMelt3D: total background melt: ',TotalBMelt
   !   END IF

   !   CALL MPI_AllReduce(MPI_IN_PLACE, TotalArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
   !   CALL MPI_AllReduce(MPI_IN_PLACE, TotalBMelt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

   !   IF(ParEnv % MyPE == 0) THEN
   !     IF(.NOT. Visited) THEN
   !       Visited = .TRUE.
   !       OPEN( UNIT=FileUnit, File=OutfileName, STATUS='UNKNOWN')
   !       WRITE(FileUnit, '(A)') "Timestep, Time, Ungrounded Area, &
   !            &Total Basal Melt (m^3)"
   !     ELSE
   !       OPEN( UNIT=FileUnit, File=OutfileName, STATUS='UNKNOWN', ACCESS='APPEND' )
   !     END IF

   !     WRITE(FileUnit, '(I0,ES20.11,ES20.11,ES20.11)') &
   !          GetTimestep(), GetTime(), TotalArea, TotalBMelt

   !     CLOSE( FileUnit )
   !   END IF

      ! IF(AverageMelt) THEN
      !   scale = Target_BMelt_Average / BMelt_Average
      !   BMeltRate = BMeltRate * scale
      !   !Also scale total value for output
      !   TotalBMelt = TotalBMelt * scale

      !   IF(Debug) PRINT *,'Plume, Scaling background melt by factor: ', scale

      !   scale = Target_PMelt_Average / PMelt_Average
      !   PMeltRate = PMeltRate * scale
      !   !Also scale total value for output
      !   TotalPMelt = TotalPMelt * scale

      !   IF(Debug) PRINT *,'Plume, Scaling plume melt by factor: ', scale
      ! END IF
   ! END IF

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
 END SUBROUTINE GMValid
