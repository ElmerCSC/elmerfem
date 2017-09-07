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
!
!/******************************************************************************
! *
! *  Solver to predict calving events in 2D based on crevasse depth calving
! *  criterion (Benn et al. 2007). Description of algorithm in 
! *  Todd & Christoffersen (2014) [doi:10.5194/tc-8-2353-2014]
! *
! *  Relies on TwoMeshes.F90 for mesh migration and interpolation 
! *  following calving, and FrontDisplacement.F90 for mesh update computation
! *  
! ******************************************************************************
! *
! *  Authors: Joe Todd
! *  Email:   jat39@st-andrews.ac.uk
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 10.11.2014
! *
! ****************************************************************************/
SUBROUTINE Find_Calving (Model, Solver, dt, TransientSimulation )

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE CrevasseGroups_t
      LOGICAL, ALLOCATABLE :: NotEmpty(:),Valid(:)
      INTEGER, ALLOCATABLE :: NodeIndexes(:,:), NoNodes(:)
      INTEGER :: CurrentGroup
      !This derived type is for storing groups of connected nodes which
      !have surface or basal crevasses.  I can't think of a sensible strategy
      !to determine the number of groups required, nor the number of nodes in
      !each group.  However, given that they only contain integers and logicals,
      !declaring them massive probably isn't an issue.
   END TYPE CrevasseGroups_t

   TYPE(Model_t) :: Model
   TYPE(Solver_t) :: Solver
   Type(Nodes_t) :: CornerElementNodes, CurrentElementNodes, &
        TargetNodes, CalvedNodes
   Type(Nodes_t), POINTER :: Nodes0
   Type(Nodes_t), TARGET :: EvalNodes
   TYPE(Matrix_t), POINTER :: Vel1Matrix !For finding neighbours
   TYPE(ValueList_t), POINTER :: Material, SolverParams
   TYPE(Element_t), POINTER :: CurrentElement, CornerElement
   TYPE(CrevasseGroups_t) :: SurfaceCrevasseGroups, BasalCrevasseGroups
   TYPE(Mesh_t), POINTER :: Mesh, EvalMesh, Mesh0 => Null()
   TYPE(Variable_t), POINTER :: DepthSol, StressSol, Stress1Sol, Stress4Sol, &
        Vel1Sol, DistanceSol, FlowSol, &
        WPressureSol, TimeVar, CSurfIndexSol, CBasalIndexSol, Var, &
        CrevasseGroupVar

   REAL (KIND=dp) :: t, dt, xcoord, ycoord, NodeDepth, NodeStress, NodeStress1, NodeStress4, &
        NodeWPressure, CalvingSurfIndex, CalvingBasalIndex, CalvingCoordHolder = 1.0E20, &
        FreeConnect, FreeConnected, Length, WaterDepth, Dw, NodeLength, RhoWF, RhoWS, Rho, &
        g, OverlapCalvingCoordinate, SeaLevel, BasalCalvingCoordinate, EvalResolution, &
        STuningParam, BTuningParam, DwStart, DwStop, season, Te, sign, Normal(3), &
        CornerNormal(3), BedSecond, BedSecondDiff, beddiff, BedToler, dx, dy, &
        LocalDist, LocalDistNode, PropDistNode, normalcond, ForceCalveSize, &
        localM(2,2), EigValues(2), EI(2), dumy, dumy2, work(8),YieldStress, &
#ifdef USE_ISO_C_BINDINGS
        rt0, rt
#else
        rt0, rt, RealTime
#endif
   REAL (KIND=DP), POINTER :: DepthValues(:), StressValues(:), Stress1Values(:), Stress4Values(:), &
        Calving1Values(:), DistanceValues(:), WPressureValues(:), FrontValues(:), &
        CSurfIndexValues(:), CBasalIndexValues(:), CIndexValues(:), CrevasseGroupValues(:), &
        Calving2Values(:), EvalPoints(:,:), Field(:), WorkReal(:)
   REAL (KIND=dp), ALLOCATABLE :: CumDist(:), PropCumDist(:), TargetCumDist(:),TargetPropDist(:)
   INTEGER :: DIM, i, j, NoNodes, MaxN, FrontNodes, BotNodes, TopNodes, &
        NofEvalPoints, NextSlot, BotCornerIndex, &
        BotSecondIndex, county, GoToNode, PrevNode, NextNode, &
        TimeSinceLast, NoNeighbours, MaxNeighbours, DOFs, ierr, StressDOFs

   INTEGER, POINTER :: DepthPerm(:), StressPerm(:), Stress1Perm(:), Stress4Perm(:),&
        Vel1Perm(:), Vel1InvPerm(:), DistancePerm(:), OrderPerm(:),&
        WPressurePerm(:), FrontPerm(:)=>NULL(), InvFrontPerm(:), &
        TopPerm(:)=>NULL(), BotPerm(:)=>NULL(), Permutation(:), &
        CSurfIndexPerm(:), CBasalIndexPerm(:), CIndexPerm(:), FieldPerm(:), &
        NodeNeighbours(:,:), NumNeighbours(:), CrevasseGroupPerm(:)

   INTEGER, ALLOCATABLE :: ThisNodeNeighbours(:)
   LOGICAL :: FirstTime = .TRUE., CalvingOccurs, RemeshOccurs, &
        BasalCalvingOccurs = .FALSE., BasalCrevasseModel, OverlapOccurs, DwSeason, &
        CornerCalving, KeepLooking, TransientSimulation, Found, Debug = .FALSE., &
        ForceCalving = .FALSE., PlaneStressOnly, Cauchy, OldWay, BasalFS, CornerBadBed, &
        CornerBadSlope
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FrontMaskName, BotMaskName, TopMaskName, &
        DwMode, CalvingModel, FlowVarName,BasalFSVarName

   SAVE :: FirstTime, DIM, Permutation, Calving1Values, Calving2Values, &
        NoNodes, Mesh, FrontMaskName, FreeConnect, WaterDepth, Rho, RhoWS, RhoWF, g, &
        FrontPerm, FrontValues, FrontNodes, BotPerm, BotNodes, TopPerm, &
        TopNodes, CSurfIndexSol, &
        CBasalIndexSol, CSurfIndexPerm, CBasalIndexPerm, CSurfIndexValues, CBasalIndexValues, &
        EvalMesh, Material, EvalNodes, NodeNeighbours, NumNeighbours, Mesh0, Nodes0, &
        TimeSinceLast, BasalCrevasseModel, PlaneStressOnly, Cauchy, SolverParams, OldWay, &
        BasalFS, BasalFSVarName

   INTERFACE
      SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
           NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
        !------------------------------------------------------------------------------
        USE Lists
        USE SParIterComm
        USE Interpolation
        USE CoordinateSystems
        !-------------------------------------------------------------------------------
        TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
        TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
        LOGICAL, OPTIONAL :: UseQuadrantTree
        LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
        TYPE(Projector_t), POINTER, OPTIONAL :: Projector
        CHARACTER(LEN=*),OPTIONAL :: MaskName
      END SUBROUTINE InterpolateMeshToMesh
   END INTERFACE

   rt0 = RealTime() !time it!

   SolverName = "Calving"

   Debug = .FALSE.

   IF(ParEnv % Pes > 1) CALL Fatal(SolverName, "Solver doesnt work in parallel!")

   Timevar => VariableGet( Model % Variables,'Time', UnfoundFatal=.TRUE.)
   t = TimeVar % Values(1)
   dt = Model % Solver % dt
   season = t - floor(t)

   RemeshOccurs = .FALSE. !For undercutting

   !*****************************Get Variables*******************************

   !This is the main solver variable, which only has a value on the calving
   !face of the glacier, and corresponds to the magnitude of retreat.
   IF(.NOT. ASSOCIATED(Solver % Variable)) CALL Fatal('Calving', 'Variable not associated!')

   !This is the Calving Solver's Surface Crevasse exported variable, which has
   !a positive value where a crevasse is predicted to exist, and a negative value elsewhere.
   CSurfIndexSol => VariableGet( Solver % Mesh % Variables, 'Calving Surface Index', &
        UnfoundFatal=.TRUE.)
   CSurfIndexPerm => CSurfIndexSol % Perm
   CSurfIndexValues => CSurfIndexSol % Values

   !This is the Calving Solver's Basal crevasse exported variable.
   !Positive value where basal crevasse
   !is predicted to exist, negative otherwise, and zero at the tip.

   !From Faezeh: C_basal = LongitudinalDevStress(Rxx)-Rho*g*(IceThickness-
   !HeightAboveGlacierBase)+RhoOcean*g*(PiezometricHeight-HeightAboveGlacierBase)

   CBasalIndexSol => VariableGet( Solver % Mesh % Variables, 'Calving Basal Index', &
        UnfoundFatal=.TRUE.)
   CBasalIndexPerm => CBasalIndexSol % Perm
   CBasalIndexValues => CBasalIndexSol % Values

   CrevasseGroupVar => VariableGet( Solver % Mesh % Variables, 'Crevasse Group ID', &
        UnfoundFatal=.TRUE.)
   CrevasseGroupPerm => CrevasseGroupVar % Perm
   CrevasseGroupValues => CrevasseGroupVar % Values

   DepthSol => VariableGet( Model % Variables, 'Depth', UnfoundFatal=.TRUE.)
   DepthPerm => DepthSol % Perm
   DepthValues => DepthSol % Values

   StressSol => VariableGet( Model % Variables, "Stress", UnfoundFatal=.TRUE.)
   StressPerm => StressSol % Perm
   StressValues => StressSol % Values
   StressDOFs = StressSol % DOFs

   Stress1Sol => VariableGet( Model % Variables, "Stress 1", &
        UnfoundFatal=.TRUE.)
   Stress1Perm => Stress1Sol % Perm
   Stress1Values => Stress1Sol % Values

   Stress4Sol => VariableGet( Model % Variables, "Stress 4", &
        UnfoundFatal=.TRUE.)
   Stress4Perm => Stress4Sol % Perm
   Stress4Values => Stress4Sol % Values

   DistanceSol => VariableGet( Model % Variables, "Distance", &
        UnfoundFatal=.TRUE.)
   DistancePerm => DistanceSol % Perm
   DistanceValues => DistanceSol % Values

   WPressureSol => VariableGet( Model % Variables, "Water Pressure", &
        UnfoundFatal=.TRUE.)
   WPressurePerm => WPressureSol % Perm
   WPressureValues => WPressureSol % Values

   IF(FirstTime) THEN

      DIM = CoordinateSystemDimension()
      IF(DIM /= 2) CALL Fatal('Calving','Solver only works in 2D!')
      MaxN = Model % Mesh % MaxElementNodes
      !To track time since last calving event. Probably not needed.
      TimeSinceLast = 0

      Mesh => Solver % Mesh
      NoNodes = SIZE( Mesh % Nodes % x )

      FrontMaskName = 'Calving Front Mask'
      TopMaskName = 'Top Surface Mask'
      BotMaskName = 'Bottom Surface Mask'
      ALLOCATE( FrontPerm(NoNodes), TopPerm(NoNodes), BotPerm(NoNodes))

      CALL MakePermUsingMask( Model,Solver,Mesh,FrontMaskName, &
           .FALSE., FrontPerm, FrontNodes )
      CALL MakePermUsingMask( Model,Solver,Mesh,TopMaskName, &
           .FALSE., TopPerm, TopNodes )
      CALL MakePermUsingMask( Model,Solver,Mesh,BotMaskName, &
           .FALSE., BotPerm, BotNodes )

      !Holds the variable values
      ALLOCATE( FrontValues(FrontNodes * 2) )

      Mesh0 => AllocateMesh()
      Mesh0 = Mesh
      Mesh0 % Name = TRIM(Mesh % Name)//'_initial'
      CALL Info('Calving','Created initial reference mesh to remap front, maintaining quality')
      ALLOCATE( Nodes0 )
      ALLOCATE( Nodes0 % x(NoNodes), Nodes0 % y(NoNodes), Nodes0 % z(NoNodes) )
      Nodes0 % x = Mesh % Nodes % x
      Nodes0 % y = Mesh % Nodes % y
      Nodes0 % z = Mesh % Nodes % z
      Mesh0 % Nodes => Nodes0

      Permutation => FrontPerm
      Calving1Values => FrontValues(1::2)
      Calving2Values => FrontValues(2::2)

      !Initialize
      Calving1Values = 0.0_dp
      Calving2Values = 0.0_dp

      !Potential to force an initial calving event
      !Specified by LOGICAL and REAL in SIF
      SolverParams => GetSolverParams()
      ForceCalving = GetLogical(SolverParams, 'Force Calving', Found)
      !TODO: This doesn't work - Calving1Values isn't associated with Solver % Variable % Values yet...
      IF(ForceCalving) THEN
         ForceCalveSize = GetConstReal (SolverParams, 'Force Calving Size', Found)
         IF(.NOT.Found) CALL Fatal(SolverName, "Force Calving Size not found.")
         Calving1Values = -ForceCalveSize
         CalvingOccurs =.TRUE.
         CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )
         RETURN
      END IF
   ENDIF

   Solver % Variable % Values => FrontValues
   Solver % Variable % Perm => FrontPerm

   Var => VariableGet(Solver % Mesh % Variables, ComponentName(Solver % Variable % Name, 1), .TRUE.)
   Var % Values => Calving1Values
   Var % Perm => Permutation
   Var => VariableGet(Solver % Mesh % Variables, ComponentName(Solver % Variable % Name, 2), .TRUE.)
   Var % Values => Calving2Values
   Var % Perm => Permutation

   IF(FirstTime .OR. Solver % Mesh % Changed) THEN
      FirstTime = .FALSE.

      !STRATEGY: Finding neighbours on the fly works fine UNLESS you are in a recursive subroutine
      !Then it messes up, because at each level, ThisNodeNeighbours is deallocate and reallocated,
      !meaning that when you jump back up, info is already overwritten. SO:
      !Keep the structure as it was with CalvingNeighbours, cycle as below to fill it, and ALSO
      !create an array to hold the number of neighbours for each node

      !Get the Matrix of the N-S Solver
      Vel1Sol => VariableGet( Solver % Mesh % Variables, 'Velocity 1', UnfoundFatal=.TRUE.)
      Vel1Perm => Vel1Sol % Perm
      Vel1Matrix => Vel1Sol % Solver % Matrix

      !Vel * DIM + Pressure...
      DOFs = DIM + 1
      MaxNeighbours = DIM * 10  !totally arbitrary...

      !Create inverse perm to lookup Matrix later
      ALLOCATE(Vel1InvPerm(MAXVAL(Vel1Perm)*DOFs)) !TODO DEALLOCATE
      !2D array to hold each nodes neighbours
      ALLOCATE(NodeNeighbours(NoNodes,MaxNeighbours))
      !1D array to hold number of neighbours for each node
      ALLOCATE(NumNeighbours(NoNodes))
      NodeNeighbours = 0
      NumNeighbours = 0
      Vel1InvPerm = 0

      j = 0
      DO i=1,SIZE(Vel1Perm)
         IF(Vel1Perm(i) == 0) CYCLE
         j = j + 1
         Vel1InvPerm( (Vel1Perm(i)*DOFs-2) : (Vel1Perm(i)*DOFs) ) = j !The 2 here is suspect...
      END DO

      DO i = 1,NoNodes
         CALL FindNodeNeighbours(i) !Updates the allocatable array 'ThisNodeNeighbours'
         NumNeighbours(i) = SIZE(ThisNodeNeighbours)
         NodeNeighbours(i,1:NumNeighbours(i)) = ThisNodeNeighbours
      END DO
   END IF

   PRINT *, '****Calving Timestep : ',t

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Calving, time taken for variable loading etc:', rt
   rt0 = RealTime()

   MaxN = Model % Mesh % MaxElementNodes
   OldWay = ListGetLogical(SolverParams, "Old Way", Found)
   IF(.NOT. Found) OldWay = .FALSE.

   !Identify the basal freesurface variable
   BasalFS = ListGetLogical(SolverParams, "Basal FreeSurface", Found)
   IF(.NOT. Found) BasalFS = .TRUE.
   IF(BasalFS) THEN
     BasalFSVarName = ListGetString(SolverParams, "Basal FreeSurface Variable Name",Found)
     IF(.NOT. Found) CALL Fatal(SolverName, &
          "Basal FreeSurface = True but no Basal FreeSurface Variable Name found!")
   END IF

   Material => GetMaterial(Model % Mesh % Elements(Solver % ActiveElements(1)))
   Cauchy = ListGetLogical( Material , 'Cauchy', Found )
   IF(.NOT. Found) THEN
     Cauchy = .FALSE.
     CALL Warn(SolverName, "Couldn't find 'cauchy' logical in material, &
          &assuming deviatoric stress.")
   END IF

   FlowVarName = ListGetString(SolverParams, "Flow Solution Variable Name",Found)
   IF(.NOT. Found) FlowVarName = "Flow Solution"

   FlowSol => VariableGet(Model % Variables, FlowVarName, .TRUE., UnfoundFatal=.TRUE.)

   !! FreeConnect: distance from calving front where free connection to proglacial
   !!water body exists. Declare in SIF Constants!
   FreeConnect = GetConstReal( Model % Constants, "FreeConnect", Found )
   IF(.NOT.Found) CALL Fatal(SolverName, "Unable to find FreeConnect")
   Dw = GetConstReal( Model % Constants, "Water Depth", Found )
   IF(.NOT.Found) THEN
     CALL Warn(SolverName, "Unable to find Water Depth, assuming zero.")
     Dw = 0.0_dp
   END IF

   Rho = GetConstReal( Model % Constants, "Rho", Found )
   IF(.NOT.Found) CALL Fatal(SolverName, "Unable to find Rho.")
   RhoWS = GetConstReal( Model % Constants, "RhoWS", Found )
   IF(.NOT.Found) CALL Fatal(SolverName, "Unable to find RhoW.")
   RhoWF = GetConstReal( Model % Constants, "RhoWF", Found )
   IF(.NOT.Found) CALL Fatal(SolverName, "Unable to find RhoWF.")
   g = GetConstReal( Model % Constants, "g", Found )
   IF(.NOT.Found) CALL Fatal(SolverName, "Unable to find Gravity.")
   IF(g < 0.0_dp) g = -1.0_dp * g

   SeaLevel = GetConstReal( Model % Constants, "Sea Level", Found )
   IF(.NOT.Found) THEN
     WRITE(Message,'(A)') 'Variable Sea level not found. &
          &Setting to 0.0'
     CALL Info('Calving', Message, level=2)
     SeaLevel = 0.0_dp
   END IF

   !This next bit stolen from MaterialModels.src
   !Only needed if calvingmodel=planestrain, could put logical?
   STuningParam = GetConstReal( Model % Constants, "Surface Crevasse Tuning Parameter", Found )
   IF(.NOT.Found) THEN
     WRITE(Message,'(A)') 'Surface Crevasse Tuning Parameter not found. &
          Setting to 1.0'
     CALL Info('Calving', Message, level=2)
     STuningParam = 1.0_dp
   END IF

   BTuningParam = GetConstReal( Model % Constants, "Basal Crevasse Tuning Parameter", Found )
   IF(.NOT.Found) THEN
     WRITE(Message,'(A)') 'Basal Crevasse Tuning Parameter not found. &
          Setting to 1.0'
     CALL Info('Calving', Message, level=2)
     BTuningParam = 1.0_dp
   END IF

   !Get stuff for variable water depth in crevasses.
   !Don't think this is the logical place for this, but whatever.
   SolverParams => GetSolverParams()

   DwMode = GetString (SolverParams, 'Dw Mode', Found)
   IF(.NOT.Found) THEN
     CALL Warn(SolverName, "Dw Mode not found, assuming off.")
     DwMode = "off"
   END IF
   DwMode = TRIM(DwMode)
   IF(DwMode /= "off") THEN
     DwStart = GetConstReal (SolverParams, 'Dw Start', Found)
     IF(.NOT.Found) CALL Fatal(SolverName, "Dw Start not found.")

     DwStop = GetConstReal (SolverParams, 'Dw Stop', Found)
     IF(.NOT.Found) CALL Fatal(SolverName, "Dw Stop not found.")
   END IF

   YieldStress = ListGetConstReal(SolverParams, "Yield Stress", Found)
   IF(.NOT. Found) THEN
     YieldStress = 0.0_dp
     CALL Warn(SolverName, "No yield stress found, assuming 0.0")
   ELSE
     WRITE(Message,'(A,f8.2)') 'Yield Stress = ',YieldStress
     CALL Info(SolverName, Message)
   END IF

   BasalCrevasseModel = ListGetLogical(SolverParams, "Basal Crevasse Model", Found)
   IF(.NOT. Found) THEN
     BasalCrevasseModel = .TRUE.
     CALL Warn(SolverName,"No 'Basal Crevasse Model' found, assuming True")
   END IF
   IF(BasalCrevasseModel) THEN
     CALL Info(SolverName,"Using Basal Crevasse Model, &
          &calving occurs when surface and basal crevasses meet", Level=2)
   ELSE
     CALL Info(SolverName,"Using Surface Crevasse Model, &
          &calving occurs when surface crevasses meet sea level",Level=2)
   END IF

   PlaneStressOnly = ListGetLogical(SolverParams, "Plane Stress Only", Found)
   IF(.NOT. Found) PlaneStressOnly = .FALSE.

   SELECT CASE(DwMode)
   CASE("off")
     WaterDepth = 0.0
   CASE ("constant")

     WaterDepth = Dw

   CASE("binary")

     DwSeason = .FALSE.
     IF(DwStart .LT. DwStop) THEN !normal case, dw season doesn't pass winter
       IF((season .GT. DwStart) .AND. (season .LT. DwStop)) DwSeason = .TRUE.
     ELSE !this is unlikely to be needed
       IF((season .GT. DwStart ) .OR. (season .LT. DwStop )) DwSeason = .TRUE.
     END IF
     IF(DwSeason) THEN
       WaterDepth = Dw
     ELSE
       WaterDepth = 0.0
     ENDIF

   CASE DEFAULT
     CALL Fatal(SolverName,"Invalid Dw mode selection")
   END SELECT

   !Finds the maximum value of Coordinate 1 (i.e. the calving face) and
   !FreeConnected (the minimum Coordinate 1 which should be checked for calving.

   !First find the length of the glacier
   Length = 0.0
   DO i = 1, NoNodes
     IF(Permutation(i) == 0) CYCLE
     NodeLength = Model % Nodes % x(i)
     IF(NodeLength > Length) Length = NodeLength
   END DO
   PRINT *, '**** Calving Glacier Length = ',Length
   CalvingCoordHolder = Length
   FreeConnected = Length - FreeConnect

   !---------------------------------------
   !
   ! Evaluate the crevasse criteria for basal and surface crevasses
   !
   !---------------------------------------

   DO i = 1, NoNodes
     xcoord = Mesh % Nodes % x(i)
     ycoord = Mesh % Nodes % y(i)

     !TODO - Tuning Parameters aren't used except in specific PlaneStressOnly case.
     !SORT THIS OUT - but what sense does the tuning parameter have in the cauchy case?
     IF(.NOT. OldWay) THEN

       !Compute maximum extensive principal stress
       localM=0.0_dp
       !xx,yy,zz,xy,yz,zx
       localM(1,1)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 1) !xx
       localM(2,2)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 2) !yy
       localM(1,2)=StressValues( StressDOFs * (StressPerm(i) - 1 ) + 4) !xy
       localM(2,1)=localM(1,2)

       IF(.NOT. Cauchy) THEN
         DO j=1,2                                       !dim+1, only runs in 3D
           localM(j,j)=localM(j,j) - FlowSol % Values(FlowSol % Perm(i)*FlowSol % DOFs)
         END DO
       END IF

       CALL DGEEV('N','N',2,localM,2,EigValues,EI,Dumy,1,Dumy2,1,Work,8,ierr )
       IF (ierr/=0) THEN
         WRITE(Message,'(A,i0)') 'Failed to compute EigenValues, error code: ',ierr
         CALL FATAL(SolverName, Message)
       END IF

       NodeStress = MAXVAL(EigValues) !Maximum principalstress

       NodeWPressure = -WPressureValues(WPressurePerm(i))
       IF(NodeWPressure < 0.0_dp) NodeWPressure = 0.0_dp

       CalvingSurfIndex = NodeStress + (WaterDepth * RhoWF * g) - YieldStress

       CalvingBasalIndex = NodeStress + NodeWPressure - YieldStress
     ELSE
       CALL Info('Calving', 'Calculating calving criterion the old way', level=2)
       NodeDepth = DepthValues(DepthPerm(i))

       NodeStress1 = Stress1Values(Stress1Perm(i))
       IF(Cauchy) &
            NodeStress1 = NodeStress1 + FlowSol % Values(FlowSol % Perm(i)*FlowSol % DOFs)


       NodeWPressure = -WPressureValues(WPressurePerm(i))
       IF(NodeWPressure < 0.0_dp) NodeWPressure = 0.0_dp

       !NodeWPressure is the water pressure in a basal crevasse (piezometric - height) basically
       !From Faezeh:
       !C_surface = LongitudinalDevStress(Rxx) - Rho*g*NodeDepth + WaterDepth*RhoW*g

       !C_basal = LongitudinalDevStress(Rxx)-Rho*g*(IceThickness-HeightAboveGlacierBase)+RhoOcean*g*(PiezometricHeight-HeightAboveGlacierBase)
       !Water pressure (Rho_w * g * Piezo) is PressureSol,PressureValues,PressurePerm

       IF(PlaneStressOnly) THEN

         CalvingSurfIndex = ( STuningParam * 2.0 * NodeStress1 ) - &
              (NodeDepth * g * Rho) + (WaterDepth * RhoWF * g) - YieldStress
         CalvingBasalIndex = (  BTuningParam * 2.0 * NodeStress1 ) - &
              (NodeDepth * g * Rho) + NodeWPressure - YieldStress
         !Last term appears if you split the brackets on the last term in C_Basal from Faezeh
       ELSE
         !Shear stress for Te calc
         NodeStress4 = Stress4Values(Stress4Perm(i))
         !2(Te^2) = Txx^2 + Tyy^2 + 2*Txy^2
         !Txx = Tyy thus:
         !Te^2 = (Txx^2) + (Txy^2)
         !Te = ((Txx^2) + (Txy^2)) ^ 0.5
         Te = ((NodeStress1**2) + (NodeStress4**2)) ** 0.5
         sign = NodeStress1/abs(NodeStress1)

         CalvingSurfIndex = ( STuningParam * 2.0 * sign * Te ) - &
              (NodeDepth * g * Rho) + (WaterDepth * RhoWF * g) - YieldStress
         CalvingBasalIndex = (  BTuningParam * 2.0 * sign * Te ) - &
              (NodeDepth * g * Rho) + NodeWPressure - YieldStress
       END IF

     END IF

     CSurfIndexValues(CSurfIndexPerm(i)) = CalvingSurfIndex
     CBasalIndexValues(CBasalIndexPerm(i)) = CalvingBasalIndex
   END DO

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Calving, time taken for evaluating CIndex:', rt
   rt0 = RealTime()

   !---------------------------------------
   !
   ! Find connected crevassed regions and check for calving
   !
   !---------------------------------------

   !TODO, allow both
   IF(BasalCrevasseModel) THEN
     CrevasseGroupValues = 0

     !Find groups of nodes which have surface crevasses
     CIndexValues => CSurfIndexValues
     CIndexPerm => CSurfIndexPerm
     CALL FindCrevasseGroups(SurfaceCrevasseGroups,.TRUE., TopPerm)

     !As above, basal crevasses
     CIndexValues => CBasalIndexValues
     CIndexPerm => CBasalIndexPerm
     CALL FindCrevasseGroups(BasalCrevasseGroups,.TRUE., BotPerm)

     !Look for touching/almost touching crevasse groups
     CALL FindCalvingBasal(SurfaceCrevasseGroups, BasalCrevasseGroups, BasalCalvingCoordinate)

   ELSE !Surface Crevasse Model

     !This dictates the resolution of the mesh for interpolating calving index
     !Only relevant for surface crevasse mode, as opposed to surf and basal
     !TODO: Unhardcode this
     EvalResolution = 1.0_dp

     !Create the points (at the waterline) at which Calving Index will be evaluated
     !and assign them to a new mesh for interpolation

     NofEvalPoints = FLOOR(FreeConnect/EvalResolution)
     EvalMesh => AllocateMesh()
     EvalMesh % Name = "Eval_Mesh"
     EvalMesh % NumberOfNodes = NofEvalPoints
     EvalMesh % Nodes => EvalNodes
     ALLOCATE(Field(NofEvalPoints), &
          FieldPerm(NofEvalPoints), &
          EvalPoints(NofEvalPoints,2), &
          EvalNodes % x(NofEvalPoints), &
          EvalNodes % y(NofEvalPoints), &
          EvalNodes % z(NofEvalPoints))

     Field = -1.0_dp
     EvalPoints(:,2) = SeaLevel
     DO i=1,NofEvalPoints
       EvalPoints(i,1) = Length - i*EvalResolution
       FieldPerm(i) = i
     END DO
     EvalMesh % Nodes % x = EvalPoints(:,1)
     EvalMesh % Nodes % y = SeaLevel
     EvalMesh % Nodes % z = 0.0_dp

     CALL VariableAdd( EvalMesh % Variables, EvalMesh, Solver, "Calving Surface Index", 1, Field, FieldPerm )
     !TEST - should be true
     CALL InterpolateMeshToMesh( Mesh, EvalMesh, Mesh % Variables, EvalMesh % Variables, .FALSE. )

     CalvingOccurs = .FALSE.
     Var => VariableGet( EvalMesh % Variables, 'Calving Surface Index', UnfoundFatal=.TRUE.)
     DO i=1,NofEvalPoints
       IF(Var % Values(Var % Perm(i)) < 0.0_dp ) CYCLE
       CalvingOccurs = .TRUE.
       PRINT *,'Surface Calving Event'
       IF(EvalMesh % Nodes % x(i) < CalvingCoordHolder) &
            CalvingCoordHolder = EvalMesh % Nodes % x(i)
     END DO

   END IF !Surface crevasse model

   IF(BasalCrevasseModel) THEN
      CalvingOccurs = .FALSE.
      IF(BasalCalvingOccurs) THEN
         CalvingOccurs = .TRUE.
         DO j = 1, NoNodes
            IF(Permutation(j) == 0) CYCLE
            Calving1Values(Permutation(j)) = BasalCalvingCoordinate - &
                 Mesh % Nodes % x(j)
            IF(Calving1Values(Permutation(j)) .GE. 0.0_dp) &
                 Calving1Values(Permutation(j)) = 0.0_dp
         END DO
         PRINT *, 'Basal Calving Event at timestep ',t
         PRINT *, 'Calving Coordinate :',BasalCalvingCoordinate
      ELSE
         Calving1Values = 0.0_dp
      END IF
   ELSE
      !Surface-only calving model
      IF(CalvingOccurs) THEN
         DO j = 1, NoNodes
            IF(Permutation(j) == 0) CYCLE
            Calving1Values(Permutation(j)) = CalvingCoordHolder - &
                 Mesh % Nodes % x(j)
            IF(Calving1Values(Permutation(j)) .GE. 0.0_dp) &
                 Calving1Values(Permutation(j)) = 0.0_dp
         END DO
         PRINT *, 'Calving Event at timestep ',t
         PRINT *, 'Calving Coordinate :',CalvingCoordHolder
      ELSE
         Calving1Values = 0.0_dp
      END IF
   END IF

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Calving, time taken for calving connectivity:', rt
   rt0 = RealTime()

   !---------------------------------------
   !   CALVING DONE
   !---------------------------------------

   !At this point, the 'calving' solution is done.  However...
   !Here we solve a problem to do with undercutting:
   !Progressive undercutting by melting of the calving front
   !can lead to a situation where 'Front' nodes start to look like
   !basal nodes.  However, they are 'officially' front nodes and so
   !don't have a friction law, grounding dynamics OR (most importantly
   !I think), a bed constraint.  Thus, it is necessary to check for this
   !occurring, and shift the bed nodes appropriately.
   !
   !Strategy:
   !-Identify the corner node by BotPerm and FrontPerm
   !
   !-Use BCelement connections to find the second to bottom
   !    node on the calving front
   !
   !-Check a condition: either
   !    second node is 'near' bed OR
   !    BCelement slope is below some critical level
   !
   !-If condition is met (i.e. we need to take action)
   !    Calving1Values @CornerNode = (X@2nd - X@corner)
   !    CalvingOccurs = .TRUE.  <-- will this have any unforeseen consequences?
   !
   !NOTE: This works in tandem with a section of TwoMeshes.f90 which does the actual
   !deformation

   !Get the node index of the bottom corner
   !NOTE: this could be 'FirstTime'd if it was also 'SAVE'd
   DO i=1,NoNodes
      IF(BotPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
         BotCornerIndex = i
      END IF
   END DO

   !Loop boundary elements, we're looking for the BCelement
   !containing BotCornerIndex and ANOTHER FrontPerm node...
   DO i=Mesh % NumberOfBulkElements+1,&
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      CurrentElement => Mesh % Elements(i)
      IF(.NOT.(ANY(CurrentElement % nodeindexes == BotCornerIndex))) CYCLE
      IF(ALL(FrontPerm(CurrentElement % nodeindexes) .GT. 0)) THEN
         !We have a winner
         CornerElement => Mesh % Elements(i)
         DO j=1,2
            IF(CurrentElement % NodeIndexes(j) .NE. BotCornerIndex) THEN
               BotSecondIndex = CurrentElement % NodeIndexes(j)
            END IF
         END DO
      END IF
   END DO

   !Check corner node isn't already calving
   IF(Calving1Values(Permutation(BotCornerIndex)) .LT. 0.0_dp) THEN
      CornerCalving = .TRUE.
   ELSE
      CornerCalving = .FALSE.
   END IF

   !Get normal vector:
   CALL GetElementNodes(CornerElementNodes, CornerElement)
   CornerNormal = NormalVector(CornerElement, CornerElementNodes)

   IF(Debug) PRINT *, 'Debug Calving, corner normal is: ' , &
        CornerNormal(1), CornerNormal(2), CornerNormal(3)

   IF(BasalFS) THEN
     BedSecond = ListGetRealAtNode( Material, "Min "//BasalFSVarName, &
          BotSecondIndex, UnfoundFatal=.TRUE. )

     IF(Debug) PRINT *, 'Debug Calving, second node bed is: ',&
          BedSecond,' and y coord is: ', Model % Nodes % y(BotSecondIndex)

     PRINT *, 'Debug Calving, second node bed is: ',&
          BedSecond,' and y coord is: ', Model % Nodes % y(BotSecondIndex)

     BedSecondDiff = Model % Nodes % y(BotSecondIndex) - BedSecond
   END IF

   !TODO - unhardcode these
   BedToler = 2.0_dp
   normalcond = 0.95_dp

   CornerBadBed = BasalFS .AND. (BedSecondDiff < BedToler)
   CornerBadSlope = ABS(CornerNormal(2)) > normalcond

   !If the slope normal is above threshold, or the second node is too close to the bed,
   !move the corner node to the second, via 'calving'
   IF((CornerBadSlope .OR. CornerBadBed) .AND. (.NOT. CornerCalving)) THEN

      IF(Debug) PRINT *,'Debug Calving, migrating mesh'
      county = 1
      GoToNode = BotSecondIndex
      PrevNode = BotCornerIndex
      KeepLooking = .TRUE.
      DO WHILE (KeepLooking)
         !Check if we should shift more than one node forward...
         IF(Debug) PRINT *, 'Debug calving: looking!'
         KeepLooking = .FALSE.
         DO i=Mesh % NumberOfBulkElements+1,&
              Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
            CurrentElement => Mesh % Elements(i)
            IF(.NOT.(ALL(FrontPerm(CurrentElement % nodeindexes) .GT. 0))) CYCLE
            !This point reached by Front BC elements, stupid way of doing it, but whatever
            IF(.NOT.(ANY(CurrentElement % nodeindexes == GoToNode))) CYCLE
            IF(ANY(CurrentElement % nodeindexes == PrevNode)) CYCLE
            !We only get here if element is the next one along from previous
            CALL GetElementNodes(CurrentElementNodes, Currentelement)
            Normal = NormalVector(CurrentElement, CurrentElementNodes)
            DO j=1,2
              IF(CurrentElement % NodeIndexes(j) .NE. GoToNode) &
                   NextNode = CurrentElement % NodeIndexes(j)
            END DO

            IF(BasalFS) THEN
              beddiff = Model % Nodes % y(NextNode) - ListGetRealAtNode( Material, &
                   "Min "//BasalFSVarName, NextNode, UnfoundFatal=.TRUE. )
            END IF

            CornerBadBed = BasalFS .AND. (beddiff < BedToler)
            CornerBadSlope = ABS(CornerNormal(2)) > normalcond

            IF(CornerBadBed .OR. CornerBadSlope) THEN
               PrevNode = GoToNode
               GoToNode = NextNode
               county = county + 1
               IF(Debug) PRINT *, 'Debug calving: Found another shift'
               KeepLooking = .TRUE.
               EXIT
            END IF
         END DO
      END DO

      RemeshOccurs = .TRUE.
   ELSE
      county = 0
   END IF

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Calving, time taken for corner problems:', rt
   rt0 = RealTime()

   IF(CalvingOccurs .OR. RemeshOccurs) THEN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Find the FrontDisplacement (= Calving 2 <-sif) for each frontal node
      !resulting from the shift in the corner node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Use FrontPerm to construct ordered node list
      !I think MakeMaskUsingPerm already ordered the nodes, from the comments:
      !! The bandwidth optimization for lines results to perfectly ordered
      !! permutations. If there is only one line the 1st node should be the
      !! lower left corner.
      ALLOCATE(InvFrontPerm(FrontNodes),&
           CumDist(FrontNodes),&
           PropCumDist(FrontNodes),&
           TargetCumDist(FrontNodes),&
           TargetPropDist(FrontNodes),&
           TargetNodes % x(FrontNodes),&
           TargetNodes % y(FrontNodes),&
           TargetNodes % z(FrontNodes),&
           CalvedNodes % x(FrontNodes))

      !InvFrontPerm(FrontNodes) points to nodeindexes in order they appear...
      !InvFrontPerm(1) points to 'lower left' corner, according to MakePermUsingMask
      DO i=1,NoNodes
         IF(FrontPerm(i) .GT. 0) THEN
            InvFrontPerm(FrontPerm(i)) = i
         END IF
      END DO

      ALLOCATE(OrderPerm(FrontNodes), WorkReal(FrontNodes))
      OrderPerm = [(i,i=1,FrontNodes)]
      DO i=1,FrontNodes
         WorkReal(i) = Mesh0 % Nodes % y(InvFrontPerm(i))
      END DO
      CALL SortD( FrontNodes, WorkReal, OrderPerm )
      DEALLOCATE(WorkReal)

      DO i=1,FrontNodes
         j = InvFrontPerm(OrderPerm(i))
         CalvedNodes % x(i) = Mesh % Nodes % x(j) + Calving1Values(Permutation(j))
         IF(Debug) PRINT *,'Debug Calving, CalvedNodes node: ',j,' is ',&
              CalvedNodes % x(i),' init coord: ',Mesh % Nodes % x(j)
      END DO

      !cycle through in order
      !First get target distribution from Mesh0
      TargetCumDist(1) = 0.0_dp
      DO i=2,FrontNodes
        dx = Mesh0 % Nodes % x(InvFrontPerm(OrderPerm(i))) - Mesh0 % Nodes % x(InvFrontPerm(OrderPerm(i-1)))
        dy = Mesh0 % Nodes % y(InvFrontPerm(OrderPerm(i))) - Mesh0 % Nodes % y(InvFrontPerm(OrderPerm(i-1)))

        TargetCumDist(i) = TargetCumDist(i-1) + (((dx**2) + (dy**2)) ** 0.5_dp)
      END DO
      TargetPropDist = TargetCumDist / MAXVAL(TargetCumDist)

      !Now find the length segments of our current line
      !If RemeshOccurs (because of bad corner node), county dictates
      !the offset from the previous bottom node of the new calving front
      CumDist(1:county+1) = 0.0_dp
      DO i=county+2,FrontNodes
        !sum coord magnitude from base upwards to give front 'length'
        !keep cumulative total
        !allocate proporitional y (and x) distances (i.e. out of 1)
        dx = CalvedNodes % x(i) - CalvedNodes % x(i-1)
        dy = Mesh % Nodes % y(InvFrontPerm(OrderPerm(i))) - Mesh % Nodes % y(InvFrontPerm(OrderPerm(i-1)))

        CumDist(i) = CumDist(i-1) + (((dx**2) + (dy**2)) ** 0.5_dp)
        !Remember first one is corner node...
        IF(Debug) PRINT *, 'Debug Calving: CumDist at node: ',&
             InvFrontPerm(OrderPerm(i)),' is ',CumDist(i)
        IF(Debug) PRINT *, 'Debug Calving: TargetDist at node: ',&
             InvFrontPerm(OrderPerm(i)),' is ',TargetCumDist(i)
      END DO
      PropCumDist = CumDist / MAXVAL(CumDist)

      !Loop each front node
      TargetNodes % x(1) = CalvedNodes % x(county+1)
      TargetNodes % y(1) = Mesh % Nodes % y(InvFrontPerm(OrderPerm(county+1)))
      TargetNodes % x(FrontNodes) = CalvedNodes % x(FrontNodes)
      TargetNodes % y(FrontNodes) = Mesh % Nodes % y(InvFrontPerm(OrderPerm(FrontNodes)))

      DO i=2,FrontNodes-1
        !and find nearest two nodes to interpolate
        DO j=county+2,FrontNodes
          IF(PropCumDist(j) .GT. TargetPropDist(i)) THEN
            !lin interp between j and j-1
            LocalDist = PropCumDist(j) - PropCumDist(j-1)
            LocalDistNode = TargetPropDist(i) - PropCumDist(j-1)

            PropDistNode = LocalDistNode / LocalDist
            IF(Debug) PRINT *, 'Debug Calving: PropDist at node: ',&
                 InvFrontPerm(OrderPerm(i)),' is ',PropDistNode

            TargetNodes % x(i) = ((1 - PropDistNode) * CalvedNodes % x(j-1))  + &
                 (PropDistNode * CalvedNodes % x(j))
            TargetNodes % y(i) = ((1 - PropDistNode) * Mesh % Nodes % y(InvFrontPerm(OrderPerm(j-1))))  + &
                 (PropDistNode * Mesh % Nodes % y(InvFrontPerm(OrderPerm(j))))
            EXIT
          END IF
        END DO
      END DO

      !At this point, we have obtained, for each FrontNode, a TargetNode % x and y
      !Thus, it simply remains to compute the two components of the displacement

      !Calving 1 = Diff X  (New % x - Old % x)
      !Calving 2 = Diff Y  (New % y - Old % y)

      DO i=1,FrontNodes
         Calving1Values(Permutation(InvFrontPerm(OrderPerm(i)))) = TargetNodes % x(i) &
              - Mesh % Nodes % x(InvFrontPerm(OrderPerm(i)))

         Calving2Values(Permutation(InvFrontPerm(OrderPerm(i)))) = TargetNodes % y(i) &
              - Mesh % Nodes % y(InvFrontPerm(OrderPerm(i)))

         IF(Debug) THEN
            PRINT *,'Debug Calving: Node: ',InvFrontPerm(OrderPerm(i)),' pos x: ',&
                 Mesh % nodes % x(InvFrontPerm(OrderPerm(i))),&
                 ' pos y: ',Mesh % nodes % y(InvFrontPerm(OrderPerm(i)))
            PRINT *,'Moving to: x: ',TargetNodes % x(i),' y: ',TargetNodes % y(i)
            PRINT *,'Displacement 1: ',Calving1Values(Permutation(InvFrontPerm(OrderPerm(i)))),&
                 'Displacement 2: ',Calving2Values(Permutation(InvFrontPerm(OrderPerm(i))))
         END IF
       END DO
     END IF

   CALL ListAddLogical( Model % Simulation, 'CalvingOccurs', CalvingOccurs )
   CALL ListAddLogical( Model % Simulation, 'RemeshOccurs', RemeshOccurs )

   rt = RealTime() - rt0
   IF(ParEnv % MyPE == 0) &
        PRINT *, 'Calving, time taken for calculating displacements (end):', rt
   rt0 = RealTime()

 CONTAINS

   SUBROUTINE FindNodeNeighbours(NodeNumber)
     INTEGER :: NodeNumber, i, count

     NoNeighbours = Vel1Matrix % Rows((Vel1Perm(NodeNumber)*DOFs)+1) - Vel1Matrix % Rows(Vel1Perm(NodeNumber)*DOFs)
     IF(MOD(NoNeighbours, DOFs).NE. 0) CALL Fatal(SolverName,"This shouldn't have happened...")
     !Each neighbour appears once per DOF, and there's also the current node thus: (x/DOFS) - 1...
     NoNeighbours = (NoNeighbours / DOFs) - 1
     IF(NoNeighbours .GT. MaxNeighbours) CALL Fatal(SolverName,"Need more neighbour slots!")

     IF(ALLOCATED(ThisNodeNeighbours)) DEALLOCATE(ThisNodeNeighbours)
     ALLOCATE(ThisNodeNeighbours(NoNeighbours))
     ThisNodeNeighbours = 0

     count = 0
     DO i=Vel1Matrix % Rows(Vel1Perm(NodeNumber)*DOFs),&
          (Vel1Matrix % Rows((Vel1Perm(NodeNumber)*DOFs)+1)-1)
        IF(MOD(i,DOFs) .NE. 0) CYCLE !Stored DOF1, DOF2, DOF3, only need every 3rd
        IF(Vel1InvPerm(Vel1Matrix % Cols(i)) == NodeNumber) CYCLE !Not our own neighbour
        count = count + 1
        ThisNodeNeighbours(count) = &
             Vel1InvPerm(Vel1Matrix % Cols(i))
     END DO

   END SUBROUTINE FindNodeNeighbours

   !After finding groups, pass them between partitions and join...
   !Can just be done in one partition and then passed back?
   SUBROUTINE FindCrevasseGroups(CrevasseGroups, CheckValidity, MaskPerm)

     TYPE(CrevasseGroups_t), INTENT(INOUT) :: CrevasseGroups
     INTEGER :: MaxGroups = 100, MaxNodesPerGroup = 10000, locali
     LOGICAL :: CheckValidity, FoundIt
     INTEGER, POINTER, OPTIONAL, INTENT(IN) :: MaskPerm(:)

     ALLOCATE( &
          CrevasseGroups % NodeIndexes(MaxGroups,MaxNodesPerGroup), &
          CrevasseGroups % NotEmpty(MaxGroups), &
          CrevasseGroups % NoNodes(MaxGroups), &
          CrevasseGroups % Valid(MaxGroups))

     CrevasseGroups % NodeIndexes = 0
     CrevasseGroups % NoNodes = 0
     CrevasseGroups % NotEmpty = .FALSE.
     CrevasseGroups % Valid = .FALSE.
     CrevasseGroups % CurrentGroup = 0

     DO i=1,NoNodes
        IF(DistanceValues(DistancePerm(i))>FreeConnect) CYCLE
        IF(CIndexValues(CIndexPerm(i))<0) CYCLE
        IF(FrontPerm(i) > 0) CYCLE

        !Check if node is already in a group
        !This aint ideal, but I put it in to prevent a bug whereby, when surface and
        !basal groups join, all surface group nodes are ALSO basal group nodes, and so
        !'calving' is detected everywhere...
        FoundIt = .FALSE.
        DO locali = 1, SurfaceCrevasseGroups % CurrentGroup
          IF(ANY(SurfaceCrevasseGroups % NodeIndexes( &
               locali,1:SurfaceCrevasseGroups % NoNodes(locali)) .EQ. i)) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        !Only check basal groups if they exist already...
        IF(ALLOCATED(BasalCrevasseGroups % NodeIndexes)) THEN
          DO locali = 1, BasalCrevasseGroups % CurrentGroup
            IF(ANY(BasalCrevasseGroups % NodeIndexes( &
                 locali,1:BasalCrevasseGroups % NoNodes(locali)) .EQ. i)) THEN
              FoundIt = .TRUE.
              EXIT
            END IF
          END DO
        END IF
        IF(FoundIt) CYCLE

        !This point is only reached by nodes which ARE crevassing
        !and not already in a group, so it takes the first slot of
        !active group
        CrevasseGroups % CurrentGroup = CrevasseGroups % CurrentGroup + 1
        CrevasseGroups % NodeIndexes(CrevasseGroups % CurrentGroup,1) = i
        CrevasseGroups % NotEmpty(CrevasseGroups % CurrentGroup) = .TRUE.
        CrevasseGroups % NoNodes(CrevasseGroups % CurrentGroup) = &
             CrevasseGroups % NoNodes(CrevasseGroups % CurrentGroup) + 1
        CrevasseGroupValues(CrevasseGroupPerm(i)) = CrevasseGroups % CurrentGroup
        NextSlot = 2
        CALL SearchNeighbours(i, CrevasseGroups)
        NextSlot = 1
     END DO

     !Now we have all the groups.  If requested, check they are valid.
     IF(CheckValidity) THEN
        DO i=1,CrevasseGroups % CurrentGroup !Cycle all groups
           DO j=1,SIZE(CrevasseGroups % NodeIndexes,2) !Cycle all nodes in group
              IF(CrevasseGroups % NodeIndexes(i,j)==0) EXIT !End of group
              IF(MaskPerm(CrevasseGroups % NodeIndexes(i,j))>0) THEN
                 !At least 1 node in group is on relevant surface
                 CrevasseGroups % Valid(i) = .TRUE.
                 EXIT
              END IF
           END DO
        END DO
     END IF

   END SUBROUTINE FindCrevasseGroups

   !Recursively looks in neighbouring nodes to construct contiguous groups of
   !crevassing nodes.
   RECURSIVE SUBROUTINE SearchNeighbours(nodenum, CrevasseGroups)
     INTEGER :: nodenum, neighbourindex
     INTEGER :: localj,locali
     TYPE(CrevasseGroups_t) :: CrevasseGroups
     LOGICAL :: FoundIt

     IF(Debug) PRINT *,'Debug, nextslot, nodenum: ', nextslot, nodenum
     DO localj = 1,NumNeighbours(nodenum)
        neighbourindex = NodeNeighbours(nodenum,localj)
        IF(DistanceValues(DistancePerm(neighbourindex))>FreeConnect) CYCLE
        IF(CIndexValues(CIndexPerm(neighbourindex))<0) CYCLE
        IF(FrontPerm(localj) > 0) CYCLE

        !Check if node is already in a group
        !NOTE: is this actually possible?
        FoundIt = .FALSE.
        DO locali = 1, SurfaceCrevasseGroups % CurrentGroup
          IF(ANY(SurfaceCrevasseGroups % NodeIndexes(&
               locali,1:SurfaceCrevasseGroups % NoNodes(locali)) &
               .EQ. neighbourindex)) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF(ALLOCATED(BasalCrevasseGroups % NodeIndexes)) THEN
          DO locali = 1, BasalCrevasseGroups % CurrentGroup
            IF(ANY(BasalCrevasseGroups % NodeIndexes(&
                 locali,1:BasalCrevasseGroups % NoNodes(locali)) &
                 .EQ. neighbourindex)) THEN
              FoundIt = .TRUE.
              EXIT
            END IF
          END DO
        END IF
        IF(FoundIt) CYCLE

        !Only get here if node is valid new calving node, not already in a group
        CrevasseGroups % NodeIndexes(CrevasseGroups % CurrentGroup, NextSlot) = neighbourindex
        CrevasseGroups % NoNodes(CrevasseGroups % CurrentGroup) = &
             CrevasseGroups % NoNodes(CrevasseGroups % CurrentGroup) + 1
        CrevasseGroupValues(CrevasseGroupPerm(neighbourindex)) = CrevasseGroups % CurrentGroup
        NextSlot = NextSlot + 1
        IF(NextSlot > SIZE(CrevasseGroups % NodeIndexes, 2)) CALL Fatal(SolverName, &
             "More than 10000 nodes in a crevasse group? Almost certainly an error in setup")

        CALL SearchNeighbours(neighbourindex,CrevasseGroups)
        !Add a check for group validity (i.e. at least one node is on the relevant boundary
     END DO
   END SUBROUTINE SearchNeighbours

   SUBROUTINE FindCalvingBasal(SurfaceCrevasseGroups, BasalCrevasseGroups, BasalCalvingCoordinate)
     TYPE(CrevasseGroups_t), INTENT(IN) :: SurfaceCrevasseGroups, BasalCrevasseGroups
     REAL (KIND=dp), INTENT(OUT) :: BasalCalvingCoordinate
     INTEGER :: i,j,k,m, BasalNode, SurfaceNode

     BasalCalvingOccurs = .FALSE.
     BasalCalvingCoordinate = HUGE(BasalCalvingCoordinate)
     !Cycle surface crevasse groups
     DO i = 1, COUNT(SurfaceCrevasseGroups % NotEmpty)
        IF(.NOT.(SurfaceCrevasseGroups % Valid(i))) CYCLE

        !Cycle basal crevasse groups
        DO j = 1, COUNT(BasalCrevasseGroups % NotEmpty)
           IF(.NOT.(BasalCrevasseGroups % Valid(j))) CYCLE

           !Cycle Nodes in Surface Crevasse Group
           DO k = 1,SIZE(SurfaceCrevasseGroups % NodeIndexes, 2)
              IF(SurfaceCrevasseGroups % NodeIndexes(i,k)==0) EXIT
              SurfaceNode = SurfaceCrevasseGroups % NodeIndexes(i,k)

              !Cycle Nodes in Basal Crevasse Group
              DO m = 1,SIZE(BasalCrevasseGroups % NodeIndexes, 2)
                 IF(BasalCrevasseGroups % NodeIndexes(j,m)==0) EXIT
                 BasalNode = BasalCrevasseGroups % NodeIndexes(j,m)
                 !Is it the same node? i.e. do the groups touch?
                 IF(SurfaceNode == BasalNode) THEN
                    !Yes, so calving coord
                    IF(Mesh % Nodes % x(SurfaceNode) .LT. BasalCalvingCoordinate) THEN
                       BasalCalvingCoordinate = Mesh % Nodes % x(SurfaceNode)
                       BasalCalvingOccurs = .TRUE.
                       !TEST
                       IF(Debug) THEN
                          PRINT *, "Debugging Calving, Surface node x: ", &
                               Mesh % Nodes % x(SurfaceNode)," y: ",&
                               Mesh % Nodes % y(SurfaceNode)
                          PRINT *, "Debugging Calving, Basal node x: ", &
                               Mesh % Nodes % x(BasalNode)," y: ",&
                               Mesh % Nodes % y(BasalNode)
                       END IF
                    END IF
                 ELSE
                    !No, but are they neighbours?
                    IF(ANY(NodeNeighbours(BasalNode,:)==SurfaceNode)) THEN
                       !yes, neighbours
                       CALL CheckCIndexOverlap(SurfaceNode, BasalNode, OverlapCalvingCoordinate, OverlapOccurs)
                       IF(OverlapOccurs .AND. (OverlapCalvingCoordinate .LT. BasalCalvingCoordinate)) THEN
                          BasalCalvingCoordinate = OverlapCalvingCoordinate
                          BasalCalvingOccurs = .TRUE.
                       END IF
                    END IF
                    !not neighbours, cycle
                 END IF

              END DO
           END DO
        END DO
     END DO
   END SUBROUTINE FindCalvingBasal

   SUBROUTINE CheckCIndexOverlap(SurfaceNode, BasalNode, OverlapCoord, OverlapOccurs)

     IMPLICIT NONE
     REAL(KIND=dp) :: CSurfSurf, CSurfBasal, CBasalSurf, &
          CBasalBasal, XSurf, YSurf, XBasal, YBasal, dx, dy,dxdy, dCSurf, &
          dCBasal, xzerobasal, yzerobasal, xzerosurf, yzerosurf, &
          dbindexdy,dsindexdy
     REAL(KIND=dp), INTENT(OUT) :: OverlapCoord
     INTEGER :: SurfaceNode, BasalNode
     LOGICAL, INTENT(OUT) :: OverlapOccurs

     OverlapOccurs = .FALSE.


     !4 values at 2 nodes
     !Format is CIndexLocation.  So CSurfBasal is the value of
     !CSurfIndex at the Basal node.

     !This used to be simple linear interpolation between two nodes, but this was problematic:
     !For surface crevassing nodes, basal crevasse index is also positive, because they are
     !governed by the same equations.  This lead to situations where both CBasalSurf and CBasalBasal
     !would be small positive, with an OBVIOUS massive gap inbetween, but because both were positive,
     !linear interp didn't see the gap.

     !Strategy to overcome the previous:
     !No problem with CSurf* because it decreases with depth as expected.
     !Need to address CBasal*: at the lower node (C*Basal), assuming negligible upward changes in
     !stress_xx (any better way?), then the y-gradient of CBasal* is predictable, as the remaining
     !components decrease linearly with depth:
     ! rho_w * g * dy   -   BTuningParam * rho_i * g * dy

     !So, we replace the part of this subroutine which used to calculate the zero level of basal
     !crevassing.  We need:
     !g, rho_i, rho_w, BTuningParam

     !the various indices
     CSurfSurf = CSurfIndexValues(CSurfIndexPerm(SurfaceNode))
     CBasalSurf = CBasalIndexValues(CBasalIndexPerm(SurfaceNode))
     CSurfBasal = CSurfIndexValues(CSurfIndexPerm(BasalNode))
     CBasalBasal = CBasalIndexValues(CBasalIndexPerm(BasalNode))

     !the coords
     XSurf = Mesh % Nodes % x(SurfaceNode)
     YSurf = Mesh % Nodes % y(SurfaceNode)
     XBasal = Mesh % Nodes % x(BasalNode)
     YBasal = Mesh % Nodes % y(BasalNode)

     !The rest of the maths assumes that the surface node is ABOVE the basal node, so check:
     !If the node in the basal crevasse field is *above* the node in the surface crevasse field
     !we have calving and we assume its x-coord is between the two.  else check...

     IF(YSurf .LE. YBasal) THEN
        IF(Debug) PRINT *, 'DEBUG Calving: Basal Node above Surf Node, calving...'
        OverlapOccurs = .TRUE.
        OverlapCoord = (XSurf + XBasal)/2.0
     ELSE
        !Gradients
        dx = XSurf - XBasal
        dy = YSurf - YBasal !+ve
        dxdy = dx/dy

        dCSurf = CSurfSurf - CSurfBasal  !+ve
        dCBasal = CBasalSurf - CBasalBasal  !-ve

        dsindexdy = dCSurf/dy !+ve
        dbindexdy = -g*(RhoWF - Rho) !-ve

        yzerobasal = YBasal - (CBasalBasal/dbindexdy)
        yzerosurf = YSurf - (CSurfSurf/dsindexdy)

        IF(yzerosurf .LT. yzerobasal) THEN
           OverlapOccurs = .TRUE.

           xzerobasal = xbasal + ((yzerobasal - YBasal)*dxdy)
           xzerosurf = xsurf + ((yzerosurf - YSurf)*dxdy)

           OverlapCoord = MIN(xzerosurf, xzerobasal)
           OverlapCoord = MAX(OverlapCoord,MIN(XSurf, XBasal))
        ELSE
           OverlapOccurs = .FALSE.
        ENDIF
     END IF

     IF(OverlapOccurs .AND. Debug) THEN
        PRINT *, "Overlap occurs!, OverlapCoord: ",OverlapCoord
        PRINT *, "Surf: ",SurfaceNode," x: ",Mesh % Nodes % x(SurfaceNode),&
             "y: ",Mesh % Nodes % y(SurfaceNode)
        PRINT *, "Base: ",BasalNode," x: ",Mesh % Nodes % x(BasalNode),&
             "y: ",Mesh % Nodes % y(BasalNode)

     END IF

   END SUBROUTINE CheckCIndexOverlap
END SUBROUTINE Find_Calving
