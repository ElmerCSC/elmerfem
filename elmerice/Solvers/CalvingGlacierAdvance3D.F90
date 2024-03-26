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
! *  Authors: Eef van Dongen, Joe Todd
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *  
! *
! *****************************************************************************
!
! A routine for getting frontal advance in 3D, given velocity and melt
! nodes on the intersection of Front and LateralMargin are also allowed to advance
! In order to assure the glacier still aligns with the valley walls,
! 'rails' are prescribed. The points on the rails need to be given in two 
! ASCII files (x y), one for left and one for right.
! x and y values need to be sorted, but the direction is not of importance. (?)
! NOTE: if the prescribed timestep is too large and 'rails' are very nonsmooth
! this implementation may not obey mass conservation.
 SUBROUTINE GlacierAdvance3D ( Model, Solver, dt, TransientSimulation )

   USE CalvingGeometry
   USE DefUtils
   IMPLICIT NONE

!-----------------------------------------------
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation
!-----------------------------------------------
   TYPE(Mesh_t), POINTER :: Mesh
   TYPE(Nodes_t), TARGET :: FrontNodes
   TYPE(ValueList_t), POINTER :: Params
   TYPE(Variable_t), POINTER :: Var, VeloVar, MeltVar, NormalVar, TangledVar, &
        DistVar
! TO DO clean all unused variables
   INTEGER :: i, j,jmin, k, m, n, DOFs, TotalNodes,&
        FaceNodeCount, DummyInt,  hits, ierr, FrontLineCount, county,&
        ShiftIdx, ShiftToIdx, PivotIdx, CornerIdx, &
        SecondIdx, FirstTangleIdx, LastTangleIdx, Nl, Nr, Naux, ok, Nrail,&
        NNodes, NBulk, NBdry, counter, LeftBCtag, RightBCtag, FrontBCtag,&
        OnRails
   INTEGER, POINTER :: Perm(:), FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        FrontNodeNums(:)=>NULL(),LeftPerm(:)=>NULL(), RightPerm(:)=>NULL(), &
        NodeIndexes(:),InflowPerm(:)=>NULL()
   INTEGER, ALLOCATABLE :: FrontLocalNodeNumbers(:), &
        NodeNumbers(:), UpdatedDirection(:)
   REAL(KIND=dp) :: NodeVelo(3), NodeMelt(3), NodeNormal(3), RailDir(2),&
        MeltRate, Displace(3), y_coord(2), epsShift, LongRangeLimit, MaxDisplacement, &
        EpsTangle,thisEps,Shift, thisY,xx,yy,TempDist,MinDist,xt,yt,t, &
        a1(2), a2(2), b1(2), b2(2), b3(2), intersect(2), DistRailNode, RDisplace(3),&
        buffer, VeloFactor
   REAL(KIND=dp), POINTER :: Advance(:)
   REAL(KIND=dp), ALLOCATABLE :: Rot_y_coords(:,:), Rot_z_coords(:,:), &
        xL(:),yL(:),xR(:),yR(:),xRail(:),yRail(:)
   LOGICAL :: Found, Debug, Parallel, Boss, ShiftLeft, LeftToRight, MovedOne, ShiftSecond, &
        Protrusion, SqueezeLeft, SqueezeRight, FirstTime=.TRUE., intersect_flag, FrontMelting, &
        MovePastRailNode, HitsRails, Reverse, ThisBC, MoveBulk
   LOGICAL, ALLOCATABLE :: DangerZone(:), WorkLogical(:), &
        DontMove(:), FrontToRight(:), FrontToLeft(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VeloVarName, MeltVarName, &
        NormalVarName, FrontMaskName, TopMaskName, &
        LeftMaskName, RightMaskName, LeftRailFName, RightRailFName,&
        InflowMaskName
   INTEGER,parameter :: io=20  
   SAVE :: FirstTime ! TO DO not actually using FirstTime here?

   !-----------------------------------------------
   ! Initialisation
   !-----------------------------------------------

   Debug = .FALSE.
   Parallel = (ParEnv % PEs > 1)
   Boss = ParEnv % MyPE == 0

   SolverName = "GlacierAdvance3D"
   Params => Solver % Values 
   Mesh => Solver % Mesh

   !The main solver var contains the magnitude of front advance
   Var => Solver % Variable
   Advance => Var % Values
   Perm => Var % Perm
   DOFs = Var % DOFs
   IF(Var % DOFs /= 3) CALL Fatal(SolverName, "Variable should have 3 DOFs...")

   !Get the flow solution
   VeloVarName = ListGetString(Params, "Flow Solution Variable Name", Found)
   IF(.NOT. Found) THEN
     CALL Info(SolverName, "Flow Solution Variable Name not found, assuming 'Flow Solution'")
     VeloVarName = "Flow Solution"
   END IF
   VeloVar => VariableGet(Mesh % Variables, VeloVarName, .TRUE., UnfoundFatal=.TRUE.)
   
   LeftRailFName = ListGetString(Params, "Left Rail File Name", Found)
   IF(.NOT. Found) THEN
      CALL Info(SolverName, "Left Rail File Name not found, assuming './LeftRail.xy'")
      LeftRailFName = "LeftRail.xy"
   END IF
   Nl = ListGetInteger(Params, "Left Rail Number Nodes", Found)
   IF(.NOT.Found) THEN
      WRITE(Message,'(A,A)') 'Left Rail Number Nodes not found'
      CALL FATAL(SolverName, Message)
   END IF
   !TO DO only do these things if firsttime=true?
   OPEN(unit = io, file = TRIM(LeftRailFName), status = 'old',iostat = ok)
   print *, ok
   IF (ok /= 0) THEN
      WRITE(message,'(A,A)') 'Unable to open file ',TRIM(LeftRailFName)
      CALL FATAL(Trim(SolverName),Trim(message))
   END IF
   ALLOCATE(xL(Nl), yL(Nl))

   ! read data
   DO i = 1, Nl
      READ(io,*,iostat = ok, end=200) xL(i), yL(i)
   END DO
200 Naux = Nl - i
   IF (Naux > 0) THEN
      WRITE(message,'(I0,A,I0,A,A)') Naux,' out of ',Nl,' datasets in file ', TRIM(LeftRailFName)
      CALL INFO(Trim(SolverName),Trim(message))
   END IF
   CLOSE(io)
   RightRailFName = ListGetString(Params, "Right Rail File Name", Found)
   IF(.NOT. Found) THEN
      CALL Info(SolverName, "Right Rail File Name not found, assuming './RightRail.xy'")
      RightRailFName = "RightRail.xy"
   END IF

   Nr = ListGetInteger(Params, "Right Rail Number Nodes", Found)
   IF(.NOT.Found) THEN
      WRITE(Message,'(A,A)') 'Right Rail Number Nodes not found'
      CALL FATAL(SolverName, Message)
   END IF
   !TO DO only do these things if firsttime=true?
   OPEN(unit = io, file = TRIM(RightRailFName), status = 'old',iostat = ok)

   IF (ok /= 0) THEN
      WRITE(message,'(A,A)') 'Unable to open file ',TRIM(RightRailFName)
      CALL FATAL(Trim(SolverName),Trim(message))
   END IF
   ALLOCATE(xR(Nr), yR(Nr))

   ! TO DO would be nice if solver warns in case there is more data in file than number Nr/Nl tells 
   ! read data
   DO i = 1, Nr
      READ(io,*,iostat = ok, end=100) xR(i), yR(i)
   END DO
100 Naux = Nr - i
   IF (Naux > 0) THEN
      WRITE(message,'(I0,A,I0,A,A)') Naux,' out of ',Nl,' datasets in file ', TRIM(RightRailFName)
      CALL INFO(Trim(SolverName),Trim(message))
   END IF
   CLOSE(io)

   !Get melt rate
   MeltVarName = ListGetString(Params, "Melt Variable Name", Found)
   IF(.NOT. Found) THEN
     CALL Info(SolverName, "Melt Variable Name not found, assuming no frontal melting")
     FrontMelting = .FALSE.
   ELSE
     FrontMelting = .TRUE.
     MeltVar => VariableGet(Mesh % Variables, MeltVarName, .TRUE., UnfoundFatal=.TRUE.)
   END IF

   !Get front normal vector
   NormalVarName = ListGetString(Params, "Normal Vector Variable Name", UnfoundFatal=.TRUE.)
   NormalVar => VariableGet(Mesh % Variables, NormalVarName, .TRUE., UnfoundFatal=.TRUE.)

   MaxDisplacement = ListGetConstReal(Params, "Maximum Node Displacement", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Maximum Node Displacement' not found, setting to 1.0E4.")
     MaxDisplacement = 1.0E4_dp
   END IF

   VeloFactor = ListGetConstReal(Params, "Lateral Margin Velocity Factor", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Lateral Margin Velocity Factor' not found, setting to 1.0.")
     VeloFactor = 1.0_dp
   END IF

   buffer = ListGetConstReal(Params, "Rail Buffer", Found, DefValue=0.1_dp)
   IF(.NOT. Found) CALL Info(SolverName, "No Rail Buffer set using default 0.1")

   MoveBulk = ListGetLogical(Params,"MoveBulk", Found, DefValue=.FALSE.)
   IF(.NOT. Found) CALL Info(SolverName, "Not moving bulk as default")

   !Get the front line
   FrontMaskName = "Calving Front Mask"
   TopMaskName = "Top Surface Mask"

   CALL MakePermUsingMask( Model, Solver, Mesh, TopMaskName, &
        .FALSE., TopPerm, dummyint)

   CALL MakePermUsingMask( Model, Solver, Mesh, FrontMaskName, &
        .FALSE., FrontPerm, FaceNodeCount)

   LeftMaskName = "Left Sidewall Mask"
   RightMaskName = "Right Sidewall Mask"

   !Generate perms to quickly get nodes on each boundary
   CALL MakePermUsingMask( Model, Solver, Mesh, LeftMaskName, &
        .FALSE., LeftPerm, dummyint)
   CALL MakePermUsingMask( Model, Solver, Mesh, RightMaskName, &
        .FALSE., RightPerm, dummyint)

   InflowMaskName = "Inflow Mask"
   CALL MakePermUsingMask( Model, Solver, Mesh, InflowMaskName, &
        .FALSE., InflowPerm, dummyint)
   
   !--------------------------------------
   ! Action: Compute lagrangian displacement for all nodes
   !         This is the main function of the Solver.
   !         TO DO: check for 'tangled' regions, see CalvingFrontAdvance3D.F90
   !--------------------------------------

   NNodes = Mesh % NumberOfNodes
   NBulk = Mesh % NumberOfBulkElements
   NBdry = Mesh % NumberOfBoundaryElements

   ALLOCATE(FrontToRight(NNodes), FrontToLeft(NNodes))
   FrontToRight = .FALSE.; FrontToLeft = .FALSE.

   Advance = 0.0_dp
   IF(MoveBulk) THEN ! for a fully Lagrangian mesh
      DO i=1, Mesh % NumberOfNodes
        IF(InflowPerm(i) > 0) CYCLE
        IF(FrontPerm(i) == 0 .AND. LeftPerm(i) == 0 .AND. RightPerm(i) == 0) THEN
          NodeVelo(1) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 1)
          NodeVelo(2) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 2)
          NodeVelo(3) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 3)

          Displace(1) =  NodeVelo(1)
          Displace(2) =  NodeVelo(2)
          Displace(3) =  NodeVelo(3)
          Displace = Displace * dt

          Advance((Perm(i)-1)*DOFs + 1) = Displace(1)
          Advance((Perm(i)-1)*DOFs + 2) = Displace(2)
          Advance((Perm(i)-1)*DOFs + 3) = Displace(3)
        END IF
      END DO
   END IF
    
   DO i=1,Mesh % NumberOfNodes
      IF(Perm(i) <= 0) CYCLE
      IF(FrontPerm(i) == 0 .AND. LeftPerm(i) == 0 .AND. RightPerm(i) == 0) CYCLE
      IF(InflowPerm(i) > 0) CYCLE

      IF(FrontMelting .AND. (FrontPerm(i) >0 )) THEN
         IF(MeltVar % Perm(i) <= 0) &
              CALL Fatal(SolverName, "Permutation error on front node!")


         !Scalar melt value from Plume solver
         MeltRate = MeltVar % Values(MeltVar % Perm(i))

         NodeNormal(1) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 1)
         NodeNormal(2) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 2)
         NodeNormal(3) = NormalVar % Values(((NormalVar % Perm(i)-1)*NormalVar % DOFs) + 3)

         NodeMelt = NodeNormal * MeltRate
      ELSE
         NodeMelt = 0.0_dp
      END IF

      !Compute front normal component of velocity
      NodeVelo(1) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 1)
      NodeVelo(2) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 2)
      NodeVelo(3) = VeloVar % Values(((VeloVar % Perm(i)-1)*VeloVar % DOFs) + 3)

      IF(FrontPerm(i) > 0) THEN
         IF(LeftPerm(i) > 0 .OR. RightPerm(i) > 0) THEN
            NodeVelo = NodeVelo * VeloFactor
         END IF
      END IF
            
      Displace = 0.0
      IF ( FrontPerm(i) > 0 ) THEN
         Displace(1) =  NodeVelo(1) - NodeMelt(1)
         Displace(2) =  NodeVelo(2) - NodeMelt(2)
         Displace(3) =  NodeVelo(3) - NodeMelt(3)
	      Displace = Displace * dt
      END IF

      OnRails = 0
      ! if not lateral node need to check we don't cross rails if fjord narrowing
      IF(LeftPerm(i) == 0 .AND. RightPerm(i) == 0) THEN
         a1(1) = Solver % Mesh % Nodes % x(i)
         a1(2) = Solver % Mesh % Nodes % y(i)
         a2(1) = a1(1) + Displace(1)
         a2(2) = a1(2) + Displace(2)

         DO j=1, Nl-1
            b1(1) = xL(j); b1(2) = yL(j)
            b2(1) = xL(j+1); b2(2) = yL(j+1)
            IF(j < NL - 1) THEN
               b3(1) = xL(j+2); b3(2) = yL(j+2)
            END IF
            tempdist = PointLineSegmDist2D(b1, b2, a1)
            IF(tempdist < buffer + AEPS) THEN
               FrontToLeft(i) = .TRUE.
               OnRails = 1 ! on rails so treat as edge
               EXIT
            END IF

            CALL LineSegmentsIntersect(a1, a2, b1, b2, intersect, HitsRails)

            IF(HitsRails) THEN
               FrontToLeft(i) = .TRUE.
               ! check direction
               Reverse = .FALSE.
               IF (PointDist2D(intersect,b2) > PointDist2D(a1, b2)) Reverse = .TRUE.

               IF(Reverse) THEN
                  b1(1) = xL(j+1); b1(2) = yL(j+1)
                  b2(1) = xL(j); b2(2) = yL(j)
                  IF(j > 1) THEN
                     b3(1) = xL(j-1); b3(2) = yL(j-1)
                  END IF
               END IF
               EXIT
            END IF
         END DO

         IF(.NOT. HitsRails) THEN ! didn't hit left rail
            DO j=1, Nr-1
               b1(1) = xR(j); b1(2) = yR(j)
               b2(1) = xR(j+1); b2(2) = yR(j+1)
               IF(j < Nr - 1) THEN
                  b3(1) = xR(j+2); b3(2) = yR(j+2)
               END IF
               tempdist = PointLineSegmDist2D(b1, b2, a1)
               IF(tempdist < buffer + AEPS) THEN
                  FrontToRight(i) = .TRUE.
                  OnRails = 2! on rails so treat as edge
                  EXIT
               END IF

               CALL LineSegmentsIntersect(a1, a2, b1, b2, intersect, HitsRails)

               IF(HitsRails) THEN
                  FrontToRight(i) = .TRUE.
                  ! check direction
                  Reverse = .FALSE.
                  IF (PointDist2D(intersect,b2) > PointDist2D(a1, b2)) Reverse = .TRUE.

                  IF(Reverse) THEN
                     b1(1) = xR(j+1); b1(2) = yR(j+1)
                     b2(1) = xR(j); b2(2) = yR(j)
                     IF(j > 1) THEN
                        b3(1) = xR(j-1); b3(2) = yR(j-1)
                     END IF
                  END IF
                  EXIT
               END IF
            END DO
         END IF

         IF(HitsRails) THEN ! hit a rail b1, b2, and b3 should be correct from above
            ! limit to rail and then move along
            Displace(1:2) = intersect-a1

            !remaining displacement
            RDisplace = 0.0_dp
            RDisplace(1) =  (NodeVelo(1) - NodeMelt(1)) * dt - Displace(1)
            RDisplace(2) =  (NodeVelo(2) - NodeMelt(2)) * dt - Displace(2)

            DistRailNode = PointDist2D(intersect, b2)
            TempDist = ( RDisplace(1)**2 + RDisplace(2)**2 ) ** 0.5

            !check if move past rail node
            MovePastRailNode = .FALSE.
            IF(TempDist > DistRailNode) MovePastRailNode=.TRUE.

            !project along rail
            IF(.NOT. MovePastRailNode) THEN

               RailDir = b1 - b2
               Displace(1:2) = Displace(1:2) + DOT_PRODUCT(RDisplace(1:2),RailDir) * &
                  RailDir / DOT_PRODUCT(RailDir,RailDir)

            ELSE

               IF((Reverse .AND. j == 1) .AND. (.NOT. Reverse .AND. j == Nr - 1)) &
                  CALL FATAL(SolverName, 'Moving past the end of the rails')

               RailDir = b1 - b2
               Displace(1:2) = Displace(1:2) + b2 - a1 ! get to new rail node
               !TempDist is total distance it would have travelled only along its original segment
               TempDist=DOT_PRODUCT(RDisplace(1:2),RailDir) * &
                  (RailDir(1)**2+RailDir(2)**2)**0.5 * dt / DOT_PRODUCT(RailDir,RailDir)
               TempDist=ABS(TempDist) ! distance no sign

               IF(TempDist < DistRailNode) THEN
                  CALL WARN(SolverName, 'Node outside rails and moving past rail node')
                  TempDist = DistRailNode
               END IF

               t=(TempDist - DistRailNode)/TempDist

               ! new rail direction
               RailDir = b2 - b3
               ! add proportion left to travel along new direction
               Displace(1:2) = Displace(1:2)+DOT_PRODUCT(RDisplace(1:2),RailDir) * &
                  RailDir * dt * t / DOT_PRODUCT(RailDir,RailDir)
            END IF
         END IF
      END IF

      ! NOTE: Displace overwritten on corner of left-front and right-front
      ! melt not taken into account on those corner nodes
      IF ((LeftPerm(i)>0) .OR. (RightPerm(i)>0) .OR. (OnRails > 0) ) THEN
         ! Strategy to project displace from x to x+v*dt along rail direction:
         ! - find closest rail node, then find cosest rail segment.
         ! - three cases to find which rail segment should be used:
         ! 1) x is rail node, then just find out which segment closest to x+dt*v
         ! 2) both x and x+dt*v are on one rail segment
         ! 3) x is on one segment and x+dt*v on another
         ! only do x,y direction; vertical update is handled by freesurf.
         ! NOTE mass is not strictly conserved but for small timesteps and
         ! a thin lateral boundary artificial mass change should be insignificant.
         xx = Solver % Mesh % Nodes % x(i)
         yy = Solver % Mesh % Nodes % y(i)
         NodeVelo(1:2) = NodeVelo(1:2) - NodeMelt(1:2)
         xt=xx+(NodeVelo(1))*dt
         yt=yy+(NodeVelo(2))*dt
         IF (LeftPerm(i)>0 .OR. OnRails == 1) THEN
            Nrail= Nl
            ALLOCATE(xRail(Nrail), yRail(Nrail))
            xRail = xL
            yRail = yL
         ELSE
            Nrail= Nr
            ALLOCATE(xRail(Nrail), yRail(Nrail))
            xRail = xR
            yRail = yR ! TO DO use pointers instead?
         END IF
         MinDist=(xRail(1)-xRail(Nrail))**2.+(yRail(1)-yRail(Nrail))**2.
         ! MinDist is actually maximum distance, needed for finding closest rail node
         DO j=1,Nrail ! Find closest point on rail
            TempDist=(xRail(j)-xx)**2.+(yRail(j)-yy)**2.
            IF(TempDist < MinDist) THEN
               MinDist=TempDist
               jmin=j
            END IF
         END DO
         IF(jmin+1>Nrail) PRINT *, 'jmin=',jmin,'ERROR: advancing/retreating beyond rail!'
         !NOTE, this assumes rails are starting at back and going towards front, should check k > 0 as well?            
         IF (MinDist < AEPS) THEN ! min distance very close to 0
            ! case 1)  mesh node is a rail point
            ! this likely happens at t=0.
            MovePastRailNode=.FALSE. !if both at a rail node and passing next rail node, rail resolution is smaller than displacement per timestep
            ! we assume that's never the case. TO DO include assert
            ! check whether x+v*dt is closest to segment jmin+1 or jmin-1  
            IF(PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                 (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                 > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                 (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
               k=jmin-1 ! x+v*dt is on jmin-1 -- jmin
            ELSE
               k=jmin+1
            END IF
            PRINT *, 'NOTE! mesh node on a rail node!'
         ELSE ! x is not a rail node, check whether x is closest to  jmin+1 or jmin-1
            IF(jmin==1) THEN
              IF(PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
              (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
              > PointLineSegmDist2D((/xRail(jmin+1),yRail(jmin+1)/), &
              (/xRail(jmin+2),yRail(jmin+2)/),(/xt,yt/))) THEN
                ! x+v*dt also closest to jmin--jmin-1
                MovePastRailNode=.TRUE. ! case 2) x+dt*v and x are on same segment
                k=jmin+1
              ELSE
                ! case 3) advancing past rail segment
                MovePastRailNode=.FALSE.
                k=jmin+1 ! x is moving from jmin-1--jmin to jmin--jmin+1
              END IF
            ELSE
              IF(PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                (/xRail(jmin+1),yRail(jmin+1)/),(/xx,yy/)) &
                > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                (/xRail(jmin-1),yRail(jmin-1)/),(/xx,yy/))) THEN
                  ! x closest to jmin-1 -- jmin segment
                  IF (PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                      (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                      > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                      (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
                      ! x+v*dt also closest to jmin--jmin-1
                      MovePastRailNode=.FALSE. ! case 2) x+dt*v and x are on same segment
                      k=jmin-1 ! x and x+v*dt are closest to jmin>jmin-1 segment
                  ELSE
                      ! case 3) advancing past rail segment
                      MovePastRailNode=.TRUE.
                      k=jmin+1 ! x is moving from jmin-1--jmin to jmin--jmin+1
                  END IF
              ELSE ! x closest to jmin -- jmin+1 segment, but closer to jmin than jmin+1
                  ! assuming not moving past jmin+1 (see above assumption on time step and rail resolution)
                  IF (PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                      (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                      > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                      (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
                      ! Case 3) x advancing past rail segment
                      MovePastRailNode=.TRUE.
                      k=jmin-1 ! x is moving from jmin+1>jmin to jmin>jmin-1
                  ELSE
                      ! Case 2) x+v*dt and x are on same segment jmin>jmin+1
                      MovePastRailNode=.FALSE.
                      k=jmin+1
                  END IF
              END IF
            END IF
         END IF
         IF(MovePastRailNode) THEN
            ! check not moving past last node
            IF((k > jmin .AND. jmin == Nrail) .OR. (k < jmin .AND. jmin == 1)) &
               CALL FATAL(SolverName, 'Moving past the end of the rails')
            ! k is defined such that the node should after displacement lay on jmin -- k segment
            ! nodes first move along the current rail segment to the nearest rail node jmin.
            ! whatever magnitude of the horizontal part of v*dt is left
            ! will be the distance the node travels on the new segment.
            RailDir=(/xRail(jmin),yRail(jmin)/)-(/xx,yy/) ! direction of current rail
            Displace(1:2)=(/xRail(jmin),yRail(jmin)/)-(/xx,yy/) ! get to new rail node
            !TempDist is total distance it would have travelled only along its original segment
            TempDist=DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
                 (RailDir(1)**2+RailDir(2)**2)**0.5 * dt / DOT_PRODUCT(RailDir,RailDir)
            TempDist=ABS(TempDist) ! distance no sign
            ! t is the proportion left to travel along new segment
            DistRailNode = ((xx-xRail(jmin))**2.+(yy-yRail(jmin))**2.)**0.5
            IF(TempDist < DistRailNode) THEN
               CALL WARN(SolverName, 'Node outside rails and moving past rail node')
               TempDist = DistRailNode
            END IF

            t=(TempDist - DistRailNode)/TempDist

            ! new rail direction
            RailDir=(/xRail(jmin),yRail(jmin)/)-(/xRail(k),yRail(k)/)
            ! add proportion left to travel along new direction
            Displace(1:2) = Displace(1:2)+DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
              RailDir * dt * t / DOT_PRODUCT(RailDir,RailDir)
         ELSE ! not moving past node, just project onto current rail
            RailDir=(/xRail(jmin),yRail(jmin)/)-(/xRail(k),yRail(k)/)
            Displace(1:2) = DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
             RailDir / DOT_PRODUCT(RailDir,RailDir)
            Displace(1:2) = Displace(1:2) * dt
         END IF
         ! TO DO write warning if distance point on margin to rail segment too large?
         ! TO DO also write warning if glacier advances out of range of rails?
         ! TO DO ASSERT that ||NewDisplace|| <= ||V||*dt
         DEALLOCATE(xRail,yRail)
     END IF


     IF(MAXVAL(Displace) > MaxDisplacement) THEN
       WRITE(Message,'(A,i0,A)') "Maximum allowable front displacement exceeded for node ",i,". Limiting..."
       CALL Warn(SolverName, Message)
       Displace = Displace * (MaxDisplacement/MAXVAL(Displace))
     END IF

     Advance((Perm(i)-1)*DOFs + 1) = Displace(1)
     Advance((Perm(i)-1)*DOFs + 2) = Displace(2)
     Advance((Perm(i)-1)*DOFs + 3) = Displace(3)

   END DO

   ! find lateral boundary tags
   DO i=1,Model % NumberOfBCs
      ThisBC = ListGetLogical(Model % BCs(i) % Values,LeftMaskName,Found)
      IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
      LeftBCtag =  Model % BCs(i) % Tag
      EXIT
   END DO
   DO i=1,Model % NumberOfBCs
      ThisBC = ListGetLogical(Model % BCs(i) % Values,RightMaskName,Found)
      IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
      RightBCtag =  Model % BCs(i) % Tag
      EXIT
   END DO
   DO i=1,Model % NumberOfBCs
      ThisBC = ListGetLogical(Model % BCs(i) % Values,FrontMaskName,Found)
      IF((.NOT. Found) .OR. (.NOT. ThisBC)) CYCLE
      FrontBCtag =  Model % BCs(i) % Tag
      EXIT
   END DO

   !reassign any elements that now only contain nodes on lateral margins
   DO i=NBulk+1, NBulk+NBdry
      IF(Mesh % Elements(i) % BoundaryInfo % constraint /= FrontBCtag) CYCLE
      NodeIndexes => Mesh % Elements(i) % NodeIndexes
      n = Mesh % Elements(i) % TYPE % NumberOfNodes
      counter = 0
      DO j=1,n
         IF(RightPerm(NodeIndexes(j)) > 0) CYCLE
         IF(FrontToRight(NodeIndexes(j))) CYCLE
         counter = counter + 1
      END DO

      IF(counter == 0) Mesh % Elements(i) % BoundaryInfo % constraint = RightBCtag
   END DO

   !reassign any elements that now only contain nodes on lateral margins
   DO i=NBulk+1, NBulk+NBdry
      IF(Mesh % Elements(i) % BoundaryInfo % constraint /= FrontBCtag) CYCLE
      NodeIndexes => Mesh % Elements(i) % NodeIndexes
      n = Mesh % Elements(i) % TYPE % NumberOfNodes
      counter = 0
      DO j=1,n
         IF(LeftPerm(NodeIndexes(j)) > 0) CYCLE
         IF(FrontToLeft(NodeIndexes(j))) CYCLE
         counter = counter + 1
      END DO

      IF(counter == 0) Mesh % Elements(i) % BoundaryInfo % constraint = LeftBCtag
   END DO

   !---------------------------------------
   !Done, just deallocations

   FirstTime = .FALSE.

  DEALLOCATE(FrontPerm, TopPerm, LeftPerm, RightPerm, InflowPerm)

 END SUBROUTINE GlacierAdvance3D
