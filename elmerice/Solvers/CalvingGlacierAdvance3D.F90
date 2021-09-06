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
        ShiftIdx, ShiftToIdx, NoTangledGroups, PivotIdx, CornerIdx, &
        SecondIdx, FirstTangleIdx, LastTangleIdx, Nl, Nr, Naux, ok, Nrail
   INTEGER, POINTER :: Perm(:), FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        FrontNodeNums(:)=>NULL(),LeftPerm(:)=>NULL(), RightPerm(:)=>NULL()
   INTEGER, ALLOCATABLE :: FrontLocalNodeNumbers(:), &
        NodeNumbers(:), TangledGroup(:), TangledPivotIdx(:), UpdatedDirection(:)
   REAL(KIND=dp) :: NodeVelo(3), NodeMelt(3), NodeNormal(3), RailDir(2),&
        MeltRate, Displace(3), y_coord(2), epsShift, LongRangeLimit, MaxDisplacement, &
        EpsTangle,thisEps,Shift, thisY,xx,yy,TempDist,MinDist,xt,yt,t
   REAL(KIND=dp), POINTER :: Advance(:)
   REAL(KIND=dp), ALLOCATABLE :: Rot_y_coords(:,:), Rot_z_coords(:,:), &
        TangledShiftTo(:), xL(:),yL(:),xR(:),yR(:),xRail(:),yRail(:)
   LOGICAL :: Found, Debug, Parallel, Boss, ShiftLeft, LeftToRight, MovedOne, ShiftSecond, &
        Protrusion, SqueezeLeft, SqueezeRight, FirstTime=.TRUE., intersect_flag, FrontMelting, &
        MovePastRailNode
   LOGICAL, ALLOCATABLE :: DangerZone(:), WorkLogical(:), &
        Tangled(:), DontMove(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VeloVarName, MeltVarName, &
        NormalVarName, FrontMaskName, TopMaskName, TangledVarName, &
        LeftMaskName, RightMaskName, LeftRailFName, RightRailFName
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

   TangledVarName = "Tangled"
   TangledVar => VariableGet(Mesh % Variables, TangledVarName, .TRUE., UnfoundFatal=.TRUE.)
   TangledVar % Values = 0.0_dp

   !Get front normal vector
   NormalVarName = ListGetString(Params, "Normal Vector Variable Name", UnfoundFatal=.TRUE.)
   NormalVar => VariableGet(Mesh % Variables, NormalVarName, .TRUE., UnfoundFatal=.TRUE.)

   MaxDisplacement = ListGetConstReal(Params, "Maximum Node Displacement", Found)
   IF(.NOT.Found) THEN
     CALL Info(SolverName, "'Maximum Node Displacement' not found, setting to 1.0E4.")
     MaxDisplacement = 1.0E4_dp
   END IF

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
   
   !--------------------------------------
   ! Action: Compute lagrangian displacement for all nodes
   !         This is the main function of the Solver.
   !         TO DO: check for 'tangled' regions, see CalvingFrontAdvance3D.F90
   !--------------------------------------
    
   DO i=1,Mesh % NumberOfNodes
      IF(Perm(i) <= 0) CYCLE

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
            
      Displace = 0.0
      IF ( FrontPerm(i) > 0 ) THEN
         Displace(1) =  NodeVelo(1) - NodeMelt(1)
         Displace(2) =  NodeVelo(2) - NodeMelt(2)
         Displace(3) =  NodeVelo(3) - NodeMelt(3)
	 Displace = Displace * dt
      END IF
      ! NOTE: Displace overwritten on corner of left-front and right-front
      ! melt not taken into account on those corner nodes
      IF ((LeftPerm(i)>0) .OR. (RightPerm(i)>0)) THEN  
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
         xt=xx+NodeVelo(1)*dt
         yt=yy+NodeVelo(2)*dt
         IF (LeftPerm(i)>0) THEN
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
         IF(MovePastRailNode) THEN
            ! k is defined such that the node should after displacement lay on jmin -- k segment
            ! nodes first move along the current rail segment to the nearest rail node jmin.
            ! whatever magnitude of the horizontal part of v*dt is left
            ! will be the distance the node travels on the new segment.
            !RailDir=(/xRail(jmin),yRail(jmin)/)-(/xx,yy/) ! direction of current rail
            !Displace(1:2)=(/xRail(jmin),yRail(jmin)/)-(/xx,yy/) ! get to new rail node
            !TempDist is total distance it would have travelled only along its original segment
            !TempDist=DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
            !     (RailDir(1)**2+RailDir(2)**2)**0.5 * dt / DOT_PRODUCT(RailDir,RailDir)
            ! t is the proportion left to travel along new segment
            !t=(TempDist - ((xx-xRail(jmin))**2.+(yy-yRail(jmin))**2.)**0.5)/TempDist
            ! new rail direction
            RailDir=(/xRail(jmin-1),yRail(jmin-1)/)-(/xRail(jmin+1),yRail(jmin+1)/)
            ! add proportion left to travel along new direction
              Displace(1:2)+DOT_PRODUCT(NodeVelo(1:2),RailDir)
              RailDir * dt * t / DOT_PRODUCT(RailDir,RailDir)
            Displace(1:2) = DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
                 RailDir * dt / DOT_PRODUCT(RailDir,RailDir)
         ELSE ! not moving past node, just project onto current rail
            RailDir=(/xRail(jmin),yRail(jmin)/)-(/xRail(k),yRail(k)/)
            Displace(1:2) = DOT_PRODUCT(NodeVelo(1:2),RailDir) * &
             RailDir / DOT_PRODUCT(RailDir,RailDir)
            Displace = Displace * dt
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


   !---------------------------------------
   !Done, just deallocations

   FirstTime = .FALSE.

  DEALLOCATE(FrontPerm, TopPerm, LeftPerm, RightPerm)

 END SUBROUTINE GlacierAdvance3D
