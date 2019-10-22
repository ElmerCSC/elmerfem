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
   REAL(KIND=dp) :: FrontOrientation(3),RotationMatrix(3,3),&
        NodeVelo(3), NodeMelt(3), NodeNormal(3), RailDir(2),&
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

   !Get the orientation of the calving front, compute rotation matrix
   !TODO: generalize and link
   FrontOrientation = GetFrontOrientation(Model)
   RotationMatrix = ComputeRotationMatrix(FrontOrientation)

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
   !         Everything following this DO loop is taking care of
   !         unprojectability/high gradient etc
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
         IF ((LeftPerm(i)>0) .OR. (RightPerm(i)>0)) THEN   ! TO DO on left or right corner
            ! Strategy to project displace from x to x+dx along rail direction:
            ! - find closest rail node, then
            ! - three cases to find which rail segment should be used:
            ! 1) x is rail node, then just find out which segment closest to x+dx
            ! 2) both x and x+dx are on one rail segment
            ! 3) x is on one segment and x+dx on another
            ! only do x,y direction; no penetration in vertical is handled by freesurf.
            ! NOTE mass is not conserved but for small timesteps and a relatively thin lateral boundary it should be insignificant..
            xx = Solver % Mesh % Nodes % x(i)
            yy = Solver % Mesh % Nodes % y(i)
            xt=xx+Displace(1)*dt
            yt=yy+Displace(2)*dt
            ! TO DO lot of inefficient coding here. better to assign RailX=xL if on left RailX=xR if on right and continue
            ! Q: need to allocate/deallocate all the time? (considering case where Nl is not Nr? and set RailN each time?
            ! IF LEFT set xrail = xL else ... etc
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
            DO j=1,Nrail
               !  TO DO save jmin,k from prev time and start checking distance from there
               ! for efficiency
               TempDist=(xRail(j)-xx)**2.+(yRail(j)-yy)**2.
               IF(TempDist < MinDist) THEN
                  MinDist=TempDist
                  jmin=j
               END IF
            END DO
            IF(jmin+1>Nrail) PRINT *, 'jmin=',jmin,'ERROR: advancing beyond left rail!'
            IF (MinDist < AEPS) THEN
               ! case 1)  mesh node is rail point'
               ! this likely happens at t=0.
               MovePastRailNode=.FALSE. !if this would be true here, rail resolution is smaller than mesh resolution
               ! we assume that's never the case. TO DO include assert  CALL Assert(1==0, SolverName, "1 does not equal zero!")
               ! check whether x+dx*dt is closest to segment jmin+1 or jmin-1  
               IF(PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                    (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                    > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                    (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
                  k=jmin-1
                  ELSE
                     k=jmin+1
                  END IF
               ELSE ! x is not a rail node
                  IF(PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin+1),yRail(jmin+1)/),(/xx,yy/)) &
                       > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin-1),yRail(jmin-1)/),(/xx,yy/))) THEN
                     IF (PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                       > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
                        MovePastRailNode=.FALSE. ! case 2) x+dx and x are on same segment'
                        k=jmin-1 ! x and x+dx are closest to jmin>jmin-1 segment
                     ELSE
                        ! case 3) advancing past rail segment
                        MovePastRailNode=.TRUE. ! TO DO should check whether projected
                        !displacement moves past rail node, not whether unprojected displacement does..?
                        k=jmin+1 ! x is moving from jmin-1>jmin to jmin>jmin+1
                     END IF
                  ELSE
                     IF (PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin+1),yRail(jmin+1)/),(/xt,yt/)) &
                       > PointLineSegmDist2D((/xRail(jmin),yRail(jmin)/), &
                       (/xRail(jmin-1),yRail(jmin-1)/),(/xt,yt/))) THEN
                        ! Case 3) x advancing past rail segment
                        MovePastRailNode=.TRUE.
                        k=jmin-1 ! x is moving from jmin+1>jmin to jmin>jmin-1
                     ELSE
                        ! Case 2) x+dx and x are on same segment jmin>jmin+1     
                        MovePastRailNode=.FALSE.
                        k=jmin+1
                     END IF
                  END IF
               END IF
            IF(MovePastRailNode) THEN
               ! We assume that the node moves from jmin>jmin-1 to jmin>jmin+1 or opposite
               ! but we exclude moving from jmin>jmin+1 to jmin+1>jmin+2 (or with minus)
               ! because then displacement in timestep is larger than rail resolution
               ! TO DO check this with CALL ASSERT
               ! we change the displacement such that the node first moves along the current rail segment
               ! to the nearest rail node. whatever magnitude of the horizontal part of the displacement is left
               ! will be the distance the node travels on the new segment.
               TempDist=( (xx-xRail(jmin))**2.+(yy-yRail(jmin))**2.)**0.5 ! distance travelled to reach new rail node
               TempDist= (Displace(1)**2.+Displace(2)**2.)**0.5 *dt - TempDist ! distance left to travel along new segment
               t=TempDist/( (xRail(k)-xRail(jmin))**2.+(yRail(k)-yRail(jmin))**2.)**0.5 ! relative distance on segment
               Displace(1)=t*(xRail(k)-xRail(jmin))
               Displace(2)=t*(yRail(k)-yRail(jmin))
               ! TO DO check if *dt should be here?! *dt is done later!!! now it's done twice!
               ! maybe move the whole IF MovePastRailNode part down after other Displacement*dt
               ! or put Displacement*dt up for all cases (front, and corner cases 1-3)
               ! TO DO above procedure can be improved to first calculate norm of displacement
               ! projected onto segm1 and then calculate how much displacement left on segm2 +
               ! project that as well!
            ELSE
               RailDir=(/xRail(jmin),yRail(jmin)/)-(/xRail(k),yRail(k)/)
               PRINT *, ' RailDir', RailDir
            Displace(1:2) = DOT_PRODUCT(Displace(1:2),RailDir) * &
                 RailDir / DOT_PRODUCT(RailDir,RailDir)
            END IF
            ! TO DO check whether x,y + displace*dt intersects with next line segment
            ! if it does, find point on next line segment with distance displ*dt from current x,y and send node there
            ! TO DO write warning if distance point to line of rail too large?
            ! TO DO also write warning if glacier advances out of range of rails?
            !  TO DO ASSERT that ||NewDisplace|| <= ||V||*dt    
            DEALLOCATE(xRail,yRail)
         END IF
         Displace = Displace * dt
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
