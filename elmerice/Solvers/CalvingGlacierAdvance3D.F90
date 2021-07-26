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
! *  Authors: Joe Todd
! *  Email:
! *  Web:     http://elmerice.elmerfem.org
! *
! *
! *****************************************************************************

!A routine for getting frontal advance in 3D, given velocity and melt

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
   INTEGER :: i, j, k, m, n, DOFs, TotalNodes,&
        FaceNodeCount, DummyInt,  hits, ierr, FrontLineCount, county,&
        ShiftIdx, ShiftToIdx, NoTangledGroups, PivotIdx, CornerIdx, &
        SecondIdx, FirstTangleIdx, LastTangleIdx
   INTEGER, POINTER :: Perm(:), FrontPerm(:)=>NULL(), TopPerm(:)=>NULL(), &
        FrontNodeNums(:)=>NULL(),LeftPerm(:)=>NULL(), RightPerm(:)=>NULL()
   INTEGER, ALLOCATABLE :: FrontLocalNodeNumbers(:), &
        NodeNumbers(:), TangledGroup(:), TangledPivotIdx(:), UpdatedDirection(:)
   REAL(KIND=dp) :: FrontOrientation(3),RotationMatrix(3,3),&
        NodeVelo(3), NodeMelt(3), NodeNormal(3), RailDir(2),&
        MeltRate, Displace(3), y_coord(2), epsShift, LongRangeLimit, MaxDisplacement, &
        EpsTangle,thisEps,Shift, thisY
   REAL(KIND=dp), POINTER :: Advance(:)
   REAL(KIND=dp), ALLOCATABLE :: Rot_y_coords(:,:), Rot_z_coords(:,:), &
        TangledShiftTo(:)
   LOGICAL :: Found, Debug, Parallel, Boss, ShiftLeft, LeftToRight, MovedOne, ShiftSecond, &
        Protrusion, SqueezeLeft, SqueezeRight, FirstTime=.TRUE., intersect_flag, FrontMelting
   LOGICAL, ALLOCATABLE :: DangerZone(:), WorkLogical(:), &
        Tangled(:), DontMove(:)
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VeloVarName, MeltVarName, &
        NormalVarName, FrontMaskName, TopMaskName, TangledVarName, &
        LeftMaskName, RightMaskName

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
            ! Strategy to project displace along rail direction:
            ! get x,y coordinate of node i, 
            ! compute direction of rail by finding prev and next rail node
            ! only do x,y direction; no penetration in vertical is handled by freesurf.
            ! NOTE mass is not conserved but for small timesteps and a relatively thin lateral boundary it should be insignificant..
            IF(LeftPerm(i)>0) THEN ! HARD CODED FOR PLANMESH, st RailDir not unitvector
               RailDir=(/-1.0,-5.0/)
            ELSE ! on right corner
               RailDir=(/1.0,-5.0/)
            END IF
            Displace(1:2) = DOT_PRODUCT(Displace(1:2),RailDir) * &
                 RailDir / DOT_PRODUCT(RailDir,RailDir)
            ! TO DO check whether x,y + displace*dt intersects with next line segment
            ! if it does, find point on next line segment with distance displ*dt from current x,y and send node there
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

!   DEALLOCATE(FrontPerm, TopPerm, LeftPerm, RightPerm)

 END SUBROUTINE GlacierAdvance3D
