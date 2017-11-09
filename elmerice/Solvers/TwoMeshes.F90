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
! *  Mesh adaptation & variable interpolation for 2D calving model. Original
! *  notes from Peter below: 
! *
! *  This solver may be used to map the field variables of a stretched mesh
! *  to the unstretched initial configuarion of the same mesh. In addition
! *  the vertical size of the mesh is accommodated so that the height of 
! *  the mesh stays fixed. All-in-all two different mapping stages and 
! *  one solution of Laplace equation is needed for the mapping.
! *  The user may also optionally nullify some field variables such as the 
! *  Mesh Update which should freshly start from the new configuration. 
! *
! *  Notes: 
! *  ----------
! *  This solver is currently only implemented for the 2D flowline case.
! *  For 3D calving, see the Calving3D test case.
! *  In places, the code would appear to permit 3D, but these haven't been 
! *  properly implemented or tested.
! *
! *  Currently InterpolateMeshToMesh will not report points for which 
! *  interpolation failed (although InterpolateVarToVarReduced *will*).
! *  When the new variable is created, it's values are set to -HUGE() 
! *  in the hope that it will be obvious to the user that a problem has
! *  occurred, though of course the long term solution should be to report
! *  failure... This is probably only a minor issue, as no points should 
! *  ever be missing, as this would indicate a larger problem with 
! *  this routine.
! *  
! *  In addition to dealing with the mesh manipulation required for a calving 
! *  event, this subroutine (in conjunction with Calving.F90), can also modify
! *  the mesh when a frontal melting overhang becomes so severe as to force 
! *  'calving front' nodes below 'bed' nodes.  This is done via the 'RemeshOccurs'
! *  switch.
! ******************************************************************************
! *
! *  Authors: Peter R�back, Joe Todd
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 19.10.2010
! *
! ****************************************************************************/

SUBROUTINE TwoMeshes( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CRSMatrix
  USE GeneralUtils
  USE ElementDescription
  USE MeshUtils  
  USE InterpVarToVar

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: MeshA => NULL(), MeshB => NULL(), Mesh0 => Null(),&
       MeshTop => Null(), MeshBot => Null()
  TYPE(Variable_t), POINTER :: Var, Var2, VarNew, VarTop, VarBot, VarDisplace1, VarDisplace2, &
       TopVar, BotVar, FrontVar, FirstBCVar, LastBodyVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), POINTER :: NodesA, NodesB, Nodes0, NodesTop, NodesBot, PreNodes, PostNodes
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, VarName, TopMaskName, BotMaskName, FrontMaskName, SolverName
  REAL(KIND=dp) :: BoundingBox(6),eps1,eps2, Norm, &
      TopDisplacement, BotDisplacement, H, B, T, y, p, &
      PreArea, PostArea, TotalLoss, BergLength, BergHeight
  INTEGER :: i,j,n,dim,istat,HeightDim, NoNodes, TopNodes, BotNodes, &
      OldBotNodes, FrontNodes, active, TopCornerIndex, BotCornerIndex, County
  INTEGER, POINTER :: TopPerm(:), BotPerm(:), FrontPerm(:), HeightPerm(:), &
       FieldPerm(:), WorkInt(:)
  REAL(KIND=dp), POINTER :: TopValues(:), BotValues(:), FrontValues(:), Height(:), ForceVector(:), Field(:),AuxReal(:)
  TYPE(Matrix_t), POINTER  :: StiffMatrix
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
  LOGICAL :: Found,RebuiltQuadrantTree,FSTop,FSBot,NeedInit=.FALSE., &
       MapCondition, RemeshCondition, Parallel, Debug = .FALSE., FirstTime = .TRUE.
  LOGICAL, POINTER :: UnfoundNodes(:)

  SAVE MeshA, MeshB, Mesh0, MeshTop, MeshBot, NodesA, NodesB, Nodes0, &
      NodesTop, NodesBot, dim, HeightDim, TopCornerIndex, BotCornerIndex, &
      TopPerm, BotPerm, FrontPerm, TopValues, BotValues, FrontValues, NoNodes,&
      TopNodes, BotNodes, TopMaskName, BotMaskName, Height, HeightPerm, &
      STIFF, FORCE, FirstTime, FSTop, FSBot, TotalLoss, PreNodes, PostNodes, &
      FirstBCVar, LastBodyVar

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

!------------------------------------------------------------------------------
  
  Debug = .FALSE.
  SolverName = "TwoMeshes"

  CALL Info( SolverName, '-----------------------------------' )
  CALL Info( SolverName, ' Adapting the calving mesh...' )
  CALL Info( SolverName, '-----------------------------------' )

  dim = CoordinateSystemDimension()
  IF(dim /= 2) THEN
     CALL Fatal(SolverName, "This solver only works in 2D, sorry!")
  END IF
  !This is somewhat redundant but it permits the possibility of future generalisation
  HeightDim = dim


  Parallel = ParEnv % Pes > 1
  IF(Parallel) CALL Fatal(SolverName, "Solver doesnt work in parallel yet!")

  IF ( FirstTime ) THEN
     FirstTime = .FALSE.

     !variable to store total mass lost through calving
     TotalLoss = 0.0_dp

     !If both top and bottom surfaces are free, we simply interpolate 
     !them from the previous (t-1) timestep.  However, repeated interpolation
     !on a fixed surface will lead to artificial smoothing, so we must store
     !the t=0 surface and interpolate from that each time.
     FSTop = ListGetLogical( Solver % Values, 'FS Top', Found)
     IF(.NOT. Found) THEN
        Call Warn(SolverName,'No value "FS Top" found, assuming Free Surface at surface.')
        FSTop = .TRUE.
     END IF

     FSBot = ListGetLogical( Solver % Values, 'FS Bottom', Found)
     IF(.NOT. Found) THEN
        Call Warn(SolverName,'No value "FS Bottom" found, assuming no Free Surface at bed.')
        FSBot = .FALSE.
     END IF

     IF((.NOT. FSBot).OR.(.NOT. FSTop))THEN
        NeedInit = .TRUE.
     END IF

     MeshA => Solver % Mesh
     NodesA => MeshA % Nodes

     MeshB => AllocateMesh()
     MeshB = MeshA
     MeshB % Name = TRIM(MeshA % Name)//'_copy'

     CALL Info(SolverName,'Allocated a new copy of the mesh')

     IF(NeedInit)THEN
        Mesh0 => AllocateMesh()
        Mesh0 = MeshA
        Mesh0 % Name = TRIM(Solver % Mesh % Name)//'_initial'
        CALL Info(SolverName,'Created initial reference mesh because at least one surface is fixed')
     END IF

     NoNodes = SIZE( MeshA % Nodes % x )
     NULLIFY( MeshB % Nodes )
     ALLOCATE( NodesB )
     ALLOCATE( NodesB % x(NoNodes), NodesB % y(NoNodes), NodesB % z(NoNodes) )
     NodesB % x = NodesA % x
     NodesB % y = NodesA % y
     NodesB % z = NodesA % z
     MeshB % Nodes => NodesB

     ALLOCATE( Nodes0 ) 
     ALLOCATE( Nodes0 % x(NoNodes), Nodes0 % y(NoNodes), Nodes0 % z(NoNodes) )
     Nodes0 % x = NodesA % x
     Nodes0 % y = NodesA % y
     Nodes0 % z = NodesA % z
     IF(NeedInit) Mesh0 % Nodes => Nodes0

     !Get the height variable for geometry interpolation
     HeightPerm => Solver % Variable % Perm
     Height => Solver % Variable % Values

     !The convenience variables MeshTop/Bot and NodesTop/Bot remove the
     !need for multiple IF statements later on.  
     !------------------------------------------------------------
     IF(FSTop) THEN
        MeshTop => MeshA
     ELSE
        MeshTop => Mesh0
     END IF

     IF(FSBot) THEN
        MeshBot => MeshA
     ELSE
        MeshBot => Mesh0
     END IF

     NodesTop => MeshTop % Nodes
     NodesBot => MeshBot % Nodes

     ! Nullify the list of variables to take use of the automatic
     ! features of the mapping
     NULLIFY( MeshB % Variables )

     ! Create the variables in MeshB
     ! These need to be there for successful mapping
     !--------------------------------------------------------------
     DO i=1,99

        WRITE( Name, '(A,I0)') 'Variable ',i
        VarName = ListGetString( Solver % Values,TRIM(Name),Found)
        IF(.NOT. Found ) EXIT

        Var => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE.) 
        IF( ASSOCIATED( Var) ) THEN
           NULLIFY( Field, FieldPerm )
           ALLOCATE( Field(SIZE(Var % Values)), FieldPerm(SIZE(Var % Perm)) )  

           !TODO InterpolateMeshToMesh should be made to report missing points in 
           !interpolation, as InterpVarToVar currently does.
           Field = -HUGE(Field)
           FieldPerm = Var % Perm

           CALL VariableAdd( MeshB % Variables, MeshB, Solver, TRIM( VarName), 1, Field, FieldPerm )

           IF(ASSOCIATED(Var % PrevValues)) THEN
              VarNew => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. )
              ALLOCATE( VarNew % PrevValues ( SIZE( Var % PrevValues,1), &
                                              SIZE( Var % PrevValues,2) ) )
           END IF
           CALL Info(SolverName,'Created variable: '//TRIM( VarName ) )
        ELSE
           WRITE(Message,'(a,a,a)') "Requested variable: ", TRIM(VarName), " but it wasn't found!"
           CALL Warn(SolverName, Message)
        END IF
     END DO

     DO i=1,99
        !Variables created here (listed in the sif Solver section) are those 
        !which are only defined on the upper and lower surface
        !They are interpolated using InterpolateVarToVarReduced

        WRITE( Name, '(A,I0)') 'Surface Variable ',i
        VarName = ListGetString( Solver % Values,TRIM(Name),Found)
        IF(.NOT. Found ) EXIT

        Var => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 

        IF( ASSOCIATED( Var) ) THEN
           NULLIFY( Field, FieldPerm)
           ALLOCATE( Field(SIZE(Var % Values)), FieldPerm(SIZE(Var % Perm)) ) 
           Field = Var % Values
           !Field = 0.0_dp
           FieldPerm = Var % Perm

           CALL VariableAdd( MeshB % Variables, MeshB, Solver, TRIM( VarName), 1, Field, FieldPerm )

           IF(i==1) THEN
              !Keep a handle on this so we can split variable list when interpolating
              FirstBCVar => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. )
           END IF

           IF(ASSOCIATED(Var % PrevValues)) THEN
              VarNew => VariableGet( MeshB % Variables, TRIM( VarName), ThisOnly = .TRUE. )
              ALLOCATE( VarNew % PrevValues ( SIZE( Var % PrevValues,1),&
                                              SIZE( Var % PrevValues,2) ) )
           END IF
           CALL Info(SolverName,'Created surface variable: '//TRIM( VarName ) )
        ELSE
           WRITE(Message,'(a,a,a)') "Requested surface variable: ", &
                TRIM(VarName), " but it wasn't found!"
           CALL Warn(SolverName, Message)        
        END IF
     END DO

     !Handle to last body variable, same as above
     Var => MeshB % Variables
     DO WHILE(ASSOCIATED(Var))
        IF(ASSOCIATED(Var % Next, FirstBCVar)) THEN
           LastBodyVar => Var
           EXIT
        END IF
        Var => Var % Next
     END DO

     ! Create variable for the top surface
     !---------------------------------------------------
     TopMaskName = 'Top Surface Mask'
     TopVar => VariableGet(MeshTop % Variables, TopMaskName, ThisOnly=.TRUE.)

     IF(ASSOCIATED(TopVar)) THEN
        TopValues => TopVar % Values
        TopPerm => TopVar % Perm
        TopNodes = COUNT(TopPerm > 0)
     ELSE
        ALLOCATE( TopPerm(NoNodes) )
        CALL MakePermUsingMask( Model,Solver,MeshTop,TopMaskName, &
             .FALSE., TopPerm, TopNodes )

        ALLOCATE( TopValues(TopNodes) )
        CALL VariableAdd( MeshTop % Variables, MeshTop, Solver, &
             TopMaskName, 1, TopValues, TopPerm )
     END IF
     
     ! Create variable for the bottom surface
     !---------------------------------------------------
     BotMaskName = 'Bottom Surface Mask'
     BotVar => VariableGet(MeshBot % Variables, BotMaskName, ThisOnly=.TRUE.)

     IF(ASSOCIATED(BotVar)) THEN
        BotValues => BotVar % Values
        BotPerm => BotVar % Perm
        BotNodes = COUNT(BotPerm > 0)
     ELSE
        ALLOCATE( BotPerm(NoNodes) )
        CALL MakePermUsingMask( Model,Solver,MeshBot,BotMaskName, &
             .FALSE., BotPerm, BotNodes )

        ALLOCATE( BotValues(BotNodes) )
        CALL VariableAdd( MeshBot % Variables, MeshBot, Solver, &
             BotMaskName, 1, BotValues, BotPerm )
     END IF

     ! Create variable for the calving front
     !---------------------------------------------------
     FrontMaskName = 'Calving Front Mask'
     FrontVar => VariableGet(MeshBot % Variables, FrontMaskName, ThisOnly=.TRUE.)

     IF(ASSOCIATED(FrontVar)) THEN
        FrontValues => FrontVar % Values
        FrontPerm => FrontVar % Perm
        FrontNodes = COUNT(FrontPerm > 0)
     ELSE
        ALLOCATE( FrontPerm(NoNodes) )
        CALL MakePermUsingMask( Model,Solver,MeshBot,FrontMaskName, &
             .FALSE., FrontPerm, FrontNodes )

        ALLOCATE( FrontValues(FrontNodes) )
        CALL VariableAdd( MeshA % Variables, MeshA, Solver, &
             FrontMaskName, 1, FrontValues, FrontPerm )
     END IF


     ! Find two corner nodes at calving front
     !---------------------------------------------------
     DO i=1,NoNodes
        IF(TopPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
           TopCornerIndex = i
        END IF
        IF(BotPerm(i) > 0 .AND. FrontPerm(i) > 0) THEN
           BotCornerIndex = i
        END IF
     END DO

     n = MeshA % MaxElementNodes 
     ALLOCATE( FORCE(n), STIFF(n,n), STAT=istat )

     !Allocate node holders for calculating calving event size later
     ALLOCATE(PreNodes,PostNodes)
     ALLOCATE(PreNodes % x(2),&
          PreNodes % y(2),&
          PostNodes % x(2),&
          PostNodes % y(2))

  END IF !FirstTime


  ! Add a suitable condition here for calving
  !---------------------------------------------
  ! This has been linked to the Calving.f90 solver by means of the "CalvingOccurs"
  ! logical in the list.
  MapCondition = ListGetLogical( Model % Simulation, 'CalvingOccurs', Found)
  IF(.NOT. Found) THEN
     CALL Warn("Two Meshes","Can't find CalvingOccurs Logical, exiting... ")
     RETURN
  END IF

  RemeshCondition = ListGetLogical( Model % Simulation, 'RemeshOccurs', Found)
  IF(.NOT. Found) THEN
     CALL Warn("Two Meshes","Can't find RemeshCondition Logical, assuming false!")
     RemeshCondition = .FALSE.
  END IF

  IF( .NOT. (MapCondition .OR. RemeshCondition)) RETURN

  PreArea = ModelArea(MeshA)

  !An array of 2 nodes, 1 being top, 2 bottom.
  !This is stored to check if berg is tabular or not...
  PreNodes % x(1) = MeshA % Nodes % x(TopCornerIndex)
  PreNodes % y(1) = MeshA % Nodes % y(TopCornerIndex)
  PreNodes % x(2) = MeshA % Nodes % x(BotCornerIndex)
  PreNodes % y(2) = MeshA % Nodes % y(BotCornerIndex)

  PRINT *,'Initial ranges'
  PRINT *,'X0:',MINVAL( Nodes0 % x), MAXVAL( Nodes0 % x)
  PRINT *,'XA:',MINVAL( NodesA % x), MAXVAL( NodesA % x)
  PRINT *,'XB:',MINVAL( NodesB % x), MAXVAL( NodesB % x)

  ! First copy the top and bottom (and front) height to variables
  !----------------------------------------------------------
  DO i=1,NoNodes
    j = TopPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        TopValues(j) = NodesTop % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        TopValues(j) = NodesTop % y(i)
      ELSE
        TopValues(j) = NodesTop % z(i)
      END IF
    END IF

    j = BotPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        BotValues(j) = NodesBot % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        BotValues(j) = NodesBot % y(i)
      ELSE
        BotValues(j) = NodesBot % z(i)
      END IF
    END IF

    j = FrontPerm(i)
    IF( j > 0 ) THEN
      IF( HeightDim == 1 ) THEN
        FrontValues(j) = NodesA % x(i)
      ELSE IF( HeightDim == 2 ) THEN
        FrontValues(j) = NodesA % y(i)
      ELSE
        FrontValues(j) = NodesA % z(i)
      END IF
    END IF
  END DO

  ! Map the top and bottom variables to the reference mesh
  ! For that purpose set the coordinates of MeshB to initial ones
  ! Stretching may be accounted for here.
  !--------------------------------------------------------------
  NodesB % x = NodesA % x
  NodesB % y = NodesA % y
  NodesB % z = NodesA % z

  !At this point, our nodes are displaced due to calving.  If we were dealing only
  !with vertical calving faces and uniform calving retreat (i.e. a perfect rectangle
  !falls off) then we can simply say NodesB % x = NodesB % x - CalvingEventSize
  !However, non-vertical calving faces require a more complex treatment.

  !FrontDisplacement returns the translation between current pre-calving node 
  !positions and the post-calving mesh whose node positions are the result of 
  !mesh displacement on the *initial* reference mesh.  This preserves mesh 
  !quality after repeated calving/advance cycles.
  VarDisplace1 => VariableGet(MeshA % Variables, 'Front Displacement 1', ThisOnly=.TRUE.)
  VarDisplace2 => VariableGet(MeshA % Variables, 'Front Displacement 2', ThisOnly=.TRUE.)

  NodesB % x = NodesA % x + VarDisplace1 % Values(VarDisplace1 % Perm)
  NodesB % y = NodesA % y + VarDisplace2 % Values(VarDisplace2 % Perm)

  !In the case where remeshing is required to deal with a severe melt overhang
  !the bottom corner node may actually 'advance' (though no physical advance of
  !the front actually takes place).
  !In this case, its height must be interpolated from some severely tilted front
  !nodes, and so we temporarily change their BotPerm here.
  IF(RemeshCondition) THEN
     IF(Debug) PRINT *, 'Debug Calving: RemeshCondition TRUE'
     !Save BotNodes so it can be restored...
     OldBotNodes = BotNodes

     !Adding all frontnodes, except the bottom corner which is already in
     BotNodes = BotNodes + FrontNodes - 1

     ALLOCATE(AuxReal(BotNodes))
     AuxReal = BotValues !(BotNodes - OldBotNodes) slots left empty

     County = 0
     DO i=1,NoNodes
        IF(FrontPerm(i) == 0) CYCLE  !if it's not a front node
        IF(BotPerm(i) /= 0) CYCLE  !if it already exists in bottom variable

        County = County + 1
        BotPerm(i) = OldBotNodes + County

        AuxReal(BotPerm(i)) = MeshA % Nodes % y(i)
        IF(Debug) THEN
           PRINT *, 'Debug ',SolverName,': Added NEW AuxReal nodenum: ',i,&
                ' botperm: ',BotPerm(i),' and AuxReal: ',AuxReal(BotPerm(i))
        END IF
     END DO

     BotVar => VariableGet(MeshBot % Variables, BotMaskName, ThisOnly=.TRUE.)
     BotVar % Values => AuxReal
     
  END IF

  !this is somewhat obtuse, but it allows InterpolateVarToVarReduced to 
  !reduce either 1 or 2 dimensions.
  ALLOCATE(WorkInt(1)); WorkInt = HeightDim;
  UnfoundNodes => NULL()
  CALL InterpolateVarToVarReduced(MeshTop, MeshB, TopMaskName, WorkInt, &
       UnfoundNodes, Variables=MeshA % Variables)

  IF(ANY(UnfoundNodes)) CALL Fatal(SolverName,&
       "Failed to find all nodes in top surface linesearch")
  CALL Info(SolverName,'Interpolated the top values to the original mesh')

  CALL InterpolateVarToVarReduced(MeshBot, MeshB, BotMaskName, WorkInt, &
       UnfoundNodes, Variables=MeshA % Variables)

  IF(ANY(UnfoundNodes)) CALL Fatal(SolverName,&
       "Failed to find all nodes in bottom surface linesearch")
  CALL Info(SolverName,'Interpolated the bottom values to the original mesh')

  IF(RemeshCondition) THEN
     VarBot => VariableGet(MeshB % Variables, BotMaskName, ThisOnly = .TRUE.)
     !Reset surplus BotPerm values
     DO i=1,NoNodes
        IF(BotPerm(i) > OldBotNodes) THEN
           BotPerm(i) = 0
           VarBot % Perm(i) = 0
        END IF
     END DO
     BotVar % Values => BotValues  !Reset variable values from AuxReal

     IF(Debug) THEN
        PRINT *, 'Debug ',SolverName,', BotVar % Values: ', BotVar % Values
        PRINT *, 'Debug ',SolverName,', BotVar % Perm: ', BotVar % Perm     
     END IF
  END IF

  ! Copy interpolated height variables back to TopValues and BotValues
  !----------------------------------------------------------

  VarTop => VariableGet(MeshB % Variables, TopMaskName, ThisOnly = .TRUE.)
  VarBot => VariableGet(MeshB % Variables, BotMaskName, ThisOnly = .TRUE.)
  DO i=1,NoNodes
    j = TopPerm(i)
    IF( j > 0 ) THEN
        TopValues(j) = VarTop % Values(VarTop % Perm(i)) 
    END IF
    j = BotPerm(i) 
    IF( j > 0 ) THEN
        BotValues(j) = VarBot % Values(VarBot % Perm(i)) 
    END IF
  END DO

  ! Find the y-shift in the calving front nodes
  ! using the two end values from MeshA and B to guide the deformation
  ! This is necessary for situations where the calving front is not 
  ! vertical, as the nodes will tend towards the outwardmost end otherwise
  !----------------------------------------------------------

  H = Nodes0 % y(TopCornerIndex) - Nodes0 % y(BotCornerIndex)
  T = Nodes0 % y(TopCornerIndex)
  B = Nodes0 % y(BotCornerIndex)

  TopDisplacement = VarTop % Values(VarTop % Perm(TopCornerIndex)) - NodesB % y(TopCornerIndex)
  BotDisplacement = VarBot % Values(VarBot % Perm(BotCornerIndex)) - NodesB % y(BotCornerIndex)

  DO i=1,NoNodes
     j = FrontPerm(i)
     IF( j > 0 ) THEN
        y = Nodes0 % y(i)
        !p is 1 at the top, 0 at the bottom
        p = ((y-B)/H)
        FrontValues(j) = NodesB % y(i) + (p*TopDisplacement) + ((1-p)*BotDisplacement) 
     END IF
  END DO


  ! Then map the mesh using Laplace equation in height direction
  !----------------------------------------------------------
  CALL Info(SolverName,'Solving the height from poisson equation')
  Solver % Mesh % Nodes => NodesB

  CALL DefaultInitialize()

  active = GetNOFActive()
  DO i=1,active
    Element => GetActiveElement(i)
    n = GetElementNOFNodes()
    CALL LocalMatrix(  STIFF, FORCE, Element, n )
    CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishAssembly()


  ! Set top and bottom Dirichlet BCs using the ones mapped from the
  ! deformed mesh. This eliminates the need for external BCs. 
  !-----------------------------------------------------------------  
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS

  DO i=1, NoNodes
    j =  TopPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, TopValues(j) )                
    END IF
    j = BotPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, BotValues(j) )                
    END IF
    j = FrontPerm(i)
    IF( j > 0 ) THEN
      CALL SetDirichtletPoint( StiffMatrix, ForceVector,1,1, &
          HeightPerm, i, FrontValues(j) )                
    END IF
  END DO

  Norm = DefaultSolve()

  ! Return the nodes to the original ones in which all other variables 
  ! have been computed in order not to case problems later on...
  Solver % Mesh % Nodes => NodesA

  ! Copy the height now to the original mesh deformed to comply
  ! with the deformed top and bottom surfaces
  !------------------------------------------------------------
  DO i=1,NoNodes
    j = HeightPerm(i)
    IF( j == 0 ) CYCLE
    IF( HeightDim == 1 ) THEN
      NodesB % x(i) = Height(j)
    ELSE IF( HeightDim == 2 ) THEN
      NodesB % y(i) = Height(j)
    ELSE
      NodesB % z(i) = Height(j)
    END IF
  END DO


  ! Map the mesh to the real positions
  ! Always built a new tree as MeshA has changed
  ! This would ideally be built inside the interpolation...
  !--------------------------------------------------------------
  RebuiltQuadrantTree = .TRUE.
  IF( RebuiltQuadrantTree ) THEN
    BoundingBox(1) = MINVAL( NodesA % x )
    BoundingBox(2) = MINVAL( NodesA % y )
    BoundingBox(3) = MINVAL( NodesA % z )
    BoundingBox(4) = MAXVAL( NodesA % x )
    BoundingBox(5) = MAXVAL( NodesA % y )
    BoundingBox(6) = MAXVAL( NodesA % z )

    eps1 = 0.1_dp
    eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
    BoundingBox(1:3) = BoundingBox(1:3) - eps2
    BoundingBox(4:6) = BoundingBox(4:6) + eps2

    IF(ASSOCIATED(MeshA % RootQuadrant)) CALL FreeQuadrantTree( MeshA % RootQuadrant )
    CALL BuildQuadrantTree( MeshA,BoundingBox,MeshA % RootQuadrant)
  END IF

  !Manipulate variable list to prevent BC variables being interpolated,
  !they're already interpolated by InterpolateVarToVarReduced
  LastBodyVar % Next => NULL()

  ! This one uses the standard interpolation routines with two meshes
  ! Note that the 4th argument (the new meshes variable) is never actually used.
  CALL InterpolateMeshToMesh( MeshA, MeshB, MeshA % Variables,MeshB % Variables, .TRUE. ) 
  CALL Info('TwoMesh','Interpolation done')

  !Restore variable list
  LastBodyVar % Next => FirstBCVar

  ! Now we still have the problem that the interpolated values sit on MeshB while
  ! we would like to work with MeshA. It seems easiest to copy all the relevant
  ! stuff back to MeshA. 
  !------------------------------------------------------------------------------

  NodesA % x = NodesB % x
  NodesA % y = NodesB % y
  NodesA % z = NodesB % z


  Var => MeshB % Variables
  DO WHILE(ASSOCIATED(Var))
     Var2 => VariableGet(MeshA % Variables, TRIM(Var % Name), ThisOnly = .TRUE.)
     IF(.NOT. ASSOCIATED(Var2)) CALL Fatal(SolverName, "Error while copying back variables")

     CALL Info(SolverName,'Copying back variable: '//TRIM( Var % Name ) )
     Var2 % Values = Var % Values
     PRINT *,'Range',MINVAL( Var2 % Values), MAXVAL( Var2 % Values )

     Var => Var % Next
  END DO

  DO i=1,99
     WRITE( Name, '(A,I0)') 'Nullify ',i
     VarName = ListGetString( Solver % Values,TRIM(Name),Found)
     IF(.NOT. Found ) EXIT
     Var2 => VariableGet( MeshA % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
     IF( ASSOCIATED( Var2) ) THEN
        CALL Info(SolverName,'Zeroing variable: '//TRIM( VarName ) )
        Var2 % Values = 0.0_dp
     ELSE
        WRITE(Message,'(a,a,a)') "Requested nullified variable: ", &
             TRIM(VarName), " but it wasn't found!"
        CALL Warn(SolverName, Message)
     END IF
  END DO

  MeshA % Changed = .TRUE.

  PRINT *,'Final ranges'
  PRINT *,'X0:',MINVAL( Nodes0 % x), MAXVAL( Nodes0 % x)
  PRINT *,'XA:',MINVAL( NodesA % x), MAXVAL( NodesA % x)
  PRINT *,'XB:',MINVAL( NodesB % x), MAXVAL( NodesB % x)

  !Calculate and print berg size
  PostArea = ModelArea(MeshA)
  TotalLoss = TotalLoss + (PreArea-Postarea)

  WRITE(Message, '(a,E10.4,a)') 'Calving event size: ',PreArea-PostArea,' m2'
  CALL Info(SolverName, Message)
  WRITE(Message, '(a,E10.4,a)') 'Total mass lost through calving: ',TotalLoss,' m2'
  CALL Info(SolverName, Message)

  !Check if tabular berg or not.
  PostNodes % x(1) = MeshA % Nodes % x(TopCornerIndex)
  PostNodes % y(1) = MeshA % Nodes % y(TopCornerIndex)
  PostNodes % x(2) = MeshA % Nodes % x(BotCornerIndex)
  PostNodes % y(2) = MeshA % Nodes % y(BotCornerIndex)

  BergLength = ((PreNodes % x(1) - PostNodes % x(1)) + (PreNodes % x(2) - PostNodes % x(2))) * 0.5
  BergHeight = ((PreNodes % y(1) - PreNodes % y(2)) + ( PostNodes % y(1) - PostNodes % y(2))) * 0.5

  IF(BergLength >= BergHeight) THEN
     WRITE(Message,'(a)') 'Calving event type: tabular'
  ELSE
     WRITE(Message,'(a)') 'Calving event type: normal'
  END IF
  CALL Info(SolverName, Message)

  DEALLOCATE(WorkInt, UnfoundNodes)
  IF(RemeshCondition) DEALLOCATE(AuxReal)

  CALL Info(SolverName,'All done')

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3),DetJ
    LOGICAL :: Stat
    INTEGER :: t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    
    FORCE = 0.0_dp
    STIFF = 0.0_dp

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
          IP % W(t),  detJ, Basis, dBasisdx )
      STIFF(1:n,1:n) = STIFF(1:n,1:n) + IP % s(t) * DetJ * &
          MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix


!------------------------------------------------------------------------------
  SUBROUTINE SetDirichtletPoint( StiffMatrix, ForceVector,DOF, NDOFs, &
      Perm, NodeIndex, NodeValue) 
!------------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(Matrix_t), POINTER :: StiffMatrix
    REAL(KIND=dp) :: ForceVector(:), NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
!------------------------------------------------------------------------------

    INTEGER :: PermIndex
    REAL(KIND=dp) :: s

!------------------------------------------------------------------------------

    PermIndex = Perm(NodeIndex)
    
    IF ( PermIndex > 0 ) THEN
      PermIndex = NDOFs * (PermIndex-1) + DOF
      
      IF ( StiffMatrix % FORMAT == MATRIX_SBAND ) THEN        
        CALL SBand_SetDirichlet( StiffMatrix,ForceVector,PermIndex,NodeValue )        
      ELSE IF ( StiffMatrix % FORMAT == MATRIX_CRS .AND. &
          StiffMatrix % Symmetric ) THEN        
        CALL CRS_SetSymmDirichlet(StiffMatrix,ForceVector,PermIndex,NodeValue)        
      ELSE                          
        s = StiffMatrix % Values(StiffMatrix % Diag(PermIndex))
        ForceVector(PermIndex) = NodeValue * s
        CALL ZeroRow( StiffMatrix,PermIndex )
        CALL SetMatrixElement( StiffMatrix,PermIndex,PermIndex,1.0d0*s )        
      END IF
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichtletPoint
!------------------------------------------------------------------------------

  FUNCTION ModelArea(Mesh) RESULT(A)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i
    REAL(KIND=dp) :: A

    A = 0.0_dp

    DO i=1,Mesh % NumberOfBulkElements

       A = A + ElementArea(Mesh, Mesh % Elements(i), &
            Mesh % Elements(i) % Type % NumberOfNodes)
    END DO

  END FUNCTION ModelArea
!------------------------------------------------------------------------------
END SUBROUTINE TwoMeshes
