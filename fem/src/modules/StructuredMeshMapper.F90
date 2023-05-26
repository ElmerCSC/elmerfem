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
! *  Authors: Peter Råback
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.3.2008
! *
! *****************************************************************************/

SUBROUTINE StructuredMeshMapper_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => GetSolverParams()  
  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)
  
END SUBROUTINE StructuredMeshMapper_init


!------------------------------------------------------------------------------
!>  Subroutine for mapping the mesh between given top and bottom surfaces.
!>  This solver assumes that the mesh is structural so that it could have 
!>  been obtained by extrusion in the direction of interest. For the given 
!>  direction the corresponding top and bottom node is computed for every node
!>  and this information is used to perform linear mapping in between.  
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE StructuredMeshMapper( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------

  USE CoordinateSystems
  USE MeshUtils
  USE ParallelUtils
  USE DefUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PSolver
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, TangledMaskVarName, MappedMeshName
  INTEGER :: i,j,k,n,dim,DOFs,itop,ibot,imid,ii,jj,Rounds,BotMode,TopMode,nsize, nnodes, &
       ActiveDirection,elem, istat, TangledCount, LimitedCount
  INTEGER, POINTER :: MaskPerm(:) => NULL(),TopPerm(:),BotPerm(:),TangledMaskPerm(:),TopPointer(:),&
       BotPointer(:),MidPointer(:),NodeIndexes(:),TmpPerm(:)
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
       DisplacementMode, MaskExists, GotVeloVar, GotUpdateVar, Tangled,&
       DeTangle, ComputeTangledMask = .FALSE., Reinitialize, &
       MidLayerExists, WriteMappedMeshToDisk = .FALSE., GotBaseVar, &
       BaseDisplaceFirst, RecompStab, RecompStabExe = .FALSE., &
       MapHeight, BotProj
  REAL(KIND=dp) :: UnitVector(3),x0loc,x0bot,x0top,x0mid,xloc,wtop,BotVal,TopVal,&
       TopVal0, BotVal0, MidVal, RefVal, ElemVector(3),DotPro,Eps,Length, MinHeight
  REAL(KIND=dp) :: at0,at1,at2,dx
  REAL(KIND=dp), POINTER :: Coord(:),BotField(:),TopField(:),TangledMask(:),CoordP(:)
  REAL(KIND=dp), ALLOCATABLE :: OrigCoord(:), Field(:), Surface(:)
  TYPE(Variable_t), POINTER :: Var, VeloVar, UpdateVar, TangledMaskVar, BaseVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t),POINTER :: BC

  INTEGER, POINTER :: FixedLayers(:),UpPointer(:),DownPointer(:),NodeLayer(:)
  INTEGER :: NumberOfLayers, NumberOfFixedLayers, RecompStabInterval, Cnt=0
  LOGICAL :: MultiLayer

  CHARACTER(*), PARAMETER :: Caller = 'StructuredMeshMapper'

  
  SAVE Visited,Initialized,UnitVector,Coord,MaskExists,MaskPerm,TopPointer,BotPointer,&
      TopMode,BotMode,TopField,BotField,TopPerm,BotPerm,Field,Surface,nsize,nnodes,OrigCoord, &
      ComputeTangledMask, MidPointer, MidLayerExists,&
      UpPointer,DownPointer,NodeLayer,NumberOfLayers, &
      RecompStabExe, Cnt

  CALL Info( Caller,'---------------------------------------',Level=4 )
  CALL Info( Caller,'Performing mapping on a structured mesh ',Level=4 )
  CALL Info( Caller,'---------------------------------------',Level=4 )

  !------------------------------------------------------------------------------
  !   Initialize the pointers to top and bottom nodes 
  !------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  PSolver => Solver
  Mesh => Solver % Mesh

  dim = Mesh % MeshDim
  
  Reinitialize = ListGetLogical(SolverParams, "Always Detect Structure", Found)
  IF( Reinitialize ) THEN
    IF( ALLOCATED(Field)) DEALLOCATE(Field)
    IF( ALLOCATED(Surface)) DEALLOCATE(Surface)
  END IF
  
  RecompStab = ListGetLogical(SolverParams, "Recompute Stabilization", Found)
  IF(.NOT. Found) THEN
    CALL Info(Caller,'Defaulting "Recompute Stabilization" to True.',Level=8)
    RecompStab = .TRUE.
  END IF
  IF (RecompStab) THEN
    RecompStabInterval = ListGetInteger(SolverParams, "Recompute Stabilization Interval", Found)
    IF (.NOT.Found) RecompStabInterval = 1
    Cnt = Cnt + 1
    IF (Cnt == RecompStabInterval) THEN
      Cnt = 0
      RecompStabExe = .TRUE.
    ELSE
      RecompStabExe = .FALSE.
    END IF
  ELSE
    RecompStabExe = .FALSE.
  END IF
    
  FixedLayers => ListGetIntegerArray( SolverParams,'Fixed Layer Indexes',MultiLayer)
  NumberOfFixedLayers = SIZE( FixedLayers )

  BotProj = ListGetLogical(SolverParams,'Project To Bottom',Found ) 
  
  IF( (.NOT. Initialized) .OR. Reinitialize ) THEN
    IF(ASSOCIATED(BotPointer)) DEALLOCATE(BotPointer)
    IF(ASSOCIATED(TopPointer)) DEALLOCATE(TopPointer)
    IF(ASSOCIATED(UpPointer)) DEALLOCATE(UpPointer)
    IF(ASSOCIATED(DownPointer)) DEALLOCATE(DownPointer)
    IF(ASSOCIATED(NodeLayer)) DEALLOCATE(NodeLayer)    
    IF( MultiLayer ) THEN
      CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
          TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
          UpNodePointer = UpPointer, DownNodePointer = DownPointer, &
          NumberOfLayers = NumberOfLayers, NodeLayer = NodeLayer )
      NumberOfLayers = NumberOfLayers + 1
      
      i = FixedLayers(1) 
      IF( i /= 1 ) THEN
        CALL Warn(Caller,'Enforcing first fixed layer to: 1 (was '//I2S(i)//')')
        FixedLayers(1) = 1
      END IF
      i = FixedLayers(NumberOfFixedLayers)
      IF( i /= NumberOfLayers ) THEN
        CALL Warn(Caller,'Enforcing last fixed layer to: '&
            //I2S(NumberOfLayers)//' (was '//I2S(i)//')')
        FixedLayers(NumberOfFixedLayers) = NumberOfLayers
      END IF
    ELSE
      CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
          TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
          MidNodePointer = MidPointer, MidLayerExists = MidLayerExists )
    END IF

    Coord => Var % Values

    ! For p-elements the number of nodes and coordinate vector differ
    ! The projection is implemented only for the true nodes
    nnodes = Mesh % NumberOfNodes
    nsize = MIN( SIZE( Coord ), Mesh % NumberOfNodes )
    Initialized = .TRUE.
   
    MaskExists = ASSOCIATED( Var % Perm )     
    IF( MaskExists ) THEN
      MaskPerm => Var % Perm
    ELSE
      IF(.NOT. ASSOCIATED(MaskPerm)) THEN
        ALLOCATE(MaskPerm(nsize))
        DO i=1,nsize
          MaskPerm(i) = i
        END DO
      END IF
    END IF

    IF(ALLOCATED(OrigCoord)) DEALLOCATE(OrigCoord)
    ALLOCATE( OrigCoord(nsize), STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL Fatal( Caller, 'Memory allocation error' )
    END IF
  END IF
  
  OrigCoord(1:nsize) = Coord(1:nsize)
  at0 = CPUTime()

  ! End of initialization
  !-------------------------------------------------------
  
  GotBaseVar = .FALSE.
  VarName = GetString( SolverParams,'Base Displacement Variable',Found)
  IF( Found ) THEN
    BaseVar => VariableGet( Mesh % Variables, VarName )
    GotBaseVar = ASSOCIATED( BaseVar )
    IF(.NOT. GotBaseVar ) THEN
      CALL Fatal(Caller,'The variable does not exist: '//TRIM(VarName))
    END IF

    BaseDisplaceFirst = GetLogical( SolverParams,'Base Displacement First',Found ) 
    IF( BaseDisplaceFirst ) THEN
      CALL Info(Caller,'Applying base displacement before structural mapping')
    END IF
  END IF
  

  ! Get the velocity variable component. 
  !-------------------------------------------------------------------
  VarName = GetString( SolverParams,'Mesh Velocity Variable',GotVeloVar)
  IF( GotVeloVar ) THEN
    VeloVar => VariableGet( Mesh % Variables, VarName )    
    IF( ASSOCIATED( VeloVar ) ) THEN
      IF( VeloVar % Dofs /= 1 ) THEN
        CALL Fatal(Caller,'The size of mesh velocity must be one')
      END IF
      IF( SIZE( VeloVar % Values ) /= nsize ) THEN
        CALL Fatal(Caller,'Sizes must agree: '//TRIM(VarName))
      END IF
    ELSE
      n = len_TRIM(VarName)
      ! Create full vector if the field is given by component
      IF( VERIFY( VarName(n:n),'123') == 0 ) THEN
        CALL DefaultVariableAdd( VarName(1:n-1), dofs = dim, Perm = MaskPerm )
        VeloVar => VariableGet( Mesh % Variables, VarName ) 
      ELSE
        CALL DefaultVariableAdd( VarName, Perm = MaskPerm, Var = VeloVar ) 
      END IF
    END IF
  END IF

  ! Get the mesh update variable component. 
  !-------------------------------------------------------------------
  VarName = GetString( SolverParams,'Mesh Update Variable',GotUpdateVar)
  IF( GotUpdateVar ) THEN
    UpdateVar => VariableGet( Mesh % Variables, VarName ) 
    IF( ASSOCIATED( UpdateVar ) ) THEN
      IF( UpdateVar % Dofs /= 1 ) THEN
        CALL Fatal(Caller,'The size of mesh update must be one')
      END IF
      IF( SIZE( UpdateVar % Values ) /= nsize ) THEN
        CALL Fatal(Caller,'Sizes must agree: '//TRIM(VarName))
      END IF
    ELSE
      n = len_TRIM(VarName)
      ! Create full vector if the field is given by component
      IF( VERIFY( VarName(n:n),'123') == 0 ) THEN
        CALL DefaultVariableAdd( VarName(1:n-1), dofs = dim, Perm = MaskPerm )
        UpdateVar => VariableGet( Mesh % Variables, VarName ) 
      ELSE
        CALL DefaultVariableAdd( VarName, Perm = MaskPerm, Var = UpdateVar )
      END IF
    END IF
  END IF
  
  DisplacementMode = GetLogical(SolverParams,'Displacement Mode',Found)

  MinHeight = GetCReal(SolverParams,'Minimum Mesh Height',GotIt)
  IF(.NOT. GotIt) MinHeight = GetCReal(SolverParams,'Minimum Height', GotIt)
  IF(.NOT. GotIt) MinHeight = EPSILON( MinHeight ) 

  WRITE(Message,'(A,E11.4)') 'Adjusting upper surface to maintain minimum height to:', MinHeight
  CALL Info(Caller,Message,Level=6)

  TangledCount = 0
  LimitedCount = 0
  
  IF( GotBaseVar ) THEN
    IF( BaseDisplaceFirst ) CALL BaseVarDisplace() 
  END IF

  IF( MultiLayer ) THEN
    CALL MultiLayerMapper()
  ELSE
    CALL BinaryLayerMapper()
  END IF

  LimitedCount = ParallelReduction(LimitedCount) 
  TangledCount = ParallelReduction(TangledCount) 
  
  IF( LimitedCount > 0 ) THEN
    CALL Info(Caller,'There seems to be '&
        //I2S(LimitedCount)//' (out of '//I2S(nsize)//&
        ') limited heights!',Level=6)
  END IF
  IF( TangledCount > 0 ) THEN
    CALL Info(Caller,'There seems to be '&
        //I2S(TangledCount)//' (out of '//I2S(nsize)//&
        ') tangled nodes!',Level=5)
  END IF    
 
  IF(ListGetLogical( SolverParams,'Mesh Mapping Passive',Found ) ) THEN
    CALL Info(Caller,'Taking back the suggested mapping!',Level=5)
    Coord(1:nsize) = OrigCoord(1:nsize)
  END IF
 
  ! If there is a mask then the coordinate is not directly linked to the real coordinate.
  ! Hence we need to do it here for the real coordinate. 
  IF( MaskExists ) THEN
    ActiveDirection = ListGetInteger( Solver % Values,'Active Coordinate')
    IF( ActiveDirection == 1 ) THEN
      CoordP => Mesh % Nodes % x
    ELSE IF( ActiveDirection == 2 ) THEN
      CoordP => Mesh % Nodes % y
    ELSE IF( ActiveDirection == 3 ) THEN
      CoordP => Mesh % Nodes % z
    ELSE
      CALL Fatal(Caller,'Unknown active coordinate!')
    END IF

    DO i=1,nnodes
      j = MaskPerm(i)
      IF( j == 0 ) CYCLE
      CoordP(i) = Coord(j)
    END DO
  END IF

  IF( GotBaseVar ) THEN
    IF (.NOT. BaseDisplaceFirst ) CALL BaseVarDisplace() 
  END IF
    
  IF( GotVeloVar .AND. .NOT. Visited ) THEN
    IF( GetLogical(SolverParams,'Mesh Velocity First Zero',Found ) ) THEN
      VeloVar % Values = 0.0_dp
    END IF
  END IF

  at1 = CPUTime()
  WRITE(Message,* ) 'Active coordinate mapping time: ',at1-at0
  CALL Info(Caller,Message)
  
  IF(.NOT. Visited ) THEN
    MappedMeshName = GetString(SolverParams,'Mapped Mesh Name', WriteMappedMeshToDisk)
    IF( WriteMappedMeshToDisk ) THEN
      CALL WriteMeshToDisk(Mesh, MappedMeshName)
    END IF
  END IF

  Visited = .TRUE.
  
  IF(RecompStabExe) CALL MeshStabParams(Mesh)

CONTAINS


  ! Perform mapping when top, bottom and possible middle layers are given.
  !-----------------------------------------------------------------------
  SUBROUTINE BinaryLayerMapper()

    !---------------- detangling stuff --------------------------------
    DeTangle = GetLogical(SolverParams,'Correct Surface',GotIt )
    IF( DeTangle ) THEN
      CALL Info(Caller,&
          '> Correct Surface < in case of intersecting upper and lower surface',Level=4)

      MapHeight = ListCheckPresent( SolverParams,'Mesh Height Map')
      IF(MapHeight) THEN
        CALL Info(Caller,'Using function to map heights',Level=5)
      END IF
      
      TangledMaskVarName = GetString(SolverParams,'Correct Surface Mask', ComputeTangledMask)
      IF (ComputeTangledMask) THEN
        TangledMaskVar => VariableGet( Mesh % Variables,  TRIM(TangledMaskVarName) )
        IF(.NOT. ASSOCIATED( TangledMaskVar ) ) THEN
          CALL Info(Caller,&
              'Given > Correct Surface Mask < variable not present, creating it.')
          ALLOCATE( TmpPerm(Mesh % NumberOfNodes) )
          DO i=1,Mesh % NumberOfNodes
            TmpPerm(i) = i
          END DO
          CALL DefaultVariableAdd( TangledMaskVarname, Perm = TmpPerm, Var = TangledMaskVar )
          NULLIFY( TmpPerm ) 
        END IF

        IF(TangledMaskVar % DOFs /= 1) THEN 
          CALL Fatal(Caller,'> Correct Surface Mask < variable should have only 1 dof')
        END IF
        TangledMask => TangledMaskVar % Values
        TangledMask = 1.0_dp
        TangledMaskPerm => TangledMaskVar % Perm
        WRITE(Message,'(A,A)') 
        CALL Info(Caller,&
            'Output of > Correct Surface Mask < to: '//TRIM(TangledMaskVarName),Level=6 )
      END IF
    END IF

    ! Get either variable or constant values for top surface
    !-------------------------------------------------------
    TopMode = 0
    TopVal0 = GetCReal(SolverParams,'Top Surface Level',GotIt)
    IF(GotIt) THEN
      TopMode = 1
    ELSE
      VarName = GetString(SolverParams,'Top Surface Variable Name',GotIt )
      IF(GotIt) THEN
        Var => VariableGet( Mesh % Variables,  VarName )
        IF(ASSOCIATED(Var)) THEN
          IF(Var % DOFs /= 1) THEN
            CALL Fatal(Caller,'Top surface variable should have only 1 dof')
          ELSE
            TopField => Var % Values
            TopPerm => Var % Perm
            TopMode = 2
          END IF
          IF( InfoActive( 20 ) ) THEN
            PRINT *,'TopField range:',MINVAL( TopField ), MAXVAL( TopField ), SUM( TopField ) / SIZE( TopField )
          END IF
        ELSE
          CALL Fatal(Caller,'Top surface variable is missing: '//TRIM(VarName))
        END IF
      END IF
    END IF

    IF(TopMode == 0) THEN
      IF( ListCheckPresentAnyBC( Model,'Top Surface') ) THEN
        TopMode = 3
        IF(.NOT. ALLOCATED(Field)) THEN
          N = Mesh % MaxElementNodes
          ALLOCATE(Field(nsize),Surface(n))
          Field = 0.0_dp
          Surface = 0.0_dp
        END IF
        DO elem = 1, Mesh % NumberOfBoundaryElements
          Element => GetBoundaryElement(elem)
          BC => GetBC()
          IF ( .NOT. ASSOCIATED( BC ) ) CYCLE
          NodeIndexes => Element % NodeIndexes
          n = GetElementNOFNodes()
          Surface(1:n) = GetReal( BC,'Top Surface',Found )
          IF(.NOT. Found) CYCLE
          
          IF( MaskExists ) THEN
            IF( ALL( MaskPerm(NodeIndexes(1:n)) > 0 ) ) THEN
              Field(MaskPerm(NodeIndexes(1:n))) = Surface(1:n) 
            END IF
          ELSE
            Field(NodeIndexes(1:n)) = Surface(1:n)              
          END IF
        END DO
      ELSE
        CALL Fatal(Caller,'Top surface BC entry is missing')
      END IF
    END IF

    ! Get either variable or constant values for bottom surface
    !----------------------------------------------------------
    BotMode = 0
    BotVal0 = GetCReal(SolverParams,'Bottom Surface Level',GotIt)
    IF(GotIt) THEN
      BotMode = 1
    ELSE
      VarName = GetString(SolverParams,'Bottom Surface Variable Name',GotIt )
      IF(GotIt) THEN
        Var => VariableGet( Mesh % Variables,  VarName )
        IF(ASSOCIATED(Var)) THEN
          IF( Var % DOFs /= 1) THEN
            CALL Fatal(Caller,'Bottom surface variable should have only 1 dof')
          ELSE
            BotField => Var % Values
            BotPerm => Var % Perm
            BotMode = 2
          END IF
          IF( InfoActive( 20 ) ) THEN
            PRINT *,'BotField range:',MINVAL( BotField ), MAXVAL( BotField ), SUM( BotField ) / SIZE( BotField )
          END IF
        ELSE
          CALL Fatal(Caller,'Bottom surface variable is missing: '//TRIM(VarName))
        END IF
      END IF
    END IF


    IF( BotMode == 0) THEN
      IF( ListCheckPresentAnyBC( Model,'Bottom Surface') ) THEN  
        BotMode = 3
        IF(.NOT. ALLOCATED(Field)) THEN
          N = Mesh % MaxElementNodes
          ALLOCATE(Field(nsize),Surface(n))
          Field = 0.0_dp
          Surface = 0.0_dp
        END IF
        DO elem = 1, Mesh % NumberOfBoundaryElements
          Element => GetBoundaryElement(elem)
          BC => GetBC()
          IF ( .NOT. ASSOCIATED( BC ) ) CYCLE

          NodeIndexes => Element % NodeIndexes
          n = GetElementNOFNodes()
          Surface(1:n) = GetReal( BC,'Bottom Surface',Found )
          IF(.NOT. Found) CYCLE

          IF( MaskExists ) THEN
            IF( ALL( MaskPerm(NodeIndexes(1:n)) > 0 ) ) THEN
              Field(MaskPerm(NodeIndexes(1:n))) = Surface(1:n) 
            END IF
          ELSE
            Field(NodeIndexes(1:n)) = Surface(1:n)              
          END IF
        END DO
      END IF
    END IF


    ! Get either variable or constant values for mid surface
    !----------------------------------------------------------
    IF( MidLayerExists ) THEN
      IF(.NOT. ALLOCATED(Field)) THEN
        N = Mesh % MaxElementNodes
        ALLOCATE(Field(nsize),Surface(n))
        Field = 0.0_dp
        Surface = 0.0_dp
      END IF
      DO elem = 1, Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(elem)
        BC => GetBC()
        IF ( ASSOCIATED( BC ) ) THEN
          NodeIndexes => Element % NodeIndexes
          n = GetElementNOFNodes()
          Surface(1:n) = GetReal( BC,'Mid Surface',Found )
          IF(Found) Field(NodeIndexes(1:n)) = Surface(1:n) 
        END IF
      END DO
    END IF
    
    ! Get the new mapping using linear interpolation from bottom and top
    !-------------------------------------------------------------------
    DO i=1,nnodes
      
      j = i
      IF( MaskExists ) THEN
        j = MaskPerm(i) 
        IF( j == 0) CYCLE
      END IF
      itop = TopPointer(j)
      ibot = BotPointer(j)

      IF( MidLayerExists ) imid = MidPointer(j)

      ! Use the previous coordinates for determining the weights
      !----------------------------------------------------------
      IF( MaskExists ) THEN
        x0top = OrigCoord(MaskPerm(itop))
        x0bot = OrigCoord(MaskPerm(ibot))
        x0loc = OrigCoord(MaskPerm(i))
      ELSE
        x0top = OrigCoord(itop)
        x0bot = OrigCoord(ibot)
        x0loc = OrigCoord(i)      
      END IF

      IF( TopMode == 1 ) THEN
        TopVal = TopVal0
      ELSE IF(TopMode == 2) THEN
        IF( TopPerm( itop ) == 0 ) THEN
          CALL Fatal(Caller,'Top surface variable perm is zero!')
        END IF
        TopVal = TopField(TopPerm(itop))
      ELSE IF(TopMode == 3) THEN
        IF( MaskExists ) THEN
          TopVal = Field(MaskPerm(itop))
        ELSE
          TopVal = Field(itop)
        END IF
      ELSE
        IF( DisplacementMode ) THEN
          TopVal = 0.0_dp
        ELSE
          TopVal = x0top
        END IF
      END IF

      IF( BotMode == 1 ) THEN
        BotVal = BotVal0 
      ELSE IF(BotMode == 2) THEN
        IF( BotPerm( ibot ) == 0 ) THEN
          CALL Fatal(Caller,'Bottom surface variable perm is zero!')
        END IF
        BotVal = BotField(BotPerm(ibot))
      ELSE IF(BotMode == 3) THEN    
        IF( MaskExists ) THEN
          BotVal = Field(MaskPerm(ibot))
        ELSE
          BotVal = Field(ibot)
        END IF
      ELSE
        IF( DisplacementMode ) THEN
          BotVal = 0.0_dp
        ELSE
          BotVal = x0bot
        END IF
      END IF

      IF( MidLayerExists ) THEN
        ! If we have a midlayer fetch the data for that too
        MidVal = Field(imid)
        x0mid = OrigCoord(imid)
      ELSE
        ! If we don't have midlayer set to to BotVal to make the tangled stuff easier
        MidVal = BotVal
        x0mid = x0bot
      END IF

      ! Check whether the mesh gets tangled
      ! Note that for midlayer existing we only check the upper part currently!
      IF( DisplacementMode ) THEN
        dx = TopVal + x0top - MidVal - x0mid
      ELSE
        dx = TopVal - MidVal 
      END IF
      Tangled = ( dx < MinHeight ) 
      
      IF( MaskExists .AND. Tangled ) THEN
        IF( DeTangle ) CALL Warn(Caller,'Cancelling tanglement when mask exists!')
        Tangled = .FALSE.
      END IF

      ! If the mesh is tangled then take some action.
      ! Here the lower surface stays intact. This is due to the main application field, 
      ! computational glaciology, where the lower surface of ice is usually nicely constrained. 
      IF( Tangled ) THEN        
        IF( dx < TINY( MinHeight ) ) THEN
          TangledCount = TangledCount + 1
        END IF
        LimitedCount = LimitedCount + 1

        IF( DeTangle ) THEN
          IF( MapHeight ) THEN
            ! Map the height from [-\infty,MinHeight] to [MinHeight/2,MinHeight], for example
            dx = ListGetFun( SolverParams,'Mesh Height Map', dx / MinHeight ) * MinHeight
          ELSE
            dx = MinHeight
          END IF

          ! Note again that if midlayer does not exist, the mid refers to bottom
          IF( DisplacementMode ) THEN
            TopVal = MidVal + x0mid + x0top + dx
          ELSE
            TopVal = MidVal + dx
          END IF
        END IF

        IF (ComputeTangledMask) THEN
          TangledMask(TangledMaskPerm(i)) = -1.0_dp 
        END IF
        IF( .FALSE. ) THEN
          WRITE(Message,'(A,E11.4,A,E11.4,A,E11.4,A,E11.4)')&
              "Corrected negative height:", TopVal - MidVal, "=",&
              TopVal ,"-", MidVal, ". New upper value:", Field(itop)
          CALL Info(Caller,Message,Level=9)
        END IF
      END IF

      ! New coordinate location
      IF( MidLayerExists ) THEN
        ! With middle layer in two parts, first the upper part
        IF( (x0top - x0mid ) * ( x0loc - x0mid ) > 0.0_dp ) THEN
          wtop = (x0loc-x0mid)/(x0top-x0mid);
          xloc = wtop * TopVal + (1.0_dp - wtop) * MidVal         
        ELSE
          wtop = (x0loc-x0bot)/(x0mid-x0bot);
          xloc = wtop * MidVal + (1.0_dp - wtop) * BotVal
        END IF
      ELSE
        ! Otherwise in one part
        wtop = (x0loc-x0bot)/(x0top-x0bot);
        xloc = wtop * TopVal + (1.0_dp - wtop) * BotVal
      END IF
      
      IF(DisplacementMode) THEN
        IF( GotVeloVar ) THEN
          IF(Velovar % Perm(i)>0) &
              VeloVar % Values( VeloVar % Perm(i) ) = xloc / dt
        END IF
        Coord(j) = OrigCoord(j) + xloc
      ELSE
        IF( GotVeloVar ) THEN
          IF(Velovar % Perm(i)>0) &
              VeloVar % Values( VeloVar % Perm(i) ) = ( xloc - OrigCoord(i) ) / dt
        END IF
        Coord(j) = xloc
      END IF
      IF( GotUpdateVar ) UpdateVar % Values ( UpdateVar % Perm(i) ) = Coord(j) - OrigCoord(j)
    END DO

    
  END SUBROUTINE BinaryLayerMapper

  
  SUBROUTINE MultiLayerMapper()
    REAL(KIND=dp), ALLOCATABLE :: Proj(:,:), StrideCoord(:),FixedCoord(:)
    INTEGER, ALLOCATABLE :: StrideInd(:),StridePerm(:)
    INTEGER :: ierr, PEs
    LOGICAL :: Hit
    REAL(KIND=dp) :: q
    TYPE(Variable_t), POINTER :: FixedVar    
    INTEGER :: status(MPI_STATUS_SIZE)
        
    ! Get the new mapping using linear interpolation from bottom and top
    !-------------------------------------------------------------------
    CALL Info(Caller,'Mapping using '//I2S(NumberOfFixedLayers)//' fixed layers',Level=6)

    IF( MaskExists ) THEN
      CALL Fatal(Caller,'Mask not available yet for multiple layers!')
    END IF
    
    DeTangle = GetLogical(SolverParams,'Correct Surface',GotIt )
    
    VarName = ListGetString( SolverParams,'Fixed Layer Variable',UnfoundFatal = .TRUE. )
    FixedVar => VariableGet( Mesh % Variables, VarName ) 
    IF(.NOT. ASSOCIATED( FixedVar ) ) THEN
      CALL Fatal(Caller,'Could not find variable: '//TRIM(VarName) )
    END IF
    IF( FixedVar % Dofs /= NumberOfFixedLayers ) THEN
      CALL Fatal(Caller,'Invalid number of components in fixed layer variable:'&
          //I2S(FixedVar % Dofs))
    END IF

    ALLOCATE( Proj(NumberOfLayers,NumberOfFixedLayers),StrideInd(NumberOfLayers),&
        StridePerm(NumberOfLayers),StrideCoord(NumberOfLayers),&
        FixedCoord(NumberOfFixedLayers))
    Proj = 0.0_dp
    
    ! Go through all 1D strides and perform mapping for mesh
    DO i=1,nnodes
      ibot = BotPointer(i)

      ! Start mapping from bottom
      IF( ibot /= i ) CYCLE

      ! Create stride for this column
      j = ibot
      StrideCoord(1) = OrigCoord(j)
      StrideInd(1) = j
      DO k = 2,NumberOfLayers
        j = UpPointer(j)
        StrideCoord(k) = OrigCoord(j)
        StrideInd(k) = j
        itop = j
      END DO

      ! Create a new projection matrix for this column
      j = 1    
      DO k = 1, NumberOfLayers
        Hit = .FALSE.
        DO j = 1, NumberOfFixedLayers+1
          IF( FixedLayers(j) == k ) THEN
            Proj(k,j) = 1.0_dp
            Hit = .TRUE.
          ELSE IF( FixedLayers(j) < k .AND. FixedLayers(j+1) > k ) THEN
            q = 1.0_dp*(StrideCoord(k)-StrideCoord(FixedLayers(j))) / &
                (StrideCoord(FixedLayers(j+1))-StrideCoord(FixedLayers(j)))
            Proj(k,j+1) = q
            Proj(k,j) = 1-q
            Hit = .TRUE.
          END IF
          IF( Hit ) EXIT
        END DO
        IF(.NOT. Hit ) THEN
          CALL Fatal(Caller,'Could not find mapping for layer: '//I2S(k))
        END IF
      END DO
      
      ! We can either have the fixed layer variable at top or bottom, not elsewhere!
      k = FixedVar % Perm(ibot)
      IF( k == 0 ) k = FixedVar % Perm(itop)
      IF( k == 0 ) CALL Fatal(Caller,'Could not find fixed fields at top or bottom!')

      ! Get the given layers
      FixedCoord = FixedVar % Values(NumberOfFixedLayers*(k-1)+1:NumberOfFixedLayers*k)
      
      IF( DisplacementMode ) THEN
        StrideCoord = StrideCoord + MATMUL( Proj, FixedCoord ) 
      ELSE
        StrideCoord = MATMUL( Proj, FixedCoord )
      END IF

      ! Fix mesh if it becomes tangled
      IF( DeTangle ) THEN
        IF( BotProj ) THEN
          k = COUNT( StrideCoord(1:NumberOfLayers-1)-StrideCoord(2:NumberOfLayers) < MinHeight )
          IF( k > 0 ) THEN
            LimitedCount = LimitedCount + k
            TangledCount = TangledCount + &
                COUNT( StrideCoord(1:NumberOfLayers-1)-StrideCoord(2:NumberOfLayers) < TINY(MinHeight) )  
            DO k = 2,NumberOfLayers
              StrideCoord(k) = MIN( StrideCoord(k), StrideCoord(k-1)-MinHeight )
            END DO
          END IF
        ELSE
          k = COUNT( StrideCoord(2:NumberOfLayers)-StrideCoord(1:NumberOfLayers-1) < MinHeight )
          IF( k > 0 ) THEN
            LimitedCount = LimitedCount + k    
            TangledCount = TangledCount + &
                COUNT( StrideCoord(2:NumberOfLayers)-StrideCoord(1:NumberOfLayers-1) < TINY(MinHeight) )
            DO k = 2,NumberOfLayers
              StrideCoord(k) = MAX( StrideCoord(k), StrideCoord(k-1)+MinHeight )
            END DO
          END IF
        END IF
      END IF
              
      Coord(StrideInd) = StrideCoord

      IF( GotVeloVar ) THEN
        StridePerm = VeloVar % Perm( StrideInd ) 
        WHERE( StridePerm > 0 ) 
          VeloVar % Values( StridePerm ) = ( Coord(StrideInd) - OrigCoord(StrideInd) ) / dt
        END WHERE
      END IF
        
      IF( GotUpdateVar ) THEN
        StridePerm = UpdateVar % Perm( StrideInd ) 
        WHERE( StridePerm > 0 ) 
          UpdateVar % Values( StridePerm ) = Coord(StrideInd) - OrigCoord(StrideInd)
        END WHERE
      END IF
    END DO

    CALL Info(Caller,'Finished multilayer mapping',Level=8)
    

  END SUBROUTINE MultiLayerMapper

  
  
  SUBROUTINE BaseVarDisplace()

    INTEGER :: dofs
    
    CALL Info(Caller,'Adding base displacement to the displacements!')

    dofs = BaseVar % Dofs
    IF( dofs /=2 .AND. dofs /= 3 ) THEN
      CALL Fatal(Caller,'Invalid base displacement dimension: '//I2S(dofs))
    END IF

    DO i=1,nsize
      j = i
      IF( MaskExists ) THEN
        j = MaskPerm(i) 
        IF( j == 0) CYCLE
      END IF
      ibot = BotPointer(i)
      
      Mesh % Nodes % x(i) = Mesh % Nodes % x(i) + &
          BaseVar % Values( dofs*(BaseVar % Perm(ibot)-1)+1 )  
      Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + &
          BaseVar % Values( dofs*(BaseVar % Perm(ibot)-1)+2 )  
      IF( dofs == 3 ) THEN
        Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + &
            BaseVar % Values( dofs*BaseVar % Perm(ibot) )  
      END IF
    END DO
       
  END SUBROUTINE BaseVarDisplace
  
  
  !------------------------------------------------------------------------------
END SUBROUTINE StructuredMeshMapper
!------------------------------------------------------------------------------


