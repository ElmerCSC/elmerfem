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
! *  Authors: Peter R�back
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.3.2008
! *
! *****************************************************************************/

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
       ActiveDirection,elem, istat, TangledCount
  INTEGER, POINTER :: MaskPerm(:),TopPerm(:),BotPerm(:),TangledMaskPerm(:),TopPointer(:),&
       BotPointer(:),MidPointer(:),NodeIndexes(:)
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
       DisplacementMode, MaskExists, GotVeloVar, GotUpdateVar, Tangled,&
       DeTangle, ComputeTangledMask = .FALSE., Reinitialize, &
       MidLayerExists, WriteMappedMeshToDisk = .FALSE., GotBaseVar, &
       BaseDisplaceFirst, RecompStab
  REAL(KIND=dp) :: UnitVector(3),x0loc,x0bot,x0top,x0mid,xloc,wtop,BotVal,TopVal,&
       TopVal0, BotVal0, MidVal, ElemVector(3),DotPro,Eps,Length, MinHeight
  REAL(KIND=dp) :: at0,at1,at2,Heps
#ifndef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: CPUTime,RealTime
#endif
  REAL(KIND=dp), POINTER :: Coord(:),BotField(:),TopField(:),TangledMask(:),CoordP(:)
  REAL(KIND=dp), ALLOCATABLE :: OrigCoord(:), Field(:), Surface(:)
  TYPE(Variable_t), POINTER :: Var, VeloVar, UpdateVar, TangledMaskVar, BaseVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t),POINTER :: BC

  SAVE Visited,Initialized,UnitVector,Coord,MaskExists,MaskPerm,TopPointer,BotPointer,&
      TopMode,BotMode,TopField,BotField,TopPerm,BotPerm,Field,Surface,nsize,nnodes,OrigCoord, &
      ComputeTangledMask, MidPointer, MidLayerExists

  CALL Info( 'StructuredMeshMapper','---------------------------------------',Level=4 )
  CALL Info( 'StructuredMeshMapper','Performing mapping on a structured mesh ',Level=4 )
  CALL Info( 'StructuredMeshMapper','---------------------------------------',Level=4 )

  !------------------------------------------------------------------------------
  !   Initialize the pointers to top and bottom nodes 
  !------------------------------------------------------------------------------

  SolverParams => GetSolverParams()
  PSolver => Solver
  Mesh => Solver % Mesh

  dim = Mesh % MeshDim
  
  Reinitialize = ListGetLogical(SolverParams, "Always Detect Structure", Found)
  IF(.NOT. Found) Reinitialize = .FALSE.

  RecompStab = ListGetLogical(SolverParams, "Recompute Stabilization", Found)
  IF(.NOT. Found) RecompStab = .FALSE.

  IF( (.NOT. Initialized) .OR. Reinitialize ) THEN
    IF(ASSOCIATED(BotPointer)) DEALLOCATE(BotPointer)
    IF(ASSOCIATED(TopPointer)) DEALLOCATE(TopPointer)
    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
        TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
        MidNodePointer = MidPointer, MidLayerExists = MidLayerExists )
    MaskExists = ASSOCIATED( Var % Perm ) 
    IF( MaskExists ) MaskPerm => Var % Perm
    Coord => Var % Values
    
    ! For p-elements the number of nodes and coordinate vector differ
    ! The projection is implemented only for the true nodes
    nnodes = Mesh % NumberOfNodes
    nsize = MIN( SIZE( Coord ), Mesh % NumberOfNodes )
    Initialized = .TRUE.

    IF(ALLOCATED(OrigCoord)) DEALLOCATE(OrigCoord)
    ALLOCATE( OrigCoord(nsize), STAT=istat)
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'StructuredMeshMapper', 'Memory allocation error' )
    END IF
  END IF
  
  OrigCoord(1:nsize) = Coord(1:nsize)
  at0 = CPUTime()

  ! End of initialization
  !-------------------------------------------------------


  MappedMeshName = GetString(SolverParams,'Mapped Mesh Name', WriteMappedMeshToDisk)

  !---------------- detangling stuff --------------------------------
  MinHeight = 0.0_dp
  DeTangle = GetLogical(SolverParams,'Correct Surface',GotIt )
  TangledCount = 0
  IF( DeTangle ) THEN
    CALL Info('StructuredMeshMapper',&
        '> Correct Surface < in case of intersecting upper and lower surface',Level=4)
    MinHeight = GetCReal(SolverParams,'Minimum Height', GotIt)
    IF (.NOT.GotIt .OR. (MinHeight <= 0.0_dp)) THEN
      CALL Fatal('StructuredMeshMapper',&
          '> Minimum Height < either set to negative/zero value or not found')
    ELSE
      WRITE(Message,'(A,E11.4)') 'Adjusting upper surface to maintain minimum height to:', MinHeight
      CALL Info('StructuredMeshMapper',Message,Level=4)
    END IF

    TangledMaskVarName = GetString(SolverParams,'Correct Surface Mask', ComputeTangledMask)
    IF (ComputeTangledMask) THEN
      TangledMaskVar => VariableGet( Mesh % Variables,  TRIM(TangledMaskVarName) )
      IF (ASSOCIATED(TangledMaskVar)) THEN
        IF(TangledMaskVar % DOFs /= 1) THEN 
          CALL Fatal('StructuredMeshMapper','> Correct Surface Mask < variable should have only 1 dof')
        END IF
        TangledMask => TangledMaskVar % Values
        TangledMask = 1.0_dp
        TangledMaskPerm => TangledMaskVar % Perm
        WRITE(Message,'(A,A)') 
        CALL Info('StructuredMeshMapper',&
            'Output of > Correct Surface Mask < to: '//TRIM(TangledMaskVarName),Level=6 )
      ELSE
        CALL Warn('StructuredMeshMapper',&
            'Ignoring '//TRIM(TangledMaskVarName)//' given as > Correct Surface Mask < variable, as not found.')
      END IF
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
          CALL Fatal('StructuredMeshMapper','Top surface variable should have only 1 dof')
        ELSE
          TopField => Var % Values
          TopPerm => Var % Perm
          TopMode = 2
        END IF
      ELSE
        CALL Fatal('StructuredMeshMapper','Top surface variable is missing: '//TRIM(VarName))
      END IF
    END IF
  END IF

  IF(TopMode == 0) THEN
    IF( ListCheckPresentAnyBC( Model,'Top Surface') ) THEN
      TopMode = 3
      IF( Reinitialize ) THEN
        IF( ALLOCATED(Field)) DEALLOCATE(Field)
        IF( ALLOCATED(Surface)) DEALLOCATE(Surface)
      END IF
      IF(.NOT. ALLOCATED(Field)) THEN
        N = Mesh % MaxElementNodes
        ALLOCATE(Field(nsize),Surface(n))
        Field = 0.0_dp
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
      CALL Fatal('StructuredMeshMapper','Top surface BC entry is missing')
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
          CALL Fatal('StructuredMeshMapper','Bottom surface variable should have only 1 dof')
        ELSE
          BotField => Var % Values
          BotPerm => Var % Perm
          BotMode = 2
        END IF
      ELSE
        CALL Fatal('StructuredMeshMapper','Bottom surface variable is missing: '//TRIM(VarName))
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

  GotBaseVar = .FALSE.
  VarName = GetString( SolverParams,'Base Displacement Variable',Found)
  IF( Found ) THEN
    BaseVar => VariableGet( Mesh % Variables, VarName )
    GotBaseVar = ASSOCIATED( BaseVar )
    IF(.NOT. GotBaseVar ) THEN
      CALL Fatal('StructuredMeshMapper','The variable does not exist: '//TRIM(VarName))
    END IF

    BaseDisplaceFirst = GetLogical( SolverParams,'Base Displacement First',Found ) 
    IF( BaseDisplaceFirst ) THEN
      CALL Info('StructuredMeshMapper','Applying base displacement before structural mapping')
    END IF
  END IF
  

  ! Get the velocity variable component. 
  !-------------------------------------------------------------------
  GotVeloVar = .FALSE.
  VarName = GetString( SolverParams,'Mesh Velocity Variable',Found)
  IF( Found ) THEN
    VeloVar => VariableGet( Mesh % Variables, VarName ) 
    IF( ASSOCIATED( VeloVar ) ) THEN
      IF( VeloVar % Dofs == 1 ) THEN
        GotVeloVar = .TRUE.
      ELSE  
        CALL Fatal('StructuredMeshMapper','The size of mesh velocity must be one')
      END IF
    ELSE
      CALL Fatal('StructuredMeshMapper','The variable does not exist: '//TRIM(VarName))
    END IF
  END IF

  ! Get the mesh update variable component. 
  !-------------------------------------------------------------------
  GotUpdateVar = .FALSE.

  VarName = GetString( SolverParams,'Mesh Update Variable',Found)
  IF( Found ) THEN
    UpdateVar => VariableGet( Mesh % Variables, VarName ) 
    IF( ASSOCIATED( UpdateVar ) ) THEN
      IF( UpdateVar % Dofs == 1 ) THEN
        GotUpdateVar = .TRUE.
      ELSE  
        CALL Fatal('StructuredMeshMapper','The size of mesh update must be one')
      END IF
    ELSE
      CALL Fatal('StructuredMeshMapper','The variable does not exist: '//TRIM(VarName))
    END IF
  END IF

  DisplacementMode = GetLogical(SolverParams,'Displacement Mode',Found)

  IF( GotBaseVar .AND. BaseDisplaceFirst ) THEN
    CALL BaseVarDisplace() 
  END IF


  ! Get the new mapping using linear interpolation from bottom and top
  !-------------------------------------------------------------------
  Heps = GetCReal(SolverParams,'Minimum Mesh Height',GotIt)
  IF(.NOT. GotIt) Heps = EPSILON( Heps )

  TangledCount = 0

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
        CALL Fatal('StructuredMeshMapper','Top surface variable perm is zero!')
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
        CALL Fatal('StructuredMeshMapper','Bottom surface variable perm is zero!')
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

    ! If we have a midlayer fetch the data for that too
    IF( MidLayerExists ) THEN
      MidVal = Field(imid)
      x0mid = OrigCoord(imid)
    END IF

    ! Check whether the mesh gets tangled
    ! Note that for midlayer existing we only check the upper part currently!
    IF( MidLayerExists ) THEN
      IF( DisplacementMode ) THEN
        Tangled = ( TopVal + x0top < MidVal + x0mid + MinHeight)
      ELSE
        Tangled = ( TopVal < MidVal + MinHeight) 
      END IF      
    ELSE      
      IF( DisplacementMode ) THEN
        Tangled = ( TopVal + x0top < BotVal + x0bot + MinHeight)
      ELSE
        Tangled = ( TopVal < BotVal + MinHeight) 
      END IF
    END IF

    IF( MaskExists .AND. Tangled ) THEN
      IF( DeTangle ) CALL Warn('StructuredMeshMapper','Cancelling tanglement when mask exists!')
      Tangled = .FALSE.
    END IF

    
    ! If the mesh is tangled then take some action.
    ! Here the lower surface stays intact. This is due to the main application field, 
    ! computational glaciology, where the lower surface of ice is usually nicely constrained. 
    IF( Tangled ) THEN
      TangledCount = TangledCount + 1

      IF( DeTangle ) THEN
        IF( MidLayerExists ) THEN
          IF( DisplacementMode ) THEN
            TopVal = MidVal + x0mid + x0top + MinHeight
          ELSE
            TopVal = MidVal + MinHeight
          END IF
        ELSE
          IF( DisplacementMode ) THEN
            TopVal = BotVal + x0bot + x0top + MinHeight
          ELSE
            TopVal = BotVal + MinHeight
          END IF
        END IF
      END IF

      IF (ComputeTangledMask) THEN
        TangledMask(TangledMaskPerm(i)) = -1.0_dp 
      END IF
      IF( .FALSE. ) THEN
        WRITE(Message,'(A,E11.4,A,E11.4,A,E11.4,A,E11.4)')&
            "Corrected negative height:", TopVal - BotVal, "=",&
            TopVal ,"-", BotVal, ". New upper value:", Field(itop)
        CALL Info('StructuredMeshMapper',Message,Level=9)
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
      CALL Fatal('StructuredMeshMapper','Unknown active coordinate!')
    END IF

    DO i=1,nnodes
      j = MaskPerm(i)
      IF( j == 0 ) CYCLE
      CoordP(i) = Coord(j)
    END DO
  END IF


  
  IF( GotBaseVar .AND. .NOT. BaseDisplaceFirst ) THEN
    CALL BaseVarDisplace() 
  END IF
  
  
  IF( GotVeloVar .AND. .NOT. Visited ) THEN
    IF( GetLogical(SolverParams,'Mesh Velocity First Zero',Found ) ) THEN
      VeloVar % Values = 0.0_dp
    END IF
  END IF

  IF( TangledCount > 0 ) THEN
    CALL Warn('StructuredMeshMapper','There seems to be '&
        //TRIM(I2S(TangledCount))//' (out of '//TRIM(I2S(nsize))//&
        ') tangled nodes!')
  END IF

  at1 = CPUTime()
  WRITE(Message,* ) 'Active coordinate mapping time: ',at1-at0
  CALL Info('StructuredMeshMapper',Message)

  IF ( (.NOT.Visited) .AND. WriteMappedMeshToDisk ) THEN
     CALL WriteMeshToDisk(Mesh, MappedMeshName)
  END IF

  Visited = .TRUE.

  IF(RecompStab) CALL MeshStabParams(Mesh)

CONTAINS

  SUBROUTINE BaseVarDisplace()

    CALL Info('StructuredMeshMapper','Adding base displacement to the displacements!')
    DO i=1,nsize
      j = i
      IF( MaskExists ) THEN
        j = MaskPerm(i) 
        IF( j == 0) CYCLE
      END IF
      ibot = BotPointer(i)

      IF( BaseVar % Dofs == 2 ) THEN
        Mesh % Nodes % x(i) = Mesh % Nodes % x(i) + &
            BaseVar % Values( 2*BaseVar % Perm(ibot)-1 )  
        Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + &
            BaseVar % Values( 2*BaseVar % Perm(ibot) )  
      ELSE
        CALL Fatal('StructuredMeshMapper','Base displacement assumed to be two!')
      END IF
    END DO
       
  END SUBROUTINE BaseVarDisplace


  
  !------------------------------------------------------------------------------
END SUBROUTINE StructuredMeshMapper
!------------------------------------------------------------------------------
