!/******************************************************************************
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! *
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write
! *  to the Free Software Foundation, Inc., 51 Franklin Street,
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *****************************************************************************/

!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Coordinate system transformations, rigid mesh mapping, extruded structure detection.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshTransform

    USE MeshBasics
    IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE DetectExtrudedStructure( Mesh, Solver, ExtVar, &
      TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, &
      MidNodePointer, MidLayerExists, NumberOfLayers, NodeLayer, &
      MaskVar )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Solver_t) :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopNodePointer(:), BotNodePointer(:), &
        UpNodePointer(:), DownNodePointer(:), MidNodePointer(:)
    INTEGER, POINTER, OPTIONAL :: NodeLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
    LOGICAL, OPTIONAL :: MidLayerExists
    TYPE(Variable_t), TARGET, OPTIONAL :: MaskVar
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    TYPE(Nodes_t), POINTER :: MeshNodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, nnodes, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
        UpHit, DownHit, bc_ind, jmin, jmax, elemmax
    INTEGER, POINTER :: NodeIndexes(:), MaskPerm(:)
    LOGICAL :: MaskExists, UpActive, DownActive, GotIt, Found, DoCoordTransform
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1, Length, UnitVector(3), Vector(3), Vector2(3), &
        ElemVector(3), DotPro, MaxDotPro, MinDotPro, Eps, MinTop, &
        MaxTop, MinBot, MaxBot
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(:), ALLOCATABLE :: VarName, CoordTransform
    CHARACTER(*), PARAMETER :: Caller="DetectExtrudedStructure"
   
    CALL Info(Caller,'Determining extruded structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim
    Params => Solver % Values
    
    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal('StructuredMeshMapper','Invalid value for Active Coordinate')
    END IF  
    UnitVector = 0.0_dp
    UnitVector(ActiveDirection) = 1.0_dp

    IF( ListGetLogical(Params,'Mapping Original Coordinates',Found ) ) THEN
      MeshNodes => Mesh % NodesOrig
    ELSE
      MeshNodes => Mesh % Nodes
    END IF
    
    IF( ListGetLogical(Params,'Project To Bottom',GotIt) ) &
        UnitVector = -1.0_dp * UnitVector

    WRITE(Message,'(A,3F8.3)') 'Unit vector of direction:',UnitVector
    CALL Info(Caller,Message,Level=8)

    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-4

    nnodes = Mesh % NumberOfNodes
    nsize = nnodes

    Var => NULL()
    IF( PRESENT(MaskVar) ) THEN
      Var => MaskVar
    ELSE          
      VarName = ListGetString(Params,'Mapping Mask Variable',GotIt )
      IF(GotIt) THEN
        Var => VariableGet( Mesh % Variables,  VarName )
      END IF
    END IF
    MaskExists = ASSOCIATED(Var)
    IF( MaskExists ) THEN
      ALLOCATE( MaskPerm( SIZE( Var % Perm ) ) )
      MaskPerm = Var % Perm 
      nsize = MAXVAL( MaskPerm ) 
      CALL Info(Caller,'Using variable as mask: '//TRIM(Var % Name),Level=8)
    ELSE
      VarName = ListGetString(Params,'Mapping Mask Name',MaskExists )
      IF( MaskExists ) THEN
        CALL Info(Caller,'Using name as mask: '//TRIM(VarName),Level=8)
        MaskPerm => NULL() 
        CALL MakePermUsingMask( CurrentModel, Solver, Mesh, VarName, &
            .FALSE., MaskPerm, nsize )
        !PRINT *,'nsize:',nsize,SIZE(MaskPerm),MAXVAL(MaskPerm(1:nnodes))
      END IF
    END IF

    IF( MaskExists ) THEN
      CALL Info(Caller,'Applying mask of size: '//I2S(nsize),Level=10)
    ELSE
      CALL Info(Caller,'Applying extrusion on the whole mesh',Level=10)
    END IF 

    CoordTransform = ListGetString(Params,'Mapping Coordinate Transformation',DoCoordTransform )
    IF( DoCoordTransform .OR. MaskExists) THEN
      Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info(Caller,'Reusing > Extruded Coordinate < variable',Level=12 )
        Values => Var % Values        
      ELSE
        NULLIFY( Values )
        ALLOCATE( Values( nsize ) )
        Values = 0.0_dp
        IF( MaskExists ) THEN
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values, MaskPerm)
        ELSE
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values)
        END IF
        Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      END IF
    ELSE IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF

    CALL Info(Caller,'Variable used to detect extrusion: '//TRIM(Var % Name),Level=10)
    IF( MaskExists .OR. DoCoordTransform) THEN
      DO i=1,Mesh % NumberOfNodes
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
        END IF
        Vector(1) = Mesh % Nodes % x(i)
        Vector(2) = Mesh % Nodes % y(i)
        Vector(3) = Mesh % Nodes % z(i)
        IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF
        Values(j) = Vector( ActiveDirection )
      END DO
    END IF
    IF( PRESENT( ExtVar ) ) ExtVar => Var
    
    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpNodePointer) .OR. PRESENT ( TopNodePointer ) 
    DownActive = PRESENT( DownNodePointer) .OR. PRESENT ( BotNodePointer ) 
    
    IF( PRESENT( NumberOfLayers) .OR. PRESENT( NodeLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nnodes
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        TopPointer(j) = i
        UpPointer(j) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nnodes        
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        BotPointer(j) = i
        DownPointer(j) = i
      END DO
    END IF
    
    CALL Info(Caller,'Determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    IF( MaskExists ) THEN
      elemmax = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
    ELSE
      elemmax = Mesh % NumberOfBulkElements 
    END IF

    
    DO elem = 1,elemmax
      
      Element => Mesh % Elements(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element
      
      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = MeshNodes % x(NodeIndexes)
      Nodes % y(1:n) = MeshNodes % y(NodeIndexes)
      Nodes % z(1:n) = MeshNodes % z(NodeIndexes)
      
      ! This is probably a copy-paste error, I comment it away for time being.   
      ! IF (.NOT. (Element % PartIndex == Parenv % Mype) ) CYCLE

      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(NodeIndexes) == 0) ) CYCLE
      END IF
      
      DO i=1,n
        ii = NodeIndexes(i)
        
        Vector(1) = Nodes % x(i)
        Vector(2) = Nodes % y(i) 
        Vector(3) = Nodes % z(i)
        
        IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF

        MaxDotPro = -1.0_dp
        MinDotPro = 1.0_dp
        
        DO j=i+1,n
          jj = NodeIndexes(j)
          
          Vector2(1) = Nodes % x(j)
          Vector2(2) = Nodes % y(j)
          Vector2(3) = Nodes % z(j)

          IF( DoCoordTransform ) THEN
            CALL CoordinateTransformationNodal( CoordTransform, Vector2 )
          END IF
          
          ElemVector = Vector2 - Vector

          Length = SQRT(SUM(ElemVector*ElemVector))
          DotPro = SUM(ElemVector * UnitVector) / Length

          IF( DotPro > MaxDotPro ) THEN
            MaxDotPro = DotPro
            jmax = jj
          END IF
          IF( DotPro < MinDotPro ) THEN
            MinDotPro = DotPro
            jmin = jj
          END IF          
        END DO
          
        IF(MaxDotPro > 1.0_dp - Eps) THEN 
          IF( MaskExists ) THEN
            IF( UpActive ) UpPointer(MaskPerm(ii)) = jmax
            IF( DownActive ) DownPointer(MaskPerm(jmax)) = ii              
          ELSE
            IF( UpActive ) UpPointer(ii) = jmax
            IF( DownActive ) DownPointer(jmax) = ii
          END IF
        END IF
            
        IF(MinDotPro < Eps - 1.0_dp) THEN
          IF( MaskExists ) THEN
            IF( DownActive ) DownPointer(MaskPerm(ii)) = jmin
            IF( UpActive ) UpPointer(MaskPerm(jmin)) = ii
          ELSE
            IF( DownActive ) DownPointer(ii) = jmin
            IF( UpActive ) UpPointer(jmin) = ii              
          END IF
        END IF

      END DO
    END DO
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      
      DO i=1,nnodes
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0) CYCLE
          IF( UpActive ) THEN
            j = UpPointer(MaskPerm(i))
            IF( TopPointer(MaskPerm(i)) /= TopPointer(MaskPerm(j)) ) THEN
              UpHit = UpHit + 1
              TopPointer(MaskPerm(i)) = TopPointer(MaskPerm(j))
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(MaskPerm(i))
            IF( BotPointer(MaskPerm(i)) /= BotPointer(MaskPerm(j)) ) THEN
              DownHit = DownHit + 1
              BotPointer(MaskPerm(i)) = BotPointer(MaskPerm(j))
            END IF
          END IF
        ELSE
          IF( UpActive ) THEN
            j = UpPointer(i)
            IF( TopPointer(i) /= TopPointer(j) ) THEN
              UpHit = UpHit + 1
              TopPointer(i) = TopPointer( j )
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(i)
            IF( BotPointer(i) /= BotPointer( j ) ) THEN
              DownHit = DownHit + 1
              BotPointer(i) = BotPointer( j )
            END IF
          END IF
        END IF
      END DO
      
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO

    ! The last round is always a check
    Rounds = Rounds - 1
    
    CALL Info(Caller,'Layered structure detected in '//I2S(Rounds)//' cycles',Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF

    ! Compute the number of layers. The Rounds above may in some cases
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        EXIT
      END DO

      j = BotPointer(1)      
      CALL Info(Caller,'Starting from node: '//I2S(j),Level=15)

      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        jj = j 
        IF( MaskExists ) THEN
          jj = MaskPerm(j)
        END IF
        k = UpPointer(jj)
        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,&
          'Extruded structure layers: '//I2S(NumberOfLayers),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( NodeLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      IF( MaskExists ) THEN
        WHERE( MaskPerm == 0 ) Layer = 0
        
        DO i=1,nnodes
          IF( MaskPerm(i) == 0 ) CYCLE
          Rounds = 1
          j = BotPointer(MaskPerm(i))
          Layer(MaskPerm(j)) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(MaskPerm(j))
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(MaskPerm(j)) = Rounds
          END DO
        END DO
      ELSE        
        DO i=1,nsize
          Rounds = 1
          j = BotPointer(i)
          Layer(j) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(j)
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(j) = Rounds
          END DO
        END DO
      END IF
        
      NodeLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

    
    IF( PRESENT( MidNodePointer ) ) THEN
      ALLOCATE( MidPointer( nsize ) )
      MidPointer = 0 
      MidLayerExists = .FALSE.

      DO elem = Mesh % NumberOfBulkElements + 1, &       
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements  
        
        Element => Mesh % Elements(elem)
        NodeIndexes => Element % NodeIndexes
        
        DO bc_ind = 1, CurrentModel % NumberOfBCs 
          IF( Element % BoundaryInfo % Constraint == &
              CurrentModel % BCs(bc_ind) % Tag ) THEN
            IF( ListCheckPresent( CurrentModel % BCs(bc_ind) % Values,'Mid Surface') ) THEN
              MidPointer( NodeIndexes ) = NodeIndexes
              MidLayerExists = .TRUE.
            END IF
            EXIT
          END IF
        END DO
      END DO

      IF( MidLayerExists ) THEN
        CALL Info(Caller,'determine mid pointers',Level=15)       
                
        DO Rounds = 1, nsize
          DownHit = 0
          UpHit = 0
          DO i=1,nsize
            IF( MaskExists ) THEN
              IF( MaskPerm(i) == 0) CYCLE
            END IF

            ! We can only start from existing mid pointer
            IF( MidPointer(i) == 0 ) CYCLE
            IF( UpActive ) THEN
              j = UpPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
            IF( DownActive ) THEN
              j = DownPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF           
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
          END DO
          IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
        END DO

        CALL Info(Caller,&
            'Mid layer structure detected in '//I2S(Rounds-1)//' cycles',Level=9)
        MidNodePointer => MidPointer
      ELSE
        DEALLOCATE( MidPointer ) 
        MidNodePointer => NULL()
      END IF
    END IF

  
    ! Count the number of top and bottom nodes, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom nodes',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i) 
          IF( j == 0 ) CYCLE
          IF(TopPointer(j) == i) THEN
            MinTop = MIN( MinTop, Var % Values(j) )
            MaxTop = MAX( MaxTop, Var % Values(j) )
            TopNodes = TopNodes + 1
          END IF
        ELSE
          IF(TopPointer(i) == i) THEN
            MinTop = MIN( MinTop, Var % Values(i) )
            MaxTop = MAX( MaxTop, Var % Values(i) )
            TopNodes = TopNodes + 1
          END IF
        END IF
      END DO
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
          IF( BotPointer(j) == i) THEN
            MinBot = MIN( MinBot, Var % Values(j))
            MaxBot = MAX( MaxBot, Var % Values(j))
            BotNodes = BotNodes + 1
          END IF
        ELSE          
          IF(BotPointer(i) == i) THEN
            MinBot = MIN( MinBot, Var % Values(i))
            MaxBot = MAX( MaxBot, Var % Values(i))
            BotNodes = BotNodes + 1
          END IF
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopNodePointer ) ) THEN
        TopNodePointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpNodePointer ) ) THEN
        UpNodePointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotNodePointer ) ) THEN
        BotNodePointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownNodePointer ) ) THEN
        DownNodePointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,* ) 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)
    CALL Info(Caller,&
        'Top and bottom pointer init rounds: '//I2S(Rounds),Level=5)
    IF( UpActive ) THEN
      CALL Info(Caller,'Number of nodes at the top: '//I2S(TopNodes),Level=6)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of nodes at the bottom: '//I2S(BotNodes),Level=6)
    END IF

    IF(DownActive .AND. UpActive ) THEN
      IF(TopNodes /= BotNodes ) THEN
        CALL Fatal(Caller, 'Something wrong: top and bottom node counts differ!')
      END IF
    END IF
    

  CONTAINS
    
    
    !---------------------------------------------------------------
    SUBROUTINE CoordinateTransformationNodal( CoordTransform, R )
      CHARACTER(LEN=*) :: CoordTransform
      REAL(KIND=dp) :: R(3)
      !---------------------------------------------------------------
      REAL(KIND=dp) :: Rtmp(3)
      REAL(KIND=dp), SAVE :: Coeff 
      LOGICAL, SAVE :: Visited = .FALSE.
      

      IF( .NOT. Visited ) THEN
        IF( ListGetLogical( Params,'Angles in Degrees') ) THEN
          Coeff = 180.0_dp / PI
        ELSE
          Coeff = 1.0_dp
        END IF
        Visited = .TRUE.
      END IF
      
      SELECT CASE ( CoordTransform )
        
      CASE('cartesian to cylindrical')
        Rtmp(1) = SQRT( R(1)**2 + R(2)**2)
        Rtmp(2) = Coeff * ATAN2( R(2), R(1)  ) 
        Rtmp(3) = R(3) 
        
      CASE('cylindrical to cartesian')
        Rtmp(1) = COS( R(2) / Coeff ) * R(1)
        Rtmp(2) = SIN( R(2) / Coeff ) * R(1)
        Rtmp(3) = R(3)
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformationNodal','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT
      
      R = Rtmp

    END SUBROUTINE CoordinateTransformationNodal
   

  END SUBROUTINE DetectExtrudedStructure
 !---------------------------------------------------------------



!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh for elements.
!> Otherwise very similar as the DetectExtrudedStructure for nodes.
!> Mesh faces may need to be created in order to determine the up and down
!> pointers.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedElements( Mesh, Solver, ExtVar, &
      TopElemPointer, BotElemPointer, UpElemPointer, DownElemPointer, &
      NumberOfLayers, ElemLayer )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Solver_t) :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopElemPointer(:), BotElemPointer(:), &
        UpElemPointer(:), DownElemPointer(:)
    INTEGER, POINTER, OPTIONAL :: ElemLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(Nodes_t) :: Nodes
    TYPE(Nodes_t), POINTER :: MeshNodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
       UpHit, DownHit, bc_ind
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: UpActive, DownActive, GotIt, Found
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1
    REAL(KIND=dp) :: FaceCenter(3),FaceDx(3),Height(2),Eps, MinTop, MaxTop, MinBot, MaxBot, Diam
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    INTEGER :: TestCounter(3),ElementIndex(2)
    CHARACTER(*),PARAMETER :: Caller="DetectExtrudedElements"
         
    CALL Info(Caller,'Determining extruded element structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim

    IF( DIM /= 3 ) THEN
      CALL Fatal(Caller,'Only implemented for 3D cases: '//I2S(dim))
    END IF

    IF( .NOT. ASSOCIATED( Mesh % Faces ) ) THEN
      CALL FindMeshFaces3D( Mesh )
    END IF

    
    Params => Solver % Values
    TestCounter = 0
    
    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal(Caller,'Invalid value for Active Coordinate')
    END IF

    IF( ListGetLogical(Params,'Mapping Original Coordinates',Found ) ) THEN
      MeshNodes => Mesh % NodesOrig
    ELSE
      MeshNodes => Mesh % Nodes
    END IF
    
    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-1

    nsize = Mesh % NumberOfBulkElements
    CALL Info(Caller,'Detecting extrusion in the mesh using coordinate: '&
        //I2S(ActiveDirection),Level=8)

    IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF      

    IF( PRESENT( ExtVar ) ) ExtVar => Var

    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpElemPointer) .OR. PRESENT ( TopElemPointer ) 
    DownActive = PRESENT( DownElemPointer) .OR. PRESENT ( BotElemPointer ) 

    IF( PRESENT( NumberOfLayers) .OR. PRESENT( ElemLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nsize
        TopPointer(i) = i
        UpPointer(i) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nsize
        BotPointer(i) = i
        DownPointer(i) = i
      END DO
    END IF

    CALL Info(Caller,'determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfFaces 

      Element => Mesh % Faces(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = MeshNodes % x(NodeIndexes)
      Nodes % y(1:n) = MeshNodes % y(NodeIndexes)
      Nodes % z(1:n) = MeshNodes % z(NodeIndexes)

      IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) CYCLE
      
      FaceCenter(1) = SUM( Nodes % x(1:n) ) / n
      FaceCenter(2) = SUM( Nodes % y(1:n) ) / n
      FaceCenter(3) = SUM( Nodes % z(1:n) ) / n

      FaceDx(1) = SUM( ABS( Nodes % x(1:n) - FaceCenter(1) ) ) 
      FaceDx(2) = SUM( ABS( Nodes % y(1:n) - FaceCenter(2) ) ) 
      FaceDx(3) = SUM( ABS( Nodes % z(1:n) - FaceCenter(3) ) ) 
      
      Diam = SQRT( SUM( FaceDx**2 ) )

      ! This is not a face that separates extruded elements
      IF( FaceDx(ActiveDirection) > Eps * Diam ) CYCLE      

      TestCounter(1) = TestCounter(1) + 1      
      
      DO k = 1, 2
        IF( k == 1 ) THEN
          Parent => Element % BoundaryInfo % Left
        ELSE
          Parent => Element % BoundaryInfo % Right
        END IF
        IF( .NOT. ASSOCIATED( Parent ) ) CYCLE
               
        n = Parent % TYPE % NumberOfNodes
        NodeIndexes => Parent % NodeIndexes        

        ElementIndex(k) = Parent % ElementIndex
        Height(k) = SUM( Var % Values(NodeIndexes) ) / n
      END DO      

      IF( Height(1) > Height(2) ) THEN
        IF( UpActive ) UpPointer(ElementIndex(2)) = ElementIndex(1)
        IF( DownActive ) DownPointer(ElementIndex(1)) = ElementIndex(2)
      ELSE
        IF( UpActive ) UpPointer(ElementIndex(1)) = ElementIndex(2)
        IF( DownActive ) DownPointer(ElementIndex(2)) = ElementIndex(1)
      END IF
    END DO  
        
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      DO i=1,nsize
        IF( UpActive ) THEN
          j = UpPointer(i)
          IF( TopPointer(i) /= TopPointer( j ) ) THEN
            UpHit = UpHit + 1
            TopPointer(i) = TopPointer( j )
          END IF
        END IF
        IF( DownActive ) THEN
          j = DownPointer(i)
          IF( BotPointer(i) /= BotPointer( j ) ) THEN
            DownHit = DownHit + 1
            BotPointer(i) = BotPointer( j )
          END IF
        END IF
      END DO
      CALL Info(Caller,'Hits in determining structure: '//I2S(UpHit+DownHit),Level=10)
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO
    ! The last round is always a check
    Rounds = Rounds - 1


    WRITE( Message,'(A,I0,A)') 'Layered elements detected in ',Rounds,' cycles'
    CALL Info(Caller,Message,Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF


    ! Compute the number of layers. The Rounds above may in some cases 
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    

      ! We start from any bottom row entry
      j = BotPointer(1)
      
      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        k = UpPointer(j)

        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO      

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,'Extruded structure layers: '//I2S(NumberOfLayers),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( ElemLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      
      DO i=1,nsize
        Rounds = 1
        j = BotPointer(i)
        Layer(j) = Rounds
        DO WHILE(.TRUE.)
          k = UpPointer(j)
          IF( k == j ) EXIT          
          Rounds = Rounds + 1
          j = k
          Layer(j) = Rounds
        END DO
      END DO
      
      ElemLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

  
    ! Count the number of top and bottom elements, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom elements',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nsize
        IF(TopPointer(i) == i) THEN
          MinTop = MIN( MinTop, Var % Values(i) )
          MaxTop = MAX( MaxTop, Var % Values(i) )
          TopNodes = TopNodes + 1
        END IF
      END DO
      CALL Info(Caller,'Number of top elements: '//I2S(TopNodes),Level=9)
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nsize
        IF(BotPointer(i) == i) THEN
          MinBot = MIN( MinBot, Var % Values(i))
          MaxBot = MAX( MaxBot, Var % Values(i))
          BotNodes = BotNodes + 1
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopElemPointer ) ) THEN
        TopElemPointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpElemPointer ) ) THEN
        UpElemPointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotElemPointer ) ) THEN
        BotElemPointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownElemPointer ) ) THEN
        DownElemPointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,'(A,ES12.3)') 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)

    CALL Info(Caller,'Top and bottom pointer init rounds: '//I2S(Rounds),Level=8)

    IF( UpActive ) THEN
      CALL Info(Caller,'Number of elements at the top: '//I2S(TopNodes),Level=8)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of elements at the bottom: '//I2S(BotNodes),Level=8)
    END IF

    IF(DownActive .AND. UpActive ) THEN
      IF(TopNodes /= BotNodes ) THEN
        CALL Fatal(Caller, 'Something wrong: top and bottom element counts differ!')
      END IF
    END IF
    

  END SUBROUTINE DetectExtrudedElements
 !---------------------------------------------------------------


  SUBROUTINE StoreOriginalCoordinates(Mesh)
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp), POINTER CONTIG :: NewCoords(:)
    INTEGER :: n

    IF( ASSOCIATED( Mesh % NodesOrig ) ) THEN
      CALL Info('StoreOriginalCoordinates','Original coordinates already stored')
    END IF

    n = SIZE( Mesh % Nodes % x )    
    ALLOCATE( NewCoords(3*n) )

    ALLOCATE( Mesh % NodesOrig ) 
    Mesh % NodesOrig % x => NewCoords(1:n)
    Mesh % NodesOrig % y => NewCoords(n+1:2*n)
    Mesh % NodesOrig % z => NewCoords(2*n+1:3*n)

    Mesh % NodesOrig % x = Mesh % Nodes % x
    Mesh % NodesOrig % y = Mesh % Nodes % y
    Mesh % NodesOrig % z = Mesh % Nodes % z

    Mesh % NodesMapped => Mesh % Nodes

    CALL Info('StoreOriginalCoordinates','Original coordinates stored',Level=6)
    
  END SUBROUTINE StoreOriginalCoordinates

    
   
  !----------------------------------------------------------------
  !> Maps coordinates from the original nodes into a new coordinate
  !> system while optionally maintaining the original coordinates. 
  !> Note that this may be called 
  !---------------------------------------------------------------
  SUBROUTINE CoordinateTransformation( Mesh, CoordTransform, Params, &
      IrreversibleTransformation )
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: CoordTransform
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL, OPTIONAL :: IrreversibleTransformation
    !---------------------------------------------------------------   
    REAL(KIND=dp) :: R0(3),R1(3),Coeff,Rad0
    LOGICAL :: Irreversible,FirstTime,Reuse,UpdateNodes,Found
    REAL(KIND=dp), POINTER :: x0(:),y0(:),z0(:),x1(:),y1(:),z1(:)
    REAL(KIND=dp), POINTER CONTIG :: NewCoords(:)
    INTEGER :: i,j,k,n,Mode
    TYPE(Variable_t), POINTER :: Var

    ! The coordinate transformation may either be global for all the solvers
    ! and this overrides the original nodes permanently. 
    ! Or it can be a solver specific transformation which saves the initial 
    ! coordinates. 
    CALL Info('CoordinateTransformation','Starting')

    IF(.NOT. ASSOCIATED(Mesh) ) THEN
      CALL Fatal('CoordinateTransformation','Mesh not associated!')
    END IF

    IF( PRESENT( IrreversibleTransformation ) ) THEN
      Irreversible = IrreversibleTransformation
    ELSE
      Irreversible = .FALSE.
    END IF

    n = Mesh % NumberOfNodes 

    x0 => Mesh % Nodes % x
    y0 => Mesh % Nodes % y
    z0 => Mesh % Nodes % z
    
    IF( Irreversible ) THEN
      UpdateNodes = .TRUE.
      ! Map to the same nodes
      x1 => Mesh % Nodes % x
      y1 => Mesh % Nodes % y
      z1 => Mesh % Nodes % z
    ELSE
      ReUse = ListGetLogical(Params,'Coordinate Transformation Reuse',Found ) 
      FirstTime = .NOT. ASSOCIATED( Mesh % NodesMapped )
      IF( FirstTime ) THEN
        ALLOCATE( Mesh % NodesMapped )
        NULLIFY( NewCoords )
        ALLOCATE( NewCoords(3*n) )
        NewCoords = 0.0_dp
        Mesh % NodesMapped % x => NewCoords(1:n)
        Mesh % NodesMapped % y => NewCoords(n+1:2*n)
        Mesh % NodesMapped % z => NewCoords(2*n+1:3*n)
        ! Mesh % NodesMapped % x => NewCoords(1::3)
        ! Mesh % NodesMapped % y => NewCoords(2::3)
        ! Mesh % NodesMapped % z => NewCoords(3::3)
      ELSE
        IF( n /= SIZE(Mesh % NodesMapped % x) ) THEN
          CALL Fatal('CoordinateTransformation','Sizes of original and mapped mesh differ!')
        END IF
      END IF

      IF( CoordTransform == 'previous' ) THEN
        IF( FirstTime ) THEN
          CALL Fatal('CoordinateTransformation','One cannot reuse unexisting transformation!')
        END IF
        ReUse = .TRUE.
      END IF

      ! Note that if many solvers reutilize the same coordinates then they must 
      ! also have the same coordinate mapping. 
      !------------------------------------------------------------------------
      UpdateNodes = FirstTime .OR. .NOT. ReUse 
      ! Map different nodes if the original ones are kept
      x1 => Mesh % NodesMapped % x
      y1 => Mesh % NodesMapped % y
      z1 => Mesh % NodesMapped % z      

      IF( FirstTime ) THEN
        IF( ListGetLogical(Params,'Coordinate Transformation Save',Found ) ) THEN
          CALL Info('CoordinateTranformation',&
              'Creating variables for > Transformed Coordinate < ')
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 1',1,x1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 2',1,y1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 3',1,z1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate',3,NewCoords)
        END IF
      END IF
    END IF
      
    IF( UpdateNodes ) THEN
      IF( ListGetLogical( Params,'Coordinate Transformation Use Degrees',Found) ) THEN
        Coeff = 180.0_dp / PI
        CALL Info('CoordinateTranformation','Using degrees for angles')
      ELSE
        Coeff = 1.0_dp
      END IF

      Rad0 = ListGetConstReal( Params,'Coordinate Transformation Radius',Found )
  
      SELECT CASE ( CoordTransform ) 
        
      CASE('cartesian to polar')
        Mode = 1
      CASE('cartesian to cylindrical')
        Mode = 1
      CASE('polar to cartesian')
        Mode = -1
      CASE('cylindrical to cartesian')
        Mode = -1
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformation','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT

      DO i=1,n    
        R0(1) = x0(i)
        R0(2) = y0(i)
        R0(3) = z0(i)
        
        IF( Mode == 1 ) THEN
          R1(1) = Rad0 + SQRT( R0(1)**2 + R0(2)**2)
          R1(2) = Coeff * ATAN2( R0(2), R0(1)  ) 
          R1(3) = R0(3)    
       
        ELSE IF( Mode == -1 ) THEN
          R1(1) = COS( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(2) = SIN( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(3) = R0(3)          
        END IF

        x1(i) = R1(1)
        y1(i) = R1(2)
        z1(i) = R1(3)

      END DO
    END IF

    IF( .NOT. Irreversible ) THEN
      Mesh % NodesOrig => Mesh % Nodes
      Mesh % Nodes => Mesh % NodesMapped

      Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
      Var % Values => Mesh % Nodes % x

      Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
      Var % Values => Mesh % Nodes % y

      Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
      Var % Values => Mesh % Nodes % z
    END IF

    CALL Info('CoordinateTransformation','All done',Level=12)

  END SUBROUTINE CoordinateTransformation
!---------------------------------------------------------------



!---------------------------------------------------------------
!> Perform rigid mesh mapping given by body force section translate,
!> Rotate and scale. This is copy-pasted from RigidMeshMapper solver
!> and intended use is when the solver cannot be used i.e. in
!> conjunction of view factor computation. 
!---------------------------------------------------------------
  SUBROUTINE RigidMeshMapping( Model, Mesh, VFMode ) 

    TYPE(Model_t) :: Model
    TYPE(Mesh_t) :: Mesh
    LOGICAL :: VFMode
    !---------------------------------------------------------------
    TYPE(ValueList_t),POINTER :: SolverParams
    INTEGER :: i,j,k,n,m,t,dim,elem,bf_id,istat,node,nodei,RotateOrder(3)
    INTEGER, TARGET :: CurrentNode(1)
    INTEGER :: NonlinIter, MaxNonlinIter, NoNodes
    INTEGER, POINTER :: NodeIndex(:), IntArray(:) => NULL()
    REAL(KIND=dp) :: x0(4), x1(4), RotMatrix(4,4), TrsMatrix(4,4), SclMatrix(4,4), &
        TrfMatrix(4,4), Identity(4,4), Origin(4), Angles(3), Scaling(3), alpha, Coord(3), &
        dCoord(3), Norm, dx(3) 
    REAL(KIND=dp) :: at0,at1,at2,Coeff,Source,relax(1),MaxDeform,AngleCoeff
    REAL(KIND=dp), POINTER :: X(:),Y(:),Z(:),PArray(:,:) => NULL()
    LOGICAL :: Found,GotMatrix,GotRotate,GotTranslate,GotScale,Visited=.FALSE.,&
        TranslateBeforeRotate, WholeMode
    LOGICAL :: AnyMeshMatrix,AnyMeshRotate,AnyMeshTranslate,AnyMeshScale,&
        AnyMeshOrigin, AnyRelax, ConstantMap, GotMap
    LOGICAL, POINTER :: NodeDone(:), FixedNode(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t),POINTER :: ValueList, PrevValueList     
    CHARACTER(LEN=MAX_NAME_LEN) :: sname
    CHARACTER(*), PARAMETER :: Caller = 'RigidMeshMapping'


    CALL Info( Caller,'---------------------------------------',Level=4 )
    CALL Info( Caller,'Performing analytic mesh mapping ',Level=4 )
    CALL Info( Caller,'---------------------------------------',Level=4 )

    at0 = CPUTime()
       
    ! First check if this is called by view factors and we use the definitions from a solver. 
    SolverParams => NULL()
    IF( VFMode ) THEN
      ! In view factor mode we usually expect the rigid mesh mapping to be given by a separate code
      ! that holds definitions for how to perform the mapping. 
      DO i=Model % NumberOfSolvers,1,-1
        sname = ListGetString(Model % Solvers(i) % Values, 'Procedure', Found)      
        j = INDEX( sname,'RigidMeshMapper') 
        IF( j > 0 ) THEN
          CALL Info(Caller,'Using rigid body mapping definitions from Solver '//I2S(i),Level=15)
          SolverParams => Model % Solvers(i) % Values
          EXIT
        END IF
      END DO
    END IF
      
    ! Do we have rigid mesh mapping requested ins simulation section.
    IF(.NOT. ASSOCIATED(SolverParams)) THEN
      SolverParams => Model % Simulation
      CALL Info(Caller,'Using rigid body mapping from Simulation section',Level=15)
    END IF
    
    IF( ListGetLogical( Model % Simulation,'Rotate in Radians',Found ) ) THEN
      AngleCoeff = 1.0_dp
    ELSE
      AngleCoeff = PI / 180.0_dp
    END IF

    TranslateBeforeRotate = ListGetLogical( SolverParams,&
        'Translate Before Rotate',Found )

    RotateOrder = (/1, 2, 3/)
    IntArray => ListGetIntegerArray( SolverParams,'Mesh Rotation Axis Order', Found)
    IF(Found) THEN
      j = SIZE(IntArray)
      IF(j /= dim) THEN
        CALL Fatal(Caller,'Size of Mesh Rotation Axis Order must match dimension.')
      END IF
      DO i=1,j
        RotateOrder(i) = IntArray(j+1-i) !reverse the order
      END DO
    END IF

    WholeMode = ListGetLogical( SolverParams,'Whole Mesh Mode',Found ) 
    IF( WholeMode ) THEN
      CALL Info(Caller,'Moving the whole mesh with keywords from same section!')
    END IF

    IF( WholeMode ) THEN
      ValueList => SolverParams      
      AnyMeshMatrix = ListCheckPresent( ValueList,'Mesh Matrix')   
      AnyMeshRotate = ListCheckPrefix( ValueList,'Mesh Rotate')
      AnyMeshTranslate = ListCheckPrefix( ValueList,'Mesh Translate')
      AnyMeshScale = ListCheckPrefix( ValueList,'Mesh Scale') 
      AnyMeshOrigin = ListCheckPrefix( ValueList,'Mesh Origin')
    ELSE
      AnyMeshMatrix = ListCheckPresentAnyBodyForce( Model,'Mesh Matrix')   
      AnyMeshRotate = ListCheckPrefixAnyBodyForce( Model,'Mesh Rotate')
      AnyMeshTranslate = ListCheckPrefixAnyBodyForce( Model,'Mesh Translate')
      AnyMeshScale = ListCheckPrefixAnyBodyForce( Model,'Mesh Scale') 
      AnyMeshOrigin = ListCheckPrefixAnyBodyForce( Model,'Mesh Origin')
    END IF

    NoNodes = Mesh % NumberOfNodes
    ALLOCATE( NodeDone(NoNodes), FixedNode(NoNodes))
    NodeDone = .FALSE.
    FixedNode = .FALSE.
    NodeIndex => CurrentNode

    Identity = 0.0d0
    DO i=1,4
      Identity(i,i) = 1.0d0
    END DO    
        
    X => Mesh % Nodes % x
    Y => Mesh % Nodes % y
    Z => Mesh % Nodes % z

    ! These are used to save time in creating the projection matrices.
    PrevValueList => NULL()
    GotMap = .FALSE.
    ConstantMap = ListGetLogical( SolverParams,'Constant Mapping',Found ) 

    GotRotate = .FALSE.
    GotTranslate = .FALSE.
    GotScale = .FALSE.
    GotMatrix = .FALSE.



    ! Implement moving and fixed BCs
    ! ------------------------------
    DO elem = 1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Element => Mesh % Elements(elem)
      Model % CurrentElement => Element
      n = Element % Type % NumberOfNodes

      IF( Element % BodyId > 0 ) THEN
        i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found )
        IF(Found) THEN
          ValueList => Model % BodyForces(i) % Values
          IF(ListGetLogical(ValueList,'Fixed Body',Found ) ) THEN
            FixedNode(Element % NodeIndexes) = .TRUE.
          END IF
        END IF
      END IF
        
      IF(elem > Mesh % NumberOfBulkElements ) THEN        
        DO i = 1, CurrentModel % NumberOfBCs 
          IF( Element % BoundaryInfo % Constraint == CurrentModel % BCs(i) % Tag ) THEN
            IF( ListGetLogical( CurrentModel % BCs(i) % Values,'Fixed Boundary',Found ) ) THEN
              FixedNode(Element % NodeIndexes) = .TRUE.
            END IF
            EXIT
          END IF            
        END DO
      END IF
    END DO

    n = COUNT(FixedNode) 
    IF(n>0) CALL Info(Caller,'Number of Fixed nodes: '//I2S(n))
    
    
    DO elem = 1,Mesh % NumberOfBulkElements      

      Element => Mesh % Elements(elem)
      Model % CurrentElement => Element
      n = Element % Type % NumberOfNodes

      IF( WholeMode ) THEN
        ! If we are doing the whole mesh then do all elements. 
        CONTINUE
      ELSE
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found )
        IF(.NOT. Found) CYCLE
        ValueList => Model % BodyForces(bf_id) % Values
        IF( ConstantMap ) THEN
          GotMap = ASSOCIATED( ValueList, PrevValueList )
          PrevValueList => ValueList
        END IF
      END IF

      DO Node=1,n
        NodeIndex(1) = Element % NodeIndexes(Node)
        NodeI = NodeIndex(1)

        IF(NodeDone(NodeI)) CYCLE

        IF(FixedNode(NodeI)) CYCLE 
        
        ! This is to save time. If we have exactly same mapping as last time then
        ! there is no use doing the same ListGet operation things again.
        !-------------------------------------------------------------------------
        IF( GotMap ) GOTO 100
        IF( WholeMode ) GotMap = .TRUE.

        ! Generic transformation matrix
        !--------------------------------
        GotMatrix = .FALSE.
        IF( AnyMeshMatrix ) THEN
          Parray => ListGetConstRealArray( ValueList,'Mesh Matrix', GotMatrix )
          IF ( GotMatrix ) THEN
            DO i=1,SIZE(Parray,1)
              DO j=1,SIZE(Parray,2)
                TrfMatrix(i,j) = Parray(j,i)
              END DO
            END DO
          END IF
        END IF

        IF(.NOT. GotMatrix ) THEN
          ! Rotations around main axis:
          !----------------------------        
          GotRotate = .FALSE.
          Angles = 0.0_dp

          IF( AnyMeshRotate ) THEN
            Parray => ListGetConstRealArray( ValueList,'Mesh Rotate', GotRotate )                
            IF ( GotRotate ) THEN
              DO i=1,SIZE(Parray,1)
                Angles(i) = Parray(i,1) 
              END DO
            ELSE 
              Angles(1:1) = ListGetReal( ValueList,'Mesh Rotate 1', 1, NodeIndex, Found )
              IF( Found ) GotRotate = .TRUE.
              Angles(2:2) = ListGetReal( ValueList,'Mesh Rotate 2', 1, NodeIndex, Found )
              IF( Found ) GotRotate = .TRUE.
              Angles(3:3) = ListGetReal( ValueList,'Mesh Rotate 3', 1, NodeIndex, Found )
              IF( Found ) GotRotate = .TRUE.
            END IF
            Angles = AngleCoeff * Angles
          END IF

          ! Scaling:
          !---------
          GotScale = .FALSE.
          IF( AnyMeshScale ) THEN
            Parray => ListGetConstRealArray( ValueList,'Mesh Scale', GotScale )
            IF ( GotScale ) THEN
              Scaling = 0.0_dp
              DO i=1,SIZE(Parray,1)
                Scaling(i) = Parray(i,1)
              END DO
            ELSE 
              Scaling(1:1) = ListGetReal( ValueList,'Mesh Scale 1',1,NodeIndex,GotScale) 
              IF(.NOT. GotScale ) Scaling(1) = 1.0_dp
              Scaling(2:2) = ListGetReal( ValueList,'Mesh Scale 2',1,NodeIndex,Found) 
              IF(.NOT. Found ) Scaling(2) = 1.0_dp
              GotScale = GotScale .OR. Found
              Scaling(3:3) = ListGetReal( ValueList,'Mesh Scale 3',1,NodeIndex,Found) 
              IF(.NOT. Found ) Scaling(3) = 1.0_dp
              GotScale = GotScale .OR. Found
            END IF
          END IF

          ! Translations:
          !---------------
          GotTranslate = .FALSE.
          IF( AnyMeshTranslate ) THEN
            Parray => ListGetConstRealArray( ValueList,'Mesh Translate', GotTranslate )
            IF ( GotTranslate ) THEN
              dCoord = 0.0_dp
              DO i=1,SIZE(Parray,1)
                dCoord(i) = Parray(i,1)
              END DO
            ELSE 
              dCoord(1:1) = ListGetReal( ValueList,'Mesh Translate 1', 1, NodeIndex, GotTranslate) 
              dCoord(2:2) = ListGetReal( ValueList,'Mesh Translate 2', 1, NodeIndex, Found) 
              GotTranslate = GotTranslate .OR. Found
              dCoord(3:3) = ListGetReal( ValueList,'Mesh Translate 3', 1, NodeIndex, Found) 
              GotTranslate = GotTranslate .OR. Found
            END IF
          END IF

          GotMatrix = GotRotate .OR. GotScale

          IF(GotMatrix) THEN
            TrsMatrix = Identity
            SclMatrix = Identity

            ! Origin:
            !---------
            IF( GotRotate ) THEN
              RotMatrix = Identity

              DO i=1,3
                j = RotateOrder(i)
                Alpha = Angles(j) 

                IF( ABS(Alpha) < TINY(Alpha) ) CYCLE
                TrfMatrix = Identity

                SELECT CASE(j)
                CASE(1)
                  TrfMatrix(2,2) =  COS(Alpha)
                  TrfMatrix(2,3) = -SIN(Alpha)
                  TrfMatrix(3,2) =  SIN(Alpha)
                  TrfMatrix(3,3) =  COS(Alpha)
                CASE(2)
                  TrfMatrix(1,1) =  COS(Alpha)
                  TrfMatrix(1,3) = -SIN(Alpha)
                  TrfMatrix(3,1) =  SIN(Alpha)
                  TrfMatrix(3,3) =  COS(Alpha)
                CASE(3)
                  TrfMatrix(1,1) =  COS(Alpha)
                  TrfMatrix(1,2) = -SIN(Alpha)
                  TrfMatrix(2,1) =  SIN(Alpha)
                  TrfMatrix(2,2) =  COS(Alpha)
                END SELECT
                RotMatrix = MATMUL( RotMatrix, TrfMatrix )
              END DO
            END IF

            IF( GotTranslate ) THEN
              DO i=1,3
                TrsMatrix(i,4) = dCoord(i)
              END DO
            END IF

            ! It may be easier to first translate the matrix to origin 
            ! and only the do the rotation than vice versa. 
            IF( TranslateBeforeRotate ) THEN
              TrfMatrix = MATMUL( RotMatrix, TrsMatrix )
            ELSE
              TrfMatrix = MATMUL( TrsMatrix, RotMatrix )
            END IF

            IF( GotScale ) THEN
              DO i=1,3
                SclMatrix(i,i) = Scaling(i)
              END DO
              TrsMatrix = TrfMatrix
              TrfMatrix = MATMUL( SClMatrix, TrsMatrix )
            END IF
          END IF
        END IF

        ! Get mesh origin
        !----------------------------------------------------
        Origin = 0.0_dp
        IF( GotMatrix .AND. AnyMeshOrigin ) THEN
          Parray => ListGetConstRealArray( ValueList,'Mesh Origin', Found )
          IF ( Found ) THEN
            DO i=1,SIZE(Parray,1)
              Origin(i) = Parray(i,1)
            END DO
          ELSE
            Origin(1:1) = ListGetReal( ValueList,'Mesh Origin 1', 1, NodeIndex, Found) 
            Origin(2:2) = ListGetReal( ValueList,'Mesh Origin 2', 1, NodeIndex, Found) 
            Origin(3:3) = ListGetReal( ValueList,'Mesh Origin 3', 1, NodeIndex, Found) 
          END IF
        END IF

100     IF( GotMatrix ) THEN                         
          x0(1) = X(NodeI)
          x0(2) = Y(NodeI)
          x0(3) = Z(NodeI)
          x0(4) = 1.0_dp
          x1 = MATMUL( TrfMatrix, x0 - Origin ) + Origin          
          dx(1:3) = x1(1:3) / x1(4) - x0(1:3)
        ELSE IF( GotTranslate ) THEN        
          dx(1:3) = dCoord(1:3)
        ELSE
          CYCLE
        END IF

        ! Now do the actual mapping!
        X(NodeI) = dx(1) + X(NodeI)
        Y(NodeI) = dx(2) + Y(NodeI)
        Z(NodeI) = dx(3) + Z(NodeI)

        ! Mark the done so that we do not revisit it again. 
        NodeDone(NodeI) = .TRUE.
      END DO
    END DO
    
    CALL Info(Caller,'Number of nodes mapped: '//I2S(COUNT(NodeDone)),Level=7)
    at1 = CPUTime()
    WRITE(Message,* ) 'Coordinate mapping time: ',at1-at0
    CALL Info(Caller,Message,Level=7)

    DEALLOCATE( NodeDone, FixedNode )

    CALL Info(Caller,'Performed internal rigid mesh mapping!',Level=10)
    
  END SUBROUTINE RigidMeshMapping
    
  

!---------------------------------------------------------------
!> Return back to the original coordinate system. 
!---------------------------------------------------------------
  SUBROUTINE BackCoordinateTransformation( Mesh, DeleteTemporalMesh )
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: DeleteTemporalMesh
!---------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var

    IF( PRESENT( DeleteTemporalMesh ) ) THEN
      IF( DeleteTemporalMesh ) THEN
        DEALLOCATE( Mesh % NodesMapped % x, &
            Mesh % NodesMapped % y, &
            Mesh % NodesMapped % z ) 
        DEALLOCATE( Mesh % NodesMapped )
      END IF
    END IF

    IF( .NOT. ASSOCIATED( Mesh % NodesOrig ) ) THEN
      CALL Fatal('BackCoordinateTransformation','NodesOrig not associated')
    END IF

    Mesh % Nodes => Mesh % NodesOrig

    Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
    Var % Values => Mesh % Nodes % x
    
    Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
    Var % Values => Mesh % Nodes % y

    Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
    Var % Values => Mesh % Nodes % z

  END SUBROUTINE BackCoordinateTransformation


END MODULE MeshTransform
