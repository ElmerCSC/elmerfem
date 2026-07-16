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
! *  Modified Data: 27.9.2012
! *
! *****************************************************************************/


!> \ingroup Solvers
!> \{

SUBROUTINE StructuredProjectToPlane_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Interpolation

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
  INTEGER :: NormInd
  LOGICAL :: GotIt

  Params => GetSolverParams()
  
  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  NormInd = ListGetInteger( Params,'Show Norm Index',GotIt)
  IF( NormInd > 0 ) THEN
    IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
      CALL ListAddString( Params,'Variable','-nooutput -global savescalars_var')
    END IF
  END IF
  
  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)
  
END SUBROUTINE StructuredProjectToPlane_init


!------------------------------------------------------------------------------
!> Subroutine for projecting results in structured 3d mesh to a 2d surface.
!>  This solver assumes that the mesh is structural so that it could have 
!>  been obtained by extrusion in the direction of interest. For the given 
!>  direction the corresponding top and bottom node is computed for every node
!>  and this information is used to perform projection to the top or bottom
!>  plane, or alternatively to the whole body. 
!------------------------------------------------------------------------------
SUBROUTINE StructuredProjectToPlane( Model,Solver,dt,Transient )
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
  TYPE(ValueList_t),POINTER :: Params
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, OldVarName, Name, Oper, Oper0, OldOper, &
      LevelsetName, TargetName
  INTEGER :: t,i,j,k,l,kk,ll,n,dim,Dofs,dof,itop,ibot,idown,iup,jup,lup,ii,jj,nsize,layer, &
      ActiveDirection,elem,TopNodes,MidNodes,NoVar,BotNodes, NormInd, nnodes, TypeIn, &
      NoLayers, rDofs, DofsIn
  INTEGER, POINTER :: MaskPerm(:),TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:),&
      MidPointer(:), NodeIndexes(:),TargetPointer(:),BotPerm(:),&
      PermOut(:),PermIn(:),LevelsetPerm(:),TopPerm(:),MidPerm(:),UnitPerm(:)=>NULL(), &
      TmpTopPointer(:), TmpBotPointer(:)
  INTEGER, ALLOCATABLE, TARGET :: InvDGPerm(:)
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
      MaskExist, GotVar, GotOldVar, GotOper, BottomTarget, ReducedDimensional, &
      MidLayerExists, UpperOper, LowerOper, ProjectEverywhere
  REAL(KIND=dp) :: dx,UnitVector(3),ElemVector(3),DotPro,Eps,Length,Level,val,q,depth,height
  REAL(KIND=dp) :: at0,at1,at2
  REAL(KIND=dp), POINTER :: FieldOut(:), FieldIn(:), Levelset(:), Coord(:),TopField(:)
  TYPE(Variable_t), POINTER :: Var, OldVar, RefVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t),POINTER :: BC
  CHARACTER(*), PARAMETER :: Caller = 'StructuredProjectToPlane'

  
  SAVE Visited,Nodes,Initialized,UnitVector,Coord,MaskExist,MaskPerm,TopPointer,&
      BotPointer,MidPointer, UpPointer,DownPointer,FieldOut,FieldIn,&
      TopNodes,MidNodes,TopPerm, MidPerm, TopField, BotNodes, BotPerm, nsize, &
      nnodes, UnitPerm, MidLayerExists, NoLayers
 
  CALL Info( Caller,'------------------------------------------',Level=4 )
  CALL Info( Caller,'Performing projection on a structured mesh ',Level=4 )
  CALL Info( Caller,'------------------------------------------',Level=4 )

!------------------------------------------------------------------------------
!   Initialize the pointers to top and bottom nodes 
!------------------------------------------------------------------------------

  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver


  IF( .NOT. Initialized ) THEN
    at0 = CPUTime()

    ! Choose active direction coordinate and set corresponding unit vector
    !---------------------------------------------------------------------
    PSolver => Solver
    CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar = Var, &
        TopNodePointer = TopPointer, BotNodePointer = BotPointer, &
        UpNodePointer = UpPointer, DownNodePointer = DownPointer, &
        MidNodePointer = MidPointer, MidLayerExists = MidLayerExists )
    MaskExist = ASSOCIATED( Var % Perm ) 
    IF( MaskExist ) MaskPerm => Var % Perm
    Coord => Var % Values
    nsize = MIN( SIZE( Coord ), Mesh % NumberOfNodes )
    nnodes = Mesh % NumberOfNodes
    Initialized = .TRUE.

    TopNodes = 0
    ALLOCATE( TopPerm( Mesh % NumberOfNodes ) )
    TopPerm = 0
    DO i=1,Mesh % NumberOfNodes
      j = i
      IF( MaskExist ) THEN
        j = MaskPerm(i)
        IF( j == 0 ) CYCLE
      END IF       
      IF(TopPointer(j) == i) THEN
        TopNodes = TopNodes + 1
        TopPerm(i) = TopNodes
      END IF
    END DO
    IF( TopNodes > 0 ) THEN
      ALLOCATE( TopField( TopNodes ) ) 
      TopField = 0.0_dp
    END IF
    CALL Info(Caller,'Number of top nodes: '//I2S(TopNodes),Level=10)

    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    NoLayers = 0
    DO i=1,Mesh % NumberOfNodes
      j = i
      IF( MaskExist ) THEN
        j = MaskPerm(i)
        IF( j == 0 ) CYCLE
      END IF
      IF(BotPointer(j) == i) THEN
        BotNodes = BotNodes + 1
        BotPerm(i) = BotNodes
      END IF
      NoLayers = NoLayers + 1
    END DO
    CALL Info(Caller,'Number of bot nodes: '//I2S(BotNodes),Level=10)
    IF(BotNodes == 0) CALL Fatal(Caller,'Cannot continue with zero BotNodes!')
    
    NoLayers = NoLayers / BotNodes 
    CALL Info(Caller,'Number of node layers: '//I2S(NoLayers),Level=10)
    IF(NoLayers < 2) THEN
      CALL Fatal(Caller,'Solver does not makse sense with '//I2S(NoLayers)//' layers!')
    END IF
    
    CALL Info(Caller,'Number of bot nodes: '//I2S(BotNodes),Level=10)
    IF(BotNodes /= TopNodes) CALL Warn(Caller,'Conflicting BotNodes vs. TopNodes: '&
        //I2S(BotNodes)//' - '//I2S(TopNodes))
    
    IF( MidLayerExists ) THEN
      MidNodes = 0
      ALLOCATE( MidPerm( Mesh % NumberOfNodes ) ) 
      MidPerm = 0
      DO i=1,Mesh % NumberOfNodes
        j = i
        IF( MaskExist ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
        END IF
        IF(MidPointer(j) == i) THEN
          MidNodes = MidNodes + 1
          MidPerm(i) = MidNodes
        END IF
      END DO
      CALL Info(Caller,'Number of mid nodes: '//I2S(MidNodes),Level=10)
      IF(BotNodes /= MidNodes) CALL Warn(Caller,'Conflicting BotNodes vs. MidNodes: '&
          //I2S(BotNodes)//' - '//I2S(MidNodes))
    END IF
  END IF
  at0 = CPUTime()

  !------------------------------------------------------------------------------
  ! Go through the variables and compute the desired projections
  !------------------------------------------------------------------------------
  GotVar  = .TRUE.
  GotOldVar = .FALSE.
  GotOper = .FALSE.
  NULLIFY(OldVar)
  NoVar = 0

  NormInd = ListGetInteger( Params,'Show Norm Index',GotIt)

  DO WHILE(.TRUE.)

    NoVar = NoVar + 1    

    WRITE (Name,'(A,I0)') 'Variable ',NoVar
    VarName = ListGetString( Params, TRIM(Name), GotVar )
    NULLIFY(Var)    
    IF(GotVar) THEN
      Var => VariableGet( Model % Variables, TRIM(VarName) )
      IF ( .NOT. ASSOCIATED( Var ) )  THEN
        CALL Fatal(Caller,'Variable does not exist: '//TRIM(VarName))
      END IF
      FieldIn => Var % Values
      DofsIn = Var % Dofs
      PermIn => Var % Perm
      TypeIn = Var % TYPE

      IF( TypeIn == Variable_on_nodes_on_elements ) THEN
        IF(.NOT. ALLOCATED( InvDGPerm ) ) THEN
          ALLOCATE( InvDGPerm( Mesh % NumberOfNodes ) )
        END IF
        InvDGPerm = 0
        DO t=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(t)
          n = Element % TYPE % NumberOfNodes

          IF( MaskExist ) THEN
            IF( ANY( MaskPerm( Element % NodeIndexes ) == 0 ) ) CYCLE
          END IF
          DO i=1,n
            k = Element % NodeIndexes(i)
            j = Element % DGIndexes(i)
            IF( InvDGPerm(k) == 0 ) THEN
              InvDGPerm(k) = PermIn(j)
            ELSE IF( InvDGPerm(k) /= PermIn(j) ) THEN
              CALL Fatal(Caller,'The DG field is not a bijection!')
            END IF
          END DO
        END DO
        PermIn => InvDGPerm
      END IF
      
      Dofs = Var % Dofs
      GotOldVar = .TRUE.
      OldVarName = VarName
    ELSE
      Dofs = 1
    END IF
   
    ! Read in the operator
    !-----------------------------------------------
    WRITE (Name,'(A,I0)') 'Operator ',NoVar
    Oper = ListGetString( Params, TRIM(Name),GotOper)
    IF(.NOT. GotOper ) THEN
      ! For the first variable the operator is definitely needed
      IF( NoVar > 1 ) THEN
        Oper = OldOper
      END IF
    END IF
 
    
    ! Either new field or new operator is needed
    !-----------------------------------------------
    IF( .NOT. (GotVar .OR. GotOper ) ) THEN
      IF( NoVar == 1 ) THEN
        CALL Warn(Caller,'Not even one field and operator to treat?')
      END IF
      EXIT
    END IF
    OldOper = Oper
    IF(GotOper .AND. .NOT. GotVar ) THEN
      GotVar = GotOldVar
      IF(GotVar) VarName = OldVarName
    END IF

    Oper0 = Oper
    UpperOper = .FALSE.
    LowerOper = .FALSE.
    TmpTopPointer => TopPointer
    TmpBotPointer => BotPointer


    IF( Oper0(1:6) == 'upper ' ) THEN
      IF( .NOT. MidLayerExists ) THEN
        CALL Fatal(Caller,'Upper operator cannot exist without midlayer')
      END IF
      UpperOper = .TRUE.
      Oper = TRIM(Oper0(7:))
      TmpBotPointer => MidPointer
      CALL Info(Caller,'Operating on the upper part with: '//TRIM(Oper),Level=10)
    ELSE IF( Oper0(1:6) == 'lower ') THEN
      IF( .NOT. MidLayerExists ) THEN
        CALL Fatal(Caller,'Lower operator cannot exist without midlayer')
      END IF
      LowerOper = .TRUE.
      Oper = TRIM(Oper0(7:))
      TmpTopPointer => MidPointer     
      CALL Info(Caller,'Operating on the lower part with: '//TRIM(Oper),Level=10)
    END IF
    

    ! Check that the variable exists for most of the operators 
    !----------------------------------------------------------
    IF( Oper == 'height' .OR. Oper == 'depth' .OR. Oper == 'index' .OR. &
        Oper == 'thickness' .OR. Oper == 'distance' ) THEN
      CONTINUE
    ELSE
      IF( .NOT. GotVar ) THEN
        CALL Fatal(Caller,'Variable required for this operator: '//TRIM(Oper))
      END IF
    END IF
        
    ! Create the projected variable if needed
    !-----------------------------------------------
    WRITE (Name,'(A,I0)') 'Target Variable ',NoVar
    TargetName = ListGetString( Params, TRIM(Name), GotIt )
    IF( .NOT. GotIt ) THEN
      WRITE (TargetName,'(A,A)') TRIM(Oper0)//' '//TRIM(VarName)
    END IF

    rdofs = dofs
    IF( Oper(1:5) == 'error' ) THEN
      Name = ListGetString( Params,'Error Variable '//I2S(NoVar),UnfoundFatal=.TRUE.)
      RefVar => VariableGet( Mesh % Variables, TRIM(Name) )
      IF(.NOT. ASSOCIATED(RefVar) ) THEN
        CALL Fatal(Caller,'Could not find required variable: '//TRIM(Name))
      END IF
      ! Error variable has always size one!
      rdofs = 1
    END IF
                
    IF( Oper == 'height' .OR. Oper == 'depth' .OR. Oper == 'index' .OR. Oper == 'distance') THEN
      ReducedDimensional = .FALSE.
    ELSE
      ReducedDimensional = .TRUE.
    END IF

    ProjectEverywhere = ListGetLogical( Params,'Project to everywhere',GotIt ) 
    IF(.NOT. GotIt) THEN
      WRITE (Name,'(A,I0,A)') 'Target Variable ',NoVar,' Everywhere' 
      ProjectEverywhere = ListGetLogical( Params, TRIM(Name), GotIt )
    END IF

    Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
    IF ( .NOT. ASSOCIATED( Var ) )  THEN      
      IF( ReducedDimensional .AND. .NOT. ProjectEverywhere ) THEN
        WRITE (Name,'(A,I0,A)') 'Target Variable ',NoVar,' At Bottom'
        IF( ListGetLogical( Params, TRIM(Name), GotIt ) ) THEN
          PermOut => BotPerm
          CALL Info(Caller,'Creating variable '//&
              TRIM(Name)//' at bottom',Level=8)
          GotIt = .TRUE.
        END IF
        IF( .NOT. GotIt .AND. MidLayerExists ) THEN
          WRITE (Name,'(A,I0,A)') 'Target Variable ',NoVar,' At Middle'
          IF( ListGetLogical( Params, TRIM(Name), GotIt ) ) THEN
            PermOut => MidPerm
            CALL Info(Caller,'Creating variable '//&
                TRIM(Name)//' at middle',Level=8)
            GotIt = .TRUE.
          END IF
        END IF
        IF(.NOT. GotIt ) THEN
          PermOut => TopPerm
          CALL Info(Caller,'Creating variable '//&
              TRIM(Name)//' at top',Level=8)
          GotIt = .TRUE.
        END IF
      ELSE
        IF( MaskExist ) THEN
          PermOut => MaskPerm
        ELSE        
          IF(.NOT. ASSOCIATED( UnitPerm ) ) THEN
            ALLOCATE( UnitPerm( nsize ) ) 
            DO i=1,nsize
              UnitPerm(i) = i
            END DO
          END IF
          PermOut => UnitPerm 
        END IF
      END IF

      CALL VariableAddVector( Mesh % Variables, Solver % Mesh, PSolver, &
          TargetName, rDofs, Perm = PermOut)           
      Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info(Caller,'Created variable: '//TRIM(TargetName),Level=9)
      ELSE
        CALL Fatal(Caller,'Could not create variable: '//TRIM(TargetName))
      END IF 
    END IF
    IF( Var % Dofs /= rDofs ) THEN
      CALL Fatal(Caller,'Mismatch in the dofs in fields!')
    END IF

    FieldOut => Var % Values
    PermOut => Var % Perm    
    FieldOut = 0.0_dp 

    IF(Oper == 'isosurface') THEN
      WRITE (Name,'(A,I0)') 'Isosurface Variable ',NoVar
      LevelsetName = ListGetString(Params,TRIM(Name),GotIt )
      IF(GotIt) THEN
        Var => VariableGet( Model % Variables, TRIM(LevelsetName) )
        Levelset => Var % Values
        LevelsetPerm => Var % Perm
      ELSE       
        Levelset => Coord
        NULLIFY(LevelsetPerm)
      END IF
      
      WRITE (Name,'(A,I0)') 'Isosurface Value ',NoVar
      Level = ListGetConstReal(Params,TRIM(Name),GotIt)
    ELSE IF ( SEQL(Oper, 'layer') ) THEN
      WRITE (Name,'(A,I0)') 'Layer Index ',NoVar
      layer = GetInteger( Params, Name, GotIt )
      IF (.NOT.GotIt) THEN
        CALL FATAL(Caller,'no "'//TRIM(Name)//'" indicated for output')
      END IF
    END IF

    ! Loop over components
    !------------------------------------------------
    DO dof = 1, rDofs

      ! Operators for dimensional reduction
      !-----------------------------------------------
      SELECT CASE(Oper)      
        
      CASE ('sum')      
        TopField = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF
                    
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          itop = TopPointer(j)
          k = i
          IF(ASSOCIATED(PermIn)) k = PermIn(k)
          IF(k == 0) CYCLE
          
          k = Dofs*(k-1)+dof
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + FieldIn(k)
        END DO
        
      CASE ('min')      
        TopField = HUGE(TopField)
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF
         
          itop = TopPointer(j)
          k = i
          IF(ASSOCIATED(PermIn)) k = PermIn(k)
            
          IF(k == 0) CYCLE
          k = Dofs*(k-1)+dof
          TopField(TopPerm(itop)) = MIN( FieldIn(k),TopField(TopPerm(itop)))
        END DO
        
      CASE ('max')      
        TopField = -HUGE(TopField)
        DO i=1,nnodes
          
          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          itop = TopPointer(j)
          k = i
          IF(ASSOCIATED(PermIn)) k = PermIn(k)
            
          IF(k == 0) CYCLE
          k = Dofs*(k-1)+dof
          TopField(TopPerm(itop)) = MAX( FieldIn(k),TopField(TopPerm(itop)))
        END DO
        
      CASE ('bottom')
        TopField = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF
         
          IF( i == TmpBotPointer(j) ) THEN
            k = i
            IF(ASSOCIATED(PermIn)) k = PermIn(k)
            k = Dofs*(k-1)+dof
            TopField(TopPerm(TopPointer(j))) = FieldIn(k)
          END IF
        END DO
        
      CASE ('top')
        TopField = 0.0_dp
        DO i=1,nnodes
          
          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
          
          IF( i == TmpTopPointer(j) ) THEN
            k = i
            IF(ASSOCIATED(PermIn)) k = PermIn(k)

            k = Dofs*(k-1)+dof
            TopField(TopPerm(TopPointer(j))) = FieldIn(k)
          END IF
        END DO

      CASE ('middle')
        TopField = 0.0_dp
               
        DO i=1,nnodes
          
          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
             
          IF( i == MidPointer(j) ) THEN
            k = i
            IF(ASSOCIATED(PermIn)) k = PermIn(k)

            k = Dofs*(k-1)+dof
            TopField(TopPerm(TopPointer(j))) = FieldIn(k)
          END IF
        END DO
        
      CASE ('layer below top')
        TopField = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
             
          IF( i == TmpTopPointer(j) ) THEN
            l = i
            DO k=1,layer
              IF( MaskExist ) THEN
                l = DownPointer(MaskPerm(l))
              ELSE
                l = DownPointer(l)
              END IF
            END DO
            IF(ASSOCIATED(PermIn)) l = PermIn(l)

            l = Dofs*(l-1)+dof
            TopField(TopPerm(i)) = FieldIn(l)
          END IF
        END DO
        
      CASE ('layer above bottom')
        TopField = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
                
          IF( i == TmpBotPointer(j) ) THEN
            l = i
            DO k=1,layer
              IF( MaskExist ) THEN
                l = UpPointer(MaskPerm(l))           
              ELSE                
                l = UpPointer(l)
              END IF
            END DO
            IF(ASSOCIATED(PermIn)) l = PermIn(l)

            l = Dofs*(l-1)+dof
            TopField(TopPerm(TopPointer(j))) = FieldIn(l)
          END IF
        END DO
        
      CASE ('isosurface')  ! not treated for mask!
        TopField = 0.0_dp
        DO i=1,nsize

          iup = UpPointer(i)
          j = i
          jup = iup
          IF(ASSOCIATED(LevelsetPerm)) THEN
            j = LevelsetPerm(j)
            jup = LevelsetPerm(jup)
          END IF
          IF(j == 0 .OR. jup == 0) CYCLE
          
          IF( (Levelset(jup) - Level) * (Levelset(j) - Level) <= 0.0_dp ) THEN
            itop = TopPointer(i)           
            dx = ABS(Levelset(jup) - Levelset(j))           
            l = i
            lup = iup
            IF(ASSOCIATED(PermIn)) THEN
              l = PermIn(l)
              lup = PermIn(lup)
            END IF
            l = Dofs*(l-1) + dof
            lup = Dofs*(lup-1) + dof
            IF( ABS(dx) < EPSILON(dx)) CALL Fatal(Caller,'dx smaller than machine Epsilon!')
            q = ABS(Levelset(jup)-Level) / dx          
            TopField(TopPerm(itop)) = q * FieldIn(l) + (1-q) * FieldIn(lup) 
          END IF
        END DO
        
      CASE ('int','int mean')

        TopField = 0.0_dp
        DO i=1,nnodes                   
          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          IF( Oper == 'int mean' ) THEN
            height = ABS(Coord(TopPointer(j)) - Coord(BotPointer(j)))
          ELSE
            height = 1.0_dp
          END IF
          
          
          ! Note for top and bottom this will automatically reduce the distance to half
          !----------------------------------------------------------------------------
          IF( i == TmpTopPointer(j) ) THEN
            iup = i
          ELSE
            iup = UpPointer(j)
          END IF

          IF( i == TmpBotPointer(j) ) THEN
            idown = i
          ELSE 
            idown = DownPointer(j)
          END IF

          IF( MaskExist ) THEN
            dx = 0.5*(Coord(MaskPerm(iup)) - Coord(MaskPerm(idown)))           
          ELSE
            dx = 0.5*(Coord(iup) - Coord(idown))
          END IF
          IF( ABS(height) < EPSILON(height)) CALL Fatal(Caller,'height smaller than machine Epsilon!')          
          dx = ABS( dx ) / height          
          k = i
          IF(ASSOCIATED(PermIn)) k = PermIn(k) 
            
          k = Dofs*(k-1) + dof
          itop = TopPointer(j)
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx * FieldIn(k)
        END DO

        
      CASE ('thickness')
        TopField = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
                    
          IF( UpperOper ) THEN  ! problem still with mask
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(j) ) THEN
            iup = i
          ELSE
            iup = UpPointer(j)
          END IF

          IF( i == TmpBotPointer(j) ) THEN
            idown = i
          ELSE 
            idown = DownPointer(j)
          END IF

          IF( MaskExist ) THEN
            dx = 0.5*(Coord(MaskPerm(iup)) - Coord(MaskPerm(idown)))            
          ELSE
            dx = 0.5*(Coord(iup) - Coord(idown))
          END IF
          dx = ABS( dx )
          itop = TopPointer(j)
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx 
        END DO

      ! Following four operators may have full dimensional results
      !--------------------------------------------------------------              
      CASE ('index')
        FieldOut = 0.0_dp
        DO i=1,nnodes
        
          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
  
          IF( UpperOper ) THEN  ! problem still with mask
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(j) ) THEN
            l = i
            
            DO k=1,nsize
              ll = l
              IF( MaskExist ) ll = MaskPerm(l)
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = 1.0_dp * k
                END IF
              ELSE
                FieldOut(l) = 1.0_dp * k
              END IF
              IF( l == TmpBotPointer(ll)) EXIT
              l = DownPointer(ll)
            END DO
          END IF
        END DO
 
      CASE ('depth') 
        FieldOut = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              

          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(j) ) THEN
            l = i            
            depth = 0.0_dp
            DO k=1,nnodes
              ll = l
              IF( MaskExist ) THEN
                ll = MaskPerm(l)
                IF( ll == 0 ) EXIT
              END IF
              IF( k > 1 ) THEN
                kk = UpPointer(ll)
                IF( MaskExist ) kk = MaskPerm(kk)
                depth = depth + ABS(Coord(kk) - Coord(ll))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = depth
                END IF
              ELSE
                FieldOut(l) = depth
              END IF
              IF( l == TmpBotPointer(ll)) EXIT
              l = DownPointer(ll)            
            END DO
          END IF
        END DO

      CASE ('height')
        FieldOut = 0.0_dp
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              

          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF
          
          IF( i == TmpBotPointer(j) ) THEN
            l = i            
            height = 0.0_dp
            DO k=1,nnodes
              ll = l
              IF( MaskExist ) THEN
                ll = MaskPerm(l)
                IF( ll == 0 ) EXIT
              END IF                
              IF( k > 1 ) THEN
                kk = DownPointer(ll)
                IF( MaskExist ) kk = MaskPerm(kk)
                height = height + ABS(Coord(ll) - Coord(kk))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = height
                END IF
              ELSE
                FieldOut(l) = height
              END IF
              IF( l == TmpTopPointer(ll)) EXIT
              l = UpPointer(ll)            
            END DO
          END IF
        END DO

 
      CASE ('distance') 
        FieldOut = 0.0_dp

        ! First check the distance to top ('depth')
        DO i=1,nnodes

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF
  
          IF( i == TmpTopPointer(j) ) THEN
            l = i
            depth = 0.0_dp
            DO k=1,nsize
              ll = l
              IF( MaskExist ) THEN
                ll = MaskPerm(l)
                IF( ll == 0 ) EXIT
              END IF
                
              IF( k > 1 ) THEN
                kk = UpPointer(ll)
                IF( MaskExist ) kk = MaskPerm(kk)
                depth = depth + ABS(Coord(kk) - Coord(ll))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = depth
                END IF
              ELSE
                FieldOut(l) = depth
              END IF
              IF( l == TmpBotPointer(ll)) EXIT
              l = DownPointer(ll)            
            END DO
          END IF
        END DO

        ! then distance to top ('height')
        DO i=1,nsize

          j = i
          IF( MaskExist ) THEN
            j = MaskPerm(i)
            IF( j == 0 ) CYCLE
          END IF              
          
          IF( UpperOper ) THEN
            IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
          END IF

          IF( i == TmpBotPointer(j) ) THEN
            l = i
            height = 0.0_dp
            DO k=1,nsize
              ll = l
              IF( MaskExist ) THEN
                ll = MaskPerm(l)
                IF( ll == 0 ) EXIT
              END IF
              IF( k > 1 ) THEN
                kk = DownPointer(ll)
                IF( MaskExist ) kk = MaskPerm(kk)
                height = height + ABS(Coord(ll) - Coord(kk))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = MIN( height, FieldOut(PermOut(l)))
                END IF
              ELSE
                FieldOut(l) = MIN( height, FieldOut(l))
              END IF
              IF( l == TmpTopPointer(ll)) EXIT
              l = UpPointer(ll)            
            END DO
          END IF
        END DO

      CASE('error norm','error projected','error max') 

        BLOCK
          INTEGER :: dofs1, dofs2, dofs0, ii 
          REAL(KIND=dp) :: s2, c, nrm1, nrm2
          REAL(KIND=dp), ALLOCATABLE :: u1(:), u2(:)
                 
          dofs1 = RefVar % Dofs
          dofs2 = DofsIn
          dofs0 = MIN(dofs1,dofs2)
          ALLOCATE(u1(NoLayers*dofs0),u2(NoLayers*dofs0))
          TopField = 0.0_dp

          DO i=1,nnodes                   
            j = i
            IF( MaskExist ) THEN
              j = MaskPerm(i)
              IF( j == 0 ) CYCLE
            END IF
            
            IF( UpperOper ) THEN
              IF( Coord(j) < Coord(MidPointer(j) ) ) CYCLE
            ELSE IF( LowerOper ) THEN
              IF( Coord(j) > Coord(MidPointer(j) ) ) CYCLE
            END IF
                        
            IF( i == TmpBotPointer(j) ) THEN
              l = i            
              u1 = 0.0_dp
              u2 = 0.0_dp
              ii = 0
              ! The end of loop should be never reached...
              DO k=1,nnodes
                ll = l
                IF( MaskExist ) THEN
                  ll = MaskPerm(l)
                  IF( ll == 0 ) EXIT
                END IF

                DO kk=1,dofs0
                  ii = ii + 1
                  u1(ii) = RefVar % Values(dofs1*(RefVar % Perm(l)-1)+kk)
                  u2(ii) = FieldIn(dofs2*(PermIn(l)-1)+kk)
                END DO
                
                IF( l == TmpTopPointer(ll)) EXIT
                l = UpPointer(ll)            
              END DO

              c = 1.0_dp
              IF( Oper == 'error projected' ) THEN
                c = SUM(u1*u2)/SUM(u2**2)
                u2 = c*u2
              END IF

              IF( Oper == 'error max' ) THEN
                nrm1 = SUM(ABS(u1))/NoLayers
                nrm2 = SUM(ABS(u2))/NoLayers
                s2 = 0.5*MAXVAL(ABS(u1-u2))/(nrm1+nrm2)
              ELSE
                nrm1 = SQRT(SUM(u1**2))
                nrm2 = SQRT(SUM(u2**2))
                s2 = 0.5*SQRT((SUM(u1-u2)**2))/(nrm1+nrm2)
              END IF
              TopField(TopPerm(l)) = s2
            END IF
          END DO
        END BLOCK
        
      CASE default
        CALL Fatal(Caller,'Unknown operator: '//TRIM(Oper))
        
      END SELECT
      

      ! Finally copy the projected values to the target variable
      ! It could be at the top, but it could also be at the bottom, or everywhere
      !----------------------------------------------------------------------------
      IF( ReducedDimensional ) THEN
        CALL Info(Caller,'Copying from surface to whole body',Level=10)        
        k = 0
        DO i=1,nnodes
          j = i
          
          IF( ASSOCIATED( PermOut ) ) THEN
            j = PermOut(i)
            IF( j == 0 ) CYCLE
          END IF
          
          IF( MaskExist ) THEN
            IF( MaskPerm(i) == 0 ) CYCLE
            k = TopPerm( TopPointer( MaskPerm(i) ) )            
          ELSE
            k = TopPerm( TopPointer(i) )
          END IF
          j = rDofs*(j-1) + dof
          FieldOut(j) = TopField(k)
        END DO
      END IF

    END DO

    IF( NormInd == NoVar ) THEN
      Solver % Variable % Values = ComputeNorm(Solver, SIZE( FieldOut ), FieldOut ) 
    END IF

  END DO


  at1 = CPUTime()

  WRITE(Message,* ) 'Projection time: ',at1-at0
  CALL Info(Caller,Message)
  CALL Info(Caller,'------------------------------------------')

!------------------------------------------------------------------------------
END SUBROUTINE StructuredProjectToPlane
!------------------------------------------------------------------------------

!> \}
