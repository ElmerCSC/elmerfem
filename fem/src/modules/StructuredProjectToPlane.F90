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

  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  NormInd = ListGetInteger( Solver % Values,'Show Norm Index',GotIt)
  IF( NormInd > 0 ) THEN
    IF( .NOT. ListCheckPresent( Solver % Values,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable','-nooutput -global savescalars_var')
    END IF
  END IF

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
  INTEGER :: i,j,k,l,n,dim,Dofs,dof,itop,ibot,idown,iup,jup,lup,ii,jj,Rounds,nsize,layer, &
      ActiveDirection,elem,TopNodes,MidNodes,NoVar,BotNodes, NormInd
  INTEGER, POINTER :: MaskPerm(:),TopPointer(:),BotPointer(:),UpPointer(:),DownPointer(:),&
      MidPointer(:), NodeIndexes(:),TargetPointer(:),BotPerm(:),&
      PermOut(:),PermIn(:),LevelsetPerm(:),TopPerm(:),MidPerm(:),UnitPerm(:)=>NULL(), &
      TmpTopPointer(:), TmpBotPointer(:)
  LOGICAL :: GotIt, Found, Visited = .FALSE., Initialized = .FALSE.,&
      Debug, MaskExist, GotVar, GotOldVar, GotOper, BottomTarget, ReducedDimensional, &
      MidLayerExists, UpperOper, LowerOper
  REAL(KIND=dp) :: dx,UnitVector(3),ElemVector(3),DotPro,Eps,Length,Level,val,q,depth,height
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at0,at1,at2
#else
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
#endif
  REAL(KIND=dp), POINTER :: FieldOut(:), FieldIn(:), Levelset(:), Coord(:),TopField(:)
  TYPE(Variable_t), POINTER :: Var, OldVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Nodes_t) :: Nodes
  TYPE(ValueList_t),POINTER :: BC

  
  SAVE Visited,Nodes,Initialized,UnitVector,Coord,MaskExist,MaskPerm,TopPointer,&
      BotPointer,MidPointer, UpPointer,DownPointer,FieldOut,FieldIn,&
      TopNodes,MidNodes,TopPerm, MidPerm, TopField, BotNodes, BotPerm, nsize, &
      UnitPerm, MidLayerExists
 
  CALL Info( 'StructuredProjectToPlane','------------------------------------------',Level=4 )
  CALL Info( 'StructuredProjectToPlane','Performing projection on a structured mesh ',Level=4 )
  CALL Info( 'StructuredProjectToPlane','------------------------------------------',Level=4 )

!------------------------------------------------------------------------------
!   Initialize the pointers to top and bottom nodes 
!------------------------------------------------------------------------------

  Debug = .FALSE.
  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver


  IF( .NOT. Initialized ) THEN

    IF(Debug) CALL Info('StructuredProjectToPlane','start init')
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
    Initialized = .TRUE.

    TopNodes = 0
    ALLOCATE( TopPerm( Mesh % NumberOfNodes ) )
    TopPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(TopPointer(i) == i) THEN
        TopNodes = TopNodes + 1
        TopPerm(i) = TopNodes
      END IF
    END DO
    IF( TopNodes > 0 ) THEN
      ALLOCATE( TopField( TopNodes ) ) 
      TopField = 0.0_dp
    END IF
    CALL Info('StructuredProjectoToPlane','Number of top nodes: '//TRIM(I2S(TopNodes)),Level=10)

    BotNodes = 0
    ALLOCATE( BotPerm( Mesh % NumberOfNodes ) )
    BotPerm = 0
    DO i=1,Mesh % NumberOfNodes
      IF(BotPointer(i) == i) THEN
        BotNodes = BotNodes + 1
        BotPerm(i) = BotNodes
      END IF
    END DO
    CALL Info('StructuredProjectoToPlane','Number of bot nodes: '//TRIM(I2S(BotNodes)),Level=10)

    IF( MidLayerExists ) THEN
      MidNodes = 0
      ALLOCATE( MidPerm( Mesh % NumberOfNodes ) ) 
      MidPerm = 0
      DO i=1,Mesh % NumberOfNodes
        IF(MidPointer(i) == i) THEN
          MidNodes = MidNodes + 1
          MidPerm(i) = MidNodes
        END IF
      END DO
      CALL Info('StructuredProjectoToPlane','Number of mid nodes: '//TRIM(I2S(MidNodes)),Level=10)
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

  NormInd = ListGetInteger( Solver % Values,'Show Norm Index',GotIt)

  debug = .FALSE.

  DO WHILE(.TRUE.)

    NoVar = NoVar + 1    
    IF(Debug) PRINT *,'NoVar',NoVar

    WRITE (Name,'(A,I0)') 'Variable ',NoVar
    VarName = ListGetString( Params, TRIM(Name), GotVar )
    NULLIFY(Var)    
    IF(GotVar) THEN
      Var => VariableGet( Model % Variables, TRIM(VarName) )
      IF ( .NOT. ASSOCIATED( Var ) )  THEN
        CALL Fatal('StructuredProjectToPlane','Variable does not exist: '//TRIM(VarName))
      END IF
      FieldIn => Var % Values
      PermIn => Var % Perm
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
        CALL Warn('StructuredProjectToPlane','Not even one field and operator to treat?')
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
        CALL Fatal('StructuredProjectoToPlane','Upper operator cannot exist without midlayer')
      END IF
      UpperOper = .TRUE.
      Oper = TRIM(Oper0(7:))
      TmpBotPointer => MidPointer
      CALL Info('StructuredProjectToPlane','Operating on the upper part with: '//TRIM(Oper),Level=10)
    ELSE IF( Oper0(1:6) == 'lower ') THEN
      IF( .NOT. MidLayerExists ) THEN
        CALL Fatal('StructuredProjectoToPlane','Lower operator cannot exist without midlayer')
      END IF
      LowerOper = .TRUE.
      Oper = TRIM(Oper0(7:))
      TmpTopPointer => MidPointer     
      CALL Info('StructuredProjectToPlane','Operating on the lower part with: '//TRIM(Oper),Level=10)
    END IF
    


    ! Check that the variable exists for most of the operators 
    !----------------------------------------------------------
    IF( Oper == 'height' .OR. Oper == 'depth' .OR. Oper == 'index' .OR. &
        Oper == 'thickness' .OR. Oper == 'distance' ) THEN
      CONTINUE
    ELSE
      IF( .NOT. GotVar ) THEN
        CALL Fatal('StructuredProjectToPlane','Variable required for this operator: '//TRIM(Oper))
      END IF
    END IF

    ! Create the projected variable if needed
    !-----------------------------------------------
    WRITE (Name,'(A,I0)') 'Target Variable ',NoVar
    TargetName = ListGetString( Params, TRIM(Name), GotIt )
    IF( .NOT. GotIt ) THEN
      WRITE (TargetName,'(A,A)') TRIM(Oper0)//' '//TRIM(VarName)
    END IF

    IF( Oper == 'height' .OR. Oper == 'depth' .OR. Oper == 'index' .OR. Oper == 'distance' ) THEN
      ReducedDimensional = .FALSE.
    ELSE
      ReducedDimensional = .TRUE.
    END IF


    Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
    IF ( .NOT. ASSOCIATED( Var ) )  THEN      
      IF( ReducedDimensional ) THEN
        WRITE (Name,'(A,I0,A)') 'Target Variable ',NoVar,' At Bottom'
        IF( ListGetLogical( Params, TRIM(Name), GotIt ) ) THEN
          PermOut => BotPerm
          CALL Info('StructuredProjectToPlane','Creating variable '//&
              TRIM(Name)//' at bottom',Level=8)
          GotIt = .TRUE.
        END IF
        IF( .NOT. GotIt .AND. MidLayerExists ) THEN
          WRITE (Name,'(A,I0,A)') 'Target Variable ',NoVar,' At Middle'
          IF( ListGetLogical( Params, TRIM(Name), GotIt ) ) THEN
            PermOut => MidPerm
            CALL Info('StructuredProjectToPlane','Creating variable '//&
                TRIM(Name)//' at middle',Level=8)
            GotIt = .TRUE.
          END IF
        END IF
        IF(.NOT. GotIt ) THEN
          PermOut => TopPerm
          CALL Info('StructuredProjectToPlane','Creating variable '//&
              TRIM(Name)//' at top',Level=8)
          GotIt = .TRUE.
        END IF
      ELSE        
        IF(.NOT. ASSOCIATED( UnitPerm ) ) THEN
          ALLOCATE( UnitPerm( nsize ) ) 
          DO i=1,nsize
            UnitPerm(i) = i
          END DO
        END IF
        PermOut => UnitPerm 
      END IF

      CALL VariableAddVector( Mesh % Variables, Solver % Mesh, PSolver, &
          TargetName, Dofs, Perm = PermOut)           
      Var => VariableGet( Mesh % Variables, TRIM(TargetName) )
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info('StructuredProjectToPlane','Created variable: '//TRIM(TargetName),Level=9)
      ELSE
        CALL Warn('StructuredProjectToPlane','Could not create variable: '//TRIM(TargetName))
      END IF 
    END IF
    IF( Var % Dofs /= Dofs ) THEN
      CALL Fatal('StructureProjectToPlane','Mismatch in the dofs in fields!')
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
        CALL FATAL('StructuredProjectToPlane','no > Layer Index < indicated for output')
      END IF
    END IF

    ! Loop over components
    !------------------------------------------------
    DO dof = 1, Dofs

      ! Operators for dimensional reduction
      !-----------------------------------------------
      SELECT CASE(Oper)      
        
      CASE ('sum')      
        TopField = 0.0_dp
        DO i=1,nsize
          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          itop = TopPointer(i)
          j = i
          IF(ASSOCIATED(PermIn)) j = PermIn(i)
          IF(j == 0) CYCLE

          j = Dofs*(j-1)+dof
          val = FieldIn(j)
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + val 
        END DO
        
      CASE ('min')      
        TopField = HUGE(TopField)
        DO i=1,nsize
          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF
         
          itop = TopPointer(i)
          j = i
          IF(ASSOCIATED(PermIn)) j = PermIn(i)

          IF(j == 0) CYCLE
          j = Dofs*(j-1)+dof
          TopField(TopPerm(itop)) = MIN( FieldIn(j),TopField(TopPerm(itop)))
        END DO
        
      CASE ('max')      
        TopField = -HUGE(TopField)
        DO i=1,nsize
          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          itop = TopPointer(i)
          j = i
          IF(ASSOCIATED(PermIn)) j = PermIn(i)

          IF(j == 0) CYCLE
          j = Dofs*(j-1)+dof
          TopField(TopPerm(itop)) = MAX( FieldIn(j),TopField(TopPerm(itop)))
        END DO
        
      CASE ('bottom')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF
         
          IF( i == TmpBotPointer(i) ) THEN
            j = i
            IF(ASSOCIATED(PermIn)) j = PermIn(i)
            j = Dofs*(j-1)+dof
            TopField(TopPerm(TopPointer(i))) = FieldIn(j)
          END IF
        END DO
        
      CASE ('top')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( i == TmpTopPointer(i) ) THEN
            j = i
            IF(ASSOCIATED(PermIn)) j = PermIn(i)
            j = Dofs*(j-1)+dof
            TopField(TopPerm(TopPointer(i))) = FieldIn(j)
          END IF
        END DO

      CASE ('middle')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( i == MidPointer(i) ) THEN
            j = i
            IF(ASSOCIATED(PermIn)) j = PermIn(i)
            j = Dofs*(j-1)+dof
            TopField(TopPerm(TopPointer(i))) = FieldIn(j)
          END IF
        END DO
        
      CASE ('layer below top')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( i == TmpTopPointer(i) ) THEN
            l = i
            DO k=1,layer
              l = DownPointer(l)
            END DO
            IF(ASSOCIATED(PermIn)) l = PermIn(l)
            l = Dofs*(l-1)+dof
            TopField(TopPerm(i)) = FieldIn(l)
          END IF
        END DO
        
      CASE ('layer above bottom')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( i == TmpBotPointer(i) ) THEN
            l = i
            DO k=1,layer
              l = UpPointer(l)
            END DO
            IF(ASSOCIATED(PermIn)) l = PermIn(l)
            l = Dofs*(l-1)+dof
            TopField(TopPerm(TopPointer(i))) = FieldIn(l)
          END IF
        END DO
        
      CASE ('isosurface')
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
            q = ABS(Levelset(jup)-Level) / dx          
            TopField(TopPerm(itop)) = q * FieldIn(l) + (1-q) * FieldIn(lup) 
          END IF
        END DO
        
      CASE ('int')
        TopField = 0.0_dp
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          ! Note for top and bottom this will automatically reduce the distance to half
          !----------------------------------------------------------------------------
          IF( i == TmpTopPointer(i) ) THEN
            iup = i
          ELSE
            iup = UpPointer(i)
          END IF

          IF( i == TmpBotPointer(i) ) THEN
            idown = i
          ELSE 
            idown = DownPointer(i)
          END IF

          dx = 0.5*(Coord(iup) - Coord(idown))
          j = i
          IF(ASSOCIATED(PermIn)) j = PermIn(i) 

          j = Dofs*(j-1) + dof
          itop = TopPointer(i)
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx * FieldIn(j)
        END DO
        
      CASE ('thickness')
        TopField = 0.0_dp
        DO i=1,nsize
          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(i) ) THEN
            iup = i
          ELSE
            iup = UpPointer(i)
          END IF

          IF( i == TmpBotPointer(i) ) THEN
            idown = i
          ELSE 
            idown = DownPointer(i)
          END IF

          dx = 0.5*(Coord(iup) - Coord(idown))
          itop = TopPointer(i)
          TopField(TopPerm(itop)) = TopField(TopPerm(itop)) + dx 
        END DO

      ! Following three operators may have full dimensional results
      !--------------------------------------------------------------      
  
      CASE ('index')
        FieldOut = 0.0_dp
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(i) ) THEN
            l = i
            DO k=1,nsize
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = 1.0_dp * k
                END IF
              ELSE
                FieldOut(l) = 1.0_dp * k
              END IF
              IF( l == TmpBotPointer(l)) EXIT
              l = DownPointer(l)            
            END DO
          END IF
        END DO
 
      CASE ('depth') 
        FieldOut = 0.0_dp
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          IF( i == TmpTopPointer(i) ) THEN
            l = i
            depth = 0.0_dp
            DO k=1,nsize
              IF( k > 1 ) THEN
                depth = depth + (Coord(UpPointer(l)) - Coord(l))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = depth
                END IF
              ELSE
                FieldOut(l) = depth
              END IF
              IF( l == TmpBotPointer(l)) EXIT
              l = DownPointer(l)            
            END DO
          END IF
        END DO

      CASE ('height')
        FieldOut = 0.0_dp
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          IF( i == TmpBotPointer(i) ) THEN
            l = i
            height = 0.0_dp
            DO k=1,nsize
              IF( k > 1 ) THEN
                height = height + (Coord(l) - Coord(DownPointer(l)))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = height
                END IF
              ELSE
                FieldOut(l) = height
              END IF
              IF( l == TmpTopPointer(l)) EXIT
              l = UpPointer(l)            
            END DO
          END IF
        END DO

 
      CASE ('distance') 
        FieldOut = 0.0_dp

        ! First check the distance to top ('depth')
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF
  
          IF( i == TmpTopPointer(i) ) THEN
            l = i
            depth = 0.0_dp
            DO k=1,nsize
              IF( k > 1 ) THEN
                depth = depth + (Coord(UpPointer(l)) - Coord(l))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = depth
                END IF
              ELSE
                FieldOut(l) = depth
              END IF
              IF( l == TmpBotPointer(l)) EXIT
              l = DownPointer(l)            
            END DO
          END IF
        END DO

        ! then distance to top ('height')
        DO i=1,nsize

          IF( UpperOper ) THEN
            IF( Coord(i) < Coord(MidPointer(i) ) ) CYCLE
          ELSE IF( LowerOper ) THEN
            IF( Coord(i) > Coord(MidPointer(i) ) ) CYCLE
          END IF

          IF( i == TmpBotPointer(i) ) THEN
            l = i
            height = 0.0_dp
            DO k=1,nsize
              IF( k > 1 ) THEN
                height = height + (Coord(l) - Coord(DownPointer(l)))
              END IF
              IF( ASSOCIATED(PermOut)) THEN
                IF( PermOut(l) > 0 ) THEN
                  FieldOut(PermOut(l)) = MIN( height, FieldOut(PermOut(l)))
                END IF
              ELSE
                FieldOut(l) = MIN( height, FieldOut(l))
              END IF
              IF( l == TmpTopPointer(l)) EXIT
              l = UpPointer(l)            
            END DO
          END IF
        END DO

        
      CASE default
        CALL Fatal('StructuredProjectToPlane','Unknown operator: '//TRIM(Oper))
        
      END SELECT
      

      ! Finally copy the projected values to the target variable
      ! It could be at the top, but it could also be at the bottom, or everywhere
      !----------------------------------------------------------------------------
      IF( ReducedDimensional ) THEN
        k = 0
        DO i=1,nsize
          j = i
          IF( ASSOCIATED( PermOut ) ) THEN
            j = PermOut(i)
            IF( j == 0 ) CYCLE
          END IF
          k = TopPerm( TopPointer(i) )
          j = Dofs*(j-1) + dof
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
  CALL Info('StructuredProjectToPlane',Message)
  CALL Info( 'StructuredProjectToPlane','------------------------------------------')

!------------------------------------------------------------------------------
END SUBROUTINE StructuredProjectToPlane
!------------------------------------------------------------------------------

!> \}
