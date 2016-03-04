!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh manipulation utilities for *Solver - routines
!------------------------------------------------------------------------------

MODULE MeshUtils

#ifdef USE_ISO_C_BINDINGS
    USE LoadMod
#ifdef HAVE_EIO
    USE EIOFortranAPI
#endif
#endif
    USE ElementUtils
    USE ElementDescription
    USE Interpolation
    USE ParallelUtils
    USE Types
    IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Allocated one single element. 
!------------------------------------------------------------------------------
   FUNCTION AllocateElement() RESULT( Element )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

     ALLOCATE( Element, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateElement', 'Unable to allocate a few bytes of memory?' )
     Element % BDOFs    =  0
     Element % NDOFs    =  0
     Element % BodyId   = -1
     Element % Splitted =  0
     Element % hK = 0
     Element % ElementIndex = 0
     Element % StabilizationMk = 0
     NULLIFY( Element % TYPE )
     NULLIFY( Element % PDefs )
     NULLIFY( Element % BubbleIndexes )
     NULLIFY( Element % DGIndexes )
     NULLIFY( Element % NodeIndexes )
     NULLIFY( Element % EdgeIndexes )
     NULLIFY( Element % FaceIndexes )
     NULLIFY( Element % BoundaryInfo )
!------------------------------------------------------------------------------
   END FUNCTION AllocateElement
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
   SUBROUTINE AllocatePDefinitions(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ALLOCATE(Element % PDefs, STAT=istat)
     IF ( istat /= 0) CALL Fatal('AllocatePDefinitions','Unable to allocate memory')

     ! Initialize fields
     Element % PDefs % P = 0 
     Element % PDefs % TetraType = 0
     Element % PDefs % isEdge = .FALSE.
     Element % PDefs % pyramidQuadEdge = .FALSE.
     Element % PDefs % localNumber = 0
     Element % PDefs % GaussPoints = 0

!------------------------------------------------------------------------------
   END SUBROUTINE AllocatePDefinitions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AllocateBoundaryInfo(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ALLOCATE(Element % BoundaryInfo, STAT=istat)
     IF ( istat /= 0) CALL Fatal('AllocateBoundaryInfo','Unable to allocate memory')

     Element % BoundaryInfo % Left => NULL()
     Element % BoundaryInfo % Right => NULL()
     Element % BoundaryInfo % GebhardtFactors => NULL()
     Element % BoundaryInfo % Constraint =  0

!------------------------------------------------------------------------------
   END SUBROUTINE AllocateBoundaryInfo
!------------------------------------------------------------------------------

!> Allocate mesh structure and return handle to it.
!------------------------------------------------------------------------------
   FUNCTION AllocateMesh() RESULT(Mesh)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     INTEGER :: istat

     ALLOCATE( Mesh, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateMesh', 'Unable to allocate a few bytes of memory?' )

!    Nothing computed on this mesh yet!
!    ----------------------------------
     Mesh % SavesDone    = 0
     Mesh % OutputActive = .FALSE.

     Mesh % AdaptiveDepth = 0
     Mesh % Changed   = .FALSE. !  TODO: Change this sometime

     Mesh % Stabilize = .FALSE.

     Mesh % Variables => NULL()
     Mesh % Parent => NULL()
     Mesh % Child => NULL()
     Mesh % Next => NULL()
     Mesh % RootQuadrant => NULL()
     Mesh % Elements => NULL()
     Mesh % Edges => NULL()
     Mesh % Faces => NULL()
     Mesh % Projector => NULL()
     Mesh % NumberOfEdges = 0
     Mesh % NumberOfFaces = 0
     Mesh % NumberOfNodes = 0
     Mesh % NumberOfBulkElements = 0
     Mesh % NumberOfBoundaryElements = 0
     Mesh % DiscontMesh = .FALSE.

     Mesh % InvPerm => NULL()

     Mesh % MaxFaceDOFs = 0
     Mesh % MaxEdgeDOFs = 0
     Mesh % MaxBDOFs = 0
     Mesh % MaxElementDOFs  = 0
     Mesh % MaxElementNodes = 0

     Mesh % ViewFactors => NULL()

     ALLOCATE( Mesh % Nodes, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateMesh', 'Unable to allocate a few bytes of memory?' )
     NULLIFY( Mesh % Nodes % x )
     NULLIFY( Mesh % Nodes % y )
     NULLIFY( Mesh % Nodes % z )
     Mesh % Nodes % NumberOfNodes = 0
     Mesh % NodesOrig => Mesh % Nodes
     NULLIFY( Mesh % NodesMapped )

     Mesh % EntityWeightsComputed = .FALSE.
     Mesh % BCWeight => NULL()
     Mesh % BodyForceWeight => NULL()
     Mesh % BodyWeight => NULL()
     Mesh % MaterialWeight => NULL()
    
     Mesh % ParallelInfo % NumberOfIfDOFs =  0
     NULLIFY( Mesh % ParallelInfo % GlobalDOFs )
     NULLIFY( Mesh % ParallelInfo % INTERFACE )
     NULLIFY( Mesh % ParallelInfo % NeighbourList )
!------------------------------------------------------------------------------
   END FUNCTION AllocateMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE GetMaxDefs(Model, Mesh, Element, ElementDef, SolverId, BodyId, Def_Dofs)
!------------------------------------------------------------------------------
     CHARACTER(*) :: ElementDef
     TYPE(Model_t) :: Model
     TYPE(MEsh_t) :: Mesh
     TYPE(Element_t) :: Element
     INTEGER :: SolverId, BodyId, Def_Dofs(:,:)

     TYPE(ValueList_t), POINTER :: Params
     INTEGER :: i, j,k,l, n, slen
     INTEGER, POINTER :: Body_Dofs(:,:)
     LOGICAL  :: stat, Found
     REAL(KIND=dp) :: x,y,z
     TYPE(Solver_t), POINTER  :: Solver
     CHARACTER(MAX_NAME_LEN) :: str, RESULT

     TYPE(ValueList_t), POINTER :: BodyParams
     CHARACTER(MAX_NAME_LEN) :: ElementDefBody
     
     BodyParams => Model % Bodies(BodyId) % Values

     ElementDefBody=ListGetString(BodyParams,'Solver '//TRIM(i2s(SolverId))//': Element',Found )
     IF (Found) THEN
       CALL Info('GetMaxDefs','Element found for body '//TRIM(i2s(BodyId))//' with solver '//TRIM(i2s(SolverId)), Level=5) 
       CALL Info('GetMaxDefs','Default element type is: '//ElementDef, Level=5)
       CALL Info('GetMaxDefs','New element type for this body is now: '//ElementDefBody, Level=5)
       ElementDef=ElementDefBody
     END IF

     Solver => Model % Solvers(SolverId)
     Params => Solver % Values

     IF ( .NOT. ALLOCATED(Solver % Def_Dofs) ) THEN
       ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
       Solver % Def_Dofs=-1
       Solver % Def_Dofs(:,:,1)=1
     END IF
     Body_Dofs => Solver % Def_Dofs(1:8,BodyId,:)

     j = INDEX(ElementDef, '-') ! FIX this to include elementtypewise defs...
     IF ( j>0 ) RETURN

     j = INDEX( ElementDef, 'n:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Body_Dofs(:,1) = l
       Def_Dofs(:,1) = MAX(Def_Dofs(:,1), l)
     END IF
          
      j = INDEX( ElementDef, 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,2) = l
        Def_Dofs(1:8,2) = MAX(Def_Dofs(1:8,2), l )
      END IF
          
      j = INDEX( ElementDef, 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,3) = l
        Def_Dofs(1:8,3) = MAX(Def_Dofs(1:8,3), l )
      END IF
          
      j = INDEX( ElementDef, 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,4) = l
        Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Body_Dofs(:,4) = 0
          Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef, 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(1:8,5) = l
        Def_Dofs(1:8,5) = MAX(Def_Dofs(1:8,5), l )
      END IF
          
      j = INDEX( ElementDef, 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2) == '%' ) THEN
          n = Element % TYPE % NumberOfNodes
          x = SUM(Mesh % Nodes % x(Element % NodeIndexes))/n
          y = SUM(Mesh % Nodes % y(Element % NodeIndexes))/n
          z = SUM(Mesh % Nodes % z(Element % NodeIndexes))/n
          WRITE( str, * ) 'cx= ',TRIM(i2s(Element % ElementIndex)),x,y,z
          str = TRIM(str) // '; ' // TRIM(ElementDef(j+3:))//'(cx)'
          slen = LEN_TRIM(str)
          CALL matc(str,RESULT,slen)
          READ(RESULT,*) x
          Body_Dofs(:,6) = 0
          Def_Dofs(1:8,6)  = MAX(Def_Dofs(1:8,6),NINT(x))
        ELSE
          READ( ElementDef(j+2:), * ) l
          Body_Dofs(:,6) = l
          Def_Dofs(1:8,6) = MAX(Def_Dofs(1:8,6), l )
        END IF
      END IF

!------------------------------------------------------------------------------
END SUBROUTINE GetMaxDefs
!------------------------------------------------------------------------------


  SUBROUTINE MarkHaloNodes( Mesh, HaloNode, FoundHaloNodes )

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, POINTER :: HaloNode(:)
    LOGICAL :: FoundHaloNodes

    INTEGER :: n,t
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: AllocDone

    ! Check whether we need to skip some elements and nodes on the halo boundary 
    ! We don't want to create additional nodes on the nodes that are on the halo only 
    ! since they just would create further need for new halo...
    FoundHaloNodes = .FALSE.
    IF( ParEnv % PEs > 1 ) THEN
      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF( ParEnv % MyPe /= Element % PartIndex ) THEN
          FoundHaloNodes = .TRUE.
          EXIT
        END IF
      END DO
    END IF


    ! If we have halo check the truly active nodes
    IF( FoundHaloNodes ) THEN
      CALL Info('MarkHaloNodes',&
          'Checking for nodes that are not really needed in bulk assembly',Level=12)

      IF( .NOT. ASSOCIATED( HaloNode ) ) THEN
        ALLOCATE( HaloNode( Mesh % NumberOfNodes ) )
        AllocDone = .TRUE.
      ELSE
        AllocDone = .FALSE.
      END IF

      ! Node is a halo node if it is not needed by any proper element
      HaloNode = .TRUE.
      DO t = 1, Mesh % NumberOfBulkElements     
        Element => Mesh % Elements(t)
        IF( ParEnv % MyPe == Element % PartIndex ) THEN
          Indexes => Element % NodeIndexes
          HaloNode( Indexes ) = .FALSE.
        END IF
      END DO

      n = COUNT( HaloNode ) 
      FoundHaloNodes = ( n > 0 ) 
      CALL Info('MarkHaloNodes','Number of passive nodes in the halo: '&
          //TRIM(I2S(n)),Level=10)

      ! If there are no halo nodes and the allocation was done within this subroutine
      ! then deallocate also. 
      IF( .NOT. FoundHaloNodes .AND. AllocDone ) THEN
        DEALLOCATE( HaloNode ) 
      END IF
    END IF

  END SUBROUTINE MarkHaloNodes

 

!> Create a discontinuous mesh over requested boundaries.
!> The nodes are duplicated in order to facilitate the discontinuity.
!> The duplicate nodes are not created by default if the connectivity 
!> of the nodes is needed by other bulk elements than those directly 
!> associated with the discontinuous boundaries. 
!------------------------------------------------------------------------------
 SUBROUTINE CreateDiscontMesh( Model, Mesh, DoAlways )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL, OPTIONAL :: DoAlways

   INTEGER, POINTER :: DisContPerm(:)
   LOGICAL, ALLOCATABLE :: DisContNode(:), DisContElem(:), ParentUsed(:), &
       MovingNode(:), StayingNode(:)
   LOGICAL :: Found, DisCont, GreedyBulk, GreedyBC, Debug, DoubleBC, UseTargetBodies, &
       UseConsistantBody, LeftHit, RightHit, Moving, Moving2, Set, Parallel
   INTEGER :: i,j,k,l,n,m,t,bc
   INTEGER :: NoNodes, NoDisContElems, NoDisContNodes, &
       NoBulkElems, NoBoundElems, NoParentElems, NoMissingElems, &
       DisContTarget, NoMoving, NoStaying, NoStayingElems, NoMovingElems, &
       NoUndecided, PrevUndecided, NoEdges, Iter, ElemFamily, DecideLimit, &
       ActiveBCs, CandA, CandB, RightBody, LeftBody, ConflictElems
   INTEGER, TARGET :: TargetBody(1)
   INTEGER, POINTER :: Indexes(:),ParentIndexes(:),TargetBodies(:)
   TYPE(Element_t), POINTER :: Element, LeftElem, RightElem, ParentElem, OtherElem
   CHARACTER(MAX_NAME_LEN) :: DiscontFlag
   LOGICAL :: CheckForHalo
   LOGICAL, POINTER :: HaloNode(:)
   TYPE(ValueList_t), POINTER :: BCList

   LOGICAL :: DoneThisAlready = .FALSE.

   IF(.NOT.PRESENT(DoAlways)) THEN
     IF (DoneThisAlready) RETURN
   ELSE 
     IF(.NOT.DoAlways) THEN
       IF (DoneThisAlready) RETURN
     END IF
   END IF
   DoneThisAlready = .TRUE.

   Discont = .FALSE.
   DoubleBC = .FALSE.
   ActiveBCs = 0
   DO bc = 1,Model % NumberOfBCs
     DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
     ! If the target boundary / periodic bc / mortar bc is zero
     ! it refers to itself. Otherwise the boundary will be doubled.
     IF( DisCont ) THEN
       i = ListGetInteger( Model % BCs(bc) % Values,'Discontinuous BC',Found )
       j = ListGetInteger( Model % BCs(bc) % Values,'Periodic BC',Found )
       k = ListGetInteger( Model % BCs(bc) % Values,'Mortar BC',Found )
       l = ListGetInteger( Model % BCs(bc) % Values,'Contact BC',Found )
       DoubleBC = ( i + j + k + l > 0 )
       ActiveBCs = ActiveBCs + 1
       BCList => Model % BCs(bc) % Values
     END IF
   END DO
   IF(ActiveBCs == 0 ) RETURN
   
   CALL Info('CreateDiscontMesh','Creating discontinuous boundaries')

   IF( ActiveBCs > 1 ) THEN
     CALL Warn('CreateDiscontMesh','Be careful when using more than one > Discontinuous Boundary < !')
   END IF

   Parallel = ( ParEnv % PEs > 1 )

   NoNodes = Mesh % NumberOfNodes
   NoBulkElems = Mesh % NumberOfBulkElements
   NoBoundElems = Mesh % NumberOfBoundaryElements
   
   ALLOCATE( DisContNode(NoNodes))
   ALLOCATE( DisContElem(NoBoundElems))
   ALLOCATE( ParentUsed(NoBulkElems))
   DisContNode = .FALSE.
   DisContElem = .FALSE.
   ParentUsed = .FALSE.
   NoDisContElems = 0
   NoMissingElems = 0


   ! Check whether we need to skip some elements and nodes on the halo boundary 
   ! We might not want to create additional nodes on the nodes that are on the halo only 
   ! since they just would create further need for new halo...
   CheckForHalo = ListGetLogical( Model % Simulation,'No Discontinuous Halo',Found ) 
   IF(.NOT. Found ) CheckForHalo = .TRUE.
   IF( CheckForHalo ) THEN
     HaloNode => NULL()
     CALL MarkHaloNodes( Mesh, HaloNode, CheckForHalo ) 
   END IF

   ! Go over all boundary elements and mark nodes that should be 
   ! discontinuous and nodes that should be continuous 
   DO t = 1, NoBoundElems
     
     Element => Mesh % Elements(NoBulkElems + t)
     Indexes => Element % NodeIndexes
     n = Element % Type % NumberOfNodes

     DisCont = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
         DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
         IF( DisCont ) EXIT
       END IF
     END DO     
     IF(.NOT. DisCont ) CYCLE
     
     DO i=1,n
       j = Indexes(i) 
       IF( CheckForHalo ) THEN
         IF( HaloNode(j) ) CYCLE
       END IF
       DisContNode(j) = .TRUE.
     END DO
     DisContElem( t ) = .TRUE.
     
     LeftElem => Element % BoundaryInfo % Left
     IF( ASSOCIATED( LeftElem ) ) THEN
       ParentUsed( LeftElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1 
     END IF
     
     RightElem => Element % BoundaryInfo % Right
     IF( ASSOCIATED( RightElem ) ) THEN
       ParentUsed( RightElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1
     END IF
   END DO
   
   IF( NoMissingElems > 0 ) THEN
     CALL Warn('CreateDiscontMesh','Missing '//TRIM(I2S(NoMissingElems))// &
     ' parent elements in partition '//TRIM(I2S(ParEnv % MyPe))) 
   END IF

   ! Calculate the number of discontinuous nodes and the number of bulk elements 
   ! associated to them. 
   NoDisContElems = COUNT( DiscontElem )
   NoDisContNodes = COUNT( DisContNode ) 
   CALL Info('CreateDiscontMesh','Number of discontinuous boundary elements: '&
       //TRIM(I2S(NoDisContElems)),Level=7)
   CALL Info('CreateDiscontMesh','Number of candicate nodes: '&
       //TRIM(I2S(NoDisContNodes)),Level=7)

   ! By default all nodes that are associated to elements immediately at the discontinuous 
   ! boundary are treated as discontinuous. However, the user may be not be greedy and release
   ! some nodes from the list that are associated also with other non-discontinuous elements.   
   ConflictElems = 0
   IF( NoDiscontNodes > 0 ) THEN
     n = NoDiscontNodes
     
     GreedyBulk = ListGetLogical( Model % Simulation,'Discontinuous Bulk Greedy',Found ) 
     IF(.NOT. Found ) GreedyBulk = .TRUE.     
     
     GreedyBC = ListGetLogical( Model % Simulation,'Discontinuous Boundary Greedy',Found ) 
     IF(.NOT. Found ) GreedyBC = .TRUE.     
     
     IF( .NOT. ( GreedyBC .AND. GreedyBulk ) ) THEN
       CALL Info('CreateDiscontMesh','Applying non-greedy strategies for Discontinuous mesh',Level=12)

       DO t = 1,NoBulkElems+NoBoundElems
         Element => Mesh % Elements(t)

         IF( t <= NoBulkElems ) THEN
           IF( GreedyBulk ) CYCLE
           IF( ParentUsed(t) ) CYCLE
         ELSE
           IF( GreedyBC ) CYCLE
           IF( DiscontElem(t-NoBulkElems) ) CYCLE
           !IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE
           ! Check that this is not an internal BC
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right) ) CYCLE
         END IF
         Indexes => Element % NodeIndexes

         IF( ANY( DisContNode( Indexes ) ) ) THEN
           !PRINT *,'t',Element % BoundaryInfo % Constraint, t,DisContElem(t), &
           !    Indexes, DisContNode( Indexes ) 
           DisContNode( Indexes ) = .FALSE.
           ConflictElems = ConflictElems + 1
         END IF
       END DO
       NoDisContNodes = COUNT( DisContNode ) 
     END IF

     IF( ConflictElems > 0 ) THEN
       CALL Info('CreateDiscontMesh','Conflicting discontinuity in elements: '&
           //TRIM(I2S(ConflictElems)))
     END IF

     IF( NoDiscontNodes < n ) THEN
       CALL Info('CreateDiscontMesh','Number of local discontinuous nodes: '&
           //TRIM(I2S(NoDisContNodes)), Level=12)
     ELSE
       CALL Info('CreateDiscontMesh','All candidate nodes used',Level=12)
     END IF
     
     IF( NoDiscontNodes == 0 ) THEN
       IF( n > 0 .AND. .NOT. GreedyBulk ) THEN
         CALL Info('CreateDiscontMesh','You might want to try the Greedy bulk strategy',Level=3)
       END IF
     END IF
   END IF
   
   i = NINT( ParallelReduction( 1.0_dp * NoDiscontNodes ) )
   CALL Info('CreateDiscontMesh','Number of discontinuous nodes: '&
       //TRIM(I2S(i)),Level=7)

   IF( i == 0 ) THEN
     CALL Warn('CreateDiscontMesh','Nothing to create, exiting...')
     IF( CheckForHalo ) DEALLOCATE( HaloNode ) 
     DEALLOCATE( DiscontNode, DiscontElem, ParentUsed )
     RETURN
   END IF

   ! Ok, we have marked discontinuous nodes, now give them an index. 
   ! This should also create the indexes in parallel.
   DisContPerm => NULL()
   ALLOCATE( DisContPerm(NoNodes) )
   DisContPerm = 0    

   ! We could end up here on an parallel case only
   ! Then we must make the parallel numbering, so jump to the end where this is done. 
   IF( NoDisContNodes == 0 ) THEN
     IF( DoubleBC ) THEN       
       Mesh % DiscontMesh = .FALSE.
       DEALLOCATE( DisContPerm ) 
     ELSE
       Mesh % DisContMesh = .TRUE.
       Mesh % DisContPerm => DisContPerm
       Mesh % DisContNodes = 0
     END IF
     GOTO 200
   END IF
   
   ! Create a table showing nodes that are related to the moving nodes by
   ! the moving elements. 
   ALLOCATE( MovingNode( NoNodes ), StayingNode( NoNodes ) ) 
   MovingNode = .FALSE.
   StayingNode = .FALSE.

   ! For historical reasons there is both single 'body' and multiple 'bodies'
   ! that define on which side of the discontinuity the new nodes will be. 
   DiscontFlag = 'Discontinuous Target Bodies'
   TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   IF(.NOT. UseTargetBodies ) THEN
     DiscontFlag = 'Discontinuous Target Body'
     TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   END IF

   ! If either parent is consistently one of the bodies then we can create a discontinuous 
   ! boundary. Note that this currently only works currently in serial!
   IF(.NOT. UseTargetBodies ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       CALL Fatal('CreateDiscontMesh','Please give > Discontinuous Target Bodies < on the BC!')
     END IF
     
     CALL Info('CreateDiscontMesh','Trying to find a dominating parent body',Level=12)

     CandA = -1
     CandB = -1
     DO t=1, NoBoundElems
       IF(.NOT. DisContElem(t) ) CYCLE
       Element => Mesh % Elements(NoBulkElems + t)

       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
         CALL Fatal('CreateDiscontMesh','Alternative strategy requires all parent elements!')
       END IF
       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
         CALL Fatal('CreateDiscontMesh','Alternative strategy requires all parent elements!')
       END IF

       LeftBody = Element % BoundaryInfo % Left % BodyId         
       RightBody = Element % BoundaryInfo % Right % BodyId

       IF( CandA == -1 ) THEN
         CandA = LeftBody 
       ELSE IF( CandA == 0 ) THEN
         CYCLE
       ELSE IF( CandA /= LeftBody .AND. CandA /= RightBody ) THEN
         CandA = 0
       END IF

       IF( CandB == -1 ) THEN
         CandB = RightBody
       ELSE IF( CandB == 0 ) THEN
         CYCLE
       ELSE IF( CandB /= LeftBody .AND. CandB /= RightBody ) THEN
         CandB = 0
       END IF
     END DO

     ! Choose the bigger one to honor the old convention
     ! This eliminates at the same time the unsuccesfull case of zero. 
     TargetBody(1) = MAX( CandA, CandB ) 

     IF( TargetBody(1) > 0 ) THEN
       CALL Info('CreateDiscontMesh',&
           'There seems to be a consistant discontinuous body: '&
           //TRIM(I2S(TargetBody(1))),Level=8)
       UseConsistantBody = .TRUE.
       TargetBodies => TargetBody
     ELSE
       CALL Fatal('CreateDiscontMesh',&
           'No simple rules available for determining discontinuous body')
     END IF
   END IF


   ! Assume we have only one active BC and we know the list of discontinuous 
   ! target bodies there. Hence we have all the info needed to set the 
   ! discontinuous elements also for other bulk elements. 
   ! This could be made more generic...
   NoUndecided = 0
   NoMovingElems = 0 
   NoStayingElems = 0

   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE

     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       NoMovingElems = NoMovingElems + 1 
       MovingNode(Indexes) = .TRUE.
     ELSE
       StayingNode(Indexes) = .TRUE.
       NoStayingElems = NoStayingElems + 1
     END IF
   END DO

   CALL Info('CreateDiscontMesh','Number of bulk elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8)
   CALL Info('CreateDiscontMesh','Number of bulk elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8)

   ! Set discontinuous nodes only if there is a real moving node associted with it
   ! Otherwise we would create a zero to the permutation vector. 
   ! If there is just a staying node then no need to create discontinuity at this node.
   DiscontNode = DiscontNode .AND. MovingNode 

   ! Create permutation numbering for the discontinuous nodes   
   ! Doubling will be done only for nodes that have both parents
   j = 0
   DO i=1,NoNodes
     IF( DisContNode(i) ) THEN
       j = j + 1
       DisContPerm(i) = j
     END IF
   END DO
   IF( j < NoDiscontNodes ) THEN
     PRINT *,'Some discontinuous nodes only needed on the other side:',&
         ParEnv % MyPe, NoDiscontNodes-j
     NoDiscontNodes = j 
   END IF


   ! Now set the new indexes for bulk elements
   ! In parallel skip the halo elements
   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE
     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       DO i=1, SIZE(Indexes) 
         j = DisContPerm(Indexes(i))
         IF( j > 0 ) Indexes(i) = NoNodes + j
       END DO
     END IF
   END DO

    
   ! Now set also the unset boundary elements by following the ownership of the parent elements
   ! or the majority opinion if this is conflicting.
   DO t=1, NoBoundElems

     Element => Mesh % Elements(NoBulkElems + t)

     ! If the element has no constraint then there is no need to treat it
     IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE

     IF( DisContElem(t) ) THEN
       LeftElem => Element % BoundaryInfo % Left
       RightElem => Element % BoundaryInfo % Right

       IF( ASSOCIATED( LeftElem ) ) THEN
         Moving = ANY( TargetBodies == LeftElem % BodyId ) 
       ELSE
         Moving = .NOT. ANY( TargetBodies == RightElem % BodyId )
       END IF
       IF( Moving ) THEN
         Element % BoundaryInfo % Left => RightElem
         Element % BoundaryInfo % Right => LeftElem 
       END IF
       CYCLE
     END IF


     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE

     ElemFamily = Element % TYPE % ElementCode / 100 
     LeftElem => Element % BoundaryInfo % Left
     RightElem => Element % BoundaryInfo % Right

     ! The boundary element follows the parent element if it is clear what to do
     Set = .TRUE.
     IF( ASSOCIATED( LeftElem ) .AND. ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( TargetBodies == LeftElem % BodyId )
       Moving2 = ANY( TargetBodies == RightElem % BodyId ) 
       IF( Moving .NEQV. Moving2) THEN
         CALL Warn('CreateDiscontMesh','Conflicting moving information')
         !PRINT *,'Moving:',t,Element % BoundaryInfo % Constraint, &
         !    Moving,Moving2,LeftElem % BodyId, RightElem % BodyId
         Set = .FALSE.
       ELSE
         IF( Moving ) THEN
           Element % BoundaryInfo % Left => RightElem
           Element % BoundaryInfo % Right => LeftElem 
         END IF
       END IF
     ELSE IF( ASSOCIATED( LeftElem ) ) THEN
       Moving = ANY( LeftElem % NodeIndexes > NoNodes ) 
     ELSE IF( ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( RightElem % NodeIndexes > NoNodes )
     ELSE
       CALL Fatal('CreateDiscontMesh','Boundary BC has no parants!')
     END IF

     ! Otherwise we follow the majority rule
     IF( .NOT. Set ) THEN
       NoMoving = COUNT( MovingNode(Indexes) ) 
       NoStaying = COUNT( StayingNode(Indexes) ) 

       IF( NoStaying /= NoMoving ) THEN
         Moving = ( NoMoving > NoStaying )
         Set = .TRUE.
       END IF
     END IF

     ! Ok, finally set whether boundary element is moving or staying
     IF( Set ) THEN
       IF( Moving ) THEN
         NoMovingElems = NoMovingElems + 1 
         DO i=1, SIZE(Indexes) 
           j = DisContPerm(Indexes(i))
           IF( j > 0 ) Indexes(i) = NoNodes + j
         END DO
       ELSE
         NoStayingElems = NoStayingElems + 1
       END IF
     ELSE
       NoUndecided = NoUndecided + 1
     END IF
   END DO

   CALL Info('CreateDiscontMesh','Number of related elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8 )
   CALL Info('CreateDiscontMesh','Number of related elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8 )
   IF( NoUndecided == 0 ) THEN
     CALL Info('CreateDiscontMesh','All elements marked either moving or staying')
   ELSE
     CALL Info('CreateDiscontMesh','Number of related undecided elements: '//TRIM(I2S(NoUndecided)) )
     CALL Warn('CreateDiscontMesh','Could not decide what to do with some boundary elements!')
   END IF


   m = COUNT( DiscontNode .AND. .NOT. MovingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not moving: ',ParEnv % MyPe, m
   END IF

   m = COUNT( DiscontNode .AND. .NOT. StayingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not staying: ',ParEnv % MyPe, m
     DO i=1,SIZE(DisContNode)
       IF( DiscontNode(i) .AND. .NOT. StayingNode(i) ) THEN
         IF( ParEnv % PEs == 1 ) THEN
           PRINT *,'Node:',ParEnv % MyPe,i
         ELSE
           PRINT *,'Node:',ParEnv % MyPe,i,Mesh % ParallelInfo % GlobalDofs(i), &
               Mesh % ParallelInfo % NeighbourList(i) % Neighbours
         END IF
         PRINT *,'Coord:',ParEnv % MyPe, Mesh % Nodes % x(i), Mesh % Nodes % y(i)
       END IF
     END DO
   END IF

   !DEALLOCATE( MovingNode, StayingNode )

   ! Now add the new nodes also to the nodes structure
   ! and give the new nodes the same coordinates as the ones
   ! that they were derived from. 
   Mesh % NumberOfNodes = NoNodes + NoDisContNodes   
   CALL EnlargeCoordinates( Mesh ) 

   CALL Info('CreateDiscontMesh','Setting new coordinate positions',Level=12)
   DO i=1, NoNodes
     j = DisContPerm(i)
     IF( j > 0 ) THEN
       k = NoNodes + j
       Mesh % Nodes % x(k) = Mesh % Nodes % x(i)
       Mesh % Nodes % y(k) = Mesh % Nodes % y(i)
       Mesh % Nodes % z(k) = Mesh % Nodes % z(i)
     END IF
   END DO


   ! If the discontinuous boundary is duplicated then no information of it 
   ! is saved. The periodic and mortar conditions now need to perform
   ! searches. On the other hand the meshes may now freely move.,
   IF( DoubleBC ) THEN
     CALL Info('CreateDiscontMesh','Creating secondary boundary for Discontinuous gap',Level=10)

     CALL EnlargeBoundaryElements( Mesh, NoDiscontElems ) 

     NoDisContElems = 0
     DO t=1, NoBoundElems

       ! Is this a boundary to be doubled?
       IF(.NOT. DisContElem(t) ) CYCLE

       Element => Mesh % Elements(NoBulkElems + t)
       IF(.NOT. ASSOCIATED(Element) ) THEN
         CALL Fatal('CreateDiscontMesh','Element '//TRIM(I2S(NoBulkElems+t))//' not associated!')
       END IF
       Indexes => Element % NodeIndexes

       DisContTarget = 0
       Found = .FALSE.
       DO bc = 1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Discontinuous BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Mortar BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Periodic BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Contact BC',Found )
           IF( Found ) EXIT
         END IF
       END DO
       IF( .NOT. Found .OR. DisContTarget == 0 ) THEN
         CALL Fatal('CreateDiscontMesh','Nonzero target boundary must be given for all, if any bc!')
       END IF

       RightElem => Element % BoundaryInfo % Right
       LeftElem => Element % BoundaryInfo % Left 

       NoDisContElems = NoDisContElems + 1              
       j = NoBulkElems + NoBoundElems + NoDisContElems 

       OtherElem => Mesh % Elements( j )
       IF(.NOT. ASSOCIATED(OtherElem) ) THEN
         CALL Fatal('CreateDiscontMesh','Other elem '//TRIM(I2S(j))//' not associated!')
       END IF

       OtherElem = Element 
       OtherElem % TYPE => Element % TYPE

       NULLIFY( OtherElem % BoundaryInfo ) 
       ALLOCATE( OtherElem % BoundaryInfo ) 
       OtherElem % BoundaryInfo % Left => Element % BoundaryInfo % Right

       ! Now both boundary elements are just one sided. Remove the associated to the other side. 
       NULLIFY( Element % BoundaryInfo % Right ) 
       NULLIFY( OtherElem % BoundaryInfo % Right )

       NULLIFY( OtherElem % NodeIndexes )
       n = SIZE( Element % NodeIndexes ) 
       ALLOCATE( OtherElem % NodeIndexes( n ) ) 

       ! Ok, we found the element to manipulate the indexes. 
       ! The new index is numbered on top of the old indexes. 
       DO i=1,n
         j = Element % NodeIndexes(i) 
         IF( DisContPerm(j) > 0 ) THEN
           OtherElem % NodeIndexes(i) = NoNodes + DisContPerm(j)
         ELSE 
           OtherElem % NodeIndexes(i) = j
         END IF
       END DO

       OtherElem % BoundaryInfo % Constraint = DisContTarget
     END DO

     CALL Info('CreateDiscontMesh','Number of original bulk elements: '&
         //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=10)
     CALL Info('CreateDiscontMesh','Number of original boundary elements: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=10)
     CALL Info('CreateDiscontMesh','Number of additional boundary elements: '&
         //TRIM(I2S(NoDisContElems)),Level=10)

     Mesh % DiscontMesh = .FALSE.
   ELSE
     Mesh % DisContMesh = .TRUE.
     Mesh % DisContPerm => DisContPerm
     Mesh % DisContNodes = NoDisContNodes 
   END IF

200 CONTINUE


   CALL EnlargeParallelInfo(Mesh, DiscontPerm )
   IF( ParEnv % PEs > 1 ) THEN
     m = COUNT( Mesh % ParallelInfo % GlobalDofs == 0) 
     IF( m > 0 ) CALL Warn('CreateDiscontMesh','There are nodes with zero global dof index: '//TRIM(I2S(m)))
   END IF

   IF( DoubleBC .AND. NoDiscontNodes > 0 ) DEALLOCATE( DisContPerm )


   DEALLOCATE( DisContNode, DiscontElem )   
  
 END SUBROUTINE CreateDiscontMesh


!> Reallocate coordinate arrays for iso-parametric p-elements,
!> or if the size of nodes has been increased due to discontinuity. 
!> This does not seem to be necessary for other types of 
!> elements (face, edge, etc.)
! -----------------------------------------------------------    
 SUBROUTINE EnlargeCoordinates(Mesh)

   TYPE(Mesh_t) :: Mesh
   INTEGER :: n0, n
   REAL(KIND=dp), POINTER :: TmpCoord(:)

   INTEGER :: i
   LOGICAL :: pelementsPresent

   n = Mesh % NumberOfNodes + &
       Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
       Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
       Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements
   n0 = SIZE( Mesh % Nodes % x )

   pelementsPresent = .FALSE.
   DO i=1,Mesh % NumberOfBulkElements
     IF(isPelement(Mesh % Elements(i))) THEN
       pelementsPresent = .TRUE.; EXIT
     END IF
   END DO

   IF ( Mesh % NumberOfNodes > n0 .OR. n > n0 .AND. pelementsPresent ) THEN
     CALL Info('EnlargeCoordinates','Increasing number of nodes from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

     TmpCoord => Mesh % Nodes % x
     ALLOCATE( Mesh % Nodes % x(n) )
     Mesh % Nodes % x(1:n0) = TmpCoord
     Mesh % Nodes % x(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % y
     ALLOCATE( Mesh % Nodes % y(n) )
     Mesh % Nodes % y(1:n0) = TmpCoord
     Mesh % Nodes % y(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % z
     ALLOCATE( Mesh % Nodes % z(n) )
     Mesh % Nodes % z(1:n0) = TmpCoord
     Mesh % Nodes % z(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )
   END IF

 END SUBROUTINE EnlargeCoordinates


 
 SUBROUTINE EnlargeBoundaryElements(Mesh, DoubleElements )

   TYPE(Mesh_t) :: Mesh
   INTEGER :: DoubleElements
   INTEGER :: n,n0,i,j
   REAL(KIND=dp), POINTER :: TmpCoord(:)
   TYPE(Element_t), POINTER :: NewElements(:),OldElements(:), Element

   IF( DoubleElements == 0 ) RETURN

   n0 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
   n = n0 + DoubleElements

   CALL Info('EnlargeBoundaryElements','Increasing number of elements from '&
       //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

   OldElements => Mesh % Elements
   CALL AllocateVector( Mesh % Elements, n, 'EnlargeBoundaryElements' )
   DO i=1,n0
     Mesh % Elements(i) = OldElements(i)
     IF(ASSOCIATED(OldElements(i) % BoundaryInfo)) THEN
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Left)) &
           Mesh % Elements(i) % BoundaryInfo % Left => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Left % ElementIndex)
       
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Right)) &
           Mesh % Elements(i) % BoundaryInfo % Right => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Right % ElementIndex)
     END IF
   END DO

   DO i=n0+1,n
     Element => Mesh % Elements(i)

     Element % DGDOFs = 0
     Element % BodyId = 0
     Element % TYPE => NULL()
     Element % BoundaryInfo => NULL()
     Element % PDefs => NULL()
     Element % DGIndexes => NULL()
     Element % EdgeIndexes => NULL()
     Element % FaceIndexes => NULL()
     Element % BubbleIndexes => NULL()
   END DO

   DEALLOCATE( OldElements ) 
   Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements + DoubleElements

 END SUBROUTINE EnlargeBoundaryElements


 SUBROUTINE EnlargeParallelInfo( Mesh, DiscontPerm )

   TYPE(Mesh_t) :: Mesh
   INTEGER, POINTER :: DiscontPerm(:)

   INTEGER :: nmax,n0,n1,i,j,istat, goffset
   INTEGER, POINTER :: TmpGlobalDofs(:) 
   INTEGER, ALLOCATABLE :: Perm(:)
   LOGICAL, POINTER :: Intf(:)
   TYPE(NeighbourList_t), POINTER :: Nlist(:)

   IF ( ParEnv % PEs <= 1 ) RETURN

   ! As index offset use the number of nodes in the whole mesh
   goffset = ParallelReduction( MAXVAL(Mesh % ParallelInfo % GlobalDofs)*1._dp,2 )

   n0 = SIZE( Mesh % ParallelInfo % GlobalDofs )
   n1 = Mesh % NumberOfNodes 
   IF( n0 >= n1 ) THEN
     CALL Info('EnlargeParallelInfo','No need to grow: '&
         //TRIM(I2S(n0))//' vs. '//TRIM(I2S(n1)),Level=10)
     RETURN
   END IF
   
   CALL Info('EnlargeParallelInfo','Increasing global numbering size from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n1)),Level=8)

   ! Create permutation table for the added nodes
   ALLOCATE(Perm(n1)); Perm  = 0
   DO i=1,n0
     IF ( DiscontPerm(i) > 0 ) THEN
       Perm(DiscontPerm(i)+n0) = i
     END IF
   END DO

   ! Create the enlarged set of global nodes indexes
   ALLOCATE( TmpGlobalDofs(n1), STAT=istat )
   IF (istat /= 0) CALL Fatal('LoadMesh', 'Unable to allocate TmpGlobalDofs array.')
   TmpGlobalDofs = 0
   DO i=1,n0
     TmpGlobalDofs(i) = Mesh % ParallelInfo % GlobalDofs(i)
   END DO
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0) THEN
       TmpGlobalDofs(i) = TmpGlobalDOfs(j) + goffset
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % GlobalDofs)
   Mesh % ParallelInfo % GlobalDOfs => TmpGlobalDofs

   ! Create the enlarged list of neighbours
   ALLOCATE(Nlist(n1))
   DO i=1,n0
     IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
       Nlist(i) % Neighbours => &
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
     ELSE 
       Nlist(i) % Neighbours => NULL()
     END IF
   END DO

   DO i=n0+1,n1
     j = Perm(i)
     IF ( j > 0 ) THEN
       IF( ASSOCIATED( Nlist(j) % Neighbours ) ) THEN
         ALLOCATE( Nlist(i) % Neighbours(SIZE(Nlist(j) % Neighbours) ) )
         Nlist(i) % Neighbours = Nlist(j) % Neighbours
       ELSE
         Nlist(i) % Neighbours => NULL()
       END IF
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
   Mesh % ParallelInfo % NeighbourList => Nlist


   ! Create logical table showing the interface nodes
   ALLOCATE( Intf(n1) )
   Intf = .FALSE.
   Intf(1:n0) = Mesh % ParallelInfo % INTERFACE(1:n0)
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0 ) THEN
       Intf(i) = Intf(j) 
     END IF
   END DO
   DEALLOCATE( Mesh % ParallelInfo % INTERFACE )
   Mesh % ParallelInfo % Interface => Intf


 END SUBROUTINE EnlargeParallelInfo




 !> Fortran reader for Elmer ascii mesh file format.
 !> This might be a Fortran replacement for the C++ eio library. 
 !------------------------------------------------------------------------
 SUBROUTINE ElmerAsciiMesh(Step, PMesh, MeshNamePar, ThisPe, IsParallel )

   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe
   LOGICAL, OPTIONAL :: IsParallel

   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: PrevStep=0, iostat
   INTEGER, PARAMETER :: FileUnit = 10
   CHARACTER(MAX_NAME_LEN) :: BaseName, FileName
   INTEGER :: i,j,k,n,BaseNameLen, SharedNodes = 0, mype
   INTEGER, POINTER :: NodeTags(:), ElementTags(:), LocalPerm(:)
   INTEGER :: MinNodeTag = 0, MaxNodeTag = 0, istat
   LOGICAL :: ElementPermutation=.FALSE., NodePermutation=.FALSE., Parallel



   SAVE PrevStep, BaseName, BaseNameLen, Mesh, mype, Parallel, &
       NodeTags, ElementTags, LocalPerm

   CALL Info('ElmerAsciiMesh','Performing step: '//TRIM(I2S(Step)),Level=8)

   IF( Step - PrevStep /= 1 ) THEN
     CALL Fatal('ElmerAsciiMesh','The routine should be called in sequence: '// &
         TRIM(I2S(PrevStep))//' : '//TRIM(I2S(Step)) )
   END IF
   PrevStep = Step
   IF( PrevStep == 6 ) PrevStep = 0 

   IF( Step == 1 ) THEN
     IF(.NOT. PRESENT( MeshNamePar ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give MeshNamePar!')
     END IF
     BaseName = TRIM( MeshNamePar ) 
     IF(.NOT. PRESENT( PMesh ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give PMesh!')
     END IF
     Mesh => PMesh
     IF(.NOT. PRESENT( ThisPe ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give ThisPe!')
     END IF
     mype = ThisPe 
     IF(.NOT. PRESENT( IsParallel ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give IsParallel!')
     END IF
     Parallel = IsParallel

     i = LEN_TRIM(MeshNamePar)
     DO WHILE(MeshNamePar(i:i) == CHAR(0))
       i=i-1
     END DO
     BaseNameLen = i
     CALL Info('LoadMesh','Base mesh name: '//TRIM(MeshNamePar(1:BaseNameLen)))
   END IF


   SELECT CASE( Step ) 

   CASE(1)       
     CALL ReadHeaderFile()

   CASE(2)
     CALL ReadNodesFile()

   CASE(3)
     CALL ReadElementsFile()

   CASE(4)
     CALL ReadBoundaryFile()
     CALL PermuteNodeNumbering()

   CASE(5)
     CALL InitParallelInfo()
     CALL ReadSharedFile()

   CASE(6) 
     IF( ASSOCIATED( LocalPerm) ) DEALLOCATE( LocalPerm ) 
     IF( ASSOCIATED( ElementTags) ) DEALLOCATE( ElementTags )

   END SELECT


 CONTAINS


   FUNCTION read_ints(s,j,halo) RESULT(n)
     INTEGER :: j(:)
     CHARACTER(LEN=*) :: s
     LOGICAL :: halo
     
     INTEGER :: i,k,l,m,n,ic
     INTEGER, PARAMETER :: ic0 = ICHAR('0'), ic9 = ICHAR('9'), icm = ICHAR('-'), &
         icd = ICHAR('/'), ics = ICHAR(' ')
     
     k = LEN_TRIM(s)
     l = 1
     n = 0
     halo = .FALSE.
     DO WHILE(l<=k.AND.n<SIZE(j))
       DO WHILE(l<=k)
         ic = ICHAR(s(l:l))
         IF( ic == ics ) THEN
           CONTINUE
         ELSE IF( ic == icd ) THEN
           halo = .TRUE.
         ELSE
           EXIT
         END IF
         l=l+1
       END DO
       IF(l>k) EXIT
       IF(.NOT.(ic==icm .OR. ic>=ic0 .AND. ic<=ic9)) EXIT
       
       m = l+1
       DO WHILE(m<=k)
         ic = ICHAR(s(m:m))
         IF(ic<ic0 .OR. ic>ic9) EXIT
         m=m+1
       END DO
       
       n = n + 1
       j(n) = s2i(s(l:m-1),m-l)
       l = m
     END DO
   END FUNCTION read_ints
   

   !---------------------------------------------------
   ! Read header file and allocate some mesh structures
   !---------------------------------------------------
   SUBROUTINE ReadHeaderFile()

     INTEGER :: TypeCount
     INTEGER :: Types(64),CountByType(64)

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(ParEnv % PEs))//&
           '/part.'//TRIM(I2S(mype+1))//'.header'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.header'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading header info from file: '//TRIM(FileName),Level=10)
     END IF

     READ(FileUnit,*,IOSTAT=iostat) Mesh % NumberOfNodes, &
         Mesh % NumberOfBulkElements,&
         Mesh % NumberOfBoundaryElements
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not read header 1st line in file: '//TRIM(FileName))
     END IF

     Types = 0
     CountByType = 0
     READ(FileUnit,*,IOSTAT=iostat) TypeCount
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not read the type count in file: '//TRIM(FileName))
     END IF
     DO i=1,TypeCount
       READ(FileUnit,*,IOSTAT=iostat) Types(i),CountByType(i)
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Could not read type count '&
             //TRIM(I2S(i))//'in file: '//TRIM(FileName))
       END IF
     END DO

     IF( Parallel ) THEN
       READ(FileUnit,*,IOSTAT=iostat) SharedNodes
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Could not read shared nodes in file: '//TRIM(FileName))
       END IF
     ELSE
       SharedNodes = 0
     END IF

     Mesh % MaxElementNodes = 0
     DO i=1,TypeCount
       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes, MODULO( Types(i), 100) )
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadHeaderFile


   !-----------------------------------------------------------------------
   ! Read nodes file and create nodal permutation if needed
   !-----------------------------------------------------------------------
   SUBROUTINE ReadNodesFile()

     REAL(KIND=dp) :: Coords(3)
     INTEGER :: NodeTag

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(ParEnv % PEs))//&
           '/part.'//TRIM(I2S(mype+1))//'.nodes'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.nodes'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ALLOCATE( NodeTags(Mesh % NumberOfNodes ) ) 
     NodeTags = 0

     NodePermutation = .FALSE.
     DO j = 1, Mesh % NumberOfNodes
       READ(FileUnit,*,IOSTAT=iostat) NodeTag, k, Coords
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Problem load node '//TRIM(I2S(j))//' in file: '//TRIM(Filename))
       END IF

       IF( NodeTags(j) /= j ) NodePermutation = .TRUE.
 
       NodeTags(j) = NodeTag
       Mesh % Nodes % x(j) = Coords(1)
       Mesh % Nodes % y(j) = Coords(2)
       Mesh % Nodes % z(j) = Coords(3)
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadNodesFile


   !------------------------------------------------------------------------------
   ! Read elements file and create elemental permutation if needed 
   !------------------------------------------------------------------------------
   SUBROUTINE ReadElementsFile()
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ElemType, Tag, Body, ElemNo, Ivals(64),nread, ioffset, partn
     CHARACTER(256) :: str
     LOGICAL :: halo


     CALL AllocateVector( ElementTags, Mesh % NumberOfBulkElements+1, 'LoadMesh')   
     ElementTags = 0
     ElementPermutation = .FALSE.

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)// &
          '/partitioning.'//TRIM(I2S(ParEnv % PEs))//&
             '/part.'//TRIM(I2S(mype+1))//'.elements'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.elements'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadElementsFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading bulk elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read start of element entry: '//TRIM(I2S(j)))
       END IF

       nread = read_ints(str,ivals,halo)

       tag = ivals(1)

       IF( halo ) THEN
         ioffset = 1
         partn = ivals(2) 
       ELSE
         ioffset = 0
         partn = 0 
       END IF
       body = ivals(ioffset+2)
       ElemType = ivals(ioffset+3)

       ElementTags(j) = tag
       IF( j /= tag ) ElementPermutation = .TRUE.             
       Element % ElementIndex = j
       Element % BodyId = body

       IF( partn > 0 ) THEN
         Element % PartIndex = partn-1
       ELSE
         Element % PartIndex = mype
       END IF

       Element % TYPE => GetElementType( ElemType )

       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadElementsFile','Element of type '&
             //TRIM(I2S(ElemType))//' could not be associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       IF( nread < n + ioffset + 3 ) THEN
         CALL Fatal('ReadElementsFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF

       CALL AllocateVector( Element % NodeIndexes, n )

       Element % NodeIndexes(1:n) = IVals(4+ioffset:nread)
     END DO
     CLOSE( FileUnit ) 

   END SUBROUTINE ReadElementsFile
   !------------------------------------------------------------------------------


   !------------------------------------------------------------------------------
   ! Read boundary elements file and remap the parents if needed.  
   !------------------------------------------------------------------------------
   SUBROUTINE ReadBoundaryFile()
     INTEGER, POINTER :: LocalEPerm(:)
     INTEGER :: MinEIndex, MaxEIndex, ElemNodes, i
     INTEGER :: Left, Right, bndry, tag, ElemType, IVals(64), nread, ioffset, partn
     TYPE(Element_t), POINTER :: Element
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(ParEnv % PEs))//&
           '/part.'//TRIM(I2S(mype+1))//'.boundary'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.boundary'
     END IF

     ! Create permutation for the elements. This is needed when the element 
     ! parents are mapped to the new order. This is needed for mapping of the 
     ! parents. Otherwise the element numbering is arbitrary. 
     !------------------------------------------------------------------------------
     IF( ElementPermutation ) THEN
       MinEIndex = MINVAL( ElementTags(1:Mesh % NumberOfBulkElements) )
       MaxEIndex = MAXVAL( ElementTags(1:Mesh % NumberOfBulkElements) )

       LocalEPerm => NULL()
       CALL AllocateVector( LocalEPerm, MaxEIndex - MinEIndex + 1, 'LoadMesh' )
       LocalEPerm = 0
       DO i=1,Mesh % NumberOfBulkElements
         LocalEPerm( ElementTags(i) - MinEIndex + 1 ) = i
       END DO
     ELSE
       MinEIndex = 1 
       MaxEIndex = Mesh % NumberOfBulkElements
     END IF


     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadBoundaryFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading boundary elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read boundary element entry: '//TRIM(I2S(j)))
       END IF
       nread = read_ints(str,ivals,halo)
       
       tag = ivals(1)

       IF( halo ) THEN
         partn = ivals(2)
         ioffset = 1
       ELSE
         partn = 0
         ioffset = 0
       END IF

       bndry = ivals(ioffset+2)
       left = ivals(ioffset+3)
       right = ivals(ioffset+4)
       ElemType = ivals(ioffset+5)
       
       Element % ElementIndex = j
       Element % TYPE => GetElementType( ElemType )
       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadBoundaryFile','Element of type '//TRIM(I2S(ElemType))//'could not be associated!')
       END IF

       ElemNodes = Element % TYPE % NumberOfNodes
       Mesh % MaxElementNodes = MAX( Mesh % MaxElementNodes, ElemNodes )

       IF( partn == 0 ) THEN
         Element % PartIndex = mype
       ELSE
         Element % PartIndex = partn-1
       END IF

       CALL AllocateBoundaryInfo( Element ) 

       Element % BoundaryInfo % Constraint = bndry
       Element % BoundaryInfo % Left => NULL()
       Element % BoundaryInfo % Right => NULL()

       IF ( Left >= MinEIndex .AND. Left <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Left  = LocalEPerm(Left - MinEIndex + 1)
         END IF
       ELSE IF ( Left > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag, Left
         CALL Error( 'ReadBoundaryFile', Message )
         Left = 0
       END IF

       IF ( Right >= MinEIndex .AND. Right <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Right = LocalEPerm(Right - MinEIndex + 1)
         END IF
       ELSE IF ( Right > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag,Right
         CALL Error( 'ReadBoundaryFile', Message )
         Right = 0
       END IF

       IF ( Left >= 1 ) THEN
         Element % BoundaryInfo % Left => Mesh % Elements(left)
       END IF

       IF ( Right >= 1 ) THEN
         Element % BoundaryInfo % Right => Mesh % Elements(right)
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( nread < 5 + n + ioffset ) THEN
         CALL Fatal('ReadBoundaryFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       Element % NodeIndexes(1:n) = Ivals(6+ioffset:nread)
     END DO
     CLOSE( FileUnit )


     IF( ElementPermutation ) THEN
       DEALLOCATE( LocalEPerm ) 
     END IF

   END SUBROUTINE ReadBoundaryFile
   !------------------------------------------------------------------------------



   ! Make a permutation for the bulk and boundary element topology if 
   ! the nodes are permuted. This is always the case in parallel.
   ! The initial numbering is needed only when the nodes are loaded and 
   ! hence this is a local subroutine. 
   !----------------------------------------------------------------------
   SUBROUTINE PermuteNodeNumbering()

     TYPE(Element_t), POINTER :: Element

     IF( NodePermutation ) THEN
       CALL Info('LoadMesh','Performing node mapping',Level=6)

       MinNodeTag = MINVAL( NodeTags )
       MaxNodeTag = MAXVAL( NodeTags )

       CALL AllocateVector( LocalPerm, MaxNodeTag-MinNodeTag+1, 'LoadMesh' )
       LocalPerm = 0
       DO i=1,Mesh % NumberOfNodes
         LocalPerm(NodeTags(i) - MinNodeTag + 1) = i
       END DO

       DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements       
         Element => Mesh % Elements(i)
         n = Element % TYPE % NumberOfNodes

         DO j=1,n
           k = Element % NodeIndexes(j) 
           Element % NodeIndexes(j) = LocalPerm(k - MinNodeTag + 1)
         END DO
       END DO
     ELSE
       CALL Info('LoadMesh','Node mapping is continuous',Level=8)
     END IF

     ! Set the for now, if the case is truly parallel we'll have to revisit these
     ! when reading the parallel information. 
     Mesh % ParallelInfo % NumberOfIfDOFs = 0
     Mesh % ParallelInfo % GlobalDOFs => NodeTags

   END SUBROUTINE PermuteNodeNumbering


   ! Initialize some parallel structures once the non-nodal 
   ! element types are known. 
   ! Currently this is here mainly because the 
   ! Elemental and Nodal tags are local
   !-------------------------------------------------------
   SUBROUTINE InitParallelInfo()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! These two have already been set, and if the case is serial
     ! case they can be as is.
     !Mesh % ParallelInfo % NumberOfIfDOFs = 0
     !Mesh % ParallelInfo % GlobalDOFs => NodeTags

     IF(.NOT. Parallel ) RETURN

     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i)
     END DO

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes)
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs

     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('LoadMesh', 'Unable to allocate NeighbourList array.')

     DO i=1,n
       NULLIFY( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % INTERFACE, n, 'LoadMesh')
     Mesh % ParallelInfo % INTERFACE = .FALSE.       

   END SUBROUTINE InitParallelInfo


   ! Read the file that shows the shared nodes.
   !------------------------------------------------------------------------
   SUBROUTINE ReadSharedFile()

     INTEGER :: Ivals(64)
     INTEGER :: npart, tag, nread
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF(.NOT. Parallel) RETURN

     FileName = BaseName(1:BaseNameLen)//&
       '/partitioning.'//TRIM(I2S(ParEnv % PEs))//&
         '/part.'//TRIM(I2S(mype+1))//'.shared'

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ! This loop could be made more effective, for example
     ! by reading tags and nparts to a temporal vector
     ! The operation using the str takes much more time.
     !-----------------------------------------------------
     DO i=1,SharedNodes          
       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read shared nodes entry: '//TRIM(I2S(i)))
       END IF
       nread = read_ints(str,ivals,halo)

       tag = ivals(1)
       npart = ivals(2)       

       k = LocalPerm( tag-MinNodeTag+1 )
       Mesh % ParallelInfo % INTERFACE(k) = .TRUE.
       CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(k) % Neighbours,npart)

       IF( nread < 2 + npart ) THEN
         CALL Fatal('ReadSharedFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       
       Mesh % ParallelInfo % NeighbourList(k) % Neighbours = ivals(3:nread) - 1

       ! this partition does not own the node
       IF ( ivals(3)-1 /= mype ) THEN
         Mesh % ParallelInfo % NumberOfIfDOFs = &
             Mesh % ParallelInfo % NumberOfIfDOFs + 1
       END IF
     END DO

     CLOSE( FileUnit )

   END SUBROUTINE ReadSharedFile

 END SUBROUTINE ElmerAsciiMesh



 !> An interface over potential mesh loading strateties. 
 !----------------------------------------------------------------- 
 SUBROUTINE LoadMeshStep( Step, PMesh, MeshNamePar, ThisPe, IsParallel ) 
   
   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe
   LOGICAL, OPTIONAL :: IsParallel

   ! Currently only one strategy to get the mesh is implemented 
   ! but there could be others.
   !
   ! This has not yet been tested in parallel and for sure
   ! it does not work for halo elements. 
   !-----------------------------------------------------------------
   CALL ElmerAsciiMesh( Step, PMesh, MeshNamePar, ThisPe, IsParallel ) 

 END SUBROUTINE LoadMeshStep



 !------------------------------------------------------------------------------
 !> Function to load mesh from disk.
 !------------------------------------------------------------------------------
 FUNCTION LoadMesh2( Model, MeshDirPar, MeshNamePar,&
     BoundariesOnly, NumProcs,MyPE, Def_Dofs ) RESULT( Mesh )
   !------------------------------------------------------------------------------
   USE PElementMaps, ONLY : GetRefPElementNodes

   IMPLICIT NONE

   CHARACTER(LEN=*) :: MeshDirPar,MeshNamePar
   LOGICAL :: BoundariesOnly    
   INTEGER, OPTIONAL :: numprocs,mype,Def_Dofs(:,:)
   TYPE(Mesh_t),  POINTER :: Mesh
   TYPE(Model_t) :: Model
   !------------------------------------------------------------------------------    
   INTEGER :: i,j,k,n
   INTEGER :: BaseNameLen, Save_Dim
   LOGICAL :: GotIt, Found
   CHARACTER(MAX_NAME_LEN) :: FileName
   TYPE(Element_t), POINTER :: Element
   TYPE(Matrix_t), POINTER :: Projector
   LOGICAL :: parallel, LoadNewMesh


   Mesh => Null()

   n = LEN_TRIM(MeshNamePar)
   DO WHILE (MeshNamePar(n:n)==CHAR(0).OR.MeshNamePar(n:n)==' ')
     n=n-1
   END DO
   IF(NumProcs<=1) THEN
     INQUIRE( FILE=MeshNamePar(1:n)//'/mesh.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Fatal('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < does not exist!')
     END IF
   ELSE
     INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// & 
         TRIM(i2s(Numprocs))//'/part.1.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Warn('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
           //TRIM(I2S(Numprocs))//' does not exist!')
       RETURN
     END IF
   END IF

   CALL Info('LoadMesh','Starting',Level=8)

   Parallel = .FALSE.
   IF ( PRESENT(numprocs) .AND. PRESENT(mype) ) THEN
     IF ( numprocs > 1 ) Parallel = .TRUE.
   END IF

   Mesh => AllocateMesh()

   ! Get sizes of mesh structures for allocation
   !--------------------------------------------------------------------
   CALL LoadMeshStep( 1, Mesh, MeshNamePar, mype, Parallel ) 

   ! Initilize and allocate mesh stuctures
   !---------------------------------------------------------------------
   CALL InitializeMesh()

   ! Get the (x,y,z) coordinates
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 2 )
   ! Permute and scale the coordinates.
   ! This also finds the mesh dimension. It is needed prior to getting the 
   ! elementtypes since wrong permutation or dimension may spoil that. 
   !-------------------------------------------------------------------
   CALL MapCoordinates()

   ! Get the bulk elements: element types, body index, topology
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 3 )

   ! Get the boundary elements: boundary types, boundary index, parents, topology
   !------------------------------------------------------------------------------
   CALL LoadMeshStep( 4 )

   ! Read elemental data - this is rarely used, parallel implementation lacking?
   !--------------------------------------------------------------------------
   i = LEN_TRIM(MeshNamePar)
   DO WHILE(MeshNamePar(i:i) == CHAR(0))
     i=i-1
   END DO
   BaseNameLen = i
   
   FileName = MeshNamePar(1:BaseNameLen)//'/mesh.elements.data'
   CALL ReadElementPropertyFile( FileName, Mesh )

   ! Read mesh.names - this could be saved by some mesh formats
   !--------------------------------------------------------------------------
   IF( ListGetLogical( Model % Simulation,'Use Mesh Names',Found ) ) THEN
     FileName = MeshNamePar(1:BaseNameLen)//'/mesh.names'
     CALL ReadTargetNames( Model, FileName )
   END IF


   ! Map bodies using Target Bodies and boundaries using Target Boundaries.
   ! This must be done before the element definitions are studied since
   ! then the pointer should be to the correct body index. 
   !------------------------------------------------------------------------
   CALL MapBodiesAndBCs()

   ! Read parallel mesh information: shared nodes
   !------------------------------------------------------------------
   CALL LoadMeshStep( 5 )

   ! Create the discontinuous mesh that accounts for the jumps in BCs
   ! This must be created after the whole mesh has been read in and 
   ! bodies and bcs have been mapped to full operation.
   ! To consider non-nodal elements it must be done before them.
   !--------------------------------------------------------------------
   CALL CreateDiscontMesh(Model,Mesh)

   ! Study the non-nodal elements (face, edge, DG, and p-elements)
   ! This must be done before parallel communication since it will 
   ! affect what needs to be communicated. 
   !-------------------------------------------------------------------
   CALL NonNodalElements()

   ! Create parallel info for the non-nodal elements
   !------------------------------------------------------------------
   CALL ParallelNonNodalElements()

   ! Deallocate some stuff no longer needed
   !------------------------------------------------------------------
   CALL LoadMeshStep( 6 )

   ! Enlarge the coordinate vectors.
   ! This must be done after the advanced elements have been detected.
   ! Currently increase is applied only for p-elements. 
   !-------------------------------------------------------------------
   CALL EnlargeCoordinates(Mesh)       

   ! Some physics related historical initializations
   !-----------------------------------------------------
   Model % FreeSurfaceNodes => NULL()
   Model % BoundaryCurvatures => NULL()

   ! If periodic BC given, compute boundary mesh projector:
   ! ------------------------------------------------------
   DO i = 1,Model % NumberOfBCs
     Model % BCs(i) % PMatrix => NULL()
     k = ListGetInteger( Model % BCs(i) % Values, 'Periodic BC', GotIt )
     IF( GotIt ) THEN
       Projector =>  PeriodicProjector( Model, Mesh, i, k )
       IF ( ASSOCIATED( Projector ) ) Model % BCs(i) % PMatrix => Projector
     END IF
   END DO

   ! Don't know why this is saved really...
   ! I guess because the could be higher dimensional meshes loaded already
   Model % DIMENSION = save_dim


   CALL Info('LoadMesh','Loading mesh done',Level=8)


 CONTAINS


   ! Initialize mesh structures after the size information has been 
   ! retrieved.
   !----------------------------------------------------------------
   SUBROUTINE InitializeMesh()

     INTEGER :: i,j,k,NoElems
     TYPE(Element_t), POINTER :: Element

     IF( Mesh % NumberOfNodes == 0 ) THEN
       CALL Fatal('LoadMesh','Mesh has zero nodes!')
     ELSE
       CALL Info('LoadMesh','Number of nodes in mesh: '&
           //TRIM(I2S(Mesh % NumberOfNodes)),Level=8)
     END IF
     IF( Mesh % NumberOfBulkElements == 0 ) THEN
       CALL Fatal('LoadMesh','Mesh has zero bulk elements!')
     ELSE
       CALL Info('LoadMesh','Number of bulk elements in mesh: '&
           //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=8)        
     END IF

     CALL Info('LoadMesh','Number of boundary elements in mesh: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=8)        

     Mesh % Nodes % NumberOfNodes = Mesh % NumberOfNodes          
     IF ( BoundariesOnly ) Mesh % NumberOfBulkElements = 0

     Mesh % MaxElementDOFs  = 0
     Mesh % MaxEdgeDOFs     = 0
     Mesh % MaxFaceDOFs     = 0
     Mesh % MaxBDOFs        = 0

     Mesh % DisContMesh = .FALSE.
     Mesh % DisContPerm => NULL()
     Mesh % DisContNodes = 0

     CALL Info('LoadMesh','Initial number of max element nodes: '&
         //TRIM(I2S(Mesh % MaxElementNodes)),Level=10) 

     ! Allocate the elements
     NoElems = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Elements, NoElems, 'LoadMesh' )

     DO j=1,NoElems        
       Element => Mesh % Elements(j)        

       Element % DGDOFs = 0
       Element % BodyId = 0
       Element % TYPE => NULL()
       Element % BoundaryInfo => NULL()
       Element % PDefs => NULL()
       Element % DGIndexes => NULL()
       Element % EdgeIndexes => NULL()
       Element % FaceIndexes => NULL()
       Element % BubbleIndexes => NULL()
     END DO

     ! Allocate the nodes
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Nodes % x, Mesh % NumberOfNodes, 'LoadMesh' )
     CALL AllocateVector( Mesh % Nodes % y, Mesh % NumberOfNodes, 'LoadMesh' )
     CALL AllocateVector( Mesh % Nodes % z, Mesh % NumberOfNodes, 'LoadMesh' )

   END SUBROUTINE InitializeMesh



   !------------------------------------------------------------------------------
   ! Map bodies and boundaries as prescirbed by the 'Target Bodies' and 
   ! 'Target Boundaries' keywords.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapBodiesAndBCs()

     TYPE(Element_t), POINTER :: Element
     INTEGER, ALLOCATABLE :: IndexMap(:), TmpIndexMap(:)
     INTEGER, POINTER :: Blist(:)
     INTEGER :: id,minid,maxid,body,bndry,DefaultTargetBC


     ! If "target bodies" is used map the bodies accordingly
     !------------------------------------------------------
     Found = .FALSE. 
     DO id=1,Model % NumberOfBodies
       IF( ListCheckPresent( Model % Bodies(id) % Values,'Target Bodies') ) THEN
         Found = .TRUE.
         EXIT
       END IF
     END DO

     IF( Found ) THEN
       CALL Info('LoadMesh','Remapping bodies',Level=8)      
       minid = HUGE( minid ) 
       maxid = -HUGE( maxid ) 
       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
         minid = MIN( id, minid ) 
         maxid = MAX( id, maxid )
       END DO
       IF( minid > maxid ) THEN
         CALL Fatal('LoadMesh','Body indexes are screwed!')
       END IF
       CALL Info('LoadMesh','Minimum initial body index: '//TRIM(I2S(minid)),Level=6 )
       CALL Info('LoadMesh','Maximum initial body index: '//TRIM(I2S(maxid)),Level=6 )

       minid = MIN( 1, minid ) 
       maxid = MAX( Model % NumberOfBodies, maxid ) 
       ALLOCATE( IndexMap(minid:maxid) )
       IndexMap = 0

       DO id=1,Model % NumberOfBodies
         BList => ListGetIntegerArray( Model % Bodies(id) % Values, &
             'Target Bodies', GotIt ) 
         IF ( Gotit ) THEN
           DO k=1,SIZE(BList)
             body = Blist(k)
             IF( body > maxid .OR. body < minid ) THEN
#if 0
               CALL Warn('LoadMesh','Unused body entry in > Target Bodies <  : '&
                   //TRIM(I2S(body)) )              
#endif
             ELSE IF( IndexMap( body ) /= 0 ) THEN
               CALL Warn('LoadMesh','Multiple bodies have same > Target Bodies < entry : '&
                   //TRIM(I2S(body)))
             ELSE
               IndexMap( body ) = id 
             END IF
           END DO
         ELSE
           IF( IndexMap( id ) /= 0 ) THEN
             CALL Warn('LoadMesh','Unset body already set by > Target Boundaries < : '&
                 //TRIM(I2S(id)) )
           ELSE 
             IndexMap( id ) = id
           END IF
         END IF

       END DO

       IF( .FALSE. ) THEN
         PRINT *,'Body mapping'
         DO id=minid,maxid
           IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
         END DO
       END IF

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
!        IF( IndexMap( id ) == 0 ) THEN
!          PRINT *,'Unmapped body: ',id
!          IndexMap(id) = id
!        END IF
         Element % BodyId = IndexMap( id ) 
       END DO

       DEALLOCATE( IndexMap )
     ELSE
       CALL Info('LoadMesh','Skipping remapping of bodies',Level=10)      
     END IF


     IF( Mesh % NumberOfBoundaryElements == 0 ) RETURN

     ! Target boundaries are usually given so this is not conditional
     !---------------------------------------------------------------
     CALL Info('LoadMesh','Remapping boundaries',Level=8)      
     minid = HUGE( minid ) 
     maxid = -HUGE( maxid ) 
     DO i=Mesh % NumberOfBulkElements+1,&
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       id = Element % BoundaryInfo % Constraint
       minid = MIN( id, minid ) 
       maxid = MAX( id, maxid )
     END DO


     CALL Info('LoadMesh','Minimum initial boundary index: '//TRIM(I2S(minid)),Level=6 )
     CALL Info('LoadMesh','Maximum initial boundary index: '//TRIM(I2S(maxid)),Level=6 )
     IF( minid > maxid ) THEN
       CALL Fatal('LoadMesh','Boundary indexes are screwed')
     END IF

     minid = MIN( minid, 1 ) 
     maxid = MAX( maxid, Model % NumberOfBCs ) 
     ALLOCATE( IndexMap(minid:maxid) )
     IndexMap = 0


     DO j=1,Model % NumberOfBoundaries
       id = ListGetInteger( Model % Boundaries(j) % Values, &
           'Boundary Condition',GotIt, minv=1, maxv=Model % NumberOFBCs )
       IF( id == 0 ) CYCLE
       bndry = Model % BoundaryId(j)
       IF( bndry > maxid ) THEN
         CALL Warn('LoadMesh','BoundaryId exceeds range')
       ELSE IF( bndry == 0 ) THEN
         CALL Warn('LoadMesh','BoundaryId is zero')
       ELSE
         IndexMap( bndry ) = id
       END IF
     END DO

     DefaultTargetBC = 0
     DO id=1,Model % NumberOfBCs
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default Target', GotIt)) DefaultTargetBC = id       
       BList => ListGetIntegerArray( Model % BCs(id) % Values, &
           'Target Boundaries', GotIt ) 
       IF ( Gotit ) THEN
         DO k=1,SIZE(BList)
           bndry = Blist(k)
           IF( bndry > maxid ) THEN
#if 0
  in my opinion, this is quite usual ... Juha
             CALL Warn('LoadMesh','Unused BC entry in > Target Boundaries <  : '&
                 //TRIM(I2S(bndry)) )              
#endif
           ELSE IF( IndexMap( bndry ) /= 0 ) THEN
             CALL Warn('LoadMesh','Multiple BCs have same > Target Boundaries < entry : '&
                 //TRIM(I2S(bndry)) )
           ELSE 
             IndexMap( bndry ) = id 
           END IF
         END DO
       ELSE
         IF( IndexMap( id ) /= 0 .AND. id /= DefaultTargetBC ) THEN
           CALL Warn('LoadMesh','Unset BC already set by > Target Boundaries < : '&
               //TRIM(I2S(id)) )
         ELSE 
           ! IndexMap( id ) = id
         END IF
       END IF
     END DO

     IF( .FALSE. ) THEN
       PRINT *,'Boundary mapping'
       DO id=minid,maxid
         IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
       END DO
     END IF

     IF( DefaultTargetBC /= 0 ) THEN
       CALL Info('LoadMesh','Default Target BC: '&
           //TRIM(I2S(DefaultTargetBC)),Level=8)
     END IF


     DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       n = Element % TYPE % NumberOfNodes
       bndry = Element % BoundaryInfo % Constraint 

       IF( bndry > maxid .OR. bndry < minid ) THEN
         CALL Warn('LoadMesh','Boundary index '//TRIM(I2S(bndry))&
             //' not in range: '//TRIM(I2S(minid))//','//TRIM(I2S(maxid)) )
       END IF

       IF( IndexMap( bndry ) < 0 ) THEN
         Element % BoundaryInfo % Constraint = 0
         CYCLE

       ELSE IF( IndexMap( bndry ) == 0 ) THEN
         IF( DefaultTargetBC /= 0 ) THEN
!          PRINT *,'Default boundary map: ',bndry,DefaultTargetBC
           IndexMap( bndry ) = DefaultTargetBC
         ELSE 
!          IF( bndry <= Model % NumberOfBCs ) THEN            
!            PRINT *,'Unmapped boundary: ',bndry
!          ELSE
!            PRINT *,'Unused boundary: ',bndry
!          END IF
           IndexMap( bndry ) = -1 
           Element % BoundaryInfo % Constraint = 0           
           CYCLE
         END IF
       END IF

       bndry = IndexMap( bndry ) 
       Element % BoundaryInfo % Constraint = bndry 

       IF( bndry <= Model % NumberOfBCs ) THEN
         Element % BodyId  = ListGetInteger( &
             Model % BCs(bndry) % Values, 'Body Id', Gotit, 1, Model % NumberOfBodies )
         Element % BoundaryInfo % OutBody = &
             ListGetInteger( Model % BCs(bndry) % Values, &
             'Normal Target Body', GotIt, maxv=Model % NumberOFBodies ) 
       END IF
     END DO

     DEALLOCATE( IndexMap ) 

   END SUBROUTINE MapBodiesAndBCs


   ! Check for the non-nodal element basis
   !--------------------------------------------------------
   SUBROUTINE NonNodalElements()

     INTEGER, POINTER :: EdgeDofs(:), FaceDofs(:)
     INTEGER :: DGIndex, body_id, body_id0, eq_id, solver_id, el_id
     LOGICAL :: NeedEdges, Found, FoundDef0, FoundDef, FoundEq, GotIt, MeshDeps, &
                FoundEqDefs, FoundSolverDefs(Model % NumberOfSolvers), FirstOrderElements
     TYPE(Element_t), POINTER :: Element
     TYPE(ValueList_t), POINTER :: Vlist
     INTEGER :: inDOFs(10,6)
     CHARACTER(MAX_NAME_LEN) :: ElementDef0, ElementDef

     EdgeDOFs => NULL()
     CALL AllocateVector( EdgeDOFs, Mesh % NumberOfBulkElements, 'LoadMesh' )
     FaceDOFs => NULL()
     CALL AllocateVector( FaceDOFs, Mesh % NumberOfBulkElements, 'LoadMesh' )

     DGIndex = 0
     NeedEdges = .FALSE.

     InDofs = 0
     InDofs(:,1) = 1
     IF ( PRESENT(Def_Dofs) ) THEN
       inDofs = Def_Dofs
     END IF

     ! P-basis only over 1st order elements:
     ! -------------------------------------
     FirstOrderElements = .TRUE.
     DO i=1,Mesh % NumberOfBulkElements
       IF (Mesh % Elements(i) % Type % BasisFunctionDegree>1) THEN
         FirstOrderElements = .FALSE.; EXIT
       END IF
     END DO

    !
    ! Check whether the "Element" definitions can depend on mesh
    ! -----------------------------------------------------------
    MeshDeps = .FALSE.; FoundEqDefs = .FALSE.;  FoundSolverDefs = .FALSE.

    DO eq_id=1,Model % NumberOFEquations
      Vlist => Model % Equations(eq_id) % Values
      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )
      FoundEqDefs = FoundEqDefs .OR. FoundDef0
      j = INDEX(ElementDef0,'p:')
      IF (j>0.AND. ElementDef0(j+2:j+2)=='%') MeshDeps = .TRUE.
    END DO

    DO solver_id=1,Model % NumberOFSolvers
      Vlist => Model % Solvers(solver_id) % Values

      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)
      FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef0

      ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef0)
      FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef0

      j = INDEX(ElementDef0,'p:')
      IF (j>0.AND. ElementDef0(j+2:j+2)=='%') meshdeps = .TRUE.
    END DO

    IF(.NOT.MeshDeps) THEN
      ElementDef = ' '
      FoundDef0 = .FALSE.
      DO body_id=1,Model % NumberOfBodies
        ElementDef0 = ' '
        Vlist => Model % Bodies(body_id) % Values
        eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
        IF( FoundEq ) THEN
          Vlist => Model % Equations(eq_id) % Values
          IF(FoundEqDefs) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

          DO solver_id=1,Model % NumberOfSolvers
            FoundDef = .FALSE.
            IF(FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)
 
            IF ( FoundDef ) THEN
              CALL GetMaxDefs( Model, Mesh, Element, ElementDef, solver_id, body_id, Indofs )
            ELSE
              IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) &
                 ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

              CALL GetMaxDefs( Model, Mesh, Element, ElementDef0, solver_id, body_id, Indofs )

              IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
            END IF
          END DO
        END IF
      END DO
    END IF

     ! non-nodal elements in bulk elements
     !------------------------------------------------------------
     body_id0 = -1; FoundDef=.FALSE.; FoundEq=.FALSE.
     ElementDef = ' '

     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       body_id = Element % BodyId
       n = Element % TYPE % NumberOfNodes

       ! Check the Solver specific element types
       IF( Meshdeps ) THEN
         IF ( body_id/=body_id0 ) THEN
           Vlist => Model % Bodies(body_id) % Values
           eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
         END IF

         ElementDef0 = ' '
         IF( FoundEq ) THEN
           Vlist => Model % Equations(eq_id) % Values
           IF( FoundEqDefs.AND.body_id/=body_id0 ) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

           DO solver_id=1,Model % NumberOfSolvers
             FoundDef = .FALSE.
             IF (FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)

             IF ( FoundDef ) THEN
               CALL GetMaxDefs( Model, Mesh, Element, ElementDef, solver_id, body_id, Indofs )
             ELSE
               IF(.NOT. FoundDef0.AND.FoundSolverDefs(solver_id)) &
                  ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

               CALL GetMaxDefs( Model, Mesh, Element, ElementDef0, solver_id, body_id, Indofs )

               IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
             END IF
           END DO
         END IF
         body_id0 = body_id
      END IF

       el_id = Element % TYPE % ElementCode / 100

       ! Apply the elementtypes
       IF ( inDOFs(el_id,1) /= 0 ) THEN
         Element % NDOFs = n
       ELSE
         Element % NDOFs = 0
       END IF

       EdgeDOFs(i) = MAX(0,inDOFs(el_id,2))
       FaceDOFs(i) = MAX(0,inDOFs(el_id,3))

       IF ( PRESENT(Def_Dofs) ) THEN
         IF ( Def_Dofs(el_id,4) == 0 ) inDOFs(el_id,4) = n
       END IF

       NULLIFY( Element % DGIndexes )
       IF ( inDOFs(el_id,4) > 0 ) THEN
         CALL AllocateVector( Element % DGIndexes, inDOFs(el_id,4))
         DO j=1,inDOFs(el_id,4)
           DGIndex = DGIndex + 1
           Element % DGIndexes(j) = DGIndex
         END DO
       ELSE
         NULLIFY( Element % DGIndexes )
       END IF
       Element % DGDOFs = MAX(0,inDOFs(el_id,4))
       NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:4)>0 )

       ! Check if given element is a p element
       IF (FirstOrderElements.AND.inDOFs(el_id,6) > 0) THEN
         CALL AllocatePDefinitions(Element)

         NeedEdges = .TRUE.

         ! Calculate element bubble dofs and set element p
         Element % PDefs % P = inDOFs(el_id,6)
         IF ( inDOFs(el_id,5) > 0 ) THEN
           Element % BDOFs = inDOFs(el_id,5)
         ELSE
           Element % BDOFs = getBubbleDOFs(Element, Element % PDefs % P)
         END IF

         ! All elements in actual mesh are not edges
         Element % PDefs % pyramidQuadEdge = .FALSE.
         Element % PDefs % isEdge = .FALSE.

         ! If element is of type tetrahedron and is a p element, 
         ! do the Ainsworth & Coyle trick
         IF (Element % TYPE % ElementCode == 504) CALL ConvertToACTetra(Element)
         CALL GetRefPElementNodes( Element,  Element % TYPE % NodeU, &
             Element % TYPE % NodeV, Element % TYPE % NodeW )
       ELSE 
         ! Clear P element definitions and set manual bubbles
         Element % PDefs => NULL()
         Element % BDOFs = MAX(0,inDOFs(el_id,5))
         ! WRITE (*,*) Element % BDOFs
       END IF

       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes,Element % TYPE % NumberOfNodes )
     END DO

     ! non-nodal elements in boundary elements
     !------------------------------------------------------------    
     DO i = Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('NonNodalElements','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
         CALL Fatal('NonNodalElements','Type in Element '//TRIM(I2S(i))//' not associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       Element % NDOFs  = n
       el_id = ELement % TYPE % ElementCode / 100

       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
         IF( Element % BoundaryInfo % Left % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Left % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Left % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
         IF ( Element % BoundaryInfo % Right % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Right % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Right % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF
     END DO

     IF ( Mesh % MaxElementDOFs <= 0 ) Mesh % MaxElementDOFs = Mesh % MaxElementNodes 

     IF ( NeedEdges ) THEN
       CALL Info('NonNodalElements','Requested elements require creation of edges',Level=8)
       CALL SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs)
     END IF

     CALL SetMeshMaxDOFs(Mesh)

     IF( ASSOCIATED(EdgeDOFs) ) DEALLOCATE(EdgeDOFs )
     IF( ASSOCIATED(FaceDOFs) ) DEALLOCATE(FaceDOFs)

   END SUBROUTINE NonNodalElements


   !------------------------------------------------------------------------------
   ! Map and scale coordinates, and increase the size of the coordinate
   ! vectors, if requested.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapCoordinates()

     REAL(KIND=dp), POINTER :: NodesX(:), NodesY(:), NodesZ(:), Wrk(:,:)
     INTEGER, POINTER :: CoordMap(:)
     REAL(KIND=dp) :: CoordScale(3)
     INTEGER :: mesh_dim, model_dim

     ! Perform coordinate mapping
     !------------------------------------------------------------
     CoordMap => ListGetIntegerArray( Model % Simulation, &
         'Coordinate Mapping',GotIt )
     IF ( GotIt ) THEN
       CALL Info('LoadMesh','Performing coordinate mapping',Level=8)

       IF ( SIZE( CoordMap ) /= 3 ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'LoadMesh', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'LoadMesh', Message )
       END IF

       IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'LoadMesh', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'LoadMesh', Message )
       END IF

       IF( CoordMap(1) == 1 ) THEN
         NodesX => Mesh % Nodes % x
       ELSE IF( CoordMap(1) == 2 ) THEN
         NodesX => Mesh % Nodes % y
       ELSE
         NodesX => Mesh % Nodes % z
       END IF

       IF( CoordMap(2) == 1 ) THEN
         NodesY => Mesh % Nodes % x
       ELSE IF( CoordMap(2) == 2 ) THEN
         NodesY => Mesh % Nodes % y
       ELSE
         NodesY => Mesh % Nodes % z
       END IF

       IF( CoordMap(3) == 1 ) THEN
         NodesZ => Mesh % Nodes % x
       ELSE IF( CoordMap(3) == 2 ) THEN
         NodesZ => Mesh % Nodes % y
       ELSE
         NodesZ => Mesh % Nodes % z
       END IF

       Mesh % Nodes % x => NodesX
       Mesh % Nodes % y => NodesY
       Mesh % Nodes % z => NodesZ
     END IF

     ! Determine the mesh dimension 
     !----------------------------------------------------------------------------
     mesh_dim = 0
     model_dim = 0
     IF ( ANY( Mesh % Nodes % x /= Mesh % Nodes % x(1) ) ) THEN
       model_dim = 1
       mesh_dim = mesh_dim + 1
     END IF
     IF ( ANY( Mesh % Nodes % y /= Mesh % Nodes % y(1) ) ) THEN
       model_dim = 2
       mesh_dim = mesh_dim + 1
     END IF
     IF ( ANY( Mesh % Nodes % z /= Mesh % Nodes % z(1) ) ) THEN
       model_dim = 3
       mesh_dim = mesh_dim + 1
     END IF

     Mesh % MeshDim = mesh_dim

     save_dim = Model % DIMENSION
     IF ( Model % DIMENSION <= 0 ) Model % DIMENSION = model_dim

     CALL Info('LoadMesh','Dimension of model is: '//TRIM(I2S(model_dim)),Level=8)
     CALL Info('LoadMesh','Dimension of mesh is: '//TRIM(I2S(mesh_dim)),Level=8)

     ! Scaling of coordinates
     !-----------------------------------------------------------------------------
     Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',GotIt )    
     IF( GotIt ) THEN            
       CoordScale = 1.0_dp
       DO i=1,mesh_dim
         j = MIN( i, SIZE(Wrk,1) )
         CoordScale(i) = Wrk(j,1)
       END DO
       WRITE(Message,'(A,3ES10.3)') 'Scaling coordinates:',CoordScale(1:mesh_dim)
       CALL Info('LoadMesh',Message) 
       Mesh % Nodes % x = CoordScale(1) * Mesh % Nodes % x
       IF( mesh_dim > 1) Mesh % Nodes % y = CoordScale(2) * Mesh % Nodes % y
       IF( mesh_dim > 2) Mesh % Nodes % z = CoordScale(3) * Mesh % Nodes % z
     END IF

   END SUBROUTINE MapCoordinates




   ! When the parallel nodal neighbours have been found 
   ! perform numbering for face and edge elements as well.
   !-------------------------------------------------------------------    
   SUBROUTINE ParallelNonNodalElements()

     TYPE(Element_t), POINTER :: Element
     IF(.NOT. Parallel ) RETURN

     n = SIZE( Mesh % ParallelInfo % NeighbourList )

     ! For unset neighbours just set the this partition to be the only owner
     DO i=1,n
       IF (.NOT.ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
         CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(i) % Neighbours,1)
         Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = mype
       END IF
     END DO

     ! Create parallel numbering of faces
     CALL SParFaceNumbering(Mesh)
     DO i=1,Mesh % NumberOfFaces
       Mesh % MaxFaceDOFs = MAX(Mesh % MaxFaceDOFs,Mesh % Faces(i) % BDOFs)
     END DO

     ! Create parallel numbering for edges
     CALL SParEdgeNumbering(Mesh)
     DO i=1,Mesh % NumberOfEdges
       Mesh % MaxEdgeDOFs = MAX(Mesh % MaxEdgeDOFs,Mesh % Edges(i) % BDOFs)
     END DO

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)        
       Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
           Element % TYPE % NumberOfNodes + &
           Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
           Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
           Element % BDOFs, &
           Element % DGDOFs )
     END DO

   END SUBROUTINE ParallelNonNodalElements

   !------------------------------------------------------------------------------
 END FUNCTION LoadMesh2
 !------------------------------------------------------------------------------




!------------------------------------------------------------------------------
  SUBROUTINE SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs)
!------------------------------------------------------------------------------
    INTEGER, OPTIONAL :: EdgeDOFs(:), FaceDOFs(:)
    TYPE(Mesh_t) :: Mesh
    INTEGER, OPTIONAL :: indofs(:,:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,el_id
    TYPE(Element_t), POINTER :: Element, Edge, Face
!------------------------------------------------------------------------------

    CALL FindMeshEdges(Mesh)

    ! Set edge and face polynomial degree and degrees of freedom for
    ! all elements
    DO i=1,Mesh % NumberOFBulkElements
       Element => Mesh % Elements(i)

       ! Iterate each edge of element
       DO j = 1,Element % TYPE % NumberOfEdges
          Edge => Mesh % Edges( Element % EdgeIndexes(j) ) 
          
          ! Set attributes of p element edges
          IF ( ASSOCIATED(Element % PDefs) ) THEN   
             ! Set edge polynomial degree and dofs
             Edge % PDefs % P = MAX( Element % PDefs % P, Edge % PDefs % P)
             Edge % BDOFs = MAX(Edge % BDOFs, Edge % PDefs % P - 1)
             Edge % PDefs % isEdge = .TRUE.
             ! Get gauss points for edge. If no dofs 2 gauss points are 
             ! still needed for integration of linear equation!
             Edge % PDefs % GaussPoints = (Edge % BDOFs+2)**Edge % TYPE % DIMENSION  

             IF (ASSOCIATED(Edge % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Right, Mesh)
             END IF
             
          ! Other element types, which need edge dofs
          ELSE IF(PRESENT(EdgeDOFs)) THEN
             Edge % BDOFs = MAX(EdgeDOFs(i), Edge % BDOFs)
          END IF

          ! Get maximum dof for edges
          Mesh % MaxEdgeDOFs = MAX(Edge % BDOFs, Mesh % MaxEdgeDOFs)
       END DO

       ! Iterate each face of element
       DO j=1,Element % TYPE % NumberOfFaces
          Face => Mesh % Faces( Element % FaceIndexes(j) )

          ! Set attibutes of p element faces
          IF ( ASSOCIATED(Element % PDefs) ) THEN
             ! Set face polynomial degree and dofs
             Face % PDefs % P = MAX(Element % PDefs % P, Face % PDefs % P)
             ! Get number of face dofs
             Face % BDOFs = MAX( Face % BDOFs, getFaceDOFs(Element, Face % PDefs % P, j) )
             Face % PDefs % isEdge = .TRUE.
             Face % PDefs % GaussPoints = getNumberOfGaussPointsFace( Face, Mesh )
             IF (ASSOCIATED(Face % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Right, Mesh)
             END IF
          ELSE IF (PRESENT(FaceDOFs)) THEN
             el_id = face % TYPE % ElementCode / 100
             Face % BDOFs = MAX(FaceDOFs(i), Face % BDOFs)
             IF ( PRESENT(inDOFs) ) Face % BDOFs = MAX(Face % BDOFs, InDOFs(el_id+6,5))
          END IF
             
          ! Get maximum dof for faces
          Mesh % MaxFaceDOFs = MAX(Face % BDOFs, Mesh % MaxFaceDOFs)
       END DO
    END DO

    ! Set local edges for boundary elements
    DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)

       ! Here set local number and copy attributes to this boundary element for left parent.
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Left % PDefs)) THEN
            CALL AllocatePDefinitions(Element)
            Element % PDefs % isEdge = .TRUE.
            CALL AssignLocalNumber(Element, Element % BoundaryInfo % Left, Mesh)
            ! CYCLE
          END IF
       END IF

       ! Here set local number and copy attributes to this boundary element for right parent
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Right % PDefs)) THEN
             CALL AllocatePDefinitions(Element)
             Element % PDefs % isEdge = .TRUE.
             CALL AssignLocalNumber(Element, Element % BoundaryInfo % Right, Mesh)
          END IF
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMeshEdgeFaceDofs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE SetMeshMaxDOFs(Mesh)
!------------------------------------------------------------------------------
   TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
   TYPE(Element_t), POINTER :: Element
   INTEGER :: i,j,n

   ! Set gauss points for each p element
   DO i=1,Mesh % NumberOfBulkElements
     Element => Mesh % Elements(i)
     IF ( ASSOCIATED(Element % PDefs) ) THEN
       Element % PDefs % GaussPoints = getNumberOfGaussPoints( Element, Mesh )
     END IF

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
          Element % TYPE % NumberOfNodes + &
          Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
          Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
          Element % BDOFs, &
          Element % DGDOFs )

     Mesh % MaxBDOFs = MAX( Element % BDOFs, Mesh % MaxBDOFs )
   END DO

   DO i=1,Mesh % NumberOFBulkElements
     Element => Mesh % Elements(i)
     IF ( Element % BDOFs > 0 ) THEN
       ALLOCATE( Element % BubbleIndexes(Element % BDOFs) )
       DO j=1,Element % BDOFs
         Element % BubbleIndexes(j) = Mesh % MaxBDOFs*(i-1)+j
       END DO
     END IF
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE SetMeshMaxDOFs
!------------------------------------------------------------------------------
 
 SUBROUTINE ReadTargetNames(Model,Filename)
     CHARACTER(LEN=*) :: FileName
     TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')
   INTEGER :: i,j,k,iostat,i1,i2,i3,n
   INTEGER :: ivals(256)
   CHARACTER(LEN=1024) :: str, name0, name1
   TYPE(ValueList_t), POINTER :: Vlist
   LOGICAL :: Found, AlreadySet

   OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT=iostat )
   IF( iostat /= 0 ) THEN
     RETURN
   ELSE
     CALL Info('ReadTargetNames','Reading names info from file: '//TRIM(FileName))
   END IF

   DO WHILE( .TRUE. ) 
     READ(FileUnit,'(A)',IOSTAT=iostat) str
     IF( iostat /= 0 ) EXIT
     i = INDEX( str,'$')     
     j = INDEX( str,'=')
     IF( i == 0 .OR. j == 0 ) CYCLE

     i = i + 1
     DO WHILE(i<=LEN_TRIM(str) .AND. str(i:i)==' ')
       i = i + 1
     END DO     
     
     i1 = i
     i2 = j-1
     i3 = j+1

     ! Move to lowercase since the "name" in sif file is also
     ! always in lowercase. 
     DO i=i1,i2
       j = i+1-i1
       k = ICHAR(str(i:i))
       IF ( k >= A .AND. k<= Z ) THEN
         name0(j:j) = CHAR(k+U2L)
       ELSE
         name0(j:j) = str(i:i)
       END IF
     END DO

     n = str2ints( str(i3:),ivals )
     IF( n == 0 ) THEN
       CALL Fatal('ReadTargetNames','Could not find arguments for: '//str(i1:i2))
     END IF

     AlreadySet = .FALSE.

     DO i=1,Model % NumberOfBCs
       Vlist => Model % BCs(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches BC '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Boundaries') ) THEN
           CALL Info('ReadTargetNames','> Target Boundaries < already defined for BC '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Boundaries',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO

     DO i=1,Model % NumberOfBodies
       Vlist => Model % Bodies(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches body '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Bodies') ) THEN
           CALL Info('ReadTargetNames','> Target Bodies < already defined for Body '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Bodies',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO
     
     IF(.NOT. AlreadySet ) THEN
       CALL Warn('ReadTargetNames','Could not map name to Body nor BC: '//name0(1:i2-i1+1) )
     END IF

   END DO

   CLOSE(FileUnit)
   
 END SUBROUTINE ReadTargetNames



!------------------------------------------------------------------------------
  SUBROUTINE ReadElementPropertyFile(FileName,Mesh)
     CHARACTER(LEN=*) :: FileName
     TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MAXLEN=1024
    CHARACTER(LEN=:), ALLOCATABLE :: str
    INTEGER :: i,j,n
    INTEGER, PARAMETER :: FileUnit = 10
    REAL(KIND=dp) :: x
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementData_t), POINTER :: PD,PD1

    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)

    OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', ERR=10 )

    DO WHILE( ReadAndTrim(FileUnit,str) )
      READ( str(9:),*) i
      IF ( i < 0 .OR. i > Mesh % NumberOFBulkElements ) THEN
        CALL Fatal( 'ReadElementProperties', 'Element id out of range.' )
      END IF
      
      IF ( SEQL( str, 'element:') ) THEN
        Element => Mesh % Elements(i)
        PD => Element % PropertyData

        DO WHILE(ReadAndTrim(FileUnit,str))
          IF ( str == 'end' ) EXIT

          i = INDEX(str, ':')
          IF ( i<=0 ) CYCLE

          IF ( .NOT.ASSOCIATED(PD)  ) THEN
            ALLOCATE( Element % PropertyData )
            PD => Element % PropertyData
            PD % Name = TRIM(str(1:i-1))
          ELSE
            DO WHILE(ASSOCIATED(PD))
              IF ( PD % Name==TRIM(str(1:i-1)) ) EXIT
              PD1 => PD
              PD => PD % Next
            END DO
            
            IF (.NOT. ASSOCIATED(PD) ) THEN
              ALLOCATE(PD1 % Next)
              PD => PD1 % Next
              PD % Name = TRIM(str(1:i-1))
            END IF
          END IF

          j = i+1
          n = 0
          DO WHILE(j<=LEN_TRIM(str))
            READ( str(j:), *, END=20,ERR=20 ) x
            n = n + 1
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
              j = j + 1
            END DO
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
              j = j + 1
            END DO
          END DO
20        CONTINUE
          IF ( n>0 ) THEN
            ALLOCATE(PD % Values(n))
            j = i+1
            n = 1
            DO WHILE(j<=LEN_TRIM(str))
              READ( str(j:), *, END=30,ERR=30 ) PD % Values(n)
              n = n + 1
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
                j = j + 1
              END DO
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
                j = j + 1
              END DO
            END DO
30          CONTINUE
          END IF
        END DO
      END IF
    END DO

    CLOSE(FileUnit)

10  CONTINUE

!------------------------------------------------------------------------------
  END SUBROUTINE ReadElementPropertyFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE MeshStabParams( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: i,n, istat
    LOGICAL :: stat
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------


    CALL Info('MeshStabParams','Computing stabilization parameters',Level=7)
    CALL ResetTimer('MeshStabParams')

    DO i=1,CurrentModel % NumberOfSolvers
       Solver => CurrentModel % Solvers(i)
       IF ( ASSOCIATED( Mesh, Solver % Mesh ) ) THEN
          Mesh % Stabilize = Mesh % Stabilize .OR. &
             ListGetLogical( Solver % Values, 'Stabilize', Stat )
          Mesh % Stabilize = Mesh % Stabilize .OR. &
             ListGetString( Solver % Values,  &
                     'Stabilization Method', Stat )=='vms'
          Mesh % Stabilize = Mesh % Stabilize .OR. &
             ListGetString( Solver % Values,  &
                     'Stabilization Method', Stat )=='stabilized'
       END IF
    END DO

    CALL AllocateVector( Nodes % x, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % y, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % z, Mesh % MaxElementNodes )

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)
       n = Element % TYPE % NumberOfNodes
       Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
       Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
       Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
       IF ( Mesh % Stabilize ) THEN
          CALL StabParam( Element, Nodes,n, &
              Element % StabilizationMK, Element % hK )
       ELSE
          Element % hK = ElementDiameter( Element, Nodes )
       END IF
    END DO
 
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

    CALL CheckTimer('MeshStabParams',Level=7,Delete=.TRUE.)
!----------------------------------------------------------------------------
  END SUBROUTINE MeshStabParams
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Given two interface meshes check the angle between them using the normal
!> vectors of the first element. Also check that all other elements are
!> aligned with the first one. Only then is it possible to determine the angle.
!------------------------------------------------------------------------------
  SUBROUTINE CheckInterfaceMeshAngle(BMesh1, BMesh2, Angles, GotAngles) 
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    REAL(KIND=dp) :: Angles(3)
    LOGICAL :: GotAngles
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: Normal(3), Normal1(3), Normal2(3), Dot1Min, Dot2Min, Alpha
    LOGICAL :: ConstantNormals

    ! Currently check of the normal direction is not enforced since at this stage 
    ! CurrentModel % Nodes may not exist!
    ! This means that there may be a 180 error in the directions. 
    ! Therefore an angle smaller than 180 is always chosen.
    !-----------------------------------------------------------------------------
    N = MAX( BMesh1 % MaxElementNodes, BMesh2 % MaxElementNodes )
    ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
 
    DO k=1,2
      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
      !-------------------------------------------------------------------------
      DO i=1, PMesh % NumberOfBoundaryElements
        Element => PMesh % Elements(i)
        
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
        
        Normal = NormalVector( Element, ElementNodes, Check = .FALSE. ) 

        ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
        !-------------------------------------------------------------------------       
        IF( i == 1 ) THEN
          Normal2 = Normal
          Dot2Min = 1.0_dp
        ELSE
          Dot2min = MIN( Dot2Min, SUM( Normal * Normal2 ) )
        END IF
      END DO

      IF( k == 1 ) THEN
        Normal1 = Normal2 
        Dot1Min = Dot2Min
      END IF
    END DO

    ConstantNormals = ( 1 - Dot1Min < 1.0d-6 ) .AND. ( 1 - Dot2Min < 1.0e-6 )     
    IF( ConstantNormals ) THEN
      WRITE(Message,'(A,3ES12.3)') 'Master normal: ',Normal1
      CALL Info('CheckInterfaceMeshAngle',Message,Level=8)    
      
      WRITE(Message,'(A,3ES12.3)') 'Initial Target normal: ',Normal2
      CALL Info('CheckInterfaceMeshAngle',Message,Level=8)    
            
      ! The full angle between the two normals
      Alpha = ACOS( SUM( Normal1 * Normal2 ) ) * 180 / PI
      WRITE(Message,'(A,ES12.3)') &
          'Suggested angle between two normals in degs (+/- 180): ',Alpha 
      CALL Info('CheckInterfaceMeshAngle',Message,Level=8)
    ELSE
      CALL Warn('CheckInterfaceMeshAngle','Could not suggest rotation angle')
    END IF


    GotAngles = .FALSE.
    Angles = 0.0_dp
    IF( .NOT. ConstantNormals ) THEN
      CALL Warn('CheckInterfaceMeshAngle','Normals are not constant, cannot test for rotation!')
    ELSE IF( Alpha > EPSILON( Alpha ) ) THEN
      ! Rotation should be performed 
      DO i=1,3
        IF( ABS ( Normal1(i) - Normal2(i) ) < EPSILON( Alpha ) ) THEN
          GotAngles = .TRUE.            
          WRITE(Message,'(A,I0,A,ES12.3)') &
              'Rotation around axis ',i,' in degs ',Alpha 
          CALL Info('CheckInterfaceMeshAngle',Message,Level=8)
          Angles(i) = Alpha
          EXIT
        END IF
      END DO
      IF(.NOT. GotAngles ) THEN
        CALL Warn('CheckInterfaceMeshAngle','could not define rotation axis, improve algorithm!')
      END IF
    END IF

    DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z )
    
  END SUBROUTINE CheckInterfaceMeshAngle
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Create master and slave mesh for the interface in order to at a later 
!> stage create projector matrix to implement periodicity or mortar elements.
!> The idea is to use a reduced set of elements and thereby speed up the 
!> mapping process. Also this gives more flexibility in transformation
!> operations since the nodes may be ereased after use. 
!------------------------------------------------------------------------------
  SUBROUTINE CreateInterfaceMeshes( Model, Mesh, This, Trgt, BMesh1, BMesh2, &
      Success ) 
!------------------------------------------------------------------------------    
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Success
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,n1,n2,k1,k2,ind,Constraint,DIM
    TYPE(Element_t), POINTER :: Element, Left, Right, Elements(:)
    LOGICAL :: ThisActive, TargetActive
    INTEGER, POINTER :: NodeIndexes(:), Perm1(:), Perm2(:), PPerm(:)
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    LOGICAL :: OnTheFlyBC, CheckForHalo, NarrowHalo, NoHalo, Found

    TYPE(Element_t), POINTER :: Parent,q
    INTEGER :: en, in, HaloCount, ActiveCount
    LOGICAL, ALLOCATABLE :: ActiveNode(:)

    CALL Info('CreateInterfaceMeshes','Making a list of elements at interface',Level=9)

    IF ( This <= 0 .OR. Trgt <= 0 ) THEN
      CALL Fatal('CreateInterfaceMeshes','Invalid target boundaries')
    END IF

    ! Interface meshes consist of boundary elements only    
    Elements => Mesh % Elements( Mesh % NumberOfBulkElements+1: )

    ! If the target is larger than number of BCs givem then 
    ! it has probably been created on-the-fly from a discontinuous boundary.
    OnTheFlyBC = ( Trgt > Model % NumberOfBCs )

    ! In parallel we may have some excess halo elements. 
    ! To eliminate them mark the nodes that are associated to elements truly owned. 
    NarrowHalo = .FALSE.
    NoHalo = .FALSE.

    IF( ParEnv % PEs > 1 ) THEN
      ! Account for halo elements that share some nodes for the master boundary
      NarrowHalo = ListGetLogical(Model % Solver % Values,'Projector Narrow Halo',Found)

      ! Do not allow for any halo elements for the master boundary
      IF( .NOT. Found ) THEN
        NoHalo = ListGetLogical(Model % Solver % Values,'Projector No Halo',Found)
      END IF
      
      IF(.NOT. Found ) THEN
        IF( ListGetLogical(Model % Solver % Values, 'Partition Local Constraints',Found) ) THEN
          NarrowHalo = .TRUE.
        ELSE
          NoHalo = .TRUE.
        END IF
      END IF
    END IF

    ! This is just temporarily set to false always until the logic has been tested. 
    CheckForHalo = NarrowHalo .OR. NoHalo

    IF( CheckForHalo ) THEN
      CALL Info('CreateInterfaceMeshes','Checking for halo elements',Level=15)
      ALLOCATE( ActiveNode( Mesh % NumberOfNodes ) )
      HaloCount = 0
      ActiveNode = .FALSE.
      DO i=1, Mesh % NumberOfBoundaryElements
        Element => Elements(i)
        IF (Element % TYPE % ElementCode<=200) CYCLE

        Left => Element % BoundaryInfo % Left 
        IF( ASSOCIATED( Left ) ) THEN
          IF( Left % PartIndex == ParEnv % MyPe ) THEN
            ActiveNode( Left % NodeIndexes ) = .TRUE.
          ELSE
            HaloCount = HaloCount + 1
          END IF
        END IF

        Right => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Right ) ) THEN
          IF( Right % PartIndex == ParEnv % MyPe ) THEN
            ActiveNode( Right % NodeIndexes ) = .TRUE.
          ELSE
            HaloCount = HaloCount + 1 
          END IF
        END IF
      END DO

      ! No halo element found on the boundary so no need to check them later
      IF( HaloCount == 0 ) THEN
        CALL Info('CreateInterfaceMeshes','Found no halo elements to eliminate',Level=15)
        DEALLOCATE( ActiveNode ) 
        CheckForHalo = .FALSE.
      ELSE
        CALL Info('CreateInterfaceMeshes','Number of halo elements to eliminate: '&
            //TRIM(I2S(HaloCount)),Level=12)
      END IF
    END IF


!   Search elements in this boundary and its periodic
!   counterpart:
!   --------------------------------------------------
    n1 = 0
    n2 = 0
    HaloCount = 0
    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      IF (Element % TYPE % ElementCode<=200) CYCLE

      Constraint = Element % BoundaryInfo % Constraint
      IF( Model % BCs(This) % Tag == Constraint ) THEN
        IF( CheckForHalo ) THEN
          IF( NarrowHalo ) THEN
            IF( ANY(ActiveNode(Element % NodeIndexes) ) ) THEN
              n1 = n1 + 1
            ELSE
              HaloCount = HaloCount + 1
            END IF
          ELSE IF( NoHalo ) THEN
            ThisActive = .FALSE.
            Left => Element % BoundaryInfo % Left 
            IF( ASSOCIATED( Left ) ) THEN
              ThisActive = ( Left % PartIndex == ParEnv % MyPe )
            END IF
            Right => Element % BoundaryInfo % Right
            IF( ASSOCIATED( Right ) ) THEN
              ThisActive = ThisActive .OR. &
                  ( Right % PartIndex == ParEnv % MyPe ) 
            END IF
            IF( ThisActive ) THEN
              n1 = n1 + 1
            ELSE
              HaloCount = HaloCount + 1
            END IF
          END IF
        ELSE
          n1 = n1 + 1
        END IF
      END IF

      IF( OnTheFlyBC ) THEN
        IF( Trgt == Constraint ) n2 = n2 + 1
      ELSE
        IF ( Model % BCs(Trgt) % Tag == Constraint ) n2 = n2 + 1
      END IF
    END DO

    IF( CheckForHalo ) THEN
      CALL Info('CreateInterfaceMeshes','Number of halo elements eliminated: '&
          //TRIM(I2S(HaloCount)),Level=12)
    END IF

    IF ( n1 <= 0 .OR. n2 <= 0 ) THEN
      ! This is too conservative in parallel
      ! CALL Warn('CreateInterfaceMeshes','There are no active boundaries!')
      Success = .FALSE.
      RETURN
    END IF


!   Initialize mesh structures for boundaries, this
!   is for getting the mesh projector:
!   ------------------------------------------------
    BMesh1 % Parent => Mesh
    BMesh2 % Parent => Mesh

    CALL AllocateVector( BMesh1 % Elements,n1 )
    CALL AllocateVector( BMesh2 % Elements,n2 )
    CALL AllocateVector( Perm1, Mesh % NumberOfNodes )
    CALL AllocateVector( Perm2, Mesh % NumberOfNodes )

 
!   Fill in the mesh element structures with the
!   boundary elements:
!   ---------------------------------------------
    n1 = 0
    n2 = 0
    Perm1 = 0
    Perm2 = 0
    BMesh1 % MaxElementNodes = 0
    BMesh2 % MaxElementNodes = 0


    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      
      Constraint = Element % BoundaryInfo % Constraint
      
      ThisActive = ( Model % BCs(This) % Tag == Constraint ) 
      IF( ThisActive .AND. CheckForHalo ) THEN
        IF( NarrowHalo ) THEN
          IF( .NOT. ANY(ActiveNode(Element % NodeIndexes) ) ) THEN
            ThisActive = .FALSE.
          END IF
        ELSE IF( NoHalo ) THEN
          ThisActive = .FALSE.
          Left => Element % BoundaryInfo % Left 
          IF( ASSOCIATED( Left ) ) THEN
            ThisActive = ( Left % PartIndex == ParEnv % MyPe )
          END IF
          Right => Element % BoundaryInfo % Right
          IF( ASSOCIATED( Right ) ) THEN
            ThisActive = ThisActive .OR. &
                ( Right % PartIndex == ParEnv % MyPe ) 
          END IF
        END IF
      END IF

      IF( OnTheFlyBC ) THEN
        TargetActive = ( Trgt == Constraint )
      ELSE
        TargetActive = ( Model % BCs(Trgt) % Tag == Constraint ) 
      END IF

      IF(.NOT. (ThisActive .OR. TargetActive ) ) CYCLE
      
      ! Set the pointers accordingly so we need to code the complex stuff
      ! only once.
      IF ( ThisActive ) THEN
        n1 = n1 + 1
        ind = n1
        PMesh => BMesh1
        PPerm => Perm1
      ELSE
        n2 = n2 + 1
        ind = n2
        PMesh => BMesh2
        PPerm => Perm2
      END IF

      n = Element % TYPE % NumberOfNodes             
      PMesh % MaxElementNodes = MAX( PMesh % MaxElementNodes, n )
      PMesh % Elements(ind) = Element

      CALL AllocateVector(PMesh % Elements(ind) % NodeIndexes,n )

      IF( Mesh % NumberOfFaces == 0 .OR. Mesh % NumberOfEdges == 0 ) THEN
        PMesh % Elements(ind) % NodeIndexes(1:n) = Element % NodeIndexes(1:n)
      ELSE
        ! If we have edge dofs we want the face element be associated with the 
        ! face list since that only has properly defined edge indexes.
        Parent => Element % BoundaryInfo % Left
        IF(.NOT. ASSOCIATED( Parent ) ) THEN
          Parent => Element % BoundaryInfo % Right
        END IF

        q => Find_Face(Parent,Element)
        PMesh % Elements(ind) % NodeIndexes(1:n) = q % NodeIndexes(1:n)

        ! set the elementindex to be faceindex as it may be needed
        ! for the edge elements.
        PMesh % Elements(ind) % ElementIndex = q % ElementIndex

        IF(ASSOCIATED(q % Pdefs)) THEN
          ALLOCATE(Pmesh % Elements(ind) % Pdefs)
          PMesh % Elements(ind) % PDefs = q % Pdefs
        END IF

        ! Set also the owner partition
!       PMesh % Elements(ind) % PartIndex = q % PartIndex

        en = q % TYPE % NumberOfEdges
        ALLOCATE(PMesh % Elements(ind) % EdgeIndexes(en))
        Pmesh % Elements(ind) % EdgeIndexes(1:en) = q % EdgeIndexes(1:en)
      END IF

      PPerm( Element % NodeIndexes(1:n) ) = 1

    END DO
  
!   Fill in the mesh node structures with the
!   boundary nodes:
!   -----------------------------------------
    BMesh1 % NumberOfBulkElements = n1
    BMesh2 % NumberOfBulkElements = n2

    BMesh2 % NumberOfNodes = COUNT(Perm2 > 0)
    BMesh1 % NumberOfNodes = COUNT(Perm1 > 0)

    ! As there were some active boundary elements this condition should 
    ! really never be possible   
    IF (BMesh1 % NumberOfNodes==0 .OR. BMesh2 % NumberOfNOdes==0) THEN
      CALL Fatal('CreateInterfaceMeshes','No active nodes on periodic boundary!')
    END IF

    WRITE(Message,'(A,I0,A,I0)') 'Number of periodic nodes: ',&
        BMesh1 % NumberOfNodes, ', ',BMesh2 % NumberOfNOdes
    CALL Info('CreateInterfaceMeshes',Message,Level=9)    
    
    ALLOCATE( BMesh1 % Nodes )
    CALL AllocateVector( BMesh1 % Nodes % x, BMesh1 % NumberOfNodes ) 
    CALL AllocateVector( BMesh1 % Nodes % y, BMesh1 % NumberOfNodes ) 
    CALL AllocateVector( BMesh1 % Nodes % z, BMesh1 % NumberOfNodes )
    
    ALLOCATE( BMesh2 % Nodes )
    CALL AllocateVector( BMesh2 % Nodes % x, BMesh2 % NumberOfNodes ) 
    CALL AllocateVector( BMesh2 % Nodes % y, BMesh2 % NumberOfNodes ) 
    CALL AllocateVector( BMesh2 % Nodes % z, BMesh2 % NumberOfNodes )
    
    CALL AllocateVector( Bmesh1 % InvPerm, BMesh1 % NumberOfNodes )
    CALL AllocateVector( Bmesh2 % InvPerm, BMesh2 % NumberOfNodes )

    ! Now, create the master and target meshes that only include the active elements
    !---------------------------------------------------------------------------
    k1 = 0; k2 = 0
    DO i=1,Mesh % NumberOfNodes

      IF ( Perm1(i) > 0 ) THEN
        k1 = k1 + 1
        Perm1(i) = k1
        BMesh1 % InvPerm(k1) = i

        BMesh1 % Nodes % x(k1) = Mesh % Nodes % x(i)
        BMesh1 % Nodes % y(k1) = Mesh % Nodes % y(i)
        BMesh1 % Nodes % z(k1) = Mesh % Nodes % z(i)
      END IF
      
      IF ( Perm2(i) > 0 ) THEN
        k2 = k2 + 1
        Perm2(i) = k2
        BMesh2 % InvPerm(k2) = i
        
        BMesh2 % Nodes % x(k2)= Mesh % Nodes % x(i)
        BMesh2 % Nodes % y(k2)= Mesh % Nodes % y(i)
        BMesh2 % Nodes % z(k2)= Mesh % Nodes % z(i)
      END IF
    END DO

!   Finally, Renumber the element node pointers to use
!   only boundary nodes:
!   ---------------------------------------------------
    DO i=1,n1
      BMesh1 % Elements(i) % NodeIndexes = Perm1(BMesh1 % Elements(i) % NodeIndexes)
    END DO

    DO i=1,n2
      BMesh2 % Elements(i) % NodeIndexes = Perm2(BMesh2 % Elements(i) % NodeIndexes)
    END DO
    DEALLOCATE( Perm1, Perm2 )

    IF( CheckForHalo ) DEALLOCATE( ActiveNode ) 

    Success = .TRUE.

  END SUBROUTINE CreateInterfaceMeshes
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !> Given two meshes that should occupy the same domain in space 
  !> use rotation, scaling and translation to achive this goal.
  !---------------------------------------------------------------------------
  SUBROUTINE OverlayIntefaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    LOGICAL :: GotIt, GotRotate
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3),x2r_min(3),x2r_max(3)
    REAL(KIND=dp) :: x(4), RotMatrix(4,4),TrsMatrix(4,4),SclMatrix(4,4), &
           TrfMatrix(4,4),Identity(4,4),Angles(3),Alpha,scl(3),s1,s2
    REAL(KIND=dp), POINTER :: PArray(:,:)
    INTEGER :: i,j,k

    ! First, check the bounding boxes
    !---------------------------------------------------------------------------
    x1_min(1) = MINVAL( BMesh1 % Nodes % x )
    x1_min(2) = MINVAL( BMesh1 % Nodes % y )
    x1_min(3) = MINVAL( BMesh1 % Nodes % z )
    
    x1_max(1) = MAXVAL( BMesh1 % Nodes % x )
    x1_max(2) = MAXVAL( BMesh1 % Nodes % y )
    x1_max(3) = MAXVAL( BMesh1 % Nodes % z )

    WRITE(Message,'(A,3ES15.6)') 'Minimum values for this periodic BC:  ',x1_min
    CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
    WRITE(Message,'(A,3ES15.6)') 'Maximum values for this periodic BC:  ',x1_max
    CALL Info('OverlayInterfaceMeshes',Message,Level=8)    

    x2_min(1) = MINVAL( BMesh2 % Nodes % x )
    x2_min(2) = MINVAL( BMesh2 % Nodes % y )
    x2_min(3) = MINVAL( BMesh2 % Nodes % z )
    
    x2_max(1) = MAXVAL( BMesh2 % Nodes % x )
    x2_max(2) = MAXVAL( BMesh2 % Nodes % y )
    x2_max(3) = MAXVAL( BMesh2 % Nodes % z )
    
    WRITE(Message,'(A,3ES15.6)') 'Minimum values for target periodic BC:',x2_min
    CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
    WRITE(Message,'(A,3ES15.6)') 'Maximum values for target periodic BC:',x2_max
    CALL Info('OverlayInterfaceMeshes',Message,Level=8)    

!    If whole transformation matrix given, it will be used directly
!    --------------------------------------------------------------
    Parray => ListGetConstRealArray( BParams,'Periodic BC Matrix', Gotit )
    IF ( GotIt ) THEN
      DO i=1,SIZE(Parray,1)
        DO j=1,SIZE(Parray,2)
          TrfMatrix(i,j) = Parray(j,i)
        END DO
      END DO
    ELSE    
      ! Otherwise check for rotation, scaling and translation
      !------------------------------------------------------

      ! Initialize the mapping matrices
      Identity = 0.0d0
      DO i=1,4
        Identity(i,i) = 1.0d0
      END DO      
      TrsMatrix = Identity
      RotMatrix = Identity
      SclMatrix = Identity
      
      !   Rotations:
      !   These are called first since they are not accounted for in the 
      !   automatic scaling and translation.
      !   ---------------------------------------------------------------      
      Angles = 0.0_dp
      Parray => ListGetConstRealArray( BParams,'Periodic BC Rotate', GotRotate )
      IF( GotRotate ) THEN
        Angles(1:3) = Parray(1:3,1)   
      ELSE
        IF( ListGetLogical( BParams,'Periodic BC Rotate Automatic', GotIt) ) THEN
          CALL CheckInterfaceMeshAngle( BMesh1, BMesh2, Angles, GotRotate ) 
        END IF
      END IF

      IF ( GotRotate ) THEN
        WRITE(Message,'(A,3ES15.6)') 'Rotating target with: ',Angles
        CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
        
        DO i=1,3
          Alpha = Angles(i) * PI / 180
          IF( ABS(Alpha) < TINY(Alpha) ) CYCLE 
          TrfMatrix = Identity
          
          SELECT CASE(i)
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
        
        DO i = 1, BMesh2 % NumberOfNodes          
          x(1) = BMesh2 % Nodes % x(i)
          x(2) = BMesh2 % Nodes % y(i)
          x(3) = BMesh2 % Nodes % z(i)
          
          x(4) = 1.0_dp
          x = MATMUL( RotMatrix, x )
          
          BMesh2 % Nodes % x(i) = x(1)
          BMesh2 % Nodes % y(i) = x(2)
          BMesh2 % Nodes % z(i) = x(3)
        END DO
        
        x2r_min(1) = MINVAL( BMesh2 % Nodes % x )
        x2r_min(2) = MINVAL( BMesh2 % Nodes % y )
        x2r_min(3) = MINVAL( BMesh2 % Nodes % z )
        
        x2r_max(1) = MAXVAL( BMesh2 % Nodes % x )
        x2r_max(2) = MAXVAL( BMesh2 % Nodes % y )
        x2r_max(3) = MAXVAL( BMesh2 % Nodes % z )
        
        WRITE(Message,'(A,3ES15.6)') 'Minimum values for rotated target:',x2r_min
        CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
        
        WRITE(Message,'(A,3ES15.6)') 'Maximum values for rotated target:',x2r_max
        CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
      ELSE
        x2r_min = x2_min
        x2r_max = x2_max
      END IF
   
!   Scaling:
!   This is either given or enforced by requiring bounding boxes to be of the same size 
!   -----------------------------------------------------------------------------------
      Parray => ListGetConstRealArray( BParams,'Periodic BC Scale', Gotit )      
      IF ( GotIt ) THEN
        DO i=1,SIZE(Parray,1)
          SclMatrix(i,i) = Parray(i,1)
        END DO
      ELSE
        ! Define scaling from the bounding boxes
        ! This assumes isotropic scaling since component-wise scaling 
        ! was prone to errors.
        !------------------------------------------------------
        s1 = SUM( ( x1_max(1:3) - x1_min(1:3) ) ** 2 )
        s2 = SUM( ( x2r_max(1:3) - x2r_min(1:3) ) ** 2 )
        IF( s2 > EPSILON( s2 ) ) THEN
          scl(1:3)  = SQRT( s1 / s2 )
        ELSE
          scl(1:3) = 1.0_dp
        END IF
        
        WRITE(Message,'(A,3ES15.6)') 'Scaling with: ',scl(1:3)
        CALL Info('OverlayInterfaceMeshes',Message)
        DO i=1,3 
          SclMatrix(i,i) = scl(i)        
        END DO
      END IF
      
!   Translations:
!   And finally define translations
!   -------------
      Parray => ListGetConstRealArray( BParams,'Periodic BC Translate', Gotit )
      IF ( gotit ) THEN
        DO i=1,SIZE(Parray,1)
          TrsMatrix(4,i) = Parray(i,1)
        END DO
      ELSE
        ! Define translations so that the lower left corner is the same
        !-------------------------------------------------------------
        DO i=1,3
          TrsMatrix(4,i) = x1_min(i) - SclMatrix(i,i) * x2r_min(i)
        END DO
      END IF
      WRITE(Message,'(A,3ES15.6)') 'Translation: ',TrsMatrix(4,1:3)
      CALL Info('OverlayInterfaceMeshes',Message)
      TrfMatrix = MATMUL( SclMatrix, TrsMatrix )
    END IF

!    Now transform the coordinates:
!    ------------------------------
    DO i=1,BMesh2 % NumberOfNodes
      x(1) = BMesh2 % Nodes % x(i)
      x(2) = BMesh2 % Nodes % y(i)
      x(3) = BMesh2 % Nodes % z(i)
      x(4) = 1.0d0
      x = MATMUL( x, TrfMatrix ) 
      BMesh2 % Nodes % x(i) = x(1) / x(4)
      BMesh2 % Nodes % y(i) = x(2) / x(4) 
      BMesh2 % Nodes % z(i) = x(3) / x(4)
    END DO

    IF(.FALSE.) THEN
      x2r_min(1) = MINVAL( BMesh2 % Nodes % x )
      x2r_min(2) = MINVAL( BMesh2 % Nodes % y )
      x2r_min(3) = MINVAL( BMesh2 % Nodes % z )
      
      x2r_max(1) = MAXVAL( BMesh2 % Nodes % x )
      x2r_max(2) = MAXVAL( BMesh2 % Nodes % y )
      x2r_max(3) = MAXVAL( BMesh2 % Nodes % z )
      
      WRITE(Message,'(A,3ES15.6)') 'Minimum values for transformed target:',x2r_min
      CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
      
      WRITE(Message,'(A,3ES15.6)') 'Maximum values for transformed target:',x2r_max
      CALL Info('OverlayInterfaceMeshes',Message,Level=8)    
    END IF

  END SUBROUTINE OverlayIntefaceMeshes
  !---------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming rotating boundaries make 
  !> a coordinate transformation to each node of the slave boundary (BMesh1) so that 
  !> they hit the master boundary (BMesh2). In case of anti-periodic projector 
  !> mark the nodes that need an odd number of periods.
  !---------------------------------------------------------------------------
  SUBROUTINE PreRotationalProjector(BMesh1, BMesh2, MirrorNode )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    !--------------------------------------------------------------------------
    LOGICAL :: AntiPeriodic
    REAL(KIND=dp) :: F2min,F2max,dFii2,Fii
    INTEGER :: i, Nfii, SectorMax
    INTEGER, ALLOCATABLE :: SectorCount(:)

    AntiPeriodic = ALLOCATED( MirrorNode )
    IF( AntiPeriodic ) MirrorNode = .FALSE.

    F2Min =  MINVAL( BMesh2 % Nodes % x )
    F2Max =  MAXVAL( BMesh2 % Nodes % x )
    dFii2 = F2Max - F2Min
    SectorMax = CEILING( 360.0 / dFii2 ) 

    WRITE( Message,'(A,I0)') 'Maximum number of sectors: ',SectorMax
    CALL Info('PreRotationalProjector',Message,Level=8)

    ALLOCATE( SectorCount(-SectorMax:SectorMax))
    SectorCount = 0

    DO i = 1, BMesh1 % NumberOfNodes
      Fii = BMesh1 % Nodes % x(i)      
      Nfii = FLOOR( (Fii-F2min) / dFii2 )
      BMesh1 % Nodes % x(i) = BMesh1 % Nodes % x(i) - Nfii * dFii2
      SectorCount(Nfii) = SectorCount(Nfii) + 1     
      IF( AntiPeriodic ) THEN
        IF( MODULO(Nfii,2) /= 0 ) THEN
          MirrorNode(i) = .TRUE.
        END IF
      END IF
    END DO

    IF( SectorCount(0) < BMesh1 % NumberOfNodes ) THEN
      CALL Info('PreRotationalProjector','Number of nodes by sectors',Level=8)
      DO i=-SectorMax,SectorMax
        IF( SectorCount(i) > 0 ) THEN
          WRITE( Message,'(A,I0,A,I0)') 'Sector:',i,'   Nodes:',SectorCount(i)
          CALL Info('MatchInterfaceNodes',Message,Level=8)
        END IF
      END DO
      IF( AntiPeriodic ) THEN
        WRITE( Message,'(A,I0)') 'Number of mirror nodes:',COUNT(MirrorNode)
        CALL Info('PreRotationalProjector',Message,Level=8)
      END IF
    ELSE
      CALL Info('PreRotationalProjector','No nodes needed mapping')
    END IF

  END SUBROUTINE PreRotationalProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Postprocess projector so that it changes the sign of the anti-periodic
!> entries as assigns by the MirrorNode flag.
!------------------------------------------------------------------------------
  SUBROUTINE PostRotationalProjector( Proj, MirrorNode )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: Proj                 !< Projection matrix
    LOGICAL, ALLOCATABLE :: MirrorNode(:)  !< Is the node a mirror node or not
!--------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)            
    REAL(KIND=dp), POINTER :: Values(:)    
    INTEGER :: i,j,n
!------------------------------------------------------------------------------

    IF( .NOT. ALLOCATED( MirrorNode ) ) RETURN
    IF( COUNT( MirrorNode ) == 0 ) RETURN

    n = Proj % NumberOfRows
    Rows => Proj % Rows
    Cols => Proj % Cols
    Values => Proj % Values

    DO i=1,n
      IF( MirrorNode(i) ) THEN
        DO j = Rows(i),Rows(i+1)-1
          Values(j) = -Values(j)
        END DO
      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE PostRotationalProjector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION Find_Face(Parent,Element) RESULT(ptr)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Ptr
    TYPE(Element_t) :: Parent, Element

    INTEGER :: i,j,k,n

    Ptr => NULL()
    DO i=1,Parent % TYPE % NumberOfFaces
      Ptr => CurrentModel % Mesh % Faces(Parent % FaceIndexes(i))
      n=0
      DO j=1,Ptr % TYPE % NumberOfNodes
        DO k=1,Element % TYPE % NumberOfNodes
          IF (Ptr % NodeIndexes(j) == Element % NodeIndexes(k)) n=n+1
        END DO
      END DO
      IF (n==Ptr % TYPE % NumberOfNodes) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION Find_Face
!------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Create a projector for mapping between interfaces using the Galerkin method
  !> A temporal mesh structure with a node for each Gaussian integration point is 
  !> created. Then this projector matrix is transferred to a projector on the nodal
  !> coordinates.   
  !---------------------------------------------------------------------------
   FUNCTION NodalProjector(BMesh2, BMesh1, &
       UseQuadrantTree, Repeating, AntiRepeating ) &
      RESULT ( Projector )
  !---------------------------------------------------------------------------
    USE Lists

    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    LOGICAL :: UseQuadrantTree, Repeating, AntiRepeating
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    INTEGER :: i,j,k,n
    INTEGER, POINTER :: Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:)

    BMesh1 % Parent => NULL()
    BMesh2 % Parent => NULL()

    InvPerm1 => BMesh1 % InvPerm
    InvPerm2 => BMesh2 % InvPerm

    ! Set the nodes of Mesh1 to be in the interval defined by Mesh2
    !-----------------------------------------------------------------
    IF( Repeating ) THEN
      IF( AntiRepeating ) THEN
        ALLOCATE( MirrorNode( BMesh1 % NumberOfNodes ) )
        MirrorNode = .FALSE.
      END IF
      CALL PreRotationalProjector(BMesh1, BMesh2, MirrorNode )
    END IF

    ! Create the projector using nodal points 
    ! This corresponds to numerical integration of the collocation method.
    !-----------------------------------------------------------------
    Projector => MeshProjector( BMesh2, BMesh1, UseQuadrantTree )    
    Projector % ProjectorType = PROJECTOR_TYPE_NODAL

    Values => Projector % Values
    Cols => Projector % Cols
    Rows => Projector % Rows

    ! One needs to change the sign of the projector for the mirror nodes
    !-----------------------------------------------------------------------------
    IF( AntiRepeating ) THEN
      CALL PostRotationalProjector( Projector, MirrorNode )
      DEALLOCATE( MirrorNode ) 
    END IF

    ! Now return from the indexes of the interface mesh system to the 
    ! original mesh system.
    !-----------------------------------------------------------------
    n = SIZE( InvPerm1 ) 
    ALLOCATE( Projector % InvPerm(n) )
    Projector % InvPerm = InvPerm1

    DO i=1,Projector % NumberOfRows
       DO j = Rows(i), Rows(i+1)-1
         k = Cols(j)    
         IF ( k > 0 ) Cols(j) = InvPerm2(k)
       END DO
    END DO

  END FUNCTION NodalProjector
!------------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !> Create a nodal projector related to discontinous interface.
  !---------------------------------------------------------------------------
   FUNCTION NodalProjectorDiscont( Mesh, bc ) RESULT ( Projector )
  !---------------------------------------------------------------------------
    USE Lists

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: bc
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    TYPE(Model_t), POINTER :: Model
    INTEGER, POINTER :: NodePerm(:)
    INTEGER :: i,j,n,m
    INTEGER, POINTER :: Rows(:),Cols(:), InvPerm(:)
    REAL(KIND=dp), POINTER :: Values(:)
    LOGICAL :: Found

    CALL Info('NodalProjectorDiscont','Creating nodal projector for discontinuous boundary',Level=7)

    Projector => Null()
    IF( .NOT. Mesh % DisContMesh ) THEN
      CALL Warn('NodalProjectorDiscont','Discontinuous mesh not created?')
      RETURN
    END IF

    Model => CurrentModel
    j = 0
    DO i=1,Model % NumberOfBCs
      IF( ListGetLogical(Model % BCs(i) % Values,'Discontinuous Boundary',Found) ) THEN
        j = j + 1
      END IF
    END DO
    ! This is a temporal limitations
    IF( j > 1 ) THEN
      CALL Warn('NodalProjectorDiscont','One BC (not '&
          //TRIM(I2S(j))//') only for discontinuous boundary!')
    END IF


    NodePerm => Mesh % DisContPerm
    n = SIZE( NodePerm ) 
    m = COUNT( NodePerm > 0 ) 

    Projector => AllocateMatrix()
    Projector % ProjectorType = PROJECTOR_TYPE_NODAL
    Projector % ProjectorBC = bc

    ALLOCATE( Projector % Cols(m) )
    ALLOCATE( Projector % Values(m) )
    ALLOCATE( Projector % Rows(m+1) )
    ALLOCATE( Projector % InvPerm(m) )

    Cols => Projector % Cols
    Values => Projector % Values
    Rows => Projector % Rows
    InvPerm => Projector % InvPerm
    Projector % NumberOfRows = m

    Values = 1.0_dp
    DO i=1,m+1
      Rows(i) = i
    END DO

    DO i=1,n
      j = NodePerm(i)
      IF( j == 0 ) CYCLE
      Cols(j) = n + j
      InvPerm(j) = i
    END DO

  END FUNCTION NodalProjectorDiscont
!------------------------------------------------------------------------------


 
  !---------------------------------------------------------------------------
  !> Create a projector for mixes nodal / edge problems assuming constant level
  !> in the 2nd direction. This kind of projector is suitable for 2D meshes where
  !> the mortar line is effectively 1D, or to 3D cases that have been created by
  !> extrusion. 
  !---------------------------------------------------------------------------
  FUNCTION LevelProjector( BMesh1, BMesh2, Repeating, AntiRepeating, &
      FullCircle, Radius, DoNodes, DoEdges, &
      NodeScale, EdgeScale, BC ) &
      RESULT ( Projector )
    !---------------------------------------------------------------------------
    USE Lists
    USE Messages
    USE Types
    USE GeneralUtils
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2, Mesh
    LOGICAL :: DoNodes, DoEdges
    LOGICAL :: Repeating, AntiRepeating, FullCircle, NotAllQuads, NotAllQuads2
    REAL(KIND=dp) :: Radius, NodeScale, EdgeScale
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Matrix_t), POINTER :: Projector    
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
    LOGICAL ::  StrongNodes, StrongLevelEdges, StrongExtrudedEdges, StrongSkewEdges
    LOGICAL :: Found, Parallel, SelfProject, EliminateUnneeded, SomethingUndone, &
        EdgeBasis, PiolaVersion, GenericIntegrator, Rotational, Cylindrical, IntGalerkin, &
        CreateDual, HaveMaxDistance
    REAL(KIND=dp) :: XmaxAll, XminAll, YminAll, YmaxAll, Xrange, Yrange, &
        RelTolX, RelTolY, XTol, YTol, RadTol, MaxSkew1, MaxSkew2, SkewTol, &
        ArcCoeff, EdgeCoeff, NodeCoeff, MaxDistance
    INTEGER :: NoNodes1, NoNodes2, MeshDim
    INTEGER :: i,j,k,n,m,Nrange,Nrange2, nrow
    INTEGER, ALLOCATABLE :: EdgePerm(:),NodePerm(:),DualNodePerm(:)
    INTEGER :: EdgeRow0, FaceRow0, EdgeCol0, FaceCol0, ProjectorRows
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Cond(:)
    TYPE(Matrix_t), POINTER :: DualProjector    
    LOGICAL :: DualMaster, DualSlave, DualLCoeff, BiorthogonalBasis

    CALL Info('LevelProjector','Creating projector for a levelized mesh',Level=7)

    IF(.NOT. (DoEdges .OR. DoNodes ) ) THEN
      CALL Warn('LevelProjector','Nothing to do, no nonodes, no edges!')
      RETURN
    END IF

    EdgeCoeff = ListGetConstReal( BC,'Projector Edge Multiplier',Found )
    IF( .NOT. Found ) EdgeCoeff = ListGetConstReal( CurrentModel % Simulation,&
        'Projector Edge Multiplier',Found )
    IF( .NOT. Found ) EdgeCoeff = 1.0_dp

    NodeCoeff = ListGetConstReal( BC,'Projector Node Multiplier',Found )
    IF( .NOT. Found ) NodeCoeff = ListGetConstReal( CurrentModel % Simulation,&
        'Projector Node Multiplier',Found )
    IF( .NOT. Found ) NodeCoeff = 1.0_dp

    Rotational = ListGetLogical( BC,'Rotational Projector',Found ) .OR. &
        ListGetLogical( BC,'Anti Rotational Projector',Found )
    Cylindrical = ListGetLogical( BC,'Cylindrical Projector',Found ) 
    
    MaxDistance = ListGetCReal( BC,'Projector Max Distance', HaveMaxDistance) 
    IF(.NOT. HaveMaxDistance ) THEN
      MaxDistance = ListGetCReal( CurrentModel % Solver % Values,&
          'Projector Max Distance', HaveMaxDistance)       
    END IF


    Parallel = ( ParEnv % PEs > 1 )
    Mesh => CurrentModel % Mesh
    BMesh1 % Parent => NULL()
    BMesh2 % Parent => NULL()

    ! Create a projector in style P=I-Q, or rather just P=Q. 
    SelfProject = .TRUE.
    
    ! Range is needed to define tolerances, and to map the angle in case 
    ! the master mesh is treated as a repeating structure. 
    XMaxAll = MAXVAL(BMesh2 % Nodes % x)
    XMinAll = MINVAL(BMesh2 % Nodes % x)
    XRange = XMaxAll - XMinAll

    YMaxAll = MAXVAL(BMesh2 % Nodes % y)
    YMinAll = MINVAL(BMesh2 % Nodes % y)
    YRange = YMaxAll - YMinAll

    ! Fix here the relative tolerance used to define the search tolerance
    RelTolY = 1.0d-4
    ! In the case of infinite target we can have tighter criteria
    IF( FullCircle .OR. Repeating ) THEN
      RelTolX = 1.0d-6
    ELSE
      RelTolX = RelTolY
    END IF
    YTol = RelTolY * YRange
    XTol = RelTolX * XRange

    ! Determine the coefficient that turns possible angles into units of
    ! ach-lenth. If this is not rotational then there are no angles. 
    IF( Rotational .OR. Cylindrical ) THEN
      ArcCoeff = (2*PI*Radius)/360.0           
    ELSE
      ArcCoeff = 1.0_dp
    END IF

    IntGalerkin = ListGetLogical( BC, 'Galerkin Projector', Found )
    MeshDIm = Mesh % MeshDim

    ! Generic integrator does not make any assumptions on the way the mesh 
    ! is constructured. Otherwise constant strides in y-direction is assumed. 
    ! For weak strategy always use the generic integrator. 
    GenericIntegrator = ListGetLogical( BC,'Level Projector Generic',Found ) 
    IF(.NOT. Found ) GenericIntegrator = IntGalerkin

    ! There is no strong strategy for skewed edges currently
    StrongSkewEdges = .FALSE.
    ! Maximum skew in degrees before treating edges as skewed
    SkewTol = 0.1  
    

    IF( MeshDim == 2 ) THEN
      ! In 2D these always yield since current only the strong method 
      ! of nodes has been implemented in 2D 
      CALL Info('LevelProjector','Initial mesh is 2D, using 1D projectors!',Level=10) 
      StrongNodes = ListGetLogical( BC,'Level Projector Nodes Strong',Found ) 
      IF(.NOT. Found) StrongNodes = ListGetLogical( BC,'Level Projector Strong',Found ) 
      IF(.NOT. Found) StrongNodes = .NOT. IntGalerkin
    ELSE ! 3D 
      ! It is assumed that that the target mesh is always un-skewed 
      ! Make a test here to be able to skip it later. No test is needed
      ! if the generic integrator is enforced. 
      IF(.NOT. GenericIntegrator ) THEN
        MaxSkew1 = CheckMeshSkew( BMesh1, NotAllQuads )
        IF( NotAllQuads ) THEN
          CALL Info('LevelProjector','This mesh has also triangles',Level=8)
        END IF
        WRITE( Message,'(A,ES12.3)') 'Maximum skew in this mesh: ',MaxSkew1
        CALL Info('LevelProjector',Message,Level=8)
        
        MaxSkew2 = CheckMeshSkew( BMesh2, NotAllQuads2 )
        IF( NotAllQuads2 ) THEN
          CALL Info('LevelProjector','Target mesh has also triangles',Level=8)
        END IF
        WRITE( Message,'(A,ES12.3)') 'Maximum skew in target mesh: ',MaxSkew2
        CALL Info('LevelProjector',Message,Level=8)
        
        IF( NotAllQuads .OR. NotAllQuads2 .OR. MaxSkew2 > SkewTol ) THEN
          IF( MaxSkew2 > MaxSkew1 .AND. MaxSkew1 < SkewTol ) THEN
            CALL Warn('LevelProjector','You could try switching the master and target BC!')
          END IF
          CALL Warn('LevelProjector','Target mesh has too much skew, using generic integrator!')
          GenericIntegrator = .TRUE. 
        END IF
      END IF
        
      ! The projectors for nodes and edges can be created either in a strong way 
      ! or weak way in the special case that the nodes are located in extruded layers. 
      ! The strong way results to a sparse projector. For constant 
      ! levels it can be quite optimal, except for the edges with a skew. 
      ! If strong projector is used for all edges then "StrideProjector" should 
      ! be recovered.
      IF( DoNodes ) THEN
        StrongNodes = .NOT. IntGalerkin
        IF( GenericIntegrator ) THEN
          StrongNodes = .FALSE.
        ELSE IF( ListGetLogical( BC,'Level Projector Strong',Found ) ) THEN
          StrongNodes = .TRUE.
        END IF
        ! The nodes could be treated with strong projector even though the edges are integrated
        ! with a weak projector. 
        IF( ListGetLogical( BC,'Level Projector Nodes Strong',Found ) ) StrongNodes = .TRUE.
      END IF
      
      IF( DoEdges ) THEN
        StrongLevelEdges = .NOT. IntGalerkin
        StrongExtrudedEdges = .NOT. IntGalerkin
        StrongSkewEdges = .FALSE.
        
        IF( GenericIntegrator ) THEN
          CALL Info('LevelProjector','Using generic weak projector for all edge dofs!')
          StrongLevelEdges = .FALSE.
          StrongExtrudedEdges = .FALSE.
        ELSE
          IF( ListGetLogical( BC,'Level Projector Strong',Found ) .OR. &
              ListGetLogical( BC,'Level Projector Edges Strong',Found ) ) THEN
            StrongLevelEdges = .TRUE.
            StrongExtrudedEdges = .TRUE.
          END IF
          IF( ListGetLogical( BC,'Level Projector Plane Edges Strong',&
              Found ) ) StrongLevelEdges = .TRUE.
          IF( ListGetLogical( BC,'Level Projector Extruded Edges Strong',&
              Found ) ) StrongExtrudedEdges = .TRUE.
          IF( ListGetLogical( BC,'Level Projector Skew Edges Strong',&
              Found ) ) StrongSkewEdges = .TRUE.
        END IF
      END IF
    END IF


    ! If the number of periods is enforced use that instead since
    ! the Xrange periodicity might not be correct if the mesh has skew.
    IF( Rotational ) THEN
      IF( FullCircle ) THEN
        Xrange = 360.0_dp
      ELSE 
        i = ListGetInteger( BC,'Rotational Projector Periods',Found,minv=1 ) 
        IF( GenericIntegrator .AND. .NOT. Found ) THEN
          CALL Fatal('LevelProjector',&
              'Generic integrator requires > Rotational Projector Periods <')
        END IF
        Xrange = 360.0_dp / i
      END IF
    END IF

    ! This is the tolerance used to define constant direction in radians
    ! For consistancy it should not be sloppier than the SkewTol
    ! but it could be equally sloppy as below. 
    RadTol = PI * SkewTol / 180.0_dp

    ! Given the inverse permutation compute the initial number of
    ! nodes in both cases. 
    NoNodes1 = BMesh1 % NumberOfNodes
    NoNodes2 = BMesh2 % NumberOfNodes

    InvPerm1 => BMesh1 % InvPerm
    InvPerm2 => BMesh2 % InvPerm

    ! Create a list matrix that allows for unspecified entries in the matrix 
    ! structure to be introduced.
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN

    CreateDual = ListGetLogical( BC,'Create Dual Projector',Found ) 
    IF( CreateDual ) THEN
      DualProjector => AllocateMatrix()
      DualProjector % FORMAT = MATRIX_LIST
      DualProjector % ProjectorType = PROJECTOR_TYPE_GALERKIN
      Projector % EMatrix => DualProjector
    END IF

    ! Check whether biorthogonal basis for projectors requested:
    ! ----------------------------------------------------------
    BiOrthogonalBasis = ListGetLogical( BC, 'Use Biorthogonal Basis', Found)

    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      BiOrthogonalBasis = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( BiOrthogonalBasis ) THEN
        CALL Info('LevelProjector',&
            'Enforcing > Use Biorthogonal Basis < to True to enable elimination',Level=8)
        CALL ListAddLogical( BC, 'Use Biorthogonal Basis',.TRUE. )
      END IF
    END IF

    IF (BiOrthogonalBasis) THEN
      IF( DoEdges ) THEN
        CALL Warn('LevelProjector','Biorthogonal basis cannot be combined with edge elements!')
      END IF

      DualSlave  = ListGetLogical(BC, 'Biorthogonal Dual Slave', Found)
      IF(.NOT.Found) DualSlave  = .TRUE.

      DualMaster = ListGetLogical(BC, 'Biorthogonal Dual Master', Found)
      IF(.NOT.Found) DualMaster = .TRUE.

      DualLCoeff = ListGetLogical(BC, 'Biorthogonal Dual Lagrange Coefficients', Found)
      IF(.NOT.Found) DualLCoeff = .FALSE.

      IF(DualLCoeff) THEN
        DualSlave  = .FALSE.
        DualMaster = .FALSE.
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.FALSE.)
      ELSE
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.TRUE.)
      END IF

      Projector % Child => AllocateMatrix()
      Projector % Child % Format = MATRIX_LIST
      CALL Info('LevelProjector','Using biorthogonal basis, as requested',Level=8)      
    END IF


    PiolaVersion = ListGetLogical( CurrentModel % Solver % Values, &
        'Use Piola Transform', Found)


    ! At the 1st stage determine the maximum size of the projector
    ! If the strong projector is used then the numbering is done as we go
    ! this way we can eliminate unneeded rows. 
    ! For the weak projector there is no need to eliminate rows. 
    IF( DoNodes ) THEN      
      ALLOCATE( NodePerm( Mesh % NumberOfNodes ) )
      NodePerm = 0

      ! in parallel only consider nodes that truly are part of this partition
      DO i=1,BMesh1 % NumberOfBulkElements
        Element => BMesh1 % Elements(i)        
        IF( Parallel ) THEN
          IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE          
        END IF        
        NodePerm( InvPerm1( Element % NodeIndexes ) ) = 1
      END DO

      n = SUM( NodePerm )
      CALL Info('LevelProjector','Initial number of periodic nodes '//TRIM(I2S(n))//&
          ' out of '//TRIM(I2S(BMesh1 % NumberOfNodes ) ), Level = 10 )

      ! Eliminate the redundant nodes by default. 
      ! These are noded that depend on themselves.
      EliminateUnneeded = ListGetLogical( BC,&
          'Level Projector Eliminate Redundant Nodes',Found ) 
      IF(.NOT. Found ) EliminateUnneeded = .TRUE.

      IF( EliminateUnneeded ) THEN
        m = 0
        n = SUM( NodePerm )
        CALL Info('LevelProjector',&
            'Number of potential nodes in projector: '//TRIM(I2S(n)),Level=10)        
        ! Now eliminate the nodes which also occur in the other mesh
        ! These must be redundant edges
        DO i=1, SIZE(InvPerm2)
          j = InvPerm2(i) 
          IF( NodePerm(j) /= 0 ) THEN
            NodePerm(j) = 0
            !PRINT *,'Removing node:',j,Mesh % Nodes % x(j), Mesh % Nodes % y(j)
            m = m + 1
          END IF
        END DO
        IF( m > 0 ) THEN
          CALL Info('LevelProjector',&
              'Eliminating redundant nodes from projector: '//TRIM(I2S(m)),Level=10)
        END IF
      END IF
      
      IF( CreateDual ) THEN
        ALLOCATE( DualNodePerm( Mesh % NumberOfNodes ) )
        DualNodePerm = 0

        DO i=1,BMesh2 % NumberOfBulkElements
          Element => BMesh2 % Elements(i)        
          IF( Parallel ) THEN
            IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE          
          END IF
          DualNodePerm( InvPerm2( Element % NodeIndexes ) ) = 1
        END DO
                
        IF( EliminateUnneeded ) THEN
          m = 0
          n = SUM( DualNodePerm )
          CALL Info('LevelProjector',&
              'Number of potential nodes in dual projector: '//TRIM(I2S(n)),Level=10)        
          ! Now eliminate the nodes which also occur in the other mesh
          ! These must be redundant edges
          DO i=1, SIZE(InvPerm1)
            j = InvPerm1(i) 
            IF( DualNodePerm(j) /= 0 ) THEN
              DualNodePerm(j) = 0
              PRINT *,'Removing dual node:',j,Mesh % Nodes % x(j), Mesh % Nodes % y(j)
              m = m + 1
            END IF
          END DO
          IF( m > 0 ) THEN
            CALL Info('LevelProjector',&
                'Eliminating redundant dual nodes from projector: '//TRIM(I2S(m)),Level=10)
          END IF
        END IF
      END IF
      
      IF( ListCheckPresent( BC,'Level Projector Condition') ) THEN
        ALLOCATE( Cond( Mesh % MaxElementNodes ) )
        Cond = 1.0_dp
        m = 0
        DO i=1, BMesh1 % NumberOfBulkElements          
          Element => Mesh % Elements( BMesh1 % Elements(i) % ElementIndex )
          CurrentModel % CurrentElement => Element
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          Cond(1:n) = ListGetReal( BC,'Level Projector Condition', n, NodeIndexes )
          DO j=1,n
            k = NodeIndexes(j)
            IF( NodePerm(k) /= 0 ) THEN
              IF( Cond(j) < 0.0 ) THEN
                m = m + 1
                NodePerm(k) = 0 
              END IF
            END IF
          END DO
        END DO
        CALL Info('LevelProjector','Eliminated nodes with negative condition: '//&
            TRIM(I2S(m)),Level=10)        
        DEALLOCATE( Cond ) 
      END IF
      
      m = 0
      DO i=1,Mesh % NumberOfNodes
        IF( NodePerm(i) > 0 ) THEN
          m = m + 1
          NodePerm(i) = m
        END IF
      END DO
      
      CALL Info('LevelProjector',&
          'Number of active nodes in projector: '//TRIM(I2S(m)),Level=8)
      EdgeRow0 = m
      
      IF( CreateDual ) THEN
        m = 0
        DO i=1,Mesh % NumberOfNodes
          IF( DualNodePerm(i) > 0 ) THEN
            m = m + 1
            DualNodePerm(i) = m
          END IF
        END DO
        ALLOCATE( DualProjector % InvPerm(m) )
        DualProjector % InvPerm = 0

        IF( DoEdges ) THEN
          CALL Fatal('LevelProjector','Dual projector cannot handle edges!')
        END IF
      END IF
    ELSE
      EdgeRow0 = 0
    END IF
    ProjectorRows = EdgeRow0


    IF( DoEdges ) THEN
      ALLOCATE( EdgePerm( Mesh % NumberOfEdges ) )
      EdgePerm = 0

      ! Mark the edges for which the projector must be created for
      DO i=1, BMesh1 % NumberOfBulkElements

        ! in parallel only consider face elements that truly are part of this partition
        IF( Parallel ) THEN
          IF( BMesh1 % Elements(i) % PartIndex /= ParEnv % MyPe ) CYCLE          
        END IF

        DO j=1, BMesh1 % Elements(i) % TYPE % NumberOfEdges
          EdgePerm( BMesh1 % Elements(i) % EdgeIndexes(j) ) = 1
        END DO
      END DO

      EliminateUnneeded = ListGetLogical( BC,&
          'Level Projector Eliminate Redundant Edges',Found )
      IF(.NOT. Found ) EliminateUnneeded = .TRUE.

      IF( EliminateUnneeded ) THEN
        n = SUM( EdgePerm )
        CALL Info('LevelProjector',&
            'Number of potential edges in projector: '//TRIM(I2S(n)),Level=10)        
        ! Now eliminate the edges which also occur in the other mesh
        ! These must be redundant edges
        DO i=1, BMesh2 % NumberOfBulkElements
          DO j=1, BMesh2 % Elements(i) % TYPE % NumberOfEdges
            EdgePerm( BMesh2 % Elements(i) % EdgeIndexes(j) ) = 0
          END DO
        END DO

        IF( DoNodes ) THEN
          IF( ListGetLogical( BC,'Level Projector Eliminate Edges Greedy',Found ) ) THEN
            DO i=1, BMesh1 % NumberOfBulkElements
              DO j=1, BMesh1 % Elements(i) % TYPE % NumberOfEdges
                k = BMesh1 % Elements(i) % EdgeIndexes(j) 
                IF( ANY( NodePerm( Mesh % Edges(k) %  NodeIndexes ) == 0 ) ) THEN
                  EdgePerm( k ) = 0
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF

      m = 0
      DO i=1,Mesh % NumberOfEdges
        IF( EdgePerm(i) > 0 ) THEN
          m = m + 1
          EdgePerm(i) = m
        END IF
      END DO

      IF( EliminateUnneeded ) THEN
        CALL Info('LevelProjector',&
            'Eliminating redundant edges from projector: '//TRIM(I2S(n-m)),Level=10)
      END IF
      CALL Info('LevelProjector',&
          'Number of active edges in projector: '//TRIM(I2S(m)),Level=8)
      FaceRow0 = EdgeRow0 + m
      ProjectorRows = FaceRow0
      
      IF( PiolaVersion ) THEN
        ! Note: this might not work in parallel with halo since some of the face elements
        ! do not then belong to the slave boundary. 
        m = BMesh1 % NumberOfBulkElements
        CALL Info('LevelProjector',&
            'Number of active faces in projector: '//TRIM(I2S(m)),Level=8)
        ! There are two additional dofs associated with each face
        ProjectorRows = FaceRow0 + 2 * m
      END IF
    END IF

    CALL Info('LevelProjector',&
        'Max number of rows in projector: '//TRIM(I2S(ProjectorRows)),Level=10)
    ALLOCATE( Projector % InvPerm(ProjectorRows) )
    Projector % InvPerm = 0

    ! If after strong projectors there are still something undone they must 
    ! be dealt with the weak projectors. 
    SomethingUndone = .FALSE.

    ! If requested, create strong mapping for node dofs
    !------------------------------------------------------------------   
    IF( DoNodes ) THEN
      IF( StrongNodes ) THEN
        IF( GenericIntegrator ) THEN 
          CALL AddNodalProjectorStrongGeneric()
        ELSE
          CALL AddNodalProjectorStrongStrides()
        END IF
      ELSE
        ! If strong projector is applied they can deal with all nodal dofs
        SomethingUndone = .TRUE.
      END IF
    END IF

    ! If requested, create strong mapping for edge dofs
    !-------------------------------------------------------------
    EdgeBasis = .FALSE.
    IF( DoEdges ) THEN
      EdgeCol0 = Mesh % NumberOfNodes
      FaceCol0 = Mesh % NumberOfNodes + Mesh % NumberOfEdges

      IF( StrongLevelEdges .OR. StrongExtrudedEdges ) THEN
        CALL AddEdgeProjectorStrongStrides()
        ! Compute the unset edge dofs. 
        ! Some of the dofs may have been set by the strong projector. 
        m = 0
        DO i=1, Mesh % NumberOfEdges
          IF( EdgePerm(i) > 0 ) m = m + 1
        END DO
        IF( m > 0 ) THEN
          CALL Info('LevelProjector',&
              'Number of weak edges in projector: '//TRIM(I2S(m)),Level=10)      
        END IF
        IF( m > 0 .OR. PiolaVersion) THEN
          SomethingUndone = .TRUE.
          EdgeBasis = .TRUE.
        END IF
      ELSE
        SomethingUndone = .TRUE.
        EdgeBasis = .TRUE.
      END IF      
    END IF

    ! And the the rest
    !-------------------------------------------------------------
    IF( SomethingUndone ) THEN      
      IF( MeshDim == 2 ) THEN
        CALL AddProjectorWeak1D()
      ELSE IF( GenericIntegrator ) THEN
        CALL AddProjectorWeakGeneric()
      ELSE
        CALL AddProjectorWeakStrides()
      END IF
    END IF

    ! Now change the matrix format to CRS from list matrix
    !--------------------------------------------------------------
    CALL List_toCRSMatrix(Projector)
    CALL CRS_SortMatrix(Projector,.TRUE.)
    CALL Info('LevelProjector','Number of rows in projector: '&
        //TRIM(I2S(Projector % NumberOfRows)),Level=12)
    CALL Info('LevelProjector','Number of entries in projector: '&
        //TRIM(I2S(SIZE(Projector % Values))),Level=12)
  

    IF(ASSOCIATED(Projector % Child)) THEN
      CALL List_toCRSMatrix(Projector % Child)
      CALL CRS_SortMatrix(Projector % Child,.TRUE.)
    END IF

    IF( CreateDual ) THEN
      CALL List_toCRSMatrix(DualProjector)
      CALL CRS_SortMatrix(DualProjector,.TRUE.)
    END IF
    
    IF( DoNodes ) DEALLOCATE( NodePerm )
    IF( CreateDual .AND. DoNodes ) DEALLOCATE( DualNodePerm )
    IF( DoEdges ) DEALLOCATE( EdgePerm )

    m = COUNT( Projector % InvPerm  == 0 ) 
    IF( m > 0 ) THEN
      CALL Warn('LevelProjector','Projector % InvPerm not set in for dofs: '//TRIM(I2S(m)))
    END IF

    CALL Info('LevelProjector','Projector created',Level=10)


  CONTAINS

    ! Currently the target mesh is assumed to be include only cartesian elements
    ! Check the angle in the elements. When we know the target mesh is cartesian
    ! we can reduce the error control in the other parts of the code. 
    !----------------------------------------------------------------------------
    FUNCTION CheckMeshSkew(BMesh, NotAllQuads) RESULT( MaxSkew )

      TYPE(Mesh_t),POINTER :: BMesh
      REAL(KIND=dp) :: MaxSkew
      LOGICAL :: NotAllQuads

      INTEGER :: i,j,n,indM,k,knext,kprev
      TYPE(Element_t), POINTER :: ElementM
      TYPE(Nodes_t) :: NodesM
      REAL(KIND=dp) :: e1(2),e2(2),DotProdM, PhiM
      INTEGER, POINTER :: IndexesM(:)

      CALL Info('LevelProjector','Checking mesh skew')

      n = 4
      ALLOCATE( NodesM % x(n), NodesM % y(n) )
      MaxSkew = 0.0_dp
      NotAllQuads = .FALSE.
      
      j = 0
      DO indM=1,BMesh % NumberOfBulkElements
        
        ElementM => BMesh % Elements(indM)        
        IF( ElementM % TYPE % ElementCode / 100 /= 4 ) THEN
          NotAllQuads = .TRUE.
        END IF
        IndexesM => ElementM % NodeIndexes
        NodesM % y(1:n) = BMesh % Nodes % y(IndexesM(1:n))
        NodesM % x(1:n) = BMesh % Nodes % x(IndexesM(1:n))
        
        ! Transfer into real length units instead of angles
        ! This gives right balance between x and y -directions. 
        NodesM % x(1:n) = ArcCoeff * NodesM % x(1:n)
        
        ! Make unit vectors of the edge
        DO k = 1, n
          knext = MODULO(k,n)+1
          kprev = MODULO(n+k-2,n)+1
          
          e1(1) = NodesM % x(knext) - NodesM % x(k) 
          e1(2) = NodesM % y(knext) - NodesM % y(k) 
          
          e2(1) = NodesM % x(kprev) - NodesM % x(k) 
          e2(2) = NodesM % y(kprev) - NodesM % y(k) 
          
          e1 = e1 / SQRT( SUM( e1**2) )
          e2 = e2 / SQRT( SUM( e2**2) )
          
          ! dot product of the unit vectors
          DotProdM = SUM( e1 * e2 )
          
          ! Cosine angle in degrees        
          PhiM = ACOS( DotProdM ) 
          MaxSkew = MAX( MaxSkew, ABS ( ABS( PhiM ) - PI/2 ) )
        END DO
      END DO

      ! Move to degrees and give the tolerance in them
      MaxSkew = MaxSkew * 180.0 / PI
        
100   DEALLOCATE( NodesM % x, NodesM % y )

    END FUNCTION CheckMeshSkew
      

    !-------------------------------------------------------------------------------------
    ! Create projector for nodes on the strides directly from a linear 
    ! combination of two nodes. This approach minimizes the size of the projector
    ! and also minimizes the need for parallel communication.
    !-------------------------------------------------------------------------------------
    SUBROUTINE AddNodalProjectorStrongStrides()

      TYPE(Element_t), POINTER :: ElementM
      INTEGER, POINTER :: IndexesM(:)
      INTEGER :: ncoeff, coeffi(2),sgn0, ind, indm, j1, j2, j3, Nundefined
      REAL(KIND=dp) :: x1, y1, x2, y2, xmin, xmax, xminm, xmaxm, Dist, MinDist
      REAL(KIND=dp) :: coeff(2), val, xm1, xm2, xm3
      INTEGER, POINTER :: EdgeMap(:,:)
      TYPE(Nodes_t) :: NodesM
      LOGICAL :: LeftCircle

      CALL Info('LevelProjector','Creating strong stride projector for nodal dofs',Level=10)

      n = Mesh % MaxElementNodes
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      NodesM % z = 0.0_dp

      ! By construction there is always two components in the projector for the nodes. 
      ncoeff = 2
      coeffi = 0
      sgn0 = 1
      Nundefined = 0

      ! This flag tells if we're working with a full circle and the problematic part of 
      ! the circle with the discontinuity in the angle. 
      LeftCircle = .FALSE.

      DO ind=1,BMesh1 % NumberOfNodes

        nrow = NodePerm( InvPerm1( ind ) )
        IF( nrow == 0 ) CYCLE
        NodePerm( InvPerm1( ind ) ) = 0
        Projector % InvPerm(nrow) = InvPerm1(ind)

        Found = .FALSE.
        x1 = BMesh1 % Nodes % x(ind)
        y1 = BMesh1 % Nodes % y(ind)
        sgn0 = 1
        coeff = 0.0_dp
        MinDist = HUGE( MinDist )

        IF( Repeating ) THEN
          Nrange = FLOOR( (x1-XMinAll) / XRange )
          x1 = x1 - Nrange * XRange
          
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
        ELSE IF( FullCircle ) THEN
          LeftCircle = ABS( x1 ) > 90.0
          IF( LeftCircle ) THEN
            IF( x1 < 0.0 ) x1 = x1 + 360.0
          END IF
        END IF

        ! If the projector is of style Px+Qx=0 then
        ! and the negative sign, otherwise let the initial sign be.
        IF( SelfProject ) sgn0 = -sgn0
        
        ! Currently a cheap n^2 loop but it could be improved
        ! Looping over master elements. Look for constant-y strides only. 
        !--------------------------------------------------------------------
        DO indM = 1, BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)
          n = ElementM % TYPE % NumberOfNodes        
          IndexesM => ElementM % NodeIndexes
          
          ! Quick tests to save time
          ! Element must have nodes at the right level
          NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))           
          IF( ALL( ABS( NodesM % y(1:n) - y1 ) > YTol ) ) CYCLE

          ! The x nodes should be in the interval
          NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))

          ! Transform the master element on-the-fly around the problematic angle
          IF( LeftCircle ) THEN
            ! The master nodes are all on right
            IF( ALL( ABS( NodesM % x(1:n) ) - 90.0_dp < Xtol ) ) CYCLE
            DO j=1,n
              IF( NodesM % x(j) < 0.0 ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          xmaxm = MAXVAL( NodesM % x(1:n) )
          xminm = MINVAL( NodesM % x(1:n) )

          ! Eliminate this special case since it could otherwise give a faulty hit
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > 180.0 ) CYCLE
          END IF

          Dist = MAX( x1-xmaxm, xminm-x1 ) 

          ! Mark the minimum distance if this would happen to be a problematic node
          MinDist = MIN( Dist, MinDist )

          IF( Dist > Xtol ) CYCLE

          ! Ok, this may be a proper element, now just find the two nodes
          ! needed for the mapping on the same stride. Basically this means 
          ! finding the correct edge but we don't need to use the data structure for that. 
          ! For 1D edge element this is trivial, note however that only 1st degree projection is used!
          j1 = 0; j2 = 0; j3 = 0
          IF( n <= 3 ) THEN
            j1 = 1 
            j2 = 2
            IF( n == 3 ) j3 = 3
          ELSE
            DO j=1,n
              IF( ABS( NodesM % y(j) - y1 ) > YTol ) CYCLE
              IF( j1 == 0 ) THEN
                j1 = j
              ELSE IF( j2 == 0 ) THEN
                j2 = j
              ELSE
                j3 = j
                ! This means that for higher order edges only three nodes are used
                EXIT
              END IF
            END DO
            IF( j2 == 0 ) CALL Warn('LevelProjector','Could not locate an edge consistently!')
          END IF

          ! The node to map must be in interval, x1 \in [xm1,xm2]
          IF( NodesM % x(j1) > NodesM % x(j2) ) THEN
             j = j2; j2 = j1; j1 = j
          END IF
          xm1 = NodesM % x(j1)
          xm2 = NodesM % x(j2)          

          ! We are at interval [xm1,xm2] now choose either [xm1,xm3] or [xm3,xm2]
          IF( j3 > 0 ) THEN
             xm3 = NodesM % x(j3)          
             IF( x1 > xm3 ) THEN
                j1 = j3; xm1 = xm3
             ELSE 
                j2 = j3; xm2 = xm3
             END IF
          END IF
          
          ! Ok, the last check, this might fail if the element had skew even though the 
          ! quick test is successfull! Then the left and right edge may have different range.
          Dist = MAX( x1-xm2, xm1-x1 )
          IF( Dist > Xtol ) CYCLE

          ! When we have the correct edge, the mapping is trivial.
          ! The sum of weights of the projectors is set to one. 
          IF( ABS(xm1-xm2) < TINY(xm1) ) THEN
            CALL Warn('LevelProjector','Degenerated edge?')
            PRINT *,'ind',ind,x1,y1,xm1,xm2,j1,j2,j3
            PRINT *,'x:',NodesM % x(1:n)
            PRINT *,'y:',NodesM % y(1:n)
            coeff(1) = 0.5_dp
          ELSE
            coeff(1) = (xm2-x1)/(xm2-xm1) 
          END IF
          coeff(2) = 1.0_dp - coeff(1)

          coeffi(1) = IndexesM(j1)
          coeffi(2) = IndexesM(j2)

          Found = .TRUE.
          
          ! If we really exactly between [xm1,xm2] then we may finish the search for good
          IF( Dist < EPSILON( Dist ) ) EXIT
        END DO

        IF(.NOT. Found ) THEN
          Nundefined = Nundefined + 1
          WRITE( Message,'(A,2I8,3ES12.3)') 'Problematic node: ',&
              ind,ParEnv % MyPe,x1,y1,MinDist
          CALL Warn('LevelProjector',Message)
          CYCLE
        END IF

        IF( SelfProject ) THEN
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPerm1(ind), NodeCoeff ) 
        END IF

        ! The scaling of the projector entries is used, for example, 
        ! to allow antiperiodic projectors. 
        Coeff(1:ncoeff) = sgn0 * Coeff(1:ncoeff)

        ! The projection weights
        DO j=1,ncoeff 

          val = Coeff(j)
          ! Skip too small projector entries
          IF( ABS( val ) < 1.0d-12 ) CYCLE

          ! Use the permutation to revert to original dofs
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPerm2(coeffi(j)), NodeScale * NodeCoeff * val ) 
        END DO

      END DO

      IF( Nundefined > 0 ) THEN
        CALL Warn('LevelProjector',&
            'Nodes could not be determined by any edge: '//TRIM(I2S(Nundefined)))          
      END IF

      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )


    END SUBROUTINE AddNodalProjectorStrongStrides
    !---------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------
    ! Adds a nodal projector assuming generic 2D mesh. 
    ! Otherwise should give same results as the one before. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddNodalProjectorStrongGeneric()

      TYPE(Element_t), POINTER :: ElementM
      INTEGER, POINTER :: IndexesM(:), coeffi(:)
      REAL(KIND=dp), POINTER :: Basis(:),coeff(:)
      INTEGER :: n, nM, ncoeff, sgn0, ind, indm, j1, j2, j3, Nundefined
      REAL(KIND=dp) :: x1, y1, z1, xmin, xmax, xminm, xmaxm, ymaxm, yminm, &
          Dist, MaxMinBasis, detJ, ArcTol, ArcRange
      REAL(KIND=dp) :: val, u, v, w
      TYPE(Nodes_t) :: NodesM
      LOGICAL :: LeftCircle, Found, Stat

      CALL Info('LevelProjector','Creating strong generic projector for nodal dofs',Level=10)

      n = Mesh % MaxElementNodes
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n), Basis(n), coeff(n), coeffi(n) )
      NodesM % z = 0.0_dp

      ncoeff = 0
      coeffi = 0
      sgn0 = 1
      Nundefined = 0
      z1 = 0.0_dp

      ArcTol = ArcCoeff * Xtol
      ArcRange = ArcCoeff * Xrange 

      ! This flag tells if we're working with a full circle and the problematic part of 
      ! the circle with the discontinuity in the angle. 
      LeftCircle = .FALSE.

      DO ind=1,BMesh1 % NumberOfNodes

        nrow = NodePerm( InvPerm1( ind ) )
        IF( nrow == 0 ) CYCLE
        NodePerm( InvPerm1( ind ) ) = 0
        Projector % InvPerm(nrow) = InvPerm1(ind)

        Found = .FALSE.
        x1 = ArcCoeff * BMesh1 % Nodes % x(ind)
        y1 = BMesh1 % Nodes % y(ind)
        IF( HaveMaxDistance ) THEN
          z1 = BMesh1 % Nodes % z(ind)
        END IF

        sgn0 = 1
        coeff = 0.0_dp
        MaxMinBasis = -HUGE(MaxMinBasis)

        IF( FullCircle ) THEN
          LeftCircle = ABS( x1 ) > ArcCoeff * 90.0
          IF( LeftCircle ) THEN
            IF( x1 < 0.0 ) x1 = x1 + ArcCoeff * 360.0
          END IF
        END IF

        ! If the projector is of style Px+Qx=0 then
        ! and the negative sign, otherwise let the initial sign be.
        IF( SelfProject ) sgn0 = -sgn0
        
        ! Currently a cheap n^2 loop but it could be improved
        ! Looping over master elements. Look for constant-y strides only. 
        !--------------------------------------------------------------------
        DO indM = 1, BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)
          nM = ElementM % TYPE % NumberOfNodes        
          IndexesM => ElementM % NodeIndexes

          IF( HaveMaxDistance ) THEN
            IF( MINVAL( ABS( BMesh2 % Nodes % z(IndexesM(1:nM)) - z1 ) ) > MaxDistance ) CYCLE          
          END IF
          
          ! Quick tests to save time
          NodesM % y(1:nM) = BMesh2 % Nodes % y(IndexesM(1:nM))           
          ymaxm = MAXVAL( NodesM % y(1:nM) )
          yminm = MINVAL( NodesM % y(1:nM) )

          Dist = MAX( y1-ymaxm, yminm-y1 ) 
          IF( Dist > Ytol ) CYCLE

          ! The x nodes should be in the interval
          NodesM % x(1:nM) = BMesh2 % Nodes % x(IndexesM(1:nM))

          ! Transform the master element on-the-fly around the problematic angle
          ! Full 2D circle is never repeating
          IF( LeftCircle ) THEN
            ! The master nodes are all on right
            IF( ALL( ABS( NodesM % x(1:nM) ) - ArcCoeff * 90.0_dp < ArcTol ) ) CYCLE
            DO j=1,nM
              IF( NodesM % x(j) < 0.0 ) NodesM % x(j) = NodesM % x(j) + ArcCoeff * 360.0_dp
            END DO
          END IF
          
          xmaxm = MAXVAL( NodesM % x(1:nM) )
          xminm = MINVAL( NodesM % x(1:nM) )

          ! Eliminate this special case since it could otherwise give a faulty hit
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0 ) CYCLE
          END IF

          IF( Repeating ) THEN
            Nrange = FLOOR( (xmaxm-x1) / XRange )
            IF( Nrange /= 0 ) THEN
              xminm = xminm - Nrange * ArcRange
              xmaxm = xmaxm - Nrange * ArcRange
              NodesM % x(1:nM) = NodesM % x(1:nM) - NRange * ArcRange 
            END IF

            ! Check whether there could be a intersection in an other interval as well
            IF( xminm + ArcRange < x1 + ArcTol ) THEN
              Nrange2 = 1
            ELSE
              Nrange2 = 0
            END IF
          END IF

100       Dist = MAX( x1-xmaxm, xminm-x1 ) 

          IF( Dist < Xtol ) THEN
            ! Integration point at the slave element
            CALL GlobalToLocal( u, v, w, x1, y1, z1, ElementM, NodesM )              
            stat = ElementInfo( ElementM, NodesM, u, v, w, detJ, Basis )
            
            IF( MINVAL( Basis(1:nM) ) > MaxMinBasis ) THEN
              MaxMinBasis = MINVAL( Basis(1:nM) )
              ncoeff = nM
              coeff(1:nM) = Basis(1:nM)
              coeffi(1:nM) = IndexesM(1:nM)
              Found = ( MaxMinBasis >= -1.0d-12 )
            END IF
         
            IF( Found ) EXIT
          END IF
          
          IF( Repeating ) THEN
            IF( NRange2 /= 0 ) THEN
              xminm = xminm + ArcCoeff * Nrange2 * ArcRange
              xmaxm = xmaxm + ArcCoeff * Nrange2 * ArcRange
              NodesM % x(1:n) = NodesM % x(1:n) + NRange2 * ArcRange 
              NRange = NRange + NRange2
              NRange2 = 0
              GOTO 100
            END IF
          END IF         

        END DO

        IF(.NOT. Found ) THEN
          IF( MaxMinBasis > -1.0e-6 ) THEN
            CALL Info('LevelProjector',Message,Level=8)
            Found = .TRUE.
          ELSE
            Nundefined = Nundefined + 1
            IF( .NOT. HaveMaxDistance ) THEN
              WRITE( Message,'(A,2I8,3ES12.3)') 'Problematic node: ',&
                  ind,ParEnv % MyPe,x1,y1,MaxMinBasis
              CALL Warn('LevelProjector',Message )
            END IF
          END IF
        END IF

        IF( Found ) THEN
          IF( SelfProject ) THEN
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                InvPerm1(ind), NodeCoeff ) 
          END IF
          
          ! The scaling of the projector entries is used, for example, 
          ! to allow antiperiodic projectors. 
          Coeff(1:ncoeff) = sgn0 * Coeff(1:ncoeff)
          
          ! Add the projection weights to the matrix
          DO j=1,ncoeff 
            
            val = Coeff(j)
            ! Skip too small projector entries
            ! These really should sum to one we now the limit quite well
            IF( ABS( val ) < 1.0d-8 ) CYCLE
            
            ! Use the permutation to revert to original dofs
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                InvPerm2(coeffi(j)), NodeScale * NodeCoeff * val ) 
          END DO
        END IF
        
      END DO

      IF( Nundefined > 0 ) THEN
        IF( HaveMaxDistance ) THEN
          CALL Info('LevelProjector',&
              'Nodes could not be found in any element: '//TRIM(I2S(Nundefined)))          
        ELSE
          CALL Warn('LevelProjector',&
              'Nodes could not be found in any element: '//TRIM(I2S(Nundefined)))          
        END IF
      END IF

      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z, Basis, coeffi, coeff )


    END SUBROUTINE AddNodalProjectorStrongGeneric
    !---------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------
    ! Create a projector for edges directly. This minmizes the size of the projector 
    ! but may result to numerically inferior projector compared to the weak projector.
    ! It seems to be ok for unskewed geometries where the simplest edge elements work 
    ! well. For skewed geometries the solution does not easily seem to be compatible
    ! with the strong projector. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddEdgeProjectorStrongStrides()

      INTEGER :: ind, indm, eind, eindm, k1, k2, km1, km2, sgn0, coeffi(100), &
          ncoeff, dncoeff, ncoeff0, i1, i2, j1, j2, Nundefined, NoSkewed, SkewPart
      TYPE(Element_t), POINTER :: Element, ElementM
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      TYPE(Nodes_t) :: NodesM, Nodes 
      INTEGER, POINTER :: EdgeMap(:,:),EdgeMapM(:,:)
      REAL(KIND=dp) :: xm1, xm2, ym1, ym2, coeff(100), signs(100), wsum, minwsum, maxwsum, val, &
          x1o, y1o, x2o, y2o, cskew, sedge
      REAL(KIND=dp) :: x1, y1, x2, y2, xmin, xmax, xminm, xmaxm, ymin, ymax, yminm, ymaxm, xmean, &
          dx,dy,Xeps
      LOGICAL :: YConst, YConstM, XConst, XConstM, EdgeReady, Repeated, LeftCircle, &
          SkewEdge, AtRangeLimit


      CALL Info('LevelProjector','Creating strong stride projector for edges assuming strides',Level=10)

      n = Mesh % NumberOfEdges
      IF( n == 0 ) RETURN      

      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      Nodes % z = 0.0_dp
      NodesM % z = 0.0_dp

      minwsum = HUGE( minwsum ) 
      maxwsum = 0.0_dp
      NoSkewed = 0
      Nundefined = 0
      LeftCircle = .FALSE.
      Xeps = EPSILON( Xeps )
      AtRangeLimit = .FALSE.

      DO ind=1,BMesh1 % NumberOfBulkElements
        
        Element => BMesh1 % Elements(ind)        
        EdgeMap => LGetEdgeMap( Element % TYPE % ElementCode / 100)

        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))

        dx = MAXVAL( Nodes % x(1:n)) - MINVAL(Nodes % x(1:n))
        dy = MAXVAL( Nodes % y(1:n)) - MINVAL(Nodes % y(1:n))

        ! Go through combinations of edges and find the edges for which the 
        ! indexes are the same. 
        DO i = 1,Element % TYPE % NumberOfEdges
          
          eind = Element % EdgeIndexes(i)
          IF( EdgePerm(eind) == 0 ) CYCLE

          nrow = EdgeRow0 + EdgePerm(eind) 
          
          ! Get the nodes of the edge
          i1 = EdgeMap(i,1) 
          i2 = EdgeMap(i,2)

          k1 = Indexes( i1 )
          k2 = Indexes( i2 )

          ! The coordinates of the edge
          x1 = Nodes % x(i1)
          y1 = Nodes % y(i1)

          x2 = Nodes % x(i2)
          y2 = Nodes % y(i2)

          YConst = ( ABS(y2-y1) < RadTol * dy )
          XConst = ( ABS(x2-x1) < RadTol * dx )

          SkewEdge = .FALSE.
          cskew = 1.0_dp
          
          IF( YConst ) THEN
            IF( .NOT. StrongLevelEdges ) CYCLE         
          ELSE IF( XConst ) THEN
            IF( .NOT. StrongExtrudedEdges ) CYCLE
          ELSE
            !print *,'skewed edge: ',ParEnv % MyPe,x1,x2,y1,y2,dx,dy
            !print *,'tol:',ABS(y2-y1)/dy,ABS(x2-x1)/dx,RadTol

            NoSkewed = NoSkewed + 1
            SkewEdge = .TRUE.
            IF(.NOT. StrongSkewEdges) CYCLE
          END IF
          

          ! Numbering of global indexes is needed to ensure correct direction 
          ! of the edge dofs. Basically the InvPerm could be used also in serial
          ! but the order of numbering is maintained when the reduced mesh is created. 
          IF(Parallel) THEN
            k1 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm1(k1))
            k2 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm1(k2))
          END IF
          ncoeff = 0 

          IF( SkewEdge ) THEN
            SkewPart = 0
            sedge = SQRT(ArcCoeff**2*(x1-x2)**2 + (y1-y2)**2)
            x1o = x1
            y1o = y1
            x2o = x2
            y2o = y2
          END IF

          ! This is mainly a test branch for skewed quadrilaters.
          ! It is based on the composition of a skewed edge into 
          ! four cartesian vectors oriented along x or y -axis. 
          ! Unfortunately the resulting projector does not seem to be 
          ! numerically favourable.           
50        IF( SkewEdge ) THEN
            IF( SkewPart < 2 ) THEN
              XConst = .TRUE.
              YConst = .FALSE.
              IF( SkewPart == 1 ) THEN
                x1 = (3.0*x1o + x2o) / 4.0_dp
              ELSE
                x1 = (x1o + 3*x2o) / 4.0_dp
              END IF
              x2 = x1
              y1 = y1o
              y2 = y2o
              cskew = 0.5 * ABS(y1-y2) / sedge
            ELSE 
              XConst = .FALSE.
              YConst = .TRUE.
              IF( SkewPart == 2 ) THEN
                x1 = x1o
                x2 = (x1o + x2o) / 2.0_dp
                y1 = y1o
                y2 = y1o
              ELSE
                x1 = (x1o + x2o) / 2.0_dp
                x2 = x2o
                y1 = y2o
                y2 = y2o
              END IF
              cskew = ArcCoeff * ABS(x1-x2) / sedge
            END IF
          END IF 

          ncoeff0 = ncoeff
          dncoeff = 0
          Repeated = .FALSE.

          ! If the edge might be treated in two periodic parts 
          ! then here study whether this is the case (Nrange2 /= 0). 
          IF( Repeating ) THEN
            Nrange = FLOOR( (x1-XMinAll) / XRange )
            x1 = x1 - Nrange * XRange
            x2 = x2 - Nrange * XRange
            
            IF( x2 > XMaxAll ) THEN
              Nrange2 = 1
            ELSE IF( x2 < XMinAll ) THEN
              Nrange2 = -1
            ELSE
              Nrange2 = 0
            END IF
          ELSE IF( FullCircle ) THEN
            ! If we have a full circle then treat the left-hand-side
            ! differently in order to circumvent the discontinuity of the
            ! angle at 180 degrees. 
            LeftCircle = ( ABS(x1) > 90.0 .AND. ABS(x2) > 90.0 )
            IF( LeftCircle ) THEN
              IF( x1 < 0.0_dp ) x1 = x1 + 360.0_dp
              IF( x2 < 0.0_dp ) x2 = x2 + 360.0_dp
            END IF
          END IF

          EdgeReady = .FALSE.
100       sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
          
          IF( SelfProject ) sgn0 = -sgn0
          
          xmin = MIN(x1,x2)
          xmax = MAX(x1,x2)
          ymin = MIN(y1,y2)
          ymax = MAX(y1,y2)
          xmean = (x1+x2) / 2.0_dp


          ! If the mesh is not repeating there is a risk that we don't exactly hit the start 
          ! or end of the range. Therefore grow the tolerance close to the ends. 
          IF(.NOT. ( Repeating .OR. FullCircle ) ) THEN
            IF ( xmax < XminAll + Xtol .OR. xmin > XmaxAll - Xtol ) THEN
              Xeps = Xtol 
            ELSE
              Xeps = EPSILON( Xeps ) 
            END IF
          END IF

          
          ! Currently a n^2 loop but it could be improved
          !--------------------------------------------------------------------
          DO indm=1,BMesh2 % NumberOfBulkElements
            
            ElementM => BMesh2 % Elements(indm)        
            n = ElementM % TYPE % NumberOfNodes        
            IndexesM => ElementM % NodeIndexes(1:n)
            
            ! Make first some coarse tests to eliminate most of the candidate elements
            ! The y nodes should always have an exact fit
            NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))           
            IF( MINVAL( ABS( ymin - NodesM % y(1:n) ) ) > YTol ) CYCLE
            IF(.NOT. YConst ) THEN
              IF( MINVAL( ABS( ymax - NodesM % y(1:n) ) ) > YTol ) CYCLE
            END IF
            
            NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))
            
            ! If we have a full circle then treat the left part differently
            IF( LeftCircle ) THEN
              IF( ALL( ABS( NodesM % x(1:n) ) - 90.0 < Xtol ) ) CYCLE
              DO j=1,n
                IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
              END DO
            END IF
            
            ! The x nodes should be in the interval
            xminm = MINVAL( NodesM % x(1:n) ) 
            xmaxm = MAXVAL( NodesM % x(1:n) ) 
            
            IF( xminm > xmax + Xeps ) CYCLE
            IF( xmaxm < xmin - Xeps ) CYCLE 
            
            ! Eliminate this special case since it could otherwise give a faulty hit
            IF( FullCircle .AND. .NOT. LeftCircle ) THEN
              IF( xmaxm - xminm > 180.0 ) CYCLE
            END IF

            yminm = MINVAL( NodesM % y(1:n) ) 
            ymaxm = MAXVAL( NodesM % y(1:n) ) 
            
            ! Ok, we have found a candicate face that will probably have some hits       
            EdgeMapM => LGetEdgeMap( ElementM % TYPE % ElementCode / 100)        
            
            ! Go through combinations of edges and find the edges for which the 
            ! indexes are the same. 
            DO j = 1,ElementM % TYPE % NumberOfEdges

              eindm = ElementM % EdgeIndexes(j)
              
              ! Eliminate the possibilitity that the same edge is accounted for twice
              ! in two different boundary elements. 
              IF( ANY( coeffi(ncoeff0+1:ncoeff) == eindm ) ) CYCLE
              
              j1 = EdgeMap(j,1)
              j2 = EdgeMap(j,2)
              
              km1 = IndexesM( j1 )
              km2 = IndexesM( j2 )
              
              ym1 = NodesM % y(j1)
              ym2 = NodesM % y(j2)
              
              xm1 = NodesM % x(j1)
              xm2 = NodesM % x(j2)
              
              ! The target mesh has already been checked that the elements are rectangular so 
              ! the edges must be have either constant y or x.
              YConstM = ( ABS(ym2-ym1) / (ymaxm-yminm) < ABS(xm2-xm1) / (xmaxm-xminm) )
              XConstM = .NOT. YConstM
              
              ! Either both are lateral edges, or both are vertical
              IF( .NOT. ( ( YConst .AND. YConstM ) .OR. ( XConst .AND. XConstM ) ) ) THEN
                CYCLE
              END IF
              
              ! sign depends on the direction and order of global numbering
              IF(Parallel) THEN
                km1 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm2(km1))
                km2 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm2(km2))
              END IF
              
              IF( YConst ) THEN
                IF( ABS( y1 - ym1 ) > YTol ) CYCLE
                
                ! Check whether the range of master x has a union with the slave x
                xmaxm = MAX( xm1, xm2 ) 
                IF( xmaxm < xmin ) CYCLE
                
                xminm = MIN( xm1, xm2 ) 
                IF( xminm > xmax ) CYCLE

                ! Ok, we have a hit register it 
                ncoeff = ncoeff + 1
                coeffi(ncoeff) = eindm

                ! weight depends on the relative fraction of overlapping
                IF( ABS( xmax-xmin) < TINY( xmax ) ) THEN
                  CALL Warn('LevelProjector','Degenerated edge 2?')
                  coeff(ncoeff) = cskew * 1.0_dp
                ELSE
                  coeff(ncoeff) = cskew * (MIN(xmaxm,xmax)-MAX(xminm,xmin))/(xmax-xmin)
                END IF

                ! this sets the sign which should be consistant 
                IF( (x1-x2)*(xm1-xm2)*(k1-k2)*(km1-km2) > 0.0_dp ) THEN
                  signs(ncoeff) = sgn0
                ELSE
                  signs(ncoeff) = -sgn0
                END IF

                ! There can be only one lateral edge hit for each element
                EXIT 
              ELSE
                dncoeff = dncoeff + 1
                ncoeff = ncoeff + 1

                IF( (y1-y2)*(ym1-ym2)*(k1-k2)*(km1-km2) > 0.0_dp ) THEN
                  signs(ncoeff) = sgn0 
                ELSE
                  signs(ncoeff) = -sgn0
                END IF

                coeffi(ncoeff) = eindm
                ! note: temporarily save the coordinate to the coefficient!
                coeff(ncoeff) = ( xm1 + xm2 ) / 2.0_dp
              END IF
            END DO

            IF( .NOT. SkewEdge ) THEN
              IF( YConst ) THEN
                ! Test whether the sum of coefficients has already reached unity
                wsum = SUM( coeff(1:ncoeff) )
                EdgeReady = ( 1.0_dp - wsum < 1.0d-12 ) 
              ELSE IF( XConst ) THEN                       
                ! If edge was found both on left and right there is no need to continue search
                EdgeReady = ( dncoeff == 2 ) 
              END IF
              IF( EdgeReady ) EXIT
            END IF
          END DO

          IF( YConst ) THEN
            ! For constant y check the 2nd part 
            ! and redo the search if it is active. 
            IF( Repeating ) THEN
              IF( NRange2 /= 0 ) THEN
                x1 = x1 - NRange2 * XRange
                x2 = x2 - NRange2 * XRange
                NRange = NRange + NRange2
                NRange2 = 0
                Repeated = .TRUE.
                GOTO 100
              END IF
            END IF
          ELSE
            ! Here there can be a second part if a proper hit was not found 
            ! due to some epsilon rules.
            IF( SkewEdge ) THEN
              IF( dncoeff == 1 ) THEN
                coeff(ncoeff) = cskew * 1.0_dp
              ELSE IF( dncoeff == 2 ) THEN
                xm1 = coeff(ncoeff-1)
                xm2 = coeff(ncoeff)
                
                IF( ABS( xm2-xm1) < TINY( xm2 ) ) THEN
                  CALL Warn('LevelProjector','Degenerated edge 3?')
                  coeff(ncoeff-1) = cskew * 0.5_dp
                ELSE
                  coeff(ncoeff-1) = cskew * ABS((xm2-xmean)/(xm2-xm1))
                END IF
                coeff(ncoeff) = cskew * 1.0_dp - coeff(1)
              END IF
            ELSE
              IF( ncoeff == 1 ) THEN
                coeff(1) = 1.0_dp
              ELSE IF( ncoeff >= 2 ) THEN
                IF( ncoeff > 2 ) THEN
                  CALL Warn('LevelProjector',&
                       'There should not be more than two target edges: '//TRIM(I2S(ncoeff))) 
                END IF
                xm1 = coeff(1)
                xm2 = coeff(2)
                IF( ABS( xm2-xm1) < TINY( xm2 ) ) THEN
                  CALL Warn('LevelProjector','Degenerated edge 3?')
                  coeff(1) = 0.5_dp
                ELSE
                  coeff(1) = ABS((xm2-xmean)/(xm2-xm1))
                END IF
                coeff(2) = 1.0_dp - coeff(1)
              END IF
            END IF

            wsum = SUM( coeff(1:ncoeff) )
          END IF

          ! Skewed edge is treated in four different parts (0,1,2,3)
          ! Go for the next part, if not finished. 
          IF( SkewEdge ) THEN
            IF( SkewPart < 3 ) THEN
              SkewPart = SkewPart + 1
              GOTO 50
            END IF
          END IF
              
          IF( ncoeff == 0 ) THEN
            Nundefined = Nundefined + 1
            WRITE( Message,'(A,2I8,4ES12.3)') 'Problematic edge: ',&
                eind,ParEnv % MyPe,x1,x2,y1,y2
            CALL Warn('LevelProjector', Message )
            WRITE( Message,'(A,I8,3L4,4ES12.3)') 'Bounding box: ',&
                eind,XConst,YConst,Repeating,XminAll,XmaxAll,YminAll,YmaxAll
            CALL Warn('LevelProjector', Message )
            CYCLE
          END IF

          wsum = SUM( ABS( coeff(1:ncoeff) ) )
          minwsum = MIN( minwsum, wsum ) 
          maxwsum = MAX( maxwsum, wsum ) 

          ! In skewed edges the sum of weights may be different from 1 but otherwise
          ! it should be very close to one. 
!          IF( ABS(wsum) < 0.999 .OR. ( ABS(wsum) > 1.001 .AND. .NOT. SkewEdge ) ) THEN
          IF(.FALSE.) THEN
            PRINT *,'*********************'
            PRINT *,'wsum',eind,ncoeff,wsum,Repeated
            PRINT *,'x coords:',x1,x2
            PRINT *,'y coords:',y1,y2
            PRINT *,'xm:',xm1,xm2
            PRINT *,'ym:',ym1,ym2
            PRINT *,'xm coords:',NodesM % x(1:4)
            PRINT *,'ym coords:',NodesM % y(1:4)
            PRINT *,'Const:',XConst,YConst,XConstM,YConstM
            PRINT *,'coeff:',ncoeff,coeff(1:ncoeff),coeffi(1:ncoeff)
          END IF

          ! Mark that this is set so it don't need to be set again
          EdgePerm(eind) = 0

          ! Ok, we found a true projector entry
          Projector % InvPerm(nrow) = EdgeCol0 + eind

          ! The reference to the edge to be projected
          IF( SelfProject ) THEN
            val = 1.0_dp
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                EdgeCol0 + eind, EdgeCoeff * val ) 
          END IF

          ! The scaling can be used to create antiperiodic projectors, for example. 
          Coeff(1:ncoeff) = signs(1:ncoeff) * Coeff(1:ncoeff)

          ! And finally add the projection weights to the projection matrix
          DO j=1,ncoeff 
            val = Coeff(j)

            IF( ABS( val ) < 1.0e-12 ) CYCLE

            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                EdgeCol0 + coeffi(j), EdgeScale * EdgeCoeff * val )
          END DO
        END DO
      END DO
         
      IF( Nundefined > 0 ) THEN
        CALL Error('LevelProjector',&
            'Number of edges could not be mapped: '//TRIM(I2S(Nundefined)))          
      END IF

      WRITE( Message,'(A,ES12.5)') 'Minimum absolute sum of edge weights: ',minwsum
      CALL Info('LevelProjector',Message,Level=10)
      
      WRITE( Message,'(A,ES12.5)') 'Maximum absolute sum of edge weights: ',maxwsum
      CALL Info('LevelProjector',Message,Level=10)
      
      IF( NoSkewed > 0 ) THEN
        CALL Info('LevelProjector','Number of skewed edge mappings: '//TRIM(I2S(NoSkewed)),Level=8)
      END IF
      CALL Info('LevelProjector','Created strong constraints for edge dofs',Level=8)      

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, &
          NodesM % x, NodesM % y, NodesM % z )

    END SUBROUTINE AddEdgeProjectorStrongStrides
    !----------------------------------------------------------------------


    !----------------------------------------------------------------------
    ! Create weak projector for the remaining nodes and edges.
    ! This uses the generic way to introduce the weights. The resulting 
    ! matrix is more dense but should be numerically favourable. 
    ! The integration is done by making an on-the-fly triangularization 
    ! into several triangles. This is not generic - it assumes constant
    ! y levels, and cartesian mesh where the search is done.  
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeakStrides()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: j1,j2,j3,j4,jj,ii,sgn0,k,kmax,ind,indM,nip,nn,ne,nf,inds(10),Ninteg,NintegGen
      TYPE(Element_t), POINTER :: Element, ElementM
      TYPE(Element_t) :: ElementT
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: RightSplit, LeftSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: x(10),y(10),xt,yt,zt,xmax,ymax,xmin,ymin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,ArcTol,u,v,w,um,vm,wm,val,Overlap,RefArea,dArea,&
          SumOverlap,SumArea,qleft, qright, qleft2, qright2, MaxErr,Err,phi(10)
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:),WBasisM(:,:),RotWbasis(:,:),dBasisdx(:,:)
      LOGICAL :: LeftCircle, Stat
      TYPE(Mesh_t), POINTER :: Mesh

      CALL Info('LevelProjector','Creating weak projector for stride mesh',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      ALLOCATE( NodesT % x(n), NodesT % y(n), NodesT % z(n) )
      ALLOCATE( Basis(n), BasisM(n) )
      ALLOCATE( dBasisdx(n,3), WBasis(n,3), WBasisM(n,3), RotWBasis(n,3) )

      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp

      MaxErr = 0.0_dp
      zt = 0.0_dp
      n = 4
      LeftCircle = .FALSE.

      ArcTol = ArcCoeff * Xtol     
      Ninteg = 0
      NintegGen = 0

      ! The temporal triangle used in the numerical integration
      ElementT % TYPE => GetElementType( 303, .FALSE. )
      ElementT % NodeIndexes => IndexesT

      DO ind=1,BMesh1 % NumberOfBulkElements

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        ne = Element % TYPE % NumberOfEdges
        IF( PiolaVersion ) THEN
          nf = 2
        ELSE
          nf = 0
        END IF
        
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))

        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))
        ymin = MINVAL(Nodes % y(1:n))
        ymax = MAXVAL(Nodes % y(1:n))

        IF( Repeating ) THEN
          Nrange = FLOOR( (xmin-XMinAll) / XRange )
          xmin = xmin - Nrange * XRange
          xmax = xmax - Nrange * XRange
          Nodes % x(1:n) = Nodes % x(1:n) - NRange * XRange 
          IF( xmax > XMaxAll ) THEN
            Nrange2 = 1
          ELSE IF( xmax < XMinAll ) THEN
            Nrange2 = -1
          ELSE
            Nrange2 = 0
          END IF
        ELSE IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Nodes % x(1:n) ) > 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Nodes % x(j) < 0.0 ) Nodes % x(j) = Nodes % x(j) + 360.0_dp
            END DO
          END IF
        END IF

        ! Transform the angle to archlength in order to have correct mapping 
        ! of skewed edges.
        Nodes % x(1:n) = ArcCoeff * Nodes % x(1:n)
        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))

        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        IP = GaussPoints( Element ) 
        RefArea = detJ * SUM( IP % s(1:IP % n) )

        SumArea = 0.0_dp
        SumOverlap = 0.0_dp
        
200     sgn0 = 1
        IF( AntiRepeating ) THEN
          IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
        END IF

        ! find an index offset such that [j1,j2,j3,j4] is ordered the as the standard
        ! nodes in bilinear elements. This could be made generic as well, but it was
        ! easier for me to fix these indexes in this way and I was feeling lazy. 
        j1 = 1; j2 = 1; j3 = 1; j4 = 1
        DO j=2,4
          ! Lower left
          IF( Nodes % x(j) + Nodes % y(j) < Nodes % x(j1) + Nodes % y(j1) ) j1 = j
          ! Lower right
          IF( Nodes % x(j) - Nodes % y(j) > Nodes % x(j2) - Nodes % y(j2) ) j2 = j
          ! Upper right
          IF( Nodes % x(j) + Nodes % y(j) > Nodes % x(j3) + Nodes % y(j3) ) j3 = j
          ! Upper left
          IF( Nodes % x(j) - Nodes % y(j) < Nodes % x(j4) - Nodes % y(j4) ) j4 = j
        END DO

        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        DO indM=1,BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes
          
          NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))
          
          ! Make the quick and dirty search first
          yminm = MINVAL( NodesM % y(1:n))
          IF( ABS( ymin - yminm ) > YTol ) CYCLE
          
          ymaxm = MAXVAL( NodesM % y(1:n))
          IF( ABS( ymax - ymaxm ) > YTol ) CYCLE
          
          NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))

          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definately on the right circle
            IF( ALL( ABS( NodesM % x(1:n) ) - 90.0 < Xtol ) ) CYCLE
            DO j=1,n
              IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          ! Transfer into real length units instead of angles
          ! This gives right balance between x and y -directions. 
          NodesM % x(1:n) = ArcCoeff * NodesM % x(1:n)

          xminm = MINVAL( NodesM % x(1:n))
          xmaxm = MAXVAL( NodesM % x(1:n))
                    
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0 ) CYCLE
          END IF
          
          Overlap = (MIN(xmax, xmaxm)- MAX(xmin,xminm))/(xmax-xmin)
          IF( Overlap < RelTolX ) CYCLE 
          
          SumOverlap = SumOverlap + Overlap
          Ninteg = Ninteg + 1
          
          ! Then if this is a possible element create a list of the corner nodes
          ! for a temporal mesh. There will be 3 to 6 corner nodes. 
          ! Check the crossings between the edges of the quadrilaters. These will
          ! be used as new points when creating the virtual triangle mesh. 
          LeftSplit = ( ( Nodes % x(j1) - xminm ) * ( xminm - Nodes % x(j4) ) > 0.0_dp )
          IF(LeftSplit) qleft =  ( Nodes % x(j1) - xminm ) / ( Nodes % x(j1) - Nodes % x(j4) )

          RightSplit = ( ( Nodes % x(j2) - xmaxm ) * ( xmaxm - Nodes % x(j3) ) > 0.0_dp )
          IF(RightSplit) qright = ( Nodes % x(j2) - xmaxm ) / ( Nodes % x(j2) - Nodes % x(j3) )

          LeftSplit2 = ( ( Nodes % x(j2) - xminm ) * ( xminm - Nodes % x(j3) ) > 0.0_dp )
          IF(LeftSplit2) qleft2 =  ( Nodes % x(j2) - xminm ) / ( Nodes % x(j2) - Nodes % x(j3) )

          RightSplit2 = ( ( Nodes % x(j1) - xmaxm ) * ( xmaxm - Nodes % x(j4) ) > 0.0_dp )
          IF(RightSplit2) qright2 = ( Nodes % x(j1) - xmaxm ) / ( Nodes % x(j1) - Nodes % x(j4) )

            ! Mark the splits on the vertical edges aligned with the y-axis
            k = 0
            IF( LeftSplit ) THEN
              k = k + 1
              x(k) = xminm
              qleft = MAX( 0.0, MIN( 1.0, qleft ) )
              y(k) = Nodes % y(j1) + qleft * ( Nodes % y(j4) - Nodes % y(j1))
            END IF
            IF( RightSplit2 ) THEN
              k = k + 1
              x(k) = xmaxm
              qright2 = MAX( 0.0, MIN( 1.0, qright2 ) )
              y(k) = Nodes % y(j1) + qright2 * ( Nodes % y(j4) - Nodes % y(j1))
            END IF
            IF( RightSplit ) THEN
              k = k + 1
              x(k) = xmaxm
              qright = MAX( 0.0, MIN( 1.0, qright ) )
              y(k) = Nodes % y(j2) + qright * ( Nodes % y(j3) - Nodes % y(j2))
            END IF
            IF( LeftSplit2 ) THEN
              k = k + 1
              x(k) = xminm
              qleft2 = MAX( 0.0, MIN( 1.0, qleft2 ) )
              y(k) = Nodes % y(j2) + qleft2 * ( Nodes % y(j3) - Nodes % y(j2))
            END IF

            ! Mark the splits on the horizontal axis
            BottomEdge = .NOT. ( ( Nodes % x(j2) < xminm ) .OR. ( Nodes % x(j1) > xmaxm ) )
            TopEdge    = .NOT. ( ( Nodes % x(j3) < xminm ) .OR. ( Nodes % x(j4) > xmaxm ) )

            IF( BottomEdge ) THEN
              k = k + 1
              x(k) = MAX( xminm, Nodes % x(j1) )
              y(k) = yminm
              k = k + 1
              x(k) = MIN( xmaxm, Nodes % x(j2) )
              y(k) = yminm
            END IF
            IF( TopEdge ) THEN
              k = k + 1
              x(k) = MIN( xmaxm, Nodes % x(j3) )
              y(k) = ymaxm
              k = k + 1
              x(k) = MAX( xminm, Nodes % x(j4) )
              y(k) = ymaxm
            END IF
            kmax = k 

            IF( kmax < 3 ) THEN
              CALL Warn('LevelProjector','Cannot integrate over '//TRIM(I2S(kmax))//' nodes')
              CYCLE
            END IF
            
            ! The polygon is convex and hence its center lies inside the polygon
            xt = SUM(x(1:kmax)) / kmax
            yt = SUM(y(1:kmax)) / kmax

            ! Set the angle from the center and order the nodes so that they 
            ! can be easily triangulated.
            DO k=1,kmax
              phi(k) = ATAN2( y(k)-yt, x(k)-xt )
              inds(k) = k
            END DO
                       
            CALL SortR(kmax,inds,phi)
            x(1:kmax) = x(inds(1:kmax))
            y(1:kmax) = y(inds(1:kmax))

            !PRINT *,'Polygon: ',ind,indm,LeftSplit, RightSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge, kmax 

          ! Deal the case with multiple corners by making 
          ! triangulariation using one corner point.
          ! This should be ok as the polygon is always convex.
          NodesT % x(1) = x(1)
          NodesT % y(1) = y(1)

          ! Use somewhat higher integration rules than the default
          IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2 ) 

          DO k=1,kmax-2                         

            ! This check over area also automatically elimiates redundant nodes
            ! that were detected twice.
            dArea = 0.5*ABS( (x(k+1)-x(1))*(y(k+2)-y(1)) -(x(k+2)-x(1))*(y(k+1)-y(1)))
            IF( dArea < RelTolY**2 * RefArea ) CYCLE

            NodesT % x(2) = x(k+1)
            NodesT % y(2) = y(k+1)
            NodesT % x(3) = x(k+2)
            NodesT % y(3) = y(k+2)
            
            ! Integration over the temporal element
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),IP % v(nip),IP % w(nip),detJ,Basis)
              
              ! We will actually only use the global coordinates and the integration weight 
              ! from the temporal mesh. 

              ! Global coordinates of the integration point
              xt = SUM( Basis(1:3) * NodesT % x(1:3) )
              yt = SUM( Basis(1:3) * NodesT % y(1:3) )
              zt = 0.0_dp

              ! Integration weight for current integration point
              Wtemp = DetJ * IP % s(nip)
              sumarea = sumarea + Wtemp
              
              ! Integration point at the slave element
              CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx,EdgeBasis=WBasis)
                ELSE
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx )
                  CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
              END IF

              ! Integration point at the master element
              CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, Basis, dBasisdx, EdgeBasis=WBasisM)
                ELSE
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, BasisM, dBasisdx )
                  CALL GetEdgeBasis(ElementM,WBasisM,RotWBasis,BasisM,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
              END IF

              ! Add the nodal dofs
              IF( DoNodes .AND. .NOT. StrongNodes ) THEN
                DO j=1,n 
                  jj = Indexes(j)                                    
                  nrow = NodePerm( InvPerm1(jj) )
                  IF( nrow == 0 ) CYCLE

                  Projector % InvPerm(nrow) = InvPerm1(jj)
                  val = Basis(j) * Wtemp
                  DO i=1,n
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val ) 

                    IF( ABS( val * BasisM(i) ) < 1.0e-10 ) CYCLE
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm2(IndexesM(i)), -NodeScale * NodeCoeff * BasisM(i) * val )   
                  END DO
                END DO
              END IF

              IF( DoEdges ) THEN
                ! Dofs are numbered as follows:
                ! 1....number of nodes
                ! + ( 1 ... number of edges )
                ! + ( 1 ... 2 x number of faces )
                !-------------------------------------------
                DO j=1,ne+nf
                  
                  IF( j <= ne ) THEN
                    jj = Element % EdgeIndexes(j) 
                    IF( EdgePerm(jj) == 0 ) CYCLE
                    nrow = EdgeRow0 + EdgePerm(jj)
                    jj = jj + EdgeCol0
                    Projector % InvPerm( nrow ) = jj
                  ELSE
                    jj = 2 * ( ind - 1 ) + ( j - 4 )
                    nrow = FaceRow0 + jj
                    jj = 2 * ( Element % ElementIndex - 1) + ( j - 4 ) 
                    Projector % InvPerm( nrow ) = FaceCol0 + jj
                  END IF
                                   
                  DO i=1,ne+nf
                    IF( i <= ne ) THEN
                      ii = Element % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( Element % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF
                    val = Wtemp * SUM( WBasis(j,:) * Wbasis(i,:) ) 
                    IF( ABS( val ) > 1.0e-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeCoeff * val ) 
                    END IF

                    IF( i <= ne ) THEN
                      ii = ElementM % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( ElementM % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF                    
                    val = -Wtemp * SUM( WBasis(j,:) * WBasisM(i,:) ) 
                    IF( ABS( val ) > 1.0e-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeScale * EdgeCoeff * val  ) 
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO
        END DO
        
        IF( Repeating ) THEN
          IF( NRange2 /= 0 ) THEN
            xmin = xmin - ArcCoeff * Nrange2 * XRange
            xmax = xmax - ArcCoeff * Nrange2 * XRange
            Nodes % x(1:n) = Nodes % x(1:n) - ArcCoeff * NRange2 * XRange 
            NRange = NRange + NRange2
            NRange2 = 0
            GOTO 200
          END IF
        END IF

        Err = SumArea/RefArea-1.0_dp
        MaxErr = MAX( MaxErr,ABS(Err))
      END DO

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )
      DEALLOCATE( NodesT % x, NodesT % y, NodesT % z )
      DEALLOCATE( Basis, BasisM )
      DEALLOCATE( dBasisdx, WBasis, WBasisM, RotWBasis )

      CALL Info('LevelProjector','Number of integration pairs: '&
          //TRIM(I2S(Ninteg)),Level=10)

      WRITE( Message,'(A,ES12.3)') 'Maximum error in area integration:',MaxErr 
      CALL Info('LevelProjector',Message,Level=8)


    END SUBROUTINE AddProjectorWeakStrides



    !----------------------------------------------------------------------
    ! Create weak projector for the remaining nodes and edges
    ! using generic algo that can deal with triangles and quadrilaterals.
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeakGeneric()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: jj,ii,sgn0,k,kmax,ind,indM,nip,nn,ne,nf,inds(10),nM,neM,nfM,iM,i2,i2M
      INTEGER :: ElemCands, TotCands, ElemHits, TotHits, EdgeHits, CornerHits, &
          MaxErrInd, MinErrInd, InitialHits, ActiveHits, TimeStep
      TYPE(Element_t), POINTER :: Element, ElementM
      TYPE(Element_t) :: ElementT
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: RightSplit, LeftSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: x(10),y(10),xt,yt,zt,xmax,ymax,xmin,ymin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,ArcTol,u,v,w,um,vm,wm,val,RefArea,dArea,&
          SumArea,TrueArea,MaxErr,MinErr,Err,phi(10),Point(3),uvw(3),ArcRange , &
          val_dual, zmin, zmax, zminm, zmaxm
      REAL(KIND=dp) :: A(2,2), B(2), C(2), absA, detA, rlen, &
          x1, x2, y1, y2, x1M, x2M, y1M, y2M, x0, y0, dist, DistTol
      REAL(KIND=dp) :: TotRefArea, TotSumArea, TotTrueArea
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:),WBasisM(:,:),RotWbasis(:,:),dBasisdx(:,:)
      LOGICAL :: LeftCircle, Stat, CornerFound(4), CornerFoundM(4)
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: TimestepVar

      ! These are used temporarely for debugging purposes
      INTEGER :: SaveInd, MaxSubElem, MaxSubTriangles, DebugInd, Nslave, Nmaster
      LOGICAL :: SaveElem, DebugElem
      CHARACTER(LEN=20) :: FileName

      REAL(KIND=dp) :: Area
      REAL(KIND=dp), ALLOCATABLE :: CoeffBasis(:), MASS(:,:)

      CALL Info('LevelProjector','Creating weak constraints using a generic integrator',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      SaveInd = ListGetInteger( BC,'Level Projector Save Element Index',Found )
      DebugInd = ListGetInteger( BC,'Level Projector Debug Element Index',Found )

      TimestepVar => VariableGet( Mesh % Variables,'Timestep',ThisOnly=.TRUE. )
      Timestep = NINT( TimestepVar % Values(1) )
 
      n = Mesh % MaxElementNodes

      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      ALLOCATE( NodesT % x(n), NodesT % y(n), NodesT % z(n) )
      ALLOCATE( Basis(n), BasisM(n) )
      ALLOCATE( dBasisdx(n,3), WBasis(n,3), WBasisM(n,3), RotWBasis(n,3) )

      IF(BiOrthogonalBasis) ALLOCATE(CoeffBasis(n), MASS(n,n))

      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp

      MaxErr = 0.0_dp
      MinErr = HUGE( MinErr )
      MaxErrInd = 0
      MinErrInd = 0
      zt = 0.0_dp
      LeftCircle = .FALSE.
      ArcTol = ArcCoeff * Xtol
      ArcRange = ArcCoeff * Xrange 
     
      DistTol = ArcTol**2 + YTol**2

      ! The temporal triangle used in the numerical integration
      ElementT % TYPE => GetElementType( 303, .FALSE. )
      ElementT % NodeIndexes => IndexesT
      TotCands = 0
      TotHits = 0
      EdgeHits = 0
      CornerHits = 0
      InitialHits = 0
      ActiveHits = 0
      TotRefArea = 0.0_dp
      TotSumArea = 0.0_dp
      TotTrueArea = 0.0_dp
      Point = 0.0_dp
      MaxSubTriangles = 0
      Nslave = 0
      Nmaster = 0


      DO ind=1,BMesh1 % NumberOfBulkElements

        ! Optionally save the submesh for specified element, for vizualization and debugging
        SaveElem = ( SaveInd == ind )
        DebugElem = ( DebugInd == ind )

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        ne = Element % TYPE % NumberOfEdges
        IF( PiolaVersion .AND. ne == 4 ) THEN
          nf = 2
        ELSE
          nf = 0
        END IF
        
        ! Transform the angle to archlength in order to have correct balance between x and y
        Nodes % x(1:n) = ArcCoeff * BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))

        IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Nodes % x(1:n) ) > ArcCoeff * 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Nodes % x(j) < 0.0 ) Nodes % x(j) = &
                  Nodes % x(j) + ArcCoeff * 360.0_dp
            END DO
          END IF
        END IF

        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))

        ymin = MINVAL(Nodes % y(1:n))
        ymax = MAXVAL(Nodes % y(1:n))
                
        IF( HaveMaxDistance ) THEN
          zmin = MINVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
          zmax = MAXVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
        END IF

        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        IP = GaussPoints( Element ) 
        RefArea = detJ * SUM( IP % s(1:IP % n) )
        SumArea = 0.0_dp
        TrueArea = 0.0_dp

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_a.dat'
          OPEN( 10,FILE=Filename)
          DO i=1,n
            WRITE( 10, * ) Nodes % x(i), Nodes % y(i)
          END DO
          CLOSE( 10 )
        END IF
        
        IF( DebugElem ) THEN
          PRINT *,'Debug Elem:',ind,n,LeftCircle
          PRINT *,'ArcTol:',ArcTol
          PRINT *,'X:',Nodes % x(1:n)
          PRINT *,'Y:',Nodes % y(1:n)
          PRINT *,'Xrange:',Xmin,Xmax
          PRINT *,'Yrange:',Ymin,Ymax
        END IF


        IF( DoNodes .AND. .NOT. StrongNodes ) THEN
          DO i=1,n
            j = InvPerm1(Indexes(i))
            nrow = NodePerm(j)
            IF( nrow == 0 ) CYCLE
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                j, 0.0_dp ) 
             IF(ASSOCIATED(Projector % Child)) &
               CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                   j, 0.0_dp ) 
          END DO
        END IF


        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        ElemCands = 0
        ElemHits = 0
        DO indM=1,BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes
          nM = ElementM % TYPE % NumberOfNodes

          IF( HaveMaxDistance ) THEN
            zminm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            zmaxm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            IF( zmaxm < zmin - MaxDistance ) CYCLE
            IF( zminm > zmax + MaxDistance ) CYCLE
          END IF
          
          NodesM % y(1:nM) = BMesh2 % Nodes % y(IndexesM(1:nM))
          
          ! Make the quick and dirty search first
          ! This requires some minimal width of the cut
          yminm = MINVAL( NodesM % y(1:nM))
          IF( yminm > ymax - Ytol ) CYCLE
          
          ymaxm = MAXVAL( NodesM % y(1:nM))
          IF( ymaxm < ymin + Ytol ) CYCLE
          
          NodesM % x(1:nM) = ArcCoeff * BMesh2 % Nodes % x(IndexesM(1:nM))
          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definately on the right circle
            IF( ALL( ABS( NodesM % x(1:nM) ) - ArcCoeff * 90.0 < ArcTol ) ) CYCLE
            DO j=1,nM
              IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = &
                  NodesM % x(j) + ArcCoeff * 360.0_dp
            END DO
          END IF
          
          xminm = MINVAL( NodesM % x(1:nM))
          xmaxm = MAXVAL( NodesM % x(1:nM))

          IF( Repeating ) THEN
            ! Enforce xmaxm to be on the same interval than xmin
            Nrange = FLOOR( (xmaxm-xmin+ArcTol) / ArcRange )
            IF( Nrange /= 0 ) THEN
              xminm = xminm - Nrange * ArcRange
              xmaxm = xmaxm - Nrange * ArcRange
              NodesM % x(1:nM) = NodesM % x(1:nM) - NRange * ArcRange 
            END IF

            ! Check whether there could be a intersection in an other interval as well
            IF( xminm + ArcRange < xmax + ArcTol ) THEN
              Nrange2 = 1
            ELSE
              Nrange2 = 0
            END IF
          END IF

          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0 ) CYCLE
          END IF

200       IF( xminm > xmax - ArcTol ) GOTO 100
          IF( xmaxm < xmin + ArcTol ) GOTO 100

          IF( DebugElem ) THEN
            PRINT *,'Candidate Elem:',indM,nM
            PRINT *,'X:',NodesM % x(1:nM)
            PRINT *,'Y:',NodesM % y(1:nM)
          END IF

          neM = ElementM % TYPE % NumberOfEdges
          IF( PiolaVersion .AND. neM == 4 ) THEN
            nfM = 2
          ELSE
            nfM = 0
          END IF
          
          k = 0
          ElemCands = ElemCands + 1
          CornerFound = .FALSE.
          CornerFoundM = .FALSE.

          ! Check through the nodes that are created in the intersections of any two edge
          DO i=1,ne
            x1 = Nodes % x(i)
            y1 = Nodes % y(i)
            i2 = i + 1 
            IF( i2 > ne ) i2 = 1
            x2 = Nodes % x(i2)
            y2 = Nodes % y(i2)

            DO iM=1,neM
              x1M = NodesM % x(iM)
              y1M = NodesM % y(iM)
              i2M = iM + 1
              IF( i2M > neM ) i2M = 1
              x2M = NodesM % x(i2M)
              y2M = NodesM % y(i2M)
              
              ! Upon solution this is tampered so it must be initialized 
              ! before each solution. 
              A(1,1) = x2 - x1
              A(2,1) = y2 - y1           
              A(1,2) = x1M - x2M
              A(2,2) = y1M - y2M

              detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
              absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))
              
              ! Lines are almost parallel => no intersection possible
              ! Check the hist at the end of the line segments.
              IF(ABS(detA) < 1.0d-8 * absA + 1.0d-20 ) CYCLE

              B(1) = x1M - x1
              B(2) = y1M - y1
              
              CALL InvertMatrix( A,2 )
              C(1:2) = MATMUL(A(1:2,1:2),B(1:2))
              
              IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) CYCLE
              
              ! We have a hit, two line segments can have only one hit
              k = k + 1
              
              x(k) = x1 + C(1) * (x2-x1)
              y(k) = y1 + C(1) * (y2-y1)

              ! If the point of intersection is at the end of a line-segment it
              ! is also a corner node.
              IF(ABS(C(1)) < 1.0d-6 ) THEN
                CornerFound(i) = .TRUE.
              ELSE IF( ABS(C(1)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFound(i2) = .TRUE.
              END IF              

              IF(ABS(C(2)) < 1.0d-6 ) THEN
                CornerFoundM(iM) = .TRUE.
              ELSE IF( ABS(C(2)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFoundM(i2M) = .TRUE.
              END IF
         
              EdgeHits = EdgeHits + 1
            END DO
          END DO

          ! Check the nodes that are one of the existing nodes i.e. corner nodes
          ! that are located inside in either element. We have to check both combinations. 
          DO i=1,n
            ! This corner was already determined active as the end of edge 
            IF( CornerFound(i) ) CYCLE

            Point(1) = Nodes % x(i)
            IF( Point(1) < xminm - ArcTol ) CYCLE
            IF( Point(1) > xmaxm + ArcTol ) CYCLE

            Point(2) = Nodes % y(i)
            IF( Point(2) < yminm - YTol ) CYCLE
            IF( Point(2) > ymaxm + YTol ) CYCLE

            ! The edge intersections should catch the sharp hits so here we can use hard criteria
            Found = PointInElement( ElementM, NodesM, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

          DO i=1,nM
            IF( CornerFoundM(i) ) CYCLE

            Point(1) = NodesM % x(i)
            IF( Point(1) < xmin - ArcTol ) CYCLE
            IF( Point(1) > xmax + ArcTol ) CYCLE

            Point(2) = NodesM % y(i)
            IF( Point(2) < ymin - YTol ) CYCLE
            IF( Point(2) > ymax + YTol ) CYCLE
         
            Found = PointInElement( Element, Nodes, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

          kmax = k          
          IF( kmax < 3 ) GOTO 100

          sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
          
          InitialHits = InitialHits + kmax

          ! The polygon is convex and hence its center lies inside the polygon
          xt = SUM(x(1:kmax)) / kmax
          yt = SUM(y(1:kmax)) / kmax
          
          ! Set the angle from the center and order the nodes so that they 
          ! can be easily triangulated.
          DO k=1,kmax
            phi(k) = ATAN2( y(k)-yt, x(k)-xt )
            inds(k) = k
          END DO

          CALL SortR(kmax,inds,phi)
          x(1:kmax) = x(inds(1:kmax))
          y(1:kmax) = y(inds(1:kmax))

          ! Eliminate redundant corners from the polygon
          j = 1
          DO k=2,kmax
            dist = (x(j)-x(k))**2 + (y(j)-y(k))**2 
            IF( dist > DistTol ) THEN
              j = j + 1
              IF( j /= k ) THEN
                x(j) = x(k)
                y(j) = y(k)
              END IF
            END IF
          END DO
          kmax = j

          IF( kmax < 3 ) GOTO 100

          IF( DebugElem ) THEN
            PRINT *,'Corners:',kmax
            PRINT *,'Center:',xt,yt
          END IF

          ElemHits = ElemHits + 1
          ActiveHits = ActiveHits + kmax

          IF( kmax > MaxSubTriangles ) THEN
            MaxSubTriangles = kmax
            MaxSubElem = ind
          END IF

          IF( SaveElem ) THEN
            FileName = 't'//TRIM(I2S(TimeStep))//'_b'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) NodesM % x(i), NodesM % y(i)
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_d'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) xt, yt
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_e'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,kmax
              WRITE( 10, * ) x(i), y(i)
            END DO
            CLOSE( 10 )           
          END IF

          
          ! Deal the case with multiple corners by making 
          ! triangulariation using one corner point.
          ! This should be ok as the polygon is always convex.
          NodesT % x(1) = x(1)
          NodesT % y(1) = y(1)
          
          ! Use somewhat higher integration rules than the default
          IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2 ) 
          
          DO k=1,kmax-2                         
            
            ! This check over area also automatically elimiates redundant nodes
            ! that were detected twice.
            dArea = 0.5*ABS( (x(k+1)-x(1))*(y(k+2)-y(1)) -(x(k+2)-x(1))*(y(k+1)-y(1)))

            IF( DebugElem ) PRINT *,'dArea:',dArea,refArea

            IF( dArea < RelTolY**2 * RefArea ) CYCLE
            
            NodesT % x(2) = x(k+1)
            NodesT % y(2) = y(k+1)
            NodesT % x(3) = x(k+2)
            NodesT % y(3) = y(k+2)

            IF(BiOrthogonalBasis) THEN
              MASS  = 0
              CoeffBasis = 0
              area = 0._dp
              DO nip=1, IP % n 
                stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                    IP % v(nip),IP % w(nip),detJ,Basis)
                IF(.NOT. Stat ) EXIT

                ! We will actually only use the global coordinates and the integration weight 
                ! from the temporal mesh. 
              
                ! Global coordinates of the integration point
                xt = SUM( Basis(1:3) * NodesT % x(1:3) )
                yt = SUM( Basis(1:3) * NodesT % y(1:3) )
                zt = 0.0_dp
              
                ! Integration weight for current integration point
                Wtemp = DetJ * IP % s(nip)
                area = area + wtemp
              
                ! Integration point at the slave element
                CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
                IF(.NOT. Stat) CYCLE

                DO i=1,n
                  DO j=1,n
                    MASS(i,j) = MASS(i,j) + wTemp * Basis(i) * Basis(j)
                  END DO
                  CoeffBasis(i) = CoeffBasis(i) + wTemp * Basis(i)
                END DO
              END DO

              IF(Area<1.d-12) GOTO 100

              CALL InvertMatrix( MASS, n )

              DO i=1,n
                DO j=1,n
                  MASS(i,j) = MASS(i,j) * CoeffBasis(i)
                END DO
              END DO
            END IF
            
            ! Integration over the temporal element
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                  IP % v(nip),IP % w(nip),detJ,Basis)
              IF(.NOT. Stat) EXIT

              ! We will actually only use the global coordinates and the integration weight 
              ! from the temporal mesh. 
              
              ! Global coordinates of the integration point
              xt = SUM( Basis(1:3) * NodesT % x(1:3) )
              yt = SUM( Basis(1:3) * NodesT % y(1:3) )
              zt = 0.0_dp
              
              ! Integration weight for current integration point
              Wtemp = DetJ * IP % s(nip)
              sumarea = sumarea + Wtemp
              
              ! Integration point at the slave element
              CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx,EdgeBasis=WBasis)
                ELSE
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx )
                  CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
              END IF

              ! Integration point at the master element
              CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, Basis, dBasisdx, EdgeBasis=WBasisM)
                ELSE
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, BasisM, dBasisdx )
                  CALL GetEdgeBasis(ElementM,WBasisM,RotWBasis,BasisM,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
              END IF
              IF(.NOT. Stat) CYCLE

              ! Add the nodal dofs
              IF( DoNodes .AND. .NOT. StrongNodes ) THEN
                IF(BiOrthogonalBasis) THEN
                  CoeffBasis = 0._dp
                  DO i=1,n
                    DO j=1,n
                      CoeffBasis(i) = CoeffBasis(i) + MASS(i,j) * Basis(j)
                    END DO
                  END DO
                END IF

                DO j=1,n 
                  jj = Indexes(j)                                    

                  nrow = NodePerm(InvPerm1(jj))
                  IF( nrow == 0 ) CYCLE

                  Projector % InvPerm(nrow) = InvPerm1(jj)
                  val = Basis(j) * Wtemp
                  IF(BiorthogonalBasis) val_dual = CoeffBasis(j) * Wtemp

                  TrueArea = TrueArea + val

                  IF( DebugElem ) PRINT *,'Vals:',val

                  DO i=1,n
                    Nslave = Nslave + 1
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val ) 

                    IF(BiOrthogonalBasis) THEN
                      CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                            InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val_dual ) 
                    END IF
                  END DO

                  DO i=1,nM
                    IF( ABS( val * BasisM(i) ) < 1.0e-10 ) CYCLE

                    Nmaster = Nmaster + 1
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )                   

                    IF(BiOrthogonalBasis) THEN
                      IF(DualMaster.OR.DualLCoeff) THEN
                        CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                              InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val_dual ) 
                      ELSE
                        CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                              InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val ) 
                      END IF
                    END IF
                  END DO
                END DO
              END IF

              IF( DoEdges ) THEN
                ! Dofs are numbered as follows:
                ! 1....number of nodes
                ! + ( 1 ... number of edges )
                ! + ( 1 ... 2 x number of faces )
                !-------------------------------------------
                DO j=1,ne+nf
                  
                  IF( j <= ne ) THEN
                    jj = Element % EdgeIndexes(j) 
                    IF( EdgePerm(jj) == 0 ) CYCLE
                    nrow = EdgeRow0 + EdgePerm(jj)
                    jj = jj + EdgeCol0
                    Projector % InvPerm( nrow ) = jj
                  ELSE
                    IF( Parallel ) THEN
                      IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
                    END IF

                    jj = 2 * ( ind - 1 ) + ( j - 4 )
                    nrow = FaceRow0 + jj
                    jj = 2 * ( Element % ElementIndex - 1) + ( j - 4 ) 
                    Projector % InvPerm( nrow ) = FaceCol0 + jj
                  END IF
                                   
                  DO i=1,neM+nfM
                    IF( i <= neM ) THEN
                      ii = Element % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( Element % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF

                    val = Wtemp * SUM( WBasis(j,:) * Wbasis(i,:) ) 
                    IF( ABS( val ) > 1.0e-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeCoeff * val ) 
                    END IF

                    IF( i <= neM ) THEN
                      ii = ElementM % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( ElementM % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF                    
                    val = -Wtemp * sgn0 * SUM( WBasis(j,:) * WBasisM(i,:) ) 
                    IF( ABS( val ) > 1.0e-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeScale * EdgeCoeff * val  ) 
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO

100       IF( Repeating ) THEN
            IF( NRange2 /= 0 ) THEN
              xminm = xminm + ArcCoeff * Nrange2 * ArcRange
              xmaxm = xmaxm + ArcCoeff * Nrange2 * ArcRange
              NodesM % x(1:n) = NodesM % x(1:n) + NRange2 * ArcRange 
              NRange = NRange + NRange2
              NRange2 = 0
              GOTO 200
            END IF
          END IF

        END DO

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_n.dat'
          OPEN( 10,FILE=Filename)
          OPEN( 10,FILE=FileName)
          WRITE( 10, * ) ElemHits 
          CLOSE( 10 )
        END IF
        
        TotCands = TotCands + ElemCands
        TotHits = TotHits + ElemHits
        TotSumArea = TotSumArea + SumArea
        TotRefArea = TotRefArea + RefArea
        TotTrueArea = TotTruearea + TrueArea

        Err = SumArea / RefArea
        IF( Err > MaxErr ) THEN
          MaxErr = Err
          MaxErrInd = Err
        END IF
        IF( Err < MinErr ) THEN
          MinErr = Err
          MinErrInd = ind
        END IF
      END DO

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )
      DEALLOCATE( NodesT % x, NodesT % y, NodesT % z )
      DEALLOCATE( Basis, BasisM )
      DEALLOCATE( dBasisdx, WBasis, WBasisM, RotWBasis )

      CALL Info('LevelProjector','Number of integration pair candidates: '&
          //TRIM(I2S(TotCands)),Level=10)
      CALL Info('LevelProjector','Number of integration pairs: '&
          //TRIM(I2S(TotHits)),Level=10)

      CALL Info('LevelProjector','Number of edge intersections: '&
          //TRIM(I2S(EdgeHits)),Level=10)
      CALL Info('LevelProjector','Number of corners inside element: '&
          //TRIM(I2S(EdgeHits)),Level=10)

      CALL Info('LevelProjector','Number of initial corners: '&
          //TRIM(I2S(InitialHits)),Level=10)
      CALL Info('LevelProjector','Number of active corners: '&
          //TRIM(I2S(ActiveHits)),Level=10)

      CALL Info('LevelProjector','Number of most subelement corners: '&
          //TRIM(I2S(MaxSubTriangles)),Level=10)
      CALL Info('LevelProjector','Element of most subelement corners: '&
          //TRIM(I2S(MaxSubElem)),Level=10)

      WRITE( Message,'(A,ES12.5)') 'Total reference area:',TotRefArea
      CALL Info('LevelProjector',Message,Level=8)
      WRITE( Message,'(A,ES12.5)') 'Total integrated area:',TotSumArea
      CALL Info('LevelProjector',Message,Level=8)

      Err = TotSumArea / TotRefArea
      WRITE( Message,'(A,ES12.3)') 'Average ratio in area integration:',Err 
      CALL Info('LevelProjector',Message,Level=8)

      WRITE( Message,'(A,ES12.5)') 'True integrated area:',TotTrueArea
      CALL Info('LevelProjector',Message,Level=8)

      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Maximum relative discrepancy in areas (element: ',MaxErrInd,'):',MaxErr-1.0_dp 
      CALL Info('LevelProjector',Message,Level=8)
      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Minimum relative discrepancy in areas (element: ',MinErrInd,'):',MinErr-1.0_dp 
      CALL Info('LevelProjector',Message,Level=8)

      CALL Info('LevelProjector','Number of slave entries: '&
          //TRIM(I2S(Nslave)),Level=10)
      CALL Info('LevelProjector','Number of master entries: '&
          //TRIM(I2S(Nmaster)),Level=10)



    END SUBROUTINE AddProjectorWeakGeneric


    !----------------------------------------------------------------------
    ! Create weak projector for the nodes in 1D mesh.
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeak1D()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: jj,ii,sgn0,k,kmax,ind,indM,nip,nn,inds(10),nM,iM,i2,i2M
      INTEGER :: ElemHits, TotHits, MaxErrInd, MinErrInd, TimeStep, AntiPeriodicHits
      TYPE(Element_t), POINTER :: Element, ElementM
      TYPE(Element_t) :: ElementT 
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: xt,yt,zt,xmax,xmin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,u,v,w,um,vm,wm,val,RefArea,dArea,&
          SumArea,MaxErr,MinErr,Err,uvw(3),val_dual,dx,dxcut, &
          zmin,zmax, zminm, zmaxm
      REAL(KIND=dp) :: TotRefArea, TotSumArea
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      LOGICAL :: LeftCircle, Stat
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: TimestepVar

      ! These are used temporarely for debugging purposes
      INTEGER :: SaveInd
      LOGICAL :: SaveElem
      CHARACTER(LEN=20) :: FileName

      REAL(KIND=dp), ALLOCATABLE :: CoeffBasis(:), MASS(:,:)

      CALL Info('LevelProjector','Creating weak constraints using a 1D integrator',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      SaveInd = ListGetInteger( BC,'Level Projector Save Element Index',Found )
      TimestepVar => VariableGet( Mesh % Variables,'Timestep',ThisOnly=.TRUE. )
      Timestep = NINT( TimestepVar % Values(1) )
 
      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      ALLOCATE( NodesT % x(n), NodesT % y(n), NodesT % z(n) )
      ALLOCATE( Basis(n), BasisM(n) )

      IF (BiOrthogonalBasis) ALLOCATE(CoeffBasis(n), MASS(n,n))

      Nodes % y  = 0.0_dp
      NodesM % y = 0.0_dp
      NodesT % y = 0.0_dp
      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp
      yt = 0.0_dp
      zt = 0.0_dp

      MaxErr = 0.0_dp
      MinErr = HUGE( MinErr )
      MaxErrInd = 0
      MinErrInd = 0
      zt = 0.0_dp
      LeftCircle = .FALSE.
     
      ! The temporal element segment used in the numerical integration
      ElementT % TYPE => GetElementType( 202, .FALSE. )
      ElementT % NodeIndexes => IndexesT
      IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2  ) 

      TotHits = 0
      AntiPeriodicHits = 0
      TotRefArea = 0.0_dp
      TotSumArea = 0.0_dp


      DO ind=1,BMesh1 % NumberOfBulkElements

        ! Optionally save the submesh for specified element, for vizualization and debugging
        SaveElem = ( SaveInd == ind )

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes
        
        n = Element % TYPE % NumberOfNodes        
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))

        ! There is a discontinuity of angle at 180 degs
        ! If we are working on left-hand-side then add 360 degs to the negative angles
        ! to remove this discontinuity.
        IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Nodes % x(1:n) ) > 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Nodes % x(j) < 0.0 ) Nodes % x(j) = &
                  Nodes % x(j) + 360.0_dp
            END DO
          END IF
        END IF

        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))
        dx = xmax - xmin 

        ! The flattened dimension is always the z-component
        IF( HaveMaxDistance ) THEN
          zmin = MINVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
          zmax = MAXVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
        END IF
                        
        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        RefArea = detJ * ArcCoeff * SUM( IP % s(1:IP % n) )
        SumArea = 0.0_dp
        
        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_a.dat'
          OPEN( 10,FILE=Filename)
          DO i=1,n
            WRITE( 10, * ) Nodes % x(i)
          END DO
          CLOSE( 10 )
        END IF

        ! Set the values to maintain the size of the matrix
        ! The size of the matrix is used when allocating for utility vectors of contact algo.
        ! This does not set the Projector % InvPerm to nonzero value that is used to 
        ! determine whether there really is a projector. 
        DO i=1,n
          j = InvPerm1(Indexes(i))
          nrow = NodePerm(j)
          IF( nrow == 0 ) CYCLE
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              j, 0.0_dp ) 
        END DO

        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        ElemHits = 0
        DO indM=1,BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes

          nM = ElementM % TYPE % NumberOfNodes

 
          NodesM % x(1:nM) = BMesh2 % Nodes % x(IndexesM(1:nM))

          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definately on the right circle
            IF( ALL( ABS( NodesM % x(1:nM) ) - 90.0 < XTol ) ) CYCLE
            DO j=1,nM
              IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = &
                  NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          xminm = MINVAL( NodesM % x(1:nM))
          xmaxm = MAXVAL( NodesM % x(1:nM))

          IF( Repeating ) THEN
            ! Enforce xmaxm to be on the same interval than xmin
            Nrange = FLOOR( (xmaxm-xmin+XTol) / XRange )
            IF( Nrange /= 0 ) THEN
              xminm = xminm - Nrange * XRange
              xmaxm = xmaxm - Nrange * XRange
              NodesM % x(1:nM) = NodesM % x(1:nM) - NRange * XRange 
            END IF

            ! Check whether there could be a intersection in an other interval as well
            IF( xminm + XRange < xmax + XTol ) THEN
              Nrange2 = 1
            ELSE
              Nrange2 = 0
            END IF
          END IF

          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > 180.0 ) CYCLE
          END IF          

200       IF( xminm >= xmax ) GOTO 100
          IF( xmaxm <= xmin ) GOTO 100

          
          ! This is a cheap test so perform that first, if requested
          IF( HaveMaxDistance ) THEN
            zminm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            zmaxm = MAXVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            IF( zmaxm < zmin - MaxDistance ) GOTO 100 
            IF( zminm > zmax + MaxDistance ) GOTO 100
          END IF
          

          NodesT % x(1) = MAX( xmin, xminm ) 
          NodesT % x(2) = MIN( xmax, xmaxm ) 
          dxcut = ABS( NodesT % x(1)-NodesT % x(2) )

          ! Too small absolute values may result to problems when inverting matrix
          IF( dxcut < 1.0d-12 ) GOTO 100

          ! Too small relative value is irrelevant
          IF( dxcut < 1.0d-8 * dx ) GOTO 100

          sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) THEN
              sgn0 = -1
              AntiPeriodicHits = AntiPeriodicHits + 1
            END IF
          END IF
          
          ElemHits = ElemHits + 1

          IF( SaveElem ) THEN
            FileName = 't'//TRIM(I2S(TimeStep))//'_b'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) NodesM % x(i)
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_e'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,2
              WRITE( 10, * ) NodesT % x(i)
            END DO
            CLOSE( 10 )           
          END IF
                   
          ! Use somewhat higher integration rules than the default
          IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2 ) 
          
          IF(BiOrthogonalBasis) THEN
            MASS  = 0
            CoeffBasis = 0
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                  IP % v(nip),IP % w(nip),detJ,Basis)

              ! Global coordinate of the integration point
              xt = SUM( Basis(1:2) * NodesT % x(1:2) )
            
              ! Integration weight for current integration point
              Wtemp = DetJ * ArcCoeff * IP % s(nip)
            
              ! Integration point at the slave element
              CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

              DO i=1,n
                DO j=1,n
                  MASS(i,j) = MASS(i,j) + wTemp * Basis(i) * Basis(j)
                END DO
                CoeffBasis(i) = CoeffBasis(i) + wTemp * Basis(i)
              END DO
            END DO

            CALL InvertMatrix( MASS, n )

            DO i=1,n
              DO j=1,n
                MASS(i,j) = MASS(i,j) * CoeffBasis(i)
              END DO
            END DO
          END IF


          DO nip=1, IP % n 
            stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                IP % v(nip),IP % w(nip),detJ,Basis)
            
            ! We will actually only use the global coordinates and the integration weight 
            ! from the temporal mesh. 
            
            ! Global coordinate of the integration point
            xt = SUM( Basis(1:2) * NodesT % x(1:2) )
            
            ! Integration weight for current integration point
            ! Use the real arc length so that this projector weights correctly 
            ! in rotational case when used with other projectors.
            Wtemp = ArcCoeff * DetJ * IP % s(nip)
            sumarea = sumarea + Wtemp

            ! Integration point at the slave element
            CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

            ! Integration point at the master element
            CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
            stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
            
            IF(BiOrthogonalBasis) THEN
              CoeffBasis = 0._dp
              DO i=1,n
                DO j=1,n
                  CoeffBasis(i) = CoeffBasis(i) + MASS(i,j) * Basis(j)
                END DO
              END DO
            END IF

            ! Add the entries to the projector
            DO j=1,n 
              jj = Indexes(j)                                    
              nrow = NodePerm(InvPerm1(jj))
              IF( nrow == 0 ) CYCLE
              
              Projector % InvPerm(nrow) = InvPerm1(jj)
              val = Basis(j) * Wtemp
              IF(BiorthogonalBasis) THEN
                val_dual = CoeffBasis(j) * Wtemp
              END IF

              DO i=1,n
                CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                      InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val )

                IF(BiorthogonalBasis ) THEN
                  CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                        InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val_dual )
                END IF
              END DO
              
              DO i=1,nM
                CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                    InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )

                IF(BiorthogonalBasis) THEN
                  IF(DualMaster .OR. DualLCoeff) THEN
                    CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val_dual )
                  ELSE
                    CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )
                  END IF
                END IF
              END DO
            END DO

            ! Add the entries to the dual projector 
            IF( CreateDual ) THEN
              DO j=1,nM 
                jj = IndexesM(j)                                    
                nrow = DualNodePerm(InvPerm2(jj))
                IF( nrow == 0 ) CYCLE
                
                DualProjector % InvPerm(nrow) = InvPerm2(jj)
                val = BasisM(j) * Wtemp

                DO i=1,nM
                  CALL List_AddToMatrixElement(DualProjector % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), sgn0 * NodeCoeff * BasisM(i) * val ) 
                END DO

                DO i=1,n
                  !IF( ABS( val * BasisM(i) ) < 1.0e-10 ) CYCLE
                  CALL List_AddToMatrixElement(DualProjector % ListMatrix, nrow, &
                      InvPerm1(Indexes(i)), -NodeScale * NodeCoeff * Basis(i) * val )                   
                END DO
              END DO
            END IF
          END DO

100       IF( Repeating ) THEN
            IF( NRange2 /= 0 ) THEN
              xminm = xminm + Nrange2 * XRange
              xmaxm = xmaxm + Nrange2 * XRange
              NodesM % x(1:n) = NodesM % x(1:n) + NRange2 * XRange 
              NRange = NRange + NRange2
              NRange2 = 0
              GOTO 200
            END IF
          END IF

        END DO

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_n.dat'
          OPEN( 10,FILE=Filename)
          WRITE( 10, * ) ElemHits 
          CLOSE( 10 )
        END IF
        
        TotHits = TotHits + ElemHits
        TotSumArea = TotSumArea + SumArea
        TotRefArea = TotRefArea + RefArea

        Err = SumArea / RefArea
        IF( Err > MaxErr ) THEN
          MaxErr = Err
          MaxErrInd = Err
        END IF
        IF( Err < MinErr ) THEN
          MinErr = Err
          MinErrInd = ind
        END IF
      END DO

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )
      DEALLOCATE( NodesT % x, NodesT % y, NodesT % z )
      DEALLOCATE( Basis, BasisM )

      CALL Info('LevelProjector','Number of integration pairs: '&
          //TRIM(I2S(TotHits)),Level=10)
      IF( AntiPeriodicHits > 0 ) THEN
        CALL Info('LevelProjector','Number of antiperiodic pairs: '&
          //TRIM(I2S(AntiPeriodicHits)),Level=10)
      END IF

      WRITE( Message,'(A,ES12.5)') 'Total reference length:',TotRefArea / ArcCoeff
      CALL Info('LevelProjector',Message,Level=8) 
      WRITE( Message,'(A,ES12.5)') 'Total integrated length:',TotSumArea / ArcCoeff
      CALL Info('LevelProjector',Message,Level=8)

      Err = TotSumArea / TotRefArea
      WRITE( Message,'(A,ES12.3)') 'Average ratio in length integration:',Err 
      CALL Info('LevelProjector',Message,Level=8)

      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Maximum relative discrepancy in length (element: ',MaxErrInd,'):',MaxErr-1.0_dp 
      CALL Info('LevelProjector',Message,Level=8)
      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Minimum relative discrepancy in length (element: ',MinErrInd,'):',MinErr-1.0_dp 
      CALL Info('LevelProjector',Message,Level=8)


    END SUBROUTINE AddProjectorWeak1D

  END FUNCTION LevelProjector
  !------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> Create a Galerkin projector related to discontinous interface.
!> This uses the information stored when the discontinuous interface 
!> was first coined. This enables simple one-to-one mapping. Integration
!> weight is used for the nodel projector to allow physical jump conditions.
!> For the edge dofs there is no such jumps and hence the projector uses
!> weights of one. 
!---------------------------------------------------------------------------
  FUNCTION WeightedProjectorDiscont(Mesh, bc ) RESULT ( Projector )
    !---------------------------------------------------------------------------
    USE Lists
    USE ListMatrix

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: bc
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: NodePerm(:)
    TYPE(Model_t), POINTER :: Model
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER :: p,q,i,j,it,nn,n,m,t,NoOrigNodes, NoDiscontNodes, indp, indq, &
        e1, e2, e12, i1, i2, j1, j2, sgn, ParentMissing, ParentFound, PosSides, ActSides, &
        InvPermSize, indpoffset
    INTEGER, POINTER :: Rows(:),Cols(:), InvPerm(:)
    REAL(KIND=dp), POINTER :: Values(:), Basis(:), WBasis(:,:), &
                 Wbasis2(:,:),RotWBasis(:,:),dBasisdx(:,:)
    REAL(KIND=dp) :: u,v,w,val,detJ,Scale,x,weight,Coeff
    INTEGER, ALLOCATABLE :: Indexes(:), DiscontIndexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element, Left, Right, OldFace, NewFace, Swap
    LOGICAL :: Stat,DisCont,Found,NodalJump,AxisSym, SetDiag, &
        SetDiagEdges, DoNodes, DoEdges, LocalConstraints, NoHalo
    LOGICAL, ALLOCATABLE :: EdgeDone(:)
    REAL(KIND=dp) :: point(3), uvw(3), DiagEps
    INTEGER, ALLOCATABLE :: EQind(:)
    INTEGER, POINTER :: OldMap(:,:), NewMap(:,:)
    TYPE(ValueList_t), POINTER :: BCParams
    LOGICAL :: CheckHaloNodes
    LOGICAL, POINTER :: HaloNode(:)

    CALL Info('WeightedProjectorDiscont','Creating projector for discontinuous boundary '&
         //TRIM(I2S(bc)),Level=7)

    Projector => NULL()
    IF( .NOT. Mesh % DisContMesh ) THEN
      CALL Warn('WeightedProjectorDiscont','Discontinuous mesh not created?')
      RETURN
    END IF

    Model => CurrentModel

    j = 0
    DO i=1,Model % NumberOfBCs
      IF( ListGetLogical(Model % BCs(i) % Values,'Discontinuous Boundary',Found) ) THEN
        j = j + 1
      END IF
    END DO
    IF( j > 1 ) THEN
      CALL Warn('WeightedProjectorDiscont','One BC (not '&
          //TRIM(I2S(j))//') only for discontinuous boundary!')
    END IF
 
    BCParams => Model % BCs(bc) % Values

    Scale = ListGetCReal( BCParams,'Mortar BC Scaling',Stat )  
    IF(.NOT. Stat) Scale = -1.0_dp

    NodalJump = ListCheckPrefix( BCParams,'Mortar BC Coefficient')
    IF(.NOT. NodalJump ) THEN
      NodalJump = ListCheckPrefix( BCParams,'Mortar BC Resistivity')
    END IF

    ! Take the full weight when creating the constraints since the values will 
    ! not be communicated
    LocalConstraints = ListGetLogical(Model % Solver % Values, &
        'Partition Local Projector',Found)
    IF(.NOT. Found ) LocalConstraints = ListGetLogical(Model % Solver % Values, &
        'Partition Local Constraints',Found)

    ! Don't consider halo when creating discontinuity
    NoHalo = ListGetLogical(Model % Solver % Values, &
        'Projector No Halo',Found)

    ! Don't consider single halo nodes when creating discontinuity
    CheckHaloNodes = ListGetLogical( Model % Solver % Values,&
        'Projector No Halo Nodes',Found ) 
    IF( CheckHaloNodes ) THEN
      CALL MarkHaloNodes( Mesh, HaloNode, CheckHaloNodes )
    END IF


    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',Found ) ) THEN
      DoEdges = .FALSE. 
    ELSE IF( ListGetLogical( BCParams,'Projector Skip Edges',Found ) ) THEN
      DoEdges = .FALSE.
    ELSE
      DoEdges = ( Mesh % NumberOfEdges > 0 )
    END IF
    IF( DoEdges .AND. Mesh % NumberOfEdges == 0 ) THEN
      CALL Warn('WeightedProjectorDiscont','Edge basis requested but mesh has no edges!')
      DoEdges = .FALSE.
    END IF

    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',Found ) ) THEN
      DoNodes = .FALSE. 
    ELSE IF( ListGetLogical( BCParams,'Projector Skip Nodes',Found ) ) THEN
      DoNodes = .FALSE.
    ELSE
      DoNodes = ( Mesh % NumberOfNodes > 0 )
    END IF

    ! Should the projector be diagonal or mass matrix type 
    SetDiag = ListGetLogical( BCParams,'Mortar BC Diag',Found ) 

    IF(.NOT. Found ) SetDiag = ListGetLogical( BCParams, 'Use Biorthogonal Basis', Found)

    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      SetDiag = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( SetDiag ) THEN
        CALL Info('WeightedProjectorDiscont',&
            'Setting > Use Biorthogonal Basis < to True to enable elimination',Level=8)
      END IF
    END IF


    SetDiagEdges = ListGetLogical( BCParams,'Mortar BC Diag Edges',Found )
    IF(.NOT. Found ) SetDiagEdges = SetDiag
    DiagEps = ListGetConstReal( BCParams,'Mortar BC Diag Eps',Found ) 

    ! Integration weights should follow the metrics if we want physical nodal jumps. 
    AxisSym = .FALSE.
    IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) THEN
      IF( NodalJump ) THEN
        AxisSym = .TRUE.
      ELSE IF (ASSOCIATED(CurrentModel % Solver)) THEN
        AxisSym = ListGetLogical(CurrentModel % Solver % Values,'Projector Metrics',Found)
      END IF
      IF( AxisSym ) CALL Info('weightedProjectorDiscont','Projector will be weighted for axi symmetry',Level=7)
    END IF


    n = Mesh % MaxElementDOFs
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    ALLOCATE( Indexes(n), DisContIndexes(n), Basis(n), Wbasis(n,3), &
            Wbasis2(n,3), dBasisdx(n,3), RotWBasis(n,3) )
    Indexes = 0
    Basis = 0.0_dp
    DiscontIndexes = 0

    NodePerm => Mesh % DisContPerm
    NoOrigNodes = SIZE( NodePerm ) 
    NoDiscontNodes = COUNT( NodePerm > 0 ) 

    IF( DoNodes ) THEN
      indpoffset = NoDiscontNodes
    ELSE
      indpoffset = 0
    END IF
    InvPerm => NULL()
    InvPermSize = indpoffset
    
    ! Compute the number of potential edges. This mimics the loop that really creates the projector 
    ! below. 
    IF( DoEdges ) THEN
      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )
      EdgeDone = .FALSE.
      indp = indpoffset

      DO t = 1, Mesh % NumberOfBoundaryElements
        
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )        
        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 
        
        IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
          CYCLE
        END IF

        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF 
        IF( NoHalo .AND. ActSides == 0 ) CYCLE
        
        ! Consistently choose the face with the old edges 
        IF( ALL( Left % NodeIndexes <= NoOrigNodes ) ) THEN
          OldFace => Left
        ELSE IF( ALL( Right % NodeIndexes <= NoOrigNodes ) ) THEN
          OldFace => Right
        ELSE
          CALL Warn('WeightedProjectorDiscont','Neither face is purely old!')
          CYCLE
        END IF

        OldMap => LGetEdgeMap( OldFace % TYPE % ElementCode / 100)

        DO i = 1,OldFace % TYPE % NumberOfEdges          
          e1 = OldFace % EdgeIndexes(i)
          IF( EdgeDone(e1) ) CYCLE

          i1 = OldFace % NodeIndexes( OldMap(i,1) )
          i2 = OldFace % NodeIndexes( OldMap(i,2) )
                    
          ! i1 and i2 were already checked to be "old" nodes
          IF( NodePerm(i1) == 0 ) CYCLE
          IF( NodePerm(i2) == 0 ) CYCLE

          indp = indp + 1
          EdgeDone(e1) = .TRUE.
        END DO
      END DO
      InvPermSize = indp
      CALL Info('WeightedProjectorDiscont',&
          'Size of InvPerm estimated to be: '//TRIM(I2S(InvPermSize)),Level=8)
    END IF

    ! Ok, nothing to do just go end tidy things up
    IF( InvPermSize == 0 ) GOTO 100

    ! Create a list matrix that allows for unspecified entries in the matrix 
    ! structure to be introduced.
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN
    Projector % ProjectorBC = bc
    
    ! Create the inverse permutation needed when the projector matrix is added to the global 
    ! matrix. 
    ALLOCATE( Projector % InvPerm( InvPermSize ) )
    InvPerm => Projector % InvPerm
    InvPerm = 0

    
    ! Projector for the nodal dofs. 
    !------------------------------------------------------------------------
    IF( DoNodes ) THEN

      ParentMissing = 0
      ParentFound = 0
      DO t = 1, Mesh % NumberOfBoundaryElements

        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        n = Element % TYPE % NumberOfNodes        
        Indexes(1:n) = Element % NodeIndexes(1:n)

        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 

        ! Here we really need both sides to be able to continue!
        !IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
        !  ParentMissing = ParentMissing + 1
        !  CYCLE
        !END IF

        PosSides = 0
        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          PosSides = PosSides + 1
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          PosSides = PosSides + 1
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF
        IF( NoHalo .AND. ActSides == 0 ) CYCLE        

        IF( LocalConstraints ) THEN
          Coeff = 1.0_dp
        ELSE
          Coeff = 1.0_dp * ActSides / PosSides 
        END IF
        IF( ABS( Coeff ) < TINY( 1.0_dp ) ) CYCLE

        ParentFound = ParentFound + 1

        ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))

        IF( ALL( NodePerm(Indexes(1:n)) == 0 ) ) CYCLE
        
        IF( CheckHaloNodes ) THEN
          IF( ALL( HaloNode(Indexes(1:n)) ) ) CYCLE
        END IF

        ! Get the indexes on the other side of the discontinuous boundary
        DO i=1,n
          j = NodePerm( Indexes(i) ) 
          IF( j == 0 ) THEN
            DiscontIndexes(i) = Indexes(i)
          ELSE
            DiscontIndexes(i) = j + NoOrigNodes
          END IF
        END DO

        IntegStuff = GaussPoints( Element )
        DO j=1,IntegStuff % n
          u = IntegStuff % u(j)
          v = IntegStuff % v(j)
          w = IntegStuff % w(j)

          Stat = ElementInfo(Element, ElementNodes, u, v, w, detJ, Basis)

          weight = Coeff * detJ * IntegStuff % s(j)
          IF( AxisSym ) THEN
            x = SUM( Basis(1:n) * ElementNodes % x(1:n) )
            weight = weight * x
          END IF

          DO p=1,n             
            indp = NodePerm( Indexes(p) )
            IF( indp == 0 ) CYCLE
            IF( CheckHaloNodes ) THEN
              IF( HaloNode( Indexes(p) ) ) CYCLE
            END IF

            val = weight * Basis(p)

            ! Only set for the nodes are are really used
            InvPerm(indp) = Indexes(p)

            IF( SetDiag ) THEN
              CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                  Indexes(p), val ) 

              CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                  DiscontIndexes(p), Scale * val )             
            ELSE
              DO q=1,n

                indq = NodePerm(Indexes(q))
                IF( indq == 0 ) CYCLE

                IF( CheckHaloNodes ) THEN
                  IF( HaloNode( Indexes(p) ) ) CYCLE
                END IF
                
                CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                    Indexes(q), Basis(q) * val ) 
                CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                    DiscontIndexes(q), Scale * Basis(q) * val ) 
              END DO
            END IF
          END DO
        END DO
      END DO
      IF( ParentMissing > 0 ) THEN
        CALL Warn('WeightedProjectorDiscont','Number of half-sided discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentMissing)) )
        CALL Warn('WeightedProjectorDiscont','Number of proper discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentFound)) )
      END IF
      CALL Info('WeightedProjectorDiscont','Created projector for '&
          //TRIM(I2S(NoDiscontNodes))//' discontinuous nodes',Level=10)
    END IF


    ! Create the projector also for edge dofs if they exist and are
    ! requested. 
    !----------------------------------------------------------------
    IF( DoEdges ) THEN
      ParentMissing = 0
      ParentFound = 0
      n = Mesh % NumberOfNodes

      val = 1.0_dp
      Scale = 1.0_dp

      indp = indpoffset
      ALLOCATE( Eqind(Mesh % NumberOfEdges) ); EQind = 0

      DO t = 1, Mesh % NumberOfBoundaryElements
        
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        
        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 
        
        ! Here we really need both sides to be able to continue!
        IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
          ParentMissing = ParentMissing + 1
          CYCLE
        END IF

        PosSides = 0
        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          PosSides = PosSides + 1
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          PosSides = PosSides + 1
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF

        IF( NoHalo .AND. ActSides == 0 ) CYCLE

        IF( LocalConstraints ) THEN
          Coeff = 1.0_dp
        ELSE          
          Coeff = (1.0_dp * ActSides) / (1.0_dp * PosSides)
        END IF

        ! Consistently choose the face with the old edges
        IF( ALL( Left % NodeIndexes <= NoOrigNodes ) ) THEN
        ELSE IF( ALL( Right % NodeIndexes <= NoOrigNodes ) ) THEN
          swap  => Left
          Left  => Right
          Right => swap
        ELSE
          ! We already complained once
          CYCLE
        END IF

        OldFace => Find_Face( Left, Element )
        nn = SIZE(Element % NodeIndexes)
        Indexes(1:nn) = Element % NodeIndexes
        Element % NodeIndexes = NodePerm(Indexes(1:nn)) + NoOrigNodes
        NewFace => Find_Face( Right, Element )
        Element % NodeIndexes = Indexes(1:nn)
 
        ParentFound = ParentFound + 1

        OldMap => LGetEdgeMap( OldFace % TYPE % ElementCode / 100 )
        NewMap => LGetEdgeMap( NewFace % TYPE % ElementCode / 100 )

        IntegStuff = GaussPoints( oldface )
        DO it = 1,IntegStuff % n
          u = integstuff % u(it)
          v = integstuff % v(it)
          w = integstuff % w(it)

          nn = OldFace % TYPE % NumberOfNodes
          ElementNodes % x(1:nn) = Mesh % Nodes % x(oldface % NodeIndexes(1:nn))
          ElementNodes % y(1:nn) = Mesh % Nodes % y(oldface % NodeIndexes(1:nn))
          ElementNodes % z(1:nn) = Mesh % Nodes % z(oldface % NodeIndexes(1:nn))

          Stat = ElementInfo( OldFace, ElementNodes,u,v,w, DetJ, Basis,dBasisdx )
          CALL GetEdgeBasis( OldFace, Wbasis, RotWbasis, Basis, dBasisdx )

          Point(1) = SUM(Basis(1:nn) * ElementNodes % x(1:nn))
          Point(2) = SUM(Basis(1:nn) * ElementNodes % y(1:nn))
          Point(3) = SUM(Basis(1:nn) * ElementNodes % z(1:nn))

          nn = NewFace % TYPE % NumberOfNodes
          ElementNodes % x(1:nn) = Mesh % Nodes % x(newface % NodeIndexes(1:nn))
          ElementNodes % y(1:nn) = Mesh % Nodes % y(newface % NodeIndexes(1:nn))
          ElementNodes % z(1:nn) = Mesh % Nodes % z(newface % NodeIndexes(1:nn))

          Found = PointInElement( NewFace, ElementNodes, Point, uvw )
          u = uvw(1); v=uvw(2); w=uvw(3)
          Stat = ElementInfo(NewFace, ElementNodes,u,v,w, detj, Basis,dbasisdx )
          CALL GetEdgeBasis( NewFace, Wbasis2, RotwBasis, Basis, dBasisdx )

          Weight = detJ * IntegStuff % s(it) * Coeff
        
          ! Go through combinations of edges and find the edges for which the 
          ! indexes are the same. 
          DO i = 1,OldFace % TYPE % NumberOfEdges
            e1 = OldFace % EdgeIndexes(i)

            IF ( EQind(e1) == 0 ) THEN
              indp = indp + 1
              EQind(e1) = indp
              InvPerm(indp) = n + e1
            END IF

            IF( SetDiagEdges ) THEN
              i1 = OldFace % NodeIndexes( OldMap(i,1) )
              i1 = NoOrigNodes + NodePerm(i1)
              i2 = OldFace % NodeIndexes( OldMap(i,2) )
              i2 = NoOrigNodes + NodePerm(i2)

              DO j = 1,NewFace % TYPE % NumberOfEdges
                j1 = NewFace % NodeIndexes( NewMap(j,1) )
                j2 = NewFace % NodeIndexes( NewMap(j,2) )
                IF (i1==j1 .AND. i2==j2 .OR. i1==j2 .AND. i2==j1 ) EXIT
              END DO
              val = Weight * SUM(WBasis(i,:) * Wbasis(i,:))
              IF ( ABS(Val)>= 10*AEPS ) &
                  CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e1, Val )
              
              e2  = NewFace % EdgeIndexes(j)
              val = Weight * SUM(WBasis(i,:) * Wbasis2(j,:))
              IF ( ABS(val) >= 10*AEPS ) &
                  CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e2, -Val )              
            ELSE
              DO j = 1,NewFace % TYPE % NumberOfEdges
                e2  = NewFace % EdgeIndexes(j)
                e12 = OldFace % EdgeIndexes(j)
                
                val = Weight * SUM(WBasis(i,:) * Wbasis(j,:))
                IF ( ABS(Val)>= 10*AEPS ) &
                    CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e12, Val )
                
                val = Weight * SUM(WBasis(i,:) * Wbasis2(j,:))
                IF ( ABS(val) >= 10*AEPS ) &
                    CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e2, -Val )
              END DO
            END IF

          END DO
        END DO
      END DO

      DEALLOCATE( EdgeDone )
      IF( .NOT. DoNodes .AND. ParentMissing > 0 ) THEN
        CALL Warn('WeightedProjectorDiscont','Number of half-sided discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentMissing)) )
        CALL Warn('WeightedProjectorDiscont','Number of proper discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentFound)) )
      END IF
      CALL Info('WeightedProjectorDiscont','Created projector for '&
          //TRIM(I2S(indp-NoDiscontNodes))//' discontinuous edges',Level=10)
    END IF

    ! Convert from list matrix to CRS matrix format
    CALL List_ToCRSMatrix(Projector)

    IF( Projector % NumberOfRows > 0) THEN
      CALL CRS_SortMatrix(Projector,.TRUE.)
      CALL Info('WeightedProjectorDiscont','Number of entries in projector matrix: '//&
          TRIM(I2S(SIZE(Projector % Cols)) ), Level=9)
    ELSE
      CALL FreeMatrix(Projector); Projector=>NULL()
    END IF

100 DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
    DEALLOCATE( Indexes, DisContIndexes, Basis, dBasisdx, WBasis, WBasis2, RotWBasis )
    IF( CheckHaloNodes ) DEALLOCATE( HaloNode )

           
  END FUNCTION WeightedProjectorDiscont
  !------------------------------------------------------------------------------
 

  !---------------------------------------------------------------------------
  ! Simply fitting of cylinder into a point cloud. This is done in two phases.
  ! 1) The axis of the cylinder is found by minimizing the \sum((n_i*t)^2)
  !    for each component of of t where n_i:s are the surface normals. 
  !    This is fully generic and assumes no positions. 
  ! 2) The radius and center point of the cylinder are found by fitting a circle
  !    in the chosen plane to three representative points. Currently the fitting
  !    can only be done in x-y plane. 
  !---------------------------------------------------------------------------
  SUBROUTINE CylinderFit(PMesh, PParams) 
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Valuelist_t), POINTER :: PParams

    INTEGER :: i,j,k,n,t,AxisI,iter
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: NiNj(3,3),A(3,3),F(3),M11,M12,M13,M14
    REAL(KIND=dp) :: d1,d2,MinDist,MaxDist,Dist,X0,Y0,Rad
    REAL(KIND=dp) :: Normal(3), AxisNormal(3), Tangent1(3), Tangent2(3), Coord(3), &
        CircleCoord(3,3)
    INTEGER :: CircleInd(3) 

    CALL Info('CylinderFit','Trying to fit a cylinder to the surface patch',Level=10)

    NiNj = 0.0_dp

    n = PMesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    ! If the initial mesh is in 2D there is really no need to figure out the 
    ! direction of the rotational axis. It can only be aligned with the z-axis. 
    IF( CurrentModel % Mesh % MeshDim == 2 ) THEN
      AxisNormal = 0.0_dp
      AxisNormal(3) = 1.0_dp
      GOTO 100 
    END IF


    ! Compute the inner product of <N*N> for the elements
    DO t=1, PMesh % NumberOfBulkElements
      Element => PMesh % Elements(t)
      
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      Nodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
      Nodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
      Nodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
      
      Normal = NormalVector( Element, Nodes, Check = .FALSE. ) 

      DO i=1,3
        DO j=1,3
          NiNj(i,j) = NiNj(i,j) + Normal(i) * Normal(j)
        END DO
      END DO      
    END DO

    ! Normalize by the number of boundary elements
    NiNj = NiNj / PMesh % NumberOfBulkElements

    ! The potential direction for the cylinder axis is the direction with 
    ! least hits for the normal.
    AxisI = 1 
    DO i=2,3
      IF( NiNj(i,i) < NiNj(AxisI,AxisI) ) AxisI = i 
    END DO

    CALL Info('CylinderFit','Axis coordinate set to be: '//TRIM(I2S(AxisI)))

    ! Keep the dominating direction fixed and iteratively solve the two other directions
    AxisNormal = 0.0_dp
    AxisNormal(AxisI) = 1.0_dp

    ! Basically we could solve from equation Ax=0 the tangent but only up to a constant.
    ! Thus we enforce the axis direction to one by manipulation the matrix equation 
    ! thereby can get a unique solution. 
    A = NiNj
    A(AxisI,1:3) = 0.0_dp
    A(AxisI,AxisI) = 1.0_dp
    CALL InvertMatrix( A, 3 )
    AxisNormal = A(1:3,AxisI)

    ! Normalize the axis normal length to one    
    AxisNormal = AxisNormal / SQRT( SUM( AxisNormal ** 2 ) )
    IF( 1.0_dp - ABS( AxisNormal(3) ) > 1.0e-5 ) THEN
      CALL Warn('CylinderFit','The cylinder axis is not aligned with z-axis!')
    END IF

100 CALL TangentDirections( AxisNormal,Tangent1,Tangent2 )

    IF(.FALSE.) THEN
      PRINT *,'Axis Normal:',AxisNormal
      PRINT *,'Axis Tangent 1:',Tangent1
      PRINT *,'Axis Tangent 2:',Tangent2
    END IF

    ! Finding three points with maximum distance in the tangent directions

    ! First, find the single extremum point in the first tangent direction
    ! Save the local coordinates in the N-T system of the cylinder
    MinDist = HUGE(MinDist) 
    DO i=1, PMesh % NumberOfNodes
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)

      d1 = SUM( Tangent1 * Coord )
      IF( d1 < MinDist ) THEN
        MinDist = d1
        CircleInd(1) = i
      END IF
    END DO

    i = CircleInd(1)
    Coord(1) = PMesh % Nodes % x(i)
    Coord(2) = PMesh % Nodes % y(i)
    Coord(3) = PMesh % Nodes % z(i)
      
    CircleCoord(1,1) = SUM( Tangent1 * Coord ) 
    CircleCoord(1,2) = SUM( Tangent2 * Coord ) 
    CircleCoord(1,3) = SUM( AxisNormal * Coord )
   

    !PRINT *,'MinDist1:',MinDist,CircleInd(1),CircleCoord(1,:)

    ! Find two more points such that their minimum distance to the previous point(s)
    ! is maximized. This takes some time but the further the nodes are apart the more 
    ! accurate it will be to fit the circle to the points. Also if there is just 
    ! a symmetric section of the cylinder it is important to find the points rigorously.
    DO j=2,3
      ! The maximum minimum distance of any node from the previously defined nodes
      MaxDist = 0.0_dp
      DO i=1, PMesh % NumberOfNodes
        Coord(1) = PMesh % Nodes % x(i)
        Coord(2) = PMesh % Nodes % y(i)
        Coord(3) = PMesh % Nodes % z(i)
        
        ! Minimum distance from the previously defined nodes
        MinDist = HUGE(MinDist)
        DO k=1,j-1
          d1 = SUM( Tangent1 * Coord )
          d2 = SUM( Tangent2 * Coord )
          Dist = ( d1 - CircleCoord(k,1) )**2 + ( d2 - CircleCoord(k,2) )**2
          MinDist = MIN( Dist, MinDist )
        END DO
        
        ! If the minimum distance is greater than in any other node, choose this
        IF( MaxDist < MinDist ) THEN
          MaxDist = MinDist 
          CircleInd(j) = i
        END IF
      END DO

      ! Ok, we have found the point now set the circle coordinates 
      i = CircleInd(j)
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(j,1) = SUM( Tangent1 * Coord ) 
      CircleCoord(j,2) = SUM( Tangent2 * Coord ) 
      CircleCoord(j,3) = SUM( AxisNormal * Coord )
    END DO
      

    !PRINT *,'Circle Indexes:',CircleInd

    ! Given three nodes it is possible to analytically compute the center point and
    ! radius of the cylinder from a 4x4 determinant equation. The matrices values
    ! m1i are the determinants of the comatrices. 

    A(1:3,1) = CircleCoord(1:3,1)  ! x
    A(1:3,2) = CircleCoord(1:3,2)  ! y
    A(1:3,3) = 1.0_dp
    m11 = Det3x3( a )

    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,2)  ! y
    A(1:3,3) = 1.0_dp
    m12 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,1)  ! x
    A(1:3,3) = 1.0_dp
    m13 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2 ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,1)  ! x
    A(1:3,3) = CircleCoord(1:3,2)  ! y
    m14 = Det3x3( a )

    !PRINT *,'determinants:',m11,m12,m13,m14

    IF( ABS( m11 ) < EPSILON( m11 ) ) THEN
      CALL Fatal('CylinderFit','Points cannot be an a circle')
    END IF

    X0 =  0.5 * m12 / m11 
    Y0 = -0.5 * m13 / m11
    rad = SQRT( x0**2 + y0**2 + m14/m11 )

    Coord = x0 * Tangent1 + y0 * Tangent2

    !PRINT *,'Center point in cartesian coordinates:',Coord
    
    CALL ListAddConstReal( PParams,'Rotational Projector Center X',Coord(1))
    CALL ListAddConstReal( PParams,'Rotational Projector Center Y',Coord(2))
    CALL ListAddConstReal( PParams,'Rotational Projector Center Z',Coord(3))

    CALL ListAddConstReal( PParams,'Rotational Projector Normal X',AxisNormal(1))
    CALL ListAddConstReal( PParams,'Rotational Projector Normal Y',AxisNormal(2))
    CALL ListAddConstReal( PParams,'Rotational Projector Normal Z',AxisNormal(3))

    
  CONTAINS
    
    ! Compute the value of 3x3 determinant
    !-------------------------------------------
    FUNCTION Det3x3( A ) RESULT ( val ) 
      
      REAL(KIND=dp) :: A(:,:)
      REAL(KIND=dp) :: val

      val = A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) &
          - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) ) &
          + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) 

    END FUNCTION Det3x3

  END SUBROUTINE CylinderFit



  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming rotating boundaries make 
  !> a coordinate transformation to (phi,z) level where the interpolation
  !> accuracy is not limited by the curvilinear coordinates. Also ensure
  !> that the master nodes manipulated so they for sure hit the target nodes.
  !---------------------------------------------------------------------------
  SUBROUTINE RotationalInterfaceMeshes(BMesh1, BMesh2, BParams, Cylindrical, &
      Radius, FullCircle )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    REAL(KIND=dp) :: Radius
    LOGICAL :: FullCircle, Cylindrical
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3),&
        x1r_min(3),x1r_max(3),x2r_min(3),x2r_max(3)
    REAL(KIND=dp) :: x(3), xcyl(3),rad2deg,F1min,F1max,F2min,F2max,dFii1,dFii2,eps_rad,&
        err1,err2,dF,Fii,Fii0,Nsymmetry,fmin,fmax,DegOffset,rad,alpha,x0(3),xtmp(3),&
        Normal(3), Tangent1(3), Tangent2(3) 
    REAL(KIND=dp), POINTER :: TmpCoord(:)
    REAL(KIND=dp),ALLOCATABLE :: Angles(:)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,ind,Nmax,Nmin,Nfii,Nnodes,MaxElemNodes,NElems
    LOGICAL :: Found, Hit0, Hit90, Hit180, Hit270, SetDegOffset
    LOGICAL :: GotNormal, GotCenter, MoveAngle

    ! We choose degrees as they are more intuitive
    rad2deg = 180.0_dp / PI
    MaxElemNodes = BMesh2 % MaxElementNodes 
    ALLOCATE( Angles(MaxElemNodes) )
    
    Nnodes = BMesh2 % NumberOfNodes
    NElems = BMesh2 % NumberOfBulkElements
    FullCircle = .FALSE.

    ! Cylindrical projector is fitted always and rotational only when requested.
    IF( ListGetLogical( BParams,'Rotational Projector Center Fit',Found ) .OR. &
       Cylindrical ) THEN
      IF( .NOT. ListCheckPresent( BParams,'Rotational Projector Center X') ) THEN
        CALL CylinderFit( BMesh1, BParams ) 
      END IF
    END IF
    
    x0(1) = ListGetCReal( BParams,'Rotational Projector Center X',GotCenter ) 
    x0(2) = ListGetCReal( BParams,'Rotational Projector Center Y',Found ) 
    GotCenter = GotCenter .OR. Found
    x0(3) = ListGetCReal( BParams,'Rotational Projector Center Z',Found ) 
    GotCenter = GotCenter .OR. Found

    Normal(1) = ListGetCReal( BParams,'Rotational Projector Normal X',GotNormal ) 
    Normal(2) = ListGetCReal( BParams,'Rotational Projector Normal Y',Found ) 
    GotNormal = GotNormal .OR. Found
    Normal(3) = ListGetCReal( BParams,'Rotational Projector Normal Z',Found ) 
    GotNormal = GotNormal .OR. Found

    IF( GotNormal ) THEN
      CALL TangentDirections( Normal,Tangent1,Tangent2 )
    END IF

    ! Go trough master (k=1) and target mesh (k=2)
    !--------------------------------------------
    DO k=1,2
     
      ! Potentially the projector may be set to rotate by just adding an offset 
      ! to the angle. This may depende on time etc. 
      IF( k == 1 ) THEN
        DegOffset = ListGetCReal(BParams,'Rotational Projector Angle Offset',SetDegOffset ) 
      ELSE
        SetDegOffset = .FALSE.
      END IF

      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      ! Check the initial bounding boxes
      !---------------------------------------------------------------------------
      x2_min(1) = MINVAL( PMesh % Nodes % x )
      x2_min(2) = MINVAL( PMesh % Nodes % y )
      x2_min(3) = MINVAL( PMesh % Nodes % z )
      
      x2_max(1) = MAXVAL( PMesh % Nodes % x )
      x2_max(2) = MAXVAL( PMesh % Nodes % y )
      x2_max(3) = MAXVAL( PMesh % Nodes % z )
      
      IF( k == 1 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Initial extrema for this boundary (x,y,z)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Initial extrema for target boundary (x,y,z)',Level=8)
      END IF
      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2_min(i),x2_max(i)
        CALL Info('RotationalInterfaceMeshes',Message,Level=8)    
      END DO

      ! Memorize the bounding box of the master mesh
      !--------------------------------------------------------------------------
      IF( k == 1 ) THEN
        x1_min = x2_min
        x1_max = x2_max
      END IF

      ! Do the actual coordinate transformation
      !---------------------------------------------------------------------------
      n = PMesh % NumberOfNodes
      DO i=1,n
        x(1) = PMesh % Nodes % x(i)
        x(2) = PMesh % Nodes % y(i)
        x(3) = PMesh % Nodes % z(i)

        ! Subtract the center of axis
        IF( GotCenter ) THEN
          x = x - x0
        END IF

        IF( GotNormal ) THEN
          xtmp = x
          x(1) = SUM( Tangent1 * xtmp ) 
          x(2) = SUM( Tangent2 * xtmp ) 
          x(3) = SUM( Normal * xtmp ) 
        END IF


        ! Set the angle to be the first coordinate as it may sometimes be the 
        ! only nonzero coordinate. Z-coordinate is always unchanged. 
        !------------------------------------------------------------------------
        alpha = rad2deg * ATAN2( x(2), x(1)  ) 
        rad = SQRT( x(1)**2 + x(2)**2)

        ! Set the offset and revert then the angle to range [-180,180] 
        IF( SetDegOffset ) THEN
          alpha = MODULO( alpha + DegOffset, 360.0 )            
          IF( alpha > 180.0 ) alpha = alpha - 360.0
        END IF

        PMesh % Nodes % x(i) = alpha
        PMesh % Nodes % y(i) = x(3)
        PMesh % Nodes % z(i) = rad      
      END DO
      

      ! For cylindrical projector follow exactly the same logic for slave and master
      !------------------------------------------------------------------------------
      IF( Cylindrical .AND. k == 2 ) THEN
        IF( MoveAngle ) THEN
          CALL Info('RotationalInterfaceMeshes','Moving the 2nd mesh discontinuity to same angle',Level=6)
          DO j=1,PMesh % NumberOfNodes
            IF( PMesh % Nodes % x(j) < Fii0 ) PMesh % Nodes % x(j) = &
                PMesh % Nodes % x(j) + 360.0_dp
          END DO
        END IF
      ELSE
        ! Let's see if we have a full angle to operate or not.
        ! If not, then make the interval continuous. 
        ! Here we check only four critical angles: (0,90,180,270) degs.
        Hit0 = .FALSE.; Hit90 = .FALSE.; Hit180 = .FALSE.; Hit270 = .FALSE.
        MoveAngle = .FALSE.; Fii = 0.0_dp; Fii0 = 0.0_dp
        
        DO i=1, PMesh % NumberOfBulkElements
          Element => PMesh % Elements(i)
          n = Element % TYPE % NumberOfNodes        
          NodeIndexes => Element % NodeIndexes
          Angles(1:n) = PMesh % Nodes % x(NodeIndexes)
          
          fmin = MINVAL( Angles(1:n) ) 
          fmax = MAXVAL( Angles(1:n) )
          
          IF( fmax - fmin > 180.0_dp ) THEN
            Hit180 = .TRUE.
          ELSE
            IF( fmax >= 0.0 .AND. fmin <= 0.0 ) Hit0 = .TRUE.
            IF( fmax >= 90.0 .AND. fmin <= 90.0 ) Hit90 = .TRUE.
            IF( fmax >= -90.0 .AND. fmin <= -90.0 ) Hit270 = .TRUE.
          END IF
        END DO
        FullCircle = Hit0 .AND. Hit90 .AND. Hit180 .AND. Hit270
        
        ! Eliminate the problematic discontinuity in case we have no full circle
        ! The discontinuity will be moved to some of angles (-90,0,90).
        IF( FullCircle ) THEN
          CALL Info('RotationalInterfaceMeshes','Cylindrical interface seems to be a full circle',&
              Level=6)
        ELSE IF( Hit180 ) THEN
          MoveAngle = .TRUE.
          IF( .NOT. Hit0 ) THEN
            Fii = 0.0_dp
          ELSE IF( .NOT. Hit270 ) THEN
            Fii = -90.0
          ELSE IF( .NOT. Hit90 ) THEN
            Fii = 90.0
          END IF

          DO j=1,PMesh % NumberOfNodes
            IF( PMesh % Nodes % x(j) < Fii ) PMesh % Nodes % x(j) = &
                PMesh % Nodes % x(j) + 360.0_dp
          END DO
          WRITE( Message,'(A,F8.3)') 'Moving discontinuity of angle to: ',Fii
          Fii0 = Fii
          CALL Info('RotationalInterfaceMesh',Message,Level=6)
        END IF
      END IF


      ! Check the transformed bounding boxes
      !---------------------------------------------------------------------------
      x2r_min(1) = MINVAL( PMesh % Nodes % x )
      x2r_min(2) = MINVAL( PMesh % Nodes % y )
      x2r_min(3) = MINVAL( PMesh % Nodes % z )
      
      x2r_max(1) = MAXVAL( PMesh % Nodes % x )
      x2r_max(2) = MAXVAL( PMesh % Nodes % y )
      x2r_max(3) = MAXVAL( PMesh % Nodes % z )
      
      IF( k == 1 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Transformed extrema for this boundary (phi,z,r)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Transformed extrema for target boundary (phi,z,r)',Level=8)
      END IF
      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2r_min(i),x2r_max(i)
        CALL Info('RotationalInterfaceMeshes',Message,Level=8)    
      END DO

      IF( x2r_min(3) < EPSILON( Radius ) ) THEN
        CALL Fatal('RotationalInterfaceMeshes','Radius cannot be almost zero!')
      END IF

      ! Memorize the bounding box for the 1st mesh
      IF( k == 1 ) THEN
        x1r_min = x2r_min
        x1r_max = x2r_max
      END IF
    END DO

    eps_rad = 1.0d-3 

    ! Choose radius to be max radius of this boundary
    Radius = x1r_max(3) 
    
    err1 = ( x1r_max(3) - x1r_min(3) ) / Radius
    err2 = ( x2r_max(3) - x2r_min(3) ) / Radius

    WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant radius:',err1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant radius:',err2
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    IF( err1 > eps_rad .OR. err2 > eps_rad ) THEN
      CALL Warn('RotationalInterfaceMeshes','Discrepancy of radius is rather large!')
    END IF

    ! Ok, so we have concluded that the interface has constant radius
    ! therefore the constant radius may be removed from the mesh description.
    ! Or perhaps we don't remove to allow more intelligent projector building 
    ! for contact mechanics. 
    !---------------------------------------------------------------------------
    !Bmesh1 % Nodes % z = 0.0_dp
    !BMesh2 % Nodes % z = 0.0_dp

    ! Check whether the z-coordinate is constant or not.
    ! Constant z-coordinate implies 1D system, otherwise 2D system.
    !---------------------------------------------------------------------------
    err1 = ( x1r_max(2) - x1r_min(2) ) / Radius
    err2 = ( x2r_max(2) - x2r_min(2) ) / Radius
    
    IF( err1 < eps_rad .AND. err2 < eps_rad ) THEN
      CALL Info('RotationalInterfaceMeshes','The effective interface meshes are 1D')
      Bmesh1 % Nodes % y = 0.0_dp
      Bmesh2 % Nodes % y = 0.0_dp
    ELSE
      CALL Info('RotationalInterfaceMeshes','The effective interface meshes are 2D')
    END IF

    ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
    Bmesh1 % MeshDim = 2
    Bmesh2 % MeshDim = 2      

    ! Cylindrical interface does not have symmetry as does the rotational!
    IF( Cylindrical .OR. FullCircle ) RETURN

    ! If were are studying a symmetric segment then anylyze further the angle 
    !-------------------------------------------------------------------------
    dFii1 = x1r_max(1)-x1r_min(1)
    dFii2 = x2r_max(1)-x2r_min(1)

    WRITE(Message,'(A,ES12.3)') 'This boundary dfii:  ',dFii1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Target boundary dfii:  ',dFii2
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    err1 = 2 * ABS( dFii1 - dFii2 ) / ( dFii1 + dFii2 )
    WRITE(Message,'(A,ES12.3)') 'Discrepancy in dfii:',err1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)        

    i = ListGetInteger(BParams,'Rotational Projector Periods',Found ) 
    IF( .NOT. Found ) THEN
      Nsymmetry = 360.0_dp / dFii2 
      WRITE(Message,'(A,ES12.3)') 'Suggested sections in target:',Nsymmetry
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
      IF( ABS( Nsymmetry - NINT( Nsymmetry ) ) > 0.01 ) THEN          
        IF( dFii1 < dFii2 ) THEN
          CALL Info('RotationalInterfaceMeshes','You might try to switch master and target!',Level=3)
        END IF
        CALL Fatal('RotationalInterfaceMeshes','Check your settings, this cannot be periodic!')
      END IF
      CALL ListAddInteger(BParams,'Rotational Projector Periods', NINT( Nsymmetry ) ) 
    ELSE
      WRITE(Message,'(A,I0)') 'Using enforced number of periods: ',i
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
      Nsymmetry = 360.0_dp / dFii2 
      WRITE(Message,'(A,ES12.3)') 'Suggested number of periods:',Nsymmetry
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
    END IF

  END SUBROUTINE RotationalInterfaceMeshes
!------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming radial boundaries make 
  !> a coordinate transformation to (r,z) level.
  !> This is always a symmetry condition and can not be a contact condition.
  !---------------------------------------------------------------------------
  SUBROUTINE RadialInterfaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3), x(3), r, phi, z, &
        err1, err2, phierr, eps_rad, rad, rad2deg
    INTEGER :: i,j,k

    ! We choose degrees as they are more intuitive
    rad2deg = 180.0_dp / PI
    
    ! Go trough master (k=1) and target mesh (k=2)
    !--------------------------------------------
    DO k=1,2
     
      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      x2_min = HUGE( x2_min )
      x2_max = -HUGE( x2_max )

      ! Loop over all nodes
      !----------------------------------------------------------------------------
      DO i=1,PMesh % NumberOfNodes
        x(1) = PMesh % Nodes % x(i)
        x(2) = PMesh % Nodes % y(i)
        x(3) = PMesh % Nodes % z(i)
        
        ! Do the actual coordinate transformation
        !---------------------------------------------------------------------------
        r = SQRT( x(1)**2 + x(2)**2 )
        phi = rad2deg * ATAN2( x(2), x(1)  )
        z = x(3)

        PMesh % Nodes % x(i) = r
        PMesh % Nodes % y(i) = z
        PMesh % Nodes % z(i) = 0.0_dp

        ! This is just to check a posteriori that the ranges are ok
        x2_min(1) = MIN(r,x2_min(1))
        IF( r > EPSILON( r ) ) THEN
          x2_min(2) = MIN(phi,x2_min(2))
        END IF
        x2_min(3) = MIN(z,x2_min(3))

        x2_max(1) = MAX(r,x2_max(1))
        IF( r > EPSILON(r) ) THEN
          x2_max(2) = MAX(phi,x2_max(2))
        END IF
        x2_max(3) = MAX(z,x2_max(3))
      END DO

      ! Memorize the bounding box of the master mesh
      !--------------------------------------------------------------------------
      IF( k == 1 ) THEN
        x1_min = x2_min
        x1_max = x2_max
      END IF

      IF( k == 1 ) THEN
        CALL Info('RadialInterfaceMeshes',&
            'Transformed extrema for this boundary (phi,r,z)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RadialInterfaceMeshes',&
            'Transformed extrema for target boundary (phi,r,z)',Level=8)
      END IF

      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2_min(i),x2_max(i)
        CALL Info('RadialInterfaceMeshes',Message,Level=8)    
      END DO

      phierr = x2_max(2) - x2_min(2)  
      WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant angle (degs):',phierr
      CALL Info('RadialInterfaceMeshes',Message,Level=8)    
    END DO

    ! Error in radius
    ! Choose radius to be max radius of either boundary
    rad = MAX( x1_max(1), x2_max(1) )    
    err1 = ABS( x1_max(1) - x2_max(1) ) / rad
    err2 = ABS( x1_min(1) - x2_min(1) ) / rad

    WRITE(Message,'(A,ES12.3)') 'Discrepancy in maximum radius:',err1
    CALL Info('RadialInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Discrepancy in minimum radius:',err2
    CALL Info('RadialInterfaceMeshes',Message,Level=8)    

    eps_rad = 1.0e-3
    IF( err1 > eps_rad .OR. err2 > eps_rad ) THEN
      CALL Warn('RadialInterfaceMeshes','Discrepancy of radius may be too large!')
    END IF

    ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
    Bmesh1 % MeshDim = 2
    Bmesh2 % MeshDim = 2      

  END SUBROUTINE RadialInterfaceMeshes
!------------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !> Given two interface meshes flatten them to (x,y) plane.
  !---------------------------------------------------------------------------
  SUBROUTINE FlatInterfaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Bmesh
    INTEGER :: FlatDim, MeshDim, MinDiffI, i, j
    REAL(KIND=dp), POINTER :: Coord(:)
    REAL(KIND=dp) :: Diff, MaxDiff, MinDiff, RelDiff, RelDiff1
    LOGICAL :: Found, ReduceDim

    CALL Info('FlatInterfaceMeshes','Flattening interface meshes to 2D',Level=8)    
    
    MeshDim = CurrentModel % DIMENSION
    FlatDim = ListGetInteger( BParams,'Flat Projector Coordinate',Found,minv=1,maxv=3) 
    ReduceDim = ListGetLogical( BParams,'Flat Projector Reduce Dimension',Found )

    IF(.NOT. Found ) THEN
      DO j=1, 2
        IF( j == 1 ) THEN
          Bmesh => BMesh1
        ELSE
          BMesh => BMesh2
        END IF
        
        MaxDiff = 0.0
        MinDiff = HUGE( MinDiff ) 
        
        DO i = 1, MeshDim
          IF( i == 1 ) THEN
            Coord => BMesh % Nodes % x 
          ELSE IF( i == 2 ) THEN
            Coord => Bmesh % Nodes % y
          ELSE
            Coord => Bmesh % Nodes % z
          END IF

          Diff = MAXVAL( Coord ) - MINVAL( Coord )
          MaxDiff = MAX( Diff, MaxDiff ) 
          IF( Diff < MinDiff ) THEN
            MinDiff = Diff
            MinDiffI = i
          END IF
        END DO

        RelDiff = MinDiff / MaxDiff
        IF( j == 1 ) THEN
          FlatDim = MinDiffI
          RelDiff1 = RelDiff 
        ELSE IF( j == 2 ) THEN
          IF( RelDiff < RelDiff1 ) FlatDim = MinDiffI
        END IF
      END DO

      CALL Info('FlatInterfaceMeshes','> Flat Projector Coordinate < set to: '//TRIM(I2S(FlatDim)))
      CALL ListAddInteger( BParams,'Flat Projector Coordinate',FlatDim )
    END IF


    DO j=1,2
      ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
      IF( j == 1 ) THEN
        Bmesh => BMesh1
      ELSE
        BMesh => BMesh2
      END IF

      ! Set the 3rd component to be the "distance" in the flat interface      
      IF( FlatDim == 3 ) THEN
        CONTINUE
      ELSE IF( FlatDim == 2 ) THEN
        Coord => BMesh % Nodes % y
        BMesh % Nodes % y => BMesh % Nodes % z
        BMesh % Nodes % z => Coord
        IF( MeshDim == 2 ) BMesh % Nodes % y = 0.0_dp
      ELSE IF( FlatDim == 1 ) THEN
        Coord => BMesh % Nodes % x
        BMesh % Nodes % x => BMesh % Nodes % y
        BMesh % Nodes % y => BMesh % Nodes % z
        Bmesh % Nodes % z => Coord
        IF( MeshDim == 2 ) BMesh % Nodes % y = 0.0_dp
      END IF

      IF( ReduceDim ) BMesh % Nodes % z = 0.0_dp

      Bmesh % MeshDim = 2
    END DO

  END SUBROUTINE FlatInterfaceMeshes
!------------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !> Given two interface meshes flatten them into the plane that 
  !> best fits either of the meshes. 
  !---------------------------------------------------------------------------
  SUBROUTINE PlaneInterfaceMeshes(BMesh1, BMesh2, BParams )
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Bmesh
    INTEGER :: i, j, n, nip, MeshDim
    REAL(KIND=dp) :: Normal(3), NormalSum(3), RefSum, Length, Planeness, &
        PlaneNormal(3,1), PlaneNormal1(3,1), Planeness1, Normal1(3), &
        Tangent(3), Tangent2(3), Coord(3), detJ, Normal0(3)
    REAL(KIND=dp), POINTER :: PNormal(:,:), Basis(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found, Stat, Normal0Set

    CALL Info('PlaneInterfaceMeshes','Flattening interface meshes to a plane',Level=8)    

    MeshDim = CurrentModel % DIMENSION
    PNormal => ListGetConstRealArray( BParams,'Plane Projector Normal',Found) 

    ! If the projector normal is not given determine it first 
    IF(.NOT. Found ) THEN     
      CALL Info('PlaneInterfaceMeshes','Could not find > Plane Projector Normal < so determining it now',Level=12)    

      n = MAX_ELEMENT_NODES
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), Basis(n) )
      ElementNodes % x = 0; ElementNodes % y = 0; ElementNodes % z = 0

      ! Fit a plane to both datasets
      DO j=1, 2
        IF( j == 1 ) THEN
          Bmesh => BMesh1
        ELSE      
          BMesh => BMesh2
        END IF

        NormalSum = 0.0_dp
        RefSum = 0.0_dp
        Normal0Set = .FALSE.

        ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
        !-------------------------------------------------------------------------
        DO i=1, BMesh % NumberOfBulkElements
          Element => BMesh % Elements(i)
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          IP = GaussPoints( Element ) 

          ElementNodes % x(1:n) = BMesh % Nodes % x(NodeIndexes(1:n))
          ElementNodes % y(1:n) = BMesh % Nodes % y(NodeIndexes(1:n))
          ElementNodes % z(1:n) = BMesh % Nodes % z(NodeIndexes(1:n))           

          DO nip=1, IP % n 
            stat = ElementInfo( Element,ElementNodes,&
                IP % u(nip),IP % v(nip),IP % w(nip),detJ,Basis)

            Normal = NormalVector( Element, ElementNodes, &
                IP % u(nip), IP % v(nip), .FALSE. ) 
            IF( .NOT. Normal0Set ) THEN
              Normal0 = Normal
              Normal0Set = .TRUE.
            END IF

            IF( SUM( Normal * Normal0 ) < 0.0 ) Normal = -Normal

            NormalSum = NormalSum + IP % S(nip) * DetJ * Normal
            RefSum = RefSum + IP % S(nip) * DetJ
          END DO
        END DO

        ! Normalize the normal to unity length
        Length = SQRT( SUM( NormalSum ** 2 ) )
        PlaneNormal(:,1) = NormalSum / Length

        ! Planeness is one if all the normals have the same direction
        Planeness = Length / RefSum 
        
        ! Save the key parameters of the first mesh
        IF( j == 1 ) THEN          
          PlaneNormal1 = PlaneNormal
          Planeness1 = Planeness
        END IF
      END DO

      ! Choose the mesh for which is close to a plane 
      IF( Planeness1 > Planeness ) THEN
        PRINT *,'PlaneNormal: Selecting slave normal'
        PlaneNormal = PlaneNormal1
      ELSE
        PRINT *,'PlaneNormal: Selecting master normal'        
        PlaneNormal = -PlaneNormal
      END IF

      PRINT *,'PlaneNormal selected:',PlaneNormal(:,1)

      CALL ListAddConstRealArray( BParams,'Plane Projector Normal',&
          3,1,PlaneNormal )
      DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, Basis )

      PNormal => ListGetConstRealArray( BParams,'Plane Projector Normal',Found) 
    END IF

    Normal = Pnormal(1:3,1)
    CALL TangentDirections( Normal, Tangent, Tangent2 )

    IF(.FALSE.) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent
      PRINT *,'Tangent2:',Tangent2
    END IF

    DO j=1,2
      IF( j == 1 ) THEN
        Bmesh => BMesh1
      ELSE
        BMesh => BMesh2
      END IF

      DO i=1,BMesh % NumberOfNodes
        Coord(1) = BMesh % Nodes % x(i)
        Coord(2) = BMesh % Nodes % y(i)
        Coord(3) = BMesh % Nodes % z(i)

        BMesh % Nodes % x(i) = SUM( Coord * Tangent )
        IF( MeshDim == 3 ) THEN
          BMesh % Nodes % y(i) = SUM( Coord * Tangent2 )
        ELSE
          BMesh % Nodes % y(i) = 0.0_dp
        END IF
        BMesh % Nodes % z(i) = SUM( Coord * Normal ) 
      END DO

      IF(.FALSE.) THEN
        PRINT *,'Range for mesh:',j
        PRINT *,'X:',MINVAL(BMesh % Nodes % x),MAXVAL(BMesh % Nodes % x)
        PRINT *,'Y:',MINVAL(BMesh % Nodes % y),MAXVAL(BMesh % Nodes % y)
        PRINT *,'Z:',MINVAL(BMesh % Nodes % z),MAXVAL(BMesh % Nodes % z)
      END IF
    END DO

    Bmesh % MeshDim = 2

  END SUBROUTINE PlaneInterfaceMeshes
  !------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given a permutation map the (x,y,z) such that the projector can better 
  !> be applied. E.g. if boundary has constant x, take that as the last coordinate.
  !---------------------------------------------------------------------------
  SUBROUTINE MapInterfaceCoordinate(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    LOGICAL :: Found
    REAL(KIND=dp), POINTER :: NodesX(:), NodesY(:), NodesZ(:), Wrk(:,:)
    INTEGER, POINTER :: CoordMap(:)
    INTEGER :: MeshNo
    TYPE(Mesh_t), POINTER :: BMesh
    
    ! Perform coordinate mapping
    !------------------------------------------------------------
    CoordMap => ListGetIntegerArray( BParams, & 
        'Projector Coordinate Mapping',Found )
    IF( .NOT. Found ) RETURN

    CALL Info('MapInterfaceCoordinates','Performing coordinate mapping',Level=8)
    
    IF ( SIZE( CoordMap ) /= 3 ) THEN
      WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
      CALL Error( 'MapInterfaceCoordinates', Message )
      WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
      CALL Fatal( 'MapInterfaceCoordinates', Message )
    END IF
    
    IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
      WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
      CALL Error( 'MapInterfaceCoordinates', Message )
      WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
      CALL Fatal( 'MapInterfaceCoordinates', Message )
    END IF

    DO MeshNo = 1,2
      IF( MeshNo == 1 ) THEN
        BMesh => BMesh1
      ELSE
        BMesh => BMesh2 
      END IF

      IF( CoordMap(1) == 1 ) THEN
        NodesX => BMesh % Nodes % x
      ELSE IF( CoordMap(1) == 2 ) THEN
        NodesX => BMesh % Nodes % y
      ELSE
        NodesX => BMesh % Nodes % z
      END IF
    
      IF( CoordMap(2) == 1 ) THEN
        NodesY => BMesh % Nodes % x
      ELSE IF( CoordMap(2) == 2 ) THEN
        NodesY => BMesh % Nodes % y
      ELSE
        NodesY => BMesh % Nodes % z
      END IF
      
      IF( CoordMap(3) == 1 ) THEN
        NodesZ => BMesh % Nodes % x
      ELSE IF( CoordMap(3) == 2 ) THEN
        NodesZ => BMesh % Nodes % y
      ELSE
        NodesZ => BMesh % Nodes % z
      END IF

      BMesh % Nodes % x => NodesX
      BMesh % Nodes % y => NodesY
      BMesh % Nodes % z => NodesZ
    END DO

  END SUBROUTINE MapInterfaceCoordinate


  ! Save projector, mainly a utility for debugging purposes
  !--------------------------------------------------------
  SUBROUTINE SaveProjector(Projector,SaveRowSum,Prefix,InvPerm,Parallel)
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: SaveRowSum 
    CHARACTER(LEN=*) :: Prefix
    INTEGER, POINTER, OPTIONAL :: InvPerm(:)
    LOGICAL, OPTIONAL :: Parallel

    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    INTEGER :: i,j,ii,jj
    REAL(KIND=dp) :: rowsum, dia, val
    INTEGER, POINTER :: IntInvPerm(:)
    LOGICAL :: GlobalInds

    IF(.NOT.ASSOCIATED(Projector)) RETURN
    
    IF( PRESENT( InvPerm ) ) THEN
      IntInvPerm => InvPerm 
    ELSE
      IntInvPerm => Projector % InvPerm
    END IF

    GlobalInds = .FALSE.
    IF(ParEnv % PEs == 1 ) THEN
      FileName = TRIM(Prefix)//'.dat'
    ELSE
      FileName = TRIM(Prefix)//'_part'//&
          TRIM(I2S(ParEnv % MyPe))//'.dat'
      IF( PRESENT( Parallel ) ) GlobalInds = Parallel
    END IF

    OPEN(1,FILE=FileName,STATUS='Unknown')    
    DO i=1,projector % numberofrows
      ii = intinvperm(i)
      IF( ii == 0 ) CYCLE
      IF( GlobalInds ) THEN
        ii = CurrentModel % Mesh % ParallelInfo % GlobalDofs(ii)
      END IF
      DO j=projector % rows(i), projector % rows(i+1)-1
        jj = projector % cols(j)
        val = projector % values(j)
        IF( GlobalInds ) THEN
          jj = CurrentModel % Mesh % ParallelInfo % GlobalDofs(jj)
          WRITE(1,*) ii,jj,ParEnv % MyPe, val
        ELSE
          WRITE(1,*) ii,jj,val
        END IF
      END DO
    END DO
    CLOSE(1)     

    IF( SaveRowSum ) THEN
      IF(ParEnv % PEs == 1 ) THEN
        FileName = TRIM(Prefix)//'_rsum.dat'
      ELSE
        FileName = TRIM(Prefix)//'_rsum_part'//&
            TRIM(I2S(ParEnv % MyPe))//'.dat'
      END IF
      
      OPEN(1,FILE=FileName,STATUS='Unknown')
      DO i=1,projector % numberofrows
        ii = intinvperm(i)
        IF( ii == 0 ) CYCLE
        rowsum = 0.0_dp
        dia = 0.0_dp

        DO j=projector % rows(i), projector % rows(i+1)-1          
          jj = projector % cols(j)
          val = projector % values(j)
          IF( ii == jj ) THEN
            dia = val
          END IF
          rowsum = rowsum + val
        END DO

        IF( GlobalInds ) THEN
          ii = CurrentModel % Mesh % ParallelInfo % GlobalDofs(ii)
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i), ParEnv % MyPe, dia, rowsum
        ELSE
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i),dia, rowsum
        END IF

      END DO
      CLOSE(1)     
    END IF

  END SUBROUTINE SaveProjector



  ! Set projector abs(rowsum) to unity
  !--------------------------------------------------------
  SUBROUTINE SetProjectorRowsum( Projector )
    TYPE(Matrix_t), POINTER :: Projector

    INTEGER :: i,j
    REAL(KIND=dp) :: rowsum

    DO i=1,projector % numberofrows
      rowsum = 0.0_dp
      DO j=projector % rows(i), projector % rows(i+1)-1
        rowsum = rowsum + ABS( projector % values(j) )
      END DO
      DO j=projector % rows(i), projector % rows(i+1)-1
        projector % values(j) = projector % values(j) / rowsum
      END DO
    END DO

  END SUBROUTINE SetProjectorRowsum


!------------------------------------------------------------------------------
!> Create a projector between Master and Target boundaries.
!> The projector may be a nodal projector x=Px or a weigted 
!> Galerking projector such that Qx=Px. In the first case the projector 
!> will be P and in the second case [Q-P]. 
!------------------------------------------------------------------------------
  FUNCTION PeriodicProjector( Model, Mesh, This, Trgt, cdim, &
      Galerkin ) RESULT(Projector)
!------------------------------------------------------------------------------   
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    INTEGER, OPTIONAL :: cdim
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL, OPTIONAL :: Galerkin
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,dim
    LOGICAL :: GotIt, UseQuadrantTree, Success, IntGalerkin, &
        Rotational, AntiRotational, Sliding, AntiSliding, Repeating, AntiRepeating, &
        Discontinuous, NodalJump, Radial, AntiRadial, DoNodes, DoEdges, &
        Flat, Plane, LevelProj, FullCircle, Cylindrical, UseExtProjector, &
        ParallelNumbering, EnforceOverlay
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    TYPE(Nodes_t), POINTER :: MeshNodes, GaussNodes
    REAL(KIND=dp) :: NodeScale, EdgeScale, Radius, Coeff
    TYPE(ValueList_t), POINTER :: BC

    INTERFACE
      FUNCTION WeightedProjector(BMesh2, BMesh1, InvPerm2, InvPerm1, &
          UseQuadrantTree, Repeating, AntiRepeating, PeriodicScale, &
          NodalJump ) &
         RESULT ( Projector )
        USE Types
        TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
        REAL(KIND=dp) :: PeriodicScale
        INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
        LOGICAL :: UseQuadrantTree, Repeating, AntiRepeating
        TYPE(Matrix_t), POINTER :: Projector
        LOGICAL :: NodalJump
      END FUNCTION WeightedProjector



    END INTERFACE
!------------------------------------------------------------------------------
    Projector => NULL()
    IF ( This <= 0  ) RETURN    


    DIM = CoordinateSystemDimension()

    CALL ResetTimer('PeriodicProjector')
    
    Projector => NULL()
    BC => Model % BCs(This) % Values
    PMesh => Mesh

    
    ! Whether to choose nodal or Galerkin projector is determined by an optional
    ! flag. The default is the nodal projector.
    !--------------------------------------------------------------------------
    IF( PRESENT( Galerkin) ) THEN
      IntGalerkin = Galerkin
    ELSE
      IntGalerkin = ListGetLogical( BC, 'Galerkin Projector', GotIt )
    END IF

    ! If the boundary is discontinuous then we have the luxury of creating the projector
    ! very cheaply using the permutation vector. This does not need the target as the 
    ! boundary is self-contained.
    !------------------------------------------------------------------------------------
    IF( ListGetLogical( BC, 'Discontinuous Boundary', GotIt ) .AND. Mesh % DisContMesh )THEN
      IF( IntGalerkin ) THEN
        Projector => WeightedProjectorDiscont( PMesh, This )
      ELSE
        Projector => NodalProjectorDiscont( PMesh, This )
      END IF
      
      IF ( .NOT. ASSOCIATED( Projector ) ) RETURN
      GOTO 100
    END IF
    
    IF ( Trgt <= 0 ) RETURN    

    ! Create the mesh projector, and if needed, also eliminate the ghost nodes
    ! There are two choices of projector: a nodal projector P in x=Px, and a 
    ! Galerkin projector [Q-P] in Qx=Px. 
    ! The projector is assumed to be either a rotational projector with no translation
    ! and rotation, or then generic one with possible coordinate mapping.
    !---------------------------------------------------------------------------------
    CALL Info('PeriodicProjector','-----------------------------------------------------',Level=8)
    WRITE( Message,'(A,I0,A,I0)') 'Creating projector between BCs ',This,' and ',Trgt
    CALL Info('PeriodicProjector',Message,Level=8)

    ! Create temporal mesh structures that are utilized when making the 
    ! projector between "This" and "Trgt" boundary.
    !--------------------------------------------------------------------------
    BMesh1 => AllocateMesh()
    BMesh2 => AllocateMesh()
    
    CALL CreateInterfaceMeshes( Model, Mesh, This, Trgt, Bmesh1, BMesh2, &
        Success ) 

    IF(.NOT. Success) THEN
      CALL ReleaseMesh(BMesh1)
      CALL ReleaseMesh(BMesh2)
      RETURN
    END IF

    ! Do we have external procedure to take care of the projection matrix creation
    UseExtProjector = ListGetLogical( BC, 'External Projector', GotIt )

    ! If requested map the interface coordinate from (x,y,z) to any permutation
    ! of these. 
    CALL MapInterfaceCoordinate( BMesh1, BMesh2, Model % BCs(This) % Values )

    ! Check whether to use (anti)rotational projector.
    ! We don't really know on which side the projector was called so 
    ! let's check both sides.
    !--------------------------------------------------------------------------
    Rotational = ListGetLogical( BC,&
        'Rotational Projector',GotIt )
    AntiRotational = ListGetLogical( BC,&
        'Anti Rotational Projector',GotIt )
    IF( AntiRotational ) Rotational = .TRUE.

    Cylindrical =  ListGetLogical( BC,&
        'Cylindrical Projector',GotIt )

    Radial = ListGetLogical( BC,&
        'Radial Projector',GotIt )
    AntiRadial = ListGetLogical( BC,&
        'Anti Radial Projector',GotIt )
    IF( AntiRadial ) Radial = .TRUE.

    Sliding = ListGetLogical( BC, &
        'Sliding Projector',GotIt )
    AntiSliding = ListGetLogical( BC, &
        'Anti Sliding Projector',GotIt )
    IF( AntiSliding ) Sliding = .TRUE. 

    Flat = ListGetLogical( BC, &
        'Flat Projector',GotIt )
    Plane = ListGetLogical( BC, &
        'Plane Projector',GotIt )

    IF( Radial ) CALL Info('PeriodicProjector','Enforcing > Radial Projector <',Level=12)
    IF( Sliding ) CALL Info('PeriodicProjector','Enforcing > Sliding Projector <',Level=12)
    IF( Cylindrical ) CALL Info('PeriodicProjector','Enforcing > Cylindrical Projector <',Level=12)
    IF( Rotational ) CALL Info('PeriodicProjector','Enforcing > Rotational Projector <',Level=12)
    IF( Flat ) CALL Info('PeriodicProjector','Enforcing > Flat Projector <',Level=12)
    IF( Plane ) CALL Info('PeriodicProjector','Enforcing > Plane Projector <',Level=12)


    NodeScale = ListGetConstReal( BC, 'Mortar BC Scaling',GotIt)
    IF(.NOT.Gotit ) THEN
      IF( AntiRadial ) THEN
        NodeScale = -1._dp
      ELSE
        NodeScale = 1.0_dp
      END IF
    END IF
    EdgeScale = NodeScale

    NodalJump = ListCheckPrefix( BC,'Mortar BC Coefficient')
    IF(.NOT. NodalJump ) THEN
      NodalJump = ListCheckPrefix( BC,'Mortar BC Resistivity')
    END IF

    ! There are tailored projectors for simplified interfaces
    !-------------------------------------------------------------

    ! Stride projector is obsolite and has been eliminated. 
    IF( ListGetLogical( BC,'Stride Projector',GotIt) ) THEN
      CALL ListAddLogical( BC,'Level Projector',.TRUE.)
      CALL ListAddLogical( BC,'Level Projector Strong',.TRUE.)
      CALL Warn('PeriodicProjector','Enforcing > Level Projector < instead of old > Stride Projector <')
    END IF

    LevelProj = ListGetLogical( BC,'Level Projector',GotIt) 
    IF( Rotational .OR. Cylindrical .OR. Radial .OR. Flat .OR. Plane ) THEN
      IF(.NOT. GotIt ) THEN
        CALL Info('PeriodicProjector','Enforcing > Level Projector = True < with dimensional reduction',&
            Level = 7 )
        LevelProj = .TRUE. 
      ELSE IF(.NOT. LevelProj ) THEN
        ! If we have dimensionally reduced projector but don't use LevelProjector 
        ! to integrate over it, then ensure that the 3rd coordinate is set to zero.
        BMesh1 % Nodes % z = 0.0_dp
        BMesh2 % Nodes % z = 0.0_dp
      END IF
    END IF


    IF( LevelProj ) THEN
      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',GotIt ) ) THEN
        DoNodes = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Nodes',GotIt) ) THEN
          DoNodes = .FALSE.
        ELSE
          DoNodes = ( Mesh % NumberOfNodes > 0 ) 
        END IF
      END IF

      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',GotIt ) ) THEN
        DoEdges = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Edges',GotIt) ) THEN
          DoEdges = .FALSE.
        ELSE
          ! We are conservative here since there may be edges in 2D which 
          ! still cannot be used for creating the projector
          DoEdges = ( Mesh % NumberOfEdges > 0 .AND. &
              Mesh % MeshDim == 3 .AND. Dim == 3 )

          ! Ensure that there is no p-elements that made us think that we have edges
          ! Here we assume that if there is any p-element then also the 1st element is such
          IF( DoEdges ) THEN
            IF(isPelement(Mesh % Elements(1))) THEN
              DoEdges = .FALSE.
              CALL Info('PeriodicProjector','Edge projector will not be created for p-element mesh',Level=10)
            END IF
          END IF
        END IF
      END IF
    END IF


    ! If the interface is rotational move to (phi,z) plane and alter the phi coordinate 
    ! so that the meshes coinside. 
    ! Otherwise make the two meshes to coinside using rotation, translation &
    ! scaling. 
    !---------------------------------------------------------------------------------
    Radius = 1.0_dp
    FullCircle = .FALSE.
    EnforceOverlay = ListGetLogical( BC, 'Mortar BC enforce overlay', GotIt )

    IF( Rotational .OR. Cylindrical ) THEN
      CALL RotationalInterfaceMeshes( BMesh1, BMesh2, BC, Cylindrical, &
          Radius, FullCircle )
    ELSE IF( Radial ) THEN
      CALL RadialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Flat ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Plane ) THEN
      CALL PlaneInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( .NOT. Sliding ) THEN
      IF( .NOT. GotIt ) EnforceOverlay = .TRUE.
    END IF

    IF( EnforceOverlay ) THEN
      CALL OverlayIntefaceMeshes( BMesh1, BMesh2, BC )
    END IF

    Repeating = ( Rotational .AND. .NOT. FullCircle ) .OR. Sliding 
    AntiRepeating = ( AntiRotational .AND. .NOT. FullCircle ) .OR. AntiSliding 

    IF( UseExtProjector ) THEN
      Projector => ExtProjectorCaller( PMesh, BMesh1, BMesh2, This )
    ELSE IF( LevelProj ) THEN 
      Projector => LevelProjector( BMesh1, BMesh2, Repeating, AntiRepeating, &
          FullCircle, Radius, DoNodes, DoEdges, &          
          NodeScale, EdgeScale, BC )
    ELSE 
      IF( FullCircle ) THEN
        CALL Fatal('PeriodicProjector','A full circle cannot be dealt with the generic projector!')
      END IF

      UseQuadrantTree = ListGetLogical(Model % Simulation,'Use Quadrant Tree',GotIt)
      IF( .NOT. GotIt ) UseQuadrantTree = .TRUE.
      IF( IntGalerkin ) THEN
        Projector => WeightedProjector( BMesh2, BMesh1, BMesh2 % InvPerm, BMesh1 % InvPerm, &
            UseQuadrantTree, Repeating, AntiRepeating, NodeScale, NodalJump )
      ELSE
        Projector => NodalProjector( BMesh2, BMesh1, &
            UseQuadrantTree, Repeating, AntiRepeating )
      END IF
    END IF


    ! Deallocate mesh structures:
    !---------------------------------------------------------------
    BMesh1 % Projector => NULL()
    BMesh1 % Parent => NULL()
    !DEALLOCATE( BMesh1 % InvPerm ) 
    CALL ReleaseMesh(BMesh1)

    BMesh2 % Projector => NULL()
    BMesh2 % Parent => NULL()
    !DEALLOCATE( BMesh2 % InvPerm ) 
    CALL ReleaseMesh(BMesh2)

100 Projector % ProjectorBC = This

    IF( ListGetLogical( BC,'Projector Set Rowsum',GotIt ) ) THEN
      CALL SetProjectorRowsum( Projector )
    END IF

    Coeff = ListGetConstReal( BC,'Projector Multiplier',GotIt) 
    IF(.NOT. GotIt) Coeff = ListGetConstReal( Model % Simulation,&
        'Projector Multiplier',GotIt) 
    IF( GotIt ) Projector % Values = Coeff * Projector % Values

    IF( ListGetLogical( BC,'Save Projector',GotIt ) ) THEN
      ParallelNumbering = ListGetLogical( BC,'Save Projector Global Numbering',GotIt )

      CALL SaveProjector( Projector, .TRUE.,'p'//TRIM(I2S(This)), Parallel = ParallelNumbering) 
      ! Dual projector if it exists
      IF( ASSOCIATED( Projector % Ematrix ) ) THEN
        CALL SaveProjector( Projector % Ematrix, .TRUE.,'pd'//TRIM(I2S(This)), &
            Projector % InvPerm, Parallel = ParallelNumbering) 
      END IF

      ! Biorthogonal projector if it exists
      IF( ASSOCIATED( Projector % Child ) ) THEN
        CALL SaveProjector( Projector % Child, .TRUE.,'pb'//TRIM(I2S(This)), & 
            Projector % InvPerm, Parallel = ParallelNumbering ) 
      END IF

      IF( ListGetLogical( BC,'Save Projector And Stop',GotIt ) ) STOP
    END IF    

    CALL CheckTimer('PeriodicProjector',Delete=.TRUE.)
    CALL Info('PeriodicProjector','Projector created, now exiting...',Level=8)

!------------------------------------------------------------------------------
  END FUNCTION PeriodicProjector
!------------------------------------------------------------------------------



  FUNCTION ExtProjectorCaller( Mesh, SlaveMesh, MasterMesh, SlaveBcInd ) RESULT ( Projector )
    !---------------------------------------------------------------------------
    USE Lists
    USE Messages
    USE Types
    USE GeneralUtils
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh, SlaveMesh, MasterMesh
    INTEGER :: SlaveBCind
    TYPE(Matrix_t), POINTER :: Projector    
    !--------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found, Parallel, CreateDual, BiorthogonalBasis
    INTEGER(KIND=AddrInt) :: ProjectorAddr

    CALL Info('ExtProjectorCaller','Creating projector using an external function',Level=7)

    Parallel = ( ParEnv % PEs > 1 )
    Mesh => CurrentModel % Mesh
    BC => CurrentModel % BCs(SlaveBcInd) % Values

    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN

    CreateDual = ListGetLogical( BC,'Create Dual Projector',Found ) 
    IF( CreateDual ) THEN
      Projector % Ematrix => AllocateMatrix()
      Projector % Ematrix % FORMAT = MATRIX_LIST
      Projector % Ematrix % ProjectorType = PROJECTOR_TYPE_GALERKIN
    ELSE
      Projector % EMatrix => NULL()
    END IF

    ! Check whether biorthogonal basis for projectors requested:
    ! ----------------------------------------------------------
    BiOrthogonalBasis = ListGetLogical( BC, 'Use Biorthogonal Basis', Found)

    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      BiOrthogonalBasis = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( BiOrthogonalBasis ) THEN
        CALL Info('ContactProjector',&
            'Setting > Use Biorthogonal Basis < to True to enable elimination',Level=8)
      END IF
    END IF

    IF (BiOrthogonalBasis) THEN
      Projector % Child => AllocateMatrix()
      Projector % Child % Format = MATRIX_LIST
      CALL Info('ContactProjector','Using biorthogonal basis, as requested',Level=8)      
    ELSE
      Projector % Child => NULL()
    END IF


    ProjectorAddr = CurrentModel % Solver % MortarProc
    IF( ProjectorAddr == 0 ) THEN
      CALL Fatal('ExtProjectorCaller','External projector requested by no > Mortar Proc < given!')
    ELSE
      CALL ExecMortarProjector( ProjectorAddr, &
          Mesh, SlaveMesh, MasterMesh, SlaveBCind, Projector )
    END IF

    ! Now change the matrix format to CRS from list matrix
    !--------------------------------------------------------------
    CALL List_toCRSMatrix(Projector)
    CALL CRS_SortMatrix(Projector,.TRUE.)

    IF( ASSOCIATED(Projector % Child) ) THEN
      CALL List_toCRSMatrix(Projector % Child)
      CALL CRS_SortMatrix(Projector % Child,.TRUE.)
    END IF

    IF( ASSOCIATED( Projector % Ematrix) ) THEN
      CALL List_toCRSMatrix(Projector % Ematrix)
      CALL CRS_SortMatrix(Projector % Ematrix,.TRUE.)
    END IF

    CALL Info('ExtProjectorCaller','Projector created',Level=10)
    
  END FUNCTION ExtProjectorCaller


!------------------------------------------------------------------------------
!> Create node distribution for a unit segment x \in [0,1] with n elements 
!> i.e. n+1 nodes. There are different options for the type of distribution.
!> 1) Even distribution 
!> 2) Geometric distribution
!> 3) Arbitrary distribution determined by a functional dependence
!> Note that the 3rd algorithm involves iterative solution of the nodal
!> positions and is therefore not bullet-proof.
!------------------------------------------------------------------------------
  SUBROUTINE UnitSegmentDivision( w, n )
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    !---------------------------------------------------------------
    INTEGER :: i,J,iter,maxiter
    REAL(KIND=dp) :: q,h1,hn,minhn,err_eps,err,xn
    REAL(KIND=dp), ALLOCATABLE :: wold(:),h(:)
    LOGICAL :: Found, GotRatio
    TYPE(Nodes_t) :: Nodes
    
    ! Linear distribution and an initial guess for the generic case
    !---------------------------------------------------------------

    ! Geometric division
    !---------------------------------------------------------------
    q = ListGetConstReal(CurrentModel % Simulation,'Extruded Mesh Ratio',GotRatio)
    IF( GotRatio ) THEN
      CALL Info('UnitSegmentDivision','Creating geometric division',Level=5)

      h1 = (1-q**(1.0_dp/n))/(1-q)
      w(0) = 0.0_dp
      hn = h1;
      DO i=1,n-1
        w(i) = w(i-1) + hn;
        hn = hn * ( q**(1.0_dp/n) )
      END DO
      w(n) = 1.0_dp


    ! Generic division given by a function
    !-----------------------------------------------------------------------
    ELSE IF( ListCheckPresent( CurrentModel % Simulation,'Extruded Mesh Density') ) THEN

      CALL Info('UnitSegmentDivision','Creating functional division',Level=5)

      ! Initial guess is an even distribtion
      DO i=0,n     
        w(i) = i/(1._dp * n)
      END DO

      ALLOCATE( wold(0:n),h(1:n))
      wold = w

      ! paramaters that determine the accuracy of the iteration
      maxiter = 10000
      err_eps = 1.0e-6

      ! Iterate to have a density distribution
      !---------------------------------------
      DO iter=1,maxiter
        
        minhn = HUGE(minhn)
        wold = w

        ! Compute the point in the local mesh xn \in [0,1]  
        ! and get the mesh parameter for that element from
        ! external function.
        !---------------------------------------------------
        DO i=1,n
          xn = (w(i)+w(i-1))/2.0_dp
          minhn = MIN( minhn, w(i)-w(i-1) )
          h(i) = ListGetFun( CurrentModel % Simulation,'Extruded Mesh Density', xn )
          IF( h(i) < EPSILON( h(i) ) ) THEN
            CALL Fatal('UnitSegmentDivision','Given value for h(i) was negative!')
          END IF
        END DO

        ! Utilize symmetric Gauss-Seidel to compute the new positions, w(i).
        ! from a weigted mean of the desired elemental densities, h(i).
        ! Note that something more clever could be applied here. 
        ! This was just a first implementation...
        !-------------------------------------------------------------
        DO i=1,n-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        DO i=n-1,1,-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        
        ! If the maximum error is small compared to the minumum elementsize then exit
        !-----------------------------------------------------------------------------
        err = MAXVAL( ABS(w-wold))/minhn

        IF( err < err_eps ) THEN
          WRITE( Message, '(A,I0,A)') 'Convergence obtained in ',iter,' iterations'
          CALL Info('UnitSegmentDivision', Message, Level=9 )
          EXIT
        END IF
      END DO

      IF( iter > maxiter ) THEN
        CALL Warn('UnitSegmentDivision','No convergence obtained for the unit mesh division!')
      END IF

    ! Uniform division 
    !--------------------------------------------------------------
    ELSE
      CALL Info('UnitSegmentDivision','Creating linear division',Level=5)
      DO i=0,n     
        w(i) = i/(1._dp * n)
      END DO
    END IF

    CALL Info('UnitSegmentDivision','Mesh division ready',Level=9)
    DO i=0,n
      WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',w(i)
      CALL Info('UnitSegmentDivision', Message, Level=9 )
    END DO


  END SUBROUTINE UnitSegmentDivision
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Given a 2D mesh extrude it to be 3D. The 3rd coordinate will always
!> be at the interval [0,1]. Therefore the adaptation for different shapes
!> must be done with StructuredMeshMapper, or some similar utility. 
!> The top and bottom surface will be assigned Boundary Condition tags
!> with indexes one larger than the maximum used on by the 2D mesh. 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrude(Mesh_in, in_levels, ExtrudedMeshName) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    INTEGER :: in_levels
    CHARACTER(LEN=MAX_NAME_LEN),INTENT(IN),OPTIONAL :: ExtrudedMeshName

!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n,cnt,cnt101,ind(8),max_baseline_bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    INTEGER :: nnodes,gnodes,gelements,ierr
    LOGICAL :: isParallel, Found, NeedEdges, PreserveBaseline
    REAL(KIND=dp)::w,MinCoord,MaxCoord,CurrCoord
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
!------------------------------------------------------------------------------
    Mesh_out => AllocateMesh()
!   Mesh_out = Mesh_in


    isParallel = ParEnv % PEs>1

    ! Generate volume nodal points:
    ! -----------------------------
    n=Mesh_in % NumberOfNodes
    nnodes=(in_levels+2)*n

    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )

    gelements = Mesh_in % NumberOfBulkElements
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo

      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % INTERFACE(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      j=0
      DO i=1,Mesh_in % NumberOfNodes
        IF (PI_in % NeighbourList(i) % &
            Neighbours(1) == ParEnv % MyPE ) j=j+1
      END DO
      CALL MPI_ALLREDUCE(j,gnodes,1, &
           MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

      j=0
      DO i=1,Mesh_in % NumberOfBulkElements
        IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
      END DO
      CALL MPI_ALLREDUCE(j,gelements,1, &
           MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    END IF


    ! Create the division for the 1D unit mesh
    !--------------------------------------------
    ALLOCATE( Wtable( 0: in_levels + 1 ) )
    CALL UnitSegmentDivision( Wtable, in_levels + 1 ) 

    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = 3 

    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF


    PreserveBaseline = ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found )
    IF(.NOT. Found) PreserveBaseline = .FALSE.

    MinCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Min Coordinate',Found )
    IF(.NOT. Found) MinCoord = 0.0_dp

    MaxCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Max Coordinate',Found )
    IF(.NOT. Found) MaxCoord = 1.0_dp

    cnt=0
    DO i=0,in_levels+1

      w = Wtable( i ) 
      CurrCoord = w * MaxCoord + (1-w) * MinCoord

      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          PI_out % INTERFACE(cnt) = PI_in % INTERFACE(j)

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(&
               SIZE(PI_in % NeighbourList(j) % Neighbours)))
          PI_out % NeighbourList(cnt) % Neighbours = &
            PI_in % NeighbourList(j) % Neighbours

          PI_out % GlobalDOFs(cnt) = PI_in % GlobalDOFs(j)+i*gnodes
        END IF

      END DO
    END DO
    Mesh_out % NumberOfNodes=cnt

    ! Count 101 elements:
    ! (these require an extra layer)
    ! -------------------

    cnt101 = 0
    DO i=Mesh_in % NumberOfBulkElements+1, &
         Mesh_in % NumberOfBulkElements+Mesh_in % NumberOfBoundaryElements
       IF(Mesh_in % Elements(i) % TYPE % ElementCode == 101) cnt101 = cnt101+1
    END DO

    n=SIZE(Mesh_in % Elements)
    IF (PreserveBaseline) THEN
        ALLOCATE(Mesh_out % Elements(n*(in_levels+3) + Mesh_in % NumberOfBoundaryElements + cnt101) )
    ELSE
	ALLOCATE(Mesh_out % Elements(n*(in_levels+3) + cnt101) )
    END IF

    ! Generate volume bulk elements:
    ! ------------------------------

    Mesh_out % MaxElementNodes = 0

    NeedEdges=.FALSE.
    n=Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt=cnt+1
        Mesh_out % Elements(cnt) = Mesh_in % Elements(j)

        l_n=0
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+i*n
        END DO
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+(i+1)*n
        END DO
        Mesh_out % Elements(cnt) % NDOFs = l_n
        Mesh_out % MaxElementNodes=MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(706)
        CASE(8)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(808)
        END SELECT

        Mesh_out % Elements(cnt) % GElementIndex = &
             Mesh_in % Elements(j) % GelementIndex + gelements*i

        Mesh_out % Elements(cnt) % ElementIndex = cnt
        ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n)) 
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % NodeIndexes = ind(1:l_n)
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()

        k = Mesh_out % Elements(cnt) % DGDOFs
        IF(k>0) THEN
          Mesh_out % Elements(cnt) % DGDOFs = &
                Mesh_out % Elements(cnt) % TYPE % NumberOFNodes
          k = Mesh_out % Elements(cnt) % DGDOFs
          ALLOCATE(Mesh_out % Elements(cnt) % DGIndexes(k))
          DO l=1,k
            dg_n = dg_n + 1
            Mesh_out % Elements(cnt) % DGIndexes(l) = dg_n
          END DO
          NeedEdges=.TRUE.
        END IF

        IF(ASSOCIATED(Mesh_in % Elements(j) % PDefs)) THEN
          NeedEdges=.TRUE.
          ALLOCATE(Mesh_out % Elements(cnt) % PDefs)
          Mesh_out % Elements(cnt) % PDefs=Mesh_in % Elements(j) % PDefs
        END IF
      END DO
    END DO
    Mesh_out % NumberOfBulkElements=cnt

    max_bid=0
    max_baseline_bid=0

    ! -------------------------------------------------------
    IF (PreserveBaseline) THEN
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        max_bid = MAX(max_bid, Mesh_in % Elements(k) % &
                BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in %  NumberOfBulkElements*(in_levels+1)+ &
	                   (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*(in_levels+1)+ &
	      (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 2
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(2)) 
          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)
          Mesh_out % Elements(cnt) % NodeIndexes = ind(1:2)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(202)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    
      IF(isParallel) THEN
        j=max_bid
        CALL MPI_ALLREDUCE(j,max_bid,1, &
            MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      END IF

      max_baseline_bid = max_bid

    END IF


    ! Add side boundaries with the bottom mesh boundary id's:
    ! (or shift ids if preserving the baseline boundary)
    ! -------------------------------------------------------
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
           Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

        max_bid = MAX(max_bid, max_baseline_bid + &
           Mesh_in % Elements(k) % BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 4
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(4)) 

          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n
          ind(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
          ind(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          Mesh_out % Elements(cnt) % NodeIndexes = ind(1:4)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes+i*n
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF 
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    END DO

    !Take care of extra 101 elements
    !-------------------------------

    IF(cnt101 > 0) THEN
       DO j=1,Mesh_in % NumberOfBoundaryElements
          k = j + Mesh_in % NumberOfBulkElements

          IF(Mesh_in % Elements(k) % TYPE % ElementCode /= 101) CYCLE
          cnt=cnt+1

          Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

          ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
          Mesh_out % Elements(cnt) % BoundaryInfo = &
               Mesh_in % Elements(k) % BoundaryInfo

          Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
               Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

          max_bid = MAX(max_bid, max_baseline_bid + &
               Mesh_in % Elements(k) % BoundaryInfo % Constraint)

          Mesh_out % Elements(cnt) % NDOFs = 1
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(1))
          Mesh_out % Elements(cnt) % NodeIndexes = &
               Mesh_in % Elements(k) % NodeIndexes+(in_levels+1)*n
          Mesh_out % Elements(cnt) % TYPE => &
               Mesh_in % Elements(k) % TYPE

          Mesh_out % Elements(cnt) % ElementIndex = cnt
          Mesh_out % Elements(cnt) % DGDOFs = 0
          Mesh_out % Elements(cnt) % DGIndexes => NULL()
          Mesh_out % Elements(cnt) % PDefs => NULL()
          Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
          Mesh_out % Elements(cnt) % FaceIndexes => NULL()
          Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
       END DO
    END IF
    
    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'First Extruded BC set to: ',max_bid+1
    CALL Info('ExtrudeMesh',Message,Level=8)

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'Number of new BCs for layers: ',max_body
    CALL Info('ExtrudeMesh',Message,Level=8)


    ! Add bottom boundary:
    ! --------------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    ! Add top boundary:
    ! -----------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(in_levels*Mesh_in % NumberOfBulkElements+i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId + max_body
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes+(in_Levels+1)*n
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    Mesh_out % NumberOfBoundaryElements=cnt-Mesh_out % NumberOfBulkElements

    Mesh_out % Name=Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs  = Mesh_out % MaxElementNodes
    Mesh_out % MeshDim = 3
    CurrentModel % DIMENSION = 3

    IF ( NeedEdges ) CALL SetMeshEdgeFaceDOFs(Mesh_out)
    CALL SetMeshMaxDOFs(Mesh_out)

    IF (PRESENT(ExtrudedMeshName)) THEN
       CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
    END IF

!------------------------------------------------------------------------------
  END FUNCTION MeshExtrude
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Writes the mesh to disk. Note that this does not include the information
!> of shared nodes needed in parallel computation. This may be used for 
!> debugging purposes and for adaptive solution, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk( NewMesh, Path )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Path
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,MaxNodes,ElmCode,Parent1,Parent2
!------------------------------------------------------------------------------

    OPEN( 1,FILE=TRIM(Path) // '/mesh.header',STATUS='UNKNOWN' )
    WRITE( 1,'(3i8)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, NewMesh % NumberOfBoundaryElements
    
    WRITE( 1,* ) 2
    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(k) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(k) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(k) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(2i8)' ) ElmCode,NewMesh % NumberOfBoundaryElements

    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(i) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(i) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(2i8)' ) ElmCode,NewMesh % NumberOfBulkElements
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.nodes', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       WRITE(1,'(i6,a,3e23.15)',ADVANCE='NO') i,' -1 ', &
            NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.elements', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       WRITE(1,'(3i7)',ADVANCE='NO') i, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          WRITE(1,'(i7)', ADVANCE='NO') &
               NewMesh % Elements(i) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.boundary', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex
       WRITE(1,'(5i7)',ADVANCE='NO') i, &
            NewMesh % Elements(k) % BoundaryInfo % Constraint, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          WRITE(1,'(i7)', ADVANCE='NO') &
               NewMesh % Elements(k) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)
!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk2(Model, NewMesh, Path, Partition )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: NewMesh
    CHARACTER(LEN=*) :: Path
    INTEGER, OPTIONAL :: Partition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeList(100),ElmCodeCounts(100),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    INTEGER, POINTER :: BList(:)
    INTEGER, ALLOCATABLE :: ElementCodes(:)
    LOGICAL :: Parallel, WarnNoTarget, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: headerFN, elementFN, nodeFN,&
         boundFN, sharedFN
!------------------------------------------------------------------------------

    IF(PRESENT(Partition)) THEN
       Parallel = .TRUE.
       WRITE(headerFN, '(A,I0,A)') '/part.',Partition+1,'.header'
       WRITE(elementFN, '(A,I0,A)') '/part.',Partition+1,'.elements'
       WRITE(nodeFN, '(A,I0,A)') '/part.',Partition+1,'.nodes'
       WRITE(boundFN, '(A,I0,A)') '/part.',Partition+1,'.boundary'
       WRITE(sharedFN, '(A,I0,A)') '/part.',Partition+1,'.shared'
    ELSE
       Parallel = .FALSE.
       headerFN = '/mesh.header'
       elementFN = '/mesh.elements'
       nodeFN = '/mesh.nodes'
       boundFN = '/mesh.boundary'
    END IF

    !Info for header file

    ElmCodeList = 0 !init array
    NumElmCodes = 0
    NumElements = NewMesh % NumberOfBoundaryElements + &
         NewMesh % NumberOfBulkElements
    ALLOCATE(ElementCodes(NumElements))

    !cycle to bring element code list into array-inquirable form
    DO i=1,NumElements
       ElementCodes(i) = NewMesh % Elements(i) % TYPE % ElementCode
    END DO

    DO i=NumElements,1,-1 !this should give element codes increasing value, which appears to be
                          !'standard' though I doubt it matters
       IF(ANY(ElmCodeList == ElementCodes(i))) CYCLE
       NumElmCodes = NumElmCodes + 1
       ElmCodeList(NumElmCodes) = ElementCodes(i)
    END DO

    DO j=1,NumElmCodes
       ElmCodeCounts(j) = COUNT(ElementCodes == ElmCodeList(j))
    END DO

    !Write header file
    OPEN( 1,FILE=TRIM(Path) // headerFN,STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, &
         NewMesh % NumberOfBoundaryElements

    WRITE( 1,'(i0)' ) NumElmCodes
    DO j=1,NumElmCodes
       WRITE( 1,'(i0,x,i0,x)' ) ElmCodeList(j),ElmCodeCounts(j)
    END DO
    IF(Parallel) THEN !need number of shared nodes
       NoShared = 0
       DO i=1,NewMesh % NumberOfNodes
          IF(SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
               Neighbours) > 1) THEN
             NoShared = NoShared + 1
          END IF
       END DO
       WRITE( 1,'(i0,x,i0)') NoShared, 0
    END IF
    CLOSE(1)

    !Write nodes file
    OPEN( 1,FILE=TRIM(Path) // nodeFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       IF (Parallel) THEN
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % ParallelInfo % GlobalDOFs(i)
       ELSE
          WRITE(1,'(i0,x)', ADVANCE='NO') i
       END IF
       WRITE(1,'(a,x,ES17.10,x,ES17.10,x,ES17.10)',ADVANCE='NO') &
            ' -1 ', NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    !Write elements file
    OPEN( 1,FILE=TRIM(Path) // elementFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       IF(Parallel) THEN
          ElemID = NewMesh % Elements(i) % GElementIndex
       ELSE
          ElemID = i
       END IF
       WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') ElemID, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(i) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(i) % NodeIndexes(j)
          END IF
          WRITE(1,'(i0,x)', ADVANCE='NO') m
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    !Write boundary file
    WarnNoTarget = .FALSE.
    OPEN( 1,FILE=TRIM(Path) // boundFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex

       IF(Parallel) THEN
          IF(parent1 /= 0) parent1 = NewMesh % Elements(parent1) % GElementIndex
          IF(parent2 /= 0) parent2 = NewMesh % Elements(parent2) % GElementIndex
       END IF

       Constraint = NewMesh % Elements(k) % BoundaryInfo % Constraint
       BList => ListGetIntegerArray( Model % BCs(Constraint) % Values, &
            'Target Boundaries', Found )
       IF(Found) THEN
          IF(SIZE(BList) > 1) THEN
             CALL WARN("WriteMeshToDisk2",&
                  "A BC has more than one Target Boundary, SaveMesh output will not match input!")
          END IF
          meshBC = BList(1)
       ELSE
          WarnNoTarget = .TRUE.
          meshBC = Constraint
       END IF

       !This meshBC stuff will *only* work if each BC has only 1 target boundary
       WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            meshBC, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(k) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(k) % NodeIndexes(j)
          END IF
          WRITE(1,'(x,i0)', ADVANCE='NO') m
       END DO
       WRITE(1,*) !blank write statement to create new line without extra space.
    END DO
    CLOSE(1)

    IF(WarnNoTarget) THEN
       CALL WARN("WriteMeshToDisk2","Couldn't find a Target Boundary, assuming mapping to self")
    END IF

    IF(.NOT. Parallel) RETURN

    !Write .shared file
    !Need to create part.n.shared from Mesh % ParallelInfo %
    !NeighbourList % Neighbours.
    OPEN( 1,FILE=TRIM(Path) // sharedFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       nneigh = SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
            Neighbours)
       IF(nneigh < 2) CYCLE
       WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') &
            NewMesh % ParallelInfo % GlobalDOFs(i),nneigh
       DO j=1,nneigh
          WRITE(1,'(I0, x)',ADVANCE='NO') NewMesh % ParallelInfo %&
               NeighbourList(i) % Neighbours(j) + 1
       END DO
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)


!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDiskPartitioned(Model, Mesh, Path, &
      ElementPart, NeighbourList )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: Path
    INTEGER, POINTER :: ElementPart(:)
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: NoBoundaryElements, NoBulkElements, NoNodes, NoPartitions, Partition
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeCounts(827),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    LOGICAL :: Found, Hit
    CHARACTER(LEN=MAX_NAME_LEN) :: DirectoryName, PrefixName
!------------------------------------------------------------------------------

    NoPartitions = MAXVAL( ElementPart ) 
    NumElmCodes = 0
    NumElements = Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements
        
    WRITE(DirectoryName, '(A,A,I0)') TRIM(PATH),'/partitioning.',NoPartitions
    CALL MakeDirectory( TRIM(DirectoryName) // CHAR(0) )
    CALL Info('WriteMeshToDiskPartitioned','Writing parallel mesh to disk: '//TRIM(DirectoryName))
   

    DO Partition = 1, NoPartitions 
      
      CALL Info('WriteMeshToDiskPartitioned','Writing piece to file: '//TRIM(I2S(Partition)),Level=12)
      
      WRITE( PrefixName,'(A,A,I0)') TRIM(DirectoryName),'/part.',Partition  

      CALL Info('WriteMeshToDiskPartitioned','Write nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.nodes', STATUS='UNKNOWN' )
      NoNodes = 0
      DO i=1,Mesh % NumberOfNodes
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          WRITE(1,'(I0,x,I0,x,3ES17.10)') i,-1, &
              Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
          NoNodes = NoNodes + 1
        END IF
      END DO
      CLOSE(1)
      

      CALL Info('WriteMeshToDiskPartitioned','Write shared nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.shared', STATUS='UNKNOWN' )
      NoShared = 0
      DO i=1,Mesh % NumberOfNodes
        nneigh = SIZE( NeighbourList(i) % Neighbours )
        IF( nneigh <= 1 ) CYCLE
        
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          NoShared = NoShared + 1
          WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') i,nneigh
          DO j=1,nneigh
            WRITE(1,'(I0, x)',ADVANCE='NO') NeighbourList(i) % Neighbours(j) 
          END DO
          WRITE( 1,* ) ''
        END IF
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write elements file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.elements', STATUS='UNKNOWN' )
      NoBulkElements = 0
      ElmCodeCounts = 0      
      DO i=1,Mesh % NumberOfBulkElements
        IF( ElementPart(i) /= Partition ) CYCLE

        Element => Mesh % Elements(i)
        WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') i, &
            Element % BodyId, Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) ''
        
        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBulkElements = NoBulkElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write boundary file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.boundary', STATUS='UNKNOWN' )
      NoBoundaryElements = 0
      DO i=Mesh % NumberOfBulkElements +1 ,&
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(i)
       
        parent1 = 0
        parent2 = 0
        Constraint = 0
        
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) &
              parent1 = Element % BoundaryInfo % Left % ElementIndex
          IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) &
              parent2 = Element % BoundaryInfo % Right % ElementIndex        
          Constraint = Element % BoundaryInfo % Constraint
        END IF

        Hit = .FALSE.
        IF( parent1 > 0 ) THEN
          IF( ElementPart( parent1 ) == Partition ) Hit = .TRUE.
        END IF
        IF( parent2 > 0 ) THEN
          IF( ElementPart( parent2 ) == Partition ) Hit = .TRUE.
        END IF

        IF( .NOT. Hit ) CYCLE

        WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            Constraint, Parent1, Parent2,&
            Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(x,i0)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) 

        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBoundaryElements = NoBoundaryElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write header file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.header',STATUS='UNKNOWN' )
      NumElmCodes = COUNT( ElmCodeCounts > 0 ) 
      WRITE( 1,'(i0,x,i0,x,i0)' ) NoNodes, &
          NoBulkElements, NoBoundaryElements      
      WRITE( 1,'(i0)' ) NumElmCodes
      DO i=SIZE(ElmCodeCounts),1,-1
        IF( ElmCodeCounts(i) == 0 ) CYCLE
        WRITE( 1,'(i0,x,i0,x)' ) i,ElmCodeCounts(i)
      END DO
      WRITE( 1,'(i0,x,i0)') NoShared, 0
      CLOSE(1)
      
      CALL Info('WriteMeshToDiskPartitioned','Done writing partition',Level=12)
    END DO

    CALL Info('WriteMeshToDiskPartitioned','Done writing parallel mesh',Level=8)

!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDiskPartitioned
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Generate element edge (faces in 3D) tables for given mesh.
!> Currently only for triangles and tetras. If mesh already
!> has edges do nothing.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges( Mesh, FindEdges)
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: FindEdges

     LOGICAL :: FindEdges3D
     INTEGER :: MeshDim, SpaceDim, MaxElemDim 

     IF(PRESENT(FindEdges)) THEN
       FindEdges3D = FindEdges
     ELSE
       FindEdges3D = .TRUE.
     END IF

!------------------------------------------------------------------------------

     SpaceDim = CoordinateSystemDimension()
     MeshDim = Mesh % MeshDim

     IF( MeshDim == 0 ) THEN
       CALL Fatal('FindMeshEdges','Mesh dimension is zero!')
     END IF
     IF( SpaceDim > MeshDim ) THEN
       CALL Warn('FindMeshEdges','Mesh dimension and space dimension differ: '&
           // TRIM(I2S(MeshDim))//' vs. '//TRIM(I2S(SpaceDim)))
     END IF

     MaxElemDim = EnsureElemDim( MeshDim ) 
     IF( MaxElemDim < MeshDim ) THEN
       CALL Warn('FindMeshEdges','Element dimension smaller than mesh dimension: '//&
           TRIM(I2S(MaxElemDim))//' vs '//TRIM(I2S(MeshDim)))
     END IF


     SELECT CASE( MaxElemDim )

     CASE(2)
       IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
         CALL Info('FindMeshEdges','Determining edges in 2D mesh',Level=8)
         CALL FindMeshEdges2D( Mesh )
       END IF

     CASE(3)
       IF ( .NOT.ASSOCIATED( Mesh % Faces) ) THEN
         CALL Info('FindMeshEdges','Determining faces in 3D mesh',Level=8)
         CALL FindMeshFaces3D( Mesh )
       END IF
       IF(FindEdges3D) THEN
         IF ( .NOT.ASSOCIATED( Mesh % Edges) ) THEN
           CALL Info('FindMeshEdges','Determining edges in 3D mesh',Level=8)
           CALL FindMeshEdges3D( Mesh )
         END IF
       END IF
     END SELECT

     CALL AssignConstraints()

CONTAINS

  ! Check that the element dimension really follows the mesh dimension
  ! The default is the MeshDim so we return immediately after that is 
  ! confirmed. 
  !--------------------------------------------------------------------
    FUNCTION EnsureElemDim(MeshDim) RESULT (MaxElemDim)

      INTEGER :: MeshDim, MaxElemDim 
      INTEGER :: i,ElemDim, ElemCode

      MaxElemDim = 0

      DO i=1,Mesh % NumberOfBulkElements
        ElemCode = Mesh % Elements(i) % Type % ElementCode
        IF( ElemCode > 500 ) THEN
          ElemDim = 3 
        ELSE IF( ElemCode > 300 ) THEN
          ElemDim = 2
        ELSE IF( ElemCode > 200 ) THEN
          ElemDim = 1
        END IF
        MaxElemDim = MAX( MaxElemDim, ElemDim ) 
        IF( MaxElemDim == MeshDim ) EXIT
      END DO
          
    END FUNCTION EnsureElemDim


    SUBROUTINE AssignConstraints()

      INTEGER, POINTER :: FaceInd(:)
      INTEGER :: i,j,k,l,n,nd,nfound
      TYPE(Element_t), POINTER :: Element, Boundary, Face, Faces(:)

      DO i=1,Mesh % NumberOfBoundaryElements
        Boundary => Mesh % Elements(Mesh % NumberOfBulkElements+i)

        Element  => Boundary % BoundaryInfo % Left
        IF (.NOT.ASSOCIATED(Element) ) &
          Element  => Boundary % BoundaryInfo % Right
        IF (.NOT.ASSOCIATED(Element) ) CYCLE

        SELECT CASE(Boundary % TYPE % DIMENSION)
        CASE(1)
          nd = Element % TYPE % NumberOfEdges
          Faces   => Mesh % Edges
          FaceInd => Element % EdgeIndexes
        CASE(2)
          nd = Element % TYPE % NumberOfFaces
          Faces   => Mesh % Faces
          FaceInd => Element % FaceIndexes
        CASE DEFAULT
          Faces => NULL()
          FaceInd => NULL()
        END SELECT

        IF ( .NOT. ASSOCIATED(Faces) .OR. .NOT. ASSOCIATED(FaceInd) ) CYCLE

        DO j=1,nd
          Face => Faces(FaceInd(j))
          IF ( .NOT.ASSOCIATED(Face % TYPE,Boundary % TYPE) ) CYCLE

          n = Boundary % TYPE % NumberOfNodes
          nfound = 0
          DO k=1,n
            DO l=1,n
              IF ( Boundary % NodeIndexes(k)==Face % NodeIndexes(l) ) &
                nfound = nfound+1
            END DO
          END DO
          IF ( nfound==n ) THEN
            Face % BoundaryInfo % Constraint = Boundary % BoundaryInfo % Constraint; EXIT
          END IF
        END DO
      END DO
    END SUBROUTINE AssignConstraints
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Find 2D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges2D( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
     
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    TYPE(Element_t), POINTER :: Element, Edges(:)

    LOGICAL :: Found
    INTEGER :: i,j,k,n,NofEdges,Edge,Swap,Node1,Node2,istat,Degree
!------------------------------------------------------------------------------
!
!   Initialize:
!   -----------
    CALL AllocateVector( Mesh % Edges, 4*Mesh % NumberOfBulkElements )
    Edges => Mesh % Edges

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector( Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       SELECT CASE( Element % TYPE % ElementCode / 100 )
         CASE(3)
            n = 3
         CASE(4)
            n = 4
       END SELECT

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n
!         We use MIN(Node1,Node2) as the hash table key:
!         ----------------------------------------------
          Node1 = Element % NodeIndexes(k)
          IF ( k<n ) THEN
             Node2 = Element % NodeIndexes(k+1)
          ELSE
             Node2 = Element % NodeIndexes(1)
          END IF

          IF ( Node2 < Node1 ) THEN
             Swap  = Node1
             Node1 = Node2
             Node2 = Swap
          END IF

!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.         
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Exisiting edge, update structures:
!         ----------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge
             Edges(Edge) % BoundaryInfo % Right => Element
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges

             Degree = Element % TYPE % BasisFunctionDegree

             Edges(Edge) % ElementIndex = Edge
             CALL AllocateVector( Edges(Edge) % NodeIndexes, Degree+1)
             ALLOCATE( Edges(Edge) % BoundaryInfo )
             Edges(Edge) % TYPE => GetElementType( 201+Degree, .FALSE. )

             Edges(Edge) % NodeIndexes(1) = Element % NodeIndexes(k)
             IF ( k < n ) THEN
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(k+1)
             ELSE
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(1)
             END IF

             DO j=2,Degree
                Edges(Edge) % NodeIndexes(j+1) = Element % NodeIndexes(k+n+j-2)
             END DO
             
             ! Create P element definitions if needed
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
               CALL AllocatePDefinitions(Edges(Edge))
               Edges(Edge) % PDefs % P = 0
             ELSE
               NULLIFY( Edges(Edge) % PDefs )
             END IF

             Edges(Edge) % NDofs = 0
             IF (Element % NDOFs /= 0 ) &
                Edges(Edge) % NDOFs  = Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             NULLIFY( Edges(Edge) % EdgeIndexes )
             NULLIFY( Edges(Edge) % FaceIndexes )

             Element % EdgeIndexes(k) = Edge

             Edges(Edge) % BoundaryInfo % Left => Element
             NULLIFY( Edges(Edge) % BoundaryInfo % Right )
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Edge = Edge
             HashPtr % Node = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfEdges = NofEdges

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshFaces3D( Mesh )
    USE PElementMaps, ONLY : GetElementFaceMap
    USE PElementBase, ONLY : isPTetra

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Node2,Face
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found
    INTEGER :: n1,n2,n3,n4
    INTEGER :: i,j,k,n,NofFaces,Face,Swap,Node1,Node2,Node3,istat,Degree
     
    TYPE(Element_t), POINTER :: Element, Faces(:)

    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,6), BrickFaceMap(6,9), &
         WedgeFaceMap(5,4), PyramidFaceMap(5,8)
    
    INTEGER :: nf(4)
!------------------------------------------------------------------------------
    

    
    TetraFaceMap(1,:) = (/ 1, 2, 3, 5, 6, 7 /)
    TetraFaceMap(2,:) = (/ 1, 2, 4, 5, 9, 8 /)
    TetraFaceMap(3,:) = (/ 2, 3, 4, 6,10, 9 /)
    TetraFaceMap(4,:) = (/ 3, 1, 4, 7, 8,10 /)

    WedgeFaceMap(1,:) = (/ 1, 2, 3,-1 /)
    WedgeFaceMap(2,:) = (/ 4, 5, 6,-1 /)
    WedgeFaceMap(3,:) = (/ 1, 2, 5, 4 /)
    WedgeFaceMap(4,:) = (/ 3, 2, 5, 6 /)
    WedgeFaceMap(5,:) = (/ 3, 1, 4, 6 /)

    PyramidFaceMap(1,:) = (/ 1, 2, 3, 4,  6,  7,  8,  9 /)
    PyramidFaceMap(2,:) = (/ 1, 2, 5, 6, 11, 10, -1, -1 /)
    PyramidFaceMap(3,:) = (/ 2, 3, 5, 7, 12, 11, -1, -1 /)
    PyramidFaceMap(4,:) = (/ 3, 4, 5, 8, 13, 12, -1, -1 /)
    PyramidFaceMap(5,:) = (/ 4, 1, 5, 9, 10, 13, -1, -1 /)

    BrickFaceMap(1,:) = (/ 1, 2, 3, 4,  9, 10, 11, 12, 25 /)
    BrickFaceMap(2,:) = (/ 5, 6, 7, 8, 17, 18, 19, 20, 26 /)
    BrickFaceMap(3,:) = (/ 1, 2, 6, 5,  9, 14, 17, 13, 21 /)
    BrickFaceMap(4,:) = (/ 2, 3, 7, 6, 10, 15, 17, 14, 22 /)
    BrickFaceMap(5,:) = (/ 3, 4, 8, 7, 11, 16, 19, 15, 23 /)
    BrickFaceMap(6,:) = (/ 4, 1, 5, 8, 12, 13, 20, 16, 24 /)

!
!   Initialize:
!   -----------
    CALL AllocateVector( Mesh % Faces, 6*Mesh % NumberOfBulkElements, 'FindMeshFaces3D' )
    Faces => Mesh % Faces

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % FaceIndexes ) ) &
          CALL AllocateVector(Element % FaceIndexes, Element % TYPE % NumberOfFaces )
       Element % FaceIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofFaces = 0
    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementFaceMap(Element, FaceMap)
          n = Element % TYPE % NumberOfFaces
       ELSE
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(5)
             n = 4
             FaceMap => TetraFaceMap
          CASE(6)
             n = 5
             FaceMap => PyramidFaceMap
          CASE(7)
             n = 5 
             FaceMap => WedgeFaceMap
          CASE(8)
             n = 6
             FaceMap => BrickFaceMap
          CASE DEFAULT
             CYCLE
             ! WRITE(Message,*) 'Element type',Element % Type % ElementCode,'not implemented.' 
             ! CALL Fatal('FindMeshFaces',Message)
          END SELECT
       END IF
 
!      Loop over every face of every element:
!      --------------------------------------
       DO k=1,n
          
          
!         We use MIN(Node1,Node2,Node3) as the hash table key:
!         ---------------------------------------------------
          SELECT CASE( Element % TYPE % ElementCode / 100 )
             CASE(5)
!
!               Tetras:
!               =======
                nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                CALL sort( 3, nf )

             CASE(6)
!
!               Pyramids:
!               =========
                IF ( k == 1 ) THEN
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                ELSE
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                END IF

             CASE(7)
!
!               Wedges:
!               =======
                IF ( k <= 2 ) THEN
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                ELSE
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                END IF
                
             CASE(8)
!
!               Bricks:
!               =======
                nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                CALL sort( 4, nf )

             CASE DEFAULT
                WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
                CALL Fatal('FindMeshFaces',Message)
          END SELECT

          Node1 = nf(1)
          Node2 = nf(2)
          Node3 = nf(3)
          
!         Look the face from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 .AND. HashPtr % Node2 == Node3) THEN
                Found = .TRUE.
                Face = HashPtr % Face
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Exisiting face, update structures:
!         ----------------------------------
          IF ( Found ) THEN       
             Element % FaceIndexes(k) = Face
             Faces(Face) % BoundaryInfo % Right => Element
          ELSE

!            Face not yet there, create:
!            ---------------------------
             NofFaces = NofFaces + 1
             Face = NofFaces
             Faces(Face) % ElementIndex = Face

             Degree = Element % TYPE % BasisFunctionDegree
             SELECT CASE( Element % TYPE % ElementCode / 100 )
                CASE(5)
!
!               for tetras:
!               -----------
                SELECT CASE( Degree ) 
                   CASE(1)
                   n1 = 3
                   CASE(2)
                   n1 = 6
                   CASE(3)
                   n1 = 10
                END SELECT
                n1 = 3
                
                Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )

                CASE(6)

!               Pyramids ( only 605 supported )
!               -------------------------------
                IF ( k == 1 ) THEN
                   n1 = 4
                   Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
                ELSE
                   n1 = 3
                   Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
                END IF
                
                CASE(7)

!               for wedges, only 706 supported:
!               -------------------------------
                IF ( k <= 2 ) THEN
                   n1 = 3
                   Faces(Face) % TYPE => GetElementType( 303, .FALSE. )
                ELSE
                   n1 = 4
                   Faces(Face) % TYPE => GetElementType( 404, .FALSE. )
                END IF

            
                CASE(8)
!
!               for bricks:
!               -----------
                SELECT CASE( Element % TYPE % NumberOfNodes ) 
                   CASE(8)
                   n1 = 4
                   CASE(20)
                   n1 = 8
                   CASE(27)
                   n1 = 9
                END SELECT

                Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE.)

                CASE DEFAULT
                   WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
                   CALL Fatal('FindMeshFaces',Message)

             END SELECT

             ! Allocate p structures for p elements
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Faces(Face))
                Faces(Face) % PDefs % P = 0
             ELSE
               NULLIFY( Faces(Face) % PDefs )
             END IF
             
             Faces(Face) % NDOFs  = 0
             IF (Element % NDOFs /= 0 ) &
                Faces(Face) % NDOFs  = Faces(Face) % TYPE % NumberOfNodes
             Faces(Face) % BDOFs  = 0
             Faces(Face) % DGDOFs = 0
             Faces(Face) % EdgeIndexes => NULL()
             Faces(Face) % FaceIndexes => NULL()

             CALL AllocateVector( Faces(Face) % NodeIndexes,n1 )
             DO n2=1,n1
                Faces(Face) % NodeIndexes(n2) = &
                         Element % NodeIndexes(FaceMap(k,n2)) 
             END DO

             Element % FaceIndexes(k) = Face

             ALLOCATE( Faces(Face) % BoundaryInfo )
             Faces(Face) % BoundaryInfo % Left => Element
             NULLIFY( Faces(Face) % BoundaryInfo % Right )
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Face = Face
             HashPtr % Node1 = Node2
             HashPtr % Node2 = Node3
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfFaces = NofFaces

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshFaces3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges3D( Mesh )
    USE PElementMaps, ONLY : GetElementEdgeMap, GetElementFaceEdgeMap
    USE PElementBase, ONLY : isPPyramid

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found
    INTEGER :: n1,n2
    INTEGER :: i,j,k,n,NofEdges,Edge,Node1,Node2,istat,Degree,ii,jj
     
    TYPE(Element_t), POINTER :: Element, Edges(:), Face

    INTEGER, POINTER :: EdgeMap(:,:), FaceEdgeMap(:,:)
    INTEGER, TARGET  :: TetraEdgeMap(6,3), BrickEdgeMap(12,3), TetraFaceMap(4,6), &
      WedgeEdgeMap(9,3), PyramidEdgeMap(8,3), TetraFaceEdgeMap(4,3), &
      BrickFaceEdgeMap(8,4), WedgeFaceEdgeMap(6,4), PyramidFaceEdgeMap(5,4)
!------------------------------------------------------------------------------
    TetraFaceMap(1,:) = (/ 1, 2, 3, 5, 6, 7 /)
    TetraFaceMap(2,:) = (/ 1, 2, 4, 5, 9, 8 /)
    TetraFaceMap(3,:) = (/ 2, 3, 4, 6,10, 9 /)
    TetraFaceMap(4,:) = (/ 3, 1, 4, 7, 8,10 /)

    TetraFaceEdgeMap(1,:) = (/ 1,2,3 /)
    TetraFaceEdgeMap(2,:) = (/ 1,5,4 /)
    TetraFaceEdgeMap(3,:) = (/ 2,6,5 /)
    TetraFaceEdgeMap(4,:) = (/ 3,4,6 /)

    TetraEdgeMap(1,:) = (/ 1,2,5 /)
    TetraEdgeMap(2,:) = (/ 2,3,6 /)
    TetraEdgeMap(3,:) = (/ 3,1,7 /)
    TetraEdgeMap(4,:) = (/ 1,4,8 /)
    TetraEdgeMap(5,:) = (/ 2,4,9 /)
    TetraEdgeMap(6,:) = (/ 3,4,10 /)

    PyramidEdgeMap(1,:) = (/ 1,2,1 /)
    PyramidEdgeMap(2,:) = (/ 2,3,1 /)
    PyramidEdgeMap(3,:) = (/ 3,4,1 /)
    PyramidEdgeMap(4,:) = (/ 4,1,1 /)
    PyramidEdgeMap(5,:) = (/ 1,5,1 /)
    PyramidEdgeMap(6,:) = (/ 2,5,1 /)
    PyramidEdgeMap(7,:) = (/ 3,5,1 /)
    PyramidEdgeMap(8,:) = (/ 4,5,1 /)

    PyramidFaceEdgeMap(1,:) = (/ 1,2,3,4 /)
    PyramidFaceEdgeMap(2,:) = (/ 1,6,5,0 /)
    PyramidFaceEdgeMap(3,:) = (/ 2,7,6,0 /)
    PyramidFaceEdgeMap(4,:) = (/ 3,8,7,0 /)
    PyramidFaceEdgeMap(5,:) = (/ 4,5,8,0 /)

    WedgeEdgeMap(1,:) = (/ 1, 2,1 /)
    WedgeEdgeMap(2,:) = (/ 2, 3,1 /)
    WedgeEdgeMap(3,:) = (/ 1, 3,1 /)
    WedgeEdgeMap(4,:) = (/ 4, 5,1 /)
    WedgeEdgeMap(5,:) = (/ 5, 6,1 /)
    WedgeEdgeMap(6,:) = (/ 6, 4,1 /)
    WedgeEdgeMap(7,:) = (/ 1, 4,1 /)
    WedgeEdgeMap(8,:) = (/ 2, 5,1 /)
    WedgeEdgeMap(9,:) = (/ 3, 6,1 /)

    WedgeFaceEdgeMap(1,:) = (/ 1,2,3,0 /)
    WedgeFaceEdgeMap(2,:) = (/ 4,5,6,0 /)
    WedgeFaceEdgeMap(3,:) = (/ 1,8,4,7 /)
    WedgeFaceEdgeMap(4,:) = (/ 2,9,5,8 /)
    WedgeFaceEdgeMap(5,:) = (/ 3,7,6,9 /)

    BrickEdgeMap(1,:) = (/ 1, 2,  9 /)
    BrickEdgeMap(2,:) = (/ 2, 3,  10 /)
    BrickEdgeMap(3,:) = (/ 4, 3,  11 /)
    BrickEdgeMap(4,:) = (/ 1, 4,  12 /)
    BrickEdgeMap(5,:) = (/ 5, 6,  13 /)
    BrickEdgeMap(6,:) = (/ 6, 7,  14 /)
    BrickEdgeMap(7,:) = (/ 8, 7,  15 /)
    BrickEdgeMap(8,:) = (/ 5, 8,  16 /)
    BrickEdgeMap(9,:) = (/ 1, 5,  17 /)
    BrickEdgeMap(10,:) = (/ 2, 6, 18 /)
    BrickEdgeMap(11,:) = (/ 3, 7, 19 /)
    BrickEdgeMap(12,:) = (/ 4, 8, 20 /)

    BrickFaceEdgeMap(1,:) = (/ 1,2,3,4   /)
    BrickFaceEdgeMap(2,:) = (/ 5,6,7,8   /)    
    BrickFaceEdgeMap(3,:) = (/ 1,10,5,9  /)
    BrickFaceEdgeMap(4,:) = (/ 2,11,6,10 /)
    BrickFaceEdgeMap(5,:) = (/ 3,12,7,11 /)
    BrickFaceEdgeMap(6,:) = (/ 4,9,8,12  /)

!
!   Initialize:
!   -----------
    CALL AllocateVector( Mesh % Edges, 12*Mesh % NumberOfBulkElements )
    Edges => Mesh % Edges

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector(Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementEdgeMap( Element, EdgeMap )
          CALL GetElementFaceEdgeMap( Element, FaceEdgeMap ) 
          n = Element % TYPE % NumberOfEdges
       ELSE 
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(5)
             n = 6
             EdgeMap => TetraEdgeMap
             FaceEdgeMap => TetraFaceEdgeMap
          CASE(6)
             n = 8
             EdgeMap => PyramidEdgeMap
             FaceEdgeMap => PyramidFaceEdgeMap
          CASE(7)
             n = 9
             EdgeMap => WedgeEdgeMap
             FaceEdgeMap => WedgeFaceEdgeMap
          CASE(8)
             n = 12
             EdgeMap => BrickEdgeMap
             FaceEdgeMap => BrickFaceEdgeMap
          CASE DEFAULT
             CYCLE
             WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
             CALL Fatal('FindMeshEdges',Message)
          END SELECT
       END IF

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n

!         Use MIN(Node1,Node2) as key to hash table:
!         ------------------------------------------
          n1 = Element % NodeIndexes(EdgeMap(k,1))
          n2 = Element % NodeIndexes(EdgeMap(k,2))
          IF ( n1 < n2 ) THEN
             Node1 = n1
             Node2 = n2
          ELSE
             Node1 = n2
             Node2 = n1
          END IF
!
!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO
!
!         Existing edge, update structures:
!         ---------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge

             ! Mark edge as an edge of pydamid square face 
             IF (isPPyramid(Element) .AND. k < 5) THEN
                Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
             END IF

             IF ( ASSOCIATED(Mesh % Faces) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces(Element % FaceIndexes(ii))
                 IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                   ALLOCATE(Face % EdgeIndexes(Face % TYPE % NumberOfEdges))
                   Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF (FaceEdgeMap(ii,jj) == k) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                          Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                          Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                       EXIT
                    END IF
                 END DO
               END DO
             END IF
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges
             Edges(Edge) % ElementIndex = Edge
             Degree = Element % TYPE % BasisFunctionDegree

!            Edge is always a line segment with deg+1 nodes:
!            -----------------------------------------------
             Edges(Edge) % TYPE => GetElementType( 201 + degree, .FALSE.)

             Edges(Edge) % NDOFs  = 0
             IF (Element % NDOFs /= 0 ) &
                Edges(Edge) % NDOFs  = Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             Edges(Edge) % EdgeIndexes => NULL()
             Edges(Edge) % FaceIndexes => NULL()

             CALL AllocateVector( Edges(Edge) % NodeIndexes, degree + 1 )
             DO n2=1,degree+1
               Edges(Edge) % NodeIndexes(n2) = &
                    Element % NodeIndexes(EdgeMap(k,n2))
             END DO

             Element % EdgeIndexes(k) = Edge
             ALLOCATE( Edges(Edge) % BoundaryInfo )
             Edges(Edge) % BoundaryInfo % Left  => NULL()
             Edges(Edge) % BoundaryInfo % Right => NULL()

             ! Allocate P element definitions 
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Edges(Edge))
             
                Edges(Edge) % PDefs % P = 0
                Edges(Edge) % PDefs % pyramidQuadEdge = .FALSE.
                ! Here mark edge as edge of pyramid if needed (or set as not)
                IF (isPPyramid(Element) .AND. k < 5) THEN
                   Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
                END IF
             ELSE
                NULLIFY( Edges(Edge) % PDefs )
             END IF

             IF ( ASSOCIATED(Mesh % Faces) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces( Element % FaceIndexes(ii) )
                 IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                    ALLOCATE( Face % EdgeIndexes( Face % TYPE % NumberOfEdges ) )
                    Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF ( FaceEdgeMap(ii,jj) == k ) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF (.NOT.ASSOCIATED( Edges(Edge) % BoundaryInfo % Left)) THEN
                          Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                          Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                    END IF
                 END DO
               END DO
             END IF

!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Edge = Edge
             HashPtr % Node1 = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfEdges = NofEdges

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    IF (ASSOCIATED(Mesh % Faces)) CALL FixFaceEdges()

CONTAINS 

    SUBROUTINE FixFaceEdges()

      INTEGER :: i,j,k,n,swap,edgeind(4),i1(2),i2(2)

      DO i=1,Mesh % NumberOfFaces
        Face => Mesh % Faces(i)
        n = Face % TYPE % NumberOfEdges
        Edgeind(1:n) = Face % EdgeIndexes(1:n)
        DO j=1,n
          i1 = Mesh % Edges(Edgeind(j)) % NodeIndexes(1:2)
          IF ( i1(1)>i1(2) ) THEN
            swap=i1(1)
            i1(1)=i1(2)
            i1(2)=swap
          END IF
          DO k=1,n
            i2(1) = k
            i2(2) = k+1
            IF ( i2(2)>n ) i2(2)=1
            i2 = Face % NodeIndexes(i2)
            IF ( i2(1)>i2(2) ) THEN
              swap=i2(1)
              i2(1)=i2(2)
              i2(2)=swap
            END IF
            IF ( ALL(i1 == i2) ) THEN
              Face % EdgeIndexes(k) = edgeind(j)
              EXIT
            END IF
          END DO
        END DO
      END DO
    END SUBROUTINE FixFaceEdges
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finds neigbours of the nodes in given direction.
!> The algorithm finds the neighbour that within 45 degrees of the 
!> given direction has the smallest distance.
!------------------------------------------------------------------------------
  SUBROUTINE FindNeighbourNodes( Mesh,Direction,Neighbours,EndNeighbours)
!------------------------------------------------------------------------------

  TYPE(Mesh_t) , POINTER :: Mesh 
  REAL(KIND=dp) :: Direction(:)
  INTEGER :: Neighbours(:)
  INTEGER, OPTIONAL :: EndNeighbours(:)

  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  REAL(KIND=dp), POINTER :: Distances(:)
  REAL(KIND=dp) :: rn(3), rs(3), ss, sn
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: i,j,k,n,t,DIM,istat

  IF(SIZE(Neighbours) < Mesh % NumberOfNodes) THEN
    CALL Warn('FindNeigbourNodes','SIZE of Neigbours should equal Number of Nodes!')
    RETURN
  END IF


  IF(PRESENT(EndNeighbours)) THEN
    IF(SIZE(EndNeighbours) < Mesh % NumberOfNodes) THEN
      CALL Warn('FindNeigbourNodes','SIZE of EndNeigbours should equal Number of Nodes!')
      RETURN
    END IF
  END IF


  DIM = CoordinateSystemDimension()
  N = Mesh % MaxElementNodes

  CALL AllocateVector( ElementNodes % x, n )
  CALL AllocateVector( ElementNodes % y, n )
  CALL AllocateVector( ElementNodes % z, n )
  CALL AllocateVector( Distances, Mesh % NumberOfNodes )

  Neighbours = 0
  Distances = HUGE(Distances)
 
  rn(1:DIM) = Direction(1:DIM)
  ss = SQRT(SUM(rn(1:DIM)**2))
  rn = rn / ss

  DO t=1,Mesh % NumberOfBulkElements

    CurrentElement => Mesh % Elements(t)
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
  
    ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
    ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
    IF(DIM == 3) THEN
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
    END IF


    DO i=1,n
      DO j=i+1,n
        rs(1) = ElementNodes % x(j) - ElementNodes % x(i)
        rs(2) = ElementNodes % y(j) - ElementNodes % y(i)
        IF (DIM == 3) THEN
          rs(3) = ElementNodes % z(j) - ElementNodes % z(i)
        END IF
        
        ss = SQRT(SUM(rs(1:DIM)**2))
        sn = SUM(rs(1:DIM)*rn(1:DIM))

        IF(ss < SQRT(2.0) * ABS(sn)) THEN
          IF(sn > 0) THEN
            IF(ss < Distances(NodeIndexes(i))) THEN
              Distances(NodeIndexes(i)) = ss
              Neighbours(NodeIndexes(i)) = NodeIndexes(j)
            END IF
          ELSE
            IF(ss < Distances(NodeIndexes(j))) THEN
              Distances(NodeIndexes(j)) = ss
              Neighbours(NodeIndexes(j)) = NodeIndexes(i)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO

  ! This loop finds the final neighbour in the end of the chain 
  IF(PRESENT(EndNeighbours)) THEN
    EndNeighbours = Neighbours

    DO t=1,Mesh%NumberOfNodes
      j = Neighbours(t)
      DO WHILE(j /= 0)
        EndNeighbours(t) = j
        j = Neighbours(j)
      END DO
    END DO
  END IF
  DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z, Distances)
!------------------------------------------------------------------------------
END SUBROUTINE FindNeighbourNodes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE UpdateSolverMesh( Solver, Mesh )
!------------------------------------------------------------------------------
     TYPE( Mesh_t ), POINTER :: Mesh
     TYPE( Solver_t ), TARGET :: Solver
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,n1,n2,DOFs
     LOGICAL :: Found, OptimizeBandwidth
     TYPE(Matrix_t), POINTER   :: Matrix
     REAL(KIND=dp), POINTER :: Work(:)
     INTEGER, POINTER :: Permutation(:)
     TYPE(Variable_t), POINTER :: TimeVar, SaveVar
!------------------------------------------------------------------------------
     SaveVar => Solver % Variable
     DOFs = SaveVar % DOFs

     Solver % Mesh => Mesh
     CALL SetCurrentMesh( CurrentModel, Mesh )
!
!    Create matrix and variable structures for
!    current equation on the new mesh:
!    -----------------------------------------
     Solver % Variable => VariableGet( Mesh % Variables, &
        Solver % Variable % Name, ThisOnly = .FALSE. )

     CALL AllocateVector( Permutation, SIZE(Solver % Variable % Perm) )

     OptimizeBandwidth = ListGetLogical( Solver % Values, 'Optimize Bandwidth', Found )
     IF ( .NOT. Found ) OptimizeBandwidth = .TRUE.

     Matrix => CreateMatrix( CurrentModel, Solver, &
        Mesh, Permutation, DOFs, MATRIX_CRS, OptimizeBandwidth, &
        ListGetString( Solver % Values, 'Equation' ) )

     Matrix % Symmetric = ListGetLogical( Solver % Values, &
             'Linear System Symmetric', Found )

     Matrix % Lumped = ListGetLogical( Solver % Values, &
             'Lumped Mass Matrix', Found )

     ALLOCATE( Work(SIZE(Solver % Variable % Values)) )
     Work = Solver % Variable % Values
     DO k=0,DOFs-1
        DO i=1,SIZE(Permutation)
           IF ( Permutation(i) > 0 ) THEN
              Solver % Variable % Values( DOFs*Permutation(i)-k ) = &
                 Work( DOFs*Solver % Variable % Perm(i)-k )
           END IF
        END DO
     END DO

     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        DO j=1,SIZE(Solver % Variable % PrevValues,2)
           Work = Solver % Variable % PrevValues(:,j)
           DO k=0,DOFs-1
              DO i=1,SIZE(Permutation)
                 IF ( Permutation(i) > 0 ) THEN
                    Solver % Variable % PrevValues( DOFs*Permutation(i) - k,j ) =  &
                        Work( DOFs * Solver % Variable % Perm(i) - k )
                  END IF
              END DO
           END DO
        END DO
     END IF
     DEALLOCATE( Work )

     Solver % Variable % Perm = Permutation
     Solver % Variable % Solver => Solver

     DEALLOCATE( Permutation )
     CALL AllocateVector( Matrix % RHS, Matrix % NumberOfRows )

     IF ( ASSOCIATED(SaveVar % EigenValues) ) THEN
        n = SIZE(SaveVar % EigenValues)

        IF ( n > 0 ) THEN
           Solver % NOFEigenValues = n
           CALL AllocateVector( Solver % Variable % EigenValues,n )
           CALL AllocateArray( Solver % Variable % EigenVectors, n, &
                    SIZE(Solver % Variable % Values) ) 

           Solver % Variable % EigenValues  = 0.0d0
           Solver % Variable % EigenVectors = 0.0d0

           CALL AllocateVector( Matrix % MassValues, SIZE(Matrix % Values) )
           Matrix % MassValues = 0.0d0
        END IF
     ELSE IF ( ASSOCIATED( Solver % Matrix ) ) THEN
        IF( ASSOCIATED( Solver % Matrix % Force) ) THEN
           n1 = Matrix % NumberOFRows
           n2 = SIZE(Solver % Matrix % Force,2)
           ALLOCATE(Matrix % Force(n1,n2))
           Matrix % Force = 0.0d0
        END IF
     END IF

     Solver % Matrix => Matrix
     Solver % Mesh % Changed = .TRUE.

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateSolverMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Split a mesh equally to smaller pieces by performing a uniform split.
!> Also known as mesh multiplication. A 2D element splits into 4 elements of
!> same form, and 3D element into 8 elements. 
!> Currently works only for linear elements.
!------------------------------------------------------------------------------
  FUNCTION SplitMeshEqual(Mesh,h) RESULT( NewMesh )
!------------------------------------------------------------------------------
    REAL(KIND=dp), OPTIONAL :: h(:)
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:),xh(:)
    INTEGER :: i, j, k, n, NewElCnt, NodeCnt, EdgeCnt, FaceCnt, Node, ParentId, Diag, NodeIt
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Eparent,Face,Faces(:)
    INTEGER, POINTER :: Child(:,:)
    INTEGER :: n1,n2,n3,EoldNodes(4),FaceNodes(4),EdgeNodes(2) ! Only linears so far
    INTEGER :: FaceNumber,Edge1,Edge2,Edge3,Edge4,Node12,Node23,Node34,Node41,Node31
    REAL(KIND=dp) :: dxyz(3,3),Dist(3),r,s,t,h1,h2
    TYPE(PElementDefs_t), POINTER :: PDefs
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER, ALLOCATABLE :: FacePerm(:), BulkPerm(:)
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN

    CALL Info( 'SplitMeshEqual', 'Mesh splitting works for first order elements 303, 404, 504, (706) and 808.', Level = 6 )

    DO i=1,Mesh % NumberOfBulkElements
      SELECT CASE(Mesh % Elements(i) % TYPE % ElementCode/100)
      CASE(6)
        CALL Fatal('SplitMeshEqual','Pyramids not supported, sorry.')
      END SELECT
    END DO

    NewMesh => AllocateMesh()

    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT.EdgesPresent) CALL FindMeshEdges( Mesh )

    CALL ResetTimer('SplitMeshEqual')

    CALL Info( 'SplitMeshEqual', '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( 'SplitMeshEqual', Message, Level=6 )
!
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = Mesh % NumberOfNodes + Mesh % NumberOfEdges
!
!   For quad faces add one node in the center:
!   ------------------------
    ALLOCATE(FacePerm(Mesh % NumberOfFaces)); FacePerm = 0
    FaceCnt = 0
    DO i = 1, Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       IF( Face % TYPE % NumberOfNodes == 4 ) THEN
         NodeCnt = NodeCnt+1
         FaceCnt = FaceCnt+1
         FacePerm(i) = NodeCnt
       END IF
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of faces : ', FaceCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   For quads and bricks, count centerpoints:
!   -----------------------------------------
    NodeIt = 0
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(4,8)
          NodeCnt = NodeCnt + 1
          NodeIt = NodeIt + 1
       END SELECT
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of bulks : ', NodeIt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   new mesh nodecoordinate arrays:
!   -------------------------------
    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:Mesh % NumberOfNodes) = u
    y(1:Mesh % NumberOfNodes) = v
    z(1:Mesh % NumberOfNodes) = w

! what is h? - pointer to nodal element size
    IF (PRESENT(h)) THEN
      ALLOCATE(xh(SIZE(x)))
      xh(1:SIZE(h)) = h
    END IF
!
!   add edge centers:
!   -----------------
    j =  Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
       j = j + 1
       Edge => Mesh % Edges(i)
       k = Edge % TYPE % NumberOfNodes
       IF (PRESENT(h)) THEN
         h1=h(Edge % NodeIndexes(1))
         h2=h(Edge % NodeIndexes(2))
         r=1._dp/(1+h1/h2)
         x(j) = r*u(Edge%NodeIndexes(1))+(1-r)*u(Edge%NodeIndexes(2))
         y(j) = r*v(Edge%NodeIndexes(1))+(1-r)*v(Edge%NodeIndexes(2))
         z(j) = r*w(Edge%NodeIndexes(1))+(1-r)*w(Edge%NodeIndexes(2))
         xh(j)=r*h1+(1-r)*h2
       ELSE
         x(j) = SUM(u(Edge % NodeIndexes))/k
         y(j) = SUM(v(Edge % NodeIndexes))/k
         z(j) = SUM(w(Edge % NodeIndexes))/k
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added edge centers to the nodes list.', Level=10 )  
!
!   add quad face centers for bricks and prisms(wedges):
!   ----------------------------
    j = Mesh % NumberOfNodes + Mesh % NumberOfEdges
    DO i=1,Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       k = Face % TYPE % NumberOfNodes
       IF( k == 4 ) THEN
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Face % EdgeIndexes(2))
            h2=xh(n+Face % EdgeIndexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Face % EdgeIndexes(3))
            h2=xh(n+Face % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Face,u(Face % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Face,v(Face % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Face,w(Face % NodeIndexes),r,s)
            xh(j) = InterpolateInElement2D(Face,h(Face % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Face % NodeIndexes))/k
            y(j) = SUM(v(Face % NodeIndexes))/k
            z(j) = SUM(w(Face % NodeIndexes))/k
          END IF
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added face centers to the nodes list.', Level=10 )
!
!   add centerpoint for quads & bricks:
!   -----------------------------------
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       k = Eold % TYPE % NumberOfNodes
       SELECT CASE( Eold % TYPE % ElementCode / 100 )

       CASE(4)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Eold % Edgeindexes(2))
            h2=xh(n+Eold % Edgeindexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Eold % EdgeIndexes(3))
            h2=xh(n+Eold % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Eold,u(Eold % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Eold,v(Eold % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Eold,w(Eold % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       CASE(8)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes+Mesh % NumberOfEdges
            h1=xh(n+Eold % FaceIndexes(4))
            h2=xh(n+Eold % FaceIndexes(6))
            r=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(5))
            h2=xh(n+Eold % FaceIndexes(3))
            s=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(2))
            h2=xh(n+Eold % FaceIndexes(1))
            t=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement3D(Eold,u(Eold % NodeIndexes),r,s,t)
            y(j) = InterpolateInElement3D(Eold,v(Eold % NodeIndexes),r,s,t)
            z(j) = InterpolateInElement3D(Eold,w(Eold % NodeIndexes),r,s,t)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       END SELECT
    END DO
!
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt
!
!   Update bulk elements:
!   =====================
!
!   First count new elements:
!   -------------------------
    NewElCnt = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode/100 )

!      Each element will be divided into 2**Dim new elements:
!      ------------------------------------------------------
       CASE(2)
          NewElCnt = NewElCnt + 2 ! lines
       CASE(3)
          NewElCnt = NewElCnt + 4 ! trias
       CASE(4)
          NewElCnt = NewElCnt + 4 ! quads
       CASE(5)
          NewElCnt = NewElCnt + 8 ! tetras
       CASE(7)
          NewElCnt = NewElCnt + 8 ! prisms (wedges)
       CASE(8)
          NewElCnt = NewElCnt + 8 ! hexas
       END SELECT
    END DO

    WRITE( Message, * ) 'Count of new elements : ', NewElCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info('SplitMeshEqual','New mesh allocated.', Level=10 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 8 )
    CALL Info('SplitMeshEqual','Array for bulk elements allocated.', Level=10 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes
    EdgeCnt = Mesh % NumberOfEdges

!
!   Index to old quad/hexa centerpoint node in the new mesh nodal arrays:
!   ---------------------------------------------------------------------
    Node = NodeCnt + EdgeCnt + FaceCnt
!
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)

       SELECT CASE( Eold % TYPE % ElementCode )
       CASE(303)
!
!         Split triangle to four triangles from
!         edge centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,3) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,4) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt

       CASE(404)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split quad to four new quads from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)


       CASE(504)
!
!         Split tetra to 8 new elements from
!         corners and edge centerpoints:
!         ----------------------------------
!
!         1st new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt

!         Then the annoying part; we still have to split the
!         remaining octahedron into four elements. This can
!         be done in three ways of which only one preserves
!         the minimum angle condition (Delaunay splitting):
!         --------------------------------------------------
          dxyz(1,1) = x(Eold % EdgeIndexes(4) + NodeCnt) &
                    - x(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(2,1) = y(Eold % EdgeIndexes(4) + NodeCnt) &
                    - y(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(3,1) = z(Eold % EdgeIndexes(4) + NodeCnt) &
                    - z(Eold % EdgeIndexes(2) + NodeCnt)

          dxyz(1,2) = x(Eold % EdgeIndexes(5) + NodeCnt) &
                    - x(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(2,2) = y(Eold % EdgeIndexes(5) + NodeCnt) &
                    - y(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(3,2) = z(Eold % EdgeIndexes(5) + NodeCnt) &
                    - z(Eold % EdgeIndexes(3) + NodeCnt)

          dxyz(1,3) = x(Eold % EdgeIndexes(6) + NodeCnt) &
                    - x(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(2,3) = y(Eold % EdgeIndexes(6) + NodeCnt) &
                    - y(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(3,3) = z(Eold % EdgeIndexes(6) + NodeCnt) &
                    - z(Eold % EdgeIndexes(1) + NodeCnt)

          Dist(1) = SQRT( dxyz(1,1)**2 + dxyz(2,1)**2 + dxyz(3,1)**2 )
          Dist(2) = SQRT( dxyz(1,2)**2 + dxyz(2,2)**2 + dxyz(3,2)**2 )
          Dist(3) = SQRT( dxyz(1,3)**2 + dxyz(2,3)**2 + dxyz(3,3)**2 )

          Diag = 1  ! The default diagonal for splitting is between edges 2-4
          IF (Dist(2) < Dist(1) .AND. Dist(2) < Dist(3)) Diag = 2 ! Edges 3-5
          IF (Dist(3) < Dist(1) .AND. Dist(3) < Dist(2)) Diag = 3 ! Edges 1-6

          SELECT CASE( Diag )
          CASE(1)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
          CASE(2)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
          CASE(3)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt

          END SELECT


       CASE(706)
!
!         Split prism to 8 new prism from edge
!         centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt 
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt 
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt 
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))

!
!         3rd new element (near node 3)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(9) + NodeCnt

!
!         4th new element (bottom center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt

!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt

!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
!
!         8th new element (top half, center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt



       CASE(808)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split brick to 8 new bricks from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 8)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(7) = Node
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(6))
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(8) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(6) = Node
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(12)+ NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(5) = Node
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(5))
!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(8) + NodeCnt
!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Node
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(2))
!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(12)+ NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(8) = Eold % NodeIndexes(8)
!
!         8th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(7) = Eold % NodeIndexes(7)
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(7) + NodeCnt

       CASE DEFAULT
          WRITE( Message,* ) 'Element type ', Eold % TYPE % ElementCode, &
              ' not supprted by the multigrid solver.'
          CALL Fatal( 'SplitMeshEqual', Message )
       END SELECT
    END DO

!
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt

!
!   Update boundary elements:
!   NOTE: Internal boundaries not taken care of...:!!!!
!   ---------------------------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements

       j = i + Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(j)
!
!      get parent of the boundary element:
!      -----------------------------------
       Eparent => Eold % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED(Eparent) ) &
          eParent => Eold % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED( Eparent ) ) CYCLE

       ParentId = Eparent % ElementIndex

       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(2)
!
!         Line segments:
!         ==============
!
!         which edge of the parent element are we ?
!         -----------------------------------------
          DO Edge1=1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             IF ( Eold % NodeIndexes(1) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(2) .OR.  &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(1) == Edge % NodeIndexes(2) ) EXIT
          END DO
!
!         index of the old edge centerpoint in the
!         new mesh nodal arrays:
!         ----------------------------------------
          Node = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr

       CASE(3)
!
!         Trias:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:3) = Eold % NodeIndexes(1:3)
          CALL sort( 3, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:3) = Face % NodeIndexes(1:3)
             CALL sort( 3, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) ) EXIT

          END DO
!
!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node31 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node31
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr

       CASE(4)
!
!         Quads:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:4) = Eold % NodeIndexes(1:4)
          CALL sort( 4, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:4) = Face % NodeIndexes(1:4)
             CALL sort( 4, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) .AND. &
                  EoldNodes(4) == FaceNodes(4) ) EXIT

          END DO

!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Fourth edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          DO Edge4 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge4) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node = FacePerm(Eparent % FaceIndexes(FaceNumber)) ! faces mid-point
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node34 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
          Node41 = Eparent % EdgeIndexes(Edge4) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Node41
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 )  CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          Enew % NodeIndexes(4) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node41
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Node34
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Node34
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
       END SELECT
    END DO

!
!   Update new mesh boundary element counts:
!   ----------------------------------------
    NewMesh % NumberOfBoundaryElements = NewElCnt - &
            NewMesh % NumberOfBulkElements
    NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
    NewMesh % MaxElementNodes = Mesh % MaxElementNodes

    j = 0
    DO i=1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
      Enew => NewMesh % Elements(i)

      IF ( Enew % DGDOFs>0 ) THEN
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF

      IF (i<=NewMesh % NumberOfBulkElements) THEN
         PDefs => Enew % PDefs

         IF(ASSOCIATED(PDefs)) THEN
           CALL AllocatePDefinitions(Enew)
           Enew % PDefs = PDefs

           ! All elements in actual mesh are not edges
           Enew % PDefs % pyramidQuadEdge = .FALSE.
           Enew % PDefs % isEdge = .FALSE.

           ! If element is of type tetrahedron and is a p element,
           ! do the Ainsworth & Coyle trick
           IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
            CALL GetRefPElementNodes( Enew,  Enew % TYPE % NodeU, &
                 Enew % TYPE % NodeV, Enew % TYPE % NodeW )
         END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO

    CALL Info( 'SplitMeshEqual', '******** New mesh ********', Level=6 )
    WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )


    ! Information of the new system size, also in parallel
    !----------------------------------------------------------------------
    ParTmp(1) = Mesh % NumberOfNodes
    ParTmp(2) = Mesh % NumberOfBulkElements
    ParTmp(3) = Mesh % NumberOfBoundaryElements
    ParTmp(4) = NewMesh % NumberOfNodes
    ParTmp(5) = NewMesh % NumberOfBulkElements
    ParTmp(6) = NewMesh % NumberOfBoundaryElements

    IF( .FALSE. .AND. ParEnv % PEs > 1 ) THEN
      CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

      CALL Info('SplitMeshEqual','Information on parallel mesh sizes')
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(3),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(4),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(5),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(6),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
    END IF


    CALL CheckTimer('SplitMeshEqual',Delete=.TRUE.)

!
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    CALL UpdateParallelMesh( Mesh, NewMesh )
!
!   If periodic BC given, compute boundary mesh projector:
!   ------------------------------------------------------
!
    DO i = 1,CurrentModel % NumberOfBCs
       k = ListGetInteger(CurrentModel % BCs(i) % Values, 'Periodic BC', Found )
       IF( Found ) THEN
         CurrentModel % BCs(i) % PMatrix => &
             PeriodicProjector( CurrentModel, Mesh, i, k )
       END IF
    END DO
!
!   Finalize:
!   ---------
    DEALLOCATE( Child )
    IF(.NOT.EdgesPresent) THEN
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    ELSE
      CALL FindMeshEdges( NewMesh )
    END IF

!call writemeshtodisk( NewMesh, "." )
!stop
CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelMesh( Mesh, NewMesh )
!------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Edge, Face, Element, BoundaryElement
       INTEGER :: i,j,k,l,m,n,p,q, istat
       INTEGER, POINTER :: IntCnts(:),IntArray(:),Reorder(:)
       INTEGER, ALLOCATABLE :: list1(:), list2(:)
       LOGICAL, ALLOCATABLE :: InterfaceTag(:)

       INTEGER :: jedges
       LOGICAL :: Found
!------------------------------------------------------------------------------

       IF ( ParEnv % PEs <= 1 ) RETURN
!
!      Update mesh interfaces for parallel execution.
!      ==============================================
!
!      Try to get an agreement about the  global numbering
!      of new mesh nodes among set of processes solving
!      this specific eq. Also allocate and generate
!      all other control information needed in parallel
!      execution:
!      ----------------------------------------------------
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) &
         CALL Fatal( 'UpdateParallelMesh', 'Allocate error.' )
       CALL AllocateVector( NewMesh % ParallelInfo % INTERFACE,n  )
       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )

       DO i=1,n
          NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % INTERFACE = .FALSE.
       NewMesh % ParallelInfo % INTERFACE(1:n) = Mesh % ParallelInfo % INTERFACE

       NewMesh % ParallelInfo % GlobalDOFs = 0
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = &
          Mesh % ParallelInfo % GlobalDOFs
!
!      My theory is, that a new node will be an
!      interface node only if all the edge or face
!      nodes which contribute to its existence are
!      interface nodes (the code immediately below
!      will only count sizes):
!      -------------------------------------------
!

       ! New version based on edges and faces (2. March 2007):
       !=====================================================
       SELECT CASE( CoordinateSystemDimension() )
          
       CASE(2)
          !
          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % INTERFACE(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'
          !
          ! Determine possible interface edges:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfEdges ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfEdges
             Edge => Mesh % Edges(i)
             IF( ASSOCIATED(Edge % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Edge % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % INTERFACE( Edge % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          !
          ! Eliminate false positives based on BoundaryElement -data:
          !----------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements( Mesh % NumberOfBulkElements + i )
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED( Element ) ) &
                  Element => BoundaryElement % BoundaryInfo % Right
             IF( .NOT.ASSOCIATED( Element ) ) CYCLE
             IF( .NOT.ASSOCIATED( Element % EdgeIndexes ) ) CYCLE
             
             ALLOCATE( list1( SIZE( BoundaryElement % NodeIndexes )))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort( SIZE(list1), list1 )
             
             DO j = 1,Element % TYPE % NumberOfEdges
                k = Element % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                IF( SIZE( Edge % NodeIndexes ) /= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Edge % NodeIndexes )))
                list2 = Edge % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO

                DEALLOCATE(list2)
                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfEdges
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Edge => Mesh % Edges(i)
             
             ! This is just for the edge count:
             !---------------------------------
             IF( NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + i) ) CYCLE
             
             ! Mark interface nodes and count edges:
             !--------------------------------------
             NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + i) = .TRUE.
             p = p+1

          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 2*p ! check
          
       CASE(3)

          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % INTERFACE(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'

          ! Determine possible interface faces:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfFaces ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( ASSOCIATED(Face % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Face % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % INTERFACE( Face % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          
          ! Eliminate false interface faces based on BoundaryElement -data:
          !----------------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED(Element) ) &
                Element => BoundaryElement % BoundaryInfo % Right
              IF( .NOT.ASSOCIATED(Element) ) CYCLE
              IF( .NOT.ASSOCIATED(Element % FaceIndexes) ) CYCLE
             
             ALLOCATE(list1(SIZE(BoundaryElement % NodeIndexes)))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort(SIZE(list1),list1)
             
             DO j = 1,Element % TYPE % NumberOfFaces
                k = Element % FaceIndexes(j)
                Face => Mesh % Faces(k)
                IF(SIZE(Face % NodeIndexes)/= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Face % NodeIndexes )))
                list2 = Face % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO
                
                DEALLOCATE(list2)

                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Count interface faces:
          !-----------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( InterfaceTag(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface faces'
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Face => Mesh % Faces(i)
             
             DO j = 1,SIZE( Face % EdgeIndexes )
                k = Face % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                
                ! This is just for the edge count:
                !---------------------------------
                IF( NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + k) ) CYCLE
                
                ! Mark interface nodes and count edges:
                !--------------------------------------
                NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + k) = .TRUE.
                p = p+1
             END DO
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 3*p ! check
          
       END SELECT

!======================================================================================================
       j = p
       jedges = p

!      For bricks, check also the faces:
!      ---------------------------------
       DO i = 1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i) 
          IF( Face % TYPE % NumberOfNodes == 4 ) THEN
             IF ( ALL( Mesh % ParallelInfo % INTERFACE( Face % NodeIndexes ) ) ) THEN
                NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes &
                     + Mesh % NumberOfEdges + i ) = .TRUE.
                j = j + 1
                k = k + Face % TYPE % NumberOfNodes
             END IF
          END IF
       END DO

!      print*,'Found',j-jedges,'interface faces'

!      CALL AllocateVector( IntCnts,  j )
!      CALL AllocateVector( IntArray, k )
!
!      Old mesh nodes were copied as is...
!      -----------------------------------
       DO i=1,Mesh % NumberOfNodes
          CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, &
                SIZE( Mesh % ParallelInfo % Neighbourlist(i) % Neighbours) )

          NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO
!
!      Take care of the new mesh internal nodes.
!      Parallel global numbering will take care
!      of the interface nodes:
!      ----------------------------------------
       DO i=Mesh % NumberOfNodes+1, NewMesh % NumberOfNodes
          IF ( .NOT. NewMesh % ParallelInfo % INTERFACE(i) ) THEN
            CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours,1 )
            NewMesh % ParallelInfo % NeighbourList(i) %  Neighbours(1) = ParEnv % MyPE
          END IF
       END DO
!
!      Copy global indices of edge and/or face nodes
!      to temporary work arrays:
!      ---------------------------------------------
!
! check also this:
!      j = 0
!      k = 0
!      DO i = 1,Mesh % NumberOfEdges
!         Edge => Mesh % Edges(i)
!         
!         ! Added check for parent elements 25.2.2007:
!         Found = .NOT.( ASSOCIATED(edge % boundaryinfo % left) &
!              .AND.  ASSOCIATED(edge % boundaryinfo % right) )
!         
!         IF ( ALL(Mesh % ParallelInfo % INTERFACE(Edge % NodeIndexes)) .AND. Found ) THEN
!            j = j + 1
!            IntCnts(j) = Edge % TYPE % NumberOfNodes
!            IntArray( k+1:k+IntCnts(j) ) = &
!                 Mesh % Parallelinfo % GlobalDOFs(Edge % NodeIndexes)
!            CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!            k = k + IntCnts(j)
!         END IF
!      END DO
!      !
!      ! For bricks, check also the faces:
!      ! ---------------------------------
!      DO i = 1,Mesh % NumberOfFaces
!         Face => Mesh % Faces(i)
!         IF( Face % TYPE % NumberOfNodes == 4 ) THEN
!            IF ( ALL( Mesh % ParallelInfo % INTERFACE(Face % NodeIndexes) ) ) THEN
!               j = j + 1
!               IntCnts(j) = Face % TYPE % NumberOfNodes
!               IntArray(k+1:k+IntCnts(j)) = &
!                    Mesh % ParallelInfo % GlobalDOFs(Face % NodeIndexes)
!               CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!               k = k + IntCnts(j)
!            END IF
!         END IF
!      END DO
!
!      Finally the beef, do the exchange of new
!      interfaces. The parallel global numbering
!      subroutine will also do reordering of the
!      nodes, hence the reorder array:
!      -------------------------------------------
       CALL AllocateVector( Reorder, NewMesh % NumberOfNodes )
       Reorder = (/ (i, i=1,NewMesh % NumberOfNodes) /)

       k = NewMesh % Nodes % NumberOfNodes - Mesh % Nodes % NumberOfNodes

       CALL ParallelGlobalNumbering( NewMesh, Mesh, k, IntCnts, IntArray, Reorder )

!      Account for the reordering of the nodes:
!      ----------------------------------------
       DO i=1,NewMesh % NumberOfBulkElements + &
            NewMesh % NumberOfBoundaryElements
          NewMesh % Elements(i) % NodeIndexes = &
              Reorder( NewMesh % Elements(i) % NodeIndexes )
       END DO

!      DEALLOCATE( IntCnts, IntArray, Reorder )
!      DEALLOCATE( Reorder )
!------------------------------------------------------------------------------
    END SUBROUTINE UpdateParallelMesh
  END FUNCTION SplitMeshEqual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMesh( Mesh )
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     TYPE(Projector_t), POINTER :: Projector
     TYPE(Projector_t), POINTER :: Projector1
     TYPE(Variable_t), POINTER  :: Var, Var1
     INTEGER :: i,j,k
     LOGICAL :: GotIt
     REAL(KIND=dp), POINTER :: ptr(:)
!------------------------------------------------------------------------------
 
!    Deallocate mesh variables:
!    --------------------------


     CALL Info('ReleaseMesh','Releasing mesh variables',Level=15)
     CALL ReleaseVariableList( Mesh % Variables )
     Mesh % Variables => NULL()

!    Deallocate mesh geometry (nodes,elements and edges):
!    ----------------------------------------------------
     IF ( ASSOCIATED( Mesh % Nodes ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh nodes',Level=15)
       IF ( ASSOCIATED( Mesh % Nodes % x ) ) DEALLOCATE( Mesh % Nodes % x )
       IF ( ASSOCIATED( Mesh % Nodes % y ) ) DEALLOCATE( Mesh % Nodes % y )
       IF ( ASSOCIATED( Mesh % Nodes % z ) ) DEALLOCATE( Mesh % Nodes % z )
       DEALLOCATE( Mesh % Nodes )

       IF ( ASSOCIATED( Mesh % ParallelInfo % GlobalDOFs ) ) &
           DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs )

       IF ( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN 
         DO i=1,Mesh % NumberOfNodes
           IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
               DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
         END DO
         DEALLOCATE( Mesh % ParallelInfo % NeighbourList )
       END IF

       IF ( ASSOCIATED( Mesh % ParallelInfo % INTERFACE ) ) &
           DEALLOCATE( Mesh % ParallelInfo % INTERFACE )
     END IF

     Mesh % Nodes => NULL()

     IF ( ASSOCIATED( Mesh % Edges ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh edges',Level=15)
       CALL ReleaseMeshEdgeTables( Mesh )
       Mesh % Edges => NULL()
     END IF

     IF ( ASSOCIATED( Mesh % Faces ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh faces',Level=15)
       CALL ReleaseMeshFaceTables( Mesh )
       Mesh % Faces => NULL()
     END IF

     IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
     CALL Info('ReleaseMesh','Releasing mesh view factors',Level=15)
       CALL ReleaseMeshFactorTables( Mesh % ViewFactors )
       Mesh % ViewFactors => NULL()
     END IF


!    Deallocate mesh to mesh projector structures:
!    ---------------------------------------------
     Projector => Mesh % Projector
     DO WHILE( ASSOCIATED( Projector ) )
       CALL Info('ReleaseMesh','Releasing mesh projector',Level=15)
       CALL FreeMatrix( Projector % Matrix )
       CALL FreeMatrix( Projector % TMatrix )
       Projector1 => Projector
       Projector => Projector % Next
       DEALLOCATE( Projector1 )
     END DO
     Mesh % Projector => NULL()


!    Deallocate quadrant tree (used in mesh to mesh interpolation):
!    --------------------------------------------------------------
     IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh quadrant tree',Level=15)
       CALL FreeQuadrantTree( Mesh % RootQuadrant )
       Mesh % RootQuadrant => NULL()
     END IF


     IF ( ASSOCIATED( Mesh % Elements ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh elements',Level=15)

        DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

!          Boundaryinfo structure for boundary elements
!          ---------------------------------------------
           IF ( Mesh % Elements(i) % Copy ) CYCLE

           IF ( i > Mesh % NumberOfBulkElements ) THEN
             IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo ) ) THEN
               IF (ASSOCIATED(Mesh % Elements(i) % BoundaryInfo % GebhardtFactors)) THEN
                 IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo % &
                     GebhardtFactors % Elements ) ) THEN
                   DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                       GebhardtFactors % Elements )
                   DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                       GebhardtFactors % Factors )
                 END IF
                 DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % GebhardtFactors )
               END IF
               DEALLOCATE( Mesh % Elements(i) % BoundaryInfo )
             END IF
           END IF

           IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
               DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
           Mesh % Elements(i) % NodeIndexes => NULL()
           
           IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
           Mesh % Elements(i) % EdgeIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
           Mesh % Elements(i) % FaceIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % DGIndexes )
           Mesh % Elements(i) % DGIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
             DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
           Mesh % Elements(i) % BubbleIndexes => NULL()

           ! This creates problems later on!!!
           !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
           !   DEALLOCATE( Mesh % Elements(i) % PDefs )

           Mesh % Elements(i) % PDefs => NULL()
 
        END DO
        DEALLOCATE( Mesh % Elements )
        Mesh % Elements => NULL()
      END IF

      CALL Info('ReleaseMesh','Releasing mesh finished',Level=15)
     
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshEdgeTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
       DO i=1,Mesh % NumberOfEdges
          Edge => Mesh % Edges(i)
          IF ( ASSOCIATED( Edge % NodeIndexes ) ) THEN
             DEALLOCATE( Edge % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Edge % BoundaryInfo ) ) THEN
             DEALLOCATE( Edge % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Edges )
    END IF
    NULLIFY( Mesh % Edges )
    Mesh % NumberOfEdges = 0

    DO i=1,Mesh % NumberOfBulkElements
       IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) THEN
          DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
          NULLIFY( Mesh % Elements(i) % EdgeIndexes )
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshEdgeTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFaceTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Face
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
       DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF ( ASSOCIATED( Face % NodeIndexes ) ) THEN
             DEALLOCATE( Face % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Face % BoundaryInfo ) ) THEN
             DEALLOCATE( Face % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Faces )
    END IF
    NULLIFY( Mesh % Faces )
    Mesh % NumberOfFaces = 0

    DO i=1,Mesh % NumberOfBulkElements
       IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) THEN
          DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
          NULLIFY( Mesh % Elements(i) % FaceIndexes )
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFaceTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFactorTables( Factors )
!------------------------------------------------------------------------------
    TYPE(Factors_t), POINTER :: Factors(:)
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Factors ) ) THEN
       DO i=1,SIZE( Factors)
          IF (ASSOCIATED(Factors(i) % Factors))  DEALLOCATE(Factors(i) % Factors)
          IF (ASSOCIATED(Factors(i) % Elements)) DEALLOCATE(Factors(i) % Elements)
       END DO
       DEALLOCATE(  Factors )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFactorTables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SetCurrentMesh( Model, Mesh )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),  POINTER :: Mesh
!------------------------------------------------------------------------------
    Model % Variables => Mesh % Variables

    Model % Mesh  => Mesh
    Model % Nodes => Mesh % Nodes
    Model % NumberOfNodes = Mesh % NumberOfNodes
    Model % Nodes % NumberOfNodes = Mesh % NumberOfNodes

    Model % Elements => Mesh % Elements
    Model % MaxElementNodes = Mesh % MaxElementNodes
    Model % NumberOfBulkElements = Mesh % NumberOfBulkElements
    Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE SetCurrentMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE DisplaceMesh( Mesh, Update, SIGN, Perm, DOFs, StabRecomp )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) , POINTER :: Mesh 
    REAL(KIND=dp) :: Update(:)
    INTEGER :: DOFs,SIGN,Perm(:)
    LOGICAL, OPTIONAL :: StabRecomp

    INTEGER :: i,k
    LOGICAL :: StabFlag

    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element

    DO i=1,MIN( SIZE(Perm), SIZE(Mesh % Nodes % x) )
       k = Perm(i)
       IF ( k > 0 ) THEN
         k = DOFs * (k-1)
         Mesh % Nodes % x(i)   = Mesh % Nodes % x(i) + SIGN * Update(k+1)
         IF ( DOFs > 1 ) &
           Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + SIGN * Update(k+2)
         IF ( DOFs > 2 ) &
           Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + SIGN * Update(k+3)
        END IF
    END DO

    StabFlag = .TRUE.
    IF ( PRESENT( StabRecomp ) ) StabFlag = StabRecomp

    IF ( SIGN == 1 .AND. StabFlag ) THEN
       k = Mesh % MaxElementDOFs
       CALL AllocateVector( ElementNodes % x,k )
       CALL AllocateVector( ElementNodes % y,k )
       CALL AllocateVector( ElementNodes % z,k )

       DO i=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          IF ( ANY( Perm( Element % NodeIndexes ) == 0 ) ) CYCLE

          k = Element % TYPE % NumberOfNodes
          ElementNodes % x(1:k) = Mesh % Nodes % x(Element % NodeIndexes)
          ElementNodes % y(1:k) = Mesh % Nodes % y(Element % NodeIndexes)
          ElementNodes % z(1:k) = Mesh % Nodes % z(Element % NodeIndexes)
          IF ( Mesh % Stabilize ) THEN
             CALL StabParam( Element,ElementNodes,k, &
                          Element % StabilizationMk, Element % Hk )
          ELSE
             Element % hK = ElementDiameter( Element, ElementNodes )
          END IF
       END DO

       DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE DisplaceMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Convert tetrahedral element to Ainsworth & Coyle type tetrahedron.
!------------------------------------------------------------------------------
  SUBROUTINE ConvertToACTetra( Tetra )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getTetraEdgeMap, getTetraFaceMap
    IMPLICIT NONE
    
    TYPE(Element_t), POINTER :: Tetra  !< Tetrahedral element to convert
!------------------------------------------------------------------------------
    INTEGER :: i, globalMin, globalMax, globalMinI
    INTEGER, DIMENSION(3) :: face, globalFace
    INTRINSIC MIN, MAX, CSHIFT

    ! Sanity check
    IF (Tetra % TYPE % ElementCode /= 504 .OR. &
         .NOT. ASSOCIATED(Tetra % PDefs)) THEN
       CALL Warn('MeshUtils::ConvertToACTetra','Element to convert not p tetrahedron!')
       RETURN
    END IF    
   
    ! Find global min and max vertices
    globalMin = Tetra % NodeIndexes(1)
    globalMinI = 1
    globalMax = Tetra % NodeIndexes(1)
    DO i=2,4
       ! Find min
       IF (globalMin > Tetra % NodeIndexes(i)) THEN
          globalMin = Tetra % NodeIndexes(i)
          globalMinI = i
       ELSE IF (globalMax < Tetra % NodeIndexes(i)) THEN
          globalMax = Tetra % NodeIndexes(i)
       END IF
    END DO
    
    ! Get face containing global min (either face 1 or 2)
    IF (globalMinI == 4) THEN
       face = getTetraFaceMap(2)
    ELSE
       face = getTetraFaceMap(1)
    END IF
    globalFace(1:3) = Tetra % NodeIndexes(face)

    ! Rotate face until first local index is min global
    DO 
       ! Check if first node matches global min node
       IF (globalMin == globalFace(1)) EXIT
       
       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Face 3 now contains global max
    face = getTetraFaceMap(3)
    globalFace(1:3) = Tetra % NodeIndexes(face)
    ! Rotate face until last local index is max global
    DO 
       ! Chech if last node matches global max node
       IF (globalMax == globalFace(3)) EXIT
       
       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Set AC tetra type
    IF (Tetra % NodeIndexes(2) < Tetra % NodeIndexes(3)) THEN
       Tetra % PDefs % TetraType = 1
    ELSE IF (Tetra % NodeIndexes(3) < Tetra % NodeIndexes(2)) THEN
       Tetra % PDefs % TetraType = 2
    ELSE 
       CALL Fatal('MeshUtils::ConvertToACTetra','Corrupt element type')
    END IF
   
  END SUBROUTINE ConvertToACTetra


!------------------------------------------------------------------------------
!>     Assign local number of edge to given boundary element. Also copies all 
!>     p element attributes from element edge to boundary edge.
!------------------------------------------------------------------------------
  SUBROUTINE AssignLocalNumber( EdgeElement, Element, Mesh )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getFaceEdgeMap 
    IMPLICIT NONE

    ! Parameters
    TYPE(Mesh_t) :: Mesh            !< Finite element mesh containing faces and edges.
    TYPE(Element_t), POINTER :: EdgeElement  !< Edge element to which assign local number
    TYPE(Element_t), POINTER :: Element      !< Bulk element with some global numbering to use to assign local number
!------------------------------------------------------------------------------
    ! Local variables

    INTEGER i,j,n,edgeNumber, numEdges, bMap(4)
    TYPE(Element_t), POINTER :: Edge

    ! Get number of points, edges or faces
    numEdges = 0
    SELECT CASE (Element % TYPE % DIMENSION)
    CASE (2)
       numEdges = Element % TYPE % NumberOfEdges
    CASE (3)   
       numEdges = Element % TYPE % NumberOfFaces
    CASE DEFAULT
       WRITE (*,*) 'MeshUtils::AssignLocalNumber, Unsupported dimension:', Element % TYPE % DIMENSION
       RETURN
    END SELECT

    ! For each edge or face in element try to find local number
    DO edgeNumber=1, numEdges
       ! If edges have not been created, stop search. This should not happen, actually.
       IF (.NOT. ASSOCIATED(Element % EdgeIndexes)) THEN
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       Edge => GetElementEntity(Element,edgeNumber,Mesh)

       ! Edge element not found. This should not be possible, unless there
       ! is an error in the mesh read in process..
       IF (.NOT. ASSOCIATED(Edge)) THEN
          CALL Warn('MeshUtils::AssignLocalNumber','Edge element not found')
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       n = 0
       ! For each element node
       DO i=1, Edge % TYPE % NumberOfNodes
          ! For each node in edge element
          DO j=1, EdgeElement % TYPE % NumberOfNodes
             ! If edge and edgeelement node match increment counter
             IF (Edge % NodeIndexes(i) == EdgeElement % NodeIndexes(j)) n = n + 1
          END DO
       END DO

       ! If all nodes are on boundary, edge was found
       IF (n == EdgeElement % TYPE % NumberOfNodes) THEN
          EdgeElement % PDefs % localNumber = edgeNumber

          ! Change ordering of global nodes to match that of element
          bMap = getElementBoundaryMap( Element, edgeNumber )
          DO j=1,n
          	EdgeElement % NodeIndexes(j) = Element % NodeIndexes(bMap(j))
	  END DO

          ! Copy attributes of edge element to boundary element
          ! Misc attributes
          EdgeElement % PDefs % isEdge = Edge % PDefs % isEdge
          
          ! Gauss points
          EdgeElement % PDefs % GaussPoints = Edge % PDefs % GaussPoints

          ! Element p (and boundary bubble dofs)
          EdgeElement % BDOFs = Edge % BDOFs
          EdgeElement % PDefs % P = Edge % PDefs % P

          ! If this boundary has edges copy edge indexes
          IF (ASSOCIATED(Edge % EdgeIndexes)) THEN
             ! Allocate element edges to element
             n = Edge % TYPE % NumberOfEdges
             bmap(1:4) = getFaceEdgeMap( Element, edgeNumber )
             
             IF ( ASSOCIATED( EdgeElement % EdgeIndexes) ) THEN
                DEALLOCATE( EdgeElement % EdgeIndexes )
             END IF
             
             CALL AllocateVector( EdgeElement % EdgeIndexes, n )
             ! Copy edges from edge to boundary edge
             DO i=1,n
                EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(bmap(i))
             !    EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(i)
             END DO
          END IF
          
          ! Edge fields copied and local edge found so return
          RETURN
       END IF
    END DO

    ! If we are here local number not found
    CALL Warn('MeshUtils::AssignLocalNumber','Unable to find local edge')
    ! EdgeElement % localNumber = 1
  CONTAINS

    FUNCTION GetElementEntity(Element, which, Mesh) RESULT(Entity)
      IMPLICIT NONE

      TYPE(Element_t), POINTER :: Element, Entity 
      INTEGER :: which
      TYPE(Mesh_t) :: Mesh

      NULLIFY(Entity)
      ! Switch by element dimension
      SELECT CASE (Element % TYPE % DIMENSION)
         CASE (2)
            Entity => Mesh % Edges( Element % EdgeIndexes(which))
         CASE (3)
            Entity => Mesh % Faces( Element % FaceIndexes(which))
         CASE DEFAULT
            WRITE (*,*) 'AssignLocalNumber::GetElementEntity: Unsupported dimension'
            RETURN
      END SELECT
    END FUNCTION GetElementEntity
  END SUBROUTINE AssignLocalNumber
    

!------------------------------------------------------------------------------
!>     Based on element degrees of freedom, return the sum of element
!>     degrees of freedom.
!------------------------------------------------------------------------------
  FUNCTION getElementMaxDOFs( Mesh, Element ) RESULT(dofs)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh        !< Finite element mesh
    TYPE(Element_t), POINTER :: Element  !< Element to get maximum dofs for
    INTEGER :: dofs                      !< maximum number of dofs for Element
!------------------------------------------------------------------------------

    TYPE(ELement_t), POINTER :: Edge, Face
    INTEGER :: i, edgeDofs, faceDofs
    
    ! Get sum of edge dofs if any
    edgeDofs = 0
    IF (ASSOCIATED(Element % EdgeIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Element % EdgeIndexes(i))
          edgeDofs = edgeDofs + Edge % BDOFs
       END DO
    END IF

    ! Get sum of face dofs if any
    faceDofs = 0
    IF (ASSOCIATED(Element % FaceIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfFaces
          Face => Mesh % Faces(Element % FaceIndexes(i))
          faceDofs = faceDofs + Face % BDOFs
       END DO
    END IF

    ! Get sum of all dofs in element
    dofs = Element % TYPE % NumberOfNodes + &
         edgeDofs + faceDofs + Element % BDOFs
  END FUNCTION getElementMaxDOFs




!------------------------------------------------------------------------------
!> Creates a permutation table for bodies or boundaries using a free chosen string
!> as mask. The resulting permutation is optimized in order, if requested. The
!> subroutine is intended to help in saving boundary data in an ordered manner,
!> but it can find other uses as well. Currently the implementation is limited
!> to normal Lagrangian elements.
!------------------------------------------------------------------------------
  SUBROUTINE MakePermUsingMask( Model,Solver,Mesh,MaskName, &
       OptimizeBW, Perm, LocalNodes, MaskOnBulk, RequireLogical )
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Mesh_t)   :: Mesh
    TYPE(SOlver_t) :: Solver
    INTEGER :: LocalNodes
    LOGICAL :: OptimizeBW
    INTEGER, POINTER :: Perm(:)
    CHARACTER(LEN=*) :: MaskName
    LOGICAL, OPTIONAL :: MaskOnBulk
    LOGICAL, OPTIONAL :: RequireLogical
!------------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm(:)
    TYPE(ListMatrix_t), POINTER :: ListMatrix(:)
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2,This,bf_id   
    LOGICAL :: Flag, Found, FirstRound, MaskIsLogical, Hit
    INTEGER :: Indexes(30), ElemStart, ElemFin, Width
    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr
    TYPE(Element_t), POINTER :: CurrentElement,Elm
    REAL(KIND=dp) :: MinDist, Dist
!------------------------------------------------------------------------------
    
    ! First check if there are active elements for this mask
    IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .FALSE.
    IF( PRESENT( RequireLogical ) ) THEN
      MaskIsLogical = RequireLogical
    ELSE
      MaskIsLogical = .FALSE.
    END IF

    IF(.NOT. ASSOCIATED( Perm ) ) THEN
      ALLOCATE( Perm( Mesh % NumberOfNodes ) )
      Perm = 0
    END IF

    ElemStart = HUGE(ElemStart) 
    ElemFin = 0     
    DO l = 1, Model % NumberOfBodyForces
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical( Model % BodyForces(l) % Values,MaskName,Found) 
       ELSE
         Hit = ListCheckPresent( Model % BodyForces(l) % Values,MaskName)
       END IF 
       IF( Hit ) THEN
          ElemStart = 1
          ElemFin = Mesh % NumberOfBulkElements
          IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .TRUE.
          EXIT
       END IF
    END DO
    DO l = 1, Model % NumberOfBCs
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical(Model % BCs(l) % Values,MaskName,Found )
       ELSE
         Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName )
       END IF
       IF( Hit ) THEN
          ElemStart = MIN( ElemStart, Mesh % NumberOfBulkElements + 1)
          ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
          EXIT
       END IF
    END DO
    
    IF( ElemFin - ElemStart <= 0) THEN
       LocalNodes = 0
       RETURN
    END IF


    k = 0
    Perm = 0
    FirstRound = .TRUE.

    ! Loop over the active elements
    ! 1st round initial numbering is given
    ! 2nd round a list matrix giving all the connections is created

100 DO t=ElemStart, ElemFin
       
       CurrentElement => Mesh % Elements(t)
       
       Hit = .FALSE.
       IF(t <= Mesh % NumberOfBulkElements) THEN
          l = CurrentElement % BodyId
	  bf_id = ListGetInteger( Model % Bodies(l) % Values, 'Body Force',Found)
	  IF( bf_id>0 ) THEN
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical( Model % BodyForces(bf_id) % Values, MaskName, Found )
            ELSE
              Hit = ListCheckPresent( Model % BodyForces(bf_id) % Values, MaskName )
            END IF
	  END IF 
       ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical(Model % BCs(l) % Values,MaskName, Found ) 
            ELSE
              Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName ) 
            END IF
            EXIT
          END DO
       END IF       
       IF( .NOT. Hit ) CYCLE       
       
       n = CurrentElement % NDOFs               
       Indexes(1:n) = CurrentElement % NodeIndexes(1:n)
       
       IF( FirstRound ) THEN
          DO i=1,n
             j = Indexes(i)
             IF ( Perm(j) == 0 ) THEN
                k = k + 1
                Perm(j) = k
             END IF
          END DO
       ELSE
          DO i=1,n
             k1 = Perm(Indexes(i))
             IF ( k1 <= 0 ) CYCLE
             DO j=1,n
                k2 = Perm(Indexes(j))
                IF ( k2 <= 0 ) CYCLE
                Lptr => List_GetMatrixIndex( ListMatrix,k1,k2 )
             END DO
          END DO
       END IF
    END DO
    LocalNodes = k

    ! Don't optimize bandwidth for parallel cases
    IF( ParEnv % PEs > 1 .OR. .NOT. OptimizeBW ) RETURN

    IF(FirstRound) THEN
       ! Allocate space 
       NULLIFY( ListMatrix )
       ALLOCATE( ListMatrix(LocalNodes) )
       DO i=1,LocalNodes
          ListMatrix(i) % Degree = 0
          NULLIFY( ListMatrix(i) % Head )
       END DO
       FirstRound = .FALSE.

       ! Find the node in the lower left corner at give it the 1st index
       ! since it will probably determine the 1st index
       MinDist = HUGE(MinDist)
       DO i=1,SIZE(Perm)
          IF( Perm(i) <= 0) CYCLE
          Dist = Mesh % Nodes % x(i) + Mesh % Nodes % y(i) + Mesh % Nodes % z(i)
          IF(Dist < MinDist) THEN
             MinDist = Dist
             j = i
          END IF
       END DO

       ! Find the 1st node and swap it with the lower corner
       DO i=1,SIZE(Perm)
          IF( Perm(i) == 1) EXIT
       END DO       
       Perm(i) = Perm(j)
       Perm(j) = 1

       GOTO 100
    END IF

!------------------------------------------------------------------------------

    ALLOCATE( InvPerm(LocalNodes) )
    InvPerm = 0
    DO i=1,SIZE(Perm)
       IF (Perm(i)>0) InvPerm(Perm(i)) = i
    END DO

    ! The bandwidth optimization for lines results to perfectly ordered 
    ! permutations. If there is only one line the 1st node should be the 
    ! lower left corner.

    Flag = .TRUE.
    Width = OptimizeBandwidth( ListMatrix, Perm, InvPerm, &

         LocalNodes, Flag, Flag, MaskName )

    ! We really only need the permutation, as there will be no matrix equation
    ! associated with it.
    DEALLOCATE( InvPerm )
    CALL List_FreeMatrix( LocalNodes, ListMatrix )

!------------------------------------------------------------------------------
  END SUBROUTINE MakePermUsingMask
!------------------------------------------------------------------------------




!------------------------------------------------------------------------
!> Find a point in the mesh structure
!> There are two strategies:
!> 1) Recursive where the same routine is repeated with sloppier criteria
!> 2) One-sweep strategy where the best hit is registered and used if of 
!>    acceptable accuracy. 
!> There are two different epsilons that control the search. One for the 
!> rough test in absolute coordinates and another one for the more accurate
!> test in local coordinates.   
!-------------------------------------------------------------------------
  FUNCTION PointInMesh(Solver, GlobalCoords, LocalCoords, HitElement, &
      CandElement, ExtInitialize ) RESULT ( Hit )
        
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: GlobalCoords(3), LocalCoords(3)
    TYPE(Element_t), POINTER :: HitElement 
    TYPE(Element_t), POINTER, OPTIONAL :: CandElement
    LOGICAL, OPTIONAL :: ExtInitialize
    LOGICAL :: Hit
!-------------------------------------------------------------------------
    LOGICAL :: Initialize, Allocated = .FALSE., Stat, DummySearch, &
        MaskExists, Found, IsRecursive
    INTEGER :: i,j,k,n,bf_id,dim,mini
    REAL(KIND=dp) :: u,v,w,dist,mindist,MinLocalCoords(3)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>NULL(), LeafQuadrant
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1 = 1e-3, GlobalEps, LocalEps
    CHARACTER(LEN=MAX_NAME_LEN) :: MaskName


    SAVE :: Allocated, ElementNodes, DummySearch, Mesh, MaskName, MaskExists, &
        GlobalEps, LocalEps, IsRecursive


    IF( PRESENT( ExtInitialize ) ) THEN
      Initialize = ExtInitialize
    ELSE
      Initialize = .NOT. Allocated 
    END IF

    IF( Initialize ) THEN
      Mesh => Solver % Mesh
      n = Mesh % MaxElementNodes
      IF( Allocated ) THEN
        DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
      END IF
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n))
      Allocated = .TRUE.

      IsRecursive = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Recursive',Stat )
!      IF(.NOT. Stat ) IsRecursive = .TRUE.

      LocalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Local Epsilon', Stat )
      IF(.NOT. stat) LocalEps = 1.0e-10

      GlobalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Global Epsilon', Stat ) 
      IF(.NOT. stat) THEN
        IF( IsRecursive ) THEN
          GlobalEps = 2.0e-10
        ELSE
          GlobalEps = 1.0e-4
        END IF
      END IF

      DummySearch = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Dummy',Stat )

      MaskName = ListGetString( CurrentModel % Simulation,&
          'Interpolation Search Mask',MaskExists )

      IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
        CALL FreeQuadrantTree( Mesh % RootQuadrant )
        Mesh % RootQuadrant => NULL()
      END IF
    END IF
      

    !-----------------------------------------------
    ! Create the octree search structure, if needed 
    !-----------------------------------------------
    IF ( .NOT. ( DummySearch .OR.  ASSOCIATED( Mesh % RootQuadrant ) ) ) THEN
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MINVAL( Mesh % Nodes % y )
      BoundingBox(3) = MINVAL( Mesh % Nodes % z )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
      
      eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
      BoundingBox(1:3) = BoundingBox(1:3) - eps2
      BoundingBox(4:6) = BoundingBox(4:6) + eps2
      
      CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)
      RootQuadrant => Mesh % RootQuadrant
      IF (.NOT. ASSOCIATED(RootQuadrant) ) THEN
        Hit = .FALSE.
        CALL Warn('PointInMesh','No RootQuadrant associated')
        RETURN
      END IF
    END IF


    Hit = .FALSE.

    ! Check that the previous hit is not hit even now
    !-------------------------------------------------
    IF( PRESENT( CandElement ) ) THEN

      IF( ASSOCIATED(CandElement)) THEN

        CurrentElement => CandElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        IF ( PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords ) ) THEN
          Hit = .TRUE.
          HitElement => CurrentElement
          RETURN
        END IF
      END IF
    END IF


    Eps1 = GlobalEps
    Eps2 = LocalEps


100 IF( DummySearch ) THEN

      mindist = HUGE( mindist ) 
      
      !----------------------------------------------------------
      ! Go through all bulk elements in a dummy search.
      ! This algorithm is mainly here for debugging purposes, or
      ! if just a few nodes need to be searched.
      !----------------------------------------------------------
      DO k=1,Mesh % NumberOfBulkElements
        CurrentElement => Mesh % Elements(k)
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        IF( MaskExists ) THEN
          bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
              'Body Force', Found )
          IF( .NOT. Found ) CYCLE
          IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
        END IF

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        Hit = PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist )
        IF( dist < mindist ) THEN
          mini = k
          mindist = dist
        END IF
        IF( Hit ) EXIT
      END DO      
    ELSE
      !-----------------------------------------------
      ! Find the right element using an octree search
      ! This is the preferred algorithms of the two.
      !-----------------------------------------------
      NULLIFY(CurrentElement)
      CALL FindLeafElements(GlobalCoords, Mesh % MeshDim, RootQuadrant, LeafQuadrant)
      IF ( ASSOCIATED(LeafQuadrant) ) THEN
        DO j=1, LeafQuadrant % NElemsInQuadrant
          k = LeafQuadrant % Elements(j)
          CurrentElement => Mesh % Elements(k)
          
          IF( MaskExists ) THEN
            bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
                'Body Force', Found )
            IF( .NOT. Found ) CYCLE
            IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
          END IF
          
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
                    
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
          
          Hit = PointInElement( CurrentElement, ElementNodes, &
              GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist ) 
          IF( dist < mindist ) THEN
            mini = k
            mindist = dist
            MinLocalCoords = LocalCoords
          END IF
          IF( Hit ) EXIT
        END DO
      END IF      
    END IF

    IF( .NOT. Hit ) THEN
      IF( IsRecursive ) THEN
        Eps1 = 10.0 * Eps1
        Eps2 = 10.0 * Eps2
        IF( Eps1 <= 1.0_dp ) GOTO 100
      ELSE
        IF( mindist < Eps1 ) THEN
          CurrentElement => Mesh % Elements(k)
          LocalCoords = MinLocalCoords
          Hit = .TRUE.
        END IF
      END IF
    END IF

    IF( Hit ) HitElement => CurrentElement
    
  END FUNCTION PointInMesh



!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh even though it is 
!> given in an unstructured format. The routine may be used by some special
!> solvers that employ the special character of the mesh.
!> The extrusion is found for a given direction and for each node the corresponding 
!> up and down, and thereafter top and bottom node is computed.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedStructure( Mesh, Solver, ExtVar, &
      TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, &
      MidNodePointer, MidLayerExists, NumberOfLayers, NodeLayer )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopNodePointer(:), BotNodePointer(:), &
        UpNodePointer(:), DownNodePointer(:), MidNodePointer(:)
    INTEGER, POINTER, OPTIONAL :: NodeLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
    LOGICAL, OPTIONAL :: MidLayerExists
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
	UpHit, DownHit, bc_ind
    INTEGER, POINTER :: NodeIndexes(:), MaskPerm(:)
    LOGICAL :: MaskExists, UpActive, DownActive, GotIt, Found, DoCoordTransform
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
#ifndef USE_ISO_C_BINDINGS
    REAL(KIND=dp) :: CPUTime
#endif
    REAL(KIND=dp) :: at0, at1, Length, UnitVector(3), Vector(3), Vector2(3), &
                 ElemVector(3), DotPro, Eps, MinTop, MaxTop, MinBot, MaxBot
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CoordTransform

   
    CALL Info('DetectExtrudedStructure','Determining extruded structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim
    Params => Solver % Values

    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal('StructuredMeshMapper','Invalid value for Active Coordinate')
    END IF  
    UnitVector = 0.0_dp
    UnitVector(ActiveDirection) = 1.0_dp


    IF( ListGetLogical(Params,'Project To Bottom',GotIt) ) &
        UnitVector = -1.0_dp * UnitVector

    WRITE(Message,'(A,3F8.3)') 'Unit vector of direction:',UnitVector
    CALL Info('DetectExtrudedStructure',Message,Level=8)

    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0e-4_dp

    VarName = ListGetString(Params,'Mapping Mask Variable',GotIt )
    MaskExists = .FALSE.
    IF(GotIt) THEN
      Var => VariableGet( Mesh % Variables,  VarName )
      IF(ASSOCIATED(Var)) THEN
        MaskExists = ASSOCIATED(Var % Perm)
        IF( MaskExists ) THEN
          ALLOCATE( MaskPerm( SIZE( Var % Perm ) ) )
          MaskPerm = Var % Perm 
          CALL Info('DetectExtrudedStructure',&
              'Using variable as mask: '//TRIM(VarName),Level=8)
        END IF
      END IF
    END IF

    IF( MaskExists ) THEN
      nsize = MAXVAL( MaskPerm ) 
      WRITE(Message,'(A,I8)') 'Applying mask of size:',nsize
      CALL Info('DetectExtrudedStructure',Message,Level=8)
    ELSE
      nsize = Mesh % NumberOfNodes
      CALL Info('DetectExtrudedStructure','Applying mask to the whole mesh',Level=8)
    END IF 

    CoordTransform = ListGetString(Params,'Mapping Coordinate Transformation',DoCoordTransform )
    IF( DoCoordTransform .OR. MaskExists) THEN
      NULLIFY( Values )
      ALLOCATE( Values( nsize ) )
      Values = 0.0_dp
      IF( MaskExists ) THEN
        CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values, MaskPerm)
      ELSE
        CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values)
      END IF
      Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
    ELSE IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF	      

    IF( MaskExists .OR. DoCoordTransform) THEN
      DO i=1,Mesh % NumberOfNodes
        j = i
	IF( MaskExists ) j = MaskPerm(i)
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
      CALL Warn('DetectExtrudedStructure','Either up or down direction should be active')
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


    CALL Info('DetectExtrudedStructure','determine up and down pointers',Level=9)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfBulkElements      

      Element => Mesh % Elements(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

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

          IF(DotPro > 1.0_dp - Eps) THEN 
            IF( UpActive ) UpPointer(ii) = jj
            IF( DownActive ) DownPointer(jj) = ii
          ELSE IF(DotPro < Eps - 1.0_dp) THEN
            IF( DownActive ) DownPointer(ii) = jj
            IF( UpActive ) UpPointer(jj) = ii
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info('DetectExtrudedStructure','determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0) CYCLE
        END IF
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
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO
    ! The last round is always a check
    Rounds = Rounds - 1

    WRITE( Message,'(A,I0,A)') 'Layered structure detected in ',Rounds,' cycles'
    CALL Info('DetectExtrudedStructure',Message,Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info('DetectExtrudedStructure','Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal('DetectExtrudedStructure','Zero rounds implies unsuccesfull operation')
    END IF


    ! Compute the number of layers. The Rounds above may in some cases 
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info('DetectExtrudedStructure','compute the number of layers',Level=9)    
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        EXIT
      END DO

      j = BotPointer(i)

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
        CALL Warn('DetectExtrudedStructure', Message )
        NumberOfLayers = Rounds
      END IF
      WRITE(Message,'(A,I0)') 'Extruded structure layers: ',NumberOfLayers
      CALL Info('DetectExtrudedStructure',Message)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( NodeLayer ) ) THEN
      CALL Info('DetectExtrudedStructure','creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      IF( MaskExists ) THEN
        WHERE( MaskPerm == 0 ) Layer = 0
      END IF
      
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
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
      
      NodeLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info('DetectExtrudedStructure',Message)
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
        CALL Info('DetectExtrudedStructure','determine mid pointers',Level=9)

        DO Rounds = 1, nsize
          DownHit = 0
          UpHit = 0
          DO i=1,nsize
            IF( MaskExists ) THEN
              IF( MaskPerm(i) == 0) CYCLE
            END IF
            IF( MidPointer(i) == 0 ) CYCLE
            IF( UpActive ) THEN
              j = UpPointer(i)
              IF( MidPointer(j) == 0 ) THEN
                UpHit = UpHit + 1
                MidPointer(j) = MidPointer(i)
              END IF
            END IF
            IF( DownActive ) THEN
              j = DownPointer(i)
              IF( MidPointer(j) == 0 ) THEN
                DownHit = DownHit + 1
                MidPointer(j) = MidPointer(i)
              END IF
            END IF
          END DO
          IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
        END DO

        CALL Info('DetectExtrudedStructure',&
            'Mid layer structure detected in '//TRIM(I2S(Rounds-1))//' cycles',Level=9)
        MidNodePointer => MidPointer
      ELSE
        DEALLOCATE( MidPointer ) 
        MidNodePointer => NULL()
      END IF
    END IF

  
    ! Count the number of top and bottom nodes, for information only
    !---------------------------------------------------------------
    CALL Info('DetectExtrudedStructure','counting top and bottom bodes',Level=9)        
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
    CALL Info('DetectExtrudedStructure','Setting pointer structures',Level=9)        
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
    CALL Info('DetectExtrudedStructure',Message)
    WRITE(Message,* ) 'Top and bottom pointer init rounds: ',Rounds
    CALL Info('DetectExtrudedStructure',Message)
    IF( UpActive ) THEN
      WRITE(Message,* ) 'Number of nodes at the top: ',TopNodes
      CALL Info('DetectExtrudedStructure',Message)
    END IF
    IF( DownActive ) THEN
      WRITE(Message,* ) 'Number of nodes at the bottom: ',BotNodes
      CALL Info('DetectExtrudedStructure',Message)
    END IF


  CONTAINS
    
    
    !---------------------------------------------------------------
    SUBROUTINE CoordinateTransformationNodal( CoordTransform, R )
      CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
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

  !----------------------------------------------------------------
  !> Maps coordinates from the original nodes into a new coordinate
  !> system while optionally maintaining the original coordinates. 
  !> Note that this may be called 
  !---------------------------------------------------------------
  SUBROUTINE CoordinateTransformation( Mesh, CoordTransform, Params, &
      IrreversibleTransformation )
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL, OPTIONAL :: IrreversibleTransformation
    !---------------------------------------------------------------   
    REAL(KIND=dp) :: R0(3),R1(3),Coeff,Rad0
    LOGICAL :: Irreversible,FirstTime,Reuse,UpdateNodes,Found
    REAL(KIND=dp), POINTER :: x0(:),y0(:),z0(:),x1(:),y1(:),z1(:),&
        NewCoords(:)
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
        Mesh % NodesMapped % x => NewCoords(1::3)
        Mesh % NodesMapped % y => NewCoords(2::3)
        Mesh % NodesMapped % z => NewCoords(3::3)
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
              'Creating variables for > Tranformed Coordinate < ')
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

    CALL Info('CoordinateTransformation','All done',Level=8)

  END SUBROUTINE CoordinateTransformation
!---------------------------------------------------------------


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
!---------------------------------------------------------------


!---------------------------------------------------------------
!> This partitions the mesh into a given number of partitions in each 
!> direction. It may be used in clustering multigrid or similar, 
!> and also to internal partitioning within ElmerSolver. 
!---------------------------------------------------------------
  SUBROUTINE ClusterNodesByDirection(Params,Mesh,Clustering,MaskActive)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, POINTER :: Clustering(:)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,GotIt,Hit
    REAL(KIND=dp), ALLOCATABLE :: Measure(:)
    INTEGER :: i,j,k,k0,l,ind,n,dim,dir,divs,nsize,elemsinpart,clusters
    INTEGER, POINTER :: Iarray(:),Order(:),NodePart(:),NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,clustersize
    REAL(KIND=dp), POINTER :: PArray(:,:), Arrange(:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), Weights(3), &
        avepart,devpart
!---------------------------------------------------------------

    ! CALL Info('ClusterNodesByDirection','')

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive )
    ELSE
      nsize = Mesh % NumberOfNodes
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('ClusterNodesByDirection','No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0e-2
      IF( dim == 3) Normal(3) = 1.0e-4
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )

    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    

    IF( .FALSE. ) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent1
      PRINT *,'Tangent2:',Tangent2
    END IF


    Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',GotIt )
    IF(.NOT. GotIt) Iarray => ListGetIntegerArray( Params,'MG Cluster Divisions',GotIt )
    Divisions = 1
    IF( GotIt ) THEN
      n = MIN( SIZE(Iarray), dim ) 
      Divisions(1:n) = Iarray(1:n)
    ELSE
      clustersize = ListGetInteger( Params,'Partitioning Size',GotIt)
      IF(.NOT. GotIt) clustersize = ListGetInteger( Params,'MG Cluster Size',GotIt)
      IF( GotIt .AND. ClusterSize > 0) THEN
        IF( dim == 2 ) THEN
          Divisions(1) = ( nsize / clustersize ) ** 0.5
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal('ClusterNodesByDirection','Clustering Divisions not given!')
      END IF
    END IF

    Clusters = Divisions(1) * Divisions(2) * Divisions(3)

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(Order(nsize),Arrange(nsize),NodePart(nsize),NoPart(Clusters))
    

    ! These are needed as an initial value for the loop over dimension
    elemsinpart = nsize
    nodepart = 1
    

    ! Go through each direction and cumulatively add to the clusters
    !-----------------------------------------------------------

    DO dir = 1,dim      
      divs = Divisions(dir)
      IF( divs <= 1 ) CYCLE
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF

      ! Initialize ordering for the current direction
      !----------------------------------------------
      DO i=1,nsize
        Order(i) = i
      END DO
      

      ! Now compute the weights for each node
      !----------------------------------------
      DO i=1,Mesh % NumberOfNodes
        j = i
        IF( MaskExists ) THEN
          IF( .NOT. MaskActive(j) ) CYCLE
        END IF
        
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)

        Arrange(j) = SUM( Weights * Coord )
      END DO

      ! Order the nodes for given direction
      !----------------------------------------------
      CALL SortR(nsize,Order,Arrange)

      ! For each direction the number of elements in cluster becomes smaller
      elemsinpart = elemsinpart / divs

      ! initialize the counter partition
      nopart = 0


      ! Go through each node and locate it to a cluster taking into consideration
      ! the previous clustering (for 1st direction all one)
      !------------------------------------------------------------------------
      j = 1
      DO i = 1,nsize
        ind = Order(i)
        
        ! the initial partition offset depends on previous partitioning
        k0 = (nodepart(ind)-1) * divs

        ! Find the correct new partitioning, this loop is just long enough
        DO l=1,divs
          Hit = .FALSE.
          
          ! test for increase of local partition
          IF( j < divs ) THEN
            IF( nopart(k0+j) >= elemsinpart ) THEN
              j = j + 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! test for decrease of local partition
          IF( j > 1 )  THEN            
            IF( nopart(k0+j-1) < elemsinpart ) THEN
              j = j - 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! If either increase or decrease is needed, this must be ok 
          IF(.NOT. Hit) EXIT
        END DO
          
        k = k0 + j
        nopart(k) = nopart(k) + 1
        nodepart(ind) = k
      END DO

    END DO


    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nsize / clusters
    devpart = 0.0_dp
    DO i=1,clusters
      minpart = MIN( minpart, nopart(i))
      maxpart = MAX( maxpart, nopart(i))
      devpart = devpart + ABS ( nopart(i) - avepart )
    END DO
    devpart = devpart / clusters

    WRITE(Message,'(A,T25,I10)') 'Min nodes in cluster:',minpart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info('ClusterNodesByDirection',Message)
    

    IF( ASSOCIATED(Clustering)) THEN
      Clustering = Nodepart 
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    DEALLOCATE(Order,Arrange,NoPart)


  END SUBROUTINE ClusterNodesByDirection



  SUBROUTINE ClusterElementsByDirection(Params,Mesh,Clustering,MaskActive)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, POINTER :: Clustering(:)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,GotIt,Hit
    REAL(KIND=dp), ALLOCATABLE :: Measure(:)
    INTEGER :: i,j,k,k0,l,ind,n,dim,dir,divs,nsize,elemsinpart,clusters
    INTEGER, POINTER :: Iarray(:),Order(:),NodePart(:),NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,clustersize
    REAL(KIND=dp), POINTER :: PArray(:,:), Arrange(:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), Weights(3), &
        avepart,devpart
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
!---------------------------------------------------------------

    ! CALL Info('ClusterElementsByDirection','')

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive ) 
    ELSE
      nsize = Mesh % NumberOfBulkElements
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('ClusterElementsByDirection','No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0e-2
      IF( dim == 3) Normal(3) = 1.0e-4
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )

    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    
    IF( .FALSE. ) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent1
      PRINT *,'Tangent2:',Tangent2
    END IF

    Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',GotIt )
    IF(.NOT. GotIt ) THEN
      Iarray => ListGetIntegerArray( Params,'MG Cluster Divisions',GotIt )
    END IF

    Divisions = 1
    IF( GotIt ) THEN
      n = MIN( SIZE(Iarray), dim ) 
      Divisions(1:n) = Iarray(1:n)
    ELSE
      clustersize = ListGetInteger( Params,'Partitioning Size',GotIt)
      IF(.NOT. GotIt) clustersize = ListGetInteger( Params,'MG Cluster Size',GotIt)
      IF( GotIt .AND. ClusterSize > 0) THEN
        IF( dim == 2 ) THEN
          Divisions(1) = ( nsize / clustersize ) ** 0.5
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal('ClusterNodesByDirection','Clustering Divisions not given!')
      END IF
    END IF

    Clusters = Divisions(1) * Divisions(2) * Divisions(3)

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(Order(nsize),Arrange(nsize),NodePart(nsize),NoPart(Clusters))
    

    ! These are needed as an initial value for the loop over dimension
    elemsinpart = nsize
    nodepart = 1
    

    ! Go through each direction and cumulatively add to the clusters
    !-----------------------------------------------------------

    DO dir = 1,dim      
      divs = Divisions(dir)
      IF( divs <= 1 ) CYCLE
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF

      ! Initialize ordering for the current direction
      !----------------------------------------------
      DO i=1,nsize
        Order(i) = i
      END DO
      

      ! Now compute the weights for each node
      !----------------------------------------
      DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        IF( MaskExists ) THEN
          IF( .NOT. MaskActive( i ) ) CYCLE
        ELSE 
          IF( i > Mesh % NumberOfBulkElements ) EXIT
        END IF
        
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes 
        n = Element % TYPE % NumberOfNodes

        Coord(1) = SUM( Mesh % Nodes % x( NodeIndexes ) ) / n
        Coord(2) = SUM( Mesh % Nodes % y( NodeIndexes ) ) / n
        Coord(3) = SUM( Mesh % Nodes % z( NodeIndexes ) ) / n

        Arrange(i) = SUM( Weights * Coord )
      END DO

      ! Order the nodes for given direction
      !----------------------------------------------
      CALL SortR(nsize,Order,Arrange)

      ! For each direction the number of elements in cluster becomes smaller
      elemsinpart = elemsinpart / divs

      ! initialize the counter partition
      nopart = 0


      ! Go through each node and locate it to a cluster taking into consideration
      ! the previous clustering (for 1st direction all one)
      !------------------------------------------------------------------------
      j = 1
      DO i = 1,nsize
        ind = Order(i)
        
        ! the initial partition offset depends on previous partitioning
        k0 = (nodepart(ind)-1) * divs

        ! Find the correct new partitioning, this loop is just long enough
        DO l=1,divs
          Hit = .FALSE.
          
          ! test for increase of local partition
          IF( j < divs ) THEN
            IF( nopart(k0+j) >= elemsinpart ) THEN
              j = j + 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! test for decrease of local partition
          IF( j > 1 )  THEN            
            IF( nopart(k0+j-1) < elemsinpart ) THEN
              j = j - 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! If either increase or decrease is needed, this must be ok 
          IF(.NOT. Hit) EXIT
        END DO
          
        k = k0 + j
        nopart(k) = nopart(k) + 1
        nodepart(ind) = k
      END DO

    END DO


    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nsize / clusters
    devpart = 0.0_dp
    DO i=1,clusters
      minpart = MIN( minpart, nopart(i))
      maxpart = MAX( maxpart, nopart(i))
      devpart = devpart + ABS ( nopart(i) - avepart )
    END DO
    devpart = devpart / clusters

    WRITE(Message,'(A,T25,I10)') 'Min nodes in cluster:',minpart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info('ClusterElementsByDirection',Message)
    
    
    IF( ASSOCIATED(Clustering)) THEN
      Clustering = Nodepart 
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    DEALLOCATE(Order,Arrange,NoPart)


  END SUBROUTINE ClusterElementsByDirection



  SUBROUTINE ClusterElementsUniform(Params,Mesh,Clustering,MaskActive)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, POINTER :: Clustering(:)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,Found
    INTEGER :: i,j,k,ind,n,dim,nsize,nmask,clusters
    INTEGER, POINTER :: Iarray(:),NodePart(:)
    INTEGER, ALLOCATABLE :: NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,Inds(3)
    REAL(KIND=dp) :: Coord(3), Weights(3), avepart,devpart
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: BoundingBox(6)
    INTEGER, ALLOCATABLE :: CellCount(:,:,:)
    

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = SIZE( MaskActive ) 
      nmask = COUNT( MaskActive ) 
      CALL Info('ClusterElementsByDirection','Using mask of size: '//TRIM(I2S(nsize)))
    ELSE
      nsize = Mesh % NumberOfBulkElements 
      nmask = nsize
      CALL Info('ClusterElementsByDirection','Applying division to all bulk elements: '//TRIM(I2S(nsize)))
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('ClusterElementsByDirection','No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    BoundingBox = 0.0_dp
    BoundingBox(1) = MINVAL( Mesh % Nodes % x )
    BoundingBox(2) = MAXVAL( Mesh % Nodes % x )
    BoundingBox(3) = MINVAL( Mesh % Nodes % y )
    BoundingBox(4) = MAXVAL( Mesh % Nodes % y )
    BoundingBox(5) = MINVAL( Mesh % Nodes % z )
    BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
    
    !PRINT *,'Bounding Box:',BoundingBox

    Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',Found)
    IF(.NOT. Found ) THEN
      CALL Fatal('ClusterNodesByDirection','> Partitioning Divisions < not given!')
    END IF

    Divisions = 1
    IF( Found ) THEN
      n = MIN( SIZE(Iarray), dim ) 
      Divisions(1:n) = Iarray(1:n)
    END IF

    ALLOCATE( CellCount(Divisions(1), Divisions(2), Divisions(3) ) )
    CellCount = 0
    Clusters = 1
    DO i=1,dim
      Clusters = Clusters * Divisions(i)
    END DO

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(NodePart(nsize),NoPart(Clusters))
    NoPart = 0
    NodePart = 0

    !----------------------------------------
    Inds = 1
    Coord = 0.0_dp

    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      IF( MaskExists ) THEN
        IF( .NOT. MaskActive( i ) ) CYCLE
      ELSE 
        IF( i > Mesh % NumberOfBulkElements ) EXIT
      END IF
      
      Element => Mesh % Elements(i)
      NodeIndexes => Element % NodeIndexes 
      n = Element % TYPE % NumberOfNodes
      
      Coord(1) = SUM( Mesh % Nodes % x( NodeIndexes ) ) / n
      Coord(2) = SUM( Mesh % Nodes % y( NodeIndexes ) ) / n
      IF( dim == 3 ) THEN
        Coord(3) = SUM( Mesh % Nodes % z( NodeIndexes ) ) / n
      END IF

      Inds = 1
      DO j=1,dim
        Inds(j) = CEILING( Divisions(j) * &
            ( Coord(j) - BoundingBox(2*j-1) ) / &
            ( BoundingBox(2*j) - BoundingBox(2*j-1) ) )
      END DO
      Inds = MAX( Inds, 1 ) 

      CellCount(Inds(1),Inds(2),Inds(3)) = &
          CellCount(Inds(1),Inds(2),Inds(3)) + 1

      ind = (Inds(1)-1)*Divisions(2)*Divisions(3) + &
          (Inds(2)-1)*Divisions(3) +  &
          Inds(3)
      NodePart(i) = ind
      NoPart(ind) = NoPart(ind) + 1
    END DO

    ! Compute statistical information of the partitioning
    n = COUNT( NoPart > 0 )    
    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nmask / n
    devpart = 0.0_dp
    DO i=1,clusters
      IF( nopart(i) > 0 ) THEN
        minpart = MIN( minpart, nopart(i))
        maxpart = MAX( maxpart, nopart(i))
        devpart = devpart + ABS ( nopart(i) - avepart )
      END IF
    END DO
    devpart = devpart / n

    WRITE(Message,'(A,T28,I0)') 'Number of partitions:',n
    CALL Info('ClusterElementsUniform',Message)
    WRITE(Message,'(A,T25,I10)') 'Min nodes in cluster:',minpart
    CALL Info('ClusterElementsUniform',Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info('ClusterElementsUniform',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info('ClusterElementsUniform',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info('ClusterElementsUniform',Message)

    ! Renumber the partitions using only the active ones
    n = 0
    DO i=1,clusters
      IF( NoPart(i) > 0 ) THEN
        n = n + 1
        NoPart(i) = n
      END IF
    END DO
    DO i=1,nsize
      j = NodePart(i)
      IF( j > 0 ) NodePart(i) = NoPart(j)
    END DO
           
    IF( ASSOCIATED(Clustering)) THEN
      Clustering(1:nsize) = Nodepart(1:nsize)
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    IF( ALLOCATED( NoPart ) ) DEALLOCATE(NoPart)

    CALL Info('ClusterElemetsUniform','Clustering finished')


  END SUBROUTINE ClusterElementsUniform

 
  !> Find the node closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-----------------------------------------------------------------
  FUNCTION ClosestNodeInMesh(Mesh,Coord,MinDist) RESULT ( NodeIndx )
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coord(3)
    REAL(KIND=dp), OPTIONAL :: MinDist
    INTEGER :: NodeIndx

    REAL(KIND=dp) :: Dist2,MinDist2,NodeCoord(3)
    INTEGER :: i

    MinDist2 = HUGE( MinDist2 ) 

    DO i=1,Mesh % NumberOfNodes
      
      NodeCoord(1) = Mesh % Nodes % x(i)
      NodeCoord(2) = Mesh % Nodes % y(i)
      NodeCoord(3) = Mesh % Nodes % z(i)
    
      Dist2 = SUM( ( Coord - NodeCoord )**2 )
      IF( Dist2 < MinDist2 ) THEN
        MinDist2 = Dist2
        NodeIndx = i  
      END IF
    END DO
    
    IF( PRESENT( MinDist ) ) MinDist = SQRT( MinDist2 ) 

  END FUNCTION ClosestNodeInMesh


  !> Find the element that owns or is closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-------------------------------------------------------------------
  FUNCTION ClosestElementInMesh(Mesh, Coords) RESULT ( ElemIndx )

    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coords(3)
    INTEGER :: ElemIndx

    REAL(KIND=dp) :: Dist,MinDist,LocalCoords(3)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: k,l,n,istat
    REAL(KIND=dp) :: ParallelHits,ParallelCands
    LOGICAL :: Hit

    n = Mesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), STAT=istat)
    IF( istat /= 0 ) CALL Fatal('ClosestElementInMesh','Memory allocation error') 	
    ElemIndx = 0
    MinDist = HUGE( MinDist ) 
    Hit = .FALSE.
    l = 0
    
    ! Go through all bulk elements and look for hit in each element.
    ! Linear search makes only sense for a small number of nodes
    DO k=1,Mesh % NumberOfBulkElements

      Element => Mesh % Elements(k)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
      
      Hit = PointInElement( Element, ElementNodes, &
          Coords, LocalCoords, LocalDistance = Dist )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        l = k
      END IF
      IF( Hit ) EXIT
    END DO
    
    ! Count the number of parallel hits
    !-----------------------------------------------------------------------
    IF( Hit ) THEN
      ParallelHits = 1.0_dp
    ELSE
      ParallelHits = 0.0_dp
    END IF
    ParallelHits = ParallelReduction( ParallelHits )
    
    ! If there was no proper hit go through the best candidates so far and 
    ! see if they would give a acceptable hit
    !----------------------------------------------------------------------
    IF( ParallelHits < 0.5_dp ) THEN	  

      ! Compute the number of parallel candidates
      !------------------------------------------
      IF( l > 0 ) THEN
        ParallelCands = 1.0_dp
      ELSE
        ParallelCands = 0.0_dp
      END IF
      ParallelCands = ParallelReduction( ParallelCands ) 

      IF( l > 0 ) THEN
        Element => Mesh % Elements(l)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        ! If there are more than two competing parallel hits then use more stringent conditions
        ! since afterwords there is no way of deciding which one was closer.
        !--------------------------------------------------------------------------------------
        IF( ParallelCands > 1.5_dp ) THEN
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0e-3_dp, LocalEps=1.0e-4_dp )	
        ELSE
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0_dp, LocalEps=0.1_dp )	
        END IF
      END IF
    END IF

    IF( Hit ) ElemIndx = l

    IF( ParallelHits < 0.5_dp ) THEN
      IF( Hit ) THEN
        ParallelHits = 1.0_dp
      ELSE
        ParallelHits = 0.0_dp
      END IF
      ParallelHits = ParallelReduction( ParallelHits )
      IF( ParallelHits < 0.5_dp ) THEN
        WRITE( Message, * ) 'Coordinate not found in any of the elements!',Coords
        CALL Warn( 'ClosestElementInMesh', Message )
      END IF
    END IF

    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
 
  END FUNCTION ClosestElementInMesh



!---------------------------------------------------------------
!> This find two fixing nodes for each coordinate direction
!> The indexes are returned in order: x1 x2 y1 y2 z1 z2.
!---------------------------------------------------------------
  SUBROUTINE FindRigidBodyFixingNodes(Solver,FixingDofs,MaskPerm)
!------------------------------------------------------------------------------
    USE GeneralUtils

    TYPE(Solver_t) :: Solver
    INTEGER, OPTIONAL :: FixingDofs(0:)
    INTEGER, OPTIONAL :: MaskPerm(:)

!---------------------------------------------------------------

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: MaskExists,FixBestDirection,FoundBetter, GotIt
    INTEGER :: i,j,k,l,ind,n,dim,dir,nsize,Sweep,MaxSweep,DirBest
    INTEGER :: PosMeasureIndex, NegMeasureIndex, FixingNodes(0:6)
    LOGICAL, ALLOCATABLE :: ForbiddenNodes(:)
    REAL(KIND=dp), POINTER :: Parray(:,:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), &
        SumCoord(3), AveCoord(3), Weights(3), RefScore, Score, &
        PosMeasure, NegMeasure, OffLineCoeff, DirDistance, &
        InLine, OffLine, Dist, MinDist, InLineMeasure, ScoreLimit
    CHARACTER(LEN=MAX_NAME_LEN) :: Method
!---------------------------------------------------------------

    CALL Info('FindRigidBodyFixingNodes','Starting',Level=6)

    Mesh => Solver % Mesh
    dim = Mesh % MeshDim 
    
    ALLOCATE( ForbiddenNodes(Mesh % NumberOfNodes) )
    CALL DetermineForbiddenNodes( )
    nsize = COUNT(.NOT. ForbiddenNodes) 

!   PRINT *,'Number of allowed Nodes:',nsize

    ! Find the center from the average of node positions
    !-----------------------------------------------------------
    SumCoord = 0.0_dp
    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      SumCoord = SumCoord + Coord
    END DO
    AveCoord = SumCoord / nsize


    ! Find the node closest to center and make that the new center
    !--------------------------------------------------------------
    MinDist = HUGE( MinDist ) 

    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      Dist = SUM( ( Coord - AveCoord )**2 )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        k = i  
      END IF
    END DO

    AveCoord(1) = Mesh % Nodes % x(k)
    AveCoord(2) = Mesh % Nodes % y(k)
    AveCoord(3) = Mesh % Nodes % z(k)
    IF(PRESENT(FixingDOFs)) FixingDOFs(0)=k
    

!   PRINT *,'AveCoord:',AveCoord

    ! Parameters of the search
    !-----------------------------------------------------------

    OffLineCoeff = ListGetConstReal( Solver % Values,'Fixing Nodes Off Line Coefficient',GotIt)
    IF(.NOT. GotIt) OffLineCoeff = 1.0_dp

    ScoreLimit = ListGetConstReal( Solver % Values,'Fixing Nodes Limit Score',GotIt)
    IF(.NOT. GotIt) ScoreLimit = 0.99_dp

    FixBestDirection = ListGetLogical( Solver % Values,'Fixing Nodes Axis Freeze',GotIt)

    Parray => ListGetConstRealArray( Solver % Values,'Fixing Nodes Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal = 0.0_dp
      Normal(1) = 1.0
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )      
    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    
    ! Find the fixing nodes by looping over all nodes
    !-----------------------------------------------------------
    DirDistance = 0.0_dp
    DirBest = 0
    DO dir = 1, dim
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF
      
      PosMeasure = 0.0_dp
      PosMeasureIndex = 0
      NegMeasure = 0.0_dp
      NegMeasureIndex = 0
      
      
      ! Choos the nodes within the cones in the given three directions
      !---------------------------------------------------------------
      DO i=1,Mesh % NumberOfNodes
        IF( ForbiddenNodes( i ) ) CYCLE
        
        Coord(1) = Mesh % Nodes % x(i) 
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        Coord = Coord - AveCoord
        Dist = SQRT( SUM( Coord ** 2 ) )
 
        ! Signed distance in in-line direction
        InLine = SUM( Coord * Weights )
        
        ! Distance in off-line direction 
        OffLine = SQRT( Dist**2 - InLine**2 )
        
        ! This defines a cone within which nodes are accepted
        InLineMeasure = ABS( InLine ) - OffLineCoeff * OffLine 
        IF( InLineMeasure < 0.0_dp ) CYCLE
        
        IF( InLine < 0.0_dp ) THEN
          IF( InLineMeasure > NegMeasure ) THEN
            NegMeasure = InLineMeasure
            NegMeasureIndex = i
          END IF
        ELSE           
          IF( InLineMeasure > PosMeasure ) THEN
            PosMeasure = InLineMeasure 
            PosMeasureIndex = i
          END IF
        END IF      
      END DO
      
      FixingNodes(2*dir-1) = NegMeasureIndex
      FixingNodes(2*dir) = PosMeasureIndex      

      IF( NegMeasureIndex > 0 .AND. PosMeasureIndex > 0 ) THEN
        IF( PosMeasure + NegMeasure > DirDistance ) THEN
          DirDistance = PosMeasure + NegMeasure
          DirBest = dir
        END IF
      END IF

    END DO


 
    ! To be on the safe side check that no node is used twice
    ! However, do not break the best direction
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      DO j=1,2*dim
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF        
        IF( FixingNodes(j) == FixingNodes(i) ) FixingNodes(j) = 0
      END DO
    END DO


    ! Go through the fixing nodes one-by-one and set the node so that the harmonic sum
    ! is minimized. This means that small distances are hopefully eliminated. 
    !-----------------------------------------------------------------------------------
    MaxSweep = ListGetInteger( Solver % Values,'Fixing Nodes Search Loops',GotIt)
    DO Sweep = 0,MaxSweep
      FoundBetter = .FALSE.
      DO j=1,2*dim 
        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        ! The first round set the unfixed nodes
        IF( Sweep == 0 ) THEN
!         PRINT *,'Initial Score:',j,RefScore
          IF( FixingNodes(j) /= 0 ) CYCLE
        END IF

        ! Fir the best direction because otherwise there are too 
        ! many moving parts.
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF

        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        DO i=1,Mesh % NumberOfNodes
          IF( ForbiddenNodes(i) ) CYCLE
          Score = FixingNodesScore(j,i)
          IF( Score < ScoreLimit * RefScore ) THEN
            RefScore = Score 
            FixingNodes(j) = i            
            FoundBetter = .TRUE.
          END IF
        END DO
      END DO
      IF(.NOT. FoundBetter ) EXIT
    END DO

    DO j=1,2*dim
      RefScore = FixingNodesScore(j,FixingNodes(j)) 
!     PRINT *,'Final Score:',j,RefScore
    END DO

    ! Output the selected nodes
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      j = FixingNodes(i)
      WRITE(Message,'(A,I0,3ES10.2)') 'Fixing Node: ',j,&
          Mesh % Nodes % x( j ), &
          Mesh % Nodes % y( j ), &
          Mesh % Nodes % z( j ) 
      CALL Info('FindRigidBodyFixingNodes',Message,Level=6)
      IF( PRESENT( FixingDofs ) ) FixingDofs(i) = j     
    END DO

    DEALLOCATE( ForbiddenNodes )


  CONTAINS

    !> Find the nodes that are either on interface, boundary or do not belong to the field.
    !-----------------------------------------------------------------------------------
    SUBROUTINE DetermineForbiddenNodes()

      TYPE(Element_t), POINTER :: Element
      LOGICAL, POINTER :: ig(:)
      INTEGER :: t
      
      ! Mark all interface nodes as forbidden nodes
      !-----------------------------------------------
      IF( ParEnv % PEs > 1 ) THEN
        ig => Mesh % ParallelInfo % INTERFACE
        ForbiddenNodes = ig(1:Mesh % NumberOfNodes)
      END IF

      ! Mark all nodes on boundary elements as forbidden nodes
      !--------------------------------------------------------
      DO t=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements

        Element => Mesh % Elements( t )
        ForbiddenNodes( Element % NodeIndexes ) = .TRUE.
      END DO

      ! If mask exists then add all nodes not in mask to forbidden nodes
      !-----------------------------------------------------------------
      IF( PRESENT( MaskPerm) ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF( MaskPerm(i) == 0 ) ForbiddenNodes(i) = .TRUE.
        END DO
      END IF
      
    END SUBROUTINE DetermineForbiddenNodes


    !> Give a value of goodness to the chosen fixing node.
    !-----------------------------------------------------------------------------------
    FUNCTION FixingNodesScore(direction,cand) RESULT ( Score )

      INTEGER :: direction, cand
      INTEGER :: i,j
      REAL(KIND=dp) :: Score

      REAL(KIND=dp) :: x0(3), x1(3), Dist

      IF( cand == 0 ) THEN
        Score = HUGE( Score ) 
        RETURN
      END IF

      Score = 0.0_dp
      x0(1) = Mesh % Nodes % x( cand )
      x0(2) = Mesh % Nodes % y( cand )
      x0(3) = Mesh % Nodes % z( cand )

      DO i=1,2*dim
        IF( i == direction ) CYCLE
        j = FixingNodes( i )

        ! Do not meausure distance to unset nodes!
        IF( j == 0 ) CYCLE

        ! This would lead to division by zero later on
        IF( cand == j ) THEN
          Score = HUGE( Score ) 
          RETURN
        END IF

        x1(1) = Mesh % Nodes % x( j )
        x1(2) = Mesh % Nodes % y( j )
        x1(3) = Mesh % Nodes % z( j )

        Dist = SQRT( SUM( (x0 - x1 ) ** 2 ) )
        Score = Score + 1 / Dist
      END DO

    END FUNCTION FixingNodesScore


!------------------------------------------------------------------------------
  END SUBROUTINE FindRigidBodyFixingNodes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Create a 1D mesh, may be used in 1D outlet conditions, for example.
!------------------------------------------------------------------------------
  FUNCTION CreateLineMesh( Params ) RESULT( Mesh )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    INTEGER :: i, j, k, n, NoNodes, NoElements, ActiveDirection, Order, BodyId
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName

!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info('CreateLineMesh','Creating 1D mesh on-the-fly')

!   Read in the parameters defining a uniform 1D mesh
!--------------------------------------------------------------    
    Order = ListGetInteger( Params,'1D Element Order',Found,minv=1,maxv=2)
    NoElements = ListGetInteger( Params,'1D Number Of Elements',minv=1)
    Length = ListGetConstReal( Params,'1D Mesh Length')    
    ActiveDirection = ListGetInteger( Params,'1D Active Direction',minv=-3,maxv=3)
    BodyId = ListGetInteger( Params,'1D Body Id',minv=1)
    MeshName = ListGetString( Params,'1D Mesh Name',Found)
    IF(.NOT. Found) MeshName = '1d_mesh'
    
    Mesh % Name = MeshName
    Mesh % OutputActive = .FALSE.

!   Compute the resulting mesh parameters
!--------------------------------------------------------------
    NoNodes = NoElements + 1 + NoElements * (Order - 1)    
    MeshVector = 0.0_dp
    MeshVector( ABS( ActiveDirection ) ) = 1.0_dp
    IF( ActiveDirection < 0 ) MeshVector = -MeshVector
    MeshVector = MeshVector * Length

!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )

    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
   
    DO i=1, NoNodes
      Coord = MeshVector * (i-1) / (NoNodes-1)

      x(i) = Coord(1)
      y(i) = Coord(2)
      z(i) = Coord(3)
    END DO
    

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    IF( Order == 1 ) THEN
      Elmt => GetElementType( 202 )
    ELSE
      Elmt => GetElementType( 203 )
    END IF

    DO i=1,NoElements
      Element => Mesh % Elements(i)      
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()     
      Element % ElementIndex = i
      IF( Order == 1 ) THEN
        CALL AllocateVector( Element % NodeIndexes, 2 )
        Element % Ndofs = 2
        Element % NodeIndexes(1) = i
        Element % NodeIndexes(2) = i + 1
      ELSE IF( Order == 2 ) THEN
        CALL AllocateVector( Element % NodeIndexes, 3 )
        Element % Ndofs = 3
        Element % NodeIndexes(1) = 2*i-1
        Element % NodeIndexes(2) = 2*i+1
        Element % NodeIndexes(3) = 2*i
      END IF
      
      Element % BodyId = BodyId
      Element % PartIndex = ParEnv % myPE
    END DO
    
!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = 1 + Order
    Mesh % MaxElementDOFs = 1 + Order
    Mesh % MeshDim = 1

    WRITE(Message,'(A,I0)') 'Number of elements created: ',NoElements
    CALL Info('CreateLineMesh',Message)

    WRITE(Message,'(A,I0)') 'Number of nodes created: ',NoNodes
    CALL Info('CreateLineMesh',Message)
 
    CALL Info('CreateLineMesh','All done')

  END FUNCTION CreateLineMesh


  !> Calcalate body average for a discontinuous galerkin field.
  !> The intended use is in conjunction of saving the results. 
  !> This tampers the field and therefore may have unwanted side effects
  !> if the solution is to be used for something else too.
  !-------------------------------------------------------------------
  SUBROUTINE CalculateBodyAverage( Mesh, Var, BodySum )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: BodySum

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: BodyAverage(:)
    INTEGER, ALLOCATABLE :: BodyCount(:)
    INTEGER :: n,i,j,k,l,nodeind,dgind
    REAL(KIND=dp) :: AveHits

    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    IF( BodySum ) THEN
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal sum for: '&
          //TRIM(Var % Name), Level=8)
    ELSE
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal average for: '&
          //TRIM(Var % Name), Level=8)
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE( BodyCount(n), BodyAverage(n) )


    DO i=1,CurrentModel % NumberOfBodies

      DO k=1,Var % Dofs
        BodyCount = 0
        BodyAverage = 0.0_dp

        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              BodyAverage( nodeind ) = BodyAverage( nodeind ) + &
                  Var % Values( Var % DOFs*( dgind-1)+k )
              BodyCount( nodeind ) = BodyCount( nodeind ) + 1 
            END IF
          END DO
        END DO

        IF( k == 1 ) THEN
          AveHits = 1.0_dp * SUM( BodyCount ) / COUNT( BodyCount > 0 )
          !PRINT *,'AveHits:',i,AveHits
        END IF

        IF(ParEnv % Pes>1) THEN
          CALL SendInterface(); CALL RecvInterface()
        END IF

        ! Do not average weighted quantities. They should only be summed, I guess... 
        
        IF( .NOT. BodySum ) THEN
          DO j=1,n
            IF( BodyCount(j) > 0 ) BodyAverage(j) = BodyAverage(j) / BodyCount(j)
          END DO
        END IF

        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              Var % Values( Var % DOFs*( dgind-1)+k ) = BodyAverage( nodeind ) 
            END IF
          END DO
        END DO
      END DO
    END DO

CONTAINS

     SUBROUTINE SendInterface()
       TYPE buf_t
         REAL(KIND=dp), ALLOCATABLE :: dval(:)
         INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       END TYPE buf_t

       INTEGER, ALLOCATABLE :: cnt(:)
       TYPE(buf_t), ALLOCATABLE :: buf(:)

       INTEGER :: i,j,k,ierr

       ALLOCATE(cnt(ParEnv % PEs), buf(ParEnv % PEs))

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % Interface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
         END DO
       END DO

       DO i=1,ParEnv % PEs
         ALLOCATE(buf(i) % gdof(cnt(i)), buf(i) % ival(cnt(i)), buf(i) % dval(cnt(i)))
       END DO

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % Interface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
           buf(k) % gdof(cnt(k)) = Mesh % ParallelInfo % GlobalDOFs(i)
           buf(k) % ival(cnt(k)) = BodyCount(i)
           buf(k) % dval(cnt(k)) = BodyAverage(i)
         END DO
       END DO

       DO i=1,ParEnv % PEs
         IF(.NOT. ParEnv % isNeighbour(i)) CYCLE

         CALL MPI_BSEND( cnt(i),1,MPI_INTEGER,i-1,1310,MPI_COMM_WORLD,ierr )
         IF(cnt(i)>0) THEN
           CALL MPI_BSEND( buf(i) % gdof,cnt(i),MPI_INTEGER,i-1,1311,MPI_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % ival,cnt(i),MPI_INTEGER,i-1,1312,MPI_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % dval,cnt(i),MPI_DOUBLE_PRECISION,i-1,1313,MPI_COMM_WORLD,ierr )
         END IF
       END DO
     END SUBROUTINE SendInterface


     SUBROUTINE RecvInterface()
       INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       REAL(KIND=dp), ALLOCATABLE :: dval(:)
       INTEGER :: i,j,k,ierr, cnt, status(MPI_STATUS_SIZE)

       DO i=1,ParEnv % PEs
         IF(.NOT. ParEnv % isNeighbour(i)) CYCLE

         CALL MPI_RECV( cnt,1,MPI_INTEGER,i-1,1310,MPI_COMM_WORLD,status,ierr )
         IF(cnt>0) THEN
           ALLOCATE( gdof(cnt), ival(cnt), dval(cnt) )
           CALL MPI_RECV( gdof,cnt,MPI_INTEGER,i-1,1311,MPI_COMM_WORLD,status,ierr )
           CALL MPI_RECV( ival,cnt,MPI_INTEGER,i-1,1312,MPI_COMM_WORLD,status,ierr )
           CALL MPI_RECV( dval,cnt,MPI_DOUBLE_PRECISION,i-1,1313,MPI_COMM_WORLD,status,ierr )

           DO j=1,cnt
             k = SearchNode(Mesh % ParallelInfo, gdof(j))
             IF (k>0) THEN
               BodyCount(k) = BodyCount(k) + ival(j)
               BodyAverage(k) = BodyAverage(k)  + dval(j)
             END IF
           END DO 
           DEALLOCATE( gdof, ival, dval )
         END IF
       END DO
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     END SUBROUTINE RecvInterface

  END SUBROUTINE CalculateBodyAverage



  !> Given an elemental DG field create a minimal reduced set of it that maintains
  !> the necessary continuities. The continuities may be requested between bodies
  !> or materials. Optionally the user may give a boundary mask which defines the 
  !> potential discontinuous nodes that may be greedy or not. 
  !-------------------------------------------------------------------------------
  FUNCTION MinimalElementalSet( Mesh, JumpMode, VarPerm, BcFlag, &
      NonGreedy ) RESULT ( SetPerm )

    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: JumpMode
    INTEGER, POINTER, OPTIONAL :: VarPerm(:)
    CHARACTER(LEN=*), OPTIONAL :: BcFlag
    LOGICAL, OPTIONAL :: NonGreedy
    INTEGER, POINTER :: SetPerm(:)

    TYPE(Element_t), POINTER :: Element, Left, Right
    INTEGER :: n,i,j,k,l,bc_id,mat_id,body_id,NoElimNodes,nodeind,JumpModeIndx,&
        LeftI,RightI,NumberOfBlocks
    LOGICAL, ALLOCATABLE :: JumpNodes(:)
    INTEGER, ALLOCATABLE :: NodeVisited(:)
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found
    

    CALL Info('MinimalDiscontSet','Creating discontinuous subset from DG field',Level=5)

    ! Calculate size of permutation vector
    ALLOCATE( NodeVisited( Mesh % NumberOfNodes ) )
    NodeVisited = 0

    NULLIFY( SetPerm ) 
    k = 0
    DO i=1,Mesh % NumberOfBulkElements         
      Element => Mesh % Elements(i)
      k = k + Element % TYPE % NumberOfNodes
    END DO
    CALL Info('MinimalElementalSet','Maximum number of dofs in DG: '//TRIM(I2S(k)),Level=12)
    ALLOCATE( SetPerm(k) )
    SetPerm = 0
    l = 0
    NoElimNodes = 0

    CALL Info('MinimalElementalSet','Reducing elemental discontinuity with mode: '//TRIM(JumpMode),Level=7)

    SELECT CASE ( JumpMode )

    CASE('db') ! discontinuous bodies
      NumberOfBlocks = CurrentModel % NumberOfBodies
      JumpModeIndx = 1

    CASE('dm') ! discontinuous materials
      NumberOfBlocks = CurrentModel % NumberOfMaterials
      JumpModeIndx = 2

    CASE DEFAULT
      CALL Fatal('MinimalElementalSet','Unknown JumpMode: '//TRIM(JumpMode))

    END SELECT
  

    IF( PRESENT( BcFlag ) ) THEN
      ALLOCATE( JumpNodes( Mesh % NumberOfNodes ) )
    END IF

    
    DO i=1,NumberOfBlocks
      
      ! Before the 1st block no numbers have been given.
      ! Also if we want discontinuous blocks on all sides initialize the whole list to zero. 
      IF( i == 1 .OR. .NOT. PRESENT( BcFlag ) ) THEN
        NodeVisited = 0

      ELSE
        ! Vector indicating the disontinuous nodes
        ! If this is not given all interface nodes are potentially discontinuous
        JumpNodes = .FALSE.
        
        DO j=Mesh % NumberOfBulkElements + 1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(j)

          DO bc_id=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
          END DO
          IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE
          IF( .NOT. ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

          Left => Element % BoundaryInfo % Left
          Right => Element % BoundaryInfo % Right
          IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

          IF( JumpModeIndx == 1 ) THEN
            LeftI = Left % BodyId
            RightI = Right % BodyId
          ELSE
            LeftI = ListGetInteger( CurrentModel % Bodies(Left % BodyId) % Values,'Material',Found)
            RightI = ListGetInteger( CurrentModel % Bodies(Right % BodyId) % Values,'Material',Found)
          END IF

          IF( LeftI /= i .AND. RightI /= i ) CYCLE
          JumpNodes( Element % NodeIndexes ) = .TRUE.
        END DO

        IF( PRESENT( NonGreedy ) ) THEN
          IF( NonGreedy ) THEN        
            DO j=Mesh % NumberOfBulkElements + 1, &
                Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
              Element => Mesh % Elements(j)

              DO bc_id=1,CurrentModel % NumberOfBCs
                IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
              END DO
              IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE

              IF( ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

              Left => Element % BoundaryInfo % Left
              Right => Element % BoundaryInfo % Right

              ! External BCs don't have a concept of jump, so no need to treat them
              IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

              JumpNodes( Element % NodeIndexes ) = .FALSE.
            END DO
          END IF
        END IF

        ! Initialize new potential nodes for the block where we found discontinuity
        WHERE( JumpNodes ) NodeVisited = 0
      END IF


      ! Now do the real thing. 
      ! Add new dofs such that minimal discontinuity is maintained 
      DO j=1,Mesh % NumberOfBulkElements         
        Element => Mesh % Elements(j)

        Body_Id = Element % BodyId 
        IF( JumpModeIndx == 1 ) THEN
          IF( Body_id /= i ) CYCLE
        ELSE
          Mat_Id = ListGetInteger( CurrentModel % Bodies(Body_Id) % Values,'Material',Found)
          IF( Mat_Id /= i ) CYCLE
        END IF

        NodeIndexes => Element % NodeIndexes
        
        DO k=1,Element % TYPE % NumberOfNodes         
          nodeind = NodeIndexes(k)
          IF( PRESENT( VarPerm ) ) THEN
            IF( VarPerm( nodeind ) == 0 ) CYCLE
          END IF
          IF( NodeVisited( nodeind ) > 0 ) THEN
            SetPerm( Element % DGIndexes(k) ) = NodeVisited( nodeind )
            NoElimNodes = NoElimNodes + 1
          ELSE
            l = l + 1
            NodeVisited(nodeind) = l
            SetPerm( Element % DGIndexes(k) ) = l
          END IF
        END DO
      END DO
    END DO

    CALL Info('MinimalElementalSet','Independent dofs in elemental field: '//TRIM(I2S(l)),Level=7)
    CALL Info('MinimalElementalSet','Redundant dofs in elemental field: '//TRIM(I2S(NoElimNodes)),Level=7)     

  END FUNCTION MinimalElementalSet


  !> Calculate the reduced DG field given the reduction permutation.
  !> The permutation must be predefined. This may be called repeatedly
  !> for different variables. Optionally one may take average, or 
  !> a plain sum over the shared nodes. 
  !-------------------------------------------------------------------
  SUBROUTINE ReduceElementalVar( Mesh, Var, SetPerm, TakeAverage )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: TakeAverage

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: SetSum(:)
    INTEGER, ALLOCATABLE :: SetCount(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp) :: AveHits

    IF(.NOT. ASSOCIATED(var)) THEN
      CALL Warn('ReduceElementalVar','Variable not associated!')
      RETURN
    END IF

    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) THEN
      CALL Warn('ReduceElementalVar','Var % Perm too small!')
      RETURN
    END IF

    IF( TakeAverage ) THEN
      CALL Info('CalculateSetAverage','Calculating reduced set average for: '&
          //TRIM(Var % Name), Level=7)
    ELSE
      CALL Info('CalculateSetAverage','Calculating reduced set sum for: '&
          //TRIM(Var % Name), Level=7)
    END IF

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    ALLOCATE( SetCount(m), SetSum(m) )
    SetCount = 0
    SetSum = 0.0_dp

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs
      SetCount = 0
      SetSum = 0.0_dp

      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        SetSum(j) = SetSum(j) + Var % Values( Var % DOFs * (l-1) + dof )
        SetCount(j) = SetCount(j) + 1
      END DO
        
      IF( TakeAverage ) THEN
        WHERE( SetCount > 0 ) SetSum = SetSum / SetCount
      END IF

      IF( dof == 1 ) THEN
        AveHits = 1.0_dp * SUM( SetCount ) / COUNT( SetCount > 0 )
        PRINT *,'AveHits:',AveHits
      END IF

      ! Copy the reduced set back to the original elemental field
      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        Var % Values( Var % DOFs * (l-1) + dof ) = SetSum(j)
      END DO
    END DO

  END SUBROUTINE ReduceElementalVar


  !> Given a elemental DG field and a reduction permutation compute the 
  !> body specific lumped sum. The DG field may be either original one
  !> or already summed up. In the latter case only one incident of the 
  !> redundant nodes is set.
  !---------------------------------------------------------------------
  SUBROUTINE LumpedElementalVar( Mesh, Var, SetPerm, AlreadySummed )
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: AlreadySummed

    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: NodeVisited(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp), ALLOCATABLE :: BodySum(:)

    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    CALL Info('LumpedElementalVar','Calculating lumped sum for: '&
        //TRIM(Var % Name), Level=8)

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    IF( AlreadySummed ) THEN
      ALLOCATE( NodeVisited(m) )
    END IF
    ALLOCATE( BodySum( CurrentModel % NumberOfBodies ) )

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs

      BodySum = 0.0_dp

      DO i=1,CurrentModel % NumberOfBodies

        IF( AlreadySummed ) THEN
          NodeVisited = .FALSE.
        END IF

        DO j=1,Mesh % NumberOfBulkElements         
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE

          DO k=1,Element % TYPE % NumberOfNodes         
            dgind = Element % DGIndexes(k)
            l = SetPerm(dgind)
            IF( l == 0 ) CYCLE

            IF( AlreadySummed ) THEN
              IF( NodeVisited(l) ) CYCLE           
              NodeVisited(l) = .TRUE.
            END IF

            BodySum(i) = BodySum(i) + &
                Var % Values( Var % Dofs * ( Var % Perm( dgind )-1) + dof )
          END DO
        END DO
      END DO

      IF( Var % Dofs > 1 ) THEN
        CALL Info('LumpedElementalVar','Lumped sum for component: '//TRIM(I2S(dof)),Level=6)
      END IF
      DO i=1,CurrentModel % NumberOfBodies
        PRINT *,'BodySum',i,BodySum(i)
      END DO

    END DO

    DEALLOCATE( NodeVisited, BodySum )

  END SUBROUTINE LumpedElementalVar

!------------------------------------------------------------------------------
END MODULE MeshUtils
!------------------------------------------------------------------------------

!> \}

